library(conflicted)
library(haven)
library(igraph)
library(magrittr)
conflicts_prefer(magrittr::extract)
conflicts_prefer(magrittr::set_names)
library(tidyverse)
conflicts_prefer(dplyr::filter)

## This function executes the complete analysis of the dataset required to create version 2 of the WHO catalogue.
## If fast = TRUE, then only the MAIN analysis (containing the WHO and ALL datasets, not CC[ATU]) gets processed.
## If correct_all = TRUE, then the analyses correct for the total number of variants considered; if it is FALSE,
## they only correct for the number of variants classified S or R, i.e. appear as a SOLO in at least one isolate.
## If skipBDQfromSA = TRUE, then a specified subset of samples are not carried forward for Bedaquiline analysis. 
## If LoF = TRUE, then pooled loss-of-function mutations are also considered hypotheses for the FDR corrections.
## If safe = TRUE (slow option), then a full conversion of variants from V1 to V2 is also performed from scratch.
## If orphanByDrug = TRUE, then the listing of orphan genotypes is stratified by the drug; this is more accurate.
## The remaining options are used to specify the relevant directories and should be left at their default values.
mainDriver = function(fast = FALSE, correct_all = TRUE, skipBDQfromSA = FALSE, LoF = TRUE, safe = TRUE, orphanByDrug = TRUE,
    EXTRACTION_ID = "2023-04-25T06_00_10.443990_jr_b741dc136e079fa8583604a4915c0dc751724ae9880f06e7c2eacc939e086536", 
    OUTPUT_DIRECTORY = paste0("Results/", EXTRACTION_ID), SCRIPT_LOCATION_DIR = "SOLOport/",
    DATA_DIRECTORY = "SOLO Algorithm Input files/DATABASE EXTRACTION files/", NON_DATABASE_DIRECTORY = "SOLO Algorithm Input files/STATA DTA files/") {
  DATA_DIRECTORY %<>%
    normalizePath()
  NON_DATABASE_DIRECTORY %<>%
    normalizePath()
  for (script in list.files(path = str_remove(SCRIPT_LOCATION_DIR, "/$"), full.names = TRUE, pattern = ".R$") %>% 
       setdiff(paste0(str_remove(SCRIPT_LOCATION_DIR, "/$"), c("/Build.R", "/Driver.R", "/LegacyCode.R")))) { 
         print(paste("Sourcing", script))
         source(script)
  }
  OUTPUT_DIRECTORY = str_replace(OUTPUT_DIRECTORY, "Results", ifelse(skipBDQfromSA, "Results_withoutSA", "Results"))
  dir.create(OUTPUT_DIRECTORY, recursive = TRUE)
  OUTPUT_DIRECTORY %<>%
    normalizePath()
  setwd(OUTPUT_DIRECTORY)
  ## Initial preprocessing of the genotypes: extract all the genotypes,  recording the drug and the tier
  allGenotypes  = extractData(inDir = paste(str_remove(DATA_DIRECTORY, "/$"), str_remove(EXTRACTION_ID, "/$"), "full_genotypes/", sep = "/"),
    drugList = DRUG_LIST, geno = TRUE)
  ## Make sure there are no missing or empty resolved symbols
  stopifnot(!any(is.na(allGenotypes$resolved_symbol)))
  ## Rename three columns for convenience and remove an unused column
  allGenotypes %<>%
    rename(gene = 'resolved_symbol', mutation = 'variant_category', effect = 'predicted_effect') %>%
    select(-neutral)
  ## Replace missing mutations or effects by "missing", then create a variant out of gene and mutation
  allGenotypes %<>%
    mutate_at(c("mutation", "effect"), ~{replace_na(., "missing")}) %>%
    mutate(variant = ifelse(mutation == "missing", "missing", paste(gene, mutation, sep = "_")))
  ## Mark any genotype entry with a low MAF or a missing variant (i.e. a sequencing defect in the gene) as a het
  allGenotypes %<>%
    mutate(het = (`max(af)` < MAF_THRESHOLD_REGULAR | variant == "missing"))
  ## Compute the positions for the mutations, first making sure that there are no more than 2 position values
  stopifnot(all(str_count(allGenotypes$mutation, "[0-9]+") <= 2L))
  allGenotypes %<>%
    mutate(allPos = map(mutation, ~{unlist(str_extract_all(., "[-]?[0-9]+"))})) %>%
    mutate(pos1 = map_int(allPos, ~{as.integer(first(.))}), pos2 = map_int(allPos, ~{as.integer(nth(., 2))})) %>% 
    select(-allPos)
  ## Initial consistency checks on the genotypes: the initial amino acid is always consistent for a given gene-position combination
  allGenotypes %<>%
    mutate(aa_start = ifelse(str_starts(mutation, "p.") & is.na(pos2), str_sub(mutation, 3, 5), NA))
  stopifnot(all(testConsistent(allGenotypes %>% filter(!is.na(aa_start)), c("gene", "pos1"), consistentVars = c("aa_start"))[[1]]))
  ## Initial consistency checks on the genotypes: the mutation is missing if and only if the effect is missing
  stopifnot(all((allGenotypes$mutation == "missing") == (allGenotypes$effect == "missing")))
  ## Initial consistency checks on the genotypes: same effect and position for each variant-drug combination
  stopifnot(all(testConsistent(allGenotypes, c("drug", "variant"), consistentVars = c("effect", "position"))[[1]]))
  ## Initial consistency checks on the genotypes: same tier for each gene-drug combination
  stopifnot(all(testConsistent(allGenotypes, c("drug", "gene"   ), consistentVars = "tier")[[1]]))
  ## Initial preprocessing of the phenotypes: extract all the phenotypes, recording the drug
  allPhenotypes = extractData(inDir = paste(str_remove(DATA_DIRECTORY, "/$"), str_remove(EXTRACTION_ID, "/$"), "phenotypes/", sep="/"),
    drugList = DRUG_LIST, geno = FALSE)
  ## Rename one column for convenience and remove an unused column
  allPhenotypes %<>%
    rename(category_phenotype = 'phenotypic_category') %>%
    select(-box)
  ## Initial consistency checks on the phenotypes: phenotype consistent for sample-drug combinations per category
  stopifnot(testConsistent(allPhenotypes, c("sample_id", "category_phenotype", "drug"), consistentVars = "phenotype")[[1]])
  ## Initial consistency checks on the phenotypes: unique isolate-drug combinations between the WHO, ALL datasets
  stopifnot(anyDuplicated(allPhenotypes %>% dplyr::filter(category_phenotype %in% c("WHO", "ALL")) %>% select(drug, sample_id)) == 0)
  ## Merge the genotypes and the phenotypes after separating the latter by phenotypic category group
  fullDataset = mergeGenoPheno(allGenotypes, allPhenotypes, phenoGroups = PHENO_GROUPS)
  ## Alternative based on a bedaquiline issue: if skipBDQfromSA is TRUE, remove the Bedaquiline entries corresponding to a list of sample_id's
  if (skipBDQfromSA) {
    curDir = getwd()
    setwd(DATA_DIRECTORY)
    badIDs = read_csv("query_result_2023-04-13T15_28_45.986784Z.csv", show_col_types = FALSE, guess_max = Inf) %>%
      pull(`Sample ID`)
    for (name in names(fullDataset)) {
      fullDataset[[name]] %<>%
        filter(!(sample_id %in% badIDs & drug == "Bedaquiline"))
    }
    setwd(curDir)
  }
  ## Quality control based on a sanity check: S isolates with one of specific variant-drug combinations fail QC
  ## For bookkeeping purposes we keep additional columns even though only the sample_ids are used for exclusion
  samplesToExclude = fullDataset[["MAIN"]] %>%
    inner_join(BAD_VAR_DRUG_PAIRS, by = c("drug", "variant")) %>%
    dplyr::filter(phenotype == "S") %>%
    mutate(neutral = FALSE) %>%
    distinct(sample_id, drug, variant, gene, phenotype, het, effect, tier, neutral, category_phenotype)
  ## Record the sample_ids that need to be excluded
  write_csv(samplesToExclude, "excluded_after_qc.csv")
  ## Set the data with WHO phenotype and add it to the dataset; it will be used to determine neutral variants
  WHOData = fullDataset[["MAIN"]] %>%
    dplyr::filter(category_phenotype == "WHO")
  fullDataset[["WHO"]] = WHOData
  ## Replace the phenotype in the MAIN dataset with "ALL"
  fullDataset[["MAIN"]] %<>%
    mutate(category_phenotype = "ALL")
  ## (Re)run the neutral algorithm to create the files containing neutral mutations (note that this is a side effect)
  if (safe) {
    WHODataMarked = WHOData %>%
      dplyr::filter(!(sample_id %in% samplesToExclude$sample_id)) %>%
      neutralAlgorithm(safe = safe, skipBDQfromSA = skipBDQfromSA, NON_DATABASE_DIRECTORY = str_remove(NON_DATABASE_DIRECTORY, "/$"), 
                       DATA_DIRECTORY = str_remove(DATA_DIRECTORY, "/$"), EXTRACTION_ID = str_remove(EXTRACTION_ID, "/$"))
  }
  ## Extract the prepared list of all neutral mutations
  allNeutrals = read_csv(paste0("neutral_mutations_WHO_", ifelse(skipBDQfromSA, "withoutSA_", ""), "F.csv"), 
                         guess_max = Inf, show_col_types = FALSE) %>%
    mutate(neutral = TRUE)
  ## Now identify the neutral variants and mark these, as well as other relevant ones, in the entire dataset
  for (name in names(fullDataset)) {
    if (!fast || name %in% c("MAIN", "WHO")) {
      print(paste("Preprocessing the", name, "dataset"))
      ## Extract the dataset corresponding to the name
      curSet = fullDataset[[name]]
      ## Compute the per-drug denominators (number of R/S isolates screened for each drug); account for exclusion
      allDenominators = curSet %>%
        dplyr::filter(!(sample_id %in% samplesToExclude$sample_id)) %>%
        group_by(drug) %>%
        distinct(sample_id, .keep_all = TRUE) %>%
        mutate(RDen = sum(phenotype == "R"), SDen = sum(phenotype == "S")) %>%
        slice(1) %>%
        ungroup() %>%
        select(drug, RDen, SDen)
      curSet %<>%
        inner_join(allDenominators, by = "drug")
      ## Adjust for the numbers of R/S isolates screened for each drug but excluded due to QC just for the bad pairs
      extraDenominators = curSet %>%
        dplyr::filter(sample_id %in% samplesToExclude$sample_id) %>%
        group_by(drug) %>%
        distinct(sample_id, .keep_all = TRUE) %>%
        mutate(RDen = sum(phenotype == "R"), SDen = sum(phenotype == "S")) %>%
        slice(1) %>%
        ungroup() %>%
        select(-variant) %>%
        inner_join(BAD_VAR_DRUG_PAIRS, by = "drug") %>%
        select(drug, variant, RDen, SDen)
      curSet %<>%
        left_join(extraDenominators, by = c("drug", "variant")) %>%
        adjustDuplicateColumns(warn = FALSE, add = TRUE)
      ## Mark the neutral mutations and convert the neutral variable to a logical one
      curSet %<>%
        left_join(allNeutrals, by = c("drug", "variant")) %>%
        mutate_at(c(paste0("set", LETTERS[1:5]), "lit_mutation", "version1", "neutral"), convertToLogical)
      ## And finally, replace into the full dataset
      fullDataset[[name]] = curSet
    }
  }
  ## Save the clean version of the input data for downstream analysis of the catalogue's performance on itself 
  write_csv(fullDataset[["MAIN"]], paste0("CompleteDataset", ifelse(skipBDQfromSA, "withoutSA_", ""), ".csv"))
  ## Now process the entire dataset in the specified order; additionally perform each of the 3 types of pooling
  for (name in c("WHO", "MAIN", "CC", "ATU")) {
    if (!fast || name %in% c("MAIN", "WHO")) {
      ## Placeholder for the primary statistics
      fullStats = tibble()
      ## Sanity check: ensures that there is only one non-missing variant of each name per drug-sample combination
      stopifnot(anyDuplicated(fullDataset[[name]] %>% dplyr::filter(variant != "missing") %>% select(drug, sample_id, variant)) == 0)
      for (pool in c(NA, names(POOLED_EFFECTS))) {
        POOLED = !is.na(pool)
        print(paste("Processing the combination of", name, "and", ifelse(POOLED, pool, "no"), "pooling"))
        ## Extract the dataset corresponding to the name
        curSet = fullDataset[[name]]
        ## Exclude the samples failing QC, after saving the rows corresponding to the problematic mutations first
        ## The samples that were among the originally excluded ones are accounted for as well for completeness
        curExcludedEntries = curSet %>%
          dplyr::filter(sample_id %in% samplesToExclude$sample_id) %>%
          bind_rows(samplesToExclude %>% dplyr::filter(name == "MAIN" | (category_phenotype == name))) %>%
          distinct(sample_id, drug, variant, gene, phenotype, het, effect, tier, neutral)
        ## Count frequency and solo frequency among the excluded samples, in the entries with "QC-relevant" mutations
        curExtraCounts = curExcludedEntries %>%
          dplyr::filter(!het) %>%
          inner_join(BAD_VAR_DRUG_PAIRS, by = c("drug", "variant")) %>%
          group_by(drug, variant) %>%
          mutate(present_R = sum(phenotype == "R"), present_S = sum(phenotype == "S")) %>%
          slice(1) %>%
          ungroup() %>%
          select(drug, variant, present_R, present_S)
        curExtraSolos = curExcludedEntries %>%
          prepMask(Syn = TRUE, Tier2 = TRUE, Neutral = TRUE, Pool = NA, SOnly = TRUE, aS = FALSE) %>%
          runSOLOPipeline(maxIter = 1L, stage = NULL, removeSOnly = FALSE) %>%
          inner_join(BAD_VAR_DRUG_PAIRS, by = c("drug", "variant")) %>%
          select(drug, variant, Rcnt, Scnt)
        curSet %<>%
          dplyr::filter(!(sample_id %in% samplesToExclude$sample_id))
        ## Add the pooling if necessary; otherwise, process "stage 0" which yields SOLO counts for neutral variants
        if (POOLED) {
          ## Update the variants being pooled accordingly; compress identical variants per drug-sample-gene combo
          ## Update the rows of the variants being pooled with missing data in the variables that are undefined
          curSet %<>%
            mutate_at("variant", ~{ ifelse(effect %in% POOLED_EFFECTS[[pool]] & !het, paste0(gene, "_", pool), .) }) %>%
            mutate_at(c("max(af)", "position"), ~{ ifelse(str_ends(variant, paste0("_", pool)), NA, .) }) %>%
            mutate_at("neutral", ~{ ifelse(str_ends(variant, paste0("_", pool)), FALSE, .) }) %>%
            distinct(drug, sample_id, gene, variant, .keep_all = TRUE)
        } else {
          ## Prepare stage 0: mask all synonymous and tier 2 variants, keep the neutral ones, and ignore the pooling
          zerothSet = prepMask(curSet, Syn = TRUE, Tier2 = TRUE, Neutral = FALSE, aS = FALSE, Pool = NA, SOnly = TRUE)
          ## Run the SOLO pipeline for 1 step of the algorithm; record stage as 0
          zerothOutputs = runSOLOPipeline(zerothSet, maxIter = 1L, stage = 0L, removeSOnly = FALSE)
          print(paste("Stage 0:", sum(zerothOutputs$class == "S"), "S and", sum(zerothOutputs$class == "R"), "R variants"))
          ## Remove the classification (we will only use SOLO counts), then filter down to only the neutral variants
          zerothOutputs %<>%
            select(-class) %>%
            inner_join(allNeutrals %>% select(drug, variant), by = c("drug", "variant"))
          ## Merge the original data with the new SOLO counts
          curSet %<>%
            left_join(zerothOutputs, by = c("drug", "variant"))
        }
        ## Prepare stage 1: mask all neutral, synonymous and tier 2 variants; SOnly is FALSE for multiple iterations!
        firstSet = prepMask(curSet, Syn = TRUE, Tier2 = TRUE, Neutral = TRUE, aS = FALSE, Pool = pool, SOnly = ifelse(POOLED, TRUE, FALSE))
        ## Run the SOLO pipeline for 2 (unpooled) or 1 (pooled) steps of the algorithm; record stage as 1
        ## New fixes introduced: removes any missing variants that get classified and adds possible het values
        firstOutputs = runSOLOPipeline(firstSet, maxIter = ifelse(POOLED, 1L, 2L), stage = 1L, removeSOnly = TRUE) %>%
          filter(variant != "missing") %>%
          inner_join(HET_TAB, by = "class", relationship = "many-to-many")
        print(paste("Stage 1:", sum(firstOutputs$class == "S")/HET_CNT["S"], "S and", sum(firstOutputs$class == "R")/HET_CNT["R"], "R variants"))
        ## Merge the original data with the new classification statuses
        curSet %<>%
          left_join(firstOutputs, by = c("drug", "variant", "het")) %>%
          adjustDuplicateColumns()
        ## The second stage of the algorithm brings in tier 2 variants for the unexplained R samples
        auxTab = findUnexplained(curSet)
        ## Prepare for stage 2: mask covers the neutral or synonymous effect variants, but not those in tier 2
        secondSet = prepMask(auxTab, Syn = TRUE, Tier2 = FALSE, Neutral = TRUE, aS = TRUE, Pool = pool, SOnly = FALSE)
        ## Run the SOLO pipeline for 1 step of the algorithm; record stage as 2
        secondOutputs = runSOLOPipeline(secondSet, maxIter = 1L, stage = 2L, removeSOnly = TRUE)
        ## Remove the classification of any pairs that were analysed in the first stage (i.e. tier 1 variants)
        ## New fixes introduced: removes any missing variants that get classified and adds possible het values
        secondOutputs %<>%
          filter(variant != "missing") %>%
          anti_join(firstOutputs, by = c("drug", "variant")) %>%
          inner_join(HET_TAB, by = "class", relationship = "many-to-many")
        print(paste("Stage 2:", sum(secondOutputs$class == "S")/HET_CNT["S"], "S and", sum(secondOutputs$class == "R")/HET_CNT["R"], "R variants"))
        ## Merge the original data with the new mutation statuses; adjust duplicated columns in order of creation
        curSet %<>%
          left_join(secondOutputs, by = c("drug", "variant", "het")) %>%
          adjustDuplicateColumns()
        if (!POOLED) {
          ##  The third stage brings in silent tier 1 variants for the still unexplained R samples
          auxTabNew = findUnexplained(curSet)
          ## Prepare for stage 3: mask the neutral variants or those in tier 2, but unmask synonymous tier 1 variants
          thirdSet = prepMask(auxTabNew, Syn = FALSE, Tier2 = TRUE, Neutral = TRUE, aS = TRUE, Pool = NA, SOnly = TRUE)
          ## Run the SOLO pipeline for 1 step of the algorithm; record stage as 3
          thirdOutputs = runSOLOPipeline(thirdSet, maxIter = 1L, stage = 3L, removeSOnly = TRUE)
          ## But remove anything that was already classified before (in tier 1) and only keep the R classifications
          thirdOutputs %<>%
            anti_join(firstOutputs, by = c("drug", "variant")) %>%
            mutate_at("class", ~{ ifelse(. != "R", "U", "R") })
          print(paste("Stage 3:", sum(thirdOutputs$class == "S"), "S and", sum(thirdOutputs$class == "R"), "R variants"))
          ## Merge the original data with the new mutation statuses; adjust duplicated columns in order of creation
          curSet %<>%
            left_join(thirdOutputs, by = c("drug", "variant")) %>%
            adjustDuplicateColumns()
        }
        ## Calculate the frequency of each variant in R isolates (present_R), in S isolates (present_S), and overall
        ## (present), and corresponding variables for the complement (absent), as well as SOLO counts per phenotype.
        ## Note that hets (ie variants in which we don't have sufficient confidence) are removed at the beginning! 
        ## Also adjust the variant and solo counts for samples excluded due to QC by adding the extra counts to each.
        curStats = curSet %>%
          dplyr::filter(!het & variant != "missing") %>%
          rename(algorithm_pass = iter, datasets = category_phenotype) %>%
          group_by(drug, variant) %>%
          mutate(present = n(), present_R = sum(phenotype == "R"), present_S = present - present_R) %>%
          left_join(curExtraCounts, by = c("drug", "variant")) %>%
          adjustDuplicateColumns(warn = FALSE, add = TRUE) %>%
          mutate(absent_R = RDen - present_R, absent_S = SDen - present_S) %>%
          left_join(curExtraSolos, by = c("drug", "variant")) %>%
          adjustDuplicateColumns(warn = FALSE, add = TRUE) %>%
          ungroup() %>%
          distinct(drug, variant, .keep_all = TRUE)
        ## Finally, complete the extraction and mark variants to FDR-correct for in SOLO and ALL modes, respectively.
        curStats %<>%
          mutate(correctAll = (!(lit_mutation | version1 | is.na(effect) | effect %in% SILENT_EFFECTS))) %>%
          rename(SOLO_R = Rcnt, SOLO_S = Scnt) %>%
          mutate(SOLO_SorR = SOLO_R + SOLO_S, correctSOLO = correctAll & !is.na(SOLO_SorR) & SOLO_SorR > 0) %>%
          select(drug, variant, class, tier, algorithm_pass, neutral, datasets, effect, position, stage, lit_mutation, version1, 
                 starts_with("present"), starts_with("absent"), starts_with("SOLO"), starts_with("correct"), starts_with("pos"), starts_with("set"))
        ## Aggregate the primary statistics with the ones computed previously, but only add new variant-drug pairs
        if (nrow(fullStats) > 0) {
          curStats %<>%
            anti_join(fullStats, by = c("drug", "variant"))
        }
        stopifnot(!POOLED || all(str_ends(curStats$variant, paste0("_", pool))))
        if (POOLED) {
          curStats %<>%
            mutate(correctAll = LoF, correctSOLO = LoF)
        }
        fullStats %<>%
          bind_rows(curStats)
        ## Save the resulting dataset
        fullDataset[[paste0(name, "_", ifelse(!POOLED, "unpooled", pool))]] = curSet
      }
      ## Save the basic mutation-wise statistics into a file
      curFilename = paste0("Stats_", name, ifelse(skipBDQfromSA, "_withoutSA", ""), ifelse(LoF, "_withLoFs", ""), ".csv")
      write_csv(fullStats, paste0("Basic_", curFilename))
      if (!fast) {
        ## Now compute the derived statistics and save them into a file
        fullStats %<>%
          computeCatalogueStats(correct_all = correct_all)
        write_csv(fullStats, curFilename)
      }
    }
  }
  print("Computing the final grades of all mutations")
  finalResult = gradeMutations(skipBDQfromSA = skipBDQfromSA, LoF = LoF, NON_DATABASE_DIRECTORY = str_remove(NON_DATABASE_DIRECTORY, "/$"))
  gradedFilename = paste0(paste("Final_graded_algorithm_catalogue", Sys.Date(), "Leonid", sep = "_"), 
                          ifelse(skipBDQfromSA, "_withoutSA", ""), ifelse(LoF, "_withLoFs", ""), ".csv")
  ## finalResult = prepareOutputForPaolo(gradedFilename, skipBDQfromSA = skipBDQfromSA, LoF = LoF)
  setwd(paste(str_remove(DATA_DIRECTORY, "/$"), str_remove(EXTRACTION_ID, "/$"), "orphan_genotypes/", sep = "/"))
  orphanTab = read_csv(list.files()[1], guess_max = Inf, show_col_types = FALSE) %>%
    rename(gene = 'resolved_symbol', mutation = 'variant_category', effect = 'predicted_effect', drug = 'drug_name') %>%
    mutate(variant = paste(gene, mutation, sep = "_"))
  if (orphanByDrug) {
    orphanTab %<>%
      select(variant, drug, sample_id) %>%
      group_by(variant, drug)
  } else {
    orphanTab %<>%
      select(variant, sample_id) %>%
      group_by(variant)
  }
  orphanTab %<>%
    mutate(Present_NoPheno = n_distinct(sample_id)) %>%
    slice(1) %>%
    ungroup() %>%
    select(-sample_id)
  print(paste("Adding in the", nrow(orphanTab), "orphan mutations"))
  setwd(OUTPUT_DIRECTORY)
  if (orphanByDrug) {
    finalResult %<>%
      full_join(orphanTab, by = c("variant", "drug"))
  } else {
    finalResult %<>%
      full_join(orphanTab, by = "variant")
  }
  ## finalResult %<>% mutate_at("Present_NoPheno", ~{replace_na(., 0)})
  write_csv(finalResult, paste0(gradedFilename, "_withOrphans", ifelse(orphanByDrug, "_Stratified", ""), ".csv"))
  setwd("../../")
  fullDataset
}
