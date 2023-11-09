## Auxiliary function to extract genotypes (if geno = TRUE) or phenotypes (if geno = FALSE) stored within inDir
extractData = function(inDir, drugList = NULL, geno = FALSE) {
  ## Record the starting directory
  startDir = getwd()
  ## Change into the input directory
  setwd(inDir)
  ## Identify sub-directories
  useDirs = list.files()
  ## Only keep the sub-directories matching the drug list, if specified
  if (!is.null(drugList)) {
    useDirs %<>%
      intersect(paste0("drug_name=", drugList))
  }
  ## Update the drug list accordingly
  drugList = str_remove_all(useDirs, "drug_name=")
  ## Iteration counter
  N = length(useDirs)
  ## Initially empty output table
  outputTab = tibble()
  ## Iterate over drugs
  for (index in 1:N) {
    ## Move into the sub-directory
    setwd(useDirs[[index]])
    if (geno) {
      ## For genotypes, list the tiers
      curLF = list.files()
      ## Extract the tier identifiers, convert to integers
      curTiers = str_remove_all(curLF, "tier=") %>%
        as.integer()
      ## Inner iteration counter
      M = length(curTiers)
      for (ind in 1:M) {
        ## Move into the sub-sub-directory
        setwd(curLF[[ind]])
        ## Identify the file, check it is unique
        curFile = list.files()
        stopifnot(length(curFile) == 1)
        ## Parse it and append it to the output table, recording the drug and the tier
        outputTab %<>%
          bind_rows(read_csv(curFile, guess_max = Inf, show_col_types = FALSE) %>% 
                      mutate(drug = drugList[index], tier = curTiers[ind]))
        ## Go back up a level
        setwd("../")
      }
    } else {
      ## Identify the file, check it is unique
      curFile = list.files()
      stopifnot(length(curFile) == 1)
      ## Parse it and append it to the output table, recording the drug
      outputTab %<>%
        bind_rows(read_csv(curFile, guess_max = Inf, show_col_types = FALSE) %>% 
                    mutate(drug = drugList[index]))
    }
    ## Go back up a level
    setwd("../")
  }
  ## Go back to the starting directory
  setwd(startDir)
  outputTab
}

## Auxiliary function for combining genotypes and phenotypes, then splitting them according to phenotype group
mergeGenoPheno = function(Genotypes, Phenotypes, phenoGroups = NULL) {
  Phenotypes %<>%
    inner_join(PHENO_GROUPS, by = "category_phenotype")
  fullDataset = right_join(Genotypes, Phenotypes, by = c("sample_id", "drug"), relationship = "many-to-many") %>%
    mutate_at("het", ~{ ifelse(is.na(variant), TRUE, .) })
  fullDataset %<>%
    split(fullDataset$group)
  fullDataset
}

## Auxiliary function for converting all Stata input files into CSV format, keeping names but changing extensions
convertAllStataFiles = function(inputDir = "SOLO Algorithm Input files/STATA DTA files/") {
  initDir = getwd()
  setwd(inputDir)
  LF = list.files(".", pattern = ".dta$")
  for (fname in LF) {
    Tab = read_dta(fname)
    print(fname)
    write_csv(Tab, paste0(initDir, "/", str_replace(fname, ".dta$", ".csv")))
  }
  setwd(initDir)
}

## Auxiliary function for testing that all the groups defined by groupingVars agree on consistentVars in a Table
testConsistent = function(Table, groupingVars, consistentVars) {
  L = length(consistentVars)
  if (L == 0) { return(TRUE) }
  cnames = colnames(Table)
  stopifnot(all(groupingVars %in% cnames))
  stopifnot(all(consistentVars %in% cnames))
  groupedTab = Table %>%
    group_by(across(all_of(groupingVars)))
  checks   = rep(FALSE, L) %>%
    magrittr::set_names(consistentVars)
  problems = vector("list", L) %>%
    magrittr::set_names(consistentVars)
  for (ind in 1:L) {
    curVar = consistentVars[ind]
    miniTab = groupedTab %>%
      select(any_of(curVar)) %>%
      distinct() %>%
      mutate(N = n()) %>%
      ungroup()
    if (!(all(miniTab$N == 1))) {
      problems[[ind]] = miniTab %>%
        dplyr::filter(N > 1)
    } else {
      checks[ind] = TRUE
    }
  }
  output = list(checks = checks, problems = problems)
}

## Auxiliary function for processing files (usually derived from Excel spreadsheets) with two header rows
## If prefix = TRUE, prepends the first header row (propagated rightwards) to the second header row
## Otherwise, directly uses the second header row whenever it is non-missing
read_csv_2headers = function(inputFile, prefix = FALSE) {
  initTab = read_csv(inputFile, guess_max = Inf, show_col_types = FALSE)
  ## Merge the second row of headers into the first one
  initHeaders = initTab %>% 
    slice(1)
  coln = colnames(initTab)
  if (prefix) { 
    for (ind in 2:length(coln)) {
      if (is.na(coln[ind])) {
        coln[ind] = coln[ind - 1]
      }
    }
    coln[!is.na(initHeaders)] = paste(coln[!is.na(initHeaders)], initHeaders[!is.na(initHeaders)], sep = "_")
  } else {
    coln[!is.na(initHeaders)] = initHeaders[!is.na(initHeaders)]
  }
  initTab %<>%
    slice(-1) %>%
    set_colnames(coln)
  initTab
}

## Auxiliary function to "manually" recode the values in initTab according to the additionally supplied manualTab
recodeValues = function(initTab, manualTab) {
  stopifnot(ncol(manualTab) == 2)
  coln = colnames(manualTab)
  stopifnot(all(str_ends(coln, "_first", negate = TRUE)) && all(str_ends(coln, "_second", negate = TRUE)))
  ## Create the completed version of the sub-table, taking the merged column's values from manualTab
  otherTab = initTab %>%
    inner_join(manualTab, by = coln[1], suffix = c("_first", "_second")) %>%
    rename_with(.fn = function(x) { str_remove_all(x, "\\_second$") }) %>%
    select(!any_of(paste0(coln[2], "_first")))
  ## Now remove the problematic entries in the original table and add them back in from the completed sub-table
  initTab %<>%
    anti_join(otherTab, by = coln[1]) %>%
    bind_rows(otherTab)
  initTab
}

## Auxiliary function for simplifying columns obtained during a merge (the second one is used if the first is NA).
## If warn = TRUE, a warning is printed out if any pair of parallel columns has non-missing elements in common.
## If add = TRUE, any pair of parallel columns is assumed to have numeric values and is added to get the result,
## and if the duplicate columns are not numeric, then their values are concatenated using a semicolon separator.
adjustDuplicateColumns = function(initTab, suffixes = c(".x", ".y"), warn = TRUE, add = FALSE) {
  coln = colnames(initTab)
  coln1 = coln[str_ends(coln, suffixes[1])]
  coln2 = coln[str_ends(coln, suffixes[2])]
  coln1Short = str_sub(coln1, end = -(nchar(suffixes[1]) + 1))
  coln2Short = str_sub(coln2, end = -(nchar(suffixes[1]) + 1))
  stopifnot(all(sort(coln1Short) == sort(coln2Short)))
  for (ind in seq_along(coln1Short)) {
    Col = coln1Short[ind]
    vec1 = initTab %>%
      pull(all_of(paste0(Col, suffixes[1])))
    vec2 = initTab %>%
      pull(all_of(paste0(Col, suffixes[2])))
    if (warn && any(vec1 != vec2, na.rm = TRUE)) {
      print(paste("Warning: column", Col, "has", sum(vec1 != vec2, na.rm = TRUE), "conflicting entries!"))
    }
    if (add) {
      if (is.numeric(vec1) && is.numeric(vec2)) {
        vec12 = ifelse(is.na(vec1), 0, vec1) + ifelse(is.na(vec2), 0, vec2)
      } else {
        vec12 = vec1
        vec12[!is.na(vec2)] = vec2[!is.na(vec2)]
        vec12[!is.na(vec1) & !is.na(vec2)] = paste(vec1[!is.na(vec1) & !is.na(vec2)], vec2[!is.na(vec1) & !is.na(vec2)], sep = "; ")
      }
    } else {
      vec12 = ifelse(!is.na(vec1), vec1, vec2)
    }
    initTab %<>%
      bind_cols(enframe(vec12, name = NULL, value = Col))
  }
  initTab %<>%
    select(!any_of(c(coln1, coln2)))
  initTab
}

## Auxiliary function for preparing a masked dataset according to 5 criteria (synonymous, tier 2, neutral, aS, Pool).
## Syn = TRUE masks synonymous variants; Tier2 = TRUE masks tier 2 variants; Neutral = TRUE masks neutral variants;
## aS = TRUE masks algorithmic S variants (classified S by the SOLO algorithm);
## Pool can either be NA (the default) or a name, in which case this masks any corresponding "poolable" variants.
## The final argument SOnly should be set to TRUE if only one iteration will be run, and to FALSE otherwise!
## When it is set, mutations that occur in S samples only will be marked rather than masked (i.e. filtered out).
## Setting it also means that samples containing at least one het will be removed from the masked dataset.
prepMask = function(inputTab, Syn = FALSE, Tier2 = TRUE, Neutral = TRUE, aS = FALSE, Pool = NA, SOnly = TRUE) {
  if (Syn) {
    inputTab %<>%
      dplyr::filter(effect == "missing" | !(effect %in% SILENT_EFFECTS))
  }
  if (Tier2) {
    inputTab %<>%
      dplyr::filter(tier != 2)
  }
  if (Neutral) {
    inputTab %<>%
      dplyr::filter(!neutral)
  }
  if (aS) {
    inputTab %<>%
      dplyr::filter(is.na(class) | class != "S")
  }
  if (!is.na(Pool)) {
    inputTab %<>%
      dplyr::filter(str_ends(variant, Pool) | effect == "missing" | !(effect %in% POOLED_EFFECTS[[Pool]]))
  }
  if (SOnly) {
    inputTab %<>%
      group_by(drug, variant) %>%
      mutate(SOnly = (sum(!het & phenotype == "S") == sum(!het))) %>%
      ungroup() %>%
      removeSamplesWithHets() %>%
      select(drug, sample_id, variant, phenotype, SOnly) %>%
      distinct()
  } else {
    inputTab %<>%
      select(drug, sample_id, variant, phenotype, het) %>%
      distinct()
  }
  inputTab
}

## Auxiliary function for identifying the samples that are not explained by the current variant classification
findUnexplained = function(inputTab) {
  inputTab %<>%
    group_by(sample_id, drug) %>%
    mutate(satisfying = all(phenotype == "S") || (all(is.na(class) | class == "S"))) %>%
    ungroup() %>%
    dplyr::filter(satisfying) %>%
    select(-satisfying)
  inputTab
}

## Auxiliary function for removing the samples that contain at least one het [eventually try to merge with the above]
removeSamplesWithHets = function(inputTab) {
  inputTab %<>%
    group_by(sample_id, drug) %>%
    mutate(satisfying = (all(!het))) %>%
    ungroup() %>%
    dplyr::filter(satisfying) %>%
    select(-satisfying)
  inputTab
}

## Auxiliary function for converting a vector into its logical form
convertToLogical = function(x) {
  output = ifelse(is.na(x), FALSE, as.logical(x))
  output
}

prepareOutputForPaolo = function(inputFile, skipBDQfromSA = FALSE, LoF = FALSE) {
  Tabs = vector("list", 2)
  Names = c("WHO", "ALL")
  inputTab = read_csv(inputFile, show_col_types = FALSE, guess_max = Inf) %>%
    rename(Datasets = datasets, Neutral_masked = neutral, OR_SOLO_FE_sig = OR_SOLO_pval_FDR_sig) %>%
    rename(PPV_conditional_SOLO = PPVc_SOLO, PPV_conditional_SOLO_lb = PPVc_SOLO_lb, PPV_conditional_SOLO_ub = PPVc_SOLO_ub) %>%
    rename(Absent_R = absent_R, Absent_S = absent_S, Present_R = present_R, Present_S = present_S, Present_SOLO_R = SOLO_R, Present_SOLO_SR = SOLO_SorR) %>%
    rename(Additionalgradingcriteria = `Additional grading criteria_v1`, Final_Confidence_Grading = Final, Initial_Confidence_Grading = Initial) %>%
    select(-class, -correctAll, -correctSOLO, -FDR_threshold, -OR_pvalue, -OR_pval_FDR_sig, -OR_pval_max, -OR_pval_rank, -OR_SOLO_pval_max, -stage) %>%
    select(-starts_with("set"))
  for (index in 1:2) {
    curTab = inputTab %>%
      filter(Datasets == Names[index])
    curColn = colnames(curTab)
    colnames(curTab) = paste0(Names[index], "_", curColn)
    Tabs[[index]] = curTab
  }
  finalTab = full_join(Tabs[[1]], Tabs[[2]], by = c("WHO_drug" = "ALL_drug", "WHO_variant" = "ALL_variant"))
  write_csv(finalTab, paste0("Merged_Graded_Stats_", Sys.Date(), ifelse(skipBDQfromSA, "_withoutSA", ""), ifelse(LoF, "_withLoFs", ""), ".csv"))
  finalTab
}
