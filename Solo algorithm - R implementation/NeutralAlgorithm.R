## This function executes the neutral algorithm, which identifies variants that should be considered as neutral
## If safe = FALSE, existing conversion files for mapping version 1 to version 2 are used directly (no checks!)
## The last option is only used for naming the (intermediate) output files, and should not normally be set
neutralAlgorithm = function(masterTab, safe = FALSE, NON_DATABASE_DIRECTORY="STATA DTA files", DATA_DIRECTORY="DATABASE EXTRACTION files", EXTRACTION_ID="2023-04-18T08_34_28.150000_jr_8b05f8aa02accfa629d731e456c23f728227d785f2352b09c1f8f098eab0d6c3") {
  for (set in c("A", "B")) {
    if (set == "B") {
      ## Goal: remove the sample-drug pairs with a URM (unconditional resistance mutation; first, import the list
      urmTab = read_csv(paste(NON_DATABASE_DIRECTORY, "urm.csv", sep="/"), show_col_types = FALSE, guess_max = Inf)
      ## Merge with the master table
      masterTab %<>%
        left_join(urmTab, by = c("drug", "variant"))
      ## Computing all the rows which are to be marked with urm = TRUE; first, convert urm to a logical variable
      masterTab %<>%
        mutate_at("urm", convertToLogical)
      ## Non-silent variants in the RRDR_GENE with exactly one of the two AAs in the RRDR are also URMs
      stopifnot(all(masterTab %>% dplyr::filter(gene == RRDR_GENE) %>% pull(drug) == "Rifampicin"))
      masterTab %<>%
        mutate(RRDR_NON_SILENT = (gene == RRDR_GENE & (pos1 %in% RRDR_INT | pos2 %in% RRDR_INT) & !effect %in% SILENT_EFFECTS))
      masterTab %<>%
        mutate(urm = or(urm, RRDR_NON_SILENT))
      ## Variants in the RRDR_GENE that are borderline mutations are also URMs; new: check they are in the URM file
      stopifnot(all(paste0(RRDR_GENE, "_", BORDERLINE_MUTATIONS) %in% (urmTab %>% filter(drug == "Rifampicin") %>% pull(variant))))
      ## Variants in tier 1 SPECIAL_URM_GENES that have an LoF mutation are also URMs
      masterTab %<>%
        mutate(urm = or(urm, gene %in% SPECIAL_URM_GENES & effect %in% POOLED_EFFECTS[["LoF"]] & tier == 1))
      ## Variants in BDQ_URM_GENE that have an effect among the BDQ_URM_EFFECTS are also URMs
      masterTab %<>%
        mutate(urm = or(urm, gene == BDQ_GENE & effect %in% BDQ_EFFECTS))
      ## Mark any sample-drug pair with at least one of the variants being an unconditional resistance mutation
      masterTab %<>%
        group_by(sample_id, drug) %>%
        mutate(urm_in_sample = any(urm)) %>%
        ungroup()
    }
    ## Compute the PPVs
    masterTab %<>%
      computePPVs(removeURM = (set == "B"))
    ## Place the variant-drug combinations with PPV_ub < threshold into the set
    masterTab %<>%
      mutate(set = (!is.na(variant) & !(variant == "missing") & (!is.na(PPV_ub)) & (PPV_ub < PPV_UB_THRESHOLD)))
    ## Change the set name according to its letter
    colnames(masterTab)[ncol(masterTab)] = paste0("set", set)
    ## Extract the neutral variants into a separate file
    curNeutral = extractNeutral(masterTab, set)
    ## Remove generic auxiliary variables from the master table
    masterTab %<>%
      select(-PPV, -PPV_ub)
  }
  ## Goal: complete set C by adding literature mutations; first, import the list of literature mutations
  litTab = read_csv(paste(NON_DATABASE_DIRECTORY, "merker_neutrals.csv", sep="/"), show_col_types = FALSE, guess_max = Inf)
  ## Merge with the master table
  masterTab %<>%
    full_join(litTab, by = c("drug", "variant"))
  ## Define and compute lit_mutation as a logical form of merker and remove merker
  masterTab %<>%
    mutate(lit_mutation = convertToLogical(merker)) %>%
    select(-merker)
  ## Define set C as the union of sets A, B and lit_mutations
  masterTab %<>%
    mutate(setC = (setA | setB | lit_mutation))
  ## Extract the neutral variants into a separate file
  curNeutral = extractNeutral(masterTab, "C")
  for (set in c("D", "E")) {
    ## filter out silent, tier 2 variants, and those in previous sets (C, D?), plus isolates containing a het
    subsetTab = masterTab %>%
      dplyr::filter(!(effect %in% SILENT_EFFECTS | tier == 2 | setC))
    if (set == "E") {
      subsetTab = subsetTab %>%
        dplyr::filter(!setD)
    }
    subsetTab %<>%
      removeSamplesWithHets()
    ## Compute the PPV_solo for this subset
    subsetTab %<>%
      computePPVs(removeURM = FALSE, solo = TRUE, restrict = TRUE)
    ## Add the computed values back to the master table 
    masterTab %<>%
      left_join(subsetTab, by = c("drug", "variant"))
    ## Place the variant-drug combinations with PPV_ub < threshold into the set
    masterTab %<>%
      mutate(set = ((!is.na(PPV_ub)) & (PPV_ub < PPV_UB_THRESHOLD)))
    ## Change the set name according to its letter
    colnames(masterTab)[ncol(masterTab)] = paste0("set", set)
    ## Extract the neutral variants into a separate file
    curNeutral = extractNeutral(masterTab, set)
    ## Remove generic auxiliary variables from the master table
    masterTab %<>%
      select(-PPV, -PPV_ub)
  }
  ## Import version 1 mutations, convert them to version 2, then mark as neutral
  convertV1ToV2(safe = safe, NON_DATABASE_DIRECTORY=NON_DATABASE_DIRECTORY, DATA_DIRECTORY=DATA_DIRECTORY, EXTRACTION_ID=EXTRACTION_ID)
  v1Tab = read_csv("neutral_mutations_catalogue_v1.csv", show_col_types = FALSE, guess_max = Inf) %>%
    mutate(version1 = TRUE) %>%
    select(-variant_category_v1)
  ## Merge with the master table
  masterTab %<>%
    full_join(v1Tab, by = c("drug", "variant"))
  ## Convert the variable to a logical one
  masterTab %<>%
    mutate_at(c("version1", paste0("set", LETTERS[1:5]), "lit_mutation"), convertToLogical)
  ## Create the overall list of mutations
  masterTab %<>%
    mutate(setF = (setC | setD | setE | version1))
  ## Extract the neutral variants into a separate file
  curNeutral = extractNeutral(masterTab, "F")
  ## Extract the complete collection of marks for each variant-drug pair
  allSets = masterTab %>%
    group_by(drug, variant) %>%
    slice(1) %>%
    ungroup() %>%
    select(drug, variant, setA, setB, setC, lit_mutation, setD, setE, version1, setF)
  ## Record the complete collection of marks
  masterTab
}

## This function computes the PPV for each variant-drug pair of a given input table
## If removeURM = TRUE, the table is first filtered down to samples without any URM
## If solo = TRUE, the table is first filtered down to only variants that are solos
## If restrict = TRUE, returns the filtered table with variants, drugs and PPV stats
computePPVs = function(inputTab, removeURM = TRUE, solo = FALSE, restrict = FALSE) {
  auxTab = inputTab
  if (removeURM && "urm_in_sample" %in% colnames(auxTab)) {
    auxTab %<>%
      dplyr::filter(!urm_in_sample)
  }
  if (solo) {
    ## Only keep the mutations that are solos within their respective sample-drug pair
    auxTab %<>%
      group_by(sample_id, drug) %>%
      mutate(N = n()) %>%
      ungroup() %>%
      mutate(solo = (N == 1)) %>%
      dplyr::filter(solo) %>%
      select(-solo, -N)
  }
  ## Compute the PPV and its confidence intervals for each variant-drug pair
  auxTab %<>%
    group_by(drug, variant) %>%
    mutate(total = n(), totalR = sum(phenotype == "R")) %>%
    slice(1) %>%
    ungroup
  PPVTab = auxTab %>%
    distinct(totalR, total, .keep_all = FALSE) %>%
    mutate(testResult = map2(totalR, total, safeBinomTest))
  auxTab %<>%
    left_join(PPVTab, by = c("totalR", "total")) %>%
    mutate(PPV = map_dbl(testResult, ~{.$estimate}), PPV_ub = map_dbl(testResult, ~{.$conf.int[2]})) %>%
    select(drug, variant, PPV, PPV_ub)
  if (restrict) {
    return(auxTab)
  } else { ## Incorporate the computed values into the original input
    inputTab %<>%
      left_join(auxTab, by = c("drug", "variant"))
    return(inputTab)
  }
}

extractNeutral = function(inputTab, setName) {
  ## Create the full variable name
  varName = paste0("set", setName)
  ## Extract the corresponding variable
  varMask = inputTab %>%
    select(one_of(varName)) %>%
    pull()
  ## Extract the variant-drug combinations in the set into a separate table
  neutralTab = inputTab %>%
    dplyr::filter(varMask) %>%
    group_by(drug, variant) %>%
    slice(1) %>%
    ungroup()
  ## Special case for set F: also keep track of all the other sets
  if (setName == "F") {
    neutralTab %<>%
      select(drug, variant, setA, setB, setC, setD, setE, version1, lit_mutation)
  } else {
    neutralTab %<>%
      select(drug, variant)
  }
  ## Save the set into a file
  outFilename = paste0("neutral_mutations_WHO_", setName, ".csv")
  write_csv(neutralTab, outFilename)
  neutralTab
}

convertV1ToV2 = function(safe = FALSE, NON_DATABASE_DIRECTORY="STATA DTA files", DATA_DIRECTORY="DATABASE EXTRACTION files", EXTRACTION_ID="2023-04-18T08_34_28.150000_jr_8b05f8aa02accfa629d731e456c23f728227d785f2352b09c1f8f098eab0d6c3") {
  ## Check if the files already exist
  if (!safe & all(file.exists(CONVERTED_FILES))) {
    print("All the required files have already been created; exiting!")
    return()
  }
  ## Save the starting directory to move back there at the end and go into the matching directory
  initDir = getwd()
  setwd(paste(DATA_DIRECTORY, EXTRACTION_ID, "v1_matching/", sep="/"))
  ## Preprocess the first file by renaming its columns and creating a new one for compatibility with version 2 
  matchTab = read_csv(list.files()[1], guess_max = Inf, show_col_types = FALSE) %>%
    rename(variant_category_v1 = 'description', gene = 'resolved_symbol', mutation = 'variant_category', effect = 'predicted_effect') %>%
    mutate(variant_category_v2 = paste(gene, mutation, sep = "_")) %>%
    select(-gene, -mutation)
  ## Keep only the non-silent mutations
  matchNonSilentTab = matchTab %>%
    dplyr::filter(!(effect %in% SILENT_EFFECTS))
  ## Extract revised list of neutral mutations
  oldTab = read_csv_2headers(paste0(NON_DATABASE_DIRECTORY, "/List_of_neutral_mutations_from_cat_ver_1_rev20Jan2023.csv"), prefix = FALSE)
  ## Rename the drugs according to their full names
  oldTab %<>%
    mutate_at("drug", ~{DRUG_LIST[match(., SHORT_NAMES)]})
  ## Remove any entries that were marked as not usable
  oldTab %<>%
    dplyr::filter(!str_starts(`FOR WHO CATALOGUE UPDATE VER. 2`, "NOT")) %>%
    select(drug, variant) %>%
    rename(variant_category_v1 = "variant")
  ## Find the common substrate between those and the matched non-silent variants
  oldTab %<>%
    inner_join(matchNonSilentTab, by = "variant_category_v1") %>%
    rename(variant = "variant_category_v2")
  ## Read in the final report from the publication (NB: the csv file contains only the Mutation_catalogue sheet)
  #### altTab = read_csv_2headers("../../../UNITAID_FIND_WHO_v1.21b_FINAL_21Apr2021.csv", prefix = TRUE) - this is to be used later on!
  altTab = read_csv_2headers(paste0(NON_DATABASE_DIRECTORY, "/WHO-UCN-GTB-PCI-2021.7-eng.csv"), prefix = FALSE) %>%
    mutate(across(starts_with("Present") | starts_with("Absent"), as.integer))
  ## Remove combo graded entries
  altTab %<>%
    dplyr::filter(`FINAL CONFIDENCE GRADING` != "combo") %>%
    rename(variant = `variant (common_name)`)
  ## Rename the drugs according to their full names
  altTab %<>%
    mutate_at("drug", ~{DRUG_LIST[match(., SHORT_NAMES)]})
  ## Remove the Genome position column, normalize variant names, remove empty ones, and rename some columns
  altTab %<>%
    select(-`Genome position`) %>%
    mutate_at("variant", ~{str_split_fixed(., " ", n = 2)[,1]}) %>%
    dplyr::filter(variant != "") %>%
    rename(variant_category_v1 = "variant", final_grading_v1 = `FINAL CONFIDENCE GRADING`) %>%
    rename(`Additional grading criteria_v1` = `Additional grading criteria`)
  ## Replace any problematic grading criteria with an empty string
  altTab %<>%
    mutate_at("Additional grading criteria_v1", ~{na_if(., PROBLEMATIC_CRITERION)})
  ## Check that all of the matched variants are present in the original table
  stopifnot(all(matchNonSilentTab$variant_category_v1 %in% altTab$variant_category_v1))
  ## Merge the two tables; variant_category_v2 is NA for the original variants that do not have a matched variant
  altTab %<>%
    full_join(matchNonSilentTab, by = "variant_category_v1")
  ## Among the entries that do not have a matched variant, exclude all except those that we want to postprocess
  ## First, drop variants in the fprA gene as it is no longer on the list
  altTab %<>%
    dplyr::filter(!(is.na(variant_category_v2) & str_starts(variant_category_v1, "fprA")))
  ## Mark the variants that are neither insertions nor deletions
  altTab %<>%
    mutate(non_indel = !str_detect(variant_category_v1, "del") & !str_detect(variant_category_v1, "ins"))
  ## Drop the variants in certain genes that are neither insertions nor deletions, and have no matched variant
  for (gene in c("rpsL", "embR", "rrs", "embC", "ubiA")) {
    altTab %<>%
      dplyr::filter(!(is.na(variant_category_v2) & str_starts(variant_category_v1, gene) & non_indel))
  }
  ## Drop all other promoter and indel variants
  altTab %<>%
    dplyr::filter(!(is.na(variant_category_v2) & (!non_indel | str_detect(variant_category_v1, "-"))))
  ## Remove the non_indel column
  altTab %<>%
    select(-non_indel)
  ## Drop all the uncertain variants
  altTab %<>%
    dplyr::filter(!(is.na(variant_category_v2) & str_detect(variant_category_v1, "[MV]1[A-Z]")))
  ## Drop the gid_V110G and inhA_T4I variants
  altTab %<>%
    dplyr::filter(!(is.na(variant_category_v2) & variant_category_v1 %in% c("gid_V110G", "inhA_T4I")))
  ## Manually recode the remaining pncA mutations
  xTab = tibble(variant_category_v1 = paste0("pncA_"  , c("H71D",     "I31T",     "L116R",     "T135N")), 
                variant_category_v2 = paste0("pncA_p.", c("His71Asp", "Ile31Thr", "Leu116Arg", "Thr135Asn")))
  altTab %<>%
    recodeValues(xTab)
  ## Double-check that every entry is annotated with a v2 variant now
  stopifnot(all(!is.na(altTab$variant_category_v2)))
  ## Keep a single representative per drug-version 2 variant combination
  ## New addition: remove any variants in v2 that correspond to more than one variant in v1
  altTab %<>%
    group_by(drug, variant_category_v2) %>%
    mutate(N = n()) %>%
    filter(N == 1) %>%
    ungroup() %>%
    select(-N)
  ## Recode the final grading values for certain specific variants
  yTab = tibble(variant_category_v2 = c("ethA_p.Met1?", "katG_p.Met1?", "pncA_p.Met1?"), final_grading_v1 = GRADES[c(1, 3, 1)])
  altTab %<>%
    recodeValues(yTab)
  ## Save a version with old and new variants and v1 confidence grading
  matchingTab = altTab %>% 
    select(drug, variant_category_v1, variant_category_v2, final_grading_v1, `Additional grading criteria_v1`)
  ## Then extract additional data from the corresponding additional dataset
  extraTab = read_csv(paste0(NON_DATABASE_DIRECTORY, "/assay_mutations_24Apr2023.csv"), guess_max = Inf, show_col_types = FALSE) %>%
    mutate_at("drug", ~{DRUG_LIST[match(., SHORT_NAMES)]}) %>%
    rename("variant_category_v2" = "variant") %>%
    mutate(`Additional grading criteria_v2` = PROBLEMATIC_CRITERION)
  matchingTab %<>%
    full_join(extraTab, by = c("variant_category_v2", "drug")) %>%
    mutate_at("Additional grading criteria_v1", ~{ifelse(`Additional grading criteria_v2` == PROBLEMATIC_CRITERION, 
                                                         PROBLEMATIC_CRITERION, .)})
  ## Save another version with new variants and additional grading criteria
  gradingTab = altTab %>% 
    rename(variant = "variant_category_v2") %>% 
    select(drug, `Additional grading criteria_v1`, final_grading_v1, variant)
  ## Go back to the starting directory and write the files
  setwd(initDir)
  ## Save the results into appropriate files
  write_csv(oldTab,      file = CONVERTED_FILES[1])
  write_csv(matchingTab, file = CONVERTED_FILES[2])
  write_csv(gradingTab,  file = CONVERTED_FILES[3])
}
