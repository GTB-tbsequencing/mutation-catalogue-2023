 gradeMutations = function(LoF = TRUE, NON_DATABASE_DIRECTORY = NULL) {
  ## Prepare the master variant table, compute genes and mutations, and conduct a basic consistency check
  Tab0 = read_csv(paste0("Stats_WHO" , ifelse(LoF, "_withLoFs", ""), ".csv"), guess_max = Inf, show_col_types = FALSE)
  Tab1 = read_csv(paste0("Stats_MAIN", ifelse(LoF, "_withLoFs", ""), ".csv"), guess_max = Inf, show_col_types = FALSE)
  stopifnot(all(Tab0$datasets == "WHO"))
  stopifnot(all(Tab1$datasets == "ALL"))
  stopifnot(nrow(Tab0 %>% anti_join(Tab1, by = c("variant", "drug"))) == 0)
  inputTab = bind_rows(Tab0, Tab1) %>%
    mutate(gene = str_split_fixed(variant, "_", n = 2)[,1], mutation = str_split_fixed(variant, "_", n = 2)[,2])
  stopifnot(all(inputTab %>% dplyr::filter(gene ==  PZA_GENE) %>% pull(drug) == "Pyrazinamide"))
  stopifnot(all(inputTab %>% dplyr::filter(gene == RRDR_GENE) %>% pull(drug) == "Rifampicin"  ))
  ## Initial confidence grading: upgrade variants matching the basic upgrade rule to grade 1; apply rules for upgrading and downgrading pncA variants
  inputTab %<>%
    mutate(Initial            = ifelse(!is.na(neutral) & neutral, 5, 3)) %>% 
    mutate(Rule_Initial       = ifelse(!is.na(neutral), ifelse(datasets == "WHO", 1, 6), ifelse(datasets == "WHO", 5, 9))) %>%
    mutate(rule_initial_up    = (SOLO_SorR >=5     & !is.na(PPVc_SOLO_lb) & PPVc_SOLO_lb >= 0.25 & OR_SOLO > 1        & OR_SOLO_pval_FDR_sig)) %>%
    mutate_at("Initial",      ~{ifelse(rule_initial_up, 1, .)}) %>%
    mutate_at("Rule_Initial", ~{ifelse(rule_initial_up, ifelse(datasets == "WHO", 2, 7), .)}) %>%
    mutate(rule_pncA_down     = (gene == PZA_GENE  & !is.na(PPV_SOLO)     & PPV_SOLO < 0.4       & PPV_SOLO_ub < 0.75 & datasets == "WHO")) %>%
    mutate_at("Rule_Initial", ~{ifelse(rule_pncA_down, 3, .)}) %>%
    mutate(rule_pncA_up       = (gene == PZA_GENE  & SOLO_R >= 2          & PPV >= 0.5)) %>%
    mutate_at("Initial",      ~{ifelse(rule_pncA_down & Initial == 3, 4, ifelse(rule_pncA_up & Initial == 3, 2, .))}) %>%
    mutate_at("Rule_Initial", ~{ifelse(rule_pncA_up, ifelse(datasets == "WHO", 4, 8), .)})
  ## Reconcile initial confidence grading between the WHO and the ALL datasets by first separating them once again
  Tab0 = inputTab %>%
    filter(datasets == "WHO")
  inputTab %<>%
    filter(datasets == "ALL") %>%
    full_join(Tab0, by = c("variant", "drug", "gene", "mutation"), suffix = c("_ALL", "_WHO")) %>%
    mutate_at("Initial_WHO", ~{ifelse(is.na(.), 3, .)}) ## Not all variants are initially WHO-graded!
  inputTab %<>%
    mutate(    Initial    = ifelse(Initial_WHO == Initial_ALL         , Initial_WHO , NA)) %>%
    mutate(    datasets   = ifelse(Initial_WHO == Initial_ALL         , "ALL+WHO"   , NA)) %>%
    mutate_at("Initial" , ~{ifelse(Initial_WHO == 3 & Initial_ALL != 3, Initial_ALL , .)}) %>%
    mutate_at("datasets", ~{ifelse(Initial_WHO == 3 & Initial_ALL != 3, "ALL"       , .)}) %>%
    mutate_at("Initial" , ~{ifelse( is.na(Initial_ALL) | (Initial_ALL == 3 & Initial_WHO != 3),             Initial_WHO , .)}) %>%
    mutate_at("datasets", ~{ifelse( is.na(Initial_ALL) | (Initial_ALL == 3 & Initial_WHO != 3),             "WHO"       , .)}) %>%
    mutate_at("Initial" , ~{ifelse( pmax(Initial_WHO, Initial_ALL) <= 2 & Initial_WHO != Initial_ALL,       Initial_WHO , .)}) %>%
    mutate_at("datasets", ~{ifelse( pmax(Initial_WHO, Initial_ALL) <= 2 & Initial_WHO != Initial_ALL,       "WHO"       , .)}) %>%
    mutate_at("Initial" , ~{ifelse( Initial_WHO == 4 & Initial_ALL <= 2,                                    INITIAL_FLAG, .)}) %>%
    mutate_at("datasets", ~{ifelse( Initial_WHO == 4 & Initial_ALL <= 2,                                    "FLAG"      , .)})
  ## Extract additional grading criteria from v1, remove unused ones, then prepare to compute the final grades and record additional grading criteria
  ## Note that only one rule is applied per variant-drug pair; those to which a rule has been applied are marked by setting anyRule to TRUE throughout
  extraTab = read_csv("v1_grades.csv", guess_max = Inf, show_col_types = FALSE) %>%
    mutate(Final_v1 = match(final_grading_v1, GRADES))
  inputTab %<>%
    left_join(extraTab, by = c("drug", "variant")) %>%
    mutate_at("Additional grading criteria_v1", ~{ifelse(. %in% IGNORED_CRITERIA, NA, .)}) %>%
    mutate(Final = Initial, `Additional grading criteria` = `Additional grading criteria_v1`, Rule_Final = 15, anyRule = FALSE)
  ## Upgrade to grade 4 variants initially graded 5 by setC only (i.e. literature only) or version 1 guidance; specify which one applies
  inputTab %<>%
    mutate(rule_COnly       = (setC_ALL     & !(setA_ALL | setB_ALL | setD_ALL | setE_ALL)                & Initial == 5)) %>%
    applyExpertRule("rule_COnly",  description = SET_C_ONLY,  finalGrade = 4, finalRule = 16) %>%
    mutate(rule_V1Only      = (version1_ALL & !(setA_ALL | setB_ALL | setD_ALL | setE_ALL | setC_ALL)     & Initial == 5)) %>%
    applyExpertRule("rule_V1Only", description = V1_NEUTRALS, finalGrade = 4, finalRule = 16.5)
  ## Downgrade to grade 2 variants initially graded 1 by the ALL dataset only
  inputTab %<>%
    mutate(rule_AllOnly     = (Initial_ALL == 1 & Initial_WHO == 3)) %>%
    applyExpertRule("rule_AllOnly",  description = ALL_ONLY,  finalGrade = 2, finalRule = 17)
  ## Mark those variants which had discrepant AwR grades in the two datasets as "Based on WHO dataset"
  inputTab %<>%
    mutate(rule_WHOBased    = (pmax(Initial_WHO, Initial_ALL) <= 2 & Initial_WHO != Initial_ALL)) %>%
    applyExpertRule("rule_WHOBased",  description = WHO_BASED, finalGrade = NA, finalRule = 18)
  ## Mark those variants which triggered a flag at the initial classification as "Manual check required"
  inputTab %<>%
    mutate(rule_FLAG        = (Initial == INITIAL_FLAG)) %>%
    applyExpertRule("rule_FLAG",      description = MANUAL_CHECK, finalGrade = NA, finalRule = 19)
  ## Upgrade rpoB borderline mutations to grade 1; note: this rule overrides any previously applied rules!
  inputTab %<>%
    mutate(rule_rpoB_borderline = (gene == RRDR_GENE & mutation %in% BORDERLINE_MUTATIONS)) %>%
    applyExpertRule("rule_rpoB_borderline", description = BORDERLINE, finalGrade = 1, finalRule = 20, applyAlways = TRUE)
  ## Downgrade any silent variant to grade 4, then set to NA their SOLO counts and related statistics!
  inputTab %<>%
    mutate(rule_silent = (effect_ALL %in% SILENT_EFFECTS                                                                    & Initial == 3)) %>%
    applyExpertRule("rule_silent",          description = SILENT,   finalGrade = 4, finalRule = 21) %>%
    mutate_at(str_subset(colnames(.), "SOLO"), ~{ifelse(Rule_Final == 21, NA, .)})
  ## Upgrade any non-silent variant in the RRDR region to grade 2
  inputTab %<>%
    mutate(rule_RRDR = (gene == RRDR_GENE & (pos1_ALL %in% RRDR_INT | pos2_ALL %in% RRDR_INT) & !effect_ALL %in% SILENT_EFFECTS & Initial == 3)) %>%
    applyExpertRule("rule_RRDR",            description = RRDR,    finalGrade = 2, finalRule = 22)
  ## If a pooled LoF is graded 1 or 2, then any grade 3 LoF mutation in the same drug-gene combination is upgraded to grade 2
  inputTab %<>%
    group_by(gene, drug) %>%
    mutate(rule_LoF_candidate = (any(mutation == "LoF" & Final <= 2))) %>%
    ungroup() %>%
    mutate(rule_LoF = (rule_LoF_candidate & mutation != "LoF" & effect_ALL %in% POOLED_EFFECTS[["LoF"]]                     & Initial == 3)) %>%
    applyExpertRule("rule_LoF",             description = LoF_MUT, finalGrade = 2, finalRule = 23)
  ## Apply the cross-resistance rules the first time around
  inputTab %<>%
    applyCrossResistanceRules(iteration = 1)
  ## Upgrade to grade 2 any variant that is "recognized as DR marker" in a WHO-endorsed assay, as well as any variant with allelic exchange evidence
  assayTab = read_csv(paste0(str_remove(NON_DATABASE_DIRECTORY, "/$"), "/assay_mutations_24Apr2023.csv"), show_col_types = FALSE, guess_max = Inf) %>%
    mutate_at("drug", ~{DRUG_LIST[match(., SHORT_NAMES)]}) %>%
    mutate(assay = TRUE)
  inputTab %<>%
    full_join(assayTab,   by = c("drug", "variant")) %>%
    mutate(assay = convertToLogical(assay),   anyRule = convertToLogical(anyRule), rule_Assay     = (assay    & (is.na(Initial) | Initial == 3))) %>%
    applyExpertRule("rule_Assay",     description = PROBLEMATIC_CRITERION, finalGrade = 2, finalRule = 27)
  allelicTab = read_csv(paste0(str_remove(NON_DATABASE_DIRECTORY, "/$"), "/allelic_exchanges.csv"),       show_col_types = FALSE, guess_max = Inf) %>%
    mutate_at("variant", ~{str_replace(., "lof$", "LoF")}) %>%
    mutate(allelic = TRUE)
  inputTab %<>%
    full_join(allelicTab, by = c("drug", "variant")) %>%
    mutate(allelic = convertToLogical(allelic), anyRule = convertToLogical(anyRule), rule_Allelic = (allelic & (is.na(Initial) | Initial == 3))) %>%
    applyExpertRule("rule_Allelic",   description = ALLELIC_EXP          , finalGrade = 2, finalRule = 28) %>%
    select(-assay, -allelic)
  ## If any of the previous 5 rules, 24-28, has upgraded a pooled LOF variant to grade 2, we apply rule 23 one more time to the corresponding variants 
  inputTab %<>%
    group_by(gene, drug) %>%
    mutate(rule_LoF_re_candidate = any(mutation == "LoF" & Final <= 2 & Rule_Final %in% 24:28)) %>%
    ungroup() %>%
    mutate(rule_LoF_re = (rule_LoF_re_candidate & mutation != "LoF" & effect_ALL %in% POOLED_EFFECTS[["LoF"]]              & Initial == 3)) %>%
    applyExpertRule("rule_LoF_re",    description = LoF_MUT,               finalGrade = 2, finalRule = 29)
  ## Apply the cross-resistance rules the second time around
  inputTab %<>%
    applyCrossResistanceRules(iteration = 2)
  ## Upgrade a specific variant to grade 2 based on Literature evidence (PMID 32571824) 
  inputTab %<>%
    mutate(rule_PZA_lit = (variant == PZA_SPEC_VAR                                                                         & Initial == 3)) %>%
    applyExpertRule("rule_PZA_lit",      description = describePMIDs(32571824), finalGrade = 2, finalRule = 30)
  ## Upgrade to grade 2 based on previous guidance to distinguish between "evidence in version 1" and "guidance before version 1" (e.g. Miotto ERJ2017)
  inputTab %<>%
    mutate(rule_prev_guidance = (!is.na(`Additional grading criteria`) & `Additional grading criteria` %in% PRE_GUIDANCE   & Initial == 3)) %>%
    applyExpertRule("rule_prev_guidance", description = PRE_GUIDANCE, finalGrade = FINAL_FLAG, finalRule = 31) %>%
    mutate(rule_v1_evidence   = (!is.na(Final_v1) & Final_v1 < 3       & !effect_ALL                  %in% INFRAME_EFFECTS & Initial == 3)) %>%
    applyExpertRule("rule_v1_evidence",   description = V1_EVIDENCE , finalGrade = FINAL_FLAG, finalRule = 31.5)
  ## Downgrade to grade 4 the potentially inflated-PPV variants:
  inputTab %<>%
    mutate(rule_inflation = (drug == "Bedaquiline" & gene == BDQ_CFZ_GENE[1] & mutation %in% INFLATED_PPV_VARS)) %>%
    applyExpertRule("rule_inflation",    description = INFLATION, finalGrade = 2, finalRule = 32)
  ## Add comment column to a specified list of mutations or mutation categories, provided they were initially graded 3 and no other rule has applied:
  commentTab = read_csv(paste0(NON_DATABASE_DIRECTORY, "/Comments_August3.csv"), show_col_types = FALSE, guess_max = Inf) %>%
    mutate_all(~{trimws(., whitespace = "[\\h\\v]")}) %>%
    set_colnames(c("drug", "gene", "mutation", "comment")) %>%
    mutate_at("mutation", ~{ifelse(. == "any AwR", 1, ifelse(. == "any AwRI", 2, ifelse(. == "any nAwRI", 4, ifelse(. == "any nAwR", 5, .))))})
  commentCategoryTab = commentTab %>%
    filter(nchar(mutation) == 1) %>%
    mutate(Final = as.integer(mutation)) %>%
    select(-mutation)
  commentLoF = commentTab %>%
    filter(str_starts(mutation, "any LoF")) %>%
    select(-mutation)
  commentSingleTab = commentTab %>%
    filter(nchar(mutation) > 1 & !str_starts(mutation, "any LoF"))
  inputTab %<>%
    full_join(commentLoF,         by = c("drug", "gene")) %>%
    mutate_at("comment",    ~{ifelse(!(effect_ALL %in% POOLED_EFFECTS[["LoF"]]), NA,  .)}) %>%
    full_join(commentCategoryTab, by = c("drug", "gene", "Final"   ), suffix = c(".x", ".z")) %>%
    adjustDuplicateColumns(suffixes = c(".x", ".z"), add = TRUE) %>%
    full_join(commentSingleTab  , by = c("drug", "gene", "mutation"), suffix = c(".x", ".z")) %>%
    adjustDuplicateColumns(suffixes = c(".x", ".z"), add = TRUE)
  ## Specify PMIDs that lead to a downgrade by rule 16:
  inputTab %<>%
    mutate_at("Additional grading criteria", ~{ifelse(drug == "Bedaquiline"      & Rule_Final == 16, describePMIDs(c(28031270, 34503982)), .)}) %>% 
    mutate_at("Additional grading criteria", ~{ifelse(drug %in% FQS              & Rule_Final == 16, describePMIDs(28137812),              .)}) %>%
    mutate_at("Additional grading criteria", ~{ifelse(drug %in% FIRST_LINE_DRUGS & Rule_Final == 16, describePMIDs(32143680),              .)})
  ## Remove variant addition variables
  inputTab %<>%
    select(-starts_with("add_"))
  ## Lastly, convert the numerical grades to their description and save the file
  inputTab %<>%
    mutate(Initial = GRADES[Initial], Final = GRADES[Final])
  write_csv(inputTab, paste0(paste("Final_graded_algorithm_catalogue", Sys.Date(), "Leonid", sep = "_"), 
                             ifelse(LoF, "_withLoFs", ""), ".csv"))
  inputTab
}

describePMIDs = function(PMIDs) {
  output = paste0("Literature evidence (PMID ", paste(PMIDs, collapse = "; "), ")")
  output
}

applyCrossResistanceRules = function(inputTab, iteration = 1) {
  ## Upgrade to grade 2 based on FQ cross-resistance (LEV-MFX on gyrA and gyrB, bilateral); add grade 1/2 non-LoF variants appearing for only one drug
  extra_FQ_variants = inputTab %>%
    group_by(variant) %>%
    mutate(add_FQ_variant = (any(sum(drug %in% FQS) == 1 & gene %in% FQ_GENE & Final < 3 & !(effect_ALL %in% POOLED_EFFECTS[["LoF"]])))) %>%
    ungroup() %>%
    filter(add_FQ_variant & drug %in% FQS) %>%
    mutate(drug = FQS[3 - match(drug, FQS)], Initial = 3) %>%
    mutate_at(as.vector(outer(outer(c("present", "absent", "SOLO"), c("_R", "_S"), paste0), c("_ALL", "_WHO"), paste0)), ~{ NA })
  inputTab %<>%
    bind_rows(extra_FQ_variants) %>%
    group_by(variant) %>%
    mutate(rule_FQ_cross_res = (gene %in% FQ_GENE      & any(drug %in% FQS & Final < 3) & drug %in% FQS                  & Initial == 3)) %>%
    applyExpertRule("rule_FQ_cross_res",     description = FQ_CROSS_RES, finalGrade = 2, finalRule = 24 + (iteration - 1)/2) %>%
    ungroup()
  ## Upgrade to grade 2 based on INH-ETH cross-resistance (INH-ETH on inhA and fabG1, bilateral); first, add non-LoF variants appearing for only one drug
  extra_IE_variants = inputTab %>%
    group_by(variant) %>%
    mutate(add_IE_variant = (any(sum(drug %in% INH_ETH) == 1 & gene %in% INH_ETH_GENE & Final < 3 & !(effect_ALL %in% POOLED_EFFECTS[["LoF"]])))) %>%
    ungroup() %>%
    filter(add_IE_variant & drug %in% INH_ETH) %>%
    mutate(drug = INH_ETH[3 - match(drug, INH_ETH)], Initial = 3) %>%
    mutate_at(as.vector(outer(outer(c("present", "absent", "SOLO"), c("_R", "_S"), paste0), c("_ALL", "_WHO"), paste0)), ~{ NA })
  inputTab %<>%
    bind_rows(extra_IE_variants) %>%
    group_by(variant) %>%
    mutate(rule_IE_cross_res = (gene %in% INH_ETH_GENE & any(drug %in% INH_ETH & Final < 3) & drug %in% INH_ETH          & Initial == 3)) %>%
    applyExpertRule("rule_IE_cross_res",      description = IE_CROSS_RES, finalGrade = 2, finalRule = 25 + (iteration - 1)/2) %>%
    ungroup()
  ## Upgrade to grade 2 based on BDQ-CFZ cross-resistance (BDQ-CFZ on Rv0678 and pepQ, BDQ over CFZ only)
  extra_BC_variants = inputTab %>%
    group_by(variant) %>%
    mutate(add_BC_variant = (any(drug == BDQ_CFZ[1]) & all(drug != BDQ_CFZ[2]) & gene %in% BDQ_CFZ_GENE & Final < 3 & !(effect_ALL %in% POOLED_EFFECTS[["LoF"]]))) %>%
    ungroup() %>%
    filter(add_BC_variant & drug == BDQ_CFZ[1]) %>%
    mutate(drug = BDQ_CFZ[2], Initial = 3) %>%
    mutate_at(as.vector(outer(outer(c("present", "absent", "SOLO"), c("_R", "_S"), paste0), c("_ALL", "_WHO"), paste0)), ~{ NA })
  inputTab %<>%
    bind_rows(extra_BC_variants) %>%
    group_by(variant) %>%
    mutate(rule_BC_cross_res = (gene %in% BDQ_CFZ_GENE & any(drug == BDQ_CFZ[1] & Final < 3) & drug == BDQ_CFZ[2]        & Initial == 3)) %>%
    applyExpertRule("rule_BC_cross_res",      description = BC_CROSS_RES, finalGrade = 2, finalRule = 26 + (iteration - 1)/2) %>%
    ungroup()
  inputTab
}

applyExpertRule = function(inputTab, ruleColumn, description = NA, finalGrade = NA, finalRule = NA, applyAlways = FALSE) {
  if (applyAlways) {
    inputTab %<>%
      mutate(applyRule = .data[[ruleColumn]])
  } else {
    inputTab %<>%
      mutate(applyRule = (!anyRule & .data[[ruleColumn]]))
  }
  if (!is.na(description)) {
    inputTab %<>%
      mutate_at("Additional grading criteria", ~{ifelse(applyRule, description, .)})
  }
  if (!is.na(finalGrade)) {
    inputTab %<>%
      mutate_at("Final",                       ~{ifelse(applyRule, finalGrade, .)})
  }
  if (!is.na(finalRule)) {
    inputTab %<>%
      mutate_at("Rule_Final",                  ~{ifelse(applyRule, finalRule, .)})
  }
  inputTab %<>%
    mutate_at("anyRule",                       ~{or(., applyRule)}) %>%
    select(-applyRule)
  inputTab
}
