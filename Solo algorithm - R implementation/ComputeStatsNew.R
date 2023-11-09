safeBinomTest = function(x, y) {
  if (y == 0) { 
    output = list(estimate = NA, conf.int = c(NA, NA))
  } else {
    output = binom.test(x, y, p = 0.5, alternative = "two.sided", conf.level = 1 - SIG_THRESHOLD)
  }
  output
}

safeFisherTest = function(x, y, z, w) {
  if (x + y == 0 || x + z == 0 || y + w == 0 || z + w == 0) {
    output = list(p.value = NA, conf.int = c(NA, NA))
  } else {
    output = fisher.test(rbind(c(x, y), c(z, w)), alternative = "two.sided", conf.level = 1 - SIG_THRESHOLD)
  }
  output
}

computeCatalogueStats = function(TabT, correct_all = TRUE) {
  TabT %<>%
    dplyr::filter(SOLO_R + SOLO_S + absent_R + absent_S > 0)
  print(paste("There are", nrow(TabT), "rows to process"))
  PPVTab = TabT %>%
    distinct(present_R, present_S, .keep_all = FALSE) %>%
    mutate(PPV = map2(present_R, present_R + present_S, safeBinomTest))
  PPV_SOLOTab = TabT %>%
    distinct(SOLO_R, SOLO_S, .keep_all = FALSE) %>%
    mutate(PPV_SOLO = map2(   SOLO_R,    SOLO_R + SOLO_S   , safeBinomTest))
  PPVc_SOLOTab = TabT %>%
    distinct(SOLO_R, present_S, .keep_all = FALSE) %>%
    mutate(PPVc_SOLO = map2(   SOLO_R,    SOLO_R + present_S, safeBinomTest))
  SensTab = TabT %>%
    distinct(present_R, absent_R, .keep_all = FALSE) %>%
    mutate(Sens = map2(present_R, present_R + absent_R , safeBinomTest))
  Sens_SOLOTab = TabT %>%
    distinct(SOLO_R, absent_R, .keep_all = FALSE) %>%
    mutate(Sens_SOLO = map2(   SOLO_R,    SOLO_R + absent_R , safeBinomTest))
  SpecTab = TabT %>%
    distinct(present_S, absent_S, .keep_all = FALSE) %>%
    mutate(Spec = map2( absent_S, present_S + absent_S , safeBinomTest))
  TabT %<>%
    left_join(PPVTab,       by = c("present_R", "present_S")) %>%
    left_join(PPV_SOLOTab,  by = c("SOLO_R"   , "SOLO_S"   )) %>%
    left_join(PPVc_SOLOTab, by = c("SOLO_R"   , "present_S")) %>%
    left_join(SensTab,      by = c("present_R", "absent_R" )) %>%
    left_join(Sens_SOLOTab, by = c("SOLO_R"   , "absent_R" )) %>%
    left_join(SpecTab,      by = c("present_S", "absent_S" )) %>%
    mutate(Spec_SOLO = Spec)
  TabT %<>%
    mutate(across(where(is.list), list(lb = ~{map_dbl(., ~{.$conf.int[1]})}, ub = ~{map_dbl(., ~{.$conf.int[2]})}))) %>%
    mutate(across(where(is.list), ~{map_dbl(., ~{.$estimate})}))
  ORTab = TabT %>%
    distinct(                    present_R,   present_S,   absent_R,   absent_S,   .keep_all = FALSE) %>%
    mutate(OR        = pmap(list(present_R,   present_S,   absent_R,   absent_S),  safeFisherTest))
  ORSoloTab = TabT %>%
    distinct(                    SOLO_R,      SOLO_S,      absent_R,   absent_S,   .keep_all = FALSE) %>%
    mutate(OR_SOLO   = pmap(list(SOLO_R,      SOLO_S,      absent_R,   absent_S),  safeFisherTest))
  TabT %<>%
    left_join(ORTab    , by = c("present_R", "present_S", "absent_R", "absent_S")) %>%
    left_join(ORSoloTab, by = c("SOLO_R",    "SOLO_S",    "absent_R", "absent_S")) %>%
    mutate(across(where(is.list), list(exact_lb = ~{map_dbl(., ~{.$conf.int[1]})}, exact_ub = ~{map_dbl(., ~{.$conf.int[2]})}, pvalue = ~{map_dbl(., ~{.$p.value})})))
  TabT %<>%
    mutate(OR = (absent_S * present_R)/(absent_R * present_S), OR_SOLO = (absent_S * SOLO_R)/(absent_R * SOLO_S), FDR_threshold  = FDR_THRESHOLD) %>%
    mutate(OR_HCI = sqrt(1/absent_R + 1/absent_S + 1/present_R + 1/present_S), OR_SOLO_HCI = sqrt(1/absent_R + 1/absent_S + 1/SOLO_R + 1/SOLO_S)) %>%
    mutate(OR_lb = exp(log(OR) - CI_COEFFICIENT * OR_HCI), OR_SOLO_lb = exp(log(OR_SOLO) - CI_COEFFICIENT * OR_SOLO_HCI)) %>%
    mutate(OR_ub = exp(log(OR) + CI_COEFFICIENT * OR_HCI), OR_SOLO_ub = exp(log(OR_SOLO) + CI_COEFFICIENT * OR_SOLO_HCI)) %>%
    mutate(correct = ifelse(rep(correct_all, nrow(.)), correctAll & !is.na(OR), correctSOLO & !is.na(OR_SOLO)), OR_HCI = NULL, OR_SOLO_HCI = NULL)
  TabToCorrect = TabT %>%
    dplyr::filter(correct) %>%
    group_by(drug) %>%
    mutate(k = n()) %>%
    mutate(OR_pval_rank      = rank(OR_pvalue)     , OR_pval_max      = max(OR_pval_rank     [OR_pvalue     / OR_pval_rank      <= FDR_threshold/k], na.rm = TRUE)) %>%
    mutate(OR_SOLO_pval_rank = rank(OR_SOLO_pvalue), OR_SOLO_pval_max = max(OR_SOLO_pval_rank[OR_SOLO_pvalue/ OR_SOLO_pval_rank <= FDR_threshold/k], na.rm = TRUE)) %>%
    mutate(OR_pval_FDR_sig = (OR_pval_rank <= OR_pval_max), OR_SOLO_pval_FDR_sig = (OR_SOLO_pval_rank <= OR_SOLO_pval_max)) %>%
    ungroup()
  TabNotToCorrect = TabT %>%
    dplyr::filter(!correct | is.na(correct)) %>%
    mutate(OR_pval_FDR_sig = FALSE, OR_SOLO_pval_FDR_sig = FALSE)
  TabT = bind_rows(TabToCorrect, TabNotToCorrect) %>%
    select(-correct)
  TabT
}

## Original version of the FDR correction: a more explicit alternative
## mutate(OR_p_rank      = rank(OR_p)     , OR_p_max      = max(OR_p_rank     [OR_p     /OR_p_rank      <= FDR_threshold/k])) %>%
## mutate(OR_SOLO_p_rank = rank(OR_SOLO_p), OR_SOLO_p_max = max(OR_SOLO_p_rank[OR_SOLO_p/OR_SOLO_p_rank <= FDR_threshold/k])) %>%
## mutate(     OR_pvalue_threshold =      OR_pvalue[     OR_pvalue_rank ==      OR_pvalue_max][1]) %>%
## mutate(OR_SOLO_pvalue_threshold = OR_SOLO_pvalue[OR_SOLO_pvalue_rank == OR_SOLO_pvalue_max][1]) %>%
