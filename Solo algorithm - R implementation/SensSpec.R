source("SOLOport/Driver.R")
source("SOLOport/ComputeStatsNew.R")
source("SOLOport/Constants.R")
source("SOLOport/Utilities.R")

computeSensSpec = function(version = 2, relaxed = FALSE, safe = TRUE, skipEpistasis = TRUE, sameRIF = TRUE,
                           DATA_DIRECTORY = "SOLO Algorithm Input files/DATABASE EXTRACTION files/",
                           EXTRACTION_ID = "2023-04-25T06_00_10.443990_jr_b741dc136e079fa8583604a4915c0dc751724ae9880f06e7c2eacc939e086536") {
  ## Extract the graded variants and the full dataset, and merge them
  gradedTab = read_csv(ifelse(version == 2, paste0("List_of_graded_variants_", ifelse(relaxed, "BDQ-CFZ-INH-DLM-relaxed_", ""), "25Apr2023.csv"), 
                                            paste0("Results/", EXTRACTION_ID, "/v1_grades.csv")), guess_max = Inf, show_col_types = FALSE)
  if (version == 2 && relaxed) {
    gradedTab %<>% 
      rename(Final_Confidence_Grading = `FINAL CONFIDENCE GRADING`)
  }
  if (version == 1) {
    gradedTab %<>% 
      rename(Final_Confidence_Grading = final_grading_v1)
  }
  gradedTab %<>% 
    mutate(Final = match(Final_Confidence_Grading, GRADES))
  fullDataset = read_csv(paste0("Results/", EXTRACTION_ID, "/CompleteDataset.csv"), guess_max = Inf, show_col_types = FALSE)
  initDir = getwd()
  setwd(paste(str_remove(DATA_DIRECTORY, "/$"), str_remove(EXTRACTION_ID, "/$"), "orphan_genotypes/", sep = "/"))
  orphanData = read_csv(list.files()[1], guess_max = Inf, show_col_types = FALSE) %>%
    rename(drug = drug_name, gene = resolved_symbol, mutation = variant_category, effect = predicted_effect) %>%
    mutate(variant = paste0(gene, "_", mutation), het = (`max(af)` < MAF_THRESHOLD_REGULAR), phenotype = "U") %>%
    mutate(allPos = map(mutation, ~{unlist(str_extract_all(., "[-]?[0-9]+"))})) %>%
    mutate(pos1 = map_int(allPos, ~{as.integer(first(.))}), pos2 = map_int(allPos, ~{as.integer(nth(., 2))}))
  setwd(initDir)
  fullDataset %<>%
    bind_rows(orphanData) %>%
    mutate(het_relaxed = (`max(af)` < MAF_THRESHOLD_RELAXED | variant == "missing"), het_strict = (`max(af)` < MAF_THRESHOLD_STRICT | variant == "missing")) %>%
    select(sample_id, drug, variant, phenotype, gene, mutation, effect, pos1, pos2, het, het_relaxed, het_strict)
  fullDataset %<>%
    left_join(gradedTab, by = c("drug", "variant")) %>%
    mutate_at(c("het_strict", "het_relaxed"), ~{replace_na(., TRUE)})
  ## Calculate the RIF resistance variable
  fullDataset %<>%
    mutate(drug_short = SHORT_NAMES[match(drug, DRUG_LIST)]) %>%
    mutate(RRDR_NON_SILENT  = (drug_short == "RIF" & gene == RRDR_GENE & (pos1 %in% RRDR_INT | pos2 %in% RRDR_INT) & !effect %in% SILENT_EFFECTS)) %>%
    group_by(sample_id) %>%
    mutate(genoRIF          = (drug_short == "RIF" & ((!is.na(Final) & Final <= 2) | RRDR_NON_SILENT))) %>%
    mutate(RIF_geno_Relaxed = any(genoRIF & !het_relaxed))
  if (sameRIF) {
    fullDataset %<>%
      mutate(RIF_geno_Regular = RIF_geno_Relaxed)
  } else {
    fullDataset %<>%
      mutate(RIF_geno_Regular = any(genoRIF & !het))
  }
  fullDataset %<>% ungroup()
  ## Implement the necessary expert rules; a previously ungraded mutation with either a non-silent RRDR variant or an LOF candidate gets grade 2
  if (version == 2) {
    fullDataset %<>%
      mutate(LOF_candidate = (effect %in% POOLED_EFFECTS[["LoF"]] & 
                                (gene %in% unlist(EXTENDED_ADD_GENES[-1]) | (drug_short == "DLM" & gene %in% EXTENDED_ADD_GENES[["DLM"]]))))
  } else {
    fullDataset %<>%
      mutate(LOF_candidate = ((drug_short == "CAP" & gene == CAP_GENE & effect %in% setdiff(POOLED_EFFECTS[["LoF"]], "start_lost")) |
                                (gene %in% unlist(SUB_EXTENDED_GENES) & effect %in% c(POOLED_EFFECTS[["LoF"]], INFRAME_EFFECTS))))
  }
  fullDataset %<>%
    mutate(Final_Relaxed = Final) %>%
    mutate_at("Final",         ~{ifelse(is.na(.)  & !het         & (RRDR_NON_SILENT | LOF_candidate), 2, .)}) %>%
    mutate_at("Final_Relaxed", ~{ifelse(is.na(.)  & !het_relaxed & (RRDR_NON_SILENT | LOF_candidate), 2, .)}) %>%
    mutate_at(c("Final", "Final_Relaxed"), ~{replace_na(., 3)}) %>%
    select(-pos1, -pos2, -RRDR_NON_SILENT, -LOF_candidate)
  fullDataset %<>%
    group_by(sample_id, drug)
  ## Implement the epistasis rules for AMI, KAN and BDQ in version 2 only
  if (version == 2) {
    fullDataset %<>%
      mutate(excludeEpi_Candidate = FALSE, excludeEpi_Regular = FALSE, excludeEpi_Relaxed = FALSE) %>%
      mutate_at("excludeEpi_Candidate", ~{ifelse(drug_short == "AMI"            , any(!het_strict & gene == "eis"                & effect %in% POOLED_EFFECTS[["LoF"]]), .)}) %>%
      mutate_at("excludeEpi_Candidate", ~{ifelse(drug_short == "KAN"            , any(!het_strict & gene == "eis"                & effect %in% POOLED_EFFECTS[["LoF"]]), .)}) %>%
      mutate_at("excludeEpi_Candidate", ~{ifelse(drug_short %in% c("BDQ", "CFZ"), any(!het_strict & gene %in% STRATIFY_BDQ_GENES & effect %in% POOLED_EFFECTS[["LoF"]]), .)}) %>%
      mutate_at("excludeEpi_Regular"  , ~{ifelse(excludeEpi_Candidate & drug_short == "AMI"             & any(!het & variant %in% EXCLUDE_SET[["AMI"]]), TRUE,
                                          ifelse(excludeEpi_Candidate & drug_short == "KAN"             & any(!het & variant %in% EXCLUDE_SET[["KAN"]]), TRUE,      
                                          ifelse(excludeEpi_Candidate & drug_short %in% c("BDQ", "CFZ") & any(!het         & gene %in% BDQ_GENE & Final <= 2), TRUE, .)))}) %>%
      mutate_at("excludeEpi_Relaxed"  , ~{ifelse(excludeEpi_Candidate & drug_short == "AMI"             & any(!het_relaxed & variant %in% EXCLUDE_SET[["AMI"]]), TRUE,
                                          ifelse(excludeEpi_Candidate & drug_short == "KAN"             & any(!het_relaxed & variant %in% EXCLUDE_SET[["KAN"]]), TRUE,          
                                          ifelse(excludeEpi_Candidate & drug_short %in% c("BDQ", "CFZ") & any(!het_relaxed & gene %in% BDQ_GENE & Final <= 2), TRUE, .)))})
    if (!skipEpistasis) {
      fullDataset %<>%
        ungroup()
      ## Compute the epistasis-specific PPV for AMI and KAN
      epiTabs = vector("list", 4) %>%
        set_names(c("AMI", "KAN", paste0("BDQ_", STRATIFY_BDQ_GENES)))
      for (drugName in c("AMI", "KAN")) {
        curEpiTab = fullDataset %>%
          filter(drug_short == drugName) %>%
          group_by(sample_id) %>%
          filter(!any(variant %in% EXCLUDE_SET[["BOTH"]]) & !(drug_short == "KAN" & sum(variant %in% EXCLUDE_SET[["KAN_EXT"]]) > 1)) %>%
          filter(any(variant %in% EXCLUDE_SET[[paste0(drugName, "_EXT")]] & !het_strict)) %>%
          mutate(Group = ifelse( any(gene == "eis" &  effect %in% POOLED_EFFECTS[["LoF"]] & !het_strict),                      "A", 
                                 ifelse(!any(gene == "eis" & !(effect %in% c(SILENT_EFFECTS, UPSTREAM_VAR)) & Final %in% 1:3), "B", NA))) %>%
          filter(variant %in% EXCLUDE_SET[[paste0(drugName, "_EXT")]] & !het_strict) %>%
          ungroup() %>%
          filter(!is.na(Group)) %>%
          group_by(Group, variant, phenotype) %>%
          mutate(N = n()) %>%
          slice(1) %>%
          ungroup() %>%
          select(Group, variant, phenotype, drug, N)
        epiTabs[[drugName]] = curEpiTab
      }
      ## Compute the epistasis-specific PPV for BDQ
      curEpiTab = fullDataset %>%
        filter(drug_short == "BDQ") %>%
        group_by(sample_id) %>%
        filter(!any(gene %in% EXCLUDE_BDQ_GENES & Final <= 2)) %>%
        filter(any(gene %in% BDQ_GENE & Final <= 2 & !het_strict))
      for (geneName in STRATIFY_BDQ_GENES) {
        subEpiTab = curEpiTab %>%
          mutate(Group = ifelse( any(gene == geneName &  effect %in% POOLED_EFFECTS[["LoF"]] & !het_strict),                    "A", 
                                 ifelse(!any(gene == geneName & !effect %in% c(SILENT_EFFECTS, UPSTREAM_VAR) & Final %in% 1:3), "B", NA))) %>%
          filter(gene %in% BDQ_GENE & Final <= 2 & !het_strict) %>%
          slice(1) %>%
          ungroup() %>%
          filter(!is.na(Group)) %>%
          group_by(Group, phenotype) %>%
          mutate(N = n()) %>%
          slice(1) %>%
          ungroup() %>%
          select(Group, gene, phenotype, drug, N)
        epiTabs[[paste0("BDQ_", geneName)]] = subEpiTab
      }
      fullDataset %<>%
        group_by(sample_id, drug)
    }
  }
  ## Assign the final "regular" and "relaxed" groups to each sample; NB: samples flagged for epistasis will be counted as not fitting the catalogue criteria
  ## The computation below adds MAX_GRADE to any het variant so that any group of interest (1, 2 or 3) can only get determined by relevant non-hets
  fullDataset %<>%
    mutate(Group_Regular = min(Final + het * MAX_GRADE), Group_Relaxed = min(Final_Relaxed + het_relaxed * MAX_GRADE)) %>%
    mutate_at(c("Group_Regular", "Group_Relaxed"), ~{ifelse(. >= 3, . + 1, .)}) %>%
    mutate(extendedCandidate = (effect %in% ADDITIONAL_EFFECTS & (gene %in% unlist(ADDITIONAL_GENES[-1]) | (drug_short == "DLM" & gene %in% EXTENDED_ADD_GENES[["DLM"]])))) %>%
    mutate_at("Group_Regular", ~{ifelse(RIF_geno_Regular & any(!het         & Final == 3 & extendedCandidate & . == 4), 3, .)}) %>%
    mutate_at("Group_Relaxed", ~{ifelse(RIF_geno_Relaxed & any(!het_relaxed & Final == 3 & extendedCandidate & . == 4), 3, .)}) %>%
    select(-extendedCandidate)
  if (version == 2) { 
    fullDataset %<>% 
      mutate_at("Group_Regular", ~{ifelse(excludeEpi_Regular, MAX_GRADE, .)}) %>%
      mutate_at("Group_Relaxed", ~{ifelse(excludeEpi_Relaxed, MAX_GRADE, .)})
  }
  if (safe) { stopifnot(all(testConsistent(fullDataset, groupingVars = c("sample_id", "drug"), consistentVars = c("Group_Regular", "Group_Relaxed"))[[1]])) }
  fullDataset %<>%
    slice(1) %>%
    ungroup()
  ## Calculate the size of groups 1, 2, 3, 1 + 2, and 1 + 2 + 3, with both regular and relaxed thresholds and each RIF_geno status, for each drug
  finalTabs = vector("list", 3 * 2 * 5) %>%
    set_names(outer(c("R", "S", "ALL"), outer(c("Regular", "Relaxed"), STAGES, function(x, y) {paste0(x, "_", y)}), function(z, w) {paste0("RIF", z, "_", w)}))
  for (relax in c(FALSE, TRUE)) {
    if (relax) {
      curDataset = fullDataset %>% 
        mutate(relevantGroup = Group_Relaxed, relevantGeno = RIF_geno_Relaxed)
    } else {
      curDataset = fullDataset %>% 
        mutate(relevantGroup = Group_Regular, relevantGeno = RIF_geno_Regular)
    }
    for (geno in c("R", "S", "ALL")) {
      if (geno != "ALL") {
        curSubset = curDataset %>%
          filter(relevantGeno == ifelse(geno == "R", TRUE, FALSE))
      } else {
        curSubset = curDataset
      }
      for (stage in STAGES) {
        if (stage <= 3) {
          curTab = curSubset %>%
            mutate(selected = (relevantGroup == stage & (stage <= 2 | (stage == 3 & drug %in% GROUP_3_DRUGS))))
        } else {
          if (stage == 12) {
            curTab = curSubset %>%
              mutate(selected = (relevantGroup <= 2))
          } else { ## stage = 123
            curTab = curSubset %>%
              mutate(selected = (relevantGroup <= 2 | (relevantGroup == 3 & drug %in% GROUP_3_DRUGS)))
          }
        }
        curTab %<>% 
          mutate_at("selected", ~{convertToLogical(.)}) %>%
          ungroup()
        curTab %<>%
          group_by(selected, phenotype, drug) %>%
          mutate(N = n_distinct(sample_id)) %>%
          slice(1) %>%
          ungroup() %>%
          select(selected, phenotype, drug, N)
        finalTabs[[paste0("RIF", geno, "_", ifelse(relax, "Relaxed", "Regular"), "_", stage)]] = curTab
      }
    }
  }
  output = list(finalTabs = finalTabs)
  if (!skipEpistasis) {
    output = c(output, list(epiTabs = epiTabs))  
  }
  result = postprocessTabs(output, version = version, relaxed = relaxed, sameRIF = sameRIF)
  result
}

postprocessTabs = function(List, version, relaxed, sameRIF) {
  fullTab = tibble()
  longTabs = List[[1]]
  for (curName in names(longTabs)) {
    curTab = longTabs[[curName]] %>%
      group_by(drug) %>%
      mutate(TP = max(N * as.integer(selected & phenotype == "R")), TN = max(N * as.integer(!selected & phenotype == "S")),
             FP = max(N * as.integer(selected & phenotype == "S")), FN = max(N * as.integer(!selected & phenotype == "R"))) %>%
      slice(1) %>%
      ungroup() %>%
      select(-selected, -phenotype, -N) %>%
      mutate(group = curName, .before = 1)
    fullTab %<>%
      bind_rows(curTab)
  }
  fullTab %<>%
    mutate(PPV  = map2(TP, TP + FP, safeBinomTest), NPV  = map2(TN, TN + FN, safeBinomTest), 
           Sens = map2(TP, TP + FN, safeBinomTest), Spec = map2(TN, FP + TN, safeBinomTest), propR = map2(TP + FN, TP + FN + FP + TN, safeBinomTest)) %>%
    mutate(across(where(is.list), list(lb = ~{map_dbl(., ~{.$conf.int[1]})}, ub = ~{map_dbl(., ~{.$conf.int[2]})}))) %>%
    mutate(across(where(is.list), ~{map_dbl(., ~{.$estimate})}))
  write_csv(fullTab, paste0("SensSpec_Leonid_Version", version, "_", ifelse(relaxed, "Relaxed_", ""), ifelse(sameRIF, "sameRIF_", ""), Sys.Date(), ".csv"))
  output = list(fullTab = fullTab)
  if (length(List) > 1) {
    shortTabs = List[[2]]
    miniTab   = tibble()
    for (curName in names(shortTabs)) {
      curTab = shortTabs[[curName]]
      if ("variant" %in% colnames(curTab)) {
        curTab %<>% rename(criterion = variant)
      } else {
        curTab %<>% rename(criterion = gene) 
      }
      curTab %<>%
        group_by(criterion) %>%
        mutate(LOF_R   = max(N * as.integer(Group == "A" & phenotype == "R")), LOF_S   = max(N * as.integer(Group == "A" & phenotype == "S")),
               nomut_R = max(N * as.integer(Group == "B" & phenotype == "R")), nomut_S = max(N * as.integer(Group == "B" & phenotype == "S"))) %>%
        slice(1) %>%
        ungroup() %>%
        select(-Group, -phenotype, -N) %>%
        mutate(group = curName, .before = 1)
      miniTab %<>%
        bind_rows(curTab)
    }
    miniTab %<>%
      mutate(PPV_LOF = map2(LOF_R, LOF_R + LOF_S, safeBinomTest), PPV_nomut  = map2(nomut_R, nomut_R + nomut_S, safeBinomTest)) %>%
      mutate(across(where(is.list), list(lb = ~{map_dbl(., ~{.$conf.int[1]})}, ub = ~{map_dbl(., ~{.$conf.int[2]})}))) %>%
      mutate(across(where(is.list), ~{map_dbl(., ~{.$estimate})}))
    write_csv(miniTab, paste0("Epistasis_Leonid_Version", version, "_", ifelse(relaxed, "Relaxed_", ""), Sys.Date(), ".csv"))
    output = c(output, list(miniTab = miniTab))
  }
  output
}
