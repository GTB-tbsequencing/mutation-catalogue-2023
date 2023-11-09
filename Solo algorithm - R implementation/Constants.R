## List of all drugs and their short names 
DRUG_LIST   = c("Amikacin", "Bedaquiline", "Capreomycin", "Delamanid", "Ethambutol", "Ethionamide", "Isoniazid", "Kanamycin", 
                "Levofloxacin", "Linezolid", "Moxifloxacin", "Pyrazinamide", "Rifampicin", "Streptomycin", "Clofazimine", "Pretomanid")
SHORT_NAMES = c("AMI", "BDQ", "CAP", "DLM", "EMB", "ETH", "INH", "KAN", 
                "LEV", "LZD", "MXF", "PZA", "RIF", "STM", "CFZ", "PTM")

## Special-status drugs
FQS              = DRUG_LIST[match(c("LEV", "MXF")                     , SHORT_NAMES)]
INH_ETH          = DRUG_LIST[match(c("INH", "ETH")                     , SHORT_NAMES)]
BDQ_CFZ          = DRUG_LIST[match(c("BDQ", "CFZ")                     , SHORT_NAMES)]
FIRST_LINE_DRUGS = DRUG_LIST[match(c("EMB", "INH", "PZA", "RIF")       , SHORT_NAMES)]
GROUP_3_DRUGS    = DRUG_LIST[match(c("BDQ", "CFZ", "DLM", "INH", "PZA"), SHORT_NAMES)]

## Special-status genes
BDQ_GENE     = "Rv0678"
CFZ_GENE     = "pepQ"
BDQ_CFZ_GENE = c(BDQ_GENE, CFZ_GENE)
CAP_GENE     = "tlyA"
DLM_GENE     = c("ddn", "fbiA", "fbiB", "fbiC", "fgd1", "Rv2983")
ETH_GENE     = "ethA"
FQ_GENE      = c("gyrA", "gyrB")
INH_GENE     = "katG"
INH_ETH_GENE = c("inhA","fabG1")
PZA_GENE     = "pncA"
RRDR_GENE    = "rpoB"
STM_GENE     = "gid"
EXTRA_GENES  = c("mshA", "panD", "Rv1979c")

## Special-status gene sets
ADDITIONAL_GENES   = list(DLM = DLM_GENE, INH = INH_GENE, PZA = PZA_GENE)
SPECIAL_URM_GENES  = c(CAP_GENE, CFZ_GENE, ETH_GENE, STM_GENE, EXTRA_GENES, unlist(ADDITIONAL_GENES)) %>% magrittr::set_names(NULL)
EXTENDED_ADD_GENES = c(ADDITIONAL_GENES, list(BDQ = BDQ_CFZ_GENE, CFZ = BDQ_CFZ_GENE, CAP = CAP_GENE, ETH = ETH_GENE, STM = STM_GENE))
ADDITIONAL_GENES   = c(ADDITIONAL_GENES, list(BDQ = BDQ_GENE, CFZ = BDQ_GENE))
SUB_EXTENDED_GENES = EXTENDED_ADD_GENES[c("ETH", "INH", "PZA", "STM")]
EXCLUDE_BDQ_GENES  = c("atpE",  "pepQ")
STRATIFY_BDQ_GENES = c("mmpL5", "mmpS5")

## Special-status variants
BORDERLINE_MUTATIONS = c("p.Leu430Pro", "p.Asp435Tyr", "p.His445Leu", "p.His445Asn", "p.His445Ser", "p.Leu452Pro", "p.Ile491Phe")
INFLATED_PPV_VARS    = c("p.Asn98fs"  , "p.Cys46Arg" , "p.Cys46fs"  , "p.Gln51fs"  , "p.Ile67Ser" , "p.Leu142fs" , "p.Met146Thr", "p.Pro48fs")
PZA_SPEC_VAR         = paste(PZA_GENE, "p.Ile31Thr", sep = "_")

## Special-status effects
POOLED_EFFECTS       = list(LoF = c("feature_ablation", "frameshift", "start_lost", "stop_gained"))
INFRAME_EFFECTS      = paste0("inframe_", c("deletion", "insertion"))
ADDITIONAL_EFFECTS   = c(POOLED_EFFECTS[["LoF"]], INFRAME_EFFECTS, "missense_variant")
BDQ_EFFECTS          = c(ADDITIONAL_EFFECTS, "stop_lost")
SILENT_EFFECTS       = c("initiator_codon_variant", "stop_retained_variant", "synonymous_variant")

## Useful constants and thresholds
RRDR_INT              = 426:452
PPV_UB_THRESHOLD      = 0.1
SIG_THRESHOLD         = 0.05
FDR_THRESHOLD         = 0.05
CI_COEFFICIENT        = 1.96
MAF_THRESHOLD_STRICT  = 0.75
MAF_THRESHOLD_REGULAR = 0.75
MAF_THRESHOLD_RELAXED = 0.25

## Grades used to prioritise variants
GRADES = paste0(c("", "not "), "assoc w ") %>% 
  str_to_sentence() %>%
  outer(c("R", "R - Interim"), function(x, y) { paste0(x, y) }) %>%
  c("Uncertain significance", "Manual check") %>%
  magrittr::extract(c(1, 3, 5, 4, 2, 6)) %>%
  enframe(value = "grading") %>%
  mutate(grading = paste0(name, ") ", grading)) %>%
  pull(grading)
INITIAL_FLAG = 6
FINAL_FLAG   = 6
MAX_GRADE    = 6

## Process-specific tables
PHENO_GROUPS       = tibble(category_phenotype = c("ALL", "WHO", "CC", "CC-ATU"), group = c("MAIN", "MAIN", "CC", "ATU"))
BAD_VAR_DRUG_PAIRS = tibble(drug = c("Isoniazid", "Rifampicin"), variant = c("katG_p.Ser315Thr", "rpoB_p.Ser450Leu"))
CONVERTED_FILES    = paste0(c("neutral_mutations_catalogue_v1", "new_variant_matched_to_old", "v1_grades"), ".csv")

## Additional information supporting the old grades
OLD_LOF    = "Indel or premature stop codon (LoF)"
ALL_ONLY   = "Evidence from ALL dataset only"
PASS_2     = "Algorithm pass 2"
RRDR       = "RRDR"
BORDERLINE = "Borderline"
PROBLEMATIC_CRITERION = "WHO-endorsed gDST assay"
IGNORED_CRITERIA = c(ALL_ONLY, PASS_2, OLD_LOF, RRDR, BORDERLINE)

## Additional information supporting the new grades
INFLATION    = "Potentially inflated PPV"
UPSTREAM_VAR = "upstream_gene_variant"
SET_C_ONLY   = "Neutrals defined by setC (literature) only"
V1_NEUTRALS  = "previous WHO guidance"
WHO_BASED    = "Interim based on WHO dataset"
MANUAL_CHECK = "Manual check required"
SILENT       = "Silent mutation"
LoF_MUT      = "Indel frameshift or premature stop codon (LoF)"
WHO_ASSAY    = "Recognized as DR marker through WHO-endorsed assay"
ALLELIC_EXP  = "Selection evidence"
PRE_GUIDANCE = "Previous WHO guidance"
V1_EVIDENCE  = "Evidence from catalogue version 1"
FQ_CROSS_RES = "FQ cross-resistance"
IE_CROSS_RES = "INH-ETH cross-resistance"
BC_CROSS_RES = "BDQ-CFZ cross-resistance"

## Constants used for SOLO algorithm results
MAX_ITER = Inf
UCODE    = 0L
RCODE    = 1L
SCODE    = -1L
CODE_KEY = c("S", "U", "R") %>% magrittr::set_names(c(SCODE, UCODE, RCODE))
HET_TAB  = tibble(class = c("S", "S", "U", "U", "R"), het = c(FALSE, TRUE, FALSE, TRUE, FALSE))
HET_CNT  = table(HET_TAB$class)
## and the final tabulation
STAGES   = c(1, 2, 3, 12, 123)

## Epistasis-related variants
EXCLUDE_SET = vector("list", 5)
EXCLUDE_SET[["AMI"]]     =                         paste0("eis_c.", "-14C>T")
EXCLUDE_SET[["AMI_EXT"]] =   EXCLUDE_SET[["AMI"]]
EXCLUDE_SET[["KAN"]]     = c(EXCLUDE_SET[["AMI"]], paste0("eis_c.", "-10G>A"))
EXCLUDE_SET[["KAN_EXT"]] = c(EXCLUDE_SET[["KAN"]], paste0("eis_c.", c("-12C>T", "-37G>T", "-8delC")))
EXCLUDE_SET[["BOTH"]]    =                         paste0("rrs_n.", c("1401A>G", "1402C>T", "1484G>T"))
