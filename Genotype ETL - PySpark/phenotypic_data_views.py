import sys, boto3, os, io, pandas, datetime
from itertools import chain

from awsglue.transforms import *
from awsglue.utils import getResolvedOptions
from pyspark.context import SparkContext
from awsglue.context import GlueContext
from awsglue.job import Job
from pyspark.sql import SQLContext
from awsglue.dynamicframe import DynamicFrame
from pyspark.sql.window import Window

import pyspark.sql.functions as F
from pyspark.sql.window import Window

# Assigning pDST data to the classification by Claudio Koser & Paolo Miotto
# On top of using the predefined categories of the phenotype_category table
# We use the concentration values of some WHO_current and WHO past
# To automatically classify resistant or susceptible results that are compatible with the WHO preferred CC

def join_phenotypes_with_categories(phenotype, phenotype_category, drug, growth_medium, assessment_method):
    # For each (drug id, medium id, method id) we select the WHO_current or WHO_past (ranked in that order) concentration if they exist
    select_who_cc = (
        phenotype_category
        .where(
            F.col("category").isin("WHO_past", "WHO_current")
        )
        .withColumn(
            "category_rank",
            F.rank().over(
                Window.partitionBy(
                    F.col("drug_id"),
                    F.col("method_id"),
                    F.col("medium_id"),
                )
                .orderBy(
                    F.when(F.col("category")=="WHO_current", 1)
                    .when(F.col("category")=="WHO_past", 2)
                )
            )
        )
        .where(
            F.col("category_rank")==1
        )
        .join(
            assessment_method,
            on="method_id",
            how="left"
        )
        .alias("")
    )

    # "Agar" on the medium column basically means uncertainty on whether "7H10" or "7H11"
    # So selecting the maximum and the miminum WHO CC for each drug on these mediums
    # Allows automatic classification of compatible R/S results into non_WHO_CC_[R|S]
    agar_cc = (
        select_who_cc
        .join(
            growth_medium,
            on="medium_id",
            how="inner"
        )
        .where(
            F.col("medium_name").isin("Middlebrook7H10", "Middlebrook7H11")
            & (F.col("method_name").isNull() | F.col("method_name").isin("Proportions"))
        )
        .groupBy(
            "drug_id",
        )
        .agg(
            F.max("concentration").alias("max(concentration)"),
            F.min("concentration").alias("min(concentration)"),
         )
        .select(
            F.col("drug_id"),
            F.col("min(concentration)"),
            F.col("max(concentration)"),
        )
        # Finally "concatenate" as a column the medium_id corresponding to Agar 
        .crossJoin(
            growth_medium.alias("agar").where(F.col("medium_name")=="Agar")
        )
    )

    # We can now perform the automatic classification
    # Very ugly with 4 successive & similar joins
    joined_phenotype_category = (
        phenotype
        .alias("phenotype")
        .join(
            growth_medium.alias("growth_medium"),
            on=F.col("phenotype.medium_id")==F.col("growth_medium.medium_id"),
            how="left"
        )
        .join(
            assessment_method.alias("assessment_method"),
            on=F.col("phenotype.method_id")==F.col("assessment_method.method_id"),
            how="left"
        )
        .join(
            other=phenotype_category.alias("phenotype_category1"),
            on=(F.col("phenotype_category1.drug_id")==F.col("phenotype.drug_id")) 
                # We need to coalesce null entries to perform an adequate join
                # So that we can join on null values
                & (F.coalesce(F.col("phenotype_category1.medium_id"), F.lit(0))==F.coalesce(F.col("phenotype.medium_id"), F.lit(0)))
                & (F.coalesce(F.col("phenotype_category1.method_id"), F.lit(0))==F.coalesce(F.col("phenotype.method_id"), F.lit(0)))
                & (F.coalesce(F.col("phenotype_category1.concentration"), F.lit(0))==F.coalesce(F.col("phenotype.concentration"), F.lit(0))),
            how="left"
        )
        .join(
            other=select_who_cc.withColumn("corrected_category", F.lit("non_WHO_CC_R")).alias("phenotype_category2"),
            on=(F.col("phenotype_category2.drug_id")==F.col("phenotype.drug_id"))
                & (F.coalesce(F.col("phenotype_category2.medium_id"), F.lit(0))==F.coalesce(F.col("phenotype.medium_id"), F.lit(0)))
                & ((F.coalesce(F.col("phenotype_category2.method_id"), F.lit(0))==F.coalesce(F.col("phenotype.method_id"), F.lit(0)))
                    | (F.col("phenotype_category2.method_name").isNull() & (F.col("assessment_method.method_name")=="Proportions"))
                    | (F.col("assessment_method.method_name").isNull() & (F.col("phenotype_category2.method_name")=="Proportions"))
                    )
                & (F.col("phenotype_category2.concentration")<F.col("phenotype.concentration"))
                & (F.col("phenotype.test_result")=="R"),
            how="left"
        )
        .join(
            other=select_who_cc.withColumn("corrected_category", F.lit("non_WHO_CC_S")).alias("phenotype_category3"),
            on=(F.col("phenotype_category3.drug_id")==F.col("phenotype.drug_id"))
                & (F.coalesce(F.col("phenotype_category3.medium_id"), F.lit(0))==F.coalesce(F.col("phenotype.medium_id"), F.lit(0)))
                & ((F.coalesce(F.col("phenotype_category3.method_id"), F.lit(0))==F.coalesce(F.col("phenotype.method_id"), F.lit(0)))
                    | (F.col("phenotype_category3.method_name").isNull() & (F.col("assessment_method.method_name")=="Proportions"))
                    | (F.col("assessment_method.method_name").isNull() & (F.col("phenotype_category3.method_name")=="Proportions"))
                    )
                & (F.col("phenotype_category3.concentration")>F.col("phenotype.concentration"))
                & (F.col("phenotype.test_result")=="S"),
            how="left"
        )
        .join(
            other=agar_cc.withColumn("corrected_category", F.lit("non_WHO_CC_S")).alias("agar1"),
            on=(F.col("agar1.drug_id")==F.col("phenotype.drug_id"))
                & (F.col("phenotype.medium_id")==F.col("agar1.medium_id"))
                & ((F.col("assessment_method.method_name")=="Proportions") | F.col("assessment_method.method_name").isNull())
                & (F.col("agar1.min(concentration)")>F.col("phenotype.concentration"))
                & (F.col("phenotype.test_result")=="S"),
            how="left"
        )
        .join(
            other=agar_cc.withColumn("corrected_category", F.lit("non_WHO_CC_R")).alias("agar2"),
            on=(F.col("agar2.drug_id")==F.col("phenotype.drug_id"))
                & (F.col("phenotype.medium_id")==F.col("agar2.medium_id"))
                & ((F.col("assessment_method.method_name")=="Proportions") | F.col("assessment_method.method_name").isNull())
                & (F.col("agar2.max(concentration)")<F.col("phenotype.concentration"))
                & (F.col("phenotype.test_result")=="R"),
            how="left"
        )
        .join(
            drug.alias("drug"),
            on=F.col("phenotype.drug_id")==F.col("drug.drug_id"),
            how="inner"
        )
        .withColumn(
            "phenotypic_category",
            # If Moxi and MGIT, we consider the concentration between new and old CC (ie 0.5 and 1)
            # almost as good as the WHO_current 
            F.when(
                (F.col("drug_name")=="Moxifloxacin")
                & (F.col("growth_medium.medium_name")=="MGIT")
                & F.col("phenotype.concentration").isin(0.5, 1),
                F.lit("non_WHO_CC_SR")
            )
            # Some categories always have priority, for safety
            .when(
                F.col("phenotype_category1.category").isin("WHO_past", "WHO_current", "to_be_excluded"),
                F.col("phenotype_category1.category")
            )
            .otherwise(
                F.coalesce(
                    F.col("phenotype_category3.corrected_category"),
                    F.col("phenotype_category2.corrected_category"),
                    F.col("agar1.corrected_category"),
                    F.col("agar2.corrected_category"),
                    F.col("phenotype_category1.category")
                )
            )
        )
        .where(
            F.col("phenotypic_category").isNotNull()
        )
        .select(
            "phenotype.*",
            F.col("phenotypic_category")
        )
    )
    return(joined_phenotype_category.alias(""))

# Ofloxacin tests are assigned to Levofloxacin, and prothionamide to ethionamide
# We need the drug information to access the actual names
def pooling_drugs(categorized_phenotypes, drug):
    pooled_phenotypes = (
        categorized_phenotypes
        .alias("categorized_phenotype")
        .join(
            drug.alias("original"),
            on="drug_id",
            how="inner"
        )
        .join(
            drug.alias("corrected")
                .withColumnRenamed("drug_id", "corrected_drug_id"),
            on=((F.col("original.drug_name")=="Ofloxacin") & (F.col("corrected.drug_name")=="Levofloxacin"))
                | ((F.col("original.drug_name")=="Prothionamide") & (F.col("corrected.drug_name")=="Ethionamide")),
            how="left"
        )
        .withColumn(
            "drug_id",
            F.coalesce(
                F.col("corrected_drug_id"),
                F.col("categorized_phenotype.drug_id"),
            )
        )
        .select(
            F.col("categorized_phenotype.*"),
            F.col("drug_id")
        )
    )
    return(pooled_phenotypes.alias(""))

# Ranking of the final pDST categories
# If a pair of (box, medium) has discordant result (i.e. more than one, and at least 1 S and 1 R)
# All the results from that pair are discarded and we move to the next (box, medium) pair
def extract_clean_phenotypes_rank_medium(categorized_phenotypes, medium):
    clean_phenotypes = (
        categorized_phenotypes
        .where(
            F.col("phenotypic_category").isNotNull()
            & ~F.col("phenotypic_category").isin("to_be_excluded", "non_WHO_CC")
        )
        .join(
            medium,
            "medium_id",
            "left"
        )
        .groupBy(
            F.col("sample_id"),
            F.col("drug_id"),
            F.col("medium_name"),
            F.col("phenotypic_category"),
        )
        # Concatenating all results from a given (sample, drug, box, medium)
        # Collect_set to fetch unique values of test_result
        .agg(
            F.concat_ws(
                "",
                F.collect_set(F.col("test_result"))
            ).alias("concat_results")
        )
        # Remove incompatible results per (sample, drug, box, medium)
        .where(
            F.col("concat_results").isin("R", "S")
        )
        # Rank the medium inside each box
        .withColumn(
            "rank_medium",
            F.rank().over(
                Window.partitionBy(
                    F.col("sample_id"),
                    F.col("drug_id"),
                    F.col("phenotypic_category")
                )
                .orderBy(
                    F.when(F.col("medium_name")=="Middlebrook7H10", 1)
                    .when(F.col("medium_name")=="LJ", 2)
                    .when(F.col("medium_name")=="Middlebrook7H11", 3)
                    .when(F.col("medium_name")=="Agar", 4)
                    .when(F.col("medium_name")=="MGIT", 5)
                    .when(F.col("medium_name")=="MODS", 6)
                    .when(F.col("medium_name")=="BACTEC460", 7)
                    .otherwise(8)
                )
            )
        )
        # Select highest ranking medium at each box
        .where(
            F.col("rank_medium")==1
        )
        .select(
            F.col("sample_id"),
            F.col("drug_id"),
            F.col("phenotypic_category"),
            F.col("concat_results").alias("test_result"),
        )
    )
    return(clean_phenotypes.alias(""))

# Perform the "box" (sub-category) ranking
def extract_clean_phenotypes_rank_category(ranked_medium_phenotypes):
    clean_phenotypes = (
        ranked_medium_phenotypes
        .where(
            F.col("phenotypic_category").isNotNull()
            & ~F.col("phenotypic_category").isin("to_be_excluded", "non_WHO_CC")
        )
        .groupBy(
            F.col("sample_id"),
            F.col("drug_id"),
            F.col("phenotypic_category"),
        )
        # Afer ranking on the medium we already have only one row left per (sample, drug, category)
        # However, MIC tests do not go through medium ranking, so we need to group again in case some samples have 2 MIC entries on the same category
        .agg(
            F.concat_ws(
                "",
                F.collect_set(F.col("test_result"))
            ).alias("concat_results")
        )
        .where(
            F.col("concat_results").isin("R", "S")
        )
        # Rank the categories for all given (sample, drug) 
        .withColumn(
            "rank_category",
            F.rank().over(
                Window.partitionBy(
                    F.col("sample_id"),
                    F.col("drug_id"),
                )
                .orderBy(
                    F.when(F.col("phenotypic_category")=="WHO_current", 1)
                    .when(F.col("phenotypic_category")=="non_WHO_CC_SR", 2)
                    .when(F.col("phenotypic_category")=="WHO_past", 3)
                    .when(F.col("phenotypic_category")=="WHO_undefined", 4)
                    .when(F.col("phenotypic_category")=="non_WHO_CC_R", 5)
                    .when(F.col("phenotypic_category")=="non_WHO_CC_S", 6)
                    .when(F.col("phenotypic_category")=="CRyPTIC_MIC", 7)
                    .when(F.col("phenotypic_category")=="MYCOTB_MIC", 8)
                )
            )
        )
        # Select highest ranking category
        .where(
            F.col("rank_category")==1
        )
        # And we are done
        .select(
            F.col("sample_id"),
            F.col("drug_id"),
            F.col("phenotypic_category"),
            F.col("concat_results").alias("phenotype"),
        )
    )
    return(clean_phenotypes)

# PySpark does not handle numeric range as Postgres does
# So we receive the data as strings and use regexps to extract the lower and upper boundaries
# Very ugly but no other way around it
# Sometimes identical numrange are cast to different string values by PySpark. Not relevant as the bounds are then cast to double anyway.
# But the output is not nicely formatted.

def binarize_mic_test(minimum_inhibitory_concentration, epidemiological_cutoff_values, plate_concentration, drug, CC_ATU=False):
    if not CC_ATU:
        ecoff_values = (
            epidemiological_cutoff_values
            .where(
                F.col("medium_name").isin("UKMYC6", "UKMYC5", "MGIT", "LJ", "7H10", "7H11", "MYCOTB")
            )
        )
    else:
        new_breakpoint = (
            plate_concentration.alias("plate_conc")
            .where(
                F.col("plate").isin("UKMYC6", "UKMYC5", "MYCOTB")
            )
            .withColumn(
                "lagging_value",
                F.lag("concentration").over(
                        Window.partitionBy("plate", "drug_id")
                        .orderBy("concentration")
                )
            )
            .fillna(0)
            .select(
                F.col("drug_id"),
                F.col("plate"),
                F.col("concentration"),
                F.col("lagging_value"),
            )
        )
        ecoff_values = (
            epidemiological_cutoff_values.alias("ecoff")
            .where(
                F.col("medium_name").isin("UKMYC6", "UKMYC5", "MGIT", "LJ", "7H10", "7H11", "MYCOTB")
            )
            .join(
                new_breakpoint.alias("new_break"),
                on=(F.col("ecoff.drug_id")==F.col("new_break.drug_id"))
                    & (F.col("ecoff.medium_name")==F.col("new_break.plate"))
                    & (F.col("ecoff.value")==F.col("new_break.concentration")),
                how="left"
            )
            .select(
                F.col("ecoff.drug_id"),
                F.col("ecoff.medium_name"),
                F.coalesce(
                    F.col("new_break.lagging_value"),
                    F.col("value")/2
                ).alias("value")
            )
        )
    binarized_plates = (
        minimum_inhibitory_concentration
        .join(
            drug,
            "drug_id",
            "inner"
        )
        # All these drugs are not on the MYCOTB plate...
        .where(
            ~(
                (F.col("plate")=="MYCOTB")
                & F.col("drug_name").isin(
                    "Capreomycin", "Cycloserine", "Para-Aminosalicylic Acid", "Rifabutin", "Streptomycin"
                )
            )
        )
        # Manually exclude impossible results for the moment...
        # I believe this is query is incorrect as testing on ranges shouldn't performed via strings !!!!!
        # I should define columns upper and lower boundaries first then perform the condition check on them
        .where(
            ~(
                (F.col("plate")=="MYCOTB")
                & (F.col("drug_name")=="Isoniazid")
                & (F.col("mic_value")=="[0.08,0.08]")
            )
        )
        .alias("mic")
        .join(
            ecoff_values.alias("ecoff"),
            on=(F.col("ecoff.drug_id")==F.col("mic.drug_id"))
                & (F.col("mic.plate")==F.col("ecoff.medium_name")),
            how="inner"
        )
        # Same regexp, extract different groups
        # Cleaner than two different regexpes for each group
        .withColumn(
            "upper_boundary",
            F.regexp_extract(
                F.col("mic_value"),
                r'[\(|\[](.*),(.*)[\]\)]',
                2
            ).cast("double")
        )
        .withColumn(
            "lower_boundary",
            F.regexp_extract(
                F.col("mic_value"),
                r'[\(|\[](.*),(.*)[\]\)]',
                1
            ).cast("double")
        )
        # When the upper boundary is equal or lower than the cutoff (say 1), it's S and has priority
        # i.e. (0.5, 1] or [1, 1] 
        # Then if the lower boundary is equal or higher than the ecoff, it's R.
        # i.e. (1, 2] or [2, 2]
        # When none of this is true, for instance (0.5, 2], it's null
        # Because the ecoff was not tested.
        .withColumn(
            "test_result",
            F.when(F.col("upper_boundary")<=F.col("ecoff.value"), 'S')
            .when(F.col("lower_boundary")>=F.col("ecoff.value"), 'R')
        )
        .where(
            F.col("test_result").isin('R', 'S')
        )
        .select(
                F.col("sample_id"),
                F.col("mic.drug_id").alias("drug_id"),
                F.col("plate"),
                F.col("mic_value"),
                F.col("test_result"),
                F.when(F.col("mic.plate")=="MYCOTB", F.lit("MYCOTB_MIC"))
                    .when(F.col("mic.plate").isin("UKMYC5", "UKMYC6"), F.lit("CRyPTIC_MIC"))
                    .when((F.col("drug_name")=="Pretomanid") & (F.col("mic.plate")=="MGIT"), F.lit("WHO_current"))
                    .otherwise(F.concat(F.lit("non_WHO_CC_"), F.col("test_result"))).alias("phenotypic_category"),
        )
    )
    return(binarized_plates.alias(""))


# These are specific filters asked by Paolo & Claudio for the CC/CC_ATU categories

def filter_tests_for_CC_CC_ATU(mic, drug):
    bin_mic_cc_cc_atu = (
        mic.alias("mic_cc_cc-atu")
        .join(
            drug,
            "drug_id",
            "inner"
        )
        .where(
            F.col("plate").isin("UKMYC5", "UKMYC6")
            | ((F.col("drug_name")=="Amikacin")
                & (F.col("plate")=="MYCOTB"))
            | ((F.col("drug_name")=="Bedaquiline")
                & F.col("plate").isin("7H11", "MGIT"))
            | ((F.col("drug_name")=="Clofazimine")
                & (F.col("plate")=="MGIT"))
            | ((F.col("drug_name")=="Ethambutol")
                & (F.col("plate")=="MYCOTB"))
            | ((F.col("drug_name")=="Isoniazid")
                & (F.col("plate")=="MYCOTB"))
            | ((F.col("drug_name")=="Linezolid")
                & F.col("plate").isin("7H10", "MGIT"))
            | ((F.col("drug_name")=="Moxifloxacin")
                & F.col("plate").isin("MGIT", "7H11", "MYCOTB"))
            | ((F.col("drug_name")=="Ofloxacin")
                & (F.col("plate")=="MYCOTB"))
            | ((F.col("drug_name")=="Rifampicin")
                & (F.col("plate")=="MYCOTB"))
            | ((F.col("drug_name")=="Streptomycin")
                & (F.col("plate")=="LJ"))
            | ((F.col("drug_name")=="Pretomanid")
                & (F.col("plate")=="MGIT"))
        )
    )

    return(bin_mic_cc_cc_atu)

    
def rename_CC_CC_ATU_phenotypic_category(binarized_mic_test, CC_ATU=False):
    column_prefix = "CC_"
    if CC_ATU:
        column_prefix = "CC-ATU_"

    return(
        binarized_mic_test
        .withColumn(
            "phenotypic_category",
            F.when(~F.col("phenotypic_category").isin("CRyPTIC_MIC", "MYCOTB_MIC"), column_prefix+"WHO_MIC")
            .otherwise(F.concat(F.lit(column_prefix), F.col("phenotypic_category")))
        )
    )

# Reciprocal filter for CC/CC_ATU
# Some samples that have more than one CRyPTIC entry can appear on each of CC/CC-ATU but not both
# For instance, two entries for AMK : (0.25,0.5] and (0.5,1]
# For CC, it's S & S
# However, for CC-ATU, it's S & R
# Thus we need to make sure that we discard those samples in both CC & CC-ATU
# No need to join on the phenotypic category because we only have one category left per sample when we perform the reciprocal filter
def reciprocal_filter_cc_cc_atu(binarized_mic_cc, binarized_mic_cc_atu):

    reciprocal_binarized_mic_cc_atu = (
        binarized_mic_cc_atu
        .join(
            binarized_mic_cc,
            on=["sample_id", "drug_id"],
            how="left_semi"
        )
    )

    reciprocal_binarized_mic_cc = (
        binarized_mic_cc
        .join(
            binarized_mic_cc_atu,
            on=["sample_id", "drug_id"],
            how="left_semi"
        )
    )

    return(reciprocal_binarized_mic_cc, reciprocal_binarized_mic_cc_atu)


def preparing_binary_data_for_final_algorithm(phenotypes, phenotypes_category, drug, growth_medium, method, mic, ecoff, tier, plate_conc):

    # Select samples based on their phenotype 
    # Only keep clean samples (ie all results either R or S)
    categorized_phenotypes = (
        join_phenotypes_with_categories(phenotypes, phenotypes_category, drug, growth_medium, method)
    )

    # Ofloxacin->levofloxacin, prothionamide->ethionamide
    pooled_lev_eth_phenotypes = (
        pooling_drugs(categorized_phenotypes, drug)
        .select(
            F.col("sample_id"),
            F.col("drug_id"),
            F.col("medium_id"),
            F.col("phenotypic_category"),
            F.col("test_result"),
        )
    )

    # Binarizing all mics
    pooled_lev_eth_mics_cryptic_mycotb = (
        pooling_drugs(
            binarize_mic_test(mic, ecoff, plate_conc, drug, CC_ATU=False),
            drug
        )
        .select(
            F.col("sample_id"),
            F.col("drug_id"),
            F.col("phenotypic_category"),
            F.col("test_result"),
        )
    )        

    pooled_pdst_mic_medium_ranked = (
        extract_clean_phenotypes_rank_medium(
            pooled_lev_eth_phenotypes, growth_medium
        )
        .unionByName(
            pooled_lev_eth_mics_cryptic_mycotb
        )
    )

    clean_phenotypes = (
        extract_clean_phenotypes_rank_category(
            pooled_pdst_mic_medium_ranked
        )
        .join(
            tier,
            on="drug_id",
            how="left_semi",
        )
    )

    return(clean_phenotypes)

def preparing_cc_cc_atu_data_for_final_algorithm(mic, drug, plate_conc, ecoff, tier):

    all_mic_cc_cc_atu = (
        filter_tests_for_CC_CC_ATU(mic, drug)
        .drop("drug_name")
    )

    bin_mic_cc = (
        extract_clean_phenotypes_rank_category(
            pooling_drugs(
                binarize_mic_test(
                    all_mic_cc_cc_atu,
                    ecoff,
                    plate_conc,
                    drug,
                    CC_ATU=False
                ),
                drug
            )
        )
        .join(
            tier,
            "drug_id",
            "left_semi"
        )
        .select(
            F.col("sample_id"),
            F.col("drug_id"),
            F.col("phenotype"),
            F.col("phenotypic_category")
        )
    )

    bin_mic_cc_atu = (
        extract_clean_phenotypes_rank_category(
            pooling_drugs(
                binarize_mic_test(
                    all_mic_cc_cc_atu,
                    ecoff,
                    plate_conc,
                    drug,
                    CC_ATU=True
                ),
                drug
            )
        )
        .join(
            tier,
            "drug_id",
            "left_semi"
        )
        .select(
            F.col("sample_id"),
            F.col("drug_id"),
            F.col("phenotype"),
            F.col("phenotypic_category")
        )
    )

    #Reciprocal filtering
    final_bin_mic_cc_cc_atu = reciprocal_filter_cc_cc_atu(
        bin_mic_cc,
        bin_mic_cc_atu
    )

    final_bin_mic_cc = rename_CC_CC_ATU_phenotypic_category(
        final_bin_mic_cc_cc_atu[0],
        CC_ATU=False
    )

    final_bin_mic_cc_atu = rename_CC_CC_ATU_phenotypic_category(
        final_bin_mic_cc_cc_atu[1],
        CC_ATU=True
    )

    return(final_bin_mic_cc, final_bin_mic_cc_atu)

# If running as main, we output a summary excel of ALL tests incerted, with minimal filterings
# Used for quality control of the binarization process for instance

if __name__ == "__main__":

    args = getResolvedOptions(sys.argv, ['JOB_NAME', "glue_database_name"])

    spark_Cont = SparkContext.getOrCreate()

    glueContext = GlueContext(spark_Cont)

    spark_session = glueContext.spark_session

    job = Job(glueContext)

    d = datetime.datetime.now().isoformat()

    # Main sources
    phenotypes = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_phenotypic_drug_susceptibility_test")
    
    phenotypes_category = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_phenotypic_drug_susceptibility_test_category")

    drug = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_drug")

    growth_medium = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_growth_medium")

    method = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_phenotypic_drug_susceptibility_assessment_method")

    categorized_phenotypes = join_phenotypes_with_categories(phenotypes, phenotypes_category, drug, growth_medium, method)

    categories_count = (
        categorized_phenotypes.alias("categorized")
        .join(
            drug,
            "drug_id",
            "inner"
        )
        .join(
            growth_medium,
            "medium_id",
            "left"
        )
        .join(
            method,
            "method_id",
            "left"
        )
        .where(
            F.col("test_result").isin('R', 'S')
        )
        .groupBy(
            F.col("drug_name"),
            F.col("medium_name"),
            F.col("method_name"),
            F.col("concentration"),
            F.col("phenotypic_category"),
            F.col("test_result"),
        )
        .count(            
        )
        .select(
            F.col("drug_name"),
            F.lit("pDST").alias("type"),
            F.when(
                F.col("method_name").isNull(), F.col("medium_name")
            )
            .otherwise(
                F.concat(F.col("medium_name"), F.lit(" ("), F.col("method_name"), F.lit(")"))
            ).alias("medium_name"),
            F.col("concentration"),
            F.col("test_result"),
            F.col("phenotypic_category"),
            F.col("count"),
        )
    )

    mic = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_minimum_inhibitory_concentration_test").alias("mic")

    ecoff = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_epidemiological_cut_off_value").alias("ecoff")

    tier = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_gene_drug_resistance_association").alias("tier")

    plate_conc = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_microdilution_plate_concentration").alias("drug")

    bin_mics = binarize_mic_test(mic, ecoff, plate_conc, drug, CC_ATU=False)

    mics_counts = (
        bin_mics.alias("categorized_mics")
        .join(
            drug,
            "drug_id",
            "inner"
        )
        .groupBy(
            F.col("drug_name"),
            F.col("plate"),
            F.col("mic_value"),
            F.col("test_result"),
            F.col("phenotypic_category"),
        )
        .count(            
        )
        .select(
            F.col("drug_name"),
            F.lit("MIC").alias("type"),
            F.col("plate"),
            F.col("mic_value"),
            F.col("test_result"),
            F.col("phenotypic_category")
        )
    )

    bin_mic_cc_cc_atu = filter_tests_for_CC_CC_ATU(mic, drug).drop("drug_name")

    bin_mic_cc = rename_CC_CC_ATU_phenotypic_category(
        binarize_mic_test(
            bin_mic_cc_cc_atu,
            ecoff,
            plate_conc,
            drug,
            CC_ATU=False
        ),
        CC_ATU=False
    )


    bin_mic_cc_counts = (
        bin_mic_cc
        .join(
            drug,
            "drug_id",
            "inner"
        )
        .groupBy(
            F.col("drug_name"),
            F.col("plate"),
            F.col("mic_value"),
            F.col("phenotypic_category"),
            F.col("test_result"),
        )
        .count()
        .select(
            F.col("drug_name"),
            F.lit("MIC").alias("type"),
            F.col("plate"),
            F.col("mic_value"),
            F.col("test_result"),
            F.col("phenotypic_category"),
            F.col("count")
        )
    )

    bin_mic_cc_atu = rename_CC_CC_ATU_phenotypic_category(
        binarize_mic_test(
            bin_mic_cc_cc_atu,
            ecoff,
            plate_conc,
            drug,
            CC_ATU=True
        ),
        CC_ATU=True
    )

    bin_mic_cc_atu_counts = (
        bin_mic_cc_atu
        .join(
            drug,
            "drug_id",
            "inner"
        )
        .groupBy(
            F.col("drug_name"),
            F.col("plate"),
            F.col("mic_value"),
            F.col("phenotypic_category"),
            F.col("test_result"),
        )
        .count()
        .select(
            F.col("drug_name"),
            F.lit("MIC").alias("type"),
            F.col("plate"),
            F.col("mic_value"),
            F.col("test_result"),
            F.col("phenotypic_category"),
            F.col("count")
        )
    )

    s3 = boto3.resource('s3')

    output = io.BytesIO()
    writer = pandas.ExcelWriter(output, engine="openpyxl")
    categories_count.toPandas().to_excel(writer, sheet_name="Phen Cat Count", index=False)
    mics_counts.toPandas().to_excel(writer, sheet_name="MIC Cat Count", index=False)
    bin_mic_cc_counts.toPandas().to_excel(writer, sheet_name="CC", index=False)
    bin_mic_cc_atu_counts.toPandas().to_excel(writer, sheet_name="CC-ATU", index=False)
    writer.save()
    data = output.getvalue()
    
    s3.Bucket('aws-glue-assets-231447170434-us-east-1').put_object(Key=args["JOB_NAME"]+"/"+d+"_"+args["JOB_RUN_ID"]+"/categories.xlsx", Body=data)