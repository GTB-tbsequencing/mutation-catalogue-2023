import pandas, openpyxl, io, boto3, datetime, sys

from awsglue.transforms import *
from awsglue.utils import getResolvedOptions
from awsglue.context import GlueContext
from awsglue.job import Job
from awsglue.dynamicframe import DynamicFrame
from pyspark.context import SparkContext
from pyspark.sql.window import Window

import pyspark.sql.functions as F

from biosql_gene_views import *

# Reimplementation of the Materialized View "ranked annotation" in spark

def ranked_annotation(variant_to_annotation, annotation, locus_tag, protein_id, resolved_symbol):
    ranked_annotation = (
        variant_to_annotation
        .join(
            annotation.alias("annotation"),
            "annotation_id",
            "inner"
        )
        .join(
            locus_tag.alias("lt"),
            F.col("lt.gene_db_crossref_id")==F.col("annotation.reference_db_crossref_id"),
            "left"
        )
        .join(
            protein_id.alias("pi"),
            F.col("pi.protein_db_crossref_id")==F.col("annotation.reference_db_crossref_id"),
            "left"
        )
        .join(
            resolved_symbol.alias("resolved"),
            F.col("resolved.gene_db_crossref_id")==F.coalesce(F.col("pi.gene_db_crossref_id"), F.col("lt.gene_db_crossref_id")),
            "inner",
        )
        .withColumn(
            "rank",
            F.rank().over(
                Window.partitionBy(
                    F.col("variant_id"),
                )
                .orderBy(
                    F.when(
                        (~F.col("predicted_effect").isin(['synonymous_variant', 'initiator_codon_variant', 'stop_retained_variant']) & ~F.col("pi.protein_db_crossref_id").isNull())
                        | (F.col("predicted_effect")=='non_coding_transcript_exon_variant'),
                        1
                    )
                    .when(
                        F.col("predicted_effect").isin(['synonymous_variant', 'initiator_codon_variant', 'stop_retained_variant']) & ~F.col("lt.rv_symbol").isNull(),
                        2
                    )
                    .when(
                        F.col("predicted_effect").isin(['upstream_gene_variant', 'downstream_gene_variant']),
                        F.col("distance_to_reference")+2
                    )
                    .asc_nulls_last()
                )
            )
        )
        .select(
            F.col("variant_id"),
            F.col("annotation_id"),
            F.col("resolved.resolved_symbol").alias("gene_name"),
            F.col("predicted_effect"),
            F.col("hgvs_value"),
            F.col("distance_to_reference"),
            F.col("rank")
        )
    )
    return(ranked_annotation)

# Reimplementation of the Materialized View "multiple variant decomposition" in spark
# We need to make sure each multiple variant (e.g. ACG -> TCA) can be decomposed to account for each individual effects
# However we can't do that for multiple variants that fall on the same codon (because the correct change can only be obtained if we consider them simultaneously)
# We make sure that does not happen in the variant categorization later, not in the MVD view
def multiple_variant_decomposition(variant):
    mvd = (
        variant.alias("variant1")
        # Selecting the variants for which the reference and alternative nucleotide both have equal and higher than one length 
        .where(
            (F.length("reference_nucleotide")>1)
            & (F.length("reference_nucleotide")==F.length("alternative_nucleotide"))
        )
        # Splitting both strings at each letter, per row
        .withColumn(
            "zipped",
            F.arrays_zip(
                F.split(F.col("reference_nucleotide"), ""),
                F.split(F.col("alternative_nucleotide"), "")
            ).alias("zip")
        )
        # Explode to new rows the splitted string with posexplode (conveniently outputs the offset as well)
        .select(
            F.col("variant1.variant_id"),
            F.col("variant1.position"),
            F.posexplode(F.col("zipped")),
        )
        # Calculate the new absolute position on the genome using the offset ("pos")
        .withColumn(
            "new_pos",
            F.col("pos")+F.col("variant1.position")
        )
        # Expand the arrays of [new_ref, new_alt] to columns
        .select(
            F.col("variant1.variant_id"),
            F.col("variant1.position"),
            F.col("new_pos"),
            F.col("col.0").alias("new_ref"),
            F.col("col.1").alias("new_alt"),
        )
        # Select the exploded rows where the letter differs (otherwise it's not a variant)
        .where(
            F.col("new_ref")!=F.col("new_alt")
        )
        # Joining back to variant to obtain the variant_id of the sub variants
        .join(
            variant.alias("variant2"),
            on=(F.col("variant2.position")==F.col("new_pos"))
            & (F.col("variant2.reference_nucleotide")==F.col("new_ref"))
            & (F.col("variant2.alternative_nucleotide")==F.col("new_alt")),
            how="inner",
        )
        .select(
            F.col("variant1.variant_id").alias("mnv_id"),
            F.col("variant2.variant_id").alias("snv_id"),
        )
        .distinct()
    )
    return(mvd)

# The infamous "FAPG"
def formatted_annotation_per_gene(variant_to_annotation, annotation, db_crossref, protein_id):
    fapg = (
        variant_to_annotation.alias("vta1")
        .join(
            annotation.alias("annotation1"),
            on="annotation_id",
            how="inner",
        )
        # First get all annotation at nucleotidic level
        .join(
            db_crossref.alias("dbxref"),
            on=F.col("dbxref.dbxref_id")==F.col("annotation1.reference_db_crossref_id"),
            how="inner",
        )
        # Get proteic annotation as well if it exists
        .join(
            # Using the variant_to_annotation table again
            other=variant_to_annotation.alias("vta2")
                .join(
                    annotation.alias("annotation2"),
                    on="annotation_id",
                    how="inner",
                )
                # Extract only proteic annotation with an inner join that associates the reference protein id to the reference gene id
                # Because annotation at a nucleotide level are linked to the gene id, whereas protein annotation are linked to protein id
                .join(
                    protein_id.alias("protein_id"),
                    on=F.col("protein_id.protein_db_crossref_id")==F.col("annotation2.reference_db_crossref_id"),
                    how="inner",
                ).alias("protein_annotation"),
            on=(F.col("protein_annotation.variant_id")==F.col("vta1.variant_id"))
            & (F.col("protein_annotation.gene_db_crossref_id")==F.col("annotation1.reference_db_crossref_id")),
            how="left",
        )
        # Selecting annotations of interest, ie either an effect on the protein, or upstream/rRNA gene variant/deletion
        # Basically everything but downstream variants
        # If a variant encompasses the upstream region + the CDS of a gene 
        # (e.g. deletion: 
        # 13557 ATGTATTCCT A upstream_gene_variant c.-8_1delAGGAATACA p.Met1fs
        # 13557 ATGTATTCCT A frameshift_variant&start_lost c.-8_1delAGGAATACA p.Met1fs)
        # We only keep the predicted effect that is different to upstream_gene_variant
        # Otherwise the same variant would be extracted twice.
        .where(
            (F.col("protein_annotation.reference_db_crossref_id").isNotNull() & (F.col("annotation1.predicted_effect")!="upstream_gene_variant"))
            | ((F.col("annotation1.predicted_effect")=="upstream_gene_variant") & F.col("protein_annotation.hgvs_value").isNull())
            | (F.col("annotation1.predicted_effect").isin(["non_coding_transcript_exon_variant", "feature_ablation"]))            
        )
        .select(
            F.col("vta1.variant_id"),
            F.col("annotation1.reference_db_crossref_id").alias("gene_db_crossref_id"),
            F.col("annotation1.predicted_effect"),
            F.col("annotation1.hgvs_value").alias("nucleotidic_annotation"),
            F.col("protein_annotation.hgvs_value").alias("proteic_annotation"),
        # Update Dec 2023
        # We needed to sanitize the distance to reference for some missense variants
        # Indeed, for MCNV, the distance given by default by snpEff is the start of the MCNV (obvious)
        # Not the location of the actual missense (eg if the MCNV starts with a synonymous and a missense on the next codon)
        # So we extract the position from the proteic annotation and simply multiply by 3.
            F.when(
                F.col("annotation1.hgvs_value").rlike("c\.[0-9]+_[0-9]+del[AGCT]+ins[AGCT]+")
                    & F.col("protein_annotation.hgvs_value").rlike("p\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}")
                    & (F.col("annotation1.predicted_effect")=="missense_variant"),
                F.regexp_extract(
                    F.col("protein_annotation.hgvs_value"),
                    r"p\.[A-Z][a-z]{2}([0-9]+)[A-Z][a-z]{2}",
                    1
                ).cast("long")*3
            )
            .otherwise(F.col("annotation1.distance_to_reference"))
            .alias("distance_to_reference"),
        )
    )
    return(fapg)

# This logic is needed to correctly handle MVD leading to global synonymous variants (i.e. no change at all & overall on the protein)
# The aim is to avoid incongruities when multiple changes on the same codon lead to a synonymous change, although each atomic change could lead to a missense 
# In that case, we must not decompose the MVD but instead choose the complex annotation in the form of HGV "delins" 
# If all atomic change also lead to synonymous change, then we can decompose.
def sanitize_synonymous_variant(formatted_annotation_per_gene, tier, multiple_variant_decomposition):
    # First select MNV syn with only atomic SNV syn
    sanitized = (
        formatted_annotation_per_gene.alias("fapg1")
        .join(
            tier,
            on="gene_db_crossref_id",  
            how="left_semi"
        )
        .where(
            F.col("fapg1.predicted_effect")=='synonymous_variant'
        )
        .join(
            multiple_variant_decomposition.alias("mvd"),
            on=F.col("fapg1.variant_id")==F.col("mvd.mnv_id"),
            how="inner"
        )
        .join(
            formatted_annotation_per_gene.alias("fapg2"),
            (F.col("mvd.snv_id")==F.col("fapg2.variant_id"))
            & (F.col("fapg1.gene_db_crossref_id")==F.col("fapg2.gene_db_crossref_id")),
            how="left",
        )
        .groupBy(
            F.col("fapg1.variant_id")
        )
        .agg(
            F.expr("bool_and(fapg2.predicted_effect='synonymous_variant')").alias("all_syn")
        )
        .where(
            F.col("all_syn")
        )
        .select("variant_id")
        .distinct()
        .alias("sanitized")
    )

    # Construct the mapping between the MVD pure-syn and the atomic nucleotidic annotations
    correct_syn_mnv = (
        multiple_variant_decomposition.alias("mvd")
        .join(
            sanitized,
            on=F.col("sanitized.variant_id")==F.col("mvd.mnv_id"),
            how="inner"
        )
        .join(
            formatted_annotation_per_gene.alias("fapg1"),
            on=F.col("fapg1.variant_id")==F.col("mvd.snv_id"),
            how="inner"
        )
        .select(
            F.col("sanitized.variant_id"),
            F.col("fapg1.gene_db_crossref_id"),
            F.col("fapg1.predicted_effect"),
            F.col("fapg1.nucleotidic_annotation"),
            F.col("fapg1.proteic_annotation"),
            F.col("fapg1.distance_to_reference"),
        )
        .alias("correct_syn_mnv")
    )
    return(correct_syn_mnv)

# I think I need to sanitize similarly MNVs
# Because if one MNVs is associated with more than one missense
# I can't count on a simple join on different distance to reference to exclude irrelevant synonymous changes
# The synonymous will be join with missense on other codons, but it still irrelevant...
def missense_codon_list(formatted_annotation_per_gene, variant, tier):
    mnvs_to_missense = (
        formatted_annotation_per_gene.alias("fapg")
        .join(
            tier,
            on="gene_db_crossref_id",  
            how="left_semi"
        )
        .join(
            variant.alias("variant"),
            on="variant_id",
            how="inner"
        )
        .where(
            (F.col("predicted_effect")=='missense_variant')
                & F.col("nucleotidic_annotation").rlike("c\.[0-9]+_[0-9]+del[AGCT]+ins[AGCT]+")
                & F.col("proteic_annotation").rlike("p\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}")
                & (F.length("reference_nucleotide")==F.length("alternative_nucleotide"))            
        )
        .select(
            F.col("variant_id"),
            F.col("position"),
            F.col("reference_nucleotide"),
            F.col("alternative_nucleotide"),
            F.col("gene_db_crossref_id"),
            F.col("proteic_annotation"),
            F.floor(1+(F.col("distance_to_reference")-1)/3).alias("codon"),
        )
        .groupBy(
            F.col("variant_id"),
            F.col("gene_db_crossref_id"),
            # F.col("position"),
            # F.col("reference_nucleotide"),
            # F.col("alternative_nucleotide"),
        )
        .agg(
            F.countDistinct("proteic_annotation").alias("count_missense"),
            F.first("gene_db_crossref_id"),
            F.concat_ws(
                "; ",
                F.collect_set(F.col("proteic_annotation"))
            ).alias("mutations"),
            F.collect_set(F.col("codon")).alias("all_codons_missense")
        )
        # .alias("mnv_to_missense")
    )

    return(
        mnvs_to_missense
        # .select(
        #     F.col("variant_id"),
        #     F.col("gene_db_crossref_id"),
        #     F.col("all_codons_missense"),
        # )
    )

    # selected_syn = (
    #     multiple_variant_decomposition.alias("mvd")
    #     .join(
    #         formatted_annotation_per_gene.alias("fapg2"),
    #         on=(F.col("mvd.snv_id")==F.col("fapg2.variant_id"))
    #             & (F.col("fapg2.predicted_effect")=="synonymous_variant"),
    #         how="inner"
    #     )
    #     .alias("synonymous")
    # )

    # mnvs_missense_syn = (
    #     mnvs_to_missense
    #     .join(
    #         selected_syn,
    #         on=(F.col("synonymous.mnv_id")==F.col("mnv_to_missense.variant_id"))
    #             & ~F.array_contains("all_codons_missense", F.floor(1+(F.col("synonymous.distance_to_reference")-1)/3)),
    #         how="left"
    #     )
    # )

    # return(mnvs_missense_syn)

# Beautiful (actually awful) query that associates variant_id to its variant category on the gene of interest first
# Variants leading to more than one missense mutation (on different codons) are already handled by decomposition of the SnpEff annotation 
# in the UpdateAnnotation ETL StepFunction pipeline, by load_annotation_variants.py
# However, multiple mutation on the same codons prevents us from using straightforwardly the MVD
def tiered_drug_variant_categories(
        formatted_annotation_per_gene,
        corrected_synonymous_mnv,
        tier,
        variant,
        multiple_variant_decomposition,
        promoter_distance,
        missense_codon_list
    ):
    variant_categories = (
        formatted_annotation_per_gene.alias("fapg1")
        .join(
            tier,
            on="gene_db_crossref_id",  
            how="inner"
        )
        .join(
            variant
            .alias("variant1")
            .select("variant_id", "position")
            .distinct(),
            on="variant_id",
            how="inner"
        )
        # Joining FAPG twice on the multiple variants decomposed into single variants
        # Joining on nucleotide only variants (syn, up, non coding) + missense
        # We still consider missense, to handle cases where a MNV leads to missense + synonymous to nearby codon
        .join(
            multiple_variant_decomposition.alias("mvd"),
            on=(F.col("fapg1.variant_id")==F.col("mvd.mnv_id"))
                & F.col("fapg1.predicted_effect").isin(
                    'upstream_gene_variant',
                    'synonymous_variant',
                    'missense_variant',
                    'non_coding_transcript_exon_variant',
                ),
            how="left",
        )
        .join(
            missense_codon_list.alias("missense"),
            on=["variant_id", "gene_db_crossref_id"],
            how="left",
        )
        # Join again with variant to get the position of the mvd
        .join(
            variant
            .alias("variant2")
            .select("variant_id", "position")
            .distinct(),
            F.col("mvd.snv_id")==F.col("variant2.variant_id"),
            how="left"
        )
        # First join the atomic synonymous variants
        # They have priority
        .join(
            corrected_synonymous_mnv.alias("cor_syn_mnv"),
            on=(F.col("cor_syn_mnv.variant_id")==F.col("fapg1.variant_id"))
                & (F.col("cor_syn_mnv.gene_db_crossref_id")==F.col("fapg1.gene_db_crossref_id")),
            how="left"
        )
        # Do not join additional missense annotations on the same codon
        # We already have the correct one !
        .join(
            formatted_annotation_per_gene.alias("fapg2"),
            (F.col("mvd.snv_id")==F.col("fapg2.variant_id"))
            & (F.col("fapg1.gene_db_crossref_id")==F.col("fapg2.gene_db_crossref_id"))
            & (
                # Extract neighboring variants only if they only have a nucleotidic impact (ie no MNV on same codon)
                (
                    F.col("fapg1.predicted_effect").isin(
                        'upstream_gene_variant',
                        'non_coding_transcript_exon_variant'
                        )
                    & (F.col("fapg1.predicted_effect")==F.col("fapg2.predicted_effect"))
                )
                # Imagine an MNV that leads to a missense in one codon and a synonymous on the next one.
                # We extract both missense and synonymous on the other 
                |
                (
                    (F.col("fapg1.predicted_effect")=="missense_variant")
                    & F.col("fapg2.predicted_effect").isin(
                        "synonymous_variant",
                        "stop_retained_variant",
                        "initiator_codon_variant"
                        )
                    & ~F.array_contains("missense.all_codons_missense", F.floor(1+(F.col("fapg2.distance_to_reference")-1)/3))
                )
            ),
            how="left",
        )
        .where(
            ~F.col("fapg1.predicted_effect").isin('downstream_gene_variant')
        )
        .withColumn(
            "variant_category",
            F.when(
                F.col("fapg1.predicted_effect")=="feature_ablation",
                F.lit("deletion")
            )
            # Sometimes more than one effect is predicted
            # Fixing frameshifts on the stop codon. This is a stop lost, not a frameshift
            .when(
                F.col("fapg1.predicted_effect").rlike(".*frameshift.*") & F.col("fapg1.proteic_annotation").rlike("p\.Ter[0-9]+fs"),
                F.regexp_replace(F.col("fapg1.proteic_annotation"), "fs", "ext*?")
            )
            # We rank then 1. stop_lost+frameshift 2. stop_lost, 3. frameshift, stop_gained, start_lost, 4. initiator_codon_variant
            .when(
                F.col("fapg1.predicted_effect").rlike(".*stop_lost.*") & F.col("fapg1.predicted_effect").rlike(".*frameshift.*"),
                F.col("fapg1.nucleotidic_annotation")
            )
            .when(
                F.col("fapg1.predicted_effect").rlike(".*stop_lost.*"),
                F.col("fapg1.proteic_annotation")
            )
            # For start lost, we sanitize all variant categories with Met
            .when(
                F.col("fapg1.predicted_effect").rlike(".*start_lost.*"),
                F.lit("p.Met1?")
            )
            # Stop gained is proteic annotation
            # Sanitizing stop gained that are also frameshifts
            # We keep the stop gained annotation
            .when(
                F.col("fapg1.predicted_effect").rlike(".*stop_gained.*"),
                F.regexp_replace(F.col("fapg1.proteic_annotation"), "fs$", "*")
            )
            # frameshift is now nucleic annotation
            .when(
                F.col("fapg1.predicted_effect").rlike(".*frameshift.*"),
                F.col("fapg1.nucleotidic_annotation")
            )
            # Sanitize initiator codon variant so that the category is different from the start_lost
            # Getting the nucleotidic annotation for initiator codon variants
            .when(
                F.col("fapg1.predicted_effect").isin(
                    "initiator_codon_variant", 
                    "initiator_codon_variant&non_canonical_start_codon",
                    "stop_retained_variant"
                ),
                F.col("fapg1.nucleotidic_annotation")
            )
            # Extract the FAPG of the multiple decomposed variant for non protein associated change
            # The corrected synonymous mnvs have priority, then undecomposed global delins if the change can't be sanitized
            .when(
                F.col("fapg1.predicted_effect").isin(
                    ['upstream_gene_variant', 'non_coding_transcript_exon_variant', 'synonymous_variant']
                ),
                F.coalesce(
                    F.col("cor_syn_mnv.nucleotidic_annotation"),
                    F.col("fapg2.nucleotidic_annotation"),
                    F.col("fapg1.nucleotidic_annotation")
                )
            )
            # Then extract FAPG for variants leading to mis+syn
            .when(
                (F.col("fapg1.predicted_effect")=="missense_variant")
                & F.col("fapg2.predicted_effect").isin(
                    "synonymous_variant",
                    "stop_retained_variant",
                    "initiator_codon_variant"
                ),
                F.col("fapg2.nucleotidic_annotation")
            )
            .otherwise(F.col("fapg1.proteic_annotation"))
        )
        .withColumn(
            "curated_effect",
            # Sanitizing inframe indel. Disruptive/conservative brings no useful information
            # For instance, the same amino acid duplication could be caused by both disruptive and conservative annotations
            F.when(
                F.col("fapg1.predicted_effect").rlike(".*start_lost.*"),
                F.lit("start_lost")
            )
            .when(
                F.col("fapg1.predicted_effect").rlike(
                    ".*((?:disruptive|conservative)_inframe_(?:deletion|insertion)).*"
                ),
                F.regexp_extract(
                    F.col("fapg1.predicted_effect"),
                    r"(inframe_(?:deletion|insertion))",
                    1
                )
            )
            # That is the case when a synonymous variant is phased with a missense variant on another codon
            .when(
                (F.col("fapg1.predicted_effect")=="missense_variant")
                & F.col("fapg2.predicted_effect").isin("synonymous_variant", "stop_retained_variant", "initiator_codon_variant"),
                F.col("fapg2.predicted_effect")
            )
            # Always favour frameshift over stop gained if both are mentionned in the predicted effect
            # Fixing frameshifts on the stop codon. This is a stop lost, not a frameshift
            .when(
                F.col("fapg1.predicted_effect").rlike(".*frameshift.*") & F.col("variant_category").rlike("p\.Ter[0-9]+ext\*\?"),
                F.lit("stop_lost")
            )
            # Stop gained has priority over frameshift
            .when(
                F.col("fapg1.predicted_effect").rlike(".*stop_gained.*"),
                F.lit("stop_gained")
            )
            .when(
                F.col("fapg1.predicted_effect").rlike(".*frameshift.*"),
                F.lit("frameshift")
            )
            .when(
                F.col("fapg1.predicted_effect").rlike(".*(stop_gained|feature_ablation).*"),
                F.regexp_extract(
                    F.col("fapg1.predicted_effect"),
                    r".*(stop_gained|feature_ablation).*",
                    1
                )
            )
            .when(
                F.col("fapg1.predicted_effect").rlike(".*stop_lost.*"),
                F.lit("stop_lost")
            )
            .when(
                F.col("fapg1.predicted_effect").rlike(".*initiator_codon_variant.*"),
                F.lit("initiator_codon_variant")
            )
            .otherwise(F.col("fapg1.predicted_effect"))
        )
        # Sometimes one (gene_db_crossref_id, variant_category) pair can have more than 1 predicted effect!
        # Think of (rpoB, p.Thr444dup) : can be "disruptive_inframe_insertion", "conservative_inframe_insertion".
        # Do not use "predicted_effect" for later groupings!
        # Update : this is now fixed by sanitizing the predicted effect
        # e.g. both "disruptive_inframe_(insertion|deletion)" and "conservative_inframe_(insertion|deletion)" are sanitized into "inframe_(insertion|deletion)"
        .select(
            F.col("drug_id"),
            F.col("fapg1.variant_id").alias("variant_id"),
            F.col("fapg1.gene_db_crossref_id"),
            F.col("curated_effect"),
            F.col("variant_category"),
            F.col("tier"),
            F.when(F.col("fapg2.predicted_effect").isNotNull(), F.col("variant2.position"))
                .otherwise(F.col("variant1.position")).alias("position"),
            F.coalesce(
                F.col("fapg2.distance_to_reference"),
                F.col("fapg1.distance_to_reference"),
            ).alias("distance_to_reference"),
        )
        .join(
            other=promoter_distance.alias("prom_dist"),
            on=(F.col("prom_dist.gene_db_crossref_id")==F.col("fapg1.gene_db_crossref_id"))
                & (F.col("curated_effect")=="upstream_gene_variant")
                & F.col("distance_to_reference").between(F.col("prom_dist.region_start"), F.col("prom_dist.region_end")),
            how="left",
        )
        # Exclude all upstream gene variant that are not located inside promoter regions for each gene
        .where(
            (F.col("curated_effect")!="upstream_gene_variant")
            | F.col("prom_dist.gene_db_crossref_id").isNotNull()
        )
        # We rank the variant categories so that each 
        # Order is identical if a variant has more than one category for a single gene i.e. :
        # If a variant leads to 2 missense changes on the same gene, we keep both
        .withColumn(
            "rank_all_categories",
            F.rank().over(
                Window.partitionBy(
                    F.col("variant_id"),
                    F.col("drug_id")
                )
                .orderBy(
                    F.when(F.col("curated_effect")=="lof", 1)
                    .when(F.col("curated_effect").isin("missense_variant", "synonymous_variant", "stop_retained_variant", "initiator_codon_variant"), 2)
                    .when(F.col("curated_effect")=="non_coding_transcript_exon_variant", 3)
                    .otherwise(4)
                )
            )
        )
        .where(F.col("rank_all_categories")==1)
        # We need to handle the case where a promoter region is shared between different gene
        # We rank the variants through tier, gene_id, 
        # Variants leading to more than one change in the promoter (think MVD) appear twice relative to the chosen gene only
        # Genes are ordered by ascending gene ids, so basically the first one along the genome will be chosen
        .withColumn(
            "rank_upstream_categories",
            F.rank().over(
                Window.partitionBy(
                    F.col("variant_id"),
                    F.col("drug_id")
                )
                .orderBy(
                    F.asc("tier"),
                    F.asc("fapg1.gene_db_crossref_id"),
                )
            )
        )
        .where(F.col("rank_upstream_categories")==1)
        .select(
            F.col("drug_id"),
            F.col("tier"),
            F.col("fapg1.gene_db_crossref_id").alias("gene_db_crossref_id"),
            F.col("variant_id"),
            F.col("position"),
            F.col("variant_category"),
            F.col("curated_effect").alias("predicted_effect"),
            F.col("distance_to_reference"),
        )
        .distinct()
    )
    return(variant_categories)


if __name__=="__main__":

    d = datetime.datetime.now().isoformat()

    args = getResolvedOptions(sys.argv, ['JOB_NAME', "glue_database_name", "sample_fraction"])

    glueContext = GlueContext(SparkContext.getOrCreate())

    spark = glueContext.spark_session

    spark._jsc.hadoopConfiguration().set('spark.sql.broadcastTimeout', '3600')

    job = Job(glueContext)

    dbxref = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_biosql_dbxref")

    sqv = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_biosql_seqfeature_qualifier_value")

    sdc = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_biosql_seqfeature_dbxref")

    protein_id = protein_id_view(dbxref, sqv, sdc).alias("protein_id")

    location = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_biosql_location")



    end_positive_genes = (
        location
        .join(
            sdc.alias("sdc"),
            "seqfeature_id",
            "inner"
        )
        .join(
            protein_id.alias("prot"),
            on=F.col("prot.gene_db_crossref_id")==F.col("sdc.dbxref_id"),
            how="inner"
        )
        .where(
            F.col("strand")==1
        )
        .select(
            F.col("gene_db_crossref_id"),
            (F.col("end_pos")-2).alias("end_pos"),
        )
        .distinct()
    )

    print(end_positive_genes.count())

    print(end_positive_genes.show())

    annot = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_annotation")

    vta = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_variant_to_annotation")

    fapg = formatted_annotation_per_gene(vta, annot, dbxref, protein_id).alias("fapg")

    promoter_distance = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_promoter_distance")

    tier = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_gene_drug_resistance_association").alias("tier")

    s3 = boto3.resource('s3')

    # csv_buffer = io.BytesIO()
    # fapg.join(tier, "gene_db_crossref_id", "left_semi").toPandas().to_csv(csv_buffer, index=False)
    # s3.Object("aws-glue-assets-231447170434-us-east-1", args["JOB_NAME"]+"/"+d+"_"+args["JOB_RUN_ID"]+"/fapg.csv").put(Body=csv_buffer.getvalue())

    variant = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_genphensql_variant")

    mvd = multiple_variant_decomposition(variant).alias("mvd")
    
    san = sanitize_synonymous_variant(fapg, tier, mvd)

    mnvs_miss = missense_codon_list(fapg, variant, tier)    

    # csv_buffer = io.BytesIO()
    # mnvs_miss.toPandas().to_csv(csv_buffer, index=False)
    # s3.Object("aws-glue-assets-231447170434-us-east-1", args["JOB_NAME"]+"/"+d+"_"+args["JOB_RUN_ID"]+"/mnvs_miss.csv").put(Body=csv_buffer.getvalue())

    var_cat = tiered_drug_variant_categories(fapg, san, tier, variant, mvd, promoter_distance, mnvs_miss)

    seqfeature = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_biosql_seqfeature").alias("seqfeature")

    term = glueContext.create_data_frame.from_catalog(database = args["glue_database_name"], table_name = "postgres_biosql_term").alias("term1")

    gene_name = gene_or_locus_tag_view(sdc, sqv, seqfeature, term, "gene_symbol").alias("gene")

    locus_tag = gene_or_locus_tag_view(sdc, sqv, seqfeature, term, "rv_symbol").alias("locus_tag")

    gene_locus_tag = merge_gene_locus_view(gene_name, locus_tag).alias("gene_locus_tag")


    variant_mapping = (
        var_cat.alias("var_cat")
        .drop(
            "position"
        )
        .join(
            variant.alias("variant1"),
            "variant_id",
            "inner"
        )
        .join(
            gene_locus_tag,
            "gene_db_crossref_id",
            "inner"
        )
        # trying to remove incorrect stop codon annotations
        # snpEff has a bug when trying to left align MCNVs on positive strand stop codons
        # it associated them with missense annotation a little upstream of the stop codons
        # which is wrong
        .join(
            end_positive_genes.alias("stop_codons"),
            on=(F.col("stop_codons.gene_db_crossref_id")==F.col("var_cat.gene_db_crossref_id"))
                & (
                    (
                        (F.col("variant1.position")==F.col("stop_codons.end_pos")) & (F.length(F.col("reference_nucleotide"))==F.length(F.col("alternative_nucleotide"))) & (F.length(F.col("reference_nucleotide"))==3)
                    )
                    |
                    (
                        (F.col("variant1.position")==(F.col("stop_codons.end_pos")+1)) & (F.length(F.col("reference_nucleotide"))==F.length(F.col("alternative_nucleotide"))) & (F.length(F.col("reference_nucleotide"))==2)
                    )
                    |
                    (
                        (F.col("variant1.position")-F.col("stop_codons.end_pos")).isin([0, 1])
                        &
                        (F.col("predicted_effect")=="frameshift")
                    )
                ),
            how="left_anti"
        )
        .select(
            F.concat_ws(
                "_",
                F.col("resolved_symbol"),
                F.col("variant_category")
            ).alias("variant"),
            F.col("resolved_symbol"),
            F.col("variant_category"),
            F.col("predicted_effect"),
            F.col("var_cat.variant_id"),
            F.col("variant1.chromosome"),
            F.col("variant1.position"),
            F.col("variant1.reference_nucleotide"),
            F.col("variant1.alternative_nucleotide"),
        )
        .distinct()
        .sample(fraction=float(args["sample_fraction"]))
        .toPandas()
    )




#    writer = pandas.ExcelWriter(output, engine="openpyxl")

    
    #writer.save()
    
    #data = output.getvalue()

    #s3.Bucket('aws-glue-assets-231447170434-us-east-1').put_object(Key=args["JOB_NAME"]+"/"+d+"_"+args["JOB_RUN_ID"]+"/variant_mapping.csv", Body=data)

    csv_buffer = io.BytesIO()
    variant_mapping.to_csv(csv_buffer, index=False)
    s3.Object("aws-glue-assets-231447170434-us-east-1", args["JOB_NAME"]+"/"+d+"_"+args["JOB_RUN_ID"]+"/variant_mapping.csv").put(Body=csv_buffer.getvalue())


    # fill_new_table = (
    #     DynamicFrame.fromDF(
    #         var_cat
    #         .select(
    #             F.col("variant_id"),
    #             F.col("gene_db_crossref_id"),
    #             F.col("predicted_effect"),
    #             F.col("variant_category.variant_category"),
    #             F.col("var_cadistance_to_reference"),
    #         )
    #         .alias("")
    #         .distinct(),
    #         glueContext,
    #         "final"
    #     )
    #     .resolveChoice(
    #         choice="match_catalog",
    #         database= args["glue_database_name"],
    #         table_name="postgres_genphensql_tiered_variant_categories"
    #     )
    # )

    # datasink5 = glueContext.write_dynamic_frame.from_catalog(frame = fill_new_table, database = args["glue_database_name"], table_name = "postgres_genphensql_tiered_variant_categories", transformation_ctx = "datasink5")


    # fill_new_table_ranked_annotation = (
    #     DynamicFrame.fromDF(
    #         ranked_annotation(vta, annot, locus_tag, protein_id, gene_locus_tag)
    #         .alias("")
    #         .distinct(),
    #         glueContext,
    #         "final"
    #     )
    #     .resolveChoice(
    #         choice="match_catalog",
    #         database= args["glue_database_name"],
    #         table_name="postgres_genphensql_ranked_annotation"
    #     )
    # )

    # datasink6 = glueContext.write_dynamic_frame.from_catalog(frame = fill_new_table_ranked_annotation, database = args["glue_database_name"], table_name = "postgres_genphensql_ranked_annotation", transformation_ctx = "datasink6")

    job.commit()
