import pyspark.sql.functions as F

# Using this view to output convenient gene name instead of gene id
def gene_or_locus_tag_view(sequence_db_crossref, seqfeature_qualifier_value, seqfeature, term, view="gene_symbol"):
    if view=="gene_symbol":
        term_value = "gene"
    elif view=="rv_symbol":
        term_value = "locus_tag"

    qualifier_name = (
        sequence_db_crossref.alias("sdc")
        .join(
            seqfeature_qualifier_value.alias("sqv"),
            on="seqfeature_id",
            how="inner",
        )
        .join(
            seqfeature.alias("seqfeature"),
            on="seqfeature_id",
            how="inner",
        )
        .join(
            term.alias("term1"),
            (F.col("term1.term_id")==F.col("sqv.term_id"))
            & (F.col("term1.name")==term_value),
            how="inner",
        )
        .join(
            term.alias("term2"),
            (F.col("term2.term_id")==F.col("seqfeature.type_term_id"))
            & (F.col("term2.name")=="gene"),
            how="inner",
        )
        .select(
            F.col("sdc.dbxref_id").alias("gene_db_crossref_id"),
            F.col("sqv.value").alias(view),
        )
    )
    return(qualifier_name)

# This view simply coalesces gene_name, locus_tag, to get the most sensible id for each gene
# Usually called the "resolved_symbol" afterwards
def merge_gene_locus_view(gene, locus_tag):
    return(
        locus_tag
        .alias("locus_tag")
        .join(
            gene
            .alias("gene"),
            on="gene_db_crossref_id",
            how="left",
        )
        .select(
            F.col("gene_db_crossref_id"),
            F.coalesce(
                F.col("gene.gene_symbol"),
                F.col("locus_tag.rv_symbol")
            ).alias("resolved_symbol")
        )
    )

# Linking the protein reference id to its associated gene reference id
# Used in the formatted annotation per gene (fapg) view
def protein_id_view (db_crossref, seqfeature_qualifier_value, seqfeature_db_crossref):
    protein_id = (
        db_crossref.alias("dbxref")
        .join(
            seqfeature_qualifier_value.alias("sqv"),
            on=F.col("sqv.value")==F.col("dbxref.accession"),
            how="inner",
        )
        .join(
            seqfeature_db_crossref.alias("sdc"),
            on=F.col("sdc.seqfeature_id")==F.col("sqv.seqfeature_id"),
            how="inner",
        )
        .where(F.col("dbname")=='Protein')
        .select(
            F.col("dbxref.dbxref_id").alias("protein_db_crossref_id"),
            F.col("sdc.dbxref_id").alias("gene_db_crossref_id"),
        )
    )
    return(protein_id)