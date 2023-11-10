# FINAL MUTATION CATALOGUE 2023 RESULT FILES
Please find in the folder *Final Result Files* the main output of the mutation catalogue initiative, 2023 version.
This includes:
1. An Excel file with the complete data generated, including variant grades and all statistical metrics per variant
2. A VCF file matching genomic coordinates to variant names only
3. A PDF documenting their use

# RAW DATA FILES
All raw data used as input for the the "SOLO" algorithm can be found in the folder *Input data files for Solo algorithms*

### full_genotypes
Contains all variants (i.e. features of the models) for each sample.  
The data is arranged in two sub-levels of folders, one for drugs and one for tiers. The csv files contain 6 columns.

```
sample_id,resolved_symbol,variant_category,predicted_effect,neutral,"max(af)",position
```

- resolved_symbol is the gene name to which the variant belongs to
- variant_category is the final variant name. It is the result of the categorization algorithm which can be found in the *Genotype ETL - PySpark* folder.
- predicted_effect can have one of the following values:
1. feature_ablation (complete deletion)
2. frameshift
3. inframe_deletion
4. inframe_insertion
5. initiator_codon_variant (start codon variant changing start codon variant)
6. missense_variant
7. non_coding_transcript_exon_variant (non coding change nucleotide changes, i.e. for *rrs/rrl*)
8. start_lost
9. stop_gained
10. stop_lost
11. stop_retained_variant (stop codon changing to another stop codon)
12. synonymous_variant
13. pstream_gene_variant


### phenotype
```
sample_id,phenotypic_category,phenotype,box
```

- phenotypic_category is the overall category the test was assigned to ("WHO" or "ALL")
- phenotype is the result (R or S)
- box is the subcategory of the test. Refer to the mutation catalogue report for more details.

# CODE FILES

All other folders each include code used to generate the mutation catalogue. This includes:

## Bioinformatic pipeline - Step Functions
The [Step Function](https://aws.amazon.com/step-functions/) definition of our bioinformatic pipeline, in JSON format. Refer to the Step Function documentation for the specifications of the format. All executed commands are found in the *Parameters.Parameters.COMMAND* entry of each Step of the Worfklow.

## Variant Annotation - Python
We used SnpEff to annote our variants with the following command: 

```
snpEff ann -config snpEff.config -upDownStreamLen 2000 -spliceSiteSize 0 -spliceRegionExonSize 0 -spliceRegionIntronMax 0 -spliceRegionIntronMin 0 NC_000962.3 tmp.vcf | grep -v "^#" | cut -f3,8 > annotated_file.txt

```

Before annotation with `snpEff ann`, we created a custom reference database with `snpEff build` and using as input the GFF file from the NCBI Genome entry ID [GCA_000195955.2, Assembly version ASM19595v2](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz).

Because of consistent errors of SnpEff, we curated, fixed and normalized its output using the scripts in the *Variant Annotation - Python* folder
The input file read by the script consists of two columns, with the first column being our internal variant id (integer) and the second column the `ANN` INFO tag of the VCF generated by `snpEff ann`.
The script can be modified to expect that column as output only if you wish to replicate our full annotation process.

## Solo algorithm - R implementation
Contains all scripts necessary to run the SOLO algorithm in R.

## Solo algorithm - Stata implementation
Contains all scripts necessary to run the SOLO algorithm in Stata.

## Database Schema - Postgres
The definition of our  database which sustains all our analysis. Partly built using the [biosql](https://github.com/biosql/biosql) specifications. Most of our ETL was implemented using PySpark scripts, found in the *Genotype ETL - PySpark* folder.
