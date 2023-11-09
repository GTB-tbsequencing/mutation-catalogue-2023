ALTER TABLE "country" ADD CONSTRAINT "country_two_letters_code_uniq" UNIQUE ("two_letters_code");
ALTER TABLE "country" ADD CONSTRAINT "country_three_letters_code_uniq" UNIQUE ("three_letters_code");

ALTER TABLE "sample" ADD CONSTRAINT "sample_biosample_id_uniq" UNIQUE("biosample_id");
ALTER TABLE "sample" ADD CONSTRAINT "sample_sample_name_uniq" UNIQUE("sample_name");
ALTER TABLE "sample" ADD CONSTRAINT "sample_sra_name_uniq" UNIQUE("sra_name");

ALTER TABLE "sequencing_data_hash" ADD CONSTRAINT "sequencing_data_hash_value_uniq" UNIQUE("sequencing_data_id", "algorithm", "value");

ALTER TABLE "dataset" ADD CONSTRAINT "dataset_name_uniq" UNIQUE("dataset_name");

ALTER TABLE "contributor" ADD CONSTRAINT "contributor_name_uniq" UNIQUE("contributor_name");
ALTER TABLE "contributor" ADD CONSTRAINT "contributor_email_uniq" UNIQUE("contributor_email");

ALTER TABLE "dataset_to_sample" ADD CONSTRAINT "dataset_to_sample_uniq" UNIQUE("dataset_id", "sample_id");

ALTER TABLE "dataset_to_contributor" ADD CONSTRAINT "dataset_to_contributor_uniq" UNIQUE("dataset_id", "contributor_id");

ALTER TABLE "sequencing_data" ADD CONSTRAINT "sequencing_data_file_path" UNIQUE("file_path");
ALTER TABLE "sequencing_data" ADD CONSTRAINT "sequencing_data_lib_name_file_path_uniq" UNIQUE("library_name", "file_path");

ALTER TABLE "additional_sample_name" ADD CONSTRAINT "additional_sample_name_uniq" UNIQUE("sample_id", "db", "db_label", "sample_name_synonym");

CREATE UNIQUE INDEX "variant_uniq" ON "variant" ("chromosome", "position", (MD5("reference_nucleotide")), (MD5("alternative_nucleotide")));

ALTER TABLE "annotation" ADD CONSTRAINT "annotation_uniq" UNIQUE ("reference_db_crossref_id", "hgvs_value", "predicted_effect");

ALTER TABLE "promoter_distance" ADD CONSTRAINT "promoter_distance_uniq" UNIQUE ("gene_db_crossref_id", "region_start", "region_end");

ALTER TABLE "variant_to_annotation" ADD CONSTRAINT "variant_to_annotation_uniq" UNIQUE ("variant_id", "annotation_id");

ALTER TABLE "genotype" ADD CONSTRAINT "genotype_uniq" UNIQUE("sample_id", "variant_id", "genotyper");

ALTER TABLE "locus_sequencing_stats" ADD CONSTRAINT "locus_sequencing_stats_uniq" UNIQUE("sample_id", "gene_db_crossref_id");

ALTER TABLE "pca_dimension_results" ADD CONSTRAINT "pca_uniq" UNIQUE("sample_id", "dimension");

ALTER TABLE "summary_sequencing_stats" ADD CONSTRAINT "summary_sequencing_stats_uniq" UNIQUE ("sample_id");

ALTER TABLE "taxonomy_stats" ADD CONSTRAINT "taxonomy_stats_uniq" UNIQUE ("sample_id", "ncbi_taxon_id");

ALTER TABLE "gene_drug_resistance_association" ADD CONSTRAINT "gene_drug_resistance_association_uniq" UNIQUE("gene_db_crossref_id", "drug_id");

ALTER TABLE "phenotypic_drug_susceptiblity_test_who_category" ADD CONSTRAINT "phenotypic_drug_susceptiblity_test_who_category_uniq_conc" UNIQUE("drug_id", "medium_id", "concentration");

ALTER TABLE "phenotypic_drug_susceptiblity_test_who_category" ADD CONSTRAINT "phenotypic_drug_susceptiblity_test_who_category_uniq_cat" UNIQUE("drug_id", "medium_id", "category");

ALTER TABLE "microdilution_plate_concentration" ADD CONSTRAINT "microdilution_plate_concentration_uniq" UNIQUE("plate", "drug_id", "concentration");

ALTER TABLE "drug_synonym" ADD CONSTRAINT "drug_synonym_uniq_name" UNIQUE("drug_name_synonym");