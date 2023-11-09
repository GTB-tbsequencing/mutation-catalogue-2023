ALTER TABLE "sample" ADD FOREIGN KEY ("ncbi_taxon_id")  REFERENCES "taxon"("ncbi_taxon_id");
ALTER TABLE "sample" ADD FOREIGN KEY ("country_id")  REFERENCES "country"("country_id");

ALTER TABLE "dataset" ADD FOREIGN KEY ("bioproject_id") REFERENCES "bioproject"("bioproject_id");

ALTER TABLE "sequencing_data_hash" ADD FOREIGN KEY ("sequencing_data_id") REFERENCES "sequencing_data"("sequencing_data_id");

ALTER TABLE "dataset_to_sample" ADD FOREIGN KEY ("dataset_id") REFERENCES "dataset"("dataset_id");
ALTER TABLE "dataset_to_sample" ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");

ALTER TABLE "dataset_to_contributor" ADD FOREIGN KEY ("dataset_id") REFERENCES "dataset"("dataset_id");
ALTER TABLE "dataset_to_contributor" ADD FOREIGN KEY ("contributor_id") REFERENCES "contributor"("contributor_id");

ALTER TABLE "sequencing_data" ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");

ALTER TABLE "pca_dimension_results" ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");

ALTER TABLE "tiered_variant_categories" ADD FOREIGN KEY ("variant_id") REFERENCES "variant"("variant_id");
ALTER TABLE "tiered_variant_categories" ADD FOREIGN KEY ("gene_db_crossref_id") REFERENCES "dbxref"("dbxref_id");

ALTER TABLE "additional_sample_name" ADD FOREIGN KEY ("sample_id")  REFERENCES "sample"("sample_id");

ALTER TABLE "annotation" ADD FOREIGN KEY ("reference_db_crossref_id") REFERENCES dbxref(dbxref_id);

ALTER TABLE "variant_to_annotation" ADD FOREIGN KEY ("variant_id") REFERENCES "variant"("variant_id");
ALTER TABLE "variant_to_annotation" ADD FOREIGN KEY ("annotation_id") REFERENCES "annotation"("annotation_id");

ALTER TABLE "genotype" ADD FOREIGN KEY ("variant_id") REFERENCES "variant"("variant_id");
ALTER TABLE "genotype" ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");

ALTER TABLE "locus_sequencing_stats" ADD FOREIGN KEY ("gene_db_crossref_id") REFERENCES dbxref(dbxref_id);
ALTER TABLE "locus_sequencing_stats" ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");

ALTER TABLE "summary_sequencing_stats" ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");

ALTER TABLE "taxonomy_stats"  ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");
ALTER TABLE "taxonomy_stats"  ADD FOREIGN KEY ("ncbi_taxon_id") REFERENCES taxon(ncbi_taxon_id);

ALTER TABLE "amplicon_target" ADD FOREIGN KEY ("gene_db_crossref_id") REFERENCES dbxref(dbxref_id);

ALTER TABLE "drug_synonym" ADD FOREIGN KEY ("drug_id") REFERENCES "drug"("drug_id");

ALTER TABLE "gene_drug_resistance_association" ADD FOREIGN KEY ("gene_db_crossref_id") REFERENCES dbxref(dbxref_id);
ALTER TABLE "gene_drug_resistance_association" ADD FOREIGN KEY ("drug_id") REFERENCES "drug"("drug_id");

ALTER TABLE "promoter_distance" ADD FOREIGN KEY ("gene_db_crossref_id") REFERENCES dbxref(dbxref_id);

ALTER TABLE "phenotypic_drug_susceptibility_test" ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");
ALTER TABLE "phenotypic_drug_susceptibility_test" ADD FOREIGN KEY ("drug_id") REFERENCES "drug"("drug_id");
ALTER TABLE "phenotypic_drug_susceptibility_test" ADD FOREIGN KEY ("medium_id") REFERENCES "growth_medium"("medium_id");
ALTER TABLE "phenotypic_drug_susceptibility_test" ADD FOREIGN KEY ("method_id") REFERENCES "phenotypic_drug_susceptibility_assessment_method"("method_id");

ALTER TABLE "minimum_inhibitory_concentration_test" ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");
ALTER TABLE "minimum_inhibitory_concentration_test" ADD FOREIGN KEY ("drug_id") REFERENCES "drug"("drug_id");

ALTER TABLE "smear_microscopy_results" ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");

ALTER TABLE "molecular_drug_resistance_test" ADD FOREIGN KEY ("sample_id") REFERENCES "sample"("sample_id");
ALTER TABLE "molecular_drug_resistance_test" ADD FOREIGN KEY ("drug_id") REFERENCES "drug"("drug_id");