CREATE INDEX "sequencing_data_sample_id_idx" ON "sequencing_data"("sample_id");

CREATE INDEX "additional_sample_name_sample_id_idx" ON "additional_sample_name"("sample_id");

CREATE INDEX "variant_pos_idx" ON "variant"("position");

CREATE INDEX "annnotation_ref_id_idx" ON "annotation"("reference_db_crossref_id");
CREATE INDEX "annnotation_value_idx" ON "annotation"("hgvs_value");

CREATE INDEX "tiered_variant_categories_var_id" ON "tiered_variant_categories"("variant_id");
CREATE INDEX "tiered_variant_categories_gen_id" ON "tiered_variant_categories"("gene_db_crossref_id");

CREATE INDEX "variant_to_annotation_variant_id_idx" ON "variant_to_annotation"("variant_id");
CREATE INDEX "variant_to_annotation_annotation_id_idx" ON "variant_to_annotation"("annotation_id");

CREATE INDEX "genotype_variant_id_idx" ON "genotype"("variant_id");
CREATE INDEX "genotype_sample_id_idx" ON "genotype"("sample_id");

CREATE INDEX "locus_sequencing_stats_sample_id_idx" ON "locus_sequencing_stats"("sample_id");
CREATE INDEX "locus_sequencing_stats_gene_db_crossref_id_idx" ON "locus_sequencing_stats"("gene_db_crossref_id");

CREATE INDEX "taxonomy_stats_sample_id_idx" ON "taxonomy_stats"("sample_id");

CREATE INDEX "summary_sequencing_stats_sample_id_idx" ON "summary_sequencing_stats"("sample_id");

CREATE INDEX "phenotypic_drug_susceptibility_test_sample_id_idx" ON "phenotypic_drug_susceptibility_test"("sample_id");
CREATE INDEX "phenotypic_drug_susceptibility_test_drug_id_idx" ON "phenotypic_drug_susceptibility_test"("drug_id");

CREATE INDEX "minimum_inhibitory_concentration_test_sample_id_idx" ON "minimum_inhibitory_concentration_test"("sample_id");
CREATE INDEX "minimum_inhibitory_concentration_test_drug_id_idx" ON "minimum_inhibitory_concentration_test"("drug_id");

CREATE INDEX "formatted_annotation_per_gene_variant_id_idx" ON "formatted_annotation_per_gene"("variant_id");
CREATE INDEX "formatted_annotation_per_gene_gene_id_idx" ON "formatted_annotation_per_gene"("gene_id");

CREATE INDEX "ranked_annotation_variant_id_idx" ON "ranked_annotation"("variant_id");
CREATE INDEX "ranked_annotation_annotation_id_idx" ON "ranked_annotation"("annotation_id");
CREATE INDEX "ranked_annotation_gene_name_idx" ON "ranked_annotation"("gene_name");

CREATE INDEX "preferred_annotation_variant_id_idx" ON "preferred_annotation"("variant_id");
CREATE INDEX "preferred_annotation_gene_name_idx" ON "preferred_annotation"("gene_name");