ALTER TABLE "country" ADD PRIMARY KEY ("country_id");

ALTER TABLE "sample" ADD PRIMARY KEY ("sample_id");

ALTER TABLE "bioproject" ADD PRIMARY KEY ("bioproject_id");

ALTER TABLE "dataset" ADD PRIMARY KEY("dataset_id");

ALTER TABLE "contributor" ADD PRIMARY KEY("contributor_id");

ALTER TABLE "patient" ADD PRIMARY KEY ("patient_id");

ALTER TABLE "sequencing_data" ADD PRIMARY KEY ("sequencing_data_id");

ALTER TABLE "variant" ADD PRIMARY KEY ("variant_id");

ALTER TABLE "annotation" ADD PRIMARY KEY ("annotation_id");

ALTER TABLE "genotype" ADD PRIMARY KEY ("genotype_id");

ALTER TABLE "amplicon_target" ADD PRIMARY KEY ("amplicon_target_id");

ALTER TABLE "growth_medium" ADD PRIMARY KEY ("medium_id");

ALTER TABLE "phenotypic_drug_susceptibility_assessment_method" ADD PRIMARY KEY ("method_id");

ALTER TABLE "drug" ADD PRIMARY KEY ("drug_id");

ALTER TABLE "phenotypic_drug_susceptibility_test" ADD PRIMARY KEY ("test_id");

ALTER TABLE "minimum_inhibitory_concentration_test" ADD PRIMARY KEY ("test_id");

ALTER TABLE "smear_microscopy_results" ADD PRIMARY KEY ("test_id");

ALTER TABLE "molecular_drug_resistance_test" ADD PRIMARY KEY ("test_id");