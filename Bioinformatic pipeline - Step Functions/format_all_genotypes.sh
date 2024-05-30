#/bin/bash

sample_name=$(basename $1)
reference_name=$2

for genotyper in freebayes gatk bcftools
do
    bcftools norm --fasta-ref $reference_name --multiallelics - $1-${genotyper}.vcf | \
      bcftools query -f "${genotyper},[%SAMPLE],%CHROM,%POS,%ID,%REF,%ALT,%QUAL,[%AD],[%DP],[%GT]\n" | \
        sed "s/,\.,/,,/g" >> \
          genotype/${sample_name}-genotype.csv
done

gzip genotype/${sample_name}-genotype.csv