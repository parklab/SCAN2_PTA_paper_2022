#!/bin/bash

echo "Converting somatic mutation RDA -> VCF"
./Convert_RDA_to_VCF.R SCAN2_PTA_sSNVs_filtered.rda SCAN2_PTA_sSNVs_filtered.vcf
./Convert_RDA_to_VCF.R SCAN2_PTA_sIndels_filtered.rda SCAN2_PTA_sIndels_filtered.vcf

echo "Running snpEff v4.3t"
java -jar snpEff/snpEff.jar -t -noStats -v hg19 SCAN2_PTA_sSNVs_filtered.vcf > SCAN2_PTA_sSNVs_filtered.snpEff.vcf

java -jar snpEff/snpEff.jar -t -noStats -v hg19 SCAN2_PTA_sIndels_filtered.vcf > SCAN2_PTA_sIndels_filtered.snpEff.vcf
