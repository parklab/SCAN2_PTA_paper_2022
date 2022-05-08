# Large result sets needed for analyses
When possible, results are made available on Github. However, some results
are either too large to host or contain potentially identifiable information
and cannot be made publicly available.

To run the scripts in this directory, the user must populate the following
directories with large data files with specific naming conventions.

## simple_chrX_genotyping
Simple chrX genotyping uses a matrix of ref/alt supporting reads from GATK
HaplotypeCaller. 

## scan2_somatic_genotype_rdas
SCAN2 v0.9 produces R data files (.rda) in snv/[sample_ID]/somatic_genotypes.rda
and indel/[sample_ID]/somatic_genotypes.pon_filter.rda. These files must be made
deposited in scan2_somatic_genotype_rdas/{sSNVs,sIndels}/[sample_ID].rda.

