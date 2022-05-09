#!/bin/bash

#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem 4000

/n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/bin/scansnv \
    --ref /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/resources/human_g1k_v37_decoy.fasta \
    --dbsnp /home/ljl11/balance/resources/dbsnp_147_b37_common_all_20160601.vcf \
    --snakefile /n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/snakemake/Snakefile \
    --scripts /n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/scripts \
    --bam 1278_ct_p1E3 1278_ct_p1E3.bam \
    --sc-sample 1278_ct_p1E3 \
    --bam 1278_ct_p1E6 1278_ct_p1E6.bam \
    --sc-sample 1278_ct_p1E6 \
    --bam 1278_ct_p1G9 1278_ct_p1G9.bam \
    --sc-sample 1278_ct_p1G9 \
    --bam 1278_ct_p2B9 1278_ct_p2B9.bam \
    --sc-sample 1278_ct_p2B9 \
    --bam 1278_ct_p2C7 1278_ct_p2C7.bam \
    --sc-sample 1278_ct_p2C7 \
    --bam 1278_ct_p2E4 1278_ct_p2E4.bam \
    --sc-sample 1278_ct_p2E4 \
    --bam 1278_ct_p2E6 1278_ct_p2E6.bam \
    --sc-sample 1278_ct_p2E6 \
    --bam 1278_ct_p2F5 1278_ct_p2F5.bam \
    --sc-sample 1278_ct_p2F5 \
    --bam 1278_ct_p2G5 1278_ct_p2G5.bam \
    --sc-sample 1278_ct_p2G5 \
    --bam 1278_heart_bulk 1278_heart_bulk.bam \
    --bulk-sample 1278_heart_bulk \
    --bam 1278BA9-A 1278BA9-A.b37.bam \
    --sc-sample 1278BA9-A  \
    --bam 1278BA9-B 1278BA9-B.b37.bam \
    --sc-sample 1278BA9-B \
    --bam 1278BA9-C 1278BA9-C.b37.bam \
    --sc-sample 1278BA9-C \
    --shapeit-panel /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/resources/1000GP_Phase3 \
    --regions-file /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/scan-snv/regions.noX.txt \
    --output-dir scansnv_fdr01_noX \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --joblimit 1000 \
    --drmaa ' -p park -A park_contrib --mem={resources.mem} -t 8:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.rda \
    --resume --somatic-indels
