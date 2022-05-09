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
    --shapeit-panel /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/resources/1000GP_Phase3 \
    --regions-file /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/scan-snv/regions.noX.txt \
    --bam 5817_ct_p1H10 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/5817_ct_p1H10.bam \
    --sc-sample 5817_ct_p1H10 \
    --bam 5817_ct_p1H2 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/5817_ct_p1H2.bam \
    --sc-sample 5817_ct_p1H2 \
    --bam 5817_ct_p1H5 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/5817_ct_p1H5.bam \
    --sc-sample 5817_ct_p1H5 \
    --bam 5817_ct_p2H6 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/5817_ct_p2H6.bam \
    --sc-sample 5817_ct_p2H6 \
    --bam 5817_liver_bulk /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/5817_liver_bulk.bam \
    --bulk-sample 5817_liver_bulk \
    --bam 5817PFC-A /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5817PFC-A.b37.bam \
    --sc-sample 5817PFC-A \
    --bam 5817PFC-B /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5817PFC-B.b37.bam \
    --sc-sample 5817PFC-B \
    --bam 5817PFC-C /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5817PFC-C.b37.bam \
    --sc-sample 5817PFC-C \
    --output-dir scansnv_fdr01_noX \
    --target-fdr 0.01 \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --joblimit 1000 \
    --drmaa ' -p park -A park_contrib --mem={resources.mem} -t 8:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.tab \
    --resume --somatic-indels
