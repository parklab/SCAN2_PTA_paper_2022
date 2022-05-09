#!/bin/bash

#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem 4000

/n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/bin/scansnv \
    --ref /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/resources/human_g1k_v37_decoy.fasta \
    --dbsnp /home/ljl11/balance/resources/dbsnp_147_b37_common_all_20160601.vcf \
    --shapeit-panel /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/resources/1000GP_Phase3 \
    --regions-file /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/scan-snv/regions.noX.txt \
    --snakefile /n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/snakemake/Snakefile \
    --scripts /n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/scripts \
    --bam 5657-bulk /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5657_0717-hrt-1b1.bam \
    --bulk-sample 5657-bulk \
    --bam 5657-pfc1D2 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5657_0902-pfc-1cp1D2.bam \
    --sc-sample 5657-pfc1D2 \
    --bam 5657-pfc1E11 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5657_0902-pfc-1cp1E11.bam \
    --sc-sample 5657-pfc1E11 \
    --bam 5657-pfc2A6 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5657_0902-pfc-1cp2A6.bam \
    --sc-sample 5657-pfc2A6 \
    --bam 5657-pfc2F1 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5657_0902-pfc-1cp2F1.bam \
    --sc-sample 5657-pfc2F1 \
    --bam 5657-pfc2G9 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5657_0902-pfc-1cp2G9.bam \
    --sc-sample 5657-pfc2G9 \
    --bam 5657PFC-A /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5657PFC-A.b37.bam \
    --sc-sample 5657PFC-A \
    --bam 5657PFC-B /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5657PFC-B.b37.bam \
    --sc-sample 5657PFC-B \
    --bam 5657PFC-C /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5657PFC-C.b37.bam \
    --sc-sample 5657PFC-C \
    --output-dir scansnv_fdr01_noX \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --joblimit 1000 \
    --drmaa ' -p park -A park_contrib --mem={resources.mem} -t 8:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.tab \
    --resume --somatic-indels
