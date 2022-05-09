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
    --bam 5559-pfc1C4 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5559_0902-pfc-1cp1C4.bam \
    --sc-sample 5559-pfc1C4 \
    --bam 5559-pfc1C7 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5559_0902-pfc-1cp1C7.bam \
    --sc-sample 5559-pfc1C7 \
    --bam 5559-pfc1E2 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5559_0902-pfc-1cp1E2.bam \
    --sc-sample 5559-pfc1E2 \
    --bam 5559-pfc1H2 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5559_0902-pfc-1cp1H2.bam \
    --sc-sample 5559-pfc1H2 \
    --bam 5559-pfc2A3 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5559_0902-pfc-1cp2A3.bam \
    --sc-sample 5559-pfc2A3 \
    --bam 5559-bulk /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20151226_Macrogen/5559_0928-hrt-1b1.bam \
    --bulk-sample 5559-bulk \
    --bam 5559PFC-A /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5559PFC-A.b37.bam \
    --sc-sample 5559PFC-A \
    --bam 5559PFC-B /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5559PFC-B.b37.bam \
    --sc-sample 5559PFC-B \
    --bam 5559PFC-C /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5559PFC-C.b37.bam \
    --sc-sample 5559PFC-C \
    --regions-file /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/scan-snv/regions.noX.txt \
    --output-dir scansnv_fdr01_noX \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --joblimit 1000 \
    --drmaa ' -p park -A park_contrib --mem={resources.mem} -t 8:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.tab \
    --resume --somatic-indels
