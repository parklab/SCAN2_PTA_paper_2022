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
    --bam 5087-hrt-1b1 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20160518_Macrogen/5087-hrt-1b1.bam \
    --bulk-sample 5087-hrt-1b1 \
    --bam 5087pfc-Lp1C5 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20170213_Alignment/5087pfc-Lp1C5.bam \
    --sc-sample 5087pfc-Lp1C5 \
    --bam 5087pfc-Rp1G4 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20170213_Alignment/5087pfc-Rp1G4.bam \
    --sc-sample 5087pfc-Rp1G4 \
    --bam 5087pfc-Rp3C5 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20170213_Alignment/5087pfc-Rp3C5.bam \
    --sc-sample 5087pfc-Rp3C5 \
    --bam 5087pfc-Rp3F4 /n/data1/bch/genetics/walsh-park/data/Hippocampus/Alignment/20170213_Alignment/5087pfc-Rp3F4.bam \
    --sc-sample 5087pfc-Rp3F4 \
    --bam 5087PFC-A /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5087PFC-A.b37.bam \
    --sc-sample 5087PFC-A \
    --bam 5087PFC-B /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5087PFC-B.b37.bam \
    --sc-sample 5087PFC-B \
    --bam 5087PFC-C /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5087PFC-C.b37.bam \
    --sc-sample 5087PFC-C \
    --regions-file /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/scan-snv/regions.noX.txt \
    --output-dir scansnv_fdr01_noX \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --joblimit 1000 \
    --drmaa ' -p park -A park_contrib --mem={resources.mem} -t 8:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.tab \
    --resume --somatic-indels
