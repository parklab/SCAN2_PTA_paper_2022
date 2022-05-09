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
    --bam 5219-Neuron-2 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/5219-Neuron-2.b37.bam \
    --sc-sample 5219-Neuron-2 \
    --bam 5219-Neuron-4 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/5219-Neuron-4.b37.bam \
    --sc-sample 5219-Neuron-4 \
    --bam 5219-Neuron-5 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/5219-Neuron-5.b37.bam \
    --sc-sample 5219-Neuron-5 \
    --bam 5219_cb_bulk /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/5219_cb_bulk.bam \
    --bulk-sample 5219_cb_bulk \
    --bam 5219_ct_p1G1 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/5219_ct_p1G1.bam \
    --sc-sample 5219_ct_p1G1 \
    --bam 5219_ct_p1G7 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/5219_ct_p1G7.bam \
    --sc-sample 5219_ct_p1G7 \
    --bam 5219_ct_p2A12 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/5219_ct_p2A12.bam \
    --sc-sample 5219_ct_p2A12 \
    --bam 5219_ct_p2C3 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/5219_ct_p2C3.bam \
    --sc-sample 5219_ct_p2C3 \
    --output-dir scansnv \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --joblimit 2000 \
    --drmaa ' -p short --mem={resources.mem} -t 12:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.tab \
    --somatic-indels \
    --resume
