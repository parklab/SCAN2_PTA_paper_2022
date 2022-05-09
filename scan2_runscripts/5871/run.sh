#!/bin/bash

#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 120:00:00
#SBATCH --mem 4000

/n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/bin/scansnv \
    --snakefile /n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/snakemake/Snakefile \
    --scripts /n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/scripts \
    --regions-file /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/scan-snv/regions.noX.txt \
    --ref /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/resources/human_g1k_v37_decoy.fasta \
    --dbsnp /home/ljl11/balance/resources/dbsnp_147_b37_common_all_20160601.vcf \
    --shapeit-panel /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/resources/1000GP_Phase3 \
    --bam 5871-BLK-liver /n/data1/bch/genetics/walsh-park/data/PTA/sourceData/5871blkliver/.PreProcessing/5871-BLK-liver.b37.bam \
    --bulk-sample 5871-BLK-liver \
    --bam 5871-Neuron-4 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/5871-Neuron-4.b37.bam \
    --sc-sample 5871-Neuron-4 \
    --bam 5871-Neuron-5 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/5871-Neuron-5.b37.bam \
    --sc-sample 5871-Neuron-5 \
    --bam 5871-Neuron-6 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/5871-Neuron-6.b37.bam \
    --sc-sample 5871-Neuron-6 \
    --output-dir scansnv \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --joblimit 5000 \
    --drmaa ' -p short --mem={resources.mem} -t 12:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.tab \
    --somatic-indels \
    --resume
