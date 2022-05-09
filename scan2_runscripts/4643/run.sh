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
    --bam 4643-Neuron-3 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/4643-Neuron-3.b37.bam \
    --sc-sample 4643-Neuron-3 \
    --bam 4643-Neuron-4 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/4643-Neuron-4.b37.bam \
    --sc-sample 4643-Neuron-4 \
    --bam 4643-Neuron-6 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/4643-Neuron-6.b37.bam \
    --sc-sample 4643-Neuron-6 \
    --bam 4643_Bulk-Liver /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_Bulk-Liver.bam \
    --bulk-sample 4643_Bulk-Liver \
    --bam 4643_MDA_1 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_MDA_1.bam \
    --sc-sample 4643_MDA_1 \
    --bam 4643_MDA_2 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_MDA_2.bam \
    --sc-sample 4643_MDA_2 \
    --bam 4643-MDA_23 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_MDA_23.bam \
    --sc-sample 4643-MDA_23 \
    --bam 4643_MDA_24 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_MDA_24.bam \
    --sc-sample 4643_MDA_24 \
    --bam 4643_MDA_26 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_MDA_26.bam \
    --sc-sample 4643_MDA_26 \
    --bam 4643_MDA_3 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_MDA_3.bam \
    --sc-sample 4643_MDA_3 \
    --bam 4643_MDA_31 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_MDA_31.bam \
    --sc-sample 4643_MDA_31 \
    --bam 4643_MDA_32 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_MDA_32.bam \
    --sc-sample 4643_MDA_32 \
    --bam 4643_MDA_4 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_MDA_4.bam \
    --sc-sample 4643_MDA_4 \
    --bam 4643_MDA_5 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4643/4643_MDA_5.bam \
    --sc-sample 4643_MDA_5 \
    --output-dir scansnv \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --joblimit 1000 \
    --drmaa ' -p short --mem={resources.mem} -t 12:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.rda \
    --resume --somatic-indels
