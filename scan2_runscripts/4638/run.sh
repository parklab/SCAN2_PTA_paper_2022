#!/bin/bash

#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 240:00:00
#SBATCH --mem 4000

/n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/bin/scansnv \
    --ref /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/resources/human_g1k_v37_decoy.fasta \
    --dbsnp /home/ljl11/balance/resources/dbsnp_147_b37_common_all_20160601.vcf \
    --snakefile /n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/snakemake/Snakefile \
    --scripts /n/data1/hms/dbmi/park/jluquette/pta/scan-snv2/scripts \
    --shapeit-panel /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/resources/1000GP_Phase3 \
    --regions-file /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/scan-snv/regions.noX.txt \
    --bam 4638-Bulk-Heart /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-Bulk-Heart.bam \
    --bulk-sample 4638-Bulk-Heart \
    --bam 4638-MDA-2 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-02.bam \
    --sc-sample 4638-MDA-2 \
    --bam 4638-MDA-03 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-03.bam \
    --sc-sample 4638-MDA-03 \
    --bam 4638-MDA-4 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-04.bam \
    --sc-sample 4638-MDA-4 \
    --bam 4638-MDA-7 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-07.bam \
    --sc-sample 4638-MDA-7 \
    --bam 4638-MDA-11 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-11.bam \
    --sc-sample 4638-MDA-11 \
    --bam 4638-MDA-12 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-12.bam \
    --sc-sample 4638-MDA-12 \
    --bam 4638-MDA-13 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-13.bam \
    --sc-sample 4638-MDA-13 \
    --bam 4638-MDA-14 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-14.bam \
    --sc-sample 4638-MDA-14 \
    --bam 4638-MDA-15 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-15.bam \
    --sc-sample 4638-MDA-15 \
    --bam 4638-MDA-20 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-20.bam \
    --sc-sample 4638-MDA-20 \
    --bam 4638-MDA-24 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/4638/4638-MDA-24.bam \
    --sc-sample 4638-MDA-24 \
    --bam 4638-Neuron-4 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/4638-Neuron-4.b37.bam \
    --sc-sample 4638-Neuron-4 \
    --bam 4638-Neuron-5 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/4638-Neuron-5.b37.bam \
    --sc-sample 4638-Neuron-5 \
    --bam 4638-Neuron-6 /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/4638-Neuron-6.b37.bam \
    --sc-sample 4638-Neuron-6 \
    --output-dir scansnv \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --joblimit 1000 \
    --drmaa ' -p short --mem={resources.mem} -t 12:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.rda \
    --resume --somatic-indels
