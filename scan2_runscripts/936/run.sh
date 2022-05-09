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
    --bam 936-hrt-1b1_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/936-hrt-1b1_20170221-WGS-final-all.bam \
    --bulk-sample 936-hrt-1b1_20170221-WGS \
    --bam 936_20141001-pfc-1cp1G11_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/936_20141001-pfc-1cp1G11_20170221-WGS-final-all.bam \
    --sc-sample 936_20141001-pfc-1cp1G11_20170221-WGS \
    --bam 936_20141001-pfc-1cp1H9_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/936_20141001-pfc-1cp1H9_20170221-WGS-final-all.bam \
    --sc-sample 936_20141001-pfc-1cp1H9_20170221-WGS \
    --bam 936_20141001-pfc-1cp2F6_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/936_20141001-pfc-1cp2F6_20170221-WGS-final-all.bam \
    --sc-sample 936_20141001-pfc-1cp2F6_20170221-WGS \
    --bam 936PFC-A /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/936PFC-A.b37.bam \
    --sc-sample 936PFC-A \
    --bam 936PFC-B /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/936PFC-B.b37.bam \
    --sc-sample 936PFC-B \
    --bam 936PFC-C /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/936PFC-C.b37.bam \
    --sc-sample 936PFC-C \
    --output-dir scansnv_fdr01_noX \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --add-chr-prefix \
    --joblimit 1000 \
    --drmaa ' -p park -A park_contrib --mem={resources.mem} -t 5:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indels \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.tab \
    --resume
