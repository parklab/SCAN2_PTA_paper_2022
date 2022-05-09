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
    --regions-file /n/data1/hms/dbmi/park/jluquette/walsh-1465/regions_1-22_1166windows.txt \
    --bam 1465-cortex_1-neuron_MDA_12 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_12.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_12 \
    --bam 1465-cortex_1-neuron_MDA_18 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_18.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_18 \
    --bam 1465-cortex_1-neuron_MDA_20 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_20.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_20 \
    --bam 1465-cortex_1-neuron_MDA_24 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_24.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_24 \
    --bam 1465-cortex_1-neuron_MDA_25 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_25.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_25 \
    --bam 1465-cortex_1-neuron_MDA_2_WGSb /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_2_WGSb.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_2_WGSb \
    --bam 1465-cortex_1-neuron_MDA_30 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_30.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_30 \
    --bam 1465-cortex_1-neuron_MDA_39 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_39.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_39 \
    --bam 1465-cortex_1-neuron_MDA_3_WGSb /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_3_WGSb.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_3_WGSb \
    --bam 1465-cortex_1-neuron_MDA_43 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_43.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_43 \
    --bam 1465-cortex_1-neuron_MDA_46 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_46.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_46 \
    --bam 1465-cortex_1-neuron_MDA_47 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_47.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_47 \
    --bam 1465-cortex_1-neuron_MDA_5 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_5.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_5 \
    --bam 1465-cortex_1-neuron_MDA_51_WGSb /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_51_WGSb.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_51_WGSb \
    --bam 1465-cortex_1-neuron_MDA_6_WGSb /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_6_WGSb.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_6_WGSb \
    --bam 1465-cortex_1-neuron_MDA_8 /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_1-neuron_MDA_8.bam \
    --sc-sample 1465-cortex_1-neuron_MDA_8 \
    --bam 1465-cortex_BulkDNA_WGSb /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-cortex_BulkDNA_WGSb.bam \
    --bulk-sample 1465-cortex_BulkDNA_WGSb \
    --bam 1465-heart_BulkDNA_WGSb /n/data1/bch/genetics/walsh-park/data/NormalSingleCell/Alignment/20160404_NYgenome_realignment/1465/1465-heart_BulkDNA_WGSb.bam \
    --bam 1465_ct_8p2h8 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/1465_ct_8p2h8.bam \
    --sc-sample 1465_ct_8p2h8 \
    --bam 1465_ct_p2B11 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/1465_ct_p2B11.bam \
    --sc-sample 1465_ct_p2B11 \
    --bam 1465_ctx_p2F06 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/1465_ctx_p2F06.bam \
    --sc-sample 1465_ctx_p2F06 \
    --bam 1465_ctx_p2g8 /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20161223_Alignment/1465_ctx_p2g8.bam \
    --sc-sample 1465_ctx_p2g8 \
    --bam 1465BA9-A /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/1465BA9-A.b37.bam \
    --sc-sample 1465BA9-A \
    --bam 1465BA9-B /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/1465BA9-B.b37.bam \
    --sc-sample 1465BA9-B \
    --bam 1465BA9-C /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/1465BA9-C.b37.bam \
    --sc-sample 1465BA9-C \
    --bam 1465BA9-D /n/data1/bch/genetics/walsh-park/PTA/.PreProcessing/1465BA9-D.b37.bam \
    --sc-sample 1465BA9-D \
    --output-dir scansnv_fdr01_noX \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --joblimit 1000 \
    --add-chr-prefix \
    --drmaa ' -p park -A park_contrib --mem={resources.mem} -t 8:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon panel_52_pta_neurons_042521.partial_100pct.mmq60.indels.rda \
    --resume --somatic-indels
