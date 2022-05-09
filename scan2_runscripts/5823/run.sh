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
    --bam 5823-tempmusc-1b1_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/5823-tempmusc-1b1_20170221-WGS-final-all.bam \
    --bulk-sample 5823-tempmusc-1b1_20170221-WGS \
    --bam 5823_20160810-dg-1cp2E1_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/5823_20160810-dg-1cp2E1_20170221-WGS-final-all.bam \
    --bam 5823_20160810-dg-1cp3D11_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/5823_20160810-dg-1cp3D11_20170221-WGS-final-all.bam \
    --bam 5823_20160810-dg-1cp3H1_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/5823_20160810-dg-1cp3H1_20170221-WGS-final-all.bam \
    --bam 5823_20160824-pfc-1cp1F11_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/5823_20160824-pfc-1cp1F11_20170221-WGS-final-all.bam \
    --sc-sample 5823_20160824-pfc-1cp1F11_20170221-WGS \
    --bam 5823_20160824-pfc-1cp2E1_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/5823_20160824-pfc-1cp2E1_20170221-WGS-final-all.bam \
    --sc-sample 5823_20160824-pfc-1cp2E1_20170221-WGS \
    --bam 5823_20160824-pfc-1cp2G5_20170221-WGS /n/data1/bch/genetics/walsh-park/data/Aging/Alignment/20170420_Alignment/5823_20160824-pfc-1cp2G5_20170221-WGS-final-all.bam \
    --sc-sample 5823_20160824-pfc-1cp2G5_20170221-WGS \
    --bam 5823PFC-A /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5823PFC-A.b37.bam \
    --sc-sample 5823PFC-A \
    --bam 5823PFC-B /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5823PFC-B.b37.bam \
    --sc-sample 5823PFC-B \
    --bam 5823PFC-C /n/data1/bch/genetics/walsh-park/PTA/NewData/.PreProcessing/5823PFC-C.b37.bam \
    --sc-sample 5823PFC-C \
    --regions-file /n/data1/hms/dbmi/park/jluquette/genotyper1/paper/scan-snv/regions.noX.txt \
    --output-dir scansnv_fdr01_noX \
    --abmodel-chunks 4 \
    --abmodel-samples-per-chunk 5000 \
    --target-fdr 0.01 \
    --joblimit 1000 \
    --drmaa ' -p park -A park_contrib --mem={resources.mem} -t 8:00:00 -o %logdir/slurm-%A.log' \
    --somatic-indel-pon /n/data1/hms/dbmi/park/jluquette/pta/panels/mda_pta_panel.partial_98pct.mmq60.indels.rda \
    --resume --somatic-indels #--snakemake-args ' --quiet --dryrun'
    #--resume --snakemake-args ' --quiet --touch --until abmodel_fit'

    #--resume --snakemake-args ' --until shapeit_gather --forcerun shapeit_prepare'
    #--resume --snakemake-args ' --restart-times 1' #' --quiet --dryrun'
