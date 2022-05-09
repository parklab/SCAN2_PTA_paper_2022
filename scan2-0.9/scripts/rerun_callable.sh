#!/bin/bash
#SBATCH -p priopark
#SBATCH -A park_contrib
#SBATCH -t 12:00:00
#SBATCH --mem=4G

echo "remember to conda activate scansnv-pta"
echo "run this script from WITHIN the scansnv_fdr01_noX directory, so"
echo "that callable_regions and config.yaml are in the working dir."


if [ ! -f "config.yaml" ]; then
    echo "no config.yaml in this directory. are you running this script in the correct dir?"
    exit 1
fi

if [ ! -d "callable_regions" ]; then
    echo "no callable_regions in this directory. are you running this script in the correct dir?"
    exit 1
fi


snakemake \
    --snakefile /n/data1/hms/dbmi/park/jluquette/pta/scan2_used_here_old/SCAN2/scripts/Snakemake.callable_rerun \
    --configfile config.yaml \
    --jobs 1 \
    --rerun-incomplete \
    --keep-going \
    --latency-wait 60 \
    `#--drmaa ' -p park -A park_contrib --mem={resources.mem} -t 30:00 -o cluster-logs/slurm-%A.log'` \
    `#--unlock` \
    #--dryrun --quiet
