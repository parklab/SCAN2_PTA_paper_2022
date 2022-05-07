# The following commands will make a conda environment suitable
# for installing the R SCAN2 0.9 package and running the scripts
# in this directory.
#
# The preferred way to create this environment is to install SCAN2
# as described on https://github.com/parklab/SCAN2.

# Create just a python 3.8 environment with the mamba package manager
conda create -n scan2_pta_paper_scripts -c conda-forge -c bioconda mamba python=3.8
conda activate scan2_pta_paper_scripts


mamba install -c conda-forge -c bioconda r-base r-devtools r-pracma r-r.utils r-argparse r-lme4 r-lmertest r-data.table r-reticulate r-plyr bioconductor-bsgenome.hsapiens.1000genomes.hs37d5


echo "You must now following the SigProfilerMatrixGenerator installation"
echo "instructions at https://github.com/parklab/SCAN2."
