#!/bin/sh

# Author: Maxwell Sherman (with input from Daniel Kwon)
# Last revision: 03/04/2016

## USAGE:
## $1: input bam file
## $2: output prefix

echo "This script requires CNV scripts from Baslan et al. Genome-wide copy number analysis of single cells. 2012"
echo "Install these scripts and set the path below to your local installation. This should contain the scripts"
echo "cbs.r, copynumber.r, cal_mapd.r, etc."

path_to_baslancnv=/home/mas138/orchestra/sandbox/Install/baslancnv/

samtools view $1 | python2 $path_to_baslancnv/varbin.50k.bam.py $2.varbin.txt $path_to_baslancnv/hg19.bin.boundaries.50k.bowtie.k50.sorted.txt $path_to_baslancnv/hg19.chrom.sizes.txt
Rscript $path_to_baslancnv/cbs.r $2.varbin.txt $2.data.txt $2.data.short.txt $path_to_baslancnv/hg19.50k.k50.bad.bins.txt $path_to_baslancnv/hg19.varbin.gc.content.50k.bowtie.k50.txt
Rscript $path_to_baslancnv/copynumber.r $2.data.txt $2.data.short.txt $2.copynumber.txt > $2.copynumber.log
Rscript $path_to_baslancnv/cal_mapd.r $2.copynumber.txt > $2.mapd
