#!/usr/bin/env Rscript

library(scan2)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 2)
    stop("usage: Convert_RDA_to_VCF.R input.rda output.vcf")

in.rda <- args[1]
out.vcf <- args[2]

muts <- get(load(in.rda))
muts$pass <- TRUE  # Write out everything

scansnv.df.to.vcf(df=muts, out.file=out.vcf, sample.name='dummy', chrs=1:22)
