#!/usr/bin/env Rscript
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 1:00:00
#SBATCH --mem=10G

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@params['nbootstraps'],
        snakemake@input['mut'],
        snakemake@input['perm'],
        snakemake@output['full'],
        snakemake@output['summary'],
        snakemake@input['beds']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

suppressMessages(library(mutenrich))



results <- command.line.analysis(function(genome, bed.files) {
        # Allow for naming each BED feature by a ":" in front of the
        # BED file name.
        bed.and.feats <- lapply(1:length(bed.files), function(i) {
            bedf <- bed.files[i]
            bedf.elts <- strsplit(bedf, ":")[[1]]
            if (length(bedf.elts) == 2) {
                bed.feat.name <- bedf.elts[1]
                bedfile <- bedf.elts[2]
            } else {
                # use the filename as a feature ID
                bed.feat.name <- bedf.elts[1]
                bedfile <- bedf.elts[1]
            }
            c(bed.feat.name, bedfile)
        })
        gbed <- read.bed(bed.and.feats[[1]][2], genome,
            feature.name=bed.and.feats[[1]][1],
            add.chr.prefix=FALSE, remove.chr.prefix=FALSE, has.metaline=TRUE)
        if (length(bed.and.feats) > 1) {
            for (i in 2:length(bed.and.feats)) {
                gbed <- read.bed(bed.and.feats[[i]][2], genome,
                    granges=gbed,
                    feature.name=bed.and.feats[[i]][1],
                    add.chr.prefix=FALSE, remove.chr.prefix=FALSE, has.metaline=TRUE)
            }
        }
        print(gbed)
        cat('Signal metadata ------------- \n')
        print(attr(gbed, 'bed.metadata'))
        enrich.data(gbed=gbed)
    },
    genome="BSgenome.Hsapiens.UCSC.hg19",
    args=commandArgs(trailingOnly=TRUE))

str(results)

e <- results$e
es <- results$es
emeta <- attr(e$edata$gbed, 'bed.metadata')

cat("Saving output to", results$fulloutfile, "and", results$summaryoutfile, "\n")
save(e, es, emeta, file=results$fulloutfile, compress=FALSE)
save(es, emeta, file=results$summaryoutfile, compress=FALSE)
