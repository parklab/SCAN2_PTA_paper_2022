#!/usr/bin/env Rscript
#SBATCH -p park
#SBATCH -A park_contrib
#SBATCH -t 1:00:00
#SBATCH --mem=16G

# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(c(
        snakemake@input['gtf'],
        snakemake@params['tissue'],
        snakemake@input['gct'],
        snakemake@output['bigwig']
    ))
    cat('Got command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- commandArgs(trailingOnly=TRUE)

if (length(args) != 4) {
    stop('usage: Convert_GTEx_expression_to_BigWig.R gencode_gene_model.gtf tissue_with_spaces tissue_median_tpm.gct out.bigwig')
}

gene.model.file <- args[1]
tissue <- args[2]
tissue.median.tpm.file <- args[3]
out.bigwig <- args[4]

tile.size=50

if (file.exists(out.bigwig))
    stop(paste('output file', out.bigwig, 'already exists, please delete it first'))

tissue.with.spaces <- tissue
tissue <- gsub(' ', '_', tissue)
tissue <- gsub('(', '_', tissue, fixed=TRUE)
tissue <- gsub(')', '_', tissue, fixed=TRUE)
cat("Mapping tissue", tissue.with.spaces, "to", tissue, "\n")

suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(data.table))
suppressMessages(library(BSgenome.Hsapiens.UCSC.hg19))


cat('Loading gene model:', gene.model.file, '\n')
# this version is the result of running GTEx's collapse script
gtex <- import.gff(gene.model.file)
# ENSEMBL IDs have underscores after them in GENCODE but not in GTEx
gtex$gene_id2 <- sapply(strsplit(gtex$gene_id,'_'),head,1)


cat('Annotating tissue median TPM expression levels:', tissue.median.tpm.file, '\n')
gtex.tpm <- fread(tissue.median.tpm.file, skip=2)
setkey(gtex.tpm, Name)
# making the entire pipeline safe for file names with spaces/parens is WAY too much work
colnames(gtex.tpm) <- gsub(' ', '_', colnames(gtex.tpm))
colnames(gtex.tpm) <- gsub('(', '_', colnames(gtex.tpm), fixed=TRUE)
colnames(gtex.tpm) <- gsub(')', '_', colnames(gtex.tpm), fixed=TRUE)
gtex.tpm <- gtex.tpm[gtex$gene_id2]
gtex$mean.expr <- rowMeans(gtex.tpm[,-(1:2)])


cat(paste0('Using tissue=', tissue, ' to compute max expression signal along genome\n'))
e <- gtex.tpm[[tissue]]

selected.genes <- !is.na(gtex$mean.expr) & apply(gtex.tpm[,-(1:2)],1,max) > 0

cat(sprintf('    %d / %d genes with non-NA expression and median expression > 0 in at least 1 of %d tissues\n',
    sum(selected.genes), length(selected.genes), ncol(gtex.tpm)-2))

gtex <- gtex[selected.genes,]
e <- e[selected.genes]

# Only retain primary contigs
gtex <- gtex[seqnames(gtex) %in% seqlevels(gtex)[1:25],]
seqlevels(gtex) <- seqlevels(gtex)[1:25]
seqinfo(gtex) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)
cat('Tiling hg19 with', tile.size,  'bp bins\n')
bins <- tileGenome(seqlengths=seqlengths(gtex), tilewidth=tile.size, cut.last.tile.in.chrom=T)

cat(paste0('Mapping GTEx to non-overlapping ', tile.size, ' bp bins and taking max expression over tissue=', tissue, ' gene expression\n'))
ols <- findOverlaps(bins, gtex)
system.time(maxexpr <- sapply(split(e[to(ols)], from(ols)), max))
emap <- data.frame(from=as.integer(names(maxexpr)), to=maxexpr)


# Make a fresh tiling GRanges with no metadata and attach the expression values for each tile.
new.bins <- GRanges(seqnames=seqnames(bins), ranges=ranges(bins), seqinfo=seqinfo(bins))
new.bins$score <- NA
new.bins$score[emap$from] <- emap$to
new.bins <- new.bins[!is.na(new.bins$score),]

cat("Writing bigwig ", out.bigwig, '\n')
export.bw(new.bins, con=out.bigwig)

if ('snakemake' %in% ls()) {
    sink()
}
