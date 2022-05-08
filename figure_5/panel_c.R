#!/usr/bin/env Rscript
library(data.table)

load('../metadata.rda')

snv.burdens <- fread('../figure_4/Collected_SCAN2_sSNV_burdens.csv')
total.snv.burden <- sum(snv.burdens[amp=='PTA']$burden)
snv.calls <- fread('../data/Collected_SCANSNV_SCAN2_sSNV_calls.csv.gz')
total.snv.calls <- sum(snv.calls[sample %in% meta[meta$amp=='PTA',]$sample, passA|passB])
est.snv.sensitivity <- total.snv.calls / total.snv.burden
snv.snpeff <- fread('SCAN2_PTA_sSNVs_filtered.snpEff.vcf')[[8]]
snv.high <- length(grep('|HIGH|', snv.snpeff, fixed=TRUE))

indel.burdens <- fread('../figure_4/Collected_SCAN2_sIndel_burdens.csv')
total.indel.burden <- sum(indel.burdens[amp=='PTA']$burden)
indel.calls <- fread('../data/Collected_SCAN2_sIndel_calls.csv.gz')
total.indel.calls <- sum(indel.calls[sample %in% meta[meta$amp=='PTA',]$sample, passA|passB])
est.indel.sensitivity <- total.indel.calls / total.indel.burden
indel.snpeff <- fread('SCAN2_PTA_sIndels_filtered.snpEff.vcf')[[8]]
indel.high <- length(grep('|HIGH|', indel.snpeff, fixed=TRUE))

cat(sprintf('sSNVs: calls=%d, total burden=%0.0f, resulting sensitivity=%0.2f, HIGH impact=%d\n',
    total.snv.calls, total.snv.burden, est.snv.sensitivity, snv.high))
cat(sprintf('sIndels: calls=%d, total burden=%0.0f, resulting sensitivity=%0.2f, HIGH impact=%d\n',
    total.indel.calls, total.indel.burden, est.indel.sensitivity, indel.high))


pdf(width=0.8*6/7, height=2.0*6/7, pointsize=5, file='panel_c.pdf')
par(mar=c(5,4,1,1))
barplot(matrix(c(indel.high, indel.high*(1/est.indel.sensitivity - 1),
                 snv.high, snv.high*(1/est.snv.sensitivity - 1)), ncol=2),
    beside=FALSE, border=NA, col=c('#4c4c4c', '#c6c6c6'), bty='n', yaxt='n', las=3,
    ylab='High impact somatic mutations', names=c('Indels', 'SNVs'))
axis(side=2, at=0:9*10, lwd=0.35, lwd.ticks=0.35)
dev.off()
