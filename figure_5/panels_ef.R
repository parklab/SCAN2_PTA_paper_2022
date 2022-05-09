#!/usr/bin/env Rscript
library(data.table)

meta <- fread('../external_data/Roadmap_Epigenomics/EID_metadata_with_H3K27ac.txt')
setkey(meta, EID)


# Add 1e-4 to p-value since 0 out of 10,000 permutations is better
# interpreted as p < 1e-4 rather than p=0.
snv.distal <- sapply(meta$EID, function(eid) {
    load(sprintf('sSNVs_vs_%s_H3K27ac_peaks.SUMMARY.rda', eid))
    ret <- unlist(es[[1]]['distal',c('enr', 'pval')])
    c(ret[1], -log10(ret[2] + 1e-4))
})
snv.proximal <- sapply(meta$EID, function(eid) {
    load(sprintf('sSNVs_vs_%s_H3K27ac_peaks.SUMMARY.rda', eid))
    ret <- unlist(es[[1]]['TSS_proximal',c('enr', 'pval')])
    c(ret[1], -log10(ret[2] + 1e-4))
})

indel.distal <- sapply(meta$EID, function(eid) {
    load(sprintf('sIndels_vs_%s_H3K27ac_peaks.SUMMARY.rda', eid))
    ret <- unlist(es[[1]]['distal',c('enr', 'pval')])
    c(ret[1], -log10(ret[2] + 1e-4))
})
indel.proximal <- sapply(meta$EID, function(eid) {
    load(sprintf('sIndels_vs_%s_H3K27ac_peaks.SUMMARY.rda', eid))
    ret <- unlist(es[[1]]['TSS_proximal',c('enr', 'pval')])
    c(ret[1], -log10(ret[2] + 1e-4))
})


pf <- function(x, xticks, no.x=FALSE, no.y=FALSE, ...) {
    anatomy <- meta[colnames(x)]$ANATOMY
    type <- meta[colnames(x)]$TYPE
    pfc <- meta[STD_NAME=='Brain_Dorsolateral_Prefrontal_Cortex']$EID

    ylab <- if (no.y) '' else 'Significance: -log10(p-value)'
    xlab <- if (no.x) '' else 'Enrichment (or depletion) ratio (obs/exp)'
    plot(t(x),
        col=ifelse(anatomy=='BRAIN' & type=='PrimaryTissue',
                   'red', ifelse(x[2,] > 1, 'black', 'grey')),
        bty='n', xaxt='n', yaxt='n', pch=20,
        xlim=range(xticks), ylim=c(0,4), cex=2,
        ylab=ylab, xlab=xlab, ...)
    # replot brain points on top of others
    points(t(x[,anatomy=='BRAIN' & type=='PrimaryTissue']), col='red',
        pch=20, cex=2)
    abline(h=1, lty='dotted', col='black', lwd=0.55)
    abline(v=1, lty='dotted', col='black', lwd=0.55)
    segments(1, 2, x[1, pfc], x[2, pfc], col='red')
    if (!no.y) axis(side=2, at=0:4, lwd=0.35, lwd.ticks=0.35)
    if (!no.x) axis(side=1, at=xticks, lwd=0.35, lwd.ticks=0.35)
}

pdf(width=5.0, height=3.75, pointsize=8, file='panel_ef.pdf')
layout(matrix(1:4,nrow=2))
par(mar=c(2,5,3,1))
pf(snv.distal, xticks=seq(0.6,1.4,0.2), main='Somatic SNVs', no.x=T)
par(mar=c(4,5,1,1))
pf(snv.proximal, xticks=seq(0.6,1.4,0.2))
par(mar=c(2,5,3,1))
pf(indel.distal, xticks=seq(0.5,2.5,0.5), main='Somatic Indels', no.x=T, no.y=T)
par(mar=c(4,5,1,1))
pf(indel.proximal, xticks=seq(0.5,2.5,0.5), no.y=T)
dev.off()
