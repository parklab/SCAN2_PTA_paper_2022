#!/usr/bin/env Rscript
library(data.table)

meta <- fread('../external_data/Roadmap_Epigenomics/EID_metadata_with_H3K27ac.txt')
setkey(meta, EID)

states <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk",
    "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk",
    "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")

snv.states <- lapply(states, function(s) {
    # Add 1e-4 to p-value since 0 out of 10,000 permutations is better
    # interpreted as p < 1e-4 rather than p=0.
    snv <- sapply(meta$EID, function(eid) {
        load(sprintf('sSNVs_vs_%s_ChromHMM15.SUMMARY.rda', eid))
        ret <- unlist(es[[1]][s,c('enr', 'pval')])
        c(ret[1], -log10(ret[2] + 1e-4))
    })
})

indel.states <- lapply(states, function(s) {
    indel <- sapply(meta$EID, function(eid) {
        load(sprintf('sIndels_vs_%s_ChromHMM15.SUMMARY.rda', eid))
        ret <- unlist(es[[1]][s,c('enr', 'pval')])
        c(ret[1], -log10(ret[2] + 1e-4))
    })
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
        ylim=c(0,4), cex=2,
        ylab=ylab, xlab=xlab, ...)
    # replot brain points on top of others
    points(t(x[,anatomy=='BRAIN' & type=='PrimaryTissue']), col='red',
        pch=20, cex=2)
    abline(h=1, lty='dotted', col='black', lwd=0.55)
    abline(v=1, lty='dotted', col='black', lwd=0.55)
    segments(1, 2, x[1, pfc], x[2, pfc], col='red')
    if (!no.y) axis(side=2, at=0:4, lwd=0.35, lwd.ticks=0.35)
    #if (!no.x) axis(side=1, at=xticks, lwd=0.35, lwd.ticks=0.35)
    if (!no.x) axis(side=1, lwd=0.35, lwd.ticks=0.35)
}

pdf(width=2*3.5, height=2*3.25, pointsize=8, file='ExtData_Fig8_sSNVs.pdf')
layout(matrix(1:16,nrow=4))
par(mar=c(4,4,2,1))
for (i in 1:length(states))
    pf(snv.states[[i]], main=states[i])
dev.off()



pdf(width=2*3.5, height=2*3.25, pointsize=8, file='ExtData_Fig8_sIndels.pdf')
layout(matrix(1:16,nrow=4))
par(mar=c(4,4,2,1))
for (i in 1:length(states))
    pf(indel.states[[i]], main=states[i])
dev.off()
