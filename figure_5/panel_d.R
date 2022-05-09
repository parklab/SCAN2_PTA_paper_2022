#!/usr/bin/env Rscript

load('Enrichment_results/sSNVs_vs_Conservation.SUMMARY.rda',verb=T)
# Reorder enrichment objects by quantiles 1..10
snvs <- lapply(es, function(e) {
    e <- e[rownames(e) %in% 1:10,]
    cbind(sort(as.integer(rownames(e))), e[order(as.integer(rownames(e))),]$enr)
})

load('Enrichment_results/sIndels_vs_Conservation.SUMMARY.rda',verb=T)
# Reorder enrichment objects by quantiles 1..10
indels <- lapply(es, function(e) {
    e <- e[rownames(e) %in% 1:10,]
    cbind(sort(as.integer(rownames(e))), e[order(as.integer(rownames(e))),]$enr)
})


pdf(width=2, height=2, pointsize=5, file='panel_d.pdf')
par(mar=c(5,4,1,1))
plot(snvs[[1]], type='l', xlim=c(1,10), ylim=c(0.8,1.6),
    xlab='Conservation (phyloP 100)\nDecile', ylab='Observed/expected',
    xaxt='n', yaxt='n', bty='n', lwd=2*0.35)
# put points over lines (R's type='b' leaves an ugly space around the points)
points(snvs[[1]], pch=16)

lines(indels[[1]], lwd=2*0.35, col="#1DAAFC")
points(indels[[1]], pch=16, col="#1DAAFC")
abline(h=1, lty='dashed', col='grey')
axis(side=1, at=1:10, lwd=0.35, lwd.ticks=0.35)
axis(side=2, at=seq(0.8,1.6,0.2), lwd=0.35, lwd.ticks=0.35)
dev.off()
