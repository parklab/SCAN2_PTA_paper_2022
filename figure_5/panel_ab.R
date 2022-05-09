#!/usr/bin/env Rscript

load('Enrichment_results/sSNVs_vs_GTEx_all_tissues.SUMMARY.rda',verb=T)
# Reorder enrichment objects by quantiles 1..10
snvs <- lapply(es, function(e) {
    e <- e[rownames(e) %in% 1:10,]
    cbind(sort(as.integer(rownames(e))), e[order(as.integer(rownames(e))),]$enr)
})

pdf(width=2, height=2, pointsize=5, file='panel_a.pdf')
par(mar=c(5,4,1,1))
for (i in 1:54) {
    if (i==1) plotf <- plot else plotf <- lines
    plotf(snvs[[i]], type='l', xlim=c(1,10), ylim=c(0.8,1.2),
        xlab='GTEx gene expression\nDecile', ylab='Somatic SNV (obs/exp)',
        xaxt='n', yaxt='n', bty='n', lwd=2*0.35)
}
# Replot the brain lines on top to add emphasis
for (i in grep('Brain', names(es))) {
    lines(snvs[[i]], type='l', col='red', lwd=2*0.35,
        xaxt='n', yaxt='n', bty='n')
}
abline(h=1, lty='dashed', col='grey')
axis(side=1, at=1:10, lwd=0.35, lwd.ticks=0.35)
axis(side=2, at=8:12/10, lwd=0.35, lwd.ticks=0.35)
dev.off()




load('Enrichment_results/sIndels_vs_GTEx_all_tissues.SUMMARY.rda',verb=T)
# Reorder enrichment objects by quantiles 1..10
indels <- lapply(es, function(e) {
    e <- e[rownames(e) %in% 1:10,]
    cbind(sort(as.integer(rownames(e))), e[order(as.integer(rownames(e))),]$enr)
})

pdf(width=2, height=2, pointsize=5, file='panel_b.pdf')
par(mar=c(5,4,1,1))
for (i in 1:54) {
    if (i==1) plotf <- plot else plotf <- lines
    plotf(indels[[i]], type='l', xlim=c(1,10), ylim=c(0.5,2.2),
        xlab='GTEx gene expression\nDecile', ylab='Somatic indel (obs/exp)',
        xaxt='n', yaxt='n', bty='n', lwd=2*0.35)
}
# Replot the brain lines on top to add emphasis
for (i in grep('Brain', names(es))) {
    lines(indels[[i]], type='l', col='red', lwd=2*0.35,
        xaxt='n', yaxt='n', bty='n')
}
abline(h=1, lty='dashed', col='grey')
axis(side=1, at=1:10, lwd=0.35, lwd.ticks=0.35)
axis(side=2, at=1:4/2, lwd=0.35, lwd.ticks=0.35)
dev.off()
