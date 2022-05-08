#!/usr/bin/env Rscript

enhprom <- lapply(c('enhancers', 'promoters'), function(regtype) {
    setNames(lapply(c('sSNVs', 'sIndels'), function(muttype) {
        sapply(c('astrocyte', 'neuron', 'microglia', 'oligo'), function(celltype) {
            load(sprintf("%s_vs_%s_%s.SUMMARY.rda", muttype, celltype, regtype))
            es[[1]]['inside',c('enr', 'enr.boot.0.95.lb', 'enr.boot.0.95.ub')]
        })
    }), c('sSNVs', 'sIndels'))
})


atac <- setNames(lapply(c('sSNVs', 'sIndels'), function(muttype) {
    sapply(c('GABA', 'GLU', 'OLIG', 'MGAS'), function(celltype) {
        load(sprintf("%s_vs_%s_ATACseq.SUMMARY.rda", muttype, celltype))
        es[[1]]['inside',c('enr', 'enr.boot.0.95.lb', 'enr.boot.0.95.ub')]
    })
}), c('sSNVs', 'sIndels'))

dnarepair <- setNames(lapply(c('sSNVs', 'sIndels'), function(muttype) {
    sapply(c('SARseq', 'Repairseq'), function(assaytype) {
        load(sprintf("%s_vs_%s.SUMMARY.rda", muttype, assaytype))
        es[[1]]['inside',c('enr', 'enr.boot.0.95.lb', 'enr.boot.0.95.ub')]
    })
}), c('sSNVs', 'sIndels'))


blue <- "#1DAAFC"
point.ests <- unlist(c(enhprom[[1]]$sSNVs['enr',], enhprom[[1]]$sIndels['enr',],
       enhprom[[2]]$sSNVs['enr',], enhprom[[2]]$sIndels['enr',],
       atac$sSNVs['enr',], atac$sIndels['enr',],
       dnarepair$sSNVs['enr',], dnarepair$sIndels['enr',]))
lbs <- unlist(c(enhprom[[1]]$sSNVs['enr.boot.0.95.lb',], enhprom[[1]]$sIndels['enr.boot.0.95.lb',],
       enhprom[[2]]$sSNVs['enr.boot.0.95.lb',], enhprom[[2]]$sIndels['enr.boot.0.95.lb',],
       atac$sSNVs['enr.boot.0.95.lb',], atac$sIndels['enr.boot.0.95.lb',],
       dnarepair$sSNVs['enr.boot.0.95.lb',], dnarepair$sIndels['enr.boot.0.95.lb',]))
ubs <- unlist(c(enhprom[[1]]$sSNVs['enr.boot.0.95.ub',], enhprom[[1]]$sIndels['enr.boot.0.95.ub',],
       enhprom[[2]]$sSNVs['enr.boot.0.95.ub',], enhprom[[2]]$sIndels['enr.boot.0.95.ub',],
       atac$sSNVs['enr.boot.0.95.ub',], atac$sIndels['enr.boot.0.95.ub',],
       dnarepair$sSNVs['enr.boot.0.95.ub',], dnarepair$sIndels['enr.boot.0.95.ub',]))

colors=rep(c('black', blue, 'black', blue, 'black', blue, 'black', blue),
        times=c(4,4,4,4,4,4,2,2))

pdf(width=4.3, height=2, pointsize=5, file='panels_ghij.pdf')
plot(point.ests,
    ylim=c(0.5, 3.5), bty='n', yaxt='n', pch=20, xaxt='n',
    ylab='Observed/expected', xlab='', col=colors)
arrows(1:length(point.ests), lbs, 1:length(point.ests), ubs,
    code=3, angle=90, col=colors, lwd=2*0.35, length=0)
axis(side=1, at=1:length(point.ests), labels=names(point.ests), las=3, lwd=0.35, lwd.ticks=0.35)
abline(v=c(4.5,12.5,20.5,26.5), col='grey', lwd=0.35)
abline(h=1, col='black', lwd=0.45, lty='dashed')
axis(side=2, at=seq(0.5,3.5,0.5), lwd=0.35, lwd.ticks=0.35)
dev.off()
