#!/usr/bin/env Rscript

segnorm <- read.table('../data/Ginkgo_100kb_varbin_SegNorm.txt.gz',
    sep='\t', header=T, stringsAsFactors=F, check.names=F)
ploidy <- read.table('../data/Ginkgo_100kb_varbin_combined_ploidies.txt',
    sep='\t', header=T, stringsAsFactors=F, check.names=F)

bulk <- '5657_0717.hrt.1b1'
bulk.ploidy <- ploidy[ploidy$sample == bulk,2]
# PTA samples are all very similar
pta <- '5657PFC.A.b37'
pta.ploidy <- ploidy[ploidy$sample == pta,2]
# MDA samples are similar except 1D2
mda <- '5657_0902.pfc.1cp1E11'
mda.ploidy <- ploidy[ploidy$sample == mda,2]


allchrs <- paste0('chr', c(1:22,c('X','Y')))
seps <- sapply(split(segnorm[,2], segnorm[,1])[allchrs], length)

plotf <- function(chrs, cns, ...) {
    par(mar=c(2,2,0.5,0.5))

    # use alternating point colors instead of background shades
    colmap <- setNames(rep(c('#888888', '#000000'), 12), allchrs)
    plot(cns, xaxt='n', yaxt='n', ylim=c(0,8), yaxs='i', xaxs='i', xpd=TRUE, bty='n',
        col=colmap[chrs], pch=20, cex=1/8, ylab='')

    axis(side=2, labels=seq(0,8,2), at=seq(0,8,2), lwd=0.35, lwd.ticks=0.35)
}

pdf(width=2.6, height=1.4, file='panel_b.pdf', pointsize=5)
layout(1:3)
plotf(segnorm[[1]], segnorm[[bulk]]*bulk.ploidy)
plotf(segnorm[[1]], segnorm[[mda]]*mda.ploidy)
plotf(segnorm[[1]], segnorm[[pta]]*pta.ploidy)
axis(side=1, labels=c(1:22, 'X', 'Y'), at=cumsum(seps) - seps/2, lwd=0, lwd.ticks=0)
dev.off()
