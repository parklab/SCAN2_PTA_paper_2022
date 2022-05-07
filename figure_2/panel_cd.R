#!/usr/bin/env Rscript
library(data.table)
library(scan2)

load('../metadata.rda')

sigb <- read.csv('../data/Lodato2018_SignatureData_Aging.csv')
sigb <- sigb$B
scf <- read.csv('../data/Petljak2019_signatures.csv')
scf <- scf$SBS.sc_F

somatic <- read.csv('../data/Collected_SCANSNV_SCAN2_sSNV_calls.csv.gz')

# There's a third infant, but it does not have MDA neurons for comparison
somatic <- somatic[somatic$donor %in% c(1278, 5817),]


# df.to.sbs96 is in library(scan2)
mda.samples <- meta[meta$amp=='MDA',]$sample
infant.mda.sig <- df.to.sbs96(somatic[somatic$sample %in% mda.samples,])
pta.samples <- meta[meta$amp=='PTA',]$sample
infant.pta.sig <- df.to.sbs96(somatic[somatic$sample %in% pta.samples,])

pdf(width=0.85, height=1.9, pointsize=5, file='panel_c.pdf')
par(mar=c(6,4,1,0.1))
barplot(c(`Infant MDA`=sum(infant.mda.sig[33:48]),
    `Infant PTA`=sum(infant.pta.sig[33:48]),
    `Sig. B`=sum(sigb[33:48]),
    `Sig. scF`=sum(scf[33:48])),
    ylim=0:1, las=3, col="#434343", yaxt='n', space=0.5, border=NA, ylab='Fraction C>T')
axis(side=2, at=0:5/5, lwd=0.35, lwd.ticks=0.35)
dev.off()

plotf <- function(v, ...) {
    p <- barplot(v, xaxt='n',col=mutsig.cols, border=NA, space=0.5, xaxt='n', yaxt='n', ...)
    abline(v = (p[seq(4, length(p) - 1, 4)] + p[seq(5, length(p), 4)])/2, col = "grey", lwd=0.35)
}

# Also from library(scansnv)
mutsig.cols[c(35,39,43,47)] <- '#FD8D7E'


pdf(width=2, height=1.8, pointsize=5, file='panel_d.pdf')
layout(1:4)
par(mar=c(1,4,1,1))
plotf(infant.mda.sig, main="MDA")
plotf(infant.pta.sig, main="PTA")
plotf(sigb, main='Signature B (Lodato et al, 2018)')
plotf(scf, main='Signature scF (Petljak et al, 2019)')
dev.off()
