#!/usr/bin/env Rscript
library(data.table)
library(pracma)
library(scan2)

load('../metadata.rda')
ni <- fread('../data/Collected_SCAN2_sIndel_calls.csv.gz')

samples <- meta[meta$amp=='PTA',]$sample

# passA also implies amp=PTA
# don't keep duplicates for signature analysis
ni <- ni[ni$passA & !ni$final.filter,]


pdf(width=3.25, height=1.48, file='panel_e.pdf', pointsize=5)
plot.indel(ni$muttype, xaxt='n', reduce.to.id83=F, yaxt='n', border=NA)
axis(side=2, at=seq(0,80,20), lwd=0.35, lwd.ticks=0.35)
dev.off()


tmp <- fread('ID83_Indel_correction_factors.csv')
sens.by.id83 <- setNames(tmp$factor, tmp$muttype)
cosmic <- fread('COSMIC_ID_signatures.csv')


# Signature analysis of somatic indels combined across samples/cells
n <- plot.indel(ni$muttype, reduce.to.id83=F, make.plot=F)
nn <- n/sens.by.id83[names(n)]
expo <- lsqnonneg(as.matrix(cosmic[,-1]), nn)$x
names(expo) <- colnames(cosmic)[-1]
expo <- expo/sum(expo)

expo.uncorrected <- lsqnonneg(as.matrix(cosmic[,-1]), n)$x
names(expo.uncorrected) <- colnames(cosmic)[-1]
expo.uncorrected <- expo.uncorrected/sum(expo.uncorrected)

pdf(width=2.0, height=1.8, file='panel_f.pdf', pointsize=5)
par(mar=c(4,4,1,1))
expo.for.plot <- expo[paste0('ID',c(1,2,5,8,3,4,6,7,9:17))]
barplot(expo.for.plot, las=3, ylab='Signature exposure (ID83 corrected)',
    col=c('#4c4c4c', '#cecece')[2 - (expo.for.plot>0.05)], border=NA, yaxt='n')
axis(side=2, at=seq(0,0.2,0.05), lwd=0.35, lwd.ticks=0.35)
abline(h=0.05, lty='dashed')
legend('topright', legend='Not detected', fill='#cecece', bty='n')
dev.off()


# Now doing COSMIC fits on each cell individually
ns <- lapply(split(ni$muttype, ni$sample),
    function(mt) {
        # plot.indel returns the spectrum
        ret <- plot.indel(mt, reduce.to.id83=F, make.plot=F)
        ret / sens.by.id83[names(ret)]
    }
)
nmat <- do.call(cbind, lapply(ns, function(n) {
    expo <- lsqnonneg(as.matrix(cosmic[,-1]), n)$x
    names(expo) <- names(cosmic)[-1]
    expo
}))

# The yellow color in the table is too bright.
meta[meta$donor==5219,]$color <- "#B1AF27"

plot.vs.age <- function(expos, signame, ...) {
    plot(meta[colnames(expos),]$age, expos, pch=24,
        bg=meta[colnames(expos),]$color,
        lwd=1/2, cex=1.5, ...)
    m <- lm(nsom ~ age,
        data=data.frame(age=meta[colnames(expos),]$age, nsom=as.vector(expos)))
    options(digits=3)
    abline(coef=coef(m), lwd=2*0.35)
    rho <- round(cor(meta[colnames(expos),]$age, as.vector(expos)),3)
    pval <- coef(summary(m))[2,4]
    legend('topleft', bty='n', title=signame, title.adj=0, inset=0.05, x.intersp=0,
        legend=bquote(rho == .(rho)~' '~ P == .(pval)))
}

pdf(width=2, height=1.8, file='panel_g.pdf', pointsize=5)
par(mar=c(4,5,1,1))
plot.vs.age(nmat["ID4",,drop=F], signame="ID4",
    xaxt='n', yaxt='n', bty='n',
    xlab='Age', ylab='Somatic indel burden\n(ID83 corrected)')
axis(side=1, at=c(0,20,40,60,80,100), lwd=0.35, lwd.ticks=0.35)
axis(side=2, at=seq(0,120,20), lwd=0.35, lwd.ticks=0.35)
dev.off()
