#!/usr/bin/env Rscript
library(data.table)
library(lme4)
library(lmerTest)

load('../metadata.rda')
ni <- fread('../data/Collected_SCAN2_sIndel_calls.csv.gz')

samples <- meta[meta$amp=='PTA',]$sample

muttab <- data.frame(meta[samples,],
    nsom=sapply(samples, function(sn) sum(ni$sample == sn & ni$status == 'A')),
    ndel=sapply(samples, function(sn) sum(ni$sample == sn & ni$status == 'A' & grepl(':Del:', ni$muttype))),
    nins=sapply(samples, function(sn) sum(ni$sample == sn & ni$status == 'A' & grepl(':Ins:', ni$muttype))),
    stringsAsFactors=F)
    


pdf(width=1.8, height=1.8, pointsize=5, file='panel_c.pdf')
par(mar=c(4,4,1,1))
plot(muttab[,c('plotage', 'ndel')],
    pch=20, col=2, bty='n', cex=1/2,
    xaxt='n', yaxt='n',
    xlab='Age', ylab='PTA somatic indels')
axis(side=1, at=c(0,20,40,60,80,100), lwd=0.35, lwd.ticks=0.35)
axis(side=2, at=seq(0,60,10), lwd=0.35, lwd.ticks=0.35)
points(muttab[,c('plotage', 'nins')],
    pch=20, col=1, cex=1/2)
inslm <- lmer(nins ~ age + (1 | donor),
    data=muttab[muttab$amp=='PTA',])
dellm <- lmer(ndel ~ age + (1 | donor),
    data=muttab[muttab$amp=='PTA',])
abline(coef=fixef(dellm), lwd=2*0.35, col=2)
abline(coef=fixef(inslm), lwd=2*0.35)
cat('------------------ PTA | insertions --------------------\n')
print(summary(inslm))
cat('------------------ PTA | insertions --------------------\n')
print(summary(dellm))
legend('topleft', col=2:1, legend=c('Deletions', 'Insertions'), pch=16, bty='n')
dev.off()




pdf(width=1.8, height=1.6, pointsize=5, file='panel_d.pdf')
par(mar=c(4,4,1,1))
plot(table((nchar(ni$altnt[ni$passA])-nchar(ni$refnt[ni$passA]))),
    lwd=1.5, lend=1,
    bty='n', xaxt='n', yaxt='n', main='',
    xlab='Indel size (deletions < 0)', ylab='PTA somatic indels')
axis(side=2, at=seq(0,500,100), lwd=0.35, lwd.ticks=0.35)
axis(side=1, at=seq(-30,17,5), lwd=0.35, lwd.ticks=0.35)
dev.off()
