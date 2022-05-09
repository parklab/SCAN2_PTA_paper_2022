#!/usr/bin/env Rscript

library(data.table)
library(lme4)
library(lmerTest)

muttab <- fread('Collected_SCAN2_sIndel_burdens.csv')

donors <- unique(muttab$donor)
colors <- sapply(donors, function(dn) muttab[donor==dn]$color[1])
ages <- sapply(donors, function(dn) muttab[donor==dn,]$age[1])

# The yellow color in the table is too bright.
muttab[muttab$donor=='5219',]$color <- "#B1AF27"


# Simple linear models for each
# Mixed-effects model: allow a donor-specific offset to remove the
# independence assumption between measurements.
plm <- lmer(burden ~ age + (1 | donor),
    data=muttab[amp=='PTA',])

# _MDA_ neurons from subjects 4638, 4643 and 5219 were not analyzed
# for indels. Do not include them in the trend line since they will
# have 0 calls.
mlm <- lmer(burden ~ age + (1 | donor),
    data=muttab[amp=='MDA' & !(donor %in% c(4638, 4643, 5219)),])
cat('------------------ PTA --------------------\n')
print(summary(plm))
cat('------------------ MDA --------------------\n')
print(summary(mlm))



muttab <- muttab[muttab$amp=='PTA',]  # Only PTA for indel panels b-d

pdf(width=1.9, height=1.8, pointsize=5, file='panel_b.pdf')
par(mar=c(4,4,1,1))
plot(muttab$plotage, muttab$burden,
    ylim=c(0,375),
    pch=ifelse(muttab$amp=='PTA', 24, 1),
    col=ifelse(muttab$amp=='PTA', 'black', muttab$color),  # color is the point outline
    lwd=ifelse(muttab$amp=='PTA', 1, 2)/2,
    bg=muttab$color, xaxt='n', yaxt='n', bty='n',
    xlab='Age', ylab='Total burden', main='', cex=1.5)
axis(side=1, at=c(0,20,40,60,80,100), lwd=0.35, lwd.ticks=0.35)
axis(side=2, at=c(0,50,100,150,200,250,300,350),
    labels=c(0,50,'',150,'',250,'',350), lwd=0.35, lwd.ticks=0.35)
# I think R decides we don't have enough space for all the labels.
# This overrides it.
for (x in c(50,150,250,350))
    axis(side=2, at=x, lwd=0.35, lwd.ticks=0.35)

abline(coef=fixef(plm), lwd=2*0.35)

dev.off()
