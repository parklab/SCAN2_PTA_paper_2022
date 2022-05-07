#!/usr/bin/env Rscript

library(data.table)
library(lme4)
library(lmerTest)

muttab <- fread('Collected_SCAN2_sSNV_burdens.csv')

donors <- unique(muttab$donor)
colors <- sapply(donors, function(dn) muttab[donor==dn]$color[1])
ages <- sapply(donors, function(dn) muttab[donor==dn,]$age[1])

# The yellow color in the table is too bright.
colors['5219'] <- "#B1AF27"  # for legend panel
muttab[muttab$donor==5219,]$color <- "#B1AF27"

pdf(width=2.40, height=1.8, pointsize=5, file='panel_a.pdf')
layout(t(1:2), widths=c(1,4))

# First plot is just the legend
par(mar=c(0,0,0,0))
plot(1, pch=NA, bty='n', xaxt='n', yaxt='n')
legend('center', legend=sub('UMB','',donors[order(ages)]),
    fill=colors[order(ages)], cex=5/6, bty='n', box.lwd=0.35)

par(mar=c(4,4,1,1))
plot(muttab$plotage, muttab$burden,
    ylim=c(0,6000),
    pch=ifelse(muttab$amp=='PTA', 24, 1),
    col=ifelse(muttab$amp=='PTA', 'black', muttab$color),  # color is the point outline
    lwd=ifelse(muttab$amp=='PTA', 1, 2)/2,
    bg=muttab$color, xaxt='n', yaxt='n', bty='n',
    xlab='Age', ylab='Total burden', main='', cex=1.5)
axis(side=1, at=c(0,20,40,60,80,100), lwd=0.35, lwd.ticks=0.35)
axis(side=2, at=1000*0:6, labels=0:6, lwd=0.35, lwd.ticks=0.35)


# Simple linear models for each
# Mixed-effects model: allow a donor-specific offset to remove the
# independence assumption between measurements.
plm <- lmer(burden ~ age + (1 | donor),
    data=muttab[amp=='PTA',])
mlm <- lmer(burden ~ age + (1 | donor),
    data=muttab[amp=='MDA',])
mlm2 <- lmer(burden ~ age + (1 | donor),
    data=muttab[amp=='MDA' & donor != 5219,])
abline(coef=fixef(plm), lwd=2*0.35)
abline(coef=fixef(mlm), lty='dashed', col='black', lwd=2*0.35)
cat('------------------ PTA --------------------\n')
print(summary(plm))
cat('------------------ MDA --------------------\n')
print(summary(mlm))
cat('------------------ MDA (no 5219) ----------\n')
print(summary(mlm2))

dev.off()
