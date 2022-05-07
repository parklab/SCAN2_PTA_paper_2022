load("../metadata.rda")

# Contains an 'abs' object of allele balances. List with one entry per
# cell/bulk. ABs are sorted and run-length encoded to save space.
load("../data/Collected_allele_balance_at_training_hSNPs.rda")

pdf(height=2, width=3.65, file='panel_d.pdf', pointsize=5)
layout(t(1:3))
par(mar=c(5,4,3,1))
for (amp in c('MDA', 'PTA', 'bulk')) {
    firstplot=TRUE
    for (i in 1:length(abs)) {
        a <- abs[[i]]
        n <- names(abs)[i]
        if (meta[n,]$amp == amp & n != '5087pfc-Rp3Cf') {
            f <- ifelse(firstplot,plot,lines)
            firstplot <- FALSE
            f(density(inverse.rle(a), from=0, to=1, na.rm=T),
                col=ifelse(amp=='PTA',"#91BFDB44",ifelse(amp=='MDA',"#FC8D5944",'#00000044')),
                type='l', main=amp, xlab='', ylab='', xaxt='n', yaxt='n', bty='n',
                ylim=c(0,6.2))
        }
    }
    axis(side=1, at=0:2/2, lwd=0.35, lwd.ticks=0.35)
    axis(side=2, at=0:6, lwd=0.35, lwd.ticks=0.35)
}
dev.off()
