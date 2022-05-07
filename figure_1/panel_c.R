load("../metadata.rda")

meta <- meta[grep('5087-1cp|-dg', meta$sample, invert=T, value=T),]
meta <- meta[order(meta$age),]

mapds <- do.call(c, lapply(unique(meta$donor), function(sn)
    list(meta$mapd[meta$donor==sn & meta$amp=='MDA'],
        meta$mapd[meta$donor==sn & meta$amp=='PTA'],
        meta$mapd[meta$donor==sn & meta$amp=='bulk'])))
names(mapds) <- as.vector(t(outer(unique(meta$donor), c('MDA','PTA','Bulk'), paste)))

pdf(width=1.75, height=2, file='panel_c.pdf', pointsize=5)
cols <- c(PTA='#91BFDB', MDA='#FC8D59', Bulk='black')
par(mar=c(4.5,4,1,1), bty='n')
for (amp in c('MDA', "PTA", "Bulk")) {
    print((strsplit(' ', names(mapds)[grep(amp, names(mapds))])))
    boxplot(mapds[grep(amp, names(mapds))],
        col=cols[amp],
        ylim=c(0,1.8), main='',
        ylab='MAPD (amplification uniformity)',
        xlab='Subject ID',
        outline=F, add=amp != 'MDA', lwd=0.35, xaxt='n', yaxt='n')
    stripchart(mapds[grep(amp, names(mapds))], vert=T, pch=20,
        method='jitter', add=TRUE, cex=1/4, xaxt='n', yaxt='n')
}
axis(side=1, at=1:17, lwd=0.35, lwd.ticks=0.35, las=3, cex.axis=5/6,
    labels=sub('UMB', '', sapply(strsplit(names(mapds)[grep(amp, names(mapds))], ' '), head, 1)))
axis(side=2, at=0:3/2, lwd=0.35, lwd.ticks=0.35)
legend('topleft', legend=c('MDA','PTA'), fill=c('#FC8D59', '#91BFDB'), box.lwd=0.35)
dev.off()
