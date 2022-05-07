#!/usr/bin/env Rscript

tab <- read.csv('../data/Collected_chrX_simple_genotyping_indels.csv')

pdf(file='panel_b.pdf', width=1.5, height=2, pointsize=5)
par(mar=c(5, 4, 1, 1), bty='n')
burdens <- do.call(c, lapply(split(tab, tab$age), function(x)
    list(MDA=x$somatic.burden[x$amp=='MDA'], PTA=x$somatic.burden[x$amp=='PTA'])))
bp <- boxplot(burdens,
    col=c('#FA8D5F', '#93BFDA'), xaxt='n', yaxt='n',
    ylab='Chromosome X sIndel burden', ylim=c(0,20),
    pars=list(boxlwd=0.35, medlwd=1, whisklwd=0.35, staplelwd=0.35))
stripchart(burdens, vertical=T, add=T, pch=c(1, 17), lwd=0.35)
legend('topleft', fill=c('#FA8D5F', '#93BFDA'),
    legend=c('MDA','PTA'), bg='white')

xat <- sapply(split(1:10, rep(1:5, each=2)), mean)
axis(side=1, at=xat, labels=c(0.4, 0.6, 17, 45, 82), lwd=0, lwd.ticks=0, line=0)
axis(side=1, at=xat, labels=c(1278, 5817, 1465, 5087, 5657), lwd=0, lwd.ticks=0, las=3, line=1.5)
# For some reason 0.6 doesn't print in the above axis
axis(side=1, at=xat[2], labels=c(0.4, 0.6, 17, 45, 82)[2], lwd=0, lwd.ticks=0, line=0)

axis(side=2, at=c(0,5,10,15,20), lwd=0.35, lwd.ticks=0.35)
dev.off()

mda.vs.pta <- sapply(split(tab, tab$age),
    function(d) sapply(split(d$somatic.burden, d$amp), mean))

cat("PLOT NOTE:\n")
cat("PLOT NOTE:\n")
cat("PLOT NOTE: Manually add point at 60.645 for 5657 MDA\n")
cat("PLOT NOTE:\n")
cat("PLOT NOTE:\n")

cat("Average sensitivity adjusted SNV burden per male\n")
print(mda.vs.pta)

cat("Excess SNV burden per male\n")
print(mda.vs.pta[1,]-mda.vs.pta[2,])

cat("Average excess SNV burden\n")
print(mean(mda.vs.pta[1,]-mda.vs.pta[2,]))
cat("Median excess SNV burden\n")
print(median(mda.vs.pta[1,]-mda.vs.pta[2,]))

cat("Simple X------------\n")
cat("Simple X Infant SNVs (raw calls):\n")
print(sapply(split(tab$n.somatic.indels[tab$age < 1], tab$amp[tab$age < 1]), mean))
cat("Simple X infant SNVs (corrected burdens):\n")
print(sapply(split(tab$somatic.burden[tab$age < 1], tab$amp[tab$age < 1]), mean))
