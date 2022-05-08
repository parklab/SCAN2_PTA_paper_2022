load("../metadata.rda")

par1=2699520
par2=154931044

# Keep male subjects only. UMB-prefix subjects and 5871 don't have
# matched MDA cells, so remove those too.
meta <- meta[!is.na(meta$sex) &
             meta$sex == 'M' &
             !grepl('UMB', meta$donor) &
             meta$donor != '5871',]
donors <- unique(meta$donor)
print(addmargins(table(meta$donor, meta$amp)))


cat("Finding germline indels\n")
indels <- do.call(rbind, lapply(donors, function(donor) {
        d <- paste0('sIndels/', donor, '.tab')
        cat('reading', d, '\n')
        g <- read.table(d, header=T, stringsAsFactors=F, check.names=F)

        # Get the approprate bulk
        bulk.sample <- meta$sample[meta$donor == donor & meta$amp == 'bulk']
        bulk.idx <- which(colnames(g) == bulk.sample)

        do.call(rbind, lapply(seq(8, ncol(g), 3), function(i) {
            x <- g[,c(1:7,i:(i+2),bulk.idx:(bulk.idx+2))]
            if (!(colnames(x)[8] %in% meta$sample) | meta[colnames(x)[8],]$amp == 'bulk')
                return(NULL)
            x$donor <- donor
            x$sample <- colnames(x)[8]
            x$sc.dp <- x[,10] + x[,9]
            x$sc.af <- x[,10] / x$sc.dp
            x$bulk.dp <- x[,12] + x[,13]
            x$bulk.alt <- x[,13]
            x$bulk.ref <- x[,12]
            x$bulk.af <- x$bulk.alt / x$bulk.dp
            x$is.called <- !is.na(x$sc.af) & x$sc.af > 0.9 &
                x$sc.dp > median(x$sc.dp) & x$bulk.dp > 10
            x$is.germline <- !is.na(x$bulk.af) & x$bulk.af > 0.9 &
                x$bulk.dp >= median(x$bulk.dp) & x$bulk.ref <= 2
            # columns 6-13 are redundant now, get rid of them
            x[x$pos >= par1 & x$pos <= par2,-(6:13)]
        }))
}))



cat("Finding somatic indels\n")
sindels <- do.call(rbind, lapply(donors, function(donor) {
        d <- paste0('/n/data1/hms/dbmi/park/jluquette/pta/',
                    donor, '/scansnv_justX/')
        print(d)
        g <- read.table(paste0(d, 'indel/mmq60.tab'),
                        header=T, stringsAsFactors=F, check.names=F)
        p <- read.table(paste0(d, 'indel/somatic_positions.chrX.tab'),
                        header=T, stringsAsFactors=F, check.names=F)
        # just subsets g by sites in p
        s <- merge(g, p)

        # Get the approprate bulk
        bulk.sample <- meta$sample[meta$donor == donor & meta$amp == 'bulk']
        bulk.idx <- which(colnames(s) == bulk.sample)

        do.call(rbind, lapply(seq(8, ncol(s), 3), function(i) {
            x <- s[,c(1:7,i:(i+2),bulk.idx:(bulk.idx+2))]
            if (!(colnames(x)[8] %in% meta$sample) | meta[colnames(x)[8],]$amp == 'bulk')
                return(NULL)
            x$donor <- donor
            x$sample <- colnames(x)[8]
            x$sc.dp <- x[,10] + x[,9]
            x$sc.af <- x[,10] / x$sc.dp
            x$bulk.dp <- x[,12] + x[,13]
            x$bulk.alt <- x[,13]
            x$is.somatic <- !is.na(x$sc.af) & x$sc.af > 0.9 & x$sc.dp > median(x$sc.dp) &
                x$bulk.dp > 10 & x$bulk.alt==0 & x$dbsnp == '.'
            x[x$pos >= par1 & x$pos <= par2,-(6:13)]
        }))
}))

write.csv(sindels, file='PROTECTED_sindels_chrX_simple_genotyping.csv') 
write.csv(indels, file='PROTECTED_germline_indels_chrX_simple_genotyping.csv')

tab <- do.call(rbind, lapply(split(indels, indels$sample), function(x) {
    N=sum(x$is.germline)
    y=sum(x$is.called[x$is.germline])
    data.frame(donor=x$donor[1], sample=x$sample[1],
        amp=meta[x$sample[1],]$amp, age=meta[x$sample[1],]$age,
        germline.sites=N, sites.called=y, sensitivity=y/N)
}))
muts <- do.call(rbind, lapply(split(sindels, sindels$sample), function(x) {
    data.frame(donor=x$donor[1], sample=x$sample[1], n.somatic.indels=sum(x$is.somatic))
}))
tab <- merge(tab, muts)
tab$somatic.burden <- tab$n.somatic.indels/tab$sensitivity

# Excluded MDA neuron
tab <- tab[tab$sample != '5087pfc-Rp3C5',]

write.csv(tab, file='Collected_chrX_simple_genotyping_indels.csv')
