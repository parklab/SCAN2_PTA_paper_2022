df.to.id83 <- function(df) {
    if (!('muttype' %in% colnames(df)))
        df <- classify.indels(df)
    # reduce.to.id83: doing this in classify.indels now
    (plot.indel(df$muttype, reduce.to.id83=F, make.plot=F)+0.1)/nrow(df)
}

classify.indels <- function(df, sample.name='dummy', save.plot=F, auto.delete=T, chrs=c(1:22,'X')) {
    ret <- classify.muts(df=df, spectype='ID', sample.name=sample.name, save.plot=save.plot,
        auto.delete=auto.delete, chrs=chrs)
    ret$muttype <- substr(ret$muttype, 3, 11)
    ret
}


classify.muts <- function(df, spectype='SNV',
    sample.name='dummy', save.plot=F, auto.delete=T, chrs=1:22)
{
    if (nrow(df) == 0)
        return(df)

    recognized.spectypes <- c('SNV', 'ID')
    if (!(spectype %in% recognized.spectypes))
        stop(sprintf("unrecognized spectype '%s', currently only supporting %s",
            spectype, paste('"', recognized.spectypes, '"', collapse=', ')))

    require(SigProfilerMatrixGeneratorR)
    td <- paste0(tempdir())
    spmgd <- paste0(td, "/spmgr/")
    if (file.exists(spmgd) & auto.delete)
        unlink(spmgd, recursive=TRUE)
    dir.create(spmgd, recursive=TRUE)

    # convert this to use scansnv.df.to.vcf at some point
    # Write out the VCF
    ### From scansnv.to.vcf
    out.file <- paste0(spmgd, sample.name, '.vcf')
    f <- file(out.file, "w")
    vcf.header <- c("##fileformat=VCFv4.0", "##source=scansnv", 
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", 
            #sprintf("##contig=<ID=%s,length=%d>", fai[, 1], fai[,2]),
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", 
            "QUAL", "FILTER", "INFO", "FORMAT", sample.name), collapse = "\t"))
    writeLines(vcf.header, con = f)
    s <- df[!is.na(df$chr),]
    s <- do.call(rbind, lapply(chrs, function(chr) {
        ss <- s[s$chr == chr, ]
        ss[order(ss$pos), ]
    }))

    # SigProfilerMatrixGenerator doesn't classify duplicated mutations
    # from the same sample, it throws errors instead. It also will not
    # detect duplicates if they are not adjacent in the file.  If the
    # duplicate is not detected in this way, it causes the bug where
    # the final newdf dataframe does not match the input df.
    # To circumvent all these headaches: just remove duplicates up front.
    mutid <- paste(s$chr, s$pos, s$refnt, s$altnt)
    dupmut <- duplicated(mutid)
    cat("Removing", sum(dupmut), "/", nrow(s), "duplicated mutations before annotating\n")
    s <- s[!dupmut,]

    # will write, eg., position=7000000 as 7e6, which will
    # confuse sigprofilermatrixgenerator
    old.opt <- options('scipen')$scipen
    options(scipen=10000)
    writeLines(paste(s$chr, s$pos, s$dbsnp, s$refnt, s$altnt, 
        ".", "PASS", ".", "GT", "0/1", sep = "\t"), con = f)
    close(f)
    options(scipen=old.opt)

    mat <- SigProfilerMatrixGeneratorR(sample.name, 'GRCh37', spmgd, seqInfo=TRUE, plot=TRUE)

    # Read in the types
    annot.files <- paste0(spmgd, '/output/vcf_files/', spectype, '/', c(1:22,'X','Y'), "_seqinfo.txt")
    if (spectype == 'ID') {
        colclasses <- c(V2='character', V5='character', V6='character')
    } else if (spectype == 'SNV') {
        colclasses <- c(V2='character')
    }
    annots <- do.call(rbind, lapply(annot.files, function(f) {
        tryCatch(x <- read.table(f, header=F, stringsAsFactors=FALSE,
                colClasses=colclasses),
            error=function(e) NULL)
    }))
    if (spectype == 'ID') {
        colnames(annots) <- c('sample', 'chr', 'pos', 'iclass', 'refnt', 'altnt', 'unknown')
        newdf <- plyr::join(df, annots[2:6], by = colnames(annots)[-c(1, 4, 7)])
    } else if (spectype == 'SNV') {
        colnames(annots) <- c('sample', 'chr', 'pos', 'iclass', 'unknown')
        newdf <- plyr::join(df, annots[2:4], by = colnames(annots)[-c(1, 4, 5)])
    }

    if (save.plot) {
        plotfiles <- list.files(paste0(spmgd, '/output/plots/'), full.names=T)
        file.copy(plotfiles, '.')
    }

    if (!all(df$chr == newdf$chr))
        stop('df and newdf do not perfectly correspond: df$chr != newdf$chr')
    if (!all(df$pos == newdf$pos))
        stop('df and newdf do not perfectly correspond: df$pos != newdf$pos')

    if (auto.delete)
        unlink(spmgd, recursive=TRUE)
    df$muttype <- newdf$iclass
    df
}


get.sig.score2 <- function (ssnvs, good.sig, lysis.sig, eps = 0.001) 
{
    require(pracma)
    #test.sig <- df.to.sbs96(ssnvs)
    nsnvs <- sum(ssnvs$filter.reason == "lysis.test")
    if (nsnvs > 0) {
        test.sig <- df.to.id83(ssnvs)
        sigs <- cbind(good.sig, lysis.sig)
        weights <- lsqnonneg(sigs, test.sig)$x
        recon <- as.vector(sigs %*% weights)
    } else {
        # doesn't matter what we return here, downstream the weights will be ignored
        weights <- c(1,1)
    }
    weights <- weights + eps
    weights <- weights/sum(weights)
    postp <- log10(lysis.sig * weights[2]) - log10(good.sig * 
        weights[1])
    list(postp = postp, nsnvs = nsnvs, weight.true = weights[1], 
        weight.artifact = weights[2])
}

rescue2 <- function (df, lysis.sig, good.sig, rescue.fdr, ...) 
{
    sigscores <- get.sig.score2(ssnvs = df[df$filter.reason == 
        "lysis.test", ], lysis.sig = lysis.sig, good.sig = good.sig, 
        ...)
    postp <- sigscores$postp
    #df$rweight <- 10^-postp[df$type.and.ctx]
    df$rweight <- 10^-postp[df$muttype]
    df$lysis.fdr2 <- df$lysis.alpha/(df$lysis.alpha + df$lysis.beta * 
        df$rweight * df$nt/df$na)
# XXX: NOTE extra filters for indels: PON and dp >= 10 are HARD CODED!!
# update this in a later version.
    df$pass2 <- !df$pass & df$lysis.fdr2 <= rescue.fdr & df$filter.reasons == 
        "lysis.test" & (df$unique.donors <= 1 | df$max.out <= 2) & df$dp >= 10
    df$filter.reasons[df$pass2] <- paste0(df$filter.reasons[df$pass2], 
        ";rescue")
    list(df = df, postp = postp, nsnvs = sigscores$nsnvs, weight.true = sigscores$weight.true, 
        weight.artifact = sigscores$weight.artifact)
}

do.rescue2 <- function (ldf, lysis.sig=pta.indel.artifact.sig.v1, min.alt = 2, calling.fdr = 0.01, 
    rescue.fdr = 0.01) 
{
    if (missing(lysis.sig)) {
        data(pta.indel.artifact.sig.v1)
        lysis.sig <- pta.indel.artifact.sig.v1
    }
    cat(sprintf("SCAN2 multi-sample rescue on %d samples, %d high confidence SNVs from target.fdr=%0.4f\n", 
        length(ldf), sum(sapply(ldf, function(d) sum(d$pass))), 
        calling.fdr))
    cat("    Step 1: annotating filter reasons and SBS96 status..\n       ")
    ldf2 <- lapply(1:length(ldf), function(i) {
        cat(sprintf(" %s", names(ldf)[i]))
        df <- ldf[[i]]
        df <- df[!is.na(df$af) & round(df$af * df$dp) >= min.alt, 
            ]
        df$filter.reasons <- get.filter.reasons(df, min.alt = min.alt, 
            fdr.threshold = calling.fdr)
        #get.3mer(df)
        classify.indels(df)
    })
    names(ldf2) <- names(ldf)
    ldf <- ldf2
    cat(".\n")

    cat("    Step 2: constructing true SNV spectrum..\n")
    x <- do.call(rbind, lapply(ldf, function(s) s[s$pass, 1:5]))
    good.sig <- df.to.id83(x)

    cat(sprintf("    Step 3: adjusting FDR for lysis artifacts using new FDR target=%0.4f..\n        ", 
        rescue.fdr))
    final <- lapply(1:length(ldf), function(i) {
        df <- ldf[[i]]
        cat(sprintf(" %s", names(ldf)[i]))
        rescue2(df, lysis.sig = lysis.sig, good.sig = good.sig, 
            rescue.fdr = rescue.fdr)
    })
    names(final) <- names(ldf)
    cat(".\n")
    list(calling.fdr = calling.fdr, rescue.fdr = rescue.fdr, 
        good.sig = good.sig, lysis.sig = lysis.sig, df = lapply(final, 
            function(l) l$df), postp = lapply(final, function(l) l$postp), 
        nsnvs = sapply(final, function(l) l$nsnvs), weight.true = sapply(final, 
            function(l) l$weight.true), weight.artifact = sapply(final, 
            function(l) l$weight.artifact))
}
