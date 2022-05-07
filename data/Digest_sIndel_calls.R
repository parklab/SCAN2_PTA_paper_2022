#!/usr/bin/env Rscript
library(scan2)

# Set all of the below to T to generate the table and save it
step1.load=F
step2.rescore=F                     # make 1% FDR calls and rescue
step3.master_table=F                # join all raw calls into a single dataframe
step4.exact_recurrence_filter=F     # SHOULDN'T have an effect on indels because of
                                    # the cross-sample filter.
step5.clustered_mut_filter=F        # removes a small number of indels
step6.finalize=F
step7.save=T

if (step1.load) {
    load('../metadata.rda')

    indel.rdas <- sprintf("../scan2_somatic_genotype_rdas/sIndels/%s.rda",
        meta$sample[meta$amp=='PTA' | (meta$amp=='MDA' & meta$type=='PFC neuron')])
    cat("PRE-FILTERING by FDR < 0.5 and cross-sample filter\n")
    ss <- lapply(indel.rdas, function(fn) {
        print(fn)
        load(fn)
        before <- nrow(somatic)
        somatic <- rescore(somatic, target.fdr=0.5, use.pon=TRUE, quiet=TRUE)
        somatic <- somatic[somatic$pass,]
        after <- nrow(somatic)
        cat('before:', before, 'after FDR < 0.5:', after, '\n')
        rescore(somatic, target.fdr=0.01, use.pon=TRUE, quiet=TRUE)
    })
    names(ss) <- meta$sample[meta$amp=='PTA' | (meta$amp=='MDA' & meta$type=='PFC neuron')]
}


if (step2.rescore) {
    ss01.mda <- lapply(ss[meta$sample[meta$amp=='MDA' & meta$type=='PFC neuron']],
        rescore, target.fdr=0.01, use.pon=T)

    # 1% FDR calls (and, implied by use.pon=T, apply a higher min. depth of 10)
    ss01.pta <- lapply(ss[meta$sample[meta$amp=='PTA']],
        rescore, target.fdr=0.01, use.pon=T)

    # remember: rescue2 is for indels
    cat('Applying indel rescue to PTA\n')
    ss.pta.resc <- do.rescue2(ldf=ss01.pta, calling.fdr=0.01, rescue.fdr=0.01)
}


# Make master table
if (step3.master_table) {
    # Convert list to table
    cols.for.tab <- c('chr', 'pos', 'refnt', 'altnt', 'dbsnp',
        'dp', 'af', 'bulk.dp', 'gp.mu', 'gp.sd', 'ab',
        paste0('pass', c('A','B','M')))
    ni.pta <- do.call(rbind, lapply(1:length(ss01.pta), function(i) {
        s <- ss01.pta[[i]]
        s$passA <- !is.na(s$pass) & s$pass
        # df[[i]]$pass should be equivalent to s$pass
        if (!all(s$pass == ss.pta.resc$df[[i]]$pass)) {
            cat(table(s$pass, ss.pta.resc$df[[i]]$pass))
            stop('sanity check: ss01 pass and ss.pta.resc pass are not identical but shuold be')
        }
        s$passB <- ss.pta.resc$df[[i]]$pass | ss.pta.resc$df[[i]]$pass2
        s$passM <- FALSE   # no MDA indels here
        sampname <- names(ss01.pta)[i]
        print(sampname)
        ret <- data.frame(sample=sampname,
            s[, cols.for.tab],
            donor=meta[sampname,]$donor,
            stringsAsFactors=F)
        # It is OK to restrict to passing mutations early for indels, but not SNVs.
        # Indels were already analyzed by the cross-sample filter, so duplicates
        # have already been removed when appropriate.
        # It will still be necessary to remove the allowed duplicates for some
        # downstream analyses (e.g., enrichment) later.
        ret <- ret[ret$passA | ret$passB | ret$passM,]
        ret$id <- paste(ret$chr, ret$pos, ret$refnt, ret$altnt)
        # This is the identical-length-but-position-off-by-1 alternative representation
        ret$id2 <- paste(ret$chr, ret$pos, nchar(ret$refnt) - nchar(ret$altnt))
        ret
    }))
    rownames(ni.pta) <- NULL

    ni.mda <- do.call(rbind, lapply(1:length(ss01.mda), function(i) {
        s <- ss01.mda[[i]]
        sampname <- names(ss01.mda)[i]
        if (nrow(s) == 0) {
            cat('sample', sampname, 'has 0 indels (MDA)\n')
            return(NULL)
        }
        s$passA <- FALSE
        s$passB <- FALSE
        s$passM <- !is.na(s$pass) & s$pass
        print(sampname)
        ret <- data.frame(sample=sampname,
            s[, cols.for.tab],
            donor=meta[sampname,]$donor,
            stringsAsFactors=F)
        # It is OK to restrict to passing mutations early for indels, but not SNVs.
        # Indels were already analyzed by the cross-sample filter, so duplicates
        # have already been removed when appropriate.
        # It will still be necessary to remove the allowed duplicates for some
        # downstream analyses (e.g., enrichment) later.
        ret <- ret[ret$passA | ret$passB | ret$passM,]
        ret$id <- paste(ret$chr, ret$pos, ret$refnt, ret$altnt)
        ret$id2 <- paste(ret$chr, ret$pos, nchar(ret$refnt) - nchar(ret$altnt))
        ret
    }))
    rownames(ni.mda) <- NULL
    
    ni <- classify.indels(rbind(ni.pta, ni.mda))
    rownames(ni) <- NULL
}



# Same recurrence filter used for SNVs, which does nothing for indels,
# because of the cross-sample filter. Still applying it for consistency, though.
if (step4.exact_recurrence_filter) {
    cat('Raw recurrence rates (PTA):\n')
    print(table(table(ni$id[ni$passA | ni$passB])))
    cat('Raw recurrence rates (MDA):\n')
    print(table(table(ni$id[ni$passM])))
    cat('Raw recurrence rates (PTA+MDA):\n')
    print(table(table(c(ni$id))))

    cat('Recurrence x donor table (all muts)\n')
    z <- split(ni$donor, ni$id)
    donors <- sapply(z, function(v) length(unique(v)))
    recs <- sapply(z, length)
    print(addmargins(table(recs,donors)))
    ni$rec.filter <- donors[ni$id] > 1

    cat('Recurrence x donor table (any passA,B,M)\n')
    z <- split(ni$donor[ni$passA | ni$passB | ni$passM],
               ni$id[ni$passA | ni$passB | ni$passM])
    donors <- sapply(z, function(v) length(unique(v)))
    recs <- sapply(z, length)
    print(addmargins(table(recs,donors)))
}


# Same clustering filter applied to SNVs
if (step5.clustered_mut_filter) {
    filter.single.sample.clusters <- function(muts, threshold=300) {
        by.sample <- split(muts, muts$sample)
        muts <- do.call(rbind, lapply(by.sample, function(d) {
            do.call(rbind, lapply(split(d, d$chr), function(dchr) {
                dchr <- dchr[order(dchr$pos),]
                up <- c(Inf, diff(dchr$pos))
                down <- c(diff(dchr$pos), Inf)
                dchr$nearest <- pmin(up, down)
                dchr
            }))
        }))
        muts <- muts[order(muts$chr, muts$pos),]
        filt.name <- paste0('clustered.filt.', threshold)
        muts[[filt.name]] <- muts$nearest < threshold
        cat(sprintf('Removing %d sites within %d bp in the same sample\n',
            sum(muts[[filt.name]]), threshold))
        print(addmargins(table(muts$sample, muts[[filt.name]])))
        rownames(muts) <- NULL
        muts
    }
    
    ni <- filter.single.sample.clusters(ni, threshold=50)
}


if (step6.finalize) {
    ni$lineage.filter <- duplicated(ni$id)
    ni$final.filter <- ni$rec.filter | ni$clustered.filt.50 | ni$lineage.filter
    ni$status <- ifelse(ni$passA, 'A', ifelse(ni$passB, 'B', ifelse(ni$passM, 'M', 'ERROR')))
}


if (step7.save) {
    outrda="Collected_SCAN2_sIndel_calls.rda"
    if (file.exists(outrda))
        stop(paste('output file', outrda, 'already exists, please delete it first'))
    save(ni, file=outrda)

    outcsv="Collected_SCAN2_sIndel_calls.csv"
    if (file.exists(outcsv))
        stop(paste('output file', outcsv, 'already exists, please delete it first'))
    write.csv(ni, file=outcsv)
}
