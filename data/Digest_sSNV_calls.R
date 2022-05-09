#!/usr/bin/env Rscript

library(scan2)
library(data.table)

# Set all values below to TRUE to run the full script
step1.load=F                         # just load data from disk. takes a while.
step2.rescore=F                      # rescore sites to FDR=1% and do SCAN2 mutation signature rescue
step3.make_initial_tables=F          # use info from (1) and (2) to make passAB determinations
step4.exact_recurrence_filter=T      # filter out SNVs called in more than one donor.
step5.clustered_mut_filter=F         # remove SNVs within 50bp in the same sample
                                     # these are mostly dinuc subs.
step6.finalize=F    
step7.save_output=F                  # write (but don't overwrite) the rda file


if (step1.load) {
    load("../metadata.rda")
    meta <- meta[meta$type == 'PFC neuron' &
                (meta$amp == 'PTA' | meta$amp=='MDA') &
                !(meta$sample %in% c('5087pfc-Rp3C5')),]

    cat('Master sample list retained in metadata:\n')
    print(addmargins(table(meta$donor, paste(meta$amp, meta$type))))

    all.ss <- lapply(meta$sample, function(sn) {
        print(sn)
        load(sprintf("../scan2_somatic_genotype_rdas/sSNVs/%s.rda", sn))
        somatic
    })
    names(all.ss) <- meta$sample
    ss.full <- all.ss
}



# For PTA, "passA" refers to VAF-based sSNV calling at target.fdr=1%. This
# is identical to SCAN-SNV calling.
# "passB" refers to SCAN2's mutation-signature rescued calls.

if (step2.rescore) {
    # First, remove all calls with FDR>50%.
    ss.pta <- lapply(ss.full[meta$sample[meta$amp=='PTA']], function(s) {
        s <- rescore(s, target.fdr=0.5, use.pon=F, quiet=F)
        s[s$pass,]
    })

    # Below, we assume ss.pta and ss.pta.resc$df data.frames contain
    # the same set of somatic SNVs. However, do.rescue trims the data.frames
    # it returns because get.filter.reasons() is extremely slow. So why
    # should the pairs of data.frames be the same? The reason is that
    # ss.pta above keeps SNVs that pass at target.fdr=0.5, which means
    # that they must meet the requisite min.alt read count and have a non-NA
    # VAF (because min.alt>0). These are exactly the criteria that do.rescue
    # filters on.

    # Rescore to FDR=1%
    ss.pta <- lapply(ss.pta, rescore, target.fdr=0.01, use.pon=F, quiet=T)
    # SCAN2 signature-based rescue
    ss.pta.resc <- do.rescue(ss.pta[meta$sample[meta$amp=='PTA' & meta$type=='PFC neuron']],
        calling.fdr=0.01, rescue.fdr=0.01)
}



if (step3.make_initial_tables) {
    # Determine confidence class for every SNV
    ss2 <- lapply(names(ss.pta), function(sn) {
        s <- ss.pta[[sn]]
        s$passA <- s$pass
        s$passB <- ss.pta.resc$df[[sn]]$pass | ss.pta.resc$df[[sn]]$pass2
        s$passM <- FALSE  # For MDA calls
        s
    })
    names(ss2) <- names(ss.pta)


    # For MDA, retain only FDR=1% calls
    # Signature-based rescue is not applied to MDA.
    ss.mda <- lapply(names(ss.full[meta$sample[meta$amp=='MDA']]), function(sn) {
        s <- ss.full[[sn]]
        s <- s[s$pass,]
        # MDA mutations don't fall under any of these categories
        s$passA <- FALSE
        s$passB <- FALSE
        s$passM <- TRUE
        s
    })
    names(ss.mda) <- names(ss.full[meta$sample[meta$amp=='MDA']])

    # Combine lists of mutations into a single table
    make.mut.df <- function(slist) {
        muts <- do.call(rbind, lapply(1:length(slist), function(i) {
            s <- slist[[i]]
            cbind(sample=names(slist)[i],
                s[,c(1:5,14:31,
                    which(colnames(s) %in% paste0('pass', c('A','B','M'))))])
        }))
        muts <- get.3mer(muts)
        muts$sample <- as.character(muts$sample)
        muts$donor <- meta[muts$sample,]$donor
        muts$chr <- as.integer(muts$chr)
        muts <- muts[order(muts$chr, muts$pos),]
        muts$id <- paste(muts$chr, muts$pos, muts$refnt, muts$altnt)
        muts
    }

    nmut <- make.mut.df(c(
        ss2[meta$sample[meta$amp=='PTA' & meta$type=='PFC neuron']],
        ss.mda[meta$sample[meta$amp=='MDA' & meta$type=='PFC neuron']]
    ))

    cat('Initial statistics:\n')
    cat('PTA:\n')
    print(sum(nmut$passA))
    print(sum(nmut$passB))
    cat('MDA:\n')
    print(sum(nmut$passM))
}



# Remove SNVs called in 2 different individuals.
# This does NOT filter SNVs called more than once in the same
# individual (=likely lineage marker). That filter is applied later.
if (step4.exact_recurrence_filter) {
    cat('Raw recurrence rates (PTA):\n')
    print(table(table(nmut$id[nmut$passA | nmut$passB])))
    cat('Raw recurrence rates (MDA):\n')
    print(table(table(nmut$id[nmut$passM])))
    cat('Raw recurrence rates (PTA+MDA):\n')
    print(table(table(c(nmut$id))))

    cat('Recurrence x donor table (all muts)\n')
    z <- split(nmut$donor, nmut$id)
    donors <- sapply(z, function(v) length(unique(v)))
    recs <- sapply(z, length)
    print(addmargins(table(recs,donors)))
    nmut$rec.filter <- donors[nmut$id] > 1

    # This is the correct" way to remove mutations called in 2
    # or more donors. However, the above filter (which does not
    # require the calls to be passed, only to be in the list of
    # FDR < 50% candidate sites) finds an additional 57 sSNVs that
    # recur exactly in another brain with fairly high read support.
    # I think it's best to regard these as likely FPs and remove them.
    cat('Recurrence x donor table (any passA,B,M)\n')
    z <- split(nmut$donor[nmut$passA | nmut$passB | nmut$passM],
               nmut$id[nmut$passA | nmut$passB | nmut$passM])
    donors <- sapply(z, function(v) length(unique(v)))
    recs <- sapply(z, length)
    print(addmargins(table(recs,donors)))
    #nmut$rec.filter <- donors[nmut$id] > 1
}


# nearby points created by a single sample are more likely to
# be artifacts. Remove the whole cluster, because it is often
# true that the entire cluster is caused by the same few reads
# that probably align poorly or are clipped.
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
        muts
    }
    
    nmut <- filter.single.sample.clusters(nmut, threshold=50)
}


if (step6.finalize) {
    nmut <- nmut[nmut$passA | nmut$passB | nmut$passM,]

    # if there are any recurrences left, be sure to only pick one. if the
    # mutation is observed in both MDA and PTA, prefer the PTA one. otherwise,
    # keep one of any copy.
    # these recurrent mutations that occur in the same subject are likely
    # lineage-related true mutations. it's important not to count them as
    # multiple independent mutations or widely-inherited lineage markers, e.g.,
    # could drive enrichment signals.
    ord.samples <- c(names(ss.pta), names(ss.mda))
    nmut <- do.call(rbind, lapply(ord.samples, function(sn)
        nmut[nmut$sample == sn,]))
    nmut$lineage.filter <- duplicated(nmut$id)
    #nmut <- nmut[!nmut$duplicated,]
    nmut <- nmut[order(nmut$chr, nmut$pos),]
    nmut$final.filter <- nmut$rec.filter | nmut$clustered.filt.50 | nmut$lineage.filter

    # give a pass status that assigns each mutation to its most confident set
    # the priority of pass status is A > B > M
    # so if a mutation is both passA and passB, it gets assigned to the passA set.
    nmut$status <- ifelse(nmut$passA, 'A',
        ifelse(nmut$passB, 'B',
                    ifelse(nmut$passM, 'M', 'ERROR')))
}

if (step7.save_output) {
    cat('step7\n')
    rdafile='Collected_SCANSNV_SCAN2_sSNV_calls.rda'
    if (file.exists(rdafile))
        stop(paste('output RDA file', rdafile, 'already exists, please delete it first'))
    save(nmut, file=rdafile)

    csvfile='Collected_SCANSNV_SCAN2_sSNV_calls.csv'
    if (file.exists(csvfile))
        stop(paste('output csv file', csvfile, 'already exists, please delete it first'))
    fwrite(nmut, file=csvfile)
}
