# General routines for calculating enrichment over a set of intervals
# defined by a GenomicRanges object.
# Usually these objects will come from BED files; gene models are the
# primary example of an exception.

# grl is a GRangesList, which allows findOverlaps to
# map all gmuts simultaneously. This is a 10x speedup,
# because it avoids the very significant overhead of
# multiple calls to findOverlaps.
# IMPORTANT: map.feature does NOT return a 1-1 mapping
# of mutations -> BED features!!! It is only useful for
# counting.
map.feature <- function(grl, edata, feat.name, verbose=TRUE) {
    thisg <- unlist(grl)
    ols <- GenomicRanges::findOverlaps(thisg, edata$gbed, type='any')

    # XXX: don't know if data.table is necessary anymore
    # data.table provides some important speed benefits for later grouping
    ret <- data.table::data.table(perm.id=thisg$perm.id[S4Vectors::from(ols)],
        to=features(edata=edata, feat.name=feat.name)[S4Vectors::to(ols)])

    ret
}

# grl is a GRangesList. This is critical for speed.
# map.and.count is a specialized implementation of edata$count.fn(edata$map.fn)
# that is MUCH more CPU and RAM efficient than calling the two functions
# separately.
# split grl into k slices and map each one individually. This
# provides a good tradeoff between performance and memory usage.
map.and.count <- function(grl, edata, verbose=1, target.k=10) {
    if (verbose > 1) {
        cat("map.and.count ----\n")
        cat('gc3\n')
        print(gc())
    }

    if (length(grl) >= target.k) {
        k = min(length(grl), target.k)
        z <- cut(1:length(grl), k, labels=FALSE)
    } else {
        k = 1
        # can't ask for 1 break in cut()
        z <- 1
    }

    if (verbose > 1) {
        cat('k =', k, '\n')
        print(table(z))
        cat('gc1a\n')
        print(gc())
    }

    # for loop: we split the GRL into k subgroups. Each subgroup contains
    # a roughly even number of GRL elements (e.g., a roughly even number
    # of permutations).
    ret <- c()
    for (i in 1:k) {
        if (verbose > 1) cat('findoverlaps', i, '-------------- \n')
        thisg <- unlist(grl[z == i])
        ols <- GenomicRanges::findOverlaps(thisg, edata$gbed, type='any')
        if (verbose > 1) { print(ols) }
        # aggregating 'mapping' like this creates a matrix with dimensions
        # (#permutations) x (#features), where #permutations is the number
        # of permutations in block 'k'.
        if (is.na(edata$use.mutclass)) {
            # new: data.table is just being used to efficiently split 'to' by 'from'
            mapping <- data.table::data.table(perm.id=thisg$perm.id[S4Vectors::from(ols)],
                to=S4Vectors::to(ols))
            if (verbose > 1) { print(mapping); print(gc()) }
            new.block <- mapping[,
                    do.call(c, lapply(feature.set(edata=edata), function(feat.name)
                        as.list(edata$count.fn(features(edata=edata, feat.name=feat.name)[to], feat.name=feat.name)))),
                by=perm.id][,-'perm.id']
        } else {
            # new: data.table is just being used to efficiently split 'to' by 'from'
            mapping <- data.table::data.table(perm.id=thisg$perm.id[S4Vectors::from(ols)],
                mutclass=mcols(thisg)[[edata$use.mutclass]][S4Vectors::from(ols)],
                to=S4Vectors::to(ols))
            if (verbose > 1) { print(mapping); print(gc()) }
            new.block <- mapping[,
                    do.call(c, lapply(feature.set(edata=edata), function(feat.name)
                        as.list(edata$count.fn(features(edata=edata, feat.name=feat.name)[to],
                                               feat.name=feat.name,
                                               mutclass=mutclass,
                                               signature.catalog=edata$signature.catalog)))),
                by=perm.id][,-'perm.id']
        }
        if (verbose > 1) print(new.block) 
        ret <- rbind(ret, new.block)
        rm(mapping)
        if (verbose > 0) cat('.')
    }
    if (verbose > 1) {
        cat('gc1b\n')
        print(gc())
    }

    ret
}

# Get counts for all possible feature types (including features
# with 0 count).
# x is a vector of the feature being measured (i.e., the return
# value of map.feature()).
# Using factors in edata$gbed allows tabulate to be used here,
# which ensures all factor levels are present and in the same order
# in the final count, even if some have 0 hits.
count.feature <- function(x, edata, feat.name) {
    setNames(tabulate(x, nbins=nlevels(x)), nm=paste0(feat.name, '|||', levels(x)))
}


# mutclass should be an ordered factor that matches the ordering
# used in the signature catalog.
# returns a vector of (feature name x feature value x signature exposure) values.
sigs.by.feature <- function(x, signature.catalog, feat.name, mutclass) {
    if (!is.factor(x)) {
        cat('x=')
        print(head(mutclass))
        stop('"x" must be a factor')
    }
    if (!is.factor(mutclass)) {
        cat('mutclass=')
        print(head(mutclass))
        stop('"mutclass" must be a factor')
    }

    # E.g., 96 x (# feature classes) table for SBS96 signatures
    # Because x and mutclass are factors, this table will include levels for
    # which no mutations are observed. This is critical, e.g., for including
    # 0s in the mutation spectra.
    tab <- table(mutclass, x)

    # Transformed to: (# signatures in catalog) x (# feature classes)
    # each value is the signature exposure estimated by least squares.
    sigs.per.feature <- apply(tab, 2, function(col) {
# Performance on a 6-signal run with a 10-signature catalog:
#> summaryRprof('Rprof.lsqnonneg.out')$sampling.time
#[1] 747.54
#> summaryRprof('Rprof.nnls.out')$sampling.time
#[1] 163.68
#> summaryRprof('Rprof.fnnls.out')$sampling.time
#[1] 567.2
#
# nnls and lsqnonneg produced essentially identical estimates:
#> summary(unlist(abs(el$perm.mat - en$perm.mat)))
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000e+00 0.000e+00 1.421e-14 2.333e-13 1.705e-13 2.092e-11 

        #z <- pracma::lsqnonneg(C=signature.catalog, d=col)
        #coefs <- z$x
        #resid <- z$resid.norm
        z <- nnls::nnls(A=signature.catalog, b=col)
        coefs <- coefficients(z)
        resid <- deviance(z)
        #z <- multiway::fnnls(XtX=crossprod(signature.catalog), Xty=crossprod(signature.catalog, col))
        #coefs <- as.vector(z)
        #resid <- sum((signature.catalog %*% coefs - col)^2)
        norm <- sum(col^2)
        pct.resid <- ifelse(norm > 0, resid/norm, 0)
#print(pct.resid)
        c(coefs, resid, norm, pct.resid)
    })

    # Unrolling a matrix in R is column-first
    setNames(as.vector(sigs.per.feature),
        nm=paste0(feat.name, '|||', rep(levels(x), each=nrow(sigs.per.feature)),
                  '|||', c(colnames(signature.catalog), c('resid.norm', 'norm', 'pct.resid'))))
}

feature.set <- function(en, edata) {
    if (!missing(en))
        edata <- en$edata

    attr(edata$gbed, 'feature.set')
}

# return an enrichment object with just 'feat.name' extracted
# strip.feat.name - remove the feature name prefix from each of the
#                   feature's levels.
get.feat <- function(en, feat.name, strip.feat.name=TRUE) {
    if (!(feat.name %in% feature.set(en)))
        stop(paste('tried to get feature', feat.name, 'from enrichment object, but no such feature is in feature.set'))

    feats <- names(en$real.obs)
    en <- subset(en, feats[sapply(strsplit(feats, '|||', fixed=TRUE), head, 1) == feat.name])
    en$feature.set <- feat.name
    if (strip.feat.name) {
        to.delete <- paste0(feat.name, '|||')
        colnames(en$perm.mat) <- sub(to.delete, '', colnames(en$perm.mat), fixed=TRUE)
        names(en$real.obs) <- sub(to.delete, '', names(en$real.obs), fixed=TRUE)
        if ('boots' %in% names(en)) {
            colnames(en$boots) <- sub(to.delete, '', colnames(en$boots), fixed=TRUE)
        }
    }
    en
}

features <- function(en, edata, feat.name) {
    if (!missing(en)) {
        # horrible hack, should use classes here
        if (is.data.frame(en))   # this is an esummary, not a full enrichment e
            return(rownames(en))
        edata <- en$edata
    }

    if (missing(feat.name)) {
        feats <- feature.set(edata=edata)
        if (length(feats) > 1)
            warn(paste('feat.name not specified, but feature set contains', length(feats), 'features. Using first feature.'))
        feat.name <- feats[1]
    }

    if (!(feat.name %in% attr(edata$gbed, 'feature.set')))
        stop(paste0('requested feature named "', feat.name, '", but no such feature in feature.set. Current feature set:', attr(edata$gbed, 'feature.set')))
    mcols(edata$gbed)[[feat.name]]
}

# Interface for enrichment analysis. At the minimum, requires
# counting and mapping functions. Additional data (such as a BED
# file of regions) can be specified in '...'.
# use.mutclass - seach mutation is annotated with a mutation class 
#     (called mutclass) and it should be passed to the count function
#     along with the mapped feature.
enrich.data <- function(count.fn=count.feature, map.fn=map.feature, use.mutclass=NA, ...) {
    list(count.fn=count.fn, map.fn=map.fn, use.mutclass=use.mutclass, ...)
}


read.bed.metadata <- function(bedfile, is.qbed=FALSE) {
    line <- readLines(bedfile, n=1)
    line <- sub('^#', '', line) # delete leading #comment character
    meta <- sapply(strsplit(line, ';')[[1]], function(x) {
        x <- strsplit(x, '=')[[1]]
        setNames(x[2], x[1])
    }, USE.NAMES=FALSE)
    if (is.qbed) {
        if (!("QBED_VERSION" %in% names(meta)))
            stop("QBED does not contain the required QBED_VERSION meta tag")
        if (meta['QBED_VERSION'] != "1")
            stop(paste("using unrecognized QBED file version", meta['QBED_VERSION']))
    }
    ret <- data.frame(t(meta[names(meta) != 'QBED_VERSION']))
    rownames(ret) <- bedfile
    ret
}


# Read a BED file with an optional 4th column indicating one or
# more interval classes. A 3 column BED will be interpreted as
# having a single state, with all intervals in the file corresponding
# to that state. Columns beyond 4 are ignored.
# The complement of the BED will be automatically determined and
# added as an additional state, which by default will be
# described by an empty string.
#
# is.qbed - is a quantile-scored BED. the important features of a
#           QBED are:
#           1. they tile the full genome
#           2. NA quantiles are reserved for bins that were excluded
#              due to either alignability or extremely high levels of
#              read mapping.
# granges - add new BED data to a previously constructed granges object.
#           Useful for constructing multi-feature enrichment analyses.
# has.metaline - the first line of the BED file is a #-prefixed single
#           line of KEY=VALUE pairs of metadata
read.bed <- function(bedfile, genome, granges.extend, feature.name=ifelse(is.qbed, 'quantile', 'feature'), add.chr.prefix=FALSE, remove.chr.prefix=TRUE, is.qbed=FALSE, has.metaline=FALSE) {
    # data.table's fread is often faster than read.table
    # make sure this new feature name won't collide with any old feature
    if (!missing(granges.extend)) {
        if (feature.name %in% attr(granges.extend, 'feature.set'))
            stop(paste('tried to add feature', feature.name, 'to feature set, but a feature of that name already exists'))
    }

    cat('Reading BED', bedfile, '\n')
    bed <- data.table::fread(bedfile, skip=ifelse(is.qbed | has.metaline, 1, 0))

    # Allow 3 or 4 column BEDs for convenience. For a 3 column
    # BED, assume that all intervals define a single class.
    # allow 3 column bed for convenience; assume there is only one state
    outside.name <- 'outside'
    if (ncol(bed) == 3) {
        bed <- cbind(bed, 'inside')
    }
    if (ncol(bed) > 4 & !is.qbed)
        cat("WARNING: BED contains >4 columns, ignoring columns 5 and greater.\n")

    colnames(bed)[1:4] <- c('chr', 'start', 'end', feature.name)

    # Remove 'chr' if requested
    if (add.chr.prefix)
        bed[['chr']] <- paste0('chr', bed[['chr']])
    if (remove.chr.prefix)
        bed[['chr']] <- sub('chr', '', bed[['chr']])

    # Make sure the BED file (with 'chr' prefix potentially stripped already)
    # is compatible with the genome the user asked for.
    # Don't allow warnings about seqinfo to get passed over.
    # Make a temporary GRanges object that will generate an error if the BED
    # and genome do not match seqlevels.
    # Can't do this in the second call below because an error will also be
    # generated if the BED file has any intervals that are longer than any
    # of the genome's contigs.
    options(warn=2)
    # Dummy construction, just to check seqinfo compatibility
    zzz <- GenomicRanges::GRanges(seqnames=bed[['chr']],
        ranges=IRanges::IRanges(start=bed[['start']], end=bed[['end']]))
    GenomeInfoDb::checkCompatibleSeqinfo(zzz, genome)
    options(warn=0)

    # Real 
    granges <- GenomicRanges::GRanges(seqnames=bed[['chr']],
        ranges=IRanges::IRanges(start=bed[['start']], end=bed[['end']]),
        seqinfo=GenomeInfoDb::seqinfo(genome))
    granges <- GenomicRanges::trim(granges)  # restrict intervals to the specified genome's contig sizes

    if (is.qbed) {
        mcols(granges)[[feature.name]] <- ifelse(is.na(bed[[feature.name]]),
            'excluded', bed[[feature.name]])
    } else {
        mcols(granges)[[feature.name]] <- bed[[feature.name]]
    }

    # Get the complement of the BED file's coverage.
    # setdiff - similar to gaps(), but does not add intervals for +/- strand
    outside <- GenomicRanges::setdiff(as(GenomeInfoDb::seqinfo(granges), 'GRanges'), granges)
    mcols(outside)[[feature.name]] <- outside.name
    granges <- sort(c(granges, outside))
    mcols(granges)[[feature.name]] <- as.factor(mcols(granges)[[feature.name]])

    # Now that the GRanges object is fully constructed, compare it to
    # the user-supplied one to ensure they are identical before merging.
    if (!missing(granges.extend)) {
        cat('Sorting new GRanges\n')
        granges <- sort(GenomeInfoDb::sortSeqlevels(granges))
        if (!all(seqnames(granges) == seqnames(granges.extend)) |
            !all(start(granges) == start(granges.extend)) |
            !all(end(granges) == end(granges.extend)) |
            !all(strand(granges) == strand(granges.extend))) {
                cat('----- old GRanges ------\n')
                print(granges.extend)
                cat('----- new GRanges ------\n')
                print(granges)

                stop(paste('tried to extend GRanges object by BED file',
                    bedfile, ', but BED file has different intervals'))
        }
        mcols(granges.extend)[[feature.name]] <- mcols(granges)[[feature.name]]
        granges <- granges.extend
    } else {
        cat("Precomputing GNCList for BED.\n")
        granges <- GenomicRanges::GNCList(sort(GenomeInfoDb::sortSeqlevels(granges)))
    }

    attr(granges, 'feature.set') <- c(attr(granges, 'feature.set'), feature.name)
    if (is.qbed | has.metaline) {
        meta <- read.bed.metadata(bedfile, is.qbed=is.qbed)
        cat("Got metadata: ")
        print(meta)
        if (!missing(granges.extend))
            meta <- rbind(attr(granges, 'bed.metadata'), meta)
        attr(granges, 'bed.metadata') <- meta
    } else {
        cat("No metadata to read.\n")
    }
    granges 
}

# Convenience function for translating mutation tables into GRanges
gr <- function(muts, seqinfo=NULL, add.chr.prefix=FALSE) {
    # Only add 'chr' if (1) the user requested it and (2) it isn't already there
    seqs <- muts$chr
    if (add.chr.prefix)
        seqs <- ifelse(grepl('chr', seqs), seqs, paste0('chr', seqs))

    GenomicRanges::GRanges(seqnames=seqs,
        ranges=IRanges::IRanges(start=muts$pos, width=1),
        seqinfo=seqinfo)
}

# grl is a GRangesList. This is critical for speed.
# verbose:
#   - 0: no messages
#   - 1: progress bar
#   - 2: extremely verbose memory profiling and map.and.count debugging
make.perm.matrix <- function(grl, edata, verbose=1) {
    if (verbose > 0) cat(sprintf("Analyzing %d permutations: [", length(grl)))
    ret <- map.and.count(grl, edata=edata, verbose=verbose)
    if (verbose > 0) cat('] 100%')
    ret
}

# bootstrap by resampling from the real mutation set. This is
# equivalent to just resampling the labels.
bootstrap <- function(gmuts, edata, n.boot, verbose=1) {
    if (verbose > 0) cat(sprintf('Bootstrapping %d samples: [', n.boot))
    if (verbose == 2) { cat('gc3a\n'); print(gc()) }

    ret <- do.call(cbind, lapply(feature.set(edata=edata), function(feat.name) {
        real.labels <- edata$map.fn(GenomicRanges::GRangesList(gmuts),
            edata=edata, feat.name=feat.name)$to

        # the combined dataset has to have an ID called 'perm.id' for
        # map.and.count to group the mapping prior to counting.
        if (verbose == 2) cat('sampling', length(gmuts)*n.boot, 'for bootstrapping\n')
        boot.obs <- data.table::data.table(t(sapply(1:n.boot, function(i) {
            if (i %% (n.boot/10) == 0 & verbose==1) cat('.')
            s <- sample(real.labels, size=length(real.labels), replace=TRUE)
            edata$count.fn(s, edata, feat.name)
        })))
        if (verbose == 2) { cat('gc3b\n'); print(gc()) }
        boot.obs
    }))
    if (verbose > 0) cat('] 100%')
    ret
}

# Either standardize measurements (as a z-score) or compute the
# enrichment of each observation (including permutations) relative
# to the mean count of permutations.
transform <- function(e, standardize=FALSE) {
    means <- colMeans(e$perm.mat)
    sds <- apply(e$perm.mat, 2, sd)
    if (standardize) {
        bmat <- apply(e$perm.mat, 2, function(col) (col-mean(col)) / sd(col))
        breal <- (e$real.obs - means) / sds
    } else {
        bmat <- apply(e$perm.mat, 2, function(col) col/mean(col))
        breal <- e$real.obs / means
    }
    ret <- list(real=as.vector(breal), pmat=bmat)
    if ('boots' %in% names(e)) {
        if (standardize)
            bs <- sapply(1:length(means), function(i)
                (e$boots[[i]] - means[i]) / sds[i])
        else
            bs <- sapply(1:length(means), function(i)
                e$boots[[i]] / means[i])
        ret <- c(ret, list(boots=bs))
    }

    ret
}

# Compute p-value of enrichment by determining how often permuted
# counts are further from 1 than the observed count.
# min.pval - for permutation tests, our p-value resolution is limited
#    to 1/(# permutations). So rather than claim something has 0 pvalue,
#    just report the minimum measurable p-value.
pval <- function(e, min.pval=1/nrow(e$perm.mat)) {
    te <- transform(e)
    lre <- log2(te$real)
    lpm <- log2(te$pmat)
    pmax(min.pval, sapply(1:length(lre), function(i) mean(abs(lpm[,i]) >= abs(lre[i]))))
}


# round.digits can make interactive use more convenient, but when writing
# to files we may want full precision.
# set round.digits=Inf to remove rounding.
esummary <- function(e, bootstrap.ci=0.95, round.digits=3) {
    pv <- pval(e)
    cmeans <- colMeans(e$perm.mat)

    stats <- t(apply(e$perm.mat, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975)))
    colnames(stats) <- paste0('perm.', c('lb.95', 'q1', 'med', 'q3', 'ub.95'))

    enr.stats <- t(apply(e$perm.mat, 2, function(col) quantile(col/mean(col), probs=c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=T)))
    colnames(enr.stats) <- paste0('enr.perm.', c('lb.95', 'q1', 'med', 'q3', 'ub.95'))

    sumdf <- data.frame(
        pval=pv,
        padj=p.adjust(pv, method='holm'),
        fdr=round(p.adjust(pv, method='fdr'),4),
        enr=round(e$real.obs / cmeans, round.digits),
        obs=e$real.obs,
        perm.mean=cmeans,
        stats,
        enr.stats
    )
    # add bootstrapping stats: number of bootstraps, 95% LB and UB
    if ('boots' %in% names(e)) {
        ci <- bootstrap.ci(e, ci=bootstrap.ci)
        sumdf$n.bootstraps <- nrow(e$boots)
        lb.name <- paste0('boot.', as.character(bootstrap.ci), '.lb')
        ub.name <- paste0('boot.', as.character(bootstrap.ci), '.ub')
        sumdf[[lb.name]] <- ci[1,]
        sumdf[[ub.name]] <- ci[2,]
        lb2.name <- paste0('enr.boot.', as.character(bootstrap.ci), '.lb')
        ub2.name <- paste0('enr.boot.', as.character(bootstrap.ci), '.ub')
        sumdf[[lb2.name]] <- round(ci[1,]/cmeans, round.digits)
        sumdf[[ub2.name]] <- round(ci[2,]/cmeans, round.digits)
    }
    if (!is.null(names(e$real.obs)))
        rownames(sumdf) <- names(e$real.obs)
    sumdf
}

bootstrap.ci <- function(e, ci=0.95) {
    alpha=1-ci
    apply(e$boots, 2, quantile, probs=c(alpha/2, 1-alpha/2))
}

# hack for plotting things in a different order
# WARNING: reordering the matrices breaks the ordering between
# levels(features) and the matrices and observations.
reorder <- function(en, new.order) {
    list(perm.mat=en$perm.mat[,..new.order],
         real.obs=en$real.obs[new.order],
         boots=en$boots[,..new.order],
         edata=en$edata)
}

# allow for subsetting an enrichment object to the set of features 'fts'
subset <- function(en, fts) {
    if ('boots' %in% names(en)) {
        list(
            perm.mat=en$perm.mat[,..fts],
            real.obs=en$real.obs[fts],
            boots=en$boots[,..fts],
            edata=en$edata
        )
    } else {
        list(
            perm.mat=en$perm.mat[,..fts],
            real.obs=en$real.obs[fts],
            edata=en$edata
        )
    }
}

# combine classes into a single class.
# this is done by simply summing up observations in permutations, real
# observations and bootstraps.
# classes.to.collapse - list where each element name denotes the new
#     class to create and each element is a character vector of classes
#     to merge together to create the new class.
#     We do not require that each new class be disjoint from other classes,
#     but some other functions do assume this (e.g., compute.outside).
# keep.old.classes - remove ALL classes in classes.to.collapse. New classes
#     are created before removing any classes in classes.to.collapse.
collapse <- function(en, classes.to.collapse=list(), keep.old.classes=FALSE) {
    if (any(names(classes.to.collapse %in% names(en$real.obs))))
        stop('new classes created by classes.to.collapse must not have the same name as any class already in e')

    # Add the new classes. Don't remove anything here because
    # classes.to.collapse may reuse the same class multiple times.
    for (i in 1:length(classes.to.collapse)) {
        newclass <- names(classes.to.collapse)[i]
        sumclasses <- classes.to.collapse[[i]]

        en$perm.mat[[newclass]] <- rowSums(en$perm.mat[,sumclasses,drop=F])
        en$real.obs[newclass] <- sum(en$real.obs[sumclasses])
        en$boots[[newclass]] <- rowSums(en$boots[,sumclasses,drop=F])
    }
    # Now remove old classes if the user specified it
    if (!keep.old.classes) {
        classes.to.remove <- unlist(classes.to.collapse)
        en$perm.mat <- en$perm.mat[,!(colnames(en$perm.mat) %in% classes.to.remove),drop=F]
        en$real.obs <- en$real.obs[!(names(en$real.obs) %in% classes.to.remove)]
        en$boots <- en$boots[,!(colnames(en$boots) %in% classes.to.remove),drop=F]
    }

    # Update the gbed with new classes
    en$last.collapse <- classes.to.collapse
    en$last.gbed <- en$edata$gbed
    new.gbed <- do.call(c, lapply(1:length(classes.to.collapse), function(i) {
        newclass <- names(classes.to.collapse)[i]
        oldclasses <- classes.to.collapse[[i]]
        gbed <- en$edata$gbed
        fname <- attr(gbed, 'feature.name')
        g <- gbed[mcols(gbed)[[fname]] %in% oldclasses,]
        g$feature <- newclass
        # Currently a GNCList. Need to coerce back to GRanges and recompute GNCList
        as(g, 'GRanges')
    }))
    new.gbed <- GenomicRanges::GNCList(sort(GenomeInfoDb::sortSeqlevels(new.gbed)))
    attr(new.gbed, 'feature.name') <- attr(en$last.gbed, 'feature.name')
    en$edata$gbed <- new.gbed
    en
}

# Plotting methods
# type:
#   'b' - plot lines and boxplots
#   'x' - plot boxplots only
#   'l' - plot lines only
# bootstrap=FALSE suppresses bootstrap plotting
plot.enrich <- function(e, es, ylim, ylab='Enrichment or depletion (%)', add=FALSE,
    bootstrap.ci=0.95, type=c('b', 'x', 'l'), las=3, lcol=1, lwd=2, ltype='b', show.asterisks=FALSE, cex.asterisk=1, ...)
{
    if (!missing(e) & !missing(es))
        stop("only one of 'e' or 'es' may be specified")
    if (missing(es))
        es <- esummary(e)

    type <- match.arg(type)

    if (missing(ylim)) {
        ylim <- get.ylim(es=es, use.enr=type == 'b' | type == 'l',
            use.boot=bootstrap.ci != FALSE,
            use.perm=type == 'b' | type == 'x',
            bootstrap.ci=bootstrap.ci)
    }

    bp <- 1:nrow(es)
    if (type == 'b' | type == 'x') {
        bxp(z=list(n=es$obs,
            names=rownames(es),
            stats=t(es[,paste0('enr.perm.', c('lb.95', 'q1', 'med', 'q3', 'ub.95'))])),
            lty='solid', ylim=ylim, ylab=ylab, border='grey', las=las,
            show.names=TRUE, ...)
    }
    if (type == 'b' | type == 'l') {
        plotf <- if (type == 'b' | add) lines else plot
        plotf(bp, es$enr, type=ltype, lwd=lwd, pch=20, ylim=ylim, col=lcol, ylab=ylab, las=las, ...)
        ast.locs <- es$enr
        if (bootstrap.ci != FALSE) {
            bs <- es[,c('enr.boot.0.95.lb','enr.boot.0.95.ub')]
            arrows(x0=bp, y0=bs[,1], x1=bp, y1=bs[,2],
                angle=90, code=3, length=0, col=lcol, lwd=2)
            ast.locs <- bs[,2]
        }

        # Asterisks are only printed if the signal is plotted (boxplots
        # only show the background distribution.)
        if (show.asterisks) {
            n.asts <- floor(-log10(es$pval))  # using pval because FDR is extremely skewed in the large joint analyses
            # e.g., when there are 0 muts mapped to a region the p-value will be NaN
            n.asts[is.na(n.asts)] <- 0
            asts <- sapply(n.asts,
                function(i) sprintf("%s", paste0(rep('*', i), collapse='')))
            text(x=bp, y=ast.locs, labels=asts, adj=c(0.5,0), cex=cex.asterisk)
        }
    }
}

# scale=100, shift=1 - for centering enrichment around 0
# use.boot - increase range to include bootstrap error bars
# use.perm - increase range to include permutation (null hypothesis) error bars
get.ylim <- function(e, es, scale=1, shift=0, use.enr=FALSE, use.boot=FALSE, use.perm=FALSE, bootstrap.ci=0.95) {
    if (missing(es) & !missing(e))
        es <- esummary(e, bootstrap.ci=bootstrap.ci)

    vals <- c()
    if (use.enr)
        vals <- c(vals, es$enr)
    if (use.boot)
        vals <- c(vals, unlist(es[,c('enr.boot.0.95.lb','enr.boot.0.95.ub')]))
    if (use.perm)
        vals <- c(vals, unlist(es[,c('enr.perm.lb.95', 'enr.perm.ub.95')]))

    r <- range(pretty(scale*(vals-shift)), na.rm=T)
    r*c(0.95,1.05)
}

# to change the order of bars in the plot, use reorder() on e first
# es - sometimes computing esummary() is expensive, so allow caller to
#      pass a precomputed summary.
barplot.enrich <- function(e, es=esummary(e), barcol='#4c4c4c', ylim, ...) {
    if (!missing(e) & !missing(es))
        stop("both 'e' and 'es' cannot be specified simultaneously")
    if (missing(es))
        es <- esummary(e)
    if (missing(ylim) & !missing(e))
        ylim <- get.ylim(e=e, scale=100, shift=1, use.boot=TRUE)
    if (missing(ylim) & !missing(es))
        ylim <- get.ylim(es=es, scale=100, shift=1, use.boot=TRUE)

    bp <- barplot(100*(es$enr-1),
        names.arg=rownames(es),
        las=3, col=barcol,
        border=NA,  ylim=ylim,
        ylab='Enrichment or depletion (%)', ...)
    abline(h=0)

    bs <- es[,c('enr.boot.0.95.lb','enr.boot.0.95.ub')]
    arrows(x0=bp, y0=100*(bs[[1]]-1), x1=bp, y1=100*(bs[[2]]-1),
        angle=90, code=3, length=0.05,lwd=1)
}

##################################################################################
# Methods for making Volcano-style plots out of multiple enrichment analyses
##################################################################################

# load an enrichment object file and do some additional things:
#    1. reorder the features
#    2. if a compute.outside > 0, then add the 'outside' interval
read.e <- function(file, check.features=NULL) {
    print(file)
    objs <- load(file)  # for FULL files, returns e, es, emeta
                        # for SUMMARY files, returns just es, emeta
    # horrible hack, but just go with it
    # attempt to automatically recognize full and summary files
    if ('e' %in% objs) {
        e <- reorder(e, order(unique(features(e))))
        return(e)
    } else {
        # also a horrible hack: assumes only 1 enrichment analysis in the file
        if (length(es) > 1)
            stop(paste('this function assumes only a single esummary exists, but this file contains', length(es)))
        e <- es[[1]]
    }

    if (!is.null(check.features)) {
        if (!all(unique(features(e)) %in% check.features))
            stop(paste('features in file', file, 'were not in the specified feature set'))
    }

    e
}


# just reads and reorders the enrichment objects and ensures
# that all enrichment objects are compatible (e.g., were run
# on the same feature sets).
read.es <- function(files) {
    # load just the first file to get the feature set.
    init.e <- read.e(files[1])
    #e <- get(load(files[1])[1])
    feature.set <- unique(features(init.e))
    es <- lapply(files, function(fn) {
        read.e(fn, check.features=feature.set)
    })
}


# map each enrichment object to a single (enrichment, p-value)
# pair. this requires choosing a specific feature on which to
# compute enrichment.
# produces a matrix suitable for evolcano()
# min.pval - should be 1 / (# of permutations used). Prevents
#    Inf p-values when none of the permutations exceeded the
#    observed values.
volcanoize <- function(es, feature='in peak', min.pval=1e-4) {
    # doing feature by name
    if (is.character(feature)) {
        if (!(feature %in% unique(features(es[[1]]))))
            stop(paste0("feature '", feature, "' not in first enrichment object"))
    } else if (is.integer(feature)) {
        if (feature > length(unique(features(es[[1]]))))
            stop(paste0("feature index ", feature, " exceeds the length of the first enrichment object"))
    }

    ret <- sapply(es, function(e) {
        # asssume data.frame means e is an esummary object
        if (is.data.frame(e))
            esum <- e
        else
            esum <- esummary(e)
        ret <- unlist(esum[feature, c('enr', 'pval')])
        ret[2] <- max(min.pval, ret[2])
        ret
    })
    ret[2,] <- -log10(ret[2,] + min.pval)
    colnames(ret) <- names(es)
    rownames(ret) <- c('enr','sig')
    t(ret)
}


evolcano <- function(emat, labels=rep('', nrow(emat)), annotate.labels=TRUE, cex=2, ...) {
    z <- emat
    # points will be red if they're labeled
    plot(z[,'enr'], z[,'sig'], pch=16,
        col=ifelse(z[,'sig'] <= 1 & labels == '', 'grey', 1 + (labels != '')),
            ylim=c(0,4.5), xlab='Enrichment (or depletion) ratio (obs/exp)',
            ylab='Significance: -log10(p-value)', cex=cex, ...)
    # overplot red points so they're always on the top layer
    points(z[,'enr'], z[,'sig'], pch=16,
        col=ifelse(labels != '',2, NA),
            ylim=c(0,4.5), xlab='Enrichment (or depletion) ratio (obs/exp)',
            ylab='Significance: -log10(p-value)', cex=cex, ...)

    if (any(labels != '') & annotate.labels) {
        z <- z[labels != '',]
        basicPlotteR::addTextLabels(xCoords=z[,'enr'], yCoords=z[,'sig'],
                    labels=labels[labels != ''], col.line=2,
                    col.label=2, cex.label=0.9)
    }
    abline(h=1, lty='dotted')
    abline(v=1, lty='dotted')
}

##################################################################################
# For making boxplots using bootstraps
##################################################################################
e.to.bootl <- function(es, feature.idx=2, labels=names(es)) {
    ret <- list(
        enrs=lapply(es, function(e) e$boots[,feature.idx]/mean(e$perm[,feature.idx])),
        pvals=pmax(1e-4,sapply(es, function(e) esummary(e)[feature.idx,]$pval)),
        fdrs=pmax(1e-4, sapply(es, function(e) esummary(e)[feature.idx,]$fdr)),
        labels=labels
    )
    ret
}

# merge bootls by:
#   1. using median enrichment values as the points for the boxplot and
#   2. using harmonic mean p-values to combine the p-values (for asterisks)
# CAUTION! CAUTION! p-values and FDRs are no longer different after applying
# this function. DO NOT interpret the 'fdrs' output of this function as FDRs.
merge.bootl <- function(bootl, group.labels) {
    harmonicmeanpvs <- sapply(split(bootl$pvals, group.labels),
        function(ps) harmonicmeanp::p.hmp(ps, L=length(ps)))
    list(
        enrs=lapply(split(bootl$enrs, group.labels), function(subenrs) sapply(subenrs, median)),
        pvals=harmonicmeanpvs,
        fdrs=harmonicmeanpvs,
        labels=group.labels[!duplicated(group.labels)])
        
}

boxp <- function(l, add.points=TRUE, asts.at=NA, ...) {
    bp <- boxplot(l$enrs, names=l$labels,
        outline=F, lty='solid', pars=list(staplewex=0),
        las=3, ...)
    abline(h=1, lty='dotted')
    if (!is.na(asts.at)) {
        asts <- sapply(floor(-log10(l$fdrs)),
            function(i) sprintf("%s", paste0(rep('*', i), collapse='')))
        print(asts)
        #text(x=1:length(l$enrs), y=asts.at, labels=asts)
        text(x=1:length(l$enrs), y=bp$stats[5,], labels=asts, adj=c(0.5,0))
    }
    if (add.points) {
        stripchart(l$enrs, vertical=T, pch=20, add=TRUE, method='jitter')
    }
}



# A standard command line analysis:
#   Given a .rda file containing mutations and a corresponding
#   set of random permutations, determine the enrichment of
#   mutations over some set of genomic regions defined by the object
#   returned by 'init.enrich'. Results are written to the output .rda
#   files.
#
#   Significance is assessed by bootstrapping with `n.boot`
#   resamplings and comparison to permutations.
#
# init.enrich - is a function that takes BSgenome object followed by
#     an arbitrary number of
#     command line arguments (starting with arg 5) and returns an
#     enrichment data object.
# args - allow the caller to override args. this allows the caller to
#     implement different command line options and then reshape the
#     args in the way command.line.analysis expects.
# genome - the name of a locally installed BSgenome
command.line.analysis <- function(init.enrich, genome, args=commandArgs(trailingOnly=TRUE)) {
    if (length(args) < 6) {
        cat("Got command args:\n")
        print(args)
        stop("usage: [enrichment_analysis.r] n.bootstraps muts.rda[:varname1] perms.rda[:use_N][:varname2] full_output.rda summary_output.rda input1 [ input2 ... inputN ]")
    }

    n.boot <- as.integer(args[1])
    mutarg <- args[2]
    permarg <- args[3]
    fulloutfile <- args[4]
    summaryoutfile <- args[5]
    inputdata <- args[6:length(args)]

    # Parse file containing mutation table. Optionally allow user to
    # specify what the name of the variable in the file is, in case it
    # contains multiple variables.
    mut.elts <- strsplit(mutarg, ":")[[1]]
    mutfile <- mut.elts[1]
    mut.varname <- NA
    if (length(mut.elts) == 2)
        mut.varname <- mut.elts[2]

    perm.elts <- strsplit(permarg, ":")[[1]]
    permfile <- perm.elts[1]
    perms.useN <- Inf
    perm.varname <- NA
    if (length(perm.elts) > 1)
        perms.useN <- as.integer(perm.elts[2])
    if (length(perm.elts) > 2)
        perm.varname <- perm.elts[3]

    if (file.exists(fulloutfile))
        stop(sprintf("output file %s already exists, please delete it first",
            fulloutfile))
    if (file.exists(summaryoutfile))
        stop(sprintf("output file %s already exists, please delete it first",
            summaryoutfile))

    cat("Getting genome", genome, "\n")
    genome <- BSgenome::getBSgenome(genome)

    eobject <- init.enrich(genome, inputdata)

    ret <- load(mutfile, verb=T)
    if (is.na(mut.varname)) {
        if (length(ret) == 1)
            mut.varname <- ret
        else
            stop(paste('no variable name provided for mut.rda and file contains multiple objects:', ret))
    }
    muts <- get(mut.varname)
    gmuts <- gr(muts, seqinfo=seqinfo(genome), add.chr.prefix=TRUE)
    gmuts$perm.id <- 1
    # Automatically recognize ID83/SBS96 signatures from our data
    if (!is.na(eobject$use.mutclass)) {
        cat('user specified mutation signature column name', eobject$use.mutclass, '\n')
        # a bunch of hackery to recognize our specific data structures
        # SNVs use this colname
        if ('type.and.ctx' %in% colnames(muts)) {
            cat('detected SNV mutation signatures in "type.and.ctx" column; adding to mutation GRanges\n')
            mcols(gmuts)[[eobject$use.mutclass]] <- scan2::sbs96(muts$type.and.ctx)
        }
        # indels use this one
        else if ('muttype' %in% colnames(muts)) {
            cat('detected indel mutation signatures in "muttype" column; adding to mutation GRanges\n')
            mcols(gmuts)[[eobject$use.mutclass]] <- scan2::id83(muts$muttype)
        } else {
            stop('expected either type.and.ctx (SNVs) or muttype (indels) in mutation object\n')
        }
    }

    ret <- load(permfile, verb=T)
    if (is.na(perm.varname)) {
        if (length(ret) == 1)
            perm.varname <- ret
        else
            stop(paste('no variable name provided for mut.rda and file contains multiple objects:', ret))
    }
    zperml <- get(perm.varname)
    seqlevels(zperml) <- paste0('chr', seqlevels(zperml))
    # remove the extra reference so memory can be easily freed later
    # corner case: unless the variable name is 'zperml', which would remove
    # the reference we just created above.
    if (perm.varname != 'zperml')
        rm(list=perm.varname)
    
    # If the user requests, trim down the permutation set
    actual.useN <- max(min(length(zperml), perms.useN), 1)
    if (actual.useN < length(zperml)) {
        cat(sprintf("Using the first %d permutations out of %d total (useN=%g)\n", 
            actual.useN, length(zperml), perms.useN))
        cat("WARNING: this creates a temporary copy in RAM of the subsetted\n")
        cat("permutations. Do not use this option if RAM usage is sensitive.\n")
        zperml <- zperml[1:actual.useN,]
    }

    # 'chr' prefixes are becoming such an issue that we need to check
    # before waiting an hour for this to run.
    # checkCompatibleSeqinfo throws a warning when there are no
    # sequence levels in common. Need to fail fast on this.
    options(warn=2)
    GenomeInfoDb::checkCompatibleSeqinfo(gmuts, eobject$gbed)
    GenomeInfoDb::checkCompatibleSeqinfo(gmuts, zperml)
    options(warn=0)

    gccheck <- function() {
        g <- gc()
        paste0('Mb used: ', sum(g[,2]), ', Peak Mb used: ', sum(g[,6]))
    }

    cat('RAM usage before analysis |', gccheck(), '\n')

    # Step 1. permutations
    perm.mat <- make.perm.matrix(zperml, edata=eobject, verbose=1)
    # remove massive permutation set to conserve memory
    rm(zperml)
    cat(' RAM usage |', gccheck(), '\n')

    # Step 2. real mutations
    real.obs <- unlist(map.and.count(GenomicRanges::GRangesList(gmuts),
                                     edata=eobject, verbose=0))

    # Step 3. bootstrapping
    if (n.boot > 0) {
        boots <- bootstrap(gmuts, edata=eobject, n.boot=n.boot)
        cat(' RAM usage |', gccheck(), '\n')
        e <- list(feature.set=feature.set(edata=eobject),
            perm.mat=perm.mat, real.obs=real.obs, edata=eobject, boots=boots)
    } else {
        cat(' RAM usage |', gccheck(), '\n')
        e <- list(feature.set=feature.set(edata=eobject),
            perm.mat=perm.mat, real.obs=real.obs, edata=eobject)
    }


    # Step 4. precompute esummary table using defaults like bootstrap CI=95%
    es <- lapply(feature.set(edata=eobject), function(feat.name) {
        this.e <- get.feat(e, feat.name)
        esummary(this.e)
    })
    names(es) <- feature.set(edata=eobject)

    list(fulloutfile=fulloutfile, e=e,
         summaryoutfile=summaryoutfile, es=es)
}
