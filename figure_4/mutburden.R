#!/usr/bin/env Rscript


library(argparse)

parser <- ArgumentParser()

parser$add_argument("mutations", type="character",
    help='CSV file containing a table of mutations annotated per sample.')
parser$add_argument("germline_control", type='character',
    help='CSV file containing heterozygous SNP sites subjected to leave-1-out SCAN2 calling. Used to estimate somatic sensitivity.')
parser$add_argument("callable_tables", type='character',
    help='RDA file containing one sequencing depth table for each sample in mutations. Each cell (i,j) in the table should contain the number of bases with depth=i+1 in the single cell and depth=j+1 in the bulk. The R object must be a list of such tables with each list element named for its sample. The R object can be named anything, but must be the only object in the RDA.')
parser$add_argument("output_summary", type='character',
    help='Write metadata joined to mutation burden estimates (in CSV format) to this file.')
parser$add_argument('--tag', type='character', metavar='STRING', default=NULL,
    help='Optional identifier tag to add as a first column in both outputs. This can be helpful when merging multiple tables (e.g., from different projects or SNV/indels).')
parser$add_argument("--metadata", type='character', metavar='FILE_PATH',
    help='Optional CSV file containing metadata to which mutation burden estimates will be joined. At minimum, must include a single "sample" column identifying each single cell. This can also specify a subset of samples present in the mutations file.  If no metadata is provided, all samples in the mutations file will be analyzed.')
parser$add_argument("--load-priors", default=FALSE, action='store_true',
    help='[NOT IMPLEMENTED] Read upper bounds on mutation burden as determined by the PRE-GENOTYPING step.')
parser$add_argument("--gbp-per-genome", default=5.845001134, type='double',
    help='Number of haploid basepairs (in billions) per genome. The default value of 5.845001134 corresponds to AUTOSOMES ONLY as determined by GRCh37.')
parser$add_argument("--min-bulk-dp", default=11, type='integer',
    help="Minimum required bulk depth for SCAN2 calling. Must match the values used for calling mutations and generating callable regions.")
parser$add_argument("--min-sc-dp", default=6, type='integer',
    help="Minimum required single cell depth for SCAN2 calling. Must match the values used for calling mutations and generating callable regions.")
parser$add_argument("--min-sc-alt", default=2, type='integer',
    help="Minimum required single cell alt read count for SCAN2 calling. Must match the values used for calling mutations and generating callable regions.")


# detect script being run by snakemake
# if so, make a mock commandArgs function
if ('snakemake' %in% ls()) {
    logfile <- snakemake@log[[1]]
    con <- file(logfile, 'w')
    sink(con, type='output')
    sink(con, type='message')

    commandArgs <- function(...) unlist(snakemake@params)
    cat('Intercepted command line arguments from snakemake:\n')
    print(commandArgs())
}

args <- parser$parse_args(commandArgs(trailingOnly=TRUE))

if (file.exists(args$output_summary))
    stop(paste('output file', args$output_summary, 'already exists, please delete it first'))


suppressMessages(library(data.table))  # read large tables like the germline control sites
suppressMessages(library(lme4))
suppressMessages(library(lmerTest)) # adds p-value estimate to lme4 output

cat("Reading mutation table..\n")
muts <- fread(args$mutations)
setkey(muts, sample)
muts$pass <- muts$passA

# Allows MDA analysis for paper. Will be removed in later versions.
if ("passM" %in% colnames(muts))
    muts$pass <- muts$pass | muts$passM

cat("Reading list of depth tables..\n")
dptabs <- get(load(args$callable_tables))

# allow the case where only one dptab (rather than a list of
# collected dptabs) is supplied for compatibility with raw SCAN2 output
if (!is.list(dptabs) & (length(unique(muts$sample)) == 1 | nrow(muts) == 0)) {
    cat("either no mutations or 1 sample in somatic mut table - assuming dptab is not a list\n")
    dptabs <- list(dptabs)
    names(dptabs) <- unique(muts$sample)[1]
}
cat(paste('Got', length(dptabs), 'depth tables with dim', nrow(dptabs[[1]]), "x", ncol(dptabs[[1]]), '\n'))


cat("Reading germline control table..\n")
germline <- fread(args$germline_control)
setkey(germline, sample)

if (nrow(germline) == 0)
    stop('got 0 germline control sites, cannot continue')

# Just allows for scenarios where no mutations are called at all.
# I.e., when analyzing a single sample (especially for indels).
if (nrow(muts) == 0) {
    if (length(unique(germline$sample)) > 1) {
        stop('0 somatic mutations in table and multiple samples in germline control')
    } else {
        warning(paste0("somatic mutation table was empty. ASSUMING sample matches the single sample ", germline$sample[1], " from the germline control"))
    }
    warning('Generating dummy row of NAs for somatic table.')
    # This makes a single row of all NA (but of the correct internal data types)
    muts <- muts[, lapply(.SD, function(x) x[NA])]
    muts$sample[1] <- germline$sample[1]
    muts$pass[1] <- FALSE # so it'll be recognized as 0 passing mutations
    muts$dp[1] <- 0       # need non-NA dp values for dptest.mid
}

if (is.null(args$metadata)) {
    meta <- data.table(sample=unique(muts[['sample']]))
} else {
    cat("Reading metadata..\n")
    meta <- fread(args$metadata)
}
setkey(meta, sample)



# load the pre-genotyping burden estimate. This should bear on
# sensitivity, since it affects N_T and N_A estimates.
# XXX: Not incorporated for now. The issue is the need to access the
# full SCAN2 output directories.
if (args$load_priors) {
    stop('Importing prior mutation burdens is not yet implemented\n')
    # XXX: donor is no longer a required argument in the metadata table
    donors <- unique(meta[samples,]$donor)
    fcs <- do.call(c, lapply(donors, function(dn) {
        this.samples <- meta$sample[meta$donor==dn & meta$sample %in% samples]
        ret <- lapply(this.samples, function(sn) {
            f <- sprintf("~/ndata1/pta/%s/scansnv_fdr01_noX/indel/%s/fdr_tuning.rda", dn, sn)
            if (!file.exists(f))
                return(list(burden=c(NA,NA)))
            load(f)
            fdr.tuning
        })
        names(ret) <- this.samples
        ret
    }))
    pre.geno.burdens <- sapply(fcs, function(fc) fc$burden[2])
}

samples <- meta[['sample']]

if (FALSE) {
# Get sensitivity estimates from germline variants
sens <- germline[, .(sens=mean(pass),
                     callable.sens=mean(pass[bulk.dp >= args$min_bulk_dp & dp >= args$min_sc_dp])),
                   by=sample]
setkey(sens, sample)

# Get counts of mutations per sample - remember, do not use signature
# rescued (passB) mutations for this.
muttab <- muts[, .(nsom=sum(passA)), by=sample]
setkey(muttab, sample)


# Mash everything together
muttab <- meta[muttab, on='sample']
}

# Add optional metadata tag
if (!is.null(args$tag)) {
    muttab <- cbind(args$tag, muttab)
}

# Break data into 4 quantiles based on depth, use the middle 2 (i.e.,
# middle 50%) to reduce noise caused by very low and very high qs.
dptest.qmid <- function(s, g, dptab, q=4) {
    gpassname <- 'pass'

    if (q != 4) stop('for qmid, q=4 is required')
    qbreaks <- quantile(g$dp, prob=0:q/q)
    if (sum(qbreaks == 0) > 1) {
        print(qbreaks)
        warning('more than one depth quartile is 0. In deep sequencing data, this usually indicates poor amplification. Consider treating this sample as an outlier.')
        qbreaks[1] <- -1
    }

    s$dpq <- cut(s$dp, qbreaks, include.lowest=T, labels=F)  # BOTH use g for the dp quantiles
    s$dpq[s$dpq==3] <- 2 # make 25-75% middle 50% a single bin
    g$dpq <- cut(g$dp, qbreaks, include.lowest=T, labels=F)
    g$dpq[g$dpq==3] <- 2
    # cut down dptab to the max value in g$dp (+1 because 1 corresponds to dp=0)
    dptab <- dptab[1:min(max(g$dp)+1, nrow(dptab)),]

    rowqs <- cut(0:(nrow(dptab)-1), qbreaks, include.lowest=T, labels=F)
    rowqs[rowqs==3] <- 2
    qstouse <- c(1,2,4)
    splits <- lapply(qstouse, function(q) s[dpq==q])
    splitg <- lapply(qstouse, function(q) g[dpq==q])
    ret <- data.frame( 
        ncalls=sapply(splits, function(d) sum(d$pass[d$bulk.dp >= args$min_bulk_dp])),
        gsens=sapply(splitg, function(d) mean(d[[gpassname]][d$bulk.dp >= args$min_bulk_dp])),
        cbp=sapply(split(rowSums(dptab[,-(1:args$min_bulk_dp)]), rowqs), sum)
    )
    ret$callable.burden <- ret$ncalls / ret$gsens
    # dividing by 2 makes it haploid gb
    ret$rate.per.gb <- ret$callable.burden / ret$cbp * 1e9/2
    ret$burden <- ret$rate.per.gb * args$gbp_per_genome
    ret$genome.sens <- ret$ncalls / ret$burden
    ret$raw.calls <- sum(ret$ncalls)
    # [2,] corresponds to middle 50%; ret should be 3 rows
    ret[2,]
}

# Specific to publication, remove in later version.
meta <- meta[meta$amp != 'bulk' & meta$sample != '5087pfc-Rp3C5',]

newburden <- do.call(rbind, lapply(meta$sample, function(sn) {
    print(sn)
    dptq <- dptest.qmid(s=muts[sample==sn], g=germline[sample==sn],
        dptab=dptabs[[sn]])
}))

cat('Writing mutation burden table..\n')
fwrite(cbind(meta, newburden), file=args$output_summary)
warnings()
