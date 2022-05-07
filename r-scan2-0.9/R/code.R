mutsig.cols <- rep(c('deepskyblue', 'black', 'firebrick2', 'grey', 'chartreuse3', 'pink2'), each=16)

muttype.map <- c(
    'A>C'="T>G",
    'A>G'="T>C",
    'A>T'="T>A",
    'C>A'="C>A",
    'C>G'="C>G",
    'C>T'="C>T",
    'G>A'='C>T',
    'G>C'='C>G',
    'G>T'='C>A',
    'T>A'='T>A',
    'T>C'='T>C',
    'T>G'='T>G')

K.func <- function (x, y, a, b, c, d) 
exp(a - (x - y)^2/b^2) + exp(c - (x - y)^2/d^2)

abc2 <- function (altreads, gp.mu, gp.sd, dp, factor = 1) 
{
    pmf <- dreads(0:dp, d = dp, gp.mu = gp.mu, gp.sd = gp.sd, 
        factor = factor)
    p.value <- sum(pmf[pmf <= pmf[altreads + 1]])
    c(p.value = p.value)
}

abmodel.approx.ctx <- function (x, y, d, hsnp.chunksize = 100) 
{
    n <- hsnp.chunksize
    list(hsnp.chunksize = as.integer(hsnp.chunksize), x = x, 
        y = y, d = d, U = numeric(n), V = numeric(n), B = numeric(n), 
        sqrtW = numeric(n), K = numeric(n * n), A = numeric(n * 
            n))
}

abmodel.approx.logp <- function (a, b, c, d, ctx, max.it = as.integer(50), verbose = FALSE) 
{
    result <- .Call("laplace_approx_chunk_cpu", ctx$hsnp.chunksize, 
        c(a, b, c, d), ctx$x, ctx$y, ctx$d, as.integer(length(ctx$x)), 
        ctx$U, ctx$V, ctx$B, ctx$sqrtW, ctx$K, ctx$A, max.it, 
        verbose, PACKAGE = "scansnv")
    return(result)
}

abmodel.sample <- function (n = 1000, alim = c(-7, 2), blim = c(2, 4), clim = c(-7, 
    2), dlim = c(2, 6), ctx, seed = 0, max.it = 50, verbose = FALSE) 
{
    set.seed(seed)
    max.it <- as.integer(max.it)
    params <- data.frame(a = runif(n, min = alim[1], max = alim[2]), 
        b = 10^runif(n, min = blim[1], max = blim[2]), c = runif(n, 
            min = clim[1], max = clim[2]), d = 10^runif(n, min = dlim[1], 
            max = dlim[2]))
    logps <- mapply(abmodel.approx.logp, a = params$a, b = params$b, 
        c = params$c, d = params$d, MoreArgs = list(ctx = ctx, 
            max.it = max.it, verbose = verbose))
    return(cbind(params, logp = logps))
}

alg3.2.2 <- function (a, b, c, d, ctx, Xnew) 
{
    mode <- ctx$B[1:length(ctx$x)]
    K <- outer(ctx$x, ctx$x, K.func, a = a, b = b, c = c, d = d)
    covK <- outer(ctx$x, Xnew, K.func, a = a, b = b, c = c, d = d)
    mean.new <- t(covK) %*% (ctx$y - ctx$d * exp(mode)/(1 + exp(mode)))
    W <- ctx$d * exp(mode)/(1 + exp(mode))^2
    sqrtW <- sqrt(W)
    L <- t(chol(diag(length(ctx$x)) + outer(sqrtW, sqrtW) * K))
    v <- forwardsolve(L, outer(sqrtW, rep(1, ncol(covK))) * covK)
    cov.new <- outer(Xnew, Xnew, K.func, a = a, b = b, c = c, 
        d = d) - t(v) %*% v
    if (any(is.na(mean.new)) | any(is.na(cov.new))) 
        stop("NA values predicted")
    list(mean = mean.new, cov = cov.new)
}

apply.fdr.tuning.parameters <- function (somatic, fdr.tuning) 
{
    somatic$popbin <- ceiling(somatic$af * fdr.tuning$bins)
    somatic$popbin[somatic$dp == 0 | somatic$popbin == 0] <- 1
    nt.na <- mapply(function(dp, popbin) {
        idx = min(dp, fdr.tuning$max.dp + 1) + 1
        if (is.null(fdr.tuning$fcs[[idx]]$pops)) 
            c(0.1, 0.1)
        else fdr.tuning$fcs[[idx]]$pops$max[popbin, ]
    }, somatic$dp, somatic$popbin)
    rownames(nt.na) <- c("nt", "na")
    nt.na
}

bin.afs <- function (afs, bins = 20) 
{
    sq <- seq(0, 1, 1/bins)
    x <- findInterval(afs, sq, left.open = TRUE)
    x[x == 0] <- 1
    tx <- tabulate(x, nbins = bins)
    names(tx) <- apply(cbind(sq[-length(sq)], sq[-1]), 1, mean)
    tx
}

df.to.sbs96 <- function (df) 
{
    if (!("type.and.ctx" %in% colnames(df))) 
        df <- get.3mer(df)
    bases <- c("A", "C", "G", "T")
    t <- rep(0, 96)
    names(t) <- paste0(rep(bases, each = 4), rep(c("C", "T"), 
        each = 48), rep(bases, times = 4), ":", rep(c("C", "T"), 
        each = 48), ">", c(rep(c("A", "G", "T"), each = 16), 
        rep(c("A", "C", "G"), each = 16)))
    t2 <- table(df$type.and.ctx)
    t[names(t2)] <- t2
    tn <- do.call(rbind, strsplit(names(t), ":"))
    t <- t[order(tn[, 2])]
    t <- t + 0.1
    t/sum(t)
}

do.rescue <- function (ldf, lysis.sig, min.alt = 2, calling.fdr = 0.01, 
    rescue.fdr = 0.01) 
{
    if (missing(lysis.sig)) {
        data(pta.snv.artifact.sig.v1)
        lysis.sig <- pta.snv.artifact.sig.v1
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
        get.3mer(df)
    })
    names(ldf2) <- names(ldf)
    ldf <- ldf2
    cat(".\n")

    cat("    Step 2: constructing true SNV spectrum..\n")
    x <- do.call(rbind, lapply(ldf, function(s) s[s$pass, 1:5]))
    good.sig <- df.to.sbs96(x)

    cat(sprintf("    Step 3: adjusting FDR for lysis artifacts using new FDR target=%0.4f..\n        ", 
        rescue.fdr))
    final <- lapply(1:length(ldf), function(i) {
        df <- ldf[[i]]
        cat(sprintf(" %s", names(ldf)[i]))
        rescue(df, lysis.sig = lysis.sig, good.sig = good.sig, 
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

dreads <- function (ys, d, gp.mu, gp.sd, factor = 1, ghd = gaussHermiteData(128)) 
{
    require(fastGHQuad)
    sapply(ys, function(y) ghQuad(function(x) {
        b <- sqrt(2) * gp.sd * x + gp.mu
        exp(dbinom(y, size = d, prob = 1/(factor * (1 + exp(-b))), 
            log = TRUE) - log(pi)/2)
    }, ghd))
}

estimate.alphabeta3 <- function (gp.mu, gp.sd, dp = 30, div = 2) 
{
    td <- data.frame(dp = 0:dp, mut = dreads(0:dp, d = dp, gp.mu = gp.mu, 
        gp.sd = gp.sd), err1 = dreads(0:dp, d = dp, gp.mu = gp.mu, 
        gp.sd = gp.sd, factor = div), err2 = dreads(0:dp, d = dp, 
        gp.mu = -gp.mu, gp.sd = gp.sd, factor = div))
    td$err <- (td$err1 + td$err2)/2
    td <- td[order(td$err), ]
    td$cumerr <- cumsum(td$err)
    return(list(td = td, alphas = td$cumerr, betas = cumsum(td$mut)))
}

estimate.somatic.burden <- function (fc, min.s = 1, max.s = 5000, n.subpops = 10, display = FALSE, 
    rough.interval = 0.99) 
{
    sim <- function(n.muts, g, s, n.samples = 1000, diagnose = FALSE) {
        samples <- rmultinom(n = n.samples, size = n.muts, prob = g)
        if (diagnose) {
            boxplot(t(samples))
            lines(s, lwd = 2, col = 2)
        }
        mean(apply(samples, 2, function(col) all(col <= s)))
    }
    srange <- c(1:100, seq(101, max.s, length.out = n.subpops))
    srange <- srange[srange < max.s]
    fraction.embedded <- sapply(srange, sim, g = fc$g, s = fc$s)
    min.burden <- max(c(0, srange[fraction.embedded >= 1 - (1 - 
        rough.interval)/2]))
    max.burden <- max(c(0, srange[fraction.embedded >= (1 - rough.interval)/2]))
    c(min = min.burden, max = max.burden)
}

fcontrol <- function (germ.df, som.df, bins = 20, rough.interval = 0.99) 
{
    germ.afs <- germ.df$af[!is.na(germ.df$af) & germ.df$af > 
        0]
    som.afs <- som.df$af[!is.na(som.df$af)]
    g <- bin.afs(germ.afs, bins = bins)
    s <- bin.afs(som.afs, bins = bins)
    if (length(s) == 0 | all(g == 0)) 
        return(list(est.somatic.burden = c(0, 0), binmids = as.numeric(names(g)), 
            g = g, s = s, pops = NULL))
    approx.ns <- estimate.somatic.burden(fc = list(g = g, s = s), 
        min.s = 1, max.s = sum(s), n.subpops = min(sum(s), 100), 
        rough.interval = rough.interval)
    cat(sprintf("fcontrol: dp=%d, max.s=%d (%d), n.subpops=%d, min=%d, max=%d\n", 
        germ.df$dp[1], nrow(som.df), sum(s), min(nrow(som.df), 
            100), as.integer(approx.ns[1]), as.integer(approx.ns[2])))
    pops <- lapply(approx.ns, function(n) {
        nt <- pmax(n * (g/sum(g)), 0.1)
        na <- pmax(s - nt, 0.1)
        cbind(nt = nt, na = na)
    })
    return(list(est.somatic.burden = approx.ns, binmids = as.numeric(names(g)), 
        g = g, s = s, pops = pops))
}

find.nearest.germline <- function (som, germ, chrs = c(1:22, "X")) 
{
    som$nearest.het <- NA
    for (chr in chrs) {
        gpos <- germ$pos[germ$chr == chr]
        spos <- som$pos[som$chr == chr]
        gidx <- findInterval(spos, gpos)
        gidx[gidx == 0] <- 1
        nearest.idx <- ifelse(abs(gpos[gidx] - spos) <= abs(gpos[gidx + 
            1] - spos), gidx, gidx + 1)
        som$nearest.het[som$chr == chr] <- gpos[nearest.idx]
    }
    som
}

genotype.somatic <- function (gatk, gatk.lowmq, sc.idx, bulk.idx, sites.with.ab, 
    sc.cigars, bulk.cigars, cigar.training, cigar.emp.score, 
    fdr.tuning, spikein = FALSE, cap.alpha = TRUE, cg.id.q = 0.05, 
    cg.hs.q = 0.05, random.seed = 0, target.fdr = 0.1, bulkref = bulk.idx + 
        1, bulkalt = bulk.idx + 2, scref = sc.idx + 1, scalt = sc.idx + 
        2, min.sc.alt = 0, min.sc.dp = 0, min.bulk.dp = 0) 
{
    call.fingerprint <- as.list(environment())
    dont.save <- c("gatk", "gatk.lowmq", "sites.with.ab", "sc.cigars", 
        "bulk.cigars")
    call.fingerprint <- call.fingerprint[!(names(call.fingerprint) %in% 
        dont.save)]
    set.seed(random.seed)
    cat("step 1: preparing data\n")
    gatk$muttype <- muttype.map[paste(gatk$refnt, gatk$altnt, 
        sep = ">")]
    gatk$dp <- gatk[, scalt] + gatk[, scref]
    gatk$af <- gatk[, scalt]/gatk$dp
    gatk$bulk.dp <- gatk[, bulkalt] + gatk[, bulkref]
    somatic <- merge(gatk, sites.with.ab, all.y = TRUE)
    cat(sprintf("        %d somatic SNV candidates\n", nrow(somatic)))
    somatic$gp.mu <- ifelse(!is.na(somatic$af) & somatic$af < 
        1/2, -abs(somatic$gp.mu), abs(somatic$gp.mu))
    somatic$ab <- 1/(1 + exp(-somatic$gp.mu))
    cat("step 2: computing p-values for filters\n")
    cat("        allele balance consistency\n")
    somatic$abc.pv <- mapply(abc2, altreads = somatic[, scalt], 
        gp.mu = somatic$gp.mu, gp.sd = somatic$gp.sd, factor = 1, 
        dp = somatic$dp)
    cat("        lysis artifacts\n")
    somatic$lysis.pv <- mapply(test2, altreads = somatic[, scalt], 
        gp.mu = somatic$gp.mu, gp.sd = somatic$gp.sd, dp = somatic$dp, 
        div = 2)
    cat("        MDA artifacts\n")
    somatic$mda.pv <- mapply(test2, altreads = somatic[, scalt], 
        gp.mu = somatic$gp.mu, gp.sd = somatic$gp.sd, dp = somatic$dp, 
        div = 4)
    cat(sprintf("step 3: tuning FDR = %0.3f\n", target.fdr))
    if (cap.alpha) 
        cat(sprintf("        cap.alpha=TRUE: alpha <= %0.3f enforced despite artifact prevalence\n", 
            target.fdr))
    nt.na <- apply.fdr.tuning.parameters(somatic, fdr.tuning)
    somatic <- cbind(somatic, t(nt.na))
    cat("        lysis artifact FDR\n")
    somatic <- cbind(somatic, lysis = t(mapply(match.fdr3, pv = somatic$lysis.pv, 
        gp.mu = somatic$gp.mu, gp.sd = somatic$gp.sd, dp = somatic$dp, 
        nt = somatic$nt, na = somatic$na, div = 2)))
    cat("        MDA artifact FDR\n")
    somatic <- cbind(somatic, mda = t(mapply(match.fdr3, pv = somatic$mda.pv, 
        gp.mu = somatic$gp.mu, gp.sd = somatic$gp.sd, dp = somatic$dp, 
        nt = somatic$nt, na = somatic$na, div = 4)))
    cat("step 5: applying optional alignment filters\n")
    if (missing(gatk.lowmq)) {
        cat("        WARNING: skipping low MQ filters. will increase FP rate\n")
        somatic$lowmq.test <- TRUE
    }
    else {
        cat(sprintf("        attaching low MQ data..\n"))
        lmq <- gatk.lowmq[, c("chr", "pos", "refnt", "altnt", 
            colnames(gatk.lowmq)[c(scref, scalt, bulkref, bulkalt)])]
        somatic <- merge(somatic, lmq, by = c("chr", "pos", "refnt", 
            "altnt"), all.x = TRUE, suffixes = c("", ".lowmq"))
        cn <- paste0(colnames(gatk.lowmq)[bulkalt], ".lowmq")
        somatic$lowmq.test <- is.na(somatic[, cn]) | somatic[, 
            cn] == 0 | spikein
    }
    somatic <- merge(somatic, sc.cigars, by = c("chr", "pos"), 
        all.x = T)
    somatic <- merge(somatic, bulk.cigars, by = c("chr", "pos"), 
        all.x = T, suffixes = c("", ".bulk"))
    somatic$id.score.y <- somatic$ID.cigars/somatic$dp.cigars
    somatic$id.score.x <- somatic$ID.cigars.bulk/somatic$dp.cigars.bulk
    somatic$id.score <- cigar.emp.score(training = cigar.training, 
        test = somatic, which = "id")
    somatic$hs.score.y <- somatic$HS.cigars/somatic$dp.cigars
    somatic$hs.score.x <- somatic$HS.cigars.bulk/somatic$dp.cigars.bulk
    somatic$hs.score <- cigar.emp.score(training = cigar.training, 
        test = somatic, which = "hs")
    cat("        Excessive indel CIGAR ops\n")
    somatic$cigar.id.test <- somatic$id.score > quantile(cigar.training$id.score, 
        prob = cg.id.q, na.rm = T)
    cat("        Excessive clipped read CIGAR ops\n")
    somatic$cigar.hs.test <- somatic$hs.score > quantile(cigar.training$hs.score, 
        prob = cg.hs.q, na.rm = T)
    somatic$dp.test <- somatic$dp >= min.sc.dp & somatic$bulk.dp >= 
        min.bulk.dp
    cat("step 6: calling somatic SNVs\n")
    somatic$hard.filter <- somatic$abc.pv > 0.05 & somatic$cigar.id.test & 
        somatic$cigar.hs.test & somatic$lowmq.test & somatic$dp.test & 
        somatic[, scalt] >= min.sc.alt
    somatic$pass <- somatic$hard.filter & somatic$lysis.fdr <= 
        target.fdr & somatic$mda.fdr <= target.fdr
    cat(sprintf("        %d passing somatic SNVs\n", sum(somatic$pass)))
    cat(sprintf("        %d filtered somatic SNVs\n", sum(!somatic$pass)))
    return(c(call.fingerprint, list(somatic = somatic)))
}

get.3mer <- function (df) 
{
    require(BSgenome)
    require(BSgenome.Hsapiens.1000genomes.hs37d5)
    comp <- c("A", "C", "G", "T")
    names(comp) <- c("T", "G", "C", "A")
    x <- df
    if (!("muttype" %in% colnames(x))) {
        cat("adding mutation types..\n")
        x$muttype <- muttype.map[paste(x$refnt, x$altnt, sep = ">")]
    }
    x$ctx <- getSeq(BSgenome.Hsapiens.1000genomes.hs37d5, names = x$chr, 
        start = x$pos - 1, end = x$pos + 1, as.character = TRUE)
    x$ctx.rc <- sapply(strsplit(x$ctx, ""), function(s) paste0(comp[s[c(3, 
        2, 1)]], collapse = ""))
    x$type.and.ctx <- ifelse(x$refnt == "C" | x$refnt == "T", 
        paste0(x$ctx, ":", x$muttype), paste0(x$ctx.rc, ":", 
            x$muttype))
    x
}

get.callable <- function (ss.dir, verbose = TRUE) 
{
    require(yaml)
    ss.config <- file.path(ss.dir, "config.yaml")
    if (!file.exists(ss.config)) 
        stop(sprintf("expected SCAN-SNV config file does not exist: %s\n", 
            ss.config))
    yaml <- read_yaml(ss.config)
    sc.samples <- yaml$sc_samples
    sapply(sc.samples, function(sn) {
        f <- sprintf("%s/callable_regions/%s/callable_regions.bed", 
            ss.dir, sn)
        if (verbose) 
            print(f)
        bed <- read.table(f, sep = "\t", header = F)
        sum(as.numeric(bed[, 3] - bed[, 2]))
    })
}

get.distance.distn <- function (d, min = 1, max = 5) 
{
    h <- hist(d[d >= min & d <= max], breaks = 50, plot = FALSE)
    h$density <- h$density/sum(h$density)
    h
}

get.fdr.tuning.parameters <- function (somatic, hsnps, bins = 20, random.seed = 0) 
{
    cat(sprintf("estimating bounds on somatic mutation rate (seed=%d)..\n", 
        random.seed))
    set.seed(random.seed)
    max.dp <- as.integer(quantile(hsnps$dp, prob = 0.8))
    fcs <- lapply(0:max.dp, function(dp) fcontrol(germ.df = hsnps[hsnps$dp == 
        dp, ], som.df = somatic[somatic$dp == dp, ], bins = bins))
    fc.max <- fcontrol(germ.df = hsnps[hsnps$dp > max.dp, ], 
        som.df = somatic[somatic$dp > max.dp, ], bins = bins)
    fcs <- c(fcs, list(fc.max))
    cat(sprintf("        profiled hSNP and somatic VAFs at depths %d .. %d\n", 
        0, max.dp))
    burden <- as.integer(c(sum(sapply(fcs, function(fc) fc$est.somatic.burden[1])), 
        sum(sapply(fcs, function(fc) fc$est.somatic.burden[2]))))
    cat(sprintf("        estimated callable somatic mutation burden range (%d, %d)\n", 
        burden[1], burden[2]))
    cat("          -> using MAXIMUM burden\n")
    list(bins = bins, burden = burden, fcs = fcs, max.dp = max.dp)
}

get.filter.reasons <- function (df, fdr.threshold = 0.01, min.alt = 2) 
{
    alt.reads <- round(df$af * df$dp)
    m <- cbind(pass = df$pass, abc.test = df$abc.pv <= 0.05, 
        lysis.test = df$lysis.fdr > fdr.threshold, mda.test = df$mda.fdr > 
            fdr.threshold, cigar.ID = !df$cigar.id.test, cigar.HS = !df$cigar.hs.test, 
        lowmq = !df$lowmq.test, dp = !df$dp.test, min.alt = is.na(alt.reads) | 
            alt.reads < min.alt)
    apply(m, 1, function(row) paste(colnames(m)[row], collapse = "&"))
}

get.scansnv <- function (ss.dir, type = "somatic", muttype = "snv", verbose = TRUE) 
{
    if (!(type %in% c("somatic", "mosaic", "hsnp_spikein"))) 
        stop(sprintf("type must be either somatic, mosaic or hsnp_spikein, not '%s'", 
            type))
    if (!(muttype %in% c("snv", "indel"))) 
        stop(spritnf("muttype must be either 'snv' or 'indel', not %s", 
            muttype))
    require(yaml)
    ss.config <- file.path(ss.dir, "config.yaml")
    if (!file.exists(ss.config)) 
        stop(sprintf("expected SCAN-SNV config file does not exist: %s\n", 
            ss.config))
    yaml <- read_yaml(ss.config)
    sc.samples <- yaml$sc_samples
    min.sc.alt <- yaml$min_sc_alt
    ret <- lapply(sc.samples, function(s) {
        path.fmt <- "%s_genotypes.rda"
        if (muttype == "indel" & type == "somatic") 
            path.fmt <- "%s_genotypes.pon_filter.rda"
        f <- file.path(ss.dir, muttype, s, sprintf(path.fmt, 
            type))
        if (verbose) 
            print(f)
        load(f)
        somatic <- get(ifelse(type == "hsnp_spikein", "spikeins", 
            type))
        scalt <- which(colnames(somatic) == make.names(s)) + 
            2
        somatic$id <- paste(somatic$chr, somatic$pos, somatic$refnt, 
            somatic$altnt)
        somatic
    })
    names(ret) <- sc.samples
    ret
}

get.sig.score <- function (ssnvs, good.sig, lysis.sig, eps = 0.001) 
{
    require(pracma)
    test.sig <- df.to.sbs96(ssnvs)
    sigs <- cbind(good.sig, lysis.sig)
    weights <- lsqnonneg(sigs, test.sig)$x
    recon <- as.vector(sigs %*% weights)
    weights <- weights + eps
    weights <- weights/sum(weights)
    nsnvs <- sum(ssnvs$filter.reason == "lysis.test")
    postp <- log10(lysis.sig * weights[2]) - log10(good.sig * 
        weights[1])
    list(postp = postp, nsnvs = nsnvs, weight.true = weights[1], 
        weight.artifact = weights[2])
}

infer.gp <- function (ssnvs, fit, hsnps, chunk = 2500, flank = 1e+05, max.hsnps = 150, 
    verbose = FALSE, spikein = FALSE) 
{
    if (verbose) 
        cat(sprintf("mode=%s\n", ifelse(spikein, "spikein", "somatic")))
    if (spikein) {
        if (chunk != 1) 
            stop("infer.gp: can only run in spikein mode with chunk=1\n")
        ssnv.is.hsnp <- ssnvs$pos %in% hsnps$pos
        cat("infer.gp: building ssnvs <-> hsnps map\n")
        hsnp.map <- sapply(1:nrow(ssnvs), function(i) {
            if (!ssnv.is.hsnp[i]) 
                NA
            else which(hsnps$pos == ssnvs$pos[i])
        })
        cat(sprintf("infer.gp: performing %d leave-1-out hSNP AB estimations\n", 
            nrow(ssnvs)))
    }
    nchunks <- ceiling(nrow(ssnvs)/chunk)
    ctx <- abmodel.approx.ctx(c(), c(), c(), hsnp.chunksize = 2 * 
        max.hsnps + 10)
    do.call(rbind, lapply(1:nchunks, function(i) {
        if (i%%100 == 0) 
            cat(sprintf("infer.gp: progress: finished %d of %d sites (%0.1f%%)\n", 
                i, nchunks, 100 * i/nchunks))
        start <- 1 + (i - 1) * chunk
        stop <- min(i * chunk, nrow(ssnvs))
        h <- hsnps
        if (spikein) {
            if (ssnv.is.hsnp[i]) 
                h <- hsnps[-hsnp.map[i], ]
        }
        infer.gp.block(ssnvs[start:stop, , drop = FALSE], fit, 
            h, ctx = ctx, flank = flank, max.hsnps = max.hsnps, 
            verbose = verbose)
    }))
}

infer.gp.block <- function (ssnvs, fit, hsnps, ctx, flank = 1e+05, max.hsnps = 150, 
    verbose = FALSE) 
{
    a <- fit$a
    b <- fit$b
    c <- fit$c
    dparam <- fit$d
    middle <- 0
    right <- findInterval(range(ssnvs$pos)[2], hsnps$pos)
    down <- findInterval(range(ssnvs$pos)[2] + flank, hsnps$pos)
    middle <- down
    down <- min(down, right + max.hsnps)
    left <- findInterval(range(ssnvs$pos)[1], hsnps$pos)
    up <- findInterval(range(ssnvs$pos)[1] - flank, hsnps$pos)
    middle <- middle - up + 1
    up <- max(up, left - max.hsnps)
    window <- c(up, down)
    d <- hsnps[max(window[1], 1):min(window[2], nrow(hsnps)), 
        ]
    if (verbose) {
        print(middle)
        print(window)
        cat(sprintf("infer.gp.block: window=%d-%d, %d nearby hets\n", 
            min(d$pos), max(d$pos), nrow(d)))
        cat(sprintf("positions:"))
        print(ssnvs$pos)
    }
    ctx$x <- d$pos
    ctx$y <- d$hap1
    ctx$d <- d$hap1 + d$hap2
    abmodel.approx.logp(a = a, b = b, c = c, d = dparam, ctx = ctx)
    z2 <- alg3.2.2(a = a, b = b, c = c, d = dparam, ctx = ctx, 
        Xnew = ssnvs$pos)
    data.frame(gp.mu = z2$mean, gp.sd = sqrt(diag(z2$cov)))
}

match.fdr3 <- function (pv, gp.mu, gp.sd, dp, nt, na, div = 2) 
{
    alphabeta <- estimate.alphabeta3(gp.mu = gp.mu, gp.sd = gp.sd, 
        dp = dp, div = div)
    x <- data.frame(alpha = alphabeta$alphas, beta = alphabeta$betas, 
        fdr = ifelse(alphabeta$alphas * na + alphabeta$betas * 
            nt > 0, alphabeta$alphas * na/(alphabeta$alphas * 
            na + alphabeta$betas * nt), 0))
    x <- rbind(c(1, 0, 1), x)
    x <- x[pv <= x$alpha, ]
    x <- x[x$fdr == min(x$fdr), , drop = F]
    unlist(x[1, ])
}

plot.3mer <- function (x, no.legend = FALSE, ...) 
{
    bases <- c("A", "C", "G", "T")
    t <- rep(0, 96)
    names(t) <- paste0(rep(bases, each = 4), rep(c("C", "T"), 
        each = 48), rep(bases, times = 4), ":", rep(c("C", "T"), 
        each = 48), ">", c(rep(c("A", "G", "T"), each = 16), 
        rep(c("A", "C", "G"), each = 16)))
    t2 <- table(x$type.and.ctx)
    t[names(t2)] <- t2
    tn <- do.call(rbind, strsplit(names(t), ":"))
    t <- t[order(tn[, 2])]
    p <- barplot(t, las = 3, col = mutsig.cols, names.arg = tn[order(tn[, 
        2]), 1], space = 0.5, border = NA, ...)
    abline(v = (p[seq(4, length(p) - 1, 4)] + p[seq(5, length(p), 
        4)])/2, col = "grey")
    if (!no.legend) 
        legend("topright", ncol = 2, legend = sort(unique(tn[, 
            2])), fill = mutsig.cols[seq(1, length(mutsig.cols), 
            16)])
    invisible(t)
}

plot.ab <- function (ab) 
{
    layout(matrix(1:4, nrow = 2, byrow = T))
    td <- ab$td[order(ab$td$dp), ]
    plot(td$dp, td$mut, type = "l", ylim = range(td[, c("mut", 
        "err1", "err2")]), xlab = "Depth", ylab = "Model probabiblity")
    lines(td$dp, pmax(td$err1, td$err2), lty = "dotted", col = 2)
    plot(x = ab$input.alphas, y = ab$alphas, xlab = "Requested alpha", 
        ylab = "Estimated alpha", log = "xy")
    plot(ab$alphas, ab$betas, xlab = "FP rate", ylab = "Power")
    abline(h = 0, lty = "dotted")
    plot(ab$alphas, ab$betas, log = "x", xlab = "log(FP rate)", 
        ylab = "Power")
    abline(h = 0, lty = "dotted")
}

plot.fcontrol <- function (fc) 
{
    layout(matrix(1:(1 + length(fc$pops)), nrow = 1))
    plot(fc$binmids, fc$g/sum(fc$g), ylim = range(c(fc$g/sum(fc$g), 
        fc$s/sum(fc$s))), type = "l", lwd = 2)
    lines(fc$binmids, fc$s/sum(fc$s), col = 2, lwd = 2)
    for (i in 1:length(fc$pops)) {
        pop <- fc$pops[[i]]
        barplot(names.arg = fc$binmids, t(pop), col = 1:2, main = sprintf("Assumption: ~%d true sSNVs", 
            sum(round(pop[, 1], 0))), las = 2)
        legend("topright", fill = 1:2, legend = c("Ntrue", "Nartifact"))
    }
}

plot.fdr <- function (fc, dps = c(10, 20, 30, 60, 100, 200), target.fdr = 0.1, 
    div = 2) 
{
    afs <- fc$binmids
    layout(matrix(1:(3 * length(fc$pops)), nrow = 3))
    for (i in 1:length(fc$pops)) {
        pop <- fc$pops[[i]]
        l <- lapply(dps, function(dp) {
            sapply(1:length(afs), function(i) match.fdr(afs[i], 
                dp, nt = pop[i, 1], na = pop[i, 2], target.fdr = target.fdr, 
                div = div))
        })
        matplot(x = afs, sapply(l, function(ll) ll[3, ]), type = "l", 
            lty = 1, main = sprintf("Assuming %d true sSNVs", 
                sum(pop[, 1])), xlab = "AF (binned)", ylab = "FDR", 
            ylim = c(0, 1.1 * target.fdr))
        abline(h = target.fdr, lty = "dotted")
        matplot(x = afs, sapply(l, function(ll) ll[1, ]), type = "l", 
            lty = 1, xlab = "AF (binned)", ylab = "log(alpha)", 
            log = "y", ylim = c(1e-05, 1))
        abline(h = 10^-(1:5), lty = 2)
        matplot(x = afs, sapply(l, function(ll) ll[2, ]), type = "l", 
            lty = 1, xlab = "AF (binned)", ylab = "Power", ylim = 0:1)
    }
}

plot.gp.confidence <- function (pos, gp.mu, gp.sd, df, sd.mult = 2, logspace = FALSE, 
    tube.col = rgb(0.9, 0.9, 0.9), line.col = "black", tube.lty = "solid", 
    line.lty = "solid", add = TRUE) 
{
    if (!missing(df)) {
        pos <- df$pos
        gp.mu <- df$gp.mu
        gp.sd <- df$gp.sd
    }
    if (!logspace) {
        cat("transforming to AF space...\n")
        sd.upper <- 1/(1 + exp(-(gp.mu + sd.mult * gp.sd)))
        sd.lower <- 1/(1 + exp(-(gp.mu - sd.mult * gp.sd)))
        gp.mu <- 1/(1 + exp(-gp.mu))
    }
    else {
        sd.upper <- gp.mu + sd.mult * gp.sd
        sd.lower <- gp.mu - sd.mult * gp.sd
    }
    if (!add) {
        plot(NA, NA, xlim = range(pos), ylim = 0:1)
    }
    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.lower)), col = tube.col, 
        border = line.col, lty = tube.lty)
    polygon(c(pos, rev(pos)), c(gp.mu, rev(sd.upper)), col = tube.col, 
        border = line.col, lty = tube.lty)
    lines(pos, gp.mu, lwd = 2, col = line.col, lty = line.lty)
}

plot.indel <- function (iclass, tsb = F, proc, reduce.to.id83 = FALSE, xaxt = "n", 
    col, border, make.plot = TRUE, ...) 
{
    iclass.order <- paste(c(rep(1, 24), rep(rep(2:5, each = 6), 
        2), c(2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5)), c(rep(c("Del", 
        "Ins"), each = 12), rep(c("Del", "Ins"), each = 24), 
        rep("Del", 11)), c(rep(rep(c("C", "T"), each = 6), 2), 
        rep("R", 48), rep("M", 11)), c(rep(0:5, 12), c(1, 1, 
        2, 1:3, 1:5)), sep = ":")
    if (tsb) 
        iclass.order <- as.vector(rbind(paste0("T:", iclass.order), 
            paste0("U:", iclass.order)))
    if (reduce.to.id83) 
        iclass <- substr(iclass, 3, 11)
    if (missing(proc)) {
        proc <- sapply(iclass.order, function(ico) sum(iclass == 
            ico))
    }
    else proc <- proc[iclass.order]
    if (missing(border)) 
        border <- rep(c("#FBBD75", "#FC7F24", "#B0DB8E", "#3B9F36", 
            "#FBC9B6", "#F8896D", "#EE453A", "#B91C22", "#C5D4E4", 
            "#8DBAD2", "#4D98C6", "#1D65A8", "#E1E1EE", "#B5B6D6", 
            "#8684BA", "#614398"), c(rep(6, 12), 1, 2, 3, 5))
    if (missing(col)) 
        col <- rep(c("#FBBD75", "#FC7F24", "#B0DB8E", "#3B9F36", 
            "#FBC9B6", "#F8896D", "#EE453A", "#B91C22", "#C5D4E4", 
            "#8DBAD2", "#4D98C6", "#1D65A8", "#E1E1EE", "#B5B6D6", 
            "#8684BA", "#614398"), c(rep(6, 12), 1, 2, 3, 5))
    if (tsb) 
        mutsig.cols <- as.vector(rbind(mutsig.cols, "#D6C2C2"))
    if (make.plot) {
        par(mar = c(4, 4, 1, 1))
        p <- barplot(proc, las = 3, col = col, names.arg = iclass.order, 
            space = 0.5, border = border, cex.names = 0.7, ...)
        abline(v = (p[seq(6, length(p) - 11, 6)] + p[seq(7, length(p) - 
            10, 6)])/2, col = "grey")
        mtext(text = c("del C", "del T", "ins C", "ins T", "del 2", 
            "del 3", "del 4", "del 5+", "ins 2", "ins 3", "ins 4", 
            "ins 5+", "microhom."), side = 1, at = c(mean(p[1:6]), 
            mean(p[7:12]), mean(p[13:18]), mean(p[19:24]), mean(p[25:30]), 
            mean(p[31:36]), mean(p[37:42]), mean(p[43:48]), mean(p[49:54]), 
            mean(p[55:60]), mean(p[61:66]), mean(p[67:72]), mean(p[73:83])))
    }
    invisible(proc)
}

plot.ssnv.region <- function (chr, pos, alt, ref, fits, fit.data, upstream = 50000, 
    downstream = 50000, gp.extend = 1e+05, n.gp.points = 100, 
    blocks = 50) 
{
    d <- fit.data[fit.data$chr == chr & fit.data$pos >= pos - 
        upstream & fit.data$pos <= pos + downstream, ]
    cat("estimating AB in region..\n")
    est.at <- c(seq(pos - upstream, pos - 1, length.out = n.gp.points/2), 
        pos, seq(pos + 1, pos + downstream, length.out = n.gp.points/2))
    fit.chr <- fits[[chr]]
    cat("WARNING: if this function produces NA predictions, try reducing the value of 'blocks'\n")
    gp <- infer.gp(ssnvs = data.frame(pos = est.at), fit = fit.chr, 
        hsnps = fit.data[fit.data$chr == chr, ], chunk = blocks, 
        flank = gp.extend, max.hsnp = 500)
    gp$pos <- est.at
    plot.gp.confidence(df = gp, add = FALSE)
    points(d$pos, d$hap1/d$dp, pch = 20, ylim = 0:1)
    af <- alt/(alt + ref)
    gp.at <- gp[gp$pos == pos, ]$gp.mu
    ab <- 1/(1 + exp(-gp.at))
    af <- ifelse(abs(af - ab) <= abs(af - (1 - ab)), af, 1 - 
        af)
    points(pos, af, pch = 20, cex = 1.25, col = 2)
}

resample.hsnps <- function (som, hsnps, chrom, M = 50) 
{
    tmpsom <- find.nearest.germline(som = som[order(som$pos), 
        ], germ = hsnps, chrs = chrom)
    spos <- log10(abs(tmpsom$pos - tmpsom$nearest.het))
    hsnps$nearest.hsnp <- pmin(diff(c(0, hsnps$pos)), diff(c(hsnps$pos, 
        max(hsnps$pos) + 1)))
    hpos <- log10(hsnps$nearest.hsnp)
    dist.s <- get.distance.distn(spos)
    dist.h <- get.distance.distn(hpos)
    ds <- dist.s$density[findInterval(hpos, dist.s$breaks, all.inside = T)]
    dh <- dist.h$density[findInterval(hpos, dist.h$breaks, all.inside = T)]
    u <- runif(n = length(ds))
    list(selection = data.frame(dist = hsnps$nearest.hsnp, ds = ds, 
        dh = dh, u = u, keep = u < ds/(M * dh)), dist.s = dist.s, 
        dist.h = dist.h)
}

rescore <- function (df, target.fdr, use.pon = FALSE, min.pon.dp = 10, quiet = TRUE) 
{
    newpass <- df$hard.filter & df$lysis.fdr <= target.fdr & 
        df$mda.fdr <= target.fdr
    if (use.pon) 
        newpass <- newpass & (df$dp >= min.pon.dp & (df$unique.donors <= 
            1 | df$max.out <= 2))
    if (!quiet) 
        cat(sprintf("rescore: %d passing -> %d passing\n", sum(df$pass), 
            sum(newpass)))
    df$pass <- newpass
    df
}

rescue <- function (df, lysis.sig, good.sig, rescue.fdr, ...) 
{
    sigscores <- get.sig.score(ssnvs = df[df$filter.reason == 
        "lysis.test", ], lysis.sig = lysis.sig, good.sig = good.sig, 
        ...)
    postp <- sigscores$postp
    df$rweight <- 10^-postp[df$type.and.ctx]
    df$lysis.fdr2 <- df$lysis.alpha/(df$lysis.alpha + df$lysis.beta * 
        df$rweight * df$nt/df$na)
    df$pass2 <- !df$pass & df$lysis.fdr2 <= rescue.fdr & df$filter.reasons == 
        "lysis.test"
    df$filter.reasons[df$pass2] <- paste0(df$filter.reasons[df$pass2], 
        ";rescue")
    list(df = df, postp = postp, nsnvs = sigscores$nsnvs, weight.true = sigscores$weight.true, 
        weight.artifact = sigscores$weight.artifact)
}

scansnv.df.to.vcf <- function (df, out.file, ss.config, yaml, sample.name, overwrite = FALSE, 
    chrs = c(1:22, "X", "Y")) 
{
    if (!missing(yaml) & !missing(ss.config)) 
        stop("only one of \"ss.config\" or \"yaml\" can be specified")
    if (!missing(ss.config)) 
        yaml <- read_yaml(ss.config)
    if (!missing(yaml) | !missing(ss.config)) {
        ref.genome <- yaml$humref
        if (missing(chrs)) 
            chrs <- yaml$chrs
        fai <- read.table(paste0(ref.genome, ".fai"), sep = "\t", 
            stringsAsFactors = F)
    }
    if (!overwrite & file.exists(out.file)) 
        stop(sprintf("output file %s already exists, please delete it first", 
            out.file))
    f <- file(out.file, "w")
    cat(sprintf("writing to %s..\n", out.file))
    vcf.header <- c("##fileformat=VCFv4.0", "##source=scansnv", 
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    if (!missing(ss.config) | !missing(yaml)) 
        vcf.header <- c(vcf.header, sprintf("##reference=%s", 
            yaml$humref), sprintf("##contig=<ID=%s,length=%d>", 
            fai[, 1], fai[, 2]))
    vcf.header <- c(vcf.header, paste(c("#CHROM", "POS", "ID", 
        "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", sample.name), 
        collapse = "\t"))
    writeLines(vcf.header, con = f)
    s <- df[!is.na(df$pass) & !is.na(df$chr) & df$pass, ]
    s <- do.call(rbind, lapply(chrs, function(chr) {
        ss <- s[s$chr == chr, ]
        ss[order(ss$pos), ]
    }))
    writeLines(paste(s$chr, s$pos, s$dbsnp, s$refnt, s$altnt, 
        ".", "PASS", ".", "GT", "0/1", sep = "\t"), con = f)
    close(f)
}

scansnv.to.vcf <- function (ss.dir, output.fmt, type = "somatic", muttype = "snv", 
    overwrite = FALSE) 
{
    if (!(type %in% c("somatic", "mosaic", "hsnp_spikein"))) 
        stop(sprintf("type must be either somatic, mosaic or hsnp_spikein, not '%s'", 
            type))
    if (!(muttype %in% c("snv", "indel"))) 
        stop(spritnf("muttype must be either 'snv' or 'indel', not %s", 
            muttype))
    require(yaml)
    ss.config <- file.path(ss.dir, "config.yaml")
    if (!file.exists(ss.config)) 
        stop(sprintf("expected SCAN-SNV config file does not exist: %s\n", 
            ss.config))
    yaml <- read_yaml(ss.config)
    sc.samples <- yaml$sc_samples
    for (s in sc.samples) {
        path.fmt <- "%s_genotypes.rda"
        if (muttype == "indel" & type == "somatic") 
            path.fmt <- "%s_genotypes.pon_filter.rda"
        f <- file.path(ss.dir, muttype, s, sprintf(path.fmt, 
            type))
        print(f)
        load(f)
        somatic <- get(ifelse(type == "hsnp_spikein", "spikeins", 
            type))
        out.file <- sprintf(output.fmt, s)
        scansnv.df.to.vcf(df = somatic, out.file = out.file, 
            yaml = yaml, sample.name = s, overwrite = overwrite)
    }
}

test2 <- function (altreads, gp.mu, gp.sd, dp, div) 
{
    err <- (dreads(0:dp, d = dp, gp.mu = gp.mu, gp.sd = gp.sd, 
        factor = div) + dreads(0:dp, d = dp, gp.mu = -gp.mu, 
        gp.sd = gp.sd, factor = div))/2
    p.value <- sum(err[err <= err[altreads + 1]])
    c(p.value = p.value)
}

