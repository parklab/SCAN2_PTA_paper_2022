#!/usr/bin/env Rscript

load("/n/data1/hms/dbmi/park/jluquette/pta/integrated_calls/metadata.rda")

library(data.table)

abs <- lapply(meta$sample, function(s) {
    donor <- meta[s,]$donor
    if (meta[s,]$amp == 'bulk') {
        f <- paste0('/home/ljl11/ndata1/pta/', donor, '/for_bulk_allelebal.tab')
        print(f)
        data <- fread(f)
    } else {
        f <- paste0('/home/ljl11/ndata1/pta/', donor, '/scansnv_fdr01_noX/',
            '/ab_model/', s, '/training.rda') # loads 'data'
        print(f)
        load(f)
    }
    # just return AFs to save memory
    data$hap1 / data$dp
})
names(abs) <- meta$sample
print(gc())

# ~3GB of data, plus perhaps somehow contains identifiable information 
abs.full <- abs  

# By sorting the VAFs, we not only destroy possibly identifiable information
# (not that I can imagine how), but also allow significant compression of
# data by run length encoding (rle).  About 100-fold size reduction. Rounding
# to 3 decimal places is another 10x size reduction. 3GB -> ~2MB.
abs <- lapply(abs.full, function(ab) rle(sort(round(ab,3))))
save(abs, file='Collected_allele_balance_at_training_hSNPs.rda')
