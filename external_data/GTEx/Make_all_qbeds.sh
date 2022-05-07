#!/bin/bash

cat GTEx_tissues.txt | while read line; do
    echo $line "-------------------------------------------------"
    filepre=$(echo $line|tr '() ' '_')
    ./Convert_GTEx_expression_to_BigWig.R \
        gencode.v26lift37.annotation.GTEX_COLLAPSED.genes_only.gtf \
        "$line" \
        GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz \
        bigwigs/$filepre.bw
    ./Make_qbed_from_BigWig.sh 10 0.8 \
        ../../data/Alignable_genome_tiles_1kb.bed.gz \
        bigwigs/$filepre.bw qbeds/$filepre.qbed \
        tmp1 tmp2 tmp3 \
        BINSIZE=1000 datasource=gtex tissue="$line"
    rm tmp1 tmp2 tmp3
done
