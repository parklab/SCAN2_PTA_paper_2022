#!/bin/bash

nbootstraps=10000

./qbedenrich.R $nbootstraps \
    ./SCAN2_PTA_sSNVs_filtered.rda \
    ./perms_10000_mutsigaware_grl.minimized.rda \
    Enrichment_results/sSNVs_vs_Conservation.FULL.rda \
    Enrichment_results/sSNVs_vs_Conservation.SUMMARY.rda \
    ../external_data/Conservation/hg19.100way.phyloP100way.qbed


./qbedenrich.R $nbootstraps \
    ./SCAN2_PTA_sIndels_filtered.rda \
    ./indel_perms_10000_mutsigaware_grl.minimized.rda \
    Enrichment_results/sIndels_vs_Conservation.FULL.rda \
    Enrichment_results/sIndels_vs_Conservation.SUMMARY.rda \
    ../external_data/Conservation/hg19.100way.phyloP100way.qbed
