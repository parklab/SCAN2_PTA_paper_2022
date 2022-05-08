#!/bin/bash

nbootstraps=10000

./qbedenrich.R $nbootstraps \
    ./SCAN2_PTA_sSNVs_filtered.rda \
    ./perms_10000_mutsigaware_grl.minimized.rda \
    sSNVs_vs_GTEx_all_tissues.FULL.rda \
    sSNVs_vs_GTEx_all_tissues.SUMMARY.rda \
    ../external_data/GTEx/qbeds/*.qbed


./qbedenrich.R $nbootstraps \
    ./SCAN2_PTA_sIndels_filtered.rda \
    ./indel_perms_10000_mutsigaware_grl.minimized.rda \
    sIndels_vs_GTEx_all_tissues.FULL.rda \
    sIndels_vs_GTEx_all_tissues.SUMMARY.rda \
    ../external_data/GTEx/qbeds/*.qbed
