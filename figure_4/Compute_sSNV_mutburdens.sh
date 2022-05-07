#!/bin/bash

# scan2 mutburden
./mutburden.R \
    --metadata ../metadata.csv \
    ../data/Collected_SCANSNV_SCAN2_sSNV_calls.csv.gz \
    ../data/PROTECTED_hSNPs_FDR_recalculated.csv.gz \
    ../data/Collected_SCAN2_sequencing_depth_tables.rda \
    Collected_SCAN2_sSNV_burdens.csv
