#!/bin/bash

# scan2 mutburden
./mutburden.R \
    --metadata ../metadata.csv \
    --min-sc-dp 10 \
    ../data/Collected_SCAN2_sIndel_calls.csv.gz \
    ../data/PROTECTED_hIndels_FDR_recalculated.csv.gz \
    ../data/Collected_SCAN2_sequencing_depth_tables.rda \
    Collected_SCAN2_sIndel_burdens.csv
