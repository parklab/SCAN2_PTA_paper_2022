#!/bin/bash

for muttype in sSNVs sIndels; do 
    perm_prefix=
    if [ $muttype == 'sIndels' ]; then
        perm_prefix='indel_'
    fi
    for celltype in GLU GABA OLIG MGAS; do
        ./bedenrich.R 10000 SCAN2_PTA_${muttype}_filtered.rda \
            ./${perm_prefix}perms_10000_mutsigaware_grl.minimized.rda \
            Enrichment_results/${muttype}_vs_${celltype}_ATACseq.FULL.rda \
            Enrichment_results/${muttype}_vs_${celltype}_ATACseq.SUMMARY.rda \
            ../external_data/ATACseq/${celltype}_DLPFC.bed
    done
done
