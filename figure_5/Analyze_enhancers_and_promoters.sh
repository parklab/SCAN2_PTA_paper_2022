#!/bin/bash

for muttype in sSNVs sIndels; do 
    perm_prefix=
    if [ $muttype == 'sIndels' ]; then
        perm_prefix='indel_'
    fi
    for celltype in astrocyte microglia neuron oligo; do
        for rtype in enhancers promoters; do
            ./bedenrich.R 10000 SCAN2_PTA_${muttype}_filtered.rda \
                ./${perm_prefix}perms_10000_mutsigaware_grl.minimized.rda \
                ${muttype}_vs_${celltype}_${rtype}.FULL.rda \
                ${muttype}_vs_${celltype}_${rtype}.SUMMARY.rda \
                ../external_data/Enhancers_and_promoters/s5_${celltype}_${rtype}.uniq.txt
        done
    done
done
