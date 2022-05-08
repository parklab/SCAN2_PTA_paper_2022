#!/bin/bash


for muttype in sSNVs sIndels; do 
    perm_prefix=
    if [ $muttype == 'sIndels' ]; then
        perm_prefix='indel_'
    fi
    cat ../external_data/Roadmap_Epigenomics/EIDs_with_H3K27ac.txt | while read eid; do
        ./bedenrich.R 10000 SCAN2_PTA_${muttype}_filtered.rda \
            ./${perm_prefix}perms_10000_mutsigaware_grl.minimized.rda \
            ${muttype}_vs_${eid}_ChromHMM15.FULL.rda \
            ${muttype}_vs_${eid}_ChromHMM15.SUMMARY.rda \
            ../external_data/Roadmap_Epigenomics/ChromHMM_15state_model/beds/${eid}.bed
    done
done
