#!/bin/bash


for muttype in sSNVs sIndels; do 
    perm_prefix=
    if [ $muttype == 'sIndels' ]; then
        perm_prefix='indel_'
    fi
    ./bedenrich.R 10000 SCAN2_PTA_${muttype}_filtered.rda \
        ./${perm_prefix}perms_10000_mutsigaware_grl.minimized.rda \
        Enrichment_results/${muttype}_vs_Repairseq.FULL.rda \
        Enrichment_results/${muttype}_vs_Repairseq.SUMMARY.rda \
        ../external_data/DNA_repair_hotspots/2021_Science_Reid_Gage_Table_S1.txt
done

for muttype in sSNVs sIndels; do 
    perm_prefix=
    if [ $muttype == 'sIndels' ]; then
        perm_prefix='indel_'
    fi
    ./bedenrich.R 10000 SCAN2_PTA_${muttype}_filtered.rda \
        ./${perm_prefix}perms_10000_mutsigaware_grl.minimized.rda \
        Enrichment_results/${muttype}_vs_SARseq.FULL.rda \
        Enrichment_results/${muttype}_vs_SARseq.SUMMARY.rda \
        ../external_data/DNA_repair_hotspots/GSE167257_SARseq_iNeuron_OverlapRep123.peaks.bed
done

