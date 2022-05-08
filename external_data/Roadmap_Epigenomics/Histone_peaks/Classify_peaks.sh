#!/bin/bash

gfile=hg19.1-22XYM.genome

echo "eid k4me3 k27ac prox dist" | tr ' ' '\t'
cat ../EIDs_with_H3K27ac.txt | while read eid; do
    k27=H3K27ac/${eid}.narrowPeak
    k4=H3K4me3/${eid}.narrowPeak
    k4slop="Distal_proximal_peaks/${eid}-H3K4me3.slop2000.bed"
    k27tss="Distal_proximal_peaks/${eid}-H3K27ac.TSS.bed"
    k27distal="Distal_proximal_peaks/${eid}-H3K27ac.distal.bed"
    k27peaks=$(wc -l < $k27 | tr -d '\t\n\r ')
    bedtools slop -b 2000 -i $k4 -g $gfile > ${k4slop}
    k4peaks=$(wc -l < $k4slop | tr -d '\t\n\r ')
    bedtools intersect -u -a $k27 -b $k4slop > $k27tss
    proximal=$(wc -l < $k27tss | tr -d '\t\n\r ')
    bedtools intersect -v -a $k27 -b $k4slop > $k27distal
    distal=$(wc -l < $k27distal | tr -d '\t\n\r ')
    k27classified=Distal_proximal_peaks/${eid}.prox_dist.bed
    awk 'BEGIN { OFS="\t"; } { print $1, $2, $3, "TSS_proximal"; }' $k27tss \
        | bedtools sort > tmp1.bed
    awk 'BEGIN { OFS="\t"; } { print $1, $2, $3, "distal"; }' $k27distal \
        | bedtools sort > tmp2.bed
    cat tmp1.bed tmp2.bed | bedtools sort > $k27classified
    echo "$eid $k4peaks $k27peaks $proximal $distal" | tr ' ' '\t'
done
