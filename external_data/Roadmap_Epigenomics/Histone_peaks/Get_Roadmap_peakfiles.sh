#!/bin/bash

for mark in H3K27ac H3K4me3; do
    cat ../EIDs_with_H3K27ac.txt | while read id; do
        wget -O $mark/${id}.narrowPeak.gz \
            "https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/${id}-${mark}.narrowPeak.gz"
    done

    echo "unzipping peak files.."
    gunzip $mark/*.gz
done
