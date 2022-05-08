#!/bin/bash

cat ../EIDs_with_H3K27ac.txt | while read id; do
    wget -O beds/${id}.bed.gz \
        "https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/${id}_15_coreMarks_mnemonics.bed.gz"
done

echo "unzipping peak files.."
gunzip beds/*.gz
