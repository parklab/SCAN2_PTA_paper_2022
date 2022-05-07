#!/bin/bash

echo "This file is >1GB in size"
echo "md5sum should be 8e2c404f8a263bc5bc5a2250f3d76970"

wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/GRCh37_mapping/gencode.v26lift37.annotation.gtf.gz
gunzip gencode.v26lift37.annotation.gtf.gz

./GTEx_collapse_annotation.py gencode.v26lift37.annotation.gtf gencode.v26lift37.annotation.GTEX_COLLAPSED.gtf

echo "GTEX_COLLAPSED file should have md5=f28688fcf6774eb360acafecf8df7ad4"
awk '$3 == "gene"' gencode.v26lift37.annotation.GTEX_COLLAPSED.gtf > gencode.v26lift37.annotation.GTEX_COLLAPSED.genes_only.gtf

echo "GTEX_COLLAPSED.genes_only.gtf should have md5=3da0557483fd3a9d20122ad560ff9de5"
