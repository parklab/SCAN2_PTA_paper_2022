#!/bin/bash

for celltype in GLU GABA OLIG MGAS; do
    wget https://s3.amazonaws.com/ggoma/${celltype}_DLPFC.bb
    bigBedToBed ${celltype}_DLPFC.bb ${celltype}_DLPFC.bed
done
