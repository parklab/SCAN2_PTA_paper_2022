#!/bin/bash

wget -O GSE167257_SARseq_iNeuron_OverlapRep123.peaks.bed.gz \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE167257&format=file&file=GSE167257%5FSARseq%5FiNeuron%5FOverlapRep123%2Epeaks%2Ebed%2Egz"

gunzip GSE167257_SARseq_iNeuron_OverlapRep123.peaks.bed.gz \
