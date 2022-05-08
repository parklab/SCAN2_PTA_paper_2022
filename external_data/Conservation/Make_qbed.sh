#!/bin/bash

../GTEx/Make_qbed_from_BigWig.sh  10 0 ../../data/Alignable_genome_tiles_1kb.bed.gz hg19.100way.phyloP100way.bw hg19.100way.phyloP100way.qbed tmp1 tmp2 tmp3 BINSIZE=1000 track=phyloP100way

rm tmp1 tmp2 tmp3
