#!/bin/bash

if [ $# -ne 6 ]; then
    echo "usage: $0 human_reference.fasta GATK_input.vcf bulk_sample_ID out.vcf out.tab [snv|indel]"
    exit 1
fi

humref=$1
gatkvcf=$2
bulksample=$3
outvcf=$4
outtab=$5
type=$6

if [ "x$type" == "xsnv" ]; then
    gatktype="SNP"
    script="totab.sh"
elif [ "x$type" != "xindel" ]; then
    gatktype="INDEL"
    script="totab.indel.sh"
else
    echo "type=snv or indel (case sensitive) is required"
    exit 1
fi

gatk3 -Xmx3G -Xms3G \
           -T SelectVariants \
           -R $humref \
           -V $gatkvcf \
           -selectType $gatktype -restrictAllelesTo BIALLELIC \
           -env -trimAlternates \
           -select 'vc.getGenotype("'"$bulksample"'").isCalled()' \
           -o $outvcf

./$script $outvcf $outtab
