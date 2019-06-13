#!/bin/bash

numTests=6

fileList=(BASW_100k.vcf.gz Medium_hmel2.5.30f4.vcf.gz OneHundredK_TotChromCloeNolan.vcf.gz RADseq_allSamples.vcf.gz hmel2.5.30f4.vcf.gz)

# following order must correspond to fileList VCFs!
popFiles=(BASW_popMap.txt popFileHmel.txt CloeNolanPopData.txt BASW_popMap.txt popFileHmel.txt)

for ((i=1; i < 3; i++))
do
    vf=${fileList[$i]}
    pf=${popFiles[$i]}
    cmd="./VCFtoSummStats -V ExampleDataFiles/${vf} -P ExampleDataFiles/${pf}"
    echo "$cmd"
    $cmd
    printf "\nDone with test $i\n\n************\n"
done
