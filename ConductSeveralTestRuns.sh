#!/bin/bash

BUFFSIZE=$1
PROGNAME=$2

numTests=6

fileList=(BASW_100k.vcf.gz Medium_hmel2.5.30f4.vcf.gz OneHundredK_TotChromCloeNolan.vcf.gz RADseq_allSamples.vcf.gz hmel2.5.30f4.vcf.gz Total_Chromosomes_Combined_with_header.vcf.gz)

# following order must correspond to fileList VCFs!
popFiles=(BASW_popMap.txt popFileHmel.txt CloeNolanPopData.txt BASW_popMap.txt popFileHmel.txt CloeNolanPopData.txt)

date > foo.out
date > bar.out

for ((i=0; i < 6; i++))
do
    vf=${fileList[$i]}
    if [ -f "ExampleDataFiles/${vf}" ]
    then
      pf=${popFiles[$i]}
      if [ $BUFFSIZE -gt 0 ]
      then
	      cmd="./${PROGNAME} -V ExampleDataFiles/${vf} -P ExampleDataFiles/${pf} -d -1.0 -B $BUFFSIZE"
      else
	      cmd="./${PROGNAME} -V ExampleDataFiles/${vf} -P ExampleDataFiles/${pf} -d -1.0"
      fi
      echo "$cmd"
      echo "$cmd" >> foo.out
      echo "$cmd" >> bar.out
      time $cmd 1>> foo.out 2>> bar.out
      printf "\nDone with test $i\n\n************\n\n" >> foo.out
      printf "\nDone with test $i\n\n************\n\n" >> bar.out
    fi
done
