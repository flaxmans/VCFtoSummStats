#!/bin/bash

mydate="8.06"

printf "\n\nOld version:\n"
./ConductSeveralTestRuns.sh 0 oldVCFprog &
wait
fname="foo.master.${mydate}.out"
mv foo.out $fname
fnames="$fname"

## 10k
printf "\n\n10k:\n"
./ConductSeveralTestRuns.sh 10000 VCFtoSummStats &
wait
fname="foo.10kbuff.${mydate}.out"
mv foo.out $fname
fnames="$fnames $fname"

## 100k
printf "\n\n100k:\n"
./ConductSeveralTestRuns.sh 100000 VCFtoSummStats &
wait
fname="foo.100kbuff.${mydate}.out"
mv foo.out $fname
fnames="$fnames $fname"

## 1M
printf "\n\n1M:\n"
./ConductSeveralTestRuns.sh 1000000 VCFtoSummStats &
wait
fname="foo.1Mbuff.${mydate}.out"
mv foo.out $fname
fnames="$fnames $fname"

## 10M
printf "\n\n10M:\n"
./ConductSeveralTestRuns.sh 10000000 VCFtoSummStats &
wait
fname="foo.10Mbuff.${mydate}.out"
mv foo.out $fname
fnames="$fnames $fname"

## 100M
printf "\n\n100M:\n"
./ConductSeveralTestRuns.sh 100000000 VCFtoSummStats &
wait
fname="foo.100Mbuff.${mydate}.out"
mv foo.out $fname
fnames="$fnames $fname"


printf "\n\nOutput files are:\n\t$fnames\n\n" >> OutputFileNames.txt

