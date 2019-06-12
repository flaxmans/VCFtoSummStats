#!/bin/bash

# filebase name is $1

mycmd="diff ${1}_Unfiltered_Summary.tsv ref_${1}_Unfiltered_Summary.tsv"
echo "$mycmd :"
printf "\n"
$mycmd

mycmd="diff ${1}_discardedLineNums.txt ref_${1}_discardedLineNums.txt"
echo "$mycmd :"
printf "\n"
$mycmd


