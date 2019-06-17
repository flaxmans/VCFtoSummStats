#!/bin/bash

cd ExampleDataFiles/

testList="Small_hmel2.5.30f4.vcf Medium_hmel2.5.30f4.vcf OneHundredK_TotChromCloeNolan.vcf BASW_100k.vcf Total_Chromosomes_Combined_with_header.vcf hmel2.5.30f4.vcf RADseq_allSamples.vcf"

endings1="_Unfiltered_Summary.tsv .bz2_Unfiltered_Summary.tsv .gz_Unfiltered_Summary.tsv"

endings2="_discardedLineNums.txt .bz2_discardedLineNums.txt .gz_discardedLineNums.txt"

for i in $testList
do
	ref="ref_${i}_Unfiltered_Summary.tsv"
	for j in $endings1
	do
		testFile="${i}${j}"
		if [ -f "$testFile" ]
		then
			cmd="diff $ref $testFile"
			echo "$cmd :"
			$cmd | wc -l
			printf "\n\t*******************\n"
		fi
	done
	ref="ref_${i}_discardedLineNums.txt"
	for j in $endings2
	do
		testFile="${i}${j}"
		if [ -f "$testFile" ]
		then
			cmd="diff $ref $testFile"
			echo "$cmd :"
			$cmd | wc -l
			printf "\n\t*******************\n"
		fi
	done
done


cd ..
