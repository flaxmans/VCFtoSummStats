#!/bin/bash

# script to invoke VCFtoSummStats program with useful inputs

# This script REQUIRES two arguments:
    # VCF file name as arg1 and
    # Population file as arg2.
# The script has one optional argument:
    # number of distinct formats in VCF FORMAT column as arg3

## --------- Step 1 DEFAULTS: --------- ##
## Arg 3:  number of different formats in VCF FORMAT column
x="foo"
if [ -z "${3+x}" ]
then
    numFormats=1
else
    numFormats=$3
fi


## ---------- STEP 2: Presence of required files -------------- ##
errMsg="\nPlease check name spelling and path, and include extension in filename.\nThe call to this script should be:\n\n\tsh WrapperForVCFtoSummStats.sh VCFfileName PopulationFileName\n\nExiting!\n"
if [ ! -f "$1" ]
then
    printf "\nError! file '$1' not found!\n"
    printf "$errMsg \n"
    exit -1
fi

if [ ! -f "$2" ]
then
    echo "Error! file '$2' not found!  Please check name spelling and path!"
    printf "$errMsg \n"
    exit -1
fi

# Communicate expected formats:
printf "\nVCF file format specification means that the VCF file's data lines should have\n"
echo "a total number of fields (columns) equal to 9 + the number of distinct samples."
echo "The first nine expected fields are:"
printf "\t#CHROM  POS  ID  REF ALT  QUAL  FILTER  INFO  FORMAT\n"
printf "\tThis is based upon:\n\thttp://samtools.github.io/hts-specs/VCFv4.3.pdf\n\t(accessed 5/31/19)\n"
printf "\tIf your VCF differs from these expectations, then the program will NOT work.\n\n"

printf "The program also expects that your population designation file ($2)\n"
printf "is a tab-separated, plain-text file in which the first column is the unique sample IDs\n"
printf "and the second column designates the population/group/taxon from which the sample\n"
printf "was derived.  By default, this file should have NO headers (see README.md)\n\n"


## --------- STEP 3: Presence of executable -------------- ##

TARGET=VCFtoSummStats;
# delete the following if the executable is installed in the PATH
if [ ! -x "${TARGET}" ]
then
	echo "Executable ${TARGET} not found.  Calling make ... "
	make
	wait
fi

## --------- STEP 4: Interpret file path ( not needed per se ) ----------- ##
# get path to VCF file from input string:
pathToVCF=$(echo $1 | rev | awk -F '/' '{ $1=""; print $0 }' | rev | tr ' ' '/')
pathLength=${#pathToVCF}
if [ $pathLength -eq 0 ]
then
    printf "\nI interpret that the data files are in the same directory as this script.\n"
    pathToVCF="./"
else
    printf "\nI interpret this as the path to the data files: $pathToVCF\n"
fi

## ------------- STEP 5: Parse population designation file ------------- ##
## calculate number of samples listed in population file
#numSamples1=$(awk '{print $1}' $2 | grep -v "^$" | wc -w | awk '{print $1}')
#numSamples2=$(awk '{print $1}' $2 | grep -v "^$" | sort | uniq | wc -l | awk '{print $1}')
#if [ $numSamples1 -eq $numSamples2 ]
#then
#    echo "I found $numSamples1 samples in the population file $2."
#else
#    echo "Error!  Samples may be repeated in population file!"
#    echo "Number of data lines = $numSamples1, but"
#    echo "Number of unique sample IDs (column 1) = $numSamples2"
#    echo "Number of unique samples will be used because VCF format"
#    echo "mandates that all sample IDs are unique."
#    echo "Please check format and accuracy of $2"
#    echo "Aborting .... "
#    exit -1
#fi
## calculate number of populations present in population file
#numPopulations=$(awk '{ print $2 }' $2 | grep -v "^$" | sort | uniq | wc -l | awk '{print $1}')
#echo "I found $numPopulations populations in the population file $2."
## create file of unique population names:
#awk '{ print $2 }' $2 | grep -v "^$" | sort | uniq > ${2}_UniquePopFile.txt
#
## expected number of fields of data in actual data rows of VCF:
#FORMATCOL=9
#numFields=$(($FORMATCOL+$numSamples2)) # 9 is the standard 9 columns listed below

# some info on VCF file:
# calculate number of fields of data:
# this should get the header row that follows the initial rows of meta-data
# numFields=$(grep -v "##" $1 | head -n 1 | awk '{print NF}')
# printf "\nI found $numFields fields of data in VCF file $1\n\n"

# read -n 1 -s -r -p "*** Press any key to continue ***"
# thanks to https://unix.stackexchange.com/questions/293940/bash-how-can-i-make-press-any-key-to-continue
# for read command

#expectedCols=$((9+$numSamples2))
#if [ ! $expectedCols -eq $numFields ]
#then
#    printf "** Error! **\n\tNumber of fields found in $1 \n"
#    printf "\tdoes not match expected number of fields based upon\n"
#    printf "\t9 + number of samples found in $2.\n"
#    printf "\tIf your popFile ($2)\n\thas a header row, call this script"
#    printf "with the third argument as 'Y'\n"
#    printf "\n\t ** Aborting execution ** \n\n"
#    # clean up:
#    if [ -f "fooDataTmp" ]
#    then
#        rm fooDataTmp
#    fi
#    exit -1
#else
#    printf "\n\n\tNumber of fields found in $1 \n"
#    printf "\tmatches expected number of fields based upon\n"
#    printf "\t9 + number of samples found in $2.\n"
#fi


## determine if all genotype formats are same:
#firstLine=$(awk ' { if ( $1 == "#CHROM" )  { print NR; exit }} ' $1)
#firstLine=$(($firstLine+1))
#printf "\nFirst line of actual data in VCF file is line no. $firstLine\n"
## get format field:
#printf "\nChecking for consistency of FORMAT (column no. $FORMATCOL) in VCF file.\n\t... This could take a few minutes ...\n"
#awk -v var="$firstLine" -v col="$FORMATCOL" '{if (NR >= var) print $col}' $1 > foofootmp
#numFormats=$(uniq foofootmp | wc -l | awk '{ print $1 }')
#rm foofootmp
#printf "\tNumber of unique formats in VCF file: $numFormats\n"

# build command for calling program:
printf "\nAbout to invoke VCFtoSummStats with following command:\n\t"
myCmd="./VCFtoSummStats -V $1 -P $2"
echo "$myCmd"
#cmdFileName="${1}_CommandToVCFtoSummStats.txt"
#printf "\tThat command will also be stored for reference in the following file:\n\t${cmdFileName}\n\n"
#echo "$myCmd" > $cmdFileName
printf "\tYou can re-use that command (without re-running this script) if\n\tyou need to re-run the program AND if nothing about your input\n\tfiles has changed.\n\n\tHere we go ...\n\n"

# echo to user:
printf "$myCmd\n\n"

# run program:
$myCmd



# clean up:
if [ -f "fooDataTmp" ]
then
    rm fooDataTmp
fi
