#!/bin/bash

# script to invoke VCFtoSummStats program with useful inputs
# this script expects two arguments: VCF file name as arg1 and
# Population file as arg2.

TARGET=VCFtoSummStats;

# check for presence of needed files:
errMsg="\nPlease check name spelling and path, and include extension in filename.\nThe call to this script should be:\n\n\tsh WrapperForVCFtoSummStats.sh VCFfileName PopulationFileName\n\nExiting!\n"
if [ ! -f "$1" ]
then
	echo "\nError! file '$1' not found!"
	echo "$errMsg"
	exit -1
fi
if [ ! -f "$2" ]
then
	echo "Error! file '$2' not found!  Please check name spelling and path!"
	echo "$errMsg"
	exit -1
fi
# delete the following if the executable is installed in the PATH
if [ ! -x "${TARGET}" ]
then
	echo "Executable ${TARGET} not found.  Calling make ... "
	make
	wait
fi

# some info on populations and samples:
# population file: does it have a header row?
printf "\nDoes the population file $2 have a header row? (Y/n):  "
read headers
echo " "
if [ "$headers" == "y" ] || [ "$headers" == "Y" ]
then
	header=1
	tail +2 $2 > fooDataTmp
	usefile=fooDataTmp
elif [ "$headers" == "n" ] || [ "$headers" == "N" ]
then
	header=0
	usefile=$2
else
	echo "Input not recognized!!  Exiting ... "
	exit -1
fi
# calculate number of samples listed in population file
numSamples1=$(wc -l $usefile | awk '{print $1}')
numSamples2=$(cut -f1 $usefile | sort | uniq | wc -l | awk '{print $1}')
if [ $numSamples1 -eq $numSamples2 ]
then
	echo "I found $numSamples1 samples in the population file $2."
else
	echo "Warning!  Samples may be repeated in population file!"
	echo "Number of data lines = $numSamples1, but"
	echo "Number of unique sample IDs (column 1) = $numSamples2"
    echo "Number of unique samples will be used because VCF format"
    echo "mandates that all sample IDs are unique."
fi
# calculate number of populations present in population file
numPopulations=$(cut -f2 $usefile | sort | uniq | wc -l | awk '{print $1}')
echo "I found $numPopulations populations in the population file $2."
if [ -f "fooDataTmp" ]
then
	rm fooDataTmp
fi

# some info on VCF file:
# calculate number of fields of data:
# this should get the header row that follows the initial rows of meta-data
numFields=$(grep -v "##" $1 | head -n 1 | awk '{print NF}')
echo "\nI found $numFields fields of data in VCF file $1\n"

# do some checking against expected format:
echo "\nVCF file format specification means that the VCF file should have"
echo "a total number of fields equal to 9 + the number of distinct samples."
echo "The first nine expected fields are:"
echo "\t#CHROM  POS  ID  REF ALT  QUAL  FILTER  INFO  FORMAT"
echo "This is based upon:\n\thttp://samtools.github.io/hts-specs/VCFv4.3.pdf\n\t(accessed 5/31/19)".
echo "If your VCF differs from these expectations, then the program will NOT work."

expectedCols=$((9+$numSamples2))
if [ ! $expectedCols -eq $numFields ]
then
    echo "** Error! **\n\tNumber of fields found in $1 "
    echo "\tdoes not match expected number of fields based upon"
    echo "\t9 + number of samples found in $2."
    echo "\t ** aborting execution ** "
    exit -1
else
    echo "\n\tNumber of fields found in $1 "
    echo "\tmatches expected number of fields based upon"
    echo "\t9 + number of samples found in $2."
fi

# count maximum line length:
maxCharPerLine=$(awk '{ print length }' $1 | sort -nr | head -n 1)

# call program:
printf "\nNow invoking VCFtoSummStats with following command:\n\t"
myCmd="./VCFtoSummStats -V $1 -P $2 -L $maxCharPerLine -H $header -N $numSamples2 -n $numPopulations -F $numFields"
echo "$myCmd"
$myCmd
