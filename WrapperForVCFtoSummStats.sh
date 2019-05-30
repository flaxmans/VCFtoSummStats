#!/bin/bash

# script to invoke VCFtoSummStats program with useful inputs
# this script expects two arguments: VCF file name as arg1 and
# Population file as arg2.

# check files for presence:
errMsg="call to this script should be:\n\n\tsh WrapperForVCFtoSummStats.sh VCFfileName PopulationFileName\n\nExiting!"
if [ ! -f "$1" ]
then
	echo "Error! file $1 not found!  Please check name spelling and path!"
	echo "$errMsg"
	exit -1
fi
if [ ! -f "$2" ]
then
	echo "Error! file $2 not found!  Please check name spelling and path!"
	echo "$errMsg"
	exit -1
fi

# count line length:
maxCharsPerLine=$(awk '{ print length }' $1 | sort -nr | head -n 1)


# call program:
./VCFtoSummStats -V $1 -P $2 -L $maxCharsPerLine
