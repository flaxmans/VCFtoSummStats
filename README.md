# VCFtoSummStats

Tools for generating a plain-text file of summary statistics from a VCF

## Building the program

**Building the program requires the g++ compiler.**

You can build (compile) the program from the UNIX/Linux-style command line by
1. cloning this repo
2. cd'ing to it on your command line
3. typing "make" at the command line and pressing return/enter


## Running the program

Use of the helper wrapper script is strongly recommended.  For example, a test case can be invoked in the UNIX/Linux-like terminal like this:

`sh WrapperForVCFtoSummStats.sh VCFtenLines.txt popMap.txt`

where `VCFtenLines.txt` is a small test snippet of a VCF file containing one header row and 9 SNP rows, and `popMap.txt` is a two-column, plain-text datafile, in which the first column is the sample IDs (corresponding to columns ID's of sample data in the VCF file) and the second column is a population designation.


