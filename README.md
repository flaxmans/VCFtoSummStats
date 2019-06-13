# VCFtoSummStats

A tool for generating a "tidy", plain-text file of SNP summaries from a VCF.

## Building the program

#### Requirements:
+ **Building the program requires the `g++` compiler and the [Boost C++ libraries](https://www.boost.org/).**  
+ Compiliing and running have only been tested in MacOS and Ubuntu 18.04 
+ See below for detailed advice about downloading the Boost libraries.  
+ MacOS users: you should already have `g++` if you've installed Xcode and the associated command line tools.

#### Compiling:
You can download the source code and then build (compile) the program from the UNIX/Linux-style command line with the following sequence of commands:
```
git clone https://github.com/flaxmans/VCFtoSummStats.git
cd VCFtoSummStats/
make
```


## Running the program

The program is invoked in the UNIX/Linux-like terminal with a command like this:

```
./VCFtoSummStats -V path/to/VCFfile.vcf -P path/to/samplesAndPopulations.txt
```

where `VCFfile.vcf` is your VCF file, and `samplesAndPopulations.txt` is a two-column, 
whitespace-delimited, plain-text datafile, in which the first column is the sample IDs 
(corresponding to columns ID's of sample data in the VCF file) and the second column 
is a population designation for each sample (more details below).
The "path/to/" should be the same for both files.
(See below for an additional example.)

The `-V` and `-P` arguments are **required**,  i.e., you must supply the names of these two files (and paths if they aren't in the current working directory)

Outputs will be written as files in the directory `path/to/`.  
Using the example invocation as given above, the two main outputs that would be generated would be `path/to/VCFfile.vcf_Unfiltered_Summary.tsv` and `path/to/VCFfile.vcf_discardedLineNums.txt`.


## Assumptions about VCF file format and Compression
The program assumes that the VCF file supplied to the program follows the VCF v4.3 format guidelines as found at [http://samtools.github.io/hts-specs/VCFv4.3.pdf](http://samtools.github.io/hts-specs/VCFv4.3.pdf), accessed 5/31/19.

The program can handle compressed VCF files as well, if they have been compressed with either `gzip` (`.gz` extension) or `bzip2` (`.bz2` extension).  *File extensions must be accurate*, i.e., an uncompressed VCF file must have a name ending with `.vcf`, and a compressed VCF file's name must end with the proper extension (`.bz2` or `.gz`).  Note that `gzip` compression is preferred to `bzip2` in the context of shortening run times. 

The population designation file (`-P` argument) must NOT be compressed.

## Requirements of file with information on samples and populations

The second argument (`-P`) to the program as shown in the example command above is the 
name of a file containing information on the samples and the population or taxonomic
group to which each sample belongs.  The sample IDs must be in the first column, and
the population/taxonomic designation must be in the second column. 
The columns must be separated by whitespace. 
The file should **NOT have a header**, but if it does, add a third argument, `-H`, when invoking the program.

## Important default assumptions
* The program's default settings assume that your VCF has a single FORMAT that applies to ALL SNPs.   If this is not true, you must invoke the program with the additional flag `-f <numFormats>` , where `<numFormats>` should be an integer > 1.  Any integer greater than 1 is sufficient to cause the program to check the format each time. 
* See note above about assumption of NO header in the population designation file.


## Example data files provided here
An example VCF and population designation file are provided in the `ExampleDataFiles/` directory here.  The VCF is a subset of a much larger file from the data archive of Schilling et al. 2018 (_Genes_ 2018, 9(6), 274).  
The original publication is freely available at: [https://doi.org/10.3390/genes9060274](https://doi.org/10.3390/genes9060274)
And the data archive is at: [http://bit.ly/2s6jeIf](http://bit.ly/2s6jeIf)

After cloning this repo, you should be able to run the program (with this directory 
as your current working directory) by using the command syntax given above.  
An example that you can use to test it quickly on your system is the following:

```
./VCFtoSummStats -V ExampleDataFiles/Small_hmel2.5.30f4.vcf.gz -P ExampleDataFiles/popFileHmel.txt
```

As that example illustrates, the program can directly read `.gz` compressed VCF files.  
Note that it relies on the file extension being accurate.

The main results of that command, assuming it runs successfully on your system, 
will be output to a file: `Small_hmel2.5.30f4.vcf.gz_Unfiltered_Summary.tsv`.


## If you need to download the boost libraries

If you are running Linux, open a terminal and use the following command:
```
sudo apt install libboost-all-dev
```

If you are running MacOS, the easiest way is using `brew` in your Terminal:
```
brew install boost
```
If you don't have [homebrew](https://brew.sh/) installed on your Mac, I recommend getting it.  
Chances are, if you are interested in `VCFtoSummStats`, you're likely to be someone who would be able 
to use `brew` to easily install lots of other command-line tools and programs that don't come with Xcode.
