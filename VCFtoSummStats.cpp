// VCFtoSummStats.cpp
// This is a program for making it easy to calculate
// population genetic summary statistics from a VCF.
// The idea is to combine aspects of
// commonly produced SNP tables with population genetic
// statistics in a way that facilitates meta-analyses,
// comparative studies, and theoretical-empirical crosstalk
// in the understanding, use, and interpretation of patterns of
// Hetergeneous Genomic Divergence (HGD)

// this program outputs allele frequencies by deme/population/taxon
// from which a wide variety of summary statistics can be calculated
// or estimated, such as FST, DXY, etc.

// produced by Samuel M. Flaxman & Chris C. R. Smith
// May 2019
// distributed without copyright

// please see accompanying README.md for more information

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <unistd.h>
using namespace std;

// global variables
const bool DEBUG = true;  	// flag for debugging options in code


// function prototypes:
void parseCommandLineInput(int argc, char *argv[], ifstream& VCFfile, ifstream& PopulationFile, int& maxCharPerLine);

int main(int argc, char *argv[])
{
	int maxCharPerLine = 0;
	
	// parse command line options and open file streams for reading:
    // VCF file and Population File
    ifstream VCFfile, PopulationFile;
    parseCommandLineInput(argc, argv, VCFfile, PopulationFile, maxCharPerLine);
	
    // parse information from files needed for figuring out what columns to use
        // (1) identify relevant columns for analysis in VCF file
        // (2) create cross referencing between population file and
    
	
	
	
	VCFfile.close();
	PopulationFile.close();
	
    return 0;
}


// --------------------- function definitions --------------------------- //
// --------------------- in alphabetical order -------------------------- //
void parseCommandLineInput(int argc, char *argv[], ifstream& VCFfile, ifstream& PopulationFile, int& maxCharPerLine)
{
	const int expectedMinArgNum = 2;
	string progname = argv[0];
	string message = "\nError!  Please supply two file names as command line arguments,\n\tin the following way:\n\t" + progname + " -V NameOfVCFfile -P NameOfPopulationFile\n\n";
	bool maxCharNotSet = true;
	if ( argc < expectedMinArgNum ) {
		cout << message;
		exit(-1);
	}
	
	string vcfName, popFileName;
	// parse command line options:
	int flag;
	while ((flag = getopt(argc, argv, "V:P:L:?")) != -1) {
		switch (flag) {
			case 'V':
				vcfName = optarg;
				break;
			case 'P':
				popFileName = optarg;
				break;
			case 'L':
				maxCharPerLine = atoi(optarg);
				maxCharNotSet = false;
				break;
			default: /* '?' */
				exit(-1);
		}
	}
	
	// open file streams and check for errors:
	VCFfile.open( vcfName );
	if ( !VCFfile.good() ) {
		cout << "\nError in parseCommandLineInput():\n\tVCF file name '" << vcfName << " 'not found!\n\t--> Check spelling and path.\n\tExiting ... \n\n";
		exit( -1 );
	}
	PopulationFile.open( popFileName );
	if ( !PopulationFile.good() ) {
		cout << "\nError in parseCommandLineInput():\n\tPopulation file name '" << popFileName << "' not found!\n\t--> Check spelling and path.\n\tExiting ... \n\n";
		exit( -1 );
	}
	if ( maxCharNotSet ) {
		maxCharPerLine = 1000;
		cout << "\nWarning: maxCharPerLine (maximum line length) is not set.  Using default of " << maxCharPerLine << endl;
	}
}
