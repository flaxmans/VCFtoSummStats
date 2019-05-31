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
#include "VCFtoSummStats.hpp"
using namespace std;

// global variables
const bool DEBUG = true;  	// flag for debugging options in code

int main(int argc, char *argv[])
{
	int maxCharPerLine, numSamples, numPopulations, numFields;
	bool popFileHeader;
	
	if ( DEBUG )
		cout << "\nI'm running!\n\n";
	
	// parse command line options and open file streams for reading:
    // VCF file and Population File
    ifstream VCFfile, PopulationFile;
    parseCommandLineInput(argc, argv, VCFfile, PopulationFile, maxCharPerLine, popFileHeader, numSamples, numPopulations, numFields);
	
    // parse information from files needed for figuring out what columns to use
	// figure get IDs and populations from population file
	sample *mySamples;
	mySamples = new sample[numSamples];
	assignSamplesToPopulations(mySamples, VCFfile, PopulationFile, numSamples, numPopulations, numFields);
	// getPopulationInformation( PopulationFile, popFileHeader, numSamples, numPopulations );
	
	// identify relevant columns for analysis in VCF file
	// parseVCFcolumns();
	
	// create cross referencing between population file and
	// calculateSNPstats();
	
	// cleanup: close files:
	VCFfile.close();
	PopulationFile.close();
	// free memory:
	//delete mySamples;
	
	// if (DEBUG)
	//	cout << "\nmaxCharPerLine = " << maxCharPerLine << "; popFileHeader = " << popFileHeader << "; numSamples = " << numSamples << endl;
	
	if ( DEBUG )
		cout << "\nI ran!!\n\n";
	
    return 0;
}


// --------------------- function definitions --------------------------- //
// --------------------- in alphabetical order -------------------------- //
void assignSamplesToPopulations(sample *mySamples, ifstream& VCFfile, ifstream& PopulationFile, int numSamples, int numPopulations, int numFields)
{
	int count = 0, firstSampleCol = (numFields - numSamples + 1);
	int infoColNum;
	string x;
	bool lookingForHeaders = true;
	while ( VCFfile >> x && count < (2 * numFields) ) {
		if ( lookingForHeaders ) {
			if ( x.substr(0,2) == "##" ) {
				VCFfile.ignore(unsigned(-1), '\n'); // move to next line
			} else if ( x.substr(0,1) == "#" ) {
				// this is the header row after the meta-data header lines
				for ( int i = 1; i < firstSampleCol ; i++ ) {
					if ( x == "INFO" ) {
						// this is the info column
						infoColNum = i;
						cout << "AF1\t"
					} else if ( x == "FORMAT" ) {
						// this is the FORMAT column
						formatColNum = i;
					}
					cout << x << "\t";
					VCFfile >> x;
					++count;
				}
				cout << endl;
				VCFfile.ignore(unsigned(-1), '\n');
				lookingForHeaders = false;
			} else {
				++count;
				if ( DEBUG ) {
					if ( x == "FORMAT" ) {
						cout << "\nFORMAT at " << count << "th column" << endl;
					}
					if ( count == (firstSampleCol - 1) ) {
						cout << "\nLast non-sample col is " << x << ", which is " << (firstSampleCol - 1) << "th field" << endl;
					} else if ( count == firstSampleCol ) {
						cout << "\nFirst sample col is " << x << ", which is " << firstSampleCol << "th field" << endl;
					}
					cout << x << "\t";
				}
				if ( count == numFields ) {
					lookingForHeaders = false;
					cout << endl;
				}
			}
		} else {
			for ( int i = 0; i < firstSampleCol; i++ ) {
				cout << x << "\t";
				VCFfile >> x;
				++count;
			}
			cout << endl;
			VCFfile.ignore(unsigned(-1), '\n');
		}
	}
	cout << endl;
}


void parseCommandLineInput(int argc, char *argv[], ifstream& VCFfile, ifstream& PopulationFile, int& maxCharPerLine, bool& popFileHeader, int& numSamples, int& numPopulations, int& numFields)
{
	const int expectedMinArgNum = 2;
	string progname = argv[0];
	string message = "\nError!  Please supply two file names as command line arguments,\n\tin the following way:\n\t" + progname + " -V NameOfVCFfile -P NameOfPopulationFile\n\n";
	bool maxCharSet = false, popFileHeaderSet = false, numFieldsSet = false;
	bool numSamplesSet = false, numPopulationsSet = false;
	if ( argc < expectedMinArgNum ) {
		cout << message;
		exit(-1);
	}
	
	string vcfName, popFileName;
	// parse command line options:
	int flag;
	while ((flag = getopt(argc, argv, "V:P:L:H:N:n:F:?")) != -1) {
		switch (flag) {
			case 'V':
				vcfName = optarg;
				break;
			case 'P':
				popFileName = optarg;
				break;
			case 'L':
				maxCharPerLine = atoi(optarg);
				maxCharSet = true;
				break;
			case 'H':
				popFileHeader = atoi(optarg);
				popFileHeaderSet = true;
				break;
			case 'N':
				numSamples = atoi(optarg);
				numSamplesSet = true;
				break;
			case 'n':
				numPopulations = atoi(optarg);
				numPopulationsSet = true;
				break;
			case 'F':
				numFields = atoi(optarg);
				numFieldsSet = true;
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
	
	// error checking on user input; some arguments are mandatory!
	if ( !maxCharSet ) {
		maxCharPerLine = 1000;
		cout << "\nWarning: maxCharPerLine (maximum line length) is not set.  Using default of " << maxCharPerLine << endl;
	}
	if ( !popFileHeaderSet ) {
		cout << "\nError! popFileHeader option (-H) needs to be set to 0 (false) or 1 (true)\nExiting ... \n\n";
		exit( -1 );
	}
	if ( !numSamplesSet ) {
		cout << "\nError! numSamples option (-N) needs to be set to total number of samples/individuals.\nExiting ... \n\n";
		exit( -1 );
	}
	if ( !numPopulationsSet ) {
		cout << "\nError! numPopulations option (-n) needs to bet set to total number of populations.\nExiting ... \n\n";
		exit( -1 );
	} else if ( numPopulations < 2 ) {
		cout << "\nError!  numPopulations = " << numPopulations << ", but it has to be >= 2 for this program.\nExiting ...\n\n";
		exit( -1 );
	}
	if ( !numFieldsSet ) {
		cout << "\nError! numFields option (-F) needs to be set to number of fields in VCF file.\nExiting ... \n\n";
		exit( -1 );
	}
}
