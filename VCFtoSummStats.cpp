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

const bool DEBUG = true;

// functions
void openFileStreams(int argc, char *argv[], ifstream& VCFfile, ifstream& PopulationFile)
{
    const int expectedArgNum = 2;
    string progname = argv[0];
    string message = "\nError!  Please supply two file names as command line arguments,\n\tin the following way:\n\t" + progname + " -V NameOfVCFfile -P NameOfPopulationFile\n\n";
    if ( argc < expectedArgNum ) {
        cout << message;
        exit(-1);
    }
    
    string vcfName, popFileName;
    // parse command line options:
    int flag;
    while ((flag = getopt(argc, argv, "V:P:?")) != -1) {
        switch (flag) {
            case 'V':
                vcfName = optarg;
                break;
            case 'P':
                popFileName = optarg;
                break;
            default: /* '?' */
                cout << message;
                exit(-1);
        }
    }
    
    if ( DEBUG ) {
        cout << "\nI found the following names: " << vcfName << "\tand\t" << popFileName << endl;
    }
    
    // open file streams:
    
}


int main(int argc, char *argv[])
{
    // open file streams for reading:
    // VCF file and Population File
    ifstream VCFfile, PopulationFile;
    openFileStreams(argc, argv, VCFfile, PopulationFile);
    
    // parse information from files needed for figuring out what columns to use
        // (1) identify relevant columns for analysis in VCF file
        // (2) create cross referencing between population file and
    
    
    
    
    
    return 0;
}


