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
#include <map>
#include <sstream>
#include "VCFtoSummStats.hpp"
using namespace std;

// global variables
const bool DEBUG = true;  	    // flag for debugging options in code
const int NUM_META_COLS = 9;    // exected number of fields of data prior to samples in VCF;

int main(int argc, char *argv[])
{
	int numSamples, numPopulations, numFields, numFormats = 1;
    long int maxCharPerLine;
	bool popFileHeader;
	
	if ( DEBUG )
		cout << "\nI'm running!\n\n";
	
	// parse command line options and open file streams for reading:
    // VCF file and Population File
    ifstream VCFfile, PopulationFile;
    parseCommandLineInput(argc, argv, VCFfile, PopulationFile, maxCharPerLine, popFileHeader, numSamples, numPopulations, numFields, numFormats);
	
    // parse information from files needed for figuring out what columns to use
	
    // create cross referencing for population membership by sample:
    map<string, int> mapOfPopulations;
    map<string, int> mapOfSamples;
    makePopulationMap( mapOfPopulations, numPopulations );
    assignPopIndexToSamples( mapOfPopulations, mapOfSamples, PopulationFile  );
    
    if ( DEBUG ) {
        string dum;
        int counter = 0;
        cout << "\nmapOfPopulations:\n";
        for( map<string, int>::const_iterator it = mapOfPopulations.begin();
            counter < numPopulations;
            counter++ ) {
            dum = it->first;
            if ( dum.length() == 1 ) {
                dum = dum + "\t\t";
            } else if ( dum.length() < 9 ) {
                dum = dum + "\t";
            }
            cout << "\t" << (counter + 1) << "\tKey:  " << dum << "\tValue:  " << it->second << endl;
            it++;
        }
        cout << "\nmapOfSamples (subset):\n";
        counter = 0;
        for( map<string, int>::const_iterator it = mapOfSamples.begin();
            counter < numSamples;
            counter++ ) {
            if ( counter % 100 == 0 || counter == (numSamples - 1) )
                cout << "\t" << (counter + 1) << "\tKey:  " << it->first << "\t\tValue:  " << it->second << endl;
            it++;
        }
    }
    
    // assign each sample column in the VCF to a population:
    int *populationReference;
    populationReference = new int[numSamples];
    bool success = assignSamplesToPopulations(VCFfile, numSamples, numFields, mapOfSamples, populationReference);
    if ( DEBUG ) {
        if ( success ) {
            cout << "\nassignSamplesToPopulations() exited cleanly\n\n";
        }
    }
    // after that function call, the  VCFfile stream has pointed
    // to the first entry of the first line of data
    
    // go through data and calculate allele frequencies:
    calculateAlleleFrequencies( VCFfile, numFormats );
    
    
    exit(0);
    
    
    // parse VCF and calculate allele frequencies by population
    
    
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
void assignPopIndexToSamples( map<string, int>& mapOfPopulations, map<string, int>& mapOfSamples, ifstream& PopulationFile )
{
    string sampleID, popMembership;
    int popIndex, count = 0;
    while ( PopulationFile >> sampleID >> popMembership ) {
        // get index from population map key-value pair:
        popIndex = mapOfPopulations[ popMembership ];
        // insert index assign that index to the sampleID in the map of samples
        mapOfSamples[ sampleID ] = popIndex;
//        if ( DEBUG ) {
//            if ( ++count % 20 == 1)
//                cout << sampleID << "\t" << popMembership << "\t" << popIndex << "\t" << mapOfSamples[ sampleID ] << endl;
//        }
    }
}


bool assignSamplesToPopulations(ifstream& VCFfile, int numSamples, int numFields, map<string, int> mapOfSamples, int *populationReference)
{
    int count = 0, firstSampleCol = (numFields - numSamples + 1);
    int popIndex;
    string x;
    
    if ( DEBUG ) {
        if ( firstSampleCol != 10 ) { // expectation based upon VCF format standards
            cout << "\nError!  First sample column was NOT estimated to be 10th field!\n";
            cout << "\t Aborting ...\n\n";
            exit(-2);
        }
    }
    
    // First: get to header row (past meta-rows) in VCF file:
    while ( VCFfile >> x ) {
        if ( x.substr(0,2) == "##" ) {
            VCFfile.ignore(unsigned(-1), '\n'); // move to next line
        } else if ( x == "#CHROM" ) {
            // this is the header row after the meta-data header lines
            for ( int i = 0; i < NUM_META_COLS ; i++ ) {
                // advance to first sample header:
                if ( DEBUG ) {
                    if ( i == (NUM_META_COLS - 1) )
                        cout << "\nYour VCF's last meta-field and some of the sample fields:\n" << x;
                }
                VCFfile >> x;
                
            }
            
            string sampleID = x;
            map<string, int>::iterator iter; // for checking existence in map
            for ( count = 0; count < numSamples; count ++ ) {
                // map sample column to population
                iter = mapOfSamples.find( sampleID );
                // check to make sure sampleID is in the map:
                if ( iter == mapOfSamples.end() ) {
                    cout << "\nError!  Sample header '" << sampleID << "' from VCF file not found in mapOfSamples!" << endl;
                    cout << "--> Please check that your population file designates\nsamples EXACTLY as they appear in the VCF." << endl;
                    cout << "\tAborting ... " << endl;
                    exit(-2);
                }
                // get popIndex from map:
                popIndex = mapOfSamples[ sampleID ];
                // store popIndex in array that maps each column to a population:
                populationReference[ count ] = popIndex;
                
                if ( DEBUG ) {
                    if ( popIndex != iter->second ) {
                        cout << "\nError!  maps aren't working like you think!\n";
                        exit(-3);
                    }
                    if ( count % 100 == 0 || count == (numSamples - 1))
                        cout << " ... " << sampleID << ", popIndex=" << populationReference[count];
                }
                
                // advance the VCF stream pointer to the next string
                if ( count < (numSamples - 1) )
                    VCFfile >> sampleID;
                else
                    VCFfile.ignore(); // don't yet read in first entry of next line, but set the stage to
            }
            
            if ( DEBUG ) {
                string oneLine;
                //oneLine = new char[55288];
                cout << "\n\nMost recently obtained string from VCFfile stream: " << sampleID << "\n\n";
//                getline( VCFfile, oneLine );
//                cout << "Results of getline(): " << oneLine << "\n\n";
//                getline( VCFfile, oneLine );
//                cout << "Results of getline(): " << oneLine << "\n\n";
                
            }
            
            if ( DEBUG ) {
                if ( mapOfSamples.find("foobar") == mapOfSamples.end() )
                    cout << "\nBogus call to mapOfSamples returned mapOfSamples.end()" << endl;
            }
            
            return true;
        } else {
            cout << "\nError!  VCF file not structured as expected!\n";
            cout << "I did NOT find a header row starting with #CHROM\n\t Aborting ...\n\n";
            exit(-2);
        }
    }
    
    return false;  // execution should never reach here unless VCF file has ONLY ## rows
}


void calculateAlleleFrequencies(ifstream& VCFfile, int numFormats )
{
    string oneLine;
    long int dumCol, SNPcount = 0;
    bool keepThis, checkFormat = true;
    int numTokensInFormat, GTtoken = -1, DPtoken = -1, GQtoken = -1;
    
    // work line by line:
    
    while ( getline( VCFfile, oneLine ) ) {
        SNPcount++;
//        if ( DEBUG ) {
//            cout << oneLine.substr(0, 20) << endl;
//        }
        // turn each line into a string stream for easy parsing by whitespace:
        stringstream lineStream( oneLine );
        // work with meta-col data:
        keepThis = parseMetaColData( lineStream, SNPcount, checkFormat, numTokensInFormat );
        if ( keepThis ) {
            // it is a biallelic SNP
        }
        if ( numFormats == 1 ) {
            checkFormat = false; // not needed after first SNP
        }
    }
}



void makePopulationMap( map<string, int>& mapOfPopulations, int numPopulations  )
{
    // helpful code for working with maps borrowed from:
    // https://thispointer.com/stdmap-tutorial-part-1-usage-detail-with-examples/
    
    ifstream uniquePopFile( "UniquePopFileTemp.txt" );
    string popID;
    int index = 0;
    
    // get unique populations
    while( uniquePopFile >> popID ) {
        mapOfPopulations[ popID ] = index++;
    }
    
    uniquePopFile.close();
    
    if ( DEBUG ) {
        if ( index != numPopulations ) {
            cout << "\nError!  numPopulations (" << numPopulations << ") != number of keys (" << index << ") in mapOfPopulations!\n\tExiting!\n\n";
            exit(-1);
        }
    }
}


void parseCommandLineInput(int argc, char *argv[], ifstream& VCFfile, ifstream& PopulationFile, long int& maxCharPerLine, bool& popFileHeader, int& numSamples, int& numPopulations, int& numFields, int& numFormats)
{
	const int expectedMinArgNum = 2;
	string progname = argv[0];
	string message = "\nError!  Please supply two file names as command line arguments,\n\tin the following way:\n\t" + progname + " -V NameOfVCFfile -P NameOfPopulationFile\n\n";
	bool maxCharSet = false, popFileHeaderSet = false, numFieldsSet = false;
	bool numSamplesSet = false, numPopulationsSet = false, numFormatsSet = false;
	if ( argc < expectedMinArgNum ) {
		cout << message;
		exit(-1);
	}
	
	string vcfName, popFileName;
	// parse command line options:
	int flag;
    while ((flag = getopt(argc, argv, "V:P:L:H:N:n:F:f:?")) != -1) {
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
            case 'f':
                numFormats = atoi(optarg);
                numFormatsSet = true;
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
    string testString;
    if ( DEBUG ) {
        cout << "\nMax length of string on this system = " << testString.max_size() << "\n\n";
    }
	if ( !maxCharSet ) {
		maxCharPerLine = testString.max_size();
		cout << "\nWarning: maxCharPerLine (maximum line length) is not set.  Using default of " << maxCharPerLine << endl;
    } else if ( maxCharPerLine > testString.max_size() ) {
        cout << "\nError!  maxCharPerLine in your VCFfile is " << maxCharPerLine << ",\n";
        cout << "but max string length on this system is " << testString.max_size() << endl;
        cout << "Please report this to flaxmans@colorado.edu\n\tAborting ... \n\n";
        exit(-1);
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
    if ( !numFormatsSet ) {
        cout << "\nWarning!! numFormats (-f) not set on command line.\nAssuming numFormats = " << numFormats << endl;
    }
}


bool parseMetaColData( stringstream& lineStream, long int SNPcount, bool checkFormat, int& numTokensInFormat )
{
    int dumi;  // field counter, starting with index of 1
    string buffer, myDelim = ":", token;
    string CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT;
    size_t pos;
    bool keepThis;
    
    if ( DEBUG ) {
        cout << "\n *** SNPcount = " << SNPcount << endl;
    }
    
    // loop over fields:
    for ( int col = 1; col <= NUM_META_COLS; col++ ) {
        lineStream >> buffer;
        switch ( col ) {
            case 1:
                CHROM = buffer;
                break;
            case 2:
                POS = buffer;
                break;
            case 3:
                ID = buffer;
                break;
            case 4:
                REF = buffer;
                break;
            case 5:
                ALT = buffer;
                break;
            case 6:
                QUAL = buffer;
                break;
            case 7:
                FILTER = buffer;
                break;
            case 8:
                INFO = buffer;
                break;
            case 9:
                FORMAT = buffer;
                break;
            default:
                cout << "\nError in calculateAlleleFrequencies():\n\t";
                cout << "case-switch not working as expected.\n\tAborting ...\n\n";
                exit(-3);
                break;
        }
    }
    
    cout <<  CHROM << "\t" << POS << "\t" << ID << "\t" << REF << "\t" << ALT << endl;
    
    if ( checkFormat ) {
        if ( DEBUG && SNPcount < 4 ) {
            dumi = 0;
            // mess with string splitting
            // thanks to https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
            cout << "pos\ttoken\n";
            while ( (pos = FORMAT.find( myDelim )) != string::npos ) {
                dumi++;
                token = FORMAT.substr(0, pos);
                cout << pos  << "\t" << token << endl;
                FORMAT.erase(0, pos + myDelim.length());
            }
            token = FORMAT;
            cout << "last\t" << token << "\twhile loop count = " << dumi << endl;
            numTokensInFormat = dumi + 1;
            
        }
    }
    // check for bi-allelic SNPs:
    if ( ALT.length() == 1 ) {
        // it should be biallelic
        keepThis = true;
    } else {
        keepThis = false;
    }
    
    if ( DEBUG ) {
        if ( ALT.length() != 1 ) {
            cout << "\nSNP #" << SNPcount << ", ID = " << ID << ", has REF = " << REF << " and ALT = " << ALT << endl;
        }
    }
    
    return keepThis;
}
