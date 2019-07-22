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

#include "VCFtoSummStats.hpp"

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <unistd.h>
#include <map>
#include <time.h>
#include <math.h>
#include <cstring>
using namespace std;

// for boost libraries for decompressing:
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/bzip2.hpp>


// global variables
const int NUM_META_COLS = 9;    // exected number of fields of data prior to samples in VCF
const char FORMAT_DELIM_DEFAULT = ':'; // expected delimiter of subfields of FORMAT column of VCF
const int MAX_SUBFIELDS_IN_FORMAT_DEFAULT = 30;
const int GT_OPS_CODE = 0, DP_OPS_CODE = 1, GQ_OPS_CODE = 2, PL_OPS_CODE = 3, SKIP_OPS_CODE = 9;
    // the latter are FORMAT parsing codes
const int ENTRIES_IN_PL = 3; // number of separate numbers in PL part of format
const string MISSING_DATA_INDICATOR = "NA";
bool VERBOSE = false;
const size_t MAX_BUFFER_SIZE = 512; // length of char arrays used as buffers
const char VCF_DELIM = '\t'; // VCF files must be tab delimited
const double OVERALL_DP_MIN_THRESHOLD_DEFAULT = 2.0;
double OVERALL_DP_MIN_THRESHOLD;
size_t MAX_TOKEN_LENGTH = 80;


int main(int argc, char *argv[])
{
    // for filtering_streambuf:
    using namespace boost::iostreams;

    clock_t startTime = clock();  // for tracking performance

    // variables for command line arguments:
    int numSamples, numPopulations, numFields, numFormats, firstDataLineNumber = -1;
    int maxSubfieldsInFormat = MAX_SUBFIELDS_IN_FORMAT_DEFAULT;
    unsigned long int VCFfileLineCount = 0;
	bool popFileHeader;
    char formatDelim = FORMAT_DELIM_DEFAULT;
    string vcfName, popFileName;
    // data file streams:
    ifstream PopulationFile;    // population and sample designations
    ofstream outputFile;
    filtering_streambuf<input> myVCFin;     // filter for VCF for dealing with compression
    ifstream vcfUnfiltered; // needed to read in unfiltered

#ifdef DEBUG
    string progname = argv[0];
    cout << "\n\t" << progname << " is running!\n\n";
#endif

    // create cross referencing for population membership by sample:
    map<string, int> mapOfPopulations;      // key = population ID, value = integer population index

	// parse command line options and open file streams for reading:
    parseCommandLineInput(argc, argv, PopulationFile, popFileHeader, numSamples, numPopulations, numFields, numFormats, formatDelim, maxSubfieldsInFormat, vcfName, popFileName, mapOfPopulations );

    createVCFfilter( myVCFin, vcfName, vcfUnfiltered );    // create filter
    istream VCFfile( &myVCFin );            // create stream from filter
    if ( !VCFfile.good() ) {
        cerr << "\nError!  istream VCFfile is not good!\n";
    }

    // create cross referencing for population membership by sample:
    //map<string, int> mapOfPopulations;      // key = population ID, value = integer population index
    map<string, int> mapOfSamples;          // key = sample ID, value = integer representing population index
    int numSamplesPerPopulation[numPopulations];    // for later frequency calculations
    //makePopulationMap( mapOfPopulations, numPopulations, popFileName );
    assignPopIndexToSamples( mapOfPopulations, mapOfSamples, PopulationFile, numSamplesPerPopulation, numPopulations, numSamples  );

    // assign each sample column in the VCF to a population:
    int *populationReference;
    populationReference = new int[numSamples];
    bool success = assignSamplesToPopulations(VCFfile, numSamples, numFields, mapOfSamples, populationReference, VCFfileLineCount, firstDataLineNumber);

#ifdef DEBUG
    if ( success ) {
        cout << "\nassignSamplesToPopulations() exited cleanly\n\n";
    }
    cout << "VCFfileLineCount after assignSamplesToPopulations() is: \t" << VCFfileLineCount << endl;
#endif

    // if all has gone well to this point, the output file can be constructed:
    setUpOutputFile( outputFile, vcfName, numPopulations, mapOfPopulations );

    // after that function call, the  VCFfile stream has pointed
    // to the first entry of the first line of data

    // go through data and calculate allele frequencies:
    parseActualData( VCFfile, numFormats, formatDelim, maxSubfieldsInFormat, VCFfileLineCount, outputFile, numSamples, numPopulations, populationReference, vcfName );

	// cleanup: close files:
	PopulationFile.close();
    outputFile.close();
	// free memory:
	//delete mySamples;


#ifdef DEBUG
		cout << "\nI ran!!\n\n";
#endif
    clock_t endTime = clock();
    int minutes;
	double seconds;
    convertTimeInterval( (endTime - startTime), minutes, seconds);
    cout << "\nIt took " << minutes << "min., " << seconds << "sec."  << " to run.\n";

    return 0;
}


// --------------------- function definitions --------------------------- //
// --------------------- in alphabetical order -------------------------- //
void assignPopIndexToSamples( map<string, int>& mapOfPopulations, map<string, int>& mapOfSamples, ifstream& PopulationFile, int numSamplesPerPopulation[], int numPopulations, int numSamples )
{
    string sampleID, popMembership;
    int popIndex, count = 0;
    for ( int i = 0; i < numPopulations; i++ )
        numSamplesPerPopulation[i] = 0; // initialize to make sure it starts at zero
    while ( PopulationFile >> sampleID >> popMembership ) {
        // get index from population map key-value pair:
        popIndex = mapOfPopulations[ popMembership ];
        // insert index assign that index to the sampleID in the map of samples
        mapOfSamples[ sampleID ] = popIndex;
        numSamplesPerPopulation[ popIndex ]++;
    }

#ifdef DEBUG
    // view number of samples per population:
    cout << "\nNumber of samples per population:\n";
    for ( int i = 0; i < numPopulations; i++ ) {
        cout << "\t" << i << ":\t" << numSamplesPerPopulation[i] << endl;
    }

    // view maps created by last two function calls:
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
#endif


}


bool assignSamplesToPopulations(istream& VCFfile, int numSamples, int numFields, map<string, int> mapOfSamples, int *populationReference, unsigned long int& VCFfileLineCount, int& firstDataLineNumber )
{
    int count = 0, firstSampleCol = (numFields - numSamples + 1);
    int popIndex;
    string x;

#ifdef DEBUG
        if ( firstSampleCol != 10 ) { // expectation based upon VCF format standards
            cout << "\nError!  First sample column was NOT estimated to be 10th field!\n";
            cout << "\t Aborting ...\n\n";
            exit(-2);
        }
#endif
    // First: get to header row (past meta-rows) in VCF file:
    while ( VCFfile >> x ) {
        VCFfileLineCount++; // incremement line counter; though the while loop goes one "word"
                            // at a time, the clauses below move through lines;
        if ( x.substr(0,2) == "##" ) {
            VCFfile.ignore(unsigned(-1), '\n'); // move to next line
        } else if ( x == "#CHROM" ) {
            firstDataLineNumber = VCFfileLineCount + 1;
            // this is the header row after the meta-data header lines
            for ( int i = 0; i < NUM_META_COLS ; i++ ) {
                // advance to first sample header:
#ifdef DEBUG
                    if ( i == (NUM_META_COLS - 1) )
                        cout << "\nYour VCF's last meta-field and some of the sample fields:\n" << x;
#endif
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

#ifdef DEBUG
                    if ( popIndex != iter->second ) {
                        cout << "\nError!  maps aren't working like you think!\n";
                        exit(-3);
                    }
                    if ( count % 100 == 0 || count == (numSamples - 1))
                        cout << " ... " << sampleID << ", popIndex=" << populationReference[count];
#endif

                // advance the VCF stream pointer to the next string
                if ( count < (numSamples - 1) )
                    VCFfile >> sampleID;
                else
                    VCFfile.ignore(unsigned(-1), '\n'); // don't yet read in first entry of next line, but set the stage to
            }

#ifdef DEBUG
                cout << "\n\nMost recently obtained string from VCFfile stream: " << sampleID << "\n\n";
                if ( mapOfSamples.find("foobar") == mapOfSamples.end() )
                    cout << "\nBogus call to mapOfSamples returned mapOfSamples.end()" << endl;
#endif
            if ( (VCFfileLineCount + 1) != firstDataLineNumber ) {
                cout << "\nError in assignSamplesToPopulations():\n\t";
                cout << "VCFfileLineCount + 1 (" << (VCFfileLineCount + 1) << ") != firstDataLineNumber (" << firstDataLineNumber << ")\n\t";
                cout << "Aborting ... \n";
                exit(-2);
            }

            return true;
        } else {
            cout << "\nError!  VCF file not structured as expected!\n";
            cout << "I did NOT find a header row starting with #CHROM\n\t Aborting ...\n\n";
            exit(-2);
        }

    }

    // execution should never reach here unless VCF file has ONLY ## rows
    cout << "\nError!  assignSamplesToPopulations() exited with status 'false'.\n\t";
    cout << "--> Please check that VCF has ## meta rows followed by one\n\t";
    cout << "header row starting with #CHROM, followed by SNP datal rows.\n\tAborting ...\n\n";
    exit(-2);

    return false;
}


inline int calculateMedian( int values[], int n, int ignoreFirst )
{
	int medianSpot = ignoreFirst + ((n - ignoreFirst)/2);
	sort( values, values + n );
	return( values[ medianSpot ] );
}


void calculateSummaryStats( istream& VCFfile, ofstream& outputFile, int numTokensInFormat, int GTtoken, int DPtoken, int GQtoken, int PLtoken, bool lookForDP, bool lookForGQ, bool lookForPL, char formatDelim, int formatOpsOrder[], int numSamples, int numPopulations, unsigned long int VCFfileLineCount, int* populationReference )
{
    int homoRefCount = 0, homoAltCount = 0, hetCount = 0, altAlleleCounts[numPopulations];
    int validSampleCounts[numPopulations], DPvalues[numSamples], GQvalues[numSamples];
    int PLvalues[ (numSamples * ENTRIES_IN_PL) ];
    // initialize all array values to zero:
    for ( int i = 0; i < numPopulations; i++ ) {
        altAlleleCounts[i] = 0;
        validSampleCounts[i] = 0;
    }
    if ( lookForDP ) {
        for ( int i = 0; i < numSamples; i++ )
            DPvalues[i] = 0;
    }
    if ( lookForGQ ) {
        for ( int i = 0; i < numSamples; i++ )
            GQvalues[i] = 0;
    }
    if ( lookForPL ) {
        for ( int i = 0; i < (numSamples * ENTRIES_IN_PL); i++ )
            PLvalues[i] = 0;
    }

    // loop over all columns of data:
    int sampleCounter = 0, operationCode, popIndex;
    //size_t pos, strStart[numTokensInFormat], strLen[numTokensInFormat];
    size_t tokenLength;
    char checkGTsep1 = '/', checkGTsep2 = '|'; // the only two expected separators
    char allele1, allele2;
	int DPnoCall = 0, GQnoCall = 0;
    char* charpt, dummyChar;
    char tokenHolder[MAX_TOKEN_LENGTH];
    for ( sampleCounter = 0; sampleCounter < numSamples; sampleCounter++ ) {
        
        popIndex = populationReference[ sampleCounter ];
        
//        VCFfile.get( dummyChar );
//#ifdef DEBUG
//        if ( dummyChar != '\t' ) {
//            cerr << "\nYo, file parsing not working way you think.\n\tIn calculateSummaryStats():\n\t";
//            cerr << "dummyChar = " << dummyChar << endl;
//        }
//#endif
        
        // parse the current sample:
        for ( int tokeni = 0; tokeni < numTokensInFormat; tokeni++ ) {
            VCFfile.get( dummyChar ); // always have to clear the delims
            if ( tokeni < ( numTokensInFormat - 1 ) )
                VCFfile.get( tokenHolder, MAX_TOKEN_LENGTH, formatDelim ); // get up to next ':'
            else if ( sampleCounter < ( numSamples - 1 ) )
                VCFfile.get( tokenHolder, MAX_TOKEN_LENGTH, VCF_DELIM ); // get up to next '\t'
            else
                VCFfile.get( tokenHolder, MAX_TOKEN_LENGTH, '\n' ); // last possible one
            // get operation code:
            operationCode = formatOpsOrder[tokeni];
            tokenLength = strlen( tokenHolder );
            if ( operationCode == GT_OPS_CODE ) {
                // parse the genotype data and add to correct population
                allele1 = tokenHolder[0];
                allele2 = tokenHolder[2]; // for biallelic SNPS, it should go like this always!

                // considering the diploid genotype, there are 9 options:
                if ( allele1 == '0' ) {
                    if ( allele2 == '0' ) {
                        homoRefCount++;
                        validSampleCounts[popIndex] += 2; // diploid; no alt alleles
                    } else if ( allele2 == '1' ) {
                        hetCount++;
                        validSampleCounts[popIndex] += 2; // diploid
                        altAlleleCounts[popIndex]++; // one alt allele
                    } else {
                        // only allele1 was valid/called:
                        validSampleCounts[popIndex]++;
                    }
                } else if ( allele1 == '1' ) {
                    if ( allele2 == '0' ) {
                        hetCount++;
                        validSampleCounts[popIndex] += 2; // diploid
                        altAlleleCounts[popIndex]++; // one alt allele
                    } else if ( allele2 == '1' ) {
                        homoAltCount++;
                        validSampleCounts[popIndex] += 2; // diploid
                        altAlleleCounts[popIndex] += 2; // two alt alleles
                    } else {
                        // allele 2 was not valid/called
                        validSampleCounts[popIndex]++;
                        altAlleleCounts[popIndex]++; // one alt allele
                    }
                } else {
                    // allele 1 was not valid/called
                    if ( allele2 == '0' || allele2 == '1' ) {
                        validSampleCounts[popIndex]++;
                        if ( allele2 == '1' )
                            altAlleleCounts[popIndex]++; // one alt allele
                    }
                }

                // now recording allele counts:

//#ifdef DEBUG

                if ( tokenHolder[1] != checkGTsep1 && tokenHolder[1] != checkGTsep2 ) {
                    cerr << "\nError in calculateSummaryStats():\n\tGT token ";
                    cerr << "does not have expected character (" << checkGTsep1 << " or " << checkGTsep2 << ") between alleles.\n\t";
                    cerr << "I found: " << tokenHolder[1] << ", and the whole token was:\n\t";
                    
                    fprintf(stderr, "[start]%s[end], length = %lu\n", tokenHolder, tokenLength);
                    fprintf(stderr, "Sample counter = %i\n", sampleCounter);
                    
                    cerr << "Aborting ... \n\n";
                    exit(-1);
                }
//                if ( sampleCounter % 100 == 0 ) {
//                    cout << "\nsample " << (sampleCounter+1) << ", loop count " << tokeni << ", GT is " << token << endl;
//                    cout << "allele1 = " << allele1 << ", allele2 = " << allele2 << "\t" << allele1.length() << "\t" << allele2.length() << endl;
//                    cout << "homoRefCount = " << homoRefCount << ", hetCount = " << hetCount << ", homoAltCount = " << homoAltCount << endl;
//                    cout << "Valid sample counts and alt allele counts by popn:\n";
//                    for ( int j = 0; j < numPopulations; j++ )
//                        cout << "\t" << validSampleCounts[j] << ", " << altAlleleCounts[j];
//                    cout << endl;
//                }
//#endif

            } else if ( operationCode == DP_OPS_CODE && lookForDP ) {
                // add the DP data to DP array
				if ( tokenHolder[0] == '.' && tokenLength == 1 ) {
					DPvalues[sampleCounter] = -1;
					DPnoCall++;
				} else {
                    //fprintf(stdout, "\nDP tokenHolder is: [%s] with length %lu\n", tokenHolder, tokenLength);
                	DPvalues[sampleCounter] = stoi(tokenHolder);
				}

//                cout << "\nsample " << (sampleCounter+1) << ", loop count " << tokeni << ", DP is " << token << endl;
            } else if ( operationCode == GQ_OPS_CODE && lookForGQ ) {
                // add the GQ data to the GQ array
				if ( tokenHolder[0] == '.' && tokenLength == 1 ) {
					GQvalues[sampleCounter] = -1;
					GQnoCall++;
				} else {
                    //fprintf(stdout, "\nGQ tokenHolder is: [%s] with length %lu\n", tokenHolder, tokenLength);
					GQvalues[sampleCounter] = stoi(tokenHolder);
				}
//                cout << "\nsample " << (sampleCounter+1) << ", loop count " << tokeni << ", GQ is " << token << endl;
            } else if ( operationCode == PL_OPS_CODE && lookForPL ) {
                
                parsePL( tokenHolder );
                
            }

            // otherwise just skip it

            // delete the substring
            //currentSample.erase(0, pos + formatDelim.length()); // old way replaced by prior for loop
        }  // end of loop over tokens in sample

    }  // end of for() loop over numSamples; used to be while() loop over lineStream

    // error checking:
    if ( sampleCounter != numSamples ) {
        cout << "\nError in calculateSummaryStats():\n\tline parsing did not give numSamples number of loops.\n\t";
        cout << "sampleCounter = " << sampleCounter << ", but numSamples = " << numSamples;
        cout << "\n\tThis suggests inconsistencies in VCF file construction\n\twith uneven numbers of samples per row";
        cout << "\n\tAborting ... ";
        exit(-5);
    }

    // calculate stats
    // record stats
    // here is the order of remaining columns to calculate and add to ofstream outputFile:
    // medianDP        medianGQ        homoRefCount    hetCount        homoAltCount
    // plus one column for each population named ALT_SNP_FREQ_popName
    int median;
    // medianDP:
    if ( lookForDP && (DPnoCall < numSamples ) ) {
        median = calculateMedian( DPvalues, numSamples, DPnoCall );
        outputFile << "\t" << median;
    } else {
        outputFile << "\t" << MISSING_DATA_INDICATOR;
    }
    // median GQ:
    if ( lookForGQ && (GQnoCall < numSamples) ) {
        median = calculateMedian( GQvalues, numSamples, GQnoCall );
        outputFile << "\t" << median;
    } else {
        outputFile << "\t" << MISSING_DATA_INDICATOR;
    }
    // diploid genotype counts:
    outputFile << "\t" << homoRefCount << "\t" << hetCount << "\t" << homoAltCount;
    double freq;
    for ( int i = 0; i < numPopulations; i++ ) {
        if ( !validSampleCounts[i] ) {
          // no div by zero
            freq = std::numeric_limits<double>::quiet_NaN();
        } else {
            freq = static_cast<double>( altAlleleCounts[i] ) / static_cast<double>( validSampleCounts[i] );
        }
        outputFile << "\t" << freq << "\t" << validSampleCounts[i];
    }

    // outputFile << endl;  not needed here; this is done in parseActualData()

}


inline void checkFormatToken( char* token, int& GTtoken, int& DPtoken, int& GQtoken, int& PLtoken, int subfieldCount  )
{
    // record sub-field:
    if ( token[0] == 'G' && token[1] == 'T' )
        GTtoken = subfieldCount;
    else if ( token[0] == 'D' && token[1] == 'P' )
        DPtoken = subfieldCount;
    else if ( token[0] == 'G' && token[1] == 'Q' )
        GQtoken = subfieldCount;
    else if ( token[0] == 'P' && token[1] == 'L' )
        PLtoken = subfieldCount;
}


void convertTimeInterval( clock_t myTimeInterval, int& minutes, double& seconds)
{
    double totalSeconds = (static_cast<double>( myTimeInterval )) / (static_cast<double>(CLOCKS_PER_SEC));
    long unsigned int totSecondsInt = static_cast<long unsigned int>( totalSeconds );

    minutes = totSecondsInt / 60;
	seconds = totalSeconds - (static_cast<double> (minutes * 60));

}


void createVCFfilter( boost::iostreams::filtering_streambuf<boost::iostreams::input>& myVCFin, string vcfName, ifstream& vcfUnfiltered )
{
    // boost libraries for filtering_streambuf
    using namespace boost::iostreams;

    // must open file:
    vcfUnfiltered.open(vcfName, ios_base::in | ios_base::binary);


    // find the file extension so we know what kind of filter, if any, to use:
    string filext;
    size_t dotPos, endPos;
    dotPos = vcfName.find_last_of( "." );
    endPos = vcfName.length();
    filext = vcfName.substr( dotPos, (endPos - dotPos));
#ifdef DEBUG
    cout << "\nfile extension on VCF file is " << filext << endl;
#endif
    //Read from the first command line argument, assume it's gzipped

    // use file extension to build filter:
    if ( filext == ".gz" ) {
        myVCFin.push(gzip_decompressor());
    } else if ( filext == ".bz2" ) {
        myVCFin.push(bzip2_decompressor());
    } else if ( filext != ".vcf" ) {
        cerr << "\nError!!  File extension '" << filext << "' not recognized!" << endl;
        cerr << "\n\tAborting ... \n\n";
        exit(-1);
    }

    // make the file the input
    myVCFin.push( vcfUnfiltered );
}


void determineFormatOpsOrder( int numTokensInFormat, int GTtoken, int DPtoken, int GQtoken, int PLtoken, bool lookForDP, bool lookForGQ, bool lookForPL, char formatDelim, int formatOpsOrder[], int maxSubfieldsInFormat )
{
    // first a safety check:
    if ( maxSubfieldsInFormat < numTokensInFormat ) {
        cout << "\nError in determineFormatOpsOrder():\n\tmaxSubfieldsInFormat (";
        cout << maxSubfieldsInFormat << ") < number of subfields in your VCF's FORMAT (";
        cout << numTokensInFormat << ")\n";
        cout << "\t--> Call program again with invocation provided by the wrapper\n\t";
        cout << "plus -S " << numTokensInFormat << "\n\tAborting ...\n";
        exit(-4);
    }

#ifdef DEBUG
    if ( lookForGQ ) {
        if ( GTtoken == GQtoken ) {
            cout << "\nError in determineFormatOpsOrder()!!\n\tGTtoken(" << GTtoken << ") == GQtoken (" << GQtoken << ")\n";
            exit(-4);
        }
    }
    if ( lookForDP ) {
        if ( GTtoken == DPtoken ) {
            cout << "\nError in determineFormatOpsOrder()!!\n\tGTtoken(" << GTtoken << ") == DPtoken (" << DPtoken << ")\n";
            exit(-4);
        }
    }
    if ( lookForGQ && lookForDP ) {
        if ( GQtoken == DPtoken ) {
            cout << "\nError in determineFormatOpsOrder()!!\n\tGQtoken(" << GQtoken << ") == DPtoken (" << DPtoken << ")\n";
            exit(-4);
        }
    }
#endif

    int index;
    for ( int i = 0; i < numTokensInFormat; i++ ) {
        index = i + 1;  // token indexes start at 1
        if ( index == GTtoken ) {
            formatOpsOrder[i] = GT_OPS_CODE;
        } else if ( index == DPtoken ) {
            formatOpsOrder[i] = DP_OPS_CODE;
        } else if ( index == GQtoken ) {
            formatOpsOrder[i] = GQ_OPS_CODE;
        } else if ( index == PLtoken ) {
            formatOpsOrder[i] = PL_OPS_CODE;
        } else {
            formatOpsOrder[i] = SKIP_OPS_CODE;
        }
    }

#ifdef DEBUG
    cout << "\nformatOpsOrder:\n";
    for ( int i = 0; i < numTokensInFormat; i++ ) {
        cout << "\t" << formatOpsOrder[i];
    }
    cout << endl;
#endif
}


inline void errorCheckTokens( int GTtoken, int DPtoken, int GQtoken, int PLtoken, bool& lookForDP, bool& lookForGQ, bool& lookForPL )
{
    if ( GTtoken == -1 ) { // -1 is flag for not set
        cout << "\nError!  GT subfield was not found in FORMAT.\nPlease double-check your format column.\n";
        cout << "If your VCF's FORMAT column uses a subfield delimiter other than the colon (:),\n";
        cout << "indicate that by using the call provided by the wrapper script with the\n";
        cout << "addition of the -D DELIM command line option, where 'DELIM' is replaced\n";
        cout << "by the delimiter your VCF uses.\n\tAborting ...\n\n";
        exit(-3);
    }
    if ( DPtoken == -1 ) { // -1 is flag for not set
        cout << "\n*** WARNING!  DP subfield was not found in FORMAT.\n";
        cout << "The medianDP column in the results file will be filled with NA.\n";
        lookForDP = false;
    } else {
        lookForDP = true;
    }
    if ( GQtoken == -1 ) { // -1 is flag for not set
        cout << "\n*** WARNING!  GQ subfield was not found in FORMAT.\n";
        cout << "The medianGQ column in the results file will be filled with NA.\n";
        lookForGQ = false;
    } else {
        lookForGQ = true;
    }
    if ( PLtoken == -1 ) { // -1 is flag for not set
        cout << "\n*** WARNING!  PL subfield was not found in FORMAT.\n";
        cout << "Any results depending upon PL scores will be filled with NA.\n";
        lookForPL = false;
    } else {
        lookForPL = true;
    }
}



//void makePopulationMap( map<string, int>& mapOfPopulations, int numPopulations, string popFileName )
//{
//    // helpful code for working with maps borrowed from:
//    // https://thispointer.com/stdmap-tutorial-part-1-usage-detail-with-examples/
//
//    string uniquePopFileName = popFileName + "_UniquePopFile.txt";
//    ifstream uniquePopFile( uniquePopFileName );
//    string popID;
//    int index = 0;
//
//    // get unique populations
//    while( uniquePopFile >> popID ) {
//        mapOfPopulations[ popID ] = index++;
//    }
//
//    uniquePopFile.close();
//
//    if ( index != numPopulations ) {
//        cout << "\nError!  numPopulations (" << numPopulations << ") != number of keys (" << index << ") in mapOfPopulations!\n\tExiting!\n\n";
//        exit(-1);
//    }
//}


double extractDPvalue( char* INFObuffer, bool& lookForDPinINFO )
{
    char dumc;
    int count = 0, buffCount, valCount;
    bool DPfound = false;
    char holdValueAsChar[80];
    double DPval;
    
    do {
        if ( INFObuffer[count] == 'D' ) {
            // could be DP
            if ( INFObuffer[(count+1)] == 'P' ) {
                if ( INFObuffer[(count+2)] == ' ' || INFObuffer[(count+2)] == '=' ) {
                    // found it;
                    DPfound = true;
                    // because spaces are allowed in VCF standard, it could be
                    // DP = num or DP= num or DP =num or DP=num
                    buffCount = count+2;
                    valCount = 0;
                    while ( INFObuffer[ buffCount ] == ' ' || INFObuffer[ buffCount ] == '=' ) {
                        ++buffCount; // get buffCount to first character that is past the equals sign
                    }
                    while ( INFObuffer[ buffCount ] != ';' && INFObuffer[ buffCount ] != '\t' && INFObuffer[ buffCount ] != ',') {
                        holdValueAsChar[ valCount++ ] = INFObuffer[ buffCount++ ];
                    }
                    holdValueAsChar[ valCount ] = '\0';  // terminate with null string
                    if ( !valCount ) {
                        cerr << "\nError in extractDPvalue():\n\tDP found in INFO but no value found following it!\n\tAborting ....\n\n";
                        exit(-5);
                    }
                    DPval = stod( holdValueAsChar );
                }
            }
        }
    } while ( (INFObuffer[++count] != '\0') && (count < MAX_BUFFER_SIZE) && !DPfound );
    
    if ( !DPfound ) {
        cout << "\nWarning!!  No DP found in INFO field...\n";
        lookForDPinINFO = false;
        DPval = std::numeric_limits<double>::quiet_NaN();
//        cout << "\n\tYo it's nan: " << DPval << endl;
    }
    
    
    return DPval;
}


inline size_t getLength( char *myCharArray )
{
    size_t totalLength = 0;
    
    while ( myCharArray[totalLength] != '\0' ) {
        totalLength++;
    }
    
    return totalLength;
}


void parseActualData(istream& VCFfile, int numFormats, char formatDelim, int maxSubfieldsInFormat, unsigned long int& VCFfileLineCount, ofstream& outputFile, int numSamples, int numPopulations, int* populationReference, string vcfName )
{
    char *CHROM, *POS, *ID, *REF, *ALT, *QUAL;
    //double QUAL;
    long int dumCol, SNPcount = 0;
    char dummyChar;
    bool keepThis, checkFormat = true, lookForDP, lookForGQ, lookForPL;
    int numTokensInFormat, GTtoken = -1, DPtoken = -1, GQtoken = -1, PLtoken = -1;
	string discardedLinesFileName = vcfName + "_discardedLineNums.txt";
	ofstream discardedLinesFile( discardedLinesFileName, ostream::out );
    // string oneLine; // old way using linestream
    // the latter ints are for parsing GT = genotype, DP = depth,
    // and GQ = quality sub-fields of the FORMAT column
    int formatOpsOrder[maxSubfieldsInFormat]; // for keeping track of how to parse FORMAT efficiently
    CHROM = new char[MAX_BUFFER_SIZE];
    POS = new char[MAX_BUFFER_SIZE];
    ID = new char[MAX_BUFFER_SIZE];
    REF = new char[MAX_BUFFER_SIZE];
    ALT = new char[MAX_BUFFER_SIZE];
    QUAL = new char[MAX_BUFFER_SIZE];

	discardedLinesFile << "VCFfileLinesNotUsed" << endl; // header row
    // work line by line:
    // stringstream lineStream( "", ios_base::in | ios_base::out ); old way
    // used to be while( getline ... )
    while ( VCFfile.get( dummyChar ) ) {
        VCFfile.putback( dummyChar ); // replace the one that get() got
        SNPcount++; // counter of how many SNP lines have been processed
        VCFfileLineCount++; // counter of how many LINES of VCF file have been processed

        // turn each line into a string stream for easy parsing by whitespace:
        // lineStream.clear();
        // lineStream.str( oneLine );

        // work with meta-col data:
        keepThis = parseMetaColData( VCFfile, SNPcount, checkFormat, numTokensInFormat, GTtoken, DPtoken, GQtoken, PLtoken, lookForDP, lookForGQ, lookForPL, formatDelim, CHROM, POS, ID, REF, ALT, QUAL);

        if ( checkFormat ) {
            determineFormatOpsOrder( numTokensInFormat, GTtoken, DPtoken, GQtoken, PLtoken, lookForDP, lookForGQ, lookForPL, formatDelim, formatOpsOrder, maxSubfieldsInFormat );
        }

        if ( keepThis ) {
            // it is a biallelic SNP
            // print out meta fields:
            outputFile << VCFfileLineCount << "\t" << CHROM << "\t" << POS << "\t" << ID << "\t" << REF << "\t" << ALT << "\t" << QUAL;

            // let's calculate and store data for one line, i.e., one SNP at a time:
            calculateSummaryStats( VCFfile, outputFile, numTokensInFormat, GTtoken, DPtoken, GQtoken, PLtoken, lookForDP, lookForGQ, lookForPL, formatDelim, formatOpsOrder, numSamples, numPopulations, VCFfileLineCount, populationReference );

            // add end of line (done with this line):
            outputFile << endl;
		} else {
			discardedLinesFile << VCFfileLineCount << endl;
            
		}
        VCFfile.ignore(unsigned(-1), '\n'); // go to end of line
        

        if ( numFormats == 1 ) {
            checkFormat = false; // not needed after first SNP
        }
//        if ( SNPcount == 2 )
//            exit(0);
    }

	discardedLinesFile.close();
    delete[] CHROM;
    delete[] POS;
    delete[] ID;
    delete[] REF;
    delete[] ALT;
    delete[] QUAL;
}


void parseCommandLineInput(int argc, char *argv[], ifstream& PopulationFile, bool& popFileHeader, int& numSamples, int& numPopulations, int& numFields, int& numFormats, char& formatDelim, int& maxSubfieldsInFormat, string& vcfName, string& popFileName, map<string, int>& mapOfPopulations )
{
	const int expectedMinArgNum = 4;
	string progname = argv[0];
    string* uniquePopulationNames;
	string message = "\nError!  Please supply two file names as command line arguments,\n\tin the following way (note flags -V and -P):\n\t" + progname + " -V NameOfVCFfile -P NameOfPopulationFile\n\n";
    bool numFormatsSet = false, vcfNameSet = false, popFileNameSet = false;
    if ( argc < expectedMinArgNum ) {
		cerr << message;
		exit(-1);
	}
    // default or automatic values:
    popFileHeader = false;  // default is NO header
    numFormats = 1;         // default is same FORMAT for every SNP
    OVERALL_DP_MIN_THRESHOLD = OVERALL_DP_MIN_THRESHOLD_DEFAULT;

	// parse command line options:
	int flag;
    while ((flag = getopt(argc, argv, "V:P:Hf:D:S:vd:")) != -1) {
		switch (flag) {
			case 'V':
				vcfName = optarg;
                vcfNameSet = true;
				break;
			case 'P':
				popFileName = optarg;
                popFileNameSet = true;
				break;
			case 'H':
				popFileHeader = true;
				break;
			case 'f':
                numFormats = atoi(optarg);
                numFormatsSet = true;
                break;
            case 'D':
                formatDelim = optarg[0];
                break;
            case 'S':
                maxSubfieldsInFormat = atoi(optarg);
                break;
            case 'v':
                VERBOSE = true;
                break;
            case 'd':
                OVERALL_DP_MIN_THRESHOLD = stod(optarg);
                break;
            default: /* '?' */
				exit(-1);
		}
	}

    if ( !popFileNameSet || !vcfNameSet ) {
        cerr << message;
        exit(-1);
    }
    
    cout << "\nOVERALL_DP_MIN_THRESHOLD is " << OVERALL_DP_MIN_THRESHOLD << endl;

    parsePopulationDesigFile( popFileName, numSamples, numPopulations, mapOfPopulations, popFileHeader );

    numFields = NUM_META_COLS + numSamples;



	// open file streams and check for errors:
//    if ( !VCFfile.good() ) {
//        cout << "\nError in parseCommandLineInput():\n\tVCF file name '" << vcfName << " 'not found!\n\t--> Check spelling and path.\n\tExiting ... \n\n";
//        exit( -1 );
//    }
	PopulationFile.open( popFileName );
	if ( !PopulationFile.good() ) {
		cout << "\nError in parseCommandLineInput():\n\tPopulation file name '" << popFileName << "' not found!\n\t--> Check spelling and path.\n\tExiting ... \n\n";
		exit( -1 );
	}

	// error checking on user input; some arguments are mandatory!
    string testString;
#ifdef DEBUG
        cout << "\nMax length of string on this system = " << testString.max_size() << "\n\n";
#endif
	if ( numPopulations < 2 ) {
		cout << "\nError!  numPopulations = " << numPopulations << ", but it has to be >= 2 for this program.\nExiting ...\n\n";
		exit( -1 );
	}
	if ( !numFormatsSet ) {
        cout << "\nWarning!! numFormats (-f) not set on command line.\nAssuming numFormats = " << numFormats << endl;
    }
}


bool parseMetaColData( istream& VCFfile, long int SNPcount, bool checkFormat, int& numTokensInFormat, int& GTtoken, int& DPtoken, int& GQtoken, int& PLtoken, bool& lookForDP, bool& lookForGQ, bool& lookForPL, char formatDelim, char* CHROM, char* POS, char* ID, char* REF, char* ALT, char* QUAL )
{
    int subfieldCount;  // field counter, starting with index of 1
    double DPval;
    char *buffer = new char[MAX_BUFFER_SIZE];
    char *token = new char[MAX_TOKEN_LENGTH];
    //char myDelim = formatDelim;
    char *FILTER, *INFO, *FORMAT;
    size_t pos, REFlength, ALTlength;
    bool keepThis = true;
    char *INFObuffer = new char[MAX_BUFFER_SIZE];
    static bool lookForDPinINFO = true;
    FILTER = new char[MAX_BUFFER_SIZE];
    INFO = new char[MAX_BUFFER_SIZE];
    FORMAT = new char[MAX_BUFFER_SIZE];


    // loop over fields:
    int col = 1;
    //while( col <= NUM_META_COLS ) {
    //    for ( int col = 1; col <= NUM_META_COLS; col++ ) {
    //VCFfile.get( buffer, MAX_BUFFER_SIZE, VCF_DELIM );
    //fprintf(stdout, "col = %i \t buffer = %s\n", col, buffer);
    //switch ( col ) {
    //            case 1:
    VCFfile.get( CHROM, MAX_BUFFER_SIZE, VCF_DELIM );
    VCFfile.ignore(unsigned(-1), VCF_DELIM);
    //                break;
    //            case 2:
    VCFfile.get( POS, MAX_BUFFER_SIZE, VCF_DELIM );
    VCFfile.ignore(unsigned(-1), VCF_DELIM);
    //                break;
    //            case 3:
    VCFfile.get( ID, MAX_BUFFER_SIZE, VCF_DELIM );
    VCFfile.ignore(unsigned(-1), VCF_DELIM);
    //                break;
    //            case 4:
    VCFfile.get( REF, MAX_BUFFER_SIZE, VCF_DELIM );
    REFlength = getLength( REF );
    //                break;
    VCFfile.ignore(unsigned(-1), VCF_DELIM);
    //            case 5:
    VCFfile.get( ALT, MAX_BUFFER_SIZE, VCF_DELIM );
    ALTlength = getLength( ALT );
    VCFfile.ignore(unsigned(-1), VCF_DELIM);
    //                break;
    //            case 6:
    VCFfile.get( QUAL, MAX_BUFFER_SIZE, VCF_DELIM );
    VCFfile.ignore(unsigned(-1), VCF_DELIM);
    //QUAL = stod(buffer);
    //                break;
    //            case 7:                 // this also handles 8!
    VCFfile.get( FILTER, MAX_BUFFER_SIZE, VCF_DELIM );
    // need to prepare to handle INFO column next
    VCFfile.ignore(unsigned(-1), VCF_DELIM); // move to tab
    VCFfile.get( INFObuffer, MAX_BUFFER_SIZE, VCF_DELIM );
    
    //                fprintf(stdout, "INFObuffer is: \t%s[end]\n", INFObuffer);
    if ( lookForDPinINFO ) {
        DPval = extractDPvalue( INFObuffer, lookForDPinINFO );
        //cout << "\tDPval extracted is: \t" << DPval << endl;
        if ( !isnan( DPval ) ) {
            if ( DPval >= OVERALL_DP_MIN_THRESHOLD )
                keepThis = true;
            else
                keepThis = false;
        }
    } else {
        DPval = std::numeric_limits<double>::quiet_NaN();
        //                    cout << "\tYo DPval is " << DPval << " bruh\n\n  **************** \n\n";
    }
    VCFfile.ignore(unsigned(-1), VCF_DELIM); // move to tab
    
    //                col++; // since this case also handles 8!
    //                break;
    
    //            case 9:
    VCFfile.get( FORMAT, MAX_BUFFER_SIZE, VCF_DELIM );
    //cout << "\nFORMAT is " << FORMAT << endl;
    
    
    //                FORMAT = buffer;
    //                break;
    //            default:
    //                cout << "\nError in parseActualData():\n\t";
    //                cout << "case-switch not working as expected.\n\tAborting ...\n\n";
    //                exit(-3);
    //                break;
    //        }
    //        if ( col < NUM_META_COLS )
    //            VCFfile.ignore(unsigned(-1), VCF_DELIM); // move to tab, but not on last one, because
    //            // calcSummaryStats() always clears a character at the beginning of its loop
    //        col++;
    //
    //    }
    //
    
    // status report:
    if ( VERBOSE ) {
        if ( SNPcount % 10000 == 0 ) {
            cout << "\nSNP lines processed so far = " << SNPcount << "; Current SNP is:\n\t";
            cout <<  CHROM << "\t" << POS << "\t" << ID << "\t" << REF << "\t" << ALT << endl;
        }
    }

    if ( checkFormat ) {
        subfieldCount = 0;
        // string splitting, with thanks to:
        // https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
#ifdef DEBUG
        cout << "** Parsing FORMAT **\n\tsubfieldcount\tpos\ttoken\n" << FORMAT << endl;
        cout << "CHROM: " << CHROM << ", INFO: " << INFO << endl;
#endif
        //
        pos = 0;
        char nextChar = FORMAT[pos];
        int subCount;
        while ( nextChar != '\0' ) {
            // get the next subfield:
            subfieldCount++; // number of subfields processed
            subCount = 0;
            do {
                //cout << "\nnextChar is " << nextChar << endl;
                token[subCount++] = nextChar; // build it char by char
                nextChar = FORMAT[++pos]; // get next char
            } while ( nextChar != formatDelim && nextChar != '\0' );
            token[subCount] = '\0'; // terminate string
            // process subfield:
            //cout << "\ntoken is " << token << endl;
            checkFormatToken( token, GTtoken, DPtoken, GQtoken, PLtoken, subfieldCount );
            
            if ( nextChar != '\0' )
                nextChar = FORMAT[++pos]; // get the next character past the delimiter
        }
        if ( !pos ) {
            cerr << "\nError in parseMetaColData():\n\tpos = 0 meaning FORMAT has length zero!\n\tFORMAT = " << FORMAT << endl;
            exit(-1);
        }
        numTokensInFormat = subfieldCount;
        errorCheckTokens( GTtoken, DPtoken, GQtoken, PLtoken, lookForDP, lookForGQ, lookForPL );
            
//        while ( (pos = FORMAT.find( myDelim )) != string::npos ) {
//            subfieldCount++;
//            token = FORMAT.substr(0, pos);
//            // record sub-field if it matches one we're looking for:
//            checkFormatToken( token, GTtoken, DPtoken, GQtoken, subfieldCount );
//#ifdef DEBUG
//                cout << "\t" << subfieldCount << "\t\t" << pos  << "\t" << token << endl;
//#endif
//            FORMAT.erase(0, pos + myDelim.length());
//        }
//        numTokensInFormat = subfieldCount + 1;
//        checkFormatToken( FORMAT, GTtoken, DPtoken, GQtoken, numTokensInFormat ); // check remaining field (after last occurrence of delimiter)
//        // error checking on occurrence of GT:
//        errorCheckTokens( GTtoken, DPtoken, GQtoken, lookForDP, lookForGQ );

#ifdef DEBUG
        cout << "\t" << numTokensInFormat << "\t\tlast\t" << FORMAT << "\twhile loop count = " << subfieldCount << endl;
        cout << "GTtoken = " << GTtoken << "; DPtoken = " << DPtoken << "; GQtoken = " << GQtoken << "; PLtoken = " << PLtoken << "; lookForPL = " << lookForPL << endl;
        
#endif


    }
    // check for bi-allelic SNPs:
    if ( keepThis ) {
        if ( REF[0] == 'N' || ALT[0] == 'N' || ALTlength != 1 || REFlength != 1 ) {
            keepThis = false;
        } else {
            keepThis = true;
        }
    }

//    if ( !keepThis ) {
//#ifdef DEBUG
//        cout << "\nSNP #" << SNPcount << ", ID = " << ID << ", has REF = " << REF << " and ALT = " << ALT << endl;
//#endif
//    }
    
    delete[] INFObuffer;
    delete[] buffer;
    delete[] FILTER;
    delete[] INFO;
    delete[] FORMAT;
    delete[] token;

    return keepThis;
}


inline void parsePL( char* tokenHolder )
{
    
    // here's the description from the VCF file specification, pp. 10-11 of
    // http://samtools.github.io/hts-specs/VCFv4.3.pdf, accessed 7/9/19:
    /*
     GL (Float): Genotype likelihoods comprised of comma
     separated floating point log10-scaled
     likelihoods for all possible genotypes given the set
     of alleles defined in the REF and ALT fields.
     In presence of the GT field the same ploidy is expected;
     without GT field, diploidy is assumed.
     
     PL (Integer): The phred-scaled genotype likelihoods rounded
     to the closest integer, and otherwise defined precisely as
     the GL field.
     */
    
#ifdef DEBUG
    cout << "\ntokenHolder = " << tokenHolder << endl;
    char* tempCharArray = strtok(tokenHolder, ",");
    int i = 0;
    while ( tempCharArray ) {
        cout << "\ttoken part " << ++i << ":\t" << tempCharArray << endl;
        tempCharArray = strtok(NULL, ",");
    }
    cout << endl;
    exit(0);
#endif
    
}


void parsePopulationDesigFile( string fname, int& numSamples, int& numPopulations, map<string,int>& mapOfPopulations, bool popFileHeader )
{
    // goal is to get numSamples, numPopulations, and uniquePopulationNames

    ifstream popMapFile( fname );
    string sampleID, popMembership;
    int countPops = 0, countSamples = 0, lineCount = 0;
    map<string, int> tempSampleMap;

    if ( !popMapFile.good() ) {
        cout << "\nError in parseCommandLineInput():\n\tPopulation file name '" << fname << "' not found!\n\t--> Check spelling and path.\n\tAborting ... \n\n";
        exit( -1 );
    }


    if ( popFileHeader ) {
        popMapFile.ignore(unsigned(-1), '\n'); // skip header line
    }

    while( popMapFile >> sampleID >> popMembership ) {
//        cout << sampleID << "\t" << popMembership << endl;
        ++lineCount;
        if ( tempSampleMap.count( sampleID ) ) {
            cerr << "\nError! Duplicate Sample ID (" << sampleID << ") found!\n\tAborting ...\n";
            exit(-1);
        } else {
            tempSampleMap.insert( pair<string, int>(sampleID, ++countSamples) );
        }
        if ( !mapOfPopulations.count( popMembership ) ) {
            mapOfPopulations.insert( pair<string, int>(popMembership, ++countPops) );
        }


    }
//    cout << "lineCount is " << lineCount << endl;

    popMapFile.close();

    numSamples = countSamples;
    numPopulations = countPops;


    string uniquePopulationNames[numPopulations];
    int counter = 0;
    //cout << "\nInitial state of mapOfPopulations:\n\tkey\tvalue\n";
    // get an alphabetical list of names:
    for( map<string, int>::const_iterator it = mapOfPopulations.begin();
        counter < numPopulations;
        counter++ )
    {
        uniquePopulationNames[counter] = it->first;
        // cout << "\t" << it->first << "\t" << it->second << endl;
        it++;
    }
    counter = 0; // reassign integer designations based upon alphabetical ordering
    string dums;
    for( map<string, int>::const_iterator it = mapOfPopulations.begin();
        counter < numPopulations;
        counter++ )
    {
        dums = uniquePopulationNames[counter];
        // change value:
        mapOfPopulations[ dums ] = counter;
        it++;
    }

    if ( VERBOSE )
        cout << "\nPopulation designations by integer ID:\n";
    for ( int i = 0; i < numPopulations; i++ ) {
        if ( VERBOSE )
            cout << "\tPopulation " << i << " is " << uniquePopulationNames[i] << endl;
        dums = uniquePopulationNames[i];
        if ( mapOfPopulations[ dums ] != i ) {
            cerr << "\nError in parsePopulationDesigFile():\n\tindexes not set up as you expect!\n\tAborting ... \n";
            exit(-1);
        }
    }

}

void setUpOutputFile (ofstream& outputFile, string vcfName, int numPopulations, map<string, int> mapOfPopulations )
{
    string filename = vcfName + "_Unfiltered_Summary" + ".tsv";
    string popHeader, popName, colHeaders, alleleCountHeader;
    int popIndex;
    map<string, int>::const_iterator it = mapOfPopulations.begin();

    // open file for output
    outputFile.open( filename, ofstream::out );
    if ( outputFile.fail() ) {
        cout << "\nError in setUpOutputFile():\n\toutputFile.fail()!\n\t--> Please make sure you have write access to the data file directory.\n\tAborting ... \n\n";
        exit(-4);
    }
    // first several column headers:
    colHeaders = "VCFlineNum\tCHROM\tPOS\tID\tREF\tALT\tQUAL\tmedianDP\tmedianGQ\thomoRefCount\thetCount\thomoAltCount";
    // put first few headers in file:
    outputFile << colHeaders;


    // loop over populations:
    popHeader = "\tALT_SNP_freq_";
    alleleCountHeader = "\trawAlleleCount_";
    for ( int i = 0; i < numPopulations; i++ ) {
        popName = it->first;
        popIndex = it->second;
        if ( popIndex != i ) {
            cout << "\nError in setUpOutputFile():\n\tmap isn't ordered as you expect!\n\tAborting ... \n\n";
            exit(-4);
        }
        outputFile << popHeader << popName << alleleCountHeader << popName;
        it++;
    }
    outputFile << endl;

    //outputFile.close();

}
