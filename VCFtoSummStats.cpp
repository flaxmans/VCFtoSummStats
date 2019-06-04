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
#include <algorithm>
using namespace std;

// global variables
const int NUM_META_COLS = 9;    // exected number of fields of data prior to samples in VCF
const string FORMAT_DELIM_DEFAULT = ":"; // expected delimiter of subfields of FORMAT column of VCF
const int MAX_SUBFIELDS_IN_FORMAT_DEFAULT = 30;
const int GT_OPS_CODE = 0, DP_OPS_CODE = 1, GQ_OPS_CODE = 2, SKIP_OPS_CODE = 3;
    // the latter are FORMAT parsing codes
const string MISSING_DATA_INDICATOR = "NA";

int main(int argc, char *argv[])
{
    // variables for command line arguments:
    int numSamples, numPopulations, numFields, numFormats = 1, firstDataLineNumber = -1;
    int maxSubfieldsInFormat = MAX_SUBFIELDS_IN_FORMAT_DEFAULT;
    unsigned long int maxCharPerLine, VCFfileLineCount = 0;
	bool popFileHeader;
    string formatDelim = FORMAT_DELIM_DEFAULT, vcfName, popFileName;
    // data file streams:
    ifstream VCFfile, PopulationFile;
    ofstream outputFile;
    
#ifdef DEBUG
        string progname = argv[0];
		cout << "\n\t" << progname << " is running!\n\n";
#endif
	
	// parse command line options and open file streams for reading:
    parseCommandLineInput(argc, argv, VCFfile, PopulationFile, maxCharPerLine, popFileHeader, numSamples, numPopulations, numFields, numFormats, formatDelim, maxSubfieldsInFormat, firstDataLineNumber, vcfName, popFileName );
	
    // create cross referencing for population membership by sample:
    map<string, int> mapOfPopulations;      // key = population ID, value = integer population index
    map<string, int> mapOfSamples;          // key = sample ID, value = integer representing population index
    int numSamplesPerPopulation[numPopulations];    // for later frequency calculations
    makePopulationMap( mapOfPopulations, numPopulations, popFileName );
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
	VCFfile.close();
	PopulationFile.close();
    outputFile.close();
	// free memory:
	//delete mySamples;
	
	
#ifdef DEBUG
		cout << "\nI ran!!\n\n";
#endif
	
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


bool assignSamplesToPopulations(ifstream& VCFfile, int numSamples, int numFields, map<string, int> mapOfSamples, int *populationReference, unsigned long int& VCFfileLineCount, int firstDataLineNumber )
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
                    VCFfile.ignore(); // don't yet read in first entry of next line, but set the stage to
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


void calculateSummaryStats( stringstream& lineStream, ofstream& outputFile, int numTokensInFormat, int GTtoken, int DPtoken, int GQtoken, bool lookForDP, bool lookForGQ, string formatDelim, int formatOpsOrder[], int numSamples, int numPopulations, unsigned long int VCFfileLineCount, int* populationReference )
{
    int homoRefCount = 0, homoAltCount = 0, hetCount = 0, altAlleleCounts[numPopulations];
    int validSampleCounts[numPopulations], DPvalues[numSamples], GQvalues[numSamples];
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
    
    // loop over all columns of data:
    int counter = 0, operationCode, popIndex;
    string currentSample, token;
    size_t pos, strStart[numTokensInFormat], strLen[numTokensInFormat];
    string checkGTsep = "/";
    string allele1, allele2;
	int DPnoCall = 0, GQnoCall = 0;
    while ( lineStream >> currentSample ) {
        // currentSample is a string with format given by FORMAT
        
        // add a delimiter for easier parsing of final subfield in loop:
        currentSample = currentSample + formatDelim;
        popIndex = populationReference[ counter ];
        pos = 0;
        for ( int i = 0; i < numTokensInFormat; i++ ) {
            strStart[i] = pos;
            pos = currentSample.find( formatDelim, pos );  // next occurrence of delimiter
            strLen[i] = pos - strStart[i]; // length of this field
            pos++; // increment for next
        }
        
        // parse the current sample:
        for ( int i = 0; i < numTokensInFormat; i++ ) {
            // get the substring:
            //pos = currentSample.find( formatDelim ); // old way replaced by for loop above
            token = currentSample.substr(strStart[i], strLen[i]);
            // get operation code:
            operationCode = formatOpsOrder[i];
            if ( operationCode == GT_OPS_CODE ) {
                // parse the genotype data and add to correct population
                allele1 = token.substr(0,1);
                allele2 = token.substr(2,1);
                
                // considering the diploid genotype, there are 9 options:
                if ( allele1 == "0" ) {
                    if ( allele2 == "0" ) {
                        homoRefCount++;
                        validSampleCounts[popIndex] += 2; // diploid; no alt alleles
                    } else if ( allele2 == "1" ) {
                        hetCount++;
                        validSampleCounts[popIndex] += 2; // diploid
                        altAlleleCounts[popIndex]++; // one alt allele
                    } else {
#ifdef DEBUG
                        cout << "\nWarning!!\n\t: Allele2 not called in VCF line " << VCFfileLineCount;
                        cout << ", sample number " << (counter + 1);
#endif
                        validSampleCounts[popIndex]++;
                    }
                } else if ( allele1 == "1" ) {
                    if ( allele2 == "0" ) {
                        hetCount++;
                        validSampleCounts[popIndex] += 2; // diploid
                        altAlleleCounts[popIndex]++; // one alt allele
                    } else if ( allele2 == "1" ) {
                        homoAltCount++;
                        validSampleCounts[popIndex] += 2; // diploid
                        altAlleleCounts[popIndex] += 2; // two alt alleles
                    } else {
#ifdef DEBUG
                        cout << "\nWarning!!\n\t: Allele2 not called in VCF line " << VCFfileLineCount;
                        cout << ", sample number " << (counter + 1);
#endif
                        validSampleCounts[popIndex]++;
                        altAlleleCounts[popIndex]++; // one alt allele
                    }
                } else {
#ifdef DEBUG
                    cout << "\nWarning!!\n\t: Allele1 not called in VCF line " << VCFfileLineCount;
                    cout << ", sample number " << (counter + 1);
#endif
                    if ( allele2 == "0" || allele2 == "1" ) {
                        validSampleCounts[popIndex] += 1; // diploid
                        if ( allele2 == "1" )
                            altAlleleCounts[popIndex]++; // one alt allele
                    } else {
#ifdef DEBUG
                        cout << "\nWarning!!\n\t: Allele2 not called in VCF line " << VCFfileLineCount;
                        cout << ", sample number " << (counter + 1);
#endif
                    }
                }
                
                // now recording allele counts:
                
//#ifdef DEBUG
                if ( token.length() != 3 ) {
                    cout << "\nError in calculateSummaryStats()\n\t:GT token (" << token << ") has length != 3\n\t";
                    cout << "Aborting ... \n\n";
                }
                if ( token.substr(1,1) != checkGTsep ) {
                    cout << "\nError in calculateSummaryStats():\n\tGT token (" << token << ")";
                    cout << "does not have expected character (" << checkGTsep << ") between alleles.\n\t";
                    cout << "I found: " << token.substr(1,1) << endl;
                    cout << "Aborting ... \n\n";
                }
//                if ( counter % 100 == 0 ) {
//                    cout << "\nsample " << (counter+1) << ", loop count " << i << ", GT is " << token << endl;
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
				if ( token == "." ) {
					DPvalues[counter] = -1;
					DPnoCall++;
				} else {
                	DPvalues[counter] = stoi(token);
				}
                
//                cout << "\nsample " << (counter+1) << ", loop count " << i << ", DP is " << token << endl;
            } else if ( operationCode == GQ_OPS_CODE && lookForGQ ) {
                // add the GQ data to the GQ array
				if ( token == "." ) {
					GQvalues[counter] = -1;
					GQnoCall++;
				} else {
					GQvalues[counter] = stoi(token);
				}
//                cout << "\nsample " << (counter+1) << ", loop count " << i << ", GQ is " << token << endl;
            }
            
            // otherwise just skip it
            
            // delete the substring
            //currentSample.erase(0, pos + formatDelim.length()); // old way replaced by prior for loop
        }  // end of loop over tokens in sample
        
        counter++;  // keep track of how many samples have been processed
    }  // end of while() loop over lineStream
    
    // error checking:
    if ( counter != numSamples ) {
        cout << "\nError in calculateSummaryStats():\n\tlineStream did not give numSamples number of loops.\n\t";
        cout << "counter = " << counter << ", but numSamples = " << numSamples;
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
        freq = static_cast<double>( altAlleleCounts[i] ) / static_cast<double>( validSampleCounts[i] );
        outputFile << "\t" << freq << "\t" << validSampleCounts[i];
    }
    
    // outputFile << endl;  not needed here; this is done in parseActualData()
    
    
    
}


inline void checkFormatToken( string token, int& GTtoken, int& DPtoken, int& GQtoken, int subfieldCount  )
{
    // record sub-field:
    if ( token == "GT" )
        GTtoken = subfieldCount;
    else if ( token == "DP" )
        DPtoken = subfieldCount;
    else if ( token == "GQ" )
        GQtoken = subfieldCount;
}


void determineFormatOpsOrder( int numTokensInFormat, int GTtoken, int DPtoken, int GQtoken, bool lookForDP, bool lookForGQ, string formatDelim, int formatOpsOrder[], int maxSubfieldsInFormat )
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


inline void errorCheckTokens( int GTtoken, int DPtoken, int GQtoken, bool& lookForDP, bool& lookForGQ )
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
}



void makePopulationMap( map<string, int>& mapOfPopulations, int numPopulations, string popFileName )
{
    // helpful code for working with maps borrowed from:
    // https://thispointer.com/stdmap-tutorial-part-1-usage-detail-with-examples/
    
    string uniquePopFileName = popFileName + "_UniquePopFile.txt";
    ifstream uniquePopFile( uniquePopFileName );
    string popID;
    int index = 0;
    
    // get unique populations
    while( uniquePopFile >> popID ) {
        mapOfPopulations[ popID ] = index++;
    }
    
    uniquePopFile.close();
    
    if ( index != numPopulations ) {
        cout << "\nError!  numPopulations (" << numPopulations << ") != number of keys (" << index << ") in mapOfPopulations!\n\tExiting!\n\n";
        exit(-1);
    }
}


void parseActualData(ifstream& VCFfile, int numFormats, string formatDelim, int maxSubfieldsInFormat, unsigned long int& VCFfileLineCount, ofstream& outputFile, int numSamples, int numPopulations, int* populationReference, string vcfName )
{
    string oneLine, CHROM, POS, ID, REF, ALT, QUAL;
    long int dumCol, SNPcount = 0;
    bool keepThis, checkFormat = true, lookForDP, lookForGQ;
    int numTokensInFormat, GTtoken = -1, DPtoken = -1, GQtoken = -1;
	string discardedLinesFileName = vcfName + "_discardedLineNums.txt";
	ofstream discardedLinesFile( discardedLinesFileName, ostream::out );
    // the latter ints are for parsing GT = genotype, DP = depth,
    // and GQ = quality sub-fields of the FORMAT column
    int formatOpsOrder[maxSubfieldsInFormat]; // for keeping track of how to parse FORMAT efficiently
	
	discardedLinesFile << "VCFfileLinesNotUsed" << endl; // header row
    // work line by line:
    stringstream lineStream( "", ios_base::in | ios_base::out );
    while ( getline( VCFfile, oneLine ) ) {
        SNPcount++; // counter of how many SNP lines have been processed
        VCFfileLineCount++; // counter of how many LINES of VCF file have been processed
        
        // turn each line into a string stream for easy parsing by whitespace:
        lineStream.clear();
        lineStream.str( oneLine );
        
        // work with meta-col data:
        keepThis = parseMetaColData( lineStream, SNPcount, checkFormat, numTokensInFormat, GTtoken, DPtoken, GQtoken, lookForDP, lookForGQ, formatDelim, CHROM, POS, ID, REF, ALT, QUAL);
        
        if ( checkFormat ) {
            determineFormatOpsOrder( numTokensInFormat, GTtoken, DPtoken, GQtoken, lookForDP, lookForGQ, formatDelim, formatOpsOrder, maxSubfieldsInFormat );
        }
        
        if ( keepThis ) {
            // it is a biallelic SNP
            // print out meta fields:
            outputFile << VCFfileLineCount << "\t" << CHROM << "\t" << POS << "\t" << ID << "\t" << REF << "\t" << ALT << "\t" << QUAL;
            
            // let's calculate and store data using lineStream:
            calculateSummaryStats( lineStream, outputFile, numTokensInFormat, GTtoken, DPtoken, GQtoken, lookForDP, lookForGQ, formatDelim, formatOpsOrder, numSamples, numPopulations, VCFfileLineCount, populationReference );
            
            // add end of line (done with this line):
            outputFile << endl;
		} else {
			discardedLinesFile << VCFfileLineCount << endl;
		}
        
        if ( numFormats == 1 ) {
            checkFormat = false; // not needed after first SNP
        }
//        if ( SNPcount == 2 )
//            exit(0);
    }
	
	discardedLinesFile.close();
}


void parseCommandLineInput(int argc, char *argv[], ifstream& VCFfile, ifstream& PopulationFile, unsigned long int& maxCharPerLine, bool& popFileHeader, int& numSamples, int& numPopulations, int& numFields, int& numFormats, string& formatDelim, int& maxSubfieldsInFormat, int& firstDataLineNumber, string& vcfName, string& popFileName )
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

	// parse command line options:
	int flag;
    while ((flag = getopt(argc, argv, "V:P:L:H:N:n:F:f:D:S:l:")) != -1) {
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
            case 'D':
                formatDelim = optarg;
                break;
            case 'S':
                maxSubfieldsInFormat = atoi(optarg);
                break;
            case 'l':
                firstDataLineNumber = atoi(optarg);
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
#ifdef DEBUG
        cout << "\nMax length of string on this system = " << testString.max_size() << "\n\n";
#endif
	if ( !maxCharSet ) {
		maxCharPerLine = static_cast<unsigned long int>(testString.max_size());
		//cout << "\nWarning: maxCharPerLine (maximum line length) is not set.  Using default of " << maxCharPerLine << endl;
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


bool parseMetaColData( stringstream& lineStream, long int SNPcount, bool checkFormat, int& numTokensInFormat, int& GTtoken, int& DPtoken, int& GQtoken, bool& lookForDP, bool& lookForGQ, string formatDelim, string& CHROM, string& POS, string& ID, string& REF, string& ALT, string& QUAL )
{
    int subfieldCount;  // field counter, starting with index of 1
    string buffer, myDelim = formatDelim, token;
    string FILTER, INFO, FORMAT;
    size_t pos;
    bool keepThis;
    
    
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
                cout << "\nError in parseActualData():\n\t";
                cout << "case-switch not working as expected.\n\tAborting ...\n\n";
                exit(-3);
                break;
        }
    }
    
    // status report:
    if ( SNPcount % 5000 == 0 ) {
        cout << "\nSNP lines processed so far = " << SNPcount << "; Current SNP is:\n\t";
        cout <<  CHROM << "\t" << POS << "\t" << ID << "\t" << REF << "\t" << ALT << endl;
    }
    
    if ( checkFormat ) {
        subfieldCount = 0;
        // string splitting, with thanks to:
        // https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
#ifdef DEBUG
            cout << "** Parsing FORMAT **\n\tsubfieldcount\tpos\ttoken\n";
#endif
        while ( (pos = FORMAT.find( myDelim )) != string::npos ) {
            subfieldCount++;
            token = FORMAT.substr(0, pos);
            // record sub-field if it matches one we're looking for:
            checkFormatToken( token, GTtoken, DPtoken, GQtoken, subfieldCount );
#ifdef DEBUG
                cout << "\t" << subfieldCount << "\t\t" << pos  << "\t" << token << endl;
#endif
            FORMAT.erase(0, pos + myDelim.length());
        }
        numTokensInFormat = subfieldCount + 1;
        checkFormatToken( FORMAT, GTtoken, DPtoken, GQtoken, numTokensInFormat ); // check remaining field (after last occurrence of delimiter)
        // error checking on occurrence of GT:
        errorCheckTokens( GTtoken, DPtoken, GQtoken, lookForDP, lookForGQ );
        
#ifdef DEBUG
            cout << "\t" << numTokensInFormat << "\t\tlast\t" << FORMAT << "\twhile loop count = " << subfieldCount << endl;
            cout << "GTtoken = " << GTtoken << "; DPtoken = " << DPtoken << "; GQtoken = " << GQtoken << endl;
#endif
        
        
    }
    // check for bi-allelic SNPs:
    if ( REF == "N" || ALT == "N" || ALT.length() != 1 ) {
        keepThis = false;
    } else {
        keepThis = true;
    }
    
    if ( !keepThis ) {
#ifdef DEBUG
        cout << "\nSNP #" << SNPcount << ", ID = " << ID << ", has REF = " << REF << " and ALT = " << ALT << endl;
#endif
    }
    
    return keepThis;
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
