// header file of function prototypes and class definitons for VCFtoSummStats
#include <iostream>
#include <fstream>
#include <string>
#include <map>
using namespace std;

#include <boost/iostreams/filtering_streambuf.hpp>



// function prototypes (in alphabetical order):
void assignPopIndexToSamples( map<string, int>& mapOfPopulations, map<string, int>& mapOfSamples, ifstream& PopulationFile, int numSamplesPerPopulation[], int numPopulations, int numSamples );

bool assignSamplesToPopulations(istream& VCFfile, int numSamples, int numFields, map<string, int> mapOfSamples, int *populationReference, unsigned long int& VCFfileLineCount, int& firstDataLineNumber );

inline int calculateMedian( int values[], int n );

void calculateSummaryStats( istream& VCFfile, ofstream& outputFile, int numTokensInFormat, int GTtoken, int DPtoken, int GQtoken, bool lookForDP, bool lookForGQ, char formatDelim, int formatOpsOrder[], int numSamples, int numPopulations, unsigned long int VCFfileLineCount, int* populationReference  );

inline void checkFormatToken( char* token, int& GTtoken, int& DPtoken, int& GQtoken, int subfieldCount  );

void convertTimeInterval( clock_t myTimeInterval, int& minutes, double& seconds);

void createVCFfilter( boost::iostreams::filtering_streambuf<boost::iostreams::input>& myVCFin, string vcfName, ifstream& vcfUnfiltered );

void determineFormatOpsOrder( int numTokensInFormat, int GTtoken, int DPtoken, int GQtoken, bool lookForDP, bool lookForGQ, char formatDelim, int formatOpsOrder[], int maxSubfieldsInFormat );

inline void errorCheckTokens( int GTtoken, int DPtoken, int GQtoken, bool& lookForDP, bool& lookForGQ );

//void makePopulationMap( map<string, int>& mapOfPopulations, int numPopulations, string popFileName );

double extractDPvalue( char* INFObuffer, bool& lookForDPinINFO );

size_t getLength( char *myCharArray );

void parseActualData(istream& VCFfile, int numFormats, char formatDelim, int maxSubfieldsInFormat, unsigned long int& VCFfileLineCount, ofstream& outputFile, int numSamples, int numPopulations, int* populationReference, string vcfName );

void parseCommandLineInput(int argc, char *argv[], ifstream& PopulationFile, bool& popFileHeader, int& numSamples, int& numPopulations, int& numFields, int& numFormats, char& formatDelim, int& maxSubfieldsInFormat, string& vcfName, string& popFileName, map<string, int>& mapOfPopulations );

bool parseMetaColData( istream& VCFfile, long int SNPcount, bool checkFormat, int& numTokensInFormat, int& GTtoken, int& DPtoken, int& GQtoken, bool& lookForDP, bool& lookForGQ, char formatDelim, char* CHROM, char* POS, char* ID, char* REF, char* ALT, char* QUAL );

void parsePopulationDesigFile( string fname, int& numSamples, int& numPopulations, map<string,int>& mapOfPopulations, bool popFileHeader );

void setUpOutputFile (ofstream& outputFile, string vcfName, int numPopulations, map<string, int> mapOfPopulations );
