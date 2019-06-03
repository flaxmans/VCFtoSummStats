// header file of function prototypes and class definitons for VCFtoSummStats
#include <iostream>
#include <fstream>
#include <string>
#include <map>
using namespace std;


// function prototypes (in alphabetical order):
void assignPopIndexToSamples( map<string, int>& mapOfPopulations, map<string, int>& mapOfSamples, ifstream& PopulationFile, int numSamplesPerPopulation[], int numPopulations, int numSamples );

bool assignSamplesToPopulations(ifstream& VCFfile, int numSamples, int numFields, map<string, int> mapOfSamples, int *populationReference, unsigned long int& VCFfileLineCount, int firstDataLineNumber );

inline int calculateMedian( int values[], int n );

void calculateSummaryStats( stringstream& lineStream, ofstream& outputFile, int numTokensInFormat, int GTtoken, int DPtoken, int GQtoken, bool lookForDP, bool lookForGQ, string formatDelim, int formatOpsOrder[], int numSamples, int numPopulations, unsigned long int VCFfileLineCount, int* populationReference  );

inline void checkFormatToken( string token, int& GTtoken, int& DPtoken, int& GQtoken, int subfieldCount  );

void determineFormatOpsOrder( int numTokensInFormat, int GTtoken, int DPtoken, int GQtoken, bool lookForDP, bool lookForGQ, string formatDelim, int formatOpsOrder[], int maxSubfieldsInFormat );

inline void errorCheckTokens( int GTtoken, int DPtoken, int GQtoken, bool& lookForDP, bool& lookForGQ );

void makePopulationMap( map<string, int>& mapOfPopulations, int numPopulations  );

void parseActualData(ifstream& VCFfile, int numFormats, string formatDelim, int maxSubfieldsInFormat, unsigned long int& VCFfileLineCount, ofstream& outputFile, int numSamples, int numPopulations, int* populationReference, string vcfName );

void parseCommandLineInput(int argc, char *argv[], ifstream& VCFfile, ifstream& PopulationFile, unsigned long int& maxCharPerLine, bool& popFileHeader, int& numSamples, int& numPopulations, int& numFields, int& numFormats, string& formatDelim, int& maxSubfieldsInFormat, int& firstDataLineNumber, string& vcfName );

bool parseMetaColData( stringstream& lineStream, long int SNPcount, bool checkFormat, int& numTokensInFormat, int& GTtoken, int& DPtoken, int& GQtoken, bool& lookForDP, bool& lookForGQ, string formatDelim, string& CHROM, string& POS, string& ID, string& REF, string& ALT, string& QUAL );

void setUpOutputFile (ofstream& outputFile, string vcfName, int numPopulations, map<string, int> mapOfPopulations );
