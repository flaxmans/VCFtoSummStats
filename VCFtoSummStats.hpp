// header file of function prototypes and class definitons for VCFtoSummStats
#include <iostream>
#include <fstream>
#include <string>
#include <map>
using namespace std;


// function prototypes (in alphabetical order):
void assignPopIndexToSamples( map<string, int>& mapOfPopulations, map<string, int>& mapOfSamples, ifstream& PopulationFile, int numSamplesPerPopulation[], int numPopulations );

bool assignSamplesToPopulations(ifstream& VCFfile, int numSamples, int numFields, map<string, int> mapOfSamples, int *populationReference);

void calculateAlleleFrequencies(ifstream& VCFfile, int numFormats, string formatDelim );

inline void checkFormatToken( string token, int& GTtoken, int& DPtoken, int& GQtoken, int subfieldCount  );

void makePopulationMap( map<string, int>& mapOfPopulations, int numPopulations  );

void parseCommandLineInput(int argc, char *argv[], ifstream& VCFfile, ifstream& PopulationFile, long int& maxCharPerLine, bool& popFileHeader, int& numSamples, int& numPopulations, int& numFields, int& numFormats, string& formatDelim );

bool parseMetaColData( stringstream& lineStream, long int SNPcount, bool checkFormat, int& numTokensInFormat, int& GTtoken, int& DPtoken, int& GQtoken, string formatDelim );

