// header file of function prototypes and class definitons for VCFtoSummStats
#include <iostream>
#include <fstream>
#include <string>
#include <map>
using namespace std;


// function prototypes (in alphabetical order):
void assignPopIndexToSamples( map<string, int>& mapOfPopulations, map<string, int>& mapOfSamples, ifstream& PopulationFile );

bool assignSamplesToPopulations(ifstream& VCFfile, int numSamples, int numFields, map<string, int> mapOfSamples, int *populationReference);

void makePopulationMap( map<string, int>& mapOfPopulations, int numPopulations  );

void parseCommandLineInput(int argc, char *argv[], ifstream& VCFfile, ifstream& PopulationFile, int& maxCharPerLine, bool& popFileHeader, int& numSamples, int& numPopulations, int& numFields);



