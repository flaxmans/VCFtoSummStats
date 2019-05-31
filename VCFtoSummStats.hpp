// header file of function prototypes and class definitons for VCFtoSummStats
#include <iostream>
#include <fstream>
#include <string>
using namespace std;


// class definitions
class sample {
	private:
		string ID; 				// name or identifier
		string population;		// population to which it is assigned
	public:
		void setID( string name ) {
			ID = name;
		}
		void setPopulation( string name ) {
			population = name;
		}
		string getID( void ) {
			return ID;
		}
		string getPopulation( void ) {
			return population;
		}
};



// function prototypes:
void parseCommandLineInput(int argc, char *argv[], ifstream& VCFfile, ifstream& PopulationFile, int& maxCharPerLine, bool& popFileHeader, int& numSamples, int& numPopulations, int& numFields);

void assignSamplesToPopulations(sample *mySamples, ifstream& VCFfile, ifstream& PopulationFile, int numSamples, int numPopulations, int numFields);



