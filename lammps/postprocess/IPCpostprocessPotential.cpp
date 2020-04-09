#include <fstream>
#include <iostream>

#include "IPCpostprocessPotential.hpp"

void IPCpotential::initialize(std::string const& potentialsDir) {
    //
    std::string fileName = potentialsDir + "/BB.table";
    double spacingBB = 0.;
    readFile(fileName, spacingBB, uBB);

    fileName = potentialsDir + "/Bs1.table";
    double spacingBS = 0.;
    readFile(fileName, spacingBS, uBs);

    fileName = potentialsDir + "/Bs1.table";
    double spacingSS = 0.;
    readFile(fileName, spacingSS, uss);

    if( !(spacingBB == spacingBS) || !(spacingBS== spacingSS) ) {
        std::cerr << __func__ << ":: something really shitty is going on in the pair potential computation.";
        exit(1);
    }
    spacing = spacingBB;
}

void IPCpotential::readFile(std::string const& fileName, double & spacing, std::vector<double> & container) {
    std::ifstream potentialFile;
    potentialFile.open(fileName);
    if(potentialFile.fail()) {
        std::cerr << "File " << fileName << " could not be opened. Aborting.\n";
        exit(1);
    }
    potentialFile.ignore(200, '\n'); // skip the comment line
    std::string dummy;
    int size;
    double dummyD;
    potentialFile >> dummy >> dummy >> size;
    container.resize(size);

    potentialFile >> dummy >> spacing >> container[0] >> dummyD;

    for (int i = 1; i < size; ++i) {
        potentialFile >> dummy >> dummyD >> container[i] >> dummyD;
    }

}
