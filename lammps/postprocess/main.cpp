#include <iostream>
#include <string>
#include <sstream>
#include "IPCpostprocess.hpp"


int main ( int argc, char *argv[] ) {

    std::stringstream helpMessage;
    helpMessage << "You need to specify a valid trajectory file, inputfile, and the directory were the potentials are stored!\n";

    std::string trajFile;
    std::string inputFile;
    std::string potDirName;
    if(argc == 4) {
        trajFile = argv[1];
        inputFile = argv[2];
        potDirName = argv[3];
    } else {
        std::cout << helpMessage.str();
        exit(1);
    }

    IPCpostprocess postProcess(trajFile, inputFile, potDirName);
    postProcess.run();
}
