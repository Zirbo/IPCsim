#include <iostream>
#include <string>
#include <sstream>
#include "IPCpostprocess.hpp"


int main ( int argc, char *argv[] ) {

    std::stringstream helpMessage;
    helpMessage << "You need to specify a valid trajectory file and an inputfile!\n";

    std::string trajFile;
    std::string inputFile;
    if(argc == 3) {
        trajFile = argv[1];
        inputFile = argv[2];
    } else {
        std::cout << helpMessage.str();
        exit(1);
    }

    IPCpostprocess postProcess(trajFile, inputFile);
    std::cout << "Starting postprocess.\n";
    postProcess.run();
    std::cout << "Postprocess finished.\n";
}
