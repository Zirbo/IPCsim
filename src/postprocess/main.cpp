#include <iostream>
#include <string>
#include <sstream>
#include "IPCpostprocess.hpp"


int main ( int argc, char *argv[] ) {

    std::stringstream helpMessage;
    helpMessage << "USAGE\nYou need to specify two arguments:\n"
                << " * the first is the number of patches in the files to be read\n"
                << " * the second is the directory where output.out and trajectory.xyz are\n"
                << "   (if you don't specify it, it will be assumed to be \"siml\"!)\n";

    std::string directoryName;
    int numberOfPatches;
    if(argc == 2) {
        numberOfPatches = std::stoi(argv[1]);
        directoryName = "siml";
    } else if(argc == 3) {
        numberOfPatches = std::stoi(argv[1]);
        directoryName = argv[2];
    } else {
        std::cout << helpMessage.str();
        exit(1);
    }

    std::cout << "Initializing postprocess with " << numberOfPatches << " patches, "
              << "the files are expected in the directory " << directoryName << ".\n";
    IPCpostprocess postProcess(numberOfPatches, directoryName);
    std::cout << "Starting postprocess.\n";
    postProcess.run();
    std::cout << "Postprocess finished.\n";
}
