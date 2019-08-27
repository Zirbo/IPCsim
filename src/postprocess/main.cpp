#include <iostream>
#include <string>
#include "IPCpostprocess.hpp"


int main ( int argc, char *argv[] ) {

    std::string directoryName;
    if(argc == 1) {
        std::cout << "No directory specified, assuming it's called \"siml\".\n";
        directoryName = "siml";
    } else if(argc == 2) {
        directoryName = argv[1];
    } else if(argc > 2) {
        std::cerr << "Too many arguments :P\n";
        exit(1);
    }

    std::cout << "Starting postprocess...\n";
    IPCpostprocess postProcess(directoryName);
    postProcess.run();
    std::cout << "Postprocess finished.\n";
}
