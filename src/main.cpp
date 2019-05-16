#include <iostream>
#include "IPCsimulation.hpp"

int main ( int argc, char *argv[] )
{
    bool restorePreviousSimulation;
    if(argc != 2) {
        std::cerr << "I only accept a single parameter, use either \"new\" or \"old\".\n";
        exit(1);
    }
    else if(std::string(argv[1])=="new")
        restorePreviousSimulation = false;
    else if(std::string(argv[1])=="old")
        restorePreviousSimulation = true;
    else {
        std::cerr << "Uncorrect parameter, use either \"new\" or \"old\".\n";
        exit(1);
    }

    IPCsimulation simulation(restorePreviousSimulation);
    simulation.run();
}
