#include <iostream>
#include "IPCsimulation.hpp"

int main ( int argc, char *argv[] )
{
    bool doWarmup;
    if(argc != 2)
    {    std::cerr<<"I only accept a single parameter, use either \"new\" or \"old\".\n";    exit(1);    }
    else if(std::string(argv[1])=="new")
        doWarmup = false;
    else if(std::string(argv[1])=="old")
        doWarmup = true;
    else
    {    std::cerr<<"Uncorrect parameter, use either \"new\" or \"old\".\n";    exit(1);    }

    IPCsimulation simulation;
    simulation.run(doWarmup);
}
