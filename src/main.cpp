#include <iostream>
#include <sstream>
#include <list>
#include "IPCsimulation.hpp"


int main ( int argc, char *argv[] ) {
    std::stringstream helpMessage;
    helpMessage << "USAGE\nYou need to specify a usage mode. You have three options:\n"
                << " * \"new\": start a new simulation, that will follow the parameters specified in input.in\n"
                << " * \"old\": resume an old simulation, that will follow the parameters specified in input.in, except the number of particles, which, together with the positions and velocities, will be read from a startingstate.xyz file; the velocities will be scaled so that the initial temperature will assume the value that you specify in the input.in --- or unchanged if this value is negative"
                << " * \"staged\": start a new simulation that will follow the parameters specified in input.in, then restarted as many times as needed ";

    if(argc != 2) {
        std::cerr << helpMessage.str();
        exit(1);
    } else if(std::string(argv[1])=="new") {
        IPCsimulation simulation(false);
        simulation.run();
    } else if(std::string(argv[1])=="old") {
        IPCsimulation simulation(true);
        simulation.run();
    } else if (std::string(argv[1])=="staged") {

        // read staging.in file
        std::ifstream stagingFile("staging.in");
        if(stagingFile.fail()) {
            std::cerr << "File staging.in could not be opened. Aborting.";
            exit(1);
        }
        std::list<std::pair<double,int>> simulationStages;
        std::pair<double,int> temp;
        while(stagingFile >> temp.first && stagingFile >> temp.second) {
            simulationStages.push_back(temp);
        }
        stagingFile.close();

        // simulate numberOfStages times :P
        bool isFirstSim = true;
        int counter = 0;
        for (auto stage: simulationStages) {
            std::cout << "Simulation stage " << ++counter << "/" << simulationStages.size() << " starting." << std::endl;
            // run the simulation stage
            IPCsimulation simulation(!isFirstSim, true, stage);
            simulation.run();
            isFirstSim = false;
            std::cout << "Simulation stage " << counter << "/" << simulationStages.size() << " is finished." << std::endl;
            // copy files
            std::stringstream command;
            command << "mv siml siml_stage-" << counter << "_T-" << stage.first;
            if(system(command.str().c_str()) != 0) {
                std::cerr << "Could not move siml/* between simulation stages. Aborting.";
                exit(1);
            }
            std::stringstream().swap(command);
            command << "cp startingstate.xyz siml_stage-" << counter << "_T-" << stage.first << "/startingstate.xyz";
            if(system(command.str().c_str()) != 0) {
                std::cerr << "Could not copy startingstate.xyz between simulation stages. Aborting.";
                exit(1);
            }
        }
    }
    else {
        std::cerr << "Uncorrect parameter, use either \"new\" or \"old\".\n";
        exit(1);
    }


}
