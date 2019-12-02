#include <iostream>
#include <string>
#include <sstream>
#include "IPEsimulation.hpp"

std::string getFlagWithOption(int argc, char* argv[], const std::string& flag);
bool getFlag(int argc, char* argv[], const std::string& flag);

int main( int argc, char *argv[] ) {
    std::stringstream helpMessage;
    helpMessage << "\n - - - - - SIMULATION OF INVERSE PATCHY COLLOIDS - - - - -\n\n"
                << "USAGE:\nYou need to specify a usage mode as the first, mandatory, argument. You have three options:\n"
                << " * \"new\": start a new simulation; you need the files input.in and stages.in;\n"
                << " * \"old\": resume an old simulation; you need the files input.in, stages.in, and a startingstate.xyz;\n"
                << "LA CACCA!\n";

   // parse number of arguments
   if(argc < 2 || argc > 4) {
        std::cerr << helpMessage.str();
        return EXIT_FAILURE;
    }

   // parse arguments with value
   std::string usageMode = argv[1];
   bool validUsageMode = (usageMode == "new") || (usageMode == "old");
   if(!validUsageMode) {
       std::cerr << helpMessage.str();
       return EXIT_FAILURE;
   }

   // open staging.in file
   std::ifstream stagesFile("stages.in");
   if(stagesFile.fail()) {
       std::cerr << "File stages.in could not be opened. Aborting.\n";
       return EXIT_FAILURE;
   }

    SimulationStage currentStage;
    currentStage.inputRestoringPreviousSimulation = (usageMode == "old");
    int simulatedStages = 0;
    std::cout << std::boolalpha;

    // read each line of the stagesFile and run a simulation with its temperature and duration
    while(stagesFile >> currentStage.inputTemperature >> currentStage.deltaTrans >> currentStage.deltaRot
                     >> currentStage.inputPrintTrajectoryAndCorrelations >> currentStage.inputStageTotalDuration) {
        ++simulatedStages;
        std::cout << "Simulation stage " << simulatedStages << " starting; set temperature: " << currentStage.inputTemperature << std::endl;
        // run the simulation stage
        IPEsimulation simulation(currentStage);
        simulation.run();
        currentStage.inputRestoringPreviousSimulation = true;
        std::cout << "Simulation stage " << simulatedStages << " is finished." << std::endl;
        // copy files
        if (currentStage.inputPrintTrajectoryAndCorrelations) {
            std::stringstream command;
            command << "mv siml siml_stage-" << simulatedStages << "_T-" << currentStage.inputTemperature;
            if(system(command.str().c_str()) != 0) {
                std::cerr << "Could not move siml/* between simulation stages. Aborting.\n";
                return EXIT_FAILURE;
            }
            std::stringstream().swap(command);
            command << "cp startingstate.xyz siml_stage-" << simulatedStages << "_T-" << currentStage.inputTemperature << "/startingstate.xyz";
            if(system(command.str().c_str()) != 0) {
                std::cerr << "Could not copy startingstate.xyz between simulation stages. Aborting.\n";
                return EXIT_FAILURE;
            }
        }
    }

    return EXIT_SUCCESS;
}


std::string getFlagWithOption(int argc, char* argv[], const std::string& flag) {
    std::string flagOption;
    for( int i = 0; i < argc; ++i) {
        const std::string argument = argv[i];
        if(0 == argument.find(flag)) {
            const std::size_t found = argument.find_last_of(flag);
            flagOption = argument.substr(found + 1);
            return flagOption;
        }
    }
    return flagOption;
}

bool getFlag(int argc, char* argv[], const std::string& flag) {
    for( int i = 0; i < argc; ++i) {
        const std::string argument = argv[i];
        if(argument == flag)
            return true;
    }
    return false;
}
