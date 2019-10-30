#include <iostream>
#include <string>
#include <cstring>
#include <sstream>
#include <list>
#include <algorithm>
#include "IPCsimulation.hpp"

std::string getFlagWithOption(int argc, char* argv[], const std::string& flag);
bool getFlag(int argc, char* argv[], const std::string& flag);

int main ( int argc, char *argv[] ) {
    std::stringstream helpMessage;
    helpMessage << "\n - - - - - SIMULATION OF INVERSE PATCHY COLLOIDS - - - - -\n\n"
                << "USAGE:\nYou need to specify a usage mode as the first, mandatory, argument. You have three options:\n"
                << " * \"new\": start a new simulation; you need the files input.in and stages.in;\n"
                << " * \"old\": resume an old simulation; you need the files input.in, stages.in, and a startingstate.xyz;\n"
                << " * \"printpot\": print potentials in lammps format in potentials.out, then exit.\n"
                << "\nOPTIONAL FLAGS:\n"
                << "   --janus        - Run a simulation with only one patch.\n"
                << "                    Information in input.in regarding the second patch will be ignored.\n"
                << "   --binary       - Run a system with a binary mixture of oppositely charged particles.\n"
                << "                    You must specify the percentage of unlike-charged particles like this:\n"
                << "                          --binary=P\n"
                << "                    where P is an integer number between 1 and 50.\n"
                << "  --printforces   - Prints also the forces in the trajectory.xyz file (only in the stages where it is created).\n"
                << "                    Be careful, the postprocess program does not support this format.\n"
                << "                    The startingstate.xyz file will not be affected."
                << "\n\nEXAMPLES:\n"
                << "./IPCsim printpot\n"
                << "./IPCsim new --binary=40\n"
                << "./IPCsim old --janus --printforces\n"
                << "\n\nIMPORTANT!!!\n"
                << "There is error checking for the flag options, but not for the flags themselves!\n"
                << "This means that, for example, something like\n"
                << "./IPCsim new --BINARY 40        or\n./IPCsim old --prntfrcs\nwould be silently ignored!\n";

    // parse number of arguments
   if(argc < 2 || argc > 4) {
        std::cerr << helpMessage.str();
        return EXIT_FAILURE;
    }

    // parse arguments with value
    std::string usageMode = argv[1];
    bool validUsageMode = (usageMode == "printpot") || (usageMode == "new") || (usageMode == "old");
    if(!validUsageMode) {
        std::cerr << helpMessage.str();
        return EXIT_FAILURE;
    }
    std::string binaryMixtureEnabled = getFlagWithOption(argc, argv, "--binary=");
    int binaryMixturePercentage = -1;
    if(!binaryMixtureEnabled.empty()) {
        binaryMixturePercentage = std::stoi(binaryMixtureEnabled);
        if(binaryMixturePercentage < 1 || binaryMixturePercentage > 50) {
            std::cerr << helpMessage.str();
            return EXIT_FAILURE;
        }
    }
    // boolean flags
    bool janusSimulation = getFlag(argc, argv, "--janus") || getFlag(argc, argv, "-j");
    bool printForces = getFlag(argc, argv, "--printforces") || getFlag(argc, argv, "-pf");


    // if we have to print the potentials, do that and not think too much
    if(usageMode == "printpot") {
        SimulationStage emptyStage;
        emptyStage.janusSimulation = janusSimulation;
        IPCsimulation simulation(emptyStage);
        simulation.printPotentials();
        return EXIT_SUCCESS;
    }

    // open staging.in file
    std::ifstream stagesFile("stages.in");
    if(stagesFile.fail()) {
        std::cerr << "File stages.in could not be opened. Aborting.\n";
        return EXIT_FAILURE;
    }

    SimulationStage currentStage;
    currentStage.janusSimulation = janusSimulation;
    currentStage.inputRestoringPreviousSimulation = (usageMode == "old");
    currentStage.binaryMixturePercentage = binaryMixturePercentage;
    currentStage.printForces = printForces;
    int simulatedStages = 0;
    double tollerance = 0.;
    std::cout << std::boolalpha;

    // read each line of the stagesFile and run a simulation with its temperature and duration
    while(stagesFile >> currentStage.inputStartingTemperature >> tollerance
                     >> currentStage.inputPrintTrajectoryAndCorrelations >> currentStage.inputStageTotalDuration) {
        // if a tollerance is given, repeat each stage until the average temperature of the run is inside the tollerance.
        bool repeatStage{false};
        do {
            ++simulatedStages;
            std::cout << "Simulation stage " << simulatedStages << " starting; set temperature: " << currentStage.inputStartingTemperature << std::endl;
            // run the simulation stage
            IPCsimulation simulation(currentStage);
            const double averageTemperatureinTheRun = simulation.run();
            currentStage.inputRestoringPreviousSimulation = true;
            std::cout << "Simulation stage " << simulatedStages << " is finished; average temperature: " << averageTemperatureinTheRun << std::endl;
            // copy files
            if (currentStage.inputPrintTrajectoryAndCorrelations) {
                std::stringstream command;
                command << "mv siml siml_stage-" << simulatedStages << "_T-" << averageTemperatureinTheRun;
                if(system(command.str().c_str()) != 0) {
                    std::cerr << "Could not move siml/* between simulation stages. Aborting.\n";
                    return EXIT_FAILURE;
                }
                std::stringstream().swap(command);
                command << "cp startingstate.xyz siml_stage-" << simulatedStages << "_T-" << averageTemperatureinTheRun << "/startingstate.xyz";
                if(system(command.str().c_str()) != 0) {
                    std::cerr << "Could not copy startingstate.xyz between simulation stages. Aborting.\n";
                    return EXIT_FAILURE;
                }
            }

            if(tollerance > 0.) {
                double relativeDifference = std::fabs(   (averageTemperatureinTheRun-currentStage.inputStartingTemperature)/averageTemperatureinTheRun   );
                repeatStage = relativeDifference > tollerance;
                std::cout << "Simulation stage " << simulatedStages << ", relative temperature difference: " << relativeDifference << std::endl;
                std::cout << "Repeat? " << repeatStage << std::endl;
            }
        } while(repeatStage);
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

