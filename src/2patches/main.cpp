#include <iostream>
#include <sstream>
#include <list>
#include "IPCsimulation.hpp"


int main ( int argc, char *argv[] ) {
    std::stringstream helpMessage;
    helpMessage << "USAGE\nYou need to specify a usage mode. You have five options:\n"
                << " * \"new\": start a new simulation; you need the files input.in and staging.in;\n"
                << " * \"old\": resume an old simulation; you need the files input.in, staging.in, and a startingstate.xyz;\n"
                << " * \"printpot\": print potentials in potentials.out, then exit.\n";

    if(argc != 2) {
        std::cerr << helpMessage.str();
        exit(1);
    } else if(std::string(argv[1]) == "printpot") {
        const SimulationStage stage;
        IPCsimulation simulation(stage);
        simulation.printPotentials();
    } else if (std::string(argv[1]) == "new" || std::string(argv[1]) == "old") {

        // open staging.in file
        std::ifstream stagingFile("staging.in");
        if(stagingFile.fail()) {
            std::cerr << "File staging.in could not be opened. Aborting.";
            exit(1);
        }

        // simulate numberOfStages times :P
        SimulationStage currentStage;
        currentStage.inputRestoringPreviousSimulation = (std::string(argv[1]) == "old");
        int counter = 0;
        double tollerance {0.}, averageTemperature{0.};
        bool repeatStage{false};

        std::cout << std::boolalpha;

        while(stagingFile >> currentStage.inputStartingTemperature >> tollerance >> currentStage.inputPrintTrajectoryAndCorrelations >> currentStage.inputStageTotalDuration) {
            do {
                std::cout << "Simulation stage " << ++counter << " starting; set temperature: " << currentStage.inputStartingTemperature << std::endl;
                // run the simulation stage
                IPCsimulation simulation(currentStage);
                averageTemperature = simulation.run();
                currentStage.inputRestoringPreviousSimulation = true;
                std::cout << "Simulation stage " << counter << " is finished; average temperature: " << averageTemperature << std::endl;
                // copy files
                if (currentStage.inputPrintTrajectoryAndCorrelations) {
                    std::stringstream command;
                    command << "mv siml siml_stage-" << counter << "_T-" << averageTemperature;
                    if(system(command.str().c_str()) != 0) {
                        std::cerr << "Could not move siml/* between simulation stages. Aborting.";
                        exit(1);
                    }
                    std::stringstream().swap(command);
                    command << "cp startingstate.xyz siml_stage-" << counter << "_T-" << averageTemperature << "/startingstate.xyz";
                    if(system(command.str().c_str()) != 0) {
                        std::cerr << "Could not copy startingstate.xyz between simulation stages. Aborting.";
                        exit(1);
                    }
                }

                if(tollerance > 0.) {
                    double relativeDifference = std::fabs(   (averageTemperature-currentStage.inputStartingTemperature)/averageTemperature   );
                    repeatStage = relativeDifference > tollerance;
                    std::cout << "Simulation stage " << counter << ", relative temperature difference: " << relativeDifference << std::endl;
                    std::cout << "Repeat? " << repeatStage << std::endl;
                }
            } while(repeatStage);
        }
    }
    else {
        std::cerr << helpMessage.str();
        exit(1);
    }


}
