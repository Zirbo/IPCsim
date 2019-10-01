#include <iostream>
#include <sstream>
#include <list>
#include "IPCsimulation.hpp"


int main ( int argc, char *argv[] ) {
    std::stringstream helpMessage;
    helpMessage << "USAGE:\nYou need to specify an IPC type and a usage mode.\n"
                << "IPC type must be either:\n * IPC\n * Janus.\n"
                << "For the usage mode there are 3 options:\n"
                << " * \"new\": start a new simulation; you need the files input.in and stages.in;\n"
                << " * \"old\": resume an old simulation; you need the files input.in, stages.in, and a startingstate.xyz;\n"
                << " * \"printpot\": print potentials in potentials.out, then exit.\n";

    // parse number of arguments
    if(argc != 3) {
        std::cerr << helpMessage.str();
        exit(1);
    }

    // parse IPC type
    bool janusSimulation;
    if(std::string(argv[2]) == "Janus") {
        janusSimulation = true;
    } else if(std::string(argv[2]) == "IPC") {
        janusSimulation = false;
    } else {
        std::cerr << helpMessage.str();
        exit(1);
    }

    // parse usage mode
    if(std::string(argv[2]) == "printpot") {
        SimulationStage emptyStage;
        emptyStage.janusSimulation = janusSimulation;
        IPCsimulation simulation(emptyStage);
        simulation.printPotentials();
    } else if (std::string(argv[2]) == "new" || std::string(argv[2]) == "old") {

        // open staging.in file
        std::ifstream stagesFile("stages.in");
        if(stagesFile.fail()) {
            std::cerr << "File stages.in could not be opened. Aborting.\n";
            exit(1);
        }

        SimulationStage currentStage;
        currentStage.janusSimulation = janusSimulation;
        currentStage.inputRestoringPreviousSimulation = (std::string(argv[2]) == "old");
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
                        exit(1);
                    }
                    std::stringstream().swap(command);
                    command << "cp startingstate.xyz siml_stage-" << simulatedStages << "_T-" << averageTemperatureinTheRun << "/startingstate.xyz";
                    if(system(command.str().c_str()) != 0) {
                        std::cerr << "Could not copy startingstate.xyz between simulation stages. Aborting.\n";
                        exit(1);
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
    }
    else {
        std::cerr << helpMessage.str();
        exit(1);
    }


}
