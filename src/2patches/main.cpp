#include <iostream>
#include <sstream>
#include <list>
#include "IPCsimulation.hpp"


int main ( int argc, char *argv[] ) {
    std::stringstream helpMessage;
    helpMessage << "USAGE\nYou need to specify a usage mode. You have four options:\n"
                << " * \"new\": start a new simulation, that will follow the parameters specified in input.in\n"
                << " * \"old\": resume an old simulation, that will follow the parameters specified in input.in,\n"
                << "   except the number of particles, which, together with the positions and velocities, will be read\n"
                << "   from a startingstate.xyz file; the velocities will be scaled so that the initial temperature will\n"
                << "   assume the value that you specify in the input.in --- or unchanged if this value is negative\n"
                << " * \"newstaged\": start a new simulation that will follow the parameters specified in input.in,\n"
                << "   expect temperature and simulation duration, and that will be restarted as many times as needed,\n"
                << "   following the temperatures and simulation durations defined in staging.in\n"
                << " * \"oldstaged\": start a new simulation that will follow the parameters specified in input.in,\n"
                << "   expect temperature and simulation duration, and that will be restarted as many times as needed,\n"
                << "   following the temperatures and simulation durations defined in staging.in\n\n";

    if(argc != 2) {
        std::cerr << helpMessage.str();
        exit(1);
    } else if(std::string(argv[1])=="new") {
        IPCsimulation simulation(false);
        simulation.run();
    } else if(std::string(argv[1])=="old") {
        IPCsimulation simulation(true);
        simulation.run();
    } else if (std::string(argv[1]) == "newstaged" || std::string(argv[1]) == "oldstaged") {

        // open staging.in file
        std::ifstream stagingFile("staging.in");
        if(stagingFile.fail()) {
            std::cerr << "File staging.in could not be opened. Aborting.";
            exit(1);
        }

        // simulate numberOfStages times :P
        std::pair<double,int> currentStage;
        bool resumingSimulation = (std::string(argv[1]) == "oldstaged");
        int counter = 0;
        double tollerance {0.}, temperature{0.};
        bool repeatStage{false}, recordStage{false};

        std::cout << std::boolalpha;

        //while(stagingFile >> currentStage.first && stagingFile >> tollerance && stagingFile >> currentStage.second) {
        while(stagingFile >> currentStage.first >> tollerance >> recordStage >> currentStage.second) {
            do {
                std::cout << "Simulation stage " << ++counter << " starting; set temperature: " << currentStage.first << std::endl;
                // run the simulation stage
                IPCsimulation simulation(resumingSimulation, true, currentStage);
                temperature = simulation.run();
                resumingSimulation = true;
                std::cout << "Simulation stage " << counter << " is finished; average temperature: " << temperature << std::endl;
                // copy files
                if (recordStage) {
                    std::stringstream command;
                    command << "mv siml siml_stage-" << counter << "_T-" << temperature;
                    if(system(command.str().c_str()) != 0) {
                        std::cerr << "Could not move siml/* between simulation stages. Aborting.";
                        exit(1);
                    }
                    std::stringstream().swap(command);
                    command << "cp startingstate.xyz siml_stage-" << counter << "_T-" << temperature << "/startingstate.xyz";
                    if(system(command.str().c_str()) != 0) {
                        std::cerr << "Could not copy startingstate.xyz between simulation stages. Aborting.";
                        exit(1);
                    }
                }

                if(tollerance > 0.) {
                    double relativeDifference = std::fabs(   (temperature-currentStage.first)/temperature   );
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
