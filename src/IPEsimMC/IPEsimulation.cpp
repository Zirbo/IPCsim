#include <cstdlib>
#include <iomanip>
#include <iostream>
#include "IPEsimulation.hpp"

//************************************************************************//
IPEsimulation::IPEsimulation(SimulationStage const& stage) {
    // clean up old data and recreate output directory
    if(system("rm -rf siml") != 0) {
        std::cerr << "Unable to delete the old 'siml/' directory with rm -rf. "
                  << "Most likely you have it open somewhere or some program is running in it.\n";
        exit(1);
    }
    if(system("mkdir siml") != 0) {
        std::cerr << "Unable to create a new 'siml/' directory. You'll never see this error message.\n";
        exit(1);
    }

    // open output files
    outputFile.open("siml/output.out");
    energyTrajectoryFile.open("siml/evolution.out");
    energyTrajectoryFile << std::scientific << std::setprecision(6);
    energyTrajectoryFile << "#t\t\t\tT\t\t\tK\t\t\tU\t\t\tE\t\t\trmin\n";

    // initialize system
    initializeSystem(stage);

    // print starting configuration and initialize output file
    outputFile << "\nPlot evolution.out to check the evolution of the system.\n";
    outputSystemEnergies(energyTrajectoryFile);


    if (printTrajectoryAndCorrelations) {
        // initialize g(r)
        pairCorrelation.initialize(20, simulationBoxSide, nIPCs);

        // initialize trajectory output file
        trajectoryFile.open("siml/trajectory.xyz");
        trajectoryFile << std::scientific << std::setprecision(24);
        outputSystemTrajectory(trajectoryFile);
    }
}

//************************************************************************//
void IPEsimulation::run() {
    time_t simulationStartTime, simulationEndTime;

    double averagePotentialEnergy = 0.;
    int sampleSizeForAverages = 0;

    // simulation begins
    time(&simulationStartTime);
    while(simulationTime < simulationTotalDuration) {
        ++simulationTime;
        computeSimulationStep();

        if( simulationTime%printingInterval == 0) {
            // compute and output energies
            computeSystemEnergy();
            outputSystemEnergies(energyTrajectoryFile);
            // output trajectory and compute g(r)
            if (printTrajectoryAndCorrelations) {
                // g(r)
                pairCorrelation.compute(particles);
                outputSystemTrajectory(trajectoryFile);
            }
            // compute averages
            ++sampleSizeForAverages;
            averagePotentialEnergy += potentialEnergy;
        }
    }
    // simulation ends
    time(&simulationEndTime);
    outputFile << "The simulation lasted " << difftime (simulationEndTime, simulationStartTime) << " seconds.\n";

    // output final state
    std::ofstream finalStateFile("startingstate.xyz");
    finalStateFile << std::scientific << std::setprecision(24);
    outputSystemTrajectory(finalStateFile);
    finalStateFile.close();

    averagePotentialEnergy /= sampleSizeForAverages;
    outputFile << "Average potential energy during the simulation run = " << averagePotentialEnergy/nIPCs << std::endl;

    // output g(r);
    if (printTrajectoryAndCorrelations) {
        const double g_r_integral = pairCorrelation.print("siml/g_r");
        outputFile << "The integral of g(r) is " << g_r_integral << " and is should be equal to the number of particles minus one, " << nIPCs-1 << std::endl;
    }
}
