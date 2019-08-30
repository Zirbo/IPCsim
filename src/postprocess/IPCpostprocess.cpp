#include "IPCpostprocess.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>


IPCpostprocess::IPCpostprocess(int inputNumberOfPatches, size_t const inputNumberOfSubSimulations, std::string const& directoryName) {
    // input checks
    if (inputNumberOfPatches != 1 && inputNumberOfPatches != 2) {
        std::cerr << "At the moment only 1 and 2 patches are supported.\n";
        exit(1);
    }
    numberOfPatches = inputNumberOfPatches;
    if (inputNumberOfSubSimulations <= 0) {
        std::cerr << "The number of subsimulations needs to be a positive integer.\n";
        exit(1);
    }
    numberOfSubSimulations = inputNumberOfSubSimulations;

    // clean up old data and recreate output directory
    if(system("rm -rf analysis") != 0) {
        std::cerr << "Unable to delete the old 'analysis/' directory with rm -rf. "
                  << "Most likely you have it open somewhere or some program is running in it.\n";
        exit(1);
    }
    if(system("mkdir analysis") != 0) {
        std::cerr << "Unable to create a new 'analysis/' directory. You'll never see this error message.\n";
        exit(1);
    }

    // open simulation output file and trajectory
    const std::string trajectoryFileName = directoryName + "/trajectory.xyz";
    trajectoryFile.open(trajectoryFileName.c_str());
    if(trajectoryFile.fail()) {
        std::cerr << "File " << trajectoryFileName << " could not be opened. Aborting.\n";
        exit(1);
    }

    readOutputFile(directoryName);

    // initialize the containers
    ipcCentersPreviousPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersCurrentPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcInitialOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCurrentOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersInitialVelocities.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersCurrentVelocities.resize(nIPCs, {0.0, 0.0, 0.0});

    displacementOfEachIPCs.resize(nIPCs, {0.0, 0.0, 0.0});
    orientationAutocorrelation.resize(subSimulationDuration, 0.);
    velocityAutocorrelation.resize(subSimulationDuration, 0.);

    // initialize output files
    autocorrelationsFile.open("analysis/autocorrelationsFile.out");
    autocorrelationsFile << std::scientific << std::setprecision(6);
    meanSquaredDisplFile.open("analysis/meanSquaredDisplacement.out");
    meanSquaredDisplFile << std::scientific << std::setprecision(6);
    //typicalOrientationsFile.open("analysis/typicalOrientations.out");
    //typicalOrientationsFile << std::scientific << std::setprecision(6);
}

void IPCpostprocess::run() {
    for (double subSymCounter = 0; subSymCounter < numberOfSubSimulations; ++subSymCounter) {
        for (double subSymSnapshot = 0; subSymSnapshot < subSimulationDuration; ++subSymSnapshot) {
            const size_t absoluteSnapshot = subSymCounter*subSimulationDuration + subSymSnapshot;
            readSnapshot(absoluteSnapshot);
            if (subSymSnapshot == 0)
                updateInitialOrientationAndVelocites();
            if (subSymCounter == 0 && subSymSnapshot == 0)
                updatePreviousPositions(); // for the initial computation of the MSD...
            computeMSD(absoluteSnapshot);
            computeAutocorrelations(subSymSnapshot);
            updatePreviousPositions();
        }
        std::cout << "Subsym " << subSymCounter << " of " << numberOfSubSimulations << " finished\n";
    }

    printAutocorrelations();

    trajectoryFile.close();
    autocorrelationsFile.close();
    meanSquaredDisplFile.close();
    //typicalOrientationsFile.close();
}


void IPCpostprocess::readSnapshot(const size_t snapshotNumber) {
    runConsistencyChecks(snapshotNumber);

    for(int i = 0; i < nIPCs; ++i) {
        char uselessChar;
        double useless;
        // IPC center
        trajectoryFile >> uselessChar
                >> ipcCentersCurrentPositions[i][0] >> ipcCentersCurrentPositions[i][1] >> ipcCentersCurrentPositions[i][2]
                >> ipcCentersCurrentVelocities[i][0] >> ipcCentersCurrentVelocities[i][1] >> ipcCentersCurrentVelocities[i][2];
        // first patch
        trajectoryFile >> uselessChar
                >> ipcCurrentOrientations[i][0] >> ipcCurrentOrientations[i][1] >> ipcCurrentOrientations[i][2]
                >> useless >> useless >> useless;
        // second patch
        if (numberOfPatches == 2) {
            trajectoryFile >> uselessChar;
            // ignore the rest of the line
            trajectoryFile.ignore(500, '\n');
        }
        for (int j: {0, 1, 2}) {
            // compute orientations
            ipcCurrentOrientations[i][j] -= ipcCentersCurrentPositions[i][j];
            relativePBC(ipcCurrentOrientations[i][j]);
            ipcCurrentOrientations[i][j] /= firstPatchEccentricity;
        }
    }
}

void IPCpostprocess::updateInitialOrientationAndVelocites() {
    // set initial state and previous state as the current
    for(int i = 0; i < nIPCs; ++i) {
        for (int j: {0, 1, 2}) {
            ipcInitialOrientations[i][j] = ipcCurrentOrientations[i][j];
            ipcCentersInitialVelocities[i][j] = ipcCentersCurrentVelocities[i][j];
        }
    }
}

void IPCpostprocess::updatePreviousPositions() {
    // set initial state and previous state as the current
    for(int i = 0; i < nIPCs; ++i) {
        for (int j: {0, 1, 2}) {
            ipcCentersPreviousPositions[i][j] = ipcCentersCurrentPositions[i][j];
        }
    }

}

void IPCpostprocess::runConsistencyChecks(size_t const snapshotNumber) {
    int nIPCsCheck;
    double simulationBoxSideCheck;
    double time;
    trajectoryFile >> nIPCsCheck >> simulationBoxSideCheck >> time;
    if ((1+numberOfPatches)*nIPCs != nIPCsCheck) {
        std::cerr << "Inconsistency in the number of particles.\n"
                  << (1+numberOfPatches)*nIPCs << " != " << nIPCsCheck << ".\n";
        exit(1);
    }
    if (simulationBoxSide != simulationBoxSideCheck) {
        std::cerr << "Inconsistency in the simulation box side.\n"
                  << simulationBoxSide << " != " << simulationBoxSideCheck << ".\n";
        exit(1);
    }
    if (snapshotNumber*printingInterval != time) {
        std::cerr << "Inconsistency in the time.\n"
                  << snapshotNumber*printingInterval << " != " << time << ".\n";
        exit(1);
    }
}

void IPCpostprocess::readOutputFile(const std::string &directoryName) {
    // open file
    const std::string outputFileName = directoryName + "/output.out";
    std::ifstream simulationOutputFile(outputFileName.c_str());
    if(simulationOutputFile.fail()) {
        std::cerr << "File " << outputFileName << " could not be opened. Aborting.\n";
        exit(1);
    }

    // read simulation parameters from the simulation output file
    simulationOutputFile.ignore(200, '\n');
    double useless;
    simulationOutputFile >> nIPCs >> density >> temperature
            >> simulationTimeStep >> printingInterval >> simulationTotalDuration
            >> useless >> useless >> useless
            >> useless >> useless >> useless
            >> useless
            >> firstPatchEccentricity >> firstPatchRadius
            >> secndPatchEccentricity >> secndPatchRadius;

    //close file
    simulationOutputFile.close();

    // process data
    simulationBoxSide = std::cbrt(nIPCs/density);
    simulationDurationInIterations = (size_t)simulationTotalDuration/printingInterval;
    subSimulationDuration = (size_t)simulationDurationInIterations/numberOfSubSimulations;

    // scale lenghts
    firstPatchEccentricity /= simulationBoxSide;
}

void IPCpostprocess::computeMSD(size_t const snapshotNumber) {
    double meanSquaredDisplacementAtThisInstant = 0.;
    for(int i = 0; i < nIPCs; ++i) {
        for (int j: {0, 1, 2}) {
            double delta_xj = ipcCentersCurrentPositions[i][j] - ipcCentersPreviousPositions[i][j];
            relativePBC(delta_xj);
            displacementOfEachIPCs[i][j] += delta_xj;
            meanSquaredDisplacementAtThisInstant += std::pow(displacementOfEachIPCs[i][j],2);
        }
    }
    meanSquaredDisplFile << snapshotNumber*printingInterval << "\t" << meanSquaredDisplacementAtThisInstant << "\n";
}

void IPCpostprocess::computeAutocorrelations(size_t const snapshotNumber) {
    for(int i = 0; i < nIPCs; ++i) {
        for (int j: {0, 1, 2}) {
            // compute autocorrelations
            orientationAutocorrelation[snapshotNumber] += ipcCurrentOrientations[i][j]*ipcInitialOrientations[i][j];
            velocityAutocorrelation[snapshotNumber] += ipcCentersCurrentVelocities[i][j]*ipcCentersInitialVelocities[i][j];
        }
    }
}

void IPCpostprocess::printAutocorrelations() {
    const double orientationScalingFactor = 1./nIPCs/numberOfSubSimulations;
    const double velocityScalingFactor = 1./velocityAutocorrelation[0];
    for (size_t i = 0; i < subSimulationDuration; ++i) {
        autocorrelationsFile << i*printingInterval << "\t"
                             << orientationAutocorrelation[i]*orientationScalingFactor << "\t"
                             << velocityAutocorrelation[i]*velocityScalingFactor << "\n";
    }
}
