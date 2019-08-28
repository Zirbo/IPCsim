#include "IPCpostprocess.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>


IPCpostprocess::IPCpostprocess(int inputNumberOfPatches, std::string const& directoryName) {
    if (inputNumberOfPatches < 1 || inputNumberOfPatches > 2) {
        std::cerr << "At the moment only 1 and 2 patches are supported.\n";
        exit(1);
    }
    numberOfPatches = inputNumberOfPatches;
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
    ipcCentersInitialPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersPreviousPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcInitialOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcPreviousOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersInitialVelocities.resize(nIPCs, {0.0, 0.0, 0.0});

    meanSquaredDisplacement.resize(simulationDurationInIterations, 0.);
    totalRotation.resize(simulationDurationInIterations, 0.);
    orientationAutocorrelation.resize(simulationDurationInIterations, 0.);
    velocityAutocorrelation.resize(simulationDurationInIterations, 0.);

    // initialize output files
    autocorrelationsFile.open("analysis/autocorrelationsFile.out");
    autocorrelationsFile << std::scientific << std::setprecision(4);
    meanSquaredDisplFile.open("analysis/meanSquaredDisplacement.out");
    meanSquaredDisplFile << std::scientific << std::setprecision(4);
    //typicalOrientationsFile.open("analysis/typicalOrientations.out");
    //typicalOrientationsFile << std::scientific << std::setprecision(4);
}

void IPCpostprocess::run() {

    trajectorySnapshot = 0;
    readFirstSnapshot();

    for (trajectorySnapshot = 1; trajectorySnapshot <= simulationDurationInIterations; ++trajectorySnapshot) {
        processSingleSnapshot();
        std::cout << trajectorySnapshot << " / " << simulationDurationInIterations << "\n";
    }

    printResults();

    trajectoryFile.close();
    autocorrelationsFile.close();
    meanSquaredDisplFile.close();
    //typicalOrientationsFile.close();
}

void IPCpostprocess::printResults() {
    // set values that were not computed
    orientationAutocorrelation[0] = 1.;
    velocityAutocorrelation[0] = 1.;
    meanSquaredDisplacement[0] = 0.;
    totalRotation[0] = 0.;

    double inverseOfNumberOfIPCs = 1./nIPCs;
    for (size_t i = 0; i <= simulationDurationInIterations; ++i) {
        autocorrelationsFile << i*printingInterval << "\t"
                             << orientationAutocorrelation[i]*inverseOfNumberOfIPCs << "\t"
                             << velocityAutocorrelation[i]*inverseOfNumberOfIPCs << "\n";
        meanSquaredDisplFile << i*printingInterval << "\t"
                             << orientationAutocorrelation[i]*inverseOfNumberOfIPCs << "\t"
                             << velocityAutocorrelation[i]*inverseOfNumberOfIPCs << "\n";
    }
}

void IPCpostprocess::readFirstSnapshot() {
    runConsistencyChecks();

    for(int i = 0; i < nIPCs; ++i) {
        char uselessChar;
        double useless;
        // IPC center
        trajectoryFile >> uselessChar
                >> ipcCentersInitialPositions[i][0] >> ipcCentersInitialPositions[i][1] >> ipcCentersInitialPositions[i][2]
                >> ipcCentersInitialVelocities[i][0] >> ipcCentersInitialVelocities[i][1] >> ipcCentersInitialVelocities[i][2];
        // first patch
        trajectoryFile >> uselessChar
                >> ipcInitialOrientations[i][0] >> ipcInitialOrientations[i][1] >> ipcInitialOrientations[i][2]
                >> useless >> useless >> useless;
        // second patch
        if (numberOfPatches == 2) {
            trajectoryFile >> uselessChar;
            // ignore the rest of the line
            trajectoryFile.ignore(500, '\n');
        }
        for (int j: {0, 1, 2}) {
            // compute orientations
            ipcInitialOrientations[i][j] -= ipcCentersInitialPositions[i][j];
            relativePBC(ipcInitialOrientations[i][j]);
            ipcInitialOrientations[i][j] /= firstPatchEccentricity;

            // set previous state as the current
            ipcPreviousOrientations[i][j] = ipcInitialOrientations[i][j];
            ipcCentersPreviousPositions[i][j] = ipcCentersInitialPositions[i][j];
        }
    }

}

void IPCpostprocess::processSingleSnapshot() {
    runConsistencyChecks();

    // start reading out IPCs
    for(int i = 0; i < nIPCs; ++i) {
        double ipcCenterPosition [3];
        double ipcCenterVelocity [3];
        double ipcOrientation [3];
        char uselessChar;
        double useless;
        // IPC center
        trajectoryFile >> uselessChar
                >> ipcCenterPosition[0] >> ipcCenterPosition[1] >> ipcCenterPosition[2]
                >> ipcCenterVelocity[0] >> ipcCenterVelocity[1] >> ipcCenterVelocity[2];
        // first patch
        trajectoryFile >> uselessChar
                >> ipcOrientation[0] >> ipcOrientation[1] >> ipcOrientation[2]
                >> useless >> useless >> useless;
        // second patch
        if (numberOfPatches == 2) {
            trajectoryFile >> uselessChar;
            trajectoryFile.ignore(500, '\n');
        }
        // compute orientation
        for (int j: {0, 1, 2}) {
            ipcOrientation[j] -= ipcCenterPosition[j];
            relativePBC(ipcOrientation[j]);
            ipcOrientation[j] /= firstPatchEccentricity;
        }

        for (int j: {0, 1, 2}) {
            // compute mean squared displacement and total rotation
            //meanSquaredDisp +=

            // compute autocorrelations
            orientationAutocorrelation[trajectorySnapshot] += ipcOrientation[j]*ipcInitialOrientations[i][j];
            velocityAutocorrelation[trajectorySnapshot] += ipcCenterVelocity[j]*ipcCentersInitialVelocities[i][j];
        }
    }
}

void IPCpostprocess::runConsistencyChecks() {
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
    if (trajectorySnapshot*printingInterval != time) {
        std::cerr << "Inconsistency in the time.\n"
                  << trajectorySnapshot*printingInterval << " != " << time << ".\n";
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

    // scale lenghts
    firstPatchEccentricity /= simulationBoxSide;
}
