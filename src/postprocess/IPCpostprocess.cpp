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
    const std::string outputFileName = directoryName + "/output.out";
    std::ifstream simulationOutputFile(outputFileName.c_str());
    if(simulationOutputFile.fail()) {
        std::cerr << "File " << outputFileName << " could not be opened. Aborting.\n";
        exit(1);
    }
    const std::string trajectoryFileName = directoryName + "/trajectory.xyz";
    trajectoryFile.open(trajectoryFileName.c_str());
    if(trajectoryFile.fail()) {
        std::cerr << "File " << trajectoryFileName << " could not be opened. Aborting.\n";
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
    simulationBoxSide = std::cbrt(nIPCs/density);
    simulationOutputFile.close();

    // double check, remove later!
    std::cout << nIPCs << "\t" << density << "\t" << temperature << "\n"
            << simulationTimeStep << "\t" << printingInterval << "\t" << simulationTotalDuration << "\n"
            << firstPatchEccentricity << "\t" << firstPatchRadius << "\n"
            << secndPatchEccentricity << "\t" << secndPatchRadius << "\n";
    // scale lenghts
    firstPatchEccentricity /= simulationBoxSide;

    // initialize the containers
    ipcCentersInitialPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersPreviousPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcInitialOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcPreviousOrientations.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersInitialVelocities.resize(nIPCs, {0.0, 0.0, 0.0});

    readFirstConfiguration();

    autocorrelationsFile.open("analysis/autocorrelationsFile.out");
    autocorrelationsFile << std::scientific << std::setprecision(4);
    meanSquaredDisplFile.open("analysis/meanSquaredDisplacement.out");
    meanSquaredDisplFile << std::scientific << std::setprecision(4);
    typicalOrientationsFile.open("analysis/typicalOrientations.out");
    typicalOrientationsFile << std::scientific << std::setprecision(4);
}

void IPCpostprocess::readFirstConfiguration() {
    // read and double check number of particles and simulation box side
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
    if (0. != time) {
        std::cerr << "Inconsistency in the time.\n"
                  << 0. << " != " << time << ".\n";
        exit(1);
    }

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

void IPCpostprocess::run() {
    const size_t simulationDurationInIterations = (size_t)simulationTotalDuration/printingInterval;

    for (simulationTime = 0; simulationTime <= simulationDurationInIterations; ++simulationTime) {
        std::cout << simulationTime << " / " << simulationDurationInIterations << "\n";
        processSingleConfiguration();
    }

    trajectoryFile.close();
    autocorrelationsFile.close();
    meanSquaredDisplFile.close();
    typicalOrientationsFile.close();
}

void IPCpostprocess::processSingleConfiguration() {
    // read and double check the simulation time
    double currentTime;
    trajectoryFile >> currentTime >> currentTime >> currentTime;

    if (simulationTime*printingInterval != currentTime) {
        std::cerr << "Inconsistency in the simulation time.\n"
                  << simulationTime*printingInterval << " != " << currentTime << ".\n";
        exit(1);
    }

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
            orientationAutocorrelation += ipcOrientation[j]*ipcInitialOrientations[i][j];
            velocityAutocorrelation += ipcCenterVelocity[j]*ipcCentersInitialVelocities[i][j];
        }
    }
    autocorrelationsFile << currentTime << "\t"
                         << orientationAutocorrelation/nIPCs << "\t"
                         << velocityAutocorrelation/nIPCs << "\n";
}
