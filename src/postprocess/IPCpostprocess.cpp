#include "IPCpostprocess.hpp"
#include <cstdlib>
#include <iostream>
#include <cmath>


IPCpostprocess::IPCpostprocess(std::string const& directoryName) {
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
    // double check, remove later!
    std::cout << nIPCs << "\t" << density << "\t" << temperature << "\n"
            << simulationTimeStep << "\t" << printingInterval << "\t" << simulationTotalDuration << "\n"
            << firstPatchEccentricity << "\t" << firstPatchRadius << "\n"
            << secndPatchEccentricity << "\t" << secndPatchRadius << "\n";

    // initialize the containers
    ipcCentersPositions.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcCentersVelocities.resize(nIPCs, {0.0, 0.0, 0.0});
    ipcOrientations.resize(nIPCs, {0.0, 0.0, 0.0});


}

void IPCpostprocess::run() {

    readSingleConfiguration();
}

void IPCpostprocess::readSingleConfiguration() {
    trajectoryFile.ignore(200, '\n');
    trajectoryFile.ignore(200, '\n');
    for(int i = 0; i < nIPCs; ++i) {

    }
}
