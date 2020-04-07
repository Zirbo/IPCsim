#include "IPCpostprocess.hpp"
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <algorithm>


IPCpostprocess::IPCpostprocess(std::string const& trajFilename, const std::string &inputFilename) {

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

    // read input file
    std::ifstream inputFile(inputFilename);
    if(inputFile.fail()) {
        std::cerr << "File " << inputFilename << " could not be opened. Aborting.\n";
        exit(1);
    }
    inputFile >> patchRadius >> patchEccentricity;
    ipcRadius = 0.5;
    interactionRange = 2*(patchEccentricity + patchRadius);


    // open simulation output file and trajectory
    trajectoryFile.open(trajFilename);
    if(trajectoryFile.fail()) {
        std::cerr << "File " << trajFilename << " could not be opened. Aborting.\n";
        exit(1);
    }
}

void IPCpostprocess::run() {
    readFirstConfiguration();
    initializeDataAnalysis();

    while (trajectoryFile.peek() != EOF) {
        readNewConfiguration();
        doDataAnalysis();
    }
}

void IPCpostprocess::readFirstConfiguration() {
    int timestep;
    trajectoryFile.ignore(200, '\n');
    trajectoryFile >> timestep;
    trajectoryFile.ignore(200, '\n');
    trajectoryFile.ignore(200, '\n');
    trajectoryFile >> nIPCs;
    nIPCs /=3;
    trajectoryFile.ignore(200, '\n');
    trajectoryFile.ignore(200, '\n');
    trajectoryFile >> boxSideX >> boxSideX;
    trajectoryFile >> boxSideY >> boxSideY;
    trajectoryFile >> boxSideZ >> boxSideZ;
    trajectoryFile.ignore(200, '\n');
    trajectoryFile.ignore(200, '\n');
    particles.resize(nIPCs);

    readIPCconfiguration();

    std::cout << timestep << "\n";
}

bool IPCpostprocess::readNewConfiguration() {
    int timestep;
    trajectoryFile.ignore(200, '\n');
    trajectoryFile.ignore(200, '\n');
    trajectoryFile >> timestep;
    if (trajectoryFile.eof())
        return false;
    // ignore the rest of the header, it must not change
    for(int i = 0; i < 8; ++i)
        trajectoryFile.ignore(200, '\n');

    readIPCconfiguration();

    std::cout << timestep << "\n";
    return true;
}

void IPCpostprocess::readIPCconfiguration() {
    int dummy;
    int ipcNumber;
    for ( IPC &ipc: particles) {
        trajectoryFile >> ipcNumber >> dummy
                >> ipc.ipcCenter.x[0] >> ipc.ipcCenter.x[1] >> ipc.ipcCenter.x[2]
                >> dummy >> dummy
                >> ipc.firstPatch.x[0] >> ipc.firstPatch.x[1] >> ipc.firstPatch.x[2]
                >> dummy >> dummy
                >> ipc.secndPatch.x[0] >> ipc.secndPatch.x[1] >> ipc.secndPatch.x[2];

        ipc.number = ipcNumber/3;
        ipc.type = 'C';
        /*for (int i: {0, 1, 2}) {
            ipc.ipcCenter.x[i] /= boxSideX;
            ipc.firstPatch.x[i] /= boxSideY;
            ipc.secndPatch.x[i] /= boxSideZ;
        }*/
    }
}

void IPCpostprocess::initializeDataAnalysis() {
    //
}

void IPCpostprocess::doDataAnalysis() {
    //
}
