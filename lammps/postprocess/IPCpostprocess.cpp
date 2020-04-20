#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "IPCpostprocess.hpp"
#include "IPCpostprocessNeighbourAnalysis.hpp"
#include "IPCpostprocessOrientationsAnalysis.hpp"
#include "IPCpairCorrelation.hpp"

IPCpostprocess::IPCpostprocess(std::string const& trajFilename, std::string const& inputFilename, std::string const& potDirName) {

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

    potential.initialize(potDirName);
}

void IPCpostprocess::run() {
    readFirstConfiguration();

    IPCneighboursAnalysis neighbourAnalysis(boxSide, interactionRange);
    IPCorientationsAnalysis orientationsAnalysis;
    IPCisotropicPairCorrelationFunction g_r(50, boxSide, nIPCs);

    neighbourAnalysis.accumulate(potential, ipcs);
    orientationsAnalysis.accumulate(ipcOrientations);
    g_r.accumulate(ipcs);

    while (trajectoryFile.peek() != EOF) {
        readNewConfiguration();
        neighbourAnalysis.accumulate(potential, ipcs);
        orientationsAnalysis.accumulate(ipcOrientations);
        g_r.accumulate(ipcs);
    }
    neighbourAnalysis.print("orientationAnalysis.out");
    orientationsAnalysis.print("orientationAnalysis.out");
    g_r.print("g_r.out");
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
    trajectoryFile >> boxSide[0] >> boxSide[0];
    trajectoryFile >> boxSide[1] >> boxSide[1];
    trajectoryFile >> boxSide[2] >> boxSide[2];
    trajectoryFile.ignore(200, '\n');
    trajectoryFile.ignore(200, '\n');
    ipcs.resize(nIPCs);
    ipcOrientations.resize(nIPCs);

    readIPCconfiguration();
    computeOrientations();

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
    computeOrientations();

    std::cout << timestep << "\n";
    return true;
}

void IPCpostprocess::readIPCconfiguration() {
    int dummy;
    int ipcNumber;
    for ( IPC &ipc: ipcs) {
        trajectoryFile >> ipcNumber >> dummy
                >> ipc.ipcCenter.x[0] >> ipc.ipcCenter.x[1] >> ipc.ipcCenter.x[2]
                >> dummy >> dummy
                >> ipc.firstPatch.x[0] >> ipc.firstPatch.x[1] >> ipc.firstPatch.x[2]
                >> dummy >> dummy
                >> ipc.secndPatch.x[0] >> ipc.secndPatch.x[1] >> ipc.secndPatch.x[2];

        ipc.number = ipcNumber/3;
        ipc.type = 'C';
    }
}

void IPCpostprocess::computePerfectOrientations() {
    const double norm[3] = {boxSide[0]/patchEccentricity, boxSide[1]/patchEccentricity, boxSide[2]/patchEccentricity};

    double checkSum;
    for (IPC ipc: ipcs) {
        if (ipc.number < 3)
            checkSum = 0.;
        for (int d: DIMENSIONS) {
            ipcOrientations[ipc.number][d] = ipc.firstPatch.x[d] - ipc.ipcCenter.x[d];
            relativePBC(ipcOrientations[ipc.number][d]);
            ipcOrientations[ipc.number][d] *= norm[d];
            if (ipc.number < 3)
                checkSum += std::pow(ipcOrientations[ipc.number][d],2);
        }
        if (ipc.number < 3) {
            checkSum = std::sqrt(checkSum);
            double diff = std::fabs(checkSum - 1.0);
            if ( diff > 5e-2) {
                std::cerr << __func__ << ":: wrong normalization of the orientation of ipc " << ipc.number << ": " << checkSum*patchEccentricity<<" !\n";
                exit(1);
            }
        }
    }
}

void IPCpostprocess::computeOrientations() {
    for (IPC ipc: ipcs) {
        double norm = 0.;
        for (int d: DIMENSIONS) {
            ipcOrientations[ipc.number][d] = ipc.firstPatch.x[d] - ipc.ipcCenter.x[d];
            relativePBC(ipcOrientations[ipc.number][d]);
            norm += std::pow(ipcOrientations[ipc.number][d],2);
        }
        norm = std::pow(norm, -0.5);
        for (int d: DIMENSIONS) {
            ipcOrientations[ipc.number][d] *= norm;
        }

    }
}
