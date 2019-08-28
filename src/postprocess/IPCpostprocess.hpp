#ifndef __IPCPOSTPROCESS_HEADER_INCLUDED__
#define __IPCPOSTPROCESS_HEADER_INCLUDED__

/*---------------------------------------------------------------------------------------
 * Postprocess for Inverse Patchy Colloids simulations with variable number of patches.
 *
 * asda
 *---------------------------------------------------------------------------------------*/


#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cmath>


class IPCpostprocess {
public:
    IPCpostprocess(int inputNumberOfPatches, std::string const& directoryName);
    void run();


private:
    IPCpostprocess();

    int numberOfPatches;
    unsigned long simulationTime;
    std::ofstream autocorrelationsFile;
    std::ofstream typicalOrientationsFile;
    std::ofstream meanSquaredDisplFile;
    std::ifstream trajectoryFile;
    // state point
    int nIPCs;
    double density, temperature;
    // simulation duration
    double simulationTotalDuration;
    double simulationTimeStep, printingInterval;
    // geometry
    double ipcRadius, firstPatchRadius, firstPatchEccentricity;
    double secndPatchRadius, secndPatchEccentricity;
    // work parameters
    double simulationBoxSide, dt;
    double squaredMinimumDistanceBetweenParticles;
    double interactionRange, squaredInteractionRange;
    double patchDistance, squaredPatchDistance;
    int nPrints;

    std::vector<std::array<double, 3>> ipcCentersInitialPositions;
    std::vector<std::array<double, 3>> ipcCentersPreviousPositions;
    std::vector<std::array<double, 3>> ipcInitialOrientations;
    std::vector<std::array<double, 3>> ipcPreviousOrientations;
    std::vector<std::array<double, 3>> ipcCentersInitialVelocities;

    std::vector<double> meanSquaredDisplacement;
    std::vector<double> totalRotation;
    std::vector<double> orientationAutocorrelation;
    std::vector<double> velocityAutocorrelation;

    void initialize(std::string const& directoryName);
    void readFirstConfiguration();
    void processSingleConfiguration();
    inline void relativePBC(double & x) {  x -= std::round(x);  }
    void runConsistencyChecks();
};

#endif //__IPCPOSTPROCESS_HEADER_INCLUDED__
