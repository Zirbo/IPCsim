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
    IPCpostprocess(int inputNumberOfPatches,
                   size_t const inputNumberOfSubSimulations,
                   std::string const& directoryName);
    void run();


private:
    IPCpostprocess();

    int numberOfPatches;
    // output files
    std::ofstream autocorrelationsFile;
    std::ofstream meanSquaredDisplFile;
//    std::ofstream typicalOrientationsFile;
    // input trajectory
    std::ifstream trajectoryFile;
    // state point
    int nIPCs;
    double density, temperature;
    // simulation duration
    double simulationTotalDuration;
    size_t simulationDurationInIterations;
    double simulationTimeStep, printingInterval;
    // subsimulation working parameters
    size_t numberOfSubSimulations, subSimulationDuration;
    // geometry
    double ipcRadius, firstPatchRadius, firstPatchEccentricity;
    double secndPatchRadius, secndPatchEccentricity;
    // work parameters
    double simulationBoxSide, dt;
    double squaredMinimumDistanceBetweenParticles;
    double interactionRange, squaredInteractionRange;
    double patchDistance, squaredPatchDistance;
    int nPrints;

    typedef std::vector<std::array<double, 3>> spaceVector;
    spaceVector ipcCentersPreviousPositions;
    spaceVector ipcCentersCurrentPositions;
    spaceVector ipcInitialOrientations;
    spaceVector ipcCurrentOrientations;
    spaceVector ipcCentersInitialVelocities;
    spaceVector ipcCentersCurrentVelocities;

    std::vector<double> orientationAutocorrelation;
    std::vector<double> velocityAutocorrelation;

    void initialize(std::string const& directoryName);
    void updateInitialOrientationAndVelocites();
    void updatePreviousPositions();
    void computeMSD(size_t const snapshotNumber);
    void computeAutocorrelations(size_t const snapshotNumber);
    inline void relativePBC(double & x) {  x -= std::round(x);  }
    void readOutputFile(std::string const& directoryName);
    void readSnapshot(size_t const snapshotNumber);
    void runConsistencyChecks(size_t const snapshotNumber);
    void printAutocorrelations(size_t const numberOfSubSimulations);
};

#endif //__IPCPOSTPROCESS_HEADER_INCLUDED__
