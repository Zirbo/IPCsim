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
    IPCpostprocess(const int inputNumberOfPatches,
                   const int inputNumberOfSubSimulations,
                   std::string const& directoryName);
    void run();


private:
    IPCpostprocess();

    int numberOfPatches;
    // input trajectory
    std::ifstream trajectoryFile;
    // state point
    int nIPCs;
    double density, temperature;
    // simulation duration
    double simulationTotalDuration;
    int simulationDurationInIterations;
    double simulationTimeStep, printingInterval;
    // subsimulation working parameters
    int numberOfSubSimulations, subSimulationDuration;
    // geometry
    double ipcRadius, firstPatchRadius, firstPatchEccentricity;
    double secndPatchRadius, secndPatchEccentricity;
    // work parameters
    double simulationBoxSide, dt;
    double squaredMinimumDistanceBetweenParticles;
    double interactionRange, squaredInteractionRange;
    double patchDistance, squaredPatchDistance;
    int nPrints;
    int orientationHistogramSize;
    long totalCollectedOrientations;

    typedef std::vector<std::array<double, 3>> spaceVector;
    spaceVector ipcCentersPreviousPositions;
    spaceVector ipcCentersCurrentPositions;
    spaceVector ipcInitialOrientations;
    spaceVector ipcCurrentOrientations;
    spaceVector ipcCentersInitialVelocities;
    spaceVector ipcCentersCurrentVelocities;

    spaceVector displacementOfEachIPCs;
    std::vector<double> orientationAutocorrelation;
    std::vector<double> velocityAutocorrelation;
    std::vector<double> meanSquaredDisplacement;
    std::vector<std::vector<double>> orientationHistogram;

    void initialize(std::string const& directoryName);
    void updateInitialOrientationAndVelocites();
    void updatePreviousPositions();
    inline void relativePBC(double & x) {  x -= std::round(x);  }
    void readOutputFile(std::string const& directoryName);
    void readSnapshot(const int snapshotNumber);
    void runConsistencyChecks(const int snapshotNumber);

    void computeMSD(const int snapshotNumber);
    void printMSD();

    void computeAutocorrelations(const int snapshotNumber);
    void printAutocorrelations();

    void accumulateTypicalOrientations();
    void printTypicalOrientations();

};

#endif //__IPCPOSTPROCESS_HEADER_INCLUDED__
