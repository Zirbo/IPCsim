#include <fstream>
#include <string>
#include <vector>
#include <array>


class IPCpostprocess {
public:
    IPCpostprocess(std::string const& directoryName);
    void run();


private:
    IPCpostprocess();

    unsigned long simulationTime;
    std::ofstream autocorrelationsFile;
    std::ofstream typicalOrientationsFile;
    std::ifstream trajectoryFile;
    // state point
    int nIPCs;
    double density, temperature;
    // simulation duration
    double simulationTotalDuration;
    double simulationTimeStep, printingInterval;
    // geometry
    double ipcRadius, firstPatchRadius, firstPatchEccentricity, secndPatchRadius, secndPatchEccentricity;
    // work parameters
    double simulationBoxSide, dt;
    double squaredMinimumDistanceBetweenParticles;
    double interactionRange, squaredInteractionRange;
    double patchDistance, squaredPatchDistance;
    int nPrints;

    std::vector<std::array<double, 3>> ipcCentersPositions;
    std::vector<std::array<double, 3>> ipcCentersVelocities;
    std::vector<std::array<double, 3>> ipcOrientations;

    void initialize(std::string const& directoryName);
    void readFirstConfiguration();
    void readSingleConfiguration();
};
