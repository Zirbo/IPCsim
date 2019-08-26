#ifndef __IPCSIMULATOR_HEADER_INCLUDED__
#define __IPCSIMULATOR_HEADER_INCLUDED__


#include <fstream>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include "../IPC.hpp"
#include "../helpers/cell_lists.hpp"
#include "../helpers/randomNumberGenerator.hpp"
#include "../helpers/isotropic_pair_correlation_function.hpp"



struct SimulationStage {
    double inputStartingTemperature;
    double inputStageTotalDuration;
    bool inputRestoringPreviousSimulation;
    bool inputPrintTrajectoryAndCorrelations;

    SimulationStage() : inputStartingTemperature{0.}, inputStageTotalDuration{0.},
                        inputRestoringPreviousSimulation{false}, inputPrintTrajectoryAndCorrelations{false}
                    {}
 };

class IPCsimulation {
public:
    IPCsimulation(SimulationStage const& stage);
    void printPotentials();
    double run();


private:
    IPCsimulation();

    unsigned long simulationTime;
    std::ofstream outputFile;
    std::ofstream energyTrajectoryFile;
    std::ofstream trajectoryFile;

    // force and potential tables computation
    std::vector<double> uBB, uBs1, uBs2, us1s2, us1s1, us2s2;
    std::vector<double> fBB, fBs1, fBs2, fs1s2, fs1s1, fs2s2;
    void compileForceAndPotentialTables();
    double computeOmega(double Ra, double Rb, double rab);
    double computeOmegaRadialDerivative(double Ra, double Rb, double rab);

    bool printTrajectoryAndCorrelations;
    // state point
    int nIPCs;
    double density, temperature;
    double initialTemperature;
    // simulation duration
    double simulationTotalDuration;
    double simulationTimeStep, printingInterval;
    // potential
    double e_BB, e_Bs1, e_Bs2, e_s1s2, e_s1s1, e_s2s2, e_min;
    double ipcRadius, firstPatchRadius, firstPatchEccentricity, secndPatchRadius, secndPatchEccentricity;
    double fakeHSexponent, fakeHScoefficient;
    double forceAndEnergySamplingStep, tollerance;
    // masses and inverse masses
    double ipcCenterMass, firstPatchMass, secndPatchMass;
    double ipcCenterInverseMass, firstPatchInverseMass, secndPatchInverseMass;
    // work parameters
    double ratioBetweenTemperatureAndKineticEnergy, totalEnergy, potentialEnergy, kineticEnergy, simulationBoxSide, dt;
    double squaredMinimumDistanceBetweenParticles;
    double interactionRange, squaredInteractionRange;
    double patchDistance, squaredPatchDistance;
    double cP11, cP12, cP1c, cP21, cP22, cP2c, alpha_1, alpha_2, alpha_sum;
    int nPrints;
    // external electric field
    bool isFieldEnabled;
    double ratioChargeFirstPatchOverIpcCenter, ratioChargeSecndPatchOverIpcCenter;
    double externalFieldIpcCenter[3], externalFieldFirstPatch[3], externalFieldSecndPatch[3];
    // particles
    std::vector<IPC> particles;
    cell_lists cells;

    IsotropicPairCorrelationFunction pairCorrelation;

    // selfexplanatory
    void initializeSystem(SimulationStage const& stage);
    void restorePreviousConfiguration();
    void initializeNewConfiguration(int N1);


    void computeTrajectoryStep();

    void computeVerletHalfStepForIPC(IPC & ipc);
    void computeVerletHalfStep();

    void finishVerletStepForIPC(IPC & ipc);
    void finishVerletStep();

    void computeFreeForces();
    struct loopVariables {
        std::vector<std::array<double, 3>> ipcCenterF;
        std::vector<std::array<double, 3>> firstPatchF;
        std::vector<std::array<double, 3>> secndPatchF;
        double U, minimumSquaredDistance;
        loopVariables(size_t nIPCs) : U{0.}, minimumSquaredDistance{1.} {
            ipcCenterF.resize(nIPCs, {0.0, 0.0, 0.0});
            firstPatchF.resize(nIPCs, {0.0, 0.0, 0.0});
            secndPatchF.resize(nIPCs, {0.0, 0.0, 0.0});
        }
    };
    void computeInteractionsWithIPCsInTheSameCell(std::list<int>::const_iterator loc, std::list<int> const& ipcsInCurrentCell, loopVariables & loopVars);
    void computeInteractionsWithIPCsInNeighbouringCells(std::list<int>::const_iterator loc, std::list<int> const& ipcsInNeighbouringCells, loopVariables & loopVars);
    void computeInteractionsBetweenTwoIPCs(const int firstIPC, const int secndIPC, loopVariables & loopVars);

    void outputSystemTrajectory(std::ofstream & outputTrajectoryFile);
    void outputSystemEnergies(std::ofstream &energyTrajectoryFile);


    // 3D boundary conditions enforcers
    void computeSystemEnergy();
    void computeSystemMomentum(double (&pcm) [3]);
    void correctTotalMomentumToZero(double (&pcm)[3], double (&pcmCorrected)[3]);
    void scaleVelocities(const double scalingFactor);
    inline void absolutePBC(double & x) {  x -= std::floor(x);  }
    inline void relativePBC(double & x) {  x -= std::round(x);  }
    // Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
    void generateRandomOrientation(double (&a)[3], RandomNumberGenerator & r);


};

#endif //__IPCSIMULATOR_HEADER_INCLUDED__
