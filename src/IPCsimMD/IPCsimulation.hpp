#ifndef __IPCSIMULATOR_HEADER_INCLUDED__
#define __IPCSIMULATOR_HEADER_INCLUDED__

/*----------------------------------------------------------------------------------
 * Simulation of Inverse Patchy Colloids with two asymmetric patches.
 *
 * It reads the system variables from a file called input.in expected in the same
 * directory of the executable (for the format of the file, peep in function
 * "initializeSystem"). Temperature and duration are to be given to the constructor
 * as parameters (peep structure "SimulationStage").
 *
 * Each IPC is represented as a center of mass plus two patches, all on a line.
 * This linear geometry makes most constraint algorithms singular, so IPCs had
 * to be reproduced using the Ciccotti reduction algorithm [1], according two which
 * only the patches are really moved according to effective forces that
 * take into account the inertia of the center of mass, which then is always
 * found in the middle point between the patches (or wherever it should be, if
 * the patches are not symmetric).
 *
 * The equation of motion for the two patches are integrated using the
 * well-known Velocity-Verlet algorithm, and the patch-patch distance is
 * constrained using the RATTLE algorithm [2]
 *
 * [1] Ciccotti, Ferrario and Ryckaert, Mol. Phys. 47-6, 1253-1264 (1982).
 * [2] Andersen, J. Comp. Phys. 52, 24-34 (1983)
 *----------------------------------------------------------------------------------*/


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
    bool janusSimulation;
    int binaryMixturePercentage;
    bool printForces;

    SimulationStage() : inputStartingTemperature{0.}, inputStageTotalDuration{0.},
                        inputRestoringPreviousSimulation{false}, inputPrintTrajectoryAndCorrelations{false},
                        janusSimulation{false}, printForces{false}
                    {}
};

class IPCsimulation {
public:
    IPCsimulation(SimulationStage const& stage);
    void printPotentials();
    double run();


private:
    IPCsimulation();

    bool isJanusSimulation;
    int binaryMixtureComposition;
    bool isNotJanusSimulation() { return !isJanusSimulation; } // needed because I am blind
    bool printForces;
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

    void printPotentialsToFile(int potentialPrintingStep, int cutoffValue);
    void printPotentialsToFileJanus(int potentialPrintingStep, int cutoffValue);

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
    double patchDistance, squaredPatchDistance, inversePatchDistance;
    double halfDtFirstPatchInverseMass, halfDtSecndPatchInverseMass;
    double inverseAlpha_sumSquaredPatchDistance;
    double cP11, cP12, cP1c, cP21, cP22, cP2c, alpha_1, alpha_2, alpha_sum;
    int nPrints;
    // external electric field
    bool isFieldEnabled;
    double ratioChargeFirstPatchOverIpcCenter, ratioChargeSecndPatchOverIpcCenter;
    double externalFieldIpcCenter[3], externalFieldFirstPatch[3], externalFieldSecndPatch[3];
    // particles
    std::vector<IPC> particles;
    std::vector<JanusIPC> janusParticles;
    cell_lists cells;

    IsotropicPairCorrelationFunction pairCorrelation;

    // selfexplanatory
    void initializeSystem(SimulationStage const& stage);
    void restorePreviousConfiguration();
    void restorePreviousJanusConfiguration();
    void initializeNewConfiguration(int N1);
    void initializeNewJanusConfiguration(int N1);


    void computeTrajectoryStep();

    void computeVerletHalfStepForIPC(IPC & ipc);
    void computeVerletHalfStepForJanusIPC(JanusIPC & ipc);
    void computeVerletHalfStep();

    void finishVerletStepForIPC(IPC & ipc);
    void finishVerletStepForJanusIPC(JanusIPC & ipc);
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
    void computeFreeJanusForces();
    struct loopVariablesJanus {
        std::vector<std::array<double, 3>> ipcCenterF;
        std::vector<std::array<double, 3>> janusPatchF;
        double U, minimumSquaredDistance;
        loopVariablesJanus(size_t nIPCs) : U{0.}, minimumSquaredDistance{1.} {
            ipcCenterF.resize(nIPCs, {0.0, 0.0, 0.0});
            janusPatchF.resize(nIPCs, {0.0, 0.0, 0.0});
        }
    };
    void computeInteractionsWithIPCsInTheSameCell(std::list<int>::const_iterator loc, std::list<int> const& ipcsInCurrentCell, loopVariables & loopVars);
    void computeInteractionsWithJanusIPCsInTheSameCell(std::list<int>::const_iterator loc, std::list<int> const& ipcsInCurrentCell, loopVariablesJanus & loopVars);
    void computeInteractionsWithIPCsInNeighbouringCells(std::list<int>::const_iterator loc, std::list<int> const& ipcsInNeighbouringCells, loopVariables & loopVars);
    void computeInteractionsWithJanusIPCsInNeighbouringCells(std::list<int>::const_iterator loc, std::list<int> const& ipcsInNeighbouringCells, loopVariablesJanus & loopVars);
    void computeInteractionsBetweenTwoIPCs(const int firstIPC, const int secndIPC, loopVariables & loopVars);
    void computeInteractionsBetweenTwoJanusIPCs(const int firstIPC, const int secndIPC, loopVariablesJanus &loopVars);

    void outputSystemTrajectory(std::ofstream & outputTrajectoryFile, const bool printTheForces);
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
