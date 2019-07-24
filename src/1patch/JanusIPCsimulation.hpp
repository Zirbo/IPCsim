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




class JanusIPCsimulation {
public:
    JanusIPCsimulation(bool restorePreviousSimulation, bool stagingEnabled = false, std::pair<double,int> stage = std::pair<double,int>{0.,0});
    double run();


private:
    JanusIPCsimulation();

    unsigned long simulationTime;
    std::ofstream outputFile;
    std::ofstream energyTrajectoryFile;
    std::ofstream trajectoryFile;

    // force and potential tables computation
    std::vector<double> uBB, uBs, uss;
    std::vector<double> fBB, fBs, fss;
    void make_table(bool printPotentials);
    double omega(double Ra, double Rb, double rab);
    double d_dr_omega(double Ra, double Rb, double rab);

    // state point
    int nIPCs;
    double density, temperature;
    double desiredTemperature;
    // simulation duration
    double simulationTotalDuration;
    double simulationTimeStep, printingInterval;
    // potential
    double e_BB, e_Bs, e_ss, e_min;
    double ipcRadius, patchRadius, patchCenterDistance;
    double fakeHSexponent, fakeHScoefficient;
    double forceAndEnergySamplingStep, tollerance;
    // masses and inverse masses
    double centerMass, patchMass;
    double centerInverseMass, patchInverseMass;
    // work parameters
    double ratioBetweenTemperatureAndKineticEnergy, totalEnergy, potentialEnergy, kineticEnergy, simulationBoxSide, dt;
    double squaredMinimumDistanceBetweenParticles;
    double squaredPatchCenterDistance;
    double interactionRange, squaredInteractionRange;
    int nPrints;
    // external electric field
    bool isFieldEnabled;
    double ratioChargePatchOverIpcCenter;
    double externalFieldIpcCenter[3], externalFieldPatch[3];
    // particles
    std::vector<JanusIPC> particles;
    cell_lists cells;

    // selfexplanatory
    void initializeSystem(bool restoreprevious, bool stagingEnabled, const std::pair<double,int> & stage);
    void restorePreviousConfiguration();
    void initializeNewConfiguration(int N1);


    void computeTrajectoryStep();

    void computeVerletHalfStepForIPC(JanusIPC & ipc);
    void computeVerletHalfStep();

    void finishVerletStepForIPC(JanusIPC & ipc);
    void finishVerletStep();

    void computeFreeForces();
    struct loopVariables {
        std::vector<std::array<double, 3>> ipcCenterF;
        std::vector<std::array<double, 3>> janusPatchF;
        double U, minimumSquaredDistance;
        loopVariables(size_t nIPCs) : U{0.}, minimumSquaredDistance{1.} {
            ipcCenterF.resize(nIPCs, {0.0, 0.0, 0.0});
            janusPatchF.resize(nIPCs, {0.0, 0.0, 0.0});
        }
    };
    void computeInteractionsWithIPCsInTheSameCell(std::list<int>::const_iterator loc, std::list<int> const& ipcsInCurrentCell, loopVariables & loopVars);
    void computeInteractionsWithIPCsInNeighbouringCells(std::list<int>::const_iterator loc, std::list<int> const& ipcsInNeighbouringCells, loopVariables & loopVars);
    void computeInteractionsBetweenTwoIPCs(const int firstIPC, const int secndIPC, loopVariables & loopVars);

    void outputSystemTrajectory(std::ofstream & outputTrajectoryFile);
    void outputSystemState(std::ofstream & outputTrajectoryFile, std::ofstream &energyTrajectoryFile);


    // 3D boundary conditions enforcers
    void computeSystemEnergy();
    void computeSystemMomentum(double (&pcm) [3]);
    void correctTotalMomentumToZero(double (&pcm)[3], double (&pcmCorrected)[3]);
    void scaleVelocities(const double scalingFactor);
    inline void absolutePBC(double & x) {  x -= std::floor(x);  }
    inline void relativePBC(double & x) {  x -= std::lround(x);  }
    // Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
    void ranor(double (&a)[3], RandomNumberGenerator & r);


};

#endif //__IPCSIMULATOR_HEADER_INCLUDED__
