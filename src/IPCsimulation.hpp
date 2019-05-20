#ifndef __IPCSIMULATOR_HEADER_INCLUDED__
#define __IPCSIMULATOR_HEADER_INCLUDED__


#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "IPC.hpp"
#include "cell_lists.hpp"
#include "randomNumberGenerator.hpp"




class IPCsimulation {
public:
    IPCsimulation(bool restorePreviousSimulation);
    void run();

private:
    IPCsimulation();

    unsigned long simulationTime;
    std::ofstream outputFile;
    std::ofstream energyTrajectoryFile;
    std::ofstream trajectoryFile;

    // force and potential tables computation
    //struct FU_table
    //{
        double *uBB, *uBs1, *uBs2, *us1s2, *us1s1, *us2s2;
        double *fBB, *fBs1, *fBs2, *fs1s2, *fs1s1, *fs2s2;
        void make_table(bool printPotentials);
    //private:
        double omega(double Ra, double Rb, double rab);
        double d_dr_omega(double Ra, double Rb, double rab);
   // } tab;

    // state point
    int nIPCs;
    double rho;
    double kT, kTimposed;
    // potential
    double e_BB, e_Bs1, e_Bs2, e_s1s2, e_s1s1, e_s2s2, e_min;
    double bigRadius, s1Radius, ecc1, s2Radius, ecc2;
    double fakeHSexp, fakeHScoef;
    double forceAndEnergySamplingStep, tollerance;
    // masses and inverse masses
    double mass[3], inverseMass[3];
    // simulation
    double dt_nonscaled;
    double SimLength;
    double PrintEvery;
    // work parameters
    double kToverK, E, U, K, L, L2, dt;
    double rmin2;
    double PotRange, PotRangeSquared;
    double PatchDistance, PatchDistanceSquared;
    double cP11, cP12, cP1c, cP21, cP22, cP2c, alpha_1, alpha_2, alpha_sum;
    int nPrints;
    int nPatc;
    // particles
    std::vector<IPC> particles;
    cell_lists cells;

    // selfexplanatory
    void initializeSystem(bool restoreprevious);
    void restorePreviousConfiguration();
    void initializeNewConfiguration(int N1);


    void computeTrajectoryStep();

    void computeVerletHalfStepForIPC(IPC & ipc);   
    void computeVerletHalfStep();

    void finishVerletStepForIPC(IPC & ipc);
    void finishVerletStep();

    void computeFreeForces();
    struct loopVariables {
        std::vector<std::vector<double>> force;
        double U, minimumSquaredDistance;
        loopVariables() : U{0.}, minimumSquaredDistance{1.} {}
    };
    void computeInteractionsWithIPCsInTheSameCell(std::list<int>::iterator loc, std::list<int> ipcsInCurrentCell, loopVariables & loopVars);
    void computeInteractionsWithIPCsInNeighbouringCells(std::list<int>::iterator loc, std::list<int> ipcsInNeighbouringCells, loopVariables & loopVars);
    void computeInteractionsBetweenTwoIPCs(int firstIPC, int secndIPC, loopVariables & loopVars);

    void outputSystemState(std::ofstream & outputTrajectoryFile, std::ofstream &energyTrajectoryFile, unsigned long simulationTime);


    // 3D boundary conditions enforcers
    void computeSystemMomentum(double (&pcm) [3]);
    inline void floorccp(double & x, double & y, double &z)  {  x-=std::floor(x);   y-=std::floor(y);   z-=std::floor(z);   }
    inline void floorccp(double & x)  {  x-=std::floor(x);  }
    inline void lroundccp(double & x, double & y, double &z) {  x-=std::lround(x);  y-=std::lround(y);  z-=std::lround(z);  }
    inline void lroundccp(double & x) {  x-=std::lround(x);  }
    // Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
    void ranor(double (&a)[3], RandomNumberGenerator & r);


};

#endif //__IPCSIMULATOR_HEADER_INCLUDED__
