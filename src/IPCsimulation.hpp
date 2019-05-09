#ifndef __IPCSIMULATOR_HEADER_INCLUDED__
#define __IPCSIMULATOR_HEADER_INCLUDED__


#include <fstream>
#include <string>
#include "zilvectors.hpp"
#include "cell_lists.hpp"
#include "randomNumberGenerator.hpp"


class IPCsimulation {
public:
    IPCsimulation(bool restorePreviousSimulation);
    void run();

private:
    IPCsimulation();
    struct Ensemble {
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
        double m1, m2, mc, im1, im2, imc;
        // simulation
        double dt_nonscaled;
        double SimLength;
        double PrintEvery;
        // work parameters
        double kToverK, E, U, K, L, L2, dt;
        double rmin2, rminbb, rminbs, rminss;
        double PotRange, PotRangeSquared;
        double PatchDistance, PatchDistanceSquared;
        double cP11, cP12, cP1c, cP21, cP22, cP2c, alpha_1, alpha_2, alpha_sum;
        int nPrints;
        int nPatc;
        // external field on cm and patches
        space::vec Ec, Ep1, Ep2;
        double qc, qp1, qp2;
    } simulationParameters;

    unsigned long simulationTime;
    space::vec *x, *v, *F;
    char *farben;
    cell_lists cells;
    std::ofstream outputFile;
    std::ofstream energyTrajectoryFile;
    std::ofstream trajectoryFile;

    // force and potential tables computation
    struct FU_table
    {
        double *uBB, *uBs1, *uBs2, *us1s2, *us1s1, *us2s2;
        double *fBB, *fBs1, *fBs2, *fs1s2, *fs1s1, *fs2s2;
        void make_table(Ensemble simulationParameters, bool printPotentials);
    } tab;

    static double omega(double Ra, double Rb, double rab);
    static double d_dr_omega(double Ra, double Rb, double rab);
    // simulation parts
    void computeFreeForce();
    void velocityVerletIteration();
    // selfexplanatory
    void initializeSystem(bool restoreprevious);
    void outputSystemState();
    void outputFINALSystemState();

private:
    // small helpers
    // 3D boundary conditions enforcers
    inline void floorccp(space::vec & a)  {  a.x-=std::floor(a.x);   a.y-=std::floor(a.y);   a.z-=std::floor(a.z);   }
    inline void lroundccp(space::vec & a) {  a.x-=std::lround(a.x);  a.y-=std::lround(a.y);  a.z-=std::lround(a.z);  }
    // Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
    void ranor(space::vec & a, RandomNumberGenerator & r);
};

#endif //__IPCSIMULATOR_HEADER_INCLUDED__
