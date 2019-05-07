#ifndef __IPCSIMULATOR_HEADER_INCLUDED__
#define __IPCSIMULATOR_HEADER_INCLUDED__


#include <fstream>
#include <string>
#include "zilrandom.hpp"
#include "zilvectors.hpp"
#include "cell_lists.hpp"


class IPCsimulation {
public:
    //IPCsimulator();
    void run(bool doWarmup);

private:
    struct Ensemble {
        // state point
        int nIPCs;
        double rho;
        double kT, kTimposed;
        // potential
        double e_BB, e_Bs1, e_Bs2, e_s1s2, e_s1s1, e_s2s2, e_min;
        double bigRadius, s1Radius, ecc1, s2Radius, ecc2;
        double FakeHSexp, FakeHScoef;
        double FUsamplingDeltar, Tollerance;
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
    } par;

    unsigned long simulationTime;
    space::vec *x, *v, *F;
    char *farben;
    cell_lists cells;
    std::ofstream outputFile;
    std::ofstream energyTrajectoryFile;


    // force and potential tables computation
    struct FU_table
    {
        double *uBB, *uBs1, *uBs2, *us1s2, *us1s1, *us2s2;
        double *fBB, *fBs1, *fBs2, *fs1s2, *fs1s1, *fs2s2;
        void make_table(Ensemble par);
    } tab;

    static double omega(double Ra, double Rb, double rab);
    static double d_dr_omega(double Ra, double Rb, double rab);
    // simulation parts
    void computeFreeForce();
    void velocityVerletIteration();
    // selfexplanatory
    void warmup(bool restoreprevious);
    void outputSystemState(std::string nome, bool append);
};

#endif //__IPCSIMULATOR_HEADER_INCLUDED__
