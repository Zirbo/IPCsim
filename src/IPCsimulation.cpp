#include <cstdlib>
#include <iomanip>
#include <iostream>
#include "IPCsimulation.hpp"
#define PI 3.1415926535897932

IPCsimulation::IPCsimulation(bool restorePreviousSimulation) {
  // clean up old data and recreate output directory
  system("rm -rf siml");
  if(system("mkdir siml") != 0) {
    std::cerr<<"Unable to create 'siml/' directory...\n";
    exit(1);
  }

  // open output files
  outputFile.open("siml/output.out");
  trajectoryFile.open("siml/trajectory.xyz");
  energyTrajectoryFile.open("siml/evolution.out");

  // initialize system
  initializeSystem(restorePreviousSimulation);

  // print starting configuration and initialize output file
  outputSystemState();
  outputFile<<"\nPlot evolution.out to check the evolution of the system.\n";
  energyTrajectoryFile<<std::scientific<<std::setprecision(10);
  energyTrajectoryFile<<"#t\t\tT\t\tK\t\tU\t\tE\t\trminbb\t\trminbs\t\trminss\t\trmin\n";
  energyTrajectoryFile<<simulationTime*simulationParameters.dt_nonscaled<<"\t"<<simulationParameters.kT<<"\t"<<simulationParameters.K/simulationParameters.nIPCs<<"\t"<<simulationParameters.U/simulationParameters.nIPCs<<"\t"<<simulationParameters.E/simulationParameters.nIPCs;
  energyTrajectoryFile<<"\t"<<simulationParameters.rminbb*simulationParameters.L<<"\t"<<simulationParameters.rminbs*simulationParameters.L<<"\t"<<simulationParameters.rminss*simulationParameters.L<<"\t"<<sqrt(simulationParameters.rmin2)*simulationParameters.L<<std::endl;
}

void IPCsimulation::run() {
  time_t simulationStartTime, simulationEndTime;

  const unsigned long simulationDurationInIterations = unsigned long(simulationParameters.SimLength/simulationParameters.dt_nonscaled);
  const unsigned long printingInterval = unsigned long(simulationParameters.PrintEvery/simulationParameters.dt_nonscaled);

  // simulation begins
  time(&simulationStartTime);
  while(simulationTime < simulationDurationInIterations) {
    velocityVerletIteration();
    if( simulationTime%printingInterval == 0)
      outputSystemState();
  }

  // check that total momentum is still zero and print final stuff
  space::vec pcm = space::vec( 0., 0., 0. );
  for(int i=0;i<3*simulationParameters.nIPCs;i++)
    pcm += v[i];
  pcm *= simulationParameters.L;

  outputFINALSystemState();
  time(&simulationEndTime);
  outputFile << "The simulation lasted " << difftime (simulationEndTime,simulationStartTime) << " seconds.\n";
  outputFile << "Residual momentum of the whole system = ( " << pcm.x << ", " << pcm.y << ", " << pcm.z << " ).\n" << std::endl;
}





//************************************************************************//
double IPCsimulation::omega(double Ra, double Rb, double rab) {
  // BKL paper, formula 18
  if ( rab > Ra+Rb )
    return 0.;
  else if ( rab <= std::fabs(Ra-Rb) )
    return 8.*std::pow(std::min(Ra,Rb),3);
  else
  {
    const double tempSum = (Ra*Ra-Rb*Rb)/(2.*rab);
    return 2.*( (2.*Ra+tempSum+rab/2.)*pow(Ra-tempSum-rab/2.,2)
            + (2.*Rb-tempSum+rab/2.)*pow(Rb+tempSum-rab/2.,2) );
  }
}

double IPCsimulation::d_dr_omega(double Ra, double Rb, double rab) {
  // BKL paper, derivative of formula 18
  if ( rab >= Ra+Rb || rab <= fabs(Ra-Rb) )    return 0.;
  else
  {
    const double tempSum = (Ra*Ra-Rb*Rb)/(2.*rab);
    const double tempSumMinus = tempSum - rab/2.;
    const double tempSumPlus = tempSum + rab/2.;
    return (6./rab) * (tempSumMinus*(Ra - tempSumPlus)*(Ra + tempSumPlus) - tempSumPlus*(Rb - tempSumMinus)*(Rb + tempSumMinus) );
  }
}
void IPCsimulation::FU_table::make_table(Ensemble par, bool printPotentials)
{
  const size_t potentialRangeSamplingSize = size_t( 2.*par.bigRadius/par.forceAndEnergySamplingStep ) + 1;

  uBB   = new double [potentialRangeSamplingSize];
  uBs1  = new double [potentialRangeSamplingSize];
  uBs2  = new double [potentialRangeSamplingSize];
  us1s2 = new double [potentialRangeSamplingSize];
  us1s1 = new double [potentialRangeSamplingSize];
  us2s2 = new double [potentialRangeSamplingSize];
  fBB   = new double [potentialRangeSamplingSize];
  fBs1  = new double [potentialRangeSamplingSize];
  fBs2  = new double [potentialRangeSamplingSize];
  fs1s2 = new double [potentialRangeSamplingSize];
  fs1s1 = new double [potentialRangeSamplingSize];
  fs2s2 = new double [potentialRangeSamplingSize];

  std::ofstream POT_OUTPUT;
  int potOutputPrintCount = 1;
  if (printPotentials) {
    POT_OUTPUT.open("siml/potentials.out");
    POT_OUTPUT<<std::scientific<<std::setprecision(6);
    POT_OUTPUT<<"#r\t\t\tpotBB\t\t\tpotBs1\t\t\tpotBs2\t\t\tpots1s2\t\t\tpots2s2\t\t\tpots1s1";
    POT_OUTPUT<<  "\t\t\tforBB\t\t\tforBs1\t\t\tforBs2\t\t\tfors1s2\t\t\tfors2s2\t\t\tfors1s1\n";
  }

  for ( size_t i = 0; i < potentialRangeSamplingSize; ++i)
  {
    double r = i*par.forceAndEnergySamplingStep;
    uBB[i]   = (par.e_BB  /par.e_min) * omega(par.bigRadius,par.bigRadius,r);
    uBs1[i]  = (par.e_Bs1 /par.e_min) * omega(par.bigRadius,par.s1Radius,r);
    uBs2[i]  = (par.e_Bs2 /par.e_min) * omega(par.bigRadius,par.s2Radius,r);
    us1s2[i] = (par.e_s1s2/par.e_min) * omega(par.s1Radius,par.s2Radius,r);
    us2s2[i] = (par.e_s2s2/par.e_min) * omega(par.s2Radius,par.s2Radius,r);
    us1s1[i] = (par.e_s1s1/par.e_min) * omega(par.s1Radius,par.s1Radius,r);

    fBB[i]   = (par.e_BB  /par.e_min) * d_dr_omega(par.bigRadius,par.bigRadius,r);
    fBs1[i]  = (par.e_Bs1 /par.e_min) * d_dr_omega(par.bigRadius,par.s1Radius,r);
    fBs2[i]  = (par.e_Bs2 /par.e_min) * d_dr_omega(par.bigRadius,par.s2Radius,r);
    fs1s2[i] = (par.e_s1s2/par.e_min) * d_dr_omega(par.s1Radius,par.s2Radius,r);
    fs2s2[i] = (par.e_s2s2/par.e_min) * d_dr_omega(par.s2Radius,par.s2Radius,r);
    fs1s1[i] = (par.e_s1s1/par.e_min) * d_dr_omega(par.s1Radius,par.s1Radius,r);

    if ( r <= 1.0 )
    {
      // setting up a Fake Hard Sphere Core
      double rm = pow(r,-par.fakeHSexp);
      uBB[i]   += par.fakeHScoef*((rm-2.)*rm+1.);
      fBB[i]   += -2.*par.fakeHSexp*par.fakeHScoef*(rm-1.)*rm/r;
    }
    if ( printPotentials && int( (1000.*i)/potentialRangeSamplingSize ) == potOutputPrintCount )
    {
      potOutputPrintCount++;
      POT_OUTPUT<<r<<"\t"<<uBB[i]<<"\t"<<uBs1[i]<<"\t"<<uBs2[i]<<"\t"<<us1s2[i]<<"\t"<<us2s2[i]<<"\t"<<us1s1[i]<<"\t";
      POT_OUTPUT<<r<<"\t"<<fBB[i]<<"\t"<<fBs1[i]<<"\t"<<fBs2[i]<<"\t"<<fs1s2[i]<<"\t"<<fs2s2[i]<<"\t"<<fs1s1[i]<<"\n";
    }
    // this division is done here to save a division during runtime;
    // it's only done now not to be seen in the plots
    double x = 1./(r);
    fBB[i] *= x;
    fBs1[i] *= x;
    fBs2[i] *= x;
    fs1s2[i] *= x;
    fs1s1[i] *= x;
    fs2s2[i] *= x;
  }
  POT_OUTPUT.close();
}





//************************************************************************//
void IPCsimulation::outputSystemState()
{
  trajectoryFile<<std::scientific<<std::setprecision(24);
  trajectoryFile<<3*simulationParameters.nIPCs<<"\n"<<simulationTime*simulationParameters.dt_nonscaled<<"\n";
  for(int i=0; i<simulationParameters.nIPCs; i++)
  {
    for(int a=0; a<3; a++)
    {
      trajectoryFile<<farben[i+a*simulationParameters.nIPCs]<<"\t";
      trajectoryFile<<x[i+a*simulationParameters.nIPCs].x<<"\t"<<x[i+a*simulationParameters.nIPCs].y<<"\t"<<x[i+a*simulationParameters.nIPCs].z<<"\t";
      trajectoryFile<<v[i+a*simulationParameters.nIPCs].x<<"\t"<<v[i+a*simulationParameters.nIPCs].y<<"\t"<<v[i+a*simulationParameters.nIPCs].z<<"\n";
    }
  }

  energyTrajectoryFile<<simulationTime*simulationParameters.dt_nonscaled<<"\t"<<simulationParameters.kT<<"\t"<<simulationParameters.K/simulationParameters.nIPCs<<"\t"<<simulationParameters.U/simulationParameters.nIPCs<<"\t"<<simulationParameters.E/simulationParameters.nIPCs;
  energyTrajectoryFile<<"\t"<<simulationParameters.rminbb*simulationParameters.L<<"\t"<<simulationParameters.rminbs*simulationParameters.L<<"\t"<<simulationParameters.rminss*simulationParameters.L<<"\t"<<sqrt(simulationParameters.rmin2)*simulationParameters.L<<std::endl;
}





//************************************************************************//
void IPCsimulation::outputFINALSystemState()
{

  std::ofstream outputFile("siml/startingstate.xyz");
  outputFile<<std::scientific<<std::setprecision(24);
  outputFile<<3*simulationParameters.nIPCs<<"\n"<<simulationTime*simulationParameters.dt_nonscaled<<"\n";
  for(int i=0; i<simulationParameters.nIPCs; i++)
  {
    for(int a=0; a<3; a++)
    {
      outputFile<<farben[i+a*simulationParameters.nIPCs]<<"\t";
      outputFile<<x[i+a*simulationParameters.nIPCs].x<<"\t"<<x[i+a*simulationParameters.nIPCs].y<<"\t"<<x[i+a*simulationParameters.nIPCs].z<<"\t";
      outputFile<<v[i+a*simulationParameters.nIPCs].x<<"\t"<<v[i+a*simulationParameters.nIPCs].y<<"\t"<<v[i+a*simulationParameters.nIPCs].z<<"\n";
    }
  }
  outputFile.close();
}





/*****************************************************************************************/
void IPCsimulation::initializeSystem(bool restoreprevious)
{
  simulationTime = 0;

  RandomNumberGenerator rand;
  int N1, N2, N3;
  // input from file
  std::fstream IN("input.in", std::ios::in);
  IN>>N1>>simulationParameters.rho>>simulationParameters.kTimposed;
  IN>>simulationParameters.dt_nonscaled>>simulationParameters.PrintEvery>>simulationParameters.SimLength;
  IN>>simulationParameters.e_BB>>simulationParameters.e_Bs1>>simulationParameters.e_Bs2;
  IN>>simulationParameters.e_s1s1>>simulationParameters.e_s2s2>>simulationParameters.e_s1s2;
  IN>>simulationParameters.e_min;
  IN>>simulationParameters.ecc1>>simulationParameters.s1Radius;
  IN>>simulationParameters.ecc2>>simulationParameters.s2Radius;
  IN>>simulationParameters.m1>>simulationParameters.m2>>simulationParameters.mc;
  IN>>simulationParameters.fakeHScoef>>simulationParameters.fakeHSexp;
  IN>>simulationParameters.forceAndEnergySamplingStep>>simulationParameters.tollerance;
  IN>>simulationParameters.Ec.x>>simulationParameters.Ec.y>>simulationParameters.Ec.z;
  IN>>simulationParameters.qc>>simulationParameters.qp1>>simulationParameters.qp2;
  IN.close();

  // processing the data
  N2 = N1*N1;        N3 = N2*N1;
  simulationParameters.nIPCs = 4*N3;
  simulationParameters.nPatc = 2*simulationParameters.nIPCs;
  if ( abs( (simulationParameters.ecc1+simulationParameters.s1Radius)-(simulationParameters.ecc2+simulationParameters.s2Radius) ) >= 1e-10 )
  {
    std::cerr<<simulationParameters.ecc1<<"+"<<simulationParameters.s1Radius<<"="<<simulationParameters.ecc1+simulationParameters.s1Radius<<"-";
    std::cerr<<simulationParameters.ecc2<<"+"<<simulationParameters.s2Radius<<"="<<simulationParameters.ecc2+simulationParameters.s2Radius<<"=\n";
    std::cerr<<(simulationParameters.ecc1+simulationParameters.s1Radius)-(simulationParameters.ecc2+simulationParameters.s2Radius)<<std::endl;
    std::cerr<<"eccentricities and radii are not consistent!\n"; exit(1);
  }
  simulationParameters.bigRadius = simulationParameters.ecc1 + simulationParameters.s1Radius;
  simulationParameters.L=cbrt(simulationParameters.nIPCs/simulationParameters.rho);
  simulationParameters.L2=simulationParameters.L*simulationParameters.L;
  simulationParameters.kToverK = 2./(5.*simulationParameters.nIPCs-3.);

  if(restoreprevious)
  {
    simulationParameters.kTimposed = sqrt(simulationParameters.kTimposed);
    char culo; double muntonarzu;
    IN.open("startingstate.xyz", std::ios::in);
    IN>>simulationParameters.nIPCs>>muntonarzu;
    simulationParameters.nIPCs /= 3;
    simulationParameters.nPatc = 2*simulationParameters.nIPCs;
    simulationParameters.L = cbrt(simulationParameters.nIPCs/simulationParameters.rho);
    simulationParameters.L2=simulationParameters.L*simulationParameters.L;
    simulationParameters.kToverK = 2./(5.*simulationParameters.nIPCs-3.);
    x = new space::vec[3*simulationParameters.nIPCs]; v = new space::vec[3*simulationParameters.nIPCs]; F = new space::vec[3*simulationParameters.nIPCs];
    farben = new char[3*simulationParameters.nIPCs];

  //  outputFile<<"Reading "<<par.nIPCs<< " particles positions and velocities from file.\n";

    for(int i=0;i<simulationParameters.nIPCs;i++)
    {
      for(int a=0; a<3; a++)
      {
        IN>>culo>>x[i+a*simulationParameters.nIPCs].x>>x[i+a*simulationParameters.nIPCs].y>>x[i+a*simulationParameters.nIPCs].z;
        IN>>v[i+a*simulationParameters.nIPCs].x>>v[i+a*simulationParameters.nIPCs].y>>v[i+a*simulationParameters.nIPCs].z;
        v[i+a*simulationParameters.nIPCs] *= simulationParameters.kTimposed;
        farben[i+a*simulationParameters.nIPCs]=culo;
      }
    }
    IN.close();
  }

  // output the data for future checks
  outputFile<<N1<<"\t"<<simulationParameters.rho<<"\t"<<simulationParameters.kTimposed<<"\n";
  outputFile<<simulationParameters.dt_nonscaled<<"\t"<<simulationParameters.PrintEvery<<"\t"<<simulationParameters.SimLength<<"\n";
  outputFile<<simulationParameters.e_BB<<"\t"<<simulationParameters.e_Bs1<<"\t"<<simulationParameters.e_Bs2<<"\n";
  outputFile<<simulationParameters.e_s1s2<<"\t"<<simulationParameters.e_s1s1<<"\t"<<simulationParameters.e_s2s2<<"\n";
  outputFile<<simulationParameters.e_min<<"\n";
  outputFile<<simulationParameters.ecc1<<"\t"<<simulationParameters.s1Radius<<"\n";
  outputFile<<simulationParameters.ecc2<<"\t"<<simulationParameters.s2Radius<<"\n";
  outputFile<<simulationParameters.m1<<"\t"<<simulationParameters.m2<<"\t"<<simulationParameters.mc;
  outputFile<<simulationParameters.fakeHScoef<<"\t"<<simulationParameters.fakeHSexp<<"\n";
  outputFile<<simulationParameters.forceAndEnergySamplingStep<<"\t"<<simulationParameters.tollerance;
  outputFile<<simulationParameters.Ec.x<<"\t"<<simulationParameters.Ec.y<<"\t"<<simulationParameters.Ec.z<<"\n";
  outputFile<<simulationParameters.qc<<"\t"<<simulationParameters.qp1<<"\t"<<simulationParameters.qp2<<"\n";

  // computing fields
  simulationParameters.Ep1 = simulationParameters.Ec*simulationParameters.qp1;
  simulationParameters.Ep2 = simulationParameters.Ec*simulationParameters.qp2;
  simulationParameters.Ec *= simulationParameters.qc;

  outputFile<<"\n*****************MD simulation in EVN ensemble for CGDH potential.********************\n";
  outputFile<<"\nDensity = "<<simulationParameters.nIPCs<<"/"<<pow(simulationParameters.L,3)<<" = ";
  outputFile<<simulationParameters.nIPCs/pow(simulationParameters.L,3)<<" = "<<simulationParameters.rho<<"\nSide = "<<simulationParameters.L<<std::endl;
  outputFile<<"Total number of simulated particles: "<<simulationParameters.nIPCs+simulationParameters.nPatc<<std::endl;

  // potential sampling
  outputFile<<"Printing potential plots in 'potentials.out'.\n";
  tab.make_table(simulationParameters, true);

  // scaling of lenghts for [0.0:1.0] simulation box
  simulationParameters.bigRadius /= simulationParameters.L;
  simulationParameters.s1Radius /= simulationParameters.L;
  simulationParameters.ecc1 /= simulationParameters.L;
  simulationParameters.s2Radius /= simulationParameters.L;
  simulationParameters.ecc2 /= simulationParameters.L;
  simulationParameters.dt = simulationParameters.dt_nonscaled/simulationParameters.L;
  simulationParameters.forceAndEnergySamplingStep /= simulationParameters.L;
  simulationParameters.PotRange = 2*simulationParameters.bigRadius;
  simulationParameters.PotRangeSquared = simulationParameters.PotRange*simulationParameters.PotRange;
  simulationParameters.PatchDistance = simulationParameters.ecc1+simulationParameters.ecc2;
  simulationParameters.PatchDistanceSquared = simulationParameters.PatchDistance*simulationParameters.PatchDistance;
  simulationParameters.im1 = 1./simulationParameters.m1;
  simulationParameters.im2 = 1./simulationParameters.m2;
  simulationParameters.imc = 1./simulationParameters.mc;
  // inverse of the I parameter from formulas!
  double iI = 1./(simulationParameters.PatchDistanceSquared*simulationParameters.imc + simulationParameters.ecc1*simulationParameters.ecc1*simulationParameters.im2 + simulationParameters.ecc2*simulationParameters.ecc2*simulationParameters.im1);
  simulationParameters.cP11 = 1.-simulationParameters.ecc2*simulationParameters.ecc2*iI*simulationParameters.im1;
  simulationParameters.cP12 = -simulationParameters.ecc1*simulationParameters.ecc2*iI*simulationParameters.im2;
  simulationParameters.cP1c = simulationParameters.PatchDistance*simulationParameters.ecc2*iI*simulationParameters.imc;
  simulationParameters.cP21 = simulationParameters.cP12;
  simulationParameters.cP22 = 1.-simulationParameters.ecc1*simulationParameters.ecc1*iI*simulationParameters.im2;
  simulationParameters.cP2c = simulationParameters.PatchDistance*simulationParameters.ecc1*iI*simulationParameters.imc;
  simulationParameters.alpha_1 = 1. - simulationParameters.ecc2*iI*(simulationParameters.ecc2*simulationParameters.im1-simulationParameters.ecc1*simulationParameters.im2);
  simulationParameters.alpha_2 = 1. + simulationParameters.ecc1*iI*(simulationParameters.ecc2*simulationParameters.im1-simulationParameters.ecc1*simulationParameters.im2);
  simulationParameters.alpha_sum = simulationParameters.alpha_1 + simulationParameters.alpha_2;
  simulationParameters.Ec  /= simulationParameters.L;
  simulationParameters.Ep1 /= simulationParameters.L;
  simulationParameters.Ep2 /= simulationParameters.L;

  // initialize positions of the particles
  if(!restoreprevious)
  {
    x = new space::vec[3*simulationParameters.nIPCs]; v = new space::vec[3*simulationParameters.nIPCs]; F = new space::vec[3*simulationParameters.nIPCs];
    farben = new char[3*simulationParameters.nIPCs];

    // scaling: sqrt(2kT/mPI) comes from boltzmann average of |v_x|
    // 4 is to compensate that the average |x| from x=rand.getRandom55() is 0.25
    double vel_scaling = 4.*sqrt(2.*simulationParameters.kTimposed/PI)/simulationParameters.L;
    // initialize IPC positions
    for(int i=0;i<N3;i++)
    {
      // FCC is obtained as 4 intersecating SC
      x[i]      = space::vec( i%N1 + .1*rand.getRandom55(), (i/N1)%N1 + .1*rand.getRandom55(), i/N2 + .1*rand.getRandom55() )/N1;
        floorccp(x[i]);
      x[i+N3]   = space::vec( .5+i%N1 + .1*rand.getRandom55(), .5+(i/N1)%N1 + .1*rand.getRandom55(), i/N2 + .1*rand.getRandom55() )/N1;
        floorccp(x[i+N3]);
      x[i+2*N3] = space::vec( i%N1 + .1*rand.getRandom55(), .5+(i/N1)%N1 + .1*rand.getRandom55(), .5+i/N2 + .1*rand.getRandom55() )/N1;
        floorccp(x[i+2*N3]);
      x[i+3*N3] = space::vec( .5+i%N1 + .1*rand.getRandom55(), (i/N1)%N1 + .1*rand.getRandom55(), .5+i/N2 + .1*rand.getRandom55() )/N1;
        floorccp(x[i+3*N3]);
      // starting from random but ONLY TRANSLATIONAL speeds, for compatibility with rattle
      v[i]          = space::vec( rand.getRandom55()*vel_scaling, rand.getRandom55()*vel_scaling, rand.getRandom55()*vel_scaling );
      v[i+N3]       = space::vec( rand.getRandom55()*vel_scaling, rand.getRandom55()*vel_scaling, rand.getRandom55()*vel_scaling );
      v[i+N3+N3]    = space::vec( rand.getRandom55()*vel_scaling, rand.getRandom55()*vel_scaling, rand.getRandom55()*vel_scaling );
      v[i+N3+N3+N3] = space::vec( rand.getRandom55()*vel_scaling, rand.getRandom55()*vel_scaling, rand.getRandom55()*vel_scaling );
    }
    //v[1]=v[2]=v[3]=space::vec(0.,0.,0.);
    // initialize patches positions
    for(int i=0;i<simulationParameters.nIPCs;i++)
    {
      space::vec a;        ranor(a,rand);

      x[i+simulationParameters.nIPCs] = x[i] + a*simulationParameters.ecc1;
        floorccp(x[i+simulationParameters.nIPCs]);
      x[i+simulationParameters.nPatc] = x[i] - a*simulationParameters.ecc2;
        floorccp(x[i+simulationParameters.nPatc]);

      space::vec surra;    ranor(surra,rand);
      surra = surra*sqrt(.5*(v[i]*v[i])/(surra*surra));
      surra = surra - a*(surra*a); // now it's orthogonal!

      v[i+simulationParameters.nPatc] = v[i] + surra;
      v[i+simulationParameters.nIPCs] = v[i] - surra;

      farben[i] = 'C';
      farben[i+simulationParameters.nIPCs] = 'P';
      farben[i+simulationParameters.nPatc] = 'G';
    }
  }

  // cell list compilation
  cells.initialize(1.,simulationParameters.PotRange,simulationParameters.nIPCs,x);
  outputFile<<"Total number of cells: "<<cells.M3<<std::endl;
  cells.compilelists(x);

  // first computation of forces
  computeFreeForce();

  // check that total momentum is zero
  space::vec pcm = space::vec( 0., 0., 0. );
  for(int i=0;i<3*simulationParameters.nIPCs;i++)
    pcm += v[i];
  pcm *= simulationParameters.L;
  outputFile<<"P whole system = ( "<<pcm.x<<", "<<pcm.y<<", "<<pcm.z<<" )."<<std::endl;

  // if not restoring, correct the total momentum to be zero
      if(!restoreprevious)
      {
        space::vec pcm_after = space::vec( 0., 0., 0. );
        pcm /= 3*simulationParameters.nIPCs*simulationParameters.L;
        for(int i=0;i<3*simulationParameters.nIPCs;i++)
        {
          v[i] -= pcm;
          pcm_after += v[i];
        }
        pcm_after *= simulationParameters.L;
        outputFile<<"P whole system corrected = ( "<<pcm_after.x<<", "<<pcm_after.y<<", "<<pcm_after.z<<" )."<<std::endl;
      }

  // compute kinetic and total energy, temperature
  simulationParameters.K = 0.;
  double kc(0.),k1(0.),k2(0.);
  for(int i=0; i<simulationParameters.nIPCs; i++)
  {
    kc += v[i]*v[i];
    k1 += v[i+simulationParameters.nIPCs]*v[i+simulationParameters.nIPCs];
    k2 += v[i+simulationParameters.nPatc]*v[i+simulationParameters.nPatc];
  }
  simulationParameters.K = (kc*simulationParameters.mc + k1*simulationParameters.m1 + k2*simulationParameters.m2)*.5*simulationParameters.L2;
  simulationParameters.kT = simulationParameters.K*simulationParameters.kToverK;
  simulationParameters.E = simulationParameters.K + simulationParameters.U;
}





/*****************************************************************************************/
void IPCsimulation::computeFreeForce()
{
  // Computes the force without accounting for constrains.
  // Force on i = sum over j of dU(r_ij)/dr * (x_j-x_i)/r_ij
  for(int i=0;i<3*simulationParameters.nIPCs;i++)
    F[i] = space::vec(0.0,0.0,0.0);
  simulationParameters.U = 0.0;  simulationParameters.rmin2 = simulationParameters.rminbb = simulationParameters.rminbs = simulationParameters.rminss = 1.;

  #pragma omp parallel
  {
    space::vec * Floop = new space::vec [3*simulationParameters.nIPCs];
    double Uloop(0.), rmin2loop(1.), rminbbloop(1.), rminbsloop(1.), rminssloop(1.);
    for(int i=0;i<3*simulationParameters.nIPCs;i++)    Floop[i] = space::vec(0.0,0.0,0.0);

    #pragma omp for
    for(int m=0; m<cells.M3; m++)  // loop over all cells
    {
      std::list<int> neighbs, local;
      cells.neighbour_cells(m,local,neighbs);

      for( std::list<int>::iterator loc = local.begin(); loc!=local.end(); loc++)
      { // loop over particles in neighbouring cells
        for( std::list<int>::iterator ext = neighbs.begin(); ext!=neighbs.end(); ext++)
        {
          space::vec ff, rji;        double r;
          // loop over cm and patches
          for (int i=0;i<3;i++)
          {
            for (int j=0;j<3;j++)
            {
              rji = x[*loc+i*simulationParameters.nIPCs]-x[*ext+j*simulationParameters.nIPCs];        lroundccp(rji);
              r = rji*rji;
              // Store the smallest distance between attraction centers:
              if( r<rmin2loop) rmin2loop=r;
              if (r <= simulationParameters.PotRangeSquared)
              {
                r=sqrt(r);                      int dist = int( r/simulationParameters.forceAndEnergySamplingStep );
                if (i==0 && j==0)        // means cm -> cm
                {
                  ff = rji*(tab.fBB[dist]);  Uloop += tab.uBB[dist];
                  if(r<rminbbloop) rminbbloop=r;    // Store the smallest distance between IPCs:
                }
                else if ( (i==0 && j==1) || (i==1 && j==0) )  // means cm -> patch1
                {
                  ff = rji*(tab.fBs1[dist]);  Uloop += tab.uBs1[dist];
                  if(r<rminbsloop) rminbsloop=r;    // Store the smallest distance between an IPC center and a patch:
                }
                else if ( (i==0 && j==2) || (i==2 && j==0) )  // means cm -> patch2
                {
                  ff = rji*(tab.fBs2[dist]);  Uloop += tab.uBs2[dist];
                  if(r<rminbsloop) rminbsloop=r;    // Store the smallest distance between an IPC center and a patch:
                }
                else if (i==2 && j==2)  // means  patch2 -> patch2
                {
                  ff = rji*(tab.fs2s2[dist]);  Uloop += tab.us2s2[dist];
                  if(r<rminssloop) rminssloop=r;    // Store the smallest distance between patches of different IPCs:
                }
                else if (i==1 && j==1)  // means  patch1 -> patch1
                {
                  ff = rji*(tab.fs1s1[dist]);  Uloop += tab.us1s1[dist];
                  if(r<rminssloop) rminssloop=r;    // Store the smallest distance between patches of different IPCs:
                }
                else if ( (i==1 && j==2) || (i==2 && j==1) )  // means patch1 -> patch2
                {
                  ff = rji*(tab.fs1s2[dist]);  Uloop += tab.us1s2[dist];
                  if(r<rminssloop) rminssloop=r;    // Store the smallest distance between patches of different IPCs:
                }
                else  { std::cerr<<"CRAAAAAAAAAAAAAAAAAAP!"<<std::endl; exit(1);}
                Floop[*loc+i*simulationParameters.nIPCs] -= ff;
                Floop[*ext+j*simulationParameters.nIPCs] += ff;
              }
            }
          }
        }
        // loop over members of the cell m (internal interactions!)
        for( std::list<int>::iterator ins = loc; ins!=local.end(); ins++)
        {
          // starts from loc which is like summing over i >= j inside the cell
          // with a list, you have to access from loc because there's no way
          // to directly access (loc+1); so you need the following
          if(ins == loc) continue;  // to skip i=j iteration
          // after that, same code
          space::vec ff, rji;        double r;
          // loop over cm and patches
          for (int i=0;i<3;i++)
          {
            for (int j=0;j<3;j++)
            {
              rji = x[*loc+i*simulationParameters.nIPCs]-x[*ins+j*simulationParameters.nIPCs];        lroundccp(rji);
              r = rji*rji;
              // Store the smallest distance between attraction centers:
              if( r<rmin2loop) rmin2loop=r;
              if (r <= simulationParameters.PotRangeSquared)
              {
                r=sqrt(r);                      int dist = int( r/simulationParameters.forceAndEnergySamplingStep );
                if (i==0 && j==0)        // means cm -> cm
                {
                  ff = rji*(tab.fBB[dist]);  Uloop += tab.uBB[dist];
                  if(r<rminbbloop) rminbbloop=r;    // Store the smallest distance between IPCs:
                }
                else if ( (i==0 && j==1) || (i==1 && j==0) )  // means cm -> patch1
                {
                  ff = rji*(tab.fBs1[dist]);  Uloop += tab.uBs1[dist];
                  if(r<rminbsloop) rminbsloop=r;    // Store the smallest distance between an IPC center and a patch:
                }
                else if ( (i==0 && j==2) || (i==2 && j==0) )  // means cm -> patch2
                {
                  ff = rji*(tab.fBs2[dist]);  Uloop += tab.uBs2[dist];
                  if(r<rminbsloop) rminbsloop=r;    // Store the smallest distance between an IPC center and a patch:
                }
                else if (i==2 && j==2)  // means  patch2 -> patch2
                {
                  ff = rji*(tab.fs2s2[dist]);  Uloop += tab.us2s2[dist];
                  if(r<rminssloop) rminssloop=r;    // Store the smallest distance between patches of different IPCs:
                }
                else if (i==1 && j==1)  // means  patch1 -> patch1
                {
                  ff = rji*(tab.fs1s1[dist]);  Uloop += tab.us1s1[dist];
                  if(r<rminssloop) rminssloop=r;    // Store the smallest distance between patches of different IPCs:
                }
                else if ( (i==1 && j==2) || (i==2 && j==1) )  // means patch1 -> patch2
                {
                  ff = rji*(tab.fs1s2[dist]);  Uloop += tab.us1s2[dist];
                  if(r<rminssloop) rminssloop=r;    // Store the smallest distance between patches of different IPCs:
                }
                else  { std::cerr<<"CRAAAAAAAAAAAAAAAAAAP!"<<std::endl; exit(1);}
                Floop[*loc+i*simulationParameters.nIPCs] -= ff;
                Floop[*ins+j*simulationParameters.nIPCs] += ff;
              }
            }
          }
        }
      }
    }
    #pragma omp critical
    {
      for(int i=0;i<3*simulationParameters.nIPCs;i++)
        F[i] += Floop[i];
      simulationParameters.U += Uloop;
      if(rmin2loop < simulationParameters.rmin2) simulationParameters.rmin2 = rmin2loop;
      if(rminbbloop < simulationParameters.rminbb) simulationParameters.rminbb = rminbbloop;
      if(rminbsloop < simulationParameters.rminbs) simulationParameters.rminbs = rminbsloop;
      if(rminssloop < simulationParameters.rminss) simulationParameters.rminss = rminssloop;
    }
    delete [] Floop;
  }
  for(int i=0;i<simulationParameters.nIPCs;i++)
    F[i] += simulationParameters.Ec;
  for(int i=simulationParameters.nIPCs;i<simulationParameters.nPatc;i++)
    F[i] += simulationParameters.Ep1;
  for(int i=simulationParameters.nPatc;i<=simulationParameters.nIPCs+simulationParameters.nPatc;i++)
    F[i] += simulationParameters.Ep1;
}





void IPCsimulation::velocityVerletIteration()
{
  /*
   * Every IPC has 3 points: the cm and two patches.
   *
   * Only the patches are really moved according to effective forces that
   * take into account the cm inertia, then the cm is moved in the mean point.
   * This is done with degree of freedom reducting algorithm from
   * Ciccotti, Ferrario and Ryckaert, Mol. Phys. 47-6, 1253-1264 (1982).
   *
   * The patches are then constrained to be at a fixed distance using RATTLE
   * Andersen, J. Comp. Phys. 52, 24-34 (1983)
   * In this procedure rij means r[i]-r[j] as in those papers;
   * some other parts of the code (i.e. free_force) have the opposite convention.
   *
   * The cm velocity is not used, but I still compute it because it's negligibly
   * expensive and I might want to add magnetic interactions sooner or later.
   *
   * The new equations of motion are derived taking into account the asymmetry.
   * They are really easy to derive following Ciccotti's paper.
   *
   */
  for ( int i=0; i<simulationParameters.nIPCs; i++ )
  {
    space::vec eFp1 = F[i+simulationParameters.nIPCs]*simulationParameters.cP11 + F[i+simulationParameters.nPatc]*simulationParameters.cP12 + F[i]*simulationParameters.cP1c;
    space::vec eFp2 = F[i+simulationParameters.nIPCs]*simulationParameters.cP21 + F[i+simulationParameters.nPatc]*simulationParameters.cP22 + F[i]*simulationParameters.cP2c;
    v[i+simulationParameters.nIPCs] += eFp1*(.5*simulationParameters.dt*simulationParameters.im1);
    v[i+simulationParameters.nPatc] += eFp2*(.5*simulationParameters.dt*simulationParameters.im2);
    space::vec r1 = x[i+simulationParameters.nIPCs] + v[i+simulationParameters.nIPCs]*simulationParameters.dt;      floorccp(r1);
    space::vec r2 = x[i+simulationParameters.nPatc] + v[i+simulationParameters.nPatc]*simulationParameters.dt;      floorccp(r2);

    // start working with constrains
    space::vec r12 = r1-r2;          lroundccp(r12);
    double diff = r12*r12-simulationParameters.PatchDistanceSquared;
    while( fabs(diff) > simulationParameters.tollerance*simulationParameters.PatchDistanceSquared )
    {
      space::vec r12old = x[i+simulationParameters.nIPCs]-x[i+simulationParameters.nPatc];      lroundccp(r12old);
      double g = diff/( 2.*(r12*r12old)*simulationParameters.alpha_sum*simulationParameters.dt );
      space::vec DX = r12old*g;
      v[i+simulationParameters.nIPCs] -= DX*simulationParameters.alpha_1;
      v[i+simulationParameters.nPatc] += DX*simulationParameters.alpha_2;
        DX *= simulationParameters.dt;
      r1 -= DX*simulationParameters.alpha_1;      floorccp(r1);
      r2 += DX*simulationParameters.alpha_2;      floorccp(r2);
      r12 = r1-r2;   lroundccp(r12);

      diff = r12*r12-simulationParameters.PatchDistanceSquared;
    }
    x[i+simulationParameters.nIPCs] = r1;
    x[i+simulationParameters.nPatc] = r2;
    x[i] = r2 + r12*simulationParameters.ecc2/simulationParameters.PatchDistance;  // with their notation r12=r1<-2
    floorccp(x[i]);
  }

  cells.compilelists(x);                // rewrite cell lists for the new iteration
  computeFreeForce();      // compute F(x[t+dt]) and the potential

  simulationParameters.K = 0.;
  double virialconstrains(0.);

  for ( int i=0; i<simulationParameters.nIPCs; i++ )
  {
    space::vec eFp1 = F[i+simulationParameters.nIPCs]*simulationParameters.cP11 + F[i+simulationParameters.nPatc]*simulationParameters.cP12 + F[i]*simulationParameters.cP1c;
    space::vec eFp2 = F[i+simulationParameters.nIPCs]*simulationParameters.cP21 + F[i+simulationParameters.nPatc]*simulationParameters.cP22 + F[i]*simulationParameters.cP2c;
    space::vec q1   = v[i+simulationParameters.nIPCs] + eFp1*(.5*simulationParameters.dt*simulationParameters.im1);
    space::vec q2   = v[i+simulationParameters.nPatc] + eFp2*(.5*simulationParameters.dt*simulationParameters.im2);

    // start working with constrains
    space::vec r12 = x[i+simulationParameters.nIPCs]-x[i+simulationParameters.nPatc];      lroundccp(r12);
    space::vec q12 = q1-q2;
    double k = (r12*q12)/( simulationParameters.alpha_sum*simulationParameters.PatchDistanceSquared );
    while( fabs(k) > simulationParameters.tollerance )
    {
      space::vec DX = r12*k;
      q1 -= DX*simulationParameters.alpha_1;
      q2 += DX*simulationParameters.alpha_2;
      q12 = q1-q2;
      k = (r12*q12)/( simulationParameters.alpha_sum*simulationParameters.PatchDistanceSquared );
    }
    virialconstrains += k;
    v[i+simulationParameters.nIPCs] = q1;
    v[i+simulationParameters.nPatc] = q2;
    v[i] = (q1*simulationParameters.ecc2+q2*simulationParameters.ecc1)/simulationParameters.PatchDistance;
    simulationParameters.K += (q1*q1)*simulationParameters.m1 + (q2*q2)*simulationParameters.m2 + (v[i]*v[i])*simulationParameters.mc;
  }
  virialconstrains *= simulationParameters.PatchDistance;
  simulationParameters.K *= .5*simulationParameters.L2;
  simulationParameters.E = simulationParameters.K + simulationParameters.U;
  simulationParameters.kT = simulationParameters.kToverK*simulationParameters.K;
  simulationTime++;
}

// Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
void IPCsimulation::ranor(space::vec & a, RandomNumberGenerator & r)
{
  double x,y,quad=2.;
  while ( quad > 1. )  {    x = r.getRandom11();    y = r.getRandom11();    quad = x*x + y*y;  }
  double norm = 2.*sqrt(1.-quad);  a.x=x*norm;  a.y=y*norm;  a.z=1.-2.*quad;
}
