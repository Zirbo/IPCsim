#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iomanip>
#include "IPCsimulation.hpp"

//************************************************************************//
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
  outputFile<<"\nPlot evolution.out to check the evolution of the system.\n";

  energyTrajectoryFile<<std::scientific<<std::setprecision(10);
  energyTrajectoryFile<<"#t\t\tT\t\tK\t\tU\t\tE\t\trmin\n";
  energyTrajectoryFile<<simulationTime*dt_nonscaled<<"\t"<<kT<<"\t"<<K/nIPCs<<"\t"<<U/nIPCs<<"\t"<<E/nIPCs;
  energyTrajectoryFile<<"\t"<<sqrt(rmin2)*L<<std::endl;

  outputSystemState(trajectoryFile, energyTrajectoryFile, 0);
}

//************************************************************************//
void IPCsimulation::run() {
    time_t simulationStartTime, simulationEndTime;

    unsigned long simulationDurationInIterations = (unsigned long)SimLength/dt_nonscaled;
    unsigned long printingInterval = (unsigned long)PrintEvery/dt_nonscaled;

    // simulation begins
    time(&simulationStartTime);
    while(simulationTime <= simulationDurationInIterations) {
        computeTrajectoryStep();
        ++simulationTime;

        if( simulationTime%printingInterval == 0)
        outputSystemState(trajectoryFile, energyTrajectoryFile, simulationTime);
    }

    // check that total momentum is still zero and print final stuff
    double pcm [3];
    computeSystemMomentum(pcm);

    std::ofstream finalStateFile("startingstate.xyz");
    outputSystemState(finalStateFile, energyTrajectoryFile, simulationTime);
    finalStateFile.close();

    time(&simulationEndTime);
    outputFile << "The simulation lasted " << difftime (simulationEndTime,simulationStartTime) << " seconds.\n";
    outputFile << "Residual momentum of the whole system = ( " << pcm[0] << ", " << pcm[1] << ", " << pcm[2] << " ).\n" << std::endl;
}



//************************************************************************//
void IPCsimulation::computeSystemMomentum(double (&pcm)[3]) {
    for (unsigned short i: {0, 1, 2})
        pcm[i] = 0.;
    for(IPC ipc: particles) {
        for (unsigned short i: {0, 1, 2}) {
            pcm[i] += mass[0]*ipc.ipcCenter.v[i] + mass[1]*ipc.firstPatch.v[i] + mass[2]*ipc.secndPatch.v[i];
        }
    }
    for (unsigned short i: {0, 1, 2})
        pcm[i] *= L;
}





/*****************************************************************************************/
void IPCsimulation::initializeSystem(bool restoreprevious)
{
  simulationTime = 0;

  int N1;
  // input from file
  std::fstream IN("input.in", std::ios::in);
  IN>>N1>>rho>>kTimposed;
  IN>>dt_nonscaled>>PrintEvery>>SimLength;
  IN>>e_BB>>e_Bs1>>e_Bs2;
  IN>>e_s1s1>>e_s2s2>>e_s1s2;
  IN>>e_min;
  IN>>ecc1>>s1Radius;
  IN>>ecc2>>s2Radius;
  IN>>mass[1]>>mass[2]>>mass[0];
  IN>>fakeHScoef>>fakeHSexp;
  IN>>forceAndEnergySamplingStep>>tollerance;
  //IN>>Ec.x>>Ec.y>>Ec.z;
  //IN>>qc>>qp1>>qp2;
  IN.close();

  // processing the data
  nIPCs = 4*N1*N1*N1;
  nPatc = 2*nIPCs;
  if ( abs( (ecc1+s1Radius)-(ecc2+s2Radius) ) >= 1e-10 )
  {
    std::cerr<<ecc1<<"+"<<s1Radius<<"="<<ecc1+s1Radius<<"-";
    std::cerr<<ecc2<<"+"<<s2Radius<<"="<<ecc2+s2Radius<<"=\n";
    std::cerr<<(ecc1+s1Radius)-(ecc2+s2Radius)<<std::endl;
    std::cerr<<"eccentricities and radii are not consistent!\n"; exit(1);
  }
  bigRadius = ecc1 + s1Radius;
  L=cbrt(nIPCs/rho);
  L2=L*L;
  kToverK = 2./(5.*nIPCs-3.);

    if(restoreprevious)
    {
        outputFile<<"Reading "<<nIPCs<< " particles positions and velocities from file.\n";
        restorePreviousConfiguration();
    }

  // output the data for future checks
  outputFile<<N1<<"\t"<<rho<<"\t"<<kTimposed<<"\n";
  outputFile<<dt_nonscaled<<"\t"<<PrintEvery<<"\t"<<SimLength<<"\n";
  outputFile<<e_BB<<"\t"<<e_Bs1<<"\t"<<e_Bs2<<"\n";
  outputFile<<e_s1s2<<"\t"<<e_s1s1<<"\t"<<e_s2s2<<"\n";
  outputFile<<e_min<<"\n";
  outputFile<<ecc1<<"\t"<<s1Radius<<"\n";
  outputFile<<ecc2<<"\t"<<s2Radius<<"\n";
  outputFile<<mass[1]<<"\t"<<mass[2]<<"\t"<<mass[0];
  outputFile<<fakeHScoef<<"\t"<<fakeHSexp<<"\n";
  outputFile<<forceAndEnergySamplingStep<<"\t"<<tollerance;
  //outputFile<<Ec.x<<"\t"<<Ec.y<<"\t"<<Ec.z<<"\n";
  //outputFile<<qc<<"\t"<<qp1<<"\t"<<qp2<<"\n";

  // computing fields
  /*Ep1 = Ec*qp1;
  Ep2 = Ec*qp2;
  Ec *= qc;*/

  outputFile<<"\n*****************MD simulation in EVN ensemble for CGDH potential.********************\n";
  outputFile<<"\nDensity = "<<nIPCs<<"/"<<pow(L,3)<<" = ";
  outputFile<<nIPCs/pow(L,3)<<" = "<<rho<<"\nSide = "<<L<<std::endl;
  outputFile<<"Total number of simulated particles: "<<nIPCs+nPatc<<std::endl;

  // potential sampling
  outputFile<<"Printing potential plots in 'potentials.out'.\n";
  make_table(true);

  // scaling of lenghts for [0.0:1.0] simulation box
  bigRadius /= L;
  s1Radius /= L;
  ecc1 /= L;
  s2Radius /= L;
  ecc2 /= L;
  dt = dt_nonscaled/L;
  forceAndEnergySamplingStep /= L;
  PotRange = 2*bigRadius;
  PotRangeSquared = PotRange*PotRange;
  PatchDistance = ecc1+ecc2;
  PatchDistanceSquared = PatchDistance*PatchDistance;
  inverseMass[1] = 1./mass[1];
  inverseMass[2] = 1./mass[2];
  inverseMass[0] = 1./mass[0];
  // inverse of the I parameter from formulas!
  double iI = 1./(PatchDistanceSquared*inverseMass[0] + ecc1*ecc1*inverseMass[2] + ecc2*ecc2*inverseMass[1]);
  cP11 = 1.-ecc2*ecc2*iI*inverseMass[1];
  cP12 = -ecc1*ecc2*iI*inverseMass[2];
  cP1c = PatchDistance*ecc2*iI*inverseMass[0];
  cP21 = cP12;
  cP22 = 1.-ecc1*ecc1*iI*inverseMass[2];
  cP2c = PatchDistance*ecc1*iI*inverseMass[0];
  alpha_1 = 1. - ecc2*iI*(ecc2*inverseMass[1]-ecc1*inverseMass[2]);
  alpha_2 = 1. + ecc1*iI*(ecc2*inverseMass[1]-ecc1*inverseMass[2]);
  alpha_sum = alpha_1 + alpha_2;
  /*Ec  /= L;
  Ep1 /= L;
  Ep2 /= L;*/

    // initialize positions of the particles
    if(!restoreprevious) {
        outputFile<<"Placing "<<nIPCs<< " IPCs on a FCC lattice.\n";
        initializeNewConfiguration(N1);
    }

    // cell list compilation
    cells.initialize(1., PotRange, nIPCs);
    outputFile<<"Total number of cells: "<<cells.M3<<std::endl;
    cells.compilelists(particles);

    // first computation of forces
    computeFreeForces();

    // check that total momentum is zero
    double pcm [3];
    computeSystemMomentum(pcm);
    outputFile << "P whole system = ( "
               << pcm[0] << ", "
               << pcm[1] << ", "
               << pcm[2] << " )." << std::endl;

    // if not restoring, correct the total momentum to be zero
    if(!restoreprevious) {
        double pcmCorrected [3] = { 0., 0., 0. };
        for (unsigned short i: {0, 1, 2})
            pcm[i] /= 3*nIPCs*L;
        for(IPC ipc: particles) {
            for (unsigned short i: {0, 1, 2}) {
                ipc.ipcCenter.v[i] -= pcm[i];
                ipc.firstPatch.v[i] -= pcm[i];
                ipc.secndPatch.v[i] -= pcm[i];

                pcmCorrected[i] += ipc.ipcCenter.v[i] + ipc.firstPatch.v[i] + ipc.secndPatch.v[i];
            }
        }
        for (unsigned short i: {0, 1, 2})
             pcmCorrected[i] *= L;
        outputFile << "P whole system corrected = ( "
                   << pcmCorrected[0] << ", "
                   << pcmCorrected[1] << ", "
                   << pcmCorrected[2] << " )." << std::endl;
    }
}





/*****************************************************************************************/


// Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
void IPCsimulation::ranor(double (&a)[3], RandomNumberGenerator & r) {
  double x,y,quad=2.;
  while ( quad > 1. )  {    x = r.getRandom11();    y = r.getRandom11();    quad = x*x + y*y;  }
  double norm = 2.*sqrt(1.-quad);  a[0]=x*norm;  a[1]=y*norm;  a[2]=1.-2.*quad;
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
//************************************************************************//
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

void IPCsimulation::make_table(bool printPotentials)
{
  const size_t potentialRangeSamplingSize = size_t( 2.*bigRadius/forceAndEnergySamplingStep ) + 1;

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
    double r = i*forceAndEnergySamplingStep;
    uBB[i]   = (e_BB  /e_min) * omega(bigRadius, bigRadius, r);
    uBs1[i]  = (e_Bs1 /e_min) * omega(bigRadius, s1Radius,  r);
    uBs2[i]  = (e_Bs2 /e_min) * omega(bigRadius, s2Radius,  r);
    us1s2[i] = (e_s1s2/e_min) * omega(s1Radius,  s2Radius,  r);
    us2s2[i] = (e_s2s2/e_min) * omega(s2Radius,  s2Radius,  r);
    us1s1[i] = (e_s1s1/e_min) * omega(s1Radius,  s1Radius,  r);

    fBB[i]   = (e_BB  /e_min) * d_dr_omega(bigRadius, bigRadius, r);
    fBs1[i]  = (e_Bs1 /e_min) * d_dr_omega(bigRadius, s1Radius,  r);
    fBs2[i]  = (e_Bs2 /e_min) * d_dr_omega(bigRadius, s2Radius,  r);
    fs1s2[i] = (e_s1s2/e_min) * d_dr_omega(s1Radius,  s2Radius,  r);
    fs2s2[i] = (e_s2s2/e_min) * d_dr_omega(s2Radius,  s2Radius,  r);
    fs1s1[i] = (e_s1s1/e_min) * d_dr_omega(s1Radius,  s1Radius,  r);

    if ( r <= 1.0 )
    {
      // setting up a Fake Hard Sphere Core
      double rm = pow(r, -fakeHSexp);
      uBB[i]   += fakeHScoef*((rm-2.)*rm+1.);
      fBB[i]   += -2.*fakeHSexp*fakeHScoef*(rm-1.)*rm/r;
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


void IPCsimulation::computeTrajectoryStep() {
    computeVerletHalfStep();
    cells.compilelists(particles);
    computeFreeForces();
    finishVerletStep();
}

void IPCsimulation::computeVerletHalfStep() {
    for(IPC &ipc: particles) {
        computeVerletHalfStepForIPC(ipc);
    }
}

void IPCsimulation::computeVerletHalfStepForIPC(IPC & ipc) {
    double x1[3], x2[3];
    double dxNew[3];
    for (unsigned short i: {0, 1, 2}) {
        // compute the half step velocities from the effective forces of the last step
        ipc.firstPatch.v[i] += ipc.eFp1[i]*(.5*dt*inverseMass[1]);
        ipc.secndPatch.v[i] += ipc.eFp2[i]*(.5*dt*inverseMass[2]);

        // compute new positions from the half step velocities
        x1[i] = ipc.firstPatch.x[i] + ipc.firstPatch.v[i]*dt;
        floorccp(x1[i]);
        x2[i] = ipc.secndPatch.x[i] + ipc.secndPatch.v[i]*dt;
        floorccp(x2[i]);

        // check how much the constraints are being violated
        dxNew[i] = x1[i] - x2[i];
        lroundccp(dxNew[i]);
    }
    double diff = (dxNew[0]*dxNew[0] + dxNew[1]*dxNew[1] + dxNew[2]*dxNew[2]) - PatchDistanceSquared;

    // correct the positions and velocities until the violation is inside the tollerance
    while( std::abs(diff) > tollerance*PatchDistanceSquared )
    {
        double dxOld[3], DX[3];
        for (unsigned short i: {0, 1, 2}) {
            dxOld[i] = ipc.firstPatch.x[i] - ipc.secndPatch.x[i];
            lroundccp(dxOld[i]);
        }
        double g = diff/( 2.*(dxOld[0]*dxNew[0] + dxOld[1]*dxNew[1] + dxOld[2]*dxNew[2]) * alpha_sum*dt );

        for (unsigned short i: {0, 1, 2}) {
            DX[i] = g*dxOld[i];

            ipc.firstPatch.v[i] -= alpha_1*DX[i];
            ipc.secndPatch.v[i] += alpha_2*DX[i];
            DX[i] *= dt;

            x1[i] -= DX[i]*alpha_1;
            floorccp(x1[i]);
            x2[i] += DX[i]*alpha_2;
            floorccp(x2[i]);

            dxNew[i] = x1[i] - x2[i];
            lroundccp(dxNew[i]);
        }

        diff = (dxNew[0]*dxNew[0] + dxNew[1]*dxNew[1] + dxNew[2]*dxNew[2]) - PatchDistanceSquared;
    }

    for (unsigned short i: {0, 1, 2}) {
        ipc.firstPatch.x[i] = x1[i];
        ipc.secndPatch.x[i] = x2[i];
        ipc.ipcCenter.x[i] = x2[i] + dxNew[i]*ecc2/PatchDistance;
        floorccp(ipc.ipcCenter.x[i]);
    }
}

void IPCsimulation::finishVerletStep() {
    K = 0.;
    for(IPC &ipc: particles) {
        finishVerletStepForIPC(ipc);
    }
    K *= .5*L2;
    E = K + U;
    kT = kToverK*K;
}

void IPCsimulation::finishVerletStepForIPC(IPC & ipc) {
    double v1[3], v2[3], dx[3], dv[3];
    for (unsigned short i: {0, 1, 2}) {
        // compute the new effective forces from the new forces
        ipc.eFp1[i] = ipc.firstPatch.F[i]*cP11 + ipc.secndPatch.F[i]*cP12 + ipc.ipcCenter.F[i]*cP1c;
        ipc.eFp2[i] = ipc.firstPatch.F[i]*cP21 + ipc.secndPatch.F[i]*cP22 + ipc.ipcCenter.F[i]*cP2c;

        // compute the the final velocities from the new effective forces
        v1[i] = ipc.firstPatch.v[i] + ipc.eFp1[i]*(.5*dt*inverseMass[1]);
        v2[i] = ipc.secndPatch.v[i] + ipc.eFp2[i]*(.5*dt*inverseMass[2]);

        // compute the patch-patch distance
        dx[i] = ipc.firstPatch.x[i] - ipc.secndPatch.x[i];
        lroundccp(dx[i]);
        dv[i] = v1[i] - v2[i];
    }
    // check how much the constraints are being violated
    double k = (dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2])/(alpha_sum*PatchDistanceSquared);
    while( fabs(k) > tollerance ) {
        // compute and apply corrections
        double DX[3];
        for (unsigned short i: {0, 1, 2}) {
            DX[i] = k*dx[i];
            v1[i] -= DX[i]*alpha_1;
            v2[i] += DX[i]*alpha_2;
            dv[i] = v1[i] - v2[i];
        }
        // recompute the violation of the constraints
        k = (dv[0]*dx[0] + dv[1]*dx[1] + dv[2]*dx[2])/(alpha_sum*PatchDistanceSquared);
    }

    for (unsigned short i: {0, 1, 2}) {
        ipc.firstPatch.v[i] = v1[i];
        ipc.secndPatch.v[i] = v2[i];
        ipc.ipcCenter.v[i] = (v1[i]*ecc2 + v2[i]*ecc1)/PatchDistance;
    }
    K += mass[1]*(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2])
       + mass[2]*(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2])
       + mass[0]*(ipc.ipcCenter.v[0]*ipc.ipcCenter.v[0] + ipc.ipcCenter.v[1]*ipc.ipcCenter.v[1] + ipc.ipcCenter.v[2]*ipc.ipcCenter.v[2]);
}



//************************************************************************//
void IPCsimulation::outputSystemState(std::ofstream & outputTrajectoryFile, std::ofstream & energyTrajectoryFile, unsigned long simulationTime)
{
    outputTrajectoryFile<<std::scientific<<std::setprecision(24);
    outputTrajectoryFile<<3*nIPCs<<"\n"<<simulationTime*dt_nonscaled<<"\n";
    for (IPC ipc: particles) {
        outputTrajectoryFile << ipc.type << "\t" << ipc.ipcCenter.x[0] << "\t" << ipc.ipcCenter.x[1] << "\t" << ipc.ipcCenter.x[2] << "\n"
                             << "P"      << "\t" << ipc.firstPatch.v[0] << "\t" << ipc.firstPatch.v[1] << "\t" << ipc.firstPatch.v[2] << "\n"
                             << "Q"      << "\t" << ipc.secndPatch.F[0] << "\t" << ipc.secndPatch.F[1] << "\t" << ipc.secndPatch.F[2] << "\n";
    }

    energyTrajectoryFile << simulationTime*dt_nonscaled << "\t" << kT << "\t"
                         << K/nIPCs << "\t" << U/nIPCs << "\t" << E/nIPCs << "\t"
                         << sqrt(rmin2)*L << "\n";
}




void IPCsimulation::initializeNewConfiguration(int N1) {
    particles.resize(nIPCs);
    RandomNumberGenerator rand;

    int N2 = N1*N1;
    int N3 = N2*N1;

    // scaling: sqrt(2kT/mPI) comes from boltzmann average of |v_x|
    // 4 is to compensate that the average |x| from x=rand.getRandom55() is 0.25
    double vel_scaling = 4.*std::sqrt(2.*kTimposed/3.1415)/L;
    // initialize IPC positions
    for(int i=0;i<N3;i++)
    {
      // FCC is obtained as 4 intersecating SC
        particles[i].ipcCenter.x[0] = (i%N1 + .1*rand.getRandom55())/N1;
        floorccp(particles[i].ipcCenter.x[0]);
        particles[i].ipcCenter.x[1] = ((i/N1)%N1 + .1*rand.getRandom55())/N1;
        floorccp(particles[i].ipcCenter.x[1]);
        particles[i].ipcCenter.x[2] = (i/N2 + .1*rand.getRandom55()) /N1;
        floorccp(particles[i].ipcCenter.x[2]);

        particles[i+N3].ipcCenter.x[0] = (.5 + i%N1 + .1*rand.getRandom55())/N1;
        floorccp(particles[i+N3].ipcCenter.x[0]);
        particles[i+N3].ipcCenter.x[1] = (.5 + (i/N1)%N1 + .1*rand.getRandom55())/N1;
        floorccp(particles[i+N3].ipcCenter.x[1]);
        particles[i+N3].ipcCenter.x[2] = (i/N2 + .1*rand.getRandom55()) /N1;
        floorccp(particles[i+N3].ipcCenter.x[2]);

        particles[i+N3+N3].ipcCenter.x[0] = (i%N1 + .1*rand.getRandom55())/N1;
        floorccp(particles[i+N3+N3].ipcCenter.x[0]);
        particles[i+N3+N3].ipcCenter.x[1] = (.5 + (i/N1)%N1 + .1*rand.getRandom55())/N1;
        floorccp(particles[i+N3+N3].ipcCenter.x[1]);
        particles[i+N3+N3].ipcCenter.x[2] = (.5 + i/N2 + .1*rand.getRandom55()) /N1;
        floorccp(particles[i+N3+N3].ipcCenter.x[2]);

        particles[i+N3+N3+N3].ipcCenter.x[0] = (.5 + i%N1 + .1*rand.getRandom55())/N1;
        floorccp(particles[i+N3+N3+N3].ipcCenter.x[0]);
        particles[i+N3+N3+N3].ipcCenter.x[1] = ((i/N1)%N1 + .1*rand.getRandom55())/N1;
        floorccp(particles[i+N3+N3+N3].ipcCenter.x[1]);
        particles[i+N3+N3+N3].ipcCenter.x[2] = (.5 + i/N2 + .1*rand.getRandom55()) /N1;
        floorccp(particles[i+N3+N3+N3].ipcCenter.x[2]);

        // starting from random but ONLY TRANSLATIONAL speeds, for compatibility with rattle
        for (unsigned short j: {0, 1, 2}) {
           particles[i].ipcCenter.v[j]          = rand.getRandom55()*vel_scaling;
           particles[i+N3].ipcCenter.v[j]       = rand.getRandom55()*vel_scaling;
           particles[i+N3+N3].ipcCenter.v[j]    = rand.getRandom55()*vel_scaling;
           particles[i+N3+N3+N3].ipcCenter.v[j] = rand.getRandom55()*vel_scaling;
        }
    }
    // initialize patches positions
    for(IPC &ipc: particles) {
        double ipcAxis[3], ipcOrthogonalAxis[3];
        ranor(ipcAxis,rand);
        ranor(ipcOrthogonalAxis, rand);
        double normOfIpcOrthogonalAxis = std::sqrt( .5*
                (std::pow(ipc.ipcCenter.v[0],2) + std::pow(ipc.ipcCenter.v[1],2) + std::pow(ipc.ipcCenter.v[2],2))
                /
                (std::pow(ipcOrthogonalAxis[0],2) + std::pow(ipcOrthogonalAxis[1],2) + std::pow(ipcOrthogonalAxis[2],2))
                );
        double scalarProductOfTheTwoAxes = ipcOrthogonalAxis[0]*ipcAxis[0] +
                ipcOrthogonalAxis[1]*ipcAxis[1] + ipcOrthogonalAxis[2]*ipcAxis[2];

        for (unsigned short i: {0, 1, 2}) {
            ipcOrthogonalAxis[i] *= normOfIpcOrthogonalAxis;
            ipcOrthogonalAxis[i] -= ipcAxis[i]*scalarProductOfTheTwoAxes;

            ipc.firstPatch.x[i] = ipc.ipcCenter.x[i] + ipcAxis[i]*ecc1;
            floorccp(ipc.firstPatch.x[i]);
            ipc.secndPatch.x[i] = ipc.ipcCenter.x[i] + ipcAxis[i]*ecc2;
            floorccp(ipc.secndPatch.x[i]);

            ipc.firstPatch.v[i] = ipc.ipcCenter.v[i] + ipcOrthogonalAxis[i]*vel_scaling;
        }
    }
}

void IPCsimulation::restorePreviousConfiguration() {
    kTimposed = sqrt(kTimposed);
    char charJunk; double doubleJunk;
    std::ifstream IN("startingstate.xyz");
    IN>>nIPCs>>doubleJunk;
    nIPCs /= 3;
    nPatc = 2*nIPCs;
    L = cbrt(nIPCs/rho);
    L2=L*L;
    kToverK = 2./(5.*nIPCs-3.);

    particles.resize(nIPCs);

    for (IPC &ipc: particles) {
        IN >> ipc.type
           >> ipc.ipcCenter.x[0] >> ipc.ipcCenter.x[1] >> ipc.ipcCenter.x[2]
           >> ipc.ipcCenter.v[0] >> ipc.ipcCenter.v[1] >> ipc.ipcCenter.v[2]
           >> ipc.ipcCenter.F[0] >> ipc.ipcCenter.F[1] >> ipc.ipcCenter.F[2];
        IN >> charJunk
           >> ipc.firstPatch.x[0] >> ipc.firstPatch.x[1] >> ipc.firstPatch.x[2]
           >> ipc.firstPatch.v[0] >> ipc.firstPatch.v[1] >> ipc.firstPatch.v[2]
           >> ipc.firstPatch.F[0] >> ipc.firstPatch.F[1] >> ipc.firstPatch.F[2];
        IN >> charJunk
           >> ipc.secndPatch.x[0] >> ipc.secndPatch.x[1] >> ipc.secndPatch.x[2]
           >> ipc.secndPatch.v[0] >> ipc.secndPatch.v[1] >> ipc.secndPatch.v[2]
           >> ipc.secndPatch.F[0] >> ipc.secndPatch.F[1] >> ipc.secndPatch.F[2];

        // scale velociteis -> move to another fct
        ipc.ipcCenter.v[0] *= kTimposed;
        ipc.ipcCenter.v[1] *= kTimposed;
        ipc.ipcCenter.v[2] *= kTimposed;
        ipc.firstPatch.v[0] *= kTimposed;
        ipc.firstPatch.v[1] *= kTimposed;
        ipc.firstPatch.v[2] *= kTimposed;
        ipc.secndPatch.v[0] *= kTimposed;
        ipc.secndPatch.v[1] *= kTimposed;
        ipc.secndPatch.v[2] *= kTimposed;
    }
    IN.close();
}



void IPCsimulation::computeFreeForces() {
    // Computes the force without accounting for constrains.
    // Force on i = sum over j of dU(r_ij)/dr * (x_j-x_i)/r_ij

    // reset all forces
    for(IPC &ipc: particles) {
        for (unsigned short i: {0, 1, 2}) {
            ipc.ipcCenter.F[i] = 0.;
            ipc.firstPatch.F[i] = 0.;
            ipc.secndPatch.F[i] = 0.;
        }
    }
    U = 0.0;  rmin2 = 1.;

    #pragma omp parallel
    {
        loopVariables loopVars;
        loopVars.force.resize(3*nIPCs, {0.0, 0.0, 0.0});

        #pragma omp for
        for(int m=0; m<cells.M3; m++)  // loop over all cells
        {
            std::list<int> ipcInNeighbouringCells, ipcInCurrentCell;
            cells.neighbour_cells(m,ipcInCurrentCell,ipcInNeighbouringCells);

            for( std::list<int>::iterator ipc = ipcInCurrentCell.begin(); ipc != ipcInCurrentCell.end(); ipc++) {
                computeInteractionsWithIPCsInNeighbouringCells(ipc, ipcInNeighbouringCells, loopVars);
                computeInteractionsWithIPCsInTheSameCell(ipc, ipcInCurrentCell, loopVars);
            }
        }
        #pragma omp critical
        {
            for (size_t j = 0; j < nIPCs; ++j) {
                for (unsigned short i: {0, 1, 2}) {
                    particles[j].ipcCenter.F[i] += loopVars.force[j][i];
                    particles[j].firstPatch.F[i] += loopVars.force[j+nIPCs][i];
                    particles[j].secndPatch.F[i] += loopVars.force[j+nIPCs+nIPCs][i];
                }
            }
            U += loopVars.U;
            if(loopVars.minimumSquaredDistance < rmin2) rmin2 = loopVars.minimumSquaredDistance;
        }
    }

/*
    for (IPC &ipc: particles) {
        for (unsigned short i: {0, 1, 2}) {
            ipc.center.F[i] += Ec[i];
            ipc.firstPatch.F[i] += Ep1[i];
            ipc.secndPatch.F[i] += Ep2[i];
          }
    }
*/
}

void IPCsimulation::computeInteractionsWithIPCsInNeighbouringCells(std::list<int>::iterator loc, std::list<int> ipcsInNeighbouringCells, loopVariables & loopVars) {
    for( auto ext = ipcsInNeighbouringCells.begin(); ext != ipcsInNeighbouringCells.end(); ++ext) {
        computeInteractionsBetweenTwoIPCs(*loc, *ext, loopVars);
    }
}



void IPCsimulation::computeInteractionsWithIPCsInTheSameCell(std::list<int>::iterator loc, std::list<int> ipcsInCurrentCell, loopVariables &loopVars) {
    // ---------------------------------------------------------------
    // SILVANOOOOO something is not working here, the loop is endless
    // ---------------------------------------------------------------
    for( auto ins = loc; ins != ipcsInCurrentCell.end(); ++ins)
    {
      // starts from loc which is like summing over i >= j inside the cell
      // with a list, you have to access from loc because there's no way
      // to directly access (loc+1); so you need the following "continue" statement
      if(ins == loc) continue;  // to skip i=j iteration

      computeInteractionsBetweenTwoIPCs(*loc, *ins, loopVars);
    }
}

void IPCsimulation::computeInteractionsBetweenTwoIPCs(int firstIPC, int secndIPC, loopVariables &loopVars) {
    IPC & first = particles[firstIPC];
    IPC & second = particles[secndIPC];
    double siteSiteSeparation[9][3];
    // center-center
    for (unsigned short i: {0, 1, 2}) {
        siteSiteSeparation[0][i] = first.ipcCenter.x[i] - second.ipcCenter.x[i];
        siteSiteSeparation[1][i] = first.ipcCenter.x[i] - second.firstPatch.x[i];
        siteSiteSeparation[2][i] = first.ipcCenter.x[i] - second.secndPatch.x[i];
        siteSiteSeparation[3][i] = first.firstPatch.x[i] - second.ipcCenter.x[i];
        siteSiteSeparation[4][i] = first.firstPatch.x[i] - second.firstPatch.x[i];
        siteSiteSeparation[5][i] = first.firstPatch.x[i] - second.secndPatch.x[i];
        siteSiteSeparation[6][i] = first.secndPatch.x[i] - second.ipcCenter.x[i];
        siteSiteSeparation[7][i] = first.secndPatch.x[i] - second.firstPatch.x[i];
        siteSiteSeparation[8][i] = first.secndPatch.x[i] - second.secndPatch.x[i];
    }
    for (unsigned short j = 0; j < 9; ++j) {
        double siteSiteSeparationModulus = siteSiteSeparation[j][0]*siteSiteSeparation[j][0]
                                         + siteSiteSeparation[j][1]*siteSiteSeparation[j][1]
                                         + siteSiteSeparation[j][2]*siteSiteSeparation[j][2];

        if (siteSiteSeparationModulus < loopVars.minimumSquaredDistance)
            loopVars.minimumSquaredDistance = siteSiteSeparationModulus;

        if (siteSiteSeparationModulus <= PotRangeSquared) {
            siteSiteSeparationModulus = std::sqrt(siteSiteSeparationModulus);
            const size_t dist = size_t( siteSiteSeparationModulus/forceAndEnergySamplingStep );
            if (j == 0) { // center - center
                loopVars.U += uBB[dist];
                for (unsigned short i: {0, 1, 2}) {
                    const double modulus = fBB[dist]*siteSiteSeparation[j][i];
                    loopVars.force[firstIPC][i] -= modulus;
                    loopVars.force[secndIPC][i] += modulus;
                }
            } else if (j == 1) { // center - patch1
                loopVars.U += uBs1[dist];
                for (unsigned short i: {0, 1, 2}) {
                    const double modulus = fBs1[dist]*siteSiteSeparation[j][i];
                    loopVars.force[firstIPC][i] -= modulus;
                    loopVars.force[secndIPC+nIPCs][i] += modulus;
                }
            } else if (j == 3) { // patch1 - center
                loopVars.U += uBs1[dist];
                for (unsigned short i: {0, 1, 2}) {
                    const double modulus = fBs1[dist]*siteSiteSeparation[j][i];
                    loopVars.force[firstIPC][i] -= modulus;
                    loopVars.force[secndIPC+nIPCs][i] += modulus;
                }
            } else if (j == 2) { // center - patch2
                loopVars.U += uBs2[dist];
                for (unsigned short i: {0, 1, 2}) {
                    const double modulus = fBs2[dist]*siteSiteSeparation[j][i];
                    loopVars.force[firstIPC][i] -= modulus;
                    loopVars.force[secndIPC+nIPCs+nIPCs][i] += modulus;
                }
            } else if (j == 6) { // patch2 - center
                loopVars.U += uBs2[dist];
                for (unsigned short i: {0, 1, 2}) {
                    const double modulus = fBs2[dist]*siteSiteSeparation[j][i];
                    loopVars.force[firstIPC+nIPCs+nIPCs][i] -= modulus;
                    loopVars.force[secndIPC][i] += modulus;
                }
            } else if (j == 4) { // patch1 - patch1
                loopVars.U += us1s1[dist];
                for (unsigned short i: {0, 1, 2}) {
                    const double modulus = fs1s1[dist]*siteSiteSeparation[j][i];
                    loopVars.force[firstIPC+nIPCs][i] -= modulus;
                    loopVars.force[secndIPC+nIPCs][i] += modulus;
                }
            } else if (j == 8) { // patch2 - patch2
                loopVars.U += us2s2[dist];
                for (unsigned short i: {0, 1, 2}) {
                    const double modulus = fBB[dist]*siteSiteSeparation[j][i];
                    loopVars.force[firstIPC+nIPCs+nIPCs][i] -= modulus;
                    loopVars.force[secndIPC+nIPCs+nIPCs][i] += modulus;
                }
            } else if (j == 5) { // patch1 - patch2
                loopVars.U += us1s2[dist];
                for (unsigned short i: {0, 1, 2}) {
                    const double modulus = fs1s2[dist]*siteSiteSeparation[j][i];
                    loopVars.force[firstIPC][i] -= modulus;
                    loopVars.force[secndIPC][i] += modulus;
                }
            } else if (j == 7) { // patch2 - patch1
                loopVars.U += us1s2[dist];
                for (unsigned short i: {0, 1, 2}) {
                    const double modulus = fs1s2[dist]*siteSiteSeparation[j][i];
                    loopVars.force[firstIPC][i] -= modulus;
                    loopVars.force[secndIPC][i] += modulus;
                }
            }
        }
    }
}
