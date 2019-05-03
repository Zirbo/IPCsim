#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include "IPCsimulator.hpp"
#define PI 3.1415926535897932

void IPCsimulator::run(bool doWarmup) {
    // declarations
    time_t program_begins, simulation_begins, simulation_ends;
    time(&program_begins);
    space::vec *x, *v, *F;
    char *farben;
    cell_lists cells;
    unsigned long t(0);
    std::fstream OUT;
    std::fstream ET;
    std::string s1="siml/trajectory.xyz";

    // clean up old data
    if(system("rm -r siml") != 0)  std::cerr<<"Just ignore this 'rm: cannot remove ...' thing\n";
    if(system("mkdir siml") != 0){ std::cerr<<"Unable to create 'siml/' directory...\n";    exit(1);  }

    // check starting flag; if ok initialize simulation
    OUT.open("siml/output.out", std::ios::out);
    warmup(x,v,F,farben,OUT,cells,doWarmup);


    // print starting configuration and initialize output file
    output(x,v,F,farben,s1,t,true);
    OUT<<"\nPlot evolution.out to check the evolution of the system.\n";
    ET.open("siml/evolution.out", std::ios::out);
    ET<<std::scientific<<std::setprecision(10);
    ET<<"#t\t\tT\t\tK\t\tU\t\tE\t\trminbb\t\trminbs\t\trminss\t\trmin\n";
    ET<<t*par.dt_nonscaled<<"\t"<<par.kT<<"\t"<<par.K/par.nIPCs<<"\t"<<par.U/par.nIPCs<<"\t"<<par.E/par.nIPCs;
    ET<<"\t"<<par.rminbb*par.L<<"\t"<<par.rminbs*par.L<<"\t"<<par.rminss*par.L<<"\t"<<sqrt(par.rmin2)*par.L<<std::endl;

    // compute printing schedule  --> CHANGE THIS AS NEEDED!
    //   par.nPrints = 2*int(log(par.SimLength)/log(2)+0.5);
    //   int *printtimes = new int [par.nPrints];
    //   printtimes[0]=1;    printtimes[1]=2;    printtimes[2]=3;    printtimes[3]=4;
    //   for(int i=4; i<par.nPrints; i++)  printtimes[i] = 2*printtimes[i-2];
    par.nPrints = int(par.SimLength/par.PrintEvery);
    int *printtimes = new int [par.nPrints];
    for(int i=0; i<par.nPrints; i++)
      printtimes[i] = int((i+1)*par.PrintEvery/par.dt_nonscaled);

    // simulation begins
    time(&simulation_begins);
    int print=0;
    while(t<int(par.SimLength/par.dt_nonscaled))
    {
      verlet(t,x,v,F,cells);

      if(t%10==0)  // print time, temperature, energies, minimum distances
      {
        ET<<t*par.dt_nonscaled<<"\t"<<par.kT<<"\t"<<par.K/par.nIPCs<<"\t"<<par.U/par.nIPCs<<"\t"<<par.E/par.nIPCs;
        ET<<"\t"<<par.rminbb*par.L<<"\t"<<par.rminbs*par.L<<"\t"<<par.rminss*par.L<<"\t"<<sqrt(par.rmin2)*par.L<<std::endl;
      }

      if( t==printtimes[print])  // print configuration
      {
        output(x,v,F,farben,s1,t,true);
        print++;
      }
    }

    // check that total momentum is still zero and print final stuff
    space::vec pcm = space::vec( 0., 0., 0. );
    for(int i=0;i<3*par.nIPCs;i++)
      pcm += v[i];
    pcm *= par.L;

    output(x,v,F,farben,"startingstate.xyz",t,false);
    time(&simulation_ends);
    double dif = difftime (simulation_ends,simulation_begins);
    OUT<<"The simulation lasted "<<dif<<" seconds.\n";
    OUT<<"Residual momentum whole system = ( "<<pcm.x<<", "<<pcm.y<<", "<<pcm.z<<" ).\n"<<std::endl;
    dif = difftime (simulation_ends,program_begins);
    OUT<<"The whole program lasted "<<dif<<" seconds.\n";

    OUT<<"\nYou divided by zero and the last monk of the only non-fake religion\n";
    OUT<<"got Syphilis while lighting a candle for the candle-making industries.\n";
}





//************************************************************************//
double IPCsimulator::omega(double Ra, double Rb, double rab)
{  // BKL paper, formula 18
  if ( rab > Ra+Rb )    return 0.;
  else if ( rab <= fabs(Ra-Rb) )   return 8.*pow(std::min(Ra,Rb),3);
  else
  {
    double cacca = (Ra*Ra-Rb*Rb)/(2.*rab);
    return 2.*( (2.*Ra+cacca+rab/2.)*pow(Ra-cacca-rab/2.,2)
            + (2.*Rb-cacca+rab/2.)*pow(Rb+cacca-rab/2.,2) );
  }
}
double IPCsimulator::d_dr_omega(double Ra, double Rb, double rab)
{  // BKL paper, derivative of formula 18
  if ( rab >= Ra+Rb || rab <= fabs(Ra-Rb) )    return 0.;
  else
  {
    double cacca_plus, cacca_mins;
    cacca_plus = (Ra*Ra-Rb*Rb)/(2.*rab);    cacca_mins=cacca_plus-rab/2.;    cacca_plus+=rab/2.;
    return (6./rab)*( cacca_mins*(Ra-cacca_plus)*(Ra+cacca_plus) - cacca_plus*(Rb-cacca_mins)*(Rb+cacca_mins) );
  }
}
void IPCsimulator::FU_table::make_table(Ensemble par)
{
  int length = int( 2.*par.bigRadius/par.FUsamplingDeltar ) + 1;
  uBB   = new double [length];
  uBs1  = new double [length];
  uBs2  = new double [length];
  us1s2 = new double [length];
  us1s1 = new double [length];
  us2s2 = new double [length];
  fBB   = new double [length];
  fBs1  = new double [length];
  fBs2  = new double [length];
  fs1s2 = new double [length];
  fs1s1 = new double [length];
  fs2s2 = new double [length];
  std::ofstream CACCA("siml/potentials.out"); CACCA<<std::scientific<<std::setprecision(6);
  CACCA<<"#r\t\t\tpotBB\t\t\tpotBs1\t\t\tpotBs2\t\t\tpots1s2\t\t\tpots2s2\t\t\tpots1s1";
  CACCA<<  "\t\t\tforBB\t\t\tforBs1\t\t\tforBs2\t\t\tfors1s2\t\t\tfors2s2\t\t\tfors1s1\n";
  int conto(1);
  for ( int i=0; i<length; i++)
  {
    double r = i*par.FUsamplingDeltar;
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
      double rm = pow(r,-par.FakeHSexp);
      uBB[i]   += par.FakeHScoef*((rm-2.)*rm+1.);
      fBB[i]   += -2.*par.FakeHSexp*par.FakeHScoef*(rm-1.)*rm/r;
    }
    if ( int( (1000.*i)/length ) == conto )
    {
      conto++;
      CACCA<<r<<"\t"<<uBB[i]<<"\t"<<uBs1[i]<<"\t"<<uBs2[i]<<"\t"<<us1s2[i]<<"\t"<<us2s2[i]<<"\t"<<us1s1[i]<<"\t";
      CACCA<<r<<"\t"<<fBB[i]<<"\t"<<fBs1[i]<<"\t"<<fBs2[i]<<"\t"<<fs1s2[i]<<"\t"<<fs2s2[i]<<"\t"<<fs1s1[i]<<"\n";
    }
    // this is to not divide by r during the execution but not having it in the plots
    double x = 1./(r);
    fBB[i] *= x;
    fBs1[i] *= x;
    fBs2[i] *= x;
    fs1s2[i] *= x;
    fs1s1[i] *= x;
    fs2s2[i] *= x;
  }
  CACCA.close();
}





//************************************************************************//
void IPCsimulator::output(space::vec x[], space::vec v[], space::vec F[], char farben[], std::string nome, unsigned long t, bool append)
{
  std::ofstream Out;
  if(append)   Out.open(nome.c_str(), std::ios::app);
  else         Out.open(nome.c_str());
  Out<<std::scientific<<std::setprecision(24);
  Out<<3*par.nIPCs<<"\n"<<t*par.dt_nonscaled<<"\n";
  for(int i=0; i<par.nIPCs; i++)
  {
    for(int a=0; a<3; a++)
    {
      Out<<farben[i+a*par.nIPCs]<<"\t";
      Out<<x[i+a*par.nIPCs].x<<"\t"<<x[i+a*par.nIPCs].y<<"\t"<<x[i+a*par.nIPCs].z<<"\t";
      Out<<v[i+a*par.nIPCs].x<<"\t"<<v[i+a*par.nIPCs].y<<"\t"<<v[i+a*par.nIPCs].z<<"\n";
    }
  }
  Out.close();
}





/*****************************************************************************************/
void IPCsimulator::warmup(space::vec *& x, space::vec *& v, space::vec *& F, char *& farben, std::ostream & OUT, cell_lists & cells, bool restoreprevious)
{
  Ran rand(483248720420);
  int N1, N2, N3;
  // input from file
  std::fstream IN("input.in", std::ios::in);
  IN>>N1>>par.rho>>par.kTimposed;
  IN>>par.dt_nonscaled>>par.PrintEvery>>par.SimLength;
  IN>>par.e_BB>>par.e_Bs1>>par.e_Bs2;
  IN>>par.e_s1s1>>par.e_s2s2>>par.e_s1s2;
  IN>>par.e_min;
  IN>>par.ecc1>>par.s1Radius;
  IN>>par.ecc2>>par.s2Radius;
  IN>>par.m1>>par.m2>>par.mc;
  IN>>par.FakeHScoef>>par.FakeHSexp;
  IN>>par.FUsamplingDeltar>>par.Tollerance;
  IN>>par.Ec.x>>par.Ec.y>>par.Ec.z;
  IN>>par.qc>>par.qp1>>par.qp2;
  IN.close();

  // processing the data
  N2 = N1*N1;        N3 = N2*N1;
  par.nIPCs = 4*N3;
  par.nPatc = 2*par.nIPCs;
  if ( abs( (par.ecc1+par.s1Radius)-(par.ecc2+par.s2Radius) ) >= 1e-10 )
  {
    std::cerr<<par.ecc1<<"+"<<par.s1Radius<<"="<<par.ecc1+par.s1Radius<<"-";
    std::cerr<<par.ecc2<<"+"<<par.s2Radius<<"="<<par.ecc2+par.s2Radius<<"=\n";
    std::cerr<<(par.ecc1+par.s1Radius)-(par.ecc2+par.s2Radius)<<std::endl;
    std::cerr<<"eccentricities and radii are not consistent!\n"; exit(1);
  }
  par.bigRadius = par.ecc1 + par.s1Radius;
  par.L=cbrt(par.nIPCs/par.rho);
  par.L2=par.L*par.L;
  par.kToverK = 2./(5.*par.nIPCs-3.);

  if(restoreprevious)
  {
    par.kTimposed = sqrt(par.kTimposed);
    char culo; double muntonarzu;
    IN.open("startingstate.xyz", std::ios::in);
    IN>>par.nIPCs>>muntonarzu;
    par.nIPCs /= 3;
    par.nPatc = 2*par.nIPCs;
    par.L = cbrt(par.nIPCs/par.rho);
    par.L2=par.L*par.L;
    par.kToverK = 2./(5.*par.nIPCs-3.);
    x = new space::vec[3*par.nIPCs]; v = new space::vec[3*par.nIPCs]; F = new space::vec[3*par.nIPCs];
    farben = new char[3*par.nIPCs];

  //  OUT<<"Reading "<<par.nIPCs<< " particles positions and velocities from file.\n";

    for(int i=0;i<par.nIPCs;i++)
    {
      for(int a=0; a<3; a++)
      {
        IN>>culo>>x[i+a*par.nIPCs].x>>x[i+a*par.nIPCs].y>>x[i+a*par.nIPCs].z;
        IN>>v[i+a*par.nIPCs].x>>v[i+a*par.nIPCs].y>>v[i+a*par.nIPCs].z;
        v[i+a*par.nIPCs] *= par.kTimposed;
        farben[i+a*par.nIPCs]=culo;
      }
    }
    IN.close();
  }

  // output the data for future checks
  OUT<<N1<<"\t"<<par.rho<<"\t"<<par.kTimposed<<"\n";
  OUT<<par.dt_nonscaled<<"\t"<<par.PrintEvery<<"\t"<<par.SimLength<<"\n";
  OUT<<par.e_BB<<"\t"<<par.e_Bs1<<"\t"<<par.e_Bs2<<"\n";
  OUT<<par.e_s1s2<<"\t"<<par.e_s1s1<<"\t"<<par.e_s2s2<<"\n";
  OUT<<par.e_min<<"\n";
  OUT<<par.ecc1<<"\t"<<par.s1Radius<<"\n";
  OUT<<par.ecc2<<"\t"<<par.s2Radius<<"\n";
  OUT<<par.m1<<"\t"<<par.m2<<"\t"<<par.mc;
  OUT<<par.FakeHScoef<<"\t"<<par.FakeHSexp<<"\n";
  OUT<<par.FUsamplingDeltar<<"\t"<<par.Tollerance;
  OUT<<par.Ec.x<<"\t"<<par.Ec.y<<"\t"<<par.Ec.z<<"\n";
  OUT<<par.qc<<"\t"<<par.qp1<<"\t"<<par.qp2<<"\n";

  // computing fields
  par.Ep1 = par.Ec*par.qp1;
  par.Ep2 = par.Ec*par.qp2;
  par.Ec *= par.qc;

  OUT<<"\n*****************MD simulation in EVN ensemble for CGDH potential.********************\n";
  OUT<<"\nDensity = "<<par.nIPCs<<"/"<<pow(par.L,3)<<" = ";
  OUT<<par.nIPCs/pow(par.L,3)<<" = "<<par.rho<<"\nSide = "<<par.L<<std::endl;
  OUT<<"Total number of simulated particles: "<<par.nIPCs+par.nPatc<<std::endl;

  // potential sampling
  OUT<<"Printing potential plots in 'potentials.out'.\n";
  tab.make_table(par);

  // scaling of lenghts for [0.0:1.0] simulation box
  par.bigRadius /= par.L;
  par.s1Radius /= par.L;
  par.ecc1 /= par.L;
  par.s2Radius /= par.L;
  par.ecc2 /= par.L;
  par.dt = par.dt_nonscaled/par.L;
  par.FUsamplingDeltar /= par.L;
  par.PotRange = 2*par.bigRadius;
  par.PotRangeSquared = par.PotRange*par.PotRange;
  par.PatchDistance = par.ecc1+par.ecc2;
  par.PatchDistanceSquared = par.PatchDistance*par.PatchDistance;
  par.im1 = 1./par.m1;
  par.im2 = 1./par.m2;
  par.imc = 1./par.mc;
  // inverse of the I parameter from formulas!
  double iI = 1./(par.PatchDistanceSquared*par.imc + par.ecc1*par.ecc1*par.im2 + par.ecc2*par.ecc2*par.im1);
  par.cP11 = 1.-par.ecc2*par.ecc2*iI*par.im1;
  par.cP12 = -par.ecc1*par.ecc2*iI*par.im2;
  par.cP1c = par.PatchDistance*par.ecc2*iI*par.imc;
  par.cP21 = par.cP12;
  par.cP22 = 1.-par.ecc1*par.ecc1*iI*par.im2;
  par.cP2c = par.PatchDistance*par.ecc1*iI*par.imc;
  par.alpha_1 = 1. - par.ecc2*iI*(par.ecc2*par.im1-par.ecc1*par.im2);
  par.alpha_2 = 1. + par.ecc1*iI*(par.ecc2*par.im1-par.ecc1*par.im2);
  par.alpha_sum = par.alpha_1 + par.alpha_2;
  par.Ec  /= par.L;
  par.Ep1 /= par.L;
  par.Ep2 /= par.L;

  // initialize positions of the particles
  if(!restoreprevious)
  {
    x = new space::vec[3*par.nIPCs]; v = new space::vec[3*par.nIPCs]; F = new space::vec[3*par.nIPCs];
    farben = new char[3*par.nIPCs];

    // scaling: sqrt(2kT/mPI) comes from boltzmann average of |v_x|
    // 4 is to compensate that the average |x| from x=rand.d55() is 0.25
    double vel_scaling = 4.*sqrt(2.*par.kTimposed/PI)/par.L;
    // initialize IPC positions
    for(int i=0;i<N3;i++)
    {
      // FCC is obtained as 4 intersecating SC
      x[i]      = space::vec( i%N1 + .1*rand.d55(), (i/N1)%N1 + .1*rand.d55(), i/N2 + .1*rand.d55() )/N1;
        floorccp(x[i]);
      x[i+N3]   = space::vec( .5+i%N1 + .1*rand.d55(), .5+(i/N1)%N1 + .1*rand.d55(), i/N2 + .1*rand.d55() )/N1;
        floorccp(x[i+N3]);
      x[i+2*N3] = space::vec( i%N1 + .1*rand.d55(), .5+(i/N1)%N1 + .1*rand.d55(), .5+i/N2 + .1*rand.d55() )/N1;
        floorccp(x[i+2*N3]);
      x[i+3*N3] = space::vec( .5+i%N1 + .1*rand.d55(), (i/N1)%N1 + .1*rand.d55(), .5+i/N2 + .1*rand.d55() )/N1;
        floorccp(x[i+3*N3]);
      // starting from random but ONLY TRANSLATIONAL speeds, for compatibility with rattle
      v[i]          = space::vec( rand.d55()*vel_scaling, rand.d55()*vel_scaling, rand.d55()*vel_scaling );
      v[i+N3]       = space::vec( rand.d55()*vel_scaling, rand.d55()*vel_scaling, rand.d55()*vel_scaling );
      v[i+N3+N3]    = space::vec( rand.d55()*vel_scaling, rand.d55()*vel_scaling, rand.d55()*vel_scaling );
      v[i+N3+N3+N3] = space::vec( rand.d55()*vel_scaling, rand.d55()*vel_scaling, rand.d55()*vel_scaling );
    }
    //v[1]=v[2]=v[3]=space::vec(0.,0.,0.);
    // initialize patches positions
    for(int i=0;i<par.nIPCs;i++)
    {
      space::vec a;        ranor(a,rand);

      x[i+par.nIPCs] = x[i] + a*par.ecc1;
        floorccp(x[i+par.nIPCs]);
      x[i+par.nPatc] = x[i] - a*par.ecc2;
        floorccp(x[i+par.nPatc]);

      space::vec surra;    ranor(surra,rand);
      surra = surra*sqrt(.5*(v[i]*v[i])/(surra*surra));
      surra = surra - a*(surra*a); // now it's orthogonal!

      v[i+par.nPatc] = v[i] + surra;
      v[i+par.nIPCs] = v[i] - surra;

      farben[i] = 'C';
      farben[i+par.nIPCs] = 'P';
      farben[i+par.nPatc] = 'G';
    }
  }

  // cell list compilation
  cells.initialize(1.,par.PotRange,par.nIPCs,x);
  OUT<<"Total number of cells: "<<cells.M3<<std::endl;
  cells.compilelists(x);

  // first computation of forces
  free_force(x,v,F,cells);

  // check that total momentum is zero
  space::vec pcm = space::vec( 0., 0., 0. );
  for(int i=0;i<3*par.nIPCs;i++)
    pcm += v[i];
  pcm *= par.L;
  OUT<<"P whole system = ( "<<pcm.x<<", "<<pcm.y<<", "<<pcm.z<<" )."<<std::endl;

  // if not restoring, correct the total momentum to be zero
      if(!restoreprevious)
      {
        space::vec pcm_after = space::vec( 0., 0., 0. );
        pcm /= 3*par.nIPCs*par.L;
        for(int i=0;i<3*par.nIPCs;i++)
        {
          v[i] -= pcm;
          pcm_after += v[i];
        }
        pcm_after *= par.L;
        OUT<<"P whole system corrected = ( "<<pcm_after.x<<", "<<pcm_after.y<<", "<<pcm_after.z<<" )."<<std::endl;
      }

  // compute kinetic and total energy, temperature
  par.K = 0.;
  double kc(0.),k1(0.),k2(0.);
  for(int i=0; i<par.nIPCs; i++)
  {
    kc += v[i]*v[i];
    k1 += v[i+par.nIPCs]*v[i+par.nIPCs];
    k2 += v[i+par.nPatc]*v[i+par.nPatc];
  }
  par.K = (kc*par.mc + k1*par.m1 + k2*par.m2)*.5*par.L2;
  par.kT = par.K*par.kToverK;
  par.E = par.K + par.U;
}





/*****************************************************************************************/
void IPCsimulator::free_force(space::vec x[], space::vec v[], space::vec F[], cell_lists cells)
{
  // Computes the force without accounting for constrains.
  // Force on i = sum over j of dU(r_ij)/dr * (x_j-x_i)/r_ij
  for(int i=0;i<3*par.nIPCs;i++)
    F[i] = space::vec(0.0,0.0,0.0);
  par.U = 0.0;  par.rmin2 = par.rminbb = par.rminbs = par.rminss = 1.;

  #pragma omp parallel
  {
    space::vec * Floop = new space::vec [3*par.nIPCs];
    double Uloop(0.), rmin2loop(1.), rminbbloop(1.), rminbsloop(1.), rminssloop(1.);
    for(int i=0;i<3*par.nIPCs;i++)    Floop[i] = space::vec(0.0,0.0,0.0);

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
              rji = x[*loc+i*par.nIPCs]-x[*ext+j*par.nIPCs];        lroundccp(rji);
              r = rji*rji;
              // Store the smallest distance between attraction centers:
              if( r<rmin2loop) rmin2loop=r;
              if (r <= par.PotRangeSquared)
              {
                r=sqrt(r);                      int dist = int( r/par.FUsamplingDeltar );
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
                Floop[*loc+i*par.nIPCs] -= ff;
                Floop[*ext+j*par.nIPCs] += ff;
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
              rji = x[*loc+i*par.nIPCs]-x[*ins+j*par.nIPCs];        lroundccp(rji);
              r = rji*rji;
              // Store the smallest distance between attraction centers:
              if( r<rmin2loop) rmin2loop=r;
              if (r <= par.PotRangeSquared)
              {
                r=sqrt(r);                      int dist = int( r/par.FUsamplingDeltar );
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
                Floop[*loc+i*par.nIPCs] -= ff;
                Floop[*ins+j*par.nIPCs] += ff;
              }
            }
          }
        }
      }
    }
    #pragma omp critical
    {
      for(int i=0;i<3*par.nIPCs;i++)
        F[i] += Floop[i];
      par.U += Uloop;
      if(rmin2loop < par.rmin2) par.rmin2 = rmin2loop;
      if(rminbbloop < par.rminbb) par.rminbb = rminbbloop;
      if(rminbsloop < par.rminbs) par.rminbs = rminbsloop;
      if(rminssloop < par.rminss) par.rminss = rminssloop;
    }
    delete [] Floop;
  }
  for(int i=0;i<par.nIPCs;i++)
    F[i] += par.Ec;
  for(int i=par.nIPCs;i<par.nPatc;i++)
    F[i] += par.Ep1;
  for(int i=par.nPatc;i<=par.nIPCs+par.nPatc;i++)
    F[i] += par.Ep1;
}





void IPCsimulator::verlet(unsigned long &t, space::vec x[], space::vec v[], space::vec F[], cell_lists cells)
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
  for ( int i=0; i<par.nIPCs; i++ )
  {
    space::vec eFp1 = F[i+par.nIPCs]*par.cP11 + F[i+par.nPatc]*par.cP12 + F[i]*par.cP1c;
    space::vec eFp2 = F[i+par.nIPCs]*par.cP21 + F[i+par.nPatc]*par.cP22 + F[i]*par.cP2c;
    v[i+par.nIPCs] += eFp1*(.5*par.dt*par.im1);
    v[i+par.nPatc] += eFp2*(.5*par.dt*par.im2);
    space::vec r1 = x[i+par.nIPCs] + v[i+par.nIPCs]*par.dt;      floorccp(r1);
    space::vec r2 = x[i+par.nPatc] + v[i+par.nPatc]*par.dt;      floorccp(r2);

    // start working with constrains
    space::vec r12 = r1-r2;          lroundccp(r12);
    double diff = r12*r12-par.PatchDistanceSquared;
    while( fabs(diff) > par.Tollerance*par.PatchDistanceSquared )
    {
      space::vec r12old = x[i+par.nIPCs]-x[i+par.nPatc];      lroundccp(r12old);
      double g = diff/( 2.*(r12*r12old)*par.alpha_sum*par.dt );
      space::vec DX = r12old*g;
      v[i+par.nIPCs] -= DX*par.alpha_1;
      v[i+par.nPatc] += DX*par.alpha_2;
        DX *= par.dt;
      r1 -= DX*par.alpha_1;      floorccp(r1);
      r2 += DX*par.alpha_2;      floorccp(r2);
      r12 = r1-r2;   lroundccp(r12);

      diff = r12*r12-par.PatchDistanceSquared;
    }
    x[i+par.nIPCs] = r1;
    x[i+par.nPatc] = r2;
    x[i] = r2 + r12*par.ecc2/par.PatchDistance;  // with their notation r12=r1<-2
    floorccp(x[i]);
  }

  cells.compilelists(x);                // rewrite cell lists for the new iteration
  free_force(x,v,F,cells);      // compute F(x[t+dt]) and the potential

  par.K = 0.;
  double virialconstrains(0.);

  for ( int i=0; i<par.nIPCs; i++ )
  {
    space::vec eFp1 = F[i+par.nIPCs]*par.cP11 + F[i+par.nPatc]*par.cP12 + F[i]*par.cP1c;
    space::vec eFp2 = F[i+par.nIPCs]*par.cP21 + F[i+par.nPatc]*par.cP22 + F[i]*par.cP2c;
    space::vec q1   = v[i+par.nIPCs] + eFp1*(.5*par.dt*par.im1);
    space::vec q2   = v[i+par.nPatc] + eFp2*(.5*par.dt*par.im2);

    // start working with constrains
    space::vec r12 = x[i+par.nIPCs]-x[i+par.nPatc];      lroundccp(r12);
    space::vec q12 = q1-q2;
    double k = (r12*q12)/( par.alpha_sum*par.PatchDistanceSquared );
    while( fabs(k) > par.Tollerance )
    {
      space::vec DX = r12*k;
      q1 -= DX*par.alpha_1;
      q2 += DX*par.alpha_2;
      q12 = q1-q2;
      k = (r12*q12)/( par.alpha_sum*par.PatchDistanceSquared );
    }
    virialconstrains += k;
    v[i+par.nIPCs] = q1;
    v[i+par.nPatc] = q2;
    v[i] = (q1*par.ecc2+q2*par.ecc1)/par.PatchDistance;
    par.K += (q1*q1)*par.m1 + (q2*q2)*par.m2 + (v[i]*v[i])*par.mc;
  }
  virialconstrains *= par.PatchDistance;
  par.K *= .5*par.L2;
  par.E = par.K + par.U;
  par.kT = par.kToverK*par.K;
  t++;
}
