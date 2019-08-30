/********************************************************************************
 * 10 June 2014
 * EVN molecular dynamics for a system with potential from
 *        Bianchi, Kahl & Likos, Soft Matter 7, 8313 (2011).
 * The length unit is the HS diameter, the energy unit is given in the input file.
 * All the masses (both the CM and the patches) are 1, which means that the total mass is 3.
 * Will use cell lists
 * 
 *******************************************************************************/
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include "zilvo.hpp"


using namespace std;
using namespace space;


struct parametri
{
  // state point
  int nIPCs;
  double rho;
  double kT, kTimposed;
  // potential
  double e_BB, e_BS, e_SS, e_min;
  double bigRadius, smallRadius, eccentricity;
  double FakeHSexp, FakeHScoef;
  double FUsamplingDeltar, Tollerance;
  // simulation
  double dt_nonscaled;
  double SimLength;
  double PrintEvery;
  // work parameters
  double kToverK, E, U, K, L, L2, dt;
  double rmin2, rminbb, rminbs, rminss;
  double PotRange, PotRangeSquared;
  double PatchDistance, PatchDistanceSquared;
  double cP1, cP2, cPcm;
  int nPrints;
  int nPatc;
};
// force and potential tables computation
struct FU_table
{
    double *uBB, *uBS, *uSS, *fBB, *fBS, *fSS;
    void make_table(parametri par);
};
double omega(double Ra, double Rb, double rab);
double d_dr_omega(double Ra, double Rb, double rab);
// simulation parts
void free_force(parametri & par, vec x[], vec v[], vec F[], FU_table tab, cell_lists cells);
void verlet(unsigned long &t, parametri & par, vec x[], vec v[], vec F[], FU_table tab, cell_lists cells);
// selfexplanatory
void warmup(parametri & par, vec *& x, vec *& v, vec *& F, char *& farben, FU_table & tab, ostream & OUT, cell_lists & cells, bool restoreprevious);
void output(vec x[], vec v[], vec F[], char farben[], string nome, parametri par, unsigned long t, bool append);

int main ( int argc, char *argv[] )
{
  // declarations
  time_t program_begins, simulation_begins, simulation_ends;
  time(&program_begins);
  parametri par;
  vec *x, *v, *F;
  char *farben;
  FU_table tab;
  cell_lists cells;
  unsigned long t(0);
  fstream OUT;
  fstream ET;
  string s1="siml/trajectory.xyz";
  
  // clean up old data
  if(system("rm -r siml") != 0)  cerr<<"Just ignore this 'rm: cannot remove ...' thing\n";
  if(system("mkdir siml") != 0){ cerr<<"Unable to create 'siml/' directory...\n";    exit(1);  }
  
  // check starting flag; if ok initialize simulation
  OUT.open("siml/output.out", ios::out); 
  if(argc != 2)
  {    cerr<<"I only accept a single parameter, use either \"new\" or \"old\".\n";    exit(1);    }
  else if(string(argv[1])=="new")
    warmup(par,x,v,F,farben,tab,OUT,cells,false);
  else if(string(argv[1])=="old")
    warmup(par,x,v,F,farben,tab,OUT,cells,true);
  else
  {    cerr<<"Uncorrect parameter, use either \"new\" or \"old\".\n";    exit(1);    }
  
  
  // print starting configuration and initialize output file
  output(x,v,F,farben,s1,par,t,true);
  OUT<<"\nPlot evolution.out to check the evolution of the system.\n";
  ET.open("siml/evolution.out", ios::out);
  ET<<scientific<<setprecision(10);
  ET<<"#t\t\tT\t\tK\t\tU\t\tE\t\trminbb\t\trminbs\t\trminss\t\trmin\n";
  ET<<t*par.dt_nonscaled<<"\t"<<par.kT<<"\t"<<par.K/par.nIPCs<<"\t"<<par.U/par.nIPCs<<"\t"<<par.E/par.nIPCs;
  ET<<"\t"<<par.rminbb*par.L<<"\t"<<par.rminbs*par.L<<"\t"<<par.rminss*par.L<<"\t"<<sqrt(par.rmin2)*par.L<<endl;
  
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
    verlet(t,par,x,v,F,tab,cells);
    
    if(t%10==0)  // print time, temperature, energies, minimum distances
    {
      ET<<t*par.dt_nonscaled<<"\t"<<par.kT<<"\t"<<par.K/par.nIPCs<<"\t"<<par.U/par.nIPCs<<"\t"<<par.E/par.nIPCs;
      ET<<"\t"<<par.rminbb*par.L<<"\t"<<par.rminbs*par.L<<"\t"<<par.rminss*par.L<<"\t"<<sqrt(par.rmin2)*par.L<<endl;
    }
    
    if( t==printtimes[print])  // print configuration
    {
      output(x,v,F,farben,s1,par,t,true);
      print++;
    }
  }
  
  // check that total momentum is still zero and print final stuff
  vec pcm = vec( 0., 0., 0. );
  for(int i=0;i<3*par.nIPCs;i++)
    pcm += v[i];
  pcm *= par.L;
  
  output(x,v,F,farben,"startingstate.xyz",par,t,false);
  time(&simulation_ends);
  double dif = difftime (simulation_ends,simulation_begins);
  OUT<<"The simulation lasted "<<dif<<" seconds.\n";
  OUT<<"Residual momentum whole system = ( "<<pcm.x<<", "<<pcm.y<<", "<<pcm.z<<" ).\n"<<endl;
  dif = difftime (simulation_ends,program_begins);
  OUT<<"The whole program lasted "<<dif<<" seconds.\n";
  
  OUT<<"\nYou divided by zero and the last monk of the only non-fake religion\n";
  OUT<<"got Syphilis while lighting a candle for the candle-making industries.\n";
}





//************************************************************************//
double omega(double Ra, double Rb, double rab)
{  // BKL paper, formula 18
  if ( rab > Ra+Rb )    return 0.;
  else if ( rab <= fabs(Ra-Rb) )   return 8.*pow(min(Ra,Rb),3);
  else
  {
    double cacca = (Ra*Ra-Rb*Rb)/(2.*rab);
    return 2.*( (2.*Ra+cacca+rab/2.)*pow(Ra-cacca-rab/2.,2)
            + (2.*Rb-cacca+rab/2.)*pow(Rb+cacca-rab/2.,2) );
  }
}
double d_dr_omega(double Ra, double Rb, double rab)
{  // BKL paper, derivative of formula 18
  if ( rab >= Ra+Rb || rab <= fabs(Ra-Rb) )    return 0.;
  else
  {
    double cacca_plus, cacca_mins;
    cacca_plus = (Ra*Ra-Rb*Rb)/(2.*rab);    cacca_mins=cacca_plus-rab/2.;    cacca_plus+=rab/2.;
    return (6./rab)*( cacca_mins*(Ra-cacca_plus)*(Ra+cacca_plus) - cacca_plus*(Rb-cacca_mins)*(Rb+cacca_mins) );
  }
}
void FU_table::make_table(parametri par)
{
  int length = int( 2.*par.bigRadius/par.FUsamplingDeltar ) + 1;
  uBB = new double [length];
  uBS = new double [length];
  uSS = new double [length];
  fBB = new double [length];
  fBS = new double [length];
  fSS = new double [length];
  ofstream CACCA("siml/potentials.out"); CACCA<<scientific<<setprecision(6);
  CACCA<<"#r\t\tpotBB\t\tpotBS\t\tpotSS\t\tforBB\t\tforBS\t\tforSS\n";
  int conto(1);
  for ( int i=0; i<length; i++)
  {
    double r = i*par.FUsamplingDeltar;
    uBB[i] = (par.e_BB/par.e_min) * omega(par.bigRadius,par.bigRadius,r);
    uBS[i] = (par.e_BS/par.e_min) * omega(par.bigRadius,par.smallRadius,r);
    uSS[i] = (par.e_SS/par.e_min) * omega(par.smallRadius,par.smallRadius,r);
    fBB[i] = (par.e_BB/par.e_min) * d_dr_omega(par.bigRadius,par.bigRadius,r);
    fBS[i] = (par.e_BS/par.e_min) * d_dr_omega(par.bigRadius,par.smallRadius,r);
    fSS[i] = (par.e_SS/par.e_min) * d_dr_omega(par.smallRadius,par.smallRadius,r);
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
      CACCA<<r<<"\t"<<uBB[i]<<"\t"<<uBS[i]<<"\t"<<uSS[i]<<"\t"<<fBB[i]<<"\t"<<fBS[i]<<"\t"<<fSS[i]<<"\n";
    }
    // this is to not divide by r during the execution but not having it in the plots
    double x = 1./(r);
    fBB[i] *= x;
    fBS[i] *= x;
    fSS[i] *= x;
  }
  CACCA.close();
}





//************************************************************************//
void output(vec x[], vec v[], vec F[], char farben[], string nome, parametri par, unsigned long t, bool append)
{
  ofstream Out;
  if(append)   Out.open(nome.c_str(), ios::app);
  else         Out.open(nome.c_str());
  Out<<scientific<<setprecision(24);
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
void warmup(parametri & par, vec *& x, vec *& v, vec *& F, char *& farben, FU_table & tab, ostream & OUT, cell_lists & cells, bool restoreprevious)
{
  Ran rand(483248720420);
  int N1, N2, N3;
  // input from file
  fstream IN("input.in", ios::in);
  IN>>N1>>par.rho>>par.kTimposed;
  IN>>par.dt_nonscaled>>par.PrintEvery>>par.SimLength;
  IN>>par.e_BB>>par.e_BS>>par.e_SS>>par.e_min;
  IN>>par.eccentricity>>par.smallRadius;
  IN>>par.FakeHScoef>>par.FakeHSexp;
  IN>>par.FUsamplingDeltar>>par.Tollerance;
  IN.close();
  
  // processing the data
  N2 = N1*N1;        N3 = N2*N1;
  par.nIPCs = 4*N3;
  par.nPatc = 2*par.nIPCs;
  par.bigRadius = par.eccentricity + par.smallRadius;
  par.L=cbrt(par.nIPCs/par.rho);
  par.L2=par.L*par.L;
  par.kToverK = 2./(5.*par.nIPCs-3.);

  if(restoreprevious)
  {
    par.kTimposed = sqrt(par.kTimposed);
    char culo; double muntonarzu;
    IN.open("startingstate.xyz", ios::in);
    IN>>par.nIPCs>>muntonarzu;
    par.nIPCs /= 3;
    par.nPatc = 2*par.nIPCs;
    par.L = cbrt(par.nIPCs/par.rho);
    par.L2=par.L*par.L;
    par.kToverK = 2./(5.*par.nIPCs-3.);
    x = new vec[3*par.nIPCs]; v = new vec[3*par.nIPCs]; F = new vec[3*par.nIPCs];
    farben = new char[3*par.nIPCs];

    OUT<<"Reading "<<par.nIPCs<< "particles positions and velocities from file.\n";
    
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
  OUT<<par.e_BB<<"\t"<<par.e_BS<<"\t"<<par.e_SS<<"\t"<<par.e_min<<"\n";
  OUT<<par.eccentricity<<"\t"<<par.smallRadius<<"\n";
  OUT<<par.FakeHScoef<<"\t"<<par.FakeHSexp<<"\n";
  OUT<<par.FUsamplingDeltar<<"\t"<<par.Tollerance;
  
  
  OUT<<"\n*****************MD simulation in EVN ensemble for CGDH potential.********************\n";
  OUT<<"\nDensity = "<<par.nIPCs<<"/"<<pow(par.L,3)<<" = ";
  OUT<<par.nIPCs/pow(par.L,3)<<" = "<<par.rho<<"\nSide = "<<par.L<<endl;
  OUT<<"Total number of simulated particles: "<<par.nIPCs+par.nPatc<<endl;

  // potential sampling
  OUT<<"Printing potential plots in 'potentials.out'.\n";
  tab.make_table(par);
  
  // scaling of lenghts for [0.0:1.0] simulation box
  par.bigRadius /= par.L;
  par.smallRadius /= par.L;
  par.eccentricity /= par.L;
  par.dt = par.dt_nonscaled/par.L;
  par.FUsamplingDeltar /= par.L;
  par.PotRange = 2*par.bigRadius;
  par.PotRangeSquared = par.PotRange*par.PotRange;
  par.PatchDistance = 2*par.eccentricity;
  par.PatchDistanceSquared = par.PatchDistance*par.PatchDistance;
  par.cP1 = 5./6.;
  par.cPcm= 1./3.;
  par.cP2= 1./6.;
  
  // initialize positions of the particles
  if(!restoreprevious)
  {
    x = new vec[3*par.nIPCs]; v = new vec[3*par.nIPCs]; F = new vec[3*par.nIPCs];
    farben = new char[3*par.nIPCs];
  
    // scaling: sqrt(2kT/mPI) comes from boltzmann average of |v_x|
    // 4 is to compensate that the average |x| from x=rand.d55() is 0.25
    double vel_scaling = 4.*sqrt(2.*par.kTimposed/PI)/par.L;
    // initialize IPC positions
    for(int i=0;i<N3;i++)
    {
      // FCC is obtained as 4 intersecating SC
      x[i]      = vec( i%N1 + .1*rand.d55(), (i/N1)%N1 + .1*rand.d55(), i/N2 + .1*rand.d55() )/N1;
        floorccp(x[i]);
      x[i+N3]   = vec( .5+i%N1 + .1*rand.d55(), .5+(i/N1)%N1 + .1*rand.d55(), i/N2 + .1*rand.d55() )/N1;
        floorccp(x[i+N3]);
      x[i+2*N3] = vec( i%N1 + .1*rand.d55(), .5+(i/N1)%N1 + .1*rand.d55(), .5+i/N2 + .1*rand.d55() )/N1;
        floorccp(x[i+2*N3]);
      x[i+3*N3] = vec( .5+i%N1 + .1*rand.d55(), (i/N1)%N1 + .1*rand.d55(), .5+i/N2 + .1*rand.d55() )/N1;
        floorccp(x[i+3*N3]);
      // starting from random but ONLY TRANSLATIONAL speeds, for compatibility with rattle
      v[i]          = vec( rand.d55()*vel_scaling, rand.d55()*vel_scaling, rand.d55()*vel_scaling );
      v[i+N3]       = vec( rand.d55()*vel_scaling, rand.d55()*vel_scaling, rand.d55()*vel_scaling );
      v[i+N3+N3]    = vec( rand.d55()*vel_scaling, rand.d55()*vel_scaling, rand.d55()*vel_scaling );
      v[i+N3+N3+N3] = vec( rand.d55()*vel_scaling, rand.d55()*vel_scaling, rand.d55()*vel_scaling );
    }
    // initialize patches positions
    for(int i=0;i<par.nIPCs;i++)
    {
      vec a;        ranor(a,rand);    a *= par.eccentricity;
      
      x[i+par.nIPCs] = x[i] + a;
        floorccp(x[i+par.nIPCs]);
      x[i+par.nPatc] = x[i] - a;
        floorccp(x[i+par.nPatc]);
      
      v[i+par.nPatc] = v[i];
      v[i+par.nIPCs] = v[i];

      farben[i] = 'C';
      farben[i+par.nIPCs] = 'P';
      farben[i+par.nPatc] = 'P';
    }
  }
  
  // cell list compilation
  cells.initialize(1.,par.PotRange,par.nIPCs,x);
  OUT<<"Total number of cells: "<<cells.M3<<endl;
  cells.compilelists(x);
  
  // first computation of forces
  free_force(par,x,v,F,tab,cells);

  // check that total momentum is zero
  vec pcm = vec( 0., 0., 0. );
  for(int i=0;i<3*par.nIPCs;i++)
    pcm += v[i];
  pcm *= par.L;
  OUT<<"P whole system = ( "<<pcm.x<<", "<<pcm.y<<", "<<pcm.z<<" )."<<endl;
  
  // if not restoring, correct the total momentum to be zero
      if(!restoreprevious)
      {
        vec pcm_after = vec( 0., 0., 0. );
        pcm /= 3*par.nIPCs*par.L;
        for(int i=0;i<3*par.nIPCs;i++)
        {
          v[i] -= pcm;
          pcm_after += v[i];
        }
        pcm_after *= par.L;
        OUT<<"P whole system corrected = ( "<<pcm_after.x<<", "<<pcm_after.y<<", "<<pcm_after.z<<" )."<<endl;
      }
  
  // compute kinetic and total energy, temperature
  par.K = 0.;
  for(int i=0; i<3*par.nIPCs; i++)
  {
    par.K += v[i]*v[i];
  }
  par.K *= .5*par.L2;
  par.kT = par.K*par.kToverK;
  par.E = par.K + par.U;
}





/*****************************************************************************************/
void free_force(parametri & par, vec x[], vec v[], vec F[], FU_table tab, cell_lists cells)
{
  // Computes the force without accounting for constrains.
  // Force on i = sum over j of dU(r_ij)/dr * (x_j-x_i)/r_ij
  for(int i=0;i<3*par.nIPCs;i++)
    F[i] = vec(0.0,0.0,0.0);
  par.U = 0.0;  par.rmin2 = par.rminbb = par.rminbs = par.rminss = 1.;
  
  #pragma omp parallel
  {
    vec * Floop = new vec [3*par.nIPCs];
    double Uloop(0.), rmin2loop(1.), rminbbloop(1.), rminbsloop(1.), rminssloop(1.);
    for(int i=0;i<3*par.nIPCs;i++)    Floop[i] = vec(0.0,0.0,0.0);
    
    #pragma omp for
    for(int m=0; m<cells.M3; m++)  // loop over all cells
    {
      std::list<int> neighbs, local;
      cells.neighbour_cells(m,local,neighbs);
      
      for( list<int>::iterator loc = local.begin(); loc!=local.end(); loc++)
      { // loop over particles in neighbouring cells
        for( list<int>::iterator ext = neighbs.begin(); ext!=neighbs.end(); ext++)
        {
          vec ff, rji;        double r;
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
                  ff = rji*(tab.fBB[dist]);
                  Uloop += tab.uBB[dist];
                  // Store the smallest distance between IPCs:
                  if(r<rminbbloop) rminbbloop=r;
                }
                else if (i==0 || j==0)   // means cm -> patch
                {
                  ff = rji*(tab.fBS[dist]);
                  Uloop += tab.uBS[dist];
                  // Store the smallest distance between an IPC center and a patch:
                  if(r<rminbsloop) rminbsloop=r;
                }
                else                     // means patch -> patch
                {
                  ff = rji*(tab.fSS[dist]);
                  Uloop += tab.uSS[dist];
                  // Store the smallest distance between patches of different IPCs:
                  if(r<rminssloop) rminssloop=r;
                }
                Floop[*loc+i*par.nIPCs] -= ff;
                Floop[*ext+j*par.nIPCs] += ff;
              }
            }
          }
        }
        // loop over members of the cell m (internal interactions!)
        for( list<int>::iterator ins = loc; ins!=local.end(); ins++)
        {
          // starts from loc which is like summing over i >= j inside the cell
          // with a list, you have to access from loc because there's no way
          // to directly access (loc+1); so you need the following
          if(ins == loc) continue;  // to skip i=j iteration
          // after that, same code
          vec ff, rji;        double r;
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
                  ff = rji*(tab.fBB[dist]);
                  Uloop += tab.uBB[dist];
                  // Store the smallest distance between IPCs:
                  if(r<rminbbloop) rminbbloop=r;
                }
                else if (i==0 || j==0)   // means cm -> patch
                {
                  ff = rji*(tab.fBS[dist]);
                  Uloop += tab.uBS[dist];
                  // Store the smallest distance between an IPC center and a patch:
                  if(r<rminbsloop) rminbsloop=r;
                }
                else                     // means patch -> patch
                {
                  ff = rji*(tab.fSS[dist]);
                  Uloop += tab.uSS[dist];
                  // Store the smallest distance between patches of different IPCs:
                  if(r<rminssloop) rminssloop=r;
                }
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
}





void verlet(unsigned long &t, parametri & par, vec x[], vec v[], vec F[], FU_table tab, cell_lists cells)
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
   */
//   bool checkconostrains = false;
  
//   if(checkconostrains) cout<<t<<"\t"<<par.Tollerance*par.PatchDistanceSquared<<endl;
  for ( int i=0; i<par.nIPCs; i++ )
  {
    vec eFp1 = F[i+par.nIPCs]*par.cP1 + F[i]*par.cPcm - F[i+par.nPatc]*par.cP2;
    vec eFp2 = F[i+par.nPatc]*par.cP1 + F[i]*par.cPcm - F[i+par.nIPCs]*par.cP2;
    v[i+par.nIPCs] += eFp1*(.5*par.dt);
    v[i+par.nPatc] += eFp2*(.5*par.dt);
    vec r1 = x[i+par.nIPCs] + v[i+par.nIPCs]*par.dt;      floorccp(r1);
    vec r2 = x[i+par.nPatc] + v[i+par.nPatc]*par.dt;      floorccp(r2);
    
    // start working with constrains
    vec r12 = r1-r2;          lroundccp(r12);
    double diff = par.PatchDistanceSquared-r12*r12;
    while( fabs(diff) > par.Tollerance*par.PatchDistanceSquared )
    {
//       if(checkconostrains)  cout<<"HAD TO because "<<fabs(diff)<<">"<<par.Tollerance*par.PatchDistanceSquared<<"!\t";
      vec r12old = x[i+par.nIPCs]-x[i+par.nPatc];      lroundccp(r12old);
      double g = diff/( (r12*r12old)*4.*par.dt );
      vec DX = r12old*g;
      v[i+par.nIPCs] += DX;
      v[i+par.nPatc] -= DX;
        DX *= par.dt;
      r1 += DX;      floorccp(r1);
      r2 -= DX;      floorccp(r2);
      r12 = r1-r2;   lroundccp(r12);
      
      diff = par.PatchDistanceSquared-r12*r12;
    }
    x[i+par.nIPCs] = r1;
    x[i+par.nPatc] = r2;
    x[i] = r2 + r12*.5;  // r12 with their notation goes from 2 to 1!
    floorccp(x[i]);
//     if(checkconostrains)  cout<<diff<<"\t";
  }
//   if(checkconostrains)  cout<<endl;
  
  cells.compilelists(x);                // rewrite cell lists for the new iteration
  free_force(par,x,v,F,tab,cells);      // compute F(x[t+dt]) and the potential
  
  par.K = 0.;
  double virialconstrains(0.);

  bool checkVconostrains = false;
  
//   if(checkVconostrains) cout<<t<<"\t"<<par.Tollerance<<endl;
  for ( int i=0; i<par.nIPCs; i++ )
  {
    vec eFp1 = F[i+par.nIPCs]*par.cP1 + F[i]*par.cPcm - F[i+par.nPatc]*par.cP2;
    vec eFp2 = F[i+par.nPatc]*par.cP1 + F[i]*par.cPcm - F[i+par.nIPCs]*par.cP2;
    vec q1   = v[i+par.nIPCs] + eFp1*(.5*par.dt);
    vec q2   = v[i+par.nPatc] + eFp2*(.5*par.dt);
    
    // start working with constrains
    vec r12 = x[i+par.nIPCs]-x[i+par.nPatc];      lroundccp(r12);
    vec q12 = q1-q2;
    double k = (r12*q12)/( 2.*par.PatchDistanceSquared );
    while( fabs(k) > par.Tollerance )
    {
//       if(checkVconostrains)  cout<<"HAD TO because "<<fabs(k)<<">"<<par.Tollerance<<"!\t";
      vec DX = r12*k;
      q1 -= DX;
      q2 += DX;
      q12 = q1-q2;
      k = (r12*q12)/( 2.*par.PatchDistanceSquared );
      virialconstrains += k;
    }
    virialconstrains += k;
    v[i+par.nIPCs] = q1;
    v[i+par.nPatc] = q2;
    v[i] = (q1+q2)*.5;
    par.K += q1*q1+q2*q2+v[i]*v[i];
//     if(checkVconostrains)  cout<<k<<"\t";
  }
//   if(checkVconostrains)  cout<<endl;
  virialconstrains *= par.PatchDistance;
  par.K *= .5*par.L2;
  par.E = par.K + par.U;
  par.kT = par.kToverK*par.K;
  t++;
}
