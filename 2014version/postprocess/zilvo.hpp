/********************************************************************************
 * ZILVO'S LIBRARY
 * Last modified  14 July 2014
 * TABLE OF CONTENTS
 * - vector apparate and random number generator
 * - correlation functions classes
 * - PBC enforcers, cell lists
 *
 * Log:
 *  -- Apr-May 2014
 * extended the higher order components for g/y (and corrected dire errors!)
 *  -- Jan 2014
 * corrected a small bug that was making it crash sometimes
 *  -- Nov 2013
 * added a list default constructor
 *  -- July 2013
 * time correlations: F_self(k,t)
 *  -- June 2013
 * cell lists using linked lists
 * time correlations: F(k,t)
 *  -- May 2013
 * spatial correlations: particles are in a box of side 1 to save multiplications
 * in the boundary conditions enforcement, added the possibility to compute third
 * order cofficients with m=0
 *  -- Feb 2013
 * totally rewritten vectors, quickened the spatial correlation computations
 *  -- Summer 2012
 * vector + spatial correlations
 * 
 *******************************************************************************/

#pragma once
#ifndef __ZILVO_HPP_INCLUDED__
#define __ZILVO_HPP_INCLUDED__ 1
#define PI 3.1415926535897932

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <list>

class Ran
{
// random number generator from NR3 for 64 bit CPUs.
public:
  Ran(long j);     // use it to initialize inserting seed
  long i64();      // 64 bit int
  double d55();    // [-0.5, +0.5)
  double d11();     // [-1.0, +1.0)
  double d01();    // [ 0.0, +1.0)
  int i32();       // 32 bit int
private:
  long u,v,w;
};
/*----------------------------------------------------------------*/
class Ran32
{
// random number generator from NR3 for 32 bit CPUs.
public:
  Ran32(int j);
  int int32();       // 32 bit int
  double d55();      // [-0.5, 0.5]
  double d01();      // [0, 1]
  double fulld55();  // more random version
  double fulld01();  // more random version
private:
  int u,v,w1,w2;
};
/*----------------------------------------------------------------*/
namespace plane
{
  // 2D vector class
  class vec
  {
  public:
    double x, y;
    vec (void)  {};
    vec (double a, double b)  {  x=a;  y=b;  }
   ~vec ()  {}
    vec operator + (vec);
    vec operator - (vec);
    vec operator * (double);
    vec operator / (double);
    void operator += (vec);
    void operator -= (vec);
    void operator *= (double);
    void operator /= (double);
    double operator * (vec);  // scalar product
  };
}
/*----------------------------------------------------------------*/
namespace space
{
  // 3D vector class
  class vec
  {
  public:
    double x, y, z;
    vec (void)  {};
    vec (double a, double b, double c)    {  x=a;  y=b;  z=c;  }
   ~vec ()  {}
    vec operator + (vec);
    vec operator - (vec);
    vec operator * (double);
    vec operator / (double);
    void operator += (vec);
    void operator -= (vec);
    void operator *= (double);
    void operator /= (double);
    double operator * (vec);  // scalar product
    vec operator ^ (vec);     // vector (cross, wedge) product
  };
}
/*----------------------------------------------------------------*
 * SPATIAL CORRELATION FUNCTIONS
 *----------------------------------------------------------------*/
// All lengths has to be in natural units, only x[] works in units in which BoxSide=1.
class G
{
// Computes the g(1,2) coefficients in the spherical harmonics expansion.
public:
  G(int Bins01, double BoxSide, int Nparticles, std::ostream & OUT, bool DoThirdOrder);
  // Bins BETWEEN 0 and 1. Limit of sampling is BoxSide/2.
  void compute(space::vec x[], space::vec w[]);
  void print (std::ostream & OUT);
private:
  short order;   int N,Q;   double L,dr;   double **g;   unsigned long M;   bool ThirdOrder;
};

class Y_KF
{
// Computes the y(1,2) coefficients for a one patch KF potential in the spherical harmonics expansion.
public:
  Y_KF(int Bins01, double BoxSide, int Nparticles, std::ofstream & OUT, double limit, double hole_depth_over_kT, double lambdaa, double cos_theta, bool DoThirdOrder);
  // Bins BETWEEN 0 and 1. Limit of sampling is lambdaa+limit. The others are the potential defining parameters.
  void compute(space::vec x[], space::vec w[], Ran & rand);
  // The last is a random numbers generator, it has to be already initialized.
  void print (void);
private:
  short order;   int N,Q;   double L,dr,e_kT, lambda, Qo, rCore;   double **y;   bool ThirdOrder;
};

class S_k
{
// Computes the (spherical) static structure factor
public:
  S_k(double kmax, double BoxSide, int Nparticles);
  // Last sampling value of k.
  void compute(space::vec x[]);
  void print (void);
private:
  double dk,L;    int nmax, ntot, N;    double *EXP;   unsigned long M;
};

class F_k
{
// Computes the intermediate scattering function in a "local" way,
// storing just the rho(k) at every instant
public:
  F_k(double MaxK, double BoxSide, int Nparticles, int SimLength);
  // ranges: k in [0,K) where K = [MaxK*L/2TT]
  // simulation is supposed in [0,SimLength)
  // the output time interval is also between [0,SimLength),
  // point i is computed from (SimLength-i) samples
  void compute(space::vec x[]);
  void print (void);
private:
  double dk,L;    int N, kmax, k3d, Mmax, M, Tmax;
  std::complex<double> **rhomk;
  
};

class Fs_k
{
// Computes the intermediate scattering function, Full+Self parts.
// storing the SINGLE PARTICLE DENSITIES. It has HUGE memory requirements
// and can ask how much will it be but only after having called the initializer
public:
  Fs_k(int nK, double Kgrid[], double BoxSide, int Nparticles, int SimLength);
  // nK is the size of kgrid, which contains the k where you approximatively want
  // the F(k,t) to be computed (they will be compared to L/2TT).
  // Simulation is supposed in [0,SimLength), there will be a check at the end.
  // point i is computed from (SimLength-i) samples
  void sketch(std::ostream & OUT);
  void compute(space::vec x[]);
  void print (std::ostream & OUT);
private:
  double dk, L;    int N, kmax, *kgrid, Mmax, M;
  std::complex<double> ***rhoikm;
  
};

class NEW_Fkt
{
// Computes the intermediate scattering function, Full+Self parts.
// storing the SINGLE PARTICLE DENSITIES. It has HUGE memory requirements
// and can ask how much will it be but only after having called the initializer
public:
  NEW_Fkt(double Kmax, double BoxSide, int Nparticles, int SimLength);
  // Size of the kgrid will be 2*int(Kmax*L/2TT)+1, many k will be coarsegrained.
  // Simulation is supposed in [0,SimLength), there will be a check at the end.
  // point i is computed from (SimLength-i) samples
  void sketch(std::ostream & OUT);
  void compute(space::vec x[]);
  void print (std::ostream & OUT);
private:
  double dk, L;    int N, kmax, kgrid, Mmax, M, *n, slots;
  std::complex<double> ***rhoikm;
  
};




/*----------------------------------------------------------------*
 * MISCELLANEOUS
 *----------------------------------------------------------------*/
void ranor(space::vec & a, Ran & r);
// 3D boundary conditions enforcers
inline void floorccp(space::vec & a)  {  a.x-=floor(a.x);   a.y-=floor(a.y);   a.z-=floor(a.z);   }
inline void lroundccp(space::vec & a) {  a.x-=lround(a.x);  a.y-=lround(a.y);  a.z-=lround(a.z);  }
// 2D boundary conditions enforcers
inline void floorccp(plane::vec & a)  {  a.x-=floor(a.x);   a.y-=floor(a.y);   }
inline void lroundccp(plane::vec & a) {  a.x-=lround(a.x);  a.y-=lround(a.y);  }
// Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
inline void ranor(space::vec & a, Ran & r)
{
  double x,y,quad=2.;
  while ( quad > 1. )  {    x = r.d11();    y = r.d11();    quad = x*x + y*y;  }
  double norm = 2.*sqrt(1.-quad);  a.x=x*norm;  a.y=y*norm;  a.z=1.-2.*quad;
}
class cell_lists
{
  // cell lists
public:
  cell_lists() {}
  cell_lists(double Side, double InteractionRange, int Nparticles, space::vec x[]);
  void initialize(double Side, double InteractionRange, int Nparticles, space::vec x[]);
  // BoxSide and IntRange have to be in the same units of the positions x,
  // it is not supposed to have a 1-side box
  void compilelists(space::vec x[]); // Only takes the first Nparticles coordinates!
  void neighbour_cells(int Cell, std::list<int> &local, std::list<int> &neigh);
  // Cell is the number of the cell you want to inquire, after the call
  // local will contain the indices of ALL the particles in cell Cell
  // neigh will contain the indices of ALL the particles in the neighbouring cells.
  int M3;
private:
  int M, M2, N;
  double l;
  std::list<int> *vicini, *lista;
  int cell(space::vec x);
  int cell(int x, int y, int z);
};









/*******************************************************************************/
// IMPLEMENTATIONS
/*******************************************************************************/

inline Ran::Ran(long j) : v(4101842887655102017LL), w(1)
{
  u = j ^ v; i64();
  v = u; i64();
  w = v; i64();
}
inline long Ran::i64()
{
  u = u * 2862933555777941757LL + 7046029254386353087LL;
  v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
  w = 4294957665U*(w & 0xffffffff) + (w >> 32);
  long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
  return (x + v) ^ w;
}
inline double Ran::d55() { return 5.42101086242752217E-20 * i64(); }
inline double Ran::d11() { return 1.084202172485504434E-19* i64(); }
inline double Ran::d01() { return d55()+0.5; } 
inline int Ran::i32() { return int( i64() ); }
/*******************************************************************************/
inline Ran32::Ran32(int j) : v(2244614371U), w1(521288629U), w2(362436069U)
{
  u = j ^ v; int32();
  v = u; int32();
}
inline int Ran32::int32()
{
  u = u * 2891336453U + 1640531513U;
  v ^= v >> 13; v ^= v << 17; v ^= v >> 5;
  w1 = 33378 * (w1 & 0xffff) + (w1 >> 16);
  w2 = 57225 * (w2 & 0xffff) + (w2 >> 16);
  int x = u ^ (u << 9); x ^= x >> 17; x ^= x << 6;
  int y = w1 ^ (w1 << 17); y ^= y >> 15; y ^= y << 5;
  return (x + v) ^ (y + w2);
}
inline double Ran32::d55() { return 2.32830643653869629E-10 * int32(); }
inline double Ran32::d01() { return 2.32830643653869629E-10 * int32()+0.5; }
inline double Ran32::fulld55()  {  return 2.32830643653869629E-10 * ( int32() + 2.32830643653869629E-10 * int32() );  }
inline double Ran32::fulld01()  {  return 2.32830643653869629E-10 * ( int32() + 2.32830643653869629E-10 * int32() ) + 0.5;  }
/*******************************************************************************/
inline plane::vec plane::vec::operator + (vec add)  {  return vec (x+add.x, y+add.y);  }
inline plane::vec plane::vec::operator - (vec add)  {  return vec (x-add.x, y-add.y);  }
inline plane::vec plane::vec::operator * (double a) {  return vec (x*a, y*a);  }
inline plane::vec plane::vec::operator / (double a) {  a=1./a; return vec (x*a, y*a);  }
inline void plane::vec::operator += (vec add) { x+=add.x; y+=add.y; }
inline void plane::vec::operator -= (vec add) { x-=add.x; y-=add.y; }
inline void plane::vec::operator *= (double a) {  x*=a; y*=a;  }
inline void plane::vec::operator /= (double a) {  a=1./a; x*=a; y*=a;  }
inline double plane::vec::operator * (vec a) {  return x*a.x + y*a.y;  }
/*******************************************************************************/
inline space::vec space::vec::operator + (vec add)  {  return vec (x+add.x, y+add.y, z+add.z);  }
inline space::vec space::vec::operator - (vec add)  {  return vec (x-add.x, y-add.y, z-add.z);  }
inline space::vec space::vec::operator * (double a) {  return vec (x*a, y*a, z*a);  }
inline space::vec space::vec::operator / (double a) {  a=1./a; return vec (x*a, y*a, z*a);  }
inline void space::vec::operator += (vec add) { x+=add.x; y+=add.y; z+=add.z; }
inline void space::vec::operator -= (vec add) { x-=add.x; y-=add.y; z-=add.z; }
inline void space::vec::operator *= (double a) {  x*=a; y*=a; z*=a;  }
inline void space::vec::operator /= (double a) {  a=1./a; x*=a; y*=a; z*=a;  }
inline double space::vec::operator * (vec a) {  return x*a.x + y*a.y + z*a.z;  }
inline space::vec space::vec::operator ^ (vec a)  {  return vec( y*a.z-z*a.y, z*a.x-x*a.z, x*a.y-y*a.x );  }
/*******************************************************************************/





/*******************************************************************************/
/*******************************************************************************/
/********************************************************************************
 * 3D correlation functions
 *******************************************************************************/


G::G(int Bins01, double BoxSide, int Nparticles, std::ostream & OUT, bool DoThirdOrder)
{
  // This computes g_ooo in the same way the isotropic g(r) is computed, and all the other
  // components like coeff*g_ooo*<spherical harmonics>. The theory can be found in the Hansen&McDonald.
  // The formula seems to come from Street&Tildesley, Proceedings of the Royal Society
  // of London. Series A, Mathematical and PhysicalSciences, Vol. 348, No. 1655 (Apr. 6, 1976), pp. 485-510
  
  ThirdOrder = DoThirdOrder;
  if (ThirdOrder) order = 30;
  else order = 14;
  L=BoxSide;   N=Nparticles;   dr=1./double(Bins01);
  M = 0;
  Q=int(    ( (sqrt(3.)/2)*L )  /dr   )+1;
  OUT<<"g(r): dr = "<<dr<<" numero suddivisioni = "<<Q<<" raggio massimo istogramma = "<<dr*Q<<std::endl;
  g = new double * [order];
  for( int i = 0 ; i < order ; i++ )
  {
    g[i] = new double [Q] ;
    for( int j = 0 ; j < Q ; j++ )
    {
      g[i][j]=0;
    }
  }
  dr/=L;
}
void G::compute(space::vec x[], space::vec w[])
{
  M++;
  for (int i=0; i<N-1; i++)
  {
    for (int k=i+1; k<N; k++)
    {
      space::vec rik = x[k] - x[i];
      lroundccp(rik);
      double r=sqrt(rik*rik);
      //questo è l'indice di distanza
      int R = int(r/dr);
      rik *= 1./r;
      double cosQ1 = rik*w[i];                //cos theta i
      double cosQ2 = rik*w[k];                //cos theta k
      double cosQ1_2 = cosQ1*cosQ1;           //cos theta i quadro
      double cosQ2_2 = cosQ2*cosQ2;           //cos theta k quadro
      //il seno tra 0 e 3.14 è sempre positivo quindi non mi serve fare distinzioni
      double senQ1_2 = 1.-cosQ1_2;
      double senQ2_2 = 1.-cosQ2_2;
      double senQ1 = sqrt(senQ1_2);
      double senQ2 = sqrt(senQ2_2);
      //come ottenere delta phi? Sottraggo alle orientazioni la componente parallela a r12;
      //il prodotto scalare dei due vettorini rimasti (normalizzati) è cos delta phi
      space::vec n1 = w[i] - rik*cosQ1;                n1 /= sqrt(n1*n1);
      space::vec n2 = w[k] - rik*cosQ2;                n2 /= sqrt(n2*n2);
      double cosp = n1*n2;
      //polinomi di Legendre pronti all'uso
      double P2Q1 = (3.*cosQ1_2-1.);
      double P2Q2 = (3.*cosQ2_2-1.);
      double cos2p = 2.*cosp*cosp-1.;

      ///////////////////////calcolo delle armoniche sferiche
      g[0][R]  += 1.;                                                 //Y00*Y00
      g[1][R]  += cosQ1;                                              //Y10*Y00
      g[2][R]  += cosQ2;                                              //Y00*Y10
      g[3][R]  += cosQ1*cosQ2;                                        //Y10*Y10
      g[4][R]  += senQ1*senQ2*cosp;                                   //Y11*Y1-1
      g[5][R]  += P2Q1;                                               //Y20*Y00
      g[6][R]  += P2Q1*cosQ2;                                         //Y20*Y10
      g[7][R]  += P2Q1*P2Q2;                                          //Y20*Y20
      g[8][R]  += cosQ1*P2Q2;                                         //Y10*Y20
      g[9][R]  += P2Q2;                                               //Y00*Y20
      g[10][R] += cosQ1*senQ1*senQ2*cosp;                             //Y21*Y1-1
      g[11][R] += cosQ1*cosQ2*senQ1*senQ2*cosp;                       //Y21*Y2-1
      g[12][R] += cosQ2*senQ1*senQ2*cosp;                             //Y11*Y2-1
      g[13][R] += senQ1_2*senQ2_2*cos2p;                              //Y22*Y2-2
    if (ThirdOrder)
    { //polinomi di Legendre pronti all'uso
      double P3Q1 = (5.*cosQ1_2-3.)*cosQ1;
      double P3Q2 = (5.*cosQ2_2-3.)*cosQ2;
      double cos3p = cosp*cos2p-sqrt((1.-cosp*cosp)*(1.-cos2p*cos2p));
      
      g[14][R] += P3Q1;                                               //Y30*Y00
      g[15][R] += P3Q1*cosQ2;                                         //Y30*Y10
      g[16][R] += P3Q1*P2Q2;                                          //Y30*Y20
      g[17][R] += P3Q1*P3Q2;                                          //Y30*Y30
      g[18][R] += P2Q1*P3Q2;                                          //Y20*Y30
      g[19][R] += cosQ1*P3Q2;                                         //Y10*Y30
      g[20][R] += P3Q2;                                               //Y00*Y30

      g[21][R] += senQ1*(5.*cosQ1_2-1.)*senQ2*cosp;                   //Y31*Y1-1
      g[22][R] += senQ1*(5.*cosQ1_2-1.)*senQ2*cosQ2*cosp;             //Y31*Y2-1
      g[23][R] += senQ1_2*cosQ1*senQ2_2*cos2p;                        //Y32*Y2-2
      g[24][R] += senQ1*(5.*cosQ1_2-1.)*senQ2*(5.*cosQ2_2-1.)*cosp;   //Y31*Y3-1
      g[25][R] += senQ1_2*cosQ1*senQ2_2*cosQ2*cos2p;                  //Y32*Y3-2
      g[26][R] += senQ1_2*senQ1*senQ2_2*senQ2*cos3p;                  //Y33*Y3-3
      g[27][R] += senQ1_2*senQ2_2*cosQ2*cos2p;                        //Y21*Y3-2
      g[28][R] += senQ1*cosQ1*senQ2*(5.*cosQ2_2-1.)*cosp;             //Y22*Y3-1
      g[29][R] += senQ1*senQ2*(5.*cosQ2_2-1.)*cosp;                   //Y11*Y3-1
    }
      ///////////////////////fine calcolo delle fottute armoniche sferiche
    }
  }
}
void G::print (std::ostream & OUT)
{
  dr*=L; //ripristina le unità naturali
  std::ofstream GG ("g");                GG<<std::scientific<<std::setprecision(8);
  double norm=(6.*pow(L,3))/(M*4.*N*PI*N*pow(dr,3));
  double integrale;
  for (int i=0; i<Q; i++)
  {
    double r       = (i+.5)*dr;
    double shell   = pow(double(i+1),3)-pow(double(i),3);
    double vol     = norm/ shell;
    g[0][i]        = g[0][i]*vol;
    g[1][i]        = sqrt(3.)*g[1][i]*vol*g[0][i];
    g[2][i]        = sqrt(3.)*g[2][i]*vol*g[0][i];
    g[3][i]        = 3.*g[3][i]*vol*g[0][i];
    g[4][i]        = -1.5*g[4][i]*vol*g[0][i];
    g[5][i]        = 0.5*sqrt(5.)*g[5][i]*vol*g[0][i];
    g[6][i]        = .5*sqrt(15)*g[6][i]*vol*g[0][i];
    g[7][i]        = 1.25*g[7][i]*vol*g[0][i];
    g[8][i]        = .5*sqrt(15)*g[8][i]*vol*g[0][i];
    g[9][i]        = 0.5*sqrt(5.)*g[9][i]*vol*g[0][i];
    g[10][i]       = -1.5*sqrt(5.)*g[10][i]*vol*g[0][i];
    g[11][i]       = -7.5*g[11][i]*vol*g[0][i];
    g[12][i]       = -1.5*sqrt(5.)*g[12][i]*vol*g[0][i];
    g[13][i]       = 1.875*g[13][i]*vol*g[0][i];
  if (ThirdOrder)
  {
    g[14][i]       = .5*sqrt(7.)*g[14][i]*vol*g[0][i];
    g[15][i]       = .5*sqrt(21.)*g[15][i]*vol*g[0][i];
    g[16][i]       = .25*sqrt(35.)*g[16][i]*vol*g[0][i];
    g[17][i]       =  1.75*g[17][i]*vol*g[0][i];
    g[18][i]       = .25*sqrt(35.)*g[18][i]*vol*g[0][i];
    g[19][i]       = .5*sqrt(21.)*g[19][i]*vol*g[0][i];
    g[20][i]       = .5*sqrt(7.)*g[20][i]*vol*g[0][i];

    g[21][i]       = -.75*sqrt(3.5)*g[21][i]*vol*g[0][i];
    g[22][i]       = -.75*sqrt(17.5)*g[22][i]*vol*g[0][i];
    g[23][i]       = 1.875*sqrt(7.)*g[23][i]*vol*g[0][i];
    g[24][i]       = -1.3125*g[24][i]*vol*g[0][i];
    g[25][i]       = 13.125*g[25][i]*vol*g[0][i];
    g[26][i]       = -2.1875*g[26][i]*vol*g[0][i];
    g[27][i]       = 1.875*sqrt(7.)*g[27][i]*vol*g[0][i];
    g[28][i]       = -.75*sqrt(17.5)*g[28][i]*vol*g[0][i];
    g[29][i]       =  -.75*sqrt(3.5)*g[29][i]*vol*g[0][i];
  }
    
    GG<<r;
    for( int j=0;j<order;j++ ) GG<<"\t"<<g[j][i];
    GG<<std::endl;
    
    integrale+=g[0][i]*shell;
  }
  GG.close();
  OUT<<"Integrale Ng000(r)/V = "<<(N/L/L/L)*integrale*( 4.*PI*pow(dr,3)/3. )<<" = "<<N-1<<" = N-1.\n";
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


Y_KF::Y_KF(int Bins01, double BoxSide, int Nparticles, std::ofstream & OUT, double limite, double hole_depth_over_kT, double lambdaa, double cos_theta, bool DoThirdOrder)
{
  // This computes the cavity function components like <(spherical harmonics)*exp(-betaDeltaPhi)> in a way similar you compute
  // the chemical potential with the Widom method, insterting a ghost particle and computing the mean of that thing. Of course
  // if the particle has any overlap it makes zero! The mothod comes from Henderson, Mol. Phys. 48, 1983, 2-389, but it is
  // unreadable; a much better explanation in LlanoRestrepo&Chapman, J. Chem. Phys. 97, 1002, 3-2046 
  ThirdOrder = DoThirdOrder;
  if (ThirdOrder) order = 30;
  else order = 14;
  e_kT = hole_depth_over_kT;
  lambda = lambdaa;
  Qo = cos_theta;
  dr = 1./double(Bins01);
  L = BoxSide;
  N = Nparticles;
  Q = int ( (lambda+limite)/dr   )+1;
  OUT<<"y(r): dr = "<<dr<<" numero suddivisioni = "<<Q<<" raggio massimo istogramma = "<<dr*Q<<std::endl;
  y = new double * [order];
  for( int i = 0 ; i < order ; i++ )
  {
    y[i] = new double [Q] ;
    for( int j = 0 ; j < Q ; j++ )
    {
      y[i][j]=0;
    }
  }
  dr/=L;
  lambda/=L;
  rCore=1./L;
}
  
void Y_KF::compute(space::vec x[], space::vec w[], Ran & rand)
{
  for (int sopra=0; sopra < N/5; sopra++)
  {
    int N1 = int(rand.d01()*N);                //sceglie la particella #1
    for (int R=0; R<Q; R++)        //cicla su tutte le distanze di campionamento della y
    {
      //ora dispongo la fantasma a distanza x*dr ed a caso sull'angolo solido.
      //rel è l'orientazione del vettore dalla fantasma alla N1; xo la posizione assoluta.
      space::vec xo, rel;                ranor(rel, rand);        xo = rel*( -dr*(0.5+double(R)) ) + x[N1];
      floorccp(xo);
      //orientazione della fantasma.
      space::vec wo;                ranor(wo, rand);
      double boltz=0;                //numero interazioni con le altre particelle
      for (int k=0; k<N; k++)
      {
        if ( k == N1 ) continue;                //nessuna interazione con la particella di riferimento
        space::vec rok = x[k] - xo;
        lroundccp(rok);
        //distanza tra la fantasma e la particella k:
        double r=sqrt(rok*rok);
        //vediamo se sono alla distanza giusta:
        if ( r < rCore )
        {
          boltz = -1;        //sovrapposizione! Allora non deve più calcolare nessuna intersezione
          break;
        }
        else if ( r < lambda )
        {        //distanza ok; vediamo se interagiscono
          rok *= 1./r;
          if ( rok*wo >= Qo && rok*w[k] <= -Qo )
            boltz++;        //c'è interazione
        }
      }
      //se boltz<0 vuol dire che c'è sovrap; exp(-beta infinito)= 0 e non vado a sommare niente a niente.
      //se boltz>0 non c'è sovrap, e ci sono boltz interazioni. il fattore è exp[-boltz*(-e/kT)]
      if ( boltz >= 0 )
      {
        boltz = exp(boltz*e_kT);
        //adesso ho a disposizione il fattore di boltzmann con cui pesare.
        //andiamo a calcolarci le medie pesate delle armoniche sferiche.
        double cosQ0 = rel*wo;                   //cos theta i
        double cosQ1 = rel*w[N1];                //cos theta k
        double cosQ0_2 = cosQ0*cosQ0;           //cos theta i quadro
        double cosQ1_2 = cosQ1*cosQ1;           //cos theta k quadro
        //il seno tra 0 e 3.14 è sempre positivo quindi non mi serve fare distinzioni
        double senQ0_2 = 1.-cosQ0_2;
        double senQ1_2 = 1.-cosQ1_2;
        double senQ0 = sqrt(senQ0_2);
        double senQ1 = sqrt(senQ1_2);
        //come ottenere delta phi? Sottraggo alle orientazioni la componente parallela a r12;
        //il prodotto scalare dei due vettorini rimasti (normalizzati) è cos delta phi
        space::vec n0 = wo - rel*cosQ0;                n0 = n0 / sqrt(n0*n0);
        space::vec n1 = w[N1] - rel*cosQ1;                n1 = n1 / sqrt(n1*n1);
        double cosp = n0*n1;
        //polinomi di Legendre pronti all'uso
        double P2Q0 = (3.*cosQ0_2-1.);
        double P2Q1 = (3.*cosQ1_2-1.);
        double cos2p = 2.*cosp*cosp-1.;

        ///////////////////////calcolo delle armoniche sferiche
        y[0][R]  += boltz;                                              //Y00*Y00
        y[1][R]  += cosQ0*boltz;                                        //Y10*Y00
        y[2][R]  += cosQ1*boltz;                                        //Y00*Y10
        y[3][R]  += cosQ0*cosQ1*boltz;                                  //Y10*Y10
        y[4][R]  += senQ0*senQ1*cosp*boltz;                             //Y11*Y1-1
        y[5][R]  += P2Q0*boltz;                                         //Y20*Y00
        y[6][R]  += P2Q0*cosQ1*boltz;                                   //Y20*Y10
        y[7][R]  += P2Q0*P2Q1*boltz;                                    //Y20*Y20
        y[8][R]  += cosQ0*P2Q1*boltz;                                   //Y10*Y20
        y[9][R]  += P2Q1*boltz;                                         //Y00*Y20
        y[10][R] += cosQ0*senQ0*senQ1*cosp*boltz;                       //Y21*Y1-1
        y[11][R] += cosQ0*cosQ1*senQ0*senQ1*cosp*boltz;                 //Y21*Y2-1
        y[12][R] += cosQ1*senQ1*senQ0*cosp*boltz;                       //Y11*Y2-1
        y[13][R] += senQ0_2*senQ1_2*cos2p*boltz;                        //Y22*Y2-2
      if (ThirdOrder)
      { //polinomi di Legendre pronti all'uso
        double P3Q0 = (5.*cosQ0_2-3.)*cosQ0;
        double P3Q1 = (5.*cosQ1_2-3.)*cosQ1;
        double cos3p = cosp*cos2p-sqrt((1.-cosp*cosp)*(1.-cos2p*cos2p));
        
        y[14][R] += P3Q0*boltz;                                         //Y30*Y00
        y[15][R] += P3Q0*cosQ1*boltz;                                   //Y30*Y10
        y[16][R] += P3Q0*P2Q1*boltz;                                    //Y30*Y20
        y[17][R] += P3Q0*P3Q1*boltz;                                    //Y30*Y30
        y[18][R] += P2Q0*P3Q1*boltz;                                    //Y20*Y30
        y[19][R] += cosQ0*P3Q1*boltz;                                   //Y10*Y30
        y[20][R] += P3Q1*boltz;                                         //Y00*Y30

        y[21][R] += senQ0*(5.*cosQ0_2-1.)*senQ1*cosp;                   //Y31*Y1-1
        y[22][R] += senQ0*(5.*cosQ0_2-1.)*senQ1*cosQ1*cosp;             //Y31*Y2-1
        y[23][R] += senQ0_2*cosQ0*senQ1_2*cos2p;                        //Y32*Y2-2
        y[24][R] += senQ0*(5.*cosQ0_2-1.)*senQ1*(5.*cosQ1_2-1.)*cosp;   //Y31*Y3-1
        y[25][R] += senQ0_2*cosQ0*senQ1_2*cosQ1*cos2p;                  //Y32*Y3-2
        y[26][R] += senQ0_2*senQ0*senQ1_2*senQ1*cos3p;                  //Y33*Y3-3
        y[27][R] += senQ0_2*senQ1_2*cosQ1*cos2p;                        //Y21*Y3-1
        y[28][R] += senQ0*cosQ0*senQ1*(5.*cosQ1_2-1.)*cosp;             //Y22*Y3-2
        y[29][R] += senQ0*senQ1*(5.*cosQ1_2-1.)*cosp;                   //Y11*Y3-1
      }
        ///////////////////////fine calcolo delle fottute armoniche sferiche
      }
    }
  }
}
  
void Y_KF::print (void)
{
  dr*=L;
  std::ofstream YY ("y");                YY<<std::scientific<<std::setprecision(8);
  for (int i=0; i<Q; i++)
  {
    y[0][i]        = y[0][i];
    y[1][i]        = sqrt(3.)*y[1][i];
    y[2][i]        = sqrt(3.)*y[2][i];
    y[3][i]        = 3.*y[3][i];
    y[4][i]        = -1.5*y[4][i];
    y[5][i]        = 0.5*sqrt(5.)*y[5][i];
    y[6][i]        = 3.*sqrt(2.5)*y[6][i];
    y[7][i]        = 1.25*y[7][i];
    y[8][i]        = .5*sqrt(15)*y[8][i];
    y[9][i]        = 0.5*sqrt(5.)*y[9][i];
    y[10][i]        = -1.5*sqrt(5.)*y[10][i];
    y[11][i]        = -7.5*y[11][i];
    y[12][i]        = -1.5*sqrt(5.)*y[12][i];
    y[13][i]        = 1.875*y[13][i];
  if (ThirdOrder)
  {
    y[14][i]       = .5*sqrt(7.)*y[14][i];
    y[15][i]       = .5*sqrt(21.)*y[15][i];
    y[16][i]       = .5*sqrt(35.)*y[16][i];
    y[17][i]       =  1.75*y[17][i];
    y[18][i]       = .5*sqrt(35.)*y[18][i];
    y[19][i]       = .5*sqrt(21.)*y[19][i];
    y[20][i]       = .5*sqrt(7.)*y[20][i];

    y[21][i]       = -.75*sqrt(3.5)*y[21][i];
    y[22][i]       = -.75*sqrt(17.5)*y[22][i];
    y[23][i]       = 1.875*sqrt(7.)*y[23][i];
    y[24][i]       = -1.3125*y[24][i];
    y[25][i]       = 13.125*y[25][i];
    y[26][i]       = -2.1875*y[26][i];
    y[27][i]       = 1.875*sqrt(7.)*y[27][i];
    y[28][i]       = -.75*sqrt(17.5)*y[28][i];
    y[29][i]       =  -.75*sqrt(3.5)*y[29][i];
  }
  
    YY<<(i+.5)*dr;
    for( int j=0;j<order;j++ ) YY<<"\t"<<y[j][i];
    YY<<std::endl;
  }
  YY.close();
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


S_k::S_k(double kmax, double BoxSide, int Nparticles)
{
  // Computes S(k) like the ensemble mean <| sum over j exp(-ik*x_j) |^2>.
  // I take the three vectors (dk,0,0),(0,dk,0),(0,0,dk) in k space and mean over them in the printing.
  // Every iteration I sample over every particles and for every particle I just compute the basic x*(dk,0,0);
  // all the other k space points are powers of this.
  // I could sample on many more k vectors but it's a bit more difficult and since near k=0 there's just no
  // way to it I do not find it useful (after all, you're more interested in the small k limit, you already
  // know that it's going to zero...). The theory is Hansen&McDonald.
  N = Nparticles;
  L = BoxSide;
  dk = 2*PI;
  nmax = int ( kmax*L/dk );
  ntot = 3*nmax+1;
  EXP = new double [ntot];
  for ( int i=0; i<ntot; i++ )        EXP[i]=0;
  M = 0;
}

void S_k::compute(space::vec x[])
{
  M++;
  std::complex<double> *temp = new std::complex<double> [ntot];
  for ( int i=0; i<ntot; i++ )        temp[i]=(0.,0.);
  for ( int i=0; i<N; i++ )
  {
    std::complex<double> expx ( cos(dk*x[i].x), -sin(dk*x[i].x) );
    std::complex<double> expy ( cos(dk*x[i].y), -sin(dk*x[i].y) );
    std::complex<double> expz ( cos(dk*x[i].z), -sin(dk*x[i].z) );
    
    temp[0] += std::complex<double>( 1., 0. );
    for ( int i=0; i<nmax; i++ )
    {
      temp[1+3*i] += pow(expx,i+1);
      temp[2+3*i] += pow(expy,i+1);
      temp[3+3*i] += pow(expz,i+1);
    }
  }
  for ( int i=0; i<ntot; i++ )
          EXP[i] += norm(temp[i]);
}

void S_k::print (void)
{
  dk/=L;
  std::ofstream AA ("S_k");                AA<<std::scientific<<std::setprecision(8);
  double norm = double(M*N);
  AA << 0. << "\t" << EXP[0]/norm << std::endl;
  for ( int i=0; i<nmax; i++ )
    AA << (i+1)*dk << "\t" << (EXP[1+3*i]+EXP[2+3*i]+EXP[3+3*i])/(3.*norm) << std::endl;
  AA.close();
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


F_k::F_k(double MaxK, double BoxSide, int Nparticles, int SimLength)
{
  // Almost same thing as S(k), but now you have to use the complete formula < rho(k,t) rho(-k,0) >/N
  // You store the rho(-k) for every configuration, and then you use it to compute
  // SimLenght different t values for F(k,t), every one with statistics i -> Mmax-i so the first
  // has Mmax and the last has just a configuration.
  // F(k,t,i) = rho(-k,0,i)*conj{rho(k,t,i)}             F(k,0,i) = |rho(-k,0,i)|^2
  // The theory is in the Hansen&McDonald.
  N = Nparticles;
  L = BoxSide;
  dk = 2*PI;
  kmax = int( MaxK*L/dk );
  k3d = 3*kmax;
  Mmax = SimLength;
  M = 0;

  rhomk = new std::complex<double> * [k3d+1];
  for( int k = 0 ; k <= k3d ; k++ )
  {
    rhomk[k] = new std::complex<double> [Mmax] ;
    for( int m = 0 ; m < Mmax ; m++ )
      rhomk[k][m] = std::complex<double> (0.,0.);
  }
}

void F_k::compute(space::vec x[])
{
  // compute rho(-k) for this configuration and store it.
  for ( int i=0; i<N; i++ )  //loop over all particles
  {
    // compute the three main exp( (dk,0,0)*x[i] ) for the three axes
    std::complex<double> expx ( cos(dk*x[i].x), sin(dk*x[i].x) );
    std::complex<double> expy ( cos(dk*x[i].y), sin(dk*x[i].y) );
    std::complex<double> expz ( cos(dk*x[i].z), sin(dk*x[i].z) );
    // now add to rho(k) the contribution of the three indipendent axes
    rhomk[0][M] += std::complex<double>(1.,0.);
    for ( int k=0; k<kmax; k++ )
    {
      rhomk[1+3*k][M] += pow(expx,k+1);
      rhomk[2+3*k][M] += pow(expy,k+1);
      rhomk[3+3*k][M] += pow(expz,k+1);
    }
  }
  M++;
}

void F_k::print (void)
{
  dk/=L;
  
  std::complex<double> ** F = new std::complex<double> * [kmax+1];
  for( int k = 0 ; k <= kmax ; k++ )
  {
    F[k] = new std::complex<double> [Mmax+1] ;
    for( int t = 0 ; t <= Mmax ; t++ )
    {
      F[k][t] = std::complex<double>(0.,0.);
      if (k==0)
      {
        for ( int i=0; i<Mmax-t; i++ )
        {
          F[0][t] += rhomk[0][i]*conj(rhomk[0][i+t]);
        }
        F[0][t] /= (Mmax-t)*N;
      }
      else
      {
        for ( int i=0; i<Mmax-t; i++ )
        {
          F[k][t] += rhomk[3*k-2][i]*conj(rhomk[3*k-2][i+t])
                   + rhomk[3*k-1][i]*conj(rhomk[3*k-1][i+t])
                   + rhomk[3*k  ][i]*conj(rhomk[3*k  ][i+t])  ;
        }
        F[k][t] /= 3.*(Mmax-t)*N;
      }
    }
  }
  // IMPRIMATUR!!! (cit. GP)
  std::ofstream AA ("F_tk");       AA<<std::scientific<<std::setprecision(8);     // kt -> kuaddu? ti...
  AA<<"#t\t\tk\t\tF.real\t\tF.imag\t\tF(k,t)/S(k)\n";
  for ( int t=0; t<Mmax; t++ )
    for ( int k=0; k<=kmax; k++ )
      AA << t << "\t" << k*dk << "\t" << F[k][t].real() << "\t" << F[k][t].imag() << "\t" << F[k][t].real()/F[k][0].real()  << "\n";
  AA.close();
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


Fs_k::Fs_k(int nK, double Kgrid[], double BoxSide, int Nparticles, int SimLength)
{
  /* A shitload of work. You store all the single particle k space densities so
   * you need a complex<double> array of length N x SimLenght x 3 x # of k points
   * pay real attention to the memory requirement or your fucking computer will
   * crush and I will fucking gut you, I've done this before (as Pat Bateman would say)
   * 
   * Just kidding!
   */
  N = Nparticles;
  L = BoxSide;
  dk = 2*PI;
  kmax = nK;
  kgrid = new int [kmax];
  for ( int i=0; i<kmax; i++ )
    kgrid[i] = int( Kgrid[i]*L/dk );
  Mmax = SimLength;
  M = 0;

  rhoikm = new std::complex<double> ** [N];
  for( int i = 0 ; i < N ; i++ )
  {
    rhoikm[i] = new std::complex<double> * [3*kmax];
    for( int k = 0 ; k < 3*kmax ; k++ )
    {
      rhoikm[i][k] = new std::complex<double> [Mmax] ;
      for( int m = 0 ; m < Mmax ; m++ )
        rhoikm[i][k][m] = std::complex<double> (0.,0.);
    }
  }
}

void Fs_k::sketch(std::ostream & OUT)
{
  OUT<<"\nYou are going to need "<<6.*sizeof(double)*N*(Mmax+1)*kmax/1048576.<<"MB = ";
  OUT<<6.*sizeof(double)*N*Mmax*kmax/1073741824.<<" GB of RAM for the one particle densities.\n"<<std::endl;
  
  OUT<<"True k grid:\n";
  for ( int k = 0 ; k < kmax ; k++ )
  {
    OUT << k <<"\t" << kgrid[k] <<"\t" << dk*double(kgrid[k])/L<<std::endl;
  }
}

void Fs_k::compute(space::vec x[])
{
  // compute exp(ik*x_i) for all the particles in this configuration and store everything.
  for ( int i=0; i<N; i++ )  //loop over all particles
  {
    // compute the three main exp( (dk,0,0)*x[i] ) for the three axes
    std::complex<double> expx ( cos(dk*x[i].x), sin(dk*x[i].x) );
    std::complex<double> expy ( cos(dk*x[i].y), sin(dk*x[i].y) );
    std::complex<double> expz ( cos(dk*x[i].z), sin(dk*x[i].z) );
    // now add to rho(i,k,m) the contribution of the three indipendent axes
//     rhoikm[i][0][M] += std::complex<double>(1.,0.); //this is a check
//     for ( int k=0; k<kmax; k++ )
//       rhoikm[i][k][M] = ( pow(expx,kgrid[k]) + pow(expy,kgrid[k]) + pow(expz,kgrid[k]) )/3.;
    for ( int k=0; k<kmax; k++ )
    {
      rhoikm[i][3*k  ][M] += pow(expx,kgrid[k]);
      rhoikm[i][3*k+1][M] += pow(expy,kgrid[k]);
      rhoikm[i][3*k+2][M] += pow(expz,kgrid[k]);
    }
  }
  M++;
}

void Fs_k::print (std::ostream & OUT)
{
  
  if (M != Mmax ) OUT<<" What the shit have you done? "<<M+1<<" is not "<<Mmax<<std::endl;
  
  dk/=L;
  
  std::ofstream AA ("Fs_tk");       AA<<std::scientific<<std::setprecision(8);
  AA<<"#t\tk\t\tF.real\t\tF.imag\t\tF(k,t)/F(k,0)\tFs.real\t\tFs.imag\t\tFs(k,t)/Fs(k,0)\n";
  double Fo[kmax], Fso[kmax];
  
  OUT<<"Computed k:\n";
  for( int k = 0 ; k < kmax ; k++ )
  {
    OUT<<k<<"/"<<kmax<<"="<<(100.*k)/kmax<<"%"<<std::endl;
    for( int t = 0 ; t < Mmax ; t++ )
    {
      std::complex<double> F  = std::complex<double> (0.,0.);
      std::complex<double> Fs = std::complex<double> (0.,0.);
      for ( int m=0; m<Mmax-t; m++ )
      {
        std::complex<double> rhokm1(0.,0.), rhokmt1(0.,0.);
        std::complex<double> rhokm2(0.,0.), rhokmt2(0.,0.);
        std::complex<double> rhokm3(0.,0.), rhokmt3(0.,0.);
        for ( int i=0; i<N; i++ )
        {
//           Fs     += rhoikm[i][k][m]*conj(rhoikm[i][k][m+t]);
//           rhokm  += rhoikm[i][k][m];
//           rhokmt += rhoikm[i][k][m+t];
          for (int j=0; j<3; j++)
            Fs     += rhoikm[i][j+3*k][m]*conj(rhoikm[i][j+3*k][m+t]);
          rhokm1  += rhoikm[i][3*k  ][m];
          rhokmt1 += rhoikm[i][3*k  ][m+t];
          rhokm2  += rhoikm[i][3*k+1][m];
          rhokmt2 += rhoikm[i][3*k+1][m+t];
          rhokm3  += rhoikm[i][3*k+2][m];
          rhokmt3 += rhoikm[i][3*k+2][m+t];
        }
//         F += rhokm*conj(rhokmt);
        F += rhokm1*conj(rhokmt1) + rhokm2*conj(rhokmt2) + rhokm3*conj(rhokmt3);
      }
//       Fs /= (Mmax-t+1)*N;
//       F  /= (Mmax-t+1)*N;
      Fs /= 3*(Mmax-t+1)*N;
      F  /= 3*(Mmax-t+1)*N;
      if (t == 0)
      {
        Fo[k]  = F.real();
        Fso[k] = Fs.real();
      }
      AA << t << "\t" << double(kgrid[k])*dk << "\t";
      AA << F.real()  << "\t" << F.imag()  << "\t" << F.real()/Fo[k]   << "\t";
      AA << Fs.real() << "\t" << Fs.imag() << "\t" << Fs.real()/Fso[k] << "\n";
    }
  }
  AA.close();
  
  //now deleting the huge array
  for( int i = 0 ; i < N ; i++ )
  {
    for( int k = 0 ; k < 3*kmax ; k++ )
      delete [] rhoikm[i][k];
    delete [] rhoikm[i];
  }
  delete [] rhoikm;
  // the end!
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


void cell_lists::initialize(double Side, double InteractionRange, int Nparticles, space::vec x[])
{
  N = Nparticles;
  M = int( Side/InteractionRange );
  M2 = M*M;    M3 = M2*M;
  l = Side/M;
  lista  = new std::list<int> [M3];     // every element is a list with the number of particles in the list
  vicini = new std::list<int> [M3];     // every element is a list with the number of nearest neighbour cells to that cell
  // fill the vicini list
  for(int x=0; x<M; x++)
  {
    for(int y=0; y<M; y++)
    {
      for(int z=0; z<M; z++)
      {
        int n      = cell(x,y,z);
        int xleft  = x-1;    if(xleft==-1)  xleft  =M-1;
        int xright = x+1;    if(xright==M)  xright =0;
        int yleft  = y-1;    if(yleft==-1)  yleft  =M-1;
        int yright = y+1;    if(yright==M)  yright =0;
        int ztop   = z+1;    if(ztop==M)    ztop   =0;
        vicini[n].push_back(cell( xright, y     , z   ));
        vicini[n].push_back(cell( xright, yright, z   ));
        vicini[n].push_back(cell( x     , yright, z   ));
        vicini[n].push_back(cell( xleft , yright, z   ));
        vicini[n].push_back(cell( xleft , yright, ztop));
        vicini[n].push_back(cell( x     , yright, ztop));
        vicini[n].push_back(cell( xright, yright, ztop));
        vicini[n].push_back(cell( xleft , y     , ztop));
        vicini[n].push_back(cell( x     , y     , ztop));
        vicini[n].push_back(cell( xright, y     , ztop));
        vicini[n].push_back(cell( xleft , yleft , ztop));
        vicini[n].push_back(cell( x     , yleft , ztop));
        vicini[n].push_back(cell( xright, yleft , ztop));
      }
    }
  }
  /* print neighbouring cells for debugging
  for(int i=0;i<M3;i++)
  {
    for(std::list<int>::iterator it = vicini[i].begin(); it!=vicini[i].end(); it++)
      std::cout<<*it<<"    ";
    std::cout<<std::endl;
  }*/
}
cell_lists::cell_lists(double Side, double InteractionRange, int Nparticles, space::vec x[])
{  // This is here only for backwards compatibility
  this->initialize(Side,InteractionRange,Nparticles,x);
}
inline int cell_lists::cell(space::vec x)
{  // gives you the list where this coordinate belong
  return int(x.x/l)+M*int(x.y/l)+M2*int(x.z/l);
}
inline int cell_lists::cell(int x, int y, int z)
{  // hope it really inlines because it's really stupid to have it
  return x + M*y + M2*z;
}
void cell_lists::compilelists(space::vec x[])
{  // empty all the lists
  for(int m=0; m<M3; m++)
    lista[m].clear();
  // put every particle in the right list
  for(int i=0; i<N; i++)
    lista[   cell(x[i])   ].push_back(i);
  // now the list contains the indices of the particles inside its volume
  
  // print out cell content for debugging
//   for(int m=0; m<M3; m++)
//   {
//     std::cout<<m<<":  ";
//     for(std::list<int>::iterator it = lista[m].begin(); it!=lista[m].end(); it++)
//       std::cout<<*it<<"  ";
//     std::cout<<std::endl;
//   }
  
}
void cell_lists::neighbour_cells(int Cell, std::list<int> &local, std::list<int> &neigh)
{
  local = lista[Cell];
  std::list<int> stampecullu;
  for(std::list<int>::iterator it = vicini[Cell].begin(); it!=vicini[Cell].end(); it++)
  {
//     std::cout<<*it<<std::endl;
    stampecullu = lista[*it];
    neigh.merge(stampecullu);
  }
}



/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/

NEW_Fkt::NEW_Fkt(double Kmax, double BoxSide, int Nparticles, int SimLength)
{
  // 30 Juni 2014
  // Takes all the possible points and coarse graines them. Can use a lot of memory.
  slots = 4;
  N = Nparticles;
  L = BoxSide;
  dk = 2*PI;
  kmax = int(Kmax*L/dk);
  kgrid = slots*kmax;
  Mmax = SimLength;
  M = 0;

  rhoikm = new std::complex<double> ** [N];
  for( int i = 0 ; i < N ; i++ )
  {
    rhoikm[i] = new std::complex<double> * [kgrid];
    for( int k = 0 ; k < kgrid ; k++ )
    {
      rhoikm[i][k] = new std::complex<double> [Mmax] ;
      for( int m = 0 ; m < Mmax ; m++ )
        rhoikm[i][k][m] = std::complex<double> (0.,0.);
    }
  }
  //std::cout<<kgrid<<" "<<kmax<<" "<<std::endl;
  n = new int [kgrid];
  for( int k = 0 ; k < kgrid ; k++ ) n[k]=0;
  for(int x=0;x<kmax;x++)
  {
    for(int y=0;y<kmax;y++)
    {
      for(int z=0;z<kmax;z++)
      {
        int k = x*x+y*y+z*z;
        k = ( k==0 ? 0 : int( slots*(sqrt(double(k))-1.) )+1 );
        if(k<kgrid)
        {
          n[k]++;
  //        std::cout<<x<<" "<<y<<" "<<z<<" "<<k<<"  "<<n[k]<<std::endl;
        }
      }
    }
  }
  //std::cout<<"Tua mamma!\n";
}

void NEW_Fkt::sketch(std::ostream & OUT)
{
  OUT<<"\nYou are going to need "<<2.*sizeof(double)*N*(Mmax)*kgrid/1048576.<<"MB = ";
  OUT<<2.*sizeof(double)*N*Mmax*kgrid/1073741824.<<" GB of RAM for the one particle densities.\n"<<std::endl;
//  for( int k = 0 ; k < kgrid ; k++ )
//  {
//    std::cout<<k<<"\t"<<n[k]<<std::endl;
//  }
}

void NEW_Fkt::compute(space::vec x[])
{
  // compute exp(ik*x_i) for all the particles in this configuration and store everything.
  for ( int i=0; i<N; i++ )  //loop over all particles
  {
    // compute the three main exp( (dk,0,0)*x[i] ) for the three axes
    std::complex<double> expx ( cos(dk*x[i].x), sin(dk*x[i].x) );
    std::complex<double> expy ( cos(dk*x[i].y), sin(dk*x[i].y) );
    std::complex<double> expz ( cos(dk*x[i].z), sin(dk*x[i].z) );
    // now add to rho(i,k,m) the contribution of the three indipendent axes
    for(int x=0;x<kmax;x++)
    {
      for(int y=0;y<kmax;y++)
      {
        for(int z=0;z<kmax;z++)
        {
          int k = x*x+y*y+z*z;
          k = ( k==0 ? 0 : int( slots*(sqrt(double(k))-1.) )+1 );
          if(k<kgrid)
            rhoikm[i][k][M] += pow(expx,x)*pow(expy,y)*pow(expz,z);
        }
      }
    }
  }
  M++;
}

void NEW_Fkt::print (std::ostream & OUT)
{
  
  if (M != Mmax ) OUT<<" What the shit have you done? "<<M+1<<" is not "<<Mmax<<std::endl;
  
  dk/=L;
  
  std::ofstream AA ("Fs_tk");       AA<<std::scientific<<std::setprecision(8);
  AA<<"#t\tk\t\tF.real\t\tF.imag\t\tF(k,t)/F(k,0)\tFs.real\t\tFs.imag\t\tFs(k,t)/Fs(k,0)\n";
  double Fo[kgrid], Fso[kgrid];
  
  OUT<<"Computed k:\n";
  for( int k = 0 ; k < kgrid ; k++ )
  {
    OUT<<k<<"/"<<kgrid<<"="<<(100.*k)/kgrid<<"%"<<std::endl;
    for( int t = 0 ; t < Mmax ; t++ )
    {
      std::complex<double> F  = std::complex<double> (0.,0.);
      std::complex<double> Fs = std::complex<double> (0.,0.);
      for ( int m=0; m<Mmax-t; m++ )
      {
        std::complex<double> rhokm(0.,0.), rhokmt(0.,0.);
        for ( int i=0; i<N; i++ )
        {
          Fs     += rhoikm[i][k][m]*conj(rhoikm[i][k][m+t]);
          rhokm  += rhoikm[i][k][m];
          rhokmt += rhoikm[i][k][m+t];
        }
        F += rhokm*conj(rhokmt);
      }
      if (n[k] != 0)
      {
        Fs /= (Mmax-t+1)*N*n[k];
        F  /= (Mmax-t+1)*N*n[k];
      }

      if (t == 0)
      {
        Fo[k]  = F.real();
        Fso[k] = Fs.real();
      }

      AA << t << "\t" << ( k==0 ? 0 : dk*(double(k-1)/slots+1.) ) << "\t";
      AA << F.real()  << "\t" << F.imag()  << "\t" << (Fo[k]!=0. ? F.real()/Fo[k] : 0. )  << "\t";
      AA << Fs.real() << "\t" << Fs.imag() << "\t" << (Fso[k]!=0. ? Fs.real()/Fso[k] : 0. ) << "\n";
    }
  }
  AA.close();
  
  //now deleting the huge array
  for( int i = 0 ; i < N ; i++ )
  {
    for( int k = 0 ; k < kgrid ; k++ )
      delete [] rhoikm[i][k];
    delete [] rhoikm[i];
  }
  delete [] rhoikm;
  // the end!
}


/*******************************************************************************/
/*******************************************************************************/
/*******************************************************************************/


class NEW_Sk
{
// Computes the (spherical) static structure factor
public:
  NEW_Sk(double Kmax, double BoxSide, int Nparticles);
  // Last sampling value of k.
  void compute(space::vec x[]);
  void print (void);
private:
  double dk,L;    int kmax, kgrid, N, slots, *n;
  double *EXP;    unsigned long M;
};


NEW_Sk::NEW_Sk(double Kmax, double BoxSide, int Nparticles)
{
  // Improved version, after the succes of NEW_Fkt
  slots = 4;
  N = Nparticles;
  L = BoxSide;
  dk = 2*PI;
  kmax = int(Kmax*L/dk);
  kgrid = slots*kmax;
  EXP = new double [kgrid];
  for ( int i=0; i<kgrid; i++ )        EXP[i]=0;
  M = 0;
  n = new int [kgrid];                                                              
  for( int k = 0 ; k < kgrid ; k++ ) n[k]=0;                                        
  for(int x=0;x<kmax;x++)                                                           
  {                                                                                 
    for(int y=0;y<kmax;y++)                                                         
    {                                                                               
      for(int z=0;z<kmax;z++)                                                       
      {                                                                             
        int k = x*x+y*y+z*z;                                                        
        k = int( slots*sqrt(double(k)) -.5 );                                       
        if(k>=kgrid)   continue;                                                    
        else                                                                        
        {                                                                           
          n[k]++;                                                                   
        }                                                                           
      }                                                                             
    }                                                                               
  }
}

void NEW_Sk::compute(space::vec x[])
{
  // compute exp(ik*x_i) for all the particles in this configuration and store everything.
  M++;
  std::complex<double> *temp = new std::complex<double> [kgrid];
  for ( int i=0; i<kgrid; i++ )        temp[i]=std::complex<double>(0.,0.);
  for ( int i=0; i<N; i++ )
  {
    // compute the three main exp( (dk,0,0)*x[i] ) for the three axes
    std::complex<double> expx ( cos(dk*x[i].x), -sin(dk*x[i].x) );
    std::complex<double> expy ( cos(dk*x[i].y), -sin(dk*x[i].y) );
    std::complex<double> expz ( cos(dk*x[i].z), -sin(dk*x[i].z) );
    // now add to rho(i,k,m) the contribution of the three indipendent axes
    for(int x=0;x<kmax;x++)                                                         
    {                                                                               
      for(int y=0;y<kmax;y++)                                                   
      {                                                                         
        for(int z=0;z<kmax;z++)                                                 
        {                                                                       
          int k = x*x+y*y+z*z;                                                  
          k = int( slots*sqrt(double(k)) -.5 );                                 
          if(k>=kgrid)   continue;                                              
          else  temp[k] += pow(expx,x)*pow(expy,y)*pow(expz,z);         
        }                                                                       
      }                                                                         
    }
  }
  for ( int k=0; k<kgrid; k++ )
          EXP[k] += norm(temp[k]);
}

void NEW_Sk::print (void)
{
  dk/=L;
  std::ofstream AA ("S_k");                AA<<std::scientific<<std::setprecision(8);
  for( int k = 0 ; k < kgrid ; k++ )
  {
    AA << ( k==0 ? dk : (.25+double(k)/slots)*dk  ) << "\t";
    AA << ( n[k]== 0 ? EXP[k] : EXP[k]/(M*N*n[k]) ) << std::endl;
  }
  AA.close();
}


#endif
