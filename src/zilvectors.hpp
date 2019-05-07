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

#ifndef __ZILVECTORS_HEADER_INCLUDED__
#define __ZILVECTORS_HEADER_INCLUDED__

#include <cmath>
#include <list>
#include "zilrandom.hpp"

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
 * MISCELLANEOUS
 *----------------------------------------------------------------*/
void ranor(space::vec & a, RandomDoubleGenerator & r);
// 3D boundary conditions enforcers
inline void floorccp(space::vec & a)  {  a.x-=std::floor(a.x);   a.y-=std::floor(a.y);   a.z-=std::floor(a.z);   }
inline void lroundccp(space::vec & a) {  a.x-=std::lround(a.x);  a.y-=std::lround(a.y);  a.z-=std::lround(a.z);  }
// 2D boundary conditions enforcers
inline void floorccp(plane::vec & a)  {  a.x-=std::floor(a.x);   a.y-=std::floor(a.y);   }
inline void lroundccp(plane::vec & a) {  a.x-=std::lround(a.x);  a.y-=std::lround(a.y);  }
// Stores in 'a' a 3D random unit vector with the (I suppose!) Marsaglia algorithm
inline void ranor(space::vec & a, RandomDoubleGenerator & r)
{
  double x,y,quad=2.;
  while ( quad > 1. )  {    x = r.getRandom11();    y = r.getRandom11();    quad = x*x + y*y;  }
  double norm = 2.*sqrt(1.-quad);  a.x=x*norm;  a.y=y*norm;  a.z=1.-2.*quad;
}


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


#endif //__ZILVECTORS_HEADER_INCLUDED__
