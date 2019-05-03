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

#ifndef __ZILRANDOM_HEADER_INCLUDED__
#define __ZILRANDOM_HEADER_INCLUDED__

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


#endif //__ZILRANDOM_HEADER_INCLUDED__
