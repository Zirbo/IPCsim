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

#include <random>
#include <chrono>

class RandomNumberGenerator {
    std::mt19937 randomNumberGenerator;
public:
    RandomNumberGenerator() {
        randomNumberGenerator.seed(std::chrono::system_clock::now().time_since_epoch().count());
    }
    double getRandomDoubleInRange(double from, double to) {
        return std::uniform_real_distribution<double>(from, to)(randomNumberGenerator);
    }
    double getRandom55() {
        return getRandomDoubleInRange(-.5, .5);
    }
    double getRandom11() {
        return getRandomDoubleInRange(-1., 1.);
    }
};

#endif //__ZILRANDOM_HEADER_INCLUDED__
