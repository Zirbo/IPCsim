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

#ifndef __CELL_LISTS_HEADER_INCLUDED__
#define __CELL_LISTS_HEADER_INCLUDED__

//#include <iostream>
#include <list>
#include "zilvectors.hpp"

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

#endif //__CELL_LISTS_HEADER_INCLUDED__