#ifndef __CELL_LISTS_HEADER_INCLUDED__
#define __CELL_LISTS_HEADER_INCLUDED__

#include <list>
#include <vector>
#include "IPC.hpp"

class cell_lists
{
  // cell lists
public:
  cell_lists() {}
  void initialize(double simulationBoxSide, double interactionRange, int howManyParticles);
  // BoxSide and IntRange have to be in the same units of the positions x, it is not supposed to have a 1-side box
  void compilelists(std::vector<IPC> const& ipcs); // Only takes the first Nparticles coordinates!

  int getNumberofCells() {
      return totalCells;
  }
  const std::list<int> & getIPCsInCell(int cell) {
      return list_of_neighbours[cell];
  }
  const std::list<int> getIPCsInNeighbouringCells(int cell);

private:
  int cellsPerSide, cellsPerSideSquared, totalCells, numberOfParticles;
  double cellSide;
  std::vector<std::list<int>> neighbouring_cells, list_of_neighbours;

  int cellNumberFromPosition(Particle const& ipc) {
      return int(ipc.x[0]/cellSide) + cellsPerSide*int(ipc.x[1]/cellSide) + cellsPerSideSquared*int(ipc.x[2]/cellSide);
  }
  int cellNumberFromPosition(int x, int y, int z) {
      return x + cellsPerSide*y + cellsPerSideSquared*z;
  }
};

#endif //__CELL_LISTS_HEADER_INCLUDED__
