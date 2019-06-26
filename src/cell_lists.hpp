#ifndef __CELL_LISTS_HEADER_INCLUDED__
#define __CELL_LISTS_HEADER_INCLUDED__

#include <list>
#include <vector>
#include "IPC.hpp"

class cell_lists
{
public:
    cell_lists() {}
    // initialize --- simulationBoxSide and interactionRange have to be in the same units of length that are used by the ipcs vector passed in compileLists
    void initialize(double simulationBoxSide, double interactionRange, int howManyParticles);
    void compileLists(std::vector<IPC> const& ipcs);

    int getNumberofCells() {
        return totalCells;
    }
    // getIPCsInCell --- returns a list containing the indices of all the particles of the cell
    const std::list<int> & getIPCsInCell(int cell) {
        return list_of_neighbours[cell];
    }
    // getIPCsInNeighbouringCells --- returns a list containing the indices of all the particles in the neighbouring cells
    const std::list<int> getIPCsInNeighbouringCells(int cell);

private:
  int cellsPerSide, cellsPerSideSquared, totalCells, numberOfParticles;
  double cellSide;
  std::vector<std::list<int>> neighbouring_cells,   // element n is a list containing the indices of the particles in cell n
                              list_of_neighbours;   // element m is a list containing the indices of the cells that are nearest neighbour of cell m

  int cellNumberFromPosition(Particle const& ipc) {
      return int(ipc.x[0]/cellSide) + cellsPerSide*int(ipc.x[1]/cellSide) + cellsPerSideSquared*int(ipc.x[2]/cellSide);
  }
  int cellNumberFromPosition(int x, int y, int z) {
      return x + cellsPerSide*y + cellsPerSideSquared*z;
  }
};

#endif //__CELL_LISTS_HEADER_INCLUDED__
