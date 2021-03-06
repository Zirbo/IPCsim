#include "cell_lists.hpp"
#include <iostream>

void cell_lists::initialize(double simulationBoxSide, double interactionRange, int howManyParticles, SimulationType simulationType) {
    numberOfParticles = howManyParticles;
    cellsPerSide = int( simulationBoxSide/interactionRange);
    cellsPerSideSquared = cellsPerSide*cellsPerSide;
    totalCells = cellsPerSideSquared*cellsPerSide;
    cellSide = simulationBoxSide/cellsPerSide;
    if(cellsPerSide < 3) {
      std::cerr << "The system is too small (only " << totalCells << " cell lists) and the list system would be degenerate. Quitting.";
      exit(1);
    }
    list_of_neighbours.resize(totalCells);
    neighbouring_cells.resize(totalCells);
    // fill the list of neighbouring cells
    for(int x=0; x<cellsPerSide; x++) {
        for(int y=0; y<cellsPerSide; y++) {
            for(int z=0; z<cellsPerSide; z++) {
                int n      = cellNumberFromPosition(x,y,z);
                int xleft  = x-1;    if(xleft==-1)             xleft  = cellsPerSide-1;
                int xright = x+1;    if(xright==cellsPerSide)  xright = 0;
                int yleft  = y-1;    if(yleft==-1)              yleft = cellsPerSide-1;
                int yright = y+1;    if(yright==cellsPerSide)  yright = 0;
                int ztop   = z+1;    if(ztop==cellsPerSide)      ztop = 0;
                neighbouring_cells[n].push_back(cellNumberFromPosition( xright, y     , z   ));
                neighbouring_cells[n].push_back(cellNumberFromPosition( xright, yright, z   ));
                neighbouring_cells[n].push_back(cellNumberFromPosition( x     , yright, z   ));
                neighbouring_cells[n].push_back(cellNumberFromPosition( xleft , yright, z   ));
                neighbouring_cells[n].push_back(cellNumberFromPosition( xleft , yright, ztop));
                neighbouring_cells[n].push_back(cellNumberFromPosition( x     , yright, ztop));
                neighbouring_cells[n].push_back(cellNumberFromPosition( xright, yright, ztop));
                neighbouring_cells[n].push_back(cellNumberFromPosition( xleft , y     , ztop));
                neighbouring_cells[n].push_back(cellNumberFromPosition( x     , y     , ztop));
                neighbouring_cells[n].push_back(cellNumberFromPosition( xright, y     , ztop));
                neighbouring_cells[n].push_back(cellNumberFromPosition( xleft , yleft , ztop));
                neighbouring_cells[n].push_back(cellNumberFromPosition( x     , yleft , ztop));
                neighbouring_cells[n].push_back(cellNumberFromPosition( xright, yleft , ztop));
                if (simulationType == SimulationType::MonteCarlo) {
                    int zdown  = z-1;    if(zdown==-1)              zdown = cellsPerSide-1;
                    neighbouring_cells[n].push_back(cellNumberFromPosition( xleft , y     , z   ));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( xright, yleft , z   ));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( x     , yleft , z   ));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( xleft , yleft , z   ));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( xleft , yright, zdown));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( x     , yright, zdown));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( xright, yright, zdown));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( xleft , y     , zdown));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( x     , y     , zdown));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( xright, y     , zdown));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( xleft , yleft , zdown));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( x     , yleft , zdown));
                    neighbouring_cells[n].push_back(cellNumberFromPosition( xright, yleft , zdown));
                }
            }
        }
    }
}

const std::list<int> cell_lists::getIPCsInNeighbouringCells(int cell) {
    std::list<int> ipcsInNeighbouringCells;
    for(auto const& it: neighbouring_cells[cell])
        std::copy(list_of_neighbours[it].begin(), list_of_neighbours[it].end(), std::back_inserter(ipcsInNeighbouringCells));
    return ipcsInNeighbouringCells;
}
