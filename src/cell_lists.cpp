#include "cell_lists.hpp"
#include <iostream>

void cell_lists::initialize(double simulationBoxSide, double interactionRange, int howManyParticles) {
    numberOfParticles = howManyParticles;
    cellsPerSide = int( simulationBoxSide/interactionRange);
    cellsPerSideSquared = cellsPerSide*cellsPerSide;
    totalCells = cellsPerSideSquared*cellsPerSide;
    cellSide = simulationBoxSide/cellsPerSide;
    list_of_neighbours.resize(totalCells);     // every element is a list with the identifying number of the particles in the list
    neighbouring_cells.resize(totalCells);     // every element is a list with the number of nearest neighbour cells to that cell
    // fill the list of neighbouring cells
    for(int x=0; x<cellsPerSide; x++) {
        for(int y=0; y<cellsPerSide; y++) {
            for(int z=0; z<cellsPerSide; z++) {
                int n      = cellNumberFromPosition(x,y,z);
                int xleft  = x-1;    if(xleft==-1)  xleft  = cellsPerSide-1;
                int xright = x+1;    if(xright==cellsPerSide)  xright = 0;
                int yleft  = y-1;    if(yleft==-1)  yleft  = cellsPerSide-1;
                int yright = y+1;    if(yright==cellsPerSide)  yright = 0;
                int ztop   = z+1;    if(ztop==cellsPerSide)    ztop   = 0;
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
            }
        }
    }
}
void cell_lists::compilelists(std::vector<IPC> const& ipcs) {
    // empty all the lists
    for(auto & m: list_of_neighbours)
        m.clear();
    // put every particle in the right list
    for(IPC const& ipc: ipcs) {
        const int ipcCell = cellNumberFromPosition(ipc.ipcCenter) ;
        list_of_neighbours[  ipcCell  ].push_back(ipc.number);
      }
    // now the list contains the indices of the particles inside its volume

    // print neighbouring cells for debugging
    std::cout << "compiling lists\n";
    for(int i=0;i<totalCells;i++)
    {
        std::cout << "list nr " << i << ":   ";
        for(auto it = list_of_neighbours[i].cbegin(); it!=list_of_neighbours[i].end(); it++)
            std::cout<<*it<<"    ";
        std::cout<<std::endl;
    }
}
const std::list<int> cell_lists::getIPCsInNeighbouringCells(int cell) {
    std::list<int> ipcsInNeighbouringCells;
    for(auto const& it: neighbouring_cells[cell])
        std::copy(list_of_neighbours[it].begin(), list_of_neighbours[it].end(), std::back_inserter(ipcsInNeighbouringCells));
    return ipcsInNeighbouringCells;
}
