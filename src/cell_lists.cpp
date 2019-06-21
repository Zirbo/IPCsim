#include "cell_lists.hpp"
#include <iostream>

void cell_lists::initialize(double Side, double InteractionRange, int Nparticles) {
    N = Nparticles;
    M = int( Side/InteractionRange);
    M2 = M*M;    M3 = M2*M;
    l = Side/M;
    list_of_neighbours.resize(M3);     // every element is a list with the identifying number of the particles in the list
    neighbouring_cells.resize(M3);     // every element is a list with the number of nearest neighbour cells to that cell
    // fill the list of neighbouring cells
    for(int x=0; x<M; x++) {
        for(int y=0; y<M; y++) {
            for(int z=0; z<M; z++) {
                int n      = cell(x,y,z);
                int xleft  = x-1;    if(xleft==-1)  xleft  = M-1;
                int xright = x+1;    if(xright==M)  xright = 0;
                int yleft  = y-1;    if(yleft==-1)  yleft  = M-1;
                int yright = y+1;    if(yright==M)  yright = 0;
                int ztop   = z+1;    if(ztop==M)    ztop   = 0;
                neighbouring_cells[n].push_back(cell( xright, y     , z   ));
                neighbouring_cells[n].push_back(cell( xright, yright, z   ));
                neighbouring_cells[n].push_back(cell( x     , yright, z   ));
                neighbouring_cells[n].push_back(cell( xleft , yright, z   ));
                neighbouring_cells[n].push_back(cell( xleft , yright, ztop));
                neighbouring_cells[n].push_back(cell( x     , yright, ztop));
                neighbouring_cells[n].push_back(cell( xright, yright, ztop));
                neighbouring_cells[n].push_back(cell( xleft , y     , ztop));
                neighbouring_cells[n].push_back(cell( x     , y     , ztop));
                neighbouring_cells[n].push_back(cell( xright, y     , ztop));
                neighbouring_cells[n].push_back(cell( xleft , yleft , ztop));
                neighbouring_cells[n].push_back(cell( x     , yleft , ztop));
                neighbouring_cells[n].push_back(cell( xright, yleft , ztop));
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
        const int ipcCell = cell(ipc.ipcCenter) ;
        list_of_neighbours[  ipcCell  ].push_back(ipc.number);
      }
    // now the list contains the indices of the particles inside its volume

    // print neighbouring cells for debugging
    std::cout << "compiling lists\n";
    for(int i=0;i<M3;i++)
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
