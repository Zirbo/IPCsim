#include "cell_lists.hpp"

void cell_lists::initialize(double Side, double InteractionRange, int Nparticles) {
    N = Nparticles;
    M = int( Side/InteractionRange);
    M2 = M*M;    M3 = M2*M;
    l = Side/M;
    lista  = new std::list<int> [M3];     // every element is a list with the number of particles in the list
    vicini = new std::list<int> [M3];     // every element is a list with the number of nearest neighbour cells to that cell
    // fill the vicini list
    for(int x=0; x<M; x++) {
        for(int y=0; y<M; y++) {
            for(int z=0; z<M; z++) {
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
}
inline int cell_lists::cell(Particle const& x)
{  // gives you the list where this coordinate belong
    return int(x.x[0]/l)+M*int(x.x[1]/l)+M2*int(x.x[2]/l);
}
inline int cell_lists::cell(int x, int y, int z)
{  // hope it really inlines because it's really stupid to have it
    return x + M*y + M2*z;
}
void cell_lists::compilelists(std::vector<IPC> const& ipcs) {  // empty all the lists
    for(int m=0; m<M3; m++)
        lista[m].clear();
    // put every particle in the right list
    for(int i=0; i<N; i++)
        lista[   cell(ipcs[i].center)   ].push_back(i);
    // now the list contains the indices of the particles inside its volume
}
void cell_lists::neighbour_cells(int Cell, std::list<int> &local, std::list<int> &neigh)
{
    local = lista[Cell];
    std::list<int> stampecullu;
    for(std::list<int>::iterator it = vicini[Cell].begin(); it!=vicini[Cell].end(); it++) {
      stampecullu = lista[*it];
      neigh.merge(stampecullu);
    }
}
