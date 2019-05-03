#include "cell_lists.hpp"

void cell_lists::initialize(double Side, double InteractionRange, int Nparticles, space::vec x[])
{
  N = Nparticles;
  M = int( Side/InteractionRange );
  M2 = M*M;    M3 = M2*M;
  l = Side/M;
  lista  = new std::list<int> [M3];     // every element is a list with the number of particles in the list
  vicini = new std::list<int> [M3];     // every element is a list with the number of nearest neighbour cells to that cell
  // fill the vicini list
  for(int x=0; x<M; x++)
  {
    for(int y=0; y<M; y++)
    {
      for(int z=0; z<M; z++)
      {
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
  /* print neighbouring cells for debugging
  for(int i=0;i<M3;i++)
  {
    for(std::list<int>::iterator it = vicini[i].begin(); it!=vicini[i].end(); it++)
      std::cout<<*it<<"    ";
    std::cout<<std::endl;
  }*/
}
cell_lists::cell_lists(double Side, double InteractionRange, int Nparticles, space::vec x[])
{  // This is here only for backwards compatibility
  this->initialize(Side,InteractionRange,Nparticles,x);
}
inline int cell_lists::cell(space::vec x)
{  // gives you the list where this coordinate belong
  return int(x.x/l)+M*int(x.y/l)+M2*int(x.z/l);
}
inline int cell_lists::cell(int x, int y, int z)
{  // hope it really inlines because it's really stupid to have it
  return x + M*y + M2*z;
}
void cell_lists::compilelists(space::vec x[])
{  // empty all the lists
  for(int m=0; m<M3; m++)
    lista[m].clear();
  // put every particle in the right list
  for(int i=0; i<N; i++)
    lista[   cell(x[i])   ].push_back(i);
  // now the list contains the indices of the particles inside its volume

  // print out cell content for debugging
//   for(int m=0; m<M3; m++)
//   {
//     std::cout<<m<<":  ";
//     for(std::list<int>::iterator it = lista[m].begin(); it!=lista[m].end(); it++)
//       std::cout<<*it<<"  ";
//     std::cout<<std::endl;
//   }

}
void cell_lists::neighbour_cells(int Cell, std::list<int> &local, std::list<int> &neigh)
{
  local = lista[Cell];
  std::list<int> stampecullu;
  for(std::list<int>::iterator it = vicini[Cell].begin(); it!=vicini[Cell].end(); it++)
  {
//     std::cout<<*it<<std::endl;
    stampecullu = lista[*it];
    neigh.merge(stampecullu);
  }
}
