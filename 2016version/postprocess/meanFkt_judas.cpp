/********************************************************************************
 * Jan 2014
 * QuasiFinal version, does the final means over the F (to be fixed!) and Fs files.
 *******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

void giuda();
void ballerino();

int main()
{
  giuda();
  cout<<"Porco chi?\n";
//   ballerino();
  cout<<"Porco giuda!\n";
}
  
void giuda()
{  
  int nFiles(10), columns(8);
  ofstream mediato("ana/Fs_tk");
  mediato<<scientific<<setprecision(8);
  string temp;
  string file = "ana/Fs_";
//   ifstream enne("inputana.dat");
//   float merdate;
//   enne>>merdate>>merdate>>merdate;  enne>>nFiles;
  ifstream *inputs = new ifstream [nFiles];
  
  for (int i = 0; i < nFiles; i++)
  {
    temp = char(i+'0');
    temp = file + temp;
    cout<<temp<<endl;
    inputs[i].open(temp.c_str());
    getline(inputs[i], temp);
  }
  mediato<<temp<<endl;
  
  while ( !inputs[0].eof() )
  {
    int pene(0), culo;
    for (int n=0; n<nFiles; n++)
    {
      inputs[n]>>culo;
      pene+=culo;
    }
    pene/=nFiles;
    mediato<<pene<<"\t";
    
    for (int c=1; c<columns; c++)
    {
    double magheggi(0.), cacca;
      for (int n=0; n<nFiles; n++)
      {
        inputs[n]>>cacca;
        magheggi+=cacca;
      }
      magheggi/=nFiles;
      mediato<<magheggi<<"\t";
    }
    mediato<<endl;
  }
}


void ballerino()
{  
  // to be revised!
  
  int nFiles(10), columns(8);
  ofstream mediato("ana/F_tk");
  mediato<<scientific<<setprecision(8);
  string temp;
  string file = "ana/F_";
//   ifstream enne("inputana.dat");
//   float merdate;
//   enne>>merdate>>merdate>>merdate;  enne>>nFiles;
  ifstream *inputs = new ifstream [nFiles];
  
  for (int i = 0; i < nFiles; i++)
  {
    temp = char(i+'0');
    temp = file + temp;
    cout<<temp<<endl;
    inputs[i].open(temp.c_str());
    getline(inputs[i], temp);
  }
  mediato<<temp<<endl;
  
  while ( !inputs[0].eof() )
  {
    int pene(0), culo;
    for (int n=0; n<nFiles; n++)
    {
      inputs[n]>>culo;
      pene+=culo;
    }
    pene/=nFiles;
    mediato<<pene<<"\t";
    
    for (int c=1; c<columns; c++)
    {
    double magheggi(0.), cacca;
      for (int n=0; n<nFiles; n++)
      {
        inputs[n]>>cacca;
        magheggi+=cacca;
      }
      magheggi/=nFiles;
      mediato<<magheggi<<"\t";
    }
    mediato<<endl;
  }
}