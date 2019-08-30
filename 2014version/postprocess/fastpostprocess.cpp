/********************************************************************************
 * 14 July 2014 (Allons enfants de la patrie, le jours de gloire est arrive'...)
 * Final version, does everything including the PORCOGGIUDA ordeal!
 *******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "zilvo.hpp"
#include <cstdlib>

using namespace space;
using namespace std;

void porcoggiudaporco(int nFiles);
void mirrorize(int tmax);

int main()
{
  if( system("rm -r ana") != 0 )    cerr<<"But that's not a problem, don't worry!\n";
  if( system("mkdir ana") != 0 )
  {
    cerr<<"Cannot create directory ana! Quitting.\n";
    exit(1);
  }
  
  int nIPCs;
  double L, rho, dt, printEvery, SimLength;
  int nSubSims = 10;
  double kmax;
  time_t start,end;
  ifstream INPUT("siml/output.out");  //read parameters from the sim
  ofstream OUTPUT("ana/minchiate.anaout");  //print state to check how far is from finishing
  ofstream MSDISP("ana/meansquareddisplacement.anaout"); //print kitties and butterflies
  MSDISP<<scientific<<setprecision(8);
  INPUT>>nIPCs>>rho>>dt;  
  INPUT>>dt>>printEvery>>SimLength;
  INPUT.close();
  int nSamples = int(SimLength/printEvery);
  
  INPUT.open("processinput.in");  //this file must contain the parameters below!
  if( !INPUT.is_open() ) { cerr<<"Giai mi processinput.in, cunn'e mongiasa!!!!\n"; exit(1); }
  
  INPUT>>nIPCs>>L>>nSubSims>>kmax; //nIPCs has to be the correct one. This way
  //the program works with all my different versions of sim. nSubSims is the number
  // of subsims in which the trajectory has to be divided. This shortens the RAM
  // requirements AND improves the statistics for large times. suggested by Jan Kurzidim!
  int LengthSubSims = nSamples/nSubSims;
  INPUT.close();
  
  time (&start);
  vec *x, *w, *v, *xo, *wo, *vo, *xold, *wold, *delta_x; double *delta_w;
  x = new vec[nIPCs];        w = new vec [nIPCs];           v = new vec [nIPCs];
                             wo = new vec[nIPCs];           vo = new vec [nIPCs];
  xold = new vec [nIPCs];    wold = new vec [nIPCs];
  delta_x = new vec[nIPCs];  delta_w = new double [nIPCs];
  double * vcorr = new double[LengthSubSims];  double * wcorr = new double[LengthSubSims];
  
  for(int i=0; i<nIPCs; i++)
  {
    delta_x[i] = vec(0.0, 0.0, 0.0);
    delta_w[i] = 0.0;
  }
  
  for(int i=0; i<LengthSubSims; i++)
  {
    vcorr[i] = 0.0; wcorr[i] = 0.0;
  }
  
  G   gigia(50,L,nIPCs,OUTPUT,false);
  NEW_Sk senno(20,L,nIPCs);
  INPUT.open("siml/trajectory.xyz", ios::in);
  
  for(int ss=0; ss<nSubSims; ss++)
  {
    OUTPUT<<"subsim number "<<ss+1<<endl;
//     NEW_Fkt trenno(kmax,L,nIPCs,LengthSubSims);
//     trenno.sketch(OUTPUT);
    for(int conf=0; conf< LengthSubSims; conf++)
    {
      if( (10*conf)%LengthSubSims == 0 )   OUTPUT<<100*conf/LengthSubSims<<"%"<<endl;
      
      int number; double time; char name;
      INPUT>>number>>time;
      for(int i=0;i<nIPCs;i++)
      {
        INPUT>>name>>x[i].x>>x[i].y>>x[i].z>>v[i].x>>v[i].y>>v[i].z;
        INPUT>>name>>w[i].x>>w[i].y>>w[i].z;
        w[i] -= x[i];      lroundccp(w[i]);        w[i] /= sqrt(w[i]*w[i]);
        INPUT>>time>>time>>time>>name>>time>>time>>time>>time>>time>>time;
      }

      // compute cfs
      gigia.compute(x,w);
      senno.compute(x);
//       trenno.compute(x);
      
      if(conf>0)
      {
        for(int i=0;i<nIPCs;i++)
        {
          vec dx = x[i] - xold[i];
          lroundccp(dx);
          delta_x[i] += dx;
          double deltatheta = wold[i]*w[i];
          if(deltatheta>1.01)
          {
            cerr<<"SOMETHING IS REALLY WRONG IN IPC "<<i<<" IN CONFIGURATION "<<conf<<"!";
            exit(1);
          } // the following is a correction for small errors in the orientation
          else if(deltatheta>1.) deltatheta = 1.;
          delta_w[i] += abs(acos(deltatheta));
          xold[i] = x[i];
          wold[i] = w[i];
        }
        // compute and print msdisplacements and total rotation
        double delta_x2(0.0), deltarot(0.0), deltaz(0.0), deltapar(0.0);
        for(int i=0; i<nIPCs; i++)
        {
          delta_x2 += delta_x[i]*delta_x[i];
          deltarot += delta_w[i];
          deltaz   += delta_x[i].z*delta_x[i].z;
          deltapar += delta_x[i].x*delta_x[i].x + delta_x[i].y*delta_x[i].y;
        }
        delta_x2 /= nIPCs;
        deltarot /= nIPCs;
        deltaz   /= nIPCs;
        deltapar /= nIPCs;
        MSDISP << conf+ss*LengthSubSims <<"\t"<< delta_x2*L*L <<"\t"<< deltarot;
        MSDISP << "\t"<< deltaz*L*L <<"\t"<< deltapar*L*L <<endl;
        // compute velocity and orientation autocorrelations
        for(int i=0; i<nIPCs; i++)
        {
          vcorr[conf] += v[i]*vo[i];
          wcorr[conf] += w[i]*wo[i];
        }
      }
      else  //conf == 0
      {     // sets the current conf as the C(0) conf
        for(int i=0;i<nIPCs;i++)
        {
          xold[i] = x[i];
          wold[i] = w[i];
          vo[i] = v[i];
          wo[i] = w[i];
          vcorr[conf] += v[i]*vo[i];
          wcorr[conf] += w[i]*wo[i];
        }
      }
    }
//     trenno.print(OUTPUT);
//     string s98, s99("mv Fs_tk ana/Fs_0");
//     stringstream cacca;   cacca<<ss;   cacca>>s98;
//     s99.replace(16,1,s98);
//     if( system(s99.c_str()) != 0 )
//     {
//       cerr<<"Cannot move Fs_"<<ss<<"...\n";
//     }
  }
  INPUT.close();
  
  gigia.print(OUTPUT);
  if( system("mv g ana/") != 0 )
    cerr<<"Cannot move g...\n";
  senno.print();
  if( system("mv S_k ana/") != 0 )
    cerr<<"Cannot move S(k)...\n";
  
  // print autocorrelations
  MSDISP.close();  MSDISP.open("ana/autocorrelations.anaout");
  MSDISP<<scientific<<setprecision(8);
  for(int c=0; c< LengthSubSims; c++)
  {
    MSDISP<<c<<"\t"<<vcorr[c]/vcorr[0]<<"\t"<<wcorr[c]/wcorr[0]<<"\n";
  }
  
  time(&end);
  double dif = difftime (end,start);
  cout<<"Porco chi?\n";
//   porcoggiudaporco(nSubSims);  //averages the 
  cout<<"Porco giuda!\n";
//   mirrorize(LengthSubSims);
  OUTPUT<<dif<<" seconds elapsed;.\n";
}

void porcoggiudaporco(int nFiles)
{
  int columns(8);
  ofstream mediato("ana/Fs_tk");
  mediato<<scientific<<setprecision(8);
  string temp;
  string file = "ana/Fs_";
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

void mirrorize(int tmax)
{
  ifstream cacca("ana/Fs_tk");
  ofstream specchiato("ana/Fskt");
  int kgrid = (std::count(std::istreambuf_iterator<char>(cacca), std::istreambuf_iterator<char>(), '\n')-1)/tmax;
  cacca.seekg(0, ios::beg);
  cout<<tmax<<" = tmax, kgrid =  "<<kgrid<<endl;
  string *outfile = new string [tmax*kgrid];
  string mistalliga;
  getline(cacca,mistalliga);
  cout<<mistalliga;
  for(int k=0;k<kgrid;k++)
  {
    for(int t=0;t<tmax;t++)
    {
      // ich lese Reihe t+tmax*k, die k+kgrid*t werden muss
      getline(cacca,outfile[k+kgrid*t]);
    }
  }
  specchiato<<mistalliga<<std::endl;
  for(int i=0; i<kgrid*tmax; i++)
  {
    specchiato<<outfile[i]<<std::endl;
  }
}
