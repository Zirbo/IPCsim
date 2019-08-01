/********************************************************************************
 * Dic 2013
 * ho visto un nero con la sciarpa di lana...
 *******************************************************************************/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "zilvo.hpp"
#include <cstdlib>

using namespace space;
using namespace std;


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
  time_t start,end;
  ifstream INPUT("siml/output.out");
  ofstream OUTPUT("ana/minchiate.anaout");
  ofstream MSDISP("ana/meansquareddisplacement.anaout");
  INPUT>>nIPCs>>rho>>dt;  
  INPUT>>dt>>printEvery>>SimLength;
  nIPCs *= 4*nIPCs*nIPCs;
  L=cbrt(nIPCs/rho);
  INPUT.close();
  int nSamples = int(SimLength/printEvery);

  time (&start);
  vec *x, *w, *v, *xo, *wo, *vo, *xold, *wold, *delta_x; double *delta_w;
  x = new vec[nIPCs];        w = new vec [nIPCs];           v = new vec [nIPCs];
//   xo = new vec [nIPCs];      wo = new vec[nIPCs];           vo = new vec [nIPCs];
  xold = new vec [nIPCs];    wold = new vec [nIPCs];
  delta_x = new vec[nIPCs];  delta_w = new double [nIPCs];
  
  for(int i=0; i++; i<nIPCs)
  {
    delta_x[i] = vec(0.0, 0.0, 0.0);
    delta_w[i] = 0.0;
  }
  
  
  G   gigia(20,L,nIPCs,OUTPUT,false);
  S_k senno(20,L,nIPCs);
  F_k nenno(20,L,nIPCs,nSamples);
  
  INPUT.open("siml/trajectory.xyz", ios::in);
  
  for(int conf=0; conf<nSamples; conf++)
  {
    if( (20*conf)%nSamples == 0 )   OUTPUT<<100*conf/nSamples<<"%"<<endl;
    
    int number; double time; char name;
    INPUT>>number>>time;
    for(int i=0;i<nIPCs;i++)
    {
      INPUT>>name>>x[i].x>>x[i].y>>x[i].z>>v[i].x>>v[i].y>>v[i].z;
      INPUT>>name>>w[i].x>>w[i].y>>w[i].z;
      w[i] -= x[i];      lroundccp(w[i]);        w[i] /= w[i]*w[i];
      INPUT>>time>>time>>time>>name>>time>>time>>time>>time>>time>>time;
    }

    // compute cfs
    gigia.compute(x,w);
    senno.compute(x);
    nenno.compute(x);
    
    if(conf>0)
    {
      for(int i=0;i<nIPCs;i++)
      {
        vec dx = x[i] - xold[i];
        lroundccp(dx);
        
        delta_x[i] += dx;
        delta_w[i] += abs(acos(wold[i]*w[i]));
        xold[i] = x[i];
        wold[i] = w[i];
      }
      double delta_x2(0.0), deltarot(0.0);
      for(int i=0; i<nIPCs; i++)
      {
        delta_x2 += delta_x[i]*delta_x[i];
        deltarot += delta_w[i]*delta_w[i];
      }
      delta_x2 /= nIPCs;
      deltarot /= nIPCs;
      MSDISP << conf <<"\t"<< delta_x2*L*L <<"\t"<< deltarot <<endl;
    }
    else
    {
      for(int i=0;i<nIPCs;i++)
      {
        xold[i] = x[i];
        wold[i] = w[i];
      }
    }
  }
  INPUT.close();
  
  gigia.print(OUTPUT);
  if( system("mv g ana/") != 0 )
    cerr<<"Cannot move g...\n";
  senno.print();
  if( system("mv S_k ana/") != 0 )
    cerr<<"Cannot move S(k)...\n";
  nenno.print();
  if( system("mv F_tk ana/") != 0 )
    cerr<<"Cannot move F(t,k)...\n";
  
  time(&end);
  double dif = difftime (end,start);
  OUTPUT<<dif<<" seconds elapsed.\n";
}
