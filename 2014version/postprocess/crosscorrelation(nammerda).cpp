/********************************************************************************
 * Jan 2014
 * Cross correlations
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
  int nIPCs;
  double L, rho, dt, printEvery, SimLength;
  int nSubSims = 10;
  ifstream INPUT("siml/output.out");
  ofstream CCOUT;
  CCOUT<<scientific<<setprecision(8);
  INPUT>>nIPCs>>rho>>dt;  
  INPUT>>dt>>printEvery>>SimLength;
  nIPCs *= 4*nIPCs*nIPCs;
  L=cbrt(nIPCs/rho);
  INPUT.close();
  int LengthSubSims = int(SimLength/printEvery)/10;
  
  vec *x, *w, *v, *xo, *wo, *vo;
  x = new vec[nIPCs];        w = new vec [nIPCs];           v = new vec [nIPCs];
  xo = new vec[nIPCs];       wo = new vec[nIPCs];           vo = new vec [nIPCs];
  double * ccx = new double[LengthSubSims];
  double * ccw = new double[LengthSubSims];
  double * ccv = new double[LengthSubSims];
  
  for(int i=0; i<LengthSubSims; i++)
  {
    ccx[i] = 0.0; ccw[i] = 0.0; ccv[i] = 0.0;
  }
  
  INPUT.open("siml/trajectory.xyz", ios::in);
  
  for(int ss=0; ss<nSubSims; ss++)
  {
    cout<<ss*10<<" %"<<endl;
    for(int conf=0; conf<LengthSubSims; conf++)
    {
      
      int number; double time; char name;
      INPUT>>number>>time;
      for(int i=0;i<nIPCs;i++)
      {
        INPUT>>name>>x[i].x>>x[i].y>>x[i].z>>v[i].x>>v[i].y>>v[i].z;
        INPUT>>name>>w[i].x>>w[i].y>>w[i].z;
        w[i] -= x[i];      lroundccp(w[i]);        w[i] /= sqrt(w[i]*w[i]);
        INPUT>>time>>time>>time>>name>>time>>time>>time>>time>>time>>time;
      }

      if(conf==0)
      {
        for(int i=0;i<nIPCs;i++)
        {
          xo[i] = x[i];
          wo[i] = w[i];
          vo[i] = v[i];
        }
      } 
      // compute cross correlations
      for(int i=0; i<nIPCs-1; i++)
      {
        for(int j=0; j<nIPCs-1; j++)
        {
          if(i==j) continue;
          ccx[conf] += xo[i]*x[j];
          ccw[conf] += vo[i]*v[j];
          ccv[conf] += wo[i]*w[j];
        }
      }
        
    }
  }
  INPUT.close();
  
  // print autocorrelations
  CCOUT.open("crosscorrelations.anaout");
  CCOUT<<scientific<<setprecision(8);
  for(int c=0; c< LengthSubSims; c++)
  {
    CCOUT<<c<<"\t"<<ccx[c]<<"\t"<<ccw[c]<<"\t"<<ccv[c]<<"\t"<<ccx[c]/ccx[0]<<"\t"<<ccw[c]/ccw[0]<<"\t"<<ccv[c]/ccv[0]<<"\n";
  }
}
