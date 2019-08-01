/********************************************************************************
 * Jan 2014
 * Computes the single particle mean squared displacement,
 * used in combination with slowfast_particles.sh
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
  int kPoints; double *kGrid;
  time_t start,end;
  ifstream INPUT("siml/output.out");
  ofstream MSDISP("msd_singleparticles.anaout");
  MSDISP<<scientific<<setprecision(8);
  INPUT>>nIPCs>>rho>>dt;  
  INPUT>>dt>>printEvery>>SimLength;
  cout<<"Number of IPCs? ";
  cin>>nIPCs;
  L=cbrt(nIPCs/rho);
  INPUT.close();
  int nSamples = int(SimLength/printEvery);
  
  time (&start);
  vec *x, *w, *v, *xo, *wo, *vo, *xold, *wold, *delta_x; double *delta_w;
  x = new vec[nIPCs];        w = new vec [nIPCs];           v = new vec [nIPCs];
                             wo = new vec[nIPCs];           vo = new vec [nIPCs];
  xold = new vec [nIPCs];    wold = new vec [nIPCs];
  delta_x = new vec[nIPCs];  delta_w = new double [nIPCs];
  
  for(int i=0; i<nIPCs; i++)
  {
    delta_x[i] = vec(0.0, 0.0, 0.0);
    delta_w[i] = 0.0;
  }
  
  INPUT.open("siml/trajectory.xyz", ios::in);
  
  for(int conf=0; conf< nSamples; conf++)
  {
    if( (10*conf)%nSamples == 0 )   cout<<100*conf/nSamples<<"%"<<endl;
    
    int number; double time; char name;
    INPUT>>number>>time;
    for(int i=0;i<nIPCs;i++)
    {
      INPUT>>name>>x[i].x>>x[i].y>>x[i].z>>v[i].x>>v[i].y>>v[i].z;
      INPUT>>name>>w[i].x>>w[i].y>>w[i].z;
      INPUT>>time>>time>>time>>name>>time>>time>>time>>time>>time>>time;
    }

    if(conf>0)
    {
      for(int i=0;i<nIPCs;i++)
      {
        vec dx = x[i] - xold[i];
        lroundccp(dx);
        delta_x[i] += dx;
        xold[i] = x[i];
      }
      // print msdisplacements and total rotation
      MSDISP << conf ;
      for(int i=0;i<nIPCs;i++)
      {
        MSDISP << "\t" << (delta_x[i]*delta_x[i])*L*L ;//<<"\t"<< delta_w[i];
      }
      MSDISP << endl;
      // compute velocity and orientation autocorrelations
    }
    else  //conf == 0
    {
      for(int i=0;i<nIPCs;i++)
      {
        xold[i] = x[i];
      }
    }
  }
  MSDISP.close();
  
  time(&end);
  double dif = difftime (end,start);
  cout<<dif<<" seconds elapsed.\n";
}
