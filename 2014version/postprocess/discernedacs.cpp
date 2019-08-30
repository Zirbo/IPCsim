/********************************************************************************
 * Jan 2014
 * Using the outputfile from slowfast_particles.sh can compute the divided autocorrelations
 * for the chocolate layers and wafer planes.
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
  INPUT>>nIPCs>>rho>>dt;  
  INPUT>>dt>>printEvery>>SimLength;
  nIPCs *= 4*nIPCs*nIPCs;
  L=cbrt(nIPCs/rho);
  INPUT.close();
  int nSamples = int(SimLength/printEvery);
  
  INPUT.open("processinput.in");
  if( !INPUT.is_open() ) { cerr<<"Giai mi processinput.in, cunn'e mongiasa!!!!\n"; exit(1); }
  INPUT>>nSubSims;
  int LengthSubSims = nSamples/nSubSims;
  INPUT.close();
  
  int Nfast;
  INPUT.open("fastpart");
  if( !INPUT.is_open() ) { cerr<<"Giai mi fastpart, cunn'e mongiasa!!!!\n"; exit(1); }
  INPUT>>Nfast;
  int * fastandfurious = new int [Nfast];
  for(int i=0; i<Nfast; i++)
    INPUT>>fastandfurious[i];
  INPUT.close();
  
  
  
  time (&start);
  vec *x, *w, *v, *xo, *wo, *vo;
  x = new vec[nIPCs];        w = new vec [nIPCs];           v = new vec [nIPCs];
                             wo = new vec[nIPCs];           vo = new vec [nIPCs];
  double * vcorrP = new double[LengthSubSims];  double * wcorrP = new double[LengthSubSims];
  double * vcorrM = new double[LengthSubSims];  double * wcorrM = new double[LengthSubSims];
  
  
  for(int i=0; i<LengthSubSims; i++)
  {
    vcorrP[i] = 0.0; wcorrP[i] = 0.0;
    vcorrM[i] = 0.0; wcorrM[i] = 0.0;
  }
  INPUT.open("siml/trajectory.xyz", ios::in);
  
  for(int ss=0; ss<nSubSims; ss++)
  {
    cout<<"subsim number "<<ss+1<<endl;
    for(int conf=0; conf< LengthSubSims; conf++)
    {
      if( (10*conf)%LengthSubSims == 0 )   cout<<100*conf/LengthSubSims<<"%"<<endl;
      
      int number; double time; char name;
      INPUT>>number>>time;
      for(int i=0;i<nIPCs;i++)
      {
        INPUT>>name>>x[i].x>>x[i].y>>x[i].z>>v[i].x>>v[i].y>>v[i].z;
        INPUT>>name>>w[i].x>>w[i].y>>w[i].z;
        w[i] -= x[i];      lroundccp(w[i]);        w[i] /= sqrt(w[i]*w[i]);
        INPUT>>time>>time>>time>>name>>time>>time>>time>>time>>time>>time;
      }
      if(conf>0)
      {
        int j=0;
        // compute velocity and orientation autocorrelations
        for(int i=0; i<nIPCs; i++)
        {
          if (i==fastandfurious[j])
          {
            j++;
            vcorrM[conf] += v[i]*vo[i];
            wcorrM[conf] += w[i]*wo[i];
          }
          else
          {
            vcorrP[conf] += v[i]*vo[i];
            wcorrP[conf] += w[i]*wo[i];
          }
        }
      }
      else  //conf == 0
      {
        int j=0;
        for(int i=0;i<nIPCs;i++)
        {
          vo[i] = v[i];
          wo[i] = w[i];
          if (i==fastandfurious[j])
          {
            j++;
            vcorrM[conf] += v[i]*vo[i];
            wcorrM[conf] += w[i]*wo[i];
          }
          else
          {
            vcorrP[conf] += v[i]*vo[i];
            wcorrP[conf] += w[i]*wo[i];
          }
        }
      }
    }
  }
  INPUT.close();
  
  // print autocorrelations
  ofstream LADRO("discerned_acs.anaout");
  LADRO<<scientific<<setprecision(8);
  for(int c=0; c< LengthSubSims; c++)
  {
    LADRO<<c<<"\t"<<vcorrP[c]/vcorrP[0]<<"\t"<<wcorrP[c]/wcorrP[0];
    LADRO<<   "\t"<<vcorrM[c]/vcorrM[0]<<"\t"<<wcorrM[c]/wcorrM[0]<<"\n";
  }
  
  time(&end);
  double dif = difftime (end,start);
  cout<<dif<<" seconds elapsed.\n";
}
