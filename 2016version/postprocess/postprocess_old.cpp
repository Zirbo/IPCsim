/********************************************************************************
 * Jan 2014
 * Final version, does everything except the final means over the F and Fs files.
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
  int nSubSims = 10;
  int kPoints; double *kGrid;
  time_t start,end;
  ifstream INPUT("siml/output.out");
  ofstream OUTPUT("ana/minchiate.anaout");
  ofstream MSDISP("ana/meansquareddisplacement.anaout");
  MSDISP<<scientific<<setprecision(8);
  INPUT>>nIPCs>>rho>>dt;  
  INPUT>>dt>>printEvery>>SimLength;
  INPUT.close();
  int nSamples = int(SimLength/printEvery);
  
  INPUT.open("processinput.in");
  if( !INPUT.is_open() ) { cerr<<"Giai mi processinput.in, cunn'e mongiasa!!!!\n"; exit(1); }
  
  INPUT>>nIPCs>>L>>nSubSims>>kPoints;
  int LengthSubSims = nSamples/nSubSims;
  kGrid = new double [kPoints];
  for (int i=0; i< kPoints; i++) INPUT>>kGrid[i];
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
  S_k senno(20,L,nIPCs);
  INPUT.open("siml/trajectory.xyz", ios::in);
  
  for(int ss=0; ss<nSubSims; ss++)
  {
    OUTPUT<<"subsim number "<<ss+1<<endl;
//     F_k nenno(20,L,nIPCs,LengthSubSims);
    Fs_k trenno(kPoints,kGrid,L,nIPCs,LengthSubSims);
    trenno.sketch(OUTPUT);
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
//       nenno.compute(x);
      trenno.compute(x);
      
//       cout<<"\t\t\t"<<vcorr[conf]<<"  "<<wcorr[conf]<<endl;
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
        double delta_x2(0.0), deltarot(0.0);
        for(int i=0; i<nIPCs; i++)
        {
          delta_x2 += delta_x[i]*delta_x[i];
          deltarot += delta_w[i];
        }
        delta_x2 /= nIPCs;
        deltarot /= nIPCs;
        MSDISP << conf+ss*LengthSubSims <<"\t"<< delta_x2*L*L <<"\t"<< deltarot <<endl;
        // compute velocity and orientation autocorrelations
        for(int i=0; i<nIPCs; i++)
        {
          vcorr[conf] += v[i]*vo[i];
          wcorr[conf] += w[i]*wo[i];
//           cout<<v[i]*vo[i]<<"  "<<w[i]*wo[i]<<"  "<<vcorr[conf]<<"  "<<wcorr[conf]<<endl;
        }
      }
      else  //conf == 0
      {
//           cout<< "zeroing number "<<ss<<endl;
        for(int i=0;i<nIPCs;i++)
        {
          xold[i] = x[i];
          wold[i] = w[i];
          vo[i] = v[i];
          wo[i] = w[i];
          vcorr[conf] += v[i]*vo[i];
          wcorr[conf] += w[i]*wo[i];
//           cout<<v[i]*vo[i]<<"  "<<w[i]*wo[i]<<"  "<<vcorr[conf]<<"  "<<wcorr[conf]<<endl;
        }
      }
//       cout<<endl<<"fineconf "<<ss*LengthSubSims+conf<<endl<<endl;
    }
    trenno.print(OUTPUT);
    string s98, s99("mv Fs_tk ana/Fs_0");
    stringstream cacca;   cacca<<ss;   cacca>>s98;
    s99.replace(16,1,s98);
    if( system(s99.c_str()) != 0 )
    {
      cerr<<"Cannot move Fs_"<<ss<<"...\n";
    }
//     nenno.print();
//     string s96, s97("mv F_tk ana/F_0");
//     cacca.clear();   cacca<<ss;   cacca>>s96;
//     s97.replace(14,1,s96);
//     if( system(s97.c_str()) != 0 )
//     {
//       cerr<<"Cannot move F_"<<ss<<"...\n";
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
  OUTPUT<<dif<<" seconds elapsed.\n";
}
