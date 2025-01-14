#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#include "taggedSF.h"
#include "GridInterp.h"
#include "LineInterp.h"
#include "constants.h"
#include "interp_arrays.h"

using namespace std;

double alpha1(double *p1perp, double p1z, double *p2perp, double p2z);
double alpha2(double *p1perp, double p1z, double *p2perp, double p2z);

int main(int argc, char** argv){
  ofstream adep_nomod_file("F2A_adep_nomod.txt");
  if (adep_nomod_file.fail()){
    cerr << "Failed to open adep_nomod_file" << endl;
    exit(-1);
  }

  ofstream xdep_nomod_file("F2A_xdep_nomod.txt");
  if (xdep_nomod_file.fail()){
    cerr << "Failed to open xdep_nomod_file" << endl;
    exit(-1);
  }

  adep_nomod_file << "for p2perp={0,0}, Q2=2" << endl;
  adep_nomod_file << std::fixed << "a2" << std::setw(20) << "xB" << std::setw(20) << "F2A" << std::setw(20) << "err" << endl;
  
  
  xdep_nomod_file << "for p2perp={0,0}, Q2=2" << endl;
  xdep_nomod_file << std::fixed << "a2" << std::setw(20) << "xB" << std::setw(20) << "F2A" << std::setw(20) << "err" << endl;

  //testing pz
  double p1z, p2z, integral_err;
  // double a1 = 0.7;
  //double a2 = 1.3;
  //double p1perp[2] = {0.2,0.2};
  //double p2perp[2] = {.1,.15};

  //double p1perp[2] = {0.,0.};
  double p2perp[2] = {0.,0.};
  double QSq = 2.;
  
  taggedSF * tag_check = new taggedSF();

  for (double xb = 0.3 ; xb <=0.7 ; xb += 0.2){
    for (double alpha2 = 1.0 ; alpha2 < 1.6 ; alpha2 += 0.01){
      double f2a = tag_check->taggedF2(xb, QSq, alpha2, p2perp, integral_err);
      adep_nomod_file << alpha2 << std::setw(20) << xb << std::setw(20) << f2a << std::setw(20) << integral_err <<endl;
    }
  }

  for (double alpha2 = 1.4 ; alpha2 <= 1.6 ; alpha2+=.1){
    for (double xb = 0.20 ; xb < 0.9 ; xb += 0.01){
      double f2a = tag_check->taggedF2(xb, QSq, alpha2, p2perp, integral_err);
      xdep_nomod_file << alpha2 << std::setw(20) << xb << std::setw(20) << f2a << std::setw(20) << integral_err <<endl;
    }
  }

  adep_nomod_file.close();
  xdep_nomod_file.close();
 


  return 0;
}


















double alpha1(double *p1perp, double p1z, double *p2perp, double p2z){

  double p1[3] = {p1perp[0],p1perp[1],p1z};
  double p2[3] = {p2perp[0],p2perp[1],p2z};
  double pcm[3] = {p1[0]+p2[0], p1[1]+p2[1], p1[2]+p2[2]};
  double p1sq = p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] ;
  double p2sq = p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2] ;
  double pcmsq = pcm[0]*pcm[0] + pcm[1]*pcm[1] + pcm[2]*pcm[2] ;

  double alpha1 = 1/mbar * ( sqrt(mN*mN + p1sq) * (mHe4 - sqrt(mH2*mH2 + pcmsq)) / (sqrt(mN*mN + p1sq) + sqrt(mN*mN + p2sq)) + p1z );

  return alpha1;
}


double alpha2(double *p1perp, double p1z, double *p2perp, double p2z){

  double p1[3] = {p1perp[0],p1perp[1],p1z};
  double p2[3] = {p2perp[0],p2perp[1],p2z};
  double pcm[3] = {p1[0]+p2[0], p1[1]+p2[1], p1[2]+p2[2]};
  double p1sq = p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] ;
  double p2sq = p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2] ;
  double pcmsq = pcm[0]*pcm[0] + pcm[1]*pcm[1] + pcm[2]*pcm[2] ;

  
  double alpha2 = 1/mbar * ( sqrt(mN*mN + p2sq) * (mHe4 - sqrt(mH2*mH2 + pcmsq)) / (sqrt(mN*mN + p1sq) + sqrt(mN*mN + p2sq)) + p2z );

  return alpha2;
}
