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
  /*  ofstream myfile;
  myfile.open("alphapz_table.txt");

  myfile << "for p1perp=p2perp={0,0}" << endl;
  myfile << std::fixed << "a_1" << std::setw(20) << "a_2" << std::setw(20) << "p1z" << std::setw(20) << "p2z" << endl;*/
  
  //testing pz
  double p1z, p2z;
  double a1 = 0.7;
  double a2 = 1.3;
  //double p1perp[2] = {0.2,0.2};
  //double p2perp[2] = {.1,.15};

  double p1perp[2] = {0.,0.};
  double p2perp[2] = {0.,0.};

  taggedSF * tag_check = new taggedSF();
  /*
  bool end = false;//pz inversion inaccurate for all a1 > a1_first_fail
  
  for (double alph2 = 1.0 ; alph2 < 2. ; alph2 += .01){
    end = false;
    for (double alph1 = 0. ; alph1 < 3. ; alph1 += .01){
      tag_check->Findpz(alph1,p1perp,alph2,p2perp,p1z,p2z);
      if ( (p1z < 10.) and (p1z > -10.) and (p2z < 10.) and (p2z > -10.) and (!end) ){
	myfile <<std::fixed << std::setprecision(2) << alph1 << setw(20) << alph2 << setw(20) << std::setprecision(6) << p1z << setw(20) << p2z << endl;
      }
      else{
	end = true;
	myfile <<std::fixed << std::setprecision(2) << alph1 << setw(20) << alph2 << setw(20) << std::setprecision(6) << 10000 << setw(20) << 10000 << endl;
      }
    }
  }
  myfile.close();*/
  
  /*
  for (int i = 0 ; i <= 400 ; i++){
    double alph1 = ((double)i)/100.;
    tag_check->Findpz(alph1,p1perp,a2,p2perp,p1z,p2z);
    double a1check = alpha1(p1perp,p1z,p2perp,p2z);
    double a2check = alpha2(p1perp,p1z,p2perp,p2z);
    cout << "\nalpha1, alpha2: " << alph1 << ", " << a2 << endl;
    cout << "calc'd p1z, p2z : " << p1z << ", " << p2z << endl;
    cout << "recalc'd alpha1, alpha2: " << a1check <<  ", " << a2check << endl;
    }
  */


  //double xb = 0.4;
  //double Q2 = 2.;


  
  double F2check = tag_check->taggedF2(0.4, 2., a2 , p2perp);

  cout << "\n tagged F2 " << F2check << endl;
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
