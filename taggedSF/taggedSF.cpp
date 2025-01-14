#include "taggedSF.h"

#include "TRandom3.h"
#include "TMath.h"

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <random>
#include <iomanip>

#include "GridInterp.h"
#include "LineInterp.h"
#include "interp_arrays.h"
#include "constants.h"
#include "Math/ProbFuncMathCore.h"

taggedSF::taggedSF()
{
  myRand = new TRandom3(0);

  int qcount = sizeof( q_list ) / sizeof( q_list[0] );
  int Qcount = sizeof( Q_list ) / sizeof( Q_list[0] );
  int pcount = sizeof( n_MF_array )/sizeof( n_MF_array[0] );


  //setting up n_np 2d interp
  n_pn_interp = new GridInterp();
  n_pn_interp->setXVals(qcount, q_list);
  n_pn_interp->setYVals(Qcount, Q_list);

  for (int iy = 0 ; iy < Qcount ; iy++){
    for (int ix = 0 ; ix < qcount ; ix++){
      int index = (iy * qcount) + ix;
      n_pn_interp->setPoint( ix , iy, n_pn_array[index]);
    }
  }

  //setting up n_nn 2d interp
  n_nn_interp = new GridInterp();
  n_nn_interp->setXVals(qcount, q_list);
  n_nn_interp->setYVals(Qcount, Q_list);

  for (int iy = 0 ; iy < Qcount ; iy++){
    for (int ix = 0 ; ix < qcount ; ix++){
      int index = (iy * qcount) + ix;
      n_nn_interp->setPoint( ix , iy, n_nn_array[index]);
    }
  }

  //setting up n_MF 1d interp
  n_MF_interp = new LineInterp();                                          
  n_MF_interp->setXVals(pcount, p_array);

  for (int i = 0 ; i < pcount ; i++){                     
    n_MF_interp->setPoint( i , n_MF_array[i]); 
  }

  //setting up n_tot 1d interp  
  n_tot_interp = new LineInterp();                          
  n_tot_interp->setXVals(pcount, p_array);
  
  for (int i = 0 ; i < pcount ; i++){                  
    n_tot_interp->setPoint( i , n_tot_array[i]);       
  }


  
}
  
taggedSF::~taggedSF()
{
  if (myRand){
    delete myRand;
  }

  delete n_pn_interp;
  delete n_nn_interp;
  delete n_MF_interp;
  delete n_tot_interp;
  
}
 


void taggedSF::Findpz(double alpha1, double *p1perp, double alpha2, double *p2perp, double &p1z, double &p2z){
  //predefine p1z and p2z to pass in to be filled

  double p1perp2 = p1perp[0]*p1perp[0] + p1perp[1]*p1perp[1];
  double p2perp2 = p2perp[0]*p2perp[0] + p2perp[1]*p2perp[1];

  double p1z_min = -10.;//in GeV
  double p1z_max = 10.;

  double tol = 1.e-10; //tolerance, adjust later
  //can go to 10^-14 if needed
  double diff = 10.;

  long int index = 0;
  //for (int i = 0 ; i < 10000 ; i++){
  while ( fabs(diff) > tol){

    if (index == 100){
      std::cerr << "\nInversion Failure" << std::endl;
      std::cerr << "\n\nalpha1 is " << alpha1 << std::endl;
      std::cerr << "p1perp is {" << p1perp[0] << " , " << p1perp[1] << "}" << std::endl;      
      std::cerr << "alpha2 is " << alpha2 << std::endl;
      std::cerr << "p2perp is {" << p2perp[0] << " , " << p2perp[1] << "}" << std::endl;
      std::cerr << "p1z is " << p1z << " and p2z is " << p2z << std::endl;
      exit(-2);}//if stuck in loop
    


    /*    if (index % 1000000000 == 0){
      std::cout << "on iteration " << index/1000000000 << " bil" << std::endl;
      }*/
    p1z = 0.5 * (p1z_max+ p1z_min);

    double pzcm = ( mH2*mH2 + ( (p1perp[0]+p2perp[0]) * (p1perp[0]+p2perp[0]) ) + ( (p1perp[1]+p2perp[1]) * (p1perp[1]+p2perp[1]) ) - ( ( mbar * mbar ) * ( 4. - alpha1 - alpha2 ) * ( 4. - alpha1 - alpha2 ) ) ) / ( 2. * mbar * ( 4. - alpha1 - alpha2 ) );

    p2z = pzcm - p1z;
    
    diff = sqrt(mN*mN+p2z*p2z+p2perp2)*(alpha1*mbar-p1z) - sqrt(mN*mN+p1z*p1z+p1perp2)*(alpha2*mbar-p2z);
    /*if (fabs(diff)<tol){
      break;
      }*/
    //std::cout << "diff is " << diff << std::endl;


    /* if (index > 44){
      std::cout << "\nIndex " << index << std::endl;
      std::cout << "p1z_min is " << p1z_min << " and p1z_max is " << p1z_max << std::endl;
      std::cout << "p1z is " << p1z << std::endl;
      std::cout << "diff is " << diff << std::endl;
      }*/
    
    if (diff > 0.0){
      p1z_min=p1z;
    }
    else{
      p1z_max=p1z;
    }

    index += 1;
  }
}

double taggedSF::n_pn(double q, double Q){
  //  std::cout << "interp n_pn" << std::endl;
  return n_pn_interp->getVal( q , Q );
  
}

double taggedSF::n_nn(double q, double Q){
  //  std::cout << "interp n_nn" << std::endl;
  return n_nn_interp->getVal( q , Q );
  
}

double taggedSF::n_MF(double p){
  //  std::cout << "interp n_MF" << std::endl;
  return  n_MF_interp->getVal( p );

}

double taggedSF::n_tot(double p){     
  //  std::cout << "interp n_tot" << std::endl;                                                
  return  n_tot_interp->getVal( p );

}


double taggedSF::da1dp1z(double *p1perp, double p1z, double *p2perp, double p2z){
  double p1sq = p1perp[0]*p1perp[0] + p1perp[1]*p1perp[1] + p1z*p1z;
  double p2sq = p2perp[0]*p2perp[0] + p2perp[1]*p2perp[1] + p2z*p2z;
  double pcmsq = (p1perp[0]+p2perp[0])*(p1perp[0]+p2perp[0]) + (p1perp[1]+p2perp[1])*(p1perp[1]+p2perp[1]) + (p1z+p2z)*(p1z+p2z);
  double da1dp1z = A/mHe4*(1+p1z/(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq)) * ( (mHe4-sqrt(mH2*mH2+pcmsq)) / sqrt(mN*mN+p1sq)-sqrt(mN*mN+p1sq)/sqrt(mH2*mH2+pcmsq)-(mHe4-sqrt(mH2*mH2+pcmsq))/(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq)) ) - p2z/(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq))*sqrt(mN*mN+p1sq)/sqrt(mH2*mH2+pcmsq));

  return da1dp1z;
}

double taggedSF::da1dp2z(double *p1perp, double p1z, double *p2perp, double p2z){
 double p1sq = p1perp[0]*p1perp[0] + p1perp[1]*p1perp[1] + p1z*p1z;
 double p2sq = p2perp[0]*p2perp[0] + p2perp[1]*p2perp[1] + p2z*p2z;
 double pcmsq = (p1perp[0]+p2perp[0])*(p1perp[0]+p2perp[0]) + (p1perp[1]+p2perp[1])*(p1perp[1]+p2perp[1]) + (p1z+p2z)*(p1z+p2z);
 
 double da1dp2z = -A/mHe4 * sqrt(mN*mN+p1sq) * ( p1z/(sqrt(mH2*mH2+pcmsq)*(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq))) + p2z/(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq)) * ( 1/sqrt(mH2*mH2+pcmsq) + (mHe4-sqrt(mH2*mH2+pcmsq))/(sqrt(mN*mN+p2sq)*(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq))) ) );

 return da1dp2z;
}

double taggedSF::da2dp1z(double *p1perp, double p1z, double *p2perp, double p2z){
 double p1sq = p1perp[0]*p1perp[0] + p1perp[1]*p1perp[1] + p1z*p1z;
 double p2sq = p2perp[0]*p2perp[0] + p2perp[1]*p2perp[1] + p2z*p2z;
 double pcmsq = (p1perp[0]+p2perp[0])*(p1perp[0]+p2perp[0]) + (p1perp[1]+p2perp[1])*(p1perp[1]+p2perp[1]) + (p1z+p2z)*(p1z+p2z);

 double da2dp1z = -A/mHe4 * sqrt(mN*mN+p2sq) * ( p2z/(sqrt(mH2*mH2+pcmsq)*(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq))) + p1z/(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq)) * ( 1/sqrt(mH2*mH2+pcmsq) + (mHe4-sqrt(mH2*mH2+pcmsq))/(sqrt(mN*mN+p1sq)*(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq))) ) );

 return da2dp1z;
}

double taggedSF::da2dp2z(double *p1perp, double p1z, double *p2perp, double p2z){
  double p1sq = p1perp[0]*p1perp[0] + p1perp[1]*p1perp[1] + p1z*p1z;
  double p2sq = p2perp[0]*p2perp[0] + p2perp[1]*p2perp[1] + p2z*p2z;
  double pcmsq = (p1perp[0]+p2perp[0])*(p1perp[0]+p2perp[0]) + (p1perp[1]+p2perp[1])*(p1perp[1]+p2perp[1]) + (p1z+p2z)*(p1z+p2z);
  double da2dp2z = A/mHe4*(1+p2z/(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq)) * ( (mHe4-sqrt(mH2*mH2+pcmsq)) / sqrt(mN*mN+p2sq)-sqrt(mN*mN+p2sq)/sqrt(mH2*mH2+pcmsq)-(mHe4-sqrt(mH2*mH2+pcmsq))/(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq)) ) - p1z/(sqrt(mN*mN+p1sq)+sqrt(mN*mN+p2sq))*sqrt(mN*mN+p2sq)/sqrt(mH2*mH2+pcmsq));

  return da2dp2z;
}


double taggedSF::J_pz_to_alpha(double alpha1, double *p1perp, double alpha2, double *p2perp, double p1z, double p2z){
    
  /*    double p1z_lo1, p1z_hi1, p2z_lo1, p2z_hi1; //used for alpha1 partials (instead of reusing variables)
  double p1z_lo2, p1z_hi2, p2z_lo2, p2z_hi2; //used for alpha2 partials

  double alpha1_lo = alpha1 - 0.005;
  double alpha1_hi = alpha1 + 0.005;

  double alpha2_lo = alpha2 - 0.005;
  double alpha2_hi = alpha2 + 0.005;//can make .05 if needed
  
  Findpz(alpha1_lo, p1perp, alpha2, p2perp, p1z_lo1, p2z_lo1);
  Findpz(alpha1_hi, p1perp, alpha2, p2perp, p1z_hi1, p2z_hi1);

  Findpz(alpha1, p1perp, alpha2_lo, p2perp, p1z_lo2, p2z_lo2);
  Findpz(alpha1, p1perp, alpha2_hi, p2perp, p1z_hi2, p2z_hi2);
  
  double d_a1_d_p1z = ( alpha1_hi - alpha1_lo ) / ( p1z_hi1 - p1z_lo1 ) ;
  double d_a1_d_p2z = ( alpha1_hi - alpha1_lo ) / ( p2z_hi1 - p2z_lo1 ) ;
  
  double d_a2_d_p1z = ( alpha2_hi - alpha2_lo ) / ( p1z_hi2 - p1z_lo2 ) ;
  double d_a2_d_p2z = ( alpha2_hi - alpha2_lo ) / ( p2z_hi2 - p2z_lo2 ) ;*/
  ////////
  
  double da1_dp1z = da1dp1z(p1perp,p1z,p2perp,p2z);
  double da1_dp2z = da1dp2z(p1perp,p1z,p2perp,p2z);
  double da2_dp1z = da2dp1z(p1perp,p1z,p2perp,p2z);
  double da2_dp2z = da2dp2z(p1perp,p1z,p2perp,p2z);
  
  double det = abs(( da1_dp1z * da2_dp2z ) - ( da1_dp2z * da2_dp1z ));


  /*
  if (det>1e9){

    std::cout << "\n\n\n\nalpha1: " << alpha1 << std::endl;
    
    std::cout << "\nAnalytic Solutions: " << std::endl;
    std::cout << "da1/dp1z is " << da1_dp1z << std::endl;
    std::cout << "da1/dp2z is " << da1_dp2z << std::endl;
    std::cout << "da2/dp1z is " << da2_dp1z << std::endl;
    std::cout << "da2/dp2z is " << da2_dp2z << std::endl;
    
    std::cout << "\n Numerical Solutions: " << std::endl;
    std::cout << "da1/dp1z is " << d_a1_d_p1z << std::endl;
    std::cout << "da1/dp2z is " << d_a1_d_p2z << std::endl;
    std::cout << "da2/dp1z is " << d_a2_d_p1z << std::endl;
    std::cout << "da2/dp2z is " << d_a2_d_p2z << std::endl;
    
    std::cout << "\n\n\nVarying alpha1:" << std::endl;
    std::cout << "lo p1z: " << p1z_lo1 << std::endl;
    std::cout << "mid p1z: " << p1z << std::endl;
    std::cout << "hi p1z: " << p1z_hi1 << std::endl;
    
    std::cout << "\nlo p2z: " << p2z_lo1 << std::endl;
    std::cout << "mid p2z: " << p2z << std::endl;
    std::cout << "hi p2z: " << p2z_hi1 << std::endl;
    
    std::cout <<"\n\nVarying alpha2: " << std::endl;
    std::cout << "lo p1z: " << p1z_lo2 << std::endl;
    std::cout << "mid p1z: " << p1z << std::endl;
    std::cout << "hi p1z: " << p1z_hi2 << std::endl;
    
    std::cout << "\nlo p2z: " << p2z_lo2 << std::endl;
    std::cout << "mid p2z: " << p2z << std::endl;
    std::cout << "hi p2z: " << p2z_hi2 << std::endl;
    
  }
  */  
  /*if (det>1e9){
    std::cout << "\n\nda1/dp1z is " << d_a1_d_p1z << std::endl;
    std::cout << "da1/dp2z is " << d_a1_d_p2z << std::endl;
    std::cout << "da2/dp1z is " << d_a2_d_p1z << std::endl;
    std::cout << "da2/dp2z is " << d_a2_d_p2z << std::endl;

    std::cout << "\nalpha1 is " << alpha1 << std::endl;
    std::cout << "p1z is " << std::setprecision (std::numeric_limits<double>::digits10 + 1) << p1z << std::endl;
    std::cout << "p1z_hi1 is " << std::setprecision (std::numeric_limits<double>::digits10 + 1) << p1z_hi1 << std::endl;
    std::cout << "p1z_lo1 is " << std::setprecision (std::numeric_limits<double>::digits10 + 1) << p1z_lo1 << std::endl;
    std::cout << "p2z is " << std::setprecision (std::numeric_limits<double>::digits10 + 1) << p2z << std::endl;
    std::cout << "p2z_hi1 is " << std::setprecision (std::numeric_limits<double>::digits10 + 1) << p2z_hi1 << std::endl;
    std::cout << "p2z_lo1 is " << p2z_lo1 << std::endl;
    
    std::cout << "\nalpha2 is " << alpha2 << std::endl;
    std::cout << "p1z is " << p1z << std::endl;
    std::cout << "p1z_hi2 is " << p1z_hi2 << std::endl;
    std::cout << "p1z_lo2 is " << p1z_lo2 << std::endl;
    std::cout << "p2z is " << p2z << std::endl;
    std::cout << "p2z_hi2 is " << p2z_hi2 << std::endl;
    std::cout << "p2z_lo2 is " << p2z_lo2 << std::endl;
    }
  */
  return det;
}



double taggedSF::rho_pn(double alpha1, double *p1perp, double alpha2, double *p2perp, double p1z, double p2z, double Jacobian){

  double p1[3] = {p1perp[0], p1perp[1], p1z};
  double p2[3] = {p2perp[0], p2perp[1], p2z};

  double p1mag = sqrt( p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] ); 
  double p2mag = sqrt( p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2] ); 
  
  double pdiff[3] = { p1[0] - p2[0] ,  p1[1] - p2[1] , p1[2] - p2[2] };
  double pdiff_mag = sqrt(pdiff[0]*pdiff[0] +  pdiff[1]*pdiff[1] + pdiff[2]*pdiff[2]);

  double q = 0.5 * pdiff_mag;

  double psum[3] = { p1[0] + p2[0] ,  p1[1] + p2[1] , p1[2] + p2[2] };

  double Q = sqrt(psum[0]*psum[0] +  psum[1]*psum[1] + psum[2]*psum[2]);

  
  double rho_np = ( 1 - ( n_MF(p1mag) / n_tot(p1mag) ) ) * ( 1 - ( n_MF(p2mag) / n_tot(p2mag) ) ) * ( 1 / (A - 1) ) * n_pn( q , Q ) * Jacobian;  

  return rho_np;
}


double taggedSF::rho_nn(double alpha1, double *p1perp, double alpha2, double *p2perp, double p1z, double p2z, double Jacobian){

  double p1[3] = {p1perp[0], p1perp[1], p1z};
  double p2[3] = {p2perp[0], p2perp[1], p2z};

  double p1mag = sqrt( p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] ); 
  double p2mag = sqrt( p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2] ); 
  
  double pdiff[3] = { p1[0] - p2[0] ,  p1[1] - p2[1] , p1[2] - p2[2] };
  double pdiff_mag = sqrt(pdiff[0]*pdiff[0] +  pdiff[1]*pdiff[1] + pdiff[2]*pdiff[2]);

  double q = 0.5 * pdiff_mag;

  double psum[3] = { p1[0] + p2[0] ,  p1[1] + p2[1] , p1[2] + p2[2] };

  double Q = sqrt(psum[0]*psum[0] +  psum[1]*psum[1] + psum[2]*psum[2]);

  
  double rho_nn = ( 1 - ( n_MF(p1mag) / n_tot(p1mag) ) ) * ( 1 - ( n_MF(p2mag) / n_tot(p2mag) ) ) * ( 2 / (A - 1) ) * n_nn( q , Q ) * Jacobian;  

  return rho_nn;
}

double taggedSF::F2p(double x, double Qsq){

  double m02 = 0.5063;
  double mP2 = 34.75;
  double mR2 = 0.0319;
  double Q02 = 1.374 ;
  double Lam02 = 0.06527;
  double aP1 = -0.11895;
  double aP2 = -0.4783;
  double aP3 = 1.353;
  double bP1 = 1.0833;
  double bP2 = 2.656;
  double bP3 = 1.771;
  double cP1 = 0.3638;
  double cP2 = 0.1211;
  double cP3 = 1.166;
  double aR1 = 0.3425;
  double aR2 = 1.0603;
  double aR3 = 0.5164;
  double bR1 = -10.408;
  double bR2 = 14.857;
  double bR3 = 0.07739;
  double cR1 = 1.3633;
  double cR2 = 2.256;
  double cR3 = 2.209;

  double t = log( log( (Qsq + Q02) / Lam02 ) / log( Q02 / Lam02 ) );

  double aP = aP1 + ( aP1 - aP2 ) * ( 1 / (1 + pow(t,aP3)) - 1 );
  double bP = bP1 + bP2 * pow(t,bP3);
  double cP = cP1 + ( cP1 - cP2 ) * ( 1 / (1+pow(t,cP3)) - 1 );

  double aR = aR1 + aR2 * pow(t,aR3);
  double bR = bR1 + bR2 * pow(t,bR3);
  double cR = cR1 + cR2 * pow(t,cR3);

  double Wsq = Qsq * ( 1/x - 1 ) + mN*mN;

  double xP = 1 / ( 1 + ( Wsq - mN*mN ) / ( Qsq + mP2 ) );
  double xR = 1 / ( 1 + ( Wsq - mN*mN ) / ( Qsq + mR2 ) );

  double F2P = cP * pow(xP,aP) * pow( (1 - x), bP); //P/p refers to Pomeron within this fn, outsid \
e of this specific fn, p refers to proton                                                          
  double F2R = cR * pow(xR,aR) * pow( (1 - x), bR);

  double F2 = ( Qsq / ( Qsq + m02 ) ) * ( F2P + F2R );

  return F2;

}

double taggedSF::Rpn(double x){

  double a=0.57;
  double b=2.2;
  double c=0.42;
  
  return a*pow(1-x,b)+c;
}


double taggedSF::taggedF2(double xB, double Qsq, double alpha2, double *p2perp, double &integral_err){

  double integral = 0;
  double weightsq = 0;
  long int n = 1.5e8; //increase once working

  long int breakcount = 0;
  const double alpha_min = xB*mN/mbar;
  // const double alpha_max = 1.8;
  //  const double alpha_max = A;
  const double alpha_max = 4.-alpha2;

  //add all (most?) variables in loop as cosntants here)
  
  for (long int i = 0 ; i < n ; i++ ){
    if (i%15000000 == 0){
      //std::cout << "i is " << i << ", n is " << n << std::endl;
      std::cerr << ((double) i)/((double)n)*100 << "% done\n";
      std::cerr << "if stopped here, F2tag would be " << integral/((double)i)  << " +/- " << sqrt(weightsq)/((double)i) << std::endl;
    }


    
    //    double alpha1 = myRand->Uniform(alpha_min,alpha_max);
    double alpha1 = 5;
    while ( (alpha1 > alpha_max) or (alpha1 < alpha_min)){
      alpha1 = myRand->Gaus(1,0.2);
    }
    double xtil = xB*mN/(alpha1*mbar);

    double norm = ROOT::Math::normal_cdf(alpha_max,0.2,1.)-ROOT::Math::normal_cdf(alpha_min,0.2,1.);
    double pdf_alpha = 1/norm * TMath::Gaus(alpha1 , 1. , .2, true);

    double p1x = myRand->Gaus(0,0.140);
    double p1y = myRand->Gaus(0,0.140);


    double p1perp[2] = {p1x,p1y};

    double p1z, p2z;
    Findpz(alpha1, p1perp, alpha2, p2perp, p1z, p2z);

    //adding this in as safeguard for when pz bisection fails
    /*if (p1z >= 10. or p1z <= -10.){
      breakcount+=1;
	continue;
      exit(-2);
      }*/

    double jacobian = J_pz_to_alpha(alpha1, p1perp, alpha2, p2perp, p1z, p2z);
    
    double p1[3] = {p1perp[0], p1perp[1], p1z};
    double p2[3] = {p2perp[0], p2perp[1], p2z};

    double p1mag = sqrt( p1[0]*p1[0] + p1[1]*p1[1] + p1[2]*p1[2] );
    double p2mag = sqrt( p2[0]*p2[0] + p2[1]*p2[1] + p2[2]*p2[2] );

    double pdiff[3] = { p1[0] - p2[0] ,  p1[1] - p2[1] , p1[2] - p2[2] };
    double pdiff_mag = sqrt(pdiff[0]*pdiff[0] +  pdiff[1]*pdiff[1] + pdiff[2]*pdiff[2]);

    double q = 0.5 * pdiff_mag;

    double psum[3] = { p1[0] + p2[0] ,  p1[1] + p2[1] , p1[2] + p2[2] };

    double Q = sqrt(psum[0]*psum[0] +  psum[1]*psum[1] + psum[2]*psum[2]);


    //adding in virtuality even tho not adding in modification yet
    double pcmperp[2] = { p1perp[0]+p2perp[0] , p1perp[1]+p2perp[1] };
    double pmin = mHe4 - ( mN*mN + sqrt( p2perp[0]*p2perp[0] + p2perp[1]*p2perp[1] )) / (mbar*alpha2) - ( mH2*mH2 + sqrt( pcmperp[0]*pcmperp[0] + pcmperp[1]*pcmperp[1] )) / (mbar*(A-alpha1-alpha2));
    double v1 = (mbar*pmin*alpha1 - sqrt( p1perp[0]*p1perp[0] + p1perp[1]*p1perp[1] ) - mN*mN) / (mN*mN);
    //std::cout << "virtuality is " << v1 << std::endl;
    
    if ( (Q < 0.000197) or (Q > 0.985) ){
      continue;
    }

    if ( (q < 0.197) or (q > 1.182) ){
      continue;
    }

    if ( (p2mag < 0.004925) or (p2mag > 1.96507) ){
      continue;
    }

    if ( (p1mag < 0.004925) or (p1mag > 1.96507) ){
      continue;
    }

    
    //double pcmperp[2] = {p1perp[0]+p2perp[0], p1perp[1]+p2perp[1]};

    double pdf_p1x = TMath::Gaus(p1x, 0  ,.140, true);
    double pdf_p1y = TMath::Gaus(p1y, 0  ,.140, true);


    double pdf_total = pdf_alpha * pdf_p1x * pdf_p1y;

    double F = 1/(A*alpha2) * F2p(xtil, Qsq) * (Z*rho_pn(alpha1, p1perp, alpha2, p2perp, p1z, p2z, jacobian) + N*Rpn(xtil)*rho_nn(alpha1, p1perp, alpha2, p2perp, p1z, p2z, jacobian)); //no modification for now

    
    double w =  F/(alpha1*pdf_total);

    if ( fabs(w) > 100000000 ){
      std::cerr << "\nalpha1 is " << alpha1 << std::endl;
      std::cerr << "p1 is { " << p1x << " , " << p1y << " , " << p1z << " } with mag " << p1mag << std::endl;
      std::cerr << "p2 is { " << p2perp[0] << " , " << p2perp[1] << " , " << p2z << " } with mag " << p2mag << std::endl;
      std::cerr << "alpha1's pdf is " << pdf_alpha << ", p1x's pdf is " << pdf_p1x << ", p1y's pdf is " << pdf_p1y << ", and the total pdf is " << pdf_total << std::endl;
      std::cerr << "virtuality is " << v1 << std::endl;
      std::cerr << "F2p(xtil, Qsq) is " << F2p(xtil, Qsq) << std::endl;
      std::cerr << "Rpn(xtil) is " << Rpn(xtil) << std::endl;
      std::cerr << "J_pz_to_alpha( alpha1, p1perp, alpha2, p2perp) is " << jacobian << std::endl;
      std::cerr << "The change to F2A is " << (integral + w)/(((double)i)+1) - integral/((double)i) << std::endl;
      }
    
    
    integral += w;
    weightsq += w*w;

    double var_sofar = weightsq/((double)i) - integral*integral/(i * i);
    /*    if ( fabs(F/(alpha1*pdf_total)) > 100000000 ){
      //if ( (i>1e5) and  (w > (integral/i + 10*sqrt(var_sofar))) ){
      std::cout << "\n1/(A*alpha2) is " << 1/(A*alpha2) << std::endl;
      std::cout << "F2p(xtil, Qsq) is " << F2p(xtil, Qsq) << std::endl;
      std::cout << "rho_pn(alpha1, p1perp, alpha2, p2perp) is " << rho_pn(alpha1, p1perp, alpha2, p2perp) << std::endl;
      std::cout << "Rpn(xtil) is " << Rpn(xtil) << std::endl;
      std::cout << "rho_nn(alpha1, p1perp, alpha2, p2perp)) is " << rho_nn(alpha1, p1perp, alpha2, p2perp) << std::endl;
      std::cout << "n_MF(p1mag) is " << n_MF(p1mag) << " and  n_MF(p2mag) is " << n_MF(p2mag) << std::endl;
      std::cout << " n_tot(p1mag) is " <<  n_tot(p1mag) << " and  n_tot(p2mag) is " << n_tot(p2mag) << std::endl;
      std::cout << "n_pn( q , Q ) is " << n_pn( q , Q ) << std::endl;
      std::cout << "n_nn( q , Q ) is " << n_nn( q , Q ) <<std::endl;
      std::cout << "integral addition " << F/(alpha1*pdf_total) << std::endl;
      std::cout << "alpha1 is " << alpha1 << std:: endl;
      std::cout << "pdf_total is " << pdf_total << std::endl;
      std::cout << "J_pz_to_alpha( alpha1, p1perp, alpha2, p2perp) is " << J_pz_to_alpha( alpha1, p1perp, alpha2, p2perp) << std::endl;
      }*/
  }
   

  integral = integral/((double)n);
  //weightsq = weightsq/((double)n);
  integral_err = sqrt(weightsq)/((double)n);
  std::cerr << "Integral +- int error: " << integral << " +/- " << integral_err << std::endl;
  std::cerr << "Findpz failed " << breakcount << "times" <<std::endl;
  return integral;
    //return 0;

}
