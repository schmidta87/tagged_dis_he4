#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

//target info
const double A = 4;
const double Z = 2;
const double N = 2;

// Masses of various nuclei and particles
const double mU = 0.9315;
const double me = 0.000510998;
const double mP = 0.938272;
const double mN = 0.939565;
const double mH2 = 2.014102*mU-me;
const double mH3 = 3.0160492*mU-me;
const double mHe3 = 3.0160293*mU - 2.*me;
const double mHe4 = 4.0026033*mU - 2.*me;
const double mbar = mHe4/4;

// Constants of nature
const double alpha_em = 0.0072973525643;

//conversions
const double invfm_to_GeV = 0.197; 

#endif
