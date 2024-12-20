#ifndef __TAGGED_SF_H__
#define __TAGGED_SF_H__

class TRandom3;
class GridInterp;
class LineInterp;

class taggedSF{

 public:
  taggedSF();
  ~taggedSF();
  double taggedF2(double xB, double Qsq, double alpha2, double *p2perp);
  void Findpz(double alpha1, double *p1perp, double alpha2, double *p2perp, double &p1z, double &p2z);
  double da1dp1z(double *p1perp, double p1z, double *p2perp, double p2z);
  double da1dp2z(double *p1perp, double p1z, double *p2perp, double p2z);
  double da2dp1z(double *p1perp, double p1z, double *p2perp, double p2z);
  double da2dp2z(double *p1perp, double p1z, double *p2perp, double p2z);
  /*double rho_nn(double alpha1, double *p1perp, double alpha2, double *p2perp);
  double rho_pn(double alpha1, double *p1perp, double alpha2, double *p2perp);
  double n_nn(double q, double Q);
  double n_pn(double q, double Q);*/
 private:
  TRandom3 * myRand;
  GridInterp * n_pn_interp;
  GridInterp * n_nn_interp;
  LineInterp * n_MF_interp;
  LineInterp * n_tot_interp;
  //void Findpz(double alpha1, double *p1perp, double alpha2, double *p2perp, double &p1z, double &p2z);
  double n_pn(double q, double Q);
  double n_nn(double q, double Q);
  double n_MF(double p);
  double n_tot(double p);
  double J_pz_to_alpha(double alpha1, double *p1perp, double alpha2, double *p2perp);
  double rho_pn(double alpha1, double *p1perp, double alpha2, double *p2perp);
  double rho_nn(double alpha1, double *p1perp, double alpha2, double *p2perp);
  double F2p(double x, double Qsq);
  double Rpn(double x);
  
};

#endif
