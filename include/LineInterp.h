#ifndef __LINEINTERP_H__
#define __LINEINTERP_H__

class LineInterp{

 private:
  int nX;
  double minX;
  double maxX;
  double *dataX;
  double *dataF;
  int findX(double val); // Get the index of the lower x-value grid point

 public:
  LineInterp();
  ~LineInterp();
  int getNX();
  void setXVals(int N,const double *data);
  void setPoint(int ix, double val);
  double getVal(double x);
};

#endif
