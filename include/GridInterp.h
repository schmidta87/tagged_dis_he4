#ifndef __GRIDINTERP_H__
#define __GRIDINTERP_H__

class GridInterp{

 private:
  int nX;
  int nY;
  double minX;
  double maxX;
  double minY;
  double maxY;
  double *dataX;
  double *dataY;
  double *dataF;
  int findX(double val); // Get the index of the lower x-value grid point
  int findY(double val); // Get the index of the lower y-value grid point

 public:
  GridInterp();
  ~GridInterp();
  int getNX();
  int getNY();
  void setXVals(int N, const double *data);
  void setYVals(int N, const double *data);
  void setPoint(int ix, int iy, double val);
  double getVal(double x, double y);
};

#endif
