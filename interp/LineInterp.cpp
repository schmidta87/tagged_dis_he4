#include "LineInterp.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>

using namespace std;

LineInterp::LineInterp()
{
  nX=-1;
  minX=0.;
  maxX=0.;
  dataX=NULL;
  dataF=NULL;
}

LineInterp::~LineInterp()
{
  if (!dataX)
    delete[] dataX;
  if (!dataF)
    delete[] dataF;
}

int LineInterp::getNX()
{
  return nX;
}

void LineInterp::setXVals(int N, const double *data)
{
  if (!dataX)
    delete[] dataX;
  
  nX=N;

  if (nX <= 0)
    return; 

  minX=data[0];
  maxX=data[N-1];
   
  dataX = new double[N];
  for (int ix=0 ; ix <N ; ix++)
    {
      dataX[ix]=data[ix];
    }

  if (nX>0)
    {
      if (!dataF)
	delete[] dataF;

      dataF = new double[nX];
    }
}

void LineInterp::setPoint(int ix, double val)
{
  if ((ix<0)||(ix>=nX))
    {
      cerr << "Asked for " << ix << " in line of size " << nX << "\n"
	   << "\t\t Aborting...\n\n";
      exit(-1);
    }

  dataF[ix]=val;
}

double LineInterp::getVal(double x)
{
  if ((x<minX)||(x>maxX))
    {
      cerr << "Asked for " << x << " in line with bounds:\n"
	   << "\t X:\t" << minX << "\t" << maxX << "\n"
	   << "\t\t Aborting...\n\n";
      exit(-2);
    }  

  int ix=findX(x);


  //  NEED TO REDO THIS!!!
  double fXLo=dataF[ix];
  double fXHi=dataF[ix+1];


  return (fXLo * (dataX[ix+1]-x) + fXHi * (x-dataX[ix]))/(dataX[ix+1]-dataX[ix]);
}

// Yes, bisection would be faster, but to do later
int LineInterp::findX(double val)
{
  int i=0;
  while (i<nX-1)
    {
      if ((dataX[i]<val)&&(dataX[i+1]>val))
	return i;

      i++;
    }
  return i;
}

