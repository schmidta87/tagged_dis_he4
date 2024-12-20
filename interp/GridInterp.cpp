#include "GridInterp.h"
#include <iostream>
#include <cstdlib>
#include <cstdio>

using namespace std;

GridInterp::GridInterp()
{
  nX=-1;
  nY=-1;
  minX=0.;
  maxX=0.;
  minY=0.;
  maxY=0.;
  dataX=NULL;
  dataY=NULL;
  dataF=NULL;
}

GridInterp::~GridInterp()
{
  if (!dataX)
    delete[] dataX;
  if (!dataY)
    delete[] dataY;
  if (!dataF)
    delete[] dataF;
}

int GridInterp::getNX()
{
  return nX;
}

int GridInterp::getNY()
{
  return nY;
}

void GridInterp::setXVals(int N, const double *data)
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

  if ((nX>0) && (nY>0))
    {
      if (!dataF)
	delete[] dataF;

      dataF = new double[nX*nY];
    }
}

void GridInterp::setYVals(int N, const double *data)
{
  if (!dataY)
    delete[] dataY;
  
  nY=N;

  if (nY <= 0)
    return; 

  minY=data[0];
  maxY=data[N-1];

  dataY = new double[N];
  for (int iy=0 ; iy <N ; iy++)
    {
      dataY[iy]=data[iy];
    }

  if ((nY>0) && (nX>0))
    {
      if (!dataF)
	delete[] dataF;

      dataF = new double[nX*nY];
    }
}

void GridInterp::setPoint(int ix, int iy, double val)
{
  if ((ix<0)||(ix>=nX)||(iy<0)||(iy>=nY))
    {
      cerr << "Asked for " << ix << ", " << iy << " in grid of size " << nX << ", " << nY << "\n"
	   << "\t\t Aborting...\n\n";
      exit(-1);
    }

  dataF[iy*nX+ix]=val;
}

double GridInterp::getVal(double x, double y)
{
  if ((x<minX)||(x>maxX)||(y<minY)||(y>maxY))
    {
      cerr << "Asked for " << x << ", " << y << " in grid with bounds:\n"
	   << "\t X:\t" << minX << "\t" << maxX << "\n"
	   << "\t Y:\t" << minY << "\t" << maxY << "\n"
	   << "\t\t Aborting...\n\n";
      exit(-2);
    }  

  int ix=findX(x);
  int iy=findY(y);

  //cout << "\n ix is " << ix << " and iy is " << iy << endl;

  double fXLo=(dataX[ix+1]-x)/(dataX[ix+1]-dataX[ix]) * dataF[iy*nX+ix] + (x-dataX[ix])/(dataX[ix+1]-dataX[ix]) * dataF[iy*nX+ix+1];
  double fXHi=(dataX[ix+1]-x)/(dataX[ix+1]-dataX[ix]) * dataF[(iy+1)*nX+ix] + (x-dataX[ix])/(dataX[ix+1]-dataX[ix]) * dataF[(iy+1)*nX+ix+1];

  //cout << "fXLo is " << fXLo << " and fXHi is " << fXHi << endl;
  
  return fXLo * (dataY[iy+1]-y)/(dataY[iy+1]-dataY[iy]) + fXHi*(y-dataY[iy])/(dataY[iy+1]-dataY[iy]);
}

// Yes, bisection would be faster, but to do later
int GridInterp::findX(double val)
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

int GridInterp::findY(double val)
{
  int i=0;
  while (i<nY-1)
    {
      if ((dataY[i]<val)&&(dataY[i+1]>val))
	return i;

      i++;
    }
  return i;
}


