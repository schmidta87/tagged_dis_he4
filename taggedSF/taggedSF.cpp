#include "taggedSF.h"

#include "TRandom3.h"

taggedSF::taggedSF()
{
  myRand = new TRandom3(0);
}

taggedSF::~taggedSF()
{
  if (myRand)
    delete myRand;
}
 
double taggedSF::taggedF2()
{
  return 0.;
}
