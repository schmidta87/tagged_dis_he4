#ifndef __TAGGED_SF_H__
#define __TAGGED_SF_H__

class TRandom3;

class taggedSF{

 public:
  taggedSF();
  ~taggedSF();
  double taggedF2();
  
 private:
  TRandom3 * myRand;
  
};


#endif
