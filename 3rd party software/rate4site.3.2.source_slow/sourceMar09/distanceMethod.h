// $Id: distanceMethod.h 759 2006-07-03 16:28:17Z ninio $

#ifndef ___DISTANCE_METHOD
#define ___DISTANCE_METHOD
#include "definitions.h"
#include "sequence.h"

/*********************************************************
Distance method is a class for computing pairwise distance 
between 2 different sequences
*******************************************************/
class distanceMethod {
public:
  virtual const MDOUBLE giveDistance(const sequence& s1,
				     const sequence& s2,
				     const vector<MDOUBLE> * weights=NULL,
				     MDOUBLE* score=NULL) const=0;
  virtual distanceMethod* clone(void) const=0;
  virtual ~distanceMethod() {}
};


#endif

