// $Id: granthamChemicalDistances.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___GRANTHAM_CHEMICAL_DISTANCES
#define ___GRANTHAM_CHEMICAL_DISTANCES

#include "definitions.h"

class granthamChemicalDistances {
public:
	explicit granthamChemicalDistances();
	MDOUBLE getGranthamDistance(const int aa1,const int aa2) const ;
	MDOUBLE getGranthamPolarityDistance(const int aa1,const int aa2) const;
	MDOUBLE getGranthamPolarity(const int aa1) const;
	virtual ~granthamChemicalDistances() {}

	MDOUBLE getHughesChargeDistance(const int aa1,const int aa2) const;// page 520
	MDOUBLE getHughesPolarityDistance(const int aa1,const int aa2) const;// page 520
	MDOUBLE getHughesHydrophobicityDistance(const int aa1,const int aa2) const;// page 520


private:

	// private members:
	MDOUBLE GranChemDist[20][20];
	MDOUBLE GranPolarityTable[20];

};


#endif


