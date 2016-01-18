// $Id: hky.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___HKY
#define ___HKY

#include "replacementModel.h"
#include <cmath>

class hky : public replacementModel {
public:
	explicit hky(const MDOUBLE inProb_a,
					const MDOUBLE inProb_c,
					const MDOUBLE inProb_g,
					const MDOUBLE inProb_t,
					const MDOUBLE TrTv);


	virtual replacementModel* clone() const { return new hky(*this); }
//	virtual nucJC* clone() const { return new nucJC(*this); } // see note down:

	const int alphabetSize() const {return 4;}


	void changeTrTv(const MDOUBLE In_TrTv);
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const;
	const MDOUBLE freq(const int i) const {return _freq[i];};
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const;

	const MDOUBLE dPij_tdBeta(const int i, const int j, const MDOUBLE t) const;

private:
	Vdouble _freq;
	MDOUBLE _a; //
	MDOUBLE _b; //

	MDOUBLE _c,_y; // relationship between probA, probC, prob G, prob T.
};

#endif

