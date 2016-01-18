// $Id: gammaDistribution.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___GAMMA_DIST
#define ___GAMMA_DIST
/************************************************************
This distribution can take several forms depending on its free parameter alpha 
(beta is assumed to be equal to alpha). For an extensive exlpanation of this distribution
see http://mathworld.wolfram.com/GammaDistribution.html.
please note that the borders of the categories are defined according to calculation of 
the gamma integral, according to numerical recipes in gammaUtilities
_globalRate represents the rate for two joint genes.
************************************************************/
#include "definitions.h"
#include "distribution.h"
class gammaDistribution : public distribution {

public:
	explicit gammaDistribution(MDOUBLE alpha=1,int in_number_of_categories=1);
	explicit gammaDistribution(const gammaDistribution& other);
	virtual ~gammaDistribution() {};
	virtual const int categories() const {return _rates.size();}
	virtual const MDOUBLE rates(const int i) const {return _rates[i]*_globalRate;}
	virtual const MDOUBLE ratesProb(const int i) const {return _ratesProb[i];}
	virtual distribution* clone() const { return new gammaDistribution(*this); }
 	virtual void setGlobalRate(const MDOUBLE x) {_globalRate = x;}
 	virtual MDOUBLE getGlobalRate()const {return _globalRate;}
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	virtual void setAlpha(MDOUBLE newAlpha);
	virtual MDOUBLE getAlpha() const {return _alpha;};
	virtual void change_number_of_categories(int in_number_of_categories);
	virtual void setGammaParameters(int numOfCategories=1 ,MDOUBLE alpha=1);
	virtual MDOUBLE getBorder(const int i) const {return _bonderi[i];}	//return the ith border. Note:  _bonderi[0] = 0, _bondery[categories()] = infinite

private:	
	int fill_mean();
	int fill_bonderi();

	MDOUBLE _alpha;
	vector<MDOUBLE> _bonderi;
	vector<MDOUBLE> _rates;
	vector<MDOUBLE> _ratesProb;
	MDOUBLE _globalRate;
};
#endif
