// $Id: generalGammaDistribution.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___GENERAL_GAMMA_DIST
#define ___GENERAL_GAMMA_DIST
/************************************************************
This distribution can take several forms depending on its free parameters alpha,beta 
(unalike gammaDist. alpha is not necessarily equal to beta). 
For an extensive exlpanation of this distribution
see http://mathworld.wolfram.com/GammaDistribution.html
************************************************************/
#include "definitions.h"
#include "distribution.h"

enum quadratureType {QUANTILE, LAGUERRE};

class generalGammaDistribution : public distribution {

public:
	explicit generalGammaDistribution(MDOUBLE alpha, MDOUBLE beta, int in_number_of_categories);
	explicit generalGammaDistribution(const generalGammaDistribution& other);
	explicit generalGammaDistribution();
	virtual ~generalGammaDistribution();
	virtual void setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta);

	virtual const int categories() const {return _rates.size();}
	virtual const MDOUBLE rates(const int i) const {return _rates[i]*_globalRate;}
	virtual const MDOUBLE ratesProb(const int i) const {return _ratesProb[i];}
	virtual distribution* clone() const { return new generalGammaDistribution(*this); }
 	virtual void setGlobalRate(const MDOUBLE x) {_globalRate = x;}
 	virtual MDOUBLE getGlobalRate()const {return _globalRate;}
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	virtual void setAlpha(MDOUBLE newAlpha);
	virtual MDOUBLE getAlpha() const {return _alpha;};
	virtual void setBeta(MDOUBLE newBeta);
	virtual MDOUBLE getBeta() const {return _beta;};
	virtual void change_number_of_categories(int in_number_of_categories);
	virtual MDOUBLE getBorder(const int i) const {return _bonderi[i];}	//return the ith border. Note:  _bonderi[0] = 0, _bondery[categories()] = infinite


private:	
	int fill_mean();
	int fill_bonderi();
	
	vector<MDOUBLE> _bonderi;
protected:
	MDOUBLE _alpha;
	MDOUBLE _beta;
	
	vector<MDOUBLE> _rates;
	vector<MDOUBLE> _ratesProb;
	MDOUBLE _globalRate;
};



#endif

