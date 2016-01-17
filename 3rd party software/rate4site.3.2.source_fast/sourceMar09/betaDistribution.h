// $Id: betaDistribution.h 779 2006-07-19 10:32:22Z eyalprivman $

#ifndef ___BETA_DIST
#define ___BETA_DIST
/************************************************************
This distribution can take several forms depending on its free parameters alpha,beta 
For an extensive exlpanation of this distribution
see http://mathworld.wolfram.com/BetaDistribution.html
************************************************************/
#include "definitions.h"
#include "distribution.h"

class betaDistribution : public distribution {

public:
	explicit betaDistribution(MDOUBLE alpha, MDOUBLE beta, int in_number_of_categories);
	explicit betaDistribution(const betaDistribution& other);
	explicit betaDistribution();
	virtual ~betaDistribution();
	virtual void setBetaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta);

	virtual const int categories() const {return _rates.size();}
	virtual const MDOUBLE rates(const int i) const {return _rates[i]*_globalRate;}
	virtual const MDOUBLE ratesProb(const int i) const {return _ratesProb[i];}
	virtual distribution* clone() const { return new betaDistribution(*this); }
 	virtual void setGlobalRate(const MDOUBLE x) {_globalRate = x;}
 	virtual MDOUBLE getGlobalRate()const {return _globalRate;}
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	virtual void setAlpha(MDOUBLE newAlpha);
	virtual MDOUBLE getAlpha() const {return _alpha;};
	virtual void setBeta(MDOUBLE newBeta);
	virtual MDOUBLE getBeta() const {return _beta;};
	virtual void change_number_of_categories(int in_number_of_categories);
	virtual MDOUBLE getBorder(const int i) const {return _boundary[i];}	//return the ith border. Note:  _bonderi[0] = 0, _bondery[categories()] = infinite


private:	
	int fill_mean();
	int fill_boundaries();
	
	vector<MDOUBLE> _boundary;
protected:
	MDOUBLE _alpha;
	MDOUBLE _beta;
	
	vector<MDOUBLE> _rates;
	vector<MDOUBLE> _ratesProb;
	MDOUBLE _globalRate;
};



#endif

