// $Id: generalGammaDistributionLaguerre.h 391 2005-06-12 12:44:57Z ninio $
// version 1.00
// last modified Sep 2004

#ifndef ___GENERAL_GAMMA_DIST_LAGUERRE
#define ___GENERAL_GAMMA_DIST_LAGUERRE
/************************************************************
This class differ from the regular generalGammaDistribution in that 
the rateCategories and their probabilities are not constructed using Yang's quantile method.
Instead the general Guass-Laguerre quadrature method is used.
For example, if we want to compute the likelihood over the rate distribution, 
then we need to solve the integral

I[0_to_infinity]{P(data|r)*P(r)} 
	= I[0_to_infinity]{P(data|r)*b^a / Gamma(a)* exp(-b*r) * r^(a-1)dr}  //a = alpha, b = beta
	= b^(a)/Gamma(a) * I[0_to_infinity]{P(data|m/b) * exp(-m) * (m/b)^(a')/bdm}  ///substitute m=b*r, a'=a-1
	= 1/Gamma(a) * I[0_to_infinity]{P(data|m/b) * exp(-m) * m^a' dm}  //
Now - we can use the Guass-Laguerre formula, to get an approximation for the Integral.
The Xj and Wj are the absicassas and weights of the Laguerre polynoms
	= 1/Gamma(a) * sum[j = 0_to_catNum]{P(data|Xj/b) * Wj}  
  
The rates are the Xj/b and their priors is Wj/Gamma(a) 
The quadrature method is explained in Numerical Recipes (Press et al.; chapter 4.5) 
and is also mentioned in Felsenstein 2001 (JME 53: 447-455).
************************************************************/
#include "definitions.h"
#include "generalGammaDistribution.h"
class generalGammaDistributionLaguerre : public generalGammaDistribution {

public:
	explicit generalGammaDistributionLaguerre(MDOUBLE alpha, MDOUBLE beta, int in_number_of_categories);
	explicit generalGammaDistributionLaguerre(const generalGammaDistributionLaguerre& other);
	virtual ~generalGammaDistributionLaguerre();
	void setGammaParameters(int numOfCategories ,MDOUBLE alpha, MDOUBLE beta);

	virtual const int categories() const {return _rates.size();}
	virtual const MDOUBLE rates(const int i) const {return _rates[i]*_globalRate;}
	virtual const MDOUBLE ratesProb(const int i) const {return _ratesProb[i];}
	virtual distribution* clone() const { return new generalGammaDistributionLaguerre(*this); }
 	virtual void setGlobalRate(const MDOUBLE x) {_globalRate = x;}
 	virtual MDOUBLE getGlobalRate()const {return _globalRate;}
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const;
	//void setAlpha(MDOUBLE newAlpha);
	//MDOUBLE getAlpha() const {return _alpha;};
	//void setBeta(MDOUBLE newBeta);
	//MDOUBLE getBeta() const {return _beta;};
	//void change_number_of_categories(int in_number_of_categories);
	MDOUBLE getBorder(const int i) const; 

protected:
	virtual void fillRatesAndProbs(int catNum);

/*private:	
	MDOUBLE _alpha;
	MDOUBLE _beta;
	vector<MDOUBLE> _rates;
	vector<MDOUBLE> _ratesProb;
	MDOUBLE _globalRate;
*/
};



#endif

