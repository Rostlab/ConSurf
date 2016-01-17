// $Id: bestAlpha.h 749 2006-06-29 11:44:32Z ninio $

#ifndef ___BEST_ALPHA
#define ___BEST_ALPHA

#include "definitions.h"

#include "likelihoodComputation.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "logFile.h"

#ifndef VERBOS
#define VERBOS
#endif

class bestAlphaFixedTree {
public:
	explicit bestAlphaFixedTree(const tree& et,
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnAlpha = 15,
					   const MDOUBLE epsilonAlphaOptimization = 0.01);
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

class bestAlphaAndBBL {
public:
	explicit bestAlphaAndBBL(tree& et, //find Best Alpha and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights=NULL,
					   const MDOUBLE initAlpha = 1.5,
					   const MDOUBLE upperBoundOnAlpha = 5.0,
					   const MDOUBLE epsilonAlphaOptimization= 0.01,
					   const MDOUBLE epsilonLikelihoodImprovment= 0.05,
					   const int maxBBLIterations=10,
					   const int maxTotalIterations=5);
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};




class C_evalAlpha{
public:
  C_evalAlpha(	const tree& et,
				const sequenceContainer& sc,
				stochasticProcess& sp,
				const Vdouble * weights = NULL)
    : _et(et),_sc(sc),_weights(weights),_sp(sp){};
private:
	const tree& _et;
	const sequenceContainer& _sc;
	const Vdouble * _weights;
	stochasticProcess& _sp;
public:
	MDOUBLE operator() (MDOUBLE alpha) {
		if (_sp.categories() == 1) {
			errorMsg::reportError(" one category when trying to optimize alpha");
		}
		(static_cast<gammaDistribution*>(_sp.distr()))->setAlpha(alpha);
		
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		//LOG(5,<<" with alpha = "<<alpha<<" logL = "<<res<<endl);
#ifdef VERBOS
		LOG(7,<<" while in brent: with alpha = "<<alpha<<" logL = "<<res<<endl);
#endif
		return -res;
	}
};

#endif


