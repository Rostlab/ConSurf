// 	$Id: bestAlphaAndNu.h 920 2006-09-21 09:26:12Z ninio $	
#ifndef ___BEST_ALPHA_AND_NU
#define ___BEST_ALPHA_AND_NU

#include "definitions.h"

#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "replacementModelSSRV.h"
#include "stochasticProcessSSRV.h"
#include "C_evalParamSSRV.h"
#include "bestAlpha.h"
#include "numRec.h"
#include "bblEM.h"
#include "logFile.h"


// Nu is fixed. The tree is fixed
class bestAlphaFixedTreeSSRV {
public:
	explicit bestAlphaFixedTreeSSRV() {}
	MDOUBLE operator()(const tree& et,
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnAlpha = 15,
					   const MDOUBLE epsilonAlphaOptimization = 0.01);
	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestL() {return _bestL;}
	
	void setAlpha(MDOUBLE alpha, ussrvModel& model) const  
	{
		model.updateAlpha(alpha);
	}

	void setBestL(MDOUBLE bestL) { _bestL = bestL;} 

private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestL;
};

// Alpha is fixed
class bestNuFixedTreeSSRV {
public:
	explicit bestNuFixedTreeSSRV(){}
	MDOUBLE operator()(const tree& et,
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnNu = 15,
					   const MDOUBLE epsilonNuOptimization = 0.01);
	MDOUBLE getBestNu() {return _bestNu;}
	MDOUBLE getBestL() {return _bestL;}
	void setNu(MDOUBLE nu, ussrvModel& model) const
	{
		model.updateNu(nu);
	}
	void setBestL(MDOUBLE bestL) { _bestL = bestL;}

private:
	MDOUBLE _bestNu;
	MDOUBLE _bestL;
};

class bestFFixedTreeSSRV {
public:
	explicit bestFFixedTreeSSRV() {}
	MDOUBLE operator()(const tree& et,
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnF = 1,
					   const MDOUBLE epsilonFOptimization = 0.01);
	MDOUBLE getBestF() {return _bestF;}
	MDOUBLE getBestL() {return _bestL;}
	void setF(MDOUBLE f, ussrvModel& model) const
	{
		if ( (f>1) || (f < 0))
		{
			LOG(5,<<"bestFFixedTreeSSRV:setF, f must be between 0 to 1. f = " << f << endl);
			return;
		}
		model.updateF(f);
	}
	void setBestL(MDOUBLE bestL) { _bestL = bestL;}

private:
	MDOUBLE _bestF;
	MDOUBLE _bestL;
};

#endif // ___BEST_ALPHA_AND_NU
