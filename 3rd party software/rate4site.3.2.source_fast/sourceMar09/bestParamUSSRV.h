// 	$Id: bestParamUSSRV.h 920 2006-09-21 09:26:12Z ninio $	
#ifndef BEST_PARAM_USSRV
#define BEST_PARAM_USSRV

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
#include "bestAlphaAndNu.h"
#include "bblEM2USSRV.h"

class bestParamUSSRV
{
public:
	explicit bestParamUSSRV(bool AlphaOptimization, bool NuOptimization,
							bool FOptimization, bool bblOptimization):
							_AlphaOptimizationFlag(AlphaOptimization),
							_NuOptimizationFlag(NuOptimization),
							_FOptimizationFlag(FOptimization),
							_bblOptimizationFlag(bblOptimization) {}
	
	MDOUBLE operator() (tree& et,
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights=NULL,
					   const MDOUBLE AlphaUpperBound = 15, 
					   const MDOUBLE NuUpperBound = 15, 
					   const MDOUBLE FUpperBound = 1, 
					   const MDOUBLE epsilonParamOptimization = 0.01,
					   const MDOUBLE epsilonLikelihoodImprovment = 0.05,
					   const int maxIterations = 10);

	MDOUBLE getBestAlpha() {return _bestAlpha;}
	MDOUBLE getBestNu() {return _bestNu;}
	MDOUBLE getBestF() {return _bestF;}
	MDOUBLE getBestL() {return _bestL;}

private:
	MDOUBLE _bestAlpha;
	MDOUBLE _bestNu;
	MDOUBLE _bestF;
	MDOUBLE _bestL;

	// flags
	bool _AlphaOptimizationFlag;
	bool _NuOptimizationFlag;
    bool _FOptimizationFlag;
	bool _bblOptimizationFlag;
};

#endif // BEST_PARAM_USSRV

