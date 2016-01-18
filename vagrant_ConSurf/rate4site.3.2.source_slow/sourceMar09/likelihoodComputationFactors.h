// $Id: likelihoodComputationFactors.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___LIKELIHOOD_COMPUTATION_FACTORS
#define ___LIKELIHOOD_COMPUTATION_FACTORS

#include "definitions.h"
#include "tree.h"
#include "computePijComponent.h"
#include "sequenceContainer.h"
#include "suffStatComponent.h"

namespace likelihoodComputation {

	MDOUBLE getLOG_LofPos(const int pos, // with a site specific rate.
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp,
					  const MDOUBLE gRate);

	// add all the other functions to use factors...


};



#endif

