// $Id: findRateOfGene.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ____FIND_RATE_OF_GENE
#define ____FIND_RATE_OF_GENE


#include "numRec.h"
#include "errorMsg.h"
#include "likelihoodComputation.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "suffStatComponent.h"
#include "definitions.h"

MDOUBLE findTheBestFactorFor(const tree &t,
							 const sequenceContainer& sc,
							stochasticProcess& sp,
							const Vdouble * weights,
							 MDOUBLE & logLresults);

void makeAverageRateEqOne(tree& et,vector<stochasticProcess> & spVec);

#endif
