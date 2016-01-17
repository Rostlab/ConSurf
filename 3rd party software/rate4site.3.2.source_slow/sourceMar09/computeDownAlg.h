// $Id: computeDownAlg.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___COMPUTE_DOWN_ALG
#define ___COMPUTE_DOWN_ALG

#include "definitions.h"
#include "tree.h"
#include "suffStatComponent.h"
#include "sequenceContainer.h"
#include "computePijComponent.h"


class computeDownAlg {
public: 
	void fillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup);

	void fillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup);

	void fillComputeDownSpecificRate(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup,
					   const MDOUBLE gRate);

};
#endif
