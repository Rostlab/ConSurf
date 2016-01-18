// $Id: bestAlpha.cpp 749 2006-06-29 11:44:32Z ninio $

#include <iostream>
using namespace std;

#include "bestAlpha.h"
#include "bblEM.h"
#include "numRec.h"
#include "logFile.h"

#ifndef VERBOS
#define VERBOS
#endif
//void bestAlpha::checkAllocation() {
//	if (_pi->stocProcessFromLabel(0)->getPijAccelerator() == NULL) {
//		errorMsg::reportError(" error in function findBestAlpha");
//	}
//}
//
// @@@@ The method works with oldL,oldA,bestA and newL,newA.
// Only when it's about to end, the members _bestAlpha and _bestL are filled.

bestAlphaAndBBL::bestAlphaAndBBL(tree& et, //find Best Alpha and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights,
					   const MDOUBLE initAlpha,
					   const MDOUBLE upperBoundOnAlpha,
					   const MDOUBLE epsilonLoglikelihoodForAlphaOptimization,
					   const MDOUBLE epsilonLoglikelihoodForBBL,
					   const int maxBBLIterations,
					   const int maxTotalIterations){
//	LOG(5,<<"find Best Alpha and best BBL"<<endl);
//	LOG(5,<<" 1. bestAlpha::findBestAlpha"<<endl);
//	brLenOpt br1(*et,*pi,weights);
	
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;
	const MDOUBLE bx=initAlpha;
	const MDOUBLE ax=0;
	const MDOUBLE cx=upperBoundOnAlpha;
//
	MDOUBLE bestA=0;
	MDOUBLE oldA=0;
	int i=0;
	for (i=0; i < maxTotalIterations; ++i) {
		newL = -brent(ax,bx,cx,
		C_evalAlpha(et,sc,sp,weights),
		epsilonLoglikelihoodForAlphaOptimization,
		&bestA);
 
#ifdef VERBOS
		LOG(5,<<"# bestAlphaAndBBL::bestAlphaAndBBL iteration " << i <<endl
		      <<"# old L = " << oldL << "\t"
		      <<"# new L = " << newL << endl
		      <<"# new Alpha = " << bestA << endl);
#endif
		if (newL > oldL+epsilonLoglikelihoodForBBL) {
		    oldL = newL;
		    oldA = bestA;
		} else {
		    oldL = newL;
		    oldA = bestA;

		    
		    _bestL = oldL;
		    _bestAlpha= oldA;
		    (static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
		    break;
		}

		(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
		bblEM bblEM1(et,sc,sp,NULL,maxBBLIterations,epsilonLoglikelihoodForBBL);//maxIterations=1000
		newL =bblEM1.getTreeLikelihood();
#ifdef VERBOS
		LOG(5,<<"# bestAlphaAndBBL::bestAlphaAndBBL iteration " << i <<endl 
		      <<"# After BBL new L = "<<newL<<" old L = "<<oldL<<endl
		      <<"# The tree:" );
		LOGDO(5,et.output(myLog::LogFile()));
#endif

		if (newL > oldL+epsilonLoglikelihoodForBBL) {
			oldL = newL;
		}
		else {
		    oldL=newL;
		    _bestL = oldL;	
			_bestAlpha= oldA; 
			(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
			break;
		}
	}
	if (i==maxTotalIterations) {
		_bestL = newL;
		_bestAlpha= bestA;
		(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
	}
}
		
bestAlphaFixedTree::bestAlphaFixedTree(const tree& et, //findBestAlphaFixedTree
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnAlpha,
					   const MDOUBLE epsilonLoglikelihoodForAlphaOptimization){
	//LOG(5,<<"findBestAlphaFixedTree"<<endl);
	MDOUBLE bestA=0;
	const MDOUBLE cx=upperBoundOnAlpha;// left, midle, right limit on alpha
	const MDOUBLE bx=static_cast<gammaDistribution*>(sp.distr())->getAlpha();
	const MDOUBLE ax=0.0;

	
	_bestL = -brent(ax,bx,cx,
		C_evalAlpha(et,sc,sp,weights),
		epsilonLoglikelihoodForAlphaOptimization,
		&bestA);
	(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(bestA);
	_bestAlpha= bestA;
}

