// $Id: bestHKYparam.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___BEST_HKY_PARAM
#define ___BEST_HKY_PARAM

#include "definitions.h"

#include "likelihoodComputation.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "hky.h"


class bestHkyParamFixedTree {
public:
	explicit bestHkyParamFixedTree(const tree& et,
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnHkyParam = 0.5,
					   const MDOUBLE epsilonHkyParamOptimization = 0.01);
	MDOUBLE getBestHkyParam() {return _bestHkyParam;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestHkyParam;
	MDOUBLE _bestL;
};

class bestHkyParamAndBBL {
public:
	explicit bestHkyParamAndBBL(tree& et, //find Best HkyParam and best BBL
					   const sequenceContainer& sc,
					   stochasticProcess& sp,
					   const Vdouble * weights=NULL,
					   const MDOUBLE upperBoundOnHkyParam = 5.0,
					   const MDOUBLE epsilonHkyParamOptimization= 0.01,
					   const MDOUBLE epsilonLikelihoodImprovment= 0.05,
					   const int maxBBLIterations=10,
					   const int maxTotalIterations=5);
	MDOUBLE getBestHkyParam() {return _bestHkyParam;}
	MDOUBLE getBestL() {return _bestL;}
private:
	MDOUBLE _bestHkyParam;
	MDOUBLE _bestL;
};




class C_evalHkyParam{
public:
  C_evalHkyParam(	const tree& et,
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
	MDOUBLE operator() (MDOUBLE HkyParam) {
		(static_cast<hky*>(_sp.getPijAccelerator()->getReplacementModel()))->changeTrTv(HkyParam);
		
		MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_sp,_weights);
		//LOG(5,<<" with HkyParam = "<<HkyParam<<" logL = "<<res<<endl);
		return -res;
	}
};









#endif


