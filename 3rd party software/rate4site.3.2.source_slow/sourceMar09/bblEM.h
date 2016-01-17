// $Id: bblEM.h 645 2006-03-20 00:50:44Z eyalprivman $
#ifndef ___BBL_EM_H
#define ___BBL_EM_H

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include "countTableComponent.h"
#include "computePijComponent.h"
#include "suffStatComponent.h"
#include <vector>
using namespace std;

class bblEM {
public:
	explicit bblEM(tree& et,
				const sequenceContainer& sc,
				const stochasticProcess& sp,
				const Vdouble * weights = NULL,
				const int maxIterations=50,
				const MDOUBLE epsilon=0.05,
				const MDOUBLE tollForPairwiseDist=0.001);
	MDOUBLE getTreeLikelihood() const {return _treeLikelihood;}

private:
	MDOUBLE compute_bblEM(const int maxIterations,
					const MDOUBLE epsilon,
					const MDOUBLE tollForPairwiseDist);
	void bblEM_it(const MDOUBLE tollForPairwiseDist);
	void computeDown(const int pos);
	void computeUp();
	void addCounts(const int pos);
	void addCounts(const int pos, tree::nodeP mynode, const doubleRep posProb, const MDOUBLE weig);
	void optimizeBranches(const MDOUBLE tollForPairwiseDist);
	void allocatePlace();


	MDOUBLE _treeLikelihood;
	tree& _et;
	const sequenceContainer& _sc;
	const stochasticProcess& _sp;
	vector<countTableComponentGam> _computeCountsV; // for each node - a table of rate*alph*alph
	computePijGam _pij;
	suffStatGlobalGam _cup;
	suffStatGlobalGamPos _cdown;
	const Vdouble * _weights;
	VdoubleRep _posLike;
};

#endif
