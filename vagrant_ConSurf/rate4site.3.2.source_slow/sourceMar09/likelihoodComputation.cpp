// $Id: likelihoodComputation.cpp 699 2006-05-27 16:39:58Z ninio $

#include "definitions.h"
#include "tree.h"
#include "computeUpAlg.h"
#include "likelihoodComputation.h"
#include "gammaUtilities.h"
#include <cmath>
#include <cassert>

using namespace likelihoodComputation;

doubleRep likelihoodComputation::getLofPos(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const computePijHom& pi,
					  const stochasticProcess& sp){

	computeUpAlg cup;
	suffStatGlobalHomPos ssc;
	cup.fillComputeUp(et,sc,pos,pi,ssc);

	doubleRep tmp = 0.0;
	for (int let = 0; let < sp.alphabetSize(); ++let) {
		doubleRep tmpLcat=
				ssc.get(et.getRoot()->id(),let)*
				sp.freq(let);
		assert(tmpLcat>=0.0);
		tmp+=tmpLcat;
	}
//	cout<<"likelihoodComputation::getLofPos: tmp = "; tmp.outputn(cout);	// DEBUG EP
	if (!(tmp>0.0)) errorMsg::reportError("likelihoodComputation::getLofPos: likelihood of pos was zero!",1);
	return tmp;
}

doubleRep likelihoodComputation::getProbOfPosWhenUpIsFilledHom(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp,
					  const suffStatGlobalHomPos& ssc){
// using the pij of stochastic process rather than pre computed pij's...
	if (ssc.size()==0) {errorMsg::reportError("error in function likelihoodComputation::getLofPosWhenUpIsFilled");}
	doubleRep tmp = 0.0;
	for (int let = 0; let < sp.alphabetSize(); ++let) {
		doubleRep tmpLcat=
				ssc.get(et.getRoot()->id(),let)*
				sp.freq(let);
		tmp+=tmpLcat;
	}
	return tmp;
}

doubleRep likelihoodComputation::getLofPosHomModelEachSiteDifferentRate(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp){
// using the pij of stochastic process rather than pre computed pij's...
	if (sp.categories()!=1) {
		  errorMsg::reportError("num of categories in function getLofPosHomModel must be one");
	}
	computeUpAlg cup;
	suffStatGlobalHomPos ssc;
	computePijHom cpij;
	cpij.fillPij(et,sp);
	cup.fillComputeUp(et,sc,pos,cpij,ssc);
	return getProbOfPosWhenUpIsFilledHom(pos,et,sc,sp,ssc);
}

doubleRep likelihoodComputation::getLofPosGamModelEachSiteDifferentRate(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp){
	computePijGam pi;
	pi.fillPij(et,sp);
	return getLofPos(pos,et,sc,pi,sp);
}

doubleRep likelihoodComputation::getLofPos(const int pos,
					  const tree& et,
					  const sequenceContainer& sc,
					  const stochasticProcess& sp,
					  const MDOUBLE gRate){ // when there is a global rate for this position
// using the pij of stochastic process rather than pre computed pij's...
	computeUpAlg cup;
	suffStatGlobalHomPos ssc;
	cup.fillComputeUpSpecificGlobalRate(et,sc,pos,sp,ssc,gRate);

	doubleRep tmp = 0.0;
	for (int let = 0; let < sp.alphabetSize(); ++let) {
		doubleRep tmpLcat=
				ssc.get(et.getRoot()->id(),let)*
				sp.freq(let);;
		assert(tmpLcat>=0.0);
		tmp+=tmpLcat;
	}
	return tmp;
}


doubleRep likelihoodComputation::getLofPos(const int pos,
										   const tree& et,
										   const sequenceContainer& sc,
										   const computePijGam& pi,
										   const stochasticProcess& sp){
//	with the pi already computed.
	doubleRep tmp=0;
	for (int i=0; i < sp.categories();++i) {
		tmp += getLofPos(pos,et,sc,pi[i],sp)*sp.ratesProb(i);
	}
	return tmp;
}

doubleRep likelihoodComputation::getLofPosAndPosteriorOfRates(const int pos,
															  const tree& et,
															  const sequenceContainer& sc,
															  const computePijGam& pi,
															  const stochasticProcess& sp,
															  VdoubleRep& postrior){
//	with the pi already computed.
	doubleRep tmp=0;
	for (int i=0; i < sp.categories();++i) {
	  postrior[i]=getLofPos(pos,et,sc,pi[i],sp)*sp.ratesProb(i);
	  tmp += postrior[i]; 
	}
	for (int i=0; i < sp.categories();++i) 
	  postrior[i] /= tmp;
	return tmp;
}

MDOUBLE likelihoodComputation::getTreeLikelihoodFromUp(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						const Vdouble * weights) {
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (int pos = 0; pos < sc.seqLen(); ++pos) {
		doubleRep tmp=0;
		for (int categor = 0; categor < sp.categories(); ++categor) {
			doubleRep veryTmp =0;
			for (int let =0; let < sc.getAlphabet()->size(); ++let) {
				veryTmp+=cup.get(pos,categor,et.getRoot()->id(),let) * sp.freq(let);
			}
			tmp += veryTmp*sp.ratesProb(categor);
		}
		like += log(tmp) * (weights?(*weights)[pos]:1);
	}
	return like;
}

MDOUBLE likelihoodComputation::getTreeLikelihoodFromUp2(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						VdoubleRep& posLike, // fill this vector with each position likelihood but without the weights.
						const Vdouble * weights) {
	posLike.clear();
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (int pos = 0; pos < sc.seqLen(); ++pos) {
		doubleRep tmp=0;
		for (int categor = 0; categor < sp.categories(); ++categor) {
			doubleRep veryTmp =0;
			for (int let =0; let < sc.alphabetSize(); ++let) {
				veryTmp+=cup.get(pos,categor,et.getRoot()->id(),let) * sp.freq(let);
			}
			tmp += veryTmp*sp.ratesProb(categor);
		}
		assert(tmp>0.0);
		like += log(tmp) * (weights?(*weights)[pos]:1);
		posLike.push_back(tmp);
	}
	return like;
}

// fill the posteriorLike matrix with each position posterior rate (p(r|D))
// but without the weights.
MDOUBLE likelihoodComputation::getPosteriorOfRates(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						VVdoubleRep& posteriorLike, 
						const Vdouble * weights) {
	suffStatGlobalGam cup;
	computeUpAlg cupAlg;
	computePijGam cpGam;
	cpGam.fillPij(et,sp);
	cupAlg.fillComputeUp(et,sc,cpGam,cup);
	return getPosteriorOfRates(et,sc,sp,cup,posteriorLike,weights);
}

// fill the posteriorLike matrix with each position posterior rate (p(r|D))
// but without the weights.
MDOUBLE likelihoodComputation::getPosteriorOfRates(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						VVdoubleRep& posteriorLike, 
						const Vdouble * weights) {
	posteriorLike.clear();
	posteriorLike.resize(sc.seqLen());
	for (int z=0; z < posteriorLike.size(); ++z) posteriorLike[z].resize(sp.categories());
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (int pos = 0; pos < sc.seqLen(); ++pos) {
		doubleRep posProb=0;
		for (int categor = 0; categor < sp.categories(); ++categor) {
			doubleRep veryTmp =0;
			for (int let =0; let < sc.getAlphabet()->size(); ++let) {
				veryTmp+=cup.get(pos,categor,et.getRoot()->id(),let) * sp.freq(let);
			}
			posProb += veryTmp*sp.ratesProb(categor);
			posteriorLike[pos][categor] += veryTmp*sp.ratesProb(categor);
		}
		like += log(posProb) * (weights?(*weights)[pos]:1);
		for (int categor1 = 0; categor1 < sp.categories(); ++categor1) {
			posteriorLike[pos][categor1] /= posProb;
		}
	}

	return like;
}


// fill the posteriorLike matrix with each position posterior rate (p(r|D))
// and the LLPP,  but without the weights.
MDOUBLE likelihoodComputation::getPosteriorOfRatesAndLLPP(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGam& cup,
						VVdoubleRep& posteriorLike, 
						VdoubleRep& LLPerPos, 
						const Vdouble * weights) {
	posteriorLike.clear();
	posteriorLike.resize(sc.seqLen());
	for (int z=0; z < posteriorLike.size(); ++z) posteriorLike[z].resize(sp.categories());
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (int pos = 0; pos < sc.seqLen(); ++pos) {
	  LLPerPos[pos] = 0.0;
		for (int categor = 0; categor < sp.categories(); ++categor) {
			doubleRep veryTmp =0;
			for (int let =0; let < sc.getAlphabet()->size(); ++let) {
				veryTmp+=cup.get(pos,categor,et.getRoot()->id(),let) * sp.freq(let);
			}
			LLPerPos[pos] += veryTmp*sp.ratesProb(categor);
			posteriorLike[pos][categor] += veryTmp*sp.ratesProb(categor);
		}
		like += log(LLPerPos[pos]) * (weights?(*weights)[pos]:1);
		for (int categor1 = 0; categor1 < sp.categories(); ++categor1) {
			posteriorLike[pos][categor1] /= LLPerPos[pos];
		}
	}

	return like;
}

MDOUBLE likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(const tree& et,
																  const sequenceContainer& sc,
																  const stochasticProcess& sp,
																  const Vdouble * const weights){
	computePijGam pi;
	pi.fillPij(et,sp);
	MDOUBLE res =0;
	doubleRep LofPos;
	int k;
	for (k=0; k < sc.seqLen(); ++k) {
	    LofPos = likelihoodComputation::getLofPos(k,//pos,
												  et,//const tree& 
												  sc,// sequenceContainer& sc,
												  pi,//const computePijGam& ,
												  sp);
		res += log(LofPos) * (weights?(*weights)[k]:1);;//const stochasticProcess& );
	}
	return res;
}


// MDOUBLE likelihoodComputation::getTreeLikelihoodFromPosteriorAndAlpha(const MDOUBLE alpha,
// 																	  const Vdouble originalBounderi,
// 																	  const VVdouble& posteriorLike,
// 																	  const VdoubleRep& LLPP,
// 																	  const Vdouble* weights)
// {
//   int nCategories = originalBounderi.size()-1;
//   Vdouble rateWeights; rateWeights.resize(nCategories);
//   for (int i=0; i<n; ++i) 
// 	rateWeights[i]=(gammp(alpha, originalBounderi[i+1]*alpha)-gammp(alpha, originalBounderi[i]*alpha))*nCategories;
  
// }

// this function forces non gamma computation of likelihoods from up.
// i.e., even if the stochastic process is really gamma - the likelihood is computed as if there's no gamma.
MDOUBLE likelihoodComputation::getTreeLikelihoodFromUpSpecifcRates(const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalHom& cup,
						VdoubleRep& posLike, // fill this vector with each position likelihood but without the weights.
						const Vdouble * weights)
{
	posLike.clear();
	MDOUBLE like = 0;
	//computing the likelihood from up:
	for (int pos = 0; pos < sc.seqLen(); ++pos) 
	{
		doubleRep tmp=0;
		for (int let =0; let < sc.getAlphabet()->size(); ++let) {
			tmp += cup.get(pos, et.getRoot()->id(), let) * sp.freq(let);
		}
		
		assert(tmp > 0);
		like += log(tmp) * (weights?(*weights)[pos]:1);
		posLike.push_back(tmp);
	}
	return like;	
}

doubleRep likelihoodComputation::getProbOfPosWhenUpIsFilledGam(const int pos,
						const tree& et,
						const sequenceContainer& sc,
						const stochasticProcess& sp,
						const suffStatGlobalGamPos& cup) {
	doubleRep tmp=0;
	for (int categor = 0; categor < sp.categories(); ++categor) {
		doubleRep veryTmp =0;
		for (int let =0; let < sc.alphabetSize(); ++let) {
			veryTmp+=cup.get(categor,et.getRoot()->id(),let) * sp.freq(let);
		}
		tmp += veryTmp*sp.ratesProb(categor);
	}
	assert(tmp>0.0);
	return tmp;
}

MDOUBLE likelihoodComputation::computeLikelihoodAndLikelihoodPerPosition(const sequenceContainer &sc, const tree &et, 
												  const stochasticProcess &sp, Vdouble &LLPerPos) {
	MDOUBLE treeLogLikelihood = 0.0;
	computePijGam cpij;
	cpij.fillPij(et, sp);
	LLPerPos.resize(sc.seqLen());
	doubleRep LofPos;
	for (int pos=0; pos < sc.seqLen() ;++pos) {
		LofPos = likelihoodComputation::getLofPos(pos, et, sc, cpij, sp);
		MDOUBLE tmpLL = log(LofPos);
		treeLogLikelihood += tmpLL;
		LLPerPos[pos] = tmpLL;
	}
	return treeLogLikelihood;
}
