// $Id: siteSpecificRate.cpp 700 2006-05-27 16:57:51Z ninio $

#include "siteSpecificRate.h"
#include "numRec.h"
#include "checkcovFanctors.h"
#include "doubleRep.h"


MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & likelihoodsV,
								   const sequenceContainer& sc,
								   const stochasticProcess& sp,
								   const tree& et,
								   const MDOUBLE maxRate,//20.0f
								   const MDOUBLE tol){//=0.0001f;

	ratesV.resize(sc.seqLen());
	likelihoodsV.resize(sc.seqLen());
	MDOUBLE Lsum = 0.0;

	for (int pos=0; pos < sc.seqLen(); ++pos) {
		computeML_siteSpecificRate(pos,sc,sp,et,ratesV[pos],likelihoodsV[pos],maxRate,tol);
		assert(likelihoodsV[pos]>0.0);
		Lsum += log(likelihoodsV[pos]);
		LOG(5,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
	}
	LOG(5,<<" number of sites: "<<sc.seqLen()<<endl);
	return Lsum;
}

MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& spAttributesVec,
						const Vint& treeAttributesVec,
						const vector<tree> & etVec,
						const vector<const stochasticProcess *> & spVec,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol){
	MDOUBLE Lsum = 0.0;
	ratesV.resize(sc.seqLen()); // the rates themselves
	likelihoodsV.resize(sc.seqLen()); // the log likelihood of each position
	
	for (int pos=0; pos < sc.seqLen(); ++pos) {
		LOG(5,<<".");
		MDOUBLE bestR=-1.0; // tree1
		//		MDOUBLE LmaxR1=0;
		
		// getting the right tree for the specific position:
		const tree*  treeForThisPosition=NULL;
		if ((etVec.size() >0 ) && (treeAttributesVec[pos]>0)) {
			treeForThisPosition = & etVec[ treeAttributesVec[pos] -1];
		} else {
			errorMsg::reportError("tree vector is empty, or treeAttribute is empty, or treeAttribute[pos] is zero (it should be one)");
		}

		// getting the right stochastic process for the specific position:

		const stochasticProcess* spForThisPosition=NULL;

		if ((spVec.size() >0 ) && (spAttributesVec[pos]>0)) {
			spForThisPosition = spVec[ spAttributesVec[pos] -1];
		} else {
			errorMsg::reportError("stochastic process vector is empty, or spAttributesVec is empty, or spAttribute[pos] is zero (it should be one)");
		}

		computeML_siteSpecificRate(pos,sc,*spForThisPosition,*treeForThisPosition,bestR,likelihoodsV[pos],maxRate,tol);
		ratesV[pos] = bestR;
		assert(likelihoodsV[pos]>0.0);
		Lsum += log(likelihoodsV[pos]);
		LOG(5,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
	}
	LOG(5,<<" number of sites: "<<sc.seqLen()<<endl);
	return Lsum;
}

// note that this places the likelihood, rather then the *log*likelihood into posL
void computeML_siteSpecificRate(int pos,
								 const sequenceContainer& sc,
								 const stochasticProcess& sp,
								 const tree &et,
								 MDOUBLE& bestRate,
								 MDOUBLE& posL,
								 const MDOUBLE maxRate,
								 const MDOUBLE tol) {
	LOG(5,<<".");
	MDOUBLE ax=0.00001f,bx=maxRate*0.25,cx=maxRate;	// MN
	posL=-brent(ax,bx,cx,Cevaluate_L_given_r(sc,et,sp,pos),tol,&bestRate);
}

MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& treeAttributesVec,
						const vector<tree> & etVec,
						const stochasticProcess& sp,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol) {
	Vint spAttributesVec(sc.seqLen(),1);
	vector<const stochasticProcess* >  spVec;
	spVec.push_back(&sp);
	return computeML_siteSpecificRate(ratesV,likelihoodsV,
		spAttributesVec,treeAttributesVec,etVec,spVec,sc,maxRate,tol);
}

MDOUBLE computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& spAttributesVec,
						const tree & et,
						const vector<const stochasticProcess* > & spVec,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol){
	Vint treeAttributesVec(sc.seqLen(),1);				
	vector<tree>  etVec;
	etVec.push_back(et);
	return computeML_siteSpecificRate(ratesV,likelihoodsV,
		spAttributesVec,treeAttributesVec,etVec,spVec,sc,maxRate,tol);
}

// THE BAYESIAN EB_EXP PART OF RATE ESTIMATION. //

void computeEB_EXP_siteSpecificRate(int pos,
								 const sequenceContainer& sc,
								 const stochasticProcess& sp,
								 const computePijGam& cpg,
								 const tree &et,
								 MDOUBLE& bestRate,
								 MDOUBLE & stdRate,
								 MDOUBLE & lowerConf,
								 MDOUBLE & upperConf,
								 const MDOUBLE alphaConf) {// alpha of 0.05 is considered 0.025 for each side.
	// here we compute P(r | data)
	VdoubleRep pGivenR(sp.categories(),0.0);
	doubleRep sum=0;
	for (int i=0; i < sp.categories(); ++i) {
		pGivenR[i] = likelihoodComputation::getLofPos(pos,et,sc,cpg[i],sp)*sp.ratesProb(i);
		sum+=pGivenR[i];
	}
	assert(sum!=0);
	
	// here we compute sigma r * P(r | data)
	bestRate = 0.0;
	MDOUBLE Ex2 = 0.0; // this is the sum of squares.
	for (int j=0; j < sp.categories(); ++j) {
		pGivenR[j]/=sum; // So that pGivenR is probability.
		                 // From here on we can convert it back
		                 // to MDOUBLE because it's not a very
		                 // small likelihood any more
		MDOUBLE tmp = convert(pGivenR[j])*sp.rates(j); 
		bestRate += tmp;
		Ex2 += (tmp*sp.rates(j));
	}

	MDOUBLE tmp2 = Ex2 - bestRate*bestRate;
	assert(tmp2>=0);
	stdRate = sqrt(tmp2);

	// detecting the confidence intervals.
	MDOUBLE oneSideConfAlpha = alphaConf/2.0; // because we are computing the two tail.
	doubleRep cdf = 0.0; // cumulative density function.
	int k=0;
	while (k < sp.categories()){
		cdf += convert(pGivenR[k]);
		if (cdf >oneSideConfAlpha) {
			lowerConf = sp.rates(k);
			break;
		} 
		k++;
	}
	while (k < sp.categories()) {
		if (cdf >(1.0-oneSideConfAlpha)) {
			upperConf = sp.rates(k);
			break;
		}
		++k;
		cdf += convert(pGivenR[k]);
	}
	if (k==sp.categories()) upperConf = sp.rates(k-1);
}

void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const sequenceContainer& sc,
								   const stochasticProcess& sp,
								   const tree& et,
								   const MDOUBLE alphaConf){
	ratesV.resize(sc.seqLen());
	stdV.resize(sc.seqLen());
	lowerBoundV.resize(sc.seqLen());
	upperBoundV.resize(sc.seqLen());

	computePijGam cpg;
	cpg.fillPij(et,sp);
	for (int pos=0; pos < sc.seqLen(); ++pos) {
		computeEB_EXP_siteSpecificRate(pos,sc,sp,cpg, et,ratesV[pos],stdV[pos],lowerBoundV[pos],upperBoundV[pos],alphaConf);
		LOG(5,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
	}
	LOG(5,<<" number of sites: "<<sc.seqLen()<<endl);
}

void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& spAttributesVec,
								   const Vint& treeAttributesVec,
							       const sequenceContainer& sc,
								   const vector<tree> & etVec,
								   const vector<const stochasticProcess *> & spVec,
								   const MDOUBLE alphaConf){
	ratesV.resize(sc.seqLen());
	stdV.resize(sc.seqLen());
	lowerBoundV.resize(sc.seqLen());
	upperBoundV.resize(sc.seqLen());
	for (int treeNum=0; treeNum<etVec.size(); ++treeNum) {
		for (int spNum = 0; spNum<spVec.size(); ++spNum) {
            computePijGam cpg;
	    cpg.fillPij(etVec[treeNum],*(spVec[spNum]));
			for (int pos=0; pos < sc.seqLen(); ++pos) {
				if (((spAttributesVec[pos]-1)!=spNum ) || ((treeAttributesVec[pos]-1)!=treeNum )) continue;
				const tree*  treeForThisPosition=NULL;
				assert ((etVec.size() >0 ) && (treeAttributesVec[pos]>0));
				treeForThisPosition = & etVec[ treeAttributesVec[pos] -1];
				const stochasticProcess* spForThisPosition=NULL;
				assert ((spVec.size() >0 ) && (spAttributesVec[pos]>0));
				spForThisPosition = spVec[ spAttributesVec[pos] -1];
				computeEB_EXP_siteSpecificRate(pos,sc,*spForThisPosition,cpg, *treeForThisPosition,ratesV[pos],stdV[pos],lowerBoundV[pos],upperBoundV[pos],alphaConf);
				LOG(5,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
			}
		}
	}
	LOG(5,<<" number of sites: "<<sc.seqLen()<<endl);
}

// one tree many sps
void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& spAttributesVec,
							       const sequenceContainer& sc,
								   const tree & et,
								   const vector<const stochasticProcess *> & spVec,
								   const MDOUBLE alphaConf){
	Vint etAttributesVec(sc.seqLen(),1);				
	vector<tree>  etVec;
	etVec.push_back(et);
	computeEB_EXP_siteSpecificRate(ratesV,stdV,lowerBoundV,upperBoundV,spAttributesVec,etAttributesVec,sc,etVec,spVec,alphaConf);
}

// one sp many trees

void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& treeAttributesVec,
							       const sequenceContainer& sc,
								   const vector<tree> & etVec,
								   const stochasticProcess & sp,
								   const MDOUBLE alphaConf){
	Vint spAttributesVec(sc.seqLen(),1);
	vector<const stochasticProcess* >  spVec;
	spVec.push_back(&sp);
	computeEB_EXP_siteSpecificRate(ratesV,stdV,lowerBoundV,upperBoundV,spAttributesVec,treeAttributesVec,sc,etVec,spVec,alphaConf);
}

