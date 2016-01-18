// 	$Id: bestAlphaAndNu.cpp 920 2006-09-21 09:26:12Z ninio $	
#include <iostream>
using namespace std;

#include "bestAlphaAndNu.h"

MDOUBLE bestFFixedTreeSSRV::operator()(const tree& et, 
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnF,
					   const MDOUBLE epsilonFOptimization){
	
	MDOUBLE bestF=0;
	const MDOUBLE cx=upperBoundOnF;// left, midle, right limit on alpha
	const MDOUBLE bx=model.getF();
	const MDOUBLE ax=0.0;
	LOG(5,<<"****    Optimizing F    **** " << endl<< "bestFFixedTreeSSRV::operator() bx is :" << bx << endl);
	LOG(9,<<"ax is :" << ax << " cx is :" << cx << endl);
	_bestL = -brent(ax,bx,cx,
		C_evalFSSRV(et,sc,baseSc,&model,weights),
		epsilonFOptimization,
		&bestF);
	setF(bestF,model);
	_bestF= bestF;
	return _bestL;
}

MDOUBLE bestAlphaFixedTreeSSRV::operator()(const tree& et, //findBestAlphaFixedTree
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnAlpha,
					   const MDOUBLE epsilonAlphaOptimization){
	
	MDOUBLE bestA=0;
	const MDOUBLE cx=upperBoundOnAlpha;// left, midle, right limit on alpha
	const MDOUBLE bx=model.getAlpha();
	const MDOUBLE ax=0.0;
	LOG(5,<<"****    Optimizing Alpha    **** " << endl<< "bestAlphaFixedTreeSSRV::operator() bx is :" << bx << endl);
	_bestL = -brent(ax,bx,cx,
		C_evalAlphaSSRV(et,sc,baseSc,&model,weights),
		epsilonAlphaOptimization,
		&bestA);
	setAlpha(bestA,model);
	_bestAlpha= bestA;
	return _bestL;
}

// Alpha is fixed
MDOUBLE bestNuFixedTreeSSRV::operator()(const tree& et, 
					   const sequenceContainer& sc,
					   const sequenceContainer& baseSc,
					   ussrvModel& model,
					   const Vdouble * weights,
					   const MDOUBLE upperBoundOnNu,
					   const MDOUBLE epsilonNuOptimization){
		
	
	MDOUBLE bestN=0;
	// define the Nu bounds
	const MDOUBLE cx=upperBoundOnNu;// left, midle, right limit on alpha
	const MDOUBLE bx= model.getNu(); 
	const MDOUBLE ax=0.0;
	LOG(5,<<"****    Optimizing Nu    **** " << endl << "bestNuFixedTreeSSRV::operator() bx is : " << bx << endl);
	_bestL = -brent(ax,bx,cx, C_evalNuSSRV(et,sc,baseSc,&model,weights), epsilonNuOptimization, &bestN);
	setNu(bestN,model);
	_bestNu= bestN;
	return _bestL;
}

