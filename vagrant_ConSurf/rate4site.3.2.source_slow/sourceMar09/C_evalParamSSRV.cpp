// 	$Id: C_evalParamSSRV.cpp 920 2006-09-21 09:26:12Z ninio $	
#include "C_evalParamSSRV.h"


MDOUBLE C_evalParamSSRV::operator() (MDOUBLE param) {
		
		setParam(param);
		MDOUBLE res = likelihoodComputation2USSRV::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_baseSc,*_pModel,_weights);
		print(param,res);
		return -res;
}

void C_evalAlphaSSRV::setParam(MDOUBLE alpha)
{
	if (_pModel->noOfCategor() == 1)
		errorMsg::reportError(" one category when trying to optimize alpha");
	_pModel->updateAlpha(alpha);	
}

void C_evalAlphaSSRV::print(MDOUBLE alpha,MDOUBLE res) {
			LOG(5,<<" with Alpha = "<<alpha<<" logL = " <<res<<endl);
}


void C_evalNuSSRV::setParam(MDOUBLE Nu)
{
	_pModel->updateNu(Nu);
}

void C_evalNuSSRV::print(MDOUBLE nu,MDOUBLE res) {
	LOG(5,<<" with Nu = "<<nu<<" logL = " <<res<<endl);
}

void C_evalFSSRV::setParam(MDOUBLE f)
{
	_pModel->updateF(f);
}

void C_evalFSSRV::print(MDOUBLE f,MDOUBLE res) {
	LOG(5,<<" with F = "<<f<<" logL = " <<res<<endl);
}
