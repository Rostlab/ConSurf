// 	$Id: C_evalParamSSRV.h 920 2006-09-21 09:26:12Z ninio $	
#ifndef ___C_EVAL_PARAM_SSRV
#define ___C_EVAL_PARAM_SSRV

#include "definitions.h"

#include "likelihoodComputation.h"
#include "likelihoodComputation2USSRV.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "replacementModelSSRV.h"
#include "stochasticProcessSSRV.h"
#include "ussrvModel.h"
#include "logFile.h"


class C_evalParamSSRV {
public:
	C_evalParamSSRV(const tree& et,
				const sequenceContainer& sc,
				const sequenceContainer& baseSc,
				ussrvModel* pModel,
				const Vdouble* weights = NULL)
    : _et(et),_sc(sc),_baseSc(baseSc),_pModel(pModel),_weights(weights){}

	MDOUBLE operator() (MDOUBLE param) ;
    virtual ~C_evalParamSSRV(){}

protected:
	const tree& _et;
	const sequenceContainer& _sc;
	const sequenceContainer& _baseSc;
	ussrvModel* _pModel;
	const Vdouble * _weights;


protected:
	virtual void setParam(MDOUBLE param) = 0;
	virtual void print(MDOUBLE param,MDOUBLE res) =0;
};


class C_evalAlphaSSRV : public C_evalParamSSRV {
public:
  C_evalAlphaSSRV(const tree& et,
				const sequenceContainer& sc,
				const sequenceContainer& baseSc,
				ussrvModel* pModel,
				const Vdouble *weights = NULL) 
				 : C_evalParamSSRV(et,sc,baseSc,pModel,weights) 
				{}
  
protected:
	virtual void setParam(MDOUBLE alpha);
	virtual void print(MDOUBLE alpha,MDOUBLE res);
};



class C_evalNuSSRV : public C_evalParamSSRV{
public:
  C_evalNuSSRV(	const tree& et,
				const sequenceContainer& sc,
				const sequenceContainer& baseSc,
				ussrvModel* pModel,
				const Vdouble * weights = NULL)
    : C_evalParamSSRV(et,sc,baseSc,pModel,weights){}

protected:
	virtual void setParam(MDOUBLE Nu);
	virtual void print(MDOUBLE nu,MDOUBLE res);
};

class C_evalFSSRV : public C_evalParamSSRV{
public:
  C_evalFSSRV(	const tree& et,
				const sequenceContainer& sc,
				const sequenceContainer& baseSc,
				ussrvModel* pModel,
				const Vdouble * weights = NULL)
    : C_evalParamSSRV(et,sc,baseSc,pModel,weights){}

protected:
	virtual void setParam(MDOUBLE F);
	virtual void print(MDOUBLE f,MDOUBLE res);
};


#endif
