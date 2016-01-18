// 	$Id: stochasticProcess.h 920 2006-09-21 09:26:12Z ninio $	

#ifndef ___STOCHASTIC_PROCESS
#define ___STOCHASTIC_PROCESS

#include "pijAccelerator.h"
#include "distribution.h"
#include <cassert>

class stochasticProcess{
public:
	explicit stochasticProcess(const distribution *in_distr,const pijAccelerator *pijAccelerator);
	explicit stochasticProcess() {
		_distr=NULL; _pijAccelerator=NULL;
	}
	
	const int alphabetSize() const {return _pijAccelerator->alphabetSize();} // The alphabet size is the same as the matrix Pij size

	stochasticProcess(const stochasticProcess& other);
	virtual const int categories() const {return _distr->categories();}
	virtual const MDOUBLE rates(const int i) const {return _distr->rates(i);}
	virtual const MDOUBLE ratesProb(const int i) const {return _distr->ratesProb(i);}

	
	virtual const MDOUBLE Pij_t(const int i, const int j, const MDOUBLE t) const {
		if (t!=0) return _pijAccelerator->Pij_t(i,j,t);
		return (i==j)? 1 : 0;
	}

	const MDOUBLE freq(const int i) const {assert(i>=0);return _pijAccelerator->freq(i);}	// P(i)
	const MDOUBLE dPij_dt(const int i,const  int j,const MDOUBLE t) const {	return _pijAccelerator->dPij_dt(i,j,t);}
	const MDOUBLE d2Pij_dt2(const int i, const int j, const MDOUBLE t) const { return _pijAccelerator->d2Pij_dt2(i,j,t);}


	virtual distribution* distr() const {return _distr;} // @@@@ this const is a lie !!!
	virtual const pijAccelerator* getPijAccelerator() const {return _pijAccelerator;}
	virtual void setDistribution(const distribution* in_distr);

	stochasticProcess& operator=(const stochasticProcess &otherStoc);
	virtual ~stochasticProcess();
 	virtual void setGlobalRate(const MDOUBLE x) {_distr->setGlobalRate(x);}
 	virtual MDOUBLE getGlobalRate() const {return _distr->getGlobalRate();}

protected:
	distribution *_distr;
	pijAccelerator *_pijAccelerator;
};



#endif


// Stochastic process is composed of two objects: a distribution of rates and a Pij accelerator.
