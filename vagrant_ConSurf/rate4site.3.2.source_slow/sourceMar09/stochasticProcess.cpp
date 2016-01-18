// $Id: stochasticProcess.cpp 391 2005-06-12 12:44:57Z ninio $

#include "stochasticProcess.h"
#include "errorMsg.h"

stochasticProcess& stochasticProcess::operator=(const stochasticProcess &otherStoc) {
	if (this != &otherStoc) {              // Check for self-assignment
		pijAccelerator* p2 = otherStoc._pijAccelerator->clone();   // Create the new one FIRST...
		delete _pijAccelerator;                   // ...THEN delete the old one
		_pijAccelerator = p2;

		distribution* d2 =  otherStoc._distr->clone();
		delete _distr;
		_distr = d2;

	}
//	if (_distr) delete _distr;
//	_distr = new distribution(*otherStoc._distr);
    return *this;
}
   
	
stochasticProcess::stochasticProcess(const distribution *in_distr,const pijAccelerator *pijAccelerator) :
	 _distr(in_distr->clone()), _pijAccelerator(pijAccelerator->clone()){
	
}
	
stochasticProcess::stochasticProcess(const stochasticProcess& other):
	 _distr(NULL), _pijAccelerator(NULL) {
		if (other._pijAccelerator != NULL) _pijAccelerator = other._pijAccelerator->clone();
		if (other._distr != NULL) _distr = other._distr->clone();
}
	
stochasticProcess::~stochasticProcess() {
		delete _distr;
		delete _pijAccelerator;
}


void stochasticProcess::setDistribution(const distribution* in_distr)
{
	if (_distr)	delete _distr;
	if (in_distr == NULL) _distr = NULL;
	else _distr = in_distr->clone();
}
