// $Id: goldmanYangModel.h 788 2006-07-19 15:20:58Z eyalprivman $

#ifndef ___GOLDMAN_YANG_MODEL
#define ___GOLDMAN_YANG_MODEL

#include "definitions.h"
#include "replacementModel.h"
#include "fromQtoPt.h"
#include "granthamChemicalDistances.h"
#include "codon.h"

class goldmanYangModel : public replacementModel {
public:
	const int alphabetSize() const {return _codonAlph.size();}
	virtual replacementModel* clone() const { return new goldmanYangModel(*this); }
	explicit goldmanYangModel(const MDOUBLE inV, const MDOUBLE inK,codon & inCodonAlph, const bool globalV=true);
	explicit goldmanYangModel(const MDOUBLE inV, const MDOUBLE inK,codon & inCodonAlph, const Vdouble& freq,const bool globalV=true);
	const MDOUBLE Pij_t(const int i,const int j, const MDOUBLE d) const {
		return _q2pt.Pij_t(i,j,d);
	}
	const MDOUBLE dPij_dt(const int i,const int j, const MDOUBLE d) const{
		return _q2pt.dPij_dt(i,j,d);
	}
	const MDOUBLE d2Pij_dt2(const int i,const int j, const MDOUBLE d) const{
		return _q2pt.d2Pij_dt2(i,j,d);
	}
	const MDOUBLE freq(const int i) const {return _freq[i];};
	void setK(const MDOUBLE newK) { _k = newK;updateQ();}
	void setV(const MDOUBLE newV) { _v = newV;updateQ();}
	void homogenousFreq(){ _freq.erase(_freq.begin(),_freq.end()),_freq.resize(_codonAlph.size(),1.0/_codonAlph.size());}

	MDOUBLE getK() {return _k;}
	MDOUBLE getV() {return _v;}

	void setGlobalV(const bool globalV){ _globalV=globalV;}
	const granthamChemicalDistances& getGCD(){return _gcd;}
	MDOUBLE getQij(const int i,const int j)const {return _Q[i][j];}
	
	VVdouble getQ() const { return _Q;}
	Vdouble getFreqs() const {return _freq;}

private:
	Vdouble _freq;
	MDOUBLE _v; //selection factor.
	MDOUBLE _k; // Tr/Tv ratio.
	void updateQ();
	q2pt _q2pt;
	granthamChemicalDistances _gcd;
	bool _globalV;   //false when compute v per site
	VVdouble _Q;
	codon & _codonAlph;

};


#endif
