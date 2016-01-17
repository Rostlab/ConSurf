// 	$Id: doubleRep.h 949 2006-10-19 10:18:12Z eyalprivman $	
#ifndef __DOUBLE_REP_H
#define __DOUBLE_REP_H

#include "definitions.h"

#ifdef DOUBLEREP

#include <iostream>
#include <cmath>
using namespace std;

/* doubleRep: enables working with much larger or smaller numbers than normally possible
by the regular double representation
 * Representation of a double x as x=_mantissa*2^_expon
 Note: Base is 2!!
 */ 

class doubleRep{
public:
	
	doubleRep(){};
	explicit doubleRep(MDOUBLE mantissa, int expon);
	doubleRep(MDOUBLE a);
	doubleRep(const doubleRep& other);
	doubleRep* clone() {return new doubleRep(*this);}

	void output(ostream &out) const{ out<<_mantissa<<string(" * 2^")<<_expon;}
    void output0x(ostream &out) const{ double e0x=_expon*0.3010299956639; // log_10(2)
	  int e=(int)(e0x)-1;
	  if (double(e) > e0x - 1.0) e--;
	  double m=_mantissa*pow(10.0,e0x-e); 
	  out<<m;
	  if (e<0)
		out<<"e"<<e;
	  else
		out<<"e+"<<e;
	}
	void outputn(ostream &out) { out<<_mantissa<<string(" * 2^")<<_expon<<endl;}
	
	friend MDOUBLE convert(const doubleRep& a);
	inline doubleRep& operator=(const doubleRep& a);
	inline doubleRep& operator+=(doubleRep a);
	friend inline doubleRep operator+(const doubleRep& a, const doubleRep& b);
	inline doubleRep& operator-=(const doubleRep& a);
	friend inline doubleRep operator-(const doubleRep& a, const doubleRep& b);
	inline doubleRep& operator*=(const doubleRep& a);
	friend inline doubleRep operator*(const doubleRep& a, const doubleRep& b);
	inline doubleRep& operator/=(const doubleRep& a);
	friend inline doubleRep operator/(const doubleRep& a, const doubleRep& b);

	friend inline bool operator==(const doubleRep& a, const doubleRep& b);
	friend inline bool operator!=(const doubleRep& a, const doubleRep& b);
	friend inline bool operator<(const doubleRep& a, const doubleRep& b);
	friend inline bool operator<=(const doubleRep& a, const doubleRep& b);
	friend inline bool operator>(const doubleRep& a, const doubleRep& b);
	friend inline bool operator>=(const doubleRep& a, const doubleRep& b);
	
    const MDOUBLE d_log() const;
//	friend ostream& operator<<(ostream &out, const doubleRep& a);

	const MDOUBLE mantissa() const {return _mantissa;}
	const int expon() const {return _expon;}

private:
	void fixParams(); 


private:
	MDOUBLE _mantissa;
	int _expon;
};
	
typedef vector<doubleRep> VdoubleRep;
typedef vector <vector<doubleRep> > VVdoubleRep;

inline doubleRep& doubleRep::operator=(const doubleRep& a){
	_mantissa=a.mantissa();
	_expon=a.expon();
	return *this;
}

// Original version by Adi Stern
inline doubleRep& doubleRep::operator+=(doubleRep a){
	//ensuring that (*this) is bigger than 'a' for sake of convenience
	if (a.expon()>_expon || ((a.expon()==_expon) && (a.mantissa()>_mantissa))){
		MDOUBLE tmpMant=0.0; int tmpExp=0;
		tmpMant=_mantissa;
		tmpExp=_expon;
		_mantissa=a.mantissa();
		a._mantissa=tmpMant;
		tmpExp=_expon;
		_expon=a.expon();
		a._expon=tmpExp;
	}
	if (a.mantissa()==0)
		return *this;
	if (_mantissa==0){
		_mantissa=a.mantissa();
		_expon=a.expon();
		return *this;
	}
	if (abs(static_cast<MDOUBLE>(_expon-a.expon()))>51.0){ //limit of epsilon difference
		return *this;
	}
	_mantissa+=a.mantissa()*pow(2.0,(a.expon()-_expon)*1.0);
	fixParams();
	return *this;
}

inline doubleRep operator+(const doubleRep& a, const doubleRep& b){
	doubleRep temp(a);
	temp+=b;
	return temp;
}

inline doubleRep& doubleRep::operator-=(const doubleRep& a){
	doubleRep b(-a.mantissa(),a.expon());
	doubleRep me(_mantissa,_expon);
	me+=b;
	_mantissa=me.mantissa();
	_expon=me.expon();
	return *this;
}

inline doubleRep operator-(const doubleRep& a, const doubleRep& b){
	doubleRep temp(a);
	temp-=b;
	return temp;
}

inline doubleRep& doubleRep::operator*=(const doubleRep& a){
	_mantissa*=a.mantissa();
	_expon+=a.expon();
	fixParams();
	return *this;
}

inline doubleRep operator*(const doubleRep& a, const doubleRep& b){
	doubleRep temp(a);
	temp*=b;
	return temp;
}

inline doubleRep& doubleRep::operator/=(const doubleRep& a){
	_mantissa/=a.mantissa();
	_expon-=a.expon();
	fixParams();
	return *this;
}

inline doubleRep operator/(const doubleRep& a, const doubleRep& b){
	doubleRep temp(a);
	temp/=b;
	return temp;
}

/************************
 * Comparison operators *
 ************************/
inline bool operator==(const doubleRep& a, const doubleRep& b){
	return (a._mantissa==b._mantissa && a._expon==b._expon);
}
inline bool operator!=(const doubleRep& a, const doubleRep& b){
	return !(a==b);
}

inline bool operator<(const doubleRep& a, const doubleRep& b){
	// if the numbers have opposite signs
    if (a._mantissa*b._mantissa<0.0){
		if (a._mantissa<b._mantissa) {return true;}
		else {return false;}
    }
	// if the expon values are different
	if (a._expon!=b._expon) {
		// special case where one number is zero
		if (a._mantissa == 0.0) {
			if (b._mantissa > 0.0) {return true;}
			else {return false;}
		}
		if (b._mantissa == 0.0) {
			if (a._mantissa < 0.0) {return true;}
			else {return false;}
		}

		if (a._expon<b._expon) {
			if (a._mantissa > 0.0) {return true;}
			else {return false;}
		} else {
			if (a._mantissa < 0.0) {return true;}
			else {return false;}
		}
		// expon values are identical
	} else {
		return (a._mantissa < b._mantissa);
	}
}

inline bool operator>(const doubleRep& a, const doubleRep& b){
	// if the numbers have opposite signs
    if (a._mantissa*b._mantissa<0.0){
		if (a._mantissa>b._mantissa) {return true;}
		else {return false;}
    }
	// if the expon values are different
	if (a._expon!=b._expon) {
		// special case where one number is zero
		if (a._mantissa == 0.0) {
			if (b._mantissa < 0.0) {return true;}
			else {return false;}
		}
		if (b._mantissa == 0.0) {
			if (a._mantissa > 0.0) {return true;}
			else {return false;}
		}

		if (a._expon>b._expon) {
			if (a._mantissa > 0.0) {return true;}
			else {return false;}
		} else {
			if (a._mantissa < 0.0) {return true;}
			else {return false;}
		}
		// expon values are identical
	} else {
		return (a._mantissa > b._mantissa);
	}
}

inline bool operator<=(const doubleRep& a, const doubleRep& b){
	return !(a>b);
}

inline bool operator>=(const doubleRep& a, const doubleRep& b){
	return !(a<b);
}




ostream& operator<<(ostream &out, const doubleRep& a);

inline MDOUBLE log(const doubleRep& d) {return d.d_log();}

inline ostream &operator<<(ostream &out, const VdoubleRep &v){
  for (int j=0;j<v.size();++j)
    out<< v[j]<<" ";
  out <<endl;
  return(out);
}

inline ostream &operator<<(ostream &out, const VVdoubleRep &m){
  for (int i=0;i<m.size();++i)
    out<<m[i];
  out <<endl;
  return(out);
}


#else

typedef MDOUBLE  doubleRep;
typedef Vdouble  VdoubleRep;
typedef VVdouble VVdoubleRep;



inline MDOUBLE convert (MDOUBLE d) {return(d);}
//inline const MDOUBLE convert (const MDOUBLE d) const  {return(d);}
#endif

#endif
