// 	$Id: replacementModelSSRV.cpp 920 2006-09-21 09:26:12Z ninio $	

#include "replacementModelSSRV.h"
#include "logFile.h"
#include <iomanip>
#include <iostream>


replacementModelSSRV::replacementModelSSRV(const distribution* dist, const replacementModel* baseRM, MDOUBLE rateOfRate /*= 1 */) :
_dist(dist->clone()),
_baseRM(baseRM->clone()),
_rateOfRate(rateOfRate)
{
	if (_rateOfRate < EPSILON) _rateOfRate = EPSILON; // Temporary - to overcome a bug in QL algorithm, when _rateOfRate == 0
	if (_dist->categories() == 0)
		errorMsg::reportError("replacementModelSSRV::replacementModelSSRV : number of categories == 0");
	
	updateFreq();
	updateQ();

	
}

//// similar to goldmanYangModel.cpp 
//replacementModelSSRV::replacementModelSSRV(const replacementModelSSRV& other) :
//_dist(other._dist->clone()),
//_baseRM(other._baseRM->clone()),
//_rateOfRate(other._rateOfRate)
//{
//	int size = alphabetSize();
//	_Q.resize(size);
//	for (int z=0; z < _Q.size();++z) 
//		_Q[z].resize(size,0);
//	updateFreq();
//	updateQ();
//}

// Instead of calling updateQ here, like in goldmanYangModel.cpp,
// this method uses the copy constructor of q2pt and also copies _freq and _Q
replacementModelSSRV::replacementModelSSRV(const replacementModelSSRV& other) :
_dist(other._dist->clone()),
_baseRM(other._baseRM->clone()),	
_rateOfRate(other._rateOfRate),
_q2pt(other._q2pt),
_freq(other._freq),
_Q(other._Q)
{
}

replacementModelSSRV::~replacementModelSSRV()
{
	if (_dist) delete (_dist);
	if (_baseRM) delete (_baseRM);
}


replacementModelSSRV& replacementModelSSRV::operator=(const replacementModelSSRV &other)
{
	if (_dist) delete (_dist);
	if (_baseRM) delete (_baseRM);
	
	_dist = other._dist->clone();
	_baseRM = other._baseRM->clone();	
	_rateOfRate = other._rateOfRate;
	_q2pt = other._q2pt; //@@@@ why doesn't this work ? explicit ?
//	_q2pt.fillFromRateMatrix(other._freq,other._Q);
	_freq = other._freq;
	_Q = other._Q;

	return (*this);
}

const int replacementModelSSRV::alphabetSize() const
{
	return (_baseRM->alphabetSize() * _dist->categories());
}



// The freq of each mulCharacter is its freq in the _baseRM divided by the mulFactor.
// mulFactor = the number of times that the _baseRM is multiplied\duplicated.
void replacementModelSSRV::updateFreq()	
{
	_freq.clear();
	int size = alphabetSize();
	int mulFactor = _dist->categories();
	_freq.resize(size);
	int idInCategory;
		
	for(idInCategory=0; idInCategory < _baseRM->alphabetSize() ; ++idInCategory)
	{
		MDOUBLE freq = _baseRM->freq(idInCategory) / mulFactor ; 
		for (int categoryNumber=0; categoryNumber < mulFactor; ++categoryNumber)
			_freq[categoryNumber*_baseRM->alphabetSize() + idInCategory] = freq;
	}

	// debug OZ
	/*cout<<"baseAlphabet freqs: " << endl;
	for (int i=0;  i < _baseRM->alphabetSize(); ++i)
		cout << _baseRM->freq(i) << '\t';
	cout << endl;
	cout << "replacementModelSSRV freqs:" << endl;
	Vdouble::iterator itr = _freq.begin();
	for (; itr != _freq.end() ; ++itr)
		cout << *itr << '\t';
	cout << endl;*/
	// end of debug
	

}


void replacementModelSSRV::updateQ()
{
	_Q.clear();
	int size = alphabetSize();
	_Q.resize(size);
	for (int z=0; z < _Q.size();++z) 
		_Q[z].resize(size,0.0);
	
	// fill Q
	int _BaseRM_alphabetSize = _baseRM->alphabetSize();
	int mulFactor = _dist->categories();
	MDOUBLE minus = (mulFactor-1) * _rateOfRate / mulFactor ;
	for (int i=0; i < _BaseRM_alphabetSize; ++i)
	{
		for (int j=0; j < _BaseRM_alphabetSize; ++j)
		{
			for (int z=0; z < mulFactor; ++z)
			{
				for (int w=0; w < mulFactor; ++w)
				{
					if (i!=j)
					{
						if (z==w)
							_Q[z*_BaseRM_alphabetSize + i][z*_BaseRM_alphabetSize+j] 
							= _dist->rates(z) * _baseRM->dPij_dt(i,j,0);
					}
					else
					{
						if (z!=w)
						{
							_Q[z*_BaseRM_alphabetSize+i][w*_BaseRM_alphabetSize+i] =  _rateOfRate / mulFactor;
						}
						else
							_Q[z*_BaseRM_alphabetSize+i][z*_BaseRM_alphabetSize+i] = _dist->rates(z) * _baseRM->dPij_dt(i,j,0) -minus;
					}

				}
			}
		}
	}
	
//	// check OZ
//	LOG(4, <<"THE Q MATRIX IS: "<<endl ) ;
//	VVdouble::iterator itr1 = _Q.begin();
//	Vdouble::iterator itr2;
//	for (; itr1 != _Q.end(); ++itr1)
//	{
//		for (itr2 = itr1->begin(); itr2 != itr1->end(); ++itr2)
//			LOG(4,<< setprecision(3) <<  setw(5) << *itr2 <<'\t');
//		LOG(4,<<endl);
//	}
//	LOG (4,<<endl);
////	 end of check

	_q2pt.fillFromRateMatrix(_freq,_Q); 
	
	/* wrong
	for (int x=0; x < _BaseRM_alphabetSize; ++x)
	{
		for (int y=0; y < _BaseRM_alphabetSize; ++y)
		{
			MDOUBLE Mxy = _baseRM->dPij_dt(x,y,0);
			if (x!=y) 
			{
				for (int z=0; z < mulFactor; ++z)
				{
				
					// Q(X,ri --> Y,ri)
					_Q[x+z*_BaseRM_alphabetSize][y+z*_BaseRM_alphabetSize] = 
						_baseRM->dPij_dt(x,y,0) * _dist->rates(z) ;
				}
			}
			else
			{
				for (int z=0; z < mulFactor; ++z)
				{
					
					// Q(X,ri --> X,ri)
					int index = x+z*_BaseRM_alphabetSize ;
					_Q[index][index] = 
						Mxy * _dist->rates(z) - minus ;
				}
				for (int z=1; z < mulFactor; ++z)
				// Q(X,ri --> X,rj)
					_Q[x][y+z*_BaseRM_alphabetSize] = _Q[x+z*_BaseRM_alphabetSize][y] = _rateOfRate / mulFactor ;
			}
		}
	}
	*/ 
}

void replacementModelSSRV::setDistribution(const distribution* dist)
 {
	 if (dist->categories() == 0)
		errorMsg::reportError("replacementModelSSRV::setDistribution : number of categories == 0");
	 if (_dist) delete (_dist);
		_dist=dist->clone();
	updateQ();
 }

MDOUBLE replacementModelSSRV::sumPijQij() const{
	MDOUBLE sum=0.0;
	for (int i=0; i < _Q.size(); ++i) {
		sum -= _Q[i][i]*_freq[i];
	}
	return sum;
}


//void replacementModelSSRV::norm(MDOUBLE scale){
//	
//	for (int i=0; i < _Q.size(); ++i) {
//		for (int j=0; j < _Q.size(); ++j) {
//			_Q[i][j]*=scale;
//		}
//	}
//	
//	_q2pt.fillFromRateMatrix(_freq,_Q);
//}


// WRONG !!!
//const MDOUBLE replacementModelSSRV::Pij_t(const int i,const int j, const MDOUBLE d) const
//{
//	int category_i = i / _baseRM->alphabetSize();
//	int category_j = j / _baseRM->alphabetSize();
//	int idInCategory_i = i % _baseRM->alphabetSize();
//	int idInCategory_j = j % _baseRM->alphabetSize();
//	
//	if (category_i == category_j)
//		return (_baseRM->Pij_t(idInCategory_i,idInCategory_j,d) * _dist->rates(category_i));
//
//	if (idInCategory_i == idInCategory_j)
//		return ((_rateOfRate / _dist->categories())*d);
//
//	return 0;
//}
//
//
//
//const MDOUBLE replacementModelSSRV::dPij_dt(const int i,const int j, const MDOUBLE d) const
//{
//	int category_i = i / _baseRM->alphabetSize();
//	int category_j = j / _baseRM->alphabetSize();
//	int idInCategory_i = i % _baseRM->alphabetSize();
//	int idInCategory_j = j % _baseRM->alphabetSize();
//
//	if (category_i == category_j)
//		return (_baseRM->dPij_dt(idInCategory_i,idInCategory_j,d) * _dist->rates(category_i));
//
//	if (idInCategory_i == idInCategory_j)
//		return (_rateOfRate / _dist->categories());
//
//	return 0;
//}
//
//const MDOUBLE replacementModelSSRV::d2Pij_dt2(const int i,const int j, const MDOUBLE d) const
//{
//	int category_i = i / _baseRM->alphabetSize();
//	int category_j = j / _baseRM->alphabetSize();
//	int idInCategory_i = i % _baseRM->alphabetSize();
//	int idInCategory_j = j % _baseRM->alphabetSize();
//
//	if (category_i == category_j)
//		return (_baseRM->d2Pij_dt2(idInCategory_i,idInCategory_j,d) * _dist->rates(category_i));
//
//	return 0;
//}






