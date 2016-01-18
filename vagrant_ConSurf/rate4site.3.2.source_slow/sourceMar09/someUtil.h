// $Id: someUtil.h 913 2006-09-20 08:42:15Z osnatzomer $

#ifndef ___SOME_UTIL_H
#define ___SOME_UTIL_H

#include "logFile.h"
#include "definitions.h"
#include <string>
using namespace std;

// STATISTICAL UTILITIES:

MDOUBLE computeAverage(const vector<int>& vec);
MDOUBLE computeAverage(const vector<MDOUBLE>& vec);
MDOUBLE computeStd(const vector<MDOUBLE>& vec);// page 60, Sokal and Rohlf
MDOUBLE computeStd(const vector<int>& vec);// page 60, Sokal and Rohlf
MDOUBLE copmutePoissonProbability(const int& k, const long double& lamda);

// TIME UTILITIES
void printTime(ostream& out);

// TEXT UTILITIES
string int2string(const int i);
string double2string(const double x, int const howManyDigitsAfterTheDot=5);
MDOUBLE string2double(const string& inString);
bool allowCharSet(const string& allowableChars, const string& string2check);
bool isCharInString(const string& stringToCheck, const char charToCheck);
void putFileIntoVectorStringArray(istream &infile,vector<string> &inseqFile);

bool fromStringIterToInt(string::const_iterator & it,
						 const string::const_iterator endOfString,
						 int& res);

string takeCharOutOfString(const string& charsToTakeOut, const string& fromString);
void toLower(string& str);
void toUpper(string& str);
// FILE UTILITIES
bool checkThatFileExist(const string& fileName); 
string* searchStringInFile(const string& string2find,
						   const int index,
						   const string& inFileName);
string* searchStringInFile(const string& string2find,
						   const string& inFileName);
bool doesWordExistInFile(const string& string2find,const string& inFileName);
void createDir(const string& curDir,const string& dirName);


//BIT UTILITIES
//void nextBit(bitset<64> &cur);

//ARITHMETIC UTILITIES
//DEQUAL: == UP TO EPSILON
//DBIG_EQUAL: >= UP TO EPSILON
//DSMALL_EQUAL: <= UP TO EPSILON
bool DEQUAL(const MDOUBLE x1, const MDOUBLE x2, const MDOUBLE epsilon = 1.192092896e-07F); // epsilon taken from WINDOW'S FILE FLOAT.H
bool DBIG_EQUAL(const MDOUBLE x1, const MDOUBLE x2, const MDOUBLE epsilon = 1.192092896e-07F); 
bool DSMALL_EQUAL(const MDOUBLE x1, const MDOUBLE x2, const MDOUBLE epsilon = 1.192092896e-07F); // {return ((x1 < x2) || DEQUAL(x1, x2));}

//swap between the 4 variables such that the first becomes the second, second becomes the third and third becomes the fourth.
//used in functoin mnbrack below.
void shift3(MDOUBLE &a, MDOUBLE &b, MDOUBLE &c, const MDOUBLE d);


// print vector and VVdoulbe util
ostream &operator<<(ostream &out, const Vdouble &v);
ostream &operator<<(ostream &out, const VVdouble &m);
void mult(Vdouble& vec, const MDOUBLE factor);
void mult(VVdouble& vec, const MDOUBLE factor);
//scale vecToScale so that its new average is AvgIn. return the scaling factor. 
MDOUBLE scaleVec(Vdouble& vecToScale, const MDOUBLE avgIn);
//determine the relative order of vecIn. The order vector is returned 
//ex: vecIn = [0.1 0.4 0.01 0.9 1.8] orderVecOut = [1 2 0 3 4] 
Vdouble orderVec(const Vdouble& vecIn); // double because sometimes ranks are 1.5 (in case of ties, for example).
MDOUBLE calcRelativeMSEDistBetweenVectors(const Vdouble& trueValues, const Vdouble& inferredValues, const MDOUBLE threshhold = 0.0);
MDOUBLE calcMSEDistBetweenVectors(const Vdouble& trueValues, const Vdouble& inferredValues);
//MAD = mean absolute deviations distance 
MDOUBLE calcMADDistBetweenVectors(const Vdouble& oneRatesVec, const Vdouble& otherRatesVec);
MDOUBLE calcRelativeMADDistBetweenVectors(const Vdouble& trueValues, const Vdouble& inferredValues, const MDOUBLE threshhold = 0.0);

//to be used for orderVec
template <class T>
class vecElem
{
public:
	vecElem();
	virtual ~vecElem() {};
	void setValue(const T val) {m_value = val;}
	T getValue() {return m_value;}
	void setPlace(const int place) {m_place = place;}
	int getPlace() {return m_place;}
	inline bool operator< (const vecElem& elemIn) const;
private:
	int m_place;
	T m_value;
};


template <class T>
vecElem< T >::vecElem()
{
	m_value = -1;
	m_place = -1;
}

//template <class T>
//vecElement< T >::~vecElement()
//{
//}
template <class T>
bool vecElem< T >::operator<(const vecElem& elemIn) const
{
	if (m_value == elemIn.m_value)
		return (m_place  < elemIn.m_place);
	else
		return (m_value < elemIn.m_value);
}
#endif

