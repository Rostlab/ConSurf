// $Id: codon.cpp 827 2006-07-24 16:51:24Z ninio $

#include "codon.h"
#include "nucleotide.h"
#include "amino.h"
#include "logFile.h"
#include "definitions.h"
#include "someUtil.h"
#include <sstream>
#include <cctype>

codon::codon(){
	geneticCodeString gcs=geneticCodeHolder::nuclearStandard;
	readMatrixFromFile(gcs.Val);
}

void codon::readMatrixFromFile(const string& matrixFileName){ //default value: "nuclearCode.txt"
  //	cout<<"in codon constructor"<<endl;
	stringstream in(matrixFileName.c_str());
	if (!in) {
		errorMsg::reportError("in codon::readMatrixFromFile: unable to open matrix data file");
	}

	int aa = -1; //initialized as -1 so in first iteration will change to 0
	int noOfCodons = 0;
	string strAmino; 
	while (!in.eof()) { //20 amino acids and stop 
		string val;
		in>>val;
		if (val.size()==1) { //amino acid
			aa++;
			strAmino=val;
			if (strAmino=="*") { _alphabetSize=noOfCodons;}
		}
		else if (val.size()==3 && val[0]!='#'){ //codon, # symbolizes a comment
			_geneticCode[val]=strAmino;
			_codon2Int[val]=noOfCodons;
			noOfCodons++;
		}
		else {
			
			if (noOfCodons!=64){
				string err="in codon::readMatrixFromFile: total number of codons = "+int2string(noOfCodons);
				errorMsg::reportError(err);
			}
			return;
		}
	}
}
codon& codon::operator=(const codon& other) {
	_geneticCode = other._geneticCode; //key - codon, value - amino acid
	_codon2Int = other._codon2Int;//key string of codon int= integer value of codon
	_alphabetSize = other._alphabetSize;


	return *this;
}
// codon::codon(const  codon& other):
// 	_geneticCode(other._geneticCode), //key - codon, value - amino acid
// 	_codon2Int(other._codon2Int),//key string of codon int= integer value of codon
// 	_alphabetSize(other._alphabetSize){}


//return -99 if not succeeds.
int codon::fromChar(const string& s, const int pos) const {
	if (s.size() <= pos+2) {
		//errorMsg::reportError("Trying to read a codon pass the end of the string. The number of nucleotide may not be divisible by three");
		string textToPrint("Trying to read a codon pass the end of the string. The number of nucleotide may not be divisible by three");
		LOG(1,<<textToPrint<<endl);
		return -99;
	}

	nucleotide nuc;
	int p1,p2,p3;
	p1 = nuc.fromChar(s[pos]);
	p2 = nuc.fromChar(s[pos+1]);
	p3 = nuc.fromChar(s[pos+2]);


	if ((p1 <0) || (p2 <0) || (p3 <0)) 
		return gap(); 
	else if ((p1 ==15) || (p2 ==15) || (p3 ==15)) return unknown(); // unknown.
	else if ((p1 >4) || (p2 >4) || (p3 >4)) return unknown(); //unknown.
	string strCodon="";
	//change U --> T
	if (p1==4) strCodon+="T";
	else  strCodon+=toupper(s[pos]);
	if (p2==4) strCodon+="T";
	else  strCodon+=toupper(s[pos+1]);
	if (p3==4) strCodon+="T";
	else  strCodon+=toupper(s[pos+2]);


	//const string strCodon = s.substr(pos,3);
	map <string,int> tmpMap=_codon2Int;
	map <string,int>::iterator it1;
	it1=tmpMap.find(strCodon);
	if (it1==tmpMap.end()){
		
		string err="error in codon::fromChar cannot find codon "+strCodon;
		errorMsg::reportError(err);
	}
	return tmpMap[strCodon];

	
}

vector<int> codon::fromString(const string &str) const {
	vector<int> vec;
	if (str.size()%3!=3) {
		errorMsg::reportError("error in function codon::fromString. String length should be a multiplication of 3");
	}
	for (int i=0;i<str.size();i+=3)
	  vec.push_back(fromChar(str,i));
	return vec;
}

string codon::fromInt(const int in_id)  const{
	if (in_id == unknown() || in_id == gap()) return "---";
	map <string, int> tmpMap = _codon2Int;
	map <string, int>::iterator it=tmpMap.begin();
	while (it!=tmpMap.end()){
		if ((*it).second==in_id){
			return (*it).first;
		}
		it++;
	}
	string err="error in function codon::fromInt: no codon found for the integer";
	errorMsg::reportError(err);
	return (string("we should never get here - the reportError above will exit"));
	
}

codonUtility::replacementType codonUtility::codonReplacement(const int c1, const int c2, codon &cod){
	if (c1 == c2) return codonUtility::sameCodon;
	else if (codonUtility::aaOf(c1,cod) == codonUtility::aaOf(c2,cod)) return codonUtility::synonymous;
	return codonUtility::non_synonymous;
}

int codonUtility::aaOf (const int c1, const codon &cod){
    amino a;
    if (c1==cod.gap()||c1==cod.unknown()) 
		return a.gap();
	string strCodon=cod.fromInt(c1);
	map <string,string> geneticCode=cod.geneticCode();
	map <string,string>::iterator pos;
	pos=geneticCode.find(strCodon);
	if (pos==geneticCode.end()){
		string err="error in codonUtility::aaOf: cannot find codon "+strCodon;
		errorMsg::reportError(err);
	}
	string strAmino=geneticCode[strCodon];
	if (strAmino.size()>1){
		errorMsg::reportError("error in codonUtility::aaOf: amino acid 1 letter code > 1");
	}
	return a.fromChar(*strAmino.c_str());

}





codonUtility::diffType codonUtility::codonDiff(const int c1, const int c2, codon &cod){
	if (c1==c2) return codonUtility::equal;
	nucleotide n;
	string s1 = cod.fromInt(c1);
	string s2 = cod.fromInt(c2);

	int pos1 = n.fromChar(s1[0])+n.fromChar(s2[0]);
	int pos2 = n.fromChar(s1[1])+n.fromChar(s2[1]);
	int pos3 = n.fromChar(s1[2])+n.fromChar(s2[2]);


	if (s1[0]!=s2[0] && s1[1]!=s2[1] && s1[2]!=s2[2])
		return  codonUtility::threesub; 

	if (s1[0]==s2[0] && s1[1]==s2[1] && s1[2]!=s2[2]) {
		if (pos3%2==0) return codonUtility::tr;
		else return codonUtility::tv;
	
	}
	if (s1[1]==s2[1] && s1[2]==s2[2] && s1[0]!=s2[0]) {
		if (pos1%2==0) return codonUtility::tr;
		else return codonUtility::tv;
	
	}
	if (s1[0]==s2[0] && s1[2]==s2[2] && s1[1]!=s2[1]) {
		if (pos2%2==0) return codonUtility::tr;
		else return codonUtility::tv;
	
	}
	
	


	if (s1[0]==s2[0] && pos2%2==0 && pos3%2==0)
		return  codonUtility::twoTrs;
	if (s1[1]==s2[1] && pos1%2==0 && pos3%2==0)
		return  codonUtility::twoTrs;
	if (s1[2]==s2[2] && pos1%2==0 && pos2%2==0)
		return  codonUtility::twoTrs;


	if (s1[0]==s2[0] && pos2%2!=0 && pos3%2!=0)
		return  codonUtility::twoTvs;
	if (s1[1]==s2[1] && pos1%2!=0 && pos3%2!=0)
		return  codonUtility::twoTvs;
	if (s1[2]==s2[2] && pos1%2!=0 && pos2%2!=0)
		return  codonUtility::twoTvs;

	return codonUtility::trtv;

}
