// $Id: codon.h 935 2006-10-17 08:38:01Z ninio $
#ifndef ____CODON
#define ____CODON

#include <cassert>
#include "definitions.h"
#include "errorMsg.h"
#include "alphabet.h"
#include "geneticCodeHolder.h"
#include <map>
class codon;


class codonUtility {
public:
	enum diffType {equal =0, tr, tv, twoTrs, twoTvs ,trtv, threesub};
	static diffType codonDiff(const int c1, const int c2, codon &cod);

	enum replacementType {sameCodon=0, synonymous, non_synonymous};
	static replacementType codonReplacement(const int c1, const int c2, codon &cod);

	static int aaOf(const int c1, const codon &cod);
};


class codon : public alphabet {
public:
	explicit codon(); //default constructor: reads "nuclearCode.txt"
	explicit codon(const geneticCodeString& matrixFileString){readMatrixFromFile(matrixFileString.Val);};
	virtual ~codon() {}
  //	explicit codon( codon& other);
	codon& operator=(const codon& other);
	virtual alphabet* clone() const { return new codon(*this); }
	void readMatrixFromFile(const string& matrixFileName);
	const map <string,string> & geneticCode()const {return _geneticCode;}
	int unknown() const  {return 64;}
	int gap() const  {return -1;}
	int size() const {return _alphabetSize;} // 3 stop codon excluded
	int stringSize() const {return 3;} // 3 letter code.
	vector<int> fromString(const string& str) const;

	int fromChar(const string& s, const int pos) const;
	string fromInt(const int in_id) const;
	// "specific" here is not unknown, nor ambiguity, nor gap (for example, for nucleotides it will true for A,C,G, or T).
	bool isSpecific(const int id) const {return (id>=0 && id < size());}


	
  int relations(const int charInSeq, const int charToCheck) const{
		if (charInSeq == -1) {
			errorMsg::reportError("gaps in the sequences. Either change gaps to ? or remove gap positions");
		}
		else if (charInSeq == unknown()) return 1;
		else if (charInSeq == charToCheck) return 1;
		assert(charInSeq < _alphabetSize);
		return 0;
	}

private:
	map <string,string> _geneticCode; //key - codon, value - amino acid
	map <string,int> _codon2Int;//key string of codon int= integer value of codon
	int _alphabetSize;
};




#endif
