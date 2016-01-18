// $Id: simulateTree.cpp 391 2005-06-12 12:44:57Z ninio $

#include "definitions.h"
#include "treeUtil.h"
#include "simulateTree.h"
#include "talRandom.h"
#include "gammaDistribution.h"

simulateTree::simulateTree(const tree&  _inEt,
						   const stochasticProcess& sp,
						   const alphabet* alph) :
  _et(_inEt), _sp(sp),_alph(alph) {};

simulateTree::~simulateTree() {}


void simulateTree::generate_seq(int seqLength) {
	sequence justAseq(_alph);
	_simulatedSequences.resize(_et.getNodesNum(),justAseq);
	for (int i=0; i < _simulatedSequences.size(); ++i) {
		_simulatedSequences[i].resize(seqLength);
	}
	generateRootSeq(seqLength); 

	vector<MDOUBLE> rateVec(seqLength);
	for (int h = 0; h < seqLength; h++)  {
		int theRanCat = getRandCategory(h);
		rateVec[h] = _sp.rates(theRanCat);
	}
	

	for (int p=0 ; p < _et.getRoot()->getNumberOfSons() ; ++p) {
	  recursiveGenerateSpecificSeq(rateVec, seqLength, _et.getRoot()->getSon(p));
	}
}


void simulateTree::generate_seq_continuous_gamma(int seqLength) {
	sequence justAseq(_alph);
	_simulatedSequences.resize(_et.getNodesNum(),justAseq);
	for (int i=0; i < _simulatedSequences.size(); ++i) {
		_simulatedSequences[i].resize(seqLength);
	}
	generateRootSeq(seqLength); 

	vector<MDOUBLE> rateVec(seqLength);
	MDOUBLE alpha= (static_cast<gammaDistribution*>(_sp.distr()))->getAlpha();
	for (int h = 0; h < seqLength; h++)  {
	  rateVec[h] = talRandom::SampleGamma(alpha);
	}
	

	for (int p=0 ; p < _et.getRoot()->getNumberOfSons() ; ++p) {
	  recursiveGenerateSpecificSeq(rateVec, seqLength, _et.getRoot()->getSon(p));
	}
}

void simulateTree::generate_seqWithRateVector(const vector<MDOUBLE>& rateVec, 
											  const int seqLength) {
	sequence justAseq(_alph);
	_simulatedSequences.resize(_et.getNodesNum(),justAseq);
	for (int i=0; i < _simulatedSequences.size(); ++i) {
		_simulatedSequences[i].resize(seqLength);
	}
	generateRootSeq(seqLength); 

	for (int p=0 ; p < _et.getRoot()->getNumberOfSons() ; ++p) {
	  recursiveGenerateSpecificSeq(rateVec,seqLength,_et.getRoot()->getSon(p));
	}
}

void simulateTree::generateRootSeq(int seqLength) {	
	for (int i = 0; i < seqLength; i++) {
		_simulatedSequences[_et.getRoot()->id()][i] =  giveRandomChar();
	}
	_simulatedSequences[_et.getRoot()->id()].setAlphabet(_alph);
	_simulatedSequences[_et.getRoot()->id()].setName(_et.getRoot()->name());
	_simulatedSequences[_et.getRoot()->id()].setID(_et.getRoot()->id());

}

void simulateTree::recursiveGenerateSpecificSeq(
							const vector<MDOUBLE> &rateVec,
							const int seqLength,
							tree::nodeP myNode) {

	for (int y = 0; y < seqLength; y++) {
		MDOUBLE lenFromFather=myNode->dis2father()*rateVec[y];
		int aaInFather = _simulatedSequences[myNode->father()->id()][y];
		int newChar = giveRandomChar(aaInFather,lenFromFather,y);
		_simulatedSequences[myNode->id()][y] = newChar;
    }
	_simulatedSequences[myNode->id()].setAlphabet(_alph);
	_simulatedSequences[myNode->id()].setName(myNode->name());
	_simulatedSequences[myNode->id()].setID(myNode->id());
	for (int x =0 ; x < myNode->getNumberOfSons(); ++x) {
	  recursiveGenerateSpecificSeq(rateVec, seqLength, myNode->getSon(x));
	}
}

int simulateTree::giveRandomChar() const {
	for (int loop =0 ;loop<100000 ;loop++) {
		MDOUBLE theRandNum = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0);
		MDOUBLE sum = 0.0;
		for (int j=0;j<_sp.alphabetSize();++j) {
			sum+=_sp.freq(j);
			if (theRandNum<sum) return j;
		}
	} 
	// there was a bug here. It was fixed on 21.3.2005.
	// the previous code returned a character with equal probabilities...
	// return(talRandom::giveIntRandomNumberBetweenZeroAndEntry(_sp.alphabetSize()));
	errorMsg::reportError("Could not give random character. The reason is probably that the P_i do not sum to one.");
	return 1;
}

int simulateTree::giveRandomChar(const int letterInFatherNode,
								 const MDOUBLE length,
								 const int pos) const {
	assert(letterInFatherNode>=0);
	assert(letterInFatherNode<_sp.alphabetSize());
	for (int loop =0 ;loop<100000 ;loop++) {
		MDOUBLE theRandNum = talRandom::giveRandomNumberBetweenZeroAndEntry(1.0);
		MDOUBLE sum = 0.0;
		for (int j=0;j<_sp.alphabetSize();++j) {
			sum+=_sp.Pij_t(letterInFatherNode,j, length);
			if (theRandNum<sum) return j;
		}
	}
	errorMsg::reportError("Could not give random character. The reason is probably that the Pij_t do not sum to one.");
	return 1;
}


int simulateTree::getRandCategory(const int pos) const {
  MDOUBLE theRandNum = talRandom::giveRandomNumberBetweenZeroAndEntry(1);
  MDOUBLE sum = 0.0;
  for (int j=0;j<_sp.categories() ;++j) {
     sum+=_sp.ratesProb(j);
   if (theRandNum<sum) return j;
  }
  errorMsg::reportError(" error in function simulateTree::getRandCategory() ");// also quit the program
  return -1;
}

sequenceContainer simulateTree::toSeqData() {
	sequenceContainer myseqData;
	for (int i=0; i < _simulatedSequences.size(); ++i) {
		myseqData.add(_simulatedSequences[i]);
	}
	return myseqData;
}

sequenceContainer simulateTree::toSeqDataWithoutInternalNodes() {
	sequenceContainer myseqData;
	for (int i=0; i < _simulatedSequences.size(); ++i) {
		tree::nodeP theCurNode = _et.findNodeByName(_simulatedSequences[i].name());
		if (theCurNode->isInternal()) continue;
		myseqData.add(_simulatedSequences[i]);
	}
	return myseqData;
}
