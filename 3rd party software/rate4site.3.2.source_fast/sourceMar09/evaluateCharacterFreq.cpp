// $Id: evaluateCharacterFreq.cpp 391 2005-06-12 12:44:57Z ninio $

#include "evaluateCharacterFreq.h"
#include <cassert>

vector<MDOUBLE> sumAlphabetCounts(const sequenceContainer & sc) {
	vector<MDOUBLE> charFreq(sc.alphabetSize(),0.0);
	sequenceContainer::constTaxaIterator tIt;
	sequenceContainer::constTaxaIterator tItEnd;
	tIt.begin(sc);
	tItEnd.end(sc);
	while (tIt!= tItEnd) {
		sequence::constIterator sIt;
		sequence::constIterator sItEnd;
		sIt.begin(*tIt);
		sItEnd.end(*tIt);
		while (sIt != sItEnd) {
			if ((*sIt >= 0) && (*sIt <charFreq.size())) ++charFreq[(*sIt)];
			++sIt;
		}
		++tIt;
	}
	return charFreq;
}

void changeCountsToFreqs(vector<MDOUBLE>& charFreq){
	MDOUBLE sumA = 0;
	int i=0;
	for (i=0; i < charFreq.size(); ++i) {
        sumA+=charFreq[i] ;
	}
	for (i=0; i < charFreq.size(); ++i) {
		charFreq[i] /= sumA;
	}
}

void makeSureNoZeroFreqs(vector<MDOUBLE> & charFreq){
	// CORRECT SO THAT THERE ARE NO ZERO FREQUENCIES.
	// ALL FREQS THAT WERE ZERO ARE CHANGED
	MDOUBLE sumB=0;
	int i=0;
	for (i=0; i < charFreq.size(); ++i) {
        if (charFreq[i]==0) charFreq[i] = 0.00001 ;
		else sumB +=charFreq[i];
	}
	for (i=0; i < charFreq.size(); ++i) {
		if (charFreq[i] != 0.00001) charFreq[i] /= sumB;
	}
}


vector<MDOUBLE> evaluateCharacterFreq(const sequenceContainer & sc) {
	vector<MDOUBLE> charFreq=sumAlphabetCounts(sc);
	changeCountsToFreqs(charFreq);
	makeSureNoZeroFreqs(charFreq);
	return charFreq;
}

VVdouble evaluateCharacterFreqOneForEachGene(const vector<sequenceContainer> & scVec){
	VVdouble charFreq;
	for (int k=0; k < scVec.size(); ++k) {
		charFreq.push_back(evaluateCharacterFreq(scVec[k]));
	}
	return charFreq;
}

		


vector<MDOUBLE> evaluateCharacterFreqBasedOnManyGenes(const vector<sequenceContainer> & scVec) {
	// note: all alphabets have to be the same!
	vector<MDOUBLE> charFreq(scVec[0].alphabetSize(),0.0);
	for (int i=0; i < scVec.size();++i) {
		assert(scVec[0].getAlphabet()->size()==scVec[i].getAlphabet()->size());
        vector<MDOUBLE> charFreqTmp=sumAlphabetCounts(scVec[i]);
		for (int z=0; z < charFreq.size();++z) charFreq[z]+=charFreqTmp[z];
	}
	changeCountsToFreqs(charFreq);
	makeSureNoZeroFreqs(charFreq);
	return charFreq;
}
