// $Id: computeDownAlg.cpp 638 2006-03-18 19:27:55Z eyalprivman $

#include "definitions.h"
#include "computeDownAlg.h"
#include "treeIt.h"
#include "doubleRep.h"

void computeDownAlg::fillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup){
	ssc.allocatePlace(et.getNodesNum(), pi.alphabetSize());
	treeIterTopDownConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int letter,letterInFather,bro,letterInSon;
		if (mynode->father()==NULL) {// if root
			for(letter=0; letter<pi.alphabetSize();letter++) {
				ssc.set(mynode->id(),letter,1.0);
			}
			mynode = tIt.next(); //continue
		}
		tree::nodeP fatherNode=mynode->father();
		const int n_bro=fatherNode->getNumberOfSons();
		for(letter=0; letter<pi.alphabetSize();letter++) {//alpha
			doubleRep totalProb=1.0;
			doubleRep fatherTerm=0;
			if (fatherNode->father()!=NULL) {
				for(letterInFather=0; letterInFather<pi.alphabetSize();letterInFather++)
					fatherTerm += pi.getPij(fatherNode->id(),letter,letterInFather)*
					ssc.get(fatherNode->id(),letterInFather);
			}
			else {
				fatherTerm=1.0;
			}
				doubleRep brotherTerm=1.0;
			for(bro = 0; bro < n_bro; bro++) {
				tree::nodeP brother = fatherNode->getSon(bro);
				if (brother != mynode) {
					doubleRep tmp_bro=0.0;
					for(letterInSon=0; letterInSon<pi.alphabetSize();letterInSon++) {
						tmp_bro+=pi.getPij(fatherNode->getSon(bro)->id(),letter,letterInSon)*
						cup.get(brother->id(),letterInSon);
					}
					brotherTerm *=tmp_bro;
				}
			}
			totalProb = fatherTerm * brotherTerm;
			ssc.set(mynode->id(),letter,totalProb);
		}
	}
}


//use Pij(t) from the stochastic process instead of precomputed probabilities (via the computePijHom class)
void computeDownAlg::fillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup){
	ssc.allocatePlace(et.getNodesNum(), sp.alphabetSize());
	treeIterTopDownConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int letter, letterInFather, bro, letterInSon;
		if (mynode->isRoot()) {// if root: set all values to 1.0
			for(letter = 0; letter < sp.alphabetSize(); letter++) {
				ssc.set(mynode->id(), letter, 1.0);
			}
			mynode = tIt.next(); //continue
		}
		tree::nodeP fatherNode = mynode->father();
		const int n_bro = fatherNode->getNumberOfSons();
		for(letter = 0; letter < sp.alphabetSize(); letter++) {
			doubleRep totalProb=1.0;
			doubleRep fatherTerm=0;
			if (fatherNode->isRoot()) 
			{
				fatherTerm = 1.0;
			}
			else
			{
				for(letterInFather = 0; letterInFather < sp.alphabetSize(); letterInFather++)
				{
					MDOUBLE dist = fatherNode->dis2father() * sp.getGlobalRate(); 
					fatherTerm += sp.Pij_t(letter, letterInFather, dist)
					* ssc.get(fatherNode->id(), letterInFather);
				}
			}
			doubleRep brotherTerm = 1.0;
			for(bro = 0; bro < n_bro; bro++) {
				tree::nodeP brother = fatherNode->getSon(bro);
				if (brother != mynode) {
					doubleRep tmp_bro=0.0;
					for(letterInSon = 0; letterInSon < sp.alphabetSize(); letterInSon++) 
					{
						MDOUBLE dist = brother->dis2father() * sp.getGlobalRate();
						tmp_bro += sp.Pij_t(letter, letterInSon, dist)
						* cup.get(brother->id(), letterInSon);
					}
					brotherTerm *= tmp_bro;
				}
			}
			totalProb = fatherTerm * brotherTerm;
			ssc.set(mynode->id(), letter, totalProb);
		}
	}
}


//compute probabilities with a site-specific rate
void computeDownAlg::fillComputeDownSpecificRate(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const stochasticProcess& sp,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup,
					   const MDOUBLE gRate){
	ssc.allocatePlace(et.getNodesNum(), sp.alphabetSize());
	treeIterTopDownConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int letter, letterInFather, bro, letterInSon;
		if (mynode->isRoot()) {// if root: set all values to 1.0
			for(letter = 0; letter < sp.alphabetSize(); letter++) {
				ssc.set(mynode->id(), letter, 1.0);
			}
			mynode = tIt.next(); //continue
		}
		tree::nodeP fatherNode = mynode->father();
		const int n_bro = fatherNode->getNumberOfSons();
		for(letter = 0; letter < sp.alphabetSize(); letter++) {
			doubleRep totalProb=1.0;
			doubleRep fatherTerm=0;
			if (fatherNode->isRoot()) 
			{
				fatherTerm = 1.0;
			}
			else
			{
				for(letterInFather = 0; letterInFather < sp.alphabetSize(); letterInFather++)
				{
					MDOUBLE dist = fatherNode->dis2father() * gRate * sp.getGlobalRate(); 
					fatherTerm += sp.Pij_t(letter, letterInFather, dist)
					* ssc.get(fatherNode->id(), letterInFather);
				}
			}
			doubleRep brotherTerm = 1.0;
			for(bro = 0; bro < n_bro; bro++) {
				tree::nodeP brother = fatherNode->getSon(bro);
				if (brother != mynode) {
					doubleRep tmp_bro=0.0;
					for(letterInSon = 0; letterInSon < sp.alphabetSize(); letterInSon++) 
					{
						MDOUBLE dist = brother->dis2father() * gRate * sp.getGlobalRate();
						tmp_bro += sp.Pij_t(letter, letterInSon, dist)
						* cup.get(brother->id(), letterInSon);
					}
					brotherTerm *= tmp_bro;
				}
			}
			totalProb = fatherTerm * brotherTerm;
			ssc.set(mynode->id(), letter, totalProb);
		}
	}
}





