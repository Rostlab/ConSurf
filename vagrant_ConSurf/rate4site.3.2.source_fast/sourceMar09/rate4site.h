#ifndef ___RATE_4_SITE____H
#define ___RATE_4_SITE____H

#include "definitions.h"
#include "rate4siteOptions.h"
#include "alphabet.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "tree.h"

class rate4site {

public:
	explicit rate4site(int argc, char* argv[]);
	virtual ~rate4site();
private:
	void computeRate4site();
	void compute_ML_Rate4site();
	void compute_EB_EXP_Rate4site();
	void getStartingStochasticProcess();
	void getStartingEvolTreeTopology(bool bCalcDistanceTable);
	void getStartingBranchLengthsAndAlpha();

	void printrate4siteInfo(ostream& out);
	void printProcessId();
	void fillOptionsParameters(int argc, char* argv[]);
	const rate4siteOptions* _options;
	void printOptionParameters();
	void getStartingSequenceData();
	void getAttributesAndValidateThatEverythingFits();
	void printOutputTree();

	// TREE SEARCH PART
	void getStartingNJtreeNjJC();
	void getStartingNJtreeNjJC_old(); // THIS IS THE OLD VERSION OF RATE4SITE.
	void getStartingTreeFromTreeFile();
	void getStartingTreeNJ_fromDistances(const VVdouble& disTab,const vector<string>& vNames);
	void getStartingNJtreeNjMLdis();
	void getStartingMLtreeFromManyNJtrees();
	
	void normalizeRates();
	void print(ostream & out, const Vdouble & rate2print);

	void printRatesBayes(ostream& out, const Vdouble & rate2print);
    void printRatesML(ostream& out, const Vdouble & rate2print);
	void printAveAndStd(ostream& out);
	void computeAveAndStd(); // fills _ave, and _std
	void fillReferenceSequence();
	void removeGapPositionAccordingToTheReferenceSequence();

	sequenceContainer _sc;
	tree _et;
	stochasticProcess* _sp;
	alphabet* _alph;

	vector<const stochasticProcess* > _spVec;
	vector<tree> _etVec;
	Vint _spAtributes;
	Vint _etAtributes;
	sequence* _refSeq; // the reference sequence

	Vdouble _rate;// the rates themselves
	Vdouble _Lrate;// the log likelihood of each position
	Vdouble _normalizedRates; // the rates when their ave = 0 and std = 1.
	MDOUBLE _ave; // the average over all rates.
	MDOUBLE _std; // the std over all rates.

	Vdouble _BayesianSTD;// the std of the Bayesian rates
	Vdouble _BayesianLowerBound;// lower bound of rate in Bayesian inference
	Vdouble _BayesianUpperBound;// upper bound of rate in Bayesian inference

	MDOUBLE _alphaConf; // the alpha confidence interval of Bayesian rates (set to 0.5).
						// the rate confidence interval will be the range of rates
						// that are in the 95% area under the curve.
};


#endif
