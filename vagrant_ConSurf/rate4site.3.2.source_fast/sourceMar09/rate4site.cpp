#include "rate4site.h"
#include "errorMsg.h"
#include "nucleotide.h"
#include "amino.h"
#include "sequenceContainer.h"
#include "uniDistribution.h"
#include "gammaDistribution.h"
#include "readDatMatrix.h"
#include "chebyshevAccelerator.h"
#include "nucJC.h"
#include "aaJC.h"
#include "trivialAccelerator.h"
#include "jcDistance.h"
#include "distanceTable.h"
#include "nj.h"
#include "numRec.h"
#include "checkcovFanctors.h"
#include "someUtil.h"
#include "likeDist.h"
#include "fastStartTree.h"
#include "checkcovFanctorsWithFactors.h"
#include "bestAlpha.h"
#include "siteSpecificRate.h"
#include "likelihoodComputation.h"
#include "recognizeFormat.h"
#include "bblEM.h"
#include <iostream>
#include <cassert>
#include <iomanip>
#include <ctime>
#include <stdio.h> 
#include <string.h>


//#ifdef unix	//for process id
//#include <unistd.h>
//#else
//#include <process.h>
//#endif

using namespace std;
using namespace likelihoodComputation;

//ofstream timingsF;

int main(int argc, char* argv[]) {
	rate4site r4s(argc,argv);
	return 0;
}

rate4site::rate4site(int argc, char* argv[]) {
	//timingsF.open("timing.dat");
	fillOptionsParameters(argc,argv);
	myLog::setLog(_options->logFile, 5);
	printProcessId();
	printrate4siteInfo(cout);
	printOptionParameters();
	getStartingSequenceData();
	getStartingStochasticProcess();// must be b4 the tree!
	getStartingEvolTreeTopology(true);// without optimizing branch lengths
	getStartingBranchLengthsAndAlpha();
	printOutputTree();
	fillReferenceSequence();
	removeGapPositionAccordingToTheReferenceSequence();
	getAttributesAndValidateThatEverythingFits();
	computeRate4site();
	
	ofstream nonNormalizedOutStream(_options->outFileNotNormalize.c_str());
	computeAveAndStd(); // put them in _ave, and _std
	print(nonNormalizedOutStream,_rate);
	nonNormalizedOutStream.close();
	normalizeRates(); // change also the _ave, the _std the quantiles, etc.
	print(_options->out(),_normalizedRates);
	//timingsF.close();
}

void rate4site::removeGapPositionAccordingToTheReferenceSequence(){
	// we remove unknown and not gaps because we changed all gaps to unknown...
	_sc.removeUnknownPositionsAccordingToAReferenceSeq(_refSeq->name());
}

void rate4site::print(ostream & out, const Vdouble & rate2print) {
	switch (_options->rateEstimationMethod){
		case (rate4siteOptions::ebExp):  
            printRatesBayes(out,rate2print);
			break;
		case (rate4siteOptions::mlRate):
           printRatesML(out,rate2print);
			break;
	}
	printAveAndStd(out);
}

void rate4site::computeAveAndStd(){
	MDOUBLE sum = 0;
	MDOUBLE sumSqr=0.0;
	for (int i=0; i < _sc.seqLen(); ++i) {
		sum+=_rate[i];
		sumSqr+=(_rate[i]*_rate[i]);
	}
	_ave = sum/_sc.seqLen();
	_std= sumSqr-(sum*sum/_sc.seqLen());
	_std /= (_sc.seqLen()-1.0);
	_std = sqrt(_std);
	if (((_ave<1e-9)) && (_ave>(-(1e-9)))) _ave=0;
	if ((_std>(1-(1e-9))) && (_std< (1.0+(1e-9)))) _std=1.0;
}

void rate4site::fillReferenceSequence(){
	if (strcmp(_options->referenceSeq.c_str(),"non")==0) {
		_refSeq = &(_sc[0]);
	}
	else {
		int id1 = _sc.getId(_options->referenceSeq,true);
		_refSeq = (&_sc[id1]);
	}
}

rate4site::~rate4site() {
	delete _alph;
	for (int i=0; i < _spVec.size(); ++i) delete _spVec[i];
	delete _options;
}

void rate4site::printrate4siteInfo(ostream& out) {
	out<<endl;
	out<<" ======================================================="<<endl;
	out<<" the rate for site project:                             "<<endl;
	out<<" Version: 2.01. Last updated 6.11.06                     "<<endl;
	out<<" Tal Pupko and his lab:     talp@post.tau.ac.il         "<<endl;
	out<<" Nir Ben-Tal and his lab:   bental@ashtoret.tau.ac.il   "<<endl;
	out<<" Itay Mayrose:  itayMay@post.tau.ac.il                  "<<endl;
	out<<" For program support, please contact Itay Mayrose       "<<endl;
	out<<" ======================================================="<<endl;
	out<<endl;
}

void rate4site::printProcessId()
{
	//int pp = getpid();
	//LOG(4,<<"#Process_id= "<<pp<<endl);
}
void rate4site::fillOptionsParameters(int argc, char* argv[]) {
	_options = new rate4siteOptions(argc, argv);
	_alphaConf = 0.5;
}

void rate4site::printOptionParameters() {
	cout<<"\n ---------------------- THE PARAMETERS ----------------------------"<<endl;
	if (_options->treefile.size()>0) cout<<"tree file is: "<<_options->treefile<<endl;
	if (_options->seqfile.size()>0) cout<<"seq file is: "<<_options->seqfile<<endl;
	if (_options->outFile.size()>0) cout<<"output file is: "<<_options->outFile<<endl;
	if 	(strcmp(_options->referenceSeq.c_str(),"non")!=0) cout<<"reference sequence is: "<<_options->referenceSeq<<endl;
	switch (_options->rateEstimationMethod){
		case (rate4siteOptions::ebExp):  {
			cout<< "rate inference method is: empirical Bayesian estimate"<<endl;
			cout<< "using a Gamma prior distribution with: "<<_options->numberOfDiscreteCategories<< " discrete categories"<<endl;
				}break;
		case (rate4siteOptions::mlRate): cout<< "rate inference method is: maximum likelihood (ML) "<<endl;break;
	}
	switch (_options->modelName){
		case (rate4siteOptions::day): cout<< "probablistic_model is: DAY" <<endl; break;
		case (rate4siteOptions::jtt): cout<< "probablistic_model is: JTT" <<endl; break;
		case (rate4siteOptions::rev): cout<< "probablistic_model is: REV" <<endl; break;
		case (rate4siteOptions::aajc): cout<< "probablistic_model is: AAJC" <<endl; break;
		case (rate4siteOptions::nucjc): cout<< "probablistic_model is: NUCJC" <<endl; break;
		case (rate4siteOptions::wag): cout<< "probablistic_model is: WAG" <<endl; break;
		case (rate4siteOptions::cprev): cout<< "probablistic_model is: CPREV" <<endl; break;
	}

	switch (_options->optimizeBranchLengths){
		case (rate4siteOptions::noBBL): cout<<"branch lengths are not optimized"<<endl; break;
		case (rate4siteOptions::mlBBLUniform): cout<<"branch lengths optimization is ML using a homogenoues model"<<endl; break;
		case (rate4siteOptions::mlAndAlphaBBL): cout<<"branch lengths optimization is ML using a gamma model"<<endl; break;
	}
	cout<<"\n -----------------------------------------------------------------"<<endl;
}

void rate4site::getStartingSequenceData(){
	if (_options->seqfile == "") {
		errorMsg::reportError("Please give a sequence file name in the command line");
	}
	ifstream in(_options->seqfile.c_str());
	int alphabetSize = _options->alphabet_size;
	if (alphabetSize==4) _alph = new nucleotide;
	else if (alphabetSize == 20) _alph = new amino;
	else errorMsg::reportError("no such alphabet in function rate4site::getStartingSequenceData");

	sequenceContainer original = recognizeFormat::read(in,_alph);;
	original.changeGaps2MissingData();
	_sc = original;
}

//	void makeSureTreeIsBiFurcatingAndUnrooted();
void rate4site::getStartingStochasticProcess() {
	if (_options->numberOfDiscreteCategories<1 || _options->numberOfDiscreteCategories>50) {
		errorMsg::reportError("number of discrete rate categories should be between 1 and 50");
	}
	distribution *dist = NULL;
	switch (_options->rateEstimationMethod){
		case (rate4siteOptions::mlRate): dist =  new uniDistribution; break;
		case (rate4siteOptions::ebExp): dist =  new gammaDistribution(1,_options->numberOfDiscreteCategories); break;
		default: dist =  new gammaDistribution(1,_options->numberOfDiscreteCategories); break;
	}

	replacementModel *probMod=NULL;
	pijAccelerator *pijAcc=NULL;
	switch (_options->modelName){
		case (rate4siteOptions::day):
			probMod=new pupAll(datMatrixHolder::dayhoff);pijAcc = new chebyshevAccelerator(probMod); break;
		case (rate4siteOptions::jtt):
			probMod=new pupAll(datMatrixHolder::jones); pijAcc = new chebyshevAccelerator(probMod); break;
		case (rate4siteOptions::rev):
			probMod=new pupAll(datMatrixHolder::mtREV24); pijAcc = new chebyshevAccelerator(probMod); break;
		case (rate4siteOptions::wag):
			probMod=new pupAll(datMatrixHolder::wag); pijAcc = new chebyshevAccelerator(probMod); break;
		case (rate4siteOptions::cprev):
			probMod=new pupAll(datMatrixHolder::cpREV45); pijAcc = new chebyshevAccelerator(probMod); break;
		case (rate4siteOptions::nucjc):
			probMod=new nucJC; pijAcc = new trivialAccelerator(probMod); break;
		case (rate4siteOptions::aajc):
			probMod=new aaJC; pijAcc = new trivialAccelerator(probMod); break;
		//case (rate4siteOptions::customQ):
		//	probMod=new pupAll(_options->userSubstitutionMatrixFile);pijAcc = new chebyshevAccelerator(probMod); break;
		default:
			errorMsg::reportError("this probablistic model is not yet available");
	}
	_sp = new stochasticProcess(dist, pijAcc);
	if (probMod) delete probMod;
	if (pijAcc) delete pijAcc;
	if (dist) delete dist;
}

//if bCalcDistanceTable==true then always calculate the distance table
void rate4site::getStartingEvolTreeTopology(bool bCalcDistanceTable){
	time_t ltime1;
	time( &ltime1 );
	LOG(4,<<"get Starting Tree Topology"<<endl);
	VVdouble disTab;
	vector<string> vNames;
	if ((_options->treefile=="") || bCalcDistanceTable) {
		distanceMethod* pDm;
		switch (_options->treeSearchAlg){
			case (rate4siteOptions::njJC):
					pDm = new jcDistance(_options->alphabet_size);
					giveDistanceTable(pDm, _sc, disTab, vNames);
				//getStartingNJtreeNjJC(); 
				break;
			case (rate4siteOptions::njJCOLD):
				pDm = new jcDistanceOLD(_options->alphabet_size);
				giveDistanceTable(pDm, _sc,disTab, vNames);

				//getStartingNJtreeNjJC_old(); 
				break;
			case (rate4siteOptions::njML): {
					uniDistribution lUni;
					const pijAccelerator* lpijAcc = _sp->getPijAccelerator();// note this is just a copy of the pointer.
					stochasticProcess lsp(&lUni,lpijAcc);
					pDm = new likeDist(lsp,0.01);
					giveDistanceTable(pDm,_sc,disTab,vNames);
					//getStartingNJtreeNjMLdis();
				}
				break;
			default:
				errorMsg::reportError("this tree search mode is not yet available");
		}
		delete pDm;

		//calc distqance table statistics
		MDOUBLE low_bound = VERYBIG;
		MDOUBLE upper_bound = VERYSMALL;
		MDOUBLE sum = 0.0;
		int count = 0;
		for (int i = 0; i < disTab.size(); ++i){
			for (int j = i+1; j < disTab[i].size(); ++j){
				sum += disTab[i][j];
				++count;
				if (disTab[i][j] < low_bound)
					low_bound = disTab[i][j];
				if (disTab[i][j] > upper_bound)
					upper_bound = disTab[i][j];
			}
		}
		MDOUBLE avg = sum / static_cast<MDOUBLE>(count);
		LOG(4,<<"#MSA diversity matrix"<<endl);
		LOG(4,<<"#Average pairwise distance= "<<avg<<endl);
		LOG(4,<<"#lower bound = "<<low_bound<<endl);
		LOG(4,<<"#upper bound = "<<upper_bound<<endl);
		LOG(4,<<"#end of MSA diversity matrix"<<endl);
	}
	if (_options->treefile=="") {
		//take ditance table from above
		getStartingTreeNJ_fromDistances(disTab, vNames);
	}
	else 
		getStartingTreeFromTreeFile();

	LOG(4,<<"After Tree Topology"<<endl);
	time_t ltime2;
	time( &ltime2 );
	int t = ltime2 - ltime1;
	//timingsF<<"time for tree topology = "<<t<<endl;
}

void rate4site::getStartingBranchLengthsAndAlpha(){
	time_t ltime1;
	time( &ltime1 );
	LOG(4,<<"get Starting Branch Lengths And Alpha"<<endl);
	int maxBBLIterations = 5;
	int maxTotalAlphaBBLIterations = 2;
	MDOUBLE epsilonForBBL= 0.1;
	MDOUBLE epsilonForAlpha= 0.2;
	MDOUBLE upperBoundAlpha = 5.0;
	MDOUBLE intitalAlpha = 1.0;
	cerr<< "Optimizing branch lengths and alpha..."<<endl;
	if (_options->rateEstimationMethod == rate4siteOptions::mlRate) {
		if (_options->optimizeBranchLengths == rate4siteOptions::noBBL) {
			return;
		} else if (_options->optimizeBranchLengths == rate4siteOptions::mlBBLUniform) {
			bblEM bblEM1(_et, _sc, *_sp, NULL, maxBBLIterations , epsilonForBBL, epsilonForBBL);
		} else {
			// Here we want to optimize branch lengths with a gamma model,
			// but sp is with a homogenoues model. Hence, we have to create a local
			// copy of a gamma stochastic process.
			if (_options->userInputAlpha != 0) intitalAlpha = _options->userInputAlpha;
			gammaDistribution localDist(intitalAlpha,_options->numberOfDiscreteCategories);
			stochasticProcess localSP(&localDist,_sp->getPijAccelerator());
			if (_options->userInputAlpha == 0) {
				// in this case we have to optimize both the alpha and the branche lengths
				bestAlphaAndBBL bbl1(_et, _sc, localSP, NULL, intitalAlpha, upperBoundAlpha, epsilonForAlpha, epsilonForBBL, maxBBLIterations, maxTotalAlphaBBLIterations);
			} else {
				// in this case we know the alpa, and we want to just optimize branch lengths with this alpha
				bestAlphaAndBBL bbl(_et, _sc, localSP, NULL, intitalAlpha, upperBoundAlpha, epsilonForAlpha, epsilonForBBL, maxBBLIterations, maxTotalAlphaBBLIterations);
			}
		}
	} else { // method for inference is Bayesian
		if (_options->optimizeBranchLengths == rate4siteOptions::noBBL) {
			//FIND BEST ALPHA, AND RETURN WITHOUT CHANING THE TREE
			if (_options->userInputAlpha == 0){
				bestAlphaFixedTree bbl2(_et, _sc, *_sp, NULL, upperBoundAlpha, epsilonForAlpha);
			} else {// in this case we just want to set the alpha to the right one
				static_cast<gammaDistribution*>(_sp->distr())->setAlpha(_options->userInputAlpha);
			}
		} else if (_options->optimizeBranchLengths == rate4siteOptions::mlBBLUniform) {
			//FIND TREE WITHOUT ALPHA with an homogenoues model. Update
			uniDistribution lUni;
			const pijAccelerator* lpijAcc = _sp->getPijAccelerator();// note this is just a copy of the pointer.
			stochasticProcess lsp(&lUni,lpijAcc);
			bestAlphaAndBBL bbl(_et, _sc, lsp, NULL, intitalAlpha,  upperBoundAlpha, epsilonForAlpha, epsilonForBBL, maxBBLIterations, maxTotalAlphaBBLIterations);
			//THEN FIND ALPHA WITHOUT OPT TREE
			if (_options->userInputAlpha == 0){
				bestAlphaFixedTree bbl3(_et,_sc,*_sp, NULL, upperBoundAlpha, epsilonForAlpha);
			} else {
				static_cast<gammaDistribution*>(_sp->distr())->setAlpha(_options->userInputAlpha);
			}
		} else {
			//ML OPT WITH GAMMA
			if (_options->userInputAlpha == 0){
				bestAlphaAndBBL bbl1(_et, _sc, *_sp, NULL, intitalAlpha, upperBoundAlpha, epsilonForAlpha, epsilonForBBL, maxBBLIterations, maxTotalAlphaBBLIterations);
			} else {// alpha is known
				static_cast<gammaDistribution*>(_sp->distr())->setAlpha(_options->userInputAlpha);
				bestAlphaAndBBL bbl1(_et, _sc, *_sp, NULL, intitalAlpha, upperBoundAlpha, epsilonForAlpha, epsilonForBBL, maxBBLIterations, maxTotalAlphaBBLIterations);
			}
		}
	}
	LOG(4,<<"After Branch Lengths And Alpha"<<endl);
	time_t ltime2;
	time( &ltime2 );
	int t = ltime2 - ltime1;
	//timingsF<<"time for alpha and branch lengths optimization = "<<t<<endl;
}


void rate4site::getStartingNJtreeNjJC() {
	jcDistance pd1(_options->alphabet_size);
	VVdouble disTab;
	vector<string> vNames;
	giveDistanceTable(&pd1,
					   _sc,
					   disTab,
					   vNames);
	getStartingTreeNJ_fromDistances(disTab,vNames);
} 

void rate4site::getStartingNJtreeNjJC_old() {
	jcDistanceOLD pd1(_options->alphabet_size);
	VVdouble disTab;
	vector<string> vNames;
	giveDistanceTable(&pd1,
					   _sc,
					   disTab,
					   vNames);
	getStartingTreeNJ_fromDistances(disTab,vNames);
} 

void rate4site::getStartingNJtreeNjMLdis() {
	// note that here ALWAYS, the ML distances are computed using
	// an homogenous rate distribution.
	uniDistribution lUni;
	const pijAccelerator* lpijAcc = _sp->getPijAccelerator();// note this is just a copy of the pointer.
	stochasticProcess lsp(&lUni,lpijAcc);

	likeDist pd1(lsp,0.01);
	VVdouble disTab;
	vector<string> vNames;
	giveDistanceTable(&pd1,
					   _sc,
					   disTab,
					   vNames);
	getStartingTreeNJ_fromDistances(disTab,vNames);
}

void rate4site::getStartingMLtreeFromManyNJtrees() {
	 int numOfNJtrees = 30;
	 if (_sc.numberOfSeqs() <4) numOfNJtrees = 1;
	 else if (_sc.numberOfSeqs() <5) numOfNJtrees = 3;
	 else if (_sc.numberOfSeqs() <6) numOfNJtrees = 15;
	 else if (_sc.numberOfSeqs() <30) numOfNJtrees = 75;
	 else if (_sc.numberOfSeqs() <50) numOfNJtrees = 15;
	 else numOfNJtrees = 5;

	 const MDOUBLE tmpForStartingTreeSearch = 1;
	 const MDOUBLE epslionWeights = 0.05;
	_et = getBestMLTreeFromManyNJtrees(_sc,
								*_sp,
								numOfNJtrees,
								tmpForStartingTreeSearch,
								epslionWeights,
								cerr);
	cerr<<"number of tree evaluated: "<<numOfNJtrees<<endl;
}

void rate4site::printOutputTree() {
	ofstream f;
	string fileName1=_options->treeOutFile;
	f.open(fileName1.c_str());
	_et.output(f);
	f.close();
	cout<<"The tree was written to a file name called "<<fileName1<<endl;
}
void rate4site::getStartingTreeNJ_fromDistances(const VVdouble& disTab,
	const vector<string>& vNames) {
	NJalg nj1;
	_et= nj1.computeTree(disTab,vNames);
	ofstream f;
	string fileName1=_options->treeOutFile;
	f.open(fileName1.c_str());
	_et.output(f);
	f.close();
}

void rate4site::getStartingTreeFromTreeFile(){
	_et= tree(_options->treefile);
	if (!_et.withBranchLength()) _et.createFlatLengthMatrix(0.05);
}

void rate4site::computeRate4site(){
	time_t ltime1;
	time( &ltime1 );
	LOG(4,<<"computing the rates"<<endl);
	cerr<< "Computing the rates..."<<endl;
	if (_options->rateEstimationMethod == rate4siteOptions::ebExp) {
		compute_EB_EXP_Rate4site();
	}
	else if (_options->rateEstimationMethod == rate4siteOptions::mlRate) {
		compute_ML_Rate4site();
	}
	else 
		errorMsg::reportError("non such method for rate inference, in function void rate4site::computeRate4site()");
	LOG(4,<<"After computing the rates"<<endl);
	time_t ltime2;
	time( &ltime2 );
	int t = ltime2 - ltime1;
	//timingsF<<"time for rates calculations= "<<t<<endl;
}

void rate4site::compute_ML_Rate4site(){
	MDOUBLE maxRate = 20.0;
	MDOUBLE tol = 0.0001f;
	if (_etVec.size() > 0) { // many trees
		if (_spVec.size()>0) {
			computeML_siteSpecificRate(_rate,_Lrate,_spAtributes,_etAtributes,_etVec,_spVec,_sc,maxRate,tol);
		} else {
			computeML_siteSpecificRate(_rate,_Lrate,_etAtributes,_etVec,*_sp,_sc,maxRate,tol);
		}
	} else {// one tree
		if (_spVec.size()>0) {
			computeML_siteSpecificRate(_rate,_Lrate,_spAtributes,_et,_spVec,_sc,maxRate,tol);
		} else {
			computeML_siteSpecificRate(_rate,_Lrate,_sc,*_sp,_et,maxRate,tol);
		}
	}
}

void rate4site::compute_EB_EXP_Rate4site(){
	if (_etVec.size() > 0) { // many trees
		if (_spVec.size()>0) {
			computeEB_EXP_siteSpecificRate(_rate,_BayesianSTD,_BayesianLowerBound,_BayesianUpperBound,_spAtributes,_etAtributes,_sc,_etVec,_spVec,_alphaConf);
		} else {// MANY TREES, 1 SP
			computeEB_EXP_siteSpecificRate(_rate,_BayesianSTD,_BayesianLowerBound,_BayesianUpperBound,_etAtributes,_sc,_etVec,*_sp,_alphaConf);
		}
	} else {// one tree
		if (_spVec.size()>0) {// many sp
			computeEB_EXP_siteSpecificRate(_rate,_BayesianSTD,_BayesianLowerBound,_BayesianUpperBound,_spAtributes,_sc,_et,_spVec,_alphaConf);
		} else {//1 tree 1 sp
			computeEB_EXP_siteSpecificRate(_rate,_BayesianSTD,_BayesianLowerBound,_BayesianUpperBound,_sc,*_sp,_et,_alphaConf);
		}
	}
}


void rate4site::normalizeRates() {
	int i=0;
	if (_std==0) errorMsg::reportError(" std = 0 in function normalizeRates",1);
	_normalizedRates.resize(_sc.seqLen(),0.0);
	for (i=0;i<_normalizedRates.size();++i) {
		_normalizedRates[i]=(_rate[i]-_ave)/_std;
	}

	if (_options->rateEstimationMethod == rate4siteOptions::ebExp) {
		for (int k=0; k < _sc.seqLen(); ++k) {
			_BayesianUpperBound[k] = (_BayesianUpperBound[k] - _ave)/_std;
			_BayesianLowerBound[k] = (_BayesianLowerBound[k] - _ave)/_std;
			_BayesianSTD[k] = (_BayesianSTD[k])/_std;
		}
	}
	_ave = 0.0;
	_std = 1.0;
}

void rate4site::printRatesML(ostream& out, const Vdouble & rate2print) {
	out<<"#Rates were calculated using Maximim Likelihood"<<endl;
	out<<"#SEQ: the amino acid in the reference sequence in one letter code."<<endl;
	out<<"#SCORE: The conservation scores. lower value = higher conservation."<<endl;
	out<<"#MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position."<<endl<<endl;
	out<<"#POS"<<" ";
	out<<"SEQ";
	out<<"  ";
	out<<"SCORE"<<"     ";
    out<<"MSA DATA"<<endl; // note position start from 1.

	#ifdef unix
		for (int pos=0; pos < _sc.seqLen(); ++pos) {
			out<<setw(5)<<pos+1<<" ";
			out<<setw(5)<<_refSeq->getAlphabet()->fromInt((*_refSeq)[pos])<<" ";
			out<<setprecision(4)<<setw(7)<<rate2print[pos]<<" ";
			out<<setw(4)<<_sc.numberOfSequencesWithoutGaps(pos)<<"/"<<_sc.numberOfSeqs()<<endl; // note position start from 1.
		}	
	#else
		for (int pos=0; pos < _sc.seqLen(); ++pos) {
			out<<left<<setw(5)<<pos+1;
			out<<left<<setw(5)<<_refSeq->getAlphabet()->fromInt((*_refSeq)[pos]);
			out<<left<<setprecision(4)<<setw(7)<<fixed<<rate2print[pos];
			out<<right<<setw(4)<<_sc.numberOfSequencesWithoutGaps(pos)<<"/"<<_sc.numberOfSeqs()<<endl; // note position start from 1.
		}
	#endif

}

void rate4site::printRatesBayes(ostream& out, const Vdouble & rate2print) {
	out<<"#Rates were calculated using the expectation of the posterior rate distribution"<<endl;
	out<<"#Prior distribution is Gamma with "<<_options->numberOfDiscreteCategories<<" discrete categories"<<endl<<endl;
	out<<"#SEQ: the amino acid in the reference sequence in one letter code."<<endl;
	out<<"#SCORE: The conservation scores. lower value = higher conservation."<<endl;
	out<<"#QQ-INTERVAL: the confidence interval for the rate estimates. The default interval is 25-75 percentiles"<<endl;
	out<<"#STD: the standard deviation of the posterior rate distribution."<<endl;
	out<<"#MSA DATA: The number of aligned sequences having an amino acid (non-gapped) from the overall number of sequences at each position."<<endl;
	
	out<<endl;
	out<<"#POS"<<" ";
	out<<"SEQ";
	out<<"  ";
	out<<"SCORE"<<"    ";
	out<<"QQ-INTERVAL"<<"     ";
	out<<"STD"<<"      ";
    out<<"MSA DATA"<<endl; // note position start from 1.
	out<<"#The alpha parameter "<<static_cast<gammaDistribution*>(_sp->distr())->getAlpha()<<endl;
	
	#ifdef unix	
		for (int pos=0; pos < _sc.seqLen(); ++pos) {
			out<<setw(5)<<pos+1<<" ";
			out<<setw(5)<<_refSeq->getAlphabet()->fromInt((*_refSeq)[pos])<<" ";
			out<<setprecision(4)<<setw(7)<<rate2print[pos]<<" ";
			out<<"  ["<<setw(6)<<setprecision(4)<<_BayesianLowerBound[pos]<<","<<setprecision(4)<<setw(6)<<_BayesianUpperBound[pos]<<"]"<<" ";
			out<<setprecision(4)<<setw(7)<<_BayesianSTD[pos]<<" ";
			out<<setw(4)<<_sc.numberOfSequencesWithoutGaps(pos)<<"/"<<_sc.numberOfSeqs()<<endl; // note position start from 1.
		}
	#else
		for (int pos=0; pos < _sc.seqLen(); ++pos) {
			out<<left<<setw(5)<<pos+1;
			out<<left<<setw(5)<<_refSeq->getAlphabet()->fromInt((*_refSeq)[pos]);
			out<<left<<setprecision(4)<<setw(7)<<fixed<<rate2print[pos];
			out<<right<<"  ["<<setw(6)<<setprecision(4)<<left<<_BayesianLowerBound[pos]<<","<<setprecision(4)<<setw(6)<<right<<_BayesianUpperBound[pos]<<"]";
			out<<right<<setprecision(4)<<setw(7)<<_BayesianSTD[pos];
			out<<right<<setw(4)<<_sc.numberOfSequencesWithoutGaps(pos)<<"/"<<_sc.numberOfSeqs()<<endl; // note position start from 1.
		}
	#endif
}


void rate4site::printAveAndStd(ostream& out) {
	out<<"#Average = "<<_ave<<endl;
	out<<"#Standard Deviation = "<<_std<<endl;
}

void rate4site::getAttributesAndValidateThatEverythingFits(){
	// This is the case when we have a file with substitution matrices attributes
	// for each position. And another file, with the substitution matrices themselves.
	// Thus, for example, if the substitution matrices attribute file contains the line:
	// 2,1,3,2
	// and the substitution matrix file contain the lines
	// alpha.dat.q
	// beta.dat.q
	// loop.dat.q
	// Then the first position will work with the beta matrix, the second 
	// with the alpha matrix, etc.

	// if we don't have attributes - we don't use this part.
	
	if ((_options->spPositionFile=="") && (_options->treePositionFile==""))
		return;

	// There are two types of attributes. One is a substitution matrix for each site
	// and the other one is a different tree for each site.
	// Here is what we do in case of different subs. matrices for each site.

	if (_options->spPositionFile!="") {
		
		// Creates the stochastic processes vector
		distribution *dist = _sp->distr();
		pijAccelerator* pijAcc;
		vector<string> vs;
		ifstream in(_options->fileListingManySubstitutionMatrices.c_str());
		putFileIntoVectorStringArray(in,vs);
		in.close();
		replacementModel *probMod=NULL;
		for (int i=0; i < vs.size(); ++i) {
			probMod=new pupAll(vs[i]);
			pijAcc = new chebyshevAccelerator(probMod);
			stochasticProcess* temp_sp = new stochasticProcess(dist, pijAcc);
			_spVec.push_back(temp_sp);
			if (probMod) delete probMod;
			if (pijAcc) delete pijAcc;
		}
	
		// Reads the stochastic process attributes
		ifstream in3;
		in3.open(_options->spPositionFile.c_str());
		if (!in3) {
			errorMsg::reportError("can't open substitution matrices attributes file");
		}
		int tmpI;
		int maxAttribute = 0;
		while (!in3.eof()) {
			in3>>tmpI;
			if (in3.fail()) break;
			_spAtributes.push_back(tmpI);
			if (tmpI > maxAttribute) maxAttribute = tmpI;
		}
		in3.close();
		//check if the maximim sthochastic process attribute doesn't exceed the number of possible sthochastic processes 
	    if (_spVec.size()< maxAttribute) {
			errorMsg::reportError("There are more attributes than stochastic processes");
		}
		//check if there is an attrute for every site
		if (_refSeq->seqLen() != _spAtributes.size()) {
			errorMsg::reportError("attributes size is different from the sequence length");
		}
	}// end of updating vector of stochastic processes.

	//	This part handles many trees
	if (_options->treePositionFile!="") {
		//read the trees
		ifstream in2;
		vector<string> treeStrings;
		putFileIntoVectorStringArray(in2,treeStrings);
		in2.close();
		for (int i=0; i < treeStrings.size();++i) {
			tree tmp(treeStrings[i]);
			if (!tmp.withBranchLength()) {
				tmp.createFlatLengthMatrix(0.05);
			}
			_etVec.push_back(tmp);
		}
		// Reads the trees attributes
		ifstream in4;
		in4.open(_options->treePositionFile.c_str());
		if (!in4) {
			errorMsg::reportError("can't open trees attributes file");
		}
		int tmpI;
		int maxAttribute = 0;
		while (!in4.eof()) {
			in4>>tmpI;
			if (in4.fail()) break;
			_etAtributes.push_back(tmpI);
			if (tmpI > maxAttribute) maxAttribute = tmpI;
		}
		in4.close();
		//check if the maximim sthochastic process attribute doesn't exceed the number of possible sthochastic processes 
	    if (_etVec.size()< maxAttribute) {
			errorMsg::reportError("There are more attributes than trees");
		}
		//check if there is a tree attrute for every site
		if (_refSeq->seqLen() != _etAtributes.size()) {
			errorMsg::reportError("tree attributes size is different from the sequence length");
		}
	}
}




