
//We change the code, so that it takes into account secondary structures.
//For this, we added the following members.
//(1) string userSubstitutionMatrixFile;
//This is used if the user want to input a single substitution file matrix.
//For example, if wag was not implemented, the user could input the wag substitution file.
//(2) string fileListingManySubstitutionMatrices;
//This file is used if we want to have several matrices (for example, alpha, beta, and loop matrices).
//Later, this file will be paired with the position file, so that position with attribute 1 will use
//matrix 1, etc.
//(3) string positionFile. This file gives the attributes of each position. For example
//position with attributes 1 will use substitution matrix 1, etc.
//(4) string manyTreeFileInput. If we want a different tree for each attribute.

#if !defined ___RATEFPRSITE__OPTION__T__
#define ___RATEFPRSITE__OPTION__T__

#ifdef SunOS
  #include <unistd.h>
#else
	#ifndef __STDC__
	#define __STDC__ 1
	#include "getopt.h"
	#undef __STDC__
	#else
	#include "getopt.h"
	#endif
#endif

#include <string>
#include <fstream>
using namespace std;

class rate4siteOptions{
public:
	enum modelNameType {rev,jtt,day,aajc,nucjc,wag,cprev};//,customQ,manyQ
	enum optimizeBranchLengthsType {noBBL,mlBBLUniform,mlAndAlphaBBL};
	enum treeSearchAlgType {njJC,njML,njJCOLD};
	enum rateEstimationMethodType {ebExp, mlRate};
public:
	explicit rate4siteOptions(int& argc, char *argv[]);
	ostream& out() const {return *outPtr;};
	//ostream& outNotNormalize() const {return *outPtrNotNormalize;};
	
	string treefile;
	string seqfile;
	string logFile;
	string referenceSeq; // the results are printed with this seq in each positions.
	int logValue;
	string outFile;
	string outFileNotNormalize;
	string treeOutFile;
	
	// different substitution matrices for different positions
	string fileListingManySubstitutionMatrices;//f
	string spPositionFile; //p
	
	// different trees for different positions
	string manyTreeFileInput;//c
	string treePositionFile; //u
	
	modelNameType modelName;
	treeSearchAlgType treeSearchAlg;
	int alphabet_size;
	
	optimizeBranchLengthsType optimizeBranchLengths;
	rateEstimationMethodType rateEstimationMethod;

	int numberOfDiscreteCategories;
	double userInputAlpha;
private:
	ostream* outPtr;
	ofstream out_f;
  //ostream* outPtrNotNormalize;
  //ofstream out_fNotNormalize;
};


#endif
