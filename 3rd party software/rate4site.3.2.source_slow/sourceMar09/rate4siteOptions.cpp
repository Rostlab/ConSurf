#include "rate4siteOptions.h"
#include "errorMsg.h"
#include <iostream>
#include <cmath>
#include <stdlib.h>  //for exit()
using namespace std;

rate4siteOptions::rate4siteOptions(int& argc, char *argv[]):

// DEFAULTS VALUES:
referenceSeq("non"),
treefile(""),
seqfile(""),
logFile(""), //log.txt
logValue(3),
outFile("r4s.res"),
outFileNotNormalize("r4sOrig.res"),
treeOutFile("TheTree.txt"),
fileListingManySubstitutionMatrices(""),
spPositionFile(""),
manyTreeFileInput(""),
treePositionFile(""),
modelName(jtt),
treeSearchAlg(njML),
alphabet_size(20), // this is 20
outPtr(&cout),
optimizeBranchLengths(mlAndAlphaBBL),
rateEstimationMethod(ebExp),
numberOfDiscreteCategories(16),
userInputAlpha(0)
{


	//static struct option long_options[] =  {{0, 0, 0, 0}};
	int option_index = 0;
	int c=0;
	bool algo_set=false;

	out_f.open(outFile.c_str()); // default output file
	outPtr=&out_f; // default output file

	while (c >= 0) {
#ifdef WIN32
		c = getopt_long(argc, argv,"A:a:b:B:c:C:d:D:f:F:Hh?i:I:k:K:l:L:m:M:O:o:p:P:s:S:T:t:u:U:v:V:x:X:y:Y:z:Z:",NULL,&option_index);
#else
		c = getopt(argc, argv,"A:a:b:B:c:C:d:D:f:F:Hh?i:I:k:K:l:L:m:M:O:o:p:P:s:S:T:t:u:U:v:V:x:X:y:Y:z:Z:");
#endif
		switch (c) {

			// tree file, seqfile 
			case 'a':case 'A': referenceSeq=optarg; break;
			case 'b':case 'B': {
				switch (optarg[0]) {
					case 'g': case 'G':  optimizeBranchLengths=mlAndAlphaBBL; break;	
					case 'h': case 'H':  optimizeBranchLengths=mlBBLUniform; break;	
					case 'n': case 'N':  optimizeBranchLengths=noBBL; break;	
					default: optimizeBranchLengths=mlAndAlphaBBL; break;
				}
			} break;
			case 'c':case 'C': manyTreeFileInput=optarg; break;
			case 'd':case 'D': userInputAlpha=atof(optarg); break;
			case 'f':case 'F': fileListingManySubstitutionMatrices=optarg; break;
			case 'h':case 'H': case '?':
				cout <<"rate4site version 2.01. Last updated: September 12, 2004"<<endl;
				cout <<"USAGE:	"<<argv[0]<<" [-options] "<<endl <<endl;
				cout <<"+----------------------------------------------+"<<endl;
				cout <<"|-t    tree file                               |"<<endl;
				cout <<"|-s    seq file (accepted formats:             |"<<endl;
				cout <<"|      Fasta, Mase, Molphy, Phylip, Clustal)   |"<<endl;
				cout <<"|-o    out file                                |"<<endl;
				cout <<"|-x    tree out file                           |"<<endl;
				cout <<"|-y    outfile for un-normalize rates          |"<<endl;
				cout <<"|----------------------------------------------|"<<endl;
				cout <<"|-M     model name                             |"<<endl;
				cout <<"|-Mj    JTT                                    |"<<endl;
				cout <<"|-Mr    REV (for mitochondrial genomes)        |"<<endl;
				cout <<"|-Md    DAY                                    |"<<endl;
				cout <<"|-Mw    WAG                                    |"<<endl;
				cout <<"|-MC    cpREV (for chloroplasts genomes)       |"<<endl;
				cout <<"|-Ma    JC amino acids                         |"<<endl;
				cout <<"|-Mn    JC nucleotides                         |"<<endl;
				cout <<"|----------------------------------------------|"<<endl;
				cout <<"|-i    rate inference method                   |"<<endl;
				cout <<"|-im   ML rate inference                       |"<<endl;
				cout <<"|-ib   empirical Bayesian (eb-Exp)             |"<<endl;
				cout <<"|-k    number of categories for Bayesian rates |"<<endl;
				cout <<"|----------------------------------------------|"<<endl;
				cout <<"|-b     branch lengths optimization            |"<<endl;
				cout <<"|-bn    no branch lengths optimization         |"<<endl;
				cout <<"|-bh    optimization using a homogenoues model |"<<endl;
				cout <<"|-bg    optimization using a gamma model       |"<<endl;
				cout <<"|default: optimization using a gamma model     |"<<endl;
				cout <<"|----------------------------------------------|"<<endl;
				cout <<"|-a		reference sequence                     |"<<endl;
				cout <<"|default: first sequence in the alignment      |"<<endl;
				cout <<"|----------------------------------------------|"<<endl;
				cout <<"|-h or -? or -H     help                       |"<<endl;
				cout <<"|capital and no captial letters are ok         |"<<endl;
				cout <<"+----------------------------------------------+"<<endl;
				cout <<"|-f many substitution matrices input file:     |"<<endl;
				cout <<"|-p attributes for substitution matrices file  |"<<endl;
				cout <<"|-c many trees input file:                     |"<<endl;
				cout <<"|-u attributes for tree file                   |"<<endl;
				cout <<"+----------------------------------------------+"<<endl;
				cout <<"|-z tree search algorithm:                     |"<<endl;
				cout <<"|-zj = JC,   -zn = NJ with ML distances        |"<<endl;
				cout <<"|-zl = ML tree based on many random NJ trees   |"<<endl;
				cout <<"|-zo = JC & gaps are counted as differences    |"<<endl;
				cout <<"+----------------------------------------------+"<<endl;
				
				// these options are supported but are for Debugging!
				//cout <<"|-d    user input alpha                        |"<<endl;
				//cout <<"|-l    logfile                                 |"<<endl;
				//cout <<"|-v    log level                               |"<<endl;
				//cout <<"|----------------------------------------------|"<<endl;
			cout<<endl;	cerr<<" please press 0 to exit "; int d; cin>>d;exit (0);
			case 'i':case 'I': {
				switch (optarg[0]) {
					case 'b': case 'B': rateEstimationMethod=ebExp; break;
					case 'm': case 'M': rateEstimationMethod=mlRate; break;
					default: rateEstimationMethod=ebExp; break;
				} break;
			}
			case 'k':case 'K': numberOfDiscreteCategories=atoi(optarg); break;
			case 'l':case 'L': logFile=optarg; break;
			case 'm':case 'M':	{
				switch (optarg[0]) {
					case 'd': case 'D':  modelName=day;alphabet_size=20; break;
					case 'j': case 'J':  modelName=jtt;alphabet_size=20; break;
					case 'r': case 'R':  modelName=rev;alphabet_size=20; break;
					case 'w': case 'W':  modelName=wag;alphabet_size=20; break;
					case 'c': case 'C':  modelName=cprev;alphabet_size=20; break;
					case 'a': case 'A':  modelName=aajc;alphabet_size=20; break;
					case 'n': case 'N':  modelName=nucjc;alphabet_size=4; break;
					//case 'q': case 'Q':  modelName=customQ;alphabet_size=20; break;
					//case 'm': case 'M':  modelName=manyQ;alphabet_size=20; break;
					default:modelName=jtt;alphabet_size=20;
					break;
				}
			} break;
			case 'o':case 'O': {
				out_f.close(); // closing the default
				outFile=optarg;
				out_f.open(outFile.c_str());
				if (out_f == NULL) errorMsg::reportError(" unable to open output file for writing. ");
				outPtr=&out_f;
			}; break;
			case 'p':case 'P': spPositionFile=optarg; break;
			case 's':case 'S': seqfile=optarg; break;
			case 't':case 'T': treefile=optarg; break;
			case 'u':case 'U': treePositionFile=optarg; break;
			case 'v':case 'V': logValue=atoi(optarg); break;
			case 'x':case 'X': treeOutFile=optarg; break;
			case 'y':case 'Y': outFileNotNormalize=optarg; break;
			case 'z':case 'Z': {
				switch (optarg[0]) {
					case 'J': case 'j':  treeSearchAlg=njJC; break;
					case 'n': case 'N':  treeSearchAlg=njML; break;
					case 'O': case 'o':  treeSearchAlg=njJCOLD; break;
					default:treeSearchAlg=njJC;
					break;
				}
			} break;

		}
	}
}

