// $Id: fastStartTree.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___FAST_START_TREE
#define ___FAST_START_TREE

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"
#include <iostream>

using namespace std;



tree getBestMLTreeFromManyNJtrees(sequenceContainer & allTogether,
								stochasticProcess& sp,
								const int numOfNJtrees,
								const MDOUBLE tmpForStartingTreeSearch,
								const MDOUBLE epslionWeights,
								ostream& out);


#endif
