// $Id: treeInference.h 391 2005-06-12 12:44:57Z ninio $
// 

// version 1.01
// last modified 23 May 2005

#ifndef ___TREE_INFERENCE
#define ___TREE_INFERENCE

#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "nj.h"
#include <vector>
using namespace std;

class treeInference {
public:
	static tree computeNJtreeWithLikeDist(const stochasticProcess &sp, const sequenceContainer &sc, 
				   const tree * const constraintTreePtr = NULL, const vector<MDOUBLE> * const weights = NULL);

};
#endif


