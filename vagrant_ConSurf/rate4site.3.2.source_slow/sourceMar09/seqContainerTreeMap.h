// $Id: seqContainerTreeMap.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___SEQUENCE_CONTAINER_TREE_MAP
#define ___SEQUENCE_CONTAINER_TREE_MAP
#include "definitions.h"
#include "tree.h"
#include "treeIt.h"
#include "sequenceContainer.h"

void checkThatNamesInTreeAreSameAsNamesInSequenceContainer(const tree& et,const sequenceContainer & sc);


class seqContainerTreeMap {
public:
	explicit seqContainerTreeMap(const sequenceContainer& sc,
								const tree& et) {
		checkThatNamesInTreeAreSameAsNamesInSequenceContainer(et,sc);
		_V.resize(et.getNodesNum());
		treeIterTopDownConst tit(et);
		for (tree::nodeP myN = tit.first();myN!=tit.end(); myN = tit.next()) {
			if (myN->isInternal()) {
				_V[myN->id()] = -1;
			} else {
				_V[myN->id()] = sc.getId(myN->name(),false);
			}
		}
	}
	int seqIdOfNodeI(const int nodeID) {
		return _V[nodeID];
	}

private:
	vector<int> _V;// _V[i] is the sequenceId of node I.
};

#endif
