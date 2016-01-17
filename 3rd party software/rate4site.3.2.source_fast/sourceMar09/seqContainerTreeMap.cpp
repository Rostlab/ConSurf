// $Id: seqContainerTreeMap.cpp 1213 2006-11-27 15:32:46Z adist $

#include "seqContainerTreeMap.h"
#include "logFile.h"


void checkThatNamesInTreeAreSameAsNamesInSequenceContainer(const tree& et,const sequenceContainer & sc){
	treeIterDownTopConst tIt(et);
	//cout<<"tree names:"<<endl;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isInternal()) 
			continue;

		bool bFound = false;
		sequenceContainer::constTaxaIterator it=sc.constTaxaBegin();
		for (;it != sc.constTaxaEnd(); ++it) 
		{
			if (it->name() == mynode->name()) 
			{
				bFound = true;
				break;
			}
		}
		if (bFound == false) 
		{
			string errMsg = "The sequence name: ";
			errMsg += mynode->name();
			errMsg += " was found in the tree file but not found in the sequence file.\n";
			LOG(4,<<errMsg<<endl);
			errorMsg::reportError(errMsg);
		}
	}
}

