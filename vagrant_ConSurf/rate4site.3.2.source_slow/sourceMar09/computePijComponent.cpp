
// $Id: computePijComponent.cpp 551 2006-02-08 17:16:59Z ninio $

#include "definitions.h"
#include "treeIt.h"
#include "computePijComponent.h"
#include "logFile.h"

void computePijHomSpec::fillPij(const MDOUBLE dis,
								const stochasticProcess& sp,
								int derivationOrder) {
  resize(sp.alphabetSize());
  int i,j;
  for (i=0; i<sp.alphabetSize(); i++) {
	switch (derivationOrder) {
	case 0:
	  _V[i][i] = sp.Pij_t(i,i,dis);
	  break;
	case 1:
	  _V[i][i] = sp.dPij_dt(i,i,dis);
	  break;
	case 2:
	  _V[i][i] = sp.d2Pij_dt2(i,i,dis);
	  break;
	default:
	  errorMsg::reportError("error in function fillPij - derivationOrder must be 0, 1 or 2");		
	}

	  
	for (j=i+1; j<sp.alphabetSize(); j++) {
	  switch (derivationOrder) {
	  case 0:
		_V[i][j] = sp.Pij_t(i,j,dis);
		break;
	  case 1:
		_V[i][j] = sp.dPij_dt(i,j,dis);
		break;
	  case 2:
		_V[i][j] = sp.d2Pij_dt2(i,j,dis);
		break;
	  default:
		errorMsg::reportError("error in function fillPij - derivationOrder must be 0, 1 or 2");		
	  }
	  if (sp.freq(j) == 0.0) {
		//if (_V[i][j]!=0.0) {
		errorMsg::reportError("error in function fillPij");
		//}
		//_V[j][i] = 0.0;
	  }
	  else {
		_V[j][i] = _V[i][j]*
		  sp.freq(i)/sp.freq(j);
	  }
	}
  }
}


void computePijHom::fillPij(const tree& et, const stochasticProcess& sp, int derivationOrder) {
	_V.resize(et.getNodesNum());
	treeIterTopDownConst tIt(et);
	tree::nodeP myNode = tIt.first();
	{// skipping the root, but allocating place for the root pij even if they are not use
	 // to maintain that all arrays have the same size.
		_V[myNode->id()].resize(sp.alphabetSize());
	}
	LOGDO(50,et.output(myLog::LogFile(),tree::ANCESTOR));
	LOGDO(50,et.output(myLog::LogFile(),tree::PHYLIP));
	for (; myNode != tIt.end(); myNode = tIt.next()) {
	  if (!(myNode->isRoot()))
		  _V[myNode->id()].fillPij(myNode->dis2father()*sp.getGlobalRate(),sp,derivationOrder);
//	  else
//	    myLog::LogFile()<<"ROOT IS "<<myNode->name()<<endl;
	}
}


void computePijGam::fillPij(const tree& et, const stochasticProcess& sp, int derivationOrder) {
	_V.resize(sp.categories());
	for (int i=0; i < _V.size(); ++i) {
		tree cp = et;
		cp.multipleAllBranchesByFactor(sp.rates(i)/sp.getGlobalRate());// the global rate is taken care of in the hom pij.
		_V[i].fillPij(cp,sp,derivationOrder);
	}
}
