// $Id: distanceTable.cpp 513 2005-12-19 13:59:07Z eyalprivman $

#include "definitions.h"
#include "distanceTable.h"

void giveDistanceTable(const distanceMethod* dis,
		       const sequenceContainer& sc,
		       VVdouble& res,
		       vector<string>& names,
		       const vector<MDOUBLE> * weights){
	res.resize(sc.numberOfSeqs());
	for (int z=0; z< sc.numberOfSeqs();++z) res[z].resize(sc.numberOfSeqs(),0.0);

	for (int i=0; i < sc.numberOfSeqs();++i) {
		for (int j=i+1; j < sc.numberOfSeqs();++j) {
			res[i][j] = dis->giveDistance(sc[i],sc[j],weights,NULL);
			//LOG(5,<<"res["<<i<<"]["<<j<<"] ="<<res[i][j]<<endl);
		}
		names.push_back(sc[i].name());
	}
}
