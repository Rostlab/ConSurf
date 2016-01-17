// $Id: distanceTable.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___DISTANCE_TABLE
#define ___DISTANCE_TABLE

#include "definitions.h"
#include "distanceMethod.h"
#include "sequenceContainer.h"

void giveDistanceTable(const distanceMethod* dis,
					   const sequenceContainer& sc,
					   VVdouble& res,
					   vector<string>& names,
					   const vector<MDOUBLE> * weights = NULL);


#endif
