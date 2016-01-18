// $Id: likeDistProp.cpp 391 2005-06-12 12:44:57Z ninio $

#include "likeDistProp.h"
#include "numRec.h"

const MDOUBLE likeDistProp::giveDistance(	const vector<countTableComponentGam>& ctc,
										   MDOUBLE& resL) const {
	const MDOUBLE MAXDISTANCE=2.0;
//	const MDOUBLE PRECISION_TOLL=0.001;
	const MDOUBLE ax=0,bx=1.0,cx=MAXDISTANCE,tol=_toll;
	MDOUBLE dist=-1.0;
	resL = -dbrent(ax,bx,cx,
		  C_evallikeDistProp(ctc,_s1),
		  C_evallikeDistProp_d(ctc,_s1),
		  tol,
		  &dist);
	return dist;
}

// the minus resL = -dbrent because C_evalDist return - value, because it is computing the min not the max...

