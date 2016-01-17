// $Id: fromCountTableComponentToDistance.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE
#define ___FROM_COUNT_TABLE_COMPONENT_TO_DISTANCE

#include "definitions.h"
#include "countTableComponent.h"
#include "stochasticProcess.h"

static const MDOUBLE startingGuessForTreeBrLen = 0.029;

class fromCountTableComponentToDistance {

public:
	explicit fromCountTableComponentToDistance(
		const countTableComponentGam& ctc,
		const stochasticProcess &sp,
		const MDOUBLE toll,
		const MDOUBLE brLenIntialGuess);// =startingGuessForTreeBrLen

	void computeDistance();// return the likelihood
	MDOUBLE getDistance() { return _distance;} // return the distance.
	MDOUBLE getLikeDistance() { return _likeDistance;} // return the distance.
private:
	const stochasticProcess & _sp;
	const countTableComponentGam& _ctc;
	MDOUBLE _toll;
	MDOUBLE _distance;
	MDOUBLE _likeDistance;
	int alphabetSize() {return _ctc.alphabetSize();}
};

#endif

