// $Id: recognizeFormat.h 391 2005-06-12 12:44:57Z ninio $

#ifndef ___RECOGNIZE_FORMAT
#define ___RECOGNIZE_FORMAT

#include "sequenceContainer.h"

class recognizeFormat{
public:
	static sequenceContainer read(istream &infile, const alphabet* alph);
	static void write(ostream &out, const sequenceContainer& sd);
	//readUnAligned: the input sequences do not need to be aligned (not all sequences are the same length).
	static sequenceContainer readUnAligned(istream &infile, const alphabet* alph);
};

#endif



