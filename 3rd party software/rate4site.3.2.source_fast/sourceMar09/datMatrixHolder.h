// $Id: datMatrixHolder.h 678 2006-05-23 11:39:04Z eyalprivman $

#ifndef ___DATMATRIXHOLDER
#define ___DATMATRIXHOLDER

#include <string>
using namespace std;

// THIS CONSTRUCT IS USED TO KEEP A STRING THAT IS THE AA SUBSTITUTION MATRIX
// THE datMatrixString IS TO BE USED WHENEVER WE USE ONE OF THE BUILD-IN AA SUBSTITUTION MATRICES.

class datMatrixString {
public:
  const string Val;
  explicit datMatrixString(const char * str): Val(str){};
};

class datMatrixHolder {
public:
  static const datMatrixString cpREV45;
  static const datMatrixString dayhoff;
  static const datMatrixString jones;	// This is JTT
  static const datMatrixString mtREV24;
  static const datMatrixString wag;
};

#endif	// ___DATMATRIXHOLDER
