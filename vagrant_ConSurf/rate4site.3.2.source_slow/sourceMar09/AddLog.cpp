// $Id: AddLog.cpp 391 2005-06-12 12:44:57Z ninio $

// version 1.00
// last modified 3 Nov 2002

#include "AddLog.h"
#include <cmath>

const int tAddLog_Precompute::G_LOGADD = 500;
const int tAddLog_Precompute::D_LOGADD = 50;

tAddLog_Precompute AddLogData;

int tAddLog_Precompute::d_logadd;

tAddLog_Precompute::tAddLog_Precompute(){
  d_logadd = int(D_LOGADD*log(10.0)*G_LOGADD);
  logaddf = new double [d_logadd+1];
  for (int i=0; i<= d_logadd; i++)
    logaddf[i] = log(1.0+exp(-static_cast<double>(i)/G_LOGADD));
}

tAddLog_Precompute::~tAddLog_Precompute(){
  delete [] logaddf;
}
