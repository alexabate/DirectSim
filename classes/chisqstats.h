#ifndef CHISQSTAT_H_SEEN
#define CHISQSTAT_H_SEEN

#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <typeinfo>
#include "fabtwriter.h"

#include "array.h"
#include "hisprof.h"
#include "histerr.h"

#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "stsrand.h"

#include "sedfilter.h"

#define PI 3.141592

// Simple chisq calculation
// For a 1D probability function finds the best fit value and the error
class ChisqStats 
{
public:
	ChisqStats(TVector<r_8>, TVector<r_8>, double dof=1);
	double BestFit(); // find best fit value
	void ErrSig(double& siglow, double& sighigh, double clevel=0.683,int npt=100); // find error
	int NearestIndex(double,TVector<r_8>);
	
	TVector<r_8> ReturnChisq(){return Chisq_;};
	TVector<r_8> ReturnXvals(){return xvals_;};
	double ReturnMinus(){return minus_;};
	double ReturnPlus(){return plus_;};
	
	/* Class variables*/
	TVector<r_8> xvals_, Chisq_;
	double dof_;
	double minus_, plus_;

};

#endif


