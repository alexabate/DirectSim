#ifndef CHISQSTAT_H_SEEN
#define CHISQSTAT_H_SEEN

#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <typeinfo>

// sophya
#include "machdefs.h"
#include "sopnamsp.h"
#include "fabtwriter.h"
#include "array.h"
#include "hisprof.h"
#include "histerr.h"
#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "stsrand.h"

// DirectSim
#include "sinterp.h"
#include "geneutils.h"
#include "constcosmo.h"

/** @class ChisqStats
  * Simple chisq calculation
  * For a 1D probability function finds the best fit value and the error
  */
class ChisqStats 
{
public:

    /** Constructor 
        @param xvals     parameter values of chi-square (must be evenly spaced)
        @param chisq     chi-square value for each parameter
        @param dof       degree of freedom                                    */
	ChisqStats(TVector<r_8> xvals, TVector<r_8> chisq, double dof=1) {};
	
	/** Find the best fit parameter value ie parameter value at the minimum 
	    chi-square                                                            */
	double BestFit();
	
	/** Estimate the parameter error 
	    @param siglow    lower error range
	    @param sighigh   upper error range
	    @param clevel    confidence level of error range
	    @param npt       number of points to interpolate chi-square func with */
	void ErrSig(double& siglow, double& sighigh, double clevel=0.683,int npt=100); // find error
	//int NearestIndex(double,TVector<r_8>);
	
	/** Return chi-square function                                            */
	TVector<r_8> ReturnChisq(){return Chisq_;};
	
	/** Return parameter values of chi-square                                 */
	TVector<r_8> ReturnXvals(){return xvals_;};
	
	/** Return lower error range                                              */
	double ReturnMinus(){return minus_;};
	
	/** Return upper error range                                              */
	double ReturnPlus(){return plus_;};
	
protected:
	TVector<r_8> xvals_;     /**< parameter values of chi-square              */
	TVector<r_8> Chisq_;     /**< chi-square value for each parameter         */
	double dof_;             /**< degree of freedom                           */
	double minus_;           /**< lower error range                           */
	double plus_;            /**< upper error range                           */

};

#endif


