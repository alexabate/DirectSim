#ifndef FITBAO_H_SEEN
#define FITBAO_H_SEEN

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

#include "pkspectrum.h"
#include "luc.h"
#include "chisqstats.h"

#define PI 3.141592
#define T_CMB_Par 2.725

//--- Class to compute Decaying sinusoid, See Blake & Glazebrook 2003 eqn 3
// The power of 1.4 in the decay term originates from the Silk damping fitting
// formula in Eisenstein & Hu 1998.  Varying this decay length as well as the
// amplitude will probably not have a significant effect on the fitted values of 
// ka
class DecaySineFunc
{
public:
	DecaySineFunc(double amp, double h, double decay = 1.4) 
    	: amp_(amp) , h_(h) , decay_(decay)
    	{
    	}
 
  	virtual double operator()(double k, double ka)
    		{
    		return ( 1+amp_*k*exp(-  (pow(k,decay_)/pow(0.1*h_,decay_))  )*sin(2*PI*k/ka) );
    		}  

   	void PrintParas()
		{ cout <<"     Amplitude = "<<amp_<<", decay length = "<<decay_<<", h = "<<h_<<endl; };
	

protected:
  double amp_, h_, decay_;
 
};


// Class to compute chisq between wiggles-only power spectrum and decaying sinosoid 
class FitBAOScale
{
public:

	// Constructor
	// Power spectrum with cosmology defined by su at redshift zref
	// is calculated in the constructor
	// Sets default ka range
	// Writes kobs_,Pobs_,sig_
	FitBAOScale(TArray<r_8> pspectrum, SimpleUniverse& su, double zref,double sig8=0.8,double n=1);
	
	// Initialise vectors
	void InitVect(TArray<r_8>);
	// Compute smooth power spectrum
	void ComputeSmoothPS(double OmegaM, double OmegaL, double OmegaB, double h, double sig8=0.8, double n=1, double R = 8);
	// Set ka range
	void Setka(double minka, double maxka, int nka)
		{ minka_=minka; maxka_=maxka; nka_=nka; };
	// Fill Pratio_ vector
	void Pratio();
	// Compute the chisq
	void ComputeChisq(double maxk = 1);
	// Write the chisq as a function of ka to a file
	void WriteChisq(string outfile);
	// Return redshift of power spectrum
	void ReturnZPS(double& zref)
		{ zref=zref_; };
	// Return power spectrum
	void ReturnPS(TVector<r_8>& ko,TVector<r_8>& Po,TVector<r_8>& sig)
		{ ko = kobs_; Po = Pobs_; sig = sig_; };
	// Return power spectrum ratio
	void ReturnPSRatio(TVector<r_8>& Pratio)
		{ Pratio = Pratio_; };
	// Write fiducial/reference power spectrum
	void WriteRefPS(string outfile);
	// Best-fit and N-sigma
	void BestfitStdDev(double& bestfit, double& siglow, double& sighigh, int nsig = 1);
	// Calculate sample variance
	void CalcSigSampVar(double VolCat);
	// Write results to a file
	void WriteResults(string outfile);
	// Write both reference power spectrum AND best-fit sinosoid to a file
	void WriteAncillaryInfo(string outfile);

	/* Class variables*/
	SimpleUniverse& su_;			  // holds cosmological parameters

	TVector<r_8> kobs_, Pobs_, sig_;
	TVector<r_8> Pref_;
	TVector<r_8> Pratio_;
	
	TVector<r_8> kavals_, Chisq_;
	double minka_,maxka_;
	int nka_;
	double zref_;
	double h_;
	double bestfit_,errup_,errdown_;
	int nsig_;
};

#endif


