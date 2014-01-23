/**
 * @file  fitkbaoscale.h
 * @brief Fit BAO scale by fitting decaying sinoid to the "wiggles only"
 *        power spectrum
 *
 * @todo implement more sophisticated fit method
 *
 * @author Alex Abate, Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2008
 * @date 2008
 *
 */
#ifndef FITBAO_H_SEEN
#define FITBAO_H_SEEN


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
#include "pkspectrum.h"
#include "luc.h"
#include "chisqstats.h"
#include "constcosmo.h"


/** @class DecaySineFunc
  *
  * Computes Decaying sinusoid, See Blake & Glazebrook 2003 eqn 3
  *
  * The power of 1.4 in the decay term originates from the Silk damping fitting
  * formula in Eisenstein & Hu 1998.  Varying this decay length as well as the
  * amplitude will probably not have a significant effect on the fitted values of 
  * ka (the BAO scale) 
  */
class DecaySineFunc
{
public:

    /** Constructor
        @param amp    amplitude of decaying sinusoid function 
        @param h      Hubble parameter in units of 100 km/s/Mpc
        @param decay  decay length                                            */
	DecaySineFunc(double amp, double h, double decay = 1.4) 
    	: amp_(amp) , h_(h) , decay_(decay) {};
    	
    /** Return decaying sinusoid function \f$ 1+Ak\exp(-(k/0.1h)^d)\sin(2\pi k/k_a) \f$
        @param k    
        @param ka                                                             */
  	virtual double operator()(double k, double ka) {
    		return ( 1+amp_*k*exp(-  (pow(k,decay_)/pow(0.1*h_,decay_))  )*sin(2*PI*k/ka) );
    		};  

    /** Print parameters of decaying sinusoid function                        */
   	void PrintParas()
		{ cout <<"     Amplitude = "<< amp_ <<", decay length = "<< decay_ <<", h = "<< h_ <<endl; };
	

protected:
    double amp_;        /**< amplitude of decaying sinusoid function            */
    double h_;          /**< Hubble parameter in units of 100 km/s/Mpc          */
    double decay_;      /**< decay length                                       */
 
};


/** @class FitBAOScale
  *
  * Computes chi-square between wiggles-only power spectrum and decaying sinosoid 
  *
  */
class FitBAOScale
{
public:

	/** Constructor
	    @param pspectrum    observed power spectrum: columns(1,2,3) = (k,P(k),sigma_P(k))
	    @param su           cosmology 
	    @param zref         redshift of power spectrum
	    @param sig8         sigma8
	    @param n            spectral index                                    */
	FitBAOScale(TArray<r_8> pspectrum, SimpleUniverse& su, double zref,
	                                   double sig8=0.8, double n=1);
	
	/** Initialise observed power spectrum variables and other vector variable 
	    sizes using pspectrum
	    @param pspectrum    observed power spectrum                           */
	void InitVect(TArray<r_8> pspectrum);
	
	/** Compute smooth (no wiggles) power spectrum
	    @param OmegaM    matter density
	    @param OmegaL    dark energy/cosmological constant density
	    @param OmegaB    baryon density
	    @param h         Hubble parameter in units of 100 km/s/Mpc
	    @param sig8      sigma8 (amplitude of fluctuations)
	    @param n         spectra index
	    @param R         scale of sigma8                                      */
	void ComputeSmoothPS(double OmegaM, double OmegaL, double OmegaB, double h, 
	                                 double sig8=0.8, double n=1, double R = 8);
	
	/** Set grid of trial ka values 
	    minka    start ka value
	    maxka    end ka value
	    nka      number of ka values to try                                   */
	void Setka(double minka, double maxka, int nka)
		{ minka_=minka; maxka_=maxka; nka_=nka; };
		
	/** Compute ratio between observed and reference power spectra            */
	void Pratio() {
	    for (int i=0; i<kobs_.Size(); i++)
		    Pratio_(i) = Pobs_(i)/Pref_(i); };
		    
	/** Compute the chi-square values as a function of ka using ratio of power
	    spectra up to some maximum k
	    @param maxk    maximum k of power spectra to use in chi-square calc   */
	void ComputeChisq(double maxk = 1);
	
	/** Write the chi-square as a function of ka to a file
	    @param outfile    file to write chi-square values to                  */
	void WriteChisq(string outfile);
	
	/** Return redshift of power spectrum                                     */
	double ReturnZPS() { return zref_; };
	
	/** Return observed power spectrum 
	    @param ko    k values
	    @param Po    observed P(k)
	    @param sig   errors on observed P(k)                                  */
	void ReturnPS(TVector<r_8>& ko, TVector<r_8>& Po, TVector<r_8>& sig)
		{ ko = kobs_; Po = Pobs_; sig = sig_; };
		
	/** Return ratio between observed and reference power spectra             */
	TVector<r_8> ReturnPSRatio() { return Pratio_; };
	
	/** Compute best-fit ka with error with given n-sigma precision 
	   @param bestfit    best-fit ka
	   @param siglow     lower error range
	   @param sighigh    upper error range
	   @param nsig       n-sigma precision                                    */
	void BestfitStdDev(double& bestfit, double& siglow, double& sighigh, int nsig = 1);
	
	/** Analytical approximation of the error on the power spectrum due to the 
	    sample variance 
	    @param VolCat    volume of galaxy catalog used to compute power spectrum */
	void CalcSigSampVar(double VolCat);
	
	/** Write results: redshift of power spectrum, best-fit ka, n-sigma of 
	    error ranges, upper error range, lower error range to a file 
	    @param outfile    file to write results to                            */
	void WriteResults(string outfile);
	
	/** Write fiducial/reference power spectrum to a file                 
	    @param outfile    file to write reference power spectrum to           */
	void WriteRefPS(string outfile);
	
	/** Write reference power spectrum AND best-fit sinosoid to a file 
	    @param outfile    file to write functions to                          */
	void WriteAncillaryInfo(string outfile);

protected:
	SimpleUniverse& su_;	 /**< cosmology                                   */
	TVector<r_8> kobs_;      /**< k values of observed power spectrum         */
	TVector<r_8> Pobs_;      /**< observed power spectrum                     */
	TVector<r_8> sig_;       /**< error on observed power spectrum            */
	TVector<r_8> Pref_;      /**< smooth reference power spectrum             */
	TVector<r_8> Pratio_;    /**< ratio between observed and reference power spectra */
	TVector<r_8> kavals_;    /**< grid of trial ka values                     */
	TVector<r_8> Chisq_;     /**< chi-square value at each ka value           */
	double minka_;           /**< min ka value in trial ka grid               */   
	double maxka_;           /**< max ka value in trial ka grid               */ 
	int nka_;                /**< number of ka values in trial ka grid        */
	double zref_;            /**< redshift of power spectrum                  */
	double h_;               /**< Hubble parameter in units of 100 km/s/Mpc   */
	double bestfit_;         /**< best fit ka                                 */
	double errup_;           /**< upper error range                           */
	double errdown_;         /**< lower error range                           */
	int nsig_;               /**< n-sigma of error ranges                     */
};

#endif


