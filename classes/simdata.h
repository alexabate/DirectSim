/**
 * @file  simdata.h
 * @brief Simulates galaxy photometry
 *
 * Could add more information here I think
 *
 * @author Alex Abate, Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2009
 * @date 2009
 *
 */
 
#ifndef  DATA_H_SEEN
#define  DATA_H_SEEN

#include <iostream>
#include <fstream>
#include <math.h>

#include "array.h"
// Sophya update v2.3 June 2013 replaces genericfunc with classfunc
//#include "genericfunc.h"
#include "classfunc.h"
#include "pexceptions.h"
#include "stsrand.h"
#include "sopnamsp.h"
#include "integ.h"
#include "fitsioserver.h"
#include "swfitsdtable.h"
#include "tarray.h"
#include "igm.h"

#include "sinterp.h"
#include "gftdist.h"
#include "cosmocalcs.h"
#include "sedfilter.h"

/** @class
  * PhotometryCalcs class
  *
  * Class holding methods that perform calculations related to photometric
  * observations
  *
  * @note a ``maggie" is not exactly flux, it is a linear measure of flux
  * whereas magnitudes are a logarithmic measure of flux.
  */
class PhotometryCalcs {
public:
    
    /** Constructor  */
    PhotometryCalcs(double lmin = 5e-8, double lmax = 2.5e-6, int npt=10000) { 
        setLminmax(lmin, lmax, npt); };
        
    /** Set the wavelength range and resolution to do the integrals over */
    void setLminmax(double lmin, double lmax, int npt)
        { lmin_ = lmin; lmax_ = lmax; npt_ = npt; };
        
        
    // Photometry calculations //
    
    /** The k-correction is the relation between the rest-frame (ie emitted frame) 
        absolute magnitude of a source in one bandpass \f$X\f$ to the observed-frame 
        apparent magnitude of the same source in another bandpass \f$Y\f$:
        \f$ K_{xy} = -2.5\log10\left( \frac{1}{1+z}
        \frac{\int f_\lambda(\lambda_o/(1+z))\lambda_oX(\lambda_o)d\lambda_o
              \int \frac{Y(\lambda_o)}{\lambda_o}d\lambda_o  }
             {\int \frac{X(\lambda_o)}{\lambda_o}d\lambda_o
              \int f_\lambda(\lambda_e)\lambda_oY(\lambda_o)d\lambda_o}\right) \f$
        @param z                redshift of object \f$z\f$
        @param sed              rest-frame SED of object \f$f_\lambda(\lambda)\f$
        @param filterX          filter object observed in \f$X(\lambda)\f$
        @param restFrameFilter  rest-frame filter \f$Y(\lambda)\f$            */
	double Kcorr(double z, SpecEnergyDist& sed, Filter& filterX, Filter& restFrameFilter);
	
	/** The k-correction calculation when the rest-frame and observed-frame bandpasses
	    are the same:  \f$ K_{xy} = -2.5\log10\left( \frac{1}{1+z} 
	    \frac{\int f_\lambda(\lambda_o/(1+z))\lambda_oX(\lambda_o)d\lambda_o}
	         {\int f_\lambda(\lambda_e)\lambda_oX(\lambda_o)d\lambda_o} \right) \f$
	    @param z                redshift of object \f$z\f$
        @param sed              rest-frame SED of object \f$f_\lambda(\lambda)\f$
        @param filterX          filter object observed in (same as rest-frame filter) \f$X(\lambda)\f$*/
	double Kcorr1Filter(double z, SpecEnergyDist& sed, Filter& filterX);

	/** Calculate galaxy color X-Y: 
	    \f$ C_{xy} = -2.5\log10\left(
        \frac{\int f_\lambda(\lambda_o/(1+z))\lambda_oX(\lambda_o)d\lambda_o
              \int \frac{Y(\lambda_o)}{\lambda_o}d\lambda_o  }
             {\int \frac{X(\lambda_o)}{\lambda_o}d\lambda_o
              \int f_\lambda(\lambda_o/(1+z))\lambda_oY(\lambda_o)d\lambda_o}\right) \f$
	    @param z                redshift of object \f$z\f$
        @param sed              rest-frame SED of object \f$f_\lambda(\lambda)\f$
        @param filterX          filter object observed in \f$X(\lambda)\f$ (should be bluer filter than Y)
        @param filterY          other filter object observed in \f$Y(\lambda)\f$  */
	double CompColor(double z, SpecEnergyDist& sed, Filter& filterX, Filter& filterY);
	
	// AA: not 100% what these two methods, restFrameFlux and restFrameFluxLambda do?
	
	/** Calculate rest-frame flux of object in band \f$X(\lambda)\f$ in FREQUENCY units: 
        \f$ F_\nu(\lambda^{eff}_e) = 
            \frac{\int f_\lambda(\lambda_e)X(\lambda_e)\lambda_e d\lambda_e}
                 {\int X(\lambda_e)/\lambda_e d\lambda_e} \f$  
        @param sed      rest-frame SED of object \f$f_\lambda(\lambda)\f$
        @param filter   filter object observed in \f$X(\lambda)\f$
        @param zs       redshift of object                                    */
	double restFrameFlux(SpecEnergyDist& sed, Filter& filter, double zs);
	            
	/** Rest-frame flux in WAVELENGTH units: 
	    \f$ F_\lambda(\lambda^{eff}_e) = \frac{F_\nu(\lambda^{eff}_e)}{\lambda_e^2} \f$*/
	double restFrameFluxLambda(SpecEnergyDist& sed, Filter& filter, double zs) {
	    double f0nu = restFrameFlux(sed, filter, zs);
	    BlueShiftFilter blueshiftFilter(filter, zs);
	    double lamEffRF = effectiveFilterWavelength(blueshiftFilter);
	    double f0l = f0nu/(lamEffRF*lamEffRF);
	    return f0l;
	    };
	    
	
	// Magnitude < - > Flux calculations //
	
	/** Return flux in FREQUENCY units given AB magnitude 
	    @note Make sure the flux you want is compatible with the AB magnitude system!
	    @param mag          AB magnitude of object
	    @param zs           redshift of object
	    @param dL           luminosity distance at zs
	    @param sed          SED of object
	    @param filterX      arbitrary filter band                             */
	double convertABMagToFlux(double mag, double zs, double dL, SpecEnergyDist& sed, Filter& filterX) {
        double magPart = pow(10,-0.4*mag);
        SEDzFilterProd sedXlambdaXfilter(sed, filterX, 0);// returns sed*lambda*filter
        FilterIntegrator integrandSED(sedXlambdaXfilter, lmin_, lmax_, npt_);
        double fluxPart = integrandSED.Value();
        double zeroPointX = getFilterZeroPointFlux(filterX);
        double fnu = magPart*dL*dL*(1.+zs)*(fluxPart/zeroPointX);
	    return fnu; };
	    
	/** Return flux in WAVELENGTH units given AB magnitude 
	    @note Make sure the flux you want is compatible with the AB magnitude system!
	    @param mag          AB magnitude of object
	    @param zs           redshift of object
	    @param dL           luminosity distance at zs
	    @param sed          SED of object
	    @param filterX      arbitrary filter band
	    @param filterY      filter band object observed in                           */
	double convertABMagToFluxLambda(double mag, double zs, double dL, SpecEnergyDist& sed, Filter& filterX, Filter& filterY) {
	    double fnu = convertABMagToFlux(mag, zs, dL, sed, filterX);
	    double lambdaEff = effectiveFilterWavelength(filterY);
        double fl = fnu/(lambdaEff*lambdaEff);
	    return fl; };
	    
	/** Return AB magnitude given maggies in FREQUENCY units:
	    \f$ m_{AB} = -2.5\log10(F_\nu) - 56.1 + zp_X \f$
	    @note Make sure the flux you want is compatible with the AB magnitude system!
	    @param flux         flux (in maggies, defined in AB system) \f$F_\nu\f$
	    @param filterX      filter band \f$X(\lambda)\f$                      */
	double convertFluxMaggiesToABMag(double flux, Filter& filterX) {
	    double zeroPoint = getFilterZeroPoint(filterX);

	    if (flux<1e-50)
	        flux = 1e-50;
	    double magAB = -2.5*log10(flux) - 56.1 + zeroPoint;
	    return magAB; 
	    };
	    
	/** Return AB magnitude given flux in maggie-style WAVELENGTH units:
	    \f$ F_\nu = F_\lambda*\lambda_{eff}^2
	        m_{AB} = -2.5\log10(F_\nu) - 56.1 + zp_X \f$
	    @note Make sure the flux you want is compatible with the AB magnitude system!
	    @param flux         flux (defined in AB system but in wavelength units) \f$ F_\lambda \f$
	    @param filterX      filter band \f$X(\lambda)\f$                      */
	double convertFluxMaggiesLambdaToABMag(double flux, Filter& filterX) {
	    double lambdaEff = effectiveFilterWavelength(filterX);
	    double fnu = lambdaEff*lambdaEff*flux;
	    double mAB = convertFluxMaggiesToABMag(fnu, filterX);
	    return mAB; };
	
	/** Return maggies in FREQUENCY units given AB magnitude
	    \f$ maggie = 10^{-0.4(m+56.1+zp_x)} \f$ */
	double convertABMagToFluxMaggies(double mag,  Filter& filterX) {
	    double zeroPoint = getFilterZeroPointFlux(filterX);
	    double flux = pow(10.,-0.4*(mag+56.1))*zeroPoint;
	    return flux;
	    }


	    
	    
    // Magnitude ERROR < - > Flux ERROR calculations //
        
    /** Convert magnitude error to a flux error: \f$ \sigma_F = 0.4\sigma_m F \ln(10) $\f
        @note units of flux error returned depend on units of flux provided
        @param sigmaM   magnitude error \f$\sigma_m \f$
        @param flux     flux of object \f$F\f$                                */
    double convertMagErrorToFluxError(double sigmaM, double flux) { 
	            double fluxError = 0.4*sigmaM*flux*log(10);   
	            return fluxError;
	            };
     
     /** Convert flux error to a magnitude error: \f$ \sigma_m = 2.5/\ln(10) \sigma_F/F \f$
        @note if AB magnitude system flux must be in FREQUENCY units
        @param sigmaFoverF   flux error \f$ \sigma_F/F \f$                    */
    double convertFluxErrorToMagError(double sigmaFoverF) { 
	            double magError = (2.5/log(10))*(sigmaFoverF);   
	            return magError;
	            };
	            	
	
	       
	// Filter calculations //
	
	/** Get "zeropoint" of filter, basically integrates dnu/nu filter(nu). 
	    Returns the zeropoint in magnitudes FREQUENCY UNITS (not flux)        */
	double getFilterZeroPoint(Filter& filterX)
	    {   double fx = getFilterZeroPointFlux(filterX);
	        double mx = 2.5*log10(fx);
	        return mx; };
	// for LSST: -3.11752 -1.65737 -1.90181 -2.06726 -2.66691 -4.33277

	/** Get "zeropoint" of filter, basically integrates dnu/nu filter(nu). 
	    Returns the zeropoint in flux FREQUENCY UNITS (not magnitudes)        */
	double getFilterZeroPointFlux(Filter& filterX);
	
	/** Calculate effective wavelength of filter.  Calculates: 
	    \f$ \frac{\int X(\lambda)\lambda d\lambda}{\int X(\lambda) d\lambda}\f$
	    where the integrals are between the filter lower and upper edges.*/
	double effectiveFilterWavelength(Filter& filterX);
	
	/** Return the value of the maximum transmission of the filter 
	    @param  filterX     filter transmission function
	    @param  lambdaAtMax wavelength at maximum transmission (returned by method)
	    @param  nStep       number of steps in the search between #lmin_, #lmax_*/
	double findFilterMax(Filter& filterX, double& lambdaAtMax, int nStep=1000);
	
	/** Find the wavelength closest to the transmission value #trans in filter
	    #iFilter between lmin and lmax
	    @param  filterX     filter transmission function
	    @param  trans       transmission value to find wavelength of 
	    @param  lmin        search for closest wavelength starting at #lmin
	    @param  lmax        search for closest wavelength ending at #lmax
	    @param  nStep       number of steps in the search between #lmin
	                        and #lmax */
	double findFilterTransValue(Filter& filterX, double trans, double lmin, double lmax, int nStep=1000);
	
	/** Find wavelength edges of filter filterX 
	    @param lmin             lower wavelength edge of filter (returned by method)
	    @param lmax             upper wavelength edge of filter (returned by method)
	    @param filterX          filter transmission function
	    @param edgeDefinition   percent of filter maximum defined as the "edge"*/
	void findFilterEdges(double& lmin, double& lmax, Filter& filterX, double edgeDefinition=0.05);
	                                            
	/** Return rest frame wavelength: \f$ \lambda_e = \frac{\lambda_o}{1+z} \f$
	    @param  lambdaObs   observed wavelength \f$ \lambda_o \f$
	    @param  zs          redshift \f$ z \f$*/                               
    double returnRestFrameWaveLength(double lambdaObs, double zs) {
            double lambdaRF = lambdaObs/(1 + zs);
            return lambdaRF;
            };
    
protected:
    double lmin_;  /**< minimum wavelength of range to do integrals over */
    double lmax_;  /**< minimum wavelength of range to do integrals over */
    int npt_;      /**< resolution of integration grid                   */
    //SimpleUniverse su_;
};


/** @class
  * SimData class
  * 
  * Class to simulate/calculate galaxy magnitudes (no photometric errors)
  *
  *
  */
class SimData : public PhotometryCalcs {
public:

    //typedef enum{None=0, Card=1, Calz=2} dustlaw;
    //typedef enum{None=0, Madau=1, Mean=2} igmmodel;
	enum dustLaw{NoDust=0, Card=1, Calz=2};
	enum igmModel{None=0, Madau=1, Mean=2};
	
	                                    
	/** Constructor: calculate magnitudes from scratch
	    @param sedArray     array holding pointers to SED objects 
	    @param filterArray  array holding pointers to Filter objects
	    @param su           object holding cosmological parameters and calculations                         
	    @param lmin         minimum wavelength
	    @param lmax         maximum wavelength
	    @param npt          wavelength resolution                                                           */
	SimData(vector<SED*> sedArray, vector<Filter*> filterArray, SimpleUniverse& su, double lmin = 5e-8, 
	        double lmax = 2.5e-6, int npt=10000)
	: sedArray_(sedArray) , filterArray_(filterArray) , su_(su) , PhotometryCalcs(lmin, lmax, npt) {
              
              nsed_ = sedArray_.size();
              nFilters_ = filterArray_.size();	
              isReadKcorr_ = false;
              //if (nsed_>=1000)
              //    throw ParmError("ERROR! Too many SEDs");
              };


	/** Constructor: read k-corrections from files and interpolate
	    @param sedNames     list of names of all SEDs
	    @param filterSet    name of filter set
	    @param filterArray  array holding pointers to Filter objects
	    @param su           object holding cosmological parameters and calculations
	    @param lmin         minimum wavelength
	    @param lmax         maximum wavelength
	    @param npt          wavelength resolution                                                           */
    SimData(vector<string> sedNames, string filterSet, vector<Filter*> filterArray, SimpleUniverse& su, 
            double lmin = 5e-8, double lmax = 2.5e-6, int npt=10000)
	: sedNames_(sedNames), filterSet_(filterSet) , filterArray_(filterArray) , su_(su) , 
	  PhotometryCalcs(lmin, lmax, npt) { isReadKcorr_=true; };
	                            
	                                    
	/** Destructor */
	virtual ~SimData(void) { };
	
	// Put at top: everything that definitely should be accessible from 
	// OUTSIDE the class
	
	// SIMULATION FUNCTIONS //

    /** Calculate magnitude of a galaxy (rest-frame filter *is* one of the filters in filterArray_)
        @param zs           redshift of galaxy
        @param absmag       absolute magnitude of galaxy in filter indexed by irestfilt
        @param sedid        id of galaxy SED in sedArray_
        @param iobsfilt     id of filter in filterArray_ to calculate magnitude in 
        @param irestfilt    id of filter in filterArray_ absmag is defined in
        @param igmtrans     IGM transmission object (returns transmission by IGM given galaxy redshift z)
        @param ext          internal dust extinction E(B-V) value
        @param law          dust law describing internal dust extinction (Cardelli or Calzetti)             */
    double getMag(double zs, double absmag, int sedid, int iobsfilt, int irestfilt, IGMTransmission& igmtrans,
                  double ext=0., dustLaw law=Card);
                  
    /** Calculate magnitude of a galaxy (rest-frame filter *is not* one of the filters in filterArray_)
        @param zs           redshift of galaxy
        @param absmag       absolute magnitude of galaxy in filter restfilter
        @param sedid        id of galaxy SED in sedArray_
        @param iobsfilt     id of filter in filterArray_ to calculate magnitude in 
        @param restfilter   filter absmag is defined in
        @param igmtrans     IGM transmission object (returns transmission by IGM given galaxy redshift z)
        @param ext          internal dust extinction E(B-V) value
        @param dustlaw      dust law describing internal dust extinction (Cardelli or Calzetti)             */
    double getMag(double zs, double absmag, int sedid, int iobsfilt, Filter& restfilter, IGMTransmission& igmtrans,
                  double ext=0., dustLaw law=Card);
                  
                  
    /** Calculate magnitude of a galaxy using pre-calculated k-correction tables 
        @param zs           redshift of galaxy
        @param absmag       absolute magnitude of galaxy in filter restfilter
        @param sedid        id of galaxy SED in sedArray_
        @param iobsfilt     column id of filter in k-correction table to calculate magnitude in 
        @param irestfilt    column id of filter in k-correction table absmag is defined in
        @param igm          IGM transmission model (could be None, Madau, Mean)
        @param ext          internal dust extinction E(B-V) value
        @param dustlaw      dust law describing internal dust extinction (Cardelli or Calzetti)            */
    double getMag(double zs, double absmag, int sedid, int iobsfilt, int irestfilt, igmModel igm, 
                  double ext, dustLaw law);
	
	
	// SETTINGS FUNCTIONS //
	
	/** Set min and max and wavelength */
	void setXminXmax(double lmin, double lmax)
			{ lmin_=lmin; lmax_=lmax; };

    /** Return fluxes of each SED in each of the rest-frame filters
        @param zs               redshift of the SED
        @param restFrameFilter  filter absolute magnitude is defined in */
    TArray<double> returnSEDFluxesInRestFrame(double zs);//, Filter& restFrameFilter);
	
	
// internal methods
protected:
	
	/** Calculate k correction (from scratch) 
	    @param sed                SED (z=0)
	    @param filter             observation filter
	    @param restFrameFilter    rest-frame filter
	    @param zs                 redshift of SED
	    @param igmtrans           transmission by IGM
	    @param ext                internal dust extinction E(B-V) value
        @param dustlaw            dust law describing internal dust extinction (Cardelli or Calzetti)       */
	double calcKcorr(SED& sed, Filter& filter, Filter& restFrameFilter, double zs, 
	                 IGMTransmission& igmtrans, double ext, dustLaw law);
	                 
	                 
	/** Interpolate k-correction from data read in from a file
	    @param filename    name of file to read in
	    @param zs          redshift of galaxy
	    @param iobsfilt    column id of filter in k-correction table to calculate magnitude in
	    @param irestfilt    column id of filter in k-correction table absmag is defined in                  */
	double interpKcorr(string filename, double zs, int iobsfilt, int irestfilt);
	
	
	/** Interpolate k correction  
	    @param sedID        SED id
	    @param iFilterObs   filter id
	    @param zs           redshift
	    @param ext          extinction                                        */
	//double interpKcorr(int sedID, int iFilterObs, double zs, double ext);
	/** Interpolate k correction  
	    @param linearIndex  SED,filter combination id
	    @param zs           redshift
	    @param ext          extinction                                        */
	//double interpKcorr(int linearIndex, double zs, double ext);
	
	/** Given a row index @param i and a column index @param j, return the single array 
	    element number if there are @param nj columns                         */
	//int returnLinearIndex(int i, int j, int nj) { return i*nj + j; };
	
     
    
protected:

	SimpleUniverse& su_;              /**< class that holds the cosmological parameters and calculations    */
    vector<SED*> sedArray_;           /**< holds the SEDs (for full calculation)                            */
    vector<Filter*> filterArray_;     /**< holds the filters (for full calculation)                         */
    vector<string> sedNames_;         /**< holds the names of all SEDs (for reading k-correction)           */
    string filterSet_;                /**< name of filter set (for reading k-correction)                    */
    bool isReadKcorr_;                /**< read k corrections from file                                     */
	int nsed_;                        /**< number of SEDs                                                   */
	int nFilters_;                    /**< number of filters                                                */
	vector<SInterp2D*> kInterpZExt_;  /**< array of pointers to k-corr interpolation                        */
	// To initialize when the references are not used                                                       */
	SimpleUniverse su_default_;       /**< to initialize #su_ when it's not used                            */
	
	
	
};


/** @class SimObservations
  *
  * Take truth photometric data and add photometric errors.
  *
  */
class SimObservations : public PhotometryCalcs {
public:

    /** Constructor 
        @param filterArray  array holding pointers to Filter objects
	    @param rg           random number generator object                                                  */
    SimObservations(vector<Filter*> filterArray, RandomGeneratorInterface& rg) 
    : filterArray_(filterArray) , rg_(rg) { setLSSTPars(); nFilters_ = filterArray_.size(); };
    
    
    /** Add GENERIC flux percentage error to magnitude. Return the observed magnitude and magnitude error in a
	    vector
	    @param mag             (true) magnitude 
	    @param percentError    error on flux in percent of flux (eg 10% is percentError=0.1) 
	    @param iFilter         index of filter (0=u, 1=g .... 5=y)            */
	vector<double> addError(double mag, double percentError, int iFilter) {
        
        Filter filter((*filterArray_[iFilter]));
        double fluxError = percentError*convertABMagToFluxMaggies(mag, filter);
        double flux = convertABMagToFluxMaggies(mag, filter);
        
        vector<double> observation; // observed magnitude and magnitude error
        observation = addFluxError(flux, fluxError, iFilter);

        double obs_err = abs(mag-observation[0]);
        observation.push_back(obs_err);
	    return observation;
        
        };
	
	
	/** Add LSST photometric error. Returns the observed magnitude and magnitude 
	    error in a vector
	    @param mag      (true) magnitude
	    @param nYear    number of LSST operation       
	    @param iF       index of filter (0=u, 1=g .... 5=y)                   */
	vector<double> addLSSTError(double mag, int nYear, int iF) {
	
	     if (iF<0 || iF>5)
	         throw ParmError("ERROR! filter index not in LSST range");
         //double m5 = calcPointSource5sigmaDepth(Cm_[iF], Msky_[iF], Theta_[iF] ,tVis_, km_[iF], airMass_);
         
         // total number of visits after nYear years
         int nvis = nYear*nVisYear_[iF];
         
         // calculate 5-sigma depth for nvis visits
         double m5nvisit = m5single_[iF] + 1.25*log10(nYear*nVisYear_[iF]);
         
         // BUT to get 1-sigma errors need to subtract -2.5*log10(5) from this
         //double m1sig = m5nvisit - 2.5*log10(5.);
         
         vector<double> obsmag = getObservedLSSTMagnitude(mag, m5nvisit, Gamma_[iF], nYear, iF);
	     return obsmag;
	     
	     };

    /** Return observed LSST magnitude and magnitude error */
    vector<double> getObservedLSSTMagnitude(double mag, double m5, double gamma, int nYear, int iFilter);
	
	
	/** Return observed magnitude and magnitude error after adding flux error 
	    in the filter indexed by #iFilter
	    @param mag          true flux
	    @param fluxError    error on flux
	    @param iFilter      index of filter  */
	vector<double> addFluxError(double flux, double fluxError, int iFilter);
	
	
	/** Return the LSST random photometric error squared
	    See equation 3.2 in LSST Science Book (divided by Nvisit)
	    @param x        \f$x=10^{0.4(m-m5)}\f$  
	    @param gamma    band-dependent parameter */
	    //@param nVis     number of visits */
	double returnLSSTRandomErrorSq(double x, double gamma) {//, double nVis) {
	
	    double sigmaRandsq = (0.04 - gamma)*x + gamma*x*x;
        // in magnitudes^2
        return sigmaRandsq; ///nVis;
        };
        
		       
    /** Return \f$10^{0.4(m-m5)}\f$ */
	double returnX(double mag, double m5) {
	    double x = pow(10.,0.4*(mag-m5));
        return x; 
        };
	
	/** Return the 5-sigma depth for point sources (within particular band) 
	    @param Cm     band dependent parameter describing system sensitivity
	    @param msky   sky brightness
	    @param theta  zenith seeing
	    @param nvis   number of visits
	    @param km     atmospheric extinction
	    @param X      airmass                                                   */
	double calcPointSource5sigmaDepth(double Cm, double msky, double theta, int nvis, double km, double X) {
	    double m5 = Cm + 0.50*(msky - 21) + 2.5*log10(0.7/theta) + 1.25*log10(nvis) - km*(X - 1);
        return m5;
        };    
        
	   
    /** Set the LSST photometric error parameters */
	void setLSSTPars();
	
                        
    /** Return effective filter restframe wavelengths if galaxy is at zs.  Basically
        returns \f$\lambda_{eff}^o = \lambda_{eff}/(1+z)\f$ for each filter   */
    vector<double> returnFilterRFWavelengths(double zs);
    
    
    
    
    
    vector<double> returnm5single() { return m5single_; };
    
protected:
    vector<Filter*> filterArray_;   /**< holds the filters                                                  */
    RandomGeneratorInterface& rg_;  /**< class that generates the random numbers                            */
    DR48RandGen rg_default_;        /**< to initialize #rg_ when it's not used                              */
    int nFilters_;                  /**< number of filters                                                  */
	vector<int> nVisYear_;    /**< number of visits per year                                     */
	double tVis_;             /**< exposure time, 2 back-to-back 15s exposures                   */
	double airMass_;          /**< median airmass                                                */
	double sigmaSys_;         /**< systematic photometric error from LSST system                 */
	vector<double> Msky_;     /**< expected median sky zenith brightness in each band            */
	vector<double> Theta_;    /**< expected delivered median zenith seeing (arcsec) in each band */
	vector<double> Gamma_;    /**< band dependent parameter                                      */
	vector<double> Cm_;       /**< band dependent parameter describing system sensitivity        */
	vector<double> km_;       /**< adopted atmospheric extinction in each band                   */
	vector<double> m5single_; /**< SINGLE VISIT 5-sigma depth for point sources                  */
	

};

/** @class BpzCorrectMags
  *
  * Correct magnitudes so flux ratios match BPZ flux ratios exactly
  *
  */
class BpzCorrectMags {
public:
    /** Constructor
        @param abFileLoc    location of BPZ .AB files
        @param filterNames  string vector containing the filenames (without extensions) of the filters      */
    BpzCorrectMags(string abFileLoc,  vector<string> filterNames)
    : abFileLoc_(abFileLoc) , filterNames_(filterNames) {
        nFilt_ = filterNames_.size();
        };
    
    
    /** Correct magnitudes, returns magnitudes corrected to match BPZ flux ratios exactly 
        @param mags          magnitude in each filter
        @param sedFileName   filename of galaxy SED (without extension)
        @param z             redshift of galaxy                                                             */
    vector<double> correct(vector<double> mags, string sedFileName, double z) {
    
        if (mags.size()!=nFilt_)
            throw ParmError("ERROR! number of magnitudes does not match number of filters");
    
        vector<double> fluxes = convertToFlux(mags);
        //cout <<"Fluxes: ";
        //for (int i=0; i<fluxes.size(); i++)
        //    cout << fluxes[i] <<"  ";
        //cout << endl;
        
        vector<double> bpzFluxes = getBpzFluxes(sedFileName, z);
        //cout <<"BPZ Fluxes: ";
        //for (int i=0; i<bpzFluxes.size(); i++)
        //    cout << bpzFluxes[i] <<"  ";
        //cout << endl;
        
        vector<double> conversions = computeConversion(fluxes, bpzFluxes);
        //cout <<"Conversions: ";
        //for (int i=0; i<conversions.size(); i++)
        //    cout << conversions[i] <<"  ";
        //cout << endl;
        
        return calcNewMags(conversions, fluxes);
         
        };
        
    /** Check magnitudes
        @param mags          (supposedly corrected) magnitude in each filter
        @param sedFileName   filename of galaxy SED (without extension)
        @param z             redshift of galaxy                                                             */
    void check(vector<double> mags, string sedFileName, double z) {
    
        vector<double> fluxes = convertToFlux(mags);
        vector<double> bpzFluxes = getBpzFluxes(sedFileName, z);

        cout <<"     Check bpz flux color divided by simulated flux color is equal to 1 in all filters:";
        for (int i=0; i<nFilt_-1; i++) {
            double bpzFluxColor = bpzFluxes[i]/bpzFluxes[i+1];
            double fluxColor = fluxes[i]/fluxes[i+1];
            cout << bpzFluxColor/fluxColor <<"  ";
            }
        cout << endl;
        }
    
    
    
protected:

    /** Convert magnitudes to flux via \f$ f=10^{-0.4m} \f$
        @param mags    magnitude in each filter                                                             */
    vector<double> convertToFlux(vector<double> mags) {
        vector<double> fluxes;
        for (int i=0; i<nFilt_; i++)
            fluxes.push_back(pow(10.,-0.4*mags[i]));
        return fluxes;
        };
        
    /** Get the BPZ fluxes of galaxy SED at redshift z
        @param sedFileName   filename of galaxy SED (without extension)
        @param z             redshift of galaxy                                                             */
    vector<double> getBpzFluxes(string sedFileName, double z) {
    
        vector<double> bpzFluxes;
        for (int i=0; i<nFilt_; i++) {
            string filename = abFileLoc_ + sedFileName + filterNames_[i] + ".AB";
            SInterp1D fz;
            fz.ReadXYFromFile(filename, 0., 12., 1024, 0, false); // 0<z<12 is BPZ flux file range
            bpzFluxes.push_back(fz(z));
            }
        return bpzFluxes;
        
        };
        
    /** Compute conversion values to make fluxes match BPZ flux ratios exactly 
        @param fluxes     flux in each filter
        @param bpzFluxes  BPZ flux in each filter                                                           */
    vector<double> computeConversion(vector<double> fluxes, vector<double> bpzFluxes) {
        vector<double> conversions;
        conversions.push_back(1.); // arbitrarily set conversion of first filter flux to 1.
        
        //cout <<"Flux colors (bpz, mine, ratio, conver): ";
        for (int i=0; i<nFilt_-1; i++) {
            double bpzFluxColor = bpzFluxes[i]/bpzFluxes[i+1];
            double fluxColor = fluxes[i]/fluxes[i+1];
            //cout << bpzFluxColor <<","<< fluxColor <<"  ";
            double ratio = bpzFluxColor/fluxColor;
            //cout << ratio <<"  "<< conversions[i]*ratio <<"  ";
            conversions.push_back( conversions[i]*(1./ratio) );
            }
        //cout << endl;
        return conversions;
        };
        
    /** Return the corrected magnitudes
        @param conversions  conversion values to make fluxes match BPZ flux ratios exactly
        @param fluxes       flux in each filter                                                             */
    vector<double> calcNewMags(vector<double> conversions, vector<double> fluxes) {
        vector<double> newmags;
        for (int i=0; i<nFilt_; i++)
            newmags.push_back(-2.5*log10(conversions[i]*fluxes[i]));
        return newmags;
    
        };

protected:
    string abFileLoc_;            /**< location of BPZ .AB files                     */
    vector<string> filterNames_;  /**< filenames (without extensions) of the filters */
    int nFilt_;                   /**< number of filters                             */

};

/** @class SEDLibColors
  *
  * Calculate colors of all SEDs in a library for a given filter set
  *
  */
class SEDLibColors : public PhotometryCalcs {
public:
    /** Constructor. The filter set must be ordered by monotonic increase in filter wavelength from blue to red
        @param sedarray   array of pointers to each SED in library
        @param filtarray  array of pointers to each filter in filter set 
        @param lmin    minimum wavelength of SED in meters
        @param lmax    maximum wavelength of SED in meters 
        @param npt     number of interpolation points for SED func                                          */
    SEDLibColors(vector<SED*> sedarray,  vector<Filter*> filtarray, double lmin=5e-8, double lmax=2.5e-6, int npt=10000):
        sedarray_(sedarray), filtarray_(filtarray) , nsed_(sedarray.size()) , nfilt_(filtarray.size()) {
        
        cout <<"     SED library contains "<< nsed_ <<" SEDs"<<endl;
        cout <<"     Filter set has "<< nfilt_ <<" filters "<<endl;
        cout << endl;
        
        setLminmax(lmin, lmax, npt);
        
        // initialize color array
        int ndim = 2;
        sa_size_t mydim[ndim];
        mydim[0] = nsed_;
        mydim[1] = nfilt_-1;
        colors_.SetSize(ndim, mydim);
        
        // call some method that calculates the colors
        for (int i=0; i<nsed_; i++) {
            for (int j=0; j<nfilt_-1; j++) {
                colors_(i,j) = calcColors(i,j);
                //cout << colors_(i,j) <<"  ";
                }
            //cout << endl;
            }
        
        };
    
    /** Return color for sed indexed by sedid and filters indexed by filtid and filtid+1. Filter index filtid
        must point to the bluer filter and filtid+1 to the redder filter
        @param sedid   index of SED to calculate color for
        @param filtid  index of redder filter */
    double calcColors(int sedid, int filtid) {
        double z=0.;
        double col = CompColor(z, (*sedarray_[sedid]), (*filtarray_[filtid]), (*filtarray_[filtid+1]));
        return col;
        };
        
        
    /** Return the rest-frame magnitude in each filter instead, useful if filter list contains more than
        one filter set.
        @param magnorm      magnitude in filter ifiltnorm
        @param ifiltnorm    ifiltnorm, magnitude set to be some value in this filter                        */
    TArray<double> getMags(double magnorm = 22., int ifiltnorm = 3) {
    
        // initialize magnitudes array
        TArray<double> mags;
        int ndim = 2;
        sa_size_t mydim[ndim];
        mydim[0] = nsed_;
        mydim[1] = nfilt_;
        mags.SetSize(ndim, mydim);
        
        // initialize filter with normalization magnitude
        for (int ised=0; ised<nsed_; ised++)
            mags(ised, ifiltnorm) = magnorm;
    
        double jmax = nfilt_ - ifiltnorm;
        for (int ised=0; ised<nsed_; ised++) {
        
            // fill from ifiltnorm upwards
            for (int j=1; j<jmax; j++) {
                mags(ised, ifiltnorm + j) = mags(ised, ifiltnorm + j - 1) - colors_(ised, ifiltnorm + j - 1);
                }
 
            // fill from ifiltnorm downwards
            for (int j=1; j<(ifiltnorm+1); j++) {
                mags(ised, ifiltnorm - j) = colors_(ised, ifiltnorm - j) + mags(ised, ifiltnorm - j + 1);
                }
            }   
        return mags;
         
        };
    
    
    /** Return index of SED with closest colors to those input. The SEDs are zero indexed and correspond
        to the index of the SED in sedarray_
        @param colors    vector of colors of a galaxy, must have size() = nfilt_-1 */
    int bestSED(vector<double> colors) {
        
        // first check size of input colors
        if (colors.size() != nfilt_-1)
            throw ParmError("ERROR! number of colors supplied does not match filter set");
        
        // initialize minimum test and current minimum logged
        double min_color_diff = 1e10;
        int closest_sed = -1; // non sensical
        
        // loop over each SED
        for (int i=0; i<nsed_; i++) {
        
            // calculate the sum of differences squared between input colors and SED colors 
            double sumsq_color_diffs=0;
            for (int j=0; j<nfilt_-1; j++)
                sumsq_color_diffs += ((colors[j] - colors_(i,j))*(colors[j] - colors_(i,j)));
                
            // test to see if this is the minimum value seen, if so re log
            if (sumsq_color_diffs < min_color_diff) {
                min_color_diff = sumsq_color_diffs; // set new minimum
                closest_sed = i; // set new closest SED
                }
            }
        
        if (closest_sed<0 || closest_sed>=nsed_)
            throw ParmError("ERROR! failed to find SED with closest colors");
        
        return closest_sed;
            
    };
    
    /** For checking */
    TArray<double> returnColorArray() { return colors_; };
    
        
    /** Output color array (each SED in row, each filter_i - filter_i+1 in column) to a file
        @param fname     name of file to output color array to                                              */
    void writeColorArray(string fname) {
    
        ofstream outp(fname.c_str(), ofstream::out);
        for (int i=0; i<nsed_; i++) {
            for (int j=0; j<nfilt_-1; j++) {
                outp << colors_(i,j) <<"  ";
                }
            outp << endl;
            }
        outp.close();
        
        };
        
    /** Output magnitude array (each SED in row, each filter in column) to a file
        @param fname     name of file to output magnitude array to  
        @param magnorm      magnitude in filter ifiltnorm
        @param ifiltnorm    ifiltnorm, magnitude set to be some value in this filter                        */
    void writeMagsArray(string fname, double magnorm = 22., int ifiltnorm = 3) {
    
        TArray<double> mags = getMags(magnorm, ifiltnorm);
    
        ofstream outp(fname.c_str(), ofstream::out);
        for (int i=0; i<nsed_; i++) {
            for (int j=0; j<nfilt_; j++) {
                outp << mags(i,j) <<"  ";
                }
            outp << endl;
            }
        outp.close();
        
        };


protected:
    vector<SED*> sedarray_;      /**< list of SED templates                                    */
    vector<Filter*> filtarray_;  /**< list of filters                                          */
    int nsed_;                   /**< number of SED templates                                  */
    int nfilt_;                  /**< number of filters                                        */
    TArray<double> colors_;      /**< color of SED for each filter pair, size nsed_*(nfilt_-1) */

};


/** @class
  * ReadKCorrections class
  * 
  * Reads in k-correction tables from files and stores them in interpolation
  * functions
  *
  */
class ReadKCorrections {
public:

    /** Constructor: read k-corrections from files and interpolate
        @param sedLib           name of SED library used to calculate k-corrections
        @param filtSet          name of relevant filter set
        @param restFrameFilt    name of rest-frame filter
	    @param zmin             min redshift k-corrections calculated from
        @param zmax             max redshift k-corrections calculated to
        @param nz               number of redshifts k-corrections calculated at
        @param emax             max extinction k-corrections calculated to
        @param ne               num of extinctions k-corrections calculated at
        @param isMadau          add Madau absorption                          */
	ReadKCorrections(string sedLib, string filtSet, string restFrameFilt, double zmin=0.,            
        double zmax=3., int nz=2000, double emax=0.3, int ne=200, bool isMadau=true)
    : sedLib_(sedLib) , filtSet_(filtSet), restFrameFilt_(restFrameFilt) ,
        zmin_(zmin) , zmax_(zmax) , nz_(nz) , emax_(emax) , ne_(ne) , isMadau_(isMadau){  
        
            double dz = (zmax_-zmin_)/(nz-1), de = emax_/(ne-1);
            for (int iz=0; iz<nz; iz++)
                zvals_.push_back(dz*iz);
            for (int ie=0; ie<ne; ie++)
                evals_.push_back(de*ie);
        };
        
    /** Read k-corrections from @param nSED x @param nFilter files and place 
        each into an array of 1D interpolation function pointers.  The array 
        will be of size: @param nSED x @param nFilter x nExt                  */
    void readInterpZ(int nSED, int nFilter);
    
    /** Read k-corrections from @param nSED x @param nFilter files and place 
        each into an array of 2D interpolation function pointers.  The array 
        will be of size: @param nSED x @param nFilter                         */
    void readInterpZExt(int nSED, int nFilter);
        
    /** Return the array of pointers to the k-correction interpolation tables */    
    vector<SInterp2D*> returnkInterpZExt() { return kInterpZExt_; };
        
    // INTERNAL FUNCS
    
    /** Return filename of k-corrections to read in
        @param iSED     id of SED
        @param iFilter  id of filter                                          */
    string getFileName(int iSED, int iFilter);
    
    
    /** Transpose the 2D array @param tab                                     */
    TArray<double> transposeTable(TArray<double> tab)
        {  
            int nDim = 2;
            sa_size_t mydim[nDim];
            mydim[0] = tab.SizeY();
            mydim[1] = tab.SizeX(); 
            
            TArray<double> tabTransposed;
            tabTransposed.SetSize(nDim,mydim);
            for (int i=0; i<mydim[0]; i++)
                for (int j=0; j<mydim[1]; j++)
                    tabTransposed(i,j) = tab(j,i);
                    
            return tabTransposed;
        };
       
protected:
    string sedLib_;         /**< name of SED library used to calculate k-corrections */
    string filtSet_;        /**< name of relevant filter set                  */
    string restFrameFilt_;  /**< name of rest-frame filter                    */
    double zmin_;           /**< min redshift k-corrections calculated from   */
    double zmax_;           /**< max redshift k-corrections calculated to     */
    int nz_;                /**< number of redshifts k-corrections calculated at*/
    double emax_;           /**< max extinction k-corrections calculated to   */
    int ne_;                /**< num of extinctions k-corrections calculated at */
    bool isMadau_;          /**< add Madau absorption                         */
    vector<double> zvals_;  /**< redshift grid of k-correction tables         */
    vector<double> evals_;  /**< extinction grid of k-correction tables       */
    vector<SInterp1D*> kInterpZ_;   /**< interpolation array (z variable only)     */
    vector<SInterp2D*> kInterpZExt_;/**< interpolation array (z and ext variables) */
};


/** @class
  * TemplateChiSquare class
  * 
  * Class to 
  *
  *
  */
class TemplateChiSquare : public PhotometryCalcs {
public:

/** Calculate chi-square distribution as a function SED type and nuisance
        @param sedArray     array holding pointers to SED objects 
        @param filterArray  array holding pointers to Filter objects
        @param su           object holding cosmological parameters and calculations
        @param lmin         minimum wavelength in meters
        @param lmax         maximum wavelength in meters                      
        @param npt          resolution of SEDs                                                              */
        TemplateChiSquare(vector<SED*> sedArray, vector<Filter*> filterArray, SimpleUniverse& su, 
            double lmin=5e-8, double lmax=2.5e-6, int npt=10000)
            : sedArray_(sedArray) , filterArray_(filterArray) , su_(su) {
        
            nsed_ = sedArray_.size();
            nFilters_ = filterArray_.size();
            
            cout <<"     "<< nFilters_ <<" filters added"<<endl;
            cout <<"     "<< nsed_ <<" templates added"<<endl;
            
            setAGrid(0.1,100.,100);
            setLminmax(lmin,lmax,npt);

            if (nsed_>=1000)
                throw ParmError("ERROR! Too many SEDs");
            };
                            
    /** Returns chi-square distribution and best fit parameters */
    TArray<double> galaxyChiSquared(vector<double> obs, vector<double> errors, 
            double zs, int& sedBestFit, double& normBestFit);
	
	/** Returns the value of the chi-square for the observations, errors and 
	    parameter values given*/
    double meritFunction(vector<double> obs, vector<double> errors, double zs, 
                                                        int iSED, double Anorm);

    /** Set the parameter grid for the normalization parameter \f$A\f$ */
    void setAGrid(double aMin, double aMax, int nA){
        aMin_ = aMin; aMax_=aMax; nA_=nA; 
        dA_ = (aMax_ - aMin_)/(nA_-1); };

    // not sure why these have to be public
	SimpleUniverse& su_;            /**< class that holds the cosmological
	                                      parameters and calculations           */
protected:
    vector<SED*> sedArray_;         /**< holds the SEDs in order of: "elliptical" 
                                    //   then "spiral" then "starburst" types */
    vector<Filter*> filterArray_;   /**< holds the filters                    */
	double lmin_;                   /**< minimum wavelength                   */
	double lmax_;                   /**< maximum wavelength                   */
	int nsed_;                      /**< number of SEDs                       */
	int nFilters_;                  /**< number of filters                    */
	double aMin_;                   /**< minimum value of normalization A     */
	double aMax_;                   /**< maximum value of normalization A     */
	double dA_;                     /**< grid step of normalization A         */
	int nA_;                        /**< number of normalization A in grid    */

};


#endif
