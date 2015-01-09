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
    PhotometryCalcs(double lmin = 5e-8, double lmax = 2.5e-6) { 
        setLminmax(lmin, lmax); };
        
    /** Set the wavelength range to do the integrals over */
    void setLminmax(double lmin, double lmax)
        { lmin_ = lmin; lmax_ = lmax; };
        
        
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
	double Kcorr(double z, ClassFunc1D& sed, Filter& filterX, Filter& restFrameFilter);
	
	/** The k-correction calculation when the rest-frame and observed-frame bandpasses
	    are the same:  \f$ K_{xy} = -2.5\log10\left( \frac{1}{1+z} 
	    \frac{\int f_\lambda(\lambda_o/(1+z))\lambda_oX(\lambda_o)d\lambda_o}
	         {\int f_\lambda(\lambda_e)\lambda_oX(\lambda_o)d\lambda_o} \right) \f$
	    @param z                redshift of object \f$z\f$
        @param sed              rest-frame SED of object \f$f_\lambda(\lambda)\f$
        @param filterX          filter object observed in (same as rest-frame filter) \f$X(\lambda)\f$*/
	double Kcorr1Filter(double z, ClassFunc1D& sed, Filter& filterX);

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
	double CompColor(double z, ClassFunc1D& sed, Filter& filterX, Filter& filterY);
	
	/** Calculate rest-frame flux of object in band \f$X(\lambda\f$ in FREQUENCY units: 
        \f$ F_\nu(\lambda^{eff}_e) = 
            \frac{\int f_\lambda(\lambda_e)X(\lambda_e)\lambda_e d\lambda_e}
                 {\int X(\lambda_e)/\lambda_e d\lambda_e} \f$  
        @param sed      rest-frame SED of object \f$f_\lambda(\lambda)\f$
        @param filter   filter object observed in \f$X(\lambda)\f$
        @param zs       redshift of object                                    */
	double restFrameFlux(ClassFunc1D& sed, Filter& filter, double zs);
	            
	/** Rest-frame flux in WAVELENGTH units: 
	    \f$ F_\lambda(\lambda^{eff}_e) = \frac{F_\nu(\lambda^{eff}_e)}{\lambda_e^2} \f$*/
	double restFrameFluxLambda(ClassFunc1D& sed, Filter& filter, double zs) {
	    double f0nu = restFrameFlux(sed,filter, zs);
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
	double convertABMagToFlux(double mag, double zs, double dL, ClassFunc1D& sed, Filter& filterX) {
        double magPart = pow(10,-0.4*mag);
        SEDzFilterProd sedXlambdaXfilter(sed, filterX, 0);// returns sed*lambda*filter
        FilterIntegrator integrandSED(sedXlambdaXfilter, lmin_, lmax_);
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
	double convertABMagToFluxLambda(double mag, double zs, double dL, ClassFunc1D& sed, Filter& filterX, Filter& filterY) {
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
	    return magAB; };
	    
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
	double effectiveFilterWavelength(ClassFunc1D& filterX);
	
	/** Return the value of the maximum transmission of the filter 
	    @param  filterX     filter transmission function
	    @param  lambdaAtMax wavelength at maximum transmission (returned by method)
	    @param  nStep       number of steps in the search between #lmin_, #lmax_*/
	double findFilterMax(ClassFunc1D& filterX, double& lambdaAtMax, int nStep=1000);
	
	/** Find the wavelength closest to the transmission value #trans in filter
	    #iFilter between lmin and lmax
	    @param  filterX     filter transmission function
	    @param  trans       transmission value to find wavelength of 
	    @param  lmin        search for closest wavelength starting at #lmin
	    @param  lmax        search for closest wavelength ending at #lmax
	    @param  nStep       number of steps in the search between #lmin
	                        and #lmax */
	double findFilterTransValue(ClassFunc1D& filterX, double trans, double lmin,
	                                            double lmax, int nStep=1000);
	
	/** Find wavelength edges of filter iFilter 
	    @param lmin             lower wavelength edge of filter (returned by method)
	    @param lmax             upper wavelength edge of filter (returned by method)
	    @param filterX          filter transmission function
	    @param edgeDefinition   percent of filter maximum defined as the "edge"*/
	void findFilterEdges(double& lmin, double& lmax, ClassFunc1D& filterX, 
	                                            double edgeDefinition=0.05);
	                                            
	/** Return rest frame wavelength: \f$ \lambda_e = \frac{\lambda_o}{1+z} \f$
	    @param  lambdaObs   observed wavelength \f$ \lambda_o \f$
	    @param  zs          redshift \f$ z \f$*/                               
    double returnRestFrameWaveLength(double lambdaObs, double zs) {
            double lambdaRF = lambdaObs/(1 + zs);
            return lambdaRF;
            };
    
protected:
    double lmin_;
    double lmax_;
    //SimpleUniverse su_;
};


/** @class
  * SimData class
  * 
  * Class to simulate galaxy observed magnitudes
  * Can also be used just to calculate colors and k-corrections
  *
  * This class should probably inherit from a base class containing the more
  * generic functions such as convertMagToFlux etc
  *
  */
class SimData : public PhotometryCalcs {
public:
	                                    
	/** Constructor: full k-correction calculation
	    @param sedArray     array holding pointers to SED objects 
	    @param filterArray  array holding pointers to Filter objects
	    @param su           object holding cosmological parameters and calculations
	    @param rg           random number generator object
	    @param nElliptcals  number of "elliptical" type galaxies in sedArray
	    @param nSpirals     number of "spiral" type galaxies in sedArray      */
	SimData(vector<SED*> sedArray, vector<Filter*> filterArray, 
	        SimpleUniverse& su, RandomGeneratorInterface& rg, 
	                int nEllipticals = 1, int nSpirals = 2);
	                        
	/** Constructor: read k-corrections from files and interpolate
	    @param su           object holding cosmological parameters and calculations
	    @param rg           random number generator object
	    @param kInterpZExt  array of pointers to k-corr interps
	    @param nFilters     number of filters
	    @param nElliptcals  number of "elliptical" type galaxies in sedArray
	    @param nSpirals     number of "spiral" type galaxies in sedArray      */
	SimData(SimpleUniverse& su, RandomGeneratorInterface& rg, vector<SInterp2D*> kInterpZExt,
	                int nFilters = 6, int nEllipticals = 1, int nSpirals = 2 );
	                            
	                                    
	/** Destructor */
	virtual ~SimData(void){};
	
	// Put at top: everything that definitely should be accessable from 
	// OUTSIDE the class
	
	// SIMULATION FUNCTIONS //

	/** Function to simulate a true magnitude given: z, sed-type, abs-mag, ext, 
	    rest-frame-filter, observed-filter. Assumes **Madau law IGM** if 
	    @param isAddMadau_ is set to true, assumes no IGM absorption if it's not
	    @param zs               redshift
	    @param sedtype          SED type of galaxy, in form of 1.XXX, 2.XXX, 3.XXX
	    @param amag             absolute magnitude in filter ifRF 
	    @param ext              extinction amount in E(B-V) magnitudes 
	    @param ifO              observation filter id
	    @param RestFrameFilter  rest-frame filter amag is defined in */
	double GetMag(double zs, double sedtype, double amag, double ext, int ifO, 
	                                                    Filter& restFrameFilter);
	                                                    
    /** Function to simulate a true magnitude given: z, abs-mag, SED, 
	    rest-frame-filter, observed-filter. Assumes **Madau law IGM** if 
	    @param isAddMadau_ is set to true, assumes no IGM absorption if it's not
	    **This is different to above because takes SED id instead of a type to convert to an id AND no
	    extinction can be added.**
	    @param zs               redshift
	    @param sedID            SED id
	    @param amag             absolute magnitude in filter RestFrameFilter 
	    @param ifO              observation filter id
	    @param RestFrameFilter  rest-frame filter amag is defined in                  */
	double GetMag(double zs, int sedID, double amag, int ifO, Filter& restFrameFilter);
	                       
	                            
	/** Function to simulate a true magnitude given: z, sed-type, abs-mag, ext, 
	    rest-frame-filter, observed-filter. Assumes Madau law IGM if 
	    @param isAddMadau_ is set to true, assumes no IGM absorption if it's not.
	    **Uses k-correction calculated from files.**
	    @param zs               redshift
	    @param sedtype          SED type of galaxy, in form of 1.XXX, 2.XXX, 3.XXX
	    @param amag             absolute magnitude in filter ifRF 
	    @param ext              extinction amount in E(B-V) magnitudes 
	    @param ifO              observation filter
	    @param RestFrameFilter  rest-frame filter amag is defined in */
	double GetMag(double zs, double sedtype, double amag, double ext, int ifO);
	                            
	                            
	/** Function to simulate a true magnitude given: z, sed-type, abs-mag, ext, 
	    rest-frame-filter, observed-filter. Uses the **line of sight IGM transmission**
	    given by igmTransmission
	    @param zs               redshift
	    @param sedtype          SED type of galaxy, in form of 1.XXX, 2.XXX, 3.XXX
	    @param amag             absolute magnitude in filter ifRF 
	    @param ext              extinction amount in E(B-V) magnitudes 
	    @param ifO              observation filter
	    @param RestFrameFilter  rest-frame filter amag is defined in 
	    @param igmTransmission  line of sight IGM transmission                */
	double GetMag(double zs, double sedtype, double amag, double ext,
	    int ifO, Filter& restFrameFilter, IGMTransmission igmTransmission);
	
	
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
		
	
	/** Add LSST u band error. Returns the observed magnitude and magnitude 
	    error in a vector
	    @param mag      (true) magnitude
	    @param nVisits  number of visits of the telescope in this band        */
	//vector<double> addLSSTuError(double mag, int nVisits);
		
	/** Add LSST g band error */
	//vector<double> addLSSTgError(double mag, int nVisits);
	
	/** Add LSST r band error */
	//vector<double> addLSSTrError(double mag, int nVisits);
	
	/** Add LSST i band error */
	//vector<double> addLSSTiError(double mag, int nVisits);
	
	/** Add LSST z band error */
	//vector<double> addLSSTzError(double mag, int nVisits);
	
	/** Add LSST y band error */
	//vector<double> addLSSTyError(double mag, int nVisits);	
	
	/** Simulate galaxy type, returns a number that corresponds to an SED in
	    #sedArray_
	    @param broad galaxy type, can be equal to 1 (elliptical), 2 (spiral) 
	           or 3 (starburst) */
	double SimSED(int gtype);
	
	/** Simulate reddening amount */
	double SimRed(double type);
	
	
	// SETTINGS FUNCTIONS //
	
	/** Set min and max and wavelength */
	void setXminXmax(double lmin,double lmax)
			{ lmin_=lmin; lmax_=lmax; };
	
	/** Set reddening amount */
	void setRed(double ebvmax,double ebvmaxEl) // if want to change from 0.3,0.1
		{ ebvmax_=ebvmax; ebvmaxEl_=ebvmaxEl; };
		
    /** Set if applying Madau absorption */
    void setMadau(bool isAddMadau, bool isLyC=true)
        { isAddMadau_ = isAddMadau; isLyC_ = isLyC; };
			
	
	// INTERNAL FUNCTIONS: THESE SHOULD BE PROTECTED? //
	
	/** Calculate k correction */
	double calcKcorr(SED& sed, Filter& filter, Filter& restFrameFilter, double zs, double ext, int law);
	
	/** Interpolate k correction  
	    @param sedID        SED id
	    @param iFilterObs   filter id
	    @param zs           redshift
	    @param ext          extinction                                        */
	double interpKcorr(int sedID, int iFilterObs, double zs, double ext);
	/** Interpolate k correction  
	    @param linearIndex  SED,filter combination id
	    @param zs           redshift
	    @param ext          extinction                                        */
	double interpKcorr(int linearIndex, double zs, double ext);
	
	/** Given a row index @param i and a column index @param j, return the single array 
	    element number if there are @param nj columns                         */
	int returnLinearIndex(int i, int j, int nj) { return i*nj + j; };
	
	/** Given an sedtype number generated by #SimSED returns the actual index 
	    of the SED in #sedArray_ 
	    @param sedtype  number generated by #SimSED */
	int returnSedId(double sedtype);
     
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
    
    /** Return fluxes of each SED in each of the rest-frame filters
        @param zs               redshift of the SED
        @param restFrameFilter  filter absolute magnitude is defined in */
    TArray<double> returnSEDFluxesInRestFrame(double zs);//, Filter& restFrameFilter);
    
    vector<double> returnm5single() { return m5single_; };

    // not sure why these have to be public
	SimpleUniverse& su_;            /**< class that holds the cosmological
	                                      parameters and calculations           */
	RandomGeneratorInterface& rg_;  /**< class that generates the random numbers*/

protected:
    vector<SED*> sedArray_;         /**< holds the SEDs in order of: "elliptical" 
                                    //   then "spiral" then "starburst" types */
    vector<Filter*> filterArray_;   /**< holds the filters                    */
    int nEllipticals_;              /**< number of ellipticals in sedArray    */
    int nSpirals_;                  /**< number of spirals in sedArray        */
    int nStarbursts_;               /**< number of starbursts in sedArray =   //
    nStarbursts_ = nsed_ - nEllipticals_ - nSpirals_                          */
    bool isAddMadau_;               /**< add Madau absorption                 */
    bool isLyC_;                    /**< include Lyman continuum in Madau absorption */
    bool isReadKcorr_;              /**< read k corrections from file         */
	int nsed_;                      /**< number of SEDs                       */
	int nFilters_;                  /**< number of filters                    */
	double ebvmax_;                 /**< max extinction to apply to galaxies  */
	double ebvmaxEl_;               /**< max extinction to apply to El galaxy */
	vector<SInterp2D*> kInterpZExt_;  /**< array of pointers to k-corr interpolation */
	// To initialize when the references are not used
	SimpleUniverse su_default_;     /**< to initialize #su_ when it's not used*/
	DR48RandGen rg_default_;        /**< to initialize #rg_ when it's not used*/
	//double uMsky_,gMsky_,rMsky_,iMsky_,zMsky_,yMsky_;
	//double uTheta_,gTheta_,rTheta_,iTheta_,zTheta_,yTheta_;
	//double uGamma_,gGamma_,rGamma_,iGamma_,zGamma_,yGamma_;
	//double uCm_,gCm_,rCm_,iCm_,zCm_,yCm_;
	//double ukm_,gkm_,rkm_,ikm_,zkm_,ykm_;
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


/** @class SEDLibColors
  *
  * Calculate colors of all SEDs in a library for a given filter set
  *
  */
class SEDLibColors : public PhotometryCalcs {
public:
    /** Constructor. The filter set must be ordered by monotonic increase in filter wavelength from blue to red
        @param sedarray   array of pointers to each SED in library
        @param filtarray  array of pointers to each filter in filter set*/
    SEDLibColors(vector<SED*> sedarray,  vector<Filter*> filtarray):
        sedarray_(sedarray), filtarray_(filtarray) , nsed_(sedarray.size()) , nfilt_(filtarray.size()) {
        
        cout <<"     SED library contains "<< nsed_ <<" SEDs"<<endl;
        cout <<"     Filter set has "<< nfilt_ <<" filters "<<endl;
        cout << endl;
        
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
    
    
    /** Return index of SED with closest colors to those input */
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
        
        if (closest_sed<0 || closest_sed>nsed_+1)
            throw ParmError("ERROR! failed to find SED with closest colors");
        
        return closest_sed;
            
    };
    
    /** For checking */
    TArray<double> returnColorArray() { return colors_; };
    
        
    /** For checking */
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
	    @param lmax         maximum wavelength in meters                      */
        TemplateChiSquare(vector<SED*> sedArray, vector<Filter*> filterArray, 
            SimpleUniverse& su, double lmin=5e-8, double lmax=2.5e-6)
            : sedArray_(sedArray) , filterArray_(filterArray) , su_(su) {
        
            nsed_ = sedArray_.size();
            nFilters_ = filterArray_.size();
            
            cout <<"     "<< nFilters_ <<" filters added"<<endl;
            cout <<"     "<< nsed_ <<" templates added"<<endl;
            
            setAGrid(0.1,100.,100);
            setLminmax(lmin,lmax);

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
