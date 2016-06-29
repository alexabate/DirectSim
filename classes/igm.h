/**
 * @file  igm.h
 * @brief Contains classes that provide functions to Monte Carlo simulate the IGM
 * 
 * LC: Contains Inoue and Iwata 2014 method, Madau method, and Meiksin method 
 *
 * Could add more information here I think
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: 20 Aug 2012
 * @date 20 Aug 2012
 *
 */

#ifndef IGM_H_SEEN
#define IGM_H_SEEN

//#include "machdefs.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm> 

// sophya
// Sophya update v2.3 June 2013 replaces genericfunc with classfunc
//#include "genericfunc.h"
#include "classfunc.h"
#include "sopnamsp.h"
#include "mydefrg.h"
#include "pexceptions.h"

// CatSim
#include "constcosmo.h"
#include "geneutils.h"
#include "sinterp.h"


/** @class AtomicCalcs
  *
  * Class that holds basic atomic calculations and constants mostly to do with
  * the Lyman series
  * 
  *
  */
class AtomicCalcs
{
public:
    /** Constructor (sets constants) */
    //AtomicCalcs() { setGammas(); setOscillatorStrength(); setConstants(); };
    AtomicCalcs();


    /** Return the Lyman Series wavelength in meters
        @param    n is the starting energy level for a Lyman transition (down to n=1) (so minimum n is 2)   */
    double returnWavelengthLymanSeries(int n) ;
        
    
    /** Return the Lyman Series wavelength in Angstroms
        @param    n is the starting energy level for a Lyman transition (down to n=1) (so minimum n is 2)   */
    double returnWavelengthAngstrLymanSeries(int n)  {
        if (n<2) throw ParmError("ERROR! Lyman transition cannot start at n<2");
        return WAVE_LYMANLIM_ANGSTR/(1.-1./(n*n));
        };
        
    
    /** Return the Lyman Series wavelength in frequency (in s^-1)
        @param  n is the starting energy level for a Lyman transition (down to n=1) (so minimum n is 2)     */
    double returnFrequencyLymanSeries(int n) {
        if (n<2) throw ParmError("ERROR! Lyman transition cannot start at n<2");
        double wl = returnWavelengthLymanSeries(n);
        return SPEED_OF_LIGHT_MS/wl;
        };


    /** Return Doppler width in meters with line center wavelength corresponding
        to Lyman series n \f$ \Delta\lambda=\lambda_i\frac{b}{c} \f$
        @param n             starting energy level of the Lyman line transition
        @param dopplerPar    doppler parameter in km/s                                                      */
    double returnDopplerWidthWL(int n, double dopplerPar) {
        double ratio = dopplerPar/SPEED_OF_LIGHT_KMS; // b/c
        return returnWavelengthLymanSeries(n)*ratio;  // lambda_i*(b/c)
        };


    /** Return Doppler width in s^-1 with line center frequency corresponding to 
        Lyman series n \f$ \Delta\nu=\nu_i\frac{b}{c} \f$
        @param n              starting energy level of the Lyman line transition
        @param dopplerPar     doppler parameter in km/s                                                     */
    double returnDopplerWidthFreq(int n, double dopplerPar) {
        double ratio = dopplerPar/SPEED_OF_LIGHT_KMS; // b/c
        return returnFrequencyLymanSeries(n)*ratio;   // freq_i*(b/c)
        };


    /** Return the (unitless) wavelength difference relative to resonant wavelength
        in Doppler units
        @param lambda        wavelength in meters                                 
        @param dopplerPar    doppler parameter in km/s                                                      */
    double returnX(double lambda, int n, double dopplerPar) {
        double dl = (lambda - returnWavelengthLymanSeries(n));
        double dw = returnDopplerWidthWL(n, dopplerPar); 
        return dl/dw; // wavelength difference relative to resonant wavelength
        };


    /** Return the damping constant for the Lyman line beginning at n     
        @param n             starting energy level of the Lyman line transition                             */
    double returnGamma(int nLine);


    /** Return the Oscillator Strength for the Lyman line beginning at n  
        @param n             starting energy level of the Lyman line transition                             */
    double returnOscillatorStrength(int nLine);


    /** Return (unitless) damping parameter a 
        @param n             starting energy level of the Lyman line transition 
        @param dopplerPar    doppler parameter in km/s                                                      */          
    double returnDampingParameter(int n, double dopplerPar) {
        double li = returnWavelengthLymanSeries(n);
        double gammai = returnGamma(n);
        double ld = returnDopplerWidthWL(n, dopplerPar);

        double numer = li*li*gammai;
        double denom = 4*PI*SPEED_OF_LIGHT_MS*ld;

        return numer/denom;
        };


    /** Print all calculations for given n and dopplerPar
        @param n             starting energy level of the Lyman line transition 
        @param dopplerPar    doppler parameter in km/s                                                      */
    void printEverything(int n, double dopplerPar);
    
    
    /** Return the maximum energy level allowed for a Lyman line transition calculation                     */
    int returnnLineMax() { return nLineMaxMax_; };

    //    void setGammas();
    //    void setOscillatorStrength();
    void setConstants();

protected:
    std::vector<double> gammaSeries_;   /**< damping constant of the Lyman series in s^-1 */
    //    std::vector<double> fSeries_;       /**< oscillator strength of the Lyman series      */
    double sigmaLymanLimitCM2_;         /**< cross-section at the lyman limit in cm^2     */
    double freqLymanLimitInvSec_;       /**< frequency of the lyman limit in s^-1         */
    int nLymanAlpha_;                   /**< starting level of a lyman alpha transition   */
    int nLineMaxMax_;                   /**< maximum Lyman series line possible to use    */
    int nGammaMax_;
	double R, s; 
};




/** HIColumnDensityLAF 
 * Inoue 2014 
 * class currently just holds values of constants
 * no normalization because not normalizing like paper for now 
 *
 */
class HIColumnDensityLAF
{
public:
    HIColumnDensityLAF();

//    void normalizeDist();

    /** Functions that pass the various members variables by reference  */
    void returnPowerLawIndex(double &betaLAF);
    void returnColDensityLimits(double &Nl, double &Nu);
    void returnColDensityBreak(double &Nc);

    void testClass();

protected:
    double betaLAF_;      /**< Power for the LAF power law          */
    double Nl_;         /**< Low bound for the column densities     */
    double Nc_;         /**< Upper bound for the column densities   */
    double Nu_;         /**< The central distrbution break density  */
//    double normB_;      /**< The normalization constant             */
};





/** HIColumnDensityDLA 
 * Inoue 2014 
 * class currently just holds values of constants
 * no normalization because not normalizing like paper for now 
 *
 */
class HIColumnDensityDLA
{
public:
    HIColumnDensityDLA();

//    void normalizeDist();

    /** Functions that pass the various members variables by reference  */
    void returnPowerLawIndex(double &betaDLA);
    void returnColDensityLimits(double &Nl, double &Nu);
    void returnColDensityBreak(double &Nc);

    void testClass();

protected:
    double betaDLA_;      /**< Power for the DLA power law          */
    double Nl_;         /**< Low bound for the column densities     */
    double Nc_;         /**< Upper bound for the column densities   */
    double Nu_;         /**< The central distrbution break density  */
//    double normB_;      /**< The normalization constant             */
};






/** AbsorberRedshiftDistributionLAF
 *
 * Class to calculate the absorber redshift distribution function for LAF
 * Based on Inoue 2014 eqn 6 f(z)
 *
 */
class AbsorberRedshiftDistributionLAF 
{
public:
    /** Constructor */
    AbsorberRedshiftDistributionLAF();

    /** Return redshift distribution at value given */
    double returnRedshiftDist(double z);
    
    /** Return power law at z value given */
    double returnFirstPowerLaw(double z);

    /** Return power law at z value given */
    double returnSecondPowerLaw(double z);

    /** Return power law at z value given */
    double returnThirdPowerLaw(double z);

    /** Return the power law indices */
    void returnPowerLawIndex(double &g1LAF, double &g2LAF, double &g3LAF);

    /** Return the redshift at the breaks */
    void returnRedshiftBreaks(double &z1LAF, double &z2LAF);

    /** Return the normalization constant */
    double returnNormalization();

    /** Return redshift distribution at z value given */
    virtual double operator()(double z)
        { return returnRedshiftDist(z); }

    void testClass();

protected:
    double ALAF_;          /**< (Normalization) Total number of absorber z=z1 with a column density 1e12<=NHI<=1e23 cm^-2 */
    double z1LAF_;         /**< Redshift value at first break */
    double z2LAF_;         /**< Redshift value at second break */
    double gamma1LAF_;     /**< Power law index */
    double gamma2LAF_;     /**< Power law index */
    double gamma3LAF_;     /**< Power law index */
};





/** AbsorberRedshiftDistributionDLA
 *
 * Class to calculate the absorber redshift distribution function for DLA
 * Based on Inoue 2014 eqn 7 f(z)
 *
 */
class AbsorberRedshiftDistributionDLA 
{
public:
    /** Constructor */
    AbsorberRedshiftDistributionDLA();

    /** Return redshift distribution at value given */
    double returnRedshiftDist(double z);
    
    /** Return power law at z value given */
    double returnFirstPowerLaw(double z);

    /** Return power law at z value given */
    double returnSecondPowerLaw(double z);

    /** Return the power law indices */
    void returnPowerLawIndex(double &g1DLA, double &g2DLA);

    /** Return the redshift at the breaks */
    void returnRedshiftBreaks(double &zDLA);

    /** Return the normalization constant */
    double returnNormalization();

    /** Return redshift distribution at z value given */
    virtual double operator()(double z)
        { return returnRedshiftDist(z); }

    void testClass();

protected:
    double ADLA_;          /**< (Normalization) Total number of absorber z=z1 with a column density 1e12<=NHI<=1e23 cm^-2 */
    double zDLA_;         /**< Redshift value at the only break */
    double gamma1DLA_;     /**< Power law index */
    double gamma2DLA_;     /**< Power law index */
};




/** DopperParDistribution
 *
 * Class to calculate the doppler parameter b distribution
 * Based on Inoue & Iwata 2008 and Hui & Rutledge 1999
 * We're assuming it's the same for Inoue 2014 
 */
class DopplerParDistribution
{
public:
    DopplerParDistribution();
    double returnDopplerDist(double b);

    virtual double operator()(double b)
        { return returnDopplerDist(b); }

    void testClass();

protected:
    double bsigma_;
};




/** ProbabilityDistAbsorbers
 *
 * Simulates a line of sight
 * Simulates some number of absorbers along this line of sight 
 * For each absorber, draws column density, redshift, and doppler parameter
 * Draws deltaZ of the next absorber 
 * Repeats until last absorber is at the set source redshift 
 * 
 */
class ProbabilityDistAbsorbers
{
public:
    ProbabilityDistAbsorbers(RandomGeneratorInterface& rg, double zS, 
 			     AbsorberRedshiftDistributionLAF& absorberZDistLAF,
                             AbsorberRedshiftDistributionDLA& absorberZDistDLA,
                             HIColumnDensityLAF& hiColumnDensityLAF,
                             HIColumnDensityDLA& hiColumnDensityDLA,
                             DopplerParDistribution& dopplerParDist); 



    // Define the min and max values for the column density and
    // the density distribution
//    void setNHiDistribution(int nStep);  //doesn't seem to do anything in old class 


    // Define the min and max values for the doppler parameter 
    // and the doppler param distribution
    void setDopplerDistribution(int nStep);

    /** Using the inverse transform method and equation 7 in 
        Inoue & Iwata 2008      */
//>>>>>>> b174eb658571fc61c81d79c7fa2a0db330092148
    double drawDeltaZLAF(double zLast);

    double drawDeltaZDLA(double zLast);


    /** Draw a column density from the distribution given in 
       Inoue & Iwata 2014 equation 5; done by numerical integration       */
    double drawHIColumnDensityLAF();
    double drawHIColumnDensityDLA();


    /** Draw a doppler parameter from the distribution given in
       Inoue & Iwata 2008 equation 6        */
    double drawDoppler();

    /** Simulate a single absorber
        @param zCurrent     the redshift of the previous absorber
        @param zNext        the redshift of the absorber being simulated
        @param bdopp        the doppler param of the absorber being simulated
        @param NHI          the HI column density of the absorber being simulated
    */
    void simulateAbsorberLAF(double zCurrent, double& zNext, double& bdopp, double& NHI);
    void simulateAbsorberDLA(double& bdopp, double& NHI);

    /** Simulate a line of sight 
        @param zStart           the redshift to begin at
        @param zMax             the redshift to end at
        @param redshifts        the vector to hold the redshift of each absorber
        @param dopplerPars      the vector to hold the doppler parameter of each absorber
        @param columnDensities  the vector to hold the column density of each absorber
    */

   void simulateLineOfSightLAF(double zStart, double zMax,
                    vector<double>& redshiftsLAF, vector<double>& dopplerParsLAF,
                    vector<double>& columnDensitiesLAF, string outfile); 

  void simulateLineOfSightDLA(double zStart, double zMax,
                    vector<double>& redshiftsDLA, vector<double>& dopplerParsDLA,
                    vector<double>& columnDensitiesDLA, string outfile); 

protected:
    RandomGeneratorInterface&       rg_;
	double zS; 
    AbsorberRedshiftDistributionLAF&   absorberZDistLAF_;
    AbsorberRedshiftDistributionDLA&   absorberZDistDLA_;
    HIColumnDensityLAF&                hiColumnDensityLAF_;
    HIColumnDensityDLA&                hiColumnDensityDLA_;
    DopplerParDistribution&         dopplerParDist_;

//    double log10Nl_;
//    double log10Nu_;
//    vector<double> log10NHIvals_;
//    vector<double> log10gvals_;
//    double log10gmin_;
//    double log10gmax_;

//    SInterp1D colDensityFunc_;

    vector<double> bvals_;
    vector<double> hvals_;
    double bmin_;
    double bmax_;
    double hmin_;
    double hmax_;
    vector<double> randomVectorLAF, columndensityVectorLAF;
    vector<double> randomVectorDLA, columndensityVectorDLA;
	vector<double> zVector, yLAFvector, yDLAvector; 
	vector<double> cppLAFvector, cppDLAvector; 

};
                                                  



class VoigtProfile:
    public AtomicCalcs
{
public:
    VoigtProfile(double dopplerPar, int nLine);

    double kFunction(double x);

    double returnHax(double lambda);

    virtual double operator()(double lambda) { 
        return returnHax(lambda); 
    };
    
protected:
    double dopplerPar_;         /*< The doppler parameter of the absorber         */
    int nLine_;                 /*< The line number to being the Lyman transition */

};




/*Class is kept exactly the same from Inoue 2008*/
class OpticalDepth:
    public AtomicCalcs
{
public:
    OpticalDepth();

    void setLymanAll();
    void setLymanContinuumOnly();
    void setLymanSeriesOnly();
    void setContribution(int option);

    /** Return the transmission in the observer's frame from an absorber with redshift zAbsorber,
        column density nhiAbsorber, and doppler parameter bAbsorber
        @param lambda observed wavelength in m
        @param zAbsorber redshift of absorber
        @param nhiAbsorber column density of absorber
        @param bAbsorber doppler parameter of absorber          */
    double returnObserverFrameTransmission(double lambda, double zAbsorber, double nhiAbsorber, double bAbsorber);

    /** Return the optical depth in the observer's frame due to an absorber with redshift zAbsorber,
        column density nhi, and doppler parameter bAbsorber
        @param lambda observed wavelength in m
        @param zAbsorber redshift of absorber
        @param nhi column density of absorber
        @param bAbsorber doppler parameter of absorber          */
    double returnObserverFrameOpticalDepth(double lambda, double zAbsorber, double nhi, double bAbsorber);

    /** Return the optical depth in the rest frame of the absorber
        @param freq rest frame frequency in 1/s
        @param bAbsorber doppler parameter of the absorber in km/s
        @param nhi the column density of the absorber in cm^-2  */
    double returnRestFrameOpticalDepth(double freq, double bAbsorber, double nhi);

    double returnLymanContinuumCrossSection(double freq);

    double returnLymanSeriesCrossSection(double freq, double bAbsorber);

    double returnLymanLineCrossSection(int n, double freq, double bAbsorber);

    double returnLineProfile(int n, double freq, double bAbsorber);

    void setMaxLine(int nLine);

protected:
    bool isLymanC_;
    bool isLymanS_;
    int nLineMax_;

};






class LineOfSightTrans
: public OpticalDepth
{
public:
    LineOfSightTrans(vector<double>& redshifts, vector<double>& dopplerPars,
                     vector<double>& columnDensities);

    /** Return the total optical depth at observed wavelength lambda
        @param lambda   observed wavelength in m
        @param zSource  redshift of source      */
    double returnOpticalDepth(double lambda, double zSource);

    double returnTransmission(double lambda, double zSource);

    int returnNumberOfAbsorbers(double zSource);

    void setReturnType(bool isOpticalDepth);

    virtual double operator()(double lambda, double zSource) {
        if(isOpticalDepth_)
            return returnOpticalDepth(lambda, zSource);
        else
            return returnTransmission(lambda, zSource);
    }; 

protected:
    vector<double> redshifts_;
    vector<double> dopplerPars_;
    vector<double> columnDensities_;
    bool isOpticalDepth_;
};








/** Madau class
  *
  */
class Madau:
    public AtomicCalcs
{
public:

    //Madau() {}; 


    /** Constructor
        @param nLineMax     Maximum Lyman-series to include (eg if nLineMax=2 will only calculate Ly-alpha) 
        @param isLyC        Include Lyman continuum absorption (not just Lyman series)                      */
    Madau(int nLineMax, bool isLyC=true)
    : nLineMax_(nLineMax) , isLyC_(isLyC) { 
        setAbsorptionStrengths(); 
        if (nLineMax_>nLineMaxMadau_) {
            cout <<"     Cannot calculate Madau absorption for n > "<< nLineMaxMadau_ <<", setting n = ";
            cout << nLineMaxMadau_ << endl;
            }
        };
    
    /** Return the transmission in the observer's frame along a line of sight to a source at some redshift z
        @param lambda       observed wavelength in meters   
        @param zSource      redshift of source                  */
    double returnObserverFrameTransmission(double lambda, double zSource)
        { double tau = returnObserverFrameOpticalDepth(lambda, zSource);
          return exp(-tau); };
          
          
    /** Return the transmission in the rest-frame along a line of sight to a source at some redshift z
        @param lambda       rest-frame wavelength in meters   
        @param zSource      redshift of source                  */
    double returnRestFrameTransmission(double lambda, double zSource) { 
            double lambdaObs = lambda*(1. + zSource);
            double tau = returnObserverFrameOpticalDepth(lambdaObs, zSource);
            return exp(-tau); };
    
    
    /** Return the optical depth in the observer's frame along a line of sight to a source at redshift z
        @param lambda       observed wavelength in meters                          
        @param zSource      redshift of source                                */
    double returnObserverFrameOpticalDepth(double lambda, double zSource);
          
    /** Return the optical depth in observer's frame along a line of sight due to Lyman continuum absorption
        only                                                                                                */
    double returnLymanContinuumOpticalDepth(double zAbsorber, double zSource); 


protected:    
    /** Set constants of Lyman series absorption strengths                    */
    void setAbsorptionStrengths();


protected:
    int nLineMax_;          /**< Maximum Lyman series to include                                            */
    int nLineMaxMadau_;     /**< Maximum Lyman series can include (due to model limits)                     */
    vector<double> Avals_;  /**< Absorption strength of Lyman-alpha to delta                                */
    bool isLyC_;            /**< Include Lyman continuum absorption                                         */


};

// Below here are classes to calculate a Meiksin 2006-style Monte Carlo simulation

/** SimulateLLS class
  *
  */
class SimulateLLS {
public:

    SimulateLLS(RandomGeneratorInterface& rg, double N0=0.25, double gamma=1.5, double beta=1.5)
    : rg_(rg) , N0_(N0) , gamma_(gamma) , beta_(beta) 
    {   zMax_ = 8.;
        opDMax_ = 1.7E+308;
        nStep_ =1000; };
    
    
    long simulateLineOfSight(double zSource, vector<double>& zAbsorbers, vector<double>& opticalDepth)
        {   long nLine = nLineofSight(zSource);
            
            if (zSource>zMax_)
                throw ParmError("ERROR! Cannot simulate out to this high a redshift");
            
            setRedshiftDistribution(nStep_);
            setOpDepthDistribution(nStep_);
            cout <<"     Set distributions "<<endl;
            cout <<"     Simulating "<< nLine <<" absorbers "<<endl;
            for (int i=0; i<nLine; i++) {
                double z = drawRedshift(zSource);
                double t = drawOpticalDepth();
                
                zAbsorbers.push_back(z);
                opticalDepth.push_back(t);
                }
            return nLine;
        };
    
    /** returns number of LLS's along lines of sight by Poisson fluctuating the
        mean number over line of sights                                       */
    long nLineofSight(double z) {   
            long nMean = calculateMeanNo(z);
            long nLine = rg_.Poisson(nMean); 
            return nLine; 
            };
    
    long calculateMeanNo(double z) { 
            // integral of N0*(1+z)^gamma  
            long num = round( (N0_*pow(1.+z, gamma_+1.))/(gamma_+1.) );
            return num;   
            };
            
    double drawRedshift(double zSource) {
        // inverted dN/dZ = (1+z)^gamma
        double random = rg_.Flat01();
		double z = pow(random * (pow((1 + zSource), gamma_ + 1) - 1) + 1,
				pow((gamma_ + 1), -1)) - 1;
		return z;
    
        /*double u1,u2;
        while (true)  {
            u1 = zMax_*rg_.Flat01();
            u2 = nzmin_ + (nzmax_ - nzmin_)*rg_.Flat01();
        
            int iElement=findClosestElement(zvals_,u1);
       
            if ( u2 <= nzvals_[iElement])
                break;
            }
        return u1;*/
        
        };
        
    double drawOpticalDepth()
        {
        // inverted dN/dtau = tau^-beta
        double random =  rg_.Flat01();
		double tau = pow(1 + random * (pow(opDMax_, 1 - beta_) - 1),
				pow((1 - beta_), -1));
		return tau;
        
         /*   double u1,u2;
            while (true)  {
                u1 = opDMax_*rg_.Flat01();
                u2 = ntmin_ + (ntmax_ - ntmin_)*rg_.Flat01();
            
                int iElement=findClosestElement(tvals_,u1);
           
                if ( u2 <= ntvals_[iElement])
                    break;
                }
            return u1;*/
            
        };
          
    void setRedshiftDistribution(long nStep)
        {  
            zvals_.clear();
            nzvals_.clear();
    
            double dz = zMax_/(nStep-1);
    
            for (int i=0; i<nStep; i++) {
                double zv = i*dz;

                zvals_.push_back(zv);
                double nz = N0_*pow(1.+zv,gamma_);
                nzvals_.push_back(nz);
                }

            nzmin_=findMinimum(zvals_);
            nzmax_=findMaximum(nzvals_);
        };
    
    void setOpDepthDistribution(long nStep)
        {  
            tvals_.clear();
            ntvals_.clear();
    
            double dt = opDMax_/(nStep-1);
    
            for (int i=0; i<nStep; i++) {
                double tv = i*dt;

                tvals_.push_back(tv);
                double nt = pow(tv,-beta_);
                ntvals_.push_back(nt);
                }

            ntmin_=findMinimum(tvals_);
            ntmax_=findMaximum(ntvals_);
        };
        
    void setZmax(double zMax){ zMax_ = zMax; };
    void setNstep(long nStep){ nStep_ = nStep; };
    
protected:
    RandomGeneratorInterface& rg_; /**< random generator                      */
    double N0_;     /**< parameter of LLS redshift distribution               */
    double gamma_;  /**< parameter of LLS redshift distribution               */
    double beta_;   /**< parameter of LLS optical depth distribution          */
    double zMax_;   /**< maximum possible redshift of an absorber             */
    long nStep_;    /**< number of steps in z, tau distributions              */
    vector<double> zvals_;
    vector<double> nzvals_;
    double nzmin_;
    double nzmax_;
    double opDMax_;
    vector<double> tvals_;
    vector<double> ntvals_;
    double ntmin_;
    double ntmax_;
};

/** LAFMeiksin class
  *
  * Calculates optical depth due to the lyman alpha forest according to Meiksin 
  * 2006. Returns a value given the observed wavelength and the redshift of the 
  * source.
  *
  */
class LAFMeiksin :
    public AtomicCalcs 
{
public:

    /** Constructor */
    LAFMeiksin(int prt=0){ prt_=prt; setCoeffs(); };
    
    /** returns transmission at observed wavelength in meters for light emitted
        from a source at zSource
        @param lambdaObs    observed wavelength in meters
        @param zSource      redshift of source                                */
    double returnTransLymanSeries(double lambdaObs, double zSource) {
        double tau = returnTauLymanSeries(lambdaObs, zSource);
        return exp(-tau);
        };
    
    /** returns optical depth at observed wavelength in meters for light emitted
        from a source at zSource
        @param lambdaObs    observed wavelength in meters
        @param zSource      redshift of source                                */
    double returnTauLymanSeries(double lambdaObs, double zSource);
    
    /** returns optical depth due to Lyman-alpha 
        @param zLyn     redshift of absorber
        @param zSource  redshift of source                                    */
    double returnTauLyA(double zLyn, double zSource) { 
        double tauLyA=0;
        if (zLyn < zSource) { // means Lyman-alpha absorption *IS OBSERVABLE*
            if (zLyn < 4.) {
                if (zLyn < 1.2)
                    tauLyA = 0.0164*pow((1. + zLyn), 1.1);  // eq ?: 0<z<1.2
                else
                    tauLyA = 0.00211*pow((1. + zLyn), 3.7); // eq 2: 1.2<z<4
				} 
            else
                tauLyA = 0.00058*pow((1. + zLyn), 4.5);     // eq 3: z>4
			}
        return tauLyA;
        };
    
    /** return nmax observable 
        @param lambdaObs    observed wavelength in meters
        @param zSource      redshift of source                                */
    int findNmax(double lambdaObs, double zSource) {
        int nmax;
        double ratio = lambdaObs / WAVE_LYMANLIM_METERS;
	    if (ratio < (1. + zSource))
		    nmax = 31; // use full series up to nMax of 31 
	    else
		    nmax = (int)pow((1. - (1. + zSource)/ratio), -0.5);
	    if (nmax > 32)
		    nmax = 31;
		return nmax;
        }; 
       
    /** returns a useful quantity, see Meiksin 2006 Table 1*/
    double returnZvalue(double z, int n){
            double zVal;
            if (z >= 3.&& n<6 )
                zVal = pow(0.25*(1. + z), (1./6.));
			else
                zVal = pow(0.25*(1. + z), (1./3.));
            return zVal;
            };
            
     void setCoeffs();
     
     double taulya(double lambda, double z);

protected:
    int prt_;                   /**< print level */
    vector<double> tauCoeff_;   /**< see Meiksin 2006 Table 1 */

};

/** PhotoElectricMeiksin class
  *
  * Calculates optical depth due to the diffuse, optically thin IGM and due to 
  * an LLS absorber
  *
  */
class PhotoElectricMeiksin  :
    public AtomicCalcs 
{
public:
    /** Constructor */
    PhotoElectricMeiksin(double A=0.07553){ A_=A; };


    /** return optical depth due to diffuse IGM (the contribution from all 
        optically thin systems)
        @param lambdaObs    observed wavelength in meters
        @param zSource      redshift of source                                */
    double returnTauDiffuse(double lambdaObs, double zSource) {

		double zLL = (lambdaObs/WAVE_LYMANLIM_METERS - 1);
		
		double tau = 0;
		if (zSource > zLL) {
		    // equation 5 in Meiksin (modified)
		    tau = A_*pow((1+zLL), 4.4)*((1/pow(1+zLL, 1.5))-(1/pow(1+zSource, 1.5)));	
		    }
		return tau;
	    };
	    
	/** return optical depth due to an LLS absorber
        @param lambdaObs    observed wavelength in meters
        @param tauVal       optical depth-like value drawn from distribution
        @param zAbsorber    redshift of LSS absorber                          */
    double returnTauLLS(double lambdaObs, double tauVal, double zAbsorber) {

        double tau = 0;
        if (lambdaObs < WAVE_LYMANLIM_METERS*(1+zAbsorber) ) {
            double k = lambdaObs/(WAVE_LYMANLIM_METERS*(1+zAbsorber));
		    tau = tauVal*k*k*k;
		    }
        return tau;
	    };
	    
	void setA(double A){ A_=A; };
	void returnA(double& A){ A=A_; };
	
protected:
    double A_; /**< diffuse IGM normalization parameter */
	    
};

/** MonteCarloMeiksin class
  *
  *
  */
class MonteCarloMeiksin :
    public LAFMeiksin, public PhotoElectricMeiksin
{
public:
    MonteCarloMeiksin(vector<double> zAbsorbers, vector<double> opticalDepths)
    : zAbsorbers_(zAbsorbers) , opticalDepths_(opticalDepths){};
            
    /** Return optical depth at observed wavelength lambda in meters of light
        from a source located at zSource, returns \f$\tau\f$
        @param lambdaObs    observed wavelength
        @param zSource      redshift of background source                     */
    double returnOpticaldepth(double lambdaObs, double zSource);
    
    /** Return transmission at observed wavelength lambda in meters of light
        from a source located at zSource, returns \f$e^{-\tau}\f$
        @param lambdaObs    observed wavelength
        @param zSource      redshift of background source                          */
    double returnTransmission(double lambdaObs, double zSource) {
        double tau = returnOpticaldepth(lambdaObs, zSource);
        return exp(-tau);
        };
        
    double returnTauLymanLimitSystems(double lambdaObs, double zSource);
        
    double returnLySeriesOnly(double lambdaObs, double zSource) {
        double tauLAF = returnTauLymanSeries(lambdaObs,zSource);// Lyman-alpha forest
        return tauLAF;
        };
        
    double returnDiffuseOnly(double lambdaObs, double zSource) {
        double taudIGM = returnTauDiffuse(lambdaObs,zSource);   // diffuse IGM
        return taudIGM;
        };
    
    double returnLLSOnly(double lambdaObs, double zSource) {
        double tauLLS = returnTauLymanLimitSystems(lambdaObs,zSource);   // Lyman-limit systems
        return tauLLS;
        };
        
protected:
    vector<double> zAbsorbers_;
    vector<double> opticalDepths_;
    double lmin_;
    double lmax_;
    int nl_;
};


/** IGMTransmission class
  *
  * Holds IGM transmission for a galaxy at fixed source redshift. Returns the transmission of a galaxy at
  * some redshift but in the REST-FRAME of the galaxy (not the observer)
  *
  * Either:
  * 
  * - Reads in the transmission values from a file. These values could be the mean transmission across all 
  *   lines of sight or for a particular line of sight
  *   @warning the IGM transmission only applies for a certain source galaxy redshift
  *
  *   File contains: observed wavelength (m)  transmission (exp[-tau])
  *   @note REALLY MUST CHANGE FILE DATA TO BE REST-FRAME WAVELENGTH
  *
  * - Calculates Madau transmission for a given source redshift
  *
  * - Returns transmission consistant with no IGM absorption, i.e. exp[-tau]=1
  *
  */
class IGMTransmission //: public SInterp1D
{
public:

    /** Default constructor                                                                                 */
    IGMTransmission() { };
    
    
    /** Constructor (if no IGM)        
        @param lmin         minimum wavelength for interpolation grid
        @param lmax         maximum wavelength for interpolation grid
        @param npt          number of points in interpolation grid                                          */
    IGMTransmission(double lmin, double lmax, int npt) {
    
        vector<double> rflam;
        vector<double> trans;
        double dl = (lmax-lmin)/(npt-1.);
        for (int i=0; i<npt; i++) {
            double lam = lmin + i*dl;
            double t = 1.;
            rflam.push_back(lam);
            trans.push_back(t);
            }
        trans_.DefinePoints(rflam, trans, rflam[0], rflam[rflam.size()-1], npt);
        };
    
    /** Constructor if reading transmission data from a file (could be mean IGM transmission or a single line
        of sight)
        @param filename     file containing IGM transmission data
        @param lmin         minimum wavelength for interpolation grid
        @param lmax         maximum wavelength for interpolation grid
        @param npt          number of points in interpolation grid
        @param nComments    number of comments at top of IGM transmission data file                         */
    IGMTransmission(string& filename, double lmin, double lmax, int npt=1024, int nComments=1) {
        cout <<"     About to read from file "<< filename << endl;
        readTransmission(filename, lmin, lmax, npt, nComments);
        };
        
    /** Constructor if using Madau law for IGM
        @param zSED    redshift of galaxy SED                                 
        @param isLyC   include Lyman continuum absorption                                                   */
    IGMTransmission(double zSED, bool isLyC=true, double lmin=1e-8, double lmax=1e-5, int npt=1000) { 
        setupMadau(zSED, isLyC, lmin, lmax, npt);
        };
        
        
    /** Return IGM transmission CHECK THIS METHOD WORKS!
        @param lambda   restframe wavelength                                  */
    virtual double operator() (double lambda) const { return trans_(lambda); };
    
    
    /** 
        @param filename     file containing IGM transmission data
        @param lmin         minimum wavelength for interpolation grid
        @param lmax         maximum wavelength for interpolation grid
        @param npt          number of points in interpolation grid
        @param nComments    number of comments at top of IGM transmission data file                         */
    void readTransmission(string& filename, double lmin, double lmax, int npt=1024, int nComments=1) {
        cout <<"     About to read from file "<< filename << endl;
        trans_.ReadXYFromFile(filename, lmin, lmax, npt, nComments, false);  
        };
    

protected:
    void setupMadau(double zSED, bool isLyC, double lmin, double lmax, int npt) {
        
        Madau madau(5, isLyC); // 5 refers to max number of Lyman lines to calculate
                               // with Madau implementation, 5 is the maximum possible
        
        vector<double> rflam;
        vector<double> trans;
        double dl = (lmax-lmin)/(npt-1.);
        for (int i=0; i<npt; i++) {
            double lam = lmin + i*dl;
            double t = madau.returnObserverFrameTransmission(lam, zSED); //returnRestFrameTransmission(lam, zSED);
            // ^^ has to be observer frame to match the fact that the IGMTransmission files are read in as a
            //    functionof observed wavelength
            rflam.push_back(lam);
            trans.push_back(t);
            }
        trans_.DefinePoints(rflam, trans, rflam[0], rflam[rflam.size()-1], npt);
        };
        
protected:
    SInterp1D trans_;             /**< Transmission look up table                                           */
};


#endif
