/**
 * @file  igm.h
 * @brief Contains classes that provide functions to Monte Carlo simulate the IGM
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

#include "machdefs.h"
#include <math.h>
#include <iostream>
#include <fstream>


// sophya
#include "genericfunc.h"
#include "sopnamsp.h"
#include "mydefrg.h"
//#include "pexceptions.h"

// CatSim
#include "constcosmo.h"
#include "geneutils.h"
#include "sinterp.h"


/** AtomicCalcs
  *
  * Class that holds basic atomic calculations and constants mostly to do with
  * the Lyman series
  * 
  *
  */
class AtomicCalcs
{
public:
    /** Constructor */
    AtomicCalcs(){ setGammas(); setConstants(); nLineMaxMax_ =  31; };

    /** Return Doppler width in meters with line center wavelength corresponding
        to Lyman series n \f$ \Delta\lambda=\lambda_i\frac{b}{c} \f$
        @param nLine starting energy level of the Lyman line transition
        @param dopplerParKMS doppler parameter in km/s                        */
    double returnDopplerWidthWL(int nLine, double dopplerParKMS)
        { return returnWavelengthLymanSeries(nLine)*(dopplerParKMS/SPEED_OF_LIGHT_KMS); };
        //THIS SEEMS RIGHT BUT DOESN'T SEEM TO WORK? 
        //{ return returnWavelengthLymanSeries(nLine)*(SPEED_OF_LIGHT_KMS/dopplerParKMS); };
        
    /** Return Doppler width in s^-1 with line center frequency corresponding to 
        Lyman series n \f$ \Delta\nu=\nu_i\frac{b}{c} \f$
        @param n starting energy level of the Lyman line transition
        @param dopplerParKMS doppler parameter in km/s                        */
    double returnDopplerWidthFreq(int nLine, double dopplerParKMS)
        { return returnFrequencyLymanSeries(nLine)*(dopplerParKMS/SPEED_OF_LIGHT_KMS); }; 
        //return SPEED_OF_LIGHT_KMS/returnDopplerWidthWL(nLine,dopplerParKMS); }; 
        //returnFrequencyLymanSeries(nLine)*(dopplerParKMS/SPEED_OF_LIGHT_KMS); };
        
    /** Return Lyman series wavelength in meters 
        @param  n is the starting energy level of the Lyman line transition
                (therefore n can have a minimum of 2: Lyman-alpha)            */
    double returnWavelengthLymanSeries(int n)
        { return WAVE_LYMANLIM_METERS/(1.-1./(n*n)); };
        
    /** Return Lyman series wavelength in angstroms 
        @param  n is the starting energy level of the Lyman line transition
                (therefore n can have a minimum of 2: Lyman-alpha)            */
    double returnWavelengthAngstrLymanSeries(int n)
        { return WAVE_LYMANLIM_ANGSTR/(1.-1./(n*n)); };
        
    /** Return Lyman series frequency in s^-1
        @param  n is the starting energy level of the Lyman line transition 
                (therefore n can have a minimum of 2: Lyman-alpha)            */
    double returnFrequencyLymanSeries(int n)
        { return SPEED_OF_LIGHT_MS/returnWavelengthLymanSeries(n); };
        
   /** Return the (unitless) wavelength difference relative to resonant wavelength
       in Doppler units
       @param lambda wavelength in meters */
    double returnX(double lambda, int nLine, double dopplerPar) {
        double dl = (lambda - returnWavelengthLymanSeries(nLine));
        double dw = returnDopplerWidthWL(nLine,dopplerPar);
        return dl/dw; };
    
    /** Return (unitless) damping parameter a */            
    double returnDampingParameter(int nLine, double dopplerPar);
        
    /** Return oscillator strength of Lyman series line nLine                 */
    double returnOscillatorStrength(int nLine) {
        int n = nLine - 2;
        return fSeries_[n]; };
  
    // Methods that set atomic constants    
    
    /** Set damping constant gamma for the Lyman series units of s^-1         */
    void setGammas();
    
    /** Set constants used in this class such as freq. of Lyman series        */
    void setConstants();
    
protected: 
    vector<double> gammaSeries_;    /**< damping constant of Lyman series s^-1*/
    vector<double> fSeries_;        /**< oscillator strength of Lyman series  */
    double freqLymanLimInvSecs_;    /**< frequency of the Lyman limit in s^-1 */
    double sigLymanLimCM2_;         /**< cross-section at the Lyman limit in cm^2*/
    int nLymanAlpha_;               /**< starting level of Lyman-a transition*/
    int nLineMaxMax_;               /**< maximum Lyman series line possible to use */
    
};

/** HIColumnDensity class
  *
  * Class to calculate the HI column density distribution function
  * Parameterization based on Inoue & Iwata 2008 eqn 4 g(N_HI)
  *
  */
class HIColumnDensity :
    public GenericFunc
{
public:
    /** Constructor 
        @param beta1 First power law index (default = 1.6)
        @param beta2 Second power law index (default = 1.3)
        @param Nc Column density value at break (default = 1.6e17/cm^2)
        @param Nl Lower bound of column density values (default = 1e12/cm^2)
        @param Nu Upper bound of column density values (default = 1e22/cm^2)
        @param Nstep Number of steps to integrate each part of power law with
    */
    HIColumnDensity(double beta1=1.6,double beta2=1.3,double Nc=1.6e17,
                            double Nl=1e12, double Nu=1e22, int Nstep=500000);
    
    /** Set normalization such that int g(NHI) dNHI = 1 */
    double setBnorm(int Nstep);
    
    /** Analytical integration of \f$ (N_{HI}/N_c)^{-\beta} \f$               */
    double integratePowerLaw(double NHI, double beta)
        {  double npower = pow(NHI,(1.-beta));
           double constant = (1.-beta)*pow(Nc_,-beta);
           return npower/constant; };
    
    /** Return column density distribution at column density value given */
    virtual double operator()(double Nh)
            { return returnColDensityDist(Nh); };
            
    /** Return column density distribution at column density value given */
    double returnColDensityDist(double Nh);
    
    /** Return power law at column density value given */
    double returnFirstPowerLaw(double Nh)
            { return pow( (Nh/Nc_),-beta1_ ); };
            
    /** Return power law at column density value given */
    double returnSecondPowerLaw(double Nh)
            { return pow( (Nh/Nc_),-beta2_ ); };
            
    /** Return the power law index values */
    void returnPowerLaws(double& beta1, double& beta2)
            { beta1=beta1_; beta2=beta2_; };
            
    /** Return the column density limits */
    void returnLowerUpperColDensValues(double& Nl, double& Nu)
            { Nl=Nl_; Nu=Nu_; };
            
    /** Return column density at power law break */
    double returnColDensAtBreak() { return Nc_; };
    
    /** Return the normalization of the power law */
    double returnPowerLawNormalization() { return Bnorm_; };
    
    /** Write the column density distribution to a file, log-spaced
        @param outfile Name of file to write to
        @param dLog    Log of step in column density values
        @param nStep   Number of column density values
     */
    void writeToFile(string outfile, double dLog, int nStep);
    
    
    /** Numerically integrate power law (incorrectly normalized) */
    double numIntegratePowerLaw(int Nstep);
    
    /** Check integration of column density function is 1 */
    double checkIntegration(int Nstep);
    

protected:
    double beta1_;          /**< Power law index */
    double beta2_;          /**< Power law index */
    double Nc_;             /**< Column density at power law break */
    double Nl_;             /**< Lower column density values */
    double Nu_;             /**< Upper column density values */
    double Bnorm_;          /**< Normalization of column density function*/

};

/** AbsorberRedshiftDistribution class
  *
  * Class to calculate the absorber (LAF, LLS, DLA) redshift distribution
  * Parameterization based on Inoue & Iwata 2008 eqn 5 f(z)
  *
  */
class AbsorberRedshiftDistribution :
    public GenericFunc
{
public:
    /** Constructor 
        @param gamma1 First power law index (default = 0.2)
        @param gamma2 Second power law index (default = 2.5)
        @param gamma3 Third power law index (default = 4.0)
        @param z1 Redshift value at first break (default = 1.2)
        @param z2 Redshift value at second break (default = 4.0)
        @param A (Normalization) Total number of absorbers at z=z1 with a 
               column density 1e12<=NHI<=1e22 /cm^2 (default = 400)
    */
    AbsorberRedshiftDistribution(double gamma1=0.2,double gamma2=2.5,
                 double gamma3=4.0, double z1=1.2, double z2=4.0, double A=400);
    
    /** Return redshift distribution at redshift value given */
    virtual double operator()(double z)
            { return returnRedshiftDist(z); };
            
    /** Return redshift distribution at z value given */
    double returnRedshiftDist(double z);
    
    /** Return power law at z value given */
    double returnFirstPowerLaw(double z)
            { return pow( ((1.+z)/(1.+z1_)),gamma1_ ); };
    /** Return power law at z given */
    double returnSecondPowerLaw(double z)
            { return pow( ((1.+z)/(1.+z1_)),gamma2_ ); };
    /** Return power law at z value given */
    double returnThirdPowerLaw(double z){
            double parta = pow( ((1.+z2_)/(1.+z1_)),gamma2_ );
            double partb = pow( ((1.+z)/(1.+z2_)),gamma3_ );
            return parta*partb; };


    /** Return the power laws */
    void returnPowerLaws(double& gamma1, double& gamma2, double& gamma3)
            { gamma1=gamma1_; gamma2=gamma2_; gamma3=gamma3_;};
    /** Return the redshifts at the breaks */
    void returnzAtBreaks(double& z1, double& z2)
            { z1=z1_; z2=z2_; };
    /** Return the normalization */
    double returnNormalization() { return A_; };
    
    /** Write the absorber redshift distribution to a file, 
        @param outfile Name of file to write to
        @param zmin    Minimum redshift
        @param dz      Step in redshifts
        @param nStep   Number of redshift values
     */
    void writeToFile(string outfile, double zmin, double dz, int nStep);
            
protected:
    double gamma1_;    /**< Power law index */
    double gamma2_;    /**< Power law index */
    double gamma3_;    /**< Power law index */
    double z1_;        /**< Redshift value at first break */
    double z2_;        /**< Redshift value at second break */
    double A_;         /**< (Normalization) Total number of absorbers at z=z1 
                        with a column density 1e12<=NHI<=1e22 /cm^2 */
};


/** DopplerParDistribution class
  *
  * Class to calculate the Doppler parameter b distribution
  * Parameterization based on Inoue & Iwata 2008, and Hui & Rutledge 1999
  *
  */
class DopplerParDistribution :
    public GenericFunc
{
public:
    /** Constructor 
        @param bsigma Single parameter of distribution (default = 23km/s)
    */
    DopplerParDistribution(double bsigma=23)
    : bsigma_(bsigma) { };
    
    /** Return Doppler distribution at b value given */
    virtual double operator()(double b)
            { return returnDopplerDist(b); };
            
    /** Return Doppler distribution at b value given */
    double returnDopplerDist(double b);
    
    /** Write the doppler parameter distribution to a file, 
        @param outfile Name of file to write to
        @param bmin    Minimum doppler parameter
        @param db      Step in doppler parameter
        @param nStep   Number of doppler parameter values
     */
    void writeToFile(string outfile, double zmin, double db, int nStep);
    
protected:
    double bsigma_; /**< Single parameter of distribution */
};

/** ProbabilityDistAbsorbers class
  *
  * Class to calculate the probability of an absorber within delta z range of
  * z', if there is already an absorber at z'
  * Parameterization based on Inoue & Iwata 2008
  *
  */
class ProbabilityDistAbsorbers //:
    //public GenericFunc
{
public:
    /** Constructor
        @param rg Class that mediates random number drawing
        @param absorberZDist Class holding the absorber redshift distribution
        @param hiColumnDensity
        @param dopplerParDist */
    ProbabilityDistAbsorbers(RandomGeneratorInterface& rg,
                            AbsorberRedshiftDistribution& absorberZDist,
                            HIColumnDensity& hiColumnDensity,
                            DopplerParDistribution& dopplerParDist
        );
    
    /** Simulate a line of sight distribution of absorbers, returns number of 
        absorbers
        @param zStart           Starting redshift of line of sight distribution
        @param zMax             Max redshift along line of sight
        @param redshifts        Vector of absorber redshifts (sorted in ascending order)
        @param dopplerPars      Vector of absorber doppler parameters 
        @param columnDensity    Vector of absorber column densities **/
    int simulateLineOfSight(double zStart,double zMax, 
                    vector<double>& redshifts, vector<double>& dopplerPars,
                               vector<double>& columnDensities, string outfile);
            
    void simulateAbsorber(double zCurrent, double& redshift, double& dopplerPar,
                                                    double& columnDensity);
    /** Draw deltaZ of next absorber (next absorber is at zLast+deltaZ).  This 
        uses the "Inverse Transformation method"
        @param zLast Redshift of last absorber                                */   
    double drawDeltaZ(double zLast);
    
    /** Draw absorber redshift given redshift of source
        @param zSource  redshift source                                       */
    double simulateAbsorberRedshift(double zSource);
    
    /** Draw HI column density of absorber.  This uses the "Rejection method" */   
    double drawHIColumnDensity();
    
    /** Draw Doppler parameter of absorber.  This uses the "Rejection method" */   
    double drawDoppler();
    
    
    /** Set the HI column density distribution to draw from */
    void setNHiDistribution(int nStep=1000);
    /** Set the doppler parameter distribution to draw from */
    void setDopplerDistribution(int nStep=1000);
    
    /** Set HI column density distribution with a new step size */
    void resetNHiDistribution(int nStep=1000)
        { setNHiDistribution(nStep); };
        
    /** Set Doppler distribution with a new min, max and step size */
    void resetDopplerDistribution(double bmin=0,double bmax=200,int nStep=1000)
        { bmin_=bmin; bmax_=bmax; 
        setDopplerDistribution(nStep); };
        
        
    /** write simulated absorbers redshifts, doppler parameters and column
        densities to a file */
    void writeToFile(string outfile, vector<double>& redshifts, 
                    vector<double>& dopplerPars, vector<double>& columnDensity);        
                    
        
    // These write**Distribution methods are really for debugging/checking
        
    /** Write redshifts to a file 
        @param outfile  File to write to
        @param zCurrent redshift of current absorber
        @param nTrial   Number of redshifts to draw
    */
    void writeZDistribution(string outfile, double zCurrent, long nTrial);
    /** Write doppler parameters to a file 
        @param outfile  File to write to
        @param nTrial   Number of doppler parameters to draw
    */
    void writeDopplerDistribution(string outfile, long nTrial);
    /** Write column densities to a file 
        @param outfile  File to write to
        @param nTrial   Number of column densities to draw
    */
    void writeNHiDistribution(string outfile, long nTrial);
    
    /** Return the grid used to draw the HI column densities from             */
    void returnColumnDensityGrid(vector<double>& log10NHIvals, vector<double>& log10gvals)
        {
            log10NHIvals = log10NHIvals_;
            log10gvals = log10gvals_;
            return;
        };
    
protected:
    RandomGeneratorInterface&       rg_; /**< For random number generation */
    /** Class holding the absorber redshift distribution */
    AbsorberRedshiftDistribution&   absorberZDist_;
    /** Class holding the HI column density distribution */
    HIColumnDensity&                hiColumnDensity_;
    /** Class holding the Doppler parameter distribution */
    DopplerParDistribution&         dopplerParDist_;
    double log10Nl_;                /**< log10 of minimum column density */
    double log10Nu_;                /**< log10 of maximum column density */
    vector<double> log10NHIvals_;   /**< log10 column density grid values */
    vector<double> log10gvals_;     /**< log10 column density grid values dist*/
    double log10gmin_;              /**< log10 of minimum of column density dist*/
    double log10gmax_;              /**< log10 of maximum of column density dist*/
    SInterp1D   colDensityFunc_;    /**< Interpolation function for column density dist*/
    vector<double> bvals_;          /**< doppler parameter b grid values */
    vector<double> hvals_;          /**< doppler parameter b grid values dist*/
    double bmin_;                   /**< minimum doppler parameter value*/
    double bmax_;                   /**< maximum doppler parameter value*/
    double hmin_;                   /**< minimum doppler parameter dist value*/
    double hmax_;                   /**< maximum doppler parameter dist value*/

};


/** OpticalDepth class
  *
  * Class to calculate the optical depth due to a single absorber at redshift z
  * with column density NHI and doppler parameter b
  *
  */
class OpticalDepth : public AtomicCalcs
{
public:
    /** Constructor
        @param zAbsorber redshift of absorber
        @param nHI column density
        @param bAbsorber doppler parameter
        @param returnOpticalDepth return optical depth or transmission */
    //OpticalDepth(double zAbsorber, double nHI, double bAbsorber, 
    //                                        bool isOpticalDepth=false)
        //: zAbsorber_(zAbsorber) , nHI_(nHI) , bAbsorber_(bAbsorber) , 
        //                        isOpticalDepth_(isOpticalDepth) { 
    OpticalDepth() {
    nLineMax_ = nLineMaxMax_; // Maximum line series to actually go up to
                              // nLineMaxMax is set in AtomicCalcs
    setLymanAll(); // by default do both contributions
    //setConstants(); push to atomic calcs class
    };
    
    /** Return optical depth or transmission at wavelength value given 
        @param lambda observed wavelength in meters
    //virtual double operator()(double lambda)
    virtual double operator()(double lambda, double zAbsorber, double nHI, double bAbsorber)
            {   double lambdaE = lambda/(1. + zAbsorber); // convert to emission frame
                double freq = SPEED_OF_LIGHT_MS/lambdaE;   // convert to frequency
                double tau = returnRestFrameOpticalDepth(freq, bAbsorber, nHI); 
                if (isOpticalDepth_)
                    return tau;
                else
                    return exp(-tau);
                    
            };*/
            
    /** Return the transmission in the observer's frame due to an absorber 
        with redshift z, column density nHI, and doppler parameter bAbsorber
        @param lambda       observed wavelength in meters                          
        @param zAbsorber    redshift of absorber
        @param nHI          HI column density of absorber (cm^-2)
        @param bAbsorber    doppler parameter of absorber (km/s)              */
    double returnObserverFrameTransmission(double lambda, double zAbsorber, double nHI, double bAbsorber)
        { double tau = returnObserverFrameOpticalDepth(lambda, zAbsorber, nHI, bAbsorber);
          return exp(-tau); };    
    
    /** Return the optical depth in the observer's frame due to an absorber 
        with redshift z, column density nHI, and doppler parameter bAbsorber
        @param lambda       observed wavelength in meters                          
        @param zAbsorber    redshift of absorber
        @param nHI          HI column density of absorber (cm^-2)
        @param bAbsorber    doppler parameter of absorber (km/s)              */
    double returnObserverFrameOpticalDepth(double lambda, double zAbsorber, double nHI, double bAbsorber)
        { double lambdaE = lambda/(1. + zAbsorber); // convert to absorption frame
          double freq = SPEED_OF_LIGHT_MS/lambdaE;  // convert to frequency
          double tau = returnRestFrameOpticalDepth(freq, bAbsorber, nHI);
          return tau; };            
    
    /** Return (unitless) optical depth in the rest-frame of the absorber
        Includes both Lyman-Limit Systems and Lyman-Alpha Forest
        @param freq         rest-frame frequency in s^-1                            
        @param bAbsorber    doppler parameter of absorber (km/s)              
        @param nHI          column density of HI in cm^-2                     */
    double returnRestFrameOpticalDepth(double freq, double bAbsorber, double nHI);
    
    /** Return the Lyman continuum cross-section in cm^2
        @param freq         rest-frame frequency in s^-1                           */
    double returnLymanContinuumCrossSection(double freq);
    
    /** Return the cross-section for the total Lyman-series in cm^2
        @param freq         rest-frame frequency in s^-1     
        @param bAbsorber    doppler parameter of absorber (km/s)              */
    double returnLymanSeriesCrossSection(double freq, double bAbsorber);
    
    /** Return the cross-section for the ith Lyman-series in cm^2
        @param nLine        starting energy level of the Lyman line transition
        @param freq         rest-frame frequency in s^-1      
        @param bAbsorber    doppler parameter of absorber (km/s)              */
    double returnLymanLineCrossSection(int nLine, double freq, double bAbsorber);
    
    /** Return the (unitless) optical depth for just the Lyman continuum contribution 
        @param freq         rest-frame frequency in s^-1 
        @param nHI          column density in cm^-2                           */
    double returnLymanContinuumOpticalDepth(double freq, double nHI) {
        return nHI*returnLymanContinuumCrossSection(freq); };
        
    /** Return the (unitless) optical depth for just the Lyman series contribution    
        @param freq         rest-frame frequency in s^-1 
        @param nHI          column density in cm^-2                           
        @param bAbsorber    doppler parameter in km/s                         */
    double returnLymanSeriesOpticalDepth(double freq, double nHI, double bAbsorber) {
        return nHI*returnLymanSeriesCrossSection(freq, bAbsorber); };
        
    /** Return the (unitless) line profile function, analytical approx by 
        Tepper-Garica 2006 
        @param nLine        starting energy level of the Lyman line transition
        @param freq         rest-frame frequency in s^-1   
        @param bAbsorber    doppler parameter  in km/s                        */
    double returnLineProf(int nLine, double freq, double bAbsorber);
        
    
        
    /** Set contribution type as both Lyman contributions                     */
    void setLymanAll() { setContribution(0); };    
    /** Set contribution type as Lyman continuum only                         */
    void setLymanSeriesOnly() { setContribution(1); };
    /** Set contribution type as Lyman series only                            */
    void setLymanContinuumOnly() { setContribution(2); };
    
                
    /** Set contribution type: Lyman continuum only, Lyman series only or both     
        @param contributionType 0=both, 1=Lyman series, 2=Ly*/
    void setContribution(int contributionType)
        { if (contributionType==0)
            { isAll_=true; isOnlyLymanC_=false; isOnlyLymanS_=false; }
          else if (contributionType==1)
            { isAll_=false; isOnlyLymanC_=false; isOnlyLymanS_=true; }
          else if (contributionType==2)
            { isAll_=false; isOnlyLymanC_=true; isOnlyLymanS_=false; }
          else
            throw ParmError("ERROR! contribution option not understood");
        };
        
    /** Set whether to return the optical depth or the transmission           
    void setReturnType(bool isOpticalDepth)
        { isOpticalDepth_ = isOpticalDepth; };*/
        
    /** Set Lyman line calculation limit */
    void setnLineMax(int nLineMax)
                {   nLineMax_ = nLineMax;
                    if (nLineMax > nLineMaxMax_) 
                        throw ParmError("ERROR! Cannot compute Lyman series this far!");
                };
                
    // generic stuff
    
    /** Return Doppler width in s^-1 with line center frequency corresponding to 
        Lyman series n 
        @param n starting energy level of the Lyman line transition
        @param dopplerParKMS doppler parameter in km/s                        */
    //double returnDopplerWidth(int nLine, double dopplerParKMS)
   //     { return returnFrequencyLymanSeries(nLine)*(dopplerParKMS/SPEED_OF_LIGHT_KMS); };
    
    /** Return oscillator strength of Lyman series line nLine                 */
    //double returnOscillatorStrength(int nLine) {
    //    int n = nLine - 2;
    //    return fSeries_[n]; };
                                                                   
    /** Return Lyman series wavelength in meters 
        @param  n is the starting energy level of the Lyman line transition
                (therefore n can have a minimum of 2: Lyman-alpha)            */
    //double returnWavelengthLymanSeries(int n)
    //    { return WAVE_LYMANLIM_METERS/(1.-1./(n*n)); };
        
    /** Return Lyman series frequency in s^-1
        @param  n is the starting energy level of the Lyman line transition 
                (therefore n can have a minimum of 2: Lyman-alpha)            */
   // double returnFrequencyLymanSeries(int n)
   //     { return SPEED_OF_LIGHT_MS/returnWavelengthLymanSeries(n); };
        
    /** Set constants used in this class such as freq. of Lyman series */
    //void setConstants();
    
    
                
protected:
    //double zAbsorber_;              /**< redshift of absorber                 */
    //double nHI_;                    /**< column density of absorber           */
    //double bAbsorber_;              /**< doppler parameter of absorber in km/s*/
    bool isAll_;                    /**< return both the Lyman contributions  */
    bool isOnlyLymanC_;             /**< return only the Lyman continuum contribution */
    bool isOnlyLymanS_;             /**< return only the Lyman series contribution    */
    //bool isOpticalDepth_;           /**< return optical depth or transmission */
    int nLineMax_;                  /**< maximum Lyman series line to use     */
    //double freqLymanLimInvSecs_;    /**< frequency of the Lyman limit in s^-1*/
    //double sigLymanLimCM2_;         /**< cross-section at the Lyman limit in cm^2*/
    //int nLymanAlpha_;               /**< starting level of Lyman-a transition*/
    //int nLineMaxMax_;               /**< maximum Lyman series line possible to use */
    //vector<double> fSeries_;        /**< oscillator strength of ith Lyman line*/
};


/** LineOfSightTrans class
  *
  * Calculates the line of sight transmission after attenuation by the IGM
  * absorber distribution as a function of observed wavelength and redshift of 
  * the source
  *
  */
class LineOfSightTrans : public GenericFunc, public OpticalDepth
{
public:
    /** Constructor
        @param redshifts of absorbers along line of sight (must be sorted in ascending order)
        @param doppler parameters of absorbers along line of sight
        @param column densities of absorbers along line of sight              */
    LineOfSightTrans(vector<double> redshifts, vector<double> dopplerPars,
                 vector<double> columnDensities, bool isOpticalDepth=false)
    : redshifts_(redshifts) , dopplerPars_(dopplerPars) , 
      columnDensities_(columnDensities) , isOpticalDepth_(isOpticalDepth) {
        bool isSorted = sortCheck(redshifts_);      
        if (!isSorted)
            throw ParmError("ERROR! redshifts vector is not sorted in ascending order");
        //setContribution(0); // by default use both contributions
        setLymanAll(); // by default use both contributions
        if (isOpticalDepth)
            cout << "     Returning optical depth"<<endl;
        else
            cout << "     Returning transmission"<<endl;
       };
    
    /** Return total optical depth (\f$\tau\f$) or transmission (\f$e^{-\tau}\f$)
        at observed wavelength value given due to all absorbers along line of sight 
        @param lambda   observed wavelength in meters                           
        @param zSource  redshift of background source                         */
    virtual double operator()(double lambda, double zSource)
        {   if (isOpticalDepth_)
                return returnOpticaldepth(lambda, zSource);
            else 
                return returnTransmission(lambda, zSource);
        };
        
    /** Return optical depth at observed wavelength lambda in meters of light
        from a source located at zSource, returns \f$e^{-\tau}\f$
        @param lambda   observed wavelength
        @param zSource  redshift of background source                         */
    double returnOpticaldepth(double lambda, double zSource);
    
    /** Return transmission at observed wavelength lambda in meters of light
        from a source located at zSource, returns \f$\tau\f$
        @param lambda observed wavelength
        @param zSource redshift of background source                          */
    double returnTransmission(double lambda, double zSource) {
        double tau = returnOpticaldepth(lambda, zSource);
        return exp(-tau);
        };
    
    /** Set whether to return the optical depth or the transmission           */
    void setReturnType(bool isOpticalDepth)
        { isOpticalDepth_ = isOpticalDepth; };
         
    /** Set contribution type: Lyman continuum only, Lyman series only or both     
        @param contributionType 0=both, 1=Lyman series, 2=Ly
    void setContribution(int contributionType)
        { if (contributionType==0)
            { isAll_=true; isOnlyLymanC_=false; isOnlyLymanS_=false; }
          else if (contributionType==1)
            { isAll_=false; isOnlyLymanC_=false; isOnlyLymanS_=true; }
          else if (contributionType==2)
            { isAll_=false; isOnlyLymanC_=true; isOnlyLymanS_=false; }
          else
            throw ParmError("ERROR! contribution option not understood");
        };*/
    
protected:
    vector<double> redshifts_;      /**< redshift of absorbers                */
    vector<double> dopplerPars_;    /**< doppler parameters of absorbers      */
    vector<double> columnDensities_;/**< column densities of absorbers        */
    bool isOpticalDepth_;           /**< return optical depth or transmission */
    //bool isOnlyLymanC_;             /**< return only the Lyman continuum contribution */
    //bool isOnlyLymanS_;             /**< return only the Lyman series contribution    */
    //bool isAll_;                    /**< return both the Lyman contributions  */
    //int nLineMax_;                  /**< maximum Lyman series line to use     */
};



/** VoigtProfile class
  *
  * Class to calculate profiles of HI absorption lines imprinted on the spectra 
  * of bright background sources by intervening absorbing systems. Uses a simple
  * analytical approximation involving the Voigt-Hjerting function and an 
  * absorption coefficient.  Follows Tepper-Garcia 2006
  *
  * Voigt profiles play an important role in the spectroscopy of stellar 
  * atmospheres where accurate measurements of line wings allow the 
  * contributions of Doppler broadening, or natural linewidth of collision line 
  * broadening to be separated. From this measurements the temperature and 
  * pressure of the emitting or absorbing layers in the stellar atmospheres can 
  * be determined.
  *
  */
class VoigtProfile :
    public GenericFunc, public AtomicCalcs
{
public:
    /** Constructor
        @param dopplerPar doppler parameter of absorber
        @param Lyman series line number (Lyman-alpha = 2, Lyman-beta=3 ...) */
    VoigtProfile(double dopplerPar, int nLine)
    : dopplerPar_(dopplerPar) , nLine_(nLine) {
        //setGammas();// set the damping constant gamma for up to Ly-series 32
        int nLinesWithGamma = gammaSeries_.size() + 1;
        if ( nLine_>nLinesWithGamma ) {
            string emsg = "ERROR! Cannot calculate profile of Lyman series n>";
            stringstream ss; ss<<(gammaSeries_.size()+1);
            emsg += ss.str();
            throw ParmError(emsg);
            }
            
        };
    
    /** Returns the (unitless) line profile at lambda, see Tepper-Garcia 2006 equation 25 
        @param lambda the wavelength in meters */
    virtual double operator()(double lambda) {
        double x = returnX(lambda, nLine_, dopplerPar_);
        double a = returnDampingParameter(nLine_, dopplerPar_);
        double K = kFunction(x);
        double H1 = exp(-x*x)*(1-a*(2./sqrt(PI))*K);
        //cout <<" x="<<x<<", a="<<a<<", K(x)="<<K<<", H1="<<H1<<endl;
        return H1; };
                
    
    /** Returns the (unitless) function defined in Tepper-Garcia 2006 equation 24 
        @param x the wavelength difference relative to resonant wavelength
               in Doppler units */
    double kFunction(double x);
    
    
        
    
    // THESE FUNCS BELOW ARE GENERIC: SHOULD THEY BE MOVED OUT SOMEWHERE?
    // probably should make a new class called AtomicCalcs and then have
    // classes inherit from that class
    
    /** Return the (unitless) wavelength difference relative to resonant wavelength
        in Doppler units
        @param lambda wavelength in meters */
    //double returnX(double lambda) {
    //    double dl = (lambda - returnWavelengthLymanSeries(nLine_));
    //    double dw = returnDopplerWidth(nLine_,dopplerPar_);
    //    return dl/dw; };
    
    /** Return (unitless) damping parameter a */            
    //double returnDampingParameter();
    
    /** Set the damping constant gamma, units s^-1. Hard-coded values for each Lyman 
        transmission as found in Table 2 of Morton 2003 and calculated from a 
        file containg values of the damping parameter a assuming b=36km/s from 
        Tepper-Garcia                                                         */
    // void setGamma();
    
    /** Return Doppler width in m with line center frequency corresponding to 
        Lyman series n 
        @param nLine starting energy level of the Lyman line transition
        @param dopplerParKMS doppler parameter in km/s 
        THIS SEEMS RIGHT BUT DOESN'T SEEM TO WORK? */
    //double returnDopplerWidth(int nLine, double dopplerParKMS)
    //    { return returnWavelengthLymanSeries(nLine_)*(SPEED_OF_LIGHT_KMS/dopplerParKMS); };
        
    /** Return Lyman series wavelength in meters 
        @param n is the starting energy level of the Lyman line transition */
    //double returnWavelengthLymanSeries(int n)
    //    { return WAVE_LYMANLIM_METERS/(1.-1./(n*n)); };
        
    /** Return Lyman series frequency in s^-1
        @param n is the starting energy level of the Lyman line transition */
    //double returnFrequencyLymanSeries(int n)
    //    { return SPEED_OF_LIGHT_MS/returnWavelengthLymanSeries(n); };
    
protected:
    double dopplerPar_;             /**< doppler parameter of absorber */
    int nLine_;                     /**< starting energy level of the Lyman line transition */
    //vector<double> gamma_;          /** damping constant */

};

/** Madau class
  *
  */
class Madau:
    public AtomicCalcs
{
public:
    /** Constructor
        @param nLineMax     Maximum Lyman-series to include                   */
    Madau(int nLineMax=5, bool isLyC=true)
    : nLineMax_(nLineMax) , isLyC_(isLyC)
        { setAbsorptionStrengths(); };
    
    /** Return the transmission in the observer's frame along a line of sight
        to a source at some redshift z
        @param lambda       observed wavelength in meters   
        @param zSource      redshift of source                  */
    double returnObserverFrameTransmission(double lambda, double zSource)
        { double tau = returnObserverFrameOpticalDepth(lambda, zSource);
          return exp(-tau); };
          
    /** Return the transmission in the rest-frame along a line of sight
        to a source at some redshift z
        @param lambda       rest-frame wavelength in meters   
        @param zSource      redshift of source                  */
    double returnRestFrameTransmission(double lambda, double zSource) { 
            double lambdaObs = lambda*(1. + zSource);
            double tau = returnObserverFrameOpticalDepth(lambdaObs, zSource);
            return exp(-tau); };
    
    /** Return the optical depth in the observer's frame along a line of sight
        to a source at redshift z
        @param lambda       observed wavelength in meters                          
        @param zSource      redshift of source                                */
    double returnObserverFrameOpticalDepth(double lambda, double zSource);
          
    /** Return the optical depth in observer's frame along a line of sight due
        to Lyman continuum absorption only                                    */
    double returnLymanContinuumOpticalDepth(double zAbsorber, double zSource); 
    
    /** Set constants of Lyman series absorption strengths                    */
    void setAbsorptionStrengths();



protected:
    int nLineMax_;          /**< Maximum Lyman series to include              */
    int nLineMaxMadau_;     /**< Maximum Lyman series can include             */
    vector<double> Avals_;  /**< Absorption strength of Lyman-alpha to delta  */
    bool isLyC_;            /**< Include Lyman continuum absorption           */


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
  * Reads in the transmission values for a particular line of sight as a 
  * function of the observed wavelength.
  *
  * File contains: observed wavelength (m)  transmission (exp[-tau])
  *
  */
class IGMTransmission : public SInterp1D
{
public:
    IGMTransmission(){};
    IGMTransmission(string& filename, int npt=1024, int nComments=1) {
        cout <<"     About to read from file "<< filename << endl;
        ReadXYFromFile(filename, 1., -1., npt, nComments, false);
        };

};


#endif
