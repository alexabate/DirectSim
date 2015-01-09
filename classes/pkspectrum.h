/**
  * @file  pkspectrum.h
  * @brief Calculates linear transfer functions and growth functions
  *
  * @todo Add non linear power spectrum calculation
  *
  * @note transfer functions, power spectra are outputed 
  * with k in units of Mpc^-1 not hMpc^-1
  *
  * @author Christophe Magneville, modified by AA
  *
  * Created on: 2008
  * @date 2008
  *
  */
#ifndef PKSPECTRUM_SEEN
#define PKSPECTRUM_SEEN

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

// sophya includes
#include "machdefs.h"
#include "pexceptions.h"
// Sophya update v2.3 June 2013 replaces genericfunc with classfunc
#include "classfunc.h" 
//#include "genericfunc.h"

// DirectSim includes
#include "constcosmo.h"
#include "geneutils.h"
#include "sinterp.h"

namespace SOPHYA {

/** InitialSpectrum class
  *
  * Generic class for inital power spectrum
  */
class InitialSpectrum : public ClassFunc1D {
public:
    /** Constructor */
    InitialSpectrum(void) {};
    /** Copy constructor */
    InitialSpectrum(InitialSpectrum& pkinf) {};
    /** Destructor */
    virtual ~InitialSpectrum(void) {};
};


/** InitialPowerLaw class
  *
  * Class to calculate the initial early universe power-law spectrum 
  */
class InitialPowerLaw : public InitialSpectrum {
public:

    /** Constructor for initial power-law spectum
        @param n    spectral index
        @param a    normalization                                             */
    InitialPowerLaw(double n, double a=1.): n_(n), A_(a) { };

    /** Copy constructor                                                      */
    InitialPowerLaw(InitialPowerLaw& pkinf) : n_(pkinf.n_), A_(pkinf.A_) { };

    /** Desctructor                                                           */
    virtual ~InitialPowerLaw(void) {};

    /** Return initial power law spectrum at wavenumber \f$k\f$ 
        @param k    wavenumber in units?                                      */
    virtual double operator() (double k) const { return A_ * pow(k,n_);};

    /** Set normalization of power spectrum                                   */
    void SetNorm(double a) {A_ = a;};

    /** Set index of power law of power spectrum                              */
    void SetSlope(double n) {n_ = n;};
    
protected:
  double n_;    /**< spectral index                                           */
  double A_;    /**< normalization                                            */
};


/** TransferFunction class 
  *
  * Generic transfer function class
  *
  */
class TransferFunction : public ClassFunc1D {
public:
    /** Constructor */
    TransferFunction(void) {};
    /** Destructor */
    virtual ~TransferFunction(void) {};
};


/** TransferEH class
  *
  * Class to calculate Eisenstein & Hu transfer function
  * k is in units of 1/Mpc (NOT h/Mpc)
  *
  * Eisenstein & Hu ApJ 496:605-614 1998 April 1 (astro-ph/9709112)
  *
  */
class TransferEH : public TransferFunction {
public:

    /** Define ReturnPart type to describe what part of transfer function is 
        returned                                                              */
    typedef enum{ALL=0, CDM=1, BARYON=2} ReturnPart;

	/** Constructor 
        @param h100	        Hubble parameter in units of 100km/s/Mpc (H0 = 100h km/s/Mpc)
        @param OmegaCDM0    cold dark matter density today
        @param OmegaBaryon0 baryon density today
        @param tcmb         CMB temperature today
        @param nobaryon     if true don't compute effect of baryons
        @param lp           print level, if lp>0 prints extra to the screen   */
    TransferEH(double h100, double OmegaCDM0, double OmegaBaryon0, double tcmb, 
                                                  bool nobaryon=false, int lp=0);
  	
  	/** Copy constructor */			
    TransferEH(TransferEH& tf);
    
    /** Destructor */
    virtual ~TransferEH(void) {};
    
    /** Set/reset cosmological parameters
        @param h100         Hubble parameter in units of 100km/s/Mpc (H0 = 100h km/s/Mpc)
        @param OmegaCDM0    cold dark matter density today
        @param OmegaBaryon0 baryon density today                              */
    bool SetParTo(double h100,double OmegaCDM0,double OmegaBaryon0);
  
    /** Return transfer function at wavenumber \f$k\f$ 
        @param k    wavenumber in units?                                      */
    virtual double operator() (double k) const;
    
    /** Return position of first acoustic peak                                */
    double KPeak(void);
  
    /** Set form of transfer function in oscillatory region (BAO region)
        @param nooscenv    nooscenv = 0 use the baryon oscillatory part of transfer 
                           function (FULL TF). nooscenv = 1 use approx. in 
                           paragraph 3.3 p610 (middle of right column). Replace 
                           \f$ j_0(k\tilde{s})  ->  [1+(k*\tilde{s})^4]^(-1/4)\f$.
                           nooscenv = 2 use formulae 29+30+31 page 612 
                           [NO WIGGLES TF]                                    */
    void SetNoOscEnv(unsigned short nooscenv=0);
    
    /** Set part of transfer function to return: only baryon or only CDM or both
        @warning only relevant for nobaryon_=false AND nooscenv!=2 (not no wiggles)
   	    @param retpart	retpart = ALL, return full transfer function, 
   	                    retpart = CDM, return only CDM part, retpart = BARYON,
   	                    return only baryon part                               */
    void SetReturnPart(ReturnPart retpart=ALL) { retpart_ = retpart; };
    
    /** Set the printing level 
    	@param lp	if lp>0 prints extra to the screen                        */
    void SetPrintLevel(int lp=0) { lp_ = lp; };
    
    /** Return Hubble parameter in units of 100km/s/Mpc (H0 = 100h km/s/Mpc)  */
    double Returnh100() { return h100_; };
    
    /** Return matter density today                                           */
    double ReturnOmegaM() { return O0_; };

protected:
    int lp_;            /**< print level if lp>0 prints extra to the screen   */
    double O0_;         /**< matter density today                             */
    double Oc_;         /**< cold dark matter density today                   */
    double Ob_;         /**< baryon density today                             */
    double h100_;       /**< Hubble parameter in units of 100km/s/Mpc         */
    double tcmb_;       /**< CMB temperature in K                             */
    double th2p7_;      /**< CMB temperature normalized to 1 at 2.7K          */
    double zeq_;        /**< redshift of matter-radiation equality            */
    double keq_;        /**< scale of particle horizon at zeq                 */
    double zd_;         /**< redshift of photon decoupling (drag epoch)       */
    double Req_;        /**< ratio of baryon to photon momentum density at zeq*/
    double Rd_;         /**< ratio of baryon to photon momentum density at zd */
    double s_;          /**< sound horizon at the drag epoch                  */
    double ksilk_;      /**< Silk damping scale                               */
    double alphac_;     /**< transfer function parameter (CDM part)           */
    double betac_;      /**< transfer function parameter (CDM part)           */
    double bnode_;      /**< for k*s<bnode shift sinosoidal nodes to higher k */
    double alphab_;     /**< transfer function parameter (baryon part)        */
    double betab_;      /**< transfer function parameter (baryon part)        */
    double alphag_;     /**< used to rescale gamma in approximate TF          */
    double sfit_;       /**< approximate value of the sound horizon           */
    double kpeak_;      /**< position of first acoustic peak                  */
    bool nobaryon_;     /**< if true effect of baryons not included           */
    unsigned short nooscenv_; /**< form of TF: full, approx or no wiggles     */
    ReturnPart retpart_;/**< return either CDM, baryon or both parts of TF    */

	/** Calculate \f$\tilde{T}_0\f$ function                                  */
    double T0tild(double k,double alphac,double betac) const;
    
    /** Initialize class transfer function parameters                         */
    void Init_(void);
    
    /** Set class transfer function parameters to zero                        */
    void zero_(void);
};


/** TransferTabulate class
  *
  * For reading transfer function from a CAMB or CMBfast output file
  *
  * CAMB/CMBFast are openly available software packages that perform the full
  * numerical solution of the Bolzmann equations to calculate the transfer 
  * function
  *
  */
class TransferTabulate : public TransferFunction {
public:

    /** Constructor */
    TransferTabulate(void);
    
    /** Copy constructor */
    TransferTabulate(TransferTabulate& tf);
  
  	/** Destructor */
    virtual ~TransferTabulate(void) {};
    
    /** Return transfer function interpolated at wavenumber \f$k\f$ from values
        read from file 
        @param k    wavenumber in units?                                      */
    virtual double operator() (double k) const 
        { return InterpTab(k, k_, tf_, interptyp_); };
    
    /** Return number of wavenumber values read from file                     */
    int NPoints(void) {return k_.size();}
    
    /** Set type of interpolation (see InterpTab in geneutils)                */
    void SetInterpTyp(int typ=0);
    
    /** Read file output from CMBFast 
    	@param filename     CMBFast filename
    	@param h100         Hubble parameter in units of 100km/s/Mpc
    	@param OmegaCDM0    cold dark matter density today
    	@param OmegaBaryon0 baryon density today                              */
    int ReadCMBFast(string filename, double h100, double OmegaCDM0, double OmegaBaryon0);
  
    /** Read file output from CAMB 
    	@param filename     CAMB filename
    	@param h100         Hubble parameter in units of 100km/s/Mpc          */
    int ReadCAMB(string filename, double h100=0.71);

protected:
    double kmin_;      /**< not actually used?                                */
    double kmax_;      /**< not actually used?                                */
    int interptyp_;    /**< type of interpolation (see InterpTab in geneutils)*/
    vector<double> k_; /**< k values of transfer function read in             */
    vector<double> tf_;/**< values of transfer function read in               */
};


/** GrowthFactor class 
  *
  * Generic growth function class
  *
  */
class GrowthFactor : public ClassFunc1D {
public:

    /** Constructor */
    GrowthFactor(void) {};
    
    /** Desctructor */
    virtual ~GrowthFactor(void) {};
    
    /** @warning not implemented */
    virtual double DsDz(double z, double);
};


/** GrowthEH class 
  *
  * Growth function, eqn A4 in E&H 1998
  *
  */
class GrowthEH : public GrowthFactor {
public:

    /** Constructor 
        @param OmegaMatter0     matter density today
        @param OmegaLambda0     dark energy/cosmological constant density today*/
    GrowthEH(double OmegaMatter0, double OmegaLambda0);
    
    /** Copy constructor */
    GrowthEH(GrowthEH& d1) : O0_(d1.O0_) , Ol_(d1.Ol_) {};
  
    /** Destructor */
    virtual ~GrowthEH(void) {};
  
    /** Return growth function at redshift z 
        @param z    redshift of growth function                               */
    virtual double operator() (double z) const;
  
    /** ? differential ? */
    virtual double DsDz(double z, double dzinc=0.01);
  
    /** Set matter and dark energy/cosmological constant density today  
        @param OmegaMatter0     matter density today
        @param OmegaLambda0     dark energy/cosmological constant density today */
    void SetParTo(double OmegaMatter0,double OmegaLambda0);
  
    /** Set matter density today                                              */
    bool SetParTo(double OmegaMatter0);
  
    /** Return matter density today                                           */
    double ReturnOmegaM() { return O0_; };

protected:
  double O0_;   /**< matter density today                                     */
  double Ol_;   /**< dark energy/cosmological constant density today          */
};


/** GrowthFN class
  *
  * Solves linear ODE to calculate growth function
  *
  * Growth functions are linear, assume Newtonian gravity, matter domination
  * and Friedmann cosmology
  *
  * @todo include DE growth function [isn't this @todone?]
  * @todo debug and fix!
  *
  */
class GrowthFN : public GrowthFactor {
public:

    /** Constructor for DE growth 
    	@param OmegaM   		matter density today
    	@param OmegaL           dark energy density today
    	@param w0               dark energy equation of state parameter
    	@param wa               dark energy equation of state parameter
    	@param na               number of steps taken by ODE integrator
    	@param prt              print level                                   */
    GrowthFN(double OmegaM, double OmegaL, double w0, double wa=0, 
                                                     int na=1000000, int prt=0);

    /** Constructor for LCDM growth 
        @param OmegaM   		matter density today
    	@param OmegaL           dark energy density today                     */
    GrowthFN(double OmegaM, double OmegaL);
    
    /** Copy constructor */
    GrowthFN(GrowthFN& gro);
  
    /** Destructor */
    virtual ~GrowthFN(void) {};

    /** Return growth function at redshift z, normalized to 1 at z=0
        @param z    redshift of growth function                               */
    virtual double operator() (double z) const { 
            if(DEGrowth_)
			    return GrowthDE(z);
            else if (LCDMGrowth_)
                return GrowthLCDM(z);
            else 
                return -999; };

    /** Calculate growth from Carroll, Press & Turner approximation (can be 
        non-flat) but must be w=-1
        @parma z    redshift of growth function                               */
    double GrowthLCDM(double z) const;
    
    /** Calculate growth by solving the 2nd Order ODE equation for the perturbations 
        @parma z    redshift of growth function                               */
    double GrowthDE(double z) const;

    /** Set value of unnormalized growth function at redshift zero            */
    void GrowthZ0();
  
    /** A fudge to get the growth at redshift zero when dark energy is included */
    double GrowthDEZ0Fit(double omegam, double w0);
  
    /** Force calculation of growth by solving the 2nd Order ODE equation 
        for the perturbations even if w0=-1 and wa=0                          */
    void ForceDEGrowth(int na=1000000){ 
            if (DEGrowth_!= true)
					{ DEGrowth_=true; LCDMGrowth_=false; CalcDE(na); } 
			};
  
    /** Integrates ODE @warning not accurate! 
        @param na   number of steps taken by ODE integrator                   */
    void CalcDE(int na);

  // densities as a function of z
  //double DensTerm(double z){ return (Omz(z)
  
    /** Return \f$\Omega_m(1+z)^3\f$
        @param z    redshift                                                  */ 
    double Omz(double z){ return OmegaM_*(1+z)*(1+z)*(1+z); };
    
    /** Return \f$\Omega_k(1+z)^2\f$
        @param z    redshift                                                  */
    double OKz(double z){ return OmegaK_*(1+z)*(1+z); };
    
    /** Return \f$\Omega_{DE}a^{-3(w_0+w_a)}\exp{-3w_a(1-a)}\f$
        @param z    redshift                                                  */
    double ODEz(double z){ double a=1/(1+z); 
			return OmegaL_*pow( a,(-3*(1+w0_+wa_)) )*exp(-3*wa_*(1-a)); };
			
  //double OmegaMz(double z){ return Omz(z)/Ez(z); };
  //double OmegaDEz(double z){ return ODEz(z)/Ez(z); };
  //double OmegaKz(double z){ return OKz(z)/Ez(z); };
  
    /** Return \f$E(z)\f$ which is Hogg notation \f$E(z)\f$ squared
        @param z    redshift                                                  */
    double Ez(double z){ return Omz(z)+ODEz(z)+OKz(z); };
    
    /** Return dE/dlog(a) 
        @param z    redshift                                                  */
    double dEdlna(double z);
    
    /** Return dlog(E)/dlog(a) 
        @param z    redshift                                                  */
    double dlnEdlna(double z){ return dEdlna(z)/Ez(z); };

    /** Set matter density and cosmological constant/dark energy density 
        @param OmegaM   		matter density today
    	@param OmegaL           dark energy density today                     */
    void SetParTo(double OmegaM, double OmegaL)
	    { OmegaM_=OmegaM; OmegaL_=OmegaL; };
  
    /** Set matter density, cosmological constant/dark energy density, w0, wa
        @param OmegaM   		matter density today
    	@param OmegaL           dark energy density today                     
    	@param w0               dark energy equation of state parameter
    	@param wa               dark energy equation of state parameter       */
    void SetParTo(double OmegaM, double OmegaL,double w0,double wa)
		{ OmegaM_=OmegaM; OmegaL_=OmegaL; w0_=w0; wa_=wa; };
		
    /** Return value of unnormalized growth function at redshift zero         */
    double ReturnGrowthZ0(){ return D1z0_; };
    
    /** Return value of incorrect unnormalized growth function at redshift zero */
    double ReturnOtherGrowthZ0(){ return dz0_; };
    
    /** For debugging */
    void ReturnDeltaZVals(vector<double>& deltavals,vector<double>& zvals)
		{ deltavals=deltavals_; zvals=zvals_; };
    
    /** For debugging */
    void ReturnDeltaddash(vector<double>& dddash1,vector<double>& dddash2)
		{ dddash1=dddash1_; dddash2=dddash2_;};


protected:
    double OmegaM_;           /**< matter density                             */
    double OmegaL_;           /**< dark energy/cosmological constant density  */
    double OmegaK_;           /**< curvature density                          */
    double w0_;               /**< dark energy equation of state parameter    */
    double wa_;               /**< dark energy equation of state parameter    */
    double D1z0_;             /**< unnormalized growth at z=0                 */
	double dz0_;              /**< incorrect unnormalized growth at z=0       */
    long na_;                 /**< number of steps taken by ODE integrator    */
    bool DEGrowth_;           /**< if true calculate DE growth (integrate ODE)*/
    bool LCDMGrowth_;         /**< if true calculate LCDM growth (use CPT approx) */
    bool CalcDEDone_;         /**< if true DE growth calculation completed    */
    SInterp1D DEgrofunc_;     /**< DE growth function interpolation           */
    vector<double> deltavals_;/**< DE growth function interp table            */
    vector<double> zvals_;    /**< z values of DE growth function interp table*/
    vector<double> dddash1_;  /**< for debugging                              */
    vector<double> dddash2_;  /**< for debugging                              */
};


/** PkSpectrum class 
  *
  * Generic power spectrum class
  *
  */
class PkSpectrum : public ClassFunc1D, ClassFunc2D {
public:

    /** Define ReturnSpectrum type to describe what format power spectrum is
        returned in: Pk(k) or Delta^2(k) = k^3*Pk(k)/2Pi^2                    */
    typedef enum {PK=0, DELTA=1} ReturnSpectrum;

    /** Constructor */
    PkSpectrum(void) : zref_(0.) , scale_(1.) , typspec_(PK) {};
    
    /** Copy constructor */
    PkSpectrum(PkSpectrum& pk)
    : zref_(pk.zref_) , scale_(pk.scale_) , typspec_(pk.typspec_) {};
  
    /** Destructor */
    virtual ~PkSpectrum(void) {};
    
    /** Set redshift of power spectrum 
        @param z    redshift of power spectrum                                */
    virtual void SetZ(double z) {zref_ = z;}
  
    /** Return redshift of power spectrum                                     */
    virtual double GetZ(void) {return zref_;}
  
    /** Set power spectrum normalisation                                      */
    virtual void SetScale(double scale=1.) {scale_ = scale;};
    
    /** Return power spectrum normalisation                                   */
    virtual double GetScale(void) {return scale_;};
    
    /** Set return spectrum type Pk(k) or Delta^2(k) = k^3*Pk(k)/2Pi^2        
        @param typspec    if =PK Pk(k) else if =DELTA Delta^2(k)              */
    virtual void SetTypSpec(ReturnSpectrum typspec=PK) {typspec_ = typspec;};
    
    /** Return return spectrum type PK if Pk(k) or DELTA if 
        Delta^2(k) = k^3*Pk(k)/2Pi^2                                          */
    virtual ReturnSpectrum GetTypSpec(void) {return typspec_ ;}

protected:
  double zref_;               /**< redshift of power spectrum                 */
  double scale_;              /**< power spectrum normalization               */
  bool hunits_;               /**< units                                      */
  ReturnSpectrum typspec_;    /**< return type of power spectrum              */
};


/** PkSpecCalc class
  *
  * For calculating power spectrum from generic initial spectrum, transfer 
  * function and growth factor
  *
  * @note that this will be arbritarily normalised: use VarianceSpectrum to normalise 
  */
class PkSpecCalc : public PkSpectrum {
public:

    /** Constructor 
        @param pkinf    initial power spectrum
        @param tf       transfer function 
        @param d1       growth function
        @param zref     redshift of power spectrum                            */
    PkSpecCalc(InitialSpectrum& pkinf, TransferFunction& tf, GrowthFactor& d1, double zref=0.)
    : pkinf_(pkinf) , tf_(tf) , d1_(d1) { zref_ = zref; }; 
    
    /** Copy constructor                                                      */
    PkSpecCalc(PkSpecCalc& pkz)
    : pkinf_(pkz.pkinf_) , tf_(pkz.tf_) , d1_(pkz.d1_) {};
    
    /** Destructor                                                            */
    virtual ~PkSpecCalc(void) {};
    
    /** Return power spectrum at wavenumber k
        @param k    wavenumber in units ?                                     */
    virtual double operator() (double k) const { return (*this)(k, zref_); };
    
    /** Return power spectrum at wavenumber k and redshift z
        @param k    wavenumber in units ?                                                
        @param z    redshift                                                  */
    virtual double operator() (double k, double z) const;
    
    /** Return initial power spectrum at wavenumber k
        @param k    wavenumber in units ?                                     */
    double ReturnPki(double k) { return pkinf_(k); };
    
    /** Return transfer function at wavenumber k
        @param k    wavenumber in units ?                                     */
    double ReturnTrans(double k) { return tf_(k); };
    
    /** Return growth function at redshift z
        @param z    redshift                                                  */
    double ReturnGrowth(double z) { return d1_(z); };
    
    /** Return initial power spectrum                                         */
    InitialSpectrum& GetPkIni(void) {return pkinf_;}
    
    /** Return transfer function                                              */
    TransferFunction& GetTransfert(void) {return tf_;}
    
    /** Return growth function                                                */
    GrowthFactor& GetGrowthFactor(void) {return d1_;}

protected:
    InitialSpectrum& pkinf_;    /**< initial power spectrum                   */
    TransferFunction& tf_;      /**< transfer function                        */
    GrowthFactor& d1_;          /**< growth function                          */
};


/** PkTabulate class
  *
  * For reading power spectrum from a CAMB or CMBfast output file
  *
  * CAMB/CMBFast are openly available software packages that perform the full
  * numerical solution of the Bolzmann equations to calculate the transfer 
  * function
  *
  */
class PkTabulate : public PkSpectrum {
public:

    /** Constructor */
    PkTabulate(void) : kmin_(1.) , kmax_(-1.) , interptyp_(0), d1_(NULL)
        { k_.resize(0); pk_.resize(0); };
  
    /** Copy constructor */
    PkTabulate(PkTabulate& pkz) 
    : kmin_(pkz.kmin_) , kmax_(pkz.kmax_) , interptyp_(pkz.interptyp_) , 
      k_(pkz.k_) , pk_(pkz.pk_) , d1_(pkz.d1_) {};
  
    /** Destructor */
    virtual ~PkTabulate(void) {};
  
    /** Return power spectrum interpolated at wavenumber \f$k\f$ from values
        read from file 
        @param k    wavenumber in units?                                      */
    virtual double operator() (double k) const;
    
    /** Return power spectrum interpolated at wavenumber \f$k\f$ and redshift z
        from values read from file @warning not implemented
        @param k    wavenumber in units?
        @param z    redshift                                                  */
    virtual double operator() (double k, double z) const;
    
    /** Set new redshift of power spectrum using growth function 
        @param z    redshift                                                  */
    virtual void SetZ(double z);
    
    /** Return number of wavenumber values read from file                     */
    int NPoints(void) {return k_.size();}
    
    /** Set type of interpolation (see InterpTab in geneutils)                */           
    void SetInterpTyp(int typ=0);
  
    /** Read file output from CAMB 
    	@param filename     CAMB filename
    	@param h100tab      Hubble parameter in units of 100km/s/Mpc
        @param zreftab      redshift                                          */
    int ReadCAMB(string filename, double h100tab=0.71, double zreftab=0.);
  
    /** Set growth function in order to be able to make power spectrum a function
        of redshift too
        @param d1    growth function                                          */
    void SetGrowthFactor(GrowthFactor* d1) {d1_ = d1;}
  
    /** Return growth function                                                */
    GrowthFactor* GetGrowthFactor(void) {return d1_;}
  
    /** Return minimum wavenumber read from file                              */
    double KMin(void) {if(k_.size()!=0) return k_[0]; else return 0.;}
    
    /** Return maximum wavenumber read from file                              */
    double KMax(void) {if(k_.size()!=0) return k_[k_.size()-1]; else return 0.;}

protected:
    double kmin_;        /**< minimum wavenumber read from file                 */
    double kmax_;        /**< maximum wavenumber read from file                 */
    int interptyp_;      /**< type of interpolation (see InterpTab in geneutils)*/
    vector<double> k_;   /**< k values of power spectrum read in                */
    vector<double> pk_;  /**< values of power spectrum read in                  */
    GrowthFactor* d1_;   /**< growth function (used to change power spectrum z  */
};

/** PkEH class
  *
  * For calculating power spectrum from initial power law spectrum, and
  * E&H 1998 transfer function and growth factor
  *
  * @note that this will be arbritarily normalised: use VarianceSpectrum to normalise
  *
  */
class PkEH : public PkSpectrum {
public:

    /** Constructor 
        @param pkinf    initlai power spectrum
        @param tf       E&H 1998 transfer function
        @param d1       E&H 1998 growth function
        @param zref     redshift of power spectrum                            */
    PkEH(InitialPowerLaw& pkinf, TransferEH& tf, GrowthEH& d1, double zref=0.)
    : pkinf_(pkinf) , tf_(tf) , d1_(d1) { zref_ = zref; };
    
    /** Copy constructor */
    PkEH(PkEH& pkz) : pkinf_(pkz.pkinf_) , tf_(pkz.tf_) , d1_(pkz.d1_) {};
    
    /** Destructor */
    virtual ~PkEH(void) {};
    
    /** Return power spectrum at wavenumber k
        @param k    wavenumber in units?                                      */
    virtual double operator() (double k) const { return (*this)(k, zref_); };
    
    /** Return power spectrum at wavenumber k, redshift z
        @param k    wavenumber in units?
        @param z    redshift                                                  */
    virtual double operator() (double k, double z) const;
    
    /** Return initial power spectrum                                         */
    InitialPowerLaw& GetPkIni(void) {return pkinf_;};
    
    /** Return E&H 1998 transfer function                                     */
    TransferEH& GetTransfert(void) {return tf_;};
    
    /** Return E&H 1998 growth function                                       */
    GrowthEH& GetGrowthFactor(void) {return d1_;};

protected:
    InitialPowerLaw& pkinf_;    /**< initial power spectrum                   */
    TransferEH& tf_;            /**< E&H 1998 transfer function               */
    GrowthEH& d1_;              /**< E&H 1998 growth function                 */
};


/** PkAltNorm class
  *
  * For calculating power spectrum normalised with the amplitude of 
  * primordial curvature fluctuations at the pivot scale k =0.05 Mpc^-1
  *
  * Power spectrum then has no need to be normalised with VarianceSpectrum
  */
class PkAltNorm : public PkSpectrum {
public:

    /** Constructor 
        @param tf       E&H 1998 transfer function
        @param d1       growth function
        @param zref     redshift of power spectrum         
        @param ns       spectral index
        @param DelRsq   amplitude of primordial curvature fluctuations at the pivot scale */
	PkAltNorm(TransferEH& tf, GrowthFN& d1, double zref=0., double ns=1, double DelRsq=2.1e-9);
	
	/** Copy constructor */
	PkAltNorm(PkAltNorm& pkz);
	
	/** Destructor */
	virtual ~PkAltNorm(void) {};
	
	/** Return power spectrum at wavenumber k
        @param k    wavenumber in units?                                      */
	virtual double operator() (double k) const { return (*this)(k, zref_); };
	
	/** Return power spectrum at wavenumber k, redshift z
        @param k    wavenumber in units?
        @param z    redshift                                                  */
  	virtual double operator() (double k, double z) const;

    /** Return ns+3 */
	double Power() const { double pwr  = ns_ +3.; return pwr; };
	
	/** Return k/k_pivot 
	    @param k    wavenumber in units? 
	    @param h    Hubble parameter in units of 100km/s/Mpc                  */
	double kNorm(double k, double h) const {
        // k is in units of Mpc^-1 not h/Mpc
        // therefore don't need to divide 
        // kpivot by h (kpivot is in units of Mpc^-1 too)
        // but - am taking it out later
        double knorm=(h*k)/kpivot_;
        return knorm; };
	
	/** Return \f$ h^4 \f$ 
	    @param h    Hubble parameter in units of 100km/s/Mpc                  */
	double hpower4(double h) const { return h*h*h*h; };
	
	/** Return growth function at redshift z
	    @param z    redshift                                                  */
	double Gro(double z) const { return d1_(z); };
	
	/** Return growth^2/\f$\Omega_m^2$/f
	    @param OmegaM    matter density                                       */
	double PowFactor(double OmegaM) const {
        double d10 = d1_.ReturnGrowthZ0();
        double powfactor = (d10/OmegaM)*(d10/OmegaM); return powfactor;
        };
	
	/** Return transfer function at wavenumber k
        @param k    wavenumber in units?                                      */
	double TransF(double k){ return tf_(k); };
	
	/** Return amplitude at pivot scale                                       */
	double Apivot(){ return Apiv_; };
	
	/** Return amplitude of primordial curvature fluctuations                 */
	double DRsq(){ return DelRsq_; };

protected:
  TransferEH& tf_;        /**< E&H 1998 transfer function                     */
  GrowthFN& d1_;          /**< growth function                                */
  double ns_;             /**< spectral index                                 */
  double DelRsq_;         /**< amplitude of primordial curvature fluctuations */
  double kpivot_;         /**< pivot scale                                    */
  double Apiv_;           /**< amplitude at pivot scale                       */
};


/** VarianceSpectrum class
  *
  * For calculating power spectrum variance -> use SetScale in the power spectrum 
  * classes in the following way:
  *
  * var = Variance(kmin,kmax);
  * scale = (sigma8*sigma8)/var
  * SetScale(scale)
  *
  * sigma8 is the variance of matter fluctuations today (z=0) on 8Mpc/h scale 
  * The scale of 8Mpc/h is a convention choice
  * Observations show sigma8 is around 0.7 to 1 
  * If h100=1 power spectrum is normalised to 8Mpc, otherwise normalised to 8Mpc/h
  * (8Mpc/h is the conventional choice)
  */
class VarianceSpectrum : public ClassFunc1D {
public:

    /** Define TypeFilter type to describe window function to compute variance
        with                                                                  */
    typedef enum {TOPHAT=0, GAUSSIAN=1, NOFILTER=2} TypeFilter;

    /** Constructor 
        @param pk            power spectrum
        @param R             size of filter 
        @param typfilter     type of filter
        @param h100          Hubble parameter in units of 100km/s/Mpc         */
    VarianceSpectrum(ClassFunc1D& pk, double R, TypeFilter typfilter, double h100=1)
    : pk_(pk)
    { SetRadius(R); SetH100(h100); SetFilter(typfilter); };
  
    /** Copy constructor */
    VarianceSpectrum(VarianceSpectrum& vpk)
    : pk_(vpk.pk_) , R_(vpk.R_) , h100_(vpk.h100_) { SetFilter(vpk.typfilter_); };
    
    /** Destructor */
    virtual ~VarianceSpectrum(void) {};

    /** Set size of the filter
        @param R    size of the filter                                        */
    void SetRadius(double R);
    
    /** Set Hubble parameter
        @param h100    Hubble parameter in units of 100km/s/Mpc               */
    void SetH100(double h100);
    
    /** Set filter type, @notethe best approx of a top-hat filter w/R is a 
        Gaussian filter with (Rg=R/sqrt(5)) 
        @param typfilter    typfilter = TOPHAT : spherical 3D top-hat,
                                      = GAUSSIAN : spherical 3D gaussian
                                      = NOFILTER : no filter                  */
    void SetFilter(TypeFilter typfilter=TOPHAT) { typfilter_ = typfilter; };
  
    /** Set Gauss-Legendre integration parameters                             */
    void SetInteg(double dperc=0.1, double dlogkinc=-1., double dlogkmax=-1.,
                                                      unsigned short glorder=4);

    /** Return variance between wavenumber kmin and kmax                      */
    double Variance(double kmin,double kmax);

    // ATTENTION: The function to integrate is : f(k)dk = k^3*Pk(k)/(2Pi^2) *filter2(k*R) *dk/k
    // input k is in 1/Mpc units NOT h/Mpc
    virtual double operator() (double k) const {  
        double kunit=k/h100_;
		return (kunit*kunit)*(pk_(k)/h100_)*Filter2(kunit*R_)/(2.*M_PI*M_PI); };
  
    /** Return value of filter function squared                               */
    double Filter2(double x) const;

    // To aid the integration
    double FindMaximum(double kmin, double kmax, double eps=1.e-3);
    int FindLimits(double high, double &kmin, double &kmax, double eps=1.e-3);

protected:

    ClassFunc1D& pk_;        /**< power spectrum                              */
    TypeFilter typfilter_;   /**< type of filter to compute variance with     */
    double R_;               /**< size of filter                              */
    double h100_;            /**< Hubble parameter in units of 100km/s/Mpc    */
    double dperc_;           /**< Gauss-Legendre integration parameter        */
    double dlogkinc_;        /**< Gauss-Legendre integration parameter        */
    double dlogkmax_;        /**< Gauss-Legendre integration parameter        */
    unsigned short glorder_; /**< Gauss-Legendre integration parameter        */

};

} // end of namespace SOPHYA

#endif
