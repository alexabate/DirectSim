#ifndef MASSFUNC_SEEN
#define MASSFUNC_SEEN

#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>

#include "genericfunc.h"
#include "geneutils.h"
#include "pkspectrum.h"
#include "cosmocalcs.h"
#include "ctimer.h"

#include "constcosmo.h"

namespace SOPHYA {

/** PSparam class
  *  
  * Class that holds power spectrum parameters 
  * @todo is this class DEFUNCT? 
  */
class PSparam {
public:
	PSparam(double spec_index = 1,  double h_power = 1, double sigma_8=0.8)
	: spec_index_(spec_index) , h_power_(h_power) , sigma_8_(sigma_8) 
	{    };	
	
	void Print()
		{ cout <<"     Printing power spectrum parameters: ";
		  cout <<"     n = "<<spec_index_<<", h = "<<h_power_<<", sigma_8 = "<<sigma_8_<<endl;}

	void Setn(double spec_index){spec_index_ = spec_index;};
	void Seth(double h_power){h_power_=h_power;};
	void Setsig8(double sigma_8){sigma_8_=sigma_8;};

	//static inline double h_power() {return h_power_;};
	//static inline double n_power() {return spec_index_;};
	//static inline double s8_power() {return sigma_8_;};

protected:
	double spec_index_, h_power_, sigma_8_;
};


/** MassFunc class
  *  
  * Class for calculating the (dark matter) mass function of universe 
  * @todo Update to take different power spectrum class PkAltNorm
  */
class MassFunc : public ClassFunc1D {
public:

    /** Constructor
        @param su       class for cosmological calculations
        @param pkz      class that calculates matter power spectrum
        @param zref     redshift of mass function
        @param sig8     \f$\sigma_8$/f normalization of power spectrum
        @param TypeLog  if true return dn/dlogM instead of dn/dM
        @param IntType  if true do integration some way?
        @param MFType   Which mass function to return, 0 = Sheth-Torman       */
	MassFunc(SimpleUniverse& su, PkSpecCalc& pkz, double zref, double sig8, 
	                    bool TypeLog=false, bool IntType = 1, int MFType = 0);

	/** Destructor */        
  	virtual ~MassFunc(void){ };
  

    /** Return mass function value at mass m 
        @param m    mass in solar masses?                                     */
	virtual double operator() (double m) const {  
				if (TypeLog_)
					return (rhobar0()/m)*fST(m)*fabs(dlnsigdlnm(m,lmstep_))*log(10);
				else
					return (rhobar0()/(m*m))*fST(m)*fabs(dlnsigdlnm(m,lmstep_))*log(10);
				};

	/** Return variance in sphere of radius @param R
	    @param R    radius of sphere in Mpc?                                  */		
	double FindVar(double R) const;

	/** Sheth-Torman function 
	    @param m    mass in solar masses?                                     */
	double fST(double mv) const;

	/** Mass variance as a function of mass
	    @param mv   mass in solar masses?                                     */
	double sigsqM(double mv) const;

	/** Return \f$ d\log\sigma / d\logm $/f: log of mass variance differentiated
	    with respect to the log of the mass
	    @param mv       mass in soloar masses?
	    @param lmstep   step in log mass                                      */
	double dlnsigdlnm(double mv, double lmstep) const;

	/** Set the step in the log of the mass           
	    @param lmstep   step size in log of the mass                          */
	void Setlmstep(double lmstep) { lmstep_ = lmstep; return; };

	/** Set redshift of mass function 
	    @param z    redshift of mass function                                 */
	void SetZ(double z) {
		zref_ = z;
		pk_.SetZ(zref_);
		return;
		}

	/** Integrate mass function between two masses
	    @param Mmin     minimum mass in solar masses?
	    @param Mmax     maximum mass in solar masses?
	    @param npt      number of points to use for integration               */
	double Integrate(double Mmin,double Mmax,int npt=100) {
		 if(npt<1) npt = 100;
		 double perc=0.01, dlxinc=(Mmax-Mmin)/npt, dlxmax=10.*dlxinc; 
		 unsigned short glorder=4;
		 double sum = IntegrateFunc(*this,Mmin,Mmax,perc,dlxinc,dlxmax,glorder);
		 return sum;
		 }

	/** Return mean matter density today (z=0)                                */
	inline double rhobar0() const { return rho_crit()*Om_; };

	/** Write out the power spectrum (for debugging)
	    @param outfile  file name to output power spectrum to                 */
	void WritePS2File(string outfile);

	// --------- Useful Constants 
	
	/** Return solar mass in KG                                               */
	static inline double SolarMassKG() { return SOLARMASS_IN_KG; }
	/** Return Gravitational constant in SI units                             */
	static inline double Gconst() { return G_NEWTON_SI; }
	/** Return number of meters in a Mpc                                      */
	static inline double MpctoM() { return NMETERS_IN_MPC; }
	/** Return Solar mass * Gravitational Constant in units s^-2 Msolar^-1 Mpc^3 */
	static inline double Gcosmo() { return (Gconst()*SolarMassKG())/pow(MpctoM(),3); }
	/** Return \f$\delta_c$\f: overdensity for collapse (spherical model)     */
	static inline double deltac() { return DELTA_C; }
	/** Return H0 in units of h/s                                             */
	static inline double H0hs() { return 100.*1.e3/MpctoM(); }
	/** Return critical density in units of h^2 M_solar / Mpc^3                */
	static inline double rho_crit() { return (3*H0hs()*H0hs()) / (8*PI*Gcosmo()); }
	
	/** Return Sheth-Torman mass function parameter \f$a$\f                   */
	static inline double aST() { return a_ST_Cst; }
	/** Return Sheth-Torman mass function parameter \f$A$\f                   */
	static inline double AST() { return A_ST_Cst; }
	/** Return Sheth-Torman mass function parameter \f$Q$\f                   */
	static inline double QST() { return Q_ST_Cst; }
	/** Return CMB temperature                                                */
	static inline double TCMB() { return T_CMB_K; }
	
protected:
	SimpleUniverse& su_;	/**< Holds cosmological parameters/calculations   */
	PkSpecCalc& pk_;	    /**< Holds power spectrum                         */
	double zref_;		    /**< redshift of mass function                    */
	bool TypeLog_; 		    /**< if true return dn/dlogm rather than dn/dm    */
	bool IntType_;          /**< if true does more accurate integration?      */
	int MFType_; 		    /**< Which mass function to return, 0 = Sheth-Torman */
	double lmstep_;		    /**< step with which to compute derivation of dlnsig/dlnm */
	double Om_;             /**< \f$\Omega_m$\f */
	

	// constants
	static double a_ST_Cst;	/**< Sheth-Torman mass function parameter         */
	static double A_ST_Cst;	/**< Sheth-Torman mass function parameter         */
	static double Q_ST_Cst;	/**< Sheth-Torman mass function parameter         */

};


} // end of namespace SOPHYA
#endif
