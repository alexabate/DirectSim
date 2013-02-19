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
#include "luc.h"
#include "ctimer.h"

#define PI 3.141592

namespace SOPHYA {

// Class that holds power spectrum parameters DEFUNCT?
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

// Class for calculating the mass function

class MassFunc : public GenericFunc {
public:
	MassFunc(SimpleUniverse&, PkSpectrumZ&, double zref, double sig8, bool TypeLog=false, bool IntType = 1, int MFType = 0);
  	virtual ~MassFunc(void);
  

	virtual double operator() (double m) {  
						if (TypeLog_)
							return (rhobar0()/m)*fST(m)*abs(dlnsigdlnm(m,lmstep_))*log(10);
						else
							return (rhobar0()/(m*m))*fST(m)*abs(dlnsigdlnm(m,lmstep_))*log(10);    }
	double FindVar(double R);
	double fST(double mv);
	double sigsqM(double mv);
	double dlnsigdlnm(double mv, double lmstep);
	void Setlmstep(double lmstep) { lmstep_ = lmstep; return; }

	// to change zref of mass function
	void SetZ(double z)
		{	
		zref_ = z;
		pk_.SetZ(zref_);
		return;
		}

	// integration function
	double Integrate(double Mmin,double Mmax,int npt=100) 
		{
		 if(npt<1) npt = 100;
		 double perc=0.01, dlxinc=(Mmax-Mmin)/npt, dlxmax=10.*dlxinc; 
		 unsigned short glorder=4;
		 double sum = IntegrateFunc(*this,Mmin,Mmax,perc,dlxinc,dlxmax,glorder);
		 return sum;
		 }

	// Mean matter density today (z=0)
	inline double rhobar0() { return rho_crit()*Om_; }

	// Output the power spectrum
	void WritePS2File(string outfile);

	// --------- Useful Constants 
	

	// Solar mass in KG
	static inline double SolarMassKG() { return MSolar_Cst; }
	// Gravitational constant
	static inline double Gconst() { return G_Newton_Cst; }
	// Mpc in meters
	static inline double MpctoM() { return MpctoMeters_Cst; }
	// Solar mass * Gravitational Constant in units s^-2 Msolar^-1 Mpc^3
	static inline double Gcosmo() { return (Gconst()*SolarMassKG())/pow(MpctoM(),3); }
	// Delta_c: overdensity for collapse (spherical model)
	static inline double deltac() { return delta_c_Cst; }
	// H0 in units of h/s
	static inline double H0hs() { return 100.*1.e3/MpctoM(); }
	// Critical density in units of h^2 M_solar / Mpc^3
	static inline double rho_crit() { return (3*H0hs()*H0hs()) / (8*PI*Gcosmo()); }
	
	// Sheth-Torman mass function parameters
	static inline double aST() { return a_ST_Cst; }
	static inline double AST() { return A_ST_Cst; }
	static inline double QST() { return Q_ST_Cst; }
	// CMB temperature
	static inline double TCMB() { return TCMB_Cst; }
protected:
	SimpleUniverse& su_;	    // Holds cosmological parameters
	PkSpectrumZ& pk_;	    // Holds power spectrum 
	double zref_;		    // redshift of mass function
	bool TypeLog_; 		    // if true return dn/dlogm rather than dn/dm
	int MFType_; 		    // Which mass function to return, 0 = Sheth-Torman
	double lmstep_;		    // step with which to compute derivation of dlnsig/dlnm
	double Om_;
	bool IntType_;

	// constants
	static double MSolar_Cst;	// Solar mass in KG
	static double delta_c_Cst;	// Overdensity for collapse (spherical model)
	static double a_ST_Cst;		// Sheth-Torman mass function parameter
	static double A_ST_Cst;		// Sheth-Torman mass function parameter
	static double Q_ST_Cst;      	// Sheth-Torman mass function parameter
	static double G_Newton_Cst; 
	static double MpctoMeters_Cst;
	static double TCMB_Cst;

};



} // end of namespace SOPHYA

#endif
