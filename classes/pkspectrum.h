/**
 * @file  pkspectrum.h
 * @brief Calculates linear transfer functions and growth functions
 *
 * @todo Add non linear power spectrum calculation
 *
 * @author Christophe Magneville
 * Contact:
 *
 * Created on: 2008
 * @date 2008
 *
 */
#ifndef PKSPECTRUM_SEEN
#define PKSPECTRUM_SEEN

#include "machdefs.h"
#include "genericfunc.h"
#include "sinterp.h"

// Written by C. Magneville
// Modified by AA

// Note that transfer functions, power spectra are outputed 
// with k in units of Mpc^-1 not hMpc^-1

namespace SOPHYA {

// Generic class for inital spectrum
class InitialSpectrum : public GenericFunc {
public:
  InitialSpectrum(void) {};
  InitialSpectrum(InitialSpectrum& pkinf) {};
  virtual ~InitialSpectrum(void) {};
};

// Class to calculate the initial early universe power-law spectrum
class InitialPowerLaw : public InitialSpectrum {
public:
  InitialPowerLaw(double n,double a=1.);
  InitialPowerLaw(InitialPowerLaw& pkinf);
  virtual ~InitialPowerLaw(void);
  virtual double operator() (double k) { return A_ * pow(k,n_);}
  void SetNorm(double a) {A_ = a;}
  void SetSlope(double n) {n_ = n;}
protected:
  double n_, A_;
};

// Generic transfer function class
class TransferFunction : public GenericFunc {
public:
  TransferFunction(void) {};
  virtual ~TransferFunction(void) {};
};

// Class to calculate Eisenstein & Hu transfer function
// k is in units of 1/Mpc (NOT h/Mpc)
class TransferEH : public TransferFunction {
public:

  typedef enum{ALL=0, CDM=1, BARYON=2} ReturnPart;

  TransferEH(double h100,double OmegaCDM0,double OmegaBaryon0,double tcmb,bool nobaryon=false,int lp=0);
  TransferEH(TransferEH& tf);
  virtual ~TransferEH(void);
  bool SetParTo(double h100,double OmegaCDM0,double OmegaBaryon0);
  virtual double operator() (double k);
  double KPeak(void);
  void SetNoOscEnv(unsigned short nooscenv=0);
  void SetReturnPart(ReturnPart retpart=ALL);
  void SetPrintLevel(int lp=0);
  double Returnh100() { return h100_; };
  double ReturnOmegaM() { return O0_; };
protected:
  int lp_;
  double O0_,Oc_,Ob_,h100_,tcmb_;
  double th2p7_;
  double zeq_,keq_,zd_,Req_,Rd_,s_,ksilk_,alphac_,betac_,bnode_,alphab_,betab_;
  double alphag_;
  double sfit_,kpeak_;

  bool nobaryon_;
  unsigned short nooscenv_;
  ReturnPart retpart_;

  double T0tild(double k,double alphac,double betac);
  void Init_(void);
  void zero_(void);
};

// For reading transfer function from a CAMB or CMBfast output file
class TransferTabulate : public TransferFunction {
public:
  TransferTabulate(void);
  TransferTabulate(TransferTabulate& tf);
  virtual ~TransferTabulate(void);
  virtual double operator() (double k);
  int NPoints(void) {return k_.size();}
  void SetInterpTyp(int typ=0);
  int ReadCMBFast(string filename,double h100,double OmegaCDM0,double OmegaBaryon0);
  int ReadCAMB(string filename, double h100=0.71);
protected:
  double kmin_,kmax_;
  int interptyp_;
  vector<double> k_, tf_;
};


// Generic growth factor class
class GrowthFactor : public GenericFunc {
public:
  GrowthFactor(void) {};
  virtual ~GrowthFactor(void) {};
  virtual double DsDz(double z, double);
};

// Growth function, eqn A4 in E&H 1998
class GrowthEH : public GrowthFactor {
public:
  GrowthEH(double OmegaMatter0,double OmegaLambda0);
  GrowthEH(GrowthEH& d1);
  virtual ~GrowthEH(void);
  virtual double operator() (double z);
  virtual double DsDz(double z,double dzinc=0.01);
  void SetParTo(double OmegaMatter0,double OmegaLambda0);
  bool SetParTo(double OmegaMatter0);
  double ReturnOmegaM() { return O0_; };
protected:
  double O0_,Ol_;
};

// New class to calculate the growth function
// Will include DE growth function at some point
// Growth functions are linear, assume Newtonian gravity, matter domination
// and Friedmann cosmology
class GrowthFN : public GrowthFactor {
public:
// constructor for DE growth
  GrowthFN(double OmegaM,double OmegaL,double w0,double wa=0,int na=1000000,int prt=0);
// constructor for LCDM growth
  GrowthFN(double OmegaM,double OmegaL);
  GrowthFN(GrowthFN& gro);
  virtual ~GrowthFN(void);

  virtual double operator() (double z) { if(DEGrowth_)
						return GrowthDE(z);
					 else if (LCDMGrowth_)
						return GrowthLCDM(z);
					 else 
						return -999; };

  double GrowthLCDM(double z);// Carroll, Press & Turner approximation (can be non-flat)
  // this growth function solve the 2nd Order ODE equation for 
  // the perturbations 
  //double GrowthDEwconst(double z);
  double GrowthDE(double z);// 


  void GrowthZ0();
  double GrowthDEZ0Fit(double omegam,double w0);
  void ForceDEGrowth(int na=1000000){ if (DEGrowth_!= true)
					{DEGrowth_=true; LCDMGrowth_=false; CalcDE(na);} };
  void CalcDE(int na);

  // densities as a function of z
  //double DensTerm(double z){ return (Omz(z)
  double Omz(double z){ return OmegaM_*(1+z)*(1+z)*(1+z); };
  double OKz(double z){ return OmegaK_*(1+z)*(1+z); };
  double ODEz(double z){ double a=1/(1+z); 
			return OmegaL_*pow( a,(-3*(1+w0_+wa_)) )*exp(-3*wa_*(1-a)); };
  //double OmegaMz(double z){ return Omz(z)/Ez(z); };
  //double OmegaDEz(double z){ return ODEz(z)/Ez(z); };
  //double OmegaKz(double z){ return OKz(z)/Ez(z); };
  double Ez(double z){ return Omz(z)+ODEz(z)+OKz(z); };
  double dEdlna(double z);
  double dlnEdlna(double z){ return dEdlna(z)/Ez(z); };

  void SetParTo(double OmegaM, double OmegaL)
		{ OmegaM_=OmegaM; OmegaL_=OmegaL; };
  void SetParTo(double OmegaM, double OmegaL,double w0,double wa)
		{ OmegaM_=OmegaM; OmegaL_=OmegaL; w0_=w0; wa_=wa; };
  double ReturnGrowthZ0(){ return D1z0_; };
  double ReturnOtherGrowthZ0(){ return dz0_; };
  void ReturnDeltaZVals(vector<double>& deltavals,vector<double>& zvals)
		{ deltavals=deltavals_; zvals=zvals_; };
  void ReturnDeltaddash(vector<double>& dddash1,vector<double>& dddash2)
		{ dddash1=dddash1_; dddash2=dddash2_;};


protected:
  double OmegaM_, OmegaL_, OmegaK_, w0_, wa_;
  double D1z0_; // growth at z=0
	double dz0_;
  long na_; //  
  bool DEGrowth_,LCDMGrowth_,CalcDEDone_;
  SInterp1D DEgrofunc_;
  vector<double> deltavals_,zvals_;
  // for debugging
  vector<double> dddash1_,dddash2_;
};


// Generic power spectrum class
class PkSpectrum : public GenericFunc {
public:
  // typsec = PK : compute Pk(k)
  //        = DELTA : compute Delta^2(k) = k^3*Pk(k)/2Pi^2
  typedef enum {PK=0, DELTA=1} ReturnSpectrum;

  PkSpectrum(void);
  PkSpectrum(PkSpectrum& pk);
  virtual ~PkSpectrum(void) {};
  virtual void   SetZ(double z) {zref_ = z;}
  virtual double GetZ(void) {return zref_;}
  // use the function SetScale below to set power spectrum normalisation
  virtual void SetScale(double scale=1.) {scale_ = scale;}
  virtual double GetScale(void) {return scale_;}
  virtual void SetTypSpec(ReturnSpectrum typspec=PK) {typspec_ = typspec;}
  virtual ReturnSpectrum GetTypSpec(void) {return typspec_ ;}
protected:
  double zref_, scale_;
  bool hunits_;
  ReturnSpectrum typspec_;
};

// For calculating power spectrum from generic initial spectrum, transfer function and growth factor
// Note that this will be arbritarily normalised: use VarianceSpectrum to normalise
class PkSpecCalc : public PkSpectrum {
public:
  PkSpecCalc(InitialSpectrum& pkinf,TransferFunction& tf,GrowthFactor& d1,double zref=0.);
  PkSpecCalc(PkSpecCalc& pkz);
  virtual ~PkSpecCalc(void);
  virtual double operator() (double k);
  virtual double operator() (double k,double z);
  double ReturnPki(double k) { return pkinf_(k); };
  double ReturnTrans(double k) { return tf_(k); };
  double ReturnGrowth(double z) { return d1_(z); };
  InitialSpectrum& GetPkIni(void) {return pkinf_;}
  TransferFunction& GetTransfert(void) {return tf_;}
  GrowthFactor& GetGrowthFactor(void) {return d1_;}
protected:
  InitialSpectrum& pkinf_;
  TransferFunction& tf_;
  GrowthFactor& d1_;
};

// For reading power spectrum from a CAMB or CMBfast output file
class PkTabulate : public PkSpectrum {
public:
  PkTabulate(void);
  PkTabulate(PkTabulate& pkz);
  virtual ~PkTabulate(void);
  virtual double operator() (double k);
  virtual double operator() (double k,double z);
  virtual void SetZ(double z);
  int NPoints(void) {return k_.size();}
  void SetInterpTyp(int typ=0);
  int ReadCAMB(string filename, double h100tab=0.71, double zreftab=0.);
  void SetGrowthFactor(GrowthFactor* d1) {d1_ = d1;}
  GrowthFactor* GetGrowthFactor(void) {return d1_;}
  double KMin(void) {if(k_.size()!=0) return k_[0]; else return 0.;}
  double KMax(void) {if(k_.size()!=0) return k_[k_.size()-1]; else return 0.;}
protected:
  double kmin_,kmax_;
  int interptyp_;
  vector<double> k_, pk_;
  GrowthFactor* d1_;
};

// For calculating power spectrum from initial power law spectrum, and
// E&H 1998 transfer function and growth factor
// Note that this will be arbritarily normalised: use VarianceSpectrum to normalise
class PkEH : public PkSpectrum {
public:
  PkEH(InitialPowerLaw& pkinf,TransferEH& tf,GrowthEH& d1,double zref=0.);
  PkEH(PkEH& pkz);
  virtual ~PkEH(void);
  virtual double operator() (double k);
  virtual double operator() (double k,double z);
  InitialPowerLaw& GetPkIni(void) {return pkinf_;}
  TransferEH& GetTransfert(void) {return tf_;}
  GrowthEH& GetGrowthFactor(void) {return d1_;}
protected:
  InitialPowerLaw& pkinf_;
  TransferEH& tf_;
  GrowthEH& d1_;
};

// For calculating power spectrum normalised with the amplitude of 
// primordial curvature fluctuations at the pivot scale k =0.05 Mpc^-1
// Power spectrum then has no need to be normalised with VarianceSpectrum
class PkAltNorm : public PkSpectrum {
public:
	PkAltNorm(TransferEH& tf,GrowthFN& d1,double zref=0.,double ns=1,double DelRsq=2.1e-9);
	PkAltNorm(PkAltNorm& pkz);
	virtual ~PkAltNorm(void);
	virtual double operator() (double k);
  	virtual double operator() (double k,double z);

	double Power();
	double kNorm(double k, double h);
	double hpower4(double h);
	double Gro(double z);
	double PowFactor(double OmegaM);
	double TransF(double k);
	double Apivot(){ return Apiv_; };
	double DRsq(){ return DelRsq_; };

protected:
  TransferEH& tf_;
  GrowthFN& d1_;
  double ns_,DelRsq_;
  double kpivot_;
  double Apiv_;
};

// For calculating power spectrum variance -> use SetScale in the power spectrum 
// classes in the following way:
// var = Variance(kmin,kmax);
// scale = (sigma8*sigma8)/var
// SetScale(scale)
// sigma8 is the variance of matter fluctuations today (z=0) on 8Mpc/h scale 
// The scale of 8Mpc/h is a convention choice
// Observations show sigma8 is around 0.7 to 1 
// If h100=1 power spectrum is normalised to 8Mpc, otherwise normalised to 8Mpc/h
// (8Mpc/h is the conventional choice)
class VarianceSpectrum : public GenericFunc {
public:

  typedef enum {TOPHAT=0, GAUSSIAN=1, NOFILTER=2} TypeFilter;

  VarianceSpectrum(GenericFunc& pk,double R,TypeFilter typfilter,double h100=1);
  VarianceSpectrum(VarianceSpectrum& vpk);
  virtual ~VarianceSpectrum(void);

  void SetRadius(double R);
  void SetH100(double h100);
  void SetFilter(TypeFilter typfilter=TOPHAT);
  void SetInteg(double dperc=0.1,double dlogkinc=-1.,double dlogkmax=-1.,unsigned short glorder=4);

  double Variance(double kmin,double kmax);

  // ATTENTION: The function to integrate is : f(k)dk = k^3*Pk(k)/(2Pi^2) *filter2(k*R) *dk/k
  // input k is in 1/Mpc units NOT h/Mpc
  virtual double operator() (double k) {  double kunit=k/h100_;
					  return (kunit*kunit)*(pk_(k)/h100_)*Filter2(kunit*R_)/(2.*M_PI*M_PI); }
  double Filter2(double x);

  // Help for the integration
  double FindMaximum(double kmin,double kmax,double eps=1.e-3);
  int FindLimits(double high,double &kmin,double &kmax,double eps=1.e-3);

protected:

  GenericFunc& pk_;
  TypeFilter typfilter_;
  double R_,h100_;

  double dperc_,dlogkinc_,dlogkmax_;
  unsigned short glorder_;

};

} // end of namespace SOPHYA

#endif
