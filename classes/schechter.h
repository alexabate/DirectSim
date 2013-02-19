/**
 * @file  schechter.h
 * @brief Calculations with Schechter functions and luminosity function parameters
 *
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2010
 * @date 2010
 *
 */

#ifndef SCHECHTER_SEEN
#define SCHECHTER_SEEN

#include "machdefs.h"
#include "genericfunc.h"
#include "geneutils.h"
#include "cosmocalcs.h"
#include "ctimer.h"

namespace SOPHYA {

class Histo;
class FunRan;

/******* LFParameters *********************************************************/

/** @class LFParameters class
  *
  * Class to read in and hold luminosity function parameters
  *
  */
class LFParameters {
public:

	/** Default constructor */
	LFParameters(){ };
	
	/** Constructor, reads in pars from LFfile and puts into table & vectors */
	LFParameters(string LFfile, int prtvl=0);
	
	/** Constructor: calculates parameter interpolation: does not have to
	    read in from file */
	LFParameters(double zmax,double dz,string type="dahlen");
	
	/** Copy constructor (if id=0 copy Dahlen-style, if id=1 copy file-read-style) 
	    @warning Not currently functional */
	LFParameters(LFParameters const& a, int id=0);
	
	virtual ~LFParameters(void){ };


    /** Returns the Schechter parameters (of ALL galaxies) at the redshift z 
        @note should this be a proper interpolation rather than be defined on a
        grid ? 
        @note also, parameters read from a file are NOT returned this way     */
	 virtual double operator()(double z,double& ps, double& ms, double& a) {  	
	 
	        double minv=1000; int idmin=-1; double diff;
		    int nz = zgrid_.size();
			for (int i=0; i<nz; i++) {
				diff=z-zgrid_[i];
				
				if (abs(diff)<minv)
					{ minv=diff; idmin=i; }
				}
			if (idmin<0)
				throw ParmError("Redshift not found within grid");
			//cout <<"minv="<<minv<<", diff="<<diff<<", idmin="<<idmin<<endl;
			ps=pstargrid_[idmin];
			ms=mstargrid_[idmin];
			a=alpgrid_[idmin];
			};


    // functions that do the file reading
    
	/** Reads from a file and stores the Schechter parameters in vectors      
	    @param LFfile   file to read from
	    @param prtvl    if prtval>0 prints parameters to the screen           */
	void File(string LFfile,int prtvl=0);
	
	/** Reads Schechter parameters from a file located at environment variable
	    value $LFLOC and stores them in an array                              
	    @param LFfile   file to read from                                     */
	TArray<r_8> ReadFile(string LFfile);
	
	/** Stores the Schecter parameters in separate vectors                    
	    @param LFTable  Array of Schechter parameters returned by ReadFile    */
	void StorePars(TArray<r_4> LFTable);
	
	/** Does the "Dahlen" interpolation of the Schecter parameters, can add more*/
	void Dahlen(double zmax,double dz);
	
	
	// functions that print stuff to the screen for checking
	
	void PrintTable(TArray<r_4> LFTable);
	
	void PrintZbins();
	
	// prints the LF pars in each z bin for galaxy of type "type"
	void PrintPars(int type=0);
	
	/** Return the parameters read from the file for the "All galaxies"       */
	void ReturnAllLF(vector<double>& z,vector<double>& p,vector<double>& a,vector<double>& m)
		{ z=zcs_; m=AMstars_; a=Aalphas_; p=Aphistars_; };
		
	/** Return the min/max redshift bounds of all the LF redshift bins read from
	    the file */
	void ReturnZbins(vector<double>& zmins, vector<double>& zmaxs, int& nb)
		{ zmins=zmins_; zmaxs=zmaxs_; nb=nZbins_; };
		
	/** Return the min and max redshifts where the LFs read from the file are 
	    defined */
	void ReturnMinMaxZbin(double& zmin, double& zmax)
		{ zmin=zmins_[0]; zmax=zmaxs_[zmaxs_.size()]; };
		
	// check over which redshift bins the galaxy survey spans
	void CheckBinRange(int& i_min, int& i_max, int& nb,double zmin, double zmax);
	
	// return the LF pars for redshift bin "ib" and for galaxy type "type"
	// type=0, All LF; type=1, Early LF; type=2, Late LF; type=3, SB LF
	void ReturnParsBini(double& Mstar,double& alpha,double& phistar,int ib,int type=0);
	
	/** Return the redshift bounds for bin "ib" read from the file            */
	void ReturnZbin(double& z1, double& z2, int ib)
		{ z1=zmins_[ib]; z2=zmaxs_[ib]; };

	vector<double>  ReturnAll()
	    { return all_;};
	
	vector<double> ReturnEarly()
		{ return early_;};

	vector<double>  ReturnLate()
		{ return late_;};

	vector<double>  ReturnSB()
		{ return sb_;};

protected:
	// Variables that help with reading in parameters from file, then storing them */
	string type_;           /**< type of LF parameter interpolation/extrapolation  */
	int ZminCol_;           /**< column in file containing redshift lower edge     */
	int ZmaxCol_;           /**< column in file containing redshift upper edge     */
	int ACol_;              /**< 1st column in file containing all galaxies pars   */
	int ECol_;              /**< 1st column in file containing early galaxies pars */
	int SCol_;              /**< 1st column in file containing spiral galaxies pars*/
	int SBCol_;             /**< 1st column in file containing starburst galaxies pars*/
	int MstarCol_;          /**< column order of M* parameters                */
	int alphaCol_;          /**< column order of alpha parameters             */
	int phistarCol_;        /**< column order of phi* parameters              */
	// These are read in from a file of LF Schechter parameters:              */
	vector<double> zcs_;        /**< redshift Schechter parameters defined at */
	vector<double> zmins_;      /**< redshift lower edge Schechter parameters defined at */
	vector<double> zmaxs_;      /**< redshift upper edge Schechter parameters defined at */
	vector<double> AMstars_;    /**< M* parameter for all galaxies              */
	vector<double> Aalphas_;    /**< alpha parameter for all galaxies           */
	vector<double> Aphistars_;  /**< phi* parameter for all galaxies            */
	vector<double> EMstars_;    /**< M* parameter for early type galaxies       */
	vector<double> Ealphas_;    /**< alpha parameter for early type galaxies    */
	vector<double> Ephistars_;  /**< phi* parameter for early type galaxies     */
	vector<double> SMstars_;    /**< M* parameter for spiral type galaxies      */
	vector<double> Salphas_;    /**< alpha parameter for spiral type galaxies   */
	vector<double> Sphistars_;  /**< phi* parameter for spiral type galaxies    */
	vector<double> SBMstars_;   /**< M* parameter for starburst type galaxies   */
	vector<double> SBalphas_;   /**< alpha parameter for starburst type galaxies*/
	vector<double> SBphistars_; /**< phi* parameter for starburst type galaxies */
	int nZbins_;                /**< number of redshift bins in file            */
	string MstarUnits_;         /**< units of M* parameters                     */
	string phistarUnits_;       /**< units of phi* parameters                   */
	// gridded Schechter parameter values                                       */
	vector<double> zgrid_;      /**< the interpolated z values                  */
	vector<double> mstargrid_;  /**< the interpolated M* values                 */
	vector<double> pstargrid_;  /**< the interpolated phi* values               */
	vector<double> alpgrid_;    /**< the interpolated alpha values              */
	vector<double> all_;
	vector<double> early_;
	vector<double> late_;
	vector<double> sb_;
};

/*********************** Schechter Function classes **************************/

/** @class Schechter class
  *
  * Returns the Schechter function \f$\phi(M)\f$ computed for a set of parameters: 
  * \f$\phi_*, M_*, \alpha \f$
  *
  * The Schechter function is: \f$ \phi(M)=0.4\log(10)\phi_*x^{\alpha+1}\exp(-x)
  * where x = 10^{-0.4(M-M_*)}
  *
  * The Schechter function describes the number density of galaxies per 
  * unit absolute magnitude (luminosity).  It can evolve with redshift and
  * galaxy type
  *
  */
class Schechter : public GenericFunc 
{
public:
    /** Constructor 
        @param phistar  Schechter parameter \f$ \phi_* \f$
        @param Mstar    Schechter parameter \f$ M_* \f$
        @param alpha    Schechter parameter \f$ \alpha \f$                    */
    Schechter(double phistar,double Mstar,double alpha);
    
    /** copy constructor */
    Schechter(Schechter& f);
                   
    Schechter(void);
    virtual ~Schechter(void);

    /** Set the output value, determine whether operator() returns dn/dm or 
        m*dn/dm */
    void SetOutValue(unsigned short outvalue=0);
      
    /** Return the output value, whether operator() returns dn/dm or 
        m*dn/dm */
    unsigned short GetOutValue(void);
      
    /** Reset the Schechter parameters to new values */
    void SetParam(double phistar,double Mstar,double alpha);
      
    /** Return current Schechter parameters */
    void GetParam(double& phistar,double& Mstar,double& alpha);
      
    /** Set the integration parameters (if want different from default) */
    void SetInteg(double Mmin,double Mmax, int npt)
        { Mmin_=Mmin; Mmax_=Mmax; npt_=npt; };

    /** Return the LF value for absolute magnitude "absmag" */
    virtual double operator() (double absmag);

    // integration functions
    double Integrate(); // an integration function
    double IntegrateQuick();// don't think this works?

    // number of gals (sky area is in steradians)
    double NGals(SimpleUniverse su,double skyarea,double z1,double z2);
    double Volume(SimpleUniverse su,double skyarea,double z1,double z2);

    // check functions
    virtual void Print(string, string);
    virtual void Print(void);
    void CheckFaintLimit(Schechter&);

protected:
    double phistar_,Mstar_,alpha_; // Schechter function parameters
    unsigned short outvalue_;
    double Mmin_,Mmax_;// for the integration
    int npt_;
    TVector<r_4> ngal_by_mpc3_, schmaxv_; //debugging/checking
};

/** @class SchechterVol class
  *
  * Returns the Schechter function \f$\phi(M|z)\f$ computed for a set of parameters: 
  * \f$\phi_*(z), M_*(z), \alpha(z) \f$ multiplied by the volume element at
  * redshift \f$z\f$.
  *
  * Returns \f$\phi(M|z)dV(z)\f$
  */
class SchechterVol : public GenericFunc
{
public:
	// constructor
	SchechterVol(LFParameters& lfpars, SimpleUniverse& su)
		: lfpars_(lfpars) , su_(su) {   
			Mmin_=-24, Mmax_=-13;// units of "M-5log10h70"
  			npt_=10000;  };

	virtual double operator() (double absmag, double z)
		{ 
			double ps,ms,a;
			lfpars_(z,ps,ms,a);
			Schechter sch(ps,ms,a);
			su_.SetEmissionRedShift(z);
			double volel=su_.VolEl();
			return sch(absmag)*volel;
		}

	// To integrate: int phi(M|z)*dV(z) dM
	double Integrate(double z);

  	// set the integration parameters (if want different from default)
  	void SetInteg(double Mmin,double Mmax, int npt=10000)
		{ Mmin_=Mmin; Mmax_=Mmax; npt_=npt; };

protected:
	LFParameters& lfpars_;
	SimpleUniverse& su_;
	double Mmin_,Mmax_;
	int npt_;
};

// returns phi(M|z)
class SchechterM : public GenericFunc
{
public:
	// constructor
	SchechterM(LFParameters& lfpars)
		: lfpars_(lfpars) {   
			Mmin_=-24, Mmax_=-13;// units of "M-5log10h70"
  			npt_=10000;  };

	virtual double operator() (double absmag, double z) { 
	
				  	double ps,ms,a;
					lfpars_(z,ps,ms,a);
					Schechter sch(ps,ms,a);
					return sch(absmag); };
					
	// To integrate: int phi(M|z) dM
	double Integrate(double z);

  	// set the integration parameters (if want different from default)
  	void SetInteg(double Mmin,double Mmax, int npt=10000)
		{ Mmin_=Mmin; Mmax_=Mmax; npt_=npt; };

protected:
	LFParameters& lfpars_;
	double Mmin_,Mmax_;
	int npt_;
};

// returns phi(z) = [int phi(M|z) dM]*dV(z)
class SchechterZVol : public GenericFunc
{
public:
	// constructor
	SchechterZVol(LFParameters& lfpars, SimpleUniverse& su)
		: lfpars_(lfpars) , su_(su) {   
			zmin_=0, zmax_=6;
  			npt_=10000; };

	virtual double operator() (double z)
					{ 
					  SchechterM schM(lfpars_);
					  su_.SetEmissionRedShift(z);
					  double volel=su_.VolEl();
					  double val=schM.Integrate(z); //tm.Split();
					  //cout << "Time to integrate = "<<tm.PartialElapsedTimems()<<endl;
					  return val*volel; };

	// To integrate: int_z int_M phi(M|z) dV(z) dM dz
	double Integrate();

  	// set the integration parameters (if want different from default)
  	void SetInteg(double zmin,double zmax, int npt=10000)
		{ zmin_=zmin; zmax_=zmax; npt_=npt; };

protected:
	LFParameters & lfpars_;
	SimpleUniverse& su_;
	double zmin_,zmax_;
	int npt_;
};

class SchechterDist : public GenericFunc
// factor to multiply luminosity function types by:
// F = schA / (schE+schL+schS)
// schA = luminosity function of ALL types, at redshift z
// schE = luminosity function of early type, at redshift z
// schL = luminosity function of late type, at redshift z
// schS = luminosity function of starburst type, at redshift z
// Schechter disribution is: schDist = schWant * F
// schWant is designated by the value of type
// if type=1: schWant=schE
// if type=2: schWant=schL
// if type=3: schWant=schS
{
public:
  SchechterDist(Schechter& schA, Schechter& schE, Schechter& schL, Schechter& schS, int type)
	: schA_(schA) , schE_(schE) , schL_(schL) , schS_(schS) , type_(type)
	{
	if ( (type_>0.5)&(type_<1.5) )
			cout << "     type=EARLY"<<endl;
	if ( (type_>1.5)&(type_<2.5) )
			cout << "     type=LATE"<<endl;
	if ( (type_>2.5)&(type_<3.5) )
			cout << "     type=STARBURST"<<endl;
	}
	
    virtual double operator()(double x) {
        if ( (type_>0.5)&(type_<1.5) )
            return (schE_(x) * (schA_(x)/(schE_(x)+schL_(x)+schS_(x))) );
        if ( (type_>1.5)&(type_<2.5) )
            return (schL_(x) * (schA_(x)/(schE_(x)+schL_(x)+schS_(x))) );
        if ( (type_>2.5)&(type_<3.5) )
            return (schS_(x) * (schA_(x)/(schE_(x)+schL_(x)+schS_(x))) );
        else {
            stringstream ss;
            ss << type_;
            string msg = "ERROR! Schechter type (="+ ss.str() + ") not understood";
            throw ParmError(msg);
            }
        };
  
  double Integrate(double Mmin,double Mmax,int npt=100) // an integration function
	{
	 if(npt<1) npt = 100;
	 double perc=0.01, dlxinc=(Mmax-Mmin)/npt, dlxmax=10.*dlxinc; unsigned short glorder=4;
	 double sum = IntegrateFunc(*this,Mmin,Mmax,perc,dlxinc,dlxmax,glorder);
	 return sum;
	 }
  
  protected:
  Schechter& schA_, schE_, schL_, schS_; 
  int type_;
};

class SchechterMassDist : public AnyDataObj {
  friend class ObjFileIO<SchechterMassDist>;
public:
  SchechterMassDist(Schechter sch,double massmin,double massmax,int nbinmass);
  SchechterMassDist(void);
  virtual ~SchechterMassDist(void);

  int GetMassLim(double& massmin,double& massmax);
  int SetNgalLim(int ngalmax,int ngalmin=1,unsigned long nalea=10000);
  int GetNgalLim(int& ngalmax,int& ngalmin);
  int GetNgalLim(void) {return nvalngal_;}
  Schechter& GetSchechter(void) {return sch_;}

  inline int IndexFrNGal(int ngal) {
    int i = ngal-ngalmin_;
    if(nvalngal_<1 || i<0) return -1;
    if(i>=nvalngal_) return -2; else return i;
  }
  inline int NGalFrIndex(int i) {
    if(nvalngal_<1 || i<0 || i>=nvalngal_) return -1;
    return ngalmin_+i;
  }

  Histo GetHmDnDm(void) const;
  FunRan GetTmDnDm(void) const;

  Histo GetHisto(int i) const;
  FunRan GetFunRan(int i) const;

  double TirMass(int ngal);

  void Print(void);
  void PrintStatus(void);

  void WritePPF(string ppfname);
  void ReadPPF(string ppfname);

protected:
 void Delete(void);

 Schechter sch_;
 unsigned short sch_outvalue_;

 double massmin_,massmax_;   int nbinmass_;
 int ngalmin_,ngalmax_,nvalngal_;
 unsigned long ntrial_dir, ntrial_tab;

 Histo* hmdndm_;
 FunRan* tirhmdndm_;

 vector<Histo> hmass_;
 vector<FunRan> tmass_;
};

// Type ratios at redshift zero
class TypeRatio0 : public GenericFunc
{
public:
    TypeRatio0(Schechter& schA, Schechter& schE, Schechter& schL,
													Schechter& schS, int type=1)
		: schA_(schA) , schE_(schE) , schL_(schL) , schS_(schS) , type_(type)
			{  
			if(type_<1||type_>3)
				throw ParmError("ERROR! Schechter type not understood");
			};
		
	/**< Return typ ratio at absolute magnitude m */
    virtual double operator()(double m) {
			
				double norm=schA_(m)/(schE_(m)+schL_(m)+schS_(m));
				double enorm=schE_(m)*norm;
				double lnorm=schL_(m)*norm;
				double snorm=schS_(m)*norm;

				if ( (type_>0.5)&(type_<1.5) )
					return ( schE_(m) / (enorm+lnorm+snorm) );
				if ( (type_>1.5)&(type_<2.5) )
					return ( schL_(m) / (enorm+lnorm+snorm) );
				if ( (type_>2.5)&(type_<3.5) )
					return ( schS_(m) / (enorm+lnorm+snorm) );
		
			};
			
    void ChangeTypeTo(int type)	{
     
			if(type_<1||type_>3)
				throw ParmError("ERROR! Schechter type not understood"); 
			else
				type_=type; 
		};
  
  
protected:
    Schechter& schA_, schE_, schL_, schS_; 
    int type_;
};
	
// Type ratios as a function of redshift 
class TypeRatio : public GenericFunc
{
public:
    /** Constructor */
    TypeRatio(TypeRatio0& tr0)
		: tr0_(tr0) { };
		
		
	virtual void operator()(double m,double z,double& Fe,double& Fs,double& Fl) {
									
									tr0_.ChangeTypeTo(1);
									Fe=tr0_(m)*pow((1+z),-0.7);
								    tr0_.ChangeTypeTo(3);
									Fs=tr0_(m)*pow((1+z),0.7);
									tr0_.ChangeTypeTo(2);
									Fl=1-(Fe+Fs);
									
									};
							
protected:
    TypeRatio0& tr0_;
};

//-----------------------------------------------------------------------------------
bool IsCompatible(Schechter& sch1,Schechter& sch2,double eps=1.e-4);

//class SchechterDist : public GenericFunc
//{
//// luminosity distribution to select from
//// schtype = luminosity function of relevant type, at redshift z
//// schF = factor to multiply LF by (calculated from SchechterFactor)
//public:
//  SchechterDist(SchechterFactor& schF, Schechter& schtype)
//	: schF_(schF) , schtype_(schtype)
//	{
//	}
//	
//  virtual double operator()(double x)
//    {
//    return (schtype_(x)*schF_(x) );
//    }
//  
//  protected:
//  SchechterFactor schF_;
//  Schechter&  schtype_;
//};


//class SchechterIntegrator : public GenericFunc
////A simple numerical summation integration class
//{
//public:
//	SchechterIntegrator(GenericFunc& f, double xmin, double xmax, int Nstep=500) 
//    : f_(f) , xmin_(xmin) , xmax_(xmax) , Nstep_(Nstep)
//    {
//    }
//	
//	double Value()
//	{
//		
//		if (Nstep_ <= 0)
//			Nstep_ = 500;
//		
//		double dx = (xmax_ - xmin_) / (Nstep_-1);
//		
//		double x = xmin_;
//		double sum=0;
//		for (int i=1; i<Nstep_; i++, x += dx)
//			sum += f_(x);
//		
//		return sum * dx;
//	}
//	
//	
//protected:
//	GenericFunc& f_; 
//	double xmin_, xmax_;
//	int Nstep_;
//};
	
} // End namespace SOPHYA
#endif
