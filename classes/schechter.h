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
 * @todo create a wrapper class the other schechter classes can inherit 
 * idential methods from, eg the integration ones Integrate() and Integrate(double)
 *
 */

#ifndef SCHECHTER_SEEN
#define SCHECHTER_SEEN

#include "machdefs.h"
#include <stdlib.h>

// Sophya update v2.3 June 2013 replaces genericfunc with classfunc
#include "classfunc.h" 
//#include "genericfunc.h"
#include "geneutils.h"
#include "cosmocalcs.h"
#include "ctimer.h"

namespace SOPHYA {

// I think I put these here because of "error: return type **** is incomplete"
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
	 virtual void operator()(double z, double& ps, double& ms, double& a) {  	
	 
	        double minv=1000; int idmin=-1; double diff;
		    int nz = zgrid_.size();
			for (int i=0; i<nz; i++) {
				diff=z-zgrid_[i];
				
				if (std::abs(diff)<minv)
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
class Schechter : public ClassFunc1D 
{
public:
    /** Constructor 
        @param phistar  Schechter parameter \f$ \phi_* \f$
        @param Mstar    Schechter parameter \f$ M_* \f$
        @param alpha    Schechter parameter \f$ \alpha \f$                    */
    Schechter(double phistar, double Mstar, double alpha)
    : phistar_(phistar) , Mstar_(Mstar) , alpha_(alpha) , outvalue_(0) {
     // default setting used in integration of the Schechter function
     Mmin_=-24, Mmax_=-13; // units of "M-5log10h70"
     npt_=10000; }
    
    /** copy constructor */
    Schechter(Schechter& f);
                   
    /** default constructor */
    Schechter(void) : phistar_(0.) , Mstar_(0.) , alpha_(0.) , outvalue_(0) {};
    
    /** destructor */
    virtual ~Schechter(void) {};

    /** Set the output value, determine the operator() returns dn/dm or 
        m*dn/dm 
        @param outvalue    if=0 returns dn/dm, if>0 returns m*dn/dm           */
    void SetOutValue(unsigned short outvalue=0);
      
    /** Return the output value, whether operator() returns dn/dm or 
        m*dn/dm                                                               */
    unsigned short GetOutValue(void){ return outvalue_; };
      
    /** Reset the Schechter parameters to new values                          */
    void SetParam(double phistar,double Mstar,double alpha) {
        phistar_ = phistar; Mstar_ = Mstar; alpha_ = alpha; };
      
    /** Return current Schechter parameters                                   */
    void GetParam(double& phistar,double& Mstar,double& alpha) {
        phistar = phistar_; Mstar = Mstar_; alpha = alpha_; };
      
    /** Set the integration parameters (if want different from default)       */
    void SetInteg(double Mmin,double Mmax, int npt)
        { Mmin_=Mmin; Mmax_=Mmax; npt_=npt; };

    /** Return the LF value for absolute magnitude "absmag"                   */
    virtual double operator() (double absmag) const;

    /** Integrate Schechter function                                          */
    double Integrate(); // an integration function
    //double IntegrateQuick();// don't think this works?

    /** Return number of galaxies within a volume defined by a sky area in 
        steradians and between redshifts z1 and z2
        @param su         cosmology 
        @param skyarea    sky area in steradians
        @param z1         redshift 1
        @param z2         redshift 2                                          */
    double NGals(SimpleUniverse su, double skyarea, double z1, double z2);
    
    /** Return volume of universe within a sky area in steradians and between 
        redshifts z1 and z2
        @param su         cosmology 
        @param skyarea    sky area in steradians
        @param z1         redshift 1
        @param z2         redshift 2                                          */
    double Volume(SimpleUniverse su, double skyarea, double z1, double z2);

    /** Print Schechter function parameters 
        @param MU    units of Mstar
        @param PU    units of phistar                                         */
    virtual void Print(string MU, string PU);
    
    /** Print Schechter function parameters                                   */
    virtual void Print(void);
    
    /** Vary faint limit of Schechter function integration to check effect on 
        galaxy number density. Results are printed to the screen              */
    void CheckFaintLimit(Schechter&);

protected:
    double phistar_;             /**< Schechter function parameter            */
    double Mstar_;               /**< Schechter function parameter            */
    double alpha_;               /**< Schechter function parameter            */
    unsigned short outvalue_;    /**< output type dn/dm or dn/logm            */
    double Mmin_;                /**< Schechter function integration limit    */
    double Mmax_;                /**< Schechter function integration limit    */
    int npt_;                    /**< integration parameter                   */
    TVector<r_4> ngal_by_mpc3_;  /**< for debugging                           */
    TVector<r_4>  schmaxv_;      /**< for debugging                           */
};


/** @class SchechterVol class
  *
  * Returns the Schechter function \f$\phi(M|z)\f$ computed for a set of parameters: 
  * \f$\phi_*(z), M_*(z), \alpha(z) \f$ multiplied by the volume element at
  * redshift \f$z\f$.
  *
  * Returns \f$\phi(M|z)dV(z)\f$
  */
class SchechterVol : public ClassFunc1D, ClassFunc2D
{
public:
	/** Constructor
	    @param lfpars    luminosity function parameters
	    @param su        cosmology                                            */
	SchechterVol(LFParameters& lfpars, SimpleUniverse& su)
		: lfpars_(lfpars) , su_(su) {   
			Mmin_=-24, Mmax_=-13;// units of "M-5log10h70"
  			npt_=10000;  };

    /** This is defined to override the pure virtual function defined in ClassFunc1D
        otherwise SchechterVol is sometimes treated as an abstract class        */
    virtual double operator() (double) const { };
    
    /** Return Schechter function multiplied by volume element at redshift z
        @param absmag    absolute magnitude Schechter func returned at
        @param z         redshift of volume element and Schechter func        */
	virtual double operator() (double absmag, double z) const
		{ 
			double ps,ms,a;
			lfpars_(z,ps,ms,a);
			Schechter sch(ps,ms,a);
			su_.SetEmissionRedShift(z);
			double volel=su_.VolEl();
			return sch(absmag)*volel;
		}

	/** Integrate \f$ \int \phi(M|z) dV(z) dM \f$
	    @param z    redshift                                                  */
	double Integrate(double z);

  	/** Set the integration parameters (if want different from default)       
  	    @param Mmin    integration lower limit
  	    @param Mmax    integration upper limit
  	    @param npt     number of points to use in integration                 */
  	void SetInteg(double Mmin, double Mmax, int npt=10000)
		{ Mmin_=Mmin; Mmax_=Mmax; npt_=npt; };

protected:
	LFParameters& lfpars_;    /**< luminosity function parameters             */
	SimpleUniverse& su_;      /**< cosmology                                  */
	double Mmin_;             /**< integration lower limit                    */
	double Mmax_;             /**< integration upper limit                    */
	int npt_;                 /**< number of points to use in integration     */
};


/** @class SchechterM class
  *
  * Returns \f$ \phi(M|z) \f$
  */
class SchechterM : public ClassFunc1D, ClassFunc2D
{
public:
	/** Constructor
	    @param lfpars    luminosity function parameters */
	SchechterM(LFParameters& lfpars)
		: lfpars_(lfpars) {   
			Mmin_=-24, Mmax_=-13;// units of "M-5log10h70"
  			npt_=10000;  };
  			
    /** This is defined to override the pure virtual function defined in ClassFunc1D
        otherwise SchechterM is sometimes treated as an abstract class        */
    virtual double operator() (double) const { };

    /** Return \f$ \phi(M|z) \f$ 
        @param absmag    absolute magnitude
        @param z         redshift                                             */
	virtual double operator() (double absmag, double z) const { 
	
				  	double ps,ms,a;
					lfpars_(z,ps,ms,a);
					Schechter sch(ps,ms,a);
					return sch(absmag); };
					
	/** Integrate \f$ \int \phi(M|z) dM \f$
	    @param z    redshift                                                  */
	double Integrate(double z);

  	/** Set the integration parameters (if want different from default)       
  	    @param Mmin    integration lower limit
  	    @param Mmax    integration upper limit
  	    @param npt     number of points to use in integration                 */
  	void SetInteg(double Mmin,double Mmax, int npt=10000)
		{ Mmin_=Mmin; Mmax_=Mmax; npt_=npt; };

protected:
	LFParameters& lfpars_;    /**< luminosity function parameters             */
	double Mmin_;             /**< integration lower limit                    */
	double Mmax_;             /**< integration upper limit                    */
	int npt_;                 /**< number of points to use in integration     */
};


/** @class SchechterZVol class
  *
  * Returns \f$  \phi(z) = [\int \phi(M|z) dM]*dV(z) \f$
  */
class SchechterZVol : public ClassFunc1D
{
public:
	/** Constructor
	    @param lfpars    luminosity function parameters
	    @param su        cosmology                                            */
	SchechterZVol(LFParameters& lfpars, SimpleUniverse& su)
		: lfpars_(lfpars) , su_(su) {   
			zmin_=0, zmax_=6;
  			npt_=10000; };

    /** Returns \f$  \phi(z) = [\int \phi(M|z) dM]*dV(z) \f$ 
        @param z    redshift                                                  */
	virtual double operator() (double z) const { 
            SchechterM schM(lfpars_);
            su_.SetEmissionRedShift(z);
            double volel=su_.VolEl();
            double val=schM.Integrate(z); //tm.Split();
            //cout << "Time to integrate = "<<tm.PartialElapsedTimems()<<endl;
            return val*volel; };

	/** Integrate \f$ \int_z \int_M \phi(M|z) dV(z) dM dz \f$                 */
	double Integrate();

  	/** Set the integration parameters (if want different from default)       
  	    @param zmin    integration lower limit
  	    @param zmax    integration upper limit
  	    @param npt     number of points to use in integration                 */
  	void SetInteg(double zmin, double zmax, int npt=10000)
		{ zmin_=zmin; zmax_=zmax; npt_=npt; };

protected:
	LFParameters& lfpars_;    /**< luminosity function parameters             */
	SimpleUniverse& su_;      /**< cosmology                                  */
	double zmin_;             /**< integration lower limit                    */
	double zmax_;             /**< integration upper limit                    */
	int npt_;                 /**< number of points to use in integration     */
};


/** @class SchechterDist class
  *
  * Returns luminosity functions normalized to correct number density 
  *
  * Factor to multiply luminosity function types by in order to preserve 
  * correct number density is:
  * F = schA / (schE+schL+schS)
  * schA = luminosity function of ALL types, at redshift z
  * schE = luminosity function of early type, at redshift z
  * schL = luminosity function of late type, at redshift z
  * schS = luminosity function of starburst type, at redshift z
  *
  * ie returns schE^correct_n = schE * F
  *
  */
class SchechterDist : public ClassFunc1D
{
public:

    /** Constructor
        @param schA    luminosity function of all galaxy types
        @param schE    luminosity function of early/elliptical galaxy type
        @param schL    luminosity function of late/spiral galaxy type
        @param schS    luminosity function of starburst galaxy type
        @param type    which luminosity function to return                    */
    SchechterDist(Schechter& schA, Schechter& schE, Schechter& schL, Schechter& schS, int type)
    : schA_(schA) , schE_(schE) , schL_(schL) , schS_(schS) , type_(type) {
        if ( (type_>0.5)&(type_<1.5) )
            cout << "     type=EARLY"<<endl;
	    if ( (type_>1.5)&(type_<2.5) )
            cout << "     type=LATE"<<endl;
        if ( (type_>2.5)&(type_<3.5) )
            cout << "     type=STARBURST"<<endl; };
	
    /** Returns luminosity function normalized to correct number density      
        @param M    absolute magnitude                                        */
    virtual double operator()(double M) const {
        if ( (type_>0.5)&(type_<1.5) )
            return (schE_(M) * (schA_(M)/(schE_(M)+schL_(M)+schS_(M))) );
        if ( (type_>1.5)&(type_<2.5) )
            return (schL_(M) * (schA_(M)/(schE_(M)+schL_(M)+schS_(M))) );
        if ( (type_>2.5)&(type_<3.5) )
            return (schS_(M) * (schA_(M)/(schE_(M)+schL_(M)+schS_(M))) );
        else {
            stringstream ss;
            ss << type_;
            string msg = "ERROR! Schechter type (="+ ss.str() + ") not understood";
            throw ParmError(msg);
            }
        };
  
    /** Integrate the Schechter function
        @param Mmin    integration lower limit
  	    @param Mmax    integration upper limit
  	    @param npt     number of points to use in integration                 */
    double Integrate(double Mmin, double Mmax, int npt=100) const {
        if(npt<1) npt = 100;
        double perc=0.01, dlxinc=(Mmax-Mmin)/npt, dlxmax=10.*dlxinc; unsigned short glorder=4;
        double sum = IntegrateFunc(*this,Mmin,Mmax,perc,dlxinc,dlxmax,glorder);
        return sum; };
  
protected:
    Schechter& schA_;    /**< luminosity function of all galaxy types         */
    Schechter& schE_;    /**< luminosity function of early/elliptical types   */
    Schechter& schL_;    /**< luminosity function of late/spiral types        */
    Schechter& schS_;    /**< luminosity function of starburst galaxy types   */
    int type_;           /**< LF to return 1=EARLY, 2=LATE, 3=STARBURST       */
};


/** @class SchechterMassDist class
  *
  * Written by Christophe Magneville
  *
  */
class SchechterMassDist : public AnyDataObj {
    friend class ObjFileIO<SchechterMassDist>;
public:
    
    /** Constructor 
        @param sch         Schechter function
        @param massmin     minimum mass bin
        @param massmax     maximum mass bin
        @param nbinmass    number of mass bins                                */
    SchechterMassDist(Schechter sch, double massmin, double massmax, int nbinmass);
    
    /** Default constructor */
    SchechterMassDist(void);
    
    virtual ~SchechterMassDist(void) { Delete(); };

    /** Return the mass bin limits and number of bins. Number of bins is returned
        by the function
        @param massmin    minimum mass bin
        @param massmax    maximum mass bin                                    */
    int GetMassLim(double& massmin,double& massmax) {
        massmin = massmin_; massmax = massmax_; return nbinmass_; };
    
    /** Create histograms for drawing from ngalmin to ngalmax galaxies 
        @param ngalmax    maximum number of galaxies to draw
        @param ngalmin    minimum number of galaxies to draw
        @param nalea      number of trials to fill histogram*/
    int SetNgalLim(int ngalmax, int ngalmin=1, unsigned long nalea=10000);
    
    /** Return min and max number of galaxies to draw
        @param ngalmax    maximum number of galaxies to draw
        @param ngalmin    minimum number of galaxies to draw                  */
    int GetNgalLim(int& ngalmax, int& ngalmin) {
        ngalmax = ngalmax_;
        ngalmin = ngalmin_;
        return nvalngal_; };

    /** Return number of histograms created used to draw galaxies?            */
    int GetNgalLim(void) { return nvalngal_; }
    
    /** Return Schechter function                                             */
    Schechter& GetSchechter(void) { return sch_; }

    
    inline int IndexFrNGal(int ngal) {
        int i = ngal-ngalmin_;
        if(nvalngal_<1 || i<0) return -1;
        if(i>=nvalngal_) return -2; else return i;
        }
  
    inline int NGalFrIndex(int i) {
        if(nvalngal_<1 || i<0 || i>=nvalngal_) return -1;
        return ngalmin_+i;
        }

    Histo GetHmDnDm(void) const { return *hmdndm_; };
    
    // Needed to comment this out to compile mass2gal
    //FunRan GetTmDnDm(void) const { return *tirhmdndm_; };
    // Was getting this error:
    ///home/alex/Software/DirectSim/classes/schechter.h: In member function ‘SOPHYA::FunRan SOPHYA::SchechterMassDist::GetTmDnDm() const’:
    ///home/alex/Software/DirectSim/classes/schechter.h:575:34: error: return type ‘struct SOPHYA::FunRan’ is incomplete


    Histo GetHisto(int i) const;
    
    FunRan GetFunRan(int i) const;

    /** Draw mass */
    double TirMass(int ngal);

    void Print(void);
    
    void PrintStatus(void);

    void WritePPF(string ppfname);
    
    void ReadPPF(string ppfname);

protected:
    void Delete(void);
    Schechter sch_;
    unsigned short sch_outvalue_;
    double massmin_;
    double massmax_;   
    int nbinmass_;
    int ngalmin_;
    int ngalmax_;
    int nvalngal_;
    unsigned long ntrial_dir, ntrial_tab;
    Histo* hmdndm_;
    FunRan* tirhmdndm_;
    vector<Histo> hmass_;
    vector<FunRan> tmass_;
};

 
/** @class TypeRatio0 class
  *
  * Galaxy type fractions at redshift zero
  *
  */
class TypeRatio0 : public ClassFunc1D
{
public:
    /** Constructor 
        @param schA    luminosity function of all galaxy types
        @param schE    luminosity function of early/elliptical galaxy type
        @param schL    luminosity function of late/spiral galaxy type
        @param schS    luminosity function of starburst galaxy type
        @param type    which galaxy type fraction to return                   */
    TypeRatio0(Schechter& schA, Schechter& schE, Schechter& schL,
													Schechter& schS, int type=1)
		: schA_(schA) , schE_(schE) , schL_(schL) , schS_(schS) , type_(type)
			{  
			if(type_<1||type_>3)
				throw ParmError("ERROR! Schechter type not understood");
			};
		
	/** Return galaxy type fraction at absolute magnitude m */
    virtual double operator()(double m) const {
			
				double norm=schA_(m)/(schE_(m)+schL_(m)+schS_(m));
				double enorm=schE_(m)*norm;
				double lnorm=schL_(m)*norm;
				double snorm=schS_(m)*norm;

				if ( (type_>0.5)&(type_<1.5) )
					return ( schE_(m) / (enorm+lnorm+snorm) );
				else if ( (type_>1.5)&(type_<2.5) )
					return ( schL_(m) / (enorm+lnorm+snorm) );
				else if ( (type_>2.5)&(type_<3.5) )
					return ( schS_(m) / (enorm+lnorm+snorm) );
		        else {
		            stringstream ss;
		            ss << type_;
		            string emsg = "TypeRatio0::operator() ERROR! Schechter type = ";
		            emsg += ss.str() + " not understood";
		            throw ParmError(emsg);
		            }
			};
			
	/** Set galaxy type fraction to return */
    void ChangeTypeTo(int type)	{
     
			if(type_<1||type_>3)
				throw ParmError("ERROR! Schechter type not understood"); 
			else
				type_=type; 
		};
  
  
protected:
    Schechter& schA_;    /**< luminosity function of all galaxy types         */
    Schechter& schE_;    /**< luminosity function of early/elliptical types   */
    Schechter& schL_;    /**< luminosity function of late/spiral types        */
    Schechter& schS_;    /**< luminosity function of starburst galaxy types   */
    int type_;           /**< type fraction to return 1=EARLY,2=LATE,3=SBURST */

};
	

/** @class TypeRatio0 class
  *
  * Galaxy type fractions as a function of redshift 
  *
  */
class TypeRatio : public ClassFunc1D
{
public:
    
    /** Constructor 
        @param tr0    galaxy type fractions at redshift zero                  */
    TypeRatio(TypeRatio0& tr0) : tr0_(tr0) { };
    
    /** This is defined to override the pure virtual function defined in ClassFunc1D
        otherwise TypeRatio is sometimes treated as an abstract class        */
    virtual double operator() (double) const { };
		
    /** Return galaxy type fractions at absolute magnitude m, redshift z
        @param m    absolute magnitude
        @param z    redshift
        @param Fe   fraction of early/elliptical types
        @param Fs   fraction of starburst types
        @param Fl   fraction of late/spiral types */
	virtual void operator()(double m, double z, double& Fe, double& Fs, double& Fl) 
	const {
             tr0_.ChangeTypeTo(1);
             Fe = tr0_(m)*pow((1+z),-0.7);
			 tr0_.ChangeTypeTo(3);
             Fs = tr0_(m)*pow((1+z),0.7);
             tr0_.ChangeTypeTo(2);
             Fl = 1-(Fe+Fs); };
							
protected:
    TypeRatio0& tr0_;    /**< galaxy type fractions at redshift zero          */
};


/** Check that the parameters of two Schecter functions are within some 
    percentage precision
    @param sch1    Schechter function 1
    @param sch2    Schechter function 2
    @param eps     percentage precision to use                                */
bool IsCompatible(Schechter& sch1, Schechter& sch2, double eps=1.e-4);

//class SchechterDist : public ClassFunc1D
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
//  virtual double operator()(double x) const
//    {
//    return (schtype_(x)*schF_(x) );
//    }
//  
//  protected:
//  SchechterFactor schF_;
//  Schechter&  schtype_;
//};


//class SchechterIntegrator : public ClassFunc1D
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
