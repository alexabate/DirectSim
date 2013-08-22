/**
 * @file  gftdist.h
 * @brief Contains a series of classes for generating galaxy distributions;
 *        galaxy redshifts, types, magnitudes
 *
 * @todo <CODE>CumulDistZ::Output2File</CODE> <BR>
 *       When the calculated redshift CDF is output to a file the cosmology and
 *       luminosity function information used to calculate it should be included
 *       
 *
 * @todo <CODE>CumulDistZ::SetUp</CODE> <BR>
 *       When the calculated redshift CDF is read from a file into an 
 *       interpolation function, that function should set any values outside of
 *       zmin, zmax to zero.  Making this fix <I>shouldn't</I> affect the 
 *       results anyway though
 *
 * @todo <CODE>DrawZ::Draw</CODE> <BR>
 *       Does it matter that redshifts are drawn from the CDF on the open interval?
 *       This means that the exact values of zmin and zmax will <I>never</I> be
 *       drawn.
 *
 * @todo <CODE>DrawM::SetArray</CODE> <BR>
 *       Find out exactly why some of the interpolated magnitude CDF values have  
 *       the value nan
 *
 * @todo <CODE>DrawM::Draw</CODE> <BR>
 *       Replace with 2D interpolation table
 *
 * @author Alex Abate, Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2009
 * @date 2009
 *
 */
#ifndef  GFTDIST_H_SEEN
#define  GFTDIST_H_SEEN


// generic stuff
#include <iostream>
#include <fstream>

// sophya stuff
#include "array.h"
#include "genericfunc.h"
#include "pexceptions.h"
#include "randinterf.h" // This define RandomGenerator interface class
#include "sopnamsp.h"
#include "ctimer.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "swfitsdtable.h"

// stuff in classes dir
#include "schechter.h"
#include "cosmocalcs.h"
#include "mydefrg.h"  // definition of default RandomGenerator
#include "sinterp.h"
#include "geneutils.h"

namespace SOPHYA {

//******* CumulDistZ *********************************************************//

/** @class
  * CumulDistZ class
  *
  * Class to create cumulative distribution to draw redshifts from 
  *
  *
  */
class CumulDistZ {
public:
	/** Default constructor */
	CumulDistZ()
	: lfpars_(lfp_default_) , su_(su_default_) , zmin_(0) , zmax_(0)
	 , nptinteg_(0) , nptinterp_(0) {    };
	
	/** Constructor when calculating table
	    @param lfpars       class holding LF parameters and their evolution
	    @param su           class that does cosmological calculations
	    @param zmin         minimum redshift of cumulative distribution function
	    @param zmax         maximum redshift of cumulative distribution function
	    @param nptinteg     number of points to use in integration of ?
	    @param nptinterp    number of points to use in interpolation of ?     */
	CumulDistZ(LFParameters& lfpars,SimpleUniverse& su,double zmin=0.,
			double zmax=6.,int nptinteg=1000,int nptinterp=100,
				int prt=0);
				
	/** Constructor when reading table 
	    @param fname    FITS file to read from
        @param prt      printing level                                        */
	CumulDistZ(string fname,int prt=0);

	//METHODS//
	
	/** Returns value of cumultive redshift distribution function at redshift z
	    \f$ F_z(z)=\frac{\int_0^z\int_M_1^M_2 \phi(M,z')dV(z')dM}
	                {\int_0^z_{max}\int_M_1^M_2  \phi(M,z')dV(z')dM}          */
	virtual double operator()(double z)
		{ return schZint_(z); };

    /** Set up for reading distribution from a FITS bintable file 
        @param fname    FITS file to read from
        @param prt      printing level                                        */
	void SetUp(string fname, int prt=0);
	
	/** Set up for doing the actual calculation                
	    @param lfpars       class holding LF parameters and their evolution
	    @param su           class that does cosmological calculations
	    @param zmin         minimum redshift of cumulative distribution function
	    @param zmax         maximum redshift of cumulative distribution function
	    @param nptinteg     number of points to use in integration of \f$\phi(z)\f$
	    @param nptinterp    number of points to use in interpolation of \f$F_z(z)\f$*/
	void SetUp(LFParameters& lfpars, SimpleUniverse& su, double zmin=0.,
			double zmax=6., int nptinteg=1000, int nptinterp=100, int prt=0);
				
	/** Print n values of the distribution between zmin and zmax+0.1 as a check */
	void PrintDist(int nvals) {
		cout << "     CumulDistZ::PrintDist Printing "<< nvals <<" values from";
		cout << " the cumulative z dist ... "<<endl;
		//int ntot=zv_.size();
		//int nskip=ntot/nvals;
		double dz = (zmax_+0.1 - zmin_)/(nvals - 1);
		for (int i=0; i<nvals; i++) { 
		    double z = zmin_ + i*dz;
		    double val = schZint_(z);
		    cout << "     "<< z <<"  "<< val <<endl; 
		    }
		    //int j=i*nskip; 
			//cout << "     "<<zv_[j] <<"  "<<scv_[j]<<endl; }
		};

    /**< Output 2D array of cumulative redshift distribution to a FITS file
         binary table.  The zmin and zmax of the calculated distribution are
         included in the filename
        @param outfileroot  root name of FITS file to write to                */
	void Output2File(string outfileroot);
	
	void returnZminZmax(double& zmin, double& zmax)
	    { zmin = zmin_; zmax = zmax_; };

protected:
	LFParameters& lfpars_;      /**< class that stores the LF parameters and evolution */
	SimpleUniverse& su_;        /**< class that calculates cosmological quantities  */
	double zmin_;               /**< min redshift of cdf                      */
	double zmax_;               /**< max redshift of cdf                      */
	int nptinteg_;              /**< number of points in integration of \f$\phi(z)\f$*/
	int nptinterp_;             /**< number of points in interpolation of \f$F_z(z)\f$*/
	vector<double> zv_;         /**< grid of z values for \f$F_z(z)\f$        */
	vector<double> scv_;        /**< grid of F values (\f$F_z(z)\f$)          */
	SInterp1D schZint_;         /**< interpolation of \f$F_z(z)\f$            */
	SimpleUniverse su_default_;
	LFParameters lfp_default_;

};

//******* CumulDistM *********************************************************//

/** @class
  * CumulDistM class
  *
  * Class to create cumulative distribution to draw magnitudes from
  *
  * Returns the value of the cumulative magnitude distribution function given a
  * magnitude M and redshift z, \f$ F_M(M,z) \f$
  *
  */
class CumulDistM {
public:

	/** Default constructor */
	CumulDistM()
	: lfpars_(lfpars_default_) , su_(su_default_) , Mmin_(0) , Mmax_(0)
	{    };

	/** Constructor 
	@param lfpars   holds LF parameters and their evolution
	@param su       does cosmological calclations
	@param Mmin     sets lower integration bound
	@param Mmax     sets upper integration bound                              */
	CumulDistM(LFParameters& lfpars, SimpleUniverse& su, double Mmin=-24, 
			double Mmax=-13)
		: lfpars_(lfpars) , su_(su) , Mmin_(Mmin) , Mmax_(Mmax){    };
		
	// copy constructor
	//CumulDistM(CumulDistM const& a);

    /** Returns value of cumultive magnitude distribution function at magnitude M
        given redshift z
	    \f$ F_M(M,z)=\frac{\int_M_1^M_2 \phi(M',z)dV(z)dM'}
	                {\int_M_1^M_2  \phi(M',z)dV(z)dM'}                        */
	virtual double operator()(double m, double z) {
		SchechterVol schvol(lfpars_,su_);
		schvol.SetInteg(Mmin_,m);
		double top=schvol.Integrate(z);
		schvol.SetInteg(Mmin_,Mmax_);
		double bot=schvol.Integrate(z);
		return top/bot; 
		}
		
	/** Return magnitude integration limits                                   */
    void returnMminMmax(double& Mmin, double& Mmax)
	    { Mmin = Mmin_; Mmax = Mmax_; };

protected:
	LFParameters& lfpars_;
	SimpleUniverse& su_;
	double Mmin_,Mmax_;
	LFParameters lfpars_default_;
	SimpleUniverse su_default_;
};



//******* DrawZ **************************************************************//

/** @class
  * DrawZ class
  *
  * Class to draw redshifts from cumulative distribution 
  *
  */
class DrawZ {
public:

    /** Constructor
        @param cumz cumulative distribution to draw from
        @param npt  number of points to use in (reverse) interpolation table  */
	DrawZ(CumulDistZ& cumz, RandomGeneratorInterface& rg, int npt=10000)
	: rg_(rg) {
	
	    // Get zmin and zmax
	    cumz.returnZminZmax(zmin_, zmax_);
	    
	
	    // Reverse interpolation table
	    vector<double> zvals, cumzvals;
		double dz = (zmax_ - zmin_)/(npt - 1);
		for (int i=0; i<npt; i++) { 
			double z=zmin_ + i*dz;
			double cv=cumz(z);	
			zvals.push_back(z);
			cumzvals.push_back(cv);
			//if ( i<10 || i>npt-10)
			//    cout << z <<"  "<< cv <<endl;
			}

        // create interpolation
        // It won't matter if cumzvals is <0 or >1 because redshifts at these 
        // CDF values will never be drawn because a flat random number generator 
        // between 0 and 1 (open interval but shouldn't matter) is used
	    revCumZ_.DefinePoints(cumzvals,zvals,0.,1.,zvals.size()*4);

		};

	double Draw(){
	    // draw a flat random number between 0 and 1 (open interval)
	    double rn=rg_.Flat01();
	    double zs = revCumZ_(rn);
	    return zs;
	    };
	    
	void returnZminZmax(double& zmin, double& zmax)
	    { zmin = zmin_; zmax = zmax_; };

protected:
	SInterp1D revCumZ_;             /**< reverse interpolation of \f$F_z(z)\f$ */
	RandomGeneratorInterface& rg_;  /**< random number generation              */
	double zmin_;
	double zmax_;
};

//******* DrawM **************************************************************//
/** @class
  * DrawM class
  *
  * Class to draw magnitudes from cumulative distribution 
  *
  */
class DrawM {
public:

	/** Default constructor 
	    @param cumm cumulative magnitude distribution \f$ F_M(M,z) \f$ (interp table) */
	DrawM(CumulDistM& cumm);
	
	/** Constructor when calculating 2D array of cumulative magnitude distribution 
	    @param cumm     cumulative magnitude distribution \f$ F_M(M,z) \f$ (interp table)
	    @param rg       draws random numbers
	    @param mmin     minimum absolute mag value in 2D array
	    @param mmax     maximum absolute mag value in 2D array
	    @param zmin     minimum z value in 2D array
	    @param zmax     maximum z value in 2D array
	    @param nptz     number of z values in 2D array
	    @param nptm     number of absolute mag values in 2D array*/
	DrawM(CumulDistM& cumm, RandomGeneratorInterface& rg, double mmin=-24,
		double mmax=-13, double zmin=0., double zmax=6., int nptz=600, int nptm=600);
				
	/** Constructor when reading 2D array of cumulative magnitude distribution 
	    from a file 
	    @param fname    FITS filename containing pre-calculated 2D array of 
	                    F_M(M,z)                                              */
	DrawM(string fname);
	
	// Destructor
	//virtual ~DrawM();
	
	/** Calculate 2D array of cumulative magnitude distribution              */
	void SetUp(RandomGeneratorInterface& rg,double mmin=-24,
		double mmax=-13,double zmin=0.,double zmax=6.,int nptz=600,
			int nptm=600);

    /** Read 2D array of cumulative magnitude distribution from a file       */
	void SetUp(string infile);
	
    /** Output 2D array of cumulative magnitude distribution to a FITS file
        The magnitude and redshift range values are written in the filename
        @param outfileroot  root name of FITS file to write to                */
	void Output2File(string outfileroot);

    /** Store cumulative magnitude distribution values in an array           */
	void SetArray() {
		if (mmin_>0)
			throw ParmError("Min magnitude NOT set");
		if (zmin_<0)
			throw ParmError("Min redshift NOT set");
			
		Timer tm("DrawM::SetArray",false);
	    tm.Split();
		cout <<"     DrawM::SetArray "<<endl;
		// loop over magnitudes and redshifts
		for (int i=0; i<mv_.Size(); i++)
			for (int j=0; j<zv_.Size(); j++) { 

				double m=mmin_+i*dm_;
				double z=zmin_+j*dz_;			
				double cv=cumm_(m,z); // this is a calculation NOT an interpolation!
				int nantest=my_isnan(cv);
				if (nantest>0) { 
					cout << "     Found nan: z = "<< z <<", m = "<< m;
					cout << ", cv = "<< cv <<" setting cv->0" <<endl;
					cv=0;
					}
				
				// filling the arrays
				mv_(i)=m;
				zv_(j)=z;
				cumval_(i,j)=cv;  // dim1=m, dim2=z
				
				if (i<1 && (j<10 && j>8 || j<20 && j>18) ) {
		            tm.Split();
		            cout <<"     10 loops took "<< tm.PartialElapsedTimems() <<"ms"<<endl;
		            }

				}

		//cout << "Checking cumval m dist ... "<<endl;
		//for (int i=0; i<npt; i++)
		//	cout << mv_(0) <<"  "<<cumval_(0)<<endl;
		};

    /** Given a redshift z draw a magnitude according the the cumulative distribution
        in the array @cumval_
        THIS SHOULD BE REPLACED WITH A 2D INTERP TABLE                        */
	double Draw(double z);
	
	void returnZminZmax(double& zmin, double& zmax)
	    { zmin = zmin_; zmax = zmax_; };

protected:
	CumulDistM& cumm_;  /**< cumulative magnitude distribution \f$ F_M(M,z) \f$ (interp table) */
	RandomGeneratorInterface& rg_;/**< draws random numbers */
	TVector<r_8> mv_;
	TVector<r_8> zv_;
	TArray<r_8> cumval_;        /**< 2D array of cumulative mag distribution      */
	double mmin_;               /**< minimum absolute magnitude value in 2D array */
	double mmax_;               /**< maximum absolute magnitude value in 2D array */
	double dm_;                 /**< step in absolute magnitude in 2D array       */
	double zmin_;               /**< minimum z value in 2D array                  */
	double zmax_;               /**< maximum z value in 2D array                  */
	double dz_;                 /**< step in z in 2D array                        */
	CumulDistM cumm_default_;
	DR48RandGen rg_default_;
};

//******* DrawType ***********************************************************//
// class to draw type
class DrawType {
public:
	DrawType(TypeRatio tr,RandomGeneratorInterface& rg)
	: tr_(tr) , rg_(rg)
	{  };
	
	virtual int operator() (double m, double z) { 
	
	    double rn=rg_.Flat01();
		double Fe,Fs,Fl;
		// return type fractions at the redshift and magnitude
		tr_(m,z,Fe,Fs,Fl);
	    if (rn>=0 && rn<Fe)
		    return 1;
		else if (rn>=Fe && rn<(Fe+Fs))
			return 3;
		else if (rn>=(Fe+Fs))
			return 2;
		 };
	
protected:
    TypeRatio tr_;
	RandomGeneratorInterface& rg_;
};

//******* SimBaseCatalog *****************************************************//

/** @class
  * SimBaseCatalog class
  *
  * Simulates the base catalog of galaxy properties from redshift, absolute
  * magnitude and type distributions
  *
  */
class SimBaseCatalog {
public:

    /** Constructor */
	SimBaseCatalog(DrawZ& drz,DrawM& drm,DrawType& drt)
	: drz_(drz) , drm_(drm) , drt_(drt)
	{  
	
	    double zminRS, zmaxRS, zminMAG, zmaxMAG;
	    drz_.returnZminZmax(zminRS, zmaxRS);
        drm_.returnZminZmax(zminMAG, zmaxMAG);
        
        if ( zminRS<zminMAG || zmaxRS>zmaxMAG )
            throw ParmError("ERROR! Magnitude and redshift CDF's don't have matching redshift ranges");
        else {
            zmin_ = zminRS;
            zmax_ = zmaxRS;
	        }
	};

    /** Simulate ngal galaxies and output to FITS file */
	void DoSim(long ngal, string outfileroot);
	
	/** Simulate one galaxy  */
	void DoSim(double& z, double& am, double& typ)
	    {   z=drz_.Draw();
		    am=drm_.Draw(z);
		    typ=(double)drt_(am,z); };

protected:
	DrawZ& drz_;
	DrawM& drm_;
	DrawType& drt_;
	double zmin_;
	double zmax_;

};

//******* NGalZ **************************************************************//

/** @class
  * NGalZ class
  *
  * Class to calculate total number of galaxies between two redshifts and over
  * some solid angle
  *
  */
class NGalZ {
public:

    /** Constructor 
        @param lfpars   LF parameters as a function of z
        @param su       class that does cosmological calculations
        @param npt      number of points to integrate LF with                 */
	NGalZ(LFParameters& lfpars,SimpleUniverse& su,int npt=1000)
       : lfpars_(lfpars) , su_(su) , npt_(npt) {       };
 
    /** Returns the number of galaxies between two redshifts and some solid angle
        @param zmin minimum redshift
        @param zmax maximum redshift
        @param sa   solid angle in steradians                                 */
	virtual long operator()(double zmin,double zmax,double sa) {
	
	    // The below class returns phi(z) = [int phi(M|z) dM]*dV(z)
		SchechterZVol schZ(lfpars_,su_);
		schZ.SetInteg(zmin,zmax,npt_);
		double val=schZ.Integrate();
		long ntot=(long)val*sa;
		return ntot; 
		
		};

    /** Reset the number of points to integrate the LF with */
	void SetInteg(int npt)
		{ npt_=npt; };

protected:
	LFParameters& lfpars_;  /**< class that holds LF parameters and their evolution */
	SimpleUniverse& su_;    /**< class that does cosmological calculations    */
	int npt_;               /**< number of points to integrate LF with        */
};

//******* GalFlxTypDist ******************************************************//

/** @class
  * GalFlxTypDist class
  *
  * This class can be used to define a galaxy distribution dist(magnitude, type) 
  * and then generate (magnitude,type) pair values according to the defined 
  * distribution 
  *
  * @warning This class has now been superceded by the DrawM, DrawType classes
  */
class GalFlxTypDist {
public:
// Constructor, needs a RandomGenerator object
explicit GalFlxTypDist(RandomGeneratorInterface& rg, int prtlev=0);
explicit GalFlxTypDist(int prtlev=0);
virtual ~GalFlxTypDist();

// To be called once magnitude distribution for all 
// types have been defined through AddGalType() call 
// return true if OK, generate an exception (error) if problem
bool Check();

// Adds a galaxy type with the corresponding magnitude distribution
// See the .cc file for more information 
int AddGalType(GenericFunc& magf, double magmin, double magmax, double fr=1, 
               int nbinmag=20, int nstep=2000);
// Once all galaxy type-magnitude functions are added using AddGalType, pick a 
// galaxy type and mag at random according to these distributions using:
int GetGalaxy(int& typ, r_8& mag) const ;
// Added by AA Feb 2011:
int GetGalMag(r_8& mag) const ;
void AddSchechter(Schechter& t1,Schechter& t2,Schechter& t3);
void GetGalType(double mag,int& type);

// Return the number of different types defined
size_t NTypes() const { return v_mag_.size(); }
// Return the number of generated galaxies
size_t NGalGen() const { return totngal_; }

// Print some information about the galaxy type list
ostream& Print(ostream& os) const;

// THESE SHOULD BE PUT IN SEPARATE CLASS - they simulate the redshift distribution 
// according to the volume element 
int GetZDist(GenericFunc& dVdz, double zmin,double zmax, int nbin=100, int nstep=2000);
int GetZ(double&) const ;

void PrintDist(string);

// .... class variables 
RandomGenerator rgdefault_;
RandomGeneratorInterface& rg_;
int prtlev_ ;						// print/debug level 
vector< TVector<r_8> > v_mag_;		// Magnitude bins for each type
vector< r_8 > v_frac_;				// Fraction of total number of galaxies, for each type
vector< r_8 > v_sumfrac_;			// Sum of fractions for each type Sum [ 0 ... type] 
r_8 tot_frac_ ;						// Sum of type fraction for all types, should be = 1 when all types are defined
mutable size_t totngal_;			// total number of galaxies generated 
TVector<r_8> dVdzbins_;
//Schechter& t1_,t2_,t3_;
bool SchTypes_;
Schechter dt1,dt2,dt3;  // Default Schecter functions
Schechter t1_,t2_,t3_;  // Schechter function objects
};

/*! operator << overloading - Prints the list on the stream \b s */
inline ostream& operator << (ostream& s, GalFlxTypDist const & gfd)
  {  gfd.Print(s);  return(s);  }
  
  
//******* VolElement *********************************************************//
class VolElement : public GenericFunc
{
// dV/dz as a function of redshift
public:
  VolElement(SimpleUniverse& su, double sa)
	: su_(su) , sa_(sa)
	{
	}
	
  virtual double operator()(double z)
    {
	su_.SetEmissionRedShift(z);
	double dVdz=su_.HubbleLengthMpc()*pow(su_.RadialCoordinateMpc(),2)*su_.Gz(z)*sa_;
    return (dVdz);
    }
  
  protected:
  SimpleUniverse& su_;
  double sa_; 
};
} // End namespace SOPHYA
#endif
