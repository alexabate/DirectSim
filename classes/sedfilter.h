/**
 * @file  sedfilter.h
 * @brief Contains classes that perform calculations with SED and filter funcs
 *
 * @todo <CODE>SED::readSED</CODE> <BR> 
 *       Check if the interpolated flux value should be zero for wavelengths
 *       outside of those read in 
 *
 * @author Alex Abate and Reza Ansari
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2008
 * @date 2008
 *
 */
 
#ifndef SEDFILTER_H_SEEN 
#define SEDFILTER_H_SEEN 

#include "machdefs.h"
#include "sopnamsp.h"
#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
// Sophya update v2.3 June 2013 replaces genericfunc with classfunc
//#include "genericfunc.h"
#include "classfunc.h"
#include "tvector.h"

#include "reddening.h"
#include "sinterp.h"
#include "igm.h"
// check the below is needed
#include "matrix.h"

// ROOT libraries
// removed these root libraries that are not needed
#include "TMinuit.h"
#include "TGraph.h"
#include "TFile.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"
#include "TPrincipal.h"

namespace SOPHYA {

/*******************************************************************************
*                                                                              *
*                                SED CLASSES                                   *
*                                                                              *
*******************************************************************************/


/** @class
  * Spectral Energy Distribution (SED) class
  * 
  * SEDs are read in from a text file with 2 columns wavelength in meters and
  * flux.  The flux value at the given wavelength is returned by the class.
  * The SED can be optionally interpolated or reddened
  * 
  */
class SED : public ClassFunc1D
{
public:
    /** Default constructor */
    SED() : sed2_(sed2init_)
    { isRead_=false; isRedden_=false; isInterp_=false; _test_=1.; };
  
    /** Constructor
        To read in wavelength,flux values from a file, and return flux values
        based on the linear interpolation of numbers in that file */
    SED(string& filename, double xmin, double xmax, int npt=1024);
    
    /** Copy constructor */
    SED(const SED& sed);
    virtual SED& Set(const SED& a);
      
    /** Read in SED from file 
        @param fname    name of file to read from */
    void readSED(string& fname, double xmin, double xmax, int npt=1024);
      
    /** @warning code doesn't seem to like this! */
    virtual double operator()(double x) const
            { return returnFlux(x); };
            
    /** Set up interpolation of the SED */
    void doInterp(SED* sed2,double a=1.,double b=0.);
      
    /** Set up reddening of the SED */
    void doRedden(double EBmV=0., int law=0, double RvCard=3.5);
    
    /** To return the relevant flux value
        What is returned depends on isRedden and isInterp values
        Either: interpSED(), addReddening() or interpAddReddening() is called 
        @param lambda   rest-frame wavelength                                 */
    double returnFlux(double lambda) const;
    
    /** To return the reddened value 
        @param lambda   rest-frame wavelength                                 */
    double addReddening(double lambda) const;
    
    /** To return the interpolated value
        @param lambda   rest-frame wavelength                                 */
    double interpSED(double lambda) const;
    
    /** To return the reddened and interpolated value
        @param lambda   rest-frame wavelength                                 */
    double interpAddReddening(double lambda) const;
    
    /** return bool holding interpolation information*/
    bool returnInterpBool(){ return isInterp_; };
    
    /** return bool holding reddening information */
    bool returnReddenBool(){ return isRedden_; };
    
    /** return bool holding file-read information */
    bool returnReadBool(){ return isRead_; };
    
    
    double returntest(){ return _test_; };
    
protected:
    SInterp1D   sed_;           /**< The SED read in and to be used           */
    bool        isRead_;        /**< true if file has been read               */
    bool        isRedden_;      /**< true if SED has been reddened            */
    bool        isInterp_;      /**< true if SED has been interpolated        */
    SED*        sed2_;          /**< 2nd SED class to interpolate between with*/
    SED*        sed2init_;      /**< a default setting of sed2                */
    double      a_,b_;          /**< how to interpolate between sed and sed2  */
    double      EBmV_,RvCard_;
    int         law_;
    double      _test_;
};


/** @class
  * SED interpolation class
  * 
  * Can linearly interpolate between 2 SED objects as specifed with the 
  * arguments. For example if you want to interpolate and make a new SED where
  * 10% of it is "sed1" and 90% of it is "sed2" set "mult1"=0.1 and "mult2"=0.9
  * 
  * @note This class is possibly redundant now this functionality is 
  * with SED class itself
  * 
  */
class SEDInterp : public ClassFunc1D
{
public:
    
    SEDInterp(SED& sed1, SED& sed2, double mult1, double mult2)
    : sed1_(sed1), sed2_(sed2), mult1_(mult1), mult2_(mult2) {
        double eps=1e-6;
        if ( std::abs(mult1+mult2-1.)>eps)
            throw ParmError("ERROR! SED factors must sum to 1");
    };
     
    /** returns interpolated SED 
        @param lambda   rest-frame wavelength                                 */
    virtual double operator()(double lambda) const
            {  return (mult1_*sed1_(lambda)+mult2_*sed2_(lambda));  };

protected:
    SED& sed1_;      /**< SED 1                                               */
    SED& sed2_;      /**< SED 2                                               */
    double mult1_;   /**< fraction of SED 1 in new SED                        */
    double mult2_;   /**< fraction of SED 2 in new SED                        */
};


/** @class SEDzFilterProd
  * Multiply SED and Filter together to create an integrand for the k-correction
  * 
  * returns: redshifted SED multiplied by the filter transmission and wavelength
  *          SED(lambda/(1+z))*Filter(lambda)*lambda
  *
  * BECAUSE for AB magnitudes, and SED in units of wavelength:
  *
  * mAB(z) = -2.5log10( int SED(lambda/(1+z))*Filter(lambda)*lambda dlambda /  
  *                                 int Filter(lambda)*lambda^-1 dlambda ) -51.6
  *
  * and when computing Kcorrection the int Filter(lambda)*lambda^-1 dlambda )-51.6
  * cancel because:
  * K12 = mAB(z) - mAB(z=0)
  */
class SEDzFilterProd : public ClassFunc1D
{
public:
    /** Constructor 
        @param f    SED
        @param g    filter
        @param z    redshift of galaxy SED                                           */
    SEDzFilterProd(ClassFunc1D& f, ClassFunc1D& g, double z=0.)
    : sed_(f) , filt_(g) , z_(z) {  };
    
    /** Returns redshifted SED multiplied by the filter transmission and wavelength,
        \f$ SED(\lambda/(1+z))*Filter(\lambda)*\lambda \f$ where \lambda is the 
        observed wavelength    
        @param lambda   observed wavelength                                   */
    virtual double operator()(double lambda) const {
        double lambdaE = lambda/(1+z_);
        return (sed_(lambdaE)*filt_(lambda)*lambda);
        }  
// *lambda if SED in units of wavelength and magnitudes have dnu/nu definition
// no lambda factor if SED in units of wavelength and magnitudes have dnu definition
// /lambda if SED in units of frequency and magnitudes have dnu/nu definition
// /lambda^2 if SED in units of frequency and magnitudes have dnu definition

protected:
  ClassFunc1D& sed_;          /**< SED                                        */
  ClassFunc1D& filt_;         /**< filter function                            */
  double z_;                  /**< redshift of galaxy SED                     */
};


/** @class SEDGOODSRedfix class
  * 
  *
  * Host galaxy reddening via Cardelli or Calzetti law: Reddens SED 
  * Compensates for the fact that LF used to simulate galaxy luminosity
  * distributions will NOT have been corrected for host galaxy reddening.
  * 
  * Explanation:
  * To "renormalize" the magnitudes so that the "correct" number of 
  * simulated galaxies is simulated in the observed GOODS R-band, make the SED 
  * brighter by the amount the galaxy is extinct at the rest-frame wavelength 
  * corresponding to the observed R-band. This will look like:
  * 
  * \f$ SED^{fix}=SED^{unred}(\lambda)/10^{-0.4k(\lambda_{eff})E(B-V)}\f$
  * where \f$ \lambda_{eff}=6510./(1+z) \f$
  * (*dividing* by \f$ 10^{-0.4k(\lambda_{eff})E(B-V)} \f$ means we are making it *brighter*)
  * 
  * here \f$\lambda\f$ is the effective wavelength of the R-band that the GOODS LF 
  * was selected in (assuming that the SEDs are in units of Angstrom). If you use
  * an LF that is derived from another observed band, for example Ks selected, 
  * then you should change the \f$\lambda_{eff}\f$ to be the effective wavelength of the Ks 
  * band instead ( 21620A for the GOODS).
  * 
  * Important note: this assumes that the SED^unred(lambda) and \f$k(\lambda)\f$ are in 
  * "rest-frame" units, i.e., not redshifted. If they are redshifted, then you 
  * should use \f$\lambda_{eff}=6510.\f$
  * 
  * Then to actually redden the SED:
  *
  * \f$ SED^{red}=SED^{fix}(\lambda)10^(-0.4k(\lambda)E(B-V)) \f$ 
  *
  */
class SEDGOODSRedfix: public ClassFunc1D
{
public:

    /** Constructor 
        @param S        SED
        @param z        redshift of SED 
        @param EBmV     E(B-V) value of galaxy
        @param law      Reddening law to use (=0 Cardelli,=1 Calzetti)
        @param RvCard   Rv parameter of Cardelli law
        @param lam_eff0 effective wavelength (for correcting GOODS reddening) */
    SEDGOODSRedfix(ClassFunc1D& S,double z,double EBmV=0.,int law=0, 
            double RvCard=3.5,double lam_eff0=6510e-10)
        : sed_(S) , z_(z) , EBmV_(EBmV) , law_(law) , RvCard_(RvCard) , 
                lam_eff0_(lam_eff0){    };
                
    /** Return reddened SED after correcting for GOODS reddening
        @param lambda   rest-frame wavelength                                 */
    virtual double operator()(double lambda) const {
        double lameff=lam_eff0_/(1+z_);
    
        Reddening red;
        double k,kfix;
        if (law_<1) {
            kfix=red.Cardelli(lameff,RvCard_);
            k=red.Cardelli(lambda,RvCard_);
            }
        if (law_>0) {
            kfix=red.Calzetti(lameff);
            k=red.Calzetti(lambda);
            }
    
        return (sed_(lambda)*pow(10,-0.4*k*EBmV_)/
                    pow(10,-0.4*kfix*EBmV_));
        };

protected:
  ClassFunc1D& sed_;       /**< SED to apply host galaxy redening to          */
  double z_;               /**< redshift of the SED                           */
  double EBmV_;            /**< E(B-V) host galaxy extinction value           */
  int law_;                /**< law_=0 uses Cardelli law, law_=1 uses Calzetti law */
  double RvCard_;          /**< Rv value                                      */
  double lam_eff0_;     //
};


/** @class SEDRedden class
  * 
  * Host galaxy reddening via Cardelli or Calzetti law: Reddens SED 
  * Reddens SED according to Reddening law and E(B-V) value
  * Default reddening law is Cardelli
  * Default Cardelli Rv=3.5
  */
class SEDRedden : public ClassFunc1D
{
public:

    /** Constructor 
        @param S        SED
        @param EBmV     E(B-V) value of galaxy
        @param law      Reddening law to use (=0 Cardelli,=1 Calzetti)
        @param RvCard   Rv parameter of Cardelli law                          */
    SEDRedden(ClassFunc1D& S, double EBmV=0., int law=0, double RvCard=3.5)
    : sed_(S) , EBmV_(EBmV) , law_(law) , RvCard_(RvCard) { }
    
    /** Return reddened SED 
        @param lambda   rest-frame wavelength                                 */
    virtual double operator()(double lambda) const {
        
        Reddening red;
    
        double k;
        if (law_<1)
            k=red.Cardelli(lambda,RvCard_);
        if (law_>0)
            k=red.Calzetti(lambda);
    
        return (sed_(lambda)*pow(10,-0.4*k*EBmV_));
        } ;

protected:
  ClassFunc1D& sed_;      /**< SED to redden                                  */
  double EBmV_;           /**< E(B-V) extinction value                        */
  int law_;               /**< if law_=0 Cardelli law, if law_=1 Calzetti law */
  double RvCard_;         /**< Rv value                                       */
};


/** @class SEDMadau class
  * 
  * Add Madau absorption
  */
class SEDMadau : public ClassFunc1D
{
public:
<<<<<<< HEAD
    /** Constructor */
    SEDMadau(GenericFunc& S, double zSED, bool isLyC=true)
    : sed_(S) , zSED_(zSED) , isLyC_(isLyC) { };
=======
    /** Constructor 
        @param S       SED to add Madau absorption to 
        @param zSED    redshift of galaxy SED                                 */
    SEDMadau(ClassFunc1D& S, double zSED)
    : sed_(S) , zSED_(zSED) { };
>>>>>>> b174eb658571fc61c81d79c7fa2a0db330092148
    
    /** Return SED with Madau absorption applied
        @param lambda   rest-frame wavelength */
<<<<<<< HEAD
    virtual double operator()(double lambda) {
        Madau madau(5,isLyC_);
=======
    virtual double operator()(double lambda) const {
        Madau madau;
>>>>>>> b174eb658571fc61c81d79c7fa2a0db330092148
        // THINK THIS SHOULD BE RETURNING THE REST FRAME TRANSMISSION
        double trans = madau.returnRestFrameTransmission(lambda, zSED_);
        return (sed_(lambda)*trans);
        };

protected:
<<<<<<< HEAD
    GenericFunc& sed_;      /**< SED to add Madau absorption to               */
    double zSED_;           /**< redshift of SED                              */ 
    bool isLyC_;            /**< include Lyman continuum                      */       
=======
    ClassFunc1D& sed_;       /**< SED to add Madau absorption to              */
    double zSED_;            /**< redshift of SED                             */        
>>>>>>> b174eb658571fc61c81d79c7fa2a0db330092148
};


/** @class SEDIGM class
  * 
  * Add IGM absorption from a particular line of sight
  */
class SEDIGM : public ClassFunc1D
{
public:
    /** Constructor 
        @param S       SED to add IGM absorption to
        @param T       IGM transmission
        @param zSED    redshift of galaxy SED                                 */
    SEDIGM(ClassFunc1D& S, ClassFunc1D& T, double zSED)
    : sed_(S) , igm_(T) , zSED_(zSED) { };
    
    /** Return SED with IGM absorption applied
        @param lambda   restframe wavelength                                  */
    virtual double operator()(double lambda) const {
        double lambdaObs = lambda*(1+zSED_);
        double trans = igm_(lambdaObs);
        return (sed_(lambda)*trans);
        };

protected:
    ClassFunc1D& sed_;            /**< SED to add IGM absorption to           */
    ClassFunc1D& igm_;            /**< IGM absorption along line of sight     */
    double zSED_;                 /**< redshift of SED                        */        
};


/** ReadSedList class
  * 
  * Class to read in SEDs from a list of filenames
  *
  */
class ReadSedList {
public:
    /** Constructor, finds the file sedFile and counts number of SEDs inside 
        @param sedFile  filename containing list of SEDs 
        @param prt      print level, if prt>0 extra statements are printed */
    ReadSedList(string sedFile, int prt=0);
    
    /** Read environment variable $SEDLOC */
    string getSedDirEnviromentVar();
    
    /** Counts SEDs in file */
    void countSeds(string sedFileFullPath);
    
    /** Main program, reads SEDs from file and arranges them in a vector of
        pointers pointing to each SED object 
        @param lmin minimum wavelength of SED in meters
        @param lmax maximum wavelength of SED in meters */
<<<<<<< HEAD
    void readSeds(double lmin=5e-8,double lmax=2.5e-6);
=======
    void readSeds(double lmin=1e-7, double lmax=1e-6);
>>>>>>> b174eb658571fc61c81d79c7fa2a0db330092148
    
    /** If interpolating between the SEDs call this method straight after
        readSeds()
        @param nInterp number of SEDs to interpolate between each SED */
    void interpSeds(int nInterp);
    
    /** If reddening all the SEDs call this method straight after #readSeds(), 
        and after #interpSeds() if interpolating.
        Currently only applies Cardelli law, and can't limit max reddening for
        elliptical galaxies 
        @param nStepRed number of times to redden 
        @param redMax   maximum limit to redden to (even steps between 0 and redMax */
    void reddenSeds(int nStepRed, double redMax);
    
    /** Write contents of sedArray to a file */
<<<<<<< HEAD
    void writeSpectra(string outFile,double lmin=5e-8,double lmax=2.5e-6,
=======
    void writeSpectra(string outFile, double lmin=1e-7, double lmax=1e-6,
>>>>>>> b174eb658571fc61c81d79c7fa2a0db330092148
                                                                 int nl=1500);
                                                                 
    /** Return sedArray */
    vector<SED*> getSedArray()
                { return sedArray_; };
    
    /** Return total number of SEDs (could be >= to nsed depending if interpolation
        was done or reddening was applied) */
    int getNTot()
                { return ntot_; };
                
    /** Return number of SEDs read in from the initial file */
    int getNSed()
                { return nsed_; };
                
    /** Reorder SEDs so that interpolated SEDs lie between the template they
        were interpolated from */
    void reorderSEDs();
    
    /** Return a vector of the SED filenames */
    vector<string> returnSedFilenames(){ return sedFiles_; };
                   
private:
    int prt_;                   /**<  print level                             */
    string sedDir_;             /**<  path of SED list file                   */
    string sedFileFullPath_;    /**<  full path and filename to SED list      */
    int nsed_;                  /**<  number of SEDs in the file sedFile      */
    int ntot_;                  /**<  total number of SEDs after interpolation /
                                      and reddening                           */
    vector<SED*> sedArray_;     /**<  pointer to ntot_ SED objects            */
    vector<string> sedFiles_;   /**<  vector of SED file names                */
    bool isInterp_;             /**<  if SEDs are interpolated                */
    bool isRedden_;             /**<  if SEDs are reddened                    */
};


/*******************************************************************************
*                                                                              *
*                              FILTER CLASSES                                  *
*                                                                              *
*******************************************************************************/

/** @class Filter class 
  * To load in Filter transmission functions:
  * two columns: observed wavelength (m) ; transmission probability (between 0 and 1) 
  */
class Filter : public SInterp1D
{
public:
    /** Default constructor */
    Filter(){};
    
    /** Constructor 
        @param filename     file to read transmission function from
        @param lmin         min observed wavelength value
        @param lmax         max observed wavelength value
        @param nComments    number of comment lines the filter file has at start
        @param zero_outside value of interpolation outside lmin,lmax is zero if true
        @param npt          number of points to do interpolations from        */
    Filter(string& filename, double lmin, double lmax, int nComments=0, 
                                    bool zero_outside = true, int npt = 1024);
                                      
   /** Read in filter from file 
       @param filename     file to read transmission function from
       @param lmin         min observed wavelength value
       @param lmax         max observed wavelength value
       @param nComments    number of comment lines the filter file has at start
       @param zero_outside value of interpolation outside lmin,lmax is zero if true
       @param npt          number of points to do interpolations from         */
    void readFilter(string& filename, double lmin, double lmax, int nComments=0, 
                                bool zero_outside = true, int npt=1024);
};


/** @class BlueShiftFilter
  *
  */
class BlueShiftFilter : public ClassFunc1D
{
public:
    
    /** Constructor */
    BlueShiftFilter(Filter& g, double z=0.)
    : filt_(g) , z_(z) { }
    
    /** returns the value of the filter transmission at the rest-frame 
        wavelength of the object */
    virtual double operator()(double lambda) const {
        double lambdaRF = lambda*(1+z_);
        return (filt_(lambdaRF));
        };

protected:
    Filter& filt_;    /**< class holding the filter function */
    double z_;        /**< redshift of the object */
};


/** @class FilterProd class
  * To multiply filter transmission by 1/lambda
  *
  * returns: Filter(lambda)*lambda^-1
  * BECAUSE for AB magnitudes, and SED in units of wavelength:
  * mAB(z) = -2.5log10( int SED(lambda/(1+z))*Filter(lambda)*lambda dlambda /  
  *                                    int Filter(lambda)*lambda^-1 dlambda )-51.6
  *
  * and when computing Kcorrection the int Filter(lambda)*lambda^-2 dlambda )-51.6
  * cancels because:
  *
  * K12 = mAB(z) - mAB(z=0)
  *
  * spectral flux density: W/m^2/Hz int across bandwidth but still want flux to 
  * be in units of W/m^2/Hz because magnitude in band X is defined flux density 
  * at a particular wavelength
  */
class FilterProd : public ClassFunc1D
{
public:
    
    /** Constructor */
    FilterProd(Filter& g)
    : filt_(g){ };
    
    /** returns the value of the filter transmission at #lambda, divided by 
        #lambda */
    virtual double operator()(double lambda) const {
        return (filt_(lambda)/lambda);//********* DOUBLE CHECK THIS *********//
        }  

protected:
  Filter& filt_;    /**< class holding the filter function */
};


/** @class FilterXLambda class
  * To multiply filter transmission value by the wavelength
  *
  */
class FilterXLambda : public ClassFunc1D
{
public:
    
    /** Constructor */
    FilterXLambda(ClassFunc1D& g)
    : filt_(g){ }
    
    /** returns the value of the filter transmission at #lambda, multiplied by 
        #lambda */
    virtual double operator()(double lambda) const {
        return (filt_(lambda)*lambda);//********* DOUBLE CHECK THIS *********//
        }  

protected:
  ClassFunc1D& filt_;          /**< class holding the filter function         */
};


/** @class FilterProdProd class
  * To multiply filter with 1/lambda^2
  *
  */
class FilterProdProd : public ClassFunc1D
{
public:
    FilterProdProd(Filter& g)
    : filt_(g) {};
    
    virtual double operator()(double lambda) {
        return (filt_(lambda)/(lambda*lambda));
        }  

protected:
    Filter& filt_; 
};


//--- Does int F(lambda)/lambda dlambda or SED(lambda)*F(lambda) or whatever
// Simple summing integration
class FilterIntegrator : public ClassFunc1D
{
public:
  FilterIntegrator(ClassFunc1D& f, double xmin, double xmax, int Nstep=500) 
    : f_(f) , xmin_(xmin) , xmax_(xmax) , Nstep_(Nstep) { }
    
   /** This is defined to override the pure virtual function defined in ClassFunc1D
        otherwise FilterIntegrator is sometimes treated as an abstract class        */
    virtual double operator() (double) const { };
    
   double Value() {
        if (xmin_>=xmax_){
            string emsg="FilterIntegrator::Value() integral limits don't make sense ";
            throw out_of_range(emsg);
            }
                    
            if (Nstep_ <= 0)
                Nstep_ = 500;

            double dx = (xmax_ - xmin_) / (Nstep_-1);
            double x = xmin_;
            double sum=0;
            for (int i=1; i<Nstep_; i++, x += dx)
                sum += f_(x);
            
            return sum * dx;
        }

protected:
  ClassFunc1D& f_; 
  double xmin_, xmax_;
  int Nstep_;
};


/** ReadFilterList class
  * 
  * Class to read in filters from a list of filenames
  *
  */
class ReadFilterList {
public:
    /** Constructor, finds the file sedFile and counts number of filters inside 
        @param sedFile  filename containing list of filters 
        @param prt      print level, if prt>0 extra statements are printed */
    ReadFilterList(string filterFile, int prt=0);
    
    /** Read environment variable $FILTLOC */
    string getFilterDirEnviromentVar();
    
    /** Counts Filters in file */
    void countFilters(string filterFileFullPath);
    
    /** Main program, reads filters from file and arranges them in a vector of
        pointers pointing to each Filter object 
        @param lmin minimum wavelength of filter in meters
        @param lmax maximum wavelength of filter in meters */
<<<<<<< HEAD
    void readFilters(double lmin=5e-8,double lmax=2.5e-6);
    
    /** Write contents of filterArray to a file */
    void writeFilters(string outFile,double lmin=5e-8,double lmax=2.5e-6,
=======
    void readFilters(double lmin=1e-7, double lmax=1e-6);
    
    /** Write contents of filterArray to a file */
    void writeFilters(string outFile, double lmin=1e-7, double lmax=1e-6,
>>>>>>> b174eb658571fc61c81d79c7fa2a0db330092148
                                                                 int nl=1500);
                                                                 
    /** Return filterArray */
    vector<Filter*> getFilterArray()
                { return filterArray_; };
    
    /** Return total number of filters */
    int getNTot() { return ntot_; };
                
private:
    int prt_;                       /**< print level */
    string filterDir_;              /**< path of filter list file */
    string filterFileFullPath_;     /**< full path and filename to filter list */
    int ntot_;                      /**< total number of filters */
    vector<Filter*> filterArray_;   /**< pointer to ntot_ Filter objects */
};


/*******************************************************************************
*                                                                              *
*                              PCA RELATED CLASSES                             *
*                                                                              *
*******************************************************************************/



/** SEDCovMat class
  * 
  * Class to make SED covariance matrix from input point to SED class array
  * Probably defunct now using ROOT TPrincipal for calculating PCA
  *
  */
class SEDCovMat
{
public:
    SEDCovMat(double lmin,double dl,int nl)
    : lmin_(lmin) , dl_(dl) , nl_(nl) {  };

    TArray<double> MakeCovMat(SED **sedarray,int nsed,
                        TVector<double>& meanSED);
    TArray<double> TransposeMult(TArray<double> SEDmatrix);
    
    void SetWavelengths(double lmin,double dl,int nl)
        { lmin_=lmin; dl_=dl; nl_=nl; };
    
protected:
    double lmin_,dl_;
    int nl_;

};



/** TemplatePCA class
  * 
  * Class to calculate principal components of SED templates 
  *
  */
class TemplatePCA {
public:

    /** Constructor for calculating eigenvectors and eigenvalues from the spectra
        in sedArray */
    TemplatePCA(vector<SED*> sedArray, double lmin=5e-8, double lmax=2.5e-6, 
                                                                 int nl=1500);
    
    /** Constructor for projecting the spectra in sedArray onto the eigenvectors 
        read from file eigVectFile */
    TemplatePCA(vector<SED*> sedArray,string eigVectFile, double lmin=5e-8, 
                                                double lmax=2.5e-6,int nl=1500);
                                                                
    /** should add argument to here to reflect different normalization choices*/
    void normalizeSpectra();
    
    /** Find mean of spectra */
    void meanSpectra();
    
    /** add each spectrum to the data matrix */
    void addSpectra();
    
    /** Calculate covariance matrix of the data */
    void calculateCovarianceMatrix();
    
    /** Normalize covariance matrix by its trace */
    void normalizeCovarianceMatrix();
    
    /* Decompose into eigenvalues and eigenvectors */
    void doPCA();
    
    /* Copy the covariance matrix into a ROOT object */
    TMatrixD copyCovMatrixToROOT();
    
    /** Project each spectrum onto eigenvectors 1:nEigKept and find the 
        eigenvalues of the projected spectrum.
        
        Basically this is a loop over reconstructSpectrum method.
        Fills class matrices: reconstructedSpectra_ and eigenvalsProjSpec_.
        
        @param nEigKept number of eigenvalues to project onto */
    void reconstructSpectra(int nEigKept);
    
    /** Project spectrum iSpectrum onto eigenvectors 1:nEigKept 
        @param nEigKept number of eigenvalues to project onto
        @param iSpectrum spectrum to project
        @param reconstructedSpectrum reconstructed spectrum */
    TVector<double> reconstructSpectrum(int nEigKept,int iSpectrum,
                                TMatrix<double>& reconstructedSpectrum);   
                                
    /** Find out how good the fit of the reconstructed to the true spectrum is */
    TVector<double> fitSpectra();                                            
                              
    // These just return stuff as a SOPHYA object
    /** Return meanValues as SOPHYA object */                
    //TVector<double> getMeanValues();
    /** Return covariance matrix as SOPHYA object */
    TMatrix<double> getCovMatrix();
    /** Return eigenvalues as SOPHYA object */
    TVector<double> getEigenValues();
    /** Return eigenvectors as SOPHYA object */
    TMatrix<double> getEigenVectors();
    
    // These just return things
    /** Return reconstructed spectra */
    TMatrix<double> returnRecSpectra(){ return reconstructedSpectra_; };
    /** Return eigenvalues of reconstructed spectra */
    TMatrix<double> returnEigValsProjSpec(){ return eigenvalsProjSpec_; };
    /** Return spectra normalization values */
    TVector<double> returnNormValues(){ return normValues_; };
    
    // Methods to write stuff to a file
    /** Write spectra normalization values to a file */
    void writeNormValues(string outFile);
    /** Write spectra mean values to a file */
    void writeMeanValues(string outFile);
    /** Write covariance matrix to a file */
    void writeCovMatrix(string outFile);
    /** Write data matrix to a file */
    void writeDataMatrix(string outFile);
    /** Write eigenvectors to a file */
    void writeEigenVectors(string outFile);
    /** Write eigenvalues to a file */
    void writeEigenValues(string outFile);
    /** Write eigenvectors and eigenvalues to a file */
    void writeEigenValVecs(string outFile);
    /** Write eigenvalues of projected spectrum to a file */
    void writeEigenValsOfProjSpec(string outFile);
    /** Write reconstructed spectra to a file */
    void writeRecSpec(string outFile);
  
    /** Read eigenvectors from a file */
    TMatrixD readEigenVectors(string inFile);
    
    // Potentially add method(s) to calculate Karnhunen-Loeve angles

protected:
    vector<SED*> sedArray_;             /**< array of SED objects */
    double lmin_;                       /**< min wavelength */
    double dl_;                         /**< step in wavelength */
    int nl_;                            /**< number of wavelengths to in SED data matrix*/
    int nsed_;                          /**< number of SEDs */
    TMatrix<double> dataMatrix_;        /**< matrix of normalized, mean subbed SEDs */
    //TPrincipal* principal_;             /**< class that does PCA calc */
    TVector<double> normValues_;        /**< normalization of each spectrum */
    TVector<double> meanValues_;        /**< mean values of each wl bin (over SEDs)*/
    TMatrix<double> covMatrix_;         /**< data covariance matrix */
    TMatrixD eigenVectors_;             /**< eigenvectors */
    TVectorD eigenValues_;              /**< eigenvalues */
    TMatrix<double> reconstructedSpectra_;/**< rows: wavelength, cols: SEDs */
    TMatrix<double> eigenvalsProjSpec_; /**< doesn't depend on # of eigenvalues kept 
                                             (well length does). Size is (nEigKept,nsed)
                                             Is basically the projection of each SED
                                             onto the first nEigKept eigenvectors */
                                             
};                 

      
}// end namespace sophya

#endif
