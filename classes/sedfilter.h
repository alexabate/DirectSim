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
 * AA: made improvements to class inheritance structure March 2015
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
#include "mydefrg.h"

// ROOT libraries
// removed these root libraries that are not needed
//#include "TMinuit.h"
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
  *
  * SED Base class
  * 
  */
class SpecEnergyDist : public ClassFunc1D
{
public:
    SpecEnergyDist() { };
};

/** @class 
  *
  * Base Filter class 
  * To load in Filter transmission functions:
  * two columns: observed wavelength (m) ; transmission probability (between 0 and 1)
  *
  * Basically just a wrapper for SInterp1D; adds level where it finds where filters are stored before 
  * trying to read file in
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
        @param npt          number of points to do interpolations with between lmin and lmax        
        @param nComments    number of comment lines the filter file has at start
        @param zero_outside value of interpolation outside lmin,lmax is zero if true */
    Filter(string& filename, double lmin, double lmax, int npt = 1024, int nComments=0, bool zero_outside=true)
    { readFilter(filename, lmin, lmax, npt, nComments, zero_outside); };
                                      
   /** Read in filter from file 
       @param filename     file to read transmission function from
       @param lmin         min observed wavelength value
       @param lmax         max observed wavelength value
       @param npt          number of points to do interpolations with between lmin and lmax   
       @param nComments    number of comment lines the filter file has at start
       @param zero_outside value of interpolation outside lmin,lmax is zero if true */
    void readFilter(string& filename, double lmin, double lmax, int npt=1024, int nComments=0, bool zero_outside=true);

};


/** @class
  *
  * Spectral Energy Distribution (SED) object. It returns flux of a single SED (in wavelength units) at
  * a given wavelength (in meters). 
  * 
  * The original SED is read from a text file that has two columns: wavelength (in meters) and
  * flux (in wavelength units).  The flux value at the given wavelength is returned by the class.
  *
  * The original SED can be optionally interpolated (with a second SED read into a class variable) or reddened.
  * Interpolation *must* be done before reddening.
  * 
  */
class SED : public SpecEnergyDist
{
public:

    // code copied in both SimData and here (should find a better way)
    //enum dustLaw{NoDust=0, Card=1, Calz=2};

    /** Default constructor */
    SED() : sed2_(sed2init_)
    { isRead_=false; isRedden_=false; isInterp_=false; _test_=1.;
      a_=1.; EBmV_=0.; RvCard_=3.5; law_=0.; };

  
  
    /** Constructor
        To read in wavelength,flux values from a file, and return flux values based on the linear 
        interpolation of numbers in that file. For the most accurate interpolation results, the interpolation
        lookup table (set by lmin, lmax, npt) must match the wavlength grid read from the file as closely as 
        possible. @warning allows wavelengths to be negative
        @param filename    file containing SED data: first column wavelength in m, second column flux_lambda
        @param lmin        minimum wavelength of SED lookup table
        @param lmax        maximum wavelength of SED lookup table
        @param npt         number of points in SED lookup table                                             */
    SED(string& filename, double lmin, double lmax, int npt=1024);
    
    
    /** Copy constructor                                                                                    */
    SED(const SED& sed);
    virtual SED& Set(const SED& a);
      
      
    /** Read in SED from file 
        @param fname    name of file to read from                                                           */
    void readSED(string& fname, double xmin, double xmax, int npt=1024);
      
      
    /** To return the relevant flux value. What is returned depends on isRedden and isInterp values
        Either: interpSED(), addReddening() or interpAddReddening() is called 
        @param lambda   rest-frame wavelength (meters)                                                      */
    virtual double operator()(double lambda) const { return returnFlux(lambda); };
    
    
    /** To return the relevant flux value. What is returned depends on isRedden and isInterp values
        Either: interpSED(), addReddening() or interpAddReddening() is called 
        @param lambda   rest-frame wavelength (meters)                                                      */
    double returnFlux(double lambda) const;
    
            
    /** Turn on interpolation of the SED with second SED (copies second SED into class variable)
        @param sed2     second SED to interpolate with (class variable now becomes defined)
        @param a        fraction of interpolation is from the original SED (the one read from a file). The
                        rest is from the second SED.                                                        */
    void doInterp(SED* sed2, double a=1.) {
        // Must take as argument a POINTER to an SED object to copy sed2 into sed2_
        sed2_ = sed2;
        a_ = a; 
        isInterp_ = true;
        };
    
    /** Turn on reddening of the SED
        @param EBmV    amount of reddening, color excess E(B-V) parameter
        @param law     Cardelli (law=0) or Calzetti (law=1) reddening law
        @param RvCard  parameter relevant to Cardelli law                                                   */
    void doRedden(double EBmV=0., int law=0, double RvCard=3.5) {
        EBmV_ = EBmV;
        law_ = law;
        RvCard_ = RvCard;
        isRedden_ = true;
        };
    
    /** Turn off interpolation */
    void stopInterp() { isInterp_ = false; a_ = 1.; };
    
    /** Turn off reddening */
    void stopRedden() { isRedden_ = false; };
    
    /** Return bool if SED was interpolated                                                                 */
    bool returnInterpBool() { return isInterp_; };
    
    /** Return bool if SED was reddened                                                                     */
    bool returnReddenBool() { return isRedden_; };
    
    /** Return bool if SED has been read from a file                                                        */
    bool returnReadBool() { return isRead_; };
    
    /** Return fraction of original SED that is in output (if interpolated will be <1)                      */ 
    double returnInterpFrac() { return a_; };
    
    /** Return reddening parameters. Returns vector of [E(B-V), law_id, RvCard*] (*if applicable). If no 
        reddening was applied E(B-V)=0 and law_id=-1                                                      */ 
    vector<double> returnReddeningPars() { 
        vector<double> redPars;
        redPars.push_back(EBmV_);
        if (!isRedden_) {
            cout <<"     No Reddening was applied! "<<endl;
            redPars.push_back(-1.);
            }
        else {
            redPars.push_back(law_);
            if (law_>0)
                cout <<"     Calzetti law applied"<<endl;
            else {
                cout <<"     Cardelli law applied"<<endl;
                redPars.push_back(RvCard_);
                }
            }
    return redPars; };
    
    /** Checks which constructor was called (1=default, 2=main, 3=copy)                                     */
    double returntest(){ return _test_; };

    
protected:
    // THESE SHOULD BE PROTECTED METHODS
    /** To return the reddened value 
        @param lambda   rest-frame wavelength                                                               */
    double addReddening(double lambda) const;
    
    /** To return the interpolated value
        @param lambda   rest-frame wavelength                                                               */
    double interpSED(double lambda) const;
    
    /** To return the reddened and interpolated value
        @param lambda   rest-frame wavelength                                                               */
    double interpAddReddening(double lambda) const;
    
    
protected:
    SInterp1D sed_;              /**< The original SED read in and to be used                               */
    bool isRead_;                /**< true if original SED file has been read                               */
    bool isRedden_;              /**< true if SED has been reddened                                         */
    bool isInterp_;              /**< true if SED has been interpolated                                     */
    SED* sed2_;                  /**< 2nd SED class to interpolate between with                             */
    SED* sed2init_;              /**< a default setting of sed2                                             */
    double a_;                   /**< fraction of original SED in final output                              */
    double EBmV_;                /**< E(B-V) reddening value                                                */
    double RvCard_;              /**< Cardelli law parameter                                                */
    int law_;                    /**< Reddening law used 0=Cardelli, 1=Calzetti                             */
    double _test_;               /**< variable for testing which constructor was called                     */
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
        @param z    redshift of galaxy SED                                                                  */
    SEDzFilterProd(SpecEnergyDist& f, Filter& g, double z=0.)
    : sed_(f) , filt_(g) , z_(z) {  };
    
    
    /** Returns redshifted SED multiplied by the filter transmission and wavelength,
        \f$ SED(\lambda/(1+z))*Filter(\lambda)*\lambda \f$ where \lambda is the 
        observed wavelength    
        @param lambda   observed wavelength (in meters)                                                     */
    virtual double operator()(double lambda) const {
        double lambdaE = lambda/(1.+z_);
        return (sed_(lambdaE)*filt_(lambda)*lambda);
        };
        
// *lambda if SED in units of wavelength and magnitudes have dnu/nu definition
// no lambda factor if SED in units of wavelength and magnitudes have dnu definition
// /lambda if SED in units of frequency and magnitudes have dnu/nu definition
// /lambda^2 if SED in units of frequency and magnitudes have dnu definition

protected:
  SpecEnergyDist& sed_;       /**< SED                                        */
  Filter& filt_;              /**< filter function                            */
  double z_;                  /**< redshift of galaxy SED                     */
};


/** @class SEDIGM class
  * 
  * Add IGM absorption
  */
class SEDIGM : public SpecEnergyDist // inherit from SED or from SpecEnergyDist?
{
public:

    /** Constructor if using Madau law for IGM
        @param S       SED to add IGM absorption to
        @param zSED    redshift of galaxy SED                                 
        @param isLyC   include Lyman continuum absorption                                                   */
    SEDIGM(SED& S, double zSED , bool isLyC=true, double lmin=1e-8, double lmax=2.5e-6, int nl=500)
    : sed_(S) , igm_(igm_default_) , zSED_(zSED) , isLyC_(isLyC) { 
        
        isMadau_ = true; 
        setupMadau(lmin, lmax, nl);
        };

    /** Constructor if applying IGM from a transmission function
        @param S       SED to add IGM absorption to
        @param T       IGM transmission function
        @param zSED    redshift of galaxy SED                                 */
    SEDIGM(SED& S, IGMTransmission& T, double zSED)
    : sed_(S) , igm_(T) , zSED_(zSED) , isLyC_(false) { isMadau_ = false; };
    
    /** Return SED with IGM absorption applied CHECK THIS METHOD WORKS!
        @param lambda   restframe wavelength                                  */
    virtual double operator()(double lambda) const {
    
        double sedWithIGM;
        if (isMadau_) {
            //Madau madau(5, isLyC_); // 5 refers to max number of Lyman lines
            // THINK THIS SHOULD BE RETURNING THE REST FRAME TRANSMISSION
            //double trans = madau.returnRestFrameTransmission(lambda, zSED_);
            
            sedWithIGM = sed_(lambda)*trans_(lambda);
            }
        else {
            double lambdaObs = lambda*(1 + zSED_); // this is because the transmission files are currently
                                                   // a function of OBS FRAME wavelength
                                                   // this undoes that so e.g. Lyman alpha effect at 1216A
                                                   // for a gal at redshift 2 will be actually found at
                                                   // 1216*(1+2) = 3648A
            double trans = igm_(lambdaObs);
            sedWithIGM = sed_(lambda)*trans;
            }
        
        return sedWithIGM;
        };
     
     
protected:   
    void setupMadau(double lmin, double lmax, int nl) {
        
        Madau madau(5, isLyC_); // 5 refers to max number of Lyman lines to calculate
                                // with Madau implementation, 5 is the maximum possible
        
        vector<double> rflam;
        vector<double> trans;
        double dl = (lmax-lmin)/(nl-1.);
        for (int i=0; i<nl; i++) {
            double lam = lmin + i*dl;
            double t = madau.returnRestFrameTransmission(lam, zSED_);
            rflam.push_back(lam);
            trans.push_back(t);
            }
        trans_.DefinePoints(rflam, trans, rflam[0], rflam[rflam.size()-1], 2*nl);
        };


protected:
    SED& sed_;                    /**< SED to add IGM absorption to                                         */
    IGMTransmission& igm_;        /**< IGM absorption along line of sight                                   */
    IGMTransmission igm_default_; /**< default class to initialise igm_                                     */
    double zSED_;                 /**< redshift of SED                                                      */
    bool isLyC_;                  /**< include Lyman continuum absorption (only can be true if Madau)       */
    bool isMadau_;                /**< true if Madau absorption being used                                  */
    SInterp1D trans_;             /**< Madau transmission look up table                                     */
};


/** ReadSedList class
  * 
  * Class to read in SEDs from a list of files
  * Can optionally interpolate new SEDs between the SEDs in the list and redden
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
    
    
    /** Main program, reads SEDs from file and arranges them in a vector of pointers pointing to each SED 
        object 
        @param lmin    minimum wavelength of SED in meters
        @param lmax    maximum wavelength of SED in meters 
        @param npt     number of interpolation points for SED func                                          */
    void readSeds(double lmin=5e-8, double lmax=2.5e-6, int npt=10000);
    
    
    /** Interpolate nInterp SEDs between each SED in list of files. If interpolating between the SEDs call 
        this method straight after readSeds() (and before reddenSeds()).
        @param nInterp number of SEDs to interpolate between each SED in file list                          */
    void interpSeds(int nInterp);
    
    
    /** Reorder SEDs so that interpolated SEDs lie between the templates they were interpolated from. This
        should be called after interpSeds() and before reddenSeds() (if reddening)                          */
    void reorderSEDs();
    
    
    /** If reddening all the SEDs call this method straight after #readSeds(), 
        and after #interpSeds() if interpolating.
        Currently only applies Cardelli law, and can't limit max reddening for
        elliptical galaxies 
        @param nStepRed number of times to redden 
        @param redMax   maximum limit to redden to (even steps between 0 and redMax */
    //void reddenSeds(int nStepRed, double redMax);
    
    /** If reddening all the SEDs call this method straight after #readSeds(), and after #interpSeds() if 
        interpolating. Assumes all ellipical SEDs are at the start of the original file list read in. Returns 
        a vector of reddening values applied to the SEDs (Ellipticals then others).
        @param nPerSED        number of reddened SEDs to make per original SED
        @param method         method for distributing reddening values: 0=uniform, 1=exponential
        @param maxidEl        maximum index of an Elliptical galaxy in the sedArray_ (-> all El must be at start)
        @param redMaxEl       maximum E(B-V) value an Elliptical galaxy is allowed 
        @param redMaxOther    maximum E(B-V) value other galaxies are allowed                               */
    vector<double> reddenSeds(int nPerSED, int method=0, int maxidEl=0, double redMaxEl=0.1, 
                                                                                      double redMaxOther=0.3);
    
    /** Write contents of sedArray to a file */
    void writeSpectra(string outFile, double lmin=5e-8, double lmax=2.5e-6, int nl=1500);
                                                                 
    /** Return sedArray */
    vector<SED*> getSedArray() { return sedArray_; };
    
    /** Return total number of SEDs (could be >= to nsed depending if interpolation
        was done or reddening was applied) */
    int getNTot() { return ntot_; };
    
                
    /** Return number of SEDs read in from the initial file */
    int getNSed() { return nsed_; };
       
    
    /** Return a vector of the SED filenames */
    vector<string> returnSedFilenames() { return sedFiles_; };
    
    
    /** Return a vector of the SED filenames without their extension (leaves dot though)                    */
    vector<string> returnSedFilenamesNoExt() { 
        vector<string> noext;
        for (int i=0; i<sedFiles_.size(); i++) {
            string fne = sedFiles_[i];
            string fn = fne.substr(0, fne.size()-3);
            noext.push_back(fn); 
            }
        return noext; };
    
    
    /** Return SED reddening values applied */
    vector<double> returnReddening() { return reds_; };
    
protected:

    /** Counts SEDs in file */
    void countSeds(string sedFileFullPath);
                   
private:
    int prt_;                   /**<  print level                             */
    string sedDir_;             /**<  path of SED list file                   */
    string sedFileFullPath_;    /**<  full path and filename to SED list      */
    int nsed_;                  /**<  number of SEDs in the file sedFile      */
    bool isRead_;               /**<  if SEDs have been read don't read again */
    int ntot_;                  /**<  total number of SEDs after interpolation /
                                      and/or reddening                        */
    vector<SED*> sedArray_;     /**<  pointer to ntot_ SED objects            */
    vector<string> sedFiles_;   /**<  vector of each SED file name (no path)  */
    bool isInterp_;             /**<  if SEDs are added to by interpolation   */
    bool isRedden_;             /**<  if SEDs are added to by reddening       */
    vector<double> reds_;       /**<  E(B-V) reddening value for each SED     */
};


/*******************************************************************************
*                                                                              *
*                              FILTER CLASSES                                  *
*                                                                              *
*******************************************************************************/




/** @class BlueShiftFilter
  *
  * Not entirely sure if this class is needed. It is used in returning the restframe flux of an SED at
  * some redshift in some filter
  */
class BlueShiftFilter : public Filter
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
    
    /** returns the value of the filter transmission at #lambda, divided by #lambda                         */
    virtual double operator()(double lambda) const {
        return (filt_(lambda)/lambda);//********* DOUBLE CHECK THIS *********//
        }  

protected:
  Filter& filt_;    /**< class holding the filter function                                                  */
};


/** @class FilterXLambda class
  * To multiply filter transmission value by the wavelength in order to calculate filter effective wavelength
  *
  */
class FilterXLambda : public ClassFunc1D
{
public:
    
    /** Constructor */
    FilterXLambda(Filter& g)
    : filt_(g){ }
    
    /** returns the value of the filter transmission at #lambda, multiplied by 
        #lambda */
    virtual double operator()(double lambda) const {
        return (filt_(lambda)*lambda);//********* DOUBLE CHECK THIS *********//
        }  

protected:
  Filter& filt_;          /**< class holding the filter function                                            */
};



//--- Does int F(lambda)/lambda dlambda or SED(lambda)*F(lambda) or whatever/any ClassFunc1D
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
        @param lmin    minimum wavelength of filter in meters
        @param lmax    maximum wavelength of filter in meters 
        @param nl      number of points between lmin and lmax        */
    void readFilters(double lmin=5e-8, double lmax=2.5e-6, int nl=500);
    
    /** Write contents of filterArray to a file */
    void writeFilters(string outFile, double lmin=5e-8, double lmax=2.5e-6, int nl=1500);
                                                                 
    /** Return filterArray */
    vector<Filter*> getFilterArray() { return filterArray_; };
    
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
