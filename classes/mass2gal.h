/**
 * @file  mass2gal.h
 * @brief Simulates galaxy distribution from input over-density distribution
 *
 * @todo Replace the old class for simulating galaxies <CODE>GalFlxTypDist</CODE> 
 *       with the newer ones 
 *
 * @author Alex Abate, Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2009
 * @date 2009
 *
 */
#ifndef MASS2GAL_H_SEEN
#define MASS2GAL_H_SEEN

#include "machdefs.h"
#include "sopnamsp.h"
#include "timing.h"
#include "fabtwriter.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "swfitsdtable.h"

#include "array.h"
#include "stsrand.h"
#include <iostream>
#include <fstream>

#include "gftdist.h"
#include "sinterp.h"
#include "schechter.h"
#include "cosmocalcs.h"
#include "simdata.h"
//#include "selectfunc.h"


/** @class Mass2Gal class 
  *
  * Class to simulate galaxy distribution from SimLSS output
  *	also simulates galaxy properties: absolute magnitude, SED type.
  *
  *	Can also be used to simulate galaxy properties without input
  *	SimLSS output: for simulating galaxy distributions without spacial clustering
  *
  */
class Mass2Gal 
{
public:

    /** Constructor: initialized and automatically removes n reduntant planes left by 
        SimLSS Fourier Transform 
        @param drho         SimLSS over-density grid
        @param su           Cosmology calculator
        @param rg           random generator
        @param nbadplanes   how many redundant planes to removed from @parah drho cube */
    Mass2Gal(TArray<r_8> drho, SimpleUniverse& su, RandomGeneratorInterface& rg,
                                                                int nbadplanes=0); 

    /** Constructor: use when NOT simulating galaxy clustering 
        @param ng           Number of galaxies
        @param su           Cosmology calculator
        @param rg           random generator                                   */
    Mass2Gal(int_8 ng, SimpleUniverse& su, RandomGeneratorInterface& rg);

    /** Constructor for when only wanting to use MaxAbsMag() function 
        @param su           Cosmology calculator
        @param rg           random generator     */
    Mass2Gal(SimpleUniverse& su,RandomGeneratorInterface& rg);
    
    /**  Copy constructor */
    Mass2Gal(Mass2Gal const& a);
    virtual ~Mass2Gal(void);
    virtual Mass2Gal& Set(const Mass2Gal& a);


    /* MAIN FUNCTIONS (these should definitely be public) */
    
    /** Read in cube properties from FITS header: zref_, Dx_, Dkx_, Nx_ etc */
    void ReadHeader(FitsInOutFile& fin);
    
    /** +1 to convert delta->rho/rho^bar and sets -ve mass pixels ->0 */
    sa_size_t CleanNegativeMassCells(); 
    
    /**  Checks there are no cells that have negative mass */
    sa_size_t CheckNegativeMassCells();
    
    /** conv=ngal_per_mpc3*pixel vol; ngal_per_mpc3 calculated by integrating 
    Schechter function of ALL gals */
    void ConvertToMeanNGal(float conv,bool Simple=false); // writes to ngals_ array  
    
    /** Poisson fluctuate ngals_ */
    sa_size_t ApplyPoisson();
    
    /** MAIN FUNCTION: simulates galaxies and puts them into a catalog        
        @param idsim    id number of simulation
        @param fitsname FITS file to write catalog to
        @param gfd      galaxy distributions to simulate with 
        @param extinct  add host galaxy extinction
        @param AMcut    throw out galaxies with too faint absolute magnitudes
        @param SkyArea  area of sky to simulate (in radians) 
        @param ZisRad   simulate cubic sky -> z dimension is exactly the radial
                        direction                                             */
    sa_size_t CreateGalCatalog(int idsim, string fitsname, GalFlxTypDist& gfd, 
        bool extinct=false, bool AMcut=false, double SkyArea=6.3, bool ZisRad=false);
        
    /** Instead of outputing whole simulation, just output a Histo object
        of redshifts to a PPF file 
        @param fitsname FITS file to write catalog to
        @param gfd      galaxy distributions to simulate wit
        @param SkyArea  area of sky to simulate (in radians)                  */
    void CreateNzHisto(string fitsname, GalFlxTypDist& gfd, bool extinct=false, 
                double SkyArea=6.3);
                
    /** Instead of outputing whole simulation, just output all the true redshifts
        to a FITS file 
        @param idsim    id number of simulation
        @param fitsname FITS file to write catalog to
        @param SkyArea  area of sky to simulate (in radians)                  */
    sa_size_t CreateTrueZFile(int idsim, string fitsname, double SkyArea=6.3);
    
    /** Instead of outputing whole simulation can just output redshifts, ra, dec
        to a FITS file 
        @param idsim    id number of simulation
        @param fitsname FITS file to write catalog to
        @param SkyArea  area of sky to simulate (in radians)                  */
    sa_size_t CreateSimpleCatalog(int idsim, string fitsname, double SkyArea=6.3);
        
    /** Calculate maximum observable absolute magnitude as a function of z    */
    void MaxAbsMag(); 
    
    /** Randomise galaxy positions within pixel*/
    void SetRandPos(bool RP){ RandPos_ = RP;
			      cout <<"    Randomise positions set"<<endl; };
			      
	/** To save memory, zero-size the large arrays */
    void ZeroSizeMassArrays(){mass_.ZeroSize(); ngals_.ZeroSize();};
			      
    //----- Added methods Feb 2011 (AA)
    //inline void SetSelectionFunction(SelectionFunctionInterface& sf) { selfuncp_=&sf;}
    // below needed by sim_mcgrids so will prob have to uncomment at some point
    //void ApplySF(SelectionFunctionInterface& sf);	
    void ApplyPZConv(string pzcf);
    void GetCellZBounds(sa_size_t i, sa_size_t j, sa_size_t k,double& zl, double& zc, double& zh);
    void NGalSmArray(TArray<r_8>& ngalssm_array) { ngalssm_array = ngalssm_; }
    void RGalSmArray(TArray<r_8>& rgalssm_array) { rgalssm_array = randgsm_; }
    void SetRandomGrid(int mean_dens);
    
    
    //----- Added methods May 2011 (AA)
    TArray<r_8> ExtractSubArray(double Z, sa_size_t nx, sa_size_t ny, sa_size_t nz);
    TArray<r_8> ExtractSubArray(sa_size_t x1, sa_size_t x2, sa_size_t y1, sa_size_t y2, sa_size_t z1, sa_size_t z2);

    // ---- Added methods Aug 2011 (AA)	
    void ResetSmGrids(){ SFApplied=false; ngalssm_=0; randgsm_=0; };
    
    /** comoving distance to redshift conversion */
    inline double dist2Redshift(double dcom) { return dist2z_(dcom); }
    	
    /** Writes LF to a fits file, for debugging */
    void WriteLFFits(Schechter, double, double, string, int npt=100); 

    /** Return Fourier space pixel size in x-dimension */
    double ReturnDKX(){return Dkx_;}
    /** Return Fourier space pixel size in y-dimension */
    double ReturnDKY(){return Dky_;}
    /** Return Fourier space pixel size in z-dimension */
    double ReturnDKZ(){return Dkz_;}
    /** Return Real space pixel size in x-dimension */
    inline double ReturnDX(){return Dx_;}
    /** Return Real space pixel size in y-dimension */
    inline double ReturnDY(){return Dy_;}
    /** Return Real space pixel size in z-dimension */
    inline double ReturnDZ(){return Dz_;}
    /** Return pixel volume */
    double ReturnPixVol(){return Dx_*Dy_*Dz_;};
    /** Return number of pixels in x-dimension */
    inline long ReturnNX(){return Nx_;}
    /** Return number of pixels in y-dimension */
    inline long ReturnNY(){return Ny_;}
    /** Return number of pixels in z-dimension */
    inline long ReturnNZ(){return Nz_;}
    /** Return total number of pixels */
    inline long ReturnNpix(){return Nx_*Ny_*Nz_;}
    /** Return reference redshift */
    inline double ReturnZref(){return zref_;}
    /** Return comoving distance at reference redshift */
    inline double ReturnDcref(){return DCref_;}
    /** Return entire cube volume */
    double ReturnCubeVol(){return Dx_*Dy_*Dz_*Nx_*Ny_*Nz_;};
    /** Return mass array */
    void MassArray(TArray<r_8>& mass_array) { mass_array = mass_; }
    /** Return galaxy array */
    void NGalArray(TArray<int_8>& ngals_array) { ngals_array = ngals_; }
    /** Returns number of pixels and pixel size in a vector */
    TVector<r_8> ReturnGridSpec();
    /** Return total number of galaxies in simulation */
    int ReturnTotNgals(){return ng_;}// set in ConvertToMeanNGal()
    /** Return maximum observable absolute magnitude as a function of z */
    void ReturnMaxAbsMag(vector<double>& zv, vector<double>& MBmax){zv = zv_; MBmax=MBmax_;};
    /** return 3D index of center pixel */
    void ReturnCenterIndex(vector<int>& idv){vector<int> idv_; idv_.push_back(idmidx_);
										    idv_.push_back(idmidy_);idv_.push_back(idmidz_); idv=idv_;}; 
										    
    /* Functions used only by class (these should probably be private???)     */
    
    //----- Added methods June & July 2010 (RA+AA)
    /** Return grid pixel coordinate */
    void GetCellCoord(sa_size_t i, sa_size_t j, sa_size_t k, double& x, 
                                                        double& y, double& z);
    /** Return grid pixel x-coordinate */                                      
    inline double GetCellX(sa_size_t i)  { return (i-(idmidx_-1))*Dx_; }
    /** Return grid pixel y-coordinate */  
    inline double GetCellY(sa_size_t j)  { return (j-(idmidy_-1))*Dy_; }
    /** Return grid pixel z-coordinate */  
    inline double GetCellZ(sa_size_t k)  { return (k-(idmidz_-1))*Dz_+DCref_; }
    
    /** Convert between Cartesian and spherical coordinates 
        NOTE : the names theta, phi are reversed compared to the usual convention */
    void Conv2SphCoord(double x, double y, double z, double& r, double& phi, double& theta);
        
    /** Draw absolute magnitude in B band MB and gal type: 0-5=early, 6-40=late,
     41-50=starburst*/
    double DrawMagType(GalFlxTypDist&, double& type); 
    
    // ---- Added methods June 2011 (AA)
    void Conv2ParaCoord(double z,double& r){ r=z; };

    //----- Added methods May 2011 (AA)
    double FindPixelAtZ(double Z, sa_size_t nx, sa_size_t ny, sa_size_t nz, sa_size_t& i, sa_size_t& j, sa_size_t& k);
    void GetRange(sa_size_t i, sa_size_t ni, sa_size_t& istart, sa_size_t& iend);

    // SurveyWindow/AddSurveyWindow PROBABLY REDUNDANT
    //void SurveyWindow(double Phi,double dbar);// fills select_ array with 0 
    //where survey does not cover, 1 where it does
    //double AddSurveyWindow(double Phi);// fills mass_ array with -1 where 
    //survey does not cover
    
    
    /* CLASS VARIABLES */
// shouldn't these be protected????

    //SelectionFunctionInterface* selfuncp_;// pointer to selection function
    int_8 ng_;					   /**< Number of galaxies */
    SimpleUniverse& su_;		   /**<  Holds cosmological parameters*/
    RandomGeneratorInterface& rg_; /**<  For random number generation */
    sa_size_t Nx_;  /**< Number of pixels in x direction: NOTE arrays not always defined as (Nx, Ny, Nz)*/
    sa_size_t Ny_;  /**< Number of pixels in y direction: NOTE arrays not always defined as (Nx, Ny, Nz)*/
    sa_size_t Nz_;	/**< Number of pixels in z direction: NOTE arrays not always defined as (Nx, Ny, Nz)*/
    double Dx_;     /**< Pixel size in x-dimension in Mpc */
    double Dy_;     /**< Pixel size in y-dimension in Mpc */
    double Dz_;		/**< Pixel size in z-dimension in Mpc */
    double Dkx_;    /**< Fourier space pixel size in x-dimension in 1/Mpc */
    double Dky_;    /**< Fourier space pixel size in y-dimension in 1/Mpc */
    double Dkz_;    /**< Fourier space pixel size in z-dimension in 1/Mpc */
    double zref_;	/**< Redshift of center pixel */
    double DCref_;  /**< Comoving distance to zref_ */
    int idmidz_;    /**< Indices of centre pixel in z-dimension (assuming 1st pixel is index 1)*/
    int idmidy_;    /**< Indices of centre pixel in y-dimension (assuming 1st pixel is index 1)*/
    int idmidx_;    /**< Indices of centre pixel in x-dimension (assuming 1st pixel is index 1)*/
    SInterp1D dist2z_;  /**< Distance to redshift converter */
    int mean_dens_; /**< mean density of random grid */
    double mean_overdensity_;   /**< mean over-density of simlss grid AFTER setting cells with <-1 to =-1 */

    bool RandPos_;      /**< If true randomise galaxy positions                                       */
    bool fg_nodrho;		/**< A flag to identify if we have drho or random generation                  */
    bool fg_readvals;	/**< A flag to identify if the SimLSS FITS header has been read in            */
    bool fg_cleancells;	/**< A flag to identify if haven't cleaned cells                              */
    bool SFApplied;		/**< Flag turns to true when selection function applied to simulated catalogs */

    /* arrays */
    TArray<r_8> mass_;          /**< 3D array holding rho/rho^bar */
    vector<double> MBmax_,zv_;	/**< Maximum observable absolute magnitude as a function of z */
    TArray<int_8> ngals_;		/**< 3D array holding total galaxy number in each pixel */
    TArray<r_8> ngalssm_,randgsm_;
};

#endif
