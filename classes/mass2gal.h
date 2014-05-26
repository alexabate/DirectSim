/**
 * @file  mass2gal.h
 * @brief Simulates galaxy distribution from input over-density distribution
 *
 * @todo Replace the old class for simulating galaxies <CODE>GalFlxTypDist</CODE> 
 *       with the newer ones 
 * 
 * @todo Add option to output more/less into galaxy catalog
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

//class FunRan;

// DirectSim classes
#include "schechter.h"
#include "gftdist.h"
#include "sinterp.h"
#include "cosmocalcs.h"
#include "simdata.h"
#include "selectfunc.h"
#include "constcosmo.h"

/** @class SimFromOverDensity
  *
  * superclass for Mass2Gal and FieldClusterGals
  */
class SimFromOverDensity
{
public:

    /** Constructor: for simulating WITH clustering information */
    SimFromOverDensity(FitsInOutFile& fin, SimpleUniverse& su, RandomGeneratorInterface& rg)
    : su_(su) , rg_(rg) { 
    
        // first read values from FITS header
        ReadHeader(fin);

        // then read in density cube
        TArray<r_8> drho;
	    fin >> drho;
	    cout << drho.Info();
	    cout << "     Print original drho array size: "<< drho.SizeX() <<"x";
	    cout << drho.SizeY() <<"x"<< drho.SizeZ() <<endl;
	    cout << endl;

        // number of junk planes to remove
        int nbadplanes = drho.SizeX()-Nz_;
        cout << "     Number of planes to remove: "<< nbadplanes << endl;

	    // remove planes
	    mass_ = drho(Range(0, drho.SizeX()-nbadplanes-1), Range::all(), Range::all()).PackElements();
  	    mass_.Show();
        cout << endl;

        // clean -ve mass cells and convert from drho/rho_bar to rho/rho_bar
        sa_size_t nbad = CleanNegativeMassCells();
        cout << "     Number of negative mass cells = "<< nbad << endl;
        cout << endl;

        // set the distance redshift relation
        setDistanceRedshiftRelation();

        fg_nodrho = false;
    };

     /** Constructor: for simulating WITHOUT clustering information */
    SimFromOverDensity(int_8 ng, SimpleUniverse& su, RandomGeneratorInterface& rg)
    : su_(su) , rg_(rg) { 

        // set the distance redshift relation
        setDistanceRedshiftRelation();
        fg_nodrho = true;
    };
    
    void returnDensMeanSigma(double& mean, double& sig) {
	    MeanSigma((mass_-=1.), mean, sig); };

protected:

    /** Read in cube properties from FITS header: zref_, Dx_, Dkx_, Nx_ etc   */
    void ReadHeader(FitsInOutFile& fin);
    
    /** +1 to convert delta->rho/rho^bar and sets -ve mass pixels ->0         */
    sa_size_t CleanNegativeMassCells(); 
    
    /** Checks there are no cells that have negative mass                    */
    sa_size_t CheckNegativeMassCells();

    /** Set distance redshift relation */
    void setDistanceRedshiftRelation() {
	    int_8 nz=1000;
	    vector<double> zrs, codist;
	    double minz=0, maxz=10;
	    double dz = (maxz-minz)/(nz-1);
	    for(int i=0; i<nz; i++) {
	
            double zs = minz + i*dz;
            su_.SetEmissionRedShift(zs);
            double cod = su_.RadialCoordinateMpc(); // radial distance
            
            if (abs(zs-zref_)<0.01)
                cout <<"     Comoving distance to redshift "<< zs <<" is "<< cod << endl; 
            
            zrs.push_back(zs);
            codist.push_back(cod); 
		    }
	    double mind = codist[0];
	    double maxd = codist[codist.size()-1];
	    dist2z_.DefinePoints(codist,zrs,mind,maxd,2*nz);
        };


    /** Convert comoving distance to redshift      
        @param dcom     comoving distance                                     */
    inline double dist2Redshift(double dcom) { return dist2z_(dcom); };

    /** Convert between Cartesian and spherical coordinates 
        NOTE : the names theta, phi are reversed compared to the usual convention 
        @param x        x coordinate
        @param y        y coordinate
        @param z        z coordinate
        @param r        r coordinate
        @param phi      phi coordinate
        @param theta    theta coordiante                                      */
    void Conv2SphCoord(double x, double y, double z, double& r, double& phi, double& theta)
        {   
        // NOTE : the names theta, phi are reversed compared to the usual (ISO 31-11) convention 
        r = sqrt(x*x + y*y + z*z);
	    theta = atan2(y,x);
	    double zoverr = z/r;
	    phi = acos(zoverr);
        };

    /** Return grid pixel comoving coordinate of the grid cell i,j,k
        @param i    index along 1st dimension of grid
        @param j    index along 2nd dimension of grid
        @param k    index along 3rd dimension of grid
        @param x    x coordinate
        @param y    y coordinate
        @param z    z coordinate           
        // i refers to 1ST dim, j refers to 2ND dim, k refers to 3RD dim      */
    void GetCellCoord(sa_size_t i, sa_size_t j, sa_size_t k, double& x, 
                                                        double& y, double& z) {
        x=GetCellX(i);
	    y=GetCellY(j);
	    z=GetCellZ(k);
        return;
        };

    /** Return grid pixel x-coordinate of grid cell i                     
        @param i    index of pixel along x-dimension                          */                         
    inline double GetCellX(sa_size_t i)  { return (i-(idmidx_-1))*Dx_; }
    /** Return grid pixel y-coordinate of grid cell j                    
        @param j    index of pixel along y-dimension                          */  
    inline double GetCellY(sa_size_t j)  { return (j-(idmidy_-1))*Dy_; }
    /** Return grid pixel z-coordinate of grid cell i                     
        @param k    index of pixel along z-dimension                          */  
    inline double GetCellZ(sa_size_t k)  { return (k-(idmidz_-1))*Dz_ + DCref_; };

public:

    /** To save memory, zero-size the large array                             */
    void ZeroSizeMassArray() { mass_.ZeroSize();  };
    
    void GetCellZBounds(sa_size_t i, sa_size_t j, sa_size_t k,
                                        double& zl, double& zc, double& zh) {

	    double x,y,z,dc;
	    GetCellCoord(i,j,k,x,y,z);
	
	    dc = sqrt(x*x+y*y+z*z);// comoving distance to each pixel center
	    zc=dist2Redshift(dc);
	    zh=dist2Redshift(dc+Dz_/2);
	    zl=dist2Redshift(dc-Dz_/2);
        };


    // Basically just 'getter' methods below

    /** Return "cleaned" mass array                                           */
    void MassArray(TArray<r_8>& mass_array) { mass_array = mass_; }
    
    /** Return "cleaned" over-density array                                   */
    void ODensArray(TArray<r_8>& odens_array) { odens_array = (mass_-= 1.); }
    
    /** Return total number of galaxies in simulation                         */
    int ReturnTotNgals() { return ng_; };

    // ---- Return various cube properties
    /** Return Fourier space pixel size in x-dimension                        */
    double ReturnDKX(){return Dkx_;}
    /** Return Fourier space pixel size in y-dimension                        */
    double ReturnDKY(){return Dky_;}
    /** Return Fourier space pixel size in z-dimension                        */
    double ReturnDKZ(){return Dkz_;}
    /** Return Real space pixel size in x-dimension                           */
    inline double ReturnDX(){return Dx_;}
    /** Return Real space pixel size in y-dimension                           */
    inline double ReturnDY(){return Dy_;}
    /** Return Real space pixel size in z-dimension                           */
    inline double ReturnDZ(){return Dz_;}
    /** Return pixel volume                                                   */
    double ReturnPixVol(){return Dx_*Dy_*Dz_;};
    /** Return number of pixels in x-dimension                                */
    inline long ReturnNX(){return Nx_;}
    /** Return number of pixels in y-dimension                                */
    inline long ReturnNY(){return Ny_;}
    /** Return number of pixels in z-dimension                                */
    inline long ReturnNZ(){return Nz_;}
    /** Return total number of pixels                                         */
    inline long ReturnNpix(){return Nx_*Ny_*Nz_;}
    /** Return reference redshift                                             */
    inline double ReturnZref(){return zref_;}
    /** Return comoving distance at reference redshift                        */
    inline double ReturnDcref(){return DCref_;}
    /** Return entire cube volume                                             */
    double ReturnCubeVol(){return Dx_*Dy_*Dz_*Nx_*Ny_*Nz_;};
    /** Returns number of pixels and pixel size in a vector                   */
    TVector<r_8> ReturnGridSpec();
    /** return 3D index of center pixel */
    void ReturnCenterIndex(vector<int>& idv){ vector<int> idv_; idv_.push_back(idmidx_);
				    idv_.push_back(idmidy_);idv_.push_back(idmidz_); idv=idv_;}; 

protected:
    SimpleUniverse& su_;		        /**<  Holds cosmological parameters       */
    RandomGeneratorInterface& rg_;  /**<  For random number generation        */
    sa_size_t Nx_;  /**< Number of pixels in x direction: NOTE arrays not always defined as (Nx, Ny, Nz)*/
    sa_size_t Ny_;  /**< Number of pixels in y direction: NOTE arrays not always defined as (Nx, Ny, Nz)*/
    sa_size_t Nz_;	/**< Number of pixels in z direction: NOTE arrays not always defined as (Nx, Ny, Nz)*/
    double Dx_;     /**< Pixel size in x-dimension in Mpc                     */
    double Dy_;     /**< Pixel size in y-dimension in Mpc                     */
    double Dz_;		/**< Pixel size in z-dimension in Mpc                     */
    double Dkx_;    /**< Fourier space pixel size in x-dimension in 1/Mpc     */
    double Dky_;    /**< Fourier space pixel size in y-dimension in 1/Mpc     */
    double Dkz_;    /**< Fourier space pixel size in z-dimension in 1/Mpc     */
    double zref_;	/**< Redshift of center pixel                             */
    double DCref_;  /**< Comoving distance to zref_                           */
    int idmidz_;    /**< Indices of centre pixel in z-dimension (assuming 1st pixel is index 1)*/
    int idmidy_;    /**< Indices of centre pixel in y-dimension (assuming 1st pixel is index 1)*/
    int idmidx_;    /**< Indices of centre pixel in x-dimension (assuming 1st pixel is index 1)*/
    TArray<r_8> mass_;          /**< 3D array holding rho/rho^bar                */
    double mean_overdensity_;   /**< mean over-density of simlss grid AFTER setting cells with <-1 to =-1 */
    SInterp1D dist2z_;          /**< Distance to redshift converter              */
    bool fg_nodrho; /**< A flag to identify if we have drho or random generation */
    int_8 ng_;				    /**< Number of galaxies                          */
};

/** @class Mass2Gal class 
  *
  * Class to simulate galaxy distribution from SimLSS output
  *	also simulates galaxy properties: absolute magnitude, SED type.
  *
  *	Can also be used to simulate galaxy properties without input
  *	SimLSS output: for simulating galaxy distributions without spacial clustering
  *
  * @todo: have taken out ApplyPZConv, ApplySF for now - need to look back at this
  *
  */
class Mass2Gal : public SimFromOverDensity
{
public:

    /** Constructor: initialized and automatically removes n reduntant planes left by 
        SimLSS Fourier Transform 
        @param fin          FITS file containing over density grid
        @param su           Cosmology calculator
        @param rg           random generator                                  */
    Mass2Gal(FitsInOutFile& fin, SimpleUniverse& su, RandomGeneratorInterface& rg)
    : SimFromOverDensity(fin, su, rg) {
	    SFApplied = false;
	    setSimProperties(2.*PI);
        }; 

    /** Constructor: use when NOT simulating galaxy clustering 
        @param ng           Number of galaxies
        @param su           Cosmology calculator
        @param rg           random generator                                  */
    Mass2Gal(int_8 ng, SimpleUniverse& su, RandomGeneratorInterface& rg)
    : SimFromOverDensity(ng, su, rg) {
        SFApplied = false;
        setSimProperties(2.*PI);
        };

    /** Constructor for when only wanting to use MaxAbsMag() function 
        @param su           Cosmology calculator
        @param rg           random generator                                  */
    Mass2Gal(SimpleUniverse& su, RandomGeneratorInterface& rg)
    : SimFromOverDensity(1, su, rg) {
	    SFApplied = false; 
        setSimProperties(2.*PI);
	    };
    
    /**  Copy constructor 
    Mass2Gal(Mass2Gal const& a)
    :  su_(a.su_) , rg_(a.rg_) {
	    cout <<"    Mass2Gal COPY constructor"<<endl;
	    Set(a); };*/

    /** Destructor */
    virtual ~Mass2Gal(void) { };


    /** Copy class variables in Mass2Gal                                      */
    //virtual Mass2Gal& Set(const Mass2Gal& a);


    /* MAIN FUNCTIONS (these should definitely be public)        
        @param SkyArea        area of sky to simulate (in radians)
        @param doAbsMagCut    throw out galaxies with too faint absolute magnitudes    
        @param isConstDens    simulation has constant number density      
        @param doRandPos      randomize galaxy positions within grid cell
        @param isZRadDir      simulate cubic sky -> z dimension is exactly the 
                              radial direction                                             
        @param doXtinct       add host galaxy extinction                      */
    void setSimProperties(double SkyArea, bool doAbsMagCut=false, bool isConstDens=false,
              bool doRandPos=true, bool isZRadDir=false, bool doXtinct=false) {
              
        SkyArea_ = SkyArea;
        doAbsMagCut_ = doAbsMagCut;
        isConstDens_ = isConstDens;
        doRandPos_ = doRandPos;
        isZRadDir_ = isZRadDir;
        doXtinct_ = doXtinct;
        
        };
    
    
    /** conv=ngal_per_mpc3*pixel vol; ngal_per_mpc3 calculated by integrating 
    Schechter function of ALL gals                                            */
    //void ConvertToMeanNGal(float conv, bool Simple=false); // writes to ngals_ array  
    
    /** Poisson fluctuate ngals_                                              */
    //sa_size_t ApplyPoisson();
    
    /** MAIN FUNCTION: simulates galaxies and puts them into a catalog        
        @param idsim         id number of simulation
        @param outfileroot   data written out to files beginning [outfileroot]
        @param conv          how to convert mass to n galaxies
        @param gfd           galaxy distributions to simulate with 
        @param doSimple      only output galaxy positions                        */
    sa_size_t CreateGalCatalog(int idsim, string outfileroot, double conv, 
                                       GalFlxTypDist& gfd, bool doSimple=false);
        
        
    /** Instead of outputing whole simulation, just output a Histo object
        of redshifts to a PPF file 
        @param fitsname FITS file to write catalog to
        @param gfd      galaxy distributions to simulate wit
        @param SkyArea  area of sky to simulate (in radians)                  */
    //void CreateNzHisto(string fitsname, GalFlxTypDist& gfd);
                
                
    /** Instead of outputing whole simulation, just output all the true redshifts
        to a FITS file 
        @param idsim    id number of simulation
        @param fitsname FITS file to write catalog to
        @param SkyArea  area of sky to simulate (in radians)                  */
    //sa_size_t CreateTrueZFile(int idsim, string fitsname);
    
    
    /** Instead of outputing whole simulation can just output redshifts, ra, dec
        to a FITS file 
        @param idsim    id number of simulation
        @param fitsname FITS file to write catalog to
        @param SkyArea  area of sky to simulate (in radians)                  */
    //sa_size_t CreateSimpleCatalog(int idsim, string fitsname);
        
        
    /** Calculate faintest observable absolute magnitude as a function of z   */
    SInterp1D MaxAbsMag(); 
   
			      
	
			      
    //----- Added methods Feb 2011 (AA)
    //inline void SetSelectionFunction(SelectionFunctionInterface& sf) { selfuncp_=&sf;}
    //void ApplySF(SelectionFunctionInterface& sf);	

    /** Apply photo-z convolution */
    //void ApplyPZConv(string pzcf);

    /** Return the redshift bounds of the grid cell i,j,k
        @param i    index along 1st dimension of grid
        @param j    index along 2nd dimension of grid
        @param k    index along 3rd dimension of grid
        @param zl   lowest redshift in grid cell
        @param zc   central redshift of grid cell
        @param zh   highest redshift in grid cell                             */
    //void GetCellZBounds(sa_size_t i, sa_size_t j, sa_size_t k,double& zl, double& zc, double& zh);

    /** Return grid of galaxies after photo-z smearing applied                */
    //void NGalSmArray(TArray<r_8>& ngalssm_array) { ngalssm_array = ngalssm_; }

    /** Return random grid after photo-z smearing applied                     */
    //void RGalSmArray(TArray<r_8>& rgalssm_array) { rgalssm_array = randgsm_; }

    /** Set the mean density of the random grid 
        @param mean_dens    mean density in the grid cells of the random grid */
    void SetRandomGrid(int mean_dens) {
        mean_dens_ = mean_dens;
        };
    
    //----- Added methods May 2011 (AA)

    /** Extract a sub grid from the full grid centered on redshift @param Z and 
        with size nx*ny*nz 
        @param Z    redshift sub-grid centered on
        @param nx   size of sub-grid along 1st dimension
        @param ny   size of sub-grid along 2nd dimension
        @param nz   size of sub-grid along 3rd dimension                      */
    TArray<r_8> ExtractSubArray(double Z, sa_size_t nx, sa_size_t ny, sa_size_t nz);

    /** Extract a sub grid from the full grid from pixel x1 to x2 along 1st 
        dimension, y1 to y2 along 2nd dimension, z1 to z2 along 3rd dimension
        @param x1   start of sub-grid along 1st dimension
        @param x2   end of sub-grid along 1st dimension
        @param y1   start of sub-grid along 2nd dimension
        @param y2   end of sub-grid along 2nd dimension                     
        @param z1   start of sub-grid along 3rd dimension
        @param z2   end of sub-grid along 3rd dimension                       */
    TArray<r_8> ExtractSubArray(sa_size_t x1, sa_size_t x2, sa_size_t y1, sa_size_t y2, 
                                sa_size_t z1, sa_size_t z2);

    // ---- Added methods Aug 2011 (AA)	
    /** Reset the photo-z smeared grids                                       */
    //void ResetSmGrids(){ SFApplied=false; ngalssm_=0; randgsm_=0; };
    
    
    	
    /** Writes LF to a fits file, for debugging                               */
    //void WriteLFFits(Schechter, double, double, string, int npt=100); 

    
    // ---- Return grids
    
    
    /** Return galaxy array                                                   */
    //void NGalArray(TArray<int_8>& ngals_array) { ngals_array = ngals_; }
    
    
    
    /** Return maximum observable absolute magnitude as a function of z */
    void ReturnMaxAbsMag(vector<double>& zv, vector<double>& MBmax)
     { zv = zv_; MBmax=MBmax_; };
    
										    
private:
        
    /** Draw absolute magnitude in B band MB and gal type: 0-5=early, 6-40=late,
        41-50=starburst*/
    double DrawMagType(GalFlxTypDist&, double& type); 
    
    // ---- Added methods June 2011 (AA)
    /** convert radial distance to be parallel to z dimension                 */
    void Conv2ParaCoord(double z,double& r){ r=z; };

    //----- Added methods May 2011 (AA)
    double FindPixelAtZ(double Z, sa_size_t nx, sa_size_t ny, sa_size_t nz, 
                        sa_size_t& i, sa_size_t& j, sa_size_t& k);
    
    void GetRange(sa_size_t i, sa_size_t ni, sa_size_t& istart, sa_size_t& iend);

    // SurveyWindow/AddSurveyWindow PROBABLY REDUNDANT
    //void SurveyWindow(double Phi,double dbar);// fills select_ array with 0 
    //where survey does not cover, 1 where it does
    //double AddSurveyWindow(double Phi);// fills mass_ array with -1 where 
    //survey does not cover
    
    
  
protected:
    bool SFApplied;		/**< true if selection function applied to simulation */
    double SkyArea_;    /**< area of sky to simulate (must be < cube area)    */
    bool doAbsMagCut_;  /**< if true don't simulate gals below some abs mag   */
    bool isConstDens_;  /**< if true simulation has constant number desity    */
    bool doRandPos_;    /**< if true randomise galaxy positions               */
    bool isZRadDir_;    /**< if true z-dimension of cube is exactly radial dir*/
    bool doXtinct_;     /**< if true simulate host galaxy extinction          */
    int mean_dens_;     /**< mean density of random grid                      */
    vector<double> MBmax_; /**< Maximum observable abs mag as a function of z */
    vector<double> zv_;	   /**< z values MBmax_ defined at                    */
    
    
    SelectionFunctionInterface* selfuncp_;/**< pointer to selection function */  
    //TArray<int_8> ngals_;		/**< 3D array holding total galaxy number in each pixel           */
    //TArray<r_8> ngalssm_;       /**< array of n galaxies per cell after applying photo-z smearing */
    //TArray<r_8> randgsm_;       /**< array of random grid after applying photo-z smearing         */
    

};


/** @class FieldClusterGals class 
  *
  * Class to simulate galaxy distribution from SimLSS output
  *	where galaxies are assigned to be in a "cluster" or in the "field"
  *
  */
class FieldClusterGals : public SimFromOverDensity
{
public:

    FieldClusterGals(FitsInOutFile& fin,  SimpleUniverse& su, RandomGeneratorInterface& rg) 
    : SimFromOverDensity(fin, su, rg)
    { };

    void simulateGalaxies(double conv, string outfile);

protected:

};

//} // End namespace SOPHYA
#endif
