/**
 * @file  cat2grid.h
 * @brief grid a galaxy catalog for power spectrum analysis
 *
 * @todo move some methods out to a selection function class eg SaveSelecFunc
 *
 * @author Alex Abate, Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2009
 * @date 2009
 *
 */
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <typeinfo>

// sophya
#include "machdefs.h"
#include "sopnamsp.h"
#include "array.h"
#include "fabtwriter.h"
#include "hisprof.h"
#include "histerr.h"
#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "stsrand.h"
#include "mydefrg.h"
#include "swfitsdtable.h"
#include "fitsmanager.h"
#include "resusage.h"
#include "timing.h"

#include "cosmocalcs.h"
#include "constcosmo.h"
#include "sinterp.h"
#include "schechter.h"
#include "selectfunc.h"
#include "mass2gal.h"

/** @class GalRecord

	GalRecord class holds galaxy variables

*/
class GalRecord {
public:

    /** Constructor */
    GalRecord() { alpha=delta=glong=glat=xcoo=ycoo=zcoo=0.; zs=zo=0.; type=1; u=g=r=i=z=y=24.; }

    /** Print galaxy properties */
	void Print()
		{ cout <<"Printing GalRecord: ";
		  cout <<"alpha="<< alpha <<", delta="<< delta <<", zs="<< zs <<", zo="<< zo;}

	double alpha;    /**< right ascension coodinate                           */
	double delta;    /**< declination coodinate                               */
	double glong;    /**< galactic longitude coordinate                       */
	double glat;     /**< galactic latitude coordinate                        */
	double xcoo;     /**< Euclidean x coordinate                              */
	double ycoo;     /**< Euclidean y coordinate                              */
	double zcoo;     /**< Euclidean z coordinate                              */
	double zs;       /**< spectroscopic redshift                              */
	double zo;       /**< observed redshift                                   */
	int type;        /**< galaxy type                                         */
	double u;        /**< u magnitude                                         */
	double g;        /**< g magnitude                                         */
	double r;        /**< r magnitude                                         */
	double i;        /**< i magnitude                                         */
	double z;        /**< z magnitude                                         */
	double y;        /**< y magnitude                                         */
};


/** @class Cat2Grid
  *
  * Reads in a galaxy catalog (ra,dec,z) from a SwFitsDataTable which must have column names:
  * "phi" and "theta", but the z column name can be set via an argument
  *
  * Lays a grid across the galaxy catalog with grid cell size (default 8 Mpc) and 
  * grid position specified in SetGrid function 
  * 
  * OPTIONALLY: adds statistical photo-z error sig = PZerr*(1+z) to z-dimension (if PZerr>0)
  * OPTIONALLY: corrects for selection effects (if SetSelectionFunction is called, which sets sfcompute_=true)
  */
class Cat2Grid
{
public:

	/** Constructor - Normal usage
	    @param dt       data table containing galaxy catalog
	    @param su       cosmology calculator
	    @param rg       random number generator
	    @param fos      FITS file containing gridded galaxy data
	    @param ZOCol    name of column in galaxy catalog containing observed z
	    @param ZSCol    name of column in galaxy catalog containing spec z
	    @param RadialZ  if true sets the z-dimension to be radial direction 
	    @param PZErr    size of Gaussian photo-z error 
	    @param Print    if true prints extra to the screen                    */
	Cat2Grid(SwFitsDataTable& dt, SimpleUniverse& su, RandomGenerator& rg,
	         FitsInOutFile& fos, string ZOCol="zp", string ZSCol="z",
	         bool RadialZ=false, double PZErr=0,bool Print=true);
	         
	/** Constructor - Takes the interp z->d functions as arguments rather than
	    calculating them every time so you can use the Row2Record, 
	    Rec2EuclidCoord functions
	    @param dt       data table containing galaxy catalog
	    @param su       cosmology calculator
	    @param rg       random number generator
	    @param dist2z   distance to redshift conversion table
	    @param z2dist   redshift to distance conversion table
	    @param ZOCol    name of column in data table containing observed z
	    @param ZSCol    name of column in data table containgin spec z
	    @param RadialZ  if true sets the z-dimension to be radial direction   */
	Cat2Grid(SwFitsDataTable& dt, SimpleUniverse& su, RandomGenerator& rg,
	         SInterp1D dist2z, SInterp1D z2dist, 
	         string ZOCol="zp", string ZSCol="zs", bool RadialZ=false);

    /** Copy constructor */
	Cat2Grid(Cat2Grid const& a);
	
	/** Destructor */
	virtual ~Cat2Grid() {};
	
	/** Copy method */
	virtual Cat2Grid& Set(const Cat2Grid& a);
	
	
	/* MAIN FUNCTIONS */
	
	// computes equivalent distance error to the photoz error (PZerr) given to the constructor
	//double CompDistErr(); 
	
	/** Find the minimum and maximum x,y,z,zs (zs=phot or spec) of galaxies in 
	    catalog                                                               */
	double FindMinMaxCoords(); 
	
	/** Set grid over the galaxy data by specifying angle and redshift range 
	    grid needs to cover @warning probably doesn't work
	    @param Phi   angular range grid to cover
	    @param zl    minimum redshift grid to cover
	    @param zh    maximum redshift grid to cover
	    @param R	 pixel size of grid cell in Mpc                           */	   
	void SetGrid(double Phi=0, double zl=-1, double zh=-1,double R=8);
	
	/** Set grid over the galaxy data by specifying redshift of central pixel, 
	    size of the pixels, and (approx) number of pixels along each dimension
	    @param Nx    (approx) number of pixels along x-dimension
	    @param Ny    (approx) number of pixels along y-dimension
	    @param Nz    (approx) number of pixels along z-dimension
	    @param R     pixel size in Mpc
	    @param zref  redshift of central pixel                                */
	void SetGrid(int_8 Nx, int_8 Ny, int_8 Nz, double R, double zref); 

	/** Create a histogram of observed redshifts divided by true redshifts 
	    and save it to a text file 
	    @param SFTextFile    selection function is saved to file called [SFTextFile]_nofz.txt 
	    @param FullCat       name of FITS file containing z of all simulated gals
	                         if >1 file to read files are named FullCat_#ofn.fits
	    @param ZCol          name of column in FullCat containing true z
	    @param Nfiles        number of files n to read in  */
	void SaveSelecFunc(string SFTextFile, string FullCat, string ZCol="z");
	
	/** Set selection function                                                */
	inline void SetSelectionFunction(SelectionFunctionInterface& sf) { 
	    selfuncp_=&sf; sfcompute_=true; 
	    cout <<"    Set selection function"<<endl; };
	    
	/** Project the galaxy distribution into the grid, fill grid arrays
	    @param SkyArea    angle covered by observation cone                   */
	void GalGrid(double SkyArea = 999);
	
	/** Count up number of grid pixels that are within observed sky area as
	    defined by SkyArea variable                                           */
	sa_size_t ObsPixels();
	
	/** Convert row from catalog data table into a GalRecord                  */
	void Row2Record(DataTableRow& rowin, GalRecord& rec);
	
	/** Return true if galaxy accepted to be used in analysis. All galaxies are
	    currently accepted                                                    */
	bool Filter(GalRecord& rec); 
	
	/** Convert galaxy position in GalRecord (theta,phi,z) into Euclidean 
	    coordinates
	    @param rec        galaxy record
	    @param x          Euclidean x coordinate
	    @param y          Euclidean y coordinate
	    @param z          Euclidean z coordinate
	    @param redshift   redshift of galaxy                                  */
	void Rec2EuclidCoord(GalRecord& rec, double& x, double& y, double& z, double& redshift);
	
	/** Add galaxy with position x,y,z and weighting 1/phi to correct pixel in
	    galaxy grid
	    @param x          galaxy Euclidean x coordinate
	    @param y          galaxy Euclidean y coordinate
	    @param z          galaxy Euclidean z coordinate
	    @param phi        weight of galaxy                                    */
	void AddToCell(double x, double y, double z,double phi=1);
	
	/** Return cartesian coord (x,y,z) of pixel cell given pixel cell index (i,j,k) 
	    @param i    pixel index in x-dimension
	    @param j    pixel index in y-dimension
	    @param k    pixel index in z-dimension
	    @param x    pixel coordinate in x-dimension
	    @param y    pixel coordinate in y-dimension
	    @param z    pixel coordinate in z-dimension                           */
	void GetCellCoord(sa_size_t i, sa_size_t j, sa_size_t k, double& x, double& y, double& z);
	
	/** Return cartesian coord of pixel cell in the x-dimenstion given a pixel 
	    cell index 
	    @param i    pixel index in x-dimension                                */
	inline double GetCellX(sa_size_t i)  
	    { return Xmin_ + cellsize_/2 + i*cellsize_; };
	
	/** Return cartesian coord of pixel cell in the y-dimenstion given a pixel 
	    cell index 
	    @param j    pixel index in y-dimension                                */
	inline double GetCellY(sa_size_t j)  
	    { return Ymin_ + cellsize_/2 + j*cellsize_; };
	
	/** Return cartesian coord of pixel cell in the z-dimenstion given a pixel 
	    cell index 
	    @param k    pixel index in z-dimension                                */
	inline double GetCellZ(sa_size_t k)  
	    { return Zmin_ + cellsize_/2 + k*cellsize_; };
	
	/** Compute random weighted grid with same selection function as data 
	    @param nc         mean density of random grid
	    @param SaveArr    if true fill an array of pixel center redshifts     */
	void RandomGrid(double nc, bool SaveArr=true);
	
	/** Normalise the galaxy number and weighted galaxy number grids          */
	void NormNArrays(){ 		
		double normm=(double)n_obs_pixels_/ng_;
		ngals_ *= normm;
		double normw = (double)n_obs_pixels_/ngw_;
		//cout << wngals_ << endl << endl;
		wngals_*= normw;
		//cout << wngals_ << endl << endl;
		cout <<"     Normalising n-gals array by "<< normm << endl; 
		cout <<"     Normalising weighted n-gals array by "<< normw <<endl; 
		cout <<"     (number of observed pixels / number of gals in array)"<< endl;
		};
	
	/** Normalise galaxy random grid                                          */
	void NormRArray(){ 		
	    double normr = (double)(n_obs_pixels_/wnrand_)*( 1/sqrt(alph_));
		cout <<"     Normalising weighted random array by "<< normr <<endl;
		wrgals_*= normr; 
		VarianceRandomGrid(); };

	/** Add Gaussian errors to redshifts 
	    @param Err   redshift error size (Err*(1+z))                          */
	void SetGaussErr(double Err) {
		AddGaussErr_=true;
		PZDerr_ = Err; 
		if (AddGaussErr_)
			cout <<"    Gaussian errors WILL be added to z-coordinate"<<endl<<endl; 
		};
		
	/** Compute variance of random grid                                       */
	void VarianceRandomGrid() {
	    double mean,sigma;
		MeanSigma(wrgals_,mean,sigma);
		cout <<"    Weighted random array has mean = "<< mean;
		cout <<", variance = "<< sigma*sigma <<endl; };

	/** Extract a subset from full sized grid                                 
	    @param nsub    sub-grid extracted from full grid
	    @param dNp     width of sub-grid in z-dimension
	    @param Np      start sub-grid extraction at this pixel number (z-dim)
	    @param theta   maximum angle the sub-grid should extend over
	    @param af      array flag <0 return ngals, ==0 wngals, >0 wrgals      */
	vector<double> XtractSubArray(TArray<r_8>& nsub, sa_size_t dNp, sa_size_t Np, 
	                                                   double theta, int af=-1);
	
	/** Extract a subset from full sized grid  
	    @param nsub    sub-grid extracted from full grid
	    @param x1      start sub-grid extraction at this x-dim pixel number
	    @param x2      end sub-grid extraction at this x-dim pixel number
	    @param y1      start sub-grid extraction at this y-dim pixel number
	    @param y2      end sub-grid extraction at this y-dim pixel number
	    @param z1      start sub-grid extraction at this z-dim pixel number
	    @param z2      end sub-grid extraction at this z-dim pixel number
	    @param af      array flag <0 return ngals, ==0 wngals, >0 wrgals      */
	vector<double> XtractSubArray(TArray<r_8>& nsub, long x1, long x2, 
	                              long y1, long y2, long z1, long z2, int af=-1);
	                              
	/** Return number of pixels in x or y dimension that cover an angle of size 
	    theta at Np pixels along the grid's z-dimension
	    @param Np     number of pixels along grid's z-dimension
	    @param theta  angle                                                   */        
	sa_size_t GetNTransPix(sa_size_t Np, double theta);

    /** Output catalog in Euclidean coords (for debugging)                    */
	void OutputEuclidCat(double SkyArea);

    /** Write FITS header */
	void WriteHeader(string IncatName) {
		 if ( ngo_<1||nrand_<1 )
			throw ParmError("ERROR! Not ready to write Fits Header yet!");

		  fos_.WriteKey("NX",Nx_," axe transverse 1");
		  fos_.WriteKey("NY",Ny_," axe transverse 2");
		  fos_.WriteKey("NZ",Nz_," axe longitudinal (redshift)");
		  fos_.WriteKey("DX",cellsize_," Mpc");
		  fos_.WriteKey("DY",cellsize_," Mpc");
		  fos_.WriteKey("DZ",cellsize_," Mpc");
		  fos_.WriteKey("ZREF",zref_," reference redshift");	
		  fos_.WriteKey("ThRAD",SkyArea_," radius of circular sky area");
		  fos_.WriteKey("MeanOD",mean_overdensity_," mean dens SimLSS delta-fudge");
		  fos_.WriteKey("InCat",IncatName," original catalog");
		  fos_.WriteKey("NGALCAT",ngall_," N gals in InCat");
		  fos_.WriteKey("NOBSPIX",n_obs_pixels_," original catalog");
		  fos_.WriteKey("NGRID",ng_," N gals in grid");
		  fos_.WriteKey("NOGRID",ngout_," N gals outside grid");
		  fos_.WriteKey("NWGRID",ngw_," N weighted gals in grid");	
		  fos_.WriteKey("NRGRID",wnrand_," N (weighted if appl.) gals in random grid");
		};
	
	/** Zero the size of the galaxy grids                                     */
	void ZeroGalArrays() { wngals_.ZeroSize(); ngals_.ZeroSize();} ;
	
	/** Write the galaxy grids to a FITS file                                 */
	void WriteGalArrays() { ngals_.PackElements(); wngals_.PackElements();
				fos_ << ngals_; fos_ << wngals_;   };

    /** Save galaxy number grid to a FITS file (for debugging)                */
    void SaveNgArray(string outfileroot, string IncatName="unknown.fits") {
    
        string outfile;
        outfile = "tmp_" + outfileroot;//outfileroot +"_ngals.fits";
        cout <<"    Writing NGALS array to "<< outfile <<endl;
        FitsInOutFile fos(outfile, FitsInOutFile::Fits_Create);
                                                        
        fos << ngals_;
        //WriteHeader(fos,IncatName);
        fos.WriteKey("NGRID",ng_," N gals in grid");
        fos.WriteKey("NOGRID",ngout_," N gals outside grid");
        };

	
	//---- SHOULD BE MOVED OUT TO SELECTION FUNCTION CLASS ----//
	/*
	void ComputeSF(string,bool MCut=false); // computes selection function and sets sfcompute_=true
	void ComputeSFfromFile(string nzfile); // computes selection function from n(z) Histo read in from file
	void ComputeSFfromLF(double phistar,double Mstar,double alpha,double mlim, double Mc); 
	// computes SF look up table and sets sfcompute_=true
		void ApplyConstantSF2ngals(double sf)
		{ 
		ngals_*=sf; 
		cout <<"    Number of galaxies left after applying selection function = "<<ngals_.Sum()<<endl;
		cout <<"    Total number of galaxies = "<<ng_<<endl;
		cout <<"    Selection function = "<<sf<<endl;
		}
	void ApplyVaryingSF();
	void ApplyConstantSF2wngals(double sf){ wngals_*=sf;};
	void ApplyVaryingSF2wngals();
	void SetConstantSF(double sf);
	void SetVaryingSF(); */
	
	/* MINOR FUNCTIONS */
	
	/** Convert comoving distance to redshift                                 */
	double ConvCoD2z(double dc) { return dist2z_(dc); };
	
	/** Set printing level (not currently used)                               */
	inline void SetPrintLevel(int lev=0) { prtlev_=lev; };
	
	/** Return the approx volume of the survey: (maxx-minx)*(maxy-miny)*(maxz-minz)*/
	r_4 ReturnSurveyVolume(){return Vol_;};
	
	/** Return the volume of the grid                                         */
	r_4 ReturnGridVolume(){return volgrid_;};
	
	/** Return number of galaxies in galaxy number grid                       */
	sa_size_t ReturnNg(){return ng_;}; 
	
	/** Return number of galaxies in WEIGHTED galaxy number grid              */
	r_8 ReturnWg(){return ngw_;}; 
	
	/** Return number of galaxies in WEIGHTED random grid                     */
	r_8 ReturnRg(){return nrand_;}; 
	
	/** Return total number of galaxies in the simulation                     */
	sa_size_t ReturnNgAll(){return ngall_;};
	
	/** Return galaxy number grid (may or may not be normalised)              */	
	void ReturnNgals(TArray<r_8>& ngals) { ngals = ngals_; };
	
	/** Return weighted galaxy number grid (may or may not be normalised)     */
	void ReturnWngals(TArray<r_8>& wngals) { wngals = wngals_; };
	
	/** Return weighted random grid (may or may not be normalised)            */
	void ReturnWrgals(TArray<r_8>& wrgals) { wrgals = wrgals_; };
	
	/** Return alpha, number of galaxies in weighted galaxy number grid divided
	    by weighted random number grid                                        */
	double ReturnAlpha(){ return alph_;};
	
	/** Return z-dimension length of galaxy number grid                       */
	sa_size_t ReturnSizeZNgals(){ return ngals_.SizeZ();};
	
	/** Return photometric redshift error in comoving distance units          */ 
	double ReturnPZDerr(){return PZDerr_;};
	
	/** Return grid specification                                             */
	TVector<r_8> ReturnGridSpec(){
	    int row=4; TVector<r_8> gridspec(row); 
		gridspec(0)=Nx_; gridspec(1)=Ny_; 
		gridspec(2)=Nz_; gridspec(3)=cellsize_;
		return gridspec; };
		
	/** Set output filename root to output debugging files to                 */
	void SetDebugOutroot(string debugoutroot)
	    { debugoutroot_=debugoutroot; DoDebug_=true; };
	
	/** Print grids to screen for checking                                    */
	void CheckArrays() {
						cout << "    CHECK NGALS"<<endl;
						cout << ngals_<<endl;	
						cout << "    CHECK WNGALS"<<endl;
						cout << wngals_<<endl;
						cout << "    CHECK WRNGALS"<<endl;
						cout << wrgals_<<endl;	
						cout << endl;		
						};
	
		
protected:
	
	
	SwFitsDataTable& dt_;			  /**< data table containing galaxy catalog */
	SimpleUniverse& su_;			  /**< cosmology                            */
	RandomGenerator& rg_;			  /**< random number generator              */
	SelectionFunctionInterface* selfuncp_;    /**< selection function           */
	SInterp1D z2dist_;                /**< redshift to distance look up table   */
	SInterp1D dist2z_;		          /**< distance to redshift look up table   */
	
	bool DoDebug_;		  /**< true if debugging                              */
	bool AddGaussErr_;    /**< true if adding Gaussian error to redshifts     */
	bool sfcompute_;	  /**< true if selection function has been set        */
	bool RadialZ_;		  /**< if true z-dimension IS radial direction        */

	FitsInOutFile& fos_; 	/**< FITS file containing gridded galaxy data     */
	string debugoutroot_; 	/**< root file name to save things to when debugging */
	
	double PZErr_;      /**< size of photo-z error                            */
	double PZDerr_;		/**< size of photo-z error in comoving distance units */
	string ZOCol_;      /**< column name containing observed redshifts        */
	string ZSCol_;		/**< column name containing spec redshifts            */
	double volgrid_;    /**< volume of grid in Mpc^3                          */
	double cellsize_;	/**< pixel size of grid in Mpc                        */
	r_4 Vol_;           /**< approx volume of galaxy catalog                  */
	double SkyArea_;    /**< sky area: not clear on what this does?           */
	double xmin_;       /**< min x-coord of galaxies                          */
	double xmax_;       /**< max x-coord of galaxies                          */
	double ymin_;       /**< min y-coord of galaxies                          */
	double ymax_;       /**< max y-coord of galaxies                          */
	double zmin_;       /**< min z-coord of galaxies                          */
	double zmax_;       /**< max z-coord of galaxies                          */
	double Xmin_;       /**< min x-coord of grid                              */
	double Xmax_;       /**< max x-coord of grid                              */
	double Ymin_;       /**< min y-coord of grid                              */
	double Ymax_;       /**< max y-coord of grid                              */
	double Zmin_;       /**< min z-coord of grid                              */
	double Zmax_;       /**< max z-coord of grid                              */
	double zsmin_;      /**< min observed redshift of gals in catalog (could be spec or photo-z) */
	double zsmax_;		/**< max observed redshift of gals in catalog (could be spec or photo-z) */
	sa_size_t Nx_;      /**< number of pixels of grid in x-dimension          */    
	sa_size_t Ny_;      /**< number of pixels of grid in y-dimension          */ 
	sa_size_t Nz_;      /**< number of pixels of grid in z-dimension          */ 
	sa_size_t Npix_;	/**< total number of pixels in grid                   */ 
	sa_size_t ngo_;     /**< observed number of gals in sim                   */
	sa_size_t ngall_;	/**< total number of gals in sim                      */
	sa_size_t ng_;      /**< number of gals inside grid                       */
	sa_size_t ngout_;	/**< number of gals outside grid                      */ 
	sa_size_t ngout2_;  /**< number of gals outside weighted grid             */
	sa_size_t n_obs_pixels_; /**< number of pixels inside survey obs cone     */
	r_8 ngw_;           /**< number of gals inside weighted grid              */
	r_8 wnrand_;        /**< number of gals in weighted random grid           */
	r_8 nrand_;         /**< number of gals in random grid                    */
	int_8 idx_;         /**< index of center pixel in x-dimension             */
	int_8 idy_;         /**< index of center pixel in y-dimension             */
	int_8 idz_;	        /**< index of center pixel in z-dimension             */		  
	double DCref_;      /**< comoving distance to center pixel                */
	double zref_;	    /**< redshift of center pixel                         */  
	vector<double> zs_; /**< selection function look up table                 */
	vector<double> dL_; /**< selection function look up table                 */
	vector<double> phi_;/**< selection function look up table                 */
	vector<double> dc_; /**< selection function look up table                 */
	double alph_;		/**< number of gals in weighted galaxy number grid/by weighted random number grid */
	sa_size_t Ic1_;     /**< index of theta? column in data table             */
	sa_size_t Ic2_;     /**< index of phi? column in data table               */
	sa_size_t Izs_;     /**< index of spec-z column in data table             */
	sa_size_t Iz_;      /**< index of obs-z column in data table              */
	sa_size_t Iid_;     /**< index of gal id column in data table             */
	sa_size_t Iug_;     /**< index of u-g color column in data table          */
	sa_size_t Igr_;     /**< index of g-r color column in data table          */
	sa_size_t Iri_;     /**< index of r-i color column in data table          */
	sa_size_t Iiz_;     /**< index of i-z color column in data table          */
	sa_size_t Izy_;     /**< index of z-y color column in data tabl           */    
	int prtlev_;		/**< print level set                                  */
	//double phistar_,Mstar_,alpha_,mlim_,Mc_; // parameters for selection function ???REDUNDANT???
	double mean_overdensity_; /**< mean over-density of simlss grid AFTER setting cells with <-1 to =-1 */
	
	TArray<r_8> ngals_;     /**< number of gals in each grid pixel, not corrected for SF */
	TArray<r_8> wngals_;	/**< weighted number of gals in each grid pixel    */
	TArray<r_8> randomcat_;	/**< random galaxy number grid                     */
	TArray<r_8> wrgals_;	/**< weighted random number grid                   */
	TArray<r_8> weights_;	/**< weight value at grid pixel centers            */
	TArray<r_8> zc_;	    /**< redshift value at grid pixel centers          */
	
	SelectionFunctionInterface defselfunc_;   // Default selection function 
	FitsInOutFile fosdefault_;
};
