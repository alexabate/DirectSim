#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <typeinfo>
#include "fabtwriter.h"

#include "array.h"
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

#include "luc.h"
#include "constcosmo.h"
#include "sinterp.h"
#include "schechter.h"
#include "selectfunc.h"
#include "mass2gal.h"
#define PI 3.141592

/*

	Class GalRecord , Cat2Grid
	
	GalRecord class holds galaxy variables

	Cat2Grid class reads in a galaxy catalog (ra,dec,z) from a SwFitsDataTable which must have column names:
	"phi" and "theta", but the z column name can be set in the program

	Lays a grid across the galaxy catalog with grid cell size (default 8 Mpc) and 
	grid position specified in SetGrid function 
	OPTIONALLY: adds statistical photo-z error sig = PZerr*(1+z) to z-dimension (if PZerr>0)
	OPTIONALLY: corrects for selection effects (if SetSelectionFunction is called, which sets sfcompute_=true)

*/

class GalRecord {
public:
GalRecord() { alpha=delta=glong=glat=xcoo=ycoo=zcoo=0.; zs=zo=0.; type=1; u=g=r=i=z=y=24.; }

	void Print()
		{ cout <<"Printing GalRecord: ";
		  cout <<"alpha="<<alpha<<", delta="<<delta<<", zs="<<zs<<", zo="<<zo;}


	double alpha,delta,glong,glat; // angular coordinates
	double xcoo,ycoo,zcoo; // Euclidean coordinates
	double zs,zo; // SPECTRO z and OBSERVED z
	int type; // galaxy type
	double u,g,r,i,z,y; // ugrizy magnitudes
};


class Cat2Grid
{
public:
	/* CONSTRUCTORS */
	// Normal usage:
	Cat2Grid(SwFitsDataTable&,SimpleUniverse&,RandomGenerator&,FitsInOutFile& fos,string ZOCol="zp",string ZSCol="z",bool RadialZ=false, double PZErr=0,bool Print=true);
	// Just so you can use the Row2Record, Rec2EuclidCoord functions
	// this constructor takes the interp z->d functions as arguments rather than
	// calculating them every time.
	Cat2Grid(SwFitsDataTable&,SimpleUniverse&,RandomGenerator&,SInterp1D dist2z, SInterp1D z2dist,string ZOCol="zp",string ZSCol="zs",bool RadialZ=false);

	Cat2Grid(Cat2Grid const& a);
	virtual ~Cat2Grid();
	virtual Cat2Grid& Set(const Cat2Grid& a);
	
	
	/* MAIN FUNCTIONS */
	// computes equivalent distance error to the photoz error (PZerr) given to the constructor
	//double CompDistErr(); 
	// finds the minimum and maximum x,y,z,zs (zs=phot or spec) of galaxies in catalog
	double FindMinMaxCoords(); 	
	// set grid over data by specifying angle and redshift range it must cover	   
	void SetGrid(double Phi=0, double zl=-1, double zh=-1,double R=8);// probably not usable
	// set grid over data by specifying redshift of central pixel and number of pixels
	void SetGrid(int_8 Nx, int_8 Ny, int_8 Nz,double R, double zref); 
	// saves selection function to a text file
	void SaveSelecFunc(string SFTextFile, string FullCat, string ZCol="z",int Nfiles=1); 
	// Set selection function
	inline void SetSelectionFunction(SelectionFunctionInterface& sf) { selfuncp_=&sf;
									   sfcompute_=true; 
									   cout <<"    Set selection function"<<endl; }	
	// projects the galaxy distribution into the grid:
	void GalGrid(double SkyArea = 999);
	// convert row in catalog table to GalRecord
	void Row2Record(DataTableRow& rowin, GalRecord& rec);
	// returns true if galaxy accepted to be used in analysis
	bool Filter(GalRecord& rec); // all currently accepted
	// convert galaxy position in GalRecord (theta,phi,z) to Euclidean coordinates
	void Rec2EuclidCoord(GalRecord& rec, double& x, double& y, double& z, double& redshift);
	// add gal with position x,y,z and weighting 1/phi to galaxy grid
	void AddToCell(double x, double y, double z,double phi=1);
	// return cartesian coord (x,y,z) of pixel cell given pixel cell index (i,j,k) 
	void GetCellCoord(sa_size_t i, sa_size_t j, sa_size_t k, double& x, double& y, double& z);
	inline double GetCellX(sa_size_t i)  { return Xmin_ + cellsize_/2 + i*cellsize_; }
	inline double GetCellY(sa_size_t j)  { return Ymin_ + cellsize_/2 + j*cellsize_; }
	inline double GetCellZ(sa_size_t k)  { return Zmin_ + cellsize_/2 + k*cellsize_; }
	// count up observed pixels
	sa_size_t ObsPixels();
	// create random grid
	void RandomGrid(double nc,bool SaveArr=true);
	
	// normalise galaxy number and number weighted arrays
	void NormNArrays(){ 		
		double normm=(double)n_obs_pixels_/ng_;
		ngals_ *= normm;
		double normw = (double)n_obs_pixels_/ngw_;
		wngals_*= normw;
		cout <<"     Normalising n-gals array by "<<normm<<endl;
		cout <<"     Normalising weighted n-gals array by "<<normw<<endl; };
	
	// normalise galaxy random array
	void NormRArray(){ 		double normr = (double)(n_obs_pixels_/wnrand_)*( 1/sqrt(alph_));
					cout <<"     Normalising weighted random array by "<<normr<<endl;
					wrgals_*= normr; 
					VarianceRandomGrid(); };

	// add Gaussian errors to redshifts
	void SetGaussErr(double Err)
		{
		AddGaussErr_=true;
		PZDerr_ = Err; 
		if (AddGaussErr_)
			cout <<"    Gaussian errors WILL be added to z-coordinate"<<endl<<endl; 
		};
	// variance of random grid
	void VarianceRandomGrid()
			{ double mean,sigma;
			  MeanSigma(wrgals_,mean,sigma);
			  cout <<"    Weighted random array has mean = "<<mean<<", variance = "<<sigma*sigma<<endl; };


	// extract a sub-array from full array
	vector<double> XtractSubArray(TArray<r_8>& nsub,sa_size_t dNp,sa_size_t Np, double theta, int af=-1);
	vector<double> XtractSubArray(TArray<r_8>& nsub,long x1,long x2,long y1,long y2,long z1,long z2,int af=-1);
	sa_size_t GetNTransPix(sa_size_t Np, double theta);

	void OutputEuclidCat(double SkyArea);

	void WriteHeader(string IncatName)
		{
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
		  fos_.WriteKey("MeanOverDensity",mean_overdensity_," mean dens SimLSS delta-fudge");
		  fos_.WriteKey("InCat",IncatName," original catalog");
		  fos_.WriteKey("NGALCAT",ngall_," N gals in InCat");
		  fos_.WriteKey("NOBSPIX",n_obs_pixels_," original catalog");
		  fos_.WriteKey("NGRID",ng_," N gals in grid");
		  fos_.WriteKey("NOGRID",ngout_," N gals outside grid");
		  fos_.WriteKey("NWGRID",ngw_," N weighted gals in grid");	
		  fos_.WriteKey("NRGRID",wnrand_," N (weighted if appl.) gals in random grid");
		};
	
	void ZeroGalArrays() { wngals_.ZeroSize(); ngals_.ZeroSize();} ;
	void WriteGalArrays() { ngals_.PackElements(); wngals_.PackElements();
				fos_ << ngals_; fos_ << wngals_;   };

        void SaveNgArray(string outfileroot,string IncatName="unknown.fits")  
                                                        {
                                                        string outfile;
                                                        outfile = "tmp_"+outfileroot;//outfileroot +"_ngals.fits";
                                                        cout <<"    Writing NGALS array to "<<outfile<<endl;
                                                        FitsInOutFile fos(outfile, FitsInOutFile::Fits_Create);
                                                        
                                                        fos << ngals_;
                                                        //WriteHeader(fos,IncatName);
                                                        fos.WriteKey("NGRID",ng_," N gals in grid");
                                                        fos.WriteKey("NOGRID",ngout_," N gals outside grid");
                                                        }

	
	// SHOULD BE MOVED OUT TO SELECTION FUNCTION CLASS
	void ComputeSF(string,bool MCut=false); // computes selection function and sets sfcompute_=true
	void ComputeSFfromFile(string nzfile); // computes selection function from n(z) Histo read in from file
	void ComputeSFfromLF(double phistar,double Mstar,double alpha,double mlim, double Mc); // computes SF look up table and sets sfcompute_=true
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
	void SetVaryingSF();
	
	/* MINOR FUNCTIONS */
	// convert comoving distance to redshift
	double ConvCoD2z(double dc) { return dist2z_(dc); };
	// set printing level
	inline void SetPrintLevel(int lev=0) { prtlev_=lev; }// has no use at this time
	// return approx volume of survey: (maxx-minx)*(maxy-miny)*(maxz-minz)
	r_4 ReturnSurveyVolume(){return Vol_;};
	// return volume of grid
	r_4 ReturnGridVolume(){return volgrid_;};
	// return number of gals in ngals grid
	sa_size_t ReturnNg(){return ng_;}; 
	// return number of gals in WEIGHTED ngals grid
	r_8 ReturnWg(){return ngw_;}; 
	// return number of gals in RANDOM WEIGHTED ngals grid
	r_8 ReturnRg(){return nrand_;}; 
	// return number of gals in simulation
	sa_size_t ReturnNgAll(){return ngall_;};
	// returned gridded data (may or may not be normalised)
	void ReturnNgals(TArray<r_8>& ngals) { ngals = ngals_; };
	void ReturnWngals(TArray<r_8>& wngals) { wngals = wngals_; };
	void ReturnWrgals(TArray<r_8>& wrgals) { wrgals = wrgals_; };
	// return alpha
	double ReturnAlpha(){ return alph_;}
	// return galaxy number grid z-dimension length
	sa_size_t ReturnSizeZNgals(){ return ngals_.SizeZ();};
	// return photometric redshift error in comoving distance 
	double ReturnPZDerr(){return PZDerr_;};
	// return grid specification
	TVector<r_8> ReturnGridSpec(){int row=4; TVector<r_8> gridspec(row); 
										gridspec(0)=Nx_; gridspec(1)=Ny_; 
										gridspec(2)=Nz_; gridspec(3)=cellsize_;
										return gridspec;
								  };
	// set output filename root to output files for debugging to
	void SetDebugOutroot(string debugoutroot){debugoutroot_=debugoutroot; DoDebug_=true;};
	// print arrays to screen
	void CheckArrays() {
						cout << "    CHECK NGALS"<<endl;
						cout << ngals_<<endl;	
						cout << "    CHECK WNGALS"<<endl;
						cout << wngals_<<endl;
						cout << "    CHECK WRNGALS"<<endl;
						cout << wrgals_<<endl;	
						cout << endl;		
						};
	
		
	/* CLASS VARIABLES */
	
	/* objects */
	SwFitsDataTable& dt_;			  // address of data table containing galaxy catalog
	SimpleUniverse& su_;			  // holds cosmological parameters
	RandomGenerator& rg_;			  // allows random number generation
	SelectionFunctionInterface defselfunc_;   // Default selection function 
	SelectionFunctionInterface* selfuncp_;    // The selection function object
	SInterp1D z2dist_,dist2z_;		  //  redshift to distance look up
	FitsInOutFile fosdefault_;
	
	/* bools */
	bool DoDebug_;		// if debugging
	bool AddGaussErr_;      // if adding Gaussian error
	bool sfcompute_;	// if selection function has been set
	bool RadialZ_;		// if true z-dimension IS radial direction

	/* outputs */
	FitsInOutFile& fos_; 	// output fits file
	string debugoutroot_; 	// root file name to save things to when debugging [DEFAULT=tmp]
	
	/* simple */
	double PZErr_, PZDerr_;			  // statistical photo-z error added, equivalent distance error
	string ZOCol_,ZSCol_;		  // column names of redshifts to read in 
	double volgrid_, cellsize_;		  // volume of GRID, pixel size of grid
	r_4 Vol_;						  // volume of SURVEY
	double SkyArea_;
	double xmin_,xmax_,ymin_,ymax_,zmin_,zmax_; // min and max coords of all gals
	double Xmin_,Xmax_,Ymin_,Ymax_,Zmin_,Zmax_; // min and max coords of grid
	double zsmin_, zsmax_;			  // min and max OBSERVED redshift of gals in cat (depends on value of ZOCol_)
	sa_size_t Nx_, Ny_, Nz_,Npix_;	  // number of pixels in grid
	sa_size_t ngo_,ngall_;			  // OBSERVED number of gals in sim, ALL gals in sim,
	sa_size_t ng_,ngout_,ngout2_;	  // Number of gals in grid, n gals outside grid, n gals outside weighted grid
	sa_size_t n_obs_pixels_;	  // Number of pixels inside survey cone
	r_8 ngw_,wnrand_,nrand_;	  // Number of gals in weighted grid,Number of gals in weighted random grid, Number of gals in random grid
	int_8 idx_,idy_,idz_;			  // indices of center pixel
	double DCref_,zref_;			  // comoving distance of center pixel
	vector<double> zs_, dL_, phi_, dc_; // selection function look up table
	double alph_;					  // = ngw_/wnrand_	
	sa_size_t Ic1_,Ic2_,Izs_,Iz_,Iid_;	  // indices of quantities in data table: theta,phi,zspec,z?
	sa_size_t Iug_,Igr_,Iri_ ,Iiz_,Izy_;// indices of quantities in data table: colors
	int prtlev_;					  // print level set by SetPrintLevel()
	double phistar_,Mstar_,alpha_,mlim_,Mc_; // parameters for selection function ???REDUNDANT???
	double mean_overdensity_; // mean over-density of simlss grid AFTER setting cells with <-1 to =-1
	
	/* arrays */
	TArray<r_8> ngals_;			  // number of gals in each grid pixel, not corrected for SF
	TArray<r_8> wngals_;			  // weighted number of gals in each grid pixel
	TArray<r_8> randomcat_;			  // random catalog grid
	TArray<r_8> wrgals_;			  // weighted random catalog grid
	TArray<r_8> weights_;			  // weight value at grid pixel center
	TArray<r_8> zc_;			  // redshift value at grid pixel center
	
	
};
