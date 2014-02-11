#include "cat2grid.h"
#include "ctimer.h"


//Improved constructor: 
Cat2Grid::Cat2Grid(SwFitsDataTable& dt, SimpleUniverse& su, RandomGenerator& rg,
                   FitsInOutFile& fos, string ZOCol, string ZSCol,
                   bool RadialZ, double PZErr, bool Print) 
: dt_(dt) , su_(su) , rg_(rg) , fos_(fos) , defselfunc_(1.), selfuncp_(&defselfunc_) 
   , ZOCol_(ZOCol) , ZSCol_(ZSCol) , RadialZ_(RadialZ) , PZErr_(PZErr)
// Constructor for this class reads in the theta, phi, redshift of all galaxies in the catalog
// from data table dt
// su:	  class which holds the cosmological parameters and can be used to 
//        calculate cosmological quantities
// rg:	  class used to generate random numbers
// ZOCol: col name of OBSERVED redshifts to read in, could be spec-z, gauss-z, phot-z
// ZSCol: col name of SPECTRO redshifts to be read in (must be spec-z)
// RadialZ: if true z-dimension IS radial direction
// PZErr: size of photometric redshift errors: sigz = PZErr(1+z); 
// Print: if true prints some extra info to screen whilst in constructor
{
    if (Print)
		cout <<endl<<"    In Cat2Grid constructor ..."<<endl;
	
	DVList dvl;
	dvl = dt.Info();
	mean_overdensity_ = dvl["MeanOverDensity"];
	cout << "    Mean over-density of fudged SimLSS grid = "<< mean_overdensity_ <<endl;
	
	DoDebug_ = false; // flag to do debug
	debugoutroot_ = "tmp"; // if want to debug, files are written to filenames with this root
	sfcompute_ = false; // flag becomes true after selection function is set
	PZDerr_=0;// set photoz error in Mpc equal to zero to start with
	AddGaussErr_ = false; // flag to add Gaussian photo-z err

	// Set counting galaxies to ZERO
	ng_=0,ngo_=0; // num gals INSIDE grid, num gals OBSERVED
	ngw_=0;// num "gals" INSIDE WEIGHTED grid
	ngout_=0; // number of galaxies OUTSIDE grid
	wnrand_=0,nrand_=0; // weighted/unweighted number of "galaxies" in random grid
	
	// Data table column names
	if (RadialZ_) {
		Ic1_ = dt_.IndexNom("x");
		Ic2_ = dt_.IndexNom("y");
		}
	else {
		Ic1_ = dt_.IndexNom("theta");
		Ic2_ = dt_.IndexNom("phi");
		}
	Izs_ = dt_.IndexNom(ZSCol_); // the SPEC-Z column
	
	if (Print) {
		cout << "    Reading SPECTRO redshifts from column name "<< ZSCol_ <<endl;
		cout << "    Reading OBSERVED redshifts from column name "<< ZOCol_ <<endl;
		}
	Iz_ = dt_.IndexNom(ZOCol_);
	Iid_= dt_.IndexNom("GalID");
	
	ngall_=dt_.NEntry();
	if (Print) {
		cout <<"    TOTAL number of galaxies in catalog ="<< ngall_ <<endl;
		cout <<"    Number of columns in catalog = "<< dt_.NCols() <<endl;
		if (RadialZ_)
			cout <<"    Coord x = column "<< Ic1_ <<", coord y = column "<< Ic2_ <<endl;
		else
			cout <<"    Angle theta = column "<< Ic1_ <<", angle phi = column "<< Ic2_ <<endl;
		cout <<"    Spec-z = column "<< Izs_ <<", z used in analysis = column "<< Iz_ <<endl;
		cout <<"    (the above will be the same column when reading in spec-z)"<<endl;
		}
		
	// Create a redshift - distance look up table to save time not using su_ class computation
	int_8 nz=1000000;
	vector<double> zrs, codist;
	double minz=0, maxz=10;
	double dz = (maxz-minz)/(nz-1);
	for(int kk=0; kk<nz; kk++) {
		
		  double zs=minz+kk*dz;
		  su_.SetEmissionRedShift(zs);
		  double cod =su_.RadialCoordinateMpc(); // radial distance 
		  zrs.push_back(zs);
		  codist.push_back(cod); 
		}
	z2dist_.DefinePoints(zrs,codist,minz,maxz,2*nz);
	dist2z_.DefinePoints(codist,zrs,codist[0],codist[codist.size()-1],2*nz);

    SetPrintLevel();
	if (Print)
		cout <<"    Exit Cat2Grid constructor ..."<<endl<<endl;
};


// make interpolation table outside Cat2Grid, just use Row2Record, Rec2EuclidCoord functions
Cat2Grid::Cat2Grid(SwFitsDataTable& dt,SimpleUniverse& su,RandomGenerator& rg,
      SInterp1D dist2z, SInterp1D z2dist,string ZOCol,string ZSCol,bool RadialZ) 
: dt_(dt) , su_(su) , rg_(rg) , defselfunc_(1.), selfuncp_(&defselfunc_) , fos_(fosdefault_) ,
   ZOCol_(ZOCol) , ZSCol_(ZSCol) , dist2z_(dist2z) , z2dist_(z2dist) , RadialZ_(RadialZ)
{

	bool Print=false;
	if (Print)
		cout <<endl<<"    In Cat2Grid constructor ..."<<endl;
		
	PZErr_=0;
	
	
	DoDebug_ = false; // flag to do debug
	debugoutroot_ = "tmp"; // if want to debug, files are written to filenames with this root
	sfcompute_ = false; // flag becomes true after selection function is set
	PZDerr_=0;// set photoz error in Mpc equal to zero to start with
	AddGaussErr_ = false; // flag to add Gaussian photo-z err
	
	// Data table column names
	if (RadialZ_) {
		Ic1_ = dt_.IndexNom("x");
		Ic2_ = dt_.IndexNom("y");
		}
	else {
		Ic1_ = dt_.IndexNom("theta");
		Ic2_ = dt_.IndexNom("phi");
		}
	Izs_ = dt_.IndexNom(ZSCol_); // the SPEC-Z column
	Iz_ = dt_.IndexNom(ZOCol_);
	if (Print)
		cout << "    Reading SPECTRO redshifts from column name "<< ZSCol_ <<endl;
	if (Print)
		cout << "    Reading OBSERVED redshifts from column name "<< ZOCol_ <<endl;

	
	ngall_=dt_.NEntry();
	if (Print) {
		cout <<"    TOTAL number of galaxies in catalog ="<< ngall_ <<endl;
		cout <<"    Number of columns in catalog = "<< dt_.NCols() <<endl;
		if (RadialZ_)
			cout <<"    Coord x = column "<< Ic1_ <<", coord y = column "<< Ic2_ <<endl;
		else
			cout <<"    Angle theta = column "<< Ic1_ <<", angle phi = column "<< Ic2_ <<endl;
		cout <<"    Spec-z = column "<< Izs_ <<", z used in analysis = column "<< Iz_ <<endl;
		cout <<"    (the above will be the same column when reading in spec-z)"<<endl;
		}
		
    	SetPrintLevel();
	if (Print)
		cout <<"    Exit Cat2Grid constructor ..."<<endl<<endl;
};


// copy constructor
Cat2Grid::Cat2Grid(Cat2Grid const& a)
: dt_(a.dt_) , su_(a.su_) , rg_(a.rg_), fos_(a.fos_) ,
  defselfunc_(a.defselfunc_.sfv_), selfuncp_(a.selfuncp_)
{
	cout <<"    Cat2Grid COPY constructor"<<endl;
	Set(a);
};


// Copy Cat2Grid 
Cat2Grid& Cat2Grid::Set(const Cat2Grid& a)
{
	//Initialised in constructor
					
	PZErr_=a.PZErr_;
	debugoutroot_=a.debugoutroot_;	
	PZDerr_=a.PZDerr_;	
	AddGaussErr_=a.AddGaussErr_;
	selfuncp_=a.selfuncp_;
	ZSCol_=a.ZSCol_;
	ZOCol_=a.ZOCol_;

	//Initialised in FindMinMaxCoords 
	zsmin_=a.zsmin_;
	zsmax_=a.zsmax_;	
	xmin_=a.xmin_;
	xmax_=a.xmax_;
	ymin_=a.ymin_;
	ymax_=a.ymax_;
	zmin_=a.zmin_;
	zmax_=a.zmax_;
		 
	ng_=a.ng_;		
	
	//Initialised in ComputeSFfromLF
	/*phistar_=a.phistar_;
	Mstar_=a.Mstar_;
	alpha_=a.alpha_;
	mlim_=a.mlim_;
	Mc_=a.Mc_;*/
	
	//Initialised in ComputeSFfromLF or ComputeSF
	if (a.dL_.size()>0)
		{
		dL_=a.dL_;
		phi_=a.phi_;
		dc_=a.dc_;
		}
	
	//Initialised in ComputeSF
	if (a.zs_.size()>0)
		zs_=a.zs_;
	
	//Initialised in ConstructGrid
	Vol_=a.Vol_;	
	Nx_=a.Nx_;
	Ny_=a.Ny_;
	Nz_=a.Nz_;
	volgrid_=a.volgrid_;
	cellsize_=a.cellsize_;
	
	Npix_=a.Npix_;
	if (a.ngals_.Size()>0)
		ngals_=a.ngals_;
	if (a.randomcat_.Size()>0)
		{
		randomcat_=a.randomcat_;
		wngals_=a.wngals_;
		weights_=a.weights_;
		wrgals_=a.wrgals_;
		}
	alph_=a.alph_;
	
		
	prtlev_ = a.prtlev_;
	
};


double Cat2Grid::FindMinMaxCoords()
// Finds min and max of cartesian x,y,z coords
// of all galaxies in simulation
// Find min/max redshift
{
	cout <<endl<<"    Cat2Grid::FindMinMaxCoords()"<<endl;
	
	if (AddGaussErr_) {
		cout <<"    Remember: adding Gaussian errors, equivalent distance error = ";
		cout << PZDerr_ <<" Mpc"<<endl;
		}
		
	long seed=1;
	rg_.SetSeed(seed);
	
	sa_size_t ng=0; // reset to zero to count up all OBSERVED galaxies
	
	// artificially set to very large min, very small max
	double minx=1e10,miny=1e10,minz=1e10;
	double maxx=0,maxy=0,maxz=0;
	double minzs=20,maxzs=0;
	
	DataTableRow rowin = dt_.EmptyRow();
	GalRecord grec;
	double x=1e8,y=1e8,z=-1e8,redshift=-10;
	
	for(long ig=0; ig<ngall_; ig++) {
		
	  	dt_.GetRow(ig, rowin);
	  	Row2Record(rowin,grec);
	  	if (!Filter(grec)) continue;
	  	Rec2EuclidCoord(grec,x,y,z,redshift);
	  	// the redshift returned here is the redshift to be used
		// in the analysis, could be spec-z, phot-z, gauss-z
	
	  	// print out first 10 gals
		if(ig<10)	
			{
			cout <<"    galid="<<ig<<": ";
			grec.Print();
			cout <<", x="<<x<<", y="<<y<<", z="<<z;
			if (AddGaussErr_)
				cout <<", zpG="<<redshift<<endl;
			else
				cout <<endl;
			}
	  	ng++; 
		
		// perform check to find minimums/maximums	
		if(x<minx)
			minx = x;
		if(y<miny)
			miny = y;
		if(z<minz)
			minz = z;
			
		if(x>maxx)
			maxx = x;
		if(y>maxy)
			maxy = y;
		if(z>maxz)
			maxz = z;
			
		if (redshift<minzs)
			minzs = redshift;
		if (redshift>maxzs)
			maxzs = redshift;

		}
		
	cout <<"    Number of observed galaxies = "<< ng;
	cout <<", Total number of galaxies = "<< ngall_ <<endl;
	
	if (ng>ngall_)
		throw ParmError("ERROR! galaxy number discrepancy");
	
	// set min and max to class variables
	xmin_ = minx;
	xmax_ = maxx;
	ymin_ = miny;
	ymax_ = maxy;
	zmin_ = minz;
	zmax_ = maxz;
	zsmin_ = minzs;
	zsmax_= maxzs;
		
	// find maximum (obs) luminosity distance in survey
	su_.SetEmissionRedShift(zsmax_);
	double maxdL=su_.LuminosityDistanceMpc();
		
	cout <<"    compute approx SURVEY volume .... "<<endl;
	cout <<"    xmin="<< xmin_ <<", xmax="<< xmax_ <<endl;
	cout <<"    ymin="<< ymin_ <<", ymax="<< ymax_ <<endl;
	cout <<"    zmin="<< zmin_ <<", zmax="<< zmax_ <<endl;
	Vol_ = (xmax_-xmin_)*(ymax_-ymin_)*(zmax_-zmin_);
	cout <<"    -> Volume = "<< Vol_ <<" Mpc"<<endl;
	cout <<"    OBSERVED redshift range (from column "<< ZOCol_ <<"): ";
	cout << zsmin_ <<"<z<"<< zsmax_ <<endl;
	cout <<"    Maximum dL="<< maxdL <<endl;
	
	cout <<"    EXIT Cat2Grid::FindMinMaxCoords()"<<endl<<endl;
	
	return maxdL;
	
};


void Cat2Grid::SetGrid(double Phi, double zl, double zh,double R)
// Computes range of grid to lay over galaxy simulation
// from specified survey volume
{
	cout <<endl<<"    Cat2Grid::SetGrid()"<<endl;
	cellsize_=R;
	
	cout << "    Using survey geometry to define grid"<<endl;
	cout << "    Grid to cover: "<< zl <<"<z<"<< zh <<", Phi = "<< Phi <<" radians"<<endl;
	cout << "    Pixel size = "<< cellsize_ <<endl;
	// To find minx,maxx,miny,maxy need to find size Theta extends
	// at the maximum redshift (zh).  Since x=y=0 at center of survey:
	su_.SetEmissionRedShift(zh);
	Xmin_ = -su_.RadialCoordinateMpc()*tan(Phi);
	Ymin_ = Xmin_;
	Xmax_ = -Xmin_;
	Ymax_ = -Ymin_;
	// Min z value is just distance to min redshift (zl) at CENTER of survey
	su_.SetEmissionRedShift(zl);
	Zmin_=su_.RadialCoordinateMpc();
	// Max z value is distance to max redshift at EDGE of survey
	su_.SetEmissionRedShift(zh);
	Zmax_=su_.RadialCoordinateMpc()/cos(Phi);
	cout <<"    CHECK Cartesian survey boundaries: "<< Xmin_ <<"<X<"<< Xmax_;
	cout <<", "<< Ymin_ <<"<Y<"<< Ymax_ <<", "<< Zmin_ <<"<Z<"<< Zmax_ <<endl;
	
	// calculate EXACT number of cells that reach from one side of sim to the other (including remainder)
	double Lx, Ly, Lz; 
	Lx = (Xmax_-Xmin_)/cellsize_;
	Ly = (Ymax_-Ymin_)/cellsize_;
	Lz = (Zmax_-Zmin_)/cellsize_;
	cout <<"    Number of cells which cover survey in each direction"<<endl;
	cout <<"    Nx="<< Lx <<", Ny="<< Ly <<", Nz="<< Lz <<endl;
	// round DOWN to get INTEGER number of cells that cover each dimension
	// round down because don't want a pixel going outside of survey area
	// integer number cells needed
	Nx_ = (sa_size_t)floor(Lx); Ny_ = (sa_size_t)floor(Ly); Nz_ = (sa_size_t)floor(Lz); 
	cout <<"    Integer number of cells: Nx="<< Nx_ <<", Ny="<< Ny_ <<", Nz="<< Nz_ <<endl;
	
	double idzd = ceil((double)Nz_/2); // index of center pixel in z direction
	double idxd = ceil((double)Nx_/2);
	double idyd = ceil((double)Ny_/2);	
	idx_ = idxd, idy_ = idyd,idz_ = idzd;
	cout <<"    Indices of center pixels: Ix="<< idx_ <<", Iy="<< idy_ <<", Iz="<< idz_ <<endl;
	DCref_ = Zmin_+cellsize_/2+(idz_-1)*cellsize_;
	cout <<"CHECK: comoving distance = "<< DCref_ <<endl;
	
	// Need to recompute Xmin_ etc to EXACTLY match grid edges
	Xmin_ = -(idx_-1)*cellsize_-cellsize_/2;
	Xmax_ = (idx_-1)*cellsize_+cellsize_/2;
	Ymin_ = -(idy_-1)*cellsize_-cellsize_/2;
	Ymax_ = (idy_-1)*cellsize_+cellsize_/2;
	Zmin_ = DCref_-(idz_-1)*cellsize_-cellsize_/2;
	Zmax_ = DCref_+(idz_-1)*cellsize_+cellsize_/2;
	
	cout <<"    Initial survey boundaries: "<< Xmin_ <<"<X<"<< Xmax_;
	cout <<", "<< Ymin_ <<"<Y<"<< Ymax_ <<", "<< Zmin_ <<"<Z<"<< Zmax_ <<endl;

//	// round down to nearest integer SHOULD REMOVE THIS
//	Xmin_=floor(Xmin_);
//	Xmax_=floor(Xmax_);
//	Ymin_=floor(Ymin_);
//	Ymax_=floor(Ymax_);
//	Zmin_=floor(Zmin_);
//	Zmax_=floor(Zmax_);
//	cout <<"    Rounded survey boundaries: "<<Xmin_<<"<X<"<<Xmax_<<", "<<Ymin_<<"<Y<"<<Ymax_<<", "<<Zmin_<<"<Z<"<<Zmax_<<endl;

	DCref_ = Zmin_+R/2+(idz_-1)*R;
	
	// convert DCref_ to a redshift
	int_8 nz=1000;
	vector<double> zrs, codist;
	double minz=0, maxz=10;
	double dz = (maxz-minz)/(nz-1);
	for(int kk=0; kk<nz; kk++) {
		  double zs=minz+kk*dz;
		  su_.SetEmissionRedShift(zs);
		  double cod =su_.RadialCoordinateMpc(); // radial distance 
		  zrs.push_back(zs);
		  codist.push_back(cod); 
		}
	double mind = codist[0];
	double maxd = codist[codist.size()-1];
	SInterp1D dist2z(codist,zrs,mind,maxd,2*nz);
	zref_ = dist2z(DCref_);
	cout <<"    Redshift of central pixel = "<< zref_;
	cout <<", found from comoving distance = "<< DCref_ <<endl;
	
	// INITIALISE ARRAYS TO GRID SIZE
	int ndim=3;
	sa_size_t mydim[ndim];
	mydim[0]=Nx_; mydim[1]=Ny_; mydim[2]=Nz_;
	ngals_.SetSize(ndim, mydim);	// number of galaxies in each cell: data
	wngals_.SetSize(ndim, mydim);	// weighted galaxy density field: data
	if (sfcompute_) {
		randomcat_.SetSize(ndim, mydim);// number of galaxies in each cell: random catalog
		weights_.SetSize(ndim, mydim);	// weights of galaxy density field (per cell)
		wrgals_.SetSize(ndim, mydim);	// weighted galaxy density field: random catalog
		}
	
	cout <<"    EXIT Cat2Grid::SetGrid()"<<endl<<endl<<endl;
};


void Cat2Grid::SetGrid(int_8 Nx, int_8 Ny, int_8 Nz, double R, double zref)
// Computes range of grid to lay over galaxy simulation
// from SimLSS cube parameters
{
	cout <<endl<<"    Cat2Grid::SetGrid()"<<endl;
	
	cout << "    Using SimLSS grid to define grid"<<endl;
	cout << "    SimLSS grid: "<< Nx <<","<< Ny <<","<< Nz <<" pixels"<<endl;
	cout << "    Centered at redshift z = "<< zref <<endl;
	cout << "    Resolution R = "<< R <<endl;
	
	cellsize_ = R;
	zref_ = zref;
	
	su_.SetEmissionRedShift(zref_);
	DCref_ = su_.RadialCoordinateMpc();
	double idzd = ceil((double)Nz/2); // index of center pixel in z direction
	double idxd = ceil((double)Nx/2);
	double idyd = ceil((double)Ny/2);	
	idx_ = idxd, idy_ = idyd,idz_ = idzd;
	cout <<"    Indices of center pixels: Ix="<< idx_ <<", Iy="<< idy_ <<", Iz="<< idz_ <<endl;
	
	Xmin_ = -(idx_-1)*cellsize_-cellsize_/2;
	Xmax_ = (idx_-1)*cellsize_+cellsize_/2;
	Ymin_ = -(idy_-1)*cellsize_-cellsize_/2;
	Ymax_ = (idy_-1)*cellsize_+cellsize_/2;
	Zmin_ = DCref_-(idz_-1)*cellsize_-cellsize_/2;
	Zmax_ = DCref_+(idz_-1)*cellsize_+cellsize_/2;
	
	// cell coords were given by:
	//X = (i-(idmidx_-1))*Dx_; }
	//(j-(idmidy_-1))*Dy_; }
	//(k-(idmidz_-1))*Dz_+DCref_; }
	
	cout <<"    Initial survey boundaries: "<< Xmin_ <<"<X<"<< Xmax_;
	cout <<", "<< Ymin_ <<"<Y<"<< Ymax_ <<", "<< Zmin_ <<"<Z<"<< Zmax_ <<endl;

//	// round down to nearest integer SHOULD REMOVE THIS
//	Xmin_=floor(Xmin_);
//	Xmax_=floor(Xmax_);
//	Ymin_=floor(Ymin_);
//	Ymax_=floor(Ymax_);
//	Zmin_=floor(Zmin_);
//	Zmax_=floor(Zmax_);
//	cout <<"    Rounded grid boundaries: "<<Xmin_<<"<X<"<<Xmax_<<", "<<Ymin_<<"<Y<"<<Ymax_<<", "<<Zmin_<<"<Z<"<<Zmax_<<endl;
	
	// calculate EXACT number of cells that reach from one side of sim to the other (including remainder)
	
	double Lx, Ly, Lz; 
	Lx = (Xmax_-Xmin_)/cellsize_;
	Ly = (Ymax_-Ymin_)/cellsize_;
	Lz = (Zmax_-Zmin_)/cellsize_;
	cout <<"    Number of cells which cover survey in each direction"<<endl;
	cout <<"    Nx="<< Lx <<", Ny="<< Ly <<", Nz="<< Lz <<endl;
	// round up to get INTEGER number of cells that cover each dimension
	// integer number cells needed
	Nx_ = (sa_size_t)ceil(Lx); Ny_ = (sa_size_t)ceil(Ly); Nz_ = (sa_size_t)ceil(Lz); 
	cout <<"    Integer number of cells: Nx="<<Nx_<<", Ny="<<Ny_<<", Nz="<<Nz_<<endl;
	
	if (Nx_!=Nx) {
		cout << "TEMP FUDGE, setting Nx to input value"<<endl;
		Nx_=Nx;
		}
	if (Ny_!=Ny) {
		cout << "TEMP FUDGE, setting Ny to input value"<<endl;
		Ny_=Ny;
		}
	if (Nz_!=Nz) {
		cout << "TEMP FUDGE, setting Nz to input value"<<endl;
		Nz_=Nz;
		}
	
	su_.SetEmissionRedShift(zref_);
	double cod = su_.RadialCoordinateMpc(); // radial distance
	cout <<"    Redshift of central pixel = "<< zref_;
	cout <<", comoving distance @zref = "<< cod <<endl;
	
	// INITIALISE ARRAYS TO GRID SIZE
	int ndim=3;
	sa_size_t mydim[ndim];
	mydim[0]=Nx_; mydim[1]=Ny_; mydim[2]=Nz_;
	ngals_.SetSize(ndim, mydim);	// number of galaxies in each cell: data
	wngals_.SetSize(ndim, mydim);	// weighted galaxy density field: data (same as ngal if no SF)
	if (sfcompute_) {
		randomcat_.SetSize(ndim, mydim);// number of galaxies in each cell: random catalog
		weights_.SetSize(ndim, mydim);	// weights of galaxy density field (per cell)
		wrgals_.SetSize(ndim, mydim);	// weighted galaxy density field: random catalog
		}
		
	cout <<"    EXIT Cat2Grid::SetGrid()"<<endl<<endl<<endl;
};


void Cat2Grid::SaveSelecFunc(string SFTextFile, string FullCat, string ZCol)
// Create Histo of observed redshifts/true redshifts and save to
// a text file called [SFTextFile]_nofz.txt
// Reads in a Fits file containing ALL the true redshifts in the
// simulation, loops over these redshifts and adds them to Histo
// object
// Loops over observed catalog and adds observed redshifts to
// Histo object
{
	cout <<"    Cat2Grid::SaveSelecFunc()"<<endl;
	Timer tm("SaveSelecFunc");
	
	// Get full catalog file name(s)
	string delim=",";
    vector<string> fnames;
    stringSplit(FullCat,delim,fnames);
    int Nfiles = fnames.size();
	
	/*// Set name of catalog containing the "full" set of z's
	string FullCatFile;
	if(Nfiles>1)
		cout <<"    Reading in "<< Nfiles <<" full catalog files"<<endl; 
	else
		FullCatFile = FullCat;*/
	
	// Initialise random generator to specified seed
	long seed=1;
	rg_.SetSeed(seed);
	
	// Set up Histogram bins
	double minz = 0, maxz=10;
	int nbin=5000;
	Histo nzFC(minz, maxz, nbin); // full catalog Histo
	Histo nzOC(minz, maxz, nbin); // obs catalog Histo
	// (will be same as above if spec-z used in analysis)
	Histo nzOCsz(minz, maxz, nbin); // obs catalog Histo (using spec-z) 
	
	// Read in Fits file containing column with all true redshifts
	// FULL CATALOG		
	double minzt, maxzt;
	
	for (int ic=0; ic<Nfiles; ic++) {
	
			string FullCatFile = fnames[ic];
			cout <<"    Read in full catalog redshifts from "<< FullCatFile;
			cout <<" from column labelled "<< ZCol <<endl;
			FitsInOutFile fin(FullCatFile,FitsInOutFile::Fits_RO);
			fin.MoveAbsToHDU(2);
			SwFitsDataTable dt(fin,512,false);
			cout <<endl;
			
			sa_size_t ncat = dt.NEntry();
			sa_size_t Izs = dt.IndexNom(ZCol);
			dt.GetMinMax(Izs,minzt,maxzt); 
			DataTableRow rowin = dt.EmptyRow();
			cout <<"    Number of galaxies in this FULL catalog = "<< ncat <<endl;
			cout <<"    Min z of this FULL catalog = "<<minzt;
			cout <<", max z of this FULL catalog = "<<maxzt<<endl;
			
			cout <<"    Add to Histogram ... "<<endl;
			for(sa_size_t i=0; i<ncat; i++) { 
				dt.GetRow(i, rowin);
				double zs=rowin[Izs];
				nzFC.Add(zs);
				}
			}

		
	sa_size_t ngFC = nzFC.NEntries();
	 
	// OBS CATALOG - loop over and add to Histo
	cout <<"    Histogram up OBSERVED catalog"<<endl;
	DataTableRow rowino = dt_.EmptyRow();
	GalRecord grec;
	double x=1e8,y=1e8,z=-1e8,redshift=-10;
	for(long ig=0; ig<ngall_; ig++) {

		dt_.GetRow(ig, rowino);
		Row2Record(rowino,grec);
		if (!Filter(grec)) continue;
		  
		Rec2EuclidCoord(grec,x,y,z,redshift);
		
		nzOC.Add(redshift); // histo up OBSERVED -z
		nzOCsz.Add(grec.zs);// histogram up spec-z no matter what
		
		
		// print out first 10 gals
        if(ig<10) {
		    cout <<"    galid="<< ig <<": ";
		    grec.Print();
		    if (AddGaussErr_)
			    cout <<", zpG="<< redshift <<endl;
		    else
			    cout <<endl;
		    }
		 }
	
	if (DoDebug_) {
		
		string outfile3 = SFTextFile+"_histo.ppf";
		// Let's keep the three histograms 
		POutPersist poh(outfile3);
		poh << PPFNameTag("nzFC") << nzFC;
		poh << PPFNameTag("nzOC") << nzOC;
		poh << PPFNameTag("nzOCsz") << nzOCsz;
		cout << "    Cat2Grid::SaveSelecFunc/Info nzFC,nzOC,nzOCsz";
		cout << " histos saved to file " << outfile3 << endl;
		}
	
	
	//READ HISTOGRAM VALUES INTO A FILE
	string outfile = SFTextFile+"_nofz.txt";
	string outfile2 = SFTextFile+"_specz_nofz.txt";
	
	cout <<"    Write n^o(z)/n^t(z) histograms to text files "<<endl;
	ifstream inp;
	ofstream outp,outp2;
	
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(ios::failbit);
		cout << "    Writing to file ..." << outfile.c_str() << " and ";
		cout << outfile2 <<endl;
		outp.open(outfile.c_str(), ofstream::out);
		outp2.open(outfile2.c_str(), ofstream::out);
		
		for(int_4 i=0;i<nbin;i++) {
			r_8 bc=nzFC.BinCenter(i);
			r_8 bc2=nzFC.BinCenter(i);
			r_8 nzt=nzFC.operator()(i);
			r_8 nzo=nzOC.operator()(i);
			r_8 nzosz=nzOCsz.operator()(i);
		
			if((bc-bc2)>0.000002)
				cout <<"DIFFERENCE BETWEEN BIN CENTERS = "<< bc-bc2 <<endl;
						
			r_8 sf1=nzo/nzt; // phot-z (if ztype is phot)
			r_8 sf2=nzosz/nzt;// spec-z
			
			if (nzo<1&&nzt<1)// ie if they are both 0
				sf1=1;
			if (nzo<1&&nzt>0)// if obs is 0
				sf1=0.000001; //vsmall not zero, otherwise weight=1/SF=INF
			if (nzo>0&&nzt<1)// if true is 0 (can get this with photo-z scattering out of sim volume)
				sf1=50; //vlarge, i.e. gal downweighted to zero weight
				
			if (nzosz<1&&nzt<1)// ie if they are both 0
				sf2=1;
			if (nzosz<1&&nzt>0)// if obs is 0
				sf2=0.000001; //vsmall not zero, otherwise weight=1/SF=INF
			if (nzosz>0&&nzt<1)// if true is 0 (should not be possible if using spec-z as ztype)
				sf2=50; //vlarge, i.e. gal downweighted to zero weight
				
//			// if bin is outside full catalog's redshift range
//			if (bc<minzt||bc>maxzt)
//				{
//				sf1=1000000;//vlarge, anthing here is downweighted to zero weight
//				sf2=1000000;
//				}
				
			outp << bc <<"      "<< sf1 <<endl;
			outp2 << bc <<"      "<< sf2 <<endl;
			}
		outp.close();
		outp2.close();
		}
	else
		cout << "Error...file """ << outfile.c_str() << """ exists" << endl;


	tm.Split();
	cout <<"    Elapsed time "<< tm.TotalElapsedTime() <<endl;
	cout <<"    EXIT Cat2Grid::SaveSelecFunc()"<<endl<<endl;
};


void Cat2Grid::GalGrid(double SkyArea)
// Lays grid over data and counts the number of galaxies in each cell
// Cartesian z-direction will be radial direction from observer to center 
// of survey area
{
	cout <<endl<<"    Cat2Grid::GalGrid()"<<endl; 
	Timer tm("GalGrid");
	// Initialise random generator so specified seed
	long seed=1;
	rg_.SetSeed(seed);

	// Set sky area
	SkyArea_ = SkyArea;
	
	// Check grid has been set
	if (!ngals_.IsAllocated())
	throw ParmError("ERROR! Grid has not been set");
	
	Npix_ = Nx_*Ny_*Nz_;
	volgrid_ = Npix_*pow(cellsize_,3);
	cout <<"    Volume of grid="<< volgrid_ <<", Volume of survey="<< Vol_ <<endl;

	// Find total number of observed pixels
	cout <<"    Finding number of observed pixels ... "<<endl;
	n_obs_pixels_ = ObsPixels();
	cout <<"    ... "<< n_obs_pixels_ <<" are observed out of "<< Npix_ <<" total pixels"<<endl;
		
	// NUMBER OF GALAXIES IN EACH CELL
	// first set array to zero
	ngals_=0.;

	cout <<"    Cells in each direction: Nx="<< Nx_ <<", Ny="<< Ny_ <<", Nz=";
	cout << Nz_ <<endl;
	cout <<"    Grid boundaries: "<< Xmin_ <<"<X<"<< Xmax_ <<", ";
	cout << Ymin_ <<"<Y<"<< Ymax_ <<", "<< Zmin_ <<"<Z<"<< Zmax_ <<endl;
	cout <<"    Grid cell size = "<< cellsize_ <<endl;
	cout <<"    Start loop over galaxies..."<<endl;
	cout <<"    applying filters ..."<<endl;
	cout <<"    Reading in OBSERVED redshifts from column "<< ZOCol_ <<endl;
	if (AddGaussErr_)
		cout <<"    Adding Gaussian photo-z error of size "<< PZDerr_ <<" Mpc to these"<<endl;
	
	
	// Take each galaxy and find which cell it lies in
	DataTableRow rowin = dt_.EmptyRow();
	GalRecord grec;
	double x=1e8,y=1e8,z=-1e8,redshift=-10;
	double selphi=1.;
	for(long ig=0; ig<ngall_; ig++) {
		
		// get row values from data table
	  	dt_.GetRow(ig, rowin); 
	  	
	  	// add row values to GalRecord class
	  	Row2Record(rowin,grec);
	  	
	  	// add possible filtering on e.g. magnitude here
	  	// no filtering is implemented yet
	  	if (!Filter(grec)) continue; 
	  
	    // count galaxy as observed 
	  	ngo_++; 
	  	
	  	// convert galaxy position into Euclid coord
	  	Rec2EuclidCoord(grec,x,y,z,redshift);
	  	
	  	// return selection function value at gal redshift
      	selphi = (*selfuncp_)(redshift);
      	
      	// basic check that selection function value isn't crazy
	  	if (selphi<0)
			cout <<"phi<0, phi="<< selphi <<", z="<< redshift <<endl;
	  	
		// print out first 10 gals, just as a basic check
	  	if(ig<10) {
			cout <<"    galid="<< ig <<": ";
			grec.Print();
			if (AddGaussErr_)
				cout <<", zpG="<< redshift <<", phi="<< selphi <<endl;
			else
				cout <<endl;
			}

        // add galaxy to correct grid pixel 
		// ng_ and ngout_ are summed up in here:
	  	AddToCell(x,y,z,selphi); 
	  
      	}

	cout <<"    Number of galaxies actually observed = "<< ngo_;
	cout <<", TOTAL number of galaxies in catalog = "<< ngall_ <<endl;
	cout <<"    Number of galaxies inside grid = "<< ng_ <<endl;
		
	cout <<"    Number of 'galaxies' inside WEIGHTED grid = "<< ngw_ <<endl;
	cout <<"    .... (should be about equal to the number of galaxies in FULL";
	cout <<" catalog - not printed here)"<<endl;
	
	cout <<"    Number of galaxies outside grid = "<< ngout_ <<endl;
	
	if (ngo_>ngall_)
		throw ParmError("ERROR! galaxy number discrepancy: ngo_!=ngall_");
	if (ng_!=ngals_.Sum())
		throw ParmError("ERROR! galaxy number discrepancy: ng_!=ngals_.Sum()");
	if (ngo_!=ngals_.Sum())
		cout <<"    Number of galaxies outside grid = "<< ngout_ <<endl;

	// normalise the arrays
	NormNArrays();

	// write the arrays here!
	WriteGalArrays();

	tm.Split();
	cout <<"    Elapsed time "<< tm.TotalElapsedTime() <<endl;
	cout <<"    EXIT Cat2Grid::GalGrid()"<<endl<<endl; 

};


sa_size_t Cat2Grid::ObsPixels()
{

	sa_size_t nobspix = 0;
	for(int i=0; i<Nx_;i++)
		for(int j=0; j<Ny_;j++)
			for(int k=0; k<Nz_;k++) {
				
				double xc,yc,zc;
				GetCellCoord(i,j,k,xc,yc,zc);
				double dcell= sqrt(xc*xc+yc*yc+zc*zc);
				double phic = acos(zc/dcell);
				if (phic<=SkyArea_)
					nobspix++;
				}
	return nobspix;
};


void Cat2Grid::Row2Record(DataTableRow& rowin, GalRecord& rec)
{

	if (RadialZ_) {
		rec.xcoo = rowin[Ic1_]; rec.ycoo= rowin[Ic2_];
		rec.alpha = -9999;   rec.delta=-9999;
		rec.zs = rowin[Izs_];  rec.zo = rowin[Iz_];
		}
	else {
		rec.xcoo = -9999 ; rec.ycoo= -9999;
		rec.alpha = rowin[Ic1_];   rec.delta=rowin[Ic2_];
		rec.zs = rowin[Izs_];  rec.zo = rowin[Iz_];
		}

	// rec.zs = SPECTRO Z
	// rec.zo = OBSERVED Z (could be spec-z, gauss-z, photo-z)
	
	// Fill in the magnitudes if available 
};


bool Cat2Grid::Filter(GalRecord& rec)
// Cut on magnitudes if applicable
{
	return true;
};


void Cat2Grid::Rec2EuclidCoord(GalRecord& rec, double& x, double& y, double& z, double& redshift)
// From GalRecord object which contains:
// double alpha,delta,glong,glat,xcoo,ycoo,zcoo;
// double zs,zo;
// int type;
// double u,g,r,i,z,y;
// Calculates: x,y,z comoving coordinates, and returns which redshift is being analysed
// Choice of 3 redshifts: photometric, spectroscopic and spectroscopic + Gaussian error
{

	redshift = rec.zo;// zo is redshift to be used in analysis (could be spec-z,phot-z,gauss-z)
	double dc = z2dist_(redshift);

	if (RadialZ_) {
		x = rec.xcoo;
		y = rec.ycoo;		
		z = dc;
		}
	else {
		double ph=rec.delta;
		double th=rec.alpha;
		if(isnan(ph)) {
			cout <<"    Cat2Grid::Rec2EuclidCoord/PB phi=nan -> set phi=0"<<endl;
			ph=0;
			}

		x=dc*cos(th)*sin(ph);
		y=dc*sin(th)*sin(ph);
		z=dc*cos(ph);
		}
		
	if (AddGaussErr_) {
		z += PZDerr_*rg_.Gaussian();
		// We have to compute back the redshift 
		redshift = dist2z_(sqrt(x*x+y*y+z*z));
		}
};


void Cat2Grid::AddToCell(double x, double y, double z,double phi)
// Find which index galaxy with (x,y,z) coordinates lies in
// Fill ngals and wngals arrays accordingly
{

	sa_size_t indexx=(sa_size_t)floor((x-Xmin_)/cellsize_);
	sa_size_t indexy=(sa_size_t)floor((y-Ymin_)/cellsize_);
	sa_size_t indexz=(sa_size_t)floor((z-Zmin_)/cellsize_);
	if (indexx<Nx_&&indexy<Ny_&&indexz<Nz_&&indexx>=0&&indexy>=0&&indexz>=0) {
	
		ngals_(indexx,indexy,indexz)++;
		if (sfcompute_) { // NOTE weight does not have to be 1/SF
			wngals_(indexx,indexy,indexz)+=1/phi;
			ngw_+=1/phi; }
		else		// AND if cat has no SF, wngals_=ngals_
			wngals_(indexx,indexy,indexz)++;
		ng_++;
		}
	else
		ngout_++;
};


void Cat2Grid::GetCellCoord(sa_size_t i, sa_size_t j, sa_size_t k, double& x, double& y, double& z)
// Given grid pixel index (i,j,k) return comoving cartesian coordinate of pixel center (x,y,z)
{

	x=GetCellX(i);
	y=GetCellY(j);
	z=GetCellZ(k);
	return;
};


void Cat2Grid::RandomGrid(double nc, bool SaveArr)
// Compute random weighted grid with same selection function as data
{
	cout <<endl<<"    Cat2Grid::RandomGrid()"<<endl;
	Timer tm("RandomGrid");
	
	cout <<"    Random catalog mean density = "<< nc <<endl;
	if (!sfcompute_)
		cout <<"    Selection function should be CONSTANT with z, check this ..."<<endl;
	else
		cout <<"    Check selection function .... "<<endl;
	cout <<"    Only filling pixels with phi <= "<< SkyArea_ <<endl;
	{
	int nz=5;
	double dz=(zsmax_-zsmin_)/(nz-1);
	for (int i=0;i<nz;i++)
		cout <<"    z="<< zsmin_+dz*i <<", phi="<< (*selfuncp_)(zsmin_+dz*i) <<endl;
	}
		
	//	throw ParmError("ERROR! selection function has NOT been computed");
		
	if (wrgals_.NbDimensions()<2) {
		int ndim=3;
		sa_size_t mydim[ndim];
		mydim[0]=Nx_; mydim[1]=Ny_; mydim[2]=Nz_;
		randomcat_.SetSize(ndim, mydim);// number of galaxies in each cell: random catalog
		weights_.SetSize(ndim, mydim);	// weights of galaxy density field (per cell)
		wrgals_.SetSize(ndim, mydim);	// weighted galaxy density field: random catalog
		}
	if (SaveArr) {
		int ndim=3;
		sa_size_t mydim[ndim];
		mydim[0]=Nx_; mydim[1]=Ny_; mydim[2]=Nz_;
		zc_.SetSize(ndim, mydim); // redshifts of pixel centers
		}
	cout <<"    Compute weights and random catalog ..."<<endl;
	
	
	for(int i=0; i<Nx_; i++) {
	
		//tm.Split();
		//cout<<"outer loop = "<<i<<", time elapsed = "<<tm.TotalElapsedTime()<<endl;
		for(int j=0; j<Ny_; j++)
			for(int k=0; k<Nz_; k++)	{
			
				double xc,yc,zc;
				GetCellCoord(i,j,k,xc,yc,zc);
				double dcell= sqrt(xc*xc+yc*yc+zc*zc);
				double phic = acos(zc/dcell);
				double redshift = dist2z_(dcell);
				double phi = (*selfuncp_)(redshift);
				
				if (phic<=SkyArea_)	{
					
					weights_(i,j,k) = 1/phi; // NOTE weight does not have to be 1/SF
					// Average number of gals expected in cell
					double mu = phi*nc; // just Poisson phi
					uint_8 npoiss = rg_.PoissonAhrens(mu); // Poisson fluctuate
					randomcat_(i,j,k) = (double)npoiss;
				
					wrgals_(i,j,k)=weights_(i,j,k)*randomcat_(i,j,k); // deleted alpha
					}
				else {
					randomcat_(i,j,k) = 0;
					weights_(i,j,k) = 0;
					wrgals_(i,j,k) = 0;
					}

				if (SaveArr)
					zc_(i,j,k) = redshift;
				
				}
		}
	nrand_=randomcat_.Sum(); 
	wnrand_=wrgals_.Sum();
	cout <<"    number of gals in in grid = "<< ng_;
	cout <<", number of gals in random grid = "<< nrand_ <<endl;
	cout <<"    Number of galaxies in RANDOM weighted grid = "<< wnrand_;
	cout <<", should be approximately input density*npixels = "<< nc*n_obs_pixels_ <<endl;
	cout <<"    Calculated weighted random grid density = "<< wnrand_/n_obs_pixels_;
	cout <<", input density = "<< nc <<" (should be approximately the same)"<<endl;
	
	// Ratio of number of observed WEIGHTED gals to WEIGHTED gals in random catalog
	if (ngw_>0)
		alph_=(double)(ngw_/wnrand_); 
	else
		alph_=(double)(ng_/wnrand_); 
	cout <<"    alpha=ngw/nsw="<< alph_ <<endl; // alpha should be ntrue/nsyn
	VarianceRandomGrid();
	cout <<"    Mean of weighted random grid should be roughly "<< nc;
	cout <<" if all pixels are observed, "<< nc*n_obs_pixels_/Npix_ <<" if not"<<endl;
	double minv,maxv;
	randomcat_.MinMax(minv,maxv);
	cout <<"    Minimum value in random cat grid = "<< minv <<endl;
	cout <<"    Maximum value in random cat grid = "<< maxv <<endl;
	weights_.MinMax(minv,maxv);
	cout <<"    Minimum value in weights grid = "<< minv <<endl;
	cout <<"    Maximum value in weights grid = "<< maxv <<endl;

	// normalise the array
	NormRArray();

	// write the arrays here!
	fos_ << wrgals_;
	fos_ << zc_;

	tm.Split();
	cout <<"    Elapsed time "<< tm.TotalElapsedTime()<<endl;
	cout <<"    EXIT Cat2Grid::RandomGrid()"<<endl<<endl; 

};


vector<double> Cat2Grid::XtractSubArray(TArray<r_8>& nsub, sa_size_t dNp, 
                                       sa_size_t Np, double theta, int arrayflag)
// Np	 : the starting pixel number of sub-array
// dNp	 : width of sub-array in z-direction
// theta : MAX angle grid should cover, if want whole grid just make sure theta
//		   is any value greater than angular extent of whole grid
// arrayflag : specifies which array: ngals,wngals,wrgals to return in nsub
{
	cout <<"    Cat2Grid::XtractSubArray() Z-axis split version"<<endl;
	
	cout <<"    Sub-array beginning at array pixel index "<< Np <<endl;
	sa_size_t Nx = GetNTransPix(Np,theta);
	
	// min and max pixel indices
	sa_size_t minxi = idx_-Nx;
	sa_size_t minyi = idy_-Nx;
	sa_size_t maxxi = idx_+Nx;
	sa_size_t maxyi = idy_+Nx;
	sa_size_t minzi = Np;
	sa_size_t maxzi = Np+dNp-1;
	
	// check min,max pixel indices
	if (maxxi>ngals_.SizeX()-1) {
		cout <<"    maxx="<< maxxi <<", so setting maxx="<< ngals_.SizeX()-1 <<endl;
		maxxi = ngals_.SizeX()-1;
		}
	if (maxyi>ngals_.SizeY()-1) {
		cout <<"    maxy="<< maxyi <<", so setting maxy="<< ngals_.SizeY()-1 <<endl;
		maxyi = ngals_.SizeY()-1;
		}
	if (maxzi>ngals_.SizeZ()-1) {
		cout <<"    maxz="<< maxzi <<", so setting maxz="<< ngals_.SizeZ()-1 <<endl;
		maxzi = ngals_.SizeZ()-1;
		}
	if (minxi<0) {
		cout <<"    minx="<< minxi <<", so setting minx=0"<<endl;
		minxi = 0;
		}
	if (minyi<0) {
		cout <<"    miny="<< minyi <<", so setting miny=0"<<endl;
		minyi=0;
		}
	if (minzi<0) {
		cout <<"    minz="<< minzi <<", so setting minz=0"<<endl;
		minzi=0;
		}
		
	if (minxi<0||minyi<0||minzi<0)
		throw ParmError("ERROR! either minx<0||miny<0||minz<0");
	if (maxxi>ngals_.SizeX()-1||maxyi>ngals_.SizeY()-1||maxzi>ngals_.SizeZ()-1)
		throw ParmError("ERROR! either maxx>ngals_.SizeX()-1||maxy>ngals_.SizeY()-1||maxz>ngals_.SizeZ()-1");

	cout <<"    Ranges of sub array: (xmin,xmax) = ("<< minxi <<","<< maxxi <<"),";
	cout <<" (ymin,ymax) = ("<< minyi <<","<< maxyi <<"), (zmin,zmax) = (";
	cout << minzi <<","<< maxzi <<")"<<endl;
	
	if(arrayflag<0)
		nsub = ngals_(Range(minxi,maxxi), Range(minyi,maxyi), Range(minzi,maxzi)).PackElements();
	else if (arrayflag==0)
		nsub = wngals_(Range(minxi,maxxi), Range(minyi,maxyi), Range(minzi,maxzi)).PackElements();
	else if (arrayflag>0)
		nsub = wrgals_(Range(minxi,maxxi), Range(minyi,maxyi), Range(minzi,maxzi)).PackElements();
		
	// need to return min and max redshifts of this sub grid
	// min redshift is redshift corresponding to: Zmin_+Np*cellsize_;
	// max redshift is redshift corresponding to:  = Nx*cellsize_/sin(theta)
	double mind = Zmin_+Np*cellsize_;
	double maxd = mind + (maxzi-minzi)*cellsize_;//(Nx*cellsize_)/sin(theta);
	double centerd = (mind + maxd)/2;
	
	// We have to initialize the distance to redshift conversion interpolator
	int_8 nz=1000;
	vector<double> zrs, codist;
	double minz=0, maxz=10;
	double dz = (maxz-minz)/(nz-1);
	for(int kk=0; kk<nz; kk++) 
		{
		  double zs=minz+kk*dz;
		  su_.SetEmissionRedShift(zs);
		  double cod =su_.RadialCoordinateMpc(); // radial distance 
		  zrs.push_back(zs);
		  codist.push_back(cod); 
		}
	SInterp1D dist2z(codist,zrs,codist[0],codist[codist.size()-1],2*nz);
	
	double minza = dist2z(mind);
	double maxza = dist2z(maxd);
	double cenza = dist2z(centerd);
	vector<double> zedge;
	zedge.push_back(minza);
	zedge.push_back(cenza);
	zedge.push_back(maxza);
	cout << "    min,center,max redshifts of array = "<< minza <<","<< cenza;
	cout <<","<< maxza <<endl;
	
	cout <<"    EXIT Cat2Grid::XtractSubArray()"<<endl<<endl;

	return zedge;
};


vector<double> Cat2Grid::XtractSubArray(TArray<r_8>& nsub, long x1, long x2,
    long y1, long y2, long z1, long z2, int arrayflag)
// Extract sub array uses pixel ranges in each dimension
{

	cout <<"    Cat2Grid::XtractSubArray() Sub-array range version"<<endl;
	
	if(arrayflag<0)
		nsub = ngals_(Range(x1,x2), Range(y1,y2), Range(z1,z2)).PackElements();
	else if (arrayflag==0)
		nsub = wngals_(Range(x1,x2), Range(y1,y2), Range(z1,z2)).PackElements();
	else if (arrayflag>0)
		nsub = wrgals_(Range(x1,x2), Range(y1,y2), Range(z1,z2)).PackElements();
	
	// need to return min and max redshifts of this sub grid
	// min redshift is redshift corresponds to distance to first pixel in z-dim
	// max redshift is redshift corresponds to distance to last pixel in z-dim
	double mind = Zmin_+z1*cellsize_;
	double maxd = mind + (z2-z1)*cellsize_;
	double centerd = (mind + maxd)/2;
	
	// We have to initialize the distance to redshift conversion interpolator
	int_8 nz=1000;
	vector<double> zrs, codist;
	double minz=0, maxz=10;
	double dz = (maxz-minz)/(nz-1);
	for(int kk=0; kk<nz; kk++) 
		{
		  double zs=minz+kk*dz;
		  su_.SetEmissionRedShift(zs);
		  double cod =su_.RadialCoordinateMpc(); // radial distance 
		  zrs.push_back(zs);
		  codist.push_back(cod); 
		}
	SInterp1D dist2z(codist,zrs,codist[0],codist[codist.size()-1],2*nz);
	
	double minza = dist2z(mind);
	double maxza = dist2z(maxd);
	double cenza = dist2z(centerd);
	vector<double> zedge;
	zedge.push_back(minza);
	zedge.push_back(cenza);
	zedge.push_back(maxza);
	cout << "    min,center,max redshifts of array = "<< minza <<","<< cenza;
	cout <<","<< maxza <<endl;
	
	cout <<"    EXIT Cat2Grid::XtractSubArray()"<<endl<<endl;

	return zedge;
};


sa_size_t Cat2Grid::GetNTransPix(sa_size_t Np, double theta)
// return number of pixels in x or y dimension that cover
// an angle of size theta at distance Dp
{

	double Dp = Zmin_ + Np*cellsize_; // min distance to center pixel on x-y plane
	double xw = Dp*tan(theta)-cellsize_/2;
	sa_size_t Nx=(sa_size_t)floor(xw/cellsize_);
	
	cout <<"    Distance to sub array = "<< Dp <<" Mpc"<<endl;
	cout <<"    Number of pixels which cover "<< theta <<" radians";
	cout <<" at distance "<< Dp <<" Mpc, = "<< Nx <<endl;
	
	return Nx;
};


void Cat2Grid::OutputEuclidCat(double SkyArea)
// Outputs catalog in Euclidean coords
{
	cout <<endl<<"    Cat2Grid::OutputEuclidCat()"<<endl; 

	string outfile = debugoutroot_+"_euclidcat.fits";
	// We create first a datatable with the fits file as the swap space 
	FitsInOutFile swf(outfile, FitsInOutFile::Fits_Create);	
	SwFitsDataTable gals(swf, 2048);
	//gals.AddLongColumn("GalID");
	gals.AddFloatColumn("theta");
	gals.AddFloatColumn("phi");
	gals.AddFloatColumn("z");
	gals.AddFloatColumn("rg");
	gals.AddFloatColumn("xg");
	gals.AddFloatColumn("yg");
	gals.AddFloatColumn("zg");
	DataTableRow row = gals.EmptyRow();

	// Set sky area
	SkyArea_ = SkyArea;
	
	Npix_= Nx_*Ny_*Nz_;
	volgrid_= Npix_*pow(cellsize_,3);
	cout <<"    Volume of grid="<< volgrid_ <<", Volume of survey="<< Vol_ <<endl;

	cout <<"    Start loop over galaxies..."<<endl;
	cout <<"    applying filters ..."<<endl;
	cout <<"    Reading in OBSERVED redshifts from column "<< ZOCol_ <<endl;
	if (AddGaussErr_)
		cout <<"    Adding Gaussian photo-z error of size "<< PZDerr_ <<" Mpc to these"<<endl;
	
	
	// Set counting galaxies to ZERO
	ngo_=0; // num gals OBSERVED

	// Take each galaxy and find which cell it lies in
	DataTableRow rowin = dt_.EmptyRow();
	GalRecord grec;
	double x=1e8,y=1e8,z=-1e8,redshift=-10;
	double selphi=1.;
	for(long ig=0; ig<ngall_; ig++) {
		
	  	dt_.GetRow(ig, rowin); // get row values from data table
	  	Row2Record(rowin,grec);// add row values to GalRecord class
	  	if (!Filter(grec)) continue; // add possible filtering on e.g. magnitude here
	  
	  	ngo_++; // count galaxy
	  	Rec2EuclidCoord(grec,x,y,z,redshift);// convert galaxy position into Euclid coord
			
		// print out first 10 gals
	  	if(ig<10) {
			
			cout <<"    galid="<< ig <<": ";
			grec.Print();
			if (AddGaussErr_)
				cout <<", zpG="<< redshift <<", phi="<< selphi <<endl;
			else
				cout <<endl;
			}

		//su_.SetEmissionRedShift(redshift);
		double dc = z2dist_(redshift);//su_.RadialCoordinateMpc();
		double ph = grec.delta;
		double th = grec.alpha;
		uint_8 gid = rowin[Iid_];

		// write to file
		/*row[0] = gid;
		row[1] = th;
		row[2] = ph;
		row[3] = redshift; 
		row[4] = dc; 
		row[5] = x; 
		row[6] = y; 
		row[7] = z;*/
		row[0] = th;
		row[1] = ph;
		row[2] = redshift; 
		row[3] = dc; 
		row[4] = x; 
		row[5] = y; 
		row[6] = z;
		gals.AddRow(row);

      	}

	cout <<"    Number of galaxies actually observed = "<< ngo_ <<endl;
	cout <<"    EXIT Cat2Grid::OutputEuclidCat()"<<endl<<endl; 

};

//double Cat2Grid::AddSurveyWindow(double Phi)
//// Sets fluctuation in all cells with an angle outside of the survey area equal to zero
//// galdens == delta == ngal/nbar - 1
//// Therefore if there are ZERO galaxies in the cell the cell will 
//// have value -1
//{
//	
//	int_8 nzero=0;
//	for(sa_size_t i=0; i<ngals_.SizeX();i++)
//		for(sa_size_t j=0; j<ngals_.SizeY();j++)
//			for(sa_size_t k=0; k<ngals_.SizeZ();k++)
//				{
//				
//				// cartesian coord of cell
//				double xc,yc,zc,phi;
//				GetCellCoord(i,j,k,xc,yc,zc);
//				phi=GetPhi(xc,yc,zc);
//				
//				if (phi>Phi)
//					{galdens_(i,j,k)=-1; nzero++;}
//				}
//	double nfrac = (double)nzero/(ngals_.SizeX()*ngals_.SizeY()*ngals_.SizeZ());
//	return nfrac;
//}
//
//
//
// junk these two
//void Cat2Grid::ApplyVaryingSF()
//{
//
//	cout <<"    Cat2Grid::ApplyVaryingSF()"<<endl;
//	
//	if (!sfcompute_)
//		throw ParmError("ERROR! selection function NOT set");
//		
//	// interpolate selection function already calculated 
//	int nint=1000;
//	SInterp1D interpSF(dc_,phi_,dc_[0],dc_[dc_.size()-1],nint);
//	
//	for(int i=0; i<Nx_;i++)
//		for(int j=0; j<Ny_;j++)
//			for(int k=0; k<Nz_;k++)
//				{
//				double dcgal=dgals_(i,j,k);
//				double phi=interpSF(dcgal);
//				//cout << "dc="<<dcgal<<", phi="<<phi<<endl;
//				ngals_(i,j,k)*=phi; // NOTE weight does not have to be 1/SF
//				if (phi==0&dcgal<dcmax_)
//					cout <<"phi=0, but dc="<<dcgal<<endl;
//				}
//				
//	cout <<"    Number of galaxies in grid AFTER selection function = "<<ngals_.Sum()<<endl;
//	cout <<"    EXIT Cat2Grid::ApplyVaryingSF()"<<endl<<endl;
//}
//
//void Cat2Grid::ApplyVaryingSF2wngals()
//{
//	cout <<"    Cat2Grid::ApplyVaryingSF2wngals()"<<endl;
//	
//	if (!sfcompute_)
//		throw ParmError("ERROR! selection function NOT set");
//		
//	// interpolate selection function already calculated 
//	int nint=1000;
//	SInterp1D interpSF(dc_,phi_,dc_[0],dc_[dc_.size()-1],nint);
//	
//	for(int i=0; i<Nx_;i++)
//		for(int j=0; j<Ny_;j++)
//			for(int k=0; k<Nz_;k++)
//				{
//				double dcgal=dgals_(i,j,k);
//				double phi=interpSF(dcgal);
//				wngals_(i,j,k)*=phi; // NOTE weight does not have to be 1/SF
//				if (phi==0&dcgal<dcmax_)
//					cout <<"phi=0, but dc="<<dcgal<<endl;
//				}
//				
//	cout <<"    Number of WEIGHTED galaxies in grid AFTER selection function = "<<wngals_.Sum()<<endl;
//	cout <<"    EXIT Cat2Grid::ApplyVaryingSF2wngals()"<<endl<<endl;
//}







// junk these too
//void Cat2Grid::ConvertGalDist2Fluct(double Phi)
//// Converts the number of galaxies in each cell to a fluctuation (delta_n) in each cell
//// delta_n = n/nbar - 1
//// n = number of galaxies in the cell
//// nbar = mean number of galaxies in a cell = total number of galaxies / total number of cells
//{ 
//	cout <<endl<<"    Cat2Grid::ConvertGalDist2Fluct()"<<endl;
//
//	double mult,neff;
//	if (sfcompute_)
//		{
//		//cout <<"     fluctuating w(r)*(n_g(r)-alpha*ns(r)) & w(r)*n_g(r)"<<endl;
//	//	galdens_=Fgals_;
//	//	neff=Fgals_.Sum();
//	//	mult=(double)totncells/neff; // 1 / average number of galaxies per cell
//	//	cout <<"    average number of galaxies in a grid cell="<<1/mult<<endl;
//		
//		cout <<"    fluctuating w(r)*n_g(r) .... "<<endl;
//		galdens_=wngals_;
//		neff=wngals_.Sum(); // neff unw = neff w*sqrt(norm)
//		mult=(double)Npix_/neff;// 1 / average number of galaxies per cell
//		cout <<"    average number of galaxies in a grid cell="<<1/mult<<endl;
//		//if (normweight)
//		//		cout <<"    Weights are normalised"<<endl;
//		}
//	else 
//		{
//		cout <<"    fluctuating n_g(r)"<<endl;
//		galdens_=ngals_;
//		
//		sa_size_t Np=0;
//		double nsum=0;
//		for(sa_size_t i=0; i<ngals_.SizeX();i++)
//			for(sa_size_t j=0; j<ngals_.SizeY();j++)
//				for(sa_size_t k=0; k<ngals_.SizeZ();k++)
//					{
//					
//					// cartesian coord of cell
//					double xc,yc,zc,phi;
//					GetCellCoord(i,j,k,xc,yc,zc);
//					phi=GetPhi(xc,yc,zc);
//					
//					double ng=ngals_(i,j,k);
//					
//					if (phi<Phi)
//						{nsum+=ng; Np++;}
//					}
//		neff=nsum;
//		mult=(double)Np/neff; // 1 / average number of galaxies per cell
//		cout <<"    sum of galaxies within survey cone = "<<neff<<endl;
//		cout <<"    number of pixels within survey cone = "<<Np<<endl;
//		cout <<"    average number of galaxies in a grid cell="<<1/mult<<endl;
//		}
//
//	galdens_*=mult;
//	galdens_-=1;
//
//	// mult1=Np/(sumN*sqrtn), mult2=Np/sumN
//	// galdens1=w * ngals, galdens2=w/sqrtn * ngals
//	// galdens1*mult1=w * ngals *Np/(sumN*sqrtn)
//	// galdens2*mult2=w/sqrtn * ngals *Np/sumN
//	// SAME!
//
//	cout <<"    EXIT Cat2Grid::ConvertGalDist2Fluct()"<<endl<<endl;
//	return;
//}
//
//void Cat2Grid::ConvertGalDist2Fluct2(TArray<r_8>& gdens)
//{
//
//	cout <<"    fluctuating n_g(r): NO -1"<<endl;
//	gdens=ngals_;
//	double neff=gdens.Sum();
//	double mult=(double)Npix_/neff; // 1 / average number of galaxies per cell
//	cout <<"    sum of galaxies within survey cone = "<<neff<<endl;
//	cout <<"    number of pixels within survey cone = "<<Npix_<<endl;
//	cout <<"    average number of galaxies in a grid cell="<<1/mult<<endl;
//	gdens*=mult;
//}
//
//void Cat2Grid::ConvertRandDist2Fluct()
//// Converts the number of galaxies in each cell to a fluctuation (delta_n) in each cell
//// delta_n = n/nbar - 1
//// n = number of galaxies in the cell
//// nbar = mean number of galaxies in a cell = total number of galaxies / total number of cells
//{ 
//	cout <<"    Cat2Grid::ConvertRandDist2Fluct()"<<endl;
//
//	double neff=wrgals_.Sum();
//	double mult=(double)Npix_/neff; // 1 / average number of galaxies per cell
//	cout <<"    average number of galaxies in a grid cell="<<1/mult<<endl;
//	
//	rgaldens_=wrgals_;
//	rgaldens_*=mult;
//	rgaldens_-=1;
//	
//	//cout << rgaldens_;
//	
//	
//	cout <<"    EXIT Cat2Grid::ConvertRandDist2Fluct()"<<endl<<endl;
//	return;
//}

// SELECTION FUNCTION STUFF ======== REMOVE
//void Cat2Grid::SetConstantSF(double sf)
//{
//
//	cout <<endl<<"    Cat2Grid::SetConstantSF()"<<endl;
//	
//	int nbin=1000; // number of bins of SF
//	
//	if (sfcompute_)
//		cout << "ERROR! selection function already exists"<<endl;
//
//    double dz = (zsmax_-zsmin_)/(nbin-1);
//	for (int i=0;i<nbin;i++)
//		{
//		// find luminosity distance at bin center
//		r_8 zg=zsmin_+i*dz;
//		su_.SetEmissionRedShift(zg);
//		dL_.push_back(su_.LuminosityDistanceMpc());
//		dc_.push_back(su_.RadialCoordinateMpc());
//		zs_.push_back(zg); // also keep redshift incase it's useful
//		phi_.push_back(sf);
//		}
//		
//	dcmin_=dc_[0];
//	dcmax_=dc_[dc_.size()-1];
//	cout <<"    SF dc min="<<dcmin_<<", SF dc max="<<dcmax_<<endl;
//	cout <<"    SF dL min="<<dmin_<<", SF dL max="<<dmax_<<endl;
//		
//	sfcompute_=true;
//	cout <<"    EXIT Cat2Grid::SetConstantSF()"<<endl<<endl;
//}
//
//void Cat2Grid::ComputeSFfromFile(string nzfile)
//{
//	cout <<"    Cat2Grid::ComputeSFfromFile()"<<endl;
//	
//	if (sfcompute_)
//		cout << "ERROR! selection function already exists"<<endl;
//		
//	PInPersist pis(nzfile);
//	Histo nzFC;
//	pis >> nzFC;
//
//	int nbin = nzFC.NBins();
//	double min=nzFC.XMin();
//	double max=nzFC.XMax();
//	
//	sa_size_t ngOC,ngFC_;
//	ngFC_=nzFC.NEntries(); // num gals in full catalog
//	ngOC=ngall_; // TOTAL num gals in obs catalog - less galaxies might be in analysied catalog
//	
//	cout<<"    number of galaxies in bins for full catalog = "<<nzFC.NEntries()<<endl;
//	cout<<"    number of galaxies in full catalog = "<<ngFC_<<endl;
//	cout<<"    Histo bin width for full catalog = "<<nzFC.BinWidth()<<endl; 
//	if (nzFC.NOver()>0||nzFC.NUnder()>0)
//		cout <<"    WARNING! galaxies in full catalog outside bin range"<<endl;
//	
//	Histo nzOC(min, max, nbin);
//	for(sa_size_t i=0; i<ngOC; i++) 
//		{
//		// set to impossible values
//		double x=1e8,y=1e8,z=-1e8,zs=-10;
//		sa_size_t ng=0;
//		FilterCat(i,x,y,z,zs);
//		if(x<1e8&&y<1e8&&z>-1e8&&zs>-10)
//			nzOC.Add(zs);
//		}
//	cout<<"    number of galaxies in bins for photoz catalog = "<<nzOC.NEntries()<<endl;
//	cout<<"    TOTAL number of galaxies in photoz catalog = "<<ngOC<<endl;
//	cout<<"    Histo bin width for photoz catalog = "<<nzOC.BinWidth()<<endl; 
//	if (nzOC.NOver()>0||nzOC.NUnder()>0)
//		cout <<"    WARNING! galaxies in photoz catalog outside bin range"<<endl;
//		
//	// GET TRUE SELECTION FUNCTION
//	// n(dL)^obs/n(dL)^true 
//	vector<double> phitmp,ntrue,nobs;
//	for (int i=0;i<nbin;i++)
//		{
//		// find luminosity distance at bin center
//		r_8 zg=nzFC.BinCenter(i);
//		su_.SetEmissionRedShift(zg);
//		dL_.push_back(su_.LuminosityDistanceMpc());
//		dc_.push_back(su_.RadialCoordinateMpc());
//		zs_.push_back(zg); // also keep redshift incase it's useful
//		
//		double nzobs=nzOC.operator()(i);
//		double nztr=nzFC.operator()(i);
//		
//		ntrue.push_back(nztr);
//		nobs.push_back(nzobs);
//		
//		phitmp.push_back(nzobs/nztr);
//		if (nzobs<1&&nztr<1)
//			phitmp[i]=1;
//		if (nzobs<1&&nztr>0)
//			phitmp[i]=0.000001; //vsmall not zero
//			
//		if(isnan(phitmp[i]))
//			cout <<"phi=no/nt is NAN: no="<<nzobs<<", nt="<<nztr<<endl;
//		}
//	dmax_=dL_[dL_.size()-1];
//	dcmax_=dc_[dc_.size()-1];
//	cout <<"    dLmax="<<dmax_<<endl;
//	cout <<"    dcmax="<<dcmax_<<endl;
//	
//	// these extra steps are totally redundant
//	// NOW STORE SF INTO phi_
//	for (int i=0;i<nbin;i++)
//		{
//		phi_.push_back(phitmp[i]);//(1./phitmp(i))*bw;// MULTIPLY BY BIN WIDTH
//		if(isnan(phi_[i]))
//			cout <<"WARNING! NAN IN SELECTION FUNCTION: DEAL WITH THIS"<<endl;
//		}
//	
//	if (nbar_>0)
//			cout <<"WARNING: LF galaxy number density is defined, nbar_="<<nbar_<<endl;
//	nbar_=0;
//	// set this to true so SF is applied, and random cat is created in ConstructGrid
//	//randcompute_=true;
//	sfcompute_=true;
//	cout <<"    EXIT Cat2Grid::ComputeSFfromFile()"<<endl<<endl;
//}
//

//
//
//void Cat2Grid::ComputeSFfromLF(double phistar,double Mstar,double alpha,double mlim, double Mc)
//// Computes selection function from LF params (phistar,Mstar,alpha) and detetion params (Mlim,Mc)
//// sets sfcompute_ to true
//{
//	cout <<"    Cat2Grid::ComputeSFfromLF()"<<endl;
//	
//	int nbin=1000; // number of bins of SF
//	
//	if (sfcompute_)
//		cout << "ERROR! selection function already exists"<<endl;
//	if (ngFC_>0)
//			cout <<"WARNING: Full catalog galaxy number is defined n_full="<<ngFC_<<endl;
//
//	// LUMINOSITY FUNCTION
//	phistar_ = phistar; Mstar_ = Mstar; alpha_ = alpha;
//	mlim_ = mlim; Mc_=Mc;
//	cout <<"    Schechter function paras: Mstar="<<Mstar_<<", phistar="<<phistar_<<", alpha="<<alpha_<<endl;
//	cout <<"    Limiting magnitude="<<mlim_<<", Complete at absolute magnitude <="<<Mc_<<endl;
//	
//	Schechter sch(phistar_,Mstar_,alpha_);
//	double schmin=-24, schmax=-13;// units of "M-5log10h70"
//	int schnpt=10000;
//	nbar_=sch.Integrate(schmin,schmax,schnpt);
//	cout <<"    Number density="<<nbar_<<" /Mpc^3"<<endl;
//	ngFC_=0; // set other equal to zero so can calcuate mean number density per cell correctly
//	
//	// compute selection function: interpolate at distance of galaxy later
//	su_.SetEmissionRedShift(zsmax_);
//	double maxdL=su_.LuminosityDistanceMpc()+50;
//	su_.SetEmissionRedShift(zsmin_);
//	double mindL=su_.LuminosityDistanceMpc()-50;
//	
//	//cout <<"    Make selection function"<<endl;
//	// create look up table - interpolate value of phi_ at
//	// each galaxy location instead of separately integrating sf
//	SelectFunc sf(phistar_,Mstar_,alpha_);
//	
//	//cout <<"    Make table"<<endl;
//	
//	//dL_.SetSize(nbin);
//	//phi_.SetSize(nbin);
//	double ddL=(maxdL-mindL)/(nbin-1);
//	for (int i=0;i<nbin;i++)
//		{
//		double x=mindL+i*ddL;
//		dL_.push_back(x);
//		double x2=sf(dL_[i],mlim_,Mc_);
//		phi_.push_back(x2);
//		}
//	dmin_=dL_[0];
//	dmax_=dL_[dL_.size()-1];
//	
//	// make look up table to convert dL's to dc's
//	// SHOULD ADD INTERP FUNCTION HERE
//	double zmin=zsmin_-0.1, zmax=zsmax_+0.1;
//	if (zmin<0)
//		zmin=0;
//	int nz=10000;
//	double dz=(zmax-zmin)/(nz-1);
//	//TVector<r_4> ztab(nz);
//	//TVector<r_4> dLtab(nz);
//	vector<double> ztab,dLtab;
//	for(int i=0;i<nz;i++)
//		{
//		double x1=zmin+i*dz;
//		ztab.push_back(x1);
//		su_.SetEmissionRedShift(ztab[i]);
//		dLtab.push_back(su_.LuminosityDistanceMpc());
//		}
//		
//		
//		// find closest codtab to comdist_(i) 
//	int idmin;
//	//dc_.SetSize(nbin);
//	for(int i=0; i<nbin; i++)
//		{
//		double dist=dL_[i];
//		double tmpmin=1e10;
//		double diff;
//		for(int m=0;m<nz;m++)
//			{
//			diff=abs(dLtab[m]-dist);
//			if(diff<tmpmin){tmpmin=diff; idmin=m;}
//			}
//		if (diff<10)
//			cout <<"difference between dL's="<<diff<<endl;
//		su_.SetEmissionRedShift(ztab[idmin]);
//		dc_.push_back(su_.RadialCoordinateMpc());
//		}
//	dcmin_=dc_[0];
//	dcmax_=dc_[dc_.size()-1];
//	cout <<"    SF dc min="<<dcmin_<<", SF dc max="<<dcmax_<<endl;
//	cout <<"    SF dL min="<<dmin_<<", SF dL max="<<dmax_<<endl;
//	
//	//cout <<"dc:"<<endl;
////	dc_.Show();
////	cout <<"dL:"<<endl;
////	dL_.Show();
////	cout <<"phi:"<<endl;
////	phi_.Show();
////	
//	//cout <<"    dcmax="<<dcmax_<<endl;
////	string fname= debugoutroot_ + "_SF.txt";
////	ifstream inp;
////	ofstream outp;
////	inp.open(fname.c_str(), ifstream::in);
////	inp.close();
////	if(inp.fail())
////	{
////		inp.clear(ios::failbit);
////		cout << "    Writing to file (1)..." << fname.c_str() << endl<< endl<< endl;
////		outp.open(fname.c_str(), ofstream::out);
////		for(int_4 i=0;i<dL_.Size();i++)
////		{
////			r_8 bcg=dL_(i);
////			r_8 bvg=phi_(i);
////
////			outp << bcg<<"      "<<bvg<<endl;
////		}
////		outp.close();
////	}
////	else cout << "Error...file """ << fname.c_str() << """ exists" << endl;
//	
//
//	
//	// set this to true so SF is applied, and random cat is created in ConstructGrid
//	sfcompute_=true;
//	//randcompute_=true;
//	cout <<"    EXIT Cat2Grid::ComputeSFfromLF()"<<endl<<endl;
//}
//void Cat2Grid::SetVaryingSF()
//{
//
//	cout <<endl<<"    Cat2Grid::SetVaryingSF()"<<endl;
//	
//	int nbin=1000; // number of bins of SF
//	
//	if (sfcompute_)
//		cout << "ERROR! selection function already exists"<<endl;
//
//    double dz = (zsmax_-zsmin_)/(nbin-1);
//	for (int i=0;i<nbin;i++)
//		{
//		// find luminosity distance at bin center
//		r_8 zg=zsmin_+i*dz;
//		su_.SetEmissionRedShift(zg);
//		dL_.push_back(su_.LuminosityDistanceMpc());
//		dc_.push_back(su_.RadialCoordinateMpc());
//		zs_.push_back(zg); // also keep redshift incase it's useful
//		double sf=(1/exp(-1))*exp(-pow(dc_[i],5)/pow(zmin_,5));
//		phi_.push_back(sf);
//		}
//		
//	dcmin_=dc_[0];
//	dcmax_=dc_[dc_.size()-1];
//	cout <<"    SF dc min="<<dcmin_<<", SF dc max="<<dcmax_<<endl;
//	cout <<"    SF dL min="<<dmin_<<", SF dL max="<<dmax_<<endl;
//		
//	sfcompute_=true;
//	cout <<"    EXIT Cat2Grid::SetVaryingSF()"<<endl<<endl;
//
//}
