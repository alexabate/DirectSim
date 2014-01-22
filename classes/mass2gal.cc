#include "mass2gal.h"
#include "ctimer.h"

/* ----------------------------------------- 
   BAO simulations , LSST & BAORadio , 
   A. Abate  LAL , 2009 
-------------------------------------------*/

//******* Constructors *******************************************************//

Mass2Gal::Mass2Gal(TArray<r_8> drho,  SimpleUniverse& su, RandomGeneratorInterface& rg, int nbadplanes) 
: su_(su) , rg_(rg)
// Reads in SimLSS output which is a 3D cube of delta values
// NOTE: SimLSS outputs cube so:
// AXIS 1 = Z (RADIAL) AXIS
// AXIS 2 = Y AXIS
// AXIS 3 = X AXIS
// delta is overdensity: (rho-rho^bar) / rho^bar
// the way the FFT is performed in SimLSS means there are 1 or 2 extra planes in the FIRST dimension
// these are removed here
// su_ holds the cosmological parameters of the SimLSS simulation
// int nbadplanes should be the difference between NX and drho.SizeX()
{

	cout <<"    Mass2Gal constructor: reads in density cube"<<endl;
	mass_ = drho(Range(0, drho.SizeX()-nbadplanes-1), Range::all(), Range::all()).PackElements();
	cout << "    Mass2Gal::Mass2Gal(TArray<r_4> drho,  int nbadplanes=" << nbadplanes << ")" << endl;
  	mass_.Show();
	
	cout <<"    Length of each dimension of mass cube: "<<endl;
	cout <<"    1st = "<< mass_.SizeX() <<", 2nd = "<< mass_.SizeY() ;
	cout <<", 3rd = "<< mass_.SizeZ() <<endl;

	// We have to initialize the distance to redshift conversion interpolator
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
	double mind = codist[0];
	double maxd = codist[codist.size()-1];
	dist2z_.DefinePoints(codist,zrs,mind,maxd,2*nz);
	
	fg_nodrho = false; // to check if we are simulating from mass distribution or not
	fg_readvals = false; // to check if header values have been read in
	fg_cleancells = false; // to check haven't cleaned cells
	SFApplied = false;
	RandPos_ = false;
}

Mass2Gal::Mass2Gal(int_8 ng,SimpleUniverse& su, RandomGeneratorInterface& rg)
: ng_(ng) , su_(su) , rg_(rg)
// Constructor for when simulating a catalog without clustering information
// used by SimData class
{
	cout <<"    Mass2Gal constructor: simulating galaxy catalog without clustering ..."<<endl;

	// We have to initialize the distance to redshift conversion interpolator
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
	dist2z_.DefinePoints(codist,zrs,mind,maxd,2*nz);
	
	fg_nodrho = true;// to check if we are simulating from mass distribution or not
	fg_readvals = false;// to check if header values have been read in
	fg_cleancells = false; // to check haven't cleaned cells
	SFApplied = false;
};


// Copy Mass2Gal 
Mass2Gal& Mass2Gal::Set(const Mass2Gal& a)
{

	Nx_ = a.Nx_;
	Ny_ = a.Ny_;
	Nz_ = a.Nz_;
	Dx_ = a.Dx_;
	Dy_ = a.Dy_;
	Dz_ = a.Dz_;	
	Dkx_ = a.Dkx_;
	Dky_ = a.Dky_;
	Dkz_ = a.Dkz_; 
	zref_ = a.zref_;		
	ng_ = a.ng_;				
	idmidz_ = a.idmidz_;
	idmidy_ = a.idmidy_;
	idmidx_ = a.idmidx_;
	mass_ = a.mass_;
	ngals_ = a.ngals_;
	//if (NbDimensions() < 1) 		
//	xvals_ = a.xvals_;
//	yvals_ = a.yvals_;
//	zvals_ = a.zvals_;
//	comdist_ = a.comdist_;
//	zdist_ = a.zdist_;
//	theta_ = a.theta_;
//	phi_ = a.phi_; 
//	MB_ = a.MB_;	
//	typeint_ = a.typeint_;		 
//	extincBmV_ = a.extincBmV_;	
}

//******* Methods  ***********************************************************//
void Mass2Gal::ReadHeader(FitsInOutFile& fin)
// reads in SimLSS fits file header so the parameters of the simulated cube can be used
// NOTE: "x" labelled quantiies refer to x DIRECTION in the simulation, and do not refer
// to dimension 1 of the FITS cube etc, because SimLSS outputs cube so:
// AXIS 1 = Z (RADIAL) AXIS
// AXIS 2 = Y AXIS
// AXIS 3 = X AXIS

{
	cout <<endl<<"    Mass2Gal::ReadHeader()"<<endl;
	
	if(fg_nodrho)
		throw ParmError("No clustering information: cannot read simulation FITS header");
	
	// read in values from fits header
	string DXs, DYs, DZs, DKXs, DKYs, DKZs, NXs, NYs, NZs, zrefs, idmids;

	// Fourier space spacing
	DKXs=fin.KeyValue("DKX");// refers to THIRD dimension
	DKYs=fin.KeyValue("DKY");// refers to SECOND dimension
	DKZs=fin.KeyValue("DKZ");// refers to FIRST dimension
	Dkx_=atof(DKXs.c_str());//(double)(drho_.Info()["DKX"];
	Dky_=atof(DKYs.c_str());//(double)(drho_.Info()["DKY"];
	Dkz_=atof(DKZs.c_str());//(double)(drho_.Info()["DKZ"];

	// Read space cube cell spacing
	DXs=fin.KeyValue("DX");// refers to THIRD dimension
	DYs=fin.KeyValue("DY");// refers to SECOND dimension
	DZs=fin.KeyValue("DZ");// refers to FIRST dimension
	Dx_=atof(DXs.c_str());//(double)(drho_.Info()["DX"];
	Dy_=atof(DYs.c_str());//(double)(drho_.Info()["DY"];
	Dz_=atof(DZs.c_str());//(double)(drho_.Info()["DZ"];

	// Number of pixels in each direction
	NXs=fin.KeyValue("NX");// refers to THIRD dimension
	NYs=fin.KeyValue("NY");// refers to SECOND dimension
	NZs=fin.KeyValue("NZ");// refers to FIRST dimension
	Nx_=(sa_size_t)atof(NXs.c_str());//(double)(drho_.Info()["NX"];
	Ny_=(sa_size_t)atof(NYs.c_str());//(double)(drho_.Info()["NY"];
	Nz_=(sa_size_t)atof(NZs.c_str());//(double)(drho_.Info()["NZ"];

	// Redshift of the centre of the cube
	zrefs=fin.KeyValue("ZREF");
	zref_=atof(zrefs.c_str()); 

	idmids=fin.KeyValue("KZREF"); // KZREF is of type DOUBLE, KZREF = Nz/2
	// -> NOT ZERO INDEXED
	// indices of the centre pixel
	idmidz_=(int)ceil(atof(idmids.c_str())); 
	idmidy_=(int)ceil(atof(NYs.c_str())/2);
	idmidx_=(int)ceil(atof(NXs.c_str())/2);
	
	su_.SetEmissionRedShift(zref_);
	DCref_=su_.RadialCoordinateMpc(); // radial distance of central pixel: Z COORD POS OF CENTRE PIXEL
	
	cout << "    check reading fits header ...."<<endl;
	cout << "    dk's: (DKX,DKY,DKZ)="<< Dkx_ <<","<< Dky_ <<","<< Dkz_ <<endl;
	cout << "    dx's: (DX,DY,DZ)="<< Dx_ <<","<< Dy_ <<","<< Dz_ <<endl;
	cout << "    N's: (NX,NY,NZ)="<< Nx_ <<","<< Ny_ <<","<< Nz_ <<endl;
	cout << "    zref:"<< zref_ <<endl;
	cout << "    index of centre pixel (z,y,x)="<< idmidz_ <<","<< idmidy_ <<",";
	cout << idmidx_ <<endl<<endl;
	
	fg_readvals = true;
	cout <<endl<<"    EXIT Mass2Gal::ReadHeader()"<<endl;
}


sa_size_t Mass2Gal::CleanNegativeMassCells()
// the way SimLSS works means there can be values of delta<-1
// this is equivalent to nonlinear structure formation
// however it is unphysical as it is equivalent to negative mass
// therefore this function does two things:
// (i)  changes delta values -> rho/rho^bar by adding 1 to all cells.
// (ii) if the cell has delta <-1 it sets rho/rho^bar=0 (equiv to first setting delta=-1, then adding 1)
{
	cout <<endl<<"    Mass2Gal::CleanNegativeMassCells()"<<endl;

	if(fg_nodrho)
		throw ParmError("No clustering information: cannot clean cells");
	
	double mean,tmp;
	MeanSigma(mass_, mean, tmp);
	cout <<"    Mean over-density BEFORE cell cleaning = "<<mean<<endl;

	sa_size_t nbad = 0;
	for(sa_size_t iz=0; iz<mass_.SizeZ(); iz++) 
		for(sa_size_t iy=0; iy<mass_.SizeY(); iy++) 
			for(sa_size_t ix=0; ix<mass_.SizeX(); ix++) 
				if (mass_(ix, iy, iz)>-1.)  
					mass_(ix, iy, iz) += (float)(1.);
				else 
					{ mass_(ix, iy, iz) = (float)(0.);  nbad++; }

	cout <<"    Mass2Gal::CleanNegativeMassCells() nbad=" << nbad << endl;

	cout <<"    Length of each dimension of mass cube: "<<endl;
	cout <<"    1st = "<< mass_.SizeX() <<", 2nd = "<< mass_.SizeY() <<", 3rd = ";
	cout << mass_.SizeZ() <<endl;

	fg_cleancells=true;
	
	MeanSigma(mass_, mean, tmp);
	mean_overdensity_ = mean - 1;
	cout <<"    Mean over-density AFTER cell cleaning = "<< mean_overdensity_ <<endl;
	
	cout <<"    EXIT Mass2Gal::CleanNegativeMassCells()"<<endl<<endl;
	return nbad;
}

sa_size_t Mass2Gal::CheckNegativeMassCells()
// double checks there are no cells with negative mass
{
	cout <<endl<<"    Mass2Gal::CheckNegativeMassCells()"<<endl;

	if(fg_nodrho)
		throw ParmError("No clustering information: cannot check cells");

	sa_size_t nbad = 0;
	for(sa_size_t iz=0; iz<mass_.SizeZ(); iz++) 
		for(sa_size_t iy=0; iy<mass_.SizeY(); iy++) 
			for(sa_size_t ix=0; ix<mass_.SizeX(); ix++) 
				if (mass_(ix, iy, iz)<0)  
					nbad++;
		
	if (nbad>0)
		cout << "PROBLEM! Mass2Gal::CheckNegativeMassCells() nbad=" << nbad << endl;

	cout <<"    EXIT Mass2Gal::CheckNegativeMassCells()"<<endl<<endl;
	return nbad;
}

void Mass2Gal::ConvertToMeanNGal(float conv, bool Simple)
// converts rho/rho^bar to a (mean) number of galaxies in the cell
// follows the following logic:
// Define quantities:
// Nc=number of gals in cell, Mc=mass in cell, (N/m)=number of gals per unit mass
// Vc=volume of cell, n=number of gals per unit vol(=number density of gals)
// rhoc=density in cell
// Nc = Mc * (N/m)  
// Mc = (rhoc/rho^bar * rho^bar)*Vc
// (N/m) = N/(rho^bar*V):  could think of N as total number of gals in universe, V as volume of universe
//       = n/rho^bar
// therefore the number of gals in a cell is:
// Nc = (rho/rho^bar * rho^bar)*Vc *  n/rho^bar = rho/rho^bar * (Vc*n)
// so all we have to do is multiply the values in the cells by Vc*n
//
// conv=Vc*n: Vc is the volume of the pixels, n=mean number density of gals calculated from 
//            integrating the Schechter luminosity function
// ALSO SWITCHES AROUND DIMENSIONS SO NOW: ngals_(Nx_,Ny_,Nz_) whereas mass_(Nz_,Ny_,Nx_)
// If Simple = true just puts conv gals in each cell
{ 
	cout <<endl<<"    Mass2Gal::ConvertToMeanNGal()"<<endl;

	if(Simple)
		cout <<"    Simulating "<< conv <<" galaxies in each cell"<<endl;
	
	if(fg_nodrho)
		throw ParmError("No clustering information: cannot convert rho/rho^bar to ngal^bar");
	
	if(!fg_readvals)
		throw ParmError("ERROR!: SimLSS cube information not read in from FITS header");
		
	cout <<"    Conversion value="<< conv <<endl;
	int ndim=3;
	sa_size_t mydim[ndim];
	mydim[0] = Nx_; mydim[1] = Ny_; mydim[2] = Nz_;
	ngals_.SetSize(ndim, mydim);
	cout <<"    Length of each dimension of ngals cube: "<<endl;
    	cout<<"    1st = "<<ngals_.SizeX()<<", 2nd = "<<ngals_.SizeY()<<", 3rd = "<<ngals_.SizeZ()<<endl;

	for(sa_size_t ix=0; ix<mass_.SizeZ(); ix++) // along sim's X-direction 
		for(sa_size_t iy=0; iy<mass_.SizeY(); iy++) // along sim's Y-direction
			for(sa_size_t iz=0; iz<mass_.SizeX(); iz++) { // along sim's Z-direction
				
				if(Simple)
					ngals_(ix,iy,iz) = conv;
				else
					ngals_(ix,iy,iz) = (int_4)floor(mass_(iz,iy,ix)*conv);
				}

	// total number of gals simulated				
	ng_ = ngals_.Sum(); // PLACE WHERE TOTAL NGALS IN SIM SET
	cout <<"    Total number of galaxies="<< ng_ <<endl;
	
	cout <<"    Length of each dimension of ngals cube: "<<endl;
    cout <<"    1st = "<< ngals_.SizeX() <<", 2nd = "<< ngals_.SizeY();
    cout <<", 3rd = "<< ngals_.SizeZ() <<endl;
  
	cout <<"    EXIT Mass2Gal::ConvertToMeanNGal()"<<endl<<endl;

};


sa_size_t Mass2Gal::ApplyPoisson()
// we can Poisson fluctuates the mean number of galaxies in each cell:
{
	cout <<endl<<"    Mass2Gal::ApplyPoisson()"<<endl;

	if(fg_nodrho)
		throw ParmError("No clustering information: cannot apply Poisson fluctuations to ngal^bar");

	sa_size_t totngals = 0;
	for(sa_size_t iz=0; iz<ngals_.SizeZ(); iz++) 
		for(sa_size_t iy=0; iy<ngals_.SizeY(); iy++) 
		  for(sa_size_t ix=0; ix<ngals_.SizeX(); ix++) {
			uint_8 ng = rg_.Poisson(ngals_(ix, iy, iz));
			ngals_(ix, iy, iz) = ng; totngals += ng;
			}
	ng_ = ngals_.Sum();
	cout <<"    EXIT Mass2Gal::ApplyPoisson()"<<endl<<endl;
	return totngals;		
};


sa_size_t Mass2Gal::CreateGalCatalog(int idsim, string fitsname, GalFlxTypDist& gfd, 
    bool extinct, bool AMcut, double SkyArea, bool ZisRad)
{
	if(fg_nodrho)
		throw ParmError("No clustering information: cannot apply Poisson fluctuations to ngal^bar");
	
	Timer tm("CreateGalCatalog");
	
	// Create swap space FITS file structure
	FitsInOutFile swf(fitsname, FitsInOutFile::Fits_Create);	
	SwFitsDataTable gals(swf, 2048);
	gals.AddLongColumn("GalID");
	if (ZisRad) {
		gals.AddFloatColumn("x");
		gals.AddFloatColumn("y");
		}
	else {
		gals.AddFloatColumn("theta");
		gals.AddFloatColumn("phi");
		}
	gals.AddFloatColumn("zs");
	gals.AddFloatColumn("type");
	gals.AddFloatColumn("MB");
	// @todo add option to output more/less info
	//gals.AddFloatColumn("d_com"); 
	//gals.AddFloatColumn("Ext");
	//gals.AddFloatColumn("cz"); // cell redshift
	DataTableRow row = gals.EmptyRow();
	
	// For true z only catalog WITHOUT ABSOLUTE MAG CUT
	string fitsname2 = "ZONLY_"+fitsname;
	FitsInOutFile swf2(fitsname2, FitsInOutFile::Fits_Create);	
	SwFitsDataTable gals2(swf2, 2048);
	gals2.AddLongColumn("GalID");
	gals2.AddFloatColumn("zs");
	DataTableRow row2 = gals2.EmptyRow();
	
	// Interpolate maximum observable absolute magnitude verses redshift
	// If not cutting on absolute magnitude the interpolation function
	// has no meaning
	if(!AMcut) {
		// If not cutting just fill with junk MB(z) function
		for (int i=0;i<10;i++)
			{zv_.push_back(i+1); MBmax_.push_back(10);}
		}
	int nint=1000;
	SInterp1D maxAM(zv_, MBmax_, zv_[0], zv_[zv_.size()-1], nint);
	
	// to create object id
	uint_8 gid0=idsim*100000000000LL; // LL seems to be needed by some compilers
	cout <<"    Simulation ID = "<< gid0 <<endl;
	uint_8 gid;
	uint_8 seq=0,seq2=0;

	cout <<  " Mass2Gal::CreateGalCatalog  start looping over cube cells NCells= " << ngals_.Size() << " NGals=" << ng_ << " ... "<<endl;
	uint_8 ngsum=0;  // counter for checking  
	uint_8 nginfile=0,nginzfile=0;  // count of galaxies going into file, total gals simulated
	
	for(sa_size_t iz=0; iz<ngals_.SizeZ(); iz++)// Z direction (~redshift direction)
        for(sa_size_t iy=0; iy<ngals_.SizeY(); iy++)// Y direction (transverse plane)
            for(sa_size_t ix=0; ix<ngals_.SizeX(); ix++) {// X direction (transverse plane) 
		  
		        // pick a cell
		        int_8 ng=ngals_(ix,iy,iz); // num galaxies in this cell
		        // comoving distance to center of cell:
		        double xc, yc, zc, dc;
		        GetCellCoord(ix,iy,iz, xc, yc, zc); // given pixel indices (ix,iy,iz) get comoving coords
		        double rx,ry,rz,rr,rtet,rphi,rreds;// note names theta,phi reverse of the usual spehric. coord convention
		        double mag,gtype,gext;  // galaxy absolute magnitude and type and internal extinction
		        dc = sqrt(xc*xc+yc*yc+zc*zc);// comoving distance to each pixel center

                for(int ing=0; ing<ng; ing++) { // from gal 1 to gal n in cell...
			
			        if(RandPos_) {
				        // We generate random positions with flat distribution inside the cell
				        rx = xc+rg_.Flatpm1()*(Dx_/2);
				        ry = yc+rg_.Flatpm1()*(Dy_/2);
				        rz = zc+rg_.Flatpm1()*(Dz_/2);
				        }
			        else {
				        // Or we use cell center position
				        rx = xc;
				        ry = yc;
				        rz = zc;
				        }
				
			    if (ZisRad)
				    Conv2ParaCoord(rz,rr);
			    else
				    Conv2SphCoord(rx,ry,rz,rr,rphi,rtet); 
			
			    // convert comoving distance into a redshift
			    rreds = dist2Redshift(rr);
			
			    // draw the galaxy properties
			    mag = DrawMagType(gfd, gtype);
			    gext=0.;// extinction is set to 0 for now
			    ngsum++;// should be same as ng_
			
			    double mAM =  maxAM(rreds); // given redshift of galaxy what's max (faintest) absolute mag possible
			    if((mag<mAM) && (rphi<SkyArea)) { // magnitude cut and sky area selection
				
				    // object id, gal id starts at 1 not 0
				    seq++; // adding up TOTAL number of gals put into file
				    gid = gid0+seq;

                    row[0] = gid;
                    row[1] = rx; 
                    row[2] = ry; 
                    row[3] = rreds; 
                    row[4] = gtype; 
                    row[5] = mag; 
                    //row[5] = nginfile;
					//row[5] = gext;
					//row[7] = creds;
					
				    if (ZisRad) {
					    row[1] = rx; 
					    row[2] = ry; 
					    }
				    else {
					    row[1] = rtet; 
					    row[2] = rphi; 
					    }

				    gals.AddRow(row);
				    nginfile++;
				    }

			    if(rphi<SkyArea) 
				    nginzfile++; // counts ALL in sim
				    
			//For the ZONLY file (only if doing magnitude cut)
			if( rphi<SkyArea&&AMcut ) { // JUST sky area selection
				
				// object id, gal id starts at 1 not 0
				seq2++; // adding up TOTAL number of gals put into file
				gid = gid0+seq2;

				row2[0] = gid;
				row2[1] = rreds; 

				gals2.AddRow(row2);
				nginzfile++; // adding number of gals in file
			    }
				
            }  // end of loop over galaxies in the cell 
        } // end of loop over ix (cells) 
	cout <<" Mass2Gal::CreateGalaxyCatalog() finished loop"<<endl;

	// want to delete swf2 file if no magnitude cut, not sure how to do this
	
	gals.Info()["NAllObject"]=nginzfile;// ALL in sim
	gals.Info()["NBrightObject"]=nginfile;// ALL in file 1
	gals.Info()["MeanOverDensity"]=mean_overdensity_;// important

	if( AMcut ){
		gals2.Info()["NAllObject"]=nginzfile;// ALL in sim
		gals2.Info()["NBrightObject"]=nginfile;// ALL in file 1
		gals2.Info()["MeanOverDensity"]=mean_overdensity_;// important
		}
	else	 {
		if( remove(fitsname2.c_str()) != 0 )
    			cout << "Error deleting ZTRUE file" << endl;
		}

	cout <<" Mass2Gal::CreateGalaxyCatalog() done - ng_ = " << ng_;
	cout << "(?=" <<  ngsum << " ) NGal in file = " << nginfile << endl;
	cout <<"                                        seq = " << seq;
	cout <<" (should == NGal in file)"<<endl;
	cout<<"                                        Total NGal simulated = " << nginzfile <<endl;
	return nginfile;
}

void Mass2Gal::CreateNzHisto(string ppfname, GalFlxTypDist& gfd, bool extinct, double SkyArea)
{
	if(fg_nodrho)
        throw ParmError("No clustering information: cannot create N(z) Histo");
	
	Timer tm("CreateNzHisto");
	
	cout <<  " Mass2Gal::CreateNzHisto start looping over cube cells NCells= ";
	cout << ngals_.Size() << " NGals=" << ng_ << " ... "<<endl;
	uint_8 ngsum=0;   // counter for checking  
	uint_8 nginfile=0;// count of galaxies in fits
	
	int nbin=10000;
	double minz=0, maxz=10;
	Histo nz(minz,maxz,nbin);
	
	for(sa_size_t iz=0; iz<ngals_.SizeZ(); iz++) // Z direction (~redshift direction)
        for(sa_size_t iy=0; iy<ngals_.SizeY(); iy++) // Y direction (transverse plane)
		    for(sa_size_t ix=0; ix<ngals_.SizeX(); ix++) { // X direction ( transverse plane) 
			  
			  // pick a cell
			  int_8 ng = ngals_(ix,iy,iz); // num galaxies in this cell
			  // comoving distance to center of cell:
			  double xc, yc, zc, dc;
			  GetCellCoord(ix,iy,iz, xc, yc, zc);
			  //cout << "(xc,yc,zc) = ("<<xc<<","<<yc<<","<<zc<<") ";
			  double rx,ry,rz,rr,rtet,rphi,rreds; // note names theta,phi reverse of the usual sph. coord convention
			  double mag,gtype,gext;  // galaxy absolute magnitude and type and internal extinction
			  dc = sqrt(xc*xc+yc*yc+zc*zc);
			  //double creds = dist2Redshift(dc);  // the cell center redshift
			  //cout <<" dc = "<<dc<<", z_c = "<<creds<<endl;
			  for(int ing=0; ing<ng; ing++) { // from gal 1 to gal n in cell...
		    	  	
				if(RandPos_) {
					// We generate random positions with flat distribution inside the cell
					rx = xc+rg_.Flatpm1()*(Dx_/2);
					ry = yc+rg_.Flatpm1()*(Dy_/2);
					rz = zc+rg_.Flatpm1()*(Dz_/2);
					}
				else {
					rx = xc;
					ry = yc;
					rz = zc;
					}
				
				Conv2SphCoord(rx,ry,rz,rr,rphi,rtet); 
			
				rreds = dist2Redshift(rr);
				mag = DrawMagType(gfd, gtype);
				gext = 0.;
				ngsum++; 
			
				if(rphi<SkyArea)  // sky area selection ONLY
					nz.Add(rreds);
				
				
                } // end of loop over galaxies in the cell 
		    } // end of loop over ix (cells) 
	cout <<" Mass2Gal::CreateNzHisto() finished loop"<<endl;
	
	POutPersist pos(ppfname);
	//FitsInOutFile fos(fitsname, FitsInOutFile::Fits_Create);
	pos << nz;
	//FitsManager::Write(fos, nz);
	
	cout<<" Mass2Gal::CreateNzHisto() done - ng_= " << ng_ << "(?=" <<  ngsum << endl; 
}

sa_size_t Mass2Gal::CreateTrueZFile(int idsim, string fitsname, double SkyArea)
{
	if(fg_nodrho)
		throw ParmError("No clustering information: cannot apply Poisson fluctuations to ngal^bar");
	
	Timer tm("CreateGalCatalog");
	
	// We create first a datatable with the fits file as the swap space 
	FitsInOutFile swf(fitsname, FitsInOutFile::Fits_Create);	
	SwFitsDataTable gals(swf, 2048);
	gals.AddLongColumn("GalID");
	gals.AddFloatColumn("z");
	DataTableRow row = gals.EmptyRow();
	
	cout << "    Mass2Gal::CreateTrueZFile start looping over cube cells NCells= " << ngals_.Size() << " NGals=" << ng_ << " ... "<<endl;
	uint_8 ngsum=0;  // counter for checking loop over whole catalog
	uint_8 nginfile=0;  // count of galaxies going into file
	
	// to create object id
	uint_8 gid0=idsim*100000000000LL; // LL seems to be needed by some compilers
	cout << "    Simulation ID = "<<gid0<<endl;
	uint_8 gid;
	uint_8 seq=0;// counter for checking ngal file1 = ngal file 2
	
    for(sa_size_t iz=0; iz<ngals_.SizeZ(); iz++) // Z direction (~redshift direction)
        for(sa_size_t iy=0; iy<ngals_.SizeY(); iy++) // Y direction (transverse plane)
            for(sa_size_t ix=0; ix<ngals_.SizeX(); ix++) { // X direction ( transverse plane) 
		  
                // pick a cell
		        int_8 ng = ngals_(ix,iy,iz); // num galaxies in this cell
		        // comoving distance to center of cell:
		        double xc, yc, zc, dc;
		        GetCellCoord(ix,iy,iz, xc, yc, zc); // given pixel indices (ix,iy,iz) get comoving coords
		        double rx,ry,rz,rr,rtet,rphi,rreds;// note names theta,phi reverse of the usual spehric. coord convention
		        //double mag,gtype,gext;  // galaxy absolute magnitude and type and internal extinction
		        dc = sqrt(xc*xc+yc*yc+zc*zc);// comoving distance to each pixel center

		        for(int ing=0; ing<ng; ing++) { // from gal 1 to gal n in cell...
			
			        if(RandPos_) {
            				// We generate random positions with flat distribution inside the cell
				        rx = xc+rg_.Flatpm1()*(Dx_/2);
				        ry = yc+rg_.Flatpm1()*(Dy_/2);
				        rz = zc+rg_.Flatpm1()*(Dz_/2);
				        }
			        else {
				        // Or we use cell center position
				        rx = xc;
				        ry = yc;
				        rz = zc;
				        }
				
			        Conv2SphCoord(rx,ry,rz,rr,rphi,rtet); 
			
			        // convert comoving distance into a redshift
			        rreds=dist2Redshift(rr);
			
			        ngsum++;// should be same as ng_
			
			        if(rphi<SkyArea) { // JUST sky area selection
				
				        // object id, gal id starts at 1 not 0
				        seq++; // adding up TOTAL number of gals put into file
				        gid = gid0+seq;

				        row[0] = gid;
				        row[1] = rreds; 

				        gals.AddRow(row);
				        nginfile++;// adding number of gals in file

			            }// end of if gal is in cat
		            }// end of loop over galaxies in the cell 
		        }// end of loop over cells
	cout <<" Mass2Gal::CreateTrueZFile() finished loop"<<endl;
	
	gals.Info()["NAllObject"]=ng_;// ALL in sim
	gals.Info()["NBrightObject"]=nginfile;// ALL in file 1 (ng_ minus angle cut)

	cout <<"    Number of galaxies in file = "<<nginfile<<endl;	
	cout <<"    Mass2Gal::CreateTrueZFile() done - ng_ = " << ng_;
	cout << "(?=" <<  ngsum << " ) NGal in file = " << nginfile << endl; 
	cout <<"                                        seq = " << seq<<" (should == NGal in file)"<<endl;
	return nginfile;
}

sa_size_t Mass2Gal::CreateSimpleCatalog(int idsim, string fitsname, double SkyArea)
{
	if(fg_nodrho)
		throw ParmError("No clustering information: cannot apply Poisson fluctuations to ngal^bar");
	
	Timer tm("CreateSimpleCatalog");
	
	// We create first a datatable with the fits file as the swap space 
	FitsInOutFile swf(fitsname, FitsInOutFile::Fits_Create);	
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
	
	cout << "    Mass2Gal::CreateSimpleCatalog start looping over cube cells NCells= " << ngals_.Size() << " NGals=" << ng_ << " ... "<<endl;
	uint_8 ngsum=0;  // counter for checking loop over whole catalog
	uint_8 nginfile=0;  // count of galaxies going into file
	
	// to create object id
	uint_8 gid0=idsim*100000000000LL; // LL seems to be needed by some compilers
	cout << "    Simulation ID = "<< gid0 <<endl;
	uint_8 gid;
	uint_8 seq=0;// counter for checking ngal file1 = ngal file 2
	
    for(sa_size_t iz=0; iz<ngals_.SizeZ(); iz++) // Z direction (~redshift direction)
        for(sa_size_t iy=0; iy<ngals_.SizeY(); iy++) // Y direction (transverse plane)
            for(sa_size_t ix=0; ix<ngals_.SizeX(); ix++) { // X direction ( transverse plane) 
		  
		        // pick a cell
		        int_8 ng = ngals_(ix,iy,iz); // num galaxies in this cell
		        // comoving distance to center of cell:
		        double xc, yc, zc, dc;
		        GetCellCoord(ix,iy,iz, xc, yc, zc); // given pixel indices (ix,iy,iz) get comoving coords
		        double rx,ry,rz,rr,rtet,rphi,rreds;// note names theta,phi reverse of the usual spehric. coord convention
		        //double mag,gtype,gext;  // galaxy absolute magnitude and type and internal extinction
		        dc = sqrt(xc*xc+yc*yc+zc*zc);// comoving distance to each pixel center

		        for(int ing=0; ing<ng; ing++) { // from gal 1 to gal n in cell...
			
			    if(RandPos_) {
				    // We generate random positions with flat distribution inside the cell
				    rx = xc+rg_.Flatpm1()*(Dx_/2);
				    ry = yc+rg_.Flatpm1()*(Dy_/2);
				    rz = zc+rg_.Flatpm1()*(Dz_/2);
				    }
			    else {
				    // Or we use cell center position
				    rx = xc;
				    ry = yc;
			        	rz = zc;
				    }
				
			    Conv2SphCoord(rx,ry,rz,rr,rphi,rtet); 
			
			    // convert comoving distance into a redshift
			    rreds=dist2Redshift(rr);
			
			    ngsum++;// should be same as ng_
			
			    if(rphi<SkyArea) { // JUST sky area selection
				
				    // object id, gal id starts at 1 not 0
				    seq++; // adding up TOTAL number of gals put into file
				    gid = gid0+seq;

				    /*row[0] = gid;
				    row[1] = rtet;
				    row[2] = rphi;
				    row[3] = rreds; 
				    row[4] = rr; 
				    row[5] = rx; 
				    row[6] = ry; 
				    row[7] = rz;*/

				    row[0] = rtet;
				    row[1] = rphi;
				    row[2] = rreds; 
				    row[3] = rr; 
				    row[4] = rx; 
				    row[5] = ry; 
				    row[6] = rz;
				    gals.AddRow(row);
				    nginfile++;// adding number of gals in file

			        }// end of if gal is in cat
		        }// end of loop over galaxies in the cell 
		    }// end of loop over cells
	cout <<" Mass2Gal::CreateSimpleCatalog() finished loop"<<endl;
	
	gals.Info()["NAllObject"]=ng_;// ALL in sim
	gals.Info()["NBrightObject"]=nginfile;// ALL in file 1 (ng_ minus angle cut)

	cout <<"    Number of galaxies in file = "<< nginfile <<endl;	
	cout <<"    Mass2Gal::CreateSimpleCatalog() done - ng_ = " << ng_;
	cout << "(?=" <<  ngsum << " ) NGal in file = " << nginfile << endl; 
	cout <<"                                        seq = " << seq<<" (should == NGal in file)"<<endl;
	return nginfile;
};


void Mass2Gal::MaxAbsMag()
// Calculate maximum absolute magnitude that 
// can be observed with LSST as a function of 
// redshift
// Needs SimData class to calculate k-correction
{

    double lmin=5e-8, lmax=2.5e-6;

    // Read in CWWK SEDs
    string sedFile = "CWWK.list";
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax);
    vector<SED*> sedArray=readSedList.getSedArray();
    int nSEDs = sedArray.size();
    
    // Read in LSST filters
    string lsstFilterFile = "LSST.filters";
	ReadFilterList readLSSTfilters(lsstFilterFile);
	readLSSTfilters.readFilters(lmin,lmax);
	vector<Filter*> LSSTfilters=readLSSTfilters.getFilterArray();
	int nFilters = LSSTfilters.size();
	
	// Read in GOODS B filter
	string goodsBFilterFile = "GOODSB.filters";
	ReadFilterList readBfilter(lsstFilterFile);
	readBfilter.readFilters(lmin,lmax);
	vector<Filter*> BFilter=readBfilter.getFilterArray();
	
	// Initialise SimData class
	string outcat="tmp";
	PhotometryCalcs photoCalcs(lmin, lmax);
	
	// 5-sig depth for point sources (10 years + 0.1 mag)
	// these should be input variables
	double udepth =27.5505+0.1;
	double gdepth =28.7324+0.1;
	double rdepth =28.45+0.1;
	double idepth =27.7536+0.1;
	double zdepth =27.0585+0.1;
	double ydepth =25.8424+0.1;
	vector<double> depths;
    depths.push_back(udepth);
    depths.push_back(gdepth);
    depths.push_back(rdepth);
    depths.push_back(idepth);
    depths.push_back(zdepth);
    depths.push_back(ydepth);
	
	// Loop over redshifts
	double zmin=0.05, zmax=4;
	int nz=500;
	double dz=(zmax-zmin)/(nz-1);
	
	for (int i=0; i<nz;i++) {
		double z=zmin+i*dz;
		zv_.push_back(z);
		
		// distance modulus
		su_.SetEmissionRedShift(z);
		double dL = su_.LuminosityDistanceMpc();
		double mu = 5.*log10(dL) + 25.;
		
		// MB = Xdepth - mu - kBX
		// To get MAX MB want Xdepth MAX and kBX MIN 
		// So want to find: max(Xdepth - kBX)		
		
		double eps=1e-10;
		double dmsMax=eps;
		for (int j=0; j<nSEDs; j++) {
		    for (int k=0; k<nFilters; k++) {
		        double dms = depths[k] - photoCalcs.Kcorr(z,(*sedArray[j]),
		                                        (*LSSTfilters[k]),(*BFilter[0]));
		        if (dms>dmsMax)
		            dmsMax=dms;
		            
		        }   
		    }
		    
		MBmax_.push_back(dmsMax - mu);
		}

};


void Mass2Gal::ApplyPZConv(string pzcf)
// Apply photometric redshift smearing to simulated 
// galaxy catalog and random catalog
{

	cout <<" Mass2Gal::ApplyPZConv()"<<endl;

	if (!SFApplied)
		throw ParmError("ERROR: have not applied selection function to simulation");
	
	int ndim=3;
	sa_size_t mydim[ndim];
	mydim[0]=ngals_.SizeX(); mydim[1]=ngals_.SizeY(); mydim[2]=ngals_.SizeZ();
	
	if (ngalssm_.NbDimensions()<2) {
		ngalssm_.SetSize(ndim, mydim);
		ngalssm_=0.; // does this definitely zero all elements?
		}
	if (randgsm_.NbDimensions()<2) {
		randgsm_.SetSize(ndim, mydim);
		randgsm_=0;
		}
	
	SInterp1D pzcfi;
	pzcfi.ReadXYFromFile(pzcf);
	double dzmin=pzcfi.XMin(), dzmax=pzcfi.XMax();
	// integrate over distribution
	double IntDist=IntegrateFunc(pzcfi,dzmin,dzmax);
	
	for(sa_size_t ix=0; ix<ngals_.SizeX(); ix++)  {
		cout << "ix="<< ix <<"/"<< ngals_.SizeX() <<endl;
		for(sa_size_t iy=0; iy<ngals_.SizeY(); iy++) 
			for(sa_size_t iz=0; iz<ngals_.SizeZ(); iz++) {

				// current pixel
				int_8 rg = randgsm_(ix,iy,iz); // random value in this pixel
				int_8 ng = ngalssm_(ix,iy,iz); // num galaxies in this pixel
				double tmp1,tmp2,zc;
				GetCellZBounds(ix,iy,iz,tmp1,zc,tmp2);// redshift of this pixel center
				
				// loop over all pixels along z dimension with x,y fixed to current
				// pixel's values
				// for all pixels along this z dim smear galaxies
				double sumg=0, sumgr=0;
				for (int izz=0; izz<ngals_.SizeZ(); izz++) {
					
					double zl,zh,tmp;
					GetCellZBounds(ix,iy,izz,zl,tmp,zh);

					// difference between pixel lowest redshift and current pixel's 
					// central redshift
					double dzl=zc-zl;
					// difference between pixel highest redshift and current pixel's 
					// central redshift
					double dzh=zc-zh;
					
					// want to integrate pzcfi between dzh to dzl
					if (dzh>dzl)
						cout <<"somthing's wrong"<<endl;
					double IntChunk=IntegrateFunc(pzcfi,dzh,dzl);

					// fraction of galaxies of current pixel that are smeared into this pixel
					double frac_sm = IntChunk/IntDist; 

					// add fraction of galaxies in current pixel to this pixel
					ngalssm_(ix,iy,izz) += round(ng*frac_sm);
					randgsm_(ix,iy,izz) += round(rg*frac_sm);
					
					sumg+=ng*frac_sm;
					sumgr+=round(ng*frac_sm);
					}
					
					//if (sumg!=sumgr)
					//	cout <<"tmp"<<endl;
						//do something
						
				
				}
		}

	double ng_left = ngalssm_.Sum();
	double rg_left = randgsm_.Sum();
	cout << "    Number of simulated galaxies left after selection and photo-z errors applied = "<< ng_left<<endl;
	cout << "    Number of random galaxies left after selection and photo-z errors applied = "<< rg_left<<endl;

	cout <<" END Mass2Gal::ApplyPZConv()"<<endl<<endl;
};


void Mass2Gal::GetCellZBounds(sa_size_t i, sa_size_t j, sa_size_t k,
                                        double& zl, double& zc, double& zh)
{

	double x,y,z,dc;
	GetCellCoord(i,j,k,x,y,z);
	
	dc = sqrt(x*x+y*y+z*z);// comoving distance to each pixel center
	zc=dist2Redshift(dc);
	zh=dist2Redshift(dc+Dz_/2);
	zl=dist2Redshift(dc-Dz_/2);
	
};


void Mass2Gal::SetRandomGrid(int mean_dens)
{
	mean_dens_ = mean_dens;
	cout <<"    Setting mean density of unclustered distribution to "<<mean_dens_<<endl<<endl;
};


TArray<r_8> Mass2Gal::ExtractSubArray(sa_size_t x1, sa_size_t x2, sa_size_t y1, 
                                        sa_size_t y2, sa_size_t z1, sa_size_t z2)
// x1, x2 refer to 3RD dimension of array
// y1, y2 refer to 2ND dimension of array
// z1, z2 refer to 1ST dimension of array
{

	cout <<"    Extracting the following pixels ... "<<endl;
	cout <<"    1st dim (z dim): "<< z1 <<","<< z2 <<endl;
	cout <<"    2nd dim (y dim): "<< y1 <<","<< y2 <<endl;
	cout <<"    3rd dim (x dim): "<< x1 <<","<< x2 <<endl;

	cout <<"    Size of mass array = "<< mass_.SizeX() <<","<< mass_.SizeY();
	cout <<","<< mass_.SizeZ() <<endl;
	TArray<r_8> submass = mass_(Range(z1,z2),Range(y1,y2),Range(x1,x2)).PackElements();
	return submass;


};


TArray<r_8> Mass2Gal::ExtractSubArray(double Z, sa_size_t nx, sa_size_t ny, sa_size_t nz)
// Assume the dimensions of the SimLSS cube are this way round: (Nz,Ny,Nx)
// BUT the nx,ny,nz refer to a cube with dimensions the OTHER way round: (nx,ny,nz)
{

	sa_size_t i,j,k;
	double zclose = FindPixelAtZ(Z,nz,ny,nx,i,j,k);
	cout <<"    Pixel closest to z = "<< Z <<" is "<< i <<","<< j <<","<< k;
	cout <<" (has z = "<<zclose<<")"<<endl;


	sa_size_t iz_start, iz_end, iy_start, iy_end, ix_start, ix_end;
	// need to find if nx,ny,nz are even or not
	GetRange(i,nx,iz_start,iz_end);
	GetRange(j,ny,iy_start,iy_end);
	GetRange(k,nz,ix_start,ix_end);
	
	cout <<"    Extracting the following pixels ... "<<endl;
	cout <<"    1st dim: "<< ix_start <<","<< ix_end <<endl;
	cout <<"    2nd dim: "<< iy_start <<","<< iy_end <<endl;
	cout <<"    3rd dim: "<< iz_start <<","<< iz_end <<endl;
	
	TArray<r_8> submass = mass_(Range(ix_start,ix_end), Range(iy_start,iy_end), Range(iz_start,iz_end)).PackElements();
	return submass;
};


void Mass2Gal::WriteLFFits(Schechter lumfunc, double schmin, double schmax, string outfile,int npt)
// Writes out luminosity function to a FITS file
// Note NOT in units of density, because *Vol
{

	if(!fg_readvals)
		throw ParmError("ERROR!: SimLSS cube information not read in from FITS header");
		
	FitsInOutFile swf(outfile, FitsInOutFile::Fits_Create);	
	SwFitsDataTable LF(swf, 2048);
	LF.AddFloatColumn("absmagB");
	LF.AddFloatColumn("phiL");
	DataTableRow row = LF.EmptyRow();

	double Vol=Dx_*Dy_*Dz_*Nx_*Ny_*Nz_;
	cout << "   print Vol"<<endl;
	cout << "   (Dx,Dy,Dz)=("<<Dx_<<","<<Dy_<<","<<Dz_<<")"<<endl;
	cout << "   (Nx,Ny,Nz)=("<<Nx_<<","<<Ny_<<","<<Nz_<<")"<<endl;

	double dm=(schmax-schmin)/(npt-1);
	for(int i=0;i<npt;i++)
		{
		double M=schmin+i*dm;
		row[0]=M;
		row[1]=lumfunc(M)*Vol;
		LF.AddRow(row);
		}
};


TVector<r_8> Mass2Gal::ReturnGridSpec()
// Returns number of pixels and pixel size in a vector
// Can use this as input into other classes
// Useful when computing power spectrum of the grid
// Can calculate Dk via: 2*pi/(N*cellsize)
{

	if(fg_nodrho)
		throw ParmError("No clustering information: cannot return simulation grid specs");
		
	if(!fg_readvals)
		throw ParmError("ERROR!: SimLSS cube information not read in from FITS header");
		
	int row=4; 
	TVector<r_8> gridspec(row); 
	gridspec(0)=Nx_; 
	gridspec(1)=Ny_; 
	gridspec(2)=Nz_; 
	gridspec(3)=Dx_;
	return gridspec;
};


void Mass2Gal::GetCellCoord(sa_size_t i, sa_size_t j, sa_size_t k, double& x, double& y, double& z)
// i refers to 1ST dim, j refers to 2ND dim, k refers to 3RD dim
{
	if(!fg_readvals)
		throw ParmError("ERROR!: SimLSS cube information not read in from FITS header");
		
	x=GetCellX(i);
	y=GetCellY(j);
	z=GetCellZ(k);

	return;
};


void Mass2Gal::Conv2SphCoord(double x, double y, double z, double& r, double& phi, double& theta)
// NOTE : the names theta, phi are reversed compared to the usual (ISO 31-11) convention 
{
	r=sqrt(x*x+y*y+z*z);
	theta=atan2(y,x);
	double zoverr=z/r;
	phi=acos(zoverr);
};


double Mass2Gal::DrawMagType(GalFlxTypDist& gfd, double& type)
{
// draw galaxy absolute magnitude and type from distribution
//
// GalFlxTypDist& gfd is 2D probability distribution of broad type and
// absolute magnitude in B band
// how to assign more precise SED type 
// DEFINE:
// type 0  = El_cww_fix2.txt
// type 10 = Sbc_cww_fix2.txt
// type 20 = Scd_cww_fix2.txt
// type 30 = Im_cww_fix2.txt
// type 40 = SB3_kin_fix2.txt
// type 50 = SB2_kin_fix2.txt

    int a1=0,  b1=5;    //early (1)
    int a2=6,  b2=40;   //late  (2)
    int a3=41, b3=50;   //starburst (3)

// Then type 1 (early) : a1-b1 (0-5)
//		type 2 (late)  : a2-b2 (6-40)
//		type 3 (SBurst): a3-b3 (41-50)
// Choice is same as Dahlen et al arXiv:0710.5532
//
// The intermediate spectra are interpolated via:
// type 0: 1.0*type0+0.0*type10
// type 1: 0.9*type0+0.1*type10
// type 2: 0.8*type0+0.2*type10
// type 3: 0.7*type0+0.3*type10
// type 4: 0.6*type0+0.4*type10
// type 5: 0.5*type0+0.5*type10
// type 6: 0.4*type0+0.6*type10
// type 7: 0.3*type0+0.7*type10
// type 8: 0.2*type0+0.8*type10
// type 9: 0.1*type0+0.9*type10
// my original choice was: a1=0,b1=8,a2=9,b2=36,a3=37,b3=50

// The type numbers are simulated as 1.XX, 2.XX, 3.XX
// where the 1,2,3 represent the original broad band type
// and the XX is the interpolated template type number
// Therefore type values will be:
// 1.00 - 1.05
// 2.06 - 2.40
// 3.41 - 3.50

// Simulation of extinction E(B-V)=extincBmV_
// Random internal galactic exinction value for each galaxy
// which extinction law to use?
// Early types 0-0.1?
// Late types  0-0.3?
// Starburst   0-0.3?
// TBD
//CHECK
//MB_.SetSize(ng_);
//typeint_.SetSize(ng_);
//extincBmV_.SetSize(ng_);
    int typ; double mag;

    int t123;
    gfd.GetGalaxy(typ,mag); 
    t123=typ; // simulates whether galaxy is early (if=1), late (if=2) or starburst (if=3)
	
    // simulate more precise type
    switch (t123) {
        case 1 : {
	        double typ1=floor(a1+(b1-a1)*rg_.Flat01()+0.5); 
	        type=1.0+typ1/100.0; }
	        break;
	    case 2 :
	        { double typ2=floor(a2+(b2-a2)*rg_.Flat01()+0.5); 
		    type=2.0+typ2/100.0; }
	        break;
	    case 3 :
	        { double typ3=floor(a3+(b3-a3)*rg_.Flat01()+0.5); 
		    type=3.0+typ3/100.0; }
	        break;
        default:
	        throw RangeCheckError("Mass2Gal::DrawMagType()/ERROR bad type < 1 or > 3 !");
	        break;
	    }
    return mag;
};


double Mass2Gal::FindPixelAtZ(double Z, sa_size_t nx, sa_size_t ny, sa_size_t nz, sa_size_t& i, sa_size_t& j, sa_size_t& k)
// Assume the dimensions of the SimLSS cube are this way round: (Nz,Ny,Nx) so 
// ix,iy,iz must be using in GetCellCoord() accordingly
//
// nz is number of sub-array pixels along z dimension IE along 1ST dimension of mass_
// ny is number of sub-array pixels along y dimension IE along 2ND dimension of mass_
// nx is number of sub-array pixels along x dimension IE along 3RD dimension of mass_
//
// syntax: x refers to 3RD dim, y refers to 2ND dim, z refers to 1ST dim
// 	   i refers to 1ST dim, j refers to 2ND dim, k refers to 3RD dim
{

	// assume sub-array is close as possible to zero-index pixel edges
	// what is the sub-array pixel center
	sa_size_t zcen_min = nz/2,ycen_min= ny/2,xcen_min= nx/2; 
	// assume sub-array is close as possible to highest-index pixel edges
	// what is the sub-array pixel center
	sa_size_t zcen_max =( (mass_.SizeX()-1)+(mass_.SizeX()-nz) )/2;// because SizeX() == 1ST dim
	sa_size_t ycen_max =( (mass_.SizeY()-1)+(mass_.SizeY()-ny) )/2;// because SizeY() == 2ND dim
	sa_size_t xcen_max =( (mass_.SizeZ()-1)+(mass_.SizeZ()-nx) )/2;// because SizeZ() == 3RD dim

	// sub-array pixel center must be >= these values, and <= these
	cout << "    Center pixel 1st (z) dimension must be >= "<< zcen_min <<" & <= "<< zcen_max <<endl;
	cout << "    Center pixel 2nd (y) dimension must be >= "<< ycen_min <<" & <= "<< ycen_max <<endl;
	cout << "    Center pixel 3rd (x) dimension must be >= "<< xcen_min <<" & <= "<< xcen_max <<endl;

	double diffz=1e10, zclose = -1; 
	for(sa_size_t ix=0; ix<mass_.SizeZ(); ix++) // X direction, 3rd dim (transverse plane) 
        	for(sa_size_t iy=0; iy<mass_.SizeY(); iy++) // Y direction, 2nd dim (transverse plane)
	    		for(sa_size_t iz=0; iz<mass_.SizeX(); iz++) { // Z direction, 3rd dim (redshift direction) 

		  		// comoving distance to center of cell:
		  		double xc, yc, zc, dc, zsc;
		  		GetCellCoord(iz,iy,ix,xc,yc,zc); // iz,ix reversed accordingly
		  		dc = sqrt(xc*xc+yc*yc+zc*zc);// comoving distance to each pixel center
				// convert comoving distance into a redshift
				zsc=dist2Redshift(dc);

				if ( (abs(zsc-Z)<diffz)&&(ix!=0)&&(iy!=0)&&(iz!=0) ) {
					// perform extra check to make sure sub-array can 
					// actually be within the array at each pixel
					
					if ( (iz>=zcen_min)&&(iz<=zcen_max)&&(iy>=ycen_min)&&(iy<=ycen_max)
					                        &&(ix>=xcen_min)&&(ix<=xcen_max) ) {
						i = iz; 
						j = iy;
						k = ix; 
						diffz = abs(zsc-Z);
						zclose = zsc;
						}
					}
				}

	return zclose;

};


void Mass2Gal::GetRange(sa_size_t i, sa_size_t ni, sa_size_t& istart, sa_size_t& iend)
{

	if (ni%2==0) { // if size is EVEN (a%b gives remainder after dividing a/b)
		cout <<"    sub-array dimension pixel number is even: ni = "<<ni<<endl;
		istart = i - (ni/2 -1);
		iend =   i  + ni/2;
		}
	else {// ni is ODD
		cout <<"    sub-array dimension pixel number is odd: ni = "<<ni<<endl;
		sa_size_t ni2 = floor((double)ni/2);
		istart = i - ni2;
		iend =   i + ni2; 
		}

};

/*void Mass2Gal::ApplySF(SelectionFunctionInterface& sf)
// Apply and correct for selection function to simulated
// galaxy catalog AND random catalog
{

	cout <<"    Mass2Gal::ApplySF()"<<endl;
	selfuncp_=&sf;

	int ndim=3;
	sa_size_t mydim[ndim];
	mydim[0]=ngals_.SizeX(); mydim[1]=ngals_.SizeY(); mydim[2]=ngals_.SizeZ();
	
	// i.e. if ngalssm not already initialised, initialise
	if (ngalssm_.NbDimensions()<2)
		{
		ngalssm_.SetSize(ndim, mydim);
		ngalssm_=0.; // zero all elements
		}
	if (randgsm_.NbDimensions()<2)
		{
		randgsm_.SetSize(ndim, mydim);
		randgsm_=0;
		}

	double ng_start = ngals_.Sum();
	cout << "    Number of simulated galaxies BEFORE selection applied = "<< ng_start<<endl;
	
	for(sa_size_t ix=0; ix<ngals_.SizeX(); ix++)
		{
		//cout <<" outer loop "<<ix<<" of "<<ngals_.SizeX()<<endl;
		for(sa_size_t iy=0; iy<ngals_.SizeY(); iy++)
			for(sa_size_t iz=0; iz<ngals_.SizeZ(); iz++)
				{
				
				int_8 ng=ngals_(ix,iy,iz); // num galaxies in this cell
				// comoving distance to center of cell:
				double xc, yc, zc, dc,redshift;
				GetCellCoord(ix,iy,iz,xc,yc,zc); // given pixel indices (ix,iy,iz) get comoving coords
				dc = sqrt(xc*xc+yc*yc+zc*zc);// comoving distance to each pixel center
				redshift=dist2Redshift(dc);
				double phi=(*selfuncp_)(redshift);
				
				double mu;
				uint_8 npoiss;				
				
				// galaxy grid
				mu = ng*phi;
				npoiss=rg_.PoissonAhrens(mu); // Poisson fluctuate
				ngalssm_(ix,iy,iz) = (double)npoiss/phi;
				
				// random grid
				mu=(phi*mean_dens_);
				npoiss=rg_.PoissonAhrens(mu); // Poisson fluctuate
				randgsm_(ix,iy,iz)=(double)npoiss/phi;
				
				}
		}
	
	SFApplied = true;

	double ng_left = ngalssm_.Sum();
	double rg_left = randgsm_.Sum();
	cout << "    Number of simulated galaxies left after selection applied = "<< ng_left<<endl;
	cout << "    Number of random galaxies left after selection applied = "<< rg_left<<endl;
	
cout <<" END Mass2Gal::ApplySF()"<<endl<<endl;
}*/





//CHECK
//void Mass2Gal::ConvertGalDist2Fluct()
//// Converts the number of galaxies in each cell to a fluctuation (delta_n) in each cell
//// delta_n = n/nbar - 1
//// n = number of galaxies in the cell
//// nbar = mean number of galaxies in a cell = total number of galaxies / total number of cells
//{ 
//cout <<endl<<"    Mass2Gal::ConvertGalDist2Fluct()"<<endl;
//int totncells=Nx_*Ny_*Nz_;
//double neff=ngals_.Sum();
//double mult=(double)totncells/neff; // 1 / average number of galaxies per cell
//cout <<"    average number of galaxies in a grid cell="<<1/mult<<endl;
//
//galdens_=ngals_;
//galdens_*=mult;
//galdens_-=1;
//
//
//cout <<"    EXIT Mass2Gal::ConvertGalDist2Fluct()"<<endl<<endl;
//return;
//}




// *********************************************** //
// EVERYTHING BELOW HERE IS PROBABLY REDUNDANT NOW //
// *********************************************** //
//
//void Mass2Gal::SurveyWindow(double Phi, double dbar)
//// Fills array selec_ with 0's where survey does not cover
//// and 1's where it does
//// This is the REAL space window function
//{
//	
//	int ndim=3;
//	sa_size_t mydim[ndim];
//	mydim[0]=mass_.SizeX(); mydim[1]=mass_.SizeY(); mydim[2]=mass_.SizeZ();
//	selec_.SetSize(ndim, mydim);
//
//	for(sa_size_t ix=0; ix<mass_.SizeZ(); ix++) // along sim's X-direction 
//		for(sa_size_t iy=0; iy<mass_.SizeY(); iy++) // along sim's Y-direction
//			for(sa_size_t iz=0; iz<mass_.SizeX(); iz++) // along sim's Z-direction
//				{
//				// find th,ph,dc value of cell
//				double x,y,z,ph,th,dc;
//				GetCellCoord(ix,iy,iz,x,y,z);
//				Conv2SphCoord(x,y,z,dc,ph,th);
//				
//				selec_(ix,iy,iz)=dbar;
//				
//				if (ph>Phi)
//					selec_(ix,iy,iz)=0;
//					
//				}
//
//}
//
//
//double Mass2Gal::AddSurveyWindow(double Phi)
//// Sets all cells with phi>Phi in over-density array mass_
//// to -1:
//// delta == rho/rhobar - 1
//// Therefore if there is ZERO mass in the cell, the cell will 
//// have value -1
//// MUST NOT HAVE CLEANED CELLS ALREADY, MUST NOT CLEAN CELLS AFTER
//// This function is just used to study real space window,
//// therefore the mass_ array must not be used afterwards
//{
//
//	if(fg_cleancells)
//		throw ParmError("ERROR!: Mass2Gal::AddSurveyWindow() Already cleaned cells");
//	
//	TArray<r_8> masstmp = mass_;
//	
//	int_8 nzero=0;
//	for(sa_size_t ix=0; ix<mass_.SizeZ(); ix++) // along sim's X-direction 
//		for(sa_size_t iy=0; iy<mass_.SizeY(); iy++) // along sim's Y-direction
//			for(sa_size_t iz=0; iz<mass_.SizeX(); iz++) // along sim's Z-direction
//				{
//				// find th,ph,dc value of cell
//				double x,y,z,ph,th,dc;
//				GetCellCoord(ix,iy,iz,x,y,z);
//				Conv2SphCoord(x,y,z,dc,ph,th);
//				
//				mass_(ix,iy,iz)=masstmp(ix,iy,iz);
//				
//				if (ph>Phi)
//					{mass_(ix,iy,iz)=-1; nzero++;}
//					
//				}
//
//	fg_apwindow=true;
//	
//	double nfrac = (double)nzero/(mass_.SizeX()*mass_.SizeY()*mass_.SizeZ());
//	return nfrac;
//}


