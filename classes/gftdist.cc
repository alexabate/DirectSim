#include "gftdist.h"

/* -----------------------------------------------------------------------------
   A. Abate , R. Ansari
------------------------------------------------------------------------------*/
// Only really CumulDistZ, DrawM, SimBaseCatalog have methods defined here

//******* CumulDistZ methods *************************************************//

// Set up for doing whole calculation
void CumulDistZ::SetUp(LFParameters& lfpars, SimpleUniverse& su, double zmin,
			double zmax, int nptinteg, int nptinterp, int prt)
{

	lfpars_ = lfpars;
	su_ = su;
	zmin_ = zmin;
	zmax_ = zmax;
	nptinteg_ = nptinteg;
	nptinterp_ = nptinterp;
	cout <<"     Setting up cumulative redshift distribution between "<< zmin_;
	cout <<" < z < "<< zmax_ <<endl;

	// if nptinteg = 1000 it will take around 28s for each nptinterp
	// (so total time in s is 28*nptinterp)
	
	// This code is calculating the equation for F_z(z) in Gorecki et al 2012 in prep
	
	// The below class returns phi(z) = [int phi(M|z) dM]*dV(z)
	SchechterZVol schZ(lfpars_, su_);

	// set up interpolation table for int phi(z) dz
	
	Timer tm("CumulDistZ::SetUp",false);
	tm.Split();
	
	double dz=(zmax-zmin)/(nptinterp_-1);
	for (int i=0; i<nptinterp_; i++) {
		if (prt>0)
			cout <<"     On "<<i+1<<" of "<<nptinterp_<<endl;                           
		double z = zmin + i*dz;
		schZ.SetInteg(0, z, nptinteg_);
		double val=schZ.Integrate();
		zv_.push_back(z);
		scv_.push_back(val);
		
		if (i<10 && i>8) {
		    tm.Split();
		    tm.PartialElapsedTime();
		    cout <<"     "<<i+1<<" redshifts took "<< tm.PartialElapsedTime() <<"s"<<endl;
		    }
		}

	// divide by value at zmax
	for (int i=0; i<nptinterp_; i++)
		scv_[i] /= scv_[nptinterp_-1];

	// check vector
	if (prt>0) {
	cout << "     Checking cumval z vector ... "<<endl;
		for (int i=0; i<nptinterp; i++)
			cout << zv_[i] <<"  "<< scv_[i] <<endl; }

	// create table
	schZint_.DefinePoints(zv_, scv_, zv_[0], zv_[zv_.size()-1], zv_.size()*4);

	/*double ztmp;
	ztmp=5.9;
	cout << "z="<<ztmp<<", cumz="<<schZint_(ztmp)<<endl;
	ztmp=6.;
	cout << "z="<<ztmp<<", cumz="<<schZint_(ztmp)<<endl;
	ztmp=6.1;
	cout << "z="<<ztmp<<", cumz="<<schZint_(ztmp)<<endl;*/

};

// Set up for reading in FITS bintable file
void CumulDistZ::SetUp(string infile, int prt)
{

	cout <<"     Reading in file "<< infile <<endl;
	FitsInOutFile fin(infile, FitsInOutFile::Fits_RO);
	fin.MoveAbsToHDU(2);
	SwFitsDataTable dt(fin, 512, false);
	DataTableRow rowin = dt.EmptyRow();
	int nv=dt.NEntry();
	for (int i=0;i<nv; i++) {
	
		dt.GetRow(i, rowin);
		double z=rowin[0];
		double c=rowin[1];
		zv_.push_back(z);
		scv_.push_back(c);
		}

	// check vector
	if (prt>0) {
	cout << "     Checking cumval z vector ... "<<endl;
		for (int i=0; i<nv; i++)
			cout << zv_[i] <<"  "<< scv_[i] <<endl; }

    zmin_ = zv_[0];
    zmax_ = zv_[nv-1];

	// create interpolation: this won't set to zero outside of zmin, zmax
	schZint_.DefinePoints(zv_, scv_, zmin_, zmax_, nv*4);

};

// Need to finish below
/*// Set up for reading in TEXT file
void CumulDistZ::SetUp(string infile,int prt)
{

	cout <<"    Reading in file "<<infile<<endl;
	ifstream ifs;
	ifs.open(infile.c_str(), ifstream::in);
	if (ifs.fail())
		throw ParmError("ERROR: failed to find cumulative z-dist file");
	sa_size_t nr, nc;
	TArray<double> tab;
	tab.ReadASCII(ifs,nr,nc);

	for (int i=0;i<nc; i++)
		{
		double z=tab(i,0);
		double c=tab(i,1);
		zv_.push_back(z);
		scv_.push_back(c);
		}

	// check vector
	if (prt>0) {
	cout << "Checking cumval z vector ... "<<endl;
		for (int i=0; i<zv_.size(); i++)
			cout << zv_[i] <<"  "<<scv_[i]<<endl; }

	// create table
	schZint_.DefinePoints(zv_,scv_,zv_[0],zv_[zv_.size()-1],zv_.size()*4);

	double ztmp;
	ztmp=5.9;
	cout << "z="<<ztmp<<", cumz="<<schZint_(ztmp)<<endl;
	ztmp=6.;
	cout << "z="<<ztmp<<", cumz="<<schZint_(ztmp)<<endl;
	ztmp=6.1;
	cout << "z="<<ztmp<<", cumz="<<schZint_(ztmp)<<endl;

}*/

// Output the calculated cumulative z function to a FITS binary table file
void CumulDistZ::Output2File(string outfileroot)
{

    // Create filename including the limits of the z in the distribution
    stringstream ss1, ss2;
    ss1 << zmin_; ss2 << zmax_;
	string outfile = outfileroot + "_" + ss1.str() + "z" + ss2.str() + "_cumz.fits";

	// Create swap space FITS file structure
	FitsInOutFile swf(outfile, FitsInOutFile::Fits_Create);	
	SwFitsDataTable func(swf, 2048);
	func.AddFloatColumn("z");
	func.AddFloatColumn("Fz");
	DataTableRow row = func.EmptyRow();

	cout << "     Writing to file ..." << outfile.c_str() << endl;
	int nz = zv_.size();
	for (int i=0; i<nz; i++) {
	
		double z = zv_[i];
		double szv = scv_[i];
		row[0] = z;
		row[1] = szv; 
		func.AddRow(row);
		}
	cout <<endl;

};

// Need to finish below
/*// Output the calculated cumulative z function to a TEXT file
void CumulDistZ::Output2File(string outfile)
{

	ifstream inp;
	ofstream outp;

	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
	  	{
		inp.clear(ios::failbit);
		cout << "     Writing to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		for (int i=0; i<zv_.size(); i++)
			{
			double z=zv_[i];
			double szv=scv_[i];
			outp << z <<"  "<<szv<<endl;
			}
		outp.close();
		}
	else
		cout <<"ERROR file "<<outfile<<" exists"<<endl;
	cout <<endl;

};*/


//******* DrawM methods ******************************************************//

// Set up for calculating table
void DrawM::SetUp(RandomGeneratorInterface& rg, double mmin, double mmax,
    double zmin, double zmax, int nptz, int nptm)
{

		rg_ = rg;
		mmin_ = mmin;
		mmax_ = mmax;
		zmin_ = zmin;
		zmax_ = zmax;
		
		// Get Mmin and Mmax from the cumulative magnitude distribution
		double Mmin, Mmax;
	    cumm_.returnMminMmax(Mmin, Mmax);
		
		if ( (mmin_<Mmin) || (mmax_>Mmax) )
		    throw ParmError("ERROR! magnitude ranges outside of integral range");

        // set array sizes
		mv_.SetSize(nptm);
		zv_.SetSize(nptz);
		int ndim=2;
		sa_size_t mydim[ndim];
		mydim[0]=nptm; mydim[1]=nptz;
		cumval_.SetSize(ndim,mydim);
		dm_=(mmax_-mmin_)/(nptm-1);
		dz_=(zmax_-zmin_)/(nptz-1);
		cout << "     DrawM::SetUp Size of array to be created: "<< nptm <<" x "<< nptz <<endl;
		cout << endl;

		SetArray();

}

// Set up for reading from FITS file 
void DrawM::SetUp(string infile)
{

    // read in table from FITS file
	FitsInFile fin(infile);
	fin >> cumval_; // dim1=m, dim2=z
	
	mmin_ = atof(fin.KeyValue("MMIN").c_str());
	dm_ = atof(fin.KeyValue("DMAG").c_str());
	int nm = atoi(fin.KeyValue("NMAG").c_str());
	mmax_ = dm_*(nm-1) + mmin_;
	
	zmin_ = atof(fin.KeyValue("ZMIN").c_str());
	dz_ = atof(fin.KeyValue("DZ").c_str());
	int nz = atoi(fin.KeyValue("NZ").c_str());
	zmax_ = dz_*(nz-1) + zmin_;
	
	cout << "     Magnitude range from "<<mmin_<<" in "<<nm<<" steps of ";
	cout << dm_ << endl;
	cout << "     Redshift range from "<<zmin_<<" in "<<nz<<" steps of ";
	cout << dz_ << endl;
	
	mv_.ZeroSize();
	zv_.ZeroSize();
	
	mv_.SetSize(nm);
	zv_.SetSize(nz);

	for (int i=0; i<nm; i++) {
		double m = mmin_ + i*dm_;
		mv_(i) = m;
		}
		
	for (int i=0; i<nz; i++) {
		double z = zmin_ + i*dz_;
		zv_(i) = z;
		}
};

// need to finish below
/*
// Set up for reading from TEXT file
void DrawM::SetUp(string infile)
{

	cout <<"    Reading in file "<<infile<<endl;
	ifstream ifs;
	ifs.open(infile.c_str(), ifstream::in);
	if (ifs.fail())
		throw ParmError("ERROR: failed to find cumulative z-dist file");
	sa_size_t nr, nc;
	TArray<double> tab;
	tab.ReadASCII(ifs,nr,nc);

	for (int i=0;i<nc; i++)
		{
		double z=tab(i,0);
		double m=tab(i,1);
		double c=tab(i,2);
		zv_.push_back(z);
		scv_.push_back(c);
		cumval_(i,j)=
		}
};*/

// To draw magnitude given redshift z
double DrawM::Draw(double z)
{
	
	double mindiff=1e10;
	int ibin=-1;

	// find closest redshift to the input z
	for (int i=0; i<zv_.Size(); i++) {
		double diff=std::abs(z-zv_(i));
		if (diff<mindiff) {
			ibin=i;
			mindiff=diff;
			}
		}
	if (ibin<0)
		throw ParmError("DrawM::Draw() No z bin found!");

	int iz=ibin;
	//cout << "DrawM::Draw z="<<z<<", izbin="<<iz<<endl;

	mindiff=1e10;
	ibin=-1;

	// find closest cumval_ in redshift bin iz to the random number val
	double rn=rg_.Flat01();

	for (int i=0; i<mv_.Size(); i++) {
		double diff=std::abs(rn-cumval_(i,iz));

		if (diff<mindiff) {
			ibin=i;
			mindiff=diff;
			}
		}
	if (ibin<0)
		throw ParmError("DrawM::Draw() No m bin found!");

	return mv_(ibin);
};

// Output the calculated cumulative m|z function to a FITS file
void DrawM::Output2File(string outfileroot)
{

    // Create filename including the limits of the z, m in the distribution
    stringstream ss1, ss2, ss3, ss4;
    ss1 << zmin_; ss2 << zmax_;
    ss3 << mmin_; ss4 << mmax_;
	string outfile = outfileroot + "_" + ss1.str() + "z" + ss2.str() + "_";
	outfile += ss3.str() + "m" + ss4.str() + "_cumm.fits";

	cout << "     Creating FITS file "<<outfile<<endl;
	int nm=mv_.Size();
	cout << "     Magnitude range from "<< mmin_ <<" in "<<nm<<" steps of ";
	cout << dm_ << endl;
	
	int nz = zv_.Size();
	cout << "     Redshift range from "<< zmin_ <<" in "<<nz<<" steps of ";
	cout << dz_ << endl;

	FitsInOutFile fos(outfile,FitsInOutFile::Fits_Create);
	fos << cumval_;  // dim1=m, dim2=z
	cout << "     Added array as FITS image"<<endl;
	cout << "     Write header data ... "<<endl;
	fos.WriteKey("MMIN",mmin_,"minimum absolute magnitude");
	fos.WriteKey("DMAG",dm_,"step in magnitudes");
	fos.WriteKey("NMAG",nm,"number of magnitudes");
	fos.WriteKey("ZMIN",zmin_,"minimum redshift");
	fos.WriteKey("DZ",dz_,"step in redshifts");
	fos.WriteKey("NZ",nz,"number of redshifts");
	cout <<"     .... done"<<endl;
	

};

// Need to finish below
/*// Output the calculated cumulative m|z function to a TEXT file
void DrawM::Output2File(string outfile)
{

	ifstream inp;
	ofstream outp;

	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
	  	{
		inp.clear(ios::failbit);
		cout << "     Writing to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);

		for (int i=0; i<mv_.Size(); i++)
			for (int j=0; j<zv_.Size(); j++)
				{ 

				double m=mv_(i);
				double z=zv_(i);
				double cv=cumval_(i,j);
	
				outp << z <<"  "<<m<<"  "<<cv<<endl;
				}
		outp.close();
		}
	else
		cout <<"ERROR file "<<outfile<<" exists"<<endl;
	cout <<endl;

};*/


//******* SimBaseCatalog methods *********************************************//

void SimBaseCatalog::DoSim(long ngal, string outfileroot)
{
	// for timing
	Timer tm("SimBaseCatalog",false);
	
	// Create filename including the limits of the z
    stringstream ss1, ss2;
    ss1 << zmin_; ss2 << zmax_;
	string outfile = outfileroot + "_" + ss1.str() + "z" + ss2.str() + "_gals.fits";

	// Create swap space FITS file structure
	FitsInOutFile swf(outfile,FitsInOutFile::Fits_Create);	
	SwFitsDataTable cat(swf, 2048);
	cat.AddFloatColumn("z");
	cat.AddFloatColumn("am");
	cat.AddFloatColumn("type");
	DataTableRow row = cat.EmptyRow();

	cout << "     Writing to file ..." << outfile.c_str() << endl;
	tm.Split(); int ii=0; int dii=floor(ngal/100);
	for (int i=0; i<ngal; i++) {

		if (i==ii) {
			cout << "     Galaxy "<<i+1<<" of ";
			cout <<ngal<<endl;
			ii+=dii;
			}
		double z=drz_.Draw();
		double am=drm_.Draw(z);
		double typ=(double)drt_(am,z);

		row[0] = z;
		row[1] = am; 
		row[2] = typ;
		cat.AddRow(row);
		}
	tm.Split();
	cout <<endl;
	cout << "     Time to simulate "<<ngal<<" gals was ";
	cout << tm.PartialElapsedTime()<<"s"<<endl;

};


//******************** Below here all is redundant **************************///

/****************************** GALFLXTYPDIST *********************************/
//-- Constructor
GalFlxTypDist::GalFlxTypDist(RandomGeneratorInterface& rg, int prtlev)
 : rg_(rg) , dt1(),  dt2() , dt3() , t1_(dt1) , t2_(dt2) , t3_(dt3)//, sa_(sa)
{
  prtlev_ = prtlev;
  tot_frac_ = 0.;
  totngal_ = 0;
  SchTypes_= false;
}

GalFlxTypDist::GalFlxTypDist(int prtlev)
 : rg_(rgdefault_)  , dt1(),  dt2() , dt3() , t1_(dt1) , t2_(dt2) , t3_(dt3)//, sa_(sa)
{
  prtlev_ = prtlev;
  tot_frac_ = 0.;
  totngal_ = 0;
  SchTypes_= false;
}

//-- Destructor 
GalFlxTypDist::~GalFlxTypDist()
{
}

//-- return true if OK, generate an exception (error) if problem
bool GalFlxTypDist::Check() 
{
  if (NTypes() < 1) {
    cout << "GalFlxTypDist::Check()/Error No galaxy type defined !" << endl;
    throw PException("GalFlxTypDist::Check()/Error No galaxy type defined");
  }
  if (fabs(tot_frac_-1.)>1.e-9) {
    cout << "GalFlxTypDist::Check()/Error - Wrong Golbal Normalization Sum_Types[Fraction] = " 
         << tot_frac_ << " != 1 " << endl;
    throw PException("GalFlxTypDist::Check()/Error - Wrong Golbal Normalization");
  }
  return true;
}

// ADD NEW GALAXY TYPE
int GalFlxTypDist::AddGalType(ClassFunc1D& magf, double magmin, double magmax, double fr, 
                              int nbinmag, int nstep)
// * magf : Magnitude distribution function (i.e. luminosity function)
// * magmin, magmax : Magnitude range
// * fr : Fraction of galaxies having the corresponding type
// * nbinmag : Number of magnitude bins for random generation
// * nstep : number of steps for computing nbinmag equi-population magnitude bins
{
	// first check fraction and mag limit values are sensible
  if ((fr<0.)||(fr>1.)||(magmin>=magmax)) { 
    cout << " GalFlxTypDist::AddGalType()/Error : Bad parameters( fr<0.)||(fr>1.)||(magmin>=magmax)" << endl;
    throw ParmError("GalFlxTypDist::AddGalType()/Error : Bad parameters");
  }
  
  // put nbinmag and nstep to sensible values if necessary
  if (nbinmag < 1)  nbinmag = 1;
  if (nstep<10*nbinmag) nstep = 10*nbinmag;
  
  // some useful vectors
  v_frac_.push_back(fr);   // v_frac_ is a vector which contains the values of the galaxy type fractions
  tot_frac_ += fr;// tot_frac_ was set to 0 in constructor, so this is current total of fractions of all added types so far
  v_sumfrac_.push_back(tot_frac_);// don't know what point of this is?
  
  // print out to help with debugging/checking
  if (prtlev_ > 1) 
    cout << " --- GalFlxTypDist::AddGalType()/Info , GalType=" << v_frac_.size() 
          << " ( NBinMag=" << nbinmag << " NStep=" << nstep << ")" << endl;
		  
  // COMPUTE APPROXIMATELY THE INTEGRAL OF magf 
  // works out bin widths so distribution has equal number of galaxies in 
  // 1) make normalised luminosity function == probability distribution
  Vector vfm(nstep);  // magnitude distribution (luminosity function) vector
  double sum = 0.;
  double delmag = (magmax-magmin)/(double)nstep; // mag spacing given magmin,magmax,nstep: -1?
  double xm = magmin;
  for(sa_size_t kk=0; kk<vfm.Size(); kk++) {
    double fdm = magf(xm);// value of luminosity function defined by magf parameters at magnitude xm
    sum += fdm; // adds up all lf values for each mag between magmin and magmax 
    vfm(kk) = fdm; // luminosity function vector
    xm += delmag; // step up mag value for next loop
  }
  vfm /= (r_8)(sum);  // renormalize the magnitude distribution vector (ie the LF) so it sums to 1
  // basically vfm is input luminosity function given by magf but defined according to magmin, magmax, delmag and normalised so it sums to 1
  
  // print out to help with debugging/checking
  if (prtlev_ > 1) 
    cout << " ... Integral[MagDist," << magmin << "," << magmax << "]=" << sum 
         << " TypeFraction=" << fr << " CurrentTotalFrac=" << tot_frac_ << endl;
		 
  xm = magmin;
  double seuil = 1./nbinmag; // because vfm sums to 1 then if have equi pop bins then each one will sum to 1/nbinmag
  sum = 0.;
  Vector magbins(nbinmag+1, BaseArray::RowVector); //define a vector with nbinmag+1 rows, and specify that it's a row vector
  sa_size_t jm=0;
  magbins(jm) = magmin;   jm++;
  
  // 2) loop over normalised luminosity function vector
  for(sa_size_t kk=0; kk<vfm.Size(); kk++) 
  {
    if (jm > nbinmag) 
	{// doesn't print this: just a check
      cout << " GalFlxTypDist::AddGalType()/Warning jm=" << jm << "/nbinmag=" 
           << nbinmag << " seuil=" << seuil << endl;
      break;
	}
	
    xm += delmag;  sum += vfm(kk); // sum will be getting larger each loop
	
	
    if (sum >= seuil) //decider if reached equipop bin yet
	{
	  // print out to help with debugging/checking
      if (prtlev_ > 2)  
        cout << " ... MagBin[" <<  jm << "] -> " << xm << " ( IntegMagDist=" 
             << seuil << " )" << endl;
	  // record mag value one LF has reached equi pop
      magbins(jm)= xm;  jm++; 
      seuil += (1./nbinmag);
    }
	
  }
  //magbins is nbinmag+1 long after this
  
  // this is to fudge the end if needed I think
  while (jm < (nbinmag+1)) 
  {
  
    if (prtlev_ > 2)  
        cout << " ... Filling MagBin[" <<  jm << "] -> " << magmax << endl;
		
    magbins(jm)= magmax;  jm++; 
  }
  
// Keep the computed magbins vector
  v_mag_.push_back(magbins); // push back a vector into a Tvector???

  return v_mag_.size();
}

// TO PICK GALAXY TYPE AND MAGNTIUDE FROM DISTRIBUTIONS IN ADDGALTYPE
int GalFlxTypDist::GetGalaxy(int& typ, r_8& mag) const 
{

  //SIMULATE TYPE
  double rndt = rg_.Flat01(); // pick a uniform random number between 0 and 1
  size_t styp = 0;
  bool badtype = true;
  
  // think this part randmonly selects styp to be 0, 1 or 2 according to the type fraction
  for(size_t kt=0; kt<NTypes(); kt++) 
	{
    if (rndt <= v_sumfrac_[kt])  {  styp = kt; badtype = false; break; }
	}
	//cout <<"print check: styp="<<styp<<", badtype="<<badtype<<endl;
	//type of gal has been simulated, it is typ=styp+1
	
  //debugging
  if (badtype) 
	{
    cout << " GalFlxTypDist::GetGalaxy()/BUG!!! Bad Type , rnd=" << rndt << endl;
    typ = 1;  mag = 0.; 
    totngal_++;
    return(totngal_);
	}
  
  //SIMULATE MAG GIVEN TYPE
  // v_mag_[0] is (luminosity function) list of magnitudes of type 1; v_mag_[1] is lf of type 2; v_mag_[2] is lf of type 3; etc
  sa_size_t maxidx = v_mag_[styp].Size()-1; //equal to nbinmag i think (max possible index of bin)
  sa_size_t idxmag = (sa_size_t)(rg_.Flat01()*(double)maxidx)+1; // to randomly select a magnitude bin 
  //cout <<"maxidx="<<maxidx<<endl;
  
  // debugging
  if ((idxmag<1)||(idxmag>maxidx)) 
	{
    cout << " GalFlxTypDist::GetGalaxy()/BUG!!! idxmag=" << idxmag << "  >=Max=" << v_mag_[styp].Size() << endl;
    typ = 1;  mag = 0.; 
    totngal_++;
    return(totngal_);
	}
	//cout<<"idxmag="<<idxmag<<endl;
	
  //v_mag_[styp]=luminosity function of simulated galaxy type (styp)
  mag = 0.5*(v_mag_[styp](idxmag)+v_mag_[styp](idxmag-1));   typ = styp+1;
  //if(mag<-22) cout<<"mag<-22, idxmag="<<idxmag<<", type="<<typ<<endl;
  totngal_++;
  return(totngal_);
}

int GalFlxTypDist::GetGalMag(r_8& mag) const 
// Use this function if you just want to draw an absolute magnitude value 
{

	size_t nt=NTypes();
	if (nt>1)
		throw ParmError("ERROR! More than one LF added ... ");
	
  
  //SIMULATE MAG GIVEN TYPE
  sa_size_t maxidx = v_mag_[0].Size()-1; //equal to nbinmag i think (max possible index of bin)
  sa_size_t idxmag = (sa_size_t)(rg_.Flat01()*(double)maxidx)+1; // to randomly select a magnitude bin 
  
  // debugging
  if ((idxmag<1)||(idxmag>maxidx)) 
	{
    cout << " GalFlxTypDist::GetGalMag()/BUG!!! idxmag=" << idxmag << "  >=Max=" << v_mag_[0].Size() << endl;
    mag = 0.; 
    totngal_++;
    return(totngal_);
	}
	
  //v_mag_[0]=luminosity function
  mag = 0.5*(v_mag_[0](idxmag)+v_mag_[0](idxmag-1));   

  totngal_++;
  return(totngal_);
}

void GalFlxTypDist::AddSchechter(Schechter& t1,Schechter& t2,Schechter& t3)
{
	double phistar,Mstar,alpha;
	t1.GetParam(phistar,Mstar,alpha);
	t1_.SetParam(phistar,Mstar,alpha);
	t2.GetParam(phistar,Mstar,alpha);
	t2_.SetParam(phistar,Mstar,alpha);
	t3.GetParam(phistar,Mstar,alpha);
	t3_.SetParam(phistar,Mstar,alpha);
	
	SchTypes_=true;
}

void GalFlxTypDist::GetGalType(double mag,int& type)
{
	
	double ft1=t1_(mag);
	double ft2=t2_(mag);
	double ft3=t3_(mag);
	double sumf=ft1+ft2+ft3;
	
	double frac1=ft1/sumf;
	double frac2=ft2/sumf;
	double frac3=ft3/sumf;
	
	
	double sumfrac=frac1+frac2+frac3;
	double eps=1e-6;
	if (sumfrac>1+eps||sumfrac<1-eps)
		{
		cout <<" sumfrac = "<<sumfrac<<endl;
		throw ParmError("ERROR fractions do not sum to 1");
		}
		
	vector<double> fracs;
	fracs.push_back(frac1);
	fracs.push_back(frac2);
	fracs.push_back(frac3);
	
	//sort(fracs.begin(),fracs.end());
	
	double rdn=rg_.Flat01();
	
	if (rdn>=0&&rdn<frac1)
		type = 1;
	else if (rdn>=frac1&&rdn<(frac1+frac2))
		type = 2;
	else if (rdn>=(frac1+frac2)&&rdn<1)
		type = 3;

}


// PRINT OUT DIST PARAS
ostream& GalFlxTypDist::Print(ostream& os) const 
{
os << "----- GalFlxTypDist::Print() N_Types= " << NTypes() << " TotalTypeFrac=" <<  tot_frac_ << endl;
for(size_t k=0; k<NTypes(); k++) {
  os << " --Type=" << k+1 << " typeFraction=" << v_frac_[k] << " SumtypeFrac=" << v_sumfrac_[k] << endl;
  v_mag_[k].Print(os, v_mag_[k].Size());
  os << endl;
}
os << " -------------------------------------------- " << endl;
return os;
}

// PRINT DISTIRBUTION TO A FILE
void GalFlxTypDist::PrintDist(string fname)
{

size_t nt=NTypes();
int index_long=0;//will be index of type with longest length of v_mag_
//int sz=0;//start with size of 0
for(size_t i=0;i<nt;i++)
	if(v_mag_[i].Size()>index_long) index_long=i;

int max=v_mag_[index_long].Size();
for(sa_size_t j=0;j<max;j++)
	{
	cout <<v_mag_[0](j)<<"    "<<v_mag_[1](j)<<"    "<<v_mag_[2](j)<<endl;
	}
};

// MAKE Z DISTRIBUTION TO FOLLOW VOLUME ELEMENT
int GalFlxTypDist::GetZDist(ClassFunc& dVdz, double zmin,double zmax, 
                              int nbin, int nstep)
{
/****************** FIRST MAKE EQUI BIN DISTRIBUTION ******************/

  // first check limit values are sensible
  if ((zmin>=zmax)) { 
    cout << " GalFlxTypDist::GetZDist()/Error : Bad parameters (zmin>=zmax)" << endl;
    throw ParmError("GalFlxTypDist::GetZDist()/Error : Bad parameters");
	}
  
  // put nbin and nstep to sensible values if necessary
  if (nbin < 1)  nbin= 1;
  if (nstep<10*nbin) nstep = 10*nbin;
  
  
  // print out to help with debugging/checking
  if (prtlev_ > 1) 
    cout  << " --- GalFlxTypDist::GetZDist()/Info "//, GalType=" << v_frac_.size() 
          << " ( NBin=" << nbin << " NStep=" << nstep << ")" << endl;
		  
  // COMPUTE APPROXIMATELY THE INTEGRAL OF dVdz 
  // works out bin widths so distribution has equal number of galaxies in 
  // 1) make normalised dVdz function == probability distribution
  Vector vfm(nstep);  // dVdz distribution (dVdz function) vector
  double sum = 0.;
  double delz = (zmax-zmin)/(double)nstep; // z spacing given zmin,zmax,nstep: -1?
  double xm = zmin;
  for(sa_size_t kk=0; kk<vfm.Size(); kk++) {
    double fdm = dVdz(xm);// value of dVdz function defined by dVdz parameters at z xm
    sum += fdm; // adds up all dVdz values for each z between zmin and zmax 
    vfm(kk) = fdm; // dVdz function vector
    xm += delz; // step up z value for next loop
  }
  vfm /= (r_8)(sum);  // renormalize the dVdz distribution vector so it sums to 1
  // basically vfm is input dVdz function given by dVdz but defined according to zmin, zmax, delz and normalised so it sums to 1
  
  // print out to help with debugging/checking
  if (prtlev_ > 1) 
    cout << " ... Integral[dVdz Dist," << zmin << "," << zmax << "]=" << sum << endl;
         //<< " TypeFraction=" << fr << " CurrentTotalFrac=" << tot_frac_ << endl;
		 
  xm = zmin;
  double seuil = 1./nbin; // because vfm sums to 1 then if have equi pop bins then each one will sum to 1/nbin
  sum = 0.;
  dVdzbins_.SetSize(nbin+1, BaseArray::RowVector); //define a vector with nbin+1 rows, and specify that it's a row vector
  sa_size_t jm=0;
  dVdzbins_(jm) = zmin;   jm++;
  
  // 2) loop over normalised dVdz function vector
  for(sa_size_t kk=0; kk<vfm.Size(); kk++) 
  {
    if (jm > nbin) 
	{// doesn't print this: just a check
      cout << " GalFlxTypDist::GetZDist()/Warning jm=" << jm << "/nbin=" 
           << nbin << " seuil=" << seuil << endl;
      break;
	}
	
    xm += delz;  sum += vfm(kk); // sum will be getting larger each loop
	
	
    if (sum >= seuil) //decider if reached equipop bin yet
	{
	  // print out to help with debugging/checking
      if (prtlev_ > 2)  
        cout << " ... dVdzBin[" <<  jm << "] -> " << xm << " ( IntegdVdzDist=" 
             << seuil << " )" << endl;
	  // record mag value one LF has reached equi pop
      dVdzbins_(jm)= xm;  jm++; 
      seuil += (1./nbin);
    }
	
  }
  //dVdzbins_ is nbin+1 long after this
  
  // this is to fudge the end if needed I think
  while (jm < (nbin+1)) 
  {
  
    if (prtlev_ > 2)  
        cout << " ... Filling dVdzBin[" <<  jm << "] -> " << zmax << endl;
		
    dVdzbins_(jm)= zmax;  jm++; 
  }

return 1; 
  /****************** HAVE MADE EQUI BIN DISTRIBUTION ******************/
  };

// PICK Z VALUE FROM Z DISTRIBUTION ABOVE
int GalFlxTypDist::GetZ(double& z) const
{
  //SIMULATE DV/DZ 
  sa_size_t maxidx = dVdzbins_.Size()-1; //equal to nbin i think (max possible index of bin)
  sa_size_t idxmag = (sa_size_t)(rg_.Flat01()*(double)maxidx)+1; // to randomly select a z bin 
  
  // debugging
  if ((idxmag<1)||(idxmag>maxidx)) {
    cout << " GalFlxTypDist::GetZ()/BUG!!! idxmag=" << idxmag << "  >=Max=" << dVdzbins_.SizeX() << endl;
    z = 0.; 
    //totngal++;
    //return(totngal);
  }
	
  //
  z = 0.5*(dVdzbins_(idxmag)+dVdzbins_(idxmag-1));  
//  totngal_++;
//  return(totngal_);
  return 0;


};
