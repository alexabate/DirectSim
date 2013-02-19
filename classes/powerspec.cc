#include "powerspec.h"

PowerSpec::PowerSpec(TArray<r_8> galdens,double Dx,double Dy,double Dz)
: drho_(galdens)
// Reads in galaxy fluctuation grid and computes Fourier transform
// The Dkx, Dky, Dkz values are computed from the array size and Dx,Dy,Dz
{
	cout << "    PowerSpec::PowerSpec():"<<endl;
	
	Dx_=Dx, Dy_=Dy, Dz_=Dz;
	ComputeFourier(drho_,four_);
	SetDKDR();

	zc_=-1;
	
	
	cout << "    Printing Fourier grid info ...."<<endl;
	cout << "    Nx,Ny,Nz = "<<Nx_<<","<<Ny_<<","<<Nz_<<endl;
	cout << "    Dkx,Dky,Dkz = "<<Dkx_<<","<<Dky_<<","<<Dkz_<<endl;
	cout << "    cell size = dx,dy,dx = "<<Dx_<<","<<Dy_<<","<<Dz_<<endl;

	cout << "    Printing drho_ array info (should match Nx,Ny,Nz above) ...."<<endl;
	cout << "    SizeX,SizeY,SizeZ = "<<drho_.SizeX()<<","<<drho_.SizeY()<<","<<drho_.SizeZ()<<endl;

	cout << "    Printing four_ array info ...."<<endl;
	cout << "    SizeX,SizeY,SizeZ = "<<four_.SizeX()<<","<<four_.SizeY()<<","<<four_.SizeZ()<<endl;
	cout << "    EXIT PowerSpec::PowerSpec():"<<endl<<endl;
	
}

PowerSpec::PowerSpec(TArray<r_8> galdens,double Dx)
: drho_(galdens)
// Reads in galaxy fluctuation grid and computes Fourier transform
// The Dkx, Dky, Dkz values are computed from the array size and Dx
{
	cout << "    PowerSpec::PowerSpec():"<<endl;
	
	Dx_=Dx, Dy_=Dx, Dz_=Dx;
	ComputeFourier(drho_,four_);
	SetDKDR();

	zc_=-1;
	
	cout << "    Printing Fourier grid info ...."<<endl;
	cout << "    Nx,Ny,Nz = "<<Nx_<<","<<Ny_<<","<<Nz_<<endl;
	cout << "    Dkx,Dky,Dkz = "<<Dkx_<<","<<Dky_<<","<<Dkz_<<endl;
	cout << "    cell size = dx,dy,dx = "<<Dx_<<","<<Dy_<<","<<Dz_<<endl;

	cout << "    Printing drho_ array info (should match Nx,Ny,Nz above) ...."<<endl;
	cout << "    SizeX,SizeY,SizeZ = "<<drho_.SizeX()<<","<<drho_.SizeY()<<","<<drho_.SizeZ()<<endl;

	cout << "    Printing four_ array info ...."<<endl;
	cout << "    SizeX,SizeY,SizeZ = "<<four_.SizeX()<<","<<four_.SizeY()<<","<<four_.SizeZ()<<endl;
	cout << "    EXIT PowerSpec::PowerSpec():"<<endl<<endl;
	
}

void PowerSpec::ComputeFourier(TArray<r_8>& real_array, TArray< complex< r_8 > >& four_array)
{
// FFT galaxy flunction grid
	FFTWServer ffts;
	ffts.FFTForward(real_array,four_array);

	cout << "    PowerSpec::ComputeFourier(): "<<endl;
	cout << "    check size of arrays ..."<<endl;
	cout << "    size of drho = [sizex sizey sizez] = ["<<real_array.SizeX()<<"  "<<real_array.SizeY()<<"  "<<real_array.SizeZ()<<"]"<<endl; 
	cout << "    size of four = [sizex sizey sizez] = ["<<four_array.SizeX()<<"  "<<four_array.SizeY()<<"  "<<four_array.SizeZ()<<"]"<<endl;
	cout << "    EXIT PowerSpec::ComputeFourier(): "<<endl<<endl;
	return;
}

void PowerSpec::ComputeFourierBack(TArray< complex< r_8 > >& four_array, TArray<r_8>& real_array)
{
	// FFT convolved galaxy grid
	FFTWServer ffts;
	ffts.FFTBackward(four_array,real_array);

	cout << "    PowerSpec::ComputeFourierBack(): "<<endl;
	cout << "    check size of arrays ..."<<endl;
	cout << "    size of drho = [sizex sizey sizez] = ["<<real_array.SizeX()<<"  "<<real_array.SizeY()<<"  "<<real_array.SizeZ()<<"]"<<endl; 
	cout << "    size of four = [sizex sizey sizez] = ["<<four_array.SizeX()<<"  "<<four_array.SizeY()<<"  "<<four_array.SizeZ()<<"]"<<endl;
	cout << "    EXIT PowerSpec::ComputeFourierBack(): "<<endl<<endl;
	return;
}


void PowerSpec::SetDKDR()
//Read in info on FT grid
{
	Nx_=drho_.SizeX();// will be short dimension
	Ny_=drho_.SizeY();
	Nz_=drho_.SizeZ();
	Dkx_=(2*PI)/(Nx_*Dx_);
	Dky_=(2*PI)/(Ny_*Dy_);
	Dkz_=(2*PI)/(Nz_*Dz_);

	NRtot_ = Nx_*Ny_*Nz_; // number of pixels in the survey
	NCx_=Nx_/2+1; // think this should be the length of the short dimension in the four_ array
	// Kny = 2/Dx // Nyquist frequency of grid, should be sufficiently greater than klin
	// klin~0.1h/Mpc at z=0, ~0.2h/Mpc at z=1

	
};

void PowerSpec::CheckPZErr()
{
	if(Err_<0)
		throw ParmError("ERROR! Photo-z error < 0");
	if(Err_>0)	
		{
		cout <<"    Accounting for photometric redshift errors"<<endl;
		//maxk = coeff / Err; 
		//cout <<"    Coefficient = "<<coeff<<", error = "<<Err<<" Mpc"<<endl;
		cout <<"    Error = "<<Err_<<" Mpc"<<endl;
		cout <<"    Maximum radial k allowed = "<<maxk_<<endl;
		cout <<"    Minimum radial k component = "<<Dkz_<<endl;
		if (maxk_<Dkz_)
			{
			// is this a good idea? In this case should just keep kz = 0 wavenumber
			//cout <<"    Changing maximum radial k allowed from "<<maxk_<<" to "<<Dkz_+1e-4<<" or power spectrum will be zero"<<endl;
			//maxk_=Dkz_+1e-4;
			cout <<"    Warning: not keeping any radial k except 0th"<<endl;
			}
		}
	else
		cout <<"    Redshift error=0"<<endl;

};

void PowerSpec::BlowUpCheck()
// If undamping power spectrum we divide by:
// exp((kz*Err)^2)
// if this number is small (<tol_corr) power spectrum will -> inf
{
	double corr_at_maxk = exp(-(maxk_*Err_)*(maxk_*Err_));
	double corr_at_k1 = exp(-(Dkz_*Err_)*(Dkz_*Err_));

	if ( (corr_at_maxk<tol_corr_)&&(Err_>0&&undamp_) )
		{
		cout <<endl<<"**************** WARNING *****************"<<endl;
		cout <<"! Keeping too many radial k's            "<<endl;
		cout <<"! Power spectrum will probably BLOW UP!  "<<endl;
		cout <<"! Correction at max radial k = "<<corr_at_maxk<<endl;
		cout <<"************** END WARNING ***************"<<endl<<endl;
		}
	if ( (corr_at_k1<tol_corr_)&&(Err_>0&&undamp_) )
		{
		cout <<"************************************** WARNING **************************************"<<endl;
		cout <<"! Array sampling in radial direction is TOO LARGE, Dkz = "<<Dkz_<<endl;
		cout <<"! Power spectrum will probably BLOW UP! "<<endl;
		cout <<"! Correction at 1st radial k wavenumber = "<<corr_at_k1<<endl;
		double resk = sqrt(-log(tol_corr_))/Err_;
		cout <<"! Increase array radial direction length from "<<Nz_<<" pixels to roughly "<<(2*PI)/(resk*Dz_) <<" pixels"<<endl; 
		cout <<"************************************ END WARNING ************************************"<<endl<<endl;
		}


}


double PowerSpec::AccumulatePowerSpectra(HProf& hp, bool pixcor, double maxk, double Err, bool undamp, double tol_corr, double snoise)
// Compute spectrum from "four_" and fill profile histogram "hp"
// (Defaults: Err=0, coeff = 1, bool undamp=true, pixcor=true, snoise=0);
// Power spectrum in hp is NOT normalised
// Radial direction is assumed to be the third dimenion (z-coord)
// four_ : in the format: four_(nkx,nky,nkz)
// FIRST dimension is short dimension
// If photometric redshifts (Err>0) only Fourier modes (kx,ky,kz) with a radial wavescale 
// |kz|< kz_max ~= 1/Err are included in the analysis (Default: maxk=1000)
// if undamp = true divides out damping: 1/exp(kz^2sigz^2)
// Err is the standard deviation of the photo-z errors in comoving coordinates
// tol_corr: if divide Fourier coefficient by number greater than this power spectrum will blow up
 { 
	cout << "    PowerSpec::AccumulatePowerSpectra():"<<endl;
	maxk_=maxk;
	Err_=Err;
	tol_corr_=tol_corr;
	undamp_=undamp;

	if(hp.NBins()<0) 
		throw ParmError("ERROR! HProf bins undefined");
	
	hp.Zero();
	
	// SHOT NOISE (NOTHING IS DONE WITH THIS AT THIS TIME)
	if(snoise<=0.) snoise = 0.;
	double snoisesq = snoise*snoise / (double)NRtot_;
  
	// PHOTO-Z ERROR
	CheckPZErr();
		
	// Blowup Check	
	BlowUpCheck();
		
	// PIXEL CORRECTION
	//TVector<r_8> vfx(NCx_);
	if(pixcor)   // kz = l*Dkz_
		{ 
		cout <<"    Pixel correction is ON"<<endl;
	//	for(long ix=0;ix<NCx_;ix++) 
	//		{vfx(ix)=pixelfilter(ix*Dkx_*Dx_/2); vfx(ix)*=vfx(ix);}
		}
	else cout <<"    Pixel correction is OFF"<<endl; 
	
	double fx, fy, fz;// filter correction
	double sum=0;
	int nkeep=0;
	// wavenumber k = n * 2pi/L where n=wave index and L = length of grid side

	for(sa_size_t iz=0; iz<four_.SizeZ(); iz++)// assumed RADIAL direction
		{
		// get wave number 3RD/RADIAL dim
		double kz = iz;// wave index for +ve freq
		if (iz > four_.SizeZ()/2) 
			kz = four_.SizeZ()-iz;// wave index for -ve freq
		kz *= Dkz_; // wave number
		
		// if pixel filtering
		if(pixcor)
			{
			fz = pixelfilter(kz*Dz_/2);
			fz *= fz;// squared
			}
		else 
			fz=1;
		
		for(sa_size_t iy=0; iy<four_.SizeY(); iy++) 
			{
			// get wave number 2ND dim
			double ky = iy;// wave index for +ve freq
			if (iy > four_.SizeY()/2) 
				ky = four_.SizeY()-iy;// wave index for -ve freq
			ky *= Dky_; // wave number

			// if pixel filtering
			if(pixcor)
				{
				fy = pixelfilter(ky*Dy_/2);
				fy *= fy;// squared
				}
			else fy=1;
			
			for(sa_size_t ix=0; ix<four_.SizeX(); ix++) // THIS IS THE SHORT DIMENSION, no neg freq
				{
				// get wave number 1ST dim (straightforward: no neg freq)
				double kx = ix*Dkx_;

				// if pixel filtering
				if(pixcor)
					{
					fx = pixelfilter(kx*Dx_/2);
					fx *= fx;// squared
					}
				else 
					fx=1;

				// k modulus
				double kmod = sqrt((kx*kx+ky*ky+kz*kz));

				// Fourier component
				complex< r_8 > za = four_(ix, iy, iz);

				// Fourier component * its complex conjugate
				double pk = za.real()*za.real()+za.imag()*za.imag();
				
				//fx = (pixcor) ? vfx(ix): 1.;
				double f = fx*fy*fz;
				
				if(Err_>0&&undamp_)// if want to undamp	
					{
					double g = exp(-(kz*Err_)*(kz*Err_));
					//if ((ix!=0)&&(iy!=0)&&(iz!=0)&&(kz<=maxk_))
					if (kz<=maxk_)
						{hp.Add(kmod, pk/(f*g)); sum+=pk; nkeep++;}
					}
				else // no undamp
					//if ((ix!=0)&&(iy!=0)&&(iz!=0)&&(kz<=maxk_))  
					if (kz<=maxk_)
						{hp.Add(kmod, pk/f); sum+=pk; nkeep++;}

				}// end loop over X dim
			}// end loop over Y dim
		}// end loop over Z dim

if(Err>0)
	cout <<"    Number of wavevectors="<<four_.Size()<<", N kept="<<nkeep<<", therefore N not included: "<<four_.Size()-nkeep<<endl;
 cout <<"    Sum of Fourier coefficients sq ="<<sum<<endl;

 return sum;
 cout << "    EXIT PowerSpec::AccumulatePowerSpectra():"<<endl<<endl;

};

void PowerSpec::WritePS(string fname,HProf& Pdata,r_4 Voldata,string PSimlssfile,double meandens)
{

	cout << "    PowerSpec::WritePS()"<<endl;
	// Read in SimLSS power spectra into an array
	ifstream ifs(PSimlssfile.c_str());
	TArray<r_8> SimLSSPSfile;
	sa_size_t nr, nc;
	SimLSSPSfile.ReadASCII(ifs,nr,nc);
	int icolk=0,icols=1,icolf=2,icolv=3;
	double Volsimlss = SimLSSPSfile(icolv,0);
	cout << "    Volume of SimLSS cube = "<<Volsimlss<<endl;
  	cout << "    Read in SimLSS power spectra, file has "<<nr<<" rows and "<<nc<<" columns"<<endl;
	cout << "    SimLSSPSfile array has "<<SimLSSPSfile.SizeX()<<" rows and "<<SimLSSPSfile.SizeY()<<" columns"<<endl;
	
	// Make an interpolation table
	
	vector<double> kvals,pvals1,pvals2;
	for(int kk=0; kk<nr; kk++) 
		{
		double kv=SimLSSPSfile(icolk,kk);
		double p1=SimLSSPSfile(icols,kk);
		double p2=SimLSSPSfile(icolf,kk);
		kvals.push_back(kv);
		pvals1.push_back(p1);
		pvals2.push_back(p2);
		//cout << kvals[kk]<<endl;
		}
	double kmin=kvals[0],kmax=kvals[nr-1];
	cout <<"    K range is "<<kmin<<"<k<"<<kmax<<endl;
	
	int_8 nk=1000000;
	SInterp1D simlss,simlssf; 
	simlss.DefinePoints(kvals,pvals1,kmin,kmax,nk);
	simlssf.DefinePoints(kvals,pvals2,kmin,kmax,nk);
	
	cout << endl;

	// Data power spectrum
	Histo Pdatah=Pdata.GetHisto();
	// number of k bins	
	int_4 nbk=Pdatah.NBins();
	cout << "    Number of bins in data power spectrum = "<<nbk<<endl;

	//read histogram values into a file
	ifstream inp;
	ofstream outp;
	inp.open(fname.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
		{
		inp.clear(ios::failbit);
		cout << "Writing to file ..." << fname.c_str() << endl<< endl<< endl;
		outp.open(fname.c_str(), ofstream::out);
		outp << "# Volume of data = "<<Voldata<<", Volume of SimLSS = "<<Volsimlss<<", mean density of fudged simlss grid = "<<meandens<<endl;
		outp << "# Nx,Ny,Nz,R(Mpc),zc = "<<Nx_<<","<<Ny_<<","<<Nz_<<","<<Dx_<<","<<zc_<<endl;
		for(int_4 i=0;i<nbk;i++)
			{
			r_8 kv = Pdatah.BinCenter(i);
			r_8 Praw = Pdatah.operator()(i);
			r_8 Pslss = simlss(kv);
			r_8 Pslssf = simlssf(kv);

			r_8 Pnorm_uncorr = Praw*Voldata*(1+meandens)*(1+meandens);
			r_8 Pnorm = Pnorm_uncorr*(Pslss/Pslssf);
			

			outp <<kv<<"    "<<Pnorm<<"    "<<Pnorm_uncorr<<"    "<<Praw*Voldata<<"    "<<Pslss*Volsimlss<<"    "<<Pslssf*Volsimlss<<endl;
			}
		outp.close();
		}
	else
	{
	cout << "Error...file """ << fname.c_str() << """ exists" << endl;
	}
	cout << "    EXIT PowerSpec::WritePS()"<<endl;
};

void PowerSpec::WritePS(string fname,HProf& Pdata,r_4 Voldata,HProf& PSimlss,HProf& PSimlssf,r_4 Volsimlss,double meandens)
// Writes 4 power spectra to a text file
// the format is:
// [k values] [normalised galaxy PS] [normalised-uncorrected galaxy PS] [raw galaxy PS] [simlss PS] [simlss PS after setting delta<-1 ->-1]
// To get undistorted galaxy power spectrum P_g(k) [normalised galaxy PS]:
// P_g(k) = P(k)*VolCat*f(k)/(1+meandens)^2
// where f(k) is found by dividing the SimLSS PS by the SimLSS PS when delta<-1 ->-1.
// Additional corrections will need to be made for shot noise, photo-z errors
{
	cout << "    PowerSpec::WritePS()"<<endl;
	Histo Pdatah=Pdata.GetHisto();
	Histo PSimlssh =PSimlss.GetHisto();
	Histo PSimlssfh=PSimlssf.GetHisto();

	// number of k bins	
	int_4 nbk=Pdatah.NBins();
	
	// check k bins are same for Simlss
	int_4 nbs= PSimlssh.NBins();
	int_4 nbf=PSimlssfh.NBins();
	double eps=1e-6;
	if(abs(nbk-nbs)>eps||abs(nbk-nbf)>eps)
		throw ParmError("ERROR! bin sizes are different between data and SimLSS");
	
	cout << "    Number of bins="<<nbk<<endl;
	
	//read histogram values into a file
	ifstream inp;
	ofstream outp;
	inp.open(fname.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
		{
		inp.clear(ios::failbit);
		cout << "Writing to file ..." << fname.c_str() << endl<< endl<< endl;
		outp.open(fname.c_str(), ofstream::out);
		outp << "# Volume of data = "<<Voldata<<", Volume of SimLSS = "<<Volsimlss<<", mean density of fudged simlss grid = "<<meandens<<endl;
		outp << "# Nx,Ny,Nz,R(Mpc),zc = "<<Nx_<<","<<Ny_<<","<<Nz_<<","<<Dx_<<","<<zc_<<endl;
		for(int_4 i=0;i<nbk;i++)
			{
			r_8 kvals=Pdatah.BinCenter(i);
			r_8 Praw=Pdatah.operator()(i);
			r_8 kv2= PSimlssh.BinCenter(i);
			r_8 Pslss= PSimlssh.operator()(i);
			r_8 kv3=PSimlssfh.BinCenter(i);
			r_8 Pslssf=PSimlssfh.operator()(i);

			if(abs(kvals-kv2)>eps||abs(kvals-kv3)>eps)
				throw ParmError("ERROR! k bins are different between data and SimLSS");

			r_8 Pnorm_uncorr = Praw*Voldata*(1+meandens)*(1+meandens);
			r_8 Pnorm = Pnorm_uncorr*(Pslss/Pslssf);
			

			outp <<kvals<<"    "<<Pnorm<<"    "<<Pnorm_uncorr<<"    "<<Praw*Voldata<<"    "<<Pslss*Volsimlss<<"    "<<Pslssf*Volsimlss<<endl;
			}
		outp.close();
		}
	else
	{
	cout << "Error...file """ << fname.c_str() << """ exists" << endl;
	}
	cout << "    EXIT PowerSpec::WritePS()"<<endl;
};

void PowerSpec::Write1PS(string fname, HProf& hpgals, double VolCat)
// writes the power spectrum to a text file
// the format is:
// [k values] [galaxy PS] [volume of galaxy catalog]
{
	cout << "    PowerSpec::Write1PS()"<<endl;
	Histo histogals=hpgals.GetHisto();
	
	int_4 nbg=histogals.NBins();
	cout << "    Number of bins="<<nbg<<endl;
	
	//read histogram values into a file
	ifstream inp;
	ofstream outp;
	inp.open(fname.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
		{
		inp.clear(ios::failbit);
		cout << "Writing to file (1)..." << fname.c_str() << endl<< endl<< endl;
		outp.open(fname.c_str(), ofstream::out);
		for(int_4 i=0;i<nbg;i++)
			{
			r_8 bcg=histogals.BinCenter(i);
			r_8 bvg=histogals.operator()(i);

			
			//r_8 be=histogals.Error(i);
			outp << bcg<<"      "<<bvg<<"      "<<VolCat<<endl;
			}
		outp.close();
		}
	else
	{
	cout << "Error...file """ << fname.c_str() << """ exists" << endl;
	}
	cout << "    EXIT PowerSpec::Write1PS()"<<endl;
}

void PowerSpec::Write2PS(string fname, HProf& hp1, HProf& hp2, double VolCat,double md,double mdf)
// writes the power spectrum to a text file
// the format is:
// [k values] [PS] [PS] [volume of galaxy catalog]
{
	cout << "    PowerSpec::Write2PS()"<<endl;
	Histo histo1=hp1.GetHisto();
	Histo histo2=hp2.GetHisto();

	int_4 nbg=histo1.NBins();
	int_4 nbg2=histo2.NBins();
	if (nbg!=nbg2)
		throw ParmError("Power spectra have different number of bins");

	cout << "    Number of bins="<<nbg<<endl;
	
	//read histogram values into a file
	ifstream inp;
	ofstream outp;
	inp.open(fname.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
		{
		inp.clear(ios::failbit);
		cout << "Writing to file ..." << fname.c_str() << endl<< endl<< endl;
		outp.open(fname.c_str(), ofstream::out);
		outp <<"# mean density of delta cube = "<<md<<", mean density of fudged";
 		outp <<" delta cube = "<<mdf<<endl;
		outp <<" cols: k values - PS_1 - PS_2 - Volume"<<endl;
		for(int_4 i=0;i<nbg;i++)
			{
			r_8 bcg=histo1.BinCenter(i);
			r_8 bvg=histo1.operator()(i);

			//r_8 bcg2=histo2.BinCenter(i);
			r_8 bvg2=histo2.operator()(i);

			outp << bcg<<"      "<<bvg<<"      "<<bvg2<<"      "<<VolCat<<endl;
			}
		outp.close();
		}
	else
	{
	cout << "Error...file """ << fname.c_str() << """ exists" << endl;
	}
	cout << "    EXIT PowerSpec::Write2PS()"<<endl;
}

//************** Stuff  below here might not work or be useful **************//

void PowerSpec::FourierSpaceWindow(HProf& hp)
// Compute spectrum from "four_" and fill profile histogram "hp"
// (Defaults: Err=0, coeff = 1, bool undamp=true, pixcor=true, snoise=0);
// Power spectrum in hp is NOT normalised
// Radial direction is the third dimenion (z-coord)
// four_ : in the format: four_(nkx,nky,nkz)
// FIRST dimension is short dimension
// If photometric redshifts (Err>0) only Fourier modes (kx,ky,kz) with a radial wavescale 
// |kz|< kz_max = coeff/Err are included in the analysis (Default: coeff=1)
// if undamp = true divides out damping: 1/exp(kz^2sigz^2)
// Err is the standard deviation of the photo-z errors in comoving coordinates
 { 
	cout << "    PowerSpec::FourierSpaceWindow():"<<endl;
	if(hp.NBins()<0) 
		throw ParmError("ERROR! HProf bins undefined");
	
	hp.Zero();
		
	// PIXEL CORRECTION
	TVector<r_8> vfx(NCx_);
	for(long ix=0;ix<NCx_;ix++) 
			{vfx(ix)=pixelfilter(ix*Dkx_*Dx_/2); vfx(ix)*=vfx(ix);}
	
	double fx, fy, fz, sum=0;
	// wavenumber k = n * 2pi/L where n=wave index and L = length of grid side
	for(sa_size_t iz=0; iz<four_.SizeZ(); iz++) 
		{
		double kz = iz;// wave index for +ve freq
		if (iz > four_.SizeZ()/2) 
			kz = four_.SizeZ()-iz;// wave index for -ve freq
		kz *= Dkz_; // wave number
		fz = pixelfilter(kz*Dz_/2);

		//cout <<"kz="<<kz<<endl;
		for(sa_size_t iy=0; iy<four_.SizeY(); iy++) 
			{
			double ky = iy;
			if (iy > four_.SizeY()/2) 
				ky = four_.SizeY()-iy;
			ky *= Dky_;
			fy = pixelfilter(ky*Dy_/2);

			for(sa_size_t ix=0; ix<four_.SizeX(); ix++) // THIS IS THE SHORT DIMENSION
				{
				double kx = ix*Dkx_;
				double kmod = sqrt((kx*kx+ky*ky+kz*kz));
				complex< r_8 > za = four_(ix, iy, iz);
				
				// modulus |W(k)|^2:
				double Wk = za.real()*za.real()+za.imag()*za.imag();
				
				fx = vfx(ix);
				double f = fx*fx*fy*fy*fz*fz;
				
				if ((ix!=0)&&(iy!=0)&&(iz!=0))
					{
					hp.Add(kmod, Wk/f); 
					sum+=Wk;
					}
	
				}
			}
		}
		
		double FourierVolEl = Dkx_*Dky_*Dkz_;
		double int_Wmodsq = sum/(2*PI*2*PI*2*PI) * FourierVolEl;
		cout << " d3k = "<<FourierVolEl<<", sum="<<sum<<", int_Wmodsq="<<int_Wmodsq<<endl;


 cout << "    EXIT PowerSpec::FourierSpaceWindow():"<<endl<<endl;

};

void PowerSpec::Conv1DGaussFT(double Gsig, double Gmean)
// Convolve Fourier Transform four_ with a 1D Gaussian mean Gmean and
// std Gsig.  Gaussian is in z-dimension direction
{
	cout << "    PowerSpec::Conv1DGaussFT():"<<endl;
	
	int ndim=3;
	sa_size_t mydim[ndim];
	mydim[0]=Nx_=four_.SizeX(); mydim[1]=four_.SizeY(); mydim[2]=four_.SizeZ();
	fourconv_.SetSize(ndim,mydim); 
		
	cout << "    Convolve Fourier space distribution ... "<<endl;
	// wavenumber k = n * 2pi/L where n=wave index and L = length of grid side
	for(sa_size_t iz=0; iz<four_.SizeZ(); iz++) 
		{
		double kz = iz;// wave index for +ve freq
		if (iz > four_.SizeZ()/2) 
			kz = four_.SizeZ()-iz;// wave index for -ve freq
		kz *= Dkz_; // wave number z dim
		
		for(sa_size_t iy=0; iy<four_.SizeY(); iy++) 
			{
			double ky = iy;
			if (iy > four_.SizeY()/2) 
				ky = four_.SizeY()-iy;
			ky *= Dky_; // wave number y dim
			for(sa_size_t ix=0; ix<four_.SizeX(); ix++) // THIS IS THE SHORT DIMENSION
				{
				double kx = ix*Dkx_; // wave number x dim
				double kmod = sqrt((kx*kx+ky*ky+kz*kz));
				complex< r_8 > za = four_(ix, iy, iz);
				
				double g = exp(-(kz*Gsig)*(kz*Gsig));
				fourconv_(ix, iy, iz) = four_(ix, iy, iz)*g;
			
				}
			}
		}

	cout << "    FT back to convolved real space distribution ... "<<endl;
	ComputeFourierBack(fourconv_,drhoconv_);

	cout << "    EXIT PowerSpec::Conv1DGaussFT():"<<endl<<endl;

};

void PowerSpec::GetObsPS(HProf& hpgals, r_4 VolCat, HProf& hpsim, HProf& hpsimf,r_4 VolSim)
// put correctly normalised observed Power spectrum into a vector
// should double check k bins from each power spectrum are the same
{
	cout << "    PowerSpec::GetObsPS()"<<endl;
	Histo histogals=hpgals.GetHisto();
	Histo histosim = hpsim.GetHisto();
	Histo histosimf=hpsimf.GetHisto();
	
	int_4 nbg=histogals.NBins();
	int_4 nbs= histosim.NBins();
	int_4 nbf=histosimf.NBins();
	double eps=1e-3;
	if(abs(nbg-nbs)>eps||abs(nbg-nbf)>eps)
		{ 
		cout << "    ERROR! bin sizes are different"<<endl;
		return;
		}
	
	kobs_.SetSize(nbg);
	Pobs_.SetSize(nbg);
		
	for(int_4 i=0;i<nbg;i++)
		{
		r_8 bcg=histogals.BinCenter(i);
		r_8 bvg=histogals.operator()(i);
		r_8 bvs= histosim.operator()(i);
		r_8 bvf=histosimf.operator()(i);
		
		kobs_(i) = bcg;
		Pobs_(i) = bvg*(bvs/bvf)*VolCat;
		
		}
	cout << "    EXIT PowerSpec::GetObsPS()"<<endl;
};




