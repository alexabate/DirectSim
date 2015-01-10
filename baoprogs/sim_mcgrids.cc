/**
  * @file  sim_mcgrids.cc
  * @brief simulate galaxy grids for computing the"multiplicative correction A(k)"
  *        and the "additive correction B(k)" to the power spectrum.
  *
  * @todo check this works!
  *
  * @author Alex Abate
  * Contact: abate@email.arizona.edu
  *
  */

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>
#include <matharr.h>

#include "machdefs.h"
#include "sopnamsp.h"
#include "timing.h"
#include "array.h"
#include "hisprof.h"
#include "histerr.h"
#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "mydefrg.h"
#include "swfitsdtable.h"
#include "resusage.h"


#include "geneutils.h"
#include "cat2grid.h"
#include "powerspec.h"
#include "mass2gal.h"
#include "schechter.h"
#include "selectfunc.h"
#include "cosmocalcs.h"
#include "gftdist.h"



void usage(void);
void usage(void) {

	cout << endl <<" Usage: sim_mcgrids [...options...]" << endl<<endl;
	
	cout << "  Program to simulate galaxy grids for computing the "<<endl;
	cout << "  'multiplicative correction A(k)' and the 'additive correction B(k)' "<<endl;
	cout << "  to the power spectrum."<<endl;
	cout << "  See section 5 of Blake et al 2007 MNRAS 374, 4, 1527-1548, "<<endl;
	cout << "  arXiv:astro-ph/0605303 for an explanation"<<endl<<endl;
			
	cout << "  The program takes the following inputs: "<<endl;
	cout << "  i)   SimLSS over-density distribution "<<endl;
	cout << "  ii)  Selection function "<<endl;
	cout << "  iii) Sub-array properties: Nx, Ny, Nz, R "<<endl;
	cout << "  iv)  Photometric redshift convolution function "<<endl;
	cout << "  v)   *_subinfo.txt file "<<endl;
	cout << "  vi)  density of random grid "<<endl;
	cout << "  vii) N Monte Carlo realisations (=1 for now) "<<endl<<endl;
			
	cout << "  The program computes the 'raw' power spectra of the simulated"<<endl;
	cout << "  galaxy grid and the random galaxy grid."<<endl;
	cout << "  (Eventually these will be averaged over N realisations)"<<endl<<endl;
			
	cout << " -C : MassDistFileName: FITS filename containing SimLSS output"<<endl;
	cout << " -O : OutputPS: Root filename to output power spectra to"<<endl;
	cout << " -d : dens: density of random grid (default=1)"<<endl;
	cout << " -a : auto initialise random number seed"<<endl;
	cout << " -e : PZCFunc: Photometric redshift convolution function"<<endl;
	cout << " -p : Nx,Ny,Nz,R: sub-array properties"<<endl;
	cout << " -s : SFunc: Survey selection function"<<endl;
	cout << " -x : xplanes: Number of SimLSS cube planes to remove (default 1)"<<endl;
	cout << " -z : SubInfoFile: *_subinfo.txt file"<<endl;
	cout << endl;
	
};


int main(int narg, char* arg[])
{
	cout << " ==== sim_mcgrids.cc program  ==== " << endl;
	// Make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;


	string infile;
	string outfileroot; 
	string SFunc;
	string subinfo;
	string pzconv;
	//double pzsig=0.03; // photometric redshift convolution function currently just Gaussian
	int meandens = 1; // must be integer
	long Nx=10,Ny=10,Nz=10;
	double Res=5;
	double zref=1.;
	double Dcref=1000;
	vector<int> idv;
	
	bool DEBUG = false;
	string debug;
	
	RandomGenerator rg;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"haC:O:d:e:p:s:z:N:")) != -1) {
	    switch (c) {
	    case 'a' :
		    rg.AutoInit(0);
		    break;
	    case 'C' :
		    infile = optarg;
		    break;
	    case 'O' :
		    outfileroot	= optarg;
		    break;
	    case 'd' :
		    sscanf(optarg,"%ld",&meandens);
		    break;
	    case 'e' :
		    pzconv = optarg;
		    break;
	    case 'p' : // CURRENTLY REDUNDANT BUT WON'T BE WHEN SIMLSS IS INSIDE
	 	    sscanf(optarg,"%ld,%ld,%ld,%lf",&Nx,&Ny,&Nz,&Res);
		    break;
	    case 's' :
		    SFunc = optarg;
		    break;
	    case 'z' :
		    subinfo = optarg;
		    break;
	    case 'N' :
		    debug = optarg;
		    DEBUG=true;
		    break;
	    case 'h' :
		    default :
		    usage(); return -1;
		   }
	    }


	cout << "     Reading command line arguments ... "<<endl;
	cout << "     SimLSS file is "<<infile<<endl;
	cout << "     Removing "<<xplanes<<" planes from SimLSS cube"<<endl;
	cout << "     Sub-array properties are Nx,Ny,Nz,R="<<Nx<<","<<Ny<<","<<Nz<<","<<Res<<endl;
	cout << "     Selection function will be read from "<<SFunc<<endl;
	//cout << "     Photometric convolution function is Gaussian with sigma = "<<pzsig<<endl;
	cout << "     Generating unclustered distribution with mean = "<<meandens<<endl;
	cout << "     Output power spectra will be written to files starting "<<outfileroot<<endl;
	cout << "     ... End reading command line arguments"<<endl<<endl;
    //-- end command line arguments
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
	ResourceUsage res;
	
	
	// READ IN SUBINFO FILE
	cout <<"     Read in subinfo file "<< subinfo <<endl;
	ifstream ifs(subinfo.c_str());
	Array B;
	sa_size_t nr, nc;
	B.ReadASCII(ifs,nr,nc);
	double minze = B(0,0);
	double zce = B(1,0);
	double maxze = B(2,0);
	cout <<"     Sub-array bounds: "<< minze <<"< z <"<< maxze <<endl;
	cout <<"     Sub-array center: zc = "<< zce <<endl;
	cout << endl;
	
	
	// READ SELECTION FUNCTION FILE
	cout << "    Read selection function file"<<endl;
	ifstream inp;
	inp.open(SFunc.c_str(), ifstream::in);
	if(inp.fail())
		throw ParmError("ERROR! Selection function file does not exist");
	ComputedSelFunc* sfp=new ComputedSelFunc(SFunc);
	cout << "     Check selection function at sub-array edges and center"<<endl;
	SelectionFunctionInterface& selfunc=*sfp;
	cout << "     SF(minz) = "<< selfunc(minze) <<", SF(zc) = "<< selfunc(zce);
	cout << ", SF(maxz) = "<< selfunc(maxze) <<endl;
	cout << endl;
	
	
	// READ PHOTO-Z CONVOLUTION FILE
	// also want to interpolate this
	cout << "     Read photo-z convolution file"<<endl;
	SInterp1D pzcfi;
	pzcfi.ReadXYFromFile(pzconv);
	cout << "     Check photo-z convolution function at sub-array edges and center"<<endl;
	cout << "     PZ(zc-minz) = "<< pzcfi(zce-minze) <<", PZ(zc-zc) = ";
	cout << pzcfi(zce-zce) <<", PZ(zc-maxz) = "<< pzcfi(zce-maxze) <<endl;
	double dzmin = pzcfi.XMin(), dzmax = pzcfi.XMax();		
	cout << "     Check photo-z convolution function integration "<<endl;
	double IntDist = IntegrateFunc(pzcfi,dzmin,dzmax);
	double IntDist2 = IntegrateFunc(pzcfi,-0.1,0.1);
	cout << "     Fraction of photo-z convolution function between -0.1<dz<0.1 = ";
	cout << IntDist2/IntDist <<endl;
	cout << endl;
	

	// INITIALISE COSMO
	cout << "     Initialise cosmology: (same as SimLSS)"<<endl;
	double h=0.71, OmegaM=0.267804, OmegaL=0.73;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout << "     OmegaK="<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
	cout << ", OmegaL="<< su.OmegaLambda() <<", OmegaB="<< su.OmegaBaryon();
	cout << ", H0="<< su.H0() <<endl;
	
	
	// READ INPUT SIMLSS FILE
	// later will change this to generate SimLSS cube here (with -p options)
	cout << "     Reading SimLSS input file= " << infile << endl;  
	FitsInOutFile fin(infile,FitsInOutFile::Fits_RO);
	TArray<r_8> dens;
	fin >> drho;
	/*TArray<r_8> drho;
	fin >> drho;
	cout << drho.Info();
	cout << "    Print original drho array size: "<<drho.SizeX()<<"x"<<drho.SizeY()<<"x"<<drho.SizeZ()<<endl<<endl<<endl;*/
	  

	cout <<"    Initialise Mass2Gal: remove planes, clean" << endl;
	Mass2Gal m2g(fin,su,rg);
	double mean, sig;
	TArray<r_8> mass;
	m2g.ODensArray(mass);
	MeanSigma(mass, mean, sig);
	cout << endl <<"    RAW DENS CUBE STATS: Mean=" << mean << " Sigma=" << sig;
	cout << endl << endl;
	long NX=m2g.ReturnNX();
	long NY=m2g.ReturnNY();
	long NZ=m2g.ReturnNZ();
	zref=m2g.ReturnZref(); // CHECK Z REF IS SAME AS ZC IN ZEDGES FILE
	double gridres = m2g.ReturnDX();
	cout << "    zref = "<< zref <<", zc = "<< zce <<endl;
	cout <<endl;
	
	
	// SIMULATE N GAL PER PIXEL
	
	cout << "    GENERATING CLUSTERED AND UNCLUSTERED GALAXY DISTRIBUTIONS"<<endl;

    // ----------------------- subsitute in LFParameters -----------------------
	cout << "    Convert rho/rho^bar To Mean NGal"<<endl;  
	cout << "    Set up Schechter functions for each galaxy type"<<endl;
	cout <<" ... GOODS B band: Early types, Late types, Starbursts"<<endl;
	cout <<" ... see Table 3 in Dahlen et al 2005"<<endl;
	string LFplace;
	char * plf=getenv("SIMBAOLF");
	if (plf==NULL) {
		cout <<"ERROR LF LOCATION ENVIRONMENT VARIABLE NOT DEFINED"<<endl;
		return 1;
		}
	else {
		LFplace=plf;
		cout <<"    Location of LF file is "<< LFplace <<endl;
		}
	string LFfile = LFplace +	"GOODS_B_LF.txt";
	ifstream ifs2;
	ifs2.open(LFfile.c_str(), ifstream::in);
	if (ifs2.fail())
		cout <<"  ERROR: failed to find luminosity function file"<<endl;

	TArray<r_4> LFTable;
	sa_size_t nr2, nc2;
	LFTable.ReadASCII(ifs2,nr2,nc2);

	
	int MstarCol=2, AlphaCol=3, PhiStarCol=4;
	// ALL GALAXIES
	double MstarAz1=LFTable(MstarCol,13),alpAz1=LFTable(AlphaCol,13),phistarAz1=LFTable(PhiStarCol,13)*1e-4;
	double MstarAz2=LFTable(MstarCol,14),alpAz2=LFTable(AlphaCol,14),phistarAz2=LFTable(PhiStarCol,14)*1e-4;
	double MstarAz3=LFTable(MstarCol,15),alpAz3=LFTable(AlphaCol,15),phistarAz3=LFTable(PhiStarCol,15)*1e-4;
	
	// EARLY TYPES
	double MstarEz1=LFTable(MstarCol,1),alpEz1=LFTable(AlphaCol,1),phistarEz1=LFTable(PhiStarCol,1)*1e-4;
	double MstarEz2=LFTable(MstarCol,6),alpEz2=LFTable(AlphaCol,6),phistarEz2=LFTable(PhiStarCol,6)*1e-4;
	double MstarEz3=LFTable(MstarCol,10),alpEz3=LFTable(AlphaCol,10),phistarEz3=LFTable(PhiStarCol,10)*1e-4;
	
	// LATE TYPES
	double MstarLz1=LFTable(MstarCol,3),alpLz1=LFTable(AlphaCol,3),phistarLz1=LFTable(PhiStarCol,3)*1e-4;
	double MstarLz2=LFTable(MstarCol,7),alpLz2=LFTable(AlphaCol,7),phistarLz2=LFTable(PhiStarCol,7)*1e-4;
	double MstarLz3=LFTable(MstarCol,11),alpLz3=LFTable(AlphaCol,11),phistarLz3=LFTable(PhiStarCol,11)*1e-4;
	
	// STARBURST TYPES
	double MstarSz1=LFTable(MstarCol,4),alpSz1=LFTable(AlphaCol,4),phistarSz1=LFTable(PhiStarCol,4)*1e-4;
	double MstarSz2=LFTable(MstarCol,8),alpSz2=LFTable(AlphaCol,8),phistarSz2=LFTable(PhiStarCol,8)*1e-4;
	double MstarSz3=LFTable(MstarCol,12),alpSz3=LFTable(AlphaCol,12),phistarSz3=LFTable(PhiStarCol,12)*1e-4;
	
	string MstarUnits="M-5log10h70";
	string phistarUnits="(Mpc/h70)^-3";
	
	cout <<" ... Schechter parameters"<<endl; 
	cout <<"     z range     Mstar     alpha     phistar     spec type"<<endl;
	cout <<"    0.75<z<1.0  "<<MstarAz3<<"     "<<alpAz3<<"       "<<phistarAz3<<"        All"<<endl;
	cout <<"                "<<MstarEz3<<"     "<<alpEz3<<"       "<<phistarEz3<<"         Early"<<endl;
	cout <<"                "<<MstarLz3<<"     "<<alpLz3<<"       "<<phistarLz3<<"          Late"<<endl;
	cout <<"                "<<MstarSz3<<"     "<<alpSz3<<"        "<<phistarSz3<<"        Starburst"<<endl<<endl;
	// ----------------------- end subsitute in LFParameters -------------------
	
	
	cout<<"    Mass to Galaxy number conversion"<<endl;
	Schechter schAz3(phistarAz3,MstarAz3,alpAz3);
	Schechter schEz3(phistarEz3,MstarEz3,alpEz3);
	Schechter schLz3(phistarLz3,MstarLz3,alpLz3);
	Schechter schSz3(phistarSz3,MstarSz3,alpSz3);
	double schmin=-24, schmax=-13;// units of "M-5log10h70"
	int schnpt=10000;
	cout <<"     ... integrating from Mbright="<< schmin <<" "<< MstarUnits;
	cout <<" to Mfaint="<< schmax <<" "<< MstarUnits <<", with step="<< schnpt <<endl;
	double nz3=schAz3.Integrate();//(schmin,schmax,schnpt);
	cout <<"     ... number density of galaxies: "<< nz3 <<" Mpc^-3"<<endl;
	
	double pixVol = m2g.ReturnPixVol();
	cout <<"     pixel volume="<< pixVol <<" Mpc^3"<<endl;
	float conv = pixVol*nz3;
	cout <<"     actual gals per pixel vol="<< pixVol*nz3 << endl;
	cout <<"     gals per pixel vol="<< conv << endl;
	//m2g.ConvertToMeanNGal(conv); //just multiplies mass_ by conv
	//m2g.setSimProperties(SkyArea, doAbsMagCut, isConstDens)
	
	// THIS WHOLE THING NEEDS TO BE RE WRITTEN
	/*
	// SET RANDOM GRID
	m2g.SetRandomGrid(meandens);
	  
	  
	// COMPUTE SIMLSS POWER SPECTRA 
	int nbin = 175;
	cout << "    Computing correction power spectra from SimLSS"<<endl<<endl;
	Mass2Gal m2g2(fin,su,rg);
	TArray<r_8> densf;
	m2g2.ODensArray(densf);
	
	double Dx=m2g2.ReturnDX(); 
	double Dy=m2g2.ReturnDY(); 
	double Dz=m2g2.ReturnDZ();
	
	PowerSpec psim(dens,Dx,Dy,Dz);
	PowerSpec psimf(densf,Dx,Dy,Dz);
	
	double kmin=0; 
	double kmax = sqrt(pow(psim.ReturnKxMax(),2)+pow(psim.ReturnKyMax(),2)+pow(psim.ReturnKzMax(),2));		
	HProf hp(kmin, kmax, nbin);
	HProf hpf(kmin, kmax, nbin);
	r_4 volsim=m2g2.ReturnCubeVol();
	double sfc=psim.AccumulatePowerSpectra(hp);
	double sfcf=psimf.AccumulatePowerSpectra(hpf);
	
	// START REALISATIONS 
	//for (int ir=0; ir<nreal; ir++)
	//	{
			
		// APPLY SELECTION FUNCTION 
		cout <<"    Applying the selection function"<<endl;
		m2g.ApplySF(*sfp);
	
		// APPLY PHOTO-Z SMEARING 
		cout <<"    Applying the photo-z convolution function"<<endl;	
		m2g.ApplyPZConv(pzconv);
	
		if (DEBUG)
			{
			TArray<int_8> ngal;
			m2g.NGalArray(ngal);
			string debugfile; 
			debugfile= debug +"_ngal.fits";
			cout <<"    Writing ngal array to "<<debugfile<<endl;
			FitsInOutFile fos(debugfile, FitsInOutFile::Fits_Create);
			fos << ngal;	
			}
	
		// RETURN ARRAYS
		cout <<"    Return the arrays"<<endl;
		TArray<r_8> ngalsm,rgalsm;
		m2g.NGalSmArray(ngalsm);	
		m2g.RGalSmArray(rgalsm);

		if (DEBUG)
			{
			string debugfile; 
			debugfile= debug +"_ngalsm.fits";
			cout <<"    Writing ngalsm array to "<<debugfile<<endl;
			FitsInOutFile fos(debugfile, FitsInOutFile::Fits_Create);
			fos << ngalsm;	
			}
	
		// NORMALISE ARRAYS 
		cout <<"    Normalise the arrays"<<endl;
		ngalsm *= ( ngalsm.Size()/ngalsm.Sum());
		rgalsm *= ( rgalsm.Size()/rgalsm.Sum());
	
		double meang, sigg, meangr, siggr;
		//MeanSigma(ngalsm, meang, sigg);
		//MeanSigma(rgalsm, meangr, siggr);
		
		// FT 
		cout <<"    FT the galaxy and random grids"<<endl;
		PowerSpec psgals(ngalsm,Res); // does FT in constructor
		PowerSpec psrand(rgalsm,Res);
	
		// COMPUTE POWER SPECTRUM 
		cout <<"    Compute power spectra"<<endl;
		cout <<"    kmin = "<<kmin<<", kmax (ngals) = "<<kmax<<", nbin="<<nbin<<endl;
		
		HProf hpgals(kmin, kmax, nbin);
		HProf hprand(kmin, kmax, nbin);
		r_4 volgrid=m2g.ReturnCubeVol(); cout <<"    volume = "<<volgrid<<endl;
		double sumFC=psgals.AccumulatePowerSpectra(hpgals);
		double sumFCr=psrand.AccumulatePowerSpectra(hprand);
	//cout <<"    Check: sum of Fourier coefficients (ngals) = "<<sumFC<<endl;
//	cout <<"           variance of real space field / 2 (ngals) = "<< sigg*sigg/2 <<endl;
//	cout <<"    Check: sum of Fourier coefficients (rgals) = "<<sumFCr<<endl;
//	cout <<"           variance of real space field / 2 (rgals) = "<< siggr*siggr/2 <<endl;
	
		// OUTPUT PS's INTO A TEXT FILE 
		//std::stringstream ss;
		//ss << ir;
		string outfile;
		outfile = outfileroot + "_ngal.txt";//"_real" + ss.str() + "_ngal.txt";
		cout <<"    Write ngal power spectrum to "<<outfile<<endl;
		psgals.WritePS(outfile,hpgals,volgrid,hp,hpf,volsim);
		outfile = outfileroot + "_rgal.txt"; //"_real" + ss.str() + "_rgal.txt";
		cout <<"    Write rgal power spectrum to "<<outfile<<endl;
		psrand.WritePS(outfile,hprand,volgrid,hp,hpf,volsim);
		
		m2g.ResetSmGrids();
	//	}*/
			
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " sim_mcgrids.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " sim_mcgrids.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " sim_mcgrids.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of sim_mcgrids.cc program  Rc= " << rc << endl;
  return rc;	
}
