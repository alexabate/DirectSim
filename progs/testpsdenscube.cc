#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <typeinfo>
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
//#include "cat2grid.h"
#include "powerspec.h"
#include "mass2gal.h"

#include "pkspectrum.h"
//#include "fitkbaoscale.h"
//#include "chisqstats.h"

#define PI 3.141592
/*
   Program to compute the power spectrum directly from the 
   over-density cube with and without fudging pixels with
   over-density<-1 (ie mass<0)
   
   This code uses the cosmology of double h=0.71, OmegaM=0.267804, OmegaL=0.73 
   (default over-density simulation cosmology)
*/

void usage(void);
void usage(void)
{
	cout << endl<<" Usage: testpsdenscube [...options...]" << endl<<endl;
	cout << " -S : cubefile : FITS filename containing over-density cube"<<endl;
	cout << " -O : outfile : filename of text file the power spectra are written to";
	cout <<endl;
	cout << " -r : Res : resolution of grid [default=6]"<<endl;
	cout << endl;
}
int main(int narg, char* arg[])
	{
	cout << " ==== testpsdenscube.cc program , compute power spectrum from";
	cout << " over-density cube  ==== " <<endl;
	// Make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	
	cout << " ==== setting defaults ===="<<endl;
	// Set defaults etc ....
	// FILES TO READ IN/OUT
	string cubefile,outfile;
	int xplanes=1;
	double res = 6;
	// POWER SPECTRUM COMPUTATION PARAMETERS
	int nbin=175;			// Number of k bins in power spectrum
	bool pixcorr = true;		// Correct for pixel size smoothing

	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hS:O:r:")) != -1) 
	{
	switch (c) 
	{
	  case 'S' :
		cubefile = optarg;
		break;
	  case 'O' :
		outfile = optarg;
		break;
	  case 'r' :
		sscanf(optarg,"%lf",&res);
		break;
	  case 'h' :
	  default :
		usage(); return -1;
	}
	}

		
	//******************** Print all command line options********************//
	cout << "     Reading command line arguments ... "<<endl;
	// IN FILE TYPE
	cout << "     SimLSS file is "<<cubefile<<endl;
	cout << "     Resolution of grid is "<<res<<" Mpc"<<endl;
	// OUTPUT FILES
	cout << "     Power spectra will be output to "<<outfile<<endl;
	cout << endl;
	//********************* end command line arguments *********************//
  
	int rc = 1;  
	try {  // exception handling try bloc at top level

	// k range
	double kmin=0, kmax=PI/res;

	/* set cosmo */
	cout << "0.1/ Initialise cosmology: (same as over-density simulation cosmology)";
	cout <<endl;
	double h=0.71, OmegaM=0.267804, OmegaL=0.73;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout <<"    OmegaK="<<su.OmegaCurv()<<", OmegaM="<<su.OmegaMatter();
	cout <<", OmegaL="<<su.OmegaLambda();
	cout <<", OmegaB="<<su.OmegaBaryon()<<", H0="<<su.H0()<<endl;
	cout << endl;

	/* read over-density cube */
	cout << "0/ Reading input file = " << cubefile << endl;  
	FitsInOutFile fin(cubefile,FitsInOutFile::Fits_RO);
	TArray<r_8> drho;
	fin >> drho;
	cout << drho.Info();
	cout << "    print original drho array size: "<<drho.SizeX()<<"x"<<drho.SizeY();
	cout <<"x"<<drho.SizeZ()<<endl<<endl<<endl;
	
	// the array to fill
	TArray<r_8> dens;
	   
	/* remove N extra planes */
	cout << "1/ Initialise Mass2Gal: remove planes" << endl;
	RandomGenerator rg;
	Mass2Gal m2g(drho,su,rg,xplanes);

	/* read in cube and pixel properties from fits header */
	cout << "1.1/ Read in cube properties from fits header" << endl;
	m2g.ReadHeader(fin);
	// xplanes should be the difference between NZ and drho.SizeX()
	int NZ=m2g.ReturnNZ();
	int diff = drho.SizeX()-NZ;
	drho.ZeroSize();
	if( xplanes!=abs(diff) )
		throw ParmError("ERROR: removed wrong number of planes from over-density cube");
		
	/* Check the std and mean */
	m2g.MassArray(dens);
	double meanc,sigc,meanfc,sigfc;
	MeanSigma(dens, meanc, sigc);
	cout << endl<<"1.2/ RAW DENS CUBE STATS: Mean = " << meanc << ", Var = ";
	cout << sigc*sigc<< ", Var/2 = " << sigc*sigc/2 << endl;
	cout<<endl;
	double volsim = m2g.ReturnCubeVol();
	cout <<"    Simulation cube volume="<<volsim<<endl;
	/* POWER SPECTRUM OF NON-FUDGED OVER-DENSITY CUBE*/
	
	cout << "2a/ Compute power spectrum of non-fudged over-density cube" << endl;

	/* ps of regular cube */
	PowerSpec psim(dens,res);
	HProf hp(kmin,kmax,nbin);
	psim.AccumulatePowerSpectra(hp,pixcorr);
	psim.ZeroSizeArrays();
	dens.ZeroSize();
	cout << endl;

	cout << "2b/ Compute power spectrum of fudged over-density cube" << endl;
	/* fudge */
	m2g.CleanNegativeMassCells(); 
	m2g.MassArray(dens);
	MeanSigma(dens, meanfc, sigfc);
	cout << endl<<"2b.2/ FUDGED DENS CUBE STATS: Mean = " << meanfc << ", Var = ";
	cout << sigfc*sigfc<< ", Var/2 = " << sigfc*sigfc/2 << endl;

	/* ps of fudge cube */
	PowerSpec pf(dens,res);
	HProf hf(kmin,kmax,nbin);
	pf.AccumulatePowerSpectra(hf,pixcorr);
	
	/* write to a file*/
	cout << "3/ Write power spectra to file " << endl;
	psim.Write2PS(outfile,hp,hf,volsim,meanc,meanfc);
	cout << endl;

	  }  // End of try bloc 
  
  
 catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testpsdenscube.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
	}
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testpsdenscube.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
	}
  catch (...) {  // catching other exceptions
    cerr << " testpsdenscube.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
	}
  cout << " ==== End of testpsdenscube.cc program  Rc= " << rc << endl;
  return rc;	
}
