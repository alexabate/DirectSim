/**
  * @file  computepsfromarray.cc
  * @brief compute power spectrum from gridded data
  *
  * @todo put cosmology in sub-grid file header
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

// sophya
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
#include "swfitsdtable.h"
#include "resusage.h"

// DirectSim
#include "mydefrg.h"
#include "geneutils.h"
#include "cat2grid.h"
#include "powerspec.h"
#include "mass2gal.h"
#include "pkspectrum.h"
#include "fitkbaoscale.h"
#include "chisqstats.h"


void usage(void);
void usage(void) {

	cout << endl<<" Usage: computepsfromarray [...options...]" << endl<<endl;
	
	cout << "  Compute power spectrum from gridded galaxy data. The     "<<endl;
	cout << "  output power spectrum is correctly normalized and the    "<<endl;
	cout << "  distortion in the simulated density distribution (from   "<<endl;
	cout << "  setting over-densities with delta<-1 equal to -1) is     "<<endl;
	cout << "  properly taken account of. This had to be done because   "<<endl;
	cout << "  delta<-1 corresponds to a negative (unphysical) density. "<<endl;
	cout << "  This can be interpreted to arise from structure formation"<<endl;   
	cout << "  on nonlinear scales not included in the simulation method."<<endl;
	cout << endl;
	
	cout << "  The file containing the gridded data for power spectrum  "<<endl;
	cout << "  analysis is supplied with option -C. This file is        "<<endl;
	cout << "  probably output from the subfromfull program.            "<<endl;
	cout << endl;
	
	cout << "  In order to do the correction described above, either the"<<endl;
	cout << "  power spectra of the undistorted and distorted over-     "<<endl;
	cout << "  density distribution must be supplied to the program, or "<<endl;
	cout << "  the original over-density distribution itself (so those  "<<endl;
	cout << "  power spectra can be computed here). The file containing "<<endl;
	cout << "  the density distribution or its power spectra are        "<<endl;
	cout << "  supplied with option -S. The density distribution file   "<<endl;
	cout << "  is probably output from the simdensity program.          "<<endl;
	cout << endl;
	
	cout << "  The shot noise power spectrum is also computing using    "<<endl;
	cout << "  gridded data made from a random catalog read in from the "<<endl;
	cout << "  same file as the gridded galaxy data. "<<endl;
	cout << endl;
	
	cout << "  The mean density of the over-density distribution is     "<<endl;
	cout << "  needed to properly normalized the power spectrum. It is  "<<endl;
	cout << "  either read from the file header, or using option -a it  "<<endl;
	cout << "  is passed to the program as an argument (overriding any  "<<endl;
	cout << "  value in the file header).                               "<<endl;
	cout << endl;
	
	cout << "  If the galaxies have approximately Gaussian photo-z errors"<<endl;
	cout << "  the magnitude of this error sigma_z (in format sigma_z*1+z)"<<endl;
	cout << "  should be supplied to the program with the -e option. Then"<<endl;
	cout << "  the power spectrum can be undampled accordingly. To turn  "<<endl;
	cout << "  off the undamping even if the photo-z error is non-zero  "<<endl;
	cout << "  use option -d.                                           "<<endl;    
	cout << endl;
	
	cout << "  The maximum k to use in the power spectrum analysis is set"<<endl;
	cout << "  with option -m.                                          "<<endl;
	cout << endl;
	
				
	cout << "  This code uses the cosmology of double h=0.71, OmegaM=0.267804,"<<endl;
	cout << "  OmegaL=0.73 (SimLSS cosmology)"<<endl;
	cout <<endl;
	
	cout << "  EXAMPLE: "<<endl;
	cout << endl;
	
	cout << "  $ computepsfromarray -C subgrids.fits -S overdensity.fits "<<endl;
	cout << "                       -O powerspectra -d -m 0.5            "<<endl;
	cout << endl;
	
	cout << " -C : infile : file containing gridded data                        "<<endl;
	cout << " -S : overdensityfile : file containing over-density distribution  "<<endl;
	cout << "                        or over-density power spectra              "<<endl;
	cout << " -O : outfile : root filename of text file the galaxy power spectra"<<endl;
	cout << "                are written to                                     "<<endl;
	cout << " -a : meandens : specify mean density of overdensity distribution  "<<endl;
	cout << " -c : compute power spectrum of simlss from same sub-grid as galaxy"<<endl;
	cout << "      data                                                         "<<endl;
	cout << " -d : Don't undamp photo-z error damping of Fourier coefficients   "<<endl;
	cout << " -e : photoZerror : size of photometric redshift error             "<<endl;
	cout << " -m : maxk_in_calc : maximum kradial used in power spectrum comp   "<<endl;
	cout << " -x : doPixCorr : turn off pixel shape correction                  "<<endl;
	cout << " -o : OutputRoot : root stem of output filename objects are written"<<endl;
	cout << "                   to if want to debug                             "<<endl;
	cout << endl;
	}
	
int main(int narg, char* arg[]) {
	cout << " ==== computepsfromarray_main.cc program , compute power spectrum";
	cout << " from grided data  ==== " <<endl;
	
	
	// Make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	

	// Set defaults etc ....
	// FILES TO READ IN/OUT
	string infile, overdensityfile, outfileroot, subinfo;
	// HOW TO COMPUTE SIMLSS POWER SPECTRUM
	bool computeOvDensityPS = false; // default don't compute, read from file
	// CATALOG-TYPE AND REDSHIFT-TYPE PARAMETERS
	double photoZerror = 0;		// Photo-z error size is 0
	// IF HAVE PHOTO-Z ERR>0 PARAMETERS
	//double coeff = 1;		// keep all kradial <= coeff / sig_r
	bool doUnDamp = true;	// undamp Fourier components
	// POWER SPECTRUM COMPUTATION PARAMETERS
	int nbin=175;			// Number of k bins in power spectrum
	bool doPixCorr = true;	// Correct for pixel size smoothing
	double maxk_in_calc = 1000; // Set maximum radial k in ps calc
	bool setMaxK=false;		// Set max k separately from z error
	bool isSameSub=false;
	double meandens;		// mean density of simlss grid after delta<-1 =-1
	bool isMeanDensitySpec = false;
	// DEBUGGING
	string debug_out = "TEST";
	bool DoDebug = false;
	
	
	// decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hdxcC:S:O:e:o:m:a:")) != -1) {
		switch (c) {
	  	    case 'C' :
			    infile = optarg;
		        break;
	        case 'S' :
		        overdensityfile = optarg;
		        break;
	        case 'O' :
		        outfileroot	= optarg;
		        break;
	        case 'd' :
		        doUnDamp = false; // don't undamp Fourier components
		        break;
	        case 'a' :
		        sscanf(optarg,"%lf",&meandens);
		        isMeanDensitySpec=true;
		        break;
	        case 'e' :
		        sscanf(optarg,"%lf",&photoZerror);
		        break;
	        case 'x' :
		        doPixCorr = false; 
		        break;
	        //case 'k' :
		    //    sscanf(optarg,"%lf",&coeff);
		    //    break;
	        case 'm' :
		        setMaxK=true;
		        sscanf(optarg,"%lf",&maxk_in_calc);
		        break;
	        case 'o' :
		        debug_out = optarg;
		        DoDebug = true;
		        break;
	        case 'c' :
		        isSameSub=true;
		        break;
	        case 'h' :
	            default :
		        usage(); return -1;
	        }
	    }


    // process whether reading in or computing over-density power spectrum
	int nc = 4; // 4 characters 
	int lf = overdensityfile.size(); // length of the filename
	string endf; // last 4 characters of overdensityfile
	for (int i=0;i<nc; i++) {
		int val=lf-nc+i;
		endf.push_back(overdensityfile[val]);
		}
	string textfile=".txt";
	string fitsfile="fits";

	if ( strcmp(textfile.c_str(),endf.c_str()) == 0 ) { // if strings ARE the same
	    cout <<" computeOvDensityPS = false"<<endl;
		computeOvDensityPS = false; // will be reading it from a file
		}
	else if ( strcmp(fitsfile.c_str(),endf.c_str()) ==0 ) {// if strings ARE the same
	    cout <<" computeOvDensityPS = true"<<endl;
		computeOvDensityPS = true;
		}
	else 
		throw ParmError("SimLSS file is of unknown file type: not fits or .txt");


	// Command line argument printing
	cout << "     Printing command line arguments ... "<<endl;
	cout << "     Galaxy sub-array file is "<< infile <<endl;

	cout << "     Maximum k_radial in analysis given by "<< maxk_in_calc <<endl;
	if (photoZerror>0) {
		if(doUnDamp)
			cout << "     Undamping Fourier coefficients"<<endl;
		else
			cout << "     NOT undamping Fourier coefficients"<<endl;
		}
	if (photoZerror)
		cout << "     Photo-z error = "<< photoZerror <<endl;
	if(doPixCorr)
		cout << "     Pixel correction is ON"<<endl;
	else
		cout << "     Pixel correction is OFF"<<endl;
	if (computeOvDensityPS) {
		cout << "     Reading simulated over-density grid from file "<< overdensityfile <<endl;
		if(isSameSub) {
			cout << "     Computing over-density power spectrum using same ";
			cout << " sub-array as galaxy data"<<endl;
			}
		}
	else
		cout << "     Reading over-density power spectrum from file "<< overdensityfile <<endl;
	cout << "     Galaxy power spectrum will be output to "<< outfileroot <<endl;
	if (DoDebug)
		cout << "     Output root filename for debugging is "<< debug_out <<endl;
	cout << endl;

  
	int rc = 1;  
	try {  // exception handling try bloc at top level
	
	// monitor memory usage
	ResourceUsage res;
	
    // Read in gridded galaxy data
	cout <<"     Read in file "<< infile <<endl;
	FitsInOutFile fin(infile, FitsInOutFile::Fits_RO);   
	TArray<r_8> ngals,wngals,wrgals;
	fin >> ngals;
	fin >> wrgals;
	fin >> wngals; // in the case where there is no selection effects on the
	               // galaxy catalog: ngals=wngals


    // Read data from file header
	sa_size_t nx = wngals.SizeX(); 
	sa_size_t ny = wngals.SizeY(); 
	sa_size_t nz = wngals.SizeZ(); 
	double z_center = atof(fin.KeyValue("ZCEN").c_str()); 
	double grid_res = atof(fin.KeyValue("R").c_str()); 
	if(!isMeanDensitySpec)
		meandens=atof(fin.KeyValue("MeanOverDensity").c_str());
	cout << "     Size of sub-array Nx,Ny,Nz = "<< nx <<","<< ny <<","<< nz;
	cout <<", resolution = "<< grid_res;
	cout <<" Mpc, mean density of distorted overdensity grid = "<< meandens <<endl;

	
	// Set cosmology (should be reading this from the header)
	cout << "     Initialise cosmology:"<<endl;
	double h = 0.71, OmegaM = 0.267804, OmegaL = 0.73;
	SimpleUniverse su(h, OmegaM, OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	double OmegaB = su.OmegaBaryon();
	cout <<"     OmegaK="<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
	cout <<", OmegaL="<< OmegaL <<", OmegaB="<< OmegaB <<", H0="<< su.H0() <<endl;
	cout << endl;
	RandomGenerator rg; // need this for cat2grid

    
    // Calculate photo-z error if needed
	double photoZerrorMpc = 0;
	if (photoZerror>0) {
		cout << "     Calculate comoving photo-z error:"<<endl;
		su.SetEmissionRedShift(z_center);
		photoZerrorMpc = su.ZErr2CoDistErr(photoZerror);
		cout <<"    Redshift of array center = "<< z_center <<", photo-z error = ";
		cout << photoZerror <<", photo-z co-distance error = "<< photoZerrorMpc;
		cout <<" Mpc"<<endl;
		cout <<endl;
		}
	
	
	// Mean and sigma of gridded galaxy data (mean should be ~0)
	double meang, sigg, meangw, siggw, meangr, siggr;
	MeanSigma(ngals, meang, sigg);
	cout << "     Mean and Sigma of (raw) galaxy fluctuation field ..."<<endl;
	cout << "     Mean="<< meang <<", Sigma="<< sigg <<endl<<endl;
	MeanSigma(wrgals, meangr, siggr);
	cout << "    ... of random galaxy grid: Mean="<< meangr <<", Sigma="<< siggr <<endl;
	MeanSigma(wngals, meangw, siggw);
	cout << "    ... of weighted galaxy grid: Mean="<< meangw <<", Sigma="<< siggw <<endl;
	cout << "    (above will be same as raw galaxy fluctuation field if "<<endl;
	cout << "     original catalog had no selection effects)"<<endl;
	cout << endl;
	
	
	r_4 volcat = wngals.SizeX()*wngals.SizeY()*wngals.SizeZ()*pow(grid_res,3); 
	cout <<"    Grid volume = "<< volcat <<endl<<endl;
	
	
	double sum_FourierCoeffs; // sum of Fourier coefficients (for checking)
	double kmax = PI/grid_res;
	double kmin = 0.; 
	
	
	// If over-density grid read in - compute its power spectrum
	HProf hp(kmin, kmax, nbin);
	HProf hpf(kmin, kmax, nbin);
	r_4 volsim;
	if (computeOvDensityPS) {
	
	    cout << "     Compute over-density power spectra "<<endl;
	    cout << "     i) with unmodified density distribution "<<endl;
	    cout << "    ii) when grid cells of delta<-1 were set to delta=-1"<<endl;
	
	
	    cout << "     Reading over-density file " << overdensityfile << endl<<endl;  
	    FitsInOutFile fsin(overdensityfile, FitsInOutFile::Fits_RO);
	    
	    // This is the unfudged over density distribution
	    TArray<r_8> denstmp; 
	    fsin >> denstmp; // delta distribution
	    cout <<"     Mean and Sigma of density fluctuations BEFORE grid cells ";
	    cout <<" of delta<-1 were set to delta=-1 ..."<<endl;
	    double meanm, sigm;
	    MeanSigma(denstmp, meanm, sigm);
	    cout <<"     Mean="<< meanm <<", Sigma="<< sigm <<", Variance="<< sigm*sigm <<endl;
	    
	    TArray<r_8> dens;
	    if(isSameSub) {
	    
		    cout <<"     Select same sub-array as galaxy catalog, ";
		    cout <<" reading pixel values from sub-array info file "<< subinfo <<endl;
		    ifstream ifs(subinfo.c_str());
		    Array B;
		    sa_size_t nr, nc;
		    B.ReadASCII(ifs,nr,nc);
		    sa_size_t x1 = B(0,1), x2 = B(2,1);
		    sa_size_t y1 = B(0,2), y2 = B(2,2);
		    sa_size_t z1 = B(0,3), z2 = B(2,3);
		    cout <<"     Read out distorted density distribution ... "<<endl;
		    cout << endl;
		    dens = denstmp(Range(z1,z2),Range(y1,y2),Range(x1,x2)).PackElements();
		    volsim = volcat;
		    }
	    else { 
	        cout <<"     Read out distorted density distribution ... "<<endl;
	        dens = denstmp;
	        }
	    denstmp.ZeroSize();
	
	
        // This is fudged over-density distribution
	    Mass2Gal m2g(fsin, su, rg);
	    cout << endl;
	    double zc = m2g.ReturnZref();
	    cout << endl;
	    
	    TArray<r_8> densf;
	    double volsim;
	    if(isSameSub) {
		
		    cout <<"     Select same sub-array as galaxy catalog, ";
		    cout <<" reading pixel values from sub-array info file "<< subinfo <<endl;
		    ifstream ifs(subinfo.c_str());
		    Array B;
		    sa_size_t nr, nc;
		    B.ReadASCII(ifs,nr,nc);
		    sa_size_t x1 = B(0,1), x2 = B(2,1);
		    sa_size_t y1 = B(0,2), y2 = B(2,2);
		    sa_size_t z1 = B(0,3), z2 = B(2,3);
		    cout <<"     Read out density distribution ... "<<endl;
		    cout << endl;
		    densf = m2g.ExtractSubArray(x1,x2,y1,y2,z1,z2);
		    }
	    else {
	        cout <<"     Read out whole array of density distribution ... "<<endl;
		    m2g.ODensArray(densf);
		    volsim = m2g.ReturnCubeVol();
		    }
	    cout << endl;
	
	    cout <<"     Mean and Sigma of density fluctuations AFTER grid cells ";
	    cout <<" of delta<-1 were set to delta=-1 ..."<<endl;
	    double meanmf, sigmf;
	    m2g.returnDensMeanSigma(meanmf, sigmf);
	    cout <<"     Mean="<< meanmf <<", Sigma="<< sigmf <<", Variance="<< sigmf*sigmf <<endl;
	    
	    
	    if ( (meanmf-1)!=meandens ) {
	        cout <<"     mean of fudged density cube (="<< meanmf-1 <<")";
	        cout <<" not equal to stored mean of fudge density cube";
	        cout <<" (="<< meandens <<")"<<endl;
		    cout <<"     probably because mean was never stored in the first";
		    cout <<" place (likely if 2nd value above is 0)"<<endl;
		    if (meandens==0)
			     meandens=meanmf-1;
		    cout <<"     Mean density using to correct power spectrum = ";
		    cout << meandens <<endl;
		    }
	    cout <<endl;


	    res.Update();
	    cout << " Memory size (KB):" << res.getMemorySize() << endl;
	    cout << " Resource usage info : \n" << res << endl;
	
	
	    cout <<"     Zero size array"<<endl;
	    m2g.ZeroSizeMassArray();
	 
	 
	    cout <<"     Computing overdensity power spectra"<<endl<<endl;
	    double Dx=m2g.ReturnDX(); double Dy=m2g.ReturnDY(); double Dz=m2g.ReturnDZ();
	    if(Dx-Dy>0||Dz-Dy>0) cout <<"   CAREFUL! Pixels are not cuboid"<<endl;
	    PowerSpec powerSpectrum_overdensity(dens,Dx,Dy,Dz);
	    PowerSpec powerSpectrum_overdensityf(densf,Dx,Dy,Dz);
	    cout <<"     Power spectrum defined using:"<<endl;
	    cout <<"     kmin="<< kmin <<", kmax="<< kmax <<", nbin="<< nbin <<endl;
    
	
	    sum_FourierCoeffs = powerSpectrum_overdensity.AccumulatePowerSpectra(hp,doPixCorr);
	    powerSpectrum_overdensity.ZeroSizeArrays();
	    cout <<"     Check: sum of Fourier coefficients = "<< sum_FourierCoeffs <<endl;
	    cout <<"            variance of real space field / 2 = "<< sigm*sigm/2 <<endl;
	    sum_FourierCoeffs = powerSpectrum_overdensityf.AccumulatePowerSpectra(hpf,doPixCorr);
	    powerSpectrum_overdensityf.ZeroSizeArrays();
	    cout <<"     Check: sum of Fourier coefficients = "<< sum_FourierCoeffs <<endl;
	    cout <<"            variance of real space field / 2 = "<< sigmf*sigmf/2 <<endl;
	    }


   cout <<"     Overdensity cube volume="<< volsim <<" Mpc";
   cout <<", survey volume="<< volcat <<" Mpc"<<endl<<endl;

	
	// Compute power spectrum
	cout <<"     Compute power spectrum of gridded galaxy data "<<endl;
	
	
	PowerSpec powerSpectrum_weighted(wngals,grid_res); // does FT in constructor
	powerSpectrum_weighted.Setzc(z_center);
	
	HProf histogram_weighted(kmin, kmax, nbin);
	sum_FourierCoeffs = powerSpectrum_weighted.AccumulatePowerSpectra
	    (histogram_weighted, doPixCorr, maxk_in_calc, photoZerrorMpc, doUnDamp);
	
	cout <<"     Check: sum of Fourier coefficients = "<< sum_FourierCoeffs <<endl;
	cout <<"            variance of real space field / 2 = "<< siggw*siggw/2 <<endl;
	
	/*if (photoZerrorMpc>0) {
	    cout <<"     If throwing out some modes the above will NOT be equal"<<endl;
		cout <<"     max k_radial = "<< maxk_in_calc <<", kmax = "<< kmax <<endl;
		if ( (coeff/photoZerrorMpc) < kmax)
		   cout <<"    Will throw out k_radial > "<< coeff/photoZerrorMpc <<endl;
		}*/
		
		
	string outfile;
	
	// Write out power spectrum
	outfile = outfileroot + "_wngal.txt";
	if(computeOvDensityPS)
	    powerSpectrum_weighted.WritePS(outfile,histogram_weighted,volcat,hp,hpf,volsim,meandens);
	else
	    powerSpectrum_weighted.WritePS(outfile,histogram_weighted,volcat,overdensityfile,meandens);
	cout << endl;
	
	
	// Compute shot noise power spectrum
	cout <<"     Compute shot noise power spectrum from random catalog grid"<<endl;
	
	PowerSpec powerSpectrum_random(wrgals, grid_res);
	powerSpectrum_random.Setzc(z_center);
		
	HProf histogram_random(kmin, kmax, nbin);
	sum_FourierCoeffs = powerSpectrum_random.AccumulatePowerSpectra
	      (histogram_random, doPixCorr, maxk_in_calc, photoZerrorMpc, doUnDamp);	
	
	cout <<"     Check: sum of Fourier coefficients = "<< sum_FourierCoeffs <<endl;
	cout <<"            variance of real space field / 2 = "<< siggr*siggr/2 <<endl;
	outfile = outfileroot + "_wrngal.txt";
	
	
	// Write out shot noise power spectrum
	if(computeOvDensityPS)
	    powerSpectrum_random.WritePS(outfile,histogram_random,volcat,hp,hpf,volsim,meandens);
	else
	    powerSpectrum_random.WritePS(outfile,histogram_random,volcat,overdensityfile,meandens);
	cout <<endl;
	

		
	  }  // End of try bloc 
  
  
 catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " computepsfromarray.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
	}
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " computepsfromarray.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
	}
  catch (...) {  // catching other exceptions
    cerr << " computepsfromarray.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
	}
  cout << " ==== End of computepsfromarray.cc program  Rc= " << rc << endl;
  return rc;	
}
