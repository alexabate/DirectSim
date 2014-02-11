/**
  * @file  fitkbao.cc
  * @brief Given an input power spectrum + errors fit the BAO scale to "wiggles only"
  *        power spectrum
  *
  * @todo read cosmology from file header
  *
  * @author Alex Abate
  * Contact: abate@email.arizona.edu
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <functional>
#include <numeric>
#include <algorithm>

// sophya
#include "sopnamsp.h"
#include "histinit.h"
#include "hisprof.h"
#include "histerr.h"
#include "histos.h"
#include "datatable.h"
#include "fitshdtable.h"
#include "swfitsdtable.h"
#include "fitsarrhand.h"
#include "fiosinit.h"
#include "tarray.h"

// DirectSim
#include "geneutils.h"
#include "cosmocalcs.h"
#include "pkspectrum.h"
#include "fitkbaoscale.h"



void usage(void);
void usage(void) {

	cout << endl<<" Usage: fitkbao [...options...]              " << endl<<endl;
	
	cout << "  Given an input power spectrum + errors fit the BAO scale."<<endl;
	cout << endl;
 
	cout << "  The input power spectrum has already been corrected for  "<<endl;
	cout << "  shot noise, selection, photo-z etc and is supplied to the"<<endl;
	cout << "  program with option -P                                   "<<endl;
	cout << endl;

	cout << "  Method: divide observed power spectrum by a reference    "<<endl;
	cout << "  power spectrum and fit a sine wave described by some     "<<endl;
	cout << "  amplitude, and a characteristic scale. To compute the    "<<endl;
	cout << "  reference power spectrum the redshift of the observed "<<endl;
	cout << "  power spectrum must be supplied with option -z, and the  "<<endl;
	cout << "  values of sigma_8 and the spectral index parameters must "<<endl;
	cout << "  supplied with option -c "<<endl;
	cout << endl;
	
	cout << "  Results are written to files starting with the root name "<<endl;
	cout << "  supplied with the -O option. The chi-square values,      "<<endl;
	cout << "  reference power spectrum, best-fit sinusoid and fit      "<<endl;
	cout << "  results are written to files "<<endl;
	cout << endl;
	
	cout << " -P : PSFile: power spectrum file to read in               "<<endl;
	cout << "              (3 columns: k (Mpc^-1), P(k) Mpc^3, err)     "<<endl;
	cout << " -O : outfile_root: file root name to write results to     "<<endl; 
	cout << " -z : zref: redshift of power spectrum                     "<<endl;
	cout << " -c : sigma8,n: sigma8 and spectral index                  "<<endl;
	cout << endl;
}



int main(int narg, char *arg[]) {

	SophyaInit();
	FitsIOServerInit();
  
	// input power spectrum
	string ps_file;
	// output file
	string outfile_root;
	double zref;
	// cosmology
	double n=1, sig8 = 0.8;
	double maxk = 1;

	//--- decoding command line arguments 
	char c;
	while ((c = getopt(narg,arg,"hP:O:z:d:c:")) != -1) {
	    switch (c) {
		    case 'P' :
			    ps_file = optarg;
			    break;
		    case 'O' :
			    outfile_root = optarg;
			    break;
		    case 'z' :
			    sscanf(optarg,"%lf",&zref);
			    break;
		    case 'c' :
			    sscanf(optarg,"%lf,%lf",&sig8,&n);
			    break;
		    case 'h' :
		        default :
			    usage(); return -1;
		    }
	    }
	
	
	cout << "     Printing command line arguments ... "<<endl<<endl;
	cout << "     Reading in observed power spectrum from :"<< ps_file <<endl;
	cout << "     Saving results to files beginning "<< outfile_root <<endl;
	cout <<endl;
	
    try {
	
	
	// Read in power spectrum file
	cout << "     Read in power spectrum file "<< ps_file <<endl;
	ifstream ifs(ps_file.c_str());
	TArray<r_8> power_spectrum;
	sa_size_t nr, nc;
	power_spectrum.ReadASCII(ifs,nr,nc);
	cout << power_spectrum ;


	// Initalise SimLSS cosmology
	cout << "     Initialise cosmology:"<<endl;
	double h = 0.71, OmegaM = 0.267804, OmegaL = 0.73, R = 8;
	SimpleUniverse su(h, OmegaM, OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	double OmegaB = su.OmegaBaryon();
	cout << "     OmegaK="<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
	cout << ", OmegaL="<< OmegaL <<", OmegaB="<< OmegaB <<", H0="<< su.H0() <<endl;
	cout << endl;
	

	// Initialise FitBAOScale
	cout << "     Compute chisq:"<<endl;
	FitBAOScale fitbao(power_spectrum, su, zref, sig8, n);
	fitbao.ComputeChisq(maxk);
	

	// Find best fit scale and 1-sig error
	cout << "     Find best-fit scale and 1-sig error:"<<endl;
	double bestfit, siglow, sighigh;
	int nsig = 1;
	fitbao.BestfitStdDev(bestfit, siglow, sighigh, nsig);
	double errup = sighigh - bestfit;
	double errdown = bestfit - siglow;
	cout <<"      ka = "<< bestfit <<"+"<< errup <<"-"<< errdown <<endl;
	cout <<endl;
	

	// print info to a file
	cout << "     Print chisq and results to files"<<endl;
	string outfile;
	outfile = outfile_root + "_chisq.txt";
	cout << "     Write chi^2 to file "<< outfile <<endl;
	fitbao.WriteChisq(outfile);
	outfile = outfile_root + "_ancillary.txt";
	cout << "     Write reference power spectrum AND best-fit sinusoid model";
	cout << " to file "<< outfile <<endl;
	fitbao.WriteAncillaryInfo(outfile);
	outfile = outfile_root + "_result.txt";
	cout << "     Write results to file "<< outfile <<endl;
	fitbao.WriteResults(outfile);
	cout << endl;


    }// end of try
  
  
catch(PThrowable exc ) {
    cerr << "fitkbao.cc , Catched exception: \n" << exc.what() << endl;
    }
catch(std::exception ex) {
    cerr << "fitkbao.cc , Catched exception ! " << (string)(ex.what()) << endl;
    }
catch(...) {
    cerr << "fitkbao.cc , Catched ... ! " << endl;
    }

cout << "--------------- fitkbao.cc / END --------------------- " << endl;
}// end of main
