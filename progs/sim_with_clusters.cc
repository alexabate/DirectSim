/**
  * @file  sim_with_clusters.cc
  * @brief Read over-density grid and simulate galaxy catalog, with clusters
  *
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
#include "resusage.h"

// DirectSim
#include "mydefrg.h"
#include "mass2gal.h"
#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "constcosmo.h"

void usage(void);
void usage(void) {

	cout << endl<<" Usage: sim_with_clusters [...options...]      "<<endl<<endl;


	cout << endl;
    }

int main(int narg, char* arg[]) {

	cout << " ==== sim_with_clusters.cc program , simulating galaxy catalog "<<endl;
	cout << " with clusters from input over-density grid  ==== " << endl;
	
	// Make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;

	string infile;            // name of over-density grid to read in
	string outfileroot;       // name of file to output to

	//--- decoding command line arguments 
	char c;
	while((c = getopt(narg,arg,"hi:o:")) != -1) {
	    switch (c) {
            case 'i' :
                infile = optarg;
                break;
            case 'o' :
                outfileroot	= optarg;
                break;
            case 'h' :
                default :
                usage(); return -1;
		    }
	    }
    cout <<"     Reading in over density grid from "<< infile << endl;
    cout <<"     Writing out galaxy catalog with clusters to "<< outfileroot << endl;
    cout << endl;
	
		
    //-- end command line arguments
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    
    // Random generator
    RandomGenerator rg;
	
  
    // Read over-density grid
	cout << "     Reading input file= " << infile << endl;  
	FitsInOutFile fin(infile,FitsInOutFile::Fits_RO);
   
   
    // Initialize cosmological parameters 
	cout << "     Initialise cosmology: "<<endl;
	double h = 0.71, OmegaM = 0.3, OmegaL = 0.7;
	SimpleUniverse su(h, OmegaM, OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout << "     OmegaK="<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
	cout << ", OmegaL="<< su.OmegaLambda() <<", OmegaB="<< su.OmegaBaryon();
	cout << ", H0="<< su.H0() <<endl;
	cout << endl;
	
	
	// Initialize FieldClusterGals and remove N extra planes 
	FieldClusterGals fieldClusterGals(fin, su, rg);
	double pixVol = fieldClusterGals.ReturnPixVol();
	double zref = fieldClusterGals.ReturnZref();
	
	/*TArray<r_8> mass;
	fieldClusterGals.MassArray(mass);
	int cnt=0;
    for (int i=0; i<mass.SizeX(); i++)
        for (int j=0; j<mass.SizeY(); j++)
            for (int k=0; k<mass.SizeZ(); k++)
                if ( (mass(i,j,k)-1)>1.6)
                        cnt++;
    cout << cnt <<" elements >1.6"<<endl;*/
	
	
	// Find mass to gal conversion
	cout << "     Convert rho/rho^bar To Mean NGal"<<endl;  
	// LF parameters: define evolution of LF as a function of redshift
	cout << "     Initializing LF parameters as a function of redshift"<<endl;
	double zmax = 5., dz = 0.01;
	LFParameters lfpars(zmax, dz);
	double ps, ms, a; // get LF parameters at redshift of cube
	lfpars(1., ps, ms, a);
	Schechter schechter(ps, ms, a);
	double schmin = -24, schmax = -13;// units of "M-5log10h70"
	int schnpt = 10000;
	schechter.SetInteg(schmin,schmax,schnpt);
	double numz = schechter.Integrate();
    double conv = numz*pixVol;
	cout <<"     gals per pixel vol = "<< conv << endl;
	cout << endl;
	
    /*TArray<r_8> odens_array = fieldClusterGals.ODensArray();
	cnt=0;
    for (int i=0; i<odens_array.SizeX(); i++)
        for (int j=0; j<odens_array.SizeY(); j++)
            for (int k=0; k<odens_array.SizeZ(); k++)
                if ( (odens_array(i,j,k))>1.6)
                        cnt++;
    cout << cnt <<" elements >1.6"<<endl;*/
    
    double ns = 1., sig8 = 0.8, m1 = 5e13, m2 = 1e15;
    int nstep = 6;
    double z = 1.;
    double delt = fieldClusterGals.findClusterDelta(ns, sig8, m1, m2, nstep, z);
	cout << "     clusters have delta > "<< delt << endl;
	
	// Convert mass in each cell to a (mean) number of galaxies
	double bias = 10.; // relative bias of cluster vs field galaxies
	fieldClusterGals.simulateGalaxies(conv, bias, outfileroot); 

	
	
	

  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " sim_with_clusters.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " sim_with_clusters.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " sim_with_clusters.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of sim_with_clusters.cc program  Rc= " << rc << endl;
  return rc;	
}
