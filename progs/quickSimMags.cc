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
#include "resusage.h"
#include "timestamp.h"
#include "ctimer.h"

#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "simdata.h"
#include "geneutils.h"

void usage(void);
void usage(void)
{
	cout << endl<<" Usage: quickSimMags -o [outfile]" << endl<<endl;
    
    cout <<"    Simple program to generate some observed magnitudes (LSST photometry) "<<endl;
    cout <<"    Output file can be read into BPZ for photo-z calculation"<<endl;
    cout <<"    Magnitudes are generated from only one galaxy type"<<endl;

	cout << " -o : OUTFILE: output filename for observed catalog (in format ready for BPZ)"<<endl;
	cout << endl;
}
int main(int narg, char* arg[])
{
	cout << " ==== quickSimMags.cc program , to simulate LSST data ==== "<<endl;

	// make sure SOPHYA modules are initialized 
    SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;

	//--- decoding command line arguments 
	string sedFile = "CWWK.list";
	string outfile;

    
    // Number of visits per year (Table 1, Ivezic et al 2008)
    /*int uVisitsPerYear = 6;
    int gVisitsPerYear = 8;
    int rVisitsPerYear = 18;
    int iVisitsPerYear = 18;
    int zVisitsPerYear = 16;
    int yVisitsPerYear = 16;*/
    int nYear = 10;
    
    // don't redden SEDs
    bool isAddRedden = false;
    // do include Madau absorption
    bool isAddMadau = true;
    
    // Number of elliptical, spiral galaxy types in CWWK template set
    int nElliptical = 1;
    int nSpiral = 2;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ho:")) != -1)  {
	    switch (c) {
	        case 'o' :
		        outfile = optarg;
		        break;
	      case 'h' :
		    default :
		    usage(); return -1;
		    }
	    }
	
	// total number of visits
	/*int uVisits = uVisitsPerYear*nYear;
    int gVisits = gVisitsPerYear*nYear;
    int rVisits = rVisitsPerYear*nYear;
    int iVisits = iVisitsPerYear*nYear;
    int zVisits = zVisitsPerYear*nYear;
    int yVisits = yVisitsPerYear*nYear;*/
	  
    cout << "     Output magnitudes to "<< outfile <<endl;
    cout << endl;
    //-- end command line arguments
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
	ResourceUsage res;
	InitTim();
	
	// this controls the drawing of random numbers
	RandomGenerator rg;

	// Set cosmology: required for calculating observed magnitudes
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout <<"     Set cosmology to: OmegaM="<<OmegaM<<", OmegaL="<<OmegaL;
	cout <<", H0="<<100*h<<endl;
	cout << endl;


	// Load in filters required
	// wavelength range of the SEDs/filters
	cout <<"     Load in LSST filters"<<endl;
	double lmin=5e-8, lmax=2.5e-6;
	string filterFile = "LSST.filters";
	ReadFilterList readFilterList(filterFile);
	readFilterList.readFilters(lmin,lmax);
	vector<Filter*> LSSTfilters=readFilterList.getFilterArray();
	int nFilter = readFilterList.getNTot();
	cout <<"     Read in "<< nFilter <<" filters"<<endl;
	//Filter restFrameFilter((*GOODSfilters[1]));
	int iU = 0;
	int iG = 1;
	int iR = 2;
	int iI = 3;
	int iZ = 4;
	int iY = 5;
	cout <<endl;
	
	
	// The rest-frame filter in which the absolute magnitude is defined
	cout <<"     Load in GOODS B filter"<<endl;
	string goodsFilterFile = "GOODSB.filters";
	ReadFilterList readGOODSBfilter(goodsFilterFile);
	readGOODSBfilter.readFilters(lmin,lmax);
	vector<Filter*> goodsBFilter=readGOODSBfilter.getFilterArray();
	
	
	// Load in SEDs
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax); // Read out SEDs into array
    int nsed=readSedList.getNSed(); // Get total number of SEDs
    int ntot = readSedList.getNTot();
    vector<SED*> sedArray=readSedList.getSedArray();
    cout <<"     Number of SEDs = "<< ntot <<endl;
	cout << endl;


	// Prepare the class which will calculate the magnitudes
	cout <<"     Initialize class to calculate magnitudes"<<endl;
	SimData simgal(sedArray, LSSTfilters, su, rg, nElliptical, nSpiral);
    cout << endl;
    
    
    // Add Madau preference
    simgal.setMadau(isAddMadau);
    
    
    // Galaxies to simulate
    int ng = 1000;
    double am = -21.;
    double zmin = 0.; 
    double zmax = 3.;
    double dz = (zmax-zmin)/(ng-1.);
    double type = 1.000; // this means elliptical
    double ext = 0.; // don't add reddening


    // Loop over all galaxies in the base catalog
	// open output file
    ofstream outp;
    cout <<"     Opening file "<< outfile << endl;
    outp.open(outfile.c_str(), ofstream::out);
    
	for (int i=0; i<ng; i++) {

		
		// Get galaxy properties: galaxy has redshift zs, absolute magnitude am,
		// SED type type and internal extinction ext.
		double zs = zmin + i*dz;


        // Calculate galaxy magnitude in observed filters ugrizy
        // Galaxy's absolute magnitude is defined in filter (*goodsBFilter[0])
		double uMagTh=simgal.GetMag(zs,type,am,ext,iU,(*goodsBFilter[0]));
		double gMagTh=simgal.GetMag(zs,type,am,ext,iG,(*goodsBFilter[0]));
		double rMagTh=simgal.GetMag(zs,type,am,ext,iR,(*goodsBFilter[0]));
		double iMagTh=simgal.GetMag(zs,type,am,ext,iI,(*goodsBFilter[0]));
		double zMagTh=simgal.GetMag(zs,type,am,ext,iZ,(*goodsBFilter[0]));
		double yMagTh=simgal.GetMag(zs,type,am,ext,iY,(*goodsBFilter[0]));
		
		//cout << zs <<"  ";
		//cout << uMagTh <<"  "<< gMagTh <<"  "<< rMagTh <<"  "<< iMagTh <<"  "<< zMagTh <<"  "<< yMagTh <<endl;
		
		// The final observations
		// The 1st element is the value of the observed magnitude
		// The 2nd element is the magnitude error
		vector<double> uObservation = simgal.addLSSTError(uMagTh, nYear, iU);
		vector<double> gObservation = simgal.addLSSTError(gMagTh, nYear, iG);
		vector<double> rObservation = simgal.addLSSTError(rMagTh, nYear, iR);
		vector<double> iObservation = simgal.addLSSTError(iMagTh, nYear, iI);
		vector<double> zObservation = simgal.addLSSTError(zMagTh, nYear, iZ);
        vector<double> yObservation = simgal.addLSSTError(yMagTh, nYear, iY);

        // Write the data to the file: zs, ugrizy, eugrizy
        outp << zs <<"  "<< uObservation[0] <<"  "<< gObservation[0] <<"  "<< rObservation[0];
        outp <<"  "<< iObservation[0] <<"  "<< zObservation[0] <<"  "<< yObservation[0];
        outp <<"  "<< uObservation[1] <<"  "<< gObservation[1] <<"  "<< rObservation[1];
        outp <<"  "<< iObservation[1] <<"  "<< zObservation[1] <<"  "<< yObservation[1] << endl;
        
        //if (my_isnan(uObservation[0]))
        //    cout <<" u is nan"<<endl;
        
        
		}
    outp.close();

	
	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " quickSimMags.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " quickSimMags.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " quickSimMags.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of quickSimMags.cc program  Rc= " << rc << endl;
  return rc;	
}
