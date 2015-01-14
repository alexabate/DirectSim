// -*- LSST-C++ -*-
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
#define PI 3.141592

void usage(void);
void usage(void)
{
	cout << endl<<" Usage: testErrors [...options...]" << endl<<endl;

	cout << " -o : OUTFILE: output filename (will be stored in testfiles/)"<<endl;
	cout << endl;
}
int main(int narg, char* arg[])
{
	cout << " ==== testErrors.cc program , to test data ==== "<<endl;

	// make sure SOPHYA modules are initialized 
    SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;

	//--- decoding command line arguments 
	string outfile="testfiles/testErrors.fits";
	string sedFile = "CWWK.list";
    int nElliptical = 1;
    int nSpiral = 2;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hi:t:z:m:o:")) != -1)  {
	    switch (c) {
	        case 'o' :
		        outfile = optarg;
		        outfile = "testfiles/" + outfile;
		        break;
	      case 'h' :
		    default :
		    usage(); return -1;
		    }
	    }
	    
	int nYear = 10;
	// Number of visits per year (Table 1, Ivezic et al 2008)
    /*int uVisitsPerYear = 6;
    int gVisitsPerYear = 8;
    int rVisitsPerYear = 18;
    int iVisitsPerYear = 18;
    int zVisitsPerYear = 16;
    int yVisitsPerYear = 16;*/
    
	// total number of visits
	/*int uVisits = uVisitsPerYear*nYear;
    int gVisits = gVisitsPerYear*nYear;
    int rVisits = rVisitsPerYear*nYear;
    int iVisits = iVisitsPerYear*nYear;
    int zVisits = zVisitsPerYear*nYear;
    int yVisits = yVisitsPerYear*nYear;*/
	
    cout << "     Number of ellipticals = "<< nElliptical <<", number of spirals = ";
    cout << nSpiral << endl;
    cout << "     Output catalog is "<< outfile <<endl;
    cout << endl;
    //-- end command line arguments
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
	ResourceUsage res;
	InitTim();
	
	// Output file
	cout <<"     Creating output file "<<outfile<<endl;
	FitsInOutFile swf(outfile,FitsInOutFile::Fits_Create);
	
	// binary data table
	SwFitsDataTable gals(swf,2048);
	gals.AddFloatColumn("mt");
	gals.AddFloatColumn("muo");
	gals.AddFloatColumn("mgo");
	gals.AddFloatColumn("mro");
	gals.AddFloatColumn("mio");
	gals.AddFloatColumn("mzo");
	gals.AddFloatColumn("myo");
	gals.AddFloatColumn("euo");
	gals.AddFloatColumn("ego");
	gals.AddFloatColumn("ero");
	gals.AddFloatColumn("eio");
	gals.AddFloatColumn("ezo");
	gals.AddFloatColumn("eyo");
	DataTableRow rowin=gals.EmptyRow();
	cout << endl;
	
	// Cosmological parameters
	cout <<"     Set cosmology ..."<<endl;
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout << "     h = "<<h<<", OmegaM = "<<OmegaM<<", OmegaL = "<<OmegaL;
	cout <<endl<<endl;

	// This controls the drawing of random numbers
	RandomGenerator rg;
	
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
	cout <<endl;
	
	// Load in SEDs
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax); // Read out SEDs into array
    int nsed=readSedList.getNSed(); // Get total number of SEDs  
    int ntot = readSedList.getNTot();
    vector<SED*> sedArray=readSedList.getSedArray();
    cout <<"     Final total number of SEDs = "<<ntot<<endl;
	cout << endl;

	
	// number of galaxies in the survey
	long ng = 10000;
	cout<<endl;
		
	
	// Prepare the class which will calculate the magnitudes
	cout <<"     Initialize class to calculate magnitudes"<<endl;
	SimData simgal(sedArray,LSSTfilters,su,rg,nElliptical,nSpiral);

    // Loop over all galaxies in the base catalog
	cout <<"     Start loop over galaxies ..."<<endl;
	
	double mMin = 15, mMax = 40;
	int nM = 10000;
	double dm = (mMax - mMin)/(nM - 1);

	for (int i=0; i<ng; i++)
		{
		cout <<"     gal "<<i+1<<" of "<<ng<<endl;

        double mag = mMin + i*dm;
		
		vector<double> uObservation = simgal.addLSSTError(mag, nYear, 0);
        vector<double> gObservation = simgal.addLSSTError(mag, nYear, 1);
        vector<double> rObservation = simgal.addLSSTError(mag, nYear, 2);
        vector<double> iObservation = simgal.addLSSTError(mag, nYear, 3);
        vector<double> zObservation = simgal.addLSSTError(mag, nYear, 4);
        vector<double> yObservation = simgal.addLSSTError(mag, nYear, 5);
		
        // Write the data to the FITS file
		rowin[0]=mag;
		rowin[1]=uObservation[0];
		rowin[2]=gObservation[0];
		rowin[3]=rObservation[0];
		rowin[4]=iObservation[0];
		rowin[5]=zObservation[0];
		rowin[6]=yObservation[0];
		rowin[7]=uObservation[1];
		rowin[8]=gObservation[1];
		rowin[9]=rObservation[1];
		rowin[10]=iObservation[1];
		rowin[11]=zObservation[1];
		rowin[12]=yObservation[1];
		gals.AddRow(rowin);
		}
	cout <<"     End loop"<<endl;
	
	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testErrors.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testErrors.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " testErrors.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of testErrors.cc program  Rc= " << rc << endl;
  return rc;	
}
