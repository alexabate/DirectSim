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

#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "simdata.h"
#include "sedfilter.h"

void usage(void);
void usage(void) {
	cout << endl<<" Usage: getSEDcolors [...options...]" << endl<<endl;
	
    cout <<"  Given a SED set and filter set calculates the rest-frame colors (or mags) of all the SEDs"<<endl;
    cout << endl;
    cout <<"  You might want to output magnitudes instead of colors if filter list supplied is a list  "<<endl;
    cout <<"  of filters from different filter sets (e.g. SDSS, LSST and CFHT). In this case a         "<<endl;
    cout <<"  'normalising' magnitude in one of the filters must be supplied, along with the index of  "<<endl;
    cout <<"  this 'normalising' filter. (filter indices are zero indexed).                            "<<endl;
    cout << endl;
	
	cout <<"  Example usage (to output colors):                                                       "<<endl;
	cout <<"  getSEDcolors -o output -t CWWKSB.list -f SDSS.filters                                   "<<endl;
	cout << endl;
	
	cout <<"  Example usage (to output mags normalized to 22 in i band which is on line 4 (index 3) of"<<endl; 
	cout <<"  the list in random.filters):                                                            "<<endl;
	cout <<"  getSEDcolors -o output -t CWWKSB.list -f random.filters -m 22,3                         "<<endl;
	cout << endl;
	
	cout << " -o: OUTFILE: write SED colors (or mags if -m option used) to this file                  "<<endl;
	cout << " -t: SEDLIB: file containing list of SED files                                           "<<endl;
	cout << " -f: FILTERS: file containing list of filters                                            "<<endl; 
	cout << " -m: MAGNORM,FILTNORM: output mags normalised to MAGNORM in filter FILTNORM              "<<endl;
	cout << " -w: LMIN,LMAX,NL: wavelength resolution of the SEDs                                     "<<endl;
	cout << endl;
    };


int main(int narg, char* arg[]) {
    
    cout << " ==== getSEDcolors.cc program , to calculate colors of SED set ====" << endl;
    cout << endl;

    // make sure SOPHYA modules are initialized 
    SophyaInit();  
    FitsIOServerInit();
    InitTim();
    cout<<endl<<endl;

    //--- decoding command line arguments 
    string outroot;
    string sedfile = "CWWKSB.list";
    string filtfile = "SDSS.filters";
    bool outColors = true;
    double magNorm;
    int iFiltNorm;
	// wavelength range of the SEDs/filters
	double lmin=1e-8, lmax=2.5e-6; // need to be careful with this given what filters are read in!
	int npt = 1000; // interpolation resolution of SEDs (need more if many fine features)

	char c;
    while((c = getopt(narg,arg,"ho:t:f:m:w:")) != -1) {
	    switch (c)  {
	        case 'o' :
	            outroot = optarg;
	            break;
	        case 't' :
	            sedfile = optarg;
	            break;
	        case 'f' :
	            filtfile = optarg;
	            break;
	        case 'm' :
	            outColors = false;
	            sscanf(optarg,"%lf,%d",&magNorm,&iFiltNorm);
	            break;
	        case 'w' :
	            sscanf(optarg,"%lf,%lf,%d",&lmin,&lmax,&npt);
	            break;
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }

    //-- end command line arguments
    cout <<"     Writing SED's and SED";
    if (outColors)
        cout <<" colors to file "<< outroot <<endl;
    else {
        cout <<" mags to file "<< outroot <<", normalized to "<< magNorm <<" in filter indexed by ";
        cout << iFiltNorm <<endl;
        }
    cout <<"     Using SED library "<< sedfile <<" with wavelength range (in m) "<< lmin <<" to "<< lmax;
    cout <<" and resolution "<< npt << endl;
    cout <<"     Using filter set "<< filtfile << endl;
    
    cout << endl;
    //-- end command line arguments
  
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
	ResourceUsage res;
	InitTim();
	string outfile;


	// GALAXY SED TEMPLATES
	ReadSedList readSedList(sedfile);
	
	// Read out SEDs into array
    readSedList.readSeds(lmin, lmax, npt);
    int nsedOrig = readSedList.getNSed();
    cout <<"     Number of original SEDs = "<< nsedOrig <<endl;
	cout << endl;
	
	// Get total number of SEDs
    vector<SED*> sedArray = readSedList.getSedArray();
    cout <<"     Number of SEDs in SED array "<< sedArray.size() << endl;
    cout << endl;
    
    // write out SEDs 
    outfile = outroot + "_SEDlib.txt";
    readSedList.writeSpectra(outfile, lmin, lmax, npt);
    
    
    // FILTERS
	
	ReadFilterList filts(filtfile);
	filts.readFilters(lmin, lmax);
	vector<Filter*> filters = filts.getFilterArray();
	int nFilters = filts.getNTot();
	cout <<"     "<< nFilters <<" filters read in "<<endl;
    cout << endl;
	
	
	// GENERATE SED MAGS or COLORS
	outfile = outroot + "_restframedata.txt";
	
	Timer tm("timer",false);
	SEDLibColors sedLibColors(sedArray, filters, lmin, lmax, npt);
	if (outColors) {
	    // for fitting SEDs to colors
	    sedLibColors.writeColorArray(outfile);
	    }
	else {
	    // for fitting SEDs to mags
	    sedLibColors.writeMagsArray(outfile, magNorm, iFiltNorm);
	    }
	tm.Split();
    cout <<"     Color computations took "<< tm.PartialElapsedTimems() <<" ms "<< endl;
	
	  
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " getSEDcolors.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " getSEDcolors.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " getSEDcolors.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of getSEDcolors.cc program  Rc= " << rc << endl;
  return rc;	
}
