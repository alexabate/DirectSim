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
	
    cout <<"  Given a SED set and filter set calculates the rest-frame colors of all the SEDs         "<<endl;
    cout << endl;
	
	cout <<"  Example usage:                                                                          "<<endl;
	cout <<"  getSEDcolors -o output -t CWWKSB.list -f SDSS.filters                                   "<<endl;
	cout << endl;
	
	cout << " -o: OUTFILE: write SED colors to file                                                   "<<endl;
	cout << " -t: SEDLIB: file containing list of SED files                                           "<<endl;
	cout << " -f: FILTERS: file containing list of filters                                            "<<endl; 
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

	char c;
    while((c = getopt(narg,arg,"ho:t:f:")) != -1) {
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
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }

    //-- end command line arguments
    cout <<"     Writing SED's and SED colors to file "<< outroot <<endl;
    cout <<"     Using SED library "<< sedfile << endl;
    cout <<"     Using filter set "<< filtfile << endl;
    cout << endl;
    //-- end command line arguments
  
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
	ResourceUsage res;
	InitTim();
	string outfile;

		
	// GALAXY SED TEMPLATES
	
	// wavelength range of the SEDs/filters
	double lmin=5e-8, lmax=2.5e-6;
	
	ReadSedList readSedList(sedfile);
	
	// Read out SEDs into array
    readSedList.readSeds(lmin,lmax);
    int nsedOrig = readSedList.getNSed();
    cout <<"     Number of original SEDs = "<< nsedOrig <<endl;
	cout << endl;
	
	// Get total number of SEDs
    vector<SED*> sedArray = readSedList.getSedArray();
    cout <<"     Number of SEDs in SED array "<< sedArray.size() << endl;
    cout << endl;
    
    // write out SEDs 
    outfile = outroot + "_SEDlib.txt";
    readSedList.writeSpectra(outfile);
    
    
    // FILTERS
	
	ReadFilterList filts(filtfile);
	filts.readFilters(lmin, lmax);
	vector<Filter*> filters = filts.getFilterArray();
	int nFilters = filts.getNTot();
	cout <<"     "<< nFilters <<" filters read in "<<endl;
    cout << endl;
	
	
	// GENERATE SED COLORS
	
	// for fitting SEDs to colors
	SEDLibColors sedLibColors(sedArray, filters);
	outfile = outroot + "_colors.txt";
	sedLibColors.writeColorArray(outfile);
	
	
	  
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
