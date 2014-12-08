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


int main(int narg, char* arg[]) {
  cout << " ==== generateTemplSet.cc program , to generate reddened template set ====" << endl << endl;
  if (narg < 2) {
        cout << " Usage: generateTemplSet outfile" << endl;
        cout << "        outfile: output filename "<<endl;
        cout << endl;
        cout <<endl;
        return 1;
        }

  // make sure SOPHYA modules are initialized 
  SophyaInit();  
  FitsIOServerInit();
  InitTim();
  cout<<endl<<endl;

  //--- decoding command line arguments 
  string outfile = arg[1];
  //-- end command line arguments
  
  int rc = 1;  
  try {  // exception handling try bloc at top level
  
	ResourceUsage res;
	InitTim();
	
	// wavelength range of the SEDs/filters
	double lmin=5e-8, lmax=2.5e-6;
	
		
	// GALAXY SED TEMPLATES
	string sedFile = "CWWKSB.list";
	ReadSedList readSedList(sedFile);
	
	// Read out SEDs into array
    readSedList.readSeds(lmin,lmax);
    
    // Redden SEDs 
    int nPerSED = 4;
    int method = 0; // uniform distribution up to max
    int maxidEl = 0;
    double redMaxEl = 0.1;
    double redMaxOther = 0.3;
    vector<double> reds = readSedList.reddenSeds(nPerSED, method, maxidEl, redMaxEl, redMaxOther);
    
    
    int nsedOrig = readSedList.getNSed();
    int nsedTot = readSedList.getNTot();
    cout <<"     Number of original SEDs = "<< nsedOrig <<endl;
    cout <<"     Number of SEDs after adding reddened ones = "<< nsedTot <<endl;
	cout << endl;
	
	
	// Get total number of SEDs
    vector<SED*> sedArray = readSedList.getSedArray();
    cout <<"     Number of SEDs in SED array "<< sedArray.size() << endl;
    cout << endl;
	
	
    // Write SEDs to a file
    readSedList.writeSpectra(outfile);
	

	  
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " generateTemplSet.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " generateTemplSet.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " generateTemplSet.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of generateTemplSet.cc program  Rc= " << rc << endl;
  return rc;	
}
