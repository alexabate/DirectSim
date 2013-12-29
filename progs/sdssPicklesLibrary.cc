#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <ostream>
#include <cstdlib>
#include <cstdio>
#include <iomanip>

#include <typeinfo>
#include "timing.h"

#include "array.h"
#include "fiosinit.h"
#include "mydefrg.h"


#include "constcosmo.h"
#include "cosmocalcs.h"
#include "simdata.h"
#include "sedfilter.h"

//using namespace std; 

int main(int narg, char* arg[])
{
  cout << " ==== sdssPicklesLibrary.cc program , compute SDSS r-i colors for all stars in";
 cout << " Pickles' Library ====" << endl << endl;
  if (narg < 2) {
        cout << " Usage: sdssPicklesLibrary outfile" << endl;
	cout << "        outfile: output filename"<<endl;
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
  string outfile=arg[1];
  
  //-- end command line arguments
  
  int rc = 1;  
  try {  // exception handling try bloc at top level
	InitTim();

	// wavelength range of the SEDs/filters
	double lmin=5e-8, lmax=2.5e-6;
	
//	std::string output_filename; 	
	
	
	// Read in the SDSS filters ugriz
	cout <<"   sdssPicklesLibrary: Loading in filters ..."<<endl;
	// these filter files must be in the directory pointed to by the 
	// environment variable FILTLOC (type echo $FILTLOC at terminal to 
	// check this)

	string filterFile = "SDSS.filters";
	ReadFilterList sdssFilters(filterFile);
	sdssFilters.readFilters(lmin,lmax);
	vector<Filter*> filterArray=sdssFilters.getFilterArray();
	int nFilters=sdssFilters.getNTot();
	cout <<"     "<<nFilters<<" SDSS filters read in "<<endl;

    // SDSS filter indices (just going to compute r-i colors)
//	int gSDSS = 1;
	int rSDSS = 2;
	int iSDSS = 3;

		
	// Read in Pickle's Library 131 star SEDs (UVILIB or UVKLIB) 
	// these SED files must be in the directory pointed to by the 
	// environment variable SEDLOC (type echo $SEDLOC at terminal to check 
	// this)
//	string sedFile = "UVILIB.list";
	string sedFile = "UVKLIB.list";
	ReadSedList picklesSEDs(sedFile);
	// Read out SEDs into array
    picklesSEDs.readSeds(lmin,lmax);
    // Get total number of SEDs
    vector<SED*> sedArray=picklesSEDs.getSedArray();
    int nsed=picklesSEDs.getNSed();

    cout <<"     Number of original SEDs = "<<nsed<<endl;
	cout << endl;

	string sedNames[nsed];
	ifstream names;
	names.open("UVKLIB.list"); 
	for (int i=0; i<nsed; i++)
	{
		getline(names, sedNames[i]);
	}
	
    // Photometry calculation class
    PhotometryCalcs photometryCalcs(lmin,lmax);
 
	//calculate SDSS r-i colors of all 131 stars at z=0
	TVector<r_8> Cri(nsed);

	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	

	for (int j=0; j<nsed; j++)
	{

			double z=0.0; 

		// "Color" in astronomy means ratio of flux in filter X compared 
		// to filter Y
		// This is equivalent to the difference in magnitudes in filter
		// X and Y 
		// because: mag ~ log10(flux)

		// mag in filter r - mag in filter i 
			Cri(j)=photometryCalcs.CompColor(z,(*sedArray[j]),
								(*filterArray[rSDSS]),(*filterArray[iSDSS]));
		}


	// Write colors to a file

	if(inp.fail()) {
	
		inp.clear(ios::failbit);
		cout << "Writing to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	for (int j=0; j<nsed; j++)
	{

			outp << sedNames[j] <<"    "<< Cri(j) <<endl;

	}
		outp.close();

	}

	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;
	  

  }  // End of try bloc 
 
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " sdssPicklesLibrary.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " sdssPicklesLibrary.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " sdssPicklesLibrary.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of sdssPicklesLibrary.cc program  Rc= " << rc << endl;
  return rc;	
}

	
