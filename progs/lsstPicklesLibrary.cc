#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <typeinfo>
#include "timing.h"

#include "array.h"
#include "fiosinit.h"
#include "mydefrg.h"

#include "constcosmo.h"
#include "cosmocalcs.h"
#include "simdata.h"
#include "sedfilter.h"


int main(int narg, char* arg[])
{
  cout << " ==== lsstPicklesLibrary.cc program , compute all LSST colors for all stars in";
 cout << " Pickles' Library ====" << endl << endl;
  if (narg < 2) {
        cout << " Usage: lsstPicklesLibrary outfile" << endl;
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
	
	
	// Read in the LSST filters ugrizy
	cout <<"   lsstPicklesLibrary: Loading in filters ..."<<endl;
	// these filter files must be in the directory pointed to by the 
	// environment variable FILTLOC (type echo $FILTLOC at terminal to 
	// check this)

	string filterFile = "LSST.filters";
	ReadFilterList lsstFilters(filterFile);
	lsstFilters.readFilters(lmin,lmax);
	vector<Filter*> filterArray=lsstFilters.getFilterArray();
	int nFilters=lsstFilters.getNTot();
	cout <<"     "<<nFilters<<" LSST filters read in "<<endl;

    // LSST filter indices (computing u-g, g-r, r-i, i-z, z-y colors)
	int uLSST = 0;
	int gLSST = 1;
	int rLSST = 2;
	int iLSST = 3;
	int zLSST = 4;
	int yLSST = 5;

		
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


	// Calculate all LSST colors of all 131 stars at z=0 

	TVector<r_8> Cug(nsed), Cgr(nsed), Cri(nsed), Ciz(nsed), Czy(nsed);
	double z=0.0; 

	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();

	for (int j=0; j<nsed; j++)
	{

		// "Color" in astronomy means ratio of flux in filter X compared 
		// to filter Y
		// This is equivalent to the difference in magnitudes in filter
		// X and Y 
		// because: mag ~ log10(flux)
		

		// mag in filter u - mag in filter g 

			Cug(j)=photometryCalcs.CompColor(z,(*sedArray[j]),
								(*filterArray[uLSST]),(*filterArray[gLSST]));

		// mag in filter g - mag in filter r 

			Cgr(j)=photometryCalcs.CompColor(z,(*sedArray[j]),
								(*filterArray[gLSST]),(*filterArray[rLSST]));

		// mag in filter r - mag in filter i 

			Cri(j)=photometryCalcs.CompColor(z,(*sedArray[j]),
								(*filterArray[rLSST]),(*filterArray[iLSST]));

		// mag in filter i - mag in filter z 

			Ciz(j)=photometryCalcs.CompColor(z,(*sedArray[j]),
								(*filterArray[iLSST]),(*filterArray[zLSST]));

		// mag in filter z - mag in filter y 

			Czy(j)=photometryCalcs.CompColor(z,(*sedArray[j]),
								(*filterArray[zLSST]),(*filterArray[yLSST]));


           }
		

	// Write colors to a file


	if(inp.fail()) {
	
		inp.clear(ios::failbit);
		cout << "Writing to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	for (int j=0; j<nsed; j++)
	{


//                      outp << j <<"    "<<Cug(j)<<"    "<<Cgr(j)<<"    "<<Cri(j)<<"    "<<Ciz(j)<<"    "<<Czy(j)<<endl;
                      outp << sedNames[j] <<"    "<<Cug(j)<<"    "<<Cgr(j)<<"    "<<Cri(j)<<"    "<<Ciz(j)<<"    "<<Czy(j)<<endl;

	}
		outp.close();
	  	
	}


	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;
	  

  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " lsstPicklesLibrary.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " lsstPicklesLibrary.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " lsstPicklesLibrary.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of lsstPicklesLibrary.cc program  Rc= " << rc << endl;
  return rc;	
}

	
