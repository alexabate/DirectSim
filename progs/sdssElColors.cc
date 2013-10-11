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
  cout << " ==== sdssElColors.cc program , compute SDSS colors of an elliptical";
  cout << " galaxy ====" << endl << endl;
  if (narg < 2) {
        cout << " Usage: sdssElColors outfile" << endl;
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
	
	
	
	// Read in the SDSS filters ugriz
	cout <<"   sdssElColors: Loading in filters ..."<<endl;
	// these filter files must be in the directory pointed to by the 
	// environment variable FILTLOC (type echo $FILTLOC at terminal to 
	// check this)

	string filterFile = "SDSS.filters";
	ReadFilterList sdssFilters(filterFile);
	sdssFilters.readFilters(lmin,lmax);
	vector<Filter*> filterArray=sdssFilters.getFilterArray();
	int nFilters=sdssFilters.getNTot();
	cout <<"     "<<nFilters<<" SDSS filters read in "<<endl;

    // SDSS filter indices (just going to compute g-r and r-i colors)
	int gSDSS = 1;
	int rSDSS = 2;
	int iSDSS = 3;
	
	
		
	// Read in galaxy SED templates (El,Sbc,Scd,Irr,SB3,SB2)
	// these SED files must be in the directory pointed to by the 
	// environment variable SEDLOC (type echo $SEDLOC at terminal to check 
	// this)
	
	// "CWW" stands for Coleman, Wu & Weedman 1980
	// in that paper they measured spectra of typical elliptical (early), 
	// spiral (late, Sbc, Scd) and irregular type galaxies
	// "K" stands for Kinney et al 1996 who measure starburst galaxies
	string sedFile = "CWWK.list";
	ReadSedList cwwSEDs(sedFile);
	// Read out SEDs into array
    cwwSEDs.readSeds(lmin,lmax);
    // Get total number of SEDs
    vector<SED*> sedArray=cwwSEDs.getSedArray();
    int nsed=cwwSEDs.getNSed();
    cout <<"     Number of original SEDs = "<<nsed<<endl;
	cout << endl;
	
	// sed index
	int iEl = 0;    // elliptical galaxy

	
	
    // Photometry calculation class
    PhotometryCalcs photometryCalcs(lmin,lmax);
    

	// Loop over redshifts, calculate color of redshifted elliptical spectrum
	double zmin=0, zmax=2.3;
	int nz=100;
	double dz=(zmax-zmin)/(nz-1);
	
	TVector<r_8> Cgr(nz), Cri(nz);
	for (int i=0; i<nz;i++){
	
		double z=zmin+i*dz;

		// "Color" in astronomy means ratio of flux in filter X compared 
		// to filter Y
		// This is equivalent to the difference in magnitudes in filter
		// X and Y 
		// because: mag ~ log10(flux)
		
		// mag in filter g - mag in filter r 
		Cgr(i)=photometryCalcs.CompColor(z,(*sedArray[iEl]),
								(*filterArray[gSDSS]),(*filterArray[rSDSS]));
		// mag in filter r - mag in filter i 
		Cri(i)=photometryCalcs.CompColor(z,(*sedArray[iEl]),
								(*filterArray[rSDSS]),(*filterArray[iSDSS]));
		
		}
		
		
		
	// Write colors to a file
	// Compare Cgr vs Cri to figure 2 in Padmanabhan et al (2005) (diamonds in the top panel)
	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
	
		inp.clear(ios::failbit);
		cout << "Writing to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (int i=0;i<nz;i++) {		
			double z=zmin+i*dz;
			outp << z <<"    "<< Cgr(i) <<"    "<< Cri(i) <<endl;
		  	}
		outp.close();
	  	}
	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;
	  

  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " sdssElColors.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " sdssElColors.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " sdssElColors.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of sdssElColors.cc program  Rc= " << rc << endl;
  return rc;	
}
