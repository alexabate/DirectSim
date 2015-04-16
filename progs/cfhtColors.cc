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

#define PI 3.141592

int main(int narg, char* arg[])
{
  cout << " ==== cfhtColors.cc program , to CFHT color distributions ====" << endl << endl;
  if (narg < 2) {
        cout << " Usage: cfhtColors outfile" << endl;
	cout << "        outfile: output filename that is stored in testing/ dir"<<endl;
	cout << endl;
	cout << " Calculates u-g and i-z color distributions with host galaxy reddening"<<endl;
	cout << " values of E(B-V)=[0, 0.1, 0.3] " << endl;
	cout <<endl;
    return 1;
  }

  // make sure SOPHYA modules are initialized 
  SophyaInit();  
  FitsIOServerInit();
  InitTim();
  cout<<endl<<endl;

  //--- decoding command line arguments 
  string outfile="testfiles/";
  string fname=arg[1];
  outfile+=fname;
  
  //-- end command line arguments
  
  int rc = 1;  
  try {  // exception handling try bloc at top level
	ResourceUsage res;
	InitTim();
	
	// wavelength range of the SEDs/filters
	double lmin=5e-8, lmax=2.5e-6;
	
	//FILTERS
	cout <<"   cfhtColors: Loading in filters ..."<<endl;
	string filterFile = "CFHT.filters";
	ReadFilterList readCFHTfilters(filterFile);
	readCFHTfilters.readFilters(lmin,lmax);
	vector<Filter*> CFHTfilters=readCFHTfilters.getFilterArray();
	int nFilter=readCFHTfilters.getNTot();
	cout <<"     "<< nFilter <<" CFHT filters read in "<<endl;
    cout << endl;
	
    // filter numbers
	int uCFHT = 0;
	int gCFHT = 1;
	int iCFHT = 3;
	int zCFHT = 4;
	
		
	// GALAXY SED TEMPLATES
	string sedFile = "CWWK.list";
	ReadSedList readSedList(sedFile);
	// Read out SEDs into array
    readSedList.readSeds(lmin,lmax);
    // Get total number of SEDs
    //vector<SED*> sedArray = readSedList.getSedArray();
    int nsed = readSedList.getNSed();
    cout <<"     Number of original SEDs = "<<nsed<<endl;
	cout << endl;
	
	// sed numbers
	int iEl = 0;    // elliptical galaxy
	int iSbc = 1;   // spiral galaxy
	int iScd = 2;   // spiral galaxy
	int iIrr = 3;   // irregular galaxy
	int iSB3 = 4;   // starburst 1
	int iSB2 = 5;   // starburst 2
	
	// Redden the SEDs
	/*SEDRedden ElRed1((*sedArray[iEl]),0.1,0);
	SEDRedden ElRed2((*sedArray[iEl]),0.3,0);
	SEDRedden SbcRed1((*sedArray[iSbc]),0.1,0);
	SEDRedden SbcRed2((*sedArray[iSbc]),0.3,0);
	SEDRedden ScdRed1((*sedArray[iScd]),0.1,0);
	SEDRedden ScdRed2((*sedArray[iScd]),0.3,0);
	SEDRedden IrrRed1((*sedArray[iIrr]),0.1,1);
	SEDRedden IrrRed2((*sedArray[iIrr]),0.3,1);
	SEDRedden SB3Red1((*sedArray[iSB3]),0.1,1);
	SEDRedden SB3Red2((*sedArray[iSB3]),0.3,1);
	SEDRedden SB2Red1((*sedArray[iSB2]),0.1,1);
	SEDRedden SB2Red2((*sedArray[iSB2]),0.3,1);*/
	
	// @todo, this will give two RANDOM reddening values to the SEDs
	// probably should code a fixed reddening version of this method too.
	// also output doesn't record what these reddening values were, and they will be different for El vs others
	int nPerSED = 2;
	int method=0;
	int maxidEl=0;
	double redMaxEl=0.3;
	double redMaxOther=0.3;
	readSedList.reddenSeds(nPerSED, method, maxidEl, redMaxEl, redMaxOther);
	vector<SED*> sedArray = readSedList.getSedArray();
	
	
    // PHOTOMETRY CALCULATION CLASS
    PhotometryCalcs photometryCalcs(lmin, lmax);
    

	// LOOP OVER REDSHIFTS
	double zmin=0, zmax=1.5;
	int nz=100;
	double dz=(zmax-zmin)/(nz-1);
	
	vector<double> Cug, Ciz;
	/*TVector<r_8> CugEllRed0(nz),CizEllRed0(nz),CugEllRed1(nz),CizEllRed1(nz),CugEllRed2(nz),CizEllRed2(nz);
	TVector<r_8> CugSbcRed0(nz),CizSbcRed0(nz),CugSbcRed1(nz),CizSbcRed1(nz),CugSbcRed2(nz),CizSbcRed2(nz);
	TVector<r_8> CugScdRed0(nz),CizScdRed0(nz),CugScdRed1(nz),CizScdRed1(nz),CugScdRed2(nz),CizScdRed2(nz);
	TVector<r_8> CugIrrRed0(nz),CizIrrRed0(nz),CugIrrRed1(nz),CizIrrRed1(nz),CugIrrRed2(nz),CizIrrRed2(nz);
	TVector<r_8> CugSB3Red0(nz),CizSB3Red0(nz),CugSB3Red1(nz),CizSB3Red1(nz),CugSB3Red2(nz),CizSB3Red2(nz);
	TVector<r_8> CugSB2Red0(nz),CizSB2Red0(nz),CugSB2Red1(nz),CizSB2Red1(nz),CugSB2Red2(nz),CizSB2Red2(nz);*/

	for (int i=0; i<nz; i++) {
	
		double z = zmin + i*dz;
		
		for (int j=0; j<sedArray.size(); j++) {
		
		    double cug = photometryCalcs.CompColor(z,(*sedArray[j]),(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		    double ciz = photometryCalcs.CompColor(z,(*sedArray[j]),(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		    
		    Cug.push_back(cug);
		    Ciz.push_back(ciz);
            }
		/*
		// El
		CugEllRed0(i)=photometryCalcs.CompColor(z,(*sedArray[iEl]),(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizEllRed0(i)=photometryCalcs.CompColor(z,(*sedArray[iEl]),(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
        CugEllRed1(i)=photometryCalcs.CompColor(z,ElRed1,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizEllRed1(i)=photometryCalcs.CompColor(z,ElRed1,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		CugEllRed2(i)=photometryCalcs.CompColor(z,ElRed2,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizEllRed2(i)=photometryCalcs.CompColor(z,ElRed2,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		
		// Sbc
		CugSbcRed0(i)=photometryCalcs.CompColor(z,(*sedArray[iSbc]),(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizSbcRed0(i)=photometryCalcs.CompColor(z,(*sedArray[iSbc]),(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
        CugSbcRed1(i)=photometryCalcs.CompColor(z,SbcRed1,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizSbcRed1(i)=photometryCalcs.CompColor(z,SbcRed1,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		CugSbcRed2(i)=photometryCalcs.CompColor(z,SbcRed2,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizSbcRed2(i)=photometryCalcs.CompColor(z,SbcRed2,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		
		// Scd
		CugScdRed0(i)=photometryCalcs.CompColor(z,(*sedArray[iScd]),(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizScdRed0(i)=photometryCalcs.CompColor(z,(*sedArray[iScd]),(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
        CugScdRed1(i)=photometryCalcs.CompColor(z,ScdRed1,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizScdRed1(i)=photometryCalcs.CompColor(z,ScdRed1,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		CugScdRed2(i)=photometryCalcs.CompColor(z,ScdRed2,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizScdRed2(i)=photometryCalcs.CompColor(z,ScdRed2,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		
		// Irr
		CugIrrRed0(i)=photometryCalcs.CompColor(z,(*sedArray[iIrr]),(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizIrrRed0(i)=photometryCalcs.CompColor(z,(*sedArray[iIrr]),(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
        CugIrrRed1(i)=photometryCalcs.CompColor(z,IrrRed1,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizIrrRed1(i)=photometryCalcs.CompColor(z,IrrRed1,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		CugIrrRed2(i)=photometryCalcs.CompColor(z,IrrRed2,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizIrrRed2(i)=photometryCalcs.CompColor(z,IrrRed2,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		
		// SB3
		CugSB3Red0(i)=photometryCalcs.CompColor(z,(*sedArray[iSB3]),(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizSB3Red0(i)=photometryCalcs.CompColor(z,(*sedArray[iSB3]),(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
        CugSB3Red1(i)=photometryCalcs.CompColor(z,SB3Red1,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizSB3Red1(i)=photometryCalcs.CompColor(z,SB3Red1,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		CugSB3Red2(i)=photometryCalcs.CompColor(z,SB3Red2,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizSB3Red2(i)=photometryCalcs.CompColor(z,SB3Red2,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		
		// SB2
		CugSB2Red0(i)=photometryCalcs.CompColor(z,(*sedArray[iSB2]),(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizSB2Red0(i)=photometryCalcs.CompColor(z,(*sedArray[iSB2]),(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
        CugSB2Red1(i)=photometryCalcs.CompColor(z,SB2Red1,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizSB2Red1(i)=photometryCalcs.CompColor(z,SB2Red1,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));
		CugSB2Red2(i)=photometryCalcs.CompColor(z,SB2Red2,(*CFHTfilters[uCFHT]),(*CFHTfilters[gCFHT]));
		CizSB2Red2(i)=photometryCalcs.CompColor(z,SB2Red2,(*CFHTfilters[iCFHT]),(*CFHTfilters[zCFHT]));*/
		}
		
	// WRITE TO A FILE

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
			
			outp << z <<"  ";
			
			for (int j=0; j<sedArray.size(); j++)
			    outp << Cug[j] <<"  "<< Ciz[j] <<"  ";
			/*// Ell
			outp <<"    "<<CugEllRed0(i)<<"    "<<CizEllRed0(i);
			outp <<"    "<<CugEllRed1(i)<<"    "<<CizEllRed1(i);
			outp <<"    "<<CugEllRed2(i)<<"    "<<CizEllRed2(i);
			// Sbc
			outp <<"    "<<CugSbcRed0(i)<<"    "<<CizSbcRed0(i);
			outp <<"    "<<CugSbcRed1(i)<<"    "<<CizSbcRed1(i);
			outp <<"    "<<CugSbcRed2(i)<<"    "<<CizSbcRed2(i);
			// Scd
			outp <<"    "<<CugScdRed0(i)<<"    "<<CizScdRed0(i);
			outp <<"    "<<CugScdRed1(i)<<"    "<<CizScdRed1(i);
			outp <<"    "<<CugScdRed2(i)<<"    "<<CizScdRed2(i);
			// Irr
			outp <<"    "<<CugIrrRed0(i)<<"    "<<CizIrrRed0(i);
			outp <<"    "<<CugIrrRed1(i)<<"    "<<CizIrrRed1(i);
			outp <<"    "<<CugIrrRed2(i)<<"    "<<CizIrrRed2(i);
			// SB3
			outp <<"    "<<CugSB3Red0(i)<<"    "<<CizSB3Red0(i);
			outp <<"    "<<CugSB3Red1(i)<<"    "<<CizSB3Red1(i);
			outp <<"    "<<CugSB3Red2(i)<<"    "<<CizSB3Red2(i);
			// SB2
			outp <<"    "<<CugSB2Red0(i)<<"    "<<CizSB2Red0(i);
			outp <<"    "<<CugSB2Red1(i)<<"    "<<CizSB2Red1(i);
			outp <<"    "<<CugSB2Red2(i)<<"    "<<CizSB2Red2(i);*/
			
			outp << endl;
		  }
		  outp.close();
	  }
	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;
	  
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " cfhtColors.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " cfhtColors.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " cfhtColors.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of cfhtColors.cc program  Rc= " << rc << endl;
  return rc;	
}
