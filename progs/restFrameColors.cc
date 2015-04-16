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


int main(int narg, char* arg[])
{
  cout << " ==== restFrameColors.cc program , calculate LSST rest frame colors ====" << endl << endl;
  if (narg < 2) {
        cout << " Usage: restFrameColors infile outfile" << endl;
        cout << "        outfile: output filename"<<endl;
        cout << "        infile: input filename"<<endl;
        cout << endl;
        cout << endl;
        return 1;
        }

  // make sure SOPHYA modules are initialized 
  SophyaInit();  
  FitsIOServerInit();
  InitTim();
  cout<<endl<<endl;

  //--- decoding command line arguments 
  string infile = arg[1];
  string outfile = arg[2];
  cout <<"     Reading in catalog from "<< infile <<endl;
  cout <<"     Writing out restframe colors to "<< outfile <<endl;
  cout << endl;
  //-- end command line arguments
  
  int rc = 1;  
  try {  // exception handling try bloc at top level
	ResourceUsage res;
	InitTim();
	
	// wavelength range of the SEDs/filters
	double lmin=5e-8, lmax=2.5e-6;
	
	//FILTERS
	cout <<"    restFrameColors: Loading in filters ..."<<endl;
	string filterFile = "LSST.filters";
	ReadFilterList readFilters(filterFile);
	readFilters.readFilters(lmin,lmax);
	vector<Filter*> filters=readFilters.getFilterArray();
	int nFilter=readFilters.getNTot();
	cout <<"     "<< nFilter <<" CFHT filters read in "<<endl;
    cout << endl;
	
    // filter numbers
	int iU = 0;
	int iG = 1;
	int iR = 2;
	int iI = 3;
	int iZ = 4;
	int iY = 5;
	
		
	// GALAXY SED TEMPLATES
	string sedFile = "CWWK.list";
	ReadSedList readSedList(sedFile);
	// Read out SEDs into array
    readSedList.readSeds(lmin,lmax);
    int nInterp = 9;
    readSedList.interpSeds(nInterp);
    readSedList.reorderSEDs();
    int ntot = readSedList.getNTot();
    vector<SED*> sedArray=readSedList.getSedArray();
    cout <<"     Final total number of SEDs = "<< ntot <<endl;
	cout << endl;
    
	
    // PHOTOMETRY CALCULATION CLASS
    PhotometryCalcs photometryCalcs(lmin,lmax);
    

	// LOOP OVER FILE
    ifstream inp;
    ofstream outp;
    
    inp.open(infile.c_str(), ifstream::in);
    outp.open(outfile.c_str(), ofstream::out);
    
    double redshift = 0.001; // effectively z=0
    int cnt=0;
    string line;
    if (inp.is_open()) {
        cout << "     File has been sucessfully opened"<<endl;
        while ( getline(inp, line) ) {
            stringstream ss;
            ss.str(line);
            double zs;
            double zp;
            int ty;
            double mu,mg,mr,mi,mz,my;
            double am;
            double ebv;
            ss >> zs >> zp >> ty >> mu >> mg >> mr >> mi >> mz >> my >> am >> ebv;
            
            cnt++;
            cout <<"     On galaxy "<< cnt << endl;
            
            // redden SED
            int law = 0; // Cardelli law  by default
            if (ty>=25)
                law = 1; // Calzetti law for more Starburst types
                
                
            sedArray[ty]->doRedden(ebv, law);
            //SEDRedden sedred((*sedArray[ty]),ebv,law);

            double cug, cgr, cri, ciz, czy;
            cug = photometryCalcs.CompColor(redshift, (*sedArray[ty]), (*filters[iU]),(*filters[iG]));
            cgr = photometryCalcs.CompColor(redshift, (*sedArray[ty]), (*filters[iG]),(*filters[iR]));
            cri = photometryCalcs.CompColor(redshift, (*sedArray[ty]), (*filters[iR]),(*filters[iI]));
            ciz = photometryCalcs.CompColor(redshift, (*sedArray[ty]), (*filters[iI]),(*filters[iZ]));
            czy = photometryCalcs.CompColor(redshift, (*sedArray[ty]), (*filters[iZ]),(*filters[iY]));
            
            outp << zs <<"  "<< ty <<"  "<< am << "  "<< mi <<"  ";
            outp << cug <<"  "<< cgr <<"  "<< cri <<"  "<< ciz <<"  "<< czy << endl; 
            }
        }
    else
        cout <<"     File failed to open"<<endl;
	
    outp.close();
    inp.close();
	  
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " restFrameColors.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " restFrameColors.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " restFrameColors.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of restFrameColors.cc program  Rc= " << rc << endl;
  return rc;	
}
