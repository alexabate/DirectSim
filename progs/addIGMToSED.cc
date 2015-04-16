#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

// sophya libraries
//#include "histinit.h"
#include "fiosinit.h"
#include "mydefrg.h"


//#include "genericfunc.h"
#include "constcosmo.h"
#include "cosmocalcs.h"
#include "igm.h"
#include "mydefrg.h"
#include "sedfilter.h"
#include "simdata.h"


// root
//#include "root_plots.h"
//#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TPrincipal.h"
#include <cmath>
//#include "iomanip.h"
#include "TRandom.h"

void usage(void);
void usage(void) {
	cout << endl<<" Usage: addIGMToSED [...options...]" << endl <<endl;
	
    cout << " -o: OUTFILEROOT: write files to filename beginning OUTROOT (saved to output/)"<<endl;
    cout << "                  [DEFAULT=addIGMToSED]"<<endl;
    cout << " -i: INFILEROOT: read line of sight transmissions from filename beginning INFILEROOT"<<endl;
	cout << " -n: NLINES: [NLINES] line of sight transmissions to read in [DEFAULT=20]"<<endl;
	cout << " -z: ZSOURCE: redshift line of sight was simulated to AND redshift of "<<endl;
	cout << "              galaxy SED [DEFAULT=3.5]"<<endl;
	cout << " -s: SEDNUM: sed number in library to add IGM absorption to [DEFAULT=7]"<<endl;
	cout << endl;
	
	cout << " Reads in [NLINES] line of sight flux fraction transmission per "<< endl;
	cout << " observed wavelength and applies it to the SED [SEDNUM] in the library"<<endl;
	cout << " taken to be at redshift [ZSOURCE] "<< endl;
	cout << endl;
	
	cout << " Outputs to files beginning [OUTFILEROOT] for each line of sight "<<endl;
	cout << " a file with the following columns:"<< endl;
    cout << " obsWaveLength, emissionWaveLength, fluxFractionTransmitted, "<<endl;
    cout << " sedFluxAtemission, sedWIGMFluxAtobsL, sedWMadauFluxAtobsL"<<endl;
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

	cout << " ==== addIGMToSED.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfileroot = "output/addIGMToSED";
    string infileroot;
    int nLines = 1000;
    double zSource = 3.5;
    int sedNo = 7;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ho:i:n:z:s:")) != -1) {
	    switch (c) {
            case 'o' :
                outfileroot = optarg;
                outfileroot = "output/" + outfileroot;
                break;
            case 'i' :
                infileroot = optarg;
                break;
            case 'n' :
                sscanf(optarg,"%d",&nLines);
                break;
            case 'z' :
                sscanf(optarg,"%lf",&zSource);
                break;
            case 's' :
                sscanf(optarg,"%d",&sedNo);
                break;
            case 'h' :
                default :
                usage(); return -1;
            }
	    }

    //-- end command line arguments
    cout <<"     Reading in "<< nLines <<" lines of sight out to z = "<< zSource <<endl;
    cout <<"     from files beginning: "<< infileroot << endl;
    cout <<"     Adding IGM absorption to SED number "<< sedNo << endl;
    cout <<"     Writing out SED data to files beginning "<< outfileroot <<endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    ifstream inp;
	ofstream outp;
	string infile, outfile;  
	
	// Load in SEDs
	double lmin=100e-10, lmax=12000e-10;
	string sedFile = "CWWKSB.list";
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax); // Read out SEDs into array
    vector<SED*> sedArray=readSedList.getSedArray();

    // Load in filters required
    // LSST
	cout <<"     Load in LSST filters"<<endl;
	string filterFile = "LSST.filters";
	ReadFilterList readFilterList(filterFile);
	readFilterList.readFilters(lmin,lmax);
	vector<Filter*> LSSTfilters=readFilterList.getFilterArray();
	// Reference
	cout <<"     Load in reference (GOODS B) filter"<<endl;
	string goodsFilterFile = "GOODSB.filters";
	ReadFilterList readGOODSBfilter(goodsFilterFile);
	readGOODSBfilter.readFilters(lmin,lmax);
	vector<Filter*> goodsBFilter=readGOODSBfilter.getFilterArray();
	
	//double lminObs = 3000e-10, lmaxObs = 1100e-9;
	// The below should really be set by the wavelength range the transmission
	// along the line of sight was calculated for
	double lminObs = 3.e-8, lmaxObs = 5.5e-7;
	int nl =10000;
    double dl = (lmaxObs-lminObs)/(nl-1);
    
    stringstream ssn;
	ssn << nLines;
    stringstream ssz;
    ssz << zSource;
    
    for (int i=0; i<nLines; i++){
    
        stringstream ss;
        ss << i+1;

        // read in file with IGM transmission
        infile = infileroot + ss.str() + "of" + ssn.str() + ".txt";
        string line;
        inp.open(infile.c_str());
        getline(inp,line);
        string nAbsString;
        stringstream num;
        for (int i=24; i<line.size(); i++)
            num << line[i];
        int nAbs = atoi(num.str().c_str());
        cout <<"     Number of absorbers = "<<nAbs <<endl;
        

        // fraction of flux transmitted per observed wavelength 
        IGMTransmission igmTransmission(infile);
        // redshift SED and apply line of sight IGM absorption
        SEDIGM sedIGM(*(sedArray[sedNo]), igmTransmission, zSource);
        // redshift SED and apply average Madau IGM absorption
        SEDIGM sedMadau(*(sedArray[sedNo]), zSource);
        
        outfile = outfileroot + "_zSource" + ssz.str() + "_SightLine";
	    outfile += ss.str() + "of" + ssn.str() + ".txt";
	    
        outp.open(outfile.c_str());
        for (int i=0; i<nl; i++) {
        
            double lamO = lminObs +i*dl;
            double lamE = lamO/(1. + zSource);
            
            // columns are:
            // obsWaveLength - emissionWaveLength - fluxFractionTransmitted - 
            // sedFluxAtemission - sedWIGMFluxAtobsL - sedWMadauFluxAtobsL
            
            outp << lamO <<"  "<< lamE <<"  "<< igmTransmission(lamO) <<"  ";
            outp << sedArray[sedNo]->returnFlux(lamO) <<"  "<< sedIGM(lamO);
            outp <<"  "<< sedMadau(lamO) << endl;
            }
        outp.close();

		
		inp.close();
        }
    cout << endl;
    
    
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " addIGMToSED.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " addIGMToSED.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " addIGMToSED.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of addIGMToSED.cc program  Rc= " << rc << endl;
  return rc;	
}

