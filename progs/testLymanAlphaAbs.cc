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


#include "constcosmo.h"
#include "cosmocalcs.h"
#include "igm.h"
#include "genericfunc.h"
#include "poly.h"
#include "mydefrg.h"
#include "geneutils.h"


void usage(void);
void usage(void) {
	cout << endl<<" Usage: testLymanAlphaAbs [...options...]" << endl<<endl;
	cout << " -o: OUTFILEROOT: write files to filename beginning OUTROOT (saved to output/)"<<endl;
	cout << endl;
    };

int main(int narg, char* arg[]) {

    cout << " ==== testLymanAlphaAbs.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;
	
	string outfileroot = "output/testLymanAlphaAbs";
    double zSource = 3.5;
    int nl = 1000;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
    while((c = getopt(narg,arg,"ho:z:w:")) != -1) {
	    switch (c)  {
	        case 'o' :
	            outfileroot = optarg;
	            outfileroot = "output/" + outfileroot;
	            break;
	        case 'z' :
	            sscanf(optarg,"%lf",&zSource);
	            break;
	        case 'w' :
	            sscanf(optarg,"%d",&nl);
	            break;
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }

    //-- end command line arguments
    cout <<"     Writing to files beginning "<< outfileroot <<endl;
    cout <<"     Redshift of source = " << zSource << endl;
    cout <<"     Wavelength resolution = "<< nl << endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    string outfile;
    ifstream inp;
	ofstream outp;
	
	// Lyman-alpha only
	int nLine = 2;
	
	// wavelength range in m
	double lmin = 1214e-10, lmax = 1218e-10;
	double dl = (lmax - lmin)/(nl - 1);
	
	double bmin = 10, bmax = 100;
	int nb = 10;
	double db = (bmax - bmin)/(nb - 1);
	
	outfile = outfileroot + "_voigt.txt";
	outp.open(outfile.c_str(), ofstream::out);
	
	// first row is lambdas
	for (int j=0; j<nl; j++) {
        double lam = lmin + j*dl;
        outp << lam << "  ";
        }
	outp << endl;
	
	for (int i=0; i<nb; i++) {
        
        double b = bmin + i*db;
        VoigtProfile voigtProf(b, nLine);
        
        for (int j=0; j<nl; j++) {
            double lam = lmin + j*dl;
            double vp = voigtProf(lam);
            outp << vp << "  ";
            }
        outp << endl;
        }
    outp.close();
    
    OpticalDepth opticalDepth;
    opticalDepth.setLymanSeriesOnly();
    opticalDepth.setMaxLine(nLine);
    
    outfile = outfileroot + "_crosssect.txt";
	outp.open(outfile.c_str(), ofstream::out);
	
	// first row is lambdas
	for (int j=0; j<nl; j++) {
        double lam = lmin + j*dl;
        outp << lam << "  ";
        }
	outp << endl;
	
	for (int i=0; i<nb; i++) {
        
        double b = bmin + i*db;
        
        for (int j=0; j<nl; j++) {
            double lam = lmin + j*dl;
            double freq = SPEED_OF_LIGHT_MS/lam;
            double cs = opticalDepth.returnLymanSeriesCrossSection(freq, b);
            outp << cs << "  ";
            }
        outp << endl;
        }
    outp.close();
    
    // new wavelength range in m
	lmin = 1200e-10, lmax = 1300e-10;
	nl *= 10;
    double nHILAF = 1e15;
    double nHIDLA = 1e21;
    outfile = outfileroot + "_optdepth.txt";
	outp.open(outfile.c_str(), ofstream::out);
	
	// first row is lambdas
	for (int j=0; j<nl; j++) {
        double lam = lmin + j*dl;
        outp << lam << "  ";
        }
	outp << endl;
	
	for (int i=0; i<nb; i++) {
        
        double b = bmin + i*db;
        
        for (int j=0; j<nl; j++) {
            double lam = lmin + j*dl;
            double freq = SPEED_OF_LIGHT_MS/lam;
            double od = opticalDepth.returnRestFrameOpticalDepth(freq, b, nHILAF);
            outp << od << "  ";
            }
        outp << endl;
        
        for (int j=0; j<nl; j++) {
            double lam = lmin + j*dl;
            double freq = SPEED_OF_LIGHT_MS/lam;
            double od = opticalDepth.returnRestFrameOpticalDepth(freq, b, nHIDLA);
            outp << od << "  ";
            }
        outp << endl;
        }
    outp.close();
    
	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testLymanAlphaAbs.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testLymanAlphaAbs.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " testLymanAlphaAbs.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of testLymanAlphaAbs.cc program  Rc= " << rc << endl;
  return rc;	
}

