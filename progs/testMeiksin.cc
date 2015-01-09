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
#define PI 3.141592
/*



*/
void usage(void);
void usage(void) {
	cout << endl<<" Usage: testMeiksin [...options...]" << endl<<endl;
	cout << " -o: OUTFILEROOT: write files to filename beginning OUTROOT (saved to testfiles/)"<<endl;
		
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

	cout << " ==== testMeiksin.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfileroot = "testfiles/testMeiksin";
    double zSource = 3.5;
    int nLines = 2000;
    int nl = 10000;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ho:z:n:w:")) != -1) {
	    switch (c) {
	        case 'o' :
	            outfileroot = optarg;
	            outfileroot = "testfiles/" + outfileroot;
	            break;
	  	    case 'z' :
	            sscanf(optarg,"%lf",&zSource);
	            break;
	        case 'n' :
	            sscanf(optarg,"%d",&nLines);
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
    cout <<"     Number of lines of sight = "<< nLines << endl;
    cout <<"     Wavelength resolution = "<< nl << endl;
    cout << endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    ifstream inp;
	ofstream outp;
	string outfile;
	
	LAFMeiksin lafMeiksin;
	
	// observed frame wavelengths in meters
    double lmin = 3.5e-7, lmax = 5.5e-7;
	double dl = (lmax - lmin)/(nl - 1);
	
	/*// open file
    outfile = outfileroot + "_lySeriesTransMeiksin.txt";
    cout << "     Writing to "<< outfile.c_str() << endl;
	outp.open(outfile.c_str(), ofstream::out);
	
	// loop over wavelengths
    for (int i=0; i<nl; i++) {
        //cout <<"     On wavelength "<<i+1<<" of "<<nl<<endl<<endl;
        double lam = lmin + i*dl;
            
        double trans = lafMeiksin.returnTransLymanSeries(lam, zSource);
        //double trans2 = lafMeiksin.taulya(lam*1e10,zSource); cout << endl;
        outp << lam <<"  "<< lam/(1+zSource) <<"  "<< trans <<endl;//<<"  "<< trans2 <<endl;
        }
    outp.close();*/
	
	
	// this controls the drawing of random numbers
	RandomGenerator rg;
	SimulateLLS simulateLLS(rg);
	
	/*// quick test of number of LLS simulated to redshift z
	for (int i=0; i<100; i++) {
	    double z = 2. + i*0.05;
	    cout << "     To z = "<< z <<" num LLS = "<< simulateLLS.nLineofSight(z)<< endl;
	    }*/
	    
	for (int l=0; l<nLines; l++) {
	
	    stringstream ss;
	    ss << l;
	    
	    // simulate a line of sight
	    vector<double> zAbsorbers, opticalDepth;
	    long nAbs = simulateLLS.simulateLineOfSight(zSource, zAbsorbers, opticalDepth);
	
	    cout <<"     Simulated line of sight "<< l <<", has "<< nAbs <<" absorbers"<<endl;
	    // open file
        outfile = outfileroot + "_lineOfSight"+ss.str()+".txt";
        cout << "     Writing to "<< outfile.c_str() << endl;
	    outp.open(outfile.c_str(), ofstream::out);
	    for ( int i=0; i<nAbs; i++) {
	        outp << "     "<< zAbsorbers[i] <<"  "<< opticalDepth[i]<<endl;
	        }
	    outp.close();
	        
	    MonteCarloMeiksin mcMeiksin(zAbsorbers, opticalDepth);
	
	    // open file
        outfile = outfileroot + "_fullTransMeiksin_lineOfSight"+ss.str()+".txt";
        cout << "     Writing to "<< outfile.c_str() << endl;
	    outp.open(outfile.c_str(), ofstream::out);
	    // loop over wavelengths
        for (int i=0; i<nl; i++) {
            //cout <<"     On wavelength "<<i+1<<" of "<<nl<<endl<<endl;
            double lam = lmin + i*dl;
                
            double trans = mcMeiksin.returnTransmission(lam,zSource);
            double transLS = exp(-mcMeiksin.returnLySeriesOnly(lam,zSource));
            double transDiff = exp(-mcMeiksin.returnDiffuseOnly(lam,zSource));
            double transLLS = exp(-mcMeiksin.returnLLSOnly(lam,zSource));

            outp << lam <<"  "<< lam/(1+zSource) <<"  "<< trans <<"  "<< transLS;
            outp <<"  "<< transDiff <<"  "<< transLLS <<endl;
            }
        outp.close();
	    }
	
	
        
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testMeiksin.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testMeiksin.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " testMeiksin.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of testMeiksin.cc program  Rc= " << rc << endl;
  return rc;	
}

