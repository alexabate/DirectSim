// -*- LSST-C++ -*-
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
	cout << endl<<" Usage: simulateAbsorberLinesOfSight [...options...]" << endl<<endl;
	cout << " -o: OUTFILEROOT: write files to filename beginning OUTROOT (saved to output/)"<<endl;
	cout << " -n: NLINES: simulate [NLINES] lines of sight "<<endl;
	cout << " -z: ZSOURCE: redshift to simulate line of sight to"<<endl;
	cout << " -w: RES: wavelength resolution "<<endl;
	cout << endl;
};

int main(int narg, char* arg[]) {

	cout << " ==== simulateAbsorberLinesOfSight.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;
	
	string outfileroot = "output/simulateAbsorberLinesOfSight";
    int nLines = 1000;
    double zSource = 3.5;
    int nl = 1000;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
    while((c = getopt(narg,arg,"ho:n:z:w:")) != -1) {
	    switch (c)  {
	        case 'o' :
	            outfileroot = optarg;
	            outfileroot = "output/" + outfileroot;
	            break;
	        case 'n' :
	            sscanf(optarg,"%d",&nLines);
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
    cout <<"     Simulating "<< nLines <<" lines of sight out to z = "<< zSource <<endl;
    cout <<"     Wavelength resolution = "<< nl << endl;
    cout <<"     Writing to files beginning "<< outfileroot <<endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    RandomGenerator rg;
  
    string outfile;
    ifstream inp;
	ofstream outp;
	
	// observed frame wavelengths in meters
    double lmin = 3.5e-7, lmax = 5.5e-7;
	double dl = (lmax - lmin)/(nl - 1);
	
	SimulateLLS simulateLLS(rg);
		    
	for (int l=0; l<nLines; l++) {
	
	    stringstream ss,ss2;
	    ss << l+1;
	    ss2 << zSource;
	    
	    // simulate a line of sight
	    vector<double> zAbsorbers, opticalDepth;
	    long nAbs = simulateLLS.simulateLineOfSight(zSource, zAbsorbers, opticalDepth);
	
	    cout <<"     Simulated line of sight "<< l+1 <<", has "<< nAbs <<" absorbers"<<endl;
	    // open file
        outfile = outfileroot + "_lineOfSight"+ss.str()+"_toRedshift" + ss2.str() +".txt";
        cout << "     Writing to "<< outfile.c_str() << endl;
	    outp.open(outfile.c_str(), ofstream::out);
	    outp <<"# number of absorbers = "<< nAbs << endl;
	    for ( int i=0; i<nAbs; i++) {
	        outp << "     "<< zAbsorbers[i] <<"  "<< opticalDepth[i]<<endl;
	        }
	    outp.close();
	    
	    MonteCarloMeiksin mcMeiksin(zAbsorbers, opticalDepth);
	
	    // open file
        outfile = outfileroot + "_fullTransMeiksin_lineOfSight"+ss.str()+"_toRedshift" + ss2.str() +".txt";
        cout << "     Writing to "<< outfile.c_str() << endl;
	    outp.open(outfile.c_str(), ofstream::out);
	    outp <<"# number of absorbers = "<< nAbs << endl;
	    // loop over wavelengths
        for (int i=0; i<nl; i++) {
            //cout <<"     On wavelength "<<i+1<<" of "<<nl<<endl<<endl;
            double lam = lmin + i*dl;
                
            double trans = mcMeiksin.returnTransmission(lam,zSource);

            outp << lam <<"  "<<"  "<< trans <<endl;
            }
        outp.close();

	    }
	
    }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " simulateAbsorberLinesOfSight.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " simulateAbsorberLinesOfSight.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " simulateAbsorberLinesOfSight.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of simulateAbsorberLinesOfSight.cc program  Rc= " << rc << endl;
  return rc;	
}

