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
	cout << " -o: OUTFILEROOT: write files to filename beginning OUTROOT"<<endl;
	cout << " -n: NLINES: simulate [NLINES] lines of sight "<<endl;
	cout << endl;
};

int main(int narg, char* arg[]) {

	cout << " ==== simulateAbsorberLinesOfSight.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;
	
	string outfileroot = "testfiles/simulateAbsorberLinesOfSight";
    int nLines = 1000;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
    while((c = getopt(narg,arg,"ho:n:")) != -1) {
	    switch (c)  {
	        case 'o' :
	            outfileroot = optarg;
	            break;
	        case 'n' :
	            sscanf(optarg,"%d",&nLines);
		        break;
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }

    //-- end command line arguments
    cout <<"     Simulating "<< nLines <<" lines of sight "<<endl;
    cout <<"     Writing to files beginning "<< outfileroot <<endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    string outfile;
    ifstream inp;
	ofstream outp;
	
	// Write HI column density distribution function to a file.
	// Calculates g(N_HI) given in eqn 4 Inoue & Iwata 2008
	cout <<"     Calculating absorber HI column density distribution"<<endl;
    double beta1=1.6, beta2=1.3;
    int nStep = 10000;
    HIColumnDensity hiColumnDensity(beta1,beta2,1.6e17,1e12,1e22, nStep);
    double Nl, Nu;
    hiColumnDensity.returnLowerUpperColDensValues(Nl,Nu);
    //double intVal = hiColumnDensity.checkIntegration(nStep);
    nStep=500000;
    double logNl=log(Nl);
    double dlogN=(log(Nu)-logNl)/(nStep-1);
    
    cout << "     dlogN = "<<dlogN << endl;
    outfile = outfileroot + "_gNHI.txt";
    hiColumnDensity.writeToFile(outfile,dlogN,nStep);
    cout << endl;

		
	// Write absorber redshift distribution to a file.
	// Calculates f(z) given in eqn 5 in Inoue & Iwata 2008
	cout <<"     Calculating absorber redshift distribution"<<endl;
    double gamma1=0.2, gamma2=2.5, gamma3=4.0;
    AbsorberRedshiftDistribution absorberZDist(gamma1,gamma2,gamma3);
    double zmin=0.2, zmax=6;
    double dz=(zmax-zmin)/(nStep-1);
    outfile= outfileroot + "_fz.txt";
    absorberZDist.writeToFile(outfile,zmin,dz,nStep);
    cout << endl;
	
	
	// Write doppler parameter distribution to a file.
	// Calculates h(b) given in eqn 6 in Inoue & Iwata 2008
	cout <<"     Calculating absorber b-parameter distribution"<<endl;
    double bsig=23;
    DopplerParDistribution dopplerParDist(bsig);
    double bmin=2, bmax=200;
    double db=(bmax-bmin)/(nStep-1);
    outfile= outfileroot + "_hb.txt";
    dopplerParDist.writeToFile(outfile,bmin,db,nStep);
    cout << endl;
	
	
    // Generate Monte Carlo distributions
    cout <<"     Generate Monte Carlo distribution of absorbers "<<endl;
    RandomGenerator rg;
	ProbabilityDistAbsorbers probDistAbsorbers(rg, absorberZDist, 
	                                           hiColumnDensity, dopplerParDist);
	double zStart = 0, zMax =6;                                       
	// Now simulate the lines of sight
	for (int i=0; i<nLines; i++) {
	    cout <<"     Simulating line of sight "<< i+1 <<" of "<< nLines <<endl;
	    
	    stringstream ss;
	    ss << i+1;
	    vector<double> redshifts, dopplers, columnDensities;
	    outfile = outfileroot + "_lineofsight"+ss.str()+".txt";
	    probDistAbsorbers.simulateLineOfSight(zStart,zMax,redshifts, dopplers, 
	                                                columnDensities, outfile);
	    cout << endl;
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

