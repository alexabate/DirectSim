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
	cout << endl<<" Usage: lineOfSightLymanAlpha [...options...]" << endl<<endl;
	cout << " -o: OUTFILEROOT: write files to filename beginning OUTROOT (saved to output/)"<<endl;
	cout << endl;
};

int main(int narg, char* arg[]) {

	cout << " ==== lineOfSightLymanAlpha.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;
	
	string outfileroot = "output/lineOfSightLymanAlpha";
    double zSource = 3.5;
    int nLines = 20;
    int nl = 10000;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
    while((c = getopt(narg,arg,"ho:z:n:w:")) != -1) {
	    switch (c)  {
	        case 'o' :
	            outfileroot = optarg;
	            outfileroot = "output/" + outfileroot;
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
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    string outfile;
    ifstream inp;
	ofstream outp;
	
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
    cout << endl;


	// Calculates f(z) given in eqn 5 in Inoue & Iwata 2008
	cout <<"     Calculating absorber redshift distribution"<<endl;
    double gamma1=0.2, gamma2=2.5, gamma3=4.0;
    AbsorberRedshiftDistribution absorberZDist(gamma1,gamma2,gamma3);
    double zmin=0.2, zmax=6;
    double dz=(zmax-zmin)/(nStep-1);
    cout << endl;
	
	
	// Calculates h(b) given in eqn 6 in Inoue & Iwata 2008
	cout <<"     Calculating absorber b-parameter distribution"<<endl;
    double bsig=23;
    DopplerParDistribution dopplerParDist(bsig);
    double bmin=2, bmax=200;
    double db=(bmax-bmin)/(nStep-1);
    cout << endl;
	
	
    // Generate Monte Carlo distributions
    cout <<"     Generate Monte Carlo distribution of absorbers "<<endl;
    RandomGenerator rg;
	ProbabilityDistAbsorbers probDistAbsorbers(rg, absorberZDist, 
	                                           hiColumnDensity, dopplerParDist);

    int nLya = 2; // lyman-alpha only
    
    // observed frame wavelengths in meters
    double lmin = 3.5e-7, lmax = 5.5e-7;
	double dl = (lmax - lmin)/(nl - 1);
	                                                
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



        // Lyman-alpha only
        LineOfSightTrans lightOfSightTransLya(redshifts, dopplers, columnDensities);
        lightOfSightTransLya.setnLineMax(nLya);
        lightOfSightTransLya.setLymanSeriesOnly();
        
        // Lyman-series only
        LineOfSightTrans lightOfSightTransSeries(redshifts, dopplers, columnDensities);
        lightOfSightTransSeries.setLymanSeriesOnly();


	    // open file
	    outfile = outfileroot + "_transLineofsight"+ss.str()+".txt";
        cout << "     Writing to "<< outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
		// loop over wavelengths
	    for (int i=0; i<nl; i++) {
	        cout <<"     On wavelength "<<i+1<<" of "<<nl<<endl;
	        double lam = lmin + i*dl;
	            
	        double transLa = lightOfSightTransLya(lam, zSource);
            double transLs = lightOfSightTransSeries(lam, zSource);

	        outp << lam <<"  "<< lam/(1+zSource) <<"  "<< transLa <<"  "<< transLs <<endl;
	        }
        outp.close();
        }
        
	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " lineOfSightLymanAlpha.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " lineOfSightLymanAlpha.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " lineOfSightLymanAlpha.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of lineOfSightLymanAlpha.cc program  Rc= " << rc << endl;
  return rc;	
}

