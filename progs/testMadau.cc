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


#include "genericfunc.h"
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
	cout << endl<<" Usage: testMadau [...options...]" << endl<<endl;
		
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

	cout << " ==== testMadau.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfileroot = "testfiles/testMadau";
    // source redshifts
    double zSource1 = 2.4;
	double zSource2 = 3.5;
	double zSource3 = 4.5;
  
	//--- decoding command line arguments 
	char c;
	while((c = getopt(narg,arg,"ho:z:")) != -1) {
	    switch (c) {
	        case 'o' :
	            outfileroot = optarg;
	            break;
	        case 'z' :
	            sscanf(optarg,"%lf,%lf,%lf",&zSource1,&zSource2,&zSource3);
	            break;
	        case 'h' :
		        default :
		        usage(); return -1;
		    }
	    }

  //-- end command line arguments
  cout << endl;
  
  int rc = 1;  
    try {  // exception handling try bloc at top level
  
    ifstream inp;
	ofstream outp;
	string outfile;
	
	// Observed wavelength ranges
	double lmin = 2500e-10, lmax = 8500e-10;
	int nl = 1000;
    double dl = (lmax - lmin)/(nl-1);
    
    // Class that calculates Madau absorption
    Madau madau;  
    
    
	
	// Write transmission as a function of observed wavelength to file
	outfile = outfileroot + "_madauTrans.txt";
    outp.open(outfile.c_str());
	cout << "     Writing transmission as a function of observed wavelength for "<<endl;
    cout << "     sources at z = "<< zSource1 <<", "<< zSource2 <<", "<< zSource3;
    cout << " to file "<< outfile << endl; 
	for (int i=0; i<nl; i++) {
	    double lambda = lmin + i*dl;
	    double trans1 = madau.returnObserverFrameTransmission(lambda, zSource1);
	    double trans2 = madau.returnObserverFrameTransmission(lambda, zSource2);
	    double trans3 = madau.returnObserverFrameTransmission(lambda, zSource3);
	    
	    outp << lambda <<"  "<< lambda/(1+zSource1) <<"  "<< lambda/(1+zSource2);
	    outp <<"  "<< lambda/(1+zSource3);
	    outp <<"  "<< trans1 <<"  "<< trans2 <<"  "<< trans3 << endl;
	    }
	outp.close();
	cout << endl;
	
	// Test effect of IGM transmission on a spectrum
	cout << "     Test the effect of IGM on a spectrum at a certain redshift"<<endl;
	
	// Load in SEDs
	lmin=5e-8, lmax=1.1e-6;
	dl = (lmax - lmin)/(nl-1);
	string sedFile = "CWWK.list";
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax); // Read out SEDs into array
    vector<SED*> sedArray=readSedList.getSedArray();
        
    // SB3 galaxy
    SED sed(*(sedArray[4]));
        
    // add Madau absorption
    SEDMadau sedMadauz0(sed,0.);
    SEDMadau sedMadauz1(sed,0.5);
    SEDMadau sedMadauz2(sed,1.0);
    SEDMadau sedMadauz3(sed,1.5);
    SEDMadau sedMadauz4(sed,2.0);
    SEDMadau sedMadauz5(sed,2.5);
	SEDMadau sedMadauz6(sed,3.0);
	SEDMadau sedMadauz7(sed,3.5);
	SEDMadau sedMadauz8(sed,4.0);
	
	outfile = outfileroot + "_SB3.txt";
	outp.open(outfile.c_str());
	cout << "     Writing SB3 flux after IGM absorption for different SB3 redshifts "<<endl;
	for (int i=0; i<nl; i++) {
	
	    double lambda = lmin + i*dl;
	    
	    double f0 = sedMadauz0(lambda);
	    double f1 = sedMadauz1(lambda);
        double f2 = sedMadauz2(lambda);
        double f3 = sedMadauz3(lambda);
        double f4 = sedMadauz4(lambda);
        double f5 = sedMadauz5(lambda);
        double f6 = sedMadauz6(lambda);
        double f7 = sedMadauz7(lambda);
        double f8 = sedMadauz8(lambda);
        
        outp << lambda <<"  "<< f0 <<"  "<< f1 <<"  "<< f2 <<"  "<< f3 <<"  ";
        outp << f4 <<"  "<< f5 <<"  "<< f6 <<"  "<< f7 <<"  "<< f8 <<endl;
        
        }
    outp.close();
    cout << endl;
    
    // Test effect of IGM on u and g mags
        
    // Load in filters required
	// wavelength range of the SEDs/filters
	cout <<"     Load in LSST filters"<<endl;
	//double lmin=5e-8, lmax=2.5e-6;
	string filterFile = "LSST.filters";
	ReadFilterList readFilterList(filterFile);
	readFilterList.readFilters(lmin,lmax);
	vector<Filter*> LSSTfilters=readFilterList.getFilterArray();
	
	cout <<"     Load in GOODS B filter"<<endl;
	string goodsFilterFile = "GOODSB.filters";
	ReadFilterList readGOODSBfilter(goodsFilterFile);
	readGOODSBfilter.readFilters(lmin,lmax);
	vector<Filter*> goodsBFilter=readGOODSBfilter.getFilterArray();
	
	// this controls the drawing of random numbers
	RandomGenerator rg;

	// Set cosmology
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout <<"     Set cosmology to: OmegaM="<<OmegaM<<", OmegaL="<<OmegaL;
	cout <<", H0="<<100*h<<endl;
	cout << endl;
	
	// Prepare the class which will calculate the magnitudes
	cout <<"     Initialize class to calculate magnitudes"<<endl;
	int nElliptical = 1;
    int nSpiral = 2;
	SimData simgal(sedArray,LSSTfilters,su,rg,nElliptical,nSpiral);
	SimData simgalNoIGM(sedArray,LSSTfilters,su,rg,nElliptical,nSpiral);
	simgalNoIGM.setMadau(false);
    cout << endl;
    
    
    double zmin = 1.5, zmax = 3.7;
    int nz = 200;
    double dz = (zmax-zmin)/(nz-1);
    
    outfile = outfileroot + "_IGMdiffs.txt";
	outp.open(outfile.c_str());
    for (int i=0; i<nz; i++){
    
        double zs = zmin + i*dz;
        double am = -18;
        double ext = 0.;
        double type = 3.004;

        // Calculate galaxy magnitude in observed filters ugrizy
        // Galaxy's absolute magnitude is defined in filter (*goodsBFilter[0])
		double uMag=simgal.GetMag(zs,type,am,ext,0,(*goodsBFilter[0]));
		double gMag=simgal.GetMag(zs,type,am,ext,1,(*goodsBFilter[0]));
		double uMagNoIGM=simgalNoIGM.GetMag(zs,type,am,ext,0,(*goodsBFilter[0]));
		double gMagNoIGM=simgalNoIGM.GetMag(zs,type,am,ext,1,(*goodsBFilter[0]));
		
		outp << zs <<"  "<< uMag <<"  "<< gMag <<"  "<< uMagNoIGM <<"  "<< gMagNoIGM <<endl;
        }
    outp.close();
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testMadau.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testMadau.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " testMadau.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of testMadau.cc program  Rc= " << rc << endl;
  return rc;	
}

