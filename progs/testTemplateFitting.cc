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
void usage(void)
{
	cout << endl<<" Usage: test [...options...]" << endl<<endl;
		
	cout << endl;
}
int main(int narg, char* arg[])
{
	cout << " ==== test.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfile = "testfiles/TESTCODE.txt";
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ho:")) != -1) 
	{
	switch (c) 
		{
	  case 'o' :
	    outfile = optarg;
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
	ofstream outp2;
	
	// this controls the drawing of random numbers
	RandomGenerator rg;

	// Set cosmology
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout <<"     Set cosmology to: OmegaM="<<OmegaM<<", OmegaL="<<OmegaL;
	cout <<", H0="<<100*h<<endl;
	cout << endl;
	
	// wavelength range of the SEDs/filters
	double lmin=5e-8, lmax=2.5e-6;
	
	// Load in SEDs
	string sedFile = "CWWK.list";
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax); // Read out SEDs into array
    int nsed=readSedList.getNSed();
    cout <<"     Number of original SEDs = "<<nsed<<endl;
    cout << endl;
    
    // Interpolate SEDs
	int nInterp = 2; 
	readSedList.interpSeds(nInterp);
    cout <<"     Interpolated SEDs "<<nInterp<<" times "<<endl;
    int ntot = readSedList.getNTot();
    cout <<"     Total number of SEDs now = "<<ntot<<endl;
	cout << endl;
	
    // Write un-reordered spectra to a file
	//outfile = "testfiles/spectraUnReordered.txt";
	//readSedList.writeSpectra(outfile);
	
	// Reorder SEDs and return them
    readSedList.reorderSEDs();
    vector<SED*> sedArray=readSedList.getSedArray();
    int nElliptical = 3;
    int nSpiral = 6;
	
	// Write reordered spectra to a file
	//outfile = "testfiles/spectraReordered.txt";
	//readSedList.writeSpectra(outfile);
	
	
	// GOODS filters
	cout <<"     GOODS filters "<<endl;
	string goodsFilterFile = "GOODS.filters";
	ReadFilterList readGOODSfilters(goodsFilterFile);
	readGOODSfilters.readFilters(lmin,lmax);
	vector<Filter*> GOODSfilters=readGOODSfilters.getFilterArray();
	int nGOODS=readGOODSfilters.getNTot();
	cout <<"     "<<nGOODS<<" GOODS filters read in "<<endl;
	int iBgoods = 1;

	
	sa_size_t ng = 10;
	double zmin = 0.01; 
	double zmax = 1; 
	double dz = (zmax - zmin)/(ng-1);
	
	
    // LSST filters
    string filterFile = "LSST.filters";
	ReadFilterList readLSSTfilters(filterFile);
	readLSSTfilters.readFilters(lmin,lmax);
	vector<Filter*> LSSTfilters=readLSSTfilters.getFilterArray();
    SimData simgal1(sedArray,LSSTfilters,su,rg,nElliptical,nSpiral);

    double lamEffu = simgal1.effectiveFilterWavelength((*LSSTfilters[0]));
    double lamEffg = simgal1.effectiveFilterWavelength((*LSSTfilters[1]));
    double lamEffr = simgal1.effectiveFilterWavelength((*LSSTfilters[2]));
    double lamEffi = simgal1.effectiveFilterWavelength((*LSSTfilters[3]));
    double lamEffz = simgal1.effectiveFilterWavelength((*LSSTfilters[4]));
    double lamEffy = simgal1.effectiveFilterWavelength((*LSSTfilters[5]));
    cout <<"     LSST filter effective wavelengths: "<<endl;
    cout <<"     u = "<< lamEffu <<endl;
	cout <<"     g = "<< lamEffg <<endl;
	cout <<"     r = "<< lamEffr <<endl;
	cout <<"     i = "<< lamEffi <<endl;
    cout <<"     z = "<< lamEffz <<endl;
    cout <<"     y = "<< lamEffy <<endl;
    
    
     // filter zero points
    double zpu = simgal1.getFilterZeroPointFlux((*LSSTfilters[0]));
    double zpg = simgal1.getFilterZeroPointFlux((*LSSTfilters[1]));
    double zpr = simgal1.getFilterZeroPointFlux((*LSSTfilters[2]));
    double zpi = simgal1.getFilterZeroPointFlux((*LSSTfilters[3]));
    double zpz = simgal1.getFilterZeroPointFlux((*LSSTfilters[4]));
    double zpy = simgal1.getFilterZeroPointFlux((*LSSTfilters[5]));
    cout <<"     LSST Filter zeropoints: "<< zpu <<"  "<< zpg <<"  "<< zpr;
    cout <<"  "<< zpi <<"  "<< zpz <<"  "<< zpy <<endl;
    cout << endl;
    

    TemplateChiSquare templateChiSq(sedArray,LSSTfilters,su);
    double aMin=0.001, aMax=0.1;
    int nA=100;
    templateChiSq.setAGrid(aMin, aMax, nA);
    
    TArray<double> chiSquare;
    //int ndim = 3;
    int ndim = 2;
    sa_size_t mydim[ndim];
    //mydim[0]=ng; mydim[1]=ntot; mydim[2] = nA;
    mydim[0]=ntot; mydim[1] = nA;
    chiSquare.SetSize(ndim, mydim);
    chiSquare = 0.;
   
    
    
    outfile = "testfiles/galFluxes.txt";
	outp.open(outfile.c_str(), ofstream::out);
	
	for (int i=0; i<ng; i++) {
	
	    cout <<"     Galaxy "<<i+1<<" of "<<ng<<endl;

        double zs = zmin + i*dz;
        //double type = 1.000;//3.012;
        double ext = 0.;
        double am = -20.;
        su.SetEmissionRedShift(zs);
        double dL = su.LuminosityDistanceMpc();
        
        // Simulate SED randomly
        int gtype = round( 0.5 + (3.5-0.5)*rg.Flat01());
        double type = simgal1.SimSED(gtype);
        int typ = simgal1.returnSedId(type);
        cout <<"     Simulated galaxy type: "<<type<<" ("<<typ<<")"<<endl;
        
        
        // Calculate galaxy magnitude in observed filters 
        // Galaxy's absolute magnitude is defined in filter (*GOODSfilters[iB])
	    double uMagTh=simgal1.GetMag(zs,type,am,ext,0,(*GOODSfilters[iBgoods]));
	    double gMagTh=simgal1.GetMag(zs,type,am,ext,1,(*GOODSfilters[iBgoods]));
	    double rMagTh=simgal1.GetMag(zs,type,am,ext,2,(*GOODSfilters[iBgoods]));
	    double iMagTh=simgal1.GetMag(zs,type,am,ext,3,(*GOODSfilters[iBgoods]));
	    double zMagTh=simgal1.GetMag(zs,type,am,ext,4,(*GOODSfilters[iBgoods]));
	    double yMagTh=simgal1.GetMag(zs,type,am,ext,5,(*GOODSfilters[iBgoods]));
	
	    // The final observations
	    // The 1st element is the value of the observed magnitude
	    // The 2nd element is the magnitude error
	    double fluxError = 0.0001;
	    vector<double> uObservation = simgal1.addError(uMagTh,fluxError,0);
	    vector<double> gObservation = simgal1.addError(gMagTh,fluxError,1);
	    vector<double> rObservation = simgal1.addError(rMagTh,fluxError,2);
	    vector<double> iObservation = simgal1.addError(iMagTh,fluxError,3);
	    vector<double> zObservation = simgal1.addError(zMagTh,fluxError,4);
	    vector<double> yObservation = simgal1.addError(yMagTh,fluxError,5);
	    double uMag = uObservation[0];
        double gMag = gObservation[0];
        double rMag = rObservation[0];
        double iMag = iObservation[0];
        double zMag = zObservation[0];
        double yMag = yObservation[0];
        double euMag = uObservation[1];
        double egMag = gObservation[1];
        double erMag = rObservation[1];
        double eiMag = iObservation[1];
        double ezMag = zObservation[1];
        double eyMag = yObservation[1];
        
        vector<double> obs;
        obs.push_back(uMag);
        obs.push_back(gMag);
        obs.push_back(rMag);
        obs.push_back(iMag);
        obs.push_back(zMag);
        obs.push_back(yMag);
        vector<double> errors;
        errors.push_back(euMag);
        errors.push_back(egMag);
        errors.push_back(erMag);
        errors.push_back(eiMag);
        errors.push_back(ezMag);
        errors.push_back(eyMag);
        
        int sedBestFit;
        double normBestFit;
        TArray<double> tmp=templateChiSq.galaxyChiSquared(obs, errors, zs, sedBestFit, normBestFit);
        stringstream ss; 
        ss << i;
        outfile = "testfiles/chkchisq"+ss.str()+".txt";
	    outp2.open(outfile.c_str(), ofstream::out);
	    for (int j=0; j<tmp.SizeX(); j++){
	        for (int k=0; k<tmp.SizeY(); k++)
	            outp2 << tmp(j,k) <<"  ";
	        outp2 << endl;
	        }
        outp2.close();
        
        // Look specifically at observed and theoretical rest-frame fluxes
        vector<double> lambdaRFs;
        lambdaRFs = simgal1.returnFilterRFWavelengths(zs);
        
        // theoretical rest-frame fluxes (these make sense!)
        //cout <<"    Calculating theoretical rest-frame fluxes"<<endl;
        double fu = simgal1.restFrameFluxLambda((*sedArray[typ]),(*LSSTfilters[0]),zs); 
        double fg = simgal1.restFrameFluxLambda((*sedArray[typ]),(*LSSTfilters[1]),zs);
        double fr = simgal1.restFrameFluxLambda((*sedArray[typ]),(*LSSTfilters[2]),zs);
        double fi = simgal1.restFrameFluxLambda((*sedArray[typ]),(*LSSTfilters[3]),zs);
        double fz = simgal1.restFrameFluxLambda((*sedArray[typ]),(*LSSTfilters[4]),zs);
        cout << endl;
        
        // observed observed-frame fluxes: (these make sense!)
        //cout <<"    Calculating observed observed-frame fluxes"<<endl;
        
        double fuo = simgal1.convertABMagToFluxLambda(uMag,zs,dL,(*sedArray[typ]),(*LSSTfilters[2]),(*LSSTfilters[0])); 
        double fgo = simgal1.convertABMagToFluxLambda(gMag,zs,dL,(*sedArray[typ]),(*LSSTfilters[2]),(*LSSTfilters[1]));
        double fro = simgal1.convertABMagToFluxLambda(rMag,zs,dL,(*sedArray[typ]),(*LSSTfilters[2]),(*LSSTfilters[2]));
        double fio = simgal1.convertABMagToFluxLambda(iMag,zs,dL,(*sedArray[typ]),(*LSSTfilters[2]),(*LSSTfilters[3]));
        double fzo = simgal1.convertABMagToFluxLambda(zMag,zs,dL,(*sedArray[typ]),(*LSSTfilters[2]),(*LSSTfilters[4]));
        cout << endl;
        
        
        outp << zs << "  " << typ << "  " << sedBestFit << "  ";
        outp << lambdaRFs[0] <<"  "<< fu <<"  "<< fuo <<"  ";
        outp << lambdaRFs[1] <<"  "<< fg <<"  "<< fgo <<"  ";
        outp << lambdaRFs[2] <<"  "<< fr <<"  "<< fro <<"  ";
        outp << lambdaRFs[3] <<"  "<< fi <<"  "<< fio <<"  ";
        outp << lambdaRFs[4] <<"  "<< fz <<"  "<< fzo <<"  ";
        outp << normBestFit << endl;
   
        }
    outp.close();        
		
	

	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " test.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " test.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " test.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of test.cc program  Rc= " << rc << endl;
  return rc;	
}

