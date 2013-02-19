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
#include "tarray.h"
#include "ctimer.h"

// my classes
#include "sinterp.h"
#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "simdata.h"

void usage(void);
void usage(void) {
	cout << endl<<" Usage: testKcorrMethod [...options...]" << endl<<endl;
	cout << " -o OUTFILEROOT: filename to write stuff out to "<<endl;
	cout << " -z ZRES: resolution of redshift grid to read in "<<endl;
	cout << " -e ERES: resolution of extinction grid to read in "<<endl;
	cout << " -c iDZT,iEZT: factor of step increase in redshift grid, extinction";
	cout << " grid to calculate with [DEFAULT=100,10] "<<endl;
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

    cout << " ==== testKcorrMethod.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfileroot = "testfiles/testKcorrMethod";
    int zres = 2000;
    int eres = 200;
    int idzt = 100;
    int idet = 10;
    
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ho:z:e:c:")) != -1) {
	    switch (c) {
            case 'o' :
                outfileroot = optarg;
                break;
            case 'z' :
                sscanf(optarg,"%d",&zres);
                break;
            case 'e' :
                sscanf(optarg,"%d",&eres);
                break;
            case 'c' :
                sscanf(optarg,"%d,%d",&idzt,&idet);
                break;
            case 'h' :
                default :
                usage(); return -1;
		    }
        }

    //-- end command line arguments
    cout <<"     Writing to file(s) beginning "<< outfileroot <<endl;
    cout <<"     Resolution of redshift grid to read in "<< zres <<endl;
    cout <<"     Resolution of extinction grid to read in "<< eres <<endl;
    cout <<"     Increase in step of redshift grid to calculate with "<< idzt <<endl;
    cout <<"     Increase in step of extinction grid to calculate with "<< idet <<endl;
    cout << endl;
  
  int rc = 1;  
    try {  // exception handling try bloc at top level
  
    ifstream inp;
	ofstream outp;
	string outfile;
	
	// For timing stuff
	Timer tm("timer",false);
	
	// this controls the drawing of random numbers
	RandomGenerator rg;
	
	// SED and filter files
	string sedFile = "CWWK.list";
	string filterFile = "LSST.filters";

	// Set cosmology
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout <<"     Set cosmology to: OmegaM="<<OmegaM<<", OmegaL="<<OmegaL;
	cout <<", H0="<<100*h<<endl;
	cout << endl;

	
	int nTrials = 10000;
	tm.Split();
	for (int i=0; i<10000; i++) {
	    double ztmp = 0.1+0.0003*i;
	    
	    su.SetEmissionRedShift(ztmp);
	    double mu=5*log10(su.LuminosityDistanceMpc())+25;
	    
	    }
	tm.Split();
    int timeMu = tm.PartialElapsedTimems();
	cout << "     Time per mu = "<< (double)timeMu/nTrials <<" ms"<< endl;


	// Load in filters required
	// wavelength range of the SEDs/filters
	cout <<"     Load in filters"<<endl;
	double lmin=5e-8, lmax=2.5e-6;
	ReadFilterList readFilterList(filterFile);
	readFilterList.readFilters(lmin,lmax);
	vector<Filter*> filters = readFilterList.getFilterArray();
	int nFilter = readFilterList.getNTot();
	cout <<"     Read in "<< nFilter <<" filters"<<endl;
	cout <<endl;
	
	// Load in the rest-frame filter
	cout <<"     Load in GOODS B filter (the rest-frame filter)"<<endl;
	string rfFilter = "GOODSB";
	string goodsFilterFile = rfFilter + ".filters";
	ReadFilterList readGOODSBfilter(goodsFilterFile);
	readGOODSBfilter.readFilters(lmin,lmax);
	vector<Filter*> goodsBFilter=readGOODSBfilter.getFilterArray();
	Filter restFrameFilter((*goodsBFilter[0]));
	
	// Load in SEDs
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax); // Read out SEDs into array  
    int nSED = readSedList.getNTot();
    vector<SED*> sedArray=readSedList.getSedArray();
    cout <<"     Final total number of SEDs = "<< nSED <<endl;
	cout << endl;
	
	// Madau absorption setting
	bool isMadau = false;
	
	// Get k-correction tables
	string sedLib="CWWK", filtSet = "LSST", restFrameFilt = "GOODSB";
	double zmin=0., zmax=3., emax=0.3;
	double dz = (zmax-zmin)/(zres-1), de = emax/(eres-1);
    ReadKCorrections readK(sedLib, filtSet, restFrameFilt, zmin, zmax, zres, emax, eres, isMadau);
    readK.readInterpZExt(nSED, nFilter);
    vector<SInterp2D*> kInterpZExt = readK.returnkInterpZExt();
    cout << " kInterpZExt.size()="<< kInterpZExt.size() << endl;
        
	// Prepare the class which will calculate the magnitudes
	cout <<"     Initialize classes to calculate magnitudes"<<endl;
	int nElliptical = 1;
    int nSpiral = 2;
	SimData simdataCalc(sedArray, filters, su, rg, nElliptical, nSpiral);
	simdataCalc.setMadau(isMadau);
	SimData simdataInterp(su, rg, kInterpZExt, nFilter, nElliptical, nSpiral);
	simdataInterp.setMadau(isMadau);
            

    // count different k corrections
    int cnt1=0;
    int cnt2=0;
    
    // statistics of worst interpolations
    double sumSED=0, sumSEDsq=0;
	double sumFilt=0, sumFiltsq=0;
	double sumz=0, sumzsq=0;
	double sume=0, sumesq=0;
	
	double ztmin = zmin + dz/2.;
	double etmin = de/2.;
	double ztmax = zmax + dz/2.;
	double etmax = emax + de/2.;
	double dzt = idzt*dz, det = idet*de;
	int nz = floor((ztmax-ztmin)/dzt);
	int ne = floor((etmax-etmin)/det);
	cout <<"     Calculating with redshift grid from "<< ztmin <<" in "<< nz;
	cout <<" steps of "<< dzt << endl;
	cout <<"     Calculating with extinction grid from "<< etmin <<" in "<< ne;
	cout <<" steps of "<< det << endl;
	cout << endl;
	
	int calcTime=0;
	int interpTime=0;
	
	string file = outfileroot + "_calcinterpDiffs.txt";
	outp.open(file.c_str());

    for (int is=0; is<nSED; is++) {
            
        // set correct reddening law
        int law = 0; // The Cardelli law
        if ( is >= (nElliptical + nSpiral))
            law = 1; // The Calzetti law

        // Current SED
        SED sed(*(sedArray[is]));
        
        // loop over filters
        for (int ift=0; ift<nFilter; ift++) {
        
            Filter filter((*filters[ift]));
	
	        for (int i=0; i<nz; i++) {
	            double zs = ztmin + i*dzt;
	            
	            for (int j=0; j<ne; j++) {
	                double ext = etmin + j*det;
	                
	                tm.Split();
                    double k1 = simdataCalc.calcKcorr(sed, filter, restFrameFilter, zs, ext, law);
                    tm.Split();
                    calcTime += tm.PartialElapsedTimems();
                   
                    
                    tm.Split();
                    double k2 = simdataInterp.interpKcorr(is, ift, zs, ext);
                    tm.Split();
                    interpTime += tm.PartialElapsedTimems();
                    
                    /*if (is<1&&ift<1&&i<1&&j<1) {
                        cout <<"     Time to calculate k-correction via calculation = "<< calcTime <<" ms,";
                        cout <<" via interpolation = "<< interpTime <<" ms"<< endl;
                        }*/
                    
                    double diff = (k2-k1)/k1;
                    
                    if (abs(diff) > 1e-4 )
                        cnt1++;
                    /*if (abs(diff) > 1e-4 ){
                        cout <<" NOTE! sed="<< is <<", filt="<< ift <<", z="<< zs;
                        cout <<", ex="<< ext<<", kdiff="<< (k2-k1)/k1 << endl; cnt1++; }*/
                    if (abs(diff) > 1e-2 ){
                        cout <<" !!!!!!!!!!!!!!!!!!!!!WARNING! sed="<< is <<", filt=";
                        cout << ift <<", z="<< zs <<", ex="<< ext <<", kdiff="<< (k2-k1)/k1 << endl;
                        cnt2++; 
                        sumSED+=is;
                        sumSEDsq+=is*is;
                        sumFilt+=ift;
                        sumFiltsq+=ift*ift;
                        sumz+=zs;
                        sumzsq+=zs*zs;
                        sume+=ext;
                        sumesq+=ext*ext;
                        
                        }
                    outp << is <<"  "<< ift <<"  "<< zs <<"  "<< ext <<"  "<< diff <<"  "<< k1 << endl;
                    }
                }
            }
        }
    outp.close();
    cout << endl;
    
    int tot = nSED*nFilter*nz*ne;
    cout << "     "<< (100.*(double)cnt1/(double)tot) <<"% >0.01% difference"<<endl;
    cout << "     "<< (100.*(double)cnt2/(double)tot) <<"% >1% difference"<<endl;
    cout << endl;
    
    double meanSED = sumSED/cnt2;
    double stdSED = sqrt(sumSEDsq/cnt2 - meanSED*meanSED);
    cout << "     Mean SED type where k-corrections were >1% different = "<< meanSED<<" std = "<< stdSED << endl;
    double meanFilt = sumFilt/cnt2;
    double stdFilt = sqrt(sumFiltsq/cnt2 - meanFilt*meanFilt);
    cout << "     Mean filter where k-corrections were >1% different = "<< meanFilt<<" std = "<< stdFilt << endl;
    double meanz = sumz/cnt2;
    double stdz = sqrt(sumzsq/cnt2 - meanz*meanz);
    cout << "     Mean redshift where k-corrections were >1% different = "<< meanz<<" std = "<< stdz << endl;
    double meane = sume/cnt2;
    double stde = sqrt(sumesq/cnt2 - meane*meane);
    cout << "     Mean extinction where k-corrections were >1% different = "<< meane<<" std = "<< stde << endl;
    cout << endl;
    
    cout << "     Time to do calculation = "<< calcTime <<" ms "<<endl;
    cout << "     Time to do interpolation = "<< interpTime <<" ms "<<endl;
    cout << endl;
    if ( interpTime> calcTime)
        cout <<"     Interpolation takes "<< (double)(interpTime/calcTime) <<" times longer"<<endl;
    else
        cout <<"     Calculation takes "<< (double)(calcTime/interpTime) <<" times longer"<<endl;
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testKcorrMethod.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testKcorrMethod.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " testKcorrMethod.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of testKcorrMethod.cc program  Rc= " << rc << endl;
  return rc;	
}

