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
#include "ctimer.h"

#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "simdata.h"

void usage(void);
void usage(void)
{
	cout << endl<<" Usage: calculateKcorrections [...options...]" << endl<<endl;

	cout << " -s : SEDFILE: reading model galaxy SEDs from file SEDFILE. Must be in order"<<endl;
	cout << "      of: elliptical types, spiral types, starburst types [DEFAULT=CWWK.list]" <<endl;
	cout << " -f : FILTFILE: reading filters from file FILTFILE [DEFAULT=LSST.filters]"<<endl;
	//cout << " -t : NELLIPTICAL,NSPIRAL: number of elliptical, spiral SEDs [DEFAULT=1,2]"<<endl;
	cout << " -z : ZMIN,ZMAX,NZ: range of k-correction calculation in redshift [DEFAULT=0,3,2000]"<<endl;
	cout << " -e : EMAX,NE: range of k-correction calculation in host galaxy extinction";
	cout << " [DEFAULT=0.3,200]"<<endl;
	cout << endl;
}
int main(int narg, char* arg[])
{
	cout << " ==== calculateKcorrections.cc program , to simulate LSST data ==== "<<endl;

	// make sure SOPHYA modules are initialized 
    SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;

	//--- decoding command line arguments 
	string outfileroot = "kCorrections/kCorr";
	string sedFile = "CWWK.list";
	string filterFile = "LSST.filters";

    //int nElliptical = 1;
    //int nSpiral = 2;
    
    double zmin=0, zmax=3.;
    double emax=0.3;
    int nz = 2000;
    int ne = 200;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hs:f:z:e:")) != -1)  {
	    switch (c) {
	        case 's' :
		        sedFile = optarg;
		        break;
            case 'f' :
		        filterFile = optarg;
		        break;
		    //case 't' :
		    //    sscanf(optarg,"%d,%d",&nElliptical,&nSpiral);
		    //    break;
		    case 'z' :
		        sscanf(optarg,"%lf,%lf,%d",&zmin,&zmax,&nz);
		        break;
		    case 'e' :
		        sscanf(optarg,"%lf,%d",&emax,&ne);
		        break;
	      case 'h' :
		    default :
		    usage(); return -1;
		    }
	    }
	    
	string delim=".";
	vector<string> result1,result2;
	stringSplit(sedFile,delim,result1);
	stringSplit(filterFile,delim,result2);
	string sedLib = result1[0];
	string filtLib = result2[0];
	
    cout << "     Reading SEDs from file "<< sedFile <<", and filters from "<< filterFile << endl;
    //cout << "     Number of ellipticals = "<< nElliptical <<", number of spirals = ";
    //cout << nSpiral << endl;    
    cout << "     Files will be output to files beginning "<< outfileroot <<endl;
    cout << "     Calculating k-corrections over redshift range "<< zmin <<"<z<"<< zmax;
    cout << " with "<< nz <<" steps"<<endl;
    cout << "     And over a host galaxy extinction range up to "<< emax <<" with "<< ne <<" steps "<<endl;
    cout << endl;
    //-- end command line arguments
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
	ResourceUsage res;
	InitTim();
	
	//enum dustLaw{NoDust=0, Card=1, Calz=2};
	
	// this controls the drawing of random numbers
	RandomGenerator rg;

	// Set cosmology
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout <<"     Set cosmology to: OmegaM="<<OmegaM<<", OmegaL="<<OmegaL;
	cout <<", H0="<<100*h<<endl;
	cout << endl;

	// Load in filters required
	// wavelength range of the SEDs/filters
	cout <<"     Load in filters"<<endl;
	double lmin=5e-8, lmax=2.5e-6;
	ReadFilterList readFilterList(filterFile);
	readFilterList.readFilters(lmin, lmax);
	vector<Filter*> filters = readFilterList.getFilterArray();
	int nFilter = readFilterList.getNTot();
	cout <<"     Read in "<< nFilter <<" filters"<<endl;
	cout <<endl;
	
	// Load in the rest-frame filter
	cout <<"     Load in GOODS B filter (the rest-frame filter)"<<endl;
	string rfFilter = "GOODSB";
	string goodsFilterFile = rfFilter + ".filters";
	ReadFilterList readGOODSBfilter(goodsFilterFile);
	readGOODSBfilter.readFilters(lmin, lmax);
	vector<Filter*> goodsBFilter=readGOODSBfilter.getFilterArray();
	Filter restFrameFilter((*goodsBFilter[0]));
	
	// Load in SEDs
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax); // Read out SEDs into array  
    int nSED = readSedList.getNTot();
    vector<SED*> sedArray=readSedList.getSedArray();
    cout <<"     Final total number of SEDs = "<< nSED <<endl;
	cout << endl;

	// Prepare the class which will calculate the magnitudes
	cout <<"     Initialize class to calculate magnitudes"<<endl;
	SimData simData(sedArray, filters, su); //, nElliptical, nSpiral);
	//SimData simgalNoMadau(sedArray, filters, su, rg); //, nElliptical, nSpiral);
	//simgalNoMadau.setMadau(false);
    cout << endl;
    
    // IGM models
    //cout <<"     Set IGM models"<<endl
    int npt = 1000;
    IGMTransmission noIGM(lmin, lmax, npt);
    
    //cout << endl;
    
    // Step in redshift and host galaxy extinction
    double dz = (zmax - zmin)/(nz - 1);
    double de = emax/(ne - 1);

	cout <<"     Start loop over SED, filter, redshift and host galaxy extinction ..."<<endl;
	Timer tm("TimingGalaxies",false);
	Timer tm2("TimingWhole",false);
	tm.Split(); tm2.Split();
    cout << endl;
    
    ofstream outp1, outp2;    
    // loop over SEDs
    for (int is=0; is<nSED; is++) {
    
        cout <<"     SED "<< is+1 <<" of "<< nSED <<endl;
            
        // set correct reddening law
        dustLaw law = Card;
        
        //int law = 0; // The Cardelli law
        //if ( is >= (nElliptical + nSpiral))
        //    law = 1; // The Calzetti law

        // Current SED
        SED sed(*(sedArray[is]));
        
        // loop over filters
        for (int ift=0; ift<nFilter; ift++) {
        
            cout <<"     Filter "<< ift+1 <<" of "<< nFilter <<endl;
        
            // Current filter
            Filter filter((*filters[ift]));
            
            // Make filename
            stringstream ss1,ss2,ss3,ss4,ss5,ss6,ss7;
            ss1 << is; 
            ss2 << ift;
            ss3 << zmin; ss4 << zmax; ss5 << nz;
            ss6 << emax; ss7 << ne;
            string fnameMadau, fnameNoMadau;
            fnameNoMadau += outfileroot + "_" + sedLib + "sed" + ss1.str();
            fnameNoMadau += "_" + filtLib + "filt" + ss2.str() + "_" + rfFilter;
            fnameNoMadau += "_zmin" + ss3.str() + "_zmax" + ss4.str() + "_nz" + ss5.str();
            fnameNoMadau += "_emax" + ss6.str() + "_ne" + ss7.str();
            
            
            fnameMadau = fnameNoMadau + "_wMadau.txt";
            fnameNoMadau = fnameNoMadau + "_woMadau.txt";
            
            // Open file
            outp1.open(fnameMadau.c_str());
            outp2.open(fnameNoMadau.c_str());
            
            // Write header
            outp1 << "# z range: zmin=" << zmin << endl;
            outp1 << "# z range: dz="<< dz << endl;
            outp1 << "# z range: zmax="<< zmax << endl;
            outp1 << "# E(B-V) range: emin=0"<< endl;
            outp1 << "# E(B-V) range: de="<< de << endl;
            outp1 << "# E(B-V) range: emax="<< emax << endl;
            outp2 << "# z range: zmin=" << zmin << endl;
            outp2 << "# z range: dz="<< dz << endl;
            outp2 << "# z range: zmax="<< zmax << endl;
            outp2 << "# E(B-V) range: emin=0"<< endl;
            outp2 << "# E(B-V) range: de="<< de << endl;
            outp2 << "# E(B-V) range: emax="<< emax << endl;
            
            double kcorrWithMadau, kcorrNoMadau;
	        for (int i=0; i<nz; i++) {
	            double zs = zmin + i*dz;
	            IGMTransmission madauIGM(zs);
	            
	            for (int j=0; j<ne; j++) {
	            
	                
                    double ex = j*de; 
                    
                    double kcorrWithMadau = simData.calcKcorr(sed, filter, restFrameFilter, zs, madauIGM, ex, law);
		            double kcorrNoMadau = simData.calcKcorr(sed, filter, restFrameFilter, zs, noIGM, ex, law);
		            
		            //double calcKcorr(SED& sed, Filter& filter, Filter& restFrameFilter, double zs, 
	                // IGMTransmission& igmtrans, double ext, dustLaw law);
	                 
                    
                    outp1 << kcorrWithMadau <<"  ";
                    outp2 << kcorrNoMadau <<"  ";
                    
                                    
                    } // end of loop over extinction values
                outp1 << endl;
                outp2 << endl;
                
                tm.Split();
                if (is<1 && i<1) {
                    if (tm.PartialElapsedTime()>0)
	                    cout <<"     Time per redshift "<< tm.PartialElapsedTime() <<" s:";
	                else
	                    cout <<"     Time per redshift "<< tm.PartialElapsedTimems() <<" ms:";
	                cout <<" number of redshift to do = "<< nz*nFilter*nSED <<endl;
	                }
                }// end of loop over redshift values
                
            // Close the files
            outp1.close();
            outp2.close();
             
            }
        cout << endl;
		}
		
	cout <<"     End loops "<<endl;
	tm2.Split();
	cout <<"     .... done, took "<< tm2.PartialElapsedTime()/60 <<" mins";
	
	
	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " calculateKcorrections.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " calculateKcorrections.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " calculateKcorrections.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of calculateKcorrections.cc program  Rc= " << rc << endl;
  return rc;	
}
