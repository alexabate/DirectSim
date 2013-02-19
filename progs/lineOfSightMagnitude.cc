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
	cout << endl<<" Usage: lineOfSightMagnitude [...options...]" << endl<<endl;
    cout << " -o: OUTFILEROOT: write files to filename beginning OUTROOT (saved to output/)"<<endl;
	cout << " -n: NLINES: simulate [NLINES] lines of sight "<<endl;
	cout << " -z: ZSOURCE: redshift to simulate line of sight to"<<endl;
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

	cout << " ==== lineOfSightMagnitude.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfileroot = "output/lineOfSightMagnitude";
    string infileroot;
    int nLines = 1000;
    double zSource = 3.5;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ho:i:n:z:")) != -1) {
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
            case 'h' :
                default :
                usage(); return -1;
            }
	    }

    //-- end command line arguments
    cout <<"     Simulating "<< nLines <<" lines of sight out to z = "<< zSource <<endl;
    cout <<"     Writing to files beginning "<< outfileroot <<endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    ifstream inp;
	ofstream outp,outpOther;
	string infile, outfile;  
	
	// Load in SEDs
	double lmin=100e-10, lmax=8000e-10;
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
	cout <<"     Calculate u and g magnitudes for a starburst galaxy at z = "<< zSource;
	cout <<" with varing IGM along the line of sight" <<endl;
	int nElliptical = 1;
    int nSpiral = 2;
	SimData simgal(sedArray,LSSTfilters,su,rg,nElliptical,nSpiral);
	
    cout << endl;

    vector<int> nAbsorbers;
    stringstream ssz;
    ssz << zSource;
    
    outfile = outfileroot + "_magsIGM_zSource" + ssz.str() + ".txt";
	outp.open(outfile.c_str());
    for (int i=0; i<nLines; i++){
    
        stringstream ss;
        ss << i+1;
        //ss2 << zSource;
        
        double zs = zSource;
        double am = -18;
        double ext = 0.;
        double type = 3.007;
        int sedNo = 7;
    
        // read in file with IGM transmission
        infile = infileroot + "_fullTransMeiksin_lineOfSight"+ss.str()+"_toRedshift" + ssz.str() +".txt";
        string line;
        inp.open(infile.c_str());
        getline(inp,line);
        string nAbsString;
        if (line.size()>25)
            nAbsString = line[line.size()-2] + line[line.size()-1];
        else
            nAbsString = line[line.size()-1];
        int nAbs = atoi(nAbsString.c_str());
        //cout <<"     Number of absorbers = "<<nAbs <<endl;

        IGMTransmission igmTransmission(infile);
        SEDIGM sedIGM(*(sedArray[sedNo]), igmTransmission, zs);
        SEDMadau sedMadau(*(sedArray[sedNo]), zs);
        
        /*//////////////// THIS PART IS FOR DEBUGGING/CHECKING 
        int nl =10000;
        double dl = (lmax-lmin)/(nl-1);
        outfile = outfileroot + "_trans" + ss.str() + ".txt";
        outpOther.open(outfile.c_str());
        for (int i=0; i<nl; i++) {
            double lamE = lmin +i*dl;
            double lamO = lamE*(1. + zs);
            outpOther << lamO <<"  "<< lamE <<"  "<< igmTransmission(lamO) <<"  ";
            outpOther << sedArray[sedNo]->returnFlux(lamE) <<"  "<< sedArray[sedNo]->returnFlux(lamO);
            outpOther << "  "<< sedIGM(lamE) << "  "<< sedMadau(lamE) << endl;
            }
        outpOther.close();
        ////////////////*/
        

		double uMag=simgal.GetMag(zs,type,am,ext,0,(*goodsBFilter[0]),igmTransmission);
		double gMag=simgal.GetMag(zs,type,am,ext,1,(*goodsBFilter[0]),igmTransmission);
		//double uMagNoIGM=simgalNoIGM.GetMag(zs,type,am,ext,0,(*goodsBFilter[0]));
		//double gMagNoIGM=simgalNoIGM.GetMag(zs,type,am,ext,1,(*goodsBFilter[0]));
		
		outp << nAbs <<"  "<< uMag <<"  "<< gMag <<"  " <<endl;//<< uMagNoIGM <<"  "<< gMagNoIGM <<endl;
		
		inp.close();
        }
    outp.close();
    cout << endl;
    
    // Now calculate u and g mags with NO igm
    cout <<"     Calculate u and g magnitudes without IGM for the starburst galaxy at ";
    double zmin = 0.01, zmax = zSource+3.;
    int nz = 10000;
    double dz = (zmax - zmin)/(nz - 1);
    cout << zmin <<"<z<"<< zmax <<endl;
    
    SimData simgalNoIGM(sedArray,LSSTfilters,su,rg,nElliptical,nSpiral);
	simgalNoIGM.setMadau(false);
    
    outfile = outfileroot + "_magsZ.txt";
	outp.open(outfile.c_str());
    for (int i=0; i<nz; i++){
    
        double zs = zmin + i*dz;
        double am = -18;
        double ext = 0.;
        double type = 3.007;
        
        double uMag=simgalNoIGM.GetMag(zs,type,am,ext,0,(*goodsBFilter[0]));
		double gMag=simgalNoIGM.GetMag(zs,type,am,ext,1,(*goodsBFilter[0]));
        outp << zs <<"  "<< uMag <<"  "<< gMag <<"  " <<endl;
        }
    outp.close();
    
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " lineOfSightMagnitude.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " lineOfSightMagnitude.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " lineOfSightMagnitude.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of lineOfSightMagnitude.cc program  Rc= " << rc << endl;
  return rc;	
}

