// -*- LSST-C++ -*-
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
#define PI 3.141592

void usage(void);
void usage(void)
{
	cout << endl<<" Usage: testSimulation [...options...]" << endl<<endl;

	cout << " -o : OUTFILE: output filename for observed LSST catalog (will be stored in testfiles/)"<<endl;
	cout << " -s : SEDFILE: reading model galaxy SEDs from file SEDFILE [DEFAULT=CWWK.list]"<<endl;
	cout << " -t : NELLIPTICAL,NSPIRAL: number of elliptical, spiral SEDs [DEFAULT=1,2]"<<endl;
	cout << " -z : zfile: read cumulative redshift distribution from file";
	cout << "  [ZFILE]"<<endl;
	cout << " -m : mfile: read cumulative mag distribution from file";
	cout << "  [MFILE]"<<endl;
	cout << endl;
}
int main(int narg, char* arg[])
{
	cout << " ==== testSimulation.cc program , to simulate test data ==== "<<endl;

	// make sure SOPHYA modules are initialized 
    SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;

	//--- decoding command line arguments 
	string outfile="testfiles/testSimulation.fits";
	string sedFile = "CWWK.list";
	string zfile,mfile;
    int nElliptical = 1;
    int nSpiral = 2;
    
    double zMax = 1.;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hi:t:z:m:o:")) != -1)  {
	    switch (c) {
	        case 'o' :
		        outfile = optarg;
		        outfile = "testfiles/" + outfile;
		        break;
	        case 's' :
		        sedFile = optarg;
		        break;
		    case 't' :
		        sscanf(optarg,"%d,%d",&nElliptical,&nSpiral);
		        break;
		    case 'z' :
			    zfile = optarg;
			    break;
		    case 'm' :
			    mfile = optarg;
			    break;
	      case 'h' :
		    default :
		    usage(); return -1;
		    }
	    }
	    
	int nYear = 10;
	// Number of visits per year (Table 1, Ivezic et al 2008)
    int uVisitsPerYear = 6;
    int gVisitsPerYear = 8;
    int rVisitsPerYear = 18;
    int iVisitsPerYear = 18;
    int zVisitsPerYear = 16;
    int yVisitsPerYear = 16;
    
	// total number of visits
	int uVisits = uVisitsPerYear*nYear;
    int gVisits = gVisitsPerYear*nYear;
    int rVisits = rVisitsPerYear*nYear;
    int iVisits = iVisitsPerYear*nYear;
    int zVisits = zVisitsPerYear*nYear;
    int yVisits = yVisitsPerYear*nYear;
	
    cout << "     Reading SEDs from file "<< endl;
    cout << "     Number of ellipticals = "<< nElliptical <<", number of spirals = ";
    cout << nSpiral << endl;
    cout << "     Output catalog is "<< outfile <<endl;
    cout << endl;
    //-- end command line arguments
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
	ResourceUsage res;
	InitTim();
	
	// Output file
	cout <<"     Creating output file "<<outfile<<endl;
	FitsInOutFile swf(outfile,FitsInOutFile::Fits_Create);
	
	// binary data table
	SwFitsDataTable gals(swf,2048);
	gals.AddFloatColumn("zs");
	gals.AddFloatColumn("am");
	gals.AddFloatColumn("type");
	gals.AddFloatColumn("ext");
	gals.AddFloatColumn("mu");
	gals.AddFloatColumn("mg");
	gals.AddFloatColumn("mr");
	gals.AddFloatColumn("mi");
	gals.AddFloatColumn("mz");
	gals.AddFloatColumn("my");
	gals.AddFloatColumn("muo");
	gals.AddFloatColumn("mgo");
	gals.AddFloatColumn("mro");
	gals.AddFloatColumn("mio");
	gals.AddFloatColumn("mzo");
	gals.AddFloatColumn("myo");
	gals.AddFloatColumn("muoo");
	gals.AddFloatColumn("mgoo");
	gals.AddFloatColumn("mroo");
	gals.AddFloatColumn("mioo");
	gals.AddFloatColumn("mzoo");
	gals.AddFloatColumn("myoo");
	gals.AddFloatColumn("muooo");
	gals.AddFloatColumn("mgooo");
	gals.AddFloatColumn("mrooo");
	gals.AddFloatColumn("miooo");
	gals.AddFloatColumn("mzooo");
	gals.AddFloatColumn("myooo");
	DataTableRow rowin=gals.EmptyRow();
	cout << endl;

	// This controls the drawing of random numbers
	RandomGenerator rg;

	// Cosmological parameters
	cout <<"     Set cosmology ..."<<endl;
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout << "     h = "<<h<<", OmegaM = "<<OmegaM<<", OmegaL = "<<OmegaL;
	cout <<endl<<endl;
	
	// Load in filters required
	// wavelength range of the SEDs/filters
	cout <<"     Load in LSST filters"<<endl;
	double lmin=5e-8, lmax=2.5e-6;
	string filterFile = "LSST.filters";
	ReadFilterList readFilterList(filterFile);
	readFilterList.readFilters(lmin,lmax);
	vector<Filter*> LSSTfilters=readFilterList.getFilterArray();
	int nFilter = readFilterList.getNTot();
	cout <<"     Read in "<< nFilter <<" filters"<<endl;
	//Filter restFrameFilter((*GOODSfilters[1]));
	int iU = 0;
	int iG = 1;
	int iR = 2;
	int iI = 3;
	int iZ = 4;
	int iY = 5;
	cout <<endl;
	
	cout <<"     Load in GOODS B filter"<<endl;
	string goodsFilterFile = "GOODSB.filters";
	ReadFilterList readGOODSBfilter(goodsFilterFile);
	readGOODSBfilter.readFilters(lmin,lmax);
	vector<Filter*> goodsBFilter=readGOODSBfilter.getFilterArray();
	
	// Load in SEDs
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax); // Read out SEDs into array
    int nsed=readSedList.getNSed(); // Get total number of SEDs  
    int ntot = readSedList.getNTot();
    vector<SED*> sedArray=readSedList.getSedArray();
    cout <<"     Final total number of SEDs = "<<ntot<<endl;
	cout << endl;

	// LF parameters
	cout << "     Initializing LF parameters as a function of redshift"<<endl;
	double zmax = 6., dz = 0.01;
	LFParameters lfpars(zmax, dz);
	cout << endl;
	
	// Redshift distributions
	cout <<"     Setting up cumulative redshift distribution"<<endl;
	CumulDistZ cumlz;
	cout <<"     ... reading cumulative redshift distribution from";
	cout <<" file"<<endl;
	cumlz.SetUp(zfile);
	cout <<endl;

	cout << "     Setting up galaxy survey simulation"<<endl;
	// number of galaxies in the survey
	long ng = 10000;
	cout<<endl;
		
	// galaxy fraction stuff
	vector<double> all=lfpars.ReturnAll();
	Schechter schA(all[0],all[1],all[2]);
	vector<double> early=lfpars.ReturnEarly();
	Schechter schE(early[0],early[1],early[2]);
	vector<double> late=lfpars.ReturnLate();
	Schechter schL(late[0],late[1],late[2]);
	vector<double> sb=lfpars.ReturnSB();
	Schechter schS(sb[0],sb[1],sb[2]);
	TypeRatio0 tr0(schA,schE,schL,schS);
	TypeRatio tr(tr0);
		
	// Redshift distribution to draw from
	cout <<" before drawz"<<endl;
	DrawZ drz(cumlz,rg); 
	cout <<" after drawz"<<endl;
	cout << endl;
	
	// Magnitude distribution to draw from
	CumulDistM cumlm(lfpars,su);
	DrawM drm(cumlm);
    drm.SetUp(mfile);
		
	// Type distribution to draw from
	cout <<"     Initialize type distribution"<<endl;
	DrawType drt(tr,rg);
	cout << endl;

	// Set up base simulation
	SimBaseCatalog sbc(drz,drm,drt);
	

	// Prepare the class which will calculate the magnitudes
	cout <<"     Initialize class to calculate magnitudes"<<endl;
	SimData simgal(sedArray,LSSTfilters,su,rg,nElliptical,nSpiral);
	SimData simgalMadau(sedArray,LSSTfilters,su,rg,nElliptical,nSpiral);
    cout << endl;
    
    // Add Madau preference
    simgal.setMadau(false);
    
    double fluxError = 0.1;

    // Loop over all galaxies in the base catalog
	cout <<"     Start loop over galaxies ..."<<endl;

	for (int i=0; i<ng; i++)
		{
		cout <<"     gal "<<i+1<<" of "<<ng<<endl;

        // Do base simulation
		double zs, am, typDouble;
		sbc.DoSim(zs,am,typDouble);
		int typ=(int)typDouble;
		double type = simgal.SimSED(typ);
        double ext = simgal.SimRed(type);

        // Calculate galaxy magnitude in observed filters ugrizy
        // Galaxy's absolute magnitude is defined in filter (*goodsBFilter[0])
		double uMagTh=simgal.GetMag(zs,type,am,ext,iU,(*goodsBFilter[0]));
		double gMagTh=simgal.GetMag(zs,type,am,ext,iG,(*goodsBFilter[0]));
		double rMagTh=simgal.GetMag(zs,type,am,ext,iR,(*goodsBFilter[0]));
		double iMagTh=simgal.GetMag(zs,type,am,ext,iI,(*goodsBFilter[0]));
		double zMagTh=simgal.GetMag(zs,type,am,ext,iZ,(*goodsBFilter[0]));
		double yMagTh=simgal.GetMag(zs,type,am,ext,iY,(*goodsBFilter[0]));
		
        double uMagThNoExt=simgal.GetMag(zs,type,am,0.,iU,(*goodsBFilter[0]));
		double gMagThNoExt=simgal.GetMag(zs,type,am,0.,iG,(*goodsBFilter[0]));
		double rMagThNoExt=simgal.GetMag(zs,type,am,0.,iR,(*goodsBFilter[0]));
		double iMagThNoExt=simgal.GetMag(zs,type,am,0.,iI,(*goodsBFilter[0]));
		double zMagThNoExt=simgal.GetMag(zs,type,am,0.,iZ,(*goodsBFilter[0]));
		double yMagThNoExt=simgal.GetMag(zs,type,am,0.,iY,(*goodsBFilter[0]));
		
		vector<double> uObservation = simgal.addError(uMagThNoExt,fluxError,0);
		vector<double> gObservation = simgal.addError(gMagThNoExt,fluxError,1);
		vector<double> rObservation = simgal.addError(rMagThNoExt,fluxError,2);
		vector<double> iObservation = simgal.addError(iMagThNoExt,fluxError,3);
		vector<double> zObservation = simgal.addError(zMagThNoExt,fluxError,4);
		vector<double> yObservation = simgal.addError(yMagThNoExt,fluxError,4);
		
		vector<double> uObservation2 = simgal.addLSSTuError(uMagThNoExt,uVisits);
        vector<double> gObservation2 = simgal.addLSSTgError(gMagThNoExt,gVisits);
        vector<double> rObservation2 = simgal.addLSSTrError(rMagThNoExt,rVisits);
        vector<double> iObservation2 = simgal.addLSSTiError(iMagThNoExt,iVisits);
        vector<double> zObservation2 = simgal.addLSSTzError(zMagThNoExt,zVisits);
        vector<double> yObservation2 = simgal.addLSSTyError(yMagThNoExt,yVisits);
		
        // Write the data to the FITS file
		rowin[0]=zs;
		rowin[1]=am;
		rowin[2]=type;
		rowin[3]=ext;
		rowin[4]=uMagTh;
		rowin[5]=gMagTh;
		rowin[6]=rMagTh;
		rowin[7]=iMagTh;
		rowin[8]=zMagTh;
		rowin[9]=yMagTh;
        rowin[10]=uMagThNoExt;
		rowin[11]=gMagThNoExt;
		rowin[12]=rMagThNoExt;
		rowin[13]=iMagThNoExt;
		rowin[14]=zMagThNoExt;
		rowin[15]=yMagThNoExt;
		rowin[16]=uObservation[0];
		rowin[17]=gObservation[0];
		rowin[18]=rObservation[0];
		rowin[19]=iObservation[0];
		rowin[20]=zObservation[0];
		rowin[21]=yObservation[0];
		rowin[22]=uObservation2[0];
		rowin[23]=gObservation2[0];
		rowin[24]=rObservation2[0];
		rowin[25]=iObservation2[0];
		rowin[26]=zObservation2[0];
		rowin[27]=yObservation2[0];
		gals.AddRow(rowin);
		}
	cout <<"     End loop"<<endl;

	
	// Write information on simulation to FITS header
	DVList  dvl;
	dvl("SedFile") = sedFile;
	dvl("NEllip") = nElliptical;
	dvl("NSpiral") = nSpiral;
	swf.WriteHeaderRecords(dvl);
	swf.MoveAbsToHDU(2);
	
	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testSimulation.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testSimulation.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " testSimulation.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of testSimulation.cc program  Rc= " << rc << endl;
  return rc;	
}
