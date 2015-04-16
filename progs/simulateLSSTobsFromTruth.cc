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

#include "hpoly.h"
//#include "genericfunc.h"
#include "constcosmo.h"
#include "cosmocalcs.h"
#include "igm.h"
#include "mydefrg.h"
#include "sedfilter.h"
#include "simdata.h"


#define PI 3.141592
/*



*/
void usage(void);
void usage(void) {
	cout << endl<<" Usage: simulateLSSTobsFromTruth [...options...]" << endl<<endl;
	cout << "  Input catalog must have no header and have the first 7 columns be:"<<endl;
	cout << "  u,g,r,i,z,y,zs "<<endl<<endl;

	cout << " -i INCAT:   name of ImSim reference catalog to read in "<<endl;
	cout << " -f NFILTER: number of filters in catalog [INCAT] [DEFAULT=6]"<<endl;
	cout << " -n NGAL:    number of galaxies in catalog [INCAT] "     <<endl;
	cout << " -o OUTCATROOT:  name of catalog to write [INCAT] data plus simulated";
	cout << "             observed magnitudes and errors to (does both .txt and .fits)"<<endl;
	cout << " -v NYEAR:   number of years of LSST operations [DEFAULT=10]"<<endl;
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

    cout << " ==== simulateLSSTobsFromTruth.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string infile, outfile;
    int nFilter = 6;
    long nGalaxy;
    int nYear = 10;
    
    // Number of visits per year (Table 1, Ivezic et al 2008)
    /*int uVisitsPerYear = 6;
    int gVisitsPerYear = 8;
    int rVisitsPerYear = 18;
    int iVisitsPerYear = 18;
    int zVisitsPerYear = 16;
    int yVisitsPerYear = 16;*/
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hi:f:n:o:v:")) != -1) {
	    switch (c) {
	        case 'i' :
	            infile = optarg;
	            break;
	        case 'f' :
	            sscanf(optarg,"%d",&nFilter);
		        break;
	        case 'n' :
	            sscanf(optarg,"%ld",&nGalaxy);
		        break;
	        case 'o' :
	            outfile = optarg;
	            outfile = "output/" + outfile;
	            break;
	        case 'v' :
	            sscanf(optarg,"%d",&nYear);
		        break;
	        case 'h' :
		        default :
		        usage(); return -1;
		    }
	    }
	    
	// output file names
	string outfileFits = outfile + ".fits";
	string outfileText = outfile + ".txt";
	    
	// total number of visits
	/*int uVisits = uVisitsPerYear*nYear;
    int gVisits = gVisitsPerYear*nYear;
    int rVisits = rVisitsPerYear*nYear;
    int iVisits = iVisitsPerYear*nYear;
    int zVisits = zVisitsPerYear*nYear;
    int yVisits = yVisitsPerYear*nYear;*/
    
    //-- end command line arguments
    cout <<"     Reading in "<< nGalaxy <<" galaxies from catalog "<< infile <<endl;
    cout <<"     Writing out data plus observed magnitudes and errors to files ";
    cout << outfileFits <<" and "<< outfileText <<endl;
    cout <<"     for "<< nYear <<" years of LSST operations"<<endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    ifstream inp;
	ofstream outp;
	
	
	SimpleUniverse su;
	RandomGenerator rg;
	
	// Load in LSST filters
	cout <<"     Load in LSST filters"<<endl;
	double lmin=5e-8, lmax=2.5e-6;
	string filterFile = "LSST.filters";
	ReadFilterList readLSSTfilters(filterFile);
	readLSSTfilters.readFilters(lmin,lmax);
	vector<Filter*> filterArray=readLSSTfilters.getFilterArray();
	
	// Load in SEDs (don't actually need this)
	string sedFile = "CWWK.list";
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lmin,lmax);
    vector<SED*> sedArray=readSedList.getSedArray();
	
	// Initialize class that simulates magnitude errors
	SimObservations simData(filterArray, rg);
	
	// Open file for reading
	cout <<"    Reading in the LSST catalog ... "<< nGalaxy <<" entries from ";
    cout << nFilter+1 <<" columns"<<endl;
    inp.open(infile.c_str());
    
    // Open text file for writing
    outp.open(outfileText.c_str());
    
    // Create swap space FITS file structure
    FitsInOutFile swf(outfileFits, FitsInOutFile::Fits_Create);	
    SwFitsDataTable catalog(swf, 2048);
    catalog.AddFloatColumn("ref_umag");
    catalog.AddFloatColumn("ref_gmag");
    catalog.AddFloatColumn("ref_rmag");
    catalog.AddFloatColumn("ref_imag");
    catalog.AddFloatColumn("ref_zmag");
    catalog.AddFloatColumn("ref_ymag");
    catalog.AddFloatColumn("z");
    catalog.AddFloatColumn("obs_umag");
    catalog.AddFloatColumn("err_umag");
    catalog.AddFloatColumn("obs_gmag");
    catalog.AddFloatColumn("err_gmag");
    catalog.AddFloatColumn("obs_rmag");
    catalog.AddFloatColumn("err_rmag");
    catalog.AddFloatColumn("obs_imag");
    catalog.AddFloatColumn("err_imag");
    catalog.AddFloatColumn("obs_zmag");
    catalog.AddFloatColumn("err_zmag");
    catalog.AddFloatColumn("obs_ymag");
    catalog.AddFloatColumn("err_ymag");
    DataTableRow row = catalog.EmptyRow();
    int NColfits=catalog.NCols();
     
    int nFilter = 6;
    
    int cntU=0., cntG=0., cntR=0., cntI=0., cntZ=0., cntY=0.;
    char a[80];
    for(int i=0; i<nGalaxy; i++) {
	    cout <<"     Entry "<< i+1 <<" of "<< nGalaxy <<endl;
        //cout << "     ";
        // read in magnitudes
        vector<double> mags;
        for(int j=0; j<nFilter; j++) {
		    inp >> a;
		    double m = atof(a);
		    //cout << m <<"  ";
		    mags.push_back(m);
		    }
		    
		// read in redshift
		inp >> a;
		double zs = atof(a);
		//cout << zs <<"  ";
		
		
		vector<double> uObs = simData.addLSSTError(mags[0], nYear, 0);
		vector<double> gObs = simData.addLSSTError(mags[1], nYear, 1);
		vector<double> rObs = simData.addLSSTError(mags[2], nYear, 2);
		vector<double> iObs = simData.addLSSTError(mags[3], nYear, 3);
		vector<double> zObs = simData.addLSSTError(mags[4], nYear, 4);
		vector<double> yObs = simData.addLSSTError(mags[5], nYear, 5);
		
		if (my_isnan(uObs[0])) {
		    uObs[0] = 99;
		    cout <<"     u band magnitude is NaN: uErr = "<< uObs[1] <<", zs = "<< zs <<endl;
		    cntU++;
		    }
		if (my_isnan(gObs[0])) {
		    gObs[0] = 99;
		    cout <<"     g band magnitude is NaN: gErr = "<< gObs[1] <<", zs = "<< zs <<endl;
		    cntG++;
		    }
		if (my_isnan(rObs[0])) {
		    rObs[0] = 99;
		    cout <<"     r band magnitude is NaN: rErr = "<< rObs[1] <<", zs = "<< zs <<endl;
		    cntR++;
		    }
		if (my_isnan(iObs[0])) {
		    iObs[0] = 99;
		    cout <<"     i band magnitude is NaN: iErr = "<< iObs[1] <<", zs = "<< zs <<endl;
		    cntI++;
		    }
		if (my_isnan(zObs[0])) {
		    zObs[0] = 99;
		    cout <<"     z band magnitude is NaN: zErr = "<< zObs[1] <<", zs = "<< zs <<endl;
		    cntZ++;
		    }
		if (my_isnan(yObs[0])) {
		    yObs[0] = 99;
		    cout <<"     y band magnitude is NaN: yErr = "<< yObs[1] <<", zs = "<< zs <<endl;
		    cntY++;
		    }
		    
		    
		//cout << uObs[0] <<"  "<< uObs[1];
		//cout << endl;
		
		// write original data
		for(int j=0; j<nFilter; j++){
		    outp << mags[j] <<"  ";
		    row[j] = mags[j];
		    }
		outp << zs <<"  "; 
		int jj = nFilter;
		row[jj] = zs;
		
		// write new data to text file
		outp << uObs[0] <<"  "<< uObs[1] <<"  ";
		outp << gObs[0] <<"  "<< gObs[1] <<"  ";
		outp << rObs[0] <<"  "<< rObs[1] <<"  ";
		outp << iObs[0] <<"  "<< iObs[1] <<"  ";
		outp << zObs[0] <<"  "<< zObs[1] <<"  ";
		outp << yObs[0] <<"  "<< yObs[1] <<"  ";
		outp << endl;
		
		// write data to FITS file
		jj++; row[jj] = uObs[0]; jj++; row[jj] = uObs[1];
		jj++; row[jj] = gObs[0]; jj++; row[jj] = gObs[1];
		jj++; row[jj] = rObs[0]; jj++; row[jj] = rObs[1];
		jj++; row[jj] = iObs[0]; jj++; row[jj] = iObs[1];
		jj++; row[jj] = zObs[0]; jj++; row[jj] = zObs[1];
		jj++; row[jj] = yObs[0]; jj++; row[jj] = yObs[1];
		
		catalog.AddRow(row);
		
		
		}
		
    inp.close();
	cout << endl;
	
	if (cntU>0)
	    cout <<"     "<< cntU <<" galaxies with zero flux in u band"<<endl;
	if (cntG>0)
	    cout <<"     "<< cntG <<" galaxies with zero flux in g band"<<endl;
	if (cntR>0)
	    cout <<"     "<< cntR <<" galaxies with zero flux in r band"<<endl;
	if (cntI>0)
	    cout <<"     "<< cntI <<" galaxies with zero flux in i band"<<endl;
	if (cntZ>0)
	    cout <<"     "<< cntZ <<" galaxies with zero flux in z band"<<endl;
	if (cntY>0)
	    cout <<"     "<< cntY <<" galaxies with zero flux in y band"<<endl;    
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " simulateLSSTobsFromTruth.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " simulateLSSTobsFromTruth.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " simulateLSSTobsFromTruth.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of simulateLSSTobsFromTruth.cc program  Rc= " << rc << endl;
  return rc;	
}

