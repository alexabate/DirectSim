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

#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "simdata.h"
#include "sedfilter.h"

void usage(void);
void usage(void) {
	cout << endl<<" Usage: produceCatalogForBPZ [...options...]" << endl<<endl;
	
    cout <<" Produces catalog with photometry that exactly matches BPZ                                "<<endl;
	
	cout <<"  Example usage:                                                                          "<<endl;
	cout <<"  produceCatalogForBPZ -i cat.fits -o out -s brown.list -l /mnt/drive2/soft/bpz-1.99.3/AB/"<<endl;
	cout <<"                      -f u_lsst_etc,g_lsst_etc,r_lsst_etc,i_lsst_etc,z_lsst_etc,y_lsst_etc"<<endl;
	cout << endl;
	
	cout << " -i: FITS: FITS file containing catalog produced by fitSEDsToColors                      "<<endl; 
	cout << " -o: OUTROOT: write files to filenames beginning OUTROOT                                 "<<endl;
	cout << " -s: SEDLIB: file containing list of SED files                                           "<<endl;
	cout << " -l: ABLOC: location of BPZ .AB files                                                    "<<endl;
	cout << " -f: FILTLIST: list of filter filenames without extensions                               "<<endl;
	cout << endl;
    };


int main(int narg, char* arg[]) {
    
    cout << " ==== produceCatalogForBPZ.cc program , produce catalog with photometry that exactly matches BPZ ====" << endl;
    cout << endl;

    // make sure SOPHYA modules are initialized 
    SophyaInit();  
    FitsIOServerInit();
    InitTim();
    cout<<endl<<endl;

    //--- decoding command line arguments 
    string outroot;
    string infile;
    string sedfile="brown.list";
    string abFileLoc="/mnt/drive2/soft/bpz-1.99.3/AB/";
    vector<string> filtlist;
    filtlist.push_back("u_lsst_etc");
    filtlist.push_back("g_lsst_etc");
    filtlist.push_back("r_lsst_etc");
    filtlist.push_back("i_lsst_etc");
    filtlist.push_back("z_lsst_etc");
    filtlist.push_back("y_lsst_etc");
    string filts;
    
    // number of years of operation
    int nYear = 10;

	char c;
    while((c = getopt(narg,arg,"ho:i:s:l:f:")) != -1) {
	    switch (c)  {
	        case 'o' :
	            outroot = optarg;
	            break;
	        case 'i' :
	            infile = optarg;
	            break;
	        case 's' :
	            sedfile = optarg;
	            break;
	        case 'l' :
	            abFileLoc = optarg;
	            break;
	        case 'f' : {
	            filts = optarg;
	            string delim=",";
	            filtlist.clear();
	            stringSplit(filts, delim, filtlist);
	            }
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }

    //-- end command line arguments
    cout <<"     Writing to files beginning "<< outroot <<endl;
    cout <<"     Reading from FITS file " << infile << endl;
    cout <<"     Using SED library "<< sedfile <<endl;
    cout <<"     BPZ .AB files are located at "<< abFileLoc << endl;
    cout <<"     LSST filters: ";
    for (int i=0; i<filtlist.size(); i++)
         cout << filtlist[i] <<", ";
    cout << endl << endl;
    //-- end command line arguments
  
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
	ResourceUsage res;
	InitTim();
	string outfile;

    // Random  number generator
	RandomGenerator rg;
	
	
	// READ SIM CATALOG
	FitsInOutFile fin(infile, FitsInOutFile::Fits_RO);
	fin.MoveAbsToHDU(2);
    SwFitsDataTable dt(fin, 512, false);
    sa_size_t ng = dt.NEntry();
    sa_size_t nc = dt.NCols();
    DataTableRow row = dt.EmptyRow();
    cout <<"     In file "<< infile <<" ... "<<endl;
    cout <<"     Number of columns = "<< nc <<", number of entries = "<< ng << endl;
    cout << endl;
    
		
	// GALAXY SED TEMPLATES
	// wavelength range of the SEDs/filters
	double lmin=5e-8, lmax=2.5e-6;
	ReadSedList readSedList(sedfile);
    readSedList.readSeds(lmin,lmax);
    vector<string> sedFiles = readSedList.returnSedFilenamesNoExt();
    vector<SED*> sedArray = readSedList.getSedArray();
    
    
    // FILTERS
    string filterFile = "LSST.filters";
	ReadFilterList lsstFilters(filterFile);
	lsstFilters.readFilters(lmin, lmax);
	vector<Filter*> lsst_filters = lsstFilters.getFilterArray();
	int nLSST = lsstFilters.getNTot();
	cout <<"     "<< nLSST <<" LSST filters read in "<<endl;
    cout << endl;
    
    
    // FOR GENERATING PHOTOMETRIC ERRORS
    //Cosmological parameters: Om = 0.25, OL = 0.75, Ob = 0.045, s8 = 0.9, h = 0.73,
    double h=0.73, OmegaM = 0.25, OmegaL = 0.75, Ob = 0.045;
	SimpleUniverse su(h, OmegaM, OmegaL);
	su.SetOmegaBaryon(Ob);
	SimData simData(sedArray, lsst_filters, su, rg);
    bool isAddMadau = true;
    simData.setMadau(isAddMadau);


    // CLASS THAT CORRECTS MAGNITUDES
	BpzCorrectMags bpzCorrectMags(abFileLoc, filtlist);
	

    // FILE TO WRITE TO
    outfile = outroot + "_bpz_lsstcat.fits";
    FitsInOutFile swf(outfile, FitsInOutFile::Fits_Create);
    //swf.MoveAbsToHDU(2);
	SwFitsDataTable gals(swf, 2048);
    gals.AddFloatColumn("Z_TRUE");
    gals.AddFloatColumn("SED_ID");
    gals.AddFloatColumn("ut_LSST");
    gals.AddFloatColumn("gt_LSST");
    gals.AddFloatColumn("rt_LSST");
    gals.AddFloatColumn("it_LSST");
    gals.AddFloatColumn("zt_LSST");
    gals.AddFloatColumn("yt_LSST");
    gals.AddFloatColumn("uo_LSST");
    gals.AddFloatColumn("go_LSST");
    gals.AddFloatColumn("ro_LSST");
    gals.AddFloatColumn("io_LSST");
    gals.AddFloatColumn("zo_LSST");
    gals.AddFloatColumn("yo_LSST");
    gals.AddFloatColumn("ERRu_LSST");
    gals.AddFloatColumn("ERRg_LSST");
    gals.AddFloatColumn("ERRr_LSST");
    gals.AddFloatColumn("ERRi_LSST");
    gals.AddFloatColumn("ERRz_LSST");
    gals.AddFloatColumn("ERRy_LSST");
    gals.AddFloatColumn("U_ABS_SDSS");
    gals.AddFloatColumn("G_ABS_SDSS");
    gals.AddFloatColumn("R_ABS_SDSS");
    gals.AddFloatColumn("I_ABS_SDSS");
    gals.AddFloatColumn("Z_ABS_SDSS");
    DataTableRow rowin = gals.EmptyRow();
    
    
    outfile = outroot + "_forbpz.txt"; 
    ofstream outp(outfile.c_str(), ofstream::out);


    // START LOOP

    // time how long
    Timer tm("timer",false);
    double maxmi=0;
    
	// loop over whole catalog
	int cnt = 0, cntlost=0;
    for (int i=0; i<ng; i++) {

        // print statement to see how far code has got
        if (i>0 && i%10000 == 0) {
            tm.Split();
            
            cout << "     On galaxy: " << i+1 <<" of "<< ng <<", took "<< tm.PartialElapsedTime() << "s ";
            cout << "to do 10,000 galaxies" << endl;
            }
            
        dt.GetRow(i,row);
        
        // data from file
        double z = row[0];
        double bf = row[1];
        double utrue = row[2];
        double gtrue = row[3];
        double rtrue = row[4];
        double itrue = row[5];
        double ztrue = row[6];
        double ytrue = row[7];
        double U = row[8];
        double G = row[9];
        double R = row[10];
        double I = row[11];
        double Z = row[12];
        
        // magnitudes
        vector<double> mags;
        mags.push_back(utrue);
        mags.push_back(gtrue);
        mags.push_back(rtrue);
        mags.push_back(itrue);
        mags.push_back(ztrue);
        mags.push_back(ytrue);
        

        // correct magnitudes to match BPZ fluxes
        vector<double> magsCorrected = bpzCorrectMags.correct(mags, sedFiles[int(bf)], z);

 
        // observed magnitudes
        vector<double> uObservation = simData.addLSSTError(magsCorrected[0], nYear, 0);
		vector<double> gObservation = simData.addLSSTError(magsCorrected[1], nYear, 1);
		vector<double> rObservation = simData.addLSSTError(magsCorrected[2], nYear, 2);
		vector<double> iObservation = simData.addLSSTError(magsCorrected[3], nYear, 3);
		vector<double> zObservation = simData.addLSSTError(magsCorrected[4], nYear, 4);
        vector<double> yObservation = simData.addLSSTError(magsCorrected[5], nYear, 5);
        
        
        // Make file ready for BPZ
        if (iObservation[0]<25.3) {
        
            // undetected mags set to m=99, e=1-sig
            // mags
            outp << uObservation[0] <<"  "<< gObservation[0] <<"  "<< rObservation[0] <<"  ";
            outp << iObservation[0] <<"  "<< zObservation[0] <<"  "<< yObservation[0] <<"  ";
            // mag errors
            outp << uObservation[1] <<"  "<< gObservation[1] <<"  "<< rObservation[1] <<"  ";
            outp << iObservation[1] <<"  "<< zObservation[1] <<"  "<< yObservation[1] <<"  ";
            // z
            outp << z <<"  "<< bf <<"  "<< U <<"  "<< G <<"  "<< R <<"  "<< I <<"  "<< Z <<endl;
            }
            
        
        // to write to FITS file
        rowin[0] = z;
        rowin[1] = float(bf);
        rowin[2] = magsCorrected[0];
        rowin[3] = magsCorrected[1];
        rowin[4] = magsCorrected[2];
        rowin[5] = magsCorrected[3];
        rowin[6] = magsCorrected[4];
        rowin[7] = magsCorrected[5];
        rowin[8] = uObservation[0];
        rowin[9] = gObservation[0];
        rowin[10] = rObservation[0];
        rowin[11] = iObservation[0];
        rowin[12] = zObservation[0];
        rowin[13] = yObservation[0];
        rowin[14] = uObservation[1];
        rowin[15] = gObservation[1];
        rowin[16] = rObservation[1];
        rowin[17] = iObservation[1];
        rowin[18] = zObservation[1];
        rowin[19] = yObservation[0];
        rowin[20] = U;
        rowin[21] = G;
        rowin[22] = R;
        rowin[23] = I;
        rowin[24] = Z;
        gals.AddRow(rowin);

        
        cnt++;

        }
    cout << endl;   
    outp.close();

	  
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " produceCatalogForBPZ.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " produceCatalogForBPZ.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " produceCatalogForBPZ.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of produceCatalogForBPZ.cc program  Rc= " << rc << endl;
  return rc;	
}
