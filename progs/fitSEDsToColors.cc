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
	cout << endl<<" Usage: fitSEDsToColors [...options...]" << endl<<endl;
	
    cout <<"  Reads absolute SDSS UGRIZ magnitudes from a FITS file and finds the best matched SED in  "<<endl;
    cout <<"  the supplied SED library to these rest-frame colors.                                     "<<endl;
    cout << endl;
    
    cout <<"  Using the redshift from the FITS file, SDSS absolute magnitude and the best fit SED,     "<<endl;
    cout <<"  observed LSST ugrizy magnitudes (with errors) are simulated.                             "<<endl;
    cout << endl;
    
    cout <<"  Three files are output:                                                                  "<<endl;
    cout <<"  a) catalog containing: LSST ugrizy+errors, SED id, z^, R^ mag, (U-G)^ color, (G-R)^ color"<<endl;
    cout <<"     ^ denotes values copied direct from FITS file                                         "<<endl;
	cout <<"  b) the SEDs that were generated to calculate the photometry                              "<<endl;
	cout <<"  c) a dictionary for matching SED id to SED in file and the E(B-V) that was added to it   "<<endl;
	cout << endl;
	
	cout <<"  Example usage:                                                                           "<<endl;
	cout <<"  fitSEDsToColors -o sim_from_file -t CWWKSB.list -i tao_xxx.fits -n 5                     "<<endl;
	cout << endl;
	
	cout << " -o: OUTROOT: write files to filenames beginning OUTROOT                                  "<<endl;
	cout << " -t: SEDLIB: file containing list of SED files                                            "<<endl;
	cout << " -i: FITS: FITS file                                                                      "<<endl; 
	cout << " -n: NRED: number of reddenings to apply to each SED in SEDLIB                            "<<endl;
	cout << endl;
    };


int main(int narg, char* arg[]) {
    
    cout << " ==== fitSEDsToColors.cc program , to generate LSST mags from TAO simulation catalog ====" << endl;
    cout << endl;

    // make sure SOPHYA modules are initialized 
    SophyaInit();  
    FitsIOServerInit();
    InitTim();
    cout<<endl<<endl;

    //--- decoding command line arguments 
    string outroot;
    string infile;
    string sedfile="CWWKSB.list";
    int nred = 5;
    bool doRedden = true;

	char c;
    while((c = getopt(narg,arg,"ho:i:t:n:")) != -1) {
	    switch (c)  {
	        case 'o' :
	            outroot = optarg;
	            break;
	        case 'i' :
	            infile = optarg;
	            break;
	        case 't' :
	            sedfile = optarg;
	            break;
	        case 'n' :
	            sscanf(optarg,"%d",&nred);
	            break;
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }
	if (nred<1)
	    doRedden = false;

    //-- end command line arguments
    cout <<"     Writing to files beginning "<< outroot <<endl;
    cout <<"     Reading from FITS file " << infile << endl;
    cout <<"     Using SED library "<< sedfile <<endl;
    cout <<"     Number of times to redden each template = "<< nred << endl;
    cout << endl;
    //-- end command line arguments
  
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
	ResourceUsage res;
	InitTim();
	string outfile;
	
	
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
	
	// Read out SEDs into array
    readSedList.readSeds(lmin,lmax);
    
    if (doRedden) {
        // Redden SEDs 
        int method = 0; // uniform distribution up to max
        int maxidEl = 0;
        double redMaxEl = 0.1;
        double redMaxOther = 0.3;
        vector<double> reds = readSedList.reddenSeds(nred, method, maxidEl, redMaxEl, redMaxOther);
        }
        
    int nsedOrig = readSedList.getNSed();
    int nsedTot = readSedList.getNTot();
    cout <<"     Number of original SEDs = "<< nsedOrig <<endl;
    if (doRedden)
        cout <<"     Number of SEDs after adding reddened ones = "<< nsedTot <<endl;
	cout << endl;
	
	// Get total number of SEDs
    vector<SED*> sedArray = readSedList.getSedArray();
    cout <<"     Number of SEDs in SED array "<< sedArray.size() << endl;
    cout << endl;
    
    // write out SEDs (this generates file b and c as described in the usage instructions)
    outfile = outroot + "_SEDlib.txt";
    readSedList.writeSpectra(outfile);
    
    
    // FILTERS
    
    // LSST
    string filterFile = "LSST.filters";
	ReadFilterList lsstFilters(filterFile);
	lsstFilters.readFilters(lmin, lmax);
	vector<Filter*> lsst_filters = lsstFilters.getFilterArray();
	int nLSST = lsstFilters.getNTot();
	cout <<"     "<< nLSST <<" LSST filters read in "<<endl;
    cout << endl;
	
	// SDSS
    filterFile = "SDSS.filters";
	ReadFilterList sdssFilters(filterFile);
	sdssFilters.readFilters(lmin, lmax);
	vector<Filter*> sdss_filters = sdssFilters.getFilterArray();
	int nSDSS = sdssFilters.getNTot();
	cout <<"     "<< nSDSS <<" SDSS filters read in "<<endl;
    cout << endl;
	
	
	// FIT SEDS AND GENERATE PHOT
	
	// for fitting SEDs to colors
	SEDLibColors sedLibColors(sedArray, sdss_filters);
	// this file also for extra checking
	outfile = outroot + "_colorarray.txt";
	sedLibColors.writeColorArray(outfile);
	
	// for calculating distance modulus
    //Cosmological parameters: Om = 0.25, OL = 0.75, Ob = 0.045, s8 = 0.9, h = 0.73,
    double h=0.73, OmegaM = 0.25, OmegaL = 0.75, Ob = 0.045;
	SimpleUniverse su(h, OmegaM, OmegaL);
	su.SetOmegaBaryon(Ob);
	RandomGenerator rg;
	
	// for generating photometry
	SimData simData(sedArray, lsst_filters, su, rg);

    // Add Madau preference
    bool isAddMadau = true;
    simData.setMadau(isAddMadau);
    
    // Number of visits per year (Table 1, Ivezic et al 2008)
    int uVisitsPerYear = 6;
    int gVisitsPerYear = 8;
    int rVisitsPerYear = 18;
    int iVisitsPerYear = 18;
    int zVisitsPerYear = 16;
    int yVisitsPerYear = 16;
    
    // total number of visits
    int nYear = 10;
	int uVisits = uVisitsPerYear*nYear;
    int gVisits = gVisitsPerYear*nYear;
    int rVisits = rVisitsPerYear*nYear;
    int iVisits = iVisitsPerYear*nYear;
    int zVisits = zVisitsPerYear*nYear;
    int yVisits = yVisitsPerYear*nYear;

    // this is a file for checking
    outfile = outroot + "_fittedSEDs.txt";
    ofstream outp(outfile.c_str(), ofstream::out);
    
    // Open FITS file for writing (file a)
    outfile = outroot + "_lsstcat.fits";
    FitsInOutFile swf(outfile, FitsInOutFile::Fits_Create);
    //swf.MoveAbsToHDU(2);
	SwFitsDataTable gals(swf, 2048);
	
    gals.AddFloatColumn("Z_TRUE");
    gals.AddFloatColumn("SED_ID");
    gals.AddFloatColumn("u_LSST");
    gals.AddFloatColumn("g_LSST");
    gals.AddFloatColumn("r_LSST");
    gals.AddFloatColumn("i_LSST");
    gals.AddFloatColumn("z_LSST");
    gals.AddFloatColumn("y_LSST");
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
    
    // for testing annz
    //ofstream outp2(outfile.c_str(), ofstream::out);
    vector<ofstream*> debug_out;
    outfile = outroot + "_annz1.txt"; // zero err, all detected
    debug_out.push_back(new ofstream());
    (*debug_out[0]).open(outfile.c_str(), ofstream::out);
    outfile = outroot + "_annz2.txt"; // zero err, not all detected
    debug_out.push_back(new ofstream());
    (*debug_out[1]).open(outfile.c_str(), ofstream::out);
    outfile = outroot + "_annz3.txt"; // i<25.3
    debug_out.push_back(new ofstream());
    (*debug_out[2]).open(outfile.c_str(), ofstream::out);
    outfile = outroot + "_annz3a.txt";// i<25.3 u->99, eu-<1-sig
    debug_out.push_back(new ofstream());
    (*debug_out[3]).open(outfile.c_str(), ofstream::out);
    outfile = outroot + "_annz4.txt";
    debug_out.push_back(new ofstream());// add absolute mag to ANNz
    (*debug_out[4]).open(outfile.c_str(), ofstream::out);
    outfile = outroot + "_annz5.txt";
    debug_out.push_back(new ofstream());// add % error to flux
    (*debug_out[5]).open(outfile.c_str(), ofstream::out);
    
    // time how long
    Timer tm("timer",false);
    //Timer tmsm("smalltimes",false);
    //int sed_time=0, thphot_time=0, obsphot_time=0;
    double maxmi=0;
    
	// loop over whole catalog
	int cnt = 0, cntlost=0;
    for (int i=0; i<ng; i++) {

        if (i>0 && i%10000 == 0) {
            tm.Split();
            
            cout << "     On galaxy: " << i+1 <<" of "<< ng <<", took "<< tm.PartialElapsedTime() << "s ";
            cout << "to do 10,000 galaxies" << endl;
            }
            
        dt.GetRow(i,row);
        //cout << row << endl;
        
        int gal_class = row[0];
        double ra = row[1];
        double dec = row[2];
        double z = row[3];
        // z obs (zcosmo+vpec)
        double U = row[5];
        double G = row[6];
        double R = row[7];
        double I = row[8];
        double mi = row[9];
        double Z = row[10];
        

        vector<double> colors;
        colors.push_back(U-G);
        colors.push_back(G-R);
        colors.push_back(R-I);
        colors.push_back(I-Z);
        
        //tmsm.Split();
        int bf = sedLibColors.bestSED(colors);
        //tmsm.Split();
        //sed_time += tmsm.PartialElapsedTimems();
        
        //cout <<"     Best fit SED = "<< bf << endl; // note this is zero indexed
        
        // for checking
        for (int j=0; j<colors.size(); j++)
            outp << colors[j] <<"  ";
        outp << bf+1 << endl;
        
        // next part, generate photometry
        //tmsm.Split();
        vector<double> mags;
        bool all_detected = true;
        for (int j=0; j<lsst_filters.size(); j++) {
            double mag = simData.GetMag(z, bf, R, j, (*sdss_filters[2]));
            mags.push_back(mag);
            if (mag>50.)
                all_detected = false;
            }
            
        // make some test files for ANNz
        
        // 1) all have to be detected: TRUTH MAG, ZERO ERROR
        if (all_detected) {
            for (int j=0; j<lsst_filters.size(); j++) 
                (*debug_out[0]) << mags[j] <<"  ";
            for (int j=0; j<lsst_filters.size(); j++) 
                (*debug_out[0]) << 0.0 <<"  ";  
            (*debug_out[0]) << z << endl;
            }
        
        // 2) none have to be detected: TRUTH MAG, ZERO ERROR
        for (int j=0; j<lsst_filters.size(); j++) 
            (*debug_out[1]) << mags[j] <<"  ";
        for (int j=0; j<lsst_filters.size(); j++) 
            (*debug_out[1]) << 0.0 <<"  ";  
        (*debug_out[1]) << z << endl;
            
            
        //tmsm.Split();
        //thphot_time += tmsm.PartialElapsedTimems();
        
        //tmsm.Split();
        vector<double> uObservation = simData.addLSSTError(mags[0],uVisits, 0);
		vector<double> gObservation = simData.addLSSTError(mags[1],gVisits, 1);
		vector<double> rObservation = simData.addLSSTError(mags[2],rVisits, 2);
		vector<double> iObservation = simData.addLSSTError(mags[3],iVisits, 3);
		vector<double> zObservation = simData.addLSSTError(mags[4],zVisits, 4);
        vector<double> yObservation = simData.addLSSTError(mags[5],yVisits, 5);
        
        double pc = 0.05;
        vector<double> uObservation2 = simData.addError(mags[0],pc, 0);
		vector<double> gObservation2 = simData.addError(mags[1],pc, 1);
		vector<double> rObservation2 = simData.addError(mags[2],pc, 2);
		vector<double> iObservation2 = simData.addError(mags[3],pc, 3);
		vector<double> zObservation2 = simData.addError(mags[4],pc, 4);
        vector<double> yObservation2 = simData.addError(mags[5],pc, 5);
        //tmsm.Split();
        //obsphot_time += tmsm.PartialElapsedTimems();
        
        // 3) none have to be detected, but iobs<25.3: OBS MAG, OBS ERR
        if (iObservation[0]<25.3) {
            // mags
            (*debug_out[2]) << uObservation[0] <<"  "<< gObservation[0] <<"  "<< rObservation[0] <<"  ";
            (*debug_out[2]) << iObservation[0] <<"  "<< zObservation[0] <<"  "<< yObservation[0] <<"  ";
            // mag errors
            (*debug_out[2]) << uObservation[1] <<"  "<< gObservation[1] <<"  "<< rObservation[1] <<"  ";
            (*debug_out[2]) << iObservation[1] <<"  "<< zObservation[1] <<"  "<< yObservation[1] <<"  ";
            // z
            (*debug_out[2]) << z << endl;
            
            // 6) add percentage error
            // mags
            (*debug_out[5]) << uObservation2[0] <<"  "<< gObservation2[0] <<"  "<< rObservation2[0] <<"  ";
            (*debug_out[5]) << iObservation2[0] <<"  "<< zObservation2[0] <<"  "<< yObservation2[0] <<"  ";
            // mag errors
            (*debug_out[5]) << uObservation2[1] <<"  "<< gObservation2[1] <<"  "<< rObservation2[1] <<"  ";
            (*debug_out[5]) << iObservation2[1] <<"  "<< zObservation2[1] <<"  "<< yObservation2[1] <<"  ";
            // z
            (*debug_out[5]) << z << endl;
            
            
            // this is to protect against non detections in u band being shown as "detected"
            double uband = uObservation[0];
            double euband = uObservation[1];
            if (uObservation[0]>28.2 || uObservation[1]>3.) {
                uband = 99.;
                euband = 28.65;
                }
            
            // 4) and set undetected u, u=99 eu=1-sig
            // mags
            (*debug_out[3]) << uband <<"  "<< gObservation[0] <<"  "<< rObservation[0] <<"  ";
            (*debug_out[3]) << iObservation[0] <<"  "<< zObservation[0] <<"  "<< yObservation[0] <<"  ";
            // mag errors
            (*debug_out[3]) << euband <<"  "<< gObservation[1] <<"  "<< rObservation[1] <<"  ";
            (*debug_out[3]) << iObservation[1] <<"  "<< zObservation[1] <<"  "<< yObservation[1] <<"  ";
            // z
            (*debug_out[3]) << z << "  ";
            for (int j=0; j<lsst_filters.size(); j++) 
                (*debug_out[3]) << mags[j] <<"  ";
            (*debug_out[3]) << R << "  "<< bf << endl;
            }
        
        // 5) add absolute magnitude as ANNz input
        if (iObservation[0]<25.3) {
            // mags
            (*debug_out[4]) << uObservation[0] <<"  "<< gObservation[0] <<"  "<< rObservation[0] <<"  ";
            (*debug_out[4]) << iObservation[0] <<"  "<< zObservation[0] <<"  "<< yObservation[0] <<"  ";
            (*debug_out[4]) << R <<"  ";
            // mag errors
            (*debug_out[4]) << uObservation[1] <<"  "<< gObservation[1] <<"  "<< rObservation[1] <<"  ";
            (*debug_out[4]) << iObservation[1] <<"  "<< zObservation[1] <<"  "<< yObservation[1] <<"  ";
            (*debug_out[4]) << 0.0 <<"  ";
            // z
            (*debug_out[4]) << z << endl;
            }
        
        // to write to FITS file
        rowin[0] = z;
        rowin[1] = float(bf);
        rowin[2] = uObservation[0];
        rowin[3] = gObservation[0];
        rowin[4] = rObservation[0];
        rowin[5] = iObservation[0];
        rowin[6] = zObservation[0];
        rowin[7] = yObservation[0];
        rowin[8] = uObservation[1];
        rowin[9] = gObservation[1];
        rowin[10] = rObservation[1];
        rowin[11] = iObservation[1];
        rowin[12] = zObservation[1];
        rowin[13] = yObservation[1];
        rowin[14] = U;
        rowin[15] = G;
        rowin[16] = R;
        rowin[17] = I;
        rowin[18] = Z;
        gals.AddRow(rowin);
        
        /*if (mi>25.5 && iObservation[0]<25.3) {
             cout <<"     LOST GAL! sdss i = "<< mi <<", observed lsst i = "<< iObservation[0] << endl;
             cntlost++;
             }*/
        
        /*if (iObservation[0]>=25.3 && iObservation[0]<=26) {
             cout <<"     FAINT GAL: observed lsst i = "<< iObservation[0] <<", truth lsst i = "<< mags[3];
             cout <<", sdss i = "<< mi;
             if (mi>iObservation[0])
                 cout <<" -> SDSS MAG FAINTER THAN LSST OBS";
             cout << endl;
             if (mi>maxmi)
                 maxmi=mi;
             }*/
        
        if (iObservation[0]<25.3 && mi>maxmi)
             maxmi = mi;
        if (iObservation[0]<25.3 && mi>26.)
            cntlost++;
        
        cnt++;
        if (cnt>1000000)
            break;
        }
    cout << endl;
        
    cout <<"     Number of lost galaxies = "<< cntlost << endl;
    cout << endl;

    /*cout <<"     Total time taken to fit SED: "<< float(sed_time)/1000. <<"s"<< endl;
    cout <<"     Average time taken to fit SED: "<< (float(sed_time)/1000.)/float(cnt) <<"s"<< endl;
    cout << endl;
    
    cout <<"     Total time taken to calculate all theoretical photometry: "<< float(thphot_time)/1000. <<"s";
    cout << endl;
    cout <<"     Average time taken to calculate single theoretical photometry: ";
    cout << (float(thphot_time)/1000.)/float(cnt) <<"s"<< endl;
    cout << endl;
    
    cout <<"     Total time taken to calculate all observed photometry: "<< float(obsphot_time)/1000. <<"s";
    cout << endl;
    cout <<"     Average time taken to calculate single observed photometry: ";
    cout << (float(obsphot_time)/1000.)/float(cnt) <<"s"<< endl;
    cout << endl;*/
    
    cout <<"     Max SDSS i for observed LSST i<25.3 was "<< maxmi << endl;
    
    outp.close();
    
    for (int i=0; i<debug_out.size(); i++)
		(*debug_out[i]).close();
	  
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " fitSEDsToColors.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " fitSEDsToColors.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " fitSEDsToColors.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of fitSEDsToColors.cc program  Rc= " << rc << endl;
  return rc;	
}
