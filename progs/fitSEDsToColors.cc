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
    cout <<"  a) catalog containing: LSST ugrizy (true), SED id, z^, SDSS UGRIZ^                       "<<endl;
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


    // this is a file for checking
    outfile = outroot + "_fittedSEDs.txt";
    ofstream outp(outfile.c_str(), ofstream::out);
    
    
    // Open FITS file for writing (file a)
    outfile = outroot + "_lsstcat.fits";
    FitsInOutFile swf(outfile, FitsInOutFile::Fits_Create);
	SwFitsDataTable gals(swf, 2048);
    gals.AddFloatColumn("Z_TRUE");
    gals.AddFloatColumn("SED_ID");
    gals.AddFloatColumn("ut_LSST");
    gals.AddFloatColumn("gt_LSST");
    gals.AddFloatColumn("rt_LSST");
    gals.AddFloatColumn("it_LSST");
    gals.AddFloatColumn("zt_LSST");
    gals.AddFloatColumn("yt_LSST");
    gals.AddFloatColumn("U_ABS_SDSS");
    gals.AddFloatColumn("G_ABS_SDSS");
    gals.AddFloatColumn("R_ABS_SDSS");
    gals.AddFloatColumn("I_ABS_SDSS");
    gals.AddFloatColumn("Z_ABS_SDSS");
    DataTableRow rowin = gals.EmptyRow();
    

    // timer
    Timer tm("timer",false);
    double maxmi=0;
    
	// loop over whole catalog
	int cnt = 0;
    for (int i=0; i<ng; i++) {

        // print statement to see how far code has got
        if (i>0 && i%10000 == 0) {
            tm.Split();
            
            cout << "     On galaxy: " << i+1 <<" of "<< ng <<", took "<< tm.PartialElapsedTime() << "s ";
            cout << "to do 10,000 galaxies" << endl;
            }
            
        dt.GetRow(i,row);

        
        // data from file
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
        

        // create vector of colors
        vector<double> colors;
        colors.push_back(U-G);
        colors.push_back(G-R);
        colors.push_back(R-I);
        colors.push_back(I-Z);
        

        // best fit SED
        int bf = sedLibColors.bestSED(colors); // this is zero indexed
        // index corresponds to index of SED in sedArray (zero=first entry)
        
        
        // for checking
        for (int j=0; j<colors.size(); j++)
            outp << colors[j] <<"  ";
        outp << bf+1 << endl;
        
        // next part, generate photometry
        vector<double> mags;
        bool all_detected = true;
        for (int j=0; j<lsst_filters.size(); j++) {
        
            double mag = simData.GetMag(z, bf, R, j, (*sdss_filters[2]));
            mags.push_back(mag);
            if (mag>50.)
                all_detected = false;
            }
            
        
        // to write to FITS file
        rowin[0] = z;
        rowin[1] = float(bf);
        rowin[2] = mags[0];
        rowin[3] = mags[1];
        rowin[4] = mags[2];
        rowin[5] = mags[3];
        rowin[6] = mags[4];
        rowin[7] = mags[5];
        rowin[8] = U;
        rowin[9] = G;
        rowin[10] = R;
        rowin[11] = I;
        rowin[12] = Z;
        gals.AddRow(rowin);
        
        
        cnt++;
        //if (cnt>1000000)
        //    break;
        }
    cout << endl;
        
    /*cout <<"     Number of lost galaxies = "<< cntlost << endl;
    cout << endl;

    cout <<"     Total time taken to fit SED: "<< float(sed_time)/1000. <<"s"<< endl;
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
    
    //cout <<"     Max SDSS i for observed LSST i<25.3 was "<< maxmi << endl;
    
    outp.close();
    
    /*for (int i=0; i<debug_out.size(); i++)
		(*debug_out[i]).close();*/
	  
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
