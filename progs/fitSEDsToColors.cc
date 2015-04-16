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
	
    cout <<"  Reads absolute magnitudes from a FITS file and finds the best matched SED in the supplied"<<endl;
    cout <<"  SED library to these rest-frame colors.                                                  "<<endl;
    cout << endl;
    
    cout <<"  Then ->                                                                                  "<<endl;
    cout <<"  Using the redshift from the FITS file, an absolute magnitude, and this best fit SED,     "<<endl;
    cout <<"  observed LSST ugrizy magnitudes *without errors* are simulated.                          "<<endl;
    cout << endl;
    
    cout <<"  IMPORTANT NOTE: a list of column names corresponding to the absolute magnitude data in   "<<endl;
    cout <<"  the FITS file must be supplied as an argument AND this list must match *exactly* the     "<<endl;
    cout <<"  contents and order in the list of filters these magnitudes are in (also supplied as an   "<<endl;
    cout <<"  argument, -f).                                                                           "<< endl;
    cout << endl;
    
    cout <<"  Files that are are output:                                                               "<<endl;
    cout <<"  MAIN OUTPUT: FITS file catalog containing: LSST ugrizy (true), SED id, z^, [ABS MAGS]^   "<<endl;
    cout <<"     ^ denotes values copied direct from input FITS file                                   "<<endl;
    cout << endl;
    cout <<"  SUPPLEMENTARY OUTPUT:                                                                    "<<endl;
	cout <<"  - the SEDs that were generated to calculate the photometry                               "<<endl;
	cout <<"  - a dictionary for matching SED in above file to the E(B-V) that was added to it (if any)"<<endl;
	cout <<"  - color of each SED in library (row order matching SED library file list; column order   "<<endl;
	cout <<"    matching filter file list)                                                             "<<endl;
	cout <<"  - all colors from input FITS file (column order matching filter file list) with last     "<<endl;
	cout <<"    column being row position of best match SED in library list                            "<<endl;
	cout << endl;
	
	cout <<"  Example usage:                                                                           "<<endl;
	cout <<"  fitSEDsToColors -i tao_xxx.fits -o sim_from_file -t CWWKSB.list -n 5 -f SDSS.filters     "<<endl;
	cout <<"         -c SDSS_u_Absolute,SDSS_g_Absolute,SDSS_r_Absolute,SDSS_i_Absolute,SDSS_z_Absolute"<<endl;
	cout <<"         -m SDSS_r_Absolute,2 -w 0.0000001,0.00002,10000                                   "<<endl;
	cout << endl;
	
	cout << " -i: FITS: input FITS file catalog                                                        "<<endl; 
	cout << " -o: OUTROOT: write files to filenames beginning OUTROOT                                  "<<endl;
	cout << " -t: SEDLIB: file containing list of SED files in library                                 "<<endl;
	cout << " -f: FILTLIB: file containing list of filters                                             "<<endl;
	cout << " -c: MCOL: column names matching order of filters in FILTLIB                              "<<endl;
	cout << " -n: NRED: number of reddenings to apply to each SED in SEDLIB                            "<<endl;
	cout << " -m: ABSMCOL,IFILT: column name containing absolute magnitude to use to define galaxy     "<<endl;
	cout << "                    luminosity, id of corresponding filter in filter list file            "<<endl;
	cout << " -w: LMIN,LMAX,NL: wavelength resolution of the SEDs                                      "<<endl;
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
    string sedfile = "CWWKSB.list";
    string filtfile = "SDSS.list";
    //string abscol = "SDSS_r_Absolute";
    //int iabscol = 2; // zero indexed
    string abs = "SDSS_r_Absolute,2";
    string magcols = "SDSS_u_Absolute,SDSS_g_Absolute,SDSS_r_Absolute,SDSS_i_Absolute,SDSS_z_Absolute"; 
    int nred = 5;
    bool doRedden = true;
    double lmin=1e-8, lmax=2.5e-6; // need to be careful with this given what filters are read in!
	int npt = 1000; // interpolation resolution of SEDs (need more if many fine features)

	char c;
    while((c = getopt(narg,arg,"ho:i:t:n:f:m:c:w:")) != -1) {
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
	        case 'f' :
	            filtfile = optarg;
	            break;
	        case 'c' :
	            magcols = optarg;
	            break;
	        case 'm' :
	            abs = optarg;
	            break;
	        case 'n' :
	            sscanf(optarg,"%d",&nred);
	            break;
	        case 'w' :
	            sscanf(optarg,"%lf,%lf,%d",&lmin,&lmax,&npt);
	            break;
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }
	if (nred<1)
	    doRedden = false;
	    
	vector<string> listOfCols;
	stringSplit(magcols,",",listOfCols);
	
	vector<string> tmp;
	stringSplit(abs,",",tmp);
	string abscol = tmp[0];
	string x = tmp[1];
	int iabscol = atoi(x.c_str());

    //-- end command line arguments
    cout <<"     Writing to files beginning "<< outroot <<endl;
    cout <<"     Reading from FITS file " << infile << endl;
    cout <<"     Using SED library "<< sedfile <<" with wavelength range (in m) "<< lmin <<" to "<< lmax;
    cout <<" and resolution "<< npt << endl;
    cout <<"     Number of times to redden each template = "<< nred << endl;
    cout <<"     Fitting to photometry in filters listed in "<< filtfile <<", and reading photometry from columns:"<<endl;
    for (int i=0; i<listOfCols.size(); i++)
        cout <<"     "<< listOfCols[i] << endl;
    cout <<"     Setting luminosity of galaxies using absolute magnitude in column "<< abscol <<", which";
    cout <<" corresponds to filter "<< iabscol <<" in "<< filtfile << endl;
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
    
    // GET COLUMN INDICES
    int idMag = dt.IndexNom(abscol);
    int idRa = dt.IndexNom("Right_Ascension");
    int idDec = dt.IndexNom("Declination");
    int idz = dt.IndexNom("Redshift_Cosmological");
    int idmapp = dt.IndexNom("SDSS_i_Apparent");
    cout << "     Column: Right_Ascension has index "<< idRa << endl;
    cout << "     Column: Declination has index "<< idDec << endl;
    cout << "     Column: Redshift_Cosmological has index "<< idz << endl;
    cout << "     Column: SDSS_i_Apparent has index "<< idmapp << endl;
    cout << "     Column: "<< abscol << " has index "<< idMag << endl;
    vector<int> idCols;
    for (int i=0; i<listOfCols.size(); i++) {
        idCols.push_back(dt.IndexNom(listOfCols[i]));
        cout << "     Column: "<< listOfCols[i] << " has index "<< idCols[i] << endl;
        }
    cout << endl;
    
		
	// GALAXY SED TEMPLATES
	ReadSedList readSedList(sedfile);
	
	// Read out SEDs into array
    readSedList.readSeds(lmin, lmax, npt);
    
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
    readSedList.writeSpectra(outfile, lmin, lmax, npt);
    
    
    // FILTERS
    
    // LSST
    string filterFile = "LSST.filters";
	ReadFilterList lsstFilters(filterFile);
	lsstFilters.readFilters(lmin, lmax);
	vector<Filter*> lsst_filters = lsstFilters.getFilterArray();
	int nLSST = lsstFilters.getNTot();
	cout <<"     "<< nLSST <<" LSST filters read in "<<endl;
    cout << endl;
	
	// PHOTOMETRY FILTERS
	ReadFilterList photFilters(filtfile);
	photFilters.readFilters(lmin, lmax);
	vector<Filter*> phot_filters = photFilters.getFilterArray();
	int nphotFilt = photFilters.getNTot();
	cout <<"     "<< nphotFilt <<" photometry filters read in "<<endl;
	if (nphotFilt != idCols.size())
	    throw ParmError("ERROR! number of photometry columns != number of photometry filters");
    cout << endl;
	
	
	// FIT SEDS AND GENERATE PHOT
	
	// for fitting SEDs to colors
	SEDLibColors sedLibColors(sedArray, phot_filters, lmin, lmax, npt);
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
	SimData simData(sedArray, lsst_filters, su); //, rg);
	//SimObservations simObs(lsst_filters, rg);

    // Add Madau preference
    //bool isAddMadau = true;
    //simData.setMadau(isAddMadau);
    

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
    for (int i=0; i<idCols.size(); i++)
        gals.AddFloatColumn(listOfCols[i]);
    DataTableRow rowin = gals.EmptyRow();
    

    // timer
    Timer tm("timer",false);
    double maxmi=0;
    
    Timer loop_time("lptm", false);
    int lt=0;
    
    Timer bestfit_time("bftm", false);
    int bft=0;
    
    Timer phot_time("phtm", false);
    int pt=0;
    
    
	// loop over whole catalog
	cout <<"     Starting loop over "<< ng <<" galaxies in the catalog"<<endl;
	int cnt = 0;
    for (int i=0; i<ng; i++) {

        // print statement to see how far code has got
        if (i>0 && i%10000 == 0) {
            tm.Split();
            
            cout << "     On galaxy: " << i+1 <<" of "<< ng <<", took "<< tm.PartialElapsedTime() << "s ";
            cout << "to do 10,000 galaxies" << endl;
            }
            
        dt.GetRow(i,row);
        //cout <<"row:"<< row << endl;

        
        // data from file
        double ra = row[idRa];
        double dec = row[idDec];
        double z = row[idz];
        double mi = row[idmapp]; // don't actually need this
        double R = row[idMag];
        //cout <<"dats:"<< ra <<" "<< dec <<" "<< z <<" "<< mi <<" "<< R << endl;
        
        // IGM model (Madau)
        IGMTransmission igm(z);

        
        // create vector of colors
        vector<double> colors;
        for (int i=0; i<idCols.size()-1; i++) {
            double c1 = row[idCols[i]];
            double c2 = row[idCols[i+1]];
            //cout << "col:"<< c1-c2 << endl;
            colors.push_back(c1-c2);
            }
        

        // best fit SED
        bestfit_time.Split();
        int bf = sedLibColors.bestSED(colors); // this is zero indexed
        bestfit_time.Split();
        bft+=bestfit_time.PartialElapsedTimems();
        
        //cout << "bf: "<< bf << endl;
        // index corresponds to index of SED in sedArray (zero=first entry)
        
        
        // for checking
        for (int j=0; j<colors.size(); j++)
            outp << colors[j] <<"  ";
        outp << bf+1 << endl;
        
        // next part, generate photometry
        phot_time.Split();
        vector<double> mags;
        bool all_detected = true;
        for (int j=0; j<lsst_filters.size(); j++) {
        
            //double mag = simData.GetMag(z, bf, R, j, (*phot_filters[iabscol]));
            
            double mag = simData.getMag(z, R, bf, j, (*phot_filters[iabscol]), igm);
            mags.push_back(mag);
            if (mag>50.)
                all_detected = false;
            }
        phot_time.Split();
        pt+=phot_time.PartialElapsedTimems();
            
        
        // to write to FITS file
        rowin[0] = z;
        rowin[1] = float(bf);
        rowin[2] = mags[0];
        rowin[3] = mags[1];
        rowin[4] = mags[2];
        rowin[5] = mags[3];
        rowin[6] = mags[4];
        rowin[7] = mags[5];
        int ii = 8;
        for (int i=0; i<idCols.size(); i++) {
            rowin[ii] = row[idCols[i]];
            ii++;
            }
        gals.AddRow(rowin);
        
        
        cnt++;
        //if (cnt>1000000)
        //    break;
        loop_time.Split();
        lt+=loop_time.PartialElapsedTimems();
        }
    cout << endl;
    cout << "Average time per loop = "<< double(lt)/ng <<" ms "<<endl;
    cout << "Average time per best-fit SED calculation = "<< double(bft)/ng <<" ms "<<endl; 
    cout << "Average time per photometry calculation = "<< double(pt)/ng <<" ms "<<endl;
    
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
