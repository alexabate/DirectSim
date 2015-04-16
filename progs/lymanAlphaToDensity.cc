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
//#include "genericfunc.h"
#include "poly.h"
#include "mydefrg.h"
#include "geneutils.h"
#include "sedfilter.h"
#include "simdata.h"


void usage(void);
void usage(void) {
	cout << endl<<" Usage: lymanAlphaToDensity [...options...]" << endl<<endl;
	cout << " -o: OUTFILEROOT: write files to filename beginning OUTROOT (saved to output/)"<<endl;
	cout << " -z: ZSOURCE: redshift of source quasar(s) "<<endl;
	cout << " -n: NLINES: number of lines of sight to quasars "<<endl;
	cout << " -w: NL: wavelength resolution "<<endl;
	cout << endl;
    };

int main(int narg, char* arg[]) {

    cout << " ==== lymanAlphaToDensity.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;
	
	string outfileroot = "output/lymanAlphaToDensity";
    double zSource = 3.5;
    int nLines = 20;
    int nl = 3500; // roughly SDSS number of pixels per spectrum
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
    while((c = getopt(narg,arg,"ho:z:n:w:")) != -1) {
	    switch (c)  {
	        case 'o' :
	            outfileroot = optarg;
	            outfileroot = "output/" + outfileroot;
	            break;
	        case 'z' :
	            sscanf(optarg,"%lf",&zSource);
	            break;
	        case 'n' :
	            sscanf(optarg,"%d",&nLines);
	            break;
	        case 'w' :
	            sscanf(optarg,"%d",&nl);
	            break;
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }

    //-- end command line arguments
    cout <<"     Writing to files beginning "<< outfileroot <<endl;
    cout <<"     Redshift of source = " << zSource << endl;
    cout <<"     Number of lines of sight = "<< nLines << endl;
    cout <<"     Wavelength resolution = "<< nl << endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    string infile, outfile;
    ifstream inp;
	ofstream outp, outp2;
	
	// Calculates g(N_HI) given in eqn 4 Inoue & Iwata 2008
	cout <<"     Calculating absorber HI column density distribution"<<endl;
    //double beta1=1.6, beta2=1.3;
    //int nStep = 10000;
    //double nHImin = 1e12, nHImax = 1e17; // want LAF's only
    HIColumnDensity hiColumnDensity;//(beta1,beta2,1.6e17,nHImin,nHImax,nStep);
    double Nl, Nu;
    hiColumnDensity.returnColDensityLimits(Nl,Nu);
    //double intVal = hiColumnDensity.checkIntegration(nStep);
    int nStep=500000;
    //double logNl=log(Nl);
    //double dlogN=(log(Nu)-logNl)/(nStep-1);
    cout << endl;


	// Calculates f(z) given in eqn 5 in Inoue & Iwata 2008
	cout <<"     Calculating absorber redshift distribution"<<endl;
    //double gamma1=0.2, gamma2=2.5, gamma3=4.0;
    AbsorberRedshiftDistribution absorberZDist;//(gamma1,gamma2,gamma3);
    //double zmin=0.2, zmax=6;
    //double dz=(zmax-zmin)/(nStep-1);
    cout << endl;
	
	
	// Calculates h(b) given in eqn 6 in Inoue & Iwata 2008
	cout <<"     Calculating absorber b-parameter distribution"<<endl;
    //double bsig=23;
    DopplerParDistribution dopplerParDist;//(bsig);
    //double bmin=2, bmax=200;
    //double db=(bmax-bmin)/(nStep-1);
    cout << endl;
	
	
    // Generate Monte Carlo distributions
    cout <<"     Generate Monte Carlo distribution of absorbers "<<endl;
    RandomGenerator rg;
	ProbabilityDistAbsorbers probDistAbsorbers(rg, absorberZDist, 
	                                           hiColumnDensity, dopplerParDist);
	                                           
	// Make CDF of Gaussian, arbitrary normalization
    double mean = 0.5, sig = 0.1;
    double xmin = 0., xmax = 1.;
    int nx = 201;
    double dx = (xmax - xmin)/(nx - 1);
    double sumGaussian = 0;
    vector<double> xvals, cdfG;
    //outfile = outfileroot + "_cdfGaussian.txt";
    //outp.open(outfile.c_str(), ofstream::out);
    for (int i=0; i<nx; i++) {
        double xval = xmin + i*dx;
        double gaussian = 1/(sqrt(2*PI*sig*sig))*exp(-(xval-mean)*(xval-mean)/(2*sig*sig));
        sumGaussian += gaussian;
        
        xvals.push_back(xval);
        cdfG.push_back(sumGaussian*dx);
        //outp << xval <<"  "<< cdfG[i] << endl;
        }
    //outp.close();
    // look up table: give it CDF value, return x value
    //cout << xvals[166] <<"  "<< cdfG[166] <<"  "<< xvals[167] <<"  " <<cdfG[167] << endl;
    SInterp1D cdfGaussian(cdfG, xvals, cdfG[0], cdfG[nx-1]);
    
    int nLya = 2; // lyman-alpha only
    
    // observed frame wavelengths in meters
    // roughly SDSS range and resolution
    double lmin = 3700e-10, lmax = 10400e-10;
	double dl = (lmax - lmin)/(nl - 1);
	                                                
	double zStart = 0, zMax =6;                                       
	// Now simulate the lines of sight
	vector<int> nAbsorbers;
	vector<double> totTrans;
	for (int i=0; i<nLines; i++) {
	    cout <<"     Simulating line of sight "<< i+1 <<" of "<< nLines <<endl;
	    
	    stringstream ss;
	    ss << i+1;
	    vector<double> redshifts, dopplers, columnDensities;
	    outfile = outfileroot + "_lineofsight"+ss.str()+".txt";
	    //int nAbs = 
	    probDistAbsorbers.simulateLineOfSight(zStart,zMax,redshifts, 
	               dopplers, columnDensities, outfile);
	    int nAbs = redshifts.size();

        // Lyman-alpha only
        LineOfSightTrans lightOfSightTransLya(redshifts, dopplers, columnDensities);
        lightOfSightTransLya.setMaxLine(nLya);
        lightOfSightTransLya.setLymanSeriesOnly();


	    outfile = outfileroot + "_lineofsightTrans"+ss.str()+".txt";
		outp.open(outfile.c_str(), ofstream::out);
		outp <<"# number of absorbers = "<< nAbs << endl;
		nAbsorbers.push_back(nAbs);
		
		
		// loop over wavelengths
		vector<double> lams, trans;
	    for (int j=0; j<nl; j++) {
	        //cout <<"     On wavelength "<<j+1<<" of "<<nl<<endl;
	        double lam = lmin + j*dl;
	            
	        double transLa = lightOfSightTransLya(lam, zSource);
	        
	        lams.push_back(lam);
	        trans.push_back(transLa);
	        
	        if (transLa<0. || transLa>1.)
	            cout <<"ERROR! transmission flux makes no sense!" << endl;

	        outp << lam <<"  "<< transLa <<endl;
	        }
        outp.close();
        
        // make flux PDF
        // create binning in order to make pdf
        double fmin = 0., fmax = 1.; // min and max values of transmitted flux
        int nf = 100;
        double df = (fmax - fmin)/(nf - 1);
        double pdf[nf];
        vector<double> fvals;
        // initialize
        for (int j=0; j<nf; j++) {
            fvals.push_back(fmin + j*df);
            pdf[j]=0.;
            }
        // fill bins, loop over all transmitted flux values
        for (int j=0; j<nl; j++) {
            double fval = trans[j];
            int ibin = (int)floor((fval-fmin)/df);
            pdf[ibin]++;
            }
        // normalize so that int pdf df = 1
        /*double sumPDF = 0.;
        for (int j=0; j<nf; j++) {
            sumPDF += pdf[j];
            }
        cout << "summing all pdf values was: "<< sumPDF << endl;*/
        //outfile = outfileroot + "_pdf"+ss.str()+".txt";
		//outp.open(outfile.c_str(), ofstream::out);
        for (int j=0; j<nf; j++) {
            pdf[j] /= nl*df;
            //outp << fvals[j] <<"  "<< pdf[j] << endl;
            }
        //outp.close();
        // check int pdf df = 1 and make flux cdf
        double sumPDF = 0.;
        vector<double> cdf;
        //outfile = outfileroot + "_cdf"+ss.str()+".txt";
		//outp.open(outfile.c_str(), ofstream::out);
        for (int j=0; j<nf; j++) {
            sumPDF += pdf[j];
            cdf.push_back(sumPDF*df);
            if (j<1 && cdf[j]>0)
                cdf[j] = 0.;
            //outp << fvals[j] <<"  "<< cdf[j] << endl;
            }
        //outp.close();
        //cout <<"int pdf df should equal 1 ====> "<< sumPDF*df <<endl;
        // look up table: give it flux value return CDF value
        SInterp1D cdfFlux(fvals, cdf, fvals[0], fvals[nf-1]);
        // loop over transmitted flux values
        // open file
	    outfile = outfileroot + "_lineofsightTransDensity"+ss.str()+".txt";
        cout << "     Writing to "<< outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
		
		double sumDens = 0., sumTransFlux = 0.;
        for (int j=0; j<nl; j++) {
            double fval = trans[j];
            double cdfval = cdfFlux(fval);
            double xval = cdfGaussian(cdfval);
            double lam = lmin + j*dl;
            outp << lam <<"  "<< lam/(1+zSource) <<"  "<< fval <<"  "<< xval << endl;
            sumDens += xval;
            sumTransFlux += fval;
            }
        totTrans.push_back(sumTransFlux);
        cout <<"     Sum of density values = "<< sumDens <<endl;
        cout <<"     Sum of transmitted flux = "<< sumTransFlux <<endl;
        outp.close();
        cout << endl;
        }
    cout << endl;
                
    cout <<"     Calculating the effect on the u-g color of a galaxy at z = "<< zSource << endl;
    // Now we are going to calculate the affect on the color of a galaxy at z=zSource
    
    // Load in SEDs
    cout <<"     Load in SEDs"<<endl;
	double lminSED=100e-10, lmaxSED=8000e-10;
	string sedFile = "CWWKSB.list";
	ReadSedList readSedList(sedFile);
    readSedList.readSeds(lminSED,lmaxSED); // Read out SEDs into array
    vector<SED*> sedArray=readSedList.getSedArray();

        
    // Load in filters required
    // LSST
	cout <<"     Load in LSST filters"<<endl;
	string filterFile = "LSST.filters";
	ReadFilterList readFilterList(filterFile);
	readFilterList.readFilters(lminSED,lmaxSED);
	vector<Filter*> LSSTfilters=readFilterList.getFilterArray();
	// Reference
	//cout <<"     Load in reference (GOODS B) filter"<<endl;
	string goodsFilterFile = "GOODSB.filters";
	ReadFilterList readGOODSBfilter(goodsFilterFile);
	readGOODSBfilter.readFilters(lminSED,lmaxSED);
	vector<Filter*> goodsBFilter=readGOODSBfilter.getFilterArray();


	// Set cosmology
	cout <<"     Set cosmology"<<endl;
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout <<"     Set cosmology to: OmegaM="<<OmegaM<<", OmegaL="<<OmegaL;
	cout <<", H0="<<100*h<<endl;
	cout << endl;
    
    int nElliptical = 1;
    int nSpiral = 2;
	SimData simgal(sedArray, LSSTfilters, su); //,rg,nElliptical,nSpiral);
	
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
        infile = outfileroot + "_lineofsightTrans"+ss.str()+".txt";
        string line;
        inp.open(infile.c_str());
        getline(inp,line);
        /*string nAbsString = line[line.size()-2] ;
        if (line.size()>25)
            nAbsString = line[line.size()-2] + line[line.size()-1];
        else
            nAbsString = line[line.size()-1];
        int nAbs = atoi(nAbsString.c_str());
        cout <<"     Number of absorbers = "<<nAbs <<endl;*/

        IGMTransmission igmTransmission(infile, lmin, lmax, nl);
        //SEDIGM sedIGM(*(sedArray[sedNo]), igmTransmission, zs);
        //SEDIGM sedMadau(*(sedArray[sedNo]), zs);
        
		//double uMag=simgal.GetMag(zs,type,am,ext,0,(*goodsBFilter[0]),igmTransmission);
		//double gMag=simgal.GetMag(zs,type,am,ext,1,(*goodsBFilter[0]),igmTransmission);
		double uMag = simgal.getMag(zs, am, sedNo, 0, (*goodsBFilter[0]), igmTransmission);
		double gMag = simgal.getMag(zs, am, sedNo, 1, (*goodsBFilter[0]), igmTransmission);
		
		outp << nAbsorbers[i] <<"  "<< uMag <<"  "<< gMag <<"  "<< totTrans[i] << endl;
		
		inp.close();
        }
    outp.close();
    cout << endl;
    
    // Now calculate u and g mags with NO igm
    cout <<"     Calculate u and g magnitudes w/Madau & w/oIGM for the starburst galaxy as a function of z";
    double zmin = 0.01, zmax = zSource+3.;
    int nz = 10000;
    double dz = (zmax - zmin)/(nz - 1);
    cout << zmin <<"<z<"<< zmax <<endl;
    
    //SimData simgalNoVaryIGM(sedArray,LSSTfilters,su,rg,nElliptical,nSpiral);
    IGMTransmission noIGM(lmin, lmax, nl);
    
    outfile = outfileroot + "_magsZ.txt";
	outp.open(outfile.c_str());
    for (int i=0; i<nz; i++){
    
        double zs = zmin + i*dz;
        double am = -18;
        double ext = 0.;
        double type = 3.007;
        
        /*double uMag = simgalNoVaryIGM.GetMag(zs,type,am,ext,0,(*goodsBFilter[0]));
		double gMag = simgalNoVaryIGM.GetMag(zs,type,am,ext,1,(*goodsBFilter[0]));
		simgalNoVaryIGM.setMadau(false);
		double uMagNoIGM = simgalNoVaryIGM.GetMag(zs,type,am,ext,0,(*goodsBFilter[0]));
		double gMagNoIGM = simgalNoVaryIGM.GetMag(zs,type,am,ext,1,(*goodsBFilter[0]));
		simgalNoVaryIGM.setMadau(true);*/
		
		int sedid = 7;
		
		// No IGM
		double uMag = simgal.getMag(zs, am, sedid, 0, (*goodsBFilter[0]), noIGM);
		double gMag = simgal.getMag(zs, am, sedid, 1, (*goodsBFilter[0]), noIGM);
		
		// with Madau
		IGMTransmission madau(zs);
		double uMagMadau = simgal.getMag(zs, am, sedid, 0, (*goodsBFilter[0]), madau);
		double gMagMadau = simgal.getMag(zs, am, sedid, 1, (*goodsBFilter[0]), madau);
		
        outp << zs <<"  "<< uMag <<"  "<< gMag <<"  "<< uMagMadau <<"  "<< gMagMadau <<endl;
        }
    outp.close();
        
	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " lymanAlphaToDensity.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " lymanAlphaToDensity.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " lymanAlphaToDensity.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of lymanAlphaToDensity.cc program  Rc= " << rc << endl;
  return rc;	
}

