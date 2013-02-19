// -*- LSST-C++ -*-
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
#include "genericfunc.h"
#include "poly.h"
#include "mydefrg.h"
#include "geneutils.h"


void usage(void);
void usage(void) {
	cout << endl<<" Usage: testsimulateIGM [...options...]" << endl<<endl;
	cout << " -o: OUTFILEROOT: write files to filename beginning OUTROOT"<<endl;
	cout << endl;
};

int main(int narg, char* arg[]) {

	cout << " ==== testsimulateIGM.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;
	
	string outfileroot = "testfiles/testsimulateIGM";
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
    while((c = getopt(narg,arg,"ho:")) != -1) {
	    switch (c)  {
	        case 'o' :
	            outfileroot = optarg;
	            break;
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }

    //-- end command line arguments
    cout <<"     Writing to files beginning "<< outfileroot <<endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    string outfile;
    ifstream inp;
	ofstream outp;
	
	// Write HI column density distribution function to a file.
	// Calculates g(N_HI) given in eqn 4 Inoue & Iwata 2008
	cout <<"     Calculating absorber HI column density distribution"<<endl;
    double beta1=1.6, beta2=1.3;
    int nStep = 10000;
    HIColumnDensity hiColumnDensity(beta1,beta2,1.6e17,1e12,1e22, nStep);
    double Nl, Nu;
    hiColumnDensity.returnLowerUpperColDensValues(Nl,Nu);
    //double intVal = hiColumnDensity.checkIntegration(nStep);
    nStep=500000;
    double logNl=log(Nl);
    double dlogN=(log(Nu)-logNl)/(nStep-1);
    
    cout << "     dlogN = "<<dlogN << endl;
    outfile = outfileroot + "_gNHI.txt";
    hiColumnDensity.writeToFile(outfile,dlogN,nStep);
    cout << endl;

		
	// Write absorber redshift distribution to a file.
	// Calculates f(z) given in eqn 5 in Inoue & Iwata 2008
	cout <<"     Calculating absorber redshift distribution"<<endl;
    double gamma1=0.2, gamma2=2.5, gamma3=4.0;
    AbsorberRedshiftDistribution absorberZDist(gamma1,gamma2,gamma3);
    double zmin=0.2, zmax=6;
    double dz=(zmax-zmin)/(nStep-1);
    outfile= outfileroot + "_fz.txt";
    absorberZDist.writeToFile(outfile,zmin,dz,nStep);
    cout << endl;
	
	
	// Write doppler parameter distribution to a file.
	// Calculates h(b) given in eqn 6 in Inoue & Iwata 2008
	cout <<"     Calculating absorber b-parameter distribution"<<endl;
    double bsig=23;
    DopplerParDistribution dopplerParDist(bsig);
    double bmin=2, bmax=200;
    double db=(bmax-bmin)/(nStep-1);
    outfile= outfileroot + "_hb.txt";
    dopplerParDist.writeToFile(outfile,bmin,db,nStep);
    cout << endl;
	
	
    // Generate Monte Carlo distributions
    cout <<"     Generate Monte Carlo distribution of absorbers "<<endl;
    RandomGenerator rg;
	ProbabilityDistAbsorbers probDistAbsorbers(rg, absorberZDist, 
	                                           hiColumnDensity, dopplerParDist);
	vector<double> log10NHIvals, log10gvals;
    probDistAbsorbers.returnColumnDensityGrid(log10NHIvals, log10gvals);
    outfile =  outfileroot + "_coldensgrid.txt";
	outp.open(outfile.c_str(), ofstream::out);
	for (int i=0; i<log10gvals.size(); i++)
	    outp <<"  "<< log10NHIvals[i] <<"  "<< log10gvals[i] << endl;
	outp.close();
    
    
    double zCurrent = 0;
    /*cout <<"     Write some files for checking ... "<<endl;
    long nTrial = 1000000;
    cout <<"     Number of absorbers in these trials = "<< nTrial <<endl;

    // Simulate NEXT absorber z distribution for CURRENT absorber being at z=0
	double zCurrent = 0;
	cout <<"     Simulate a distribution of the NEXT absorber z for CURRENT";
    cout <<" absorber being at "<< zCurrent <<endl;
	outfile="testfiles/prz.txt";
	probDistAbsorbers.writeZDistribution(outfile,zCurrent,nTrial);
	
	// Simulate doppler distribution
	cout <<"     Simulate a distribution of dopper parameters"<<endl;
	outfile="testfiles/prb.txt";
	probDistAbsorbers.writeDopplerDistribution(outfile,nTrial);
	
	// Simulate column density distribution
	cout <<"     Simulate a distribution of column densities"<<endl;
	outfile="testfiles/prc.txt";
	probDistAbsorbers.writeNHiDistribution(outfile,nTrial);
	cout <<"     End writing check files "<<endl;
	cout << endl;*/
	
	
	// Now simulate a line of sight
	cout <<"     Simulate one line of sight starting at z = 0"<<endl;
	double zStart = 0;
	zCurrent = zStart;
	double zMax =6;
	vector<double> redshifts, dopplers, columnDensities, redshifts2;
	outfile = outfileroot + "_lineofsight.txt";
	probDistAbsorbers.simulateLineOfSight(zStart,zMax,redshifts, dopplers, 
	                                                columnDensities, outfile);
	/*for (int i=0; i<redshifts.size(); i++) {
	    double z = probDistAbsorbers.simulateAbsorberRedshift(zMax);
	    redshifts2.push_back(z);
	    }
    cout << endl;
    outfile = outfileroot + "_lineofsight2.txt";
    probDistAbsorbers.writeToFile(outfile, redshifts2, dopplers, columnDensities);
    
    // One absorber, look at transmission as a function of observed z
    double bval = 30;
    double NHI = 1e14;
    double zA = 3.;
    
    // emission frame wavelengths
    double lmin = 9.e-8, lmax = 1.25e-7;
	int nl = 1000;
	double dl = (lmax - lmin)/(nl - 1);
    
    // Study Voigt profiles
    outfile =  outfileroot + "_voigtProfiles.txt";
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
	    inp.clear(ios::failbit);
	    cout << "     Writing to "<< outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (int i=0; i<nl; i++) {

	        double lam = lmin + i*dl;
	        outp << lam <<"  ";
	        for (int n=2; n<=31; n++) {
	            VoigtProfile voigtProf(bval, n);
	            double vprof = voigtProf(lam);
	            outp << vprof <<"  ";
	            }
            outp << endl;
            
	        }
        outp.close();
	    }
    else
        cout << "Error...file " << outfile.c_str() << " exists" << endl;
    
    // build up lyman-series cross-section calculation
    OpticalDepth optDepth;
    outfile =  outfileroot + "_lyScrossSec.txt";
	outp.open(outfile.c_str(), ofstream::out);
	
	double constants = (sqrt(PI)*ELECTRON_CHARGE_STATC*ELECTRON_CHARGE_STATC)/
                                        (ELECTRON_MASS_G*SPEED_OF_LIGHT_CMS);// in cm^2/s
    cout << "     Contant = "<< constants <<endl;

    for (int i=0; i<nl; i++) {

        double lam = lmin + i*dl;
        outp << lam <<"  "<< constants <<"  ";
        for (int n=2; n<=31; n++) {
        
            // doppler width
            double dW =  optDepth.returnDopplerWidthFreq(n,bval);
        
            // line strength
            double fi = optDepth.returnOscillatorStrength(n);
            
            // voigt profile
            VoigtProfile voigtProf(bval, n);
            double vprof = voigtProf(lam);
            
            double val = constants*(fi/dW)*vprof;
            
            outp << fi <<"  "<< dW << "  "<< vprof << "  "<< val <<"  ";
            }
        outp << endl;
        
        }
    outp.close();
    
    
    
	outfile =  outfileroot + "_crossSections.txt";
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
	    inp.clear(ios::failbit);
	    cout << "     Writing to "<< outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (int i=0; i<nl; i++) {
	        //cout <<"     On wavelength "<<i+1<<" of "<<nl<<endl;
	        double lam = lmin + i*dl;

	        double freq = SPEED_OF_LIGHT_MS/lam;
	        double transLC = optDepth.returnLymanContinuumOpticalDepth(freq,NHI);
	        double transLS = optDepth.returnLymanSeriesOpticalDepth(freq,NHI,bval);
	        
	        //double trans = optDepth(lam);
	        outp << lam <<"  "<< transLC <<"  "<< transLS <<endl;
	        }
        outp.close();
	    }
    else
        cout << "Error...file " << outfile.c_str() << " exists" << endl;*/
    
    
    /*double zAMin = 0.1, zAMax = 3;
    int nA = 100;
    double dzA = (zAMax - zAMin)/(nA - 1);
    
    for (int i=0; i<nA; i++)
        double zA = zAMin + i*dzA;
        OpticalDepth optDepth(zA,NHI,bval);*/
        
    /*// make some fake data, needed for the below
    vector<double> ztest, nlow, btest, nhigh;
    int nAbsorber = redshifts.size();
    double zamin = 1.5, zamax =2.99;
    double dza = (zamax - zamin)/(nAbsorber - 1);
    for (int i=0; i<nAbsorber; i++) {
        double zA = zamin + dza*i;
        //cout <<"     zAbs = "<<zA<<endl;
        ztest.push_back(0.006*i);
        nlow.push_back(1e14);
        btest.push_back(bval);
        nhigh.push_back(1e18);
        }
    cout <<"     N absorbers = "<< nAbsorber <<", max z = "<< ztest[nAbsorber-1]<<endl;
        
    int nLya = 2; // lyman-alpha only
   
    // NOT REAL DISTRIBUTION, LOW NHI
    // Lyman-alpha only
    LineOfSightTrans lightOfSightTrans1(ztest, btest, nlow);
    lightOfSightTrans1.setnLineMax(nLya);
    lightOfSightTrans1.setLymanSeriesOnly();
    // Lyman series only  
    LineOfSightTrans lightOfSightTrans2(ztest, btest, nlow); 
    lightOfSightTrans2.setLymanSeriesOnly();
    // Lyman continuum only
    LineOfSightTrans lightOfSightTrans3(ztest, btest, nlow);
    lightOfSightTrans3.setLymanContinuumOnly();
    // All
    LineOfSightTrans lightOfSightTrans4(ztest, btest, nlow); 
    
    // NOT REAL DISTRIBUTION, HIGH NHI
    // Lyman-alpha only
    LineOfSightTrans lightOfSightTrans5(ztest, btest, nhigh);
    lightOfSightTrans5.setnLineMax(nLya);
    lightOfSightTrans5.setLymanSeriesOnly();
    // Lyman series only  
    LineOfSightTrans lightOfSightTrans6(ztest, btest, nhigh); 
    lightOfSightTrans6.setLymanSeriesOnly();
    // Lyman continuum only
    LineOfSightTrans lightOfSightTrans7(ztest, btest, nhigh);
    lightOfSightTrans7.setLymanContinuumOnly();
    // All
    LineOfSightTrans lightOfSightTrans8(ztest, btest, nhigh); 

    // REAL-Z DISTRIBUTION, LOW NHI
    // Lyman-alpha only
    LineOfSightTrans lightOfSightTrans9(redshifts, btest, nlow);
    lightOfSightTrans9.setnLineMax(nLya);
    lightOfSightTrans9.setLymanSeriesOnly();
    // Lyman series only  
    LineOfSightTrans lightOfSightTrans10(redshifts, btest, nlow); 
    lightOfSightTrans10.setLymanSeriesOnly();
    // Lyman continuum only
    LineOfSightTrans lightOfSightTrans11(redshifts, btest, nlow);
    lightOfSightTrans11.setLymanContinuumOnly();
    // All
    LineOfSightTrans lightOfSightTrans12(redshifts, btest, nlow); 
    
    // REAL-Z DISTRIBUTION, HIGH NHI
    // Lyman-alpha only
    LineOfSightTrans lightOfSightTrans13(redshifts, btest, nhigh);
    lightOfSightTrans13.setnLineMax(nLya);
    lightOfSightTrans13.setLymanSeriesOnly();
    // Lyman series only  
    LineOfSightTrans lightOfSightTrans14(redshifts, btest, nhigh); 
    lightOfSightTrans14.setLymanSeriesOnly();
    // Lyman continuum only
    LineOfSightTrans lightOfSightTrans15(redshifts, btest, nhigh);
    lightOfSightTrans15.setLymanContinuumOnly();
    // All
    LineOfSightTrans lightOfSightTrans16(redshifts, btest, nhigh);
    
    // REAL-Z DISTRIBUTION, REAL NHI
    // Lyman-alpha only
    LineOfSightTrans lightOfSightTrans17(redshifts, btest, columnDensities);
    lightOfSightTrans17.setnLineMax(nLya);
    lightOfSightTrans17.setLymanSeriesOnly();
    // Lyman series only  
    LineOfSightTrans lightOfSightTrans18(redshifts, btest, columnDensities); 
    lightOfSightTrans18.setLymanSeriesOnly();
    // Lyman continuum only
    LineOfSightTrans lightOfSightTrans19(redshifts, btest, columnDensities);
    lightOfSightTrans19.setLymanContinuumOnly();
    // All
    LineOfSightTrans lightOfSightTrans20(redshifts, btest, columnDensities);
    
    // ALL REAL
    // Lyman-alpha only
    LineOfSightTrans lightOfSightTrans21(redshifts, dopplers, columnDensities);
    lightOfSightTrans21.setnLineMax(nLya);
    lightOfSightTrans21.setLymanSeriesOnly();
    // Lyman series only  
    LineOfSightTrans lightOfSightTrans22(redshifts, dopplers, columnDensities); 
    lightOfSightTrans22.setLymanSeriesOnly();
    // Lyman continuum only
    LineOfSightTrans lightOfSightTrans23(redshifts, dopplers, columnDensities);
    lightOfSightTrans23.setLymanContinuumOnly();
    // All
    LineOfSightTrans lightOfSightTrans24(redshifts, dopplers, columnDensities);*/
    
    /*double zwant = 2.78, nwant = 6.5e17;
    int cnt=0;
    vector<double> redshiftsFudge, dopplersFudge, columnDensitiesFudge;
    for (int i=0; i<redshifts.size(); i++) {
    
        //cout << i <<endl;
        double zv = redshifts[i];
        double bv = dopplers[i];
        double nv = columnDensities[i];
        
        redshiftsFudge.push_back(zv);
        dopplersFudge.push_back(bv);
        columnDensitiesFudge.push_back(nv);
        
        if ( (zv<zwant)&&(redshifts[i+1]>zwant) ){
            cout << "     Adding absorber between index "<<i<<" and "<<i+1<<endl;
            redshiftsFudge.push_back(zwant);
            dopplersFudge.push_back(bval);
            columnDensitiesFudge.push_back(nwant);
            cnt++;
            }
        }
    if (cnt>1)
        cout <<"EEEP!"<<endl;
    if (cnt<1)
        cout <<"     Did not add absorber!"<< endl;
    LineOfSightTrans lightOfSightTrans10(redshiftsFudge, dopplersFudge, columnDensitiesFudge);*/
    
    
    /*
    double zSource = 3.; //
    
    VoigtProfile voigtProfile(bval, 2);
    
    // observed frame wavelengths
    lmin = 3e-7, lmax = 5.5e-7;
	nl = 100000;
	dl = (lmax - lmin)/(nl - 1);
	
	outfile = outfileroot + "_transvary.txt";
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
	    inp.clear(ios::failbit);
	    cout << "     Writing to "<< outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (int i=0; i<nl; i++) {
	        cout <<"     On wavelength "<<i+1<<" of "<<nl<<endl;
	        double lam = lmin + i*dl;
	        
	        double transLS1 = lightOfSightTrans1(lam, zSource);
	        double transLS2 = lightOfSightTrans2(lam, zSource);
	        double transLS3 = lightOfSightTrans3(lam, zSource);
	        double transLS4 = lightOfSightTrans4(lam, zSource);
	        //double transLS5 = lightOfSightTrans5(lam, zSource);
            //double transLS6 = lightOfSightTrans6(lam, zSource);
            //double transLS7 = lightOfSightTrans7(lam, zSource);
            //double transLS8 = lightOfSightTrans8(lam, zSource);
            //double transLS9 = lightOfSightTrans9(lam, zSource);
            //double transLS10 = lightOfSightTrans10(lam, zSource);
            //double transLS11 = lightOfSightTrans11(lam, zSource);
            //double transLS12 = lightOfSightTrans12(lam, zSource);
            //double transLS13 = lightOfSightTrans13(lam, zSource);
            //double transLS14 = lightOfSightTrans14(lam, zSource);
            //double transLS15 = lightOfSightTrans15(lam, zSource);
            //double transLS16 = lightOfSightTrans16(lam, zSource);
            //double transLS17 = lightOfSightTrans17(lam, zSource);
            //double transLS18 = lightOfSightTrans18(lam, zSource);
            //double transLS19 = lightOfSightTrans19(lam, zSource);
            //double transLS20 = lightOfSightTrans20(lam, zSource);
            double transLS21 = lightOfSightTrans21(lam, zSource);
            double transLS22 = lightOfSightTrans22(lam, zSource);
            double transLS23 = lightOfSightTrans23(lam, zSource);
            double transLS24 = lightOfSightTrans24(lam, zSource);
            
            //int iDistant = findClosestElement(ztest,zSource);
            //cout << iDistant <<endl;*/
            
            /*double crossSec = 0., crossSec2 = 0., profSum=0.;
            vector<double> profs, lambdas;
            for (int j=0; j<=iDistant; j++) {
                double lambdaE = lam/(1. + ztest[j]);
                lambdas.push_back(lambdaE);
                //double freq = SPEED_OF_LIGHT_MS/lambdaE;
                double vf = voigtProfile(lambdaE);
                profs.push_back(vf);
                profSum +=vf;
                //crossSec += lightOfSightTrans1.returnLymanLineCrossSection(2, freq, bval);
                crossSec2 += 2.5255e-14*vf;
                }
            if (profs.size() >2)
                cout << profs.size() <<"?"<<endl;
            //double transTest1 = exp(-1e14*crossSec);
            //double transTest1a = exp(-1e15*crossSec);
            //double transTest1b = exp(-1e16*crossSec);
            //double transTest1c = exp(-1e17*crossSec);
            //double transTest5 = exp(-1e18*crossSec);
            double transTest = exp(-1e18*crossSec2);*/
            /*
	        outp << lam <<"  "<< lam/(1+zSource) <<"  "<< transLS1  <<"  ";
	        outp << transLS2  <<"  "<< transLS3  <<"  "<< transLS4  <<"  ";
	        //outp << transLS5 ;
	        //outp <<"  "<< transLS6  <<"  "<< transLS7  <<"  ";
	        //outp << transLS8  <<"  "<< transLS9  <<"  "<< transLS10 <<"  ";
	        //outp << transLS11 <<"  "<< transLS12 <<"  "<< transLS13 <<"  ";
	        //outp << transLS14 <<"  "<< transLS15 <<"  "<< transLS16 <<"  "; 
	        //outp << transLS17 <<"  "<< transLS18 <<"  "<< transLS19 <<"  ";
	        //outp << transLS20 <<"  ";
	        outp << transLS21 <<"  "<< transLS22 <<"  ";
	        outp << transLS23 <<"  "<< transLS24 <<"  ";
	        outp << endl;
	        }
        outp.close();
	    }
    else
        cout << "Error...file " << outfile.c_str() << " exists" << endl;*/
    
	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testsimulateIGM.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testsimulateIGM.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " testsimulateIGM.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of testsimulateIGM.cc program  Rc= " << rc << endl;
  return rc;	
}

