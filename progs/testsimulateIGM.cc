#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

#include "fiosinit.h"
#include "mydefrg.h"

#include "constcosmo.h"
#include "cosmocalcs.h"
#include "igm.h"
#include "poly.h"
#include "geneutils.h"
#include "time.h"
#include <vector>

using namespace std;

void findOpticalDepthAlongLineOfSight(double lambda, vector<double>& redshifts, vector<double>& dopplerPars,
                                           vector<double>& colDens);
void printDopplerParameters(ProbabilityDistAbsorbers &pdist, int n);
void printColDens(ProbabilityDistAbsorbers &pdist, int n);
void printRedshift(ProbabilityDistAbsorbers &pdist, int numTests, double z);

int main()
{
//    cout << endl << endl << endl << endl;
    double versionNum = 1.21;
    cout << endl << endl << "****           Version " << versionNum << "            ****" << endl << endl;



    // Initialize RNG
    cout << "       Initializing the RNG" << endl;
    RandomGenerator rg;
    rg.SetSeed(time(NULL));

    // Initialize the distributions
    cout << "       Initializing the Distributions" << endl;
    HIColumnDensity colDensityDist;
    AbsorberRedshiftDistribution zDist;
    DopplerParDistribution bDist;

    ProbabilityDistAbsorbers pdist(rg, zDist, colDensityDist, bDist);

    // Simulate the line of sight
    cout << "       Simulating the line of sight" << endl;
    vector<double> redshifts;
    vector<double> columnDensities;
    vector<double> dopplerPars;
    double zMax = 6.0;
    double zStart = 0.0;
    string outfile = "";

    pdist.simulateLineOfSight(zStart, zMax, redshifts, dopplerPars, 
                                columnDensities, outfile);

    // Define the the problem
    int nmax = 10;                  //Highest Lyman Line to go to
    double zSource = 5.0;           //Redshift of the background source
    double lambdaMinA = 4200.;      //Minimum WL in Angstroms
    double lambdaMaxA = 7800.;      //Maximum WL in Angstroms
    double lambdaResA = 0.1;        //WL Resolution in Angstroms
    string transFileName = "./testfiles/testsimulateIGM_IandI2008_zSource5";



/*
            // TEST CDFS

    int numTests = 1000000;

    // Print n doppler parameters
    printDopplerParameters(pdist, numTests);

    // Print n column densities
    printColDens(pdist, numTests);

    // Print n redshifts at each functional range
    double z1 = 0.5;
    double z2 = 2.0;
    double z3 = 4.5;
    printRedshift(pdist, numTests, z1);
    printRedshift(pdist, numTests, z2);
    printRedshift(pdist, numTests, z3);
// */


/* 
            // FIND OPTICAL DEPTH THROUGH EACH ABSORBER AND PRINT

    // Find the optical Depth of each absorber along a line of sight
    double lambda = 4000*pow(10.,-10);
    findOpticalDepthAlongLineOfSight(lambda, redshifts, dopplerPars, columnDensities);
// */


// /*

            // FIND THE TRANSMISSION ON A LINE OF SIGHT FOR A VECTOR OF WAVELENGTHS

    // Set up the Line of Sight Transmission Object
    cout << "       Creating LineOfSightTransmission Object" << endl;
    LineOfSightTrans lineOfSight(redshifts, dopplerPars, columnDensities);
    lineOfSight.setMaxLine(nmax);
    lineOfSight.setLymanAll();


    // Make the vector for the wavelength range
    cout << "       Defining the Wavelength Range" << endl;
    vector<double> WL;
    double lambdaCA = lambdaMinA;
    while(lambdaCA < lambdaMaxA) {
        WL.push_back(lambdaCA*pow(10.,-10));
        lambdaCA += lambdaResA;
    }


    // Calculate the transmission at each wavelength
    cout << "       Beginning transmission calculation" << endl;
    vector<double> transmission;
    for(int i=0; i<WL.size(); i++) {
        double transWL = lineOfSight.returnTransmission(WL[i], zSource);
        transmission.push_back(transWL);
    }


    // Output the Transmission to a file
    cout << "       Outputting transmission" << endl;
    ofstream ofile;
    ofile.open(transFileName.c_str());
    for(int i=0; i<transmission.size(); i++) {
        ofile << WL[i] << " " << transmission[i] << endl;
    }
    ofile.close();
// */



    return 0;
};

// Given a line of sight and wavelength, print the optical depth through every absorber
void findOpticalDepthAlongLineOfSight(double lambda, vector<double>& redshifts, vector<double>& dopplerPars,
                                           vector<double>& colDens)
{
    OpticalDepth opti;
    opti.setMaxLine(10);
    opti.setLymanSeriesOnly();

    ofstream ofileOptiDepth;
    ofileOptiDepth.open("./testfiles/testsimulateIGM_absorberOpticalDepth");

    int numAbsorbers = redshifts.size();
    for(int i=0; i<numAbsorbers; i++) {
        double tau = opti.returnObserverFrameOpticalDepth(lambda, redshifts[i],
                                    colDens[i], dopplerPars[i]);
//        if(tau > 5)
//            tau = 0;

        ofileOptiDepth << tau << endl;
    }

    return;
};


void printDopplerParameters(ProbabilityDistAbsorbers &pdist, int n)
{
    string fname = "./testfiles/testsimulateIGM_prb.txt";
    ofstream ofile;
    ofile.open(fname.c_str());

    cout << "Drawing doppler parameters." << endl;

    for(int i=0; i<n; i++) {
        ofile << pdist.drawDoppler() << endl;
    }

    ofile.close();

    return;
};

void printColDens(ProbabilityDistAbsorbers &pdist, int n)
{
    string fname = "./testfiles/testsimulateIGM_prn.txt";
    ofstream ofile;
    ofile.open(fname.c_str());

    cout << "Drawing column densities." << endl;

    for(int i=0; i<n; i++) {
        ofile << pdist.drawHIColumnDensity() << endl;
    }

    return;
};


void printRedshift(ProbabilityDistAbsorbers &pdist, int n, double z)
{
    string fname = "";

    if(z <= 1.2) 
        fname = "./testfiles/testsimulateIGM_prz1";
    else if(z <= 4.0)
        fname = "./testfiles/testsimulateIGM_prz2";
    else
        fname = "./testfiles/testsimulateIGM_prz3";

    ofstream ofile;
    ofile.open(fname.c_str());

    cout << "Drawing redshifts." << endl;

    for(int i=0; i<n; i++) {
        ofile << pdist.drawDeltaZ(z) << endl;
    }

    return;
};


















