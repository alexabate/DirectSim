#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>
#include <sstream>

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

int main()
{
    double versionNum = 2.01;
    cout << endl << endl << "****           Version " << versionNum << "            ****" << endl << endl;

    // ********************************
    // Initialize RNG
    // ********************************

    cout << "       Initializing the RNG" << endl;
    RandomGenerator rg;
    rg.SetSeed(time(NULL));


    // ********************************
    // Initialize the distributions
    // ********************************

    cout << "       Initializing the Distributions" << endl;
    HIColumnDensity colDensityDist;
    AbsorberRedshiftDistribution zDist;
    DopplerParDistribution bDist;

    ProbabilityDistAbsorbers pdist(rg, zDist, colDensityDist, bDist);


    // ********************************
    // Create some vectors and other Variables
    // ********************************

    vector<double> redshifts;
    vector<double> columnDensities;
    vector<double> dopplerPars;
    vector<LineOfSightTrans> LoSvector;
    string outfile = "";


    // ********************************
    // Define the the problem
    // ********************************

    int nmax = 10;                                  //Highest Lyman Line to go to
    int nLoS = 400;                                 //The number of lines of sight to examine for each source redshift
    double zSource = 1.8;                           //Starting redshift of the background source


    double lambdaMinA = 3000.;                    //Minimum WL in Angstroms (LSST minimum wavelength)
    double lambdaMaxA = 100. + 1250.*(1+zSource); //Maximum WL in Angstroms (About 100A past Lya from an absorber at zSource) (<11,000A)

// Incorrect way for min and max wavelength 
//    double lambdaMinA = 3000.*(1.+zSource-0.1)/(zSource+1.);         //Minimum WL in Angstroms (LSST minimum wavelength)
//    double lambdaMaxA = 3000.;                      //Maximum WL in Angstroms (About 100A past Lya from an absorber at zSource) (<11,000A)


    double lambdaResA = 0.1;                        //WL Resolution in Angstroms
//    double lambdaResA = 0.125;                        //WL Resolution in Angstroms

    double zAbsorberMax = zSource;// + zSource*0.1;    //Maximum absorber redshift
    double zAbsorberStart = 0.0;                    //Starting absorber redshift   

    string transFileName = "";                      //File name for the Transmission file
    string LoSFileName = "";                        //File name for file containing LoS info
    string outputPath = "./testfiles/";             //Path from pwd to output files

    stringstream ss;
    ss << zSource << "zSource_" << nLoS << "nLos_" << nmax << "LynMax.dat";
    transFileName = outputPath + "Transmission_" + ss.str();
    LoSFileName = outputPath + "LineOfSightData_" + ss.str();

    string avgNAbsorbsOutFile = outputPath + "avgNAbsorbers.out";

    // ********************************
    // Initialize Clock
    // ********************************

/*
    clock_t starttime, endtime;
    starttime = clock();
*/

    // ********************************
    // Count number of Absorbers per redshift
    // ********************************
/*
    ofstream countAbsorbs;
    countAbsorbs.open(avgNAbsorbsOutFile.c_str());

    int totalAbs = 0;

    for(int j=0; j<16; j++) {
        for(int i=0; i<nLoS; i++) {
            pdist.simulateLineOfSight(zAbsorberStart, zAbsorberMax+j*0.1, redshifts, dopplerPars,
                                      columnDensities, outfile);

            totalAbs += columnDensities.size();

            redshifts.clear();
            dopplerPars.clear();
            columnDensities.clear();
        }

        countAbsorbs << zAbsorberMax+j*0.1 << "     " << totalAbs/400. << endl;
        totalAbs = 0;
    }

    countAbsorbs.close();
// */

    // ********************************
    // Simulate nLoS lines of sight and record them to a file
    // ********************************

// /*
    ofstream LoSData;
    LoSData.open(LoSFileName.c_str());

    LoSData << "Each line of sight goes redshift, doppler parameter, column density" << endl
            << "zSource: " << zSource << endl
            << "nLoS: "    << nLoS << endl;

    for(int i=0; i<nLoS; i++) {
        pdist.simulateLineOfSight(zAbsorberStart, zAbsorberMax, redshifts, dopplerPars,
                                  columnDensities, outfile);

        LineOfSightTrans lineOfSight(redshifts, dopplerPars, columnDensities);
        lineOfSight.setMaxLine(nmax);
        lineOfSight.setLymanAll();

        for(int j=0; j<redshifts.size(); j++)
            LoSData << redshifts[j] << " ";
        LoSData << endl;
        for(int j=0; j<dopplerPars.size(); j++)
            LoSData << dopplerPars[j] << " ";
        LoSData << endl;
        for(int j=0; j<columnDensities.size(); j++)
            LoSData << columnDensities[j] << " ";
        LoSData << endl;

        LoSvector.push_back(lineOfSight);
        redshifts.clear();
        dopplerPars.clear();
        columnDensities.clear();
    }

    LoSData.close();

// */

    // ********************************
    // Make and write the wavelength vector
    // ********************************

// /*
    ofstream transmissionData;
    transmissionData.open(transFileName.c_str());

    transmissionData << "First line is wavelength in Angstroms each successive line is transmission along a line of sight" << endl
            << "zSource: " << zSource << endl
            << "nLoS: "    << nLoS << endl;

    vector<double> wavelengths;
    double lambdaCA = lambdaMinA;
    while(lambdaCA < lambdaMaxA) {
        double wavelengthTemp = lambdaCA*pow(10.,-10);
        wavelengths.push_back(wavelengthTemp);

        transmissionData << wavelengthTemp << " ";
//        transmissionData << wavelengthTemp << endl;

        lambdaCA += lambdaResA;
    }

    transmissionData << endl;
//    transmissionData << " ";

// */

    // ********************************
    // Calculate the Transmission
    // ********************************

// /*

    int lineLength = 0;

    for(int k=0; k<nLoS; k++) {
        vector<double> transmission;
        for(int i=0; i<wavelengths.size(); i++) {
            double transWL = LoSvector[k].returnTransmission(wavelengths[i], zSource);
            transmission.push_back(transWL);
        }

        for(int i=0; i<transmission.size(); i++) {
            transmissionData << transmission[i] << " ";
//            transmissionData << transmission[i] << endl;


        }
        transmissionData << endl;
//	transmissionData << " ";

        cout << transmission.size() << endl;
    }

    transmissionData.close();
// */


    // ********************************
    // End clock and output time
    // ********************************
/*
    endtime = clock();
    float runTimeTicks = ((float)endtime - (float)starttime);
    float runTimeSeconds = runTimeTicks / CLOCKS_PER_SEC;

    string clockFileName = outputPath + "clock_" + ss.str();
    ofstream clockOutput;
    clockOutput.open(clockFileName.c_str());
    clockOutput << "Run time was " << runTimeSeconds << " for " << nLoS
                << " lines of sight." << endl;
    clockOutput.close();
    

*/
    return 0;
};

