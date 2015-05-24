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
#include "poly.h"
#include "geneutils.h"
#include "time.h"
#include "igm.h"
#include <vector>


//class ACalcs {
//    static double sigmaLymanLimitCM2;         /**< cross-section at the lyman limit in cm^2     */
//    static double freqLymanLimitInvSec;       /**< frequency of the lyman limit in s^-1         */
//    static int nLymanAlpha;                   /**< starting level of a lyman alpha transition   */
//    static int nLineMaxMax;                   /**< maximum Lyman series line possible to use    */
//};




using namespace std;

int main(int narg, char* arg[]) {
    double versionNum = 2.01;
    cout << endl << endl << "****           Version " << versionNum << "            ****" << endl << endl;

    // ********************************
    // Initialize RNG
    // ********************************

    cout << "       Initializing the RNG" << endl;
    cout << "       Time: " << time(NULL) << endl;
    RandomGenerator rg;
    ///rg.SetSeed(time(NULL));
    rg.AutoInit();
    uint_2 seeds[3];
    rg.GetSeed(seeds);
    cout << "Random Seed: " << seeds[0] << " " << seeds[1] << " " << seeds[2]<< endl;


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
    // Define the problem
    // ********************************

    int nmax = 10;                                  //Highest Lyman Line to go to
    int nLoS = 400;                                 //The number of lines of sight to examine for each source redshift
    double zSource = 1.8;                           //Starting redshift of the background source
    int nSequence = -1;                             //Sequence number for determining multiple runs

    string transFileName = "";                      //File name for the Transmission file
    string LoSFileName = "";                        //File name for file containing LoS info
    string outputPath = "./testfiles/";             //Path from pwd to output files

    //--- decoding command line arguments 
    cout << " ==== decoding command line arguments ===="<<endl;
    char c;
    while((c = getopt(narg,arg,"n:z:s:h")) != -1) 
      {
	switch (c) 
	  {
	  case 'n' :
	    nLoS = atoi(optarg);
	    break;
	  case 'z' :
	    zSource = atof(optarg);
	    break;
	  case 's' :
	    nSequence = atoi(optarg);
	    break;
	  case 'h' :
	  default :
	    cout << "Usage: testsimulateIGM [-n N lines of sight] [-z redshift] [-h] [-s sequence]";
	    return -1;
	  }
      }

  //-- end command line arguments

    cout << "nSequence: " << nSequence << " nLoS: " << nLoS << " zSource: " << zSource << endl;

    // Modify the wavelength range based upon the galaxy redshift
    double lambdaMinA = 3000.;                    //Minimum WL in Angstroms (LSST minimum wavelength)
    double lambdaMaxA = 100. + 1250.*(1+zSource); //Maximum WL in Angstroms (About 100A past Lya from an absorber at zSource) (<11,000A)

// Incorrect way for min and max wavelength 
//    double lambdaMinA = 3000.*(1.+zSource-0.1)/(zSource+1.);         //Minimum WL in Angstroms (LSST minimum wavelength)
//    double lambdaMaxA = 3000.;                      //Maximum WL in Angstroms (About 100A past Lya from an absorber at zSource) (<11,000A)

    double lambdaResA = 0.1;                        //WL Resolution in Angstroms
//    double lambdaResA = 0.125;                        //WL Resolution in Angstroms

    double zAbsorberMax = zSource;// + zSource*0.1;    //Maximum absorber redshift
    double zAbsorberStart = 0.0;                    //Starting absorber redshift   

    stringstream ss;
    ss << zSource << "zSource_" << nLoS << "nLos_" << nmax << "LynMax";
    if (nSequence>0) ss << "_" << nSequence;
    ss << ".dat";
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
    LoSData.open(LoSFileName.c_str(),ios::out | ios::binary);

    LoSData.write((char*)&zSource, sizeof(double));
    LoSData.write((char*)&nLoS, sizeof(int));

    for(int i=0; i<nLoS; i++) {

      if (i%100 == 0 ) 	cout << "nSequence: " << nSequence << ". Line of Sight: " << i << endl;

        pdist.simulateLineOfSight(zAbsorberStart, zAbsorberMax, redshifts, dopplerPars,
                                  columnDensities, outfile);

        LineOfSightTrans lineOfSight(redshifts, dopplerPars, columnDensities);
        lineOfSight.setMaxLine(nmax);
        lineOfSight.setLymanAll();

	int nredshifts = redshifts.size();
	LoSData.write((char*)&nredshifts, sizeof(int));
        for(int j=0; j<redshifts.size(); j++)
	  LoSData.write((char*)&redshifts[j],sizeof(double));

	int ndopplers = dopplerPars.size();
	LoSData.write((char*)&ndopplers, sizeof(int));
        for(int j=0; j<dopplerPars.size(); j++)
	  LoSData.write((char*)&dopplerPars[j],sizeof(double));

	int ncoldens = columnDensities.size();
	LoSData.write((char*)&ncoldens, sizeof(int));
        for(int j=0; j<columnDensities.size(); j++)
	  LoSData.write((char*)&columnDensities[j], sizeof(double));


        LoSvector.push_back(lineOfSight);
        redshifts.clear();
        dopplerPars.clear();
        columnDensities.clear();
    }
    cout << "nSequence: " << nSequence << ". End Line of Sight simulations." << endl;

    LoSData.close();

// */

    // ********************************
    // Make and write the wavelength vector
    // ********************************

// /*

    cout << "nSequence: " << nSequence <<  ". Preparing transmission data." << endl;
    ofstream transmissionData;
    transmissionData.open(transFileName.c_str(), ios::out | ios::binary);

    //    transmissionData << "First line is wavelength in Angstroms each successive line is transmission along a line of sight" << endl
    //            << "zSource: " << zSource << endl
    //            << "nLoS: "    << nLoS << endl;
    transmissionData.write((char*)&zSource,sizeof(double));
    transmissionData.write((char*)&nLoS,sizeof(int));

    vector<double> wavelengths;
    double lambdaCA = lambdaMinA;
    while(lambdaCA < lambdaMaxA) {
        double wavelengthTemp = lambdaCA*pow(10.,-10);
        wavelengths.push_back(wavelengthTemp);

	//        transmissionData << wavelengthTemp << " ";
	transmissionData.write((char*)&wavelengthTemp, sizeof(double));
//        transmissionData << wavelengthTemp << endl;

        lambdaCA += lambdaResA;
    }
    double endmark = -100.0;
    transmissionData.write((char*)&endmark, sizeof(double));

    //    transmissionData << endl;
//    transmissionData << " ";

// */

    // ********************************
    // Calculate the Transmission
    // ********************************

// /*

    cout << "nSequence: " << nSequence << ". Writing transmission data." << endl;
    int lineLength = 0;

    for(int k=0; k<nLoS; k++) {
      //      cout << "Line of Sight: " << k << endl;
      if (k%100 == 0 ) 	cout << "nSequence: " << nSequence << ". Line of Sight: " << k << endl;

      vector<double> transmission;
      for(int i=0; i<wavelengths.size(); i++) {
	double transWL = LoSvector[k].returnTransmission(wavelengths[i], zSource);
	transmission.push_back(transWL);
      }
      //      cout << "  Obtained transmission." << endl;

      for(int i=0; i<transmission.size(); i++) {
	//	transmissionData << transmission[i] << " ";
	transmissionData.write((char*)&transmission[i], sizeof(double));
//            transmissionData << transmission[i] << endl;
      }

      // write out end mark
      double endmark = -100;
      transmissionData.write((char*)&endmark, sizeof(double));

      //      transmissionData << endl;
      //      cout << "  Wrote transmission." << endl;
//	transmissionData << " ";
      
//      cout << "Transmission size: " << transmission.size() << endl;
    }

    transmissionData.close();

    cout << "nSequence: " << nSequence << ". Finished IGM transmission simulations." << endl;
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

