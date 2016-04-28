#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>
#include <sstream>
#include <algorithm>

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

    HIColumnDensityLAF colDensityDistLAF;
    HIColumnDensityDLA colDensityDistDLA;

    AbsorberRedshiftDistributionLAF zDistLAF;
    AbsorberRedshiftDistributionDLA zDistDLA;

    DopplerParDistribution bDist;

    ProbabilityDistAbsorbers pdist(rg, zDistLAF, zDistDLA, colDensityDistLAF, colDensityDistDLA, bDist);


    // ********************************
    // Create some vectors and other Variables
    // ********************************

    vector<double> redshiftsLAF, redshiftsDLA;
    vector<double> columnDensitiesLAF, columnDensitiesDLA;
    vector<double> dopplerParsLAF, dopplerParsDLA;
    vector<LineOfSightTrans> LoSvectorLAF, LoSvectorDLA;
    string outfile = "";


    // ********************************
    // Define the problem
    // ********************************

    int nmax = 10;                                  //Highest Lyman Line to go to
    int nLoS = 200;                                 //The number of lines of sight to examine for each source redshift
    double zSource;                           	    //Source galaxy redshift 
    int nSequence = -1;                             //Sequence number for determining multiple runs

    string transFileName = "";                      //File name for the Transmission file
    string LoSFileNameLAF = "";                     //File name for file containing LAF LoS (z, NHI, b) 
    string LoSFileNameDLA = "";                     //File name for file containing DLA LoS (z, NHI, b)
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

    // Modify the wavelength range based upon the source galaxy redshift
    double lambdaMinA = 3000.;                    //Minimum WL in Angstroms (LSST minimum wavelength)
    double lambdaMaxA = 100. + 1250.*(1+zSource); //Maximum WL in Angstroms (About 100A past Lya from an absorber at zSource) (<11,000A)

    double lambdaResA = 0.1;                      //WL Resolution in Angstroms

    double zAbsorberMax = zSource;                //Maximum absorber redshift
    double zAbsorberStart = 0.0;                  //Starting absorber redshift   

    stringstream ss;
    ss << zSource << "zSource_" << nLoS << "nLos_" << nmax << "LynMax";
    if (nSequence>0) ss << "_" << nSequence;
    ss << ".dat";
    transFileName = outputPath + "Transmission_" + ss.str();
    LoSFileNameLAF = outputPath + "LineOfSightData_LAF_" + ss.str();
    LoSFileNameDLA = outputPath + "LineOfSightData_DLA_" + ss.str();

    string avgNAbsorbsOutFile = outputPath + "avgNAbsorbers.out";


    // ********************************
    // Simulate nLoS for LAF absorbers AND nLoS for DLA absorbers
    // simulate as 2 sets of lines of sight
    // record absorbers info to 2 files  
    // ********************************

    ofstream LoSDataLAF, LoSDataDLA;

//    LoSDataLAF.open(LoSFileNameLAF.c_str(),ios::out | ios::binary);
//    LoSDataDLA.open(LoSFileNameDLA.c_str(),ios::out | ios::binary);

    LoSDataLAF.open(LoSFileNameLAF.c_str(), ofstream::out);
    LoSDataDLA.open(LoSFileNameDLA.c_str(), ofstream::out);

    LoSDataLAF.write((char*)&zSource, sizeof(double));
    LoSDataDLA.write((char*)&zSource, sizeof(double));

    LoSDataLAF.write((char*)&nLoS, sizeof(int));
    LoSDataDLA.write((char*)&nLoS, sizeof(int));

  

    // ********************************
    //loop to simulate LAF absorbers for nLoS
    // ********************************

    for(int i=0; i<nLoS; i++) {

//      if (i%100 == 0 ) 	cout << "nSequence: " << nSequence << ". Line of Sight: " << i << endl;

        pdist.simulateLineOfSightLAF(zAbsorberStart, zAbsorberMax, redshiftsLAF, dopplerParsLAF,
                                  columnDensitiesLAF, outfile);

        LineOfSightTrans lineOfSight(redshiftsLAF, dopplerParsLAF, columnDensitiesLAF);
        lineOfSight.setMaxLine(nmax);
        lineOfSight.setLymanAll();

	int nredshifts = redshiftsLAF.size();
	LoSDataLAF.write((char*)&nredshifts, sizeof(int));
        for(int j=0; j<redshiftsLAF.size(); j++)
	  LoSDataLAF.write((char*)&redshiftsLAF[j],sizeof(double));

	int ndopplers = dopplerParsLAF.size();
	LoSDataLAF.write((char*)&ndopplers, sizeof(int));
        for(int j=0; j<dopplerParsLAF.size(); j++)
	  LoSDataLAF.write((char*)&dopplerParsLAF[j],sizeof(double));

	int ncoldens = columnDensitiesLAF.size();
	LoSDataLAF.write((char*)&ncoldens, sizeof(int));
        for(int j=0; j<columnDensitiesLAF.size(); j++)
	  LoSDataLAF.write((char*)&columnDensitiesLAF[j], sizeof(double));


        LoSvectorLAF.push_back(lineOfSight);
        redshiftsLAF.clear();
        dopplerParsLAF.clear();
        columnDensitiesLAF.clear();
    }
    cout << "nSequence: " << nSequence << ". End Line of Sight simulations." << endl;

    LoSDataLAF.close();



    // ********************************
    //loop to simulate DLA absorbers for nLoS
    // ********************************


    for(int i=0; i<nLoS; i++) {

//      if (i%100 == 0 ) 	cout << "nSequence: " << nSequence << ". Line of Sight: " << i << endl;

        pdist.simulateLineOfSightDLA(zAbsorberStart, zAbsorberMax, redshiftsDLA, dopplerParsDLA,
                                  columnDensitiesDLA, outfile);

        LineOfSightTrans lineOfSight(redshiftsDLA, dopplerParsDLA, columnDensitiesDLA);
        lineOfSight.setMaxLine(nmax);
        lineOfSight.setLymanAll();

	int nredshifts = redshiftsDLA.size();
	LoSDataDLA.write((char*)&nredshifts, sizeof(int));
        for(int j=0; j<redshiftsDLA.size(); j++)
	  LoSDataDLA.write((char*)&redshiftsDLA[j],sizeof(double));

	int ndopplers = dopplerParsDLA.size();
	LoSDataDLA.write((char*)&ndopplers, sizeof(int));
        for(int j=0; j<dopplerParsDLA.size(); j++)
	  LoSDataDLA.write((char*)&dopplerParsDLA[j],sizeof(double));

	int ncoldens = columnDensitiesDLA.size();
	LoSDataDLA.write((char*)&ncoldens, sizeof(int));
        for(int j=0; j<columnDensitiesDLA.size(); j++)
	  LoSDataDLA.write((char*)&columnDensitiesDLA[j], sizeof(double));


        LoSvectorDLA.push_back(lineOfSight);
        redshiftsDLA.clear();
        dopplerParsDLA.clear();
        columnDensitiesDLA.clear();
    }
    cout << "nSequence: " << nSequence << ". End Line of Sight simulations." << endl;

    LoSDataDLA.close();



    // ********************************
    // Make and write the wavelength vector
    // ********************************

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


    // ********************************
    // Calculate the Transmission
    // Add each LAF line of sight to each DLA line of sight
    // transmission should then be exp(-tau_LAF)*exp(-tau_DLA)
    // ********************************


    cout << "nSequence: " << nSequence << ". Writing transmission data." << endl;
    int lineLength = 0;

    for(int k=0; k<nLoS; k++) {
//      if (k%100 == 0 ) 	cout << "nSequence: " << nSequence << ". Line of Sight: " << k << endl;

      vector<double> transmission;
      for(int i=0; i<wavelengths.size(); i++) {

	double transWL = LoSvectorLAF[k].returnTransmission(wavelengths[i], zSource) * LoSvectorDLA[k].returnTransmission(wavelengths[i], zSource);
	transmission.push_back(transWL);
      }

      for(int i=0; i<transmission.size(); i++) {
	//	transmissionData << transmission[i] << " ";
	transmissionData.write((char*)&transmission[i], sizeof(double));
//            transmissionData << transmission[i] << endl;
      }

      // write out end mark
      double endmark = -100;
      transmissionData.write((char*)&endmark, sizeof(double));

    }

    transmissionData.close();

    cout << "nSequence: " << nSequence << ". Finished IGM transmission simulations." << endl;

    return 0;
};


