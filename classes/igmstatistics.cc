
#include "igmstatistics.h"


/************************************************************
*                                                           *
*                                                           *
*                       IGMVariance                         *
*                                                           *
*                                                           *
************************************************************/

IGMVariance::IGMVariance(double zshift)
{
    if(zshift > 2.9)
        throw ParmError("Redshift must be less than or equal to 2.9");

    char * pt = getenv("IGMLUTLOC");
    if(pt == NULL)
        throw ParmError("LOCAITON OF IGM LOOKUP TABLES ENVIRONMENT VARIABLE -IGMLUTLOC- NOT DEFINED ");

    redshift_ = zshift;
    numAbsorbers_ = calcNumberAbsorbers(redshift_);
    findLookupTables();
    
    // Open file
    ifstream ifile;
    string filename = pt+lookupTable_;
    ifile.open(filename.c_str());
    if(!ifile.is_open())
        throw ParmError("Lookup Table file was not opened");

    // Read in the file
    double tlam = 0.0, tmu = 0.0, tsig = 0.0;
    vector<double> readlam, readmu, readsig;
    while(!ifile.eof()) {
        ifile >> tlam >> tmu >> tsig;
        readlam.push_back(tlam);
        readsig.push_back(tsig);
    }
    ifile.close();

    // Shift the wavelength to the appropriate redshift
    for(int i=0; i<readlam.size(); i++)
        wavelengthVar_.push_back(readlam[i]*(1.0+redshift_)/(1.0+lookupZ_));

    // Shift the variance by the number of absorbers
    for(int i=0; i<readsig.size(); i++)
        variance_.push_back(readsig[i]*sqrt(numAbsorbers_)/sqrt(calcNumberAbsorbers(lookupZ_)));

    // Calculate the mean transmission
    calcMean();

};

double IGMVariance::calcNumberAbsorbers(double z)
{
    double temp = 390.7*z*z - 677.*z + 729.1;
    return temp;
};

void IGMVariance::findLookupTables()
{
    double lookupZ [15] = {1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9};
    stringstream ss;

    for(int i=0; i<14; i++) {
        if(lookupZ[i+1] >= redshift_) {
            lookupZ_ = lookupZ[i+1];
            lowZ_ = lookupZ[i];
            break;
        }
    }

    ss << "LookupTable_zSource" << lookupZ_ << ".tbl";
    lookupTable_ = ss.str();

    ss.str("");
    ss << "LookupTable_zSource" << lowZ_ << ".tbl";
    lookupTableLow_ = ss.str();

    highZ_ = lookupZ_;

    return;
};

void IGMVariance::calcMean()
{
    vector<double> wlL, wlU;

    // Load the lower file
    char * pt = getenv("IGMLUTLOC");
    ifstream ifileL;
    string filename1 = pt+lookupTableLow_;
    ifileL.open(filename1.c_str());
    if(!ifileL.is_open())
        throw ParmError("Lookup Table file was not opened");

    // Read in the lower file
    double tlam = 0.0, tmu = 0.0, tsig = 0.0;
    vector<double> readlamL, readmuL;
    while(!ifileL.eof()) {
        ifileL >> tlam >> tmu >> tsig;
        readlamL.push_back(tlam);
        readmuL.push_back(tmu);
    }
    ifileL.close();

    // Load the upper file
    ifstream ifileU;
    string filename2 = pt+lookupTable_;
    ifileU.open(filename2.c_str());
    if(!ifileU.is_open())
        throw ParmError("Lookup Table file was not opened");

    // Read in the upper file
    vector<double> readlamU, readmuU;
    while(!ifileU.eof()) {
        ifileU >> tlam >> tmu >> tsig;
        readlamU.push_back(tlam);
        readmuU.push_back(tmu);
    }
    ifileU.close();

    // Shift the lower wavelengths to the appropriate redshift
    for(int i=0; i<readlamL.size(); i++)
        wlL.push_back(readlamL[i]*(1.0+redshift_)/(1.0+lowZ_));

    // Shift the upper wavelengths to the appropriate redshift
    for(int i=0; i<readlamU.size(); i++)
        wlU.push_back(readlamU[i]*(1.0+redshift_)/(1.0+highZ_));

    // Perform the two interpolations to make the upper and lower data have the same grid locations
    vector<double> upperWL, lowerWL, upperT, lowerT;
    interpVectors1024(lowerWL, lowerT, wlL, readmuL, 3000., (1.+lowZ_)*1250.+0.1);
    interpVectors1024(upperWL, upperT, wlU, readmuU, 3000., (1.+lowZ_)*1250.+0.1);

    // Interpolate the new values as a weighted mean between the upper and lower curves
    combineCurves(lowerT, upperT);
    wavelengthMean_ = lowerWL;

    return;
};


void IGMVariance::interpVectors1024(vector<double> &newWL, vector<double> &newT, vector<double> &oldWL, vector<double> &oldT, double lamMin, double lamMax)
{
    int N = 1024;
    double dl = (lamMax-lamMin)/(N-1);

    // Construct newWL
    for(int i=0; i<N; i++)
        newWL.push_back(lamMin + dl*i);

    // Perform the interpolation
    double T0 = 0.0, T1 = 0.0, T = 0.0, lam0 = 0.0, lam1 = 0.0, lam = 0.0;
    for(int i=0; i<newWL.size(); i++) {
        lam = newWL[i];

        // Find the point before lam
        for(int j=0; j<oldWL.size(); j++) {
            if(oldWL[j] < lam) {
                lam0 = oldWL[j];
                T0 = oldT[j];
            }
            if(oldWL[j] > lam)
                break;
        }

        // Find the point after lam
        for(int j=0; j<oldWL.size(); j++) {
            if(oldWL[j] > lam) {
                lam1 = oldWL[j];
                T1 = oldT[j];
                break;
            }
        }

        // Calculate the new transmission value
        T = T0 + ((lam-lam0)*(T1-T0))/(lam1-lam0);
        newT.push_back(T);
    }

    return;
};

void IGMVariance::combineCurves(vector<double> &lowerT, vector<double> &upperT)
{
    double Tmid = 0.0, TL = 0.0, TU = 0.0;

    if(lowerT.size() != upperT.size())
        throw ParmError("Interpolated vectors are not the same length");

    for(int i=0; i<lowerT.size(); i++) {
        TL = lowerT[i];
        TU = upperT[i];

        Tmid = TL*(1.0 - ((redshift_ - lowZ_)/0.1)) + TU*(1.0 - ((highZ_ - redshift_)/0.1));
        mean_.push_back(Tmid);
    }

    return;
};


void IGMVariance::testPrint(vector<double> &lam, vector<double> &data)
{
    ofstream ofile;
    ofile.open("./testfiles/test.test");

    if(lam.size() != data.size())
        cout << "Trying to print 2 vectors of different sizes" << endl;
    for(int i=0; i<lam.size(); i++)
        ofile << lam[i] << "\t" << data[i] << endl;

    ofile.close();

    return;
};


