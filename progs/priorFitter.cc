#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>

#include <typeinfo>
#include "timing.h"

// sophya classes
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

//root classes
//#include "TH2.h"
#include "TMinuit.h"
#include "TGraph.h"
#include "TFile.h"

// my classes
#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "simdata.h"
#define PI 3.141592

// Global variables needed for the minimization function
double mLow, zLow, dM, dZ;
int nM, nZ;
TMatrixD magZData;

// Minimization function
void fcn(int &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {


    Double_t mc = 18;

    // the fit parameters
    Double_t A = par[0]; // nuisance parameter
    Double_t alpha = par[1];
    Double_t z0 = fabs(par[2]);// always keep positive
    Double_t km = par[3];

    Double_t sumsq = 0.;
    for (int im=0; im<nM; im++) {
        Double_t mv = mLow + dM*im;
        
        for (int iz=0; iz<nZ; iz++) {
            Double_t zv = zLow + dZ*iz;
                    
            // data
            Double_t data = magZData(im,iz);
            Double_t error2 = data; // Poisson error on bin content
            
            
            // model
            Double_t zm = z0 + km*(mv-mc);
	        Double_t prior1;
	        if (zv == 0) prior1 = 0;
	        else prior1 = pow(zv, alpha);

	        Double_t arg1;
	        if (zm != 0) arg1 = zv/zm;
	        else arg1 = 0;

            Double_t prior2 = exp(-pow(arg1,alpha));

            Double_t prior = A*prior1*prior2;

	        if (iflag == 3) {
	            // fill histogram values
                }

            if (error2 > 0) 
	            sumsq += (data-prior)*(data-prior)/error2;
	        else
	            sumsq += (data-prior)*(data-prior);
    
            }
    }
    f = sumsq;
};


void usage(void);
void usage(void) {
	cout << endl<<" Usage: priorFitter [...options...]" << endl<<endl;
		
	cout << "  Fit a prior using a spectroscopic sub-sample of galaxies"<< endl;
	cout <<endl;


	cout << " -t : INFILE: input filename (if text format) "<<endl;
	cout << " -f : INFILE: input filename (if FITS format) "<<endl;
	cout << " -o : OUTFILEROOT: output filename for prior fit parameter values etc"<<endl;
	cout << " -z : ZSCOL: column number of INFILE containing the redshifts [DEFAULT ZSCOL=0]"<<endl;
	cout << " -m : MAGCOL: column number of INFILE containing the magnitudes [DEFAULT MAGCOL=4]"<<endl;
	cout << "      (note column numbers will be zero indexed) "<<endl;  
	cout << " -l : ZL,ML: set the lower limits of the redshift (ZL) and "<<endl;
	cout << "             magnitude (ML) bins [DEFAULT ZL=0,ML=10]"<<endl;
	cout << " -u : ZU,MU: set the upper limits of the redshift (ZU) and "<<endl;
	cout << "             magnitude (MU) bins [DEFAULT ZU=1,MU=30]"<<endl;
	cout << " -n : NZ,MZ: set the number of redshift (NZ) and magnitude (MZ) "<<endl;
	cout <<"              bins [DEFAULT NZ=20,MZ=20]"<<endl;
	cout << endl;
    }
int main(int narg, char* arg[]) {

	cout << " ==== priorFitter.cc program , to fit a prior function to a ";
	cout << " data set  ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;

	//--- decoding command line arguments 
	string outfileroot="output/priorfit";
	string infile;
	int mCol=4,zCol=0;
	bool isText=false;
	bool isFITS=false;
	double zLower=0,mLower=18;
	double zUpper=4,mUpper=24;
	int nZbins=50,nMbins=50;
  
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ht:f:o:z:m:l:u:n:")) != -1) {
	    switch (c) {
	    
            case 't' :
		        infile = optarg;
		        isText = true;
		        break;
		    case 'f' :
		        infile = optarg;
		        isFITS = true;
		        break;
	        case 'o' :
		        outfileroot = optarg;
		        outfileroot = "output/"+outfileroot;
		        break;
	        case 'z' :
		        sscanf(optarg,"%d",&zCol);
		        break;
		    case 'm' :
		        sscanf(optarg,"%d",&mCol);
		        break;
		    case 'l' :
		        sscanf(optarg,"%lf,%lf",&zLower,&mLower);
		        break;
		    case 'u' :
		        sscanf(optarg,"%lf,%lf",&zUpper,&mUpper);
		        break;
		    case 'n' :
		        sscanf(optarg,"%d,%d",&nZbins,&nMbins);
		        break;
	        case 'h' :
		        default :
		        usage(); return -1;
		    }
	    }
	  if (isText&&isFITS)
	    throw ParmError("ERROR! Cannot give both a text and FITS file");
	  if (!isText&&!isFITS)
	    throw ParmError("ERROR! Must give either a text or FITS file");
	  if (isText)
	    throw ParmError("ERROR! Text format not supported yet soz");

      cout << "     Input catalog file is "<<infile; 
      if (isText)
        cout <<" (text format)"<<endl;
      if (isFITS)
        cout <<" (FITS format)"<<endl;
      cout << "     Output filename root is "<<outfileroot<<endl;
      cout << "     Redshifts to be read from column "<<zCol<<endl;
      cout << "     Magnitudes to be read from column "<<mCol<<endl;
      cout << "     Binning data into "<< nZbins <<" bins over "<<zLower<<"<z<"<<zUpper;
      cout << " and "<< nMbins <<" bins over "<<mLower<<"<m<"<<mUpper<<endl;
      cout << " ==== end decoding command line arguments ===="<<endl;
      cout << endl;
  //-- end command line arguments
  
  int rc = 1;  
  try {  // exception handling try bloc at top level
	ResourceUsage res;
    InitTim();
    
    ifstream inp;
	ofstream outp;
	string outfile;
	

	// Read in observed magnitudes and true redshifts from FITS file
	FitsInOutFile fin(infile,FitsInOutFile::Fits_RO);
	fin.MoveAbsToHDU(2);
	SwFitsDataTable dt(fin,512,false);
	sa_size_t ng=dt.NEntry();
	sa_size_t nc=dt.NCols();
	DataTableRow row=dt.EmptyRow();
	cout <<"     In file "<<infile<<" ... "<<endl;
	cout <<"     Number of columns = "<<nc<<", number of entries = "<<ng;
	cout << endl;
	cout <<"     Columns in the file are:"<<endl;
	cout <<"     #    Name "<<endl;
	for (int i=0; i<nc; i++)
	    cout << "     "<<i<<"    "<<dt.NomIndex(i)<<endl;
	cout << endl;
    double mMin,mMax;
    dt.GetMinMax(mCol,mMin,mMax);
    cout <<"     Range of "<<dt.NomIndex(mCol)<<" column is "<<mMin<<"<m<"<<mMax<<endl;
    double zMin,zMax;
    dt.GetMinMax(zCol,zMin,zMax);
    cout <<"     Range of "<<dt.NomIndex(zCol)<<" column is "<<zMin<<"<z<"<<zMax<<endl;
    cout << endl;
    
    // Create 2D histogram of i band magnitude and redshifts
    cout <<"     Creating 2D histogram of i band magnitude and redshifts ... "<<endl;
    
    // Define bins
    double dZbin = (zUpper - zLower)/(nZbins - 1);
    double dMbin = (mUpper - mLower)/(nMbins - 1);
    
    TArray<int> magRedshiftData;
    int nDim = 2;
    sa_size_t myDim[nDim];
    myDim[0]=nMbins; myDim[1]=nZbins;
    magRedshiftData.SetSize(nDim,myDim);
    
    for (long ig=0; ig<ng; ig++) {
        dt.GetRow(ig,row);
        double mag = row[mCol];
        double zs = row[zCol];
        
        sa_size_t mIndex = (sa_size_t)floor((mag - mLower)/dMbin);
        sa_size_t zIndex = (sa_size_t)floor((zs - zLower)/dZbin);
        
        if ( mIndex<nMbins && zIndex<nZbins && mIndex>0 && zIndex>0)// if within bins
             magRedshiftData(mIndex,zIndex)++;
             
        }
        
    outfile = outfileroot + "_magdist.txt";
    outp.open(outfile.c_str());
    cout <<"     Writing histogram to file "<< outfile << endl;
    for (int i=0; i<nMbins; i++){
        for (int j=0; j<nZbins; j++)
            outp << magRedshiftData(i,j) << "  ";
        outp << endl;
        }
    outp.close();
    cout << endl;
    

    // Brut force grid method for parameter estimation

    cout <<"     Estimating parameters using grid method "<<endl;
    
    // Trial parameter grids    
    cout <<"     Setting parameter grids "<<endl;
    // alpha
    int nalpha = 40;
    double amin = 0.2, amax = 4.;
    double da = (amax - amin)/(nalpha - 1);
    cout <<"     "<< amin <<" < alpha < "<< amax <<" with step size = "<< da <<endl;

    // k
    int nk = 20;
    double kmin = 0.01, kmax=0.2;
    double dk = (kmax - kmin)/(nk - 1);
    cout <<"     "<< kmin <<" < k < "<< kmax <<" with step size = "<< dk <<endl;
    
    // z0
    int nz0 = 50;
    double z0min = 0.01, z0max = 0.5;
    double dz0 = (z0max - z0min)/(nz0 - 1);
    cout <<"     "<< z0min <<" < z0 < "<< z0max <<" with step size = "<< dz0 <<endl;
    cout << endl;
    
    cout <<"     Beginning loop over parameters ...."<<endl;

    double mc = 18; // maybe should be 18
    TArray<double> chisq, prob;
    int ndim = 3;
    sa_size_t mydim[ndim];
    mydim[0] = nalpha; mydim[1] = nz0; mydim[2] = nk;
    chisq.SetSize(ndim,mydim);
    prob.SetSize(ndim,mydim);
    // loop over the trial parameters
    for (int ia = 0; ia<nalpha; ia++) {
        cout << "     Global loop "<< ia+1 <<" of "<< nalpha <<endl;
        double alpha = amin + da*ia;
        
        for (int iz0=0; iz0<nz0; iz0++) {
            double z0 = z0min + dz0*iz0;
            
            for (int ik=0; ik<nk; ik++) {
                double k = kmin + dk*ik;
            
                //vector<double> diff;
                // loop over each magnitude bin
                double sumsq = 0;
                for (int im=0; im<nMbins; im++) {
                    double mv = mLower + dMbin*im;
                    
                    vector<double> data, prior;
                    for (int iz=0; iz<nZbins; iz++) {
                        double zv = zLower + dZbin*iz;
                    
                        data.push_back(magRedshiftData(im,iz));
                    
                        double zm = z0+k*(mv-mc);
                        double pr = pow(zv,alpha)*exp(-pow(zv/zm,alpha));
                        prior.push_back(pr);
                        

                        
                        }
                    int tmp;
                    double maxData = findMaximum(data,tmp);
                    double maxPrior = findMaximum(prior,tmp);
                    //cout << "max values: "<<maxData <<"  "<<maxPrior <<endl;
                    
                    for (int iz=0; iz<nZbins; iz++) {
                       
                        double d = data[iz], p = prior[iz];
                        if (d>0)
                            d/=maxData;
                        if (p>0)
                            p/=maxPrior;  
                        
                        
                        double dmp = d - p;
                        int nantest=my_isnan(dmp);
                        double zm = z0+k*(mv-mc);
                        if (nantest>0)
                            cout <<" is nan!: mv = "<< mv <<", data[iz] = "<< data[iz] <<", prior[iz] = "<< prior[iz] <<endl;
                        
                        sumsq += dmp*dmp;
                        }
                        
                        
                    }
                chisq(ia,iz0,ik) = sumsq;
                prob(ia,iz0,ik) = exp(-sumsq/2);
                }
            }
        }
    cout << endl;
    
    cout <<"     Marginalize over the parameters ... "<<endl;
    cout << endl;
    
    cout <<"     ... k "<<endl;
    vector<double> probk;
    for (int ik=0; ik<nk; ik++) {
        
        double sum = 0;
        for (int ia = 0; ia<nalpha; ia++)
            for (int iz0=0; iz0<nz0; iz0++)
                sum += prob(ia,iz0,ik);
        
        probk.push_back(sum);
        
        }
    outfile = outfileroot + "_margprobk.txt";
    outp.open(outfile.c_str());
    for (int ik=0; ik<nk; ik++) {
        double k = kmin + dk*ik;
        outp << k <<"  "<< probk[ik] <<endl;
        }
    outp.close();
    cout << endl;

    cout <<"     ... z0 "<<endl;
    vector<double> probz;
    for (int iz0=0; iz0<nz0; iz0++) {
        
        double sum = 0;
        for (int ia = 0; ia<nalpha; ia++)
            for (int ik=0; ik<nk; ik++)
                sum += prob(ia,iz0,ik);
        
        probz.push_back(sum);
        
        }
    outfile = outfileroot + "_margprobz.txt";
    outp.open(outfile.c_str());
    for (int iz0=0; iz0<nz0; iz0++) {
        double z0 = z0min + dz0*iz0;
        outp << z0 <<"  "<< probz[iz0] <<endl;
        }
    outp.close();
    cout << endl;
        
    cout <<"     ... alpha "<<endl;
    vector<double> proba;
    for (int ia = 0; ia<nalpha; ia++){
        
        double sum = 0;
        for (int iz0=0; iz0<nz0; iz0++) 
            for (int ik=0; ik<nk; ik++)
                sum += prob(ia,iz0,ik);
        
        proba.push_back(sum);
        
        }
    outfile = outfileroot + "_margproba.txt";
    outp.open(outfile.c_str());
    for (int ia = 0; ia<nalpha; ia++) {
        double alpha = amin + da*ia;
        outp << alpha <<"  "<< proba[ia] <<endl;
        }
    outp.close();
    cout << endl;
    
    int iMax;
    findMaximum(probk,iMax);
    double kbest = kmin + dk*iMax;
    findMaximum(proba,iMax);
    double abest = amin + da*iMax;
    findMaximum(probz,iMax);
    double z0best = z0min + dz0*iMax;
    cout <<"     Results from grid method:"<<endl;
    cout <<"  NO.   NAME      VALUE "<<endl;
    cout <<"   1    alpha     "<< abest <<endl;
    cout <<"   2    z0        "<< z0best <<endl;
    cout <<"   3    km        "<< kbest <<endl;
    cout << endl;
    
    // Now do proper minimization technique
    
    cout <<"     Now estimating parameters using MINUIT "<<endl;
    // Set the global variables needed by the minimization function
    cout <<"     Setting global variables needed by minimization function "<<endl;
    mLow = mLower;
    zLow = zLower;
    dM = dMbin;
    dZ = dZbin;
    nM = nMbins;
    nZ = nZbins;
    magZData.ResizeTo(nMbins,nZbins);
    // convert to ROOT object
    for (int im=0; im<nMbins; im++)
        for (int iz=0; iz<nZbins; iz++) 
	        magZData(im,iz)=magRedshiftData(im,iz);
	cout << endl;
	        
	cout <<"     Initializing MINUIT "<<endl;
	//initialize TMinuit with a maximum of 4 params
    TMinuit *gMinuit = new TMinuit(4);  
    gMinuit->SetFCN(fcn);
    cout << endl << endl;
    
    // The mnexcm method "Minuit execute command" executes whatever command is
    // given in the first argument.
    // The 2nd argument to mnexcm supplies the actual arguments to give to the 
    // command given in the first argument 
    // The 3rd argument specifies how many arguments are supplied in that 2nd
    // argument
    
    
    Double_t arglist[10]; // for use with mnexcm
    Int_t ierflg = 0; // error flag
    // This is the list of parameters for the fit
    
    
    // Do a chi2 fit
    // setting the SET ERR command to equal 1 means you are doing a chi-square
    // minimization.  Probably this is related to 1-sig = dchisq = 1
    // For likelihood this number should be 0.5
    arglist[0] = 1;
    cout <<"     Set chi-square fit ... "<<endl;
    gMinuit->mnexcm("SET ERR", arglist, 1, ierflg);
    cout << endl << endl;
    
    // Set starting and step values for the parameters
    cout <<"     Set starting parameters and step values "<<endl;
    gMinuit->mnparm(0, "A", 1e5, 1., 0, 1e10, ierflg);
    gMinuit->mnparm(1, "alpha", 1.3, 0.1, 0, 0, ierflg);
    gMinuit->mnparm(2, "z0", 0.1, 0.1, 0, 0, ierflg);
    gMinuit->mnparm(3, "km", 0.1, 0.001, 0, 5, ierflg);
    cout << endl << endl;
    
    
    
    // Minimize the chi2
    
    arglist[0] = 500; // max number of calls
    arglist[1] = 1.;  // tolerance
    
    cout <<"     Minimizing chi-squared with MIGRAD"<<endl;
    double AMigrad, alphaMigrad, z0Migrad, kmMigrad;
    double dAMigrad, dalphaMigrad, dz0Migrad, dkmMigrad;
    gMinuit->mnexcm("MIGRAD", arglist, 2, ierflg);
    gMinuit->GetParameter(0,AMigrad,dAMigrad );
    gMinuit->GetParameter(1,alphaMigrad,dalphaMigrad);
    gMinuit->GetParameter(2,z0Migrad,dz0Migrad);
    gMinuit->GetParameter(3,kmMigrad,dkmMigrad);
    cout << endl << endl;
    
    cout <<"     Minimizing chi-squared with MINI"<<endl;
    double AMini, alphaMini, z0Mini, kmMini;
    double dAMini, dalphaMini, dz0Mini, dkmMini;
    gMinuit->mnexcm("MINI",   arglist, 2, ierflg);
    gMinuit->GetParameter(0,AMini,dAMini );
    gMinuit->GetParameter(1,alphaMini,dalphaMini);
    gMinuit->GetParameter(2,z0Mini,dz0Mini);
    gMinuit->GetParameter(3,kmMini,dkmMini);
    cout << endl << endl;
    
    cout <<"     Minimizing chi-squared with MINOS"<<endl;
    double AMinos, alphaMinos, z0Minos, kmMinos;
    double dAMinos, dalphaMinos, dz0Minos, dkmMinos;
    gMinuit->mnexcm("MINOS",  arglist, 2, ierflg);
    gMinuit->GetParameter(0,AMinos,dAMinos );
    gMinuit->GetParameter(1,alphaMinos,dalphaMinos);
    gMinuit->GetParameter(2,z0Minos,dz0Minos);
    gMinuit->GetParameter(3,kmMinos,dkmMinos);
    cout << endl << endl;
    
    // Fix parameters
    cout <<"     Minimizing chi-squared with MIGRAD but z0 is fixed to 0.1"<<endl;
    TMinuit *gMinuitFix = new TMinuit(4);  
    gMinuitFix->SetFCN(fcn);
    arglist[0] = 1;
    gMinuitFix->mnexcm("SET ERR", arglist, 1, ierflg);
    gMinuitFix->mnparm(0, "A", 1e5, 1., 0, 1e10, ierflg);
    gMinuitFix->mnparm(1, "alpha", 1.3, 0.1, 0, 0, ierflg);
    gMinuitFix->mnparm(2, "z0", 0.1, 0.1, 0, 0, ierflg);
    gMinuitFix->mnparm(3, "km", 0.1, 0.001, 0, 5, ierflg);
    gMinuitFix->FixParameter(2);
    arglist[0] = 500; // max number of calls
    arglist[1] = 1.;  // tolerance
    double AMigradFix, alphaMigradFix, z0MigradFix, kmMigradFix;
    double dAMigradFix, dalphaMigradFix, dz0MigradFix, dkmMigradFix;
    gMinuitFix->mnexcm("MIGRAD", arglist, 2, ierflg);
    gMinuitFix->GetParameter(0,AMigradFix,dAMigradFix );
    gMinuitFix->GetParameter(1,alphaMigradFix,dalphaMigradFix);
    gMinuitFix->GetParameter(2,z0MigradFix,dz0MigradFix);
    gMinuitFix->GetParameter(3,kmMigradFix,dkmMigradFix);
    cout << endl << endl;

    // Print results
    cout <<"    Printing results .... "<<endl;
    cout <<"            MIGRAD       MINI         MINOS "<<endl;
    cout <<"         A: "<< AMigrad <<"       "<< AMini <<"       "<< AMinos <<endl;
    cout <<"     alpha: "<< alphaMigrad <<"      "<< alphaMini <<"      "<< alphaMinos <<endl;
    cout <<"        z0: "<< z0Migrad <<"  "<< z0Mini <<"  "<< z0Minos <<endl;
    cout <<"        km: "<< kmMigrad <<"     "<< kmMini <<"     "<< kmMinos <<endl;
    cout << endl << endl;
    
    cout <<"    Printing results fixed z0 .... "<<endl;
    cout <<"            MIGRAD     "<<endl;
    cout <<"         A: "<< AMigradFix     << endl;
    cout <<"     alpha: "<< alphaMigradFix << endl;
    cout <<"        z0: "<< z0MigradFix << endl;
    cout <<"        km: "<< kmMigradFix << endl;
    
    /*Double_t chi2,edm,errdef;
    Int_t nvpar,nparx,istat;
    
    gMinuit->mnstat(chi2,edm,errdef,nvpar,nparx,istat);
    
    cout <<"    Printing values, errors and limits of parameters "<<endl;
    gMinuit->mnprin(1,chi2); // prints values, errors, limits of parameters when fcn=chi2
    
    Double_t a0, da0;
    gMinuit->GetParameter(0, a0, da0);*/

    if (ierflg) {
      Printf("Error in minimization.");
      //return ierflg;
    }


	
	
	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " priorFitter.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " priorFitter.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " priorFitter.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of priorFitter.cc program  Rc= " << rc << endl;
  return rc;	
}
