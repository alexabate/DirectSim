
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

// sophya libraries
//#include "histinit.h"
#include "machdefs.h"
#include "sopnamsp.h"
#include "fiosinit.h"
#include "mydefrg.h"

#include "hpoly.h"
//#include "genericfunc.h"
#include "constcosmo.h"
#include "cosmocalcs.h"
#include "igm.h"
#include "mydefrg.h"
#include "sedfilter.h"
#include "simdata.h"


void usage(void);
void usage(void) {
	cout << endl<<" Usage: fitLSSTspectra [...options...]" << endl<<endl;
    cout << " -o OUTFILEROOT "<<endl;
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

    cout << " ==== fitLSSTspectra.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfileroot = "output/fitLSSTspectra";
    
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ho:")) != -1) {
	    switch (c) {
	        case 'o' :
	            outfileroot = optarg;
	            outfileroot = "output/" + outfileroot;
	            break;
	        case 'h' :
		        default :
		        usage(); return -1;
		    }
	    }
	    

    
    //-- end command line arguments
    cout << "     Output to be written to files beginning "<< outfileroot <<endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    ifstream inp;
	ofstream outp;
	string outfile;

    cout <<"     Reading in SED libraries "<<endl;	
    string sedFile;
    // Wavelength range and resolution of SEDs
	double lmin=1e-7, lmax=1e-6;
	int nl=3000;
	double dl = (lmax - lmin)/(nl -1);
	
	// Load in CWWK SEDs
	sedFile = "CWWKSB.list";
	ReadSedList readCWWK(sedFile);
    readCWWK.readSeds(lmin,lmax);
    vector<SED*> sedCWWK=readCWWK.getSedArray();
    string sedLoc = readCWWK.getSedDirEnviromentVar();
    vector<string> sedFilenames=readCWWK.returnSedFilenames();
    int nCWWK = readCWWK.getNSed();
    cout << "     "<< nCWWK <<" CWWK seds "<<endl;
    
    // Load in LSST SEDs
	sedFile = "LSST.list";
	ReadSedList readLSST(sedFile);
    readLSST.readSeds(lmin,lmax);
    vector<SED*> sedLSST=readLSST.getSedArray();
	int nLSST = readLSST.getNSed();
    cout << "     "<< nLSST <<" LSST seds "<<endl;
    cout << endl;
    
    // Find normalization values of the SEDs
    double lamNorm = 5500e-10;
    cout << "     Normalizing spectra to 1 at "<< lamNorm*1e10 <<" angstroms "<<endl;
    vector<double> normCWWK;
    for (int i=0; i<nCWWK; i++) {
        double val = sedCWWK[i]->returnFlux(lamNorm);
        normCWWK.push_back(val);
        }
    vector<double> normLSST;
    for (int i=0; i<nLSST; i++) {
        double val = sedLSST[i]->returnFlux(lamNorm);
        normLSST.push_back(val);
        }
    cout << endl;
    
    // Open file to write to
    outfile = outfileroot + "_fitvalues.txt";
    outp.open(outfile.c_str());
      
    // Fit SEDs to each other
    vector<int> bestFitSpectra;
    for (int i=0; i<nLSST; i++) { // loop over LSST SEDs
        cout <<"     On LSST sed "<< i+1 <<" of "<< nLSST;
    
        vector<double> fitValues;
        for (int j=0; j<nCWWK; j++) { // for each CWWKSB SED
        
            double sum = 0;
            for (int il=0; il<nl; il++) {
                double lam = lmin + il*dl;
                double sCWWK = (sedCWWK[j]->returnFlux(lam))/normCWWK[j];
                double sLSST = (sedLSST[i]->returnFlux(lam))/normLSST[i];
                sum += (sCWWK-sLSST)*(sCWWK-sLSST); // Chi-square type fit
                }
                
            fitValues.push_back(sum);
            }
            
        for (int j=0; j<nCWWK; j++)
            outp << fitValues[j] <<"  ";
        outp << endl;
            
        // find CWWKSB SED with smallest chi-square value (this is the best fit!)
        int iElement;
        findMinimumPosition(fitValues, iElement);
        bestFitSpectra.push_back(iElement);
        cout <<", best fit CWWK spectrum: "<< iElement <<endl;
        }
    outp.close();
    cout << endl;   
        
    double eps = 1e-6;
    TArray<double> newSEDs;
    sa_size_t ndim = 2;
    sa_size_t mydim[ndim];
    mydim[0] = nl;
    mydim[1] = nCWWK + 1;
    newSEDs.SetSize(ndim,mydim);
    
    outfile = outfileroot + "_newSEDs.txt";
    outp.open(outfile.c_str());
    for (int il=0; il<nl; il++) {
        double lam = lmin + il*dl;
        outp << lam <<"  ";
        newSEDs(il,0) = lam;
        
        int cntTot = 0;
        for (int i=0; i<nCWWK; i++) { // loop over each CWWKSB SED
            int iWant = i;
        
            int cnt = 0;
            double sumSEDvals = 0, sumSEDvalsSq = 0;
	        for (int j=0; j<nLSST; j++) {
	            int iValue = bestFitSpectra[j];
	                
	            // Find all LSST SEDs that had this CWWKSB SED as the best fit,
	            // and average them together
	            if ( abs(iWant - iValue)<eps ) {
	                double val = (sedLSST[j]->returnFlux(lam));
	                sumSEDvals += val;
	                sumSEDvalsSq += val*val;
	                cnt++;
	                }
	            }
	        //cout <<"     "<< cnt <<" LSST spectra have a best-fit CWWK = "<<iWant<<endl;
	        cntTot += cnt;
	        
	        double sedVal = sumSEDvals/cnt;
	        double sedVar = sumSEDvalsSq/cnt - sedVal*sedVal;
	        
	        outp << sedVal <<"  "<< sedVar << "  ";
	        newSEDs(il,i+1) = sedVal;
	        
	        }// end loop over CWWK SEDs
	    outp << endl;
	    
	    if (cntTot!=nLSST)
	        throw ParmError("ERROR! best-fit SED numbers don't add up");
	        
        }// end loop over wavelengths
    outp.close();  
    cout << endl;
     
    // Write these new spectra to files
    cout <<"     Writing the new spectra to files: "<<endl;
    for (int i=0; i<nCWWK; i++) {
        outfile = sedLoc+"/"+sedFilenames[i]+".lsst";
        cout <<"     "<<outfile<<endl;
        outp.open(outfile.c_str());
    
	    for (int il=0; il<nl; il++) 
	        outp << newSEDs(il,0) <<"  "<<newSEDs(il,i+1) << endl;
	    
	    outp.close();
	    }
	cout << endl;
        
    // Write these new spectra to files
    cout <<"     Writing the new spectra to files BPZ style: "<<endl;
    cout <<"     (wavelength is in angstroms and filename ends in .sed)"<<endl;
    for (int i=0; i<nCWWK; i++) {
        outfile = sedLoc+"/"+sedFilenames[i]+".lsst.sed";
        cout <<"     "<<outfile<<endl;
        outp.open(outfile.c_str());
    
	    for (int il=0; il<nl; il++) 
	        outp << newSEDs(il,0)*1e10 <<"  "<<newSEDs(il,i+1) << endl;
	    
	    outp.close();
	    }
	cout << endl;
	
    outfile = sedLoc+"/LSSTshort.list";
    cout <<"     Writing a new spectra list file "<< outfile << endl;
    outp.open(outfile.c_str());
    for (int i=0; i<nCWWK; i++)
        outp << sedFilenames[i]+".lsst" << endl;
    outp.close();
    cout << endl;
    
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " fitLSSTspectra.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " fitLSSTspectra.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " fitLSSTspectra.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of fitLSSTspectra.cc program  Rc= " << rc << endl;
  return rc;	
}

