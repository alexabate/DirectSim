// -*- LSST-C++ -*-
#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

// sophya libraries
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

// my codes
#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "simdata.h"
#include "matrix.h"

// root
#include "TMinuit.h"
#include "TGraph.h"
#include "TFile.h"
#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TPrincipal.h"

#define PI 3.141592
/*



*/
void usage(void);
void usage(void)
{
	cout << endl<<" Usage: projectTemplates [...options...]" << endl<<endl;
	
	cout << " Reads in spectra listed in file INFILE "<<endl;
	cout << " Reads in eigenvectors to project spectra onto from EIGVFILE\n";
	cout << " Projects spectra onto NEIG of the eigenvectors"<<endl;
	cout << endl;
	cout << " Outputs eigenvalues of the reconstructed spectra to \n";
	cout << " OUTFILEROOT_eigenvals.txt "<<endl;
	cout << endl;
	cout << " Optionally can add reddening or interpolation between the templates";
	cout << endl << endl;
	
		
	cout << " -i : INFILE: filename of template spectra list "<<endl;
	cout << " -v : EIGVFILE: filename of eigenvectors "<<endl;
	cout << " -o : OUTFILEROOT: base output file name"<<endl;
	cout << " -n : NEIG: Number of eigenvectors to project onto [DEFAULT=6]"<<endl;
	cout << " -r : NSTEP,RMAX add reddening up to value of RMAX, in NSTEP\n";
	cout << "      steps [DEFAULT=NO]"<<endl;
	cout << " -p : NINTERP interpolate NINTERP templates linearly between\n";
	cout << "      the templates [DEFAULT=0]"<<endl;
	cout << endl;
}
int main(int narg, char* arg[])
{
	cout << " ==== projectTemplates.cc program , to calculate PC's of a";
	cout << " template set  ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;

	//--- decoding command line arguments 
	string inFile="CWWK.list";
	string eigVectFile;
	string outFileRoot="testfiles/test";
	int nEigKept=6;
	int nStepRed=0;
	double redMax=0.3;
	int nInterp=0;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hi:v:o:n:r:p:")) != -1) 
	{
	switch (c) 
		{
	  case 'i' :
		inFile = optarg;
		break;
	  case 'v' :
	    eigVectFile = optarg;
	    break;
	  case 'o' :
		outFileRoot = optarg;
		break;
	  case 'n' :
		sscanf(optarg,"%d",&nEigKept);
		break;
	  case 'r' :
		sscanf(optarg,"%d,%lf",&nStepRed,&redMax);
		break;
	  case 'p' :
		sscanf(optarg,"%d",&nInterp);
		break;
	  case 'h' :
		default :
		usage(); return -1;
		}
	}
  //-- end command line arguments
  cout <<"     Reading in templates listed in file "<<inFile<<endl;
  cout <<"     Reading eigenvectors from file "<<eigVectFile<<endl;
  cout <<"     Outputs go to file beginning: "<<outFileRoot<<endl;
  cout <<"     Projecting onto "<<nEigKept<<" eigenvectors"<<endl;
  if (nStepRed>0){
    cout <<"     Applying reddening to templates up to 0<E(B-V)<"<<redMax;
    cout <<"     in "<<nStepRed<<" steps"<<endl;
    }
  else
    cout <<"     Not applying reddening to templates"<<endl;
  if (nInterp>0)
    cout <<"     Interpolating between templates "<<nInterp<<" times"<<endl;
  cout << endl;
  
  int rc = 1;  
  try {  // exception handling try bloc at top level
	ResourceUsage res;
	InitTim();
	
	// Input/output streams
	ifstream ifs;
	ofstream outp;
	string outFile;
	
    // Read in SEDs from file
	ReadSedList readSedList(inFile,1);
	int nsed=readSedList.getNSed();
	cout <<"     Number of original SEDs = "<<nsed<<endl;
	cout << endl;
	
	// Wavelength range and resolution of SEDs
	double lmin=1e-7, lmax=1e-6;
	int nl=1500;
	double dl=(lmax-lmin)/(nl-1);
	
	// Read out SEDs into array
    readSedList.readSeds(lmin,lmax);
    cout << endl;
        
    // Interpolate?
    if (nInterp>0) { 
        readSedList.interpSeds(nInterp);
        cout <<"     Interpolated SEDs"<<endl;
        cout <<endl;
        }
    
    
    // Add reddening?
    if (nStepRed>0) { 
        readSedList.reddenSeds(nStepRed,redMax);
        cout <<"     Reddened SEDs"<<endl;
        cout <<endl;
        }
    
    //cout <<"     Write spectra to a file"<<endl;
    //outFile = outFileRoot+"_spectra.txt";
    //readSedList.writeSpectra(outFile,lmin,lmax);
    //cout << endl;
    
    // Get total number of SEDs
    vector<SED*> sedArray=readSedList.getSedArray();
    int ntot=readSedList.getNTot();
    
    // Project spectra
    TemplatePCA templatePca(sedArray,eigVectFile,lmin,lmax,nl);
    templatePca.reconstructSpectra(nEigKept);
    TMatrix<double> eigenvalsReconstruct=templatePca.returnEigValsProjSpec();

    // Write eigenvalues of reconstructed spectra to a file
	cout <<"     Write eigenvalues of spectra to a file"<<endl;
	outFile = outFileRoot+"_eigenvals.txt";
	templatePca.writeEigenValsOfProjSpec(outFile);

	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " projectTemplates.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " projectTemplates.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " projectTemplates.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of projectTemplates.cc program  Rc= " << rc << endl;
  return rc;	
}
