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
void usage(void){
	cout << endl<<" Usage: pcaTemplates [...options...]" << endl<<endl;
	
	cout << " Reads in spectra listed in file INFILE. "<<endl;
	cout << " Finds the eigenvalues and eigenvectors of the set of spectra.\n";
	cout << " Reconstructs the spectra using only the first NEIG eigenvalues\n";
	cout << " or each of 1:NEIGMAX eigenvalues.\n";
	cout << endl;
	cout << " Optionally the spectra can be interpolated and/or reddened."<<endl;
	cout << endl;
	
	cout << " Outputs:"<<endl;
	cout << "  i) lots of useful things for checking "<<endl;
	cout << " ii) the reconstructed spectra"<<endl;
	cout << "iii) the eigenvalues of the reconstructed spectra"<<endl;
	cout << endl;
		
	cout << " -i : INFILE: filename of template list (stored in $SEDLOC)"<<endl;
	cout << " -o : OUTFILEROOT: base output file name"<<endl;
	cout << " -n : NEIG: number of eigenvalues to keep [DEFAULT=6]"<<endl;
	cout << " -m : NEIGMAX: Keep 1:NEIGMAX eigenvalues via looping [DEFAULT=NO]";
	cout << endl;
	cout << " -r : NSTEP,RMAX add reddening up to value of RMAX, in NSTEP\n";
	cout << "      steps [DEFAULT=NO]"<<endl;
	cout << " -p : NINTERP interpolate NINTERP templates linearly between\n";
	cout << "      the templates [DEFAULT=0]"<<endl;
	cout << endl;
};

int main(int narg, char* arg[]) {
	cout << " ==== pcaTemplates.cc program , to calculate PC's of a";
	cout << " template set  ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	cout<<endl<<endl;

	//--- decoding command line arguments 
	string inFile="CWWK.list";
	string outFileRoot="testfiles/test";
	bool isLoop=false;
	int nEigKept=6;
	int nEigKeptMax=0;
	int nStepRed=0;
	double redMax=0.3;
	int nInterp=0;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hi:o:n:m:r:p:")) != -1) {
	switch (c) {
	  case 'i' :
		inFile = optarg;
		break;
	  case 'o' :
		outFileRoot = optarg;
		break;
	  case 'n' :
		sscanf(optarg,"%d",&nEigKept);
		break;
	  case 'm' :
		sscanf(optarg,"%d",&nEigKeptMax);
		isLoop=true;
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
  cout <<"     Outputs go to file beginning: "<<outFileRoot<<endl;
  if (isLoop)
    cout <<"     Keeping 1:"<<nEigKeptMax<<" eigenvalues"<<endl;
  else
    cout <<"     Keeping "<<nEigKept<<" eigenvalues"<<endl;
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
	ReadSedList readSedList(inFile);
	int nsed=readSedList.getNSed();
	cout <<"     Number of original SEDs = "<<nsed<<endl;
	cout << endl;
	
	// Wavelength range and resolution of SEDs
	double lmin=1e-7, lmax=1e-6;
	int nl=1500;
	//double dl=(lmax-lmin)/(nl-1);
	
	// Read out SEDs into array
    readSedList.readSeds(lmin,lmax);
    
        
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
    
    cout <<"     Write spectra to a file"<<endl;
    outFile = outFileRoot+"_spectra.txt";
    readSedList.writeSpectra(outFile,lmin,lmax);
    cout << endl;
    
    
    // Get total number of SEDs
    vector<SED*> sedArray=readSedList.getSedArray();
    //int ntot=readSedList.getNTot();
    
	
	// Perform PCA on SED array
	cout <<"     Perform PCA on spectra"<<endl<<endl;
	TemplatePCA templatePca(sedArray,lmin,lmax,nl);
	
	// Write out loads of things (useful for checking)
	// write data normalization
	outFile = outFileRoot+"_normvals.txt";
	cout <<"     Writing data normalization values to "<<outFile<<endl;
	templatePca.writeNormValues(outFile);
	// write data mean values
	outFile = outFileRoot+"_meanvals.txt";
	cout <<"     Writing data mean values to "<<outFile<<endl;
	templatePca.writeMeanValues(outFile);
	// write data matrix
	outFile = outFileRoot+"_datamatrix.txt";
	cout <<"     Writing data matrix to "<<outFile<<endl;
	templatePca.writeDataMatrix(outFile);
	// write cov matrix
	outFile = outFileRoot+"_covmat.txt";
	cout <<"     Writing covariance matrix and mean values to "<<outFile<<endl;
	templatePca.writeCovMatrix(outFile);// also writes mean values of data
	// write eigenvectors
	outFile = outFileRoot+"_eigenvectors.txt";
	cout <<"     Writing eigenvectors to "<<outFile<<endl;
	TMatrix<double> eigv=templatePca.getEigenVectors();
	//cout <<"     Size of eigenvector matrix: "<<eigv.SizeX()<<"x"<<eigv.SizeY()<<endl;
	templatePca.writeEigenVectors(outFile);
	// write eigenvalues
	outFile = outFileRoot+"_eigenvalues.txt";
	cout <<"     Writing eigenvectors to "<<outFile<<endl;
	templatePca.writeEigenValues(outFile);
	cout << endl;
	
	cout <<"     Reconstruct spectra .... "<<endl;
	cout << endl;
	
	if (isLoop) { // if chose option m, reconstructs spectra after keeping
	              // between 1 and nEigKeptMax eigenvalues
	    cout <<"     Looping over number of eigenvalues to keep ... "<<endl;
	
	    // File to put mean fit to spectrum in
	    string outFileFit = outFileRoot+"_meanfit.txt";
        ifs.open(outFileFit.c_str(),ifstream::in);
	    ifs.close();
	    if(ifs.fail()) {
	    
		    ifs.clear(ios::failbit);
		    outp.open(outFileFit.c_str(),ofstream::out);
	
	        for (int i=0; i<nEigKeptMax; i++) {
	        
	            // reconstruct spectra using nek eigenvalues
	            int nek=i+1;
	            cout <<"     Reconstructing spectra with only "<<nek<<" eigenvalues\n";
	            cout << endl;
	            templatePca.reconstructSpectra(nek);
	                
	            // Write reconstructed spectra to a file
	            // each row is a different wavelength, each column a different spectrum
	            stringstream ss;
	            ss<<nek;
	            cout <<"     Write reconstructed spectrum to a file"<<endl;
	            outFile = outFileRoot+"_recspecEigMax"+ss.str()+".txt";
	            templatePca.writeRecSpec(outFile);

                // Write eigenvalues of reconstructed spectra to a file
                // each row is a different eigenvalue, each column a different spectrum
                // the number of rows will match the value of nek
	            cout <<"     Write eigenvalues of reconstructed spectrum to a file"<<endl;
	            outFile = outFileRoot+"_receigvEigMax"+ss.str()+".txt";
	            templatePca.writeEigenValsOfProjSpec(outFile);
	            
	            // Fit the rec-spectra to find best number of eigenvalues to keep
	            TVector<double> fits=templatePca.fitSpectra();
	            double mean=0; 
	            for (int j=0; j<fits.Size(); j++)
	                mean+=fits(j);
	            mean/=fits.Size();
	            outp << nek<<"  "<<mean<<endl;

	            }
	         outp.close();
            }
        else
		    cout <<"     ERROR! file "<<outFileFit<<" exists"<<endl;
	    cout << endl;
	    }
	else { // If chose option n or neither n nor m.
	       // Only keep nEigKept eigenvalues   
	            
	    cout <<"     Reconstructing spectra with only "<<nEigKept<<" eigenvalues\n";
	    cout << endl;
	    templatePca.reconstructSpectra(nEigKept);

        // Write reconstructed spectra to a file
        // each row is a different wavelength, each column a different spectrum
        stringstream ss;
	    ss<<nEigKept;
	    cout <<"     Write reconstructed spectrum to a file"<<endl;
	    outFile = outFileRoot+"_recspec_"+ss.str()+"Eig.txt";
	    templatePca.writeRecSpec(outFile);

        // Write eigenvalues of reconstructed spectra to a file
        // each row is a different eigenvalue, each column a different spectrum
        // the number of rows will match the value of nEigKept
	    cout <<"     Write eigenvalues of reconstructed spectrum to a file"<<endl;
	    outFile = outFileRoot+"_receigv_"+ss.str()+"Eig.txt";
	    templatePca.writeEigenValsOfProjSpec(outFile);
	        
	    }

	}  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " pcaTemplates.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " pcaTemplates.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " pcaTemplates.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of pcaTemplates.cc program  Rc= " << rc << endl;
  return rc;	
}
