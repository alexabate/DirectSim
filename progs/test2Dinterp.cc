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
#include "tarray.h"
#include "ctimer.h"

// my classes
#include "sinterp.h"

double func(double a, double b)
{
    //return 5.*a + b;
    //return 5.*a*a + b;
    return 5.*a*a*a + b;
};

void usage(void);
void usage(void) {
	cout << endl<<" Usage: test2Dinterp [...options...]" << endl<<endl;
		
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

	cout << " ==== test2Dinterp.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfileroot = "testfiles/test2Dinterp";
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ho:")) != -1) 
	{
	switch (c) 
		{
	  case 'o' :
	    outfileroot = optarg;
	    break;
	  case 'h' :
		default :
		usage(); return -1;
		}
	}

  //-- end command line arguments
  cout << endl;
  
  int rc = 1;  
    try {  // exception handling try bloc at top level
  
    ifstream inp;
	ofstream outp;
	string outfile;
	
	// variables xa and xb
	vector<double> xa, xb;
	int na = 10, nb = 10;
	
	// fill xa
	for (int i=0; i<na; i++)
	    xa.push_back(i);
	    
	// fill xb 
	for (int i=0; i<nb; i++)
	    xb.push_back(i);
	
	TArray<double> y;
	int ndim=2;
	sa_size_t mydim[ndim];
	mydim[0]=na;
	mydim[1]=nb;
	y.SetSize(ndim,mydim);
	
	for (int i=0; i<na; i++) {
	    double xaval = xa[i];
	    for (int j=0; j<nb; j++) {
	        double xbval = xb[j];
	        y(i,j) = func(xaval, xbval);
	        }    
	    }
	    
	for (int i=0; i<na; i++) {
	    for (int j=0; j<nb; j++) 
	        cout << y(i,j) <<"  ";
	    cout << endl;
	    }
	
	SInterp2D interpTest(xa, xb, y);
	
	Timer tm("timer",false);
	
	// Check the interpolation
	double xmin = 0.1, xmax = 8.9;
	int nTest = 100;
	double dx = (xmax-xmin)/(nTest-1);
	
	for (int i=0; i<nTest; i++)
	    for (int j=0; j<nTest; j++) {
	    
	        double x1 = xmin + i*dx;
	        double x2 = xmin + j*dx;
	
	        tm.Split();
	        double y1 = interpTest.biLinear(x1,x2);
	        tm.Split();
	        if (i<1&&j<1)
	            cout <<"Binear takes: "<< tm.PartialElapsedTimems() <<" ms"<<endl;
            tm.Split();
	        double y1o = interpTest.biLinearAccurate(x1,x2);
	        tm.Split();
	        if (i<1&&j<1)
	            cout <<"Binear-accurate takes: "<< tm.PartialElapsedTimems() <<" ms"<<endl;
	        double yTrue = func(x1, x2);
	        double perc = 100.*(y1-yTrue)/yTrue;
	        double perco = 100.*(y1o-yTrue)/yTrue;
	        cout << yTrue <<"  "<< y1 << "  "<< y1o <<"  "<< perc <<"  "<< perco << endl;
	        }
	        
	// get filename of relevant k-correction
    stringstream ss1,ss2,ss3,ss4,ss5,ss6,ss7;
    ss1 << 0; 
    ss2 << 0;
    ss3 << 0; ss4 << 3; ss5 << 2000;
    ss6 << 0.3; ss7 << 200;
    string fname;
    fname  = "kCorrections/kCorr_CWWKsed" + ss1.str();
    fname += "_LSSTfilt" + ss2.str() + "_GOODSB";
    fname += "_zmin" + ss3.str() + "_zmax" + ss4.str() + "_nz" + ss5.str();
    fname += "_emax" + ss6.str() + "_ne" + ss7.str();
    fname = fname + "_woMadau.txt";
    
    // read data from file into an array
    ifstream ifs;
	ifs.open(fname.c_str(), ifstream::in);
	if (ifs.fail()) { 
	    string emsg = "ERROR: failed to find k-correction file " + fname;
		throw ParmError(emsg);
		}
	sa_size_t nr, nc;
	TArray<double> tab;
	tab.ReadASCII(ifs,nr,nc);
	cout <<" nr = "<< nr <<", tab.SizeX() = "<< tab.SizeX() <<endl;
    cout <<" nc = "<< nc <<", tab.SizeY() = "<< tab.SizeY() <<endl;
    
    vector<double> zvals, evals;
    double zmax = 3, emax = 0.3;
    int nz = 2000, ne = 200;
    double dz = zmax/(nz-1), de=emax/(ne-1);
    for (int iz=0; iz<nz; iz++)
        zvals.push_back(dz*iz);
    for (int ie=0; ie<ne; ie++)
        evals.push_back(de*ie);
    
    
    SInterp2D interpkCorrTest(evals, zvals, tab);
    string outFile = "testfiles/tmp.txt";
    outp.open(outFile.c_str(),ofstream::out);
    for (int ie=0; ie<1; ie++)
        for (int iz=0; iz<nz; iz++) {
            double zv = 0.01 + dz*iz;
            double ev = 0.001 + de*ie;
            tm.Split();
            double kint = interpkCorrTest.biLinearAccurate(ev,zv);
            tm.Split();
	        if (iz<1&&ie<1)
	            cout <<"K correction interpolation takes: "<< tm.PartialElapsedTimems() <<" ms"<<endl;
            tm.Split();
	        double kint2 = interpkCorrTest.biLinear(ev,zv);
            tm.Split();
	        if (iz<1&&ie<1)
	            cout <<"K correction inaccurate interpolation takes: "<< tm.PartialElapsedTimems() <<" ms"<<endl;
            outp << zv <<"  "<< ev <<"  "<< kint << endl;
            }
    outp.close();
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " test2Dinterp.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " test2Dinterp.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " test2Dinterp.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of test2Dinterp.cc program  Rc= " << rc << endl;
  return rc;	
}

