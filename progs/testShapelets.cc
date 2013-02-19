// testing branching: can i just index the change?
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

#include "hpoly.h"
#include "genericfunc.h"
#include "constcosmo.h"
#include "cosmocalcs.h"
#include "igm.h"
#include "mydefrg.h"
#include "sedfilter.h"
#include "simdata.h"
#include "shapelets.h"


// root
//#include "root_plots.h"
//#include "TDecompSVD.h"
#include "TMatrixDEigen.h"
#include "TPrincipal.h"
#include <cmath>
//#include "iomanip.h"
#include "TRandom.h"
#define PI 3.141592
/*



*/
void usage(void);
void usage(void) {
	cout << endl<<" Usage: test [...options...]" << endl<<endl;
		
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

	cout << " ==== test.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string outfile = "testfiles/TESTCODE.txt";
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"ho:")) != -1) 
	{
	switch (c) 
		{
	  case 'o' :
	    outfile = optarg;
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
	//string outfile;
	
	outfile = "testfiles/hermitebasis.txt";

	Hermite hmite;
	BasisFuncs basisFuncs;
	
	double xmin=-3, xmax=3;
	int nx = 1000;
	double dx = (xmax-xmin)/(nx-1);
	
	outp.open(outfile.c_str(), ofstream::out);
	
	double beta = 0.1;
	int order;
	for (int i=0; i<nx; i++){
	    double x = xmin+i*dx;
	    	    
	    order=0;
	    double h0=hmite.returnHermiteN(order,x);
	    double p0=basisFuncs.phiBasisFunc(order,x);
	    double b0=basisFuncs.bBasisFunc(order,x,beta);
	    order=1;
	    double h1=hmite.returnHermiteN(order,x);
	    double p1=basisFuncs.phiBasisFunc(order,x);
	    double b1=basisFuncs.bBasisFunc(order,x,beta);
		order=2;
	    double h2=hmite.returnHermiteN(order,x);
	    double p2=basisFuncs.phiBasisFunc(order,x);
	    double b2=basisFuncs.bBasisFunc(order,x,beta);
	    order=3;
	    double h3=hmite.returnHermiteN(order,x);
	    double p3=basisFuncs.phiBasisFunc(order,x);
	    double b3=basisFuncs.bBasisFunc(order,x,beta);
	    
	    //cout << x <<"  "<< h0 <<"  "<< h1 << "  "<< h2 << "  "<< h3 <<endl;
	    outp << x  <<"  "<< h0 <<"  "<< h1 <<"  "<< h2 <<"  "<< h3 <<"  ";
	    outp << p0 <<"  "<< p1 <<"  "<< p2 <<"  "<< p3 <<"  ";
	    outp << b0 <<"  "<< b1 <<"  "<< b2 <<"  "<< b3 << endl;
	    }
	outp.close();
	
	// Read in photo-z distribution file
    ifstream infs;
    string infile = "testfiles/testDZcat.txt";
    infs.open(infile.c_str(),ifstream::in);
    TArray<r_4> photoZdist;
    sa_size_t nr,nc;
    photoZdist.ReadASCII(infs,nr,nc);
    cout <<"      ... done "<<endl;
    cout <<"     Catalog has "<< nr <<" rows and "<< nc <<" columns "<<endl;
    infs.close();
	cout << endl;
	
	vector<double> xv, yv, yv1;
	for (int i=0; i<nr; i++) {
	    xv.push_back(photoZdist(0,i));
	    yv.push_back(photoZdist(1,i));
	    yv1.push_back(1.);
	    }
	    
	
	SInterp1D func(xv,yv);
	SInterp1D func1(xmin,xmax, yv1);
	Shapelets shapelets(func);
	//Shapelets shapelets1(func1);
	int nmax = 6;
	int nmax2 = 15;
	
	outfile = "testfiles/shapelets.txt";
	outp.open(outfile.c_str(), ofstream::out);
	for (int i=0; i<nx; i++){
	    cout <<"     On "<< i+1 <<" of "<< nx <<endl;
	    double x = xmin+i*dx;
	    
	    double fx = func(x);
	    double fxr = shapelets.functionRep(x,nmax,beta);
        double fxr1 = shapelets.functionRep(x,nmax2,beta);
	    outp << x <<"  "<< fx <<"  "<< fxr <<"  "<< fxr1 << endl;
	    }
	outp.close();
	

	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " test.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " test.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " test.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of test.cc program  Rc= " << rc << endl;
  return rc;	
}

