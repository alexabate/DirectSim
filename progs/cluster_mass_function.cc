#include <iostream>
#include <fstream>

#include "machdefs.h"
#include "sopnamsp.h"
#include "fiosinit.h"


#include <massfunc.h>
#include <constcosmo.h>

void usage(void);
void usage(void) {

	cout << endl<<"  Usage: cluster_mass_function [...options...]   "<<endl<<endl;
	cout << "  Program to calculate the number of clusters within         "<<endl;
    cout << "  some volume at some redshift between some mass limits"<<endl<<endl;

	cout << " -z [redshift=1]                          "<<endl;
	cout << " -v [volume=11e6] (in Mpc^3)              "<<endl;
	cout << " -m [m1=5e13,m2=1e15] (in solar masses)   "<<endl;
	cout << " -o [outfile]                             "<<endl; 
	cout << endl;
};


int main(int narg, char* arg[]){

  	cout << " ==== cluster_mass_function.cc program ====" << endl;

  	// make sure SOPHYA modules are initialized 
  	SophyaInit();  
  	cout<<endl<<endl;

    string outfile;
  	double redshift = 1.;
  	double volume = 11e6; // Mpc^3
  	double m1 = 5e13; // solar masses
  	double m2 = 1e15; // solar masses

  	//--- decoding command line arguments 
    cout << " decoding command line arguments ..."<<endl;
    char c;
    while((c = getopt(narg,arg,"hz:v:m:o:")) != -1) {
	    switch (c) {
	        case 'o' :
	            outfile = optarg;
	            break;
	    	case 'z' :
                sscanf(optarg,"%lf",&redshift);
			    break;
		    case 'v' :
                sscanf(optarg,"%lf",&volume);
			    break;
			case 'm' :
                sscanf(optarg,"%lf,%lf",&m1,&m2);
			    break;
		    case 'h' :
		    default :
			    usage(); return -1;
		    }
	    }
	    
	cout <<"     Returning the expected number of clusters at redshift ";
	cout << redshift <<" in volume "<< volume <<" Mpc^3 between mass limits ";
	cout << m1 <<" and "<< m2 <<" solar masses"<<endl;
	cout <<"     Mass function will be written to file "<< outfile << endl;
	cout << endl;
	
	int rc = 1;  
  	try {  // exception handling try bloc at top level

	// Output streams
	ofstream outp;
	
	
	// cosmological parameters
	double oc = 0.253;
	double ob = 0.048;
	double om = oc + ob;
	double ol = 1. - om;
	double h = 0.7;
	double n = 1;
	double w0 = -1.;
	double sig8 = 0.9;
	
	
	// cosmological calculations
	SimpleUniverse su(h,om,ol);
	su.SetOmegaPhoton(0);
	su.SetOmegaRadiation(0);
	
	
	// power spectrum calculations
	InitialPowerLaw ipl(n);
	TransferEH tf(h, oc, ob, T_CMB_K);
	GrowthFN gro(om, ol, w0);
	PkSpecCalc pkz(ipl, tf, gro, redshift);
	
	
	// mass function calculations
	bool typeLog = false;
	MassFunc massFunc(su, pkz, redshift, sig8, typeLog);
	double nc = volume*massFunc.Integrate(m1, m2, 1000);
	cout <<"     Number of clusters in volume "<< volume <<" at z = "<< redshift;
	cout <<" between "<< m1 <<" and "<< m2 <<" solar masses is "<< nc << endl;
	
	
	// write out
	outp.open(outfile.c_str(), ofstream::out);
	
	// masses, log spaced
    int nm = 997; // ACCURACY
    double Mmin = 1.06e-7, Mmax=2.8301e16;
    double logMmin = log(Mmin), logMmax = log(Mmax);
    double dlogM = (logMmax-logMmin)/(nm-1);
    
    for (int i=0; i<nm; i++) {
    
        double logM = logMmin + i*dlogM;
        double M = exp(logM);
        double mf = massFunc(M);
        outp << M <<"  "<< mf << endl;
        
        }
    outp.close();

	
  }  // End of try bloc 
  
 
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " cluster_mass_function.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " cluster_mass_function.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " cluster_mass_function.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of cluster_mass_function.cc program  Rc= " << rc << endl;
  return rc;	
}
	   
