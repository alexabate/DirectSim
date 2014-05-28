/**
  * @file  baseSimulation.cc
  * @brief Program that simulates type, absolute magnitude values from an input
  *        luminosity function       
  *
  * @author Alex Abate
  * @date January 2012 
  * Contact: abate@email.arizona.edu
  *
  */
#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

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
#include "ctimer.h"

// stuff in classes/ dir
#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "sinterp.h"


void usage(void);
void usage(void)
{
	cout << endl<<"  Usage: baseSimulation [...options...]" << endl << endl;

	cout << "  Program that draws redshift, type, absolute magnitude\n";
	cout << "   (i.e. luminosity) values from input luminosity\n";
	cout << "  functions" <<endl<<endl;
	
	cout << "  If options -z and/or -m are not selected the cumulative distributions\n";
	cout << "  are calculated instead of being read from files: this can take"<<endl;
	cout << "  a while"<<endl<<endl;
	
	cout << "  If reading in one or both cumulative distributions, make sure that\n";
	cout << "  the zmin and zmax of the cumulative distributions match the survey\n";
	cout << "  that you wish to simulate "<<endl<<endl;

	cout << "  Output: (to output/)"<<endl;
	cout << "	- a FITS file containing redshift, absolute magnitude and broad";
	cout << " galaxy type" <<endl;
	cout << "  Optionally:"<<endl;
	cout << "   - a FITS binary table containing the cumulative redshift distribution"<<endl;
	cout << "   - a FITS file containting the 2D cumulative magnitude distribution at\n";
	cout << "     a series of redshifts"<<endl<<endl;
                                
 	cout << "  AA April 2010"<<endl<<endl;

	cout << " -o : outfileroot: root filename of outputs, stored in output/ dir (DEFAULT=baseSimulation)"<<endl;
	cout << " -z : zfile: read cumulative redshift distribution from FITS binary table file";
	cout << "  [ZFILE]"<<endl;
	cout << " -m : mfile: read cumulative mag distribution array from FITS file [MFILE]"<<endl;
	cout << " -a : skyarea: area of sky in square degrees (DEFAULT=2)"<<endl;
	cout << " -Z : zMin,zMax: minimum,maximum redshift to simulate from, up to (DEFAULT=0.01,3)"<<endl;
	cout << " -M : mMin,mMAx,zmMin,zmMax: minimum,maximum magnitude, redshift used to calculate"<<endl;
	cout << "      the magnitude, redshift CDF with (DEFAULT=-24,-13,0,6) "<< endl;
	cout << "      Do not need this if cumulative mag distribution is read from file with option -m"<<endl;
	cout << endl;
}

int main(int narg, char* arg[])
{
    cout << " ==== baseSimulation.cc program , to make simulation of galaxy data";
    cout << " from LF ====" << endl;

    // make sure SOPHYA modules are initialized 
    SophyaInit();  
    FitsIOServerInit();
    InitTim();
    cout<<endl<<endl;

    //--- decoding command line arguments 
    string outfileroot="baseSimulation";
    string zfile; bool ReadZ=false;
    string mfile; bool ReadM=false;
    double skyArea = 2; // in square degrees
    double zMin=0.01, zMax = 3.;
    double mMin=-24., mMax=-13.;
    double zMinMag = 0., zMaxMag = 6.;
  
    cout << " decoding command line arguments ..."<<endl;
    char c;
    while((c = getopt(narg,arg,"ho:z:m:a:Z:M:")) != -1)  {
        switch (c)  {
		    case 'o' :
			    outfileroot = optarg;
			    break;
		    case 'z' :
			    zfile = optarg;
			    ReadZ = true;
			    break;
		    case 'm' :
			    mfile = optarg;
			    ReadM = true;
			    break;
		    case 'a' :
		        sscanf(optarg,"%lf",&skyArea);
		        break;
		    case 'Z' :
		        sscanf(optarg,"%lf,%lf",&zMin,&zMax);
		        break;
		    case 'M' :
		        sscanf(optarg,"%lf,%lf,%lf,%lf",&mMin,&mMax,&zMinMag,&zMaxMag);
		        break;
		    case 'h' :
		    default :
			    usage(); return -1;
		    }
	}
    //-- end command line arguments

  cout << "     Output will be saved to files beginning "<<outfileroot<<endl;
  if (ReadZ) {
  	cout << "     Cumulative redshift distribution will be read from ";
  	cout << zfile;
  	}
  else
	cout << "     Explicitly calculating cumulative redshift distribution"; 
  cout<<endl;
  if (ReadM)
  	{
  	cout << "     Cumulative magnitude distribution will be read from ";
  	cout <<mfile;
  	}
  else
	cout << "     Explicitly calculating cumulative magnitude distribution"; 
  cout<<endl;
    
  cout << "     Simulating galaxies over redshifts "<< zMin <<" to "<< zMax;
  cout << " with absolute magnitudes "<< mMin <<" to "<< mMax << endl; 
  cout << "     And over "<< skyArea <<" square degrees "<<endl;
  cout << " ... finished decoding command line arguments "<<endl;
  cout <<endl;

  int rc = 1;  
  try {  // exception handling try bloc at top level
	
    // to output stuff to a file	
	string outfile;
	ifstream inp;
	ofstream outp;

	// this controls the drawing of random numbers
	RandomGenerator rg;

	// cosmological parameters
	cout <<"     Set cosmology ..."<<endl;
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout << "     h = "<<h<<", OmegaM = "<<OmegaM<<", OmegaL = "<<OmegaL;
	cout <<endl<<endl;

	// LF parameters: define evolution of LF as a function of redshift
	cout << "     Initializing LF parameters as a function of redshift"<<endl;
	double zmax = 6., dz = 0.01;
	if (zMax>zmax)
	    throw ParmError("ERROR! LF interpolation does not go up to max redshift");
	LFParameters lfpars(zmax, dz);
	cout << endl;
	
	// Redshift distributions: pre-calculate or read in an interpolation table
	// of F_z(z)
	Timer tm("TimingCumulDistZ",false);
	cout <<"     Setting up cumulative redshift distribution"<<endl;
	CumulDistZ cumlz;
    if (ReadZ) {
		cout <<"     ... reading cumulative redshift distribution from";
		cout <<" file"<<endl;
		cumlz.SetUp(zfile);
		cumlz.PrintDist(10);
		}
	else {
	    // this takes around an hour or two
		tm.Split();
		cout <<"     ... explicitly calculating cumulative redshift";
		cout <<" distribution"<<endl;
		cumlz.SetUp(lfpars, su, zMin, zMax); tm.Split();
		cout <<"     .... done, took "<<tm.PartialElapsedTime()/60<<" mins";
		cout <<endl;
		outfile = "output/" + outfileroot;
		cumlz.Output2File(outfile);
		cumlz.PrintDist(10);
		}
	cout <<endl;
	
	// Prepare redshift distribution to draw from
	DrawZ drz(cumlz,rg); 
	cout << endl;
	
	
	// Magnitude distribution to draw from
	cout <<"     Setting up cumulative magnitude distribution"<<endl;
	tm.Split();
	CumulDistM cumlm(lfpars,su); // calculates magnitude cdf as a function of z
	tm.Split();
	//cout<<"     Time to initialize CumulDistM "<<tm.PartialElapsedTimems();
	//cout <<" ms"<<endl;

    // Magnitude distributions: pre-calculate or read in an interpolation table
	// of F_m(m,z)
	DrawM drm(cumlm);
	if (ReadM) {
		cout <<"     ... reading cumulative magnitude distribution";
		cout <<" from file"<<endl;
		drm.SetUp(mfile);
		}
	else {
	    // this takes a few hours
		tm.Split();
		cout <<"     ... explicitly calculating magnitude distribution";
		cout <<endl<<endl;
		drm.SetUp(rg, mMin, mMax, zMinMag, zMaxMag); tm.Split();
		cout <<"     .... done, took "<<tm.PartialElapsedTime()/60<<" mins";
		cout <<endl<<endl;
		outfile = "output/"+outfileroot;
		drm.Output2File(outfile);
		}
	cout <<endl;
	
    // Galaxy type distributions
	cout <<"     Initialize type distribution"<<endl;
	// get low redshift galaxy type LF parameters
	vector<double> all=lfpars.ReturnAll();
	Schechter schA(all[0],all[1],all[2]);
	vector<double> early=lfpars.ReturnEarly();
	Schechter schE(early[0],early[1],early[2]);
	vector<double> late=lfpars.ReturnLate();
	Schechter schL(late[0],late[1],late[2]);
	vector<double> sb=lfpars.ReturnSB();
	Schechter schS(sb[0],sb[1],sb[2]);
	// calculate type fraction at low redshift
	TypeRatio0 tr0(schA,schE,schL,schS);
	// calculate type fraction redshift evolution
	TypeRatio tr(tr0);
	// set up galaxy type drawing
	DrawType drt(tr,rg);
	cout <<"     .... done"<<endl;
	cout << endl;
	
	cout << "     Setting up galaxy survey simulation"<<endl;
	// number of galaxies in the survey
	NGalZ ng(lfpars,su);
	double sa = NSTERADIANS_IN_SQDEG*skyArea; // area of sky in steradians
	long ntot=ng(zMin, zMax, sa);
	cout <<"     "<< ntot <<" galaxies over "<< skyArea <<" degrees or "<< sa;
	cout <<" steradians between "<< zMin <<"<z<"<< zMax <<endl;
	cout<<endl;

	// full simulation
	cout <<"     Run galaxy simulation "<<endl;
	stringstream ss;
    ss << skyArea;
	outfile="output/" + outfileroot + "_" + ss.str() + "SQDEG";
	SimBaseCatalog sbc(drz,drm,drt);
	sbc.DoSim(ntot, outfile); // this part is fairly fast
	
}

  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " baseSimulation.cc: Catched Exception (PThrowable)" 
    	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " baseSimulation.cc: Catched std::exception "  << " - what()= " 
         << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " baseSimulation.cc: some other exception (...) was caught ! "<< endl;
    rc = 97;
  }
  cout << " ==== End of baseSimulation.cc program  Rc= " << rc << endl;
  return rc;	
}
