// -*- LSST-C++ -*-
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

#define PI 3.141592

/*******************************************************************************
  Program that simulates type, absolute magnitude values 
  from an input luminosity function     
             
  A Abate January 2012             
                
*******************************************************************************/

void usage(void);
void usage(void)
{
	cout << endl<<"  Usage: testbasesim [...options...]" << endl<<endl;

	cout << "  Testing program for drawing redshift, type, absolute\n";
	cout << "  magnitude (i.e. luminosity) values from input luminosity\n";
	cout << "  functions" <<endl<<endl;

	cout << "  Output: (to testing/)"<<endl;
	cout << "	- stuff "<<endl<<endl;
                                
 	cout << "  AA April 2010"<<endl<<endl;

	cout << " -o : outfileroot: root filename of outputs (DEFAULT:";
	cout << " testbasesim)"<<endl;
	cout << " -z : zfile: read cumulative redshift distribution from file";
	cout << "  [ZFILE]"<<endl;
	cout << " -m : mfile: read cumulative mag distribution from file";
	cout << "  [MFILE]"<<endl;
	cout << " -a : skyarea: area of sky in square degrees (DEFAULT=2)"<<endl;
	cout << endl;
}

int main(int narg, char* arg[])
{
  cout << " ==== testbasesim.cc program , to test simulation of galaxy data";
  cout << " from LF ====" << endl;

  // make sure SOPHYA modules are initialized 
  SophyaInit();  
  FitsIOServerInit();
  InitTim();
  cout<<endl<<endl;

  //--- decoding command line arguments 
  string outfileroot="testbasesim";
  string zfile; bool ReadZ=false;
  string mfile; bool ReadM=false;
  double skyArea = 2; // in square degrees
  cout << " decoding command line arguments ..."<<endl;
  char c;
  	while((c = getopt(narg,arg,"ho:z:m:a:")) != -1) 
	{
	switch (c) 
		{
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
		case 'h' :
		default :
			usage(); return -1;
		}
	}
    //-- end command line arguments

  cout << "     Something will be saved to "<<outfileroot<<endl;
  if (ReadZ)
  	{
  	cout << "     Cumulative redshift distribution will be read from ";
  	cout <<zfile;
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

	// redshift range
	double zmax=6, dz=0.05;
	int nz=zmax/dz;

	// magnitude range
	double mmin=-24,mmax=-13;
	int nmag=100;
	double dm=(mmax-mmin)/(nmag-1);

	// LF parameters
	cout << "     LF parameters as a function of redshift"<<endl;
	LFParameters lfpars(zmax,dz); cout <<"made it through lfpars"<<endl;
	string fname="GOODS_B_LF.txt";
	LFParameters lfpars_table(fname,1);
	double ms_table,a_table,ps_table;
	lfpars_table.ReturnParsBini(ms_table,a_table,ps_table,0,0);
	cout <<endl;
	

	// LF parameters as a function of z
	cout << "     Output LF parameters as a function of redshift"<<endl;
	outfile="testfiles/lfpars.txt";
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
	  	{
		inp.clear(ios::failbit);
		cout << "     Writing to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		for (int i=0; i<nz; i++)
			{
			double z=i*dz;
			double ps,ms,a;
			lfpars(z,ps,ms,a);
			outp << z <<"  "<<ps<<"  "<<ms<<"  "<<a<<endl;
			}
		outp.close();
		}
	else
		cout <<"ERROR file "<<outfile<<" exists"<<endl;
	cout <<endl;


	cout <<"     Set z=0.3, output LF from i) extrapolated LF pars,";
	cout <<"  ii) measured LF pars"<<endl;
	double z=0.3;
	double ms,ps,a;
	lfpars(z,ps,ms,a);
	cout << "     Extrapolated LF parameters at z=0.3 are: ps = "<<ps;
	cout << ", ms = "<<ms<<", a = "<<a<<endl;
	cout << "     Measured LF parameters at z=0.3 are: ps = "<<ps_table;
	cout << ", ms = "<<ms_table<<", a = "<<a_table<<endl;
	Schechter sch(ps,ms,a);
	Schechter schtable(ps_table,ms_table,a_table);
	outfile="testfiles/lf.txt";
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
	  	{
		inp.clear(ios::failbit);
		cout << "     Writing to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		for (int i=0; i<nmag; i++)
			{
			double m=mmin+i*dm;
			double lf=sch(m);
			double lft=schtable(m);
			outp << m <<"  "<<lf<<"  "<<lft<<endl;
			}
		outp.close();
		}
	else
		cout <<"ERROR file "<<outfile<<" exists"<<endl;
	cout << endl;

	cout << "     Write phi(z) = int phi(M|z) dM dV(z) to a file"<<endl;
	// phi(z) = int phi(M|z) dM dV(z)
	SchechterZVol schZ(lfpars,su);
	outfile="testfiles/schzvol.txt";
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
	  	{
		inp.clear(ios::failbit);
		cout << "     Writing to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		for (int i=0; i<nz; i++)
			{
			double z=i*dz;
			double szv=schZ(z);
			outp << z <<"  "<<szv<<endl;
			}
		outp.close();
		}
	else
		cout <<"ERROR file "<<outfile<<" exists"<<endl;
	cout <<endl;
	
	// redshift dists
	Timer tm("TimingCumulDistZ",false);
	cout <<"     Setting up cumulative redshift distribution"<<endl;
	CumulDistZ cumlz;
        if (ReadZ)
		{
		cout <<"     ... reading cumulative redshift distribution from";
		cout <<" file"<<endl;
		cumlz.SetUp(zfile);
		cumlz.PrintDist(10);
		}
	else
		{
		tm.Split();
		cout <<"     ... explicitly calculating cumulative redshift";
		cout <<" distribution"<<endl;
		cumlz.SetUp(lfpars,su); tm.Split();
		cout <<" .... done, took "<<tm.PartialElapsedTime()/60<<" mins";
		cout <<endl;
		outfile="testfiles/"+outfileroot;
		cumlz.Output2File(outfile);
		cumlz.PrintDist(10);
		}
	cout <<endl;

	cout << "     Setting up galaxy survey simulation"<<endl;
	// number of galaxies in the survey
	NGalZ ng(lfpars,su);
	double sa = NSTERADIANS_IN_SQDEG*skyArea; // area of sky in steradians
	long ntot=ng(0.,zmax,sa);
	cout <<"     "<<ntot<<" galaxies over "<<sa<<" steradians between 0<z<";
	cout <<zmax<<endl;
	cout<<endl;
		
	// galaxy fraction stuff
	vector<double> all=lfpars.ReturnAll();
	Schechter schA(all[0],all[1],all[2]);
	vector<double> early=lfpars.ReturnEarly();
	Schechter schE(early[0],early[1],early[2]);
	vector<double> late=lfpars.ReturnLate();
	Schechter schL(late[0],late[1],late[2]);
	vector<double> sb=lfpars.ReturnSB();
	Schechter schS(sb[0],sb[1],sb[2]);
	TypeRatio0 tr0(schA,schE,schL,schS);
	TypeRatio tr(tr0);
		
	// Redshift distribution to draw from
	DrawZ drz(cumlz,rg); cout << endl;
	cout << "     Check drawing redshifts"<<endl;
	for (int i=0; i<100;i++)
		{
		double z=drz.Draw();
		cout << "z = "<<z<<endl;
		}
	
	
	
	// Magnitude distribution to draw from
	tm.Split();
	CumulDistM cumlm(lfpars,su);
	tm.Split();
	cout<<"     Time to initialize CumulDistM "<<tm.PartialElapsedTimems();
	cout <<" ms"<<endl;

	DrawM drm(cumlm);
	if (ReadM)
		{
		cout <<"     ... reading cumulative magnitude distribution";
		cout <<" from file"<<endl;
		drm.SetUp(mfile);
		}
	else
		{
		tm.Split();
		cout <<"     ... explicitly calculating magnitude distribution";
		cout <<endl;
		drm.SetUp(rg); tm.Split();
		cout <<" .... done, took "<<tm.PartialElapsedTime()/60<<" mins";
		cout <<endl;
		outfile="testfiles/"+outfileroot;
		drm.Output2File(outfile);
		}
	cout <<endl;
		
	// types
	cout <<"     Initialize type distribution"<<endl;
	DrawType drt(tr,rg);
	cout << endl;

	// full simu
	cout <<"     Run galaxy simulation "<<endl;
	outfile="testfiles/"+outfileroot+".fits";
	SimBaseCatalog sbc(drz,drm,drt);
	sbc.DoSim(ntot,outfile);
	
}

  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testbasesim.cc: Catched Exception (PThrowable)" 
    	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testbasesim.cc: Catched std::exception "  << " - what()= " 
         << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " testbasesim.cc: some other exception (...) was caught ! "<< endl;
    rc = 97;
  }
  cout << " ==== End of testbasesim.cc program  Rc= " << rc << endl;
  return rc;	
}
