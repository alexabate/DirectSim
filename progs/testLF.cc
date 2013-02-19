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

// stuff in classes/ dir
#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"

#define PI 3.141592

/**************************************************
  Testing program for drawing type,abs mag values 
  from input luminosity function     
             
  A Abate January 2012             
                
**************************************************/

void usage(void);
void usage(void)
{
	cout << endl<<"  Usage: testLF [...options...]" << endl<<endl;

	cout << "  Testing program for drawing redshift, type, absolute mag\n";
	cout << "  (i.e.luminosity) values from input luminosity functions\n\n";

	cout << "  This is the old method: uses the GOODS LFs without tweaking";
	cout << "\n  to get a better match to the data"<<endl<<endl;

	cout << "  Simulates broad galaxy type and absolute magnitude for a \n";
	cout << "  galaxy survey 0<z<3 over 0.3sq.deg"<<endl;
	cout << "  Output: (to testing/) simulated gals:"<<endl;
	cout << "	- 0<z<0.5 (*bin1.txt) "<<endl;
	cout << "	- 0.5<z<0.75 (*bin2.txt) "<<endl;
	cout << "	- z>0.75 (*bin3.txt) "<<endl;
	cout <<endl;
                                
 	cout << "  AA April 2010"<<endl<<endl;

	cout << " -o : outfileroot: root filename of outputs (default testlf)";
	cout << endl;
	cout << endl;
}

int main(int narg, char* arg[])
{
  cout << " ==== testLF.cc program , to test simulation of galaxy data from LF"; 
  cout << " ==== " << endl;

  // make sure SOPHYA modules are initialized 
  SophyaInit();  
  FitsIOServerInit();
  InitTim();
  cout<<endl<<endl;

  //--- decoding command line arguments 
  string outfileroot="testlf";
  cout << " decoding command line arguments ..."<<endl;
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

  cout << "     Something will be saved to testing/"<<outfileroot<<endl;
  cout << " ... finished decoding command line arguments "<<endl;
  cout <<endl;

  int rc = 1;  
  try {  // exception handling try bloc at top level
	

	// this controls the drawing of random numbers
	RandomGenerator rg;

	cout <<"     Set cosmology ..."<<endl;
	double h=0.7, OmegaM=0.3, OmegaL=0.7;
	SimpleUniverse su(h,OmegaM,OmegaL);
	su.SetFlatUniverse_OmegaMatter();
	cout << "     h = "<<h<<", OmegaM = "<<OmegaM<<", OmegaL = "<<OmegaL<<endl;
	cout <<endl;

	// Read in LF parameters from file "fname"
	// this file must be in the directory pointed to by the 
	// environment variable LFLOC (type echo $LFLOC at terminal to 
	// check this)
	string fname="GOODS_B_LF.txt";
	LFParameters lfpars(fname,1);
	lfpars.PrintZbins();
	//lfpars.PrintPars(0);
	//lfpars.PrintPars(1);
	//lfpars.PrintPars(2);
	//lfpars.PrintPars(3);

	// Galaxy survey characteristics
	// between zmin and zmax over a sky area of sqdeg
	double zmin=0, zmax=3;
	double sqdeg=0.306;
	double sqrad=sqdeg*(PI/180)*(PI/180);
	cout <<"     Galaxy survey from redshift "<<zmin<<" to redshift "<< zmax<<endl;
	cout <<"     Over an area of "<<sqdeg<<" deg^2 or "<<sqrad<<" steradians"<<endl;
	cout <<endl;

	// Work out how many LF bins the galaxy survey (defined by zmin and zmax 
	// above) spans across
	int imin,imax,nb;
	lfpars.CheckBinRange(imin,imax,nb,zmin,zmax);

	// Calculate number of galaxies in each of these redshift bins
	vector<double> ngbin;
	for (int i=0; i<nb; i++)
		{
		int ib = imin+i;
		double Mstar,alpha,phistar;
		// return the LF parameters for All galaxies LF in bin "ib"
		lfpars.ReturnParsBini(Mstar,alpha,phistar,ib,0);
		cout << "     Parameters of bin "<<i+1<<" are:"<<endl;
		cout << "     Mstar="<<Mstar<<", alpha="<<alpha<<", phistar=";
		cout <<phistar<<endl;
		// set up Schechter function
		Schechter sch(phistar,Mstar,alpha);
		double z1,z2;
		lfpars.ReturnZbin(z1,z2,ib);
		if (ib==0&&zmin<z1)//if 1st bin and min survey z is < lower bin 
			z1=zmin;   //edge
		if (ib>(nb-2)&&zmax>z2)//if last bin and max survey z is > upper
			z2=zmax;       // bin edge
		// integrate Schechter function to return number of gals in this 
		// bin
		double ng=sch.NGals(su,sqrad,z1,z2);
		cout << "     Number of galaxies in this redshift bin = "<<ng;
		cout << endl;
		ngbin.push_back(ng);
		cout <<endl;
		}
	cout << endl;

	// for integration below
	double schmin=-24, schmax=-13;// units of "M-5log10h70"
	int schnpt=10000;
	int magbin=100000;// for MB-type dist

	// Loop over redshift bins and draw galaxies
	for (int ib=0; ib<nb; ib++)
		{
		cout <<"     Redshift bin "<<ib+1<<" of "<<nb<<endl;
		// 1) set up each separate Schechter function
		double Mstar,alpha,phistar;
		lfpars.ReturnParsBini(Mstar,alpha,phistar,ib,0);
		//cout << "     Mstar="<<Mstar<<", alpha="<<alpha<<", phistar=";
		//cout <<phistar<<endl;
		Schechter schA(phistar,Mstar,alpha);
		lfpars.ReturnParsBini(Mstar,alpha,phistar,ib,1);
		//cout << "     Mstar="<<Mstar<<", alpha="<<alpha<<", phistar=";
		//cout <<phistar<<endl;
		Schechter schE(phistar,Mstar,alpha);
		lfpars.ReturnParsBini(Mstar,alpha,phistar,ib,2);
		//cout << "     Mstar="<<Mstar<<", alpha="<<alpha<<", phistar=";
		//cout <<phistar<<endl;
		Schechter schL(phistar,Mstar,alpha);
		lfpars.ReturnParsBini(Mstar,alpha,phistar,ib,3);
		//cout << "     Mstar="<<Mstar<<", alpha="<<alpha<<", phistar=";
		//cout <<phistar<<endl;
		Schechter schSB(phistar,Mstar,alpha);
		// 2) Re-normalize the type=specific distributions 
		// so that they give the correct total amount of 
		// galaxies
		int type;
		type=1;
		SchechterDist schDistE(schA,schE,schL,schSB,type);
		type=2;
		SchechterDist schDistL(schA,schE,schL,schSB,type);
		type=3;
		SchechterDist schDistS(schA,schE,schL,schSB,type);
		cout <<endl;
		// 3) integrate work out type fractions
		double nfE=schDistE.Integrate(schmin,schmax,schnpt);
		double nfL=schDistL.Integrate(schmin,schmax,schnpt);
		double nfS=schDistS.Integrate(schmin,schmax,schnpt);
		double totalnr=nfE+nfL+nfS;
		double fE=nfE/totalnr;
		double fL=nfL/totalnr;
		double fS=nfS/totalnr;
		// 4) Make 2D distribution of galaxy type vs magnitude
		int PrtLevel = 0;
		GalFlxTypDist gfd(rg, PrtLevel); 	
		gfd.AddGalType(schDistE,schmin,schmax,fE,magbin,schnpt); 
		gfd.AddGalType(schDistL,schmin,schmax,fL,magbin,schnpt); 
		gfd.AddGalType(schDistS,schmin,schmax,fS,magbin,schnpt);
		// 5) draw galaxies and write them to a file
		stringstream ss;
		ss << ib+1;
		string outfile="testfiles/"+outfileroot+"_bin"+ss.str()+".txt";
		ifstream inp;
    		ofstream outp;
		inp.open(outfile.c_str(), ifstream::in);
		inp.close();
		if(inp.fail())
			{
			inp.clear(ios::failbit);
			cout << "     Writing to file ..." << outfile.c_str();
			cout << endl << endl;
			outp.open(outfile.c_str(), ofstream::out);
			double z1,z2;
			lfpars.ReturnZbin(z1,z2,ib);
			if (ib==0&&zmin<z1)//if 1st bin and min survey z < lower
				z1=zmin;   // bin edge
			if (ib>(nb-2)&&zmax>z2)//if last bin and max survey z >
				z2=zmax;       //upper bin edge
			double vol=schA.Volume(su,sqrad,z1,z2);
			int ng=int(ngbin[ib]);
			cout <<"     Writing "<<ng<<" galaxies to the file\n";
			outp <<" # volume = "<<vol<<" Mpc^3"<<endl;
			for(int i=0; i<ng; i++) 
				{
				double mag;
				int typ;
				gfd.GetGalaxy(typ,mag); 
				outp<<typ<<"    "<<mag<<endl;
				}
			outp.close();
			cout << endl;
	  		}
	  	else
			cout <<"Error...file "<<outfile.c_str()<<" exists\n";
		}
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testLF.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testLF.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " testLF.cc: some other exception (...) was caught ! " << endl;
    rc = 97;
  }
  cout << " ==== End of testLF.cc program  Rc= " << rc << endl;
  return rc;	
}
