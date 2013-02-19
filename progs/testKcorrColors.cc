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

#include "schechter.h"
#include "cosmocalcs.h"
#include "gftdist.h"
#include "simdata.h"
#include "sedfilter.h"

#define PI 3.141592
/******************************************************
*                                                     *
* PROGRAM TO TEST K CORRECTION AND COLOR COMPUTATION  *
* FUNCTIONS IN SIMDATA CLASS                          *
*                                                     *
* CALCULATES AS A FUNCTION OF Z:                      *
* a) mg-mr and mr-mi for SDSS-like/LSST filters       *
* b) k correction IN B,V,R FILTERS                    *
*                                                     *
******************************************************/
int main(int narg, char* arg[])
{
  cout << " ==== testKcorrColor.cc program , to check k correction & color";
  cout << " computation  ====" << endl << endl;
  if (narg < 2) {
        cout << " Usage: testKcorrColor outfile" << endl;
	cout << "        outfile: output filename"<<endl;
	cout << endl;
	cout << " Calculates:"<<endl;
	cout << "          i) the g-r and r-i color (LSST filters)"<<endl;
	cout << "         ii) B-R and R-I color (GOODS filters)"<<endl;
	cout << "        iii) k-correction for an elliptical galaxy in B,V,R";
	cout << " (GOODS filters)"<<endl;
	cout << "         iv) k-correction for an Sbc galaxy in B,V,R (GOODS";
	cout << " filters)"<<endl;
	cout << "          v) k-correction for an Scd galaxy in B,V,R (GOODS";
	cout << " filters)"<<endl;
	cout << "         vi) k-correction for an Im galaxy in B,V,R (GOODS";
	cout << " filters)"<<endl;
	cout << "        vii) k-correction for an elliptical galaxy in g (LSST";
	cout << " filter)"<<endl<<endl;
	cout << " Outputs a file with the following columns:"<<endl;
	cout << " redshift C_gr C_ri C_BR CRI kBe kVe kRe kBs1 kVs1 kRs1 kBs2";
	cout << " kVs2 kRs2 kBi kVi kRi kge"<<endl<<endl;
	cout << " * Compare kBe(z),kVe(z), kRe(z), kBs1(z), kVs1(z), kRs1(z),";
	cout << " kBs2(z), kVs2(z), kRs2(z),"<<endl;
	cout << " kBi(z), kVi(z), kRi(z) to figure 9 in deLapparent et al";
	cout << " (2004)"<<endl<<endl;
	cout << " * Compare C_gr vs C_ri to figure 2 in Padmanabhan et al";
	cout << " (2005) (diamonds in the top panel)"<<endl;
	cout <<endl;
    return 1;
  }

  // make sure SOPHYA modules are initialized 
  SophyaInit();  
  FitsIOServerInit();
  InitTim();
  cout<<endl<<endl;

  //--- decoding command line arguments 
  string outfile="testfiles/";
  string fname=arg[1];
  outfile+=fname;
  
  //-- end command line arguments
  
  int rc = 1;  
  try {  // exception handling try bloc at top level
	ResourceUsage res;
	InitTim();

	// wavelength range of the SEDs/filters
	double lmin=5e-8, lmax=2.5e-6;
	
	//FILTERS
	cout <<"   testKcorrColor: Loading in filters ..."<<endl;
	// these filter files must be in the directory pointed to by the 
	// environment variable FILTLOC (type echo $FILTLOC at terminal to 
	// check this)

	// GOODS (Great Observatories Origins Deep Survey) uses a number of 
	// telescopes including Spitzer, HST and Chandra.  One of the telescopes
	// which did the imaging was the 2.2m La Silla telescope.The UBVRI 
	// filters are on the "Johnson & Morgan" system
	cout <<"   ... GOODS"<<endl;
	string goodsFilterFile = "GOODS.filters";
	ReadFilterList readGOODSfilters(goodsFilterFile);
	readGOODSfilters.readFilters(lmin,lmax);
	vector<Filter*> GOODSfilters=readGOODSfilters.getFilterArray();
	int nGOODS=readGOODSfilters.getNTot();
	cout <<"     "<<nGOODS<<" GOODS filters read in "<<endl;
		cout <<"   Check GOODs filters ..."<<endl;
	cout <<" ... printing B filter to the screen"<<endl;
	int nlam=100;
	//double lmin2=3.935e-07, lmax2=5.135e-07;
	double dlam = (lmax-lmin)/(nlam-1);
	for (int i=0;i<nlam;i++) {
	
		double lam = lmin+i*dlam;
		cout << lam <<"  "<< GOODSfilters[1]->operator()(lam) << endl;
		
		}
	cout << endl;
	
    // GOODS filter numbers
	int Bgoods = 1;
	int Vgoods = 2;
	int Rgoods = 3;
	int Igoods = 3;
	
	// LSST filters ugrizy are on the "Thuan-Gunn" system
	// (optimized for faint galaxies)
	cout <<"   ... LSST"<<endl;
	string lsstFilterFile = "LSST.filters";
	ReadFilterList readLSSTfilters(lsstFilterFile);
	readLSSTfilters.readFilters(lmin,lmax);
	vector<Filter*> LSSTfilters=readLSSTfilters.getFilterArray();
	int nLSST=readLSSTfilters.getNTot();
	cout <<"     "<<nLSST<<" LSST filters read in "<<endl;
		cout <<"   Check LSST filters ..."<<endl;
	cout <<" ... printing LSST g filter and GOODS B filters to the screen";
	cout <<endl;
	//lmin2=3.83e-7, lmax2=5.79e-7;
	dlam = (lmax-lmin)/(nlam-1);
	for (int i=0;i<nlam;i++) {
	
		double lam = lmin+i*dlam;
		cout << "lam="<<lam<<", g trans = "<< LSSTfilters[1]->operator()(lam);
		cout <<", B trans = "<< GOODSfilters[1]->operator()(lam) << endl;

		}
	cout << endl;
	
	// LSST filter numbers
	int gLSST = 1;
	int rLSST = 2;
	int iLSST = 3;
		
	// GALAXY SED TEMPLATES
	// these SED files must be in the directory pointed to by the 
	// environment variable SEDLOC (type echo $SEDLOC at terminal to check 
	// this)
	
	// "CWW" stands for Coleman, Wu & Weedman 1980
	// in that paper they measured spectra of typical elliptical (early), 
	// spiral (late, Sbc, Scd) and irregular type galaxies
	// "K" stands for Kinney et al 1996 who measure starburst galaxies
	string sedFile = "CWWK.list";
	ReadSedList readSedList(sedFile);
	// Read out SEDs into array
    readSedList.readSeds(lmin,lmax);
    // Get total number of SEDs
    vector<SED*> sedArray=readSedList.getSedArray();
    int nsed=readSedList.getNSed();
    cout <<"     Number of original SEDs = "<<nsed<<endl;
	cout << endl;
	
	// sed numbers
	int iEl = 0;    // elliptical galaxy
	int iSbc = 1;   // spiral galaxy
	int iScd = 2;   // spiral galaxy
	int iIrr = 3;    // irregular galaxy
	
    // PHOTOMETRY CALCULATION CLASS
    PhotometryCalcs photometryCalcs(lmin,lmax);
    

	// LOOP OVER REDSHIFTS
	double zmin=0, zmax=2.3;
	int nz=100;
	double dz=(zmax-zmin)/(nz-1);
	
	TVector<r_8> kBe(nz),kVe(nz),kRe(nz);
	TVector<r_8> kBs1(nz),kVs1(nz),kRs1(nz);
	TVector<r_8> kBs2(nz),kVs2(nz),kRs2(nz);
	TVector<r_8> kBi(nz),kVi(nz),kRi(nz);
	TVector<r_8> kge(nz);
	TVector<r_8> Cgr(nz), Cri(nz);
	TVector<r_8> Cgrb(nz), Crib(nz);
	TVector<r_8> Cgrc(nz), Cric(nz);
	TVector<r_8> CBR(nz), CRI(nz);
	for (int i=0; i<nz;i++)
		{
		double z=zmin+i*dz;

		// "Color" in astronomy means ratio of flux in filter X compared 
		// to filter Y
		// This is equivalent to the difference in magnitudes in filter
		// X and Y 
		// because: mag ~ log10(flux)
		
		// mag in filter g - mag in filter r 
		Cgr(i)=photometryCalcs.CompColor(z,(*sedArray[iEl]),(*LSSTfilters[gLSST]),(*LSSTfilters[rLSST]));
		// mag in filter r - mag in filter i 
		Cri(i)=photometryCalcs.CompColor(z,(*sedArray[iEl]),(*LSSTfilters[rLSST]),(*LSSTfilters[iLSST]));
		// mag in filter B - mag in filter R 
		CBR(i)=photometryCalcs.CompColor(z,(*sedArray[iEl]),(*GOODSfilters[Bgoods]),(*GOODSfilters[Rgoods]));
		// mag in filter R - mag in filter I
		CRI(i)=photometryCalcs.CompColor(z,(*sedArray[iEl]),(*GOODSfilters[Rgoods]),(*GOODSfilters[Igoods]));

		// The k-correction corrects for the fact that the flux observed
		// through filter X is not the "same" flux that was emitted in 
		// the rest frame of the galaxy within the wavelength region of 
		// filter X

		// k-correction for the elliptical galaxy in the B,V,R bands
		kBe(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iEl]),(*GOODSfilters[Bgoods])); 
		kVe(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iEl]),(*GOODSfilters[Vgoods])); 
		kRe(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iEl]),(*GOODSfilters[Rgoods])); 

		// for the Sbc galaxy in BVR
		kBs1(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iSbc]),(*GOODSfilters[Bgoods])); 
		kVs1(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iSbc]),(*GOODSfilters[Vgoods])); 
		kRs1(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iSbc]),(*GOODSfilters[Rgoods]));
		
		// for the Scd galaxy in BVR
		kBs2(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iScd]),(*GOODSfilters[Bgoods])); 
		kVs2(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iScd]),(*GOODSfilters[Vgoods])); 
		kRs2(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iScd]),(*GOODSfilters[Rgoods]));

		// for the irregular galaxy in BVR
		kBi(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iIrr]),(*GOODSfilters[Bgoods])); 
		kVi(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iIrr]),(*GOODSfilters[Vgoods])); 
		kRi(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iIrr]),(*GOODSfilters[Rgoods]));

		// for the elliptical galaxy in g band
		kge(i)=photometryCalcs.Kcorr1Filter(z,(*sedArray[iEl]),(*LSSTfilters[gLSST]));

		cout <<"Bfilter kEl(z="<<z<<")="<<kBe(i)<<endl;
		}
		
	// WRITE TO A FILE

	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
	  {
		inp.clear(ios::failbit);
		cout << "Writing to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
	    for (int i=0;i<nz;i++)
		  {		
			double z=zmin+i*dz;
			outp << z <<"    "<<Cgr(i)<<"    "<<Cri(i)<<"    ";
			outp <<CBR(i)<<"    "<<CRI(i)<<"    "<<kBe(i)<<"    ";
			outp <<kVe(i)<<"    "<<kRe(i)<<"    "<<kBs1(i)<<"    ";
			outp <<kVs1(i)<<"    "<<kRs1(i)<<"    "<<kBs2(i);
			outp <<"    "<<kVs2(i)<<"    "<<kRs2(i)<<"    "<<kBi(i);
			outp <<"    "<<kVi(i)<<"    "<<kRi(i)<<"    "<<kge(i);
			outp <<"    "<<Cgrb(i)<<"    "<<Crib(i)<<"    ";
			outp <<Cgrc(i)<<"    "<<Cric(i)<<endl;
		  }
		  outp.close();
	  }
	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;
	  

	/*// Loop over redshifts again to check Kxy = Cxc + Kcy
	SimData data2(lmin,lmax);
	string outfile2="testfiles/check.txt";//+outfile;
	inp.open(outfile2.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
	 	{
		inp.clear(ios::failbit);
		cout << "Writing to file ..." << outfile2.c_str() << endl;
		outp.open(outfile2.c_str(), ofstream::out);
		for (int i=0; i<nz;i++)
			{
			double z=zmin+i*dz;

			double kgr=data2.Kcorr(z,sed0,lsstg,lsstr);
			double cgi=data2.CompColor(z,sed0,lsstg,lssti);
			double kir=data2.Kcorr(z,sed0,lssti,lsstr);

			// calculate k-correction parts separately
			// 1) kgr
			SEDzFilterProd szfx(sed0,lsstg,z);
			FilterProd  fprody(lsstr);          
			SEDzFilterProd szfB0(sed0,lsstr,0.);
			FilterProd  fprodx(lsstg);             
			FilterIntegrator trpzx(szfx,lmin,lmax);	
			FilterIntegrator trpzfy(fprody,lmin,lmax);	
			FilterIntegrator trpzy(szfB0,lmin,lmax);
			FilterIntegrator trpzfx(fprodx,lmin,lmax);	
			double kgr2=-2.5*log10(pow((1+z),-1)*(trpzx.Value()/trpzfx.Value())
					*(trpzfy.Value()/trpzy.Value())); 

			outp<< z <<"  "<<kgr<<"  "<<cgi<<"  "<<kir<<"  "<<kgr2<<endl;
			}
		outp.close();
		}
	else
		cout << "Error...file """ << outfile2.c_str() << """ exists" << endl;*/
  }  // End of try bloc 
  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " testKcorrColor.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " testKcorrColor.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " testKcorrColor.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of testKcorrColor.cc program  Rc= " << rc << endl;
  return rc;	
}
