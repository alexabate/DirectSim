/**
 * @file  getsf.cc
 * @brief Compute selection functions \f$ n_{pz}(z)/n_{true}(z) \f$ and
 *        \f$ n_{sz}(z)/n_{true}(z) \f$
 *
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: Aug 2010
 * @date Aug 2010
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

// sophya
#include "sopnamsp.h"
#include "histinit.h"
#include "hisprof.h"
#include "histerr.h"
#include "histos.h"
#include "datatable.h"
#include "fitshdtable.h"
#include "swfitsdtable.h"
#include "fitsarrhand.h"
#include "fiosinit.h"
#include "tarray.h"

// DirectSim
#include "geneutils.h"
#include "cat2grid.h"
#include "cosmocalcs.h"
#include "constcosmo.h"

void usage(void);
void usage(void) {

	cout << endl<<" Usage: getsf [...options...]" << endl<<endl;
	
    cout << "  Compute selection functions n_{pz}(z)/n_{true}(z) and    "<<endl;
    cout << "  n_{sz}(z)/n_{true}(z)                                    "<<endl;
    cout << endl;
	
	cout << "  Read in a catalog containing all the simulated galaxies  "<<endl;
	cout << "  and a catalog of the observed galaxies. The catalog that "<<endl;
	cout << "  that contains all the simulated galaxies may be stored in"<<endl;
	cout << "  more than one file. Supply all these filenames separated "<<endl;
	cout << "  by commas to the program with the -F option. Supply the  "<<endl;
    cout << "  filename of the catalog of observed galaxies with the -O "<<endl;
    cout << "  option.                                                  "<<endl;
	cout << endl;
	
	cout << "  The -z option passes the names of the columns to read the"<<endl;
	cout << "  redshifts from, separated by commas, in the following    "<<endl;
	cout << "  order: observed photometric redshifts, observed spec     "<<endl;
	cout << "  redshifts, simulated redshifts. The default names are:   "<<endl;
	cout << "  'zp', 'z' and 'z'                                        "<<endl;
	cout << endl;

    cout << "  The selection function n_{sz}/n_{true} will be written to"<<endl;
    cout << "  the file [SFTextFile]_specz_nofz.txt and the selection   "<<endl;
    cout << "  function n_{pz}/n_{true} will be written to the file     "<<endl;
    cout << "  [SFTextFile]_nofz.txt, where SFTextFile is supplied to   "<<endl;
    cout << "  the program using the -o option.                         "<<endl;
	
	cout << "  If the -d option is used ppf files containing the n(z)   "<<endl;
	cout << "  are output to a file                                     "<<endl;
	cout << endl;
	
	cout << "  EXAMPLES: " <<endl;
	
	cout << "  $ getsf -F full.fits -O obs.fits -o selectfunc -z zp,z,z "<<endl;
    cout << "  $ getsf -F full1.fits,full2.fits,full3.fits -O obs.fits  "<<endl;
    cout << "                                  -o selectfunc -z zp,z,z  "<<endl;
    cout << endl;
	
	cout << " -F : FullCat    FITS filename containing simulated catalog(s)  "<<endl;
	cout << " -O : ObsCat     FITS file containing observed catalog          "<<endl;
	cout << " -o : sfunc      root name of selection function file           "<<endl;
	cout << " -z : zp,zs,zz   column names to read redshifts from (see above)"<<endl;
	cout << " -d : [noarg]    Save n(z) Histos to a ppf file                 "<<endl;
	//cout << " -N : nFiles: Number of FITS files containing RDLSS output [default=1]"<<endl;
	cout << endl;
}


int main(int narg, char *arg[])	{

	SophyaInit();
	FitsIOServerInit();
  
	string FullCat, ObsCat, SFFileName;
	//int nFiles = 1;
	bool DoDebug = false;
	bool isZColSpecified = false;
	string ZCol;
	string ZOCol = "zp"; // read in photo-z's as OBSERVED redshifts from OBSERVED catalog
	string ZSCol = "z";  // read in SPECTRO-z from OBSERVED catalog
	string ZFCol = "z";  // read in SPECTRO-z from FULL catalog
	
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hdF:O:o:z:")) != -1) { 
	    switch (c) {
		    case 'F' :
			    FullCat = optarg;
			    break;
//			case 'N' :
//			    sscanf(optarg,"%d",&nFiles);
//			    break;
		    case 'O' :
			    ObsCat = optarg;
			    break;
		    case 'o' :
			    SFFileName = optarg;
			    break;
		    case 'z' :
			    ZCol = optarg;
			    isZColSpecified = true;
			    break;
		    case 'd' :
			    DoDebug=true;
			    break;
		    case 'h' :
		        default :
			    usage(); return -1;
		        }
	        }

    // read in names of the redshift columns
	if (isZColSpecified) { 
		string delim=",";
		vector<string> results;
		stringSplit(ZCol,delim,results);
		vector<string>::iterator i;
		i = results.begin();
		ZOCol=*i; // column name of observed redshifts
		i++;
		if (i!=results.end())
			ZSCol=*i; // column name of spec redshifts
		i++;
		if (i!=results.end())
			ZFCol=*i; // column name of spec redshifts in full catalog
		}
	
	cout << "     Printing command line arguments ... "<<endl<<endl;
	cout << "     Observed catalog in file "<< ObsCat <<endl;
	cout << "     OBSERVED redshifts to be read from "<< ZOCol <<" column"<<endl;
	cout << "     SPECTRO redshifts to be read from "<< ZSCol <<" column"<<endl;
	cout << "     True redshifts in file(s) " << FullCat << endl;
	cout << "     SPECTRO redshifts to be read from "<< ZFCol <<endl;
	cout << "     Selection function will be written to "<< SFFileName <<"_nofz.txt and ";
	cout << SFFileName <<"_specz_nofz.txt"<<endl;	
	if (DoDebug)
		cout << "     Saving n(z)'s to ppf file "<< SFFileName <<"_histo.ppf"<< endl;
	cout <<endl;
	
  try {
  
		// Read in redshifts from observed catalog
		cout <<"0/ Read in observed catalog "<< ObsCat <<endl;
		FitsInOutFile fin(ObsCat, FitsInOutFile::Fits_RO);
		fin.MoveAbsToHDU(2);
		SwFitsDataTable dt(fin,512,false);
		sa_size_t Izs = dt.IndexNom(ZSCol);
		sa_size_t Izp = dt.IndexNom(ZOCol);
		double minzo,maxzo,minzop,maxzop;
		dt.GetMinMax(Izs, minzo, maxzo); 
		dt.GetMinMax(Izp, minzop, maxzop); 
		cout <<"    Min SPECTRO z of OBS catalog = "<< minzo;
		cout <<", max SPECTRO z of OBS catalog = "<< maxzo <<endl;
		cout <<"    Min OBSERVED z of OBS catalog = "<< minzop;
		cout <<", max OBSERVED z of OBS catalog = "<< maxzop <<endl;
		cout << endl;
	
	
		// Set cosmology
		cout << "0.1/ Initialise cosmology:"<<endl;
		double h = 0.71, OmegaM = 0.267804, OmegaL = 0.73;
		SimpleUniverse su(h, OmegaM, OmegaL);
		su.SetFlatUniverse_OmegaMatter();
		double OmegaB = su.OmegaBaryon();
		cout <<"    OmegaK = "<< su.OmegaCurv() <<", OmegaM="<< su.OmegaMatter();
		cout <<", OmegaL = "<< OmegaL <<", OmegaB = "<< OmegaB;
		cout <<", H0 = "<< su.H0() <<endl;
		
		
		// Calculate selection function
		RandomGenerator rg;
		string tmp="tmptmp";
		FitsInOutFile fos(tmp,FitsInOutFile::Fits_Create);
		Cat2Grid cat(dt, su, rg, fos, ZOCol, ZSCol);
		
		string OutRoot = "tmp";
		if (DoDebug)
			cat.SetDebugOutroot(OutRoot);
		
		string sffile;
		sffile = SFFileName+"_nofz.txt";
	
		ifstream inp;
		inp.open(sffile.c_str(), ifstream::in);
		inp.close();
		if(inp.fail()) {
			inp.clear(ios::failbit);
			cat.SaveSelecFunc(SFFileName, FullCat, ZFCol);
			// both SPECTRO-z and OBSERVED-z sf's are computed here
			}

		if( remove(tmp.c_str()) != 0 )
            cout << "Error deleting temporary file" << endl;
  
  
			
    }
    catch(PThrowable exc ) {
        cerr << "getsf.cc , Catched exception: \n" << exc.what() << endl;
        }
    catch(std::exception ex) {
        cerr << "getsf.cc , Catched exception ! " << (string)(ex.what()) << endl;
        }
    catch(...) {
        cerr << "getsf.cc , Catched ... ! " << endl;
        }

    cout << "--------------- getsf.cc / END --------------------- " << endl;
}
