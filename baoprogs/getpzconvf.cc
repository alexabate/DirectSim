/**
 * @file  getpzconvf.cc
 * @brief Compute the photo-z convolution function given the photo-z - true-z 
 *        distribution within a given redshift range
 *
 * @todo histogram bins should be option to program
 * @todo not sure if zp-zs or (zp-zs)/(1+zs) should be added to histogram?
 * 
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: Feb 2011
 * @date 2011
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

// SOPHYA
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


void usage(void);
void usage(void) {
	
	cout << endl<<" Usage: getpzconvf [...options...]              "<<endl<<endl;
		
	cout << "  Compute photo-z convolution function from the redshift   "<<endl;
	cout << "  range given in a *_subinfo.txt file (produced by the     "<<endl;
	cout << "  program subfromfull)                                     "<<endl;
	cout << endl;
		
	cout << " -O : obs_cat_file: FITS file containing observed catalog(s) "<<endl;
	cout << "                    (up to 2 separated by commas)            "<<endl;
	cout << " -Z : z_range_file: *_subinfo.txt file to read in            "<<endl;
	cout << " -o : out_file: file to write photo-z convolution function to"<<endl;
	cout << " -z : obs_z_col,true_z_col: read OBSERVED redshifts from     "<<endl;
	cout << "      column labeled obs_z_col, true redshifts from column   "<<endl;
	cout << "      labeled true_z_col                                     "<<endl;
	cout << endl;
	}



int main(int narg, char *arg[]) {

	SophyaInit();
	FitsIOServerInit();
  
	string obscats;        // observed catalog(s) to read in (up to 2)
	string z_range_file;   // *_subinfo.txt file containing redshift range
	string outfile;        // output file
	string obs_cat1;       // first catalog to read in
	string obs_cat2;       // second catalog to read in (if specified)
	int num_cats = 0;      // number of catalogs to read in
	string zcols;          // string to get column names from 'z' argument
	bool getColNames = false; // if true new redshift column names were specifed
	string obs_z_col = "zp";  // by default OBSERVED z column called "zp" is read 
	string true_z_col = "zs"; // by default SPECTRO z column called "zs" is read 

	//--- decoding command line arguments 
	char c;
	while((c = getopt(narg,arg,"hO:Z:o:z:")) != -1) {
	    switch (c) {
		    case 'O' :
			    obscats = optarg;
			    break;
		    case 'Z' :
			    z_range_file = optarg;
			    break;
		    case 'o' :
			    outfile = optarg;
			    break;
		    case 'z' :
			    zcols = optarg;
			    getColNames = true;
		 	    break;
		    case 'h' :
		    default :
			    usage(); return -1;
		    }
	    }
	

    // get catalog(s) names
	string delim=",";
	vector<string> cats;
	stringSplit(obscats, delim, cats);
	obs_cat1 = cats[0];
	if (cats.size()>1) {
	    obs_cat2 = cats[1];
	    num_cats++;
	    }

		
	// get two z column names
	if (getColNames) {
		vector<string> results;
		stringSplit(zcols, delim, results);
		obs_z_col = results[0];
		if (results.size()>1)
		    true_z_col = results[1];
		}
		
	
	cout << "     Printing command line arguments ...             "<<endl<<endl;
	
	cout << "     Reading in "<< num_cats <<" observed catalog(s), from file";
	if (num_cats<2)
		cout << obs_cat1 <<endl;
	else
		cout <<"s "<< obs_cat1 <<" and "<< obs_cat2 <<endl;

	cout << "     Photometric redshifts will be read from "<< obs_z_col;
	cout << " column, true redshifts from "<< true_z_col <<" column     "<<endl;
	cout << "     Redshift range will be read from "<< z_range_file      <<endl;
	cout << "     Photo-z convolution function will be saved to "<< outfile <<endl;
	cout << endl;
	
	
    try {
		
		
		// Read in redshift range from file
		cout <<"     Read in redshift range from file "<< z_range_file <<endl;
		ifstream ifs(z_range_file.c_str());
		Array B;
		sa_size_t nr, nc;
		B.ReadASCII(ifs,nr,nc);
		double minz = B(0,0);
		double zc = B(1,0);
		double maxz = B(2,0);
		cout <<"    Redshift bounds: "<< minz <<"< z <"<< maxz <<endl;
		cout <<"    Redshift center: zc = "<< zc <<endl;


		// Initialize histogram
		// @todo histogram bins should be option to program
		int nbin = 200;
		double mindz = -2, maxdz = 2;
		Histo deltaz(mindz, maxdz, nbin);
		cout << endl;
  
  
		// Read in observed catalog
		cout <<"     Read in observed catalog "<< obs_cat1;
		if (num_cats>1)
		    cout << " 1 of 2" <<endl;
		else
		    cout << endl;
		FitsInOutFile fin(obs_cat1, FitsInOutFile::Fits_RO);
		fin.MoveAbsToHDU(2);
		SwFitsDataTable dt(fin,512,false);
		DataTableRow row = dt.EmptyRow();
		sa_size_t ng = dt.NEntry();
		sa_size_t Izs = dt.IndexNom(true_z_col);
		sa_size_t Izp = dt.IndexNom(obs_z_col);
		double minzo,maxzo,minzop,maxzop;
		dt.GetMinMax(Izs,minzo,maxzo); 
		dt.GetMinMax(Izp,minzop,maxzop); 
		cout <<"     "<< ng <<" galaxies in the catalog"<<endl;
		cout <<"     Min true z of OBS catalog = "<< minzo;
		cout <<", max true z of OBS catalog = "<< maxzo <<endl;
		cout <<"     Min photo-z of OBS catalog = "<< minzop;
		cout <<", max photo-z of OBS catalog = "<< maxzop <<endl;
		cout << endl;
		
		
		// Loop over catalog and add redshift difference to histogram
		cout <<"     Loop over catalog and add to Histo "<<endl;
		for (sa_size_t i=0; i<ng; i++) {
			
			dt.GetRow(i,row);
			double zs = row[Izs];
			double zp = row[Izp];
			double dz = zp - zs; // should this be /(1+zs)??
			
			// if redshift is within redshift range that was read in
			if ( zs>=minz && zs<maxz )
				deltaz.Add(dz);
			}
			
			
	    // Now check if there's a second catalog and repeat
		if (num_cats>1) {
			
			// Read in 2nd catalog
			cout <<"     Read in observed catalog "<< obs_cat2 <<" 2 of 2"<<endl;
			FitsInOutFile fin2(obs_cat2, FitsInOutFile::Fits_RO);
			fin2.MoveAbsToHDU(2);
			SwFitsDataTable dt2(fin2,512,false);
			DataTableRow row2 = dt2.EmptyRow();
			sa_size_t ng2 = dt2.NEntry();
			sa_size_t Izs2 = dt2.IndexNom(true_z_col);
		    sa_size_t Izp2 = dt2.IndexNom(obs_z_col);
			double minzo2,maxzo2,minzop2,maxzop2;
			dt2.GetMinMax(Izs2,minzo2,maxzo2); 
			dt2.GetMinMax(Izp2,minzop2,maxzop2); 
			cout <<"     "<< ng2 <<" galaxies in catalog 2"<<endl;
			cout <<"     Min  true z of OBS catalog = "<< minzo2;
			cout <<", max true z of OBS catalog = "<< maxzo2 <<endl;
			cout <<"     Min photo-z of OBS catalog = "<< minzop2;
			cout <<", max photo-z of OBS catalog = "<< maxzop2 <<endl;
			cout << endl;
			
			
			// Loop over catalog and add redshift difference to histogram
			cout <<"     Loop over catalog 2 and add to Histo "<<endl;
			for (sa_size_t i=0; i<ng2; i++) {
				
				dt2.GetRow(i,row2);
				double zs=row2[Izs2];
				double zp=row2[Izp2];
				double dz = zp - zs;// should this be /(1+zs)??
			
			    // if redshift is within redshift range that was read in
				if ( zs>=minz && zs<maxz )
					deltaz.Add(dz);
				}
			
			}
		
		
		// Write histogram to a file
		ifstream inp;
		ofstream outp;
		inp.open(outfile.c_str(), ifstream::in);
		inp.close();
		if (inp.fail()) {
			
			inp.clear(ios::failbit);
			cout << "     Writing to file ..." << outfile.c_str() << endl;
			outp.open(outfile.c_str(), ofstream::out);
			
			for (int_4 i=0; i<nbin; i++) {
				
			    r_8 bc = deltaz.BinCenter(i);
			    r_8 nz = deltaz.operator()(i);
		
			    outp << bc <<"  "<< nz <<endl;
			    }
		    outp.close();
		    }
	    else
		    cout << "Error...file """ << outfile.c_str() << """ exists" << endl;
		
	    cout << endl;
      }
      

catch(PThrowable exc ) {
    cerr << "getpzconvf.cc , Catched exception: \n" << exc.what() << endl;
    }
catch(std::exception ex) {
    cerr << "getpzconvf.cc , Catched exception ! " << (string)(ex.what()) << endl;
    }
catch(...) {
    cerr << "getpzconvf.cc , Catched ... ! " << endl;
    }

cout << "--------------- getpzconvf.cc / END --------------------- " << endl;
}
