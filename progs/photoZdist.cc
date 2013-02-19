// -*- LSST-C++ -*-
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

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

#include "geneutils.h"

/*  

 Output a distribution of photo-z for each bin in spectroscopic redshift
 
  
 A Abate - Oct 12
 
*/

void usage(void);
void usage(void) {
	cout << endl<<" Usage: photoZdist [...options...]" << endl<<endl;
	cout << " Output a distribution of photo-z for each bin in redshift, by"<<endl;
	cout << " default this is photometric redshift, choose -s to make it by"<<endl;
	cout << " spectroscopic redshift instead "<< endl << endl;
	
	cout << " Works with .fit files output from SDSS server: ZSCol=z and ZPCol=Column1"<<endl;
	cout << endl;
	cout << " -i : ObsCat: FITS file containing observed catalog"<<endl;
	cout << " -o : OutFile: Output filename"<<endl;
	cout << " -z : ZSCol,ZPCol: read spec-z,photo-z from columns labeled [ZSCOL,ZPCOL] resp."<<endl;
	cout << " -b : zMin,zMax,nz: [NZ] redshift bins between [ZMIN] and [ZMAX]. [Default: 0.,3.,5]"<<endl;
	cout << " -s : Bin by spectroscopic redshift instead "<<endl;
	cout << " -d : dzMin,dzMax,ndz: define histogram of (zp-zs)/(1+zs). [Default: -0.5,0.5,200]"<<endl;
	cout <<endl;
	}


int main(int narg, char *arg[]) {

	SophyaInit();
	FitsIOServerInit();
  
	cout << " ==== setting defaults ===="<<endl;
	string infile, outfile, cols;
	string ZSCol = "zs", ZPCol = "zp";
	bool isColSpec = false;
	double zMin = 0., zMax = 3.;
	int nBin = 5;
	bool isSpec = false;
	int nbin=200;
	double mindz=-0.5, maxdz=0.5;

	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hsi:o:z:b:d:")) != -1) {
	switch (c) {
		case 'i' :
			infile = optarg;
			break;
		case 'o' :
			outfile = optarg;
			break;
		case 'b' :
			sscanf(optarg,"%lf,%lf,%d",&zMin,&zMax,&nBin);
			break;
		case 'z' :
			cols = optarg;
			isColSpec = true;
			break;
	    case 's' :
	        isSpec = true;
	        break;
	    case 'd' :
			sscanf(optarg,"%lf,%lf,%d",&mindz,&maxdz,&nbin);
			break;
		case 'h' :
		    default :
			usage(); return -1;
		}
	}

    if (isColSpec) {
        // decode column names
       	string delim=",";
	    vector<string> results;
	    stringSplit(cols,delim,results);
	    ZSCol = results[0];
	    ZPCol = results[1];
    }

	cout << " ==== finished decoding command line arguments ==== "<<endl<<endl;
	
	cout << "     Printing command line arguments ... "<<endl<<endl;

    cout <<"     Reading catalog from "<< infile << endl;
    cout <<"     and writing photo-z distributions to "<< outfile << endl;
    cout <<"     Reading spectroscopic redshifts from column labeled "<< ZSCol;
    cout <<", reading photometeric redshifts from column labeled "<< ZPCol <<endl;
    cout <<"     Binning z into "<< nBin <<" bins between z = "<< zMin;
    cout <<" and z = "<< zMax <<endl;
    if (isSpec)
        cout <<"     using spectroscopic redshifts "<<endl;
	cout <<endl;
	
    try {
  
        // Histogram bin width
        double dZbin = (zMax - zMin)/(nBin - 1);

  
		// Read in catalog FITS file
		cout <<"     Read in observed catalog ..."<< infile <<endl;
		FitsInOutFile fin(infile, FitsInOutFile::Fits_RO);
		fin.MoveAbsToHDU(2);
		SwFitsDataTable dt(fin,512,false);
		DataTableRow row = dt.EmptyRow();
		sa_size_t ng = dt.NEntry();
		sa_size_t Izs = dt.IndexNom(ZSCol);
		sa_size_t Izp = dt.IndexNom(ZPCol);
		double minzo,maxzo,minzop,maxzop;
		dt.GetMinMax(Izs,minzo,maxzo); 
		dt.GetMinMax(Izp,minzop,maxzop); 
		cout <<"     "<<ng<<" galaxies in the catalog"<<endl;
		cout <<"     Min z of OBS catalog = "<<minzo<<", max z of OBS catalog = "<<maxzo<<endl;
		cout <<"     Min zp of OBS catalog = "<<minzop<<", max zp of OBS catalog = "<<maxzop<<endl;
		cout <<endl;
		
		// Set up histograms
		cout <<"     Set up histograms:"<<endl;
		vector<Histo*> histograms;
		
		double bW;
		cout <<"     "<< nbin <<" bins with "<< mindz <<" < dz < "<< maxdz<<endl;
		for (int i=0; i<nBin; i++) {
		    histograms.push_back(new Histo());
		    histograms[i]->ReSize(mindz,maxdz,nbin);
		    bW = histograms[i]->BinWidth();
		    }
		cout <<"     Bin width = "<< bW << endl;
		cout << endl;

		
		// Loop over catalog and add to Histogram 
		cout <<"     Loop over catalog and add to Histo "<<endl;
		for (sa_size_t i=0; i<ng; i++){
		
		    cout <<"     On galaxy "<< i+1 <<" of "<< ng;
			dt.GetRow(i,row);
			double zs=row[Izs];
			double zp=row[Izp];
			double dz = (zp - zs)/(1. + zs);
			//cout << dz << endl;
			
			// Find which histogram dz belongs in
			int index;
			if (isSpec)
			    index=(int)floor((zs-zMin)/dZbin);
			else 
			    index=(int)floor((zp-zMin)/dZbin);
		    //cout <<", index = "<< index <<" (max="<<nBin-1<<")"<< endl;
		    
		    if (index<nBin && index>=0)
			    histograms[index]->Add(dz);
			else {
			    cout <<", rejecting galaxy: z = "<<zs;
			    if (index<0)
			        cout <<" and WARNING! z<0!"<<endl;
			    }
			cout << endl;

            //cout <<" zs = "<<zs<<", index = "<<index<<endl;
			}
		cout << endl;
		
		cout <<"     Total number of galaxies was = "<< ng <<endl;
		long cnt=0;
		for(int j=0; j<nBin; j++){
	        long nInBin = histograms[j]->NEntries();
	        double bc = zMin + j*dZbin + dZbin/2.;
		    cout << "     Number in bin "<<j+1<<" with bin center = "<< bc <<" is "<< nInBin << endl;
		    cnt+=nInBin;
		    }
		cout <<"     Total in all redshift bins = "<<cnt<<endl;
	    cout << endl;
				
    // Write to text file
		
    ifstream inp;
    ofstream outp;
    inp.open(outfile.c_str(), ifstream::in);
    inp.close();
    if(inp.fail()) {
    
        inp.clear(ios::failbit);
		cout << "    Writing to file ..." << outfile.c_str() <<endl;
		outp.open(outfile.c_str(), ofstream::out);
			
		for (int i=0; i<nbin; i++) {
			
		    outp << histograms[0]->BinCenter(i) <<"  ";
			    
		    for(int j=0; j<nBin; j++)
			    outp << histograms[j]->operator()(i) <<" ";
		    outp << endl;
		    }
		outp.close();
		}
	else
		cout << "Error...file """ << outfile.c_str() << """ exists" << endl;
		
	cout << endl;
  }
  catch(PThrowable exc ) {
    cerr << "photoZdist.cc , Catched exception: \n" << exc.what() << endl;
  }
  catch(std::exception ex) {
    cerr << "photoZdist.cc , Catched exception ! " << (string)(ex.what()) << endl;
  }
  catch(...) {
    cerr << "photoZdist.cc , Catched ... ! " << endl;
  }

cout << "--------------- photoZdist.cc / END --------------------- " << endl;
}
