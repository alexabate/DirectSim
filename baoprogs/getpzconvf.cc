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
#include "cat2grid.h"

/*  

 Compute photo-z convolution function from the z-range
 given in a *_subinfo.txt file (produced by subfromfull)
 
 A Abate - Feb 2011
 
*/

void usage(void);
void usage(void)
	{
	cout << endl<<" Usage: getpzconf [...options...]" << endl<<endl;
		
	cout << "  Compute photo-z convolution function from the z-range"<<endl;
	cout << "  given in a *_subinfo.txt file (produced by subfromfull)"<<endl;
	cout << "  "<<endl;
		
	cout << " -O : ObsCat: FITS file containing observed catalog(s) (up to 2)"<<endl;
	cout << " -Z : SubInfoFile: *_subinfo.txt file to read in"<<endl;
	cout << " -o : PZConvFile: photo-z convolution function text filename"<<endl;
	cout << " -z : ZOCol,ZSCol: read OBSERVED redshifts from column labeled ZOCol, SPECTRO redshifts from column labeled ZSCol"<<endl;
	cout <<endl;
	}

void StringSplit(string str, string delim, vector<string>& results)
	{
	int cutAt;
	while( (cutAt = str.find_first_of(delim)) != str.npos )
		{
		if(cutAt > 0)
			results.push_back(str.substr(0,cutAt));
		str = str.substr(cutAt+1);
		}
		if(str.length() > 0)
			results.push_back(str);
		
	}

int main(int narg, char *arg[])
	{

	SophyaInit();
	FitsIOServerInit();
  
	cout << " ==== setting defaults ===="<<endl;
	string obscats,subinfo,pzconv;
	string ObCat1,ObCat2;
	int NC=0;
	string ZCol; 
	bool Zsp=false;
	string ZOCol = "zp"; // by default OBSERVED redshift column labelled "zp" is read in
	string ZSCol = "zs"; // by default SPECTRO redshift column labelled "zs" is read in

	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hO:Z:o:z:")) != -1) 
	{
	switch (c) 
		{
		case 'O' :
			obscats = optarg;
			break;
		case 'Z' :
			subinfo = optarg;
			break;
		case 'o' :
			pzconv = optarg;
			break;
		case 'z' :
			ZCol = optarg;
			break;
		case 'h' :
		default :
			usage(); return -1;
		}
	}
	

	string delim=",";
	vector<string> results;
	StringSplit(obscats,delim,results);
	vector<string>::iterator i;
	i = results.begin();
	ObCat1=*i;
	i++;
	NC++;
	if (i!=results.end())
		{
		ObCat2=*i;
		i++;
		NC++;
		}
		
	// get two z column names
	if (Zsp)	
		{ 
		string delim=",";
		vector<string> results;
		StringSplit(ZCol,delim,results);
		vector<string>::iterator i;
		i = results.begin();
		ZOCol=*i;
		i++;
		if (i!=results.end())
			ZSCol=*i;
		}
		
	cout << " ==== finished decoding command line arguments ==== "<<endl<<endl;
	
	cout << "     Printing command line arguments ... "<<endl<<endl;
	
	cout << "     Reading in "<<NC<<" observed catalogs, filename";
	if (NC<2)
		cout <<" is "<<ObCat1<<endl;
	else
		cout <<"s are "<<ObCat1<<" and "<<ObCat2<<endl;

	cout << "     Photometric redshifts to be read from "<<ZOCol<<" column, spectroscopic from "<<ZSCol<<" column"<<endl;
	cout << "     Z range will be read from "<<subinfo<<endl;
	cout << "     Photo-z convolution function will be saved to "<<pzconv<<endl;
	cout <<endl;
	
  try {
		/* READ IN SUBINFO FILE*/
		cout <<"0/ Read in subinfo file "<<subinfo<<endl;
		ifstream ifs(subinfo.c_str());
		Array B;
		sa_size_t nr, nc;
		B.ReadASCII(ifs,nr,nc);
		double minz = B(0,0);
		double zc = B(1,0);
		double maxz = B(2,0);
		cout <<"    Sub-array bounds: "<<minz<<"< z <"<<maxz<<endl;
		cout <<"    Sub-array center: zc = "<<zc<<endl;

		/* INITIALISE HISTOGRAM */
		int nbin=200;
		double mindz=-2, maxdz=2;
		Histo deltaz(mindz, maxdz, nbin);
		cout <<endl;
  
		/* READ IN CATALOG FITS FILE */
		cout <<"1/ Read in observed catalog "<<ObCat1<<endl;
		FitsInOutFile fin(ObCat1,FitsInOutFile::Fits_RO);
		fin.MoveAbsToHDU(2);
		SwFitsDataTable dt(fin,512,false);
		DataTableRow row = dt.EmptyRow();
		sa_size_t ng = dt.NEntry();
		sa_size_t Izs = dt.IndexNom(ZSCol);
		sa_size_t Izp = dt.IndexNom(ZOCol);
		double minzo,maxzo,minzop,maxzop;
		dt.GetMinMax(Izs,minzo,maxzo); 
		dt.GetMinMax(Izp,minzop,maxzop); 
		cout <<"    "<<ng<<" galaxies in the catalog"<<endl;
		cout <<"    Min z of OBS catalog = "<<minzo<<", max z of OBS catalog = "<<maxzo<<endl;
		cout <<"    Min zp of OBS catalog = "<<minzop<<", max zp of OBS catalog = "<<maxzop<<endl;
		cout <<endl;
		
		/* LOOP OVER CATALOG FITS FILE AND ADD TO HISTO */
		cout <<"2/ Loop over catalog and add to Histo "<<endl;
		for (sa_size_t i=0; i<ng; i++)
			{
			dt.GetRow(i,row);
			double zs=row[Izs];
			double zp=row[Izp];
			double dz = zp - zs;
			if (zs>=minz&&zs<maxz)
				deltaz.Add(dz);
			}
			
		if (NC>1)
			{
			/* READ IN 2nd CATALOG FITS FILE */
			cout <<"1a/ Read in 2nd observed catalog "<<ObCat2<<endl;
			FitsInOutFile fin2(ObCat2,FitsInOutFile::Fits_RO);
			fin2.MoveAbsToHDU(2);
			SwFitsDataTable dt2(fin2,512,false);
			DataTableRow row2 = dt2.EmptyRow();
			sa_size_t ng2 = dt2.NEntry();
			sa_size_t Izs2 = dt2.IndexNom("z");
			sa_size_t Izp2 = dt2.IndexNom(ZCol);
			double minzo2,maxzo2,minzop2,maxzop2;
			dt2.GetMinMax(Izs2,minzo2,maxzo2); 
			dt2.GetMinMax(Izp2,minzop2,maxzop2); 
			cout <<"    "<<ng2<<" galaxies in the catalog"<<endl;
			cout <<"    Min z of OBS catalog = "<<minzo2<<", max z of OBS catalog = "<<maxzo2<<endl;
			cout <<"    Min zp of OBS catalog = "<<minzop2<<", max zp of OBS catalog = "<<maxzop2<<endl;
			cout <<endl;
			
			cout <<"1b/ Loop over 2nd catalog and add to Histo "<<endl;
			/* LOOP OVER CATALOG FITS FILE AND ADD TO HISTO */
			for (sa_size_t i=0; i<ng2; i++)
				{
				dt2.GetRow(i,row2);
				double zs=row2[Izs2];
				double zp=row2[Izp2];
				double dz = zp - zs;
				if (zs>=minz&&zs<maxz)
					deltaz.Add(dz);
				}
			
			}
		
		/* WRITE HISTO TO TEXT FILE*/
		
		ifstream inp;
		ofstream outp;
		inp.open(pzconv.c_str(), ifstream::in);
		inp.close();
		if(inp.fail())
			{
			inp.clear(ios::failbit);
			cout << "    Writing to file ..." << pzconv.c_str() << " and "<<pzconv<<endl;
			outp.open(pzconv.c_str(), ofstream::out);
			for(int_4 i=0;i<nbin;i++)
				{
			r_8 bc=deltaz.BinCenter(i);
			r_8 nz=deltaz.operator()(i);
		
			outp <<bc<<"      "<<nz<<endl;
			}
		outp.close();
		}
	else
		cout << "Error...file """ << pzconv.c_str() << """ exists" << endl;
		
	cout << endl;
  }
  catch(PThrowable exc ) {
    cerr << "getpzconf.cc , Catched exception: \n" << exc.what() << endl;
  }
  catch(std::exception ex) {
    cerr << "getpzconf.cc , Catched exception ! " << (string)(ex.what()) << endl;
  }
  catch(...) {
    cerr << "getpzconf.cc , Catched ... ! " << endl;
  }

cout << "--------------- getpzconf.cc / END --------------------- " << endl;
}
