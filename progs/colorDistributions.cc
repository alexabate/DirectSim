#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

// sophya libraries
//#include "histinit.h"
#include "fitsioserver.h"
#include "swfitsdtable.h"
#include "fiosinit.h"
#include "mydefrg.h"

// CatSim classes
#include "simdata.h"
#include "sedfilter.h"

#define PI 3.141592

void usage(void);
void usage(void) {
	cout << endl<<" Usage: colorDistributions [...options...]" << endl<<endl;
	
	cout << " Makes histograms of colors vs redshift from an input catalog of "<<endl;
	cout << " of observational data.  The observational data is ugriz[y] AB "<<endl;
	cout << " magnitudes plus redshifts. "<<endl<<endl;
	
	cout << " Color vs redshift tracks for the spectra given in the supplied SED "<<endl;
	cout << " list file are also calculated"<<endl<<endl;
	
	cout << " The input catalog file must be a FITS bintable and you must know the "<<endl;
	cout << " column names of the ugriz[y] and redshift columns.  If a column name "<<endl;
	cout << " for the y-filter is NOT supplied it is assumed there are only ugriz "<<endl;
	cout << " magnitudes in the catalog "<<endl<<endl;

	cout << " -i INCAT:       Name of file containing observed catalog (must be FITS) "<<endl;
	cout << " -o OUTFILE:     Name of file to output distributions to (will be saved in output/) "<<endl;
	cout << " -s SPECTRAFILE: Name of SED list file [DEFAULT=CWWKSB.list]"<<endl;
	cout << " -f FILTERFILE:  Name of filter list file [DEFAULT=LSST.filters]"<<endl;
	cout << " -c COLS:        Column names of zs, ugriz or ugrizy magnitudes"<<endl;
	cout << " -l ZL,CL:       set the lower limits of the redshift (ZL) and "<<endl;
	cout << "                 color (CL) bins [DEFAULT ZL=0,CL=-2]"<<endl;
	cout << " -u ZU,CU:       set the upper limits of the redshift (ZU) and "<<endl;
	cout << "                 color (CU) bins [DEFAULT ZU=4,CU=3]"<<endl;
	cout << " -n NZ,NC:       set the number of redshift (NZ) and color (NC)"<<endl;
	cout <<"                  bins [DEFAULT NZ=50,NC=50]"<<endl;
	cout << endl;
    }
    
int main(int narg, char* arg[]) {

    cout << " ==== colorDistributions.cc program ==== "<<endl;

	// make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	cout<<endl<<endl;

	//--- decoding command line arguments 
    string infile, outfileroot;
    string spectraFile = "CWWKSB.list";
    string filterFile = "LSST.filters";
    string zsCol="z";
    string uCol="obs_umag",gCol="obs_gmag",rCol="obs_rmag",iCol="obs_imag";
    string zCol="obs_zmag",yCol="obs_ymag";
    string cols;
    double zLower=0., zUpper=4.;
    double cLower=-2., cUpper=3.;
    int nZbins=50, nCbins=50;
    int nFilter = 6;
  
	//--- decoding command line arguments 
	cout << " ==== decoding command line arguments ===="<<endl;
	char c;
	while((c = getopt(narg,arg,"hi:o:s:f:c:l:u:n:")) != -1) {
	    switch (c) {
	        case 'i' :
	            infile = optarg;
	            break;
	        case 'o' :
	            outfileroot = optarg;
	            outfileroot = "output/" + outfileroot;
	            break;
	        case 's' :
	            spectraFile = optarg;
	            break;
	        case 'f' :
	            filterFile = optarg;
	            break;
	        case 'c' :
	            cols = optarg;
	            break;
	        case 'l' :
		        sscanf(optarg,"%lf,%lf",&zLower,&cLower);
		        break;
		    case 'u' :
		        sscanf(optarg,"%lf,%lf",&zUpper,&cUpper);
		        break;
		    case 'n' :
		        sscanf(optarg,"%d,%d",&nZbins,&nCbins);
		        break;
	        case 'h' :
		        default :
		        usage(); return -1;
		    }
	    }
	    
	string delim=",";
	vector<string> results;
	stringSplit(cols,delim,results);
    nFilter = results.size()-1;

    zsCol=results[0];
    uCol=results[1];
    gCol=results[2];
    rCol=results[3];
    iCol=results[4];
    zCol=results[5];
	if (nFilter>5)
	    yCol = results[6];
	       
    //-- end command line arguments
    cout <<"     Reading observations catalog "<< infile <<endl;
    cout <<"     Reading magnitudes from columns named: "<< uCol <<", "<< gCol;
    cout <<", "<< rCol <<", "<< iCol <<", "<< zCol;
    if (nFilter>5)
        cout <<", "<< yCol<<endl;
    else
        cout << endl;
    cout <<"     Writing out color distributions and SED color-redshift tracks "<<endl;
    cout <<"     to files beginning "<< outfileroot << endl;
    cout <<"     Binning colors by "<< nCbins <<" bins between "<< cLower <<" and "<< cUpper <<endl;
    cout <<"     Binning redshifts by "<< nZbins <<" bins between "<< zLower <<" and "<< zUpper <<endl;
    cout <<"     Reading spectra from files listed in "<< spectraFile <<endl;
    cout << endl;
  
    int rc = 1;  
    try {  // exception handling try bloc at top level
  
    ifstream inp;
	ofstream outp;
	string outfile;
	
	// Read in the observed catalog file
	cout <<"     Reading in file "<< infile <<endl;
	FitsInOutFile fin(infile, FitsInOutFile::Fits_RO);
	fin.MoveAbsToHDU(2);
	SwFitsDataTable dt(fin,512,false);
	cout <<endl;
	DataTableRow rowIn = dt.EmptyRow();
	sa_size_t ng=dt.NEntry();
	sa_size_t nc=dt.NCols();
	DataTableRow row=dt.EmptyRow();
	
	// Look at file contents
	cout <<"     In file "<< infile <<" ... "<<endl;
	cout <<"     Number of columns = "<<nc<<", number of entries = "<< ng << endl;
	cout <<"     Columns in the file are:"<<endl;
	cout <<"     #    Name "<<endl;
	for (int i=0; i<nc; i++)
	    cout << "     "<<i<<"    "<<dt.NomIndex(i)<<endl;
	cout << endl;
	
	// Get column indices
	sa_size_t iZs = dt.IndexNom(zsCol);
	sa_size_t iU = dt.IndexNom(uCol);
	sa_size_t iG = dt.IndexNom(gCol);
	sa_size_t iR = dt.IndexNom(rCol);
	sa_size_t iI = dt.IndexNom(iCol);
	sa_size_t iZ = dt.IndexNom(zCol);
	sa_size_t iY;
	if (nFilter>5)
	    iY = dt.IndexNom(yCol);
	
	// Look at ranges of magnitude and redshift columns
    double mMin,mMax;
    dt.GetMinMax(dt.ColumnIndex(uCol),mMin,mMax);
    cout <<"     Range of "<<dt.NomIndex(iU)<<" ("<<iU<<"th) column is "<<mMin<<"<m<"<<mMax<<endl;
    dt.GetMinMax(dt.ColumnIndex(gCol),mMin,mMax);
    cout <<"     Range of "<<dt.NomIndex(iG)<<" ("<<iG<<"th) column is "<<mMin<<"<m<"<<mMax<<endl;
    dt.GetMinMax(dt.ColumnIndex(rCol),mMin,mMax);
    cout <<"     Range of "<<dt.NomIndex(iR)<<" ("<<iR<<"th) column is "<<mMin<<"<m<"<<mMax<<endl;
    dt.GetMinMax(dt.ColumnIndex(iCol),mMin,mMax);
    cout <<"     Range of "<<dt.NomIndex(iI)<<" ("<<iI<<"th) column is "<<mMin<<"<m<"<<mMax<<endl;
    dt.GetMinMax(dt.ColumnIndex(zCol),mMin,mMax);
    cout <<"     Range of "<<dt.NomIndex(iZ)<<" ("<<iZ<<"th) column is "<<mMin<<"<m<"<<mMax<<endl;
    if (nFilter>5) {
        dt.GetMinMax(dt.ColumnIndex(yCol),mMin,mMax);
        cout <<"     Range of "<<dt.NomIndex(iY)<<" ("<<iY<<"th) column is "<<mMin<<"<m<"<<mMax<<endl;
        }
    double zMin,zMax;
    dt.GetMinMax(dt.ColumnIndex(zsCol),zMin,zMax);
    cout <<"     Range of "<<dt.NomIndex(iZs)<<" ("<<iZs<<"th) column is "<<zMin<<"<z<"<<zMax<<endl;
    cout << endl;
    
    // Create 2D histogram of i band magnitude and redshifts
    cout <<"     Creating 2D histograms of colors vs redshifts ... "<<endl;
    
    // Define bins
    double dZbin = (zUpper - zLower)/(nZbins - 1);
    double dCbin = (cUpper - cLower)/(nCbins - 1);
    
    // Set up histograms
	cout <<"     Set up histograms: "<<endl;
	vector<TArray<int>*> colorZ;
	for (int i=0; i<nFilter-1; i++){
    
        int nDim = 2;
        sa_size_t myDim[nDim];
        myDim[0]=nCbins; myDim[1]=nZbins;
        colorZ.push_back(new TArray<int>());
        colorZ[i]->SetSize(nDim,myDim);
        
        }
    
    // Loop over catalog putting each data point in correct histogram
    for (long ig=0; ig<ng; ig++) {
    
        // Get the data
        dt.GetRow(ig,row);
        double umag = row[iU];
        double gmag = row[iG];
        double rmag = row[iR];
        double imag = row[iI];
        double zmag = row[iZ];
        double ymag;
        if (nFilter>5)
            ymag = row[iY];
        double zs = row[iZs];
        
        // The colors
        double cUG = umag - gmag;
        double cGR = gmag - rmag;
        double cRI = rmag - imag;
        double cIZ = imag - zmag;
        double cZY;
        if (nFilter>5)
            cZY = zmag - ymag;
        
        //cout << cUG <<"  "<< cGR <<"  "<< cRI <<"  "<< cIZ <<"  ";
        
        // The histo bins
        sa_size_t c1Index = (sa_size_t)floor((cUG - cLower)/dCbin);
        sa_size_t c2Index = (sa_size_t)floor((cGR - cLower)/dCbin);
        sa_size_t c3Index = (sa_size_t)floor((cRI - cLower)/dCbin);
        sa_size_t c4Index = (sa_size_t)floor((cIZ - cLower)/dCbin);
        sa_size_t c5Index;
        if (nFilter>5)
            c5Index = (sa_size_t)floor((cZY - cLower)/dCbin);
        sa_size_t zIndex = (sa_size_t)floor((zs - zLower)/dZbin);
        
        // Add to correct histo
        if ( c1Index<nCbins && zIndex<nZbins && c1Index>0 && zIndex>0)// if within bins
             colorZ[0]->operator()(c1Index,zIndex)++;
        if ( c2Index<nCbins && zIndex<nZbins && c2Index>0 && zIndex>0)// if within bins
             colorZ[1]->operator()(c2Index,zIndex)++;
        if ( c3Index<nCbins && zIndex<nZbins && c3Index>0 && zIndex>0)// if within bins
             colorZ[2]->operator()(c3Index,zIndex)++;
        if ( c4Index<nCbins && zIndex<nZbins && c4Index>0 && zIndex>0)// if within bins
             colorZ[3]->operator()(c4Index,zIndex)++;
        if (nFilter>5) {
            if ( c5Index<nCbins && zIndex<nZbins && c5Index>0 && zIndex>0)// if within bins
                 colorZ[4]->operator()(c5Index,zIndex)++;
            }

                
        }
        
    for (int ic=0; ic<nFilter-1; ic++){    
        stringstream ss;
        ss << ic+1;
        outfile = outfileroot + "_colordist"+ss.str()+".txt";
        outp.open(outfile.c_str());
        cout <<"     Writing color histogram "<<ic+1<<" to file "<< outfile << endl;
        for (int i=0; i<nCbins; i++){
            for (int j=0; j<nZbins; j++)
                outp << colorZ[ic]->operator()(i,j) << "  ";
            outp << endl;
            }
        outp.close();
        cout << endl;
        }
        
    // wavelength range of the SEDs/filters
	double lmin=5e-8, lmax=2.5e-6;
	
	// Photometry calculation class
    PhotometryCalcs photometryCalcs(lmin,lmax);
        
    // Read in galaxy spectra
    ReadSedList readSedList(spectraFile);
    readSedList.readSeds(lmin,lmax);
    vector<SED*> sedArray=readSedList.getSedArray();
    int nsed=readSedList.getNSed();
    cout <<"     Number of spectra = "<<nsed<<endl;
	cout << endl;
	
	// Read in filters
	ReadFilterList readFilterList(filterFile);
	readFilterList.readFilters(lmin,lmax);
	vector<Filter*> filters = readFilterList.getFilterArray();
	int nFilterRead=readFilterList.getNTot();
	cout <<"     "<< nFilterRead <<" filters read in "<<endl;
        
    // Calculate template color tracks
    int nz = 1000;
    double zmin = 0., zmax = zUpper;
    double dz = (zmax-zmin)/(nz-1);
    
    for (int is=0; is<nsed; is++) {
    
        int nDim = 2;
        sa_size_t myDim[nDim];
        myDim[0]=5; myDim[1]=nz;
        TArray<double> colorTracks;
        colorTracks.SetSize(nDim,myDim);
        
        for (int i=0; i<nz;i++) {
		    double z=zmin+i*dz;

            for (int j=0; j<nFilter-1; j++) {
		        colorTracks(j,i) = photometryCalcs.CompColor(z,(*sedArray[is]),(*filters[j]),(*filters[j+1]));
		        }
		        
		    }
		
		stringstream ss;
		ss << is;
		outfile = outfileroot + "_colortrackSED"+ss.str()+".txt";
        outp.open(outfile.c_str());
        cout <<"     Writing color track for spectrum "<< is <<" to file "<< outfile << endl;
        for (int i=0; i<nz; i++){
            for (int j=0; j<5; j++)
                outp << colorTracks(j,i) << "  ";
            outp << endl;
            }
        outp.close();
        cout << endl;
		    
		}
  
	}  // End of try bloc 
   
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " colorDistributions.cc: Catched Exception (PThrowable)" 
	 << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " colorDistributions.cc: Catched std::exception "  << " - what()= " 
	 << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " colorDistributions.cc: some other exception (...) was caught ! " 
	 << endl;
    rc = 97;
  }
  cout << " ==== End of colorDistributions.cc program  Rc= " << rc << endl;
  return rc;	
}

