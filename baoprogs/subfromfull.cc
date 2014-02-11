/**
  * @file  subfromfull.cc
  * @brief Read in either gridded galaxy data or a SimLSS grid and output a sub grid
  *
  * @todo add cosmology into header (read from read in file's header)
  * @todo include in header if catalog has selection effects (both in file read
  *       in and file written out
  * @todo remove debug part? doesn't really do anything
  *
  * @author Alex Abate
  * Contact: abate@email.arizona.edu
  *
  */

#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

// sophya
#include "machdefs.h"
#include "sopnamsp.h"
#include "timing.h"
#include "array.h"
#include "hisprof.h"
#include "histerr.h"
#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "swfitsdtable.h"
#include "resusage.h"

// DirectSim/
//#include "geneutils.h"
//#include "cat2grid.h"
//#include "powerspec.h"
//#include "mass2gal.h"
//#include "pkspectrum.h"
//#include "fitkbaoscale.h"
//#include "chisqstats.h"


void usage(void);
void usage(void) {

	cout << endl<<" Usage: subfromfull [...options...]           " <<endl<<endl;
	
	cout << "  Read in either a file containing gridded galaxy data     "<<endl;
	cout << "  (output from grid_data program) or a SimLSS grid of      "<<endl;
	cout << "  overdensities (output from simdensity program) using     "<<endl;
	cout << "  option -F, and write out a sub-portion of these grid(s)  "<<endl;
	cout << "  to a new file specified via the -O option.               "<<endl;
	cout << endl;
	
	cout << "  To specify the sub-portion of these grid(s) give its     "<<endl;
	cout << "  pixel ranges in the following way via option -x (range of"<<endl;
	cout << "  sub-grid pixels in x-dimension), -y (range of sub-grid   "<<endl;
	cout << "  pixels in y-dimension), -z (range of sub-grid pixels in  "<<endl;
	cout << "  z-dimension).                                            "<<endl;
	cout << endl;
		
	cout << " -F : infile : file containing gridded galaxy data or SimLSS grid"<<endl;
	cout << " -O : outfile : filename to write sub-grids to             "<<endl;
	cout << " -x : x1,x2 : range of sub-grid pixels in x-dimension      "<<endl;
	cout << " -y : y1,y2 : range of sub-grid pixels in y-dimension      "<<endl;
	cout << " -z : z1,z2 : range of sub-grid pixels in z-dimension      "<<endl;
	cout << " -s : [noarg] : file is SimLSS density grid                "<<endl;
	}


int main(int narg, char* arg[]) {

	cout << " ==== subfromfull.cc program  ==== " <<endl;
	// Make sure SOPHYA modules are initialized 
	SophyaInit();  
	FitsIOServerInit();
	InitTim();
	
	// Set defaults etc ....
	// FILES TO READ IN/OUT
	string infile, outfile;
	// SUB ARRAY PARAMETERS
	sa_size_t x1=0, x2=0;	// Range of sub array in x-dim
	sa_size_t y1=0, y2=0;	// Range of sub array in y-dim
	sa_size_t z1=0, z2=0;	// Range of sub array in z-dim
	bool isSimLSS=false;
	// DEBUG
	bool DeBug=false;
	
	string dbfile="tmp";
	
	
	//--- decoding command line arguments 
	char c;
	while((c = getopt(narg,arg,"hsF:O:x:y:z:d:")) != -1) {
	    switch (c) {
	        case 'F' :
		        infile = optarg;
		        break;
	        case 'O' :
		        outfile	= optarg; // filename of subfromfulls
		        break;
	        case 'x' :
		        sscanf(optarg,"%ld,%ld",&x1,&x2);
		        break;
	        case 'y' :
		        sscanf(optarg,"%ld,%ld",&y1,&y2);
		        break;
	        case 'z' :
		        sscanf(optarg,"%ld,%ld",&z1,&z2);
		        break;
	        case 's' :
		        isSimLSS=true;
		        break;
	        case 'd' :
		        dbfile = optarg;
		        DeBug = true;
		        break;
	        case 'h' :
	            default :
		        usage(); return -1;
	        }
        }

		
	if ( x1<1 && x2<1 && y1<1 && y2<1 && z1<1 && z2<1)
		throw ParmError("ERROR! Must specify all sub-array pixel ranges");


	cout << "     Printing command line arguments ... "<<endl<<endl;
	
	// IN FILES
	cout << "     *GRID DETAILS*"<<endl;
	if(!isSimLSS)
		cout << "     Gridded data arrays read from "<< infile <<endl;
	else
		cout << "     Reading SimLSS cube from "<< infile;
	cout << endl;
	
	
	// SUB ARRAY STUFF
	cout << "     *SUB-ARRAY DETAILS*"<<endl;
	cout << "     Sub array has pixel ranges: (x1,x2)=("<< x1 <<","<< x2 <<")"<<endl;
	cout << "                                 (y1,y2)=("<< y1 <<","<< y2 <<")"<<endl;
	cout << "                                 (z1,z2)=("<< z1 <<","<< z2 <<")"<<endl;
	cout << endl;
	
	// OUTPUT FILES
	cout << "     *OUTPUT DETAILS*"<<endl;
	cout << "     Sub-array will be output to "<< outfile <<endl;
	cout << endl;
	//-- end command line arguments
  
	int rc = 1;  
	try {  // exception handling try bloc at top level
	
	ResourceUsage res;
	res.Update();
	cout << " Memory size (KB):" << res.getMemorySize() << endl;
	cout << " Max memory size (KB):" << res.getMaxMemorySize() << endl;
	cout << " Maximum allowed data segment size (KB):"<<res.getMaxDataSize()<<endl;
	cout << " Resource usage info : \n" << res << endl;
	
	
	// IF input file is NOT a SimLSS density grid
	if(!isSimLSS){
   	
   	    // Read in gridded galaxy data 
	    TArray<r_8> ngal;
	    TArray<r_8> wrgal, wngal, zc;
	
	    cout <<"     Read in gridded galaxy data from file "<< infile <<endl;
	    FitsInOutFile fin(infile, FitsInOutFile::Fits_RO);
	    fin >> ngal;
	    fin >> wngal; // this is same as ngal if no selection function effects
	    fin >> wrgal;
	    fin >> zc; // array of redshifts at each pixel value

        
        // Read key values from the file header
	    //string Rs = fin.KeyValue("DX");
	    double R = atof(fin.KeyValue("DX").c_str()); // size of grid cell
	    //string mean_denss=fin.KeyValue("MeanOverDensity");
	    double mean_dens = atof(fin.KeyValue("MeanOverDensity").c_str()); 
	    cout <<"     Mean density of fudged SimLSS cube = "<< mean_dens <<endl;


        // Print min and max values in the grids
	    double min, max;
	    
	    cout <<"     Min and max of galaxy grid:"<<endl;
	    ngal.MinMax(min, max);
	    cout <<"     min = "<< min <<", max = "<< max <<endl;
	    
	    cout <<"     Min and max of random galaxy grid:"<<endl;
	    wrgal.MinMax(min, max);
	    cout <<"     min = "<< min <<", max = "<< max <<endl;
	 
	    cout <<"     Min and max of redshifts grid:"<<endl;
	    zc.MinMax(min, max);
	    cout <<"     min = "<< min <<", max = "<< max <<endl;
	
		cout <<"     Min and max of weighted galaxy grid:"<<endl;
		cout <<"     (same as galaxy grid if original catalog had no";
		cout <<" selection effects) "<<endl;
		wngal.MinMax(min, max);
		cout <<"     min = "<< min <<", max = "<< max <<endl;

	
	
	    res.Update();
	    cout << " Memory size (KB):" << res.getMemorySize() << endl;
	    cout << " Max memory size (KB):" << res.getMaxMemorySize() << endl;
	    cout << " Maximum allowed data segment size (KB):"<<res.getMaxDataSize()<<endl;
	    cout << " Resource usage info : \n" << res << endl;
	
	
	    // Extract the sub arrays
	    cout <<"     Extract sub arrays"<<endl;
	    TArray<r_8> nsub;
	    TArray<r_8> wsub, wrsub;
			
	    cout << "     Extract NGALS"<<endl;
	    nsub = ngal(Range(x1,x2), Range(y1,y2), Range(z1,z2)).PackElements();

	    res.Update();
	    cout << "     Extracted ngals sub array "<<endl;
	    cout << " Memory size (KB):" << res.getMemorySize() << endl;
	    cout << " Resource usage info : \n" << res << endl;
	    
	    cout << "     Check ngals sub array: "<<endl;
	    cout << nsub(0,0,0) <<"   "<< nsub(0,1,0) <<"   "<< nsub(0,2,0) << endl;
	    cout << nsub(1,0,0) <<"   "<< nsub(1,1,0) <<"   "<< nsub(1,2,0) << endl;
	    cout << nsub(2,0,0) <<"   "<< nsub(2,1,0) <<"   "<< nsub(2,2,0) << endl;
		
	    cout << "     Extract RANDOM"<<endl;
	    wrsub = wrgal(Range(x1,x2), Range(y1,y2), Range(z1,z2)).PackElements();

	    res.Update();
	    cout << "     Extracted random sub array "<<endl;
	    cout << " Memory size (KB):" << res.getMemorySize() << endl;
	    cout << " Resource usage info : \n" << res << endl;
	
	    cout << "     Check wrngals sub array: "<<endl;
	    cout << wrsub(0,0,0) <<"   "<< wrsub(0,1,0) <<"   "<< wrsub(0,2,0) << endl;
	    cout << wrsub(1,0,0) <<"   "<< wrsub(1,1,0) <<"   "<< wrsub(1,2,0) << endl;
	    cout << wrsub(2,0,0) <<"   "<< wrsub(2,1,0) <<"   "<< wrsub(2,2,0) << endl;

		cout << "     Extract WEIGHTED"<<endl;
		wsub = wngal(Range(x1,x2), Range(y1,y2), Range(z1,z2)).PackElements();

		res.Update();
		cout << "     Extracted weighted sub array "<<endl;
		cout << " Memory size (KB):" << res.getMemorySize() << endl;
		cout << " Resource usage info : \n" << res << endl;
		
		cout << "     Check wngals sub array: "<<endl;
		cout << wsub(0,0,0) <<"   "<< wsub(0,1,0) <<"   "<< wsub(0,2,0) << endl;
		cout << wsub(1,0,0) <<"   "<< wsub(1,1,0) <<"   "<< wsub(1,2,0) << endl;
		cout << wsub(2,0,0) <<"   "<< wsub(2,1,0) <<"   "<< wsub(2,2,0) << endl;	
		
		
		// Write sub arrays to a file
	    cout <<"     Writing sub-arrays to "<< outfile <<endl;
	    FitsInOutFile fos(outfile, FitsInOutFile::Fits_Create);	

	    fos << nsub;
	    fos << wrsub;
		fos << wsub;
		
	    if (DeBug) {
		
		    string outfiledb; 
		
		    outfiledb= dbfile +"_ngals.fits";
		    cout <<"    Writing ngals subarray to "<< outfiledb <<endl;
		    FitsInOutFile fosdb1(outfiledb, FitsInOutFile::Fits_Create);	
		    fosdb1 << nsub;
		
		    outfiledb= dbfile +"_wrgals.fits";
		    cout <<"    Writing wrgals subarray to "<<outfiledb<<endl;
		    FitsInOutFile fosdb2(outfiledb, FitsInOutFile::Fits_Create);	
		    fosdb2 << wrsub;
		
			outfiledb= dbfile +"_wngals.fits";
			cout <<"    Writing wngals subarray to "<<outfiledb<<endl;
			FitsInOutFile fosdb3(outfiledb, FitsInOutFile::Fits_Create);	
			fosdb3 << wsub;
		
		    }
		
			
	    // Find redshift bounds of arrays		
	    cout <<"    Find redshift bounds of arrays"<<endl;
	    sa_size_t xc,yc; // center of sub-array face in x,y coords
	    sa_size_t zl,zcen,zh;// front,center,back of subarray in z coord
	
	    zl = z1;
	    zh = z2;
	    zcen = (z2+z1)/2;
	    xc = (x2+x1)/2;
	    yc = (y2+y1)/2;
	
	    cout <<"     Index of center front pixel: ("<< xc <<","<< yc <<","<< zl;
	    cout <<")"<<endl;
	    cout <<"     Index of center pixel: ("<< xc <<","<< yc <<","<< zcen;
	    cout <<")"<<endl;
	    cout <<"     Index of center back pixel: ("<< xc <<","<< yc <<","<< zh;
	    cout <<")"<<endl;
	
	    cout <<"    Redshift of center front pixel: "<< zc(xc,yc,zl) <<endl;
	    cout <<"    Redshift of center pixel: "<< zc(xc,yc,zcen) <<endl;
	    cout <<"    Redshift of center back pixel: "<<zc(xc,yc,zh) <<endl;

	
	    // write data to FITS header
	    fos.WriteKey("MeanOverDensity",mean_dens," mean dens SimLSS delta-fudge");
	    fos.WriteKey("ZMIN",zc(xc,yc,zl),"z of front,center pixel");
	    fos.WriteKey("ZCEN",zc(xc,yc,zcen),"z of central pixel");
	    fos.WriteKey("ZMAX",zc(xc,yc,zh),"z of back,center pixel");
	    fos.WriteKey("X1",x1,"x-dim sub-array pixel start"); 
	    fos.WriteKey("X2",x2,"x-dim sub-array pixel end"); 	
	    fos.WriteKey("XC",xc,"x-dim sub-array pixel center"); 
	    fos.WriteKey("Y1",y1,"y-dim sub-array pixel start"); 
	    fos.WriteKey("Y2",y2,"y-dim sub-array pixel end"); 	
	    fos.WriteKey("YC",yc,"y-dim sub-array pixel center"); 
	    fos.WriteKey("Z1",z1,"z-dim sub-array pixel start"); 
	    fos.WriteKey("Z2",z2,"z-dim sub-array pixel end"); 	
	    fos.WriteKey("ZC",zcen,"z-dim sub-array pixel center"); 
	    fos.WriteKey("R",R,"pixel size Mpc"); 
		
	    res.Update();
	    cout << "     Written arrays to file "<<endl;
	    cout << " Memory size (KB):" << res.getMemorySize() << endl;
	    cout << " Resource usage info : \n" << res << endl;
	    }
	else { // file is a SimLSS density grid

	    FitsInOutFile fsin(infile,FitsInOutFile::Fits_RO);
	    TArray<r_8> drho; fsin >> drho; 
	    TArray<r_8> dsub;
			
	    cout << "    Extracting SimLSS sub cube"<<endl;
	    dsub = drho(Range(x1,x2), Range(y1,y2), Range(z1,z2)).PackElements();
	    FitsInOutFile fos(outfile, FitsInOutFile::Fits_Create);	
	    fos << dsub;
	    }
	

    }  // End of try bloc 
  
  
    catch (PThrowable & exc) { // catching SOPHYA exceptions
        cerr << " subfromfull.cc: Catched Exception (PThrowable)" << (string)typeid(exc).name() 
             << "\n...exc.Msg= " << exc.Msg() << endl;
        rc = 99;
	    }
    catch (std::exception & e) { // catching standard C++ exceptions
        cerr << " subfromfull.cc: Catched std::exception "  << " - what()= " << e.what() << endl;
        rc = 98;
	    }
    catch (...) { // catching other exceptions
        cerr << " subfromfull.cc: some other exception (...) was caught ! " << endl;
        rc = 97;
	    }
    cout << " ==== End of subfromfull.cc program  Rc= " << rc << endl;
    return rc;	
}
