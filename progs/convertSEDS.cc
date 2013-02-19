// -*- LSST-C++ -*-
#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <typeinfo>

// sophya libraries
//#include "histinit.h"
#include "fiosinit.h"
//#include "mydefrg.h"
#include "tarray.h"
#include <pexceptions.h>


#define PI 3.141592
/*



*/
void usage(void);
void usage(void)
{
	cout << endl<<" Usage: convertSEDs [...options...]" << endl<<endl;
	cout << " -i [SEDLIST]: File containing list of SEDs to read in, must be in same directory as SEDs.\n";
	cout << "               Must have full path to file (unless in current dir) [DEFAULT=LSST.list]"<<endl;
	cout << " -d [DIRNAME]: Directory to save to [DEFAULT=SEDs/]"<<endl;
	cout << " -u [UNITS]: Units the SEDs are in, either NM (nanometers) or AA (angstroms) [DEFAULT=NM]"<<endl;	
	cout << endl;
}


string splitFilename (const string& str)
{
  size_t found;
  //cout << "Splitting: " << str << endl;
  found=str.find_last_of("/\\");
  string directory=str.substr(0,found);
  return directory;
  //cout << " folder: " << str.substr(0,found) << endl;
  //cout << " file: " << str.substr(found+1) << endl;
}

int main(int narg, char* arg[])
{
    cout << " ==== convertSEDs.cc program ==== "<<endl;

    // make sure SOPHYA modules are initialized 
    SophyaInit();  
    FitsIOServerInit();
    cout<<endl<<endl;


    //--- decoding command line arguments 
    string inFile=" /raid00/catalogs/LSST/SEDS/LSST.list",dirOut="SEDs/";
    string units="NM";
    string nanoMeters="NM",angstroms="AA";
    bool isNano=true;

    cout << " ==== decoding command line arguments ===="<<endl;
    char c;
    while((c = getopt(narg,arg,"hi:d:u:")) != -1) {
        switch (c) {
            case 'i' :
                inFile = optarg;
                break;
            case 'd' :
                dirOut = optarg;
                break;
            case 'u' :
                units = optarg;
                break;
            case 'h' :
                default :
                usage(); return -1;
	        }
        }

    if (units==angstroms)
        isNano = false;
    else if (units==nanoMeters)
        isNano = true;
    else
        throw ParmError("ERROR! Don't understand spectra units");

    //-- end command line arguments
    cout <<"     Reading in files listed in "<<inFile<<endl;
    cout <<"     Converting from";
    if (isNano)
        cout <<" nanometers";
    else
        cout <<" angstroms";
    cout <<" to meters"<<endl;
    cout <<"     Saving to "<<dirOut<<" directory"<<endl;
    cout << endl;

    int rc = 1;  
    try {  // exception handling try bloc at top level

    ifstream ifs;
    ofstream outp;
    string line;

    // Unit convertor
    double unitConversion;
    if (isNano)
        unitConversion=1e-9;
    else
        unitConversion=1e-10;
        

    // Separate list filename and path to list file
    string dirIn=splitFilename(inFile);
    cout <<"     SEDs are in directory "<<dirIn<<endl;

    // Count number of SEDS in the file
    int nsed=0;
    ifs.open(inFile.c_str(),ifstream::in);
    if (ifs.fail()){
	    string emsg="error: Unable to open file ";
	    emsg+=inFile;
	    throw ParmError(emsg);
	    }

    while ( ifs.good() ){
	    getline(ifs,line);
        nsed++;
	    }
    nsed-=1;
    ifs.close();
    cout <<"     Number of SEDs to read in "<<nsed<<endl;

    // read in all the SED filenames
    int prt=1;
    ifs.open(inFile.c_str(),ifstream::in);
    string fileNames[nsed];
    if (prt > 0)
        cout <<"     File contains the following SEDs:"<<endl;
    for (int i=0; i<nsed; i++)
	    {
	    getline(ifs,line);
	    fileNames[i]=line;//dir +"/"+line;
	    if (prt > 0)
	        cout <<"     "<<fileNames[i]<<endl;
	    }
    ifs.close();
    cout <<endl;
	
    // Loop over all files
    for (int i=0; i<nsed; i++) {
	
        // read in SED files
        string file=dirIn+"/"+fileNames[i];
        cout <<"     Reading in file "<<file;
        ifs.open(file.c_str(), ifstream::in);
        sa_size_t nr, nc;
        TArray<r_4> spectrum;
        spectrum.ReadASCII(ifs,nr,nc);
        
        // write out a new
        string outFile=dirOut+fileNames[i];
        cout <<", writing out file "<<outFile<<endl;
        ifs.open(outFile.c_str(),ifstream::in);
        ifs.close();
        if(ifs.fail())
	        {
	        ifs.clear(ios::failbit);
	        outp.open(outFile.c_str(),ofstream::out);
	        for (int j=0; j<nr; j++)
	            outp <<spectrum(0,j)*unitConversion<<"  "<<spectrum(1,j)<<endl;   
		
	        outp.close();
	        }
        else
	        cout <<"     ERROR! file "<<outFile<<" exists"<<endl;
        cout << endl;
        
        }


    }  // End of try bloc 


    catch (PThrowable & exc) {  // catching SOPHYA exceptions
        cerr << " convertSEDs.cc: Catched Exception (PThrowable)" 
        << (string)typeid(exc).name() << "\n...exc.Msg= " << exc.Msg() << endl;
        rc = 99;
        }
    catch (std::exception & e) {  // catching standard C++ exceptions
        cerr << " convertSEDs.cc: Catched std::exception "  << " - what()= " 
        << e.what() << endl;
        rc = 98;
        }
    catch (...) {  // catching other exceptions
        cerr << " convertSEDs.cc: some other exception (...) was caught ! " 
        << endl;
        rc = 97;
        }
    cout << " ==== End of convertSEDs.cc program  Rc= " << rc << endl;
        return rc;	
    }

