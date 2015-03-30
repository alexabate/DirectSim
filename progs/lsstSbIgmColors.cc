#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cstdlib>

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


#include "constcosmo.h"
#include "cosmocalcs.h"
#include "simdata.h"
#include "sedfilter.h"
//#include "schechter.h"

// AA: updated old Sophya version's GenericFunc to ClassFunc1D

// @todo's below
// AA: put the below function into the igm.h, igm.cc source files
vector<string> returnFileList(string fname);
// AA: add these two functions as methods in PhotometryCalcs class
double computeColorIGM(double z, ClassFunc1D& sed, Filter& filterX, Filter& filterY, ClassFunc1D& transmission, double lmin, double lmax);
//double restFrameFluxIGM(ClassFunc1D &sed, Filter& filter, double zs, ClassFunc1D& transmission, double lmin, double lmax);
double restFrameFluxIGM(ClassFunc1D& sed, Filter& filter, double zs, ClassFunc1D& transmission, double lmin, double lmax);
// AA: then just can use the getFilterZeroPointFlux method in PhotometryCalcs instead of this
double getFilterZeroPointFlux(Filter& filter, double lmin, double lmax);


// AA: made a usage function
void usage(void);
void usage(void) {
	cout << endl<<" Usage: lsstSbIgmColors [...options...]" << endl<<endl;
	
    cout <<"  Calculates observed magnitudes for a galaxy at the redshift supplied (if not supplied  "<<endl;
	cout <<"  default is z=2.45) with and with absorption by a series of IGM lines of sight. For     "<<endl;
	cout <<"  galaxies with IGM absorption, magnitudes are calculated with and without LSST          "<<endl;
	cout <<"  photometric errors added.                                                              "<<endl;
	cout << endl;
	
	cout << " -o: OUTLOC: write files to location OUTLOC                                              "<<endl;
	cout << " -l: FILTLOC,SEDLOC,TRANSLOC: set locations of filters, SEDs and IGM transmissions       "<<endl;
	cout << " -z: ZSOURCE: redshift of source galaxy                                                  "<<endl; 
	cout << " -n: NLINES: number of IGM lines of sight to process                                     "<<endl;
	cout << endl;
    };


int main(int narg, char* arg[])
{

    /*************************************
    *  Make Sure Sophya is Initialized   *
    *************************************/

    SophyaInit();  
    FitsIOServerInit();
    InitTim();
    cout<<endl<<endl;
    
    //--- decoding command line arguments 

//for general use 
//    string outloc; 
//    string FILTLOC="./filters/";
//    string SEDLOC="./SEDs/";
//    string TRANSLOC="./transmission/";

//hard coded in the subdirectories as I'm working through this 
    string outloc= "/home/lidenscheng/MK_DirectSim/testfiles/StarburstLOS1Yr/";
//    string outloc= "/home/lidenscheng/MK_DirectSim/testfiles/StarburstLOS10Yr/";
    string FILTLOC= "/home/lidenscheng/MK_DirectSim/filters/";
//    string SEDLOC= "/home/lidenscheng/MK_DirectSim/SEDs/";
    string SEDLOC= "/home/lidenscheng/DirectSim/SEDs/"; 
    string TRANSLOC= "/home/lidenscheng/MK_DirectSim/transmission/";
//    string TRANSLOC= "/home/lidenscheng/MK_DirectSim/meanIGMTransmissions/";

    string locs;
    int nz = 400;    //number of lines of sight 
    double z = 1.8; //constant redshift 

	char c;
    while((c = getopt(narg,arg,"ho:l:z:n:")) != -1) {
	    switch (c)  {
	        case 'o' :
	            outloc = optarg;
	            break;
	        case 'l' : {
	            locs = optarg;
	            string delim=",";
	            vector<string> results;
	            stringSplit(locs, delim, results);
	            FILTLOC=results[0];
	            SEDLOC=results[1];
	            TRANSLOC=results[2];
	            }
	            break;
	        case 'z' :
	            sscanf(optarg,"%lf",&z);
	            break;
	        case 'n' :
	            sscanf(optarg,"%d",&nz);
	            break;
	        case 'h' :
		        default :
		    usage(); return -1;
		    }
	    }
	// AA: @note could add option of SED list to read in (set sedFile)
	// AA: @note could add option of number of years of LSST observation (set nYear)
	// AA: @note could add option to set different i band mag (set imag_fixed)
	    
    //-- end command line arguments
    cout <<"     Writing data to files located at "<< outloc <<endl;
    cout <<"     Filters located at "<< FILTLOC << endl;
    cout <<"     SEDs located at "<< SEDLOC << endl;
    cout <<"     Transmissions located at "<< TRANSLOC <<endl;
    cout <<"     Redshift of source galaxies "<< z <<endl;
    cout <<"     Calculating using "<< nz <<" IGM lines of sight"<<endl;
    cout << endl;
    //-- end command line arguments
  
  int rc = 1;  
  try {  // exception handling try bloc at top level
    InitTim();


    /*************************************
    *       Set Parameters               *
    *************************************/

    // wavelength range of the SEDs/filters
    //double lmin=3000.e-10, lmax=11000.e-10;
    double lmin = 5e-8, lmax = 2.5e-6;
	
	
    // Set environment variables
    // AA: removed hard coding so others can use
    string tmp1,tmp2,tmp3;
    char *c1,*c2,*c3;
    
    tmp1 = "FILTLOC=" + FILTLOC; 
    c1 = &tmp1[0];
    putenv(c1);
    
    tmp2 = "SEDLOC=" + SEDLOC; 
    c2 = &tmp2[0];
    putenv(c2);
    
    tmp3 = "TRANSLOC=" + TRANSLOC; 
    c3 = &tmp3[0];
    putenv(c3);
	
    // Set output files
    // AA: removed hard coding, made vector of strings to prevent code copy
    // AA: can further improve this by reading first part of filename from SED list read in
    stringstream ss;
    ss << z;
    vector<string> outfiles;

    outfiles.push_back(outloc + "SB3_B2004a_" + ss.str() + "z_Mags.txt");
    outfiles.push_back(outloc + "SB2_B2004a_" + ss.str() + "z_Mags.txt");
    outfiles.push_back(outloc + "ssp_25Myr_z008_" + ss.str() + "z_Mags.txt");
    outfiles.push_back(outloc + "ssp_5Myr_z008_" + ss.str() + "z_Mags.txt");

//meanIGM files use constant IGM (IGM averaged over 400LOS)
//    outfiles.push_back(outloc + "meanIGM_SB3_B2004a_1%_" + ss.str() + "z_Mags.txt");
//    outfiles.push_back(outloc + "meanIGM_SB2_B2004a_1%_" + ss.str() + "z_Mags.txt");
//    outfiles.push_back(outloc + "meanIGM_ssp_25Myr_z008_1%_" + ss.str() + "z_Mags.txt");
//    outfiles.push_back(outloc + "meanIGM_ssp_5Myr_z008_1%_" + ss.str() + "z_Mags.txt");


    /*************************************
    *       Read in Filters              *
    *************************************/


    string filterFile = "LSST.filters";
    ReadFilterList lsstFilters(filterFile);
    lsstFilters.readFilters(lmin,lmax);
    vector<Filter*> filterArray = lsstFilters.getFilterArray();
    int nFilters = lsstFilters.getNTot();
    cout <<"     "<< nFilters <<" LSST filters read in "<<endl;

    // LSST filter indices
    int uLSST = 0;
    int gLSST = 1;
    int rLSST = 2;
    int iLSST = 3;
    int zLSST = 4;
    int yLSST = 5;


    /*************************************
    *           Read in SEDs             *
    *************************************/

//    string sedFile = "LSSTBurst1.list";//"LSST.list";
    string sedFile = "igmSBGalaxies.list";

    ReadSedList cwwSEDs(sedFile);
    // Read out SEDs into array
    cwwSEDs.readSeds(lmin,lmax);
    // Get total number of SEDs
    vector<SED*> sedArray=cwwSEDs.getSedArray();
    int nsed=cwwSEDs.getNSed();
    cout <<"     Number of original SEDs = "<<nsed<<endl;
    cout << endl;
    
    
    
    /*****************************************************************
    *     Read in mean transmission files or LOS transmission files  *
    ******************************************************************/

    string transFile = "Transmission_" + ss.str() + "z.list";
//    string transFile = "meanTransmission_" + ss.str() + "z.list";
    vector<string> transFileList = returnFileList(transFile);
    vector<SInterp1D> transIgm;

    for(int i=0; i<transFileList.size()-1; i++) {
        SInterp1D trans;
        trans.ReadXYFromFile(transFileList[i], lmin, lmax, 1024, 0, false);
        transIgm.push_back(trans);
        }


// AA: vv THIS SECTION CALCS NOT ACTUALLY USED
    /*************************************
    *     Calculate Flux (No IGM)        *
    *************************************/
    int iSb = 0;

    PhotometryCalcs photometryCalcs(lmin,lmax);

//for(iSb=0; iSb<nSb; iSb)
//{
	
    TVector<r_8> uFlux(nz), gFlux(nz), rFlux(nz), iFlux(nz), zFlux(nz), yFlux(nz);
    for (int i=0; i<nz; i++) {
	
//        double z=zmin+i*dz;

        // flux in u filter
        uFlux(i)=photometryCalcs.restFrameFlux((*sedArray[iSb]),
                                (*filterArray[uLSST]), z);
        // flux in g filter
        gFlux(i)=photometryCalcs.restFrameFlux((*sedArray[iSb]),
                                (*filterArray[gLSST]), z);
        // flux in r filter
        rFlux(i)=photometryCalcs.restFrameFlux((*sedArray[iSb]),
                                (*filterArray[rLSST]), z);
        // flux in i filter
        iFlux(i)=photometryCalcs.restFrameFlux((*sedArray[iSb]),
                                (*filterArray[iLSST]), z);
        // flux in z filter
        zFlux(i)=photometryCalcs.restFrameFlux((*sedArray[iSb]),
                                (*filterArray[zLSST]), z);
        // flux in y filter
        yFlux(i)=photometryCalcs.restFrameFlux((*sedArray[iSb]),
                                (*filterArray[yLSST]), z);
    }



    /*************************************
    *      Calculate Mag (No IGM)        *
    *************************************/

    TVector<r_8> uMag(nz), gMag(nz), rMag(nz), iMag(nz), zMag(nz), yMag(nz);
    for (int i=0; i<nz; i++) {

//      double z = zmin + i*dz;

        // magnitude in u filter
        uMag(i) = photometryCalcs.convertFluxMaggiesToABMag(uFlux(i), 
                                (*filterArray[uLSST]));
        // magnitude in g filter
        gMag(i) = photometryCalcs.convertFluxMaggiesToABMag(gFlux(i), 
                                (*filterArray[gLSST]));
        // magnitude in r filter
        rMag(i) = photometryCalcs.convertFluxMaggiesToABMag(rFlux(i), 
                                (*filterArray[rLSST]));
        // magnitude in i filter
        iMag(i) = photometryCalcs.convertFluxMaggiesToABMag(iFlux(i), 
                                (*filterArray[iLSST]));
        // magnitude in z filter
        zMag(i) = photometryCalcs.convertFluxMaggiesToABMag(zFlux(i), 
                                (*filterArray[zLSST]));
        // magnitude in y filter
        yMag(i) = photometryCalcs.convertFluxMaggiesToABMag(yFlux(i), 
                                (*filterArray[yLSST]));
    }







    /************************************
    *  Calculate Flux - IGM Absorption  *
    ************************************/

    TVector<r_8> uFluxR(nz), gFluxR(nz), rFluxR(nz), iFluxR(nz), zFluxR(nz), yFluxR(nz);
    for (int i=0; i<nz; i++) {

	
//        double z=zmin+i*dz;


        // flux in u filter
        uFluxR(i)=restFrameFluxIGM((*sedArray[iSb]), 
                                (*filterArray[uLSST]), z,
                                transIgm[i], lmin, lmax);
        // flux in g filter
        gFluxR(i)=restFrameFluxIGM((*sedArray[iSb]),
                                (*filterArray[gLSST]), z,
                                transIgm[i], lmin, lmax);
        // flux in r filter
        rFluxR(i)=restFrameFluxIGM((*sedArray[iSb]),
                                (*filterArray[rLSST]), z,
                                transIgm[i], lmin, lmax);
        // flux in i filter
        iFluxR(i)=restFrameFluxIGM((*sedArray[iSb]),
                                (*filterArray[iLSST]), z,
                                transIgm[i], lmin, lmax);
        // flux in z filter
        zFluxR(i)=restFrameFluxIGM((*sedArray[iSb]),
                                (*filterArray[zLSST]), z,
                                transIgm[i], lmin, lmax);
        // flux in y filter
        yFluxR(i)=restFrameFluxIGM((*sedArray[iSb]),
                                (*filterArray[yLSST]), z,
                                transIgm[i], lmin, lmax);
    }


    /*************************************
    *      Calculate Mag (IGM Absorption)        *
    *************************************/

    TVector<r_8> uMagR(nz), gMagR(nz), rMagR(nz), iMagR(nz), zMagR(nz), yMagR(nz);
    for (int i=0; i<nz; i++) {

//        double z = zmin + i*dz;

        // magnitude in u filter
        uMagR(i) = photometryCalcs.convertFluxMaggiesToABMag(uFluxR(i), 
                                (*filterArray[uLSST]));
        // magnitude in g filter
        gMagR(i) = photometryCalcs.convertFluxMaggiesToABMag(gFluxR(i), 
                                (*filterArray[gLSST]));
        // magnitude in r filter
        rMagR(i) = photometryCalcs.convertFluxMaggiesToABMag(rFluxR(i), 
                                (*filterArray[rLSST]));
        // magnitude in i filter
        iMagR(i) = photometryCalcs.convertFluxMaggiesToABMag(iFluxR(i), 
                                (*filterArray[iLSST]));
        // magnitude in z filter
        zMagR(i) = photometryCalcs.convertFluxMaggiesToABMag(zFluxR(i), 
                                (*filterArray[zLSST]));
        // magnitude in y filter
        yMagR(i) = photometryCalcs.convertFluxMaggiesToABMag(yFluxR(i), 
                                (*filterArray[yLSST]));
    }
// AA: ^^THIS SECTION CALCS NOT ACTUALLY USED ----- END


    /*************************************
    *        Write mags to a file        *
    *************************************/
	
    double m_u, m_g, m_r, m_i, m_z, m_y; 
    double m_uR, m_gR, m_rR, m_iR, m_zR, m_yR; 
    vector<double> error_U, error_G, error_R, error_I, error_Z, error_Y; //photometric errors  
    TVector<r_8> Cug(nz), Cgr(nz), Cri(nz), Ciz(nz), Czy(nz), Cyu(nz);
    TVector<r_8> CugR(nz), CgrR(nz), CriR(nz), CizR(nz), CzyR(nz), CyuR(nz);

// Either 1Yr or 10Yr LSST photometric errors
//    int nYear = 10; // could as as prog arg
    int nYear = 1;
    
    // AA: No longer needs these numbers (they are coded within SimData class)
    // Number of visits per year (Table 1, Ivezic et al 2008)
    //int uVisitsPerYear = 6;
    //int gVisitsPerYear = 8;
    //int rVisitsPerYear = 18;
    //int iVisitsPerYear = 18;
    //int zVisitsPerYear = 16;
    //int yVisitsPerYear = 16;

    // total number of visits
    //int uVisits = uVisitsPerYear*nYear;
    //int gVisits = gVisitsPerYear*nYear;
    //int rVisits = rVisitsPerYear*nYear;
    //int iVisits = iVisitsPerYear*nYear;
    //int zVisits = zVisitsPerYear*nYear;
    //int yVisits = yVisitsPerYear*nYear;

    SimpleUniverse su;
    RandomGenerator rg;
    
    // AA: code m_i=24. so we know what it is
    double imag_fixed = 24.; // can make this prog arg

    // Initialize class that simulates magnitude errors
    SimData simData(sedArray, filterArray, su, rg);


    // AA: add loop over SEDs to save code copying
    for (int iSb=0; iSb<sedArray.size(); iSb++) {

    
        ifstream inp;
        ofstream outp;
        inp.open(outfiles[iSb].c_str(), ifstream::in);
        inp.close();
        if(inp.fail()) {
            inp.clear(ios::failbit);
            cout << "Writing to file ..." << outfiles[iSb].c_str() << endl;
            outp.open(outfiles[iSb].c_str(), ofstream::out);

            // Faking magnitudes (no IGM and IGM) by using colors. Manually set m_i=24. 


		
	        for (int i=0; i<nz; i++) {
		

//            double z=zmin+i*dz;

		        // mag in filter u - mag in filter g 
		        Cug(i) = photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[uLSST]), (*filterArray[gLSST]));
		
		        // mag in filter g - mag in filter r 
		        Cgr(i) = photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[gLSST]), (*filterArray[rLSST]));
		
		        // mag in filter r - mag in filter i 
		        Cri(i)=photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[rLSST]), (*filterArray[iLSST]));
								
		        // mag in filter i - mag in filter z 
		        Ciz(i) = photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[iLSST]), (*filterArray[zLSST]));
		
		        // mag in filter z - mag in filter y 
		        Czy(i) = photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[zLSST]), (*filterArray[yLSST]));

		        // mag in filter y - mag in filter z 
		        Cyu(i) = photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[yLSST]), (*filterArray[uLSST]));
								
			    // no IGM mags
	            m_r = Cri(i) + imag_fixed;
	            m_i = imag_fixed;
	            m_z = imag_fixed - Ciz(i);
	            m_y = m_z - Czy(i);
	            m_g = Cgr(i) + m_r;
	            m_u = Cug(i)+m_g; 


                // Colors with IGM absorption 
               CugR(i) = computeColorIGM(z, (*sedArray[iSb]), (*filterArray[uLSST]), (*filterArray[gLSST]), transIgm[i], lmin, lmax); 
               CgrR(i) = computeColorIGM(z, (*sedArray[iSb]), (*filterArray[gLSST]), (*filterArray[rLSST]), transIgm[i], lmin, lmax); 
               CriR(i) = computeColorIGM(z, (*sedArray[iSb]), (*filterArray[rLSST]), (*filterArray[iLSST]), transIgm[i], lmin, lmax); 
               CizR(i) = computeColorIGM(z, (*sedArray[iSb]), (*filterArray[iLSST]), (*filterArray[zLSST]), transIgm[i], lmin, lmax); 
               CzyR(i) = computeColorIGM(z, (*sedArray[iSb]), (*filterArray[zLSST]), (*filterArray[yLSST]), transIgm[i], lmin, lmax); 

               // w/IGM mags
	           m_rR = CriR(i) + imag_fixed;
	           m_iR = imag_fixed;
	           m_zR = imag_fixed - CizR(i);
	           m_yR = m_zR - CzyR(i);
	           m_gR = CgrR(i) + m_rR;
	           m_uR = CugR(i) + m_gR;

//	for(int k=0; k<400; k++){

              // AA: updated to new method that adds LSST errors
              // AA: @ note photometric errors only added to IGM absorbed magnitudes
	          error_U = simData.addLSSTError(m_uR, nYear, uLSST); 
	          error_G = simData.addLSSTError(m_gR, nYear, gLSST); 
	          error_R = simData.addLSSTError(m_rR, nYear, rLSST); 
	          error_I = simData.addLSSTError(m_iR, nYear, iLSST); 
	          error_Z = simData.addLSSTError(m_zR, nYear, zLSST); 
	          error_Y = simData.addLSSTError(m_yR, nYear, yLSST); 

//           	outp << z 
//                << "    " << m_u << "    " << m_g<< "    " << m_r << "    " 
//		 << 24.0 << "    " << m_z << "    " << m_y << "	" << 0.0 << endl;

            // write to the file: igm line of sight index, mags_noerror_noigm, mags_noerror_wigm, mags_werror_wigm, errors
            // (redshift is in the filename)

           	outp << i+1 <<"  "; // IGM line of sight index
           	outp << m_u <<"  "<< m_g <<"  "<< m_r <<"  "<< m_i << "  "<< m_z <<"  "<< m_y <<"  "; // no error no IGM (same each row) 
           	outp << m_uR <<"  "<< m_gR <<"  "<< m_rR <<"  "<< m_iR <<"  "<< m_zR <<"  "<< m_yR <<"  "<< 0.0 << " "; // no error w/IGM
           	outp << error_U[0] <<"  "<< error_G[0] <<"  "<< error_R[0] <<"  "<< error_I[0] <<"  "<< error_Z[0] <<"  "<< error_Y[0] <<"  "; // w/error w/IGM
           	outp << error_U[1] <<"  "<< error_G[1] <<"  "<< error_R[1] <<"  "<< error_I[1] <<"  "<< error_Z[1] <<"  "<< error_Y[1] <<"  "; // errors

                outp << endl;

//           	outp2 << i+1 
//                 << "    " << error_U[1] << "    " << error_G[1] << "    " << error_R[1] << "    " 
//		 << error_I[1] << "    " << error_Z[1] << "    " << error_Y[1] << endl;

//only calculate and print u-g, g-r colors for 400 LOS at some particular z 
//              outp << i+1 << "  " << Cug(i) << "	" << Cgr(i) << endl;
//              outp << i+1 << "	" << CugR(i) << "	" << CgrR(i) << endl; 
  
//only calculate and print u-g, g-r colors for different z with no IGM and with IGM 
//              outp << z << "  " << Cug(i) << "	" << Cgr(i) << endl; 
//              outp << z << "  " << CugR(i) << " 	" << CgrR(i) << endl;  

//            }
         }

        outp.close();
        //outp2.close();
	    }

    else {
        cout << "Error...file " << outfiles[iSb].c_str() << " exists" << endl;
        //cout << "Error...file " << outfile2.c_str() << " exists" << endl;
	    }
    }
 
  }  // End of try bloc 

  
  
  catch (PThrowable & exc) {  // catching SOPHYA exceptions
    cerr << " lsstSbIgmColors.cc: Catched Exception (PThrowable)" 
     << (string)typeid(exc).name() 
         << "\n...exc.Msg= " << exc.Msg() << endl;
    rc = 99;
  }
  catch (std::exception & e) {  // catching standard C++ exceptions
    cerr << " lsstSbIgmColors.cc: Catched std::exception "  << " - what()= " 
     << e.what() << endl;
    rc = 98;
  }
  catch (...) {  // catching other exceptions
    cerr << " lsstSbIgmColors.cc: some other exception (...) was caught ! " 
     << endl;
    rc = 97;
  }
  cout << " ==== End of lsstSbIgmColors.cc program  Rc= " << rc << endl;
  return rc;	
}


/*****************************************************
* Opens a file of filenames and makes it into a list *
*****************************************************/

vector<string> returnFileList(string fname) {
    vector<string> flist;

    char * pt = getenv("TRANSLOC");
    string fpath = "";
    if(pt==NULL)
        throw ParmError("ERROR TRANSMISSON LOCATION ENVIRONMENT VARIABLE -TRANSLOC- NOT DEFINED");
    else
        fpath = pt+fname;

    ifstream ifile;
    ifile.open(fpath.c_str());

    string tmp;
    while(!ifile.eof()) {
        getline(ifile, tmp);
        flist.push_back(tmp);
    }

    ifile.close();

    return flist;
}


/**************************************************
*   Computer Galaxy Color with IGM Transmission   *
**************************************************/

double computeColorIGM(double z, ClassFunc1D& sed, Filter& filterX, Filter& filterY, ClassFunc1D& transmission, double lmin, double lmax) {

    SEDzFilterProdIGM SEDFX(sed, filterX, transmission, z);
    SEDzFilterProdIGM SEDFY(sed, filterY, transmission, z);

    // Integrate SED across filter
    FilterIntegrator intsedX(SEDFX, lmin, lmax);
    FilterIntegrator intsedY(SEDFY, lmin, lmax);

    // Calculate filter zero points
    FilterProd FX(filterX);
    FilterProd FY(filterY);
    FilterIntegrator intFX(FX, lmin, lmax);
    FilterIntegrator intFY(FY, lmin, lmax);
    double zpCxy = -2.5*log10(intFY.Value()/intFX.Value());

    // Calculate color
    double Cxy = -2.5*log10(intsedX.Value()/intsedY.Value()) + zpCxy;

    return Cxy;
}

/***************************************************
*   Compute Galaxy Magntude with IGM Transmission  *
***************************************************/

//double restFrameFluxIGM(ClassFunc1D &sed, Filter& filter, double zs, ClassFunc1D& transmission, double lmin, double lmax) {

double restFrameFluxIGM(ClassFunc1D& sed, Filter& filter, double zs, ClassFunc1D& transmission, double lmin, double lmax) {

    SEDzFilterProdIGM sedXfilter(sed, filter, transmission, zs);
    FilterIntegrator intSED(sedXfilter, lmin, lmax);

    double zpFluxFilter = getFilterZeroPointFlux(filter, lmin, lmax);
    double f0 = intSED.Value()/zpFluxFilter;

    return f0;
}

/********************************
*   Compute Filter Zero Point   *
********************************/

double getFilterZeroPointFlux(Filter& filter, double lmin, double lmax) {
    FilterProd FX(filter);
    FilterIntegrator intFX(FX, lmin, lmax);
    double zp = intFX.Value();

  
  return zp;
}


