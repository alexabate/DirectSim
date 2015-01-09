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
#include "fiosinit.h"
#include "mydefrg.h"


#include "constcosmo.h"
#include "cosmocalcs.h"
#include "simdata.h"
#include "sedfilter.h"

vector<string> returnFileList(string fname);
double computeColorIGM(double z, GenericFunc& sed, Filter& filterX, Filter& filterY, GenericFunc& transmission, double lmin, double lmax);
//double restFrameFluxIGM(GenericFunc &sed, Filter& filter, double zs, GenericFunc& transmission, double lmin, double lmax);
double restFrameFluxIGM(GenericFunc& sed, Filter& filter, double zs, GenericFunc& transmission, double lmin, double lmax);
double getFilterZeroPointFlux(Filter& filter, double lmin, double lmax);


int main(int narg, char* arg[])
{

    /*************************************
    *  Make Sure Sophya is Initialized   *
    *************************************/

    SophyaInit();  
    FitsIOServerInit();
    InitTim();
    cout<<endl<<endl;
  
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
    putenv("FILTLOC=/home/lidenscheng/MK_DirectSim/filters/");

    putenv("SEDLOC=/home/lidenscheng/MK_DirectSim/SEDs/");
//    putenv("SEDLOC=/home/lidenscheng/bpz-1.99.3/SED/LSST/");

    putenv("TRANSLOC=/home/lidenscheng/MK_DirectSim/transmission/");
//    putenv("TRANSLOC=/home/lidenscheng/MK_DirectSim/meanIGMTransmissions/");
	
    // Set output files

    string outfile = "/home/lidenscheng/MK_DirectSim/testfiles/StarburstLOS10Yr/SB3_B2004a_z008_2.45z_Mags.txt";
    string outfile3 = "/home/lidenscheng/MK_DirectSim/testfiles/StarburstLOS10Yr/SB2_B2004a_z008_2.45z_Mags.txt";
    string outfile5 = "/home/lidenscheng/MK_DirectSim/testfiles/StarburstLOS10Yr/ssp_25Myr_z008_2.45z_Mags.txt";
    string outfile7 = "/home/lidenscheng/MK_DirectSim/testfiles/StarburstLOS10Yr/ssp_5Myr_z008_2.45z_Mags.txt";

    string outfile2= "/home/lidenscheng/MK_DirectSim/testfiles/StarburstLOS10Yr/SB3_B2004a_z008_2.45z_Errors.txt";
    string outfile4= "/home/lidenscheng/MK_DirectSim/testfiles/StarburstLOS10Yr/SB2_B2004a_z008_2.45z_Errors.txt";
    string outfile6= "/home/lidenscheng/MK_DirectSim/testfiles/StarburstLOS10Yr/ssp_25Myr_z008_2.45z_Errors.txt";
    string outfile8= "/home/lidenscheng/MK_DirectSim/testfiles/StarburstLOS10Yr/ssp_5Myr_z008_2.45z_Errors.txt";

//    string redOutfile = "/home/lidenscheng/MK_DirectSim/testfiles/lsstColorsExtinction.txt";

    // Redshifts
/*    double zmin=1.4, zmax=2.9;
    int nz=16;
    double dz=(zmax-zmin)/(nz-1);
*/

      int nz=400; //number of lines of sight 
      double z=2.45; //constant redshift 

    // Which SED to examine (0-3 for SB3_B2004a, SB2_B2004a, ssp_25Myr_z008, ssp_5Myr_z008)

    int iSb = 0;
//      int iSb; 
//      int nSb= 6; 

    /*************************************
    *       Read in Filters              *
    *************************************/


    string filterFile = "LSST.filters";
    ReadFilterList lsstFilters(filterFile);
    lsstFilters.readFilters(lmin,lmax);
    vector<Filter*> filterArray=lsstFilters.getFilterArray();
    int nFilters=lsstFilters.getNTot();
    cout <<"     "<<nFilters<<" LSST filters read in "<<endl;

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


    /*************************************
    *     Calculate Flux (No IGM)        *
    *************************************/

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



    /*************************************
    *     Read in mean transmission files or LOS transmission files  *
    *************************************/

    string transFile = "Transmission_2.45z.list";
//    string transFile = "Transmission.list";
    vector<string> transFileList= returnFileList(transFile);
    vector<SInterp1D> transIgm;

    for(int i=0; i<transFileList.size()-1; i++) {
        SInterp1D trans;
        trans.ReadXYFromFile(transFileList[i], lmin, lmax, 1024, 0, false);
        transIgm.push_back(trans);
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


    /*************************************
    *        Write mags to a file        *
    *************************************/
	
    double m_u, m_g, m_r, m_i, m_z, m_y; 
    double m_uR, m_gR, m_rR, m_iR, m_zR, m_yR; 
    vector<double> error_U, error_G, error_R, error_I, error_Z, error_Y; //photometric errors  
    TVector<r_8> Cug(nz), Cgr(nz), Cri(nz), Ciz(nz), Czy(nz), Cyu(nz);
    TVector<r_8> CugR(nz), CgrR(nz), CriR(nz), CizR(nz), CzyR(nz), CyuR(nz);

    int nYear = 10;
//    int nYear = 1;
    
    // Number of visits per year (Table 1, Ivezic et al 2008)
    int uVisitsPerYear = 6;
    int gVisitsPerYear = 8;
    int rVisitsPerYear = 18;
    int iVisitsPerYear = 18;
    int zVisitsPerYear = 16;
    int yVisitsPerYear = 16;

    // total number of visits
    int uVisits = uVisitsPerYear*nYear;
    int gVisits = gVisitsPerYear*nYear;
    int rVisits = rVisitsPerYear*nYear;
    int iVisits = iVisitsPerYear*nYear;
    int zVisits = zVisitsPerYear*nYear;
    int yVisits = yVisitsPerYear*nYear;

    SimpleUniverse su;
    RandomGenerator rg;

    // Initialize class that simulates magnitude errors
    SimData simData(sedArray, filterArray, su, rg);

    iSb=0; 

    ifstream inp, inp2;
    ofstream outp, outp2;
    inp.open(outfile.c_str(), ifstream::in);
    inp2.open(outfile2.c_str(), ifstream::in);
    inp.close();
    inp2.close();
    if(inp.fail() && inp2.fail()) {

        inp.clear(ios::failbit);
        inp2.clear(ios::failbit);
        cout << "Writing to file ..." << outfile.c_str() << endl;
        cout << "Writing to file ..." << outfile2.c_str() << endl;
        outp.open(outfile.c_str(), ofstream::out);
        outp2.open(outfile2.c_str(), ofstream::out);

//Faking magnitudes (no IGM and IGM) by using colors. Manually set m_i=24. 
		
	    for (int i=0;i<nz;i++) {		

//            double z=zmin+i*dz;

		// mag in filter u - mag in filter g 
		Cug(i)=photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[uLSST]),(*filterArray[gLSST]));
		
		// mag in filter g - mag in filter r 
		Cgr(i)=photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[gLSST]),(*filterArray[rLSST]));
		// mag in filter r - mag in filter i 
		Cri(i)=photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[rLSST]),(*filterArray[iLSST]));
		// mag in filter i - mag in filter z 
		Ciz(i)=photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[iLSST]),(*filterArray[zLSST]));
		// mag in filter z - mag in filter y 
		Czy(i)=photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[zLSST]),(*filterArray[yLSST]));

		// mag in filter y - mag in filter z 
		Cyu(i)=photometryCalcs.CompColor(z,(*sedArray[iSb]),
								(*filterArray[yLSST]),(*filterArray[uLSST]));
	m_r= Cri(i)+24.0;
	m_i= 24.0;
	m_z= 24.0-Ciz(i);
	m_y= m_z-Czy(i);
	m_g= Cgr(i)+m_r;
	m_u= Cug(i)+m_g; 


        // Colors with IGM absorption 
        CugR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[uLSST]), (*filterArray[gLSST]), transIgm[i], lmin, lmax); 
        CgrR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[gLSST]), (*filterArray[rLSST]), transIgm[i], lmin, lmax); 
        CriR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[rLSST]), (*filterArray[iLSST]), transIgm[i], lmin, lmax); 
        CizR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[iLSST]), (*filterArray[zLSST]), transIgm[i], lmin, lmax); 
        CzyR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[zLSST]), (*filterArray[yLSST]), transIgm[i], lmin, lmax); 

	m_rR= CriR(i)+24.0;
	m_iR= 24.0;
	m_zR= 24.0-CizR(i);
	m_yR= m_zR-CzyR(i);
	m_gR= CgrR(i)+m_rR;
	m_uR= CugR(i)+m_gR;

	error_U= simData.addLSSTuError(m_uR, uVisits); 
	error_G= simData.addLSSTgError(m_gR, gVisits); 
	error_R= simData.addLSSTrError(m_rR, rVisits); 
	error_I= simData.addLSSTiError(m_iR, iVisits); 
	error_Z= simData.addLSSTzError(m_zR, zVisits); 
	error_Y= simData.addLSSTyError(m_yR, yVisits); 

//           	outp << z 
//                << "    " << m_u << "    " << m_g<< "    " << m_r << "    " 
//		 << 24.0 << "    " << m_z << "    " << m_y << "	" << 0.0 << endl;

           	outp << i+1 
                 << "    " << m_uR << "    " << m_gR << "    " << m_rR << "    " 
		 << 24.0 << "    " << m_zR << "    " << m_yR << "	" << 0.0 << endl;

           	outp2 << i+1 
                 << "    " << error_U[1] << "    " << error_G[1] << "    " << error_R[1] << "    " 
		 << error_I[1] << "    " << error_Z[1] << "    " << error_Y[1] << endl;

//only calculate and print u-g, g-r colors for 400 LOS at some particular z 
//              outp << i+1 << "  " << Cug(i) << "	" << Cgr(i) << endl;
//              outp << i+1 << "	" << CugR(i) << "	" << CgrR(i) << endl; 
  
//only calculate and print u-g, g-r colors for different z with no IGM and with IGM 
//              outp << z << "  " << Cug(i) << "	" << Cgr(i) << endl; 
//              outp << z << "  " << CugR(i) << " 	" << CgrR(i) << endl;  

}


        outp.close();
        outp2.close();
	    }

    else{
        cout << "Error...file " << outfile.c_str() << " exists" << endl;
        cout << "Error...file " << outfile2.c_str() << " exists" << endl;
	}



    iSb=1; 

    ifstream inp3, inp4;
    ofstream outp3, outp4;
    inp3.open(outfile3.c_str(), ifstream::in);
    inp4.open(outfile4.c_str(), ifstream::in);
    inp3.close();
    inp4.close();
    if(inp3.fail() && inp4.fail()) {

        inp3.clear(ios::failbit);
        inp4.clear(ios::failbit);
        cout << "Writing to file ..." << outfile3.c_str() << endl;
        cout << "Writing to file ..." << outfile4.c_str() << endl;
        outp3.open(outfile3.c_str(), ofstream::out);
        outp4.open(outfile4.c_str(), ofstream::out);

//Faking magnitudes (no IGM and IGM) by using colors. Manually set m_i=24. 
		
	    for (int i=0;i<nz;i++) {		

        // Colors with IGM absorption 
        CugR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[uLSST]), (*filterArray[gLSST]), transIgm[i], lmin, lmax); 
        CgrR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[gLSST]), (*filterArray[rLSST]), transIgm[i], lmin, lmax); 
        CriR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[rLSST]), (*filterArray[iLSST]), transIgm[i], lmin, lmax); 
        CizR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[iLSST]), (*filterArray[zLSST]), transIgm[i], lmin, lmax); 
        CzyR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[zLSST]), (*filterArray[yLSST]), transIgm[i], lmin, lmax); 

	m_rR= CriR(i)+24.0;
	m_iR= 24.0;
	m_zR= 24.0-CizR(i);
	m_yR= m_zR-CzyR(i);
	m_gR= CgrR(i)+m_rR;
	m_uR= CugR(i)+m_gR;

	error_U= simData.addLSSTuError(m_uR, uVisits); 
	error_G= simData.addLSSTgError(m_gR, gVisits); 
	error_R= simData.addLSSTrError(m_rR, rVisits); 
	error_I= simData.addLSSTiError(m_iR, iVisits); 
	error_Z= simData.addLSSTzError(m_zR, zVisits); 
	error_Y= simData.addLSSTyError(m_yR, yVisits); 

           	outp3 << i+1 
                 << "    " << m_uR << "    " << m_gR << "    " << m_rR << "    " 
		 << 24.0 << "    " << m_zR << "    " << m_yR << "	" << 0.0 << endl;

           	outp4 << i+1 
                 << "    " << error_U[1] << "    " << error_G[1] << "    " << error_R[1] << "    " 
		 << error_I[1] << "    " << error_Z[1] << "    " << error_Y[1] << endl;

}


        outp3.close();
        outp4.close();
	    }

    else{
        cout << "Error...file " << outfile3.c_str() << " exists" << endl;
        cout << "Error...file " << outfile4.c_str() << " exists" << endl;
	}


    iSb=2; 

    ifstream inp5, inp6;
    ofstream outp5, outp6;
    inp5.open(outfile5.c_str(), ifstream::in);
    inp6.open(outfile6.c_str(), ifstream::in);
    inp5.close();
    inp6.close();
    if(inp5.fail() && inp6.fail()) {

        inp5.clear(ios::failbit);
        inp6.clear(ios::failbit);
        cout << "Writing to file ..." << outfile5.c_str() << endl;
        cout << "Writing to file ..." << outfile6.c_str() << endl;
        outp5.open(outfile5.c_str(), ofstream::out);
        outp6.open(outfile6.c_str(), ofstream::out);

//Faking magnitudes (no IGM and IGM) by using colors. Manually set m_i=24. 
		
	    for (int i=0;i<nz;i++) {		

        // Colors with IGM absorption 
        CugR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[uLSST]), (*filterArray[gLSST]), transIgm[i], lmin, lmax); 
        CgrR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[gLSST]), (*filterArray[rLSST]), transIgm[i], lmin, lmax); 
        CriR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[rLSST]), (*filterArray[iLSST]), transIgm[i], lmin, lmax); 
        CizR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[iLSST]), (*filterArray[zLSST]), transIgm[i], lmin, lmax); 
        CzyR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[zLSST]), (*filterArray[yLSST]), transIgm[i], lmin, lmax); 

	m_rR= CriR(i)+24.0;
	m_iR= 24.0;
	m_zR= 24.0-CizR(i);
	m_yR= m_zR-CzyR(i);
	m_gR= CgrR(i)+m_rR;
	m_uR= CugR(i)+m_gR;

	error_U= simData.addLSSTuError(m_uR, uVisits); 
	error_G= simData.addLSSTgError(m_gR, gVisits); 
	error_R= simData.addLSSTrError(m_rR, rVisits); 
	error_I= simData.addLSSTiError(m_iR, iVisits); 
	error_Z= simData.addLSSTzError(m_zR, zVisits); 
	error_Y= simData.addLSSTyError(m_yR, yVisits); 

           	outp5 << i+1 
                 << "    " << m_uR << "    " << m_gR << "    " << m_rR << "    " 
		 << 24.0 << "    " << m_zR << "    " << m_yR << "	" << 0.0 << endl;

           	outp6 << i+1 
                 << "    " << error_U[1] << "    " << error_G[1] << "    " << error_R[1] << "    " 
		 << error_I[1] << "    " << error_Z[1] << "    " << error_Y[1] << endl;

}


        outp5.close();
        outp6.close();
	    }

    else{
        cout << "Error...file " << outfile5.c_str() << " exists" << endl;
        cout << "Error...file " << outfile6.c_str() << " exists" << endl;
	}


    iSb=3; 

    ifstream inp7, inp8;
    ofstream outp7, outp8;
    inp7.open(outfile7.c_str(), ifstream::in);
    inp8.open(outfile8.c_str(), ifstream::in);
    inp7.close();
    inp8.close();
    if(inp7.fail() && inp8.fail()) {

        inp7.clear(ios::failbit);
        inp8.clear(ios::failbit);
        cout << "Writing to file ..." << outfile7.c_str() << endl;
        cout << "Writing to file ..." << outfile8.c_str() << endl;
        outp7.open(outfile7.c_str(), ofstream::out);
        outp8.open(outfile8.c_str(), ofstream::out);

//Faking magnitudes (no IGM and IGM) by using colors. Manually set m_i=24. 
		
	    for (int i=0;i<nz;i++) {		

        // Colors with IGM absorption 
        CugR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[uLSST]), (*filterArray[gLSST]), transIgm[i], lmin, lmax); 
        CgrR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[gLSST]), (*filterArray[rLSST]), transIgm[i], lmin, lmax); 
        CriR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[rLSST]), (*filterArray[iLSST]), transIgm[i], lmin, lmax); 
        CizR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[iLSST]), (*filterArray[zLSST]), transIgm[i], lmin, lmax); 
        CzyR(i)=computeColorIGM(z, (*sedArray[iSb]), (*filterArray[zLSST]), (*filterArray[yLSST]), transIgm[i], lmin, lmax); 

	m_rR= CriR(i)+24.0;
	m_iR= 24.0;
	m_zR= 24.0-CizR(i);
	m_yR= m_zR-CzyR(i);
	m_gR= CgrR(i)+m_rR;
	m_uR= CugR(i)+m_gR;

	error_U= simData.addLSSTuError(m_uR, uVisits); 
	error_G= simData.addLSSTgError(m_gR, gVisits); 
	error_R= simData.addLSSTrError(m_rR, rVisits); 
	error_I= simData.addLSSTiError(m_iR, iVisits); 
	error_Z= simData.addLSSTzError(m_zR, zVisits); 
	error_Y= simData.addLSSTyError(m_yR, yVisits); 

           	outp7 << i+1 
                 << "    " << m_uR << "    " << m_gR << "    " << m_rR << "    " 
		 << 24.0 << "    " << m_zR << "    " << m_yR << "	" << 0.0 << endl;

           	outp8 << i+1 
                 << "    " << error_U[1] << "    " << error_G[1] << "    " << error_R[1] << "    " 
		 << error_I[1] << "    " << error_Z[1] << "    " << error_Y[1] << endl;

}


        outp7.close();
        outp8.close();
	    }

    else{
        cout << "Error...file " << outfile7.c_str() << " exists" << endl;
        cout << "Error...file " << outfile8.c_str() << " exists" << endl;
	}



//        outp.close();


//} // end of for loop of some number of SB galaxies


 
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
        cout << "ERROR - TRANSLOC NOT DEFINED!" << endl;
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

double computeColorIGM(double z, GenericFunc& sed, Filter& filterX, Filter& filterY, GenericFunc& transmission, double lmin, double lmax) {

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

//double restFrameFluxIGM(GenericFunc &sed, Filter& filter, double zs, GenericFunc& transmission, double lmin, double lmax) {

double restFrameFluxIGM(GenericFunc& sed, Filter& filter, double zs, GenericFunc& transmission, double lmin, double lmax) {

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
