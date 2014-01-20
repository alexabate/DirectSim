#include "machdefs.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pexceptions.h"
#include "histos.h"
#include "perandom.h"
#include "tvector.h"
#include "cspline.h"
#include "fioarr.h"
#include "genericfunc.h"

#include "constcosmo.h"
#include "geneutils.h"
#include "schechter.h"

namespace SOPHYA {

//******* LFParameters *******************************************************//

// Constructor: reads in the LF parameters from file "LFfile"
// and stores these parameters in vectors
LFParameters::LFParameters(string LFfile, int prtvl)
: type_("file")
{
    File(LFfile, prtvl);
};


// Constructor: calculates parameter interpolation: does not have to read in 
// from file
LFParameters::LFParameters(double zmax,double dz, string type)
: type_(type)
{

    if (strcmp(type_.c_str(),"dahlen")==0)
        Dahlen(zmax,dz);
    else
        throw ParmError("Don't understand LF type!");

};


// Copy constructor: THIS DOESN'T WORK
LFParameters::LFParameters(LFParameters const& a, int id)
{

    cout << "id="<<id<<endl;
    // Dahlen constructor
    if (id<1) {
        cout <<" in copy constructor"<<endl;
        /*int ng=a.zgrid_.size();
        
        for (int i=0; i<ng; i++)
            {
            zgrid_.push_back((a.zgrid_[i]));
            mstargrid_.push_back((a.mstargrid_[i]));
            pstargrid_.push_back((a.pstargrid_[i]));
            alpgrid_.push_back((a.alpgrid_[i]));
            all_.push_back((a.all_[i]));
            early_.push_back((a.early_[i]));
            late_.push_back((a.late_[i]));
            sb_.push_back((a.sb_[i]));
            }*/
        }
    else // Read-file constructor
        {
        // stuff to be able to read in parameters from file, then store them

        /*int nc=a.zcs_.size();

        for (int i=0; i<nc; i++)
            {
            zcs_.push_back(a.zcs_[i]);
            zmins_.push_back(a.zmins_[i]);
            zmaxs_.push_back(a.zmaxs_[i]);
            AMstars_.push_back(a.AMstars_[i]);
            Aalphas_.push_back(a.Aalphas_[i]);
            Aphistars_.push_back(a.Aphistars_[i]);
            EMstars_.push_back(a.EMstars_[i]);
            Ealphas_.push_back(a.Ealphas_[i]);
            Ephistars_.push_back(a.Ephistars_[i]);
            SMstars_.push_back(a.SMstars_[i]);
            Salphas_.push_back(a.Salphas_[i]);
            Sphistars_.push_back(a.Sphistars_[i]);
            SBMstars_.push_back(a.SBMstars_[i]);
            SBalphas_.push_back(a.SBalphas_[i]);
            SBphistars_.push_back(a.SBphistars_[i]);
            }

        ZminCol_=a.ZminCol_;
        ZmaxCol_=a.ZmaxCol_;
        ACol_=a.ACol_;
        ECol_=a.ECol_;
        SCol_=a.SCol_;
        SBCol_=a.SBCol_;
        MstarCol_=a.MstarCol_;
        alphaCol_=a.alphaCol_;
        phistarCol_=a.phistarCol_;
        nZbins_=a.nZbins_;
        MstarUnits_=a.MstarUnits_;
        phistarUnits_=a.phistarUnits_;*/
        }


};


//******* LFParameters methods ***********************************************//

void LFParameters::File(string LFfile,int prtvl)
{

    // read the file into an array
    TArray<r_4> LFTable=ReadFile(LFfile);

    // if choose to print
    if (prtvl>0)
        PrintTable(LFTable);

    // store the redshift bins and Schechter parameters in vectors
    StorePars(LFTable);

};


// How the LF is interpolated and extrapolated in the paper: Dahlen et el 2008
void LFParameters::Dahlen(double zmax, double dz)
{
    // See Table 3 in Dahlen et al 2005
    // LF pars measured in 3 redshift bins (0.1<z<0.5, 0.5<z<0.75, 0.75<z<1)
    // and for 3 different galaxy types: early, late, starburst
    
    // The LF pars for the FULL SAMPLE are the only ones interpolated (ie not
    // the parameters measured per galaxy type)

    // PHISTAR
    // Take phiStar from Table 3 in Dahlen at lowest redshift bin, extrapolate to z=0
    double phistarz0=29.7e-4; // 28.1@z=0.3 extrapolated to z=0
    double phistarz2=0.63*phistarz0;    // evolution to z=2 
    double phistarz6=0.42*phistarz0;    // evolution to z=6
    // then we assume linear evolution in between
    double slope1_pstar=(phistarz2-phistarz0)/2;
    double slope2_pstar=(phistarz6-phistarz2)/(6-2);
    double c1_pstar=phistarz0;
    double c2_pstar=phistarz6-slope2_pstar*6;

    // ALPHA
    // Alpha is fixed for simplicity
    double alpha=-1.37; // value for "fixed faint end slope fit" in Table 3
    
    // MSTAR
    // Mstar is fixed to M_B=-21.22 at z=0, and then it is assumed that this 
    // value evolves to become brighter by 1 mag by z=1.5
    // Above z>1.5, Mstar is kept constant
    double mstarz0=-21.22;      // this is the value of Mstar for all gals in the lowest z bin
    double mstarz1p5=mstarz0-1.;// 1 mag brighter (the value at z=1.5)
    double mstarmaxz=mstarz1p5; // M* is constant from z>=1.5
    // then we assume linear evolution in between
    double slope1_mstar=(mstarz1p5-mstarz0)/1.5;
    double slope2_mstar=0;
    double c1_mstar=mstarz0;
    double c2_mstar=mstarmaxz;

    // build grid
    int nz=ceil(zmax/dz);
    for (int i=0; i<nz; i++) {

        double zv=dz*i; // starting at z = 0
        zgrid_.push_back(zv);

        // do Mstar interpolation
        double ms;
        if (zv<=1.5)
            ms=slope1_mstar*zv+c1_mstar;
        else 
            ms=slope2_mstar*zv+c2_mstar;
            
        // do phiStar interpolation
        double ps;
        if (zv<=2)
            ps=slope1_pstar*zv+c1_pstar;
        else 
            ps=slope2_pstar*zv+c2_pstar;

        mstargrid_.push_back(ms);
        pstargrid_.push_back(ps);
        alpgrid_.push_back(alpha);
        }

    // All galaxies LF parameters from low redshift bin
    all_.push_back(28.1e-4); all_.push_back(-21.22); all_.push_back(-1.37);
    // early type LF parameters from low redshift bin (excluding faint up turn)
    early_.push_back(15.2e-4); early_.push_back(-20.99); early_.push_back(-0.64);
    // late type LF parameters from low redshift bin
    late_.push_back(21.2e-4 ); late_.push_back(-21.00); late_.push_back(-1.35); 
    // starburst type LF parameters from low redshift bin
    sb_.push_back(42.6e-4); sb_.push_back(-18.72); sb_.push_back(-1.02); 

    // These "per type" parameters will be used to calculate the relative 
    // abundance of each of the types.

    // Mstar evolves as per the above for all galaxies

    // Basically these vectors will be returned and used in a class called 
    // TypeRatio to calculate the type fraction and its evolution with z

};


// Function that reads in the parameters from the file into an array "LFTable"
TArray<r_8> LFParameters::ReadFile(string LFfile)
{

    // get lf file location from enviroment variable
    string LFplace;
    char * plf=getenv("LFLOC");
    if (plf==NULL) {
        string emsg="ERROR LF LOCATION ENVIRONMENT VARIABLE -LFLOC-";
        emsg+=" NOT DEFINED";
        throw ParmError(emsg);
        }
    else {
        LFplace=plf;
        cout <<"    Location of LF file is "<<LFplace<<endl;
        }
    string fname=LFplace + LFfile;
    
    cout <<"    Reading in file "<<LFfile<<endl;
    ifstream ifs;
    ifs.open(fname.c_str(), ifstream::in);
    if (ifs.fail()) {
        string emsg="ERROR: failed to find luminosity function file";
        throw ParmError(emsg);
        }
    sa_size_t nr, nc;
    TArray<r_4> LFTable;
    LFTable.ReadASCII(ifs,nr,nc);
    // note that ReadASCII puts the rows as columns (ie in the 2nd dim) and
        // columns as rows (ie in the 1st dim)
    cout <<"     Size of LFTable: number of rows = "<<LFTable.SizeX();
    cout <<" (nc="<<nc<<")"<<endl;
    cout <<"                      number of cols = "<<LFTable.SizeY();
    cout <<" (nr="<<nr<<")"<<endl;
    cout <<endl;

    return LFTable;
};


// Function that stores the LF parameters in array "LFTable" in separate vectors
// this includes two vectors containing the max and min z of each redshift bin
void LFParameters::StorePars(TArray<r_4> LFTable)
{

    // "Column" (actually row in the array) numbers of values in LFTable    
    
    // Redshift range
    ZminCol_=0, ZmaxCol_=1;
    
    // Starting column number for each galaxy type: All, Early, Spiral, Starburst 
    ACol_=2, ECol_=5, SCol_=8, SBCol_=11;
    
    // Column order of Schechter parameters
    MstarCol_=0, alphaCol_=1, phistarCol_=2;
    
    // Eg. if want column number of alpha for the early types it will be:
    // ECol_ + alphaCol_ = 6

    int nr=LFTable.SizeX();
    int nc=LFTable.SizeY();
    cout <<"     Parameter table has "<< nr <<" rows and "<< nc <<" columns"<<endl;

    // storing the parameters (+z bounds) in vectors
    // the storing order is from values in the lowest z bin to the highest
    // Loop is over rows in the file (columns in the array)
    // One row for each redshift bin
    for (int i=0; i<nc; i++) {
    
        // z bounds
        zmins_.push_back(LFTable(ZminCol_,i));
        zmaxs_.push_back(LFTable(ZmaxCol_,i));
        zcs_.push_back(0.5*(zmins_[i]+zmaxs_[i]));
        
        // all galaxies pars
        AMstars_.push_back(LFTable(ACol_+MstarCol_,i));
        Aalphas_.push_back(LFTable(ACol_+alphaCol_,i));
        Aphistars_.push_back(LFTable(ACol_+phistarCol_,i));
        
        // early galaxies pars
        EMstars_.push_back(LFTable(ECol_+MstarCol_,i));
        Ealphas_.push_back(LFTable(ECol_+alphaCol_,i));
        Ephistars_.push_back(LFTable(ECol_+phistarCol_,i));
        
        // late galaxies pars
        SMstars_.push_back(LFTable(SCol_+MstarCol_,i));
        Salphas_.push_back(LFTable(SCol_+alphaCol_,i));
        Sphistars_.push_back(LFTable(SCol_+phistarCol_,i));
        
        // starburst galaxies pars
        SBMstars_.push_back(LFTable(SBCol_+MstarCol_,i));
        SBalphas_.push_back(LFTable(SBCol_+alphaCol_,i));
        SBphistars_.push_back(LFTable(SBCol_+phistarCol_,i));
        }

    // This is the number of z bins that have an LF read in from the file
    nZbins_=zmins_.size();

    // this is the actual units for these pars (alpha has no units)
    MstarUnits_="M-5log10h70";
    phistarUnits_="(Mpc/h70)^-3";
};


// function to print the parameters read in (useful for checking)
// called in the constructor if prtvl>0 (not by default)
void LFParameters::PrintTable(TArray<r_4> LFTable)
{
    int nr=LFTable.SizeX();
    int nc=LFTable.SizeY();

    cout << "     Printing LF parameters ..."<<endl<<endl;
    cout << "                         ALL               EARLY            ";
    cout << "LATE                STARBURST"<<endl;
    cout << "     zmin  zmax  [M*  alpha  phi*]  [M*  alpha  phi*]  [M*  ";
    cout << "alpha  phi*]  [M*  alpha  phi*] "<<endl;
    for (int i=0; i<nc; i++) {
        cout <<"     ";
        for (int j=0; j<nr; j++)
            cout << LFTable(j,i) << "  ";
        cout <<endl;
        }   
    cout <<endl<<endl;

};


// function to print the z bin ranges of the LF's read in (useful for checking)
void LFParameters::PrintZbins()
{
    cout << "     Redshift bins"<<endl;
    for (int i=0; i<nZbins_; i++)
        cout <<"     bin "<<i+1<<": "<<zmins_[i]<<" to "<<zmaxs_[i]<<endl;    
    cout <<endl;
};


// function to print the LF pars for a particular galaxy type, once they have
// been stored in the vectors. similar to PrintTable() (useful for checking)
void LFParameters::PrintPars(int type)
{
    if (type<0 || type >3)
        throw ParmError("Galaxy type not understood");

    cout << "     Parameters for ";
    if (type<1 && type>-1)
        cout <<"All";
    else if (type>0 && type<2)
        cout <<"Early";
    else if (type>1 && type<3)
        cout <<"Late";
    else if (type>2 && type<4)
        cout <<"Starburst";
    cout << " galaxy types in each redshift bin"<<endl;
    for (int i=0; i<nZbins_; i++) {
    
        if (type<1 && type>-1) {
            cout <<"     bin "<<i+1<<": Mstar="<<AMstars_[i];
            cout <<", alpha="<<Aalphas_[i]<<", phistar=";
            cout <<Aphistars_[i]<<endl;
            }
        else if (type>0 && type<2) {
            cout <<"     bin "<<i+1<<": Mstar="<<EMstars_[i];
            cout <<", alpha="<<Ealphas_[i]<<", phistar=";
            cout <<Ephistars_[i]<<endl;
            }
        else if (type>1 && type<3) {
            cout <<"     bin "<<i+1<<": Mstar="<<EMstars_[i];
            cout <<", alpha="<<Salphas_[i]<<", phistar=";
            cout <<Sphistars_[i]<<endl;
            }
        else if (type>2 && type<4) {
            cout <<"     bin "<<i+1<<": Mstar="<<SBMstars_[i];
            cout <<", alpha="<<SBalphas_[i]<<", phistar=";
            cout <<SBphistars_[i]<<endl;
            }
        }
    cout <<endl;
};


// function that works out how many of the LF redshift bins the galaxy survey 
// spans probably this will be all of them since LF's are usually defined only 
// 0<z<1
void LFParameters::CheckBinRange(int& i_min, int& i_max, int& nb,double zmin, double zmax)
{
// LF is defined in "nb" redshift bins as read in from the file in the 
// constructor i_min/i_max are the LF bin indices that the galaxy survey covers
// if the galaxy survey covers ALL of the LF redshift bins: i_min=0 and 
// i_max=nZbins_-1 note that it is fine if zmin<min z bin of the LF, and 
// zmax>max zbin of the LF we just have to extrapolate

    i_min=-1,i_max=-1;
    
    for (int i=0; i<nZbins_; i++) {
        if ( (zmin<zmaxs_[i])&&(i_min<0) )
            i_min=i;
        if ( zmax>zmins_[i] )
            i_max=i;
        }
    
    if (i_min<0 || i_max<0)
        throw ParmError("Bin range not set");

    if (i_min<1 && i_max>=nZbins_-1 ) { 
        cout <<"     Galaxy survey spans all LF bins"<<endl; 
        nb=nZbins_; 
        }
    else {
        cout << "     Galaxy survey spans from LF bin "<<i_min+1;
        cout << " to LF bin "<< i_max+1<<" of the "<<nZbins_<<" bins\n";
        nb=(i_max-i_min)+1;
        }
    cout <<endl;

};


// function that returns the LF pars for a given redshift bin, for a given gal 
// type (default is for "all" galaxies) 
void LFParameters::ReturnParsBini(double& Mstar, double& alpha, double& phistar,
                                int ib, int type)
{
    // ib is the LF redshift bin number
    // type refers to the galaxy type luminosity function:
    // type=0: all
    // type=1: early/elliptical
    // type=2: late/spiral
    // type=3: starburst

    if (type<0 || type >3)
        throw ParmError("Galaxy type not understood");

    // if type=0 return parameters for "all" galaxies LF
    if (type<1 && type>-1) {
        Mstar=AMstars_[ib];
        alpha=Aalphas_[ib];
        phistar=Aphistars_[ib];
        }
    // if type=1 return parameters for "early" galaxies LF
    else if (type>0 && type<2) {
        Mstar=EMstars_[ib];
        alpha=Ealphas_[ib];
        phistar=Ephistars_[ib];
        }
    // if type=2 retrun parameters for "late" galaxies LF
    else if (type>1 && type<3) {
        Mstar=SMstars_[ib];
        alpha=Salphas_[ib];
        phistar=Sphistars_[ib];
        }
    // if type=3 return parameters for "starburst" galaxies LF
    else if (type>2 && type<4) {
        Mstar=SBMstars_[ib];
        alpha=SBalphas_[ib];
        phistar=SBphistars_[ib];
        }
    else {
        stringstream ss;
        ss << type;
        string emsg = "ERROR! type value = " + ss.str() + " unknown";
        throw ParmError(emsg);
        }

};

//-----------------------CLASSES FOR CALCULATING SCHECHTER FUNCS -------------//

//******* Schechter **********************************************************//
    
// Copy constructor
Schechter::Schechter(Schechter& f)
  : phistar_(f.phistar_) , Mstar_(f.Mstar_) , alpha_(f.alpha_) , outvalue_(f.outvalue_) ,
    Mmin_(f.Mmin_) , Mmax_(f.Mmax_) , npt_(f.npt_)
{ };


void Schechter::SetOutValue(unsigned short outvalue)
// outvalue = 0 : give dn/dm
//          = 1 : give m*dn/dm
{
    if(outvalue>1) {
        cout<<"Schechter::SetOutValue: Error bad outvalue: "<<outvalue<<endl;
        throw ParmError("Schechter::SetOutValue: Error bad outvalue");
        }
    outvalue_ = outvalue;
};


double Schechter::operator() (double M) const
// Return : "phi(M) = f(M)" or "M*phi(M) = f(M)"
// think Mstar_ will alway be Mstar-5log10h units
{
    double x = pow(10,-0.4*(M-Mstar_));
    double phiM =0.4*log(10)*phistar_* pow(x,alpha_+1) * exp(-x);
    if(outvalue_) phiM *= M;  // if you want M*phi(M)
    return phiM;
};


double Schechter::Integrate()
// Integrate from Mmin to Mmax with at least npt points linear spaced
// Adaptive integration using Gauss-Legendre method
// Function is defined in geneutils 
{
     if(npt_<1) npt_ = 100;
     //double lmassmin = log10(Mmin), lmassmax = log10(Mmax);
     //double perc=0.01, dlxinc=(lmassmax-lmassmin)/npt, dlxmax=10.*dlxinc; 
     // unsigned short glorder=4;
     //double sum = IntegrateFuncLog(*this,lmassmin,lmassmax,perc,dlxinc,dlxmax,glorder);
     double perc=0.01, dlxinc=(Mmax_-Mmin_)/npt_, dlxmax=10.*dlxinc; 
     unsigned short glorder=4;
     double sum = IntegrateFunc(*this,Mmin_,Mmax_,perc,dlxinc,dlxmax,glorder);
     return sum;
};


// LUMINOSITY SCHECHTER FUNCTION: (Lmax approaching infinity)
// For Min approaching zero, integrating over dL gives:
// Ntot = phistar_*GAMMA(alpha_+1,M/Mstar_) which is useful for normalizations. 
// Gamma is the incomplete gamma function
// integrating L*dL gives:
// Ltot = phistar_*Mstar_*Gamma(alpha_+2) 

// To calculate maximum luminosity observable in survey:
// Lmax: m_limit=m_solar-2.5*log10((Lmax/L_solar)*(r_solar/r_limit)^2)
// m_limit=apparent magnitude limit of survey
// r_solar=earth-sun distance
// r_limit=distance to max survey z

/*double Schechter::IntegrateQuick()
// Integrate from Mmin to Mmax with at least npt points linear spaced
// Use simple numerical summation
{
    if (npt_ <= 0)
        npt_ = 500;
    
    Schechter sch;
    sch=*this;// THIS DOES NOT WORK
    double dM = (Mmax_ - Mmin_) / (npt_-1);
    
    double M = Mmin_;
    double sum=0;
    
    for (int i=1; i<npt_; i++, M += dM)
        sum += sch(M);
    
    return sum * dM;
}*/


// Number of galaxies over "skyarea" between "z1" and "z2"
double Schechter::NGals(SimpleUniverse su, double skyarea, double z1, double z2)
{
    // first get number density by simply integrating the Schechter function
    double ndens=Integrate();
    cout <<"     Number density = "<< ndens <<" gals/Mpc^3"<<endl;

    // calculate volume between z1 and z2 over sky area of skyarea
    double vol = Volume(su, skyarea, z1, z2);
    return vol*ndens;
};


// given the cosmological model (stored in "su") and the galaxy survey coverage ("skyarea" 
// between "z1" and "z2") this calculates the volume of the universe in comoving coordinates
double Schechter::Volume(SimpleUniverse su, double skyarea, double z1, double z2)
{

    // fraction of whole sky
    double skyfrac=skyarea/(4*PI);

    double dvdOmega;
    int npt = 1000;
    su.NumIntegVolElz(z1, z2, npt, dvdOmega);
    double DH = su.HubbleLengthMpc();
    double vol = skyfrac*4*PI*dvdOmega*DH;
    //cout <<"     Volume between z1(="<<z1<<") and z2(="<<z2<<") is "<<vol<<endl;
    return vol;
};


void Schechter::Print(string MU, string PU)
//MU=Mstar units, PU=phistar units
{
    cout <<"    Schechter::Print: phistar="<< phistar_ <<" "<<PU
         <<"  Mstar="<< Mstar_ <<" "<< MU
         <<"  alpha="<< alpha_
         <<"  (outvalue="<< outvalue_ <<" -> return ";
    if (outvalue_) cout <<"m*dn/dm)"; else cout <<"dn/dm)";
    cout << endl;
};


void Schechter::Print(void)
{
    cout <<"Schechter::Print: phistar="<< phistar_ <<" Mpc^-3"
         <<"  mstar="<< Mstar_ <<" MSol"
         <<"  alpha="<< alpha_
         <<"  (outvalue="<< outvalue_ <<" -> return ";
    if (outvalue_) cout <<"m*dn/dm)"; else cout <<"dn/dm)";
    cout << endl;
};


// vary faint limit of Schechter function integration to check effect on 
// galaxy number density
void Schechter::CheckFaintLimit(Schechter& sch)
{
    //brightest quasar ~-30 (ie brightest thing ever)
    //faintest dwarf galaxy ~-6 (ie faintest thing ever)
    
    double schmin=-30; // units of "M-5log10h70"
    TVector<r_4> schmaxv_(20);// units of "M-5log10h70"
    int schnpt=200;

    double faint_lim_start = -20;
    TVector<r_4> ngal_by_mp3_(20);
    for(int i=0;i<ngal_by_mp3_.Size();i++) {
    
        schmaxv_(i)=faint_lim_start+i*1; // this will go from -20 to -1 
        SetInteg(schmin, schmaxv_(i), schnpt);
        // get num dens of gals between schmin and current schmaxv
        ngal_by_mp3_(i)=sch.Integrate(); 
        cout <<"Faint limit="<< schmaxv_(i) <<" -> Int="<< ngal_by_mp3_(i) <<endl;
        }

    // set integration back to correct parameters
    SetInteg(Mmin_, Mmax_, npt_);
};


//******* SchechterVol *******************************************************//

double SchechterVol::Integrate(double z)
// Integrate from Mmin to Mmax with at least npt points linear spaced
// Adaptive integration using Gauss-Legendre method
// Function is defined in geneutils 
{
     if(npt_<1) npt_ = 100;

     double perc=0.01, dlxinc=(Mmax_-Mmin_)/npt_, dlxmax=10.*dlxinc; 
     unsigned short glorder=4;
     double sum = IntegrateFunc(*this,z,Mmin_,Mmax_,perc,dlxinc,dlxmax,glorder);
     return sum;
};


//******* SchechterM *********************************************************//

double SchechterM::Integrate(double z)
// Integrate from Mmin to Mmax with at least npt points linear spaced
// Adaptive integration using Gauss-Legendre method
// Function is defined in geneutils 
{
     if(npt_<1) npt_ = 100;

     double perc=0.01, dlxinc=(Mmax_-Mmin_)/npt_, dlxmax=10.*dlxinc; 
     unsigned short glorder=4;
     double sum = IntegrateFunc(*this,z,Mmin_,Mmax_,perc,dlxinc,dlxmax,glorder);
     return sum;
};


//******* SchechterZVol ******************************************************//

double SchechterZVol::Integrate()
// Integrate from zmin to zmax with at least npt points linear spaced
// Adaptive integration using Gauss-Legendre method
// Function is defined in geneutils 
{
     if(npt_<1) npt_ = 100;

     double perc=0.01, dlxinc=(zmax_-zmin_)/npt_, dlxmax=10.*dlxinc; 
     unsigned short glorder=4;
     double sum = IntegrateFunc(*this,zmin_,zmax_,perc,dlxinc,dlxmax,glorder);
     return sum;
};


//******* SchechterMassDist **************************************************//

SchechterMassDist::SchechterMassDist(Schechter sch,double massmin,double massmax,int nbinmass)
// Need a function to draw randomly N galaxies according to a Schechter function
// To optimize this (eventually) some histograms
// ATTENTION: The mass limits for the histos are given in mass but their abscissa
// are in log10(mass): histo [log10(massmin),log10(massmax)] with nbinmass bins
// If nbinmass<0 then it's the number of points by decade
  : sch_(sch) , sch_outvalue_(sch.GetOutValue())
  , massmin_(massmin) , massmax_(massmax) , nbinmass_(nbinmass)
  , ngalmin_(0) , ngalmax_(0) , nvalngal_(0)
  , ntrial_dir(0) , ntrial_tab(0)
  , hmdndm_(NULL) , tirhmdndm_(NULL)
{
    if(massmin_>massmax_  || massmin_<=0. || massmax_<=0.|| nbinmass_==0) {
        cout<<"SchechterMassDist::SchechterMassDist: error in input values"<<endl;
        throw ParmError("SchechterMassDist::SchechterMassDist: error in input values");
        }

    // Draw randomly according to Schechter function
    sch_outvalue_ = sch.GetOutValue();  // get original form of Schechter (dn/dm or m*dn/dm)
    sch.SetOutValue(1);  // want m*dN/dm
    // Bin the function in log10 of the mass
    double lnx1 = log10(massmin_), lnx2 = log10(massmax_); 
    if(nbinmass_<0) 
        nbinmass_ = int((-nbinmass_)*(lnx2-lnx1+1.));
  
    // create histogram according to mass limits and number of bins
    hmdndm_ = new Histo(lnx1, lnx2, nbinmass_); hmdndm_->ReCenterBin();
    // convert Schechter function to a histogram
    FuncToHisto(sch,*hmdndm_, true);  // true -> bin in log10(x)
    tirhmdndm_ = new FunRan(*hmdndm_,true);  // true -> histo is pdf
    sch.SetOutValue(sch_outvalue_);  // return to the initial value
};


// default constructor
SchechterMassDist::SchechterMassDist(void)
  : sch_outvalue_(0)
  , massmin_(0.) , massmax_(0.) , nbinmass_(0)
  , ngalmin_(0) , ngalmax_(0) , nvalngal_(0)
  , ntrial_dir(0) , ntrial_tab(0)
  , hmdndm_(NULL) , tirhmdndm_(NULL)
{
};


/*// delete class variables
void SchechterMassDist::Delete(void)
{
   if(hmdndm_) delete hmdndm_;
   if(tirhmdndm_) delete tirhmdndm_;
   hmass_.resize(0);
   tmass_.resize(0);
   ntrial_dir = ntrial_tab = 0;
};*/


int SchechterMassDist::SetNgalLim(int ngalmax, int ngalmin, unsigned long nalea)
// Creation of Histos for drawing from ngalmin to ngalmax galaxies
{
    int lp=2;
    ngalmin_=ngalmax_=nvalngal_=0;
    if(ngalmin<=0) ngalmin=1;
    if(ngalmax<ngalmin || ngalmax==1) return 0;
    ngalmin_ = ngalmin;
    ngalmax_ = ngalmax;
    nvalngal_ = ngalmax-ngalmin+1;
  
    if(nalea<1) nalea = 100000;

    if(lp>0) cout <<"SchechterMassDist::SetNgalLim: ngal=["
                  << ngalmin_ <<","<< ngalmax_ <<"] n="<< nvalngal_
                  <<" filling with "<< nalea <<" trials"<<endl;

    //------- Construct histo
    double lnx1 = log10(massmin_), lnx2 = log10(massmax_);
    if(lp>0) cout <<"> Creating "<< nvalngal_ <<" histos ["<< lnx1 
                  <<","<< lnx2 <<"] n="<<nbinmass_<<endl;
 
    for(int i=ngalmin_;i<=ngalmax_;i++) {
        Histo h(*hmdndm_); h.Zero();
        hmass_.push_back(h);
        }
    if(lp>1) cout<<"...number of histos is "<<hmass_.size()<<endl;

    //------- Random filling
    sch_.SetOutValue(1);  // want m*dN/dm
    int lpmod = nalea/20; if(lpmod<=0) lpmod=1;
    double s1=0., sc1=0.; unsigned long ns1=0;
    double sax=0., scax=0.; unsigned long nsax=0;
    // loop over number of histograms
    for(unsigned long ia=0; ia<nalea; ia++) {
        
        if(lp>1 && ia%lpmod==0) cout<<"...drawing "<< ia <<endl;
        double sum = 0.;
        
        for(int i=1; i<=ngalmax_; i++) {
     
            //double l10m = tirhmdndm_->Random();
            // draw log10 mass from Schechter function
            double l10m = tirhmdndm_->RandomInterp(); 
            double m = pow(10.,l10m);
            sum += m;
            s1 += m; sc1 += m*m; ns1++;
            int ipo = i-ngalmin_;
     
            if(ipo<0) continue; // don't add to histogram, keep doing sum
            
            // ATTENTION: the hist for drawing stores the log10(mean=sum_masses/ngal).
            //            Helps to have the same binning whatever ngal is
            double v = log10(sum/(double)i); // mean mass of gals drawn
            
            hmass_[ipo].Add(v);
            if(i==ngalmax) {
                sax += sum/(double)i; // mean of all
                scax += sum/(double)i*sum/(double)i; 
                nsax++;
                }
            }
        }
    sch_.SetOutValue(sch_outvalue_);  // Return ot the initial value

   if(ns1>1) {
       s1 /= ns1; sc1 = sc1/ns1 - s1*s1;
       cout <<"...Mean mass for ngal=1: "<< s1 <<" ("<< log10(fabs(s1)) <<")"
            <<" s="<< sqrt(fabs(sc1)) <<" (ntrials "<< ns1 <<")"<<endl;
       }
    if(nsax>1) {
        sax /= nsax; scax = scax/nsax - sax*sax;
        cout <<"...Mean mass for ngal="<< ngalmax_ <<": "<< sax 
             <<" ("<< log10(fabs(sax)) <<")"
             <<" s="<< sqrt(fabs(scax)) <<" (ntrials "<< nsax <<")"<<endl;
        }

    //------- Generation of random drawing classes and histos for checking
    if(lp>0) cout <<"> Creating "<< nvalngal_ <<" FunRan"<<endl;
 
    for(unsigned int i=0;i<hmass_.size();i++) {
        FunRan t(hmass_[i],true);
        tmass_.push_back(t);
        }
    if(lp>1) cout <<"...number of funran is "<< tmass_.size() <<endl;

    return nvalngal_;
};


Histo SchechterMassDist::GetHisto(int i) const
{
    if(i<0 || i>=nvalngal_) {
        cout<<"SchechterMassDist::GetHisto: error in input values"<<endl;
        throw ParmError("SchechterMassDist::GetHisto: error in input values");
        }
    return hmass_[i];
};


FunRan SchechterMassDist::GetFunRan(int i) const
{
    if(i<0 || i>=nvalngal_) {
        cout<<"SchechterMassDist::GetFunRan: error in input values"<<endl;
        throw ParmError("SchechterMassDist::GetFunRan: error in input values");
        }
    return tmass_[i];
};


double SchechterMassDist::TirMass(int ngal)
{
    if(ngal<1) return 0.;

    int ipo = IndexFrNGal(ngal);
    double masse_des_ngal = 0.;
    if(ipo<0) {  // No drawing from the histos for this number of galaxies
        for(long i=0; i<ngal; i++) {  // Draw ngal times from the Schechter
            // double lm = tirhmdndm_->Random();
            double lm = tirhmdndm_->RandomInterp();
            masse_des_ngal += pow(10.,lm);  // ATTENTION draws log10(mass)
            }
       ntrial_dir++;
       } 
   else {
       // ATTENTION the histos store the log10(mean=sum_masses/ngal)
       //double lmngal = tmass_[ipo].Random();
       double lmngal = tmass_[ipo].RandomInterp();
       masse_des_ngal = pow(10.,lmngal) * ngal;
       ntrial_tab++;
       }

    return masse_des_ngal;
};


void SchechterMassDist::Print(void)
{
    cout <<"SchechterMassDist::Print: mass=["<< massmin_ <<","<< massmax_ 
         <<"] n="<<nbinmass_<<endl;
    cout <<"ngal=["<< ngalmin_ <<","<< ngalmax_ <<"] n="<< nvalngal_ <<endl;
    sch_.Print();
};


void SchechterMassDist::PrintStatus(void)
{
    cout <<"SchechterMassDist::PrintStatus: number of trials: direct="<< ntrial_dir
         <<" tabulated="<< ntrial_tab <<endl;
};


void SchechterMassDist::WritePPF(string ppfname)
{
    char str[64];
    cout <<"SchechterMassDist::WritePPF into "<< ppfname <<endl;
    POutPersist pos(ppfname.c_str());

    double nstar,mstar,alpha;
    sch_.GetParam(nstar,mstar,alpha);
    TVector<r_8> tdum(20); tdum = 0.;
    tdum(0) = nstar;
    tdum(1) = mstar;
    tdum(2) = alpha;
    tdum(3) = sch_outvalue_;
    tdum(4) = massmin_;
    tdum(5) = massmax_;
    tdum(6) = nbinmass_;
    tdum(7) = ngalmin_;
    tdum(8) = ngalmax_;
    tdum(9) = nvalngal_;
    tdum(10) = hmass_.size();
    tdum(11) = tmass_.size();
    pos << PPFNameTag("SMDparam") << tdum;

    pos << PPFNameTag("SMDhmdndm") << *hmdndm_;
    pos << PPFNameTag("SMDtirhmdndm") << *tirhmdndm_;

    if(hmass_.size()>0) {
        for(unsigned int i=0;i<hmass_.size();i++) {
            sprintf(str,"SMDh%d",NGalFrIndex(i));
            pos << PPFNameTag(str) << hmass_[i];
            }
        }

    if(tmass_.size()>0) {
        for(unsigned int i=0;i<tmass_.size();i++) {
            sprintf(str,"SMDt%d",NGalFrIndex(i));
            Histo hdum(tmass_[i]);
            pos << PPFNameTag(str) << hdum;
            }
        }

};


void SchechterMassDist::ReadPPF(string ppfname)
{
    //Delete(); // De-allocate if already filled!

    char str[64];
    cout <<"SchechterMassDist::ReadPPF from "<< ppfname <<endl;
    PInPersist pis(ppfname.c_str());

    TVector<r_8> tdum;
    pis >> PPFNameTag("SMDparam") >> tdum;
    sch_.SetParam(tdum(0),tdum(1),tdum(2));
    sch_.SetOutValue((unsigned short)(tdum(3)+0.1));
    massmin_ = tdum(4);
    massmax_ = tdum(5);
    nbinmass_ = int(tdum(6)+0.1);
    ngalmin_ = int(tdum(7)+0.1);
    ngalmax_ = int(tdum(8)+0.1);
    nvalngal_ = int(tdum(9)+0.1);
    unsigned int nhmass = (unsigned int)(tdum(10)+0.1);
    unsigned int ntmass = (unsigned int)(tdum(11)+0.1);

    {
    Histo hdum;
    pis >> PPFNameTag("SMDhmdndm") >> hdum;
    hmdndm_ = new Histo(hdum);
    pis >> PPFNameTag("SMDtirhmdndm") >> hdum;
    tirhmdndm_ = new FunRan(hdum,false);
    }

    if(nhmass>0) {
        for(unsigned int i=0;i<nhmass;i++) {
            sprintf(str,"SMDh%d",NGalFrIndex(i));
            Histo hdum;
            pis >> PPFNameTag(str) >> hdum;
            hmass_.push_back(hdum);
            }
        }

    if(ntmass>0) {
        for(unsigned int i=0;i<ntmass;i++) {
            sprintf(str,"SMDt%d",NGalFrIndex(i));
            Histo hdum;
            pis >> PPFNameTag(str) >> hdum;
            FunRan fdum(hdum,false);
            tmass_.push_back(fdum);
            }
        }

};


////////////////////////////////////////////////////////////////////////////////
//******************* Check Function *****************************************//
////////////////////////////////////////////////////////////////////////////////

bool IsCompatible(Schechter& sch1, Schechter& sch2, double eps)
// compare the differences to eps pres
{
    if(eps<=0.) eps=1.e-4;
    double nstar1,mstar1,alpha1;
    sch1.GetParam(nstar1,mstar1,alpha1);
    double nstar2,mstar2,alpha2;
    sch2.GetParam(nstar2,mstar2,alpha2);

    // nstar and mstar are never nul
    if(fabs(nstar1-nstar2)>fabs(nstar1+nstar2)/2.*eps) return false;
    if(fabs(mstar1-mstar2)>fabs(mstar1+mstar2)/2.*eps) return false;

    // alpha maybe nul
    if(fabs(alpha1)<1.e-100 && fabs(alpha2)<1.e-100 && fabs(alpha1-alpha2)>eps) return false;
    if(fabs(alpha1-alpha2)>fabs(alpha1+alpha2)/2.*eps) return false;
    
    return true;
};


}  // End namespace SOPHYA
