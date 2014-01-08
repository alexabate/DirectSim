#include "machdefs.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "pexceptions.h"

#include "constcosmo.h"
#include "geneutils.h"
#include "pkspectrum.h"

namespace SOPHYA {


//****************** TransferEH **********************************************//

// From Eisenstein & Hu ApJ 496:605-614 1998 April 1 (or astro-ph/9709112)
TransferEH::TransferEH(double h100, double OmegaCDM0, double OmegaBaryon0,
                                             double tcmb, bool nobaryon, int lp)
: lp_(lp) ,
  Oc_(OmegaCDM0) , Ob_(OmegaBaryon0) , h100_(h100) , tcmb_(tcmb) ,
  nobaryon_(nobaryon) , nooscenv_(0), retpart_(ALL)
{
    zero_();
    Init_();
};


TransferEH::TransferEH(TransferEH& tf)
: lp_(tf.lp_) ,
  Oc_(tf.Oc_) , Ob_(tf.Ob_) , h100_(tf.h100_) , tcmb_(tf.tcmb_) ,
  nobaryon_(tf.nobaryon_) , nooscenv_(tf.nooscenv_), retpart_(tf.retpart_)
{
    zero_();
    Init_();
};


void TransferEH::zero_(void)
{
    th2p7_=zeq_=keq_=zd_=Req_=Rd_=s_=ksilk_=alphac_=betac_=bnode_
          =alphab_=betab_=alphag_=sfit_=kpeak_=1.e99;
};


void TransferEH::Init_(void)
{

    O0_ = Oc_ + Ob_;
    //if(nobaryon_) {O0_ = Oc_; Ob_ = 0.;} //AA: don't want to set Ob=0 here, fixed Sept2011
    double H0 = 100. * h100_, h2 = h100_*h100_;
    if(lp_) { cout <<"h100="<< h100_ <<" H0="<< H0 <<") Omatter="<< O0_;
              cout <<" Ocdm="<< Oc_ <<" Ob="<< Ob_ <<endl; }

    if(tcmb_<0.) tcmb_ = T_CMB_K;
    th2p7_ = tcmb_/2.7;
    double th2p7P4 = th2p7_*th2p7_*th2p7_*th2p7_;
    if(lp_) cout <<"tcmb = "<< tcmb_ <<" K = "<< th2p7_ <<" *2.7K "<<endl;

    // formula 2 p 606
    zeq_ = 2.50e4 * O0_ * h2 / th2p7P4;
    if(lp_) cout <<"zeq = "<< zeq_ <<" (redshift of matter-radiation equality)"<<endl;

    // formula 3 p 607
    // (attention here C=1 : H0 -> H0/C if we use the first formula)
    //  keq_ = sqrt(2.*O0_*H0*H0*zeq_) / SpeedOfLight_Cst;
    keq_ = 7.46e-2 * O0_ * h2 / (th2p7_*th2p7_);
    if(lp_) cout <<"keq = "<< keq_ <<" Mpc^-1 (scale of equality)"<<endl;
    
    //cout <<"Oc="<<Oc_<<", Ob="<<Ob_<<endl;
    //cout << "Oc+Ob="<<O0_<<",  h^2="<< h2 <<",Omh^2="<<O0_ * h2<<", TCMB/2.7="<<th2p7_<<endl;
    //cout<<"keq = "<<keq_<<" Mpc^-1 (scale of equality)"<<endl;
    
    // Stops here if there's no baryons
    if(nobaryon_) return;

    // formula 4 p 607
    double b1_eq4 = 0.313*pow(O0_*h2,-0.419)*(1. + 0.607*pow(O0_*h2,0.674));
    double b2_eq4 = 0.238*pow(O0_*h2,0.223);
    zd_ = 1291. * pow(O0_*h2,0.251) / (1.+0.659* pow(O0_*h2,0.828))
              * (1. + b1_eq4*pow(Ob_*h2,b2_eq4));
    if(lp_) cout <<"zd = "<< zd_ <<" (Redshift of drag epoch)"<<endl;

    // formula 5 page 607    (R = 3*rho_baryon/4*rho_gamma)
    Req_ = 31.5*Ob_*h2 / th2p7P4 * (1.e3/zeq_);
    //WARNING: W.Hu's code (tf_fit.c) is in disagreement with the article: zd -> (1+zd)
    Rd_  = 31.5*Ob_*h2 / th2p7P4 * (1.e3/zd_);
    //in tf_fit.c: Rd_  = 31.5*Ob_*h2 / th2p7P4 * (1.e3/(1.+zd_));
    if(lp_) {
        cout <<"Req = "<< Req_ <<" Rd = "<< Rd_
             <<" (Photon-baryon ratio at equality/drag epoch)"<<endl;
        cout <<"Sound speed at equality "<< 1./sqrt(3.*(1.+Req_))
             <<", at drag "<< 1./sqrt(3.*(1.+Rd_)) <<" in unit of C"<<endl;
        }

    // formula 6 p 607
    s_ = 2./(3.*keq_) * sqrt(6./Req_)
         * log( (sqrt(1.+Rd_) + sqrt(Rd_+Req_)) / (1.+sqrt(Req_)) );
    if(lp_) cout <<"s = "<< s_ <<" Mpc (sound horizon at drag epoch)"<<endl;

    // formula 7 page 607
    ksilk_ = 1.6*pow(Ob_*h2,0.52)*pow(O0_*h2,0.73) * (1. + pow(10.4*O0_*h2,-0.95));
    if(lp_) cout <<"ksilk = "<< ksilk_ <<" Mpc^-1 (silk damping scale)"<<endl;

    // formulas 10 page 608
    double a1 = pow(46.9*O0_*h2,0.670) * (1. + pow(32.1*O0_*h2,-0.532));
    double a2 = pow(12.0*O0_*h2,0.424) * (1. + pow(45.0*O0_*h2,-0.582));
    alphac_ = pow(a1,-Ob_/O0_) * pow(a2,-pow(Ob_/O0_,3.));
    double b1 = 0.944 / (1. + pow(458.*O0_*h2,-0.708));
    double b2 = pow(0.395*O0_*h2,-0.0266);
    betac_ = 1 / ( 1. + b1*(pow(Oc_/O0_,b2) - 1.) );
    if(lp_) cout <<"alphac = "<< alphac_ <<" betac = "<< betac_
                 <<" (CDM suppression/log shift)"<<endl;

    // formula 23 page 610
    bnode_ = 8.41 * pow(O0_*h2,0.435);
    if(lp_) cout<<"bnode = "<< bnode_ <<" (sound horizon shift)"<<endl;

    // formula 14 page 608
    //WARNING: W.Hu's code (tf_fit.c) is in disagreement with the article: (1+zeq) -> zeq
    double y = (1.+zeq_)/(1.+zd_);
    //in tf_fit.c: double y = zeq_/(1.+zd_);
    double s1py = sqrt(1.+y);
    double Gy = y*( -6.*s1py + (2.+3.*y)*log((s1py+1.)/(s1py-1.)) );
    alphab_ = 2.07*keq_*s_*pow(1.+Rd_,-3./4.)*Gy;

    // formula 24 page 610
    betab_ = 0.5 + Ob_/O0_
           + (3.-2.*Ob_/O0_) * sqrt(pow(17.2*O0_*h2,2.) + 1.);
    if(lp_) cout <<"alphab = "<< alphab_ <<" betab = "<< betab_
                 <<" (Baryon suppression/envelope shift)"<<endl;

    // formula 31 page 612
    alphag_ = 1.
            - 0.328*log(431.*O0_*h2)*Ob_/O0_
            + 0.38*log(22.3*O0_*h2)*pow(Ob_/O0_,2.);
    if(lp_) cout <<"alphag = "<< alphag_ <<" (gamma suppression in approximate TF)"<<endl;

    // The approximate value of the sound horizon, formula 26 page 611
    sfit_ = 44.5*log(9.83/(O0_*h2)) / sqrt(1.+10.*pow(Ob_*h2,3./4.));  // Mpc
    if(lp_) cout <<"sfit="<< sfit_ <<" Mpc (fit to sound horizon)"<<endl;

    // Position of the 1st acoustic peak, formula 25 page 611
    kpeak_ = 5*M_PI/(2.*sfit_) * (1.+0.217*O0_*h2);  // 1/Mpc
    if(lp_) cout <<"kpeak="<< kpeak_ <<" Mpc^-1 (fit to wavenumber of first peak)"<<endl;

    return;
};


bool TransferEH::SetParTo(double h100,double OmegaCDM0,double OmegaBaryon0)
// To change the values of the parameters (a possible re-initialisation follows)
// If h100,Omega...<=0. there's no change however, keeps the old value
{
    bool haschanged = false;

    if(h100>0.) { h100_ = h100; haschanged = true; }
    if(OmegaCDM0>0.) { Oc_ = OmegaCDM0; haschanged = true; }
    if(OmegaBaryon0>0.) { Ob_ = OmegaBaryon0; haschanged = true; }

    // recalculate the initialisations
    if(haschanged) Init_();

    return haschanged;
};


void TransferEH::SetNoOscEnv(unsigned short nooscenv)
// To obtain an approximate form of the non-oscillatory part of the transfer function
// nooscenv = 0 : use the baryon oscillatory part of transfert function (FULL TF)
// nooscenv = 1 : use approx. paragraph 3.3 p610 (middle of right column)
//                Replace  j0(k*stilde)  ->  [1+(k*stilde)^4]^(-1/4)
// nooscenv = 2 : use formulae 29+30+31 page 612 [NO WIGGLES TRANSFER FUNCTION]
//                The value of an approximate transfer function that captures
//                the non-oscillatory part of a partial baryon transfer function.
//                In other words, the baryon oscillations are left out,
//                but the suppression of power below the sound horizon is included.
{
    if(nooscenv!=1 && nooscenv!=2) nooscenv = 0;
    nooscenv_ = nooscenv;
};


double TransferEH::T0tild(double k, double alphac, double betac) const
{
    // formula 10 p 608
    //double q = k*th2p7_*th2p7_/(O0_*h100_*h100_);
    double q = k/(13.41*keq_);
    // formula 20 p 610
    double C = (14.2/alphac) + 386./(1.+69.9*pow(q,1.08));
    // formula 19 p 610
    double x = log(M_E+1.8*betac*q);
    return x / (x + C*q*q);
}

double TransferEH::operator() (double k) const
{

    // --- For zero baryon
    //  or for smooth function without BAO
    if(nobaryon_  || nooscenv_ == 2) {
        double gamma = O0_*h100_;
        // Calculation of Gamma_eff, formula 30 page 612 (for smooth function)
        if( nobaryon_==false && nooscenv_ == 2 )
            gamma = O0_*h100_*(alphag_ + (1.-alphag_)/(1.+pow(0.43*k*sfit_,4.))); // Gamma_eff
    
        // formula 28 page 612 : which is equivalent to:
        //    q = k / h100_ * th2p7_*th2p7_ / gamma;
        // which is equivalent to:
        //    q = k / (13.41 * keq)                   for Ob=0
        //    q = k / (13.41 * keq) * (O0*h/Gamma)    for the smooth spectrum
        // The results are slightly different because of the approximate values
        // of the numerical constants: take like W.Hu's code (tf_fit.c)
        //double q = k / h100_ * th2p7_*th2p7_ / gamma;  // Mpc^-1
        
        // AA: fixed q calculation for ZERO BARYONS
        double q;
        if (nobaryon_)
            q = k/(13.41*keq_);  // Mpc^-1
        else
	        q = k/(13.41*keq_) * (O0_*h100_/gamma);  // Mpc^-1
	        
        // formulas 29 page 612
        double l0 = log(2.*M_E + 1.8*q);
        double c0 = 14.2 + 731./(1.+62.5*q);
        return l0 / (l0 + c0*q*q);
        }

    // --- For CDM + Baryons
    // --- CDM
    double f = 1. / (1. + pow(k*s_/5.4,4.));
    double Tc = f*T0tild(k,1.,betac_) + (1.-f)*T0tild(k, alphac_, betac_);
    if(retpart_ == CDM) return Tc; 

    // --- Baryons
    // formula 22 page 610
    double stilde, ksbnode = k*s_/bnode_;
    if(ksbnode<0.001) 
        stilde =s_ * ksbnode;
    else   
        stilde = s_ / pow(1. + pow(1./ksbnode,3.), 1./3.);
  
    // formula 21 page 610
    double j0kst = 0.;
    if(nooscenv_ == 1) 
        j0kst = pow(1.+pow(k*stilde,4.) , -1./4.);  //smooth without BAO
    else {
        double x = k*stilde;
        if(x<0.01) 
            j0kst = 1. - x*x/6.*(1.-x*x/20.);
        else 
            j0kst = sin(x)/x;
        //cout<<"DEBUG: k="<<k<<" stilde="<<stilde<<" x="<<x<<" j0kst="<<j0kst<<endl;
        }
        
    double Tb = T0tild(k,1.,1.) / (1. + pow(k*s_/5.2,2.));
    Tb += alphab_/(1.+pow(betab_/(k*s_),3.)) * exp(-pow(k/ksilk_,1.4));
    Tb *= j0kst;
  
    if(retpart_ == BARYON) 
        return Tb;

    // --- Total
    double T = (Ob_/O0_)*Tb + (Oc_/O0_)*Tc;

    return T;
};


double TransferEH::KPeak(void)
// Position of 1st acoustic peak
{
    if(nobaryon_) return -1.;
    return kpeak_;
};


//******************* TransferTabulate ***************************************//

// For reading transfer function from a CAMB output file
TransferTabulate::TransferTabulate(void)
: kmin_(1.) , kmax_(-1.) , interptyp_(0)
{
    k_.resize(0);
    tf_.resize(0);
};


TransferTabulate::TransferTabulate(TransferTabulate& tf)
: kmin_(tf.kmin_) , kmax_(tf.kmax_) , interptyp_(tf.interptyp_) , k_(tf.k_) , tf_(tf.tf_)
{
};


void TransferTabulate::SetInterpTyp(int typ)
// see comment in InterpTab
{
    if(typ<0) typ=0; 
    else if(typ>2) typ=2;
    
    interptyp_ = typ;
};


int TransferTabulate::ReadCMBFast(string filename, double h100, 
                                          double OmegaCDM0, double OmegaBaryon0)
{
    FILE *file = fopen(filename.c_str(),"r");
    if(file==NULL) return -1;
    
    cout <<"TransferTabulate::ReadCMBFast: fn="<< filename <<" h100="<< h100
         <<" OmegaCDM0="<< OmegaCDM0 <<" OmegaBaryon0="<< OmegaBaryon0 <<endl;

    const int lenline = 512;
    char *line = new char[lenline];

    k_.resize(0); tf_.resize(0);
    double tmax = -1.;
    while ( fgets(line,lenline,file) != NULL ) {
     
        double k,tc,tb,tf;
        sscanf(line,"%lf %lf %lf",&k,&tc,&tb);
        k *= h100;     // convert h Mpc^-1  ->  Mpc^-1
        tf = (OmegaCDM0*tc+OmegaBaryon0*tb)/(OmegaCDM0+OmegaBaryon0);
        if(tf>tmax) 
            tmax = tf;
        k_.push_back(k);
        tf_.push_back(tf);
        }

    cout<<"TransferTabulate::ReadCMBFast: nread="<<tf_.size()<<" tf_max="<<tmax<<endl;
    delete [] line;
    if(tf_.size()==0) 
        return (int)tf_.size();

    for (unsigned int i=0; i<tf_.size(); i++) 
        tf_[i] /= tmax;

    return (int)tf_.size();
};


int TransferTabulate::ReadCAMB(string filename, double h100)
{
    FILE *file = fopen(filename.c_str(),"r");
    if(file==NULL) return -1;
    
    cout <<"TransferTabulate::ReadCAMB: fn="<< filename <<" h100="<< h100 <<endl;

    const int lenline = 512;
    char *line = new char[lenline];

    k_.resize(0); tf_.resize(0);
    double tmax = -1.;
    while ( fgets(line,lenline,file) != NULL ) {
        double k,tcdm,tbar,tph,trel,tnu,ttot, tf;
        sscanf(line,"%lf %lf %lf %lf %lf %lf %lf",&k,&tcdm,&tbar,&tph,&trel,&tnu,&ttot);
        k *= h100;     // convert h Mpc^-1  ->  Mpc^-1
        tf = ttot;
        if(tf>tmax) 
            tmax = tf;
        k_.push_back(k);
        tf_.push_back(tf);
        }

    cout<<"TransferTabulate::ReadCAMB nread="<<tf_.size()<<" tf_max="<<tmax<<endl;
    delete [] line;
    if(tf_.size()==0) 
        return (int)tf_.size();

    for(unsigned int i=0; i<tf_.size(); i++) 
        tf_[i] /= tmax;

    return (int)tf_.size();
};


//********************* GrowthFactor *****************************************//
double GrowthFactor::DsDz(double z, double)
{
    cout <<"GrowthFactor::DsDz_Error not implemented"<<endl;
    throw AllocationError(" GrowthFactor::DsDz_Error not implemented");
}


//********************* GrowthEH *********************************************//
// From Eisenstein & Hu ApJ 496:605-614 1998 April 1
// To have D(z) = 1/(1+z) make: OmegaMatter0=1 OmegaLambda0=0
GrowthEH::GrowthEH(double OmegaMatter0,double OmegaLambda0)
  : O0_(OmegaMatter0) , Ol_(OmegaLambda0)
{
    if(OmegaMatter0==0.) {
        cout <<"GrowthEH::GrowthEH_Error: bad OmegaMatter0  value : "
             << OmegaMatter0 <<endl;
        throw ParmError("GrowthEH::GrowthEH_Error:  bad OmegaMatter0  value");
        }
};


double GrowthEH::operator() (double z) const
// see Formulae A4 + A5 + A6 page 614
{
    z += 1.;
    double z2 = z*z, z3 = z2*z;

    // Calculation of the normalisation (for z=0 -> growth=1.)
    double D1z0 = pow(O0_,4./7.) - Ol_ + (1.+O0_/2.)*(1.+Ol_/70.);
    D1z0 = 2.5*O0_ / D1z0;

    // Calculation of the growthfactor vs z
    double Ok = 1. - O0_ - Ol_;
    double den = Ol_ + Ok*z2 + O0_*z3;
    double o0z = O0_ *z3 / den;
    double olz = Ol_ / den;

    double D1z = pow(o0z,4./7.) - olz + (1.+o0z/2.)*(1.+olz/70.);
    D1z = 2.5*o0z / z / D1z;

    return D1z / D1z0;
};


double GrowthEH::DsDz(double z, double dzinc)
// y-y0 = a*(x-x0)^2 + b*(x-x0)
// dy = a*dx^2 + b*dx   with  dx = x-x0
// dy'(dx) = 2*a*dx + b -> for x=x0 on a dy'(0) = b
//
// x-x0 == z1 or z2 - z
// y-y0 == g(z1) or g(z2) - g(z)
// no idea what is happening here
{
    if(z<0. || dzinc<=0.) {
        cout<<"GrowthEH::DsDz_Error: z<0 or dzinc<=0. !"<<endl;
        throw ParmError("GrowthEH::DsDz_Error: z<0 or dzinc<=0. !");
        }

    double z1, z2;
    if(z>dzinc/2.) {
        // case where z is sufficiently far from zero
        // resolution with 2 points manages x0=z
        z1 = z - dzinc; if(z1<0.) z1 = 0.;
        z2 = z + dzinc;
        } 
    else {
        // case where z is close to zero
        // resolution with 2 points above, x0=z1
        z1 = z + dzinc;
        z2 = z + 2.*dzinc;
        }

    double gz  = (*this)(z); // growth at redshift z
    double dgz1 = (*this)(z1) - gz; // growth at z1 - growth at z
    double dgz2 = (*this)(z2) - gz; // growth at z2 - growth at z

    z1 -= z;
    z2 -= z;
    return (dgz2*z1*z1 - dgz1*z2*z2)/(z1*z2*(z1-z2));
};


void GrowthEH::SetParTo(double OmegaMatter0, double OmegaLambda0)
{
    if(OmegaMatter0>0.) O0_ = OmegaMatter0;
    Ol_ = OmegaLambda0;
};


bool GrowthEH::SetParTo(double OmegaMatter0)
// identical to the above without changing OmegaLambda0
{
    if(OmegaMatter0<=0.) return false;
    O0_ = OmegaMatter0;
    return true;
};

//********************* GrowthFN *********************************************//
// new growth function, more general, also probably inaccurate

// Constructor, for dark energy model
GrowthFN::GrowthFN(double OmegaM, double OmegaL, double w0, double wa, int na, int prt)
: OmegaM_(OmegaM) , OmegaL_(OmegaL) , w0_(w0) , wa_(wa)
{
    if(OmegaM==0.) {
        cout<<"GrowthFN::GrowthFN_Error: bad OmegaM  value : "<<OmegaM<<endl;
        throw ParmError("GrowthFN::GrowthEH_Error:  bad OmegaM value");
        }  

    CalcDEDone_=false;

    // Curvature
    OmegaK_ = 1-OmegaM_-OmegaL_;

    if ( abs(w0_+1)<1e-6 && abs(wa_)<1e-6) { // If parameters==LCDM use LCDM growth func
	    DEGrowth_=false;
  	    LCDMGrowth_=true;
	    if (prt>0) { 
	        cout <<"     w0 = "<<w0_<<", wa_ = "<<wa_<<" so using LCDM growth";
		    cout <<" function"<<endl; 
		    }
	    }
    else { // for now just do w=const
	    DEGrowth_=true;
  	    LCDMGrowth_=false;
	    if (prt>0)
		    cout <<"     Using DE growth function"<<endl;
	    // This function makes an interpolation table
	    // AND calculates growth at z=0 too
  	    CalcDE(na);
	    }
	
    // Currently if DEGrowth=true this calculation cheats with a polyfit
    // the the growth at z=0 calculated by WLCode
    GrowthZ0();

};


// Constructor for LCDM model
GrowthFN::GrowthFN(double OmegaM, double OmegaL)
: OmegaM_(OmegaM) , OmegaL_(OmegaL)
{
    if(OmegaM==0.) {
        cout<<"GrowthFN::GrowthFN_Error: bad OmegaM  value : "<< OmegaM <<endl;
        throw ParmError("GrowthFN::GrowthEH_Error:  bad OmegaM value");
        }

    w0_=-1;
    wa_=0;

    CalcDEDone_=false;

    // Curvature
    OmegaK_ = 1-OmegaM_-OmegaL_;

    DEGrowth_=false;
    LCDMGrowth_=true;
    GrowthZ0(); // calculate growth at z=0

};


// Copy constructor
GrowthFN::GrowthFN(GrowthFN& d1)
: OmegaM_(d1.OmegaM_) , OmegaL_(d1.OmegaL_) , w0_(d1.w0_) , wa_(d1.wa_)
{
    DEGrowth_=d1.DEGrowth_;
    LCDMGrowth_=d1.LCDMGrowth_;
    GrowthZ0(); // calculate growth at z=0
};


// LCDM calculation: w0=-1, wa=0
double GrowthFN::GrowthLCDM(double z) const
// see Formulae A4 + A5 + A6 page 614
// in Eisenstein and Hu 1998, the
// original ref is Carroll, Press & Turner 1992
{
    z += 1.;
    double z2 = z*z, z3 = z2*z;

    // Calculation of the growthfactor vs z
    double Ok = 1. - OmegaM_ - OmegaL_;
    double den = OmegaL_ + Ok*z2 + OmegaM_*z3;
    double o0z = OmegaM_ *z3 / den;
    double olz = OmegaL_ / den;

    double D1z = pow(o0z,4./7.) - olz + (1.+o0z/2.)*(1.+olz/70.);
    D1z = 2.5*o0z / z / D1z;

    // normalised to 1 at z=0
    // D1z0_ was already calculated in constructor
    return D1z / D1z0_;
};


// DE calculation
double GrowthFN::GrowthDE(double z) const
{
    if (!CalcDEDone_) {
        string emsg = "GrowthFN::GrowthDE Have not done Dark energy growth calculation!"; 
	    throw ParmError(emsg);
	    }
    double gro = DEgrofunc_(z);
    return gro;
};


// Growth at redshift zero
void GrowthFN::GrowthZ0()
{

    if(LCDMGrowth_) {
	    double D1z0= pow(OmegaM_,4./7.) - OmegaL_ + (1.+OmegaM_/2.)*(1.+OmegaL_/70.);
  	    D1z0_ = 2.5*OmegaM_ / D1z0;
	    }
    else if (DEGrowth_)
	    D1z0_=GrowthDEZ0Fit(OmegaM_,w0_);
    else
	    D1z0_ =-999;

};

// this comes from a 5 degree 2D polynomial fit to the distribution of growth at z=0
// calculated by WLCcode for varying OmegaM,w0
double GrowthFN::GrowthDEZ0Fit(double omegam, double w0)
{

    vector<double> coeffs;
    coeffs.push_back(-0.4919);
    coeffs.push_back(-1.3517);
    coeffs.push_back(-0.5338);
    coeffs.push_back(0.2570);
    coeffs.push_back(0.3374);
    coeffs.push_back(0.0937);
    coeffs.push_back(4.3718);
    coeffs.push_back(0.5802);
    coeffs.push_back(0.0836);
    coeffs.push_back(0.2434);
    coeffs.push_back(0.0890);
    coeffs.push_back(-11.7486);
    coeffs.push_back(2.0054);
    coeffs.push_back(1.7116);
    coeffs.push_back(0.2848);
    coeffs.push_back(22.1383);
    coeffs.push_back(-0.9427);
    coeffs.push_back(-0.5363);
    coeffs.push_back(-20.0254);
    coeffs.push_back(0.0403);
    coeffs.push_back(6.8244);

    double order = 0.5 * (sqrt(8*coeffs.size()+1) - 3);// should equal 5
    double dz0 = 0;
    int column = 0;
    for (int xpower=0; xpower<(order+1); xpower++)// from 0 to 5
        for (int ypower = 0; ypower<(order-xpower+1); ypower++) {
	   
            dz0 += (coeffs[column]*pow(omegam,xpower)*pow(w0,ypower));
            column++;
	        }

    return dz0;

};


// not convinced this is accurate
void GrowthFN::CalcDE(int na)
{
    // value for na is set in constructor
    // the default value is 1,000,000
    // this gives a 5% accuracy to the calculation
    // need to implement a better ODE integrator 
    // to get a higher accuracy with less steps

    //vector<double> deltavals,zvals;
    vector<double> deltmp,ztmp,dddashtmp1,dddashtmp2;

    // set up initial conditions
    double astart=1e-4, astop=1;// start value same as WLCode
    double da=(astop-astart)/(na-1);
    int ia=0;
    // delta \propto a at matter domination
    // so set this as initial condition
    // d delta /da then = 1
    double delta=astart; 
    double deltadash=1;// equivalent value to WLCode, but 4.6 gives "better" results
    double a=astart;
    // set index zero values
    deltmp.push_back(delta);
    ztmp.push_back(1/astart-1);
    dddashtmp1.push_back(-1);
    dddashtmp2.push_back(-1);

    double deltaold,deltadashold,aold;

    for (ia=1;ia<na;ia++) { // start at index 1, index 0 already set
		
        // save previous redshift delta,deltadash,a values
        deltaold=delta;
        deltadashold=deltadash;
        //aold = a;
    
        // current z, a values
        a=astart+ia*da;
        double z = 1/a-1;
        ztmp.push_back(z);
    
        double deltaddash1=-(3.+0.5*dlnEdlna(z))*(deltadashold/a);
        double deltaddash2=(3./2.)*OmegaM_/(a*a*a*Ez(z))*(deltaold/(a*a));
        double deltaddash = deltaddash1+deltaddash2;
						
        dddashtmp1.push_back(deltaddash1);
        dddashtmp2.push_back(deltaddash2);

        deltadash=deltadashold+da*deltaddash;
        delta=deltaold + da*deltadash;
        deltmp.push_back(delta);
        }
    // ************* WARNING *****************//
    // This is not the correct growth for w0!=-1 and wa!=0
    //GrowthZ0(); // for now going to cheat with polyfit to WLCode

    // get value of delta at z=0
    dz0_=deltmp[na-1];
    //cout << "OmegaM = "<<OmegaM_<<", w0 = "<<w0_<<", delta at z=0, "<<dz0_<<", na="<<na<<endl;
    //cout << endl;

    // reverse the order, and divide so normalised to 1 at z=0
    for (ia=0;ia<na;ia++) {
        zvals_.push_back(ztmp[na-(ia+1)]);
        deltavals_.push_back(deltmp[na-(ia+1)]/dz0_);
        dddash1_.push_back(dddashtmp1[na-(ia+1)]);
        dddash2_.push_back(dddashtmp2[na-(ia+1)]);
        }

    // make the interp table
    DEgrofunc_.DefinePoints(zvals_,deltavals_,zvals_[0],zvals_[zvals_.size()-1],2*na);
    CalcDEDone_=true;

};


double GrowthFN::dEdlna(double z)
{
    // valid for dark energy e.o.s parameterised as: w=w0+(1-a)wa
    // and non-flat universes

    double aa=1/(1+z);

    double p1 = -3.0*(1.0+w0_+wa_);
    double p2 = -3.0*(1.0-aa)*wa_;
    double p3 = -3.0*((2./3.) + w0_ + wa_);
      
    double dEEdlna = -3.0*OmegaM_/(aa*aa*aa) -2.0*OmegaK_/(aa*aa)
                     + p1*OmegaL_*(pow(aa,p1))*exp(p2)
                     + 3.0*wa_*OmegaL_*(pow(aa,p3))*exp(p2);

    return dEEdlna;

};


//********************** PkSpecCalc ******************************************//

double PkSpecCalc::operator() (double k,double z) const
{
    double tf = tf_(k);
    double pkinf = pkinf_(k);
    double d1 = d1_(z);

    double v = pkinf * (tf*tf) * (d1*d1);
    if(typspec_ == DELTA) 
        v *= k*k*k/(2.*M_PI*M_PI);

    return scale_ * v;
};


//********************** PkTabulate ******************************************//

void PkTabulate ::SetInterpTyp(int typ)
// see comment in InterpTab
{
    if(typ<0) typ=0; else if(typ>2) typ=2;
    interptyp_ = typ;
};


double PkTabulate::operator() (double k) const
{
    double v = InterpTab(k,k_,pk_,interptyp_);
    if(typspec_ == DELTA) 
        v *= k*k*k/(2.*M_PI*M_PI);
    return scale_ * v;
};


double PkTabulate::operator() (double k, double z) const
{
    cout<<"PkTabulate::operator(double k,double z)_Error: not implemented"<<endl;
    throw AllocationError("PkTabulate::operator(double k,double z)_Error: not implemented");
}

void PkTabulate::SetZ(double z)
{
    if(d1_ == NULL) {
        string emsg = "PkTabulate::SetZ_Error: d1==NULL, no possible redshift change";
        emsg += " for tabulated Pk";
        throw ParmError(emsg);
        }
    if(pk_.size() == 0) {
        string emsg = "PkTabulate::SetZ_Error: pk_.size()==0, no possible redshift";
        emsg += " change for tabulated Pk";
        throw ParmError(emsg);
        }

    double zold = zref_;
    if(fabs(z-zold)<1.e-4) return;

    zref_ = z;
    double d0 = (*d1_)(zold);
    double d1 = (*d1_)(zref_);
    double conv = d1*d1 / (d0*d0);
    cout <<"PkTabulate::SetZ: change redshift from "<< zold <<" (d="<< d0
         <<") to "<< zref_ <<" (d="<< d1 <<") conv="<< conv <<endl;
    
    for(unsigned int i=0;i<pk_.size();i++) 
        pk_[i] *= conv;
};


int PkTabulate::ReadCAMB(string filename, double h100tab, double zreftab)
{
    FILE *file = fopen(filename.c_str(),"r");
    if(file==NULL) return -1;
    cout <<"PkTabulate::ReadCAMB: fn="<< filename;
    cout <<" h100="<< h100tab <<" zreftab = "<< zreftab <<endl;

    const int lenline = 512;
    char *line = new char[lenline];
    double h = h100tab, h3 = pow(h100tab,3.);

    k_.resize(0); pk_.resize(0);
    double kmax = 0., pkmax = 0.;
    while ( fgets(line,lenline,file) != NULL ) {
   
        double k, pk;
        sscanf(line,"%lf %lf",&k,&pk);
        k *= h;      // convert h    Mpc^-1  ->  Mpc^-1
        pk /= h3;    // convert h^-3 Mpc^3   ->  Mpc^3
        
        if(pk>pkmax) { pkmax = pk; kmax = k; }
        
        k_.push_back(k);
        pk_.push_back(pk);
        }

    zref_ = zreftab;
    cout <<"  nread="<< pk_.size() <<" zref="<< GetZ();
    cout <<" , k,pk: max="<< kmax <<","<< pkmax;
    if(pk_.size()>0) {
        cout<<" [0]="<< k_[0] <<","<< pk_[0]
            <<" [n]="<<k_[pk_.size()-1]<<","<<pk_[pk_.size()-1];
        }
    cout<<endl;

    delete [] line;

    return (int)pk_.size();
};


//********************* PkEH *************************************************//

// P(k) calculated from E&H 1998 approximations

double PkEH::operator() (double k, double z) const
{
    double tf = tf_(k);
    double pkinf = pkinf_(k);
    double d1 = d1_(z);

    double v = pkinf * (tf*tf) * (d1*d1);
    if(typspec_ == DELTA) 
        v *= k*k*k/(2.*M_PI*M_PI);

    return scale_ * v;
};


//********************* PkAltNorm ********************************************//

// P(k) calculated from E&H 1998 approximations AND different normalisation
PkAltNorm::PkAltNorm(TransferEH& tf,GrowthFN& d1,double zref, double ns, double DelRsq)
  : tf_(tf) , d1_(d1) , ns_(ns) , DelRsq_(DelRsq)
{
    zref_ = zref;
    kpivot_ = 0.05;
    Apiv_ = 8.0967605e7;
};


PkAltNorm::PkAltNorm(PkAltNorm& pkz)
  : tf_(pkz.tf_) , d1_(pkz.d1_)
{

    ns_=pkz.ns_;
    DelRsq_=pkz.DelRsq_;
    kpivot_=pkz.kpivot_;
    Apiv_=pkz.Apiv_;

};


double PkAltNorm::operator() (double k, double z) const
{

    double OmegaM=tf_.ReturnOmegaM();
    double h = tf_.Returnh100();

    // transfer function
	double tf = tf_(k);

    // ns+3
    double pwr  = Power();

    // k in units of Mpc/h / kpivot (??)
    double knorm = kNorm(k,h);

    // hubble para^4
    double h4 = hpower4(h);

    //  growth function, see Carroll,Press,Turner 92
    double d1 = Gro(z);

    // growth function^2/OmegaM^2
    double powfactor = PowFactor(OmegaM);

    double DELTAK2 = Apiv_*pow(knorm,pwr)*tf*tf*powfactor*DelRsq_*d1*d1/h4;

    double v;
    if(typspec_ == DELTA)
        v = DELTAK2;
    else
        { v = (2.*M_PI*M_PI)*DELTAK2/(k*k*k); v*=(h*h*h)/pow(h,pwr); }
    // need to multiply the above by h^3 because k is in units
    // of Mpc^-1, not h/Mpc, another factor of h is also cancelled

    return v;
};


//******************* VarianceSpectrum ***************************************//

void VarianceSpectrum::SetRadius(double R)
// R = size of filter, top-hat or Gaussian
{
    if(R<=0.) {
        cout<<"VarianceSpectrum::SetRadius_Error: R<=0"<<endl;
        throw ParmError("VarianceSpectrum::SetRadius_Error: R<=0");
        }
    R_ = R;
};


void VarianceSpectrum::SetH100(double h100)
{
    if(h100<=0.) {
        cout<<"VarianceSpectrum::SetH100_Error: h100<=0"<<endl;
        throw ParmError("VarianceSpectrum::SetH100_Error: h100<=0");
        }
    h100_ = h100;
};


void VarianceSpectrum::SetInteg(double dperc, double dlogkinc, double dlogkmax, 
                                                         unsigned short glorder)
// WARNING: we don't integrate f(k)*dk but k*f(k)*d(log10(k))
// see argument details in function IntegrateFuncLog (geneutils.cc)
{
    dperc_ = dperc;  if(dperc_<=0.) dperc_ = 0.1;
    dlogkinc_ = dlogkinc;
    dlogkmax_ = dlogkmax;
    glorder_ = glorder;
};


double VarianceSpectrum::Filter2(double x) const
// WANRING: the filter-squared is returned
{
    // Just integrate the spectrum without filtering
    if(typfilter_ == NOFILTER) return 1.;

    double x2 = x*x;
    // Gaussian filter G(x) = exp(-x^2/2)
    //          note G(x)^2 = exp(-x^2)
    // Take the DL of G(x)^2 for x->0 to the order O(x^6)
    //             DL(x) = 1-x^2*(1-x^2/2)
    //             for x<0.01  |DL(x)-G(X)^2|<2.0e-13
    if(typfilter_ == GAUSSIAN)
        { if(x<0.01) return 1.-x2*(1.-x2/2.); else return exp(-x2); }

    // Top-hat filter T(x) = 3*(sin(x)-x*cos(x))/x^3
    // --- Take care of the pseudo-divergence for x->0
    // Take the DL of T(x)^2 for x->0 to the order O(x^7)
    //             DL(x) = 1-x^2/5*(1-3*x^2/35*(1-4*x^2/81))
    //             for x<0.1  |DL(x)-T(X)^2|<2.5e-13
    double f2=0.;
    if(x<0.1)
        f2 = 1.-x2/5.*(1.-3.*x2/35.*(1.-4.*x2/81.));
    else {
        f2 = 3.*(sin(x)-x*cos(x))/(x2*x);
        f2 *= f2;
        }
    return f2;
  
};


double VarianceSpectrum::Variance(double kmin, double kmax)
// Compute variance of spectrum pk_ by integration
// Input:
//     kmin,kmax = integration limits of k for calculating the variance
// Return:
//     value of the variance (sigma^2)
// Note:
//   The mass variance is returned
{
    if(kmin<=0 || kmax<=0. || kmin>=kmax) {
        cout<<"VarianceSpectrum::Variance_Error: kmin<=0 or kmax<=0 or kmin>=kmax"<<endl;
        throw ParmError("VarianceSpectrum::Variance_Error: kmin<=0 or kmax<=0 or kmin>=kmax");
        }

    double lkmin = log10(kmin), lkmax = log10(kmax);
    double var = IntegrateFuncLog(*this,lkmin,lkmax,dperc_,dlogkinc_,dlogkmax_,glorder_);

    return var;
};


double VarianceSpectrum::FindMaximum(double kmin, double kmax, double eps)
// Return the maximum of the function to be integrated
// Search between [kmin,kmax] but not logarithmically
// Input:
//     kmin,kmax : search interval
//     eps : precision required on the values
// Return:
//     position (in k) of the maximum
{
    if(kmin<=0 || kmax<=0. || kmin>=kmax) {
        string emsg = "VarianceSpectrum::FindMaximum_Error: kmin<=0 or kmax<=0";
        emsg +=" or kmin>=kmax || eps<=0";
        throw ParmError(emsg);
        }

    int n = 10; // always >2
    double lkmin = log10(kmin), lkmax = log10(kmax), dlk = (lkmax-lkmin)/n;

    double lkfind=lkmin, pkfind=-1.;
    while(1) {
        for(int i=0; i<=n; i++) {
            double lk = lkmin  + i*dlk;
            double v = (*this)(pow(10.,lk));
            if(v<pkfind) continue;
            pkfind = v; lkfind = lk;
            }
    
        //cout<<"VarianceSpectrum::FindMaximum: lkfind="<<lkfind<<" pkfind="<<pkfind
        //    <<" lkmin,max="<<lkmin<<","<<lkmax<<" dlk="<<dlk<<endl;
        
        // --- Convergence if "kfind" is such that "dk/kfind<eps"
        // Have dk = 10^(lkfind+dlk) - 10^(lkfind-dlk) = kfind * (10^(dlk) - 10^(-dlk))
        if( pow(10.,dlk)-pow(10.,-dlk) < eps ) 
            break;
    
        if(lkfind-dlk>lkmin) lkmin = lkfind-dlk;
        if(lkfind+dlk<lkmax) lkmax = lkfind+dlk;
        dlk = (lkmax-lkmin)/n;
        }

    return pow(10.,lkfind);
};


int VarianceSpectrum::FindLimits(double high, double &kmin, double &kmax, double eps)
// Return "[kmin,kmax]" such that the function to be integrated has "f(k) <= high"
// Search between [kmin,kmax] not logarithmically
// Input:
//     kmin,kmax : search interval
//     eps : precision required on the values kmin and kmax
// Output:
//     kmin,kmax such that "f(k) <= high"
// Return:
//     rc  = 0 if OK
//     rc |= 1 "f(kmin) >= high"   (bit0 =1)
//     rc |= 2 "f(kmax) >= high"   (bit1 =1)
//     rc |= 4 "f(k) < high for all k"   (bit2 =1)
{
    if(kmin<=0 || kmax<=0. || kmin>=kmax  || eps<=0.) {
        string emsg = "VarianceSpectrum::FindLimits_Error: kmin<=0 or kmax<=0";
        emsg +=" or kmin>=kmax or eps<=0";
        throw ParmError(emsg);
        }

    int n = 10; // always >2

    int rc = 0;
    double lkmin,lkmax,dlk,lkfind;

    // --- Find kmin
    lkmin=log10(kmin); lkmax=log10(kmax); dlk=(lkmax-lkmin)/n;
    while(1) {
        
        lkfind = lkmin;
        for(int i=0;i<=n;i++) {
            if( (*this)(pow(10,lkfind)) >= high ) 
                break;
            lkfind = lkmin + i*dlk;
            }
        //cout<<"VarianceSpectrum::FindLimits[kmin]: lkfind="<<lkfind
        //    <<" lkmin,max="<<lkmin<<","<<lkmax<<" dlk="<<dlk<<endl;
    
        // protect against f(k)<high for all k
        if(fabs(lkfind-lkmax)<dlk/2.) 
            {rc |= 4; return rc;} 
            
        if( pow(10.,dlk)-pow(10.,-dlk) < eps ) 
            break;
        
        if(lkfind-dlk>lkmin) lkmin = lkfind-dlk;
        if(lkfind+dlk<lkmax) lkmax = lkfind+dlk;
        dlk = (lkmax-lkmin)/n;
        }
        
    if(lkfind-lkmin<dlk/2.) rc |= 1;  // f(kmin) >= high
    else kmin = pow(10.,lkmin);
    //cout<<"rc="<<rc<<" lkmin="<<lkmin<<"  pk="<<(*this)(pow(10.,lkmin))<<endl;

    // --- Find kmax
    lkmin=log10(kmin); lkmax=log10(kmax); dlk=(lkmax-lkmin)/n;
    while(1) {
    
        lkfind=lkmax;
        for(int i=0;i<=n;i++) {
            if( (*this)(pow(10,lkfind)) >= high ) 
                break;
            lkfind -= dlk;
            lkfind = lkmax - i*dlk;
            }
        //cout<<"VarianceSpectrum::FindLimits[kmax]: lkfind="<<lkfind
        //    <<" lkmin,max="<<lkmin<<","<<lkmax<<" dlk="<<dlk<<endl;
    
        if( pow(10.,dlk)-pow(10.,-dlk) < eps ) 
            break;
    
        if(lkfind-dlk>lkmin) lkmin = lkfind-dlk;
        if(lkfind+dlk<lkmax) lkmax = lkfind+dlk;
        dlk = (lkmax-lkmin)/n;
        }
  
    if(lkmax-lkfind<dlk/2.) rc |= 2;  // f(kmax) >= high
    else kmax = pow(10.,lkmax);
    //cout<<"rc="<<rc<<" lkmax="<<lkmax<<"  pk="<<(*this)(pow(10.,lkmax))<<endl;

    return rc;
};

}  // end of namespace SOPHYA
