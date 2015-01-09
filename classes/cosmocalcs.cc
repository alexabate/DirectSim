#include "cosmocalcs.h"

/******* SimpleUniverse constructors ******************************************/

/* Default parameters */
SimpleUniverse::SimpleUniverse()
{
    Initialize();
    h_ = 0.719;
    omegarad_ = PhotonNuDensityKgm3() / CriticalDensityKgm3();
    omegaphot_ = PhotonDensityKgm3() / CriticalDensityKgm3();;
    omegabaryon_ = 0.0441;
    omegamat_ = 0.214+omegabaryon_;
    // We make it a flat universe using OmegaLambda
    SetFlatUniverse_OmegaLambda();
}

/* Universe only with Omega_M */
SimpleUniverse::SimpleUniverse(double h, double omega0)
{
    Initialize();
    h_ = h;
    omegarad_ = PhotonNuDensityKgm3() / CriticalDensityKgm3();
    omegaphot_ = PhotonDensityKgm3() / CriticalDensityKgm3();;
    omegabaryon_ = 0.044;
    omegamat_ = omega0;
    ComputeCurvature();
}

/* Universe with Omega_M and Omega_L */
SimpleUniverse::SimpleUniverse(double h, double omega0, double omegaL)
{
    Initialize();
    h_ = h;
    omegarad_ = PhotonNuDensityKgm3() / CriticalDensityKgm3();
    omegaphot_ = PhotonDensityKgm3() / CriticalDensityKgm3();;
    omegabaryon_ = 0.044;
    omegamat_ = omega0;
    if (omegaL > 1.e-79) {
        omegaL_ = omegaL;
        hasL_ = true;
        }
    ComputeCurvature();
}

/******* SimpleUniverse methods ***********************************************/


/* --Method-- */
void SimpleUniverse::SetOmegaMatter(double omega0)
{
    omegamat_ = omega0;
    if (omegamat_ < omegabaryon_) {
        cout << " SimpleUniverse::SetOmegaMatter()/Warning  Setting " 
	    << " OmegaBaryon=OmegaMatter=" <<  omega0 << endl;
        omegabaryon_ = omegamat_;
        }
    ComputeCurvature();
}

/* --Method-- */
void SimpleUniverse::SetOmegaLambda(double omegaL)
{
    if (omegaL > 1.e-39) {
        omegaL_ = omegaL;
        hasL_ = true;
        hasX_ = false;
        }
    else {
        omegaL_ = 0.;
        hasL_ = false;
        }
    ComputeCurvature();
}

/* --Method-- */
void SimpleUniverse::SetOmegaRadiation(double omegarad)
{
    omegarad_ = omegarad;
    if (omegarad_ < omegaphot_) {
        cout << " SimpleUniverse::SetOmegaRadiation()/Warning  Setting " 
	    << " OmegaPhoton=OmegaRad=" <<  omegarad << endl;
        omegaphot_ = omegarad_;
        }
    ComputeCurvature();
}

/* --Method-- */
void SimpleUniverse::SetOmegaBaryon(double omegabar)
{
    omegabaryon_ = omegabar;
    if (omegabaryon_ > omegamat_) {
        cout << " SimpleUniverse::SetOmegaBaryon()/Warning  Setting " 
	    << " OmegaMatter=OmegaBaryon=" <<  omegabar << endl;
        omegamat_ = omegabaryon_;
        ComputeCurvature();
        }
}

/* --Method-- */
void SimpleUniverse::SetOmegaPhoton(double omegaphot)
{
    omegaphot_ = omegaphot;
    if (omegaphot_ > omegarad_) {
        cout << " SimpleUniverse::SetOmegaPhoton()/Warning  Setting " 
	    << " OmegaRad=OmegaPhoton=" <<  omegaphot << endl;
        omegarad_ = omegaphot_;
        ComputeCurvature();
        }
  
}

/* --Method-- */
// for some crazy reason omegaL_ is *still* a quantity even 
// if dark energy is set.  Added omegaL_=0 to the below
// suppose you could have *both* dark energy AND cosmological constant
// but it's too upsetting to consider
void SimpleUniverse::SetDarkEnergy(double omegaX, double wX)
{
    if (omegaX > 1.e-39) {
        omegaX_ = omegaX;
        wX_ = wX;
        omegaL_ = 0;
        hasX_ = true;
        hasL_ = false;
        }
    else {
        omegaX_ = 0.;
        wX_ = -1.;
        omegaL_ = 0;
        hasX_ = false;
        hasL_ = false;
        }
    ComputeCurvature();
}

/* --Method-- */
void SimpleUniverse::SetDarkEnergyEoS(double wX)
{

    if (omegaX_ > 0.) {
        wX_ = wX;
        hasX_ = true;
        hasL_ = false;
        }
    else {
        cout <<" SimpleUniverse::SetDarkEnergyEoS() Warning: setting dark";
        cout <<" energy density using cosmological constant density "<< endl;
        omegaX_ = omegaL_;
        omegaL_ = 0;
        wX_ = wX;
        hasX_ = true;
        hasL_ = false;
        }
    ComputeCurvature();
}

/* --Method-- */
void SimpleUniverse::SetFlatUniverse_OmegaMatter()
{
    omegamat_ = 1. - (omegarad_ + omegaL_ + omegaX_);
    if (omegamat_ < 0.) {
        cerr << " SimpleUniverse::SetFlatUniverse_OmegaMatter()/Error -> " 
	    << " OmegaMatter=" << omegamat_ << " < 0" << endl;
        omegamat_ = 0.;
        throw ParmError("SimpleUniverse::SetFlatUniverse_OmegaMatter() OmegaMatter<0");
        }
    ComputeCurvature();
}

/* --Method-- */
void SimpleUniverse::SetFlatUniverse_OmegaLambda()
{
    omegaL_ = 1. - (omegamat_ + omegarad_ + omegaX_);
    if (omegaL_ < 0.) {
        cerr << " SimpleUniverse::SetFlatUniverse_OmegaLambda()/Error -> " 
	    << " OmegaL=" << omegaL_ << " < 0" << endl;
        omegaL_ = 0.;
        throw ParmError("SimpleUniverse::SetFlatUniverse_OmegaLambda() OmegaL<0");
        }
    ComputeCurvature();
}


/* --Method-- */
void SimpleUniverse::Initialize()
{
    h_ = 0.5;
    omegamat_ = 1.;
    omegarad_ = 0;
    omegabaryon_ = 0;
    omegaphot_ = 0;
    omegaL_ = 0;
    omegaX_ = 0.;
    wX_ = -1.;
    hasL_ = hasX_ = false;

    // We initilalize ze_, ...  to today
    ze_ = 0.;
    te_ = 0.;
    chie_ = 0.;
    da_ = 0.;
    dl_ = 0.;
    EofZe_ = 1.;
    integGz_ = 0;
    integGTz_ = 0;

    prt_dbg_level = 0;
}

/* --Method-- */
void SimpleUniverse::ComputeCurvature()
{

    omegaCurv_ = (1. - (omegamat_ + omegarad_ + omegaL_ + omegaX_));
    if (fabs(omegaCurv_) < 1.e-39) {
        omegaCurv_ = 0.;
        kcurvature_ = 0;
        }
    else { 
        //    if (omegaCurv_ < 0) { kcurvature_ = 1;  omegaCurv_ = -omegaCurv_; }
        if (omegaCurv_ < 0) kcurvature_ = 1;  
        else kcurvature_ = -1;
        }

    // We initilalize ze_, ...  to today
    ze_ = 0.;
    te_ = 0.;
    chie_ = 0.;
    da_ = 0.;
    dl_ = 0.;
    EofZe_ = 1.;
}


/* --Method-- */
double SimpleUniverse::OmegaMatterZE()  const
{
    double zz = (1+ze_);
    return (OmegaMatter()*zz*zz*zz/(EofZe_*EofZe_));
}

/* --Method-- */
double SimpleUniverse::OmegaRadiationZE() const
{
    double zz = (1+ze_);
    return (OmegaRadiation()*zz*zz*zz*zz/(EofZe_*EofZe_));
}

/* --Method-- */
double SimpleUniverse::OmegaCurvZE() const
{
    double zz = (1+ze_);
    return (OmegaCurv()*zz*zz/(EofZe_*EofZe_));
}

/* --Method-- */
double SimpleUniverse::OmegaMatter(double z)  const
{
    double zz = (1+z);
    return (OmegaMatter()*zz*zz*zz/(Ez(z)*Ez(z)));
}

/* --Method-- */
double SimpleUniverse::OmegaRadiation(double z) const
{
    double zz = (1+z);
    return (OmegaRadiation()*zz*zz*zz*zz/(Ez(z)*Ez(z)));
}

/* --Method-- */
double SimpleUniverse::PhotonDensitycm3(double z) const
{
// s/k = g 2 pi^2 / 45 ( k T / hbar c )^3 = 3.602 ngamma 
// g = 2 for photons, ngamma = nbgamma / m^3
// Obtain however ngamma = 410 photons / cm^3 = 4.1 e8 photons / m^3
    if (z>1.e-38) return PhotonNumberDensitycm3()*(1+z)*(1+z)*(1+z) ;
    else return PhotonNumberDensitycm3();
}

/* --Method-- */
double SimpleUniverse::BaryonDensitycm3(double z) const
{
    double nbar = OmegaBaryon()*CriticalDensityGeVcm3()/PROTON_MASS_IN_GEV;
//  cout << "\n BaryonDensitycm3-DEBUG Ob= " << OmegaBaryon() 
//       << " CritDens=" << CriticalDensityGeVcm3()
//       << " pmassGeV= " << ProtonMassGeV_Cst 
//       << " nbar= " << nbar << endl; 
    if (z>1.e-38) return nbar*(1+z)*(1+z)*(1+z) ;
    else return nbar;
}

/* --Method-- */
double SimpleUniverse::EtaBaryonPhoton() const
{
    double nbar = OmegaBaryon()*CriticalDensityGeVcm3()/PROTON_MASS_IN_GEV;
    return ( nbar / PhotonNumberDensitycm3() ); 
}


/******* --MAIN Method-- ******************************************************/
void SimpleUniverse::SetEmissionRedShift(double ze, double prec, bool fginc)
{
    if (ze < 0.) { // z is unphysical
        cout << " SimpleUniverse::SetEmissionRedShift/Error ze = " << ze 
	    << " less than 0 ! " << endl;
        throw ParmError("SimpleUniverse::SetEmissionRedShift(ze < 0.)");
        }

    if (prt_dbg_level > 1) // if printing
        Timer tm("SimpleUniverse::SetEmissionRedShift(...)");

    if (ze < 1.e-39) { // if ~today
        ze_ = 0.;
        te_ = 0.;
        chie_ = 0.;
        da_ = 0.;
        dl_ = 0.;
        EofZe_ = 1.;
        integGz_ = 0;
        integGTz_ = 0;
        return;
        }

    double zelast = ze_; // set previous emission redshift to be last 
    ze_ = ze;            // set current emission redshift
    EofZe_ = Ez(ze);

    // Reference: Principles of Physical Cosmology - P.J.E. Peebles  
    //            Princeton University Press - 1993
    //              ( See Chapter 13)
    // We have to integrate Integral(dz / E(z)) from 0 to ze  (cf 13.29)
    //      E(z) = Sqrt(Omega0*(1+z)^3 + OmegaCurv*(1+z)^2 + OmegaL) (cf 13.3)
    // G(z) = 1/E(z) : integration element for the calculation of D_A  
    // GT(z) = 1/[(1+z)*E(z)] : integration element for the calculation of age/time (dt = GTz(z)*dz)

    double idzoez = 0.;
    double idzToez = 0.;
    if (fginc) {  // Incremental calculation
        double inGz = 0.;
        double inGTz = 0.;
        
        if (zelast < ze_) {
            // Numerical integration of i) inGz: G(z)dz = 1/E(z) dz between [zelast,ze], precision prec
            // Numerical integration of ii) inGTz: GT(z)dz = 1/[(1+z)*E(z)] dz between [zelast,ze], precision prec
            NumIntegrateGzGTz(zelast, ze_, prec, inGz, inGTz);
            integGz_ += inGz;
            integGTz_ += inGTz;
            }
        else {
            NumIntegrateGzGTz(ze_, zelast, prec, inGz, inGTz);
            integGz_ -= inGz;
            integGTz_ -= inGTz;
            }
        }
    else NumIntegrateGzGTz(0., ze_, prec, integGz_, integGTz_);

    idzoez = integGz_;
    idzToez = integGTz_;
    if ((prt_dbg_level > 1) && !hasL_ && !hasX_ && (kcurvature_ != 0) ) {
        double idzoez_ana = IntegrateGz_NoLambda(ze_) - IntegrateGz_NoLambda(0.);
        cout << " SimpleUniverse::SetEmissionRedShift(...)/Check idzoez= " 
	    << idzoez << "  idzoez_analyt =" << idzoez_ana 
	    << " RelDiff= " << 100.*(idzoez-idzoez_ana)/idzoez_ana << " %" << endl;
        }
    // 1/(H0a0R) = Sqrt(OmegaCurv)   (cf 13.4)

    dcl_=idzoez;
    // Notation of Hogg 1999:
    // Hubble distance:
    // dh = c/H0
    // Line of sight comoving distance:
    // dc = dh int_0^z dz/E(z')
    // Transverse comoving distance:
    // i)     OPEN UNI dm = dh 1/sqrt(OmK)sinh[sqrt(OmK) dc/dh]
    // ii)    FLAT UNI dm = dc
    // iii) CLOSED UNI dm = dh 1/sqrt(OmK)sin[sqrt(|OmK|) dc/dh]
    if ( IsFlat() ) {
    
	    // here chie_= int_z1^z2 dz/E(z')
	    // therefore it is equal to dc/dh=dm/dh
	    dct_ = chie_ = da_ = idzoez; } 
    else {
	    double sov;
  	    if (IsOpen())  
		    sov = sqrt(OmegaCurv());
    	    else 
		    sov = sqrt(-OmegaCurv());

	    // here chie_= sqrt(OmK) int_z1^z2 dz/E(z') = sqrt(OmK) dc/dh
	    // **and is NOT equal to dc/dh or dm/dh**
    	    chie_ = sov * idzoez;
    	    if (IsOpen()) 
		    da_ = sinh(chie_);
    	    else  
		    da_ = sin(chie_);
    	    da_ /= sov;
	    dct_=da_;
        }
        
    // September 2003 : Correction of a bug with da_ and dl_  (Reza)
    // The value of da_ for the angular diameter distance should be taken
    // at the emission redshift (ze) so that delta_theta = Proper_distance / Distance_da_
    // Then must divide da_ calculated below by (1+ze_) to get
    // to the distance corresponding to ze_
    // See Peebles formulae (13.46 - 13.47) 
    da_ /= (1+ze_);

    // There is a factor (1+ze_)^2 between da_ and dl_
    // > (1+ze_) to get da_ to the value z=0 (today) (1)
    // > (1+ze_) to compensate for the redshift of the photons (2)
    // See for example Weinberg, Gravitation and Cosmology, page 423 (14.4.22)
    dl_ = da_*(1+ze_)*(1+ze_);   

    // The look back time
    te_ = idzToez;

    // volume element, see Hogg 1999, added by AA June 2011
    //  dvdzdo_ = da_*da_*(1+ze_)*(1+ze_);
}

/* --Method-- */
double SimpleUniverse::DeltaTime(double z1, double z2, double prec)
{
    double xx, deltatime;
    if (z1 > z2)  
        NumIntegrateGzGTz(z2, z1, prec, xx, deltatime);
    else
        NumIntegrateGzGTz(z1, z2, prec, xx, deltatime);
    return deltatime;
}

/* --Method-- */
double SimpleUniverse::DeltaConformalTime(double z1, double z2, double prec)
{
    double xx, deltaeta;
    if (z1 > z2)  
        NumIntegrateGzGTz(z2, z1, prec, deltaeta, xx);
    else
        NumIntegrateGzGTz(z1, z2, prec, deltaeta, xx);
    return deltaeta;
}

/* --Method-- */
void SimpleUniverse::Print(ostream & os) const 
{
    os  << "--------------------------------------------------------------------------" 
        << endl;

    os  << " Time=0 (now)" << " OmegaMatter= " << OmegaMatter()
        << " OmegaRad= " << OmegaRadiation()  << endl;
    os  << " OmegaBaryon= " << OmegaBaryon() << " OmegaCDM= " << OmegaCDM() 
        << " OmegaPhoton= " << OmegaPhoton();
        
    if (OmegaPhoton() > 1.e-39) os << " Zequ= " << Z_Equality() << endl;
    else os << " Zequ= ???? " << endl;

    if (hasX_) os << " OmegaDE= " << OmegaDE();
    else os << " OmegaL= " << OmegaLambda();
    os << " OmegaCurv=" << OmegaCurv() << " KCurv=" << KCurvature();

    if (KCurvature() == 0) os << " Flat " << endl;
    else if (KCurvature() < 0) os << " Open " << endl;
    else  os << " Closed " << endl;

    os  << " H0= " << H0() << " km/sec/Mpc (h=" << h() 
        << ") HubbleTimeYear= " << HubbleTimeYear()
        << " HubbleLengthMpc= " << HubbleLengthMpc() << endl;
    os  << " >> Emission RedShift= ze= " << ZE() << " ChiE= " << chie_
        << " H(ze)= " << HZE() << " km/s/Mpc " << endl;
    os  << " >> @ ze: Time=" << LookBackTime() 
        << " (HubbleTimeUnit)  TimeYear=" << LookBackTimeYear()  << endl;
    os  << " >> OmegaMatterZE= " << OmegaMatterZE()
        << " OmegaRadiationZE= " << OmegaRadiationZE()  
        << " OmegaCurvZE=" << OmegaCurvZE() << endl; 
  //  os << " >> AngularDiameterDistance= " << AngularDiameterDistance()
  //   << " LuminosityDistance= " << LuminosityDistance() << endl;
    os  << " >> AngDiamDistanceMpc= " << AngularDiameterDistanceMpc()
        << " (Comov=" << ComovAngDiamDistanceMpc() << ")" 
        << " LumDistanceMpc= " << LuminosityDistanceMpc() << endl;
    os  << "--------------------------------------------------------------------------" 
        << endl;
}



/* --Method-- */
double SimpleUniverse::Z_Equality() const 
{
    if (omegarad_ < 1.e-79) {
        cerr << "SimpleUniverse::Z_Equality()/Error OmegaRad = 0 !" << endl;
        throw ParmError("SimpleUniverse::Z_Equality()/Error OmegaRad = 0");
        }
    return(omegamat_/omegarad_-1.);
}



/* --Method-- */
double SimpleUniverse::Ez(double z) const
{
    double zz = 1+z;
    double ez2 = zz*zz* ( omegarad_*zz*zz + omegamat_*zz + omegaCurv_);
    ez2 += omegaL_;
    if (hasX_) ez2 += omegaX_*pow(zz,3.*(1.+wX_));
    return sqrt(ez2);
}

/* --Method-- added by AA June 2011 -- integrate Volume element dV/(dzdO) */
void SimpleUniverse::NumIntegVolElz(double z1, double z2, int nstep, 
				       double& dVdO)
// Simple trapesium integration with respect to z of
// (1+z)^2*DA^2/E(z): see Eqn 28 in Hogg 1999
// returned in hubble length units (not Mpc^3) therefore need to multiply
// by c/H0 to get Mpc^3 units
// multiply by solid angle X to get volume between z1 and z2 over solid angle X
{
 
    if ((z1 < 0.) || (z2 < 0.) || (z2 < z1)) {
        cout << " SimpleUniverse::NumIntegVolElz()/Error invalid values for z1,z2: " 
	         << z1 << "," << z2 << endl;
        throw ParmError("SimpleUniverse::NumIntegVolElz() invalid z1/z2 values");
        }

    double step = (z2-z1)/(nstep-1);
    double sum =0;
    for (int i=0;i<nstep;i++) {
    
	    double z = z1 + i*step;
	    SetEmissionRedShift(z);
	    sum+=(1+z)*(1+z)*AngularDiameterDistanceMpc()*AngularDiameterDistanceMpc()/Ez(z);
	    }

    dVdO=sum*step;

    return;
}

/* --Method-- */
void SimpleUniverse::NumIntegrateGzGTz(double z1, double z2, double prec, 
				       double& resG, double& resGT)
{
    NumIntegGzGTz(z1, z2, prec, resG, resGT);
    return;
}

  
/******* --MAIN CALCULATION method-- ******************************************/
void SimpleUniverse::NumIntegGzGTz(double z1, double z2, double prec, 
						double& resG, double& resGT)
{
    if ((z1 < 0.) || (z2 < 0.) || (z2 < z1)) {
        cout << " SimpleUniverse::NumIntegGzGTz()/Error invalid values for z1,z2: " 
	    << z1 << "," << z2 << endl;
        throw ParmError("SimpleUniverse::NumIntegGzGTz() invalid z1/z2 values");
        }

   /* -----  How the precision is used to calculate the step - 
   We use the relation :
      | Integ_{a,b} f(x) dx - (b-a)f(a) | <= (b-a)^2/2 Sup_{a,b} |f'(x)| 
   The function to integrate G(z) is a decreasing positive function
   We approximate Sup |f'(x)| par |f'(a)| = -f'(a) 
   I = Integ_{a,b} f(x) dx  >= (b-a)f(b)    , delta I / I <= prec
     ===> delta I/I <= (b-a)/2 | f'(a) / f(b) |  < prec 
   step = (b-a) ===> b-a <= 2*prec* | f(b)/f'(a) | 
   We will use : step = prec* | f(b)/f'(a) |
   */

    resG = resGT = 0.;

    double accG, accGT; 

    double z1orig = z1;
    double z2orig = z2;

    // Integration with a fixed step up to z=zbrk (=10.)
    double zbrk = 10.;
    // for z < zbrk, we take a=0 , G'(0) ~ 1 G(1) ~ 1 ===> step~prec 
    double cfac = fabs(DerivateGz(0.,1.e-3));
    if (cfac < 1.) cfac = 1.;
    double dz = prec*Gz(1)/cfac;
    if (dz > 0.1) dz = 0.1;
    double dzini = dz;
 
    //DBG  cout << " **DBG*A**  prec=" << prec << " zbrk=" << zbrk
    //DBG << " G(1.)= " << Gz(1.) << " dG(0)=" << DerivateGz(0.,dz) << " ==>dz " << dz << endl;

    double zc = z1;
    accG = accGT = 0.;
    if (z1<zbrk) {    // Integration up to zbrk 
        if (z2 > zbrk)  z2 = zbrk;
     
        while (zc < z2) {   // Calculate time and distance
            accG += Gz(zc);
            accGT += GTz(zc);
            zc += dz;   
            }
            
        zc -= dz;
        accG -= 0.5*(Gz(z1)+Gz(zc));
        accGT -= 0.5*(GTz(z1)+GTz(zc));
        accG *= dz;   accGT *= dz;
        dz = z2-zc;
        
        if (dz > 1.e-69) {       // Correction for the last step 
            accG += 0.5*dz*(Gz(zc)+Gz(z2));
            accGT += 0.5*dz*(GTz(zc)+GTz(z2));
            }
        }

    resG = accG;    resGT = accGT;

    if (z2orig <= zbrk)  return;
   
    z1 = zbrk;
    // if integration after z1 > zbrk 
    if (z1orig > zbrk)   { 
        z1 = z1orig;
        dz = prec*fabs(Gz(z1*4.)/DerivateGz(z1,0.1));
        //DBG    cout << " **DBG*C**  dz=" << dz << endl; 
        }
    else dz = dzini;

    z2 = z2orig;


    long kkMAX = 100;
    bool fgencore = true;
    double MAXRAPDZ = 10.;
     
    zc = z1;
    while ((zc < z2)&&fgencore) {

        double Delz = (double)kkMAX*dz;
        if (Delz < 1) Delz = 1.;

        // Calculation of step :
        double lastdz = dz;
        dz = prec*fabs(Gz(zc+Delz)/DerivateGz(zc,dz));
        //     double dzs = dz;

        if (dz/lastdz > MAXRAPDZ)  dz = lastdz*MAXRAPDZ;
        long kkmx = (long)((z2-zc)/dz)-1;
        if (kkmx > kkMAX) kkmx = kkMAX;

        //DBG cout << " **DBG*B**  prec=" << prec << " zc=" << zc << " Delz=" << Delz 
        //DBG  << " G(zc)= " << Gz(zc+Delz) << " dG=" << DerivateGz(zc,dz) << " ==>dz " 
        //DBG  << dz << " dzs=" << dzs << " kkmx=" << kkmx << endl; 

        double za = zc;
        accG = accGT = 0.;
        for (long kkp=0; kkp<kkmx; kkp++) {
            accG += Gz(zc);
            accGT += GTz(zc);
            zc += dz;   
            }
        zc -= dz;

        accG -= 0.5*(Gz(za)+Gz(zc));
        accGT -= 0.5*(GTz(za)+GTz(zc));
        accG *= dz;   accGT *= dz;
        resG += accG;    resGT += accGT;      
        if (kkmx < kkMAX) fgencore = false;
        }
   
    dz = z2-zc;
    if (dz > 1.e-69) {       // Correction for the last step 
        resG += 0.5*dz*(Gz(zc)+Gz(z2));
        resGT += 0.5*dz*(GTz(zc)+GTz(z2));
        }

    return;
}

/* 
1/ Integral[1/E(z) dz] is analytic if Lambda = 0
   ----> Result from mathematica :
   Eoz[z_] := 
              Sqrt[omr*(z + 1)^4 + omm*(z + 1)^3 + omk*(z + 1)^2]
   IEz = Integrate[1/Eoz[z], z]
   CForm[IEz]

   -(((1 + z)*Sqrt(omk + omm + omr + omm*z + 2*omr*z + 
   omr*Power(z,2))*
   Log((2*omk + omm + omm*z)/(Sqrt(omk)*(1 + z)) + 
   (2*Sqrt(omk + omm + omr + omm*z + 2*omr*z + 
   omr*Power(z,2)))/(1 + z)))/
   (Sqrt(omk)*Sqrt(omk*Power(1 + z,2) + 
   omm*Power(1 + z,3) + omr*Power(1 + z,4))))


   IETz = Integrate[1/((z + 1)Eoz[z]), z]
   CForm[IETz]

   -((omk + omm + omr + omm*z + 2*omr*z + omr*Power(z,2))/
      (omk*Sqrt(omk*Power(1 + z,2) + omm*Power(1 + z,3) + omr*Power(1 + z,4)))) + 
   (omm*(1 + z)*Sqrt(omk + omm + omr + omm*z + 2*omr*z + omr*Power(z,2))*
      Log((-2*(2*Power(omk,2) + omk*omm + omk*omm*z))/(Sqrt(omk)*omm*(1 + z)) - 
        (4*omk*Sqrt(omk + omm + omr + omm*z + 2*omr*z + omr*Power(z,2)))/(omm*(1 + z)))
      )/(2.*Power(omk,1.5)*Sqrt(omk*Power(1 + z,2) + omm*Power(1 + z,3) + 
        omr*Power(1 + z,4)))

   ---> Derivative of 1/E(z)  Lambda=0    
   dEoz = D[1/Eoz[z], z]
   CForm[dEoz] 
        -(2*omk*(1 + z) + 3*omm*Power(1 + z,2) + 4*omr*Power(1 + z,3))/
	(2.*Power(omk*Power(1 + z,2) + omm*Power(1 + z,3) + omr*Power(1 + z,4),1.5))

   ---> Derivative of 1/E(z) avec Lambda    
   Eozl[z_] := 
      Sqrt[omr*(z + 1)^4 + omm*(z + 1)^3 + omk*(z + 1)^2 + oml]
   dEozl = D[1/Eozl[z], z]
   CForm[dEozl]
      -(2*omk*(1 + z) + 3*omm*Power(1 + z,2) + 4*omr*Power(1 + z,3))/
      (2.*Power(oml + omk*Power(1 + z,2) + omm*Power(1 + z,3) + omr*Power(1 + z,4),1.5))

   ---> Derivative of 1/E(z) with Lambda and X (Dark energy)  
   EozX[z_] := 
              Sqrt[omr*(z + 1)^4 + omm*(z + 1)^3 + omk*(z + 1)^2 
	      + oml + omX*(1 + z)^(3(1 + w))]
   dEozX := D[1/EozX[z], z]

   CForm[dEozX]
         -(2*omk*(1 + z) + 3*omm*Power(1 + z,2) + 4*omr*Power(1 + z,3) + 
         3*omX*(1 + w)*Power(1 + z,-1 + 3*(1 + w)))/
	 (2.*Power(oml + omk*Power(1 + z,2) + omm*Power(1 + z,3) + omr*Power(1 + z,4) + 
	 omX*Power(1 + z,3*(1 + w)),1.5))

*/


inline double power_2(double x) { return x*x; }
inline double power_3(double x) { return x*x*x; }
inline double power_4(double x) { return x*x*x*x; }
inline double power_5(double x) { return x*x*x*x*x; }

/* --Method-- */
// Return the analytic integral Integ[1/E[z] dz] if OmegaLambda = 0
double SimpleUniverse::IntegrateGz_NoLambda(double z)
{
    return (
        -(((1 + z)*sqrt(omegaCurv_ + omegamat_ + omegarad_ + omegamat_*z + 2*omegarad_*z + 
        omegarad_*power_2(z))*
        log((2*omegaCurv_ + omegamat_ + omegamat_*z)/(sqrt(omegaCurv_)*(1 + z)) + 
        (2*sqrt(omegaCurv_ + omegamat_ + omegarad_ + omegamat_*z + 2*omegarad_*z + 
        omegarad_*power_2(z)))/(1 + z)))/
        (sqrt(omegaCurv_)*sqrt(omegaCurv_*power_2(1 + z) + 
        omegamat_*power_3(1 + z) + omegarad_*power_4(1 + z))))
        );
}

/* --Method-- */
// Return the analytic integral Integ[1/((1+z)E[z]) dz] if OmegaLambda = 0
double SimpleUniverse::IntegrateGTz_NoLambda(double z)
{
    return (
	    -((omegaCurv_ + omegamat_ + omegarad_ + omegamat_*z + 
	     2*omegarad_*z + omegarad_*power_2(z)) /
	      ( omegaCurv_*sqrt(omegaCurv_*power_2(1+z) + 
			      omegamat_*power_3(1+z) + omegarad_*power_4(1+z)) ) ) + 
	  (omegamat_*(1+z)*sqrt(omegaCurv_+omegamat_ + omegarad_+omegamat_*z + 
				2*omegarad_*z + omegarad_*power_2(z)) *
	   log((-2.*(2.*power_2(omegaCurv_) + omegaCurv_*omegamat_ + 
		     omegaCurv_*omegamat_*z)) / (sqrt(omegaCurv_)*omegamat_*(1+z)) - 
	       (4*omegaCurv_*sqrt(omegaCurv_+omegamat_+omegarad_+omegamat_*z + 
				  2*omegarad_*z + omegarad_*power_2(z))) / 
	       (omegamat_*(1+z)) )  ) / 
	  (2.*pow(omegaCurv_,1.5)*sqrt( omegaCurv_*power_2(1+z) + 
				        omegamat_*power_3(1+z) + 
				        omegarad_*power_4(1+z) ) )

        );  
}


/* --Method-- */
// Return the derivation d G(z) / dz
double SimpleUniverse::DerivateGz(double z, double dz)
{
    if (dz > 1.e-19) {  // Calculate a numerical derivative
        return((Gz(z+dz)-Gz(z))/dz);
        }
    else {   // Using the analytic form, problem ??
        if (hasX_) 
            return DerivateGz_Lambda_X(z);
        else { 
            if (hasL_) return DerivateGz_Lambda(z);
            else return DerivateGz_NoLambda(z);
            }
        }
}

/* --Method-- */
// Return the derivation d G(z) / dz Lambda = 0, OmegaX = 0
double SimpleUniverse::DerivateGz_NoLambda(double z)
{
    return ( -(2*omegaCurv_*(1+z) + 3*omegamat_*power_2(1+z) 
	     + 4*omegarad_*power_3(1+z)) / 
	   (2.*pow(omegaCurv_*power_2(1+z) 
		   + omegamat_*power_3(1+z) 
		   + omegarad_*power_4(1+z),1.5)) ) ;
}    

/* --Method-- */
// Return the derivation d G(z) / dz Lambda non zero
double SimpleUniverse::DerivateGz_Lambda(double z)
{
    return ( -(2.*omegaCurv_*(1+z) + 3.*omegamat_*power_2(1+z) + 
	     4.*omegarad_*power_3(1+z)) /
	   (2.*pow(omegaL_ + omegaCurv_*power_2(1+z) + 
		   omegamat_*power_3(1+z) + omegarad_*power_4(1+z), 1.5) ) );
}

/* --Method-- */
// Return the derivation d G(z) / dz Lambda and OmegaX (dark energy) non zero
double SimpleUniverse::DerivateGz_Lambda_X(double z)
{

    return ( -(2*omegaCurv_*(1+z) + 3*omegamat_*power_2(1+z) 
	     + 4*omegarad_*power_3(1+z) + 
	     3*omegaX_*(1+wX_)*pow(1+z, -1.+3.*(1.+wX_)) ) /
	   (2.*pow(omegaL_+omegaCurv_*power_2(1+z) + omegamat_*power_3(1+z) 
		   + omegarad_*power_4(1+z) + 
		   omegaX_*pow(1+z, 3.*(1.+wX_)),1.5) ) );

}
