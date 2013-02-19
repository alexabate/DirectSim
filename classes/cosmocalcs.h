/**
 * @file  cosmocalcs.h
 * @brief Contains methods that perform cosmological calculations
 *
 * Could add more information here I think
 *
 * @author Alex Abate, Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2008
 * @date 2008
 *
 */

#ifndef LUC_H_SEEN
#define LUC_H_SEEN

#include "machdefs.h"
#include <math.h>
#include <iostream>

// sophya
#include "sopnamsp.h"
#include "ctimer.h"
#include "pexceptions.h"

#include "constcosmo.h"


/** SimpleUniverse class
  * 
  * Class for calculating the different coordinates and parameters of the 
  * universe. Formulae and notation according to:
  * Principles of Physical Cosmology - P.J.E. Peebles - Chapter 13
  * 
  */
class SimpleUniverse {
public:
    // Constructor function of the cosmological parameters
    // h: Hubble constant in units of 100 km/s/Mpc (H0 = h x 100 km/s/Mpc)
    // omega0: Total matter density, in units of the critical density 
    // omegaL: Cosmological constant in units of the critical density

    // Default density of the photons  -> corresponds to T_CMB = 2.725 K
    // Default density of the neutrinos -> corresponds to T_nu ~= 1.7 K
    // Default density of the baryons ~ 0.044 * RhoCrit
    // Universe standard flat LambdaCDM WMAP 
  
    /** Default constructor */
    SimpleUniverse();
  
    /** Constructor for a universe with just one matter component determined
        by omega0. Omega_B is set to 0.044. */
    SimpleUniverse(double h, double omega0);
  
    /** Constructor for a universe composed of matter component (omega0) and 
        lambda (omegaL). Omega_B is set to 0.044. */
    SimpleUniverse(double h, double omega0, double omegaL);

    virtual ~SimpleUniverse();

    /** Set print debug level */
    inline void   SetPrtDbgLevel(int lev=0) { prt_dbg_level = lev; }

    /** Set the Hubble parameter h using H0 the Hubble constant */
    inline  void  SetH0(double H0) { Seth(H0/100.); }
    
    /** Set the Hubble parameter h */
    void  Seth(double h);

    // To define the density of the different components,
    // in units of critical density
    void  SetOmegaMatter(double omega0);        /**< Set Omega_M */
    void  SetOmegaLambda(double omegaL);        /**< Set Omega_Lambda */
    void  SetOmegaRadiation(double omegarad);   /**< Set Omega_Rad */
    void  SetOmegaBaryon(double omegabar);      /**< Set Omega_B */
    void  SetOmegaPhoton(double omegaphot);     /**< Set Omega_Photon */
    
    // Dark energy component, with equation of state
    // of type pressure = w rho  (rho = density)
    /** Set dark energy, Omega_DE and equation of state */
    void  SetDarkEnergy(double omegaX, double wX);
    
    /** Set dark energy equation of state */
    void  SetDarkEnergyEoS(double wX);
  
    /** Set flat universe using Omega_M */
    void  SetFlatUniverse_OmegaMatter();
    /** Set flat universe using Omega_Lambda */
    void  SetFlatUniverse_OmegaLambda();

    //---- Hubble constant and related quantities
    /** Return Hubble parameter (Hubble constant value normalized to 100 km/s/Mpc) */
    inline double h() const { return(h_); }
    /** Return Hubble constant value in km/s/Mpc */
    inline double H0() const { return(h_*100.); }
    /** Return Hubble constant value in s^-1 */
    inline double H0InvSec() const 
        { return(h_*100.*1.e3/MpctoMeters()); }
    /** Return Hubble constant value in h/s */
    inline double H0hperSec() const 
        { return(100.*1.e3/MpctoMeters()); }
    /** Return 1/Hubble constant in s */
    inline double HubbleTimeSec() const 
        { return(1./H0InvSec()); }
    /** Return 1/Hubble constant in years */
    inline double HubbleTimeYear() const 
        { return(HubbleTimeSec()/YeartoSec()); }
    /** Return Hubble length = 1/H0 x SpeedOfLight , in meters */
    inline double HubbleLengthMeter() const 
        { return(SpeedOfLight()*HubbleTimeSec()); }
    /** Return Hubble length = 1/H0 x SpeedOfLight , in Mpc */
    inline double HubbleLengthMpc() const 
        { return(HubbleLengthMeter()/MpctoMeters()); }
    /** Return Hubble length = 1/H0 x SpeedOfLight , in Mpc (Hogg 1999 notation) */
    // factor of 1000 to convert m/s to km/s
    inline double DH() const { return(SpeedOfLight()/(1000*H0())); }
    /** Return Hubble length = 1/H0 x SpeedOfLight , in Mpc/h (Hogg 1999 notation) */
    // Hubble length (Hogg 1999 notation) in h units
    inline double DHh() const { return(SpeedOfLight()/(1000*100)); }
    /** Expansion rate at emission redshift  H(tE) */
    inline double HZE() const { return(H0()*EofZe_); }


    // Matter and various energy densities today (z=0)
    inline double OmegaMatter() const { return(omegamat_); }/**< Return Omega_M */
    inline double OmegaRadiation() const { return(omegarad_); }/**< Return Omega_Rad */
    inline double OmegaPhoton() const { return(omegaphot_); }/**< Return Omega_Photon */
    inline double OmegaBaryon() const { return(omegabaryon_); }/**< Return Omega_B */
    inline double OmegaCDM() const { return(omegamat_-omegabaryon_); }/**< Return Omega_CDM */
    inline double OmegaLambda() const { return(omegaL_); }/**< Return Omega_Lambda */
    inline double OmegaDE() const { return(omegaX_); }/**< Return Omega_DE */
    inline double OmegaCurv() const { return(omegaCurv_); }/**< Return Omega_K */

    // Dark energy equation of state
    /** Return dark energy equation of state */
    inline double DEw() const { return(wX_); }

    // Matter and various energy densities
    /** Return Omega_M at redshift z=ze_*/
    double OmegaMatterZE() const ;
    /** Return Omega_Rad at redshift z=ze_*/
    double OmegaRadiationZE() const ;
    /** Return Omega_K at redshift z=ze_*/
    double OmegaCurvZE() const ; 

    /** Return Omega_M at redshift z */
    double OmegaMatter(double z) const ; 
    /** Return Omega_Rad at redshift z */
    double OmegaRadiation(double z) const ; 


    /** Critical density (in Kg/cm^3) */
    inline double CriticalDensityKgm3() const 
        { return( h_*h_*RHOCRIT_h2KGM3); }
    /** Critical density (in GeV/cm^3)*/
    inline double CriticalDensityGeVcm3() const 
        { return( h_*h_*RHOCRIT_h2GEVCM3); }
    /** Critical density (in SolarMass/Mpc^3)*/
    inline double CriticalDensityMSolMpc3() const 
        { return( h_*h_*RHOCRIT_h2MSOLMPC3); }

    /** Critical density as a function of h (in Kg/cm^3)*/
    static inline double CriticalDensityKgm3(double h)
        { return( h*h*RHOCRIT_h2KGM3); }

    /** Density in photons / cm^3  = n_gamma */
    double PhotonDensitycm3(double z=0.) const ;
    /** Baryon number density / cm^3  = n_baryons*/
    double BaryonDensitycm3(double z=0.) const ;
    /** Ratio n_baryons / n_photons */
    double EtaBaryonPhoton() const ;

    // -------------- SetEmissionRedShift() ----------------
    /** Defines the Emission redshift then computes the corresponding 
        coordinates: time, radial coordinates etc. Most of the computation, and 
        integration (ie int 1/E(z) dz) is done in this routine, which may take a 
        while. 
        @param ze The emission redshift
        @param prec The target precision used to compute the adaptative 
               integration step
        @param fginc If fginc = true, does incremental calculation of integrals
               */
    void SetEmissionRedShift(double ze,double prec=0.001,bool fginc=true);
  
    //--------- Emission redshift and coordinates
    /** Return the current emission redshift*/
    inline double GetEmissionRedShift() const { return ze_; }
    /** Return the current emission redshift*/
    inline double ZE() const { return ze_; }
    /** Return the emission point radial coordinate Xi in Mpc (Observer->0.)
    @warning DON'T USE THIS! NOT SURE WHAT IT IS!!!*/
    inline double RadialCoordinateMpc() const { return chie_*HubbleLengthMpc(); }
    /** Return the emission point radial coordinate Xi in Hubble length units 
    (Observer->0.) */
    inline double ChiCoordinate() const { return chie_; }
  
    //--------- Comoving z error
    // Added by AA Nov 2010
    /** Equivalent error in comoving distance from redshift error */
    inline double ZErr2CoDistErr(double zerr)
				{ return zerr * (1+ZE()) * ( (SPEED_OF_LIGHT_MS)/(HZE()*1e3) );}

    //--------- The distances
    //--------- Comoving distances
    /** Line of sight comoving distance, Hubble length units*/
    inline double LineOfSightComovDistance() const { return dcl_; };
     /** Line of sight comoving distance in Mpc*/
    inline double LineOfSightComovDistanceMpc() const { return dcl_*HubbleLengthMpc(); };
    /** Transverse comoving distance, Hubble length units*/
    inline double TransComovDistance() const { return dct_; };
    /** Transverse comoving distance, Hubble length units*/
    inline double TransComovDistanceMpc() const { return dct_*HubbleLengthMpc(); };

    //--------- Angular Diameter Distance D_A
    /** Return the angular diameter distance D_A at ze_ in Hubble length units
    @note the D_A (Ang.Diam.Distance) IS NOT a comoving length !*/
    inline double AngularDiameterDistance() const { return da_; }
    /** Return the angular diameter distance D_A at ze_ in Mpc
    @note the D_A (Ang.Diam.Distance) IS NOT a comoving length !*/
    inline double AngularDiameterDistanceMpc() const 
        { return da_*HubbleLengthMpc(); }
    /** Return the comoving ang. diam. distance (1+z)*D_A */
    inline double ComovAngDiamDistanceMpc() const 
        { return da_*HubbleLengthMpc()*(1+ze_); }

    //--------- Luminosity Distance D_L
    /** Luminosity distance = D_A*(1+z)*(1+z) -
        Return the luminosity distance D_A at ze_ in Hubble length units
        @note D_L IS a comoving length */
    inline double LuminosityDistance() const { return dl_; }
    /** Luminosity distance = D_A*(1+z)*(1+z) -
        Return the luminosity distance D_A at ze_ in mega parsec (Mpc)
        @note D_L IS a comoving length */
    inline double LuminosityDistanceMpc() const 
        { return dl_*HubbleLengthMpc(); }

    //--------- Volume element
    // Added by AA Nov 2010, see Hogg 1999
    /** Return the volume element: dV/dzdOmega */
    double VolEl() { return (DH()*DH()*DH()*da_*da_*(1+ze_)*(1+ze_)*Gz(ze_)); }

  
    //--------- Time 
    /** Return the look back time, from emission to today (in Hubble Time units)*/
    inline double LookBackTime()  const { return te_; }
    /** Return the look back time, in years */
    inline double LookBackTimeYear()  const { return te_*HubbleTimeYear(); }
    /** Return the time in Hubble time units between the two specified redshifts*/
    double	DeltaTime(double z1, double z2, double prec=0.001);
    /** Return the time in years between the two specified redshits*/
    inline double DeltaTimeYear(double z1, double z2, double prec=0.001)
        { return DeltaTime(z1, z2, prec)*HubbleTimeYear(); }
    /** Return the age of the universe */
    inline double AgeYearToday(double zmax=50000.)
        { return DeltaTime(0., zmax)*HubbleTimeYear(); }
    /** Return the age of the universe at the emission redshift */
    double AgeYearAtZE(double zmax=50000.)
        { return DeltaTime(ze_, zmax)*HubbleTimeYear(); }
    /** Conformal time eta, d eta = dt*(1+z)
        Return the conformal time (eta) in Hubble time units between the two 
        specified redshifts*/
    double	DeltaConformalTime(double z1, double z2, double prec=0.001);
    /** Return the conformal time (eta) in years between the two specified 
        redshifts */
    inline double DeltaConformalTimeYear(double z1, double z2, double prec=0.001)
        { return DeltaConformalTime(z1, z2, prec)*HubbleTimeYear(); }

    //--------- The universe curvature
    /** Return curvature */
    inline int    KCurvature() const { return( kcurvature_ ); }
    /** Return bool if flat universe */
    inline bool   IsFlat() const { return ( kcurvature_ == 0 ); }
    /** Return bool if open universe */
    inline bool   IsOpen() const { return ( kcurvature_ == -1 ); } 
    /** Return bool if closed universe */
    inline bool   IsClosed() const { return ( kcurvature_ == 1 ); } 

    //--------- Printing
    /** Print cosmological parameters and quantities evalulated at the current
        emission redshift */
    void          Print(ostream & os) const ;
    /** Print cosmological parameters and quantities evalulated at the current
        emission redshift */
    inline void   Print() const { Print(cout); } ;
    /** Print cosmological parameters and quantities evalulated at the current
        emission redshift */
    inline void   print(ostream & os) const { Print(os); }
    /** Print cosmological parameters and quantities evalulated at the current
        emission redshift */
    inline void   print() const { Print(cout); } ;

    //--------- Universe parameters
    /** Statistical value of the z-decoupling (data) */
    static inline double Z_DecPar()  { return Z_RECOMBINATION; }
    /** Matter-radiation equality epoch */
    double        Z_Equality() const ;


    //--------- CMB temperature and photon and neutrino densities (in Kg/m^3)
    // (today z=0) , etc ...
    /** Return CMB temperature in Kelvin */
    static inline double T_CMB() { return T_CMB_K; }
    /** Return photon density today in kg/m^3 */
    static inline double PhotonDensityKgm3() 
            { return PHOTON_DENSITY_TODAY_KGM3; }
    /** Return photon number density today per cm^3 */
    static inline double PhotonNumberDensitycm3() 
            { return PHOTON_NUMDENSITY_PER_CM3_TODAY; }
    /** Return neutrino density today in kg/m^3 */
    static inline double NeutrinoDensityKgm3() 
          { return PHOTON_AND_NU_DENSITY_TODAY_KGM3-PHOTON_DENSITY_TODAY_KGM3; }
    /** Return photon+neutrino density today in kg/m^3 */
    static inline double PhotonNuDensityKgm3() 
        { return PHOTON_AND_NU_DENSITY_TODAY_KGM3; }

    // --------- Useful Constants 
    /** Return Mpc to meter conversion factor */
    static inline double MpctoMeters() { return NMETERS_IN_MPC; }
    /** Return lightyear to meter conversion factor  */
    static inline double LightYear() { return NMETERS_IN_LIGHTYEAR; }       
    /** Return speed of light  m/sec  */
    static inline double SpeedOfLight() { return SPEED_OF_LIGHT_MS; }    
    /** Return 1 year in seconds */
    static inline double YeartoSec() { return NSECS_IN_YEAR; }       
    /** Return Newton's constant (SI units, m^3/kg/s^2)  */
    static inline double G_Newton() { return G_NEWTON_SI; }
    /** Return Planck constant (SI units, Js) */
    static inline double h_Planck() { return HPLANCK_SI; }
    /** Return Planck constant/2pi (SI units, Js) */
    static inline double hbar() { return HPLANCK_SI/2./PI; }
    /** Return Boltzmann constant (SI units, J/K)  */
    static inline double k_Boltzman() { return KBOLZMANN_SI; }
    /** Return electron charge (SI units, As or Coulomb) */
    static inline double ElectronCharge() { return ELECTRON_CHARGE_SI; }
    /** Return Kg -> GeV conversion factor  */
    static inline double KgtoGeV() 
        { return SpeedOfLight()*SpeedOfLight()/ElectronCharge()/1.e9; }
    /** Return Planck Mass (in Kg)  */
    static inline double PlanckMassKg() 
        { return sqrt(hbar()*SpeedOfLight()/G_Newton()); }
    /** Return Planck Mass (in GeV) */
    static inline double PlanckMassGeV() { return PlanckMassKg()*KgtoGeV(); }


    // ----------------------------------------------------------------------
    // Main methods to calculate cosmological quantities
    
    /** Return E(z) : Notation of Peebles / Principle of Physical Cosmology 
        (Chap 13) */
    double         Ez(double z) const;
    /** Return G(z) = 1/E(z) : integration element for the calculation of D_A*/  
    inline double  Gz(double z)  const { return ( 1./Ez(z) ) ; }
    /** Return GT(z) = 1/[(1+z)*E(z)] : integration element for the calculation 
        of age/time (dt = GTz(z)*dz) */
    inline double  GTz(double z)  const { return ( 1./((1+z)*Ez(z)) ) ; }

    /** Numerical integration of G(z) dz = 1/E(z) dz and numerical
        integration of GT(z) dz = 1/[(1+z)*E(z)] dz between [z1,z2]
    @param z1 Lower redshift bound
    @param z2 Upper redshift bound
    @param prec The target precision used to compute the adaptative 
               integration step
    @param resG Resulting integration of G(z)
    @param resGT Resulting integration of GT(z)
    */
    void NumIntegrateGzGTz(double z1, double z2, double prec, 
				double& resG, double& resGT);
    /** Numerical integration of volume element over redshifts (not solid
        angle)
    @param z1 Lower redshift bound
    @param z2 Upper redshift bound
    @param nstep Number of integration steps
    @param dVdO Resulting integration of volume element over z
    */
    void NumIntegVolElz(double z1, double z2, int nstep, 
				       double& dVdO);

    /** Value of the analytical primitive of H(z) = Integ[G(z) dz] with Lambda=0
    */
    double      IntegrateGz_NoLambda(double z);
    /** Value of the analytical primitive of Integ[GT(z) dz] with Lambda=0
    */
    double      IntegrateGTz_NoLambda(double z);
    /** Value of d G(z) / dz*/
    double      DerivateGz(double z, double dz=-1.);
    /** Value of d G(z) / dz Lambda = 0, OmegaX = 0*/
    double      DerivateGz_NoLambda(double z);
    /** Value of d G(z) / dz with Lambda */
    double      DerivateGz_Lambda(double z);
    /**  Value of d G(z) / dz with Lambda and Dark energy (X)*/
    double      DerivateGz_Lambda_X(double z);

protected:
    /** Set parameters to some default values */
    void          Initialize();
    /** Compute curvature from other component densities */
    void          ComputeCurvature();

    /** Numerical integration of G(z) dz = 1/E(z) dz and numerical
        integration of GT(z) dz = 1/[(1+z)*E(z)] dz between [z1,z2] 
    @param z1 Lower redshift bound
    @param z2 Upper redshift bound
    @param prec The target precision used to compute the adaptative 
               integration step
    @param resG Resulting integration of G(z)
    @param resGT Resulting integration of GT(z)
    */
    void          NumIntegGzGTz(double z1, double z2, double prec, 
			                        double& resG, double& resGT);

    // data members
    double h_;              ///< H0 in units of 100 km/s/Mpc : 0.5 -> 50 km/s/MPc
    double omegamat_;       ///< Total matter density in universe @z=0 (inc baryons)
    double omegarad_;       ///< Total radiation density in universe at z=0 
    double omegabaryon_;    ///< Baryon density (included in omegamat_); 
    double omegaphot_;      ///< Photon energy density, (included in omegarad_) 
    double omegaL_;         ///< cosmological constant
    double omegaX_;         ///< The dark energy density 
    double wX_;             ///< Dark energy density equation of state
    bool hasL_;             ///< True-> Non zero cosmological constant
    bool hasX_;             ///< True-> Non zero dark energy density
    double omegaCurv_;      ///< curvature density (1-( omega0_ + omegaL_) )
    int kcurvature_;        ///< 0 Flat , +1 closed, -1 Open
    double ze_;             ///< Emission redshift
    double te_;             ///< Emission time (0=today)
    //double re_;                Radial coordinate (emission)
    double chie_;           ///< Chi coordinate
    double dcl_;            ///< Line of sight comoving distance (units c/H0)
    double dct_;            ///< Transverse comoving distance (units c/H0)
    double da_;             ///< Angular diameter distance (units c/H0)
    double dl_;             ///< Luminosity distance (units c/H0)
    double EofZe_;          ///< E(ze)  
    double integGz_;        ///< Integral[G(z) dz]_[ze_ ... 0]
    double integGTz_;       ///< Integral[GT(z) dz]_[ze_ ... 0]

    int prt_dbg_level;      ///< Debug/printing level

  /*// ---- LUC version number -----
  static double LUC_version;

  // ---- Useful constants ----
  static double MpctoMeters_Cst;        // Mpc to meter conversion factor
  static double LightYear_Cst;          // LightYear to meter conversion factor
  static double SpeedOfLight_Cst;       // Speed of light  m/sec
  static double YeartoSec_Cst;          // 1 Year in second
  static double G_Newton_Cst;           // Newton's constant G (SI Units) 
  static double h_Planck_Cst;           // Planck constant     (SI Units)
  static double k_Boltzman_Cst;         // Boltzmann constant  (SI Units)
  static double k_BoltzmaneV_Cst;       // Boltzmann constant  (in eV)
  static double ElectronCharge_Cst;     // Electron charge 
  static double ElectronMassGeV_Cst;    // Electron mass in GeV
  static double ProtonMassGeV_Cst;      // Electron mass in GeV
  static double KgtoGeV_Cst;            // Kg to GeV conversion factor

  // --- Critical density 
  static double RhoCrit_h1_Cst;          // Critical density for h=1 Kg/m^3
  static double RhoCrit_h1_MSolMpc3_Cst; //   in SolarMass/Mpc^3
  static double RhoCrit_h1_GeVcm3_Cst;   //   in GeV/cm^3

  // ---- Universe parameters 
  static double T_CMB_Par;              // CMB temperature
  static double rho_gamma_Par;          // Photon density today
  static double numberdens_gamma_Par;   // Photon number density today
  static double rho_gamma_nu_Par;       // Photon+neutrino density today
  static double Z_Rec_Par;              // Recombination redshift*/
};

// ostream << operator definition (print)
inline ostream & operator << (ostream & os, SimpleUniverse const & su)
{ su.Print(os); return os; }




#endif
