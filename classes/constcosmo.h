/**
 * @file  constcosmo.h
 * @brief Holds physical, unit conversion and cosmology constants
 *
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: 17th August 2012
 * @date 17th August 2012
 *
 */
 
#ifndef CONSTCOSMO_SEEN
#define CONSTCOSMO_SEEN


// :::::::::::::::::::::::: CONSTANTS ::::::::::::::::::::::::::

// Speed of light
static const double SPEED_OF_LIGHT_CMS = 29979245800; // speed of light cm/sec
static const double SPEED_OF_LIGHT_KMS = 299792.458; // speed of light km/sec
static const double SPEED_OF_LIGHT_MS = 2.99792458e8;// speed of light m/sec

// Time
static const double NSECS_IN_YEAR = 31556925.2;      // 1 year into seconds

// Distances and angles
static const double NMETERS_IN_MPC = 3.0856e+22;     // 1 Mpc into meters
static const double NMETERS_IN_LIGHTYEAR = 9.46073042e+15;// 1 lightyear in m
static const double NRADIANS_IN_DEG = (M_PI/180.); // 1 degree in radians
static const double NSTERADIANS_IN_SQDEG = (M_PI/180.)*(M_PI/180.); // 1 square deg in steradians

// Physics
static const double G_NEWTON_SI = 6.6742e-11; // G in SI units m^3/kg/s^2
static const double HPLANCK_SI = 6.6260693e-34;  // Planck constant SI units Js
static const double KBOLZMANN_SI = 1.3806503e-23;// Boltzmann constant SI J/K
static const double KBOLZMANN_EV = 8.617342e-5;  // Boltzmann constant in eV
static const double STEFAN_BOLZ_SI = 5.670400e-8;// Stefan-Boltzman in W/m^2/K^4
static const double ZETA_3 = 1.2020569031595942854; // Zeta(3)
static const double PROTON_MASS_IN_KG = 1.67262171e-27; // Proton mass in Kg
static const double PROTON_MASS_IN_G = 1.67262171e-24; // Proton mass in g
static const double PROTON_MASS_IN_GEV = 0.938271998;  // Proton mass in GeV
static const double ELECTRON_CHARGE_SI = 1.602176462e-19;// Electric charge in SI units As/coulomb
static const double ELECTRON_CHARGE_STATC = 4.80320425e-10; // Electric charge in cgs units, statcoulombs
static const double ELECTRON_MASS_GEV = 0.510998902e-3;// Electron mass in GeV 
static const double ELECTRON_MASS_SI = 9.10938215e-31; // Electron mass SI (Kg)
static const double ELECTRON_MASS_G = 9.10938215e-28; // Electron mass in g

// Sigma = 2*PI^5*K^4/(15.*C^2*H^3) = PI^2*K^4/(60.*Hbar^3*C2)

// Atomic
static const double WAVE_LYMANLIM_METERS = 91.175e-9;// in m
static const double WAVE_LYMANLIM_ANGSTR = 911.75;  // in angstroms

// Astronomy units
static const double JANSKY_IN_WM2HZ = 1.e-26;  // 1 Jansky in Watts/m^2/Hz
static const double SOLARMASS_IN_G = 1.98844e+33; // Solar mass in grams
static const double SOLARMASS_IN_KG = 1.98844e+30; // Solar mass in kilograms
static const double GRAMPERCM3_IN_MSOLPERMPC3 = 1.477428208143872e+40; 
                                                    // 1 g/cm^3 in 1 Msol/Mpc^3
// GCm3toMsolMpc3_Cst = (MpctoMeters_Cst*100)^3 / SolarMass_Cst

// Cosmology 
static const double T_CMB_K = 2.725;  // temperature of CMB in K
static const double T_NU_K = 1.9;  // temperature of neutrinos in K
static const double FREQ_21CM_HI_IN_GHZ = 1.420405751786;
                                        // frequency of 21cm radiation in GHz

static const double RHOCRIT_h2KGM3 = 1.879e-26;// critical density in h^2 Kg/m^3
static const double RHOCRIT_h2MSOLMPC3 = 2.77536627e11;//critical density in MSol/Mpc^3
static const double RHOCRIT_h2GEVCM3 = 1.0539e-5;//critical density in in GeV/cm^3
static const double PHOTON_DENSITY_TODAY_KGM3 = 4.6417e-31;// Photon density today (kg/m^3) 
static const double PHOTON_NUMDENSITY_PER_CM3_TODAY = 410.4;// Photon number density today /cm^3 
static const double PHOTON_AND_NU_DENSITY_TODAY_KGM3 = 7.8042e-31; // Photon+neutrino density (kg/m^3) 
static const double Z_RECOMBINATION = 1088.; // Value in WMAP ApJ paper
static const double DELTA_C = 1.686;    // Overdensity for collapse (spherical model)

// PI!
static const double PI = 3.141592;



#endif
