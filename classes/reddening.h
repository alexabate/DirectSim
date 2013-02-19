/**
 * @file  reddening.h
 * @brief Calculate Cardelli and Calzetti reddening laws
 *
 * Could add more information here I think
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2010
 * @date 2010
 *
 */
 
#ifndef REDDENING_H_SEEN 
#define REDDENING_H_SEEN 

#include "machdefs.h"
#include "sopnamsp.h"
#include "pexceptions.h"
#include <math.h>
#include <iostream>

//#include "genericfunc.h"

/** Reddening class
  *
  * Class that calculates the Cardelli and Calzetti reddening laws
  *
  */
class Reddening //: public GenericFunc
{
public :
    /** Constructor */
	Reddening(); 
	//virtual Reddening() { }

    /** Returns the Calzetti reddening law \f$ k(\lambda) \f$ for starburst 
        galaxies. Reddening can then be applied to a galaxy spectrum, \f$ S(\lambda) \f$,
        via: \f$ S^{red}(\lambda) = S(\lambda)10^{-0.4k(\lambda)E(B-v)} \f$
        @param lambda   wavelength in meters (in rest-frame of galaxy)        */
	double Calzetti(double lambda);
	
	/** Returns the Cardelli reddening law \f$ k(\lambda) \f$ for elliptical and 
        spiral galaxies. Reddening can then be applied to a galaxy spectrum, \f$ S(\lambda) \f$,
        via: \f$ S^{red}(\lambda) = S(\lambda)10^{-0.4k(\lambda)E(B-v)} \f$
        @param lambda   wavelength in meters (in rest-frame of galaxy)        */
	double Cardelli(double,double);
		
};

//The Calzetti law is given by :
//k(λ) = 2.659(−1.857 + 1.040/λ) + Rv,	if 0.63 μm ≤ λ ≤ 2.20 μm
//k(λ) = 2.659(−2.156 + 1.509/λ − 0.198/λ^2+ 0.011/λ^3 ) + Rv if 0.12 μm ≤ λ ≤ 0.63 μm
//with Rv =4.05 and E(B−V) varies from 0 to 0.3.

//The Cardelli law:
//k(λ) = a(x) + b(x)/Rv
//where x = 1/λ, and a and b take different forms depending on λ :
//a(x) = 0.574x^1.61 
//b(x) = −0.527x^1.61
//if 0.3μm−1 ≤ x ≤ 1.1μm−1
//a(x) = 1 + 0.17699y − 0.50447y^2 − 0.02427y^3 + 0.72085y^4 + 0.01979y^5 − 0.77530y^6 + 0.32999y^7
//b(x) = 1.41338y + 2.28305y^2 + 1.07233y^3 − 5.38434y^4 − 0.62251y^5 + 5.30260y^6 − 2.09002y^7 
//if 1.1μm−1 ≤ x ≤ 3.3μm−1
//a(x) = 1.752 − 0.31x − 0.104/[(x − 4.67)^2 + 0.341] + Fa(x) 
//b(x) = −3.090 + 1.825x + 1.206/[(x − 4.62)^2 + 0.263] + Fb(x)
//Fa(x) = −0.04473(x − 5.9)^2 − 0.009779(x − 5.9)^3	8μm−1 ≥ x ≥ 5.9μm−1
//0 x ≤ 5.9μm−1
//Fb(x) = 􏰞0.2130(x − 5.9)2 + 0.1207(x − 5.9)3 8μm−1 ≥ x ≥ 5.9μm−1
// 0 x ≤ 5.9μm−1
//The value of Rv is usually taken to be equal to 3.1 and 
//E(B − V ) varies from 0 to 0.1 for early type galaxies and from 0 to 0.3 for late types.
//--- SED class 
// to load in SEDs: 
// two columns: wavelength (m) ;  flux (W/m^2/Hz)




#endif

