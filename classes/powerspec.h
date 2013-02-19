/**
 * @file  powerspec.h
 * @brief Compute power spectrum from an array of over-densities
 *
 * Could add more information here I think
 *
 * @author Alex Abate and Reza Ansari
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2008
 * @date 2008
 *
 */

#ifndef POWERSPEC_H_SEEN
#define POWERSPEC_H_SEEN

#include "machdefs.h"
#include "sopnamsp.h"
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>

#include <typeinfo>
#include "fabtwriter.h"

#include "array.h"
#include "hisprof.h"
#include "histerr.h"

#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "stsrand.h"

#include "sinterp.h"

#define PI 3.141592

// Reads in galaxy fluctuation grid and computes power spectrum then writes it to a file
class PowerSpec
{
public:
	PowerSpec(TArray<r_8>,double,double,double); // Reads in galaxy fluctuation grid 
	PowerSpec(TArray<r_8>,double); // Reads in galaxy fluctuation grid 
	void ComputeFourier(TArray<r_8>&, TArray< complex< r_8 > >&);  // Fourier transforms galaxy fluctuation grid 
	void ComputeFourierBack(TArray< complex<r_8> >& four_array, TArray<r_8>& real_array);
	void SetDKDR(); // Calc specs of FT grid
	double AccumulatePowerSpectra(HProf& hp,bool pixcor=true,double maxk=1000,double Err=0,bool undamp=true,double tol_corr=0.15, double snoise=0); 
	void Conv1DGaussFT(double Gsig, double Gmean=0);
	void FourierSpaceWindow(HProf& hpW);
	void GetObsPS(HProf&, r_4, HProf&, HProf&, r_4);
	// writes power spectra to a file
	void WritePS(string fname,HProf& Pdata,r_4 Voldata,HProf& PSimlss,HProf& PSimlssf,r_4 Volsimlss,double meandens=0); 
	void WritePS(string fname,HProf& Pdata,r_4 Voldata,string PSimlssfile,double meandens=0);// supercedes the above
	void Write1PS(string fname,HProf&,double Vol);
	void Write2PS(string fname,HProf&,HProf&,double Vol,double md=0,double mdf=0);
			
	/* Minor functions */
	//void SetMaxKrad(double maxk){maxk_=maxk;};
	double ReturnKxMax(){kxmax=four_.SizeX()*Dkx_; return kxmax;};
	double ReturnKyMax(){kymax=four_.SizeY()*Dky_/2; return kymax;};
	double ReturnKzMax(){kzmax=four_.SizeZ()*Dkz_/2; return kzmax;};
	inline TArray<r_8>& ReturnRealDist() { return drho_; };
	inline TArray< complex< r_8 > >& ReturnFour() { return four_; };
	void ReturnConvDist(TArray<r_8>& drhoconv) { drhoconv = drhoconv_; };
	TVector<r_8> ReturnPk(){return Pobs_;};
	TVector<r_8> Returnk(){return kobs_;};
	void ZeroSizeArrays(){drho_.ZeroSize(); four_.ZeroSize();};
	void ZeroFourArray(){ four_.ZeroSize(); };
	void Setzc(double zc){ zc_=zc; };

	/* Added functions AA June 2011 */
	void CheckPZErr();
	void BlowUpCheck();

	/* Class Variables */
	TArray<r_8> drho_;	// Galaxy fluctuation array
	TArray<r_8> drhoconv_;	// Convolved galaxy fluctuation array
	TArray< complex< r_8 > > four_; // FT of galaxy fluctuation array
	TArray< complex< r_8 > > fourconv_; // FT of galaxy fluctuation array convolved with Gaussian
	TVector<r_8> Pobs_, kobs_;
	r_4 Vol_;			// Volume of survey
	sa_size_t Nx_,Ny_,Nz_;   // Number of pixels of FT grid
	int_8 NRtot_;		// total number of pixels
	long NCx_;			// Nz_/2+1
	double Dx_,Dy_,Dz_;	// Pixel size of FT grid
	double Dkx_,Dky_,Dkz_, kxmax, kymax, kzmax; // These are calculated from the above
	double maxk_,Err_;     // maximum k radial used in power spectrum computation, photo-z error
	double tol_corr_;      // tolerance on undamping correction
	bool undamp_;  // radial coordinate of real array has error/no error
	double zc_;
	
// Filter by the pixel window function (parallelepiped/Top Hat cube)
// FT = 1/(dx*dy*dz)*Int[{-dx/2,dx/2},{-dy/2,dy/2},{-dz/2,dz/2}]
//                   e^(ik_x*x) e^(ik_y*y) e^(ik_z*z) dxdydz
//    = 2/(k_x*dx) * sin(k_x*dx/2)  * (same y) * (same z)
// Avoid divergence at kx = 0 with the approximation: 
//								sin(y)/y = 1 - y^2/6*(1-y^2/20)
//								with y = k_x*dx/2
// Approximation works well until y~1
	protected:
		inline double pixelfilter(double x)// x = kR (0.5*wavevector*pixel size)
			{return (x<0.025) ? 1.-x*x/6.*(1.-x*x/20.): sin(x)/x;}
};

#endif


