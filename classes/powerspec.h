/**
 * @file  powerspec.h
 * @brief Compute power spectrum from an array of over-densities
 *
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


#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <typeinfo>

// sophya
#include "machdefs.h"
#include "sopnamsp.h"
#include "fabtwriter.h"
#include "array.h"
#include "hisprof.h"
#include "histerr.h"
#include "histinit.h"
#include "fiosinit.h"
#include "fitsioserver.h"
#include "fftwserver.h"
#include "stsrand.h"

// DirectSim
#include "sinterp.h"
#include "constcosmo.h"

/** @class PowerSpec
  *
  * Reads in galaxy fluctuation grid and computes power spectrum 
  */
class PowerSpec
{
public:
    /** Constructor 
        @param delta_gal    galaxy fluctuation grid
        @param dx           grid pixel size in x-dimension (Mpc)
        @param dy           grid pixel size in y-dimension (Mpc)
        @param dz           grid pixel size in z-dimension (Mpc)              */
	PowerSpec(TArray<r_8> delta_gal, double dx, double dy, double dz); 
	
	/** Constructor 
        @param delta_gal    galaxy fluctuation grid
        @param dx           grid pixel size in x-dimension (Mpc)              */
	PowerSpec(TArray<r_8> delta_gal, double dr); 
	
	/** Compute Fourier transform of galaxy fluctuation grid
	    @param real_grid      real value array 
	    @param fourier_grid   Fourier space array                             */
	void ComputeFourier(TArray<r_8>& real_grid, TArray< complex< r_8 > >& fourier_grid);  
	
	/** Compute Fourier transform of Fourier space grid back to real space 
	    @param fourier_grid   Fourier space array
	    @param real_grid      real value array                                */
	void ComputeFourierBack(TArray< complex<r_8> >& fourier_grid, TArray<r_8>& real_grid);
	
	/** Compute spacing of Fourier transformed grid ie \f$ dk_x =  2\pi/(N_xdx) \f$*/
	void SetDKDR(); 
	
	/** Compute power spectrum from Fourier transformed galaxy fluctuation grid.
	    and fill histogram. Power spectrum is a one-dimensional function of
	    wavevector length (universe is isotropic)
	    @param hp        histogram of power spectrum values (to be filled) 
	    @param pixcor    correct for the grid pixel shape 
	    @param maxk      max wavevector length
	    @param Err       photo-z error
	    @param undamp    if true undamp power spectrum given photo-z error size
	    @param tol_corr  maximum value can to divide power spectrum by 
	                     (else numerically unstable)
	    @param snoise    shot noise: not implemented currently                */
	double AccumulatePowerSpectra(HProf& hp, bool pixcor=true, double maxk=1000,
	     double Err=0, bool undamp=true, double tol_corr=0.15, double snoise=0); 
	
	/** Convolve Fourier transformed grid with a 1D Gaussian in z-dimension 
	    direction. @warning may not work or be useful
	    @param Gsig    standard deviation of Gaussian
	    @param Gmean   mean of Gaussian                                       */
	void Conv1DGaussFT(double Gsig, double Gmean=0);
	
	/** Compute Fourier space window function? Not sure what this does 
	    @warning may not work or be useful                                    */
	void FourierSpaceWindow(HProf& hpW);
	
	/** Remove distortion effect of setting \f$ \delta<-1 \f$ grid cells to 
	    equal -1 in estimated power spectrum. @warning may not work or be useful
	    @param hp      raw estimated power spectrum
	    @param volCat  volume of catalog power spectrum estimated from
	    @param hpsim   power spectrum estimated from simulation wo/distortion
	    @param hpsimf  power spectrum estimated from simulation w/distortion
	    @param volSim  volume of simulation                                   */
	void GetObsPS(HProf& hp, r_4 volCat, HProf& hpsim, HProf& hpsimf, r_4 volSim);
	
	/** Write power spectra to a file. None of the power spectra read in as
	    arguments are normalised yet. Columns written to the file are (1) k-values, 
	    (2) estimated P(k), (3) raw estimated P(k), (4) raw P(k), (5) simulation P(k),
	    (6) simulation with distortion P(k)
	    @param fname       file to write power spectra to
	    @param Pdata       raw estimated power spectrum
	    @param volData     volume of catalog used to estimate power spectrum
	    @param Psim        power spectrum estimated from simulation wo/distortion
	    @param Psimf       power spectrum estimated from simulation w/distortion
	    @param volSim      volume of simulation
	    @param meandens    mean density of simulation with distortion         */
	void WritePS(string fname, HProf& Pdata, r_4 Voldata, HProf& PSimlss, 
	                         HProf& PSimlssf, r_4 Volsimlss, double meandens=0); 
	
	/** Write power spectra to a file. Columns written to the file are (1) k-values, 
	    (2) estimated P(k), (3) raw estimated P(k), (4) raw P(k), (5) simulation P(k),
	    (6) simulation with distortion P(k) @note supercedes the other WritePS
	    @param fname       file to write power spectra to
	    @param Pdata       raw estimated power spectrum
	    @param volData     volume of catalog used to estimate power spectrum 
	    @param Psimfile    file to read power spectra estimated from simulation 
	                       w/ and wo/distortion from
	    @param meandens    mean density of simulation with distortion         */
	void WritePS(string fname, HProf& Pdata, r_4 Voldata, 
	                                     string PSimlssfile, double meandens=0);
	
	/** Write power spectrum to a file 
	    @param fname    file to write power spectrum to
	    @param hp       unnormalized power spectrum
	    @param Vol      volume power spectrum estimated from                  */
	void Write1PS(string fname, HProf& hp, double Vol);
	
	/** Write two power spectra to a file 
	    @param fname    file to write power spectrum to
	    @param hp1      unnormalized power spectrum 1
	    @param hp2      unnormalized power spectrum 2
	    @param Vol      volume power spectrum estimated from           
	    @param md1      mean density of distribution power spectrum 1 estimated from
	    @param md2      mean density of distribution power spectrum 2 estimated from    */
	void Write2PS(string fname, HProf& hp, HProf& hp2, double Vol, double md1=0, double md2=0);
			
	// Minor functions
	
	//void SetMaxKrad(double maxk){maxk_=maxk;};
	
	/** Return maximum wavevector in the x-dimension                          */
	double ReturnKxMax(){kxmax=four_.SizeX()*Dkx_; return kxmax;};
	
	/** Return maximum wavevector in the y-dimension                          */
	double ReturnKyMax(){kymax=four_.SizeY()*Dky_/2; return kymax;};
	
	/** Return maximum wavevector in the z-dimension                          */
	double ReturnKzMax(){kzmax=four_.SizeZ()*Dkz_/2; return kzmax;};
	
	/** Return galaxy fluctuation grid                                        */
	inline TArray<r_8>& ReturnRealDist() { return drho_; };
	
	/** Return Fourier transformed galaxy fluctuation grid                    */
	inline TArray< complex< r_8 > >& ReturnFour() { return four_; };
	
	/** Return galaxy fluctuation grid convolved with 1D Gaussian             */
	void ReturnConvDist(TArray<r_8>& drhoconv) { drhoconv = drhoconv_; };
	
	/** Return estimated power spectrum                                       */
	TVector<r_8> ReturnPk() { return Pobs_; };
	
	/** Return k values of estimated power spectrum                           */
	TVector<r_8> Returnk() { return kobs_; };
	
	/** Zero the size of the galaxy fluctuation grid and FT of galaxy 
	    fluctuation grid                                                      */
	void ZeroSizeArrays() { drho_.ZeroSize(); four_.ZeroSize(); };
	
	/** Zero the size of the FT of galaxy fluctuation grid                    */
	void ZeroFourArray() { four_.ZeroSize(); };
	
	/** Set redshift of power spectrum */
	void Setzc(double zc) { zc_=zc; };


	/* Added functions AA June 2011 */
	void CheckPZErr();

	void BlowUpCheck();

protected:
	TArray<r_8> drho_;	            /**< galaxy fluctuation grid              */
	TArray<r_8> drhoconv_;	        /**< galaxy fluctuation grid convolved with 1D Gaussian  */
	TArray< complex< r_8 > > four_; /**< FT of galaxy fluctuation grid                       */
	TArray< complex< r_8 > > fourconv_; /**< FT of galaxy fluc. grid convolved with Gaussian */
	TVector<r_8> Pobs_;             /**< estimated power spectrum             */
	TVector<r_8> kobs_;             /**< k values of estimated power spectrum */
	r_4 Vol_;			            /**< volume of survey                     */
	sa_size_t Nx_;                  /**< number of grid pixels in x-dimension */
	sa_size_t Ny_;                  /**< number of grid pixels in y-dimension */
	sa_size_t Nz_;                  /**< number of grid pixels in z-dimension */
	int_8 NRtot_;		            /**< total number of pixels               */
	long NCx_;			            /**< Nz_/2+1                              */
	double Dx_;                     /**< size of grid pixel in x-dimension    */
	double Dy_;                     /**< size of grid pixel in y-dimension    */
	double Dz_;	                    /**< size of grid pixel in z-dimension    */
	double Dkx_;         /**< size of grid pixel in Fourier space x-dimension */
	double Dky_;         /**< size of grid pixel in Fourier space y-dimension */
	double Dkz_;         /**< size of grid pixel in Fourier space z-dimension */
	double kxmax;        /**< (2PI)/(Nx_*Dx_)                                 */
	double kymax;        /**< (2PI)/(Ny_*Dy_)                                 */
	double kzmax;        /**< (2PI)/(Nz_*Dz_)                                 */
	double maxk_;    /**< maximum k radial used in power spectrum computation */
	double Err_;         /**< photo-z error                                   */
	double tol_corr_;    /** tolerance on undamping correction                */
	bool undamp_;        /**< if true do undamping correction                 */
	double zc_;          /**< redshift of power spectrum                      */
	
// Filter by the pixel window function (parallelepiped/Top Hat cube)
// FT = 1/(dx*dy*dz)*Int[{-dx/2,dx/2},{-dy/2,dy/2},{-dz/2,dz/2}]
//                   e^(ik_x*x) e^(ik_y*y) e^(ik_z*z) dxdydz
//    = 2/(k_x*dx) * sin(k_x*dx/2)  * (same y) * (same z)
// Avoid divergence at kx = 0 with the approximation: 
//								sin(y)/y = 1 - y^2/6*(1-y^2/20)
//								with y = k_x*dx/2
// Approximation works well until y~1
	
    inline double pixelfilter(double x) // x = kR (0.5*wavevector*pixel size)
			{ return (x<0.025) ? 1.-x*x/6.*(1.-x*x/20.): sin(x)/x; };
};

#endif
