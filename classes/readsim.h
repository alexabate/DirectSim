/**
 * @file  readsim.h
 * @brief 
 *
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2010
 * @date 2010
 *
 */
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
#include "swfitsdtable.h"

// DirectSim
#include "cosmocalcs.h"

/** @class ReadSim
  *
  * Read in the SimLSS mass distribution and the header of its FITS file
  * Optionally can also fudge the delta<-1 cells
  */
class ReadSim
{
public:

    /** Constructor 
        @param drho          mass distribution
        @param nbadplanes    number of junk planes of cube                    */
	ReadSim(TArray<r_8> drho, int nbadplanes=0);
	
	/** Read header of FITS file                                              */
	void ReadHeader(FitsInOutFile& fin);
	
	/** Set grid cells with negative mass to zero                             */
	sa_size_t CleanNegativeMassCells();
	
	/** Return Fourier space step in x-dimension                              */
	double ReturnDKX() { return Dkx_; };
	
	/** Return Fourier space step in y-dimension                              */
	double ReturnDKY() { return Dky_; };
	
	/** Return Fourier space step in z-dimension                              */
	double ReturnDKZ() { return Dkz_;};
	
	/** Return pixel size in x-dimension                                      */
	double ReturnDX() { return Dx_; };
	
	/** Return pixel size in y-dimension                                      */
	double ReturnDY() { return Dy_; };
	
	/** Return pixel size in z-dimension                                      */
	double ReturnDZ() { return Dz_; };
	
	/** Return number of pixels in x-dimension                                */
	long ReturnNX() { return Nx_; };
	
	/** Return number of pixels in y-dimension                                */
	long ReturnNY() { return Ny_; };
	
	/** Return number of pixels in z-dimension                                */
	long ReturnNZ() { return Nz_; };
	
	/** Return pixel volume                                                   */
	double ReturnPixVol() { return Dx_*Dy_*Dz_; };
	
	/** Return grid volume                                                    */
	double ReturnCubeVol() { return Nx_*Ny_*Nz_*Dx_*Dy_*Dz_; };
	
	/** Return mass distribution                                              */
	TArray<r_8> MassArray() { return drho_; };

protected:
	sa_size_t Nx_;    /**< number of pixels in x-dimension                    */
	sa_size_t Ny_;    /**< number of pixels in y-dimension                    */
	sa_size_t Nz_;    /**< number of pixels in z-dimension                    */
	double Dx_;       /**< pixel size in x-dimension                          */
	double Dy_;       /**< pixel size in y-dimension                          */
	double Dz_;       /**< pixel size in z-dimension                          */
	double Dkx_;      /**< pixel size in Fourier space x-dimension            */
	double Dky_;      /**< pixel size in Fourier space y-dimension            */
	double Dkz_;      /**< pixel size in Fourier space z-dimension            */
	double zref_;     /**< redshift of center pixel                           */
	long idmidz_;     /**< index of center x-dim pixel (1st pixel is index 1) */
	double idmidy_;   /**< index of center y-dim pixel (1st pixel is index 1) */
	double idmidx_;   /**< index of center z-dim pixel (1st pixel is index 1) */
	TArray<r_8> drho_;/**< arrays defined as Nz,Ny,Nx */
	TArray<r_8> d_;   /**< */

};
