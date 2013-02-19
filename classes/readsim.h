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
#include "swfitsdtable.h"
#include "luc.h"


class ReadSim
//Function of this class is really just to read in the SimLSS mass distribution
//and the header of its FITS file
//Optionally can also fudge the delta<-1 cells
{
public:
	ReadSim(TArray<r_8>, int nbadplanes=0);
	void ReadHeader(FitsInOutFile& fin);
	sa_size_t CleanNegativeMassCells();
	
	double ReturnDKX(){return Dkx_;}
	double ReturnDKY(){return Dky_;}
	double ReturnDKZ(){return Dkz_;}
	double ReturnDX(){return Dx_;}
	double ReturnDY(){return Dy_;}
	double ReturnDZ(){return Dz_;}
	long ReturnNX(){return Nx_;}
	long ReturnNY(){return Ny_;}
	long ReturnNZ(){return Nz_;}
	double ReturnPixVol(){return Dx_*Dy_*Dz_;}
	double ReturnCubeVol(){return Nx_*Ny_*Nz_*Dx_*Dy_*Dz_;}
	TArray<r_8> MassArray() {return drho_; }

protected:
	sa_size_t Nx_,Ny_,Nz_; // number of pixels in each direction: NOTE arrays defined as (Nz, Ny, Nx)
	double Dx_,Dy_,Dz_; // pixel size
	double Dkx_,Dky_,Dkz_; // Fourier space pixel
	double zref_; // redshift of centre pixel
	long idmidz_, idmidy_, idmidx_; // indexs of centre pixel (assuming 1st pixel is index 1)
	TArray<r_8> drho_, d_;

};
