/**
 * @file  genefluct3d.h
 * @brief Generates over-densities on a grid
 *
 * Could add more information here I think
 *
 * @author Reza Ansari, Christophe Magneville
 * Contact:
 *
 * Created on: 2008
 * @date 2008
 *
 */

#ifndef GENEFLUCT3D_SEEN
#define GENEFLUCT3D_SEEN

#include "machdefs.h"
#include <math.h>
#include "genericfunc.h"
#include "tarray.h"
#include "histerr.h"
#include "hist2err.h"
#include "perandom.h"

#include "FFTW/fftw3.h"
#include "FitsIO/fitsio.h"

#include <vector>
#include <algorithm>

#include "cosmocalcs.h"
#include "pkspectrum.h"

#define WITH_FFTW_THREAD

//DON'T UNCOMMENT, USE THE MAKEFILE #define GEN3D_FLOAT
//#define GEN3D_FLOAT 1

#if defined(GEN3D_FLOAT)
#define GEN3D_TYPE r_4
#define GEN3D_FFTW_PLAN fftwf_plan
#define GEN3D_FFTW_COMPLEX fftwf_complex
#else
#define GEN3D_TYPE r_8
#define GEN3D_FFTW_PLAN fftw_plan
#define GEN3D_FFTW_COMPLEX fftw_complex
#endif

// written by Christophe Magneville
// edited by AA

namespace SOPHYA {

//-----------------------------------------------------------------------------------
class GeneFluct3D {
public:
  GeneFluct3D(long nx,long ny,long nz,double dx,double dy,double dz,unsigned short nthread=0,int lp=0);  // Mpc
  GeneFluct3D(unsigned short nthread=0);
  virtual ~GeneFluct3D(void);

  // Comoving distance to the observer
  void SetObservator(double redshref=0.,double kredshref=0.);
    inline double DXcom(long i) {return i*Dx_ - xobs_[0];}
    inline double DYcom(long j) {return j*Dy_ - xobs_[1];}
    inline double DZcom(long k) {return k*Dz_ - xobs_[2];}
    inline double Dcom(long i,long j,long k) {
      double dx=DXcom(i), dy=DYcom(j), dz=DZcom(k); 
      return sqrt(dx*dx+dy*dy+dz*dz);
    }
  void SetCosmology(SimpleUniverse& cosmo);
  void SetGrowthFactor(GrowthFactor& growth);
  long LosComRedshift(double zinc=0.001,long npoints=-1);

  TArray< complex<GEN3D_TYPE> >& GetComplexArray(void) {return T_;}
  GEN3D_FFTW_COMPLEX * GetComplexPointer(void) {return fdata_;}
  TArray<GEN3D_TYPE>& GetRealArray(void) {return R_;}
  GEN3D_TYPE* GetRealPointer(void) {return data_;}

  // For index data_[ip]
  inline int_8 IndexR(long i,long j,long k) {return (int_8)(k+NTz_*(j+Ny_*i));}
  // For index fdata_[ip][0-1]
  inline int_8 IndexC(long i,long j,long k) {return (int_8)(k+NCz_*(j+Ny_*i));}
  // Could also index:
  // TArray< complex<r_8> >& pk = gf3d.GetComplexArray();
  //    pk(k,j,i) avec k=[0,NCz_[  j=[0,Ny_[   i=[0,Nx_[
  //    pk[IndexC(i,j,k)]
  // TArray<r_8>& rgen = gf3d.GetRealArray();
  //    rgen(k,j,i) avec k=[0,NTz_[  j=[0,Ny_[   i=[0,Nx_[
  //                mais seul k=[0,Nz_[ est utile
  //    rgen[IndexR(i,j,k)]
  // ATTENTION: TArray adresse en memoire a l'envers du C !
  //            Tarray(n1,n2,n3) == Carray[n3][n2][n1]

  vector<long> GetNpix(void) {return N_;}
  int_8 NPix(void) {return NRtot_;}
  long GetNx(void) {return Nx_;}
  long GetNy(void) {return Ny_;}
  long GetNz(void) {return Nz_;}

  // Return |K_i| module relative to pixel indices
  inline r_8 Kx(long i) {long ii=(i>Nx_/2)? Nx_-i :i; return ii*Dkx_;}
  inline r_8 Ky(long j) {long jj=(j>Ny_/2)? Ny_-j :j; return jj*Dky_;}
  inline r_8 Kz(long l) {return l*Dkz_;}

  vector<r_8> GetDinc(void) {return D_;}
  double GetDVol(void) {return dVol_;}
  double GetVol(void) {return Vol_;}

  vector<r_8> GetKinc(void) {return Dk_;}
  vector<r_8> GetKnyq(void) {return Knyq_;}
  double GetKmax(void) {return sqrt(Knyqx_*Knyqx_+Knyqy_*Knyqy_+Knyqz_*Knyqz_);}
  double GetKTmax(void) {return sqrt(Knyqx_*Knyqx_+Knyqy_*Knyqy_);}
  double GetKincMin(void)
    {vector<r_8>::const_iterator it = min_element(Dk_.begin(), Dk_.end()); return *it;}
  double GetKincMax(void)
    {vector<r_8>::const_iterator it = max_element(Dk_.begin(), Dk_.end()); return *it;}
  double GetKTincMin(void) {return min(Dk_[0],Dk_[1]);}
  double GetKTincMax(void) {return max(Dk_[0],Dk_[1]);}

  void ComputeFourier0(GenericFunc& pk_at_z);
  void ComputeFourier(GenericFunc& pk_at_z);
  void FilterByPixel(void);

  void ComputeReal(void);
  void ApplyGrowthFactor(int type_evol=1);

  void ReComputeFourier(void);

  int  ComputeSpectrum(HistoErr& herr);
  int  ComputeSpectrum2D(Histo2DErr& herr);
  int  ComputeSpectrum(HistoErr& herr,double sigma,bool pixcor);
  int  ComputeSpectrum2D(Histo2DErr& herr,double sigma,bool pixcor);

  int_8 VarianceFrReal(double R,double& var);
  int_8 MeanSigma2(double& rm,double& rs2,double vmin=1.,double vmax=-1.
                  ,bool useout=false,double vout=0.);
  int_8 MinMax(double& xmin,double& xmax,double vmin=1.,double vmax=-1.);
  int_8 NumberOfBad(double vmin=-1.e+150,double vmax=1.e+150);
  int_8 SetToVal(double vmin, double vmax,double val0=0.);
  void  ScaleOffset(double scalecube=1.,double offsetcube=0.);

  void TurnFluct2Mass(void);
  double TurnFluct2MeanNumber(double val_by_mpc3);
  double ApplyPoisson(void);
  double TurnNGal2Mass(FunRan& massdist,bool axeslog=false);
  double TurnNGal2MassQuick(SchechterMassDist& schmdist);
  double TurnMass2Flux(void);
  //void AddAGN(double lfjy,double lsigma,double powlaw=0.);
  void AddNoise2Real(double snoise,int type_evol=0);

  void WriteFits(string cfname,int bitpix=FLOAT_IMG);
  void ReadFits(string cfname);

  void WritePPF(string cfname,bool write_real=true);
  void ReadPPF(string cfname);
  void WriteSlicePPF(string cfname);
  void NTupleCheck(POutPersist &pos,string ntname,unsigned long nent);

  void SetPrtLevel(int lp=0) {lp_ = lp;}
  void Print(void);

//-------------------------------------------------------------------

protected:
  void init_default(void);
  void setsize(long nx,long ny,long nz,double dx,double dy,double dz);
  void setalloc(void);
  void setpointers(bool from_real);
  void init_fftw(void);
  void delete_fftw(void);
  long manage_coefficients(void);
  double compute_power_carte(void);
  void check_array_alloc(void);
  inline double pixelfilter(double x)
    {return (x<0.025) ? 1.-x*x/6.*(1.-x*x/20.): sin(x)/x;}

  // real space values
  long Nx_,Ny_,Nz_;  vector<long> N_;
  long NCz_,NTz_;
  int_8 NRtot_;

  double Dx_,Dy_,Dz_;  vector<double> D_;

  // k space values
  double Dkx_,Dky_,Dkz_;  vector<double> Dk_;
  double Knyqx_,Knyqy_,Knyqz_;  vector<double> Knyq_;
  double Dk3_;
  double dVol_, Vol_;

  // management of FFT
  bool is_set_fft_plan;
  GEN3D_FFTW_PLAN pf_,pb_;
  unsigned short nthread_;
  int lp_;

  // storage of data cube and pointers
  bool array_allocated_;  // true if array has been allocated
  unsigned short array_type; // 0=empty, 1=real, 2=complex
  TArray< complex<GEN3D_TYPE> > T_;
  GEN3D_FFTW_COMPLEX *fdata_;
  TArray<GEN3D_TYPE> R_;
  GEN3D_TYPE *data_;

  // the observer
  SimpleUniverse *cosmo_;
  GrowthFactor *growth_;
  double redsh_ref_,kredsh_ref_,dred_ref_;
  double loscom_ref_,dtrc_ref_, dlum_ref_, dang_ref_;
  double nu_ref_, dnu_ref_ ;
  double xobs_[3];
  double loscom_min_, loscom_max_;
  vector<double> zred_, loscom_;
  double loscom2zred_min_, loscom2zred_max_;
  vector<double> loscom2zred_;

};

} // end of namespace SOPHYA

#endif
