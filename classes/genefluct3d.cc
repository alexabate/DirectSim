#include "sopnamsp.h"
#include "machdefs.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "tarray.h"
#include "pexceptions.h"
#include "perandom.h"
#include "srandgen.h"

#include "fabtcolread.h"
#include "fabtwriter.h"
#include "fioarr.h"
#include "ntuple.h"

#include "arrctcast.h"

#include "constcosmo.h"
#include "geneutils.h"
#include "schechter.h"

#include "genefluct3d.h"

// written by Christophe Magneville
// edited by AA

//#define GEN3D_FLOAT

#if defined(GEN3D_FLOAT)
#define GEN3D_FFTW_INIT_THREADS       fftwf_init_threads
#define GEN3D_FFTW_CLEANUP_THREADS    fftwf_cleanup_threads
#define GEN3D_FFTW_PLAN_WITH_NTHREADS fftwf_plan_with_nthreads
#define GEN3D_FFTW_PLAN_DFT_R2C_3D    fftwf_plan_dft_r2c_3d
#define GEN3D_FFTW_PLAN_DFT_C2R_3D    fftwf_plan_dft_c2r_3d
#define GEN3D_FFTW_DESTROY_PLAN       fftwf_destroy_plan
#define GEN3D_FFTW_EXECUTE            fftwf_execute
#else
#define GEN3D_FFTW_INIT_THREADS       fftw_init_threads
#define GEN3D_FFTW_CLEANUP_THREADS    fftw_cleanup_threads
#define GEN3D_FFTW_PLAN_WITH_NTHREADS fftw_plan_with_nthreads
#define GEN3D_FFTW_PLAN_DFT_R2C_3D    fftw_plan_dft_r2c_3d
#define GEN3D_FFTW_PLAN_DFT_C2R_3D    fftw_plan_dft_c2r_3d
#define GEN3D_FFTW_DESTROY_PLAN       fftw_destroy_plan
#define GEN3D_FFTW_EXECUTE            fftw_execute
#endif

#define MODULE2(_x_) ((double)((_x_).real()*(_x_).real() + (_x_).imag()*(_x_).imag()))

namespace SOPHYA {

// Constructor
GeneFluct3D::GeneFluct3D(long nx, long ny, long nz, double dx, double dy, double dz,
                                                 unsigned short nthread, int lp)
{
    init_default();

    lp_ = lp;
    nthread_ = nthread;

    setsize(nx,ny,nz,dx,dy,dz);
    setalloc();
    setpointers(false);
    init_fftw();
};


// default constructor
GeneFluct3D::GeneFluct3D(unsigned short nthread)
{
    init_default();
    setsize(2,2,2,1.,1.,1.);
    nthread_ = nthread;
    setalloc();
    setpointers(false);
    init_fftw();
};


void GeneFluct3D::init_default(void)
{
 Nx_ = Ny_ = Nz_ = 0;
 is_set_fft_plan = false;
 nthread_ = 0;
 lp_ = 0;
 array_allocated_ = false; array_type = 0;
 cosmo_ = NULL;
 growth_ = NULL;
 redsh_ref_ = -999.;
 kredsh_ref_ = 0.;
 dred_ref_ = -999.;
 loscom_ref_ = -999.;
 dtrc_ref_ = dlum_ref_ = dang_ref_ = -999.;
 nu_ref_ = dnu_ref_ = -999.;
 loscom_min_ = loscom_max_ = -999.;
 loscom2zred_min_ = loscom2zred_max_ = 0.;
 xobs_[0] = xobs_[1] = xobs_[2] = 0.;
 zred_.resize(0);
 loscom_.resize(0);
 loscom2zred_.resize(0);
}

void GeneFluct3D::setsize(long nx,long ny,long nz,double dx,double dy,double dz)
{
 if(lp_>1) cout<<"--- GeneFluct3D::setsize: N="<<nx<<","<<ny<<","<<nz
                <<" D="<<dx<<","<<dy<<","<<dz<<endl;
 if(nx<=0 || dx<=0.) {
   const char *bla = "GeneFluct3D::setsize_Error: bad value(s) for nn/dx";
   cout<<bla<<endl; throw ParmError(bla);
 }

 // Sizes of the tables
 Nx_ = nx;
 Ny_ = ny;  if(Ny_ <= 0) Ny_ = Nx_;
 Nz_ = nz;  if(Nz_ <= 0) Nz_ = Nx_;
 N_.resize(0); N_.push_back(Nx_); N_.push_back(Ny_); N_.push_back(Nz_);
 NRtot_ = Nx_*Ny_*Nz_; // number of pixels in the survey
 NCz_ =  Nz_/2 +1;
 NTz_ = 2*NCz_;

 // step in space (Mpc)
 Dx_ = dx;
 Dy_ = dy; if(Dy_ <= 0.) Dy_ = Dx_;
 Dz_ = dz; if(Dz_ <= 0.) Dz_ = Dx_;
 D_.resize(0); D_.push_back(Dx_); D_.push_back(Dy_); D_.push_back(Dz_);
 dVol_ = Dx_*Dy_*Dz_;
 Vol_ = (Nx_*Dx_)*(Ny_*Dy_)*(Nz_*Dz_);

 // step in Fourier space (Mpc^-1)
 Dkx_ = 2.*M_PI/(Nx_*Dx_);
 Dky_ = 2.*M_PI/(Ny_*Dy_);
 Dkz_ = 2.*M_PI/(Nz_*Dz_);
 Dk_.resize(0); Dk_.push_back(Dkx_); Dk_.push_back(Dky_); Dk_.push_back(Dkz_);
 Dk3_ = Dkx_*Dky_*Dkz_;
 
 // Nyquist frequency in k (Mpc^-1)
 Knyqx_ = M_PI/Dx_;
 Knyqy_ = M_PI/Dy_;
 Knyqz_ = M_PI/Dz_;
 Knyq_.resize(0); Knyq_.push_back(Knyqx_); Knyq_.push_back(Knyqy_); Knyq_.push_back(Knyqz_);
}

void GeneFluct3D::setalloc(void)
{
#if defined(GEN3D_FLOAT)
 if(lp_>1) cout<<"--- GeneFluct3D::setalloc FLOAT ---"<<endl;
#else
 if(lp_>1) cout<<"--- GeneFluct3D::setalloc DOUBLE ---"<<endl;
#endif
 // Sizing of complex table<r_8>
 // ATTENTION: TArray is arranged in memory the opposite to C
 //            Tarray(n1,n2,n3) == Carray[n3][n2][n1]
 sa_size_t SzK_[3] = {NCz_,Ny_,Nx_}; // a l'envers
 try {
   T_.ReSize(3,SzK_);
   array_allocated_ = true; array_type=0;
   if(lp_>1) cout<<"  allocating: "<<T_.Size()*sizeof(complex<GEN3D_TYPE>)/1.e6<<" Mo"<<endl;
 } catch (...) {
   cout<<"GeneFluct3D::setalloc_Error: Problem allocating T_"<<endl;
 }
 T_.SetMemoryMapping(BaseArray::CMemoryMapping);
}

void GeneFluct3D::setpointers(bool from_real)
{
 if(lp_>1) cout<<"--- GeneFluct3D::setpointers ---"<<endl;
 if(from_real) T_ = ArrCastR2C(R_);
   else        R_ = ArrCastC2R(T_);
 // Fill the pointers
 fdata_ = (GEN3D_FFTW_COMPLEX *) (&T_(0,0,0));
 data_ = (GEN3D_TYPE *) (&R_(0,0,0));
}

void GeneFluct3D::init_fftw(void)
{
 if( is_set_fft_plan ) delete_fftw();

 // --- Initialisation of fftw3 (attention data is over-written on initalisation)
 if(lp_>1) cout<<"--- GeneFluct3D::init_fftw ---"<<endl;
#ifdef WITH_FFTW_THREAD
 if(nthread_>0) {
   cout<<"...Computing with "<<nthread_<<" threads"<<endl;
   GEN3D_FFTW_INIT_THREADS();
   GEN3D_FFTW_PLAN_WITH_NTHREADS(nthread_);
 }
#endif
 if(lp_>1) cout<<"...forward plan"<<endl;
 pf_ = GEN3D_FFTW_PLAN_DFT_R2C_3D(Nx_,Ny_,Nz_,data_,fdata_,FFTW_ESTIMATE);
 if(lp_>1) cout<<"...backward plan"<<endl;
 pb_ = GEN3D_FFTW_PLAN_DFT_C2R_3D(Nx_,Ny_,Nz_,fdata_,data_,FFTW_ESTIMATE);
 is_set_fft_plan = true;
}

void GeneFluct3D::delete_fftw(void)
{
 if( !is_set_fft_plan ) return;
 GEN3D_FFTW_DESTROY_PLAN(pf_);
 GEN3D_FFTW_DESTROY_PLAN(pb_);
#ifdef WITH_FFTW_THREAD
 if(nthread_>0) GEN3D_FFTW_CLEANUP_THREADS();
#endif
 is_set_fft_plan = false;
}

void GeneFluct3D::check_array_alloc(void)
// To test if the table T_ is allocated
{
 if(array_allocated_) return;
 char bla[90];
 sprintf(bla,"GeneFluct3D::check_array_alloc_Error: array is not allocated");
 cout<<bla<<endl; throw ParmError(bla);
}

//-------------------------------------------------------
void GeneFluct3D::SetObservator(double redshref,double kredshref)
// The observer is at redshift z=0
//               is situated "perpendicular" to the x,y face
//               on an axis that goes through the center of the face
// It is necessary to put the cube z-axis, ie the redshifts:
//     redshref  = reference redshift
//                 If redshref<0 then redshref=0
//     kredshref = indice (double) corresponding to this redshift
//                 If kredshref<0 then kredshref=nz/2 (middle of cube)
// Example: redshref=1.5 kredshref=250.75
//    -> The pixel i=nx/2 j=ny/2 k=250.75 is at redshift 1.5
{
    if(redshref<0.) redshref = 0.;
    if(kredshref<0.) {
        if(Nz_<=0) {
            const char *bla = "GeneFluct3D::SetObservator_Error: for kredsh_ref<0 define cube geometry first";
            cout<<bla<<endl; throw ParmError(bla);
            }
        kredshref = Nz_/2.;
        }
    redsh_ref_  = redshref;
    kredsh_ref_ = kredshref;
    if(lp_>0)
        cout <<"--- GeneFluct3D::SetObservator zref="<< redsh_ref_ <<" kref="<< kredsh_ref_ << endl;
};


long GeneFluct3D::LosComRedshift(double zinc, long npoints)
// Given a position of the cube relative to the observer
// and a cosmology
// (SetObservator() and SetCosmology() should have been called !)
// This routine filled:
//   the vector "zred_" of scanned redshift (by zinc increments)
//   the vector "loscom_" of corresponding comoving distance
// -- Input:
// zinc : redshift increment for computation
// npoints : number of points required for inverting loscom -> zred
// 
{
    if(lp_>0) cout<<"--- LosComRedshift: zinc="<< zinc <<" , npoints="<< npoints <<endl;

    if(cosmo_ == NULL || redsh_ref_<0.) {
        const char *bla = "GeneFluct3D::LosComRedshift_Error: set Observator and Cosmology first";
        cout<< bla <<endl; throw ParmError(bla);
        }

    // get reference pixel quantities
    cosmo_->SetEmissionRedShift(redsh_ref_);
    // angular/luminosity/Dnu distance from reference pixel
    //dred_ref_ = Dz_/(cosmo_->Dhubble()/cosmo_->E(redsh_ref_));
    dred_ref_ = Dz_/(cosmo_->DH()/cosmo_->Ez(redsh_ref_));
    //loscom_ref_ = cosmo_->Dloscom(redsh_ref_);
    loscom_ref_ = cosmo_->LineOfSightComovDistanceMpc();
    //dtrc_ref_ = cosmo_->Dtrcom(redsh_ref_);
    dtrc_ref_ = cosmo_->TransComovDistanceMpc();
    //dlum_ref_ = cosmo_->Dlum(redsh_ref_);
    dlum_ref_ = cosmo_->LuminosityDistanceMpc();
    //dang_ref_ = cosmo_->Dang(redsh_ref_);
	dang_ref_ = cosmo_->AngularDiameterDistanceMpc();
    nu_ref_   = FREQ_21CM_HI_IN_GHZ/(1.+redsh_ref_); // GHz
    dnu_ref_  = FREQ_21CM_HI_IN_GHZ*dred_ref_/pow(1.+redsh_ref_,2.); // GHz
 
    if(lp_>0) {
        cout <<"...reference pixel redshref="<< redsh_ref_
             <<", dredref="<< dred_ref_
             <<", nuref="<< nu_ref_ <<" GHz"
             <<", dnuref="<< dnu_ref_ <<" GHz" <<endl
             <<"   dlosc="<< loscom_ref_ <<" Mpc com"
             <<", dtrc="<< dtrc_ref_ <<" Mpc com"
             <<", dlum="<< dlum_ref_ <<" Mpc"
             <<", dang="<< dang_ref_ <<" Mpc" <<endl;
        }

    // Calculate the observer coordinates at the reference point of the cube,
    // ie the origin is the center of pixel i=j=l=0.
    // The observer is on an axis centered on the middle of the Oxy face
    xobs_[0] = Nx_/2.*Dx_;
    xobs_[1] = Ny_/2.*Dy_;
    xobs_[2] = kredsh_ref_*Dz_ - loscom_ref_;

    // Is the observer in the cube?
    bool obs_in_cube = false;
    if(xobs_[2]>=0. && xobs_[2]<=Nz_*Dz_) obs_in_cube = true;

    // Find MINIMUM los com distance to the observer:
    // It's the centre of the face at k=0
    // (or zero if the observer is in the cube)
    loscom_min_ = 0.;
    if(!obs_in_cube) loscom_min_ = -xobs_[2];

    // TO BE FIXED TO BE FIXED TO BE FIXED TO BE FIXED TO BE FIXED TO BE FIXED
    if(loscom_min_<=1.e-50) {
        string emsg = "WARNING! CODE DOES NOT WORK FOR AN OBSERVER INSIDE THE CUBE";
        for(int i=0;i<50;i++)
            cout<< emsg << endl;
        }
    // TO BE FIXED TO BE FIXED TO BE FIXED TO BE FIXED TO BE FIXED TO BE FIXED


    // Find MAXIMUM los com distance to the observer:
    // where that positions the observer, the max
    // distance is one of the corners of the cube
    loscom_max_ = 0.;
    for(long i=0;i<=1;i++) {
        double dx2 = DXcom(i*(Nx_-1)); dx2 *= dx2;
        for(long j=0;j<=1;j++) {
            double dy2 = DYcom(j*(Ny_-1)); dy2 *= dy2;
            for(long k=0;k<=1;k++) {
                double dz2 = DZcom(k*(Nz_-1)); dz2 *= dz2;
                dz2 = sqrt(dx2+dy2+dz2);
                if(dz2>loscom_max_) loscom_max_ = dz2;
                }
            }
        }
        
    if(lp_>0) {
        cout <<"...zref="<< redsh_ref_ <<" kzref="<< kredsh_ref_ <<" losref="<< loscom_ref_ <<" Mpc\n"
             <<"   xobs="<< xobs_[0] <<" , "<< xobs_[1] <<" , "<< xobs_[2] <<" Mpc "
             <<" in_cube="<< obs_in_cube
             <<" loscom_min="<< loscom_min_ <<" loscom_max="<< loscom_max_ <<" Mpc (com)"<<endl;
        }

    // Fill the corresponding vectors for loscom and zred
    // Be sure to have one dlc <loscom_min and one >loscom_max
    if(zinc<=0.) zinc = 0.01;
 
    for(double z=0.; ; z+=zinc) {
        //double dlc = cosmo_->Dloscom(z);
        
        // calculate line of sight comoving distance to redshift z
	    cosmo_->SetEmissionRedShift(z);
	    double dlc = cosmo_->LineOfSightComovDistanceMpc();
	    
	    // if this distance is outside of cube (lower distance than the cube)
	    // resize the loscom(z) table to zero
        if(dlc<loscom_min_) { zred_.resize(0); loscom_.resize(0); }
        
        zred_.push_back(z);
        loscom_.push_back(dlc);
        z += zinc;
        
        // quit loop when distance is outside far side of cube
        if(dlc>loscom_max_) break; // sort after to have stored dlc>dlcmax
        }

    if(lp_>0) {
        long n = zred_.size();
        cout <<"...zred/loscom tables[zinc="<< zinc <<"]: n="<< n;
        if(n>0) cout <<" z="<< zred_[0] <<" -> d="<< loscom_[0];
        if(n>1) cout <<" , z="<< zred_[n-1] <<" -> d="<< loscom_[n-1];
        cout<<endl;
        }

    // Compute the parameters and tables needed for inversion loscom->zred
    if(npoints<3) npoints = zred_.size();
    InverseFunc invfun(zred_,loscom_);
    invfun.ComputeParab(npoints,loscom2zred_);
    loscom2zred_min_ = invfun.YMin();
    loscom2zred_max_ = invfun.YMax();

    if(lp_>0) {
        long n = loscom2zred_.size();
        cout <<"...loscom -> zred[npoints="<< npoints <<"]: n="<< n
             <<" los_min="<< loscom2zred_min_
             <<" los_max="<< loscom2zred_max_
             <<" -> zred=[";
   
        if(n>0) cout << loscom2zred_[0];
        cout <<",";
        if(n>1) cout << loscom2zred_[n-1];
        cout <<"]"<<endl;
        
        if(lp_>1 && n>0)
            for(int i=0;i<n;i++)
                if(i<2 || abs(i-n/2)<2 || i>=n-2) {
                    cout <<"    i="<< i
                         <<"  d="<< loscom2zred_min_+i*(loscom2zred_max_-loscom2zred_min_)/(n-1.)
                         <<" Mpc   z="<<loscom2zred_[i] <<endl;
                    }
        }

    return zred_.size();
};


//-------------------------------------------------------
void GeneFluct3D::WriteFits(string cfname,int bitpix)
{
 cout<<"--- GeneFluct3D::WriteFits: Writing Cube to "<<cfname<<endl;
 try {
   FitsImg3DWriter fwrt(cfname.c_str(),bitpix,5);
   fwrt.WriteKey("NX",Nx_," axe transverse 1");
   fwrt.WriteKey("NY",Ny_," axe transverse 2");
   fwrt.WriteKey("NZ",Nz_," axe longitudinal (redshift)");
   fwrt.WriteKey("DX",Dx_," Mpc");
   fwrt.WriteKey("DY",Dy_," Mpc");
   fwrt.WriteKey("DZ",Dz_," Mpc");
   fwrt.WriteKey("DKX",Dkx_," Mpc^-1");
   fwrt.WriteKey("DKY",Dky_," Mpc^-1");
   fwrt.WriteKey("DKZ",Dkz_," Mpc^-1");
   fwrt.WriteKey("ZREF",redsh_ref_," reference redshift");
   fwrt.WriteKey("KZREF",kredsh_ref_," reference redshift on z axe");
   fwrt.Write(R_);
 } catch (PThrowable & exc) {
   cout<<"Exception : "<<(string)typeid(exc).name()
       <<" - Msg= "<<exc.Msg()<<endl;
   return;
 } catch (...) {
   cout<<" some other exception was caught !"<<endl;
   return;
 }
}

void GeneFluct3D::ReadFits(string cfname)
{
 cout<<"--- GeneFluct3D::ReadFits: Reading Cube from "<<cfname<<endl;
 try {
   FitsImg3DRead fimg(cfname.c_str(),0,5);
   fimg.Read(R_);
   long nx = fimg.ReadKeyL("NX");
   long ny = fimg.ReadKeyL("NY");
   long nz = fimg.ReadKeyL("NZ");
   double dx = fimg.ReadKey("DX");
   double dy = fimg.ReadKey("DY");
   double dz = fimg.ReadKey("DZ");
   double zref = fimg.ReadKey("ZREF");
   double kzref = fimg.ReadKey("KZREF");
   setsize(nx,ny,nz,dx,dy,dz);
   setpointers(true);
   init_fftw();
   SetObservator(zref,kzref);
   array_allocated_ = true;
 } catch (PThrowable & exc) {
   cout<<"Exception : "<<(string)typeid(exc).name()
       <<" - Msg= "<<exc.Msg()<<endl;
   return;
 } catch (...) {
   cout<<" some other exception was caught !"<<endl;
   return;
 }
}

void GeneFluct3D::WritePPF(string cfname,bool write_real)
// Write either TArray<r_8> or TArray<complex <r_8> >
{
 cout<<"--- GeneFluct3D::WritePPF: Writing Cube (real="<<write_real<<") to "<<cfname<<endl;
 try {
   R_.Info()["NX"] = (int_8)Nx_;
   R_.Info()["NY"] = (int_8)Ny_;
   R_.Info()["NZ"] = (int_8)Nz_;
   R_.Info()["DX"] = (r_8)Dx_;
   R_.Info()["DY"] = (r_8)Dy_;
   R_.Info()["DZ"] = (r_8)Dz_;
   R_.Info()["ZREF"] = (r_8)redsh_ref_;
   R_.Info()["KZREF"] = (r_8)kredsh_ref_;
   POutPersist pos(cfname.c_str());
   if(write_real) pos << PPFNameTag("rgen")  << R_;
     else         pos << PPFNameTag("pkgen") << T_;
 } catch (PThrowable & exc) {
   cout<<"Exception : "<<(string)typeid(exc).name()
       <<" - Msg= "<<exc.Msg()<<endl;
   return;
 } catch (...) {
   cout<<" some other exception was caught !"<<endl;
   return;
 }
}

void GeneFluct3D::ReadPPF(string cfname)
{
 cout<<"--- GeneFluct3D::ReadPPF: Reading Cube from "<<cfname<<endl;
 try {
   bool from_real = true;
   PInPersist pis(cfname.c_str());
   string name_tag_k = "pkgen";
   bool found_tag_k = pis.GotoNameTag("pkgen");
   if(found_tag_k) {
     cout<<"           ...reading spectrum into TArray<complex <r_8> >"<<endl;
     pis >> PPFNameTag("pkgen")  >> T_;
     from_real = false;
   } else {
     cout<<"           ...reading space into TArray<r_8>"<<endl;
     pis >> PPFNameTag("rgen")  >> R_;
   }
   setpointers(from_real);  // put here to re-read DVInfo
   int_8 nx = R_.Info()["NX"];
   int_8 ny = R_.Info()["NY"];
   int_8 nz = R_.Info()["NZ"];
   r_8 dx = R_.Info()["DX"];
   r_8 dy = R_.Info()["DY"];
   r_8 dz = R_.Info()["DZ"];
   r_8 zref = R_.Info()["ZREF"];
   r_8 kzref = R_.Info()["KZREF"];
   setsize(nx,ny,nz,dx,dy,dz);
   init_fftw();
   SetObservator(zref,kzref);
   array_allocated_ = true;
 } catch (PThrowable & exc) {
   cout<<"Exception : "<<(string)typeid(exc).name()
       <<" - Msg= "<<exc.Msg()<<endl;
   return;
 } catch (...) {
   cout<<" some other exception was caught !"<<endl;
   return;
 }
}

void GeneFluct3D::WriteSlicePPF(string cfname)
// Write 3 slices of the cube according to each axis
{
 cout<<"--- GeneFluct3D::WriteSlicePPF: Writing Cube Slices "<<cfname<<endl;
 try {

   POutPersist pos(cfname.c_str());
   TMatrix<r_4> S;
   char str[16];
   long i,j,l;

   // slices in Z
   for(int s=0;s<3;s++) {
     S.ReSize(Nx_,Ny_);
     if(s==0) l=0; else if(s==1) l=(Nz_+1)/2; else  l=Nz_-1;
     sprintf(str,"z%ld",l);
     for(i=0;i<Nx_;i++) for(j=0;j<Ny_;j++) S(i,j)=data_[IndexR(i,j,l)];
     pos<<PPFNameTag(str)<<S; S.RenewObjId();
   }

   // slices in Y
   for(int s=0;s<3;s++) {
     S.ReSize(Nz_,Nx_);
     if(s==0) j=0; else if(s==1) j=(Ny_+1)/2; else  j=Ny_-1;
     sprintf(str,"y%ld",j);
     for(i=0;i<Nx_;i++) for(l=0;l<Nz_;l++) S(l,i)=data_[IndexR(i,j,l)];
     pos<<PPFNameTag(str)<<S; S.RenewObjId();
   }

   // slices in X
   for(int s=0;s<3;s++) {
     S.ReSize(Nz_,Ny_);
     if(s==0) i=0; else if(s==1) i=(Nx_+1)/2; else  i=Nx_-1;
     sprintf(str,"x%ld",i);
     for(j=0;j<Ny_;j++) for(l=0;l<Nz_;l++) S(l,j)=data_[IndexR(i,j,l)];
     pos<<PPFNameTag(str)<<S; S.RenewObjId();
   }

 } catch (PThrowable & exc) {
   cout<<"Exception : "<<(string)typeid(exc).name()
       <<" - Msg= "<<exc.Msg()<<endl;
   return;
 } catch (...) {
   cout<<" some other exception was caught !"<<endl;
   return;
 }
}

//-------------------------------------------------------
void GeneFluct3D::NTupleCheck(POutPersist &pos,string ntname,unsigned long nent)
// Fill the NTuple "ntname" with "nent" values of the cube (real or complex) and write it in "pos"
{
  if(ntname.size()<=0 || nent==0) return;
  int nvar = 0;
  if(array_type==1) nvar = 3;
  else if(array_type==2) nvar = 4;
  else return;
  const char *vname[4] = {"t","z","re","im"};
  float xnt[4];
  NTuple nt(nvar, const_cast<char **>(vname));

  if(array_type==1) {
    unsigned long nmod = Nx_*Ny_*Nz_/nent; if(nmod==0) nmod=1;
    unsigned long n=0;
    for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
      if(n==nmod) {
        int_8 ip = IndexR(i,j,l);
        xnt[0]=sqrt(i*i+j*j); xnt[1]=l; xnt[2]=data_[ip];
        nt.Fill(xnt);
        n=0;
      }
      n++;
    }
  } else {
    unsigned long nmod = Nx_*Ny_*NCz_/nent; if(nmod==0) nmod=1;
    unsigned long n=0;
    for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<NCz_;l++) {
      if(n==nmod) {
        xnt[0]=sqrt(i*i+j*j); xnt[1]=l; xnt[2]=T_(l,j,i).real(); xnt[3]=T_(l,j,i).imag();
        nt.Fill(xnt);
        n=0;
      }
      n++;
    }
  }

  pos.PutObject(nt,ntname);
}

//-------------------------------------------------------
void GeneFluct3D::Print(void)
{
 cout<<"GeneFluct3D(T_alloc="<<array_allocated_<<"):"<<endl;
 cout<<"Space Size : nx="<<Nx_<<" ny="<<Ny_<<" nz="<<Nz_<<" ("<<NTz_<<")  size="
     <<NRtot_<<endl;
 cout<<"      Resol: dx="<<Dx_<<" dy="<<Dy_<<" dz="<<Dz_<<" Mpc"
     <<", dVol="<<dVol_<<", Vol="<<Vol_<<" Mpc^3"<<endl;
 cout<<"Fourier Size : nx="<<Nx_<<" ny="<<Ny_<<" nz="<<NCz_<<endl;
 cout<<"        Resol: dkx="<<Dkx_<<" dky="<<Dky_<<" dkz="<<Dkz_<<" Mpc^-1"
     <<", Dk3="<<Dk3_<<" Mpc^-3"<<endl;
 cout<<"          (2Pi/k: "<<2.*M_PI/Dkx_<<" "<<2.*M_PI/Dky_<<" "<<2.*M_PI/Dkz_<<" Mpc)"<<endl;
 cout<<"      Nyquist: kx="<<Knyqx_<<" ky="<<Knyqy_<<" kz="<<Knyqz_<<" Mpc^-1"
     <<", Kmax="<<GetKmax()<<" Mpc^-1"<<endl;
 cout<<"          (2Pi/k: "<<2.*M_PI/Knyqx_<<" "<<2.*M_PI/Knyqy_<<" "<<2.*M_PI/Knyqz_<<" Mpc)"<<endl;
 cout<<"Redshift "<<redsh_ref_<<" for z axe at k="<<kredsh_ref_<<endl;
}

//-------------------------------------------------------
void GeneFluct3D::ComputeFourier0(ClassFunc1D& pk_at_z)
// cf ComputeFourier() but with another method of generating the spectrum
//    (attention make one fft to generate the spectrum)
{

 // --- realization of a table Gaussian drawing
 if(lp_>0) cout<<"--- ComputeFourier0: before gaussian filling ---"<<endl;
 // On tient compte du pb de normalisation de FFTW3
 double sntot = sqrt((double)NRtot_);
 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   data_[ip] = NorRand()/sntot;
 }

 // --- realization of a table Gaussian drawing
 if(lp_>0) cout<<"...before fft real ---"<<endl;
 GEN3D_FFTW_EXECUTE(pf_);

 // --- Fill with one realization
 if(lp_>0) cout<<"...before Fourier realization filling"<<endl;
 T_(0,0,0) = complex<GEN3D_TYPE>(0.);  // initialize
 long lmod = Nx_/10; if(lmod<1) lmod=1;
 for(long i=0;i<Nx_;i++) {
   long ii = (i>Nx_/2) ? Nx_-i : i;
   double kx = ii*Dkx_;  kx *= kx;
   if(lp_>0 && i%lmod==0) cout<<"i="<<i<<" ii="<<ii<<endl;
   for(long j=0;j<Ny_;j++) {
     long jj = (j>Ny_/2) ? Ny_-j : j;
     double ky = jj*Dky_;  ky *= ky;
     for(long l=0;l<NCz_;l++) {
       double kz = l*Dkz_;  kz *= kz;
       if(i==0 && j==0 && l==0) continue; 
       double k = sqrt(kx+ky+kz);
       // cf normalisation: Peacock, Cosmology, formula 16.38 p504
       double pk = pk_at_z(k)/Vol_;
       // here no "/2" because of the remark below
       T_(l,j,i) *= sqrt(pk);
     }
   }
 }

 array_type = 2;

 if(lp_>0) cout<<"...computing power"<<endl;
 double p = compute_power_carte();
 if(lp_>0) cout<<"Puissance dans la realisation: "<<p<<endl;

}

//-------------------------------------------------------
void GeneFluct3D::ComputeFourier(ClassFunc1D& pk_at_z)
// Calculate a realisation of spectrum "pk_at_z"
// Attention: in TArray first index varies fastest
// Normalisation explanation: see Coles & Lucchin, Cosmology, p264-265
// FFTW3: note N=Nx*Ny*Nz
// f  --(FFT)-->  F = TF(f)  --(FFT^-1)-->  fb = TF^-1(F) = TF^-1(TF(f))
// sum(f(x_i)^2) = S
//                sum(F(nu_i)^2) = S*N
//                                          sum(fb(x_i)^2) = S*N^2
{
 // --- RaZ of table
 T_ = complex<GEN3D_TYPE>(0.);

 // --- Fill with a realization
 if(lp_>0) cout<<"--- ComputeFourier ---"<<endl;
 long lmod = Nx_/10; if(lmod<1) lmod=1;
 for(long i=0;i<Nx_;i++) 
	{
	long ii = (i>Nx_/2) ? Nx_-i : i;
	double kx = ii*Dkx_;  kx *= kx;
	if(lp_>0 && i%lmod==0) cout<<"i="<<i<<" ii="<<ii<<endl;
	for(long j=0;j<Ny_;j++) 
		{
		long jj = (j>Ny_/2) ? Ny_-j : j;
		double ky = jj*Dky_;  ky *= ky;
		for(long l=0;l<NCz_;l++) 
			{
			double kz = l*Dkz_;  kz *= kz;
			if(i==0 && j==0 && l==0) continue;
			double k = sqrt(kx+ky+kz);
			// cf normalisation: Peacock, Cosmology, formula 16.38 p504
			double pk = pk_at_z(k)/Vol_;
			// Explanation of the /2: see perandom.cc
			// or equally Coles & Lucchin, Cosmology formula 13.7.2 p279
			T_(l,j,i) = ComplexGaussianRand(sqrt(pk/2.));
			}
		}
	}

 array_type = 2;
 manage_coefficients();   // big effect for the spectra that use it !

 if(lp_>0) cout<<"...computing power"<<endl;
 double p = compute_power_carte();
 if(lp_>0) cout<<"Power in the realisation: "<<p<<endl;

}

long GeneFluct3D::manage_coefficients(void)
// Take into account the real and complexe conjugate coefficients
// because we want a realization of a real data in real space
{
 if(lp_>1) cout<<"...managing coefficients"<<endl;
 check_array_alloc();

 // 1./ Le Continu and Nyquist are real
 long nreal = 0;
 for(long kk=0;kk<2;kk++) {
   long k=0;  // continu
   if(kk==1) {if(Nz_%2!=0) continue; else k = Nz_/2;}  // Nyquist
   for(long jj=0;jj<2;jj++) {
     long j=0;
     if(jj==1) {if( Ny_%2!=0) continue; else j = Ny_/2;}
     for(long ii=0;ii<2;ii++) {
       long i=0;
       if(ii==1) {if( Nx_%2!=0) continue; else i = Nx_/2;}
       int_8 ip = IndexC(i,j,k);
       //cout<<"i="<<i<<" j="<<j<<" k="<<k<<" = ("<<fdata_[ip][0]<<","<<fdata_[ip][1]<<")"<<endl;
       fdata_[ip][1] = 0.; fdata_[ip][0] *= M_SQRT2;
       nreal++;
     }
   }
 }
 if(lp_>1) cout<<"Number of forced real number ="<<nreal<<endl;

 // 2./ The elements complex conjugates (all in the plan k=0,Nyquist)

 // a./ rows and columns of continu and nyquist
 long nconj1 = 0;
 for(long kk=0;kk<2;kk++) {
   long k=0;  // continu
   if(kk==1) {if(Nz_%2!=0) continue; else k = Nz_/2;}  // Nyquist
   for(long jj=0;jj<2;jj++) { // according to j
     long j=0;
     if(jj==1) {if( Ny_%2!=0) continue; else j = Ny_/2;}
     for(long i=1;i<(Nx_+1)/2;i++) {
       int_8 ip = IndexC(i,j,k);
       int_8 ip1 = IndexC(Nx_-i,j,k);
       fdata_[ip1][0] = fdata_[ip][0]; fdata_[ip1][1] = -fdata_[ip][1];
       nconj1++;
     }
   }
   for(long ii=0;ii<2;ii++) {
     long i=0;
     if(ii==1) {if( Nx_%2!=0) continue; else i = Nx_/2;}
     for(long j=1;j<(Ny_+1)/2;j++) {
       int_8 ip = IndexC(i,j,k);
       int_8 ip1 = IndexC(i,Ny_-j,k);
       fdata_[ip1][0] = fdata_[ip][0]; fdata_[ip1][1] = -fdata_[ip][1];
       nconj1++;
     }
   }
 }
 if(lp_>1) cout<<"Number of forced conjugate on cont+nyq ="<<nconj1<<endl;

 // b./ lines and columns outside continu and nyquist
 long nconj2 = 0;
 for(long kk=0;kk<2;kk++) {
   long k=0;  // continu
   if(kk==1) {if(Nz_%2!=0) continue; else k = Nz_/2;}  // Nyquist
   for(long j=1;j<(Ny_+1)/2;j++) {
     if(Ny_%2==0 && j==Ny_/2) continue; // don't restate nyquist in j
     for(long i=1;i<Nx_;i++) {
       if(Nx_%2==0 && i==Nx_/2) continue; // don't restate nyquist in i
       int_8 ip = IndexC(i,j,k);
       int_8 ip1 = IndexC(Nx_-i,Ny_-j,k);
       fdata_[ip1][0] = fdata_[ip][0]; fdata_[ip1][1] = -fdata_[ip][1];
       nconj2++;
     }
   }
 }
 if(lp_>1) cout<<"Number of forced conjugate hors cont+nyq ="<<nconj2<<endl;

 if(lp_>1) cout<<"Check: ddl= "<<NRtot_<<" =?= "<<2*(Nx_*Ny_*NCz_-nconj1-nconj2)-8<<endl;

 return nreal+nconj1+nconj2;
}

double GeneFluct3D::compute_power_carte(void)
// Calculation of the power of the realisation of spectrum Pk
{
 check_array_alloc();

 double s2 = 0.;
 for(long l=0;l<NCz_;l++)
   for(long j=0;j<Ny_;j++)
     for(long i=0;i<Nx_;i++) s2 += MODULE2(T_(l,j,i));

 double s20 = 0.;
 for(long j=0;j<Ny_;j++)
   for(long i=0;i<Nx_;i++) s20 += MODULE2(T_(0,j,i));

 double s2n = 0.;
 if(Nz_%2==0)
   for(long j=0;j<Ny_;j++)
     for(long i=0;i<Nx_;i++) s2n += MODULE2(T_(NCz_-1,j,i));

 return 2.*s2 -s20 -s2n;
}

//-------------------------------------------------------------------
void GeneFluct3D::FilterByPixel(void)
// Filtering by pixel window function (parallelepipede)
// TF = 1/(dx*dy*dz)*Int[{-dx/2,dx/2},{-dy/2,dy/2},{-dz/2,dz/2}]
//                   e^(ik_x*x) e^(ik_y*y) e^(ik_z*z) dxdydz
//    = 2/(k_x*dx) * sin(k_x*dx/2)  * (idem y) * (idem z)
// Manage divergence at 0: sin(y)/y = 1 - y^2/6*(1-y^2/20)
//                          with y = k_x*dx/2
{
 if(lp_>0) cout<<"--- FilterByPixel ---"<<endl;
 check_array_alloc();

 for(long i=0;i<Nx_;i++) {
   long ii = (i>Nx_/2) ? Nx_-i : i;
   double kx = ii*Dkx_ *Dx_/2;
   double pk_x = pixelfilter(kx);
   for(long j=0;j<Ny_;j++) {
     long jj = (j>Ny_/2) ? Ny_-j : j;
     double ky = jj*Dky_ *Dy_/2;
     double pk_y = pixelfilter(ky);
     for(long l=0;l<NCz_;l++) {
       double kz = l*Dkz_ *Dz_/2;
       double pk_z =  pixelfilter(kz);
       T_(l,j,i) *= pk_x*pk_y*pk_z;
     }
   }
 }

}

//-------------------------------------------------------------------
void GeneFluct3D::ApplyGrowthFactor(int type_evol)
// Apply Growth to real space
// Using the correspondance between redshift and los comoving distance
// describe in vector "zred_" "loscom_"
// type_evol = 1 : evolution with distance to the observer
//             2 : evolution with distance to plane Z
//             (all pixels in plane Z are put at the same redshift z, as that in middle)
{
 if(lp_>0) cout<<"--- ApplyGrowthFactor: evol="<<type_evol<<endl;
 check_array_alloc();

 if(growth_ == NULL) {
   const char *bla = "GeneFluct3D::ApplyGrowthFactor_Error: set GrowthFactor first";
   cout<<bla<<endl; throw ParmError(bla);
 }
 if(type_evol<1 || type_evol>2) {
   const char *bla = "GeneFluct3D::ApplyGrowthFactor_Error: bad type_evol value";
   cout<<bla<<endl; throw ParmError(bla);
 }

 InterpFunc interpinv(loscom2zred_min_,loscom2zred_max_,loscom2zred_);
 //unsigned short ok;

 //CHECK: Histo hgr(0.9*zred_[0],1.1*zred_[n-1],1000);
 for(long i=0;i<Nx_;i++) {
   double dx2 = DXcom(i); dx2 *= dx2;
   for(long j=0;j<Ny_;j++) {
     double dy2 = DYcom(j); dy2 *= dy2;
     for(long l=0;l<Nz_;l++) {
       double dz = DZcom(l);
       if(type_evol==1) dz = sqrt(dx2+dy2+dz*dz);
         else dz = fabs(dz); // all the Z planes at same redshift
       double z = interpinv(dz);
       //CHECK: hgr.Add(z);
       double dzgr = (*growth_)(z);   // interpolation by parts
       //double dzgr = growth_->Linear(z,ok);  // linear interpolation
       //double dzgr = growth_->Parab(z,ok);  // parabolic interpolation
       int_8 ip = IndexR(i,j,l);
       data_[ip] *= dzgr;
     }
   }
 }

 //CHECK: {POutPersist pos("applygrowth.ppf"); string tag="hgr"; pos.PutObject(hgr,tag);}

}

//-------------------------------------------------------------------
void GeneFluct3D::ComputeReal(void)
// Calculate a realisation in real space
{
 if(lp_>0) cout<<"--- ComputeReal ---"<<endl;
 check_array_alloc();

 // On fait la FFT
 GEN3D_FFTW_EXECUTE(pb_);
 array_type = 1;
}

//-------------------------------------------------------------------
void GeneFluct3D::ReComputeFourier(void)
{
 if(lp_>0) cout<<"--- ReComputeFourier ---"<<endl;
 check_array_alloc();

 // Do the FFT
 GEN3D_FFTW_EXECUTE(pf_);
 array_type = 2;

 // Correct the problem of the normalisation of FFTW3
 complex<r_8> v((r_8)NRtot_,0.);
 for(long i=0;i<Nx_;i++)
   for(long j=0;j<Ny_;j++)
     for(long l=0;l<NCz_;l++) T_(l,j,i) /= v;
}

//-------------------------------------------------------------------
int GeneFluct3D::ComputeSpectrum(HistoErr& herr)
// Compute spectrum from "T" and fill HistoErr "herr"
// T : in the standard format of GeneFuct3D: T(nz,ny,nx)
// ie T(kz,ky,kx) with  0<kz<kz_nyq  -ky_nyq<ky<ky_nyq  -kx_nyq<kx<kx_nyq
{
 if(lp_>0) cout<<"--- ComputeSpectrum changed?---"<<endl;
 check_array_alloc();

 if(herr.NBins()<0) return -1;
 herr.Zero();

 // Attention on the order
	for(long i=0;i<Nx_;i++) 
		{
		long ii = (i>Nx_/2) ? Nx_-i : i;
		double kx = ii*Dkx_;  
		kx *= kx;
		for(long j=0;j<Ny_;j++) 
			{
			long jj = (j>Ny_/2) ? Ny_-j : j;
			double ky = jj*Dky_;  
			ky *= ky;
			for(long l=0;l<NCz_;l++) 
				{
				double kz = l*Dkz_;
				double k = sqrt(kx+ky+kz*kz);
				double pk = MODULE2(T_(l,j,i));
				herr.Add(k,pk);
				}
			}	
		}
	herr.ToVariance();

 // renormalize to directly compare to original spectrum
 double norm = Vol_;
 herr *= norm;
 cout <<"print spec norm: "<<norm<<endl;
 cout <<"print (dkx,dky,dkz)=("<<Dkx_<<","<<Dky_<<","<<Dkz_<<")"<<endl;
 cout <<"print (Nx,Ny,NCz)=("<<Nx_<<","<<Ny_<<","<<NCz_<<")"<<endl;
 return 0;
}

int GeneFluct3D::ComputeSpectrum2D(Histo2DErr& herr)
{
 if(lp_>0) cout<<"--- ComputeSpectrum2D ---"<<endl;
 check_array_alloc();

 if(herr.NBinX()<0 || herr.NBinY()<0) return -1;
 herr.Zero();

 // Attention  to the order
 for(long i=0;i<Nx_;i++) {
   long ii = (i>Nx_/2) ? Nx_-i : i;
   double kx = ii*Dkx_;  kx *= kx;
   for(long j=0;j<Ny_;j++) {
     long jj = (j>Ny_/2) ? Ny_-j : j;
     double ky = jj*Dky_;  ky *= ky;
     double kt = sqrt(kx+ky);
     for(long l=0;l<NCz_;l++) {
       double kz = l*Dkz_;
       double pk = MODULE2(T_(l,j,i));
       herr.Add(kt,kz,pk);
     }
   }
 }
 herr.ToVariance();

 // renormalize to directly compare to original spectrum
 double norm = Vol_;
 herr *= norm;

 return 0;
}

//-------------------------------------------------------------------
int GeneFluct3D::ComputeSpectrum(HistoErr& herr,double sigma,bool pixcor)
// Compute spectrum from "T" and fill HistoErr "herr"
// WITH the substraction of noise and the correction by filterpixel()
// If this isn't done, end up with a non-isotropic spectrum!
//
// T : in the standard format of GeneFuct3D: T(nz,ny,nx)
// ie T(kz,ky,kx) with  0<kz<kz_nyq  -ky_nyq<ky<ky_nyq  -kx_nyq<kx<kx_nyq
{
  if(lp_>0) cout<<"--- ComputeSpectrum: sigma="<<sigma<<endl;
 check_array_alloc();

 if(sigma<=0.) sigma = 0.;
 double sigma2 = sigma*sigma / (double)NRtot_;

 if(herr.NBins()<0) return -1;
 herr.Zero();

 TVector<r_8> vfz(NCz_);
 if(pixcor)  // kz = l*Dkz_
   for(long l=0;l<NCz_;l++) {vfz(l)=pixelfilter(l*Dkz_ *Dz_/2); vfz(l)*=vfz(l);}

 // Attention to the order
 for(long i=0;i<Nx_;i++) {
   long ii = (i>Nx_/2) ? Nx_-i : i;
   double kx = ii*Dkx_;
   double fx = (pixcor) ? pixelfilter(kx*Dx_/2): 1.;
   kx *= kx; fx *= fx;
   for(long j=0;j<Ny_;j++) {
     long jj = (j>Ny_/2) ? Ny_-j : j;
     double ky = jj*Dky_;
     double fy = (pixcor) ? pixelfilter(ky*Dy_/2): 1.;
     ky *= ky; fy *= fy;
     for(long l=0;l<NCz_;l++) {
       double kz = l*Dkz_;
       double k = sqrt(kx+ky+kz*kz);
       double pk = MODULE2(T_(l,j,i)) - sigma2;
       double fz = (pixcor) ? vfz(l): 1.;
       double f = fx*fy*fz;
       if(f>0.) herr.Add(k,pk/f);
     }
   }
 }
 herr.ToVariance();
 for(int i=0;i<herr.NBins();i++) herr(i) += sigma2;

 // renormalize to directly compare to original spectrum
 double norm = Vol_;
 herr *= norm;

 return 0;
}

int GeneFluct3D::ComputeSpectrum2D(Histo2DErr& herr,double sigma,bool pixcor)
// WTIH the substraction of noise and the correction by filterpixel()
{
 if(lp_>0) cout<<"--- ComputeSpectrum2D: sigma="<<sigma<<endl;
 check_array_alloc();

 if(sigma<=0.) sigma = 0.;
 double sigma2 = sigma*sigma / (double)NRtot_;

 if(herr.NBinX()<0 || herr.NBinY()<0) return -1;
 herr.Zero();

 TVector<r_8> vfz(NCz_);
 if(pixcor)  // kz = l*Dkz_
   for(long l=0;l<NCz_;l++) {vfz(l)=pixelfilter(l*Dkz_ *Dz_/2); vfz(l)*=vfz(l);}

 // Attention to the order
 for(long i=0;i<Nx_;i++) {
   long ii = (i>Nx_/2) ? Nx_-i : i;
   double kx = ii*Dkx_;
   double fx = (pixcor) ? pixelfilter(kx*Dx_/2) : 1.;
   kx *= kx; fx *= fx;
   for(long j=0;j<Ny_;j++) {
     long jj = (j>Ny_/2) ? Ny_-j : j;
     double ky = jj*Dky_;
     double fy = (pixcor) ? pixelfilter(ky*Dy_/2) : 1.;
     ky *= ky; fy *= fy;
     double kt = sqrt(kx+ky);
     for(long l=0;l<NCz_;l++) {
       double kz = l*Dkz_;
       double pk = MODULE2(T_(l,j,i)) - sigma2;
       double fz = (pixcor) ? vfz(l): 1.;
       double f = fx*fy*fz;
       if(f>0.) herr.Add(kt,kz,pk/f);
     }
   }
 }
 herr.ToVariance();
 for(int i=0;i<herr.NBinX();i++)
    for(int j=0;j<herr.NBinY();j++) herr(i,j) += sigma2;

 // renormalize to directly compare to original spectrum
 double norm = Vol_;
 herr *= norm;

 return 0;
}

//-------------------------------------------------------
int_8 GeneFluct3D::VarianceFrReal(double R,double& var)
// Recompute MASS variance in spherical top-hat (radius=R)
// By definition: SigmaR^2 = <(M-<M>)^2>/<M>^2
//                 where M = mass in a sphere of radius R
// --- ATTENTION: the variance calculated has a very big dispersion
//  (especially all if the volume of the cube is small). To check 
//  that the sigmaR calculed  by this method agrees with the
//  sigmaR input, must do many simulations (~100) and look at the 
//  mean sigmaR recovered
{
  if(lp_>0) cout<<"--- VarianceFrReal R="<<R<<endl;
 check_array_alloc();

 long dnx = long(R/Dx_)+1; if(dnx<=0) dnx = 1;
 long dny = long(R/Dy_)+1; if(dny<=0) dny = 1;
 long dnz = long(R/Dz_)+1; if(dnz<=0) dnz = 1;
 if(lp_>0) cout<<"dnx="<<dnx<<" dny="<<dny<<" dnz="<<dnz<<endl;

 double sum=0., sum2=0., sn=0., r2 = R*R;
 int_8 nsum=0;

 for(long i=dnx;i<Nx_-dnx;i+=2*dnx) {
   for(long j=dny;j<Ny_-dny;j+=2*dny) {
     for(long l=dnz;l<Nz_-dnz;l+=2*dnz) {
       double m=0.; int_8 n=0;
       for(long ii=i-dnx;ii<=i+dnx;ii++) {
         double x = (ii-i)*Dx_; x *= x;
         for(long jj=j-dny;jj<=j+dny;jj++) {
           double y = (jj-j)*Dy_; y *= y;
           for(long ll=l-dnz;ll<=l+dnz;ll++) {
             double z = (ll-l)*Dz_; z *= z;
             if(x+y+z>r2) continue;
             int_8 ip = IndexR(ii,jj,ll);
             m += 1.+data_[ip];  // 1+drho/rho
             n++;
	   }
	 }
       }
       if(n>0) {sum += m; sum2 += m*m; nsum++; sn += n;}
       //cout<<i<<","<<j<<","<<l<<" n="<<n<<" m="<<m<<" sum="<<sum<<" sum2="<<sum2<<endl;
     }
   }
 }

 if(nsum<=1) {var=0.; return nsum;}
 sum /= nsum;
 sum2 = sum2/nsum - sum*sum;
 sn /= nsum;
 if(lp_>0) cout<<"...<n>="<<sn<<", nsum="<<nsum<<" <M>="<<sum<<" <(M-<M>)^2>="<<sum2<<endl;
 var = sum2/(sum*sum);  // <dM^2>/<M>^2
 if(lp_>0) cout<<"...sigmaR^2 = <(M-<M>)^2>/<M>^2 = "<<var
               <<" -> sigmaR = "<<sqrt(var)<<endl;

 return nsum;
}

//-------------------------------------------------------
int_8 GeneFluct3D::NumberOfBad(double vmin,double vmax)
// number of pixels outside of ]vmin,vmax[ extremities excluded
//     ->  vmin and vmax are considered as bad
{
 check_array_alloc();

 int_8 nbad = 0;
 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   double v = data_[ip];
   if(v<=vmin || v>=vmax) nbad++;
 }

 if(lp_>0) cout<<"--- NumberOfBad "<<nbad<<" px out of ]"<<vmin<<","<<vmax
               <<"[ i.e. frac="<<nbad/(double)NRtot_<<endl;
 return nbad;
}

int_8 GeneFluct3D::MinMax(double& xmin,double& xmax,double vmin,double vmax)
// Calculation of the values xmin et xmax in the real cube with values ]vmin,vmax[ extremities excluded
{
 bool tstval = (vmax>vmin)? true: false;
 if(lp_>0) {
   cout<<"--- MinMax";
   if(tstval) cout<<"  range=]"<<vmin<<","<<vmax<<"[";
   cout<<endl;
 }
 check_array_alloc();

 int_8 n = 0;
 xmin = xmax = data_[0];

 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   double x = data_[ip];
   if(tstval && (x<=vmin || x>=vmax)) continue;
   if(x<xmin) xmin = x;
   if(x>xmax) xmax = x;
   n++;
 }

 if(lp_>0) cout<<"  n="<<n<<" min="<<xmin<<" max="<<xmax<<endl;

 return n;
}

int_8 GeneFluct3D::MeanSigma2(double& rm,double& rs2,double vmin,double vmax
                             ,bool useout,double vout)
// Calculation of mean,sigma2 in the real cube with values ]vmin,vmax[ extremities excluded
// useout = false: don't use the pixels outside the limits to calculate mean,sigma2
//          true : use the pixels outside the limits to calculate mean,sigma2
//                 replace their values by "vout"
{
 bool tstval = (vmax>vmin)? true: false;
 if(lp_>0) {
   cout<<"--- MeanSigma2";
   if(tstval) cout<<"  range=]"<<vmin<<","<<vmax<<"[";
   if(useout) cout<<", useout="<<useout<<" vout="<<vout;
   cout<<endl;
 }
 check_array_alloc();

 int_8 n = 0;
 rm = rs2 = 0.;

 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   double v = data_[ip];
   if(tstval) {
     if(v<=vmin || v>=vmax) {if(useout) v=vout; else continue;}
   }
   rm += v;
   rs2 += v*v;
   n++;
 }

 if(n>1) {
   rm /= (double)n;
   rs2 = rs2/(double)n - rm*rm;
 }

 if(lp_>0) cout<<"  n="<<n<<" m="<<rm<<" s2="<<rs2<<" s="<<sqrt(fabs(rs2))<<endl;

 return n;
}

int_8 GeneFluct3D::SetToVal(double vmin, double vmax,double val0)
// set to "val0" if out of range ]vmin,vmax[ extremities excluded
// ie set to "val0" if in [vmin,vmax] -> vmin and vmax are set to val0
{
 check_array_alloc();

 int_8 nbad = 0;
 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   double v = data_[ip];
   if(v<=vmin || v>=vmax) {data_[ip] = val0; nbad++;}
 }

 if(lp_>0) cout<<"--- SetToVal "<<nbad<<" px set to="<<val0
               <<" because out of range=]"<<vmin<<","<<vmax<<"["<<endl;
 return nbad;
}

void GeneFluct3D::ScaleOffset(double scalecube,double offsetcube)
// Replace  "V"  by  "scalecube * ( V + offsetcube )"
{
 if(lp_>0) cout<<"--- ScaleCube scale="<<scalecube<<" offset="<<offsetcube<<endl;

 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   data_[ip] = scalecube * ( data_[ip] + offsetcube );
 }

 return;
}

//-------------------------------------------------------
void GeneFluct3D::TurnFluct2Mass(void)
// d_rho/rho -> Mass  (add one!)
{
 if(lp_>0) cout<<"--- TurnFluct2Mass ---"<<endl;
 check_array_alloc();


 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   data_[ip] += 1.;
 }
}

double GeneFluct3D::TurnFluct2MeanNumber(double val_by_mpc3)
// ATTENTION: managment of pixels<0 proposed here induces a loss of variance,
// the spectrum Pk reconstructed will be lower!
// The effect will be greater if the number of pixels<0 is big
{
 if(lp_>0) cout<<"--- TurnFluct2MeanNumber : "<<val_by_mpc3<<" quantity (gal or mass)/Mpc^3"<<endl;

 // First convert dRho/Rho into 1+dRho/Rho
 int_8  nball = 0; double sumall = 0., sumall2 = 0.;
 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   data_[ip] += 1.;
   nball++; sumall += data_[ip]; sumall2 += data_[ip]*data_[ip];
 }
 if(nball>2) {
   sumall /= (double)nball;
   sumall2 = sumall2/(double)nball - sumall*sumall;
   if(lp_>0) cout<<"1+dRho/Rho: mean="<<sumall<<" variance="<<sumall2
                 <<" -> "<<sqrt(fabs(sumall2))<<endl;
 }

 // Find contribution for positive pixels
 int_8  nbpos = 0; double sumpos = 0. , sumpos2 = 0.;
 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   double v = data_[ip];
   if(data_[ip]>0.) {nbpos++; sumpos += v; sumpos2 += v*v;}
 }
 if(nbpos<1) {
   cout<<"TurnFluct2MeanNumber_Error: nbpos<1"<<endl;
   throw RangeCheckError("TurnFluct2MeanNumber_Error: nbpos<1");
 }
 sumpos2 = sumpos2/nball - sumpos*sumpos/(nball*nball);
 if(lp_>0)
   cout<<"1+dRho/Rho with v<0 set to zero: mean="<<sumpos/nball
       <<" variance="<<sumpos2<<" -> "<<sqrt(fabs(sumpos2))<<endl;
   cout<<"Sum of positive values: sumpos="<<sumpos
       <<" (n(v>0) = "<<nbpos<<" frac(v>0)="<<nbpos/(double)NRtot_<<")"<<endl;

 // - Put exactly val_by_mpc3*Vol galaxies (or Msol) in our survey
 // - Only in the pixels with mass >0.
 // - Put to zero the pixels <0
 double dn = val_by_mpc3 * Vol_ / sumpos;
 if(lp_>0) cout<<"...density move from "
               <<val_by_mpc3*dVol_<<" to "<<dn<<" / pixel"<<endl;

 double sum = 0.;
 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   if(data_[ip]<=0.) data_[ip] = 0.;
   else {
     data_[ip] *= dn;
     sum += data_[ip];
   }
 }

 if(lp_>0) cout<<"...quantity put into survey "<<sum<<" / "<<val_by_mpc3*Vol_<<endl;

 return sum;
}

double GeneFluct3D::ApplyPoisson(void)
// do NOT treate negative or nul mass  -> leave it as it is
{
 if(lp_>0) cout<<"--- ApplyPoisson ---"<<endl;
 check_array_alloc();

 double sum = 0.;
 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   double v = data_[ip];
   if(v>0.) {
     uint_8 dn = PoissonRand(v,10.);
     data_[ip] = (double)dn;
     sum += (double)dn;
   }
 }
 if(lp_>0) cout<<sum<<" galaxies put into survey"<<endl;

 return sum;
}

double GeneFluct3D::TurnNGal2Mass(FunRan& massdist,bool axeslog)
// do NOT treat negative or nul mass  -> let it as it is
// INPUT:
//   massdist : distribution of mass (m*dn/dm)
//   axeslog = false : return the mass
//           = true  : return log10(mass)
// RETURN the total mass
{
 if(lp_>0) cout<<"--- TurnNGal2Mass ---"<<endl;
 check_array_alloc();

 double sum = 0.;
 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   double v = data_[ip];
   if(v>0.) {
     long ngal = long(v+0.1);
     data_[ip] = 0.;
     for(long i=0;i<ngal;i++) {
       double m = massdist.RandomInterp();  // massdist.Random();
       if(axeslog) m = pow(10.,m);
       data_[ip] += m;
     }
     sum += data_[ip];
   }
 }
 if(lp_>0) cout<<sum<<" MSol HI mass put into survey"<<endl;

 return sum;
}

double GeneFluct3D::TurnNGal2MassQuick(SchechterMassDist& schmdist)
// identical to TurnNGal2Mass but a lot faster
{
 if(lp_>0) cout<<"--- TurnNGal2MassQuick ---"<<endl;
 check_array_alloc();

 double sum = 0.;
 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) for(long l=0;l<Nz_;l++) {
   int_8 ip = IndexR(i,j,l);
   double v = data_[ip];
   if(v>0.) {
     long ngal = long(v+0.1);
     data_[ip] = schmdist.TirMass(ngal);
     sum += data_[ip];
   }
 }
 if(lp_>0) cout<<sum<<" MSol HI mass put into survey"<<endl;

 return sum;
}

void GeneFluct3D::AddNoise2Real(double snoise,int type_evol)
// add noise to every pixels (meme les <=0 !)
// type_evol = 0 : no evolution of noise power spectrum
//             1 : evolution of noise power spectrum with the distance to the observer
//             2 : evolution of noise power spectrum with the distance to plane Z
//                 (all the Z plane are at the same redshift as their center)
{
 if(lp_>0) cout<<"--- AddNoise2Real: snoise = "<<snoise<<" evol="<<type_evol<<endl;
 check_array_alloc();

 if(type_evol<0) type_evol = 0;
 if(type_evol>2) {
   const char *bla = "GeneFluct3D::AddNoise2Real_Error: bad type_evol value";
   cout<<bla<<endl; throw ParmError(bla);
 }

 vector<double> correction;
 InterpFunc *intercor = NULL;

 if(type_evol>0) {
   // Sigma_Noise(en mass) :
   //      Slim ~ 1/sqrt(DNu) * sqrt(nlobe)   en W/m^2Hz
   //      Flim ~ sqrt(DNu) * sqrt(nlobe)   en W/m^2
   //      Mlim ~ sqrt(DNu) * (Dlum)^2 * sqrt(nlobe)  en Msol
   //                                    nlobe ~ 1/Dtrcom^2
   //      Mlim ~ sqrt(DNu) * (Dlum)^2 / Dtrcom
   if( (cosmo_ == NULL) || (redsh_ref_<0.| loscom2zred_.size()<1 )) {
     const char *bla = "GeneFluct3D::AddNoise2Real_Error: set Observator and Cosmology first";
     cout<<bla<<endl; throw ParmError(bla);
   }
   InterpFunc interpinv(loscom2zred_min_,loscom2zred_max_,loscom2zred_);
   long nsz = loscom2zred_.size(), nszmod=((nsz>10)? nsz/10: 1);
   for(long i=0;i<nsz;i++) {
     double d = interpinv.X(i);
     double zred = interpinv(d);
		 cosmo_->SetEmissionRedShift(zred);
     //double dtrc = cosmo_->Dtrcom(zred);  // for varying the solid angle
		 double dtrc = cosmo_->TransComovDistanceMpc();  // for varying the solid angle
     //double dlum = cosmo_->Dlum(zred);  // for varying the conversion to mass HI
     double dlum = cosmo_->LuminosityDistanceMpc();
		 //double dred = Dz_/(cosmo_->Dhubble()/cosmo_->E(zred));
		 double dred = Dz_/(cosmo_->DH()/cosmo_->Ez(zred));
     double dnu  = FREQ_21CM_HI_IN_GHZ*dred/pow(1.+zred,2.); // to vary dNu
     double corr = sqrt(dnu/dnu_ref_) * pow(dlum/dlum_ref_,2.) * dtrc_ref_/dtrc;
     if(lp_>0 && (i==0 || i==nsz-1 || i%nszmod==0))
       cout<<"i="<<i<<" d="<<d<<" red="<<zred<<" dred="<<dred<<" dnu="<<dnu
           <<" dtrc="<<dtrc<<" dlum="<<dlum<<" -> cor="<<corr<<endl;
     correction.push_back(corr);
   }
   intercor = new InterpFunc(loscom2zred_min_,loscom2zred_max_,correction);
 } 

 double corrlim[2] = {1.,1.};
 for(long i=0;i<Nx_;i++) {
   double dx2 = DXcom(i); dx2 *= dx2;
   for(long j=0;j<Ny_;j++) {
     double dy2 = DYcom(j); dy2 *= dy2;
     for(long l=0;l<Nz_;l++) {
       double corr = 1.;
       if(type_evol>0) {
         double dz = DZcom(l);
         if(type_evol==1) dz = sqrt(dx2+dy2+dz*dz);
           else dz = fabs(dz); // all the Z planes are at the same redshift
         corr = (*intercor)(dz);
         if(corr<corrlim[0]) corrlim[0]=corr; else if(corr>corrlim[1]) corrlim[1]=corr;
       }
       int_8 ip = IndexR(i,j,l);
       data_[ip] += snoise*corr*NorRand();
     }
   }
 }
 if(type_evol>0)
   cout<<"correction factor range: ["<<corrlim[0]<<","<<corrlim[1]<<"]"<<endl;

 if(intercor!=NULL) delete intercor;
}

}  // End namespace SOPHYA




/*********************************************************************
void GeneFluct3D::AddAGN(double lfjy,double lsigma,double powlaw)
// Add AGN flux into simulation:
// --- Procedure:
// 1. lancer "cmvdefsurv" avec les parametres du survey
//        (au redshift de reference du survey)
//    et recuperer l'angle solide "angsol sr" du pixel elementaire
//    au centre du cube.
// 2. lancer "cmvtstagn" pour cet angle solide -> cmvtstagn.ppf
// 3. regarder l'histo "hlfang" et en deduire un equivalent gaussienne
//    cad une moyenne <log10(S)> et un sigma "sig"
//    Attention: la distribution n'est pas gaussienne les "mean,sigma"
//               de l'histo ne sont pas vraiment ce que l'on veut
// --- Limitations actuelle du code:
// . les AGN sont supposes evoluer avec la meme loi de puissance pour tout theta,phi
// . le flux des AGN est mis dans une colonne Oz (indice k) et pas sur la ligne de visee
// . la distribution est approximee a une gaussienne
// ... C'est une approximation pour un observateur loin du centre du cube
//     et pour un cube peu epais / distance observateur 
// --- Parametres de la routine:
// llfy : c'est le <log10(S)> du flux depose par les AGN
//        dans l'angle solide du pixel elementaire de reference du cube
// lsigma : c'est le sigma de la distribution des log10(S)
// powlaw : c'est la pente de la distribution cad que le flux "lmsol"
//          et considere comme le flux a 1.4GHz et qu'on suppose une loi
//             F(nu) = (1.4GHz/nu)^powlaw * F(1.4GHz)
// - Comme on est en echelle log10():
//   on tire log10(Msol) + X
//   ou X est une realisation sur une gaussienne de variance "sig^2"
//   La masse realisee est donc: Msol*10^X
// - Pas de probleme de pixel negatif car on a une multiplication!
{
  if(lp_>0) cout<<"--- AddAGN: <log10(S Jy)> = "<<lfjy<<" , sigma = "<<lsigma<<endl;
  check_array_alloc();

 if(cosmo_ == NULL || redsh_ref_<0.| loscom2zred_.size()<1) {
   char *bla = "GeneFluct3D::AddAGN_Error: set Observator and Cosmology first";
   cout<<bla<<endl; throw ParmError(bla);
 }

 // Le flux des AGN en Jy et en mass solaire
 double fagnref = pow(10.,lfjy)*(dnu_ref_*1.e9); // Jy.Hz = W/m^2
 double magnref = FluxHI2Msol(fagnref*Jansky2Watt_cst,dlum_ref_); // Msol
 if(lp_>0)
   cout<<"Au pixel de ref: fagnref="<<fagnref
       <<" Jy.Hz (a 1.4GHz), magnref="<<magnref<<" Msol"<<endl;

 if(powlaw!=0.) {
   // F(nu) = F(1.4GHz)*(nu GHz/1.4 Ghz)^p = F(1.4GHz)*(1/(1+z))^p , car nu = 1.4 GHz/(1+z)
   magnref *= pow(1/(1.+redsh_ref_),powlaw);
   if(lp_>0) cout<<" powlaw="<<powlaw<<"  -> change magnref to "<<magnref<<" Msol"<<endl;
 }

 // Les infos en fonction de l'indice "l" selon Oz
 vector<double> correction;
 InterpFunc interpinv(loscom2zred_min_,loscom2zred_max_,loscom2zred_);
 long nzmod = ((Nz_>10)?Nz_/10:1);
 for(long l=0;l<Nz_;l++) {
   double z = fabs(DZcom(l));
   double zred = interpinv(z);
   double dtrc = cosmo_->Dtrcom(zred);  // pour variation angle solide
   double dlum = cosmo_->Dlum(zred);  // pour variation conversion mass HI
   double dred = Dz_/(cosmo_->Dhubble()/cosmo_->E(zred));
   double dnu  = Fr_HyperFin_Par *dred/pow(1.+zred,2.); // pour variation dNu
   // on a: Mass ~ DNu * Dlum^2 / Dtrcom^2
   double corr = dnu/dnu_ref_*pow(dtrc_ref_/dtrc*dlum/dlum_ref_,2.);
   // F(nu) = F(1.4GHz)*(nu GHz/1.4 Ghz)^p = F(1.4GHz)*(1/(1+z))^p , car nu = 1.4 GHz/(1+z)
   if(powlaw!=0.) corr *= pow((1.+redsh_ref_)/(1.+zred),powlaw);
   correction.push_back(corr);
   if(lp_>0 && (l==0 || l==Nz_-1 || l%nzmod==0)) {
     cout<<"l="<<l<<" z="<<z<<" red="<<zred<<" dred="<<dred<<" dnu="<<dnu
         <<" dtrc="<<dtrc<<" dlum="<<dlum
         <<" -> cor="<<corr<<endl;
   }
 }

 double sum=0., sum2=0., nsum=0.;
 for(long i=0;i<Nx_;i++) for(long j=0;j<Ny_;j++) {
   double a = lsigma*NorRand();
   a = magnref*pow(10.,a);
   // On met le meme tirage le long de Oz (indice k)
   for(long l=0;l<Nz_;l++) {
     int_8 ip = IndexR(i,j,l);
     data_[ip] += a*correction[l];
   }
   sum += a; sum2 += a*a; nsum += 1.;
 }

 if(lp_>0 && nsum>1.) {
   sum /= nsum;
   sum2 = sum2/nsum - sum*sum;
   cout<<"...Mean mass="<<sum<<" Msol , s^2="<<sum2<<" s="<<sqrt(fabs(sum2))<<endl;
 }
 
}
*********************************************************************/
