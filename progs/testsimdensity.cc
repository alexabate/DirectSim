#include "sopnamsp.h"
#include "machdefs.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "sophyainit.h"
#include "timing.h"
#include "ntuple.h"
#include "matharr.h"
#include "randfmt.h"
//#include "randr48.h"
#include "srandgen.h"

#include "constcosmo.h"
#include "cosmocalcs.h"
#include "schechter.h"
#include "geneutils.h"
#include "genefluct3d.h"

// written by Christophe Magneville
// edited by AA

void usage(void);
void usage(void)
{
 cout<<"testsimdensity [...options...]"<<endl
     <<" -a : auto init random seed (needed for multiple simul)"<<endl
     <<" -0 : use ComputeFourier0 method (defaut: no, use normal way)"<<endl
     <<" -G typevol: compute Pk(z=0) and apply growth factor in real space"<<endl
     <<"       typevol=1 evolved with distance / observer (def)"<<endl
     <<"       typevol=2 evolved with distance to middle of Z planes"<<endl
     <<"       else : no evol, spectrum Pk(z=z_median) for all cube (def)"<<endl
     <<" -F : filter spectrum by pixel shape (0=no 1=yes(default)"<<endl
//     <<" -U typ =0 compute <NGal>, poisson fluctuate then convert to HI mass (def)"<<endl
//     <<"        >0 compute directly <HI mass>, do NOT poisson fluctuate with <Ngal>"<<endl
//     <<"        <0 just multiply dRho/Rho by mass_by_pixel (possible negative pixel"
//     << " values)"<<endl
     <<" -x nx,dx : size along x axis (npix,Mpc)"<<endl
     <<" -y ny,dy : size along y axis (npix,Mpc)"<<endl
     <<"            if ny or dy <=0 take same value as for x"<<endl
     <<" -z nz,dz : size along z axis (redshift axis, npix,Mpc)"<<endl
     <<" -Z zref : redshift for the center of the simulation cube"<<endl
//     <<" -s snoise,typevol : Gaussian noise sigma in equivalent Msol"<<endl
//     <<"           typevol=0 no evolution (def)"<<endl
//     <<"           typevol=1 noise evolved with distance / observer"<<endl
//     <<"           typevol=2 noise evolved with distance to middle of Z planes"<<endl
     <<" -2 : compute also 2D spectrum (default: no)"<<endl
//     <<" -N scalemass: facteur d\'unite pour la masse (default: -1)"<<endl
//     <<"               ex: si on veut unites 10^8 Msol -> scalemass=1.e-8"<<endl
//     <<"               si <0 alors facteur=-scalemass*Mpix"<<endl
     <<" -M schmin,schmax,nsch : min,max mass and nb points for schechter HI"<<endl
     <<"                         If nsch<0 in fact number of points per decade"<<endl
     <<" -Q naleagal : use quick method for turning ngal to mass"<<endl
     <<" -R schmassdist.ppf : read mass distribution for trials from file"<<endl
     <<"               instead of computing it (ONLY if \"-Q\" option is activated)"<<endl
     <<" -8 sigmaR,R : normalisation of power spectrum, R in Mpc"<<endl
     <<"               (default sigmaR=1, R=8/h100 Mpc)"<<endl
     <<" -W : write cube in FITS format (complex cube is coded as real cube)"<<endl
     <<" -P : write cube in PPF format"<<endl
//     <<" -O a,b,c,d,e : tell what you want to write (with the -W and -P options)"<<endl
//     <<"              a=1 : write generated fourier cube (_k0)"<<endl
//     <<"              b=1 : write real space cube dRho/Rho at z (_r0)"<<endl
//     <<"              c=1 : write real space cube conv to gal,poiss,mass etc... (_r)"<<endl
//     <<"              d=1 : write previous real space cube + noise (_rf)"<<endl
//     <<"              e=1 : write fourier transform of previous cube (_k)"<<endl
     <<" -O a,b : tell what you want to write (with the -W and -P options)"<<endl
     <<"              a=1 : write generated fourier cube (_k0)"<<endl
     <<"              b=1 : write real space cube dRho/Rho at z (_r0)"<<endl
     <<" -S : write cube slices in PPF format"<<endl
     <<" -o root_name_out : root string for output file name (def: testsimdensity)"<<endl
//     <<" -V : compute variance from real space (for check, default: no)"<<endl
     <<" -T nth : number of threads (if compiled multi-thread, default: 0)"<<endl
     <<endl;
}

int main(int narg,char *arg[])
{
 SophyaInit();
 InitTim();

 //-----------------------------------------------------------------
 // *** Survey definition
 long nx=360, ny=-1, nz=64; double dx=1., dy=-1., dz=-1.;

 // *** Cosmography definition   (WMAP)
 unsigned short flat = 0;
 double ob0 = 0.0444356;
 double h100=0.71, om0=0.267804, or0=7.9e-05, ol0=0.73,w0=-1.;
 double zref = 0.5;
 double perc=0.01,dzinc=-1.,dzmax=-1.; unsigned short glorder=4;

 // *** Spectrum and variance definition
 double ns = 1., as = 1.;
 double R=8./h100, Rg=R/sqrt(5.);
 double sigmaR = 1.;

 double kmin=1e-5,kmax=1000.;
 int npt = 10000;
 double lkmin=log10(kmin), lkmax=log10(kmax);
 double eps=1.e-3;

 // *** Schechter mass function definition
 double h75 = h100 / 0.75;
 double nstar = 0.006*pow(h75,3.);
 double mstar = pow(10.,9.8);  // MSol
 double alpha = -1.37;

 double schmin=1.e7, schmax=1.e13;
 int schnpt = -100;
 bool use_schmassdist = false;
 long naleagal = 100000;
 bool recompute_schmassdist = true;
 string schmassdistfile = "";
// int no_poisson_type = 0;

// double scalemass = -1.;

// // *** Noise level
// double snoise= 0.;   // in equivalent MSol
// int noise_evol = 0;

 // *** type of generation
 bool computefourier0=false;
 int use_growth_factor = 0;
 unsigned short nthread=0;
 int filter_by_pixel = 1;

 // *** What to do
 bool comp2dspec = false;
 bool wfits = false;
 bool wppf = false;
 bool wslice = false;
// bool compvarreal = false;
 //unsigned short whattowrt[5] = {1,1,1,1,1};
 unsigned short whattowrt[2] = {1,1};
 string rootnameout = "testsimdensity";

 unsigned long ntnent = 10000;  // 0 = do not fill NTuple

 // --- Decode arguments
 if(narg>0) {
   cout<<"\n--- Arguments: "<<endl;
   for(int i=0;i<narg;i++) cout<<arg[i]<<" ";
   cout<<endl;
 }
 system("date -u");

 // --- Choice of random generator (done here before AutoInitRand)
 FMTRandGen *RandGen = new FMTRandGen;
 RandGen->SelectGaussianAlgo(C_Gaussian_RandLibSNorm);
 RandGen->SelectPoissonAlgo(C_Poisson_Ahrens);
 RandomGeneratorInterface::SetGlobalRandGenP(RandGen);
 
 // --- Decode arguments
 char c;
 while((c = getopt(narg,arg,"ha0PWS2:G:F:x:y:z:Z:M:A:T:Q:R:8:O:o:")) != -1) {
  int nth = 0;
  switch (c) {
  case 'a' :
    AutoInitRand(5);
    break;
  case '0' :
    computefourier0 = true;
    break;
  case 'G' :
    sscanf(optarg,"%d",&use_growth_factor);
    break;
//  case 'U' :
//    sscanf(optarg,"%d",&no_poisson_type);
//    break;
  case 'F' :
    sscanf(optarg,"%d",&filter_by_pixel);
    break;
  case 'x' :
    sscanf(optarg,"%ld,%lf",&nx,&dx);
    break;
  case 'y' :
    sscanf(optarg,"%ld,%lf",&ny,&dy);
    break;
  case 'z' :
    sscanf(optarg,"%ld,%lf",&nz,&dz);
    break;
//  case 's' :
//    sscanf(optarg,"%lf,%d",&snoise,&noise_evol);
//    break;
  case 'Z' :
    sscanf(optarg,"%lf",&zref);
    break;
  case '2' :
    comp2dspec = true;
    break;
//  case 'N' :
//    sscanf(optarg,"%lf",&scalemass);
//    if(scalemass==0.) scalemass = 1.;
//    break;
  case 'M' :
    sscanf(optarg,"%lf,%lf,%d",&schmin,&schmax,&schnpt);
    break;
  case 'Q' :
    use_schmassdist = true;
    sscanf(optarg,"%ld",&naleagal);
    break;
  case 'R' :
    schmassdistfile = optarg;
    break;
    ////  case 'A' :
    ////do_agn = true;
    ////sscanf(optarg,"%lf,%lf,%lf",&lfjy_agn,&lsigma_agn,&powlaw_agn);
    ////break;
  case '8' :
    sscanf(optarg,"%lf,%lf",&sigmaR,&R);
    break;
//  case 'V' :
//    compvarreal = true;
//    break;
  case 'W' :
    wfits = true;
    break;
  case 'P' :
    wppf = true;
    break;
  case 'O' :
//    sscanf(optarg,"%hu,%hu,%hu,%hu,%hu"
//          ,&whattowrt[0],&whattowrt[1],&whattowrt[2],&whattowrt[3],&whattowrt[4]);
    sscanf(optarg,"%hu,%hu",&whattowrt[0],&whattowrt[1]);
    break;
  case 'S' :
    wslice = true;
    break;
  case 'o' :
    rootnameout = optarg;
    break;
  case 'T' :
    sscanf(optarg,"%d",&nth);
    nthread = (nth<1)? 0: nth;
    break;
  case 'h' :
  default :
    usage(); return -1;
  }
 }

 double lschmin=log10(schmin), lschmax=log10(schmax);
 if(schnpt<=0) {  // in fact this is the number of pointer per decade
   schnpt = long( (-schnpt)*(lschmax-lschmin+1.) + 0.5 );
   if(schnpt<=2) schnpt = 1000;
 }
 if(naleagal<=2) naleagal = 100000;

 cout<<"zref="<<zref<<endl;
 cout<<"nx="<<nx<<" dx="<<dx<<" ny="<<ny<<" dy="<<dy<<" nz="<<nz<<" dz="<<dz<<endl;
 cout<<"kmin="<<kmin<<" ("<<lkmin<<"), kmax="<<kmax<<" ("<<lkmax<<") Mpc^-1"
     <<", npt="<<npt<<endl;
 cout<<"Filter by pixel = "<<filter_by_pixel<<endl;
 cout<<"R="<<R<<" Rg="<<Rg<<" Mpc, sigmaR="<<sigmaR<<endl;
 cout<<"Use_growth_factor = "<<use_growth_factor<<endl;
 cout<<"nstar= "<<nstar<<"  mstar="<<mstar<<"  alpha="<<alpha<<endl;
 cout<<"schmin="<<schmin<<" ("<<lschmin
     <<"), schmax="<<schmax<<" ("<<lschmax<<") Msol"
     <<", schnpt="<<schnpt<<endl;
// if(no_poisson_type!=0) cout<<"No poisson fluctuation, direct conversion to HI mass, typ="
//                            <<no_poisson_type<<endl;
// cout<<"snoise="<<snoise<<" equivalent Msol, evolution="<<noise_evol<<endl;
// cout<<"scalemass="<<scalemass<<endl;
 cout<<"wfits="<<wfits<<" wppf="<<wppf<<" wslice="<<wslice<<" what?="
//     <<whattowrt[0]<<","<<whattowrt[1]<<","<<whattowrt[2]<<","<<whattowrt[3]<<","
//     <<whattowrt[4]<<endl;
     <<whattowrt[0]<<","<<whattowrt[1]<<endl;
 cout<<"rootnameout="<<rootnameout<<endl;
 ShowRandom();
 cout<<"   First random is: "<<drand01()<<endl;

 string tagobs = rootnameout + ".ppf";
 POutPersist posobs(tagobs);

 //-----------------------------------------------------------------
 cout<<endl<<"\n--- Create Cosmology"<<endl;

 //CosmoCalc univ(flat,true,zref+1.);
 SimpleUniverse univ(h100,om0,ol0);//(flat,true,zref+1.);
 //univ.SetInteg(perc,dzinc,dzmax,glorder);
 //univ.SetDynParam(h100,om0,or0,ol0,w0);
 //univ.PrtInteg();
 //univ.Print();
 univ.SetOmegaBaryon(ob0);// AA added
 univ.SetOmegaPhoton(0);// AA added
 univ.SetOmegaRadiation(0);// AA added
 univ.SetDarkEnergy(ol0,w0);// AA added
 univ.SetEmissionRedShift(zref);// AA added
 //double loscomref = univ.Dloscom(zref);
 double loscomref = univ.LineOfSightComovDistanceMpc();
 cout<<"\nzref = "<<zref<<" -> dloscom = "<<loscomref<<" Mpc"<<endl;
 univ.Print();

 //-----------------------------------------------------------------
 cout<<endl<<"\n--- Create Spectrum"<<endl;

 //InitialSpectrum pkini(ns,as);
 InitialPowerLaw Pkinit(ns,as);

 //TransfertEisenstein tf(h100,om0-ob0,ob0,T_CMB_Par,false);
 TransferEH tf(h100,om0-ob0,ob0,T_CMB_K,false);

 //GrowthFactor growth(om0,ol0);
 GrowthFN growth(om0,ol0);
 double growth_at_z = growth(zref);
 cout<<"...Growth factor at z="<<zref<<" = "<<growth_at_z<<endl;

 //PkSpectrum0 pk0(pkini,tf);
 //PkSpecCalc pk0(Pkinit,tf,growth,0.);

 //PkSpectrumZ pkz(pk0,growth,zref);
 PkSpecCalc pkz(Pkinit,tf,growth,zref);

 //-----------------------------------------------------------------
 pkz.SetZ(0.);
 cout<<endl<<"\n--- Compute variance for top-hat R="<<R
     <<" at z="<<pkz.GetZ()<<endl;
 VarianceSpectrum varpk_th(pkz,R,VarianceSpectrum::TOPHAT);
 double kfind_th = varpk_th.FindMaximum(kmin,kmax,eps);
 double pkmax_th = varpk_th(kfind_th);
 cout<<"kfind_th = "<<kfind_th<<" ("<<log10(kfind_th)<<"), integrand="<<pkmax_th<<endl;
 double k1=kmin, k2=kmax;
 int rc = varpk_th.FindLimits(pkmax_th/1.e4,k1,k2,eps);
 cout<<"limit_th: rc="<<rc<<" : "<<k1<<" ("<<log10(k1)<<") , "
     <<k2<<" ("<<log10(k2)<<")"<<endl;

 double ldlk = (log10(k2)-log10(k1))/npt;
 varpk_th.SetInteg(0.01,ldlk,-1.,4);
 double sr2 = varpk_th.Variance(k1,k2);
 cout<<"varpk_th="<<sr2<<"  ->  sigma="<<sqrt(sr2)<<endl;

 // normalise the power spectrum to the values given for sigmaR,R
 double normpkz = sigmaR*sigmaR/sr2;
 pkz.SetScale(normpkz);
 cout<<"Spectrum normalisation = "<<pkz.GetScale()<<endl;

 {
 Histo hpkz0(lkmin,lkmax,npt); hpkz0.ReCenterBin();
 FuncToHisto(pkz,hpkz0,true);
 tagobs = "hpkz0"; posobs.PutObject(hpkz0,tagobs);
 }

 // set power spectrum back to be at zref
 pkz.SetZ(zref);

 {
 Histo hpkz(lkmin,lkmax,npt); hpkz.ReCenterBin();
 FuncToHisto(pkz,hpkz,true);
 tagobs = "hpkz"; posobs.PutObject(hpkz,tagobs);
 }

 //-----------------------------------------------------------------
 cout<<endl<<"\n--- Compute variance for Pk at z="<<pkz.GetZ()<<endl;
 VarianceSpectrum varpk_int(pkz,R,VarianceSpectrum::NOFILTER);

 double kfind_int = varpk_int.FindMaximum(kmin,kmax,eps);
 double pkmax_int = varpk_int(kfind_int);
 cout<<"kfind_int = "<<kfind_int<<" ("<<log10(kfind_int)<<"), integrand="<<pkmax_int<<endl;
 double k1int=kmin, k2int=kmax;
 int rcint = varpk_int.FindLimits(pkmax_int/1.e4,k1int,k2int,eps);
 cout<<"limit_int: rc="<<rcint<<" : "<<k1int<<" ("<<log10(k1int)<<") , "
     <<k2int<<" ("<<log10(k2int)<<")"<<endl;

 double ldlkint = (log10(k2int)-log10(k1int))/npt;
 varpk_int.SetInteg(0.01,ldlkint,-1.,4);
 double sr2int = varpk_int.Variance(k1int,k2int);
 cout<<"varpk_int="<<sr2int<<"  ->  sigma="<<sqrt(sr2int)<<endl;
 
 //-----------------------------------------------------------------
 cout<<endl<<"\n--- Create mass function, compute number/mass density, init mass trials";
 cout<<endl;

 Schechter sch(nstar,mstar,alpha);
 sch.Print();
 sch.SetInteg(schmin,schmax,schnpt); // AA added

 sch.SetOutValue(0);
 cout<<"sch(mstar) = "<<sch(mstar)<<" /Mpc^3/Msol"<<endl;
 //double ngal_by_mpc3 = sch.Integrate(schmin,schmax,schnpt);
 double ngal_by_mpc3 = sch.Integrate();
 cout<<"Galaxy density number = "<<ngal_by_mpc3<<" /Mpc^3 between limits"<<endl;

 sch.SetOutValue(1);
 cout<<"mstar*sch(mstar) = "<<sch(mstar)<<" Msol/Mpc^3/Msol"<<endl;
 //double mass_by_mpc3 = sch.Integrate(schmin,schmax,schnpt);
 double mass_by_mpc3 = sch.Integrate();
 cout<<"Galaxy mass density= "<<mass_by_mpc3<<" Msol/Mpc^3 between limits"<<endl;
 //cout<<"Omega_HI at z=0 is "<<mass_by_mpc3/(univ.Rhoc(0.)*GCm3toMsolMpc3_Cst)<<endl
 //    <<"         at z="<<zref<<" is "<<mass_by_mpc3/(univ.Rhoc(zref)*GCm3toMsolMpc3_Cst)<<endl;
 cout<<"Omega_HI at z=0 is "<<mass_by_mpc3/(univ.CriticalDensityMSolMpc3())<<endl
    <<"         at z="<<zref<<" is "<<mass_by_mpc3/(univ.CriticalDensityMSolMpc3())<<endl;


 SchechterMassDist schmdist(sch,schmin,schmax,schnpt);
 if(use_schmassdist && schmassdistfile.size()>0) {
   cout<<"\nSchechterMassDist read from "<<schmassdistfile<<endl;
   schmdist.ReadPPF(schmassdistfile);
   bool compsch = IsCompatible(sch,schmdist.GetSchechter());
   double m1,m2; schmdist.GetMassLim(m1,m2);
   if( !compsch || fabs(m1-schmin)>1e-4*schmin || fabs(m2-schmax)>1e-4*schmax) {
     cout<<"FATAL_ERROR: INCONSISTENT SchechterMassDist file / Schechter or limits"<<endl;
     cout<<"Schechter:"<<endl; sch.Print();
     cout<<"schmin="<<schmin<<" schmax="<<schmax<<endl;
     cout<<"SchechterMassDist:"<<endl; schmdist.Print();
     return -10;
   }
   recompute_schmassdist = false;
 }
 schmdist.Print();
 Histo hmdndm = schmdist.GetHmDnDm();
 FunRan tirhmdndm = schmdist.GetTmDnDm();
 {
 tagobs = "hmdndm"; posobs.PutObject(hmdndm,tagobs);
 Histo hdum1(tirhmdndm);
 tagobs = "tirhmdndm"; posobs.PutObject(hdum1,tagobs);
 }

 PrtTim(">>>> End of definition");

 //-----------------------------------------------------------------
 // FFTW3 (p26): faster if sizes 2^a 3^b 5^c 7^d 11^e 13^f  with e+f=0 ou 1
 cout<<endl<<"\n--- Initialisation of GeneFluct3D"<<endl;

 GeneFluct3D fluct3d(nx,ny,nz,dx,dy,dz,nthread,2);
 fluct3d.SetObservator(zref,-nz/2.);
 fluct3d.SetCosmology(univ);
 fluct3d.SetGrowthFactor(growth);
 fluct3d.LosComRedshift(0.001,-1);
 //TArray< complex<GEN3D_TYPE> >& pkgen = fluct3d.GetComplexArray();
 //TArray<GEN3D_TYPE>& rgen = fluct3d.GetRealArray();
 cout<<endl; fluct3d.Print();
 cout<<"\nMean number of galaxies per pixel = "<<ngal_by_mpc3*fluct3d.GetDVol()<<endl;
 double mass_by_pixel = mass_by_mpc3 * fluct3d.GetDVol();
 cout<<"Mean mass per pixel = "<<mass_by_pixel<<endl;

 double dkmin = fluct3d.GetKincMin();
 double knyqmax = fluct3d.GetKmax();
 long nherr = long(knyqmax/dkmin+0.5);
 cout<<"\nFor HistoErr: d="<<dkmin<<" max="<<knyqmax<<" n="<<nherr<<endl;

 double dktmin = fluct3d.GetKTincMin();
 double ktnyqmax = fluct3d.GetKTmax();
 long nherrt = long(ktnyqmax/dktmin+0.5);
 double dkzmin = fluct3d.GetKinc()[2];
 double kznyqmax = fluct3d.GetKnyq()[2];
 long nherrz = long(kznyqmax/dkzmin+0.5);
 cout<<"For Histo2DErr: d="<<dktmin<<","<<dkzmin
     <<" max="<<ktnyqmax<<","<<kznyqmax<<" n="<<nherrt<<","<<nherrz<<endl;

 //-----------------------------------------------------------------
 cout<<"\n--- Computing spectra variance up to Kmax at z="<<pkz.GetZ()<<endl;
 // In fact, works on a cube inside a sphere of radius kmax:
 // sphere: Vs = 4Pi/3 k^3 , cube inside (side = k*sqrt(2)): Vc = (k*sqrt(2))^3
 // Vc/Vs = 0.675   ->  keff = kmax * (0.675)^(1/3) = kmax * 0.877
 double knyqmax_mod = 0.877*knyqmax;
 ldlkint = (log10(knyqmax_mod)-log10(k1int))/npt;
 varpk_int.SetInteg(0.01,ldlkint,-1.,4);
 double sr2int_kmax = varpk_int.Variance(k1int,knyqmax_mod);
 cout<<"varpk_int(<"<<knyqmax_mod<<")="<<sr2int_kmax<<"  ->  sigma="<<sqrt(sr2int_kmax)<<endl;

 PrtTim(">>>> End Initialisation de GeneFluct3D");

 //-----------------------------------------------------------------
 cout<<"\n--- Computing a realization in Fourier space"<<endl;
 if(use_growth_factor>0) pkz.SetZ(0.); else pkz.SetZ(zref);
 cout<<"Power spectrum set at redshift: "<<pkz.GetZ()<<endl;
 if(computefourier0) fluct3d.ComputeFourier0(pkz);
   else              fluct3d.ComputeFourier(pkz);
 fluct3d.NTupleCheck(posobs,string("ntpkgen"),ntnent);
 PrtTim(">>>> End Computing a realization in Fourier space");

 cout<<"\n--- Checking realization spectra"<<endl;
 HistoErr hpkgen(0.,knyqmax,nherr);
 hpkgen.ReCenterBin(); hpkgen.Zero();
 hpkgen.Show();
 fluct3d.ComputeSpectrum(hpkgen);
 {
 tagobs = "hpkgen"; posobs.PutObject(hpkgen,tagobs);
 }
 PrtTim(">>>> End Checking realization spectra");

 if(comp2dspec) {
   cout<<"\n--- Checking realization 2D spectra"<<endl;
   Histo2DErr hpkgen2(0.,ktnyqmax,nherrt,0.,kznyqmax,nherrz);
   hpkgen2.ReCenterBin(); hpkgen2.Zero();
   hpkgen2.Show();
   fluct3d.ComputeSpectrum2D(hpkgen2);
   {
   tagobs = "hpkgen2"; posobs.PutObject(hpkgen2,tagobs);
   }
   PrtTim(">>>> End Checking realization 2D spectra");
 }

 if(filter_by_pixel!=0) {
   cout<<"\n--- Computing convolution by pixel shape"<<endl;
   fluct3d.FilterByPixel();
   fluct3d.NTupleCheck(posobs,string("ntpkgenf"),ntnent);
   PrtTim(">>>> End Computing convolution by pixel shape");

   cout<<"\n--- Checking realization spectra after pixel shape convol."<<endl;
   HistoErr hpkgenfb(0.,knyqmax,nherr);
   hpkgenfb.ReCenterBin(); hpkgenfb.Zero();
   hpkgenfb.Show();
   fluct3d.ComputeSpectrum(hpkgenfb);
   {
   tagobs = "hpkgenfb"; posobs.PutObject(hpkgenfb,tagobs);
   }
   PrtTim(">>>> End Checking realization spectra");

   cout<<"\n--- Checking realization spectra after pixel shape convol. with pixel correc."<<endl;
   HistoErr hpkgenf(hpkgenfb); hpkgenf.Zero();
   fluct3d.ComputeSpectrum(hpkgenf,0.,filter_by_pixel);
   {
   tagobs = "hpkgenf"; posobs.PutObject(hpkgenf,tagobs);
   }
   PrtTim(">>>> End Checking realization spectra with pixel correc.");

   if(comp2dspec) {
     cout<<"\n--- Checking realization 2D spectra after pixel shape convol."<<endl;
     Histo2DErr hpkgenfb2(0.,ktnyqmax,nherrt,0.,kznyqmax,nherrz);
     hpkgenfb2.ReCenterBin(); hpkgenfb2.Zero();
     hpkgenfb2.Show();
     fluct3d.ComputeSpectrum2D(hpkgenfb2);
     {
     tagobs = "hpkgenfb2"; posobs.PutObject(hpkgenfb2,tagobs);
     }
     PrtTim(">>>> End Checking realization 2D spectra");

     cout<<"\n--- Checking realization 2D spectra after pixel shape convol. with pixel correc."<<endl;
     Histo2DErr hpkgenf2(hpkgenfb2); hpkgenf2.Zero();
     fluct3d.ComputeSpectrum2D(hpkgenf2,0.,filter_by_pixel);
     {
     tagobs = "hpkgenf2"; posobs.PutObject(hpkgenf2,tagobs);
     }
     PrtTim(">>>> End Checking realization 2D spectra with pixel correc.");
   }
 }

 if(whattowrt[0]==1) {
 if(wfits) {
   tagobs = "!" + rootnameout + "_k0.fits";
   fluct3d.WriteFits(tagobs);
   PrtTim(">>>> End WriteFits");
 }
 if(wppf) {
   tagobs = rootnameout + "_k0.ppf";
   fluct3d.WritePPF(tagobs,false);
   PrtTim(">>>> End WritePPF");
 }
 }

 //-----------------------------------------------------------------
 cout<<"\n--- Computing a realization in real space"<<endl;
 fluct3d.ComputeReal();
 double rmin,rmax; fluct3d.MinMax(rmin,rmax);
 cout<<"rgen.Min = "<<rmin<<" , Max="<<rmax<<endl;
 fluct3d.NTupleCheck(posobs,string("ntreal"),ntnent);
 PrtTim(">>>> End Computing a realization in real space");

 if(use_growth_factor>0) {
   cout<<"\n--- Apply Growth factor"<<endl;
   cout<<"...D(z=0)="<<growth(0.)<<"  D(z="<<zref<<")="<<growth(zref)<<endl;
   fluct3d.ApplyGrowthFactor(use_growth_factor);
   fluct3d.MinMax(rmin,rmax);
   cout<<"rgen.Min = "<<rmin<<" , Max="<<rmax<<endl;
   fluct3d.NTupleCheck(posobs,string("ntgrow"),ntnent);
   PrtTim(">>>> End Applying growth factor");
 }

 int_8 nm;
 double rmref,rs2ref;
 cout<<"\n--- Computing reference variance in real space"<<endl;
 nm = fluct3d.MeanSigma2(rmref,rs2ref);
 cout<<" rs2ref= "<<rs2ref<<" , rmref="<<rmref<<" ("<<nm<<")"<<endl;
 PrtTim(">>>> End Computing reference variance in real space");

 if(whattowrt[1]==1) {
 if(wfits) {
   tagobs = "!" + rootnameout + "_r0.fits";
   fluct3d.WriteFits(tagobs);
   PrtTim(">>>> End WriteFits");
 }
 if(wppf) {
   tagobs = rootnameout + "_r0.ppf";
   fluct3d.WritePPF(tagobs,true);
   PrtTim(">>>> End WritePPF");
 }
 if(wslice) {
   tagobs = rootnameout + "_s_r0.ppf";
   fluct3d.WriteSlicePPF(tagobs);
   PrtTim(">>>> End WriteSlicePPF");
 }
 }

 //-----------------------------------------------------------------
 delete RandGen;
 PrtTim(">>>> End Of Job");

 return 0;
}

