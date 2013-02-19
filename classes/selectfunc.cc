#include "selectfunc.h"
#include "schechter.h"
#include "geneutils.h"

/* --Methode-- */
ComputedSelFunc::ComputedSelFunc(string filename)
{
  myinterp_.ReadXYFromFile(filename);
}

/* --Methode-- */
double ComputedSelFunc::SFValue(int type, double redshift, double glong, double glat)
{
  return myinterp_.YInterp(redshift);
}

/* --Methode-- */
AnalyticalSelFunc::AnalyticalSelFunc(SimpleUniverse& su, SelectFunc& sf, double mlim, double Mc)
: su_(su), sf_(sf), mlim_(mlim), Mc_(Mc)
{
}

/* --Methode-- */
double AnalyticalSelFunc::SFValue(int type, double redshift, double glong, double glat)
{
   // we can make it faster by interpolating DLum(redshift)  in the constructor
  su_.SetEmissionRedShift(redshift); 
  double dlum = su_.LuminosityDistanceMpc();
  return sf_(dlum, mlim_, Mc_);
}

///////////////////////////////////////////////////////////
//*********** Compute selection Function from LF *********//
///////////////////////////////////////////////////////////

SelectFunc::SelectFunc(double phistar,double Mstar,double alpha)
  : phistar_(phistar) , Mstar_(Mstar) , alpha_(alpha)
{
	
	// Create Schechter luminosity function
	// See eqn 1.7 in Statistics of the Galaxy Distribution by Martinez and Saar
	sch_.SetParam(phistar_,Mstar_,alpha_);
}

//SelectFunc::SelectFunc(Schechter& f, SimpleUniverse& su)
//  : phistar_(f.phistar_) , Mstar_(f.Mstar_) , alpha_(f.alpha_)
//{
//}
//
//SelectFunc::SelectFunc(void)
//  : phistar_(0.) , Mstar_(0.) , alpha_(0.)
//{
//}
//
//SelectFunc::SelectFunc(void)
//{
//}


void SelectFunc::SetParam(double phistar,double Mstar,double alpha)
{
  phistar_ = phistar; Mstar_ = Mstar; alpha_ = alpha;
}

void SelectFunc::GetParam(double& phistar,double& Mstar,double& alpha)
{
  phistar = phistar_; Mstar = Mstar_; alpha = alpha_;
}

double SelectFunc::operator() (double LumDist,double mlim, double Mc)
// Computes selection function assuming a luminosity function which 
// has a Schechter function form with parameters: phistar,alpha,Mstar
// Computes eqn 1.8 in Statistics of the Galaxy Distribution by Martinez and Saar
// LumDist is the luminosity distance
// Mmax = min(Mlim,Mc) - i.e. which is brighter?
//	where Mc is the absolute magnitude for which the catalog is COMPLETE
//		- i.e. the faintest luminosity a galaxy can have to be included at any distance
//		  Mlim = mlim - 25 - 5log10(LumDist) - calculated in this function
//		  Therefore Mc = mlim - 25 - 5log10(max(LumDist))
// mlim is the magnitude limit of the survey
// not taking into account 
{
	// Choose really bright min abs mag so LF~0 anyway (essentially -inf)
	double Mmin=-26;
	int npt=50000; // number of points to perform integration with
	
	// Limiting abs mag for observation at LumDist given apparent mag limit is mlim
	double Mlim=mlim-25-5*log10(LumDist);
	
	// Eqn .18 is
	// phi(x) = int_-inf^Mlim phi(M)dM / int_-inf^Mc phi(M)dM
	
	// top integral:
	double no=sch_.Integrate(Mmin,Mlim,npt);
	
	// bottom integral
	// first check that Mlim is not brighter than Mc 
	// if it is then galaxy DEFINITELY must be observed so top integral should=bottom
	// (will happen if galaxy is at a distance further than most distant galaxy in survey)
	double Mb;
	if (Mlim<Mc) // if Mlim is fainter than Mc
		 Mb=Mc;  // then set upper integral limit to Mc
	else 
		Mb=Mlim; // set to same limit as other integral, then selection func=1
	double nt=sch_.Integrate(Mmin,Mb,npt);
	//double nt=sch_.IntegrateQuick(Mmin,Mb,npt);
	
	// however i don't get this method so instead:
	double Mmax=-13; // some reasonable maximum luminosity
	nt=sch_.Integrate(Mmin,Mmax,npt);

	double sf=no/nt;
	return sf;
}

double SelectFunc::operator() (double Mlim, double Mc)
// Computes selection function assuming a luminosity function which 
// has a Schechter function form with parameters: phistar,alpha,Mstar
// Computes eqn 1.8 in Statistics of the Galaxy Distribution by Martinez and Saar
// LumDist is the luminosity distance
// Mmax = min(Mlim,Mc) - i.e. which is brighter?
//	where Mc is the absolute magnitude for which the catalog is COMPLETE
//		- i.e. the faintest luminosity a galaxy can have to be included at any distance
//		  Mlim = mlim - 25 - 5log10(LumDist) - calculated in this function
//		  Therefore Mc = mlim - 25 - 5log10(max(LumDist))
// mlim is the magnitude limit of the survey
// not taking into account 
{
	// Choose really bright min abs mag so LF~0 anyway (essentially -inf)
	double Mmin=-26;
	int npt=50000; // number of points to perform integration with
	
	// Eqn .18 is
	// phi(x) = int_-inf^Mlim phi(M)dM / int_-inf^Mc phi(M)dM
	
	// top integral:
	double no=sch_.Integrate(Mmin,Mlim,npt);
	
	// bottom integral
	// first check that Mlim is not brighter than Mc 
	// if it is then galaxy DEFINITELY must be observed so top integral should=bottom
	// (will happen if galaxy is at a distance further than most distant galaxy in survey)
	double Mb;
	if (Mlim<Mc) // if Mlim is fainter than Mc
		 Mb=Mc;  // then set upper integral limit to Mc
	else 
		Mb=Mlim; // set to same limit as other integral, then selection func=1
	double nt=sch_.Integrate(Mmin,Mb,npt);
	//double nt=sch_.IntegrateQuick(Mmin,Mb,npt);
	
	// however i don't get this method so instead:
	double Mmax=-13; // some reasonable maximum luminosity
	nt=sch_.Integrate(Mmin,Mmax,npt);

	double sf=no/nt;
	return sf;
}


void SelectFunc::Print()
{
  cout<<"    SelectFunc::Print: phistar="<<phistar_<<", Mstar="<<Mstar_<<", alpha="<<alpha_;
  cout<<endl;
}
