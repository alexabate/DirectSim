#include "selectfunc.h"

/******* SelectFunc ***********************************************************/

SelectFunc::SelectFunc(double phistar, double Mstar, double alpha)
  : phistar_(phistar) , Mstar_(Mstar) , alpha_(alpha)
{
	// Create Schechter luminosity function
	// See eqn 1.7 in Statistics of the Galaxy Distribution by Martinez and Saar
	//sch_.SetParam(phistar_,Mstar_,alpha_);
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


double SelectFunc::operator() (double LumDist,double mlim, double Mc) const
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
{
	Schechter sch(phistar_,Mstar_,alpha_);

	// Choose really bright min abs mag so LF~0 anyway (essentially -inf)
	double Mmin = -26;
	int npt=50000; // number of points to perform integration with
	
	// Limiting abs mag for observation at LumDist given apparent mag limit is mlim
	// @warning: this assumes bolometric magnitudes!
	double Mlim = mlim - 25 - 5*log10(LumDist);
	
	// Eqn 1.8 is
	// phi(x) = int_-inf^Mlim phi(M)dM / int_-inf^Mc phi(M)dM
	
	// top integral of eqn 1.8
	sch.SetInteg(Mmin, Mlim, npt);
	double no = sch.Integrate();
	
	// bottom integral
	// first check that Mlim is not brighter than Mc 
	// if it is then galaxy DEFINITELY must be observed so top integral should=bottom
	// (will happen if galaxy is at a distance further than most distant galaxy in survey)
	double Mb;
	if (Mlim<Mc) // if Mlim is fainter than Mc
		 Mb = Mc;  // then set upper integral limit to Mc
	else 
		Mb = Mlim; // set to same limit as other integral, then selection func=1
	
	// bottom integral of eqn 1.8
	sch.SetInteg(Mmin, Mb, npt);
	double nt=sch.Integrate();
	//double nt=sch.IntegrateQuick(Mmin,Mb,npt);
	
	// WTF? Why is Mb replaced with Mmax?
	// however i don't get this method so instead
	double Mmax = -13; // some reasonable maximum luminosity
	sch.SetInteg(Mmin,Mmax,npt);
	nt=sch.Integrate();

	double sf=no/nt;
	return sf;
};


double SelectFunc::operator() (double Mlim, double Mc) const
// Computes selection function assuming a luminosity function which 
// has a Schechter function form with parameters: phistar,alpha,Mstar
// Computes eqn 1.8 in Statistics of the Galaxy Distribution by Martinez and Saar
// Mmax = min(Mlim,Mc) - i.e. which is brighter?
{
    Schechter sch(phistar_,Mstar_,alpha_);
 
	// Choose really bright min abs mag so LF~0 anyway (essentially -inf)
	double Mmin=-26;
	int npt=50000; // number of points to perform integration with
	
	// Eqn .18 is
	// phi(x) = int_-inf^Mlim phi(M)dM / int_-inf^Mc phi(M)dM
	
	// top integral:
	sch.SetInteg(Mmin, Mlim, npt);
	double no=sch.Integrate();
	
	// bottom integral
	// first check that Mlim is not brighter than Mc 
	// if it is then galaxy DEFINITELY must be observed so top integral should=bottom
	// (will happen if galaxy is at a distance further than most distant galaxy in survey)
	double Mb;
	if (Mlim<Mc) // if Mlim is fainter than Mc
		 Mb=Mc;  // then set upper integral limit to Mc
	else 
		Mb=Mlim; // set to same limit as other integral, then selection func=1
		
    sch.SetInteg(Mmin,Mb,npt);
	double nt = sch.Integrate();
	//double nt=sch.IntegrateQuick(Mmin,Mb,npt);
	
	// WTF? Why is Mb replaced with Mmax?
	// however i don't get this method so instead:
	double Mmax=-13; // some reasonable maximum luminosity
	sch.SetInteg(Mmin,Mmax,npt);
	nt=sch.Integrate();

	double sf=no/nt;
	return sf;
};


/******* AnalyticalSelFunc ****************************************************/

double AnalyticalSelFunc::SFValue(int type, double redshift, double glong, double glat)
{
    // we can make it faster by interpolating DLum(redshift)  in the constructor
    su_.SetEmissionRedShift(redshift); 
    double dlum = su_.LuminosityDistanceMpc();
    return sf_(dlum, mlim_, Mc_);
};
