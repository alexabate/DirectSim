#include "massfunc.h"

// ::::::::::::  Constants and parameters ::::::::::::
double MassFunc::a_ST_Cst = 0.707;       // Sheth-Torman mass function parameter
double MassFunc::A_ST_Cst = 0.3222;      // Sheth-Torman mass function parameter
double MassFunc::Q_ST_Cst = 0.3;         // Sheth-Torman mass function parameter

// :::::::::::: Constructor ::::::::::::: 
MassFunc::MassFunc(SimpleUniverse& su, PkSpecCalc& pk, double zref, double sig8, 
                                        bool TypeLog, bool IntType, int MFType)
: su_(su) , pk_(pk) , zref_(zref) , TypeLog_(TypeLog) , IntType_(IntType) , MFType_(MFType)
// Because zref_ is an input to this class it shouldn't matter what z pk_ is initially
// calculated at
// MFType_ currently not used: if other mass functions than ST are included then 
// this will be used.
{     

	cout << "     MassFunc::Constructor"<<endl;
	if (TypeLog_)
		cout << "     Outputing mass function dn/dlogM"<<endl;
	else
		cout << "     Outputing mass function dn/dM"<<endl;
	if (IntType_)
		cout << "     Integrating power spectrum using adaptive technique"<<endl;
	else
		cout << "     Integrating power spectrum using quick trapezium technique"<<endl;

	// set lmstep to default value
	lmstep_ = 0.05;

	// cosmo params
	Om_= su_.OmegaMatter();
	double R = 8; // sigma8 definition

	// Normalise to sigma8
	cout << "     Normalise power spectrum to sigma"<< R <<" = "<< sig8 << endl;
	
	// Compute variance at z=0 
	pk_.SetZ(0.);
 	cout << "     Compute variance for top-hat R = "<< R <<" (sigma"<< R;
 	cout << ") at z = "<< pk_.GetZ() <<endl;
 	double sr2 = FindVar(R);
	cout <<"     Current variance = "<< sr2 <<"  ->  sigma"<< R <<" = "<< sqrt(sr2) <<endl;

	// renormalise
	double normpk = sig8*sig8/sr2;
	pk_.SetScale(normpk);
	cout <<"     Spectrum normalisation now = "<< pk_.GetScale() <<endl;
	pk_.SetZ(zref_); // set back to z
	cout << "     ... Finished initialising power spectrum "<<endl;

	cout << "     EXIT MassFunc::Constructor"<<endl<<endl;

}; 


double MassFunc::FindVar(double R) const
{

	VarianceSpectrum varpk(pk_, R, VarianceSpectrum::TOPHAT, IntType_);

	// set k range to integrate spectrum
	double kmin=1e-6,kmax=1000.;
	int npt = 1000;

	double ldlk,k1,k2;
	if (IntType_) {
		// this makes the adaptive integration faster
		double eps=1.e-3;
	 	double kfind_th = varpk.FindMaximum(kmin,kmax,eps);// return k at maximum of P(k)
	 	double pkmax_th = varpk(kfind_th);// returns k^2*P(k)*W^2(kR) -> the integrand of variance integration
	 	k1=kmin, k2=kmax; // not sure need to set these to anything here (these will be integration limits)
	 	varpk.FindLimits(pkmax_th/1.e4,k1,k2,eps); // think this stops the integration doing ranges where P(k) -> 0
		double ldlk = (log10(k2)-log10(k1))/npt; // should be log10 because integration done on log10
		varpk.SetInteg(0.01,ldlk,log10(kmax),4);
		}
	else {
		varpk.SetInteg(npt);
		k1 = kmin; k2 = kmax;
		}

	double sr2 = varpk.Variance(k1,k2);
	return sr2;

};


double MassFunc::fST(double mv) const
// Sheth-Torman function
{
	// the parameters
	double deltc = deltac();
	
	double sigsq = sigsqM(mv); // variance at R = cuberoot[(3M)/(4PI*RHO)]
	double epart = exp( -0.5*aST()*( (deltc*deltc)/sigsq ) );
	double endpart = 1 + 1/( aST()*deltc*deltc/pow(sigsq,QST()) );

	double fst = AST()*sqrt(2*aST()/PI)*(deltc/sqrt(sigsq))*epart*endpart;

	return fst;

};


double MassFunc::sigsqM(double mv) const
// Variance as a function of M
{

	// equivalent radius to mass mv
	double Rcubed = (3*mv)/(4*PI*rhobar0());
	double R = pow(Rcubed,(1./3.));

	// variance of power spectrum at this scale
	double sigsq = FindVar(R);

	return sigsq;

};


double MassFunc::dlnsigdlnm(double mv, double lmstep) const
// Differentiate ln(variance) by ln(mass)
{

    double logm2 = log(mv)+lmstep;
	double m2 = exp(logm2);

	double sigm1 = sqrt(sigsqM(mv));
	double sigm2 = sqrt(sigsqM(m2));
	double diff=( log(sigm2)-log(sigm1) )/( logm2-log(mv) );
	
	return diff;
};


void MassFunc::WritePS2File(string outfile)
{

	cout << endl << "     MassFunc::WritePS2File"<<endl;
	double kmin=1e-6, kmax=1000.;
	int npt = 1000;
	double dk = (kmax-kmin)/(npt-1);

	cout << "     Writing power spectrum to "<< outfile <<endl;
	ifstream inp;
	ofstream outp;

	// check file does not already exist
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) { // file does not exist
	  	
		inp.clear(ios::failbit);
		cout << "     Writing to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
		for (int i=0; i<npt; i++) {
	
			double k = kmin + i*dk;
			double Pk = pk_(k);
			outp << k <<"    "<< Pk <<endl;
			}

	 	outp.close();
	  	}
	else // file *does* exist
		cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

	cout  << "     EXIT MassFunc::WritePS2File"<<endl<< endl;

};

