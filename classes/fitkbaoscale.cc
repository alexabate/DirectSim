#include "fitkbaoscale.h"

//FitBAOScale::FitBAOScale(TVector<r_8> kobs, TVector<r_8> Pobs, TVector<r_8> Pref, double h)
//: kobs_(kobs) , Pobs_(Pobs) , Pref_(Pref) , h_(h)
//{
//
//SetRatioObsSpec();
//
//};

FitBAOScale::FitBAOScale(TArray<r_8> pspectrum, SimpleUniverse& su, double zref,double sig8,double n)
: su_(su) , zref_(zref)
{

	
	// fill power spectrum vectors
	InitVect(pspectrum);

	// Smooth power spectrum
 	double OmegaM=su_.OmegaMatter();
	double OmegaB=su_.OmegaBaryon();
	double OmegaL=su_.OmegaLambda();
	double R=8; // sigma8 definition
	h_=su_.h();
	cout << "     h = "<<h_<<endl;
	ComputeSmoothPS(OmegaM, OmegaL, OmegaB, h_, sig8, n, R);

	// Compute Pobs/Pref
	Pratio();
	
	// set default ka range
	minka_=0.03; maxka_=0.05; nka_=1000;

	bestfit_=-10; // uninitialised
}

void FitBAOScale::InitVect(TArray<r_8> pspectrum)
{
	int nk = pspectrum.SizeY(); // 2nd dim is k value direction
	cout << "     Power spectrum has "<<nk<<" k values"<<endl;
	
	// set sizes of vectors
	kobs_.SetSize(nk); Pobs_.SetSize(nk); sig_.SetSize(nk);
	Pref_.SetSize(nk);
	Pratio_.SetSize(nk);

	for (int i=0; i<nk; i++)
		{
		kobs_(i) = pspectrum(0,i);
		Pobs_(i) = pspectrum(1,i);
		sig_(i) = pspectrum(2,i);
		//if (i<100)
		//	cout <<	kobs_(i) << "    "<<Pobs_(i) << "    "<<sig_(i) <<endl;
		}
	

}

// Calculate SMOOTH fiducial power spectrum
void FitBAOScale::ComputeSmoothPS(double OmegaM, double OmegaL, double OmegaB, double h, double sig8, double n, double R)
{
	
	cout << "     Computing smooth fiducial power spectrum ..."<<endl;
	InitialSpectrum Pkinit(n);
	TransferEH tf(h,OmegaM-OmegaB,OmegaB,T_CMB_Par,false);
	tf.SetNoOscEnv(1); // want smooth power spectrum
	GrowthFactor growth(OmegaM,OmegaL);
 	double growth_at_z = growth(zref_);
 	cout<< "     Growth factor at z="<<zref_<<" = "<<growth_at_z<<endl;
 	PkSpectrum0 pk0(Pkinit,tf);
 	PkSpectrumZ pkz(pk0,growth,zref_);
	pkz.SetZ(0.);
 	cout<<endl<<"     Compute variance for top-hat R="<<R<<" (sigma"<<R<<") at z="<<pkz.GetZ()<<endl;
 	VarianceSpectrum varpk_th(pkz,R,VarianceSpectrum::TOPHAT);
	double kmin=1e-5,kmax=1000.;
	int npt = 10000;
	double lkmin=log10(kmin), lkmax=log10(kmax);
	double eps=1.e-3;
 	double kfind_th = varpk_th.FindMaximum(kmin,kmax,eps);
 	double pkmax_th = varpk_th(kfind_th);
 	cout<<"     kfind_th = "<<kfind_th<<" ("<<log10(kfind_th)<<"), integrand="<<pkmax_th<<endl;
 	double k1=kmin, k2=kmax;
 	int rc = varpk_th.FindLimits(pkmax_th/1.e4,k1,k2,eps);
 	cout<<"     limit_th: rc="<<rc<<" : "<<k1<<" ("<<log10(k1)<<") , "<<k2<<" ("<<log10(k2)<<")"<<endl;
	double ldlk = (log10(k2)-log10(k1))/npt;
	varpk_th.SetInteg(0.01,ldlk,-1.,4);
	double sr2 = varpk_th.Variance(k1,k2);
	cout<<"     varpk_th="<<sr2<<"  ->  sigma="<<sqrt(sr2)<<endl;
	double normpkz = sig8*sig8/sr2;
	pkz.SetScale(normpkz);
	cout<<"     Spectrum normalisation = "<<pkz.GetScale()<<endl;
	pkz.SetZ(zref_);
	 cout<<endl<<"     Compute variance for Pk at z="<<pkz.GetZ()<<endl;
 	VarianceSpectrum varpk_int(pkz,R,VarianceSpectrum::NOFILTER);
	double kfind_int = varpk_int.FindMaximum(kmin,kmax,eps);
 	double pkmax_int = varpk_int(kfind_int);
 	cout<<"     kfind_int = "<<kfind_int<<" ("<<log10(kfind_int)<<"), integrand="<<pkmax_int<<endl;
 	double k1int=kmin, k2int=kmax;
 	int rcint = varpk_int.FindLimits(pkmax_int/1.e4,k1int,k2int,eps);
 	cout<<"     limit_int: rc="<<rcint<<" : "<<k1int<<" ("<<log10(k1int)<<") , "<<k2int<<" ("<<log10(k2int)<<")"<<endl;
	double ldlkint = (log10(k2int)-log10(k1int))/npt;
	varpk_int.SetInteg(0.01,ldlkint,-1.,4);
	double sr2int = varpk_int.Variance(k1int,k2int);
	cout<<"     varpk_int="<<sr2int<<"  ->  sigma="<<sqrt(sr2int)<<endl;
	cout << "     ... Finished computing smooth fiducial power spectrum "<<endl;
	cout << endl;

	for (int i=0; i<kobs_.Size(); i++)
		Pref_(i) = pkz(kobs_(i));

}

void FitBAOScale::Pratio()
{
	for (int i=0; i<kobs_.Size(); i++)
		{Pratio_(i) = Pobs_(i)/Pref_(i);
		 //cout << Pratio_(i) <<"   "<<Pobs_(i)<<"   "<<Pref_(i)<<endl;
		}

}
void FitBAOScale::CalcSigSampVar(double VolCat)
// Analytical approximation of the error on the power spectrum
// due to the sample variance
// Won't be accurate for photo-z power spectrum I expect
{

sig_.SetSize(kobs_.Size());
double dk = kobs_(2)-kobs_(1);
for (int i=0; i<kobs_.Size(); i++)
	{
	
	sig_(i)=sqrt(2*pow(2*PI,3)/VolCat*(1./(4*PI*pow(kobs_(i),2)*dk)))*Pratio_(i);
	
	}
};

void FitBAOScale::ComputeChisq(double maxk)
{

	double dka = (maxka_-minka_)/(nka_-1);
	double amp=2.5;// amplitude is fixed
	kavals_.SetSize(nka_);
	Chisq_.SetSize(nka_);
	cout << "    Computing chisq from ka = "<<minka_<<" to ka = "<<maxka_<<" in steps of "<<dka<<endl;
	

	DecaySineFunc decs(amp, h_);
	cout << "    Check decaying sine function"<<endl;
	//double ka =0.04;
	//for (int ik=0; ik<kobs_.Size(); ik++)
	//		{
	//		double pred=decs(kobs_(ik),ka);
	//		cout << kobs_(ik) <<"   "<< pred<<endl;
	//		}
	decs.PrintParas();

	for (int ika=0; ika<nka_; ika++)
		{
		 double chisq=0;
		 kavals_(ika) = minka_ + ika*dka;
	 
		for (int ik=0; ik<kobs_.Size(); ik++)
			{
			double pred=decs(kobs_(ik),kavals_(ika));
			double diff= Pratio_(ik)-pred;
			
			if (sig_(ik)==0)
				throw ParmError("Error is zero!");


			if (kobs_(ik)<maxk)
				chisq += diff*diff/(sig_(ik)*sig_(ik));
			//if (isnan(pow(diff,2.)/(sig_(ik)*sig_(ik))))
				//cout << "sig_(ik)="<<sig_(ik)<<", Pref^o="<<Pratio_(ika)<<", Pref^p="<<pred<<endl;
			}
	
		Chisq_(ika)=chisq;
		}

}

void FitBAOScale::WriteChisq(string outfile)
{

	if (bestfit_<0)
		throw ParmError("Have not calculated best fit ka yet!");

	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
		{
		inp.clear(ios::failbit);
		cout << "    Writing chisq to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		outp << "Redshift of power spectrum = "<<zref_<<", best fit ka = "<<bestfit_<<", "<<nsig_<<"-sigma errors: + "<<errup_<<" - "<<errdown_<<endl;
		for (int i=0;i<Chisq_.Size();i++)// loop over ka values
			outp << kavals_(i) << "   "<< Chisq_(i)<< endl;
		outp.close();
		} // end of write
		else
			cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

}

void FitBAOScale::WriteRefPS(string outfile)
{

	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
		{
		inp.clear(ios::failbit);
		cout << "    Writing reference power spectrum to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
		for (int i=0;i<Pref_.Size();i++)// loop over k values
			outp << kobs_(i) << "   "<< Pref_(i)<< endl;
		outp.close();
		} // end of write
		else
			cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

}

void FitBAOScale::WriteAncillaryInfo(string outfile)
{

	if (bestfit_<0)
		throw ParmError("Have not calculated best fit ka yet!");

	double amp=2.5;
	DecaySineFunc decs(amp,h_);

	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
		{
		inp.clear(ios::failbit);
		cout << "    Writing reference power spectrum and decaying sinosoid to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
		for (int i=0;i<Pref_.Size();i++)// loop over k values
			{
			double pred=decs(kobs_(i),bestfit_);
			outp << kobs_(i) << "   "<< Pref_(i)<< "   "<< pred <<endl;
			}
		outp.close();
		} // end of write
		else
			cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

}


void FitBAOScale::BestfitStdDev(double& bestfit, double& siglow, double& sighigh, int nsig)
{
	nsig_ = nsig;
	if (nsig_>3)
		throw ParmError("Can't compute above 3-sigma");

	vector<double> clevels;
	clevels.push_back(0.683); clevels.push_back(0.954); clevels.push_back(0.997);

	double clevel = clevels[nsig_-1];
	cout <<"    Confidence level = "<<clevel<<endl;
	ChisqStats chisqstat(kavals_,Chisq_);
	bestfit = chisqstat.BestFit();
	chisqstat.ErrSig(siglow,sighigh,clevel,100);

	bestfit_=bestfit;
	errup_  = sighigh - bestfit_;
	errdown_= bestfit_ - siglow;
}

void FitBAOScale::WriteResults(string outfile)
{

	ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail())
		{
		inp.clear(ios::failbit);
		cout << "    Writing results to file ..." << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
		cout << "    Redshift of power spectrum = "<<zref_<<", best fit ka = "<<bestfit_<<", "<<nsig_<<"-sigma errors: + "<<errup_<<" - "<<errdown_<<endl;
		outp << "z_ps : best fit ka : n-sigma of errors : +error : -error "<<endl;
		outp <<zref_<<"   "<<bestfit_<<"   "<<nsig_<<"   "<<errup_<<"   "<<errdown_<<endl;
		outp.close();
		} // end of write
		else
			cout << "Error...file """ << outfile.c_str() << """ exists" << endl;

}
