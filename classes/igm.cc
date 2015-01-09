#include "igm.h"
//#include <string>
//#include <cmath>

/*  */

/******* AtomicCalcs **********************************************************/

AtomicCalcs::AtomicCalcs()
{
    setGammas();
    setOscillatorStrength();
    setConsants();
};

/******* AtomicCalcs methods **************************************************/

void AtomicCalcs::setGammas()
{
        // Found in Table 2 of Morton 2003
    gammaSeries_.push_back(6.265e8);// 2
    gammaSeries_.push_back(1.897e8);// 3
    gammaSeries_.push_back(8.127e7);// 4
    gammaSeries_.push_back(4.204e7);// 5
    gammaSeries_.push_back(2.450e7);// 6
    
    // Backed out of a file containg values of the damping parameter a assuming 
    // b=36km/s from Tepper-Garcia.  See Fig 1 in Tepper-Garcia 2006
    gammaSeries_.push_back(1.236e7);// 7
    gammaSeries_.push_back(8.252e6);// 8
    gammaSeries_.push_back(5.781e6);// 9
    gammaSeries_.push_back(4.209e6);// 10
    gammaSeries_.push_back(3.159e6);// 11
    
    gammaSeries_.push_back(2.431e6);// 12
    gammaSeries_.push_back(1.91e6); 
    gammaSeries_.push_back(4.522e5);// 21 
    
    gammaSeries_.push_back(3.932e5);// 22
    gammaSeries_.push_back(3.442e5);// 23
    gammaSeries_.push_back(3.029e5);// 24
    gammaSeries_.push_back(2.678e5);// 25
    gammaSeries_.push_back(2.381e5);// 26
    
    gammaSeries_.push_back(2.127e5);// 27
    gammaSeries_.push_back(1.906e5);// 28
    gammaSeries_.push_back(1.716e5);// 29
    gammaSeries_.push_back(1.549e5);// 30
    gammaSeries_.push_back(1.405e5);// 31

    gammaSeries_.push_back(1.277e5);// 32
    // size of gammaSeries =  31

    return;
};

void AtomicCalcs::setOscillatorStrength()
{
        // Absorption oscillator strength of ith Lyman line (unitless)
    // Taken from Wiese et al 1966
    fSeries_.push_back(0.4162); 
    fSeries_.push_back(7.910e-2);  
    fSeries_.push_back(2.899e-2);  
    fSeries_.push_back(1.394e-2);
    fSeries_.push_back(7.799e-3);
    
    fSeries_.push_back(4.814e-3);
    fSeries_.push_back(3.183e-3);
    fSeries_.push_back(2.216e-3);
    fSeries_.push_back(1.605e-3);
    fSeries_.push_back(1.201e-3);
    
    fSeries_.push_back(9.214e-4);
    fSeries_.push_back(7.227e-4);
    fSeries_.push_back(5.774e-4);
    fSeries_.push_back(4.686e-4);
    fSeries_.push_back(3.856e-4);
    
    fSeries_.push_back(3.211e-4);
    fSeries_.push_back(2.702e-4);
    fSeries_.push_back(2.296e-4);
    fSeries_.push_back(1.967e-4);
    fSeries_.push_back(1.698e-4);
    
    fSeries_.push_back(1.476e-4);
    fSeries_.push_back(1.291e-4);
    fSeries_.push_back(1.136e-4);
    fSeries_.push_back(1.005e-4);
    fSeries_.push_back(8.928e-5);
    
    fSeries_.push_back(7.970e-5);
    fSeries_.push_back(7.144e-5);
    fSeries_.push_back(6.429e-5);
    fSeries_.push_back(5.806e-5);
    fSeries_.push_back(5.261e-5);
    
    fSeries_.push_back(4.782e-5);
    fSeries_.push_back(4.360e-5);
    fSeries_.push_back(3.986e-5);
    fSeries_.push_back(3.653e-5);
    fSeries_.push_back(3.357e-5);
    
    fSeries_.push_back(3.092e-5);
    fSeries_.push_back(2.854e-5);
    fSeries_.push_back(2.640e-5);
    fSeries_.push_back(2.446e-5);
    
    // size of fSeries = 39

    return;
};

void AtomicCalcs::setConsants()
{
    sigmaLymanLimitCM2_ = 6.30e-18;              //In cm^2
    freqLymanLimitInvSec_ = SPEED_OF_LIGHT_MS/WAVE_LYMANLIM_METERS;
    nLymanAlpha_ = 2;
    nLineMaxMax_ = 31;

    return;
};

double AtomicCalcs::returnWavelengthLymanSeries(int n)
{
    return WAVE_LYMANLIM_METERS/(1.-1./(n*n));
};

double AtomicCalcs::returnWavelengthAngstrLymanSeries(int n)
{
    return WAVE_LYMANLIM_ANGSTR/(1.-1./(n*n));
};

double AtomicCalcs::returnFrequencyLymanSeries(int n)
{
    double wl = returnWavelengthLymanSeries(n);
    return SPEED_OF_LIGHT_MS/wl;
};

double AtomicCalcs::returnDopplerWidthWL(int nLine, double dopplerParamKMS)
{
    double ratio = dopplerParamKMS/SPEED_OF_LIGHT_KMS;
    return returnWavelengthLymanSeries(nLine)*ratio;
};

double AtomicCalcs::returnDopplerWidthFreq(int nLine, double dopplerParamKMS)
{
    double ratio = dopplerParamKMS/SPEED_OF_LIGHT_KMS;
    return returnFrequencyLymanSeries(nLine)*ratio;
};

double AtomicCalcs::returnX(double lambda, int nLine, double dopplerParam)
{
    double dl = (lambda - returnWavelengthLymanSeries(nLine));
    double dw = returnDopplerWidthWL(nLine, dopplerParam);

    return dl/dw;
};

double AtomicCalcs::returnGamma(int nLine)
{
    int iSeries = nLine-2;
    return gammaSeries_[iSeries];
};

double AtomicCalcs::returnOscillatorStrength(int nLine)
{
    int iSeries = nLine-2;
    return fSeries_[iSeries];
};

double AtomicCalcs::returnDampingParameter(int nLine, double dopplerParamKMS)
{
    double li = returnWavelengthLymanSeries(nLine);
    double gammai = returnGamma(nLine);
    double ld = returnDopplerWidthWL(nLine, dopplerParamKMS);

    double numer = li*li*gammai;
    double denom = 4*PI*SPEED_OF_LIGHT_MS*ld;

    return numer/denom;
};

void AtomicCalcs::printEverything(int nLine, double dopplerPar)
{
    std::cout << "Starting energy level:     " << nLine << std::endl;
    std::cout << "Doppler Param:             " << dopplerPar << std::endl;

    std::cout << "Lyman Series WL in m       " << returnWavelengthLymanSeries(nLine) << std::endl;
    std::cout << "Lyman Series WL in A       " << returnWavelengthAngstrLymanSeries(nLine) << std::endl;
    std::cout << "Lyman Series freq in s^-1  " << returnFrequencyLymanSeries(nLine) << std::endl;

    std::cout << "Doppler Width WL:          " << returnDopplerWidthWL(nLine, dopplerPar) << std::endl;
    std::cout << "Doppler Width freq:        " << returnDopplerWidthFreq(nLine,dopplerPar) << std::endl;

    std::cout << "X for lambda = 500nm       " << returnX(0.0000005,nLine,dopplerPar) << std::endl;
    std::cout << "gamma                      " << returnGamma(nLine) << std::endl;
    std::cout << "Oscillator Strength        " << returnOscillatorStrength(nLine) << std::endl;
    std::cout << "Return Damping Param       " << returnDampingParameter(nLine,dopplerPar) << std::endl;

    return;
};


/******* HIColumnDensity **********************************************************/

HIColumnDensity::HIColumnDensity() 
: beta1_(1.6), beta2_(1.3), Nl_(1e12), Nc_(1.6e17), Nu_(1e22)
{
    normB_ = 1;
    normalizeDist();
    std::cout << "Normalization constant:    " << normB_ << std::endl;
    //checkNormalize();

};

/******* HIColumnDensity Methods **************************************************/

void HIColumnDensity::normalizeDist()
{
    double term1 = pow(Nc_, beta1_)*integratePowerLaw(Nl_, Nc_, -beta1_);
    double term2 = pow(Nc_, beta2_)*integratePowerLaw(Nc_, Nu_, -beta2_);
    double sumTerms = term1+term2;
    normB_ = 1/sumTerms;
    return;
};

double HIColumnDensity::integratePowerLaw(double low, double high, double power)
{
    double constant = 1/(1+power);
    double integratedTerm = pow(high, 1+power) - pow(low, 1+power);
    return constant*integratedTerm;
};

double HIColumnDensity::returnColDensityDist(double NHI)
{
    double dist = 0.0; 
    if(NHI < Nc_) {
        dist = returnFirstPowerLaw(NHI);
    }
    else if(NHI >= Nc_) {
        dist = returnSecondPowerLaw(NHI);        
    }
    else {
        std::cout << "Column Density value not recognized!" << std::endl;
    }

    return dist;
};

double HIColumnDensity::returnFirstPowerLaw(double NHI)
{
    double ratio = NHI/Nc_;
    return returnNormB()*pow(ratio,-beta1_);
};

double HIColumnDensity::returnSecondPowerLaw(double NHI)
{
        double logb, logN, logNc, logg;
        logb = log(returnNormB());
        logN = log(NHI);
        logNc = log(Nc_);
        logg = logb-(beta2_*logN)+(beta2_*logNc);
        return exp(logg);
};

double HIColumnDensity::returnNormB()
{
    return normB_;
};

void HIColumnDensity::returnPowerLawIndex(double &beta1, double &beta2)
{
    beta1 = beta1_;
    beta2 = beta2_;
    return;
};

void HIColumnDensity::returnColDensityLimits(double &Nl, double &Nu)
{
    Nl = Nl_;
    Nu = Nu_;
    return;
};

void HIColumnDensity::returnColDensityBreak(double &Nc)
{
    Nc = Nc_;
    return;
};

void HIColumnDensity::writeToFile(std::string outfile, double dLog, int nStep)
{

    double logNl = log(Nl_);

    std::ifstream inp;
	std::ofstream outp;
	inp.open(outfile.c_str(), std::ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(std::ios::failbit);
		std::cout << "     Writing HI column density distribution to file ...";
		std::cout << outfile.c_str() << std::endl;
		outp.open(outfile.c_str(), std::ofstream::out);
		
	    for (int i=0; i<nStep; i++) {
	    
			double logNh=logNl+i*dLog;
			double Nh=exp(logNh);
			outp << Nh <<"  "<< returnColDensityDist(Nh) << std::endl;
            }
            
		outp.close();
	    }
	else
		std::cout << "Error...file " << outfile.c_str() << " exists" << std::endl;

};

void HIColumnDensity::testClass()
{
    std::cout << "Confirming calculations in HIColumnDensity..." << std::endl;
    std::cout << "Normalization constant should be:         " << "2.82654e-21" << std::endl;

    double b1,b2,nl,nc,nu;
    returnPowerLawIndex(b1,b2);
    returnColDensityLimits(nl,nu);
    returnColDensityBreak(nc);

    std::cout << "beta1: " << b1 << std::endl << "beta2: " << b2 << std::endl
                << "Nl: " << nl << std::endl << "Nc: " << nc << std::endl << "Nu: " << nu << std::endl;

    return;
};

// AA: I think the below conflict should be removed? Left in just in case
/*<<<<<<< HEAD
int ProbabilityDistAbsorbers::simulateLineOfSight(double zStart,double zMax,
    vector<double>& redshifts, vector<double>& dopplerPars,
                                  vector<double>& columnDensities, string outfile)
{

    double zCurrent = zStart;
    int iAbsorber = 0;
    
    while ( zCurrent < zMax ) {
	    //cout <<"     On absorber "<<iAbsorber+1<<endl;
	    double zDraw, bdopp, NHI;
	    simulateAbsorber(zCurrent,zDraw,bdopp,NHI);
	    
	    redshifts.push_back(zDraw);
	    dopplerPars.push_back(bdopp);
	    columnDensities.push_back(NHI);	    
	    
	    zCurrent = zDraw;
	    iAbsorber++;
	    }
	int nAbsorbers = columnDensities.size();
    cout <<"     Total number of absorbers along line of sight = "<<nAbsorbers<<".";
	cout <<" Max redshift = "<<redshifts[redshifts.size()-1]<<endl;
	
	writeToFile(outfile, redshifts, dopplerPars, columnDensities);

    return nAbsorbers;
=======*/
/******* AbsorberRedshiftDistribution *********************************************/
AbsorberRedshiftDistribution::AbsorberRedshiftDistribution() 
: A_(400.), z1_(1.2), z2_(4.0), gamma1_(0.2), gamma2_(2.5), gamma3_(4.0)
{

//>>>>>>> b174eb658571fc61c81d79c7fa2a0db330092148
};

/******* AbsorberRedshiftDistribution methods *************************************/
double AbsorberRedshiftDistribution::returnRedshiftDist(double z)
{
    double dist = 0.0;
    if(z <= z1_) {
        dist = returnNormalization()*returnFirstPowerLaw(z);
    }
    else if(z <= z2_ && z > z1_) {
        dist = returnNormalization()*returnSecondPowerLaw(z);
    }
    else if(z > z2_) {
        dist = returnNormalization()*returnThirdPowerLaw(z);
    }
    else {
        std::cout << "Redshift not recognized!" << std::endl;
    }

    return dist;
};

double AbsorberRedshiftDistribution::returnFirstPowerLaw(double z)
{
    double ratio = 0.0;
    ratio = (1+z)/(1+z1_);
    return pow(ratio, gamma1_);
};

double AbsorberRedshiftDistribution::returnSecondPowerLaw(double z)
{
    double ratio = 0.0;
    ratio = (1+z)/(1+z1_);
    return pow(ratio, gamma2_);
};

double AbsorberRedshiftDistribution::returnThirdPowerLaw(double z)
{
    double ratio2 = 0.0, ratio3 = 0.0;
    ratio2 = (1+z2_)/(1+z1_);
    ratio3 = (1+z)/(1+z2_);
    return pow(ratio2,gamma2_)*pow(ratio3,gamma3_);
};

void AbsorberRedshiftDistribution::returnPowerLawIndex(double &g1, double &g2, double &g3)
{
    g1 = gamma1_;
    g2 = gamma2_;
    g3 = gamma3_;
    return;
};

void AbsorberRedshiftDistribution::returnRedshiftBreaks(double &z1, double &z2)
{
    z1 = z1_;
    z2 = z2_;
    return;
};

double AbsorberRedshiftDistribution::returnNormalization()
{
    return A_;
};

void AbsorberRedshiftDistribution::testClass()
{
    double g1, g2, g3, z1, z2, A;
    returnPowerLawIndex(g1,g2,g3);
    returnRedshiftBreaks(z1,z2);
    A = returnNormalization();

    std::cout << "gamma1 " << g1 << std::endl
                << "gamma2 " << g2 << std::endl
                << "gamma3 " << g3 << std::endl
                << "z1     " << z1 << std::endl
                << "z2     " << z2 << std::endl
                << "A      " << A  << std::endl;
    std::cout << std::endl;

    return;
};

/******* DopplerParDistribution ***************************************************/
DopplerParDistribution::DopplerParDistribution()
: bsigma_(23)
{

};

/******* DopplerParDistribution methods *******************************************/

double DopplerParDistribution::returnDopplerDist(double b)
{
    double b4, bsig4, ratio4, p1, p2;
    b4 = b*b*b*b;
    bsig4 = bsigma_*bsigma_*bsigma_*bsigma_;
    ratio4 = bsig4/b4;
    p1 = 4*ratio4/b;
    p2 = exp(-ratio4);
    return p1*p2;
};


void DopplerParDistribution::testClass()
{
    std::cout << "bsigma " << bsigma_ << std::endl;
    return;
};




/******* ProbabilityDistAbsorbers *************************************************/
ProbabilityDistAbsorbers::ProbabilityDistAbsorbers(RandomGeneratorInterface& rg,
                                                   AbsorberRedshiftDistribution& absorberZDist,
                                                   HIColumnDensity& hiColumnDensity,
                                                   DopplerParDistribution& dopplerParDist)
: rg_(rg), absorberZDist_(absorberZDist), hiColumnDensity_(hiColumnDensity), dopplerParDist_(dopplerParDist)

{
    int nStep = 1000;
    setNHiDistribution(nStep);
    bmin_ = 0.0;
    bmax_ = 200.;
    setDopplerDistribution(nStep);

};

/******* ProbabilityDistAbsorbers methods *****************************************/

void ProbabilityDistAbsorbers::setNHiDistribution(int nStep)
{
    log10NHIvals_.clear();
    log10gvals_.clear();

    double Nl, Nu;
    hiColumnDensity_.returnColDensityLimits(Nl,Nu);
    log10Nl_ = log10(Nl);
    log10Nu_ = log10(Nu);

    double dlogN = (log10Nu_ - log10Nl_)/(nStep-1);

    for(int i=0; i<nStep; i++) {
        double log10NHI = log10Nl_ + i*dlogN;
        log10NHIvals_.push_back(log10NHI);

        double NHI = pow(10., log10NHI);
        double g = hiColumnDensity_(NHI);
        log10gvals_.push_back(log10(g));
    }

    log10gmin_ = findMinimum(log10gvals_);
    log10gmax_ = findMaximum(log10gvals_);

    // Create interpolation
//    colDensityFunc_.DefinePoints(log10NHIvals_, log10gvals_);

    return;
  
};

void ProbabilityDistAbsorbers::setDopplerDistribution(int nStep)
{
    bvals_.clear();
    hvals_.clear();

    double db = (bmax_-bmin_)/(nStep-1);

    for(int i=0; i<nStep; i++) {
        double bv = bmin_ + i*db;
        bvals_.push_back(bv);

        double hv = dopplerParDist_(bv);
        hvals_.push_back(hv);
    }

    hmin_ = findMinimum(hvals_);
    hmax_ = findMaximum(hvals_);

    return;
};

double ProbabilityDistAbsorbers::drawDeltaZ(double zLast)
{
    double fz = absorberZDist_(zLast);
    double rn = rg_.Flat01();

    // p(\DeltaZ;z) = f(z)*exp( -f(z)*\DeltaZ )
    double diff = -log(rn);
    double deltaZ = diff/fz;

    return deltaZ;
};

double ProbabilityDistAbsorbers::drawHIColumnDensity()
{
/*    double u1, u2, nhi;

    while(true) {
        u1 = log10Nl_ + (log10Nu_ - log10Nl_)*rg_.Flat01();
        u2 = log10gmin_ + (log10gmax_ - log10gmin_)*rg_.Flat01();

        nhi = pow(10.,u1);
        double gval = hiColumnDensity_(nhi);
        double log10gval = log10(gval);

        // Use interpolation to draw the number
        //double log10gval = colDensityFunc_(u1);

        if(u2 <= log10gval) {
            break;
        }
    }

    return pow(10., u1);        */

    //IMPLEMENTING THE TRANSFORMATION PERFORMED IN CODE BY
    //INOOUE & IWATA FOR THEIR 2008 PUBLICATION

    double beta1, beta2, Nl, Nc, Nu;
    hiColumnDensity_.returnPowerLawIndex(beta1, beta2);
    hiColumnDensity_.returnColDensityLimits(Nl, Nu);
    hiColumnDensity_.returnColDensityBreak(Nc);

    double nhi = 0.0;
    double index1 = 1.0 / (1.0 - beta1);
    double index2 = 1.0 / (1.0 - beta2);
    double Rlc = pow((Nl/Nc),(1.0-beta1));
    double Ruc = pow((Nu/Nc),(1.0-beta2));
    double xc = (Rlc - 1.0) / (Rlc - 1.0 + (1.0 - Ruc) * ((beta1 - 1.0) / (beta2 - 1.0)));

    double ran = rg_.Flat01();

    if(ran < xc) {
        nhi = Nc * pow((((xc-ran)*Rlc + ran) / (xc)), index1);
    } else {
        nhi = Nc * pow(((1.0 - ran + (ran-xc)*Ruc) / (1.0 - xc)), index2);
    }

    return nhi;

};

double ProbabilityDistAbsorbers::drawDoppler()
{
    double u1, u2;

    while(true) {
        u1 = bmin_ + (bmax_ - bmin_)*rg_.Flat01();
        u2 = hmin_ + (hmax_ - hmin_)*rg_.Flat01();

        double hval = dopplerParDist_(u1);

        if(u2 <= hval) {
            break;
        }
    }

    return u1;
};

void ProbabilityDistAbsorbers::simulateLineOfSight(double zStart, double zMax,
                    vector<double>& redshifts, vector<double>& dopplerPars,
                    vector<double>& columnDensities, string outfile)
{
    double zCurrent = zStart;
    int iAbsorber = 0;

    while(zCurrent < zMax) {
        double zNext, bdopp, NHI;
        simulateAbsorber(zCurrent,zNext,bdopp,NHI);

        redshifts.push_back(zNext);
        dopplerPars.push_back(bdopp);
        columnDensities.push_back(NHI);

        zCurrent = zNext;
        iAbsorber++;
    }

    cout << "       The total number of absorbers along the line of sight = "
            << iAbsorber << endl;
    cout << "       Max redshift = " << redshifts[iAbsorber-1] << endl;

    return;
};

void ProbabilityDistAbsorbers::simulateAbsorber(double zCurrent, double& zNext,
                    double& bdopp, double& NHI)
{
    zNext = zCurrent + drawDeltaZ(zCurrent);
    bdopp = drawDoppler();
    NHI = drawHIColumnDensity();

    return;
};


/******* VoigtProfile *************************************************************/
VoigtProfile::VoigtProfile(double dopplerPar, int nLine)
: dopplerPar_(dopplerPar), nLine_(nLine)
{
    double nGammaMax = gammaSeries_.size() + 1;
    if(nLine > nGammaMax) {
        cout << "STARTING LINE IS TOO LARGE, NO GAMMA VALUES FOR THAT TRANSITION"
                << endl;
    }

};

/******* VoigtProfile Methods *****************************************************/
double VoigtProfile::kFunction(double x)
{
    double prefact = 1./(2*x*x);
    double t1 = (4*x*x + 3.) * (x*x + 1.) * exp(-x*x);
    double t2 = (1./(x*x)) * (2*x*x + 3.) * sinh(x*x);
    double k = prefact*(t1 - t2);

    int isInf = my_isinf(k);
    if(isInf != 0) {
//        k = 0.;
        k = 1000000.;
    }

    return k;
};

double VoigtProfile::returnHax(double lambda)
{
    double x = returnX(lambda, nLine_, dopplerPar_);
    double K = kFunction(x);
    double a = returnDampingParameter(nLine_, dopplerPar_);

    double H1 = exp(-x*x) * (1. - a*2.*K/sqrt(PI));

    return H1;
};


/******* OpticalDepth *************************************************************/

OpticalDepth::OpticalDepth()
{
    nLineMax_ = nLineMaxMax_;
    setLymanAll();

};


/******* OpticalDepth Methods *****************************************************/

void OpticalDepth::setLymanAll()
{
    setContribution(0);
    return;
};

void OpticalDepth::setLymanContinuumOnly()
{
    setContribution(1);
    return;
};

void OpticalDepth::setLymanSeriesOnly()
{
    setContribution(2);
    return;
};

void OpticalDepth::setContribution(int option)
{
    isLymanC_ = true;
    isLymanS_ = true;

    if(option == 1) {
        isLymanS_ = false;
    } else if(option == 2) {
        isLymanC_ = false;
    }

    return;
};

double OpticalDepth::returnObserverFrameTransmission(double lambda, double zAbsorber,
                    double nhiAbsorber, double bAbsorber)
{
    double tau = returnObserverFrameOpticalDepth(lambda, zAbsorber, nhiAbsorber, bAbsorber);
    return exp(-tau);
};

double OpticalDepth::returnObserverFrameOpticalDepth(double lambda, double zAbsorber,
                    double nhi, double bAbsorber)
{
    double lambdaE = lambda/(1.0+zAbsorber);
    double freq = SPEED_OF_LIGHT_MS/lambdaE;
    double tau = returnRestFrameOpticalDepth(freq, bAbsorber, nhi);

    return tau;
};

double OpticalDepth::returnRestFrameOpticalDepth(double freq, double bAbsorber, double nhi)
{
    double sigmaLC = 0.0, sumSigmai = 0.0, sigmaT = 0.0, restFrameOpticalDepth = 0.0;

    sigmaLC = returnLymanContinuumCrossSection(freq);
    sumSigmai = returnLymanSeriesCrossSection(freq, bAbsorber);

    if(isLymanC_)
        sigmaT += sigmaLC;
    if(isLymanS_)
        sigmaT += sumSigmai;

    restFrameOpticalDepth = nhi*sigmaT;

    return restFrameOpticalDepth;
};

double OpticalDepth::returnLymanContinuumCrossSection(double freq)
{
    double sigmaLC = 0.0;

    if(freq >= freqLymanLimitInvSec_) {
        double ratio = freqLymanLimitInvSec_/freq;
        sigmaLC = sigmaLymanLimitCM2_*ratio*ratio*ratio;
    } else {
        sigmaLC = 0.0;
    }

    return sigmaLC;
};

double OpticalDepth::returnLymanSeriesCrossSection(double freq, double bAbsorber)
{
    double sigmaSumLymanSeries = 0.0;
    for(int n=nLymanAlpha_; n<=nLineMax_; n++) {
        sigmaSumLymanSeries += returnLymanLineCrossSection(n, freq, bAbsorber);
    }

    return sigmaSumLymanSeries;
};

double OpticalDepth::returnLymanLineCrossSection(int n, double freq, double bAbsorber)
{
    double constNumer = sqrt(PI)*ELECTRON_CHARGE_STATC*ELECTRON_CHARGE_STATC;
    double constDenom = ELECTRON_MASS_G*SPEED_OF_LIGHT_CMS;

    double vd = returnDopplerWidthFreq(n, bAbsorber);
    double fi = returnOscillatorStrength(n);
    double phii = returnLineProfile(n, freq, bAbsorber);

    double numer = constNumer*fi*phii;
    double denom = constDenom*vd;

    double sigmai = numer/denom;

    return sigmai;
};

double OpticalDepth::returnLineProfile(int n, double freq, double bAbsorber)
{
    VoigtProfile voigtProf(bAbsorber, n);
    double lambda = SPEED_OF_LIGHT_MS/freq;

    double vprof = voigtProf(lambda);
    
    return vprof;
};

void OpticalDepth::setMaxLine(int nLine)
{
    nLineMax_ = nLine;
    return;
};


/******* LineOfSightTrans *********************************************************/
LineOfSightTrans::LineOfSightTrans(vector<double>& redshifts, vector<double>& dopplerPars,
                                   vector<double>& columnDensities)
: redshifts_(redshifts), dopplerPars_(dopplerPars), columnDensities_(columnDensities)
{
    setLymanAll();
    isOpticalDepth_=false;
};

/******* LineOfSightTrans Methods *************************************************/

double LineOfSightTrans::returnOpticalDepth(double lambda, double zSource)
{
//    int numAbsorbers = returnNumberOfAbsorbers(zSource);
    int numAbsorbers = redshifts_.size();
    double tau = 0.0;

    for(int i=0; i<numAbsorbers; i++) {
        double z = redshifts_[i];
        double nhi = columnDensities_[i];
        double b = dopplerPars_[i];

        if(z > zSource)
            break;

        double tauAbsorber = returnObserverFrameOpticalDepth(lambda, z, nhi, b);
        
        tau += tauAbsorber;
    }

    return tau; 
};

double LineOfSightTrans::returnTransmission(double lambda, double zSource)
{
    double tau = returnOpticalDepth(lambda, zSource);
    return exp(-tau);
};

int LineOfSightTrans::returnNumberOfAbsorbers(double zSource)
{
    int numAbsorbers = 0;
    int nz = redshifts_.size();
    for(int i=0; i<nz; i++) {
        if(redshifts_[i] < zSource)
            numAbsorbers++;
    }

    return numAbsorbers;
};

void LineOfSightTrans::setReturnType(bool isOpticalDepth)
{
    isOpticalDepth_ = isOpticalDepth;
    return;
};

















/******* Madau methods ********************************************************/

double Madau::returnObserverFrameOpticalDepth(double lambdaObs, double zSource)
{
	// Put the wavelength into the emission frame
	double lambdaEm = lambdaObs / (1 + zSource);
	
	// this variable is set but not used
	int cas = 0;  // Case whether pure Lyman series or Lyman series+continuum
	
	double lambdaAlpha = returnWavelengthLymanSeries(2);
	double lambdaLineMin = returnWavelengthLymanSeries(nLineMax_);

	// Initialize optical depth value as transparent
	double tau = 0;

	if (lambdaObs > lambdaAlpha)
		tau = 0;//return_value = 1;
		
    // If emission wavelength is between the lyman limit and lyman-alpha
	if (WAVE_LYMANLIM_METERS <= lambdaEm && lambdaEm <= lambdaAlpha) {// eq 15
	
	    // If emission wavelength is between lyman limit and lyman-highest
		if (lambdaEm < lambdaLineMin ) {
		
		    for (int n=2; n<=nLineMax_; n++) { // do whole series absorption
		    
		        double lambdan = returnWavelengthLymanSeries(n);
				tau += Avals_[n-2]*pow((lambdaObs/lambdan), 3.46);
				
			    }
		    }
		// If emission wavelength is between lyman-highest and lyman-alpha
        else  {
        
			for (int n=3; n<=nLineMax_; n++) {
			
			    // If emission wavelength is between line n and line n-1 
				if (returnWavelengthLymanSeries(n)<= lambdaEm && lambdaEm <= returnWavelengthLymanSeries(n-1)) {
				
				    // loop from lyman-alpha up to line n
					for (int j=2; j<n; j++) {
					    double lambdan = returnWavelengthLymanSeries(j);
						tau += Avals_[j-2] * pow((lambdaObs/lambdan), 3.46);

					    }
				    }
			    }
		    }

		cas = 1;
	}
	// If emission wavelength is BELOW lyman limit
	else if (lambdaEm <= WAVE_LYMANLIM_METERS && isLyC_) { //eq 16

        // Do Lyman series first
		for (int n = 2; n<=nLineMax_; n++) {
		    double lambdan = returnWavelengthLymanSeries(n);
		    tau += Avals_[n-2]*pow((lambdaObs/lambdan), 3.46);
		    }

        // Now do photoelectric absorption   
		
		// Minimum z an absorber could have
		double zmin = lambdaObs / WAVE_LYMANLIM_METERS - 1;
		if (zmin < 0)
			zmin = 0.;
        
        tau += returnLymanContinuumOpticalDepth(zmin,zSource);

		cas = 2;
	    }

	return tau;
 
};

double Madau::returnLymanContinuumOpticalDepth(double zAbsorber, double zSource) 
{
// The approximate integration of equation 16 in Madau 1995

    double tau = 0;
    double xe = 1 + zSource;
    double xc = 1 + zAbsorber;
    double term1 = 0.25 * pow(xc, 3) * (pow(xe, 0.46) - pow(xc, 0.46));
    double term2 = 9.4 * pow(xc, 1.5) * (pow(xe, 0.18) - pow(xc, 0.18));
    double term3 = -0.7 * pow(xc, 3) * (pow(xc, -1.32) - pow(xe, -1.32));
    double term4 = -0.023 * (pow(xe, 1.68) - pow(xc, 1.68));
    tau = term1 + term2 + term3 + term4;

	return tau;
}

void Madau::setAbsorptionStrengths()
{

    Avals_.push_back(0.0036);
    Avals_.push_back(1.7e-3);
    Avals_.push_back(1.2e-3);
    Avals_.push_back(9.3e-4);
    nLineMaxMadau_=Avals_.size()+1;

};

/******* LAFMeiksin methods ********************************************************/

double LAFMeiksin::returnTauLymanSeries(double lambdaObs, double zSource)
{
    // see Meiksin 2006 
    // The Lyman series is the series of transitions and resulting 
    // ultraviolet emission lines of the hydrogen atom as an electron goes 
    // from n ≥ 2 to n = 1 
    // Therefore Lyman-alpha is the transition from n=2:

    // total transmission due to all lines
	double tauTotal = 0; // begin by setting to zero

    // getting Lyman-alpha line contribution
    // zSource: redshift *LIGHT* is emitted from
    // zLyAlpha: redshift of Lyman-alpha if observed at wavelength lambdaObs
    double lambdaLyAlpha = returnWavelengthLymanSeries(2);
	double zLyAlpha = lambdaObs/lambdaLyAlpha - 1.;
	double tauLyA = 0;
    tauLyA = returnTauLyA(zLyAlpha, zSource);
	tauTotal += tauLyA;

	//if (prt_>0)
	//	cout << tauTotal <<"  "<< ratio << endl;
	
    // work out how many Lyman series are needed
    int nmax = findNmax (lambdaObs, zSource);
	
	// getting Lyman-series line contributions
	// zLyn: redshift of Lyman-line if observed at wavelength lambdaObs
    double zLyn, lambdaLyn, tauLyn;

    // The next transition is from n=3+
	for (int n = 3; n <= nmax; n++) {
	
	    // wavelength in meters of current Lyman series line
	    lambdaLyn = returnWavelengthLymanSeries(n);
	    
	    // redshift of Lyman-n absorber that would be observable at lambdaObs
        zLyn = lambdaObs/lambdaLyn - 1.;
        
        // find Lyman alpha transmission at the redshift of this Lyman line absorber
	    tauLyA = 0;
        tauLyA = returnTauLyA(zLyn, zSource);
        
        // useful quantity
        double zVal = returnZvalue(zLyn,n);
        
	    double tauRatio;
	    if (n <=9) {
	    
	        // get ratio of tauLyn/tauLyA
	        tauRatio = tauCoeff_[n-3]*zVal;
	        
	        if (prt_>0)
                cout << tauTotal <<"  "<< n << endl;

            }
        else {

            // get ratio of tauLyn/tauLyA
            tauRatio = 0.028266 * (720./(n*(n*n - 1.)))*zVal; // eq. 4

            if (prt_>0)
                cout << tauTotal <<"  "<< n << endl;
			}
	    // multiply by tauLyA 
        tauLyn = tauRatio*tauLyA;
			
        // add to total
        tauTotal += tauLyn;

	    }

	return tauTotal;

};

void LAFMeiksin::setCoeffs()
{
    // see Meiksin 2006 Table 1
    tauCoeff_.push_back(0.348);
    tauCoeff_.push_back(0.179);
    tauCoeff_.push_back(0.109);
    tauCoeff_.push_back(0.0722);
    tauCoeff_.push_back(0.0508);
    tauCoeff_.push_back(0.0373);
    tauCoeff_.push_back(0.0283);
    
};

double LAFMeiksin::taulya(double lambda, double z) 
{
		/* add Lyman series resonant scattering */

		//int debug = 0;
		int nmax, n;
		double lambdaL = 911.75;
		double tauLyn;
		double tautot = 0;

		double lambdaLyn = lambdaL / (1. - 1. / 4.);
		double zLyn = lambda / lambdaLyn - 1.;
		double tauLyA = 0;
		double ratio = lambda / lambdaL;

		if (zLyn < z) {
			if (zLyn < 4.) {
				if (zLyn < 1.2)
					tauLyA = 0.0164 * pow((1. + zLyn), 1.1);
				else
					tauLyA = 0.00211 * pow((1. + zLyn), 3.7);
			} else
				tauLyA = 0.00058 * pow((1. + zLyn), 4.5);
		}
		tautot += tauLyA;
		cout <<"     lambda = "<< lambda <<", lambdaLyn = "<< lambdaLyn <<endl;
		cout <<"     zLyAlpha = " << zLyn << " zSource = "<< z <<endl;
		cout <<"     Value added to tautot = " << tauLyA <<", n = 2"<<endl;

		//if (debug == 1) {
		//	System.out.println(tautot + " " + ratio);
		//}

		if (ratio < (1. + z)) {
			nmax = 31; /* add full Lyman series, cut at n = 31 */
		} else
			nmax = (int) pow((1. - (1. + z) / ratio), -0.5);
		if (nmax > 32) {
			nmax = 31;
		}
        cout <<"     nmax = "<< nmax <<endl;
		for (n = 3; n <= nmax; n++) {
		
		
			switch (n) {

			case 3:
				lambdaLyn = lambdaL / (1. - 1. / 9.);
				zLyn = lambda / lambdaLyn - 1.;
				tauLyA = 0;
				if (zLyn < z) {
					if (zLyn < 4.) {
						if (zLyn < 1.2)
							tauLyA = 0.0164 * pow((1. + zLyn), 1.1);
						else
							tauLyA = 0.00211 * pow((1. + zLyn), 3.7);
					} else
						tauLyA = 0.00058 * pow((1. + zLyn), 4.5);
				}
				if (zLyn < 3) {
					tauLyn = 0.34827 * pow(0.25 * (1. + zLyn), (1. / 3.));
				} else
					tauLyn = 0.34827 * pow(0.25 * (1. + zLyn), (1. / 6.));
				tautot += tauLyn * tauLyA;
				cout <<"     Value added to tautot = " << (tauLyn * tauLyA) <<", n = "<< n <<", tauLyA = "<<tauLyA << endl;
				//if (debug == 1) {
				//	System.out.println(tautot + " " + n);
				//}
				

				break;

			case 4:
				lambdaLyn = lambdaL / (1. - 1. / 16.);
				zLyn = lambda / lambdaLyn - 1.;
				tauLyA = 0;
				if (zLyn < z) {
					if (zLyn < 4.) {
						if (zLyn < 1.2)
							tauLyA = 0.0164 * pow((1. + zLyn), 1.1);
						else
							tauLyA = 0.00211 * pow((1. + zLyn), 3.7);
					} else
						tauLyA = 0.00058 * pow((1. + zLyn), 4.5);
				}
				if (zLyn < 3) {
					tauLyn = 0.17932 * pow(0.25 * (1. + zLyn), (1. / 3.));
				} else
					tauLyn = 0.17932 * pow(0.25 * (1. + zLyn), (1. / 6.));
				tautot += tauLyn * tauLyA;
				cout <<"     Value added to tautot = " << (tauLyn * tauLyA) <<", n = "<< n <<", tauLyA = "<<tauLyA << endl;
				//if (debug == 1) {
				//	System.out.println(tautot + " " + n);
				//}
				break;

			case 5:
				lambdaLyn = lambdaL / (1. - 1. / 25.);
				zLyn = lambda / lambdaLyn - 1.;
				tauLyA = 0;
				if (zLyn < z) {
					if (zLyn < 4.) {
						if (zLyn < 1.2)
							tauLyA = 0.0164 * pow((1. + zLyn), 1.1);
						else
							tauLyA = 0.00211 * pow((1. + zLyn), 3.7);
					} else
						tauLyA = 0.00058 * pow((1. + zLyn), 4.5);
				}
				tauLyn = 0.10872 * pow(0.25 * (1. + zLyn), (1. / 3.));
				tautot += tauLyn * tauLyA;
				cout <<"     Value added to tautot = " << (tauLyn * tauLyA) <<", n = "<< n <<", tauLyA = "<<tauLyA << endl;
				//if (debug == 1) {
				//	System.out.println(tautot + " " + n);
				//}
				break;

			case 6:
				lambdaLyn = lambdaL / (1. - 1. / 36.);
				zLyn = lambda / lambdaLyn - 1.;
				tauLyA = 0;
				if (zLyn < z) {
					if (zLyn < 4.) {
						if (zLyn < 1.2)
							tauLyA = 0.0164 * pow((1. + zLyn), 1.1);
						else
							tauLyA = 0.00211 * pow((1. + zLyn), 3.7);
					} else
						tauLyA = 0.00058 * pow((1. + zLyn), 4.5);
				}
				tauLyn = 0.072242 * pow(0.25 * (1. + zLyn), (1. / 3.));
				tautot += tauLyn * tauLyA;
				cout <<"     Value added to tautot = " << (tauLyn * tauLyA) <<", n = "<< n <<", tauLyA = "<<tauLyA << endl;
				//if (debug == 1) {
				//	System.out.println(tautot + " " + n);
				//}
				break;

			case 7:
				lambdaLyn = lambdaL / (1. - 1. / 49.);
				zLyn = lambda / lambdaLyn - 1.;
				tauLyA = 0;
				if (zLyn < z) {
					if (zLyn < 4.) {
						if (zLyn < 1.2)
							tauLyA = 0.0164 * pow((1. + zLyn), 1.1);
						else
							tauLyA = 0.00211 * pow((1. + zLyn), 3.7);
					} else
						tauLyA = 0.00058 * pow((1. + zLyn), 4.5);
				}
				tauLyn = 0.050852 * pow(0.25 * (1. + zLyn), (1. / 3.));
				tautot += tauLyn * tauLyA;
				cout <<"     Value added to tautot = " << (tauLyn * tauLyA) <<", n = "<< n <<", tauLyA = "<<tauLyA << endl;
				//if (debug == 1) {
				//	System.out.println(tautot + " " + n);
				//}
				break;

			case 8:
				lambdaLyn = lambdaL / (1. - 1. / 64.);
				zLyn = lambda / lambdaLyn - 1.;
				tauLyA = 0;
				if (zLyn < z) {
					if (zLyn < 4.) {
						if (zLyn < 1.2)
							tauLyA = 0.0164 * pow((1. + zLyn), 1.1);
						else
							tauLyA = 0.00211 * pow((1. + zLyn), 3.7);
					} else
						tauLyA = 0.00058 * pow((1. + zLyn), 4.5);
				}
				tauLyn = 0.037336 * pow(0.25 * (1. + zLyn), (1. / 3.));
				tautot += tauLyn * tauLyA;
				cout <<"     Value added to tautot = " << (tauLyn * tauLyA) <<", n = "<< n <<", tauLyA = "<<tauLyA << endl;
				//if (debug == 1) {
				//	System.out.println(tautot + " " + n);
				//}
				break;

			case 9:
				lambdaLyn = lambdaL / (1. - 1. / 81.);
				zLyn = lambda / lambdaLyn - 1.;
				tauLyA = 0;
				if (zLyn < z) {
					if (zLyn < 4.) {
						if (zLyn < 1.2)
							tauLyA = 0.0164 * pow((1. + zLyn), 1.1);
						else
							tauLyA = 0.00211 * pow((1. + zLyn), 3.7);
					} else
						tauLyA = 0.00058 * pow((1. + zLyn), 4.5);
				}
				tauLyn = 0.028266 * pow(0.25 * (1. + zLyn), (1. / 3.));
				tautot += tauLyn * tauLyA;
				cout <<"     Value added to tautot = " << (tauLyn * tauLyA) <<", n = "<< n <<", tauLyA = "<<tauLyA << endl;
				//if (debug == 1) {
				//	System.out.println(tautot + " " + n);
				//}
				break;

			default:
				lambdaLyn = lambdaL / (1. - 1. / (n * n));
				zLyn = lambda / lambdaLyn - 1.;
				tauLyA = 0;
				if (zLyn < z) {
					if (zLyn < 4.) {
						if (zLyn < 1.2)
							tauLyA = 0.0164 * pow((1. + zLyn), 1.1);
						else
							tauLyA = 0.00211 * pow((1. + zLyn), 3.7);
					} else
						tauLyA = 0.00058 * pow((1. + zLyn), 4.5);
				}
				tauLyn = 0.028266 * (720. / (n * (n * n - 1.)))
						* pow(0.25 * (1. + zLyn), (1. / 3.));
				tautot += tauLyn * tauLyA;
				cout <<"     Value added to tautot = " << (tauLyn * tauLyA) <<", n = "<< n <<", tauLyA = "<<tauLyA << endl;
				//if (debug == 1) {
				//	System.out.println(tautot + " " + n);
				//}
				break;
			}
		}

		// System.out.println(nmax);
		return exp(-tautot);

}
/******* MonteCarloMeiksin methods ********************************************/

double MonteCarloMeiksin::returnOpticaldepth(double lambdaObs, double zSource)
{

    // Transmission that isn't Monte Carlo
    double tauLAF = returnTauLymanSeries(lambdaObs,zSource);// Lyman-alpha forest
    double taudIGM = returnTauDiffuse(lambdaObs,zSource);   // diffuse IGM
    
    // Transmission calculated from line of sight distribution
    double tauLLS=returnTauLymanLimitSystems(lambdaObs, zSource);

    double tauTotal = tauLLS + taudIGM + tauLAF;
    return tauTotal; 
}

double MonteCarloMeiksin::returnTauLymanLimitSystems(double lambdaObs, double zSource)
{

    int nAbsorber = zAbsorbers_.size();
    
    double tauLLS=0.;
    for (int i=0; i<nAbsorber; i++) {
    
        double zA = zAbsorbers_[i];
        double tauVal = opticalDepths_[i];
        tauLLS += returnTauLLS(lambdaObs, tauVal, zA);
    
        }
    return tauLLS;

};

