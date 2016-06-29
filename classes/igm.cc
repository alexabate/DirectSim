#include "igm.h"


/*  */

/******* AtomicCalcs **********************************************************/

AtomicCalcs::AtomicCalcs()
{
  //    setGammas();
    //    setOscillatorStrength();
    setConstants();
};

void AtomicCalcs::setConstants()
{
    sigmaLymanLimitCM2_ = 6.30e-18;              //In cm^2
    freqLymanLimitInvSec_ = SPEED_OF_LIGHT_MS/WAVE_LYMANLIM_METERS;
    nLymanAlpha_ = 2;
    nLineMaxMax_ = 40;   //changed to 40 from 32
    nGammaMax_ = 40;     //changed to 40 from 24
    R= 1.34;
    s= 2.99;

    return;
};


double AtomicCalcs::returnWavelengthLymanSeries(int n)
{
  static std::vector<double> wavelength_lyman_;
  static int init = 1;

  // initialize lookup table
  if (init) {
    init = 0;
    for (int i=0; i<nLineMaxMax_; i++) {
      int k = i+1;
      if (i==0) 
	wavelength_lyman_.push_back(0);
      else
	wavelength_lyman_.push_back(WAVE_LYMANLIM_METERS/(1.-1./(k*k)));
    }
  }
  
  //  return WAVE_LYMANLIM_METERS/(1.-1./(n*n));
  return wavelength_lyman_[n-1];
};


double AtomicCalcs::returnGamma(int nLine)
{
  static double gammaSeries_[] = {
	4.67e8,   //Taken from Inoue's 2008 MC
	9.93e7,
	3.00e7,
	1.15e7,
	5.17e6,
	2.60e6,
	1.43e6,
	8.40e5,
	5.21e5,
	3.37e5,
	2.26e5,
	1.57e5,
	1.11e5,
	8.3e4,
	6.0e4,
	4.4e4,
	3.3e4,
	2.5e4,
	2.0e4,
	1.5e4,
	1.2e4,
	9.8e3,
	7.9e3,
	6.4e3,
	5.2e3,
	4.4e3,
	3.6e3,
	3.0e3,
	2.5e3,
	2.1e3,
	1.9e3,
	1.6e3,
	1.4e3,
	1.2e3,
	1.0e3,
	9.1e2,
	7.9e2,
	7.0e2,
	6.1e2,
};
    int iSeries = nLine-2;
    return gammaSeries_[iSeries];
};

double AtomicCalcs::returnOscillatorStrength(int nLine)
{
  static double fSeries_[] = {
    0.4162, 
    7.910e-2,  
    2.899e-2,  
    1.394e-2,
    7.799e-3,
    4.814e-3,
    3.183e-3,
    2.216e-3,
    1.605e-3,
    1.201e-3,
    9.214e-4,
    7.227e-4,
    5.774e-4,
    4.686e-4,
    3.856e-4,
    3.211e-4,
    2.702e-4,
    2.296e-4,
    1.967e-4,
    1.698e-4,
    1.476e-4,
    1.291e-4,
    1.136e-4,
    1.005e-4,
    8.928e-5,
    7.970e-5,
    7.144e-5,
    6.429e-5,
    5.806e-5,
    5.261e-5,
    4.782e-5,
    4.360e-5,
    3.986e-5,
    3.653e-5,
    3.357e-5,
    3.092e-5,
    2.854e-5,
    2.640e-5,
    2.446e-5};

    int iSeries = nLine-2;
    return fSeries_[iSeries];
};

/****************
double AtomicCalcs::returnDampingParameter(int nLine, double dopplerParamKMS)
{
    double li = returnWavelengthLymanSeries(nLine);
    double gammai = returnGamma(nLine);
    double ld = returnDopplerWidthWL(nLine, dopplerParamKMS);

    double numer = li*li*gammai;
    double denom = 4*PI*SPEED_OF_LIGHT_MS*ld;

    return numer/denom;
};
*************/

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


/******* HIColumnDensityLAF **********************************************************/

HIColumnDensityLAF::HIColumnDensityLAF() 
: betaLAF_(1.7), Nl_(1e12), Nc_(1e21), Nu_(1e23)
{
//put normalization here if needed
};

void HIColumnDensityLAF::returnPowerLawIndex(double &betaLAF)
{
    betaLAF = betaLAF_;
    return;
};

void HIColumnDensityLAF::returnColDensityLimits(double &Nl, double &Nu)
{
    Nl = Nl_;
    Nu = Nu_;
    return;
};

void HIColumnDensityLAF::returnColDensityBreak(double &Nc)
{
    Nc = Nc_;
    return;
};

void HIColumnDensityLAF::testClass()
{
    std::cout << "Confirming calculations in HIColumnDensity..." << std::endl;
    std::cout << "Normalization constant should be:         " << "something I've haven't yet decided to include or not" << std::endl;

    double bLAF, nl,nc,nu;
    returnPowerLawIndex(bLAF);
    returnColDensityLimits(nl,nu);
    returnColDensityBreak(nc);

    std::cout << "betaLAF: " << bLAF << std::endl
                << "Nl: " << nl << std::endl << "Nc: " << nc << std::endl << "Nu: " << nu << std::endl;

    return;
};




/******* HIColumnDensityDLA **********************************************************/

HIColumnDensityDLA::HIColumnDensityDLA() 
: betaDLA_(0.9), Nl_(1e12), Nc_(1e21), Nu_(1e23)
{
//put normalization here if needed
};

void HIColumnDensityDLA::returnPowerLawIndex(double &betaDLA)
{
    betaDLA = betaDLA_;
    return;
};

void HIColumnDensityDLA::returnColDensityLimits(double &Nl, double &Nu)
{
    Nl = Nl_;
    Nu = Nu_;
    return;
};

void HIColumnDensityDLA::returnColDensityBreak(double &Nc)
{
    Nc = Nc_;
    return;
};

void HIColumnDensityDLA::testClass()
{
    std::cout << "Confirming calculations in HIColumnDensity..." << std::endl;
    std::cout << "Normalization constant should be:         " << "something I've haven't yet decided to include or not" << std::endl;

    double bDLA, nl,nc,nu;
    returnPowerLawIndex(bDLA);
    returnColDensityLimits(nl,nu);
    returnColDensityBreak(nc);

    std::cout << "betaDLA: " << bDLA << std::endl
                << "Nl: " << nl << std::endl << "Nc: " << nc << std::endl << "Nu: " << nu << std::endl;

    return;
};





/******* AbsorberRedshiftDistributionLAF *********************************************/
AbsorberRedshiftDistributionLAF::AbsorberRedshiftDistributionLAF() 
: ALAF_(500.), z1LAF_(1.2), z2LAF_(4.7), gamma1LAF_(0.2), gamma2LAF_(2.7), gamma3LAF_(4.5)
{

//>>>>>>> b174eb658571fc61c81d79c7fa2a0db330092148
};

/******* AbsorberRedshiftDistributionLAF methods *************************************/
double AbsorberRedshiftDistributionLAF::returnRedshiftDist(double z)
{
    double dist = 0.0;
    if(z < z1LAF_) {
        dist = returnNormalization()*returnFirstPowerLaw(z);
    }
    else if(z < z2LAF_ && z >= z1LAF_) {
        dist = returnNormalization()*returnSecondPowerLaw(z);
    }
    else if(z >= z2LAF_) {
        dist = returnNormalization()*returnThirdPowerLaw(z);
    }
    else {
        std::cout << "Redshift not recognized!" << std::endl;
    }

    return dist;
};

double AbsorberRedshiftDistributionLAF::returnFirstPowerLaw(double z)
{
    double ratio = 0.0;
    ratio = (1+z)/(1+z1LAF_);
    return pow(ratio, gamma1LAF_);
};

double AbsorberRedshiftDistributionLAF::returnSecondPowerLaw(double z)
{
    double ratio = 0.0;
    ratio = (1+z)/(1+z1LAF_);
    return pow(ratio, gamma2LAF_);
};

double AbsorberRedshiftDistributionLAF::returnThirdPowerLaw(double z)
{
    double ratio2 = 0.0, ratio3 = 0.0;
    ratio2 = (1+z2LAF_)/(1+z1LAF_);
    ratio3 = (1+z)/(1+z2LAF_);
    return pow(ratio2,gamma2LAF_)*pow(ratio3,gamma3LAF_);
};


void AbsorberRedshiftDistributionLAF::returnPowerLawIndex(double &g1LAF, double &g2LAF, double &g3LAF)
{
    g1LAF = gamma1LAF_;
    g2LAF = gamma2LAF_;
    g3LAF = gamma3LAF_;
    return;
};

void AbsorberRedshiftDistributionLAF::returnRedshiftBreaks(double &z1LAF, double &z2LAF)
{
    z1LAF = z1LAF_;
    z2LAF = z2LAF_;
    return;
};

double AbsorberRedshiftDistributionLAF::returnNormalization()
{
    return ALAF_;
};

void AbsorberRedshiftDistributionLAF::testClass()
{
    double g1LAF, g2LAF, g3LAF, z1LAF, z2LAF, ALAF;
    returnPowerLawIndex(g1LAF,g2LAF,g3LAF);
    returnRedshiftBreaks(z1LAF,z2LAF);
    ALAF = returnNormalization();

    std::cout << "gamma1LAF " << g1LAF << std::endl
                << "gamma2LAF " << g2LAF << std::endl
                << "gamma3LAF " << g3LAF << std::endl
                << "z1LAF     " << z1LAF << std::endl
                << "z2LAF     " << z2LAF << std::endl
                << "ALAF      " << ALAF  << std::endl;
    std::cout << std::endl;

    return;
};



/******* AbsorberRedshiftDistributionDLA *********************************************/
AbsorberRedshiftDistributionDLA::AbsorberRedshiftDistributionDLA() 
: ADLA_(1.1), zDLA_(2.0), gamma1DLA_(1.0), gamma2DLA_(2.0)
{

//>>>>>>> b174eb658571fc61c81d79c7fa2a0db330092148
};

/******* AbsorberRedshiftDistributionDLA methods *************************************/
double AbsorberRedshiftDistributionDLA::returnRedshiftDist(double z)
{
    double dist = 0.0;
    if(z < zDLA_) {
        dist = returnNormalization()*returnFirstPowerLaw(z);
    }

    else if(z >= zDLA_) {
        dist = returnNormalization()*returnSecondPowerLaw(z);
    }
    else {
        std::cout << "Redshift not recognized!" << std::endl;
    }

    return dist;
};

double AbsorberRedshiftDistributionDLA::returnFirstPowerLaw(double z)
{
    double ratio = 0.0;
    ratio = (1+z)/(1+zDLA_);
    return pow(ratio, gamma1DLA_);
};

double AbsorberRedshiftDistributionDLA::returnSecondPowerLaw(double z)
{
    double ratio = 0.0;
    ratio = (1+z)/(1+zDLA_);
    return pow(ratio, gamma2DLA_);
};

void AbsorberRedshiftDistributionDLA::returnPowerLawIndex(double &g1DLA, double &g2DLA)
{
    g1DLA = gamma1DLA_;
    g2DLA = gamma2DLA_;
    return;
};

void AbsorberRedshiftDistributionDLA::returnRedshiftBreaks(double &zDLA)
{
    zDLA = zDLA_;
    return;
};

double AbsorberRedshiftDistributionDLA::returnNormalization()
{
    return ADLA_;
};

void AbsorberRedshiftDistributionDLA::testClass()
{
    double g1DLA, g2DLA, zDLA, ADLA;
    returnPowerLawIndex(g1DLA,g2DLA);
    returnRedshiftBreaks(zDLA);
    ADLA = returnNormalization();

    std::cout << "gamma1DLA " << g1DLA << std::endl
                << "gamma2DLA " << g2DLA << std::endl
                << "zDLA     " << zDLA << std::endl
                << "ADLA      " << ADLA  << std::endl;
    std::cout << std::endl;

    return;
};




/******* DopplerParDistribution ***************************************************/
DopplerParDistribution::DopplerParDistribution()
: bsigma_(23.0) 
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
ProbabilityDistAbsorbers::ProbabilityDistAbsorbers(RandomGeneratorInterface& rg, double zS, 
                                                   AbsorberRedshiftDistributionLAF& absorberZDistLAF,
                                                   AbsorberRedshiftDistributionDLA& absorberZDistDLA,
                                                   HIColumnDensityLAF& hiColumnDensityLAF,
                                                   HIColumnDensityDLA& hiColumnDensityDLA,
                                                   DopplerParDistribution& dopplerParDist)
: rg_(rg), zS(zS), absorberZDistLAF_(absorberZDistLAF), absorberZDistDLA_(absorberZDistDLA), hiColumnDensityLAF_(hiColumnDensityLAF), hiColumnDensityDLA_(hiColumnDensityDLA), dopplerParDist_(dopplerParDist)

{
    int nStep = 1000, nMaxLAF=3000, nMaxDLA= 50;

// preparing for MC LAF and DLA generation, redshift distribution part
	double dz= zS/nStep; 
	double z, yLAF, yDLA, pLAF, pDLA;
	zVector.push_back(0.0);
//	yLAFvector.push_back(0.0);
	yDLAvector.push_back(0.0);
	
	for (int i=1; i<nStep; i++)
	{
		z = zVector[i-1] + dz;
		zVector.push_back(z); 

//		cout << i << "	" << zVector[i] << endl; 

//		yLAF = yLAFvector[i-1] + 0.5*dz*absorberZDistLAF_(zVector[i-1]) + 0.5*dz*absorberZDistLAF_(zVector[i]); 
//		yLAFvector.push_back(yLAF);

		yDLA = yDLAvector[i-1] + 0.5*dz*absorberZDistDLA_(zVector[i-1]) + 0.5*dz*absorberZDistDLA_(zVector[i]); 
		yDLAvector.push_back(yDLA); 
	}

//	double NavgLAF= yLAFvector[nStep-1]; 
	double NavgDLA= yDLAvector[nStep-1];

//	cout << "LAF avg number " << NavgLAF << endl;
//	cout << "DLA avg number " << NavgDLA << endl;

//	transform(yLAFvector.begin(), yLAFvector.end(), yLAFvector.begin(), bind1st(divides<double>(), NavgLAF)); //divides the whole vector by const
//	transform(yDLAvector.begin(), yDLAvector.end(), yDLAvector.begin(), bind1st(divides<double>(), NavgDLA));  //using transform, bind1st, divides

	for(int i=0; i<nStep; i++)
	{
//		yLAFvector[i] = yLAFvector[i]/NavgLAF;
		yDLAvector[i] = yDLAvector[i]/NavgDLA;
//		cout << i << "	" << yLAFvector[i] << "	" << yDLAvector[i] << endl; 
	}

/*	for (int i=0; i<nMaxLAF; i++)
	{
		if(i==0)
		{
			pLAF = exp(-NavgLAF);
			cppLAFvector.push_back(pLAF);

		}
		
		else
		{
			pLAF = (pLAF*NavgLAF)/i;
			cppLAFvector.push_back(cppLAFvector[i-1]+pLAF);

		}

		cout << i << "	" << cppLAFvector[i] << endl; 
	} 
*/

	for (int i=0; i<nMaxDLA; i++)
	{
		if(i==0)
		{
			pDLA = exp(-NavgDLA);
			cppDLAvector.push_back(pDLA);
		}
		
		else
		{
		
			pDLA = (pDLA*NavgDLA)/i;
			cppDLAvector.push_back(cppDLAvector[i-1]+pDLA);
		}


//		cout << i << "	" << cppDLAvector[i] << endl; 
	}
			
	
	//doppler distribution part
    bmin_ = 0.0;
    bmax_ = 200.;
    setDopplerDistribution(nStep);


	//column density distribution part 
    string infileLAF= "LookupTableLAF.tbl";  	//Lookup table for LAF column density
    string infileDLA= "LookupTableDLA.tbl";	//Lookup table for DLA column density 

    ifstream inputLAF; 
    ifstream inputDLA; 

    inputLAF.open(infileLAF.c_str(), ifstream::in);
    inputDLA.open(infileDLA.c_str(), ifstream::in);

    double randomLAF, columnDensityLAF; 
    double randomDLA, columnDensityDLA; 

    while(inputLAF >> randomLAF >> columnDensityLAF)
    {
	randomVectorLAF.push_back(randomLAF);
	columndensityVectorLAF.push_back(columnDensityLAF);
    }

    while(inputDLA >> randomDLA >> columnDensityDLA)
    {
	randomVectorDLA.push_back(randomDLA);
	columndensityVectorDLA.push_back(columnDensityDLA);
    }

    inputLAF.close();
    inputDLA.close(); 

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

double ProbabilityDistAbsorbers::drawDeltaZLAF(double zLast) 
{
    double fz = absorberZDistLAF_(zLast);
    double rn = rg_.Flat01();

    // p(\DeltaZ;z) = f(z)*exp( -f(z)*\DeltaZ )

    double diff = -log(rn); 
    double deltaZ = diff/fz;

//    cout << fz << endl; 	 
 
    return deltaZ;
};

double ProbabilityDistAbsorbers::drawDeltaZDLA(double zLast)  
{
    double fz = absorberZDistDLA_(zLast);
    double rn = rg_.Flat01();

    // p(\DeltaZ;z) = f(z)*exp( -f(z)*\DeltaZ )

    double diff = -log(rn); 
    double deltaZ= diff/fz; 

    return deltaZ;
};


double ProbabilityDistAbsorbers::drawHIColumnDensityLAF()
{
    double nhi = 0.0;
    double ran = rg_.Flat01();

    vector<double>::iterator low;		
    low= lower_bound (randomVectorLAF.begin(), randomVectorLAF.end(), ran);   //search through lookup table

    nhi= columndensityVectorLAF[low-randomVectorLAF.begin()]; 
    return nhi;

}


double ProbabilityDistAbsorbers::drawHIColumnDensityDLA()
{
    double nhi=0.0;
    double ran = rg_.Flat01();

    vector<double>::iterator low;
    low= lower_bound (randomVectorDLA.begin(), randomVectorDLA.end(), ran);	//search through lookup table

    nhi= columndensityVectorDLA[low-randomVectorDLA.begin()]; 
    return nhi;

}



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




void ProbabilityDistAbsorbers::simulateLineOfSightLAF(double zStart, double zMax,
                    vector<double>& redshiftsLAF, vector<double>& dopplerParsLAF,
                    vector<double>& columnDensitiesLAF, string outfile)
{
    double zCurrent = zStart;
    int iAbsorber = 0;

    while(zCurrent < zMax) {
        double zNext, bdopp, NHI;
        simulateAbsorberLAF(zCurrent,zNext,bdopp,NHI);

        redshiftsLAF.push_back(zNext);
        dopplerParsLAF.push_back(bdopp);
        columnDensitiesLAF.push_back(NHI);

        zCurrent = zNext;
        iAbsorber++;
    }

      cout << iAbsorber << endl;  //counting number of LAF absorbers per line of sight 

    return;
};

void ProbabilityDistAbsorbers::simulateAbsorberLAF(double zCurrent, double& zNext,
                    double& bdopp, double& NHI)
{
    zNext = zCurrent + drawDeltaZLAF(zCurrent);
    bdopp = drawDoppler();
    NHI = drawHIColumnDensityLAF();
    return;
};





void ProbabilityDistAbsorbers::simulateLineOfSightDLA(double zStart, double zMax,
                    vector<double>& redshiftsDLA, vector<double>& dopplerParsDLA,
                    vector<double>& columnDensitiesDLA, string outfile)
{
//    double zCurrent = zStart;
//    int iAbsorber = 0;
	int nMax= 50, NN= 1000;
	int numDLA, k; 
    double zNext, bdopp, NHI;

    double ran = rg_.Flat01();

/*	for(int k=0; k<nMax; k++)
	{
		if(ran <= cppDLAvector[k]) 
			numDLA = k; 
	}
*/

    vector<double>::iterator lowNum, lowZ;
    lowNum= lower_bound(cppDLAvector.begin(), cppDLAvector.end(), ran);	//search through vector 
    numDLA= (lowNum-cppDLAvector.begin()); 

//	cout << "Number of DLAs in LoS: " << numDLA << endl; 

	for(int j=1; j<=numDLA; j++)
	{
		ran = rg_.Flat01();

		lowZ= lower_bound(yDLAvector.begin(), yDLAvector.end(), ran);

		k= (lowZ-yDLAvector.begin());

		zNext = (zVector[k-1]*(ran-yDLAvector[k-1]) + zVector[k]*(yDLAvector[k]-ran))/(yDLAvector[k]-yDLAvector[k-1]); 
		redshiftsDLA.push_back(zNext); 

/*		for(int k=1; k<NN; k++)
		{
			if(ran < yDLAvector[k])
			{
				zNext = (zVector[k-1]*(ran-yDLAvector[k-1]) + zVector[k]*(yDLAvector[k]-ran))/(yDLAvector[k]-yDLAvector[k-1]); 
				redshiftsDLA.push_back(zNext); 
			}

		}
*/
//	   cout << "DLA redshift is " << zNext << endl; 

	   simulateAbsorberDLA(bdopp,NHI);
       dopplerParsDLA.push_back(bdopp);
       columnDensitiesDLA.push_back(NHI);
	}


/*    while(zCurrent < zMax) {
        double zNext, bdopp, NHI;
        simulateAbsorberLAF(zCurrent,zNext,bdopp,NHI);

        redshiftsLAF.push_back(zNext);
        dopplerParsLAF.push_back(bdopp);
        columnDensitiesLAF.push_back(NHI);

        zCurrent = zNext;
        iAbsorber++;
    }
*/

      cout << numDLA << endl;  //counting number of DLA absorbers per line of sight 

    return;
};

void ProbabilityDistAbsorbers::simulateAbsorberDLA(double& bdopp, double& NHI)
{
//    zNext = zCurrent + drawDeltaZDLA(zCurrent);
    bdopp = drawDoppler();
    NHI = drawHIColumnDensityDLA();

    return;
};





/******* VoigtProfile *************************************************************/
VoigtProfile::VoigtProfile(double dopplerPar, int nLine)
: dopplerPar_(dopplerPar), nLine_(nLine)
{
  //    double nGammaMax = gammaSeries_.size() + 1;
    double nGammaMax = nGammaMax_ + 1;
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

    if(freq >= freqLymanLimitInvSec_) 
	{
        double ratio = freqLymanLimitInvSec_/freq;
//        sigmaLC = sigmaLymanLimitCM2_*ratio*ratio*ratio;
		sigmaLC= sigmaLymanLimitCM2_*(R*pow(ratio,s)+(1.0-R)*pow(ratio,s+1.0)); 
    } 
	else if(freq > returnFrequencyLymanSeries(35))
	{
		sigmaLC = sigmaLymanLimitCM2_;
	}
	else 
	{
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
    nLineMaxMadau_ = Avals_.size() + 1;

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

