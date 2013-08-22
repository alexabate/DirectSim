#include "igm.h"

/******* AtomicCalcs methods **************************************************/

double AtomicCalcs::returnDampingParameter(int nLine, double dopplerPar)
{   
    
    double lambdaN = returnWavelengthLymanSeries(nLine);
    
    int iSeries = nLine - 2;
    double gammaLine = gammaSeries_[iSeries];
    
    double a = (lambdaN*lambdaN*gammaLine)/
                    (4*PI*SPEED_OF_LIGHT_MS*returnDopplerWidthWL(nLine, dopplerPar));
                    
    return a; 

};

// To set the distributions and the constants
void AtomicCalcs::setConstants()
{

    freqLymanLimInvSecs_ = SPEED_OF_LIGHT_MS/WAVE_LYMANLIM_METERS; // in s^-1
    sigLymanLimCM2_  = 6.30e-18; // in cm^2
    nLymanAlpha_ = 2;
    
   
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
    
};

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
    gammaSeries_.push_back(1.91e6); // 13
    gammaSeries_.push_back(1.528e6);// 14
    gammaSeries_.push_back(1.242e6);// 15
    gammaSeries_.push_back(1.024e6);// 16
    
    gammaSeries_.push_back(8.532e5);// 17
    gammaSeries_.push_back(7.185e5);// 18
    gammaSeries_.push_back(6.109e5);// 19 
    gammaSeries_.push_back(5.235e5);// 20
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
   
};

/******* HIColumnDensity ******************************************************/

HIColumnDensity::HIColumnDensity(double beta1, double beta2, double Nc, 
                                                double Nl, double Nu, int Nstep)
: beta1_(beta1) , beta2_(beta2) , Nc_(Nc) , Nl_(Nl) , Nu_(Nu)
{

    Bnorm_ = 1; // temporarily set for now 
    Bnorm_ = setBnorm(Nstep);
    cout <<"     Normalization value B = "<< Bnorm_ <<endl;
    double intVal = checkIntegration(Nstep);
    double numIntVal = numIntegratePowerLaw(Nstep);
    cout <<"     After normalization column density function numercially integrates to "<< intVal <<endl;

};

/******* HIColumnDensity methods **********************************************/

double HIColumnDensity::setBnorm(int Nstep)
{

    double beta1part = integratePowerLaw(Nc_,beta1_) - integratePowerLaw(Nl_,beta1_);                            
    double beta2part = integratePowerLaw(Nu_,beta2_) - integratePowerLaw(Nc_,beta2_);
                                                
    double Bnorm = 1./(beta1part + beta2part);
    return Bnorm;


};


double HIColumnDensity::returnColDensityDist(double Nh)
{
   
    double dist;
    if ( Nh < Nc_ ) {
        dist = returnPowerLawNormalization()*returnFirstPowerLaw(Nh);
        }
    else if ( Nh >= Nc_ ) {
        dist = returnPowerLawNormalization()*returnSecondPowerLaw(Nh);
        }
    else {
        string emsg = "ERROR! column density value not recognized";
        throw emsg;
        }
      
    return dist;
     
};

void HIColumnDensity::writeToFile(string outfile, double dLog, int nStep)
{

    double logNl = log(Nl_);

    ifstream inp;
	ofstream outp;
	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(ios::failbit);
		cout << "     Writing HI column density distribution to file ...";
		cout << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (int i=0; i<nStep; i++) {
	    
			double logNh=logNl+i*dLog;
			double Nh=exp(logNh);
			outp << Nh <<"  "<< returnColDensityDist(Nh) << endl;
            }
            
		outp.close();
	    }
	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;

};

double HIColumnDensity::checkIntegration(int Nstep)
{

    double numIntVal = returnPowerLawNormalization()*numIntegratePowerLaw(Nstep);
    return numIntVal;
};

double HIColumnDensity::numIntegratePowerLaw(int Nstep)
{
    // do integration in log steps
    double logNl = log(Nl_);
    double logNu = log(Nu_);
    double dlogNh = (logNu - logNl )/(Nstep - 1);

    double sum=0;
    for ( int i=0; i<Nstep; i++ ) {
        double logNh = logNl + i*dlogNh;
        double Nh = exp(logNh);
 
        if ( Nh < Nc_ )
            sum += returnFirstPowerLaw(Nh)*Nh*dlogNh;
        else if ( Nh >= Nc_ )
            sum += returnSecondPowerLaw(Nh)*Nh*dlogNh;
    
        }
        
    return sum;
};

/******* AbsorberRedshiftDistribution *****************************************/

AbsorberRedshiftDistribution::AbsorberRedshiftDistribution(double gamma1,
double gamma2, double gamma3, double z1, double z2, double A)
: gamma1_(gamma1) , gamma2_(gamma2) , gamma3_(gamma3) , z1_(z1) , z2_(z2) ,A_(A)
{
             
};

/******* AbsorberRedshiftDistribution methods *********************************/

double AbsorberRedshiftDistribution::returnRedshiftDist(double z)
{

    double dist;
    if ( z <= z1_ ) {
        dist = A_*returnFirstPowerLaw(z);
        }
    else if ( z <= z2_ && z > z1_ ) {
        dist = A_*returnSecondPowerLaw(z);
        }
    else if ( z > z2_ )
        dist = A_*returnThirdPowerLaw(z);
    else
        throw ParmError("ERROR! redshift value not recognized");
        
    return dist;


};

void AbsorberRedshiftDistribution::writeToFile(string outfile, double zmin, 
                                                           double dz, int nStep)
{

    ifstream inp;
	ofstream outp;

    inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(ios::failbit);
		cout << "     Writing to absorber redshift distribution to file ...";
		cout << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (int i=0; i<nStep; i++) {		
			double zv=zmin+i*dz;
			outp << zv <<"  "<< returnRedshiftDist(zv) << endl;
            }
            
		outp.close();
	    }
	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;

};



/******* DopplerParDistribution methods ***************************************/

double DopplerParDistribution::returnDopplerDist(double b)
{

    double bsig4=bsigma_*bsigma_*bsigma_*bsigma_;
    double b4=b*b*b*b;
    double bratio=bsig4/b4;
    double bratio2=bratio/b;
    
    return 4*bratio2*exp(-bratio);

};

void DopplerParDistribution::writeToFile(string outfile, double bmin, 
                                                           double db, int nStep)
{

    ifstream inp;
	ofstream outp;

    inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(ios::failbit);
		cout << "     Writing Doppler parameter distribution to file ...";
		cout << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (int i=0; i<nStep; i++) {		
			double bv=bmin+i*db;
			outp << bv <<"  "<< returnDopplerDist(bv) << endl;
            }
            
		outp.close();
	    }
	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;

};

/******* ProbabilityDistAbsorbers  ********************************************/

ProbabilityDistAbsorbers::ProbabilityDistAbsorbers(RandomGeneratorInterface& rg,
                            AbsorberRedshiftDistribution& absorberZDist,
                            HIColumnDensity& hiColumnDensity,
                            DopplerParDistribution& dopplerParDist
        )
: rg_(rg) , absorberZDist_(absorberZDist) , hiColumnDensity_(hiColumnDensity) ,
                                                dopplerParDist_(dopplerParDist)
{ 

    int nStep = 1000;
    setNHiDistribution(nStep);
    bmin_ = 0, bmax_ = 200;
    setDopplerDistribution(nStep);
    



};

/******* ProbabilityDistAbsorbers methods *************************************/

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
};

void ProbabilityDistAbsorbers::simulateAbsorber(double zCurrent, 
                    double& redshift, double& dopplerPar, double& columnDensity)
{
    redshift = zCurrent + drawDeltaZ(zCurrent);
	dopplerPar = drawDoppler();
	columnDensity = drawHIColumnDensity();

};

/*double ProbabilityDistAbsorbers::simulateAbsorberRedshift(double zSource)
{
    double redshift = findzInoue(zSource);
    return redshift;
};*/


// To draw from the distributions
double ProbabilityDistAbsorbers::drawDeltaZ(double zLast)
{
    
    // returns number of absorbers expected at z
    double fz=absorberZDist_(zLast);

    // Probability the next absorber on a line of sight is at zLast + deltaZ:
    // P(deltaZ | zLast) = fz*exp(-fz*deltaZ)
    
    // Normalize this distribution so it is properly cumulative (i.e. tends
    // to Prob = 1 ->
    // P(deltaZ | zLast) = exp(-fz*deltaZ)
    
    // Use inverse transformation method to generate correct distribution
    double rn = rg_.Flat01();
    
    double deltaZ=-1/fz*log(rn); // this is inverse of above function
    
    return deltaZ;

};

// Code below is from IGMtransmission code - works same as drawDeltaZ above
double ProbabilityDistAbsorbers::simulateAbsorberRedshift(double z0) 
{
		double random = rg_.Flat01();
		double gamma1, gamma2, gamma3;
		absorberZDist_.returnPowerLaws(gamma1, gamma2, gamma3);
        double z1, z2;
        absorberZDist_.returnzAtBreaks(z1, z2);

		double z = 0;
		if (z0 <= z1) {
			double normalisation1 = ((1 + z1) / (gamma1 + 1))
					* pow(((1 + z0) / (z1 + 1)), (gamma1 + 1))
					- ((1 + z1) / (gamma1 + 1))
					* pow((1 / (z1 + 1)), (gamma1 + 1));
			z = pow(
					(random * normalisation1 * ((gamma1 + 1) / (1 + z1)) + 
					pow((1 / (1 + z1)), (gamma1 + 1))),
					(1 / (gamma1 + 1)))
					* (1 + z1) - 1;
		} else if (z1 < z0 && z0 <= z2) {
			double normalisation2 = ((1 + z1) / (1 + gamma1))
					* (1 - pow((1 / (1 + z1)), (gamma1 + 1)))
					+ ((1 + z1) / (1 + gamma2))
					* (pow(((1 + z0) / (1 + z1)), (gamma2 + 1)) - 1);
			double Fz1 = ((1 + z1) / (gamma1 + 1))
					* (1 - pow((1 / (1 + z1)), (gamma1 + 1)))
					/ normalisation2;
			if (0 < random && random <= Fz1) {
				z = pow((random * normalisation2
						* ((gamma1 + 1) / (1 + z1)) + pow((1 / (1 + z1)),
						(gamma1 + 1))), (1 / (gamma1 + 1)))
						* (1 + z1) - 1;
			} else if (Fz1 < random && random <= 1) {
				z = pow(((random * normalisation2 - Fz1 * normalisation2)
						* ((gamma2 + 1) / (1 + z1)) + 1), (1 / (gamma2 + 1)))
						* (1 + z1) - 1;
			} else {
				throw ParmError("ERROR!");

			}
		} else if (z0 > z2) {
			double normalisation3 = ((1 + z1) / (gamma1 + 1))
					* (1 - pow((1 / (1 + z1)), (gamma1 + 1)))
					+ ((1 + z1) / (gamma2 + 1))
					* (pow(((1 + z2) / (1 + z1)), (gamma2 + 1)) - 1)
					+ ((1 + z2) / (gamma3 + 1))
					* (pow(((1 + z2) / (1 + z1)), gamma2))
					* (pow(((1 + z0) / (1 + z2)), (gamma3 + 1)) - 1);
			double Fz1 = (((1 + z1) / (gamma1 + 1)) * (1 - pow(
					(1 / (1 + z1)), (gamma1 + 1))))
					/ normalisation3;
			double Fz2 = Fz1
					+ (((1 + z1) / (gamma2 + 1)) * (pow(
							((1 + z2) / (1 + z1)), (gamma2 + 1)) - 1))
					/ normalisation3;
			// System.out.println("Fz1: " + Fz1 + " Fz2: " + Fz2);
			if (0 < random && random <= Fz1) {
				z = pow((random * normalisation3
						* ((gamma1 + 1) / (1 + z1)) + pow((1 / (1 + z1)),
						(gamma1 + 1))), (1 / (gamma1 + 1)))
						* (1 + z1) - 1;
			} else if (Fz1 < random && random <= Fz2) {
				z = pow(
						((normalisation3 * random - ((1 + z1) / (gamma1 + 1))
								* (1 - pow((1 / (1 + z1)), (gamma1 + 1))))
								* ((gamma2 + 1) / (1 + z1)) + 1),
						(1 / (gamma2 + 1)))
						* (1 + z1) - 1;
			} else if (random > Fz2) {

				z = pow(
								(((gamma3 + 1) / (1 + z2))
										* pow(((1 + z1) / (1 + z2)),
												gamma2)
										* (random
												* normalisation3
												- ((1 + z1) / (gamma1 + 1))
												* (1 - pow((1 / (1 + z1)),
														(gamma1 + 1))) - ((1 + z1) / (gamma2 + 1))
												* (pow(
														((1 + z2) / (1 + z1)),
														(gamma2 + 1)) - 1)) + 1),
								(1 / (gamma3 + 1)))
						* (1 + z2) - 1;
			} else {
				throw ParmError("ERROR!");
			}
		} else {
			throw ParmError("ERROR!");
		}
		return z;
	};


double ProbabilityDistAbsorbers::drawHIColumnDensity()
{

    //cout <<"     log10(Nl) = "<<log10Nl_<<",  log10(Nu) = "<<log10Nu_<< endl;
    //cout <<"     log10(gin) = "<<log10gmin_<<",  log10(gmax) = "<<log10gmax_<< endl;
    double u1,u2;
    while (true) {
        u1 = log10Nl_ + (log10Nu_ - log10Nl_)*rg_.Flat01();
        u2 = log10gmin_ + (log10gmax_ - log10gmin_)*rg_.Flat01();
       // u2 = log10gmin_*rg_.Flat01();// between 0 and log10gmax
        
        double log10gval = colDensityFunc_(u1);
        //int iElement=findClosestElement(log10NHIvals_,u1);
   
        if ( u2 <= log10gval)//log10gvals_[iElement])
            break;
        }
    return pow(10.,u1);

};


double ProbabilityDistAbsorbers::drawDoppler()
{

    double u1,u2;
    while (true)  {
        u1 = bmin_ + (bmax_ - bmin_)*rg_.Flat01();
        u2 = hmin_ + (hmax_ - hmin_)*rg_.Flat01();
    
        int iElement=findClosestElement(bvals_,u1);
   
        if ( u2 <= hvals_[iElement])
            break;
        }
    return u1;

};


void ProbabilityDistAbsorbers::setNHiDistribution(int nStep)
{

    log10NHIvals_.clear();
    log10gvals_.clear();
    
    double Nl,Nu;
    hiColumnDensity_.returnLowerUpperColDensValues(Nl,Nu);
    log10Nl_ = log10(Nl);
    log10Nu_ = log10(Nu);
    
    double dlogN=(log10Nu_ - log10Nl_)/(nStep-1);
    
    for (int i=0; i<nStep; i++) {
        double log10NHI = log10Nl_ + i*dlogN;
        double NHI = pow(10.,log10NHI);
        log10NHIvals_.push_back(log10NHI);
        double g=hiColumnDensity_.returnColDensityDist(NHI);
        log10gvals_.push_back(log10(g));
        }
        
    // Want to subtract the minimum values of log10(g) from all the elements in
    // the log10gvals_ vector because log10gvals is negative
    double log10gminTmp = findMinimum(log10gvals_);
    //for (int i=0; i<nStep; i++)
    //    log10gvals_[i] -= log10gminTmp;

        
    // check the vectors log10gvals and log10NHIvals_
    // check the maximum values here
    log10gmin_=findMinimum(log10gvals_); // new minimum must be zero
    log10gmax_=findMaximum(log10gvals_);
    cout << " log10gmax = " << log10gmax_ << " log10gmin = " << log10gmin_ << endl;
    
    // Make interpolation
    colDensityFunc_.DefinePoints(log10NHIvals_,log10gvals_);

    double eps = 1e-4;
    //if (abs(log10gmin_)>eps)
    //    throw ParmError("ERROR! Column density distribution minimum not zero!");
    
};


void ProbabilityDistAbsorbers::setDopplerDistribution(int nStep)
{

    bvals_.clear();
    hvals_.clear();
    
    double db = (bmax_-bmin_)/(nStep-1);
    
    for (int i=0; i<nStep; i++) {
        double bv = bmin_ + i*db;

        bvals_.push_back(bv);
        double hv=dopplerParDist_(bv);
        hvals_.push_back(hv);
        }

    hmin_=findMinimum(hvals_);
    hmax_=findMaximum(hvals_);
};





// Write the line of sight absorbers to a file 
void ProbabilityDistAbsorbers::writeToFile(string outfile, 
    vector<double>& redshifts, vector<double>& dopplerPars, 
                                                vector<double>& columnDensity)
{

    int nAbsorbers = columnDensity.size();

    ifstream inp;
	ofstream outp;

	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(ios::failbit);
		cout << "     Writing line of sight to file ...";
		cout << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (long i=0; i<nAbsorbers; i++) {
	        
	        double zv = redshifts[i];
	        double bv = dopplerPars[i];
	        double cv = columnDensity[i];
	        
	        outp << zv <<"  "<< bv <<"  "<< cv <<endl;
	        
	        }
	    outp.close();
	    }
	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;



};


// Mostly for debugging below here
void ProbabilityDistAbsorbers::writeZDistribution(string outfile, 
                                            double zCurrent, long nTrial)
{

    ifstream inp;
	ofstream outp;

	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(ios::failbit);
		cout << "     Writing delta z list to file ...";
		cout << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (long i=0; i<nTrial; i++) {
	        
	        double deltaZ=drawDeltaZ(zCurrent);
	        
	        outp << deltaZ <<endl;
	        
	        }
	    outp.close();
	    }
	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;


};


void ProbabilityDistAbsorbers::writeDopplerDistribution(string outfile, 
                                                                   long nTrial)
{

    ifstream inp;
	ofstream outp;

	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(ios::failbit);
		cout << "     Writing doppler list to file ...";
		cout << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (long i=0; i<nTrial; i++) {
	        
	        double doppler=drawDoppler();
	        
	        outp << doppler <<endl;
	        
	        }
	    outp.close();
	    }
	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;


};

void ProbabilityDistAbsorbers::writeNHiDistribution(string outfile, long nTrial)
{


    ifstream inp;
	ofstream outp;

	inp.open(outfile.c_str(), ifstream::in);
	inp.close();
	if(inp.fail()) {
		
		inp.clear(ios::failbit);
		cout << "     Writing column density list to file ...";
		cout << outfile.c_str() << endl;
		outp.open(outfile.c_str(), ofstream::out);
		
	    for (long i=0; i<nTrial; i++) {
	        
	        double NHI=drawHIColumnDensity();
	        
	        outp << NHI <<endl;
	        
	        }
	    outp.close();
	    }
	else
		cout << "Error...file " << outfile.c_str() << " exists" << endl;


};




/******* OpticalDepth methods *************************************************/

// All frequencies are REST FRAME (of the absorber) in these calculations
double OpticalDepth::returnRestFrameOpticalDepth(double freq, double bAbsorber, double nHI)
{

    double sigmaLC = 0., sumSigmai = 0., opticalDepthRestFrame = 0.;
    
    // cross-section for the Lyman continuum 
    sigmaLC = returnLymanContinuumCrossSection(freq);
    
    // cross-section for the Lyman series
    sumSigmai = returnLymanSeriesCrossSection(freq, bAbsorber);
    
    // total optical depth
    if (isAll_)
        opticalDepthRestFrame = nHI*(sigmaLC + sumSigmai);
    else if (isOnlyLymanC_)
        opticalDepthRestFrame = nHI*sigmaLC;
    else if (isOnlyLymanS_)
        opticalDepthRestFrame = nHI*sumSigmai;
    else
        throw ParmError("ERROR! Unknown optical depth contributions");
    
    return opticalDepthRestFrame;

};

// All frequencies are REST FRAME (of the absorber) in these calculations
// Continuum
double OpticalDepth::returnLymanContinuumCrossSection(double freq)
{
    
    double sigLC;
    
    if ( freq >= freqLymanLimInvSecs_) {
        double freqRatio = freqLymanLimInvSecs_/freq;
        sigLC = sigLymanLimCM2_*freqRatio*freqRatio*freqRatio;
        }
    else
        sigLC = 0.;
        
    return sigLC; // in cm^2
};

// All frequencies are REST FRAME (of the absorber) in these calculations
// Series
double OpticalDepth::returnLymanSeriesCrossSection(double freq, double bAbsorber)
{

    double sumSigmai=0.;
    for (int n=nLymanAlpha_; n<=nLineMax_; n++)
        sumSigmai += returnLymanLineCrossSection(n, freq, bAbsorber);
    return sumSigmai; // in cm^2
    
};

// All frequencies are REST FRAME (of the absorber) in these calculations
double OpticalDepth::returnLymanLineCrossSection(int nLine, double freq, double bAbsorber)
{

    // Constant: in cgs units
    double constants = (sqrt(PI)*ELECTRON_CHARGE_STATC*ELECTRON_CHARGE_STATC)/
                                        (ELECTRON_MASS_G*SPEED_OF_LIGHT_CMS);// in cm^2/s
                            
    // Doppler width of line i
    double dopplerWidth = returnDopplerWidthFreq(nLine, bAbsorber); // in s^-1
    
    // Oscillator strength of line i
    double fi = returnOscillatorStrength(nLine); // unitless
    
    // Line profile function of line i
    double phii = returnLineProf(nLine, freq, bAbsorber); // unitless
    
    // Lyman series cross-section
    double sigmai = constants*(fi/dopplerWidth)*phii; // units of cm^2
    
    //cout << constants*(fi/dopplerWidth) <<endl;
    //cout << "Nline="<<nLine<<"  "<<constants <<"  "<< fi <<"  "<< dopplerWidth <<"  "<<phii<<endl;
                       
    return sigmai;

};

// All frequencies are REST FRAME (of the absorber) in these calculations
double OpticalDepth::returnLineProf(int nLine,double freq, double bAbsorber)
{
    VoigtProfile voigtProf(bAbsorber, nLine);
    double lambda = SPEED_OF_LIGHT_MS/freq;
    
    double vprof = voigtProf(lambda);
    return vprof;

};

/*// To set the distributions and the constants
void OpticalDepth::setConstants()
{

    freqLymanLimInvSecs_ = SPEED_OF_LIGHT_MS/WAVE_LYMANLIM_METERS; // in s^-1
    sigLymanLimCM2_  = 6.30e-18; // in cm^2
    nLymanAlpha_ = 2;
    nLineMaxMax_ =  31; // Maximum line series possible to go up to
    nLineMax_ = nLineMaxMax_; // Maximum line series to actually go up to
   
    // Absorption oscillator strength of ith Lyman line
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
    
};*/

/******* LineOfSightTrans methods *********************************************/

double LineOfSightTrans::returnOpticaldepth(double lambda, double zSource)
{

    // Set contributions of absorbers
    if (isAll_)
        setLymanAll();
    else if (isOnlyLymanC_)
        setLymanContinuumOnly();
    else if (isOnlyLymanS_)
        setLymanSeriesOnly();
    else
        throw ParmError("ERROR! Unknown optical depth contributions");

    // find most distant absorber 
    int iDistant = findClosestElement(redshifts_,zSource);
    
    // loop over all the absorbers, summing up the optical depth of each
    //bool isOpticalDepth=true;
    double tau=0.;
    for (int i=0; i<=iDistant; i++) {
        
        double zAbsorber = redshifts_[i];
        double nHI = columnDensities_[i];
        double bAbsorber = dopplerPars_[i];
        
        double tauAbsorber = returnObserverFrameOpticalDepth(lambda, zAbsorber, nHI, bAbsorber);

        tau += tauAbsorber;
        }
        
    return tau;

};

/******* VoigtProfile methods *************************************************/

/*double VoigtProfile::returnDampingParameter()
{   
    
    double lambdaN = returnWavelengthLymanSeries(nLine_);
    
    int iSeries = nLine_-2;
    double gammaLine = gammaSeries_[iSeries];
    
    double a = (lambdaN*lambdaN*gammaLine)/
                    (4*PI*SPEED_OF_LIGHT_MS*returnDopplerWidth(nLine_,dopplerPar_));
                    
    return a; 

};*/


double VoigtProfile::kFunction(double x)
{

    double term1 = (4.*x*x + 3.)*(x*x + 1.)*exp(-x*x);
    double term2 = (1./(x*x))*(2.*x*x + 3.)*sinh(x*x);
    double k = (1./(2*x*x))*(term1 - term2);
    
    // check if inf
    int isInf = my_isinf(k);
    
    if (isInf != 0)// if k IS infinite set it to zero
        k = 0.;
    
    return k;

};

/*void VoigtProfile::setGamma()
{

    // Found in Table 2 of Morton 2003
    gamma_.push_back(6.265e8);// 2
    gamma_.push_back(1.897e8);// 3
    gamma_.push_back(8.127e7);// 4
    gamma_.push_back(4.204e7);// 5
    gamma_.push_back(2.450e7);// 6
    
    // Backed out of a file containg values of the damping parameter a assuming 
    // b=36km/s from Tepper-Garcia.  See Fig 1 in Tepper-Garcia 2006
    gamma_.push_back(1.236e7);// 7
    gamma_.push_back(8.252e6);// 8
    gamma_.push_back(5.781e6);// 9
    gamma_.push_back(4.209e6);// 10
    gamma_.push_back(3.159e6);// 11
    
    gamma_.push_back(2.431e6);// 12
    gamma_.push_back(1.91e6); // 13
    gamma_.push_back(1.528e6);// 14
    gamma_.push_back(1.242e6);// 15
    gamma_.push_back(1.024e6);// 16
    
    gamma_.push_back(8.532e5);// 17
    gamma_.push_back(7.185e5);// 18
    gamma_.push_back(6.109e5);// 19 
    gamma_.push_back(5.235e5);// 20
    gamma_.push_back(4.522e5);// 21 
    
    gamma_.push_back(3.932e5);// 22
    gamma_.push_back(3.442e5);// 23
    gamma_.push_back(3.029e5);// 24
    gamma_.push_back(2.678e5);// 25
    gamma_.push_back(2.381e5);// 26
    
    gamma_.push_back(2.127e5);// 27
    gamma_.push_back(1.906e5);// 28
    gamma_.push_back(1.716e5);// 29
    gamma_.push_back(1.549e5);// 30
    gamma_.push_back(1.405e5);// 31

    gamma_.push_back(1.277e5);// 32
    
    // size of gamma =  31
   
};*/

/******* Madau methods ********************************************************/

double Madau::returnObserverFrameOpticalDepth(double lambdaObs, double zSource)
{
	// Put the wavelength into the emission frame
	double lambdaEm = lambdaObs / (1 + zSource);
	
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

		int debug = 0;
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

