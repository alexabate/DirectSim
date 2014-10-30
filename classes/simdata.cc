#include "simdata.h"

/******* PhotometryCalcs methods **********************************************/


// K-correction calculations 

double PhotometryCalcs::Kcorr(double z, ClassFunc1D& sed, Filter& filterX, 
                                                        Filter& restFrameFilter)
{
// rest-frame and observed filters need NOT be the SAME
// compute general K correction for a galaxy with:
// -redshift z
// -spectral energy distribution sed (sed might have extinction applied already)
// -from REST FRAME bandpass Y
// -to OBSERVED bandpass X

	double kxy;
		
	// returns sedred(lambda/(1+z))*filtX(lambda)*lambda;
	SEDzFilterProd szfx(sed, filterX, z);
	// returns filtY(lambda)/lambda  
	FilterProd  fprody(restFrameFilter);          
	// returns sed(lambda)*filtY(lambda)*lambda;    
	SEDzFilterProd szfB0(sed, restFrameFilter, 0.);
	 // returns filtX(lambda)/lambda
	FilterProd  fprodx(filterX);             
	
    // To integrate the above
	FilterIntegrator trpzx(szfx, lmin_, lmax_);	
	FilterIntegrator trpzfy(fprody, lmin_, lmax_);	
	FilterIntegrator trpzy(szfB0, lmin_, lmax_);
	FilterIntegrator trpzfx(fprodx, lmin_, lmax_);	
		
	kxy=-2.5*log10(pow((1+z),-1)*(trpzx.Value()/trpzfx.Value())*
	                                    (trpzfy.Value()/trpzy.Value()));
	if (z<0.0001)
	    cout << z <<"  "<< kxy <<endl;
	return kxy;
};


double PhotometryCalcs::Kcorr1Filter(double z, ClassFunc1D& sed, Filter& filterX)
{
// 1 FILTER: rest-frame and observed filters are the SAME
// compute general K correction for a galaxy with:
// -redshift z
// -spectral energy distribution sed (sed might have extinction applied already)

	double kxx;
	
	// returns sedred(lambda/(1+z))*filtX(lambda)*lambda;
	SEDzFilterProd szfx(sed, filterX, z);  
	// returns sedred(lambda)*filtX(lambda)*lambda;
	SEDzFilterProd szfx0(sed, filterX, 0);  

    // integrate the above
	FilterIntegrator trpzx(szfx, lmin_, lmax_);	
	FilterIntegrator trpzx0(szfx0, lmin_, lmax_);	
		
	kxx=-2.5*log10( pow((1+z),-1)*(trpzx.Value()/trpzx0.Value()) ); 
	return kxx;
};


double PhotometryCalcs::CompColor(double z, ClassFunc1D& sed, Filter& filterX, Filter& filterY)
{
// compute color: mag in filter X - filter Y, for a galaxy with:
// -redshift z
// -spectral energy distribution sed (sed might have extinction applied already)

	//SEDzFilterProd returns: (sedred^zshift*filtX*lambda);
	SEDzFilterProd SEDFX(sed, filterX, z); 
	SEDzFilterProd SEDFY(sed, filterY, z); 
	// integrate SED across each filter 
	FilterIntegrator intsedX(SEDFX, lmin_, lmax_);
	FilterIntegrator intsedY(SEDFY, lmin_, lmax_);

	// To get correct color zeropoints: should use getFilterZeroPoint() here
	//FilterProd returns (filtX/lambda)
	FilterProd  FX(filterX); 
	FilterProd  FY(filterY); 
	FilterIntegrator intFX(FX, lmin_, lmax_);
	FilterIntegrator intFY(FY, lmin_, lmax_);
	double zpCxy = -2.5*log10(intFY.Value()/intFX.Value());

	// Color
	double Cxy = -2.5*log10(intsedX.Value()/intsedY.Value())+zpCxy;

	return Cxy;
};


double PhotometryCalcs::restFrameFlux(ClassFunc1D& sed, Filter& filter, double zs)
{
    BlueShiftFilter blueshiftFilter(filter, zs);// returns filter(lam*(1+zs))
    SEDzFilterProd sedXlambdaXfilter(sed, blueshiftFilter, 0.);// returns sed*lambda*filter
    // if last arg !=0 then sed would be redshifted (sed(lam/(1+z)))
	
	// integrate
	FilterIntegrator integrandSED(sedXlambdaXfilter, lmin_, lmax_);

	double zpFluxFilter = getFilterZeroPointFlux(filter);//3631e-26*
	//cout << zpFluxFilter <<"  ";
	
	double f0 = integrandSED.Value()/zpFluxFilter;
	//cout << integrandSED.Value() <<"  ";
	
	return f0;

};



// Filter methods


double PhotometryCalcs::getFilterZeroPointFlux(Filter& filterX)
{

    //FilterProd returns (filtX/lambda)
    FilterProd  FX(filterX);
    FilterIntegrator intFX(FX, lmin_, lmax_);
    double zp = intFX.Value();
    return zp;
};


double PhotometryCalcs::effectiveFilterWavelength(ClassFunc1D& filter)
{

    double x, y, lambdaEff, lmin, lmax;
    
    // Find the edges of the filter so can integrate across correct wavelength
    // range
    findFilterEdges(lmin, lmax, filter);
    //cout <<"     Filter edges are "<<lmin<<" to "<<lmax<<endl;

    // filter multiplied by lambda
    FilterXLambda  filterXLambda(filter);
    
    FilterIntegrator intFilterXLambda(filterXLambda,lmin,lmax);
    FilterIntegrator intFilter(filter,lmin,lmax);
    
    x = intFilterXLambda.Value();
    y = intFilter.Value();
  
    lambdaEff = x/y;
    
    //cout <<"     Filter "<<iFilter<<" effective wavelength = "<<lambdaEff<<endl;

    return lambdaEff;
};


double PhotometryCalcs::findFilterMax(ClassFunc1D& filter, double& lambdaAtMax, int nStep)
{

    // over whole range lmin_:lmax_
    double dStep = (lmax_ - lmin_)/(nStep - 1);
    
    double maxVal = 1e-6;
    int iRecord = -1;
    for (int i=0; i<nStep; i++) {
        double lam = lmin_ + i*dStep;
        double tmp = filter(lam);
        if (tmp > maxVal){
            iRecord = i;
            maxVal = tmp;
            }
        }
        
    if (iRecord<0)
        throw ParmError("Failed to find maximum of filter transmission");
        
    lambdaAtMax = lmin_ + iRecord*dStep;
    //cout <<"     Max transmission of filter is "<< maxVal <<endl; 
    //cout <<"     Wavelength at filter "<<iFilter<<" max is "<<lambdaAtMax<<endl;
    //cout << endl;
    
    return maxVal;

};


double PhotometryCalcs::findFilterTransValue(ClassFunc1D& filter, double trans, double lmin, 
                                                         double lmax, int nStep)
{

    // over given range lmin:lmax
    double dStep = (lmax - lmin)/(nStep - 1);
    double minDiff = 1e6;
    int iRecord = -1;
    double lambdaAtTrans = -999.;
    
    
    for (int i=0; i<nStep; i++) {
    
        double lam = lmin + i*dStep;
        double tmp = filter(lam);
        double diff = abs(tmp - trans);
            
        if (diff < minDiff){
            iRecord = i;
            minDiff = diff;
            }
        }
            
    lambdaAtTrans = lmin + iRecord*dStep;

          
    if (iRecord<0)
        throw ParmError("Failed to find filter transmission value");
        
    //cout <<"     Wavelength at filter "<<iFilter<<" transmission = "<<trans;
    //cout <<" is "<<lambdaAtTrans<<endl;
    //cout << endl;
    
    return lambdaAtTrans;

};


void PhotometryCalcs::findFilterEdges(double& lmin, double& lmax, ClassFunc1D& filter, 
	                                                double edgeDefinition)
{

    double lTMax;
    double filterMax = findFilterMax(filter, lTMax);
    double minTranmsission = edgeDefinition*filterMax;
    //cout <<"     wavelength of filter maximum is "<<lTMax<<endl;
    
    double lLow, lHigh;
    lLow = lmin_;
    lHigh = lTMax;
    //cout <<"     Searching between "<<lLow<<" and "<<lHigh<<endl;
    lmin = findFilterTransValue(filter, minTranmsission, lLow, lHigh);
    lLow = lTMax;
    lHigh = lmax_;
    //cout <<"     Now searching between "<<lLow<<" and "<<lHigh<<endl;
    lmax = findFilterTransValue(filter, minTranmsission, lLow, lHigh);

};


/******* SimData **************************************************************/

//constructor for when want to use this class to **calculate** k-corrections
SimData::SimData(vector<SED*> sedArray, vector<Filter*> filterArray, 
            SimpleUniverse& su, RandomGeneratorInterface& rg,
                                    int nEllipticals, int nSpirals)
: sedArray_(sedArray) , filterArray_(filterArray) , nEllipticals_(nEllipticals),
    nSpirals_(nSpirals) , su_(su) , rg_(rg)
{
	ebvmax_=0.3; ebvmaxEl_=0.1;
	nsed_ = sedArray_.size();
	nFilters_ = filterArray_.size();
	nStarbursts_ = nsed_ - nEllipticals_ - nSpirals_;
	
	isAddMadau_ = true;
	isLyC_ = true;
	isReadKcorr_ = false;
	
	setLSSTPars();
	
	if (nsed_>=1000)
	    throw ParmError("ERROR! Too many SEDs");
};

//constructor for when want to use this class to **read** k-corrections
SimData::SimData(SimpleUniverse& su, RandomGeneratorInterface& rg, 
    vector<SInterp2D*> kInterpZExt, int nFilters, int nEllipticals, int nSpirals)
: su_(su) , rg_(rg) , kInterpZExt_(kInterpZExt) , nFilters_(nFilters) , nEllipticals_(nEllipticals), nSpirals_(nSpirals)
{

	isReadKcorr_ = true;
	setLSSTPars();
	
};


/******* SimData methods ******************************************************/
                                               
// Simulation methods

// Simulate a "true" magnitude with or without Madau IGM absorption
double SimData::GetMag(double zs, double sedtype, double amag, double ext,
                                       int iFilterObs, Filter& restFrameFilter)
{

    Timer tm("timer",false);

    if (isReadKcorr_)
        throw ParmError("ERROR! Not set up to read k-correction from a file");

	// calculate distance modulus
	// @warning galaxy magnitudes for redshifts z<<0.01 don't have much meaning
	double mu;
	if (zs<1e-5) 
	    mu = 25.;
	else {
	    su_.SetEmissionRedShift(zs);
	    mu=5.*log10(su_.LuminosityDistanceMpc())+25;
	    }
	    
	// retrieve SED id of galaxy
	int sedID = returnSedId(sedtype);// Check this is doing the correct job

    // set correct reddening law
    int law = 0; // The Cardelli law
    if ( sedtype >= (nEllipticals_+nSpirals_) )
        law = 1; // The Calzetti law
    
    tm.Split();
    // copy contents of sedArray_[sedID] into new SED object
    SED sed(*(sedArray_[sedID]));
    
    // copy contents of filterArray_[iFilter] into a new Filter object
    Filter filter((*filterArray_[iFilterObs]));
    
    double kcorr = calcKcorr(sed, filter, restFrameFilter, zs, ext, law);
    /*if (isAddMadau_) {
        // add Madau absorption
        SEDMadau sedMadau(sed, zs);
    
        // redden SED
        SEDGOODSRedfix sedReddened(sedMadau,zs,ext,law);
		
        // k correction from REST FRAME filter to OBS FRAME filter
	    kcorr = Kcorr(zs,sedReddened,filter,restFrameFilter);    
	    }
	else{ 
	    // redden SED
        SEDGOODSRedfix sedReddened(sed,zs,ext,law);
		
        // k correction from REST FRAME filter to OBS FRAME filter
	    kcorr = Kcorr(zs,sedReddened,filter,restFrameFilter); 
        }*/
    tm.Split();
    //cout <<"Time to do k-corr part = "<< tm.PartialElapsedTimems() <<" ms"<<endl;
    
    tm.Split();
	// magnitude
	double mag =  amag + mu + kcorr;
	
	int isInf = my_isinf(kcorr);
    if (isInf != 0) // if kcorrection is infinite because filter has
        mag = 99;   // shifted out of where galaxy has flux values set 
                    // magnitude to 99 (undetected)	
                    
    int isMagInf = my_isinf(mag);
    if (isMagInf != 0) {
        cout <<"     Warning magnitude is infinite, z = "<< zs <<", k = "<< kcorr <<", amag = "<< amag;
        cout <<", mu = "<< mu <<endl;
        }
        
	return mag;
};

// Simulate a "true" magnitude with or without Madau IGM absorption
// READ K-CORRECTION FROM A FILE
double SimData::GetMag(double zs, double sedtype, double amag, double ext, int iFilterObs)
{
    Timer tm("timer",false);

    if (!isReadKcorr_)
        throw ParmError("ERROR! Not set up to read k-correction from a file");

	// calculate distance modulus
	su_.SetEmissionRedShift(zs);
	double mu=5*log10(su_.LuminosityDistanceMpc())+25;
	
	// retrieve SED id of galaxy
	int sedID = returnSedId(sedtype);// Check this is doing the correct job

    // get k correction
    tm.Split();
    double kcorr = interpKcorr(sedID, iFilterObs, zs, ext);
    tm.Split();
    cout <<"Time to do k-corr part = "<< tm.PartialElapsedTimems() <<" ms"<<endl;
    
	// magnitude
	double mag =  amag + mu + kcorr;
	
	int isInf = my_isinf(kcorr);
    if (isInf != 0) // if kcorrection is infinite because filter has
        mag = 99;   // shifted out of where galaxy has flux values set 
                    // magnitude to 99 (undetected)	
                    
    int isMagInf = my_isinf(mag);
    if (isMagInf != 0)
        cout <<"     Warning magnitude is infinite"<<endl;

	return mag;
};


// Simulate a "true" magnitude with a particular line of sight IGM absorption
double SimData::GetMag(double zs, double sedtype, double amag, double ext,
    int iFilterObs, Filter& restFrameFilter, IGMTransmission igmTransmission)
{
    if (isReadKcorr_)
        throw ParmError("ERROR! Not set up to read k-correction from a file");

	// calculate distance modulus
	su_.SetEmissionRedShift(zs);
	double mu=5*log10(su_.LuminosityDistanceMpc())+25;
	
	// retrieve SED id of galaxy
	int sedID = returnSedId(sedtype);// Check this is doing the correct job

    // set correct reddening law
    int law = 0; // The Cardelli law
    if ( sedtype >= (nEllipticals_+nSpirals_) )
        law = 1; // The Calzetti law
        
    // copy contents of sedArray_[sedID] into new SED object
    SED sed(*(sedArray_[sedID]));
    
    // copy contents of filterArray_[iFilter] into a new Filter object
    Filter filter((*filterArray_[iFilterObs]));
    
    double kcorr;
    
    // add IGM absorption
    SEDIGM sedIGM(sed, igmTransmission, zs);
    
    // redden SED
    SEDGOODSRedfix sedReddened(sedIGM,zs,ext,law);
		
    // k correction from REST FRAME filter to OBS FRAME filter
    kcorr = Kcorr(zs,sedReddened,filter,restFrameFilter);    

        
	// magnitude
	double mag =  amag + mu + kcorr;
	
	int isInf = my_isinf(kcorr);
    if (isInf != 0) // if kcorrection is infinite because filter has
        mag = 99;   // shifted out of where galaxy has flux values set 
                    // magnitude to 99 (undetected)	
                    
    int isMagInf = my_isinf(mag);
    if (isMagInf != 0)
        cout <<"     Warning magnitude is infinite"<<endl;
	
	return mag;
};


// add percentage magnitude error
vector<double> SimData::addError(double mag, double percentError, int iFilter)
{
    Filter filter((*filterArray_[iFilter]));
    double fluxError = percentError*convertABMagToFluxMaggies(mag, filter);
    vector<double> observation; // observed magnitude and magnitude error
    observation = addFluxError(mag, fluxError, iFilter);
	return observation;
};


// add LSST u band error
vector<double> SimData::addLSSTuError(double mag, int nVisits)
{   
    int iFilter = 0;
    double m5 = returnPointSource5sigmaDepth(uCm_,uMsky_,uTheta_,tVis_,
                                                            ukm_,airMass_);

    vector<double> obsmag = getObservedLSSTMagnitude(mag, m5, uGamma_, nVisits, iFilter);
	return obsmag;
};


// add LSST g band error
vector<double> SimData::addLSSTgError(double mag, int nVisits)
{   
    int iFilter = 1;
    double m5 = returnPointSource5sigmaDepth(gCm_,gMsky_,gTheta_,tVis_,
                                                            gkm_,airMass_);
    vector<double> obsmag = getObservedLSSTMagnitude(mag, m5, gGamma_, nVisits, iFilter);
	return obsmag;
};

	
// add LSST r band error
vector<double> SimData::addLSSTrError(double mag, int nVisits)
{   
    int iFilter = 2;
    double m5 = returnPointSource5sigmaDepth(rCm_,rMsky_,rTheta_,tVis_,
                                                            rkm_,airMass_);
    vector<double> obsmag = getObservedLSSTMagnitude(mag, m5, rGamma_, nVisits, iFilter);
	return obsmag;
};


// add LSST i band error
vector<double> SimData::addLSSTiError(double mag, int nVisits)
{   
    int iFilter = 3;
    double m5 = returnPointSource5sigmaDepth(iCm_,iMsky_,iTheta_,tVis_,
                                                            ikm_,airMass_);
    vector<double> obsmag = getObservedLSSTMagnitude(mag, m5, iGamma_, nVisits, iFilter);
	return obsmag;
};


// add LSST z band error
vector<double> SimData::addLSSTzError(double mag, int nVisits)
{   
    int iFilter = 4;
    double m5 = returnPointSource5sigmaDepth(zCm_,zMsky_,zTheta_,tVis_,
                                                            zkm_,airMass_);
    vector<double> obsmag = getObservedLSSTMagnitude(mag, m5, zGamma_, nVisits, iFilter);
	return obsmag;
};


// add LSST y band error
vector<double> SimData::addLSSTyError(double mag, int nVisits)
{   
    int iFilter = 5;
    double m5 = returnPointSource5sigmaDepth(yCm_,yMsky_,yTheta_,tVis_,
                                                            ykm_,airMass_);
    vector<double> obsmag = getObservedLSSTMagnitude(mag, m5, yGamma_, nVisits, iFilter);
	return obsmag;
};


// simulate an SED
double SimData::SimSED(int gtype)
// The 1st nEllipticals_ elements in sedArray correspond to elliptical galaxies
// The next (nEllipticals_+1) to (nSpirals_+nEllipticals_) correspond to spiral galaxies
// The next (nSpirals_+nEllipticals_+1) to nsed_ correspond to starbursts
// Also be aware: elliptical == early type; spiral == late type
// Also remember the zero indexing!
// WARNING: only can deal with <1000 SEDs
{

    if ( (gtype<1)||(gtype>3) )
        throw ParmError("ERROR! Unknown galaxy type");

    // to make the flat random number draw over a CLOSED interval
    // i.e. including the edges
	double fudge=0.5; 
			 
    // bounds of elliptical galaxy SEDs in sedArray_
	double E1 = 0 - fudge;
	double E2 = nEllipticals_-1 + fudge;
	// bounds of spiral galaxy SEDs in sedArray_
	double L1 = nEllipticals_ - fudge;
	double L2 = nEllipticals_ + nSpirals_-1 + fudge;
	// bounds of starburst galaxy SEDs in sedArray_
	double S1 = nEllipticals_ + nSpirals_ - fudge;
	double S2 = nsed_-1 + fudge;

	double sedtype;
	if (gtype<2) // if gal is elliptical type
		sedtype=1+round(E1+(E2-E1)*rg_.Flat01())/1000;
	else if (gtype>1 && gtype<3)// if gal is late type
		sedtype=2+round(L1+(L2-L1)*rg_.Flat01())/1000;
	else if (gtype>2)// if gal is starburst type
		sedtype=3+round(S1+(S2-S1)*rg_.Flat01())/1000;
    else
        throw ParmError("ERROR! Galaxy type not understood!");

    // sedtype is a number in the following form:
    // 1.XXX or 2.XXX or 3.XXX
    // The 1., 2. or 3. correspond to whether the galaxy was originally an 
    // elliptical, spiral or starburst galaxy type.
    // The XXX part is an integer from 0 to 999 (though always with 3 signficant
    // figures, padded with zeros when necessary) and corresponds to the exact
    // SED of the galaxy, found in sedArray e.g sedArray[0], sedArray[1] ....
	return sedtype;
};


// simulate reddening amount
double SimData::SimRed(double type)
{
		double ebv;
	if (type<2)
		ebv=ebvmaxEl_*rg_.Flat01();
	else
		ebv=ebvmax_*rg_.Flat01();
			
	return ebv;
};


// Calculation methods

vector<double> SimData::returnFilterRFWavelengths(double zs)
{

    vector<double> lambdaRFs;
    for (int i=0; i<nFilters_; i++) {
    
        // current filter
        Filter filter((*filterArray_[i]));
        
        // find wavelength of filter in galaxy's restframe
        double lamObs = effectiveFilterWavelength(filter);
        //cout << lamObs <<"  ";
        double lamRF = returnRestFrameWaveLength(lamObs,zs);
        lambdaRFs.push_back(lamRF);
        
        }
    //cout << endl;
        
    return lambdaRFs;

};


TArray<double> SimData::returnSEDFluxesInRestFrame(double zs)//, Filter& restFrameFilter)
{
    TArray<double> templateFluxes;
    int ndim = 2;
    sa_size_t mydim[ndim];
    mydim[0]=nsed_; mydim[1] = nFilters_;
    templateFluxes.SetSize(ndim, mydim);
    
    for (int i=0; i<nsed_; i++) {
        SED sed((*sedArray_[i]));
    
        for (int j=0; j<nFilters_; j++) {
        
            Filter filter((*filterArray_[j]));
            //double k = Kcorr(zs,i,j,restFrameFilter);
            //double flux = pow(10.,-0.4*k);
            
            double flux = restFrameFluxLambda(sed,filter,zs);
            templateFluxes(i,j) = flux;
            }
        }
    /*for (int i=0; i<nsed_; i++) {
        for (int j=0; j<nFilters_; j++) {
            double lam = lambdaRFs[j];
            double val = sedArray_[i]->returnFlux(lam);
            templateFluxes(i,j) = val;
            }
         }*/
         
    return templateFluxes;
};

// INTERNAL FUNCTIONS

double SimData::calcKcorr(SED& sed, Filter& filter, Filter& restFrameFilter, double zs, double ext, int law)
{
    double kcorr;
    if (isAddMadau_) {
        // add Madau absorption
        SEDMadau sedMadau(sed, zs, isLyC_);
    
        // redden SED
        SEDGOODSRedfix sedReddened(sedMadau,zs,ext,law);
		
        // k correction from REST FRAME filter to OBS FRAME filter
	    kcorr = Kcorr(zs,sedReddened,filter,restFrameFilter);    
	    }
	else{ 
	    // redden SED
        SEDGOODSRedfix sedReddened(sed,zs,ext,law);
		
        // k correction from REST FRAME filter to OBS FRAME filter
	    kcorr = Kcorr(zs,sedReddened,filter,restFrameFilter); 
        }

    return kcorr;
};

double SimData::interpKcorr(int sedID, int iFilterObs, double zs, double ext)
{
    int id = returnLinearIndex(sedID,iFilterObs,nFilters_);
    double kcorr = kInterpZExt_[id]->operator()(zs,ext);
    return kcorr;

};

double SimData::interpKcorr(int linearIndex, double zs, double ext)
{
    double kcorr = kInterpZExt_[linearIndex]->operator()(zs,ext);
    return kcorr;
};

int SimData::returnSedId(double sedtype)
{
    double eps = 1e-6;
    
    int sedID;
    
    if( (floor(sedtype)-1)<eps )        // If galaxy is elliptical type
			sedID=(int)round(1000*sedtype-1000);
    else if ( (floor(sedtype)-2)<eps )  // If galaxy is spiral type
            sedID=(int)round(1000*sedtype-2000);
    else if ( (floor(sedtype)-3)<eps )  // If galaxy is starburst type
            sedID=(int)round(1000*sedtype-3000);
    else
        throw ParmError("Unknown galaxy SED");
        
    return sedID;
};


vector<double> SimData::addFluxError(double mag, double fluxError, int iFilter)
{
    Filter filter((*filterArray_[iFilter]));
	double flux = convertABMagToFluxMaggies(mag, filter); // flux in FREQ units
	double fluxobs = flux+fluxError*rg_.Gaussian();       // flux in FREQ units
	
	
	double obsmag = convertFluxMaggiesToABMag(fluxobs, filter);
	if (my_isnan(obsmag))
	    cout <<"obs mag is nan, fluxobs = "<< fluxobs << endl;
	double fE = fluxError/flux;
	double magError = convertFluxErrorToMagError(fE);// should be fluxError/flux?????
	
	vector<double> observation;
	observation.push_back(obsmag);
	observation.push_back(magError);
	
	return observation;
};


vector<double> SimData::getObservedLSSTMagnitude(double mag, double m5, double gamma, int nVis, int iFilter)
{
    Filter filter((*filterArray_[iFilter]));

    // get random component of photometric error
    double x = returnX(mag, m5);
    double sigSq = returnLSSTRandomErrorSq(x,gamma,nVis);
    
    // total photometric error (in magnitudes)
    double sigmaM = sqrt(sigSq + sigmaSys_*sigmaSys_);
    
    // convert photometric error to flux units
    double flux = convertABMagToFluxMaggies(mag, filter);
	double fluxError = convertMagErrorToFluxError(sigmaM, flux); 
	
	// add the error
    vector<double> observation;
    observation = addFluxError(mag, fluxError, iFilter);
    
	return observation;
};


double SimData::returnLSSTRandomErrorSq(double x, double gamma, double nVis)
{

    double sigmaRandsq = (0.04 - gamma)*x + gamma*x*x;
    // in magnitudes^2
    return sigmaRandsq/nVis;

};


double SimData::returnPointSource5sigmaDepth(double Cm, double msky,
                                double theta, double tvis, double km, double X)
{

    double m5 = Cm + 0.50*(msky - 21) + 2.5*log10(0.7/theta) +
                                        1.25*log10(tvis/30) - km*(X - 1);
    return m5;                              
};


void SimData::setLSSTPars()
{
    // expected median sky zenith brightness at Cerro Pachon, assuming mean solar cycle
    // and three-day old Moon (mag/arcsec^2)
    uMsky_ = 21.8; gMsky_ = 22.0; rMsky_ = 21.3; iMsky_ = 20.0; zMsky_ = 19.1; 
    yMsky_ = 17.5;
    
    // expected delivered median zenith seeing (arcsec). For larger airmass X seeing
    // is proportional to X^0.6
	uTheta_ = 0.77; gTheta_ = 0.73; rTheta_ = 0.70; iTheta_ = 0.67; zTheta_ = 0.65; 
	yTheta_ = 0.63;
	
	// band dependent parameter
	uGamma_ = 0.037; gGamma_ = 0.038; rGamma_ = 0.039; iGamma_ = 0.039; zGamma_ = 0.040;
	yGamma_ = 0.040;
	
	// band dependent parameter
	uCm_ = 23.60; gCm_ = 24.57; rCm_ = 24.57; iCm_ = 24.47; zCm_ = 24.19; yCm_ = 23.74;
	
	// adopted atmospheric extinction
	ukm_ = 0.48; gkm_ = 0.21; rkm_ = 0.1; ikm_ = 0.07; zkm_ = 0.06; ykm_ = 0.06;
	
	// exposure time, 2 back-to-back 15s exposures
	tVis_ = 2*15;
	
	// median air mass
	airMass_ = 1.2;

    // systematic error
    sigmaSys_ = 0.005;
};

/******* ReadKCorrections methods ********************************************/



// When interpolating across redshifts and extinctions
void ReadKCorrections::readInterpZExt(int nSED, int nFilter)
{
            
    int index = 0;
    for (int i=0; i<nSED; i++) {
        for (int j=0; j<nFilter; j++) {
    
            string fname = getFileName(i,j);

            ifstream ifs;
	        ifs.open(fname.c_str(), ifstream::in);
	        if (ifs.fail()) { 
	            string emsg = "ERROR: failed to find k-correction file " + fname;
		        throw ParmError(emsg);
		        }
	        sa_size_t nr, nc;
	        TArray<double> tab;
	        tab.ReadASCII(ifs,nr,nc);
	        ifs.close();

	        // row direction of tab corresponds to extinction dimension
	        // col direction of tab corresponds to redshift dimension
	        // THIS IS THE OPPOSITE DIRECTION TO THAT IN THE FILE!
	
	        // and: doing transpose tab and swapping evals, zvals order will speed this up
	        // (only because ne<nz)
	        TArray<double> tabT=transposeTable(tab);

            // Store interp table
		    kInterpZExt_.push_back( new SInterp2D() );
            kInterpZExt_[index] -> definePoints(zvals_, evals_, tabT, true);
            
            index++;
            }
        }
};

// When just want to interpolate across redshifts
void ReadKCorrections::readInterpZ(int nSED, int nFilter)
{

    int index = 0;
    for (int i=0; i<nSED; i++) {
        for (int j=0; j<nFilter; j++) {
    
            string fname = getFileName(i,j);

            ifstream ifs;
	        ifs.open(fname.c_str(), ifstream::in);
	        if (ifs.fail()) { 
	            string emsg = "ERROR: failed to find k-correction file " + fname;
		        throw ParmError(emsg);
		        }
	        sa_size_t nr, nc;
	        TArray<double> tab;
	        tab.ReadASCII(ifs,nr,nc);
	        ifs.close();

	        // row direction of tab corresponds to extinction dimension
	        // col direction of tab corresponds to redshift dimension
	        // THIS IS THE OPPOSITE DIRECTION TO THAT IN THE FILE!
	
	        // and: doing transpose tab and swapping evals, zvals order will speed this up
	        TArray<double> tabT=transposeTable(tab);
	        
	        
	        for (int k=0; k<ne_; k++) {
	        
	            // turn array column into a vector
	            vector<double> vec;
	            for (int a=0; a<nz_; a++)
	                vec.push_back(tabT(a,k));

                // Store interp table
                kInterpZ_.push_back( new SInterp1D() );
                kInterpZ_[index] -> DefinePoints(zvals_, vec);
                
                index++;
                }
        }
    }   
};

string ReadKCorrections::getFileName(int iSED, int iFilter)
{
    stringstream ss1,ss2,ss3,ss4,ss5,ss6,ss7;
    ss1 << iSED; 
    ss2 << iFilter;
    ss3 << zmin_; ss4 << zmax_; ss5 << nz_;
    ss6 << emax_; ss7 << ne_;
    string fname;
    fname += "kCorrections/kCorr_" + sedLib_ + "sed" + ss1.str();
    fname += "_" + filtSet_ + "filt" + ss2.str() + "_" + restFrameFilt_;
    fname += "_zmin" + ss3.str() + "_zmax" + ss4.str() + "_nz" + ss5.str();
    fname += "_emax" + ss6.str() + "_ne" + ss7.str();
    // madau absorption?
    if (isMadau_)
        fname = fname + "_wMadau.txt";
    else
        fname = fname + "_woMadau.txt";
        
    return fname;

};



/******* TemplateChiSquare methods ********************************************/

TArray<double> TemplateChiSquare::galaxyChiSquared(vector<double> obs, 
         vector<double> errors, double zs, int& sedBestFit, double& normBestFit)
{

    int no = obs.size();
    int ne = errors.size();
    if ( (no!=nFilters_)||(ne!=nFilters_) )
        throw ParmError("ERROR! Number of filters does not match amount of data");
        
    TArray<double> chiSquare;
    int ndim = 2;
    sa_size_t mydim[ndim];
    mydim[0]=nsed_; mydim[1] = nA_;
    chiSquare.SetSize(ndim, mydim);

    for (int i=0; i<nsed_; i++) {
        for (int j=0; j<nA_; j++) {
        
            double a = aMin_ + dA_*j;
            //cout << a<<endl;
            chiSquare(i,j) = meritFunction(obs, errors, zs, i, a);
            }
        }
        
    // Get best-fit parameters
    int iSED, iA;
    double minChisq = findMinimumPosition(chiSquare,iSED,iA);
    sedBestFit = iSED;
    normBestFit = aMin_ + dA_*iA;
    
    cout <<"     Min chi-square/dof = "<<minChisq/nFilters_<<", SED = "<<sedBestFit;
    cout <<", A = "<<normBestFit<<endl;

    return chiSquare;
};



double TemplateChiSquare::meritFunction(vector<double> obs, vector<double> errors, double zs, int iSED, double Anorm)
{

    su_.SetEmissionRedShift(zs);
    double dL = su_.LuminosityDistanceMpc();
    
    SED sed((*sedArray_[iSED]));
    Filter arbitraryFilter((*filterArray_[2]));
    
    double chisq = 0.;
    for (int i=0; i<nFilters_; i++) {
    
        Filter filter((*filterArray_[i]));
    
        double fluxPred = Anorm*restFrameFluxLambda(sed,filter,zs);
        double fluxObs = convertABMagToFluxLambda(obs[i],zs,dL,sed,arbitraryFilter,filter); 
        double fluxError = convertMagErrorToFluxError(errors[i],fluxObs);
        
        double val = (fluxObs - fluxPred)/fluxError;
        chisq+=(val*val);
        }

    return chisq;

};
                                                                

