#include "simdata.h"

/******* PhotometryCalcs methods **********************************************/


// K-correction calculations 
// Must make SED wrapper class that is NOT ClassFunc1D!

double PhotometryCalcs::Kcorr(double z, SpecEnergyDist& sed, Filter& filterX, Filter& restFrameFilter)
{
// rest-frame and observed filters need NOT be the SAME
// compute general K correction for a galaxy with:
// -redshift z
// -spectral energy distribution sed (sed might have extinction applied already)
// -from REST FRAME bandpass Y
// -to OBSERVED bandpass X

    Timer tm("timer",false);
    
    //cout << filterX.Npoints() <<"  "<< restFrameFilter.Npoints() << endl;

	double kxy;
		
   // tm.Split();
	// returns sedred(lambda/(1+z))*filtX(lambda)*lambda;
	SEDzFilterProd szfx(sed, filterX, z);
	// returns filtY(lambda)/lambda  
	FilterProd  fprody(restFrameFilter);          
	// returns sed(lambda)*filtY(lambda)*lambda;    
	SEDzFilterProd szfB0(sed, restFrameFilter, 0.);
	 // returns filtX(lambda)/lambda
	FilterProd  fprodx(filterX);
	//tm.Split();
    //cout <<"Time to set up k-correction = "<< tm.PartialElapsedTimems() <<" ms, ";
	
	tm.Split();
	for (int i=0; i<100000; i++)
	    double tmp = szfx(4e-7);
	tm.Split();
	cout <<"Time of test1 = "<< tm.PartialElapsedTimems() <<" ms, ";
	
	tm.Split();
	for (int i=0; i<100000; i++)
	    double tmp = szfB0(4e-7);
	tm.Split();
	cout <<"Time of test1a = "<< tm.PartialElapsedTimems() <<" ms, ";
	
	//tm.Split();
	for (int i=0; i<100000; i++)
        double tmp = fprody(4e-7);
	//tm.Split();
	//cout <<"Time of test2 = "<< tm.PartialElapsedTimems() <<" ms, "<< endl;
	
    // To integrate the above
    
	FilterIntegrator trpzx(szfx, lmin_, lmax_, npt_); //tm.Split(); //
	//tm.Split();
	//cout <<"Time to int1 = "<< tm.PartialElapsedTimems() <<" ms, "<< endl;
	FilterIntegrator trpzfy(fprody, lmin_, lmax_, npt_); //tm.Split(); //cout <<"Time to int2 = "<< tm.PartialElapsedTimems() <<" ms, ";
	//tm.Split();
	FilterIntegrator trpzy(szfB0, lmin_, lmax_, npt_); //tm.Split(); //cout <<"Time to int3 = "<< tm.PartialElapsedTimems() <<" ms, ";
	//tm.Split();
	FilterIntegrator trpzfx(fprodx, lmin_, lmax_, npt_); //tm.Split(); //cout <<"Time to int4 = "<< tm.PartialElapsedTimems() <<" ms, ";
		
    //tm.Split();
    double a = trpzx.Value(); // this takes a long time
    //tm.Split();
	//cout <<"Time to do a = "<< tm.PartialElapsedTimems() <<" ms, ";
	//tm.Split();
    double b = trpzfx.Value();
    //tm.Split();
	//cout <<"Time to do b = "<< tm.PartialElapsedTimems() <<" ms, ";
	//tm.Split();
    double c = trpzfy.Value();
    //tm.Split();
	//cout <<"Time to do c = "<< tm.PartialElapsedTimems() <<" ms, ";
	//tm.Split();
    double d = trpzy.Value(); // this takes a long time
    //tm.Split();
	//cout <<"Time to do d = "<< tm.PartialElapsedTimems() <<" ms, ";
	kxy=-2.5*log10(pow((1.+z),-1)*(a/b)*(c/d));
	//tm.Split();
	//cout <<"Time to do end = "<< tm.PartialElapsedTimems() <<" ms"<<endl;
	//if (z<0.0001)
	//    cout << z <<"  "<< kxy <<endl;
	return kxy;
};


double PhotometryCalcs::Kcorr1Filter(double z, SpecEnergyDist& sed, Filter& filterX)
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
	FilterIntegrator trpzx(szfx, lmin_, lmax_, npt_);	
	FilterIntegrator trpzx0(szfx0, lmin_, lmax_, npt_);	
		
	kxx=-2.5*log10( pow((1+z),-1)*(trpzx.Value()/trpzx0.Value()) ); 
	return kxx;
};


double PhotometryCalcs::CompColor(double z, SpecEnergyDist& sed, Filter& filterX, Filter& filterY)
{
// compute color: mag in filter X - filter Y, for a galaxy with:
// -redshift z
// -spectral energy distribution sed (sed might have extinction applied already)

	//SEDzFilterProd returns: (sedred^zshift*filtX*lambda);
	SEDzFilterProd SEDFX(sed, filterX, z); 
	SEDzFilterProd SEDFY(sed, filterY, z); 
	// integrate SED across each filter 
	FilterIntegrator intsedX(SEDFX, lmin_, lmax_, npt_);
	FilterIntegrator intsedY(SEDFY, lmin_, lmax_, npt_);

	// To get correct color zeropoints: should use getFilterZeroPoint() here
	//FilterProd returns (filtX/lambda)
	FilterProd  FX(filterX); 
	FilterProd  FY(filterY); 
	FilterIntegrator intFX(FX, lmin_, lmax_, npt_);
	FilterIntegrator intFY(FY, lmin_, lmax_, npt_);
	double zpCxy = -2.5*log10(intFY.Value()/intFX.Value());
	//cout << intsedY.Value() <<"  "<< intsedX.Value() <<"  "<<-2.5*log10(intsedX.Value()/intsedY.Value()) <<"  "<<zpCxy << endl;

	// Color
	double Cxy = -2.5*log10(intsedX.Value()/intsedY.Value())+zpCxy;

	return Cxy;
};


double PhotometryCalcs::restFrameFlux(SpecEnergyDist& sed, Filter& filter, double zs)
{
 /* I think there is a problem with this - Matt Kirby 
     
    BlueShiftFilter blueshiftFilter(filter, zs);// returns filter(lam*(1+zs))
    SEDzFilterProd sedXlambdaXfilter(sed, blueshiftFilter, 0.);// returns sed*lambda*filter
    // if last arg !=0 then sed would be redshifted (sed(lam/(1+z)))
	
	// integrate
	FilterIntegrator integrandSED(sedXlambdaXfilter, lmin_, lmax_);

	double zpFluxFilter = getFilterZeroPointFlux(filter);//3631e-26*
	//cout << zpFluxFilter <<"  ";
	
	double f0 = integrandSED.Value()/zpFluxFilter;
	//cout << integrandSED.Value() <<"  ";
// */

// /* Matt's code

    SEDzFilterProd sedXlambdaXfilter(sed, filter, zs);

    //integrate
    FilterIntegrator integrandSED(sedXlambdaXfilter, lmin_, lmax_, npt_);

    double zpFluxFilter = getFilterZeroPointFlux(filter);

    double f0 = integrandSED.Value()/zpFluxFilter;
	
	return f0;
// */
};



// Filter methods


double PhotometryCalcs::getFilterZeroPointFlux(Filter& filterX)
{

    //FilterProd returns (filtX/lambda)
    FilterProd  FX(filterX);
    FilterIntegrator intFX(FX, lmin_, lmax_, npt_);
    double zp = intFX.Value();
    return zp;
};


double PhotometryCalcs::effectiveFilterWavelength(Filter& filter)
{

    double x, y, lambdaEff, lmin, lmax;
    
    // Find the edges of the filter so can integrate across correct wavelength
    // range
    findFilterEdges(lmin, lmax, filter);
    //cout <<"     Filter edges are "<<lmin<<" to "<<lmax<<endl;

    // filter multiplied by lambda
    FilterXLambda  filterXLambda(filter);
    
    FilterIntegrator intFilterXLambda(filterXLambda, lmin, lmax, npt_);
    FilterIntegrator intFilter(filter, lmin, lmax, npt_);
    
    x = intFilterXLambda.Value();
    y = intFilter.Value();
  
    lambdaEff = x/y;
    
    //cout <<"     Filter "<<iFilter<<" effective wavelength = "<<lambdaEff<<endl;

    return lambdaEff;
};


double PhotometryCalcs::findFilterMax(Filter& filter, double& lambdaAtMax, int nStep)
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


double PhotometryCalcs::findFilterTransValue(Filter& filter, double trans, double lmin, 
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


void PhotometryCalcs::findFilterEdges(double& lmin, double& lmax, Filter& filter, double edgeDefinition)
{

    // first find wavelength at filter maximum
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



/******* SimData methods ******************************************************/
                                               
// Simulation methods

// Full calculation of galaxy magnitude (rest-frame filter *is* one of the filters in filterArray_)
double SimData::getMag(double zs, double absmag, int sedid, int iobsfilt, int irestfilt, 
                       IGMTransmission& igmtrans, double ext, dustLaw law)
{

    // copy contents of filterArray_[irestfilt] into a new Filter object
    Filter rfFilter((*filterArray_[irestfilt])); // restframe filter of absmag
    
    return getMag(zs, absmag, sedid, iobsfilt, rfFilter, igmtrans, ext, law);

};


// Full calculation of galaxy magnitude (rest-frame filter *is not* one of the filters in filterArray_)
double SimData::getMag(double zs, double absmag, int sedid, int iobsfilt, Filter& rfFilter, 
                       IGMTransmission& igmtrans, double ext, dustLaw law)
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
	    mu = 5.*log10(su_.LuminosityDistanceMpc()) + 25.;
	    }
	    
	// copy contents of sedArray_[sedid] into new SED object
    SED sed(*(sedArray_[sedid])); // SED (no internal dust or IGM applied yet)
    
    // copy contents of filterArray_[iobsfilt] into a new Filter object
    Filter obsFilter((*filterArray_[iobsfilt])); // observation filter
    
    tm.Split();
    double kcorr = calcKcorr(sed, obsFilter, rfFilter, zs, igmtrans, ext, law);
    tm.Split();
    cout <<"Time to do k-corr part = "<< tm.PartialElapsedTimems() <<" ms"<<endl;
    
	// magnitude
	double mag =  absmag + mu + kcorr;
	
	// Check for infinities
	int isInf = my_isinf(kcorr);
    if (isInf != 0) // if kcorrection is infinite because filter has
        mag = 99;   // shifted out of where galaxy has flux values set 
                    // magnitude to 99 (undetected)	
                    
    int isMagInf = my_isinf(mag);
    if (isMagInf != 0) {
        cout <<"     Warning magnitude is infinite, z = "<< zs <<", k = "<< kcorr <<", amag = "<< absmag;
        cout <<", mu = "<< mu <<endl;
        }
        
	return mag;

};


// Calculate magnitude of a galaxy using pre-calculated k-correction tables
double SimData::getMag(double zs, double absmag, int sedid, int iobsfilt, int irestfilt, igmModel igm, 
                  double ext, dustLaw law)
{

    Timer tm("timer",false);

    if (!isReadKcorr_)
        throw ParmError("ERROR! Not set up to read k-correction from a file");

	// calculate distance modulus
	// @warning galaxy magnitudes for redshifts z<<0.01 don't have much meaning
	double mu;
	if (zs<1e-5) 
	    mu = 25.;
	else {
	    su_.SetEmissionRedShift(zs);
	    mu = 5.*log10(su_.LuminosityDistanceMpc()) + 25.;
	    }
	    
	// build k-correction filename to read in
	
	// SED name
	string sedname = sedNames_[sedid];
	
	// IGM type
	string igmType;
	if (igm == None)
	    igmType = "None";
	else if (igm == Mad)
	    igmType = "Madau";
	else if (igm == Aver)
	    igmType = "Mean";
	else {
	    cout <<"IGM model type = "<< igm << endl;
	    throw ParmError("ERROR! IGM type not understood");
	    }
	    
	// dust type
	string dustType;
	if (law == NoDust)
	    dustType = "NoDust";
	else if (law == Card)
	    dustType = "Card";
	else if (law == Calz)
	    dustType = "Calz";
	else {
	    cout <<"Dust model type = "<< law << endl;
	    throw ParmError("ERROR! Dust law not understood");
	    }
	    
	// dust amount
    char* dchar = new char[8];
    float dval=ext;
    sprintf(dchar, "%0.2f", dval);
    stringstream ss;
    ss<<dchar;
    
	string filename = sedname + "_" + filterSet_ + "_igm" + igmType + "_ext" + dustType + ss.str() + ".dat";
	    
	double kcorr;

    // magnitude
	double mag =  absmag + mu + kcorr;
	
	// Check for infinities
	int isInf = my_isinf(kcorr);
    if (isInf != 0) // if kcorrection is infinite because filter has
        mag = 99;   // shifted out of where galaxy has flux values set 
                    // magnitude to 99 (undetected)	
                    
    int isMagInf = my_isinf(mag);
    if (isMagInf != 0) {
        cout <<"     Warning magnitude is infinite, z = "<< zs <<", k = "<< kcorr <<", amag = "<< absmag;
        cout <<", mu = "<< mu <<endl;
        }
        
	return mag;
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

double SimData::calcKcorr(SED& sed, Filter& filter, Filter& restFrameFilter, double zs, 
                          IGMTransmission& igmtrans, double ext, dustLaw law)
{
    Timer tm("timer",false);
    
    // Check if need to add extinction
    if (ext>0) {
    
        // redden SED   
        int lawtype;
        if (law == Card)
            lawtype = 0;
	    else if (law == Calz)
	        lawtype = 1;
        else {
	        cout <<"Dust model type = "<< law << endl;
	        throw ParmError("ERROR! Dust law not understood");
	        }
        
        sed.doRedden(ext, lawtype);
        }
    
    
    // add IGM absorption
    SEDIGM sedIGM(sed, igmtrans, zs);


    // k correction from REST FRAME filter to OBS FRAME filter
    tm.Split();
	double kcorr = Kcorr(zs, sedIGM, filter, restFrameFilter);
	tm.Split();
    cout <<"Time to calc k-correction = "<< tm.PartialElapsedTimems() <<" ms"<<endl;


    return kcorr;
};


double SimData::interpKcorr(string filename, double zs, int iobsfilt, int irestfilt)
{

    // interpolate kcorr value 
    SInterp1D kCorrData;
    kCorrData.ReadXYFromFile(filename, 0., zs+1., npt_);
    double sedpart = kCorrData(zs);
    
    
    //------------------------------------------------------------------
    // integrate obs and restframe filters -> make a PhotometryCalcs method
    Filter obsFilter((*filterArray_[iobsfilt])); // observation filter
    Filter rfFilter((*filterArray_[irestfilt])); // rest-frame filter
    
    FilterProd  obsFiltIntegrand(obsFilter);
    FilterProd  rfFiltIntegrand(rfFilter);
    
    FilterIntegrator intObsFilt(obsFiltIntegrand, lmin_, lmax_, npt_);
    FilterIntegrator intRfFilt(rfFiltIntegrand, lmin_, lmax_, npt_);
    
    double filterpart = intRfFilt.Value()/intObsFilt.Value();
    //------------------------------------------------------------------
    
    // calc final k correction
    double kcorr = -2.5*log10((1./(1.+zs)) * sedpart * filterpart);

    return kcorr;
};

/*double SimData::interpKcorr(int sedID, int iFilterObs, double zs, double ext)
{
    int id = returnLinearIndex(sedID,iFilterObs,nFilters_);
    double kcorr = kInterpZExt_[id]->operator()(zs,ext);
    return kcorr;

};

double SimData::interpKcorr(int linearIndex, double zs, double ext)
{
    double kcorr = kInterpZExt_[linearIndex]->operator()(zs,ext);
    return kcorr;
};*/

/****** SimObservations methods *****************************************************************************/

vector<double> SimObservations::returnFilterRFWavelengths(double zs)
{

    vector<double> lambdaRFs;
    for (int i=0; i<nFilters_; i++) {
    
        // current filter
        Filter filter((*filterArray_[i]));
        
        // find wavelength of filter in galaxy's restframe
        double lamObs = effectiveFilterWavelength(filter);
        //cout << lamObs <<"  ";
        double lamRF = returnRestFrameWaveLength(lamObs, zs);
        lambdaRFs.push_back(lamRF);
        
        }
    //cout << endl;
        
    return lambdaRFs;

};


vector<double> SimObservations::addFluxError(double flux, double fluxError, int iFilter)
{
    
    Filter filter((*filterArray_[iFilter]));
    
    // convert magnitude to flux
	//double flux = convertABMagToFluxMaggies(mag, filter); // flux in FREQ units
	//cout <<", flux="<< flux;
	
	// scatter flux according to flux error size
	double fluxobs = flux + fluxError*rg_.Gaussian();       // flux in FREQ units
	//if (fluxobs<0)
	//    cout <<"     WARNING! negative flux "<< fluxobs << endl;
	
	// convert flux back to magnitude
	double obsmag = convertFluxMaggiesToABMag(fluxobs, filter);
	if (my_isnan(obsmag))
	    cout <<"obs mag is nan, fluxobs = "<< fluxobs << endl;
	//cout <<", obsmag="<< obsmag << endl;
	
	// convert flux error back to magnitude error
	// don't actually need to do this because the flux error was probably calculated *from* the mag error
	double fE = fluxError/flux;
	double magError = convertFluxErrorToMagError(fE);
	
	// package ready to return
	vector<double> observation;
	observation.push_back(obsmag);
	observation.push_back(magError);
	
	return observation;
};


vector<double> SimObservations::getObservedLSSTMagnitude(double mag, double m5, double gamma, int nYear, int iFilter)
{
    Filter filter((*filterArray_[iFilter]));

    // get random component of photometric error
    double x = returnX(mag, m5);
    double sigSq = returnLSSTRandomErrorSq(x, gamma); // eqn 3.2 in Science Book (not /nVis)
    
    // total photometric error (in magnitudes)
    double sigmaM = sqrt(sigSq + sigmaSys_*sigmaSys_); // eqn 3.1 in Science Book

    // convert photometric error to flux units
    double flux = convertABMagToFluxMaggies(mag, filter);
	double fluxError = convertMagErrorToFluxError(sigmaM, flux); 

	// add the error
    vector<double> observation;
    observation = addFluxError(flux, fluxError, iFilter);
    
    // don't do this here - post-process elsewhere ????
    // if very fainter than full visit depth m5 count as non-detection 
    // (mag=99, error on mag=1-sigma detection limit, full visit depth)
    double m5nvisit = m5single_[iFilter] + 1.25*log10(nYear*nVisYear_[iFilter]);
    double m1sig = m5nvisit + 2.5*log10(5.);
    
    if (observation[0]>m1sig) {
        observation[0]=99.;
        observation[1]=m1sig;
        }
    
    // also return the 1-sigma limiting mag
    observation.push_back(m1sig);
    
	return observation;
};



void SimObservations::setLSSTPars()
{
    // Number of visits per year (Table 1, Ivezic et al 2008)
    nVisYear_.push_back(6);
    nVisYear_.push_back(8);
    nVisYear_.push_back(18);
    nVisYear_.push_back(18);
    nVisYear_.push_back(16);
    nVisYear_.push_back(16);

    // Numbers below are from Table 3.2 in the Science Book

    // expected median sky zenith brightness at Cerro Pachon, assuming mean solar cycle
    // and three-day old Moon (mag/arcsec^2)
    Msky_.push_back(21.8);  // 22.9
    Msky_.push_back(22.0);  // 22.3
    Msky_.push_back(21.3);  // 21.2
    Msky_.push_back(20.0);  // 20.5
    Msky_.push_back(19.1);  // 19.6
    Msky_.push_back(17.5);  // 18.6
    //uMsky_ = 21.8; gMsky_ = 22.0; rMsky_ = 21.3; iMsky_ = 20.0; zMsky_ = 19.1; 
    //yMsky_ = 17.5;
    
    // expected delivered median zenith seeing (arcsec). For larger airmass X seeing
    // is proportional to X^0.6
    Theta_.push_back(0.77); 
    Theta_.push_back(0.73); 
    Theta_.push_back(0.70); 
    Theta_.push_back(0.67); 
    Theta_.push_back(0.65); 
	Theta_.push_back(0.63);
	//uTheta_ = 0.77; gTheta_ = 0.73; rTheta_ = 0.70; iTheta_ = 0.67; zTheta_ = 0.65; 
	//yTheta_ = 0.63;
	
	// band dependent parameter
	Gamma_.push_back(0.037); 
	Gamma_.push_back(0.038); 
	Gamma_.push_back(0.039); 
	Gamma_.push_back(0.039); 
	Gamma_.push_back(0.040);
	Gamma_.push_back(0.040);
	//uGamma_ = 0.037; gGamma_ = 0.038; rGamma_ = 0.039; iGamma_ = 0.039; zGamma_ = 0.040;
	//yGamma_ = 0.040;
	
	// band dependent parameter
	Cm_.push_back(23.60); // 22.92
	Cm_.push_back(24.57); // 24.29
	Cm_.push_back(24.57); // 24.33
	Cm_.push_back(24.47); // 24.20
	Cm_.push_back(24.19); // 24.07
	Cm_.push_back(23.74); // 23.69
	//uCm_ = 23.60; gCm_ = 24.57; rCm_ = 24.57; iCm_ = 24.47; zCm_ = 24.19; yCm_ = 23.74;
	
	// adopted atmospheric extinction
	km_.push_back(0.48); // 0.451
	km_.push_back(0.21); // 0.163
	km_.push_back(0.1);  // 0.087
	km_.push_back(0.07); // 0.065
	km_.push_back(0.06); // 0.043
	km_.push_back(0.06); // 0.138
	//ukm_ = 0.48; gkm_ = 0.21; rkm_ = 0.1; ikm_ = 0.07; zkm_ = 0.06; ykm_ = 0.06;
	
	// exposure time, 2 back-to-back 15s exposures
	tVis_ = 2.*15.;
	
	// median air mass
	airMass_ = 1.2;
	
	// pre-calculate SINGLE VISIT m5 here
	int nvis = 1;
	for (int i=0; i<Cm_.size(); i++) {
	    double m5 = calcPointSource5sigmaDepth(Cm_[i], Msky_[i], Theta_[i], nvis, km_[i], airMass_);
	    m5single_.push_back(m5);
	    }

    // systematic error
    sigmaSys_ = 0.0025;
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
                                                                

