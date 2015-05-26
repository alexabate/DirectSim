#include "sedfilter.h"

namespace SOPHYA {

/*******************************************************************************
*                                                                              *
*                                SED CLASSES                                   *
*                                                                              *
*******************************************************************************/


//******* SED ****************************************************************//

// Main Constructor
SED::SED(string& fname, double xmin, double xmax, int npt)
: sed2_(sed2init_)
{

	// read the SED
	readSED(fname, xmin, xmax, npt);

	// set some defaults
	a_ = 1.;    // fraction of SED is 100% of one read in
	EBmV_ = 0.; // no reddening applied yet
	RvCard_=3.5;
	law_=0;

    // flags: SED not interpolated or reddened
	isInterp_ = false;
	isRedden_ = false;
	_test_=2.;
};


// Copy Constructor
SED::SED(SED const& a)
:  sed_(a.sed_)
{
	Set(a);
	_test_=3.;

};


// Copy SED method: this does a DEEP copy
SED& SED::Set(const SED& a)
{
    sed2_=a.sed2_;// missing this line cost me a whole live long day!
	isRead_ = a.isRead_;
    isRedden_=a.isRedden_;
    isInterp_=a.isInterp_;
    a_=a.a_;
    EBmV_=a.EBmV_;
    RvCard_=a.RvCard_;
    law_=a.law_;
	
};


//******* SED methods **************************************************************************************//

void SED::readSED(string& fname, double xmin, double xmax, int npt)
{

    // find location of SED file
  	string Tplace;
	char * pt=getenv("SEDLOC");
	if (pt==NULL)
		throw ParmError("ERROR SED LOCATION ENVIRONMENT VARIABLE -SEDLOC- NOT DEFINED");
	else
		{
		Tplace=pt;
		cout <<"     Location of template file is "<< Tplace <<endl;
		}
	string filename = Tplace + fname;

    // Set interpolation to zero outside of x-range in file 
    bool isSetZero = false; // definitely set to false!
    int nComments = 0; // no comments at top of file (could add as argument to generalise code)

    //cout << xmin <<"  "<< xmax <<" "<< npt << endl;
	sed_.ReadXYFromFile(filename, xmin, xmax, npt, nComments, isSetZero);
	isRead_=true;

};


// Main return function
double SED::returnFlux(double lambda) const
{

    if (isInterp_&&!isRedden_)// if SED is to be interpolated
        return interpSED(lambda);
    else if(isRedden_&&!isInterp_)// if SED is to be reddened
        return addReddening(lambda);
    else if (isRedden_&&isInterp_)// if SED is to be interpolated AND reddened
        return interpAddReddening(lambda);
    else // if SED will be au natural
        return sed_(lambda);

};


double SED::addReddening(double lambda) const
{
    if (!isRedden_)
        throw ParmError("ERROR! Cannot redden spectrum");
        
    if (EBmV_<0.00001)
        return sed_(lambda);
        
    Reddening red;
	double k;
	if (law_<1)
		k = red.Cardelli(lambda, RvCard_);
	if (law_>0)
		k = red.Calzetti(lambda);
	
	return  (sed_(lambda)*pow(10,-0.4*k*EBmV_));
};


double SED::interpAddReddening(double lambda) const
{

    if (!isRedden_)
        throw ParmError("ERROR! Cannot redden spectrum");
    if (!isInterp_)
        throw ParmError("ERROR! Cannot interpolate spectrum");
        
    if (EBmV_<0.00001)
        return interpSED(lambda);
    
    Reddening red;
	double k;
	if (law_<1)
		k = red.Cardelli(lambda,RvCard_);
	if (law_>0)
		k = red.Calzetti(lambda);
	
	// redden interpolated SED
	double redPart = pow(10,-0.4*k*EBmV_);
	double interpPart = interpSED(lambda);
	return  (interpPart*redPart);

};


double SED::interpSED(double lambda) const
{

    if (!isInterp_)
        throw ParmError("ERROR! Cannot interpolate spectrum");
        
    double b = 1. - a_;
        
    // sed2_ is a POINTER to a SED object
    double sed1=a_*(sed_(lambda));
    double sed2=b*(sed2_->returnFlux(lambda));
    return (sed1+sed2);
};


//******* ReadSedList ********************************************************//


ReadSedList::ReadSedList(string sedFile, int prt)
: prt_(prt) {

    // Initialize these to zero to start
    ntot_ = nsed_ = 0;

	// first get location of SED files
	sedDir_ = getSedDirEnviromentVar();
	sedFileFullPath_ = sedDir_ + sedFile;
		
	// open file, read all lines and count number of SEDs within
	countSeds(sedFileFullPath_);
    cout <<"     There are "<< nsed_ <<" SEDs to read in"<<endl;
    cout << endl;
    
    // Initialize total number of SEDs to number of SEDs to be read in
    ntot_ = nsed_; 
    
    // Initialize if SEDs added by interpolation or reddening to false
    isInterp_ = isRedden_ = false;
    
    // SEDs have not been read yet
    isRead_ = false;

};

//******* ReadSedList methods ************************************************//


string ReadSedList::getSedDirEnviromentVar()
{

    string tPlace;
    char * pt=getenv("SEDLOC");
	if (pt==NULL) {
		string emsg="ERROR SED LOCATION ENVIRONMENT VARIABLE";
		emsg+=" -SEDLOC- NOT DEFINED";
		throw ParmError(emsg);
		}
	else {
		tPlace=pt;
		cout <<"     Location of SED template files is "<< tPlace <<endl;
		}
	cout << endl;
	return tPlace;

};


void ReadSedList::countSeds(string sedFileFullPath)
{
    // open up file of list of SEDs
    ifstream ifs;
    ifs.open(sedFileFullPath.c_str(),ifstream::in);
	
	// if it fails throw an error
	if (ifs.fail()) {
		string emsg="error: Unable to open file ";
		emsg += sedFileFullPath;
		throw ParmError(emsg);
		}
		
    // otherise count number of lines in file
    string line;
    nsed_ = 0;
   	while ( ifs.good() ) {
   		getline(ifs,line);
	    nsed_++;
   		}
   	nsed_ -= 1; // not sure why this is here?
   	ifs.close();
};


void ReadSedList::readSeds(double lmin, double lmax, int npt)
{
    if (isRead_)
        throw ParmError("ERROR! SEDs have already been read");
    
    // open up file of list of SEDs
    ifstream ifs;
	ifs.open(sedFileFullPath_.c_str(),ifstream::in);
	
	// read in SED file names
	string fileNames[nsed_];
	string line;
	if (prt_ > 0)
	    cout <<"     File contains the following SEDs:"<<endl;
	for (int i=0; i<nsed_; i++)	{
		getline(ifs,line);
		sedFiles_.push_back(line);
		fileNames[i] = line; //sedDir_ + line; becase SED class can find location!
		if (prt_ > 0)
		    cout <<"     "<< fileNames[i] <<endl;
		}
		
	// read in all SEDs
	if (prt_ > 0)
	    cout <<"     Reading in all the SED files ... "<<endl;

	for (int i=0; i<nsed_; i++) {
		sedArray_.push_back(new SED()); // assigned memory for SED pointer
		sedArray_[i]->readSED(fileNames[i], lmin, lmax, npt);
		reds_.push_back(0.);
		}
		
    // close file of list of SEDs
    ifs.close();
    
    isRead_ = true; // SEDs have now been read
};

void ReadSedList::interpSeds(int nInterp)
{

    if (isRedden_)
        throw ParmError("ERROR! Make sure interpolation is done first!");

    // If we have nsed's we have nsed-1 spaces inbetween seds
    // Therefore number of SEDs being added is nInterp*(nsed-1)

    if (prt_>0) {
	    cout <<"     Before interpolation, size of SED array is ";
	    cout << sedArray_.size()<<endl;
	    }
	    
	int iSED = ntot_; // starting index of next SED
    for (int i=0; i<(nsed_-1); i++) {
    
		double a=0;

		// loop over number of interps to do
		for (int j=0; j<nInterp; j++) {
		
		    // need to find values of a and b depending on what nInterp is
		    a+=1./(nInterp+1);
		    double b=1-a;
		
		    // push new SED to end of array
		    sedArray_.push_back( new SED(*(sedArray_[i])) );// copies SED in i
		    sedArray_[iSED]->doInterp(sedArray_[i+1], a);
		    reds_.push_back(0.); // hmm this could be bad if interpolated two reddened SEDs, BUT! interp
		                         // must be done first!
		                        
            iSED++;
		    }
		 
		}
		
	int nAdd = nInterp*(nsed_-1);
    ntot_+=nAdd;
    
    if (prt_>0){
	    cout <<"     Size of SED array is now "<< sedArray_.size() <<endl;
	    cout <<"     Total number of SEDs is "<< ntot_ <<endl;
	    cout << endl; 
        }
        
     isInterp_ = true; 
};


void ReadSedList::reorderSEDs()
{
    if (isRedden_)
        throw ParmError("ERROR! Cannot reorder SEDs if reddening has been done!");
    if (!isInterp_)
        throw ParmError("ERROR! No need to reorder if not interpolated!");

    vector<SED*> tmpsedArray;
    
    // If we have nsed's we have nsed-1 spaces inbetween seds
    // Therefore number of SEDs being added is: nInterped = nInterp*(nsed-1)
    // nInterp is the number of interpolations done between each original SED
    
    
    int nInterpd = ntot_ - nsed_;
    int nInterp = nInterpd/(nsed_ - 1);
    
    //int iBeginInt = nsed_; // index of sedArray_ where the interpolated SEDS begin
    
    for (int i=0; i<(nsed_-1); i++) {
    
        // push back original SED from file first
        tmpsedArray.push_back(sedArray_[i]);
        
        // then get the nInterp SEDs between that SED and the one at i+1
        int ii = i*nInterp;
         
        for (int j=0; j<nInterp; j++) {
            int ik = nsed_ + ii + j;
            tmpsedArray.push_back(sedArray_[ik]); 
            }
        }
    
    // push back last SED in original list
    tmpsedArray.push_back(sedArray_[nsed_-1]); 
    
    int tmp = tmpsedArray.size();
    if ( tmp!=ntot_ )
        throw ParmError("ERROR! Wrong number of SEDs!");
            
    sedArray_.clear();
    sedArray_ = tmpsedArray;
};


//void ReadSedList::reddenSeds(int nStepRed, double redMax)
vector<double> ReadSedList::reddenSeds(int nPerSED, int method, int maxidEl, double redMaxEl, double redMaxOther)
{
// This is not properly implemented, the value of redMax and law will change
// depending on the SED involved
// Will need to define another array that encodes these properties e.g. an
// integer array ..
// broadtype_[i] = 0,1,2
// where 0 == elliptical type: redMax cannot be >0.1, law=0 (Cardelli)
//       1 == late type: no limit on redMax, law=0 (Cardelli)
//       1 == starburst type: no limit on redMax, law=1 (Calzetti)
    
    // Check method type makes sense
    if (method<0 || method>1) {
        stringstream ss;
        ss << method;
        string emsg = "ERROR! Reddening method = " + ss.str() + " not understood";
        throw ParmError(emsg);
        }
          
    // Check maximum index of an El makes sense
    if (maxidEl<0 || maxidEl>(sedArray_.size()-1)) {
        stringstream ss;
        ss << maxidEl;
        string emsg = "ERROR! Max index of an elliptical galaxy, " + ss.str() + ", outside array range of SEDs";
        throw ParmError(emsg);
        }
        
    if (prt_>0) {
	    cout <<"     Before reddening, size of SED array is ";
	    cout << sedArray_.size() <<endl;
	    }
	    
	// get reddening values, same for all SEDs except elliptical
	RandomGenerator rg;
	vector<double> EBV_El;    // five reddening values to apply to El's (or first maxidEl+1 SEDs)
	vector<double> EBV_Other; // five reddening values to apply to any other SED
	
	// if uniform distribution of reddening
	if (method==0) {
	
	    double EBmV;
	    for (int j=0; j<nPerSED; j++) {
            
            double r = rg.Flat01();
	        
	        EBmV = redMaxEl*rg.Flat01();
	        EBV_El.push_back(EBmV);
	        
	        EBmV = redMaxOther*rg.Flat01();
	        EBV_Other.push_back(EBmV);
	        
	        }
	    }
	// if exponential distribution of reddening
	else if (method==1) {
	
	    // not yet implemented
	    throw ParmError("ERROR! exponential reddening not yet implemented");
	    
	    }
	    

	// sort into ascending order of redness
	sort(EBV_El.begin(), EBV_El.end());
	sort(EBV_Other.begin(), EBV_Other.end());
	
	// store in a vector El then others.
	vector<double> redvals;
	for (int j=0; j<nPerSED; j++)
	    redvals.push_back(EBV_El[j]);
	for (int j=0; j<nPerSED; j++)
	    redvals.push_back(EBV_Other[j]);
	   

    // loop over each original SED
	for (int i=0; i<ntot_; i++) {
	    
	    cout <<"     On SED "<< i+1 <<" of "<< ntot_ <<endl;
		    
		// loop over number of times reddening to be done
		for (int j=0; j<nPerSED; j++) {
		    
		    // default expect galaxy is Elliptical:
		    int law = 0; // Cardelli
		    double EBmV = EBV_El[j]; // Elliptical reddening value

		    // if not Elliptical template, change law and reddening value
		    if (i>maxidEl) {
		        law = 1;
		        EBmV = EBV_Other[j];
		        }
		        
		    cout <<"     On redden "<< j+1 <<" of "<< nPerSED <<", E(B-V)="<< EBmV << endl;
		    
		    // push new SED to end of array
		    sedArray_.push_back(new SED(*(sedArray_[i])));
		    reds_.push_back(EBmV); // actual order of reddening values applied
		    
		    int lastIndex = sedArray_.size()-1; // because size just grew by 1
		    sedArray_[lastIndex]->doRedden(EBmV, law);
		    }
		    
	    }
	
    ntot_=(ntot_)*nPerSED+ntot_; // reddening all seds nStepRead times

    if (prt_>0){
	    cout <<"     Size of SED array is now "<< sedArray_.size() <<endl;
	    cout <<"     Total number of SEDs is "<< ntot_ <<endl; 
        }
        
    isRedden_ = true;

    return redvals;
};


void ReadSedList::writeSpectra(string outFile, double lmin, double lmax, int nl) {
    
    ifstream ifs;
	ofstream outp;
	
	cout <<"     Number of spectra to write = "<< sedArray_.size() <<endl;
	
	// wavelength spacing
    double dl = (lmax-lmin)/(nl-1);
    
	ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail()) {
	
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		
		for (int j=0; j<nl; j++) {
		
			double lam = lmin + j*dl;
			outp << lam <<"  ";
			
			for (int i=0; i<sedArray_.size(); i++) {
			
				double val = sedArray_[i]->returnFlux(lam);
				outp << val <<"  ";
				
				}
			outp << endl;
			
			// also print to screen for first 100 lambdas
			if ( prt_ & (j<100) ) {
				cout << lam <<"  ";
				for (int k=0; k<ntot_; k++)
					cout<< sedArray_[k]->returnFlux(lam) <<" ";
				cout << endl;
				}
				
			}
		outp.close();
		}
	else
		cout <<"     ERROR! file "<< outFile <<" exists"<<endl;
	cout << endl;


    // now write accompanying reddening value file
    string redfile;
    
    // make filename
    string s = outFile;
    
    // remove extension
    string delimiter = ".";
    size_t pos = s.find(delimiter);
    string token = s.substr(0, pos);
    
    // add new extension
    redfile = token + "_redvals.txt";

    
    ifs.open(redfile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail()) {
	
		ifs.clear(ios::failbit);
		outp.open(redfile.c_str(),ofstream::out);
		
		// WRITE
		outp <<"# SEDs in file "<< outFile << endl;
		outp <<"# column 1: column number of SED in above file"<<endl;
		outp <<"# column 2: E(B-V) value that was applied to that SED"<<endl;
		for (int i=0; i<ntot_; i++) 
		    outp << i+1 <<"  "<< reds_[i] << endl;

        outp.close();
		}
	else
		cout <<"     ERROR! file "<< redfile <<" exists"<<endl;
	cout << endl;
};


/*******************************************************************************
*                                                                              *
*                              FILTER CLASSES                                  *
*                                                                              *
*******************************************************************************/



//******* Filter methods *****************************************************//

void Filter::readFilter(string& fname, double lmin, double lmax, int npt, int nComments, bool zero_outside)
{

    string Fplace;
	char * pf=getenv("FILTLOC");
	if (pf==NULL) {
	    string emsg = "ERROR FILTER LOCATION ENVIRONMENT VARIABLE -FILTLOC- ";
	    emsg += "NOT DEFINED";
		throw ParmError(emsg);
        }
	else {
		Fplace=pf;
		cout <<"    Location of filter files are "<< Fplace <<endl;
		}

//	string filename = Fplace + fname;
      string filename = fname;

	// zero_outside = true: then interpolation sets stuff to zero outside xmin,xmax
	// definitely want this to be true for the filter
    ReadXYFromFile(filename, lmin, lmax, npt, nComments, zero_outside);

};

//******* ReadFilterList *****************************************************//

ReadFilterList::ReadFilterList(string filterFile, int prt)
{
    // Initialize these to zero to start
    ntot_ = 0;

	// first get location of filter files
	filterDir_ = getFilterDirEnviromentVar();
	filterFileFullPath_ = filterDir_ + filterFile;
		
	// open file, read all lines and count number of filters within
	countFilters(filterFileFullPath_);
    cout <<"     There are "<< ntot_ <<" filters to read in"<<endl;
    cout <<endl;
    
};
    
//******* ReadFilterList methods *********************************************//

string ReadFilterList::getFilterDirEnviromentVar()
{

    string tPlace;
    char * pt=getenv("FILTLOC");
	if (pt==NULL) {
		string emsg="ERROR FILTER LOCATION ENVIRONMENT VARIABLE";
		emsg+=" -FILTLOC- NOT DEFINED";
		throw ParmError(emsg);
		}
	else {
		tPlace=pt;
		cout <<"     Location of Filter template files is "<< tPlace <<endl;
		}
	cout << endl;
	return tPlace;

};
    

void ReadFilterList::countFilters(string filterFileFullPath)
{
    ifstream ifs;
    ifs.open(filterFileFullPath.c_str(),ifstream::in);
	if (ifs.fail())
		{
		string emsg="error: Unable to open file ";
		emsg+=filterFileFullPath;
		throw ParmError(emsg);
		}
    string line;
    ntot_=0;
   	while ( ifs.good() )
    	{
   		getline(ifs,line);
	    ntot_++;
   		}
   	ntot_-=1;
   	ifs.close();

};
    

void ReadFilterList::readFilters(double lmin, double lmax, int npt)
{
    ifstream ifs;
    
    // read in all the filter filenames
	ifs.open(filterFileFullPath_.c_str(),ifstream::in);
	string fileNames[ntot_];
	string line;
	if (prt_ > 0)
	    cout <<"     File contains the following filters:"<<endl;
	for (int i=0; i<ntot_; i++) {
		getline(ifs,line);
		fileNames[i]=filterDir_+line;
		if (prt_ > 0)
		    cout <<"     "<< fileNames[i] <<endl;
		}
		
	// read in all sed files
	if (prt_ > 0)
	    cout <<"     Reading in all the filter files ... "<<endl;

	for (int i=0; i<ntot_; i++) {
		filterArray_.push_back(new Filter());// assigned memory for SED pointer
		filterArray_[i]->readFilter(fileNames[i], lmin, lmax, npt); 
		}

};
    

void ReadFilterList::writeFilters(string outFile,double lmin,double lmax,int nl)
{
    ifstream ifs;
	ofstream outp;
	
	cout <<"     Filter number = "<<filterArray_.size()<<endl;
	
    double dl=(lmax-lmin)/(nl-1);
	ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int j=0; j<nl; j++)
			{
			double lam=lmin+j*dl;
			//cout <<"     Writing wavelength "<<j+1<<endl;
			outp <<lam<<"  ";
			for (int i=0; i<ntot_; i++)
				{
				//cout <<"     Writing spectrum "<<i+1<<endl;
				double val=filterArray_[i]->operator()(lam);
				outp << val <<"  ";
				}
			outp <<endl;
			if ( prt_ & (j<100) ) {
				cout <<lam<<"  ";
				for (int k=0; k<ntot_; k++)
					cout<<filterArray_[k]->operator()(lam)<<" ";
				cout <<endl;
				}
			}
		outp.close();
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;
	cout << endl;


};



/*******************************************************************************
*                                                                              *
*                              PCA RELATED CLASSES                             *
*                                                                              *
*******************************************************************************/

//******* SEDCovMat methods **************************************************//
TArray<double> SEDCovMat::MakeCovMat(SED **sedarray, int nsed,
						TVector<double>& meanSED)
{

	// Initialize size of SED matrix
	TArray<double> SEDmatrix;
	int ndim=2;
	sa_size_t mydim[ndim];
	mydim[0]=nl_; mydim[1]=nsed;
	SEDmatrix.SetSize(ndim, mydim);
	
	// Vector of mean SEDs
	meanSED.SetSize(nsed);
	meanSED=0; // set vector elements to equal zero
	
	// loop over wavelengths and SEDs and fill matrix, and calculate mean
	// SED values
	for (int i=0; i<nl_; i++)
		{
		double lam=lmin_+i*dl_;
		for (int j=0; j<nsed; j++)
			{
			double val=sedarray[j]->returnFlux(lam);
			//cout <<val<<endl;
			SEDmatrix(i,j)=val;
			meanSED(j)+=SEDmatrix(i,j);
			}
		}
	meanSED/=nl_;
	for (int j=0; j<nsed; j++)
		cout <<"mean of SED "<<j<<" is "<<meanSED(j)<<endl;
	
	// Subtract mean
	for (int i=0; i<nl_; i++)
		for (int j=0; j<nsed; j++)
			SEDmatrix(i,j)-=meanSED(j);
			
	/*cout <<"print SEDmatrix part to screen M(0:10,1:6) "<<endl;
	for (int i=0; i<100; i++)
		{
		cout <<" lam = "<<lmin_+i*dl_<<"  ";
		for (int j=0; j<6; j++)
			cout <<SEDmatrix(i,j)<<"  ";
		cout <<endl;
		}*/

	TArray<double> CovMat=TransposeMult(SEDmatrix);
	CovMat/=(nl_-1);
	
	return CovMat;
};


TArray<double> SEDCovMat::TransposeMult(TArray<double> SEDmatrix)
{

	//SVD matinv(SEDmatrix);
	// Do Transpose
	TArray<double> TSEDmatrix=Transpose(SEDmatrix);
	cout <<"     Size of original matrix = "<<SEDmatrix.SizeX()<<"x";
	cout <<SEDmatrix.SizeY()<<endl;
	cout <<"     Size of transposed matrix = "<<TSEDmatrix.SizeX()<<"x";
	cout <<TSEDmatrix.SizeY()<<endl;
	/*for (int i=0; i<6; i++)
		{
		for (int j=100; j<110; j++)
			cout <<TSEDmatrix(i,j)<<"  ";
		cout <<endl;
		}*/
	
	// Do multiplication
	TArray<double> TransMult=Mult(SEDmatrix,TSEDmatrix);

	return TransMult;

};

//******* TemplatePCA ********************************************************//

TemplatePCA::TemplatePCA(vector<SED*> sedArray,double lmin,double lmax,int nl)
: sedArray_(sedArray) , lmin_(lmin) , nl_(nl)
{

    dl_=(lmax-lmin_)/(nl_-1);
    nsed_=sedArray_.size();
    
    dataMatrix_.SetSize(nsed_,nl_);
    cout <<"     Set data matrix size: "<<dataMatrix_.NRows()<<"x";
    cout << dataMatrix_.NCols() << endl;
    cout << endl;

    // Normalize each spectrum such that fmean=f/(sqrt(sum(f^2)))
    // This is not the only option for normalization
    cout <<"     Find normalizations of spectra "<<endl;
    normalizeSpectra();
    cout << endl;
    
    // Find mean value of each spectrum
    cout <<"     Find means of spectra "<<endl;
    meanSpectra();
    cout << endl;
    
    // Add each normalized, mean subtracted spectrum to the data matrix
    cout <<"     Add spectra to data matrix "<<endl;
    addSpectra();
    cout << endl;
    
    // Calculate covariance matrix
    cout <<"     Calculate covariance matrix of data"<<endl;
    calculateCovarianceMatrix();
    cout << endl;

    // Do PCA
    cout <<"     Find eigenvalues and eigenvectors "<<endl;
    doPCA();
    cout << endl;
   
};

TemplatePCA::TemplatePCA(vector<SED*> sedArray, string eigVectFile, double lmin,
                                                        double lmax, int nl)
: sedArray_(sedArray) , lmin_(lmin) , nl_(nl)
{

    dl_=(lmax-lmin_)/(nl_-1);
    nsed_=sedArray_.size();
    
    dataMatrix_.SetSize(nsed_,nl_);
    cout <<"     Set data matrix size "<<dataMatrix_.NRows()<<"x";
    cout << dataMatrix_.NCols() << endl;
    cout << endl;

    // Normalize each spectrum such that fmean=f/(sqrt(sum(f^2)))
    // This is not the only option for normalization
    cout <<"     Find normalization of spectra "<<endl;
    normalizeSpectra();
    cout << endl;
    
    // Find mean value of each spectrum
    cout <<"     Find mean of spectra "<<endl;
    meanSpectra();
    cout << endl;
    
    // Add each normalized, mean subtracted spectrum to the data matrix
    cout <<"     Add spectra to data matrix "<<endl;
    addSpectra();
    cout << endl;
    
    // Don't need to calculate covariance matrix or do PCA
    
    // Read eigenvalues from a file
    cout <<"     Read eigenvectors from file "<<endl;
    eigenVectors_.ResizeTo(nl_,nl_);
    eigenVectors_=readEigenVectors(eigVectFile);
    
    cout <<"    Check eigenvectors ... "<<endl;
    cout <<"    Size of eigenvector matrix = "<< eigenVectors_.GetNrows() <<"x";
    cout << eigenVectors_.GetNcols() <<endl;
    cout <<"    Printing first 10x10 ... "<<endl;
    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++)
            cout <<eigenVectors_(i,j) <<" ";
        cout << endl;
        }
            
    
    cout << endl;
    
};

//******* TemplatePCA methods ************************************************//

void TemplatePCA::normalizeSpectra()
{
    normValues_.SetSize(nsed_);
    
    for (int i=0; i<nsed_; i++)
        {
        double sumSquared=0;
		for (int j=0; j<nl_; j++)
			{
			double lam=lmin_+j*dl_;
			double val=sedArray_[i]->returnFlux(lam);
			sumSquared+=(val*val);
			}
		normValues_(i)=sqrt(sumSquared);
		}
		
};


void TemplatePCA::meanSpectra()
{
    meanValues_.SetSize(nsed_);

    for (int i=0; i<nsed_; i++)
        {
        double sum=0;
		for (int j=0; j<nl_; j++)
			{
			double lam=lmin_+j*dl_;
			double val=sedArray_[i]->returnFlux(lam);
			sum+=(val/normValues_(i));
			}
		meanValues_(i)=sum/nl_;
		}

};


void TemplatePCA::addSpectra()
{

    for (int i=0; i<nsed_; i++) {
        //cout << "     Adding spectrum "<<i+1<<" of "<<nsed_<<endl;
        for (int j=0; j<nl_; j++) {
            
            //cout <<"     On lambda "<<j+1<<" of "<<nl_<<endl;
			double lam=lmin_+j*dl_;
			double val=sedArray_[i]->returnFlux(lam);
			dataMatrix_(i,j)=(val/normValues_(i))-meanValues_(i);
			
			}
        }

};


void TemplatePCA::calculateCovarianceMatrix()
{

    // Transpose data matrix
    //cout <<"     Transpose data matrix "<<endl;
    TMatrix<double> dataMatrixTransposed;
	dataMatrixTransposed=Transpose(dataMatrix_);
	
	// Make covariance matrix: D^T*D 
	//cout <<"     Make covariance matrix "<<endl;
	covMatrix_=Mult(dataMatrixTransposed,dataMatrix_);
	
	// Normalize covariance matrix by its trace (sum(D_ii))
	//cout <<"     Normalize covariance matrix "<<endl;
    normalizeCovarianceMatrix();
};
	        

void TemplatePCA::normalizeCovarianceMatrix()
{
    double traceCovMatrix=0;
	for (int i=0; i<nl_; i++)
	    traceCovMatrix += covMatrix_(i,i);
	cout <<"     -Trace of the covariance matrix is = "<< traceCovMatrix <<endl;
	cout <<"     -Normalize by this value"<<endl;
	covMatrix_/=traceCovMatrix;

};


void TemplatePCA::doPCA()
{
    // Copy covariance matrix to ROOT object 
    //cout <<"     Copy covariance matrix to ROOT object"<<endl;
    TMatrixD fCovarianceMatrix=copyCovMatrixToROOT();

    // Decompose into eigenvectors and eigenvalues
    //TVectorD fEigenValues(nl_);
    //TMatrixD fEigenVectors(nl_,nl_);
    eigenValues_.ResizeTo(nl_);
    eigenVectors_.ResizeTo(nl_,nl_);
    TMatrixDSym sym; sym.Use(fCovarianceMatrix.GetNrows(),
	                            fCovarianceMatrix.GetMatrixArray());
    TMatrixDSymEigen eigen(sym);
    eigenVectors_ = eigen.GetEigenVectors();
    eigenValues_ = eigen.GetEigenValues();
    
};


TMatrixD TemplatePCA::copyCovMatrixToROOT()
{
       
	TMatrixD fCovarianceMatrix(nl_,nl_);

	// convert to ROOT object
	for (int i=0; i<nl_; i++)
        for (int j=0; j<nl_; j++)
	        fCovarianceMatrix(i,j)=covMatrix_(i,j);
	        
    return fCovarianceMatrix;
    
};


void TemplatePCA::reconstructSpectra(int nEigKept)
{

    // first zero size incase have called method before with a different
    // nEigKept
    reconstructedSpectra_.ZeroSize();
	eigenvalsProjSpec_.ZeroSize();

    reconstructedSpectra_.SetSize(nl_,nsed_);
	eigenvalsProjSpec_.SetSize(nEigKept,nsed_);
    for (int i=0; i<nsed_; i++)
        {
	    int iSpectrum=i;
	    cout <<"     On spectrum "<<iSpectrum+1<<" of "<<nsed_<<endl;
	    TMatrix<double> reconstructedSpectrum;
	    TVector<double> zCut=reconstructSpectrum(nEigKept,iSpectrum,
	                                                     reconstructedSpectrum);
	    cout << "     Size of projected eigenvalues = "<<zCut.Size()<<endl;
	                                                     
	    for (int j=0; j<nl_; j++)
	        reconstructedSpectra_(j,i)=reconstructedSpectrum(j,1);
	        
	    for (int j=0; j<nEigKept; j++) // think this should be nl!!!!
	        eigenvalsProjSpec_(j,i)=zCut(j);
	        
	    cout << endl;
	    }
	cout << endl;

};

TVector<double> TemplatePCA::reconstructSpectrum(int nEigKept, int iSpectrum, 
                                        TMatrix<double>& reconstructedSpectrum)
{

    // Keep only 1st nEigKept
    //cout <<"     Only taking first "<<nEigKept<<" eigenvalues "<<endl;
    
    // Chop the eigenvector matrix at the maximum eigenvalue to keep
	TMatrixT<double> eigenVectorsCut(eigenVectors_.GetNrows(),nEigKept);
	for (int i=0; i<eigenVectors_.GetNrows(); i++)
		for (int j=0; j<nEigKept; j++)
            eigenVectorsCut(i,j)=eigenVectors_(i,j);
            
    /*cout <<"    Printing first nEigKeptxnEigKept of chopped eigenvector matrix "<<endl;
    for (int i=0; i<nEigKept; i++) {
        for (int j=0; j<nEigKept; j++)
            cout <<eigenVectorsCut(i,j) <<" ";
        cout << endl;
        }*/
            
    // Transpose the chopped eigenvector matrix 
	//cout <<"     Transpose cut eigenvector matrix"<<endl;
	cout <<"     Size of eigenvector matrix:";
	int nr=eigenVectorsCut.GetNrows();
	int nc=eigenVectorsCut.GetNcols();
	cout <<nr<<"x"<<nc<<endl;
	TMatrixT<double> eigenVectorsCutTransposed(nc,nr);
	eigenVectorsCutTransposed.Transpose(eigenVectorsCut);
	cout <<"     Size of transposed eigenvector matrix:";
	cout <<eigenVectorsCutTransposed.GetNrows()<<"x";
	cout <<eigenVectorsCutTransposed.GetNcols()<<endl;
	
	/*cout <<"    Printing first nEigKeptxnEigKept of chopped eigenvector matrix";
	cout <<" transposed "<<endl;
    for (int i=0; i<nEigKept; i++) {
        for (int j=0; j<nEigKept; j++)
            cout <<eigenVectorsCutTransposed(i,j) <<" ";
        cout << endl;
        }*/
	
	// check that the size of eigenvector matrix is ok
	if (nr!=nl_) {
	    stringstream ss1,ss2;
	    ss1<<nr; ss2<<nl_;
	    string emsg="ERROR! Length of eigenvectors "+ss1.str()+" is not equal to the";
	    emsg+=" number of wavelengths "+ss2.str();
	    throw ParmError(emsg);
	    }
	
    // Get mean subtracted data 
	TMatrixT<double> spectrumMeanSubtracted(nr,1);// needs to be a TMatrix type
	for (int i=0; i<nl_; i++){
			//double lam=lmin_+i*dl_;
			//double val=sedArray_[iSpectrum]->returnFlux(lam);
			spectrumMeanSubtracted(i,0)=dataMatrix_(iSpectrum,i);
			                            //(val/normValues_(iSpectrum))-
			                            //                      meanValues_(i);
			}

    // Matrix multiply the transposed eigenvectors and the mean subtracted spectrum
	// zCut are the eigenvalues of the spectrum projected onto the nEigKept
	// eigenvectors
	TMatrixT<double> zCut(nEigKept,1);
	zCut.Mult(eigenVectorsCutTransposed,spectrumMeanSubtracted);
	cout <<"     Result size = "<<zCut.GetNrows()<<"x"<<zCut.GetNcols()<<endl;
    
    
	// Rotate data back
	// Multiply each column of the truncated eigenvector matrix by each
	// eigenvalue of the projected spectrum
	TMatrixT<double> rot(nl_,nEigKept);
	for (int i=0; i<nl_; i++)
	    for (int j=0; j<nEigKept; j++)
	        rot(i,j)=eigenVectorsCut(i,j)*zCut(j,0);
	        
    // Sum over eigenvalues (the columns)
	TMatrixT<double> recSpec(nl_,1);
	recSpec=0;
	for (int i=0; i<nl_; i++)
	    for (int j=0; j<nEigKept; j++)
	        recSpec(i,0)+=( rot(i,j) ); //+meanValues_(i)/nEigKept );
	     
	     
	// Add back mean spectrum
	for (int i=0; i<nl_; i++)
	       recSpec(i,0) = recSpec(i,0) + meanValues_(iSpectrum);
	        
	
	// Basically done now!
	
	// Convert to Sophya objects
	TVector<double> eigenvaluesSpectrum(nEigKept);
	for (int i=0; i<nEigKept; i++)
	    eigenvaluesSpectrum(i)=zCut(i,0);
	
	reconstructedSpectrum.SetSize(nl_,2);
	for (int i=0; i<nl_; i++)
	        {
	        double lam=lmin_+i*dl_;
	        reconstructedSpectrum(i,0)=lam; // add in lambda for good measure
	        reconstructedSpectrum(i,1)=recSpec(i,0);
	        }
	        
	return eigenvaluesSpectrum;

};


TVector<double> TemplatePCA::fitSpectra()
{

    TVector<double> fitValues(nsed_);
    for (int i=0; i<nsed_; i++)
        {
        double chisq=0.;
        for (int j=0; j<nl_; j++)
			{
			double lam=lmin_+j*dl_;
			double sp1=sedArray_[i]->returnFlux(lam);
			double sp2=reconstructedSpectra_(j,i)*normValues_(i);
			
            chisq+=((sp1-sp2)*(sp1-sp2));

			}
		fitValues(i)=chisq/(double)nl_;
		}
	return fitValues;

};


TMatrix<double> TemplatePCA::getCovMatrix()
{


    int nr=covMatrix_.NRows();
    int nc=covMatrix_.NCols();
    TMatrix<double> covMatrixSophyaMatrix(nr,nc);
    for (int i=0; i<nr; i++)
			for (int j=0; j<nc; j++)
				covMatrixSophyaMatrix(i,j)=covMatrix_(i,j);
				
	return covMatrixSophyaMatrix;

};

TVector<double> TemplatePCA::getEigenValues()
{

    int nv=eigenValues_.GetNoElements();
    TVector<double> eigenValuesSophyaVector(nv);
    for (int i=0; i<nv; i++)
        eigenValuesSophyaVector(i)=eigenValues_(i);
        
    return eigenValuesSophyaVector;
};

TMatrix<double> TemplatePCA::getEigenVectors()
{

    int nr=eigenVectors_.GetNrows();
    int nc=eigenVectors_.GetNcols();

    TMatrix<double> eigenVectorsSophyaMatrix(nr,nc);
    for (int i=0; i<nr; i++)
			for (int j=0; j<nc; j++)
				eigenVectorsSophyaMatrix(i,j)=eigenVectors_(i,j);
	return eigenVectorsSophyaMatrix;

};

// WRITE TO FILES

// Data normalization values
void TemplatePCA::writeNormValues(string outFile)
{

    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<normValues_.Size(); i++)
				outp<< normValues_(i) <<endl;
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;   


};


// Data mean values
void TemplatePCA::writeMeanValues(string outFile)
{

    int nMeanValues = meanValues_.Size();
    //cout <<"     Number of mean values = "<<nMeanValues<<endl;

    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<nMeanValues; i++)
				outp<< meanValues_(i) <<endl;
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;   


};


// Covariance matrix
void TemplatePCA::writeCovMatrix(string outFile)
{

    int nr = covMatrix_.NRows();
    int nc = covMatrix_.NCols();
    //cout <<"     Size of covariance matrix = "<< nr << "x" << nc <<endl;
    
    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<nr; i++)
			{
			for (int j=0; j<nc; j++)
				outp<< covMatrix_(i,j)<<"  ";
			outp << endl; 
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;

};

// data matrix
void TemplatePCA::writeDataMatrix(string outFile)
{

    
    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<nsed_; i++) {
			for (int j=0; j<nl_; j++)
                outp<< dataMatrix_(i,j)<<"  ";
			outp << endl;
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;

};

// eigenvectors
void TemplatePCA::writeEigenVectors(string outFile)
{

    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<eigenVectors_.GetNrows(); i++)
			{
			for (int j=0; j<eigenVectors_.GetNcols(); j++)
				outp<< eigenVectors_(i,j)<<"  ";
			outp << endl;
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;

};

// eigenvalues
void TemplatePCA::writeEigenValues(string outFile)
{

    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<eigenValues_.GetNrows(); i++)
				outp<< eigenValues_(i)<<endl;
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;



};

// both eigenvectors and eigenvalues
void TemplatePCA::writeEigenValVecs(string outFile)
{

    //const TVectorD* eigenValues=principal_->GetEigenValues();
    
    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<eigenVectors_.GetNrows(); i++)
			{
			for (int j=0; j<eigenVectors_.GetNcols(); j++)
				outp<< eigenVectors_(i,j)<<"  ";
			outp <<eigenValues_(i)<< endl;
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;

};


void TemplatePCA::writeEigenValsOfProjSpec(string outFile)
{

    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<eigenvalsProjSpec_.NRows(); i++)
			{
			for (int j=0; j<nsed_; j++)
			    outp <<eigenvalsProjSpec_(i,j) << "  ";
			outp << endl;
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;
	cout << endl;


};

void TemplatePCA::writeRecSpec(string outFile)
{
    ifstream ifs;
	ofstream outp;
    ifs.open(outFile.c_str(),ifstream::in);
	ifs.close();
	if(ifs.fail())
		{
		ifs.clear(ios::failbit);
		outp.open(outFile.c_str(),ofstream::out);
		for (int i=0; i<nl_; i++)
			{
			double lam=lmin_+i*dl_;
			outp <<lam<<"  ";
			for (int j=0; j<nsed_; j++)
			    outp << reconstructedSpectra_(i,j)*normValues_(j) << "  ";
			outp << endl;
			}
		outp.close();	
		}
	else
		cout <<"     ERROR! file "<<outFile<<" exists"<<endl;
	cout << endl;

};

// Read from file
TMatrixD TemplatePCA::readEigenVectors(string inFile)
{

        ifstream ifs;
        cout <<"     Reading in eigenvectors from file "<<inFile<<endl;
	    ifs.open(inFile.c_str(), ifstream::in);
	    sa_size_t nr, nc;
	    TArray<r_4> spectrum;
        spectrum.ReadASCII(ifs,nr,nc);
        cout <<"    Number of eigenvectors = "<<nc<<endl;
        cout <<"    Length of eigenvectors = "<<nr<<endl;

        //TMatrixT<double> eigenvals(nr,nc);
        TMatrixD eigenvals(nr,nc);
        for (int j=0; j<nr; j++)
            for (int i=0; i<nc; i++)
                eigenvals(j,i)=spectrum(i,j);
                
        return eigenvals;


};

}// end namespace SOPHYA
