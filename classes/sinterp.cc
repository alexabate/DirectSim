#include "sinterp.h"

//-------------------------------------------
// --- Class SInterp1D : 
//  Simple linear 1D interpolation class 
//-------------------------------------------

// Default constructor, function is y=x
SInterp1D::SInterp1D()
  : xmin_(0.), xmax_(1.), dx_(1.), ksmx_(1), npoints_(0)
{
    setzero_=false;
    xs_.push_back(0.);
    xs_.push_back(1.); 
    ys_.push_back(0.);
    ys_.push_back(1.); 
};


// constructor with regularly spaced x
SInterp1D::SInterp1D(double xmin, double xmax, vector<double>& yreg)
  : xmin_(0.), xmax_(1.), dx_(1.), ksmx_(1), npoints_(0)
{
    setzero_=false;
    DefinePoints(xmin, xmax, yreg);
};


// constructor without regularly spaced x
SInterp1D::SInterp1D(vector<double>& xs, vector<double>& ys, double xmin, double xmax, size_t npt)
  : xmin_(0.), xmax_(1.), dx_(1.), ksmx_(1), npoints_(0)
{
    setzero_=false;
    DefinePoints(xs, ys, xmin, xmax, npt);
};


// constructor without regularly spaced x and TVector
SInterp1D::SInterp1D(TVector<r_8>& xs, TVector<r_8>& ys,  double xmin, double xmax, size_t npt)
  : xmin_(0.), xmax_(1.), dx_(1.), ksmx_(1), npoints_(0)
{
    vector<double> xs_std, ys_std;
    for (int i=0; i<xs.Size(); i++) {
        xs_std.push_back(xs(i));
        ys_std.push_back(ys(i));
        }
        
    setzero_=false;
    DefinePoints(xs_std, ys_std, xmin, xmax, npt);
};


// Interpolate y value
double SInterp1D::YInterp(double x) const
{
    if (npoints_>0) { // use regularly spaced points (computed in DefinePoints)
        
        long i = (long)((x-xmin_)/dx_); // find closest index to given x value
    
        // if index is negative (i.e. x < xmin-dx)
        if (i<0) {
            
            if(!setzero_) { // if not automatically setting zero ...
                // linearly extrapolate outside range:
				return ( yreg_[0]+(x-xmin_)*(yreg_[1]-yreg_[0])/dx_ );
                } 
            else if (setzero_) // if automatically set to zero
                return 0;
            }
            
        // if index outside yreg_ (x>xmax_)
		if (i>=(int)npoints_) {
		
            if(!setzero_) // if not automatically setting zero ...
                return ( yreg_[npoints_]+(x-xmax_)*(yreg_[npoints_]-yreg_[npoints_-1])/dx_ );
            else if (setzero_) // if automatically set to zero
                return 0;
            }

        // otherwise, if within range
        return (yreg_[i]+(x-X(i))*(yreg_[i+1]-yreg_[i])/dx_);
	
		}
  	else { // use the points xs_,ys_ directly 
		
		// if x outside x,y pair values linearly extrapolate outside range
		
		// if x<xs
        if (x<=xs_[0]) 
    	    return ( ys_[0]+(x-xs_[0])*(ys_[1]-ys_[0])/(xs_[1]-xs_[0]) );
    		
        // if x>xs
        if (x>=xs_[ksmx_]) 
    	    return ( ys_[ksmx_]+(x-xs_[ksmx_])*(ys_[ksmx_]-ys_[ksmx_-1])/
    	                                           (xs_[ksmx_]-xs_[ksmx_-1]) );
        // if x within range
        size_t k=1;
        while(x>xs_[k]) k++; {
        
            if (k>=xs_.size()) {
	            // this should not happen ...
	            string emsg = " SInterp1D::YInterp() out of range k -> BUG in code ";
                throw out_of_range(emsg);
                }
            }

        double rv = ys_[k-1]+(x-xs_[k-1])*(ys_[k]-ys_[k-1])/(xs_[k]-xs_[k-1]);
//    cout << " DBG- x=" << x << " k=" << k << " xs[k]=" << xs_[k] << " ys[k]" << ys_[k] 
//	 << " rv=" << rv << endl;
        return rv;
  		}

};


// Define interpolation table with regularly spaced y values
void SInterp1D::DefinePoints(double xmin, double xmax, vector<double>& yreg)
{
  
    if (yreg.size()<2) {
        string emsg = "SInterp1D::DefinePoints(xmin,xmax,yreg) Bad parameters yreg.size()<2 ";
        throw range_error(emsg);
        }
  
    // not sure xmin and xmax are ever checked?
    xmin_ = xmin; 
    xmax_ = xmax;
    npoints_ = yreg.size()-1; // because x max is not included
    dx_ = (xmax_-xmin_)/(double)npoints_;
    yreg_ = yreg;
};


// Define interpolation table with x,y pairs
void SInterp1D::DefinePoints(vector<double>& xs, vector<double>& ys, 
                                           double xmin, double xmax, size_t npt)
{

    // check sizes
    if ((xs.size() != ys.size())||(xs.size()<2)) {
        string emsg = "SInterp1D::DefinePoints() Bad parameters ";
        emsg+="(xs.size() != ys.size())||(xs.size()<2) ";
        throw range_error(emsg);
        }
        
    // check x is sorted
    for(size_t k=1; k<xs.size(); k++) {
	
        if (xs[k-1]>=xs[k]) { 
            string emsg =  "SInterp1D::DefinePoints() unsorted xs";
            throw range_error(emsg);
            }
        }

		
    xs_ = xs;
    ys_ = ys;
    ksmx_ = xs_.size()-1;
    npoints_ = npt; // npt could be equal to zero
	
    // check xmin,xmax,npt make sense
    if (xmin>=xmax) { // use xs range to set xmin,xmax
        xmin_ = xs_[0]; 
        xmax_ = xs_[ksmx_];
        }
    else {
        xmin_ = xmin; 
        xmax_ = xmax;
        if (xmin_<xs_[0])
            xmin_ = xs_[0];
        if (xmax_>xs_[ksmx_])
            xmax_ = xs_[ksmx_];
        }
    if (npoints_<1) {
        // does not defined regularly spaced points 
        // - interpolation done directly on points given
        dx_ = (xmax_-xmin_)/(double)(xs_.size()-1);
        return;
        }
        
    // interpolation table variables
    dx_ = (xmax_-xmin_)/(double)npoints_;
    yreg_.resize(npoints_+1); // npoints+1 y vals from regularly spaced x

    // Compute the the y values for regularly spaced x xmin <= x <= xmax 
    // and keep values in the yreg vector
    // yreg defined between input y vector's limits
    yreg_[0] = ys_[0]; // set first element
    yreg_[npoints_] = ys_[ksmx_];
    size_t k=1; // k=1 because already set k=0
    for (size_t i=0; i<npoints_; i++) {
    
        double x = X(i); // X is a function that returns xmin_ + i*dx_
        
        while (x>xs_[k]) k++; // to get to the first index where xs_[k]>x
		
        if (k>=xs_.size()) { 
            // will happen if xmax given is greater than max of xs_
            string emsg="SInterp1D::DefinePoints()  out of range k -> BUG in code ";
            throw out_of_range(emsg);
            //k = ksmx_;
            }
            
        // xs[k-1] < x < xs[k]
        // yreg_[i] lies between ys_[k-1] and ys_[k]
        
        double m_grad = (ys_[k]-ys_[k-1])/(xs_[k]-xs_[k-1]);
        
        // intercept = ys[k-1] - m_grad*xs[k-1]
        // yreg[i] = m_grad*x + intercept
        // yreg[i] = m_grad*x + ys[k-1] - m_grad*xs[k-1] = ys[k-1] + m_grad*(x-xs[k-1])
		
        yreg_[i] = ys_[k-1] + (x-xs_[k-1])*m_grad;
		//DBG cout << " DBG* i=" << i << " X(i)=" << X(i) << " yreg_[i]= ";
		//cout << yreg_[i] << " X^2= " << X(i)*X(i) 
		//DBG << " k=" << k << " xs[k]=" << xs_[k] << endl;
		//cout << " DBG* i=" << i << " X(i)=" << X(i) << " yreg_[i]= " << yreg_[i] 
		//     << " k=" << k << " xs[k]=" << xs_[k] << endl;
		}
    return;
};


// Read regularly spaced y values from a file
size_t SInterp1D::ReadYFromFile(string const& filename, double xmin, 
                                                    double xmax, int nComments)
{
    ifstream inputFile;
    inputFile.open(filename.c_str(), ifstream::in);  
    if ( !inputFile.is_open() ) {
        string emsg = "  SInterp1D::ReadYFromFile() problem opening file ";
        emsg += filename;
        throw runtime_error(emsg);
        }
	    
    for (int i=0; i<nComments; i++) {
        string line;
        getline(inputFile,line);
        }
    vector<double> xsv, ysv;

    size_t cnt=0;
    double cola;
    while(!inputFile.eof())  { 
        inputFile.clear(); 
        inputFile >> cola; 
        if ( (!inputFile.good()) || inputFile.eof()) break;
        //cout << cola<< "    "<<colb<<endl;
        ysv.push_back(cola);
        cnt++;
        }
    inputFile.close();
    cout << " SInterp1D::ReadYFromFile()/Info: " << cnt;
    cout << " Y-values read from file " << filename << endl;  
    DefinePoints(xmin, xmax, ysv);
    return cnt;
};


// Read x,y pairs from a file
size_t SInterp1D::ReadXYFromFile(string const& filename,double xmin,double xmax, 
		size_t npt, int nComments, bool setzero)
{
    //so outside xsv[0] and xsv[xsv.size()-1] is forced to be set to zero
    setzero_ = setzero; 

    ifstream inputFile;
    inputFile.open(filename.c_str(), ifstream::in);  
    if( !inputFile.is_open() ) {
        string emsg = "SInterp1D::ReadXYFromFile() problem opening file ";
        emsg += filename;
        throw runtime_error(emsg);
        }
		
    for (int i=0; i<nComments; i++) {
        string line;
        getline(inputFile,line);
        }
    vector<double> xsv, ysv;// the two columns to be read in

    size_t cnt=0; // count number of lines
	//cout << "nComments = "<< nComments << ", cnt = "<< cnt <<endl;
    double cola, colb;
	while(!inputFile.eof()) {  
        inputFile.clear();
        inputFile >> cola >> colb; // read each line into cola and colb
        //cout <<" cola = "<< cola <<", colb = "<< colb << endl;
        
        if ( (!inputFile.good()) || inputFile.eof() ) {   
            break; // read until end of file
            }

        xsv.push_back(cola);
        ysv.push_back(colb);

        cnt++;
        }
    inputFile.close();
    size_t nLines = cnt - nComments;
    //cout << " SInterp1D::ReadXYFromFile()/Info: " << nLines;
    //cout << " (x,y) pairs read from file " << filename << endl;  
    DefinePoints(xsv, ysv, xmin, xmax, npt);
	
    return nLines;
};


// print to output stream
ostream& SInterp1D::Print(ostream& os, int lev)  const
{
    os << " ---- SInterp1D::Print() XMin=" << XMin() << " XMax=" << XMax();
    os << " NPoints=" << npoints_ << endl;
    os << "  xs_.size()= " << xs_.size() << "  ys_.size()= " << ys_.size();
    os << "  yreg_.size()= " << yreg_.size() << endl;
  
    if ((lev>0)&&(xs_.size()>0)) {
        for(size_t i=0; i<xs_.size(); i++) 
            os << " xs[" << i << " ]=" << xs_[i] << " -> ys[" << i << "]=" << ys_[i] << endl;
        }
         
    if ((lev>0)&&(yreg_.size()>0)) {
        for(size_t i=0; i<yreg_.size(); i++) { 
            os << " Regularly Spaced X(" << i << " )=" << X(i) << " -> yreg_[";
            os << i << "]=" << yreg_[i] << endl;
            }
        } 
    os << " ----------------------------------------------------------" << endl;
};


/******* SInterp2D methods ****************************************************/

void SInterp2D::rangeChecks(vector<double>& xa, vector<double>& xb, TArray<double>& y)
{
    int na = xa.size();
    int nb = xb.size();
    int nya = y.SizeX();
    int nyb = y.SizeY();
    
    // Check dimensions match up 
    if ((na != nya) || (nb != nyb)) {
		string emsg = "ERROR! variable vectors not the same size as array dimensions! ";
		throw range_error(emsg);
		}
	
	// Check that the array has more elements than 1!
	if ( nya*nyb < 2) {
	    string emsg = "ERROR! array has size=1";
	    throw range_error(emsg);
		}
		
    // Check variables are sorted
	for (int k=1; k<na; k++) {
		if ( xa[k-1]>=xa[k] ) { 
			string emsg =  "ERROR! xa is unsorted";
			throw range_error(emsg);
			}
		}
    for (int k=1; k<nb; k++) {
		if ( xb[k-1]>=xb[k] ) { 
			string emsg =  "ERROR! xb is unsorted";
			throw range_error(emsg);
			}
		}

};


// 1 == a, 2 == b
double SInterp2D::biLinear(double x1, double x2) 
{
    // x1 corresponds to rows direction of y (y-axis direction)
    // x2 corresponds to columns direction of y (x-axis direction)

    // first find closest xa element to x1[x2], where xa<x1[xb<x2]
    int ia = findxaElement(x1); // x1 y-axis == rows dir == 1st dim
    int ib = findxbElement(x2); // x2 x-axis == cols dir == 2nd dim
    
    double t = (x2 - xb_[ib])/(xb_[ib+1] - xb_[ib]);
    double u = (x1 - xa_[ia])/(xa_[ia+1] - xa_[ia]);
    
    // values of the array at the grid square within which the point (x1,x2) falls
    double y1 = y_(ia+1,ib);
    double y2 = y_(ia+1,ib+1);
    double y3 = y_(ia,ib+1);
    double y4 = y_(ia,ib);
    
    //cout << "Corners: " << y1 <<"  "<< y2 <<"  "<< y3 <<"  "<< y4 <<endl;
    
    // Bilinear interpolation
    double y = (1 - t)*(1 - u)*y1 + t*(1 - u)*y2 + t*u*y3 + (1 - t)*u*y4;
    
    return y;
};


double SInterp2D::biLinearAccurate(double x1, double x2) 
// Given arrays xa_[1..na_] and xb_[1..nb_] of independent variables, and a submatrix of function
// values y_[1..na_][1..nb_], tabulated at the grid points defined by xa and xb; and given values
// x1 and x2 of the independent variables; this routine returns an interpolated function value y,
// and an accuracy indication dy (based only on the interpolation in the x1 direction, however).
{

    vector<double> ymtmp;
    
    for (int j=0; j<nb_; j++) { // Loop over columns
        vector<double> ya = getColumn(j);
        SInterp1D interpa(xa_,ya);
        //Interpolate answer into temporary storage
        double y1=interpa(x1);
        ymtmp.push_back(y1);
        } 
    SInterp1D interpa(xb_,ymtmp); // Do the final interpolation.
    double y = interpa(x2);
    
    return y;
};
