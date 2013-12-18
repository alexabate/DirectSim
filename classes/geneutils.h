/**
 * @file  geneutils.h
 * @brief Contains a series of useful, generic methods
 *
 * @note Better interpolations class is in sinterp.h
 *
 * @author Alex Abate, Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2008
 * @date 2008
 *
 */
#ifndef GENEUTILS_SEEN
#define GENEUTILS_SEEN

#include "machdefs.h"
#include <math.h>
// Sophya update v2.3 June 2013 replaces genericfunc with classfunc
//#include "genericfunc.h"
#include "classfunc.h"
#include "histos.h"
#include "tvector.h"
#include "cspline.h"

#include <vector>
#include <algorithm>


namespace SOPHYA {


/** @class
  * InterpFunc class
  * 
  * Simple linear 1D interpolation class 
  * @note Don't use this! use SInterp1D
  *
  */
class InterpFunc {
public:
    /** Constructor */
    InterpFunc(double xmin, double xmax, vector<double>& y);
    virtual ~InterpFunc(void) { }

    double XMin(void) {return _xmin;}
    double XMax(void) {return _xmax;}
    inline double X(long i) {return _xmin + i*_dx;}

    /** Return the nearest element to f giving y=f(x) */
    inline double operator()(double x) {
        x -= _xmin;
        long i = long(x/_dx+0.5);  // to take the nearest "i"
        if(i<0) i=0; else if(i>=_nm1) i=_nm1-1;
        return _y[i];
        }

    /** idem operator(double) and return
        ok==0 if the value is found, 1 if x<xmin, 2 if x>xmax */
    inline double operator()(double x,unsigned short& ok)
        { ok=0; if(x<_xmin) ok=1; else if(x>_xmax) ok=2; return (*this)(x); }

    /** Return the linear interpolation of f giving y=f(x)
        ok==0 if the value is found, 1 if x<xmin, 2 if x>xmax */
    double Linear(double x,unsigned short& ok);

    /**  Return the parabolic interpolation of f giving y=f(x)
        ok==0 if the value is found, 1 if x<xmin, 2 if x>xmax */
    double Parab(double x,unsigned short& ok);

protected:
    double _xmin;   /**< min of interpolation table                           */
    double _xmax;   /**< max of interpolation table                           */
    double _dx;     /**< x spacing of original vector used to create table    */
    long _nm1;      /**< n-1                                                  */
    vector<double>& _y; /**< y values of original vector used to create table */
};


/** @class
  * InverseFunc class
  *
  * Invert 1D function. The function MUST be growing monotonically and be 
  * regularly spaced
  */
class InverseFunc {
public:
    /** Constructor, invert function @param y(@param x). The function @param y
        must be monotonically increasing, where x(i)<x(i+1) and y(i)<y(i+1)
        @param x    x values of function 
        @param y    y values of function                                      */
    InverseFunc(vector<double>& x, vector<double>& y);

    /** Destructor */
    virtual ~InverseFunc(void){ };

    /** computation of the inverted function via linear interpolation 
        @param n        number of points of inverted function
        @param xfcty    inverted function                                     */
    int ComputeLinear(long n,vector<double>& xfcty);

    /** computation of the inverted function via parabolic interpolation 
        @param n        number of points of inverted function
        @param xfcty    inverted function                                     */
    int ComputeParab(long n,vector<double>& xfcty);

    double YMin(void) {return _ymin;}
    double YMax(void) {return _ymax;}
protected:

    /** find value position in table of y values 
        @param x        y(?) value to find in _y (check)
        @param klo      index just below y value
        @param khi      index just above y value                             */
    inline void find_in_y(double x, long& klo, long& khi) {
        long k;
        klo=0, khi=_y.size()-1;
        while (khi-klo > 1) {
            k = (khi+klo) >> 1; // find index ~exactly in between klo and khi
	        if (_y[k] > x) khi=k; else klo=k;
            }
        }

    double _ymin;       /**< minimum value of function to invert              */
    double _ymax;       /**< maximum value of function to invert              */
    vector<double>& _x; /**< x values of function to invert                   */
    vector<double>& _y; /**< y values of function to invert                   */
};


/** Interp table */
double InterpTab(double x0, vector<double> const& X, vector<double> const& Y, unsigned short typint=0);

//-------- Convert functions ---------------------------------------------------

/** Convert function to a histogram 
    @param func     function to convert to histogram
    @param h        histogram to fill
    @param logaxex  if true:  abscisses "x" of the bins are filled with f(10^x)
                    if false: abscisses "x" of the bins are filled with f(x)  */
int FuncToHisto(ClassFunc1D& func, Histo& h, bool logaxex=false);

/** Convert function to a vector 
    @param func     function to convert to histogram
    @param h        vector to fill
    @param xmin     minimum x in vector
    @param xmax     maximum x in vector
    @param logaxex  if true:  abscisses "x" of the bins are filled with f(10^x)
                    if false: abscisses "x" of the bins are filled with f(x)  */
int FuncToVec(ClassFunc1D& func, TVector<r_8>& h, double xmin, double xmax, bool logaxex=false);

//------- Angle area functions -------------------------------------------------

/** Return the solid angle of a "rectangle" @param dtheta x @param dphi centered
    on @param theta_0
    @param dtheta   delta-angle in theta [0,pi] coord direction
    @param dphi     delta-angle in phi [0,2pi] coord direction
    @param theta0   theta coord value rectangle is centered on                */
double AngSol(double dtheta, double dphi, double theta0=M_PI/2.);

/** Return solid angle of a cap of opening angle @param dtheta
    @param dtheta   opening angle                                             */
double AngSol(double dtheta);

/** Return the opening angle of a spherical cap with solid angle @param angsol
    @param angsol   solid angle                                               */
double FrAngSol(double angsol);

//------- sin functions --------------------------------------------------------

/** Calculate sin(x)/x. If x^2 < 1.7e-4 or @param app=true returns: 
    1 - x^2/6 * ( 1 - x^2/20 * (1 - x^2/42) )
    @param x    x
    @param app  if true always use approximation                              */
double SinXsX(double x, bool app=false);

/** Calculate (sin(x)/x)^2. If x^2 < 6.8e-5 or @param app=true returns:
    1. - x^2/3 * ( 1 - 2*x^2/15 * (1 - x^2/14) )
    @param x    x
    @param app  if true always use approximation                              */
double SinXsX_Sqr(double x, bool app=false);

/** Calculate sin(Nx)/sin(x). Approximation used if sin(x)^2 < 3.5e-6/N^2
    @param x    x
    @param N    N
    @param app  if true always use approximation                              */
double SinNXsX(double x,unsigned long N,bool app=false);

/** Calculate (sin(Nx)/sin(x))^2. Approximation used if sin(x)^2 ~ 1.5e-6/N^2
    @param x    x
    @param N    N
    @param app  if true always use approximation   */  
double SinNXsX_Sqr(double x, unsigned long N, bool app=false);

//------- Integration functions ------------------------------------------------

/** Integrate function with ONE variable 
    @param func     function to integrate
    @param xmin     lower integral bound
    @param xmax     upper integral bound
    @param perc     integral precision
    @param dxinc    search increment
    @param dxmax    max possible increment
    @param glorder  Gauss-Legendre order                                      */
double IntegrateFunc(ClassFunc1D const& func, double xmin, double xmax,
    double perc=0.1, double dxinc=-1., double dxmax=-1., unsigned short glorder=4);

/** Integrate function with ONE variable and ONE parameter/ or one variable set to a constant
    parameter/constant variable must be the SECOND argument to func.          
    @param func     function to integrate
    @param par      parameter to set in function
    @param xmin     lower integral bound
    @param xmax     upper integral bound
    @param perc     integral precision
    @param dxinc    search increment
    @param dxmax    max possible increment
    @param glorder  Gauss-Legendre order                                      */
double IntegrateFunc(ClassFunc2D const&  func, double par, double xmin, double xmax,
    double perc=0.1, double dxinc=-1., double dxmax=-1., unsigned short glorder=4);

/** Integrate function with ONE variable in log space
    @param func     function to integrate
    @param lxmin    lower integral bound (log)
    @param lxmax    upper integral bound (log)
    @param perc     integral precision
    @param dxinc    search increment
    @param dxmax    max possible increment
    @param glorder  Gauss-Legendre order                                      */
double IntegrateFuncLog(ClassFunc1D const& func, double lxmin, double lxmax,
    double perc=0.1, double dlxinc=-1., double dlxmax=-1., unsigned short glorder=4);

/** Integrate function w/trapezium method in log space
    @param func     function to integrate 
    @param xmin     lower integral bound
    @param xmax     upper integral bound
    @param nx       integration precision (number of trapesiums)              */
double IntegrateFuncTrapLogSpace(ClassFunc1D const& func, double xmin, double xmax, int nx);

/** Compute Gauss-Legendre                               
    @param glorder  Gauss-Legendre order
    @param x        x values where coefficients are calculated
    @param w        value of associated coefficients
    @param x1       lower interval limit
    @param x2       upper interval limit                                      */
void Compute_GaussLeg(unsigned short glorder, vector<double>& x, vector<double>& w,
    double x1=0., double x2=1.);



//------- Some special functions -----------------------------------------------

/** Log of \f$ \Gamma $/f function
    @param x    x                                                             */
double LogGammaFunc(double x);

/** Error function
    @param x    x                                                             */
double ERF(double x);

/** Incomplete \f$ \Gamma $/f function. Uses either GCF or GSER method depending
    on values of x and a
    @param a    a
    @param x    x                                                             */
double GAMMP(double a, double x);

/** Returns incomplete \f$ \Gamma $/f function evaluated by its continued
    fraction representation, also returns the log of the \f$ \Gamma $/f function
    @param a
    @param x
    @param gammcf   incomplete \f$ \Gamma $/f function
    @param gln      log of \f$ \Gamma $/f function                            */
double GCF(double a, double x, double& gammcf, double& gln);

/** Returns the incomplete \f$ \Gamma $/f function evaluated by its series 
    representation. Also returns log of \f$ \Gamma $/f function
    @param a
    @param x
    @param gamser   incomplete \f$ \Gamma $/f function
    @param gln      log of \f$ \Gamma $/f function                            */
double GSER(double a, double x, double& gamser, double& gln);
 
//------ Root finder -----------------------------------------------------------

/** RootFinder Using Ridder's Method From Numerical Recipes In Fortran.
    Given the function @param func which is a function of x and some 
    parameters, calculate the root of the function known to lie between @param x1
    and @param x2. The value of the root  is refined until its accuracy is 
    +/- @param xacc.
    Translated into C++ by AA 
    @param func     function to find root of
    @param x1       root of function is at x>x1
    @param x2       root of function is at x<x2
    @param xacc     accuracy of root                                          */
double  RootRidder(ClassFunc1D const& func, double x1, double x2, double xacc=1e-4);

//------ NaN/Inf checker ------------------------------------------------------//

/** Check if x is not-a-number 
    @param x    value to check if NaN                                         */
inline int my_isnan(double x)
	{ return x != x;}

/** Check if x is infinite 
    @param x    value to check if infinite                                    */
inline int my_isinf(double x) { 
	if ((x == x) && ((x - x) != 0.0)) 
		return (x < 0.0 ? -1 : 1);
	else 
		return 0; }

//double MAX(double a,double b)
//	{ return (a+b + abs(a-b))/2; };
	
//	double MIN(double a,double b)
//	{ return (a+b - abs(a-b))/2; };
	
//int	SIGN(double x)
//	{ return (int)x/abs(x); };

//------ Useful vector stuff -------------------------------------------------//

/** Find maximum value within a vector
    @param values is a vector of doubles                                      */
double findMaximum(vector<double> values);

/** Find maximum value within a vector
    @param values is a vector of doubles
    @param iMax is the index of the maximum value                             */
double findMaximum(vector<double> values, int& iMax);

/** Find position of maximum value within a 2-dimensional array
    @param array    TArray of doubles
    @param iRow     index (zero-indexed) of row containing max value
    @param jCol     index (zero-indexed of column containing max value        */
double findMaximumPosition(TArray<double> array, int& iRow, int& jCol);
	
/** Find minimum value within a vector
    @param values is a vector of doubles                                      */
double findMinimum(vector<double> values);

/** Find position of minimum value within a 2-dimensional array
    @param array    TArray of doubles
    @param iRow     index (zero-indexed) of row containing min value
    @param jCol     index (zero-indexed of column containing min value        */
double findMinimumPosition(TArray<double> array, int& iRow, int& jCol);

/** Find position of minimum value within a 1-dimensional array
    @param array    vector of doubles
    @param iRow     index (zero-indexed) of element containing min value      */
double findMinimumPosition(vector<double> array, int& iElement);

/** Find closest value within a vector to some given value
    @param values is a vector of doubles
    @param val is a double                                                    */
int findClosestElement(vector<double> values,double val);

/** Check if a vector is sorted in ascending order                            */
bool sortCheck(vector<double> values);

/** Check non-decreasing */
bool checkNonDecreasing(vector<double> xs);

/** Remove all repeats from a vector */
vector<int> uniqueVector(vector<int> vect,vector<int>& ids);


/** Template class for help with sorting vectors */
template<class T> 
class indexCompare {
        const T arr;
        public:
            indexCompare(const T arr) : arr(arr) {};
            bool operator()(const size_t a, const size_t b) const
                { return arr[a] < arr[b]; }
        
        };

vector<double> sortAndGetIndices(vector<double> unsorted, vector<int>& indices);


//------- Random useful stuff --------------------------------------------------

/** Split up a string according to delimeter given */
void stringSplit(string str, string delim, vector<string>& results);

/** Compute factorial */
int factorial(int n);

/** Gaussian distribution with unit mean, sigma and normalization */
double unitGaussian(double x);

/** Return a row from a data file, double format*/
vector<double> getDataFileRow(char* a, string deliminator);


}  // End namespace SOPHYA




#endif
