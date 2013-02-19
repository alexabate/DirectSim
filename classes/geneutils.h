/**
 * @file  geneutils.h
 * @brief Contains a series of useful, generic methods
 *
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
#include "genericfunc.h"
#include "histos.h"
#include "tvector.h"
#include "cspline.h"

#include <vector>
#include <algorithm>

// Some useful classes for interpolation and integration etc
// Better interpolations class is in sinterp.h

namespace SOPHYA {


//----------------------------------------------------
class InterpFunc {
public:
  InterpFunc(double xmin,double xmax,vector<double>& y);
  virtual ~InterpFunc(void) { }

  double XMin(void) {return _xmin;}
  double XMax(void) {return _xmax;}
  inline double X(long i) {return _xmin + i*_dx;}

  //! Return the nearest element to f giving y=f(x)
  inline double operator()(double x)
         {
         x -= _xmin;
         long i = long(x/_dx+0.5);  // to take the nearest "i"
         if(i<0) i=0; else if(i>=_nm1) i=_nm1-1;
         return _y[i];
         }

  // idem operator(double) and return
  // ok==0 if the value is found, 1 if x<xmin, 2 if x>xmax
  inline double operator()(double x,unsigned short& ok)
    {ok=0; if(x<_xmin) ok=1; else if(x>_xmax) ok=2; return (*this)(x);}

  //! Return the linear interpolation of f giving y=f(x)
  // ok==0 if the value is found, 1 if x<xmin, 2 if x>xmax
  double Linear(double x,unsigned short& ok);

  //!  Return the parabolic interpolation of f giving y=f(x)
  // ok==0 if the value is found, 1 if x<xmin, 2 if x>xmax
  double Parab(double x,unsigned short& ok);

protected:
  double _xmin,_xmax,_dx; // _dx is spacing of original vector
  long _nm1;  // n-1
  vector<double>& _y;
};

//----------------------------------------------------
class InverseFunc {
public:
  InverseFunc(vector<double>& x,vector<double>& y);
  virtual ~InverseFunc(void);
  int ComputeLinear(long n,vector<double>& xfcty);
  int ComputeParab(long n,vector<double>& xfcty);
  double YMin(void) {return _ymin;}
  double YMax(void) {return _ymax;}
protected:
  inline void find_in_y(double x,long& klo,long& khi)
    {
      long k;
      klo=0, khi=_y.size()-1;
      while (khi-klo > 1) {
	k = (khi+klo) >> 1;
	if (_y[k] > x) khi=k; else klo=k;
      }
    }

  double _ymin,_ymax;
  vector<double>& _x;
  vector<double>& _y;
};

//----------------------------------------------------
double InterpTab(double x0,vector<double>& X,vector<double>& Y,unsigned short typint=0);

//----------------------------------------------------
int FuncToHisto(GenericFunc& func,Histo& h,bool logaxex=false);
int FuncToVec(GenericFunc& func,TVector<r_8>& h,double xmin,double xmax,bool logaxex=false);

//----------------------------------------------------
double AngSol(double dtheta,double dphi,double theta0=M_PI/2.);
double AngSol(double dtheta);
double FrAngSol(double angsol);

double SinXsX(double x,bool app=false);
double SinXsX_Sqr(double x,bool app=false);

double SinNXsX(double x,unsigned long N,bool app=false);
double SinNXsX_Sqr(double x,unsigned long N,bool app=false);

//------ Integration functions

// for function with ONE variable
double IntegrateFunc(GenericFunc& func,double xmin,double xmax
         ,double perc=0.1,double dxinc=-1.,double dxmax=-1.,unsigned short glorder=4);
// for function with ONE variable and ONE parameter/ or one variable set to a constant
// parameter/constant variable must be the SECOND argument to func.
double IntegrateFunc(GenericFunc& func,double par, double xmin,double xmax
         ,double perc=0.1,double dxinc=-1.,double dxmax=-1.,unsigned short glorder=4);

double IntegrateFuncLog(GenericFunc& func,double lxmin,double lxmax
         ,double perc=0.1,double dlxinc=-1.,double dlxmax=-1.,unsigned short glorder=4);

double IntegrateFuncTrapLogSpace(GenericFunc& func,double xmin, double xmax, int nx);

void Compute_GaussLeg(unsigned short glorder,vector<double>& x,vector<double>& w,double x1=0.,double x2=1.);

//------ Some special functions -----------------------

double LogGammaFunc(double x);
double ERF(double x);
double GAMMP(double a, double x);
double GCF(double a, double x, double& gammcf,double& gln);
double GSER(double a,double x,double& gamser,double& gln);

//------ Root finder ---------------

// RootFinder Using Ridder's Method From Numerical Recipes In Fortran.
// Given the function 'func' and which is a function of x and the 
// parameters fp(np), calculate the root of the function.
// Translated into C++ by AA
double  RootRidder(GenericFunc& func,double x1, double x2, double xacc=1e-4);

//------ NaN/Inf checker ----------//

/** Check if x is not-a-number */
inline int my_isnan(double x)
	{ return x != x;}

/** Check if x is infinite */
inline int my_isinf(double x)
	{ 
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
    @param values is a vector of doubles
*/
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
    @param val is a double
*/
int findClosestElement(vector<double> values,double val);

/** Check if a vector is sorted in ascending order                            */
bool sortCheck(vector<double> values);

/** Split up a string according to delimeter given */
void stringSplit(string str, string delim, vector<string>& results);

/** Compute factorial */
int factorial(int n);

/** Gaussian distribution with unit mean, sigma and normalization */
double unitGaussian(double x);

/** Return a row from a data file, double format*/
vector<double> getDataFileRow(char* a, string deliminator);

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

}  // End namespace SOPHYA




#endif
