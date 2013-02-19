/**
 * @file  sinterp.h
 * @brief Interpolation class
 *
 *
 * @author Reza Ansari, ...
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2009
 * @date 2009
 *
 */
 
#ifndef SINTERP_H_SEEN
#define SINTERP_H_SEEN

#include <iostream>
#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

using namespace std;

// To compile without SOPHYA : Comment the following two line 
#include "genericfunc.h"
#include "tarray.h"
using namespace SOPHYA; 



//-------------------------------------------

//-------------------------------------------

/** @class
  * SInterp1D class
  * 
  * Simple linear 1D interpolation class 
  *
  */
// To compile without SOPHYA : exchange the comment between the following two lines 
// class SInterp1D
class SInterp1D : public GenericFunc
{
public :
  /** Default constructor - represent the function y=x */
  SInterp1D(); 
  
  /** Constructor: Regularly spaced points in X with Y values defined by yreg */
  SInterp1D(double xmin, double xmax, vector<double>& yreg);
  
  // Interpolate to a finer regularly spaced grid, from xmin to xmax with npt points 
  // DOES NOT interpolate IF npt=0 , use xs limits if xmax < xmin 
  SInterp1D(vector<double>& xs, vector<double>& ys, double xmin=1., double xmax=-1., size_t npt=0); 

  virtual ~SInterp1D() { }
        
  double XMin() const { return xmin_; }
  double XMax() const { return xmax_; }
  double DeltaX()  { return dx_; }
  inline double X(long i) const {return xmin_ + i*dx_;}  // returns x value for index i

  // --------------------------------------------------------------
  //  Interpolated Y value as a function of X 
  double YInterp(double x) const ;
// To compile without SOPHYA : Comment the following line 
  virtual inline double operator()(double x) {  return YInterp(x); }
  // --------------------------------------------------------------
        
  // Define the interpolation points through a set of regularly spaced points on X
  void DefinePoints(double xmin, double xmax, vector<double>& yreg);
  // Interpolate to a finer regularly spaced grid, from xmin to xmax with npt points 
  // DOES NOT interpolate IF npt=0 , use xs limits if xmax < xmin 
  void DefinePoints(vector<double>& xs, vector<double>& ys, double xmin=1., double xmax=-1., size_t npt=0); 

  // Read  Y's  ( one  / line) for regularly spaced X's from file and call DefinePoints(xmin, xmax, yreg)
  size_t ReadYFromFile(string const& filename, double xmin, double xmax, int nComments=0);
  // Read pairs of X Y ( one pair / line) from file and call DefinePoints(xs, ys ...)
  size_t ReadXYFromFile(string const& filename, double xmin=1., double xmax=-1.,
         size_t npt=0, int nComments=0, bool setzero=false);
  
  vector<double>& GetVX()  { return xs_; }
  vector<double>& GetVY()  { return ys_; }

  ostream& Print(ostream& os,int lev=0) const ;
  inline ostream& Print(int lev=0) const { return Print(cout, lev); } 

protected:
  vector<double> yreg_, xs_, ys_;  // interpolated y value for regularly spaced x 
  double xmin_, xmax_, dx_;        // dx is spacing of finer grid of x's
  size_t ksmx_;                    // Maximum index value in xs_, ys_
  size_t npoints_;                 // Number of regularly spaced points, xmax not included 
  bool setzero_;				   // Set Y to zero if outside given x-vector range
};

inline ostream& operator << (ostream& s, SInterp1D const& a) 
{ a.Print(s,0);  return s; }


// Super un-sophisticated 2D interpolation
class SInterp2D : public GenericFunc
{
public :

    /** Default constructor */
    SInterp2D(){ };

    /** Constructor: linearly interpolate 2D array in @param y 
        @param xa           variable corresponding to y array rows
        @param xb           variable corresonding to y array columns 
        @param y            2D array to interpolate  
        @param isAccurate   choose interpolation method                       */
    SInterp2D(vector<double> xa, vector<double> xb, TArray<double> y, bool isAccurate=true) { 
            definePoints(xa, xb, y, isAccurate);  
            };
        
    void definePoints(vector<double> xa, vector<double> xb, TArray<double> y, bool isAccurate=true)
        {
            rangeChecks(xa, xb, y);
            
            xa_=xa;
            xb_=xb;
            y_.Set(y); // not sure this is right
        
            na_=xa_.size(); // number of rows of y
            nb_=xb_.size(); // number of columns of y
            
            dxa_ = (xa_[na_-1] - xa_[0])/(na_ - 1);
            dxb_ = (xb_[nb_-1] - xb_[0])/(nb_ - 1);
            
            isAccurate_=isAccurate;
        };
        
     /** Check the tables to interpolate make sense                           */
     void rangeChecks(vector<double>& xa, vector<double>& xb, TArray<double>& y);

    /** Interpolated Y value as a function of x1, x2. Simple bilinear interpolation, faster.                 
        @param x1   corresponds to rows direction of y (y-axis direction)
        @param x2   corresponds to columns direction of y (x-axis direction)  */
    double biLinear(double x1, double x2);
    
    /** Interpolated Y value as a function of x1, x2. More accurate bilinear interpolation                       
        @param x1   corresponds to rows direction of y (y-axis direction)
        @param x2   corresponds to columns direction of y (x-axis direction)  */
    double biLinearAccurate(double x1, double x2);

    virtual inline double operator()(double x1, double x2) {  
        if (isAccurate_)
            return biLinearAccurate(x1,x2); 
        else
            { return biLinear(x1,x2);} };
    
    /** Find closest xa element to x1, where xa<x1                            */
    int findxaElement(double x1) { 
            int ia = (int)floor((x1-xa_[0])/dxa_); 
            if (ia < 0)
                { cout << "WARNING! x1 outside of grid!" <<endl; ia = 0;}
            if (ia > na_-1 )
                { cout << "WARNING! x1 outside of grid!" <<endl; ia = na_-1;}
            return ia;
            };
            
    /** Find closest xb element to x2, where xb<x2                            */
    int findxbElement(double x2)
        { 
            int ib = (int)floor((x2-xb_[0])/dxb_); 
            if (ib < 0)
                { cout << "WARNING! x2 outside of grid!" <<endl; ib = 0;}
            if (ib > nb_-1 )
                { cout << "WARNING! x2 outside of grid!" <<endl; ib = nb_-1;}
            return ib;
            };
            
     vector<double> getColumn(int colID) {
            vector<double> yColumn;
            for (int i=0; i<na_; i++)
                yColumn.push_back(y_(i,colID));
            return yColumn;
            };
    
    
protected:
    vector<double> xa_;
    vector<double> xb_;
    TArray<double> y_;
    bool isAccurate_;
    int na_;
    int nb_;
    double dxa_;
    double dxb_;
};

typedef SInterp2D* InterpPtr;


#endif
