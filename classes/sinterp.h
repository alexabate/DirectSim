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

// To compile without SOPHYA : Comment the following two lines 
// Sophya update v2.3 June 2013 replaces genericfunc with classfunc
#include "classfunc.h" 
//#include "genericfunc.h"
#include "tarray.h"
using namespace SOPHYA; 

/** @class
  * SInterp1D class
  * 
  * Simple linear 1D interpolation class 
  *
  */
// To compile without SOPHYA : exchange the comment between the following two lines 
// class SInterp1D
class SInterp1D : public ClassFunc1D
{
public :
    /** Default constructor - represents the function y=x */
    SInterp1D(); 
  
    /** Constructor: Regularly spaced points in X with Y values defined by 
        @param yreg
        @param xmin     x value corresponding to yreg[0]
        @param xmax     x value corresponding to last yreg value              
        @param yreg     y values of function                                  */
    SInterp1D(double xmin, double xmax, vector<double>& yreg);
  
    /** Constructor: pairs of (x,y) points, no need to be regularly spaced in x,
        but MUST be sorted. Interpolates to a finer regularly spaced grid, from 
        xmin to xmax with npt points. If npt<1 interpolation is done on 
        irregular spaced grid. xs limits are used if xmax < xmin
        @param xs    x values of the points
        @param ys    y values of the points
        @param xmin  min x value to do interpolation from
        @param xmax  max x value to do interpolation from
        @param npt   number of points to use in interpolation                 */
    SInterp1D(vector<double>& xs, vector<double>& ys, 
                                 double xmin=1., double xmax=-1., size_t npt=0); 

    /** Destructor */
    virtual ~SInterp1D() { };

    /** Return min x value of interpolation table                             */
    double XMin() const { return xmin_; };
    
    /** Return max x value of interpolation table                             */
    double XMax() const { return xmax_; }
  
    /** Return dx step of interpolation table                                 */
    double DeltaX()  { return dx_; }
  
    /** Returns x value for index i in interpolation table                    */
    inline double X(long i) const {return xmin_ + i*dx_;}  

    // --------------------------------------------------------------
    /** Return interpolated Y value as a function of X                        */ 
    double YInterp(double x) const ;
  
    // To compile without SOPHYA : Comment the following line 
    /** Return interpolated Y value as a function of X                        */
    virtual inline double operator()(double x) const {  return YInterp(x); }
    // --------------------------------------------------------------
        
    /** Define the interpolation table through a set of regularly spaced points
        in X
        @param xmin    min x value of interpolation table
        @param xmax    max x value of interpolation table
        @param yreg    y values of interpolation table                        */             
    void DefinePoints(double xmin, double xmax, vector<double>& yreg);
    
    /** Define interpolation table through pairs of (x,y) points, no need to be 
        regularly spaced in x as it interpolates them to a finer regularly 
        spaced grid, from xmin to xmax with npt points. If npt<1 interpolation 
        is done on irregular spaced grid. xs limits are used if xmax < xmin. 
        xs must be sorted.
        @param xs    x values of the points
        @param ys    y values of the points
        @param xmin  min x value of interpolation table
        @param xmax  max x value of interpolation table
        @param npt   number of points to use in interpolation (or npt<1 for no
                     table)                                                   */
    void DefinePoints(vector<double>& xs, vector<double>& ys, 
                                 double xmin=1., double xmax=-1., size_t npt=0); 

    /** Read Y's (one per line) (for regularly spaced X's) from a file and call 
        DefinePoints(xmin, xmax, yreg) 
        @param filename    file to read y values from
        @param xmin        x value corresponding to first y value read
        @param xmax        x value corresponding to last y value read
        @param nComments   number of comment lines at the top of the file     */
    size_t ReadYFromFile(string const& filename, double xmin, double xmax, int nComments=0);
  
    /** Read pairs of X Y (one pair per line) from a file and call 
        DefinePoints(xs, ys ...) 
        @param filename    file to read x,y values from
        @param xmin        min x value to do interpolation from
        @param xmax        max x value to do interpolation from
        @param npt         number of points to use in interpolation
        @param nComments   number of comment lines at the top of the file     
        @param setzero     if true and x outside xmin to xmax range interpolated
                           y value = zero                                     */
    size_t ReadXYFromFile(string const& filename, double xmin=1., double xmax=-1.,
                             size_t npt=0, int nComments=0, bool setzero=false);
  
    /** Return x values of interpolation table                                */
    vector<double>& GetVX()  { return xs_; };
    
    /** Return y values of interpolation table                                */
    vector<double>& GetVY()  { return ys_; };

    /** Print properties of interpolation to the output stream
        @param os    output stream to print to                                
        @param lev   level of printing                                        */
    ostream& Print(ostream& os, int lev=0) const ;
  
    /** Print properties of interpolation to the screen
        @param lev   level of printing                                        */
    inline ostream& Print(int lev=0) const { return Print(cout, lev); };

protected:
  vector<double> yreg_;  /**< interp table y value for regularly spaced x     */
  vector<double> xs_;    /**< x value of x,y pair used to make interp table   */
  vector<double> ys_;    /**< y value of x,y pair used to make interp table   */
  double xmin_;          /**< min x of interpolation table                    */
  double xmax_;          /**< max x of interpolation table                    */
  double dx_;            /**< dx spacing of interpolation table               */
  size_t ksmx_;          /**< maximum index value in xs_, ys_                 */
  size_t npoints_;       /**< num of regularly spaced points, xmax not incl.  */ 
  bool setzero_;		 /**< set Y to zero if outside given x range          */
};


inline ostream& operator << (ostream& s, SInterp1D const& a) 
{ a.Print(s,0);  return s; };


// Super un-sophisticated 2D interpolation
// Also NOT properly checked since sophya update
class SInterp2D : public ClassFunc1D
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
            
    /** This is defined to override the pure virtual function defined in ClassFunc1D
        otherwise SInterp2D is sometimes treated as an abstract class        */
    virtual double operator() (double) const { };
        
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
