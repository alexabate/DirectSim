/**
 * @file  shapelets.h
 * @brief 
 *
 *
 * @author Alex Abate
 * Contact: abate@email.arizona.edu
 *
 * Created on: 18 Oct 2012
 * @date 18 Oct 2012
 *
 */
 
#ifndef SHAPELETS_SEEN
#define SHAPELETS_SEEN

#include "machdefs.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>

// sophya
#include "sopnamsp.h"
// Sophya update v2.3 June 2013 replaces genericfunc with classfunc
//#include "genericfunc.h"
#include "classfunc.h"
#include "pexceptions.h"
#include "mydefrg.h"

// DirectSim
#include "constcosmo.h"
#include "geneutils.h"
#include "hpoly.h"


/** @class BasisFuncs 
  *
  * Calculates dimensionless basis functions involving Hermite polynomials
  * 
  */
class BasisFuncs : public ClassFunc1D, public Hermite
{
public:
    /** Constructor */
    BasisFuncs(){ };
    
    /** This is defined to override the pure virtual function defined in ClassFunc1D
        otherwise BasisFuncs is sometimes treated as an abstract class        */
    virtual double operator() (double) const { };
    
    /** Return the basis function \f$\phi_n(x)\f$
        @param n    order of Hermite polynomial
        @param x    value basis function evalulated at                        */
    virtual double operator()(int n, double x) const
        { return phiBasisFunc(n, x); };
        
    /** Return the rescaled basis function \f$B_n(x)\f$
        @param n     order of Hermite polynomial
        @param x     value basis function evalulated at                        
        @param beta  scale size                                               */
    virtual double operator()(int n, double x, double beta) const
        { return bBasisFunc(n, x, beta); };
    
    /** Return the basis function \f$\phi_n(x)\f$                                
        @param n    order of Hermite polynomial
        @param x    value basis function evalulated at                        */
    double phiBasisFunc(int n, double x) const;
    
    /** Return the rescaled basis function \f$B_n\f$
        @param n     order of Hermite polynomial
        @param x     value basis function evalulated at                        
        @param beta  scale size                                               */
    double bBasisFunc(int n, double x, double beta) const {
        double xb = x/beta;
        return phiBasisFunc(n, xb)/sqrt(beta); };
    

};


/** @class FunctionXBasisFunction 
  *
  * Multiplies generic function with basis function \f$B_n\f$, returns as a 
  * function
  */
class FunctionXBasisFunction : public ClassFunc1D
{
public:
    /** Constructor 
        @param func    function to multiply with basis function
        @param n       order of Hermite polynomial
        @param beta    shapelet scale size                                    */
    FunctionXBasisFunction(ClassFunc1D& func, int n, double beta)
    : func_(func) , n_(n) , beta_(beta) { };
    
    virtual double operator()(double x) const
        {   BasisFuncs basisFunc;
            return func_(x)*basisFunc.bBasisFunc(n_, x, beta_); };
    
protected:
    ClassFunc1D& func_; /**< function to multiply with basis function         */
    int n_;             /**< order of Hermite polynomial                      */
    double beta_;       /**< shapelet scale size                              */
};


/** @class Shapelets
  *
  * Expand a function into shapelets
  */
class Shapelets : public ClassFunc1D
{
public:

    /** Constructor 
        @param func    function to expand into shapelets                      */
    Shapelets(ClassFunc1D& func)
    : func_(func) { 
        // for integration purposes	
	    npt_=1000;
	    perc_=0.1, dxinc_=0.005, dxmax_=10*dxinc_;
	    glorder_=4;
	    xmin_=-1, xmax_=1;// No idea what these should be!
	    };
    
    /** Return function expanded into shapelets
        @param x
        @param nmax  maximum order of Hermite polynomial
        @param beta  shapelet scale size                                      */
    double functionRep(double x, int nmax, double beta) {
            BasisFuncs basisfuncs;
            double fx=0;
            for (int n=0; n<=nmax; n++) {
                double fn = shapeletCoefficient(n, beta);
                double bn = basisfuncs(n, x, beta);
                fx+=(fn*bn);
                }
            return fx;
            };
             
    /** Return shapelet coefficient for order n                               
        @param n     order of Hermite polynomial
        @param beta  shapelet scale size                                      */
    double shapeletCoefficient(int n, double beta) {
            FunctionXBasisFunction funcXbasisFunc(func_, n, beta);
            double fn = IntegrateFunc(funcXbasisFunc,xmin_,xmax_,perc_,dxinc_,dxmax_,glorder_);
            return fn; 
            };
            
    /** Set integration parameters                                            */
    void setIntegration(int npt, double perc, double dxinc, unsigned short glorder)
        { npt_ = npt; perc_ = perc; dxinc_ = dxinc; glorder_ = glorder; };
        
    /** Set integration limits                                                */
    void setIntegLimts(double xmin, double xmax)
        { xmin_ = xmin; xmax_ = xmax; };
          
protected:
    ClassFunc1D&    func_;      /**< function to expand into shapelets        */
    // integration parameters   
	int npt_;			        /**< number of points to use in integration   */
	double perc_;               /**< integration parameter                    */
	double dxinc_;              /**< integration parameter                    */
	double dxmax_;              /**< integration parameter                    */
	unsigned short glorder_;    /**< integration parameter                    */
    double xmin_;               /**< integration lower limit                  */
    double xmax_;               /**< integration upper limit                  */
};

#endif
