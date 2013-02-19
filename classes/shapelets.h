// -*- LSST-C++ -*-

/**
 * @file  shapelets.h
 * @brief 
 *
 * Could add more information here I think
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
#include "genericfunc.h"
#include "pexceptions.h"
#include "mydefrg.h"

// CatSim
#include "constcosmo.h"
#include "geneutils.h"
#include "hpoly.h"


/** BasisFuncs class
  *
  * Calculates dimensionless basis functions involving Hermite polynomials
  * 
  */
class BasisFuncs : public GenericFunc, public Hermite
{
public:
    /** Constructor */
    BasisFuncs(){ };
    
    /** Return the basis function \f$\phi_n\f$                                */
    virtual double operator()(int n, double x)
        { return phiBasisFunc(n, x); };
        
    /** Return the rescaled basis function \f$B_n\f$                          */
    virtual double operator()(int n, double x, double beta)
        { return bBasisFunc(n, x, beta); };
    
    /** Return the basis function \f$\phi_n\f$                                */
    double phiBasisFunc(int n, double x);
    
    /** Return the rescaled basis function \f$B_n\f$                          */
    double bBasisFunc(int n, double x, double beta)
            {   double xb=x/beta;
                return phiBasisFunc(n, xb)/sqrt(beta); };
    

};

/** FunctionXBasisFunction class
  *
  * Multiplies generic function with basis function \f$B_n\f$, returns as a 
  * function
  */
class FunctionXBasisFunction : public GenericFunc
{
public:
    FunctionXBasisFunction(GenericFunc& func, int n, double beta)
    : func_(func) , n_(n) , beta_(beta) { };
    
    virtual double operator()(double x)
        {   BasisFuncs basisFunc;
            return func_(x)*basisFunc.bBasisFunc(n_, x, beta_); };
    
protected:
    GenericFunc&    func_;
    int             n_;   
    double          beta_;      
};

/** Shapelets class
  *
  *
  */
class Shapelets : public GenericFunc
{
public:
    /** Constructor */
    Shapelets(GenericFunc& func)
    : func_(func) { 
        // for integration purposes	
	    npt_=1000;
	    perc_=0.1, dxinc_=0.005, dxmax_=10*dxinc_;
	    glorder_=4;
	    xmin_=-1, xmax_=1;// No idea what these should be!
	    };
    
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
             
    /** Return shapelet coefficient for order n                               */
    double shapeletCoefficient(int n, double beta)
            { 
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
    GenericFunc&    func_;      /**< function to expand into shapelets        */
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
