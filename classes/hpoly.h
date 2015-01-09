/**
 * @file  hpoly.h
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
 
#ifndef HPOLY_SEEN
#define HPOLY_SEEN

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

// CatSim
#include "constcosmo.h"
#include "geneutils.h"

//using namespace std;

//namespace SOPHYA {

/** Hermite class
  *
  * Class to calculate the 1D Hermite polynomials Hn(x) as a function of abscissa x
  * Faster than the astlib HERMITE function because it stores a
  * pre-compiled set of low-order hermite functions
  *
  */
class Hermite
{
public:
    /** Constructor */
    Hermite(){ setLowOrderCoeffs(); };
    
   // /** Returns the 1D Hermite polynomial of order n Hn(x)                    */
   // virtual double operator()(int order, double x)
   //     { return returnHermiteN(order, x); };

    /** Returns the 1D Hermite polynomial of order n Hn(x)                    */
    double returnHermiteN(int order, double x) const;
    
    /** Returns the Hermite coefficients for 1D polynomial of order n         */
    vector<double> getHermiteCoeffs(int order) const;
    
    /** Calculates the Hermite coefficients for 1D polynomial of order n      */
    vector<double> calculateHermiteCoeffs(int order) const;
    
    /** Set the pre-calculated lower order coefficients */
    void setLowOrderCoeffs();

protected:
    vector<double> c0_;    /**< pre-calculated coefficients                      */
    vector<double> c1_;    /**< pre-calculated coefficients                      */
    vector<double> c2_;    /**< pre-calculated coefficients                      */
    vector<double> c3_;    /**< pre-calculated coefficients                      */
    vector<double> c4_;    /**< pre-calculated coefficients                      */
    vector<double> c5_;    /**< pre-calculated coefficients                      */
    vector<double> c6_;    /**< pre-calculated coefficients                      */
    vector<double> c7_;    /**< pre-calculated coefficients                      */
    vector<double> c8_;    /**< pre-calculated coefficients                      */
    vector<double> c9_;    /**< pre-calculated coefficients                      */
    vector<double> c10_;    /**< pre-calculated coefficients                      */
    vector<double> c11_;    /**< pre-calculated coefficients                      */
    vector<double> c12_;    /**< pre-calculated coefficients                      */
    vector<double> c13_;    /**< pre-calculated coefficients                      */
    vector<double> c14_;    /**< pre-calculated coefficients                      */
    vector<double> c15_;    /**< pre-calculated coefficients                      */
    vector<double> c16_;    /**< pre-calculated coefficients                      */
    vector<double> c17_;    /**< pre-calculated coefficients                      */
    vector<double> c18_;    /**< pre-calculated coefficients                      */
    vector<double> c19_;    /**< pre-calculated coefficients                      */
    vector<double> c20_;    /**< pre-calculated coefficients                      */
    vector<double> c21_;    /**< pre-calculated coefficients                      */
    vector<double> c22_;    /**< pre-calculated coefficients                      */
    vector<double> c23_;    /**< pre-calculated coefficients                      */
    vector<double> c24_;    /**< pre-calculated coefficients                      */
    vector<double> c25_;    /**< pre-calculated coefficients                      */
    
};
#endif
