#ifndef SELECTFUNC_SEEN
#define SELECTFUNC_SEEN

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// sophya
#include "sopnamsp.h"
#include "machdefs.h"
#include "pexceptions.h"
#include "histos.h"
#include "perandom.h"
#include "tvector.h"
#include "cspline.h"
#include "fioarr.h"
// Sophya update v2.3 June 2013 replaces genericfunc with classfunc
//#include "genericfunc.h"
#include "classfunc.h"

// DirectSim
#include "sinterp.h"
#include "cosmocalcs.h"
#include "constcosmo.h"
#include "geneutils.h"
#include "schechter.h"

/** @class SelectionFunctionInterface
  *
  * This super class define the interface for a survey selection function
  */
class SelectionFunctionInterface {
public:
    
    /** Constructor 
        @param sfv    probabilty of finding an object                        */
    SelectionFunctionInterface(double sfv=1.) { sfv_ = sfv; };
    
    /** Desctructor */
    virtual ~SelectionFunctionInterface() { }
    
    /** Return the probability of finding an object (galaxy) with the following
        properties
        @param type        galaxy type
        @param redshift    redshift
        @param glong       galactic longitude
        @param glat        galactic latitiude                                 */
        virtual double SFValue(int type, double redshift, double glong, double glat)
            { return sfv_; };   
	
	/** Return the probability of finding an object (galaxy) with the following
        properties
        @param type        galaxy type
        @param redshift    redshift
        @param glong       galactic longitude
        @param glat        galactic latitiude                                 */
    inline double operator()(int type, double redshift, double glong, double glat)
        { return SFValue(type, redshift, glong, glat); };
        
    /** Return the probability of finding an object (galaxy) with the following
        properties
        @param type        galaxy type
        @param redshift    redshift                                           */
    inline double operator()(int type, double redshift) 
        { return SFValue(type, redshift, 0., 0.); }
  
    /** Return the probability of finding an object (galaxy) with the following
        properties
        @param redshift    redshift                                           */
    inline double operator()(double redshift)
        { return SFValue(1, redshift, 0., 0.); }

  double sfv_;    /**< probabilty of finding an object                        */
};


/** @class ComputedSelFunc
  *
  * Reads in a n(z) ascii file and returns the value of the selection 
  * function according the n(z) distribution
  */
class ComputedSelFunc : public SelectionFunctionInterface {
public:

    /** Constructor 
        @param filename    text (ascii) file with pair of [z,n(z)] on each line */
    ComputedSelFunc(string filename) { selectionfunc_.ReadXYFromFile(filename); };
    
    /** Return the selection function given the following properties (@note only
        redshift is implemented currently)
        @param type        galaxy type
        @param redshift    redshift
        @param glong       galactic longitude
        @param glat        galactic latitiude                                 */
    virtual double SFValue(int type, double redshift, double glong, double glat)
        { return selectionfunc_.YInterp(redshift); };

    SInterp1D selectionfunc_;    /**< selection function                      */
};


/** @class SelectFunc
  *
  * Compute selection function from luminosity function
  *
  */
class SelectFunc : public ClassFunc2D, ClassFunc3D  {
public:

    /** Constructor - Initialise Schechter luminosity function. See eqn 1.7 in 
        Statistics of the Galaxy Distribution by Martinez and Saar
        @param phistar    Schechter function parameter  
        @param Mstar      Schechter function parameter  
        @param alpha      Schechter function parameter                        */
	SelectFunc(double phistar, double Mstar, double alpha);
	
    //SelectFunc(Schechter& f); Try to include LF params directly from a Schechter class
    //SelectFunc(void);
    //virtual SelectFunc(void);

    /** Set Schechter function parameters
        @param phistar    Schechter function parameter  
        @param Mstar      Schechter function parameter  
        @param alpha      Schechter function parameter                        */
	void SetParam(double phistar, double Mstar, double alpha) {
        phistar_ = phistar; Mstar_ = Mstar; alpha_ = alpha; };

    /** Get Schechter function parameters
        @param phistar    Schechter function parameter  
        @param Mstar      Schechter function parameter  
        @param alpha      Schechter function parameter                        */
	void GetParam(double& phistar,double& Mstar,double& alpha) {
        phistar = phistar_; Mstar = Mstar_; alpha = alpha_; };

    /** Compute selection function given luminosity distance assuming the 
        Schechter luminosity function. Computes eqn 1.8 in Statistics of the 
        Galaxy Distribution by Martinez and Saar
        @param LumDist    luminosity distance
        @param mlim       (aparent) magnitude limit of the survey
        @param Mc         absolute magnitude survey is complete down to       */
	virtual double operator() (double LumDist,double mlim, double Mc) const;
	
	/** Compute selection function given luminosity distance assuming the 
        Schechter luminosity function. Computes eqn 1.8 in Statistics of the 
        Galaxy Distribution by Martinez and Saar
        @param Mlim       min absolute magnitude observable (eg at some distance)
        @param Mc         absolute magnitude survey is complete down to       */
	virtual double operator() (double Mlim, double Mc) const;
	
	/** Print Schechter function parameters to the screen                     */
	virtual void Print() {
        cout <<"    SelectFunc::Print: phistar="<< phistar_ <<", Mstar=";
        cout << Mstar_ <<", alpha="<< alpha_ <<endl; };

protected:
	double phistar_;    /**< Schechter function parameter                     */
	double Mstar_;      /**< Schechter function parameter                     */
	double alpha_;      /**< Schechter function parameter                     */
	//Schechter sch_;     /**< Schechter function                               */
};


/** @class  AnalyticalSelFunc
  *
  * This class computes a selection function using the luminosity 
  * distribution function, the survey mag limit and cosmology */
class AnalyticalSelFunc : public SelectionFunctionInterface {
public:

    /** Constructor
        @param su    cosmology
        @param sf    selection function given by Schechter function
        @param mlim  magnitude limit of survey
        @param Mc    absolute magnitude survey is complete down to            */
    AnalyticalSelFunc(SimpleUniverse& su, SelectFunc& sf, double mlim, double Mc)
    : su_(su), sf_(sf), mlim_(mlim), Mc_(Mc) { };
  
    /** Return selection function (@note only redshift is implemented currently)
        @param type        galaxy type
        @param redshift    redshift
        @param glong       galactic longitude
        @param glat        galactic latitiude                                 */
    virtual double SFValue(int type, double redshift, double glong, double glat);
  
    SimpleUniverse& su_;    /**< cosmology                                    */
    SelectFunc& sf_;        /**< selection function                           */
    double mlim_;           /**< magnitude limit of survey                    */
    double Mc_;             /**< absolute magnitude survey is complete down to*/
};


#endif
