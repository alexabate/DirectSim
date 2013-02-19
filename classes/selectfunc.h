#ifndef SELECTFUNC_SEEN
#define SELECTFUNC_SEEN

#include <iostream>
#include "machdefs.h"
#include "genericfunc.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "sopnamsp.h"


#include "pexceptions.h"
#include "histos.h"
#include "perandom.h"
#include "tvector.h"
#include "cspline.h"
#include "fioarr.h"
#include "genericfunc.h"

#include "sinterp.h"
#include "luc.h"

#include "constcosmo.h"
#include "geneutils.h"
#include "schechter.h"

//------------------------------
// This super class define the interface for a survey selection function
class SelectionFunctionInterface {
public:
  SelectionFunctionInterface(double sfv=1.) { sfv_=sfv;}
  virtual ~SelectionFunctionInterface() { }
  // This method returns the probability of finding an object (galaxy) at the 
  // corresponding coordinates (type, redshift, glong, glat)
  virtual double SFValue(int type, double redshift, double glong, double glat)
    { return sfv_; }   
	
  inline double operator()(int type, double redshift, double glong, double glat)
    { return SFValue(type, redshift, glong, glat); }
  inline double operator()(int type, double redshift) 
    { return SFValue(type, redshift, 0., 0.); }
  inline double operator()(double redshift)
    { return SFValue(1, redshift, 0., 0.); }

  double sfv_;
};

// This class reads in an n(z) fits file and returns the value of the selection 
// function according the the corresponding distribution
class ComputedSelFunc : public SelectionFunctionInterface {
public:
// filename is a text (ascii) file with pair of (z value) on each line 
  ComputedSelFunc(string filename);
  virtual double SFValue(int type, double redshift, double glong, double glat);
  
  SInterp1D myinterp_;
};


class SelectFunc : public GenericFunc {
public:
	SelectFunc(double phistar,double Mstar,double alpha);
  //SelectFunc(Schechter& f); Try to include LF params directly from a Schechter class
  //SelectFunc(void);
  //virtual SelectFunc(void);

	void SetParam(double phistar,double Mstar,double alpha);
	void GetParam(double& phistar,double& Mstar,double& alpha);

	virtual double operator() (double LumDist,double mlim, double Mc);
	virtual double operator() (double mlim, double Mc);
	
	virtual void Print();

protected:
	double phistar_,Mstar_,alpha_; // Schechter function parameters
	Schechter sch_;
};

// This class computes a selection function using the luminosity 
// distribution function, the survey mag limite and cosmology 
// function according the the corresponding distribution
class AnalyticalSelFunc : public SelectionFunctionInterface {
public:
// filename is a text (ascii) file with pair of (z value) on each line 
  AnalyticalSelFunc(SimpleUniverse& su, SelectFunc& sf, double mlim, double Mc);
  virtual double SFValue(int type, double redshift, double glong, double glat);
  
  SimpleUniverse& su_;
  SelectFunc& sf_;
  double mlim_, Mc_;
};


#endif
