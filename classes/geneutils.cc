#include "machdefs.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <vector>

#include "pexceptions.h"

#include "histos.h"
#include "hisprof.h"
#include "srandgen.h"

#include "geneutils.h"

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAXSTP 10000
#define TINY 1.0e-30

namespace SOPHYA {
//-------------------------------------------------------------------
// InterpFuncXY:
// Linear interpolation class when the elements are NOT regularly spaced
/*
InterpFuncXY::InterpFuncXY(vector<double>& xs,vector<double>& ys, double xmin, double xmax, int npt)
{
  DefinePoints(xs, ys, xmin, xmax, npt);
}

InterpFuncXY::DefinePoints(vector<double>& xs,vector<double>& ys, double xmin, double xmax, int npt)
{
  xmin_ = xmin; 
  xmax_ = xmax;
  npoints_ = npt;
 
  
  if(xmin_>=_xmax_|| npoints_()<=2) {
   cout<<"InterpFuncXY::InterpFuncXY : bad argument values"<<endl;
   throw ParmError("InterpFuncX::InterpFuncXY : bad argument values");
	}
   
  yreg_.SetSize(npoints_+1);
  dx_ = (xmax_-xmin_)/(double)npoints_;
  for(int k=0; k<npoints_+1; k++)  
	{
    double x = k*dx_;
	yreg_[k]= ys[k] + (ys[k+1]-ys[k])*(x-xs[k])/(xs[k+1]-xs[k]);
	}
}


*/


/******* InterpFunc **********************************************************/
// InterpFunc:
// Class of linear interpolation
// The vector y has n elements y_i such as y_i = f(x_i)
//   where thee x_i are regularly spaced
//   and x_0=xmin and x_{n-1}=xmax   (xmax included!)
InterpFunc::InterpFunc(double xmin,double xmax,vector<double>& y)
  : _xmin(xmin), _xmax(xmax), _y(y)
{
  if(_xmin>=_xmax || _y.size()<=2) {  // less than 3 points!
   cout<<"InterpFunc::InterpFunc : bad arguments values"<<endl;
   throw ParmError("InterpFunc::InterpFunc : bad arguments values");
  }
  _nm1   = _y.size()-1;
  _dx    = (_xmax-_xmin)/(double)_nm1;
}

double InterpFunc::Linear(double x,unsigned short& ok)
{
  ok=0; if(x<_xmin) ok=1; else if(x>_xmax) ok=2;
  x -= _xmin;
  long i = long(x/_dx);  // to take the "i" just below
  if(i<0) i=0; else if(i>=_nm1) i=_nm1-1;
  return _y[i] + (_y[i+1]-_y[i])/_dx*(x-i*_dx);
}

double InterpFunc::Parab(double x,unsigned short& ok)
{
  ok=0; if(x<_xmin) ok=1; else if(x>_xmax) ok=2;
  x -= _xmin;
  long i = long(x/_dx+0.5);  // to take the nearest "i"
  if(i<1) i=1; else if(i>=_nm1-1) i=_nm1-2;
  double a = (_y[i+1]-2.*_y[i]+_y[i-1])/(2.*_dx*_dx);
  double b = (_y[i+1]-_y[i-1])/(2.*_dx);
  return _y[i] + (x-i*_dx)*(a*(x-i*_dx)+b);
}


/******* InverseFunc **********************************************************/
InverseFunc::InverseFunc(vector<double>& x, vector<double>& y)
  : _ymin(0.) , _ymax(0.) , _x(x) , _y(y)
{
    uint_4 ns = _x.size();
    if(ns<3 || _y.size()<=0 || ns!=_y.size())
        throw ParmError("InverseFunc::InverseFunc_Error: bad array size");

    // Check "x" is growing monotonically
    for(uint_4 i=0;i<ns-1;i++)
        if(_x[i+1]<=_x[i]) {
            cout<<"InverseFunc::InverseFunc_Error: _x array not stricly growing"<<endl;
            throw ParmError("InverseFunc::InverseFunc_Error: _x array not stricly growing");
            }

    // Check "y" is growing monotonically
    for(uint_4 i=0;i<ns-1;i++)
        if(_y[i+1]<_y[i]) {
            cout<<"InverseFunc::InverseFunc_Error: _y array not growing"<<endl;
            throw ParmError("InverseFunc::InverseFunc_Error: _y array not growing");
            }

    // define limits
    _ymin = _y[0];
    _ymax = _y[ns-1];

}

// Compute table "xfcty" by linear interpolation of "x" versus "y"
// on "nout" points from "ymin" to "ymax":
// xfcty[i] = interpolation of function "x" for "ymin+i*(ymax-ymin)/(nout-1)"
int InverseFunc::ComputeLinear(long nout,vector<double>& xfcty)
{
    if(nout<3) return -1;

    xfcty.resize(nout);

    long i1,i2;
    double x;
    for(int_4 i=0;i<nout;i++) {
    
        double y = _ymin + i*(_ymax-_ymin)/(nout-1.);
        find_in_y(y, i1, i2);
        double dy = _y[i2]-_y[i1];
        
        if(dy==0.) 
            x = (_x[i2]+_x[i1])/2.; // the function to invert is flat!
        else 
            x = _x[i1] + (_x[i2]-_x[i1])/dy * (y-_y[i1]);
   
        xfcty[i] = x;
        }

    return 0;
}

// Compute table "xfcty" by parabolic interpolation of "x" versus "y"
// on "nout" points from "ymin" to "ymax":
// xfcty[i] = interpolation of function "x" for "ymin+i*(ymax-ymin)/(nout-1)"
int InverseFunc::ComputeParab(long nout,vector<double>& xfcty)
{
    if(nout<3) return -1;

    xfcty.resize(nout);

    long i1,i2,i3;
    double x;
    for(int_4 i=0;i<nout;i++) {
    
        double y = _ymin + i*(_ymax-_ymin)/(nout-1.);
        find_in_y(y,i1,i2);

        // Find the 3rd point according to the position of y / 2 next to it
        double my = (_y[i1]+_y[i2])/2.;

        if(y<my) {i3=i2; i2=i1; i1--;} else {i3=i2+1;}

        // Protection
        if(i1<0) {i1++; i2++; i3++;}
        if(i3==(long)_y.size()) {i1--; i2--; i3--;}
        
        // parabolic interpolation 
        double dy = _y[i3]-_y[i1];
        if(dy==0.) 
            x = (_x[i3]+_x[i1])/2.; // the function to invert is flat!
        else {
            double X1=_x[i1]-_x[i2], X3=_x[i3]-_x[i2];
            double Y1=_y[i1]-_y[i2], Y3=_y[i3]-_y[i2];
            double den = Y1*Y3*dy;
            double a = (X3*Y1-X1*Y3)/den;
            double b = (X1*Y3*Y3-X3*Y1*Y1)/den;
            y -= _y[i2];
            x = (a*y+b)*y + _x[i2];
            }
        xfcty[i] = x;
        }

    return 0;
}

/******* InterpTab **********************************************************/
double InterpTab(double x0, vector<double>const& X, vector<double> const& Y, unsigned short typint)
// Interpole in x0 the table Y = f(X)
//           X doit etre ordonne par ordre croissant (strictement)
// typint = 0 : nearest value
//          1 : linear interpolation
//          2 : parabolique interpolation
{
 long n = X.size();
 if(n>(long)Y.size() || n<2)
   throw ParmError("InterpTab_Error :  incompatible size between X and Y tables!");

 if(x0<X[0] || x0>X[n-1]) return 0.;
 if(typint>2) typint = 0;

 // Recherche des indices encadrants par dichotomie
 long k, klo=0, khi=n-1;
 while (khi-klo > 1) {
   k = (khi+klo) >> 1;
   if (X[k] > x0) khi=k; else klo=k;
 }

 // Quel est le plus proche?
 k = (x0-X[klo]<X[khi]-x0) ? klo: khi;

 // On retourne le plus proche
 if(typint==0) return Y[k];

 // On retourne l'extrapolation lineaire
 if(typint==1 || n<3)
   return Y[klo] + (Y[khi]-Y[klo])/(X[khi]-X[klo])*(x0-X[klo]);

 // On retourne l'extrapolation parabolique
 if(k==0) k++; else if(k==n-1) k--;
 klo = k-1; khi = k+1;
 double x1 = X[klo]-X[k], x2 = X[khi]-X[k];
 double y1 = Y[klo]-Y[k], y2 = Y[khi]-Y[k];
 double den = x1*x2*(x1-x2);
 double a = (y1*x2-y2*x1)/den;
 double b = (y2*x1*x1-y1*x2*x2)/den;
 x0 -= X[k];
 return Y[k] + (a*x0+b)*x0;;

}

//-------------------------------------------------------------------
int FuncToHisto(ClassFunc1D& func, Histo& h, bool logaxex)
// Fill the histo 1D "h" with the function "func"
// INPUT:
// logaxex = false : fill linearly
//        the abscisses "x" of the bins are filled with f(x)
// logaxex = true : fill logarithmically (base 10)
//        the abscisses "x" of the bins are filled f(10^x)
// RETURN:
//       0 = OK
//       1 = error
{
    if(h.NBins()<=0) return 1;

    h.Zero();

    for(int_4 i=0;i<h.NBins();i++) {
    
        double x = h.BinCenter(i);
        if(logaxex) x = pow(10.,x);
        h.SetBin(i,func(x));
        }

    return 0;
}

int FuncToVec(ClassFunc1D& func, TVector<r_8>& v, double xmin, double xmax, bool logaxex)
// Fill the TVector with the function "func"
// INPUT:
// logaxex = false : fill linearly
//        the abscisses "x" of the bins are filled with f(x)
// logaxex = true : fill logarithmically (base 10)
//        the abscisses "x" of the bins are filled with f(10^x)
// RETURN:
//       0 = OK
//       1 = error
// Note:
//  v(i) = f(xmin+i*dx) with dx = (xmax-xmin)/v.NEts()
{
    if(v.NElts()<=0 || xmax<=xmin) return 1;

    v = 0.;
    double dx = (xmax-xmin)/v.NElts();
   
    for(int_4 i=0;i<v.NElts();i++) {
        double x = xmin + i * dx;;
        if(logaxex) x = pow(10.,x);
        v(i) = func(x);
        }

    return 0;
}

//-------------------------------------------------------------------
double AngSol(double dtheta, double dphi, double theta0)
// Return the solid angle of a "rectangle" and spherical coordinates
// of half-angles "dtheta" x "dphi" centered on "theta0"
// Note: the "theta0" at the equator is Pi/2  (and not zero)
//          The units are in radians:
//          theta0 in [0,Pi]
//          dtheta in [0,Pi]
//          dphi   in [0,2Pi]
// Return: solid angle in steradians
{
    double theta1 = theta0-dtheta, theta2 = theta0+dtheta;
    if(theta1<0.) theta1=0.;
    if(theta2>M_PI) theta2=M_PI;

    return   2.*dphi * (cos(theta1)-cos(theta2));
}

double AngSol(double dtheta)
// Return solid angle of a spherical cap with angle opening "dtheta" 
// Note: The units are in radians:
//               dtheta in [0,Pi]
// Return: solid angle in steradians
// Approximation for small theta: PI * theta^2
{
    return   2.*M_PI * (1.-cos(dtheta));
}

double FrAngSol(double angsol)
// Return the opening angle "dtheta" of a spherical cap with solid angle "angsol"
// Input: solid angle in steradians
// Return: opening angle in radians
{
    angsol = 1. - angsol/(2.*M_PI);
    if(angsol<-1. || angsol>1.) return -1.;
    return acos(angsol);
}

//-------------------------------------------------------------------
double SinXsX(double x, bool app)
// Calculation of sin(x)/x
// if x^2 < 1.7e-4 or app=true returns:
// 1 - x^2/6 * ( 1 - x^2/20 * (1 - x^2/42) )
{
    double x2 = x*x;
    if(app || x2<1.7e-4) return 1. - x2/6.*(1. - x2/20.*(1. - x2/42.));
    return sin(x)/x;
}

double SinXsX_Sqr(double x, bool app)
// Calculation of (sin(x)/x)^2
// Approximation: if x^2 < 6.8e-5 returns:
// 1. - x^2/3 * ( 1 - 2*x^2/15 * (1 - x^2/14) )
{
    double x2 = x*x;
    if(app || x2<6.8e-5) return 1. - x2/3.*(1. - 2.*x2/15.*(1. - x2/14.));
    x2 = sin(x)/x;
    return x2*x2;
}

double SinNXsX(double x, unsigned long N, bool app)
// Calculation of sin(N*x)/sin(x) where N is a positive integer
// ATTENTION: N is integer 
//  1. In case sin(N*x) and sin(x) are both zero
//     (do DL popur sin(N*x) and sin(x))
//  2. To treat DL correctly in x = p*Pi+e  with e<<1 and p an integer
//     sin(N*x)/sin(x) = sin(N*p*Pi+N*e)/sin(p*Pi+e)
//                     = [sin(N*p*Pi)*cos(N*e)+cos(N*p*Pi)*sin(N*e)]
//                         / [sin(p*Pi)*cos(e)+cos(p*Pi)*sin(e)]
//                     when sin(N*p*Pi)=0
//                     = [cos(N*p*Pi)*sin(N*e)] / [cos(p*Pi)*sin(e)]
//                     = [sin(N*e)/sin(e)] * [cos(N*p*Pi)/cos(p*Pi)]
//                     = [DL autour de x=0] * (+1 ou -1)
// The only case where we have "-1" is when "p=odd" (cos(p*Pi)=-1) and "N=pair"
// Approx: if sin(x)^2 < 3.5e-6/N^2
{
    if(N==0) return 0;
    
    double sx = sin(x), N2 = N*N;
    if (app || sx*sx<3.5e-6/N2) {
        double x2 = asin(sx); x2 *= x2;
        double s = 1.;
        if (N%2==0 && cos(x)<0.) s = -1.; // when x ~ (2p+1)*Pi and N is even
        return s*N*(1.-(N-1.)*(N+1.)/6.*x2*(1.-(3.*N2-7.)/60.*x2));
        }
    return sin((double)N*x)/sx;
}

double SinNXsX_Sqr(double x, unsigned long N, bool app)
// Calculation of [sin(N*x)/sin(x)]^2  where N is positive integer
// ATTENTION: cf remark as N integer in SinNXsX method
// Approx: if sin(x)^2 ~ 1.5e-6/N^2
{
    if(N==0) return 0;

    double sx = sin(x), N2 = N*N;
    if(app || sx*sx<1.5e-6/N2) {
        double x2 = asin(sx); x2 *= x2;
        return N2*(1. - (N-1.)*(N+1.)/3.*x2*(1. - (2.*N2-3.)/15.*x2));
        }
    sx = sin((double)N*x)/sx;
    return sx*sx;
}

//-------------------------------------------------------------------

static unsigned short IntegrateFunc_GlOrder = 0;
static vector<double> IntegrateFunc_x;
static vector<double> IntegrateFunc_w;

///////////////////////////////////////////////////////////
//************** Integration of Functions ***************//
///////////////////////////////////////////////////////////


///// ONE VARIABLE //////////////////////////////////////////////////
double IntegrateFunc(ClassFunc1D const& func, double xmin, double xmax,
    double perc, double dxinc, double dxmax, unsigned short glorder)
// --- Adaptative integration ---
//     Integrate[func[x], {x,xmin,xmax}]
// ..xmin,xmax are the integration limits
// ..dxinc is the searching increment x = xmin+i*dxinc
//         if <0  npt = int(|dxinc|)   (if<2 then npt=100)
//                and dxinc = (xmax-xmin)/npt
// ..dxmax is the maximum possible increment (if dxmax<=0 no test)
// ---
// Split [xmin,xmax] in intervals [x(i),x(i+1)]:
// Traverse [xmin,xmax] in steps of "dxinc" : x(i) = xmin + i*dxinc
// After creating interval [x(i),x(i+1)]
//     - if |f[x(i+1)] - f[x(i)]| > perc*|f[[x(i)]|
//     - if |x(i+1)-x(i)| >= dxmax  (if dxmax>0.)
// In the interval [x(i),x(i+1)] the function is integrated
// by the Gauss-Legendre method with the order: "glorder"
{

	// perform checks
	double signe = 1.;
	if (xmin>xmax) 
		{ double tmp=xmax; xmax=xmin; xmin=tmp; signe=-1.; }
	if (glorder==0) 
		glorder = 4;
	if (perc<=0.) 
		perc=0.1;
	if (dxinc<=0.) 
		{ int n=int(-dxinc); if(n<2) n=100; dxinc=(xmax-xmin)/n; }
	if (glorder != IntegrateFunc_GlOrder) {
		IntegrateFunc_GlOrder = glorder;
		Compute_GaussLeg(glorder, IntegrateFunc_x, IntegrateFunc_w, 0., 1.);
  		}

	// Search the interval: [x(i),x(i+1)]
	int_4 ninter = 0;
	double sum = 0., xbas=xmin, fbas=func(xbas), fabsfbas=fabs(fbas);
	for (double x=xmin+dxinc; x<xmax+dxinc/2.; x += dxinc) {
	
		double f = func(x);
		double dx = x-xbas;
		bool doit = false;
		// checks whether to do integral
    		if ( x>xmax ) 
			{doit = true; x=xmax;}
		else if ( dxmax>0. && dx>dxmax ) 
			doit = true;
    		else if ( fabs(f-fbas)>perc*fabsfbas ) 
			doit = true;
		if ( !doit ) 
			continue;
    		double s = 0.;
		for (unsigned short i=0;i<IntegrateFunc_GlOrder;i++)//default up to 4
      			s += IntegrateFunc_w[i]*func(xbas+IntegrateFunc_x[i]*dx);
		sum += s*dx;
		xbas = x; fbas = f; fabsfbas=fabs(fbas); ninter++;
		}
    //cout<<"Ninter="<<ninter<<endl;

    return sum*signe;
}

///// ONE VARIABLE, ONE PARAMETER //////////////////////////////////
// PARAMETER "PAR" MUST BE THE SECOND ARGUMENT TO FUNC
double IntegrateFunc(ClassFunc2D const& func, double par, double xmin, double xmax,
    double perc, double dxinc, double dxmax, unsigned short glorder)
// --- Adaptative integration ---
//     Integrate[func[x], {x,xmin,xmax}]
// ..xmin,xmax are the integration limits
// ..dxinc is the searching increment x = xmin+i*dxinc
//         if <0  npt = int(|dxinc|)   (if<2 then npt=100)
//                and dxinc = (xmax-xmin)/npt
// ..dxmax is the maximum possible increment (if dxmax<=0 no test)
// ---
// Split [xmin,xmax] in intervals [x(i),x(i+1)]:
// Traverse [xmin,xmax] in steps of "dxinc" : x(i) = xmin + i*dxinc
// After creating interval [x(i),x(i+1)]
//     - if |f[x(i+1)] - f[x(i)]| > perc*|f[[x(i)]|
//     - if |x(i+1)-x(i)| >= dxmax  (if dxmax>0.)
// In the interval [x(i),x(i+1)] the function is integrated
// by the Gauss-Legendre method with the order: "glorder"
{
    // perform checks
    double signe = 1.;
    if (xmin>xmax) 
        { double tmp=xmax; xmax=xmin; xmin=tmp; signe=-1.; }
    if (glorder==0) 
        glorder = 4;
    if (perc<=0.) 
        perc=0.1;
    if (dxinc<=0.) 
        { int n=int(-dxinc); if(n<2) n=100; dxinc=(xmax-xmin)/n; }
    if (glorder != IntegrateFunc_GlOrder) {
        IntegrateFunc_GlOrder = glorder;
        Compute_GaussLeg(glorder, IntegrateFunc_x, IntegrateFunc_w, 0., 1.);
        }

    // Search the interval: [x(i),x(i+1)]
    int_4 ninter = 0;
    double sum = 0., xbas=xmin, fbas=func(xbas,par), fabsfbas=fabs(fbas);
    for(double x=xmin+dxinc; x<xmax+dxinc/2.; x += dxinc) {
    
        double f = func(x,par);
        double dx = x-xbas;
        bool doit = false;
        // checks whether to do integral
        if ( x>xmax ) 
            {doit = true; x=xmax;}
        else if ( dxmax>0. && dx>dxmax ) 
            doit = true;
        else if ( fabs(f-fbas)>perc*fabsfbas ) 
            doit = true;
        if ( !doit ) 
            continue;
        double s = 0.;
        for (unsigned short i=0;i<IntegrateFunc_GlOrder;i++)
            s += IntegrateFunc_w[i]*func(xbas+IntegrateFunc_x[i]*dx,par);
        sum += s*dx;
        xbas = x; fbas = f; fabsfbas=fabs(fbas); ninter++;
        }
    //cout<<"Ninter="<<ninter<<endl;

    return sum*signe;
}

////////////////////////////////////////////////////////////////////////////////////
double IntegrateFuncLog(ClassFunc1D const& func, double lxmin, double lxmax,
    double perc, double dlxinc, double dlxmax, unsigned short glorder)
// --- Adaptative integration ---
// Idem IntegrateFunc but integrate on logarithmic (base 10) intervals:
//    Int[ f(x) dx ] = Int[ x*f(x) dlog10(x) ] * log(10)
// ..lxmin,lxmax are the log10() of the integration limits
// ..dlxinc is the searching logarithmic (base 10) increment lx = lxmin+i*dlxinc
// ..dlxmax is the maximum possible logarithmic (base 10) increment (if dlxmax<=0 no test)
// Remark: to be used if "x*f(x) versus log10(x)" looks like a polynomial
//            better than "f(x) versus x"
// ATTENTION: the function func that is passed as an argument
//            is "func(x)" and not "func(log10(x))"
{

    // perform checks
    double signe = 1.;
    if (lxmin>lxmax) 
        { double tmp=lxmax; lxmax=lxmin; lxmin=tmp; signe=-1.; }
    if (glorder==0) 
        glorder = 4;
    if (perc<=0.) 
        perc=0.1;
    if (dlxinc<=0.) 
        { int n=int(-dlxinc); if(n<2) n=100; dlxinc=(lxmax-lxmin)/n; }
    if (glorder != IntegrateFunc_GlOrder) {
        IntegrateFunc_GlOrder = glorder;
        Compute_GaussLeg(glorder, IntegrateFunc_x, IntegrateFunc_w, 0., 1.);
        }

    // Search some intervals [lx(i),lx(i+1)]
    int_4 ninter = 0;
    double sum = 0., lxbas=lxmin, fbas=func(pow(10.,lxbas)), fabsfbas=fabs(fbas);
    for (double lx=lxmin+dlxinc; lx<lxmax+dlxinc/2.; lx += dlxinc) {
    
        double f = func(pow(10.,lx));
        double dlx = lx-lxbas;
        bool doit = false;
        // checks whether to do integral
        if ( lx>lxmax ) 
            {doit = true; lx=lxmax;}
        else if( dlxmax>0. && dlx>dlxmax ) 
            doit = true;
        else if( fabs(f-fbas)>perc*fabsfbas ) 
            doit = true;
        if( ! doit ) 
            continue;
        double s = 0.;
        for (unsigned short i=0;i<IntegrateFunc_GlOrder;i++) {
            double y = pow(10.,lxbas+IntegrateFunc_x[i]*dlx);
            s += IntegrateFunc_w[i]*y*func(y);
            }
        sum += s*dlx;
        lxbas = lx; fbas = f; fabsfbas=fabs(fbas); ninter++;
         }
    //cout<<"Ninter="<<ninter<<endl;

    return M_LN10*sum*signe;
}

////////////////////////////////////////////////////////////////////////////////////
/*

 Simple trapesium integration function: function is log_e spaced

*/
double IntegrateFuncTrapLogSpace(ClassFunc1D const& func,double xmin, double xmax, int nx)
{

	double logxmin = log(xmin);
	double dlogx = (log(xmax) - logxmin) / (nx-1);
	
	
	double sum = 0;
	for(int i=0; i<nx; i++)
		{
	 	double logx = logxmin + i*dlogx;
    		double f = func(exp(logx));
		sum += f*dlogx;
		}
 
	return sum;
}

////////////////////////////////////////////////////////////////////////////////
/*
Integration of one dimensional functions y=f(x) via Gauss-Legendre method.
  --> Calculation of coefficients for [x1,x2]
| INPUT:
|  x1,x2 : interval limits (in nbinteg.h -> x1=-0.5 x2=0.5)
|  glorder = degree n of Gauss-Legendre
| OUTPUT:
|  x[] = x values where coefficients were calculated (dim=n)
|  w[] = value of associated coefficients
| REMARKS:
|  - x and w have dimensions up to n.
|  - the integration is rigorous if on the interval of integration
|    the function f(x) can be approximated by a polynomial
|    of degree 2*m (monotomic the + haut x**(2*n-1) )
*/
void Compute_GaussLeg(unsigned short glorder, vector<double>& x, vector<double>& w, double x1, double x2)
{
    if(glorder==0) return;
    int n = (int)glorder;
    x.resize(n,0.); w.resize(n,0.);

    int m,j,i;
    double z1,z,xm,xl,pp,p3,p2,p1;

    m=(n+1)/2;
    xm=0.5*(x2+x1);
    xl=0.5*(x2-x1);
    for (i=1;i<=m;i++)  {
    
        z=cos(M_PI*(i-0.25)/(n+0.5));
        do {
            p1=1.0;
            p2=0.0;
            for (j=1;j<=n;j++) {
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
            }  while (fabs(z-z1) > 3.0e-11);  // epsilon
        x[i-1]=xm-xl*z;
        x[n-i]=xm+xl*z;
        w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
        w[n-i]=w[i-1];
        }
}

//********************* Special functions ************************//
// the functions below are (fairly) blindly copied from Andrew Zentner's fortran 
// code

// Log of gamma function
double LogGammaFunc(double x)
{
    double ser=1.000000000190015;
	const double stp = 2.5066282746310005;
	double gammln, y, tmp;

	/****
	vector<double> cof;
	cof.push_back(76.18009172947146);
	cof.push_back(-86.50532032941677);
	cof.push_back(24.01409824083091);
	cof.push_back(-1.231739572450155);
	cof.push_back(0.1208650973866179e-2);
	cof.push_back(-0.5395239384953e-5);
	****/

	double cof[6] = {76.18009172947146,
			 -86.50532032941677,
			 24.01409824083091,
			 -1.231739572450155,
			 0.1208650973866179e-2,
			 -0.5395239384953e-5};


    y=x;
    tmp=x+5.5;
    tmp = (x+0.5)*log(tmp)-tmp;
    for (int j=0;j<6;j++) {
        y += 1.0;
        double cof_i = cof[j];
        ser+= cof_i/y;
		}
    gammln=tmp+log(stp*ser/x);
        
    return gammln;

}

// Error function
double ERF(double x)
{

	double erf;
	if (x<0.0)
        	erf=-GAMMP(0.5,x*x);
	else
        erf=GAMMP(0.5,x*x);
      
	return erf;
}

// Incomplete gamma function
double GAMMP(double a, double x)
{    
    double gammp;
    double gammcf,gamser,gln;

	if ( (x<0.0) || (a<=0.0) )
        cout <<"BAD ARGUMENTS IN GAMMP"<<endl;
	if ( x < (a+1.0) ) {// use the series representation
        	GSER(a,x,gamser,gln);
        	gammp=gamser;
		}
	else {// use the continued fraction representation
        	GCF(a,x,gammcf,gln);
        	gammp=1.0-gammcf;// and take its complement
    		}

    return gammp;
};

// Returns the incomplete gamma function evaluated by its continued
// fraction representation, also returns the log of the gamma ftn
double GCF(double a, double x, double& gammcf, double& gln)
{

	int itmax=200;     // maximum number of iterations
	double eps =3e-7;  // relative accuracy
	double fpmin=1e-30;// number near the smallest representable floating point number

	gln = LogGammaFunc(a);// checked a=0.5, gln=0.572365
	//cout << "gln="<<gln<<endl;


	// Set up for evaluating continued fraction
	// by modified Lentz's method with b0=0
	double b=x+1.0-a;
	double c=1.0/fpmin;
	double d=1.0/b;
	double h=d;

	for (int i=0;i<itmax;i++) // iterate to convergence
		{
 	        double an=-(i+1)*((i+1)-a);
		b=b+2.0;
		d=an*d+b;
		if( abs(d) < fpmin )
			d=fpmin;
        	c=b+an/c;
        	if( abs(c) < fpmin )
			c=fpmin;
        	d=1.0/d;
        	double del=d*c;
        	h=h*del;
        	if( abs(del-1.0) < eps )
			{ 
			gammcf=exp(-x+a*log(x)-gln)*h; 
			return gammcf;
			}
		}

	cout <<"A TOO LARGE, ITMAX TOO SMALL IN GCF"<<endl; 
	gammcf=-99;

	return gammcf;
}

// Returns the incomplete gamma function evaluated by its series 
// representation
// Also returns log of gamma function
double GSER(double a,double x,double& gamser,double& gln)
{

 
	int itmax=200; // max number of iterations
	double eps =3e-7;// relative accuracy, fpmin=1e-30;

	gln = LogGammaFunc(a);// checked a=0.5, gln=0.572365

	if( x <= 0.0 )
		{
	        if( x<0.0 )
			cout <<"X < 0 IN GSER'"<<endl;
		gamser=0.0;	
        	return gamser;
		}
	double ap=a;
	double sum=1.0/a;
	double del=sum;

	for (int n=0;n<itmax; n++)
		{
		ap=ap+1.0;
		del=del*x/ap;
		sum=sum+del;
		if( abs(del)< abs(sum)*eps )
			{
			gamser=sum*exp(-x+a*log(x)-gln); 
			return gamser; 
			}
		}

     	cout <<"A TOO LARGE, ITMAX TOO SMALL IN GSER"<<endl ; 
	gamser=-99;

	return gamser;
}

//******* Root finder ********************************************************//

// From Numerical recipes, Ridder's method.
// Using the false position method, find the root of a function func, known to lie
// between x1 and x2.  The root rtf is refined until its accuracy is +/- xacc
double  RootRidder(ClassFunc1D const& func,double x1, double x2, double xacc)
{

	int maxit = 30; // maximum number of iterations

	//double xl,xh,swap;
	//double dx, del, f, rtf;

	double fl=func(x1);
	double fh=func(x2);

	if ( (fl>0.0 && fh<0.0) || (fl<0.0 &&fh>0.0) ) {
	
		double xl=x1;
		double xh=x2;
		double ans=-9.99e99; // any highly unlikely value 
		
		for (int j=0; j<maxit;j++) {
		
			double xm=0.5*(xl+xh);
			double fm=func(xm);// 1st of 2 function evaluations per iteration
			double s=sqrt(fm*fm-fl*fh);
			if (s==0.0)
				return ans;
			double xnew = xm+(xm-xl)*((fl>=fh ? 1.0 : -1.0)*fm/s); // updating formula
			if (abs(xnew-ans) <= xacc)
				return ans;
			ans=xnew;
			double fnew=func(ans); // 2nd of 2 function evaluations per iteration
			if(fnew == 0.0)
				return ans;
			if(SIGN(fm,fnew)!=fm) {// bookeeping to keep the root bracketed on next iteration
				xl = xm;
				fl = fm;
				xh = ans;
				fh = fnew;
				}
			else if (SIGN(fl,fnew) !=fl) {
				xh = ans;
				fh = fnew;
				}
			else if (SIGN(fh,fnew) !=fh) {
				xl=ans;
				fl=fnew;
				}
			else
				throw("RootRidder: never get here.");
			if (abs(xh-xl)<=xacc)
				return ans;
			}
		throw("RootRidder exceeded maximum iterations");
		}
	else {
		if (fl == 0.0)
			return x1;
		if (fh == 0.0)
			return x2;
		throw("RootRidder: root must be bracketed in RootRidder");
		}
};

//******* Useful *************************************************************//

double findMaximum(vector<double> values)
{
    double maxval=-1e12;
    int nv = values.size();
    int iMax = -1;
    for (int i=0; i<nv; i++) {
        
        if ( maxval<values[i] ) {
            maxval = values[i];
            iMax = i;
            }
            
        }
    if (iMax<0)
        throw ParmError("ERROR! Did not find maximum in vector");
        
    return maxval;

};


double findMaximum(vector<double> values, int& iMax)
{
    double maxval=-1e12;
    int nv = values.size();
    iMax = -1;
    for (int i=0; i<nv; i++) {
        
        if ( maxval<values[i] ) {
            maxval = values[i];
            iMax = i;
            }
            
        }
    if (iMax<0)
        throw ParmError("ERROR! Did not find maximum in vector");
        
    return maxval;

};

double findMaximumPosition(TArray<double> array, int& iRow, int& jCol)
{

    // set i,j to unphysical indices
    iRow = -1;
    jCol = -1;

    double maxval=-1e12;
    int nRows = array.SizeX();
    int nCols = array.SizeY();
    
    for (int i=0; i<nRows; i++) {
        for (int j=0; j<nCols; j++) {
        
            if ( maxval<array(i,j) ){
                maxval = array(i,j);
                iRow = i;
                jCol = j;
                }
            }
        }
        
        
    if ( (iRow<0)||(jCol<0) )
        throw ParmError("ERROR! Did not find maximum value in array");
        
    return maxval;
    

};



double findMinimum(vector<double> values)
{
    double minval=1e12;
    int nv = values.size();
    for (int i=0; i<nv; i++) {
        
        if ( minval>values[i] )
            minval = values[i];
            
        }
        
    return minval;

};


double findMinimumPosition(TArray<double> array, int& iRow, int& jCol)
{

    // set i,j to unphysical indices
    iRow = -1;
    jCol = -1;

    double minval=1e12;
    int nRows = array.SizeX();
    int nCols = array.SizeY();
    
    for (int i=0; i<nRows; i++) {
        for (int j=0; j<nCols; j++) {
        
            if ( minval>array(i,j) ){
                minval = array(i,j);
                iRow = i;
                jCol = j;
                }
            
            }
        }
        
    if ( (iRow<0)||(jCol<0) )
        throw ParmError("ERROR! Did not find minimum value in array");
        
    return minval;
    

};


double findMinimumPosition(vector<double> array, int& iElement)
{

 // set iElement to unphysical index
    iElement = -1;

    double minval=1e12;
    int nElement = array.size();
    
    for (int i=0; i<nElement; i++) {
        
        if ( minval>array[i] ){
            minval = array[i];
            iElement = i;
            }

        }
        
    if ( (iElement<0) )
        throw ParmError("ERROR! Did not find minimum value in array");
        
    return minval;

};


int findClosestElement(vector<double> values,double val)
{

    double minval=1e100;
    int iElement = -1;
    int nv = values.size();
    for (int i=0; i<nv; i++) {
        
        double diff = abs(val-values[i]);
        
        if ( minval>diff ) {
            iElement = i;
            minval = diff;
            }
            
        }
        
    if ( iElement < 0)
        throw ParmError("findClosestElement failed");
        
    return iElement;


};

bool sortCheck(vector<double> values)
{
    int nVals = values.size();
    bool isSorted = true; // assume sorted
    for (int i=0; i<(nVals-1); i++) {
        if ( values[i+1]<values[i] )
            isSorted = false;
        if (!isSorted)// if fails sort test immediately return isSorted=false
            return isSorted;
        }
       
    return isSorted;
    
};

bool checkNonDecreasing(vector<double> xs)
{

    int nx = xs.size();
    
    bool isIncreasing = true;
    for (int i=0; i<(nx-1); i++) {
        if ( xs[i+1]<=xs[i] )
            isIncreasing = false;
        }
    return isIncreasing;

};

vector<int> uniqueVector(vector<int> vect,vector<int>& ids)
{

    int n = vect.size();
    vector<int> uniqueValues;

    int cnt = 0;
    for (int i=0; i<n; i++) {
        int val = vect[i];
        //cout <<" vector value = "<< val << endl;
        if (i<1) {
            //cout <<" Adding 1st element "<<endl;
            uniqueValues.push_back(val);
            ids.push_back(cnt); cnt++;
            
            }
        else {
            
            bool isRepeat = false;
            for (int j=0; j<uniqueValues.size(); j++){
                //cout <<" checking if value is repeated "<< endl;
                int uval = uniqueValues[j];
                //cout << uval << endl;
                if (val == uval ) {
                    //cout <<" REPEAT!"<<endl;
                    ids.push_back(j);
                    isRepeat = true; break;
                    }
                }
            
            if (!isRepeat){
                //cout <<" NO REPEAT!"<<endl;
                uniqueValues.push_back(val);
                ids.push_back(cnt); cnt++;
                }
            //cout << endl << endl;
            }
        
        }
    return uniqueValues;

};

/******* for sorting **********************************************************/

vector<double> sortAndGetIndices(vector<double> unsorted, vector<int>& indices)
{
    for (unsigned i = 0; i<unsorted.size(); ++i)
        indices.push_back(i); // b = [0, 1, 2]
    
    sort(indices.begin(), indices.end(), indexCompare<vector<double>&>(unsorted));

    vector<double> sorted;
    for (unsigned i = 0; i<unsorted.size(); ++i)
        sorted.push_back(unsorted[indices[i]]);
        
    return sorted;
}

/******* Random useful ********************************************************/

void stringSplit(string str, string delim, vector<string>& results)
{
    int cutAt;
	while( (cutAt = str.find_first_of(delim)) != str.npos )	{
	
		if(cutAt > 0)
			results.push_back(str.substr(0,cutAt));
		str = str.substr(cutAt+1);
		
        }
        
	if(str.length() > 0)
		results.push_back(str);
		
};



double unitGaussian(double x)
{

    return exp(-(x*x)/2.);

};


vector<double> getDataFileRow(char* a, string deliminator)
{


    string row = a;
    vector<string> rowArray;
    stringSplit(a, deliminator, rowArray);
    
    vector<double> rowArrayValues;
    int nRow = rowArray.size();
    for (int j=0; j<nRow; j++)
        rowArrayValues.push_back(atof(rowArray[j].c_str()));
    
    return rowArrayValues;

};

int factorial(int n)
{

    int i=0, fact=1;
    if( n<= 1)
        return fact;
    else {

        for(i=1; i<=n; i++)
            fact*=i;

        return(fact);

        }

};


}  // End namespace SOPHYA
