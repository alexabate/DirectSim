// -*- LSST-C++ -*-
#include "shapelets.h"

/******* BasisFuncs methods ***************************************************/

double BasisFuncs::phiBasisFunc(int n, double x)
{

    if (n<0)
        throw ParmError("ERROR! Order of polynomial is negative!");
        
    double A = 1./sqrt(pow(2.,n)*sqrt(PI)*factorial(n));
    double Hn = returnHermiteN(n,x);
    double G = unitGaussian(x);
    return A*Hn*G;

};
