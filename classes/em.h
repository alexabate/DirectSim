/**
 * @file  em.h
 * @brief Expectation-maximization algorithm for estimating a mixture of n
 *        Gaussians
 *
 * @todo This algorithm seems to fail sometimes when the Gaussians are
 *       spaced far apart. Probably should code up a metric that flags 
 *       up when this happens - a simple (data-fit)^2 value should do
 *
 * @author Alex Abate, Michael Chen (Matlab version)
 * Contact: abate@email.arizona.edu
 *
 * Created on: 29/Jan/2013
 * @date 29/Jan/2013
 *
 */
#ifndef  EM_H_SEEN
#define  EM_H_SEEN


// generic stuff
#include <iostream>
#include <fstream>

// sophya stuff
#include "tarray.h"
#include "randinterf.h"
#include "sopnamsp.h"

// my stuff
#include "mydefrg.h"
#include "geneutils.h"
#include "matrix.h"
#include "constcosmo.h"

struct modelGaussian {

    // length nGaussian_
    vector<double> mus;
    vector<double> sigmas;
    vector<double> weights;
    // length of number of iterations performed
    vector<double> llhs; 
    // length of data
    vector<int> labels;
    int convNum;
    
};

class EMAlgorithm 
{
public:

    /** Default constructor */
    //EMAlgorithm(){};
    
    /** Constructor 
        @param data         data to fit Gaussians to 
        @param rg           random number generator 
        @param nGaussian    number of Gaussians to fit data to 
        @param tol          tolerance */
    EMAlgorithm(vector<double> data, RandomGeneratorInterface& rg, int nGaussian=3, 
            double tol=1e-5, int maxIter=500)
    : data_(data) , rg_(rg) , nGaussian_(nGaussian) , tol_(tol) , maxIter_(maxIter)
    {};
    
    ~EMAlgorithm(){};
    
    /** Expectation-maximization algorithm */
    modelGaussian doExpectationMaximization();
    
    
    // LOCAL METHODS
    
    /** Initialize the algorithm */
    TArray<double> initialization();
    
    void maximization(TArray<double> R);
    
    TArray<double> expectation();
    
    double logGaussPDF(double datum, double mu, double sigma);
    
    /** Initialize the Gaussian model */
    void initializeModel(int n) { 
        for (int i=0; i<nGaussian_; i++) {
            model_.mus.push_back(1e10); // infinity
            model_.sigmas.push_back(-1);
            model_.weights.push_back(-1);
            }
        for (int j=0; j<n; j++)
            model_.labels.push_back(-1);
        model_.convNum = -1;
        };
        
    /** Compute log(sum(exp(x),dim)) while avoiding numerical underflow */
    TVector<double> computeLogSumExp(TArray<double> logRho, int n);

    void getLabels(TArray<double> R, int n);
    
protected:
    vector<double> data_;           /**< data to fit Gaussians to */ 
    RandomGeneratorInterface& rg_;  /**< random generator */
    int nGaussian_;                 /**< number of Gaussians to fit to the data */
    double tol_;                    /**< tolerance */
    double maxIter_;                /**< maximum number of iterations */
    modelGaussian model_;           /**< structure holding the current model parameters */
};

#endif
