/**
 * @file    igmstatistics.h
 * @brief   Contains classes to return the variance and mean of the transmission through IGM for a given redshift
 *
 * @todo    Write classes to return both the variance in and the average of the IGM transmission
 *
 * @author  Matthew Kirby
 * Contact: matthewkirby@email.arizona.edu
 *
 * Created: 2014
 * @date    2014
 *
 */


#ifndef IGMSTATISTICS_H_SEEN
#define IGMSTATISTICS_H_SEEN

#include "machdefs.h"
#include "sopnamsp.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "sinterp.h"



class IGMVariance
{
public: 
    /** Default constructor */
    IGMVariance()
        { redshift_ = 0.0; }

    /** Constructor
        Takes in a redshift and gets transmission variance for
        a galaxy at that redshift 
        @param zshift   redshift of galaxy */
    IGMVariance(double zshift); 

    /** calcNumberAbsorbers
        Calculate the average number of absorbers along a line
        of sight of length redshift_. Uses a fit that was done
        over data averaged over many lines of sight. */
    double calcNumberAbsorbers(double z);

    /** findClosestLookupTable
        Determines which lookup table shoudl be used to model
        the variance and the mean */
    void findLookupTables();

    /** calcMean
        Calculates the mean IGM Transmission given the two
        closest lookup tables */
    void calcMean();

    /** interpVectors1024
        Takes a set of vectors and interpolates to 1250 points 
        over a constant grid size between lamMin and lamMax
        @param newWL    The new wavelength vector
        @param newT     The transmission corresponding to newWL
        @param oldWL    The initial wavelength vector
        @param oldT     The initial transmission vector
        @param lamMin   The wavelength to begin the interpolation from
        @param lamMax   The wavelength to end the interpolation
      */
    void interpVectors1024(vector<double> &newWL, vector<double> &newT, 
                           vector<double> &oldWL, vector<double> &oldT,
                           double lamMin, double lamMax);

    /** combineCurves
        Takes 2 vectors and combines them using the lowZ_ and highZ_
        @param lowerT   The transmission curve for lowZ_
        @param upperT   The transmission curve for highZ_ */
    void combineCurves(vector<double> &lowerT, vector<double> &upperT);

    void returnMean(vector<double> & wl, vector<double> & mu)
        { wl = wavelengthMean_; mu = mean_; return; };
    void returnVariance(vector<double> & wl, vector<double> & sig)
        { wl = wavelengthVar_; sig = variance_; return; };

    /** testPrint
        Debugging function to print two vectors to a file
        @param lam      The x-values to be printed
        @param data     The y-values to be printed */
    void testPrint(vector<double> &lam, vector<double> &data);

protected:
    double      redshift_;          /**< The redshift of the galaxy */
    double      lookupZ_;           /**< Redshift of the lookup table being used for variance */
    double      lowZ_;              /**< Redshift of the lookup table being shifted up for mean */
    double      highZ_;             /**< Redshift of the lookup table being shifted down for mean */
    int         numAbsorbers_;      /**< The number of IGM absorbers between us and redshift_ */
    string      lookupTable_;       /**< Filename for the upper lookup table being used */
    string      lookupTableLow_;    /**< Filename for the lower lookup table being used  */

    vector<double>  wavelengthMean_;/**< The vector containing the wavelength for the mean */
    vector<double>  wavelengthVar_; /**< The vector containing the wavelength for the variance */
    vector<double>  mean_;          /**< The vector containing the mean transmission */
    vector<double>  variance_;      /**< The vector containing the variance of the transmission */
};


#endif

