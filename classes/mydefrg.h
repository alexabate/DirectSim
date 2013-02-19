/**
 * @file  mydefrg.h
 * @brief Wrapper for new random generator class in Sophya
 *
 *
 * @author Reza Ansari 
 * Contact: abate@email.arizona.edu
 *
 * Created on: 2009
 * @date 2009
 *
 */
#ifndef  MYDEFRG_H_SEEN
#define  MYDEFRG_H_SEEN


// We define the old name RandomGenerator as the new class ThSR48RandGen
#include "randr48.h"
//typedef  ThSDR48RandGen RandomGenerator ;  
typedef  DR48RandGen RandomGenerator ; 
#endif
