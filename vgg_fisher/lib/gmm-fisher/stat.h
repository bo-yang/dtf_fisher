///  
/// \brief    Statistics
///
/// \version  1.0
/// \author   Jorge A. Sanchez
/// \date     28/02/2008
///

#ifndef __STAT_H
#define __STAT_H

#include <iostream>
#include <string.h>
#include <math.h>
#include <vector>
#include <assert.h>
#include <limits>

template<class T> void sample_mean( const std::vector<T*> &samples, T *mean, int dim );
template<class T> void sample_variance( const std::vector<T*> &samples, const T *mean, T *variance, int dim );
template<class T> void sample_variance( const std::vector<T*> &samples, T *variance, int dim );
template<class T> void standardize( std::vector<T*> &samples, T *mean, T *variance, int dim );
template<class T> void destandardize( std::vector<T*> &samples, const T *mean, const T *variance, int dim );

#endif
