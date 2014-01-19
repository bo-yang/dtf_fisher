
#ifndef __GMM_H
#define __GMM_H

#include <vector>
#include <assert.h>
#include <iostream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <limits>

#include "stat.h"
#include "simd_math.h"

struct em_param 
{
  em_param(): 
    max_iter(100), 
    alpha(100.0),
    llh_diff_thr(0.001),
    min_gamma(1e-4),
    variance_floor(1.0e-9),
    variance_floor_factor(0.01) {}
  int max_iter;         // max. number of EM iterations
  float alpha;          // Dirichlet prior on mixture weights: alpha_k=alpha (for all k)
  float llh_diff_thr;   // average Log-Likelihood difference threshold
  float grow_factor;    // growing factor for split
  float min_gamma;      // min. posterior prob. for a sample
  float variance_floor; // hard variance floor
  float variance_floor_factor; // factor for the adaptive flooring
  void print();
};

/// \class    gaussian_mixture gmm.h "gmm.h"
///  
/// \brief    Gaussian (diagonal covariances) Mixture Model using EM-algorithm
///
/// \author   Jorge Sanchez
/// \date     29/07/2009

template<class T>
class gaussian_mixture 
{
public:
  template<class TT> friend class fisher;

  gaussian_mixture( const char* modelfile );

  gaussian_mixture( int n_gauss, int n_dim );

  gaussian_mixture( em_param &p, int n_gauss, int n_dim );

  ~gaussian_mixture();

  void set( std::vector<T*> &_mean, std::vector<T*> &_var, std::vector<T>  &_coef );
  void random_init( std::vector<T*> &samples, int seed=-1 );

  void em( std::vector<T*> &samples );

  T posterior( T* sample, T *pst );

  T log_likelihood( std::vector<T*> &samples );

  int n_dim(){ return ndim; }
  int n_gauss(){ return ngauss; }

  int load( const char* filename );
  int save( const char* filename );
  void print( bool _coef=true, bool _mean=false, bool _var=false );
  
  inline T* get_mean( int idx ) { return mean[idx]; }
  inline T* get_variance( int idx ) { return var[idx]; }
  inline T get_mixing_coefficients( int idx ) { return coef[idx]; }

  inline em_param get_params() { return param; }

  void set_mean( std::vector<T*> &_mean );
  void set_variance( std::vector<T*> &_var );
  void set_mixing_coefficients( std::vector< T > &_coef ); 

protected:

  T **mean, **var, *coef;

  void reset_stat_acc();
  T accumulate_statistics( T* sample, bool _s0=true, bool _s1=true, bool _s2=true,
			   T* s0_ext=0, T** s1_ext=0, T** s2_ext=0 );
  T *s0, **s1, **s2;

  em_param param;

  void init();

  void clean();

  void compute_variance_floor( std::vector<T*> &x );

  void precompute_aux_var();

  void update_model();

  T log_p( std::vector<T*> &samples );

  T **i_var, *var_floor, *log_coef;

  double *log_var_sum; // accumulate as double

  int ngauss, ndim, nsamples;

  T ndim_log_2pi;

  T log_p( T* x, T *log_pst=0 );

  T log_gauss( int k, T* x );

private:

};

#endif
