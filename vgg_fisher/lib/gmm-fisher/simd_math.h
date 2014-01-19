///
/// \file	simd_math.h
/// \brief	Math routines using SIMD vector extensions
///
/// \author	Jorge Sanchez
/// \date	31/08/10 
///

#ifndef __SIMD_MATH_H
#define __SIMD_MATH_H

#include <math.h>
#include <stdio.h>

class simd
{
#if defined(__GNUC__) && defined(__SSE2__)
  typedef float  v4sf __attribute__ ((vector_size(16)));
  typedef double v2df __attribute__ ((vector_size(16)));
#endif

public:
  simd(){};

  // sum_i a_i*b_i
  static float dot_product( int ndim, const float *a, const float *b );
  static double dot_product( int ndim, const double *a, const double *b);

  // a_i *= c
  static void scale( int ndim, float *a, const float c );
  static void scale( int ndim, double *a, const double c );

  // a_i += c
  static void offset( int ndim, float *a, const float c );
  static void offset( int ndim, double *a, const double c );

  // a_i -= b_i
  static void sub( int ndim, float *a, const float *b );
  static void sub( int ndim, double *a, const double *b );

  // a_i -= b_i^2
  static void sub2( int ndim, float *a, const float *b );
  static void sub2( int ndim, double *a, const double *b );

  // a_i += b_i
  static void add( int ndim, float *a, const float *b );
  static void add( int ndim, double *a, const double *b );

  // a_i += c * b_i
  static void add( int ndim, float *a, const float *b, const float c );
  static void add( int ndim, double *a, const double *b, const double c );

  // a_i += c * b_i^2
  static void add2( int ndim, float *a, const float *b, const float c );
  static void add2( int ndim, double *a, const double *b, const double c );

  // a_i *= b_i
  static void mult( int ndim, float *a, const float *b );  
  static void mult( int ndim, double *a, const double *b );

  // sum_i (a_i-b_i)^2
  static float l2_sq( int ndim, const float *a, const float *b );
  static double l2_sq( int ndim, const double *a, const double *b );

  // sum_i c_i*(a_i-b_i)^2
  static float weighted_l2_sq( int ndim, const float *a, const float *b, const float *c );
  static double weighted_l2_sq( int ndim, const double *a, const double *b, const double *c );

  // s1_i += x_i
  // s2_i += x_i * x_i
  static void accumulate_stat( int ndim, float *s1, float *s2, const float *x );
  static void accumulate_stat( int ndim, double *s1, double *s2, const double *x );

  // s1_i += weight * x_i
  // s2_i += weight * x_i * x_i
  static void accumulate_stat( int ndim, float *s1, float *s2, const float *x, const float weight );
  static void accumulate_stat( int ndim, double *s1, double *s2, const double *x, const double weight );

  // s1_i += x_i
  static void accumulate_stat( int ndim, float *s1, const float *x );
  static void accumulate_stat( int ndim, double *s1, const double *x );

  // s1_i += weight * x_i
  static void accumulate_stat( int ndim, float *s1, const float *x, const float weight );
  static void accumulate_stat( int ndim, double *s1, const double *x, const double weight );

  // s2_i += weight * ( x_i - mu_i )^2
  static void accumulate_stat_centered( int ndim, float *s2, const float *x, const float *mu, const float weight );
  static void accumulate_stat_centered( int ndim, double *s2, const double *x, const double *mu, const double weight );

  // y_i = (x_i-m_i)*istd_i
  static void standardize( int ndim, float *y, const float *x, const float *m, const float *istd );
  static void standardize( int ndim, double *y, const double *x, const double *m, const double *istd );


};

#endif
