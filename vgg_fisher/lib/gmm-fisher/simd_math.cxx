
#include "simd_math.h"

/// \brief Dot-product
///
/// \param ndim dimensionality
/// \param a input vector
/// \param b input vector
///
/// \return sum_i a_i*b_i
///
/// \author Jorge Sánchez
/// \date Nov. 2010

float simd::dot_product( int ndim, const float *a, const float *b )
{
  float sum = 0.0;
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) ) 
  {
    float sums[4] = {0,0,0,0};
    v4sf acc = *(v4sf*)sums;
    while (ndim >= 4) 
    {
      acc = acc + *(v4sf*)a * *(v4sf*)b;
      a += 4;
      b += 4;
      ndim -= 4;
    };
    *(v4sf*)sums = acc;
    sum = sums[0] + sums[1] + sums[2] + sums[3];
  }
#endif
  while( ndim-- ) 
  {
    sum += (*(a++)) * (*(b++));
  };
  return sum;
}

double simd::dot_product( int ndim, const double *a, const double *b)
{
  double sum = 0.0;
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    double sums[2] = {0,0};
    v2df acc = *(v2df*)sums;
    while( ndim >= 2 )
    {
      acc = acc + *(v2df*)a * *(v2df*)b;
      a += 2;
      b += 2;
      ndim -= 2;
    };
    *(v2df*)sums = acc;
    sum = sums[0] + sums[1];
  }
#endif
  while( ndim-- ) 
  {
    sum += (*(a++)) * (*(b++));
  };
  return sum;
}

/// \brief Inplace scaling of dimensions: a_i *= c
///
/// \param ndim dimensionality
/// \param a vector
/// \param c constant
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::scale( int ndim, float *a, const float c )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    float _c[4];
    _c[0] = _c[1] = _c[2] = _c[3] = c;
    v4sf cv = *(v4sf*)_c;
    while( ndim >= 4 )
    {
      *(v4sf*)a = *(v4sf*)a * cv;
      a += 4;
      ndim -= 4;
    };
  }
#endif
  while( ndim-- ) 
  {
    *(a++) *= c;
  };
}

void simd::scale( int ndim, double *a, const double c )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    double _c[2];
    _c[0] = _c[1] = c;
    v2df cv = *(v2df*)_c;
    while( ndim >= 2 )
    {
      *(v2df*)a = *(v2df*)a * cv;
      a += 2;
      ndim -= 2;
    };
  }
#endif
  while( ndim-- ) 
  {
    *(a++) *= c;
  };
}

/// \brief Inplace offset of dimensions: a_i += c
///
/// \param ndim dimensionality
/// \param a vector
/// \param c constant
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::offset( int ndim, float *a, const float c )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    float _c[4];
    _c[0] = _c[1] = _c[2] = _c[3] = c;
    v4sf cv = *(v4sf*)_c;
    while( ndim >= 4 )
    {
      *(v4sf*)a = *(v4sf*)a + cv;
      a += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) += c;
  };
}

void simd::offset( int ndim, double *a, const double c )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    double _c[2];
    _c[0] = _c[1] = c;
    v2df cv = *(v2df*)_c;
    while( ndim >= 2 )
    {
      *(v2df*)a = *(v2df*)a + cv;
      a += 2;
      ndim -= 2;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) += c;
  };
}

/// \brief Inplace difference: a_i -= b_i
///
/// \param ndim dimensionality
/// \param a minuend
/// \param b substraend
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::sub( int ndim, float *a, const float *b )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    while( ndim >= 4 )
    {
      *(v4sf*)a = *(v4sf*)a - *(v4sf*)b;
      a += 4;
      b += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) -= *(b++);
  };
}

void simd::sub( int ndim, double *a, const double *b )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    while( ndim >= 2 )
    {
      *(v2df*)a = *(v2df*)a - *(v2df*)b;
      a += 2;
      b += 2;
      ndim -= 2;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) -= *(b++);
  };
}

/// \brief Inplace difference (squared substraend): a_i -= b_i^2
///
/// \param ndim dimensionality
/// \param a minuend
/// \param b substraend
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::sub2( int ndim, float *a, const float *b )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    while( ndim >= 4 )
    {
      v4sf bv = *(v4sf*)b;
      *(v4sf*)a = *(v4sf*)a - bv * bv;
      a += 4;
      b += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) -= (*b) * (*b);
    b++;
  };
}

void simd::sub2( int ndim, double *a, const double *b )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    while( ndim >= 2 )
    {
      v2df bv = *(v2df*)b;
      *(v2df*)a = *(v2df*)a - bv * bv;
      a += 2;
      b += 2;
      ndim -= 2;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) -= (*b) * (*b);
    b++;
  };
}

/// \brief Inplace sum: a_i += b_i
///
/// \param ndim dimensionality
/// \param a addend
/// \param b addend
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::add( int ndim, float *a, const float *b )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    while( ndim >= 4 )
    {
      *(v4sf*)a = *(v4sf*)a + *(v4sf*)b;
      a += 4;
      b += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) += *(b++);
  };
}

void simd::add( int ndim, double *a, const double *b )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    while( ndim >= 2 )
    {
      *(v2df*)a = *(v2df*)a + *(v2df*)b;
      a += 2;
      b += 2;
      ndim -= 2;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) += *(b++);
  };
}

/// \brief Inplace scale+sum: a_i += c * b_i
///
/// \param ndim dimensionality
/// \param a addend
/// \param b addend
/// \param c constant
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::add( int ndim, float *a, const float *b, const float c )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    float _c[4]; 
    _c[0] = _c[1] = _c[2] = _c[3] = c;
    v4sf cv = *(v4sf*)_c;
    while( ndim >= 4 )
    {
      *(v4sf*)a = *(v4sf*)a + cv * *(v4sf*)b;
      a += 4;
      b += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) += c * (*(b++));
  };
}

void simd::add( int ndim, double *a, const double *b, const double c )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    double _c[2]; _c[0]=c; _c[1]=c;
    v2df cv = *(v2df*)_c;
    while( ndim >= 2 )
    {
      *(v2df*)a = *(v2df*)a + cv * *(v2df*)b;
      a += 2;
      b += 2;
      ndim -= 2;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) += c * (*(b++));
  };
}

/// \brief Inplace scale+squared sum: a_i += c * b_i^2
///
/// \param ndim dimensionality
/// \param a sum
/// \param b addend
/// \param c constant
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::add2( int ndim, float *a, const float *b, const float c )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    float _c[4]; 
    _c[0] = _c[1] = _c[2] = _c[3] = c;
    v4sf cv = *(v4sf*)_c;
    while( ndim >= 4 )
    {
      v4sf bv = *(v4sf*)b;
      *(v4sf*)a = *(v4sf*)a + cv * bv * bv;
      a += 4;
      b += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) += c * (*b) * (*b);
    b++;
  };
}

void simd::add2( int ndim, double *a, const double *b, const double c )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    double _c[2]; _c[0]=c; _c[1]=c;
    v2df cv = *(v2df*)_c;
    while( ndim >= 2 )
    {
      v2df bv = *(v2df*)b;
      *(v2df*)a = *(v2df*)a + cv * bv * bv;
      a += 2;
      b += 2;
      ndim -= 2;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) += c * (*b) * (*b);
    b++;
  };
}

/// \brief Inplace term-by-term product: a_i *= b_i
///
/// \param ndim dimensionality
/// \param a vector
/// \param b vector
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::mult( int ndim, float *a, const float *b )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    while( ndim >= 4 )
    {
      *(v4sf*)a = *(v4sf*)a * *(v4sf*)b;
      a += 4;
      b += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) *= *(b++);
  };
}

void simd::mult( int ndim, double *a, const double *b )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    while( ndim >= 2 )
    {
      *(v2df*)a = *(v2df*)a * *(v2df*)b;
      a += 2;
      b += 2;
      ndim -= 2;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(a++) *= *(b++);
  };
}

/// \brief Squared L2-distance
///
/// \param ndim dimensionality
/// \param a vector
/// \param b vector
///
/// \return sum_i (a_i-b_i)^2
///
/// \author Jorge Sánchez
/// \date Nov. 2010

float simd::l2_sq( int ndim, const float *a, const float *b )
{
  float sum = 0.0;
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    float sums[4] = {0,0,0,0};
    v4sf acc = *(v4sf*)sums;
    while( ndim >= 4 )
    {
      v4sf dv = *(v4sf*)a - *(v4sf*)b;
      acc = acc + dv * dv;
      a += 4;
      b += 4;
      ndim -= 4;
    };    
    *(v4sf*)sums = acc;
    sum = sums[0] + sums[1] + sums[2] + sums[3];
  }
#endif
  while( ndim-- ) 
  {
    float d = (*(a++)) - (*(b++));
    sum += d * d;
  }
  return sum;
}

double simd::l2_sq( int ndim, const double *a, const double *b )
{
  double sum = 0.0;
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    double sums[2] = {0,0};
    v2df acc = *(v2df*)sums;
    while( ndim >= 2 )
    {
      v2df dv = *(v2df*)a - *(v2df*)b;
      acc = acc + dv * dv;
      a += 2;
      b += 2;
      ndim -= 2;
    };    
    *(v2df*)sums = acc;
    sum = sums[0] + sums[1];
  }
#endif
  while( ndim-- ) 
  {
    double d = (*(a++)) - (*(b++));
    sum += d * d;
  }
  return sum;
}

/// \brief Squared Weighted L2-distance
///
/// \param ndim dimensionality
/// \param a vector
/// \param b vector
/// \param c per-dimension weights
///
/// \return sum_i c_i*(a_i-b_i)^2
///
/// \author Jorge Sánchez
/// \date Nov. 2010

float simd::weighted_l2_sq( int ndim, const float *a, const float *b, const float *c )
{
  float sum = 0.0;
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    float sums[4] = {0,0,0,0};
    v4sf acc = *(v4sf*)sums;
    while( ndim >= 4 )
    {
      v4sf dv = *(v4sf*)a - *(v4sf*)b;
      acc = acc + dv * dv * *(v4sf*)c;
      a += 4;
      b += 4;
      c += 4;
      ndim -= 4;
    };    
    *(v4sf*)sums = acc;
    sum = sums[0] + sums[1] + sums[2] + sums[3];
  }
#endif
  while( ndim-- ) 
  {
    float d = (*(a++)) - (*(b++));
    sum += d * d * (*(c++));
  }
  return sum;
}

double simd::weighted_l2_sq( int ndim, const double *a, const double *b, const double *c )
{
  double sum = 0.0;
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    double sums[2] = {0,0};
    v2df acc = *(v2df*)sums;
    while( ndim >= 2 )
    {
      v2df dv = *(v2df*)a - *(v2df*)b;
      acc = acc + dv * dv * *(v2df*)c;
      a += 2;
      b += 2;
      c += 2;
      ndim -= 2;
    };    
    *(v2df*)sums = acc;
    sum = sums[0] + sums[1];
  }
#endif
  while( ndim-- ) 
  {
    double d = (*(a++)) - (*(b++));
    sum += d * d * (*(c++));
  }
  return sum;
}

/// \brief Accumulate 1st and 2nd order statistics:
///          s1_i += x_i
///          s2_i += x_i * x_i
///
/// \param ndim dimensionality
/// \param s1 1st order stat. vector
/// \param s2 2nd order stat. vector
/// \param x input vector
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::accumulate_stat( int ndim, float *s1, float *s2, const float *x )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    while( ndim >= 4 )
    {
      v4sf xv = *(v4sf*)x;
      *(v4sf*)s1 = *(v4sf*)s1 + xv;
      *(v4sf*)s2 = *(v4sf*)s2 + xv * xv;

      s1 += 4;
      s2 += 4;
      x += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(s1++) += (*x);
    *(s2++) += (*x) * (*x);
    x++;
  };
}

void simd::accumulate_stat( int ndim, double *s1, double *s2, const double *x )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    while( ndim >= 2 )
    {
      v2df xv = *(v2df*)x;
      *(v2df*)s1 = *(v2df*)s1 + xv;
      *(v2df*)s2 = *(v2df*)s2 + xv * xv;

      s1 += 2;
      s2 += 2;
      x += 2;
      ndim -= 2;
    };    
  }
#endif
  while( ndim-- ) 
  {
    *(s1++) += (*x);
    *(s2++) += (*x) * (*x);
    x++;
  };
}

/// \brief Accumulate (weighted) 1st and 2nd order statistics:
///          s1_i += weight * x_i
///          s2_i += weight * x_i * x_i
///
/// \param ndim dimensionality
/// \param s1 1st order stat. vector
/// \param s2 2nd order stat. vector
/// \param x input vector
/// \param weight sample weight
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::accumulate_stat( int ndim, float *s1, float *s2, const float *x, const float weight )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    float w[4]; 
    w[0] = w[1] = w[2] = w[3] = weight;
    v4sf wv = *(v4sf*)w;
    while( ndim >= 4 )
    {
      v4sf xv = *(v4sf*)x;
      v4sf wv_xv = wv * xv;

      *(v4sf*)s1 = *(v4sf*)s1 + wv_xv;
      *(v4sf*)s2 = *(v4sf*)s2 + wv_xv * xv;

      s1 += 4;
      s2 += 4;
      x += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    float w_x = weight * (*x);
    *(s1++) += w_x;
    *(s2++) += w_x * (*(x++));    
  };
}

void simd::accumulate_stat( int ndim, double *s1, double *s2, const double *x, const double weight )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    double w[2]; w[0]=weight; w[1]=weight;
    v2df wv = *(v2df*)w;
    while( ndim >= 2 )
    {
      v2df xv = *(v2df*)x;
      v2df wv_xv = wv * xv;

      *(v2df*)s1 = *(v2df*)s1 + wv_xv;
      *(v2df*)s2 = *(v2df*)s2 + wv_xv * xv;

      s1 += 2;
      s2 += 2;
      x += 2;
      ndim -= 2;
    };    
  }
#endif
  while( ndim-- ) 
  {
    double w_x = weight * (*x);
    *(s1++) += w_x;
    *(s2++) += w_x * (*(x++));
  };
}

/// \brief s1_i += x_i
void simd::accumulate_stat( int ndim, float *s1, const float *x )
{
  return simd::add( ndim, s1, x );
}
void simd::accumulate_stat( int ndim, double *s1, const double *x )
{
  return simd::add( ndim, s1, x );
}

/// \brief s1_i += weight * x_i
void simd::accumulate_stat( int ndim, float *s1, const float *x, const float weight )
{
  return simd::add( ndim, s1, x, weight );
}
void simd::accumulate_stat( int ndim, double *s1, const double *x, const double weight )
{
  return simd::add( ndim, s1, x, weight );
}


/// \brief Accumulate (weighted) 2nd order statistics by centering samples:
///          s1_i += weight * ( x_i - mu_i )^2
///
/// \param ndim dimensionality
/// \param s1 1st order stat. vector
/// \param x input vector
/// \param mu mean vector
/// \param weight sample weight
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Jan 2011

void simd::accumulate_stat_centered( int ndim, float *s2, const float *x, const float *mu, const float weight )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    float w[4]; 
    w[0] = w[1] = w[2] = w[3] = weight;
    v4sf wv = *(v4sf*)w;
    while( ndim >= 4 )
    {
      v4sf xmv = *(v4sf*)x - *(v4sf*)mu;
      *(v4sf*)s2 = *(v4sf*)s2 + wv * xmv * xmv;

      s2 += 4;
      x += 4;
      mu += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    float w_x = weight * (*x - *(mu++));
    *(s2++) += w_x * (*(x++));    
  };
}

void simd::accumulate_stat_centered( int ndim, double *s2, const double *x, const double *mu, const double weight )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    double w[4]; 
    w[0] = w[1] = w[2] = w[3] = weight;
    v2df wv = *(v2df*)w;
    while( ndim >= 2 )
    {
      v2df xmv = *(v2df*)x - *(v2df*)mu;
      *(v2df*)s2 = *(v2df*)s2 + wv * xmv * xmv;

      s2 += 4;
      x += 4;
      mu += 4;
      ndim -= 4;
    };    
  }
#endif
  while( ndim-- ) 
  {
    double w_x = weight * (*x - *(mu++));
    *(s2++) += w_x * (*(x++));    
  };
}

/// \brief Vector standarization
///
/// \param ndim dimensionality
/// \param y output vector
/// \param x input vector
/// \param m mean vector
/// \param istd inverse std. dev. vector
///
/// \return none
///
/// \author Jorge Sánchez
/// \date Nov. 2010

void simd::standardize( int ndim, float *y, const float *x, const float *m, const float *istd )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(float)==4) && (ndim%4==0) )
  {
    while( ndim >= 4 )
    {
      *(v4sf*)y = ( *(v4sf*)x - *(v4sf*)m ) * *(v4sf*)istd;
      y += 4;
      x += 4;
      m += 4;
      istd += 4;
      ndim -= 4;
    };
  }
#endif
  while( ndim-- ) 
  {
    *(y++) = ( *(x++) - *(m++) ) * *(istd++);
  };
}

void simd::standardize( int ndim, double *y, const double *x, const double *m, const double *istd )
{
#if defined(__GNUC__) && defined(__SSE2__)
  if( (sizeof(double)==8) && (ndim%2==0) )
  {
    while( ndim >= 2 )
    {
      *(v2df*)y = ( *(v2df*)x - *(v2df*)m ) * *(v2df*)istd;
      y += 2;
      x += 2;
      m += 2;
      istd += 2;
      ndim -= 2;
    };
  }
#endif
  while( ndim-- ) 
  {
    *(y++) = ( *(x++) - *(m++) ) * *(istd++);
  };
}
