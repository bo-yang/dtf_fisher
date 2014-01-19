
/// \Class    fisher fisher.h "fisher.h"
///  
/// \brief
///
/// \version  1.0
/// \author   Jorge A. Sanchez
/// \date     02/08/2010

#ifndef __FISHER_H
#define __FISHER_H

#include <limits>

#include "simd_math.h"
#include "gmm.h"

// -------------------------------------------------------------------------
// Fisher Vector

struct fisher_param {
  fisher_param() :
    grad_weights(false), 
    grad_means(true), 
    grad_variances(true),
    alpha(0.5), 
    pnorm(2.0) { }
  bool grad_weights;
  bool grad_means;
  bool grad_variances;
  float alpha;
  float pnorm;
  float gamma_eps;
  void print();
};

// -------------------------------------------------------------------------

template<class T>
class fisher
{

public:
  
  fisher( fisher_param &_param );
  ~fisher( );
  
  void set_model( gaussian_mixture<T> &_gmm );

  // unweighted
  int compute( std::vector<T*> &x, T *fk );

  // weighted
  int compute( std::vector<T*> &x, std::vector<T> &wgh, T *fk);

  int dim(){ return fkdim; }

private:

  bool equal( T a, T b )
  { 
    if( fabs((T)a-(T)b)<std::numeric_limits<T>::epsilon() )
      return true;
    return false;
  }

  void alpha_and_lp_normalization( T *fk );

protected:

  fisher_param param;

  int ndim, ngauss, ngrad, fkdim;

  gaussian_mixture<T> *gmm;

  T *iwgh, *istd;
};

#endif
