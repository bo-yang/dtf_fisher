
#include "gmm.h"

/// \brief constructor
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
gaussian_mixture<T>::gaussian_mixture( const char* modelfile )
  : mean(0), var(0), coef(0), s0(0), s1(0), s2(0), i_var(0), var_floor(0), log_coef(0), log_var_sum(0), ngauss(0), ndim(0)
{  
  int err = this->load( modelfile );
  assert( err==0 );
}

template<class T>
gaussian_mixture<T>::gaussian_mixture( int n_gauss, int n_dim )
  : ngauss(n_gauss), ndim(n_dim)
{
  init();
}

template<class T>
gaussian_mixture<T>::gaussian_mixture( em_param &p, int n_gauss, int n_dim )
  : param(p), ngauss(n_gauss), ndim(n_dim)
{
  init();
}

template<class T>
void 
gaussian_mixture<T>::init()
{
  if( (ngauss==0) || (ndim==0) )
  {
    mean=0;
    var=0;
    coef=0;

    log_coef=0;
    log_var_sum=0;
    i_var=0;

    var_floor = 0;
  }
  else
  {
    mean = new T*[ngauss];
    for( int k=ngauss; k--; )
      mean[k] = new T[ndim];

    var = new T*[ngauss];
    for( int k=ngauss; k--; )
      var[k] = new T[ndim];

    coef = new T[ngauss];

    log_coef = new T[ngauss];

    log_var_sum = new double[ngauss];

    i_var = new T*[ngauss];
    for( int k=ngauss; k--; )
      i_var[k] = new T[ndim];

    // initial default values  
    for( int k=ngauss; k--; )
    {
      for( int i=ndim; i--; )
      {
        mean[k][i] = 0.0;
        var[k][i] = (T)1.0;
      }
      coef[k] = (T)1.0/(T)ngauss;
    }

    ndim_log_2pi = (T)(double(ndim)*log(2.0*M_PI));

    var_floor = 0;
  }

  s0 = 0;
  s1 = 0;
  s2 = 0;
}

/// \brief destructor
/// 
/// \param none
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
gaussian_mixture<T>::~gaussian_mixture()
{
  clean();
}

template<class T>
void
gaussian_mixture<T>::clean()
{

  delete[] coef; 
  coef=0;

  delete[] log_coef; 
  log_coef=0;

  delete[] log_var_sum; 
  log_var_sum=0;

  delete[] var_floor; 
  var_floor=0;

  delete[] s0; 
  s0=0;

  if( mean )
  {
    for( int k=ngauss; k--; )
      delete[] mean[k];
    delete[] mean;
    mean = 0;
  }

  if( var )
  {
    for( int k=ngauss; k--; )
      delete[] var[k];
    delete[] var;
    var = 0;
  }

  if( i_var )
  {
    for( int k=ngauss; k--; )
      delete[] i_var[k];
    delete[] i_var;
    i_var = 0;
  }


  if( s1 )
  {
    for( int k=ngauss; k--; ) 
    {
      if( s1[k] )
        delete[] s1[k];      
    }    
    delete[] s1;
    s1=0;
  }

  if( s2 )
  {
    for( int k=ngauss; k--; ) 
    {
      if( s2[k] )
        delete[] s2[k];      
    }    
    delete[] s2;
    s2=0;
  }

}

/// \brief set mean for the GMM 
/// 
/// \param mean new means (of the same dim)
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void
gaussian_mixture<T>::set_mean( std::vector<T*> &_mean )
{
  assert( (int)_mean.size()==ngauss );
  for( int k=0; k<ngauss; k++ )
  {
    for( int i=0; i<ndim; i++ )  
    {
      mean[k][i] = _mean[k][i];
    }
  }
}

/// \brief set variance for the GMM 
/// 
/// \param var new variances (of the same dim)
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void
gaussian_mixture<T>::set_variance( std::vector<T*> &_var )
{
  assert( (int)_var.size()==ngauss );

  for( int k=ngauss; k--; )
  {
    for( int i=ndim; i--; )
    {
      var[k][i] = std::max( (T)param.variance_floor, (T)_var[k][i] );
    }
  }

}

/// \brief set mixing coefficients for the GMM 
/// 
/// \param coef new mixing coefficients (of the same dim)
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void
gaussian_mixture<T>::set_mixing_coefficients( std::vector<T> &_coef )
{
  assert( (int)_coef.size()==ngauss );
  T sum=0.0;
  for( int k=ngauss; k--; )
    sum += _coef[k];
  assert(sum>0.0);
  for( int k=ngauss; k--; )
    coef[k] = _coef[k]/sum;
}

/// \brief set parameters for the GMM 
/// 
/// \param mean new means
/// \param var new variances
/// \param coef new mixing coefficients
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void
gaussian_mixture<T>::set( std::vector<T*> &_mean, 
                          std::vector<T*> &_var, 
                          std::vector<T>  &_coef )
{
  set_mean( _mean );
  set_variance( _var );
  set_mixing_coefficients( _coef );

  // precompute constants
  precompute_aux_var();

  // reset accumulators
  reset_stat_acc();

}

/// \bief random initialization of mean vectors
/// 
/// \param samples samples list
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void 
gaussian_mixture<T>::random_init( std::vector<T*> &samples, int seed )
{
  int N=samples.size();
  assert( N>0 );

  T dmin[ndim], dmax[ndim];
  for( int i=0; i<ndim; ++i )
    dmin[i] = dmax[i] = samples[0][i];

  for( int n=1; n<N; ++n )
  {
    for( int i=0; i<ndim; ++i )
    {
      if( samples[n][i]<dmin[i] )
        dmin[i] = samples[n][i];
      else if( samples[n][i]>dmax[i] )
        dmax[i] = samples[n][i];
    }
  }
  T m[ndim], v[ndim];
  sample_mean( samples, m, ndim );
  sample_variance( samples, m, v, ndim );
  
  srand( seed );

  for( int k=ngauss; k--; )
  {
    for( int i=ndim; i--; )
    {
      T drange = dmax[i]-dmin[i];
      mean[k][i] = dmin[i]+drange*T(rand())/T(RAND_MAX);
      var[k][i] = std::max( (T)param.variance_floor, (T)0.1*drange*drange );
    }
    coef[k] = 1.0/T(ngauss);
  }

  // precompute constants
  precompute_aux_var();

  // reset accumulators
  reset_stat_acc();
}

/// \brief Adaptive variance floornig
/// 
/// \param samples samples list
///
/// \return none
///
/// \author Jorge Sanchez
/// \date   October 2010

template<class T>
void
gaussian_mixture<T>::compute_variance_floor( std::vector<T*> &x )
{
  if( var_floor )
    delete[] var_floor;
  var_floor = new T[ndim];

  // variance of the sample set
  sample_variance( x, var_floor, ndim );

  for( int i=ndim; i--;  ) 
  {
    // a small factor of the sample variance
    var_floor[i] = param.variance_floor_factor * var_floor[i];

    // and floored
    var_floor[i] = std::max( var_floor[i], (T)param.variance_floor );
  }
}

/// \brief EM iterations
/// 
/// \param x samples
///
/// \return none
///
/// \author Jorge Sanchez
/// \date   Feb. 2012

template<class T>
void
gaussian_mixture<T>::em( std::vector<T*> &samples )
{

  nsamples = (int)samples.size();

  std::cout << "Computing variance floor..." << std::endl;
  compute_variance_floor(samples);

  std::cout << "  Number of Gaussians: " << ngauss   << std::endl;
  std::cout << "  Number of samples: "   << nsamples << std::endl;
  std::cout << "  Sample dimensions: "   << ndim     << std::endl;
  std::cout << std::endl;  

  double llh_init=0, llh_prev=0, llh_curr=0, llh_diff=0;
  for( int iter=0; iter<param.max_iter; ++iter )
  {

    // precompute constants
    precompute_aux_var();

    // reset accumulators
    reset_stat_acc();

    //  update sample statistics (E-step)
    llh_curr=0.0;
    for( int n=0; n<nsamples; ++n )
    {
      if (n % 5000 == 0) std::cout << "     (e-step): updating statistics for sample " << n << " of " << nsamples << "..." << std::endl;
      llh_curr += (double)accumulate_statistics( samples[n] );
    }
    llh_curr /= double(nsamples);

    // check for convergence
    if( iter==0 )
    {
      llh_init = llh_curr;
      std::cout << "  iter 0, avg. llh = " << llh_init << std::endl;  
    }
    else
    {
      llh_diff = (llh_curr-llh_prev)/(llh_curr-llh_init);     
      std::cout << "  iter " << std::noshowpos << iter << ", avg. llh = " << llh_curr << " (" << std::showpos << llh_diff << ")" << std::endl;      
      if( llh_diff<(double)param.llh_diff_thr )
        break;
    }
    llh_prev = llh_curr;

    // update model parameters (M-step)
    std::cout << "     (m-step): updating model..." << std::endl;
    update_model();
  }
}

/// \brief Clean accumulators
/// 
/// \param none
///
/// \return none
///
/// \author Jorge Sanchez
/// \date   Feb. 2012

template<class T>
void
gaussian_mixture<T>::reset_stat_acc()
{

  if( s0 )
    delete[] s0;
  s0 = new T[ngauss];
  memset( s0, 0, ngauss*sizeof(T));

  if( s1 )
  {
    for( int k=ngauss; k--; ) 
    {
      if( s1[k] )
        delete[] s1[k];      
    }    
    delete[] s1;
  }
  s1 = new T*[ngauss];
  for( int k=ngauss; k--; ) 
  {
    s1[k] = new T[ndim];
    memset( s1[k], 0, ndim*sizeof(T));
  }

  if( s2 )
  {
    for( int k=ngauss; k--; ) 
    {
      if( s2[k] )
        delete[] s2[k];      
    }    
    delete[] s2;
  }
  s2 = new T*[ngauss];
  for( int k=ngauss; k--; ) 
  {
    s2[k] = new T[ndim];
    memset( s2[k], 0, ndim*sizeof(T));
  }
}

/// \brief precompute auxiliary variables
/// 
/// \param none
///
/// \return none
///
/// \author Jorge Sanchez
/// \date   Feb. 2012

template<class T>
void
gaussian_mixture<T>::precompute_aux_var()
{
#pragma omp parallel for
  for( int k=0; k<ngauss; k++ ) 
  {
    log_coef[k] = log( coef[k] );

    log_var_sum[k] = 0.0;
    for( int i=ndim; i--; ) 
    {
      i_var[k][i] = (T)1.0/var[k][i];
      log_var_sum[k] += (double)log( var[k][i] );
    }
  }
}

/// \brief Accumulate statistics
/// 
/// \param x sample
///
/// \return none
///
/// \author Jorge Sanchez
/// \date   Feb. 2012

template<class T>
T
gaussian_mixture<T>::accumulate_statistics( T* x, bool _s0, bool _s1, bool _s2,
					    T* s0_ext, T** s1_ext, T** s2_ext )
{
  T* s0_active;
  T** s1_active;
  T** s2_active;

  if (s0_ext)
  {
    s0_active = s0_ext;
  } else {
    s0_active = s0;
  }
  if (s1_ext)
  {
    s1_active = s1_ext;
  } else {
    s1_active = s1;
  }
  if (s2_ext)
  {
    s2_active = s2_ext;
  } else {
    s2_active = s2;
  }
  
  T *pst=new T[ngauss];
  T llh = posterior( x, pst );

  // s0_active
  if( _s0 )
  {
    simd::add( ngauss, s0_active, pst ); 
  }
  if( _s1 )
  {
    // s1 and s2
    if( _s2 )
    {
#pragma omp parallel for
      for( int k=0; k<ngauss; ++k )
      {
        if( pst[k]<param.min_gamma )
          continue;

        simd::accumulate_stat( ndim, s1_active[k], s2_active[k], x, pst[k] );
      }
    }
    // s1 only
    else
    {
#pragma omp parallel for
      for( int k=0; k<ngauss; ++k )
      {
        if( pst[k]<param.min_gamma )
          continue;
        
        simd::add( ndim, s1_active[k], x, pst[k] );
      }
    }
  }
  // s2 only
  else if( _s2 )
  {
#pragma omp parallel for
    for( int k=0; k<ngauss; ++k )
    {
      if( pst[k]<param.min_gamma )
        continue;
        
      simd::add2( ndim, s2_active[k], x, pst[k] );
    }    
  }
  delete[] pst; pst=0;
  return llh;
}

/// \brief Update model parameters
/// 
/// \param none
///
/// \return none
///
/// \author Jorge Sanchez
/// \date   Feb. 2012

template<class T>
void
gaussian_mixture<T>::update_model()
{
  // mean and variances  
#pragma omp parallel for
  for( int k=0; k<ngauss; ++k )
  {
    if( mean[k] )
      delete[] mean[k];
    mean[k] = s1[k];
    s1[k] = 0;
    simd::scale( ndim, mean[k], (T)1.0/s0[k] );
    
    if( var[k] )
      delete[] var[k];
    var[k] = s2[k];
    s2[k] = 0;
    simd::scale( ndim, var[k], (T)1.0/s0[k] );
    simd::sub2( ndim, var[k], mean[k] );
    for( int i=ndim; i--; )
    {
      var[k][i] = std::max( var[k][i], var_floor[i] );   
    }
  }
  s1=0;
  s2=0;

  // Dirichlet prior on the mxture weights
  simd::offset( ngauss, s0, param.alpha*(T)nsamples );

  // mixing coeficients
  if( coef )
    delete[] coef;
  coef = s0;
  s0 = 0;
  T psum=0;
#pragma omp parallel for reduction(+:psum)
  for( int k=0; k<ngauss; ++k )
  {
    psum += coef[k];
  }
  simd::scale( ngauss, coef, (T)1.0/psum );  
}

/// \brief log-likelihood (and log-posterior) of a sample
///        avoiding numerical underflow.
///
/// The likelihood of a sample x can be written as:
///
///    p(x) = sum_k ( pi_k * N_k(x) )
///         = sum_k exp( log(pi_k) + log(N_k(x)) )
///
/// and the log-likelihood:
///
///   llh = log[ sum_k( exp(y_k) ) ]
///
/// where y_k = log(pi_k) + log(N_k(x)). Let's write y_max = max_k(y_k).
/// Defining y'_k = y_k-y_max, one can write:
///
///   log[ sum_k( exp(y_k) ) = log[ sum_k( exp(y'_k+y_max) ) ]
///                          = y_max + log[ sum_k( exp(y'_k) ) ]
///
/// \param x sample
/// \param log_pst k-dimensional vector (output)
///
/// \return log-likelihood of the sample
///
/// \author Jorge Sanchez
/// \date    August 2010

template<class T>
T
gaussian_mixture<T>::log_p( T *x, T *log_pst )
{

  T *lp=0;
  if( log_pst )
  {
    lp = log_pst;
  }
  else
  {
    lp = new T[ngauss];
  }

#pragma omp parallel for
  for( int k=0; k<ngauss; ++k )
  {
    lp[k] = log_gauss( k, x );
  }
  simd::add( ngauss, lp, log_coef );

  T lp_max = lp[0];
  for( int k=1; k<ngauss; ++k )
  {
    if( lp[k] > lp_max )
      lp_max = lp[k];
  }

  T log_p_sum = 0.0;
#pragma omp parallel for reduction(+:log_p_sum)
  for( int k=0; k<ngauss; ++k )
  {
    log_p_sum += (T)exp( lp[k]-lp_max );
  }
  log_p_sum = lp_max + log(log_p_sum);

  if( log_pst ) // Compute Log-Posteriors
  {
    simd::offset( ngauss, lp, -log_p_sum );
  }
  else // Just interested on the Log-likelihood of the sample
  {
    delete[] lp;
    lp = 0;
  }
  
  return log_p_sum;
}

/// \brief Gaussian Log-probability
/// 
/// \param k Gaussian component
/// \param x sample
///
/// \return log(p_k)
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
T
gaussian_mixture<T>::log_gauss( int k, T* x )
{
  T log_p = ndim_log_2pi + log_var_sum[k] 
    + simd::weighted_l2_sq( ndim, x, mean[k], i_var[k] );
  return -(T)0.5 * log_p;
}

/// \brief Posterior of a sample (responsibilities)
/// 
/// \param x sample
/// \param pst k-dimensional vector (output)
///
/// \return log-likelihood for the sample
///
/// \author Jorge Sanchez
/// \date   August 2010

template<class T>
T
gaussian_mixture<T>::posterior( T *x, T *pst )
{
  T llh=log_p( x, pst );

  T pst_sum=0.0;
#pragma omp parallel for reduction(+:pst_sum)
  for( int k=0; k<ngauss; ++k )
  {
    pst[k] = exp( pst[k] );
    pst_sum += pst[k];
  }

  if( pst_sum>0.0 )
  {
    simd::scale( ngauss, pst, (T)1.0/pst_sum );
  }

  return llh;
}

/// \brief sample log likelihoog
/// 
/// \param sample list of samples
///
/// \return average log-likelihood
///
/// \author Jorge Sanchez
/// \date   August 2009

template<class T>
T
gaussian_mixture<T>::log_likelihood( std::vector<T*> &samples )
{
  T llh=0.0;
  for( int n=0; n<(int)samples.size(); ++n )
  {
    llh += (double)log_p( samples[n] );
  }
  llh /= double(samples.size());

  return (T)llh;
}

/// \brief print model
/// 
/// \param none
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void
gaussian_mixture<T>::print( bool print_coef, bool print_mean, bool print_var )
{
  assert( print_coef || print_mean || print_var );

  for( int k=0; k<ngauss; ++k )
  {    
    if( print_coef )
      std::cout << " coef[" << k << "] = " << coef[k] << std::endl;

    if( print_mean )
    {
      std::cout << " mean[" << k << "] = [";
      for( int i=0; i<ndim; ++i )
        std::cout << " " << mean[k][i];
      std::cout << " ]" << std::endl;
    }

    if( print_var )
    {
      std::cout << " variance[" << k << "] = [";
      for( int i=0; i<ndim; ++i )
        std::cout << " " << var[k][i];
      std::cout << " ]" << std::endl;
    }
  }
}

/// \brief load model from file
/// 
/// \param none
///
/// \return -1 if error
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
int 
gaussian_mixture<T>::load( const char* filename )
{
  clean();

  std::ifstream fin( filename, std::ios::in | std::ios::binary );
  if( fin.fail() )
    return -1;

  fin.read( (char*)&ngauss, sizeof(int) );
  fin.read( (char*)&ndim,   sizeof(int) );

  init();
  
  for( int k=0; k<ngauss; k++ )
  {
    fin.read( (char*)&coef[k], sizeof(T) );
    fin.read( (char*)&mean[k][0], ndim*sizeof(T) );
    fin.read( (char*)&var[k][0], ndim*sizeof(T) );
  }

  fin.close();
  if( fin.fail() )
    return -1;

  // precompute constants
  precompute_aux_var();

  // reset accumulators
  reset_stat_acc();

  return 0;
}

/// \brief save model to file
/// 
/// \param none
///
/// \return -1 if error
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
int 
gaussian_mixture<T>::save( const char* filename )
{
  std::ofstream fout( filename, std::ios::out | std::ios::binary );
  if( fout.fail() )
    return -1;

  fout.write( (char*)&ngauss, sizeof(int) );
  fout.write( (char*)&ndim, sizeof(int) );
  
  for( int k=0; k<ngauss; k++ )
  {
    fout.write( (char*)&coef[k], sizeof(T) );
    fout.write( (char*)&mean[k][0], ndim*sizeof(T) );
    fout.write( (char*)&var[k][0], ndim*sizeof(T) );
  }

  fout.close();
  if( fout.fail() )
    return -1;

  return 0;
}

/// \brief print
/// 
/// \param none
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

void
em_param::print()
{
  std::cout << "  max_iter = " << max_iter << std::endl;
  std::cout << "  alpha = " << alpha << std::endl;
  std::cout << "  llh_diff_thr = " << llh_diff_thr << std::endl;
  std::cout << "  min_gamma = " << min_gamma << std::endl;
  std::cout << "  variance_floor = " << variance_floor << std::endl;
  std::cout << "  variance_floor_factor = " << variance_floor_factor << std::endl;
}

template class gaussian_mixture<float>;
template class gaussian_mixture<double>;

