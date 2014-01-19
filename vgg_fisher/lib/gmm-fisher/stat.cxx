
#include "stat.h"

/// \bief sample mean
/// 
/// \param samples sample list (input)
/// \param mean mean (output)
/// \param dim dimension of the input samples
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void
sample_mean( const std::vector<T*> &samples, T *mean, int dim )
{
  int N=samples.size();
  assert(N>0);

  double *acc = new double[dim];
  memset( acc, 0, dim*sizeof(double) );

  for( int n=N; n--; )  
  {
    for( int i=dim; i--; )
    {
      acc[i] += (double)samples[n][i];
    }
  }

  double iN = 1.0/double(N);
  for( int i=dim; i--; )
  {
    mean[i] = (T)(acc[i]*iN);
  } 

  delete[] acc;
}


/// \bief sample variance (over dimensions)
/// 
/// \param samples sample list (input)
/// \param mean precomputed sample mean (input)
/// \param variance variance (output)
/// \param dim dimension of the input samples
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void
sample_variance( const std::vector<T*> &samples, const T *mean, T *variance, int dim )
{
  int N=samples.size();
  assert(N>1);

  double *acc = new double[dim];
  memset( acc, 0, dim*sizeof(double) );

  for( int n=N; n--; )  
  {
    for( int i=dim; i--; )
    {
      T dm = samples[n][i]-mean[i];
      acc[i] += (double)(dm*dm);
    }
  }

  double iN1 = 1.0/double(N-1);
  for( int i=dim; i--; )
  {
    variance[i] = (T)(acc[i]*iN1);
  }

  delete[] acc;
}

/// \bief sample variance (over dimensions)
/// 
/// \param samples sample list (input)
/// \param variance variance (output)
/// \param dim dimension of the input samples
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void
sample_variance( const std::vector<T*> &samples, T *variance, int dim )
{
  T *mean = new T[dim];
  sample_mean( samples, mean, dim );
  sample_variance( samples, mean, variance, dim );
  delete[] mean;
  return;
}

/// \bief Sample standardization (-> zero mean, unit variance)
/// 
/// \param samples sample list (input) -> standarized samples (output)
/// \param mean sample mean (output)
/// \param variance sample variance (output)
/// \param dim dimension of the input samples
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void 
standardize( std::vector<T*> &samples, T *mean, T *variance, int dim )
{

  sample_mean( samples, mean, dim );
  sample_variance( samples, mean, variance, dim );  

  T i_stddev[dim];
  for( int i=dim; i--; )
  {
    if( variance[i]>0.0 )
      i_stddev[i] = 1.0/sqrt(variance[i]);
    else
      i_stddev[i] = 0.0;
  }

#pragma omp parallel for
  for( int n=0; n<(int)samples.size(); n++ )
  {
    for( int i=0; i<dim; i++ )   
    {
      samples[n][i] = (samples[n][i]-mean[i])*i_stddev[i];
    }
  }

}

/// \bief Samples de-standardization (zero mean, unit variance -> original mean and var)
/// 
/// \param samples standardized sample list (input) -> original standarized samples (output)
/// \param mean original sample mean (input)
/// \param variance original sample variance (input)
/// \param dim dimension of the input samples
///
/// \return none
///
/// \author Jorge Sanchez
/// \date    August 2009

template<class T>
void
destandardize( std::vector<T*> &samples, const T *mean, const T *variance, int dim )
{
  assert(mean && variance);

  T stddev[dim];
  for( int i=0; i<dim; i++ )
  {
    stddev[i] = sqrt(variance[i]);
  }

#pragma omp parallel for
  for( int n=0; n<(int)samples.size(); n++ )
  {
    for( int i=0; i<dim; i++ )   
    {
      samples[n][i] = samples[n][i]*stddev[i]+mean[i];
    }
  }

}




template void sample_mean<double>( const std::vector<double*> &samples, double *mean, int dim );
template void sample_variance<double>( const std::vector<double*> &samples, const double *mean, double *variance, int dim );
template void sample_variance<double>( const std::vector<double*> &samples, double *variance, int dim );
template void standardize<double>( std::vector<double*> &samples, double *mean, double *variance, int dim );
template void destandardize<double>( std::vector<double*> &samples, const double *mean, const double *variance, int dim );



template void sample_mean<float>( const std::vector<float*> &samples, float *mean, int dim );
template void sample_variance<float>( const std::vector<float*> &samples, const float *mean, float *variance, int dim );
template void sample_variance<float>( const std::vector<float*> &samples, float *variance, int dim );
template void standardize<float>( std::vector<float*> &samples, float *mean, float *variance, int dim );
template void destandardize<float>( std::vector<float*> &samples, const float *mean, const float *variance, int dim );
