
/*
 * mexGmmTrainSP.cxx
 *
 * gmm = mexGmmTrainSP(n_gauss,n_dim,samples,em_param,init_mean,init_var,init_coef)
 *
 * Author: Ken Chatfield
 * July 2011
 *
 * This is a MEX-file for MATLAB.
*/

#include "mex.h"
#include "../gmm.h"
#include <stdint.h>
#include <time.h>
#include <set>
#include <vector>
#include <string>

void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray *prhs[])
{        
    mwSize scalar_dims[]={1,1};
    
    if (nrhs < 2)
        mexErrMsgTxt("Function called as : gmm = mexGmmTrainSP(samples,n_gauss,[em_param],[init_mean,init_var,init_coef])");
    
    // input parameters
    const mxArray *pmxData = prhs[0];
    const mxArray *pmxNGauss = prhs[1];
    
    if (!mxIsClass(pmxData, "single"))
        mexErrMsgTxt("samples must be single");
    
    const mxArray *pmxParams = 0;
    
    if (nrhs >= 3)
    {
        pmxParams = prhs[2];
        
        if (!mxIsStruct(pmxParams))
            mexErrMsgTxt("em_param must be a structure");
    }
    
    const mxArray *pmxInitMean = 0;
    const mxArray *pmxInitVar = 0;
    const mxArray *pmxInitCoef = 0;
    
    if (nrhs >= 4)
    {
        if (nrhs == 6)
        {
            pmxInitMean = prhs[3];
            pmxInitVar = prhs[4];
            pmxInitCoef = prhs[5];
            
            if (!mxIsClass(pmxInitMean, "single"))
                mexErrMsgTxt("init_mean must be single");
            
            if (!mxIsClass(pmxInitVar, "single"))
                mexErrMsgTxt("init_var must be single");
            
            if (!mxIsClass(pmxInitCoef, "single"))
                mexErrMsgTxt("init_coef must be single");             
        }
        else
            mexErrMsgTxt("you must specify init_mean, init_var, and init_coef");
    }
        
    // number of Gaussians
    int n_gauss = (int)mxGetScalar(pmxNGauss);
        
    // data to cluster
    float* samples_arr = (float*)mxGetData(pmxData);
    
    int n_dim = (int)mxGetM(pmxData); //size(samples,1) - dimension of samples
    int numsamples = (int)mxGetN(pmxData); //size(samples,2) - number of samples i.e. samples stored as columns
    // copy across to a vector
    std::vector<float*> samples(numsamples);
    for (int si = 0; si < numsamples; ++si) {
        samples[si] = &samples_arr[n_dim*si];
    }
    
    // construct a c++ struct with default parameter values
    em_param gmm_params;
    
    // load in optional parameter struct if specified
    if (pmxParams) {
        // first load in a list of existing fields
        int nfields = mxGetNumberOfFields(pmxParams);
        
        // now construct a list of fieldnames
        std::set<std::string> fnames;
        for (int ifield=0; ifield < nfields; ++ifield) {
            fnames.insert(std::string(mxGetFieldNameByNumber(pmxParams,ifield)));
        }
        
        // for each field that exists in the input matlab struct, replace the default
        // value in the c++ struct with its value
        mxArray *tmp;
        if (fnames.count("max_iter") > 0) {
            tmp = mxGetField(pmxParams,0,"max_iter");
            gmm_params.max_iter = (int)mxGetScalar(tmp);
        }
        if (fnames.count("alpha") > 0) {
            tmp = mxGetField(pmxParams,0,"alpha");
            gmm_params.alpha = (float)mxGetScalar(tmp);
        }
        if (fnames.count("llh_diff_thr") > 0) {
            tmp = mxGetField(pmxParams,0,"llh_diff_thr");
            gmm_params.llh_diff_thr = (float)mxGetScalar(tmp);
        }
        if (fnames.count("min_gamma") > 0) {
            tmp = mxGetField(pmxParams,0,"min_gamma");
            gmm_params.min_gamma = (float)mxGetScalar(tmp);
        }
        if (fnames.count("variance_floor") > 0) {
            tmp = mxGetField(pmxParams,0,"variance_floor");
            gmm_params.variance_floor = (float)mxGetScalar(tmp);
        }
        if (fnames.count("variance_floor_factor") > 0) {
            tmp = mxGetField(pmxParams,0,"variance_floor_factor");
            gmm_params.variance_floor_factor = (float)mxGetScalar(tmp);
        }
    }
    
    std::vector<float*> init_mean(n_gauss), init_var(n_gauss);
    std::vector<float> init_coef;
    // load in initial mean, variance and mixing coefficients if specified
    if (pmxInitMean) {
        float* init_mean_arr = (float*)mxGetData(pmxInitMean);
        // init_mean is a matrix of column means
        ptrdiff_t init_mean_n_dim = mxGetM(pmxInitMean); //size(init_mean,1) - should always be equal to n_dim
        ptrdiff_t init_mean_n_gauss = mxGetN(pmxInitMean); //size(init_mean,2) - should always be equal to n_gauss
        if ((init_mean_n_dim != n_dim) || (init_mean_n_gauss != n_gauss)) {
            mexErrMsgTxt("init_mean must be a matrix of [n_dim x n_gauss] size");
        }
        // copy across to a vector
        for (int mi = 0; mi < n_gauss; ++mi) {
            init_mean[mi] = &init_mean_arr[n_dim*mi];
        } 
    }
    
    if (pmxInitVar) {
        float* init_var_arr = (float*)mxGetData(pmxInitVar);
        // init_mean is a matrix of column means
        ptrdiff_t init_var_n_dim = mxGetM(pmxInitVar); //size(init_var,1) - should always be equal to n_dim
        ptrdiff_t init_var_n_gauss = mxGetN(pmxInitVar); //size(init_var,2) - should always be equal to n_gauss
        if ((init_var_n_dim != n_dim) || (init_var_n_gauss != n_gauss)) {
            mexErrMsgTxt("init_var must be a matrix of [n_dim x n_gauss] size");
        }
        // copy across to a vector
        for (int mi = 0; mi < n_gauss; ++mi) {
            init_var[mi] = &init_var_arr[n_dim*mi];
        } 
    }
    
    if (pmxInitCoef) {
        float* init_coef_arr = (float*)mxGetData(pmxInitCoef);
        // init_mean is a matrix of column means
        ptrdiff_t init_coef_n_gauss = mxGetM(pmxInitCoef); //size(init_coef,1) - should always be equal to n_gauss
        ptrdiff_t init_coef_cols = mxGetN(pmxInitCoef); //size(init_coef,2) - should always be equal to 1 (column vector)
        if ((init_coef_n_gauss != n_gauss) || (init_coef_cols != 1)) {
            mexErrMsgTxt("init_coef must be a column vector of length n_gauss");
        }
        // copy across to a vector
        init_coef = std::vector<float>(init_coef_arr, init_coef_arr + n_gauss);
    }
    
    // prepare output struct
    const char *keys[] = { "mean", "variance", "coef", "n_gauss", "n_dim", "log_likelihood" };
    plhs[0] = mxCreateStructMatrix(1, 1, 6,  keys);
    
    // call gmm library
    
    mexPrintf("Initialising GMM Builder...");
    gaussian_mixture<float> gmmproc(gmm_params, n_gauss, n_dim);
    mexPrintf("DONE\n");
    
    if (pmxInitMean && pmxInitVar && pmxInitCoef) 
    {
        mexPrintf("Setting up the model to specified mean, var, and coef...");
        gmmproc.set(init_mean, init_var, init_coef);
    }
    else
    {
        //int seed = 42;
        int seed = time(NULL);
	
        mexPrintf("Setting up the model to random mean, var, and coef...");
        gmmproc.random_init(samples, seed);
    }
        
    mexPrintf("DONE\n");
    
    mexPrintf("Running EM...");
    gmmproc.em(samples);
    mexPrintf("DONE\n");
    
    mexPrintf("Computing final log likelihood...");
    float llh_val = gmmproc.log_likelihood(samples);
    mexPrintf("DONE\n");
    
    // load model into output structure
    mwSize meanvar_dims[]={n_dim,n_gauss};
    mxArray *mean_mat = mxCreateNumericArray(2, meanvar_dims, mxSINGLE_CLASS, mxREAL);
    float* mean = (float*)mxGetData(mean_mat);
    mxArray *variance_mat = mxCreateNumericArray(2, meanvar_dims, mxSINGLE_CLASS, mxREAL);
    float* variance = (float*)mxGetData(variance_mat);
    for (int j = 0; j < n_gauss; ++j) {
        float* componentmean = gmmproc.get_mean(j);
        float* componentvariance = gmmproc.get_variance(j);
        for (int i = 0; i < n_dim; ++i) {
            mean[i+j*n_dim] = componentmean[i];
            variance[i+j*n_dim] = componentvariance[i];
        }
    }
    mxArray *llh_mat = mxCreateNumericArray(2, scalar_dims, mxSINGLE_CLASS, mxREAL);
    float* llh = (float*)mxGetData(llh_mat);
    (*llh) = llh_val;
    
    mxSetField(plhs[0], 0, "mean", mean_mat);
    mxSetField(plhs[0], 0, "variance", variance_mat);
    mxSetField(plhs[0], 0, "log_likelihood", llh_mat);
    
    mwSize coef_dims[]={n_gauss,1};
    mxArray *coef_mat = mxCreateNumericArray(2, coef_dims, mxSINGLE_CLASS, mxREAL);
    float* coef = (float*)mxGetData(coef_mat);
    for (int i = 0; i < n_gauss; ++i) {
        coef[i] = gmmproc.get_mixing_coefficients(i);
    }
    mxSetField(plhs[0], 0, "coef", coef_mat);
    
    mxArray *n_gauss_mat = mxCreateNumericArray(2, scalar_dims, mxSINGLE_CLASS, mxREAL);
    float *n_gauss_arr = (float*)mxGetData(n_gauss_mat);
    n_gauss_arr[0] = (float)n_gauss;
    mxSetField(plhs[0], 0, "n_gauss", n_gauss_mat);
    
    mxArray *n_dim_mat = mxCreateNumericArray(2, scalar_dims, mxSINGLE_CLASS, mxREAL);
    float *n_dim_arr = (float*)mxGetData(n_dim_mat);
    n_dim_arr[0] = (float)n_dim;
    mxSetField(plhs[0], 0, "n_dim", n_dim_mat);
    
}

















