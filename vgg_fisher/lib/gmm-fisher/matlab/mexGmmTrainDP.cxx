
/*
 * mexGmmTrainDP.cxx
 *
 * gmm = mexGmmTrainDP(n_gauss,n_dim,samples,em_param,init_mean,init_var,init_coef)
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
    // parameter validation
    if ((nrhs < 3) || !mxIsClass(prhs[0], "double") ||
            !mxIsClass(prhs[1], "double") ||
            !mxIsClass(prhs[2], "double") ||
            (nrhs >= 4 && !mxIsStruct(prhs[3])) ||
            (nrhs >= 5 && !mxIsClass(prhs[4], "double")) ||
            (nrhs >= 6 && !mxIsClass(prhs[5], "double")) ||
            (nrhs >= 7 && !mxIsClass(prhs[6], "double")) ||
            (nrhs > 5 && nrhs < 7))
        mexErrMsgTxt("Function called as : gmm = mexGmmTrainDP(n_gauss,n_dim,samples,em_param,init_mean,init_var,init_coef)");
    
    // load in data
    int n_gauss = (int)mxGetScalar(prhs[0]);
    int n_dim = (int)mxGetScalar(prhs[1]);
    
    double* samples_arr = mxGetPr(prhs[2]);
    ptrdiff_t samplesdim = mxGetM(prhs[2]); //size(samples,1) - dimension of samples
    ptrdiff_t numsamples = mxGetN(prhs[2]); //size(samples,2) - number of samples i.e. samples stored as columns
    // copy across to a vector
    std::vector<double*> samples(numsamples);
    for (int si = 0; si < numsamples; ++si) {
        samples[si] = &samples_arr[samplesdim*si];
    }
    
    // construct a c++ struct with default parameter values
    em_param gmm_params;
    
    // load in optional parameter struct if specified
    if (nrhs >= 4) {
        // first load in a list of existing fields
        int nfields = mxGetNumberOfFields(prhs[3]);
        
        // now construct a list of fieldnames
        std::set<std::string> fnames;
        for (int ifield=0; ifield < nfields; ++ifield) {
            fnames.insert(std::string(mxGetFieldNameByNumber(prhs[3],ifield)));
        }
        
        // for each field that exists in the input matlab struct, replace the default
        // value in the c++ struct with its value
        mxArray *tmp;
        if (fnames.count("max_iter") > 0) {
            tmp = mxGetField(prhs[3],0,"max_iter");
            gmm_params.max_iter = (int)mxGetScalar(tmp);
        }
        if (fnames.count("alpha") > 0) {
            tmp = mxGetField(prhs[3],0,"alpha");
            gmm_params.alpha = (float)mxGetScalar(tmp);
        }
        if (fnames.count("llh_diff_thr") > 0) {
            tmp = mxGetField(prhs[3],0,"llh_diff_thr");
            gmm_params.llh_diff_thr = (float)mxGetScalar(tmp);
        }
        if (fnames.count("min_gamma") > 0) {
            tmp = mxGetField(prhs[3],0,"min_gamma");
            gmm_params.min_gamma = (float)mxGetScalar(tmp);
        }
        if (fnames.count("variance_floor") > 0) {
            tmp = mxGetField(prhs[3],0,"variance_floor");
            gmm_params.variance_floor = (float)mxGetScalar(tmp);
        }
        if (fnames.count("variance_floor_factor") > 0) {
            tmp = mxGetField(prhs[3],0,"variance_floor_factor");
            gmm_params.variance_floor_factor = (float)mxGetScalar(tmp);
        }
    }
    
    std::vector<double*> init_mean(n_gauss), init_var(n_gauss);
    std::vector<double> init_coef;
    // load in initial mean, variance and mixing coefficients if specified
    if (nrhs >= 5) {
        double* init_mean_arr = mxGetPr(prhs[4]);
        // init_mean is a matrix of column means
        ptrdiff_t init_mean_n_dim = mxGetM(prhs[4]); //size(init_mean,1) - should always be equal to n_dim
        ptrdiff_t init_mean_n_gauss = mxGetN(prhs[4]); //size(init_mean,2) - should always be equal to n_gauss
        if ((init_mean_n_dim != n_dim) || (init_mean_n_gauss != n_gauss)) {
            mexErrMsgTxt("init_mean must be a matrix of [n_dim x n_gauss] size");
        }
        // copy across to a vector
        for (int mi = 0; mi < n_gauss; ++mi) {
            init_mean[mi] = &init_mean_arr[n_dim*mi];
        } 
    }
    if (nrhs >= 6) {
        double* init_var_arr = mxGetPr(prhs[5]);
        // init_mean is a matrix of column means
        ptrdiff_t init_var_n_dim = mxGetM(prhs[5]); //size(init_var,1) - should always be equal to n_dim
        ptrdiff_t init_var_n_gauss = mxGetN(prhs[5]); //size(init_var,2) - should always be equal to n_gauss
        if ((init_var_n_dim != n_dim) || (init_var_n_gauss != n_gauss)) {
            mexErrMsgTxt("init_var must be a matrix of [n_dim x n_gauss] size");
        }
        // copy across to a vector
        for (int mi = 0; mi < n_gauss; ++mi) {
            init_var[mi] = &init_var_arr[n_dim*mi];
        } 
    }
    if (nrhs >= 7) {
        double* init_coef_arr = mxGetPr(prhs[6]);
        // init_mean is a matrix of column means
        ptrdiff_t init_coef_n_gauss = mxGetM(prhs[6]); //size(init_coef,1) - should always be equal to n_gauss
        ptrdiff_t init_coef_cols = mxGetN(prhs[6]); //size(init_coef,2) - should always be equal to 1 (column vector)
        if ((init_coef_n_gauss != n_gauss) || (init_coef_cols != 1)) {
            mexErrMsgTxt("init_coef must be a column vector of length n_gauss");
        }
        // copy across to a vector
        init_coef = std::vector<double>(init_coef_arr, init_coef_arr + n_gauss);
    }
    
    // prepare output struct
    const char *keys[] = { "mean", "variance", "coef", "n_gauss", "n_dim", "log_likelihood" };
    plhs[0] = mxCreateStructMatrix(1, 1, 6,  keys);
    
    // call gmm library
    
    mexPrintf("Initialising GMM Builder...");
    gaussian_mixture<double> gmmproc(gmm_params, n_gauss, n_dim);
    mexPrintf("DONE\n");
    
    mexPrintf("Setting up Model...");
    if (nrhs >= 5) {
        gmmproc.set(init_mean, init_var, init_coef);
    }
    mexPrintf("DONE\n");
    
    mexPrintf("Running EM...");
    gmmproc.em(samples);
    mexPrintf("DONE\n");
    
    mexPrintf("Computing final log likelihood...");
    double llh_val = gmmproc.log_likelihood(samples);
    mexPrintf("DONE\n");
    
    // load model into output structure
    
    mxArray *mean_mat = mxCreateDoubleMatrix(n_dim,n_gauss,mxREAL);
    double* mean = mxGetPr(mean_mat);
    mxArray *variance_mat = mxCreateDoubleMatrix(n_dim,n_gauss,mxREAL);
    double* variance = mxGetPr(variance_mat);
    for (int j = 0; j < n_gauss; ++j) {
        double* componentmean = gmmproc.get_mean(j);
        double* componentvariance = gmmproc.get_variance(j);
        for (int i = 0; i < n_dim; ++i) {
            mean[i+j*n_dim] = componentmean[i];
            variance[i+j*n_dim] = componentvariance[i];
        }
    }
    mxArray *llh_mat = mxCreateDoubleMatrix(1,1,mxREAL);
    double* llh = (double*)mxGetData(llh_mat);
    (*llh) = llh_val;
    
    mxSetField(plhs[0], 0, "mean", mean_mat);
    mxSetField(plhs[0], 0, "variance", variance_mat);
    mxSetField(plhs[0], 0, "log_likelihood", llh_mat);
    
    mxArray *coef_mat = mxCreateDoubleMatrix(n_gauss,1,mxREAL);
    double* coef = mxGetPr(coef_mat);
    for (int i = 0; i < n_gauss; ++i) {
        coef[i] = gmmproc.get_mixing_coefficients(i);
    }
    mxSetField(plhs[0], 0, "coef", coef_mat);
    
    mwSize scalar_dims[]={1,1};
    mxArray *n_gauss_mat = mxCreateNumericArray(2, scalar_dims, mxDOUBLE_CLASS, mxREAL);
    double *n_gauss_arr = (double*)mxGetData(n_gauss_mat);
    n_gauss_arr[0] = (double)n_gauss;
    mxSetField(plhs[0], 0, "n_gauss", n_gauss_mat);
    
    mxArray *n_dim_mat = mxCreateNumericArray(2, scalar_dims, mxDOUBLE_CLASS, mxREAL);
    double *n_dim_arr = (double*)mxGetData(n_dim_mat);
    n_dim_arr[0] = (double)n_dim;
    mxSetField(plhs[0], 0, "n_dim", n_dim_mat);
    
}

















