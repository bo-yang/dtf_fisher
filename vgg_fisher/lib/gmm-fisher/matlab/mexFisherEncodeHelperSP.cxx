
/*
 * mexFisherEncodeHelper.cxx
 *
 * handle = mexFisherEncodeHelper('init',gmm_mat,<fisher_params>)
 * FV = mexFisherEncodeHelper('encode',handle,x,<weights>)
 * mexFisherEncodeHelper('clear',handle)
 *
 * Author: Ken Chatfield
 * July 2011
 *
 * This is a MEX-file for MATLAB.
*/

#include "mex.h"
#include "../fisher.h"
#include "../gmm.h"
#include "mlclasshandler/fisher_handle.h"
#include <set>
#include <vector>
#include <string>
#include <string.h>

void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray *prhs[])
{    
    // first parameter validation
    if (nrhs < 1)
        mexErrMsgTxt("No command name passed");
    
    // get string name of command from first input parameter
    char *command = mxArrayToString(prhs[0]);
    
    if (strcmp(command,"init") == 0) {
        // subfunction parameter validation
        if (nrhs < 2)
            mexErrMsgTxt("At least a GMM structure must be passed");
        if (!mxIsStruct(prhs[1]))
            mexErrMsgTxt("second parameter should be a GMM structure");
        // load in GMM fields (to mxArrays)
        mxArray *mean_mat = mxGetField(prhs[1], 0, "mean");
        mxArray *variance_mat = mxGetField(prhs[1], 0, "variance");
        mxArray *coef_mat = mxGetField(prhs[1], 0, "coef");
        mxArray *n_gauss_mat = mxGetField(prhs[1], 0, "n_gauss");
        mxArray *n_dim_mat = mxGetField(prhs[1], 0, "n_dim");
        // validate GMM field names
        if ((mean_mat == 0) || (variance_mat == 0) || (coef_mat == 0) ||
                (n_gauss_mat == 0) || (n_dim_mat == 0))
            mexErrMsgTxt("Invalid GMM Structure - not all fields exist");
        // load in data from GMM fields
        float *mean_arr = (float*)mxGetData(mean_mat);
        float *variance_arr = (float*)mxGetData(variance_mat);
        float *coef_arr = (float*)mxGetData(coef_mat);
        int n_gauss = (int)mxGetScalar(n_gauss_mat);
        int n_dim = (int)mxGetScalar(n_dim_mat);
        //validate GMM field dimensions
        size_t mean_mat_m = mxGetM(mean_mat);
        size_t mean_mat_n = mxGetN(mean_mat);
        size_t variance_mat_m = mxGetM(variance_mat);
        size_t variance_mat_n = mxGetN(variance_mat);
        size_t coef_mat_m = mxGetM(coef_mat);
        size_t coef_mat_n = mxGetN(coef_mat);
        if ((mean_mat_m != n_dim) || (variance_mat_m != n_dim) ||
                (coef_mat_m != n_gauss) || (mean_mat_n != n_gauss) ||
                (variance_mat_n != n_gauss) || (coef_mat_n != 1)) {
            mexPrintf(" coef_mat: [%d %d]\n mean_mat: [%d %d]\n variance_mat: [%d %d]\n n_gauss: %d\n n_dim: %d\n",
                      (int)coef_mat_m, (int)coef_mat_n, (int)mean_mat_m, (int)mean_mat_n,
                      (int)variance_mat_m, (int)variance_mat_n, (int)n_gauss, (int)n_dim);
            mexErrMsgTxt("Invalid GMM Structure - dimension error");
        }
        // convert GMM field arrays to vectors
        std::vector<float*> mean(n_gauss);
        std::vector<float*> variance(n_gauss);
        for (int j = 0; j < n_gauss; ++j) {
            mean[j] = &mean_arr[j*n_dim];
            variance[j] = &variance_arr[j*n_dim];
        }
        std::vector<float> coef(coef_arr, coef_arr + n_gauss);
        // prepare a GMM model with data from the structure
        gaussian_mixture<float> *gmmproc = new gaussian_mixture<float>(n_gauss,n_dim);
        gmmproc->set(mean, variance, coef);
        
        // construct a c++ struct with default parameter values
        fisher_param fisher_encoder_params;
        
        // check if a set of fisher parameters were passed as the third
        // parameter - if present will be used to initialise the fisher
        // encoding class 'fisher'
        if (nrhs >= 3) {
            if (!mxIsStruct(prhs[2]))
                mexErrMsgTxt("Third parameter, if passed, should be a structure of fisher encoding options");
            // first load in a list of existing fields
            int nfields = mxGetNumberOfFields(prhs[2]);
            // now construct a list of fieldnames
            std::set<std::string> fnames;
            for (int ifield=0; ifield < nfields; ++ifield) {
                fnames.insert(std::string(mxGetFieldNameByNumber(prhs[2],ifield)));
            }
            // for each field that exists in the input matlab struct, replace the default
            // value in the c++ struct with its value
            mxArray *tmp;
            if (fnames.count("grad_weights") > 0) {
                tmp = mxGetField(prhs[2],0,"grad_weights");
                fisher_encoder_params.grad_weights = (bool)mxGetScalar(tmp);
            }
            if (fnames.count("grad_means") > 0) {
                tmp = mxGetField(prhs[2],0,"grad_means");
                fisher_encoder_params.grad_means = (bool)mxGetScalar(tmp);
            }
            if (fnames.count("grad_variances") > 0) {
                tmp = mxGetField(prhs[2],0,"grad_variances");
                fisher_encoder_params.grad_variances = (bool)mxGetScalar(tmp);
            }
            if (fnames.count("alpha") > 0) {
                tmp = mxGetField(prhs[2],0,"alpha");
                fisher_encoder_params.alpha = (float)mxGetScalar(tmp);
            }
            if (fnames.count("pnorm") > 0) {
                tmp = mxGetField(prhs[2],0,"pnorm");
                fisher_encoder_params.pnorm = (float)mxGetScalar(tmp);
            }
            if (fnames.count("gamma_eps") > 0) {
                tmp = mxGetField(prhs[2],0,"gamma_eps");
                fisher_encoder_params.gamma_eps = (float)mxGetScalar(tmp);
            }
        }            
        
        fisher_handle<float> *fisher_encoder =
        new fisher_handle<float>(*gmmproc,fisher_encoder_params);
        // initialise encoder with a GMM model (vocabulary)
        fisher_encoder->set_model(*gmmproc);
        // return handle to fisher encoder class instance
        plhs[0] = convertPtr2Mat<float>(fisher_encoder);
    } else if (strcmp(command,"encode") == 0) {
        // subfunction parameter validation
        if (nrhs < 3)
            mexErrMsgTxt("At least a handle and matrix of vectors to encode must be passed");
        if (nlhs < 1)
            mexErrMsgTxt("No output array specified");
        // get encoder from handle
        fisher_handle<float> *fisher_encoder =
                convertMat2Ptr<float>(prhs[1]);
        // get matrix of vectors to encode
        float *x_arr = (float*)mxGetData(prhs[2]);
        size_t x_m = mxGetM(prhs[2]);
        size_t x_n = mxGetN(prhs[2]);
        // convert input vectors to c++ std::vector format
        std::vector<float*> x(x_n);
        for (int j = 0; j < x_n; ++j) {
            x[j] = &x_arr[j*x_m];
        }
        // load in weights if specified
        std::vector<float> wgh;
        if (nrhs >= 4) {
            float *wgh_arr = (float*)mxGetData(prhs[3]);
            size_t wgh_m = mxGetM(prhs[2]);
            size_t wgh_n = mxGetN(prhs[2]);
            //validate dimensions
            if ((wgh_n != 1)  || (wgh_m != x_n))
                mexErrMsgTxt("Weight vector of mismatched dimensions");
            // convert to std::vector
            wgh = std::vector<float>(wgh_arr, wgh_arr + x_n);
        }
        // do encoding
        mwSize fk_dims[]={fisher_encoder->dim(),1};
        plhs[0] = mxCreateNumericArray(2, fk_dims, mxSINGLE_CLASS, mxREAL);
        float* fk = (float*)mxGetData(plhs[0]);
        if (nrhs < 4) {
            fisher_encoder->compute(x, fk);
        } else {
            fisher_encoder->compute(x, wgh, fk);
        }
    } else if (strcmp(command,"getdim") == 0) {
        // subfunction parameter validation
        if (nrhs < 2)
            mexErrMsgTxt("At least a handle must be passed");
        if (nlhs < 1)
            mexErrMsgTxt("No output array specified");
        // get encoder from handle
        fisher_handle<float> *fisher_encoder =
                convertMat2Ptr<float>(prhs[1]);
        // get dimensionality of fisher vector and move to output array
        mwSize float_dims[]={1,1};
        plhs[0] = mxCreateNumericArray(2, float_dims, mxSINGLE_CLASS, mxREAL);
        float* fkdim = (float*)mxGetData(plhs[0]);
        (*fkdim) = (float)fisher_encoder->dim();
    } else if (strcmp(command,"clear") == 0) {
        // get encoder from handle
        fisher_handle<float> *fisher_encoder =
                convertMat2Ptr<float>(prhs[1]);
        gaussian_mixture<float>* gmmproc = fisher_encoder->getGmmPtr();
        if (gmmproc)
            delete gmmproc;
        delete fisher_encoder;
    } else {
        mexErrMsgTxt("Command name not recognised");
    }
    
}

















