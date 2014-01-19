
/*
 * LLCEncodeHelper.cpp
 *
 * encoding = LLCEncodeHelper(dictWords,imWords,nearestNeighborResult,beta,outputFullMat)
 *
 * where size(dictWords) = [m n], size(imWords) = [m l] 
 * and size(nearestNeighborResult) = [n l]
 *
 * Author: Yuning Chai 
 * Modified By: Ken Chatfield
 * May 11 - added sparse encoding matrix output, type checking
 * Sep 11 - added option to either output matrix or histogram
 *
 * This is a MEX-file for MATLAB.
*/

#include "mex.h"
#include "lapack.h"
#include "blas.h"

void mult_help(const double* A, ptrdiff_t height, ptrdiff_t width, double* C){
    
    double tmp;
    for(ptrdiff_t i=0;i<width;i++){
        for(ptrdiff_t j=i;j<width;j++){
            tmp = 0;
            for(ptrdiff_t k=0;k<height;k++){
                tmp = tmp + A[i*height+k]*A[j*height+k];
            }
            C[i*width+j] = tmp;
            C[j*width+i] = tmp;
        }
    }

}


void encode(const double* dict_words, const double* im_words, 
        const double* nn_ind, ptrdiff_t dim, ptrdiff_t dict_size, 
        ptrdiff_t im_size, ptrdiff_t nn_size, double beta, 
        double* encoding_sr, mwIndex* encoding_ir, mwIndex* encoding_jc, 
        double* hist, bool output_full_mat){
    
    char upper_triangle = 'U';
    ptrdiff_t INFO;
    ptrdiff_t int_one = 1;
    double sum;

    double* z = new double[dim*nn_size];
    double* C = new double[nn_size*nn_size];
    double* b = new double[nn_size];
    ptrdiff_t* IPIV = new ptrdiff_t[nn_size];
    
    if (output_full_mat) {
        // first column with non-zero value in encoding has index 0
        encoding_jc[0] = 0;
    }

    for(mwIndex i=0; i<im_size; i++){
        
        mwIndex tmp_ind;
        for(mwIndex n=0;n<nn_size;n++){
            for(mwIndex m=0;m<dim;m++){
                tmp_ind = (mwIndex)nn_ind[i*nn_size+n];
                z[n*dim+m] = dict_words[tmp_ind*dim+m]-im_words[i*dim+m];
            }
        }
        
        mult_help(z,dim,nn_size,C);
        
//         dgemm(&transpose,&non_trans,&nn_size,&nn_size,&dim,&one,z,&nn_size,z,&dim,&zero,C,&nn_size);
        
        sum = 0;
        for(mwIndex m=0;m<nn_size;m++) sum += C[m*nn_size+m];
        sum = sum*beta;
        for(mwIndex m=0;m<nn_size;m++) C[m*nn_size+m] += sum;
            
        for(mwIndex m=0;m<nn_size;m++) b[m] = 1;
        
        dposv(&upper_triangle,&nn_size,&int_one,C,&nn_size,b,&nn_size,&INFO);
//         dgesv(&nn_size,&int_one,C,&nn_size,IPIV,b,&nn_size,&INFO);
                
        sum = 0;
        for(mwIndex m=0;m<nn_size;m++) sum += b[m];
        for(mwIndex m=0;m<nn_size;m++) b[m] /= sum;
        
        for(mwIndex m=0;m<nn_size;m++){
            tmp_ind = (mwIndex)nn_ind[i*nn_size+m];
            if (output_full_mat) {
                encoding_sr[encoding_jc[i]+m] = b[m];
                encoding_ir[encoding_jc[i]+m] = tmp_ind;
            } else {
                if(hist[tmp_ind]<b[m]) hist[tmp_ind] = b[m];
            }
        }
        if (output_full_mat) {
            encoding_jc[i+1] = encoding_jc[i] + nn_size;
        }
    }
    
    delete [] z;
    delete [] C;
    delete [] b;
    delete [] IPIV;    
    
    return;
}

//encoding = ...
//  LLCEncodeHelper(dictWords,imWords,nearestNeighborResult,beta,outputFullMat)
// * dictWords, imWords and nearestNeighborResult must be of type double
// * encoding is returned as a sparse matrix
// * hist is returned as a vector of type double
void mexFunction( int nlhs, mxArray *plhs[], 
                  int nrhs, const mxArray *prhs[])
{
    // parameter validation
    if (!(mxIsClass(prhs[0], "double")) || !(mxIsClass(prhs[1], "double")) ||
            !(mxIsClass(prhs[2], "double"))){
        mexErrMsgTxt("dictWords, imWords and nearestNeighborResult "
                "should all be of type DOUBLE.");
    }
    
    // load in data
    double* dict_words = mxGetPr(prhs[0]);
    double* im_words = mxGetPr(prhs[1]);
    double* nn_ind = mxGetPr(prhs[2]);
    double beta = mxGetScalar(prhs[3]); 
    bool output_full_mat = mxIsLogicalScalarTrue(prhs[4]);
    ptrdiff_t dim = mxGetM(prhs[0]);       //size(dictWords,1) - dims of descs
    ptrdiff_t dict_size = mxGetN(prhs[0]); //size(dictWords,2) - # visual words
    ptrdiff_t im_size = mxGetN(prhs[1]);   //size(imWords,2) - # descs in img
    ptrdiff_t nn_size = mxGetM(prhs[2]);   //number of nearest neighbours   
    
    // matlab index starts with 1 instead 0
    for(mwIndex i=0;i<nn_size*im_size;i++) nn_ind[i]--;
    
    double* encoding_sr;
    mwIndex* encoding_ir;
    mwIndex* encoding_jc;
    double* hist;
    
    if (output_full_mat) {
        /* encoding matrix should be sparse */
        plhs[0] = mxCreateSparse(dict_size,im_size,nn_size*im_size,mxREAL);

        /* encoding matrix indexed as (val,row,col) tuple */
        encoding_sr = mxGetPr(plhs[0]);
        encoding_ir = mxGetIr(plhs[0]);
        encoding_jc = mxGetJc(plhs[0]);
        
        /* create empty vars for encoding hist */
        hist = 0;
    } else {
        /* histogram should be double */
        plhs[0] = mxCreateDoubleMatrix(dict_size,1,mxREAL);
        hist = mxGetPr(plhs[0]);
        
        /* create empty vars for encoding matrix */
        encoding_sr = 0;
        encoding_ir = 0;
        encoding_jc = 0;
    }
    
    encode(dict_words, im_words, nn_ind, 
        dim, dict_size, im_size, nn_size, beta, 
        encoding_sr, encoding_ir, encoding_jc, 
        hist, output_full_mat);
    
}

















