
/// \Class    fisher_handle fisher_handle.h "fisher_handle.h"
///  
/// \brief
///
/// \version  1.0
/// \author   Ken Chatfield
/// \date     08/07/2011

#ifndef __FISHER_HANDLE_H
#define __FISHER_HANDLE_H

#include "../../fisher.h"
#include "../../gmm.h"
#include "mex.h"
#include <stdint.h>

#define CLASS_HANDLE_SIGNATURE 0xa5a50f0f
template<class T> class fisher_handle: public fisher<T>
{
public:
    fisher_handle(gaussian_mixture<T> &gmm, fisher_param params): fisher<T>(params) { signature = CLASS_HANDLE_SIGNATURE; gmmproc = &gmm; }
    ~fisher_handle() { signature = 0; }
    bool isValid() { return (signature == CLASS_HANDLE_SIGNATURE); }
    gaussian_mixture<T>* getGmmPtr() { return gmmproc; }
private:
    uint32_t signature;
    gaussian_mixture<T> *gmmproc;
};

template<class T> inline mxArray *convertPtr2Mat(fisher_handle<T> *ptr)
{
    mxArray *out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *((uint64_t *)mxGetData(out)) = reinterpret_cast<uint64_t>(ptr);
    return out;
}

template<class T> inline fisher_handle<T> *convertMat2Ptr(const mxArray *in)
{
    if (mxGetNumberOfElements(in) != 1 || mxGetClassID(in) != mxUINT64_CLASS ||
            mxIsComplex(in))
        mexErrMsgTxt("Input must be a real uint64 scalar");
    fisher_handle<T> *ptr =
            reinterpret_cast<fisher_handle<T> *>(*((uint64_t *)mxGetData(in)));
    if (!ptr->isValid())
        mexErrMsgTxt("Handle not valid");
    return ptr;
}

#endif
