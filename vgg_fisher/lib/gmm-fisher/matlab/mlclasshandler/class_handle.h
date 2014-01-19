
/// \Class    class_handle class_handle.h "class_handle.h"
///  
/// \brief
///
/// \version  1.0
/// \author   Ken Chatfield
/// \date     08/07/2011

#ifndef __CLASS_HANDLE_H
#define __CLASS_HANDLE_H

#include "mex.h"
#include <stdint.h>

#define CLASS_HANDLE_SIGNATURE 0xa5a50f0f
template<class base> class class_handle: public base
{
public:
    class_handle(): base() { signature = CLASS_HANDLE_SIGNATURE; }
    ~class_handle() { signature = 0; }
    bool isValid { return (signature == CLASS_HANDLE_SIGNATURE); }
private:
    uint32_t signature;
};

template<class base> inline mxArray *convertPtr2Mat(class_handle<base> *ptr)
{
    mxArray *out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *((uint64_t *)mxGetData(out)) = reinterpret_cast<uint64_t>(ptr);
    return out;
}

template<class base> inline class_handle<base> *convertMat2Ptr(const mxArray *in)
{
    if (mxGetNumberOfElements(in) != 1 || mxGetClassID(in) != mxUINT64_CLASS ||
            mxIsComplex(in))
        mexErrMsgTxt("Input must be a real uint64 scalar");
    class_handle<base> *ptr =
        reinterpret_cast<class_handle<base> *>(*((uint64_t *)mxGetData(in)));
    if (!ptr->isValid())
        mexErrMsgTxt("Handle not valid");
    return ptr;
}

#endif
