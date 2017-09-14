#ifndef COMMON_H
#define COMMON_H

#include "mex.h"
#include "math.h"
#include <omp.h>

typedef mxArray* M_Ptr;
typedef double   scalar_t;

#define C_CAST(CONST_M_PTR) const_cast<M_Ptr>(CONST_M_PTR)
#define dMat(N, M) mxCreateNumericMatrix(N,M, mxDOUBLE_CLASS, mxREAL)
#define iMat(N, M) mxCreateNumericMatrix(N,M, mxINT32_CLASS, mxREAL)
#define MEX_EPS 1e-12

template <typename T> T* M_Cast(M_Ptr _ptr){ return static_cast<T*> (mxGetData(_ptr)); }

extern "C" bool mxUnshareArray(M_Ptr array_ptr, bool noDeepCopy);

#endif //COMMON_H
