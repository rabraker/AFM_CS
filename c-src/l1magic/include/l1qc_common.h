#ifndef _L1QC_COMMON_
#define _L1QC_COMMON_

extern int N_aligned;


#define ALIGNMENT_BYTES 64

#ifdef __MATLAB__
#include "mex.h"
#define PRINT(...) mexPrintf(__VA_ARGS__);
#else
#include <stdio.h>
#define PRINT(...) printf(__VA_ARGS__);
#endif


#ifdef _USEMKL_
#include "mkl.h"
#define free_double(x) mkl_free(x)
#define malloc_double(N) (double*)mkl_malloc(N*sizeof(double), 64)
#else
#include "cblas.h"

#define free_double(x) free(x)
#define malloc_double(N) calloc(N, sizeof(double))
#endif

#define min(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
    __typeof__ (b) _b = (b);                    \
    _a < _b ? _a : _b; })

#define max(a,b)                                \
  ({ __typeof__ (a) _a = (a);                   \
    __typeof__ (b) _b = (b);                    \
    _a > _b ? _a : _b; })


#endif
