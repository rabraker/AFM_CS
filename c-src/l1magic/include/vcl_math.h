#ifndef __VCL_MATH__
#define __VCL_MATH__
#include "l1c_common.h"

#ifdef __cplusplus
extern "C" {
#endif
  double vcl_logsum(const l1c_int N, const double alpha, const double *x);

  double vcl_sum(const l1c_int N, const double *x);

  void vcl_dxMy_pz(const int N, const double *x, const double *y, double *z);
#ifdef __cplusplus
}
#endif

#endif
