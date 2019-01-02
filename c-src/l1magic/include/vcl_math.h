#ifndef __VCL_MATH__
#define __VCL_MATH__

#ifdef __cplusplus
extern "C" {
#endif
  double vcl_logsum(const int N, const double alpha, const double *x);

  double vcl_sum(const int N, const double *x);

  double vcl_ddot(const int N, const double *x, const double *y);

  double vcl_dnorm2(const int N, const double *x);

  void vcl_daxpy(const int N, const double alpha, const double *x, double *y);

  void vcl_daxpby(const int N, const double alpha, const double *x, const double beta, double *y);

  void vcl_dxpby(const int N, const double *x, const double beta, double *y);
#ifdef __cplusplus
}
#endif

#endif
