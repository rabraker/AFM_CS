#ifndef _FIT_ZPK_
#define _FIT_ZPK_

#include "config.h"

#include <complex.h>
#include "mpfit.h"

struct  FitSOSResult {
  int stable_num;
  int stable_den;
  double bestnorm;
  double orignorm;
  int niter;
  int nfev;
  int mpfit_status;
  int nfree;
  int npegged;
  int nfunc;

};




DLL_PUBLIC int fit_sos_st(int N_omegas, double *omegas, double *resp_real,
                          double *resp_imag, int N_theta, double *theta,
                          struct  FitSOSResult *fit_sos_result);



#endif
