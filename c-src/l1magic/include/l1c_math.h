#ifndef __L1C_MATH__
#define __L1C_MATH__

#include "l1c_common.h"


void l1c_dxmuly_z(l1c_int N, double *x, double *y, double *z);

void l1c_daxpy_z(l1c_int N, double alpha, double *x, double *y, double *z);

double l1c_dnorm1(l1c_int N, double *x);

double l1c_dsum(l1c_int N, double *x);

double l1c_dlogsum(l1c_int N,  double alpha, double *x);



#endif
