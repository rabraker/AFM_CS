#ifndef _SOS_
#define _SOS_
#include "config.h"

/* Use this to wrap function arguments that arent used, but are required by the callers
 callback spec.
*/
#define UNUSED(x) (void)(x)

#include <complex.h>

DLL_PUBLIC void get_frf(int N_w, double *omega, int m, double *theta, double Ts, double complex *Resp);

double complex freqresp_SOS(double omega, double *theta, double Ts);

typedef struct FRF_Data_ {
  double *omegas;
  double complex *Resp;
  double Ts;
  int N;
} FRF_Data;

#endif
