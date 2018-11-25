#include <complex.h>
#include <math.h>

#include "mpfit.h"
// #include "fit_ZPK.h"
#include "SOS.h"



int absfunc(int m, int n, double *theta, double *resid, double **dvec, void *user_data)
{
  /*
    int m         - number of data points
    int n         - number of parameters (length(theta) )
    double *theta - array of n parameters
    double *resid - array of m residuals to be returned by myfunct()
    double **derivs - used for user-computed derivatives.
    (= 0  when automatic finite differences are computed)
  */


  int i;
  FRF_Data *model_data = ( FRF_Data *) user_data;
  double Ts = model_data->Ts;

  double *omegas, *ey, f;
  double complex *Resp;
  double complex G_est;
  omegas = model_data->omegas;
  Resp = model_data->Resp;

  for (i=0; i<m; i++) {
    G_est = freqresp_SOS(omegas[i], theta, Ts);
    resid[i] = log(cabs(G_est)) - log(cabs(Resp[i]));
  }

  return 0;
}








// typedef struct SoS_{
//   // assume z^2 + a1*z a2. Do not store first coef.
//   double Den[2];
//   double Num[2];
//   double Ts;
//   double k;

// } SoS;
// double complex eval_frf(double omega, SoS sos){

//   double complex z1;
//   double complex z2;
//   double complex den;
//   double complex num;
//   double complex G;
//   // e^{j * Ts * w}
//   z1 = cexp( I * omega * sos.Ts);
//   z2 = cpow(z1, 2);

//   num = z2 + z1 * sos.Num[0] + sos.Num[1];
//   den = z2 + z1 * sos.Den[0] + sos.Den[1];

//   G = sos.k * num / den;
//   return G;
// }
