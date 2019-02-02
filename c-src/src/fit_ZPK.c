#include <complex.h>
#include <stdlib.h>
#include <math.h>

#include "mpfit.h"
// #include "fit_ZPK.h"
#include "SOS.h"
#include "plotting.h"


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

  UNUSED(n);
  UNUSED(dvec);
  int i;
  FRF_Data *model_data = ( FRF_Data *) user_data;
  double Ts = model_data->Ts;

  double *omegas;
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



int fit_sos(int N_omegas, double *omegas, double *resp_real,
             double *resp_imag, int N_theta, double *theta, double result[9]){
  int status=0, i=0;

  mp_result *mpres = malloc(sizeof(mp_result));
  double complex *resp_fit = malloc(N_omegas * sizeof( double complex));
  double complex *resp = malloc(N_omegas * sizeof( double complex));

  mpres->bestnorm=0;
  mpres->orignorm=0;
  mpres->niter=0;
  mpres->nfev=0;
  mpres->status=0;
  mpres->npar=0;
  mpres->nfree=0;
  mpres->npegged=0;
  mpres->nfunc=0;
  mpres->resid=0;
  mpres->xerror=0;
  mpres->covar=0;


  for (i=0; i < N_omegas; i++){
    resp[i] = resp_real[i] + resp_imag[i]*I;
  }


  FRF_Data frd = {.N = N_omegas,
                  .Resp=resp,
                  .omegas = omegas,
                  .Ts=40e-6};

  status = mpfit(&absfunc, N_omegas, N_theta, theta, 0, 0, &frd, mpres);

  result[0] = mpres->bestnorm;
  result[1] = mpres->orignorm;
  result[2] = (double) mpres->niter;
  result[3] = (double) mpres->nfev;
  result[4] = (double) mpres->status;
  result[5] = (double) mpres->npar;
  result[6] = (double) mpres->nfree;
  result[7] = (double) mpres->npegged;
  result[8] = (double) mpres->nfunc;

  free(mpres);
  free(resp_fit);
  free(resp);

  return status;

}



