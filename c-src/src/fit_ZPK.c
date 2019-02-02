#include <complex.h>
#include <stdlib.h>
#include <math.h>

#include "mpfit.h"
#include "fit_ZPK.h"
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


/* Return 1 of roots of a2*z^2 + a1*z + a0 are within the
 unit circle.

 Return 0 otherwise.
*/
int poly2_is_z_stable(double a, double b, double c){

  double root1 = 0.0, root2 = 0.0;

  root1 =( -b + csqrt(b * b - 4*a * c)) / (2 * a);
  root2 =( -b - csqrt(b * b - 4*a * c)) / (2 * a);

  if( ( cabs(root1) < 1.0) && ( cabs(root2) < 1) ){
    return 1;
  }

  return 0;

}


/*

  Given an initial guess in parameter vector theta, fits the
  a SOS transfer function of the form

  k = theta[0]
  b1 = theta[1]
  b0 = theta[2]
  a1 = theta[3]
  a0 = theta[4]

           z^2 + b1*z + b0
  G(z) =k ----------------------
        z^2 + a1*z + a0

  to the FRF data contained in resp_real + j*resp_imag, defined
  at frequencies omega.

*/
int fit_sos(int N_omegas, double *omegas, double *resp_real,
             double *resp_imag, int N_theta, double *theta, double result[11]){
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

  result[9] = poly2_is_z_stable(1.0, theta[1], theta[2]); // Numerator is stable?
  result[10] = poly2_is_z_stable(1.0, theta[3], theta[4]); // Denominator is stable?
  free(mpres);
  free(resp_fit);
  free(resp);

  return status;

}




int fit_sos_st(int N_omegas, double *omegas, double *resp_real,
               double *resp_imag, int N_theta, double *theta,  struct FitSOSResult *fit_sos_result){
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

  fit_sos_result->stable_num =  poly2_is_z_stable(1.0, theta[1], theta[2]);  // Numerator is stable?
  fit_sos_result->stable_den =  poly2_is_z_stable(1.0, theta[3], theta[4]); // Denominator is stable?
  fit_sos_result->bestnorm = mpres->bestnorm;
  fit_sos_result->orignorm = mpres->orignorm;
  fit_sos_result->niter = mpres->niter;
  fit_sos_result->nfev = mpres->nfev;
  fit_sos_result->mpfit_status = mpres->status;
  fit_sos_result->npegged = mpres->npegged;
  fit_sos_result->nfunc = mpres->nfunc;

  // result[9] = poly2_is_z_stable(1.0, theta[1], theta[2]); // Numerator is stable?
  // result[10] = poly2_is_z_stable(1.0, theta[3], theta[4]); // Denominator is stable?
  // free(mpres);
  free(resp_fit);
  free(resp);

  return status;

}
