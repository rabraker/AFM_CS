#include <complex.h>
#include <stdlib.h>
#include <math.h>

#include "mpfit.h"
#include "fit_ZPK.h"
#include "SOS.h"
#include "plotting.h"



static int absfunc(int m, int n, double *theta, double *resid, double **dvec, void *user_data)
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


/* Return 1 if roots of a2*z^2 + a1*z + a0 are within the
   unit circle.

   Return 0 otherwise.
*/
int poly2_is_z_stable(double a, double b, double c){

  double complex root1;
  double complex root2;
  double complex a_c = a + 0*I;
  double complex b_c = b + 0*I;
  double complex c_c = c + 0*I;

  root1 =( -b_c + csqrt(b_c * b_c - 4*a_c * c_c)) / (2 * a_c);
  root2 =( -b_c - csqrt(b_c * b_c - 4*a_c * c_c)) / (2 * a_c);

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
DLL_PUBLIC int fit_sos_st(int N_omegas, double *omegas, double *resp_real,
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
  mp_par pars[N_theta];
  double ub[5] = {20, 2, 1, 2, 1};
  double lb[5] = {-20, -2, -1, -2, -1};

  for (i=0; i<N_theta; i++){
    pars[i].fixed = 0;

    pars[i].limited[0] = 1;
    pars[i].limited[1] = 1;
    pars[i].limits[0] = lb[i];
    pars[i].limits[1] = ub[i];

    pars[i].parname = 0;
    pars[i].step = 0;
    pars[i].relstep = 0;
    pars[i].side = 0;
    pars[i].deriv_debug = 0;

  }
  status = mpfit(&absfunc, N_omegas, N_theta, theta, pars, 0, &frd, mpres);

  fit_sos_result->stable_num =  poly2_is_z_stable(1.0, theta[1], theta[2]);  // Numerator is stable?
  fit_sos_result->stable_den =  poly2_is_z_stable(1.0, theta[3], theta[4]); // Denominator is stable?
  fit_sos_result->bestnorm = mpres->bestnorm;
  fit_sos_result->orignorm = mpres->orignorm;
  fit_sos_result->niter = mpres->niter;
  fit_sos_result->nfev = mpres->nfev;
  fit_sos_result->mpfit_status = mpres->status;
  fit_sos_result->npegged = mpres->npegged;
  fit_sos_result->nfunc = mpres->nfunc;

  // free(mpres);
  free(resp_fit);
  free(resp);

  return status;

}
