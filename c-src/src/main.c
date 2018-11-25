/* http://cow.physics.wisc.edu/~craigm/idl/cmpfit.html */
// gcc -Icmpfit/ -Lcmpfit/ test_data.c main.c -lc -lm -lmpfit

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "test_data.h"

#include "SOS.h"
#include "mpfit.h"
#include "fit_ZPK.h"
#include "plotting.h"
// #include "plConfig.h"
// #include "plplot.h"

void print_status(int status);

void printresult(mp_result *result);

void print_frd(FRF_Data frd);



void main(){
  int i, status, len_theta;
  // double complex G;
  // double omega = 10;
  double Ts = 40e-6;
  double theta[5] = {0.678648984004927, -1.995999123944854, 0.998919373769841, -1.996142567427762, 0.998941969049303};
  len_theta = 5;
  FRF_Data frd;


  double complex *resp;
  double complex *resp_fit;
  resp = malloc( LEN_omegas * sizeof( double complex));
  resp_fit = malloc( LEN_omegas * sizeof( double complex));

  for (i=0; i < LEN_omegas; i++){
    resp[i] = resp_real[i] + resp_imag[i]*I;
  }

  frd.N = LEN_omegas;
  frd.Resp = resp;
  frd.omegas = omegas;
  frd.Ts = Ts;

  mp_result *result;
  result = malloc(sizeof(mp_result));

  status = mpfit(&absfunc, LEN_omegas, len_theta, theta, 0, 0, &frd, result);

  printresult(result);

  print_status(status);



  for (i=0; i<len_theta; i++){
    printf("Theta[%d] = %.12f\n", i, theta[i]);
  }
  get_frf(LEN_omegas, frd.omegas, 5, theta, Ts, resp_fit);
  // for (i=0; i<LEN_omegas; i++){
  //   printf("resp[i] = %.6f + %.6fj\n", creal(resp_fit[i]), cimag(resp_fit[i]));
  // }

  bode(LEN_omegas, omegas, frd.Resp, resp_fit);


  free(resp);
  free(result);
  }


/*Print out the frf, for debugging */
void print_frd(FRF_Data frd){
  int i;
  printf("omega   |    resp    \n");

  for (i=0; i<frd.N; i++){
    printf("%.6f     %.6f+ %.6f\n", frd.omegas[i], creal(frd.Resp[i]), cimag(frd.Resp[i]));
  }
}

/* Simple routine to print the fit results */
void printresult(mp_result *result)
{
int i;

printf("  CHI-SQUARE = %f    (%d DOF)\n",
         result->bestnorm, result->nfunc-result->nfree);
printf("        NPAR = %d\n", result->npar);
printf("       NFREE = %d\n", result->nfree);
printf("     NPEGGED = %d\n", result->npegged);
printf("     NITER = %d\n", result->niter);
printf("      NFEV = %d\n", result->nfev);
printf("\n");

}


/* Routine to parse the error code and print it */
void print_status(int status){
  printf("Status code: %d\n", status);

  switch (status){
  case MP_ERR_INPUT     :
    printf("General input parameter error              \n");
    break;
  case MP_ERR_NAN       :
    printf("User function produced non-finite values   \n");
    break;
  case MP_ERR_FUNC      :
    printf("* No user function was supplied            \n");
    break;
  case MP_ERR_NPOINTS   :
    printf("* No user data points were supplied        \n");
    break;
  case MP_ERR_NFREE     :
    printf("* No free parameters                       \n");
    break;
  case MP_ERR_MEMORY    :
    printf("* Memory allocation error                  \n");
    break;
  case MP_ERR_INITBOUNDS:
    printf("* Initial values inconsistent w constraints* \n");
    break;
  case MP_ERR_BOUNDS    :
    printf("* Initial constraints inconsistent         \n");
    break;
  case MP_ERR_PARAM     :
    printf("* General input parameter error            \n");
    break;
  case MP_ERR_DOF       :
    printf("* Not enough degrees of freedom            \n");
    break;
  case    MP_OK_CHI   :
    printf("* Convergence in chi-square value           \n");
    break;
  case    MP_OK_PAR   :
    printf("* Convergence in parameter value            \n");
    break;
  case    MP_OK_BOTH  :
    printf("* Both MP_OK_PAR and MP_OK_CHI hold         \n");
    break;
  case    MP_OK_DIR   :
      printf("* Convergence in orthogonality              \n");
    break;
  case    MP_MAXITER  :
    printf("* Maximum number of iterations reached      \n");
    break;
  case    MP_FTOL     :
    printf("* ftol is too small; no further improvement \n");
    break;
  case    MP_XTOL     :
    printf("* xtol is too small; no further improvement \n");
    break;
  case    MP_GTOL     :
    printf("* gtol is too small; no further improvement \n");
    break;
  }
}
