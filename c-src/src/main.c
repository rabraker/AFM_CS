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

#ifndef __MINGW32__
#include "plotting.h"
#endif


void print_status(int status);

void printresult(mp_result *result);

void print_frd(FRF_Data frd);


int main(){
#define len_theta  5

  int i=0, status=0;

  double Ts = 40e-6;
  double result[9];
  double theta[len_theta] = {0.678648984004927, -1.995999123944854, 0.998919373769841,
                             -1.996142567427762, 0.998941969049303};

  double complex *resp = malloc( LEN_omegas * sizeof( double complex));
  double complex *resp_fit =malloc( LEN_omegas * sizeof( double complex));


  for (i=0; i< LEN_omegas; i++){
    resp[i] = resp_real[i] + resp_imag[i]*I;
  }
  status = fit_sos(LEN_omegas, omegas, resp_real, resp_imag, len_theta, theta, result);


  get_frf(LEN_omegas, omegas, 5, theta, Ts, resp_fit);
  // printresult(result);

  print_status(status);



  for (i=0; i<len_theta; i++){
    printf("Theta[%d] = %.12f\n", i, theta[i]);
  }


  printf("idx  |   resp_real   | resp_fit_real   | resp_imag   |resp_fit_imag\n");
  for (i=0; i<LEN_omegas; i++){
    printf("%d    %.6f    %.6f    %.6fj    %.6fj\n", i, resp_real[i], creal(resp_fit[i]),
           resp_imag[i], cimag(resp_fit[i]));
  }



#ifndef __MINGW32__
  bode(LEN_omegas, omegas, resp, resp_fit);
#endif

  free(resp);
  free(resp_fit);

  return 0;
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
