#include <math.h>
#include <complex.h>

#include "SOS.h"



/*Evaluate the frequency response at a single omega. Assumes theta  has size 5 */
double complex freqresp_SOS(double omega, double *theta, double Ts){

  double complex z1;
  double complex z2;
  double complex den;
  double complex num;
  double complex G;


  // e^{j * Ts * w}
  z1 = cexp( I * omega * Ts);
  z2 = cpow(z1, 2);

  // Theta shall look like:
  // [K, b1, b2, a1, a2]
  num = z2 + z1 *theta[1]  + theta[2];
  den = z2 + z1 * theta[3] + theta[4];

  G = theta[0] * num / den;
  return G;

}


/*Get the frequency response of a grid of omegas. */
void get_frf(int N_w, double *omega, int m, double *theta, double Ts, double complex *Resp){
  int i;
  for (i=0; i<N_w; i++){
    Resp[i] =  freqresp_SOS(omega[i], theta,  Ts);
  }
}
