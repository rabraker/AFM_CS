#include <math.h>
#include <complex.h>
#include "plplot.h"
#include "plConfig.h"

#include "plotting.h"

#define PI 3.141592653589793
/* Forward declaration of functions in this module */
static double min_vec(int N, double *x);
static double max_vec(int N, double *x);

// void bode( int N, double *omegas, double complex *resp, double complex *resp_fit);


/* Plot the magnitude bode */
void bode( int N, double *omegas, double complex *resp, double complex *resp_fit)
{

  PLFLT *x, *y, *y_fit;
  x = malloc( N*sizeof(PLFLT));
  y_fit = malloc( N*sizeof(PLFLT));
  y = malloc( N*sizeof(PLFLT));
  int   i;

  PLFLT xmin, xmax, ymin, ymax;

  // Prepare data to be plotted.

  for ( i = 0; i < N; i++ )
    {
      x[i] = log10(omegas[i] / (2 * PI));
      y[i] = 20*log10(cabs(resp[i]));
      y_fit[i] = 20*log10(cabs(resp_fit[i]));
    }

  xmin = x[0];
  xmax = x[N-1];
  ymin = min_vec(N, (double *) y)*1.5;
  ymax = max_vec(N, (double *) y)*1.5;

  printf("xmin: %.6f, xmax: %.6f\n", xmin, xmax);
  //


  plsdev("xcairo");
  // Initialize plplot
  plscol0(15, 0, 0, 0);
  plscolbg(255, 255, 255);

  plinit();
  plfont( 2 );

  plcol0(15);

  // plenv(xmin, xmax, ymin, ymax, 0, 13);
  pladv( 0 );

  plvpor( 0.1, 0.95, 0.1, 0.9 ); //Basically, the axes limits.
  plwind( xmin, xmax, ymin, ymax);
  // // plwind( xmin, xmax, ymin, ymax );
  // // Create a labelled box to hold the plot.
  // // plenv(xmin, xmax, ymin, ymax, 0, 0 );
  // plcol0(15);
  // plbox("bcgl", .0, 0, "bcg", 0.0, 0 );

  double m, b;
  m = 1 / log10( omegas[N-1]/omegas[0]);
  b = - m * log10( omegas[0]/(2*PI) );
  char buf[20];
  // plschr(0, .75);

  double x_space = 10;
  double x_minor_space = 2;
  double x_minor_pt;
  double xtick_ht = (ymax - ymin)/50;
  double x_pt = pow(10, xmin);
  while (x_pt < pow(10, xmax)){
    // Reset linewidth, color, linestyle.
    plwidth(0.95);
    plcol0(7);
    pllsty(1);
    //Draw xticks
    pljoin(log10(x_pt), ymin+xtick_ht, log10(x_pt), ymax-xtick_ht);

    // Minor ticks. Starting point is major tick.
    x_minor_pt = x_pt;

    pllsty(2);
    while (x_minor_pt < x_pt + x_space){
      if (x_minor_pt == x_pt){
        pljoin(log10(x_minor_pt), ymin+xtick_ht, log10(x_minor_pt), ymax-xtick_ht);
      }
      else{
        pljoin(log10(x_minor_pt), ymin, log10(x_minor_pt), ymax);
      }
      x_minor_pt = x_minor_pt + x_minor_space;
    }
    plwidth(1.5);
    plcol0(15);
    pllsty(1);
    pljoin(log10(x_pt), ymin, log10(x_pt), ymin + xtick_ht);
    pljoin(log10(x_pt), ymax, log10(x_pt), ymax - xtick_ht);

    snprintf(buf, 6, "%.0f", x_pt);
    plmtex("b", 1., m*log10(x_pt)+b, 0.5, buf);

    x_pt = x_pt + x_space; //Update for next round
  }

  plschr(0, 1);

  pllsty(1);
  plcol0(7);
  plbox("", 0, 0, "g", 2.0, 0 );
  plcol0(15);
  plbox("bcn", 0, 0, "bnstv", 2.0, 0 );


  plcol0(15);
  pllab( "frequency", "Mag [dB]", "Bode Magnitude" );

  plcol0( 1 );
  // Plot the data that was prepared above.
  plwidth(2);
  pllsty(1);
  // plcol0(9);
  plline( N, x, y );

  plwidth(1.5);
  plcol0(3);
  pllsty(2);
  plline( N, x, y_fit );

  // Close PLplot library
  plend();

}




/*return the min of a vector */
static double min_vec(int N, double *x){
  int i;
  double min_x = x[0];
  for (i=1; i<N; i++){
    min_x = fmin(min_x, x[i]);
  }
  return min_x;
}

/*return the max of a vector */
static double max_vec(int N, double *x){
  int i;
  double max_x = x[0];
  for (i=1; i<N; i++){
    max_x = fmax(max_x, x[i]);
  }
  return max_x;
}
