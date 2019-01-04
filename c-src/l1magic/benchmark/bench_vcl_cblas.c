#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "omp.h"

#include "mkl.h"
#include "vcl_math.h"
#include "l1qc_common.h"

#include <stdint.h>
#include <sys/timeb.h>

// needs -lrt (real-time lib)
// 1970-01-01 epoch UTC time, 1 mcs resolution (divide by 1M to get time_t)
// uint64_t ClockGetTime()
// {
//   struct timespec ts;
//   clock_gettime(CLOCK_REALTIME, &ts);
//   return (uint64_t)ts.tv_sec * 1000000LL + (long)ts.tv_nsec / 1000LL;
//   //return (uint64_t)ts.tv_nsec;
// }

double cddot(int N, double *x, double *y){
  return cblas_ddot(N, x, 1, y, 1);
}

inline double dot(const int N, double * restrict x, double * restrict y){
  double *x_ = __builtin_assume_aligned(x, 64);
  double *y_ = __builtin_assume_aligned(y, 64);
  double total = 0.0;

#pragma omp parallel
  {
#pragma omp for reduction(+:total)
  for (int i =0; i<N; i++){
    total += x_[i] * y_[i];
  }
  }
  return total;
}

// void axpy(const int N, const double alpha, double * restrict x, double * restrict y){
//   double *x_ = __builtin_assume_aligned(x, 64);
//   double *y_ = __builtin_assume_aligned(y, 64);

//   for (int i =0; i<N; i++){
//     y[i] = x_[i] * alpha + y_[i];
//   }

// }



// double norm2(const int N, double * restrict x){
//   double *x_ = __builtin_assume_aligned(x, 64);

//   double total = 0.0;
//   for (int i =0; i<N; i++){
//     total += x_[i] * x[i];
//   }

//   return total;
// }

void load_xy(int N, double *x, double *y){
  for (int i=0; i<N; i++){
    x[i] = ((double) rand())/ (double)32767;
    y[i] = ((double) rand())/ (double)32767;
  }
}
int main(){
  int i=0, N = 512*512;
  double *x, *y;
  int N_trials = 200*100; //200*163;
  x = malloc_double(N);
  y = malloc_double(N);

  double d_vcl, d_mkl, d_naive;
  int nProcessors = omp_get_max_threads();
  printf("number of processors: %d \n", nProcessors);
  // omp_set_num_threads(nProcessors);/
  omp_set_num_threads(8);
  int mkln = mkl_get_max_threads();
  printf("mkl-nthreads: %d\n", mkln);

  load_xy(N, x, y);
  d_mkl = cblas_ddot(N, x, 1, y, 1);

  /*---------------------------- DOT-PRODUCT------------------ */
  double begin_mkl =0, t_mkl=0;
  double begin_naive=0, t_naive=0;
  double begin_vcl=0, t_vcl=0;
  printf("------------- dot-product -----------------------------------------\n");
  printf("Method | trial  | Total Time     |  Mean Time     |     value     |\n");

  for(int k=0; k<10; k++){
    double Ntr = (double) N_trials;
    if (1){
      t_naive=0;
      for (i=1; i<=N_trials; i++){
        begin_naive = omp_get_wtime();
        d_naive = dot(N-i, x, y);
        t_naive += omp_get_wtime() - begin_naive;
      }
      double time_naive = (double)(t_naive );
      printf(" naive     %d      %f       %.10e    %.6e\n", k, time_naive, time_naive/Ntr, d_naive);
    }

    if(1){
      if (1){
        t_vcl = 0;
        for (i=1; i<=N_trials; i++){
          begin_vcl = omp_get_wtime();
          d_vcl = vcl_ddot(N-i, x, y);
          t_vcl += omp_get_wtime() - begin_vcl;
        }
        double time_vcl = (double)(t_vcl);
        printf("  VCL      %d      %f       %.10e    %.6e\n", k, time_vcl, time_vcl/Ntr, d_vcl);

      t_mkl = 0;
      for (i=1; i<=N_trials; i++){
        begin_mkl = omp_get_wtime();
        d_mkl = cblas_ddot(N-i, x, 1, y, 1);
        t_mkl += omp_get_wtime() - begin_mkl;
      }
      double time_mkl = (double)(t_mkl);// / 1e6;
      printf("  MKL      %d      %f       %.10e    %.6e\n", k, time_mkl, time_mkl/Ntr, d_mkl);
    }

    }
    mkl_free_buffers();;
  }

  free_double(x);
  free_double(y);
}
