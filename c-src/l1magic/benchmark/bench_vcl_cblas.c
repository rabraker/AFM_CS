#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "mkl.h"
#include "vcl_math.h"
#include "l1qc_common.h"


double dot(const int N, double * restrict x, double * restrict y){
  double *x_ = __builtin_assume_aligned(x, 64);
  double *y_ = __builtin_assume_aligned(y, 64);

  double total = 0.0;
  for (int i =0; i<N; i++){
    total += x_[i] * y_[i];
  }

  return total;
}

void axpy(const int N, const double alpha, double * restrict x, double * restrict y){
  double *x_ = __builtin_assume_aligned(x, 64);
  double *y_ = __builtin_assume_aligned(y, 64);

  for (int i =0; i<N; i++){
    y[i] = x_[i] * alpha + y_[i];
  }

}


double norm2(const int N, double * restrict x){
  double *x_ = __builtin_assume_aligned(x, 64);

  double total = 0.0;
  for (int i =0; i<N; i++){
    total += x_[i] * x[i];
  }

  return total;
}


int main(){
  int i=0, N = 50;
  double *x, *y, *z;
  int N_trials = 100;
  double alpha = 3.0;
  x = malloc_double(N);
  y = malloc_double(N);
  z = malloc_double(N);

  double d_vcl, d_mkl, d_naive;

  for (i=0; i<N; i++){
    x[i] = (double) (rand()/32767);
    y[i] = (double) (rand()/32767);
  }

  /*---------------------------- DOT-PRODUCT------------------ */
  /* MKL dot*/
  clock_t begin_mkl = clock();
  for (i=1; i<=N_trials; i++){
    d_mkl = cblas_ddot(N, x, 1, y, 1);
  }
  clock_t end_mkl = clock();

  /* niave dot*/
  clock_t begin_naive = clock();
  for (i=1; i<=N_trials; i++){
    d_naive = dot(N, x, y);
      }
  clock_t end_naive = clock();

  /* VCL dot*/
  clock_t begin_vcl = clock();
  for (i=1; i<=N_trials; i++){
    d_vcl = vcl_ddot(N, x, y);
  }
  clock_t end_vcl = clock();


  double time_vcl = (double)(end_vcl - begin_vcl) / CLOCKS_PER_SEC;
  double time_mkl = (double)(end_mkl - begin_mkl) / CLOCKS_PER_SEC;
  double time_naive = (double)(end_naive - begin_naive) / CLOCKS_PER_SEC;
  double Ntr = (double) N_trials;

  printf("------------- dot-product-------------------------\n");
  printf("VCL-dot: %f,   MKL-dot: %f,  naive-dot:%f\n", d_vcl, d_mkl, d_naive);
  printf("Method | Total Time     |  Mean Time     |\n");
  printf("  VCL     %f       %.10e\n", time_vcl, time_vcl/Ntr);
  printf("  MKL     %f       %.10e\n", time_mkl, time_mkl/Ntr);
  printf(" naive    %f       %.10e\n", time_naive, time_naive/Ntr);


  /*-------------------------2NORM------------------ */
  /* MKL dot*/
  begin_mkl = clock();
  for (i=1; i<=N_trials; i++){
    d_mkl = cblas_ddot(N, x, 1, x, 1);
  }
  end_mkl = clock();

  /* niave dot*/
  begin_naive = clock();
  for (i=1; i<=N_trials; i++){
    d_naive = norm2(N, x);
      }
  end_naive = clock();

  /* VCL dot*/
  begin_vcl = clock();
  for (i=1; i<=N_trials; i++){
    d_vcl = vcl_dnorm2(N, x);
  }
  end_vcl = clock();


  time_vcl = (double)(end_vcl - begin_vcl) / CLOCKS_PER_SEC;
  time_mkl = (double)(end_mkl - begin_mkl) / CLOCKS_PER_SEC;
  time_naive = (double)(end_naive - begin_naive) / CLOCKS_PER_SEC;
  Ntr = (double) N_trials;
  printf("------------- 2-norm-------------------------\n");
  printf("VCL-norm: %f,   MKL-norm: %f,  naive-norm:%f\n", d_vcl, d_mkl, d_naive);
  printf("Method | Total Time     |  Mean Time     |\n");
  printf("  VCL     %f       %.10e\n", time_vcl, time_vcl/Ntr);
  printf("  MKL     %f       %.10e\n", time_mkl, time_mkl/Ntr);
  printf(" naive    %f       %.10e\n", time_naive, time_naive/Ntr);

  /*---------------------------- a*X + Y = Y------------------ */
  /* MKL dot*/
  begin_mkl = clock();
  for (i=1; i<=N_trials; i++){
    cblas_daxpy(N, alpha, x, 1, y, 1);
  }
  end_mkl = clock();

  /* niave dot*/
  begin_naive = clock();
  for (i=1; i<=N_trials; i++){
    axpy(N, alpha, x, y);
      }
  end_naive = clock();

  /* VCL dot*/
  begin_vcl = clock();
  for (i=1; i<=N_trials; i++){
    vcl_daxpy(N, alpha, x, y);
  }
  end_vcl = clock();


  time_vcl = (double)(end_vcl - begin_vcl) / CLOCKS_PER_SEC;
  time_mkl = (double)(end_mkl - begin_mkl) / CLOCKS_PER_SEC;
  time_naive = (double)(end_naive - begin_naive) / CLOCKS_PER_SEC;
  Ntr = (double) N_trials;
  printf("------------- axpy-------------------------\n");

  printf("Method | Total Time     |  Mean Time     |\n");
  printf("  VCL     %f       %.10e\n", time_vcl, time_vcl/Ntr);
  printf("  MKL     %f       %.10e\n", time_mkl, time_mkl/Ntr);
  printf(" naive    %f       %.10e\n", time_naive, time_naive/Ntr);

}
