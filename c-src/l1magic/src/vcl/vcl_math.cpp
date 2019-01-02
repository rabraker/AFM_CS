#define VCL_NAMESPACE vcl
#include "vectorclass.h"
using namespace vcl;

#include "vectormath_exp.h"
#include <math.h>
#include <stdio.h>
#define VECTORSIZE 4

extern "C" double vcl_logsum(const int N, const double alpha, const double *x){

  const int regularpart = N & (~(VECTORSIZE-1));
  // (AND-ing with -vectorsize-1 will round down to nearest
  // lower multiple of vectorsize. This works only if
  // vectorsize is a power of 2)

  int i=0;
  double total = 0.0;
  Vec4d xvec, sum_vec(0);
  // if (N>4){
    for(i=0; i<regularpart; i+=VECTORSIZE){
      xvec.load_a(x+i);
      sum_vec += vcl::log(xvec * alpha);
    }
  // }
  for(; i<N; i++){
    total += log(x[i]*alpha);
  }
  total += vcl::horizontal_add(sum_vec);
  return total;
}


extern "C" double vcl_sum(const int N, const double *x){
  const int regularpart = N & (~(VECTORSIZE-1));
  // (AND-ing with -vectorsize-1 will round down to nearest
  // lower multiple of vectorsize. This works only if
  // vectorsize is a power of 2)

  int i=0;
  double total = 0.0;
  Vec4d xvec, sum_vec(0.0);
  for(i=0; i<regularpart; i+=VECTORSIZE){
    xvec.load_a(x+i);
    sum_vec += xvec;
  }

  for(; i<N; i++){
    total += x[i];
  }
  total += vcl::horizontal_add(sum_vec);
  return total;
}

extern "C" double vcl_ddot(const int N, const double *x, const double *y){
  const int regularpart = N & (~(VECTORSIZE-1));
  // (AND-ing with -vectorsize-1 will round down to nearest
  // lower multiple of vectorsize. This works only if
  // vectorsize is a power of 2)

  int i=0;
  double total = 0.0;
  Vec4d xvec, yvec, sum_vec(0.0);
  for(i=0; i<regularpart; i+=VECTORSIZE){
    xvec.load_a(x+i);
    yvec.load_a(y+i);
    sum_vec += xvec * yvec;
  }

  for(; i<N; i++){
    total += x[i] * y[i];
  }
  total += vcl::horizontal_add(sum_vec);
  return total;
}

extern "C" double vcl_dnorm2(const int N, const double *x){
  const int regularpart = N & (~(VECTORSIZE-1));
  // (AND-ing with -vectorsize-1 will round down to nearest
  // lower multiple of vectorsize. This works only if
  // vectorsize is a power of 2)

  int i=0;
  double total = 0.0;
  Vec4d xvec, sum_vec(0.0);
  for(i=0; i<regularpart; i+=VECTORSIZE){
    xvec.load_a(x+i);
    sum_vec += vcl::square(xvec) ;
  }

  for(; i<N; i++){
    total += x[i] * x[i];
  }
  total += vcl::horizontal_add(sum_vec);
  return total;
}


extern "C" void vcl_daxpy(const int N, const double alpha, const double *x, double *y){
  const int regularpart = N & (~(VECTORSIZE-1));
  // (AND-ing with -vectorsize-1 will round down to nearest
  // lower multiple of vectorsize. This works only if
  // vectorsize is a power of 2)

  int i=0;
  Vec4d xvec, yvec, zvec;
  for(i=0; i<regularpart; i+=VECTORSIZE){
    xvec.load_a(x+i);
    yvec.load_a(y+i);

    zvec = xvec * alpha + yvec;
    zvec.store(y+i);
  }

  for(; i<N; i++){
    y[i] = x[i] * alpha + y[i];
  }
}

extern "C" void vcl_daxpby(const int N, const double alpha, const double *x, const double beta, double *y){
  const int regularpart = N & (~(VECTORSIZE-1));
  // (AND-ing with -vectorsize-1 will round down to nearest
  // lower multiple of vectorsize. This works only if
  // vectorsize is a power of 2)

  int i=0;
  Vec4d xvec, yvec, sum_vec(0.0);
  for(i=0; i<regularpart; i+=VECTORSIZE){
    xvec.load_a(x+i);
    yvec.load_a(y+i);

    yvec = xvec * alpha + beta * yvec;
    yvec.store(y+i);
  }

  for(; i<N; i++){
    y[i] = x[i] * alpha + beta * y[i];
  }
}

extern "C" void vcl_dxpby(const int N, const double *x, const double beta, double *y){
  const int regularpart = N & (~(VECTORSIZE-1));
  // (AND-ing with -vectorsize-1 will round down to nearest
  // lower multiple of vectorsize. This works only if
  // vectorsize is a power of 2)

  int i=0;
  Vec4d xvec, yvec, sum_vec(0.0);
  for(i=0; i<regularpart; i+=VECTORSIZE){
    xvec.load_a(x+i);
    yvec.load_a(y+i);

    yvec = xvec + beta * yvec;
    yvec.store(y+i);
  }

  for(; i<N; i++){
    y[i] = x[i] + beta * y[i];
  }
}
