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
  for(i=0; i<N-4+1; i+=4){
    xvec.load_a(x+i);
    sum_vec += xvec;
  }

  for(; i<N; i++){
    total += x[i];
  }
  total += vcl::horizontal_add(sum_vec);
  return total;
}
