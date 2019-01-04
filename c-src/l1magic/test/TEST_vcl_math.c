/*
This is a test suite for the fourier integration functions contained in
../src/ss_fourier_functions.c. The data used to test these functions is
contained in the header file test_data_ss_ff.h, which defines several global
variables. The header file is generated from the matlab script called
generate_test_data.m

This test suite uses the libcheck framework. On my computer, this got installed
into /usr/local/lib, which was not by default found by the system. Thus, I have
to do
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/local/lib


  https://libcheck.github.io/check/doc/check_html/check_4.html#No-Fork-Mode

 */


#define CK_FLOATING_DIG 20
#include <stdlib.h>
#include <stdio.h>
#include <math.h> //Constants
#include <check.h>

#include "l1qc_common.h"

#include "cJSON.h"
#include "json_utils.h"
#include "vcl_math.h"
#include "check_utils.h"

/* Tolerances and things */
#include "test_constants.h"
#include "check_utils.h"
#include "l1qc_common.h"



START_TEST(test_vcl_sum)
{
  double x_[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
  double *x = malloc_double(12);
  for (int i =0; i<12; i++){
    x[i] = x_[i];
  }

  double sum_exp0 = 6.0;
  double sum_exp1 = 36.0;
  double sum_exp2 = 45.0;
  double sum_exp3 = 55.0;
  double sum_exp4 = 66.0;
  double sum_exp5 = 78.0;
  double sum_x = 0.0;

  sum_x = vcl_sum(3, x);
  ck_assert_double_eq_tol(sum_exp0, sum_x, TOL_DOUBLE);

  sum_x = vcl_sum(8, x);
  ck_assert_double_eq_tol(sum_exp1, sum_x, TOL_DOUBLE);

  sum_x = vcl_sum(9, x);
  ck_assert_double_eq_tol(sum_exp2, sum_x, TOL_DOUBLE);

  sum_x = vcl_sum(10, x);
  ck_assert_double_eq_tol(sum_exp3, sum_x, TOL_DOUBLE);

  sum_x = vcl_sum(11, x);
  ck_assert_double_eq_tol(sum_exp4, sum_x, TOL_DOUBLE);

  sum_x = vcl_sum(12, x);
  ck_assert_double_eq_tol(sum_exp5, sum_x, TOL_DOUBLE);

  free_double(x);
}
END_TEST


START_TEST(test_vcl_logsum)
{
  double x_[] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0};
  double *x = malloc_double(12);
  for (int i =0; i<12; i++){
    x[i] = x_[i];
  }
  double logsum_exp0 = 5.087596335232384;
  double logsum_exp1 = 19.393501212090126;
  double logsum_exp2 = 22.689338078094455;
  double logsum_exp3 = 26.090535459756609;
  double logsum_exp4 = 29.587043021223089;
  double logsum_exp5 = 33.170561959679198;

  double logsum_x = 0.0;

  double alpha = 3.0;

  logsum_x = vcl_logsum(3, alpha, x);
  ck_assert_double_eq_tol(logsum_exp0, logsum_x, TOL_DOUBLE);

  logsum_x = vcl_logsum(8, alpha, x);
  ck_assert_double_eq_tol(logsum_exp1, logsum_x, TOL_DOUBLE);

  logsum_x = vcl_logsum(9, alpha, x);
  ck_assert_double_eq_tol(logsum_exp2, logsum_x, TOL_DOUBLE);

  logsum_x = vcl_logsum(10, alpha, x);
  ck_assert_double_eq_tol(logsum_exp3, logsum_x, TOL_DOUBLE);

  logsum_x = vcl_logsum(11, alpha, x);
  ck_assert_double_eq_tol(logsum_exp4, logsum_x, TOL_DOUBLE);

  logsum_x = vcl_logsum(12, alpha, x);
  ck_assert_double_eq_tol(logsum_exp5, logsum_x, TOL_DOUBLE);

  free_double(x);
}
END_TEST

START_TEST(test_vcl_daxpy_simple)
{
  double alpha = 2.0;
  double *x = malloc_double(12);
  double *y = malloc_double(12);
  double *y_tmp = malloc_double(12);
  double *z_exp = malloc_double(12);
  for (int i =0; i<12; i++){
    x[i] = 1.0;
    y[i] = 2.0;
    z_exp[i] = 4.0; //1*2 + *2 = 4;
  }


  cblas_dcopy(12, y, 1, y_tmp, 1);
  vcl_daxpy(3, alpha, x, y_tmp);
  ck_assert_double_array_eq_tol(3, z_exp, y_tmp, TOL_DOUBLE);

  cblas_dcopy(12, y, 1, y_tmp, 1);
  vcl_daxpy(8, alpha, x, y_tmp);
  ck_assert_double_array_eq_tol(8, z_exp, y_tmp, TOL_DOUBLE);

  cblas_dcopy(12, y, 1, y_tmp, 1);
  vcl_daxpy(9, alpha, x, y_tmp);
  ck_assert_double_array_eq_tol(9, z_exp, y_tmp, TOL_DOUBLE);

  cblas_dcopy(12, y, 1, y_tmp, 1);
  vcl_daxpy(10, alpha, x, y_tmp);
  ck_assert_double_array_eq_tol(10, z_exp, y_tmp, TOL_DOUBLE);

  cblas_dcopy(12, y, 1, y_tmp, 1);
  vcl_daxpy(11, alpha, x, y_tmp);
  ck_assert_double_array_eq_tol(11, z_exp, y_tmp, TOL_DOUBLE);

  cblas_dcopy(12, y, 1, y_tmp, 1);
  vcl_daxpy(12, alpha, x, y_tmp);
  ck_assert_double_array_eq_tol(12, z_exp, y_tmp, TOL_DOUBLE);


  free_double(x);
  free_double(y);
  free_double(z_exp);
  free_double(y_tmp);
}
END_TEST

START_TEST(test_vcl_daxpy_large)
{
  int N = 512*512;
  double alpha = 2.0;
  double *x = malloc_double(N);
  double *y = malloc_double(N);
  double *y_tmp = malloc_double(N);
  double *z_exp = malloc_double(N);
  for (int i =0; i<N; i++){
    x[i] = ( (double)rand()) / (double)65534;
    y[i] = ( (double)rand()) / (double)65534;
  }


  cblas_dcopy(N, y, 1, z_exp, 1);
  cblas_daxpy(N, alpha, x, 1, z_exp, 1);

  cblas_dcopy(N, y, 1, y_tmp, 1);
  vcl_daxpy(N, alpha, x, y_tmp);
  ck_assert_double_array_eq_tol(N, z_exp, y_tmp, TOL_DOUBLE_SUPER);

  free_double(x);
  free_double(y);
  free_double(z_exp);
  free_double(y_tmp);
}
END_TEST


START_TEST(test_vcl_dxpby_large)
{
  int N = 512*512;
  double beta = 2.0;
  double *x = malloc_double(N);
  double *y = malloc_double(N);
  double *y_tmp = malloc_double(N);
  double *z_exp = malloc_double(N);
  for (int i =0; i<N; i++){
    x[i] = ( (double)rand()) / (double)65534;
    y[i] = ( (double)rand()) / (double)65534;
  }


  cblas_dcopy(N, y, 1, z_exp, 1);
  cblas_daxpby(N, 1.0, x, 1, beta, z_exp, 1);

  cblas_dcopy(N, y, 1, y_tmp, 1);
  vcl_dxpby(N, x, beta, y_tmp);
  ck_assert_double_array_eq_tol(N, z_exp, y_tmp, TOL_DOUBLE_SUPER);

  free_double(x);
  free_double(y);
  free_double(z_exp);
  free_double(y_tmp);
}
END_TEST


START_TEST(test_vcl_ddot)
{
  int N = 512*512;
  double *x = malloc_double(N);
  double *y = malloc_double(N);
  double dd_exp = 0, dd_act=0;
  for (int i =0; i<N; i++){
    x[i] = ( (double)rand()) / (double)65534;
    y[i] = ( (double)rand()) / (double)65534;
  }


  dd_exp = cblas_ddot(N, x, 1, y, 1);
  dd_act = vcl_ddot(N, x, y);

  ck_assert_double_eq_tol(dd_exp, dd_act, TOL_DOUBLE_SUPER);

  free_double(x);
  free_double(y);
}
END_TEST


START_TEST(test_vcl_dnorm2)
{
  int N = 512*512;
  double *x = malloc_double(N);
  double dd_exp = 0, dd_act=0;
  for (int i =0; i<N; i++){
    x[i] = ( (double)rand() ) / (double)65534;
  }


  dd_exp = cblas_ddot(N, x, 1, x, 1);
  dd_act = vcl_dnorm2(N, x);

  ck_assert_double_eq_tol(dd_exp, dd_act, TOL_DOUBLE_SUPER);

  free_double(x);
}
END_TEST


Suite *vcl_math_suite(void)
{
  Suite *s;

  TCase  *tc_mathfuns;
  s = suite_create("vcl_math");


  tc_mathfuns = tcase_create("vcl_math");
  tcase_add_test(tc_mathfuns, test_vcl_sum);
  tcase_add_test(tc_mathfuns, test_vcl_logsum);
  tcase_add_test(tc_mathfuns, test_vcl_daxpy_simple);
  tcase_add_test(tc_mathfuns, test_vcl_daxpy_large);
  tcase_add_test(tc_mathfuns, test_vcl_dxpby_large);
  tcase_add_test(tc_mathfuns, test_vcl_ddot);
  tcase_add_test(tc_mathfuns, test_vcl_dnorm2);

  /*Add test cases to the suite */
  suite_add_tcase(s, tc_mathfuns);

  return s;

}
