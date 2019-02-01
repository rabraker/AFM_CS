#ifndef _FIT_ZPK_
#define _FIT_ZPK_


#if defined _WIN32 || defined __CYGWIN__ || defined __MINGW32__
#ifdef BUILDING_DLL
#ifdef __GNUC__
#define DLL_PUBLIC __attribute__ ((dllexport))
#else
#define DLL_PUBLIC __declspec(dllexport) // Note: actually gcc seems to also supports this syntax.
#endif
#else
#ifdef __GNUC__
#define DLL_PUBLIC __attribute__ ((dllimport))
#else
#define DLL_PUBLIC __declspec(dllimport) // Note: actually gcc seems to also supports this syntax.
#endif
#endif
#define DLL_LOCAL
#else
#if __GNUC__ >= 4
#define DLL_PUBLIC __attribute__ ((visibility ("default")))
#define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#else
#define DLL_PUBLIC
#define DLL_LOCAL
#endif
#endif

#include <complex.h>

int absfunc(int m, int n, double *theta, double *resid, double **dvec, void *user_data);

DLL_PUBLIC  int fit_sos(int Len_omegas, double *omegas, double *resp_real,
                        double *resp_imag, int len_theta, double *theta,
                        double *result);

#endif
