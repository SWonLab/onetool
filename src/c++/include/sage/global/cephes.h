#ifndef CEPHES_H
#define CEPHES_H

#include "sage/global/definition.h"

namespace SAGE
{

//  Prototypes for probability functions from the cephes library

// Combinatorial functions
double factorial(int n);
double log_factorial(int n);
double choose(int n, int k);
double log_choose(int n, int k);

double bin_prob(size_t q, double p, size_t x);

// Binomial distribution
//double bdtr(int k, int n, double p);
//double bdtrc(int k, int n, double p); // complemented

// Negative binomial distribution
//double nbdtr(int k, int n, double p);
//double nbdtrc(int k, int n, double p); // complemented

// Beta distribution
//double incbet(double a, double b, double x); // incomplete beta intgral
//double incbi(double a, double b, double p);  // inv incomplete beta integral

// Chi-square distribution
double chdtr(double df, double x);
double chdtrc(double df, double x);
double chdtri(double df, double p);

// F distribution
//double fdtr(int df1, int df2, double x);
//double fdtrc(int df1, int df2, double x); // complemented
//double fdtri(int df1, int df2, double p); // inverse of complemented

// Gamma function
double gammafn(double x);         // gamma function
double lgam(double x);            // log of gamma function
double igam(double a, double x);  // incomplete gamma integral
double igamc(double a, double x); // complemented
double igami(double a, double p); // inverse incomplete gamma integral

// Gamma distribution
//double gdtr(double a, double b, double x);
//double gdtrc(double a, double b, double x); // complemented

// Normal distribution
double ndtr(double x);
double ndtri(double p); // inverse

// Poisson distribution
//double pdtr(int k, double m);
//double pdtrc(int k, double m); // complemented
//double pdtri(int k, double p); // inverse

// Student's t distribution
//double stdtr(int k, double t);
//double stdtri(int k, double p); // inverse
//double stdtr_mean_test(int k, double t);    // density from -t to t
//double stdtr_mean_test_inv(int k, double p); // t for density -t to t for p

#ifdef NEEDS_LOG1P
extern "C" double log1p(double x);
#endif

#ifdef NEEDS_EXPM1
extern "C" double expm1(double x);
#endif

#ifdef NEEDS_COSM1
extern "C" double cosm1(double x);
#endif

// inlined poly eval functions

inline double polevl(double x, double coef[], int N)
{
  double ans;
  int i;
  double *p;

  p = coef;
  ans = *p++;
  i = N;

  do
    ans = ans * x  +  *p++;
  while( --i );
  
  return( ans );
}

inline double p1evl(double x, double coef[], int N)
{
  double ans;
  double *p;
  int i;

  p = coef;
  ans = x + *p++;
  i = N-1;

  do
    ans = ans * x  + *p++;
  while( --i );

  return( ans );
}

}

#endif
