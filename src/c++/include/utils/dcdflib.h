#pragma once
#ifndef __WISARD_DCDFLIB_H__
#define __WISARD_DCDFLIB_H__

#include "global/common.h"

namespace ONETOOL {

	extern void cdfbet(int*, double*, double*, double*, double*, double*, double*,
		int*, double*);
	extern void cdfbin(int*, double*, double*, double*, double*, double*, double*,
		int*, double*);
	extern void cdfchi(int*, double*, double*, double*, double*, int*, double*);
	extern void cdfchn(int*, double*, double*, double*, double*, double*, int*, double*);
	extern void cdff(int*, double*, double*, double*, double*, double*, int*, double*);
	extern void cdffnc(int*, double*, double*, double*, double*, double*, double*,
		int*s, double*);
	extern void cdfgam(int*, double*, double*, double*, double*, double*, int*, double*);
	extern void cdfnbn(int*, double*, double*, double*, double*, double*, double*,
		int*, double*);
	extern void cdfpoi(int*, double*, double*, double*, double*, int*, double*);
	extern void cdft(int*, double*, double*, double*, double*, int*, double*);

	double dchisq(double x, double df, int give_log);
	double dnorm(double x, double mu = 0.0, double sigma = 1.0, int give_log = 0);
	wsReal ptdist(wsReal *Rp_tStat, wsReal *Rp_tDF);

	double dbinom(double x, double n, double p, int give_log = 0);
	double stirlerr(double n);
	double pgamma(double x, double alph, double scale, int lower_tail = 1, int log_p = 0);
	double qgamma(double p, double alpha, double scale, int lower_tail = 1, int log_p = 0);
	double pf(double x, double df1, double df2, int lower_tail);
	double pT(double t, double df);

	double dbeta(double x, double a, double b, int log_p);
	double phyper(wsUint _x, wsUint NR, wsUint NB, wsUint n,
		int lower_tail = 1, int log_p = 0);

} // End namespace ONETOOL

#endif
