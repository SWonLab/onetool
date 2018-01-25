#pragma once
#ifndef __CDFLIB_H__
#define __CDFLIB_H__

namespace ONETOOL {

	void cdfbet(int*, double*, double*, double*, double*, double*, double*,
		int*, double*);
	void cdfbin(int*, double*, double*, double*, double*, double*, double*,
		int*, double*);
	void cdfchi(int*, double*, double*, double*, double*, int*, double*);
	void cdfchn(int*, double*, double*, double*, double*, double*, int*, double*);
	void cdff(int*, double*, double*, double*, double*, double*, int*, double*);
	void cdffnc(int*, double*, double*, double*, double*, double*, double*,
		int*s, double*);
	void cdfgam(int*, double*, double*, double*, double*, double*, int*, double*);
	void cdfnbn(int*, double*, double*, double*, double*, double*, double*,
		int*, double*);
	void cdfpoi(int*, double*, double*, double*, double*, int*, double*);
	void cdft(int*, double*, double*, double*, double*, int*, double*);
	void cumchn(double, double, double, double*, double*);
	void cumfnc(double*, double*, double*, double*, double*, double*);
	void x_cumpoi(double, double, double*, double*);
	void dinvr(int*, double*, double, unsigned long*, unsigned long*);
	void dstinv(double, double, double, double, double, double,
		double);
	void dzror(int*, double*, double, double*, double *,
		unsigned long*qleft = NULL, unsigned long*qhi = NULL);
	void dstzr(double zxlo, double zxhi, double zabstl, double zreltl);

} // End namespace ONETOOL

#endif
