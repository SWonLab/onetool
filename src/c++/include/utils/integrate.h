#pragma once
#ifndef __WISARD_INTEGRATE_H__
#define __WISARD_INTEGRATE_H__

#include <limits>
#include "global/common.h"
#include "utils/util.h"

using namespace std;

namespace ONETOOL {

struct _xIntgInp;
typedef void (* funcIntg)(double *, int, struct _xIntgInp *);
typedef struct _xIntgInp
{
	funcIntg H_func;
	double R_lower, R_upper;
	int N_subDiv;
	double R_absTol, R_relTol;
	double *Ra_X;
	void *Vp_param;
	_xIntgInp(funcIntg H_inpFunc, double R_inpLower, double R_inpUpper,
		void *Vp_inpParam, int N_sDiv=100, double R_inpAtol=1e-6) {
			H_func = H_inpFunc;
			R_lower = R_inpLower;
			R_upper = R_inpUpper;
			Vp_param = Vp_inpParam;
			N_subDiv = N_sDiv;
			R_absTol = R_inpAtol;//pow(numeric_limits<double>::epsilon(), 0.25);
			R_relTol = R_absTol;
			Ra_X = NULL;
	}
	~_xIntgInp() {
		DEALLOC(Ra_X);
	}
} xIntgInp;

typedef struct _xIntgRes
{
	double R_value;
	double R_absErr;
	int N_subDiv;
	wsUint N_err;
} xIntgRes;

xIntgRes* integrate(xIntgInp *Xp_inp, char B_stopOnErr=0);

} // End namespace ONETOOL

#endif
