#pragma once
#ifndef __NEW_STAT_H__
#define __NEW_STAT_H__
namespace ONETOOL {
	double	x_alngam(double *n);
	void	x_bratio(double, double, double, double, double*, double*, int*);
	double	x_algdiv(double, double);
	double	x_alnrel(double);
	double	x_apser(double, double, double, double);
	double	x_basym(double*, double*, double*, double*);
	double	x_bcorr(double, double);
	double	x_betaln(double, double);
	double	x_bfrac(double*, double*, double*, double*, double*, double*);
	void	x_bgrat(double*, double*, double*, double*, double*, double*, int*i);
	double	x_bpser(double, double, double, double);
	double	x_brcmp1(int*, double*, double*, double*, double*);
	double	x_brcomp(double*, double*, double*, double*);
	double	x_bup(double*, double*, double*, double*, int*, double*);

	double	x_dlanor(double*);
	double	x_dln1mx(double);
	double	x_dln1px(double);
	double	x_dlnbet(double, double);
	double	x_dlngam(double);
	double	x_dstrem(double);
	double	x_dt1(double, double, double);
	double	x_erf1(double*);
	double	x_erfc1(int*, double*);
	double	x_esum(int, double);
	double	x_exparg(int*);
	double	x_fpser(double, double, double, double);
	double	x_gam1(double);
	void	x_gaminv(double*, double*, double*, double*, double*, int*);
	double	x_gamln(double);
	double	x_gamln1(double);
	double	x_Xgamm(double*);
	void	x_grat1(double*, double*, double*, double*, double*, double*);
	void	x_gratio(double*, double*, double*, double*, int*);
	double	x_gsumln(double*, double*);
	double	x_psi(double);
	double	x_rcomp(double*, double*);
	double	x_rexp(double*);
	double	x_rlog(double*);
	double	x_rlog1(double*);
	double	x_spmpar(const int);
	double	x_stvaln(double*);
	double	x_fifdint(double);
	double	x_fifdmax1(double, double);
	double	x_fifdmin1(double, double);
	double	x_fifdsign(double, double);
	long	x_fifidint(double);
	long	x_fifmod(long, long);
	double	x_stirlerr(double n);
	double	x_truncate(double x);

	double	x_dblSterlingRmder4CompBetaFun(double, double);
	double	x_dblEvalExpXminus1(double);
	double	x_dblEvalPolyX(double[], int, double);
	double	x_dblNormalDistInv(double p, double q);

	void	x_cdfnor(int, double*, double*, double*, double*, double*, int*, double*);
	void	x_cumbet(double, double, double, double, double*, double*);
	void	x_cumbin(double, double, double, double, double*, double*);
	void	x_cumchi(double, double, double*, double*);
	void	x_cumf(double, double, double, double*, double*);
	void	x_cumgam(double*, double*, double*, double*);
	void	x_cumnbn(double, double, double, double, double*, double*);
	void	x_cumnor(double*, double*, double*);
	void	x_cumt(double*, double*, double*, double*);

	double	x_betaPraw(double x, double pin, double qin, int lower_tail, int log_p);
	double	lbeta(double a, double b);
	double	gammafn(double x);
	double	lgammafn_sign(double x, int *sgn);
	double	lgammafn(double x);
	double	lgammacor(double x);
	double	qbeta(double alpha, double p, double q, int lower_tail, int log_p);
} // End namespace ONETOOL

#endif
