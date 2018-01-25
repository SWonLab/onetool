#ifndef __WISARD_STAT_H__
#define __WISARD_STAT_H__
#pragma once

namespace ONETOOL {

/* Fisher's Exact Test (2-side) */
wsRealCst		fisher(int N_11, int N_12, int N_21, int N_22);
wsRealCst		fisher(int *Na_2x2);
double		fisher(wsReal R_11, wsReal R_12, wsReal R_13, wsReal R_21, wsReal R_22, wsReal R_23);

/* ����-�� ���� */
double	ranksum(int _cx, double *x, int _cy, double *y);

struct RANK {		/* ranksum() �Լ��� ���� ����ü */
	double R_score;
	double rank;
	int B_isY;
};

double	betaP(double x, double pin, double qin, int lower_tail, int log_p);
double	binoP(double x, double n, double p, int lower_tail=1, int log_p=0);

double	davies(double R_Q, double *Ra_lambda, wsUint N_lambda);
double	farebrother(double R_Q, double *Ra_lambda, wsUint N_lambda);
double	ruben(double *lambda, int *mult, double *delta, int n, double c,
	double *mode, int *maxit, double *eps, double *dnsty, int *ifault);

double	metaFisher(wsUint N_pv, wsVec Ra_pv);
double	metaStouffer(wsUint N_pv, wsVec Ra_pv, wsVec Ra_weight=NULL);
double	metaKost(wsUint N_pv, wsVec Ra_pv, wsSym Ra_cor);

/* bonferroni correction */
wsRealCst	bonferroni(wsReal p, int n);
double		qnorm(double p, double mu=0.0, double sigma=1.0, int lower_tail=1, int log_p=0);
double		PVchisq(double x, double df, double *nc=NULL, char B_quiet=0);
double		norprobP(double x, double mean=W0, double sd=W1);
double		pnorm(double x, double mu=0.0, double sigma=1.0, int lower_tail=1, int log_p=0);

wsReal	test2x2(wsUint N_11, wsUint N_12, wsUint N_21, wsUint N_22,
	char B_corr);
wsReal	test2x2_PEDCMC(wsUint N_11, wsUint N_12, wsUint N_21, wsUint N_22,
	char B_corr, wsUint *Np_typeTest);

} // End namespace ONETOOL

#endif
