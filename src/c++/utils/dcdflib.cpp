#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global/common.h"
#include "utils/cdflib.h"
#include "utils/util.h"
#include "utils/stat.h"
#include "utils/new_stat.h"
//#define lgammafn(x) alngam(&x)
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-util_expm1(x)) : util_log1p(-exp(x)))
#define M_LN_SQRT_2PI 0.918938533204672741780329736406 
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */
#define R_P_bounds_01(x, x_min, x_max) 	\
	if(x <= x_min) return R_DT_0;		\
	if(x >= x_max) return R_DT_1

namespace ONETOOL {

	pthread_mutex_t	H_ctrl;
	double x_stirlerr(double n);

	inline int ipmpar_v2(const int i)
	{
		if (i < 1 || i > 10) halt("Invalid range");
		static const int imach[11] = { 0, 2, 31, 2147483647, 2, 24, -125,
			128, 53, -1021, 1024 };
		return imach[i];
	}

	inline double exparg_v2(int l)
		/*
		--------------------------------------------------------------------
		IF L = 0 THEN  EXPARG(L) = THE LARGEST POSITIVE W FOR WHICH
		EXP(W) CAN BE COMPUTED.

		IF L IS NONZERO THEN  EXPARG(L) = THE LARGEST NEGATIVE W FOR
		WHICH THE COMPUTED VALUE OF EXP(W) IS NONZERO.

		NOTE... ONLY AN APPROXIMATE VALUE FOR EXPARG(L) IS NEEDED.
		--------------------------------------------------------------------
		*/
	{
		double exparg,lnb;
		int b,m;
		//	  0  1  2   3           4  5   6     7    8   9      10
		//	{ 0, 2, 31, 2147483647, 2, 24, -125, 128, 53, -1021, 1024 }
		b = 2;
		if(b == 2)
			lnb = 0.693147180559945309417232121458176568e0;
		else if(b == 8)
			lnb = 2.0794415416798359282516963643745297e0;
		else if(b == 16)
			lnb = 2.7725887222397812376689284858327062723e0;
		else
			lnb = log((double)b);
		if(l == 0) {
			m = 1024;
			exparg = 0.99999e0*((double)m*lnb);
		} else {
			m = (-1021)-1;
			exparg = 0.99999e0*((double)m*lnb);
		}
		return exparg;
	}

	/* EVALUATION OF THE REAL ERROR FUNCTION */
	inline double erf1_v2(double x)
	{
		static const double c = .564189583547756e0;
		static const double a[5] = {
			.771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
			.479137145607681e-01,.128379167095513e+00
		};
		static const double b[3] = {
			.301048631703895e-02,.538971687740286e-01,.375795757275549e+00
		};
		static const double p[8] = {
			-1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
			4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
			4.51918953711873e+02,3.00459261020162e+02
		};
		static const double q[8] = {
			1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
			2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
			7.90950925327898e+02,3.00459260956983e+02
		};
		static const double r[5] = {
			2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
			4.65807828718470e+00,2.82094791773523e-01
		};
		static const double s[4] = {
			9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
			1.80124575948747e+01
		};
		double erf1,ax,bot,t,top;

		ax = fabs(x);

		if(ax <= 0.5e0) {
			t = x*x;
			top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0e0;
			bot = ((b[0]*t+b[1])*t+b[2])*t+1.0e0;
			erf1 = x*(top/bot);
		} else if (ax <= 4.0e0) {
			top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax+p[6])*ax+p[
				7];
			bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax+q[6])*ax+q[
				7];
			erf1 = 0.5e0+(0.5e0-exp(-(x*x))*top/bot);
			if(x < 0.0e0)
				erf1 = -erf1;
		} else if(ax >= 5.8e0)
			erf1 = x_fifdsign(1.0e0,x);
		else {
			double x2 = x*x;
			t = 1.0e0/x2;
			top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
			bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0e0;
			erf1 = (c-top/(x2*bot))/ax;
			erf1 = 0.5e0+(0.5e0-exp(-x2)*erf1);
			if (x < 0.0e0)
				erf1 = -erf1;
		}

		return erf1;
	}

	inline double erfc1_v2(int ind, double x)
		/*
		-----------------------------------------------------------------------
		EVALUATION OF THE COMPLEMENTARY ERROR FUNCTION

		ERFC1(IND,X) = ERFC(X)            IF IND = 0
		ERFC1(IND,X) = EXP(X*X)*ERFC(X)   OTHERWISE
		-----------------------------------------------------------------------
		*/
	{
		static const double c = .564189583547756e0;
		static const double a[5] = {
			.771058495001320e-04,-.133733772997339e-02,.323076579225834e-01,
			.479137145607681e-01,.128379167095513e+00
		};
		static const double b[3] = {
			.301048631703895e-02,.538971687740286e-01,.375795757275549e+00
		};
		static const double p[8] = {
			-1.36864857382717e-07,5.64195517478974e-01,7.21175825088309e+00,
			4.31622272220567e+01,1.52989285046940e+02,3.39320816734344e+02,
			4.51918953711873e+02,3.00459261020162e+02
		};
		static const double q[8] = {
			1.00000000000000e+00,1.27827273196294e+01,7.70001529352295e+01,
			2.77585444743988e+02,6.38980264465631e+02,9.31354094850610e+02,
			7.90950925327898e+02,3.00459260956983e+02
		};
		static const double r[5] = {
			2.10144126479064e+00,2.62370141675169e+01,2.13688200555087e+01,
			4.65807828718470e+00,2.82094791773523e-01
		};
		static const double s[4] = {
			9.41537750555460e+01,1.87114811799590e+02,9.90191814623914e+01,
			1.80124575948747e+01
		};
		static const int K1 = 1;
		double ax, t, top, bot, erfc1;

		/*
		ABS(X) .LE. 0.5
		*/
		ax = fabs(x);
		if(ax <= 0.5e0) {
			t = x*x;
			top = (((a[0]*t+a[1])*t+a[2])*t+a[3])*t+a[4]+1.0e0;
			bot = ((b[0]*t+b[1])*t+b[2])*t+1.0e0;
			erfc1 = 0.5e0+(0.5e0-x*(top/bot));
			if(ind != 0)
				erfc1 = exp(t)*erfc1;
			return erfc1;
		}
		/*
		0.5 .LT. ABS(X) .LE. 4
		*/
		if(ax <= 4.0e0) {
			top = ((((((p[0]*ax+p[1])*ax+p[2])*ax+p[3])*ax+p[4])*ax+p[5])*ax
				+p[6])*ax+p[7];
			bot = ((((((q[0]*ax+q[1])*ax+q[2])*ax+q[3])*ax+q[4])*ax+q[5])*ax
				+q[6])*ax+q[7];
			erfc1 = top/bot;
		} else {
			/*
			ABS(X) .GT. 4
			*/
			if(x <= -5.6e0) {
				/*
				LIMIT VALUE FOR LARGE NEGATIVE X
				*/
				erfc1 = 2.0e0;
				if(ind != 0)
					erfc1 = 2.0e0*exp(x*x);
				return erfc1;
			}
			if(ind == 0) {
				if(x > 100.0e0 || (x*x) > -exparg_v2(K1)) {
					/*
					LIMIT VALUE FOR LARGE POSITIVE X
					WHEN IND = 0
					*/
					erfc1 = 0.0e0;
					return erfc1;
				}
			}
			t = SQR(1.0e0/x);
			top = (((r[0]*t+r[1])*t+r[2])*t+r[3])*t+r[4];
			bot = (((s[0]*t+s[1])*t+s[2])*t+s[3])*t+1.0e0;
			erfc1 = (c-t*top/bot)/ax;
		}
		/*
		FINAL ASSEMBLY
		*/
		if(ind == 0) {
			double w = x*x;
			t = w;
			double e = w-t;
			erfc1 = (0.5e0+(0.5e0-e))*exp(-t)*erfc1;
			if(x < 0.0e0) erfc1 = 2.0e0-erfc1;
			return erfc1;
		}
		if(x < 0.0e0)
			erfc1 = 2.0e0*exp(x*x)-erfc1;
		return erfc1;
	}

	/* EVALUATION OF THE FUNCTION EXP(X) - 1 */
	inline double rexp_v2(double x)
	{
		static const double p1 = .914041914819518e-09;
		static const double p2 = .238082361044469e-01;
		static const double q1 = -.499999999085958e+00;
		static const double q2 = .107141568980644e+00;
		static const double q3 = -.119041179760821e-01;
		static const double q4 = .595130811860248e-03;
		double rexp, w;
		/*
		..
		.. Executable Statements ..
		*/
		if(fabs(x) > 0.15e0) {
			w = exp(x);
			if(x > 0.0e0)
				rexp = w*(0.5e0+(0.5e0-1.0e0/w));
			else
				rexp = w-0.5e0-0.5e0;
		} else
			rexp = x*(((p2*x+p1)*x+1.0e0)/((((q4*x+q3)*x+q2)*x+q1)*x+1.0e0));

		return rexp;
	}

	/* COMPUTATION OF  X - 1 - LN(X) */
	inline double x_rlog2(double x)
	{
		static const double a = .566749439387324e-01;
		static const double b = .456512608815524e-01;
		static const double p0 = .333333333333333e+00;
		static const double p1 = -.224696413112536e+00;
		static const double p2 = .620886815375787e-02;
		static const double q1 = -.127408923933623e+01;
		static const double q2 = .354508718369557e+00;
		double rlog, w, w1;
		double u, r, t;

		if(x < 0.61e0 || x > 1.57e0) {
			r = x-0.5e0-0.5e0;
			rlog = r-log(x);
		} else {
			if(x < 0.82e0) {
				u = x-0.7e0;
				u /= 0.7e0;
				w1 = a-u*0.3e0;
			} else if(x > 1.18e0) {
				u = 0.75e0*x-1.e0;
				w1 = b+u/3.0e0;
			} else {
				/*
				ARGUMENT REDUCTION
				*/
				u = x-0.5e0-0.5e0;
				w1 = 0.0e0;
			}
			/*
			SERIES EXPANSION
			*/
			r = u/(u+2.0e0);
			t = r*r;
			w = ((p2*t+p1)*t+p0)/((q2*t+q1)*t+1.0e0);
			rlog = 2.0e0*t*(1.0e0/(1.0e0-r)-r*w)+w1;
		}

		return rlog;
	}


	void cdfbet(int *which,double *p,double *q,double *x,double *y,
		double *a,double *b,int *status,double *bound)
		/**********************************************************************

		void cdfbet(int *which,double *p,double *q,double *x,double *y,
		double *a,double *b,int *status,double *bound)

		Cumulative Distribution Function
		BETa Distribution


		Function


		Calculates any one parameter of the beta distribution given
		values for the others.


		Arguments


		WHICH --> Integer indicating which of the next four argument
		values is to be calculated from the others.
		Legal range: 1..4
		iwhich = 1 : Calculate P and Q from X,Y,A and B
		iwhich = 2 : Calculate X and Y from P,Q,A and B
		iwhich = 3 : Calculate A from P,Q,X,Y and B
		iwhich = 4 : Calculate B from P,Q,X,Y and A

		P <--> The integral from 0 to X of the chi-square
		distribution.
		Input range: [0, 1].

		Q <--> 1-P.
		Input range: [0, 1].
		P + Q = 1.0.

		X <--> Upper limit of integration of beta density.
		Input range: [0,1].
		Search range: [0,1]

		Y <--> 1-X.
		Input range: [0,1].
		Search range: [0,1]
		X + Y = 1.0.

		A <--> The first parameter of the beta density.
		Input range: (0, +infinity).
		Search range: [1D-300,1D300]

		B <--> The second parameter of the beta density.
		Input range: (0, +infinity).
		Search range: [1D-300,1D300]

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound
		3 if P + Q .ne. 1
		4 if X + Y .ne. 1

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method


		Cumulative distribution function  (P)  is calculated directly by
		code associated with the following reference.

		DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
		Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
		Trans. Math.  Softw. 18 (1993), 360-373.

		Computation of other parameters involve a seach for a value that
		produces  the desired  value  of P.   The search relies  on  the
		monotinicity of P with the other parameter.


		Note


		The beta density is proportional to
		t^(A-1) * (1-t)^(B-1)

		**********************************************************************/
	{
#define tol (1.0e-8)
#define atol (1.0e-50)
#define zero (1.0e-300)
#define inf 1.0e300
#define one 1.0e0
		static int K1 = 1;
		static double K2 = 0.0e0;
		static double K3 = 1.0e0;
		static double K8 = 0.5e0;
		static double K9 = 5.0e0;
		double v_fx,xhi,xlo,cum,ccum,xy,pq;
		unsigned long qhi,qleft,qporq;
		double T12,T13,T14,T15;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		if(!(*which < 1 || *which > 4)) goto S30;
		if(!(*which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 4.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(*which == 1) goto S70;
		/*
		P
		*/
		if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
		if(!(*p < 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = 1.0e0;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(*which == 1) goto S110;
		/*
		Q
		*/
		if(!(*q < 0.0e0 || *q > 1.0e0)) goto S100;
		if(!(*q < 0.0e0)) goto S80;
		*bound = 0.0e0;
		goto S90;
	S80:
		*bound = 1.0e0;
	S90:
		*status = -3;
		return;
	S110:
	S100:
		if(*which == 2) goto S150;
		/*
		X
		*/
		if(!(*x < 0.0e0 || *x > 1.0e0)) goto S140;
		if(!(*x < 0.0e0)) goto S120;
		*bound = 0.0e0;
		goto S130;
	S120:
		*bound = 1.0e0;
	S130:
		*status = -4;
		return;
	S150:
	S140:
		if(*which == 2) goto S190;
		/*
		Y
		*/
		if(!(*y < 0.0e0 || *y > 1.0e0)) goto S180;
		if(!(*y < 0.0e0)) goto S160;
		*bound = 0.0e0;
		goto S170;
	S160:
		*bound = 1.0e0;
	S170:
		*status = -5;
		return;
	S190:
	S180:
		if(*which == 3) goto S210;
		/*
		A
		*/
		if(!(*a <= 0.0e0)) goto S200;
		*bound = 0.0e0;
		*status = -6;
		return;
	S210:
	S200:
		if(*which == 4) goto S230;
		/*
		B
		*/
		if(!(*b <= 0.0e0)) goto S220;
		*bound = 0.0e0;
		*status = -7;
		return;
	S230:
	S220:
		if(*which == 1) goto S270;
		/*
		P + Q
		*/
		pq = *p+*q;
		if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S260;
		if(!(pq < 0.0e0)) goto S240;
		*bound = 0.0e0;
		goto S250;
	S240:
		*bound = 1.0e0;
	S250:
		*status = 3;
		return;
	S270:
	S260:
		if(*which == 2) goto S310;
		/*
		X + Y
		*/
		xy = *x+*y;
		if(!(fabs(xy-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S300;
		if(!(xy < 0.0e0)) goto S280;
		*bound = 0.0e0;
		goto S290;
	S280:
		*bound = 1.0e0;
	S290:
		*status = 4;
		return;
	S310:
	S300:
		if(!(*which == 1)) qporq = *p <= *q;
		/*
		Select the minimum of P or Q
		Calculate ANSWERS
		*/
		if(1 == *which) {
			/*
			Calculating P and Q
			*/
			x_cumbet(*x,*y,*a,*b,p,q);
			*status = 0;
		} else if(2 == *which) {
			/*
			Calculating X and Y
			*/
			pthread_mutex_lock(&H_ctrl);
			dstzr(0.0, 1.0, atol, tol);
			if(!qporq) goto S340;
			*status = 0;
			dzror(status,x,v_fx,&xlo,&xhi,&qleft,&qhi);
			*y = one-*x;
		S320:
			if(!(*status == 1)) goto S330;
			x_cumbet(*x,*y,*a,*b,&cum,&ccum);
			v_fx = cum-*p;
			dzror(status,x,v_fx,&xlo,&xhi,&qleft,&qhi);
			*y = one-*x;
			goto S320;
		S330:
			goto S370;
		S340:
			*status = 0;
			dzror(status,y,v_fx,&xlo,&xhi,&qleft,&qhi);
			*x = one-*y;
		S350:
			if(!(*status == 1)) goto S360;
			x_cumbet(*x,*y,*a,*b,&cum,&ccum);
			v_fx = ccum-*q;
			dzror(status,y,v_fx,&xlo,&xhi,&qleft,&qhi);
			*x = one-*y;
			goto S350;
		S370:
		S360:
			if(!(*status == -1)) goto S400;
			if(!qleft) goto S380;
			*status = 1;
			*bound = 0.0e0;
			goto S390;
		S380:
			*status = 2;
			*bound = 1.0e0;
		S400:
		S390:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if (3 == *which) {
			/*
			Computing A
			*/
			pthread_mutex_lock(&H_ctrl);
			*a = 5.0e0;
			dstinv(zero,inf,K8,K8,K9,atol,tol);
			*status = 0;
			dinvr(status,a,v_fx,&qleft,&qhi);
		S410:
			if(!(*status == 1)) goto S440;
			x_cumbet(*x,*y,*a,*b,&cum,&ccum);
			if(!qporq) goto S420;
			v_fx = cum-*p;
			goto S430;
		S420:
			v_fx = ccum-*q;
		S430:
			dinvr(status,a,v_fx,&qleft,&qhi);
			goto S410;
		S440:
			if(!(*status == -1)) goto S470;
			if(!qleft) goto S450;
			*status = 1;
			*bound = zero;
			goto S460;
		S450:
			*status = 2;
			*bound = inf;
		S470:
		S460:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(4 == *which) {
			/*
			Computing B
			*/
			pthread_mutex_lock(&H_ctrl);
			*b = 5.0e0;
			T12 = zero;
			T13 = inf;
			T14 = atol;
			T15 = tol;
			dstinv(zero,inf,K8,K8,K9,atol,tol);
			*status = 0;
			dinvr(status,b,v_fx,&qleft,&qhi);
		S480:
			if(!(*status == 1)) goto S510;
			x_cumbet(*x,*y,*a,*b,&cum,&ccum);
			if(!qporq) goto S490;
			v_fx = cum-*p;
			goto S500;
		S490:
			v_fx = ccum-*q;
		S500:
			dinvr(status,b,v_fx,&qleft,&qhi);
			goto S480;
		S510:
			if(!(*status == -1)) goto S540;
			if(!qleft) goto S520;
			*status = 1;
			*bound = zero;
			goto S530;
		S520:
			*status = 2;
			*bound = inf;
		S530:
			;
			pthread_mutex_unlock(&H_ctrl);\
		}
	S540:
		return;
#undef tol
#undef atol
#undef zero
#undef inf
#undef one
	}
	void cdfbin(int *which,double *p,double *q,double *s,double *xn,
		double *pr,double *ompr,int *status,double *bound)
		/**********************************************************************

		void cdfbin(int *which,double *p,double *q,double *s,double *xn,
		double *pr,double *ompr,int *status,double *bound)

		Cumulative Distribution Function
		BINomial distribution


		Function


		Calculates any one parameter of the binomial
		distribution given values for the others.


		Arguments


		WHICH --> Integer indicating which of the next four argument
		values is to be calculated from the others.
		Legal range: 1..4
		iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
		iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
		iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
		iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN

		P <--> The cumulation from 0 to S of the binomial distribution.
		(Probablility of S or fewer successes in XN trials each
		with probability of success PR.)
		Input range: [0,1].

		Q <--> 1-P.
		Input range: [0, 1].
		P + Q = 1.0.

		S <--> The number of successes observed.
		Input range: [0, XN]
		Search range: [0, XN]

		XN  <--> The number of binomial trials.
		Input range: (0, +infinity).
		Search range: [1E-300, 1E300]

		PR  <--> The probability of success in each binomial trial.
		Input range: [0,1].
		Search range: [0,1]

		OMPR  <--> 1-PR
		Input range: [0,1].
		Search range: [0,1]
		PR + OMPR = 1.0

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound
		3 if P + Q .ne. 1
		4 if PR + OMPR .ne. 1

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method


		Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of
		Mathematical   Functions (1966) is   used  to reduce the  binomial
		distribution  to  the  cumulative incomplete    beta distribution.

		Computation of other parameters involve a seach for a value that
		produces  the desired  value  of P.   The search relies  on  the
		monotinicity of P with the other parameter.


		**********************************************************************/
	{
#define atol (1.0e-50)
#define tol (1.0e-8)
#define zero (1.0e-300)
#define inf 1.0e300
#define one 1.0e0
		static int K1 = 1;
		static double K2 = 0.0e0;
		static double K3 = 0.5e0;
		static double K4 = 5.0e0;
		static double K11 = 1.0e0;
		double v_fx,xhi,xlo,cum,ccum,pq,prompr;
		unsigned long qhi,qleft,qporq;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		if(!(*which < 1 && *which > 4)) goto S30;
		if(!(*which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 4.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(*which == 1) goto S70;
		/*
		P
		*/
		if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
		if(!(*p < 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = 1.0e0;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(*which == 1) goto S110;
		/*
		Q
		*/
		if(!(*q < 0.0e0 || *q > 1.0e0)) goto S100;
		if(!(*q < 0.0e0)) goto S80;
		*bound = 0.0e0;
		goto S90;
	S80:
		*bound = 1.0e0;
	S90:
		*status = -3;
		return;
	S110:
	S100:
		if(*which == 3) goto S130;
		/*
		XN
		*/
		if(!(*xn <= 0.0e0)) goto S120;
		*bound = 0.0e0;
		*status = -5;
		return;
	S130:
	S120:
		if(*which == 2) goto S170;
		/*
		S
		*/
		if(!(*s < 0.0e0 || (*which != 3 && *s > *xn))) goto S160;
		if(!(*s < 0.0e0)) goto S140;
		*bound = 0.0e0;
		goto S150;
	S140:
		*bound = *xn;
	S150:
		*status = -4;
		return;
	S170:
	S160:
		if(*which == 4) goto S210;
		/*
		PR
		*/
		if(!(*pr < 0.0e0 || *pr > 1.0e0)) goto S200;
		if(!(*pr < 0.0e0)) goto S180;
		*bound = 0.0e0;
		goto S190;
	S180:
		*bound = 1.0e0;
	S190:
		*status = -6;
		return;
	S210:
	S200:
		if(*which == 4) goto S250;
		/*
		OMPR
		*/
		if(!(*ompr < 0.0e0 || *ompr > 1.0e0)) goto S240;
		if(!(*ompr < 0.0e0)) goto S220;
		*bound = 0.0e0;
		goto S230;
	S220:
		*bound = 1.0e0;
	S230:
		*status = -7;
		return;
	S250:
	S240:
		if(*which == 1) goto S290;
		/*
		P + Q
		*/
		pq = *p+*q;
		if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S280;
		if(!(pq < 0.0e0)) goto S260;
		*bound = 0.0e0;
		goto S270;
	S260:
		*bound = 1.0e0;
	S270:
		*status = 3;
		return;
	S290:
	S280:
		if(*which == 4) goto S330;
		/*
		PR + OMPR
		*/
		prompr = *pr+*ompr;
		if(!(fabs(prompr-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S320;
		if(!(prompr < 0.0e0)) goto S300;
		*bound = 0.0e0;
		goto S310;
	S300:
		*bound = 1.0e0;
	S310:
		*status = 4;
		return;
	S330:
	S320:
		if(!(*which == 1)) qporq = *p <= *q;
		/*
		Select the minimum of P or Q
		Calculate ANSWERS
		*/
		if(1 == *which) {
			/*
			Calculating P
			*/
			x_cumbin(*s,*xn,*pr,*ompr,p,q);
			*status = 0;
		} else if(2 == *which) {
			/*
			Calculating S
			*/
			pthread_mutex_lock(&H_ctrl);
			*s = 5.0e0;
			dstinv(K2,*xn,K3,K3,K4,atol,tol);
			*status = 0;
			dinvr(status,s,v_fx,&qleft,&qhi);
		S340:
			if(!(*status == 1)) goto S370;
			x_cumbin(*s,*xn,*pr,*ompr,&cum,&ccum);
			if(!qporq) goto S350;
			v_fx = cum-*p;
			goto S360;
		S350:
			v_fx = ccum-*q;
		S360:
			dinvr(status,s,v_fx,&qleft,&qhi);
			goto S340;
		S370:
			if(!(*status == -1)) goto S400;
			if(!qleft) goto S380;
			*status = 1;
			*bound = 0.0e0;
			goto S390;
		S380:
			*status = 2;
			*bound = *xn;
		S400:
		S390:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(3 == *which) {
			/*
			Calculating XN
			*/
			pthread_mutex_lock(&H_ctrl);
			*xn = 5.0e0;
			dstinv(zero,inf,K3,K3,K4,atol,tol);
			*status = 0;
			dinvr(status,xn,v_fx,&qleft,&qhi);
		S410:
			if(!(*status == 1)) goto S440;
			x_cumbin(*s,*xn,*pr,*ompr,&cum,&ccum);
			if(!qporq) goto S420;
			v_fx = cum-*p;
			goto S430;
		S420:
			v_fx = ccum-*q;
		S430:
			dinvr(status,xn,v_fx,&qleft,&qhi);
			goto S410;
		S440:
			if(!(*status == -1)) goto S470;
			if(!qleft) goto S450;
			*status = 1;
			*bound = zero;
			goto S460;
		S450:
			*status = 2;
			*bound = inf;
		S470:
		S460:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(4 == *which) {
			/*
			Calculating PR and OMPR
			*/
			pthread_mutex_lock(&H_ctrl);
			dstzr(K2,K11,atol,tol);
			if(!qporq) goto S500;
			*status = 0;
			dzror(status,pr,v_fx,&xlo,&xhi,&qleft,&qhi);
			*ompr = one-*pr;
		S480:
			if(!(*status == 1)) goto S490;
			x_cumbin(*s,*xn,*pr,*ompr,&cum,&ccum);
			v_fx = cum-*p;
			dzror(status,pr,v_fx,&xlo,&xhi,&qleft,&qhi);
			*ompr = one-*pr;
			goto S480;
		S490:
			goto S530;
		S500:
			*status = 0;
			dzror(status,ompr,v_fx,&xlo,&xhi,&qleft,&qhi);
			*pr = one-*ompr;
		S510:
			if(!(*status == 1)) goto S520;
			x_cumbin(*s,*xn,*pr,*ompr,&cum,&ccum);
			v_fx = ccum-*q;
			dzror(status,ompr,v_fx,&xlo,&xhi,&qleft,&qhi);
			*pr = one-*ompr;
			goto S510;
		S530:
		S520:
			if(!(*status == -1)) goto S560;
			if(!qleft) goto S540;
			*status = 1;
			*bound = 0.0e0;
			goto S550;
		S540:
			*status = 2;
			*bound = 1.0e0;
		S550:
			;
			pthread_mutex_unlock(&H_ctrl);
		}
	S560:
		return;
#undef atol
#undef tol
#undef zero
#undef inf
#undef one
	}
	void cdfchi(int *which,double *p,double *q,double *x,double *df,
		int *status,double *bound)
		/**********************************************************************

		void cdfchi(int *which,double *p,double *q,double *x,double *df,
		int *status,double *bound)

		Cumulative Distribution Function
		CHI-Square distribution


		Function


		Calculates any one parameter of the chi-square
		distribution given values for the others.


		Arguments


		WHICH --> Integer indicating which of the next three argument
		values is to be calculated from the others.
		Legal range: 1..3
		iwhich = 1 : Calculate P and Q from X and DF
		iwhich = 2 : Calculate X from P,Q and DF
		iwhich = 3 : Calculate DF from P,Q and X

		P <--> The integral from 0 to X of the chi-square
		distribution.
		Input range: [0, 1].

		Q <--> 1-P.
		Input range: (0, 1].
		P + Q = 1.0.

		X <--> Upper limit of integration of the non-central
		chi-square distribution.
		Input range: [0, +infinity).
		Search range: [0,1E300]

		DF <--> Degrees of freedom of the
		chi-square distribution.
		Input range: (0, +infinity).
		Search range: [ 1E-300, 1E300]

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound
		3 if P + Q .ne. 1
		10 indicates error returned from cumgam.  See
		references in cdfgam

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method


		Formula    26.4.19   of Abramowitz  and     Stegun, Handbook  of
		Mathematical Functions   (1966) is used   to reduce the chisqure
		distribution to the incomplete distribution.

		Computation of other parameters involve a seach for a value that
		produces  the desired  value  of P.   The search relies  on  the
		monotinicity of P with the other parameter.

		**********************************************************************/
	{
#define tol (1.0e-8)
#define atol (1.0e-50)
#define zero (1.0e-300)
#define inf 1.0e300
		int K1 = 1;
		double K2 = 0.0e0;
		double K4 = 0.5e0;
		double K5 = 5.0e0;
		double v_fx = 0,cum = 0,ccum = 0,pq = 0,porq = 0;
		unsigned long qhi = 0,qleft = 0,qporq = 0;
		double T3 = 0,T6 = 0,T7 = 0,T8 = 0,T9 = 0,T10 = 0,T11 = 0;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		if(!(*which < 1 || *which > 3)) goto S30;
		if(!(*which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 3.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(*which == 1) goto S70;
		/*
		P
		*/
		if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
		if(!(*p < 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = 1.0e0;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(*which == 1) goto S110;
		/*
		Q
		*/
		if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
		if(!(*q <= 0.0e0)) goto S80;
		*bound = 0.0e0;
		goto S90;
	S80:
		*bound = 1.0e0;
	S90:
		*status = -3;
		return;
	S110:
	S100:
		if(*which == 2) goto S130;
		/*
		X
		*/
		if(!(*x < 0.0e0)) goto S120;
		*bound = 0.0e0;
		*status = -4;
		return;
	S130:
	S120:
		if(*which == 3) goto S150;
		/*
		DF
		*/
		if(!(*df <= 0.0e0)) goto S140;
		*bound = 0.0e0;
		*status = -5;
		return;
	S150:
	S140:
		if(*which == 1) goto S190;
		/*
		P + Q
		*/
		pq = *p+*q;
		if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S180;
		if(!(pq < 0.0e0)) goto S160;
		*bound = 0.0e0;
		goto S170;
	S160:
		*bound = 1.0e0;
	S170:
		*status = 3;
		return;
	S190:
	S180:
		if(*which == 1) goto S220;
		/*
		Select the minimum of P or Q
		*/
		qporq = *p <= *q;
		if(!qporq) goto S200;
		porq = *p;
		goto S210;
	S200:
		porq = *q;
	S220:
	S210:
		/*
		Calculate ANSWERS
		*/
		if(1 == *which) {
			/*
			Calculating P and Q
			*/
			*status = 0;
			x_cumchi(*x,*df,p,q);
			if(porq > 1.5e0) {
				*status = 10;
				return;
			}
		} else if(2 == *which) {
			/*
			Calculating X
			*/
			pthread_mutex_lock(&H_ctrl);
			*x = 5.0e0;
			dstinv(0.0,inf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,x,v_fx,&qleft,&qhi);
		S230:
			if(!(*status == 1)) goto S270;
			x_cumchi(*x,*df,&cum,&ccum);
			if(!qporq) goto S240;
			v_fx = cum-*p;
			goto S250;
		S240:
			v_fx = ccum-*q;
		S250:
			if(!(v_fx+porq > 1.5e0)) goto S260;
			*status = 10;
			return;
		S260:
			dinvr(status,x,v_fx,&qleft,&qhi);
			goto S230;
		S270:
			if(!(*status == -1)) goto S300;
			if(!qleft) goto S280;
			*status = 1;
			*bound = 0.0e0;
			goto S290;
		S280:
			*status = 2;
			*bound = inf;
		S300:
		S290:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(3 == *which) {
			/*
			Calculating DF
			*/
			pthread_mutex_lock(&H_ctrl);
			*df = 5.0e0;
			dstinv(zero,inf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,df,v_fx,&qleft,&qhi);
		S310:
			if(!(*status == 1)) goto S350;
			x_cumchi(*x,*df,&cum,&ccum);
			if(!qporq) goto S320;
			v_fx = cum-*p;
			goto S330;
		S320:
			v_fx = ccum-*q;
		S330:
			if(!(v_fx+porq > 1.5e0)) goto S340;
			*status = 10;
			return;
		S340:
			dinvr(status,df,v_fx,&qleft,&qhi);
			goto S310;
		S350:
			if(!(*status == -1)) goto S380;
			if(!qleft) goto S360;
			*status = 1;
			*bound = zero;
			goto S370;
		S360:
			*status = 2;
			*bound = inf;
		S370:
			;
			pthread_mutex_unlock(&H_ctrl);
		}
	S380:
		return;
#undef tol
#undef atol
#undef zero
#undef inf
	}
	void cdfchn(int *which,double *p,double *q,double *x,double *df,
		double *pnonc,int *status,double *bound)
		/**********************************************************************

		void cdfchn(int *which,double *p,double *q,double *x,double *df,
		double *pnonc,int *status,double *bound)

		Cumulative Distribution Function
		Non-central Chi-Square


		Function


		Calculates any one parameter of the non-central chi-square
		distribution given values for the others.


		Arguments


		WHICH --> Integer indicating which of the next three argument
		values is to be calculated from the others.
		Input range: 1..4
		iwhich = 1 : Calculate P and Q from X and DF
		iwhich = 2 : Calculate X from P,DF and PNONC
		iwhich = 3 : Calculate DF from P,X and PNONC
		iwhich = 3 : Calculate PNONC from P,X and DF

		P <--> The integral from 0 to X of the non-central chi-square
		distribution.
		Input range: [0, 1-1E-16).

		Q <--> 1-P.
		Q is not used by this subroutine and is only included
		for similarity with other cdf* routines.

		X <--> Upper limit of integration of the non-central
		chi-square distribution.
		Input range: [0, +infinity).
		Search range: [0,1E300]

		DF <--> Degrees of freedom of the non-central
		chi-square distribution.
		Input range: (0, +infinity).
		Search range: [ 1E-300, 1E300]

		PNONC <--> Non-centrality parameter of the non-central
		chi-square distribution.
		Input range: [0, +infinity).
		Search range: [0,1E4]

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method


		Formula  26.4.25   of   Abramowitz   and   Stegun,  Handbook  of
		Mathematical  Functions (1966) is used to compute the cumulative
		distribution function.

		Computation of other parameters involve a seach for a value that
		produces  the desired  value  of P.   The search relies  on  the
		monotinicity of P with the other parameter.


		WARNING

		The computation time  required for this  routine is proportional
		to the noncentrality  parameter  (PNONC).  Very large  values of
		this parameter can consume immense  computer resources.  This is
		why the search range is bounded by 10,000.

		**********************************************************************/
	{
#define tent4 1.0e4
#define tol (1.0e-8)
#define atol (1.0e-50)
#define zero (1.0e-300)
#define one (1.0e0-1.0e-16)
#define inf 1.0e300
		static double K1 = 0.0e0;
		static double K3 = 0.5e0;
		static double K4 = 5.0e0;
		double v_fx,cum,ccum;
		unsigned long qhi,qleft;
		double T7,T8,T9,T10;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		if(!(*which < 1 || *which > 4)) goto S30;
		if(!(*which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 4.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(*which == 1) goto S70;
		/*
		P
		*/
		if(!(*p < 0.0e0 || *p > one)) goto S60;
		if(!(*p < 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = one;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(*which == 2) goto S90;
		/*
		X
		*/
		if(!(*x < 0.0e0)) goto S80;
		*bound = 0.0e0;
		*status = -4;
		return;
	S90:
	S80:
		if(*which == 3) goto S110;
		/*
		DF
		*/
		if(!(*df <= 0.0e0)) goto S100;
		*bound = 0.0e0;
		*status = -5;
		return;
	S110:
	S100:
		if(*which == 4) goto S130;
		/*
		PNONC
		*/
		if(!(*pnonc < 0.0e0)) goto S120;
		*bound = 0.0e0;
		*status = -6;
		return;
	S130:
	S120:
		/*
		Calculate ANSWERS
		*/
		if(1 == *which) {
			/*
			Calculating P and Q
			*/
			cumchn(*x,*df,*pnonc,p,q);
			*status = 0;
		} else if(2 == *which) {
			/*
			Calculating X
			*/
			pthread_mutex_lock(&H_ctrl);
			*x = 5.0e0;
			dstinv(K1,inf,K3,K3,K4,atol,tol);
			*status = 0;
			dinvr(status,x,v_fx,&qleft,&qhi);
		S140:
			if(!(*status == 1)) goto S150;
			cumchn(*x,*df,*pnonc,&cum,&ccum);
			v_fx = cum-*p;
			dinvr(status,x,v_fx,&qleft,&qhi);
			goto S140;
		S150:
			if(!(*status == -1)) goto S180;
			if(!qleft) goto S160;
			*status = 1;
			*bound = 0.0e0;
			goto S170;
		S160:
			*status = 2;
			*bound = inf;
		S180:
		S170:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(3 == *which) {
			/*
			Calculating DF
			*/
			pthread_mutex_lock(&H_ctrl);
			*df = 5.0e0;
			T7 = zero;
			T8 = inf;
			T9 = atol;
			T10 = tol;
			dstinv(zero,inf,K3,K3,K4,atol,tol);
			*status = 0;
			dinvr(status,df,v_fx,&qleft,&qhi);
		S190:
			if(!(*status == 1)) goto S200;
			cumchn(*x,*df,*pnonc,&cum,&ccum);
			v_fx = cum-*p;
			dinvr(status,df,v_fx,&qleft,&qhi);
			goto S190;
		S200:
			if(!(*status == -1)) goto S230;
			if(!qleft) goto S210;
			*status = 1;
			*bound = zero;
			goto S220;
		S210:
			*status = 2;
			*bound = inf;
		S230:
		S220:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(4 == *which) {
			/*
			Calculating PNONC
			*/
			pthread_mutex_lock(&H_ctrl);
			*pnonc = 5.0;
			dstinv(K1,tent4,K3,K3,K4,atol,tol);
			*status = 0;
			dinvr(status,pnonc,v_fx,&qleft,&qhi);
		S240:
			if(!(*status == 1)) goto S250;
			cumchn(*x,*df,*pnonc,&cum,&ccum);
			v_fx = cum-*p;
			dinvr(status,pnonc,v_fx,&qleft,&qhi);
			goto S240;
		S250:
			if(!(*status == -1)) goto S280;
			if(!qleft) goto S260;
			*status = 1;
			*bound = zero;
			goto S270;
		S260:
			*status = 2;
			*bound = tent4;
		S270:
			;
			pthread_mutex_unlock(&H_ctrl);
		}
	S280:
		return;
#undef tent4
#undef tol
#undef atol
#undef zero
#undef one
#undef inf
	}
	void cdff(int *which,double *p,double *q,double *f,double *dfn,
		double *dfd,int *status,double *bound)
		/**********************************************************************

		void cdff(int *which,double *p,double *q,double *f,double *dfn,
		double *dfd,int *status,double *bound)

		Cumulative Distribution Function
		F distribution


		Function


		Calculates any one parameter of the F distribution
		given values for the others.


		Arguments


		WHICH --> Integer indicating which of the next four argument
		values is to be calculated from the others.
		Legal range: 1..4
		iwhich = 1 : Calculate P and Q from F,DFN and DFD
		iwhich = 2 : Calculate F from P,Q,DFN and DFD
		iwhich = 3 : Calculate DFN from P,Q,F and DFD
		iwhich = 4 : Calculate DFD from P,Q,F and DFN

		P <--> The integral from 0 to F of the f-density.
		Input range: [0,1].

		Q <--> 1-P.
		Input range: (0, 1].
		P + Q = 1.0.

		F <--> Upper limit of integration of the f-density.
		Input range: [0, +infinity).
		Search range: [0,1E300]

		DFN < --> Degrees of freedom of the numerator sum of squares.
		Input range: (0, +infinity).
		Search range: [ 1E-300, 1E300]

		DFD < --> Degrees of freedom of the denominator sum of squares.
		Input range: (0, +infinity).
		Search range: [ 1E-300, 1E300]

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound
		3 if P + Q .ne. 1

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method


		Formula   26.6.2   of   Abramowitz   and   Stegun,  Handbook  of
		Mathematical  Functions (1966) is used to reduce the computation
		of the  cumulative  distribution function for the  F  variate to
		that of an incomplete beta.

		Computation of other parameters involve a seach for a value that
		produces  the desired  value  of P.   The search relies  on  the
		monotinicity of P with the other parameter.

		WARNING

		The value of the  cumulative  F distribution is  not necessarily
		monotone in  either degrees of freedom.  There  thus may  be two
		values  that  provide a given CDF  value.   This routine assumes
		monotonicity and will find an arbitrary one of the two values.

		**********************************************************************/
	{
#define tol (1.0e-8)
#define atol (1.0e-50)
#define zero (1.0e-300)
#define inf 1.0e300
		static int K1 = 1;
		static double K2 = 0.0e0;
		static double K4 = 0.5e0;
		static double K5 = 5.0e0;
		double pq,v_fx,cum,ccum;
		unsigned long qhi,qleft,qporq;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		if(!(*which < 1 || *which > 4)) goto S30;
		if(!(*which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 4.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(*which == 1) goto S70;
		/*
		P
		*/
		if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
		if(!(*p < 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = 1.0e0;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(*which == 1) goto S110;
		/*
		Q
		*/
		if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
		if(!(*q <= 0.0e0)) goto S80;
		*bound = 0.0e0;
		goto S90;
	S80:
		*bound = 1.0e0;
	S90:
		*status = -3;
		return;
	S110:
	S100:
		if(*which == 2) goto S130;
		/*
		F
		*/
		if(!(*f < 0.0e0)) goto S120;
		*bound = 0.0e0;
		*status = -4;
		return;
	S130:
	S120:
		if(*which == 3) goto S150;
		/*
		DFN
		*/
		if(!(*dfn <= 0.0e0)) goto S140;
		*bound = 0.0e0;
		*status = -5;
		return;
	S150:
	S140:
		if(*which == 4) goto S170;
		/*
		DFD
		*/
		if(!(*dfd <= 0.0e0)) goto S160;
		*bound = 0.0e0;
		*status = -6;
		return;
	S170:
	S160:
		if(*which == 1) goto S210;
		/*
		P + Q
		*/
		pq = *p+*q;
		if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S200;
		if(!(pq < 0.0e0)) goto S180;
		*bound = 0.0e0;
		goto S190;
	S180:
		*bound = 1.0e0;
	S190:
		*status = 3;
		return;
	S210:
	S200:
		if(!(*which == 1)) qporq = *p <= *q;
		/*
		Select the minimum of P or Q
		Calculate ANSWERS
		*/
		if(1 == *which) {
			/*
			Calculating P
			*/
			x_cumf(*f,*dfn,*dfd,p,q);
			*status = 0;
		} else if(2 == *which) {
			/*
			Calculating F
			*/
			pthread_mutex_lock(&H_ctrl);
			*f = 5.0;
			dstinv(K2,inf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,f,v_fx,&qleft,&qhi);
		S220:
			if(!(*status == 1)) goto S250;
			x_cumf(*f,*dfn,*dfd,&cum,&ccum);
			if(!qporq) goto S230;
			v_fx = cum-*p;
			goto S240;
		S230:
			v_fx = ccum-*q;
		S240:
			dinvr(status,f,v_fx,&qleft,&qhi);
			goto S220;
		S250:
			if(!(*status == -1)) goto S280;
			if(!qleft) goto S260;
			*status = 1;
			*bound = 0.0e0;
			goto S270;
		S260:
			*status = 2;
			*bound = inf;
		S280:
		S270:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(3 == *which) {
			/*
			Calculating DFN
			*/
			pthread_mutex_lock(&H_ctrl);
			*dfn = 5.0;
			dstinv(zero,inf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,dfn,v_fx,&qleft,&qhi);
		S290:
			if(!(*status == 1)) goto S320;
			x_cumf(*f,*dfn,*dfd,&cum,&ccum);
			if(!qporq) goto S300;
			v_fx = cum-*p;
			goto S310;
		S300:
			v_fx = ccum-*q;
		S310:
			dinvr(status,dfn,v_fx,&qleft,&qhi);
			goto S290;
		S320:
			if(!(*status == -1)) goto S350;
			if(!qleft) goto S330;
			*status = 1;
			*bound = zero;
			goto S340;
		S330:
			*status = 2;
			*bound = inf;
		S350:
		S340:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(4 == *which) {
			/*
			Calculating DFD
			*/
			pthread_mutex_lock(&H_ctrl);
			*dfd = 5.0;
			dstinv(zero,inf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,dfd,v_fx,&qleft,&qhi);
		S360:
			if(!(*status == 1)) goto S390;
			x_cumf(*f,*dfn,*dfd,&cum,&ccum);
			if(!qporq) goto S370;
			v_fx = cum-*p;
			goto S380;
		S370:
			v_fx = ccum-*q;
		S380:
			dinvr(status,dfd,v_fx,&qleft,&qhi);
			goto S360;
		S390:
			if(!(*status == -1)) goto S420;
			if(!qleft) goto S400;
			*status = 1;
			*bound = zero;
			goto S410;
		S400:
			*status = 2;
			*bound = inf;
		S410:
			;
			pthread_mutex_unlock(&H_ctrl);
		}
	S420:
		return;
#undef tol
#undef atol
#undef zero
#undef inf
	}
	void cdffnc(int *which,double *p,double *q,double *f,double *dfn,
		double *dfd,double *phonc,int *status,double *bound)
		/**********************************************************************

		void cdffnc(int *which,double *p,double *q,double *f,double *dfn,
		double *dfd,double *phonc,int *status,double *bound)

		Cumulative Distribution Function
		Non-central F distribution


		Function


		Calculates any one parameter of the Non-central F
		distribution given values for the others.


		Arguments


		WHICH --> Integer indicating which of the next five argument
		values is to be calculated from the others.
		Legal range: 1..5
		iwhich = 1 : Calculate P and Q from F,DFN,DFD and PNONC
		iwhich = 2 : Calculate F from P,Q,DFN,DFD and PNONC
		iwhich = 3 : Calculate DFN from P,Q,F,DFD and PNONC
		iwhich = 4 : Calculate DFD from P,Q,F,DFN and PNONC
		iwhich = 5 : Calculate PNONC from P,Q,F,DFN and DFD

		P <--> The integral from 0 to F of the non-central f-density.
		Input range: [0,1-1E-16).

		Q <--> 1-P.
		Q is not used by this subroutine and is only included
		for similarity with other cdf* routines.

		F <--> Upper limit of integration of the non-central f-density.
		Input range: [0, +infinity).
		Search range: [0,1E300]

		DFN < --> Degrees of freedom of the numerator sum of squares.
		Input range: (0, +infinity).
		Search range: [ 1E-300, 1E300]

		DFD < --> Degrees of freedom of the denominator sum of squares.
		Must be in range: (0, +infinity).
		Input range: (0, +infinity).
		Search range: [ 1E-300, 1E300]

		PNONC <-> The non-centrality parameter
		Input range: [0,infinity)
		Search range: [0,1E4]

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound
		3 if P + Q .ne. 1

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method


		Formula  26.6.20   of   Abramowitz   and   Stegun,  Handbook  of
		Mathematical  Functions (1966) is used to compute the cumulative
		distribution function.

		Computation of other parameters involve a seach for a value that
		produces  the desired  value  of P.   The search relies  on  the
		monotinicity of P with the other parameter.

		WARNING

		The computation time  required for this  routine is proportional
		to the noncentrality  parameter  (PNONC).  Very large  values of
		this parameter can consume immense  computer resources.  This is
		why the search range is bounded by 10,000.

		WARNING

		The  value  of the  cumulative  noncentral F distribution is not
		necessarily monotone in either degrees  of freedom.  There  thus
		may be two values that provide a given  CDF value.  This routine
		assumes monotonicity  and will find  an arbitrary one of the two
		values.

		**********************************************************************/
	{
#define tent4 1.0e4
#define tol (1.0e-8)
#define atol (1.0e-50)
#define zero (1.0e-300)
#define one (1.0e0-1.0e-16)
#define inf 1.0e300
		static double K1 = 0.0e0;
		static double K3 = 0.5e0;
		static double K4 = 5.0e0;
		double v_fx,cum,ccum;
		unsigned long qhi,qleft;
		double T7,T8,T9,T10;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		if(!(*which < 1 || *which > 5)) goto S30;
		if(!(*which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 5.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(*which == 1) goto S70;
		/*
		P
		*/
		if(!(*p < 0.0e0 || *p > one)) goto S60;
		if(!(*p < 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = one;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(*which == 2) goto S90;
		/*
		F
		*/
		if(!(*f < 0.0e0)) goto S80;
		*bound = 0.0e0;
		*status = -4;
		return;
	S90:
	S80:
		if(*which == 3) goto S110;
		/*
		DFN
		*/
		if(!(*dfn <= 0.0e0)) goto S100;
		*bound = 0.0e0;
		*status = -5;
		return;
	S110:
	S100:
		if(*which == 4) goto S130;
		/*
		DFD
		*/
		if(!(*dfd <= 0.0e0)) goto S120;
		*bound = 0.0e0;
		*status = -6;
		return;
	S130:
	S120:
		if(*which == 5) goto S150;
		/*
		PHONC
		*/
		if(!(*phonc < 0.0e0)) goto S140;
		*bound = 0.0e0;
		*status = -7;
		return;
	S150:
	S140:
		/*
		Calculate ANSWERS
		*/
		if(1 == *which) {
			/*
			Calculating P
			*/
			cumfnc(f,dfn,dfd,phonc,p,q);
			*status = 0;
		} else if(2 == *which) {
			/*
			Calculating F
			*/
			pthread_mutex_lock(&H_ctrl);
			*f = 5.0;
			dstinv(K1,inf,K3,K3,K4,atol,tol);
			*status = 0;
			dinvr(status,f,v_fx,&qleft,&qhi);
		S160:
			if(!(*status == 1)) goto S170;
			cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
			v_fx = cum-*p;
			dinvr(status,f,v_fx,&qleft,&qhi);
			goto S160;
		S170:
			if(!(*status == -1)) goto S200;
			if(!qleft) goto S180;
			*status = 1;
			*bound = 0.0e0;
			goto S190;
		S180:
			*status = 2;
			*bound = inf;
		S200:
		S190:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(3 == *which) {
			/*
			Calculating DFN
			*/
			pthread_mutex_lock(&H_ctrl);
			*dfn = 5.0e0;
			T7 = zero;
			T8 = inf;
			T9 = atol;
			T10 = tol;
			dstinv(zero,inf,K3,K3,K4,atol,tol);
			*status = 0;
			dinvr(status,dfn,v_fx,&qleft,&qhi);
		S210:
			if(!(*status == 1)) goto S220;
			cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
			v_fx = cum-*p;
			dinvr(status,dfn,v_fx,&qleft,&qhi);
			goto S210;
		S220:
			if(!(*status == -1)) goto S250;
			if(!qleft) goto S230;
			*status = 1;
			*bound = zero;
			goto S240;
		S230:
			*status = 2;
			*bound = inf;
		S250:
		S240:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(4 == *which) {
			/*
			Calculating DFD
			*/
			pthread_mutex_lock(&H_ctrl);
			*dfd = 5.0;
			dstinv(zero,inf,K3,K3,K4,atol,tol);
			*status = 0;
			dinvr(status,dfd,v_fx,&qleft,&qhi);
		S260:
			if(!(*status == 1)) goto S270;
			cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
			v_fx = cum-*p;
			dinvr(status,dfd,v_fx,&qleft,&qhi);
			goto S260;
		S270:
			if(!(*status == -1)) goto S300;
			if(!qleft) goto S280;
			*status = 1;
			*bound = zero;
			goto S290;
		S280:
			*status = 2;
			*bound = inf;
		S300:
		S290:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(5 == *which) {
			/*
			Calculating PHONC
			*/
			pthread_mutex_lock(&H_ctrl);
			*phonc = 5.0;
			dstinv(K1,tent4,K3,K3,K4,atol,tol);
			*status = 0;
			dinvr(status,phonc,v_fx,&qleft,&qhi);
		S310:
			if(!(*status == 1)) goto S320;
			cumfnc(f,dfn,dfd,phonc,&cum,&ccum);
			v_fx = cum-*p;
			dinvr(status,phonc,v_fx,&qleft,&qhi);
			goto S310;
		S320:
			if(!(*status == -1)) goto S350;
			if(!qleft) goto S330;
			*status = 1;
			*bound = 0.0e0;
			goto S340;
		S330:
			*status = 2;
			*bound = tent4;
		S340:
			;
			pthread_mutex_unlock(&H_ctrl);
		}
	S350:
		return;
#undef tent4
#undef tol
#undef atol
#undef zero
#undef one
#undef inf
	}
	void cdfgam(int *which,double *p,double *q,double *x,double *shape,
		double *scale,int *status,double *bound)
		/**********************************************************************

		void cdfgam(int *which,double *p,double *q,double *x,double *shape,
		double *scale,int *status,double *bound)

		Cumulative Distribution Function
		GAMma Distribution


		Function


		Calculates any one parameter of the gamma
		distribution given values for the others.


		Arguments


		WHICH --> Integer indicating which of the next four argument
		values is to be calculated from the others.
		Legal range: 1..4
		iwhich = 1 : Calculate P and Q from X,SHAPE and SCALE
		iwhich = 2 : Calculate X from P,Q,SHAPE and SCALE
		iwhich = 3 : Calculate SHAPE from P,Q,X and SCALE
		iwhich = 4 : Calculate SCALE from P,Q,X and SHAPE

		P <--> The integral from 0 to X of the gamma density.
		Input range: [0,1].

		Q <--> 1-P.
		Input range: (0, 1].
		P + Q = 1.0.

		X <--> The upper limit of integration of the gamma density.
		Input range: [0, +infinity).
		Search range: [0,1E300]

		SHAPE <--> The shape parameter of the gamma density.
		Input range: (0, +infinity).
		Search range: [1E-300,1E300]

		SCALE <--> The scale parameter of the gamma density.
		Input range: (0, +infinity).
		Search range: (1E-300,1E300]

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound
		3 if P + Q .ne. 1
		10 if the gamma or inverse gamma routine cannot
		compute the answer.  Usually happens only for
		X and SHAPE very large (gt 1E10 or more)

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method


		Cumulative distribution function (P) is calculated directly by
		the code associated with:

		DiDinato, A. R. and Morris, A. H. Computation of the  incomplete
		gamma function  ratios  and their  inverse.   ACM  Trans.  Math.
		Softw. 12 (1986), 377-393.

		Computation of other parameters involve a seach for a value that
		produces  the desired  value  of P.   The search relies  on  the
		monotinicity of P with the other parameter.


		Note



		The gamma density is proportional to
		T**(SHAPE - 1) * EXP(- SCALE * T)

		**********************************************************************/
	{
#define tol (1.0e-8)
#define atol (1.0e-50)
#define zero (1.0e-300)
#define inf 1.0e300
		static int K1 = 1;
		static double K5 = 0.5e0;
		static double K6 = 5.0e0;
		double xx,v_fx,xscale,cum,ccum,pq,porq;
		int ierr;
		unsigned long qhi,qleft,qporq;
		double T2,T9;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		if(!(*which < 1 || *which > 4)) goto S30;
		if(!(*which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 4.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(*which == 1) goto S70;
		/*
		P
		*/
		if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
		if(!(*p < 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = 1.0e0;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(*which == 1) goto S110;
		/*
		Q
		*/
		if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
		if(!(*q <= 0.0e0)) goto S80;
		*bound = 0.0e0;
		goto S90;
	S80:
		*bound = 1.0e0;
	S90:
		*status = -3;
		return;
	S110:
	S100:
		if(*which == 2) goto S130;
		/*
		X
		*/
		if(!(*x < 0.0e0)) goto S120;
		*bound = 0.0e0;
		*status = -4;
		return;
	S130:
	S120:
		if(*which == 3) goto S150;
		/*
		SHAPE
		*/
		if(!(*shape <= 0.0e0)) goto S140;
		*bound = 0.0e0;
		*status = -5;
		return;
	S150:
	S140:
		if(*which == 4) goto S170;
		/*
		SCALE
		*/
		if(!(*scale <= 0.0e0)) goto S160;
		*bound = 0.0e0;
		*status = -6;
		return;
	S170:
	S160:
		if(*which == 1) goto S210;
		/*
		P + Q
		*/
		pq = *p+*q;
		if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S200;
		if(!(pq < 0.0e0)) goto S180;
		*bound = 0.0e0;
		goto S190;
	S180:
		*bound = 1.0e0;
	S190:
		*status = 3;
		return;
	S210:
	S200:
		if(*which == 1) goto S240;
		/*
		Select the minimum of P or Q
		*/
		qporq = *p <= *q;
		if(!qporq) goto S220;
		porq = *p;
		goto S230;
	S220:
		porq = *q;
	S240:
	S230:
		/*
		Calculate ANSWERS
		*/
		if(1 == *which) {
			/*
			Calculating P
			*/
			*status = 0;
			xscale = *x**scale;
			x_cumgam(&xscale,shape,p,q);
			if(porq > 1.5e0) *status = 10;
		} else if(2 == *which) {
			/*
			Computing X
			*/
			pthread_mutex_lock(&H_ctrl);
			T2 = -1.0e0;
			x_gaminv(shape,&xx,&T2,p,q,&ierr);
			if(ierr < 0.0e0) {
				*status = 10;
				return;
			}
			else  {
				*x = xx/ *scale;
				*status = 0;
			}
			pthread_mutex_unlock(&H_ctrl);
		} else if(3 == *which) {
			/*
			Computing SHAPE
			*/
			pthread_mutex_lock(&H_ctrl);
			*shape = 5.0;
			xscale = *x**scale;
			dstinv(zero,inf,K5,K5,K6,atol,tol);
			*status = 0;
			dinvr(status,shape,v_fx,&qleft,&qhi);
		S250:
			if(!(*status == 1)) goto S290;
			x_cumgam(&xscale,shape,&cum,&ccum);
			if(!qporq) goto S260;
			v_fx = cum-*p;
			goto S270;
		S260:
			v_fx = ccum-*q;
		S270:
			if(!((qporq && cum > 1.5e0) || (!qporq && ccum > 1.5e0))) goto S280;
			*status = 10;
			return;
		S280:
			dinvr(status,shape,v_fx,&qleft,&qhi);
			goto S250;
		S290:
			if(!(*status == -1)) goto S320;
			if(!qleft) goto S300;
			*status = 1;
			*bound = zero;
			goto S310;
		S300:
			*status = 2;
			*bound = inf;
		S320:
		S310:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(4 == *which) {
			/*
			Computing SCALE
			*/
			T9 = -1.0e0;
			x_gaminv(shape,&xx,&T9,p,q,&ierr);
			if(ierr < 0.0e0) {
				*status = 10;
				return;
			}
			else  {
				*scale = xx/ *x;
				*status = 0;
			}
		}
		return;
#undef tol
#undef atol
#undef zero
#undef inf
	}
	void cdfnbn(int *which,double *p,double *q,double *s,double *xn,
		double *pr,double *ompr,int *status,double *bound)
		/**********************************************************************

		void cdfnbn(int *which,double *p,double *q,double *s,double *xn,
		double *pr,double *ompr,int *status,double *bound)

		Cumulative Distribution Function
		Negative BiNomial distribution


		Function


		Calculates any one parameter of the negative binomial
		distribution given values for the others.

		The  cumulative  negative   binomial  distribution  returns  the
		probability that there  will be  F or fewer failures before  the
		XNth success in binomial trials each of which has probability of
		success PR.

		The individual term of the negative binomial is the probability of
		S failures before XN successes and is
		Choose( S, XN+S-1 ) * PR^(XN) * (1-PR)^S


		Arguments


		WHICH --> Integer indicating which of the next four argument
		values is to be calculated from the others.
		Legal range: 1..4
		iwhich = 1 : Calculate P and Q from S,XN,PR and OMPR
		iwhich = 2 : Calculate S from P,Q,XN,PR and OMPR
		iwhich = 3 : Calculate XN from P,Q,S,PR and OMPR
		iwhich = 4 : Calculate PR and OMPR from P,Q,S and XN

		P <--> The cumulation from 0 to S of the  negative
		binomial distribution.
		Input range: [0,1].

		Q <--> 1-P.
		Input range: (0, 1].
		P + Q = 1.0.

		S <--> The upper limit of cumulation of the binomial distribution.
		There are F or fewer failures before the XNth success.
		Input range: [0, +infinity).
		Search range: [0, 1E300]

		XN  <--> The number of successes.
		Input range: [0, +infinity).
		Search range: [0, 1E300]

		PR  <--> The probability of success in each binomial trial.
		Input range: [0,1].
		Search range: [0,1].

		OMPR  <--> 1-PR
		Input range: [0,1].
		Search range: [0,1]
		PR + OMPR = 1.0

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound
		3 if P + Q .ne. 1
		4 if PR + OMPR .ne. 1

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method


		Formula   26.5.26   of   Abramowitz  and  Stegun,  Handbook   of
		Mathematical Functions (1966) is used  to  reduce calculation of
		the cumulative distribution  function to that of  an  incomplete
		beta.

		Computation of other parameters involve a seach for a value that
		produces  the desired  value  of P.   The search relies  on  the
		monotinicity of P with the other parameter.

		**********************************************************************/
	{
#define tol (1.0e-8)
#define atol (1.0e-50)
#define inf 1.0e300
#define one 1.0e0
		static int K1 = 1;
		static double K2 = 0.0e0;
		static double K4 = 0.5e0;
		static double K5 = 5.0e0;
		static double K11 = 1.0e0;
		double v_fx,xhi,xlo,pq,prompr,cum,ccum;
		unsigned long qhi,qleft,qporq;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		if(!(*which < 1 || *which > 4)) goto S30;
		if(!(*which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 4.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(*which == 1) goto S70;
		/*
		P
		*/
		if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
		if(!(*p < 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = 1.0e0;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(*which == 1) goto S110;
		/*
		Q
		*/
		if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
		if(!(*q <= 0.0e0)) goto S80;
		*bound = 0.0e0;
		goto S90;
	S80:
		*bound = 1.0e0;
	S90:
		*status = -3;
		return;
	S110:
	S100:
		if(*which == 2) goto S130;
		/*
		S
		*/
		if(!(*s < 0.0e0)) goto S120;
		*bound = 0.0e0;
		*status = -4;
		return;
	S130:
	S120:
		if(*which == 3) goto S150;
		/*
		XN
		*/
		if(!(*xn < 0.0e0)) goto S140;
		*bound = 0.0e0;
		*status = -5;
		return;
	S150:
	S140:
		if(*which == 4) goto S190;
		/*
		PR
		*/
		if(!(*pr < 0.0e0 || *pr > 1.0e0)) goto S180;
		if(!(*pr < 0.0e0)) goto S160;
		*bound = 0.0e0;
		goto S170;
	S160:
		*bound = 1.0e0;
	S170:
		*status = -6;
		return;
	S190:
	S180:
		if(*which == 4) goto S230;
		/*
		OMPR
		*/
		if(!(*ompr < 0.0e0 || *ompr > 1.0e0)) goto S220;
		if(!(*ompr < 0.0e0)) goto S200;
		*bound = 0.0e0;
		goto S210;
	S200:
		*bound = 1.0e0;
	S210:
		*status = -7;
		return;
	S230:
	S220:
		if(*which == 1) goto S270;
		/*
		P + Q
		*/
		pq = *p+*q;
		if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S260;
		if(!(pq < 0.0e0)) goto S240;
		*bound = 0.0e0;
		goto S250;
	S240:
		*bound = 1.0e0;
	S250:
		*status = 3;
		return;
	S270:
	S260:
		if(*which == 4) goto S310;
		/*
		PR + OMPR
		*/
		prompr = *pr+*ompr;
		if(!(fabs(prompr-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S300;
		if(!(prompr < 0.0e0)) goto S280;
		*bound = 0.0e0;
		goto S290;
	S280:
		*bound = 1.0e0;
	S290:
		*status = 4;
		return;
	S310:
	S300:
		if(!(*which == 1)) qporq = *p <= *q;
		/*
		Select the minimum of P or Q
		Calculate ANSWERS
		*/
		if(1 == *which) {
			/*
			Calculating P
			*/
			x_cumnbn(*s,*xn,*pr,*ompr,p,q);
			*status = 0;
		} else if(2 == *which) {
			/*
			Calculating S
			*/
			pthread_mutex_lock(&H_ctrl);
			*s = 5.0;
			dstinv(K2,inf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,s,v_fx,&qleft,&qhi);
		S320:
			if(!(*status == 1)) goto S350;
			x_cumnbn(*s,*xn,*pr,*ompr,&cum,&ccum);
			if(!qporq) goto S330;
			v_fx = cum-*p;
			goto S340;
		S330:
			v_fx = ccum-*q;
		S340:
			dinvr(status,s,v_fx,&qleft,&qhi);
			goto S320;
		S350:
			if(!(*status == -1)) goto S380;
			if(!qleft) goto S360;
			*status = 1;
			*bound = 0.0e0;
			goto S370;
		S360:
			*status = 2;
			*bound = inf;
		S380:
		S370:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(3 == *which) {
			/*
			Calculating XN
			*/
			pthread_mutex_lock(&H_ctrl);
			*xn = 5.0;
			dstinv(K2,inf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,xn,v_fx,&qleft,&qhi);
		S390:
			if(!(*status == 1)) goto S420;
			x_cumnbn(*s,*xn,*pr,*ompr,&cum,&ccum);
			if(!qporq) goto S400;
			v_fx = cum-*p;
			goto S410;
		S400:
			v_fx = ccum-*q;
		S410:
			dinvr(status,xn,v_fx,&qleft,&qhi);
			goto S390;
		S420:
			if(!(*status == -1)) goto S450;
			if(!qleft) goto S430;
			*status = 1;
			*bound = 0.0e0;
			goto S440;
		S430:
			*status = 2;
			*bound = inf;
		S450:
		S440:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(4 == *which) {
			/*
			Calculating PR and OMPR
			*/
			pthread_mutex_lock(&H_ctrl);
			dstzr(K2, K11, atol, tol);
			if(!qporq) goto S480;
			*status = 0;
			dzror(status,pr,v_fx,&xlo,&xhi,&qleft,&qhi);
			*ompr = one-*pr;
		S460:
			if(!(*status == 1)) goto S470;
			x_cumnbn(*s,*xn,*pr,*ompr,&cum,&ccum);
			v_fx = cum-*p;
			dzror(status,pr,v_fx,&xlo,&xhi,&qleft,&qhi);
			*ompr = one-*pr;
			goto S460;
		S470:
			goto S510;
		S480:
			*status = 0;
			dzror(status,ompr,v_fx,&xlo,&xhi,&qleft,&qhi);
			*pr = one-*ompr;
		S490:
			if(!(*status == 1)) goto S500;
			x_cumnbn(*s,*xn,*pr,*ompr,&cum,&ccum);
			v_fx = ccum-*q;
			dzror(status,ompr,v_fx,&xlo,&xhi,&qleft,&qhi);
			*pr = one-*ompr;
			goto S490;
		S510:
		S500:
			if(!(*status == -1)) goto S540;
			if(!qleft) goto S520;
			*status = 1;
			*bound = 0.0e0;
			goto S530;
		S520:
			*status = 2;
			*bound = 1.0e0;
		S530:
			;
			pthread_mutex_unlock(&H_ctrl);
		}
	S540:
		return;
#undef tol
#undef atol
#undef inf
#undef one
	}
	void x_cdfnor(int which,double *p,double *q,double *x,double *mean,
		double *sd,int *status,double *bound)
		/**********************************************************************

		void cdfnor(int which,double p,double q,double *x,double *mean,
		double *sd,int *status,double *bound)

		Cumulative Distribution Function
		NORmal distribution


		Function


		Calculates any one parameter of the normal
		distribution given values for the others.


		Arguments


		WHICH  --> Integer indicating  which of the  next  parameter
		values is to be calculated using values  of the others.
		Legal range: 1..4
		iwhich = 1 : Calculate P and Q from X,MEAN and SD
		iwhich = 2 : Calculate X from P,Q,MEAN and SD
		iwhich = 3 : Calculate MEAN from P,Q,X and SD
		iwhich = 4 : Calculate SD from P,Q,X and MEAN

		P <--> The integral from -infinity to X of the normal density.
		Input range: (0,1].

		Q <--> 1-P.
		Input range: (0, 1].
		P + Q = 1.0.

		X < --> Upper limit of integration of the normal-density.
		Input range: ( -infinity, +infinity)

		MEAN <--> The mean of the normal density.
		Input range: (-infinity, +infinity)

		SD <--> Standard Deviation of the normal density.
		Input range: (0, +infinity).

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound
		3 if P + Q .ne. 1

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method




		A slightly modified version of ANORM from

		Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
		Package of Special Function Routines and Test Drivers"
		acm Transactions on Mathematical Software. 19, 22-32.

		is used to calulate the  cumulative standard normal distribution.

		The rational functions from pages  90-95  of Kennedy and Gentle,
		Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
		starting values to Newton's Iterations which compute the inverse
		standard normal.  Therefore no  searches  are necessary for  any
		parameter.

		For X < -15, the asymptotic expansion for the normal is used  as
		the starting value in finding the inverse standard normal.
		This is formula 26.2.12 of Abramowitz and Stegun.


		Note


		The normal density is proportional to
		exp( - 0.5 * (( X - MEAN)/SD)**2)

		**********************************************************************/
	{
		double z;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		*status = 0;
		if(!(which < 1 || which > 4)) goto S30;
		if(!(which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 4.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(which == 1) goto S70;
		/*
		P
		*/
		if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
		if(!(*p <= 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = 1.0e0;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(which == 1) goto S110;
		/*
		Q
		*/
		if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
		if(!(*q <= 0.0e0)) goto S80;
		*bound = 0.0e0;
		goto S90;
	S80:
		*bound = 1.0e0;
	S90:
		*status = -3;
		return;
	S110:
	S100:
		double pq = W0;
		if(which == 1) goto S150;
		/*
		P + Q
		*/
		pq = *p + *q;
		if(!(fabs(pq-1.0) > 3.0* x_spmpar(1))) goto S140;
		if(!(pq < 0.0)) goto S120;
		*bound = 0.0;
		goto S130;
	S120:
		*bound = 1.0;
	S130:
		*status = 3;
		return;
	S150:
	S140:
		if(which == 4) goto S170;
		/*
		SD
		*/
		if(!(*sd <= 0.0)) goto S160;
		*bound = 0.0;
		*status = -6;
		return;
	S170:
	S160:
		/*
		Calculate ANSWERS
		*/
		if(1 == which) {
			/*
			Computing P
			*/
			z = (*x-*mean)/ *sd;
			x_cumnor(&z,p,q);
		} else if(2 == which) {
			/*
			Computing X
			*/
			z = x_dblNormalDistInv(*p,*q);
			*x = *sd*z+*mean;
		} else if(3 == which) {
			/*
			Computing the MEAN
			*/
			z = x_dblNormalDistInv(*p,*q);
			*mean = *x-*sd*z;
		} else if(4 == which) {
			/*
			Computing SD
			*/
			z = x_dblNormalDistInv(*p,*q);
			*sd = (*x-*mean)/z;
		}
		return;
	}
	void cdfpoi(int *which,double *p,double *q,double *s,double *xlam,
		int *status,double *bound)
		/**********************************************************************

		void cdfpoi(int *which,double *p,double *q,double *s,double *xlam,
		int *status,double *bound)

		Cumulative Distribution Function
		POIsson distribution


		Function


		Calculates any one parameter of the Poisson
		distribution given values for the others.


		Arguments


		WHICH --> Integer indicating which  argument
		value is to be calculated from the others.
		Legal range: 1..3
		iwhich = 1 : Calculate P and Q from S and XLAM
		iwhich = 2 : Calculate A from P,Q and XLAM
		iwhich = 3 : Calculate XLAM from P,Q and S

		P <--> The cumulation from 0 to S of the poisson density.
		Input range: [0,1].

		Q <--> 1-P.
		Input range: (0, 1].
		P + Q = 1.0.

		S <--> Upper limit of cumulation of the Poisson.
		Input range: [0, +infinity).
		Search range: [0,1E300]

		XLAM <--> Mean of the Poisson distribution.
		Input range: [0, +infinity).
		Search range: [0,1E300]

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound
		3 if P + Q .ne. 1

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method


		Formula   26.4.21  of   Abramowitz  and   Stegun,   Handbook  of
		Mathematical Functions (1966) is used  to reduce the computation
		of  the cumulative distribution function to that  of computing a
		chi-square, hence an incomplete gamma function.

		Cumulative  distribution function  (P) is  calculated  directly.
		Computation of other parameters involve a seach for a value that
		produces  the desired value of  P.   The  search relies  on  the
		monotinicity of P with the other parameter.

		**********************************************************************/
	{
#define tol (1.0e-8)
#define atol (1.0e-50)
#define inf 1.0e300
		static int K1 = 1;
		static double K2 = 0.0e0;
		static double K4 = 0.5e0;
		static double K5 = 5.0e0;
		double v_fx,cum,ccum,pq;
		unsigned long qhi,qleft,qporq;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		if(!(*which < 1 || *which > 3)) goto S30;
		if(!(*which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 3.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(*which == 1) goto S70;
		/*
		P
		*/
		if(!(*p < 0.0e0 || *p > 1.0e0)) goto S60;
		if(!(*p < 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = 1.0e0;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(*which == 1) goto S110;
		/*
		Q
		*/
		if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
		if(!(*q <= 0.0e0)) goto S80;
		*bound = 0.0e0;
		goto S90;
	S80:
		*bound = 1.0e0;
	S90:
		*status = -3;
		return;
	S110:
	S100:
		if(*which == 2) goto S130;
		/*
		S
		*/
		if(!(*s < 0.0e0)) goto S120;
		*bound = 0.0e0;
		*status = -4;
		return;
	S130:
	S120:
		if(*which == 3) goto S150;
		/*
		XLAM
		*/
		if(!(*xlam < 0.0e0)) goto S140;
		*bound = 0.0e0;
		*status = -5;
		return;
	S150:
	S140:
		if(*which == 1) goto S190;
		/*
		P + Q
		*/
		pq = *p+*q;
		if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S180;
		if(!(pq < 0.0e0)) goto S160;
		*bound = 0.0e0;
		goto S170;
	S160:
		*bound = 1.0e0;
	S170:
		*status = 3;
		return;
	S190:
	S180:
		if(!(*which == 1)) qporq = *p <= *q;
		/*
		Select the minimum of P or Q
		Calculate ANSWERS
		*/
		if(1 == *which) {
			/*
			Calculating P
			*/
			x_cumpoi(*s,*xlam,p,q);
			*status = 0;
		} else if(2 == *which) {
			/*
			Calculating S
			*/
			pthread_mutex_lock(&H_ctrl);
			*s = 5.0;
			dstinv(K2,inf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,s,v_fx,&qleft,&qhi);
		S200:
			if(!(*status == 1)) goto S230;
			x_cumpoi(*s,*xlam,&cum,&ccum);
			if(!qporq) goto S210;
			v_fx = cum-*p;
			goto S220;
		S210:
			v_fx = ccum-*q;
		S220:
			dinvr(status,s,v_fx,&qleft,&qhi);
			goto S200;
		S230:
			if(!(*status == -1)) goto S260;
			if(!qleft) goto S240;
			*status = 1;
			*bound = 0.0e0;
			goto S250;
		S240:
			*status = 2;
			*bound = inf;
		S260:
		S250:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(3 == *which) {
			/*
			Calculating XLAM
			*/
			pthread_mutex_lock(&H_ctrl);
			*xlam = 5.0;
			dstinv(K2,inf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,xlam,v_fx,&qleft,&qhi);
		S270:
			if(!(*status == 1)) goto S300;
			x_cumpoi(*s,*xlam,&cum,&ccum);
			if(!qporq) goto S280;
			v_fx = cum-*p;
			goto S290;
		S280:
			v_fx = ccum-*q;
		S290:
			dinvr(status,xlam,v_fx,&qleft,&qhi);
			goto S270;
		S300:
			if(!(*status == -1)) goto S330;
			if(!qleft) goto S310;
			*status = 1;
			*bound = 0.0e0;
			goto S320;
		S310:
			*status = 2;
			*bound = inf;
		S320:
			;
			pthread_mutex_unlock(&H_ctrl);
		}
	S330:
		return;
#undef tol
#undef atol
#undef inf
	}

	void cdft(int *which,double *p,double *q,double *t,double *df,
		int *status,double *bound)
		/**********************************************************************

		void cdft(int *which,double *p,double *q,double *t,double *df,
		int *status,double *bound)

		Cumulative Distribution Function
		T distribution


		Function


		Calculates any one parameter of the t distribution given
		values for the others.


		Arguments


		WHICH --> Integer indicating which  argument
		values is to be calculated from the others.
		Legal range: 1..3
		iwhich = 1 : Calculate P and Q from T and DF
		iwhich = 2 : Calculate T from P,Q and DF
		iwhich = 3 : Calculate DF from P,Q and T

		P <--> The integral from -infinity to t of the t-density.
		Input range: (0,1].

		Q <--> 1-P.
		Input range: (0, 1].
		P + Q = 1.0.

		T <--> Upper limit of integration of the t-density.
		Input range: ( -infinity, +infinity).
		Search range: [ -1E300, 1E300 ]

		DF <--> Degrees of freedom of the t-distribution.
		Input range: (0 , +infinity).
		Search range: [1e-300, 1E10]

		STATUS <-- 0 if calculation completed correctly
		-I if input parameter number I is out of range
		1 if answer appears to be lower than lowest
		search bound
		2 if answer appears to be higher than greatest
		search bound
		3 if P + Q .ne. 1

		BOUND <-- Undefined if STATUS is 0

		Bound exceeded by parameter number I if STATUS
		is negative.

		Lower search bound if STATUS is 1.

		Upper search bound if STATUS is 2.


		Method


		Formula  26.5.27  of   Abramowitz   and  Stegun,   Handbook   of
		Mathematical Functions  (1966) is used to reduce the computation
		of the cumulative distribution function to that of an incomplete
		beta.

		Computation of other parameters involve a seach for a value that
		produces  the desired  value  of P.   The search relies  on  the
		monotinicity of P with the other parameter.

		**********************************************************************/
	{
#define tol (1.0e-8)
#define atol (1.0e-50)
#define zero (1.0e-300)
#define inf 1.0e300
#define maxdf 1.0e10
		static int K1 = 1;
		static double K4 = 0.5e0;
		static double K5 = 5.0e0;
		double v_fx,cum,ccum,pq;
		unsigned long qhi,qleft,qporq;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Check arguments
		*/
		if(!(*which < 1 || *which > 3)) goto S30;
		if(!(*which < 1)) goto S10;
		*bound = 1.0e0;
		goto S20;
	S10:
		*bound = 3.0e0;
	S20:
		*status = -1;
		return;
	S30:
		if(*which == 1) goto S70;
		/*
		P
		*/
		if(!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
		if(!(*p <= 0.0e0)) goto S40;
		*bound = 0.0e0;
		goto S50;
	S40:
		*bound = 1.0e0;
	S50:
		*status = -2;
		return;
	S70:
	S60:
		if(*which == 1) goto S110;
		/*
		Q
		*/
		if(!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
		if(!(*q <= 0.0e0)) goto S80;
		*bound = 0.0e0;
		goto S90;
	S80:
		*bound = 1.0e0;
	S90:
		*status = -3;
		return;
	S110:
	S100:
		if(*which == 3) goto S130;
		/*
		DF
		*/
		if(!(*df <= 0.0e0)) goto S120;
		*bound = 0.0e0;
		*status = -5;
		return;
	S130:
	S120:
		if(*which == 1) goto S170;
		/*
		P + Q
		*/
		pq = *p+*q;
		if(!(fabs(pq-0.5e0-0.5e0) > 3.0e0*x_spmpar(K1))) goto S160;
		if(!(pq < 0.0e0)) goto S140;
		*bound = 0.0e0;
		goto S150;
	S140:
		*bound = 1.0e0;
	S150:
		*status = 3;
		return;
	S170:
	S160:
		if(!(*which == 1)) qporq = *p <= *q;
		/*
		Select the minimum of P or Q
		Calculate ANSWERS
		*/
		if(1 == *which) {
			/*
			Computing P and Q
			*/
			x_cumt(t,df,p,q);
			*status = 0;
		} else if(2 == *which) {
			/*
			Computing T
			.. Get initial approximation for T
			*/
			pthread_mutex_lock(&H_ctrl);
			*t = x_dt1(*p, *q, *df);
			dstinv(-inf,inf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,t,v_fx,&qleft,&qhi);
		S180:
			if(!(*status == 1)) goto S210;
			x_cumt(t,df,&cum,&ccum);
			if(!qporq) goto S190;
			v_fx = cum-*p;
			goto S200;
		S190:
			v_fx = ccum-*q;
		S200:
			dinvr(status,t,v_fx,&qleft,&qhi);
			goto S180;
		S210:
			if(!(*status == -1)) goto S240;
			if(!qleft) goto S220;
			*status = 1;
			*bound = -inf;
			goto S230;
		S220:
			*status = 2;
			*bound = inf;
		S240:
		S230:
			;
			pthread_mutex_unlock(&H_ctrl);
		} else if(3 == *which) {
			/*
			Computing DF
			*/
			pthread_mutex_lock(&H_ctrl);
			*df = 5.0;
			dstinv(zero,maxdf,K4,K4,K5,atol,tol);
			*status = 0;
			dinvr(status,df,v_fx,&qleft,&qhi);
		S250:
			if(!(*status == 1)) goto S280;
			x_cumt(t,df,&cum,&ccum);
			if(!qporq) goto S260;
			v_fx = cum-*p;
			goto S270;
		S260:
			v_fx = ccum-*q;
		S270:
			dinvr(status,df,v_fx,&qleft,&qhi);
			goto S250;
		S280:
			if(!(*status == -1)) goto S310;
			if(!qleft) goto S290;
			*status = 1;
			*bound = zero;
			goto S300;
		S290:
			*status = 2;
			*bound = maxdf;
		S300:
			;
			pthread_mutex_unlock(&H_ctrl);
		}
	S310:
		return;
#undef tol
#undef atol
#undef zero
#undef inf
#undef maxdf
	}
	void x_cumbet(double x,double y,double a,double b,double *cum,
		double *ccum)
		/*
		**********************************************************************

		void cumbet(double *x,double *y,double *a,double *b,double *cum,
		double *ccum)

		Double precision cUMulative incomplete BETa distribution


		Function


		Calculates the cdf to X of the incomplete beta distribution
		with parameters a and b.  This is the integral from 0 to x
		of (1/B(a,b))*f(t)) where f(t) = t**(a-1) * (1-t)**(b-1)


		Arguments


		X --> Upper limit of integration.
		X is DOUBLE PRECISION

		Y --> 1 - X.
		Y is DOUBLE PRECISION

		A --> First parameter of the beta distribution.
		A is DOUBLE PRECISION

		B --> Second parameter of the beta distribution.
		B is DOUBLE PRECISION

		CUM <-- Cumulative incomplete beta distribution.
		CUM is DOUBLE PRECISION

		CCUM <-- Compliment of Cumulative incomplete beta distribution.
		CCUM is DOUBLE PRECISION


		Method


		Calls the routine BRATIO.

		References

		Didonato, Armido R. and Morris, Alfred H. Jr. (1992) Algorithim
		708 Significant Digit Computation of the Incomplete Beta Function
		Ratios. ACM ToMS, Vol.18, No. 3, Sept. 1992, 360-373.

		**********************************************************************
		*/
	{
		int ierr;
		/*
		..
		.. Executable Statements ..
		*/
		if(!(x <= 0.0)) goto S10;
		*cum = 0.0;
		*ccum = 1.0;
		return;
	S10:
		if(!(y <= 0.0)) goto S20;
		*cum = 1.0;
		*ccum = 0.0;
		return;
	S20:
		// Call bratio routine
		x_bratio(a,b,x,y,cum,ccum,&ierr);
		return;
	}
	void x_cumbin(double s,double xn,double pr,double ompr,
		double *cum,double *ccum)
		/*
		**********************************************************************

		void cumbin(double *s,double *xn,double *pr,double *ompr,
		double *cum,double *ccum)

		CUmulative BINomial distribution


		Function


		Returns the probability   of 0  to  S  successes in  XN   binomial
		trials, each of which has a probability of success, PBIN.


		Arguments


		S --> The upper limit of cumulation of the binomial distribution.
		S is DOUBLE PRECISION

		XN --> The number of binomial trials.
		XN is DOUBLE PRECISIO

		PBIN --> The probability of success in each binomial trial.
		PBIN is DOUBLE PRECIS

		OMPR --> 1 - PBIN
		OMPR is DOUBLE PRECIS

		CUM <-- Cumulative binomial distribution.
		CUM is DOUBLE PRECISI

		CCUM <-- Compliment of Cumulative binomial distribution.
		CCUM is DOUBLE PRECIS


		Method


		Formula  26.5.24    of   Abramowitz  and    Stegun,  Handbook   of
		Mathematical   Functions (1966) is   used  to reduce the  binomial
		distribution  to  the  cumulative    beta distribution.

		**********************************************************************
		*/
	{
		/*
		..
		.. Executable Statements ..
		*/
		if(!(s < xn)) goto S10;
		x_cumbet(pr,ompr,s+1.0,xn-s,ccum,cum);
		goto S20;
	S10:
		*cum = 1.0;
		*ccum = 0.0;
	S20:
		return;
	}
	void x_cumchi(double x,double df,double *cum,double *ccum)
		/*
		**********************************************************************

		void cumchi(double *x,double *df,double *cum,double *ccum)
		CUMulative of the CHi-square distribution


		Function


		Calculates the cumulative chi-square distribution.


		Arguments


		X       --> Upper limit of integration of the
		chi-square distribution.
		X is DOUBLE PRECISION

		DF      --> Degrees of freedom of the
		chi-square distribution.
		DF is DOUBLE PRECISION

		CUM <-- Cumulative chi-square distribution.
		CUM is DOUBLE PRECISIO

		CCUM <-- Compliment of Cumulative chi-square distribution.
		CCUM is DOUBLE PRECISI


		Method


		Calls incomplete gamma function (CUMGAM)

		**********************************************************************
		*/
	{
		double a,xx;
		/*
		..
		.. Executable Statements ..
		*/
		a = df*0.5;
		xx = x*0.5;
		x_cumgam(&xx,&a,cum,ccum);
		return;
	}
	void cumchn(double x,double df,double pnonc,double *cum,
		double *ccum)
		/*
		**********************************************************************

		void cumchn(double *x,double *df,double *pnonc,double *cum,
		double *ccum)

		CUMulative of the Non-central CHi-square distribution


		Function


		Calculates     the       cumulative      non-central    chi-square
		distribution, i.e.,  the probability   that  a   random   variable
		which    follows  the  non-central chi-square  distribution,  with
		non-centrality  parameter    PNONC  and   continuous  degrees   of
		freedom DF, is less than or equal to X.


		Arguments


		X       --> Upper limit of integration of the non-central
		chi-square distribution.
		X is DOUBLE PRECISION

		DF      --> Degrees of freedom of the non-central
		chi-square distribution.
		DF is DOUBLE PRECISION

		PNONC   --> Non-centrality parameter of the non-central
		chi-square distribution.
		PNONC is DOUBLE PRECIS

		CUM <-- Cumulative non-central chi-square distribution.
		CUM is DOUBLE PRECISIO

		CCUM <-- Compliment of Cumulative non-central chi-square distribut
		CCUM is DOUBLE PRECISI


		Method


		Uses  formula  26.4.25   of  Abramowitz  and  Stegun, Handbook  of
		Mathematical    Functions,  US   NBS   (1966)    to calculate  the
		non-central chi-square.


		Variables


		EPS     --- Convergence criterion.  The sum stops when a
		term is less than EPS*SUM.
		EPS is DOUBLE PRECISIO

		NTIRED  --- Maximum number of terms to be evaluated
		in each sum.
		NTIRED is INTEGER

		QCONV   --- .TRUE. if convergence achieved -
		i.e., program did not stop on NTIRED criterion.
		QCONV is LOGICAL

		CCUM <-- Compliment of Cumulative non-central
		chi-square distribution.
		CCUM is DOUBLE PRECISI

		**********************************************************************
		*/
	{
#define dg(i) (df+2.0*(double)(i))
#define qsmall(xx) (int)(sum < 1.0e-20 || (xx) < eps*sum)
#define qtired(i) (int)((i) > ntired)
		static double eps = 1.0e-5;
		int ntired = 1000;
		double adj,centaj,centwt,chid2,dfd2,lcntaj,lcntwt,lfact,pcent,pterm,sum,
			sumadj,term,wt,xnonc;
		int i,icent,iterb,iterf;
		double T1,T3;
		/*
		..
		.. Executable Statements ..
		*/
		if(!(x <= 0.0)) goto S10;
		*cum = 0.0e0;
		*ccum = 1.0e0;
		return;
	S10:
		if(!(pnonc <= 1.0e-10)) goto S20;
		/*
		When non-centrality parameter is (essentially) zero,
		use cumulative chi-square distribution
		*/
		x_cumchi(x,df,cum,ccum);
		return;
	S20:
		xnonc = pnonc/2.0;
		/*
		**********************************************************************
		The following code calcualtes the weight, chi-square, and
		adjustment term for the central term in the infinite series.
		The central term is the one in which the poisson weight is
		greatest.  The adjustment term is the amount that must
		be subtracted from the chi-square to move up two degrees
		of freedom.
		**********************************************************************
		*/
		icent = x_fifidint(xnonc);
		if(icent == 0) icent = 1;
		chid2 = x/2.0;
		/*
		Calculate central weight term
		*/
		T1 = (double)(icent+1);
		lfact = x_alngam(&T1);
		lcntwt = -xnonc+(double)icent*log(xnonc)-lfact;
		centwt = exp(lcntwt);
		/*
		Calculate central chi-square
		*/
		x_cumchi(x,dg(icent),&pcent,ccum);
		/*
		Calculate central adjustment term
		*/
		dfd2 = dg(icent)/2.0;
		T3 = 1.0e0+dfd2;
		lfact = x_alngam(&T3);
		lcntaj = dfd2*log(chid2)-chid2-lfact;
		centaj = exp(lcntaj);
		sum = centwt*pcent;
		/*
		**********************************************************************
		Sum backwards from the central term towards zero.
		Quit whenever either
		(1) the zero term is reached, or
		(2) the term gets small relative to the sum, or
		(3) More than NTIRED terms are totaled.
		**********************************************************************
		*/
		iterb = 0;
		sumadj = 0.0e0;
		adj = centaj;
		wt = centwt;
		i = icent;
		goto S40;
	S30:
		if(qtired(iterb) || qsmall(term) || i == 0) goto S50;
	S40:
		dfd2 = dg(i)/2.0e0;
		/*
		Adjust chi-square for two fewer degrees of freedom.
		The adjusted value ends up in PTERM.
		*/
		adj = adj*dfd2/chid2;
		sumadj += adj;
		pterm = pcent+sumadj;
		/*
		Adjust poisson weight for J decreased by one
		*/
		wt *= ((double)i/xnonc);
		term = wt*pterm;
		sum += term;
		i -= 1;
		iterb += 1;
		goto S30;
	S50:
		iterf = 0;
		/*
		**********************************************************************
		Now sum forward from the central term towards infinity.
		Quit when either
		(1) the term gets small relative to the sum, or
		(2) More than NTIRED terms are totaled.
		**********************************************************************
		*/
		sumadj = adj = centaj;
		wt = centwt;
		i = icent;
		goto S70;
	S60:
		if(qtired(iterf) || qsmall(term)) goto S80;
	S70:
		/*
		Update weights for next higher J
		*/
		wt *= (xnonc/(double)(i+1));
		/*
		Calculate PTERM and add term to sum
		*/
		pterm = pcent-sumadj;
		term = wt*pterm;
		sum += term;
		/*
		Update adjustment term for DF for next iteration
		*/
		i += 1;
		dfd2 = dg(i)/2.0e0;
		adj = adj*chid2/dfd2;
		sumadj += adj;
		iterf += 1;
		goto S60;
	S80:
		*cum = sum;
		*ccum = 0.5e0+(0.5e0-*cum);
		return;
#undef dg
#undef qsmall
#undef qtired
	}
	void x_cumf(double f,double dfn,double dfd,double *cum,double *ccum)
		/*
		**********************************************************************

		void cumf(double *f,double *dfn,double *dfd,double *cum,double *ccum)
		CUMulative F distribution


		Function


		Computes  the  integral from  0  to  F of  the f-density  with DFN
		and DFD degrees of freedom.


		Arguments


		F --> Upper limit of integration of the f-density.
		F is DOUBLE PRECISION

		DFN --> Degrees of freedom of the numerator sum of squares.
		DFN is DOUBLE PRECISI

		DFD --> Degrees of freedom of the denominator sum of squares.
		DFD is DOUBLE PRECISI

		CUM <-- Cumulative f distribution.
		CUM is DOUBLE PRECISI

		CCUM <-- Compliment of Cumulative f distribution.
		CCUM is DOUBLE PRECIS


		Method


		Formula  26.5.28 of  Abramowitz and   Stegun   is  used to  reduce
		the cumulative F to a cumulative beta distribution.


		Note


		If F is less than or equal to 0, 0 is returned.

		**********************************************************************
		*/
	{
#define half 0.5e0
#define done 1.0e0
		int ierr;
		/*
		..
		.. Executable Statements ..
		*/
		if(!(f <= 0.0)) goto S10;
		*cum = 0.0;
		*ccum = 1.0;
		return;
	S10:
		double prod = dfn*f;
		/*
		XX is such that the incomplete beta with parameters
		DFD/2 and DFN/2 evaluated at XX is 1 - CUM or CCUM
		YY is 1 - XX
		Calculate the smaller of XX and YY accurately
		*/
		double dsum = dfd+prod;
		double xx = dfd/dsum, yy;
		if(xx > half) {
			yy = prod/dsum;
			xx = done-yy;
		} else yy = done-xx;
		x_bratio(dfd*half,dfn*half,xx,yy,ccum,cum,&ierr);
		return;
#undef half
#undef done
	}
	void cumfnc(double *f,double *dfn,double *dfd,double *pnonc,
		double *cum,double *ccum)
		/*
		**********************************************************************

		F -NON- -C-ENTRAL F DISTRIBUTION



		Function


		COMPUTES NONCENTRAL F DISTRIBUTION WITH DFN AND DFD
		DEGREES OF FREEDOM AND NONCENTRALITY PARAMETER PNONC


		Arguments


		X --> UPPER LIMIT OF INTEGRATION OF NONCENTRAL F IN EQUATION

		DFN --> DEGREES OF FREEDOM OF NUMERATOR

		DFD -->  DEGREES OF FREEDOM OF DENOMINATOR

		PNONC --> NONCENTRALITY PARAMETER.

		CUM <-- CUMULATIVE NONCENTRAL F DISTRIBUTION

		CCUM <-- COMPLIMENT OF CUMMULATIVE


		Method


		USES FORMULA 26.6.20 OF REFERENCE FOR INFINITE SERIES.
		SERIES IS CALCULATED BACKWARD AND FORWARD FROM J = LAMBDA/2
		(THIS IS THE TERM WITH THE LARGEST POISSON WEIGHT) UNTIL
		THE CONVERGENCE CRITERION IS MET.

		FOR SPEED, THE INCOMPLETE BETA FUNCTIONS ARE EVALUATED
		BY FORMULA 26.5.16.


		REFERENCE


		HANDBOOD OF MATHEMATICAL FUNCTIONS
		EDITED BY MILTON ABRAMOWITZ AND IRENE A. STEGUN
		NATIONAL BUREAU OF STANDARDS APPLIED MATEMATICS SERIES - 55
		MARCH 1965
		P 947, EQUATIONS 26.6.17, 26.6.18


		Note


		THE SUM CONTINUES UNTIL A SUCCEEDING TERM IS LESS THAN EPS
		TIMES THE SUM (OR THE SUM IS LESS THAN 1.0E-20).  EPS IS
		SET TO 1.0E-4 IN A DATA STATEMENT WHICH CAN BE CHANGED.

		**********************************************************************
		*/
	{
#define qsmall(x) (int)(sum < 1.0e-20 || (x) < eps*sum)
#define half 0.5e0
#define done 1.0e0
		static double eps = 1.0e-4;
		double dsum,dummy,prod,xx,yy,adn,aup,b,betdn,betup,centwt,dnterm,sum,
			upterm,xmult,xnonc;
		int i,icent,ierr;
		double T1,T4,T5,T6;
		/*
		..
		.. Executable Statements ..
		*/
		if(!(*f <= 0.0e0)) goto S10;
		*cum = 0.0e0;
		*ccum = 1.0e0;
		return;
	S10:
		if(!(*pnonc < 1.0e-10)) goto S20;
		/*
		Handle case in which the non-centrality parameter is
		(essentially) zero.
		*/
		x_cumf(*f,*dfn,*dfd,cum,ccum);
		return;
	S20:
		xnonc = *pnonc/2.0e0;
		/*
		Calculate the central term of the poisson weighting factor.
		*/
		icent = (int)xnonc;
		if(icent == 0) icent = 1;
		/*
		Compute central weight term
		*/
		T1 = (double)(icent+1);
		centwt = exp(-xnonc+(double)icent*log(xnonc)-x_alngam(&T1));
		/*
		Compute central incomplete beta term
		Assure that minimum of arg to beta and 1 - arg is computed
		accurately.
		*/
		prod = *dfn**f;
		dsum = *dfd+prod;
		yy = *dfd/dsum;
		if(yy > half) {
			xx = prod/dsum;
			yy = done-xx;
		}
		else  xx = done-yy;
		x_bratio(*dfn*half+(double)icent,
			*dfd*half,
			xx,yy,&betdn,&dummy,&ierr);
		adn = *dfn/2.0e0+(double)icent;
		aup = adn;
		b = *dfd/2.0e0;
		betup = betdn;
		sum = centwt*betdn;
		/*
		Now sum terms backward from icent until convergence or all done
		*/
		xmult = centwt;
		i = icent;
		T4 = adn+b;
		T5 = adn+1.0e0;
		dnterm = exp(x_alngam(&T4)-x_alngam(&T5)-x_alngam(&b)+adn*log(xx)+b*log(yy));
	S30:
		if(qsmall(xmult*betdn) || i <= 0) goto S40;
		xmult *= ((double)i/xnonc);
		i -= 1;
		adn -= 1.0;
		dnterm = (adn+1.0)/((adn+b)*xx)*dnterm;
		betdn += dnterm;
		sum += (xmult*betdn);
		goto S30;
	S40:
		i = icent+1;
		/*
		Now sum forwards until convergence
		*/
		xmult = centwt;
		if(aup-1.0+b == 0) upterm = exp(-x_alngam(&aup)-x_alngam(&b)+(aup-1.0)*log(xx)+
			b*log(yy));
		else  {
			T6 = aup-1.0+b;
			upterm = exp(x_alngam(&T6)-x_alngam(&aup)-x_alngam(&b)+(aup-1.0)*log(xx)+b*
				log(yy));
		}
		goto S60;
	S50:
		if(qsmall(xmult*betup)) goto S70;
	S60:
		xmult *= (xnonc/(double)i);
		i += 1;
		aup += 1.0;
		upterm = (aup+b-2.0e0)*xx/(aup-1.0)*upterm;
		betup -= upterm;
		sum += (xmult*betup);
		goto S50;
	S70:
		*cum = sum;
		*ccum = 0.5e0+(0.5e0-*cum);
		return;
#undef qsmall
#undef half
#undef done
	}
	void x_cumgam(double *x,double *a,double *cum,double *ccum)
		/*
		**********************************************************************

		void cumgam(double *x,double *a,double *cum,double *ccum)
		Double precision cUMulative incomplete GAMma distribution


		Function


		Computes   the  cumulative        of    the     incomplete   gamma
		distribution, i.e., the integral from 0 to X of
		(1/GAM(A))*EXP(-T)*T**(A-1) DT
		where GAM(A) is the complete gamma function of A, i.e.,
		GAM(A) = integral from 0 to infinity of
		EXP(-T)*T**(A-1) DT


		Arguments


		X --> The upper limit of integration of the incomplete gamma.
		X is DOUBLE PRECISION

		A --> The shape parameter of the incomplete gamma.
		A is DOUBLE PRECISION

		CUM <-- Cumulative incomplete gamma distribution.
		CUM is DOUBLE PRECISION

		CCUM <-- Compliment of Cumulative incomplete gamma distribution.
		CCUM is DOUBLE PRECISIO


		Method


		Calls the routine GRATIO.

		**********************************************************************
		*/
	{
		int K1 = 0;
		/*
		..
		.. Executable Statements ..
		*/
		if(!(*x <= 0.0e0)) goto S10;
		*cum = 0.0e0;
		*ccum = 1.0e0;
		return;
	S10:
		x_gratio(a,x,cum,ccum, &K1);
		/*
		Call gratio routine
		*/
		return;
	}
	void x_cumnbn(double s,double xn,double pr,double ompr,
		double *cum,double *ccum)
		/*
		**********************************************************************

		void cumnbn(double *s,double *xn,double *pr,double *ompr,
		double *cum,double *ccum)

		CUmulative Negative BINomial distribution


		Function


		Returns the probability that it there will be S or fewer failures
		before there are XN successes, with each binomial trial having
		a probability of success PR.

		Prob(# failures = S | XN successes, PR)  =
		( XN + S - 1 )
		(            ) * PR^XN * (1-PR)^S
		(      S     )


		Arguments


		S --> The number of failures
		S is DOUBLE PRECISION

		XN --> The number of successes
		XN is DOUBLE PRECISIO

		PR --> The probability of success in each binomial trial.
		PR is DOUBLE PRECISIO

		OMPR --> 1 - PR
		OMPR is DOUBLE PRECIS

		CUM <-- Cumulative negative binomial distribution.
		CUM is DOUBLE PRECISI

		CCUM <-- Compliment of Cumulative negative binomial distribution.
		CCUM is DOUBLE PRECIS


		Method


		Formula  26.5.26    of   Abramowitz  and    Stegun,  Handbook   of
		Mathematical   Functions (1966) is   used  to reduce the  negative
		binomial distribution to the cumulative beta distribution.

		**********************************************************************
		*/
	{
		/*
		..
		.. Executable Statements ..
		*/
		x_cumbet(pr,ompr,xn,s+1.,cum,ccum);
		return;
	}
	void x_cumnor(double *arg,double *result,double *ccum)
		/*
		**********************************************************************

		void cumnor(double *arg,double *result,double *ccum)


		Function


		Computes the cumulative  of    the  normal   distribution,   i.e.,
		the integral from -infinity to x of
		(1/sqrt(2*pi)) exp(-u*u/2) du

		X --> Upper limit of integration.
		X is DOUBLE PRECISION

		RESULT <-- Cumulative normal distribution.
		RESULT is DOUBLE PRECISION

		CCUM <-- Compliment of Cumulative normal distribution.
		CCUM is DOUBLE PRECISION

		Renaming of function ANORM from:

		Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
		Package of Special Function Routines and Test Drivers"
		acm Transactions on Mathematical Software. 19, 22-32.

		with slight modifications to return ccum and to deal with
		machine constants.

		**********************************************************************
		Original Comments:
		------------------------------------------------------------------

		This function evaluates the normal distribution function:

		/ x
		1       |       -t*t/2
		P(x) = ----------- |      e       dt
		sqrt(2 pi)  |
		/-oo

		The main computation evaluates near-minimax approximations
		derived from those in "Rational Chebyshev approximations for
		the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
		This transportable program uses rational functions that
		theoretically approximate the normal distribution function to
		at least 18 significant decimal digits.  The accuracy achieved
		depends on the arithmetic system, the compiler, the intrinsic
		functions, and proper selection of the machine-dependent
		constants.

		*******************************************************************
		*******************************************************************

		Explanation of machine-dependent constants.

		MIN   = smallest machine representable number.

		EPS   = argument below which anorm(x) may be represented by
		0.5  and above which  x*x  will not underflow.
		A conservative value is the largest machine number X
		such that   1.0 + X = 1.0   to machine precision.
		*******************************************************************
		*******************************************************************

		Error returns

		The program returns  ANORM = 0     for  ARG .LE. XLOW.


		Intrinsic functions required are:

		ABS, AINT, EXP


		Author: W. J. Cody
		Mathematics and Computer Science Division
		Argonne National Laboratory
		Argonne, IL 60439

		Latest modification: March 15, 1992

		------------------------------------------------------------------
		*/
	{
		static double a[5] = {
			2.2352520354606839287e00,1.6102823106855587881e02,1.0676894854603709582e03,
			1.8154981253343561249e04,6.5682337918207449113e-2
		};
		static double b[4] = {
			4.7202581904688241870e01,9.7609855173777669322e02,1.0260932208618978205e04,
			4.5507789335026729956e04
		};
		static double c[9] = {
			3.9894151208813466764e-1,8.8831497943883759412e00,9.3506656132177855979e01,
			5.9727027639480026226e02,2.4945375852903726711e03,6.8481904505362823326e03,
			1.1602651437647350124e04,9.8427148383839780218e03,1.0765576773720192317e-8
		};
		static double d[8] = {
			2.2266688044328115691e01,2.3538790178262499861e02,1.5193775994075548050e03,
			6.4855582982667607550e03,1.8615571640885098091e04,3.4900952721145977266e04,
			3.8912003286093271411e04,1.9685429676859990727e04
		};
		static double half = 0.5e0;
		static double p[6] = {
			2.1589853405795699e-1,1.274011611602473639e-1,2.2235277870649807e-2,
			1.421619193227893466e-3,2.9112874951168792e-5,2.307344176494017303e-2
		};
		static double one = 1.0e0;
		static double q[5] = {
			1.28426009614491121e00,4.68238212480865118e-1,6.59881378689285515e-2,
			3.78239633202758244e-3,7.29751555083966205e-5
		};
		static double sixten = 1.60e0;
		static double sqrpi = 3.9894228040143267794e-1;
		static double thrsh = 0.66291e0;
		static double root32 = 5.656854248e0;
		static double zero = 0.0e0;
		static int K1 = 1;
		static int K2 = 2;
		static int i;
		static double del,eps,temp,x,xden,xnum,y,xsq,min;
		/*
		------------------------------------------------------------------
		Machine dependent constants
		------------------------------------------------------------------
		*/
		eps = x_spmpar(K1)*0.5e0;
		min = x_spmpar(K2);
		x = *arg;
		y = fabs(x);
		if(y <= thrsh) {
			/*
			------------------------------------------------------------------
			Evaluate  anorm  for  |X| <= 0.66291
			------------------------------------------------------------------
			*/
			xsq = zero;
			if(y > eps) xsq = x*x;
			xnum = a[4]*xsq;
			xden = xsq;
			for(i=0; i<3; i++) {
				xnum = (xnum+a[i])*xsq;
				xden = (xden+b[i])*xsq;
			}
			*result = x*(xnum+a[3])/(xden+b[3]);
			temp = *result;
			*result = half+temp;
			*ccum = half-temp;
		}
		/*
		------------------------------------------------------------------
		Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
		------------------------------------------------------------------
		*/
		else if(y <= root32) {
			xnum = c[8]*y;
			xden = y;
			for(i=0; i<7; i++) {
				xnum = (xnum+c[i])*y;
				xden = (xden+d[i])*y;
			}
			*result = (xnum+c[7])/(xden+d[7]);
			xsq = x_fifdint(y*sixten)/sixten;
			del = (y-xsq)*(y+xsq);
			*result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
			*ccum = one-*result;
			if(x > zero) {
				temp = *result;
				*result = *ccum;
				*ccum = temp;
			}
		}
		/*
		------------------------------------------------------------------
		Evaluate  anorm  for |X| > sqrt(32)
		------------------------------------------------------------------
		*/
		else  {
			*result = zero;
			xsq = one/(x*x);
			xnum = p[5]*xsq;
			xden = xsq;
			for(i=0; i<4; i++) {
				xnum = (xnum+p[i])*xsq;
				xden = (xden+q[i])*xsq;
			}
			*result = xsq*(xnum+p[4])/(xden+q[4]);
			*result = (sqrpi-*result)/y;
			xsq = x_fifdint(x*sixten)/sixten;
			del = (x-xsq)*(x+xsq);
			*result = exp(-(xsq*xsq*half))*exp(-(del*half))**result;
			*ccum = one-*result;
			if(x > zero) {
				temp = *result;
				*result = *ccum;
				*ccum = temp;
			}
		}
		if(*result < min) *result = 0.0e0;
		/*
		------------------------------------------------------------------
		Fix up for negative argument, erf, etc.
		------------------------------------------------------------------
		----------Last card of ANORM ----------
		*/
		if(*ccum < min) *ccum = 0.0e0;
	}

	void x_cumpoi(double s,double xlam,double *cum,double *ccum)
		/*
		**********************************************************************

		void cumpoi(double *s,double *xlam,double *cum,double *ccum)
		CUMulative POIsson distribution


		Function


		Returns the  probability  of  S   or  fewer events in  a   Poisson
		distribution with mean XLAM.


		Arguments


		S --> Upper limit of cumulation of the Poisson.
		S is DOUBLE PRECISION

		XLAM --> Mean of the Poisson distribution.
		XLAM is DOUBLE PRECIS

		CUM <-- Cumulative poisson distribution.
		CUM is DOUBLE PRECISION

		CCUM <-- Compliment of Cumulative poisson distribution.
		CCUM is DOUBLE PRECIS


		Method


		Uses formula  26.4.21   of   Abramowitz and  Stegun,  Handbook  of
		Mathematical   Functions  to reduce   the   cumulative Poisson  to
		the cumulative chi-square distribution.

		**********************************************************************
		*/
	{
		double chi,df;
		/*
		..
		.. Executable Statements ..
		*/
		df = 2.0*(s+1.0);
		chi = 2.0*xlam;
		x_cumchi(chi,df,ccum,cum);
	}
	void x_cumt(double *t,double *df,double *cum,double *ccum)
		/*
		**********************************************************************

		void cumt(double *t,double *df,double *cum,double *ccum)
		CUMulative T-distribution


		Function


		Computes the integral from -infinity to T of the t-density.


		Arguments


		T --> Upper limit of integration of the t-density.
		T is DOUBLE PRECISION

		DF --> Degrees of freedom of the t-distribution.
		DF is DOUBLE PRECISIO

		CUM <-- Cumulative t-distribution.
		CCUM is DOUBLE PRECIS

		CCUM <-- Compliment of Cumulative t-distribution.
		CCUM is DOUBLE PRECIS


		Method


		Formula 26.5.27   of     Abramowitz  and   Stegun,    Handbook  of
		Mathematical Functions  is   used   to  reduce the  t-distribution
		to an incomplete beta.

		**********************************************************************
		*/
	{
		static double K2 = 0.5e0;
		double xx,a,oma,tt,yy,dfptt;
		/*
		..
		.. Executable Statements ..
		*/
		tt = *t**t;
		dfptt = *df+tt;
		xx = *df/dfptt;
		yy = tt/dfptt;
		x_cumbet(xx,yy,0.5**df,0.5,&a,&oma);
		if(!(*t <= 0.0e0)) goto S10;
		*cum = 0.5e0*a;
		*ccum = oma+*cum;
		goto S20;
	S10:
		*ccum = 0.5e0*a;
		*cum = oma+*ccum;
	S20:
		return;
	}
	double x_dblSterlingRmder4CompBetaFun(double a,double b)
		/*
		**********************************************************************

		double dbetrm(double *a,double *b)
		Double Precision Sterling Remainder for Complete
		Beta Function


		Function


		Log(Beta(A,B)) = Lgamma(A) + Lgamma(B) - Lgamma(A+B)
		where Lgamma is the log of the (complete) gamma function

		Let ZZ be approximation obtained if each log gamma is approximated
		by Sterling's formula, i.e.,
		Sterling(Z) = LOG( SQRT( 2*PI ) ) + ( Z-0.5 ) * LOG( Z ) - Z

		Returns Log(Beta(A,B)) - ZZ


		Arguments


		A --> One argument of the Beta
		DOUBLE PRECISION A

		B --> The other argument of the Beta
		DOUBLE PRECISION B

		**********************************************************************
		*/
	{
		double _dbetrm;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		Try to sum from smallest to largest
		*/
		_dbetrm = -x_dstrem(a+b);
		_dbetrm += x_dstrem(x_fifdmax1(a,b));
		_dbetrm += x_dstrem(x_fifdmin1(a,b));
		return _dbetrm;
	}
	double x_dblEvalPolyX(double a[],int n,double x)
		/*
		**********************************************************************

		double devlpl(double a[],int *n,double *x)
		Double precision EVALuate a PoLynomial at X


		Function


		returns
		A(1) + A(2)*X + ... + A(N)*X**(N-1)


		Arguments


		A --> Array of coefficients of the polynomial.
		A is DOUBLE PRECISION(N)

		N --> Length of A, also degree of polynomial - 1.
		N is INTEGER

		X --> Point at which the polynomial is to be evaluated.
		X is DOUBLE PRECISION

		**********************************************************************
		*/
	{
		/*
		..
		.. Executable Statements ..
		*/
		double term = a[n-1];
		for(int i=n-2 ; i>=0 ; i--) term = a[i] + term*x;
		return term;
	}
	double x_dblEvalExpXminus1(double x)
		/*
		**********************************************************************

		double dexpm1(double *x)
		Evaluation of the function EXP(X) - 1


		Arguments


		X --> Argument at which exp(x)-1 desired
		DOUBLE PRECISION X


		Method


		Renaming of function rexp from code of:

		DiDinato, A. R. and Morris,  A.   H.  Algorithm 708: Significant
		Digit Computation of the Incomplete  Beta  Function Ratios.  ACM
		Trans. Math.  Softw. 18 (1993), 360-373.

		**********************************************************************
		*/
	{
		static double p1 = .914041914819518e-09;
		static double p2 = .238082361044469e-01;
		static double q1 = -.499999999085958e+00;
		static double q2 = .107141568980644e+00;
		static double q3 = -.119041179760821e-01;
		static double q4 = .595130811860248e-03;
		double w;
		/*
		..
		.. Executable Statements ..
		*/
		if(fabs(x) > 0.15) goto S10;
		return x*(((p2*x+p1)*x+1.0)/((((q4*x+q3)*x+q2)*x+q1)*x+1.0));
	S10:
		w = exp(x);
		if(x > 0.0) goto S20;
		return w-1.0;
	S20:
		return w*(0.5+(0.5-1.0/w));
	}
	double x_dblNormalDistInv(double p,double q)
		/*
		**********************************************************************

		double dinvnr(double *p,double *q)
		Double precision NoRmal distribution INVerse


		Function


		Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
		infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P


		Arguments


		P --> The probability whose normal deviate is sought.
		P is DOUBLE PRECISION

		Q --> 1-P
		P is DOUBLE PRECISION


		Method


		The  rational   function   on  page 95    of Kennedy  and  Gentle,
		Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
		value for the Newton method of finding roots.


		Note


		If P or Q .lt. machine EPS returns +/- DINVNR(EPS)

		**********************************************************************
		*/
	{
#define maxit 100
#define eps (1.0e-13)
#define r2pi 0.3989422804014326e0
#define nhalf (-0.5e0)
#define dennor(x) (r2pi*exp(nhalf*(x)*(x)))
		double _dinvnr,strtx,xcur,cum,ccum,pp,dx;
		int i;
		unsigned long qporq;
		/*
		..
		.. Executable Statements ..
		*/
		/*
		FIND MINIMUM OF P AND Q
		*/
		qporq = p <= q;
		if(!qporq) goto S10;
		pp = p;
		goto S20;
	S10:
		pp = q;
	S20:
		/*
		INITIALIZATION STEP
		*/
		strtx = x_stvaln(&pp);
		xcur = strtx;
		/*
		NEWTON INTERATIONS
		*/
		for(i=1; i<=maxit; i++) {
			x_cumnor(&xcur,&cum,&ccum);
			dx = (cum-pp)/dennor(xcur);
			xcur -= dx;
			if(fabs(dx/xcur) < eps) goto S40;
		}
		_dinvnr = strtx;
		/*
		IF WE GET HERE, NEWTON HAS FAILED
		*/
		if(!qporq) _dinvnr = -_dinvnr;
		return _dinvnr;
	S40:
		/*
		IF WE GET HERE, NEWTON HAS SUCCEDED
		*/
		_dinvnr = xcur;
		if(!qporq) _dinvnr = -_dinvnr;
		return _dinvnr;
#undef maxit
#undef eps
#undef r2pi
#undef nhalf
#undef dennor
	}
	/* DEFINE DINVR */
	static void E0000(int IENTRY,int *status,double *x,double fx,
		unsigned long *qleft,unsigned long *qhi,double zabsst=W0,
		double zabsto=W0,double zbig=W0,double zrelst=W0,
		double zrelto=W0,double zsmall=W0,double zstpmu=W0)
	{
#define qxmon(zx,zy,zz) (int)((zx) <= (zy) && (zy) <= (zz))
		static double absstp,abstol,big,fbig,fsmall,relstp,reltol,sm,step,stpmul,xhi,
			xlb,xlo,xsave,xub,yy;
		static int i99999;
		static unsigned long qbdd,qcond,qdum1,qdum2,qincr,qlim,qup;
		switch(IENTRY){case 0: goto DINVR; case 1: goto DSTINV;}
	DINVR:
		if(*status > 0) goto S310;
		qcond = !qxmon(sm,*x,big);
		if(qcond) halt(" SMALL, X, BIG not monotone in INVR");
		xsave = *x;
		/*
		See that SMALL and BIG bound the zero and set QINCR
		*/
		*x = sm;
		/*
		GET-FUNCTION-VALUE
		*/
		i99999 = 1;
		goto S300;
	S10:
		fsmall = fx;
		*x = big;
		/*
		GET-FUNCTION-VALUE
		*/
		i99999 = 2;
		goto S300;
	S20:
		fbig = fx;
		qincr = fbig > fsmall;
		if(!qincr) goto S50;
		if(fsmall <= 0.0e0) goto S30;
		*status = -1;
		*qleft = *qhi = 1;
		return;
	S30:
		if(fbig >= 0.0e0) goto S40;
		*status = -1;
		*qleft = *qhi = 0;
		return;
	S40:
		goto S80;
	S50:
		if(fsmall >= 0.0e0) goto S60;
		*status = -1;
		*qleft = 1;
		*qhi = 0;
		return;
	S60:
		if(fbig <= 0.0e0) goto S70;
		*status = -1;
		*qleft = 0;
		*qhi = 1;
		return;
	S80:
	S70:
		*x = xsave;
		step = x_fifdmax1(absstp,relstp*fabs(*x));
		/*
		YY = F(X) - Y
		GET-FUNCTION-VALUE
		*/
		i99999 = 3;
		goto S300;
	S90:
		yy = fx;
		if(!(yy == 0.0e0)) goto S100;
		*status = 0;
		return;
	S100:
		qup = (qincr && yy < 0.0e0) || (!qincr && yy > 0.0e0);
		/*
		++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		HANDLE CASE IN WHICH WE MUST STEP HIGHER
		++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		*/
		if(!qup) goto S170;
		xlb = xsave;
		xub = x_fifdmin1(xlb+step,big);
		goto S120;
	S110:
		if(qcond) goto S150;
	S120:
		/*
		YY = F(XUB) - Y
		*/
		*x = xub;
		/*
		GET-FUNCTION-VALUE
		*/
		i99999 = 4;
		goto S300;
	S130:
		yy = fx;
		qbdd = (qincr && yy >= 0.0e0) || (!qincr && yy <= 0.0e0);
		qlim = xub >= big;
		qcond = qbdd || qlim;
		if(qcond) goto S140;
		step = stpmul*step;
		xlb = xub;
		xub = x_fifdmin1(xlb+step,big);
	S140:
		goto S110;
	S150:
		if(!(qlim && !qbdd)) goto S160;
		*status = -1;
		*qleft = 0;
		*qhi = !qincr;
		*x = big;
		return;
	S160:
		goto S240;
	S170:
		/*
		++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		HANDLE CASE IN WHICH WE MUST STEP LOWER
		++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		*/
		xub = xsave;
		xlb = x_fifdmax1(xub-step,sm);
		goto S190;
	S180:
		if(qcond) goto S220;
	S190:
		/*
		YY = F(XLB) - Y
		*/
		*x = xlb;
		/*
		GET-FUNCTION-VALUE
		*/
		i99999 = 5;
		goto S300;
	S200:
		yy = fx;
		qbdd = (qincr && yy <= 0.0e0) || (!qincr && yy >= 0.0e0);
		qlim = xlb <= sm;
		qcond = qbdd || qlim;
		if(qcond) goto S210;
		step = stpmul*step;
		xub = xlb;
		xlb = x_fifdmax1(xub-step,sm);
	S210:
		goto S180;
	S220:
		if(!(qlim && !qbdd)) goto S230;
		*status = -1;
		*qleft = 1;
		*qhi = qincr;
		*x = sm;
		return;
	S240:
	S230:
		dstzr(xlb,xub,abstol,reltol);
		/*
		++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
		++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
		*/
		*status = 0;
		goto S260;
	S250:
		if(!(*status == 1)) goto S290;
	S260:
		dzror(status,x,fx,&xlo,&xhi,&qdum1,&qdum2);
		if(!(*status == 1)) goto S280;
		/*
		GET-FUNCTION-VALUE
		*/
		i99999 = 6;
		goto S300;
	S280:
	S270:
		goto S250;
	S290:
		*x = xlo;
		*status = 0;
		return;
	DSTINV:
		sm = zsmall;
		big = zbig;
		absstp = zabsst;
		relstp = zrelst;
		stpmul = zstpmu;
		abstol = zabsto;
		reltol = zrelto;
		return;
	S300:
		/*
		TO GET-FUNCTION-VALUE
		*/
		*status = 1;
		return;
	S310:
		switch((int)i99999){case 1: goto S10;case 2: goto S20;case 3: goto S90;case 
			4: goto S130;case 5: goto S200;case 6: goto S270;default: break;}
#undef qxmon
	}
	static void E0000_v2(int IENTRY, int *status, double *x, double fx,
		unsigned long *qleft,
		unsigned long *qhi,
		/* PARAM */ double zabsst=W0,
		/* PARAM */ double zabsto=W0,
		/* PARAM */ double zbig=W0,
		/* PARAM */ double zrelst=W0,
		/* PARAM */ double zrelto=W0,
		/* PARAM */ double zsmall=W0,
		/* PARAM */ double zstpmu=W0)
	{
		static double sm, big, absstp, relstp, stpmul, abstol, reltol;
		if (IENTRY == 1) {
			sm = zsmall;
			big = zbig;
			absstp = zabsst;
			relstp = zrelst;
			stpmul = zstpmu;
			abstol = zabsto;
			reltol = zrelto;
			return;
		}
		else if (IENTRY != 0) {
			halt("E0000 error");
			return;
		}

#define qxmon(zx,zy,zz) (int)((zx) <= (zy) && (zy) <= (zz))
		static double			fsmall, step, xhi, xlb, xlo, xsave, xub, yy;
		static int				i99999;
		static unsigned long	qincr;

		unsigned long qbdd, qcond, qlim;

		// if *status <= 0
		if (*status <= 0) {
			qcond = !qxmon(sm, *x, big);
			if (qcond) halt(" SMALL, X, BIG not monotone in INVR");
			xsave = *x;
			// See that SMALL and BIG bound the zero and set QINCR
			*x = sm;
			// GET-FUNCTION-VALUE
			i99999 = 1;
			*status = 1;
			return;
		}

		switch ((int)i99999) {
		case 1:
			fsmall = fx;
			*x = big;
			// GET-FUNCTION-VALUE
			i99999 = 2;
			*status = 1;
			return;
		case 2: {
			double fbig = fx;
			qincr = fbig > fsmall;
			if (!qincr) {
				if (fsmall >= 0.0e0) {
					if (fbig > 0.0e0) {
						*status = -1;
						*qleft = 0;
						*qhi = 1;
						return;
					}
					// goto S80
				} else {
					*status = -1;
					*qleft = 1;
					*qhi = 0;
					return;
				}
			} else if (fsmall <= 0.0e0) {
				if (fbig < 0.0e0) {
					*status = -1;
					*qleft = *qhi = 0;
					return;
				}
				// goto S80
			} else {
				*status = -1;
				*qleft = *qhi = 1;
				return;
			}
			// S80
			*x = xsave;
			step = x_fifdmax1(absstp, relstp*fabs(*x));
			// YY = F(X) - Y ==> GET-FUNCTION-VALUE
			i99999 = 3;
			*status = 1;
		}
				return;
		case 3:
			yy = fx;
			if (!(yy == 0.0e0)) {
				/* HANDLE CASE IN WHICH WE MUST STEP HIGHER */
				if ((qincr && yy < 0.0e0) || (!qincr && yy > 0.0e0)) {
					xlb = xsave;
					xub = x_fifdmin1(xlb + step, big);
					// YY = F(XUB) - Y
					*x = xub;
					// GET-FUNCTION-VALUE
					i99999 = 4;
					*status = 1;
					return;
				}
				// HANDLE CASE IN WHICH WE MUST STEP LOWER
				xub = xsave;
				xlb = x_fifdmax1(xub - step, sm);
				// YY = F(XLB) - Y
				*x = xlb;
				// GET-FUNCTION-VALUE
				i99999 = 5;
				*status = 1;
				return;
			}
			*status = 0;
			return;
		case 4:
			// S130
			do {
				yy = fx;
				qbdd = (qincr && yy >= 0.0e0) || (!qincr && yy <= 0.0e0);
				qlim = xub >= big;
				qcond = qbdd || qlim;
				if (qcond) break;
				step = stpmul*step;
				xlb = xub;
				xub = x_fifdmin1(xlb + step, big);
			} while (!qcond);

			if (!(qlim && !qbdd)) goto S240;

			*status = -1;
			*qleft = 0;
			*qhi = !qincr;
			*x = big;
			return;
		case 5:
			// S200
			do {
				yy = fx;
				qbdd = (qincr && yy <= 0.0e0) || (!qincr && yy >= 0.0e0);
				qlim = xlb <= sm;
				qcond = qbdd || qlim;
				if (!qcond) {
					step = stpmul*step;
					xub = xlb;
					xlb = x_fifdmax1(xub - step, sm);
				}
			} while (!qcond);

			// S220
			if (qlim && !qbdd) {
				*status = -1;
				*qleft = 1;
				*qhi = qincr;
				*x = sm;
				return;
			}
		S240:
			dstzr(xlb, xub, abstol, reltol);
			// IF WE REACH HERE, XLB AND XUB BOUND THE ZERO OF F.
			*status = 0;
		S260:
			dzror(status, x, fx, &xlo, &xhi);
			if (*status == 1) {
				// GET-FUNCTION-VALUE
				i99999 = 6;
				*status = 1;
				return;
			}
		case 6:
			if (!(*status == 1)) {
				*x = xlo;
				*status = 0;
				return;
			}
			goto S260;
		default: break;
		}
		// CANNOT BE REACHED
		halt("E0000 cannot be reached to here");
#undef qxmon
	}
	void dinvr(int *status,double *x,double fx,
		unsigned long *qleft,unsigned long *qhi)
		/*
		**********************************************************************

		void dinvr(int *status,double *x,double *fx,
		unsigned long *qleft,unsigned long *qhi)

		Double precision
		bounds the zero of the function and invokes zror
		Reverse Communication


		Function


		Bounds the    function  and  invokes  ZROR   to perform the   zero
		finding.  STINVR  must  have   been  called  before this   routine
		in order to set its parameters.


		Arguments


		STATUS <--> At the beginning of a zero finding problem, STATUS
		should be set to 0 and INVR invoked.  (The value
		of parameters other than X will be ignored on this cal

		When INVR needs the function evaluated, it will set
		STATUS to 1 and return.  The value of the function
		should be set in FX and INVR again called without
		changing any of its other parameters.

		When INVR has finished without error, it will return
		with STATUS 0.  In that case X is approximately a root
		of F(X).

		If INVR cannot bound the function, it returns status
		-1 and sets QLEFT and QHI.
		INTEGER STATUS

		X <-- The value of X at which F(X) is to be evaluated.
		DOUBLE PRECISION X

		FX --> The value of F(X) calculated when INVR returns with
		STATUS = 1.
		DOUBLE PRECISION FX

		QLEFT <-- Defined only if QMFINV returns .FALSE.  In that
		case it is .TRUE. If the stepping search terminated
		unsucessfully at SMALL.  If it is .FALSE. the search
		terminated unsucessfully at BIG.
		QLEFT is LOGICAL

		QHI <-- Defined only if QMFINV returns .FALSE.  In that
		case it is .TRUE. if F(X) .GT. Y at the termination
		of the search and .FALSE. if F(X) .LT. Y at the
		termination of the search.
		QHI is LOGICAL

		**********************************************************************
		*/
	{
		//	E0000_v2(0, status, x, fx, qleft, qhi);
		E0000(0, status, x, fx, qleft, qhi);
	}
	void dstinv(double zsmall,double zbig,double zabsst,
		double zrelst,double zstpmu,double zabsto,
		double zrelto)
		/*
		**********************************************************************
		void dstinv(double *zsmall,double *zbig,double *zabsst,
		double *zrelst,double *zstpmu,double *zabsto,
		double *zrelto)

		Double Precision - SeT INverse finder - Reverse Communication
		Function
		Concise Description - Given a monotone function F finds X
		such that F(X) = Y.  Uses Reverse communication -- see invr.
		This routine sets quantities needed by INVR.
		More Precise Description of INVR -
		F must be a monotone function, the results of QMFINV are
		otherwise undefined.  QINCR must be .TRUE. if F is non-
		decreasing and .FALSE. if F is non-increasing.
		QMFINV will return .TRUE. if and only if F(SMALL) and
		F(BIG) bracket Y, i. e.,
		QINCR is .TRUE. and F(SMALL).LE.Y.LE.F(BIG) or
		QINCR is .FALSE. and F(BIG).LE.Y.LE.F(SMALL)
		if QMFINV returns .TRUE., then the X returned satisfies
		the following condition.  let
		TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
		then if QINCR is .TRUE.,
		F(X-TOL(X)) .LE. Y .LE. F(X+TOL(X))
		and if QINCR is .FALSE.
		F(X-TOL(X)) .GE. Y .GE. F(X+TOL(X))
		Arguments
		SMALL --> The left endpoint of the interval to be
		searched for a solution.
		SMALL is DOUBLE PRECISION
		BIG --> The right endpoint of the interval to be
		searched for a solution.
		BIG is DOUBLE PRECISION
		ABSSTP, RELSTP --> The initial step size in the search
		is MAX(ABSSTP,RELSTP*ABS(X)). See algorithm.
		ABSSTP is DOUBLE PRECISION
		RELSTP is DOUBLE PRECISION
		STPMUL --> When a step doesn't bound the zero, the step
		size is multiplied by STPMUL and another step
		taken.  A popular value is 2.0
		DOUBLE PRECISION STPMUL
		ABSTOL, RELTOL --> Two numbers that determine the accuracy
		of the solution.  See function for a precise definition.
		ABSTOL is DOUBLE PRECISION
		RELTOL is DOUBLE PRECISION
		Method
		Compares F(X) with Y for the input value of X then uses QINCR
		to determine whether to step left or right to bound the
		desired x.  the initial step size is
		MAX(ABSSTP,RELSTP*ABS(S)) for the input value of X.
		Iteratively steps right or left until it bounds X.
		At each step which doesn't bound X, the step size is doubled.
		The routine is careful never to step beyond SMALL or BIG.  If
		it hasn't bounded X at SMALL or BIG, QMFINV returns .FALSE.
		after setting QLEFT and QHI.
		If X is successfully bounded then Algorithm R of the paper
		'Two Efficient Algorithms with Guaranteed Convergence for
		Finding a Zero of a Function' by J. C. P. Bus and
		T. J. Dekker in ACM Transactions on Mathematical
		Software, Volume 1, No. 4 page 330 (DEC. '75) is employed
		to find the zero of the function F(X)-Y. This is routine
		QRZERO.
		**********************************************************************
		*/
	{
		//	E0000_v2(1, NULL, NULL, NULL, NULL, NULL, zabsst, zabsto, zbig, zrelst, zrelto, zsmall,
		//		zstpmu);
		E0000(1, NULL, NULL, W0, NULL, NULL, zabsst, zabsto, zbig, zrelst, zrelto, zsmall,
			zstpmu);
	}
	/* DEFINE DZROR */
	static void E0001(int IENTRY,int *status,double *x,double fx,
		double *xlo,double *xhi,unsigned long *qleft,
		unsigned long *qhi,double zabstl=W0,double zreltl=W0,
		double zxhi=W0,double zxlo=W0)
	{
		static double a,abstol,b,c,d,fa,fb,fc,fd,fda,fdb,m,mb,p,q,reltol,tol,w,xxhi,xxlo;
		static 	int ext,i99999;
		static unsigned long first,qrzero;
		switch(IENTRY){case 0: goto DZROR; case 1: goto DSTZR;}
	DZROR:
		if(*status > 0) goto S280;
		*xlo = xxlo;
		*xhi = xxhi;
		b = *x = *xlo;
		/*
		GET-FUNCTION-VALUE
		*/
		i99999 = 1;
		*status = 1;
		return;
	S10:
		fb = fx;
		*xlo = *xhi;
		a = *x = *xlo;
		/*
		GET-FUNCTION-VALUE
		*/
		i99999 = 2;
		*status = 1;
		return;
	S20:
		/*
		Check that F(ZXLO) < 0 < F(ZXHI)  or
		F(ZXLO) > 0 > F(ZXHI)
		*/
		if(!(fb < 0.0)) goto S40;
		if(!(fx < 0.0)) goto S30;
		*status = -1;
		*qleft = fx < fb;
		*qhi = 0;
		return;
	S40:
	S30:
		if(!(fb > 0.0)) goto S60;
		if(!(fx > 0.0)) goto S50;
		*status = -1;
		*qleft = fx > fb;
		*qhi = 1;
		return;
	S60:
	S50:
		fa = fx;
		first = 1;
	S70:
		c = a;
		fc = fa;
		ext = 0;
	S80:
		if(!(fabs(fc) < fabs(fb))) goto S100;
		if(!(c != a)) goto S90;
		d = a;
		fd = fa;
	S90:
		a = b;
		fa = fb;
		*xlo = c;
		b = *xlo;
		fb = fc;
		c = a;
		fc = fa;
	S100:
		tol = 0.5 * x_fifdmax1(abstol,reltol*fabs(*xlo));
		m = (c+b)*.5;
		mb = m-b;
		if(!(fabs(mb) > tol)) goto S240;
		if(!(ext > 3)) goto S110;
		w = mb;
		goto S190;
	S110:
		tol = x_fifdsign(tol,mb);
		p = (b-a)*fb;
		if(!first) goto S120;
		q = fa-fb;
		first = 0;
		goto S130;
	S120:
		fdb = (fd-fb)/(d-b);
		fda = (fd-fa)/(d-a);
		p = fda*p;
		q = fdb*fa-fda*fb;
	S130:
		if(!(p < 0.0)) goto S140;
		p = -p;
		q = -q;
	S140:
		if(ext == 3) p *= 2.0;
		if(!(p*1.0 == 0.0 || p <= q*tol)) goto S150;
		w = tol;
		goto S180;
	S150:
		if(!(p < mb*q)) goto S160;
		w = p/q;
		goto S170;
	S160:
		w = mb;
	S190:
	S180:
	S170:
		d = a;
		fd = fa;
		a = b;
		fa = fb;
		b += w;
		*xlo = b;
		*x = *xlo;
		/*
		GET-FUNCTION-VALUE
		*/
		i99999 = 3;
		*status = 1;
		return;
	S200:
		fb = fx;
		if(!(fc*fb >= 0.0)) goto S210;
		goto S70;
	S210:
		if(!(w == mb)) goto S220;
		ext = 0;
		goto S230;
	S220:
		ext += 1;
	S230:
		goto S80;
	S240:
		*xhi = c;
		qrzero = (fc >= 0.0 && fb <= 0.0) || (fc < 0.0 && fb >= 0.0);
		if(!qrzero) goto S250;
		*status = 0;
		goto S260;
	S250:
		*status = -1;
	S260:
		return;
	DSTZR:
		xxlo = zxlo;
		xxhi = zxhi;
		abstol = zabstl;
		reltol = zreltl;
		return;
	S280:
		switch((int)i99999){case 1: goto S10;case 2: goto S20;case 3: goto S200;
		default: break;}
	}
	/* DEFINE DZROR */
	static void E0001_v2(int IENTRY,int *status,double *x,double fx,
		double *xlo,double *xhi,
		/* RETRN */ unsigned long *qleft=NULL,
		/* RETRN */ unsigned long *qhi=NULL,

		/* PARAM */ double zabstl=W0,
		/* PARAM */ double zreltl=W0,
		/* PARAM */ double zxhi=W0,
		/* PARAM */ double zxlo=W0)
	{
		static double			xxlo, xxhi, abstol, reltol;
		static int				ext, i99999;
		static unsigned long	first;
		static double			a, b, c, d, w, fa, fb, fc, fd, mb;

		if (IENTRY == 1) {
			xxlo = zxlo;
			xxhi = zxhi;
			abstol = zabstl;
			reltol = zreltl;
			return;
		} else if (IENTRY != 0) {
			halt("E0001 error");
			return;
		}

		// if *status <= 0
		if(*status <= 0) {
			*xlo = xxlo;
			*xhi = xxhi;
			b = *x = *xlo;
			// GET-FUNCTION-VALUE
			i99999 = 1;
			*status = 1;
			return;
		}

		switch ((int)i99999) {
		case 1:
			fb = fx;
			*xlo = *xhi;
			a = *x = *xlo;
			// GET-FUNCTION-VALUE
			i99999 = 2;
			*status = 1;
			return;
		case 2:
			/*
			Check that F(ZXLO) < 0 < F(ZXHI)  or
			F(ZXLO) > 0 > F(ZXHI)
			*/
			if((fb < 0.0) && (fx < 0.0)) {
				*status = -1;
				if (qleft) *qleft = fx < fb;
				if (qhi) *qhi = 0;
				return;
			}
			if((fb > 0.0) && (fx > 0.0)) {
				*status = -1;
				if (qleft) *qleft = fx > fb;
				if (qhi) *qhi = 1;
				return;
			}
			fa = fx;
			first = 1;
			c = a;
			fc = fa;
			ext = 0;
			goto S80;
		case 3:
			fb = fx;
			if(!(fc*fb >= 0.0)) {
				if (w == mb)
					ext = 0;
				else
					ext++;
			} else {
				c = a;
				fc = fa;
				ext = 0;
			}
			goto S80;
		default: break;
		}
		return;
		////////////////////////
	S80:
		if(fabs(fc) < fabs(fb)) {
			if(c != a) {
				d = a;
				fd = fa;
			}
			a = b;
			fa = fb;
			*xlo = c;
			b = *xlo;
			fb = fc;
			c = a;
			fc = fa;
		}
		double tol = 0.5 * x_fifdmax1(abstol,reltol*fabs(*xlo));
		double m = (c+b)*.5;
		mb = m-b;
		if (fabs(mb) <= tol) {
			*xhi = c;
			*status = (fc >= 0.0 && fb <= 0.0) || (fc < 0.0 && fb >= 0.0) ? 0 : -1;
			return;
		}

		if (ext > 3) w = mb;
		else {
			tol = x_fifdsign(tol,mb);
			double p = (b-a)*fb;
			double q;
			if (!first) {
				double fdb = (fd-fb)/(d-b);
				double fda = (fd-fa)/(d-a);
				p = fda*p;
				q = fdb*fa-fda*fb;
			} else {
				q = fa-fb;
				first = 0;
			}
			if (p < 0.0) {
				p = -p;
				q = -q;
			}
			if (ext == 3) p *= 2.0;

			if (p*1.0 == 0.0 || p <= q*tol)
				w = tol;
			else
				w = p < mb*q ? p/q : mb;
		}

		d = a;
		fd = fa;
		a = b;
		fa = fb;
		b += w;
		*xlo = b;
		*x = *xlo;
		// GET-FUNCTION-VALUE
		i99999 = 3;
		*status = 1;
		return;
	}

	void dzror(int *status,double *x,double fx,double *xlo,
		double *xhi,unsigned long *qleft/*=NULL*/,unsigned long *qhi/*=NULL*/)
		/*
		**********************************************************************

		void dzror(int *status,double *x,double fx,double *xlo,
		double *xhi,unsigned long *qleft,unsigned long *qhi)

		Double precision ZeRo of a function -- Reverse Communication


		Function


		Performs the zero finding.  STZROR must have been called before
		this routine in order to set its parameters.


		Arguments


		STATUS <--> At the beginning of a zero finding problem, STATUS
		should be set to 0 and ZROR invoked.  (The value
		of other parameters will be ignored on this call.)

		When ZROR needs the function evaluated, it will set
		STATUS to 1 and return.  The value of the function
		should be set in FX and ZROR again called without
		changing any of its other parameters.

		When ZROR has finished without error, it will return
		with STATUS 0.  In that case (XLO,XHI) bound the answe

		If ZROR finds an error (which implies that F(XLO)-Y an
		F(XHI)-Y have the same sign, it returns STATUS -1.  In
		this case, XLO and XHI are undefined.
		INTEGER STATUS

		X <-- The value of X at which F(X) is to be evaluated.
		DOUBLE PRECISION X

		FX --> The value of F(X) calculated when ZROR returns with
		STATUS = 1.
		DOUBLE PRECISION FX

		XLO <-- When ZROR returns with STATUS = 0, XLO bounds the
		inverval in X containing the solution below.
		DOUBLE PRECISION XLO

		XHI <-- When ZROR returns with STATUS = 0, XHI bounds the
		inverval in X containing the solution above.
		DOUBLE PRECISION XHI

		QLEFT <-- .TRUE. if the stepping search terminated unsucessfully
		at XLO.  If it is .FALSE. the search terminated
		unsucessfully at XHI.
		QLEFT is LOGICAL

		QHI <-- .TRUE. if F(X) .GT. Y at the termination of the
		search and .FALSE. if F(X) .LT. Y at the
		termination of the search.
		QHI is LOGICAL

		**********************************************************************
		*/
	{
		E0001(0, status, x, fx, xlo, xhi, qleft, qhi);
		//	E0001_v2(0, status, x, fx, xlo, xhi, qleft, qhi);
	}
	void dstzr(double zxlo,double zxhi,double zabstl,double zreltl)
		/*
		**********************************************************************
		void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl)
		Double precision SeT ZeRo finder - Reverse communication version
		Function
		Sets quantities needed by ZROR.  The function of ZROR
		and the quantities set is given here.
		Concise Description - Given a function F
		find XLO such that F(XLO) = 0.
		More Precise Description -
		Input condition. F is a double precision function of a single
		double precision argument and XLO and XHI are such that
		F(XLO)*F(XHI)  .LE.  0.0
		If the input condition is met, QRZERO returns .TRUE.
		and output values of XLO and XHI satisfy the following
		F(XLO)*F(XHI)  .LE. 0.
		ABS(F(XLO)  .LE. ABS(F(XHI)
		ABS(XLO-XHI)  .LE. TOL(X)
		where
		TOL(X) = MAX(ABSTOL,RELTOL*ABS(X))
		If this algorithm does not find XLO and XHI satisfying
		these conditions then QRZERO returns .FALSE.  This
		implies that the input condition was not met.
		Arguments
		XLO --> The left endpoint of the interval to be
		searched for a solution.
		XLO is DOUBLE PRECISION
		XHI --> The right endpoint of the interval to be
		for a solution.
		XHI is DOUBLE PRECISION
		ABSTOL, RELTOL --> Two numbers that determine the accuracy
		of the solution.  See function for a
		precise definition.
		ABSTOL is DOUBLE PRECISION
		RELTOL is DOUBLE PRECISION
		Method
		Algorithm R of the paper 'Two Efficient Algorithms with
		Guaranteed Convergence for Finding a Zero of a Function'
		by J. C. P. Bus and T. J. Dekker in ACM Transactions on
		Mathematical Software, Volume 1, no. 4 page 330
		(Dec. '75) is employed to find the zero of F(X)-Y.
		**********************************************************************
		*/
	{
		//	E0001_v2(1, NULL, NULL, NULL, NULL, NULL, NULL, NULL, zabstl, zreltl, zxhi, zxlo);
		E0001(1, NULL, NULL, W0, NULL, NULL, NULL, NULL);
	}

	/* stirlerr(n) = log(n!) - log( sqrt(2*pi*n)*(n/e)^n )
	*             = log Gamma(n+1) - 1/2 * [log(2*pi) + log(n)] - n*[log(n) - 1]
	*             = log Gamma(n+1) - (n + 1/2) * log(n) + n - log(2*pi)/2
	*
	* see also lgammacor() in ./lgammacor.c  which computes almost the same!
	*/

#include <limits>
#include <float.h>
	using namespace std;

	double bd0(double x, double np)
	{
		double ej, s, s1, v;
		int j;

		if(x == numeric_limits<double>::infinity() ||
			np == numeric_limits<double>::infinity() || np == 0.0)
			return numeric_limits<double>::quiet_NaN();

		if (fabs(x-np) < 0.1*(x+np)) {
			v = (x-np)/(x+np);
			s = (x-np)*v;/* s using v -- change by MM */
			ej = 2*x*v;
			v = v*v;
			for (j=1; ; j++) { /* Taylor series */
				ej *= v;
				s1 = s+ej/((j<<1)+1);
				if (s1==s) /* last term was effectively 0 */
					return(s1);
				s = s1;
			}
		}
		/* else:  | x - np |  is not too small */
		return(x*log(x/np)+np-x);
	}

	double dpois_raw(double x, double lambda, int give_log)
	{
#define M_2PI 6.283185307179586476925286766559 
#define R_D_fexp(f,x)     (give_log ? -0.5*log(f)+(x) : exp(x)/sqrt(f))
		/*       x >= 0 ; integer for dpois(), but not e.g. for pgamma()!
		lambda >= 0
		*/
		if (lambda == 0) return( (x == 0) ? 1.0 : 0.0 );
		if (lambda == numeric_limits<double>::infinity()) return 0.0;
		if (x < 0) return 0.0;
		if (x <= lambda * DBL_MIN) return(exp(-lambda) );

		if (lambda < x * DBL_MIN) {
			double xp1 = x+1;
			return(exp(-lambda + x*log(lambda) - x_alngam(&xp1)));
		}
		return(R_D_fexp( M_2PI*x, -x_stirlerr(x)-bd0(x,lambda) ));
	}

#define ML_POSINF numeric_limits<double>::infinity()
#define ML_NEGINF -numeric_limits<double>::infinity()

	double dgamma(double x, double shape, double scale, int give_log)
	{
		double pr;

		if (shape < 0 || scale <= 0)
			return numeric_limits<double>::quiet_NaN();
		if (x < 0)
			return 0.0;
		if (shape == 0) /* point mass at 0 */
			return (x == 0)? ML_POSINF : 0.0;
		if (x == 0) {
			if (shape < 1) return ML_POSINF;
			if (shape > 1) return 0.0;
			/* else */
			return give_log ? -log(scale) : 1 / scale;
		}

		if (shape < 1) {
			pr = dpois_raw(shape, x/scale, give_log);
			return give_log ?  pr + log(shape/x) : pr*shape/x;
		}
		/* else  shape >= 1 */
		pr = dpois_raw(shape-1, x/scale, give_log);
		return give_log ? pr - log(scale) : pr/scale;
	}

	double dchisq(double x, double df, int give_log)
	{
		return dgamma(x, df / 2., 2., give_log);
	}

	double dnorm(double x, double mu, double sigma, int give_log)
	{
#define LN_SQRT_2PI 0.918938533204672741780329736406
#define ONE_SQRT_2PI 0.398942280401432677939946059934
		if (sigma == numeric_limits<double>::infinity())
			return 0.0;
		if (x == numeric_limits<double>::infinity() && mu == x)
			return numeric_limits<double>::quiet_NaN(); /* x-mu is NaN */
		if (sigma <= 0) {
			if (sigma < 0) return numeric_limits<double>::signaling_NaN();
			/* sigma == 0 */
			return (x == mu) ? numeric_limits<double>::infinity() : 0.0;
		}
		x = (x - mu) / sigma;

		if(x == numeric_limits<double>::infinity()) return 0.0;
		return (give_log ?
			-(LN_SQRT_2PI  +	0.5 * x * x + log(sigma)) :
			ONE_SQRT_2PI * exp(-0.5 * x * x)  /	  sigma);
	}

#define EPS double(1.0e-10)
#define FPMIN double(1.0e-60)

	double betacf(double a, double b, 
		double x, int MAXIT=100);

	double betai(double a,double b,double x)
	{
		double bt;

		if (x < 0.0 || x > 1.0) halt("Bad x in routine betai");
		if (x == 0.0 || x == 1.0) bt=double(0.0);
		else {
			double ab = a+b;
			bt=exp(x_alngam(&ab)-x_alngam(&a)-x_alngam(&b)+a*log(x)+b*log(1.0-x));
		}
		if (x < (a+1.0)/(a+b+2.0))
			return bt*betacf(a,b,x)/a;
		else
			return 1.0-bt*betacf(b,a,1.0-x)/b;
	}


	double betacf(double a, double b,double x,int MAXIT)
	{
		int m,m2;
		double aa,c,d,del,h,qab,qam,qap;

		qab=a+b;
		qap=a+1.0;
		qam=a-1.0;
		c=1.0;
		d=1.0-qab*x/qap;
		if (fabs(d) < FPMIN) d=FPMIN;
		d=1.0/d;
		h=d;
		for (m=1;m<=MAXIT;m++) 
		{
			m2=2*m;
			aa=m*(b-m)*x/((qam+m2)*(a+m2));
			d=1.0+aa*d;
			if (fabs(d) < FPMIN) d=FPMIN;
			c=1.0+aa/c;
			if (fabs(c) < FPMIN) c=FPMIN;
			d=1.0/d;
			h *= d*c;
			aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
			d=1.0+aa*d;
			if (fabs(d) < FPMIN) d=FPMIN;
			c=1.0+aa/c;
			if (fabs(c) < FPMIN) c=FPMIN;
			d=1.0/d;

			del=d*c;
			h *= del;
			if (fabs(del-1.0) < EPS) break;
		}
		if (m > MAXIT)
			halt("a or b too big, or MAXIT too small in betacf");
		return h;
	}

	// double betacf(double a, double b, double x)
	// {
	// 	int which = 1, st;
	// 	double p, q, bd;
	// 	double y = 1-x;
	// 	cdfbet(&which, &p, &q, &x, &y, &a, &b, &st, &bd);
	// 
	// 	return p;
	// }

#undef MAXIT
#undef EPS
#undef FPMIN

#define R_D_nonint(x)   (fabs((x) - floor((x)+0.5)) > 1e-7)
#define R_D_negInonint(x)   (x < 0. || R_D_nonint(x))
#define R_D_nonint_check(x)
#define R_D_forceint(x)   floor((x) + 0.5)
#define	R_D_exp(x)	(give_log	?  (x)	 : exp(x))
#define	R_D_exp2(x)	(log_p	?  (x)	 : exp(x))

	wsReal ptdist(wsReal *Rp_tStat, wsReal *Rp_tDF)
	{
		double R_tDF = (double)(*Rp_tDF);
		double R_tStat = (double)(*Rp_tStat);
		return (wsReal)betai(0.5*R_tDF, 0.5, R_tDF/(R_tDF+SQR(R_tStat)));
	}

	double dbinom_raw(double x, double n, double p, double q, int give_log)
	{
		double lf, lc;

		if (p == 0) return((x == 0) ? 1.0 : 0.0);
		if (q == 0) return((x == n) ? 1.0 : 0.0);

		if (x == 0) {
			if(n == 0) return 1.0;
			lc = (p < 0.1) ? -bd0(n,n*q) - n*p : n*log(q);
			return( R_D_exp(lc) );
		}
		if (x == n) {
			lc = (q < 0.1) ? -bd0(n,n*p) - n*q : n*log(p);
			return( R_D_exp(lc) );
		}
		if (x < 0 || x > n) return( 0.0 );

		/* n*p or n*q can underflow to zero if n and p or q are small.  This
		used to occur in dbeta, and gives NaN as from R 2.3.0.  */
		lc = x_stirlerr(n) - x_stirlerr(x) - x_stirlerr(n-x) - bd0(x,n*p) - bd0(n-x,n*q);

		/* f = (M_2PI*x*(n-x))/n; could overflow or underflow */
		/* Upto R 2.7.1:
		* lf = log(M_2PI) + log(x) + log(n-x) - log(n);
		* -- following is much better for  x << n : */
		lf = log(M_2PI) + log(x) + util_log1p(- x/n);

		return R_D_exp(lc - 0.5*lf);
	}

	double dbinom(double x, double n, double p, int give_log/*=0*/)
	{
		if (p < 0 || p > 1 || R_D_negInonint(n))
			return numeric_limits<double>::quiet_NaN();
		R_D_nonint_check(x);
		if (x < 0 || numeric_limits<double>::infinity()==x) return 0.0;

		n = R_D_forceint(n);
		x = R_D_forceint(x);

		return dbinom_raw(x, n, p, 1-p, give_log);
	}


	/*
	*  Mathlib : A C Library of Special Functions
	*  Copyright (C) 2005-6 Morten Welinder <terra@gnome.org>
	*  Copyright (C) 2005-10 The R Foundation
	*  Copyright (C) 2006-10 The R Core Team
	*
	*  This program is free software; you can redistribute it and/or modify
	*  it under the terms of the GNU General Public License as published by
	*  the Free Software Foundation; either version 2 of the License, or
	*  (at your option) any later version.
	*
	*  This program is distributed in the hope that it will be useful,
	*  but WITHOUT ANY WARRANTY; without even the implied warranty of
	*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	*  GNU General Public License for more details.
	*
	*  You should have received a copy of the GNU General Public License
	*  along with this program; if not, a copy is available at
	*  http://www.r-project.org/Licenses/
	*
	*  SYNOPSIS
	*
	*	#include <Rmath.h>
	*
	*	double pgamma (double x, double alph, double scale,
	*		       int lower_tail, int log_p)
	*
	*	double __log1pmx	(double x)
	*	double lgamma1p (double a)
	*
	*	double logspace_add (double logx, double logy)
	*	double logspace_sub (double logx, double logy)
	*
	*
	*  DESCRIPTION
	*
	*	This function computes the distribution function for the
	*	gamma distribution with shape parameter alph and scale parameter
	*	scale.	This is also known as the incomplete gamma function.
	*	See Abramowitz and Stegun (6.5.1) for example.
	*
	*  NOTES
	*
	*	Complete redesign by Morten Welinder, originally for Gnumeric.
	*	Improvements (e.g. "while NEEDED_SCALE") by Martin Maechler
	*
	*  REFERENCES
	*
	*/


	/*----------- DEBUGGING -------------
	*	make CFLAGS='-DDEBUG_p -g -I/usr/local/include -I../include'
	* (cd ~/R/D/r-devel/Linux-inst/src/nmath; gcc -std=gnu99 -I. -I../../src/include -I../../../R/src/include -I/usr/local/include -DDEBUG_p -g -O2 -c ../../../R/src/nmath/pgamma.c -o pgamma.o)
	*/

	/* Scalefactor:= (2^32)^8 = 2^256 = 1.157921e+77 */
	//#define SQR(x) ((x)*(x))
	static const double scalefactor = SQR(SQR(SQR(4294967296.0)));
#undef SQR

	/* If |x| > |k| * M_cutoff,  then  log[ exp(-x) * k^x ]	 =~=  -x */
#ifndef M_LN2
#define M_LN2      0.693147180559945309417
#endif
	static const double M_cutoff = M_LN2 * DBL_MAX_EXP / DBL_EPSILON;/*=3.196577e18*/

																	 /* Continued fraction for calculation of
																	 *    1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3*d) + ... = sum_{k=0}^Inf x^k/(i+k*d)
																	 *
																	 * auxilary in __log1pmx() and lgamma1p()
																	 */
	static double
		logcf (double x, double i, double d,
			double eps /* ~ relative tolerance */)
	{
		double c1 = 2 * d;
		double c2 = i + d;
		double c4 = c2 + d;
		double a1 = c2;
		double b1 = i * (c2 - i * x);
		double b2 = d * d * x;
		double a2 = c4 * c2 - b2;

#if 0
		assert (i > 0);
		assert (d >= 0);
#endif

		b2 = c4 * b1 - i * b2;

		while (fabs(a2 * b1 - a1 * b2) > fabs(eps * b1 * b2)) {
			double c3 = c2*c2*x;
			c2 += d;
			c4 += d;
			a1 = c4 * a2 - c3 * a1;
			b1 = c4 * b2 - c3 * b1;

			c3 = c1 * c1 * x;
			c1 += d;
			c4 += d;
			a2 = c4 * a1 - c3 * a2;
			b2 = c4 * b1 - c3 * b2;

			if (fabs (b2) > scalefactor) {
				a1 /= scalefactor;
				b1 /= scalefactor;
				a2 /= scalefactor;
				b2 /= scalefactor;
			} else if (fabs (b2) < 1 / scalefactor) {
				a1 *= scalefactor;
				b1 *= scalefactor;
				a2 *= scalefactor;
				b2 *= scalefactor;
			}
		}

		return a2 / b2;
	}

	/* Accurate calculation of log(1+x)-x, particularly for small x.  */
	double __log1pmx (double x)
	{
		static const double minLog1Value = -0.79149064;

		if (x > 1 || x < minLog1Value)
			return util_log1p(x) - x;
		else { /* -.791 <=  x <= 1  -- expand in  [x/(2+x)]^2 =: y :
			   * log(1+x) - x =  x/(2+x) * [ 2 * y * S(y) - x],  with
			   * ---------------------------------------------
			   * S(y) = 1/3 + y/5 + y^2/7 + ... = \sum_{k=0}^\infty  y^k / (2k + 3)
			   */
			double r = x / (2 + x), y = r * r;
			if (fabs(x) < 1e-2) {
				static const double two = 2;
				return r * ((((two / 9 * y + two / 7) * y + two / 5) * y +
					two / 3) * y - x);
			} else {
				static const double tol_logcf = 1e-14;
				return r * (2 * y * logcf (y, 3, 2, tol_logcf) - x);
			}
		}
	}


	/* Compute  log(gamma(a+1))  accurately also for small a (0 < a < 0.5). */
	double lgamma1p (double a)
	{
		const double eulers_const =	 0.5772156649015328606065120900824024;

		/* coeffs[i] holds (zeta(i+2)-1)/(i+2) , i = 0:(N-1), N = 40 : */
		const int N = 40;
		static const double coeffs[40] = {
			0.3224670334241132182362075833230126e-0,/* = (zeta(2)-1)/2 */
			0.6735230105319809513324605383715000e-1,/* = (zeta(3)-1)/3 */
			0.2058080842778454787900092413529198e-1,
			0.7385551028673985266273097291406834e-2,
			0.2890510330741523285752988298486755e-2,
			0.1192753911703260977113935692828109e-2,
			0.5096695247430424223356548135815582e-3,
			0.2231547584535793797614188036013401e-3,
			0.9945751278180853371459589003190170e-4,
			0.4492623673813314170020750240635786e-4,
			0.2050721277567069155316650397830591e-4,
			0.9439488275268395903987425104415055e-5,
			0.4374866789907487804181793223952411e-5,
			0.2039215753801366236781900709670839e-5,
			0.9551412130407419832857179772951265e-6,
			0.4492469198764566043294290331193655e-6,
			0.2120718480555466586923135901077628e-6,
			0.1004322482396809960872083050053344e-6,
			0.4769810169363980565760193417246730e-7,
			0.2271109460894316491031998116062124e-7,
			0.1083865921489695409107491757968159e-7,
			0.5183475041970046655121248647057669e-8,
			0.2483674543802478317185008663991718e-8,
			0.1192140140586091207442548202774640e-8,
			0.5731367241678862013330194857961011e-9,
			0.2759522885124233145178149692816341e-9,
			0.1330476437424448948149715720858008e-9,
			0.6422964563838100022082448087644648e-10,
			0.3104424774732227276239215783404066e-10,
			0.1502138408075414217093301048780668e-10,
			0.7275974480239079662504549924814047e-11,
			0.3527742476575915083615072228655483e-11,
			0.1711991790559617908601084114443031e-11,
			0.8315385841420284819798357793954418e-12,
			0.4042200525289440065536008957032895e-12,
			0.1966475631096616490411045679010286e-12,
			0.9573630387838555763782200936508615e-13,
			0.4664076026428374224576492565974577e-13,
			0.2273736960065972320633279596737272e-13,
			0.1109139947083452201658320007192334e-13/* = (zeta(40+1)-1)/(40+1) */
		};

		const double c = 0.2273736845824652515226821577978691e-12;/* zeta(N+2)-1 */
		const double tol_logcf = 1e-14;
		double lgam;
		int i;

		if (fabs (a) >= 0.5)
			return lgammafn (a + 1);

		/* Abramowitz & Stegun 6.1.33 : for |x| < 2,
		* <==> log(gamma(1+x)) = -(log(1+x) - x) - gamma*x + x^2 * \sum_{n=0}^\infty c_n (-x)^n
		* where c_n := (Zeta(n+2) - 1)/(n+2)  = coeffs[n]
		*
		* Here, another convergence acceleration trick is used to compute
		* lgam(x) :=  sum_{n=0..Inf} c_n (-x)^n
		*/
		lgam = c * logcf(-a / 2, N + 2, 1, tol_logcf);
		for (i = N - 1; i >= 0; i--)
			lgam = coeffs[i] - a * lgam;

		return (a * lgam - eulers_const) * a - __log1pmx (a);
	} /* lgamma1p */



	  /*
	  * Compute the log of a sum from logs of terms, i.e.,
	  *
	  *     log (exp (logx) + exp (logy))
	  *
	  * without causing overflows and without throwing away large handfuls
	  * of accuracy.
	  */
	double logspace_add (double logx, double logy)
	{
		return fmax2 (logx, logy) + util_log1p (exp (-fabs (logx - logy)));
	}


	/*
	* Compute the log of a difference from logs of terms, i.e.,
	*
	*     log (exp (logx) - exp (logy))
	*
	* without causing overflows and without throwing away large handfuls
	* of accuracy.
	*/
	double logspace_sub (double logx, double logy)
	{
		return logx + R_Log1_Exp(logy - logx);
	}


	/* dpois_wrap (x_P_1,  lambda, g_log) ==
	*   dpois (x_P_1 - 1, lambda, g_log) :=  exp(-L)  L^k / gamma(k+1) ,  k := x_P_1 - 1
	*/
	static double
		dpois_wrap (double x_plus_1, double lambda, int give_log)
	{
#ifdef DEBUG_p
		REprintf (" dpois_wrap(x+1=%.14g, lambda=%.14g, log=%d)\n",
			x_plus_1, lambda, give_log);
#endif
		if (lambda == numeric_limits<double>::infinity() ||
			-lambda == numeric_limits<double>::infinity())
			return give_log ? ML_NEGINF : 0.;
		//return R_D__0;
		if (x_plus_1 > 1)
			return dpois_raw (x_plus_1 - 1, lambda, give_log);
		if (lambda > fabs(x_plus_1 - 1) * M_cutoff)
			return R_D_exp(-lambda - lgammafn(x_plus_1));
		else {
			double d = dpois_raw (x_plus_1, lambda, give_log);
#ifdef DEBUG_p
			REprintf ("  -> d=dpois_raw(..)=%.14g\n", d);
#endif
			return give_log
				? d + log (x_plus_1 / lambda)
				: d * (x_plus_1 / lambda);
		}
	}

	/*
	* Abramowitz and Stegun 6.5.29 [right]
	*/
	static double
		pgamma_smallx (double x, double alph, int lower_tail, int log_p)
	{
		double sum = 0, c = alph, n = 0, term;

#ifdef DEBUG_p
		REprintf (" pg_smallx(x=%.12g, alph=%.12g): ", x, alph);
#endif

		/*
		* Relative to 6.5.29 all terms have been multiplied by alph
		* and the first, thus being 1, is omitted.
		*/

		do {
			n++;
			c *= -x / n;
			term = c / (alph + n);
			sum += term;
		} while (fabs (term) > DBL_EPSILON * fabs (sum));

#ifdef DEBUG_p
		REprintf (" %d terms --> conv.sum=%g;", n, sum);
#endif
		if (lower_tail) {
			double f1 = log_p ? util_log1p (sum) : 1 + sum;
			double f2;
			if (alph > 1) {
				f2 = dpois_raw (alph, x, log_p);
				f2 = log_p ? f2 + x : f2 * exp (x);
			} else if (log_p)
				f2 = alph * log (x) - lgamma1p (alph);
			else
				f2 = pow (x, alph) / exp (lgamma1p (alph));
#ifdef DEBUG_p
			REprintf (" (f1,f2)= (%g,%g)\n", f1,f2);
#endif
			return log_p ? f1 + f2 : f1 * f2;
		} else {
			double lf2 = alph * log (x) - lgamma1p (alph);
#ifdef DEBUG_p
			REprintf (" 1:%.14g  2:%.14g\n", alph * log (x), lgamma1p (alph));
			REprintf (" sum=%.14g  log(1+sum)=%.14g	 lf2=%.14g\n",
				sum, log1p (sum), lf2);
#endif
			if (log_p)
				return R_Log1_Exp (util_log1p (sum) + lf2);
			else {
				double f1m1 = sum;
				double f2m1 = util_expm1 (lf2);
				return -(f1m1 + f2m1 + f1m1 * f2m1);
			}
		}
	} /* pgamma_smallx() */

	static double
		pd_upper_series (double x, double y, int log_p)
	{
		double term = x / y;
		double sum = term;

		do {
			y++;
			term *= x / y;
			sum += term;
		} while (term > sum * DBL_EPSILON);

		/* sum =  \sum_{n=1}^ oo  x^n     / (y*(y+1)*...*(y+n-1))
		*	   =  \sum_{n=0}^ oo  x^(n+1) / (y*(y+1)*...*(y+n))
		*	   =  x/y * (1 + \sum_{n=1}^oo	x^n / ((y+1)*...*(y+n)))
		*	   ~  x/y +  o(x/y)   {which happens when alph -> Inf}
		*/
		return log_p ? log (sum) : sum;
	}

	/* Continued fraction for calculation of
	*    scaled upper-tail F_{gamma}
	*  ~=  (y / d) * [1 +  (1-y)/d +  O( ((1-y)/d)^2 ) ]
	*/
	static double
		pd_lower_cf (double y, double d)
	{
		double f= 0.0 /* -Wall */, of, f0;
		double i, c2, c3, c4,  a1, b1,  a2, b2;

#define	NEEDED_SCALE				\
	  (b2 > scalefactor) {			\
	    a1 /= scalefactor;			\
	    b1 /= scalefactor;			\
	    a2 /= scalefactor;			\
	    b2 /= scalefactor;			\
	}

#define max_it 200000

#ifdef DEBUG_p
		REprintf("pd_lower_cf(y=%.14g, d=%.14g)", y, d);
#endif
		if (y == 0) return 0;

		f0 = y/d;
		/* Needed, e.g. for  pgamma(10^c(100,295), shape= 1.1, log=TRUE): */
		if(fabs(y - 1) < fabs(d) * DBL_EPSILON) { /* includes y < d = Inf */
#ifdef DEBUG_p
			REprintf(" very small 'y' -> returning (y/d)\n");
#endif
			return (f0);
		}

		if(f0 > 1.) f0 = 1.;
		c2 = y;
		c4 = d; /* original (y,d), *not* potentially scaled ones!*/

		a1 = 0; b1 = 1;
		a2 = y; b2 = d;

		while NEEDED_SCALE

			i = 0; of = -1.; /* far away */
		while (i < max_it) {

			i++;	c2--;	c3 = i * c2;	c4 += 2;
			/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i odd */
			a1 = c4 * a2 + c3 * a1;
			b1 = c4 * b2 + c3 * b1;

			i++;	c2--;	c3 = i * c2;	c4 += 2;
			/* c2 = y - i,  c3 = i(y - i),  c4 = d + 2i,  for i even */
			a2 = c4 * a1 + c3 * a2;
			b2 = c4 * b1 + c3 * b2;

			if NEEDED_SCALE

				if (b2 != 0) {
					f = a2 / b2;
					/* convergence check: relative; "absolute" for very small f : */
					if (fabs (f - of) <= DBL_EPSILON * fmax2(f0, fabs(f))) {
#ifdef DEBUG_p
						REprintf(" %g iter.\n", i);
#endif
						return f;
					}
					of = f;
				}
		}

		//    MATHLIB_WARNING(" ** NON-convergence in pgamma()'s pd_lower_cf() f= %g.\n",
		//		    f);
		return f;/* should not happen ... */
	} /* pd_lower_cf() */
#undef NEEDED_SCALE


	static double
		pd_lower_series (double lambda, double y)
	{
		double term = 1, sum = 0;

#ifdef DEBUG_p
		REprintf("pd_lower_series(lam=%.14g, y=%.14g) ...", lambda, y);
#endif
		while (y >= 1 && term > sum * DBL_EPSILON) {
			term *= y / lambda;
			sum += term;
			y--;
		}
		/* sum =  \sum_{n=0}^ oo  y*(y-1)*...*(y - n) / lambda^(n+1)
		*	   =  y/lambda * (1 + \sum_{n=1}^Inf  (y-1)*...*(y-n) / lambda^n)
		*	   ~  y/lambda + o(y/lambda)
		*/
#ifdef DEBUG_p
		REprintf(" done: term=%g, sum=%g, y= %g\n", term, sum, y);
#endif

		if (y != floor (y)) {
			/*
			* The series does not converge as the terms start getting
			* bigger (besides flipping sign) for y < -lambda.
			*/
			double f;
#ifdef DEBUG_p
			REprintf(" y not int: add another term ");
#endif
			/* FIXME: in quite few cases, adding  term*f  has no effect (f too small)
			*	  and is unnecessary e.g. for pgamma(4e12, 121.1) */
			f = pd_lower_cf (y, lambda + 1 - y);
#ifdef DEBUG_p
			REprintf("  (= %.14g) * term = %.14g to sum %g\n", f, term * f, sum);
#endif
			sum += term * f;
		}

		return sum;
	} /* pd_lower_series() */

	  /*
	  * Compute the following ratio with higher accuracy that would be had
	  * from doing it directly.
	  *
	  *		 dnorm (x, 0, 1, FALSE)
	  *	   ----------------------------------
	  *	   pnorm (x, 0, 1, lower_tail, FALSE)
	  *
	  * Abramowitz & Stegun 26.2.12
	  */
	static double
		dpnorm (double x, int lower_tail, double lp)
	{
		/*
		* So as not to repeat a pnorm call, we expect
		*
		*	 lp == pnorm (x, 0, 1, lower_tail, TRUE)
		*
		* but use it only in the non-critical case where either x is small
		* or p==exp(lp) is close to 1.
		*/

		if (x < 0) {
			x = -x;
			lower_tail = !lower_tail;
		}

		if (x > 10 && !lower_tail) {
			double term = 1 / x;
			double sum = term;
			double x2 = x * x;
			double i = 1;

			do {
				term *= -i / x2;
				sum += term;
				i += 2;
			} while (fabs (term) > DBL_EPSILON * sum);

			return 1 / sum;
		} else {
			double d = dnorm (x, 0., 1., 0);
			return d / exp (lp);
		}
	}

	/*
	* Asymptotic expansion to calculate the probability that Poisson variate
	* has value <= x.
	* Various assertions about this are made (without proof) at
	* http://members.aol.com/iandjmsmith/PoissonApprox.htm
	*/
	static double
		ppois_asymp (double x, double lambda, int lower_tail, int log_p)
	{
		static const double coefs_a[8] = {
			-1e99, /* placeholder used for 1-indexing */
			2/3.,
			-4/135.,
			8/2835.,
			16/8505.,
			-8992/12629925.,
			-334144/492567075.,
			698752/1477701225.
		};

		static const double coefs_b[8] = {
			-1e99, /* placeholder */
			1/12.,
			1/288.,
			-139/51840.,
			-571/2488320.,
			163879/209018880.,
			5246819/75246796800.,
			-534703531/902961561600.
		};

		double elfb, elfb_term;
		double res12, res1_term, res1_ig, res2_term, res2_ig;
		double dfm, pt_, s2pt, f, np;
		int i;

		dfm = lambda - x;
		/* If lambda is large, the distribution is highly concentrated
		about lambda.  So representation error in x or lambda can lead
		to arbitrarily large values of pt_ and hence divergence of the
		coefficients of this approximation.
		*/
		pt_ = - __log1pmx (dfm / x);
		s2pt = sqrt (2 * x * pt_);
		if (dfm < 0) s2pt = -s2pt;

		res12 = 0;
		res1_ig = res1_term = sqrt (x);
		res2_ig = res2_term = s2pt;
		for (i = 1; i < 8; i++) {
			res12 += res1_ig * coefs_a[i];
			res12 += res2_ig * coefs_b[i];
			res1_term *= pt_ / i ;
			res2_term *= 2 * pt_ / (2 * i + 1);
			res1_ig = res1_ig / x + res1_term;
			res2_ig = res2_ig / x + res2_term;
		}

		elfb = x;
		elfb_term = 1;
		for (i = 1; i < 8; i++) {
			elfb += elfb_term * coefs_b[i];
			elfb_term /= x;
		}
		if (!lower_tail) elfb = -elfb;
#ifdef DEBUG_p
		REprintf ("res12 = %.14g   elfb=%.14g\n", elfb, res12);
#endif

		f = res12 / elfb;

		np = pnorm (s2pt, 0.0, 1.0, !lower_tail, log_p);

		if (log_p) {
			double n_d_over_p = dpnorm (s2pt, !lower_tail, np);
#ifdef DEBUG_p
			REprintf ("pp*_asymp(): f=%.14g	 np=e^%.14g  nd/np=%.14g  f*nd/np=%.14g\n",
				f, np, n_d_over_p, f * n_d_over_p);
#endif
			return np + util_log1p (f * n_d_over_p);
		} else {
			double nd = dnorm (s2pt, 0., 1., log_p);

#ifdef DEBUG_p
			REprintf ("pp*_asymp(): f=%.14g	 np=%.14g  nd=%.14g  f*nd=%.14g\n",
				f, np, nd, f * nd);
#endif
			return np + f * nd;
		}
	} /* ppois_asymp() */

	double pgamma_raw (double x, double alph, int lower_tail, int log_p)
	{
		/* Here, assume that  (x,alph) are not NA  &  alph > 0 . */

		double res;

#ifdef DEBUG_p
		REprintf("pgamma_raw(x=%.14g, alph=%.14g, low=%d, log=%d)\n",
			x, alph, lower_tail, log_p);
#endif
		R_P_bounds_01(x, 0., ML_POSINF);

		if (x < 1) {
			res = pgamma_smallx (x, alph, lower_tail, log_p);
		} else if (x <= alph - 1 && x < 0.8 * (alph + 50)) {
			/* incl. large alph compared to x */
			double sum = pd_upper_series (x, alph, log_p);/* = x/alph + o(x/alph) */
			double d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
			REprintf(" alph 'large': sum=pd_upper*()= %.12g, d=dpois_w(*)= %.12g\n",
				sum, d);
#endif
			if (!lower_tail)
				res = log_p
				? R_Log1_Exp (d + sum)
				: 1 - d * sum;
			else
				res = log_p ? sum + d : sum * d;
		} else if (alph - 1 < x && alph < 0.8 * (x + 50)) {
			/* incl. large x compared to alph */
			double sum;
			double d = dpois_wrap (alph, x, log_p);
#ifdef DEBUG_p
			REprintf(" x 'large': d=dpois_w(*)= %.14g ", d);
#endif
			if (alph < 1) {
				if (x * DBL_EPSILON > 1 - alph)
					sum = R_D__1;
				else {
					double f = pd_lower_cf (alph, x - (alph - 1)) * x / alph;
					/* = [alph/(x - alph+1) + o(alph/(x-alph+1))] * x/alph = 1 + o(1) */
					sum = log_p ? log (f) : f;
				}
			} else {
				sum = pd_lower_series (x, alph - 1);/* = (alph-1)/x + o((alph-1)/x) */
				sum = log_p ? util_log1p (sum) : 1 + sum;
			}
#ifdef DEBUG_p
			REprintf(", sum= %.14g\n", sum);
#endif
			if (!lower_tail)
				res = log_p ? sum + d : sum * d;
			else
				res = log_p
				? R_Log1_Exp (d + sum)
				: 1 - d * sum;
		} else { /* x >= 1 and x fairly near alph. */
#ifdef DEBUG_p
			REprintf(" using ppois_asymp()\n");
#endif
			res = ppois_asymp (alph - 1, x, !lower_tail, log_p);
		}

		/*
		* We lose a fair amount of accuracy to underflow in the cases
		* where the final result is very close to DBL_MIN.	 In those
		* cases, simply redo via log space.
		*/
		if (!log_p && res < DBL_MIN / DBL_EPSILON) {
			/* with(.Machine, double.xmin / double.eps) #|-> 1.002084e-292 */
#ifdef DEBUG_p
			REprintf(" very small res=%.14g; -> recompute via log\n", res);
#endif
			return exp (pgamma_raw (x, alph, lower_tail, 1));
		} else
			return res;
	}


	double pgamma(double x, double alph, double scale, int lower_tail, int log_p)
	{
#ifdef IEEE_754
		if (ISNAN(x) || ISNAN(alph) || ISNAN(scale))
			return x + alph + scale;
#endif
		if(alph < 0. || scale <= 0.)
			return numeric_limits<double>::quiet_NaN();
		x /= scale;
#ifdef IEEE_754
		if (ISNAN(x)) /* eg. original x = scale = +Inf */
			return x;
#endif
		if(alph == 0.) /* limit case; useful e.g. in pnchisq() */
			return (x <= 0) ? R_DT_0: R_DT_1; /* <= assert  pgamma(0,0) ==> 0 */
		return pgamma_raw (x, alph, lower_tail, log_p);
	}
	/* From: terra@gnome.org (Morten Welinder)
	* To: R-bugs@biostat.ku.dk
	* Cc: maechler@stat.math.ethz.ch
	* Subject: Re: [Rd] pgamma discontinuity (PR#7307)
	* Date: Tue, 11 Jan 2005 13:57:26 -0500 (EST)

	* this version of pgamma appears to be quite good and certainly a vast
	* improvement over current R code.  (I last looked at 2.0.1)  Apart from
	* type naming, this is what I have been using for Gnumeric 1.4.1.

	* This could be included into R as-is, but you might want to benefit from
	* making logcf, __log1pmx, lgamma1p, and possibly logspace_add/logspace_sub
	* available to other parts of R.

	* MM: I've not (yet?) taken  logcf(), but the other four
	*/

	/*
	*  Mathlib : A C Library of Special Functions
	*  Copyright (C) 1998 Ross Ihaka
	*  Copyright (C) 2000--2011 The R Core Team
	*  Copyright (C) 2004--2009 The R Foundation
	*  based on AS 91 (C) 1979 Royal Statistical Society
	*
	*  This program is free software; you can redistribute it and/or modify
	*  it under the terms of the GNU General Public License as published by
	*  the Free Software Foundation; either version 2 of the License, or
	*  (at your option) any later version.
	*
	*  This program is distributed in the hope that it will be useful, but
	*  WITHOUT ANY WARRANTY; without even the implied warranty of
	*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	*  General Public License for more details.
	*
	*  You should have received a copy of the GNU General Public License
	*  along with this program; if not, a copy is available at
	*  http://www.r-project.org/Licenses/
	*
	*  DESCRIPTION
	*
	*	Compute the quantile function of the gamma distribution.
	*
	*  NOTES
	*
	*	This function is based on the Applied Statistics
	*	Algorithm AS 91 ("ppchi2") and via pgamma(.) AS 239.
	*
	*	R core improvements:
	*	o  lower_tail, log_p
	*      o  non-trivial result for p outside [0.000002, 0.999998]
	*	o  p ~ 1 no longer gives +Inf; final Newton step(s)
	*
	*  REFERENCES
	*
	*	Best, D. J. and D. E. Roberts (1975).
	*	Percentage Points of the Chi-Squared Distribution.
	*	Applied Statistics 24, page 385.  */

#define R_Q_P01_check(p)			\
	if ((log_p	&& p > 0) ||			\
	(!log_p && (p < 0 || p > 1)) )		\
	return numeric_limits<double>::quiet_NaN()

#ifdef DEBUG_qgamma
# define DEBUG_q
#endif

	/* log(1-exp(x)):  R_D_LExp(x) == (log1p(- R_D_qIv(x))) but even more stable:*/
#define R_D_LExp(x)     (log_p ? R_Log1_Exp(x) : util_log1p(-x))
#define R_D_log(p)	(log_p	?  (p)	 : log(p))	/* log(p) */
#define R_DT_log(p)	(lower_tail? R_D_log(p) : R_D_LExp(p))/* log(p) in qF */
#define R_DT_Clog(p)	(lower_tail? R_D_LExp(p): R_D_log(p))/* log(1-p) in qF*/

	double qchisq_appr(double p, double nu, double g /* = log Gamma(nu/2) */,
		int lower_tail, int log_p, double tol /* EPS1 */)
	{
#define C7	4.67
#define C8	6.66
#define C9	6.73
#define C10	13.32

		double alpha, a, c, ch, p1;
		double p2, q, t, x;

		/* test arguments and initialise */

#ifdef IEEE_754
		if (ISNAN(p) || ISNAN(nu))
			return p + nu;
#endif
		R_Q_P01_check(p);
		if (nu <= 0) return numeric_limits<double>::quiet_NaN();

		alpha = 0.5 * nu;/* = [pq]gamma() shape */
		c = alpha-1;

		if (nu < (-1.24)*(p1 = R_DT_log(p))) {	/* for small chi-squared */
												/* log(alpha) + g = log(alpha) + log(gamma(alpha)) =
												*        = log(alpha*gamma(alpha)) = lgamma(alpha+1) suffers from
												*  catastrophic cancellation when alpha << 1
												*/
			double lgam1pa = (alpha < 0.5) ? lgamma1p(alpha) : (log(alpha) + g);
			ch = exp((lgam1pa + p1)/alpha + M_LN2);
#ifdef DEBUG_qgamma
			REprintf(" small chi-sq., ch0 = %g\n", ch);
#endif

		} else if(nu > 0.32) {	/*  using Wilson and Hilferty estimate */

			x = qnorm(p, 0, 1, lower_tail, log_p);
			p1 = 2./(9*nu);
			ch = nu*pow(x*sqrt(p1) + 1-p1, 3);

#ifdef DEBUG_qgamma
			REprintf(" nu > .32: Wilson-Hilferty; x = %7g\n", x);
#endif
			/* approximation for p tending to 1: */
			if( ch > 2.2*nu + 6 )
				ch = -2*(R_DT_Clog(p) - c*log(0.5*ch) + g);

		} else { /* "small nu" : 1.24*(-log(p)) <= nu <= 0.32 */
			ch = 0.4;
			a = R_DT_Clog(p) + g + c*M_LN2;
#ifdef DEBUG_qgamma
			REprintf(" nu <= .32: a = %7g\n", a);
#endif
			do {
				q = ch;
				p1 = 1. / (1+ch*(C7+ch));
				p2 = ch*(C9+ch*(C8+ch));
				t = -0.5 +(C7+2*ch)*p1 - (C9+ch*(C10+3*ch))/p2;
				ch -= (1- exp(a+0.5*ch)*p2*p1)/t;
			} while(fabs(q - ch) > tol * fabs(ch));
		}

		return ch;
	}

	double qgamma(double p, double alpha, double scale, int lower_tail, int log_p)
		/*			shape = alpha */
	{
#define EPS1 1e-2
#define EPS2 5e-7/* final precision of AS 91 */
#define EPS_N 1e-15/* precision of Newton step / iterations */
#define LN_EPS -36.043653389117156 /* = log(.Machine$double.eps) iff IEEE_754 */

#define MAXIT 1000/* was 20 */

#define pMIN 1e-100   /* was 0.000002 = 2e-6 */
#define pMAX (1-1e-14)/* was (1-1e-12) and 0.999998 = 1 - 2e-6 */

		const static double
			i420  = 1./ 420.,
			i2520 = 1./ 2520.,
			i5040 = 1./ 5040;

		double p_, a, b, c, g, ch, ch0, p1;
		double p2, q, s1, s2, s3, s4, s5, s6, t, x;
		int i, max_it_Newton = 1;

		/* test arguments and initialise */

#ifdef IEEE_754
		if (ISNAN(p) || ISNAN(alpha) || ISNAN(scale))
			return p + alpha + scale;
#endif
		R_Q_P01_boundaries(p, 0., ML_POSINF);

		if (alpha < 0 || scale <= 0) return numeric_limits<double>::quiet_NaN();

		if (alpha == 0) /* all mass at 0 : */ return 0.;

		if (alpha < 1e-10) {
			/* Warning seems unnecessary now: */
#ifdef _DO_WARN_qgamma_
			MATHLIB_WARNING("value of shape (%g) is extremely small: results may be unreliable", alpha);
#endif
			max_it_Newton = 7;/* may still be increased below */
		}

		p_ = R_DT_qIv(p);/* lower_tail prob (in any case) */

#ifdef DEBUG_qgamma
		REprintf("qgamma(p=%7g, alpha=%7g, scale=%7g, l.t.=%2d, log_p=%2d): ",
			p,alpha,scale, lower_tail, log_p);
#endif
		g = lgammafn(alpha);/* log Gamma(v/2) */

							/*----- Phase I : Starting Approximation */
		ch = qchisq_appr(p, /* nu= 'df' =  */ 2*alpha, /* lgamma(nu/2)= */ g,
			lower_tail, log_p, /* tol= */ EPS1);
		if (ch == numeric_limits<double>::infinity() ||
			-ch == numeric_limits<double>::infinity()) {
			//    if(!R_FINITE(ch)) {
			/* forget about all iterations! */
			max_it_Newton = 0; goto END;
		}
		if(ch < EPS2) {/* Corrected according to AS 91; MM, May 25, 1999 */
			max_it_Newton = 20;
			goto END;/* and do Newton steps */
		}

		/* FIXME: This (cutoff to {0, +Inf}) is far from optimal
		* -----  when log_p or !lower_tail, but NOT doing it can be even worse */
		if(p_ > pMAX || p_ < pMIN) {
			/* did return ML_POSINF or 0.;	much better: */
			max_it_Newton = 20;
			goto END;/* and do Newton steps */
		}

#ifdef DEBUG_qgamma
		REprintf("\t==> ch = %10g:", ch);
#endif

		/*----- Phase II: Iteration
		*	Call pgamma() [AS 239]	and calculate seven term taylor series
		*/
		c = alpha-1;
		s6 = (120+c*(346+127*c)) * i5040; /* used below, is "const" */

		ch0 = ch;/* save initial approx. */
		for(i=1; i <= MAXIT; i++ ) {
			q = ch;
			p1 = 0.5*ch;
			p2 = p_ - pgamma_raw(p1, alpha, /*lower_tail*/TRUE, /*log_p*/0);
#ifdef DEBUG_qgamma
			if(i == 1) REprintf(" Ph.II iter; ch=%g, p2=%g\n", ch, p2);
			if(i >= 2) REprintf("     it=%d,  ch=%g, p2=%g\n", i, ch, p2);
#endif
#ifdef IEEE_754
			if(!R_FINITE(p2) || ch <= 0)
#else
			if(ch <= 0)
#endif
			{ ch = ch0; max_it_Newton = 27; goto END; }/*was  return ML_NAN;*/

			t = p2*exp(alpha*M_LN2+g+p1-c*log(ch));
			b = t/ch;
			a = 0.5*t - b*c;
			s1 = (210+ a*(140+a*(105+a*(84+a*(70+60*a))))) * i420;
			s2 = (420+ a*(735+a*(966+a*(1141+1278*a)))) * i2520;
			s3 = (210+ a*(462+a*(707+932*a))) * i2520;
			s4 = (252+ a*(672+1182*a) + c*(294+a*(889+1740*a))) * i5040;
			s5 = (84+2264*a + c*(1175+606*a)) * i2520;

			ch += t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
			if(fabs(q - ch) < EPS2*ch)
				goto END;
			if(fabs(q - ch) > 0.1*ch) {/* diverging? -- also forces ch > 0 */
				if(ch < q) ch = 0.9 * q; else ch = 1.1 * q;
			}
		}
		/* no convergence in MAXIT iterations -- but we add Newton now... */
#ifdef DEBUG_q
		MATHLIB_WARNING3("qgamma(%g) not converged in %d iterations; rel.ch=%g\n",
			p, MAXIT, ch/fabs(q - ch));
#endif
		/* was
		*    ML_ERROR(ME_PRECISION, "qgamma");
		* does nothing in R !*/

	END:
		/* PR# 2214 :	 From: Morten Welinder <terra@diku.dk>, Fri, 25 Oct 2002 16:50
		--------	 To: R-bugs@biostat.ku.dk     Subject: qgamma precision

		* With a final Newton step, double accuracy, e.g. for (p= 7e-4; nu= 0.9)
		*
		* Improved (MM): - only if rel.Err > EPS_N (= 1e-15);
		*		    - also for lower_tail = FALSE	 or log_p = TRUE
		* 		    - optionally *iterate* Newton
		*/
		x = 0.5*scale*ch;
		if(max_it_Newton) {
			/* always use log scale */
			if (!log_p) {
				p = log(p);
				log_p = TRUE;
			}
			if(x == 0) {
				const double _1_p = 1. + 1e-7;
				const double _1_m = 1. - 1e-7;
				x = DBL_MIN;
				p_ = pgamma(x, alpha, scale, lower_tail, log_p);
				if(( lower_tail && p_ > p * _1_p) ||
					(!lower_tail && p_ < p * _1_m))
					return(0.);
				/* else:  continue, using x = DBL_MIN instead of  0  */
			}
			else
				p_ = pgamma(x, alpha, scale, lower_tail, log_p);
			if(p_ == ML_NEGINF) return 0; /* PR#14710 */
			for(i = 1; i <= max_it_Newton; i++) {
				p1 = p_ - p;
#ifdef DEBUG_qgamma
				if(i == 1) REprintf("\n it=%d: p=%g, x = %g, p.=%g; p1=d{p}=%g\n",
					i, p, x, p_, p1);
				if(i >= 2) REprintf("          x{it= %d} = %g, p.=%g, p1=d{p}=%g\n",
					i,    x, p_, p1);
#endif
				if(fabs(p1) < fabs(EPS_N * p))
					break;
				/* else */
				if((g = dgamma(x, alpha, scale, log_p)) == R_D__0) {
#ifdef DEBUG_q
					if(i == 1) REprintf("no final Newton step because dgamma(*)== 0!\n");
#endif
					break;
				}
				/* else :
				* delta x = f(x)/f'(x);
				* if(log_p) f(x) := log P(x) - p; f'(x) = d/dx log P(x) = P' / P
				* ==> f(x)/f'(x) = f*P / P' = f*exp(p_) / P' (since p_ = log P(x))
				*/
				t = log_p ? p1*exp(p_ - g) : p1/g ;/* = "delta x" */
				t = lower_tail ? x - t : x + t;
				p_ = pgamma (t, alpha, scale, lower_tail, log_p);
				if (fabs(p_ - p) > fabs(p1) ||
					(i > 1 && fabs(p_ - p) == fabs(p1)) /* <- against flip-flop */) {
					/* no improvement */
#ifdef DEBUG_qgamma
					if(i == 1 && max_it_Newton > 1)
						REprintf("no Newton step done since delta{p} >= last delta\n");
#endif
					break;
				} /* else : */
#ifdef Harmful_notably_if_max_it_Newton_is_1
				  /* control step length: this could have started at
				  the initial approximation */
				if(t > 1.1*x) t = 1.1*x;
				else if(t < 0.9*x) t = 0.9*x;
#endif
				x = t;
			}
		}

		return x;
	}

	double pf(double x, double df1, double df2, int lower_tail)
	{
		int log_p = 0;
#ifdef IEEE_754
		if (ISNAN(x) || ISNAN(df1) || ISNAN(df2))
			return x + df2 + df1;
#endif
		if (df1 <= 0. || df2 <= 0.) return numeric_limits<double>::quiet_NaN();

		R_P_bounds_01(x, 0., ML_POSINF);

		/* move to pchisq for very large values - was 'df1 > 4e5' in 2.0.x,
		now only needed for df1 = Inf or df2 = Inf {since pbeta(0,*)=0} : */
		if (df2 == ML_POSINF) {
			if (df1 == ML_POSINF) {
				if(x <  1.) return R_DT_0;
				if(x == 1.) return 0.5;
				if(x >  1.) return R_DT_1;
			}
			double p, q;
			int st = 0; // error variable
			int w = 1; // function variable
			double bnd = 1; // boundary function
			double _x = x*df1;


			cdfchi(&w,&p,&q,&_x,&df1,&st,&bnd);
			return lower_tail?p:q;
			//return pchisq(x * df1, df1, lower_tail, log_p);
		}

		if (df1 == ML_POSINF) {/* was "fudge"	'df1 > 4e5' in 2.0.x */
			double p, q;
			int st = 0; // error variable
			int w = 1; // function variable
			double bnd = 1; // boundary function
			double _x = df2/x;


			cdfchi(&w,&p,&q,&_x,&df2,&st,&bnd);
			return lower_tail?q:p;
			//return pchisq(df2 / x , df2, !lower_tail, log_p);
		}

		/* Avoid squeezing pbeta's first parameter against 1 :  */
		if (df1 * x > df2)
			x = betaP(df2 / (df2 + df1 * x), df2 / 2., df1 / 2., 
				!lower_tail, 0);
		else
			x = betaP(df1 * x / (df2 + df1 * x), df1 / 2., df2 / 2.,
				lower_tail, 0);

		return x==x ? x : numeric_limits<double>::quiet_NaN();
	}

	double pT(double t, double df)
	{
		int which = 1, stat = 0;
		double p = WISARD_NAN;
		double q = WISARD_NAN;
		double bound = WISARD_NAN;
		cdft(&which, &p, &q, &t, &df, &stat, &bound);

		return p;
	}

	double dbeta(double x, double a, double b, int log_p)
	{
		if (a <= 0 || b <= 0) return numeric_limits<double>::quiet_NaN();
		if (x < 0 || x > 1) return(R_D__0);
		if (x == 0) {
			if(a > 1) return(R_D__0);
			if(a < 1) return(ML_POSINF);
			/* a == 1 : */ return(R_D_val(b));
		}
		if (x == 1) {
			if(b > 1) return(R_D__0);
			if(b < 1) return(ML_POSINF);
			/* b == 1 : */ return(R_D_val(a));
		}
		double lval;
		if (a <= 2 || b <= 2)
			lval = (a-1)*log(x) + (b-1)*util_log1p(-x) - lbeta(a, b);
		else if(!_finite(a) || !_finite(b))
			lval = ML_NEGINF;
		else
			lval = log(a+b-1) + dbinom_raw(a-1, a+b-2, x, 1-x, TRUE);

		return R_D_exp2(lval);
	}

	double pdhyper (double x, double NR, double NB, double n, int log_p)
	{
		/*
		* Calculate
		*
		*	    phyper (x, NR, NB, n, TRUE, FALSE)
		*   [log]  ----------------------------------
		*	       dhyper (x, NR, NB, n, FALSE)
		*
		* without actually calling phyper.  This assumes that
		*
		*     x * (NR + NB) <= n * NR
		*
		*/
		double sum = 0;
		double term = 1;

		while (x > 0 && term >= DBL_EPSILON * sum) {
			term *= x * (NB - n + x) / (n + 1 - x) / (NR + 1 - x);
			sum += term;
			x--;
		}

		double ss = (double) sum;
		return log_p ? util_log1p(ss) : 1 + ss;
	}


	double dhyper(double x, double r, double b, double n, int give_log)
	{
		double p, q, p1, p2, p3;

		if (R_D_negInonint(r) || R_D_negInonint(b) || R_D_negInonint(n) || n > r+b)
			return WISARD_NAN;
		if(x < 0) return W0;
		R_D_nonint_check(x);// incl warning

		if (n < x || r < x || n - x > b) return(W0);
		if (n == 0) return((x == 0) ? W1 : W0);

		p = ((double)n)/((double)(r+b));
		q = ((double)(r+b-n))/((double)(r+b));

		p1 = dbinom_raw(x,	r, p,q,give_log);
		p2 = dbinom_raw(n-x,b, p,q,give_log);
		p3 = dbinom_raw(n,r+b, p,q,give_log);

		return( (give_log) ? p1 + p2 - p3 : p1*p2/p3 );
	}

#define R_DT_Log(p)	(lower_tail? (p) : R_Log1_Exp(p))

	/* FIXME: The old phyper() code was basically used in ./qhyper.c as well
	* -----  We need to sync this again!
	*/
	double phyper (wsUint _x, wsUint NR, wsUint NB, wsUint n,
		int lower_tail, int log_p)
	{
		/* Sample of  n balls from  NR red  and	 NB black ones;	 x are red */
		double x = (double)_x;
		double d, pd;

		x = floor (x + 1e-7);

		if ((NR + NB) == WISARD_INF || (n > NR + NB))
			return WISARD_NAN;

		if (x * (NR + NB) > n * NR) {
			/* Swap tails.	*/
			wsUint oldNB = NB;
			NB = NR;
			NR = oldNB;
			x = n - x - 1;
			lower_tail = !lower_tail;
		}

		if (x < 0) return W0;
		if (x >= NR || x >= n) return W1;

		d  = dhyper (x, NR, NB, n, log_p);
		pd = pdhyper(x, NR, NB, n, log_p);

		return log_p ? R_DT_Log(d + pd) : R_D_Lval(d * pd);
	}

} // End namespace ONETOOL
