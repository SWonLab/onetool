#define SSEINLINE

/*
 * Global include
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <setjmp.h>
#include <limits>
#ifdef _WIN32
#	include <windows.h>
#	include <tchar.h>
#else
#	include <unistd.h>
#	include <limits.h>
#endif

/*
 * Local include
 */
#include "utils/dcdflib.h"
#include "global/option.h"
#include "global/common.h"
#include "utils/util.h"
#include "global/worker.h"
#include "global/io.h"
#include "utils/data.h"
#include "utils/matrix.h"
#include "input/stream.h"
#include "utils/stat.h"
#include "utils/new_stat.h"

/*
 * BSD-specific dependency for system function access
 */
#if defined(__FreeBSD__) || defined(__APPLE__)
#	ifdef __FreeBSD__
#		include <kvm.h>
#	else
#		include <sys/ptrace.h>
#	endif
#	include <sys/types.h>
#	include <sys/sysctl.h>
#	include <sys/user.h>
#endif

//#define CHECK_SSEIGEN
//#define USE_SSEIGEN
#define CHECK_SSEIGEN_THR 1e-12

namespace ONETOOL {

/************************************************************************
FIFDINT:
Truncates a double precision number to an integer and returns the
value in a double.4
************************************************************************/
inline double fifdint(double a)
/* a     -     number to be truncated */
{
	return (double) ((int) a);
}

// template<typename T> inline T SQR(T a)
// {
// 	return a*a;
// }

double pythag(const double a, const double b)
{
	double absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

/*
IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
THAT IS USED. 
*/
int ipmpar(int *i)
{
	static int imach[11];
	static int ipmpar;
	imach[1] = 2;
	imach[2] = 31;
	imach[3] = 2147483647;
	imach[4] = 2;
	imach[5] = 24;
	imach[6] = -125;
	imach[7] = 128;
	imach[8] = 53;
	imach[9] = -1021;
	imach[10] = 1024;
	ipmpar = imach[*i];
	return ipmpar;
}
inline int ipmpar_v2(const int i)
{
	if (i < 1 || i > 10) halt("Invalid range");
	static const int imach[11] = { 0, 2, 31, 2147483647, 2, 24, -125,
		128, 53, -1021, 1024 };
	return imach[i];
}
/*
SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
THE COMPUTER BEING USED.
*/
double spmpar(int i)
{
	static int K1 = 4;
	static int K2 = 8;
	static int K3 = 9;
	static int K4 = 10;
	static double spmpar,b,binv,bm1,one,w,z;
	static int emax,emin,ibeta,m;
	if(i > 1) goto S10;
	b = ipmpar(&K1);
	m = ipmpar(&K2);
	spmpar = pow(b,(double)(1-m));
	return spmpar;
S10:
	if(i > 2) goto S20;
	b = ipmpar(&K1);
	emin = ipmpar(&K3);
	one = 1.0;
	binv = one/b;
	w = pow(b,(double)(emin+2));
	spmpar = w*binv*binv*binv;
	return spmpar;
S20:
	ibeta = ipmpar(&K1);
	m = ipmpar(&K2);
	emax = ipmpar(&K4);
	b = ibeta;
	bm1 = ibeta-1;
	one = 1.0;
	z = pow(b,(double)(m-1));
	w = ((z-one)*b+bm1)/(b*z);
	z = pow(b,(double)(emax-2));
	spmpar = w*z*b*b;
	return spmpar;
}

void cumnor(double *arg,double *result,double *ccum)
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
	eps = spmpar(K1)*0.5e0;
	min = spmpar(K2);
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
		xsq = fifdint(y*sixten)/sixten;
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
		xsq = fifdint(x*sixten)/sixten;
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

void float_fmt(wsReal v, char *p)
{
	int e = 0;
	wsReal x = 1.0;
	if (v == 0.0)
		sprintf(p, "0.0");
	else {
		while (v < x) {
			e++;
			v *= 10;
		}
		sprintf(p, "%.2lfe-%d", v, e);
	}
}

/* LOG2���� ���� */
double LOG2(double R_number)
{
	const double R_log_e_2 = log(2.0);

	// log(n)/log(2) is log2
	return log(R_number)  / R_log_e_2;
}

/* Bonferroni correction
 *
 * -------------+----------------------------------------------------
 * wsReal	R_val	| A value to apply Bonferroni correction
 * int	N_count	| # of multiple cases
 * -------------+----------------------------------------------------
 *  return		| Corrected p-value, return 1 if saturated
 */
wsRealCst bonferroni(wsReal R_val, int N_count)
{
	wsReal R_adjP = R_val * N_count;
	return R_adjP > W1 ? W1 : R_adjP;
} /* END OF bonferroni() */

/* Launch R with specific R script path
 *
 * ------------------+----------------------------------------------------
 * char	*S_Rfn		 | A file path of R script to be executed
 * ------------------+---------------------------------------------------------
 *  return			 | Return code from R
 */
int		launchR(const char *S_Rfn)
{
	int		N_res;
	char	S_buf[512];

	sprintf(S_buf, "%s\\Rscript %s", OPT_STRING(rpath), S_Rfn);
#ifdef _WIN32
	STARTUPINFO si;
	::ZeroMemory (&si, sizeof (STARTUPINFO));
	si.cb = sizeof (STARTUPINFO);
	PROCESS_INFORMATION pi;

	if (::CreateProcess(NULL, S_buf, NULL,
		NULL, FALSE, NORMAL_PRIORITY_CLASS, NULL, NULL, &si, &pi)) {
		DWORD N_retCode;
		::CloseHandle (pi.hThread);
		::WaitForSingleObject (pi.hProcess, INFINITE);
		::CloseHandle (pi.hProcess);
		GetExitCodeProcess(pi.hProcess, &N_retCode);

		N_res = (int)N_retCode;
	} else
		halt_fmt(WISARD_FAIL_LAUNCH_R, S_buf, GetLastError());
#else
	FILE *fp = popen(S_buf, "r");
	if (fp == NULL) {
		perror("Failed to launch");
		halt("See above error message");
	}
	while (fgets(S_buf, 512, fp));
	N_res = (int)pclose(fp);
#endif
	if (N_res != 0)
		halt("An error from R [Return code : %d]", N_res);

	return N_res;
}

double PVchisq(double x, double df, double *nc/*=NULL*/, char B_quiet/*=0*/)
{
//	if ( ! realnum(x) ) return -
 	if (x != x || df != df) {
//		if (B_quiet == 0)
//			LOG("Chi-distribution probability is invalid\n");
		return numeric_limits<double>::quiet_NaN();
	}

	double p, q;
	int st = 0; // error variable
	int w = 1; // function variable
	double bnd = 1; // boundary function

	if (x == numeric_limits<double>::infinity())
		return 0.0;

	// NCP is set to 0
	if (nc)
		cdfchn(&w, &p, &q, &x, &df, nc, &st, &bnd);
	else
		cdfchi(&w,&p,&q,&x,&df,&st,&bnd);

	// Check status
	if (st != 0)
		return 1.0;

	// Return p-value
	return q;
}

double norprobP(double x, double mean, double sd)
{
	double p, q;
	double bnd = 1;
	int st = 0;
	x_cdfnor(1, &p, &q, &x, &mean, &sd, &st, &bnd);

	if (st != 0)
		return 1.0;

	return p;
}

double chebyshev_eval(double x, const double *a, const int n)
{
	double b0, b1, b2, twox;
	int i;

	if (n < 1 || n > 1000) return numeric_limits<double>::quiet_NaN();

	if (x < -1.1 || x > 1.1) return numeric_limits<double>::quiet_NaN();

	twox = x * 2;
	b2 = b1 = 0;
	b0 = 0;
	for (i = 1; i <= n; i++) {
		b2 = b1;
		b1 = b0;
		b0 = twox * b1 - b2 + a[n - i];
	}
	return (b0 - b2) * 0.5;
}

double util_log1p(double x)
{
    /* series for log1p on the interval -.375 to .375
     *				     with weighted error   6.35e-32
     *				      log weighted error  31.20
     *			    significant figures required  30.93
     *				 decimal places required  32.01
     */
    const static double alnrcs[43] = {
	+.10378693562743769800686267719098e+1,
	-.13364301504908918098766041553133e+0,
	+.19408249135520563357926199374750e-1,
	-.30107551127535777690376537776592e-2,
	+.48694614797154850090456366509137e-3,
	-.81054881893175356066809943008622e-4,
	+.13778847799559524782938251496059e-4,
	-.23802210894358970251369992914935e-5,
	+.41640416213865183476391859901989e-6,
	-.73595828378075994984266837031998e-7,
	+.13117611876241674949152294345011e-7,
	-.23546709317742425136696092330175e-8,
	+.42522773276034997775638052962567e-9,
	-.77190894134840796826108107493300e-10,
	+.14075746481359069909215356472191e-10,
	-.25769072058024680627537078627584e-11,
	+.47342406666294421849154395005938e-12,
	-.87249012674742641745301263292675e-13,
	+.16124614902740551465739833119115e-13,
	-.29875652015665773006710792416815e-14,
	+.55480701209082887983041321697279e-15,
	-.10324619158271569595141333961932e-15,
	+.19250239203049851177878503244868e-16,
	-.35955073465265150011189707844266e-17,
	+.67264542537876857892194574226773e-18,
	-.12602624168735219252082425637546e-18,
	+.23644884408606210044916158955519e-19,
	-.44419377050807936898878389179733e-20,
	+.83546594464034259016241293994666e-21,
	-.15731559416479562574899253521066e-21,
	+.29653128740247422686154369706666e-22,
	-.55949583481815947292156013226666e-23,
	+.10566354268835681048187284138666e-23,
	-.19972483680670204548314999466666e-24,
	+.37782977818839361421049855999999e-25,
	-.71531586889081740345038165333333e-26,
	+.13552488463674213646502024533333e-26,
	-.25694673048487567430079829333333e-27,
	+.48747756066216949076459519999999e-28,
	-.92542112530849715321132373333333e-29,
	+.17578597841760239233269760000000e-29,
	-.33410026677731010351377066666666e-30,
	+.63533936180236187354180266666666e-31,
    };


# define nlnrel 22
//    const static double xmin = -0.999999985;
/* 22: for IEEE double precision where DBL_EPSILON =  2.22044604925031e-16 */

    if (x == 0.) return 0.;/* speed */
    if (x == -1) return(ML_NEGINF);
    if (x  < -1) return numeric_limits<double>::quiet_NaN();

    if (fabs(x) <= .375) {
        /* Improve on speed (only);
	   again give result accurate to IEEE double precision: */
	if(fabs(x) < .5 * DBL_EPSILON)
	    return x;

	if( (0 < x && x < 1e-8) || (-1e-9 < x && x < 0))
	    return x * (1 - .5 * x);
	/* else */
	return x * (1 - x * chebyshev_eval(x / .375, alnrcs, nlnrel));
    }
    /* else */
//    if (x < xmin) {
	/* answer less than half precision because x too near -1 */
//	ML_ERROR(ME_PRECISION, "log1p");
//    }
    return log(1 + x);
}

double qnorm(double p, double mu, double sigma, int lower_tail, int log_p)
{
    double p_, q, r, val;

#ifdef IEEE_754
    if (ISNAN(p) || ISNAN(mu) || ISNAN(sigma))
	return p + mu + sigma;
#endif
    R_Q_P01_boundaries(p, ML_NEGINF, ML_POSINF);

    if(sigma  < 0)	return numeric_limits<double>::quiet_NaN();
    if(sigma == 0)	return mu;

    p_ = R_DT_qIv(p);/* real lower_tail prob. p */
    q = p_ - 0.5;

#ifdef DEBUG_qnorm
    REprintf("qnorm(p=%10.7g, m=%g, s=%g, l.t.= %d, log= %d): q = %g\n",
	     p,mu,sigma, lower_tail, log_p, q);
#endif


/*-- use AS 241 --- */
/* double ppnd16_(double *p, long *ifault)*/
/*      ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3

        Produces the normal deviate Z corresponding to a given lower
        tail area of P; Z is accurate to about 1 part in 10**16.

        (original fortran code used PARAMETER(..) for the coefficients
         and provided hash codes for checking them...)
*/
    if (fabs(q) <= .425) {/* 0.075 <= p <= 0.925 */
        r = .180625 - q * q;
	val =
            q * (((((((r * 2509.0809287301226727 +
                       33430.575583588128105) * r + 67265.770927008700853) * r +
                     45921.953931549871457) * r + 13731.693765509461125) * r +
                   1971.5909503065514427) * r + 133.14166789178437745) * r +
                 3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                     28729.085735721942674) * r + 39307.89580009271061) * r +
                   21213.794301586595867) * r + 5394.1960214247511077) * r +
                 687.1870074920579083) * r + 42.313330701600911252) * r + 1.);
    }
    else { /* closer than 0.075 from {0,1} boundary */

	/* r = min(p, 1-p) < 0.075 */
	if (q > 0)
	    r = R_DT_CIv(p);/* 1-p */
	else
	    r = p_;/* = R_DT_Iv(p) ^=  p */

	r = sqrt(- ((log_p &&
		     ((lower_tail && q <= 0) || (!lower_tail && q > 0))) ?
		    p : /* else */ log(r)));
        /* r = sqrt(-log(r))  <==>  min(p, 1-p) = exp( - r^2 ) */
#ifdef DEBUG_qnorm
	REprintf("\t close to 0 or 1: r = %7g\n", r);
#endif

        if (r <= 5.) { /* <==> min(p,1-p) >= exp(-25) ~= 1.3888e-11 */
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                       .0227238449892691845833) * r + .24178072517745061177) *
                     r + 1.27045825245236838258) * r +
                    3.64784832476320460504) * r + 5.7694972214606914055) *
                  r + 4.6303378461565452959) * r +
                 1.42343711074968357734)
                / (((((((r *
                         1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                        r + .0151986665636164571966) * r +
                       .14810397642748007459) * r + .68976733498510000455) *
                     r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1.);
        }
        else { /* very close to  0 or 1 */
            r += -5.;
            val = (((((((r * 2.01033439929228813265e-7 +
                       2.71155556874348757815e-5) * r +
                      .0012426609473880784386) * r + .026532189526576123093) *
                    r + .29656057182850489123) * r +
                   1.7848265399172913358) * r + 5.4637849111641143699) *
                 r + 6.6579046435011037772)
                / (((((((r *
                         2.04426310338993978564e-15 + 1.4215117583164458887e-7)*
                        r + 1.8463183175100546818e-5) * r +
                       7.868691311456132591e-4) * r + .0148753612908506148525)
                     * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1.);
        }

	if(q < 0.0)
	    val = -val;
        /* return (q >= 0.)? r : -r ;*/
    }
    return mu + sigma * val;
}

double util_expm1(double x)
{
	double y, a = fabs(x);

	if (a < DBL_EPSILON) return x;
	if (a > 0.697) return exp(x) - 1;  /* negligible cancellation */

	if (a > 1e-8)
		y = exp(x) - 1;
	else /* Taylor expansion, more accurate in this range */
		y = (x / 2 + 1) * x;

	/* Newton step for solving   log(1 + y) = x   for y : */
	/* WARNING: does not work for y ~ -1: bug in 1.5.0 */
	y -= (1 + y) * (util_log1p (y) - x);
	return y;
}

wsReal** sseVarCovRow(wsReal **Ra_mat, wsUint N_row, wsUint N_col)
{
	wsReal **Ra_ret = sseMatrix(N_row, N_row);
	wsUint i, j, k;

	/*** SHOULD copy Ra_ret[i][j] to Ra_ret[j][i] to make return matrix
	 *** SYMMETRIC!! i and j will cover the upper triangular part of
	 *** return matrix!! */

	for (i=0 ; i<N_row ; i++) {
		for (j=i ; j<N_row ; j++) {
#ifdef USE_SSE
			wsUint N_med = getMed(N_row);

			for (k=0 ; k<N_med ; k+=sseJmp) {
				/* SSE-enabled part */
			}
			for (k=N_med ; k<N_row ; k++) {
				/* Rest part */
			}
#else
			for (k=0 ; k<N_row ; k++) {
				/* Non-sse enabled full part */
			}
#endif
		}
	}

	return Ra_ret;
}

#ifdef USE_MEMDBG
wsReal** _sseEmptyMatrix(wsStrCst S_file, wsUint N_line, wsUint N_row, wsUint N_col)
#else
wsReal** _sseEmptyMatrix(wsUint N_row, wsUint N_col)
#endif
{
	wsReal **Ra_ret = NULL;
	MULTI_MALLOC_SYM(Ra_ret, wsReal*, N_row);
	for (wsUint i=0 ; i<N_row ; i++) {
		sseCallocSym(Ra_ret[i], wsReal, N_col);
	}

	return Ra_ret;
}

#ifdef USE_MEMDBG
wsReal** _sseEmptySymMat(wsStrCst S_file, wsUint N_line, wsUint N_sz)
#else
wsReal** _sseEmptySymMat(wsUint N_sz)
#endif
{
	wsReal **Ra_ret = NULL;
	MULTI_MALLOC_SYM(Ra_ret, wsReal*, N_sz);
	for (wsUint i=0 ; i<N_sz ; i++) {
#ifdef USE_SYM
		sseCallocSym(Ra_ret[i], wsReal, (i+1));
#else
		sseCallocSym(Ra_ret[i], wsReal, N_sz);
#endif
	}

	return Ra_ret;
}

wsRealCst sseVpPpV(wsReal *Ra_V, wsUint N_sz, wsReal *Ra_p, wsUint N_c,
	wsUint *Na_rEnd, wsUint *Na_cIdx, wsReal *Ra_W/*=NULL*/,
	wsReal *Ra_mask/*=NULL*/)
{
	if (Ra_W == NULL) Ra_W = Ra_V;

	wsReal *Ra_vM = sseVpP(Ra_V, N_sz, Ra_p, N_c, Na_rEnd, Na_cIdx, Ra_mask);
	wsReal R;
	if (Ra_mask != NULL)
		R = sseVVsel(Ra_mask, Ra_vM, N_sz, Ra_W);
	else
		R = sseVV(Ra_vM, N_sz, Ra_W);
	sseFree(Ra_vM);

	return R;
}

void LRtest(wsReal R_v1, wsReal R_v2, wsReal *Rp_lrt, wsReal *Rp_pValue,
	wsReal R_df/*=W1*/)
{
	wsReal R_lrt = REAL_CONST(-2.0)*(R_v1-R_v2);
	wsReal R_pValue = (wsReal)PVchisq(R_lrt, R_df);

	if (Rp_lrt) *Rp_lrt = R_lrt;
	if (Rp_pValue) *Rp_pValue = R_pValue;
}

#define SIGN(Ra_mat,b)	((b)>=0?((Ra_mat)>=0?(Ra_mat):-(Ra_mat)):((Ra_mat)>=0)?-(Ra_mat):(Ra_mat))
#ifndef MAX
#	define MAX(Ra_mat,b)	((b)>(Ra_mat)?(b):(Ra_mat))
#	define MIN(Ra_mat,b)	((b)<(Ra_mat)?(b):(Ra_mat))
#endif

inline wsReal SVDpythag(const wsReal a, const wsReal b)
{
	double absa,absb;

	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return (wsReal)(absa*sqrt(1.0+SQR(absb/absa)));
	else return (absb == 0.0 ? W0 : (wsReal)(absb*sqrt(1.0+SQR(absa/absb))));
}

/* 
 * svdcomp - SVD decomposition routine. 
 * Takes an mxn matrix a and decomposes it into udv, where u,v are
 * left and right orthogonal transformation matrices, and d is a 
 * diagonal matrix of singular values.
 *
 * This routine is adapted from svdecomp.c in XLISP-STAT 2.1 which is 
 * code from Numerical Recipes adapted by Luke Tierney and David Betz.
 *
 * Input to dsvd is as follows:
 *   a = mxn matrix to be decomposed, gets overwritten with u
 *   m = row dimension of a
 *   n = column dimension of a
 *   w = returns the vector of singular values of a
 *   v = returns the right orthogonal transformation matrix
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
 
int dsvd(wsVec a, int m, wsReal* w)
{
    int its, j, l;
    double f, s;
    double anorm = 0.0, g = 0.0, scale = 0.0;
    double *rv1;
  
    if (m < 1) 
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }
  
	rv1 = (double *)malloc((unsigned int)1*sizeof(double));

	/* Householder reduction to bidiagonal form */ {
		/* left-hand reduction */
		l = 1;
		rv1[0] = scale * g;
		g = s = scale = 0.0;

		for (int k = 0; k < m; k++)
			scale += fabs((double)a[k]);
		if (scale)
		{
			for (int k = 0; k < m; k++)
			{
				a[k] = (float)((double)a[k]/scale);
				s += ((double)a[k] * (double)a[k]);
			}
			f = (double)a[0];
			g = -SIGN(sqrt(s), f);
			a[0] = (float)(f - g);
			for (int k = 0; k < m; k++)
				a[k] = (float)((double)a[k]*scale);
		}
		w[0] = (float)(scale * g);

		/* right-hand reduction */
		g = s = scale = 0.0;
		anorm = MAX(anorm, (fabs((double)w[0]) + fabs(rv1[0])));
	}
  
	/* accumulate the right-hand transformation */ {
		wsUint i = 0;
		g = rv1[i];
		l = i;
	}
  
	/* accumulate the left-hand transformation */ {
		l = 1;
		g = (double)w[0];
		if (g) {
			g = 1.0 / g;
			for (j = 0; j < m; j++)
				a[j] = (float)((double)a[j]*g);
		} else {
			for (j = 0; j < m; j++)
				a[j] = 0.0;
		}
		++a[0];
	}

	/* diagonalize the bidiagonal form */ {
		/* loop over singular values */
		for (its = 0; its < 30; its++)
		{                         /* loop over allowed iterations */
			l = 0; /* test for splitting */
			/* Always true */
			if (fabs(rv1[l]) + anorm == anorm)
			{
				break;
			}
		}
	}
    free((void*) rv1);
    return(1);
}

#define sign3(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
bool SingularValueDecomposition(wsReal **a, int n, wsReal *w, wsReal **v)
{
	bool flag;
	int i,its,j,jj,k,l=0,nm=0;
	double anorm,c,f,g,h,s,scale,x,y,z;
	double volatile temp;

	vector<double> rv1(n);
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < n) {
			for (k=i;k<n;k++) scale += fabs(a[k][i]);
			if (scale != 0.0) {
				for (k=i;k<n;k++) {
					a[k][i] /= (wsReal)scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=(wsReal)(f-g);
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<n;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<n;k++) a[k][j] += (wsReal)f*a[k][i];
				}
				for (k=i;k<n;k++) a[k][i] *= (wsReal)scale;
			}
		}
		w[i]=(wsReal)(scale *g);
		g=s=scale=0.0;
		if (i+1 <= n && i+1 != n) {
			for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					a[i][k] /= (wsReal)scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l-1]=(wsReal)(f-g);
				for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l-1;k<n;k++) a[j][k] += (wsReal)(s*rv1[k]);
				}
				for (k=l-1;k<n;k++) a[i][k] *= (wsReal)scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++)
					v[j][i]=(a[i][j]/a[i][l])/(wsReal)g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += (wsReal)s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(n,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<n;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<n;k++) a[k][j] += (wsReal)f*a[k][i];
			}
			for (j=i;j<n;j++) a[j][i] *= (wsReal)g;
		} else for (j=i;j<n;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
				nm=l-1;
				temp=fabs(rv1[l])+anorm;
				if (temp == anorm) {
					flag=false;
					break;
				}
				temp=fabs(w[nm])+anorm;
				if (temp == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					temp = fabs(f)+anorm;
					if (temp == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=(wsReal)h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<n;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=(wsReal)(y*c+z*s);
						a[j][i]=(wsReal)(z*c-y*s);
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = (wsReal)(-z);
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) 
				return false; // cannot converge: multi-collinearity?
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=(wsReal)(x*c+z*s);
					v[jj][i]=(wsReal)(z*c-x*s);
				}
				z=pythag(f,h);
				w[j]=(wsReal)z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<n;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=(wsReal)(y*c+z*s);
					a[jj][i]=(wsReal)(z*c-y*s);
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=(wsReal)x;
		}
	}
	return true;
}

bool SVDcmp(wsReal **a/*RR*/, wsReal *w/*R*/, wsReal **v/*RR*/, int N)
{
	bool flag;
	int i,its,j,jj,k,l=0,nm=0;
	wsReal anorm,c,f,g,h,s,scale,x,y,z;
	wsReal volatile temp;

	int m=N;
	if (m == 0)
		halt("Internal problem in SVD function (no observations left?)");
	int n=N;

	vector<wsReal> rv1(n);
	g=scale=anorm=W0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=W0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a[k][i]);
			if (scale != W0) {
				for (k=i;k<m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=W0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
			if (scale != W0) {
				for (k=l-1;k<n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l-1];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l-1]=f-g;
				for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l-1;k<n;k++) a[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) a[i][k] *= scale;
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != W0) {
				for (j=l;j<n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g != W0) {
			g=W1/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += f*a[k][i];
			}
			for (j=i;j<m;j++) a[j][i] *= g;
		} else for (j=i;j<m;j++) a[j][i]=W0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
				nm=l-1;
				temp=fabs(rv1[l])+anorm;
				if (temp == anorm) {
					flag=false;
					break;
				}
				temp=fabs(w[nm])+anorm;
				if (temp == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					temp = fabs(f)+anorm;
					if (temp == anorm) break;
					g=w[i];
					h=(wsReal)pythag(f,g);
					w[i]=h;
					h=W1/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < W0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) 
				return false; // cannot converge: multi-collinearity?
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(W2*h*y);
			g=(wsReal)pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=(wsReal)pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=(wsReal)pythag(f,h);
				w[j]=z;
				if (z) {
					z=W1/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=W0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	return true;
}

bool SVDcmp2(wsReal **a/*RR*/, wsReal *w/*R*/, wsReal **v/*RR*/, int r)
{
	bool flag;
	int i,it,j,jj,k,l=0,nm=0;
	wsReal anorm,c,f,g=W0,h,s,scale,x,y,z;
	wsReal volatile temp;
	wsReal rv1[10];

	int m=r;
	int n=m;

	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		wsReal *pa = a[i];

		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a[k][i]);
			if (scale != 0.0) {
				for (k=i ; k<m ; k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=pa[i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				pa[i]=f-g;
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += f*a[k][i];
				}
				for (k=i;k<m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale*g;
		g=s=scale=0.0;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1 ; k<n ; k++) scale += fabs(pa[k]);
			if (scale != 0.0) {
				for (k=l-1 ; k<n ; k++) {
					a[i][k]	/= scale;
					s		+= SQR(pa[k]);
				}
				f		= pa[l-1];
				g		= -SIGN(sqrt(s),f);
				h		= f*g-s;
				pa[l-1]	= f-g;
				for (k=l-1 ; k<n ; k++) rv1[k]=pa[k]/h;
				for (j=l-1 ; j<m ; j++) {
					for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*pa[k];
					for (k=l-1 ; k<n ; k++) a[j][k] += s*rv1[k];
				}
				for (k=l-1;k<n;k++) pa[k] *= scale;
			}
		}
		anorm = MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		wsReal *pa = a[i];
		if (i < n-1) {
			if (g != 0.0) {
				for	(j=l ; j<n ; j++)
					v[j][i] = (pa[j]/pa[l])/g;
				for (j=l ; j<n ; j++) {
					for (s=0.0,k=l ; k<n ; k++) s += pa[k]*v[k][j];
					for (k=l ; k<n ; k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l ; j<n ; j++) v[i][j] = v[j][i] = 0.0;
		}
		v[i][i]	= 1.0;
		g		= rv1[i];
		l		= i;
	}
	for (i=MIN(m,n)-1 ; i>=0 ; i--) {
		wsReal *pa = a[i];
		l = i+1;
		g = w[i];
		for (j=l ; j<n ; j++) pa[j]=0.0;
		if (g != 0.0) {
			g=W1/g;
			for (j=l ; j<n ; j++) {
				for (s=0.0,k=l ; k<m ; k++) s += a[k][i]*a[k][j];
				f = (s/pa[i])*g;
				for (k=i ; k<m ; k++) a[k][j] += f*a[k][i];
			}
			for (j=i ; j<m ; j++) a[j][i] *= g;
		} else for (j=i  ;j<m ; j++) a[j][i] = 0.0;
		++pa[i];
	}
	for (k=n-1 ; k>=0 ; k--) {
		for (it=0 ; it<30 ; it++) {
			flag = true;
			for (l=k ; l>=0 ; l--) {
				nm		= l-1;
				temp	= fabs(rv1[l])+anorm;
				if (temp == anorm) {
					flag=false;
					break;
				}
				temp	= fabs(w[nm])+anorm;
				if (temp == anorm) break;
			}
			if (flag) {
				c = 0.0;
				s = 1.0;
				for (i=l ; i<k+1 ; i++) {
					f		= s*rv1[i];
					rv1[i]	= c*rv1[i];
					temp	= fabs(f)+anorm;
					if (temp == anorm) break;
					g		= w[i];
					h		= SVDpythag(f,g);
					w[i]	= h;
					h		= W1/h;
					c		= g*h;
					s		= -f*h;
					for (j=0 ; j<m ; j++) {
						wsReal *pa = a[j];
						y		= pa[nm];
						z		= pa[i];
						pa[nm]	= y*c+z*s;
						pa[i]	= z*c-y*s;
					}
				}
			}
			z = w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (it == 29) 
				return false; // cannot converge: multi-collinearity?
			x	= w[l];
			nm	= k-1;
			y	= w[nm];
			g	= rv1[nm];
			h	= rv1[k];
			f	= ((y-z)*(y+z)+(g-h)*(g+h))/(W2*h*y);
			g	= SVDpythag(f,1.0);
			f	= ((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c = s = 1.0;
			for (j=l ; j<=nm ; j++) {
				i		= j+1;
				g		= rv1[i];
				y		= w[i];
				h		= s*g;
				g		= c*g;
				z		= SVDpythag(f,h);
				rv1[j]	= z;
				c		= f/z;
				s		= h/z;
				f		= x*c+g*s;
				g		= g*c-x*s;
				h		= y*s;
				y		*= c;
				for (jj=0 ; jj<n ; jj++) {
					wsReal *pv=v[jj];
					x		= pv[j];
					z		= pv[i];
					pv[j]	= x*c+z*s;
					pv[i]	= z*c-x*s;
				}
				z		= SVDpythag(f,h);
				w[j]	= z;
				if (z) {
					z = W1/z;
					c = f*z;
					s = h*z;
				}
				f = c*g+s*y;
				x = c*y-s*g;
				for (jj=0 ; jj<m ; jj++) {
					wsReal *pa=a[jj];
					y		= pa[j];
					z		= pa[i];
					pa[j]	= y*c+z*s;
					pa[i]	= z*c-y*s;
				}
			}
			rv1[l]	= 0.0;
			rv1[k]	= f;
			w[k]	= x;
		}
	}

	return true;
}

bool EDinverse(wsReal **Ra_matOrig, wsUint N_sz, wsReal ***Rp_res,
	wsUint *Np_rank/*=NULL*/)
{
	wsUint			i;//, j, k;
	const wsReal	R_eps		= numeric_limits<wsReal>::epsilon();
//	wsReal			**Ra_v		= sseEmptyMatrix(N_sz, N_sz); /* checked */
// 	wsReal			*Ra_wgt		= NULL; /* checked */
// 	sseCalloc(Ra_wgt, wsReal, N_sz);
//	wsReal			**Ra_mat	= sseMatrixP(N_sz, N_sz, Ra_matOrig); /* checked */
	wsMat			Ra_mat		= NULL; // P, not P^t
	char			B_failed	= 0;
	wsReal			*Ra_wgt		= EIGENDECOMPOSITION(Ra_matOrig, N_sz, &Ra_mat, 1);
//	bool			B_cvged		= eigenDecomp2(Ra_mat, N_sz, Ra_wgt, Ra_v);

	// Look for singular values
	wsReal wmax = W0;
	for (i=0; i<N_sz; i++)
		wmax = Ra_wgt[i] > wmax ? Ra_wgt[i] : wmax;
	wsReal wmin = wmax * R_eps;
	for (i=0; i<N_sz; i++)
		Ra_wgt[i] = Ra_wgt[i] < wmin ? 0 : 1/Ra_wgt[i];

	// Get rank of Ra_mat if assigned
	if (Np_rank) {
		*Np_rank = 0;
		for (i=0; i<N_sz; i++)
			if (Ra_wgt[i]) (*Np_rank)++;
	}

	// u w t(v)
	// row U * 1/w
	// results matrix
	wsReal **Ra_ret = sseMDMt(Ra_mat, N_sz, N_sz, Ra_wgt);
// 	for (i=0 ; i<N_sz ; i++)
// 		sseVpV(Ra_mat[i], Ra_wgt, Ra_mat[i], N_sz);
	sseFree(Ra_wgt);
	// 	for (i=0; i<N_sz; i++)
	// 		for (j=0; j<N_sz; j++)
	// 			Ra_mat[i][j] *= Ra_wgt[j];

	// [nxn].[t(v)] 
	//	printf("O matrix[%x]\n", O);
// 	wsReal **Ra_ret = sseMpMt(Ra_mat, N_sz, N_sz,
// 		Ra_v, N_sz, N_sz);
	// 	for (i=0 ; i<N_sz ; i++)
	// 		for (j=0 ; j<N_sz ; j++)
	// 			for (k=0 ; k<N_sz ; k++)
	// 				Ra_ret[i][j] += Ra_mat[i][k] * Ra_v[j][k];

	sseUnmat(Ra_mat, N_sz);
//	sseUnmat(Ra_v, N_sz);
	*Rp_res = Ra_ret;

	return !B_failed;
}

bool SVDinverse(wsReal **Ra_matOrig, wsUint N_sz, wsReal ***Rp_res,
				wsUint *Np_rank/*=NULL*/)
{
	wsUint			i;//, j, k;
	const wsReal	R_eps		= numeric_limits<wsReal>::epsilon();
	wsReal			**Ra_v		= sseEmptyMatrix(N_sz, N_sz); /* checked */
	wsReal			*Ra_wgt		= NULL; /* checked */
	sseCalloc(Ra_wgt, wsReal, N_sz);
	wsReal			**Ra_mat	= sseMatrixP(N_sz, N_sz, Ra_matOrig); /* checked */
	bool			B_cvged		= SingularValueDecomposition(Ra_mat, N_sz, Ra_wgt, Ra_v);

	// Look for singular values
	wsReal wmax = W0;
	for (i=0; i<N_sz; i++)
		wmax = Ra_wgt[i] > wmax ? Ra_wgt[i] : wmax;
	wsReal wmin = wmax * R_eps;
	for (i=0; i<N_sz; i++)
		Ra_wgt[i] = Ra_wgt[i] < wmin ? 0 : 1/Ra_wgt[i];

	// Get rank of Ra_mat if assigned
	if (Np_rank) {
		*Np_rank = 0;
		for (i=0; i<N_sz; i++)
			if (Ra_wgt[i]) (*Np_rank)++;
	}

	// u w t(v)
	// row U * 1/w
	// results matrix
	for (i=0 ; i<N_sz ; i++)
		sseVpV(Ra_mat[i], Ra_wgt, Ra_mat[i], N_sz);
	sseFree(Ra_wgt);
// 	for (i=0; i<N_sz; i++)
// 		for (j=0; j<N_sz; j++)
// 			Ra_mat[i][j] *= Ra_wgt[j];

	// [nxn].[t(v)] 
	//	printf("O matrix[%x]\n", O);
	wsReal **Ra_ret = sseMpMt(Ra_mat, N_sz, N_sz,
		Ra_v, N_sz, N_sz);
// 	for (i=0 ; i<N_sz ; i++)
// 		for (j=0 ; j<N_sz ; j++)
// 			for (k=0 ; k<N_sz ; k++)
// 				Ra_ret[i][j] += Ra_mat[i][k] * Ra_v[j][k];

	sseUnmat(Ra_mat, N_sz);
	sseUnmat(Ra_v, N_sz);
	*Rp_res = Ra_ret;
	return B_cvged;
}

wsReal** transpose(wsReal **Ra_mat, wsUint N_matR, wsUint N_matC)
{
	wsReal **Ra_ret = NULL;
	wsAlloc(Ra_ret, wsReal*, N_matC);
	for (wsUint i=0 ; i<N_matC ; i++) {
		sseMalloc(Ra_ret[i], wsReal, N_matR);
		for (wsUint j=0 ; j<N_matR ; j++)
			Ra_ret[i][j] = Ra_mat[j][i];
	}
	return Ra_ret;
}

void transposeSelf(wsReal **Ra_mat, wsUint N_matR, wsUint N_matC)
{
	if (N_matR != N_matC) halt("Can't do self-transpose[%d x %d]", N_matR, N_matC);
	wsUint N = N_matR;

	for (wsUint i=0 ; i<N ; i++)
		for (wsUint j=0 ; j<i ; j++) {
			wsReal R_tmp = Ra_mat[i][j];
			Ra_mat[i][j] = Ra_mat[j][i];
			Ra_mat[j][i] = R_tmp;
		}
}

wsReal** makeMatrix(const char *S_path, wsUint *Np_r, wsUint *Np_c)
{
	FILE	*H_fp = fopen(S_path, "r");
	char	*S_buf = NULL;
	wsUint	N_sz, N_c, N_r, i, j;
	if (H_fp == NULL)
		halt_fmt(WISARD_FAIL_OPEN_FILEWITHDESC, "source matrix", S_path);

	/* Get the size of columns */
	for (N_sz=0 ; !feof(H_fp) ; N_sz++) {
		int ch = fgetc(H_fp);
		if (ch == '\r' || ch == '\n') break;
	}
	rewind(H_fp);
	N_sz <<= 1;
	wsAlloc(S_buf, char, N_sz);

	/* Get the number of columns */
	if (fgets(S_buf, N_sz, H_fp) == NULL)
		halt_fmt(WISARD_NULL_FILEWITHDESC, "source matrix", S_path);
	char *a, *b;
	for (N_c=0,a=b=S_buf ; b ; N_c++,a=b)
		getString(&a, &b);

	/* Get the number of rows */
	for (N_r=1 ; fgets(S_buf, N_sz , H_fp) ; N_r++) {
		getString(&S_buf, &a);
		if (!S_buf[0]) break;
	}
	rewind(H_fp);

	LOG("[%d*%d matrix]\n", N_r, N_c);

	/* Allocate matrix */
	wsReal **Ra_ret = sseMatrix(N_r, N_c);
	for (i=0 ; fgets(S_buf, N_sz , H_fp) ; i++) {
		for (j=0,a=b=S_buf ; b ; j++,a=b) {
			getString(&a, &b);
			Ra_ret[i][j] = (wsReal)atof(a);
		}
		if (j != N_c)
			halt("Count does not match at %d line[%d get, %d Expected]",
			i+1, j, N_c);
	}

	if (Np_r) *Np_r = N_r;
	if (Np_c) *Np_c = N_c;

	DEALLOC(S_buf);
	return Ra_ret;
}

wsReal** makeMatrixSSE(const char *S_path, wsUint *Np_r, wsUint *Np_c)
{
	FILE	*H_fp = fopen(S_path, "r");
	char	*S_buf = NULL;
	wsUint	N_sz, N_c, N_r, i, j;
	if (H_fp == NULL)
		halt("Can't open matrix file [%s]", S_path);

	/* Get the size of columns */
	for (N_sz=0 ; !feof(H_fp) ; N_sz++) {
		int ch = fgetc(H_fp);
		if (ch == '\r' || ch == '\n') break;
	}
	rewind(H_fp);
	N_sz <<= 1;
	wsAlloc(S_buf, char, N_sz);

	/* Get the number of columns */
	if (fgets(S_buf, N_sz, H_fp) == NULL)
		halt("Failed to read matrix content, maybe this file is empty?");
	char *a, *b;
	for (N_c=0,a=b=S_buf ; a ; N_c++,a=b) {
		getString(&a, &b);
		if (b && !b[0])
			break;
	}

	/* Get the number of rows */
	for (N_r=1 ; fgets(S_buf, N_sz , H_fp) ; N_r++) {
		getString(&S_buf, &a);
		if (!S_buf[0]) {
			N_r-- ;
			break;
		}
	}
	rewind(H_fp);

	//LOG("[%d*%d matrix]\n", N_r, N_c);

	/* Allocate matrix */
	wsReal **Ra_ret = sseMatrix(N_r, N_c);
	for (i=0 ; fgets(S_buf, N_sz , H_fp) ; i++) {
		for (j=0,a=b=S_buf ; b ; j++,a=b) {
			getString(&a, &b);
			if (!strcmp("NA", a))
				Ra_ret[i][j] = WISARD_NA_REAL;
			else
				Ra_ret[i][j] = (wsReal)atof(a);
		}
		if (j != N_c)
			halt("Count does not match at %d line[%d get, %d Expected]",
			i+1, j, N_c);
	}

	if (Np_r) *Np_r = N_r;
	if (Np_c) *Np_c = N_c;

	DEALLOC(S_buf);
	return Ra_ret;
}

wsReal** makeSymMatSSE(const char *S_path, wsUint *Np_r)
{
	FILE	*H_fp = fopen(S_path, "r");
	char	*S_buf = NULL;
	wsUint	N_sz, N_c, N_r, i, j;
	if (H_fp == NULL)
		halt("Can't open matrix file [%s]", S_path);

	/* Get the size of columns */
	for (N_sz=0 ; !feof(H_fp) ; N_sz++) {
		int ch = fgetc(H_fp);
		if (ch == '\r' || ch == '\n') break;
	}
	rewind(H_fp);
	N_sz <<= 1;
	wsAlloc(S_buf, char, N_sz);

	/* Get the number of columns */
	if (fgets(S_buf, N_sz, H_fp) == NULL)
		halt("Failed to read matrix content, maybe this file is empty?");
	char *a, *b;
	for (N_c=0,a=b=S_buf ; a ; N_c++,a=b) {
		getString(&a, &b);
		if (b && !b[0])
			break;
	}

	/* Get the number of rows */
	for (N_r=1 ; fgets(S_buf, N_sz , H_fp) ; N_r++) {
		getString(&S_buf, &a);
		if (!S_buf[0]) {
			N_r-- ;
			break;
		}
	}
	rewind(H_fp);

	//LOG("[%d*%d matrix]\n", N_r, N_c);

	/* Allocate matrix */
	wsReal **Ra_ret = sseMatrix(N_r, N_c);
	for (i=0 ; fgets(S_buf, N_sz , H_fp) ; i++) {
		for (j=0,a=b=S_buf ; j<=i ; j++,a=b) {
			getString(&a, &b);
			if (!strcmp("NA", a))
				Ra_ret[i][j] = WISARD_NA_REAL;
			else
				Ra_ret[i][j] = (wsReal)atof(a);
		}
	}

	if (Np_r) *Np_r = N_r;

	DEALLOC(S_buf);
	return Ra_ret;
}

void compareMatrix(wsReal **Rp_mat, wsUint N_r, wsUint N_c, const char *S_path
	, const char *S_prefix)
{
	wsUint	N_sr, N_sc;
	wsReal	**Ra_mat = makeMatrix(S_path, &N_sr, &N_sc);

	if (N_sr != N_r) halt("Unmatch row [mem %d, file %d]", N_r, N_sr);
	if (N_sc != N_c) halt("Unmatch column [mem %d, file %d]", N_c, N_sc);

	if (S_prefix) {
		puts("");
		puts(S_prefix);
	}

	for (wsUint i=0 ; i<N_sr ; i++) {
		for (wsUint j=0 ; j<N_sc ; j++)
			printf("%g	", Rp_mat[i][j]-Ra_mat[i][j]);
		puts("");
	}
	deallocMatrix(Ra_mat, N_sr);
}

void compareMatrix(wsReal **Rp_mat1, wsUint N_r, wsUint N_c, wsReal **Rp_mat2, const char *S_path)
{
	double xxx = 0.0;
	FILE *fp = fopen(S_path, "w+");
	for (wsUint i=0 ; i<N_r ; i++) {
		for (wsUint j=0 ; j<N_c ; j++) {
			wsReal v = Rp_mat1[i][j]-Rp_mat2[i][j];
			if (v)
				fprintf(fp, "%g(%d,%d)	", v, i, j);
			xxx += v;
		}
		fprintf(fp, "\n");
	}
	printf("Error sum : %g\n", xxx);
	printf("Error mean : %g\n", xxx/(double)(N_r*N_c));
	fclose(fp);
}

char getEvEx(wsMat Ra_retEV, wsReal *realEigenvalues, wsReal *e, wsUint n,
	char B_negCorrect, char B_noEVX=0)
{
	//transposeSelf(Ra_retEV, n, n);
	wsMat z = Ra_retEV;

	/*
	realEigenvalues = new double[n];
	final double[] e = new double[n];
	for (int i = 0; i < n - 1; i++) {
		realEigenvalues[i] = main[i];
		e[i] = secondary[i];
	}
	realEigenvalues[n - 1] = main[n - 1];
	e[n - 1] = 0;*/

	// Determine the largest main and secondary value in absolute term.
	double maxAbsoluteValue = 0;
	for (wsUint i=0 ; i<n ; i++) {
		if (fabs(realEigenvalues[i]) > maxAbsoluteValue)
			maxAbsoluteValue = fabs(realEigenvalues[i]);
		if (fabs(e[i]) > maxAbsoluteValue)
			maxAbsoluteValue = fabs(e[i]);
	}
	// Make null any main and secondary value too small to be significant
	if (maxAbsoluteValue != 0) {
		for (wsUint i=0 ; i<n ; i++) {
			if (fabs(realEigenvalues[i]) <= DBL_EPSILON * maxAbsoluteValue) {
				realEigenvalues[i] = 0;
			}
			if (fabs(e[i]) <= DBL_EPSILON * maxAbsoluteValue) {
				e[i]=0;
			}
		}
	}

	for (wsUint j=0 ; j<n ; j++) {
		wsUint its = 0;
		wsUint m;
		do {
			for (m=j ; m<n-1 ; m++) {
				double delta = fabs(realEigenvalues[m]) +
					fabs(realEigenvalues[m + 1]);
				if (fabs(e[m]) + delta == delta) {
					break;
				}
			}
			if (m != j) {
				/* STOP */
				if (its == 200)
					return 1;

				its++;
				double q = (realEigenvalues[j + 1] - realEigenvalues[j]) / (2 * e[j]);
				double t = sqrt(1 + q * q);
				if (q < 0.0) {
					q = realEigenvalues[m] - realEigenvalues[j] + e[j] / (q - t);
				} else {
					q = realEigenvalues[m] - realEigenvalues[j] + e[j] / (q + t);
				}
				double u = 0.0;
				double s = 1.0;
				double c = 1.0;
				int i;
				for (i = (int)m-1 ; i>=(int)j ; i--) {
					double p = s * e[i];
					double h = c * e[i];
					if (fabs(p) >= fabs(q)) {
						c = q / p;
						t = sqrt(c * c + 1.0);
						e[i + 1] = p * t;
						s = 1.0 / t;
						c = c * s;
					} else {
						s = p / q;
						t = sqrt(s * s + 1.0);
						e[i + 1] = q * t;
						c = 1.0 / t;
						s = s * c;
					}
					if (e[i + 1] == 0.0) {
						realEigenvalues[i + 1] -= u;
						e[m] = 0.0;
						break;
					}
					q = realEigenvalues[i + 1] - u;
					t = (realEigenvalues[i] - q) * s + 2.0 * c * h;
					u = s * t;
					realEigenvalues[i + 1] = q + u;
					q = c * t - h;
#if 0
					for (UINT_t ia=0 ; ia<n ; ia++) {
						p = z[i + 1][ia];
						z[i + 1][ia] = s * z[i][ia] + c * p;
						z[i][ia] = c * z[i][ia] - s * p;
					}
#else

					if (!B_noEVX) {
						wsUint N_med = 0;
// #ifdef USE_SSE
// 						N_med = getMed(n);
// 						sse_t sse_s = sseSet(s);
// 						sse_t sse_c = sseSet(c);
// 						for (wsUint ia=0 ; ia<N_med ; ia+=sseJmp) {
// 							sse_t *sse_i0 = (sse_t *)(z[i] + ia);
// 							sse_t *sse_i1 = (sse_t *)(z[i+1] + ia);
// 							sse_t sse_v1 = *sse_i1;
// 
// 							*sse_i1 = sseAdd(sseMul(sse_s, *sse_i0), sseMul(sse_c, sse_v1));
// 							*sse_i0 = sseSub(sseMul(sse_c, *sse_i0), sseMul(sse_s, sse_v1));
// 	// 						p = z[i+1][ia];
// 	// 						z[i + 1][ia] = s * z[i][ia] + c * p;
// 	// 						z[i][ia] = c * z[i][ia] - s * p;
// 						}
//#endif
						for (wsUint ia=N_med ; ia<n ; ia++) {
							p = z[ia][i + 1];
							z[ia][i + 1] = s * z[ia][i] + c * p;
							z[ia][i] = c * z[ia][i] - s * p;

// 							p = z[i + 1][ia];
// 							z[i + 1][ia] = s * z[i][ia] + c * p;
// 							z[i][ia] = c * z[i][ia] - s * p;
						}
#endif
					}
				}
				if (t == 0.0 && i >= (int)j) {
					continue;
				}
				realEigenvalues[j] -= u;
				e[j] = q;
				e[m] = 0.0;
			}
		} while (m != j);
	}

	//Sort the eigen values (and vectors) in increase order
//	wsReal *Ra_p = sseVector(n);
// 	for (wsUint i=0 ; i<n ; i++) {
// 		wsUint k = i;
// 		double p = realEigenvalues[i];
// 		for (wsUint j=i+1 ; j<n ; j++) {
// 			if (realEigenvalues[j] > p) {
// 				k = j;
// 				p = realEigenvalues[j];
// 			}
// 		}
// 		if (!B_noEVX && k != i) {
//  			wsReal *Ra_p = z[i];
// // 			z[i] = z[k];
// // 			z[k] = Ra_p;
// 			realEigenvalues[k] = realEigenvalues[i];
// 			realEigenvalues[i] = p;
// 			memcpy(Ra_p, z[i], sizeof(wsReal)*n);
// 			memcpy(z[i], z[k], sizeof(wsReal)*n);
// 			memcpy(z[k], Ra_p, sizeof(wsReal)*n);
// // 			for (wsUint j=0 ; j<n ; j++) {
// // 				p = z[j][i];
// // 				z[j][i] = z[j][k];
// // 				z[j][k] = p;
// // 			}
// 		}
// 	}
//	sseFree(Ra_p);

	// Determine the largest eigen value in absolute term.
	maxAbsoluteValue = 0;
	for (wsUint i = 0; i < n; i++) {
		if (fabs(realEigenvalues[i]) > maxAbsoluteValue) {
			maxAbsoluteValue=fabs(realEigenvalues[i]);
		}
	}
	// Make null any eigen value too small to be significant
	if (maxAbsoluteValue!=0.0) {
		if (B_negCorrect) for (wsUint i=0; i < n; i++) {
			if (fabs(realEigenvalues[i]) < DBL_EPSILON * maxAbsoluteValue ||
				realEigenvalues[i] < W0) {
				realEigenvalues[i] = 0;
			}
		} else for (wsUint i=0; i < n; i++) {
			if (fabs(realEigenvalues[i]) < DBL_EPSILON * maxAbsoluteValue) {
				realEigenvalues[i] = 0;
			}
		}
	}
// 	eigenvectors = new ArrayRealVector[n];
// 	final double[] tmp = new double[n];
// 	for (int i = 0; i < n; i++) {
// 		for (int j = 0; j < n; j++) {
// 			tmp[j] = z[j][i];
// 		}
// 		eigenvectors[i] = new ArrayRealVector(tmp);
// 	}

	return 0;
}

void trdg(wsMat Ra_evx, wsUint N_sz, wsReal *Ra_ev, wsReal *Ra_e)
{
	wsReal	*z = sseVector(N_sz);

#if 0 /* DOES NOT WORK PROPERLY */
//	double	*main = new double[N_sz];
//	double	*secondary = new double[N_sz];
//	Ra_e[N_sz-1] = 0.0;
	Ra_e[0] = 0.0;

//	for (wsUint k=0 ; k<N_sz-1 ; k++) {
	for (int k=N_sz-1 ; k>=1 ; k--) {
		// zero-out a row and a column simultaneously
		double	*hK = Ra_evx[k];
		double	xNormSqr = 0;

		/* Set diagonal */
		Ra_ev[k] = hK[k];
//		for (wsUint j=k+1 ; j<N_sz ; j++)
		for (UINT_t j=0 ; j<k ; j++)
			xNormSqr += SQR(hK[j]);

//		wsReal *hKQ = hK + (k+1);
		REAL_t *hKQ = hK + (k-1);
		double a = (*hKQ > 0) ? -sqrt(xNormSqr) : sqrt(xNormSqr);
		Ra_e[k] = a;

		if (a != 0.0) {
			// apply Householder transform from left and right simultaneously
			*hKQ -= a;
			REAL_t	beta = -1 / (a * *hKQ);

			// compute a = beta A v, where v is the Householder vector
			// this loop is written in such a way
			//   1) only the upper triangular part of the matrix is accessed
			//   2) access is cache-friendly for a matrix stored in rows
//			memset(z+k+1, 0x00, sizeof(double)*(N_sz-k-1));
			memset(z, 0x00, sizeof(double)*k);
			//Arrays.fill(z, k + 1, N_sz, 0);
//			for (UINT_t i=k+1 ; i<N_sz ; ++i) {
			for (UINT_t i=0 ; i<k ; ++i) {
				REAL_t	*hI	= Ra_evx[i];
				REAL_t	hKI	= hK[i];
				double zI = hI[i] * hKI;
				for (UINT_t j=i+1 ; j<k ; ++j) {
//				for (UINT_t j=i+1 ; j<N_sz ; ++j) {
					double hIJ = hI[j];
					zI   += hIJ * hK[j];
					z[j] += hIJ * hKI;
				}
				z[i] = beta * (z[i] + zI);
			}

			// compute gamma = beta vT z / 2
			double gamma = 0;
//			for (UINT_t i=k+1 ; i<N_sz ; ++i) {
			for (UINT_t i=0 ; i<k ; ++i) {
				gamma += z[i] * hK[i];
			}
			gamma *= beta / 2;

			// compute z = z - gamma v
//			for (UINT_t i=k+1 ; i<N_sz ; ++i) {
			for (UINT_t i=0 ; i<k ; ++i) {
				z[i] -= gamma * hK[i];
			}

			// update matrix: A = A - v zT - z vT
			// only the upper triangular part of the matrix is updated
//			for (UINT_t i=k+1 ; i<N_sz ; ++i) {
			for (UINT_t i=0 ; i<k ; ++i) {
				REAL_t *hI = Ra_evx[i];
				for (UINT_t j=i ; j<k ; ++j)
//				for (UINT_t j=i ; j<N_sz ; ++j)
					hI[j] -= hK[i] * z[j] + z[i] * hK[j];
			}
		}
	}
//	Ra_ev[N_sz-1] = Ra_evx[N_sz-1][N_sz-1];
	Ra_ev[0] = Ra_evx[0][0];
#else
	//	double	*main = new double[N_sz];
	//	double	*secondary = new double[N_sz];
	Ra_e[N_sz-1] = 0.0;

	cTimer u;
	u.start();
	for (wsUint k=0 ; k<N_sz-1 ; k++) {
		if (u.get() > REAL_CONST(1000.0)) {
			notice("Tridiagonalization [dim %d]...%g%% over\r", N_sz,
				k/(wsReal)(N_sz-1)*REAL_CONST(100.0));
			u.start();
		}

		// zero-out a row and a column simultaneously
		wsReal	*hK = Ra_evx[k];
		wsReal	xNormSqr = 0;

		/* Set diagonal */
		Ra_ev[k] = hK[k];

// #ifdef USE_SSE
// 		/* Minimum starting point of SSE */
// 		wsUint N1 = k+1 + sseJmp-1;
// 		N1 = getMed(N1);
// 		/* Maximum reach point of SSE */
// 		wsUint N2 = getMed(N_sz);
// 		/* Use SSE method ONLY IF k<=N1<N2<=N_sz */
// 		if (N1 < N2) {
// 			for (wsUint j=k+1 ; j<N1 ; j++)
// 				xNormSqr += SQR(hK[j]);
// 			sse_t sse_sqr = sseSet(0.0);
// 			for (wsUint j=N1 ; j<N2 ; j+=sseJmp) {
// 				sse_t *sse_hK = (sse_t *)(hK + j);
// 				sse_sqr = sseAdd(sse_sqr, sseMul(*sse_hK, *sse_hK));
// 			}
// 			sseAppendSum(sse_sqr, xNormSqr);
// 			for (wsUint j=N2 ; j<N_sz ; j++)
// 				xNormSqr += SQR(hK[j]);
// 		} else
// #endif
			for (wsUint j=k+1 ; j<N_sz ; j++)
				xNormSqr += SQR(hK[j]);

		wsReal *hKQ = hK + (k+1);
		wsReal a = (*hKQ > 0) ? -sqrt(xNormSqr) : sqrt(xNormSqr);
		Ra_e[k] = a;

		if (a != 0.0) {
			// apply Householder transform from left and right simultaneously
			*hKQ -= a;
			wsReal	beta = -1 / (a * *hKQ);

			// compute a = beta A v, where v is the Householder vector
			// this loop is written in such a way
			//   1) only the upper triangular part of the matrix is accessed
			//   2) access is cache-friendly for a matrix stored in rows
			memset(z+k+1, 0x00, sizeof(double)*(N_sz-k-1));
			//Arrays.fill(z, k + 1, N_sz, 0);
			for (wsUint i=k+1 ; i<N_sz ; ++i) {
				wsReal	*hI	= Ra_evx[i];
				wsReal	hKI	= hK[i];
				wsReal zI = hI[i] * hKI;
// #ifdef USE_SSE
// 				N1 = i+1 + sseJmp-1;
// 				N1 = getMed(N1);
// 
// 				if (N1 < N2) {
// 					for (wsUint j=i+1 ; j<N1 ; ++j) {
// 						wsReal hIJ = hI[j];
// 						zI   += hIJ * hK[j];
// 						z[j] += hIJ * hKI;
// 					}
// 					sse_t sse_zI = sseSet(0.0);
// 					sse_t sse_hKI = sseSet(hKI);
// 					for (wsUint j=N1 ; j<N2 ; j+=sseJmp) {
// 						sse_t *sse_hI = (sse_t *)(hI + j);
// 						sse_t *sse_hK = (sse_t *)(hK + j);
// 						sse_t *sse_zJ = (sse_t *)(z + j);
// 
// 						sse_zI = sseAdd(sse_zI, sseMul(*sse_hI, *sse_hK));
// 						*sse_zJ = sseAdd(*sse_zJ, sseMul(*sse_hI, sse_hKI));
// 					}
// 					sseAppendSum(sse_zI, zI);
// 					for (wsUint j=N2 ; j<N_sz ; ++j) {
// 						wsReal hIJ = hI[j];
// 						zI   += hIJ * hK[j];
// 						z[j] += hIJ * hKI;
// 					}
// 				} else
// #endif
					for (wsUint j=i+1 ; j<N_sz ; ++j) {
						wsReal hIJ = hI[j];
						zI   += hIJ * hK[j];
						z[j] += hIJ * hKI;
					}
				z[i] = beta * (z[i] + zI);
			}

			// compute gamma = beta vT z / 2
			wsReal gamma = 0;
			for (wsUint i=k+1 ; i<N_sz ; ++i) {
				gamma += z[i] * hK[i];
			}
			gamma *= beta / 2;

			// compute z = z - gamma v
			for (wsUint i=k+1 ; i<N_sz ; ++i) {
				z[i] -= gamma * hK[i];
			}

			// update matrix: A = A - v zT - z vT
			// only the upper triangular part of the matrix is updated
			for (wsUint i=k+1 ; i<N_sz ; ++i) {
				wsReal *hI = Ra_evx[i];
// #ifdef USE_SSE
// 				N1 = i + sseJmp-1;
// 				N1 = getMed(N1);
// 				//notice("do %d\n", i);
// 				if (N1 < N2) {
// 					for (wsUint j=i ; j<N1 ; ++j)
// 						hI[j] -= hK[i] * z[j] + z[i] * hK[j];
// 					sse_t zI	= sseSet(z[i]);
// 					sse_t hkI	= sseSet(hK[i]);
// 					for (wsUint j=N1 ; j<N2 ; j+=sseJmp) {
// 						sse_t *sse_zJ = (sse_t *)(z + j);
// 						sse_t *sse_hkJ = (sse_t *)(hK + j);
// 						sse_t *sse_hiJ = (sse_t *)(hI + j);
// 
// 						*sse_hiJ = sseSub(*sse_hiJ, sseAdd(
// 							sseMul(hkI, *sse_zJ),
// 							sseMul(zI, *sse_hkJ)
// 						));
// 					}
// 					for (wsUint j=N2 ; j<N_sz ; ++j)
// 						hI[j] -= hK[i] * z[j] + z[i] * hK[j];
// 				} else
// #endif
					for (wsUint j=i ; j<N_sz ; ++j)
						hI[j] -= hK[i] * z[j] + z[i] * hK[j];
			}
		}
	}
	Ra_ev[N_sz-1] = Ra_evx[N_sz-1][N_sz-1];
#endif

	sseFree(z);
}

void tred2(wsReal **Ra_mat, int N_elem, wsReal *R_d, wsReal *R_e)
{
	int l, k, j, i;
	wsReal scale, hh, h, g, f;

	for (i=N_elem-1 ; i>=1 ; i--) {
		l = i-1;
		h = scale = 0.0;

		if (l > 0) {
			/*
			 * scale = sum(Ra_mat[i][1:(i-1)])
			 */
			scale = sseVasum(Ra_mat[i], i);
// 			for (k=0 ; k<i ; k++)
// 				scale += fabs(Ra_mat[i][k]);

			if (scale == 0.0)
				R_e[i] = Ra_mat[i][l];
			else {
				/*
				 * Ra_mat[i][1:(i-1)] /= scale
				 * h = sum(Ra_mat[i][1:(i-1)] * Ra_mat[i][1:(i-1)])
				 */
				sseVpC(Ra_mat[i], W1/scale, Ra_mat[i], i, &h);
// 				for (k=0 ; k<i ; k++) {
// 					Ra_mat[i][k] /= scale;
// 					h += SQR(Ra_mat[i][k]);
// 				}
				f = Ra_mat[i][l];
				g = f>0 ? -sqrt(h) : sqrt(h);
				R_e[i] = scale * g;
				h -= f * g;
				Ra_mat[i][l] = f - g;
				f = 0.0;

				/*
				 * Ra_mat[1:(i-1)][i] = Ra_mat[i][1:(i-1)]/h
				 * g =  sum(Ra_mat[ j in 1:(i-1) ][1:j] * Ra_mat[i][1:j])
				 * g += sum(Ra_mat[ j in 1:(i-1) ][j+1:(i-1)] * Ra_mat[i][j+1:(i-1)])
				 * R_e[1:(i-1)] = g / h
				 * f = sum(R_e[1:(i-1)] * Ra_mat[i][1:(i-1)])
				 */
				for (j=0 ; j<i ; j++) {
					Ra_mat[j][i] = Ra_mat[i][j]/h;
					g = 0.0;

					wsUint N_med = 0;
#ifdef USE_SSE
					N_med = getMed(j+1);
					sse_t sse_G = sseSet(0.0);
					for (wsUint l=0 ; l<N_med ; l+=sseJmp) {
						sse_t *sse_ik	= (sse_t *)(Ra_mat[i] + l);
						sse_t *sse_jk	= (sse_t *)(Ra_mat[j] + l);

						sse_G = sseAdd(sse_G, sseMul(*sse_ik, *sse_jk));
//						g += Ra_mat[j][k] * Ra_mat[i][k];
					}
					sseSum(sse_G, g);
#endif
					for (k=N_med ; k<=j ; k++)
						g += Ra_mat[j][k] * Ra_mat[i][k];
					for (k=j+1 ; k<i ; k++)
						g += Ra_mat[k][j] * Ra_mat[i][k];

					R_e[j] = g / h;
					f	+= R_e[j] * Ra_mat[i][j];
				}
				hh = f / (h + h);
				for (j=0 ; j<i ; j++) {
					f = Ra_mat[i][j];
					R_e[j] = g = R_e[j] - hh * f;

					wsUint N_med = 0;
#ifdef USE_SSE
					N_med = getMed(j+1);
					sse_t sse_F = sseSet(f);
					sse_t sse_G = sseSet(g);
					for (wsUint l=0 ; l<N_med ; l+=sseJmp) {
						sse_t *sse_ik	= (sse_t *)(Ra_mat[i] + l);
						sse_t *sse_jk	= (sse_t *)(Ra_mat[j] + l);
						sse_t *sse_k	= (sse_t *)(R_e + l);

						*sse_jk = sseSub(*sse_jk,
							sseAdd(sseMul(sse_F, *sse_k),
								sseMul(sse_G, *sse_ik)));
//						Ra_mat[j][k] -= f*R_e[k] + g*Ra_mat[i][k];
					}
#endif
					for (k=N_med ; k<=j ; k++)
						Ra_mat[j][k] -= f*R_e[k] + g*Ra_mat[i][k];
				}
			}
		} else
			R_e[i] = Ra_mat[i][l];
		R_d[i] = h;
	}
	R_d[0] = 0.0;
	R_e[0] = 0.0;
	for (i=0 ; i<N_elem ; i++) {
		if (R_d[i]) {
			for (j=0 ; j<i ; j++) {
				g = 0.0;
				for (k=0 ; k<i ; k++)	g += Ra_mat[i][k] * Ra_mat[k][j];
				for (k=0 ; k<i ; k++)	Ra_mat[k][j] -= g * Ra_mat[k][i];
			}
		}
		R_d[i]	= Ra_mat[i][i];
		Ra_mat[i][i]	= 1.0;
		for (j=0 ; j<i ; j++)
			Ra_mat[j][i] = Ra_mat[i][j] = 0.0;
	}
}

void tred2_dbl(double **Ra_mat, int N_elem, wsReal *R_d, wsReal *R_e)
{
	int l, k, j, i;
	double scale, hh, h, g, f;

	for (i=N_elem-1 ; i>=1 ; i--) {
		l = i-1;
		h = scale = 0.0;

		if (l > 0) {
			/*
			 * scale = sum(Ra_mat[i][1:(i-1)])
			 */
			for (k=0 ; k<i ; k++)
				scale += fabs(Ra_mat[i][k]);

			if (scale == 0.0)
				R_e[i] = Ra_mat[i][l];
			else {
				/*
				 * Ra_mat[i][1:(i-1)] /= scale
				 * h = sum(Ra_mat[i][1:(i-1)] * Ra_mat[i][1:(i-1)])
				 */
				for (k=0 ; k<i ; k++) {
					Ra_mat[i][k] /= scale;
					h += SQR(Ra_mat[i][k]);
				}
				f = Ra_mat[i][l];
				g = f>0 ? -sqrt(h) : sqrt(h);
				R_e[i] = scale * g;
				h -= f * g;
				Ra_mat[i][l] = f - g;
				f = 0.0;

				/*
				 * Ra_mat[1:(i-1)][i] = Ra_mat[i][1:(i-1)]/h
				 * g =  sum(Ra_mat[ j in 1:(i-1) ][1:j] * Ra_mat[i][1:j])
				 * g += sum(Ra_mat[ j in 1:(i-1) ][j+1:(i-1)] * Ra_mat[i][j+1:(i-1)])
				 * R_e[1:(i-1)] = g / h
				 * f = sum(R_e[1:(i-1)] * Ra_mat[i][1:(i-1)])
				 */
				for (j=0 ; j<i ; j++) {
					Ra_mat[j][i] = Ra_mat[i][j]/h;
					g = 0.0;

					wsUint N_med = 0;
#ifdef USE_SSE
					N_med = getMed(j+1);
					sse_t sse_G = sseSet(0.0);
					for (wsUint l=0 ; l<N_med ; l+=sseJmp) {
						sse_t *sse_ik	= (sse_t *)(Ra_mat[i] + l);
						sse_t *sse_jk	= (sse_t *)(Ra_mat[j] + l);

						sse_G = sseAdd(sse_G, sseMul(*sse_ik, *sse_jk));
//						g += Ra_mat[j][k] * Ra_mat[i][k];
					}
					sseSum(sse_G, g);
#endif
					for (k=N_med ; k<=j ; k++)
						g += Ra_mat[j][k] * Ra_mat[i][k];
					for (k=j+1 ; k<i ; k++)
						g += Ra_mat[k][j] * Ra_mat[i][k];

					R_e[j] = g / h;
					f	+= R_e[j] * Ra_mat[i][j];
				}
				hh = f / (h + h);
				for (j=0 ; j<i ; j++) {
					f = Ra_mat[i][j];
					R_e[j] = g = R_e[j] - hh * f;

					wsUint N_med = 0;
#ifdef USE_SSE
					N_med = getMed(j+1);
					sse_t sse_F = sseSet(f);
					sse_t sse_G = sseSet(g);
					for (wsUint l=0 ; l<N_med ; l+=sseJmp) {
						sse_t *sse_ik	= (sse_t *)(Ra_mat[i] + l);
						sse_t *sse_jk	= (sse_t *)(Ra_mat[j] + l);
						sse_t *sse_k	= (sse_t *)(R_e + l);

						*sse_jk = sseSub(*sse_jk,
							sseAdd(sseMul(sse_F, *sse_k),
								sseMul(sse_G, *sse_ik)));
//						Ra_mat[j][k] -= f*R_e[k] + g*Ra_mat[i][k];
					}
#endif
					for (k=N_med ; k<=j ; k++)
						Ra_mat[j][k] -= f*R_e[k] + g*Ra_mat[i][k];
				}
			}
		} else
			R_e[i] = Ra_mat[i][l];
		R_d[i] = h;
	}
	R_d[0] = 0.0;
	R_e[0] = 0.0;
	for (i=0 ; i<N_elem ; i++) {
		if (R_d[i]) {
			for (j=0 ; j<i ; j++) {
				g = 0.0;
				for (k=0 ; k<i ; k++)	g += Ra_mat[i][k] * Ra_mat[k][j];
				for (k=0 ; k<i ; k++)	Ra_mat[k][j] -= g * Ra_mat[k][i];
			}
		}
		R_d[i]	= Ra_mat[i][i];
		Ra_mat[i][i]	= 1.0;
		for (j=0 ; j<i ; j++)
			Ra_mat[j][i] = Ra_mat[i][j] = 0.0;
	}
}

void sseTred2(wsReal **Ra_mat, int N_elem, wsReal *R_d, wsReal *R_e)
{
	int l, k, j, i;
	wsReal scale, hh, h, g, f;

	for (i=N_elem-1 ; i>=1 ; i--) {
#ifdef USE_SSE
		int N_med = getMed(i);
#endif
		l = i-1;
		h = scale = 0.0;

		if (l > 0) {
			/*
			 * scale = sum(Ra_mat[i][1:(i-1)])
			 */
#ifdef USE_SSE
			sse_t sse_scale = sseSet(0.0);
			for (k=0 ; k<N_med ; k+=sseJmp)
				sse_scale = sseAdd(sse_scale,
					sseAbs(*(sse_t *)(Ra_mat[i] + k)));
			sseSum(sse_scale, scale);
			for (k=N_med ; k<i ; k++)
				scale += fabs(Ra_mat[i][k]);
#else
			for (k=0 ; k<i ; k++)
				scale += fabs(Ra_mat[i][k]);
#endif

			if (scale == W0)
				R_e[i] = Ra_mat[i][l];
			else {
				/*
				 * Ra_mat[i][1:(i-1)] /= scale
				 * h = sum(Ra_mat[i][1:(i-1)] * Ra_mat[i][1:(i-1)])
				 */
#ifdef USE_SSE
				sse_t sse_h = sseSet(0.0);
// 				sse_scale = sseSet(scale);
// 				for (k=0 ; k<N_med ; k+=sseJmp) {
// 					sse_t *sse_mat = (sse_t *)(Ra_mat[i] + k);
// 
// 					*sse_mat = sseDiv(*sse_mat, sse_scale);
// 					sse_h = sseAdd(sse_h, sseMul(*sse_mat, *sse_mat));
// 				}
// 				sseSum(sse_h, h);
// 				for (k=N_med ; k<i ; k++) {
// 					Ra_mat[i][k] /= scale;
// 					h += Ra_mat[i][k]*Ra_mat[i][k];
// 				}
				sse_scale = sseSet(W1/scale);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_mat = (sse_t *)(Ra_mat[i] + k);

					*sse_mat = sseMul(*sse_mat, sse_scale);
					sse_h = sseAdd(sse_h, sseMul(*sse_mat, *sse_mat));
				}
				sseSum(sse_h, h);
				for (k=N_med ; k<i ; k++) {
					Ra_mat[i][k] /= scale;
					h += Ra_mat[i][k]*Ra_mat[i][k];
				}
#else
				for (k=0 ; k<i ; k++) {
					Ra_mat[i][k] /= scale;
					h += Ra_mat[i][k]*Ra_mat[i][k];
				}
#endif
				f = Ra_mat[i][l];
				g = f>0 ? -sqrt(h) : sqrt(h);
				R_e[i] = scale * g;
				h -= f * g;
				Ra_mat[i][l] = f - g;
				f = 0.0;

				/*
				 * Ra_mat[1:(i-1)][i] = Ra_mat[i][1:(i-1)]/h
				 * g =  sum(Ra_mat[ j in 1:(i-1) ][1:j] * Ra_mat[i][1:j])
				 * g += sum(Ra_mat[ j in 1:(i-1) ][j+1:(i-1)] * Ra_mat[i][j+1:(i-1)])
				 * R_e[1:(i-1)] = g / h
				 * f = sum(R_e[1:(i-1)] * Ra_mat[i][1:(i-1)])
				 */
				for (j=0 ; j<i ; j++) {
					Ra_mat[j][i] = Ra_mat[i][j]/h;
					g = W0;

#ifdef USE_SSE
					int N_jMed = getMed(j+1);
					sse_t sse_g = sseSet(0.0);
					for (k=0 ; k<N_jMed ; k+=sseJmp) {
						sse_t *sse_jMat = (sse_t *)(Ra_mat[j] + k);
						sse_t *sse_iMat = (sse_t *)(Ra_mat[i] + k);

						// g += Ra_mat[j][k] * Ra_mat[i][k];
						sse_g = sseAdd(sse_g, sseMul(*sse_jMat, *sse_iMat));
					}
					sseSum(sse_g, g);
					for (k=N_jMed ; k<=j ; k++)
						g += Ra_mat[j][k] * Ra_mat[i][k];
#else
					for (k=0 ; k<=j ; k++)
						g += Ra_mat[j][k] * Ra_mat[i][k];
#endif
					for (k=j+1 ; k<=l ; k++)
						g += Ra_mat[k][j] * Ra_mat[i][k];

					R_e[j] = g / h;
					f	+= R_e[j] * Ra_mat[i][j];
				}
				hh = f / (h + h);
				for (j=0 ; j<i ; j++) {
					f = Ra_mat[i][j];
					R_e[j] = g = R_e[j] - hh * f;
#ifdef USE_SSE
					int N_jMed = getMed(j+1);
					sse_t sse_f = sseSet(f);
					sse_t sse_g = sseSet(g);
					for (k=0 ; k<N_jMed ; k+=sseJmp) {
						sse_t *sse_jMat = (sse_t *)(Ra_mat[j] + k);

						*sse_jMat = sseSub(*sse_jMat, sseAdd(
							sseMul(sse_f, *(sse_t *)(R_e + k)),
							sseMul(sse_g, *(sse_t *)(Ra_mat[i] + k))
							)
						);
					} 
					for (k=N_jMed ; k<=j ; k++)
						Ra_mat[j][k] -= f*R_e[k] + g*Ra_mat[i][k];
#else
					for (k=0 ; k<=j ; k++)
						Ra_mat[j][k] -= f*R_e[k] + g*Ra_mat[i][k];
#endif
				}
			}
		} else
			R_e[i] = Ra_mat[i][l];
		R_d[i] = h;
	}
	R_d[0] = 0.0;
	R_e[0] = 0.0;
	for (i=0 ; i<N_elem ; i++) {
		int N_med = getMed(i);
		if (R_d[i]) {
			for (j=0 ; j<i ; j++) {
				g = 0.0;
#ifdef USE_SSE
				sse_t sse_g = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_jMat = (sse_t *)(Ra_mat[j] + k);
					sse_t *sse_iMat = (sse_t *)(Ra_mat[i] + k);
					sse_g = sseAdd(sse_g, sseMul(*sse_jMat, *sse_iMat));
				}
				sseSum(sse_g, g);
				for (k=N_med ; k<i ; k++)
					g += Ra_mat[i][k] * Ra_mat[j][k];
#else
				for (k=0 ; k<i ; k++)
					g += Ra_mat[i][k] * Ra_mat[j][k];
#endif
#ifdef USE_SSE
				sse_g = sseSet(g);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_jMat = (sse_t *)(Ra_mat[j] + k);
					sse_t *sse_iMat = (sse_t *)(Ra_mat[i] + k);
					*sse_jMat = sseSub(*sse_jMat, sseMul(sse_g, *sse_iMat));
				}
				for (k=N_med ; k<i ; k++)
					Ra_mat[j][k] -= g * Ra_mat[i][k];
#else
				for (k=0 ; k<i ; k++)
					Ra_mat[j][k] -= g * Ra_mat[i][k];
#endif
			}
		}
		R_d[i]	= Ra_mat[i][i];
		Ra_mat[i][i]	= 1.0;
		for (j=0 ; j<i ; j++)
			Ra_mat[j][i] = Ra_mat[i][j] = 0.0;
	}
}

char tqli(wsReal *R_d, wsReal *R_e, int N_elem, wsReal **Ra_mat)
{
	int m, l, iter, i, k;
	wsReal s, r, p, g, f, dd, c, b;

	for (i=1 ; i<N_elem ; i++)
		R_e[i-1] = R_e[i];
	R_e[N_elem-1] = 0.0;
	for (l=0 ; l<N_elem ; l++) {
		iter = 0;
		do {
			for (m=l ; m<N_elem-1 ; m++) {
				dd = fabs(R_d[m]) + fabs(R_d[m+1]);
				if (fabs(R_e[m]) + dd == dd) break;
			} if (m!=l) {
				if (iter++ == 200) {
		//			LOG("No convergence in TLQI!");
					return 1;
				}

				g = (R_d[l+1] - R_d[l]) / (2.0f * R_e[l]);
				r = sqrt((g * g) + 1.0f);
				g = R_d[m] - R_d[l] + R_e[l] / (g + SIGN(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i=m-1 ; i>=l ; i--) {
					f = s*R_e[i];
					b = c*R_e[i];

					if (fabs(f) >= fabs(g)) {
						c = g / f;
						r = sqrt((c * c) + 1.0f);
						R_e[i+1] = f * r;
						c *= (s = 1.0f/r);
					} else {
						s = f / g;
						r = sqrt((s * s) + 1.0f);
						R_e[i+1] = g * r;
						s *= (c = 1.0f/r);
					}
					g = R_d[i+1] - p;
					r = (R_d[i] - g) * s + 2.0f * c * b;
					p = s * r;
					R_d[i+1] = g + p;
					g = c * r - b;
					for (k=0 ; k<N_elem ; k++) {
						f = Ra_mat[k][i+1];
						Ra_mat[k][i+1]	= s*Ra_mat[k][i] + c*f;
						Ra_mat[k][i]	= c*Ra_mat[k][i] - s*f;
					}
				}
				R_d[l] = R_d[l] - p;
				R_e[l] = g;
				R_e[m] = 0.0;
			}
		} while (m != l);
	}

	return 0;
}

char tqli_dbl(wsReal *R_d, wsReal *R_e, int N_elem, double **Ra_mat)
{
	int m, l, iter, i, k;
	wsReal s, r, p, g, f, dd, c, b;

	for (i=1 ; i<N_elem ; i++)
		R_e[i-1] = R_e[i];
	R_e[N_elem-1] = 0.0;
	for (l=0 ; l<N_elem ; l++) {
		iter = 0;
		do {
			for (m=l ; m<N_elem-1 ; m++) {
				dd = fabs(R_d[m]) + fabs(R_d[m+1]);
				if (fabs(R_e[m]) + dd == dd) break;
			} if (m!=l) {
				if (iter++ == 200) {
					//			LOG("No convergence in TLQI!");
					return 1;
				}

				g = (R_d[l+1] - R_d[l]) / (2.0f * R_e[l]);
				r = sqrt((g * g) + 1.0f);
				g = R_d[m] - R_d[l] + R_e[l] / (g + SIGN(r, g));
				s = c = 1.0;
				p = 0.0;
				for (i=m-1 ; i>=l ; i--) {
					f = s*R_e[i];
					b = c*R_e[i];

					if (fabs(f) >= fabs(g)) {
						c = g / f;
						r = sqrt((c * c) + 1.0f);
						R_e[i+1] = f * r;
						c *= (s = 1.0f/r);
					} else {
						s = f / g;
						r = sqrt((s * s) + 1.0f);
						R_e[i+1] = g * r;
						s *= (c = 1.0f/r);
					}
					g = R_d[i+1] - p;
					r = (R_d[i] - g) * s + 2.0f * c * b;
					p = s * r;
					R_d[i+1] = g + p;
					g = c * r - b;
					for (k=0 ; k<N_elem ; k++) {
						f = Ra_mat[k][i+1];
						Ra_mat[k][i+1]	= s*Ra_mat[k][i] + c*f;
						Ra_mat[k][i]	= c*Ra_mat[k][i] - s*f;
					}
				}
				R_d[l] = R_d[l] - p;
				R_e[l] = g;
				R_e[m] = 0.0;
			}
		} while (m != l);
	}

	return 0;
}

int jcbi(wsReal **Ra_mat, wsReal **Ra_eV, int N_elem, wsReal eps, int jt)
{
	int i,j,p=-1,q=-1,l;
	wsReal fm,cn,sn,omega,x,y,d;
	wsReal *Ra_colQ = NULL;
	wsReal *Ra_colP = NULL;

	wsAlloc(Ra_colQ, wsReal, N_elem);
	wsAlloc(Ra_colP, wsReal, N_elem);

	l=1;
	int nP = -1, nQ = -1;
	int nRecP = 0, nRecQ = 0;
	while (1) {
		fm=0.0;
		for (i=1 ; i<N_elem ; i++) {
			for (j=0 ; j<i ; j++) {
				d=fabs(Ra_mat[i][j]);
				if (d>fm) {
					fm	= d;
					nP	= i;
					nQ	= j;
				}
			}
		}
		if ((l%1000) == 0)
			printf("%d loops recP/Q %d/%d, %g\r", l, nRecP, nRecQ, fm);
		//printf("fm %g\n", fm);
		if (fm<eps) break;
		if (l>jt)
			break;
		wsReal *Ra_rowQ = Ra_mat[nQ];
		wsReal *Ra_rowP = Ra_mat[nP];
		if (nP != p) {
			for (i=0 ; i<N_elem ; i++)
				Ra_colP[i] = Ra_mat[i][nP];
		} else
			nRecP++;
		if (nQ != q) {
			for (i=0 ; i<N_elem ; i++)
				Ra_colQ[i] = Ra_mat[i][nQ];
		} else
			nRecQ++;
		p = nP;
		q = nQ;
		l++;
		x = -Ra_rowP[q];
		y = (Ra_rowQ[q]-Ra_rowP[p])/W2;
		omega=x/sqrt(x*x+y*y);

		if (y<W0)
			omega=-omega;
		sn = W1+sqrt(W1-omega*omega);
		sn = omega/sqrt(W2*sn);
		cn = sqrt(W1-sn*sn);
		fm = Ra_rowP[p];
		Ra_colP[p]=Ra_rowP[p]=fm*cn*cn+Ra_rowQ[q]*sn*sn+Ra_rowP[q]*omega;
		Ra_colQ[q]=Ra_rowQ[q]=fm*sn*sn+Ra_rowQ[q]*cn*cn-Ra_rowP[q]*omega;
		Ra_rowP[q]=W0;
		Ra_rowQ[p]=W0;
		Ra_colP[q]=W0;
		Ra_colQ[p]=W0;
		for (j=0 ; j<q; j++) {
			fm		= Ra_rowP[j];
			Ra_rowP[j]	= fm*cn+Ra_rowQ[j]*sn;
			Ra_rowQ[j]	= Ra_rowQ[j]*cn-fm*sn;
			fm		= Ra_colP[j];
			Ra_colP[j]	= fm*cn+Ra_colQ[j]*sn;
			Ra_colQ[j]	= -fm*sn+Ra_colQ[j]*cn;
		}
		for (j=q+1 ; j<p ; j++) {
			fm		= Ra_rowP[j];
			Ra_rowP[j]	= fm*cn+Ra_rowQ[j]*sn;
			Ra_rowQ[j]	= Ra_rowQ[j]*cn-fm*sn;
			fm		= Ra_colP[j];
			Ra_colP[j]	= fm*cn+Ra_colQ[j]*sn;
			Ra_colQ[j]	= -fm*sn+Ra_colQ[j]*cn;
		}
		for (j=p+1 ; j<N_elem ; j++) {
			fm		= Ra_rowP[j];
			Ra_rowP[j]	= fm*cn+Ra_rowQ[j]*sn;
			Ra_rowQ[j]	= Ra_rowQ[j]*cn-fm*sn;
			fm		= Ra_colP[j];
			Ra_colP[j]	= fm*cn+Ra_colQ[j]*sn;
			Ra_colQ[j]	= Ra_colQ[j]*cn-fm*sn;
		}
		for (i=0 ; i<q ; i++) {
			Ra_mat[i][p] = Ra_colP[i];
			Ra_mat[i][q] = Ra_colQ[i];
		}
		for (i=q+1 ; i<p ; i++) {
			Ra_mat[i][p] = Ra_colP[i];
			Ra_mat[i][q] = Ra_colQ[i];
		}
		for (i=p+1 ; i<N_elem ; i++) {
			Ra_mat[i][p] = Ra_colP[i];
			Ra_mat[i][q] = Ra_colQ[i];
		}
	}
	sseMalloc(*Ra_eV, wsReal, N_elem);
	for (i=0 ; i<N_elem ; i++)
		(*Ra_eV)[i] = Ra_mat[i][i];
	//printf("%d loops recP/Q %d/%d\n", l, nRecP, nRecQ);

	DEALLOC(Ra_colQ);
	DEALLOC(Ra_colP);

	return(1);
}

wsReal* EVpower(wsReal **Ra_mat, wsUint N_sz, wsUint N_pc, double R_thr)
{
	wsUint i, j, k;
	wsReal *Ra_ev	= NULL;
	sseMalloc(Ra_ev, wsReal, N_pc);
	wsReal *ttt		= NULL;

	/* Matrix x (initial guess) */
	wsReal *xk = NULL;	///< Matrix xk
	wsReal m0;	///< x * x variable
	wsReal m1;	///< x * xk variable
	wsReal m2;	///< xk * xk variable
	wsReal q = W0;	///< Calculated Raleigh quotient
	wsReal e;	///< Calculated error bound
	wsReal max;	///< Maximum magnitude of xk

	wsAlloc(xk, wsReal, N_sz);
	sseMalloc(ttt, wsReal, N_sz);

	/* Calculate xk vector, Raleigh Quotient, and Error. */
	for (wsUint N_iter=0 ; N_iter<N_pc ; N_iter++) {
		// x0<-matrix(rep(1,ncol(dat)),ncol=1)
		for (i=0 ; i<N_sz ; i++)
			ttt[i] = 1.0;
		double pe = -999;

		for (k=0 ; k<100 ; k++) {
			/* Load matrix 'x' w/ zeros. */
			memset(xk, 0x00, sizeof(wsReal)*N_sz);

			/* Set coeffs to zero */
			m0 = m1 = m2 = 0;
			max = 0;

			/* Get Ax */
			for(i=0; i<N_sz ; i++) {
				wsReal *Ra_mt = sseMpV(&Ra_mat[i], 1, N_sz, ttt);
				xk[i] += Ra_mt[0];
				sseFree(Ra_mt);
//				for(j=0; j<N_sz; j++) /* Calculate A * x. */
//					xk[i] += Ra_mat[i][j]*ttt[j];
				m0 += ttt[i]*ttt[i];   /* x * x vector   */
				m1 += ttt[i]*xk[i];  /* x * xk vector  */
				m2 += xk[i]*xk[i]; /* xk * xk vector */

				/* Find max magnitude of xk vector */
				if (xk[i]>max) max=xk[i];
				else if (0-xk[i]>max) max=0-xk[i];
			}

			/* Display x(k) vector */
			for(i=0; i<N_sz ; i++)
				ttt[i] = xk[i]/max;

			q = m1/m0; /* Calculate Raleigh quotient */
			e = sqrt((m2/m0)-(q*q)); /* Calculate error bound */

			if (pe != -999) {
				//				LOG("Loop %d e %g...\r", k, e);
				if (e < R_thr) {
					// 					LOG("\nx%d=[",k);
					// 					for(i=0; i<N ; i++)
					// 						LOG("%lf ",ttt[i]);
					// 					if (N_sz!=N)
					// 						LOG("...");
					// 					LOG("] q%d=%g e=%g\n",q,N_iter,e);
					break;
				}
			} else
				pe = e;
		}
// 		LOG("Iteration %d : (%s:%d) e.value=%g, error=%g\n",
// 			N_iter+1, k==100?"STOP":"OK", k, q, e);

		/* Substract */
		wsReal xxx = 0;
		for (j=0 ; j<N_sz ; j++)
			xxx += ttt[j]*ttt[j];
		xxx = sqrt(xxx);
		Ra_ev[N_iter] = q;

		for (j=0 ; j<N_sz ; j++)
			ttt[j] /= xxx;
		//		LOG("Extract matrix\n");
		for (i=0 ; i<N_sz ; i++) {
			wsReal xx = 0;
			for (j=0 ; j<N_sz ; j++)
				xx += Ra_mat[i][j]*ttt[j];

			/* Get vec */
			for (j=0 ; j<N_sz ; j++) {
				Ra_mat[i][j] -= xx*ttt[j];
				//				printf("%g ", A[i][j]);
			}
			//			printf("\n");
		}
	}
	sseFree(ttt);
	DEALLOC(xk);

	qsort(Ra_ev, N_pc, sizeof(wsReal), sortabs_real_unorder);
	return Ra_ev;
}

inline wsReal r_sign(wsReal *a, wsReal *b)
{
	wsReal x;
	x = (*a >= 0 ? *a : - *a);
	return( *b >= 0 ? x : -x);
}

void tred1(wsUint N_sz, wsReal **Ra_mat, wsReal *Ra_diag, wsReal *e)
{
	/* System generated locals */
	wsReal r__1;

	/* Local variables */
	static wsReal f, g, h__;
	wsUint i, ii, j, k, l, jp1;
	static wsReal scale;

	/*     this subroutine is a translation of the algol procedure tred1, */
	/*     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson. */
	/*     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971). */

	/*     this subroutine reduces a real symmetric matrix */
	/*     to a symmetric tridiagonal matrix using */
	/*     orthogonal similarity transformations. */

	/*     on input */

	/*        N_sz is the order of the matrix. */

	/*        Ra_mat contains the real symmetric input matrix.  only the */
	/*          lower triangle of the matrix need be supplied. */

	/*     on output */

	/*        Ra_mat contains information about the orthogonal trans- */
	/*          formations used in the reduction in its strict lower */
	/*          triangle.  the full upper triangle of a is unaltered. */

	/*        Ra_diag contains the diagonal elements of the tridiagonal matrix. */

	/*        e contains the subdiagonal elements of the tridiagonal */
	/*          matrix in its last n-1 positions.  e(1) is set to zero. */

	/*        e2 contains the squares of the corresponding elements of e. */
	/*          e2 may coincide with e if the squares are not needed. */

	/*     questions and comments should be directed to burton s. garbow, */
	/*     mathematics and computer science div, argonne national laboratory */

	/*     this version dated august 1983. */

	/*     ------------------------------------------------------------------ */

	/* Parameter adjustments */
	--e;
	--Ra_diag;
	--Ra_mat;

	/* Function Body */
	for (i = 1; i <= N_sz; ++i) {
		Ra_diag[i] = Ra_mat[i][N_sz-1];
		Ra_mat[i][N_sz-1] = Ra_mat[i][i-1];
		/* L100: */
	}
	/*     .......... for i=n step -1 until 1 do -- .......... */
	for (ii = 1; ii <= N_sz ; ++ii) {
		i = N_sz + 1 - ii;
		l = i - 1;
		h__ = W0;
		scale = W0;
		if (l < 1) {
			goto L130;
		}
		/*     .......... scale row (algol tol then not needed) .......... */
#if 0
		static const __m128 sse_absMask = 
			_mm_castsi128_ps(_mm_set1_epi32(0x80000000));
		UINT_t	N_med		= getMed(l+1);
		__m128	sse_scale	= _mm_set1_ps(0.0f);
		for (k=1 ; k<=N_med ; k+=sseJmp)
			sse_scale = _mm_add_ps(sse_scale,
				_mm_andnot_ps(sse_absMask, *((__m128 *)Ra_diag+k)));
		sseSum(sse_scale, scale);
		for (k=N_med+1 ; k<=l ; ++k)
			scale += fabs(Ra_diag[k]);
#else
		for (k = 1; k <= l; ++k)
			scale += fabs(Ra_diag[k]);
#endif

		if (scale != 0.f) {
			goto L140;
		}

		for (j = 1; j <= l; ++j) {
			Ra_diag[j] = Ra_mat[j][l-1];
			Ra_mat[j][l-1] = Ra_mat[j][i-1];
			Ra_mat[j][i-1] = W0;
			/* L125: */
		}

L130:
		e[i] = W0;
		goto L300;

L140:
		for (k = 1; k <= l; ++k) {
			Ra_diag[k] /= scale;
			h__ += Ra_diag[k] * Ra_diag[k];
			/* L150: */
		}

		f = Ra_diag[l];
		r__1 = sqrt(h__);
		g = -r_sign(&r__1, &f);
		e[i] = scale * g;
		h__ -= f * g;
		Ra_diag[l] = f - g;
		if (l == 1) {
			goto L285;
		}
		/*     .......... form a*u .......... */
		for (j = 1; j <= l; ++j) {
			/* L170: */
			e[j] = W0;
		}

		for (j = 1; j <= l; ++j) {
			f = Ra_diag[j];
			g = e[j] + Ra_mat[j][j-1] * f;
			jp1 = j + 1;
			if (l < jp1) {
				goto L220;
			}

			for (k = jp1; k <= l; ++k) {
				g += Ra_mat[j][k-1] * Ra_diag[k];
				e[k] += Ra_mat[j][k-1] * f;
				/* L200: */
			}

L220:
			e[j] = g;
			/* L240: */
		}
		/*     .......... form p .......... */
		f = W0;

		for (j = 1; j <= l; ++j) {
			e[j] /= h__;
			f += e[j] * Ra_diag[j];
			/* L245: */
		}

		h__ = f / (h__ + h__);
		/*     .......... form q .......... */
		for (j = 1; j <= l; ++j) {
			/* L250: */
			e[j] -= h__ * Ra_diag[j];
		}
		/*     .......... form reduced a .......... */
		for (j = 1; j <= l; ++j) {
			f = Ra_diag[j];
			g = e[j];

			for (k = j; k <= l; ++k)
				/* L260: */
				Ra_mat[j][k-1] -= f * e[k] - g * Ra_diag[k];
		}

L285:
		for (j = 1; j <= l; ++j) {
			f = Ra_diag[j];
			Ra_diag[j] = Ra_mat[j][l-1];
			Ra_mat[j][l-1] = Ra_mat[j][i-1];
			Ra_mat[j][i-1] = f * scale;
			/* L290: */
		}

L300:
		;
	}
} /* tred1_ */
	
//void tqli1(wsUint N_sz, wsReal *d__, wsReal *e, wsUint *ierr)
int tqli1(int n, wsReal *d, wsReal *e)
{
	d--;
	e--;
	int m,l,iter,i;
	wsReal s,r,p,g,f,dd,c,b;

	for (i=2;i<=n;i++) e[i-1]=e[i];
	e[n]=W0;
	for (l=1;l<=n;l++) {
		iter=0;
		do {
			for (m=l ; m<=n-1 ; m++) {
				dd=fabs(d[m])+fabs(d[m+1]);
				if ((double)(fabs(e[m])+dd) == dd) break;
			}
			if (m != l) {
				if (iter++ == 100)
					return 1;
				g=(d[l+1]-d[l])/(W2*e[l]);
				r=(wsReal)pythag(g, W1);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=W1;
				p=W0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=(wsReal)pythag(f,g));
					if (r == W0) {
						d[i+1] -= p;
						e[m]=W0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+W2*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
				}
				if (r == W0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=W0;
			}
		} while (m != l);
	}

	qsort(d+1, n, sizeof(wsReal), sortabs_real_unorder);
	return 0;
}

int tqli2(wsReal *d, wsReal *e, int n, 
		  wsReal **z, char B_getEvec)
{
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;

	for (l=0;l<n;l++) {
		iter=0;
		do { 
			for (m=l;m<n-1;m++) { 
				dd=fabs(d[m])+fabs(d[m+1]);
				if (fabs(e[m])+dd == dd) break;
			}
			if (m!=l) { 
				if (iter++ == 30) { 
					exportMatrix("d.tqli2", &d, 1, n);
					exportMatrix("e.tqli2", &e, 1, n);
					halt("Too many iterations in tqli()\n");
					return 0;
				}
				g=(d[l+1]-d[l])/(W2*e[l]);
				r=sqrt((g*g)+W1);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) { 
					f=s*e[i];
					b=c*e[i];
					if (fabs(f) >= fabs(g)) { 
						c=g/f;r=sqrt((c*c)+W1);
						e[i+1]=(wsReal)(f*r);
						c *= (s=1.0/r);
					} else { 
						s=f/g;r=sqrt((s*s)+W1);
						e[i+1]=(wsReal)(g*r);
						s *= (c=1.0/r);
					}
					g=d[i+1]-p;
					r=(d[i]-g)*s+W2*c*b;
					p=s*r;
					d[i+1]=(wsReal)(g+p);
					g=c*r-b;
					/*EVECTS*/
					if (B_getEvec == 1) {
						for (k=0;k<n;k++) { 
							f=z[k][i+1];
							z[k][i+1]=(wsReal)(s*z[k][i]+c*f);
							z[k][i]=(wsReal)(c*z[k][i]-s*f);
						}
					}//Evects
				}
				d[l] -=	(wsReal)p;
				e[l] =	(wsReal)g;
				e[m] =	W0;
			}
		} while (m!=l);
	}

	return 1;
}

wsReal* eigen(wsReal **Ra_mat, wsUint N_sz, char *Bp_stat)
{
	wsUint	N_ret;
	wsReal	* Ra_e, *Ra_diag = NULL;
	sseMalloc(Ra_e, wsReal, N_sz+1);
	sseMalloc(Ra_diag, wsReal, N_sz+1);

	cTimer c;
	//	c.start();
#if 1
	sseTred2(Ra_mat, N_sz, Ra_diag, Ra_e);
#else
	tred2(Ra_mat, N_sz, Ra_diag, Ra_e);
#endif
	//	LOG("Tred2 %s\n", c.getReadable());
	//	exportMatrix("eig.e", &Ra_e, 1, N_sz+1);
	//	exportMatrix("eig.diag", &Ra_diag, 1, N_sz);
	//	c.start();
	N_ret = tqli1(N_sz, Ra_diag, Ra_e);
	//	LOG("tqli1 %s\n", c.getReadable());

	//	exportMatrix("eig.e2", &Ra_e, 1, N_sz+1);
	//	exportMatrix("eig.diag2", &Ra_diag, 1, N_sz);

	*Bp_stat = (char)N_ret;
	sseFree(Ra_e);

	return Ra_diag;
}

wsReal** sqrtMatrix(wsReal **Ra_mat, wsUint N_sz, char B_sym/*=0*/)
{
	wsUint	i, j, k;
//	char	B_isFailed = 0;
	wsReal	**Ra_eVec	= NULL; /* P, not P^t */
	wsReal	*Ra_eVal	= NULL;
	
#ifdef USE_ED2
	Ra_eVal = EIGENDECOMPOSITION(Ra_mat, N_sz, &Ra_eVec, 1);
#else
	if (B_sym) {
		wsReal **Ra_m = sseSym2Mat(Ra_mat, N_sz);
		Ra_eVal = eigenDecomp(Ra_m, N_sz, &B_isFailed, &Ra_eVec);
		sseUnmat(Ra_m, N_sz);
	} else
		Ra_eVal = eigenDecomp(Ra_mat, N_sz, &B_isFailed, &Ra_eVec);
#endif

	wsReal	*Ra_b		= NULL;
	wsReal	**Ra_ret	= sseEmptyMatrix(N_sz, N_sz);
	sseMalloc(Ra_b, wsReal, N_sz);

	/* Take sqrt to eVal */
	// Look for singular values
	const wsReal	R_eps		= numeric_limits<wsReal>::epsilon();
	wsReal wmax = W0;
	for (wsUint i=0; i<N_sz; i++)
		wmax = Ra_eVal[i] > wmax ? Ra_eVal[i] : wmax;
	wsReal wmin = wmax * R_eps;
	for (wsUint i=0; i<N_sz; i++)
		Ra_eVal[i] = Ra_eVal[i] < wmin ? 0 : sqrt(Ra_eVal[i]);
//	sseVsqrt(Ra_eVal, N_sz);

	/* Do evec%*%sqrt(eVal)%*%t(evec)
	 *   [ e11v11e11+...+e1nvnne1n e11v11e21+...+e1nvnne2n ... e11v11en1+...+e1nvnnenn ]
	 * = [ e21v11e11+...+e2nvnne1n e21v11e21+...+e2nvnne2n ... e21v11en1+...+e2nvnnenn ]
	 *   [          . . .                   . . .                       . . .
	 *   [ eikvkkejk
	 *   [ en1v11e11+...+ennvnne1n en1v11e21+...+ennvnne2n ... en1v11en1+...+ennvnnenn ] */
#ifdef USE_SSE
	wsUint N_med = getMed(N_sz);
#endif
	for (i=0 ; i<N_sz ; i++) {
#ifdef USE_SSE
		for (k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *sse_b = (sse_t *)(Ra_b + k);
			sse_t *sse_v = (sse_t *)(Ra_eVal + k);
			sse_t *sse_e = (sse_t *)(Ra_eVec[i] + k);

			*sse_b = sseMul(*sse_e, *sse_v);
		}
		for (k=N_med ; k<N_sz ; k++)
			Ra_b[k] = Ra_eVec[i][k]*Ra_eVal[k];
#else
		for (k=0 ; k<N_sz ; k++)
			Ra_b[k] = Ra_eVec[i][k]*Ra_eVal[k];
#endif
		for (j=0 ; j<N_sz ; j++) {
#ifdef USE_SSE
			sse_t sse_sum = sseSet(0.0);
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_b = (sse_t *)(Ra_b + k);
				sse_t *sse_v = (sse_t *)(Ra_eVec[j] + k);
				sse_sum = sseAdd(sse_sum, sseMul(*sse_b, *sse_v));
			}
			sseSum(sse_sum, Ra_ret[i][j]);

			for (k=N_med ; k<N_sz ; k++)
				Ra_ret[i][j] += Ra_b[k]*Ra_eVec[j][k];
#else
			for (k=0 ; k<N_sz ; k++)
				Ra_ret[i][j] += Ra_b[k]*Ra_eVec[j][k];
#endif
		}
	}
	sseFree(Ra_b);
	sseFree(Ra_eVal);
	sseUnmat(Ra_eVec, N_sz);

	return Ra_ret;
}

wsReal** sqrtMatrixEV(wsReal **Ra_eVec, wsReal *Ra_eVal, wsUint N_sz,
	char B_sym/*=0*/)
{
	wsUint	i, j, k;
	
	wsReal	*Ra_b		= NULL;
	wsReal	**Ra_ret	= sseEmptyMatrix(N_sz, N_sz);
	sseMalloc(Ra_b, wsReal, N_sz);

	/* Take sqrt to eVal */
	const wsReal	R_eps		= numeric_limits<wsReal>::epsilon();
	wsReal wmax = W0;
	for (i=0; i<N_sz; i++)
		wmax = Ra_eVal[i] > wmax ? Ra_eVal[i] : wmax;
	wsReal wmin = wmax * R_eps;
	for (i=0; i<N_sz; i++)
		Ra_eVal[i] = Ra_eVal[i] < wmin ? 0 : sqrt(Ra_eVal[i]);

	/* Do evec%*%sqrt(eVal)%*%t(evec)
	 *   [ e11v11e11+...+e1nvnne1n e11v11e21+...+e1nvnne2n ... e11v11en1+...+e1nvnnenn ]
	 * = [ e21v11e11+...+e2nvnne1n e21v11e21+...+e2nvnne2n ... e21v11en1+...+e2nvnnenn ]
	 *   [          . . .                   . . .                       . . .
	 *   [ eikvkkejk
	 *   [ en1v11e11+...+ennvnne1n en1v11e21+...+ennvnne2n ... en1v11en1+...+ennvnnenn ] */
#ifdef USE_SSE
	wsUint N_med = getMed(N_sz);
#endif
	for (i=0 ; i<N_sz ; i++) {
#ifdef USE_SSE
		for (k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *sse_b = (sse_t *)(Ra_b + k);
			sse_t *sse_v = (sse_t *)(Ra_eVal + k);
			sse_t *sse_e = (sse_t *)(Ra_eVec[i] + k);

			*sse_b = sseMul(*sse_e, *sse_v);
		}
		for (k=N_med ; k<N_sz ; k++)
			Ra_b[k] = Ra_eVec[i][k]*Ra_eVal[k];
#else
		for (k=0 ; k<N_sz ; k++)
			Ra_b[k] = Ra_eVec[i][k]*Ra_eVal[k];
#endif
		for (j=0 ; j<N_sz ; j++) {
#ifdef USE_SSE
			sse_t sse_sum = sseSet(0.0);
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_b = (sse_t *)(Ra_b + k);
				sse_t *sse_v = (sse_t *)(Ra_eVec[j] + k);
				sse_sum = sseAdd(sse_sum, sseMul(*sse_b, *sse_v));
			}
			sseSum(sse_sum, Ra_ret[i][j]);

			for (k=N_med ; k<N_sz ; k++)
				Ra_ret[i][j] += Ra_b[k]*Ra_eVec[j][k];
#else
			for (k=0 ; k<N_sz ; k++)
				Ra_ret[i][j] += Ra_b[k]*Ra_eVec[j][k];
#endif
		}
	}
	sseFree(Ra_b);
// 	sseFree(Ra_eVal);
// 	sseUnmat(Ra_eVec, N_sz);

	return Ra_ret;
}

wsReal** eigenInvMatrix(wsReal **Ra_mat, wsUint N_sz)
{
	wsUint	i, j, k;
	char	B_isFailed = 0;
	wsReal	**Ra_eVec	= NULL;
	wsReal	*Ra_eVal	= eigenDecomp(Ra_mat, N_sz, &B_isFailed, &Ra_eVec);
	wsReal	*Ra_b		= NULL;
	wsReal	**Ra_ret	= sseEmptyMatrix(N_sz, N_sz);
	sseMalloc(Ra_b, wsReal, N_sz);

	/* Take sqrt to eVal */
	sseVsqrt(Ra_eVal, N_sz);

	/* Do evec%*%sqrt(eVal)%*%t(evec)
	 *   [ e11v11e11+...+e1nvnne1n e11v11e21+...+e1nvnne2n ... e11v11en1+...+e1nvnnenn ]
	 * = [ e21v11e11+...+e2nvnne1n e21v11e21+...+e2nvnne2n ... e21v11en1+...+e2nvnnenn ]
	 *   [          . . .                   . . .                       . . .
	 *   [ eikvkkejk
	 *   [ en1v11e11+...+ennvnne1n en1v11e21+...+ennvnne2n ... en1v11en1+...+ennvnnenn ] */
#ifdef USE_SSE
	wsUint N_med = getMed(N_sz);
#endif
	for (i=0 ; i<N_sz ; i++) {
#ifdef USE_SSE
		for (k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *sse_b = (sse_t *)(Ra_b + k);
			sse_t *sse_v = (sse_t *)(Ra_eVal + k);
			sse_t *sse_e = (sse_t *)(Ra_eVec[i] + k);

			*sse_b = sseMul(*sse_e, *sse_v);
		}
		for (k=N_med ; k<N_sz ; k++)
			Ra_b[k] = Ra_eVec[i][k]*Ra_eVal[k];
#else
		for (k=0 ; k<N_sz ; k++)
			Ra_b[k] = Ra_eVec[i][k]*Ra_eVal[k];
#endif
		for (j=0 ; j<N_sz ; j++) {
#ifdef USE_SSE
			sse_t sse_sum = sseSet(0.0);
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_b = (sse_t *)(Ra_b + k);
				sse_t *sse_v = (sse_t *)(Ra_eVec[j] + k);
				sse_sum = sseAdd(sse_sum, sseMul(*sse_b, *sse_v));
			}
			sseSum(sse_sum, Ra_ret[i][j]);

			for (k=N_med ; k<N_sz ; k++)
				Ra_ret[i][j] += Ra_b[k]*Ra_eVec[j][k];
#else
			for (k=0 ; k<N_sz ; k++)
				Ra_ret[i][j] += Ra_b[k]*Ra_eVec[j][k];
#endif
		}
	}
	sseFree(Ra_b);
	sseFree(Ra_eVal);
	sseUnmat(Ra_eVec, N_sz);

	return Ra_ret;
}

wsMat eigenInvSymMatrix(wsMat Ra_eVec, wsReal *Ra_eVal, wsUint N_sz)
{
	wsUint	i, j, k;
	wsReal	*Ra_b		= NULL;
	wsReal	**Ra_ret	= sseEmptyMatrix(N_sz, N_sz);
	sseMalloc(Ra_b, wsReal, N_sz);
	wsReal	*Ra_wgt		= Ra_eVal;
	const wsReal	R_eps		= numeric_limits<wsReal>::epsilon();

	// Look for singular values
	wsReal wmax = W0;
	for (i=0; i<N_sz; i++)
		wmax = Ra_eVal[i] > wmax ? Ra_eVal[i] : wmax;
	wsReal wmin = wmax * R_eps;
	for (i=0; i<N_sz; i++)
		Ra_wgt[i] = Ra_eVal[i] < wmin ? 0 : 1/Ra_eVal[i];


	/* Do evec%*%sqrt(eVal)%*%t(evec)
	 *   [ e11v11e11+...+e1nvnne1n e11v11e21+...+e1nvnne2n ... e11v11en1+...+e1nvnnenn ]
	 * = [ e21v11e11+...+e2nvnne1n e21v11e21+...+e2nvnne2n ... e21v11en1+...+e2nvnnenn ]
	 *   [          . . .                   . . .                       . . .
	 *   [ eikvkkejk
	 *   [ en1v11e11+...+ennvnne1n en1v11e21+...+ennvnne2n ... en1v11en1+...+ennvnnenn ] */
#ifdef USE_SSE
	wsUint N_med = getMed(N_sz);
#endif
	for (i=0 ; i<N_sz ; i++) {
#ifdef USE_SSE
		for (k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *sse_b = (sse_t *)(Ra_b + k);
			sse_t *sse_v = (sse_t *)(Ra_eVal + k);
			sse_t *sse_e = (sse_t *)(Ra_eVec[i] + k);

			*sse_b = sseMul(*sse_e, *sse_v);
		}
		for (k=N_med ; k<N_sz ; k++)
			Ra_b[k] = Ra_eVec[i][k]*Ra_eVal[k];
#else
		for (k=0 ; k<N_sz ; k++)
			Ra_b[k] = Ra_eVec[i][k]*Ra_eVal[k];
#endif
		for (j=0 ; j<=i ; j++) {
#ifdef USE_SSE
			sse_t sse_sum = sseSet(0.0);
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_b = (sse_t *)(Ra_b + k);
				sse_t *sse_v = (sse_t *)(Ra_eVec[j] + k);
				sse_sum = sseAdd(sse_sum, sseMul(*sse_b, *sse_v));
			}
			sseSum(sse_sum, Ra_ret[i][j]);

			for (k=N_med ; k<N_sz ; k++)
				Ra_ret[i][j] += Ra_b[k]*Ra_eVec[j][k];
#else
			for (k=0 ; k<N_sz ; k++)
				Ra_ret[i][j] += Ra_b[k]*Ra_eVec[j][k];
#endif
		}
	}
	sseFree(Ra_b);
	//sseFree(Ra_eVal);
	//sseUnmat(Ra_eVec, N_sz);

	return Ra_ret;
}

wsReal* eigenDecomp2(wsReal **Ra_mat, wsUint N_sz, char *Bp_stat,
	wsReal ***Rp_eVec, char B_noTranspose/*=0*/, char B_negCorrect/*=1*/)
{
	char	N_ret;
	wsReal	*Ra_e, *Ra_diag = NULL;

	/* Size check */
	if (N_sz == 0)
		return NULL;

	sseMalloc(Ra_e, wsReal, N_sz+1);
	sseMalloc(Ra_diag, wsReal, N_sz+1);
	wsReal	**Ra_eVec = NULL;
	Ra_eVec = sseMatrixP(N_sz, N_sz, Ra_mat);

	cTimer xxx;
	xxx.start();
	trdg(Ra_eVec, N_sz, Ra_diag, Ra_e);
	/* Print only it takes more than 1sec */
	char B_prt = false;
	if (xxx.get() > REAL_CONST(1000.0)) {
		pverbose("TRDG %s\n", xxx.getReadable());
		B_prt = true;
	}
	//exportMatrix("EVX", Ra_eVec, N_sz, N_sz);
	//wsReal	**Ra_eVec2 = Ra_eVec;
	wsReal **Ra_eVec2 = sseEmptyMatrix(N_sz, N_sz);

	// build up first part of the matrix by applying Householder transforms
	xxx.start();
	if (Rp_eVec) {
#if 1
		for (int k=N_sz-1 ; k>=1 ; k--) {
#ifdef USE_SSE
			/* Minimum starting point of SSE */
			wsUint N1 = k+1 + sseJmp-1;
			N1 = getMed(N1);
			/* Maximum reach point of SSE */
			//wsUint N2 = getMed(N_sz);
			/* Use SSE method ONLY IF k<=N1<N2<=N_sz */
#endif
			wsReal	*hK = Ra_eVec[k - 1];
			wsReal	*cK = Ra_eVec2[k];
			//memset(cK, 0x00, sizeof(wsReal)*N_sz);
	//		Ra_eVec2[k][k] = 1;
			cK[k] = 1;
			if (hK[k] != 0.0) {
				wsReal inv = 1.0 / (Ra_e[k - 1] * hK[k]);
				wsReal beta = 1.0 / Ra_e[k - 1];
				cK[k] = 1 + beta * hK[k];

// #ifdef USE_SSE
// 				if (N1 < N2) {
// 					/* k -> N1 */
// 					for (wsUint i=(wsUint)k+1 ; i<N1 ; i++)
// 						cK[i] = beta * hK[i];
// 					sse_t sse_beta = sseSet(beta);
// 					for (wsUint i=N1 ; i<N2 ; i+=sseJmp) {
// 						*((sse_t *)(cK + i)) = sseMul(sse_beta,
// 							*((sse_t *)(hK + i)));
// 	//					cK[i] = beta * hK[i];
// 					}
// 					for (wsUint i=N2 ; i<N_sz ; i++)
// 						cK[i] = beta * hK[i];
// 				} else 
// #endif
					for (wsUint i=(wsUint)k+1 ; i<N_sz ; ++i)
						cK[i] = beta * hK[i];

	//			cK[k] = W1;
				for (wsUint j=(wsUint)(k+1) ; j<N_sz ; ++j) {
					wsReal *cJ = Ra_eVec2[j];
					beta = 0;
// #ifdef USE_SSE
// 					if (N1 < N2) {
// 						/* k -> N1 */
// 						for (wsUint i=(wsUint)k+1 ; i<N1 ; i++)
// 							beta += cJ[i] * hK[i];
// 						sse_t sse_beta = sseSet(0.0);
// 						for (wsUint i=N1 ; i<N2 ; i+=sseJmp) {
// 							sse_beta = sseAdd(sse_beta, sseMul(
// 								*((sse_t *)(cJ + i)), *((sse_t *)(hK + i))));
// 	//						beta += cJ[i] * hK[i];
// 						}
// 						sseAppendSum(sse_beta, beta);
// 						for (wsUint i=N2 ; i<N_sz ; i++)
// 							beta += cJ[i] * hK[i];
// 					} else
// #endif
						for (wsUint i=(wsUint)k+1 ; i<N_sz ; ++i)
							beta += cJ[i] * hK[i];

					beta *= inv;
					cJ[k] = beta * hK[k];
// #ifdef USE_SSE
// 					if (N1 < N2) {
// 						/* k -> N1 */
// 						for (wsUint i=(wsUint)k+1 ; i<N1 ; i++)
// 							cJ[i] += beta * hK[i];
// 						sse_t sse_beta = sseSet(beta);
// 						for (wsUint i=N1 ; i<N2 ; i+=sseJmp) {
// 							sse_t *sse_cJ = (sse_t *)(cJ + i);
// 							*sse_cJ = sseAdd(*sse_cJ, sseMul(
// 								sse_beta, *((sse_t *)(hK + i))));
// 	//						cJ[i] += beta * hK[i];
// 						}
// 						//sseAppendSum(sse_beta, beta);
// 						for (wsUint i=N2 ; i<N_sz ; i++)
// 							cJ[i] += beta * hK[i];
// 					} else
// #endif
						for (wsUint i=(wsUint)k+1 ; i<N_sz ; ++i)
							cJ[i] += beta * hK[i];
				}
			}
		}
		memset(Ra_eVec2[0], 0x00, sizeof(wsReal)*N_sz);
		//exportMatrix("evt", Ra_eVec2, N_sz, N_sz);
		Ra_eVec2[0][0] = 1;
		if (B_prt) pverbose("EVCM %s\n", xxx.getReadable());
	} else
		sseUnmat(Ra_eVec2, N_sz);

#else /* DOES NOT WORK PROPERLY */
//	for (int k = N_sz-1 ; k>=1 ; --k) {
	for (int k=0 ; k<(N_sz-1) ; k++) {
//		wsUint _k = k-1;
		wsUint _k = k+1;
		wsReal *hK = Ra_eVec[_k];
		Ra_eVec2[k][k] = 1;
		if (hK[k] != 0.0) {
			wsReal inv = 1.0 / (Ra_e[_k] * hK[k]);
			wsReal beta = 1.0 / Ra_e[_k];
			Ra_eVec2[k][k] = 1 + beta * hK[k];
//			for (wsUint i=(wsUint)k+1 ; i<N_sz ; ++i) {
			for (wsUint i=0 ; i<k ; ++i) {
				Ra_eVec2[k][i] = beta * hK[i];
			}
//			for (int j=k+1 ; j<N_sz ; ++j) {
			for (int j=0 ; j<k ; ++j) {
				beta = 0;
//				for (wsUint i=(wsUint)k+1 ; i<N_sz ; ++i) {
				for (wsUint i=0 ; i<k ; ++i) {
					beta += Ra_eVec2[j][i] * hK[i];
				}
				beta *= inv;
				Ra_eVec2[j][k] = beta * hK[k];
				for (wsUint i=0 ; i<k ; ++i) {
//				for (wsUint i=(wsUint)k+1 ; i<N_sz ; ++i) {
					Ra_eVec2[j][i] += beta * hK[i];
				}
			}
		}
	}
	Ra_eVec2[N_sz-1][N_sz-1] = 1;
#endif
	xxx.start();
//	sseUnmat(Ra_eVec, N_sz);
	N_ret = getEvEx(Ra_eVec2, Ra_diag, Ra_e, N_sz, B_negCorrect, !Rp_eVec);
	if (B_noTranspose)
		transposeSelf(Ra_eVec2, N_sz, N_sz);
	if (B_prt) pverbose("EVEX %s\n", xxx.getReadable());

	*Bp_stat = N_ret;
	if (Rp_eVec) *Rp_eVec = Ra_eVec2;
	sseFree(Ra_e);
	return Ra_diag;
}

wsReal* EIGENDECOMPOSITION(wsSym Ra_m, int N_sz, wsMat *Ra_evec/*=NULL*/,
	char B_noTranspose/*=0*/, char B_noSort/*=0*/)
{
	return EIGENDECOMPOSITION(Ra_m, N_sz, NULL, Ra_evec, B_noTranspose, B_noSort);
}

wsReal* EIGENDECOMPOSITION(wsSym Ra_m, int N_sz, char *Bp_isFailed, 
	wsMat *Ra_evec/*=NULL*/, char B_noTranspose/*=0*/, char B_noSort/*=0*/)
{
	wsMat householderVectors = sseSym2Mat(Ra_m, N_sz);
#ifdef CHECK_SSEIGEN
	__int64	N_sse		= 0;
	__int64	N_nosse		= 0;
#endif

	wsReal *main		= sseEmptyVec(N_sz);
	wsReal *secondary	= sseEmptyVec(N_sz);
#if defined(USE_SSE) && defined(USE_SSEIGEN)
	/* Maximum ending point of SSE */
	wsUint Ne0 = getMed(N_sz);
#endif

	{
		wsReal *zz			= sseEmptyVec(N_sz);
		for (int k=0 ; k<N_sz-1 ; k++) {
			//zero-out a row and a column simultaneously
			wsReal *hK = householderVectors[k];
			main[k] = hK[k];
			wsReal xNormSqr = 0;
#if defined(USE_SSE) && defined(USE_SSEIGEN)
			/* Minimum starting point of SSE */
			wsUint Nk1 = k+1 + sseJmp-1;
			Nk1 = getMed(Nk1);
#endif

#if defined(USE_SSE) && defined(USE_SSEIGEN)
			if (Nk1 < Ne0) {
#	ifdef CHECK_SSEIGEN
				if ((k+1) > Nk1) halt("ERRx");
#	endif
				for (wsUint i=(wsUint)k+1 ; i<Nk1 ; i++) {
#	ifdef CHECK_SSEIGEN
					N_nosse++;
#	endif
					wsReal c = hK[i];
					xNormSqr += c * c;
				}
				sse_t sse_ns = sseSet(0.0);
				for (wsUint i=Nk1 ; i<Ne0 ; i+=sseJmp) {
#	ifdef CHECK_SSEIGEN
					N_sse += sseJmp;
#	endif
					sse_t sse_hki = *(sse_t *)(hK + i);
					sse_ns = sseAdd(sse_ns, sseMul(sse_hki, sse_hki));
				}
				sseAppendSum(sse_ns, xNormSqr);
				for (int i=Ne0 ; i<N_sz ; i++) {
#	ifdef CHECK_SSEIGEN
					N_nosse++;
#	endif
					wsReal c = hK[i];
					xNormSqr += c * c;
				}
#	ifdef CHECK_SSEIGEN
				wsReal xx = 0.0;
				for (int j=k+1 ; j<N_sz ; ++j) {
					wsReal c = hK[j];
					xx += c * c;
				}
				if (fabs(xx-xNormSqr) > CHECK_SSEIGEN_THR)
					halt("ERR1 [%g != %g]", xNormSqr, xx);
#	endif
			} else
#endif
			for (int j=k+1 ; j<N_sz ; ++j) {
#	ifdef CHECK_SSEIGEN
				N_nosse++;
#	endif
				wsReal c = hK[j];
				xNormSqr += c * c;
			}
			wsReal a = (hK[k + 1] > 0) ? -sqrt(xNormSqr) : sqrt(xNormSqr);
			secondary[k] = a;
			if (a != 0.0) {
				// apply Householder transform from left and right simultaneously

				hK[k + 1] -= a;
				wsReal beta = -1 / (a * hK[k + 1]);

				// compute a = beta A v, where v is the Householder vector
				// this loop is written in such a way
				//   1) only the upper triangular part of the matrix is accessed
				//   2) access is cache-friendly for a matrix stored in rows
				//Arrays.fill(z, k + 1, m, 0);
				memset(zz+k+1, 0x00, sizeof(wsReal)*(N_sz-k-1));
				for (int i=k+1 ; i<N_sz ; ++i) {
					wsReal *hI = householderVectors[i];
					wsReal hKI = hK[i];
					wsReal zI = hI[i] * hKI;

#if defined(USE_SSE) && defined(USE_SSEIGEN)
					/* Minimum starting point of SSE */
					wsUint Ni1 = i+1 + sseJmp-1;
					Ni1 = getMed(Ni1);
					if (Ni1 < Ne0) {
#	ifdef CHECK_SSEIGEN
						wsReal R_zzI = zI;
						wsReal *Ra_zz = sseVector(N_sz);
						memcpy(Ra_zz, zz, N_sz*sizeof(wsReal));
						
						for (int j=i+1 ; j<N_sz ; ++j) {
							wsReal hIJ = hI[j];
							R_zzI		+= hIJ * hK[j];
							Ra_zz[j]	+= hIJ * hKI;
						}
						if ((i+1) > Ni1) halt("ERRx");
#	endif
						for (wsUint j=(wsUint)i+1 ; j<Ni1 ; j++) {
#	ifdef CHECK_SSEIGEN
							N_nosse++;
#	endif
							wsReal hIJ = hI[j];
							zI		+= hIJ * hK[j];
							zz[j]	+= hIJ * hKI;
						}
						sse_t sse_hKI = sseSet(hKI);
						sse_t sse_zI = sseSet(0.0);
						for (wsUint j=Ni1 ; j<Ne0 ; j+=sseJmp) {
#	ifdef CHECK_SSEIGEN
							N_sse += sseJmp;
#	endif
							sse_t *sse_hij = (sse_t *)(hI + j);
							sse_t *sse_hkj = (sse_t *)(hK + j);
							sse_t *sse_zzj = (sse_t *)(zz + j);

							sse_zI = sseAdd(sse_zI, sseMul(*sse_hij, *sse_hkj));
							*sse_zzj = sseAdd(*sse_zzj, sseMul(*sse_hij, sse_hKI));
						}
						sseAppendSum(sse_zI, zI);
						for (int j=Ne0 ; j<N_sz ; j++) {
#	ifdef CHECK_SSEIGEN
							N_nosse++;
#	endif
							wsReal hIJ = hI[j];
							zI		+= hIJ * hK[j];
							zz[j]	+= hIJ * hKI;
						}
#	ifdef CHECK_SSEIGEN
						if (fabs(R_zzI-zI) > CHECK_SSEIGEN_THR)
							halt("ERR2 [%g != %g]", R_zzI, zI);
						for (int j=i+1 ; j<N_sz ; ++j) {
							if (fabs(Ra_zz[j]-zz[j]) > CHECK_SSEIGEN_THR)
								halt("ERR3 [%g != %g]", Ra_zz[j], zz[j]);
						}
						sseFree(Ra_zz);
#endif
					} else
#endif
					for (int j=i+1 ; j<N_sz ; ++j) {
#	ifdef CHECK_SSEIGEN
						N_nosse++;
#	endif
						wsReal hIJ = hI[j];
						zI		+= hIJ * hK[j];
						zz[j]	+= hIJ * hKI;
					}
					zz[i] = beta * (zz[i] + zI);
				}

				// compute gamma = beta vT z / 2
				wsReal gamma = 0;
#if defined(USE_SSE) && defined(USE_SSEIGEN)
				if (Nk1 < Ne0) {
#	ifdef CHECK_SSEIGEN
					if ((k+1) > Nk1) halt("ERRx");
#	endif
					for (wsUint i=(wsUint)k+1 ; i<Nk1 ; i++) {
#	ifdef CHECK_SSEIGEN
						N_nosse++;
#	endif
						gamma += zz[i] * hK[i];
					}

					sse_t sse_ga = sseSet(0.0);
					for (wsUint i=Nk1 ; i<Ne0 ; i+=sseJmp) {
#	ifdef CHECK_SSEIGEN
						N_sse += sseJmp;
#	endif
						sse_t *sse_zzi = (sse_t *)(zz + i);
						sse_t *sse_hki = (sse_t *)(hK + i);
						sse_ga = sseAdd(sse_ga, sseMul(*sse_zzi, *sse_hki));
					}
					sseAppendSum(sse_ga, gamma);
					for (int i=Ne0 ; i<N_sz ; i++) {
#	ifdef CHECK_SSEIGEN
						N_nosse++;
#	endif
						gamma += zz[i] * hK[i];
					}

#	ifdef CHECK_SSEIGEN
					wsReal R_gg = 0;
					for (int i=k+1 ; i<N_sz ; ++i)
						R_gg += zz[i] * hK[i];
					if (fabs(R_gg-gamma) > CHECK_SSEIGEN_THR)
						halt("ERR4 [%g != %g]", R_gg, gamma);
#	endif
				} else
#endif
				for (int i=k+1 ; i<N_sz ; ++i) {
#	ifdef CHECK_SSEIGEN
					N_nosse++;
#	endif
					gamma += zz[i] * hK[i];
				}
				gamma *= beta / 2;

				// compute z = z - gamma v
#if defined(USE_SSE) && defined(USE_SSEIGEN)
				if (Nk1 < Ne0) {
#	ifdef CHECK_SSEIGEN
					if ((k+1) > Nk1) halt("ERRx");
					wsReal *Ra_zz = sseVector(N_sz);
					memcpy(Ra_zz, zz, sizeof(wsReal)*N_sz);
					for (int i=k+1 ; i<N_sz ; ++i)
						Ra_zz[i] -= gamma * hK[i];
#	endif
					for (wsUint i=(wsUint)k+1 ; i<Nk1 ; i++) {
#	ifdef CHECK_SSEIGEN
						N_nosse++;
#	endif
						zz[i] -= gamma * hK[i];
					}

					sse_t sse_ga = sseSet(gamma);
					for (wsUint i=Nk1 ; i<Ne0 ; i+=sseJmp) {
#	ifdef CHECK_SSEIGEN
						N_sse += sseJmp;
#	endif
						sse_t *sse_zzi = (sse_t *)(zz + i);
						sse_t *sse_hki = (sse_t *)(hK + i);
						*sse_zzi = sseSub(*sse_zzi, sseMul(sse_ga, *sse_hki));
					}

					for (int i=Ne0 ; i<N_sz ; i++) {
#	ifdef CHECK_SSEIGEN
						N_nosse++;
#	endif
						zz[i] -= gamma * hK[i];
					}
#	ifdef CHECK_SSEIGEN
					for (int i=k+1 ; i<N_sz ; ++i)
						if (fabs(Ra_zz[i]-zz[i]) > CHECK_SSEIGEN_THR)
							halt("ERR5 [%g != %g]", Ra_zz[i], zz[i]);
					sseFree(Ra_zz);
#	endif
				} else
#endif
				for (int i=k+1 ; i<N_sz ; ++i) {
#	ifdef CHECK_SSEIGEN
					N_nosse++;
#	endif
					zz[i] -= gamma * hK[i];
				}

				// update matrix: A = A - v zT - z vT
				// only the upper triangular part of the matrix is updated
				for (int i=k+1 ; i<N_sz ; ++i) {
					wsReal zzi = zz[i];
					wsReal hKi = hK[i];
					wsReal *hI = householderVectors[i];

#if defined(USE_SSE) && defined(USE_SSEIGEN)
					/* Minimum starting point of SSE */
					wsUint Ni0 = i + sseJmp-1;
					Ni0 = getMed(Ni0);
					if (Ni0 < Ne0) {
#	ifdef CHECK_SSEIGEN
						if (i > Ni0) halt("ERRx");
						wsReal *Ra_hi = sseVector(N_sz);
						memcpy(Ra_hi, hI, sizeof(wsReal)*N_sz);
						for (int j=i ; j<N_sz ; ++j)
							Ra_hi[j] -= hKi * zz[j] + zzi * hK[j];
#	endif
						for (wsUint j=(wsUint)i ; j<Ni0 ; j++) {
#	ifdef CHECK_SSEIGEN
							N_nosse++;
#	endif
							hI[j] -= hKi * zz[j] + zzi * hK[j];
						}

						sse_t sse_hKI = sseSet(hKi);
						sse_t sse_zzI = sseSet(zzi);
						for (wsUint j=Ni0 ; j<Ne0 ; j+=sseJmp) {
#	ifdef CHECK_SSEIGEN
							N_sse += sseJmp;
#	endif
							sse_t *sse_hij = (sse_t *)(hI + j);
							sse_t *sse_hkj = (sse_t *)(hK + j);
							sse_t *sse_zzj = (sse_t *)(zz + j);

							*sse_hij = sseSub(*sse_hij, sseAdd(
								sseMul(sse_hKI, *sse_zzj),
								sseMul(sse_zzI, *sse_hkj)));
						}
						for (int j=Ne0 ; j<N_sz ; j++) {
#	ifdef CHECK_SSEIGEN
							N_nosse++;
#	endif
							hI[j] -= hKi * zz[j] + zzi * hK[j];
						}
#	ifdef CHECK_SSEIGEN
						for (int j=i ; j<N_sz ; ++j)
							if (fabs(Ra_hi[j]-hI[j]) > CHECK_SSEIGEN_THR)
								halt("ERR6 [%g != %g]", Ra_hi[j], hI[j]);
						sseFree(Ra_hi);
#	endif
					} else
#endif
					for (int j=i ; j<N_sz ; ++j) {
#	ifdef CHECK_SSEIGEN
						N_nosse++;
#	endif
						hI[j] -= hKi * zz[j] + zzi * hK[j];
					}
				}
			}
		}
		main[N_sz - 1] = householderVectors[N_sz - 1][N_sz - 1];
		sseFree(zz);
	}

	// build up first part of the matrix by applying Householder transforms
	wsMat z = NULL;
	if (Ra_evec) {
		z = sseEmptyMatrix(N_sz, N_sz);
		for (int k=N_sz-1 ; k>=1 ; --k) {
#if defined(USE_SSE) && defined(USE_SSEIGEN)
			/* Minimum starting point of SSE */
			wsUint N1 = k+1 + sseJmp-1;
			N1 = getMed(N1);
#endif
			wsReal *hK = householderVectors[k - 1];
			wsReal *zk = z[k];
			zk[k] = 1;
			if (hK[k] != 0.0) {
				wsReal inv = 1.0 / (secondary[k - 1] * hK[k]);
				wsReal beta = 1.0 / secondary[k - 1];

				zk[k] = 1 + beta * hK[k];
#if defined(USE_SSE) && defined(USE_SSEIGEN)
				if (N1 < Ne0) {
#	ifdef CHECK_SSEIGEN
					if ((k+1) > N1) halt("ERRx");
					wsReal *Ra_zk = sseVector(N_sz);
					for (int i=k+1 ; i<N_sz ; ++i)
						Ra_zk[i] = beta * hK[i];
#	endif
					/* k -> N1 */
					for (wsUint i=(wsUint)k+1 ; i<N1 ; i++) {
#	ifdef CHECK_SSEIGEN
						N_nosse++;
#	endif
						zk[i] = beta * hK[i];
					}

					sse_t sse_beta = sseSet(beta);
					for (wsUint i=N1 ; i<Ne0 ; i+=sseJmp) {
#	ifdef CHECK_SSEIGEN
						N_sse += sseJmp;
#	endif
						*((sse_t *)(zk + i)) = sseMul(sse_beta,
							*((sse_t *)(hK + i)));
					}

					for (int i=Ne0 ; i<N_sz ; i++) {
#	ifdef CHECK_SSEIGEN
						N_nosse++;
#	endif
						zk[i] = beta * hK[i];
					}
#	ifdef CHECK_SSEIGEN
					for (int i=k+1 ; i<N_sz ; ++i)
						if (fabs(Ra_zk[i]-zk[i]) > CHECK_SSEIGEN_THR)
							halt("ERR7 [%g != %g]", Ra_zk[i], zk[i]);
					sseFree(Ra_zk);
#	endif
				} else 
#endif
				for (int i=k+1 ; i<N_sz ; ++i) {
#	ifdef CHECK_SSEIGEN
					N_nosse++;
#	endif
					zk[i] = beta * hK[i];
				}

				for (int j=k+1 ; j<N_sz ; ++j) {
					wsReal *zj = z[j];
					beta = 0;

#if defined(USE_SSE) && defined(USE_SSEIGEN)
					if (N1 < Ne0) {
						/* k -> N1 */
						for (wsUint i=(wsUint)k+1 ; i<N1 ; i++) {
#	ifdef CHECK_SSEIGEN
							N_nosse++;
#	endif
							beta += zj[i] * hK[i];
						}

						sse_t sse_beta = sseSet(0.0);
						for (wsUint i=N1 ; i<Ne0 ; i+=sseJmp) {
#	ifdef CHECK_SSEIGEN
							N_sse += sseJmp;
#	endif
							sse_beta = sseAdd(sse_beta, sseMul(
								*((sse_t *)(zj + i)), *((sse_t *)(hK + i))));
						}
						sseAppendSum(sse_beta, beta);
						for (int i=Ne0 ; i<N_sz ; i++) {
#	ifdef CHECK_SSEIGEN
							N_nosse++;
#	endif
							beta += zj[i] * hK[i];
						}
#	ifdef CHECK_SSEIGEN
						wsReal R_gg = 0;
						for (int i=k+1 ; i<N_sz ; ++i)
							R_gg += zj[i] * hK[i];
						if (fabs(R_gg-beta) > CHECK_SSEIGEN_THR)
							halt("ERR8 [%g != %g]", R_gg, beta);
#	endif
					} else
#endif
					for (int i=k+1 ; i<N_sz ; ++i) {
#	ifdef CHECK_SSEIGEN
						N_nosse++;
#	endif
						beta += zj[i] * hK[i];
					}

					beta *= inv;
					zj[k] = beta * hK[k];

#if defined(USE_SSE) && defined(USE_SSEIGEN)
					if (N1 < Ne0) {
#	ifdef CHECK_SSEIGEN
						wsReal *Ra_zz = sseVector(N_sz);
						memcpy(Ra_zz, zj, sizeof(wsReal)*N_sz);
						for (int i=k+1 ; i<N_sz ; ++i)
							Ra_zz[i] += beta * hK[i];
#	endif
						/* k -> N1 */
						for (wsUint i=(wsUint)k+1 ; i<N1 ; i++) {
#	ifdef CHECK_SSEIGEN
							N_nosse++;
#	endif
							zj[i] += beta * hK[i];
						}

						sse_t sse_beta = sseSet(beta);
						for (wsUint i=N1 ; i<Ne0 ; i+=sseJmp) {
#	ifdef CHECK_SSEIGEN
							N_sse += sseJmp;
#	endif
							sse_t *sse_cJ = (sse_t *)(zj + i);
							*sse_cJ = sseAdd(*sse_cJ, sseMul(
								sse_beta, *((sse_t *)(hK + i))));
						}
						for (int i=Ne0 ; i<N_sz ; i++) {
#	ifdef CHECK_SSEIGEN
							N_nosse++;
#	endif
							zj[i] += beta * hK[i];
						}

#	ifdef CHECK_SSEIGEN
						for (int i=k+1 ; i<N_sz ; ++i)
							if (fabs(Ra_zz[i]-zj[i]) > CHECK_SSEIGEN_THR)
								halt("ERR9 [%g != %g]", Ra_zz[i], zj[i]);
						sseFree(Ra_zz);
#	endif
					} else
#endif
					for (int i=k+1 ; i<N_sz ; ++i) {
#	ifdef CHECK_SSEIGEN
						N_nosse++;
#	endif
						zj[i] += beta * hK[i];
					}
				}
			}
		}
		z[0][0] = 1;
	}
	sseUnmat(householderVectors, N_sz);

//	594        final double[][]z = householderMatrix.clone();
	int n = N_sz;
	wsReal *realEigenvalues = sseVector(n);
	wsReal *e = sseVector(n);
	for (int i = 0; i < n - 1; i++) {
		realEigenvalues[i] = main[i];
		e[i] = secondary[i];
	}
	realEigenvalues[n - 1] = main[n - 1];
	e[n - 1] = 0;
	sseFree(main);
	sseFree(secondary);

	// Determine the largest main and secondary value in absolute term.
	wsReal maxAbsoluteValue = 0;
	for (int i = 0; i < n; i++) {
		if (fabs(realEigenvalues[i]) > maxAbsoluteValue) {
			maxAbsoluteValue = fabs(realEigenvalues[i]);
		}
		if (fabs(e[i]) > maxAbsoluteValue) {
			maxAbsoluteValue = fabs(e[i]);
		}
	}
	// Make null any main and secondary value too small to be significant
	if (maxAbsoluteValue != 0) {
		for (int i=0; i < n; i++) {
			if (fabs(realEigenvalues[i]) <= DBL_EPSILON* maxAbsoluteValue) {
				realEigenvalues[i] = 0;
			}
			if (fabs(e[i]) <= DBL_EPSILON* maxAbsoluteValue) {
				e[i]=0;
			}
		}
	}

	for (int j = 0; j < n; j++) {
		int its = 0;
		int m;
		do {
//			wsReal R_del = WISARD_NAN;
			for (m = j; m < n - 1; m++) {
				wsReal delta = fabs(realEigenvalues[m]) +
					fabs(realEigenvalues[m + 1]);
				if (fabs(e[m]) + delta == delta) {
					break;
				}
//				R_del = e[m];
			}
			if (m != j) {
				if (its == 500) {
					if (Bp_isFailed) {
						*Bp_isFailed = 1;
						return NULL;
					} else {
						exportSymMat("eigen.failed.mat", Ra_m, N_sz);
						halt("Eigendecomposition failed of dimension [%d] matrix", n);
					}
				}
				its++;
				wsReal q = (realEigenvalues[j + 1] - realEigenvalues[j]) / (2 * e[j]);
				wsReal t = sqrt(1 + q * q);
				if (q < 0.0) {
					q = realEigenvalues[m] - realEigenvalues[j] + e[j] / (q - t);
				} else {
					q = realEigenvalues[m] - realEigenvalues[j] + e[j] / (q + t);
				}
				wsReal u = 0.0;
				wsReal s = 1.0;
				wsReal c = 1.0;
				int i;
				for (i = m - 1; i >= j; i--) {
					wsReal p = s * e[i];
					wsReal h = c * e[i];

					if (fabs(p) >= fabs(q)) {
						c = q / p;
						t = sqrt(c * c + 1.0);
						e[i + 1] = p * t;
						s = 1.0 / t;
						c = c * s;
					} else {
						s = p / q;
						t = sqrt(s * s + 1.0);
						e[i + 1] = q * t;
						c = 1.0 / t;
						s = s * c;
					}
					if (e[i + 1] == 0.0) {
						realEigenvalues[i + 1] -= u;
						e[m] = 0.0;
						break;
					}
					q = realEigenvalues[i + 1] - u;
					t = (realEigenvalues[i] - q) * s + 2.0 * c * h;
					u = s * t;
					realEigenvalues[i + 1] = q + u;
					q = c * t - h;
					if (Ra_evec) {
// 						for (int ia = 0; ia < n; ia++) {
// 							p = z[i + 1][ia];
// 							z[i + 1][ia] = s * z[i][ia] + c * p;
// 							z[i][ia] = c * z[i][ia] - s * p;
// 						}
						int N_med = 0;
#if defined(USE_SSE) && defined(USE_SSEIGEN)
						N_med = getMed(n);
						sse_t sse_s = sseSet(s);
						sse_t sse_c = sseSet(c);
#	ifdef CHECK_SSEIGEN
						wsReal *z1 = sseVector(N_sz);
						wsReal *z2 = sseVector(N_sz);
						memcpy(z1, z[i], sizeof(wsReal)*N_sz);
						memcpy(z2, z[i+1], sizeof(wsReal)*N_sz);
						for (int ia=0 ; ia<N_med ; ia++) {
							p = z2[ia];
							z2[ia] = s * z1[ia] + c * p;
							z1[ia] = c * z1[ia] - s * p;
						}
#	endif
						for (wsUint ia=0 ; ia<N_med ; ia+=sseJmp) {
							sse_t *sse_i0 = (sse_t *)(z[i] + ia);
							sse_t *sse_i1 = (sse_t *)(z[i+1] + ia);
							sse_t sse_v1 = sseLoad(z[i+1] + ia);
#	ifdef CHECK_SSEIGEN
							N_sse += sseJmp;
#	endif

							*sse_i1 = sseAdd(sseMul(sse_s, *sse_i0), sseMul(sse_c, sse_v1));
							*sse_i0 = sseSub(sseMul(sse_c, *sse_i0), sseMul(sse_s, sse_v1));
						}
#	ifdef CHECK_SSEIGEN
						for (int ia=0 ; ia<N_med ; ++ia)
							if (fabs(z1[ia]-z[i][ia]) > CHECK_SSEIGEN_THR ||
								fabs(z2[ia]-z[i+1][ia]) > CHECK_SSEIGEN_THR)
								halt("ERR10 [%g != %g] or [%g != %g]", z1[ia], z[i][ia],
									z2[ia], z[i+1][ia]);
						sseFree(z1);
						sseFree(z2);
#	endif
#endif
						for (int ia=N_med ; ia<n ; ia++) {
#	ifdef CHECK_SSEIGEN
							N_nosse++;
#	endif
							p = z[i + 1][ia];
							z[i + 1][ia] = s * z[i][ia] + c * p;
							z[i][ia] = c * z[i][ia] - s * p;
						}
					}
				}
				if (t == 0.0 && i >= j) {
					continue;
				}
				realEigenvalues[j] -= u;
				e[j] = q;
				e[m] = 0.0;
			}
		} while (m != j);
	}
	sseFree(e);

#ifdef CHECK_SSEIGEN
	LOG("SSE [" FMT_INT64 "] NoSSE[" FMT_INT64 "] Prop [%g]\n", N_sse, N_nosse, N_sse/(wsReal)N_nosse);
#endif

	if (!B_noSort) {
		//Sort the eigen values (and vectors) in increase order
		xRealSort *Xa_r = buildRealSort(realEigenvalues, n);
		qsort(Xa_r, n, sizeof(xRealSort), sort_real_desc);
		if (Ra_evec) {
			wsMat R = NULL;
			wsCalloc(R, wsReal*, n);
			for (int i = 0; i < n; i++) {
				realEigenvalues[i] = Xa_r[i].V;
				R[i] = z[Xa_r[i].i];
			}
			memcpy(z, R, sizeof(wsReal*)*n);
			DEALLOC(R);
		} else for (int i = 0; i < n; i++) {
			realEigenvalues[i] = Xa_r[i].V;
		}
		DEALLOC(Xa_r);
	}

	// Determine the largest eigen value in absolute term.
	maxAbsoluteValue = 0;
	for (int i = 0; i < n; i++) {
		if (fabs(realEigenvalues[i]) > maxAbsoluteValue) {
			maxAbsoluteValue=fabs(realEigenvalues[i]);
		}
	}
	// Make null any eigen value too small to be significant
	if (maxAbsoluteValue != 0.0) {
		for (int i=0; i < n; i++) {
			if (realEigenvalues[i] < DBL_EPSILON * maxAbsoluteValue) {
				realEigenvalues[i] = 0;
			}
		}
	}
	if (Ra_evec) {
		if (B_noTranspose)
			transposeSelf(z, n, n);
		*Ra_evec = z;
	}
	if (Bp_isFailed) *Bp_isFailed = 1;
	return realEigenvalues;
}

wsReal* eigenDecomp(wsReal **Ra_mat, wsUint N_sz, char *Bp_stat,
	wsReal ***Rp_eVec, char B_transpose/*=0*/, char B_negCorrect/*=1*/)
{
	char	N_ret;
	wsReal	*Ra_e, *Ra_diag = NULL;
	wsUint	i, j;

	sseMalloc(Ra_e, wsReal, N_sz+1);
	sseMalloc(Ra_diag, wsReal, N_sz+1);
	wsReal	**Ra_eVec;
	Ra_eVec = sseMatrixP(N_sz, N_sz, Ra_mat);

	cTimer c;
	//	c.start();
#if 0
	sseTred2(Ra_eVec, N_sz, Ra_diag, Ra_e);
#else
	tred2(Ra_eVec, N_sz, Ra_diag, Ra_e);
#endif
	//	LOG("Tred2 %s\n", c.getReadable());
	//	exportMatrix("eig.e", &Ra_e, 1, N_sz+1);
	//	exportMatrix("eig.diag", &Ra_diag, 1, N_sz);
	//	c.start();
	N_ret = tqli(Ra_diag, Ra_e, N_sz, Ra_eVec);

	/* Make order index */
	xRealSort *Xa_r = NULL;
	wsAlloc(Xa_r, xRealSort, N_sz);
	if (B_negCorrect) for (i=0 ; i<N_sz ; i++) {
		Xa_r[i].i = i;
		Xa_r[i].V = Ra_diag[i] < 0 ? 0 : Ra_diag[i];
		/* Should be all positive */
	} else for (i=0 ; i<N_sz ; i++) {
		Xa_r[i].i = i;
		Xa_r[i].V = Ra_diag[i];
		/* Should be all positive */
	}
	qsort(Xa_r, N_sz, sizeof(xRealSort), sortabs_real);
	for (i=0 ; i<N_sz ; i++)
		Ra_diag[i] = Xa_r[i].V;

	/* Re-order EVecs */
	wsReal **Ra_eVecNew = sseMatrix(N_sz, N_sz);
	if (B_transpose) {
		for (i=0 ; i<N_sz ; i++) {
			wsUint N_idx = Xa_r[i].i;
			for (j=0 ; j<N_sz ; j++)
				Ra_eVecNew[i][j] = Ra_eVec[j][N_idx];
		}
	} else {
		for (i=0 ; i<N_sz ; i++) {
			wsUint N_idx = Xa_r[i].i;
			for (j=0 ; j<N_sz ; j++)
				Ra_eVecNew[j][i] = Ra_eVec[j][N_idx];
		}
	}
	sseUnmat(Ra_eVec, N_sz);
	*Rp_eVec = Ra_eVecNew;
	DEALLOC(Xa_r);
	//	LOG("tqli1 %s\n", c.getReadable());

	//	exportMatrix("eig.e2", &Ra_e, 1, N_sz+1);
	//	exportMatrix("eig.diag2", &Ra_diag, 1, N_sz);

	*Bp_stat = N_ret;
	sseFree(Ra_e);

	return Ra_diag;
}

wsReal detCholesky(wsReal **Ra_m, wsUint N_sz)
	//	Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
	//	decomposition, A = L �� L'. On input, only the upper triangle of a need be given; it is not
	//	modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
	//	elements which are returned in p[1..n].
{
	wsUint	i, j, k;
	wsReal	R_sum, R_det;

	R_sum = Ra_m[0][0];
	if (R_sum <= W0) // a, with rounding errors, is not positive definite.
		halt("choldc failed");
	R_det = sqrt(R_sum);

	for (j=1 ; j<N_sz ; j++)
		Ra_m[j][0] = Ra_m[j][0]/R_det;

	R_det = R_sum;

	for (i=1 ; i<N_sz ; i++) {
		wsReal R_sqSum;
		R_sum = Ra_m[i][i];
#ifdef USE_SSE
		sse_t sse_s = sseSet(0.0);
		wsUint N_med = getMed(i);

		for (k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *sse_ik = (sse_t *)(Ra_m[i] + k);
			sse_s = sseSub(sse_s, sseMul(*sse_ik, *sse_ik));
		}
		sseAppendSum(sse_s, R_sum);
		for (k=N_med ; k<i ; k++)
			R_sum -= Ra_m[i][k]*Ra_m[i][k];
#else
		for (k=i-1 ; ; k--) {
			R_sum -= Ra_m[i][k]*Ra_m[i][k];
			if (k == 0) break;
		}
#endif
		if (R_sum <= W0) // a, with rounding errors, is not positive definite.
			return numeric_limits<wsReal>::quiet_NaN();
			//halt("choldc failed");

		R_sqSum	= sqrt(R_sum);
		R_det	*= R_sum;
		if (R_det == numeric_limits<wsReal>::infinity())
			halt("Too large determinant for det()");

		for (j=i+1 ; j<N_sz ; j++) {
			R_sum=Ra_m[i][j];
#ifdef USE_SSE
			sse_s = sseSet(0.0);
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_ik = (sse_t *)(Ra_m[i] + k);
				sse_t *sse_jk = (sse_t *)(Ra_m[j] + k);

				sse_s = sseSub(sse_s, sseMul(*sse_ik, *sse_jk));
			}
			sseAppendSum(sse_s, R_sum);
			for (k=N_med ; k<i ; k++)
				R_sum -= Ra_m[i][k]*Ra_m[j][k];
#else
			for (k=i-1 ; ; k--) {
				R_sum -= Ra_m[i][k]*Ra_m[j][k];
				if (k == 0) break;
			}
#endif
			Ra_m[j][i] = R_sum/R_sqSum;
		}
	}

	return R_det;
}

wsReal determinantCholesky(wsReal **Ra_m, wsUint N_sz)
	//	Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
	//	decomposition, A = L �� L'. On input, only the upper triangle of a need be given; it is not
	//	modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
	//	elements which are returned in p[1..n].
{
	wsUint	i, j, k;
	wsReal	R_sum;
	wsReal	R_det;

	R_sum = Ra_m[0][0];
	if (R_sum <= W0) // a, with rounding errors, is not positive definite.
		halt("choldc failed");
	R_det = sqrt(R_sum);

	for (j=1 ; j<N_sz ; j++)
		Ra_m[j][0] = Ra_m[j][0]/R_det;

	R_det = log(fabs(R_sum));

	for (i=1 ; i<N_sz ; i++) {
		wsReal R_sqSum;
		R_sum = Ra_m[i][i];
#ifdef USE_SSE
		sse_t sse_s = sseSet(0.0);
		wsUint N_med = getMed(i);

		for (k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *sse_ik = (sse_t *)(Ra_m[i] + k);
			sse_s = sseSub(sse_s, sseMul(*sse_ik, *sse_ik));
		}
		sseAppendSum(sse_s, R_sum);
		for (k=N_med ; k<i ; k++)
			R_sum -= Ra_m[i][k]*Ra_m[i][k];
#else
		for (k=i-1 ; ; k--) {
			R_sum -= Ra_m[i][k]*Ra_m[i][k];
			if (k == 0) break;
		}
#endif
		if (R_sum <= W0) // a, with rounding errors, is not positive definite.
			halt("choldc failed");
		R_sqSum	= sqrt(R_sum);
		R_det	+= log(fabs(R_sum));

		for (j=i+1 ; j<N_sz ; j++) {
			R_sum = Ra_m[i][j];
#ifdef USE_SSE
			sse_s = sseSet(0.0);
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_ik = (sse_t *)(Ra_m[i] + k);
				sse_t *sse_jk = (sse_t *)(Ra_m[j] + k);

				sse_s = sseSub(sse_s, sseMul(*sse_ik, *sse_jk));
			}
			sseAppendSum(sse_s, R_sum);
			for (k=N_med ; k<i ; k++)
				R_sum -= Ra_m[i][k]*Ra_m[j][k];
#else
			for (k=i-1 ; ; k--) {
				R_sum -= Ra_m[i][k]*Ra_m[j][k];
				if (k == 0) break;
			}
#endif
			Ra_m[j][i] = R_sum/R_sqSum;
		}
	}

	return R_det;
}

void chol(wsReal **Ra_m, wsUint N_sz)
	//	Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
	//	decomposition, A = L �� L'. On input, only the upper triangle of a need be given; it is not
	//	modified. The Cholesky factor L is returned in the lower triangle of a, except for its diagonal
	//	elements which are returned in p[1..n].
{
	wsUint	i, j, k;
	wsReal	R_sum;
	wsReal	*Ra_p = NULL;
	sseMalloc(Ra_p, wsReal, N_sz);

	R_sum = Ra_m[0][0];
	if (R_sum <= W0) { // a, with rounding errors, is not positive definite.
		exportMatrix("chol.fail.mat", Ra_m, N_sz, N_sz);
		halt("choldc failed");
	}
	Ra_p[0] = sqrt(R_sum);

	for (j=1 ; j<N_sz ; j++)
		Ra_m[j][0] = Ra_m[j][0]/Ra_p[0];

	for (i=1 ; i<N_sz ; i++) {
		R_sum = Ra_m[i][i];
#ifdef USE_SSE
		sse_t sse_s = sseSet(0.0);
		wsUint N_med = getMed(i);

		for (k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *sse_ik = (sse_t *)(Ra_m[i] + k);
			sse_s = sseSub(sse_s, sseMul(*sse_ik, *sse_ik));
		}
		sseAppendSum(sse_s, R_sum);
		for (k=N_med ; k<i ; k++)
			R_sum -= Ra_m[i][k]*Ra_m[i][k];
#else
		for (k=i-1 ; ; k--) {
			R_sum -= Ra_m[i][k]*Ra_m[i][k];
			if (k == 0) break;
		}
#endif
		if (R_sum <= W0) { // a, with rounding errors, is not positive definite.
			exportMatrix("chol.fail.mat", Ra_m, N_sz, N_sz);
			halt("choldc failed");
		}
		Ra_p[i]	= sqrt(R_sum);

		for (j=i+1 ; j<N_sz ; j++) {
			R_sum = Ra_m[i][j];
#ifdef USE_SSE
			sse_s = sseSet(0.0);
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_ik = (sse_t *)(Ra_m[i] + k);
				sse_t *sse_jk = (sse_t *)(Ra_m[j] + k);

				sse_s = sseSub(sse_s, sseMul(*sse_ik, *sse_jk));
			}
			sseAppendSum(sse_s, R_sum);
			for (k=N_med ; k<i ; k++)
				R_sum -= Ra_m[i][k]*Ra_m[j][k];
#else
			for (k=i-1 ; ; k--) {
				R_sum -= Ra_m[i][k]*Ra_m[j][k];
				if (k == 0) break;
			}
#endif
			Ra_m[j][i] = R_sum/Ra_p[i];
		}
	}
	for (i=0 ; i<N_sz ; i++) {
		Ra_m[i][i] = Ra_p[i];
		memset(Ra_m[i]+i+1, 0x00, sizeof(wsReal)*(N_sz-i-1));
	}
	sseFree(Ra_p);
}

/*detCholesky <- function(a, n)
{
	n <- dim(a)[2]
	p <- c()

	i <- 1
	sum <- a[i,i]
	cat("Loop ", i, " : ", sum, "\n")
	if (sum <= 0) # Ra_m, with rounding errors, is not positive definite.
		stop(paste("choldc failed", i, i, sum))
	p[i] <- sqrt(sum)

	for (j in (i+1):n) {
		sum <- a[i,j]
		a[j,i] <- sum/p[i]
	}
	for (i in 1:n)
		cat(a[i,], "\n")

	for (i in 2:n) {
		sum <- a[i,i]
		sum <- sum-sum(a[i,(i-1):1]*a[i,(i-1):1])
		cat("Loop ", i, " : ", sum, "\n")
		if (sum <= 0) # Ra_m, with rounding errors, is not positive definite.
			stop(paste("choldc failed", i, i, sum))
		p[i] <- sqrt(sum)
		
		if (i != n) {
			for (j in (i+1):n) {
				sum <- a[i,j]
				sum <- sum-sum(a[i,(i-1):1]*a[j,(i-1):1])
				a[j,i] <- sum/p[i]
			}
		}
	}
	det <- sum(log(p^2))
	det
}*/

wsReal* EVlanczos(wsReal **Ra_mat, wsUint N_sz, wsUint N_iter,
			wsUint *Np_getE, wsUint N_szE, wsReal R_eps)
{
	cTimer x;
	if (N_szE == 0xffffffff)
		N_szE = N_sz;
	wsReal *Ra_a, *Ra_b;

	wsReal *Ra_vp = NULL;	/* v_j */
	wsReal *Ra_vpp = NULL;	/* v_(j-1) */

	sseMalloc(Ra_vp, wsReal, N_sz);
	sseMalloc(Ra_vpp, wsReal, N_sz);

	sseMalloc(Ra_a, wsReal, N_iter);
	sseMalloc(Ra_b, wsReal, N_iter+1);
	Ra_b[0] = W0;

	memset(Ra_vpp, 0x00, sizeof(wsReal)*N_sz);
	sseInit(Ra_vp, N_sz, W1/(wsReal)sqrt((double)N_sz));
	/* Tridiagonalize */
	x.start();
	for (wsUint i=0 ; i<N_iter ; i++) {
		/*    v <- list()
		v[[1]] <- matrix(0, nrow=nr, ncol=1)    # define v_0 <- 0
		v[[2]] <- matrix(1, nrow=nr, ncol=1)
		v[[2]] <- v[[2]] / sqrt(sum(v[[2]]^2))    # define v_1, ||v_1|| = 1

		# v_i == v_(j+1) because R cannot define v[[0]]

		b <- c(0)    # Adjacent diagonal element in tridiagonal
		a <- c(0)    # Diagonal element in tridiagonal

		w <- list()
		for (j in 1:ni) {*/

		/* Av <- A%*%vp */
		wsReal *Ra_Av	= sseMpV(Ra_mat, N_sz, N_sz, Ra_vp);
		wsReal *Ra_w	= Ra_Av;
		/* Bv <- b[j]*vpp */
		wsReal *Ra_Bv	= Ra_vpp;
		sseVpC(Ra_vpp, Ra_b[i], Ra_Bv, N_sz);
		/* vpp <- vp */
		memcpy(Ra_vpp, Ra_vp, sizeof(wsReal)*N_sz);

		/* w[[j]] <- Av - Bv    # w_j <- A %*% v_j  - b_j v_(j-1) */
		sseVsV(Ra_Av, Ra_Bv, Ra_w, N_sz);
		/* a[j] <- w[[j]] (dot) vp           # a_j <- w_j^T %*% v_j */
		wsReal *Ra_ai = sseMpV(&Ra_w, 1, N_sz, Ra_vp);
		Ra_a[i] = Ra_ai[0];
		/* av <- a[j] * vp */
		sseVpC(Ra_vp, Ra_ai[0], Ra_vp, N_sz);
		/* w[[j]] <- w[[j]] - av       # w_j <- w_j - a_j * v_j */
		sseVsV(Ra_w, Ra_vp, Ra_w, N_sz);

		/* b[j+1] <- sqrt(sum(w[[j]]^2))        # b_(j+1) <- || w_j || */
		wsReal R_b = W0;
		for (wsUint j=0 ; j<N_sz ; j++)
			R_b += Ra_Av[j] * Ra_Av[j];
		R_b = sqrt(R_b);

		/* v[[j+2]] <- w[[j]] / b[j+1]        # v_(j+1) <- w_j / b_(j+1) */
		for (wsUint j=0 ; j<N_sz ; j++)
			Ra_vp[j] = Ra_Av[j] / R_b;
		Ra_b[i+1] = R_b;
		sseFree(Ra_Av);
		sseFree(Ra_ai);
	}
	//LOG("EVlanczos loop %s\n", x.getReadable());
	sseFree(Ra_vp);
	sseFree(Ra_vpp);
	/* Get eigenvalues */
	x.start();
	tqli1(N_iter, Ra_a, Ra_b);
	//LOG("EVlanczos tqli1 %s\n", x.getReadable());
	sseFree(Ra_b);

	/* Tryout eigenvalues */
	wsReal *Rp_v = NULL;
	char B_p = 0;

	wsReal *Ra_eig = NULL;
	sseMalloc(Ra_eig, wsReal, N_szE);
	wsUint N_idxEig = 0;

	for (wsUint i=0 ; i<N_iter ; i++) {
		if (Rp_v == NULL) {
			Rp_v = Ra_a + i;
			B_p = 0;
		} else {
			if ((*Rp_v-Ra_a[i]) < R_eps) {
				if (B_p == 0) {
					printf("%g ", *Rp_v);
					Ra_eig[N_idxEig++] = *Rp_v;
					if (N_idxEig == N_szE)
						break;
					B_p = 1;
				}
			} else {
				B_p = 0;
			}
		}
		Rp_v = Ra_a+i;
	}
	//LOG("Extracted %d evs\n", N_idxEig);
	*Np_getE = N_idxEig;
	sseFree(Ra_a);
	return Ra_eig;
	/* 	E <- eigen(M)$values
	v <- NULL
	p <- F
	r <- c()
	for (i in 1:length(E)) {
	if (is.null(v)) {
	v <- E[i]
	p <- F
	} else {
	if ((v-E[i]) < eps) {
	if (p == F) {
	r <- c(r,v)
	p <- T
	}
	} else {
	p <- F
	#			cat(v,"\n")
	}
	v <- E[i]
	}
	}
	r*/
}

wsUint* kmeans(wsReal **Ra_val, wsUintCst N_dim, wsUintCst N_szVal, wsUintCst N_clus)
{
	wsReal	**Ra_cent	= sseEmptyMatrix(N_clus, N_dim);
	wsUint	i, j, k;
	wsUint	*Na_pt2cl	= NULL;
	vector<wsUint>
		*Xa_cl2pt	= NULL;
	bool move;
	bool some_point_is_moving = true;
//	unsigned int num_iterations = 0;

	wsAlloc(Xa_cl2pt, vector<wsUint>, N_clus);
	sseMalloc(Na_pt2cl, wsUint, N_szVal);

	//
	// Initial partition of points
	//
	for (i=0 ; i<N_szVal ; i++) {
		wsUint N_idClus = i%N_clus;
		Na_pt2cl[i] = N_idClus;
		Xa_cl2pt[N_idClus].push_back(i);
	}

	//
	// Until not converge
	//
	while (some_point_is_moving) {
		some_point_is_moving = false;

		// For each centroid
		for (i=0 ; i<N_clus ; i++) {
			for (j=0 ; j<Xa_cl2pt[i].size() ; j++) {
				wsReal *Rp_val = Ra_val[Xa_cl2pt[i][j]];

				for (k=0 ; k<N_dim ; k++)
					Ra_cent[j][k] += Rp_val[k];
			}
			//
			// if no point in the clusters, this goes to inf (correct!)
			//
			for (j=0 ; j<N_dim ; j++)
				Ra_cent[i][j] /= (wsReal)(Xa_cl2pt[i].size());
		}
		//      std::cout << "Centroids" << std::endl << centroids__;      

		//
		// for each point
		//
		for (i=0 ; i<N_szVal ; i++) {
			// distance from current cluster
			wsUint N_myClus = Na_pt2cl[i];
			wsUint N_toClus = 0xffffffff;
			wsReal R_dist = W0;

			for (k=0 ; k<N_dim ; k++)
				R_dist += SQR(Ra_cent[N_myClus][k]-Ra_val[i][k]);
			//
			// foreach centroid
			//
			move = false;
			for (j=0 ; j<N_clus ; j++) {
				wsReal R_nDist = W0;

				for (k=0 ; k<N_dim ; k++)
					R_nDist += SQR(Ra_cent[j][k]-Ra_val[i][k]);

				if (R_nDist < R_dist){
					R_dist = R_nDist;
					move = true;
					N_toClus = j;

					// remove from current cluster
					FOREACH (vector<wsUint>::iterator, Xa_cl2pt[N_myClus],
						it) {
							if (*it == i) {
								Xa_cl2pt[N_myClus].erase(it);
								break;
							}
					}

					some_point_is_moving = true;
				}
			}

			//
			// move towards a closer centroid 
			//
			if (move) {
				Na_pt2cl[i] = N_toClus;
				Xa_cl2pt[N_toClus].push_back(i);
				// insert
				//				std::cout << "\t\tmove to cluster=" << to_cluster << std::endl;
			}
		}      
		//		num_iterations++;
	} // end while (some_point_is_moving)

	sseUnmat(Ra_cent, N_clus);
	DEALLOC(Xa_cl2pt);

	//std::cout << std::endl << "Final clusters" << std::endl;
	//std::cout << clusters_to_points__;
	return Na_pt2cl;
}

wsUint* kmeans(wsReal *Ra_val, wsUintCst N_szVal, wsUintCst N_clus)
{
	return kmeans(&Ra_val, 1, N_szVal, N_clus);
}

double Brent_fmin(double ax, double bx, double (*f)(double, void *),
				  void *info, double tol)
{
	/*  c is the squared inverse of the golden ratio */
	const double c = (3. - sqrt(5.)) * .5;

	/* Local variables */
	double a, b, d, e, p, q, r, u, v, w, x;
	double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;

	/*  eps is approximately the square root of the relative machine precision. */
	eps = DBL_EPSILON;
	tol1 = eps + 1.;/* the smallest 1.000... > 1 */
	eps = sqrt(eps);

	a = ax;
	b = bx;
	v = a + c * (b - a);
	w = v;
	x = v;

	d = 0.;/* -Wall */
	e = 0.;
	fx = (*f)(x, info);
	fv = fx;
	fw = fx;
	tol3 = tol / 3.;

	/*  main loop starts here ----------------------------------- */

	for(;;) {
		xm = (a + b) * .5;
		tol1 = eps * fabs(x) + tol3;
		t2 = tol1 * 2.;

		/* check stopping criterion */

		if (fabs(x - xm) <= t2 - (b - a) * .5) break;
		p = 0.;
		q = 0.;
		r = 0.;
		if (fabs(e) > tol1) { /* fit parabola */

			r = (x - w) * (fx - fv);
			q = (x - v) * (fx - fw);
			p = (x - v) * q - (x - w) * r;
			q = (q - r) * 2.;
			if (q > 0.) p = -p; else q = -q;
			r = e;
			e = d;
		}

		if (fabs(p) >= fabs(q * .5 * r) ||
			p <= q * (a - x) || p >= q * (b - x)) { /* a golden-section step */

				if (x < xm) e = b - x; else e = a - x;
				d = c * e;
		} else { /* a parabolic-interpolation step */

			d = p / q;
			u = x + d;

			/* f must not be evaluated too close to ax or bx */

			if (u - a < t2 || b - u < t2) {
				d = tol1;
				if (x >= xm) d = -d;
			}
		}

		/* f must not be evaluated too close to x */

		if (fabs(d) >= tol1)
			u = x + d;
		else if (d > 0.)
			u = x + tol1;
		else
			u = x - tol1;

		fu = (*f)(u, info);

		/*  update  a, b, v, w, and x */
//		if(1) {
//			LOG(" tol  = %g, u %g, v %g\n", tol, u, v);
//			LOG(" %g <= %g\n", fabs(x - xm), t2 - (b - a) * .5);
//		}


		if (fu <= fx) {
			if (u < x) b = x; else a = x;
			v = w;    w = x;   x = u;
			fv = fw; fw = fx; fx = fu;
		} else {
			if (u < x) a = u; else b = u;
			if (fu <= fw || w == x) {
				v = w; fv = fw;
				w = u; fw = fu;
			} else if (fu <= fv || v == x || v == w) {
				v = u; fv = fu;
			}
		}
	}
	/* end of main loop */

	return x;
}

char** loadStringValues(wsStrCst Sp_val, wsUint *Np_var)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev		= 0;
	char		**Sp_ret	= NULL;
	if (Sp_val == NULL) {
		*Np_var = 0;
		return NULL;
	}
	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') {
		*Np_var = 0;
		return NULL;
	}

	/* Get the number of given prevalences */
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		N_prev++;
	}
	// Since --mqls now requires the number of prevalence>2, it removed (130423)
	// 	if (N_prev > 2)
	// 		halt("Too much of prevalences (>2)");

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
	wsAlloc(Sp_ret, char*, N_prev);
	N_prev = 0;
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		if (b) *b = '\0';
		Sp_ret[N_prev++] = strdup(a+1);
	}
	free(Sp_prev);

	*Np_var = N_prev;
	return Sp_ret;
}

char** loadStringValues2(wsStrCst Sp_val, wsUint *Np_var, char S_sep=',')
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev		= 0;
	char		**Sp_ret	= NULL;
	if (Sp_val == NULL) {
		*Np_var = 0;
		return NULL;
	}
	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') {
		*Np_var = 0;
		return NULL;
	}

	/* Get the number of given prevalences */
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, S_sep);
		N_prev++;
	}
	// Since --mqls now requires the number of prevalence>2, it removed (130423)
	// 	if (N_prev > 2)
	// 		halt("Too much of prevalences (>2)");

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
	wsAlloc(Sp_ret, char*, N_prev);
	N_prev = 0;
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, S_sep);
		if (b) *b = '\0';
		Sp_ret[N_prev++] = strdup(a+1);
	}
	free(Sp_prev);

	*Np_var = N_prev;
	return Sp_ret;
}

char** loadStringValuesByWS(wsStrCst Sp_val, wsUint *Np_var)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev		= 0;
	char		**Sp_ret	= NULL;
	if (Sp_val == NULL) {
		*Np_var = 0;
		return NULL;
	}
	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') {
		*Np_var = 0;
		return NULL;
	}

	/* Get the number of given prevalences */
	for (a=Sp_prev ; a && a[0] ; a=b) {
		getString(&a, &b);
		N_prev++;
	}
	wsUint N_exp = N_prev;
	free(Sp_prev);
	Sp_prev = strdup(Sp_val);

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
	wsAlloc(Sp_ret, char*, N_prev);
	N_prev = 0;
	for (a=Sp_prev ; a && a[0] ; a=b) {
		getString(&a, &b);
		Sp_ret[N_prev++] = strdup(a);
	}
	if (N_exp != N_prev)
		halt("This function have an error! [%s]", Sp_val);
	free(Sp_prev);

	*Np_var = N_prev;
	return Sp_ret;
}

void loadGroupValues(wsStrCst Sp_val, mStrReal &Xm_map)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev		= 0;
	if (Sp_val == NULL)
		return;

	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0')
		return;

	/* Get the number of given prevalences */
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		N_prev++;
	}

	/* Process prevalence string into real */
	N_prev = 0;
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		if (b) *b = '\0';
		/* Find = */
		char *Sp_val = strchr(a+1, '=');
		if (Sp_val == NULL)
			halt("[%d]th value [%s] does not have equal sign(=)!", N_prev+1, a+1);
		*(Sp_val++) = '\0';
		Xm_map[a+1] = (wsReal)atof(a+1);
	}
}

void loadFlags(wsStrCst Sp_val, mStr &Xm_map)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev		= 0;
	if (Sp_val == NULL)
		return;

	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0')
		return;

	/* Get the number of given prevalences */
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		N_prev++;
	}

	/* Process prevalence string into real */
	N_prev = 0;
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		if (b) *b = '\0';
		/* Find = */
		char *Sp_val = strchr(a+1, '=');
		if (Sp_val == NULL) {
			Xm_map[a+1].assign("");
		} else {
			*(Sp_val++) = '\0';
			Xm_map[a+1].assign(Sp_val);
		}
	}

	free(Sp_prev);
}

wsVec loadRealValues(wsStrCst Sp_val, wsUint *Np_var)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev	= 0;
	wsVec		Ra_prev	= NULL;
	if (Sp_val == NULL) {
		*Np_var = 0;
		return NULL;
	}
	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') {
		*Np_var = 0;
		return NULL;
	}

	/* Get the number of given prevalences */
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		N_prev++;
	}
	// Since --mqls now requires the number of prevalence>2, it removed (130423)
// 	if (N_prev > 2)
// 		halt("Too much of prevalences (>2)");

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
	Ra_prev = sseVector(N_prev);
	N_prev = 0;
	for (a=Sp_prev-1 ; a ; a=b) {
		char *Sp_shouldNull = NULL;

		b = strchr(a+1, ',');
		if (b) *b = '\0';
		wsReal R_val = (wsReal)str2dbl(a+1, &Sp_shouldNull);
		/* If invalid character */
		if (Sp_shouldNull && Sp_shouldNull[0]) {
			*Np_var = 0;
			sseFree(Ra_prev);
			return NULL;
		}
		Ra_prev[N_prev++] = (wsReal)R_val;
	}
	free(Sp_prev);

	*Np_var = N_prev;
	return Ra_prev;
}
 
wsVec loadRealValuesByWS(wsStrCst Sp_val, wsUint *Np_var)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev	= 0;
	wsVec		Ra_prev	= NULL;
	if (Sp_val == NULL) {
		*Np_var = 0;
		return NULL;
	}
	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') {
		*Np_var = 0;
		return NULL;
	}

	/* Get the number of given prevalences */
	for (a=Sp_prev ; a && a[0] ; a=b) {
		getString(&a, &b);
		N_prev++;
	}
	wsUint N_exp = N_prev;
	free(Sp_prev);
	Sp_prev = strdup(Sp_val);

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
	Ra_prev = sseVector(N_prev);
	N_prev = 0;
	for (a=Sp_prev ; a && a[0] ; a=b) {
		char *Sp_shouldNull = NULL;

		getString(&a, &b);
		wsReal R_val = (wsReal)str2dbl(a, &Sp_shouldNull);
		/* If invalid character */
		if (Sp_shouldNull && Sp_shouldNull[0]) {
			*Np_var = 0;
			sseFree(Ra_prev);
			return NULL;
		}
		Ra_prev[N_prev++] = (wsReal)R_val;
	}
	if (N_exp != N_prev)
		halt("This function have an error!");
	free(Sp_prev);

	*Np_var = N_prev;
	return Ra_prev;
}

wsUint* loadIntValues(char *Sp_val, wsUint *Np_var)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev	= 0;
	wsUint*		Na_prev	= NULL;
	if (Sp_val == NULL) {
		*Np_var = 0;
		return NULL;
	}
	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') {
		*Np_var = 0;
		return NULL;
	}

	/* Get the number of given prevalences */
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		N_prev++;
	}
	// Since --mqls now requires the number of prevalence>2, it removed (130423)
	// 	if (N_prev > 2)
	// 		halt("Too much of prevalences (>2)");

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
	sseMalloc(Na_prev, wsUint, N_prev);
	N_prev = 0;
	for (a=Sp_prev-1 ; a ; a=b) {
		char *Sp_shouldNull = NULL;

		b = strchr(a+1, ',');
		if (b) *b = '\0';
		wsUint N_val = (wsUint)strtol(a+1, &Sp_shouldNull, 10);
		/* If invalid character */
		if (Sp_shouldNull && Sp_shouldNull[0]) {
			*Np_var = 0;
			sseFree(Na_prev);
			return NULL;
		}
		Na_prev[N_prev++] = N_val;
	}

	*Np_var = N_prev;
	return Na_prev;
}

wsUint* loadIntValuesByWS(char* Sp_val, wsUint *Np_var)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev	= 0;
	wsUint*		Na_prev	= NULL;
	if (Sp_val == NULL) {
		*Np_var = 0;
		return NULL;
	}
	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') {
		*Np_var = 0;
		return NULL;
	}

	/* Get the number of given prevalences */
	for (a=Sp_prev ; a && a[0] ; a=b) {
		getString(&a, &b);
		N_prev++;
	}
	wsUint N_exp = N_prev;
	free(Sp_prev);
	Sp_prev = strdup(Sp_val);

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
	sseMalloc(Na_prev, wsUint, N_prev);
	N_prev = 0;
	for (a=Sp_prev ; a && a[0] ; a=b) {
		char *Sp_shouldNull = NULL;

		getString(&a, &b);
		wsUint N_val = strtol(a, &Sp_shouldNull, 10);
		/* If invalid character */
		if (Sp_shouldNull && Sp_shouldNull[0]) {
			*Np_var = 0;
			sseFree(Na_prev);
			return NULL;
		}
		Na_prev[N_prev++] = N_val;
	}
	if (N_exp != N_prev)
		halt("This function have an error! Expect [%d] get [%d] [%s]", N_exp, N_prev, Sp_val);
	free(Sp_prev);

	*Np_var = N_prev;
	return Na_prev;
}

/* Calculate L %*% diag(D) %*% t(L) */
/* # R code
 * # lt <- lower triangular matrix
 * # d <- diagonal term
 * LDLt <- function(lt, d) {
 *   n <- length(d)
 *   res <- matrix(NA, nrow=dim(lt)[1], ncol=dim(lt)[2])
 *   for (i in 1:n) {
 *     intm <- d[1:i]*lt[i,1:i]
 *     for (j in 1:i) {
 *       res[i,j] <- sum(intm * lt[j,1:i])
 *       res[j,i] <- res[i,j]
 *     }
 *   }
 *   res
 * }
 */
void sseLDLt(wsReal **Ra_mLowTri, wsUint N_sz, wsReal *Ra_diag,
	wsReal **Ra_dest)
{
	wsReal *Ra_intm = NULL;
	sseMalloc(Ra_intm, wsReal, N_sz);
	// for (i in 1:n) {
	for (wsUint i=0 ; i<N_sz ; i++) {
		// b_k * g_ik
		// intm <- d[1:i]*lt[i,1:i]
		sseVpV(Ra_diag, Ra_mLowTri[i], Ra_intm, i+1);
// 		for (wsUint j=0 ; j<=i ; j++)
// 			Ra_intm[j] = Ra_diag[j]*Ra_mLowTri[i][j];

		// for (j in 1:i) {
		for (wsUint j=0 ; j<=i ; j++) {
			// res[i,j] <- intm * lt[j,1:1]
			Ra_dest[i][j] = sseVV(Ra_intm, i+1, Ra_mLowTri[j]);

// 			Ra_dest[i][j] = W0;
// 			for (wsUint k=0 ; k<=i ; k++)
// 				Ra_dest[i][j] += Ra_intm[k]*Ra_mLowTri[j][k];

			/* Since it is symmetric */
			// res[j,i] <- res[i,j]
			Ra_dest[j][i] = Ra_dest[i][j];
		}
	}
	sseFree(Ra_intm);
}
void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p);

double pnorm(double x, double mu, double sigma, int lower_tail, int log_p)
{
    double p, cp;

    /* Note: The structure of these checks has been carefully thought through.
     * For example, if x == mu and sigma == 0, we get the correct answer 1.
     */
#ifdef IEEE_754
    if(ISNAN(x) || ISNAN(mu) || ISNAN(sigma))
	return x + mu + sigma;
#endif
    if(x==numeric_limits<double>::infinity() && mu == x)
		return numeric_limits<double>::quiet_NaN();/* x-mu is NaN */
    if (sigma <= 0) {
	if(sigma < 0) return numeric_limits<double>::signaling_NaN();
	/* sigma = 0 : */
	return (x < mu) ? 0.0 : 1.0;
    }
    p = (x - mu) / sigma;
    if(p == numeric_limits<double>::infinity())
	return (x < mu) ? 0.0 : 1.0;
    x = p;

    pnorm_both(x, &p, &cp, (lower_tail ? 0 : 1), log_p);

    return(lower_tail ? p : cp);
}

double truncate(double x)
{
	if(x >= 0) return floor(x);
	else return ceil(x);
}

#define SIXTEN	16 /* Cutoff allowing exact "*" and "/" */

void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p)
{
/* i_tail in {0,1,2} means: "lower", "upper", or "both" :
   if(lower) return  *cum := P[X <= x]
   if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
*/
    const static double a[5] = {
	2.2352520354606839287,
	161.02823106855587881,
	1067.6894854603709582,
	18154.981253343561249,
	0.065682337918207449113
    };
    const static double b[4] = {
	47.20258190468824187,
	976.09855173777669322,
	10260.932208618978205,
	45507.789335026729956
    };
    const static double c[9] = {
	0.39894151208813466764,
	8.8831497943883759412,
	93.506656132177855979,
	597.27027639480026226,
	2494.5375852903726711,
	6848.1904505362823326,
	11602.651437647350124,
	9842.7148383839780218,
	1.0765576773720192317e-8
    };
    const static double d[8] = {
	22.266688044328115691,
	235.38790178262499861,
	1519.377599407554805,
	6485.558298266760755,
	18615.571640885098091,
	34900.952721145977266,
	38912.003286093271411,
	19685.429676859990727
    };
    const static double p[6] = {
	0.21589853405795699,
	0.1274011611602473639,
	0.022235277870649807,
	0.001421619193227893466,
	2.9112874951168792e-5,
	0.02307344176494017303
    };
    const static double q[5] = {
	1.28426009614491121,
	0.468238212480865118,
	0.0659881378689285515,
	0.00378239633202758244,
	7.29751555083966205e-5
    };
#define M_SQRT_32	5.656854249492380195206754896838
#define ONE_SQRT_2PI 0.398942280401432677939946059934

    double xden, xnum, temp, del, eps, xsq, y;
#ifdef NO_DENORMS
    double min = DBL_MIN;
#endif
    int i, lower, upper;

#ifdef IEEE_754
    if(ISNAN(x)) { *cum = *ccum = x; return; }
#endif

    /* Consider changing these : */
    eps = DBL_EPSILON * 0.5;

    /* i_tail in {0,1,2} =^= {lower, upper, both} */
    lower = i_tail != 1;
    upper = i_tail != 0;

    y = fabs(x);
    if (y <= 0.67448975) { /* qnorm(3/4) = .6744.... -- earlier had 0.66291 */
	if (y > eps) {
	    xsq = x * x;
	    xnum = a[4] * xsq;
	    xden = xsq;
	    for (i = 0; i < 3; ++i) {
		xnum = (xnum + a[i]) * xsq;
		xden = (xden + b[i]) * xsq;
	    }
	} else xnum = xden = 0.0;

	temp = x * (xnum + a[3]) / (xden + b[3]);
	if(lower)  *cum = 0.5 + temp;
	if(upper) *ccum = 0.5 - temp;
	if(log_p) {
	    if(lower)  *cum = log(*cum);
	    if(upper) *ccum = log(*ccum);
	}
    }
    else if (y <= M_SQRT_32) {

	/* Evaluate pnorm for 0.674.. = qnorm(3/4) < |x| <= sqrt(32) ~= 5.657 */

	xnum = c[8] * y;
	xden = y;
	for (i = 0; i < 7; ++i) {
	    xnum = (xnum + c[i]) * y;
	    xden = (xden + d[i]) * y;
	}
	temp = (xnum + c[7]) / (xden + d[7]);

#define do_del(X)							\
	xsq = truncate(X * SIXTEN) / SIXTEN;				\
	del = (X - xsq) * (X + xsq);					\
	if(log_p) {							\
	    *cum = (-xsq * xsq * 0.5) + (-del * 0.5) + log(temp);	\
	    if((lower && x > 0.) || (upper && x <= 0.))			\
		  *ccum = util_log1p(-exp(-xsq * xsq * 0.5) *		\
				exp(-del * 0.5) * temp);		\
	}								\
	else {								\
	    *cum = exp(-xsq * xsq * 0.5) * exp(-del * 0.5) * temp;	\
	    *ccum = 1.0 - *cum;						\
	}

#define swap_tail						\
	if (x > 0.) {/* swap  ccum <--> cum */			\
	    temp = *cum; if(lower) *cum = *ccum; *ccum = temp;	\
	}

	do_del(y);
	swap_tail;
    }

/* else	  |x| > sqrt(32) = 5.657 :
 * the next two case differentiations were really for lower=T, log=F
 * Particularly	 *not*	for  log_p !

 * Cody had (-37.5193 < x  &&  x < 8.2924) ; R originally had y < 50
 *
 * Note that we do want symmetry(0), lower/upper -> hence use y
 */
    else if((log_p && y < 1e170) /* avoid underflow below */
	/*  ^^^^^ MM FIXME: can speedup for log_p and much larger |x| !
	 * Then, make use of  Abramowitz & Stegun, 26.2.13, something like

	 xsq = x*x;

	 if(xsq * DBL_EPSILON < 1.)
	    del = (1. - (1. - 5./(xsq+6.)) / (xsq+4.)) / (xsq+2.);
	 else
	    del = 0.;
	 *cum = -.5*xsq - M_LN_SQRT_2PI - log(x) + util_log1p(-del);
	 *ccum = util_log1p(-exp(*cum)); /.* ~ log(1) = 0 *./

 	 swap_tail;

	 [Yes, but xsq might be infinite.]

	*/
	    || (lower && -37.5193 < x  &&  x < 8.2924)
	    || (upper && -8.2924  < x  &&  x < 37.5193)
	) {

	/* Evaluate pnorm for x in (-37.5, -5.657) union (5.657, 37.5) */
	xsq = 1.0 / (x * x); /* (1./x)*(1./x) might be better */
	xnum = p[5] * xsq;
	xden = xsq;
	for (i = 0; i < 4; ++i) {
	    xnum = (xnum + p[i]) * xsq;
	    xden = (xden + q[i]) * xsq;
	}
	temp = xsq * (xnum + p[4]) / (xden + q[4]);
	temp = (ONE_SQRT_2PI - temp) / y;

	do_del(x);
	swap_tail;
    } else { /* large x such that probs are 0 or 1 */
	if(x > 0) {	*cum = 1.0; *ccum = 0.0;	}
	else {	        *cum = 0.0; *ccum = 1.0;	}
    }


#ifdef NO_DENORMS
    /* do not return "denormalized" -- we do in R */
    if(log_p) {
	if(*cum > -min)	 *cum = -0.;
	if(*ccum > -min)*ccum = -0.;
    }
    else {
	if(*cum < min)	 *cum = 0.;
	if(*ccum < min)	*ccum = 0.;
    }
#endif
    return;
}

wsMat multMPt(wsMat Ra_mat, wsUint N_row, wsUint N_col, wsReal *Rp_data,
	wsUint N_pr, wsUint N_pc, wsUint *Na_rEnd, wsUint *Na_cIdx)
{
	if (N_col != N_pc)
		halt("Column size of first matrix [%d] is not match with column size "
			"of second sparase matrix [%d]", N_col, N_pc);

	wsMat		Ra_ret	= sseEmptyMatrix(N_row, N_pr);
	wsUint		i		= 0;
	wsUint		N_s		= 0;
	wsUint		N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e) {
			/* Sum up them */
			for (wsUint r=0 ; r<N_row ; r++)
				for (wsUint j=N_s ; j<N_e ; j++) {
					wsUint TC = Na_cIdx[j];
					Ra_ret[r][i] += Ra_mat[r][TC] * Rp_data[j];
				}
		}

		/* Update start point */
		i++;
		if (i < N_pr)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	return Ra_ret;
}

wsMat multMtP(wsMat Ra_mat, wsUint N_row, wsUint N_col, wsReal *Rp_data,
	wsUint N_pr, wsUint N_pc, wsUint *Na_rEnd, wsUint *Na_cIdx)
{
	if (N_row != N_pr)
		halt("Row size of first matrix [%d] is not match with row size "
			"of second sparase matrix [%d]", N_row, N_pr);

	wsMat		Ra_ret	= sseEmptyMatrix(N_col, N_pc);
	wsUint		i		= 0;
	wsUint		N_s		= 0;
	wsUint		N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e) {
			/* Sum up them */
			for (wsUint r=0 ; r<N_col ; r++)
				for (wsUint j=N_s ; j<N_e ; j++) {
					wsUint TC = Na_cIdx[j];
					Ra_ret[TC][r] += Ra_mat[i][r] * Rp_data[j];
				}
		}

		/* Update start point */
		i++;
		if (i < N_pr)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	return Ra_ret;
}

// wsReal** sseSD(SYM_t Ra_symMat, wsUint N_s1, DIAG_t Ra_diagMat)
// {
// 	wsUint i, j;
// 	wsReal **Ra_ret = sseMatrix(N_s1, N_s1);
// 
// #if 1
// 	for (i=0 ; i<N_s1 ; i++) {
// 		wsUint N_med = 0;
// #ifdef USE_SSE
// 		N_med = getMed(i);
// 		sse_t sse_di = sseSet(Ra_diagMat[i]);
// 		for (j=0 ; j<N_med ; j+=sseJmp) {
// 			sse_t *sse_sij = (sse_t *)(Ra_symMat[i] + j);
// 			sse_t *sse_rij = (sse_t *)(Ra_ret[i] + j);
// 
// 			*sse_rij = sseMul(*sse_sij, *(sse_t *)(Ra_diagMat + j));
// //			Ra_ret[i][j] = Ra_symMat[i][j] * Ra_diagMat[j];
// 			sse_t sse_rji = sseMul(*sse_sij, sse_di);
// 			Ra_ret[j][i] = sse_rji.m128d_f64[0];
// 			Ra_ret[j+1][i] = sse_rji.m128d_f64[1];
// //			Ra_ret[j][i] = Ra_symMat[i][j] * Ra_diagMat[i];
// 		}
// #endif
// 		for (j=N_med ; j<i ; j++) {
// 			Ra_ret[i][j] = Ra_symMat[i][j] * Ra_diagMat[j];
// 			Ra_ret[j][i] = Ra_symMat[i][j] * Ra_diagMat[i];
// 		}
// 		Ra_ret[i][i] = Ra_symMat[i][i] * Ra_diagMat[i];
// 	}
// #else
// 	for (i=0 ; i<N_s1 ; i++) {
// 		for (j=0 ; j<i ; j++) {
// 			Ra_ret[i][j] = Ra_symMat[i][j] * Ra_diagMat[j];
// 		}
// 		Ra_ret[i][i] = Ra_symMat[i][i] * Ra_diagMat[i];
// 		for (j=i+1 ; j<N_s1 ; j++)
// 			Ra_ret[i][j] = Ra_symMat[j][i] * Ra_diagMat[j];
// 	}
// #endif
// 
// 	return Ra_ret;
// }
// 
// wsReal** sseDS(DIAG_t Ra_diagMat, wsUint N_s1, SYM_t Ra_symMat)
// {
// 	wsUint i, j;
// 	wsReal **Ra_ret = sseMatrix(N_s1, N_s1);
// 
// #if 1
// 	for (i=0 ; i<N_s1 ; i++) {
// 		wsUint N_med = 0;
// #ifdef USE_SSE
// 		N_med = getMed(i);
// 		sse_t sse_di = sseSet(Ra_diagMat[i]);
// 		for (j=0 ; j<N_med ; j+=sseJmp) {
// 			sse_t *sse_sij = (sse_t *)(Ra_symMat[i] + j);
// 			sse_t *sse_rij = (sse_t *)(Ra_ret[i] + j);
// 
// 			*sse_rij = sseMul(*sse_sij, sse_di);
// 			//			Ra_ret[i][j] = Ra_symMat[i][j] * Ra_diagMat[j];
// 			sse_t sse_rji = sseMul(*sse_sij, *(sse_t *)(Ra_diagMat + j));
// 			Ra_ret[j][i] = sse_rji.m128d_f64[0];
// 			Ra_ret[j+1][i] = sse_rji.m128d_f64[1];
// 			//			Ra_ret[j][i] = Ra_symMat[i][j] * Ra_diagMat[i];
// 		}
// #endif
// 		for (j=N_med ; j<i ; j++) {
// 			Ra_ret[i][j] = Ra_symMat[i][j] * Ra_diagMat[i];
// 			Ra_ret[j][i] = Ra_symMat[i][j] * Ra_diagMat[j];
// 		}
// 		Ra_ret[i][i] = Ra_symMat[i][i] * Ra_diagMat[i];
// 	}
// #else
// 	for (i=0 ; i<N_s1 ; i++) {
// 		for (j=0 ; j<i ; j++) {
// 			Ra_ret[i][j] = Ra_symMat[i][j] * Ra_diagMat[j];
// 		}
// 		Ra_ret[i][i] = Ra_symMat[i][i] * Ra_diagMat[i];
// 		for (j=i+1 ; j<N_s1 ; j++)
// 			Ra_ret[i][j] = Ra_symMat[j][i] * Ra_diagMat[j];
// 	}
// #endif
// 
// 	return Ra_ret;
// }
// 
// wsReal** sseSD(DIAG_t Ra_diagMat, wsUint N_s1, SYM_t Ra_symMat)
// {
// 	wsUint i, j;
// 	wsReal **Ra_ret = sseMatrix(N_s1, N_s1);
// 
// 	for (i=0 ; i<N_s1 ; i++) {
// 		wsUint N_med = 0;
// #ifdef USE_SSE
// 		N_med = getMed(i);
// 		sse_t sse_di = sseSet(Ra_diagMat[i]);
// 		for (j=0 ; j<N_med ; j+=sseJmp) {
// 			sse_t *sse_sij = (sse_t *)(Ra_symMat[i] + j);
// 			sse_t *sse_rij = (sse_t *)(Ra_ret[i] + j);
// 
// 			*sse_rij = sseMul(*sse_sij, sse_di);
// 			//			Ra_ret[i][j] = Ra_symMat[i][j] * Ra_diagMat[j];
// 			sse_t sse_rji = sseMul(*sse_sij, *(sse_t *)(Ra_diagMat + j));
// 			Ra_ret[j][i] = sse_rij->m128d_f64[0];
// 			Ra_ret[j+1][i] = sse_rij->m128d_f64[1];
// 			//			Ra_ret[j][i] = Ra_symMat[i][j] * Ra_diagMat[i];
// 		}
// #endif
// 		for (j=N_med ; j<i ; j++) {
// 			Ra_ret[i][j] = Ra_symMat[i][j] * Ra_diagMat[i];
// 			Ra_ret[j][i] = Ra_symMat[i][j] * Ra_diagMat[j];
// 		}
// 		Ra_ret[i][i] = Ra_symMat[i][i] * Ra_diagMat[i];
// 	}
// 
// 	return Ra_ret;
// }

wsReal distEuclidean(wsRealCst *Ra_1, wsRealCst *Ra_2, wsUint N_sz)
{
	wsUint N_med = 0, i;
	wsReal R_ret = W0;
#ifdef USE_SSE
	sse_t sse_sum = sseSet(0.0);
	for (i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_1 = (sse_t *)(Ra_1 + i);
		sse_t *sse_2 = (sse_t *)(Ra_2 + i);
		sse_t sse_int = sseSub(*sse_1, *sse_2);;
		sse_sum = sseAdd(sse_sum, sseMul(sse_int, sse_int));
	}
	sseSum(sse_sum, R_ret);
#endif
	for (i=N_med ; i<N_sz ; i++) {
		wsReal R_v = Ra_1[i]-Ra_2[i];
		R_ret += SQR(R_v);
	}

	return R_ret;
}

wsSym sseSymSbsMatrix(wsReal **Rp_inpMatrix, wsUint N_sz,
	wsUint *Na_newSeq, wsUint N_newSz/*=0xffffffff*/)
{
	if (N_newSz == 0xffffffff)
		N_newSz = N_sz;
	wsSym Ra_ret = sseSymMat(N_newSz);

	for (wsUint i=0 ; i<N_newSz ; i++)
		for (wsUint j=0 ; j<=i ; j++)
			Ra_ret[i][j] = Rp_inpMatrix[Na_newSeq[i]][Na_newSeq[j]];

	return Ra_ret;
}

wsReal* sseRowMatrix(wsMat Rp_mat, wsUint N_sz, wsUint N_idx,
	wsUint N_rem/*=0xffffffff*/)
{
	if (N_rem != 0xffffffff) {
		wsReal *Ra_ret = sseVector(N_sz-1);
		memcpy(Ra_ret, Rp_mat[N_idx], N_rem*sizeof(wsReal));
		memcpy(Ra_ret+N_rem, Rp_mat[N_idx]+N_rem+1, sizeof(wsReal)*(N_sz-N_rem-1));
		return Ra_ret;
	}
	wsReal *Ra_ret = sseVector(N_sz);
	memcpy(Ra_ret, Rp_mat[N_idx], N_sz*sizeof(wsReal));
	return Ra_ret;
}

wsReal normFrobenius(wsReal **Ra_mat, wsUint N_r, wsUint N_c)
{
	wsReal R_norm = W0;

	for (wsUint i=0 ; i<N_r ; i++)
		for (wsUint j=0 ; j<N_c ; j++)
			R_norm += SQR(Ra_mat[i][j]);

	return R_norm;
}

wsReal normFrobenius(wsSym Ra_mat, wsUint N_sz)
{
	wsReal R_norm = W0;

	for (wsUint i=0 ; i<N_sz ; i++) {
		for (wsUint j=0 ; j<i ; j++)
			R_norm += W2*SQR(Ra_mat[i][j]);
		R_norm += SQR(Ra_mat[i][i]);
	}

	return R_norm;
}

wsReal compMMfrobenius(wsReal **Ra_m1, wsUint N_r1, wsUint N_c1,
	wsReal **Ra_m2, wsUint N_r2, wsUint N_c2)
{
	wsReal R_norm = W0;

	/* Sanity check */
	if (N_r1 != N_r2 || N_c1 != N_c2)
		halt("Dimension mismatch [%d*%d] and [%d*%d]", N_r1, N_c1, N_r2,
			N_c2);

	for (wsUint i=0 ; i<N_r1 ; i++)
		for (wsUint j=0 ; j<N_c1 ; j++) {
			wsReal R_diff = Ra_m1[i][j] - Ra_m2[i][j];
			R_norm += SQR(R_diff);
		}

	return R_norm;
}

wsReal compSMfrobenius(wsReal **Ra_s1, wsUint N_sz, wsReal **Ra_m2,
	wsUint N_r2, wsUint N_c2)
{
	wsReal R_norm = W0;

	/* Sanity check */
	if (N_sz != N_r2 || N_sz != N_c2)
		halt("Dimension mismatch [%d*%d] and [%d*%d]", N_sz, N_sz, N_r2,
			N_c2);

	for (wsUint i=0 ; i<N_sz ; i++) {
		for (wsUint j=0 ; j<=i ; j++) {
			wsReal R_diff = Ra_s1[i][j] - Ra_m2[i][j];
			R_norm += SQR(R_diff);
		}
		for (wsUint j=i+1 ; j<N_sz ; j++) {
			wsReal R_diff = Ra_s1[j][i] - Ra_m2[i][j];
			R_norm += SQR(R_diff);
		}
	}

	return R_norm;
}

wsReal compSSfrobenius(wsReal **Ra_s1, wsUint N_s1, wsReal **Ra_m2,
	wsUint N_s2)
{
	wsReal R_norm = W0;

	/* Sanity check */
	if (N_s1 != N_s2)
		halt("Dimension mismatch [%d*%d] and [%d*%d]", N_s1, N_s1,
			N_s2, N_s2);

	for (wsUint i=0 ; i<N_s1 ; i++) {
		for (wsUint j=0 ; j<i ; j++) {
			wsReal R_diff = Ra_s1[i][j] - Ra_m2[i][j];
			R_norm += W2 * SQR(R_diff);
		}
		wsReal R_diff = Ra_s1[i][i] - Ra_m2[i][i];
		R_norm += SQR(R_diff);
	}

	return R_norm;
}






int imin2(int x, int y)
{
	return (x < y) ? x : y;
}

double fmin2(double x, double y)
{
	return (x < y) ? x : y;
}

double fmax2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? y : x;
}



#define M_LN_SQRT_2PI 0.918938533204672741780329736406 

double lgammacor(double x);

/* liu c
SKAT_liu <- function(q, lambda, h = rep(1,length(lambda)), delta = rep(0,length(lambda))) {
*/
double liuEV(double q, double *Ra_lambda, wsUint r)
{
	/* h <- 1 */
	/* delta <- 0 */
	double c1 = W0; // c1 <- sum(lambda*h) + sum(lambda*delta)
	double c2 = W0; // c2 <- sum(lambda^2*h) + 2*sum(lambda^2*delta)
	double c3 = W0; // c3 <- sum(lambda^3*h) + 3*sum(lambda^3*delta)
	double c4 = W0; // c4 <- sum(lambda^4*h) + 4*sum(lambda^4*delta)

	for (wsUint i=0 ; i<r ; i++) {
		double x = Ra_lambda[i];
		c1 += x;
		x *= Ra_lambda[i];
		c2 += x;
		x *= Ra_lambda[i];
		c3 += x;
		x *= Ra_lambda[i];
		c4 += x;
	}

	// s1 <- c3/(c2^(3/2))
	double s1 = c3 / pow(c2, 1.5);
	// s2 <- c4/c2^2
	double s2 = c4/SQR(c2);
	// muQ <- c1
	double muQ = c1;
	// sigmaQ <- sqrt(2*c2)
	double sigQ = sqrt(2.0*c2);
	// tstar
	double tstar = (q-muQ)/sigQ;

	double a, delt, l;
	if ((s1*s1) > s2) {
		// a <- 1/(s1-sqrt(s1^2-s2))
		a = 1.0 / (s1 - sqrt((s1*s1)-s2));
		// delta <- s1*a^3-a^2
		delt = s1*(a*a*a) - a*a;
		// l <- a^2-2*delta
		l = a*a - 2.0*delt;
	} else {
		// a <- 1/s1
		a = 1.0/s1;
		// delta <- 0
		delt = 0.0;
		// l <- c2^3/c3^2
		l = CUBE(c2) / SQR(c3);
	}

	// muX <- l+delta
	double muX = l+delt;

	// sigmaX <- sqrt(2)*a
	double sigX = sqrt(2.0)*a;

	return 1.0-PVchisq(tstar*sigX+muX, l, &delt);
	// Qq <- pchisq(tstar*sigmaX+muX,df=l,ncp=delta,lower.tail=FALSE)
	// return(Qq)
}

double wsUnifrand()
{
	return cWisardRandom::getInstance()->getDouble();
// #if RAND_MAX == INT_MAX
// 	double r = wsRand()/(double)INT_MAX;
// #else
// 	double r = (double)wsRand()/((double)RAND_MAX*(double)RAND_MAX+1.0);
// #endif
// //	LOG("[%d (%d,%d)]\n", wsRand(), RAND_MAX, INT_MAX);
// 	return r;
	//return (rand())/(double)((wsUint)RAND_MAX+1);
}

#ifndef M_LN2
#	define M_LN2 0.693147180559945309417
#endif
#define expmax	(DBL_MAX_EXP * M_LN2)/* = log(DBL_MAX) */

double genBeta(double aa, double bb)
{
	double a, b, alpha;
	double r, s, t, u1, u2, v, w, y, z;

	int qsame;
	/* FIXME:  Keep Globals (properly) for threading */
	/* Uses these GLOBALS to save time when many rv's are generated : */
	static double beta, gamma, delta, k1, k2;
	static double olda = -1.0;
	static double oldb = -1.0;

	if (aa <= 0. || bb <= 0. || (!_finite(aa) && !_finite(bb)))
		return WISARD_NAN;

	if (!_finite(aa))
		return 1.0;

	if (!_finite(bb))
		return 0.0;

	/* Test if we need new "initializing" */
	qsame = (olda == aa) && (oldb == bb);
	if (!qsame) { olda = aa; oldb = bb; }

	a = fmin2(aa, bb);
	b = fmax2(aa, bb); /* a <= b */
	alpha = a + b;

#define v_w_from__u1_bet(AA) 			\
	v = beta * log(u1 / (1.0 - u1));	\
	if (v <= expmax) {			\
	w = AA * exp(v);		\
	if(!_finite(w)) w = DBL_MAX;	\
	} else				\
	w = DBL_MAX


	if (a <= 1.0) {	/* --- Algorithm BC --- */

		/* changed notation, now also a <= b (was reversed) */

		if (!qsame) { /* initialize */
			beta = 1.0 / a;
			delta = 1.0 + b - a;
			k1 = delta * (0.0138889 + 0.0416667 * a) / (b * beta - 0.777778);
			k2 = 0.25 + (0.5 + 0.25 / delta) * a;
		}
		/* FIXME: "do { } while()", but not trivially because of "continue"s:*/
		for(;;) {
			u1 = wsUnifrand();
			u2 = wsUnifrand();
			if (u1 < 0.5) {
				y = u1 * u2;
				z = u1 * y;
				if (0.25 * u2 + z - y >= k1)
					continue;
			} else {
				z = u1 * u1 * u2;
				if (z <= 0.25) {
					v_w_from__u1_bet(b);
					break;
				}
				if (z >= k2)
					continue;
			}

			v_w_from__u1_bet(b);

			if (alpha * (log(alpha / (a + w)) + v) - 1.3862944 >= log(z))
				break;
		}
		return (aa == a) ? a / (a + w) : w / (a + w);

	}
	else {		/* Algorithm BB */

		if (!qsame) { /* initialize */
			beta = sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
			gamma = a + 1.0 / beta;
		}
		do {
			u1 = wsUnifrand();
			u2 = wsUnifrand();

			v_w_from__u1_bet(a);

			z = u1 * u1 * u2;
			r = gamma * v - 1.3862944;
			s = a + r - w;
			if (s + 2.609438 >= 5.0 * z)
				break;
			t = log(z);
			if (s > t)
				break;
		}
		while (r + alpha * log(alpha / (b + w)) < t);

		return (aa != a) ? b / (b + w) : w / (b + w);
	}
}

#define repeat for(;;)

double genBinomial(double nin, double pp)
{
    /* FIXME: These should become THREAD_specific globals : */

    static double c, fm, npq, p1, p2, p3, p4, qn;
    static double xl, xll, xlr, xm, xr;

    static double psave = -1.0;
    static int nsave = -1;
    static int m;

    double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
    double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
    int i, ix, k, n;

    if (!_finite(nin)) return WISARD_NAN;
    r = floor(nin + 0.5);
    if (r != nin) return WISARD_NAN;
    if (!_finite(pp) ||
	/* n=0, p=0, p=1 are not errors <TSL>*/
	r < 0 || pp < 0. || pp > 1.)
	return WISARD_NAN;

    if (r == 0 || pp == 0.) return 0;
    if (pp == 1.) return r;

    if (r >= INT_MAX)
		halt("Too large nin");

    /* else */
    n = (int) r;

    p = fmin2(pp, 1. - pp);
    q = 1. - p;
    np = n * p;
    r = p / q;
    g = r * (n + 1);

    /* Setup, perform only when parameters change [using static (globals): */

    /* FIXING: Want this thread safe
       -- use as little (thread globals) as possible
    */
    if (pp != psave || n != nsave) {
	psave = pp;
	nsave = n;
	if (np < 30.0) {
	    /* inverse cdf logic for mean less than 30 */
	    qn = pow(q, (double) n);
	    goto L_np_small;
	} else {
	    ffm = np + p;
	    m = (int) ffm;
	    fm = m;
	    npq = np * q;
	    p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
	    xm = fm + 0.5;
	    xl = xm - p1;
	    xr = xm + p1;
	    c = 0.134 + 20.5 / (15.3 + fm);
	    al = (ffm - xl) / (ffm - xl * p);
	    xll = al * (1.0 + 0.5 * al);
	    al = (xr - ffm) / (xr * q);
	    xlr = al * (1.0 + 0.5 * al);
	    p2 = p1 * (1.0 + c + c);
	    p3 = p2 + c / xll;
	    p4 = p3 + c / xlr;
	}
    } else if (n == nsave) {
	if (np < 30.0)
	    goto L_np_small;
    }

    /*-------------------------- np = n*p >= 30 : ------------------- */
    repeat {
      u = wsUnifrand() * p4;
      v = wsUnifrand();
      /* triangular region */
      if (u <= p1) {
	  ix = (int)(xm - p1 * v + u);
	  goto finis;
      }
      /* parallelogram region */
      if (u <= p2) {
	  x = xl + (u - p1) / c;
	  v = v * c + 1.0 - fabs(xm - x) / p1;
	  if (v > 1.0 || v <= 0.)
	      continue;
	  ix = (int) x;
      } else {
	  if (u > p3) {	/* right tail */
	      ix = (int)(xr - log(v) / xlr);
	      if (ix > n)
		  continue;
	      v = v * (u - p3) * xlr;
	  } else {/* left tail */
	      ix = (int)(xl + log(v) / xll);
	      if (ix < 0)
		  continue;
	      v = v * (u - p2) * xll;
	  }
      }
      /* determine appropriate way to perform accept/reject test */
      k = abs(ix - m);
      if (k <= 20 || k >= npq / 2 - 1) {
	  /* explicit evaluation */
	  f = 1.0;
	  if (m < ix) {
	      for (i = m + 1; i <= ix; i++)
		  f *= (g / i - r);
	  } else if (m != ix) {
	      for (i = ix + 1; i <= m; i++)
		  f /= (g / i - r);
	  }
	  if (v <= f)
	      goto finis;
      } else {
	  /* squeezing using upper and lower bounds on log(f(x)) */
	  amaxp = (k / npq) * ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
	  ynorm = -k * k / (2.0 * npq);
	  alv = log(v);
	  if (alv < ynorm - amaxp)
	      goto finis;
	  if (alv <= ynorm + amaxp) {
	      /* stirling's formula to machine accuracy */
	      /* for the final acceptance/rejection test */
	      x1 = ix + 1;
	      f1 = fm + 1.0;
	      z = n + 1 - fm;
	      w = n - ix + 1.0;
	      z2 = z * z;
	      x2 = x1 * x1;
	      f2 = f1 * f1;
	      w2 = w * w;
	      if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
		  goto finis;
	  }
      }
  }

 L_np_small:
    /*---------------------- np = n*p < 30 : ------------------------- */

  repeat {
     ix = 0;
     f = qn;
     u = wsUnifrand();
     repeat {
	 if (u < f)
	     goto finis;
	 if (ix > 110)
	     break;
	 u -= f;
	 ix++;
	 f *= (g / ix - r);
     }
  }
 finis:
    if (psave > 0.5)
	 ix = n - ix;
  return (double)ix;
}

#if defined(__FreeBSD__) || defined(__APPLE__)
void getcmdline(pid_t pid, char *S)
{
	int mib[4];
	size_t len;

#if defined(__FreeBSD__)
	mib[0] = CTL_KERN;
	mib[1] = KERN_PROC;
	mib[2] = KERN_PROC_ARGS;
	mib[3] = pid;

	len = MAX_PATH;
	if (sysctl(mib, 4, S, &len, NULL, 0) == -1)
		return;
#else
	mib[0] = CTL_KERN;
	mib[1] = KERN_PROCARGS2;
	mib[2] = pid;

	len = MAX_PATH;
	if (sysctl(mib, 3, S, &len, NULL, 0) == -1)
		return;
#endif

	return;
}
#endif

char getItselfPath(char *Sp_path)
{
#ifdef _WIN32
	Sp_path[0] = '\0';
	DWORD r = GetModuleFileName(NULL, Sp_path, MAX_PATH);
	/* Remove file itself */
	for (char *q=Sp_path+r-1 ; q>=Sp_path && *q != '\\' ; q--)
		*q = '\0';
	if (!Sp_path[0]) return 0;
#else
	char	S_tmp[PATH_MAX];
	pid_t	pid = getpid();
#	if defined(__FreeBSD__) || defined(__APPLE__)
	getcmdline(pid, S_tmp);
	//strcpy(Sp_path, S_tmp);
	::realpath(S_tmp, Sp_path);
	for (char *q=Sp_path+strlen(Sp_path)-1 ; q>=Sp_path && *q != '/' ; q--)
		*q = '\0';
#	elif defined(__sun)
	strcpy(Sp_path, getexecname());
#	else
	sprintf(S_tmp, "/proc/%d/exe", pid);
	if (readlink(S_tmp, Sp_path, PATH_MAX) == -1) {
		perror("readlink");
		return 0;
	}
	/* Remove file itself */
	for (char *q=Sp_path+strlen(Sp_path)-1 ; q>=Sp_path && *q != '/' ; q--)
		*q = '\0';
#	endif
#endif
	return 1;
}

#ifdef WIN32
// gettimeofday in windows
int gettimeofday(struct timeval *tv, struct timezone *tz)
{
	FILETIME ft;
	unsigned __int64 tmpres = 0;
	static int tzflag;
	if (NULL != tv) {
		// system time�� ���ϱ�
		GetSystemTimeAsFileTime(&ft);

		// unsigned 64 bit�� �����
		tmpres |= ft.dwHighDateTime;
		tmpres <<= 32;
		tmpres |= ft.dwLowDateTime;
		// 100nano�� 1micro�� ��ȯ�ϱ�
		tmpres /= 10;
		// epoch time���� ��ȯ�ϱ�
		tmpres -= DELTA_EPOCH_IN_MICROSECS;    

		// sec�� micorsec���� ���߱�
		tv->tv_sec = (long)(tmpres / 1000000UL);
		tv->tv_usec = (tmpres % 1000000UL);
	}
	// timezone ó��
	if (NULL != tz) {
		if (!tzflag) {
			_tzset();
			tzflag++;
		}
		tz->tz_minuteswest = _timezone / 60;
		tz->tz_dsttime = _daylight;
	}
	return 0;
}
#endif

void multcomp(cIO *Cp_IO, vector<double>& chi, string title)
{
	wsUint N_mkr = Cp_IO->sizeVariant();
	if (N_mkr == 0) return;
	
	vInt tcnt;
	if (N_mkr != chi.size())
		halt("SYSERR: Internal problem in multiple comparison routine");
	LOG("Computing corrected significance values (fdr, Sidak, etc)\n");
	bool altern_pval = tcnt.size() > 0;

	wsUint		N_avail	= 0;
	xRealSort*	Xa_p	= NULL;
	xRealSort*	Xa_chi	= NULL;
	wsAlloc(Xa_p, xRealSort, N_mkr);
	wsAlloc(Xa_chi, xRealSort, N_mkr);
	for (wsUint i=0 ; i<N_mkr ; i++) {
		if (chi[i] < 0 || chi[i] != chi[i]) continue;
		double p = altern_pval ? pT(sqrt(chi[i]), tcnt[i]) :
			PVchisq(chi[i], 1);
		if (p <= -1) continue;

		/* Insert */
		Xa_p[N_avail].V	= p;
		Xa_p[N_avail].i	= i;
		Xa_chi[N_avail].V	= chi[i];
		Xa_chi[N_avail].i	= altern_pval ? (int)tcnt[i] : 0;

		/* Increase index */
		N_avail++;
	}

	if (N_avail == 0) {
		LOG("Zero valid tests computed -- no adjusted values calculated\n");
		return;
	}

	/* Sort p-values */
	qsort(Xa_p, N_avail, sizeof(xRealSort), sort_real_desc);
	qsort(Xa_chi, N_avail, sizeof(xRealSort), sort_real_desc);
	wsReal R_avail = (wsReal)N_avail;

	/* Genomic control */
	double R_lambda;
	double lambda_mean = 0;
	if (N_avail % 2 == 0 )
		R_lambda = ( Xa_chi[N_avail/ 2 - 1 ].V + Xa_chi[ N_avail / 2 ].V ) / 2 ;
	else
		R_lambda = Xa_chi[ (N_avail-1) / 2 ].V ;
	for (wsUint i=0; i<N_avail ; i++)
		lambda_mean += Xa_chi[i].V;
	lambda_mean /= R_avail;

	if (IS_ASSIGNED(usergc)) {
		R_lambda = OPT_REAL(usergc);
	} else {
		R_lambda /= 0.456;
		if (R_lambda < 1) R_lambda = 1.00;
	}

	/* Print */
	LOG("Genomic control 'lambda' is [%g] (%s)\n", R_lambda,
		IS_ASSIGNED(usergc)?"Given as parameter":
		"Computed from median of chi-squared statistics");

	LOG("Mean chi-squared statistic is %g\n", R_lambda);
	LOG("Correcting for %d tests\n", N_avail);

	// Consider each test
	// Bonferroni correction

	vector<double> pv_GC(N_avail);
	vector<double> pv_sidakSS(N_avail);
	vector<double> pv_sidakSD(N_avail);
	vector<double> pv_holm(N_avail);
	vector<double> pv_BH(N_avail);
	vector<double> pv_BY(N_avail);

	// Genomic control (reverse order)
	int i2=0;
	for (int i=(int)N_avail-1 ; i>=0 ; i--)
	{
		pv_GC[i2++] = altern_pval ? 
			pT(sqrt(Xa_chi[i].V / R_lambda), Xa_chi[i].i) : 
			PVchisq(Xa_chi[i].V / R_lambda, 1 );      
	}

	// Base adjust values on GC p-values?
	if (OPT_ENABLED(genoctrl)) {
		LOG("Using genomic-controlled p-values for adjusted p-values\n"); 
		for (wsUint i=0 ; i<N_avail ; i++)
			Xa_p[i].V = pv_GC[i];
	}

	// Holm 
	pv_holm[0] = Xa_p[0].V*R_avail > 1 ? 1 : Xa_p[0].V*R_avail;
	for (wsUint i=1 ; i<N_avail ; i++)
	{
		double x = (N_avail-i)*Xa_p[i].V < 1 ? (N_avail-i)*Xa_p[i].V : 1;
		pv_holm[i] = pv_holm[i-1] > x ? pv_holm[i-1] : x;
	}

	// Sidak SS
	for (wsUint i=0 ; i<N_avail ; i++)
		pv_sidakSS[i] = 1 - pow(W1 - Xa_p[i].V, R_avail);


	// Sidak SD
	pv_sidakSD[0] = 1 - pow(W1 - Xa_p[0].V , R_avail);
	for (wsUint i=1 ; i<N_avail ; i++) {
		double x = W1 - pow(W1 - Xa_p[i].V , R_avail-i);
		pv_sidakSD[i] = pv_sidakSD[i-1] > x ? pv_sidakSD[i-1] : x ; 
	}

	// BH
	pv_BH[N_avail-1] = Xa_p[N_avail-1].V;
	for (int i=(int)N_avail-2 ; i>=0 ; i--) {
		double u = (R_avail/(double)(i+1))*Xa_p[i].V;
		double x = u < 1 ? u : W1;
		pv_BH[i] = pv_BH[i+1] < x ? pv_BH[i+1] : x;
	}

	// BY
	double a = 0;
	for (double i=1 ; i<=R_avail ; i++)
		a += 1/i;

	pv_BY[N_avail-1] = a * Xa_p[N_avail-1].V < 1 ? a * Xa_p[N_avail-1].V : 1 ; 

	for (int i=N_avail-2 ; i>=0 ; i--) {      
		double u = ((R_avail*a)/(double)(i+1))*Xa_p[i].V;
		double x = u < 1 ? u : W1;
		pv_BY[i] = pv_BY[i+1] < x ? pv_BY[i+1] : x;
	}

	// Output
	
	cExporter *Cp_ex = cExporter::summon(title.c_str());
 	LOGoutput("Multiple-test corrected p-values are exported to [%s.%s]\n",
		OPT_STRING(out), title.c_str());
	headerVariant(Cp_ex);
	Cp_ex->put("	P_ORIG	GC	P_BON	P_HOLM	P_SIDAKSS P_SIDAKSD	P_FDRBH	P_FDRBY\n");
// 	MT << setw(4) << "CHR" << " "
// 		<< setw(par::pp_maxsnp) << "SNP" << " "
// 		<< setw(10) << "UNADJ" << " "
// 		<< setw(10) << "GC" << " ";
// 	if ( par::qq_plot ) 
// 		MT << setw(10) << "QQ" << " ";
// 	MT << setw(10) << "BONF" << " "
// 		<< setw(10) << "HOLM" << " "
// 		<< setw(10) << "SIDAK_SS" << " "
// 		<< setw(10) << "SIDAK_SD" << " "
// 		<< setw(10) << "fdr_BH" << " "
// 		<< setw(10) << "fdr_BY" << "\n";
// 
	vVariant& Xa_mkr = Cp_IO->getVariant();
	for (wsUint l=0; l<N_avail; l++) {
		entryVariant(Cp_ex, Xa_mkr[Xa_p[l].i]);

		double bonferroni = Xa_p[l].V*R_avail > 1 ? 1 : Xa_p[l].V*R_avail;
		Cp_ex->fmt("	%g	%g	%g	%g	%g	%g	%g	%g\n", Xa_p[l].V,
			pv_GC[l], bonferroni, pv_holm[l], pv_sidakSS[l], pv_sidakSD[l],
			pv_BH[l], pv_BY[l]);
	}
	delete Cp_ex;
// 		MT << setw(4) << locus[sp[l].l]->chr << " "
// 			<< setw(par::pp_maxsnp) << locus[sp[l].l]->name << " ";
// 
// 		// Unadjusted
// 		pprint(MT,sp[l].p);
// 
// 		// Genomic control
// 		pprint(MT,pv_GC[l]);
// 
// 		// Q-Q plot?
// 		if ( par::qq_plot ) 
// 		{
// 			pprint(MT,(l+0.5)/(double)t);
// 		}
// 
// 		// Bonferroni, etc
// 		double bonferroni = sp[l].p*t > 1 ? 1 : sp[l].p*t;
// 		pprint(MT,bonferroni);
// 		pprint(MT,pv_holm[l]);
// 		pprint(MT,pv_sidakSS[l]);
// 		pprint(MT,pv_sidakSD[l]);
// 		pprint(MT,pv_BH[l]);
// 		pprint(MT,pv_BY[l]);
// 		MT << "\n";
// 	}
// 
// 	MT.close();  
}

void multcomp(cIO *Cp_IO, double *Ra_stat, wsUint N_stat, string title)
{
	wsUint N_mkr = Cp_IO->sizeVariant();
	if (N_mkr == 0) return;

	vInt tcnt;
	if (N_mkr != N_stat)
		halt("SYSERR: Internal problem in multiple comparison routine");
	LOG("Computing corrected significance values (fdr, Sidak, etc)\n");
	bool altern_pval = tcnt.size() > 0;

	wsUint		N_avail	= 0;
	xRealSort*	Xa_p	= NULL;
	xRealSort*	Xa_chi	= NULL;
	wsAlloc(Xa_p, xRealSort, N_mkr);
	wsAlloc(Xa_chi, xRealSort, N_mkr);
	for (wsUint i=0 ; i<N_mkr ; i++) {
		if (Ra_stat[i] < 0 || NA(Ra_stat[i])) continue;
		double p = altern_pval ? pT(sqrt(Ra_stat[i]), tcnt[i]) :
			PVchisq(Ra_stat[i], 1);
		if (p <= -1) continue;

		/* Insert */
		Xa_p[N_avail].V		= p;
		Xa_p[N_avail].i	= i;
		Xa_chi[N_avail].V	= Ra_stat[i];
		Xa_chi[N_avail].i	= altern_pval ? (int)tcnt[i] : 0;

		/* Increase index */
		N_avail++;
	}

	if (N_avail == 0) {
		LOG("Zero valid tests computed -- no adjusted values calculated\n");
		return;
	}

	/* Sort p-values */
	qsort(Xa_p, N_avail, sizeof(xRealSort), sort_real_desc);
	qsort(Xa_chi, N_avail, sizeof(xRealSort), sort_real_desc);
	wsReal R_avail = (wsReal)N_avail;

	/* Genomic control */
	double R_lambda;
	double lambda_mean = 0;
	if (N_avail % 2 == 0 )
		R_lambda = ( Xa_chi[N_avail/ 2 - 1 ].V + Xa_chi[ N_avail / 2 ].V ) / 2 ;
	else
		R_lambda = Xa_chi[ (N_avail-1) / 2 ].V ;
	for (wsUint i=0; i<N_avail ; i++)
		lambda_mean += Xa_chi[i].V;
	lambda_mean /= R_avail;

	if (IS_ASSIGNED(usergc)) {
		R_lambda = OPT_REAL(usergc);
	} else {
		R_lambda /= 0.456;
		if (R_lambda < 1) R_lambda = 1.00;
	}

	/* Print */
	LOG("Genomic control 'lambda' is [%g] (%s)\n", R_lambda,
		IS_ASSIGNED(usergc)?"Given as parameter":
		"Computed from median of chi-squared statistics");

	LOG("Mean chi-squared statistic is %g\n", R_lambda);
	LOG("Correcting for %d tests\n", N_avail);

	// Consider each test
	// Bonferroni correction

	vector<double> pv_GC(N_avail);
	vector<double> pv_sidakSS(N_avail);
	vector<double> pv_sidakSD(N_avail);
	vector<double> pv_holm(N_avail);
	vector<double> pv_BH(N_avail);
	vector<double> pv_BY(N_avail);

	// Genomic control (reverse order)
	int i2=0;
	for (int i=(int)N_avail-1 ; i>=0 ; i--)
	{
		pv_GC[i2++] = altern_pval ? 
			pT(sqrt(Xa_chi[i].V / R_lambda), Xa_chi[i].i) : 
			PVchisq(Xa_chi[i].V / R_lambda, 1 );      
	}

	// Base adjust values on GC p-values?
	if (OPT_ENABLED(genoctrl)) {
		LOG("Using genomic-controlled p-values for adjusted p-values\n"); 
		for (wsUint i=0 ; i<N_avail ; i++)
			Xa_p[i].V = pv_GC[i];
	}

	// Holm 
	pv_holm[0] = Xa_p[0].V*R_avail > 1 ? 1 : Xa_p[0].V*R_avail;
	for (wsUint i=1 ; i<N_avail ; i++)
	{
		double x = (N_avail-i)*Xa_p[i].V < 1 ? (N_avail-i)*Xa_p[i].V : 1;
		pv_holm[i] = pv_holm[i-1] > x ? pv_holm[i-1] : x;
	}

	// Sidak SS
	for (wsUint i=0 ; i<N_avail ; i++)
		pv_sidakSS[i] = 1 - pow(W1 - Xa_p[i].V, R_avail);


	// Sidak SD
	pv_sidakSD[0] = 1 - pow(W1 - Xa_p[0].V , R_avail);
	for (wsUint i=1 ; i<N_avail ; i++) {
		double x = W1 - pow(W1 - Xa_p[i].V , R_avail-i);
		pv_sidakSD[i] = pv_sidakSD[i-1] > x ? pv_sidakSD[i-1] : x ; 
	}

	// BH
	pv_BH[N_avail-1] = Xa_p[N_avail-1].V;
	for (int i=(int)N_avail-2 ; i>=0 ; i--) {
		double u = (R_avail/(double)(i+1))*Xa_p[i].V;
		double x = u < 1 ? u : W1;
		pv_BH[i] = pv_BH[i+1] < x ? pv_BH[i+1] : x;
	}

	// BY
	double a = 0;
	for (double i=1 ; i<=R_avail ; i++)
		a += 1/i;

	pv_BY[N_avail-1] = a * Xa_p[N_avail-1].V < 1 ? a * Xa_p[N_avail-1].V : 1 ; 

	for (int i=N_avail-2 ; i>=0 ; i--) {      
		double u = ((R_avail*a)/(double)(i+1))*Xa_p[i].V;
		double x = u < 1 ? u : W1;
		pv_BY[i] = pv_BY[i+1] < x ? pv_BY[i+1] : x;
	}

	// Output

	cExporter *Cp_ex = cExporter::summon(title.c_str());
	LOGoutput("Multiple-test corrected p-values are exported to [%s.%s]\n",
		OPT_STRING(out), title.c_str());
	headerVariant(Cp_ex);
	Cp_ex->put("	P_ORIG	GC	P_BON	P_HOLM	P_SIDAKSS P_SIDAKSD	P_FDRBH	P_FDRBY\n");
	// 	MT << setw(4) << "CHR" << " "
	// 		<< setw(par::pp_maxsnp) << "SNP" << " "
	// 		<< setw(10) << "UNADJ" << " "
	// 		<< setw(10) << "GC" << " ";
	// 	if ( par::qq_plot ) 
	// 		MT << setw(10) << "QQ" << " ";
	// 	MT << setw(10) << "BONF" << " "
	// 		<< setw(10) << "HOLM" << " "
	// 		<< setw(10) << "SIDAK_SS" << " "
	// 		<< setw(10) << "SIDAK_SD" << " "
	// 		<< setw(10) << "fdr_BH" << " "
	// 		<< setw(10) << "fdr_BY" << "\n";
	// 
	vVariant& Xa_mkr = Cp_IO->getVariant();
	for (wsUint l=0; l<N_avail; l++) {
		entryVariant(Cp_ex, Xa_mkr[Xa_p[l].i]);

		double bonferroni = Xa_p[l].V*R_avail > 1 ? 1 : Xa_p[l].V*R_avail;
		Cp_ex->fmt("	%g	%g	%g	%g	%g	%g	%g	%g\n", Xa_p[l].V,
			pv_GC[l], bonferroni, pv_holm[l], pv_sidakSS[l], pv_sidakSD[l],
			pv_BH[l], pv_BY[l]);
	}
	delete Cp_ex;
	// 		MT << setw(4) << locus[sp[l].l]->chr << " "
	// 			<< setw(par::pp_maxsnp) << locus[sp[l].l]->name << " ";
	// 
	// 		// Unadjusted
	// 		pprint(MT,sp[l].p);
	// 
	// 		// Genomic control
	// 		pprint(MT,pv_GC[l]);
	// 
	// 		// Q-Q plot?
	// 		if ( par::qq_plot ) 
	// 		{
	// 			pprint(MT,(l+0.5)/(double)t);
	// 		}
	// 
	// 		// Bonferroni, etc
	// 		double bonferroni = sp[l].p*t > 1 ? 1 : sp[l].p*t;
	// 		pprint(MT,bonferroni);
	// 		pprint(MT,pv_holm[l]);
	// 		pprint(MT,pv_sidakSS[l]);
	// 		pprint(MT,pv_sidakSD[l]);
	// 		pprint(MT,pv_BH[l]);
	// 		pprint(MT,pv_BY[l]);
	// 		MT << "\n";
	// 	}
	// 
	// 	MT.close();  
}

char fetchMarkerList(mDataIdx* &Xm_mmap, wsStrCst S_fn, vInt &Xv_ret, wsStrCst S_desc/*=NULL*/)
{
	/* Build marker map */
// 	wsUint		j = 0;
// 	FOREACHDO (vMarker_it, Xv_snp, i, j++) {
// 		if (i->chr > 0 && i->chr <= MAX_CHR)
// 			Xm_mmap[i->chr].insert(make_pair((string)i->name, j));
// 		else
// 			Xm_mmap[0].insert(make_pair((string)i->name, j));
// 	}
	char	N_ret	= 0;
	char*	Sp_buf	= NULL;
	char*	a		= NULL;
	wsAlloc(Sp_buf, char, 8192);

	/* Try to open file */
	cStrFile C_in(S_fn, S_desc, 1);
	if (C_in.isFailed()) goto _end;

	/* Read the first line */
	if (C_in.gets(Sp_buf, 8192) == NULL)
		goto _end;

	/* Check the first element is 'CHR' */
	getString(&Sp_buf, &a);
	if (!stricmp(Sp_buf, "CHR")) while (C_in.gets(Sp_buf, 8192)) {
		char *b = NULL;
		a = NULL;
		getString(&Sp_buf, &a);
		/* Read second column */
		getString(&a, &b);
		/* a = MARKER */
		string	S_key	= (string)a;
		wsUint	N_idx	= 0;
		wsUint	j		= 0;
		for (j=0 ; j<=NCHR_SPECIES ; j++) {
			mDataIdx_it X_find = Xm_mmap[j].find(S_key);
			if (X_find != Xm_mmap[j].end()) {
				N_idx = X_find->second;
				break;
			}
		}
		/* Failed to found */
		if (j > NCHR_SPECIES) continue;
		Xv_ret.push_back(N_idx);
	} else do {
		if (!a) getString(&Sp_buf, &a);
		/* Sp_buf = MARKER */
		string	S_key	= (string)Sp_buf;
		wsUint	N_idx	= 0;
		wsUint	j		= 0;
		for (j=0 ; j<=NCHR_SPECIES ; j++) {
			mDataIdx_it X_find = Xm_mmap[j].find(S_key);
			if (X_find != Xm_mmap[j].end()) {
				N_idx = X_find->second;
				break;
			}
		}
		a = NULL;
		/* Failed to found */
		if (j > NCHR_SPECIES) continue;
		Xv_ret.push_back(N_idx);
	} while (C_in.gets(Sp_buf, 8192));

	/* Mark as success */
	N_ret = 1;
_end:
	DEALLOC(Sp_buf);
	return N_ret;
}

char fetchMarkerList(wsStrCst S_fn, vStr &Xv_ret, wsStrCst S_desc/*=NULL*/)
{
	/* Build marker map */
	// 	wsUint		j = 0;
	// 	FOREACHDO (vMarker_it, Xv_snp, i, j++) {
	// 		if (i->chr > 0 && i->chr <= MAX_CHR)
	// 			Xm_mmap[i->chr].insert(make_pair((string)i->name, j));
	// 		else
	// 			Xm_mmap[0].insert(make_pair((string)i->name, j));
	// 	}
	char	N_ret	= 0;
	char*	Sp_buf	= NULL;
	char*	a		= NULL;
	wsAlloc(Sp_buf, char, 8192);

	/* Try to open file */
	cStrFile C_in(S_fn, S_desc, 1);
	if (C_in.isFailed()) goto _end;

	/* Read the first line */
	if (C_in.gets(Sp_buf, 8192) == NULL)
		goto _end;

	/* Check the first element is 'CHR' */
	getString(&Sp_buf, &a);
	if (!stricmp(Sp_buf, "CHR")) while (C_in.gets(Sp_buf, 8192)) {
		char *b = NULL;
		a = NULL;
		getString(&Sp_buf, &a);
		/* Read second column */
		getString(&a, &b);
		Xv_ret.push_back(a);
	} else do {
		if (!a) getString(&Sp_buf, &a);
		Xv_ret.push_back(Sp_buf);
	} while (C_in.gets(Sp_buf, 8192));

	/* Mark as success */
	N_ret = 1;
_end:
	DEALLOC(Sp_buf);
	return N_ret;
}

char fetchMarkerList(wsStrCst S_fn, eStr &Xe_ret, wsStrCst S_desc/*=NULL*/)
{
	char	N_ret	= 0;
	char*	Sp_buf	= NULL;
	char*	a		= NULL;
	char*	b		= NULL;
	wsAlloc(Sp_buf, char, 8192);

	/* Try to open file */
	cStrFile C_in(S_fn, S_desc, 1);
	if (C_in.isFailed()) goto _end2;

	/* Read the first line */
	if (C_in.gets(Sp_buf, 8192) == NULL)
		goto _end;

	/* Check the first element is 'CHR' */
	getString(&Sp_buf, &a);
	if (!stricmp(Sp_buf, "CHR")) while (C_in.gets(Sp_buf, 8192)) {
		b = Sp_buf;
		a = NULL;
		getString(&b, &a);
		/* Read second column */
		getString(&a, &b);
		Xe_ret.insert(a);
	} else do {
		b = Sp_buf;
		if (!a) getString(&b, &a);
		if (b[0])
			Xe_ret.insert(b);
	} while (C_in.gets(Sp_buf, 8192));

	/* Mark as success */
	N_ret = 1;
_end:
	DEALLOC(Sp_buf);
	return N_ret;
_end2:
	LOG("Marker list [%s] looks like not a file, try to get as list\n", S_fn);
	/* Regarding it is not a file */
	wsUint N_fil = 0;
	char ** Sa_fils = loadStringValues(S_fn, &N_fil);
	for (wsUint i=0 ; i<N_fil ; i++)
		Xe_ret.insert(Sa_fils[i]);
	DEALLOC(Sp_buf);
	N_ret = 1;
	return N_ret;
}

bool fileExists(wsStrCst S_fn)
{
	FILE *fp = fopen(S_fn, "r");
	if (fp == NULL) return false;
	fclose(fp);
	return true;
}

void getEigenPhi(wsSym Ra_phi, wsUint N_anaSamp, wsMat *Rp_eVec,
	wsReal **Rp_eVal, wsUint *Na_testFamIdx/*=NULL*/)
{
	char	S_fn[512];
	char	B_ret	= 0;
	wsMat	Ra_eVec	= NULL;
	wsReal*	Ra_eVal	= NULL;
	wsStrCst	S_fx = IS_ASSIGNED(ev) ? OPT_STRING(ev) : (OPT_ENABLED(makeev) ?
		OPT_STRING(out) : NULL);
	
	/* Pass if null */
	if (S_fx == NULL) goto _calc;

	sprintf(S_fn, "%s.eigen.val", S_fx);
	/* Do calc if not exists */
	if (!fileExists(S_fn)) goto _calc;

	/* Find file */
#ifdef USE_ED2
	sprintf(S_fn, "%s.eigen.vect", S_fx);
#else
	sprintf(S_fn, "%s.eigen.vec", S_fx);
#endif

	/* Do calc if not exists */
	if (!fileExists(S_fn)) goto _calc;

	/* Try to init from file */ {
		wsUint N_r = 0, N_c = 0;
		Ra_eVec = makeMatrixSSE(S_fn, &N_r, &N_c);
		if (N_r != N_c || N_r != N_anaSamp) {
			sseUnmat(Ra_eVec, N_r);
			goto _calc;
	// 		halt("Sample size[%d * %d] is not match with analysis size[%d]",
	// 			N_r, N_c, N_anaSamp);
		}

		sprintf(S_fn, "%s.eigen.val", S_fx);
		wsMat Ra_eval = makeMatrixSSE(S_fn, &N_r, &N_c);
		if (N_r != 1 || N_c != N_anaSamp) {
			sseUnmat(Ra_eVec, N_anaSamp);
			sseUnmat(Ra_eval, N_r);
			goto _calc;
// 			halt("Sample size[%d * %d] is not match with analysis size[%d]",
// 				N_r, N_c, N_anaSamp);
		}
		Ra_eVal = Ra_eval[0];
		DEALLOC(Ra_eval);
		LOG("Cached eigendecomposition [%s.eigen.vec] and [%s.eigen.val] loaded\n",
			S_fx, S_fx);
	}
#ifndef USE_ED2
	/* Need to be transposed in this part
	 * P -> P^t */
	transposeSelf(Ra_eVec, N_anaSamp, N_anaSamp);
#endif
	goto _ret;

_calc:
	if (Na_testFamIdx)
		Ra_eVal = famEIGENDECOMPOSITION(Ra_phi, N_anaSamp,
			Na_testFamIdx, &Ra_eVec);
	else
		Ra_eVal = EIGENDECOMPOSITION(Ra_phi, N_anaSamp, &Ra_eVec);
	if (B_ret == 1)
		halt("Eigendecomposition failed");
	if (OPT_ENABLED(makeev)) {
		exportMatrix("eigen.vect", Ra_eVec, N_anaSamp, N_anaSamp);
		exportVector("eigen.val", Ra_eVal, N_anaSamp);
	}
_ret:
	*Rp_eVec = Ra_eVec;
	*Rp_eVal = Ra_eVal;
}

/* Table of constant values */

typedef double (*D_fp)(...), (*E_fp)(...);
static int c__8000 = 8000;
static int c__0 = 0;
static int c__2 = 2;
static int c__4 = 4;
static double c_b19 = 1.;

double mvbvtl_(int *nu, double *dh, double *dk, double *
	r__);
double mvphi_(double *z__);
double mvbvu_(double *sh, double *sk, double *r__);
double mvbvt_(int *nu, double *_lower, double *_upper, int *
	infin, double *correl);
double stdjac_(int *nu, double *t);
double studnt_(int *nu, double *t);
/* Subroutine */ int mvtsrt_(int *n, int *nu, double *_lower, 
	double *_upper, double *correl, int *infin, double *y, 
	int *infis, double *a, double *b, int *infi, 
	double *cov, double *d__, double *e);
/* Subroutine */ int mvtlms_(int *nu, double *a, double *b, 
	int *infin, double *_lower, double *_upper);
double stdinv_(int *n, double *z__);
double mvtnit_(int *n, int *nuin, double *correl, double *
	_lower, double *_upper, int *infin, int *infis, double *
	d__, double *e);
double fncmvt_(int *n, double *w);
double fulsum_(int *s, double *center, double *hwidth, 
	double *x, double *g, D_fp f);
/* Subroutine */ int rcswap_(int *p, int *q, double *a, 
	double *b, int *infin, int *n, double *c__);
/* Subroutine */ int limits_(double *a, double *b, int *infin, 
	double *_lower, double *_upper);
double bvnu_(double *sh, double *sk, double *r__);
double phi_(double *z__);
double bvn_(double *_lower, double *_upper, int *infin, 
	double *correl);
/* Subroutine */ int ncvsrt_(int *n, double *_lower, double *_upper,
	double *correl, int *infin, double *y, int *infis, 
	double *a, double *b, int *infi, double *cov, 
	double *d__, double *e);
double phinv_(double *p);
/* Subroutine */ int rulnrm_(int *lenrul, int *numnul, int *
	rulpts, double *w, double *rulcon);
/* Subroutine */ int trestr_(int *pointr, int *sbrgns, double *
	pontrs, double *rgners);
/* Subroutine */ int basrul_(int *ndim, double *a, double *b, 
	double *width, D_fp functn, double *w, int *lenrul, 
	double *g, double *center, double *z__, double *
	rgnert, double *basest);
/* Subroutine */ int differ_(int *ndim, double *a, double *b, 
	double *width, double *z__, double *dif, D_fp functn, 
	int *divaxn, int *difcls);
/* Subroutine */ int bsinit_(int *ndim, double *w, int *lenrul, 
	double *g);
double mvnnit_(int *n, double *correl, double *_lower, 
	double *_upper, int *infin, int *infis, double *d__, 
	double *e);
double mvnfnc_(int *n, double *w);
int adapt_(int *ndim, int *mincls, int *maxcls, 
	D_fp functn, double *absreq, double *relreq, int *lenwrk, 
	double *work, double *absest, double *finest, int *
	inform__);
/* Subroutine */ int adbase_(int *ndim, int *mincls, int *maxcls, 
	D_fp functn, double *absreq, double *relreq, double *
	absest, double *finest, int *sbrgns, int *mxrgns, int 
	*rulcls, int *lenrul, double *errors, double *values, 
	double *pontrs, double *lowers, double *uppers, 
	double *meshes, double *weghts, double *points, 
	double *_lower, double *_upper, double *width, double *
	mesh, double *work, int *inform__);

inline double d_sign(double *a, double *b)
{
	double x;
	x = (*a >= 0 ? *a : - *a);
	return( *b >= 0 ? x : -x);
}

int pow_ii(int *ap, int *bp)
{
	int pow, x, n;
	unsigned long u;

	x = *ap;
	n = *bp;

	if (n <= 0) {
		if (n == 0 || x == 1)
			return 1;
		if (x != -1)
			return x == 0 ? 1/x : 0;
		n = -n;
	}
	u = n;
	for(pow = 1; ; )
	{
		if(u & 01)
			pow *= x;
		if(u >>= 1)
			x *= x;
		else
			break;
	}
	return(pow);
}

/* Selected portions of code taken from: */
/*    http://www.math.wsu.edu/math/faculty/genz/software/mvn.f */
/*    http://www.math.wsu.edu/math/faculty/genz/software/mvt.f */
/* with a few minor modifications (search for 'AA' to find them) */

/* Author: */
/*          Alan Genz */
/*          Department of Mathematics */
/*          Washington State University */
/*          Pullman, WA 99164-3113 */
/*          Email : alangenz@wsu.edu */
/* except for some auxiliary functions whose authors are indicated */
/* in the respective code below. */
/* ----------------------------------------------------------------------- */
/* Subroutine */ int sadmvn_(int *n, double *_lower, double *_upper,
	 int *infin, double *correl, int *maxpts, double *
	abseps, double *releps, double *error, double *value, 
	int *inform__)
{
    /* System generated locals */
    int i__1, i__2, i__3, i__4;
    double d__1, d__2;

    /* Local variables */
	static double d__ = 0, e = 0;
    static int m = 0;
    static int infis = 0;
    static double oldval = 0;
    static int maxcls = 0, newcls = 0, rulcls = 0, totcls = 0;
	double *work = NULL;
	wsCalloc(work, double, *maxpts);
	//memset(work, 0x00, sizeof(double)*8000);


/*     A subroutine for computing multivariate normal probabilities. */
/*     This subroutine uses an algorithm given in the paper */
/*     "Numerical Computation of Multivariate Normal Probabilities", in */
/*     J. of Computational and Graphical Stat., 1(1992), pp. 141-149, by */
/*          Alan Genz */
/*          Department of Mathematics */
/*          Washington State University */
/*          Pullman, WA 99164-3113 */
/*          Email : alangenz@wsu.edu */

/*  Parameters */

/*     N      INTEGER, the number of variables. */
/*     LOWER  REAL, array of lower integration limits. */
/*     UPPER  REAL, array of upper integration limits. */
/*     INFIN  INTEGER, array of integration limits flags: */
/*            if INFIN(I) < 0, Ith limits are (-infinity, infinity); */
/*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
/*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
/*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
/*     CORREL REAL, array of correlation coefficients; the correlation */
/*            coefficient in row I column J of the correlation matrix */
/*            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I. */
/*     MAXPTS INTEGER, maximum number of function values allowed. This */
/*            parameter can be used to limit the time taken. A */
/*            sensible strategy is to start with MAXPTS = 1000*N, and then */
/*            increase MAXPTS if ERROR is too large. */
/*     ABSEPS REAL absolute error tolerance. */
/*     RELEPS REAL relative error tolerance. */
/*     ERROR  REAL estimated absolute error, with 99% confidence level. */
/*     VALUE  REAL estimated value for the integral */
/*     INFORM INTEGER, termination status parameter: */
/*            if INFORM = 0, normal completion with ERROR < EPS; */
/*            if INFORM = 1, completion with ERROR > EPS and MAXPTS */
/*                           function vaules used; increase MAXPTS to */
/*                           decrease ERROR; */
/*            if INFORM = 2, N > 20 or N < 1. */

    /* Parameter adjustments */
    --correl;
    --infin;
    --_upper;
    --_lower;

    /* Function Body */
    if (*n > 100 || *n < 1) {
	*inform__ = 2;
	*value = 0.;
	*error = 1.;
	return 0;
    }
    *inform__ = (int) mvnnit_(n, &correl[1], &_lower[1], &_upper[1], &infin[
	    1], &infis, &d__, &e);
    m = *n - infis;
    if (m == 0) {
	*value = 1.;
	*error = 0.;
    } else if (m == 1) {
	*value = e - d__;
	*error = (float)2e-16;
    } else {

/*        Call the subregion adaptive integration subroutine */

	--m;
	rulcls = 1;
	adapt_(&m, &rulcls, &c__0, (D_fp)mvnfnc_, abseps, releps, &c__8000, 
		work, error, value, inform__);
/* Computing MIN */
	i__1 = rulcls * 10;
	maxcls = min(i__1,*maxpts);
	totcls = 0;
	adapt_(&m, &totcls, &maxcls, (D_fp)mvnfnc_, abseps, releps, &c__8000, 
		work, error, value, inform__);
/* Computing MAX */
	d__1 = *abseps, d__2 = *releps * fabs(*value);
	if (*error > max(d__1,d__2)) {
L10:
	    oldval = *value;
/* Computing MAX */
/* Computing MIN */
	    i__3 = maxcls * 3 / 2, i__4 = *maxpts - totcls;
	    i__1 = rulcls << 1, i__2 = min(i__3,i__4);
	    maxcls = max(i__1,i__2);
	    newcls = -1;
	    adapt_(&m, &newcls, &maxcls, (D_fp)mvnfnc_, abseps, releps, &
		    c__8000, work, error, value, inform__);
	    totcls += newcls;
/* Computing 2nd power */
	    d__2 = *error;
	    *error = (d__1 = *value - oldval, fabs(d__1)) + sqrt(rulcls * (
		    d__2 * d__2) / totcls);
/* Computing MAX */
	    d__1 = *abseps, d__2 = *releps * fabs(*value);
	    if (*error > max(d__1,d__2)) {
		if (*maxpts - totcls > rulcls << 1) {
		    goto L10;
		}
	    } else {
		*inform__ = 0;
	    }
	}
    }
    return 0;
} /* sadmvn_ */

/* -------------------------------------------------------------------------- */

/* Subroutine */ int adapt_(int *ndim, int *mincls, int *maxcls, 
	D_fp functn, double *absreq, double *relreq, int *lenwrk, 
	double *work, double *absest, double *finest, int *
	inform__)
{
    /* Local variables */
    static int inmesh, invals, sbrgns, inwdth, lenrul, inerrs, rulcls, 
	    inmshs, inwork, inlowr, inpnts, inwgts, inuppr, mxrgns, inptrs, 
	    inlwrs, inuprs;


/*   Adaptive Multidimensional Integration Subroutine */

/*   Author: Alan Genz */
/*           Department of Mathematics */
/*           Washington State University */
/*           Pullman, WA 99164-3113 USA */

/*  This subroutine computes an approximation to the integral */

/*      1 1     1 */
/*     I I ... I       FUNCTN(NDIM,X)  dx(NDIM)...dx(2)dx(1) */
/*      0 0     0 */

/* **************  Parameters for ADAPT  ******************************** */

/* ***** Input Parameters */

/*  NDIM    Integer number of integration variables. */
/*  MINCLS  Integer minimum number of FUNCTN calls to be allowed; MINCLS */
/*          must not exceed MAXCLS. If MINCLS < 0, then ADAPT assumes */
/*          that a previous call of ADAPT has been made with the same */
/*          integrand and continues that calculation. */
/*  MAXCLS  Integer maximum number of FUNCTN calls to be used; MAXCLS */
/*          must be >= RULCLS, the number of function calls required for */
/*          one application of the basic integration rule. */
/*           IF ( NDIM .EQ. 1 ) THEN */
/*              RULCLS = 11 */
/*           ELSE IF ( NDIM .LT. 15 ) THEN */
/*              RULCLS = 2**NDIM + 2*NDIM*(NDIM+3) + 1 */
/*           ELSE */
/*              RULCLS = 1 + NDIM*(24-NDIM*(6-NDIM*4))/3 */
/*           ENDIF */
/*  FUNCTN  Externally declared real user defined integrand. Its */
/*          parameters must be (NDIM, Z), where Z is a real array of */
/*          length NDIM. */
/*  ABSREQ  Real required absolute accuracy. */
/*  RELREQ  Real required relative accuracy. */
/*  LENWRK  Integer length of real array WORK (working storage); ADAPT */
/*          needs LENWRK >= 16*NDIM + 27. For maximum efficiency LENWRK */
/*          should be about 2*NDIM*MAXCLS/RULCLS if MAXCLS FUNCTN */
/*          calls are needed. If LENWRK is significantly less than this, */
/*          ADAPT may be less efficient. */

/* ***** Output Parameters */

/*  MINCLS  Actual number of FUNCTN calls used by ADAPT. */
/*  WORK    Real array (length LENWRK) of working storage. This contains */
/*          information that is needed for additional calls of ADAPT */
/*          using the same integrand (input MINCLS < 0). */
/*  ABSEST  Real estimated absolute accuracy. */
/*  FINEST  Real estimated value of integral. */
/*  INFORM  INFORM = 0 for normal exit, when ABSEST <= ABSREQ or */
/*                     ABSEST <= |FINEST|*RELREQ with MINCLS <= MAXCLS. */
/*          INFORM = 1 if MAXCLS was too small for ADAPT to obtain the */
/*                     result FINEST to within the requested accuracy. */
/*          INFORM = 2 if MINCLS > MAXCLS, LENWRK < 16*NDIM + 27 or */
/*                     RULCLS > MAXCLS. */

/* *********************************************************************** */

/*     Begin driver routine. This routine partitions the working storage */
/*      array and then calls the main subroutine ADBASE. */

    /* Parameter adjustments */
    --work;

    /* Function Body */
    if (*ndim == 1) {
	lenrul = 5;
	rulcls = 9;
    } else if (*ndim < 12) {
	lenrul = 6;
	rulcls = pow_ii(&c__2, ndim) + (*ndim << 1) * (*ndim + 2) + 1;
    } else {
	lenrul = 6;
	rulcls = (*ndim << 1) * ((*ndim << 1) + 1) + 1;
    }
    if (*lenwrk >= lenrul * (*ndim + 4) + *ndim * 10 + 3 && rulcls <= *maxcls 
	    && *mincls <= *maxcls) {
	mxrgns = (*lenwrk - lenrul * (*ndim + 4) - *ndim * 7) / (*ndim * 3 + 
		3);
	inerrs = 1;
	invals = inerrs + mxrgns;
	inptrs = invals + mxrgns;
	inlwrs = inptrs + mxrgns;
	inuprs = inlwrs + mxrgns * *ndim;
	inmshs = inuprs + mxrgns * *ndim;
	inwgts = inmshs + mxrgns * *ndim;
	inpnts = inwgts + (lenrul << 2);
	inlowr = inpnts + lenrul * *ndim;
	inuppr = inlowr + *ndim;
	inwdth = inuppr + *ndim;
	inmesh = inwdth + *ndim;
	inwork = inmesh + *ndim;
	if (*mincls < 0) {
	    sbrgns = (int) work[*lenwrk];
	}
	adbase_(ndim, mincls, maxcls, (D_fp)functn, absreq, relreq, absest, 
		finest, &sbrgns, &mxrgns, &rulcls, &lenrul, &work[inerrs], &
		work[invals], &work[inptrs], &work[inlwrs], &work[inuprs], &
		work[inmshs], &work[inwgts], &work[inpnts], &work[inlowr], &
		work[inuppr], &work[inwdth], &work[inmesh], &work[inwork], 
		inform__);
	work[*lenwrk] = (double) sbrgns;
    } else {
	*inform__ = 2;
	*mincls = rulcls;
    }
    return 0;
} /* adapt_ */


/* Subroutine */ int adbase_(int *ndim, int *mincls, int *maxcls, 
	D_fp functn, double *absreq, double *relreq, double *
	absest, double *finest, int *sbrgns, int *mxrgns, int 
	*rulcls, int *lenrul, double *errors, double *values, 
	double *pontrs, double *lowers, double *uppers, 
	double *meshes, double *weghts, double *points, 
	double *_lower, double *_upper, double *width, double *
	mesh, double *work, int *inform__)
{
    /* System generated locals */
    int lowers_dim1, lowers_offset, uppers_dim1, uppers_offset, 
	    meshes_dim1, meshes_offset, i__1, i__2;
    double d__1, d__2;

    /* Local variables */
    static int i__ = 0, j = 0;
    static int difcls = 0, divaxn = 0, rgncls = 0, funcls = 0;
    static int nwrgns = 0;
    static int top = 0;


/*        Main adaptive integration subroutine */


/*     Initialization of subroutine */

    /* Parameter adjustments */
    meshes_dim1 = *ndim;
    meshes_offset = 1 + meshes_dim1 * 1;
    meshes -= meshes_offset;
    uppers_dim1 = *ndim;
    uppers_offset = 1 + uppers_dim1 * 1;
    uppers -= uppers_offset;
    lowers_dim1 = *ndim;
    lowers_offset = 1 + lowers_dim1 * 1;
    lowers -= lowers_offset;
    --errors;
    --values;
    --pontrs;
    --weghts;
    --points;
    --_lower;
    --_upper;
    --width;
    --mesh;
    --work;

    /* Function Body */
    *inform__ = 2;
    funcls = 0;
    bsinit_(ndim, &weghts[1], lenrul, &points[1]);
    if (*mincls >= 0) {

/*       When MINCLS >= 0 determine initial subdivision of the */
/*       integration region and apply basic rule to each subregion. */

	*sbrgns = 0;
	i__1 = *ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    _lower[i__] = 0.;
	    mesh[i__] = 1.;
	    width[i__] = 1 / (mesh[i__] * 2);
	    _upper[i__] = 1.;
	}
	divaxn = 0;
	rgncls = *rulcls;
	nwrgns = 1;
L10:
	differ_(ndim, &_lower[1], &_upper[1], &width[1], &work[1], &work[*ndim 
		+ 1], (D_fp)functn, &divaxn, &difcls);
	funcls += difcls;
	if (funcls + rgncls * (mesh[divaxn] + 1) / mesh[divaxn] <= (
		double) (*mincls)) {
	    rgncls = (int) (rgncls * (mesh[divaxn] + 1) / mesh[divaxn]);
	    nwrgns = (int) (nwrgns * (mesh[divaxn] + 1) / mesh[divaxn]);
	    ++mesh[divaxn];
	    width[divaxn] = 1 / (mesh[divaxn] * 2);
	    goto L10;
	}
	if (nwrgns <= *mxrgns) {
	    i__1 = *ndim;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		_upper[i__] = _lower[i__] + width[i__] * 2;
		mesh[i__] = 1.;
	    }
	}

/*     Apply basic rule to subregions and store results in heap. */

L20:
	++(*sbrgns);
	basrul_(ndim, &_lower[1], &_upper[1], &width[1], (D_fp)functn, &weghts[
		1], lenrul, &points[1], &work[1], &work[*ndim + 1], &errors[*
		sbrgns], &values[*sbrgns]);
	trestr_(sbrgns, sbrgns, &pontrs[1], &errors[1]);
	i__1 = *ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    lowers[i__ + *sbrgns * lowers_dim1] = _lower[i__];
	    uppers[i__ + *sbrgns * uppers_dim1] = _upper[i__];
	    meshes[i__ + *sbrgns * meshes_dim1] = mesh[i__];
	}
	i__1 = *ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    _lower[i__] = _upper[i__];
	    _upper[i__] = _lower[i__] + width[i__] * 2;
	    if (_lower[i__] + width[i__] < 1.) {
		goto L20;
	    }
	    _lower[i__] = 0.;
	    _upper[i__] = _lower[i__] + width[i__] * 2;
	}
	funcls += *sbrgns * *rulcls;
    }

/*     Check for termination */

L30:
    *finest = 0.;
    *absest = 0.;
    i__1 = *sbrgns;
    for (i__ = 1; i__ <= i__1; ++i__) {
	*finest += values[i__];
	*absest += errors[i__];
    }
/* Computing MAX */
    d__1 = *absreq, d__2 = *relreq * fabs(*finest);
    if (*absest > max(d__1,d__2) || funcls < *mincls) {

/*     Prepare to apply basic rule in (parts of) subregion with */
/*     largest error. */

	top = (int) pontrs[1];
	rgncls = *rulcls;
	i__1 = *ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    _lower[i__] = lowers[i__ + top * lowers_dim1];
	    _upper[i__] = uppers[i__ + top * uppers_dim1];
	    mesh[i__] = meshes[i__ + top * meshes_dim1];
	    width[i__] = (_upper[i__] - _lower[i__]) / (mesh[i__] * 2);
	    rgncls = (int) (rgncls * mesh[i__]);
	}
	differ_(ndim, &_lower[1], &_upper[1], &width[1], &work[1], &work[*ndim 
		+ 1], (D_fp)functn, &divaxn, &difcls);
	funcls += difcls;
	rgncls = (int) (rgncls * (mesh[divaxn] + 1) / mesh[divaxn]);
	if (funcls + rgncls <= *maxcls) {
	    if (*sbrgns + 1 <= *mxrgns) {

/*     Prepare to subdivide into two pieces. */

		nwrgns = 1;
		width[divaxn] /= 2;
	    } else {
		nwrgns = 0;
		width[divaxn] = width[divaxn] * mesh[divaxn] / (mesh[divaxn] 
			+ 1);
		meshes[divaxn + top * meshes_dim1] = mesh[divaxn] + 1;
	    }
	    if (nwrgns > 0) {

/*     Only allow local subdivision when space is available. */

		i__1 = *sbrgns + nwrgns;
		for (j = *sbrgns + 1; j <= i__1; ++j) {
		    i__2 = *ndim;
		    for (i__ = 1; i__ <= i__2; ++i__) {
			lowers[i__ + j * lowers_dim1] = _lower[i__];
			uppers[i__ + j * uppers_dim1] = _upper[i__];
			meshes[i__ + j * meshes_dim1] = mesh[i__];
		    }
		}
		uppers[divaxn + top * uppers_dim1] = _lower[divaxn] + width[
			divaxn] * 2;
		lowers[divaxn + (*sbrgns + 1) * lowers_dim1] = uppers[divaxn 
			+ top * uppers_dim1];
	    }
	    funcls += rgncls;
	    basrul_(ndim, &lowers[top * lowers_dim1 + 1], &uppers[top * 
		    uppers_dim1 + 1], &width[1], (D_fp)functn, &weghts[1], 
		    lenrul, &points[1], &work[1], &work[*ndim + 1], &errors[
		    top], &values[top]);
	    trestr_(&top, sbrgns, &pontrs[1], &errors[1]);
	    i__1 = *sbrgns + nwrgns;
	    for (i__ = *sbrgns + 1; i__ <= i__1; ++i__) {

/*     Apply basic rule and store results in heap. */

		basrul_(ndim, &lowers[i__ * lowers_dim1 + 1], &uppers[i__ * 
			uppers_dim1 + 1], &width[1], (D_fp)functn, &weghts[1],
			 lenrul, &points[1], &work[1], &work[*ndim + 1], &
			errors[i__], &values[i__]);
		trestr_(&i__, &i__, &pontrs[1], &errors[1]);
	    }
	    *sbrgns += nwrgns;
	    goto L30;
	} else {
	    *inform__ = 1;
	}
    } else {
	*inform__ = 0;
    }
    *mincls = funcls;
    return 0;
} /* adbase_ */

/* Subroutine */ int bsinit_(int *ndim, double *w, int *lenrul, 
	double *g)
{
    /* System generated locals */
    int w_dim1, w_offset, g_dim1, g_offset, i__1, i__2;
    double d__1;

    /* Builtin functions */
    int pow_ii(int *, int *);

    /* Local variables */
    static double lamp = 0;
    static int i__ = 0, j = 0;
    static double rulcon = 0;
	static int rulpts[6] = { 0, };
    static double lam1 = 0, lam2 = 0, lam3 = 0;


/*     For initializing basic rule weights and symmetric sum parameters. */


/*     The following code determines rule parameters and weights for a */
/*      degree 7 rule (W(1,1),...,W(5,1)), two degree 5 comparison rules */
/*      (W(1,2),...,W(5,2) and W(1,3),...,W(5,3)) and a degree 3 */
/*      comparison rule (W(1,4),...W(5,4)). */

/*       If NDIM = 1, then LENRUL = 5 and total points = 9. */
/*       If NDIM < SDIM, then LENRUL = 6 and */
/*                      total points = 1+2*NDIM*(NDIM+2)+2**NDIM. */
/*       If NDIM > = SDIM, then LENRUL = 6 and */
/*                      total points = 1+2*NDIM*(1+2*NDIM). */

    /* Parameter adjustments */
    g_dim1 = *ndim;
    g_offset = 1 + g_dim1 * 1;
    g -= g_offset;
    w_dim1 = *lenrul;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;

    /* Function Body */
    i__1 = *lenrul;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *ndim;
	for (j = 1; j <= i__2; ++j) {
	    g[j + i__ * g_dim1] = 0.;
	}
	for (j = 1; j <= 4; ++j) {
	    w[i__ + j * w_dim1] = 0.;
	}
    }
    rulpts[4] = (*ndim << 1) * (*ndim - 1);
    rulpts[3] = *ndim << 1;
    rulpts[2] = *ndim << 1;
    rulpts[1] = *ndim << 1;
    rulpts[0] = 1;
    lamp = (float).85;
    lam3 = (float).4707;
    lam2 = 4 / (15 - 5 / lam3);
/* Computing 2nd power */
    d__1 = lam2;
    w[w_dim1 + 5] = (3 - lam3 * 5) / ((lam2 - lam3) * 180 * (d__1 * d__1));
    if (*ndim < 12) {
	lam1 = lam3 * 8 * (lam3 * 31 - 15) / ((lam3 * 3 - 1) * (lam3 * 5 - 3) 
		* 35);
/* Computing 3rd power */
	d__1 = lam3 * 3;
	w[*lenrul + w_dim1] = 1 / (d__1 * (d__1 * d__1)) / pow_ii(&c__2, ndim)
		;
    } else {
	lam1 = (lam3 * (15 - lam2 * 21) + (*ndim - 1) * 35 * (lam2 - lam3) / 
		9) / (lam3 * (21 - lam2 * 35) + (*ndim - 1) * 35 * (lam2 / 
		lam3 - 1) / 9);
/* Computing 3rd power */
	d__1 = lam3 * 3;
	w[w_dim1 + 6] = 1 / (d__1 * (d__1 * d__1) * 4);
    }
    w[w_dim1 + 3] = (15 - (lam3 + lam1) * 21 + lam3 * 35 * lam1) / (lam2 * 
	    210 * (lam2 - lam3) * (lam2 - lam1)) - ((*ndim - 1) << 1) * w[
	    w_dim1 + 5];
    w[w_dim1 + 2] = (15 - (lam3 + lam2) * 21 + lam3 * 35 * lam2) / (lam1 * 
	    210 * (lam1 - lam3) * (lam1 - lam2));
    if (*ndim < 12) {
	rulpts[*lenrul - 1] = pow_ii(&c__2, ndim);
	lam3 = sqrt(lam3);
	i__1 = *ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    g[i__ + *lenrul * g_dim1] = lam3;
	}
    } else {
/* Computing 3rd power */
	d__1 = lam3 * 3;
	w[w_dim1 + 6] = 1 / (d__1 * (d__1 * d__1) * 4);
	rulpts[5] = (*ndim << 1) * (*ndim - 1);
	lam3 = sqrt(lam3);
	for (i__ = 1; i__ <= 2; ++i__) {
	    g[i__ + g_dim1 * 6] = lam3;
	}
    }
    if (*ndim > 1) {
/* Computing 2nd power */
	d__1 = lam2 * 6;
	w[(w_dim1 << 1) + 5] = 1 / (d__1 * d__1);
/* Computing 2nd power */
	d__1 = lam2 * 6;
	w[w_dim1 * 3 + 5] = 1 / (d__1 * d__1);
    }
    w[(w_dim1 << 1) + 3] = (3 - lam1 * 5) / (lam2 * 30 * (lam2 - lam1)) -
		((*ndim - 1) << 1) * w[(w_dim1 << 1) + 5];
    w[(w_dim1 << 1) + 2] = (3 - lam2 * 5) / (lam1 * 30 * (lam1 - lam2));
    w[w_dim1 * 3 + 4] = (3 - lam2 * 5) / (lamp * 30 * (lamp - lam2));
    w[w_dim1 * 3 + 3] = (3 - lamp * 5) / (lam2 * 30 * (lam2 - lamp)) -
		((*ndim - 1) << 1) * w[w_dim1 * 3 + 5];
    w[(w_dim1 << 2) + 2] = 1 / (lam1 * 6);
    lamp = sqrt(lamp);
    lam2 = sqrt(lam2);
    lam1 = sqrt(lam1);
    g[(g_dim1 << 1) + 1] = lam1;
    g[g_dim1 * 3 + 1] = lam2;
    g[(g_dim1 << 2) + 1] = lamp;
    if (*ndim > 1) {
	g[g_dim1 * 5 + 1] = lam2;
	g[g_dim1 * 5 + 2] = lam2;
    }
    for (j = 1; j <= 4; ++j) {
	w[j * w_dim1 + 1] = 1.;
	i__1 = *lenrul;
	for (i__ = 2; i__ <= i__1; ++i__) {
	    w[j * w_dim1 + 1] -= rulpts[i__ - 1] * w[i__ + j * w_dim1];
	}
    }
    rulcon = 2.;
    rulnrm_(lenrul, &c__4, rulpts, &w[w_offset], &rulcon);
    return 0;
} /* bsinit_ */


/* Subroutine */ int rulnrm_(int *lenrul, int *numnul, int *
	rulpts, double *w, double *rulcon)
{
    /* System generated locals */
    int w_dim1, w_offset, i__1, i__2, i__3;

    /* Local variables */
    static int i__ = 0, j = 0, k = 0;
    static double alpha = 0, normcf = 0, normnl = 0;


/*     Compute orthonormalized null rules. */

    /* Parameter adjustments */
    w_dim1 = *lenrul;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;
    --rulpts;

    /* Function Body */
    normcf = 0.;
    i__1 = *lenrul;
    for (i__ = 1; i__ <= i__1; ++i__) {
	normcf += rulpts[i__] * w[i__ + w_dim1] * w[i__ + w_dim1];
    }
    i__1 = *numnul;
    for (k = 2; k <= i__1; ++k) {
	i__2 = *lenrul;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w[i__ + k * w_dim1] -= w[i__ + w_dim1];
	}
	i__2 = k - 1;
	for (j = 2; j <= i__2; ++j) {
	    alpha = 0.;
	    i__3 = *lenrul;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		alpha += rulpts[i__] * w[i__ + j * w_dim1] * w[i__ + k * 
			w_dim1];
	    }
	    alpha = -alpha / normcf;
	    i__3 = *lenrul;
	    for (i__ = 1; i__ <= i__3; ++i__) {
		w[i__ + k * w_dim1] += alpha * w[i__ + j * w_dim1];
	    }
	}
	normnl = 0.;
	i__2 = *lenrul;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    normnl += rulpts[i__] * w[i__ + k * w_dim1] * w[i__ + k * w_dim1];
	}
	alpha = sqrt(normcf / normnl);
	i__2 = *lenrul;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w[i__ + k * w_dim1] = alpha * w[i__ + k * w_dim1];
	}
    }
    i__1 = *numnul;
    for (j = 2; j <= i__1; ++j) {
	i__2 = *lenrul;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    w[i__ + j * w_dim1] /= *rulcon;
	}
    }
    return 0;
} /* rulnrm_ */

/* -------------------------------------------------------------------------- */
double mvnfnc_0_(int n__, int *n, double *w, double *correl, 
	double *_lower, double *_upper, int *infin, int *infis, 
	double *d__, double *e)
{
    /* System generated locals */
    int i__1, i__2;
    double ret_val, d__1, d__2;

    /* Local variables */
	static int infi[100] = { 0, };
	static double prod, a[100] = { 0, }, b[100] = { 0, };
    static int i__, j;
	static double y[100] = { 0, };
    static double d1 = 0, e1 = 0, di, ei;
    static int ij;
	static double cov[5050] = { 0, }, sum;


/*     Integrand subroutine */

    /* Parameter adjustments */
    if (w) {
	--w;
	}
    if (correl) {
	--correl;
	}
    if (_lower) {
	--_lower;
	}
    if (_upper) {
	--_upper;
	}
    if (infin) {
	--infin;
	}

    /* Function Body */
    switch(n__) {
	case 1: goto L_mvnnit;
	}

    di = d1;
    ei = e1;
    prod = e1 - d1;
    ij = 1;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	d__1 = di + w[i__] * (ei - di);
	y[i__ - 1] = phinv_(&d__1);
	sum = 0.;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    ++ij;
	    sum += cov[ij - 1] * y[j - 1];
	}
	++ij;
	if (cov[ij - 1] > 0.) {
	    d__1 = a[i__] - sum;
	    d__2 = b[i__] - sum;
	    limits_(&d__1, &d__2, &infi[i__], &di, &ei);
	} else {
	    d__1 = a[i__] - sum;
	    di = (d_sign(&c_b19, &d__1) + 1) / 2;
	    d__1 = b[i__] - sum;
	    ei = (d_sign(&c_b19, &d__1) + 1) / 2;
	}
	prod *= ei - di;
    }
    ret_val = prod;
    return ret_val;

/*     Entry point for intialization. */


L_mvnnit:
    ret_val = 0.;

/*     Initialization and computation of covariance Cholesky factor. */

    ncvsrt_(n, &_lower[1], &_upper[1], &correl[1], &infin[1], y, infis, a, b, 
	    infi, cov, d__, e);
    d1 = *d__;
    e1 = *e;
    if (*n - *infis == 2) {
/* Computing 2nd power */
	d__1 = cov[1];
	*d__ = sqrt(d__1 * d__1 + 1);
	a[1] /= *d__;
	b[1] /= *d__;
	d__1 = cov[1] / *d__;
	*e = bvn_(a, b, infi, &d__1);
	*d__ = 0.;
	++(*infis);
    }
    return ret_val;
} /* mvnfnc_ */

double mvnfnc_(int *n, double *w)
{
    return mvnfnc_0_(0, n, w, (double *)0, (double *)0, (double *)
	    0, (int *)0, (int *)0, (double *)0, (double *)0);
    }

double mvnnit_(int *n, double *correl, double *_lower, 
	double *_upper, int *infin, int *infis, double *d__, 
	double *e)
{
    return mvnfnc_0_(1, n, (double *)0, correl, _lower, _upper, infin, 
	    infis, d__, e);
    }

/* Subroutine */ int limits_(double *a, double *b, int *infin, 
	double *_lower, double *_upper)
{
    *_lower = 0.;
    *_upper = 1.;
    if (*infin >= 0) {
	if (*infin != 0) {
	    *_lower = phi_(a);
	}
	if (*infin != 1) {
	    *_upper = phi_(b);
	}
    }
    return 0;
} /* limits_ */

/* -------------------------------------------------------------------------- */
double phi_(double *z__)
{
    /* System generated locals */
    double ret_val, d__1;

    /* Local variables */
    static double zabs = 0, p = 0, expntl = 0;


/*     Normal distribution probabilities accurate to 1.e-15. */
/*     Z = no. of standard deviations from the mean. */

/*     Based upon algorithm 5666 for the error function, from: */
/*     Hart, J.F. et al, 'Computer Approximations', Wiley 1968 */

/*     Programmer: Alan Miller */

/*     Latest revision - 30 March 1986 */


    zabs = fabs(*z__);

/*     |Z| > 37 */

    if (zabs > 37.) {
	p = 0.;
    } else {

/*     |Z| <= 37 */

/* Computing 2nd power */
	d__1 = zabs;
	expntl = exp(-(d__1 * d__1) / 2);

/*     |Z| < CUTOFF = 10/SQRT(2) */

	if (zabs < 7.071067811865475) {
	    p = expntl * ((((((zabs * .03526249659989109 + .7003830644436881) 
		    * zabs + 6.37396220353165) * zabs + 33.912866078383) * 
		    zabs + 112.0792914978709) * zabs + 221.2135961699311) * 
		    zabs + 220.2068679123761) / (((((((zabs * 
		    .08838834764831844 + 1.755667163182642) * zabs + 
		    16.06417757920695) * zabs + 86.78073220294608) * zabs + 
		    296.5642487796737) * zabs + 637.3336333788311) * zabs + 
		    793.8265125199484) * zabs + 440.4137358247522);

/*     |Z| >= CUTOFF. */

	} else {
	    p = expntl / (zabs + 1 / (zabs + 2 / (zabs + 3 / (zabs + 4 / (
		    zabs + .65))))) / 2.506628274631001;
	}
    }
    if (*z__ > 0.) {
	p = 1 - p;
    }
    ret_val = p;
    return ret_val;
} /* phi_ */


double phinv_(double *p)
{
    /* System generated locals */
    double ret_val, d__1, d__2;

    /* Local variables */
    static double q, r__;


/* 	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3 */

/* 	Produces the normal deviate Z corresponding to a given lower */
/* 	tail area of P. */

/* 	The hash sums below are the sums of the mantissas of the */
/* 	coefficients.   They are included for use in checking */
/* 	transcription. */


/*     Coefficients for P close to 0.5 */

/*     HASH SUM AB    55.88319 28806 14901 4439 */

/*     Coefficients for P not close to 0, 0.5 or 1. */

/*     HASH SUM CD    49.33206 50330 16102 89036 */

/* 	Coefficients for P near 0 or 1. */

/*     HASH SUM EF    47.52583 31754 92896 71629 */

    q = (*p * 2 - 1) / 2;
    if (fabs(q) <= .425) {
	r__ = .180625 - q * q;
	ret_val = q * (((((((r__ * 2509.0809287301226727 + 
		33430.575583588128105) * r__ + 67265.770927008700853) * r__ + 
		45921.953931549871457) * r__ + 13731.693765509461125) * r__ + 
		1971.5909503065514427) * r__ + 133.14166789178437745) * r__ + 
		3.387132872796366608) / (((((((r__ * 5226.495278852854561 + 
		28729.085735721942674) * r__ + 39307.89580009271061) * r__ + 
		21213.794301586595867) * r__ + 5394.1960214247511077) * r__ + 
		687.1870074920579083) * r__ + 42.313330701600911252) * r__ + 
		1);
    } else {
/* Computing MIN */
	d__1 = *p, d__2 = 1 - *p;
	r__ = min(d__1,d__2);
	if (r__ > 0.) {
	    r__ = sqrt(-log(r__));
	    if (r__ <= 5.) {
		r__ += -1.6;
		ret_val = (((((((r__ * 7.7454501427834140764e-4 + 
			.0227238449892691845833) * r__ + 
			.24178072517745061177) * r__ + 1.27045825245236838258)
			 * r__ + 3.64784832476320460504) * r__ + 
			5.7694972214606914055) * r__ + 4.6303378461565452959) 
			* r__ + 1.42343711074968357734) / (((((((r__ * 
			1.05075007164441684324e-9 + 5.475938084995344946e-4) *
			 r__ + .0151986665636164571966) * r__ + 
			.14810397642748007459) * r__ + .68976733498510000455) 
			* r__ + 1.6763848301838038494) * r__ + 
			2.05319162663775882187) * r__ + 1);
	    } else {
		r__ += -5.;
		ret_val = (((((((r__ * 2.01033439929228813265e-7 + 
			2.71155556874348757815e-5) * r__ + 
			.0012426609473880784386) * r__ + 
			.026532189526576123093) * r__ + .29656057182850489123)
			 * r__ + 1.7848265399172913358) * r__ + 
			5.4637849111641143699) * r__ + 6.6579046435011037772) 
			/ (((((((r__ * 2.04426310338993978564e-15 + 
			1.4215117583164458887e-7) * r__ + 
			1.8463183175100546818e-5) * r__ + 
			7.868691311456132591e-4) * r__ + 
			.0148753612908506148525) * r__ + 
			.13692988092273580531) * r__ + .59983220655588793769) 
			* r__ + 1);
	    }
	} else {
	    ret_val = 9.;
	}
	if (q < 0.) {
	    ret_val = -ret_val;
	}
    }
    return ret_val;
} /* phinv_ */

double bvn_(double *_lower, double *_upper, int *infin, 
	double *correl)
{
    /* System generated locals */
    double ret_val = 0.0, d__1, d__2, d__3, d__4;

/*     A function for computing bivariate normal probabilities. */

/*  Parameters */

/*     LOWER  REAL, array of _lower integration limits. */
/*     UPPER  REAL, array of _upper integration limits. */
/*     INFIN  INTEGER, array of integration limits flags: */
/*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
/*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
/*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
/*     CORREL REAL, correlation coefficient. */

    /* Parameter adjustments */
    --infin;
    --_upper;
    --_lower;

    /* Function Body */
    if (infin[1] == 2 && infin[2] == 2) {
	ret_val = bvnu_(&_lower[1], &_lower[2], correl) - bvnu_(&_upper[1], &
		_lower[2], correl) - bvnu_(&_lower[1], &_upper[2], correl) + 
		bvnu_(&_upper[1], &_upper[2], correl);
    } else if (infin[1] == 2 && infin[2] == 1) {
	ret_val = bvnu_(&_lower[1], &_lower[2], correl) - bvnu_(&_upper[1], &
		_lower[2], correl);
    } else if (infin[1] == 1 && infin[2] == 2) {
	ret_val = bvnu_(&_lower[1], &_lower[2], correl) - bvnu_(&_lower[1], &
		_upper[2], correl);
    } else if (infin[1] == 2 && infin[2] == 0) {
	d__1 = -_upper[1];
	d__2 = -_upper[2];
	d__3 = -_lower[1];
	d__4 = -_upper[2];
	ret_val = bvnu_(&d__1, &d__2, correl) - bvnu_(&d__3, &d__4, correl);
    } else if (infin[1] == 0 && infin[2] == 2) {
	d__1 = -_upper[1];
	d__2 = -_upper[2];
	d__3 = -_upper[1];
	d__4 = -_lower[2];
	ret_val = bvnu_(&d__1, &d__2, correl) - bvnu_(&d__3, &d__4, correl);
    } else if (infin[1] == 1 && infin[2] == 0) {
	d__1 = -_upper[2];
	d__2 = -(*correl);
	ret_val = bvnu_(&_lower[1], &d__1, &d__2);
    } else if (infin[1] == 0 && infin[2] == 1) {
	d__1 = -_upper[1];
	d__2 = -(*correl);
	ret_val = bvnu_(&d__1, &_lower[2], &d__2);
    } else if (infin[1] == 1 && infin[2] == 1) {
	ret_val = bvnu_(&_lower[1], &_lower[2], correl);
    } else if (infin[1] == 0 && infin[2] == 0) {
	d__1 = -_upper[1];
	d__2 = -_upper[2];
	ret_val = bvnu_(&d__1, &d__2, correl);
    }
    return ret_val;
} /* bvn_ */

double bvnu_(double *sh, double *sk, double *r__)
{
    /* Initialized data */

    static struct {
	double e_1[3];
	double fill_2[7];
	double e_3[6];
	double fill_4[4];
	double e_5[10];
	} equiv_89 = {
		{ .1713244923791705, .3607615730481384, .4679139345726904 },
		{ 0, },
		{ .04717533638651177, .1069393259953183, .1600783285433464, .2031674267230659, .2334925365383547, .2491470458134029 },
		{ 0, },
		{ .01761400713915212, 
		.04060142980038694, .06267204833410906, .08327674157670475, 
		.1019301198172404, .1181945319615184, .1316886384491766, 
		.1420961093183821, .1491729864726037, .1527533871307259 } };

#define w ((double *)&equiv_89)

    static struct {
	double e_1[3];
	double fill_2[7];
	double e_3[6];
	double fill_4[4];
	double e_5[10];
	} equiv_90 = {
		{ -.9324695142031522, -.6612093864662647, -.238619186083197 },
		{ 0, },
		{ -.9815606342467191, -.904117256370475, -.769902674194305, -.5873179542866171, -.3678314989981802,  -.1252334085114692 },
		{ 0, },
		{ -.9931285991850949, 
		-.9639719272779138, -.9122344282513259, -.8391169718222188, 
		-.7463319064601508, -.636053680726515, -.5108670019508271, 
		-.3737060887154196, -.2277858511416451, -.07652652113349733 } };

#define x ((double *)&equiv_90)


    /* System generated locals */
    int i__1;
    double ret_val, d__1, d__2, d__3, d__4;

    /* Local variables */
    static double a, b, c__, d__, h__;
    static int i__;
    static double k;
    static int lg;
    static double as;
    static int ng;
    static double bs, hk, hs, sn, rs, xs;
    static double bvn, asr;


/*     A function for computing bivariate normal probabilities. */

/*       Yihong Ge */
/*       Department of Computer Science and Electrical Engineering */
/*       Washington State University */
/*       Pullman, WA 99164-2752 */
/*       Email : yge@eecs.wsu.edu */
/*     and */
/*       Alan Genz */
/*       Department of Mathematics */
/*       Washington State University */
/*       Pullman, WA 99164-3113 */
/*       Email : alangenz@wsu.edu */

/* BVN - calculate the probability that X is larger than SH and Y is */
/*       larger than SK. */

/* Parameters */

/*   SH  REAL, integration limit */
/*   SK  REAL, integration limit */
/*   R   REAL, correlation coefficient */
/*   LG  INTEGER, number of Gauss Rule Points and Weights */

/*     Gauss Legendre Points and Weights, N =  6 */
/*     Gauss Legendre Points and Weights, N = 12 */
/*     Gauss Legendre Points and Weights, N = 20 */
    if (fabs(*r__) < (float).3) {
	ng = 1;
	lg = 3;
    } else if (fabs(*r__) < (float).75) {
	ng = 2;
	lg = 6;
    } else {
	ng = 3;
	lg = 10;
    }
    h__ = *sh;
    k = *sk;
    hk = h__ * k;
    bvn = 0.;
    if (fabs(*r__) < (float).925) {
	hs = (h__ * h__ + k * k) / 2;
	asr = asin(*r__);
	i__1 = lg;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    sn = sin(asr * (x[i__ + ng * 10 - 11] + 1) / 2);
	    bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn))
		    ;
	    sn = sin(asr * (-x[i__ + ng * 10 - 11] + 1) / 2);
	    bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn))
		    ;
/* L10: */
	}
	d__1 = -h__;
	d__2 = -k;
	bvn = bvn * asr / 12.566370614359172 + phi_(&d__1) * phi_(&d__2);
    } else {
	if (*r__ < 0.) {
	    k = -k;
	    hk = -hk;
	}
	if (fabs(*r__) < 1.) {
	    as = (1 - *r__) * (*r__ + 1);
	    a = sqrt(as);
/* Computing 2nd power */
	    d__1 = h__ - k;
	    bs = d__1 * d__1;
	    c__ = (4 - hk) / 8;
	    d__ = (12 - hk) / 16;
	    bvn = a * exp(-(bs / as + hk) / 2) * (1 - c__ * (bs - as) * (1 - 
		    d__ * bs / 5) / 3 + c__ * d__ * as * as / 5);
	    if (hk > -160.) {
		b = sqrt(bs);
		d__1 = -b / a;
		bvn -= exp(-hk / 2) * sqrt(6.283185307179586) * phi_(&d__1) * 
			b * (1 - c__ * bs * (1 - d__ * bs / 5) / 3);
	    }
	    a /= 2;
	    i__1 = lg;
	    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
		d__1 = a * (x[i__ + ng * 10 - 11] + 1);
		xs = d__1 * d__1;
		rs = sqrt(1 - xs);
		bvn += a * w[i__ + ng * 10 - 11] * (exp(-bs / (xs * 2) - hk / 
			(rs + 1)) / rs - exp(-(bs / xs + hk) / 2) * (c__ * xs 
			* (d__ * xs + 1) + 1));
/* Computing 2nd power */
		d__1 = -x[i__ + ng * 10 - 11] + 1;
		xs = as * (d__1 * d__1) / 4;
		rs = sqrt(1 - xs);
		bvn += a * w[i__ + ng * 10 - 11] * exp(-(bs / xs + hk) / 2) * 
			(exp(-hk * (1 - rs) / ((rs + 1) * 2)) / rs - (c__ * 
			xs * (d__ * xs + 1) + 1));
/* L20: */
	    }
	    bvn = -bvn / 6.283185307179586;
	}
	if (*r__ > 0.) {
	    d__1 = -max(h__,k);
	    bvn += phi_(&d__1);
	}
	if (*r__ < 0.) {
/* Computing MAX */
	    d__3 = -h__;
	    d__4 = -k;
	    d__1 = 0., d__2 = phi_(&d__3) - phi_(&d__4);
	    bvn = -bvn + max(d__1,d__2);
	}
    }
    ret_val = bvn;
    return ret_val;
} /* bvnu_ */

#undef x
#undef w



/* -------------------------------------------------------------------------- */
/* Subroutine */ int ncvsrt_(int *n, double *_lower, double *_upper,
	 double *correl, int *infin, double *y, int *infis, 
	double *a, double *b, int *infi, double *cov, 
	double *d__, double *e)
{
    /* System generated locals */
    int i__1, i__2, i__3;
    double d__1;

    /* Local variables */
    static double amin, bmin, dmin__, emin;
    static int jmin, i__, j, k, iflag;
    static double sumsq, aj, bj;
    static int ii, ij;
    static double cvdiag, yl, yu;
    static double sum;


/*     Subroutine to sort integration limits. */

    /* Parameter adjustments */
    --cov;
    --infi;
    --b;
    --a;
    --y;
    --infin;
    --correl;
    --_upper;
    --_lower;

    /* Function Body */
    ij = 0;
    ii = 0;
    *infis = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	infi[i__] = infin[i__];
	if (infi[i__] < 0) {
	    ++(*infis);
	} else {
	    a[i__] = 0.;
	    b[i__] = 0.;
	    if (infi[i__] != 0) {
		a[i__] = _lower[i__];
	    }
	    if (infi[i__] != 1) {
		b[i__] = _upper[i__];
	    }
	}
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    ++ij;
	    ++ii;
	    cov[ij] = correl[ii];
	}
	++ij;
	cov[ij] = 1.;
    }

/*     First move any doubly infinite limits to innermost positions */
/*     [AA, recoded to avoid GOTO jump outside IF block] */

    if (*infis < *n) {
	i__1 = *n - *infis + 1;
	for (i__ = *n; i__ >= i__1; --i__) {
	    iflag = 0;
	    if (infi[i__] >= 0) {
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
		    if (infi[j] < 0 && iflag == 0) {
			rcswap_(&j, &i__, &a[1], &b[1], &infi[1], n, &cov[1]);
			iflag = 1;
		    }
		}
	    }
/* L10: */
	}

/*     Sort remaining limits and determine Cholesky decomposition */

	ii = 0;
	i__1 = *n - *infis;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*     Determine the integration limits for variable with minimum */
/*      expected probability and interchange that variable with Ith. */

	    emin = 1.;
	    dmin__ = 0.;
	    jmin = i__;
	    cvdiag = 0.;
	    ij = ii;
	    i__2 = *n - *infis;
	    for (j = i__; j <= i__2; ++j) {
		sum = 0.;
		sumsq = 0.;
		i__3 = i__ - 1;
		for (k = 1; k <= i__3; ++k) {
		    sum += cov[ij + k] * y[k];
/* Computing 2nd power */
		    d__1 = cov[ij + k];
		    sumsq += d__1 * d__1;
		}
		ij += j;
/* Computing MAX */
		d__1 = cov[ij] - sumsq;
		sumsq = sqrt((max(d__1,0.)));
		if (sumsq > 0.) {
		    if (infi[j] != 0) {
			aj = (a[j] - sum) / sumsq;
		    }
		    if (infi[j] != 1) {
			bj = (b[j] - sum) / sumsq;
		    }
		    limits_(&aj, &bj, &infi[j], d__, e);
		    if (emin - dmin__ >= *e - *d__) {
			jmin = j;
			if (infi[j] != 0) {
			    amin = aj;
			}
			if (infi[j] != 1) {
			    bmin = bj;
			}
			dmin__ = *d__;
			emin = *e;
			cvdiag = sumsq;
		    }
		}
	    }
	    if (jmin != i__) {
		rcswap_(&i__, &jmin, &a[1], &b[1], &infi[1], n, &cov[1]);
	    }

/*     Compute Ith column of Cholesky factor. */

	    ij = ii + i__;
	    cov[ij] = cvdiag;
	    i__2 = *n - *infis;
	    for (j = i__ + 1; j <= i__2; ++j) {
		if (cvdiag > 0.) {
		    sum = cov[ij + i__];
		    i__3 = i__ - 1;
		    for (k = 1; k <= i__3; ++k) {
			sum -= cov[ii + k] * cov[ij + k];
		    }
		    cov[ij + i__] = sum / cvdiag;
		} else {
		    cov[ij + i__] = 0.;
		}
		ij += j;
	    }

/*     Compute expected value for Ith integration variable and */
/*     scale Ith covariance matrix row and limits. */

	    if (cvdiag > 0.) {
		if (emin > dmin__ + 1e-8) {
		    yl = 0.;
		    yu = 0.;
		    if (infi[i__] != 0) {
/* Computing 2nd power */
			d__1 = amin;
			yl = -exp(-(d__1 * d__1) / 2) / 2.5066282746310005024;
		    }
		    if (infi[i__] != 1) {
/* Computing 2nd power */
			d__1 = bmin;
			yu = -exp(-(d__1 * d__1) / 2) / 2.5066282746310005024;
		    }
		    y[i__] = (yu - yl) / (emin - dmin__);
		} else {
		    if (infi[i__] == 0) {
			y[i__] = bmin;
		    }
		    if (infi[i__] == 1) {
			y[i__] = amin;
		    }
		    if (infi[i__] == 2) {
			y[i__] = (amin + bmin) / 2;
		    }
		}
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    ++ii;
		    cov[ii] /= cvdiag;
		}
		if (infi[i__] != 0) {
		    a[i__] /= cvdiag;
		}
		if (infi[i__] != 1) {
		    b[i__] /= cvdiag;
		}
	    } else {
		y[i__] = 0.;
		ii += i__;
	    }
	}
	limits_(&a[1], &b[1], &infi[1], d__, e);
    }
    return 0;
} /* ncvsrt_ */

/* -------------------------------------------------------------------------- */
/* Subroutine */ int basrul_(int *ndim, double *a, double *b, 
	double *width, D_fp functn, double *w, int *lenrul, 
	double *g, double *center, double *z__, double *
	rgnert, double *basest)
{
    /* System generated locals */
    int w_dim1, w_offset, g_dim1, g_offset, i__1;
    double d__1, d__2;

    /* Local variables */
    static int i__;
    static double rgncmp, rgnval, rgncpt, rgnerr, rgnvol;
    static double fsymsm;


/*     For application of basic integration rule */


/*     Compute Volume and Center of Subregion */

    /* Parameter adjustments */
    --z__;
    --center;
    --width;
    --b;
    --a;
    g_dim1 = *ndim;
    g_offset = 1 + g_dim1 * 1;
    g -= g_offset;
    w_dim1 = *lenrul;
    w_offset = 1 + w_dim1 * 1;
    w -= w_offset;

    /* Function Body */
    rgnvol = 1.;
    i__1 = *ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rgnvol = rgnvol * 2 * width[i__];
	center[i__] = a[i__] + width[i__];
    }
    *basest = 0.;
    *rgnert = 0.;

/*     Compute basic rule and error */

L10:
    rgnval = 0.;
    rgnerr = 0.;
    rgncmp = 0.;
    rgncpt = 0.;
    i__1 = *lenrul;
    for (i__ = 1; i__ <= i__1; ++i__) {
	fsymsm = fulsum_(ndim, &center[1], &width[1], &z__[1], &g[i__ * 
		g_dim1 + 1], (D_fp)functn);
/*     Basic Rule */
	rgnval += w[i__ + w_dim1] * fsymsm;
/*     First comparison rule */
	rgnerr += w[i__ + (w_dim1 << 1)] * fsymsm;
/*     Second comparison rule */
	rgncmp += w[i__ + w_dim1 * 3] * fsymsm;
/*     Third Comparison rule */
	rgncpt += w[i__ + (w_dim1 << 2)] * fsymsm;
    }

/*     Error estimation */

/* Computing 2nd power */
    d__1 = rgncmp;
/* Computing 2nd power */
    d__2 = rgnerr;
    rgnerr = sqrt(d__1 * d__1 + d__2 * d__2);
/* Computing 2nd power */
    d__1 = rgncpt;
/* Computing 2nd power */
    d__2 = rgncmp;
    rgncmp = sqrt(d__1 * d__1 + d__2 * d__2);
    if (rgnerr * 4 < rgncmp) {
	rgnerr /= 2;
    }
    if (rgnerr * 2 > rgncmp) {
	rgnerr = max(rgnerr,rgncmp);
    }
    *rgnert += rgnvol * rgnerr;
    *basest += rgnvol * rgnval;

/*     When subregion has more than one piece, determine next piece and */
/*      loop back to apply basic rule. */

    i__1 = *ndim;
    for (i__ = 1; i__ <= i__1; ++i__) {
	center[i__] += width[i__] * 2;
	if (center[i__] < b[i__]) {
	    goto L10;
	}
	center[i__] = a[i__] + width[i__];
    }
    return 0;
} /* basrul_ */

double fulsum_(int *s, double *center, double *hwidth, 
	double *x, double *g, D_fp f)
{
    /* System generated locals */
    int i__1, i__2;
    double ret_val;

    /* Local variables */
    static int i__, l;
    static double gi, gl;
    static int ixchng, lxchng;
    static double intsum;


/* ***  To compute fully symmetric basic rule sum */

    /* Parameter adjustments */
    --g;
    --x;
    --hwidth;
    --center;

    /* Function Body */
    ret_val = 0.;

/*     Compute centrally symmetric sum for permutation of G */

L10:
    intsum = 0.;
    i__1 = *s;
    for (i__ = 1; i__ <= i__1; ++i__) {
	x[i__] = center[i__] + g[i__] * hwidth[i__];
    }
L20:
    intsum += (*f)(s, &x[1]);
    i__1 = *s;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g[i__] = -g[i__];
	x[i__] = center[i__] + g[i__] * hwidth[i__];
	if (g[i__] < 0.) {
	    goto L20;
	}
    }
    ret_val += intsum;

/*     Find next distinct permuation of G and loop back for next sum */

    i__1 = *s;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (g[i__ - 1] > g[i__]) {
	    gi = g[i__];
	    ixchng = i__ - 1;
	    i__2 = (i__ - 1) / 2;
	    for (l = 1; l <= i__2; ++l) {
		gl = g[l];
		g[l] = g[i__ - l];
		g[i__ - l] = gl;
		if (gl <= gi) {
		    --ixchng;
		}
		if (g[l] > gi) {
		    lxchng = l;
		}
	    }
	    if (g[ixchng] <= gi) {
		ixchng = lxchng;
	    }
	    g[i__] = g[ixchng];
	    g[ixchng] = gi;
	    goto L10;
	}
    }

/*     End loop for permutations of G and associated sums */

/*     Restore original order to G's */

    i__1 = *s / 2;
    for (i__ = 1; i__ <= i__1; ++i__) {
	gi = g[i__];
	g[i__] = g[*s + 1 - i__];
	g[*s + 1 - i__] = gi;
    }
    return ret_val;
} /* fulsum_ */

/* Subroutine */ int differ_(int *ndim, double *a, double *b, 
	double *width, double *z__, double *dif, D_fp functn, 
	int *divaxn, int *difcls)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__;
    static double frthdf, funcen, widthi;


/*     Compute fourth differences and subdivision axes */

    /* Parameter adjustments */
    --dif;
    --z__;
    --width;
    --b;
    --a;

    /* Function Body */
    *difcls = 0;
    *divaxn = *divaxn % *ndim + 1;
    if (*ndim > 1) {
	i__1 = *ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    dif[i__] = 0.;
	    z__[i__] = a[i__] + width[i__];
	}
L10:
	funcen = (*functn)(ndim, &z__[1]);
	i__1 = *ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    widthi = width[i__] / 5;
	    frthdf = funcen * 6;
	    z__[i__] -= widthi * 4;
	    frthdf += (*functn)(ndim, &z__[1]);
	    z__[i__] += widthi * 2;
	    frthdf -= (*functn)(ndim, &z__[1]) * 4;
	    z__[i__] += widthi * 4;
	    frthdf -= (*functn)(ndim, &z__[1]) * 4;
	    z__[i__] += widthi * 2;
	    frthdf += (*functn)(ndim, &z__[1]);
/*     Do not include differences below roundoff */
	    if (funcen + frthdf / 8 != funcen) {
		dif[i__] += fabs(frthdf) * width[i__];
	    }
	    z__[i__] -= widthi * 4;
	}
	*difcls = *difcls + (*ndim << 2) + 1;
	i__1 = *ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    z__[i__] += width[i__] * 2;
	    if (z__[i__] < b[i__]) {
		goto L10;
	    }
	    z__[i__] = a[i__] + width[i__];
	}
	i__1 = *ndim;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (dif[*divaxn] < dif[i__]) {
		*divaxn = i__;
	    }
	}
    }
    return 0;
} /* differ_ */

/* -------- */
/* Subroutine */ int trestr_(int *pointr, int *sbrgns, double *
	pontrs, double *rgners)
{
    double rgnerr;
    int subrgn, pointp, subtmp, points;

/* ***BEGIN PROLOGUE TRESTR */
/* ***PURPOSE TRESTR maintains a heap for subregions. */
/* ***DESCRIPTION TRESTR maintains a heap for subregions. */
/*            The subregions are ordered according to the size of the */
/*            greatest error estimates of each subregion (RGNERS). */

/*   PARAMETERS */

/*     POINTR Integer. */
/*            The index for the subregion to be inserted in the heap. */
/*     SBRGNS Integer. */
/*            Number of subregions in the heap. */
/*     PONTRS Real array of dimension SBRGNS. */
/*            Used to store the indices for the greatest estimated errors */
/*            for each subregion. */
/*     RGNERS Real array of dimension SBRGNS. */
/*            Used to store the greatest estimated errors for each */
/*            subregion. */

/* ***ROUTINES CALLED NONE */
/* ***END PROLOGUE TRESTR */

/*   Global variables. */


/*   Local variables. */

/*   RGNERR Intermediate storage for the greatest error of a subregion. */
/*   SUBRGN Position of child/parent subregion in the heap. */
/*   SUBTMP Position of parent/child subregion in the heap. */


/* ***FIRST PROCESSING STATEMENT TRESTR */

    /* Parameter adjustments */
    --rgners;
    --pontrs;

    /* Function Body */
    rgnerr = rgners[*pointr];
    if ((double) (*pointr) == pontrs[1]) {

/*        Move the new subregion inserted at the top of the heap */
/*        to its correct position in the heap. */

	subrgn = 1;
L10:
	subtmp = subrgn << 1;
	if (subtmp <= *sbrgns) {
	    if (subtmp != *sbrgns) {

/*              Find maximum of left and right child. */

		points = (int) pontrs[subtmp];
		pointp = (int) pontrs[subtmp + 1];
		if (rgners[points] < rgners[pointp]) {
		    ++subtmp;
		}
	    }

/*           Compare maximum child with parent. */
/*           If parent is maximum, then done. */

	    points = (int) pontrs[subtmp];
	    if (rgnerr < rgners[points]) {

/*              Move the pointer at position subtmp up the heap. */

		pontrs[subrgn] = pontrs[subtmp];
		subrgn = subtmp;
		goto L10;
	    }
	}
    } else {

/*        Insert new subregion in the heap. */

	subrgn = *sbrgns;
L20:
	subtmp = subrgn / 2;
	if (subtmp >= 1) {

/*           Compare child with parent. If parent is maximum, then done. */

	    points = (int) pontrs[subtmp];
	    if (rgnerr > rgners[points]) {

/*              Move the pointer at position subtmp down the heap. */

		pontrs[subrgn] = pontrs[subtmp];
		subrgn = subtmp;
		goto L20;
	    }
	}
    }
    pontrs[subrgn] = (double) (*pointr);

/* ***END TRESTR */

    return 0;
} /* trestr_ */

/* -------------------------------------------------------------------------- */
/* Subroutine */ int rcswap_(int *p, int *q, double *a, 
	double *b, int *infin, int *n, double *c__)
{
    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i__, j;
    static double t;
    static int ii, jj;


/*     Swaps rows and columns P and Q in situ. */

    /* Parameter adjustments */
    --c__;
    --infin;
    --b;
    --a;

    /* Function Body */
    t = a[*p];
    a[*p] = a[*q];
    a[*q] = t;
    t = b[*p];
    b[*p] = b[*q];
    b[*q] = t;
    j = infin[*p];
    infin[*p] = infin[*q];
    infin[*q] = j;
    jj = *p * (*p - 1) / 2;
    ii = *q * (*q - 1) / 2;
    t = c__[jj + *p];
    c__[jj + *p] = c__[ii + *q];
    c__[ii + *q] = t;
    i__1 = *p - 1;
    for (j = 1; j <= i__1; ++j) {
	t = c__[jj + j];
	c__[jj + j] = c__[ii + j];
	c__[ii + j] = t;
    }
    jj += *p;
    i__1 = *q - 1;
    for (i__ = *p + 1; i__ <= i__1; ++i__) {
	t = c__[jj + *p];
	c__[jj + *p] = c__[ii + i__];
	c__[ii + i__] = t;
	jj += i__;
    }
    ii += *q;
    i__1 = *n;
    for (i__ = *q + 1; i__ <= i__1; ++i__) {
	t = c__[ii + *p];
	c__[ii + *p] = c__[ii + *q];
	c__[ii + *q] = t;
	ii += i__;
    }
    return 0;
} /* rcswap_ */

/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */
/* -------------------------------------------------------------------------- */

/* Subroutine */ int sadmvt_(int *n, int *nu, double *_lower, 
	double *_upper, int *infin, double *correl, int *
	maxpts, double *abseps, double *releps, double *error, 
	double *value, int *inform__)
{
    /* System generated locals */
    int i__1, i__2, i__3, i__4;
    double d__1, d__2;

    /* Local variables */
	static double work[8000] = { 0, }, d__ = 0, e = 0;
    static int m = 0;
    static int infis = 0;
    static double oldval = 0;
    static int maxcls, newcls;
    static int rulcls, totcls;


/*     A subroutine for computing multivariate t probabilities. */
/*          Alan Genz */
/*          Department of Mathematics */
/*          Washington State University */
/*          Pullman, WA 99164-3113 */
/*          Email : AlanGenz@wsu.edu */

/*  Parameters */

/*     N      INTEGER, the number of variables. */
/*     NU     INTEGER, the number of degrees of freedom. */
/*     LOWER  REAL, array of _lower integration limits. */
/*     UPPER  REAL, array of _upper integration limits. */
/*     INFIN  INTEGER, array of integration limits flags: */
/*            if INFIN(I) < 0, Ith limits are (-infinity, infinity); */
/*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
/*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
/*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
/*     CORREL REAL, array of correlation coefficients; the correlation */
/*            coefficient in row I column J of the correlation matrix */
/*            should be stored in CORREL( J + ((I-2)*(I-1))/2 ), for J < I. */
/*     MAXPTS INTEGER, maximum number of function values allowed. This */
/*            parameter can be used to limit the time taken. A sensible */
/*            strategy is to start with MAXPTS = 1000*N, and then */
/*            increase MAXPTS if ERROR is too large. */
/*     ABSEPS REAL absolute error tolerance. */
/*     RELEPS REAL relative error tolerance. */
/*     ERROR  REAL, estimated absolute error, with 99% confidence level. */
/*     VALUE  REAL, estimated value for the integral */
/*     INFORM INTEGER, termination status parameter: */
/*            if INFORM = 0, normal completion with ERROR < EPS; */
/*            if INFORM = 1, completion with ERROR > EPS and MAXPTS */
/*                           function vaules used; increase MAXPTS to */
/*                           decrease ERROR; */
/*            if INFORM = 2, N > 20 or N < 1. */

    /* Parameter adjustments */
    --correl;
    --infin;
    --_upper;
    --_lower;

    /* Function Body */
    if (*n > 20 || *n < 1) {
	*inform__ = 2;
	*value = 0.;
	*error = 1.;
	return 0;
    }
    *inform__ = (int) mvtnit_(n, nu, &correl[1], &_lower[1], &_upper[1], &
	    infin[1], &infis, &d__, &e);
    m = *n - infis;
    if (m == 0) {
	*value = 1.;
	*error = 0.;
    } else if (m == 1) {
	*value = e - d__;
	*error = (float)2e-16;
    } else {

/*        Call the subregion adaptive integration subroutine */

	--m;
	rulcls = 1;
	adapt_(&m, &rulcls, &c__0, (D_fp)fncmvt_, abseps, releps, &c__8000, 
		work, error, value, inform__);
/* Computing MIN */
	i__1 = rulcls * 10;
	maxcls = min(i__1,*maxpts);
	totcls = 0;
	adapt_(&m, &totcls, &maxcls, (D_fp)fncmvt_, abseps, releps, &c__8000, 
		work, error, value, inform__);
/* Computing MAX */
	d__1 = *abseps, d__2 = *releps * fabs(*value);
	if (*error > max(d__1,d__2)) {
L10:
	    oldval = *value;
/* Computing MAX */
/* Computing MIN */
	    i__3 = maxcls * 3 / 2, i__4 = *maxpts - totcls;
	    i__1 = rulcls << 1, i__2 = min(i__3,i__4);
	    maxcls = max(i__1,i__2);
	    newcls = -1;
	    adapt_(&m, &newcls, &maxcls, (D_fp)fncmvt_, abseps, releps, &
		    c__8000, work, error, value, inform__);
	    totcls += newcls;
/* Computing 2nd power */
	    d__2 = *error;
	    *error = (d__1 = *value - oldval, fabs(d__1)) + sqrt(rulcls * (
		    d__2 * d__2) / totcls);
/* Computing MAX */
	    d__1 = *abseps, d__2 = *releps * fabs(*value);
	    if (*error > max(d__1,d__2)) {
		if (*maxpts - totcls > rulcls << 1) {
		    goto L10;
		}
	    } else {
		*inform__ = 0;
	    }
	}
    }
    return 0;
} /* sadmvt_ */


double fncmvt_0_(int n__, int *n, double *w, int *nuin, 
	double *correl, double *__lower, double *__upper, int *
	infin, int *infis, double *d__, double *e)
{
    /* System generated locals */
    int i__1, i__2;
    double ret_val, d__1, d__2;

    /* Local variables */
	static int infi[20] = { 0, };
    static double prod = 0, a[20] = { 0, }, b[20] = { 0, };
    static int i__ = 0, j = 0;
    static double y[20] = { 0, }, d1 = 0, e1 = 0, di, ei;
    static int ij = 0;
    static double yd, ui = 0;
    static int nu = 0;
    static double cov[210] = { 0, }, sum = 0;


/*     Integrand subroutine */

    /* Parameter adjustments */
    if (w) {
	--w;
	}
    if (correl) {
	--correl;
	}
    if (__lower) {
	--__lower;
	}
    if (__upper) {
	--__upper;
	}
    if (infin) {
	--infin;
	}

    /* Function Body */
    switch(n__) {
	case 1: goto L_mvtnit;
	}

    di = d1;
    ei = e1;
    prod = ei - di;
    ij = 1;
    yd = 1.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = nu + i__ - 1;
	d__1 = di + w[i__] * (ei - di);
	ui = stdinv_(&i__2, &d__1);
	y[i__ - 1] = ui / yd;
	yd /= sqrt((ui - 1) * (ui + 1) / (nu + i__) + 1);
	sum = 0.;
	i__2 = i__;
	for (j = 1; j <= i__2; ++j) {
	    ++ij;
	    sum += cov[ij - 1] * y[j - 1];
	}
	++ij;
	i__2 = nu + i__;
	d__1 = (a[i__] - sum) * yd;
	d__2 = (b[i__] - sum) * yd;
	mvtlms_(&i__2, &d__1, &d__2, &infi[i__], &di, &ei);
	prod *= ei - di;
    }
    ret_val = prod;
    return ret_val;

/*     Entry point for intialization */


L_mvtnit:
    ret_val = 0.;

/*     Initialization and computation of covariance matrix Cholesky factor */

    mvtsrt_(n, nuin, &__lower[1], &__upper[1], &correl[1], &infin[1], y, infis, a,
	     b, infi, cov, d__, e);
    nu = *nuin;
    d1 = *d__;
    e1 = *e;
    return ret_val;
} /* fncmvt_ */

double fncmvt_(int *n, double *w)
{
    return fncmvt_0_(0, n, w, (int *)0, (double *)0, (double *)0, 
	    (double *)0, (int *)0, (int *)0, (double *)0, (
	    double *)0);
    }

double mvtnit_(int *n, int *nuin, double *correl, double *
	_lower, double *_upper, int *infin, int *infis, double *
	d__, double *e)
{
    return fncmvt_0_(1, n, (double *)0, nuin, correl, _lower, _upper, infin,
	     infis, d__, e);
    }

/* Subroutine */ int mvtlms_(int *nu, double *a, double *b, 
	int *infin, double *_lower, double *_upper)
{
    *_lower = 0.;
    *_upper = 1.;
    if (*infin >= 0) {
	if (*infin != 0) {
	    *_lower = studnt_(nu, a);
	}
	if (*infin != 1) {
	    *_upper = studnt_(nu, b);
	}
    }
    return 0;
} /* mvtlms_ */


/* Subroutine */ int mvtsrt_(int *n, int *nu, double *_lower, 
	double *_upper, double *correl, int *infin, double *y, 
	int *infis, double *a, double *b, int *infi, 
	double *cov, double *d__, double *e)
{
    /* System generated locals */
    int i__1, i__2, i__3;
    double d__1, d__2, d__3;

    /* Local variables */
    static double amin, bmin, dmin__, emin;
    static int jmin, i__, j, k, iflag;
    static double sumsq, ai, bi;
    static int ii, ij;
    static double cvdiag, yd, yl, conodd, yu, conevn;
    static double con, sum;


/*     Sort limits */

    /* Parameter adjustments */
    --cov;
    --infi;
    --b;
    --a;
    --y;
    --infin;
    --correl;
    --_upper;
    --_lower;

    /* Function Body */
    ij = 0;
    ii = 0;
    *infis = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	infi[i__] = infin[i__];
	if (infi[i__] < 0) {
	    ++(*infis);
	} else {
	    a[i__] = 0.;
	    b[i__] = 0.;
	    if (infi[i__] != 0) {
		a[i__] = _lower[i__];
	    }
	    if (infi[i__] != 1) {
		b[i__] = _upper[i__];
	    }
	}
	i__2 = i__ - 1;
	for (j = 1; j <= i__2; ++j) {
	    ++ij;
	    ++ii;
	    cov[ij] = correl[ii];
	}
	++ij;
	cov[ij] = 1.;
    }
    conodd = .31830988618379069;
    conevn = .5;
    i__1 = *nu - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (i__ % 2 == 0) {
	    if (i__ > 2) {
		conevn = conevn * (i__ - 1) / (i__ - 2);
	    }
	} else {
	    if (i__ > 2) {
		conodd = conodd * (i__ - 1) / (i__ - 2);
	    }
	}
    }

/*     First move any doubly infinite limits to innermost positions */
/*     [AA, recoded to avoid GOTO jump outside IF block] */

    if (*infis < *n) {
	i__1 = *n - *infis + 1;
	for (i__ = *n; i__ >= i__1; --i__) {
	    iflag = 0;
	    if (infi[i__] >= 0) {
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
		    if (infi[j] < 0 && iflag == 0) {
			rcswap_(&j, &i__, &a[1], &b[1], &infi[1], n, &cov[1]);
			iflag = 1;
		    }
		}
	    }
/* L10: */
	}

/*     Sort remaining limits and determine Cholesky decomposition */

	ii = 0;
	yd = 1.;
	i__1 = *n - *infis;
	for (i__ = 1; i__ <= i__1; ++i__) {

/*     Determine the integration limits for variable with minimum */
/*      expected probability and interchange that variable with Ith. */

	    emin = 1.;
	    dmin__ = 0.;
	    jmin = i__;
	    cvdiag = 0.;
	    ij = ii;
	    i__2 = *n - *infis;
	    for (j = i__; j <= i__2; ++j) {
		sum = 0.;
		sumsq = 0.;
		i__3 = i__ - 1;
		for (k = 1; k <= i__3; ++k) {
		    sum += cov[ij + k] * y[k];
/* Computing 2nd power */
		    d__1 = cov[ij + k];
		    sumsq += d__1 * d__1;
		}
		ij += j;
/* Computing MAX */
		d__1 = cov[ij] - sumsq;
		sumsq = sqrt((max(d__1,0.)));
		if (sumsq > 0.) {
		    ai = yd * (a[j] - sum) / sumsq;
		    bi = yd * (b[j] - sum) / sumsq;
		    i__3 = *nu + j - 1;
		    mvtlms_(&i__3, &ai, &bi, &infi[j], d__, e);
		    if (emin - dmin__ >= *e - *d__) {
			jmin = j;
			amin = ai;
			bmin = bi;
			dmin__ = *d__;
			emin = *e;
			cvdiag = sumsq;
		    }
		}
	    }
	    if (jmin != i__) {
		rcswap_(&i__, &jmin, &a[1], &b[1], &infi[1], n, &cov[1]);
	    }

/*     Compute Ith column of Cholesky factor. */

	    ij = ii + i__;
	    cov[ij] = cvdiag;
	    i__2 = *n - *infis;
	    for (j = i__ + 1; j <= i__2; ++j) {
		if (cvdiag > 0.) {
		    sum = cov[ij + i__];
		    i__3 = i__ - 1;
		    for (k = 1; k <= i__3; ++k) {
			sum -= cov[ii + k] * cov[ij + k];
		    }
		    cov[ij + i__] = sum / cvdiag;
		} else {
		    cov[ij + i__] = 0.;
		}
		ij += j;
	    }

/*     Compute expected value for Ith integration variable and */
/*     scale Ith covariance matrix row and limits. */

	    if ((*nu + i__ - 1) % 2 == 0) {
		if (*nu + i__ - 3 > 0) {
		    conevn = conevn * (*nu + i__ - 2) / (*nu + i__ - 3);
		}
		con = conevn;
	    } else {
		if (*nu + i__ - 3 > 0) {
		    conodd = conodd * (*nu + i__ - 2) / (*nu + i__ - 3);
		}
		con = conodd;
	    }
	    if (cvdiag > 0.) {
		yl = 0.;
		yu = 0.;
		if (infi[i__] != 0 && *nu + i__ - 2 > 0) {
/* Computing 2nd power */
		    d__2 = amin;
		    d__1 = d__2 * d__2 / (*nu + i__ - 1) + 1;
		    d__3 = (*nu + i__ - 2) / 2.;
		    yl = -con * (*nu + i__ - 1) / (*nu + i__ - 2) / pow(d__1, d__3);
		}
		if (infi[i__] != 1 && *nu + i__ - 2 > 0) {
/* Computing 2nd power */
		    d__2 = bmin;
		    d__1 = d__2 * d__2 / (*nu + i__ - 1) + 1;
		    d__3 = (*nu + i__ - 2) / 2.;
		    yu = -con * (*nu + i__ - 1) / (*nu + i__ - 2) / pow(d__1, d__3);
		}
		y[i__] = (yu - yl) / (emin - dmin__) / yd;
		i__2 = i__;
		for (j = 1; j <= i__2; ++j) {
		    ++ii;
		    cov[ii] /= cvdiag;
		}
		if (infi[i__] != 0) {
		    a[i__] /= cvdiag;
		}
		if (infi[i__] != 1) {
		    b[i__] /= cvdiag;
		}
	    } else {
		y[i__] = 0.;
		ii += i__;
	    }
	    yd /= sqrt((y[i__] * yd + 1) * (y[i__] * yd - 1) / (*nu + i__) + 
		    1);
	}
	mvtlms_(nu, &a[1], &b[1], &infi[1], d__, e);
    }
    return 0;
} /* mvtsrt_ */

/* -- */
double studnt_(int *nu, double *t)
{
    /* System generated locals */
    double ret_val;

    /* Local variables */
    static int j;
    static double snthe, polyn, rn, ts, tt, cssthe;


/*     Student t Distribution Function */

/*                       T */
/*         STUDNT = C   I  ( 1 + y*y/NU )**( -(NU+1)/2 ) dy */
/*                   NU -INF */

    if (*nu == 1) {
	ret_val = (atan(*t) * 2 / 3.141592653589793 + 1) / 2;
    } else if (*nu == 2) {
	ret_val = (*t / sqrt(*t * *t + 2) + 1) / 2;
    } else {
	tt = *t * *t;
	cssthe = 1 / (tt / *nu + 1);
	polyn = 1.;
	for (j = *nu - 2; j >= 2; j += -2) {
	    polyn = (j - 1) * cssthe * polyn / j + 1;
	}
	if (*nu % 2 == 1) {
	    rn = (double) (*nu);
	    ts = *t / sqrt(rn);
	    ret_val = ((atan(ts) + ts * cssthe * polyn) * 2 / 
		    3.141592653589793 + 1) / 2;
	} else {
	    snthe = *t / sqrt(*nu + tt);
	    ret_val = (snthe * polyn + 1) / 2;
	}
	ret_val = max(0.,ret_val);
    }
    return ret_val;
} /* studnt_ */

/* -- */
double stdinv_(int *n, double *z__)
{
    /* System generated locals */
    double ret_val, d__1, d__2;

    /* Local variables */
    static double a, b, c__, d__, p, x, y;


/*     Inverse Student t Distribution Function */

/*                    STDINV */
/*           Z = C   I      (1 + y*y/N)**(-(N+1)/2) dy */
/*                N  -INF */

/*      Reference: G.W. Hill, Comm. ACM Algorithm 395 */
/*                 Comm. ACM 13 (1970), pp. 619-620. */

/*      Conversions to double precision and other modifications by */
/*                 Alan Genz, 1993-4. */

    if (0. < *z__ && *z__ < 1.) {
	if (*n == 1) {
	    ret_val = tan((*z__ * 2 - 1) * 3.141592653589793 / 2);
	} else if (*n == 2) {
	    ret_val = (*z__ * 2 - 1) / sqrt(*z__ * 2 * (1 - *z__));
	} else {
	    if (*z__ * 2 >= 1.) {
		p = (1 - *z__) * 2;
	    } else {
		p = *z__ * 2;
	    }
	    a = 1 / (*n - (float).5);
	    b = 48 / (a * a);
	    c__ = ((a * 20700 / b - 98) * a - 16) * a + (float)96.36;
	    d__ = (((float)94.5 / (b + c__) - 3) / b + 1) * sqrt(a * 
		    3.141592653589793 / 2) * *n;
	    x = d__ * p;
	    d__1 = 2. / *n;
	    y = pow(x, d__1);
	    if (y > a + (float).05) {
		d__1 = p / 2;
		x = phinv_(&d__1);
		y = x * x;
		if (*n < 5) {
		    c__ += (*n - (float)4.5) * 3 * (x * 10 + 6) / 100;
		}
		c__ = (((d__ * x - 100) * x / 20 - 7) * x - 2) * x + b + c__;
		y = (((((y * 4 + 63) * y / 10 + 36) * y + (float)94.5) / c__ 
			- y - 3) / b + 1) * x;
		y = a * y * y;
		if (y > (float).002) {
		    y = exp(y) - 1;
		} else {
		    y *= y / 2 + 1;
		}
	    } else {
		y = ((1 / (((*n + 6) / (*n * y) - d__ * (float).089 - (float)
			.822) * (*n * 3 + 6)) + (float).5 / (*n + 4)) * y - 1)
			 * (*n + 1) / (*n + 2) + 1 / y;
	    }
	    ret_val = sqrt(*n * y);
	    if (*z__ * 2 < 1.) {
		ret_val = -ret_val;
	    }
	    if (fabs(ret_val) > 0.) {

/*     Use one third order correction to the single precision result */

		x = ret_val;
		d__ = *z__ - studnt_(n, &x);
		ret_val = x + d__ * 2 / (2 / stdjac_(n, &x) - d__ * (*n + 1) /
			 (*n / x + x));
	    }
	}
    } else {

/*     Use cutoff values for Z near 0 or 1. */

	d__1 = sqrt(*n * 6.2831853071795862) * 2e-16;
	d__2 = 2. / *n;
	ret_val = sqrt(*n / pow(d__1, d__2));
	if (*z__ * 2 < 1.) {
	    ret_val = -ret_val;
	}
    }
    return ret_val;
} /* stdinv_ */

/* -- */
double stdjac_(int *nu, double *t)
{
    /* Initialized data */

    double nuold = 0.;

    /* System generated locals */
    int i__1;
    double ret_val, d__1;

    /* Local variables */
    static int j;
    static double const__, tt;


/*     Student t Distribution Transformation Jacobean */

/*          T            STDINV(NU,T) */
/*         I  f(y) dy = I   f(STDINV(NU,Z) STDJAC(NU,STDINV(NU,Z)) dZ */
/*         -INF          0 */

    if (*nu == 1) {
	ret_val = (*t * *t + 1) * 3.141592653589793;
    } else if (*nu == 2) {
/* Computing 3rd power */
	d__1 = sqrt(*t * *t + 2);
	ret_val = d__1 * (d__1 * d__1);
    } else {
	if ((double) (*nu) != nuold) {
	    nuold = (double) (*nu);
	    if (*nu % 2 == 0) {
		const__ = sqrt(nuold) * 2;
	    } else {
		const__ = sqrt(nuold) * 3.141592653589793;
	    }
	    for (j = *nu - 2; j >= 1; j += -2) {
		const__ = j * const__ / (j + 1);
	    }
	}
	tt = *t * *t / *nu + 1;
	i__1 = (*nu + 1) / 2;
	ret_val = const__ * pow(tt, i__1);
	if (*nu % 2 == 0) {
	    ret_val *= sqrt(tt);
	}
    }
    return ret_val;
} /* stdjac_ */

	/* mnormt/src/biv-nt.f -- translated by f2c (version 20000121).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

/* Table of constant values */

static double c_b8 = 1.;

/* Selected portion of code taken from: */
/*    http://www.math.wsu.edu/faculty/genz/software/mvtdstpack.f */
/* to compute bivariate normal and Student's t distribution functions. */

/* Author: */
/*          Alan Genz */
/*          Department of Mathematics */
/*          Washington State University */
/*          Pullman, WA 99164-3113 */
/*          Email : alangenz@wsu.edu */

/* except for some auxiliary functions whose authors are indicated */
/* in the respective code below. */

/* In addition the dummy subroutine below has been added, needed to interface */
/* R and Fortran.  Adelchi Azzalini, 2008-12-06 */
/* Subroutine */ int smvbvt_(double *prob, int *nu, double *_lower,
	 double *_upper, int *infin, double *correl)
{
    /* Parameter adjustments */
    --infin;
    --_upper;
    --_lower;

    /* Function Body */
    *prob = mvbvt_(nu, &_lower[1], &_upper[1], &infin[1], correl);
    return 0;
} /* smvbvt_ */

/* *********************************************************************** */

double mvbvn_(double *_lower, double *_upper, int *infin, 
			  double *correl)
{
	/* System generated locals */
	double ret_val, d__1, d__2, d__3, d__4;
	/*     A function for computing bivariate normal probabilities. */

	/*  Parameters */

	/*     LOWER  REAL, array of _lower integration limits. */
	/*     UPPER  REAL, array of _upper integration limits. */
	/*     INFIN  INTEGER, array of integration limits flags: */
	/*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
	/*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
	/*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
	/*     CORREL REAL, correlation coefficient. */

	/* Parameter adjustments */
	--infin;
	--_upper;
	--_lower;

	/* Function Body */
	if (infin[1] == 2 && infin[2] == 2) {
		ret_val = mvbvu_(&_lower[1], &_lower[2], correl) - mvbvu_(&_upper[1], &
			_lower[2], correl) - mvbvu_(&_lower[1], &_upper[2], correl) + 
			mvbvu_(&_upper[1], &_upper[2], correl);
	} else if (infin[1] == 2 && infin[2] == 1) {
		ret_val = mvbvu_(&_lower[1], &_lower[2], correl) - mvbvu_(&_upper[1], &
			_lower[2], correl);
	} else if (infin[1] == 1 && infin[2] == 2) {
		ret_val = mvbvu_(&_lower[1], &_lower[2], correl) - mvbvu_(&_lower[1], &
			_upper[2], correl);
	} else if (infin[1] == 2 && infin[2] == 0) {
		d__1 = -_upper[1];
		d__2 = -_upper[2];
		d__3 = -_lower[1];
		d__4 = -_upper[2];
		ret_val = mvbvu_(&d__1, &d__2, correl) - mvbvu_(&d__3, &d__4, correl);
	} else if (infin[1] == 0 && infin[2] == 2) {
		d__1 = -_upper[1];
		d__2 = -_upper[2];
		d__3 = -_upper[1];
		d__4 = -_lower[2];
		ret_val = mvbvu_(&d__1, &d__2, correl) - mvbvu_(&d__3, &d__4, correl);
	} else if (infin[1] == 1 && infin[2] == 0) {
		d__1 = -_upper[2];
		d__2 = -(*correl);
		ret_val = mvbvu_(&_lower[1], &d__1, &d__2);
	} else if (infin[1] == 0 && infin[2] == 1) {
		d__1 = -_upper[1];
		d__2 = -(*correl);
		ret_val = mvbvu_(&d__1, &_lower[2], &d__2);
	} else if (infin[1] == 1 && infin[2] == 1) {
		ret_val = mvbvu_(&_lower[1], &_lower[2], correl);
	} else if (infin[1] == 0 && infin[2] == 0) {
		d__1 = -_upper[1];
		d__2 = -_upper[2];
		ret_val = mvbvu_(&d__1, &d__2, correl);
	} else {
		ret_val = 1.;
	}
	return ret_val;
} /* mvbvn_ */

// double mvbvu_(double *sh, double *sk, double *r__)
// {
// 	/* Initialized data */
// 
// 	static struct {
// 		double e_1[3];
// 		double fill_2[7];
// 		double e_3[6];
// 		double fill_4[4];
// 		double e_5[10];
// 	} equiv_21 = {
// 		{ .1713244923791705, .3607615730481384, .4679139345726904 },
// 		{ 0, },
// 		{ .04717533638651177, .1069393259953183, .1600783285433464, .2031674267230659, .2334925365383547, .2491470458134029 },
// 		{ 0, },
// 		{ .01761400713915212, .04060142980038694, .06267204833410906, .08327674157670475,
// 		.1019301198172404, .1181945319615184, .1316886384491766, 
// 		.1420961093183821, .1491729864726037, .1527533871307259 } };
// 
// #define w ((double *)&equiv_21)
// 
// 	static struct {
// 		double e_1[3];
// 		double fill_2[7];
// 		double e_3[6];
// 		double fill_4[4];
// 		double e_5[10];
// 	} equiv_22 = {
// 		{ -.9324695142031522, -.6612093864662647, -.238619186083197 },
// 		{ 0, },
// 		{ -.9815606342467191, -.904117256370475, -.769902674194305, -.5873179542866171, -.3678314989981802, -.1252334085114692 },
// 		{ 0, },
// 		{ -.9931285991850949, 
// 		-.9639719272779138, -.9122344282513259, -.8391169718222188, 
// 		-.7463319064601508, -.636053680726515, -.5108670019508271, 
// 		-.3737060887154196, -.2277858511416451, -.07652652113349733 } };
// 
// #define x ((double *)&equiv_22)
// 
// 
// 	/* System generated locals */
// 	int i__1;
// 	double ret_val, d__1, d__2, d__3, d__4;
// 
// 	/* Local variables */
// 	static double a, b, c__, d__, h__;
// 	static int i__;
// 	static double k;
// 	static int lg;
// 	static double as;
// 	static int ng;
// 	static double bs, hk, hs, sn, rs, xs, bvn, asr;
// 
// 
// 	/*     A function for computing bivariate normal probabilities; */
// 	/*       developed using */
// 	/*         Drezner, Z. and Wesolowsky, G. O. (1989), */
// 	/*         On the Computation of the Bivariate Normal Integral, */
// 	/*         J. Stat. Comput. Simul.. 35 pp. 101-107. */
// 	/*       with extensive modications for double precisions by */
// 	/*         Alan Genz and Yihong Ge */
// 	/*         Department of Mathematics */
// 	/*         Washington State University */
// 	/*         Pullman, WA 99164-3113 */
// 	/*         Email : alangenz@wsu.edu */
// 
// 	/* BVN - calculate the probability that X is larger than SH and Y is */
// 	/*       larger than SK. */
// 
// 	/* Parameters */
// 
// 	/*   SH  REAL, integration limit */
// 	/*   SK  REAL, integration limit */
// 	/*   R   REAL, correlation coefficient */
// 	/*   LG  INTEGER, number of Gauss Rule Points and Weights */
// 
// 	/*     Gauss Legendre Points and Weights, N =  6 */
// 	/*     Gauss Legendre Points and Weights, N = 12 */
// 	/*     Gauss Legendre Points and Weights, N = 20 */
// 	if (fabs(*r__) < (float).3) {
// 		ng = 1;
// 		lg = 3;
// 	} else if (fabs(*r__) < (float).75) {
// 		ng = 2;
// 		lg = 6;
// 	} else {
// 		ng = 3;
// 		lg = 10;
// 	}
// 	h__ = *sh;
// 	k = *sk;
// 	hk = h__ * k;
// 	bvn = 0.;
// 	if (fabs(*r__) < (float).925) {
// 		hs = (h__ * h__ + k * k) / 2;
// 		asr = asin(*r__);
// 		i__1 = lg;
// 		for (i__ = 1; i__ <= i__1; ++i__) {
// 			sn = sin(asr * (x[i__ + ng * 10 - 11] + 1) / 2);
// 			bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn))
// 				;
// 			sn = sin(asr * (-x[i__ + ng * 10 - 11] + 1) / 2);
// 			bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn))
// 				;
// 		}
// 		d__1 = -h__;
// 		d__2 = -k;
// 		bvn = bvn * asr / 12.566370614359172 + mvphi_(&d__1) * mvphi_(&d__2);
// 	} else {
// 		if (*r__ < 0.) {
// 			k = -k;
// 			hk = -hk;
// 		}
// 		if (fabs(*r__) < 1.) {
// 			as = (1 - *r__) * (*r__ + 1);
// 			a = sqrt(as);
// 			/* Computing 2nd power */
// 			d__1 = h__ - k;
// 			bs = d__1 * d__1;
// 			c__ = (4 - hk) / 8;
// 			d__ = (12 - hk) / 16;
// 			bvn = a * exp(-(bs / as + hk) / 2) * (1 - c__ * (bs - as) * (1 - 
// 				d__ * bs / 5) / 3 + c__ * d__ * as * as / 5);
// 			if (hk > -160.) {
// 				b = sqrt(bs);
// 				d__1 = -b / a;
// 				bvn -= exp(-hk / 2) * sqrt(6.283185307179586) * mvphi_(&d__1) 
// 					* b * (1 - c__ * bs * (1 - d__ * bs / 5) / 3);
// 			}
// 			a /= 2;
// 			i__1 = lg;
// 			for (i__ = 1; i__ <= i__1; ++i__) {
// 				/* Computing 2nd power */
// 				d__1 = a * (x[i__ + ng * 10 - 11] + 1);
// 				xs = d__1 * d__1;
// 				rs = sqrt(1 - xs);
// 				bvn += a * w[i__ + ng * 10 - 11] * (exp(-bs / (xs * 2) - hk / 
// 					(rs + 1)) / rs - exp(-(bs / xs + hk) / 2) * (c__ * xs 
// 					* (d__ * xs + 1) + 1));
// 				/* Computing 2nd power */
// 				d__1 = -x[i__ + ng * 10 - 11] + 1;
// 				xs = as * (d__1 * d__1) / 4;
// 				rs = sqrt(1 - xs);
// 				bvn += a * w[i__ + ng * 10 - 11] * exp(-(bs / xs + hk) / 2) * 
// 					(exp(-hk * (1 - rs) / ((rs + 1) * 2)) / rs - (c__ * 
// 					xs * (d__ * xs + 1) + 1));
// 			}
// 			bvn = -bvn / 6.283185307179586;
// 		}
// 		if (*r__ > 0.) {
// 			d__1 = -max(h__,k);
// 			bvn += mvphi_(&d__1);
// 		}
// 		if (*r__ < 0.) {
// 			/* Computing MAX */
// 			d__3 = -h__;
// 			d__4 = -k;
// 			d__1 = 0., d__2 = mvphi_(&d__3) - mvphi_(&d__4);
// 			bvn = -bvn + max(d__1,d__2);
// 		}
// 	}
// 	ret_val = bvn;
// 	return ret_val;
// } /* mvbvu_ */

#undef x
#undef w



double mvstdt_(int *nu, double *t)
{
	/* System generated locals */
	double ret_val;

	/* Local variables */
	static int j;
	static double csthe, snthe;
	static double polyn, rn, ts, tt;


	/*     Student t Distribution Function */

	/*                       T */
	/*         TSTDNT = C   I  ( 1 + y*y/NU )**( -(NU+1)/2 ) dy */
	/*                   NU -INF */

	if (*nu < 1) {
		ret_val = mvphi_(t);
	} else if (*nu == 1) {
		ret_val = (atan(*t) * 2 / 3.141592653589793 + 1) / 2;
	} else if (*nu == 2) {
		ret_val = (*t / sqrt(*t * *t + 2) + 1) / 2;
	} else {
		tt = *t * *t;
		csthe = *nu / (*nu + tt);
		polyn = 1.;
		for (j = *nu - 2; j >= 2; j += -2) {
			polyn = (j - 1) * csthe * polyn / j + 1;
		}
		if (*nu % 2 == 1) {
			rn = (double) (*nu);
			ts = *t / sqrt(rn);
			ret_val = ((atan(ts) + ts * csthe * polyn) * 2 / 
				3.141592653589793 + 1) / 2;
		} else {
			snthe = *t / sqrt(*nu + tt);
			ret_val = (snthe * polyn + 1) / 2;
		}
		if (ret_val < 0.) {
			ret_val = 0.;
		}
	}
	return ret_val;
} /* mvstdt_ */


double mvbvt_(int *nu, double *_lower, double *_upper, int *
			  infin, double *correl)
{
	/* System generated locals */
	double ret_val, d__1, d__2, d__3, d__4;

	/*     A function for computing bivariate normal and t probabilities. */

	/*  Parameters */

	/*     NU     INTEGER degrees of freedom parameter; NU < 1 gives normal case. */
	/*     LOWER  REAL, array of _lower integration limits. */
	/*     UPPER  REAL, array of _upper integration limits. */
	/*     INFIN  INTEGER, array of integration limits flags: */
	/*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
	/*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
	/*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
	/*     CORREL REAL, correlation coefficient. */

	/* Parameter adjustments */
	--infin;
	--_upper;
	--_lower;

	/* Function Body */
	if (*nu < 1) {
		ret_val = mvbvn_(&_lower[1], &_upper[1], &infin[1], correl);
	} else {
		if (infin[1] == 2 && infin[2] == 2) {
			ret_val = mvbvtl_(nu, &_upper[1], &_upper[2], correl) - mvbvtl_(nu, 
				&_upper[1], &_lower[2], correl) - mvbvtl_(nu, &_lower[1], &
				_upper[2], correl) + mvbvtl_(nu, &_lower[1], &_lower[2], 
				correl);
		} else if (infin[1] == 2 && infin[2] == 1) {
			d__1 = -_lower[1];
			d__2 = -_lower[2];
			d__3 = -_upper[1];
			d__4 = -_lower[2];
			ret_val = mvbvtl_(nu, &d__1, &d__2, correl) - mvbvtl_(nu, &d__3, &
				d__4, correl);
		} else if (infin[1] == 1 && infin[2] == 2) {
			d__1 = -_lower[1];
			d__2 = -_lower[2];
			d__3 = -_lower[1];
			d__4 = -_upper[2];
			ret_val = mvbvtl_(nu, &d__1, &d__2, correl) - mvbvtl_(nu, &d__3, &
				d__4, correl);
		} else if (infin[1] == 2 && infin[2] == 0) {
			ret_val = mvbvtl_(nu, &_upper[1], &_upper[2], correl) - mvbvtl_(nu, 
				&_lower[1], &_upper[2], correl);
		} else if (infin[1] == 0 && infin[2] == 2) {
			ret_val = mvbvtl_(nu, &_upper[1], &_upper[2], correl) - mvbvtl_(nu, 
				&_upper[1], &_lower[2], correl);
		} else if (infin[1] == 1 && infin[2] == 0) {
			d__1 = -_lower[1];
			d__2 = -(*correl);
			ret_val = mvbvtl_(nu, &d__1, &_upper[2], &d__2);
		} else if (infin[1] == 0 && infin[2] == 1) {
			d__1 = -_lower[2];
			d__2 = -(*correl);
			ret_val = mvbvtl_(nu, &_upper[1], &d__1, &d__2);
		} else if (infin[1] == 1 && infin[2] == 1) {
			d__1 = -_lower[1];
			d__2 = -_lower[2];
			ret_val = mvbvtl_(nu, &d__1, &d__2, correl);
		} else if (infin[1] == 0 && infin[2] == 0) {
			ret_val = mvbvtl_(nu, &_upper[1], &_upper[2], correl);
		} else {
			ret_val = 1.;
		}
	}
	return ret_val;
} /* mvbvt_ */


double mvbvtc_(int *nu, double *l, double *u, int *infin, 
			   double *rho)
{
	/* System generated locals */
	double ret_val;

	/* Local variables */
	static double b;
	static int i__;
	static double lw[2], up[2];
	static int inf[2];


	/*     A function for computing complementary bivariate normal and t */
	/*       probabilities. */

	/*  Parameters */

	/*     NU     INTEGER degrees of freedom parameter. */
	/*     L      REAL, array of lower integration limits. */
	/*     U      REAL, array of upper integration limits. */
	/*     INFIN  INTEGER, array of integration limits flags: */
	/*            if INFIN(1) INFIN(2),        then MVBVTC computes */
	/*                 0         0              P( X>U(1), Y>U(2) ) */
	/*                 1         0              P( X<L(1), Y>U(2) ) */
	/*                 0         1              P( X>U(1), Y<L(2) ) */
	/*                 1         1              P( X<L(1), Y<L(2) ) */
	/*                 2         0      P( X>U(1), Y>U(2) ) + P( X<L(1), Y>U(2) ) */
	/*                 2         1      P( X>U(1), Y<L(2) ) + P( X<L(1), Y<L(2) ) */
	/*                 0         2      P( X>U(1), Y>U(2) ) + P( X>U(1), Y<L(2) ) */
	/*                 1         2      P( X<L(1), Y>U(2) ) + P( X<L(1), Y<L(2) ) */
	/*                 2         2      P( X>U(1), Y<L(2) ) + P( X<L(1), Y<L(2) ) */
	/*                               +  P( X>U(1), Y>U(2) ) + P( X<L(1), Y>U(2) ) */

	/*     RHO    REAL, correlation coefficient. */


	/* Parameter adjustments */
	--infin;
	--u;
	--l;

	/* Function Body */
	for (i__ = 1; i__ <= 2; ++i__) {
		if (infin[i__] % 2 == 0) {
			inf[i__ - 1] = 1;
			lw[i__ - 1] = u[i__];
		} else {
			inf[i__ - 1] = 0;
			up[i__ - 1] = l[i__];
		}
	}
	b = mvbvt_(nu, lw, up, inf, rho);
	for (i__ = 1; i__ <= 2; ++i__) {
		if (infin[i__] == 2) {
			inf[i__ - 1] = 0;
			up[i__ - 1] = l[i__];
			b += mvbvt_(nu, lw, up, inf, rho);
		}
	}
	if (infin[1] == 2 && infin[2] == 2) {
		inf[0] = 1;
		lw[0] = u[1];
		b += mvbvt_(nu, lw, up, inf, rho);
	}
	ret_val = b;
	return ret_val;
} /* mvbvtc_ */


double mvbvtl_(int *nu, double *dh, double *dk, double *
			   r__)
{
	/* System generated locals */
	int i__1;
	double ret_val, d__1, d__2, d__3;

	/* Builtin functions */
	double d_sign(double *, double *);

	/* Local variables */
	static double gmph, gmpk, hkrn, qhrk, xnkh, xnhk;
	static int j, hs, ks;
	static double btnckh, btnchk, btpdkh, btpdhk, hkn, hpk, hrk, krh, bvt,
		ors, snu;


	/*     a function for computing bivariate t probabilities. */

	/*       Alan Genz */
	/*       Department of Mathematics */
	/*       Washington State University */
	/*       Pullman, Wa 99164-3113 */
	/*       Email : alangenz@wsu.edu */

	/*    this function is based on the method described by */
	/*        Dunnett, C.W. and M. Sobel, (1954), */
	/*        A bivariate generalization of Student's t-distribution */
	/*        with tables for certain special cases, */
	/*        Biometrika 41, pp. 153-169. */

	/* mvbvtl - calculate the probability that x < dh and y < dk. */

	/* parameters */

	/*   nu number of degrees of freedom */
	/*   dh 1st lower integration limit */
	/*   dk 2nd lower integration limit */
	/*   r   correlation coefficient */

	snu = sqrt((double) (*nu));
	ors = 1 - *r__ * *r__;
	hrk = *dh - *r__ * *dk;
	krh = *dk - *r__ * *dh;
	if (fabs(hrk) + ors > 0.) {
		/* Computing 2nd power */
		d__1 = hrk;
		/* Computing 2nd power */
		d__2 = hrk;
		/* Computing 2nd power */
		d__3 = *dk;
		xnhk = d__1 * d__1 / (d__2 * d__2 + ors * (*nu + d__3 * d__3));
		/* Computing 2nd power */
		d__1 = krh;
		/* Computing 2nd power */
		d__2 = krh;
		/* Computing 2nd power */
		d__3 = *dh;
		xnkh = d__1 * d__1 / (d__2 * d__2 + ors * (*nu + d__3 * d__3));
	} else {
		xnhk = 0.;
		xnkh = 0.;
	}
	d__1 = *dh - *r__ * *dk;
	hs = (int) d_sign(&c_b8, &d__1);
	d__1 = *dk - *r__ * *dh;
	ks = (int) d_sign(&c_b8, &d__1);
	if (*nu % 2 == 0) {
		bvt = atan2(sqrt(ors), -(*r__)) / 6.2831853071795862;
		/* Computing 2nd power */
		d__1 = *dh;
		gmph = *dh / sqrt((*nu + d__1 * d__1) * 16);
		/* Computing 2nd power */
		d__1 = *dk;
		gmpk = *dk / sqrt((*nu + d__1 * d__1) * 16);
		btnckh = atan2(sqrt(xnkh), sqrt(1 - xnkh)) * 2 / 
			3.14159265358979323844;
		btpdkh = sqrt(xnkh * (1 - xnkh)) * 2 / 3.14159265358979323844;
		btnchk = atan2(sqrt(xnhk), sqrt(1 - xnhk)) * 2 / 
			3.14159265358979323844;
		btpdhk = sqrt(xnhk * (1 - xnhk)) * 2 / 3.14159265358979323844;
		i__1 = *nu / 2;
		for (j = 1; j <= i__1; ++j) {
			bvt += gmph * (ks * btnckh + 1);
			bvt += gmpk * (hs * btnchk + 1);
			btnckh += btpdkh;
			btpdkh = (j << 1) * btpdkh * (1 - xnkh) / ((j << 1) + 1);
			btnchk += btpdhk;
			btpdhk = (j << 1) * btpdhk * (1 - xnhk) / ((j << 1) + 1);
			/* Computing 2nd power */
			d__1 = *dh;
			gmph = gmph * ((j << 1) - 1) / ((j << 1) * (d__1 * d__1 / *nu + 1)
				);
			/* Computing 2nd power */
			d__1 = *dk;
			gmpk = gmpk * ((j << 1) - 1) / ((j << 1) * (d__1 * d__1 / *nu + 1)
				);
		}
	} else {
		/* Computing 2nd power */
		d__1 = *dh;
		/* Computing 2nd power */
		d__2 = *dk;
		qhrk = sqrt(d__1 * d__1 + d__2 * d__2 - *r__ * 2 * *dh * *dk + *nu * 
			ors);
		hkrn = *dh * *dk + *r__ * *nu;
		hkn = *dh * *dk - *nu;
		hpk = *dh + *dk;
		bvt = atan2(-snu * (hkn * qhrk + hpk * hkrn), hkn * hkrn - *nu * hpk *
			qhrk) / 6.2831853071795862;
		if (bvt < -1e-15) {
			bvt += 1;
		}
		/* Computing 2nd power */
		d__1 = *dh;
		gmph = *dh / (snu * 6.2831853071795862 * (d__1 * d__1 / *nu + 1));
		/* Computing 2nd power */
		d__1 = *dk;
		gmpk = *dk / (snu * 6.2831853071795862 * (d__1 * d__1 / *nu + 1));
		btnckh = sqrt(xnkh);
		btpdkh = btnckh;
		btnchk = sqrt(xnhk);
		btpdhk = btnchk;
		i__1 = (*nu - 1) / 2;
		for (j = 1; j <= i__1; ++j) {
			bvt += gmph * (ks * btnckh + 1);
			bvt += gmpk * (hs * btnchk + 1);
			btpdkh = ((j << 1) - 1) * btpdkh * (1 - xnkh) / (j << 1);
			btnckh += btpdkh;
			btpdhk = ((j << 1) - 1) * btpdhk * (1 - xnhk) / (j << 1);
			btnchk += btpdhk;
			/* Computing 2nd power */
			d__1 = *dh;
			gmph = (j << 1) * gmph / (((j << 1) + 1) * (d__1 * d__1 / *nu + 1)
				);
			/* Computing 2nd power */
			d__1 = *dk;
			gmpk = (j << 1) * gmpk / (((j << 1) + 1) * (d__1 * d__1 / *nu + 1)
				);
		}
	}
	ret_val = bvt;

	/*     end mvbvtl */

	return ret_val;
} /* mvbvtl_ */



double mvphi_(double *z__)
{
	/* Initialized data */

	static double a[44] = { .610143081923200417926465815756,
		-.434841272712577471828182820888,.176351193643605501125840298123,
		-.060710795609249414860051215825,.017712068995694114486147141191,
		-.004321119385567293818599864968,8.54216676887098678819832055e-4,
		-1.2715509060916274262889394e-4,1.1248167243671189468847072e-5,
		3.13063885421820972630152e-7,-2.70988068537762022009086e-7,
		3.0737622701407688440959e-8,2.515620384817622937314e-9,
		-1.02892992132031912759e-9,2.9944052119949939363e-11,
		2.605178968726693629e-11,-2.634839924171969386e-12,
		-6.43404509890636443e-13,1.12457401801663447e-13,
		1.7281533389986098e-14,-4.264101694942375e-15,
		-5.45371977880191e-16,1.58697607761671e-16,2.0899837844334e-17,
		-5.900526869409e-18,-9.41893387554e-19,2.1497735647e-19,
		4.6660985008e-20,-7.243011862e-21,-2.387966824e-21,1.91177535e-22,
		1.20482568e-22,-6.72377e-25,-5.747997e-24,-4.28493e-25,
		2.44856e-25,4.3793e-26,-8.151e-27,-3.089e-27,9.3e-29,1.74e-28,
		1.6e-29,-8e-30,-2e-30 };

	/* System generated locals */
	double ret_val;

	/* Local variables */
	static double b;
	static int i__;
	static double p, t, bm, bp, xa;


	/*     Normal distribution probabilities accurate to 1d-15. */
	/*     Reference: J.L. Schonfelder, Math Comp 32(1978), pp 1232-1240. */


	xa = fabs(*z__) / 1.414213562373095048801688724209;
	if (xa > 100.) {
		p = 0.;
	} else {
		t = (xa * 8 - 30) / (xa * 4 + 15);
		bm = 0.;
		b = 0.;
		for (i__ = 24; i__ >= 0; --i__) {
			bp = b;
			b = bm;
			bm = t * b - bp + a[i__];
		}
		p = exp(-xa * xa) * (bm - bp) / 4;
	}
	if (*z__ > 0.) {
		p = 1 - p;
	}
	ret_val = p;
	return ret_val;
} /* mvphi_ */

/* pmrnorm */
double pmnorm(wsReal *Ra_X, wsUint N_sz, wsSym Ra_vcov, wsReal *Ra_mean)
{
	double*		Ra_inX		= NULL;
	double**	Ra_inVcov	= NULL;
	double*		Ra_inMean	= NULL;
	sseMalloc(Ra_inX, double, N_sz);
	if (Ra_mean)
		sseMalloc(Ra_inMean, double, N_sz);
#ifndef USE_DBL
	wsAlloc(Ra_inVcov, double*, N_sz);
	for (wsUint i=0 ; i<N_sz ; i++)
		sseMalloc(Ra_inVcov[i], double, i+1);
	for (wsUint i=0 ; i<N_sz ; i++) {
		for (wsUint j=0 ; j<=i ; j++)
			Ra_inVcov[i][j] = (double)Ra_vcov[i][j];
		Ra_inX[i] = (double)Ra_X[i];
		if (Ra_mean)
			Ra_inMean[i] = (double)Ra_mean[i];
	}
#else
	Ra_inVcov = sseSymMatP(N_sz, Ra_vcov);
	memcpy(Ra_inX, Ra_X, sizeof(wsReal)*N_sz);
	if (Ra_mean)
		memcpy(Ra_inMean, Ra_mean, sizeof(wsReal)*N_sz);
#endif
	double *Ra_lb = NULL;
	wsAlloc(Ra_lb, double, N_sz);
	for (wsUint i=0 ; i<N_sz ; i++) Ra_lb[i] = -numeric_limits<double>::infinity();
	int *Na_inf = NULL;
	wsAlloc(Na_inf, int, N_sz);

	/* get sd of varcov matrix */
	double *Ra_sd = NULL;
	wsAlloc(Ra_sd, double, N_sz);
	for (wsUint i=0 ; i<N_sz ; i++)
		Ra_sd[i] = W1 / sqrt(Ra_inVcov[i][i]);

	if (Ra_inMean == NULL) {
		/* upper <- as.double(upper/sd) */
		for (wsUint i=0 ; i<N_sz ; i++) {
			Ra_inX[i] *= Ra_sd[i];
			Na_inf[i] = Ra_inX[i] > -9.0 ? 2 : 1;
		}
	} else {
		/* upper <- as.double((upper-mean)/sd) */
		for (wsUint i=0 ; i<N_sz ; i++) {
			Ra_inX[i] = (Ra_inX[i] - Ra_inMean[i]) * Ra_sd[i];
			Na_inf[i] = Ra_inX[i] > -9.0 ? 2 : 1;
		}
	}

	wsUint N_sig = (N_sz * (N_sz-1)) >> 1;
	double* Ra_sig = NULL;
	wsAlloc(Ra_sig, double, N_sig);
	for (wsUint i=1,k=0 ; i<N_sz ; i++)
		for (wsUint j=0 ; j<i ; j++,k++)
			Ra_sig[k] = Ra_inVcov[i][j] * Ra_sd[i] * Ra_sd[j];
	double R_prob = 1e-6;
	double R_releps = 1e-5;
	double R_error = 1e-6;
	double R_bound = 0.5 * 1e-6;
	int N_fault = 0;
	int N_tsz = N_sz;
	int N_maxpts = 2000 * N_sz;
	sadmvn_(&N_tsz, Ra_lb, Ra_inX, Na_inf, Ra_sig, &N_maxpts, &R_bound, &R_releps,
		&R_error, &R_prob, &N_fault);
// 	mulnor(Ra_inX, Ra_lb, Ra_sig, 1e-6, N_sz, Na_inf, &R_prob, &R_bound, 
// 		   &N_fault);

	sseFree(Ra_inX);
	sseFree(Ra_mean);
	for (wsUint i=0 ; i<N_sz ; i++)
		sseFree(Ra_inVcov[i]);
	DEALLOC(Ra_inVcov);

	/* Return probability if N_fault is valid */
	return N_fault < 2 ? R_prob : 0.0;
}

/* Subroutine */ int mvtdst_(int *n, int *nu, double *_lower, 
	double *_upper, int *infin, double *correl, double *
	delta, int *maxpts, double *abseps, double *releps, 
	double *error, double *value, int *inform__);

/* pmrnorm */
double pmnorm(wsReal R_X, wsUint N_sz, wsSym Ra_cor)
{
	double*		Ra_inX		= NULL;
	double**	Ra_inCor	= NULL;
	sseMalloc(Ra_inX, double, N_sz);
#ifndef USE_DBL
	wsAlloc(Ra_inCor, double*, N_sz);
	for (wsUint i=0 ; i<N_sz ; i++)
		sseMalloc(Ra_inCor[i], double, i+1);
	for (wsUint i=0 ; i<N_sz ; i++) {
		for (wsUint j=0 ; j<=i ; j++)
			Ra_inCor[i][j] = (double)Ra_cor[i][j];
		Ra_inX[i] = (double)R_X;
	}
#else
	Ra_inCor = sseSymMatP(N_sz, Ra_cor);
	sseVinit(Ra_inX, N_sz, R_X);
#endif
	double *Ra_lb = NULL;
	wsAlloc(Ra_lb, double, N_sz);
	for (wsUint i=0 ; i<N_sz ; i++) Ra_lb[i] = -R_X;
	int *Na_inf = NULL;
	wsAlloc(Na_inf, int, N_sz);

	/* upper <- as.double(upper/sd) */
	for (wsUint i=0 ; i<N_sz ; i++)
		Na_inf[i] = 2;

	wsUint N_sig = (N_sz * (N_sz-1)) >> 1;
	double* Ra_sig = NULL;
	wsAlloc(Ra_sig, double, N_sig);
	for (wsUint i=1,k=0 ; i<N_sz ; i++)
		for (wsUint j=0 ; j<i ; j++,k++)
			Ra_sig[k] = Ra_inCor[i][j];
	double R_prob = 1e-6;
	double R_releps = 0;
	double R_error = 1e-6;
	double R_abseps = 0.001;
	int N_fault = 0;
	int N_tsz = N_sz;
	int N_maxpts = 25000;
#if 0
	sadmvn_(&N_tsz, Ra_lb, Ra_inX, Na_inf, Ra_sig, &N_maxpts, &R_bound, &R_releps,
		&R_error, &R_prob, &N_fault);
#else
	int df = 0;
	double *delta = new double[N_sz];
	for (wsUint i=0 ; i<N_sz ; i++) delta[i] = 0;

	mvtdst_(&N_tsz, &df, Ra_lb, Ra_inX, Na_inf, Ra_sig, delta, &N_maxpts,
		&R_abseps, &R_releps, &R_error, &R_prob, &N_fault);
#endif
	DEALLOC(Ra_lb);
	// 	mulnor(Ra_inX, Ra_lb, Ra_sig, 1e-6, N_sz, Na_inf, &R_prob, &R_bound, 
	// 		   &N_fault);
	DEALLOC(Na_inf);

	sseFree(Ra_inX);
	for (wsUint i=0 ; i<N_sz ; i++)
		sseFree(Ra_inCor[i]);
	DEALLOC(Ra_inCor);
	DEALLOC(Ra_sig);

	/* Return probability if N_fault is valid */
	return N_fault < 2 ? R_prob : 0.0;
}



/*** mvtnorm ***/

double d_mod(double *x, double *y)
{
	double quotient;
	if( (quotient = *x / *y) >= 0)
		quotient = floor(quotient);
	else
		quotient = -floor(-quotient);
	return(*x - (*y) * quotient );
}

/* Subroutine */ int mvints_(int *n, int *nuin, double *correl, 
							 double *_lower, double *_upper, double *delta, int *
							 infin, int *nd, double *vl, double *er, int *inform__);
/* Subroutine */ int mvkbrv_(int *ndim, int *minvls, int *maxvls, 
							 int *nf, D_fp funsub, double *abseps, double *releps, 
							 double *abserr, double *finest, int *inform__);
/* Subroutine */ int mvkbrv_(int *ndim, int *minvls, int *maxvls, 
							 int *nf, D_fp funsub, double *abseps, double *releps, 
							 double *abserr, double *finest, int *inform__);
/* Subroutine */ int mvsubr_(int *n, double *w, int *nf, 
							 double *f);
/* Subroutine */ int mvkrsv_(int *ndim, int *kl, double *values, 
							 int *prime, double *vk, int *nf, D_fp funsub, double *
							 x, double *r__, int *pr, double *fs);
double mvchnc_(double *lkn, int *n, double *p, double *
			   r__);
/* Subroutine */ int mvvlsb_(int *n, double *w, double *r__, 
							 double *dl, int *infi, double *a, double *b, 
							 double *cov, double *y, double *di, double *ei, 
							 int *nd, double *value);
double mvchnv_(int *n, double *p);
/* Subroutine */ int mvsort_(int *n, double *_lower, double *_upper,
							 double *delta, double *correl, int *infin, double *y,
							 bool *pivot, int *nd, double *a, double *b, 
							 double *dl, double *cov, int *infi, int *inform__);
/* Subroutine */ int mvspcl_(int *nd, int *nu, double *a, 
							 double *b, double *dl, double *cov, int *infi, 
							 double *snu, double *vl, double *er, int *inform__);
/* Subroutine */ int mvlims_(double *a, double *b, int *infin, 
							 double *_lower, double *_upper);
double mvphnv_(double *p);
/* Subroutine */ int mvswap_(int *p, int *q, double *a, 
							 double *b, double *d__, int *infin, int *n, 
							 double *c__);
double mvtdns_(int *nu, double *x);
/* Subroutine */ int mvsswp_(double *x, double *y);

struct _WISARD_XXXXX_1_ {
	int ivls;
} ptblck_;

#define ptblck_1 ptblck_

/* Table of constant values */

static int c__1 = 1;
static bool c_true = true;
static double c_b24 = 1.;
static double c_b32 = 2.;
static int c__100 = 100;


/*    $Id: mvt.f 245 2013-05-29 14:10:53Z thothorn $ */

/* Subroutine */ int mvtdst_(int *n, int *nu, double *_lower, 
							 double *_upper, int *infin, double *correl, double *
							 delta, int *maxpts, double *abseps, double *releps, 
							 double *error, double *value, int *inform__)
{
	static double e[1], v[1];
	static int nd;

	/*     A subroutine for computing non-central multivariate t probabilities. */
	/*     This subroutine uses an algorithm (QRSVN) described in the paper */
	/*     "Comparison of Methods for the Computation of Multivariate */
	/*         t-Probabilities", by Alan Genz and Frank Bretz */
	/*         J. Comp. Graph. Stat. 11 (2002), pp. 950-971. */

	/*          Alan Genz */
	/*          Department of Mathematics */
	/*          Washington State University */
	/*          Pullman, WA 99164-3113 */
	/*          Email : AlanGenz@wsu.edu */

	/* 	Original source available from */
	/* 	http://www.math.wsu.edu/faculty/genz/software/fort77/mvtdstpack.f */

	/* 	This is version 28/05/2013 */

	/*  Parameters */

	/*     N      INTEGER, the number of variables. */
	/*     NU     INTEGER, the number of degrees of freedom. */
	/*            If NU < 1, then an MVN probability is computed. */
	/*     LOWER  DOUBLE PRECISION, array of _lower integration limits. */
	/*     UPPER  DOUBLE PRECISION, array of _upper integration limits. */
	/*     INFIN  INTEGER, array of integration limits flags: */
	/*             if INFIN(I) < 0, Ith limits are (-infinity, infinity); */
	/*             if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)]; */
	/*             if INFIN(I) = 1, Ith limits are [LOWER(I), infinity); */
	/*             if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)]. */
	/*     CORREL DOUBLE PRECISION, array of correlation coefficients; */
	/*            the correlation coefficient in row I column J of the */
	/*            correlation matrixshould be stored in */
	/*               CORREL( J + ((I-2)*(I-1))/2 ), for J < I. */
	/*            The correlation matrix must be positive semi-definite. */
	/*     DELTA  DOUBLE PRECISION, array of non-centrality parameters. */
	/*     MAXPTS INTEGER, maximum number of function values allowed. This */
	/*            parameter can be used to limit the time. A sensible */
	/*            strategy is to start with MAXPTS = 1000*N, and then */
	/*            increase MAXPTS if ERROR is too large. */
	/*     ABSEPS DOUBLE PRECISION absolute error tolerance. */
	/*     RELEPS DOUBLE PRECISION relative error tolerance. */
	/*     ERROR  DOUBLE PRECISION estimated absolute error, */
	/*            with 99% confidence level. */
	/*     VALUE  DOUBLE PRECISION estimated value for the integral */
	/*     INFORM INTEGER, termination status parameter: */
	/*            if INFORM = 0, normal completion with ERROR < EPS; */
	/*            if INFORM = 1, completion with ERROR > EPS and MAXPTS */
	/*                           function vaules used; increase MAXPTS to */
	/*                           decrease ERROR; */
	/*            if INFORM = 2, N > 1000 or N < 1. */
	/*            if INFORM = 3, correlation matrix not positive semi-definite. */

	/* Parameter adjustments */
	--delta;
	--correl;
	--infin;
	--_upper;
	--_lower;

	/* Function Body */
	ptblck_1.ivls = 0;
	//rndstart_();
	if (*n > 1000 || *n < 1) {
		*value = 0.;
		*error = 1.;
		*inform__ = 2;
	} else {
		mvints_(n, nu, &correl[1], &_lower[1], &_upper[1], &delta[1], &infin[1],
			&nd, value, error, inform__);
		if (*inform__ == 0 && nd > 0) {

			/*           Call the lattice rule integration subroutine */

			mvkbrv_(&nd, &ptblck_1.ivls, maxpts, &c__1, (D_fp)mvsubr_, abseps,
				releps, e, v, inform__);
			*error = e[0];
			*value = v[0];
		}
	}
	//rndend_();
	return 0;
} /* mvtdst_ */


/* Subroutine */ int mvsubr_0_(int n__, int *n, double *w, int *
							   nf, double *f, int *nuin, double *correl, double *
							   _lower, double *_upper, double *delta, int *infin, int *
							   nd, double *vl, double *er, int *inform__)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	static double a[1000], b[1000], r__, y[1000], di, ei, dl[1000];
	static int nu, ny;
	static double cov[500500], snu;
	static int infi[1000];


	/*     Integrand subroutine */

	/* Parameter adjustments */
	if (w) {
		--w;
	}
	if (f) {
		--f;
	}
	if (correl) {
		--correl;
	}
	if (_lower) {
		--_lower;
	}
	if (_upper) {
		--_upper;
	}
	if (delta) {
		--delta;
	}
	if (infin) {
		--infin;
	}

	/* Function Body */
	switch(n__) {
	case 1: goto L_mvints;
	}

	if (nu <= 0) {
		r__ = 1.;
		i__1 = *n + 1;
		mvvlsb_(&i__1, &w[1], &r__, dl, infi, a, b, cov, y, &di, &ei, &ny, &f[
			1]);
	} else {
		r__ = mvchnv_(&nu, &w[*n]) / snu;
		mvvlsb_(n, &w[1], &r__, dl, infi, a, b, cov, y, &di, &ei, &ny, &f[1]);
	}
	return 0;

	/*     Entry point for intialization. */


L_mvints:

	/*     Initialization and computation of covariance Cholesky factor. */

	mvsort_(n, &_lower[1], &_upper[1], &delta[1], &correl[1], &infin[1], y, &
		c_true, nd, a, b, dl, cov, infi, inform__);
	nu = *nuin;
	mvspcl_(nd, &nu, a, b, dl, cov, infi, &snu, vl, er, inform__);
	return 0;
} /* mvsubr_ */

/* Subroutine */ int mvsubr_(int *n, double *w, int *nf, 
							 double *f)
{
	return mvsubr_0_(0, n, w, nf, f, (int *)0, (double *)0, (
		double *)0, (double *)0, (double *)0, (int *)0, (
		int *)0, (double *)0, (double *)0, (int *)0);
}

/* Subroutine */ int mvints_(int *n, int *nuin, double *correl, 
							 double *_lower, double *_upper, double *delta, int *
							 infin, int *nd, double *vl, double *er, int *inform__)
{
	return mvsubr_0_(1, n, (double *)0, (int *)0, (double *)0, 
		nuin, correl, _lower, _upper, delta, infin, nd, vl, er, inform__);
}


/* Subroutine */ int mvspcl_(int *nd, int *nu, double *a, 
							 double *b, double *dl, double *cov, int *infi, 
							 double *snu, double *vl, double *er, int *inform__)
{
	/* System generated locals */
	double d__1;

	/* Local variables */
	static double r__;


	/*     Special cases subroutine */

	/* Parameter adjustments */
	--infi;
	--cov;
	--dl;
	--b;
	--a;

	/* Function Body */
	if (*inform__ > 0) {
		*vl = 0.;
		*er = 1.;
	} else {

		/*        Special cases */

		if (*nd == 0) {
			*er = 0.;
			/*  Code added to fix ND = 0 bug, 24/03/2009 -> */
			*vl = 1.;
			/*  <- Code added to fix ND = 0 bug, 24/03/2009 */
		} else if (*nd == 1 && (*nu < 1 || fabs(dl[1]) == 0.)) {

			/*           1-d case for normal or central t */

			*vl = 1.;
			if (infi[1] != 1) {
				d__1 = b[1] - dl[1];
				*vl = mvstdt_(nu, &d__1);
			}
			if (infi[1] != 0) {
				d__1 = a[1] - dl[1];
				*vl -= mvstdt_(nu, &d__1);
			}
			if (*vl < 0.) {
				*vl = 0.;
			}
			*er = 2e-16;
			*nd = 0;
		} else if (*nd == 2 && (*nu < 1 || fabs(dl[1]) + fabs(dl[2]) == 0.)) {

			/*           2-d case for normal or central t */

			if (infi[1] != 0) {
				a[1] -= dl[1];
			}
			if (infi[1] != 1) {
				b[1] -= dl[1];
			}
			if (infi[2] != 0) {
				a[2] -= dl[2];
			}
			if (infi[2] != 1) {
				b[2] -= dl[2];
			}
			if (fabs(cov[3]) > 0.) {

				/*              2-d nonsingular case */

				/* Computing 2nd power */
				d__1 = cov[2];
				r__ = sqrt(d__1 * d__1 + 1);
				if (infi[2] != 0) {
					a[2] /= r__;
				}
				if (infi[2] != 1) {
					b[2] /= r__;
				}
				cov[2] /= r__;
				*vl = mvbvt_(nu, &a[1], &b[1], &infi[1], &cov[2]);
				*er = 1e-15;
			} else {

				/*              2-d singular case */

				if (infi[1] != 0) {
					if (infi[2] != 0) {
						a[1] = max(a[1],a[2]);
					}
				} else {
					if (infi[2] != 0) {
						a[1] = a[2];
					}
				}
				if (infi[1] != 1) {
					if (infi[2] != 1) {
						b[1] = min(b[1],b[2]);
					}
				} else {
					if (infi[2] != 1) {
						b[1] = b[2];
					}
				}
				if (infi[1] != infi[2]) {
					infi[1] = 2;
				}
				*vl = 1.;
				/*  A(1), B(1) Bug Fixed, 28/05/2013 */
				if (infi[1] != 1) {
					*vl = mvstdt_(nu, &b[1]);
				}
				if (infi[1] != 0) {
					*vl -= mvstdt_(nu, &a[1]);
				}
				if (*vl < 0.) {
					*vl = 0.;
				}
				*er = 2e-16;
			}
			*nd = 0;
		} else {
			if (*nu > 0) {
				*snu = sqrt((double) (*nu));
			} else {
				--(*nd);
			}
		}
	}
	return 0;
} /* mvspcl_ */


/* Subroutine */ int mvvlsb_(int *n, double *w, double *r__, 
							 double *dl, int *infi, double *a, double *b, 
							 double *cov, double *y, double *di, double *ei, 
							 int *nd, double *value)
{
	/* System generated locals */
	int i__1, i__2;
	double d__1, d__2;

	/* Local variables */
	static int i__, j;
	static double ai, bi;
	static int ij;
	static double sum;
	static int infa, infb;


	/*     Integrand subroutine */

	/* Parameter adjustments */
	--y;
	--cov;
	--b;
	--a;
	--infi;
	--dl;
	--w;

	/* Function Body */
	*value = 1.;
	infa = 0;
	infb = 0;
	*nd = 0;
	ij = 0;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		sum = dl[i__];
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
			++ij;
			if (j <= *nd) {
				sum += cov[ij] * y[j];
			}
		}
		if (infi[i__] != 0) {
			if (infa == 1) {
				/* Computing MAX */
				d__1 = ai, d__2 = *r__ * a[i__] - sum;
				ai = max(d__1,d__2);
			} else {
				ai = *r__ * a[i__] - sum;
				infa = 1;
			}
		}
		if (infi[i__] != 1) {
			if (infb == 1) {
				/* Computing MIN */
				d__1 = bi, d__2 = *r__ * b[i__] - sum;
				bi = min(d__1,d__2);
			} else {
				bi = *r__ * b[i__] - sum;
				infb = 1;
			}
		}
		++ij;
		if (i__ == *n || cov[ij + *nd + 2] > 0.) {
			i__2 = infa + infa + infb - 1;
			mvlims_(&ai, &bi, &i__2, di, ei);
			if (*di >= *ei) {
				*value = 0.;
				return 0;
			} else {
				*value *= *ei - *di;
				++(*nd);
				if (i__ < *n) {
					d__1 = *di + w[*nd] * (*ei - *di);
					y[*nd] = mvphnv_(&d__1);
				}
				infa = 0;
				infb = 0;
			}
		}
	}
	return 0;
} /* mvvlsb_ */


/* Subroutine */ int mvsort_(int *n, double *_lower, double *_upper,
							 double *delta, double *correl, int *infin, double *y,
							 bool *pivot, int *nd, double *a, double *b, 
							 double *dl, double *cov, int *infi, int *inform__)
{
	/* System generated locals */
	int i__1, i__2, i__3, i__4;
	double d__1;

	/* Local variables */
	static double d__, e;
	static int i__, j, k, l, m;
	static double aj, bj;
	static int ii, ij, il, jl;
	static double sum, amin, bmin;
	static int jmin;
	static double epsi, demin, sumsq, cvdiag;


	/*     Subroutine to sort integration limits and determine Cholesky factor. */

	/* Parameter adjustments */
	--infi;
	--cov;
	--dl;
	--b;
	--a;
	--y;
	--infin;
	--correl;
	--delta;
	--_upper;
	--_lower;

	/* Function Body */
	*inform__ = 0;
	ij = 0;
	ii = 0;
	*nd = *n;
	i__1 = *n;
	for (i__ = 1; i__ <= i__1; ++i__) {
		a[i__] = 0.;
		b[i__] = 0.;
		dl[i__] = 0.;
		infi[i__] = infin[i__];
		if (infi[i__] < 0) {
			--(*nd);
		} else {
			if (infi[i__] != 0) {
				a[i__] = _lower[i__];
			}
			if (infi[i__] != 1) {
				b[i__] = _upper[i__];
			}
			dl[i__] = delta[i__];
		}
		i__2 = i__ - 1;
		for (j = 1; j <= i__2; ++j) {
			++ij;
			++ii;
			cov[ij] = correl[ii];
		}
		++ij;
		cov[ij] = 1.;
	}

	/*     First move any doubly infinite limits to innermost positions. */

	if (*nd > 0) {
		i__1 = *nd + 1;
		for (i__ = *n; i__ >= i__1; --i__) {
			if (infi[i__] >= 0) {
				i__2 = i__ - 1;
				for (j = 1; j <= i__2; ++j) {
					if (infi[j] < 0) {
						mvswap_(&j, &i__, &a[1], &b[1], &dl[1], &infi[1], n, &
							cov[1]);
						goto L10;
					}
				}
			}
L10:
			;
		}

		/*     Sort remaining limits and determine Cholesky factor. */

		ii = 0;
		jl = *nd;
		i__1 = *nd;
		for (i__ = 1; i__ <= i__1; ++i__) {

			/*        Determine the integration limits for variable with minimum */
			/*        expected probability and interchange that variable with Ith. */

			demin = 1.;
			jmin = i__;
			cvdiag = 0.;
			ij = ii;
			epsi = i__ * 1e-10;
			if (! (*pivot)) {
				jl = i__;
			}
			i__2 = jl;
			for (j = i__; j <= i__2; ++j) {
				if (cov[ij + j] > epsi) {
					sumsq = sqrt(cov[ij + j]);
					sum = dl[j];
					i__3 = i__ - 1;
					for (k = 1; k <= i__3; ++k) {
						sum += cov[ij + k] * y[k];
					}
					aj = (a[j] - sum) / sumsq;
					bj = (b[j] - sum) / sumsq;
					mvlims_(&aj, &bj, &infi[j], &d__, &e);
					if (demin >= e - d__) {
						jmin = j;
						amin = aj;
						bmin = bj;
						demin = e - d__;
						cvdiag = sumsq;
					}
				}
				ij += j;
			}
			if (jmin > i__) {
				mvswap_(&i__, &jmin, &a[1], &b[1], &dl[1], &infi[1], n, &cov[
					1]);
			}
			if (cov[ii + i__] < -epsi) {
				*inform__ = 3;
			}
			cov[ii + i__] = cvdiag;

			/*        Compute Ith column of Cholesky factor. */
			/*        Compute expected value for Ith integration variable and */
			/*         scale Ith covariance matrix row and limits. */

			if (cvdiag > 0.) {
				il = ii + i__;
				i__2 = *nd;
				for (l = i__ + 1; l <= i__2; ++l) {
					cov[il + i__] /= cvdiag;
					ij = ii + i__;
					i__3 = l;
					for (j = i__ + 1; j <= i__3; ++j) {
						cov[il + j] -= cov[il + i__] * cov[ij + i__];
						ij += j;
					}
					il += l;
				}

				/*              Expected Y = -( density(b) - density(a) )/( b - a ) */

				if (demin > epsi) {
					y[i__] = 0.;
					if (infi[i__] != 0) {
						y[i__] = mvtdns_(&c__0, &amin);
					}
					if (infi[i__] != 1) {
						y[i__] -= mvtdns_(&c__0, &bmin);
					}
					y[i__] /= demin;
				} else {
					if (infi[i__] == 0) {
						y[i__] = bmin;
					}
					if (infi[i__] == 1) {
						y[i__] = amin;
					}
					if (infi[i__] == 2) {
						y[i__] = (amin + bmin) / 2;
					}
				}
				i__2 = i__;
				for (j = 1; j <= i__2; ++j) {
					++ii;
					cov[ii] /= cvdiag;
				}
				a[i__] /= cvdiag;
				b[i__] /= cvdiag;
				dl[i__] /= cvdiag;
			} else {
				il = ii + i__;
				i__2 = *nd;
				for (l = i__ + 1; l <= i__2; ++l) {
					cov[il + i__] = 0.;
					il += l;
				}

				/*        If the covariance matrix diagonal entry is zero, */
				/*         permute limits and rows, if necessary. */


				for (j = i__ - 1; j >= 1; --j) {
					if ((d__1 = cov[ii + j], fabs(d__1)) > epsi) {
						a[i__] /= cov[ii + j];
						b[i__] /= cov[ii + j];
						dl[i__] /= cov[ii + j];
						if (cov[ii + j] < 0.) {
							mvsswp_(&a[i__], &b[i__]);
							if (infi[i__] != 2) {
								infi[i__] = 1 - infi[i__];
							}
						}
						i__2 = j;
						for (l = 1; l <= i__2; ++l) {
							cov[ii + l] /= cov[ii + j];
						}
						i__2 = i__ - 1;
						for (l = j + 1; l <= i__2; ++l) {
							if (cov[(l - 1) * l / 2 + j + 1] > 0.) {
								ij = ii;
								i__3 = l;
								for (k = i__ - 1; k >= i__3; --k) {
									i__4 = k;
									for (m = 1; m <= i__4; ++m) {
										mvsswp_(&cov[ij - k + m], &cov[ij + m]
										);
									}
									mvsswp_(&a[k], &a[k + 1]);
									mvsswp_(&b[k], &b[k + 1]);
									mvsswp_(&dl[k], &dl[k + 1]);
									m = infi[k];
									infi[k] = infi[k + 1];
									infi[k + 1] = m;
									ij -= k;
								}
								goto L20;
							}
						}
						goto L20;
					}
					cov[ii + j] = 0.;
				}
L20:
				ii += i__;
				y[i__] = 0.;
			}
		}
	}
	return 0;
} /* mvsort_ */


double mvtdns_(int *nu, double *x)
{
	/* System generated locals */
	int i__1;
	double ret_val, d__1;

	/* Local variables */
	static int i__;
	static double prod;

	ret_val = 0.;
	if (*nu > 0) {
		prod = 1 / sqrt((double) (*nu));
		for (i__ = *nu - 2; i__ >= 1; i__ += -2) {
			prod = prod * (i__ + 1) / i__;
		}
		if (*nu % 2 == 0) {
			prod /= 2;
		} else {
			prod /= 3.141592653589793;
		}
		d__1 = sqrt(*x * *x / *nu + 1);
		i__1 = *nu + 1;
		ret_val = prod / pow(d__1, i__1);
	} else {
		if (fabs(*x) < 10.) {
			ret_val = exp(-(*x) * *x / 2) / 2.506628274631001;
		}
	}
	return ret_val;
} /* mvtdns_ */


/* Subroutine */ int mvlims_(double *a, double *b, int *infin, 
							 double *_lower, double *_upper)
{
	extern double mvphi_(double *);

	*_lower = 0.;
	*_upper = 1.;
	if (*infin >= 0) {
		if (*infin != 0) {
			*_lower = mvphi_(a);
		}
		if (*infin != 1) {
			*_upper = mvphi_(b);
		}
	}
	*_upper = max(*_upper,*_lower);
	return 0;
} /* mvlims_ */


/* Subroutine */ int mvsswp_(double *x, double *y)
{
	static double t;

	t = *x;
	*x = *y;
	*y = t;
	return 0;
} /* mvsswp_ */


/* Subroutine */ int mvswap_(int *p, int *q, double *a, 
							 double *b, double *d__, int *infin, int *n, 
							 double *c__)
{
	/* System generated locals */
	int i__1;

	/* Local variables */
	static int i__, j, ii, jj;


	/*     Swaps rows and columns P and Q in situ, with P <= Q. */

	/* Parameter adjustments */
	--c__;
	--infin;
	--d__;
	--b;
	--a;

	/* Function Body */
	mvsswp_(&a[*p], &a[*q]);
	mvsswp_(&b[*p], &b[*q]);
	mvsswp_(&d__[*p], &d__[*q]);
	j = infin[*p];
	infin[*p] = infin[*q];
	infin[*q] = j;
	jj = *p * (*p - 1) / 2;
	ii = *q * (*q - 1) / 2;
	mvsswp_(&c__[jj + *p], &c__[ii + *q]);
	i__1 = *p - 1;
	for (j = 1; j <= i__1; ++j) {
		mvsswp_(&c__[jj + j], &c__[ii + j]);
	}
	jj += *p;
	i__1 = *q - 1;
	for (i__ = *p + 1; i__ <= i__1; ++i__) {
		mvsswp_(&c__[jj + *p], &c__[ii + i__]);
		jj += i__;
	}
	ii += *q;
	i__1 = *n;
	for (i__ = *q + 1; i__ <= i__1; ++i__) {
		mvsswp_(&c__[ii + *p], &c__[ii + *q]);
		ii += i__;
	}
	return 0;
} /* mvswap_ */

double mvphnv_(double *p)
{
	/* System generated locals */
	double ret_val, d__1, d__2;

	/* Local variables */
	static double q, r__;


	/* 	ALGORITHM AS241  APPL. STATIST. (1988) VOL. 37, NO. 3 */

	/* 	Produces the normal deviate Z corresponding to a given _lower */
	/* 	tail area of P. */

	/* 	The hash sums below are the sums of the mantissas of the */
	/* 	coefficients.   They are included for use in checking */
	/* 	transcription. */


	/*     Coefficients for P close to 0.5 */

	/*     HASH SUM AB    55.88319 28806 14901 4439 */

	/*     Coefficients for P not close to 0, 0.5 or 1. */

	/*     HASH SUM CD    49.33206 50330 16102 89036 */

	/* 	Coefficients for P near 0 or 1. */

	/*     HASH SUM EF    47.52583 31754 92896 71629 */

	q = (*p * 2 - 1) / 2;
	if (fabs(q) <= .425) {
		r__ = .180625 - q * q;
		ret_val = q * (((((((r__ * 2509.0809287301226727 + 
			33430.575583588128105) * r__ + 67265.770927008700853) * r__ + 
			45921.953931549871457) * r__ + 13731.693765509461125) * r__ + 
			1971.5909503065514427) * r__ + 133.14166789178437745) * r__ + 
			3.387132872796366608) / (((((((r__ * 5226.495278852854561 + 
			28729.085735721942674) * r__ + 39307.89580009271061) * r__ + 
			21213.794301586595867) * r__ + 5394.1960214247511077) * r__ + 
			687.1870074920579083) * r__ + 42.313330701600911252) * r__ + 
			1);
	} else {
		/* Computing MIN */
		d__1 = *p, d__2 = 1 - *p;
		r__ = min(d__1,d__2);
		if (r__ > 0.) {
			r__ = sqrt(-log(r__));
			if (r__ <= 5.) {
				r__ += -1.6;
				ret_val = (((((((r__ * 7.7454501427834140764e-4 + 
					.0227238449892691845833) * r__ + 
					.24178072517745061177) * r__ + 1.27045825245236838258)
					* r__ + 3.64784832476320460504) * r__ + 
					5.7694972214606914055) * r__ + 4.6303378461565452959) 
					* r__ + 1.42343711074968357734) / (((((((r__ * 
					1.05075007164441684324e-9 + 5.475938084995344946e-4) *
					r__ + .0151986665636164571966) * r__ + 
					.14810397642748007459) * r__ + .68976733498510000455) 
					* r__ + 1.6763848301838038494) * r__ + 
					2.05319162663775882187) * r__ + 1);
			} else {
				r__ += -5.;
				ret_val = (((((((r__ * 2.01033439929228813265e-7 + 
					2.71155556874348757815e-5) * r__ + 
					.0012426609473880784386) * r__ + 
					.026532189526576123093) * r__ + .29656057182850489123)
					* r__ + 1.7848265399172913358) * r__ + 
					5.4637849111641143699) * r__ + 6.6579046435011037772) 
					/ (((((((r__ * 2.04426310338993978564e-15 + 
					1.4215117583164458887e-7) * r__ + 
					1.8463183175100546818e-5) * r__ + 
					7.868691311456132591e-4) * r__ + 
					.0148753612908506148525) * r__ + 
					.13692988092273580531) * r__ + .59983220655588793769) 
					* r__ + 1);
			}
		} else {
			ret_val = 9.;
		}
		if (q < 0.) {
			ret_val = -ret_val;
		}
	}
	return ret_val;
} /* mvphnv_ */

double mvbvu_(double *sh, double *sk, double *r__)
{
	/* Initialized data */

	static struct {
		double e_1[3];
		double fill_2[7];
		double e_3[6];
		double fill_4[4];
		double e_5[10];
	} equiv_83 = {
		{ .1713244923791705, .3607615730481384, .4679139345726904 },
		{ 0, },
		{ .04717533638651177, .1069393259953183, .1600783285433464, .2031674267230659, .2334925365383547, .2491470458134029 },
		{ 0, },
		{ .01761400713915212, 
		.04060142980038694, .06267204833410906, .08327674157670475, 
		.1019301198172404, .1181945319615184, .1316886384491766, 
		.1420961093183821, .1491729864726037, .1527533871307259 } };

#define w ((double *)&equiv_83)

	static struct {
		double e_1[3];
		double fill_2[7];
		double e_3[6];
		double fill_4[4];
		double e_5[10];
	} equiv_84 = {
		{ -.9324695142031522, -.6612093864662647, -.238619186083197 },
		{ 0, },
		{ -.9815606342467191, -.904117256370475, -.769902674194305, -.5873179542866171, -.3678314989981802, -.1252334085114692 },
		{ 0, },
		{ -.9931285991850949, 
		-.9639719272779138, -.9122344282513259, -.8391169718222188, 
		-.7463319064601508, -.636053680726515, -.5108670019508271, 
		-.3737060887154196, -.2277858511416451, -.07652652113349733 } };

#define x ((double *)&equiv_84)


	/* System generated locals */
	int i__1;
	double ret_val, d__1, d__2;

	/* Local variables */
	static double a, b, c__, d__, h__;
	static int i__;
	static double k;
	static int lg;
	static double as;
	static int ng;
	static double bs, hk, hs, sn, rs, xs, bvn, asr;


	/*     A function for computing bivariate normal probabilities; */
	/*       developed using */
	/*         Drezner, Z. and Wesolowsky, G. O. (1989), */
	/*         On the Computation of the Bivariate Normal Integral, */
	/*         J. Stat. Comput. Simul.. 35 pp. 101-107. */
	/*       with extensive modications for double precisions by */
	/*         Alan Genz and Yihong Ge */
	/*         Department of Mathematics */
	/*         Washington State University */
	/*         Pullman, WA 99164-3113 */
	/*         Email : alangenz@wsu.edu */

	/* BVN - calculate the probability that X is larger than SH and Y is */
	/*       larger than SK. */

	/* Parameters */

	/*   SH  REAL, integration limit */
	/*   SK  REAL, integration limit */
	/*   R   REAL, correlation coefficient */
	/*   LG  INTEGER, number of Gauss Rule Points and Weights */

	/*     Gauss Legendre Points and Weights, N =  6 */
	/*     Gauss Legendre Points and Weights, N = 12 */
	/*     Gauss Legendre Points and Weights, N = 20 */
	if (fabs(*r__) < .3f) {
		ng = 1;
		lg = 3;
	} else if (fabs(*r__) < .75f) {
		ng = 2;
		lg = 6;
	} else {
		ng = 3;
		lg = 10;
	}
	h__ = *sh;
	k = *sk;
	hk = h__ * k;
	bvn = 0.;
	if (fabs(*r__) < .925f) {
		hs = (h__ * h__ + k * k) / 2;
		asr = asin(*r__);
		i__1 = lg;
		for (i__ = 1; i__ <= i__1; ++i__) {
			sn = sin(asr * (x[i__ + ng * 10 - 11] + 1) / 2);
			bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn))
				;
			sn = sin(asr * (-x[i__ + ng * 10 - 11] + 1) / 2);
			bvn += w[i__ + ng * 10 - 11] * exp((sn * hk - hs) / (1 - sn * sn))
				;
		}
		d__1 = -h__;
		d__2 = -k;
		bvn = bvn * asr / 12.566370614359172 + mvphi_(&d__1) * mvphi_(&d__2);
	} else {
		if (*r__ < 0.) {
			k = -k;
			hk = -hk;
		}
		if (fabs(*r__) < 1.) {
			as = (1 - *r__) * (*r__ + 1);
			a = sqrt(as);
			/* Computing 2nd power */
			d__1 = h__ - k;
			bs = d__1 * d__1;
			c__ = (4 - hk) / 8;
			d__ = (12 - hk) / 16;
			bvn = a * exp(-(bs / as + hk) / 2) * (1 - c__ * (bs - as) * (1 - 
				d__ * bs / 5) / 3 + c__ * d__ * as * as / 5);
			if (hk > -160.) {
				b = sqrt(bs);
				d__1 = -b / a;
				bvn -= exp(-hk / 2) * sqrt(6.283185307179586) * mvphi_(&d__1) 
					* b * (1 - c__ * bs * (1 - d__ * bs / 5) / 3);
			}
			a /= 2;
			i__1 = lg;
			for (i__ = 1; i__ <= i__1; ++i__) {
				/* Computing 2nd power */
				d__1 = a * (x[i__ + ng * 10 - 11] + 1);
				xs = d__1 * d__1;
				rs = sqrt(1 - xs);
				bvn += a * w[i__ + ng * 10 - 11] * (exp(-bs / (xs * 2) - hk / 
					(rs + 1)) / rs - exp(-(bs / xs + hk) / 2) * (c__ * xs 
					* (d__ * xs + 1) + 1));
				/* Computing 2nd power */
				d__1 = -x[i__ + ng * 10 - 11] + 1;
				xs = as * (d__1 * d__1) / 4;
				rs = sqrt(1 - xs);
				/* Computing 2nd power */
				d__1 = rs + 1;
				bvn += a * w[i__ + ng * 10 - 11] * exp(-(bs / xs + hk) / 2) * 
					(exp(-hk * xs / (d__1 * d__1 * 2)) / rs - (c__ * xs * 
					(d__ * xs + 1) + 1));
			}
			bvn = -bvn / 6.283185307179586;
		}
		if (*r__ > 0.) {
			d__1 = -max(h__,k);
			bvn += mvphi_(&d__1);
		} else {
			bvn = -bvn;
			if (k > h__) {
				if (h__ < 0.) {
					bvn = bvn + mvphi_(&k) - mvphi_(&h__);
				} else {
					d__1 = -h__;
					d__2 = -k;
					bvn = bvn + mvphi_(&d__1) - mvphi_(&d__2);
				}
			}
		}
	}
	ret_val = bvn;
	return ret_val;
} /* mvbvu_ */

#undef x
#undef w

double mvchnv_(int *n, double *p)
{
	/* Initialized data */

	static int no = 0;

	/* System generated locals */
	double ret_val, d__1;

	/* Local variables */
	static int i__;
	static double r__, ro, lkn;


	/*                  MVCHNV */
	/*     P =  1 - K  I     exp(-t*t/2) t**(N-1) dt, for N >= 1. */
	/*               N  0 */

	/*                 LRP =   LOG( SQRT( 2/PI ) ) */
	if (*n <= 1) {
		d__1 = *p / 2;
		r__ = -mvphnv_(&d__1);
	} else if (*p < 1.) {
		if (*n == 2) {
			r__ = sqrt(log(*p) * -2);
		} else {
			if (*n != no) {
				no = *n;
				lkn = 0.;
				for (i__ = *n - 2; i__ >= 2; i__ += -2) {
					lkn -= log((double) i__);
				}
				if (*n % 2 == 1) {
					lkn += -.22579135264472743235;
				}
			}
			if ((double) (*n) >= log(1 - *p) * -5 / 4) {
				r__ = 2. / (*n * 9);
				/* Computing 3rd power */
				d__1 = -mvphnv_(p) * sqrt(r__) + 1 - r__;
				r__ = *n * (d__1 * (d__1 * d__1));
				if (r__ > (double) ((*n << 1) + 6)) {
					r__ = (lkn - log(*p)) * 2 + (*n - 2) * log(r__);
				}
			} else {
				r__ = exp((log((1 - *p) * *n) - lkn) * 2. / *n);
			}
			r__ = sqrt(r__);
			ro = r__;
			r__ = mvchnc_(&lkn, n, p, &r__);
			if ((d__1 = r__ - ro, fabs(d__1)) > 1e-6) {
				ro = r__;
				r__ = mvchnc_(&lkn, n, p, &r__);
				if ((d__1 = r__ - ro, fabs(d__1)) > 1e-6) {
					r__ = mvchnc_(&lkn, n, p, &r__);
				}
			}
		}
	} else {
		r__ = 0.;
	}
	ret_val = r__;
	return ret_val;
} /* mvchnv_ */


double mvchnc_(double *lkn, int *n, double *p, double *
			   r__)
{
	/* System generated locals */
	double ret_val, d__1;

	/* Local variables */
	static int i__;
	static double df, ai, bi, al, ci, di, dl, rn, rr, chi;


	/*     Third order Schroeder correction to R for MVCHNV */

	/*                 LRP =   LOG( SQRT( 2/PI ) ) */
	rr = *r__ * *r__;
	if (*n < 2) {
		d__1 = -(*r__);
		chi = mvphi_(&d__1) * 2;
	} else if (*n < 100) {

		/*        Use standard Chi series */

		rn = 1.;
		for (i__ = *n - 2; i__ >= 2; i__ += -2) {
			rn = rr * rn / i__ + 1;
		}
		rr /= 2;
		if (*n % 2 == 0) {
			chi = exp(log(rn) - rr);
		} else {
			d__1 = -(*r__);
			chi = exp(log(*r__ * rn) - .22579135264472743235 - rr) + mvphi_(&
				d__1) * 2;
		}
	} else {
		rr /= 2;
		al = *n / 2.;
		chi = exp(-rr + al * log(rr) + *lkn + log(2.) * (*n - 2) / 2);
		if (rr < al + 1) {

			/*           Use Incomplete Gamma series */

			dl = chi;
			for (i__ = 1; i__ <= 1000; ++i__) {
				dl = dl * rr / (al + i__);
				chi += dl;
				if ((d__1 = dl * rr / (al + i__ + 1 - rr), fabs(d__1)) < 1e-14)
				{
					goto L10;
				}
			}
L10:
			chi = 1 - chi / al;
		} else {

			/*           Use Incomplete Gamma continued fraction */

			bi = rr + 1 - al;
			ci = 1e14;
			di = bi;
			chi /= bi;
			for (i__ = 1; i__ <= 250; ++i__) {
				ai = i__ * (al - i__);
				bi += 2;
				ci = bi + ai / ci;
				if (ci == 0.) {
					ci = 1e-14;
				}
				di = bi + ai / di;
				if (di == 0.) {
					di = 1e-14;
				}
				dl = ci / di;
				chi *= dl;
				if ((d__1 = dl - 1, fabs(d__1)) < 1e-14) {
					goto L20;
				}
			}
		}
	}
L20:
	df = (*p - chi) / exp(*lkn + (*n - 1) * log(*r__) - rr);
	ret_val = *r__ - df * (1 - df * (*r__ - (*n - 1) / *r__) / 2);
	return ret_val;
} /* mvchnc_ */


/* Subroutine */ int mvkbrv_(int *ndim, int *minvls, int *maxvls, 
							 int *nf, D_fp funsub, double *abseps, double *releps, 
							 double *abserr, double *finest, int *inform__)
{
	/* Initialized data */

	static int p[28] = { 31,47,73,113,173,263,397,593,907,1361,2053,3079,
		4621,6947,10427,15641,23473,35221,52837,79259,118891,178349,
		267523,401287,601943,902933,1354471,2031713 };
	static int c__[2772]	/* was [28][99] */ = { 12,13,27,35,64,111,163,
		246,347,505,794,1189,1763,2872,4309,6610,9861,10327,19540,34566,
		31929,40701,103650,165843,130365,333459,500884,858339,9,11,28,27,
		66,42,154,189,402,220,325,888,1018,3233,3758,6977,3647,7582,19926,
		9579,49367,69087,125480,90647,236711,375354,566009,918142,9,17,10,
		27,28,54,83,242,322,601,960,259,1500,1534,4034,1686,4073,7124,
		11582,12654,10982,77576,59978,59925,110235,102417,399251,501970,
		13,10,11,36,28,118,43,102,418,644,528,1082,432,2941,1963,3819,
		2535,8214,11113,26856,3527,64590,46875,189541,125699,383544,
		652979,234813,12,15,11,22,44,20,82,250,215,612,247,725,1332,2910,
		730,2314,3430,9600,24585,37873,27066,39397,77172,67647,56483,
		292630,355008,460565,12,15,20,29,44,31,92,250,220,160,247,811,
		2203,393,642,5647,9865,10271,8726,38806,13226,33179,83021,74795,
		93735,41147,430235,31996,12,15,11,29,55,31,150,102,339,206,338,
		636,126,1796,1502,3953,2830,10193,17218,29501,56010,10858,126904,
		68365,234469,374614,328722,753018,12,15,11,20,67,72,59,250,339,
		206,366,965,2240,919,2246,3614,9328,10800,419,17271,18911,38935,
		14541,167485,60549,48032,670680,256150,12,15,28,45,10,17,76,280,
		339,206,847,497,1719,446,3834,5115,4320,9086,4918,3663,40574,
		43129,56299,143918,1291,435453,405585,199809,12,15,13,5,10,94,76,
		118,337,422,753,497,1284,919,1511,423,5913,2365,4918,10763,20767,
		35468,43636,74912,93937,281493,405585,993599,12,22,13,5,10,14,47,
		196,218,134,753,1490,878,919,1102,423,10365,4409,4918,18955,20767,
		35468,11655,167289,245291,358168,424646,245149,12,15,28,5,10,14,
		11,118,315,518,236,1490,1983,1117,1102,5408,8272,13812,15701,1298,
		9686,5279,52680,75517,196061,114121,670180,794183,3,15,13,21,10,
		11,11,191,315,134,334,392,266,103,1522,7426,3706,5661,17710,26560,
		47603,61518,88549,8148,258647,346892,670180,121349,3,6,13,21,10,
		14,100,215,315,134,334,1291,266,103,1522,423,6186,9344,4037,17132,
		47603,61518,29804,172106,162489,238990,641587,150619,3,6,13,21,38,
		14,131,121,315,518,461,508,266,103,3427,423,7806,9344,4037,17132,
		11736,27945,101894,126159,176631,317313,215580,376952,12,6,14,21,
		38,14,116,121,167,652,711,508,266,103,3427,487,7806,10362,15808,
		4753,11736,70975,113675,35867,204895,164158,59048,809123,7,15,14,
		21,10,94,116,49,167,382,652,1291,747,103,3928,6227,7806,9344,
		11401,4753,41601,70975,48040,35867,73353,35497,633320,809123,7,15,
		14,21,10,10,116,49,167,206,381,1291,747,103,915,2660,8610,9344,
		19398,8713,12888,86478,113675,35867,172319,70530,81010,804319,12,
		9,14,21,10,10,116,49,167,158,381,508,127,103,915,6227,2563,8585,
		25950,18624,32948,86478,34987,121694,28881,70530,20789,67352,12,
		13,14,21,10,10,116,49,361,441,381,1291,127,2311,3818,1221,11558,
		11114,25950,13082,30801,20514,48308,52171,136787,434839,389250,
		969594,12,2,14,21,10,10,116,49,201,179,652,508,2074,3117,3818,
		3811,11558,13080,4454,6791,44243,20514,97926,95354,122081,24754,
		389250,434796,12,2,14,21,49,14,138,49,124,441,381,508,127,1101,
		3818,197,9421,13080,24987,1122,53351,73178,5475,113969,122081,
		24754,638764,969594,12,2,14,21,49,14,138,49,124,56,381,867,2074,
		3117,3818,4367,1181,13080,11719,19363,53351,73178,49449,113969,
		275993,24754,638764,804319,12,13,14,21,49,14,138,49,124,559,381,
		867,1400,3117,4782,351,9421,6949,8697,34695,16016,43098,6850,
		76304,64673,393656,389250,391368,12,11,14,21,49,14,138,49,124,559,
		381,867,1383,1101,4782,1281,1181,3436,1452,18770,35086,43098,
		62545,123709,211587,118711,389250,761041,12,11,14,21,49,14,138,49,
		124,56,381,867,1383,1101,4782,1221,1181,3436,1452,18770,35086,
		4701,62545,123709,211587,118711,398094,754049,12,10,14,21,49,14,
		138,49,124,56,381,934,1383,1101,3818,351,1181,3436,1452,18770,
		32581,59979,9440,144615,211587,148227,80846,466264,3,15,14,21,49,
		14,138,49,124,56,381,867,1383,1101,4782,351,9421,13213,1452,18770,
		2464,59979,33242,123709,282859,271087,147776,754049,3,15,14,29,49,
		11,138,171,124,56,226,867,1383,1101,3818,351,1181,6130,1452,15628,
		2464,58556,9440,64958,282859,355831,147776,754049,3,15,14,17,49,
		11,138,171,124,56,326,867,1383,2503,3818,7245,1181,6130,8697,
		18770,49554,69916,33242,64958,211587,91034,296177,466264,12,15,14,
		17,49,11,101,171,124,56,326,867,1383,2503,1327,1984,10574,8159,
		8697,18770,2464,15170,9440,32377,242821,417029,398094,754049,7,15,
		31,17,49,8,101,171,124,56,326,867,1383,2503,1327,2999,10574,8159,
		6436,18770,2464,15170,33242,193002,256865,417029,398094,754049,7,
		15,31,17,49,8,101,171,231,56,326,867,1383,2503,1327,2999,3534,
		11595,21475,18770,49554,4832,9440,193002,256865,91034,147776,
		282852,12,15,5,17,38,8,101,171,231,56,326,867,1383,2503,1327,2999,
		3534,8159,6436,33766,49554,4832,62850,25023,256865,91034,147776,
		429907,12,15,5,17,38,8,101,171,90,56,326,1284,1400,2503,1327,2999,
		3534,3436,22913,20837,2464,43064,9440,40017,122203,417029,396313,
		390017,12,15,5,17,31,8,101,171,90,56,326,1284,1383,2503,1327,2999,
		3534,7096,6434,20837,81,71685,9440,141605,291915,91034,578233,
		276645,12,6,31,17,4,8,101,171,90,56,126,1284,1383,2503,1327,2999,
		3534,7096,18497,20837,27260,4832,9440,189165,122203,299843,578233,
		994856,12,6,13,17,4,8,101,171,90,56,326,1284,1383,429,1387,3995,
		2898,7096,11089,20837,10681,15170,90308,189165,291915,299843,
		578233,250142,12,6,11,17,31,18,101,171,90,56,326,1284,1383,429,
		1387,2063,2898,7096,11089,20837,2185,15170,90308,141605,291915,
		413548,19482,144595,12,15,11,23,64,18,101,171,90,101,326,1284,
		1383,429,1387,2063,2898,7096,11089,20837,2185,15170,90308,189165,
		122203,413548,620706,907454,12,15,11,23,4,18,101,171,90,101,326,
		1284,1383,429,1387,2063,3450,7096,11089,6545,2185,27679,47904,
		189165,25639,308300,187095,689648,12,9,11,23,4,18,101,171,90,56,
		326,1284,1383,429,1387,2063,2141,7096,3036,6545,2185,27679,47904,
		141605,25639,413548,620706,687580,3,13,11,23,4,18,101,171,90,101,
		326,1284,507,429,1387,1644,2141,7096,3036,6545,2185,27679,47904,
		141605,291803,413548,187095,687580,3,2,11,23,64,113,101,171,90,
		101,326,563,1073,429,1387,2063,2141,7096,14208,6545,2185,60826,
		47904,141605,245397,413548,126467,687580,3,2,13,23,45,62,101,171,
		90,101,326,563,1073,1702,1387,2077,2141,7096,14208,6545,2185,
		60826,47904,189165,284047,308300,241663,687580,12,2,13,23,45,62,
		101,171,90,101,326,563,1073,1702,1387,2512,2141,7096,14208,12138,
		18086,6187,47904,127047,245397,308300,241663,978368,7,13,13,23,45,
		45,101,171,90,101,326,563,1073,1702,2339,2512,2141,7096,14208,
		12138,18086,6187,47904,127047,245397,308300,241663,687580,7,11,13,
		23,45,45,101,171,90,101,195,1010,1990,184,2339,2512,2141,7096,
		12906,12138,18086,4264,47904,127047,245397,413548,241663,552742,
		12,11,13,23,45,113,101,171,48,101,195,1010,1990,184,2339,2077,
		7055,7096,12906,12138,18086,4264,47904,127047,245397,308300,
		241663,105195,12,10,13,23,45,113,101,171,48,101,55,1010,1990,184,
		2339,2077,7055,7096,12906,12138,18086,4264,41143,127047,245397,
		308300,241663,942843,12,15,13,23,66,113,101,171,48,193,55,208,
		1990,184,2339,2077,7055,7096,12906,12138,17631,4264,41143,127047,
		245397,308300,241663,768249,12,15,14,21,66,113,116,171,48,193,55,
		838,1990,184,2339,2077,7055,7096,12906,12138,17631,4264,41143,
		127047,245397,308300,241663,307142,12,15,14,27,66,113,116,171,90,
		193,55,563,507,105,2339,754,7055,7096,12906,12138,18086,45567,
		41143,127047,94241,308300,241663,307142,12,15,14,3,66,113,116,171,
		90,193,55,563,507,105,2339,754,7055,4377,12906,12138,18086,32269,
		41143,127047,66575,15311,241663,307142,12,15,14,3,66,113,116,171,
		90,193,55,563,507,105,2339,754,7055,7096,12906,12138,18086,32269,
		41143,127047,66575,15311,241663,307142,12,15,14,3,66,113,116,171,
		90,193,55,759,507,105,2339,754,7055,4377,7614,12138,37335,32269,
		41143,127047,217673,15311,241663,880619,12,15,14,24,66,113,116,
		171,90,193,55,759,507,105,2339,754,7055,4377,7614,12138,37774,
		32269,36114,127047,217673,15311,321632,880619,3,15,14,27,66,113,
		100,171,90,101,55,564,507,105,2339,754,7055,4377,7614,12138,37774,
		62060,36114,127047,217673,176255,23210,880619,3,15,14,27,66,113,
		100,171,90,101,55,759,507,105,2339,754,7055,4377,7614,12138,37774,
		62060,36114,127047,217673,176255,23210,880619,3,6,14,17,66,113,
		100,171,90,101,55,759,507,105,3148,754,7055,4377,5021,30483,26401,
		62060,36114,127047,217673,23613,394484,880619,12,6,14,29,66,113,
		100,171,90,101,55,801,507,105,3148,754,7055,5410,5021,30483,26401,
		62060,36114,127047,217673,23613,394484,880619,7,6,14,29,66,113,
		100,171,90,101,55,801,1073,105,3148,754,7055,5410,5021,30483,
		26401,62060,24997,127047,217673,23613,394484,880619,7,15,14,29,66,
		113,138,161,90,101,55,801,1073,105,3148,754,7055,4377,5021,30483,
		26401,62060,65162,127047,217673,23613,78101,117185,12,15,14,17,66,
		113,138,161,90,101,55,801,1073,105,3148,754,2831,4377,5021,30483,
		26401,62060,65162,127047,217673,23613,78101,117185,12,9,14,5,66,
		113,138,161,90,101,55,759,1073,105,3148,754,8204,4377,5021,12138,
		26401,62060,65162,127047,217673,23613,78101,117185,12,13,14,5,66,
		63,138,161,90,101,55,759,1073,105,3148,754,8204,4377,10145,12138,
		26401,62060,65162,127785,217673,172210,542095,117185,12,2,14,5,66,
		63,138,161,90,101,55,759,1073,105,3148,754,8204,4377,10145,12138,
		26401,1803,65162,127785,217673,204328,542095,117185,12,2,31,5,66,
		53,101,161,90,101,55,759,1073,105,3148,754,8204,4377,10145,12138,
		26401,1803,65162,127785,217673,204328,542095,117185,12,2,31,21,66,
		63,101,161,90,101,195,759,1073,105,3148,754,8204,4377,10145,12138,
		26401,1803,65162,127785,217673,204328,542095,117185,12,13,5,21,11,
		67,101,161,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
		26401,1803,65162,127785,217673,204328,542095,117185,12,11,5,21,66,
		67,101,14,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
		26401,1803,65162,127785,217673,121626,542095,117185,7,11,5,21,66,
		67,101,14,90,101,195,563,1073,105,3148,1097,8204,4377,10145,12138,
		26401,1803,65162,127785,217673,121626,542095,117185,3,10,11,21,66,
		67,101,14,90,101,195,563,1073,105,3148,1097,8204,4377,10145,12138,
		12982,1803,65162,127785,217673,121626,542095,117185,3,10,13,21,66,
		67,101,14,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
		40398,1803,65162,127785,217673,121626,542095,60731,3,15,11,21,66,
		67,101,14,90,101,195,563,1073,105,3148,754,8204,4377,10145,12138,
		40398,1803,65162,127785,210249,121626,542095,60731,7,15,11,21,66,
		67,101,14,243,101,132,563,1073,105,3148,754,8204,4377,10145,12138,
		40398,1803,65162,80822,210249,200187,542095,60731,7,15,11,21,66,
		67,101,14,243,101,132,563,1073,105,3148,754,8204,4377,10145,12138,
		40398,1803,47650,80822,210249,200187,542095,60731,7,15,11,21,66,
		67,101,14,243,101,132,226,1073,105,1776,248,8204,4377,10145,12138,
		40398,1803,47650,80822,210249,200187,542095,60731,3,15,11,21,66,
		67,101,14,243,122,132,226,22,105,1776,754,8204,4377,10145,12138,
		40398,1803,47650,80822,210249,200187,542095,60731,3,15,11,21,45,
		67,101,14,243,122,132,226,22,105,1776,1097,8204,4377,10145,12138,
		3518,51108,47650,80822,210249,200187,542095,60731,3,15,11,21,11,
		67,101,14,243,122,132,226,22,105,3354,1097,8204,4377,10145,12138,
		3518,51108,47650,80822,210249,121551,542095,60731,3,15,13,21,7,67,
		101,14,243,122,132,226,22,105,3354,1097,8204,4377,10145,12138,
		3518,51108,47650,131661,210249,121551,542095,60731,3,6,13,21,3,67,
		101,14,243,122,132,226,22,105,3354,1097,8204,4377,10145,12138,
		37799,51108,47650,131661,210249,248492,542095,60731,3,2,11,21,2,
		67,101,14,243,122,132,226,22,105,925,222,8204,4377,10145,9305,
		37799,51108,40586,131661,210249,248492,542095,60731,3,3,13,17,2,
		51,101,14,243,122,132,226,1073,105,3354,222,8204,4377,10145,11107,
		37799,51108,40586,131661,94453,248492,277743,178309,3,2,5,17,2,51,
		101,14,283,122,132,226,452,105,3354,222,8204,4377,10145,11107,
		37799,51108,40586,131661,94453,248492,277743,178309,3,3,5,17,27,
		51,38,14,283,122,387,226,452,784,925,222,8204,4377,10145,11107,
		37799,51108,40586,131661,94453,248492,277743,178309,3,2,5,6,5,51,
		38,10,283,122,387,226,452,784,925,754,8204,4377,10145,11107,37799,
		51108,40586,131661,94453,248492,457259,178309,3,2,5,17,3,51,38,10,
		283,122,387,226,452,784,925,1982,4688,4377,10145,11107,37799,
		51108,40586,131661,94453,248492,457259,74373,3,2,14,17,3,12,38,10,
		283,122,387,226,452,784,925,1982,4688,4377,4544,11107,37799,51108,
		40586,131661,94453,248492,457259,74373,3,2,13,6,5,51,38,10,283,
		122,387,226,452,784,925,1982,4688,4377,4544,11107,37799,51108,
		38725,131661,94453,248492,457259,74373,3,2,5,3,5,12,38,10,283,122,
		387,226,318,784,2133,1982,2831,4377,4544,11107,4721,55315,38725,
		131661,94453,248492,457259,74373,3,2,5,6,2,51,38,10,283,122,387,
		226,301,784,2133,1982,2831,4377,4544,11107,4721,55315,38725,
		131661,94453,248492,457259,74373,3,2,5,6,2,5,38,103,283,122,387,
		226,301,784,2133,1982,2831,4377,4544,11107,4721,54140,38725,
		131661,94453,248492,457259,74373,3,2,5,3,2,3,3,10,16,122,387,226,
		301,784,2133,1982,2831,440,4544,11107,4721,54140,88329,131661,
		94453,13942,457259,74373,3,2,5,3,2,3,3,10,283,101,387,226,301,784,
		2133,1982,2831,440,8394,11107,7067,54140,88329,131661,94453,13942,
		457259,74373,3,2,5,3,2,2,3,10,16,101,387,226,86,784,2133,1982,
		2831,1199,8394,11107,7067,54140,88329,131661,94453,13942,457259,
		214965,3,2,5,3,2,2,3,10,283,101,387,226,86,784,2133,1982,2831,
		1199,8394,9305,7067,54140,88329,7114,94453,13942,457259,214965,3,
		2,5,3,2,5,3,5,283,101,387,226,15,784,2133,1982,2831,1199,8394,
		9305,7067,13134,88329,131661,94453,13942,457259,214965 };

	/* System generated locals */
	int i__1, i__2;
	double d__1, d__2, d__3;

	/* Local variables */
	static int i__, k;
	static double r__[1000], x[1000], fs[5000];
	static int np;
	static double vk[1000];
	static int pr[1000], kmx;
	static double difint, finval[5000], varprd;
	static int sampls;
	static double values[5000], varest[5000], varsqr[5000];
	static int intvls;


	/*  Automatic Multidimensional Integration Subroutine */

	/*         AUTHOR: Alan Genz */
	/*                 Department of Mathematics */
	/*                 Washington State University */
	/*                 Pulman, WA 99164-3113 */
	/*                 Email: AlanGenz@wsu.edu */

	/*         Last Change: 12/15/00 */

	/*  MVKBRV computes an approximation to the integral */

	/*      1  1     1 */
	/*     I  I ... I       F(X)  dx(NDIM)...dx(2)dx(1) */
	/*      0  0     0 */

	/*    F(X) is a real NF-vector of integrands. */

	/*  It uses randomized Korobov rules. The primary references are */
	/*   "Randomization of Number Theoretic Methods for Multiple Integration" */
	/*    R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13, pp. 904-14, */
	/*  and */
	/*   "Optimal Parameters for Multidimensional Integration", */
	/*    P. Keast, SIAM J Numer Anal, 10, pp.831-838. */
	/*  If there are more than 100 variables, the remaining variables are */
	/*  integrated using the rules described in the reference */
	/*   "On a Number-Theoretical Integration Method" */
	/*   H. Niederreiter, Aequationes Mathematicae, 8(1972), pp. 304-11. */

	/* **************  Parameters ******************************************** */
	/* ***** Input parameters */
	/*  NDIM    Number of variables, must exceed 1, but not exceed 100 */
	/*  MINVLS  Integer minimum number of function evaluations allowed. */
	/*          MINVLS must not exceed MAXVLS.  If MINVLS < 0 then the */
	/*          routine assumes a previous call has been made with */
	/*          the same integrands and continues that calculation. */
	/*  MAXVLS  Integer maximum number of function evaluations allowed. */
	/*  NF      Number of integrands, must exceed 1, but not exceed 5000 */
	/*  FUNSUB  EXTERNALly declared user defined integrand subroutine. */
	/*          It must have parameters ( NDIM, Z, NF, FUNVLS ), where */
	/*          Z is a real NDIM-vector and FUNVLS is a real NF-vector. */

	/*  ABSEPS  Required absolute accuracy. */
	/*  RELEPS  Required relative accuracy. */
	/* ***** Output parameters */
	/*  MINVLS  Actual number of function evaluations used. */
	/*  ABSERR  Maximum norm of estimated absolute accuracy of FINEST. */
	/*  FINEST  Estimated NF-vector of values of the integrals. */
	/*  INFORM  INFORM = 0 for normal exit, when */
	/*                     ABSERR <= MAX(ABSEPS, RELEPS*||FINEST||) */
	/*                  and */
	/*                     INTVLS <= MAXCLS. */
	/*          INFORM = 1 If MAXVLS was too small to obtain the required */
	/*          accuracy. In this case a value FINEST is returned with */
	/*          estimated absolute accuracy ABSERR. */
	/* *********************************************************************** */
	/* Parameter adjustments */
	--finest;

	/* Function Body */
	*inform__ = 1;
	intvls = 0;
	varprd = 0.;
	if (*minvls >= 0) {
		i__1 = *nf;
		for (k = 1; k <= i__1; ++k) {
			finest[k] = 0.;
			varest[k - 1] = 0.;
		}
		sampls = 8;
		for (i__ = min(*ndim,10); i__ <= 28; ++i__) {
			np = i__;
			if (*minvls < (sampls << 1) * p[i__ - 1]) {
				goto L10;
			}
		}
		/* Computing MAX */
		i__1 = 8, i__2 = *minvls / (p[np - 1] << 1);
		sampls = max(i__1,i__2);
	}
L10:
	vk[0] = 1. / p[np - 1];
	k = 1;
	i__1 = *ndim;
	for (i__ = 2; i__ <= i__1; ++i__) {
		if (i__ <= 100) {
			/* Computing MIN */
			i__2 = *ndim - 1;
			d__1 = c__[np + min(i__2,99) * 28 - 29] * (double) k;
			d__2 = (double) p[np - 1];
			k = (int) d_mod(&d__1, &d__2);
			vk[i__ - 1] = k * vk[0];
		} else {
			d__1 = (double) (i__ - 100) / (*ndim - 99);
			vk[i__ - 1] = (double) ((int) (p[np - 1] * pow(c_b32, d__1)));
			d__1 = vk[i__ - 1] / p[np - 1];
			vk[i__ - 1] = d_mod(&d__1, &c_b24);
		}
	}
	i__1 = *nf;
	for (k = 1; k <= i__1; ++k) {
		finval[k - 1] = 0.;
		varsqr[k - 1] = 0.;
	}

	i__1 = sampls;
	for (i__ = 1; i__ <= i__1; ++i__) {
		mvkrsv_(ndim, &c__100, values, &p[np - 1], vk, nf, (D_fp)funsub, x, 
			r__, pr, fs);
		i__2 = *nf;
		for (k = 1; k <= i__2; ++k) {
			difint = (values[k - 1] - finval[k - 1]) / i__;
			finval[k - 1] += difint;
			/* Computing 2nd power */
			d__1 = difint;
			varsqr[k - 1] = (i__ - 2) * varsqr[k - 1] / i__ + d__1 * d__1;
		}
	}

	intvls += (sampls << 1) * p[np - 1];
	kmx = 1;
	i__1 = *nf;
	for (k = 1; k <= i__1; ++k) {
		varprd = varest[k - 1] * varsqr[k - 1];
		finest[k] += (finval[k - 1] - finest[k]) / (varprd + 1);
		if (varsqr[k - 1] > 0.) {
			varest[k - 1] = (varprd + 1) / varsqr[k - 1];
		}
		if ((d__1 = finest[k], fabs(d__1)) > (d__2 = finest[kmx], fabs(d__2))) {
			kmx = k;
		}
	}
	*abserr = sqrt(varsqr[kmx - 1] / (varprd + 1)) * 7 / 2;
	/* Computing MAX */
	d__2 = *abseps, d__3 = (d__1 = finest[kmx], fabs(d__1)) * *releps;
	if (*abserr > max(d__2,d__3)) {
		if (np < 28) {
			++np;
		} else {
			/* Computing MIN */
			i__1 = sampls * 3 / 2, i__2 = (*maxvls - intvls) / (p[np - 1] << 
				1);
			sampls = min(i__1,i__2);
			sampls = max(8,sampls);
		}
		if (intvls + (sampls << 1) * p[np - 1] <= *maxvls) {
			goto L10;
		}
	} else {
		*inform__ = 0;
	}
	*minvls = intvls;

	/*    Optimal Parameters for Lattice Rules */


	return 0;
} /* mvkbrv_ */


/* Subroutine */ int mvkrsv_(int *ndim, int *kl, double *values, 
							 int *prime, double *vk, int *nf, D_fp funsub, double *
							 x, double *r__, int *pr, double *fs)
{
	/* System generated locals */
	int i__1, i__2;
	double d__1;

	/* Local variables */
	static int j, k, jp;


	/*     For lattice rule sums */

	/* Parameter adjustments */
	--fs;
	--pr;
	--r__;
	--x;
	--vk;
	--values;

	/* Function Body */
	i__1 = *nf;
	for (j = 1; j <= i__1; ++j) {
		values[j] = 0.;
	}

	/*     Determine random shifts for each variable; scramble lattice rule */

	i__1 = *ndim;
	for (j = 1; j <= i__1; ++j) {
		r__[j] = wsUnifrand();
		if (j < *kl) {
			jp = (int) (j * r__[j] + 1);
			if (jp < j) {
				pr[j] = pr[jp];
			}
			pr[jp] = j;
		} else {
			pr[j] = j;
		}
	}

	/*     Compute latice rule sums */

	i__1 = *prime;
	for (k = 1; k <= i__1; ++k) {
		i__2 = *ndim;
		for (j = 1; j <= i__2; ++j) {
			r__[j] += vk[pr[j]];
			if (r__[j] > 1.) {
				--r__[j];
			}
			x[j] = (d__1 = r__[j] * 2 - 1, fabs(d__1));
		}
		(*funsub)(ndim, &x[1], nf, &fs[1]);
		i__2 = *nf;
		for (j = 1; j <= i__2; ++j) {
			values[j] += (fs[j] - values[j]) / ((k << 1) - 1);
		}
		i__2 = *ndim;
		for (j = 1; j <= i__2; ++j) {
			x[j] = 1 - x[j];
		}
		(*funsub)(ndim, &x[1], nf, &fs[1]);
		i__2 = *nf;
		for (j = 1; j <= i__2; ++j) {
			values[j] += (fs[j] - values[j]) / (k << 1);
		}
	}

	return 0;
} /* mvkrsv_ */

#define LOW 0.02425
#define HIGH 0.97575

/* Coefficients in rational approximations. */

static const double a[] =
{
	-3.969683028665376e+01,
	2.209460984245205e+02,
	-2.759285104469687e+02,
	1.383577518672690e+02,
	-3.066479806614716e+01,
	2.506628277459239e+00
};

static const double b[] =
{
	-5.447609879822406e+01,
	1.615858368580409e+02,
	-1.556989798598866e+02,
	6.680131188771972e+01,
	-1.328068155288572e+01
};

static const double c[] =
{
	-7.784894002430293e-03,
	-3.223964580411365e-01,
	-2.400758277161838e+00,
	-2.549732539343734e+00,
	4.374664141464968e+00,
	2.938163982698783e+00
};

static const double d[] =
{
	7.784695709041462e-03,
	3.224671290700398e-01,
	2.445134137142996e+00,
	3.754408661907416e+00
};

double ltqnorm(double p)
{
	double q, r;

	errno = 0;

	if (p < 0 || p > 1) {
		return 0.0;
	} else if (p == 0) {
		return -HUGE_VAL /* minus "infinity" */;
	} else if (p == 1) {
		return HUGE_VAL /* "infinity" */;
	} else if (p < LOW) {
		/* Rational approximation for lower region */
		q = sqrt(-2*log(p));
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	} else if (p > HIGH) {
		/* Rational approximation for upper region */
		q  = sqrt(-2*log(1-p));
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
			((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
	} else {
		/* Rational approximation for central region */
		q = p - 0.5;
		r = q*q;
		return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
			(((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
	}
}

wsReal satterthwaite(wsVecCst Ra_a, wsUint N_len, wsUint N_df, wsReal* Rp_df)
{
	wsUint	L_all	= N_len * N_df;
	wsVec	Ra_res	= sseVector(L_all);

	LOOP (i, N_df)
		memcpy(Ra_res + i*N_len, Ra_a, sizeof(wsReal)*N_len);

	wsReal	R_tr2	= W0;
	wsReal	R_tr	= sseVsum(Ra_res, L_all, &R_tr2);

	R_tr	/= (wsReal)L_all;
	R_tr2	/= (wsReal)L_all * SQR(R_tr);

	*Rp_df	= L_all / R_tr2;
	return R_tr * R_tr2;
}

wsReal sseVmax(wsVecCst Ra_v, wsUint N_sz)
{
	if (N_sz == 0) return WISARD_NAN;
	wsReal R_ret = Ra_v[0];
	for (wsUint i=1 ; i<N_sz ; i++)
		if (R_ret < Ra_v[i]) R_ret = Ra_v[i];
	return R_ret;
}

inline wsReal k0(wsVec Ra_zeta, wsVec Ra_lambda, wsUint N_sz)
{
	wsReal	R_ret	= W0;
	wsUint	N_med	= 0;
// #ifdef USE_SSE
// 	N_med = getMed(N_sz);
// 	sse_t sse_1 = sseSet(W1);
// 	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
// 		sse_t* sse_z = (sse_t *)(Ra_zeta + i);
// 		sse_t* sse_l = (sse_t *)(Ra_lambda + i);
// 		sse_t sse_T = sseSet(W2);
// 		sse_T = sseMul(sse_T, sseMul(*sse_z, *sse_l));
// 		sse_T = sseSub(sse_1, sse_T);
// 		
// 		
// 	}
// #endif
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_ret -= log(W1 - W2*Ra_zeta[i]*Ra_lambda[i]);
	return R_ret / W2;
}

//kprime0<-function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
inline wsReal k0(wsReal R_zeta, wsVec Ra_lambda, wsUint N_sz)
{
	wsReal	R_ret	= W0;
	wsUint	N_med	= 0;
// #ifdef USE_SSE
// 	N_med = getMed(N_sz);
// 	sse_t sse_1 = sseSet(W1);
// 	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
// 		sse_t* sse_z = (sse_t *)(Ra_zeta + i);
// 		sse_t* sse_l = (sse_t *)(Ra_lambda + i);
// 		sse_t sse_T = sseSet(W2);
// 		sse_T = sseMul(sse_T, sseMul(*sse_z, *sse_l));
// 		sse_T = sseSub(sse_1, sse_T);
// 		
// 		
// 	}
// #endif
	wsReal R_z2 = W2 * R_zeta;
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_ret -= log(W1 - R_z2*Ra_lambda[i]);
	return R_ret / W2;
}

inline wsReal kprime0(wsVec Ra_zeta, wsVec Ra_lambda, wsUint N_sz)
{
	wsReal	R_ret	= W0;
	wsUint	N_med	= 0;
	// #ifdef USE_SSE
	// 	N_med = getMed(N_sz);
	// 	sse_t sse_1 = sseSet(W1);
	// 	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
	// 		sse_t* sse_z = (sse_t *)(Ra_zeta + i);
	// 		sse_t* sse_l = (sse_t *)(Ra_lambda + i);
	// 		sse_t sse_T = sseSet(W2);
	// 		sse_T = sseMul(sse_T, sseMul(*sse_z, *sse_l));
	// 		sse_T = sseSub(sse_1, sse_T);
	// 		
	// 		
	// 	}
	// #endif
	for (wsUint i=N_med ; i<N_sz ; i++)
		R_ret -= log(W1 - W2*Ra_zeta[i]*Ra_lambda[i]);
	return R_ret / W2;
} 

// kprime0<-function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
inline wsReal kprime0(wsReal R_zeta, wsVec Ra_lambda, wsUint N_sz)
{
	wsReal	R_ret	= W0;
	wsUint	N_med	= 0;
	for (wsUint i=N_med ; i<N_sz ; i++) {
		wsReal R_l2 = Ra_lambda[i];
		wsReal R_denom = W1 - W2*R_zeta*Ra_lambda[i];

		R_ret += R_l2 / R_denom;
	}
	return R_ret;
}

inline wsReal kpprime0(wsVec Ra_zeta, wsVec Ra_lambda, wsUint N_sz)
{
	wsReal	R_ret	= W0;
	wsUint	N_med	= 0;
	for (wsUint i=N_med ; i<N_sz ; i++) {
		wsReal R_l2 = SQR(Ra_lambda[i]);
		wsReal R_denom = W1 - W2*Ra_zeta[i]*Ra_lambda[i];

		R_ret += R_l2 / SQR(R_denom);
	}
	return R_ret * W2;
} 

inline wsReal kpprime0(wsReal R_zeta, wsVec Ra_lambda, wsUint N_sz)
{
	wsReal	R_ret	= W0;
	wsUint	N_med	= 0;
	for (wsUint i=N_med ; i<N_sz ; i++) {
		wsReal R_l2 = SQR(Ra_lambda[i]);
		wsReal R_denom = W1 - W2*R_zeta*Ra_lambda[i];

		R_ret += R_l2 / SQR(R_denom);
	}
	return R_ret * W2;
}

double zeroin2(			/* An estimate of the root */
    double ax,				/* Left border | of the range	*/
    double bx,				/* Right border| the root is seeked*/
    double fa, double fb,		/* f(a), f(b) */
    double (*f)(double x, void *info),	/* Function under investigation	*/
    void *info,				/* Add'l info passed on to f	*/
    double *Tol,			/* Acceptable tolerance		*/
    int *Maxit)				/* Max # of iterations */
{
    double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
    double tol;
    int maxit;

    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = * Tol;

    /* First test if we have found a root at an endpoint */
    if(fa == 0.0) {
		*Tol = 0.0;
		*Maxit = 0;
		return a;
    }
    if(fb ==  0.0) {
		*Tol = 0.0;
		*Maxit = 0;
		return b;
    }

    while(maxit--) {		/* Main iteration loop	*/
		double prev_step = b-a;		/* Distance from the last but one
									   to the last approximation	*/
		double tol_act;			/* Actual tolerance		*/
		double p;			/* Interpolation step is calcu- */
		double q;			/* lated in the form p/q; divi-
						 * sion operations is delayed
						 * until the last moment	*/
		double new_step;		/* Step at this iteration	*/

		if (fabs(fc) < fabs(fb)) {
			/* Swap data for b to be the	*/
			a = b;  b = c;  c = a;	/* best approximation		*/
			fa=fb;  fb=fc;  fc=fa;
		}
		tol_act = 2*numeric_limits<double>::epsilon()*fabs(b) + tol/2;
		new_step = (c-b)/2;

		if (fabs(new_step) <= tol_act || fb == (double)0) {
			*Maxit -= maxit;
			*Tol = fabs(c-b);
			return b;			/* Acceptable approx. is found	*/
		}

		/* Decide if the interpolation can be tried	*/
		if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
			&& fabs(fa) > fabs(fb) ) {	/* and was in true direction,
						 * Interpolation may be tried	*/
			register double t1,cb,t2;
			cb = c-b;
			if( a==c ) {		/* If we have only two distinct	*/
						/* points linear interpolation	*/
			t1 = fb/fa;		/* can only be applied		*/
			p = cb*t1;
			q = 1.0 - t1;
			}
			else {			/* Quadric inverse interpolation*/

			q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
			p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
			q = (q-1.0) * (t1-1.0) * (t2-1.0);
			}
			if( p>(double)0 )		/* p was calculated with the */
			q = -q;			/* opposite sign; make p positive */
			else			/* and assign possible minus to	*/
			p = -p;			/* q				*/

			if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
				&& p < fabs(prev_step*q/2) )	/* and isn't too large	*/
			new_step = p/q;			/* it is accepted
							 * If p/q is too large then the
							 * bisection procedure can
							 * reduce [b,c] range to more
							 * extent */
		}

		if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
			if( new_step > (double)0 )	/* than tolerance		*/
				new_step = tol_act;
			else
				new_step = -tol_act;
		}
		a = b;	fa = fb;			/* Save the previous approx. */
		b += new_step;	fb = (*f)(b, info);	/* Do step to a new approxim. */
		if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
			/* Adjust c for it to have a sign opposite to that of b */
			c = a;  fc = fa;
		}
    }

    /* failed! */
    *Tol = fabs(c-b);
    *Maxit = -1;
    return b;
}

double uniroot(double (*H_func)(double x, void *info), wsReal R_lower,
	wsReal R_upper, void *Vp_info, double R_tol=1e-8)
{
	if (NA(R_lower) || NA(R_upper) || R_lower >= R_upper) 
		halt("lower < upper  is not fulfilled");

	wsReal R_fLower = H_func(R_lower, Vp_info);
	wsReal R_fUpper = H_func(R_upper, Vp_info);

    if (NA(R_fLower)) 
        halt("f.lower = f(lower) is NA");
    if (NA(R_fUpper)) 
        halt("f.upper = f(upper) is NA");

    if (R_fLower * R_fUpper > 0) 
        halt("f() values at end points not of opposite sign");
	int N_maxIt = 1000;
	double R_ret = zeroin2(R_lower, R_upper, R_fLower, R_fUpper, H_func,
		Vp_info, &R_tol, &N_maxIt);

	return R_ret;
}

typedef struct _xZet {
	wsVec Ra_lambda;
	wsUint N_sz;
} xZet;

double zet(double x, void *info)
{
	xZet* Xp_i = (xZet *)info;
	return kprime0(x, Xp_i->Ra_lambda, Xp_i->N_sz) - x;
}

//saddle<-function(x,lambda){
wsReal saddlepoint(wsReal R_x, wsVec Ra_lambda, wsUint N_sz)
{
	// d<-max(lambda)
	wsReal d = sseVmax(Ra_lambda, N_sz);
	// lambda<-lambda/d
	sseVpC(Ra_lambda, W1/d, Ra_lambda, N_sz);
	// x<-x/d
	R_x /= d;
	//	k0<-function(zeta) -sum(log(1-2*zeta*lambda))/2
	//	kprime0<-function(zeta) sapply(zeta, function(zz) sum(lambda/(1-2*zz*lambda)))
	//	kpprime0<-function(zeta) 2*sum(lambda^2/(1-2*zeta*lambda)^2)
	// n<-length(lambda)

	wsReal	R_lmin	= WISARD_NAN;
	wsReal	R_v		= W0;
	char	B_negL	= 0;
	LOOP (i, N_sz) {
		if (Ra_lambda[i] < W0) B_negL = 1;
		R_v += Ra_lambda[i];
	}

	// if (any(lambda < 0)) {
	if (B_negL) {
		wsReal R_curLmax = -numeric_limits<wsReal>::infinity();
		LOOP (i, N_sz) if (Ra_lambda[i] < W0) {
			wsReal R_curV = W1 / (W2 * Ra_lambda[i]);
			if (R_curLmax < R_curV) R_curLmax = R_curV;
		}
		R_lmin = R_curLmax * REAL_CONST(.99999);
		// lmin <- max(1/(2 * lambda[lambda < 0])) * 0.99999
	// } else if (x>sum(lambda)){
	} else if (R_x > R_v) {
		R_lmin = -0.01;
		// lmin <- -0.01
	} else {
		R_lmin = -(N_sz / (W2 * R_x));
		// lmin<- -length(lambda)/(2*x)
	}

	wsReal R_curLmin = numeric_limits<wsReal>::infinity();
	LOOP (i, N_sz) if (Ra_lambda[i] > W0) {
		wsReal R_curV = W1 / (W2 * Ra_lambda[i]);
		if (R_curLmin > R_curV) R_curLmin = R_curV;
	}
	wsReal R_lmax = R_curLmin * REAL_CONST(.99999);
	// lmax<-min(1/(2*lambda[lambda>0]))*0.99999

	// hatzeta <- uniroot(function(zeta) kprime0(zeta) - x, 
	//	lower = lmin, upper = lmax, tol = 1e-08)$root
	xZet X = { Ra_lambda, N_sz };
	wsReal R_hz = uniroot(zet, R_lmin, R_lmax, &X, 1e-8);

	//	k0<-function(zeta) -sum(log(1-2*zeta*lambda))/2
#define sign(a) (R_hz >= W0 ? 1 : -1)
	wsReal R_w = sign(R_hz)*sqrt(W2*(R_hz*R_x-k0(R_hz, Ra_lambda, N_sz)));
	wsReal R_vv = R_hz * sqrt(kpprime0(R_hz, Ra_lambda, N_sz));
//	if (abs(hatzeta)<1e-4)
//		NA
//	else
//		pnorm(w+log(v/w)/w, lower.tail=FALSE)
	if (0)
		return WISARD_NAN;
	else
		return pnorm(R_w+log(R_vv/R_w)/R_w, W0, W1, 0, 0);
}

wsRealCst pchisqsum(double x, double df, wsVec a, wsUint N_sz, char lower_tail/*=TRUE*/,
	char method/*=0*/)
/*	0 = "satterthwaite"
	1 = "integration"
	2 = "saddlepoint"*/
{
	wsReal R_ret	= WISARD_NAN;
	wsReal R_satDF = WISARD_NAN;
	// sat<-satterthwaite(a,df)
	wsReal R_satScale = satterthwaite(a, N_sz, (wsUint)df, &R_satDF);
	// guess<-pchisq(x/sat$scale,sat$df,lower.tail=lower.tail)
	wsReal R_guess = lower_tail ? W1 - PVchisq(x/R_satScale, R_satDF)
		: PVchisq(x/R_satScale, R_satDF);

	// if (method=="satterthwaite") return(guess)

	wsReal R_abstol = R_guess / REAL_CONST(1000.0);
	R_abstol = max(1e-9, R_abstol);
//	wsReal R_reltol = REAL_CONST(.001);

	if (method == 1) {
		char B_useDavies = 0;
		LOOP (i, N_sz) if (a[i] <= W0) {
			B_useDavies = 1;
			break;
		}

		wsReal R_itGuess = WISARD_NAN;
		if (B_useDavies) {
			R_itGuess = davies(x, a, N_sz);
//						if (f$ifault>0) warning("Probable loss of accuracy ")
//							guess[i]<-f$Qq
//				}
//				if(any(guess<1e-6)) warning("Probable loss of accuracy ")
		} else {
			R_itGuess = farebrother(x, a, N_sz);
//			if(any(guess<1e-9)) warning("Probable loss of accuracy ")
		}
		if (lower_tail)
			R_itGuess = W1 - R_itGuess;

		/* Set the return p-value */
		if (NA(R_itGuess)) R_ret = R_guess;
		else R_ret = R_itGuess;
	} else if (method == 2) {
		wsUint	N_df		= (wsUint)df;
		wsVec	Ra_lambda	= sseVector(N_df*N_sz);
		LOOP (i, N_df)
			memcpy(Ra_lambda+(i*N_sz), a, sizeof(wsReal)*N_sz);
		wsReal	R_spGuess	= saddlepoint(x, Ra_lambda, N_df*N_sz);

		if (lower_tail)
			R_spGuess = W1 - R_spGuess;

		/* Set the return p-value */
		if (NA(R_spGuess)) R_ret = R_guess;
		else R_ret = R_spGuess;
	} else
		R_ret = R_guess;

	return R_ret;
}

#define white_space(c) ((c) == ' ' || (c) == '\t')
#define valid_digit(c) ((c) >= '0' && (c) <= '9')

double ws_str2dbl(const char *S_str, char **Sp_end)
{
	int frac;
	double R_ssgn, R_value, R_scale;

	// Skip leading white space, if any.

	while (white_space(*S_str)) S_str++;

	// Get sign, if any.

	R_ssgn = 1.0;
	if (*S_str == '-') {
		R_ssgn = -1.0;
		S_str += 1;
	} else if (*S_str == '+')
		S_str += 1;

	// Get digits before decimal point or exponent, if any.

	for (R_value = 0.0; valid_digit(*S_str); S_str += 1) {
		R_value = R_value * 10.0 + (*S_str - '0');
	}

	// Get digits after decimal point, if any.

	if (*S_str == '.') {
		double pow10 = 10.0;
		S_str += 1;
		while (valid_digit(*S_str)) {
			R_value += (*S_str - '0') / pow10;
			pow10 *= 10.0;
			S_str += 1;
		}
	}

	// Handle exponent, if any.

	frac = 0;
	R_scale = 1.0;
	if ((*S_str == 'e') || (*S_str == 'E')) {
		unsigned int expon;

		// Get sign of exponent, if any.

		S_str += 1;
		if (*S_str == '-') {
			frac = 1;
			S_str += 1;

		} else if (*S_str == '+') {
			S_str += 1;
		}

		// Get digits of exponent, if any.

		for (expon = 0; valid_digit(*S_str); S_str += 1) {
			expon = expon * 10 + (*S_str - '0');
		}
		if (expon > 308) expon = 308;

		// Calculate scaling factor.

		while (expon >= 50) { R_scale *= 1E50; expon -= 50; }
		while (expon >=  8) { R_scale *= 1E8;  expon -=  8; }
		while (expon >   0) { R_scale *= 10.0; expon -=  1; }
	}
	// 
	if (Sp_end) *Sp_end = (char *)S_str;

	// Return signed and scaled floating point result.
	return R_ssgn * (frac ? (R_value / R_scale) : (R_value * R_scale));
}

wsReal randInRange(wsReal R_min, wsReal R_max)
{
	return R_min + ((wsReal)rand() / (wsReal)RAND_MAX ) * (R_max - R_min);
}


wsVec sseVrank(wsVecCst Ra_v, wsUint N_sz, xRankTie X_tie/*=RT_AVERAGE*/)
{
#define  MAX_LEVELS  300
	wsVec Ra_ret = sseVector(N_sz);
	xRealSort* Xa_sort = buildRealSort(Ra_v, N_sz);
	wsReal R_rank = W0;
	wsUint N_tie = 0;
	{
		int  beg[MAX_LEVELS], end[MAX_LEVELS], i = 0, L, R, swap;

		/* Build realsort */
		xRealSort piv;

		beg[0] = 0; end[0] = N_sz;
		while (i >= 0) {
			L = beg[i]; R = end[i] - 1;

			if (L < R) {
				piv.i = Xa_sort[L].i;
				piv.V = Xa_sort[L].V;

				while (L < R) {
					while (Xa_sort[R].V >= piv.V && L < R) R--;
					if (L < R) {
						Xa_sort[L].i = Xa_sort[R].i;
						Xa_sort[L].V = Xa_sort[R].V;
						L++;
					}
					while (Xa_sort[L].V <= piv.V && L < R) L++;
					if (L < R) {
						Xa_sort[R].i = Xa_sort[L].i;
						Xa_sort[R].V = Xa_sort[L].V;
						R--;
					}
				}
				Xa_sort[L].i = piv.i;
				Xa_sort[L].V = piv.V;
				beg[i + 1] = L + 1; end[i + 1] = end[i]; end[i++] = L;
				if (end[i] - beg[i]>end[i - 1] - beg[i - 1]) {
					swap = beg[i]; beg[i] = beg[i - 1]; beg[i - 1] = swap;
					swap = end[i]; end[i] = end[i - 1]; end[i - 1] = swap;
				}
			} else i--;
		}

		wsReal R_lastV = WISARD_NAN;
		for (wsUint i=0 ; i<N_sz ; i++) {
			wsUint I = Xa_sort[i].i;
			wsRealCst V = Xa_sort[i].V;

			/* Tie broken */
			if (R_lastV != V) {
				/* If there was more than one ties */
				if (N_tie > 1) {
					wsReal R_finRank = 0;
					switch (X_tie) {
					case RT_MIN: R_finRank = R_rank; break;
					case RT_MAX: R_finRank = R_rank + N_tie - 1; break;
					case RT_AVERAGE:
						R_finRank = ((R_rank + N_tie - W1)*(R_rank + N_tie) - R_rank*(R_rank - W1)) / W2; break;
					default:
						halt("Invalid rank tie method");
					}
					/* (i-1) ~ (i-N_tie), fill value */
					for (wsUint j=i - N_tie ; j<i ; j++)
						Ra_ret[Xa_sort[j].i] = R_finRank;
				}
			} else {
				/* TIe keep */
				N_tie++;
				continue;
			}

			/* Normal rank change */
			N_tie = 1;
			R_lastV = V;
			R_rank = (wsReal)i+1;
			Ra_ret[I] = R_rank;
		}
	}

	/* If there was more than one ties */
	if (N_tie > 1) {
		wsReal R_finRank = 0;
		switch (X_tie) {
		case RT_MIN: R_finRank = R_rank; break;
		case RT_MAX: R_finRank = R_rank + N_tie - 1; break;
		case RT_AVERAGE:
			R_finRank = ((R_rank + N_tie - W1)*(R_rank + N_tie) - R_rank*(R_rank - W1)) / W2; break;
		default:
			halt("Invalid rank tie method");
		}
		/* (i-1) ~ (i-N_tie), fill value */
		for (wsUint j=N_sz-N_tie ; j<N_sz ; j++)
			Ra_ret[Xa_sort[j].i] = R_finRank;
	}

	DEALLOC(Xa_sort);
	return Ra_ret;
}

wsVec sseVarank(wsVecCst Ra_v, wsUint N_sz, xRankTie X_tie/*=RT_AVERAGE*/)
{
#define  MAX_LEVELS  300
	wsVec Ra_ret = sseVector(N_sz);
	xRealSort* Xa_sort = buildArealSort(Ra_v, N_sz);
	wsReal R_rank = W0;
	wsUint N_tie = 0;
	{
		int  beg[MAX_LEVELS], end[MAX_LEVELS], i = 0, L, R, swap;

		/* Build realsort */
		xRealSort piv;

		beg[0] = 0; end[0] = N_sz;
		while (i >= 0) {
			L = beg[i]; R = end[i] - 1;

			if (L < R) {
				piv.i = Xa_sort[L].i;
				piv.V = Xa_sort[L].V;

				while (L < R) {
					while (Xa_sort[R].V >= piv.V && L < R) R--;
					if (L < R) {
						Xa_sort[L].i = Xa_sort[R].i;
						Xa_sort[L].V = Xa_sort[R].V;
						L++;
					}
					while (Xa_sort[L].V <= piv.V && L < R) L++;
					if (L < R) {
						Xa_sort[R].i = Xa_sort[L].i;
						Xa_sort[R].V = Xa_sort[L].V;
						R--;
					}
				}
				Xa_sort[L].i = piv.i;
				Xa_sort[L].V = piv.V;
				beg[i + 1] = L + 1; end[i + 1] = end[i]; end[i++] = L;
				if (end[i] - beg[i]>end[i - 1] - beg[i - 1]) {
					swap = beg[i]; beg[i] = beg[i - 1]; beg[i - 1] = swap;
					swap = end[i]; end[i] = end[i - 1]; end[i - 1] = swap;
				}
			} else i--;
		}

		wsReal R_lastV = WISARD_NAN;
		for (wsUint i=0 ; i<N_sz ; i++) {
			wsUint I = Xa_sort[i].i;
			wsRealCst V = Xa_sort[i].V;

			/* Tie broken */
			if (R_lastV != V) {
				/* If there was more than one ties */
				if (N_tie > 1) {
					wsReal R_finRank = 0;
					switch (X_tie) {
					case RT_MIN: R_finRank = R_rank; break;
					case RT_MAX: R_finRank = R_rank + N_tie - 1; break;
					case RT_AVERAGE:
						R_finRank = ((R_rank + N_tie - W1)*(R_rank + N_tie) - R_rank*(R_rank - W1)) / W2 / (wsReal)N_tie; break;
					default:
						halt("Invalid rank tie method");
					}
					if (R_finRank > N_sz)
						halt("Errornous rank [%g], should not over [%d]!\n", R_finRank, N_sz);
					/* (i-1) ~ (i-N_tie), fill value */
					for (wsUint j=i - N_tie ; j<i ; j++)
						Ra_ret[Xa_sort[j].i] = R_finRank;
				}
			} else {
				/* TIe keep */
				N_tie++;
				continue;
			}

			/* Normal rank change */
			N_tie = 1;
			R_lastV = V;
			R_rank = (wsReal)i+1;
			Ra_ret[I] = R_rank;
		}
	}

	/* If there was more than one ties */
	if (N_tie > 1) {
		wsReal R_finRank = 0;
		switch (X_tie) {
		case RT_MIN: R_finRank = R_rank; break;
		case RT_MAX: R_finRank = R_rank + N_tie - 1; break;
		case RT_AVERAGE:
			R_finRank = ((R_rank + N_tie - W1)*(R_rank + N_tie) - R_rank*(R_rank - W1)) / W2 / (wsReal)N_tie; break;
		default:
			halt("Invalid rank tie method");
		}
		if (R_finRank > N_sz)
			halt("Errornous rank [%g], should not over [%d]!\n", R_finRank, N_sz);
		/* (i-1) ~ (i-N_tie), fill value */
		for (wsUint j=N_sz-N_tie ; j<N_sz ; j++)
			Ra_ret[Xa_sort[j].i] = R_finRank;
	}

	DEALLOC(Xa_sort);
	return Ra_ret;
}

void getMachineName(char* Sp_mname)
{
	char S_name[150];

#ifdef WIN32
	int i = 0;
	TCHAR S_infoBuf[150];
	DWORD N_szBuf = 150;
	memset(S_name, 0, 150);
	if (GetComputerName(S_infoBuf, &N_szBuf)) {
		for (i=0; i < 150; i++)
			S_name[i] = S_infoBuf[i];
	} else
		strcpy(S_name, "Unknown_Host_Name");
#else
	memset(S_name, 0, 150);
	gethostname(S_name, 150);
#endif
	strncpy(Sp_mname, S_name, 150);
}

void Qbenjhoch(xRealSort* Xa_rs, wsUintCst N_sz, wsVec Ra_ret)
{
	wsVec Ra_tmp = sseVector(N_sz);

	// test.p<-length(pvals)/(1:length(pvals))*pvals
	LOOP (i, N_sz) Ra_tmp[i] = (N_sz/(wsReal)(i+W1)) * Xa_rs[i].V;
	// for { adj.p[i] <- min(test.p[i:length(test.p)]); ifelse(adj.p[i]>1, 1, adj.p[i])
	LOOP (i, N_sz) {
		wsReal R_min = Ra_tmp[i];
		for (wsUint j=i + 1; j<N_sz; j++) if (R_min > Ra_tmp[j]) R_min = Ra_tmp[j];
		Ra_ret[Xa_rs[i].i] = min(W1, R_min);
	}
	sseFree(Ra_tmp);
}

} // End namespace ONETOOL
