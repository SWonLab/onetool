#include <math.h>
#include <setjmp.h>
#include "utils/dcdflib.h"
#include "utils/cdflib.h"
#include "utils/util.h"
#include "utils/stat.h"
#include "utils/new_stat.h"

#ifdef __cplusplus
#define W_MAX(Ra_mat,b) (((Ra_mat) > (b)) ? (Ra_mat) : (b))
#define min(Ra_mat,b) (((Ra_mat) < (b)) ? (Ra_mat) : (b))
#endif

namespace ONETOOL {

/******* FET ********/

/* Global values for FET calculation */
//double		R_sProb;
//static int	N_s1X, N_s11, N_sX1, N_sXX;

/* Macros for EFT calculation */
#define FET_B(N_count,k)	(FET_F(N_count)-FET_F(k)-FET_F((N_count)-(k)))
#define FET_X(Ra_mat,b,c,d)		exp(FET_B(b,Ra_mat)+FET_B(d-b,c-Ra_mat)-FET_B(d,c))
#define FET_H(Ra_mat)			FET_HY(Ra_prob, N_s11, N_s1X, N_sX1, N_sXX, Ra_mat, 0, 0, 0)

double FET_G(double z)
{
	double x = 0;
	x += 0.1659470187408462e-06/(z+7);
	x += 0.9934937113930748e-05/(z+6);
	x -= 0.1385710331296526    /(z+5);
	x += 12.50734324009056     /(z+4);
	x -= 176.6150291498386     /(z+3);
	x += 771.3234287757674     /(z+2);
	x -= 1259.139216722289     /(z+1);
	x += 676.5203681218835     /(z);
	x += 0.9999999999995183;

	return(log(x)-5.58106146679532777-z+(z-0.5)*log(z+6.5));
}

inline double FET_F(double n)
{
	if (n<=1.0) return 0.0;
	return FET_G(n+1);
}

double FET_HY(double& R_sProb, int& N_s11, int& N_s1X, int& N_sX1, int& N_sXX, int N_11i, int N_1Xi, int N_X1i, int N_XXi)
{
	/* For hyper1 */
	if (N_1Xi|N_X1i|N_XXi) {
		N_s11 = N_11i;
		N_s1X = N_1Xi;
		N_sX1 = N_X1i;
		N_sXX = N_XXi;
	} else {
		if (!(N_11i % 10 == 0)) {
			if(N_11i == (N_s11+1)) {
				R_sProb *= ((N_s1X-N_s11)/(double)(N_11i))*((N_sX1-N_s11)/(double)(N_11i+N_sXX-N_s1X-N_sX1));
				N_s11 = N_11i;
				return R_sProb;
			}
			if(N_11i==N_s11-1) {
				R_sProb *= ((N_s11)/(double)(N_s1X-N_11i))*((N_s11+N_sXX-N_s1X-N_sX1)/(double)(N_sX1-N_11i));
				N_s11 = N_11i;
				return R_sProb;
			}
		}
		N_s11 = N_11i;
	}

	R_sProb = FET_X(N_s11, N_s1X, N_sX1, N_sXX);
	return R_sProb;
}

/* Fisher's Exact Test (two-side)
 *
 * -------------+----------------------------------------------------
 * int	N_11	| Count of (1,1) cell
 * int	N_12	| Count of (1,2) cell
 * int	N_21	| Count of (2,1) cell
 * int	N_22	| Count of (2,2) cell
 * -------------+----------------------------------------------------
 * int	*Na_2x2	| 4 int array represents 2x2 contingency table
 * -------------+----------------------------------------------------
 *  return		| p-value of FET
 */
wsRealCst	fisher(int *Na_2x2)
{
	return fisher(Na_2x2[0], Na_2x2[1], Na_2x2[2], Na_2x2[3]);
}
wsRealCst	fisher(int N_11, int N_12, int N_21, int N_22)
{
	int		i;
	/* For abbreviation */
	int		N_1X = N_11 + N_12;
	int		N_X1 = N_11 + N_21;
	int		N_XX = N_11 + N_12 +N_21 +N_22;
	/* Values for calculation */
	double	R_R = 0.0, R_L = 0.0;
	double	R_p, R_baseP;
	int		N_mx = N_1X >= N_X1 ? N_X1 : N_1X;
	int		N_mi = N_1X + N_X1 - N_XX;
	double	R_sProb = 0.0;
	int		N_s11 = 0;
	int		N_s1X = 0;
	int		N_sX1 = 0;
	int		N_sXX = 0;

	/* Invalid value */
	if (N_11<0 || N_12<0 || N_21<0 || N_22<0) {
		LOG("Invalid cell value for Fisher's exact test [%d,%d ; %d,%d]!\n", N_11, N_12, N_21, N_22);
		exit(1);
	}
	if (N_mi < 0) N_mi = 0;

	/* Tie case */
	if (N_mi == N_mx)
		return 1.0;

	/* Calculate hypergeometric prob. */
	R_baseP	= FET_HY(R_sProb, N_s11, N_s1X, N_sX1, N_sXX, N_11, N_1X, N_X1, N_XX);
	R_p		= FET_HY(R_sProb, N_s11, N_s1X, N_sX1, N_sXX, N_mi, 0, 0, 0);
	//printf("prob %g p %g\n", R_baseP, R_p);
	for(i=N_mi+1 ; R_p<0.99999999*R_baseP ; i++) {
		R_L += R_p;
		R_p	=  FET_HY(R_sProb, N_s11, N_s1X, N_sX1, N_sXX, i, 0, 0, 0);
	}
	if (R_p<1.00000001*R_baseP) R_L += R_p;

	R_p		= FET_HY(R_sProb, N_s11, N_s1X, N_sX1, N_sXX, N_mx, 0, 0, 0);
	for(i=N_mx-1 ; R_p<0.99999999*R_baseP ; i--) {
		R_R	+= R_p;
		R_p	=  FET_HY(R_sProb, N_s11, N_s1X, N_sX1, N_sXX, i, 0, 0, 0);
	}
	if (R_p<1.00000001*R_baseP) R_R += R_p;

	//printf("prob %g %g\n", R_L+R_R, (wsReal)(R_L+R_R));
	R_p = R_L+R_R;
	if(R_p > 1.0)
		R_p = 1.0;

	return (wsReal)R_p;
} /* END OF fisher() */

double fisher(wsReal R_11, wsReal R_12, wsReal R_13, wsReal R_21, wsReal R_22, wsReal R_23)
{
// 	int nrow = 2;
// 	int ncol = 3;
// 	double table[200];
// 
// 	int c=0;
// 	table[0] = R_11;
// 	table[1] = R_12;
// 	table[2] = R_13;
// 	table[3] = R_21;
// 	table[4] = R_22;
// 	table[5] = R_23;
// 	double expect = -1.0;
// 	double percnt = 100.0;
// 	double emin = 0;
 	double pre = 0;
//	double prt = 0;
// 	int ws = 300000;

//	fexact(&nrow, &ncol, table, &nrow, &expect, &percnt, &emin, &prt, &pre, &ws);

	return pre;
}

/*
 * Komogorov-Smirnov test
 */

//______________________________________________________________________________
int Nint(double x)
{
	// Round to nearest integer. Rounds half integers to the nearest
	// even integer.

	int i;
	if (x >= 0) {
		i = int(x + 0.5);
		if (x + 0.5 == double(i) && i & 1) i--;
	}
	else {
		i = int(x - 0.5);
		if (x - 0.5 == double(i) && i & 1) i++;

	}
	return i;
}

//______________________________________________________________________________
double PVkolmogorov(double z)
{
	// Calculates the Kolmogorov distribution function,
	// which gives the probability that Kolmogorov's test statistic will exceed
	// the value z assuming the null hypothesis. This gives a very powerful
	// test for comparing two one-dimensional distributions.
	// see, for example, Eadie et al, "statistocal Methods in Experimental
	// Physics', pp 269-270).
	//
	// This function returns the confidence level for the null hypothesis, where:
	//   z = dn*sqrt(n), and
	//   dn  is the maximum deviation between a hypothetical distribution
	//       function and an experimental distribution with
	//   n    events
	//
	// NOTE: To compare two experimental distributions with m and n events,
	//       use z = sqrt(m*n/(m+n))*dn
	//
	// Accuracy: The function is far too accurate for any imaginable application.
	//           Probabilities less than 10^-15 are returned as zero.
	//           However, remember that the formula is only valid for "large" n.
	// Theta function inversion formula is used for z <= 1
	//
	// This function was translated by Rene Brun from PROBKL in CERNLIB.

	double fj[4] = { -2, -8, -18, -32 }, r[4];
	const double w = 2.50662827;
	// c1 - -pi**2/8, c2 = 9*c1, c3 = 25*c1
	const double c1 = -1.2337005501361697;
	const double c2 = -11.103304951225528;
	const double c3 = -30.842513753404244;

	double u = fabs(z);
	double p;
	if (u < 0.2) {
		p = 1;
	}
	else if (u < 0.755) {
		double v = 1. / (u*u);
		p = 1 - w*(exp(c1*v) + exp(c2*v) + exp(c3*v)) / u;
	}
	else if (u < 6.8116) {
		r[1] = 0;
		r[2] = 0;
		r[3] = 0;
		double v = u*u;
		int maxj = W_MAX(1, Nint(3. / u));
		for (int j = 0; j < maxj; j++) {
			r[j] = exp(fj[j] * v);
		}
		p = 2 * (r[0] - r[1] + r[2] - r[3]);
	}
	else {
		p = 0;
	}
	return p;
}

//______________________________________________________________________________
wsReal kstest(int N_aCount, const wsReal *Ra_aData, int N_bCount, const wsReal *Ra_bData)
{
	//  Statistical test whether two one-dimensional sets of points are compatible
	//  with coming from the same parent distribution, using the Kolmogorov test.
	//  That is, it is used to compare two experimental distributions of unbinned data.
	//
	//  Input:
	//  a,b: One-dimensional arrays of length na, nb, respectively.
	//       The elements of a and b must be given in ascending order.
	//  option is a character string to specify options
	//         "D" Put out a line of "Debug" printout
	//         "M" Return the Maximum Kolmogorov distance instead of prob
	//
	//  Output:
	// The returned value prob is a calculated confidence level which gives a
	// statistical test for compatibility of a and b.
	// Values of prob close to zero are taken as indicating a small probability
	// of compatibility. For two point sets drawn randomly from the same parent
	// distribution, the value of prob should be uniformly distributed between
	// zero and one.
	//   in case of error the function return -1
	//   If the 2 sets have a different number of points, the minimum of
	//   the two sets is used.
	//
	// Method:
	// The Kolmogorov test is used. The test statistic is the maximum deviation
	// between the two integrated distribution functions, multiplied by the
	// normalizing factor (rdmax*sqrt(na*nb/(na+nb)).
	//
	//  Code adapted by Rene Brun from CERNLIB routine TKOLMO (Fred James)
	//   (W.T. Eadie, D. Drijard, F.E. James, M. Roos and B. Sadoulet,
	//      Statistical Methods in Experimental Physics, (North-Holland,
	//      Amsterdam 1971) 269-271)
	//
	//  Method Improvement by Jason A Detwiler (JADetwiler@lbl.gov)
	//  -----------------------------------------------------------
	//   The nuts-and-bolts of the TMath::KolmogorovTest() algorithm is a for-loop
	//   over the two sorted arrays a and b representing empirical distribution
	//   functions. The for-loop handles 3 cases: when the next points to be
	//   evaluated satisfy a>b, a<b, or a=b:
	//
	//      for (int i=0;i<na+nb;i++) {
	//         if (a[ia-1] < b[ib-1]) {
	//            rdiff -= sa;
	//            ia++;
	//            if (ia > na) {ok = kTRUE; break;}
	//         } else if (a[ia-1] > b[ib-1]) {
	//            rdiff += sb;
	//            ib++;
	//            if (ib > nb) {ok = kTRUE; break;}
	//         } else {
	//            rdiff += sb - sa;
	//            ia++;
	//            ib++;
	//            if (ia > na) {ok = kTRUE; break;}
	//            if (ib > nb) {ok = kTRUE; break;}
	//        }
	//         rdmax = TMath::Max(rdmax,TMath::Abs(rdiff));
	//      }
	//
	//   For the last case, a=b, the algorithm advances each array by one index in an
	//   attempt to move through the equality. However, this is incorrect when one or
	//   the other of a or b (or both) have a repeated value, call it x. For the KS
	//   statistic to be computed properly, rdiff needs to be calculated after all of
	//   the a and b at x have been tallied (this is due to the definition of the
	//   empirical distribution function; another way to convince yourself that the
	//   old CERNLIB method is wrong is that it implies that the function defined as the
	//   difference between a and b is multi-valued at x -- besides being ugly, this
	//   would invalidate Kolmogorov's theorem).
	//
	//   The solution is to just add while-loops into the equality-case handling to
	//   perform the tally:
	//
	//         } else {
	//            double x = a[ia-1];
	//            while(a[ia-1] == x && ia <= na) {
	//              rdiff -= sa;
	//              ia++;
	//            }
	//            while(b[ib-1] == x && ib <= nb) {
	//              rdiff += sb;
	//              ib++;
	//            }
	//            if (ia > na) {ok = kTRUE; break;}
	//            if (ib > nb) {ok = kTRUE; break;}
	//         }
	//
	//
	//  NOTE1
	//  A good description of the Kolmogorov test can be seen at:
	//    http://www.itl.nist.gov/div898/handbook/eda/section3/eda35g.htm

	//  LM: Nov 2010: clean up and returns now a zero distance when vectors are the same 

	double prob = -1;
	//      Require at least two points in each graph
	if (!Ra_aData || !Ra_bData || N_aCount <= 2 || N_bCount <= 2) {
		LOG("KolmogorovTest :: Sets must have more than 2 points");
		return (wsReal)prob;
	}
	//     Constants needed
	double rna = (double)N_aCount;
	double rnb = (double)N_bCount;
	double sa = 1. / rna;
	double sb = 1. / rnb;
	double rdiff = 0;
	double rdmax = 0;
	int ia = 0;
	int ib = 0;

	//    Main loop over point sets to find max distance
	//    rdiff is the running difference, and rdmax the max.
	bool ok = false;
	for (int i = 0; i < N_aCount + N_bCount; i++) {
		if (Ra_aData[ia] < Ra_bData[ib]) {
			rdiff -= sa;
			ia++;
			if (ia >= N_aCount) { ok = true; break; }
		}
		else if (Ra_aData[ia] > Ra_bData[ib]) {
			rdiff += sb;
			ib++;
			if (ib >= N_bCount) { ok = true; break; }
		}
		else {
			// special cases for the ties 
			double x = Ra_aData[ia];
			while (Ra_aData[ia] == x && ia < N_aCount) {
				rdiff -= sa;
				ia++;
			}
			while (Ra_bData[ib] == x && ib < N_bCount) {
				rdiff += sb;
				ib++;
			}
			if (ia >= N_aCount) { ok = true; break; }
			if (ib >= N_bCount) { ok = true; break; }
		}
		rdmax = W_MAX(rdmax, fabs(rdiff));
		//rdmax = TMath::Max(rdmax,TMath::Abs(rdiff));
	}
	if (ok == false) {
		LOG("ERROR!\n");
		exit(1);
	}

	if (ok) {
		rdmax = W_MAX(rdmax, fabs(rdiff));
		double z = rdmax * sqrt(rna*rnb / (rna + rnb));
		prob = PVkolmogorov(z);
	}
	// debug printout
	//	if (opt.Contains("D")) {
	printf(" Kolmogorov Probability = %g, Max Dist = %g\n", prob, rdmax);
	//	}
	return (wsReal)prob;
	//	if(opt.Contains("M")) return rdmax;
	//	else                  return prob;
}

/*
 * Rank-sum test
 */

/* ascend-order sorting function for phylop score */
int asc_sort_phyloP(const void *X_l, const void *X_r)
{
	double _a = ((RANK *)X_l)->R_score;
	double _b = ((RANK *)X_r)->R_score;

	return _a > _b ? 1 : (_a < _b ? -1 : 0);
}

/* ����-�� �׽�Ʈ */
double ranksum(int N_cntX, double *Ra_X, int N_cntY, double *Ra_Y)
{
	int		N_cntSum		= N_cntX+N_cntY;
	double	R_cntX			= (double)N_cntX;
	double	R_cntY			= (double)N_cntY;
	double	R_cntSum		= R_cntX+R_cntY;
	RANK	*Sa_ranks		= NULL;
	int		*Na_tieSizes	= NULL;
	double	R_Wval			= 0.0;
	int		N_tieCount		= 0;
	int		N_sameRankCount;
	double	R_p, R_q;
	//	double	R_Ua;
	int		i, j;
	/* Z-norm ���� */
	double	R_TjSum = 0.0;
	double	R_varW, R_Znorm;
	double	R_clipVal = 0.5;

	wsAlloc(Sa_ranks, RANK, N_cntSum);
	wsAlloc(Na_tieSizes, int, MAX_TIES_CNT);

	for (i=j=0 ; i<N_cntX ; i++,j++) {
		Sa_ranks[i].R_score = Ra_X[i];
		Sa_ranks[j].B_isY = 0;
	}
	for (i=0 ; i<N_cntY ; i++,j++) {
		Sa_ranks[j].R_score = Ra_Y[i];
		Sa_ranks[j].B_isY = 1;
	}
	qsort(Sa_ranks, N_cntSum, sizeof(RANK), asc_sort_phyloP);

	//	printf("Sorted result : ");
	//	for (i=0 ; i<cs ; i++)
	//		printf("%lf ", s[i].v);
	//	printf("\n");

	/* ������ ����Ͽ� ��ũ �ű�� */
	for (i=1 ; i<=N_cntSum ; i+=N_sameRankCount) {
		/* �ϴ� �ּ��� �ϳ��� �����Ƿ� */
		N_sameRankCount = 1;

		/* �ټ�����? */
		for (j=i-1 ; Sa_ranks[j].R_score == Sa_ranks[j+1].R_score ; j++)
			N_sameRankCount++;

		if (N_sameRankCount > 1) { /* var ���ϱ� ���� �����׷� ��� */
			if (N_tieCount == MAX_TIES_CNT) {
				printf("Too many ties!\n");
				exit(1);
			}

			Na_tieSizes[N_tieCount++] = N_sameRankCount;
		}

		if (N_sameRankCount == 1)
			Sa_ranks[i-1].rank = (double)i;
		else for (j=0 ; j<N_sameRankCount ; j++)
			Sa_ranks[i+j-1].rank = (i+(i+N_sameRankCount-1))/2.0;
	}

	//	printf("Ranked result : \n");
	//	for (i=0 ; i<cs ; i++)
	//		printf("\t%c %lf (%lf)\n", s[i].isy?'Y':'X', s[i].v, s[i].rank);
	//	printf("\n");

	/* W ��� */
	for (i=0 ; i<N_cntSum ; i++)
		if (!Sa_ranks[i].B_isY) R_Wval += Sa_ranks[i].rank;

	//	R_Ua = (R_cntX*R_cntY) + (R_cntX*(R_cntX+1.0) / (double)2.0) - R_Wval;

	/* R_Znorm ��� */
	for (i=0 ; i<N_tieCount ; i++)
		R_TjSum += Na_tieSizes[i] * (Na_tieSizes[i] * Na_tieSizes[i] - 1.0);
	R_varW = (R_cntX*R_cntY)/12.0 * (R_cntSum+1.0 - R_TjSum/(R_cntSum*(R_cntSum-1.0)));
	//R_Znorm = (W-(cx*(cs+1))/2.0) / sqrt(R_varW);
	if (R_Wval-(R_cntX*(R_cntSum+1)/2.0) > 0) R_clipVal *= -1;
	R_Znorm = (R_Wval+R_clipVal-(R_cntX*(R_cntSum+1))/2.0) / sqrt(R_varW);
	//R_Znorm = 0.63;
	//	printf("R_Ua %lf W %lf, R_Znorm %lf\n", R_Ua, R_Wval, R_Znorm);

	/* p-value ��� */
	x_cumnor(&R_Znorm, &R_p, &R_q);

	DEALLOC(Na_tieSizes);
	DEALLOC(Sa_ranks);

	return R_p>R_q?R_q*2:R_p*2;
}

double betaP(double x, double pin, double qin, int lower_tail, int log_p)
{
	if (pin <= 0 || qin <= 0) numeric_limits<double>::quiet_NaN();

	if (x <= 0)
		return (lower_tail ? (log_p ? ML_NEGINF : 0.) : (log_p ? 0. : 1.));
	if (x >= 1)
		return (lower_tail ? (log_p ? 0. : 1.) : (log_p ? ML_NEGINF : 0.));
	return x_betaPraw(x, pin, qin, lower_tail, log_p);
}

/*
 *
 * Mixed chi-square p-value
 * 
 */

static double  log1(double x, char first)
	/* if (first) log(1 + x) ; else  log(1 + x) - x */
{
	if (fabs(x) > 0.1)
	{
		return (first ? log(1.0 + x) : (log(1.0 + x) - x));
	}
	else
	{
		double s, s1, term, y, k;
		y = x / (2.0 + x);  term = 2.0 * (y*y*y);  k = 3.0;
		s = (first ? 2.0 : - x) * y;
		y = SQR(y);
		for (s1 = s + term / k; s1 != s; s1 = s + term / k)
		{ k = k + 2.0; term = term * y; s = s1; }
		return s;
	}
}

static double exp1(double x)               /* to avoid underflows  */
{ return x < -50.0 ? 0.0 : exp(x); }

static jmp_buf	davies_env;
static char		davies_fail;
static int		davies_count;
static char		davies_ndtsrt;
static double davies_sigsq, davies_lmax, davies_lmin, davies_mean;
static double davies_intl,davies_ersm;
static int *davies_th;
static int davies_lim;
static double davies_c;
static int *davies_n, davies_r;
static double *davies_lb,*davies_nc;

#define pi		3.1415926535897932384
#define log28	0.08664339756999316077835  /*  log(2.0) / 8.0  */

static void davies_order(void)
	/* find order of absolute values of lb */
{
	int		j, k;
	double	lj;

	for (j=0 ; j<davies_r ; j++) {
		lj = fabs(davies_lb[j]);

		for (k=j-1 ; k>=0 ; k--) {
			if ( lj > fabs(davies_lb[davies_th[k]]) )
				davies_th[k + 1] = davies_th[k];
			else goto l1;
		}
		k = -1;
l1 :
		davies_th[k+1] = j;
	}
	davies_ndtsrt = 0;
}

static void davies_counter(void)
	/*  count number of calls to errbd, truncation, cfe */
{
	if ( ++davies_count > davies_lim ) longjmp(davies_env, 1);
}

static double   davies_errbd(double u, double* cx)
	/*  find bound on tail probability using mgf, cutoff
	point returned to *cx */
{
	double	sum1, lj, ncj, x, y, xconst;
	int		j, nj;

	davies_counter();
	xconst	= u*davies_sigsq;
	sum1	= u*xconst;
	u		= 2.0 * u;

	for (j=davies_r-1 ; j>=0 ; j--) {
		nj		= davies_n[j];
		lj		= davies_lb[j];
		ncj		= davies_nc[j];
		x		= u * lj;
		y		= 1.0 - x;
		xconst	+= lj*(ncj/y + nj)/y;
		sum1	+= ncj*SQR(x/y) + nj*(SQR(x)/y + log1(-x, 0));
	}
	*cx = xconst;
	
	return exp1(-0.5*sum1);
}

static double  davies_ctff(double accx, double* upn)
	/*  find ctff so that p(qf > ctff) < accx  if (upn > 0,
	p(qf < ctff) < accx otherwise */
{
	double u1, u2, u, rb, xconst, c1, c2;
	u2 = *upn;
	u1 = 0.0;
	c1 = davies_mean;
	rb = 2.0 * ((u2 > 0.0) ? davies_lmax : davies_lmin);
	for (u=u2/(1.0 + u2*rb) ; davies_errbd(u, &c2)>accx ; u=u2/(1.0 + u2*rb)) {
		u1 = u2;
		c1 = c2;
		u2 = 2.0*u2;
	}
	for (u=(c1-davies_mean) / (c2-davies_mean) ; u<0.9 ; u=(c1-davies_mean) / (c2-davies_mean)) {
		u = (u1+u2) / 2.0;
		if (davies_errbd(u / (1.0 + u * rb), &xconst) > accx) {
			u1 = u;
			c1 = xconst;
		} else {
			u2 = u;
			c2 = xconst;
		}
	}
	*upn = u2;
	return c2;
}

static double davies_truncation(double u, double tausq)
	/* bound integration error due to truncation at u */
{
	double sum1, sum2, prod1, prod2, prod3, lj, ncj,
		x, y, err1, err2;
	int j, nj, s;

	davies_counter();
	sum1	= 0.0;
	prod2	= 0.0;
	prod3	= 0.0;
	s		= 0;
	sum2	= (davies_sigsq+tausq) * SQR(u);
	prod1	= 2.0 * sum2;
	u		*= 2.0;

	for (j=0 ; j<davies_r ; j++) {
		lj		= davies_lb[j];
		ncj		= davies_nc[j];
		nj		= davies_n[j];
		x		= SQR(u * lj);
		sum1	+= ncj*x / (1.0+x);

		if (x > 1.0) {
			prod2	+= nj * log(x);
			prod3	+= nj * log1(x, TRUE);
			s		+= nj;
		} else
			prod1	+= nj * log1(x, TRUE);
	}
	sum1	*= 0.5;
	prod2	+= prod1;
	prod3	+= prod1;

	x		= exp1(-sum1 - 0.25*prod2) / pi;
	y		= exp1(-sum1 - 0.25*prod3) / pi;
	err1	= (s == 0) ? 1.0 : x*2.0 / s;
	err2	= (prod3 > 1.0) ? 2.5*y : 1.0;
	if (err2 < err1)
		err1 = err2;
	x		= 0.5 * sum2;
	err2	= x<=y ? 1.0 : y/x;

	return err1<err2 ? err1 : err2;
}

static void davies_findu(double* utx, double accx)
	/*  find u such that truncation(u) < accx and truncation(u / 1.2) > accx */
{
	int i;
	double u, ut;
	static double divis[] = { 2.0, 1.4, 1.2, 1.1 };

	ut	= *utx;
	u	= ut / 4.0;
	if (davies_truncation(u, 0.0) > accx) {
		for (u=ut ; davies_truncation(u, 0.0)>accx ; u=ut)
			ut *= 4.0;
	} else {
		ut = u;
		for (u/=4.0 ; davies_truncation(u, 0.0)<=accx ; u/=4.0 )
			ut = u;
	}
	for (i=0 ; i<4 ; i++) {
		u = ut/divis[i];
		if (davies_truncation(u, 0.0) <= accx )
			ut = u;
	}
	*utx = ut;
}


static void davies_integrate(int nterm, double interv, double tausq, char mainx)
	/*  carry out integration with nterm terms, at stepsize
	interv.  if (! mainx) multiply integrand by
	1.0-exp(-0.5*tausq*u^2) */
{
	double inpi, u, sum1, sum2, sum3, x, y, z;
	int k, j, nj;
	inpi = interv / pi;

	for (k=nterm ; k>=0 ; k--) {
		u		= (k+0.5) * interv;
		sum1	= - 2.0 * u * davies_c;
		sum2	= fabs(sum1);
		sum3	= - 0.5 * davies_sigsq * SQR(u);

		for (j=davies_r-1 ; j>=0 ; j--) {
			nj		= davies_n[j];
			x		= 2.0 * davies_lb[j] * u;
			y		= SQR(x);
			sum3	-= 0.25 * nj * log1(y, TRUE );
			y		= davies_nc[j]*x / (1.0+y);
			z		= nj*atan(x) + y;
			sum1	+= z;
			sum2	+= fabs(z);
			sum3	-= 0.5 * x * y;
		}

		x = inpi * exp1(sum3) / u;
		if (!mainx)
			x *= 1.0 - exp1(-0.5 * tausq * SQR(u));
		sum1 = sin(0.5 * sum1) * x;
		sum2 *= 0.5*x;
		davies_intl += sum1;
		davies_ersm += sum2;
	}
}

static double davies_cfe(double x)
	/*  coef of tausq in error when convergence factor of
	exp1(-0.5*tausq*u^2) is used when df is evaluated at x */
{
	double axl, axl1, axl2, sxl, sum1, lj;
	int j, k, t;

	davies_counter();
	if (davies_ndtsrt)
		davies_order();

	axl		= fabs(x);
	sxl		= x>0.0 ? 1.0 : -1.0;
	sum1	= 0.0;

	for (j=davies_r-1 ; j>=0 ; j--) {
		t = davies_th[j];
		if (davies_lb[t]*sxl > 0.0) {
			lj		= fabs(davies_lb[t]);
			axl1	= axl - lj * (davies_n[t] + davies_nc[t]);
			axl2	= lj / log28;

			if (axl1 > axl2)
				axl = axl1;
			else {
				if (axl > axl2)
					axl = axl2;
				sum1 = (axl-axl1) / lj;
				for (k=j-1 ; k>=0 ; k--)
					sum1 += davies_n[davies_th[k]] + davies_nc[davies_th[k]];
				goto l;
			}
		}
	}
l:
	if (sum1 > 100.0) {
		davies_fail = TRUE;
		return 1.0;
	} else
		return pow(2.0, sum1/4.0) / (pi*SQR(axl));
}

double qfc(double* lb1, double* nc1, int* Na_df, int N_lambda, double R_sigma,
	double *c1, int N_lim, double R_acc, double* trace, int* ifault)

	/*  distribution function of a linear combination of non-central
	chi-squared random variables :

	input:
	lb[j]            coefficient of j-th chi-squared variable
	nc[j]            non-centrality parameter
	Na_df[j]         degrees of freedom
	j = 0, 2 ... N_lambda-1
	R_sigma          coefficient of standard normal variable
	c                point at which df is to be evaluated
	lim              maximum number of terms in integration
	acc              maximum error

	output:
	ifault = 1       required accuracy NOT achieved
	2       round-off error possibly significant
	3       invalid parameters
	4       unable to locate integration parameters
	5       out of memory

	trace[0]         absolute sum
	trace[1]         total number of integration terms
	trace[2]         number of integrations
	trace[3]         integration interval in final integration
	trace[4]         truncation point in initial integration
	trace[5]         s.d. of initial convergence factor
	trace[6]         cycles to locate integration parameters     */

{
	int			j, nj, nt, ntm;
	double		acc1, almx, xlim, xnt, xntm;
	double		utx, tausq, sd, intv, intv1, x, up, un, d1, d2, lj, ncj;
	double		qfval = -1.0;
	static int	rats[]={ 1, 2, 4, 8 };

	if (setjmp(davies_env) != 0)
		halt("Failed to set long jump in Davies method");

	/* Externs */
	davies_lim	= N_lim;
	davies_r	= N_lambda;
	davies_c	= c1[0];
	davies_n	= Na_df;
	davies_lb	= lb1;
	davies_nc	= nc1;

	for (j=0 ; j<7 ; j++)
		trace[j] = 0.0;

	*ifault = 0;
	davies_count = 0;
	davies_intl = 0.0;
	davies_ersm = 0.0;
	qfval = -1.0;
	acc1 = R_acc;
	davies_ndtsrt = 1;
	davies_fail = 0;
	xlim = (double)davies_lim;

	wsAlloc(davies_th, int, N_lambda);

	/* find mean, sd, max and min of lb,
	check that parameter values are valid */
	davies_sigsq = SQR(R_sigma);
	sd = davies_sigsq;
	davies_lmax = 0.0; davies_lmin = 0.0; davies_mean = 0.0;

	for (j=0 ; j<N_lambda ; j++) {
		nj	= davies_n[j];
		lj	= davies_lb[j];
		ncj	= davies_nc[j];
		if ( nj < 0  ||  ncj < 0.0 ) { *ifault = 3;  goto  endofproc;  }
		sd  = sd  + SQR(lj) * (2 * nj + 4.0 * ncj);
		davies_mean = davies_mean + lj * (nj + ncj);
		if (davies_lmax < lj) davies_lmax = lj ;
		else if (davies_lmin > lj) davies_lmin = lj;
	}

	if ( sd == 0.0  ) {
		qfval = (davies_c > 0.0) ? 1.0 : 0.0;
		goto  endofproc;
	}

	/* Check parameter */
	if (davies_lmin==0.0 && davies_lmax==0.0 && R_sigma==0.0) {
		*ifault = 3;
		goto  endofproc;
	}

	sd		= sqrt(sd);
	almx	= (davies_lmax < -davies_lmin) ? -davies_lmin : davies_lmax;

	/* starting values for findu, ctff */
	utx	= 16.0 / sd;
	up	= 4.5 / sd;
	un	= -up;

	/* truncation point with no convergence factor */
	davies_findu(&utx, .5*acc1);

	/* does convergence factor help */
	if (davies_c!=0.0 && (almx > 0.07*sd)) {
		tausq = .25 * acc1 / davies_cfe(davies_c);

		if (davies_fail)
			davies_fail = 0;
		else if (davies_truncation(utx, tausq) < .2*acc1) {
			davies_sigsq += tausq;
			davies_findu(&utx, .25*acc1);
			trace[5] = sqrt(tausq);
		}
	}
	trace[4] = utx;  acc1 = 0.5 * acc1;

	/* find RANGE of distribution, quit if outside this */
l1:
	d1 = davies_ctff(acc1, &up) - davies_c;
	if (d1 < 0.0) {
		qfval = 1.0;
		goto endofproc;
	}

	d2 = davies_c - davies_ctff(acc1, &un);
	if (d2 < 0.0) {
		qfval = 0.0;
		goto endofproc;
	}

	/* find integration interval */
	intv = 2.0 * pi / ((d1 > d2) ? d1 : d2);

	/* calculate number of terms required for main and auxilary integrations */
	xnt		= utx / intv;
	xntm	= 3.0 / sqrt(acc1);

	if (xnt > xntm*1.5) {
		/* parameters for auxillary integration */
		if (xntm > xlim) { *ifault = 1; goto endofproc; }
		ntm		= (int)floor(xntm+0.5);
		intv1	= utx / ntm;  x = 2.0 * pi / intv1;
		if (x <= fabs(davies_c))
			goto l2;

		/* calculate convergence factor */
		tausq = .33 * acc1 / (1.1 * (davies_cfe(davies_c - x) + davies_cfe(davies_c + x)));
		if (davies_fail)
			goto l2;
		acc1 = .67 * acc1;

		/* auxillary integration */
		davies_integrate(ntm, intv1, tausq, 0 );
		xlim -= xntm;
		davies_sigsq += tausq;
		trace[2] += 1;
		trace[1] += ntm + 1;

		/* find truncation point with new convergence factor */
		davies_findu(&utx, .25 * acc1);
		acc1 = 0.75*acc1;

		goto l1;
	}

	/* main integration */
l2:
	trace[3] = intv;
	if (xnt > xlim) { *ifault = 1; goto endofproc; }
	nt = (int)floor(xnt+0.5);
	davies_integrate(nt, intv, 0.0, TRUE );
	trace[2]	+= 1;
	trace[1]	+= nt + 1;
	qfval		= 0.5 - davies_intl;
	trace[0]	= davies_ersm;

	/* test whether round-off error could be significant
	allow for radix 8 or 16 machines */
	up	= davies_ersm;
	x	= up + R_acc / 10.0;
	for (j=0 ; j<4 ; j++) {
		if (rats[j]*x == rats[j]*up)
			*ifault = 2;
	}

endofproc:
	DEALLOC(davies_th);
	trace[6]	= (double)davies_count;

	if (qfval < 0.0)
		return WISARD_NAN;
	return 1.0 - qfval;
}

double davies(double R_Q, double *Ra_lambda, wsUint N_lambda)
{
	int		*Na_df	= NULL; wsAlloc(Na_df, int, N_lambda);
	double	*Ra_nc	= NULL; wsAlloc(Ra_nc, double, N_lambda);
	double	*Ra_Q	= NULL; wsAlloc(Ra_Q, double, N_lambda);
	double	Ra_trace[7];
	int		N_fault;
	for (wsUint i=0 ; i<N_lambda ; i++) {
		Ra_Q[i] = R_Q;
		Na_df[i] = 1;
		Ra_nc[i] = 0.0;
	}

	double	R_Qq = qfc(Ra_lambda, Ra_nc, Na_df, N_lambda, 0.0,
		Ra_Q, 10000, 1e-04, Ra_trace, &N_fault);

	DEALLOC(Na_df);
	DEALLOC(Ra_nc);
	DEALLOC(Ra_Q);

	return R_Qq!=R_Qq?-1.0:R_Qq;
}

double farebrother(double R_Q, double *Ra_lambda, wsUint N_lambda)
{
	int		*Na_df	= NULL; wsAlloc(Na_df, int, N_lambda);
	double	*Ra_nc	= NULL; wsAlloc(Ra_nc, double, N_lambda);
	for (wsUint i=0 ; i<N_lambda ; i++) {
		Na_df[i] = 1;
		Ra_nc[i] = 0.0;
	}
	double	mode	= 1.0;
	int		maxit	= 100000;
	double	eps		= 10e-10;
	double	dnsty	= 0.0;
	int		ifault	= 0;

	// 	h = rep(1, length(lambda)), delta = rep(0, 
	// 		length(lambda)), maxit = 1e+05, eps = 10^(-10), mode = 1) 
	// 	{
	// 		r <- length(lambda)
	// 			if (length(h) != r) 
	// 				stop("lambda and h should have the same length!")
	// 				if (length(delta) != r) 
	// 					stop("lambda and delta should have the same length!")
	// 					dnsty <- 0
	// 					ifault <- 0
	// 					res <- 0
	wsReal out = ruben(Ra_lambda, Na_df, Ra_nc, N_lambda, R_Q, &mode,
		&maxit, &eps, &dnsty, &ifault);
	// 	lambda = as.double(lambda), h = as.integer(h), 
	// 					delta = as.double(delta), r = as.integer(r), q = as.double(q), 
	// 					mode = as.double(mode), maxit = as.integer(maxit), eps = as.double(eps), 
	// 					dnsty = as.double(dnsty), ifault = as.integer(ifault), 
	// 					res = as.double(res), PACKAGE = "CompQuadForm")
	return 1.0 - out;
	//		out$res <- 1 - out$res
	//		return(out)
}

double ruben(double *lambda, int *mult, double *delta, int n, double c,
			 double *mode, int *maxit, double *eps, double *dnsty, int *ifault)
{
	double res = WISARD_NAN;
	int i,k,m,j;
	double ao, aoinv, z, bbeta, eps2, hold, hold2, sum, sum1, dans, lans, pans, prbty, tol;
	double *gamma, *theta, *a, *b;
	gamma	= new double[n];
	theta	= new double[n];
	a		= new double[maxit[0]];
	b		= new double[maxit[0]];



	if ((n<1) || (c<=0) || (maxit[0] <1) || (eps[0]<=0.0)) {
		/* ERROR */
		ifault[0] = 2;
		delete[] gamma;
		delete[] theta;
		delete[] a;
		delete[] b;

		return res;
	} 
	tol = -200.0;

	// Preliminaries
	sum = lambda[0];
	bbeta = sum;

	for (i=0 ; i<n ; i++) {
		hold = lambda[i];
		if ((hold<=0.0) || (mult[i]<1) || (delta[i]<0.0)) {
			/* ERROR */
			ifault[0] = -(int)(i+1);
			delete[] gamma;
			delete[] theta;
			delete[] a;
			delete[] b;
			return res;
		}	
		if (bbeta > hold) bbeta = hold; // calcul du max des lambdas
		if (sum < hold) sum = hold;    // calcul du min des lambdas
	}


	if (mode[0] > 0.0) {
		// if ((2.0/(1.0/bbeta+1.0/sum))>1.8*sum) bbeta = sum; // comme dans NAG : methode avec betaA
		bbeta = mode[0]*bbeta;
	} else {
		bbeta = 2.0/(1.0/bbeta+1.0/sum);  // methode avec betaB
	}

	k = 0;
	sum = 1.0;
	sum1 = 0.0;
	for (i=0 ; i<n ; i++) {
		hold = bbeta/lambda[i];
		gamma[i] = 1.0 - hold;
		sum = sum * pow(hold, mult[i]); //???? pas sur ..
		sum1 = sum1 + delta[i];
		k = k + mult[i];
		theta[i] = 1.0;
	}

	ao = exp(0.5*(log(sum)-sum1));
	if (ao <= 0.0) {
		res = 0.0;
		dnsty[0] = 0.0;
		ifault[0] = 1;
	} else { // evaluate probability and density of chi-squared on k degrees of freedom. The constant 0.22579135264473 is ln(sqrt(pi/2))
		z = c/bbeta;

		if ((k%2) == 0) { // k est un entier donc on regarde si k est divisible par 2: k == (k/2)*k 
			i = 2;
			lans = -0.5*z;
			dans = exp(lans);
			pans = 1.0 - dans;
		} else {
			i = 1;
			lans = -0.5*(z+log(z)) - 0.22579135264473;
			dans = exp(lans);
			pans = pnorm(sqrt(z),0.0,1.0,1,0) - pnorm(-sqrt(z),0.0,1.0,1,0); 
		}

		k = k-2;
		for (j=i;j<=k;j=j+2) {
			if (lans < tol) {
				lans = lans + log(z/(double)j);
				dans = exp(lans);
			} else {
				dans = dans*z/(double)j;
			}
			pans = pans -dans;
		}

		// evaluate successive terms of expansion

		prbty = pans;
		dnsty[0] = dans;
		eps2 = eps[0]/ao;
		aoinv = 1.0/ao;
		sum = aoinv - 1.0;


		for (m=1;m<=maxit[0];m++) {
			sum1 = 0.0;
			for (i=0 ; i<n ; i++) {
				hold = theta[i];
				hold2 = hold*gamma[i];
				theta[i] = hold2;
				sum1 += hold2*mult[i]+m*delta[i]*(hold-hold2);
			}
			sum1 = 0.5*sum1;
			b[m-1] = sum1;
			for (i=m-1 ; i>=1 ; i--) {
				sum1 += b[i-1]*a[m-i-1]; 
			}
			sum1 = sum1/(double)m;
			a[m-1] = sum1;
			k = k + 2;
			if (lans < tol) {
				lans = lans + log(z/(double)k);
				dans = exp(lans);
			} else {
				dans = dans*z/(double)k;
			}
			pans = pans - dans;
			sum = sum - sum1;
			dnsty[0] = dnsty[0] + dans*sum1;
			sum1 = pans*sum1;
			prbty = prbty + sum1;
			if (prbty<(-aoinv)) {
				/* ERROR */
				res = WISARD_NAN;
				delete[] gamma;
				delete[] theta;
				delete[] a;
				delete[] b;
				ifault[0] = 3;
				return res;
			}
			if (fabs(pans*sum) < eps2) {
				if (fabs(sum1) < eps2) {
					ifault[0] = 0;

					m = maxit[0]+1;
					break;

				}
			}
		}

		ifault[0] = 4;
		dnsty[0] = ao*dnsty[0]/(bbeta+bbeta);
		prbty = ao*prbty;
		if (prbty<0.0 || prbty>1.0) {
			ifault[0] = ifault[0] + 5;
		} else {
			if (dnsty[0]<0.0) ifault[0] = ifault[0] + 6;
		}
		res = prbty;
	}

	delete[] gamma;
	delete[] theta;
	delete[] a;
	delete[] b;

	return res;
}

double	metaFisher(wsUint N_pv, wsVec Ra_pv)
{
	wsReal R_sum = 0;
	LOOP(i, N_pv) R_sum += log(Ra_pv[i]);
	R_sum *= -W2;

	return PVchisq(R_sum, 2 * N_pv);
}

double	metaStouffer(wsUint N_pv, wsVec Ra_pv, wsVec Ra_weight/*=NULL*/)
{
	wsReal R_sum = W0;
	// zp <- (qnorm(p[keep], lower.tail = FALSE) %*% weights[keep])/sqrt(sum(weights[keep]^2))
	if (!Ra_weight) {
		LOOP(i, N_pv) R_sum += qnorm(Ra_pv[i], 0.0, 1.0, 0);
		R_sum /= sqrt((wsReal)N_pv);
	} else {
		wsReal R_denom = W0;
		LOOP(i, N_pv) {
			R_sum += qnorm(Ra_pv[i], 0.0, 1.0, 0) * Ra_weight[i];
			R_denom += Ra_weight[i] * Ra_weight[i];
		}
		R_sum /= sqrt(R_denom);
	}
	return pnorm(R_sum, 0.0, 1.0, 0);
}

double	metaKost(wsUint N_pv, wsVec Ra_pv, wsSym Ra_cor)
{
	wsReal R_cov = W0;
	for (wsUint i=1 ; i<N_pv ; i++) for (wsUint j=0 ; j<i ; j++) {
		// 3.263*cor.offdiag + 0.71*cor.offdiag ^ 2 + 0.027*cor.offdiag ^ 3
		wsReal R_cor = Ra_cor[i][j];
		R_cov += 3.263*R_cor + 0.71*R_cor*R_cor + 0.027*R_cor*R_cor*R_cor;
	}
	// var.psi <- 4*k + 2 * sum(3.263*cor.offdiag + 0.71*cor.offdiag^2 + 0.027*cor.offdiag^3)
	wsReal varPsi = (wsReal)(4*N_pv) + W2*R_cov;
	// e.psi <- 2*k
	wsReal ePsi = (wsReal)(2 * N_pv);
	// f <- e.psi^2 / var.psi
	wsReal f = ePsi*ePsi / varPsi;
	// c <- var.psi / (2*e.psi)
	wsReal c = varPsi / (W2*ePsi);

	// s.psi <--2 * sum(log(cur.pval))
	wsReal sPsi = W0;
	LOOP(i, N_pv) sPsi += log(Ra_pv[i]);
	sPsi *= -W2;

	// stat <- s.psi / c
	wsReal stat = sPsi / c;

	return PVchisq(stat, W2*f);
}

wsReal test2x2(wsUint N_11, wsUint N_12, wsUint N_21, wsUint N_22,
			   char B_corr)
{
	wsUint N_r1 = N_11 + N_12;
	wsUint N_r2 = N_21 + N_22;
	wsUint N = N_r1+N_r2;
	//wsReal NNN = N*N*N;
	wsUint N_c1 = N_11 + N_21;
	wsUint N_c2 = N_12 + N_22;

	wsReal E_11 = (wsReal)(N_r1*N_c1)/(wsReal)N;
	wsReal E_12 = (wsReal)(N_r1*N_c2)/(wsReal)N;
	wsReal E_21 = (wsReal)(N_r2*N_c1)/(wsReal)N;
	wsReal E_22 = (wsReal)(N_r2*N_c2)/(wsReal)N;

	if (E_11 < 5 || E_12 < 5 || E_21 < 5 || E_22 < 5)
		return fisher(N_11, N_12, N_21, N_22);

	wsReal R_Yates = B_corr ? REAL_CONST(0.5) : W0;
	// 	wsReal v_11 = N_c1*N_r1*(N-N_r1)*(N-N_c1)/NNN;
	// 	wsReal v_12 = N_c2*N_r1*(N-N_r1)*(N-N_c2)/NNN;
	// 	wsReal v_21 = N_c1*N_r2*(N-N_r2)*(N-N_c1)/NNN;
	// 	wsReal v_22 = N_c2*N_r2*(N-N_r2)*(N-N_c2)/NNN;

	// 	if (correct && nrow(x) == 2 && ncol(x) == 2) {
	// 		YATES <- 0.5
	// 			METHOD <- paste(METHOD, "with Yates' continuity correction")
	// 	}
	// 	else YATES <- 0
	// 		STATISTIC <- sum((abs(x - E) - YATES)^2/E)
	// 		PARAMETER <- (nr - 1) * (nc - 1)
	// 		PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)

	wsReal R_stChi =
		SQR(fabs(N_11-E_11) - R_Yates) / E_11 +
		SQR(fabs(N_12-E_12) - R_Yates) / E_12 +
		SQR(fabs(N_21-E_21) - R_Yates) / E_21 +
		SQR(fabs(N_22-E_22) - R_Yates) / E_22;

	return (wsReal)PVchisq(R_stChi, 1.0);
}

wsReal test2x2_PEDCMC(wsUint N_11, wsUint N_12, wsUint N_21, wsUint N_22,
			   char B_corr, wsUint *Np_typeTest)
{
	wsUint N_r1 = N_11 + N_12;
	wsUint N_r2 = N_21 + N_22;
	wsUint N = N_r1+N_r2;
	//wsReal NNN = N*N*N;
	wsUint N_c1 = N_11 + N_21;
	wsUint N_c2 = N_12 + N_22;

	wsReal E_11 = (wsReal)(N_r1*N_c1)/(wsReal)N;
	wsReal E_12 = (wsReal)(N_r1*N_c2)/(wsReal)N;
	wsReal E_21 = (wsReal)(N_r2*N_c1)/(wsReal)N;
	wsReal E_22 = (wsReal)(N_r2*N_c2)/(wsReal)N;

	if (E_11 < 5 || E_12 < 5 || E_21 < 5 || E_22 < 5) {
		*Np_typeTest = 1;
		return fisher(N_11, N_12, N_21, N_22);
	}
	*Np_typeTest = 0;

	wsReal R_Yates = B_corr ? REAL_CONST(0.5) : W0;
	// 	wsReal v_11 = N_c1*N_r1*(N-N_r1)*(N-N_c1)/NNN;
	// 	wsReal v_12 = N_c2*N_r1*(N-N_r1)*(N-N_c2)/NNN;
	// 	wsReal v_21 = N_c1*N_r2*(N-N_r2)*(N-N_c1)/NNN;
	// 	wsReal v_22 = N_c2*N_r2*(N-N_r2)*(N-N_c2)/NNN;

	// 	if (correct && nrow(x) == 2 && ncol(x) == 2) {
	// 		YATES <- 0.5
	// 			METHOD <- paste(METHOD, "with Yates' continuity correction")
	// 	}
	// 	else YATES <- 0
	// 		STATISTIC <- sum((abs(x - E) - YATES)^2/E)
	// 		PARAMETER <- (nr - 1) * (nc - 1)
	// 		PVAL <- pchisq(STATISTIC, PARAMETER, lower.tail = FALSE)

	wsReal R_stChi =
		SQR(fabs(N_11-E_11) - R_Yates) / E_11 +
		SQR(fabs(N_12-E_12) - R_Yates) / E_12 +
		SQR(fabs(N_21-E_21) - R_Yates) / E_21 +
		SQR(fabs(N_22-E_22) - R_Yates) / E_22;

	return R_stChi;
}

} // End namespace ONETOOL
