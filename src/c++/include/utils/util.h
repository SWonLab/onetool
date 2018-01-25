#ifndef __WISARD_UTIL_H__
#define __WISARD_UTIL_H__
#pragma once

#ifndef _WIN32
#	include <sys/time.h>
#	include <pthread.h>
#	include <unistd.h>
#	define pthread_mutex_init(a)	pthread_mutex_init(&a, NULL)
#	define pthread_mutex_uninit(a)
#else
#	include <windows.h>
typedef CRITICAL_SECTION pthread_mutex_t;
#	define pthread_mutex_init(a)	InitializeCriticalSection(&a)
#	define pthread_mutex_uninit(a)	DeleteCriticalSection(&a)
#	define pthread_mutex_lock		EnterCriticalSection
#	define pthread_mutex_unlock		LeaveCriticalSection
#endif

#define str2dbl	ws_str2dbl

#define ML_POSINF numeric_limits<double>::infinity()
#define ML_NEGINF -numeric_limits<double>::infinity()

#ifndef SSEINLINE
#	if __GNUC__ && !__GNUC_STDC_INLINE__
#		define SSEINLINE extern inline
#	else
#		define SSEINLINE inline
#	endif
#endif

#define SET_AVAIL() (IS_ASSIGNED(setconsec) || IS_ASSIGNED(set) || IS_ASSIGNED(setrandom))

#include "global/common.h"
#include "global/datatype.h"
#ifdef _WIN32
#	include <process.h>
#else
#	include <unistd.h>
#	define _execlp execlp
#	define _aligned_free ::free
#	define MAX_PATH 512
#	define _finite isfinite
#endif
#include <string>
#include <vector>
#include <map>
#include <unordered_set>
#include <exception>
#include <limits>
#ifndef __APPLE__
#	ifdef __FreeBSD__
#		include <stdlib.h>
#	else
#		include <malloc.h>
#	endif
#endif
#include <stdlib.h>
#include <limits>
#include <float.h>
#include <limits.h>
#include "output/exporter.h"
#include "utils/comb.h"
#include "utils/log.h"
#include "utils/simd.h"
#include "utils/rand.h"
#include "utils/timer.h"
#include "utils/genet.h"
#include "utils/sort.h"
#include "utils/memmgr.h"

#ifdef USE_MEMDBG
#	include "memmgr.h"
#endif

#define NO_MULTITHREAD (OPT_NUMBER(thread) == 1)
#define IS_MULTITHREAD (OPT_NUMBER(thread) > 1)

namespace ONETOOL {

/* Defines
 *
 */

typedef vector<std::string>	vStr;
typedef vStr::iterator		vStr_it;

/* iterator macro for STL foreach */
#define	FOREACH(t,c,i) for (t i=(c).begin() ; i!=(c).end() ; i++)
#define	RFOREACH(t,c,i) for (t i=(c).rbegin() ; i!=(c).rend() ; i++)
#define	FOREACHDO(t,c,i,r) for (t i=(c).begin() ; i!=(c).end() ; i++,r)

/* loop */
#define	LOOP(i,to) for (wsUint i=0 ; i<to ; i++)
#define	LOOPi(i,to) for (int i=0 ; i<to ; i++)
#define LOOPV(i, v, to) for (wsUint i=v ; i<to ; i++)

#define PVAL_FAIL(p) (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), (p)))
#define PVAL_FAIL2(p1, p2) (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), (p1)) && !isInRange(OPT_RANGE(pvalrange), (p2)))
#define PVAL_FAIL3(p1, p2, p3) (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), (p1)) && !isInRange(OPT_RANGE(pvalrange), (p2)) && !isInRange(OPT_RANGE(pvalrange), (p3)))
 
/* t ????? ???? sz?? ????? ????? ??? ?? p?? ??????? ??? */
#ifdef USE_MEMDBG
#	define wsCalloc(R_p, t, sz)			R_p = (t *)(MEM().alloc(#R_p, __FILE__, __LINE__, sizeof(t)*(sz), 1))
#	define wsAlloc(R_p, t, sz)			R_p = (t *)(MEM().alloc(#R_p, __FILE__, __LINE__, sizeof(t)*(sz)))
#	define MULTI_MALLOC_SYM(R_p, t, sz)		R_p = (t *)(MEM().alloc(#R_p, S_file, N_line, sizeof(t)*(sz)))
#	define MULTI_ALIGNED_MALLOC(R_p, t, sz)	R_p = (t *)(MEM().alloc16(#R_p, __FILE__, __LINE__, sizeof(t)*(sz)))
#	define MULTI_ALIGNED_CALLOC(R_p, t, sz)	R_p = (t *)(MEM().alloc16(#R_p, __FILE__, __LINE__, sizeof(t)*(sz), 1))
#	define MULTI_ALN_MALLOC_SYM(R_p, t, sz)	R_p = (t *)(MEM().alloc16(#R_p, S_file, N_line, sizeof(t)*(sz)))
#	define MULTI_ALN_CALLOC_SYM(R_p, t, sz)	R_p = (t *)(MEM().alloc16(#R_p, S_file, N_line, sizeof(t)*(sz), 1))
#	define DEALLOC(R_p)						MEM().memfree(__FILE__, __LINE__, R_p)
#	define DEALLOC_ALIGNED(R_p)				MEM().memfree(__FILE__, __LINE__, R_p)
#else
#	ifndef DEALLOC
#		define wsCalloc					CALLOC
#		define wsAlloc					MALLOC
#		define MULTI_MALLOC_SYM				MALLOC
#		define MULTI_ALIGNED_MALLOC			MALLOC_ALIGNED
#		define MULTI_ALIGNED_CALLOC			CALLOC_ALIGNED
#		define MULTI_ALN_MALLOC_SYM			MALLOC_ALIGNED
#		define MULTI_ALN_CALLOC_SYM			CALLOC_ALIGNED
#		define DEALLOC						MFREE
#		define DEALLOC_ALIGNED				MFREE_ALIGNED
#	endif
#endif

/* Memmgr-related */
#ifdef USE_MEMDBG
	class cMemManager;
	cMemManager& MEM();
#endif

class cAnalysis;
class cIO;
cAnalysis*	getCorrelation(cIO *io);

/* Baseline memory allocation function that wrapped from various macros */
#define CALLOC(R_p, t, sz)			try { R_p = new t[sz]; } \
	catch (exception &e) { \
		halt_fmt(WISARD_FAIL_MEMALLOC_W_SYSMSG, sizeof(t)*(sz), __FILE__, __LINE__, e.what()); \
	} memset(R_p, 0x00, sizeof(t)*(sz))
#define MALLOC(R_p, t, sz)			try { R_p = new t[sz]; } \
	catch (exception &e) { \
		halt_fmt(WISARD_FAIL_MEMALLOC_W_SYSMSG, sizeof(t)*(sz), __FILE__, __LINE__, e.what()); \
	}

#ifdef USE_AVX
#	define SSE_SZALN		32
#	define SSE_ALIGN(t,sz)	((((sizeof(t)*(sz))+31)>>5)<<5)
#else
#	define SSE_SZALN		16
#	define SSE_ALIGN(t,sz)	((((sizeof(t)*(sz))+15)>>4)<<4)
#endif

#ifdef __FreeBSD__
#	define MALLOC_ALIGNED(R_p, t, sz)	posix_memalign((void **)&(R_p), SSE_SZALN, sizeof(t)*(sz)); \
	if ((R_p) == NULL) { \
	halt_fmt(WISARD_FAIL_MEMALLOC, sizeof(t)*(sz), __FILE__, __LINE__); \
	}
#	define CALLOC_ALIGNED(R_p, t, sz)	posix_memalign((void **)&(R_p), SSE_SZALN, sizeof(t)*(sz)); \
	if ((R_p) == NULL) { \
	halt_fmt(WISARD_FAIL_MEMALLOC, sizeof(t)*(sz), __FILE__, __LINE__); \
	} memset(R_p, 0x00, sizeof(t)*(sz))
#else
#	define MALLOC_ALIGNED(R_p, t, sz)	if ((R_p = (t *)_aligned_malloc( \
	SSE_ALIGN(t,sz), SSE_SZALN)) == NULL) { \
	halt_fmt(WISARD_FAIL_MEMALLOC, sizeof(t)*(sz), __FILE__, __LINE__); \
	}
#	define CALLOC_ALIGNED(R_p, t, sz)	if ((R_p = (t *)_aligned_malloc( \
	SSE_ALIGN(t,sz), SSE_SZALN)) == NULL) { \
	halt_fmt(WISARD_FAIL_MEMALLOC, sizeof(t)*(sz), __FILE__, __LINE__); \
	} memset(R_p, 0x00, sizeof(t)*(sz))
#endif

#define MFREE(R_p)					if (R_p) { delete [] (R_p); R_p=NULL; }
#define MFREE_ALIGNED(R_p)			if (R_p) { _aligned_free(R_p); R_p=NULL; }

#ifdef USE_SSE
#	define sseMalloc		MULTI_ALIGNED_MALLOC
#	define sseCalloc		MULTI_ALIGNED_CALLOC
#	define sseMallocSym		MULTI_ALN_MALLOC_SYM
#	define sseCallocSym		MULTI_ALN_CALLOC_SYM
#	define sseFree			DEALLOC_ALIGNED
#	define SSE_LOOP(i,to)	for (wsUint i=0 ; i<to ; i+=sseJmp)	
/* sseInit */
#	define sseInit(p, s, v) { \
		sse_t vvv=sseSet(v); \
		wsUint vv=getMed(s); \
		for (wsUint q=0 ; q<vv ; q+=sseJmp) { \
			*((sse_t *)(p+q)) = vvv; \
		} \
		for (wsUint q=vv ; q<s ; q++) p[q] = (wsReal)(v); \
	}
#else
#	define sseMalloc			MULTI_MALLOC
#	define sseCalloc			MULTI_CALLOC
#	define sseFree				DEALLOC
#	define sseInit(p, s, v) { \
		for (wsUint q=0 ; q<s ; q++) p[q] = REAL_CONST(v); \
	}
#	define SSE_LOOP(i,to)		for (wsUint i=0 ; i<to ; i++)
#endif

/* nucleotide char?? ?????? ????????? */
#define BASEtoBIN(s) ((s)&0x04 ? ((s)&0x01 ? 2 : 3) : ((s)&0x02 ? 1 : 0))
/* Is the specified genotype is missing? */
#define isMissing(Ra_mat)		((Ra_mat) == WISARD_NA)
#define isAvailable(Ra_mat)		((Ra_mat) != WISARD_NA)
#define isMissingReal(Ra_mat)	((Ra_mat) == WISARD_NA_REAL)
#define isAvailableReal(Ra_mat)	((Ra_mat) != WISARD_NA_REAL)

#define mapIsInexist(obj, elem)	(obj).find(elem) == (obj).end()
#define mapIsExist(obj, elem)	(obj).find(elem) != (obj).end()

template<class T>
inline const T SIGN(const T &a, const T &b)
{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}
template<class T>
inline const T MAX(const T &a, const T &b)
{return b > a ? (b) : (a);}
template<class T>
inline const T MIN(const T &a, const T &b)
{return b < a ? (b) : (a);}

#ifndef TRUE
	#define TRUE 1
#endif

/* Launch R with specific file */
int			launchR(const char *S_Rfn);
void		float_fmt(wsReal v, char *p);

/* ????? seq?��??? [start, end] ?? ????? ?????? ret?? ??????
 * ??? ret?? ?????? ?????? ????? ??? ???? ??
 */
int getPartialSeq(char *S_srcSeq, int N_startOffset, int N_endOffset, char *Sp_seqRet);

/* LOG2???? ???? */
double	LOG2(double R_number);

/* Similar as strtok but more complicated */
void	getString(char **Sp_str, char **Sp_nextStr);
void	getString(char **Sp_str, char **Sp_nextStr, char S_delim);

#ifndef USE_DBL
#error "DONOT"
#endif

#ifdef USE_DBL
void	deallocMatrix(wsFloat **Rp_mat, wsUint N_row, void *Vp_supp=NULL);
#endif
void	deallocMatrix(wsUmat Np_mat, wsUint N_row, void *Vp_supp=NULL);
void	deallocMatrix(wsMat Rp_mat, wsUint N_row, void *Vp_supp=NULL);
void	deallocMatrix(USHORT_t **Np_mat, wsUint N_row, void *Vp_supp=NULL);
wsMat	makeMatrix(wsStrCst S_path, wsUint *Np_r=NULL, wsUint *Np_c=NULL);
wsMat	makeMatrixSSE(wsStrCst S_path, wsUint *Np_r=NULL, wsUint *Np_c=NULL);
wsMat	makeSymMatSSE(wsStrCst S_path, wsUint *Np_r=NULL);
void	compareMatrix(wsMat Rp_mat, wsUint N_r, wsUint N_c, wsStrCst S_path,
	wsStrCst S_prefix=NULL);
void	compareMatrix(wsMat Rp_mat1, wsUint N_r, wsUint N_c, wsMat Rp_mat2,
	wsStrCst S_path);
wsMat	allocateMatrix(wsUint N_row, wsUint N_col);

/*
 * Structure definition
 */

typedef	map<string,int>				mDataIdx;
typedef	mDataIdx::iterator			mDataIdx_it;
typedef unordered_set<string>		eStr;
typedef	eStr::iterator				eStr_it;
typedef map<string,wsReal>			mStrReal;
typedef mStrReal::iterator			mStrReal_it;
typedef vector<mStrReal>			vmStrReal;
typedef vmStrReal::iterator			vmStrReal_it;
typedef vmStrReal::reverse_iterator	vmStrReal_rit;
typedef map<string,vInt>			mvInt;
typedef mvInt::iterator				mvInt_it;
typedef map<string,vReal>			mvReal;
typedef mvReal::iterator			mvReal_it;
typedef map<wsUint, int>			mInt;
typedef mInt::iterator				mInt_it;


double pythag(const double a, const double b);
double util_log1p(double x);

void	LRtest(wsReal R_v1, wsReal R_v2, wsReal *Rp_lrt,
	wsReal *Rp_pValue, wsReal R_df=W1);

bool	SVDinverse(wsReal **Ra_mat, wsUint N_sz, wsReal ***Rp_res,
	wsUint *Np_rank=NULL);
bool	EDinverse(wsReal **Ra_mat, wsUint N_sz, wsReal ***Rp_res,
	wsUint *Np_rank=NULL);
wsMat	eigenInvSymMatrix(wsMat Ra_eVec, wsReal *Ra_eVal, wsUint N_sz);

wsMat	transpose(wsReal **Ra_mat, wsUint N_matR, wsUint N_matC);
void	transposeSelf(wsReal **Ra_mat, wsUint N_matR, wsUint N_matC);

int		tqli1(int n, wsReal *d, wsReal *e);
char	tqli(wsReal *d, wsReal *e, int n, wsReal **z);
char	tqli_dbl(wsReal *d, wsReal *e, int n, double **z);
void	tred2(wsReal **a, int n, wsReal *d, wsReal *e);
void	tred2_dbl(double **a, int n, wsReal *d, wsReal *e);
int		jcbi(wsReal **Ra_mat, wsReal **Ra_eV, int N_elem, wsReal eps, int jt);
wsReal*	EVpower(double **Ra_mat, wsUint N_sz, wsUint N_pc, double R_thr);

wsReal*	eigen(wsReal **Ra_mat, wsUint N_sz, char *Bp_stat);
wsReal*	eigenDecomp(wsReal **Ra_mat, wsUint N_sz, char *Bp_stat,
	wsReal ***Ra_eVec, char B_transpose=0, char B_negCorrect=1);
wsReal*	eigenDecomp2(wsReal **Ra_mat, wsUint N_sz, char *Bp_stat,
	wsReal ***Ra_eVec, char B_noTranspose=0, char B_negCorrect=1);
wsReal* EIGENDECOMPOSITION(wsSym Ra_m, int N_sz, char *Bp_isFailed,
	wsMat *Ra_evec=NULL, char B_noTranspose=0, char B_noSort=0);	
wsReal* EIGENDECOMPOSITION(wsSym Ra_m, int N_sz, wsMat *Ra_evec=NULL,
	char B_noTranspose=0, char B_noSort=0);	
wsMat	sqrtMatrix(wsReal **Ra_mat, wsUint N_sz, char B_sym=0);
wsMat	sqrtMatrixEV(wsMat Ra_eVec, wsVec Ra_eVal, wsUint N_sz,
	char B_sym=0);

wsReal	detCholesky(wsReal **Ra_m, wsUint N_sz);
wsReal	determinantCholesky(wsReal **Ra_m, wsUint N_sz);
void	chol(wsMat Ra_m, wsUint N_sz);

wsReal*	EVlanczos(wsReal **Ra_mat, wsUint N_sz, wsUint N_iter,
				wsUint *Np_getE,
				wsUint N_szE=0xffffffff, wsReal R_eps=5e-15);

wsUint*	kmeans(wsReal *Ra_val, wsUintCst N_szVal, wsUintCst N_clus);
wsUint*	kmeans(wsReal **Ra_val, wsUintCst N_dim, wsUintCst N_szVal, wsUintCst N_clus);

double	Brent_fmin(double ax, double bx, double (*f)(double, void *),
	void *info, double tol);

bool	SingularValueDecomposition(wsMat a, int n, wsReal *w, wsReal **v);
bool	SingularValueDecomposition(wsVec a, int n, wsReal *d, wsReal *v);
int		dsvd(wsVec a, int m, wsReal* w);

char**	loadStringValues(wsStrCst Sp_val, wsUint *Np_var);
char**	loadStringValues2(wsStrCst Sp_val, wsUint *Np_var, char S_sep);
char**	loadStringValuesByWS(wsStrCst Sp_val, wsUint *Np_var);
void	loadGroupValues(wsStrCst Sp_val, mStrReal &Xm_map);
void	loadFlags(wsStrCst Sp_val, mStr &Xm_map);

/**
 * loadRealValues extracts a set of real values from string
 *
 * @param     Sp_val		A pointer to the string includes real values to be extracted
 * @param     Np_var		A pointer to the integer that the number of extracted real values be stored
 * @return    (wsReal*)		A pointer to the array of real values extracted
 */
wsVec	loadRealValues(wsStrCst Sp_val, wsUint *Np_var);
wsVec	loadRealValuesByWS(wsStrCst Sp_val, wsUint *Np_var);

/**
 * loadIntValues extracts a set of integer values from string
 *
 * @param     Sp_val		A pointer to the string includes integer values to be extracted
 * @param     Np_var		A pointer to the integer that the number of extracted integer values be stored
 * @return    (wsUint*)		A pointer to the array of integer values extracted
 */
wsUint*	loadIntValues(char *Sp_val, wsUint *Np_var);
wsUint*	loadIntValuesByWS(char *Sp_val, wsUint *Np_var);

void	sseLDLt(wsReal **Ra_mLowTri, wsUint N_sz, wsReal *Ra_diag,
	 wsReal **Ra_dest);

/**
 * sseVP	Performs VP multiplication
 *
 * @param	Ra_v	Former vector (V) [n*1] to be multiplicated
 * @param	N_sz	size of V
 * @param	Ra_p	Value vector of sparse matrix (P) [n*q] to be multiplicated
 * @param	N_c		# column of P
 * @param	Na_rEnd	Offset vector of row [size N_pr] end
 * @param	Na_cIdx	Index vector of column [size Na_rEnd[END]] for items
 * @param	Ra_mask	(optional) Computation mask [size N_sz],
 *					0 then skip, otherwise compute
 * @return	(VEC_t)	Return of VP [q*1] = [1*n] %*% [n*q]
 */
wsVec	sseVpP(wsVec Ra_v, wsUint N_sz, wsReal *Ra_p, wsUint N_c,
	wsUint *Na_rEnd, wsUint *Na_cIdx, wsReal *Ra_mask=NULL);

/**
 * sseVpPpV	Performs VPW multiplication, return will be scalar
 *
 * @param	Ra_V		Former vector (V) [n*1] to be multiplicated
 * @param	N_sz		size of V
 * @param	Ra_p		Value vector of sparse matrix (P) [n*q] to be multiplicated
 * @param	N_c			# column of P
 * @param	Na_rEnd		Offset vector of row [size N_pr] end
 * @param	Na_cIdx		Index vector of column [size Na_rEnd[END]] for items
 * @param	Ra_W		(optional) Latter vector (W) [q*1], V if not assigned
 * @param	Ra_mask		(optional) Computation mask [size N_sz],
 *						0 then skip, otherwise compute
 * @return	(REAL_c)	Return of VPW [1] = [1*n] %*% [n*q] %*% [q*1]
 */
wsRealCst	sseVpPpV(wsVec Ra_V, wsUint N_sz, wsReal *Ra_p, wsUint N_c,
	wsUint *Na_rEnd, wsUint *Na_cIdx, wsVec Ra_W=NULL, wsReal *Ra_mask=NULL);

// wsReal**	sseSD(SYM_t Ra_symMat, wsUint N_s1, DIAG_t Ra_diagMat);
// wsReal**	sseDS(DIAG_t Ra_diagMat, wsUint N_s1, SYM_t Ra_symMat);

wsReal	distEuclidean(wsRealCst *Ra_1, wsRealCst *Ra_2, wsUint N_sz);

wsSym	sseSymSbsMatrix(wsMat Rp_inpMatrix, wsUint N_sz, wsUint *Na_newSeq,
	wsUint N_newSz=0xffffffff);
wsReal*	sseRowMatrix(wsMat Rp_mat, wsUint N_sz, wsUint N_idx,
	wsUint N_rem=0xffffffff);

/* Frobenius norm calculation */
wsReal	normFrobenius(wsMat Ra_mat, wsUint N_r, wsUint N_c);
wsReal	normFrobenius(wsSym Ra_mat, wsUint N_sz);

// Ra_m is difference itself
wsReal compMMfrobenius(wsMat Ra_m1, wsUint N_r1, wsUint N_c1, wsMat Ra_m2,
	wsUint N_r2, wsUint N_c2);
wsReal compSMfrobenius(wsSym Ra_s1, wsUint N_sz, wsMat Ra_m2, wsUint N_r2,
	wsUint N_c2);
wsReal compSSfrobenius(wsSym Ra_s1, wsUint N_s1, wsSym Ra_s2, wsUint N_s2);


/* Do the boundaries exactly for q*() functions :
 * Often  _LEFT_ = ML_NEGINF , and very often _RIGHT_ = ML_POSINF;
 *
 * R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)  :<==>
 *
 *     R_Q_P01_check(p);
 *     if (p == R_DT_0) return _LEFT_ ;
 *     if (p == R_DT_1) return _RIGHT_;
 *
 * the following implementation should be more efficient (less tests):
 */
#define R_Q_P01_boundaries(p, _LEFT_, _RIGHT_)		\
    if (log_p) {					\
	if(p > 0)					\
	    return numeric_limits<double>::quiet_NaN();				\
	if(p == 0) /* upper bound*/			\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
	if(p == ML_NEGINF)				\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
    }							\
    else { /* !log_p */					\
	if(p < 0 || p > 1)				\
	    return numeric_limits<double>::quiet_NaN();				\
	if(p == 0)					\
	    return lower_tail ? _LEFT_ : _RIGHT_;	\
	if(p == 1)					\
	    return lower_tail ? _RIGHT_ : _LEFT_;	\
    }
#define R_D_Lval(p)		(lower_tail ? (p) : (0.5 - (p) + 0.5))	/*  p  */
#define R_DT_qIv(p)		(log_p ? (lower_tail ? exp(p) : - util_expm1(p)) \
							: R_D_Lval(p))
#define R_D_Clog(p)		(log_p	? util_log1p(-(p)) : (0.5 - (p) + 0.5)) /* [log](1-p) */
#define R_D_val(x)		(log_p	? log(x) : (x))		/*  x  in pF(x,..) */
#define R_DT_Cval(x)	(lower_tail ? R_D_Clog(x) : R_D_val(x))
#define R_D_Cval(p)		(lower_tail ? (0.5 - (p) + 0.5) : (p))	/*  1 - p */
#define R_DT_CIv(p)		(log_p ? (lower_tail ? -util_expm1(p) : exp(p)) \
							: R_D_Cval(p))
#define P1P(v) ((v)/(W1-(v))) // P / (1-P)

int		imin2(int x, int y);
double	fmin2(double x, double y);
double	fmax2(double x, double y);
double	qbeta(double alpha, double p, double q, int lower_tail, int log_p);
double	util_expm1(double x);

/* liu c
SKAT_liu <- function(q, lambda, h = rep(1,length(lambda)), delta = rep(0,length(lambda))) {
*/
double	liuEV(double q, double *Ra_lambda, wsUint r);
double	chebyshev_eval(double x, const double *a, const int n);
double	lgammafn(double x);
double	genBeta(double aa, double bb);
double	genBinomial(double nin, double pp);

/* Split string by white-spaces
 *
 * ------------------+----------------------------------------------------
 * char	**Sp_str	 | A pointer that points to the beginning of string
 * char **Sp_nextStr | A pointer that be pointed to the beginning of
 *                   |   next string
 * ------------------+---------------------------------------------------------
 *  return			 | Nothing
 */
#ifdef _DEBUG
#	define getString(a, b) _getString((a), (b), __LONGFUNC__)
#	define getStringDelim(a, b, c) _getStringDelim((a), (b), (c), __LONGFUNC__)
#else
#	define getString(a, b) _getString((a), (b))
#	define getStringDelim(a, b, c) _getStringDelim((a), (b), (c))
#endif
//#define getString(a, b, c) _getString((a), (b), (c), __LONGFUNC__)
inline void _getString(char **Sp_str, char **Sp_nextStr, wsStrCst S_fname=NULL)
{
	char *p = *Sp_str;
	if (p == NULL) halt("Null pointer given to make string in function [%s]", S_fname?S_fname:"");

	/* Skip whitespace */
	for ( ; *p == '\t' || *p == '\r' || *p == '\n' || *p == ' ' ; p++);
	*Sp_str = p;
	for ( ; *p != '\t' && *p != '\r' && *p != '\n' && *p != ' ' && *p != '\0' ; p++);

	/* If Sp_delim */
	if (*p != '\t' && *p != ' ') {
		*p = '\0';
		*Sp_nextStr = NULL;
	} else {
		*p = '\0';
		/* Skip whitespace */
		for (p++ ; *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' ; p++);
		if (*p == '\0')
			*Sp_nextStr = NULL;
		else
			*Sp_nextStr = p;
	}
}

inline void _getStringDelim(char **Sp_str, char **Sp_nextStr, char S_delim,
	wsStrCst S_fname=NULL)
{
	char *p = *Sp_str;
	if (p == NULL) halt("Null pointer given to make string in function [%s]",
		S_fname?S_fname:"");

	for ( ; *p != '\t' && *p != '\r' && *p != '\n' && *p != ' ' && *p != '\0' && *p != S_delim ; p++);

	/* If Sp_delim */
	if (*p != '\t' && *p != ' ' && *p != S_delim) {
		*p = '\0';
		*Sp_nextStr = NULL;
	} else {
		*p = '\0';
		/* Skip whitespace */
		for (p++ ; *p == ' ' || *p == '\t' || *p == '\r' || *p == '\n' ; p++);
		if (*p == '\0')
			*Sp_nextStr = NULL;
		else
			*Sp_nextStr = p;
	}
}

char getItselfPath(char *Sp_path);

/* ??????/?????? ?????? ???? */
#ifdef _WIN32
// epoch time???? ????? ???
#	if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#		define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#	else
#		define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

	// for timezone
	struct timezone
	{
		int  tz_minuteswest; /* minutes W of Greenwich */
		int  tz_dsttime;     /* X_logType of dst correction */
	};

	// gettimeofday in windows
	int gettimeofday(struct timeval *tv, struct timezone *tz);
#endif

inline int findVstr(vStr &Xv_str, string S_val)
{
	int k = 0, ret = -1;
	FOREACHDO (vStr_it, Xv_str, K, k++)
		if (S_val.compare(*K) == 0) {
			ret = k;
			break;
		}
		return ret;
}

void multcomp(cIO *Cp_IO, vector<double>& chi, string title);
void multcomp(cIO *Cp_IO, double *Ra_stat, wsUint N_stat, string title);

char fetchMarkerList(mDataIdx* &Xm_mmap, wsStrCst S_fn, vInt &Xv_ret,
	wsStrCst S_desc=NULL);
char fetchMarkerList(wsStrCst S_fn, vStr &Xv_ret, wsStrCst S_desc=NULL);
char fetchMarkerList(wsStrCst S_fn, eStr &Xe_ret, wsStrCst S_desc=NULL);
void getEigenPhi(wsSym Ra_phi, wsUint N_anaSamp, wsMat *Rp_eVec,
	wsReal **Rp_eVal, wsUint *Na_testFamIdx=NULL);

#define MAX_TIES_CNT				65536		/* ranksum() ??????? ??? ??? ???? */

/* Prob. of multinormal dist */
double pmnorm(wsReal *Ra_X, wsUint N_sz, wsSym Ra_vcov, wsReal *Ra_mean=NULL);
double pmnorm(wsReal R_X, wsUint N_sz, wsSym Ra_cor);

double lbeta(double a, double b);

/* time */
#define WISARD_NOW(r) struct tm r; { time_t t = ::time(NULL); r = *localtime(&t); }

double ltqnorm(double q);
double FET_HY(int N_11i, int N_1Xi, int N_X1i, int N_XXi);

typedef struct _xRange {
	wsUint N_start;
	wsUint N_end;
} xRange;

typedef vector<xRange> vRange;
typedef vRange::iterator	vRange_it;

/* Satterthwaite DF function */
wsReal satterthwaite(wsVecCst Ra_a, wsUint N_len, wsUint N_df, wsReal* Rp_df);

wsRealCst pchisqsum(double x, double df, wsVec a, wsUint N_sz, char lower_tail=TRUE,
	char method=0);

double ws_str2dbl(const char *S_str, char **Sp_end);

#define WSSTR_INITSIZE		1024
#define WSSTR_EXPANDSIZE	1024
class wsString
{
	size_t	N_len;
	size_t	N_buf;
	char*	S_buf;
	char*	S_cur;
public:
	wsString() {
		S_buf = NULL;
		S_cur = NULL;
		N_len = N_buf = 0;
	}
	~wsString() {
		DEALLOC(S_buf);
	}
	void rem() {
		DEALLOC(S_buf);
		S_buf = NULL;
		S_cur = NULL;
		N_len = N_buf = 0;
	}
	wsStrCst get() {
		return S_buf?S_buf:"";
	}
	void request(size_t N_req) {
		size_t N_remain = 0;

		/* Have buffer? */
		if (N_buf)
			/* Enough space left? */
				N_remain = N_buf - N_len - 1;
		/* Need space? */
		if (N_remain < N_req) {
			size_t N_new = (size_t)((N_req - N_remain + WSSTR_EXPANDSIZE - 1) / (wsReal)WSSTR_EXPANDSIZE)
				* WSSTR_EXPANDSIZE + N_buf;
			char* S_newBuf = NULL;
			wsAlloc(S_newBuf, char, N_new);
			if (N_buf) {
				strcpy(S_newBuf, S_buf);
				DEALLOC(S_buf);
			} else
				S_newBuf[0] = '\0';
			S_buf = S_newBuf;
			S_cur = S_buf + N_len;
			N_buf = N_new;
		}
	}
	wsString(wsStrCst S_str) {
		rem();
		size_t N_str = strlen(S_str);
		request(N_len);
		strcpy(S_cur, S_str);
		S_cur += N_str;
		N_len += N_str;
	}
	wsString& operator+=(string& S_str) {
		size_t N_sz = S_str.size();
		request(N_sz);
		strcpy(S_cur, S_str.c_str());
		S_cur += N_sz;
		N_len += N_sz;
		return *this;
	}
	wsString& operator+=(wsUint N_val) {
		request(11); /* Maximum possible size */
		sprintf(S_cur, "%d", N_val);
		wsUint N_sz = (wsUint)strlen(S_cur);
		S_cur += N_sz;
		N_len += N_sz; 
		return *this;
	}
	wsString& append(char S_chr, wsUint N_times) {
		request(N_times); /* Maximum possible size */
		while (N_times--)
			*(S_cur++) = S_chr;
		*S_cur = '\0';
		N_len += N_times;
		return *this;
	}
	wsString& operator+=(char S_chr) {
		request(1); /* Maximum possible size */
		*(S_cur++) = S_chr;
		*S_cur = '\0';
		N_len++;
		return *this;
	}
	wsString& operator+=(wsStrCst S_str) {
		size_t N_str = strlen(S_str);
		request(N_len);
		strcpy(S_cur, S_str);
		S_cur += N_str;
		N_len += N_str;
		return *this;
	}
	wsString& operator+=(wsReal R_val) {
		request(32); /* Maximum possible size */
		sprintf(S_cur, "%g", R_val);
		size_t N_sz = strlen(S_cur);;
		S_cur += N_sz;
		N_len += N_sz;
		return *this;
	}
};

wsReal randInRange(wsReal R_min, wsReal R_max);

typedef enum _xRankTie {
	RT_MIN,
	RT_AVERAGE,
	RT_MAX
} xRankTie;
wsVec sseVrank(wsVecCst Ra_arr, wsUint N_sz, xRankTie X_tie=RT_AVERAGE);
wsVec sseVarank(wsVecCst Ra_arr, wsUint N_sz, xRankTie X_tie=RT_AVERAGE);

void getMachineName(char* Sp_mname);

// Benjamini-Hochberg Q value
void Qbenjhoch(xRealSort* Xa_rs, wsUintCst N_sz, wsVec Ra_ret);

} // End namespace ONETOOL

#endif
