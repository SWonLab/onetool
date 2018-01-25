#include <limits>
#include <algorithm>
#include "utils/cdflib.h"
#include "analyses/gstest.h"
#include "utils/data.h"
#include "analyses/emai.h"
#include "analyses/setmgr.h"
#include "analyses/ppp.h"
#include "analyses/regr.h"
#include "utils/integrate.h"
#include "utils/dcdflib.h"
#include "utils/matrix.h"
#include "analyses/gemma.h"
#include "analyses/fqls.h"
#include "utils/stat.h"
#define MIN_PERM_VT 200

//#define GSvalidate

namespace ONETOOL {

/*
 *
 * DEFINITIONS
 *
 */
wsRealCst _getPEDCMC(cIO* Cp_IO, wsUintCst N_samp, cMatrix* Mp_phi, cVector* Vp_yGSunAdj,
	const char *S_gs, vInt &Xa_set,
	wsReal **Ra_PEDCMCca, wsReal **Ra_PEDCMCct, wsReal *Rp_stat,
	wsUint N_ctrl, wsUint N_case, wsUint N_common, wsUint N_rare,
	cMask &V_mask, char *Ba_filt, wsUint *Np_typeTest, wsUint *Np_rank);
xLiuCov** _getSKATo(wsUintCst N_phe,
	const char *S_gs, vInt &Xa_set,
	cStdMatrix&	M_U,		/* Checked = GW */
	cVector&	V_wGScurr,
	//	VEC_t		Ra_wGScurr,	/* Checked */
	cStdMatrix&	M_genes,
	wsVec		Ra_xxG,			/* Checked */
	wsUintCst	N_curSamp,
	cVector&	V_ySKATo,
	wsReal*		Rp_rhoSKATO,
	wsReal*		Rp_Popt,
	cSymMatrix*	Mp_covs,
	xGeneTestParam& X);
wsRealCst _getCollapsing(cMask &V_mask, wsRealCst R_vars, cVector& V_aV,
	cVector &V_pldClp,
	wsReal *Rp_statClp, wsUint N_samp, wsReal *Rp_aVy);
wsRealCst _getCMC(cMask &V_mask, cVector& V_hat, wsRealCst R_vars,
	cVector& V_aV,
	cVector &V_pldCMC, wsReal *Rp_statCMC, wsUintCst N_samp);
wsUint _getSKAT(wsRealCst R_pppMult, const char *S_gs, vInt &Xa_set,
	cStdMatrix& M_adjGeno, wsRealCst R_stat, wsUintCst N_actualSamp, wsReal *Rp_pVal);

/* Running only once
 * Divide SNPs ALMOST equally to each thread */
int forGene_equal(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	int			*Np_idx	= (int *)Vp_thrData;
	xAnaThread	*Cp_ana	= (xAnaThread *)Vp_shareData;
	cSetManagerAnalysis*
				Cp_gs	= (cSetManagerAnalysis *)(Cp_ana->Vp_data);
	mGeneDef&	Xm_gene	= Cp_gs->getGeneDef();
	wsUintCst		N_gene	= (wsUint)(Xm_gene.size());
	static wsUint
				N_proc	= 0;
	wsReal		R_szDiv	= (wsReal)N_gene/(wsReal)N_thread;
	int			N_ret	= N_proc != 0 ? 0 : 1;
	int			i=0;
	wsReal		j=W0;

	switch (N_mode) {
	case DISTFUN_INIT:
		i = 0;
		j = W0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		if (N_ret == 0)
			return 0;
		for ( ; i<N_thread ; i++,j+=R_szDiv) {
			wsUint E = (wsUint)(j+R_szDiv+REAL_CONST(0.5));
			if (E > N_gene)
				E = N_gene;
			Np_idx[3*i] = (wsUint)(j+REAL_CONST(0.5));
			Np_idx[3*i+1] = E;
			Np_idx[3*i+2] = 0;
		}
		N_proc = N_ret = N_thread;
		break;
	case DISTFUN_UNINIT:
		LOG("%d/%d variants processed\n", N_gene, N_gene);
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}

	return N_ret;
}

inline float   abs1(float x)    {return (x>0.0 ? x :-x);}
inline int      abs1(int x)       {return (x>0 ? x :-x);}

typedef enum _xRankerMethod
{
	RM_AVERAGE,
	RM_MIN,
	RM_MAX,
	RM_DEFAULT
} xRankerMethod;

template <class T>
class lt { public: static int compare(T a, T b) { return(a < b); } };
template <class T>
class gt { public: static int compare(T a, T b) { return(a > b); } };

template <class T, class C>
class ranker
{
private:
	const T* p;
	size_t sz;

public:
	ranker(const vector<T>& v) : p(&v[0]), sz(v.size()) { }
	ranker(const T* tp, size_t s) : p(tp), sz(s) { }

	int operator()(size_t i1, size_t i2) const { return(C::compare(p[i1],p[i2])); }

	template <class S>
	void get_orders(vector<S>& w) const {
		w.resize(sz);
		w.front() = 0;
		for (typename vector<S>::iterator i = w.begin(); i != w.end() - 1; ++i)
			*(i + 1) = *i + 1;
		std::sort(w.begin(), w.end(), *this);
	}

	template <class S>
	void get_partial_orders(vector<S>& w, size_t num) const {
		if (num > sz) num = sz;
		w.resize(sz);
		w.front() = 0;
		for (typename vector<S>::iterator i = w.begin(); i != w.end() - 1; ++i)
			*(i + 1) = *i + 1;
		std::partial_sort(w.begin(), w.begin() + num, w.end(), *this);
		w.resize(num);
	}

	template <class S>
	void get_ranks(vector<S>& w, const xRankerMethod method) const {
		w.resize(sz);
		vInt tmp(w.size());
		get_orders(tmp);
		switch (method) {
		case RM_AVERAGE:
			for (size_t c=0,reps ; c<w.size(); c+=reps) { reps = 1;
			while (c + reps < w.size() && p[tmp[c]] == p[tmp[c + reps]]) ++reps ;
			for (size_t k=0 ; k<reps ; ++k)
				w[tmp[c + k]] = S(2 * c + reps - 1) / 2 + 1;
			} break;
		case RM_MIN:
			for (size_t c=0,reps ; c<w.size(); c+=reps) { reps = 1;
			while (c + reps < w.size() && p[tmp[c]] == p[tmp[c + reps]]) ++reps ;
			for (size_t k=0 ; k<reps ; ++k) w[tmp[c + k]] = S(c + 1);
			} break;
		case RM_MAX:
			for (size_t c=0,reps ; c<w.size(); c+=reps) { reps = 1;
			while (c + reps < w.size() && p[tmp[c]] == p[tmp[c + reps]]) ++reps ;
			for (size_t k=0 ; k<reps ; ++k) w[tmp[c + k]] = S(c + reps);
			} break;
		case RM_DEFAULT: {
			// default
			// for (uint c=0; c<num ; ++c) w[tmp[c]] = c + 1;
			wsUint counter = 0;
			for (size_t c=0,reps ; c<w.size(); c+=reps) { reps = 1; 
			while (c + reps < w.size() && p[tmp[c]] == p[tmp[c + reps]]) ++reps ;
			for (size_t k=0 ; k<reps ; ++k) w[tmp[c + k]] = counter;
			counter++;
			}
						 } break;
		}
	}

	template <class S>
	void get_partial_ranks(vector<S>& w, const xRankerMethod method, size_t num) const {
		if (num > sz) num = sz;
		vInt tmp(sz);
		get_partial_orders(tmp, num);
		w.resize(sz);
		fill(w.begin(), w.end(), 0);
		switch (method) {
		case RM_AVERAGE:
			for (size_t c=0,reps ; c<num ; c+=reps) { reps = 1;
			while (c + reps < num && p[tmp[c]] == p[tmp[c + reps]]) ++reps ;
			for (size_t k=0 ; k<reps ; ++k)
				w[tmp[c + k]] = S(2 * c + reps - 1) / 2 + 1;
			} break;
		case RM_MIN:
			for (size_t c=0,reps ; c<num ; c+=reps) { reps = 1;
			while (c + reps < num && p[tmp[c]] == p[tmp[c + reps]]) ++reps ;
			for (size_t k=0 ; k<reps ; ++k) w[tmp[c + k]] = c + 1;
			} break;
		case RM_MAX:
			for (size_t c=0,reps ; c<num ; c+=reps) { reps = 1;
			while (c + reps < num && p[tmp[c]] == p[tmp[c + reps]]) ++reps ;
			for (size_t k=0 ; k<reps ; ++k) w[tmp[c + k]] = c + reps;
			} break;
		case RM_DEFAULT: {
			// default
			// for (wsUint c=0; c<num ; ++c) w[tmp[c]] = c + 1;
			wsUint counter = 0;
			for (size_t c=0,reps ; c<w.size() ; c+=reps) { reps = 1;
			while (c + reps < w.size() && p[tmp[c]] == p[tmp[c + reps]]) ++reps;
			for (size_t k=0 ; k<reps ; ++k) w[tmp[c + k]] = counter;
			counter++;
			}
						 } break;
		}
	}
};

template <class T, class S>
inline void _rank(const vector<T>& v, vector<S>& w,
				  xRankerMethod method = RM_AVERAGE)
{ ranker<T, lt<T> > r(v); r.get_ranks(w, method); }


/* !!! FARVAT-specific analysis !!! */
#if ((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE & FUNCTION_FARVAT) || \
	(TOOLSET_TYPE & FUNCTION_MFARVAT) || (TOOLSET_TYPE & FUNCTION_FARVATX)

typedef struct _xGStest
{
	wsReal	R_statSKAT, R_statCMC, R_statClp;
	wsReal	R_pvSKAT, R_pvCMC, R_pvClp;
	wsReal	R_statPEDCMC, R_pvPEDCMC;
} xGStest;

typedef struct _xGStestInput
{
	mGeneDef_it	&X_it;
	cIO			*Cp_IO;
	wsUint		N_storeIdx;
	wsUint		N_testSamp;
} xGStestInput;

wsVec getOptimalRhos(wsUint& N_rho)
{
	wsVec	Ra_rhos	= NULL;
	char**	Sa_rhos	= NULL; 

	if (IS_ASSIGNED(skatondiv))
		N_rho = OPT_NUMBER(skatondiv);
	else if (IS_ASSIGNED(skatodivs))
		Sa_rhos = loadStringValues(OPT_STRING(skatodivs), &N_rho);
	else
		N_rho = 8;
	Ra_rhos = sseVector(N_rho);

	if (IS_ASSIGNED(skatondiv)) {
		wsReal R_div = W1 / (wsReal)(N_rho-1);
		wsReal R_val = W0;
		for (wsUint i=0 ; i<N_rho; i++,R_val+=R_div)
			Ra_rhos[i] = R_val;
	} else if (IS_ASSIGNED(skatodivs)) LOOP (i, N_rho) {
		/* Load rhos */
		char *Sp_buf = NULL;
		Ra_rhos[i] = str2dbl(Sa_rhos[i], &Sp_buf);
		if (Sp_buf)
			halt("[%d]th optimal weight [%s] is invalid", i+1,
			Sa_rhos[i]);
	} else {
		Ra_rhos[0] = 0;
		Ra_rhos[1] = 0.01;
		Ra_rhos[2] = 0.04;
		Ra_rhos[3] = 0.09;
		Ra_rhos[4] = 0.16;
		Ra_rhos[5] = 0.25;
		Ra_rhos[6] = 0.5;
		Ra_rhos[7] = 1;
	}
	return Ra_rhos;
}

wsRealCst quadfactor(wsSym kinmat, wsUint N_sz, char chrom, wsVec resid)
{
	wsReal sex_factor = W2;
	if (chrom) {
// 		sex_factor <- 2*(sex==2) %o% (sex==2) +
// 			male.dose^2*(sex==1) %o% (sex==1) + 
// 			sqrt(2)*male.dose*((sex==2) %o% (sex==1)+(sex==1) %o% (sex==2)) 
	}

//	genet.var <- 2*kinmat
//## if no inbreeding, d = 1 so genet.var = genet.cor, but
//## below allows for inbred subjects
	wsVec d = sseSdiag(kinmat, N_sz);
	sseVsqrt(d, N_sz);
	//sqrt(diag(genet.var))

	wsSym genet_cor = sseSymMat(N_sz);
	LOOP (i, N_sz)
		sseVdV(kinmat[i], d, genet_cor[i], i+1, sex_factor/d[i]);
	//	<- sex.factor * (genet.var / (d%o%d))

	wsReal R_ret = sseVpSpV(resid, genet_cor, N_sz);
	sseUnmat(genet_cor, N_sz);
	return R_ret;
//	c.factor <-  t(resid) %*% genet.cor %*% resid 
//	return(c.factor)
}

xLiuCov* liuCov(wsSym Ra_cov, wsUint N_sz, wsSym Ra_R,
	wsReal **Ra_ttt/*=NULL*/, wsReal **Ra_ttt2/*=NULL*/,
	wsReal R_muQ/*=WISARD_NAN*/)
{
	cTimer	x;
	xLiuCov	*R = NULL;
	char	B_deAllocTTT	= Ra_ttt == NULL;
	char	B_deAllocTTT2	= Ra_ttt2 == NULL;
	wsCalloc(R, xLiuCov, 1);

	if (Ra_ttt == NULL) {
		x.start();
		Ra_ttt = sseThrRowSS(Ra_R, N_sz, Ra_cov, N_sz);
		pverbose("ttt1 %s\n", x.getReadable());
	}

	if (R_muQ != R_muQ) {
		R->R_muQ = W0;
		for (wsUint i=0 ; i<N_sz ; i++)
			R->R_muQ += Ra_ttt[i][i];
	} else
		R->R_muQ = R_muQ;

	if (Ra_ttt2 == NULL) {
		x.start();
		Ra_ttt2 = sseThrRowMM(Ra_ttt, N_sz, N_sz, Ra_ttt, N_sz, N_sz);
		pverbose("ttt2 %s\n", x.getReadable());
	}
	for (wsUint i=0 ; i<N_sz ; i++)
		R->R_c2 += Ra_ttt2[i][i];

	/* Compute c3 and c4 */
	x.start();
	R->R_c4 = diagSum2matrix(Ra_ttt2, N_sz, N_sz, Ra_ttt2, N_sz, N_sz);
	R->R_c3 = diagSum2matrix(Ra_ttt2, N_sz, N_sz, Ra_ttt, N_sz, N_sz);

	/* Release memory */
	if (B_deAllocTTT2)	sseUnmat(Ra_ttt2, N_sz);
	if (B_deAllocTTT)	sseUnmat(Ra_ttt, N_sz);

	R->R_sigQ		= sqrt(W2*R->R_c2);
	wsReal R_s1		= W0;
	wsReal R_s2		= W0;
	if (R->R_c2 != W0) {
		R_s1		= R->R_c3 / (R->R_c2*sqrt(R->R_c2));
		R_s2		= R->R_c4 / SQR(R->R_c2);
	}
	wsReal R_ss1	= SQR(R_s1);

	wsReal R_a		= W0;
	if (R_s2 != W0)
		R_a			= R_ss1>R_s2 ? W1/(R_s1-sqrt(R_ss1-R_s2)) :
			W1 / sqrt(R_s2);
	R->R_delta		= R_ss1>R_s2 ? R_s1*CUBE(R_a) - SQR(R_a) : W0;
//	LOG("ss1 %g s2 %g\n", R_ss1, R_s2);
	R->R_l			= R_a*R_a - (R->R_delta+R->R_delta);
	R->R_muX		= R->R_l + R->R_delta;
	R->R_sigX		= sqrt(W2 * (R->R_l + W2*R->R_delta));
	pverbose("rest %s\n", x.getReadable());

	return R;
}

xLiuCov* liuCovSS(wsSym Ra_ttt, wsSym Ra_ttt2, wsUint N_sz, wsReal R_muQ=WISARD_NAN)
{
	cTimer	x;
	xLiuCov	*R = NULL;
	wsCalloc(R, xLiuCov, 1);

	if (NA(R_muQ))
		R->R_muQ = diagSum(Ra_ttt, N_sz);
	else
		R->R_muQ = R_muQ;
	R->R_c2 = diagSum(Ra_ttt2, N_sz);

	/* Compute c3 and c4 */
	x.start();
	R->R_c4 = sseTrSS(Ra_ttt2, N_sz, Ra_ttt2);
	R->R_c3 = sseTrSS(Ra_ttt2, N_sz, Ra_ttt);

	R->R_sigQ		= sqrt(W2*R->R_c2);
	wsReal R_s1		= W0;
	wsReal R_s2		= W0;
	if (R->R_c2 != W0) {
		R_s1		= R->R_c3 / (R->R_c2*sqrt(R->R_c2));
		R_s2		= R->R_c4 / SQR(R->R_c2);
	}
	wsReal R_ss1	= SQR(R_s1);

	wsReal R_a		= W0;
	if (R_s2 != W0)
		R_a			= R_ss1>R_s2 ? W1/(R_s1-sqrt(R_ss1-R_s2)) :
		W1 / sqrt(R_s2);
	R->R_delta		= R_ss1>R_s2 ? R_s1*CUBE(R_a) - SQR(R_a) : W0;
	R->R_l			= R_a*R_a - (R->R_delta+R->R_delta);
	R->R_muX		= R->R_l + R->R_delta;
	R->R_sigX		= sqrt(W2 * (R->R_l + W2*R->R_delta));
	pverbose("rest %s\n", x.getReadable());

	return R;
}

// 	integrateSKATo<-function(x,qmins,rRho,rhos,muQ,sigQ,sigXi,lius){
void integrateOptimal(double *Ra_X, int N, xIntgInp *Vp_i)
{
	xIntgSKATO*	Xp_s	= (xIntgSKATO *)(Vp_i->Vp_param);
	xLiuCov*	Xp_liu	= Xp_s->Xp_liuFinal;
	double *dx = NULL;
	wsAlloc(dx, double, N);

	for (int i=0 ; i<N ; i++) {
		// 		fff <- function(xxx,qmins,rRho,rhos) min( (qmins - rRho*xxx)/(1-rhos) )
		wsReal R_min = (Xp_s->Ra_qLiu[0] - Xp_s->Ra_rRho[0]*Ra_X[i]) /
			(1.0 - Xp_s->Ra_rho[0]);
		for (wsUint j=1 ; j<Xp_s->N_rho ; j++) {
			wsReal R_res = (Xp_s->Ra_qLiu[j] - Xp_s->Ra_rRho[j]*Ra_X[i])
				/ (1.0 - Xp_s->Ra_rho[j]);
			if (R_res < R_min)
				R_min = R_res;
		}
		// 			dx  <- unlist( lapply(x,fff,qmins=qmins,rRho=rRho,rhos=rhos) )
		dx[i] = R_min;
	}
	wsReal R_mul = sqrt(SQR(Xp_s->R_sigQ) - SQR(Xp_s->R_sigXi)) / Xp_s->R_sigQ;
	for (int i=0 ; i<N ; i++) {
		// 			deltax<-(dx - muQ)*sqrt(sigQ^2 - sigXi^2)/sigQ + muQ
		dx[i] = (dx[i] - Xp_s->R_muQ) * R_mul + Xp_s->R_muQ;

		// pliu(deltax,lius,Lower=TRUE)*dchisq(x,df=1)
		//	tt <- (quan-lius$muQ)/lius$sigQ*lius$sigX + lius$muX
		wsReal R_t = (dx[i]-Xp_liu->R_muQ)/Xp_liu->R_sigQ*Xp_liu->R_sigX
			+ Xp_liu->R_muX;
		//	pchisq(tt,df=lius$l,ncp=lius$delta,lower=Lower)

		// 			pliu(deltax,lius,Lower=TRUE)*dchisq(x,df=1)
		double R_nc = (double)(Xp_liu->R_delta);
		Ra_X[i] = (1.0 - PVchisq(R_t, Xp_liu->R_l, &R_nc, 1))
			* dchisq(Ra_X[i], 1.0, 0);
	}

	DEALLOC(dx);
}
// 	}


// qliu   <- function(pval,lius,Lower=FALSE){
// 	(qchisq(pval,df=lius$l,ncp=lius$delta,lower=Lower)-lius$muX)*lius$sigQ/lius$sigX + lius$muQ
// }
wsReal qLiu(wsReal R_pVal, xLiuCov *Xp_liuCov)
{
	int N_which = 2, N_stat;
	double R_pValInp = 1.0 - (double)R_pVal;
	double R_X, R_bound;
	double R_df = (double)(Xp_liuCov->R_l);
	double R_nc = (double)(Xp_liuCov->R_delta);
	if (R_pValInp == 1.0) /* In case of invalid range */
		R_X = WISARD_INF;
	else
		cdfchn(&N_which, &R_pValInp, NULL, &R_X, &R_df,
			&R_nc, &N_stat, &R_bound);
	return (R_X - Xp_liuCov->R_muX) *
		Xp_liuCov->R_sigQ/Xp_liuCov->R_sigX + Xp_liuCov->R_muQ;
}

int thr_GStest(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	// 	cTestAnalysis
	// 			*Cp_anaGS	= (cGSTestAnalysis *)Vp_shareData;
	// 	wsUint	N_rIdx		= ((xGStestInput *)Vp_data)->N_storeIdx;
	// 	mGeneSet_it
	// 			it			= ((xGStestInput *)Vp_data)->X_it;
	// 	wsUint	N_samp		= ((xGStestInput *)Vp_data)->N_testSamp;
	// 	xGStest	*Xp_res		= (xGStest *)Vp_result + N_rIdx;
	// 	vSnp	Xa_SNP		= ((xGStestInput *)Vp_data)->Cp_IO->getSNPs();
	// 
	// 	wsReal	*Ra_pldCMC	= NULL;
	// 	wsReal	*Ra_pldClp	= NULL;
	// 	sseMalloc(Ra_pldCMC , wsReal, N_samp);
	// 	sseMalloc(Ra_pldClp, wsReal, N_samp);
	// 
	// 	wsUint N_rare = 0, N_common = 0;
	// 	FOREACH (vector<int>::iterator, it->second, iit) {
	// 		wsUint N_idx = *iit;
	// 		if (Xa_SNP[N_idx].maf < OPT_REAL(raremaf))
	// 			N_rare++;
	// 		else
	// 			N_common++;
	// 	}
	// 	wsUint N_col = N_common+(N_rare > 0);
	// 	wsReal	**Ra_PEDCMC = sseMatrix(N_samp, N_col);
	// 
	// 	/* Fill to MISSING_REAL if RV exists */
	// 	if (N_rare) for (wsUint i=0 ; i<N_samp ; i++)
	// 		Ra_PEDCMC[i][N_common] = MISSING_REAL;
	// 
	// 	/* Get pooled data */
	// 	wsReal **Ra_adjGeno = Cp_anaGS->_getGeneWiseStatV2(it,
	// 		Ba_
	// 		Ra_pldCMC, Ra_pldClp, &(Xp_res->R_statSKAT),
	// 		N_common, Ra_PEDCMC);
	// 
	// 	/* Perform SKAT */
	// 	Xp_res->R_pvSKAT = Cp_anaGS->_getSKAT(it, Ra_adjGeno,
	// 		Xp_res->R_statSKAT, &Xp_res->R_pvSKAT);
	// 
	// 	/* Perform CMC */
	// 	Xp_res->R_pvCMC = Cp_anaGS->_getCMC(Ra_pldCMC, &(Xp_res->R_statCMC),
	// 		N_samp);
	// 
	// 	/* Perform Collapsing */
	// 	Xp_res->R_pvClp = Cp_anaGS->_getCollapsing(Ra_pldClp,
	// 		&(Xp_res->R_statClp), N_samp);
	// 
	// 	if (OPT_ENABLED(pedcmc)) {
	// 		Xp_res->R_pvPEDCMC = Cp_anaGS->_getPEDCMC(it,
	// 			Ra_PEDCMC, &(Xp_res->R_statPEDCMC), N_ctrl, N_case,
	// 			N_common, N_rare);
	// 	}
	// 
	// 	sseFree(Ra_pldCMC);
	// 	sseFree(Ra_pldClp);
	// 
	// 	/* For async thread */
	// 	WORKER().setWait(N_idx);
	return 0;
}

int dist_GStest(int N_mode, int N_thread, void *Vp_data,
				void *Vp_shareData, wsUint *Na_waitQueue)
{
	cGSTestAnalysis
		*Cp_anaGS	= (cGSTestAnalysis *)Vp_shareData;
	mGeneDef	&Xa_gs		= Cp_anaGS->getGeneDef();
	xGStestInput
		*Xp_data	= (xGStestInput *)Vp_data;
	static mGeneDef_it
		X_currIdx;
	static wsUint	N_idx = 0;
	int i = 0;

	switch (N_mode) {
	case DISTFUN_INIT:
		X_currIdx = Xa_gs.begin();
		break;
	case DISTFUN_DISTRIBUTE:
		for (i=0 ; i<N_thread && X_currIdx!=Xa_gs.end() ;
			i++,N_idx++,X_currIdx++) {
				Xp_data[i].N_storeIdx	= N_idx;
				Xp_data[i].X_it		= X_currIdx;
				Xp_data[i].Cp_IO		= Cp_anaGS->getIO();
				Xp_data[i].N_testSamp	= Cp_anaGS->getAvailableSampleSize();

				if ((N_idx%10) == 0)
					notice("%d/%d gene-sets processed\r", N_idx, Xa_gs.size());
		}
		break;
	case DISTFUN_UNINIT:
		break;

		/* For async run */
	case DISTFUN_CHECKEND:
		return X_currIdx != Xa_gs.end();
		break;
	case DISTFUN_DISTASYNC:
		for (i=0 ; i<N_thread ; i++) {
			if (Na_waitQueue[i]) {
				if (Xa_gs.end() == X_currIdx)
					Na_waitQueue[i] = 0;
				else {
					Xp_data[i].N_storeIdx = N_idx++;
					Xp_data[i].X_it		= X_currIdx++;
					Xp_data[i].N_testSamp	=
						Cp_anaGS->getAvailableSampleSize();
				}
			}
		}
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt_fmt(WISARD_SYST_INVL_DISTFUNC_CMD, N_mode);
		break;
	}
	return i;
}

int forAllSample_cmc(int N_mode, int N_thread, void *Vp_data,
	void *Vp_shareData)
{
	int		*Np_idx = (int *)Vp_data;
	static wsUint N_idxSample;
	wsUint N_availSamp = ((cGSTestAnalysis *)Vp_shareData)->getAvailableSampleSize();
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		N_idxSample = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++) {
			if ((N_idxSample+i) >= N_availSamp)
				break;
			Np_idx[i] = N_idxSample+i;
			if (((N_idxSample+i)%10) == 0)
				notice("%d/%d samples processed\r", N_idxSample+i, N_availSamp);
		}
		N_idxSample += N_thread;
		break;
	case DISTFUN_UNINIT:
		LOG("%d/%d samples processed\n", N_availSamp, N_availSamp);
		break;
	}
	return i;
}

/* Get residual via linear regression model */
wsReal* getResidual(wsUint N_obs, wsUint N_var, wsVecCst Ra_Y, wsReal **Ra_X)
{
	/* Calculate XtX */
/**/wsSym	Ra_XtX		= sseMtM(Ra_X, N_obs, N_var, Ra_X, N_obs, N_var);
	/* Inv */
/**/wsSym	Ra_XtX_inv	= invSymMat(Ra_XtX, N_var, 1);
	sseUnmat(Ra_XtX, N_var);
	/* Calc X(XtX)^-1 */
/**/wsMat	Ra_XXtXinv	= sseMpMt(Ra_X, N_obs, N_var, Ra_XtX_inv, N_var, N_var);
	sseUnmat(Ra_XtX_inv, N_var);
	/* Calc XtY */
/**/wsVec	Ra_XtY		= sseMtV((wsMatCst)Ra_X, N_obs, N_var, Ra_Y, N_obs);
	/* Get ^y */
/**/wsVec	Ra_yHat		= sseMpV(Ra_XXtXinv, N_obs, N_var, Ra_XtY);
	sseUnmat(Ra_XXtXinv, N_obs);
	sseFree(Ra_XtY);
	/* Get residual */
	wsReal *Ra_res = sseVector(N_obs);
	sseVsV(Ra_Y, Ra_yHat, Ra_res, N_obs);
	sseFree(Ra_yHat);

	return Ra_res;
}

cGSTestAnalysis::cGSTestAnalysis(cIO *Cp_IO, cPPPAnalysisV2 *Cp_inpAnaPPP,
	cAnalysis *Cp_inpAnaPhi,
	cSetManagerAnalysis *Cp_inpGsmAna,
	cEmAiAnalysisV2 *Cp_inpAnaEmai/*=NULL*/) : cAnalysis(Cp_IO)
{
	wsUint	i;
	wsUint	N_gsVrt		= 0;					///< # of variants that included in gene-set
	wsUint	N_sample	= Cp_IO->sizeSample();	///< # of samples in dataset
	wsVecCst	Ra_pheno	= Cp_IO->getPheno();	///< Non-adjusted phenotype
	wsUint	N_pheno		= Cp_IO->sizePheno();

	Mp_phiInv = NULL;

	if (OPT_ENABLED(farvatx)) {
		/* --farvatx && --farvat is m.e.
		 * (well, actually, all of non-X-chromosome methods should not be run w/ this...) */
		if (OPT_ENABLED(farvat))
			halt_fmt(WISARD_CANT_EXCL_OPT, "--farvatx", "--farvat");
		/* Can't do w/ --indep */
		if (OPT_ENABLED(indep))
			halt_fmt(WISARD_CANT_EXCL_OPT, "--farvatx", "--indep");
		/* --farvatx && --autoonly is m.e. */
		if (OPT_ENABLED(autoonly))
			halt_fmt(WISARD_CANT_EXCL_OPT, "--farvatx", "--autoonly");
		/* Check X chromosome is available */
		const char* Ba_chrAllowed = getIO()->getChrAllowance();
		wsUint N_chrXpos = OPT_NUMBER(maxNumAutoChr);
		for (wsStrCst	Sa_nautoChrSeq = OPT_STRING(nonAutoChrSeq) ;
			*Sa_nautoChrSeq != 'X' ; Sa_nautoChrSeq++,N_chrXpos++);
		if (!Ba_chrAllowed[N_chrXpos])
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--farvatx",
				"X chromosome is excluded");

		LOGnote("FARVATX will be performed\n");

		/* w/o --x2, do warning */
		if (!OPT_ENABLED(x2))
			LOGwarn("--farvat is running without --x2, might be incorrect");
		/* Assign --farvatxndiv if not */
		if (!IS_ASSIGNED(farvatxndiv)) {
			OPTION().assign("farvatxndiv", OPTION().getDefVal("farvatxndiv"));
			OPTION().FORCE_OPT_NUMBER(farvatxndiv);
		}
		/* --farvatxd auto on */
		if (!IS_ASSIGNED(farvatxd)) {
			OPTION().assign("farvatxd", OPTION().getDefVal("farvatxd"));
			OPTION().FORCE_OPT_REAL(farvatxd);
		}
	}

	B_isCont	= Cp_IO->isContinuous();
	Vp_yGS		= NULL;
	Mp_Psq		= NULL;
	Vp_yGSunAdj	= NULL;
	Ra_pheVT	= NULL;
	Ba_misPheno	= NULL;
	Ra_phiHat	= NULL;

	/* SKAT and SKAT-o must be not permitted with binary phenotype when
	there is an explicit assignment */
	if (OPT_ENABLED(skat) || OPT_ENABLED(skato)) {
		if (Cp_IO->sizePheno() > 1 || !Cp_IO->isContinuous())
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--skat / --skato",
				"binary phenotype or multiple phenotypes");
	}
	
	if (OPT_ENABLED(farvat) || OPT_ENABLED(farvatx)) {
		/* --skato auto on */
		if (!OPT_ENABLED(skato)) {
			OPTION().assign("skato", "1");
			OPTION().FORCE_OPT_NUMBER(skato);
		}
		B_isCont = 0;
	}

	/* Currently, this analysis only applicable to single phenotype */
	if (Cp_IO->sizePheno() != 1) {
		if (OPT_ENABLED(adjf1) || OPT_ENABLED(adjf2))
			halt("Phenotype adjustment with FQLS only can be applied to single phenotype!");
		if (OPT_ENABLED(logistic))
			halt("Gene-set test with logistic regression fitting cannot applied to multi-phenotype");
		if (IS_ASSIGNED(longitudinal)) {
			if (OPT_ENABLED(mfhom) || OPT_ENABLED(mfhet))
				halt_fmt(WISARD_CANT_EXCL_OPT, "--longitudinal", "--mfhet/--mfhom");
			for (wsUint x=0 ; x<N_pheno ; x++)
				if (!Cp_IO->isContinuous(x))
					halt_fmt(WISARD_CANT_BINPHENO_W_LONGIANA, x+1);
			LOG("Multiple phenotype found, longitudinal analysis will be performed\n");
		} else if (!OPT_ENABLED(mfhom) && !OPT_ENABLED(mfhet))
			halt("Gene-set analysis with multiple phenotype is currently unavailable");

		if (OPT_ENABLED(pedcmc))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--pedcmc", "multiple phenotypes");
		if (OPT_ENABLED(kbac))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--kbac", "multiple phenotypes");
		if (OPT_ENABLED(wsum))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--wsum", "multiple phenotypes");
		if (OPT_ENABLED(ggemma))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--ggemma", "multiple phenotypes");
		if (OPT_ENABLED(imputepheno))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--imputepheno",
				"Multiple phenotype assignment");
		if ((OPT_ENABLED(mfhet) || OPT_ENABLED(mfhom)) && !OPT_ENABLED(skato)) {
			OPTION().assign("skato", "1");
			OPTION().FORCE_OPT_NUMBER(skato);
		}
	} else {
		/* --adjf1 and --adjf2 and --logistic are m.e. */
		if (OPT_ENABLED(adjf1)+OPT_ENABLED(adjf2)+OPT_ENABLED(logistic) > 1)
			halt("--logistic/--adjf1/--adjf2 are mutually exclusive!");
//		if (OPT_ENABLED(mfhet) || OPT_ENABLED(mfhom))
//			halt("--mfhet/--mfhom cannot be used under single phenotype");
		if (IS_ASSIGNED(longitudinal))
			LOGwarn("--longitudinal will be ignored by # of assigned phenotypes[%d]\n", Cp_IO->sizePheno());
		else if (OPT_ENABLED(logistic) || OPT_ENABLED(adjf1) || OPT_ENABLED(adjf2)) {
			if (OPT_ENABLED(adjf1) || OPT_ENABLED(adjf2)) {
				if (!IS_ASSIGNED(heri) || !IS_ASSIGNED(prevalence))
					halt("--adjf1/--adjf2 requires both --heri and --prevalence!");
				if (OPT_ENABLED(adjf2) && !getIO()->isProbandAssigned())
					halt("--adjf2 requires proband information!");
			}
			if (OPT_ENABLED(pedcmc) || OPT_ENABLED(fbskat))
				halt("Cannot perform PEDCMC/FB-SKAT with binary trait via --logistic/--adjf1/--adjf2");
			else if (B_isCont)
				halt("--logistic/--adjf1/--adjf2 option only can be applied to dichotomous phenotype");
			else if ((OPT_ENABLED(logistic) || OPT_ENABLED(adjf1)) && getIO()->isContinuous())
				halt("--adjf1/--adjf2 option only can be applied to dichotomous phenotype");
			else if (OPT_ENABLED(imputepheno))
				halt("--logistic/--adjf1/--adjf2 option cannot used with --imputepheno");
			else {
				LOG("Gene-set test will be performed after logistic regression\n");
				LOG("Note that samples will be considered independent\n");
			}
		} else if (OPT_ENABLED(imputepheno) && B_isCont)
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--imputepheno",
				"Continuous phenotype");
		else if (IS_ASSIGNED(blup)) {
			if (OPT_ENABLED(pedcmc))
				halt_fmt(WISARD_CANT_OPT_W_OTHEROPT, "--pedcmc", "--blup");
		}
	}
	if (!IS_ASSIGNED(genemiss)) {
		OPTION().assign("genemiss", OPTION().getDefVal("genemiss"));
		OPTION().FORCE_OPT_REAL(genemiss);
	}

	Cp_anaGsm	= Cp_inpGsmAna;
	Cp_anaPPP	= Cp_inpAnaPPP;
	Cp_anaCorr	= Cp_inpAnaPhi;
	Cp_anaEmai	= Cp_inpAnaEmai;

	/* Check gene set */
	mGeneDef &Xa_gs = Cp_anaGsm->getGeneDef();
	FOREACH (mGeneDef_it, Xa_gs, it)
		N_gsVrt += (wsUint)(it->second.size());
	LOG("[%d] variants included in [%d] gene-sets\n", N_gsVrt, Xa_gs.size());

	/* Check type of phenotype */
	if (B_isCont) {
		/* If the phenotype is continuous,
		 * pedcmc will not allowed, and skato must be assigned */
		if (OPT_ENABLED(pedcmc))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--pedcmc", "continuous phenotype");
		if (OPT_ENABLED(wsum))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--wsum", "continuous phenotype");
		if (OPT_ENABLED(kbac))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--kbac", "continuous phenotype");
	} else {
		/* Can't do --ggemma with binary phenotype */
		if (OPT_ENABLED(ggemma))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--ggemma", "binary phenotype");
		/* Can't do --famvt with binary phenotype */
		if (OPT_ENABLED(famvt))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--famvt", "binary phenotype");
	}

	/* Check missingness of phenotype */
	if (OPT_ENABLED(mfhom) || OPT_ENABLED(mfhet)) {
		/* Remove NAs and rebuild */
		Ba_misPheno	= Cp_IO->getPheCovMissing(&N_anaSamp);
		N_availPhenoSamp = N_anaSamp;
	} else if (!IS_ASSIGNED(prevalence)) {
		N_anaSamp = 0;
		N_availPhenoSamp = 0;
		wsCalloc(Ba_misPheno, char, N_sample); /* 120917 */
		for (i=0 ; i<N_sample ; i++) {
			if (isMissingReal(Ra_pheno[i])) {
				if (!OPT_ENABLED(imputepheno))
					Ba_misPheno[i] = 1; /* Set to missing only --imputepheno is disabled */
				else
					N_anaSamp++;
			} else {
				N_anaSamp++;
				N_availPhenoSamp++;
			}
		}
	} else {
		N_anaSamp = N_availPhenoSamp = N_sample;
	}
	if (N_anaSamp == 0) {
		LOG("No available sample left, can't do this analysis\n");
		return;
	}

	wsReal *Ra_yGS = _makePheno(B_isCont);		///< Adjusted phenotype
	if (B_isCont)
		LOGnote("Gene-level test with continuous phenotype will be performed\n");
	else
		LOGnote("Gene-level test with dichotomous phenotype will be performed\n");

	/* Make phenotype vector */
	if (IS_ASSIGNED(longitudinal))
		Vp_yGS = new cVector(Ra_yGS, N_anaSamp*N_pheno);
	else
		Vp_yGS = new cVector(Ra_yGS, N_anaSamp);

	/* Create phenotype mask */ {
		wsUint _i=0;
		V_availPhenoMask.init(N_anaSamp, 1);
		wsReal *Ra_mask = V_availPhenoMask.get();
		for (i=0 ; i<N_sample ; i++) {
			if (Ba_misPheno && Ba_misPheno[i]) continue;
			if (isMissingReal(Ra_pheno[i]))
				Ra_mask[_i] = W0;
			_i++;
		}
	}

	/* Build unadjusted phenotype vector */
	if (!B_isCont && !IS_ASSIGNED(blup)) {
		if (Cp_IO->isContinuous() && IS_ASSIGNED(farvat)) {
			Vp_yGSunAdj = NULL;
		} else {
			wsUint _i = 0;
			i = 0;
			vSampPtr Xa_samp = Cp_IO->getSample();

			Vp_yGSunAdj = new cVector(N_anaSamp);
			wsReal *Ra_yUnadjGS = Vp_yGSunAdj->get();
	//		sseMalloc(Ra_yUnadjGS, wsReal, N_sampGS);
			FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
				if (Ba_misPheno && Ba_misPheno[i]) continue;

				if (Ra_pheno[i] == WISARD_UNAFFECTED)
					Ra_yUnadjGS[_i] = W0;
				else if (Ra_pheno[i] == WISARD_AFFECTED)
					Ra_yUnadjGS[_i] = W1;
				else if (isMissingReal(Ra_pheno[i]))
					Ra_yUnadjGS[_i] = WISARD_NA_REAL;
				else
					halt("SYSERR: Should not reach to here");

				_i++;
			}
		}
	} else
		Vp_yGSunAdj = NULL;
//		Ra_yUnadjGS = NULL;

	/* Make VT */
	if (OPT_ENABLED(vt)) {
		if (!IS_ASSIGNED(nperm)) {
			OPTION().assign("nperm", OPTION().getDefVal("nperm"));
			OPTION().FORCE_OPT_NUMBER(nperm);
			LOG("# of permutations in VT set to %d\n", OPT_NUMBER(nperm));
		}

		vSampPtr &Xa_samp	= Cp_IO->getSample();
		/* Build X vector */
		wsMat	Ra_covs		= Cp_IO->getCovariates();
		wsUint	N_origCov	= Cp_IO->sizeCovar();
		wsMat	Ra_regrCovs	= sseMatrix(N_anaSamp, N_origCov+1);

		/* Set intercept */
		for (wsUint j=0,_j=0 ; j<N_sample ; j++) {
			if (Ba_misPheno && Ba_misPheno[j]) continue;
			Ra_regrCovs[_j][0] = W1;
		}
		/* Fill covariates */
		for (i=1 ; i<=N_origCov ; i++) {
			for (wsUint j=0,_j=0 ; j<N_sample ; j++) {
				xSample *Xp_sa = Xa_samp[j];

				if (Ba_misPheno && Ba_misPheno[j]) continue;

				/* Missing check */
				if (isMissingReal(Ra_covs[i][j]))
					halt_fmt(WISARD_NULL_VTCOVAR, i,
						Xp_sa->S_FID.c_str(), Xp_sa->S_IID.c_str());

				Ra_regrCovs[_j][i] = Ra_covs[i][j];
				_j++;
			}
		}

		/* Do regression */
		Ra_pheVT = sseMatrix(N_anaSamp, OPT_NUMBER(nperm));
//		MULTI_MALLOC(Ra_pheVT, wsReal*, N_sampGS);
		wsReal *Ra_oriRes = getResidual(N_anaSamp, N_origCov+1, Ra_yGS, Ra_regrCovs);
		for (wsUint i=0 ; i<N_anaSamp ; i++)
			Ra_pheVT[i][0] = Ra_oriRes[i];
		R_pheVTmean = sseVsum(Ra_oriRes, N_anaSamp) / (wsReal)N_anaSamp;
 		for (int q=1 ; q<OPT_NUMBER(nperm) ; q++) {
			/* Permute Y */
			wsReal *Ra_vals		= NULL;
			sseMalloc(Ra_vals, wsReal, N_anaSamp);
			memcpy(Ra_vals, Ra_oriRes, sizeof(wsReal)*N_anaSamp);

			for (wsUint k=0,j=N_anaSamp ; k<N_anaSamp ; k++,j--) {
				wsUint N_selIdx = rand()%j;

				Ra_pheVT[k][q] = Ra_vals[N_selIdx];
				/* Do erase */
				memmove(Ra_vals+N_selIdx, Ra_vals+N_selIdx+1,
					sizeof(wsReal)*(j-N_selIdx-1));
			}
			sseFree(Ra_vals);
		}
		sseFree(Ra_oriRes);

		/* Cleanup */
		sseUnmat(Ra_regrCovs, N_anaSamp);
	} else {
		Ra_pheVT	= NULL;
		R_pheVTmean	= WISARD_NAN;
	}

#ifdef GSvalidate
	exportMatrix("gs.Y", &Ra_yGS, 1, N_anaSamp*N_pheno);
#endif

	/* Get the phi matrix */
	Ra_phiGS		= NULL;
	Ra_origPhi		= _makeCorr(&N_szPhi);
	if (Ra_origPhi == NULL)
		halt("Failed to get correlation matrix");

	if (OPT_ENABLED(imputepheno)) {
		/* If --imputepheno is on, use all samples even its phenotype is missing */
		N_anaSamp	= N_sample;
		Ra_phiGS	= Ra_origPhi;

		LOG("Sample relatedness matrix will be retained by --imputepheno\n");
	} else if (!IS_ASSIGNED(longitudinal)) {
		/* Otherwise, include samples that are having its phenotype */
		if (N_anaSamp == N_sample)
			Ra_phiGS = Ra_origPhi;
		else {
			LOG("Sample relatedness matrix will reduced from [%d] samples "
				"to [%d] samples from missingness\n", N_sample, N_anaSamp);
		//	Ra_phiGS = sseMatrix(N_anaSamp, N_anaSamp); /* 120917 */

			/* Resizing phi matrix */
			Ra_phiGS = sseMsubsetRect(Ra_origPhi, N_sample, Ba_misPheno,
				0, N_anaSamp);
		}
	} else {
		Ra_phiGS = Ra_origPhi;
	}

	wsMat	Ra_eVec		= NULL;
	wsReal*	Ra_eVal		= NULL;
	Mp_phiInv			= NULL;
	if (IS_ASSIGNED(longitudinal)) {
		Mp_phi = new cIdtMatrix(N_szPhi*N_anaSamp);
		Mp_phiInv = &(Mp_phi->inv());
// 	} else if (OPT_ENABLED(logistic)) {
// 		/* Inverse of phihat(1-phihat) should be inserted */
// 		cVector V_ph(Ra_phiHat, N_anaSamp, 1);
// 		cVector V_ph2 = V_ph.p1p();
// 		cVector V_ph3 = V_ph2.sqrt();
// 		V_ph3.setDontDealloc();
// 		Mp_phi = new cDiagMatrix(N_anaSamp, V_ph3.get());
// 		Mp_phiInv = &(Mp_phi->inv());
	} else if (OPT_ENABLED(logistic) || OPT_ENABLED(indep)) {
		Mp_phi = new cIdtMatrix(N_anaSamp);
		Mp_phiInv = &(Mp_phi->inv());
	} else {
		Mp_phi = new cSymMatrix(Ra_phiGS, N_anaSamp);

		if (OPT_ENABLED(makeev) || IS_ASSIGNED(ev)) {
			wsUint N_r = 0, N_c = 0;
			char S_fn[512];
			wsStrCst	Sp_vv		= IS_ASSIGNED(ev) ? OPT_STRING(ev) : OPT_STRING(out);

#ifdef USE_ED2
			sprintf(S_fn, "%s.eigen.vect", Sp_vv);
#else
			sprintf(S_fn, "%s.eigen.vec", Sp_vv);
#endif
			Ra_eVec = makeMatrixSSE(S_fn, &N_r, &N_c);
#ifdef USE_ED2
			/* Need to be transposed in this part */
			transposeSelf(Ra_eVec, N_anaSamp, N_anaSamp);
#endif
			if (N_r != N_c || N_r != N_anaSamp)
				halt("Sample size[%d * %d] is not match with analysis size[%d]",
				N_r, N_c, N_anaSamp);

			sprintf(S_fn, "%s.eigen.val", Sp_vv);
			wsMat Ra_eval = makeMatrixSSE(S_fn, &N_r, &N_c);
			if (N_r != 1 || N_c != N_anaSamp)
				halt("Sample size[%d * %d] is not match with analysis size[%d]",
				N_r, N_c, N_anaSamp);
			Ra_eVal = Ra_eval[0];
			DEALLOC(Ra_eval);
			LOG("Cached eigendecomposition [%s.eigen.vec] and [%s.eigen.val] loaded\n",
				Sp_vv, Sp_vv);

			Mp_phiInv = &(((cSymMatrix *)Mp_phi)->inv(Ra_eVec, Ra_eVal, N_anaSamp));
		} else
			Mp_phiInv = &(Mp_phi->inv());
	}
#ifdef GSvalidate
	Mp_phi->file("gs.phi");
#endif

	/*** MULTI-FARVAT ***/
	if (OPT_ENABLED(mfhom) || OPT_ENABLED(mfhet)) {
//		wsMat		Ra_yOri	= Cp_IO->getPhenos();
		wsMat		Ra_oCov	= Cp_IO->getCovariates();
		wsUint		N_cv	= Cp_IO->sizeCovar() + 1;
		wsMat		Ra_yRep	= sseMatrix(N_pheno, N_anaSamp);
		cStdMatrix	M_cov(N_cv, N_anaSamp);
		wsMat		Ra_cov	= M_cov.get();
		M_cov.setRow(0, W1);
		for (wsUint i=1 ; i<N_cv ; i++) for (wsUint j=0,J=0 ; j<N_sample ; j++) {
			if (Ba_misPheno && Ba_misPheno[j]) continue;
			Ra_cov[i][J] = Ra_oCov[i-1][j];
			J++;
		}
		for (wsUint j=0,J=0 ; j<N_sample ; j++) {
			if (Ba_misPheno && Ba_misPheno[j]) continue;
			for (wsUint i=0 ; i<N_pheno ; i++)
				Ra_yRep[i][J] = Ra_yGS[j*N_pheno+i];
			J++;
		}
		M_Tt.init(N_pheno, N_anaSamp, Ra_yRep);

		/* Eigendecomposition */
		/* e.ks   <- eigen(akinship)
			* e.P    <- e.ks$vectors
			* e.D    <- e.ks$values */
		cStdMatrix	M_eVec;
		cVector		V_eVal;
		//cSymMatrix	M_phi(Ra_phiGS, N_anaSamp);
		if (Ra_eVal) {
			V_eVal.init(N_anaSamp, Ra_eVal, NULL, 0, 1);
			M_eVec.init(MATOPT_DELNONE, N_anaSamp, N_anaSamp, 0,
				Ra_eVec);
		} else {
			V_eVal = Mp_phi->eigen(M_eVec);
			V_eVal.setDontDealloc();
		}
		/* Z'P
			* Zp     <- t(P)%*%Z */
		cStdMatrix	M_ZtP	= M_cov * M_eVec;
		/* PD
			* e.PD   <- e.P %*% e.D */
		cDiagMatrix	M_D(N_anaSamp, V_eVal.get());
		cStdMatrix	M_PD	= M_eVec * M_D;

		/* Adjust phenotype */
		wsMat		Ra_est	= NULL;
		if (!Cp_anaEmai && B_isCont) {
			cStdMatrix	M_XtY_t		= M_Tt.Mt(M_cov);
			cSymMatrix	M_XtX		= M_cov.Mt();
			cSymMatrix&	M_XtXi		= M_XtX.inv();
			cStdMatrix	M_bHat_t	= M_XtY_t * M_XtXi;
			delete &M_XtXi;
			cVector		V_YtY		= M_Tt.rSS();

			cStdMatrix	M_yHat		= M_bHat_t * M_cov;	/* yHat <- X%*%matrix(solve(XX)%*%Xy,ncol=1) */
			cVector		V_btXtY(N_pheno);
			wsReal*		Ra_btXtY	= V_btXtY.get();
			wsReal		R_nk		= (wsReal)(N_anaSamp-M_cov.row());
			for (wsUint i=0 ; i<N_pheno ; i++) {
				cVector V_curBhat	= M_bHat_t.r2v(i);
				cVector V_curXtY	= M_XtY_t.r2v(i);
				Ra_btXtY[i] = V_curBhat.sum(V_curXtY);
			}
			cVector		V_var		= (V_YtY-V_btXtY) / (R_nk - W1);
			V_var.get();
			Ra_est = sseMatrix(V_var.size(), 2);
			for (wsUint i=0 ; i<V_var.size() ; i++) {
				Ra_est[i][0] = V_var.get()[i];
				Ra_est[i][1] = W0;
			}
		} else if (B_isCont)
			Ra_est = Cp_anaEmai->getEstimates();
		for (wsUint i=0 ; i<N_pheno ; i++) {
			/* Dinv   <- 1 / (sigg2*e.D + sige2) */
			cVector V_y(Ra_yRep[i], N_anaSamp, 1);
			if (B_isCont) {
				/* T is already adjusted 150925 */
#if 0
				wsReal*	Ra_e	= Ra_est[i];

				// Vinv <- solve(sigg2*phi + sige2*diag(totind))
				cMatrix& M_phi = *Mp_phi * Ra_e[1];
				M_phi.addDiag(Ra_e[0]);
				cMatrix& M_Vinv = M_phi.inv();
				// Q <- solve(t(Z)%*%Vinv%*%Z)
				cSymMatrix Qo = M_cov.MMt(M_Vinv);
				cSymMatrix& Q = Qo.inv();
				Qo.rem();
				// P <- Vinv - Vinv%*%Z%*%Q%*%t(Z)%*%Vinv
				cStdMatrix VZt = M_cov * ((cSymMatrix &)M_Vinv);
				cStdMatrix VZ = VZt.transpose();
				cSymMatrix VZQZV = VZ.MMt(Q);
				cSymMatrix P = ((cSymMatrix &)M_Vinv) - VZQZV;
				cStdMatrix QZ = Q * M_cov;
				delete& Q;
				// Py <- P%*%phens
				cVector Py = P * V_y;
				// b <- phig2 %*% Py
				M_phi.addDiag(-Ra_e[0]);
				cVector beta = M_phi * Py;
				delete& M_phi;
				// Vy <- Vinv%*%phens
				cVector ViY = M_Vinv * V_y;
				cVector alpha = QZ * ViY;
				cVector V_offset = alpha * M_cov;
				V_offset += beta;				
#endif	
// /**/			cVector	V_sgD_sI	= V_eVal * Ra_e[1];
// 				V_sgD_sI += Ra_e[0];
// /**/			cDiagMatrix	M_dL(N_anaSamp, V_sgD_sI.get(), 1);
// /**/			cDiagMatrix	M_dLinv = M_dL.inv2();
// 				V_sgD_sI.rem();
// 
// 				/* Check */
// //				cSymMatrix M_RR = M_eVec.MMt(V_eVal.inv());
// //				LOG("Diff %g\n", compSSfrobenius(Mp_phiInv->get(), Mp_phi->row(),
// //					M_RR.get(), M_RR.row()));
// 
// 				/* ZD     <- t(Zp) %*% Dinv */
// /**/			cStdMatrix	M_ZD	= M_ZtP * M_dLinv;
// /**/			MAT_t		Ra_ZD	= M_ZD.get();
// 				/* ZDZi   <- solve(ZD%*%Zp) */
// /**/			SYM_t		Ra_ZDZ	= sym_sseMpMt(Ra_ZD, N_cv, N_anaSamp,
// 					M_ZtP.get(), N_cv, N_anaSamp);
// /**/			cSymMatrix	M_ZDZ(Ra_ZDZ, N_cv);
// /**/			cSymMatrix&	M_ZDZi	= M_ZDZ.inv();
// 				M_ZDZ.rem();
// 
// 				/* QZp */
// 				cStdMatrix	M_QZ	= M_ZDZi * M_ZtP;
// 				cStdMatrix	M_PtZ	= M_ZtP.transpose();
// 				MAT_t		Ra_ZQZ = sym_sseMpM(M_PtZ.get(), M_PtZ.row(), M_PtZ.col(),
// 					M_QZ.get(), M_QZ.row(), M_QZ.col());
// 				/* t(diag(r) - Ki %*% Zp %*% Qi %*% t(Zp))) */
// 				cSymMatrix	M_ZQZ(Ra_ZQZ, N_anaSamp);
// 				cStdMatrix	M_kZQZt	= M_ZQZ * M_dLinv;
// 				cStdMatrix	M_I_ZQZk	= M_kZQZt.transpose();
// 				M_I_ZQZk.subIself();
// 
// 				/* sigg2 * P %*% D */
// 				cStdMatrix&	M_PDx	= M_PD * Ra_e[1];
// 				cStdMatrix	M_v2	= M_PDx.Mt(M_I_ZQZk);
// 				delete &M_PDx;
// 				/* ZQZp */
// 				cStdMatrix	M_v1	= M_cov.tM(M_QZ);
// 				M_v1 += M_v2;
// 
//  				/* Pty
// 				 * Yp     <- t(P) %*% phens */
// 				cVector		V_KYp		= M_eVec.tV(V_y);
// 				V_KYp *= M_dLinv;
// 				cVector		V_offset	= M_v1 * V_KYp;

				/* T <- phens-offset */
				/* T is already adjusted 150925 */
				//V_y -= V_offset;
			}// else
				/* Binary case, offset should be mean */
			//	V_y -= V_y.mean();
		}
		if (!Cp_anaEmai) {
			sseUnmat(Ra_est, N_pheno);
		}
	}
	sseUnmat(Ra_eVec, N_anaSamp);
	sseFree(Ra_eVal);

	/* Get the inverse of phi matrix */

// 	wsReal **Ra_phiInv = NULL; /* 120917 */
// 	if (Cp_IO->getTwinSize() || OPT_ENABLED(ginv)) {
// 		LOG("Use generalized inverse matrix to get phiInv\n");
// 		SVDinverse(Ra_phiGS, N_sampGS, &Ra_phiInv);
// 	} else {
// 		Ra_phiInv = invSymMat(Ra_phiGS, N_sampGS);
// 		if (Ra_phiInv == NULL)
// 			halt("Can't get phiInv matrix");
// 	}
#ifdef GSvalidate
	Mp_phiInv->file("gs.phiInv");
#endif


	// [#samp*#pheno]
	// II         <- matrix(1, nrow=totind, ncol=1)
	// [#samp*#pheno] = [#samp*#samp]%*%[#samp*#pheno]
	// phiInvII   <- akinshipInv%*%II
	// [1]
	// IIphiInvII <- sum(phiInvII)
// 	wsReal *Ra_phiInvII = sseMsumR(Ra_phiInv, N_sampGS, N_sampGS,
// 		&R_IIphiInvII); /* 120917 */
	V_phiInv1 = Mp_phiInv->sumR(&R_1phiInv1);

	// [1]
	// IIphiInvII_Inv <- 1/IIphiInvII
	if (R_1phiInv1 == W0)
		halt("1' %*% phi^-1 %*% 1 = 0, cannot proceed!");
	R_1phiInv1_inv = W1 / R_1phiInv1;
	//LOG("R_IIphiInvII : %g, %g(inv)\n", R_1phiInv1, R_1phiInv1_inv);

	if (OPT_ENABLED(farvatx)) {
		/* IIphiInvII_Inv <- solve(t(II)%*%akinshipInv%*%II) where II = 1 male 2 female */
		cVector		II(W1, N_anaSamp);
		wsVec		Ra_II	= II.get();
		vSampPtr&	Xv_samp	= getIO()->getSample();
		const char*	Ba_fnd	= getIO()->getIsFounder();
		wsUintCst	N_fnd	= getIO()->sizeFounder();
		wsVec		Ra_IIf	= sseVector(N_fnd);
		sseVinit(Ra_IIf, N_fnd, W1);
		wsRealCst	R_v		= OPT_REAL(farvatxd) * W2;

		/* FARVATX */
		wsUint		N_anaFnd = 0;
		wsVec		Ra_phiFnd = sseVector(N_fnd);
		for (wsUint i = 0, I = 0 ; i < N_sample; i++) {
			if (Ba_misPheno && Ba_misPheno[i]) continue;
			/* Check ith sample sex */
			if (Xv_samp[i]->N_sex == 2) {
				if (Ba_fnd[i])
					Ra_IIf[N_anaFnd] = R_v;
				Nv_idxFemale.push_back(I);
				Ra_II[I] = W2;
			}
			if (Ba_fnd[i]) {
				Ra_phiFnd[N_anaFnd] = Ra_origPhi[i][i];
				N_anaFnd++;
				Nv_idxFnder.push_back(I);
			}
			I++;
		}
		cVector		IIf(Ra_IIf, N_anaFnd);
		cDiagMatrix M_piv(N_anaFnd, Ra_phiFnd);
		cDiagMatrix M_pivInv = M_piv.inv();

		wsStrCst S_freq = IS_ASSIGNED(freq) ? OPT_STRING(freq) : OPTION().getDefVal("freq");
		V_phiInv1_Xall = II * *((cSymMatrix *)Mp_phiInv);
		R_1phiInv1_inv_Xall = W1 / II.qf(*((cSymMatrix *)Mp_phiInv));
		if (!stricmp(S_freq, "founder")) {
			V_phiInv1_X = IIf * M_pivInv;
			R_1phiInv1_inv_X = W1 / IIf.qf(M_pivInv);
		} else {
			V_phiInv1_X = V_phiInv1_Xall;
			R_1phiInv1_inv_X = R_1phiInv1_inv_Xall;
		}
		cVector V_11pi1i = II * R_1phiInv1_inv_Xall;
		M_I_pi1pi1i = V_11pi1i.tV(V_phiInv1_Xall);

		/* I - akinshipInv*IIphiInvII_Inv */
		M_I_pi1pi1i.subIself();
	}

	if (OPT_ENABLED(mfhet) || OPT_ENABLED(mfhom)) {
		cVector		V_1T	= M_Tt.sumR();
		cStdMatrix	M_Tsub	= V_1T.tV(V_phiInv1, R_1phiInv1_inv);
		M_TPA		= M_Tt - M_Tsub;

		if (OPT_ENABLED(indep)) {
			cIdtMatrix	M_phi(N_anaSamp);
			cMatrix&	M_PAP	= M_phi - R_1phiInv1_inv;
			M_TPAPT	= M_Tt.MMt(M_PAP);
			delete &M_PAP;
		} else {
			cSymMatrix&	M_PAP	= *(cSymMatrix *)Mp_phi - R_1phiInv1_inv;
			M_TPAPT	= M_Tt.MMt(M_PAP);
			delete &M_PAP;
		}
		R_sumTPAPT	= M_TPAPT.sum();
		V_1TPA		= M_TPA.sumC();
	}

#ifdef GSvalidate
	wsReal *Rp_ipi = &R_1phiInv1;
	exportMatrix("gs.IIphiInvII", &Rp_ipi, 1, 1);
#endif
	/* SKAT-o dedicated part */
//	if (OPT_ENABLED(skato)) {
		/* HH2     <- diag(totind) - matrix(IIphiInvII_Inv, totind,totind)%*%akinshipInv
		 *       [ 1/(1 p^-1 1) ... 1/(1 p^-1 1) ]     [ p^-1_11 ... p^-1_1n ]
		 * I_n - [     ...              ...      ] %*% [   ...         ...   ]
		 *       [ 1/(1 p^-1 1) ... 1/(1 p^-1 1) ]     [ p^-1_n1 ... p^-1_nn ]
		 *
		 *       [ p^-1_+1 / (1 p^-1 1) ... p^-1_+n / (1 p^-1 1) ]
		 * I_n - [ p^-1_+1 / (1 p^-1 1) ... p^-1_+n / (1 p^-1 1) ]
		 *       [ p^-1_+1 / (1 p^-1 1) ... p^-1_+n / (1 p^-1 1) ] */
// 
// 		Ra_HHskato = sseMatrix(N_sampGS, N_sampGS, Ra_phiGS);
// 		for (i=0 ; i<N_sampGS ; i++) {
// 			wsReal R_sumPhiInv = sseVsum(Ra_phiInv[i], N_sampGS) / R_1phiInv1;
// 			for (j=0 ; j<N_sampGS ; j++)
// 				Ra_HHskato[j][i] = -R_sumPhiInv;
// 			Ra_HHskato[i][i] += W1;
// 		}

#ifdef GSvalidate
//		exportMatrix("skato.HH", Ra_HHskato, N_sampGS, N_sampGS);
#endif
//	} else
//		Ra_HHskato = NULL;

	/* --metafarvat */
	if (OPT_ENABLED(makefarvat)) {
		cVector csum = ((cSymMatrix *)Mp_phiInv)->sumC();
		csum *= R_1phiInv1_inv;
		csum.neg();
		wsUint n = Mp_phiInv->row();
		cStdMatrix cc(n, n);
		wsMat ccc = cc.get();
		LOOP (i, n) {
			memcpy(ccc[i], csum.get(), sizeof(wsReal)*n);
			ccc[i][i] += W1;
		}
		V_metaFarvat = *Vp_yGS * cc;
	}

//	sseUnmat(Ra_phiInv, N_sampGS);

	// [#pheno*#samp] = [#pheno*#samp]*1
	// hat    <- t(phiInvII)*IIphiInvII_Inv	# [1,#i]%*%[1]
	V_hat = V_phiInv1 * R_1phiInv1_inv;
//	Ra_hat = Ra_phiInvII;
//	sseVpC(Ra_phiInvII, R_1phiInv1_inv, Ra_hat, N_sampGS);
#ifdef GSvalidate
	V_hat.file("gs.hat");
//	exportMatrix("gs.hat", &Ra_hat, 1, N_sampGS);
#endif

	// [#samp*#samp] = [#samp*#pheno]%*%[#pheno*#samp]
	// hatMat <- II%*%hat
	//  - Actually it copies one 'hat' vector repeats II times
	// [#samp*#samp] = [#samp*#samp]-[#samp*#samp]
	// V      <- diag(totind) - hatMat
// 	wsReal **Ra_V = sseMatrix(N_sampGS, N_sampGS); /* 120917 */
// 	for (i=0 ; i<N_sampGS ; i++) {
// 		for (j=0 ; j<N_sampGS ; j++)
// 			Ra_V[i][j] = -Ra_hat[j];
// 		Ra_V[i][i] += W1;
// 	}
#ifdef GSvalidate
//	exportMatrix("gs.V", Ra_V, N_sampGS, N_sampGS);
#endif

	/* Get the sum of Y*hat */
	//wsReal R_sumYH = (*Vp_yGS * V_hat).sum();
	wsReal R_ySum = Vp_yGS->sum();

	// [#pheno*#samp] = [#pheno*#samp]%*%[#samp*#samp]
	// aV     <- t(YY)%*%V
	/* aV[i] <- YY[i] - sumYH */
	cVector V_tmp = V_hat * R_ySum;
	V_aV = *Vp_yGS - V_tmp;

#ifdef GSvalidate
	V_aV.file("gs.aV");
//	exportMatrix("gs.aV", Ra_aV, 1, N_sampGS);
#endif
//	sseUnmat(Ra_V, N_sampGS);

	/* Build weight */
	wsReal*		Rp_wGS	= getIO()->getWeight(Cp_anaPPP);
	if (Rp_wGS)
		Vp_wGS = new cVector(Rp_wGS, getIO()->sizeVariant(), 1);
	else
		Vp_wGS = new cUnitVector(getIO()->sizeVariant());
//	Ra_wGS = _makeWeight(); /* 120917 */	///< Weight for variants
#ifdef GSvalidate
//	Vp_wGS->file("gs.weight");
//	exportMatrix("gs.weight", &Ra_wGS, 1, Cp_IO->getVariantSize());
#endif

	/* These variables will be determined in run() when --skato */
	Mp_P = NULL;
	Mp_Psq = NULL;
}

cGSTestAnalysis::~cGSTestAnalysis()
{
	if (Mp_phiInv) delete Mp_phiInv;
	if (Vp_yGS) delete Vp_yGS;
//	if (Ra_aV) sseUnmat(Ra_aV, 1);
	/* Deallocate matrix only if it is newly created in this analysis */
//	if (Ra_phiGS != Ra_origPhi)
//		sseUnmat(Ra_phiGS, N_sampGS);

	/* SKAT-o dedicated part */
//	if (Ra_HHskato)	sseUnmat(Ra_HHskato, N_sampGS);
//	if (Ra_P)		sseUnmat(Ra_P, N_sampGS);
//	if (Cp_P)		delete Cp_P;
	if (Mp_Psq)		delete Mp_Psq;

	/* VT */
	if (OPT_ENABLED(vt))
		sseUnmat(Ra_pheVT, N_anaSamp);
	if (Vp_yGSunAdj) delete Vp_yGSunAdj;

//	sseFree(Ra_yGS);
//	sseFree(Ra_yUnadjGS);
//	sseFree(Ra_hat);
	DEALLOC(Ba_misPheno);
// 	if (Ra_logisticV)
// 		sseFree(Ra_logisticV);
	if (Ra_phiHat)
		sseFree(Ra_phiHat);
}

wsRealCst _qlsUniv(cMask &V_mask, cVector& V_aV, cVector &V_pooled, wsReal R_vars,
	wsReal *Rp_stat, wsReal *Rp_aVy=NULL)
{
	wsReal R_aVy = V_mask.vv(V_aV, V_pooled);
	/* Store result */
	if (Rp_aVy) *Rp_aVy = R_aVy;
	// stat <- aVy^2/vars
	*Rp_stat = SQR(R_aVy)/R_vars;
	/* Check and return */
	if (*Rp_stat != *Rp_stat) return WISARD_NA_REAL;
	else
		// pchisq(stat, df=1, lower=F)
		return (wsReal)PVchisq(*Rp_stat, 1.0);
}

// REAL_c cGSTestAnalysis::_qlsUniv(wsReal *Ra_maskTest,
// 	wsReal *Ra_pooled, wsReal R_vars, wsReal *Rp_stat, wsReal *Rp_aVy/*=NULL*/)
// {
// 	// aVy  <- (aV%*%genes)[1,1]
// 	wsReal	R_aVy	= sseVVsel(Ra_maskTest, Ra_aV[0], N_sampGS,
// 		Ra_pooled);
// 	if (Rp_aVy)
// 		*Rp_aVy = R_aVy;
// 
// 	// stat <- aVy^2/vars
// 	*Rp_stat = R_aVy*R_aVy/R_vars;
// 
// 	if (*Rp_stat != *Rp_stat)
// 		return WISARD_NA_REAL;
// 	else
// 		// pchisq(stat, df=1, lower=F)
// 		return (wsReal)chiprobP(*Rp_stat, 1.0);
// }

wsUintCst cGSTestAnalysis::_filterByGenotypingRate(const char *S_gs,
	vInt &X_set, char *Ba_filt, wsReal *Ra_mask)//, wsUint *Np_case)
{
	char	**Na_data	= Cp_IO->getGenotype();
	wsUintCst	N_vrt		= Cp_IO->sizeVariant();
	wsUintCst	N_sample	= Cp_IO->sizeSample();
	wsReal	R_filtThr	= OPT_REAL(genemiss);
	wsUint	N_filter	= 0;

//	*Np_case = 0;
	if (OPT_ENABLED(imputepheno)) {
		for (wsUint i=0 ; i<N_sample ; i++) {
			wsUint N = 0;

			FOREACH (vInt_it, X_set, it) {
				wsUint N_idx = (wsUint)*it;	
				if (N_idx >= N_vrt)
					halt("Invalid gene set input %d in %s", N_idx, S_gs);
				char N_geno = Na_data[i][N_idx];

				if (isAvailable(N_geno))
					N++;
			}
			/* Mark as FILTERED if genotyping rate is NOT ENOUGH */
			wsReal R_currSampGpRate = (wsReal)N / (wsReal)X_set.size();
			if (R_currSampGpRate < R_filtThr) {
				Ba_filt[i] = 1;
				N_filter++;
				Ra_mask[i] = MASK_OFF;
 			}// else {
// 				if (Ra_Y[i] == 1)
// 					(*Np_case)++;
// 			}
		}
	} else {
		for (wsUint i=0,_i=0 ; i<N_sample ; i++) {
			if (Ba_misPheno[i]) continue;
			if (_i >= N_anaSamp)
				halt_fmt(WISARD_SYST_INVL_GS_SZSAMP, N_anaSamp);
			wsUint N = 0;

			FOREACH (vInt_it, X_set, it) {
				wsUint N_idx = (wsUint)*it;	
				if (N_idx >= N_vrt)
					halt("Invalid gene set input %d in %s", N_idx, S_gs);
				char N_geno = Na_data[i][N_idx];

				if (isAvailable(N_geno))
					N++;
			}
			/* Mark as FILTERED if genotyping rate is NOT ENOUGH */
			wsReal R_currSampGpRate = (wsReal)N / (wsReal)X_set.size();
			if (R_currSampGpRate < R_filtThr) {
				Ba_filt[i] = 1;
				N_filter++;
				Ra_mask[_i] = MASK_OFF;
 			}// else {
// 				if (Ra_Y[_i] == 1)
// 					(*Np_case)++;
// 			}
			_i++;
		}
	}

	return N_filter;
}

template <typename T>
cStdMatrix _getGeneWiseStatV2(
	xGeneTestParam& X,
	const char *S_gs,
	T**		Na_data,
	vInt&	Xa_set,
	char *Ba_filt,
	wsUintCst N_anaSamp,
	wsUintCst N_actualSamp,
	cVector &V_pldCMC, cVector &V_pldClp, wsUint *Na_clp,
	wsReal *Rp_statSKAT,
	wsUintCst N_cmIdx, wsReal **Ra_PEDCMCca, wsReal **Ra_PEDCMCct,
	wsUint *Np_Ica, wsUint *Np_Ict, wsUint *Np_imp,
	wsReal *Ra_xx, wsReal *Ra_xxG, wsReal **Ra_Wg_t, cStdMatrix& M_genes,
	wsReal *Ra_wGScurr, wsUint N_case, wsUint N_ctrl,
	wsUint *Np_caMi, wsUint *Np_ctMi,
	wsUint *Np_caMa, wsUint *Np_ctMa, 
	wsUint N_colPEDCMC, wsReal **Ra_VTsum1, wsReal *Ra_VTsum2,
	wsReal **Ra_VTgeno, wsMat Ra_Wsum, char **Na_permPhe)
{
	wsUint		N_sample	= X.Cp_IO->sizeSample();
	wsUintCst	N_gs		= (wsUint)Xa_set.size();
	wsMat		Ra_adjG		= NULL;
	wsVec		Ra_yGS		= X.Vp_yGS->get();
	wsReal		R_allSum	= W0;
	wsUint		N_pheno		= IS_ASSIGNED(longitudinal) ? X.Cp_IO->sizePheno() : 1;
	wsMat		Ra_genes	= M_genes.get();
	wsUint		k			= 0;
	wsUint		i, _i, __i, j;

	if (OPT_ENABLED(skat))
		Ra_adjG = sseMatrix(N_gs, N_anaSamp*N_pheno); /* 120917 */

#ifdef GSvalidate
	char t[256];
	sprintf(t, "gs.X.%s", S_gs);
/**/cExporter* ce = cExporter::summon(t);

	sprintf(t, "gs.weight.%s", S_gs);
/**/cExporter* ce2 = cExporter::summon(t);
#endif

	/* For all variants included in current gene set */
	wsReal *R_hd	= NULL;
	wsReal *Ra_hat	= X.V_hat.get();
	if (OPT_ENABLED(skat)) {
		sseMalloc(R_hd, wsReal, N_actualSamp);
	}

	/* PEDCMC */
	cVector V_pcCase(N_case), V_pcCtrl(N_ctrl);
	V_pcCase.set0();
	V_pcCtrl.set0();
	wsReal *Ra_pcCase = V_pcCase.get();
	wsReal *Ra_pcCtrl = V_pcCtrl.get();

	/* Initialize form */
//	sseInit(Ra_pldClp, N_testSamp, MISSING_REAL);
	V_pldCMC.set0();
	V_pldClp.set0();
//	memset(Ra_pldCMC, 0, sizeof(wsReal)*N_sampGS);
//	memset(Ra_pldClp, 0, sizeof(wsReal)*N_sampGS);

	/* Wsum */
	wsUint	N_permVT	= IS_ASSIGNED(nperm) ? min(MIN_PERM_VT, OPT_NUMBER(nperm)) : MIN_PERM_VT;
	wsUint	N_perm		= IS_ASSIGNED(nperm) ? OPT_NUMBER(nperm) : 1000;
	wsUint	N_ni		= 0;
	wsUint*	Na_miu		= NULL; ///< [#perm]
	wsUint*	Na_niu		= NULL; ///< [#perm]
	if (OPT_ENABLED(wsum)) {
		wsAlloc(Na_miu, wsUint, N_perm+1);
		wsAlloc(Na_niu, wsUint, N_perm+1);

		for (i=0 ; i<N_actualSamp ; i++)
			memset(Na_permPhe[i], 0x00, sizeof(char)*N_perm);

		/* Make permuted phenotypes */
		for (i=0 ; i<N_perm ; i++) {
			for (_i=0 ; _i<N_case ; _i++) {
				wsUint N_idx = 0;
				do {
					N_idx = wsRand()%N_actualSamp;
				} while (Na_permPhe[N_idx][i]);
				Na_permPhe[N_idx][i] = 1;
			}
		}
	}

	wsVecCst	Ra_pheno	= X.Cp_IO->getPheno();
	wsVec		Ra_pldCMC	= V_pldCMC.get();
	wsVec		Ra_pldClp	= V_pldClp.get();
	wsVec		Ra_wGS		= X.Vp_wGS->get();
	xMaf*		Xp_maf		= X.Cp_IO->getMAF();
	//wsReal *Ra_yGSunadj	= Vp_yGSunAdj->get();

	wsUint		N_currColPEDCMC = 0;
 	*Np_imp = 0;
	FOREACH (vInt_it, Xa_set, it) {
		wsUint	N_idx	= (wsUint)*it;
		xMaf&	X_maf	= Xp_maf[N_idx];
		wsReal	R_wgt	= Ra_wGS ? Ra_wGS[N_idx] : W1;

		/* Wsum */
		if (OPT_ENABLED(wsum)) {
			N_ni = 0;
			for (i=0 ; i<=N_perm ; i++) {
				memset(Na_miu, 0x00, sizeof(wsUint)*N_actualSamp);
				memset(Na_niu, 0x00, sizeof(wsUint)*N_actualSamp);
			}
		}

		/* SKAT statistics for current gene set */
		wsReal	R_colSum = W0;

		/* Calculate current N_idx's result */
		if (OPT_ENABLED(skat)) {
			for (i=_i=__i=0 ; i<N_sample ; i++) {
				/* [N_sample]	Ba_filt     
				 * [N_sample]	Ba_misPheno */
				if (!Ba_filt[i]) {
					R_hd[__i] = isAvailable(Na_data[i][N_idx]) ?
						Ra_hat[_i]*(wsReal)Na_data[i][N_idx] : W0;
					__i++;
				}
				if (!X.Ba_misPheno[i])
					_i++;
			}
		}

		wsUint N_idxPEDCMC;
		if (X_maf.R_maf < OPT_REAL(raremaf))
			N_idxPEDCMC = N_cmIdx;
		else
			N_idxPEDCMC = N_currColPEDCMC++;

		/* For PEDCMC */
		wsUint Ica = 0;
		wsUint Ico = 0;

		/* SKAT : Get column-wise sum */
		wsReal R_skatSum = W0;
		if (OPT_ENABLED(skat)) {
			for (j=0 ; j<N_actualSamp ; j++)
				R_skatSum += R_hd[j];
		}

		/* For all sample */
		for (i=_i=__i=0 ; i<N_sample ; i++) {
			/* [N_sample]	Ba_filt     
			 * [N_sample]	Ba_misPheno
			 * [N_actSamp]	R_hd
			 * [N_testSamp]	Ra_pldCMC
			 * [N_testSamp]	Ra_pldClp
			 * [N_actSamp]	Ra_adjG
			 * 
			 *  - skat-o
			 *  
			 * [N_gs]		Ra_xxG
			 * [N_gs*N_actSamp]
			 *				Ra_genes
			 * [N_actSamp*N_gs]
			 *				Ra_Wg_t */
			if (X.Ba_misPheno && X.Ba_misPheno[i]) continue;
			if (!Ba_filt[i]) {
				T		T_geno	= Na_data[i][N_idx];
				bool	B_miss	= isMissing(T_geno);
				wsReal	R_geno;

				/* If missing, impute its genotype into MAF*2 */
				if (B_miss) {
					R_geno = X_maf.R_maf*W2;
					(*Np_imp)++;
				} else {
					R_geno = (wsReal)T_geno;
					Na_clp[k] += (wsUint)(T_geno * N_pheno);
				}

				/* SKAT-o */
				if (OPT_ENABLED(skato)) {
					if (X.B_isCont) {
						wsReal R_val = R_wgt*R_geno;
						/* For continuous trait */
						Ra_Wg_t[__i][k] = R_val;
						if (OPT_ENABLED(mfhom) || OPT_ENABLED(mfhet))
							LOOP (l, N_pheno) Ra_genes[k][__i + l] = R_geno;
					} else {
						/* For dichotomous trait */
						if (Ra_xxG)
							Ra_xxG[k] += (wsReal)N_pheno * Ra_xx[__i] * R_geno;
						LOOP (l, N_pheno) Ra_genes[k][__i + l] = R_geno;
// Ra_xxG	[vrt]
// Ra_genes [snp * obs]
					}
				}
				/* Make geno for other cases */
				if ((OPT_ENABLED(kbac) || OPT_ENABLED(makegeno) || OPT_ENABLED(pedgene)) &&
					(!OPT_ENABLED(skato) || X.B_isCont))
					for (wsUint l=0 ; l<N_pheno ; l++)
						Ra_genes[k][__i + l] = R_geno;

				/* PEDCMC */
				if (OPT_ENABLED(pedcmc)) {
					/* Index overflow check */
					if (N_idxPEDCMC >= N_colPEDCMC)
						halt_fmt(WISARD_SYST_INVL_PEDCMCIDX, N_idxPEDCMC,
							N_colPEDCMC);

					/* Set store index of PEDCMC */
					wsReal	*Rp_PEDCMC = NULL;	///<
					wsReal	*Rp_pc = NULL;
					if (Ra_pheno[i] == WISARD_AFFECTED) {
						if (Ica >= N_case)
							halt("SYSERR: Too large rowIdx on PEDCMC case[%d>%d] or ctrl[%d>%d]", Ica, N_case, Ico, N_ctrl);
						Rp_PEDCMC	= Ra_PEDCMCca[N_idxPEDCMC] + Ica;
						Rp_pc		= Ra_pcCase + (Ica++);
					} else if (Ra_pheno[i] == WISARD_UNAFFECTED) {
						if (Ico >= N_ctrl)
							halt("SYSERR: Too large rowIdx on PEDCMC case[%d>%d] or ctrl[%d>%d]", Ica, N_case, Ico, N_ctrl);
						Rp_PEDCMC	= Ra_PEDCMCct[N_idxPEDCMC] + Ico;
						Rp_pc		= Ra_pcCtrl + (Ico++);
					}

					/* PEDCMC :
					 * If this variant is rare, fill RARE index
					 * otherwise, fill common index */
					if (Rp_PEDCMC) {
						if (X_maf.R_maf < OPT_REAL(raremaf)) {
							if (Na_data[i][N_idx] && isAvailable(Na_data[i][N_idx])) {
								*Rp_PEDCMC = 1;
								*Rp_pc = 1;
							}
						} else {
							*Rp_PEDCMC = R_geno;
							if (R_geno > W0)
								*Rp_pc = 1;
						}
					}
				}

				/* Weighted sum */
				if (OPT_ENABLED(wsum) && !B_miss) {
					/* Set ni */
					N_ni++;

					/* If unaffected */
					if (Ra_pheno[i] == WISARD_UNAFFECTED) {
						Na_niu[0]++;
						Na_miu[0] += (char)T_geno;
					}
					/* For permutation */
					for (wsUint q=1 ; q<=N_perm ; q++) if (Na_permPhe[i][q] == WISARD_UNAFFECTED) {
						Na_niu[q]++;
						Na_miu[q] += (char)T_geno;
					}						
				}

				/* CMC and collapsing */
				if (R_geno) {
                    if (!B_miss) for (wsUint q=0 ; q<N_pheno ; q++)
						Ra_pldCMC[_i+q] = 1;
					for (wsUint q=0 ; q<N_pheno ; q++)
						Ra_pldClp[_i+q] += (wsReal)R_geno * R_wgt;
				}

				/* VT */
				if (OPT_ENABLED(vt)) {
					/* For each variant, sum(geno * (PHENO - mean(PHENO))) */
					/* For each variant, geno^2 */
					LOOP (l, N_permVT)
						Ra_VTsum1[k][l] += R_geno * (X.Ra_pheVT[_i][l] - X.R_pheVTmean);
					Ra_VTsum2[k] += SQR(R_geno);
					//Ra_VTgeno[k][_i] = R_geno;
				} else if (OPT_ENABLED(famvt))
					Ra_VTgeno[k][_i] = R_geno;

				/* SKAT */
				if (OPT_ENABLED(skat)) {
					wsReal R_val = ((wsReal)R_geno-R_skatSum)*R_wgt;
					for (wsUint l=0 ; l<N_pheno ; l++)
						Ra_adjG[k][__i + l]	= R_val;
					R_colSum += N_pheno * R_val*Ra_yGS[_i];
// Ra_adjG [vrt * obs]
				}

#ifdef GSvalidate
				ce->fmt("%g ", R_geno);
#endif
				__i += N_pheno;
			} 
			_i += N_pheno;
		} /* END OF sample loop */

		/* SKAT-o */
		if (Ra_wGScurr)
			Ra_wGScurr[k] = R_wgt;

		/* Wsum = compute weight */
		if (OPT_ENABLED(wsum)) {
			wsReal R_qi = (wsReal)(Na_miu[0]+1) / (wsReal)((Na_niu[0]<<1) + 2);
			Ra_Wsum[k][0] = sqrt((wsReal)N_ni * R_qi * (W1-R_qi));

			for (wsUint q=1 ; q<=N_perm ; q++) {
				wsReal R_qi = (wsReal)(Na_miu[q]+1) / (wsReal)((Na_niu[q]<<1) + 2);
				Ra_Wsum[k][q] = sqrt((wsReal)N_ni * R_qi * (W1-R_qi));
			}
		}

		*Np_Ica = Ica;
		*Np_Ict = Ico;
#ifdef GSvalidate
		ce2->fmt("%g ", R_wgt);
		ce->put("\n");
#endif
		if (OPT_ENABLED(skat))
			R_allSum += R_colSum*R_colSum;
		k++;
	}

	/* PEDCMC */
	if (OPT_ENABLED(pedcmc)) {
		*Np_caMi = *Np_ctMi = 0;
		*Np_caMa = *Np_ctMa = 0;
		for (i=0 ; i<N_case ; i++) {
			if (Ra_pcCase[i] > W0) (*Np_caMi)++;
			else (*Np_caMa)++;
		}
		for (i=0 ; i<N_ctrl ; i++) {
			if (Ra_pcCtrl[i] > W0) (*Np_ctMi)++;
			else (*Np_ctMa)++;
		}
	}

	/* Store SKAT statistics */
	if (OPT_ENABLED(skat)) {
		*Rp_statSKAT = R_allSum;
		sseFree(R_hd);
	}

	DEALLOC(Na_miu);
	DEALLOC(Na_niu);
#ifdef GSvalidate
	delete ce;
	delete ce2;
#endif

	return cStdMatrix(N_gs, N_anaSamp*N_pheno, Ra_adjG);
}

// wsReal* cGSTestAnalysis::_makeWeight()
// {
// 	wsUint	i;
// 	REAL_c	*Ra_ppp_sqrt2pq
// 					= Cp_anaPPP->getPPPsq();
// 	///< Number of given variants
// 	wsUint	N_vrt	= Cp_IO->getVariantSize();
// 	wsReal	*Ra_w	= NULL;
// 
// 	sseMalloc(Ra_w, wsReal, N_vrt);
// 
// 	if (IS_ASSIGNED(weight) || IS_ASSIGNED(betaweight)) {
// 		memcpy(Ra_w, Cp_IO->getVariantWeight(), sizeof(wsReal)*N_vrt);
// 	} else if (OPT_ENABLED(noweight)) {
// #ifdef USE_SSE
// 		wsUint N_med = getMed(N_vrt);
// 		for (i=0 ; i<N_med ; i+=sseJmp)
// 			*((sse_t *)(Ra_w+i)) = sseSet(1.0);
// 		for (i=N_med ; i<N_vrt ; i++)
// 			Ra_w[i] = W1;
// #else
// 		for (i=0 ; i<N_vrt ; i++)
// 			Ra_w[i] = W1/Ra_ppp_sqrt2pq[i];
// #endif
// 	} else {
// #ifdef USE_SSE
// 		sse_t	sse_H	= sseSet(1.0);
// 		wsUint	N_med	= getMed(N_vrt);
// 		for (i=0 ; i<N_med ; i+=sseJmp) {
// 			sse_t*	sse_pW	= (sse_t *)(Ra_w + i);
// 			sse_t*	sse_pP	= (sse_t *)(Ra_ppp_sqrt2pq + i);
// 			*sse_pW = sseDiv(sse_H, *sse_pP);
// 		}
// 		for (i=N_med ; i<N_vrt ; i++)
// 			Ra_w[i] = W1/Ra_ppp_sqrt2pq[i];
// #else
// 		for (i=0 ; i<N_vrt ; i++)
// 			Ra_w[i] = W1/Ra_ppp_sqrt2pq[i];
// #endif
// 	}
// 
// 	return Ra_w;
// }

wsReal*	cGSTestAnalysis::_makePheno(char &B_isContinuous)
{
	wsUint		i, _i;
	wsReal		**Ra_pheno	= Cp_IO->getPhenos();
	wsUint		N_pheno		= Cp_IO->sizePheno();
	vSampPtr	&Xa_samp	= Cp_IO->getSample();

	/* Get the structure of --prevalence */
	wsUint		N_prev		= 0;
	wsReal*		Ra_prev		= loadRealValues(OPT_STRING(prevalence), &N_prev);

	/* With --logistic */
	//wsReal *Ra_phiHat = NULL;
	if (OPT_ENABLED(logistic)) {
		if (Cp_IO->sizePheno() > 1 || B_isCont)
			halt("Cannot apply --logistic since #pheno[%d] > 1, or phenotype is continuous",
				Cp_IO->sizePheno());
		cTimer t;

		t.start();
		wsUint N_cov	= Cp_IO->sizeCovar();
		wsUint N_samp	= Cp_IO->sizeSample();
/**/	wsReal **Ra_cov	= Cp_IO->getCovariates();
/**/	wsReal *Ra_coef	= NULL;
/**/	wsReal *Ra_anaY	= NULL;

		/* Update phenotype missingness to incorportate covariate missing */
		for (wsUint j=0 ; j<N_cov ; j++)
			for (i=0 ; i<N_samp ; i++)
				if (isMissingReal(Ra_cov[j][i]))
					Ba_misPheno[i] = 1;
		/* Update N_sampGS */
		N_anaSamp = 0;
		for (i=0 ; i<N_samp ; i++)
			if (!Ba_misPheno[i]) N_anaSamp++;

		sseMalloc(Ra_anaY, wsReal, N_anaSamp);
/**/	wsMat Ra_anaXt = sseMatrix(N_cov+1, N_anaSamp);
		for (i=0 ; i<N_anaSamp ; i++)
			Ra_anaXt[0][i] = W1;
		for (wsUint j=1 ; j<=N_cov ; j++) {
			for (i=_i=0 ; i<N_samp ; i++) {
				if (Ba_misPheno[i]) continue;
				Ra_anaXt[j][_i] = Ra_cov[j-1][i];
				_i++;
			}
		}
		for (i=_i=0 ; i<N_samp ; i++) {
			if (Ba_misPheno[i]) continue;
			Ra_anaY[_i] = Ra_pheno[0][i];
			_i++;
		}

		/* Initialize */
		sseMalloc(Ra_logisticV, wsReal, N_anaSamp);
		sseCalloc(Ra_coef, wsReal, N_cov+1);
		for (i=0 ; i<=N_cov ; i++)
			Ra_coef[i] = W0;

		bool B_conv = false;
		/*int N_iter = */cRegrAnalysis::fitLogisticByNR(N_anaSamp, N_cov+1, Ra_anaXt,
			Ra_anaY, Ra_coef, Ra_logisticV, 20, B_conv, &Ra_phiHat, &R_logisticS);
		//sseVsqrt(Ra_logisticV, N_anaSamp);
		/* Take sqrt */

		/* Change to continuous phenotype */
		B_isContinuous = 1;
		deallocMatrix(Ra_anaXt, N_cov+1, (void *)1);
		sseFree(Ra_coef);
		sseFree(Ra_anaY);
		LOG("[%s] Logistic regression fitting for gene-level test\n",
			t.getReadable());
	} else if (OPT_ENABLED(adjf1) || OPT_ENABLED(adjf2)) {
		N_anaSamp = Cp_IO->sizeSample();
		Ra_logisticV	= NULL;
		Ra_phiHat		= NULL;
	} else {
		Ra_logisticV	= NULL;
		Ra_phiHat		= NULL;
	}
	/* Allocate buffer for adjusted phenotype */
	wsReal		*Ra_Y		= sseEmptyVec(N_anaSamp*N_pheno);


	/* Phenotype adjustment by --prevalence
	 * 
	 * 1. In case of one prevalence
	 *  CTRL -> -P/(1-P)
	 *  CASE -> 1
	 *  MISS -> 0
	 * 
	 * 2. In case of two prevalence
	 *  CTRL -> -P(SEX)/(1-P(SEX)) if SEX is available
	 *       -> 0                  if SEX is missing
	 *  CASE -> 1
	 *  MISS -> 0
	 */

	i	= 0; /* Sample index in DATA     */
	_i	= 0; /* Sample index w/o MISSING */

	/* CASE 1 : BLUP assigned, DO NOT adjust phenotype */
	if (IS_ASSIGNED(blup) || OPT_ENABLED(nolmm)) {
		lverbose("Gene-test phenotype adjustment: --blup or --nolmm\n");
		FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
			if (Ba_misPheno && Ba_misPheno[i]) continue;	

			for (wsUint j=0 ; j<N_pheno ; j++)
				Ra_Y[_i + j] = Ra_pheno[j][i];
			_i += N_pheno;
		}
	}
	/* CASE 2 : FQLS adjustment */
	else if (OPT_ENABLED(adjf1) || OPT_ENABLED(adjf2)) {
		lverbose("Gene-test phenotype adjustment: --adjf1 or --adjf2\n");

		/* Prevalence */
		wsUint		N_prev;
		wsVec		Ra_prev			= loadRealValues(OPT_STRING(prevalence), &N_prev);
		wsSym		Ra_offsetCor	= getFullCorMat(Cp_anaCorr);
		cVector		V_prev(Ra_prev, N_prev);

		wsVec		Ra_adjFullY_fqls	= genAdjPheno(getIO(), Cp_anaEmai,
			V_prev, NULL);
		memset(Ra_Y, 0x00, sizeof(wsReal)*Xa_samp.size()); /* Default set to 0 */
		if (OPT_ENABLED(adjf1))
			genFqlsPheno(getIO(), Ra_adjFullY_fqls, Ra_offsetCor,
				&Ra_Y, NULL);
		else
			genFqlsPheno(getIO(), Ra_adjFullY_fqls, Ra_offsetCor,
				NULL, &Ra_Y);
		sseFree(Ra_adjFullY_fqls);

		/* Use all samples */
		_i = (wsUint)Xa_samp.size();
	}
	/* CASE 2 : --logistic */
	else if (OPT_ENABLED(logistic)) {
		lverbose("Gene-test phenotype adjustment: --logistic\n");

		/* Adjust by p */
		FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
			if (Ba_misPheno[i]) continue;

			/* Get deviance residual */
// 			for (wsUint j=0 ; j<N_pheno ; j++)
// 				if (Ra_pheno[j][i] == WISARD_AFFECTED)
// 					Ra_Y[_i + j] = sqrt(REAL_CONST(-2.0) * log(Ra_phiHat[_i]));
// 				else
// 					Ra_Y[_i + j] = -sqrt(REAL_CONST(-2.0) * log(W1 - Ra_phiHat[_i]));
			/* 140111 fixed to normal residual (y-p) */
			for (wsUint j=0 ; j<N_pheno ; j++)
				Ra_Y[_i + j] = Ra_pheno[j][i] - Ra_phiHat[_i];
			_i += N_pheno;
		}
	}
	/* CASE 2 : LMM fitted */
	else if (Cp_anaEmai) { // && Cp_anaEmai->getBLUP()) {
		lverbose("Gene-test phenotype adjustment: Linear mixed model\n");

		if (!getIO()->isBLUPapplied()) LOOP (j, N_pheno) {
 			wsRealCst *Ra_blup = Cp_anaEmai->getBLUP(j);
 			wsRealCst *Ra_pred = Cp_anaEmai->getPred(j);

			wsUint I = 0;
			i = 0;
			FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
				if (Ba_misPheno[i]) continue;	

				Ra_Y[j + I*N_pheno] = Ra_pheno[j][i] - Ra_blup[I] - Ra_pred[I];
				I++;
				_i++;
			}
		} else LOOP (j, N_pheno) {
			wsUint I = 0;
			i = 0;
			FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
				if (Ba_misPheno[i]) continue;	

				Ra_Y[j + I*N_pheno] = Ra_pheno[j][i];// - Ra_blup[I] - Ra_pred[I];
				I++;
				_i++;
			}
		}
	}
	/* CASE 3 : Binary phenotype */
	else if (!B_isCont) {
		/* CASE 3-2 : --prevalence */
		if (IS_ASSIGNED(prevalence)) {
			lverbose("Gene-test phenotype adjustment: --prevalence\n");

			/* CASE 3-2-1 : Single phenotype */
			if (Cp_IO->sizePheno() == 1) {
				switch (N_prev) {
				case 0xffffffff: { /* Invalid character? */
					vStr &Xv_grp = Cp_IO->getSampGrpIdx2grp();
					LOGnote("Population-wise prevalence will be applied to FARVAT\n");
					if (Xv_grp.size() == 0)
						halt("--prevalence value [%s] is invalid!",
							OPT_STRING(prevalence));

						/* Try to load prevalence by group */
						mStrReal Xm_grpPrev;
						loadGroupValues(OPT_STRING(prevalence), Xm_grpPrev);

						if (OPT_ENABLED(verbose)) FOREACH (mStrReal_it, Xm_grpPrev, it)
							LOG("Population [%s] = %g\n", it->first.c_str(), it->second);

						/* Assign group-wise prevalence */
						FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
							if (Ba_misPheno && Ba_misPheno[i]) continue;

							wsUint N_grp = (*it)->N_grpFst;
							string &S_grp = Xv_grp[N_grp];
							/* Try to find prevalence of group */
							mStrReal_it X_find = Xm_grpPrev.find(S_grp);
							if (X_find == Xm_grpPrev.end())
								halt("Prevalence of group [%s] is undefined!",
									S_grp.c_str());
							wsReal R_curPrev = X_find->second;

							if (Ra_pheno[0][i] == WISARD_UNAFFECTED)
								Ra_Y[_i] = -R_curPrev;
							else if (Ra_pheno[0][i] == WISARD_AFFECTED)
								Ra_Y[_i] = W1 - R_curPrev;
							else
								Ra_Y[_i] = W0;
							_i++;
						}
					} break;
				case 1: /* Equal to male/female */
					LOGnote("Single prevalence will be applied to FARVAT\n");
					FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
//						if (Ba_misPheno && Ba_misPheno[i]) continue;

						if (Ra_pheno[0][i] == WISARD_UNAFFECTED)
							Ra_Y[_i] = -Ra_prev[0];
						else if (Ra_pheno[0][i] == WISARD_AFFECTED)
							Ra_Y[_i] = W1 - Ra_prev[0];
						else
							Ra_Y[_i] = W0;
						_i++;
					}
					break;
				case 2: /* Non-equal to male/female */
					LOGnote("Sex-wise prevalence will be applied to FARVAT\n");
					FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
//						if (Ba_misPheno[i]) continue;

						if (Ra_pheno[0][i] == WISARD_UNAFFECTED) {
							switch ((*it)->N_sex) {
							case 1:
								Ra_Y[_i] = -Ra_prev[0];
								break;
							case 2:
								Ra_Y[_i] = -Ra_prev[1];
								break;
							default: /* FIXME : How to treat in case of missing sex? */
								Ra_Y[_i] = W0;
							}
						}
						else if (Ra_pheno[0][i] == WISARD_AFFECTED) switch ((*it)->N_sex) {
						case 1:
							Ra_Y[_i] = W1 - Ra_prev[0];
							break;
						case 2:
							Ra_Y[_i] = W1 - Ra_prev[1];
							break;
						default: /* FIXME : How to treat in case of missing sex? */
							halt("ERROR of sex");
						} else
							Ra_Y[_i] = W0;

						_i++;
					}
					break;
				}
			}
			/* CASE 3-2-2 : Multiple phenotype */
//			else if (N_prev != Cp_IO->sizePheno())
//				halt("# of prevalences[%d] does not match with # of phenotypes[%d]",
//					N_prev, Cp_IO->sizePheno());
			else FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
				if (Ba_misPheno && Ba_misPheno[i]) continue;

				if (N_prev == 1) {
					for (wsUint j=0 ; j<N_pheno ; j++)
						if (Ra_pheno[j][i] == WISARD_UNAFFECTED)
							Ra_Y[_i+j] = -Ra_prev[0];
						else if (Ra_pheno[j][i] == WISARD_AFFECTED)
							Ra_Y[_i+j] = W1 - Ra_prev[0];
						else
							Ra_Y[_i+j] = W0;
				}
				else for (wsUint j = 0; j<N_pheno; j++)
					if (Ra_pheno[j][i] == WISARD_UNAFFECTED)
						Ra_Y[_i+j] = -Ra_prev[j];
					else if (Ra_pheno[j][i] == WISARD_AFFECTED)
						Ra_Y[_i+j] = W1 - Ra_prev[j];
					else
						Ra_Y[_i+j] = W0;

				_i += N_pheno;
			}
		} /* END OF CASE 3-2 */
		else if (OPT_ENABLED(farvat)) {
			lverbose("Gene-test phenotype adjustment: --farvat\n");

			/* Use 'as is' */
			FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
				if (Ba_misPheno && Ba_misPheno[i]) continue;	

				for (wsUint j=0 ; j<N_pheno ; j++)
					Ra_Y[_i + j] = Ra_pheno[j][i];
				_i += N_pheno;
			}
		}
		/* CASE 3-3 : ERROR : Should not reach to here */
		else halt("Error : No matched case for determining binary Y, should not reach to here");
	}
	/* CASE 4 : Continuous phenotype, but should not reach to here
	 * FIXME : Anyway use 'as is' */
	else {
		lverbose("Gene-test phenotype adjustment: No adjustment\n");

		//	halt("Error : No matched case for determining continuous Y, should not reach to here");
		FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
			if (Ba_misPheno && Ba_misPheno[i]) continue;	

			for (wsUint j=0 ; j<N_pheno ; j++)
				Ra_Y[_i + j] = Ra_pheno[j][i];
			_i += N_pheno;
		}
	}


	/* Check number */
	if (_i != N_anaSamp*N_pheno)
		halt("SYSERR: Unmatched N_availPhenoSamp %d, expected to be %d * %d",
			_i, N_anaSamp, N_pheno);

	sseFree(Ra_prev);
	return Ra_Y;
}

void getVT(wsUint N_curPerm, xMaf *Xp_maf, vInt& Xv_set, wsVec Ra_mafBins, wsUint N_mafBins, wsMat Ra_VTsum1, wsVec Ra_VTsum2,
	wsReal& R_pvalVT, wsReal& R_thrVT, wsReal& R_statVT)
{
	/* Memory for z(T) */
	wsVec	Ra_numerator	= sseEmptyVec(N_curPerm);
	wsVec	Ra_zMax			= sseVector(N_curPerm);
	wsVec	Ra_mafMax		= sseVector(N_curPerm);
	sseVinit(Ra_zMax, N_curPerm, -numeric_limits<wsReal>::infinity());

	wsReal	R_denominator = W0;

	/* For each MAF bin, get stat */
	LOOP(i, N_mafBins) {
		wsReal R_curMAF = Ra_mafBins[i];

		/* Find variants having its MAF to Ra_mafBins[i] */
		wsUint j=0;
		FOREACHDO(vInt_it, Xv_set, x, j++)
			if (Xp_maf[*x].R_maf == R_curMAF) {
				sseVaV(Ra_numerator, Ra_VTsum1[j], Ra_numerator, N_curPerm);
				// 						for (int r=0 ; r<OPT_NUMBER(nperm) ; r++)
				// 							Ra_numerator[r]	+= Ra_VTsum1[j][r];
				R_denominator	+= Ra_VTsum2[j];
			}
		LOOP(r, N_curPerm) {
			wsReal R_curZ = Ra_numerator[r] / sqrt(R_denominator);
			if (R_curZ > Ra_zMax[r]) {
				Ra_zMax[r] = R_curZ;
				Ra_mafMax[r] = R_curMAF;
			}
		}
	}
	/* z(T)=
	*		��_(i=1)^m ��_i(T) ��_(j=1)^n C_ij (��_j-��bar) /
	*		[��_(i=1)^m ��_i(T) ��_(j=1)^n (C_ij)^2 ]^(1/2)
	*/

	/* Find z(max) */
	wsUint x = 0;
	for (wsUint r=1 ; r<N_curPerm ; r++)
		if (Ra_zMax[0] < Ra_zMax[r]) x++;
	R_pvalVT = (wsReal)(x+1) / (wsReal)(N_curPerm+1);
	R_thrVT = Ra_mafMax[0];
	R_statVT = Ra_zMax[0];

	sseFree(Ra_numerator);
	sseFree(Ra_zMax);
	sseFree(Ra_mafMax);
}

void _getPedgene(xMaf* Xp_maf, wsStrCst S_gs, vInt& Xa_set,
	cVector& V_wgt, cVector& V_resid, cStdMatrix& M_genes, wsUintCst N_curSamp,
	wsRealCst R_cFtr, wsReal& R_Tbdn, wsReal& R_Pbdn, wsReal& R_Tkrn, wsReal& R_Pkrn)
{
	wsUint		N_gs	= (wsUint)Xa_set.size();
	wsSym		Ra_R	= sseMcorRow(M_genes.get(), N_gs, N_curSamp);

//	f <- wt * sqrt(maf*(1-maf))
	wsVec Ra_f = sseVectorP(N_gs, V_wgt.get());
	LOOP (i, N_gs) {
		wsReal R_maf = Xp_maf[Xa_set[i]].R_maf;
		Ra_f[i] *= sqrt(R_maf * (W1 - R_maf));
	}

//	fRmat <- (f %o% f) * r.mat
//	var.z <- fRmat * c.factor
	LOOP (i, N_gs)
		sseVpV(Ra_f, Ra_R[i], Ra_R[i], i+1, Ra_f[i]*R_cFtr);
	wsReal R_stDen = sseSsum(Ra_R, N_gs) / R_cFtr;

//# score males according to male.dose
	cSymMatrix	M_kmat;
	cVector		V_bScSj;
	if (0) {//x.chrom & (male.dose !=1) )  {
// 		geno.score <- geno     
// 		geno.score[sex==1,] <- geno[sex==1,]*male.dose       
// ## kernel stat info
// 		kmat <- geno.score %*% diag(wt^2,nrow=length(wt),ncol=length(wt)) %*% t(geno.score)
// ## Burden stat info
// 		burden.score.subject <- as.vector(geno.score %*% wt)        
	} else {
		cStdMatrix M_gt = M_genes.transpose();
		/* kernel stat info */
		M_kmat	= M_gt.MMt(V_wgt);
		// kmat <- geno %*% diag(wt^2,nrow=length(wt),ncol=length(wt)) %*% t(geno)

		/* Burden stat info */
		V_bScSj	= M_gt * V_wgt;
		// burden.score.subject <- as.vector(geno %*% wt)
	}

//########### Burden stat (2-sided)

//	stat.num <- (resid %*% burden.score.subject)
//	factor.sum <- sum(fRmat)
//	stat.denom <- sqrt(factor.sum * c.factor)    
	R_Tbdn = V_resid.sum(V_bScSj) / R_stDen;
//	burden.stat <- stat.num / stat.denom

	R_Pbdn = PVchisq(SQR(R_Tbdn), W1);
//	burden.pval <- pchisq(burden.stat^2, 1,lower.tail=FALSE)

//########## Quadratic kernel stat
//## if only 1 variant, reduces to burden stat, o/w do kernel test
				
	if (N_gs > 1) {
		R_Tkrn = V_resid.qf(M_kmat);
//		stat.kernel <- as.vector(resid %*% kmat %*% resid)
//		eig <-eigen(var.z, symmetric=T, only.values=T)
//		evals <-eig$values[eig$values>1e-6*eig$values[1]]
		wsVec	Ra_eVals = EIGENDECOMPOSITION(Ra_R, N_gs);
				
		// Try Davies first
		R_Pkrn = davies(R_Tkrn, Ra_eVals, N_gs);
		if (NA(R_Pkrn))
			R_Pkrn = pchisqsum(R_Tkrn, 1.0, Ra_eVals, N_gs, 0, 2);
		sseFree(Ra_eVals);
//		if(method=="davies") {
//			pval.kernel <- davies(stat.kernel , evals, acc=acc.davies)$Qq        
//## Davies' method sometimes instable, returns out-of-range p-value
//## Help file suggests lerge acc.davies of 1e-4 or so to fix
//		}
//		if (method=="kounen") {
//## in main method, require(survery) for kounen
//			pval.kernel <- pchisqsum(stat.kernel, rep(1, length(evals)),
//			evals, lower.tail=FALSE, method="saddle")
//		}

		/* davies and kounen sometimes return out of range p-values.
		 * fix to 1 or 0 on either side */
		if (R_Pkrn < 0 || R_Pkrn > 1)
			R_Pkrn = WISARD_NAN;
//		pval.kernel <- min(pval.kernel, 1)
//		pval.kernel <- max(pval.kernel,0)
	} else {
		R_Tkrn = SQR(R_Tbdn);
		R_Pkrn	= R_Pbdn;
//		stat.kernel <- burden.stat^2  
//		pval.kernel <- burden.pval
	}
	sseUnmat(Ra_R, N_gs);
}

double _getFamVT(xMaf* Xp_maf, wsUintCst N_anaSamp, cMatrix* Mp_P,
	wsMat &Ra_VTgeno, wsReal *Ra_wGScurr,
	vVariant &Xv_vrt, vInt& Xa_set, wsUint N_gs, wsReal *Ra_ySKATo)
{
	/* Reorder VTgeno */
	xRealSort*	Xa_sort = NULL;
	wsAlloc(Xa_sort, xRealSort, N_gs);
	wsUint j = 0;
	FOREACHDO (vInt_it, Xa_set, iit, j++) {
		Xa_sort[j].i = j;
		Xa_sort[j].V = Xp_maf[*iit].R_maf;
	}
	/* Sort set */
	qsort(Xa_sort, N_gs, sizeof(xRealSort), sort_real);

	wsMat Ra_newVTgeno = NULL;
	wsReal *Ra_wGSvt = sseVector(N_gs);
	wsAlloc(Ra_newVTgeno, wsReal*, N_gs);
	for (j=0 ; j<N_gs ; j++) {
		Ra_newVTgeno[j] = Ra_VTgeno[Xa_sort[j].i];
		Ra_wGSvt[j] = Ra_wGScurr[Xa_sort[j].i];
	}
	DEALLOC(Ra_VTgeno);
	Ra_VTgeno = Ra_newVTgeno;

	/*
	for(i in 1:nSNP) for(j in 1:i){
	aa[i,j]<-aa[j,i]<-W[i]*W[j]*integrateSKATo(genes[,i],P,genes[,j])
	}
	*/
	wsSym Ra_aa = sseSymMat(N_gs);
	for (wsUint i=0 ; i<N_gs ; i++) {
		wsReal R_wi = Ra_wGSvt[i];
		cVector V_gi(Ra_VTgeno[i], N_anaSamp, 1);
		cVector V_iP = V_gi * *(cSymMatrix *)Mp_P;
		for (wsUint j=0 ; j<=i ; j++)
			Ra_aa[i][j] = R_wi * Ra_wGSvt[j] *
				sseVV(V_iP.get(), N_anaSamp, Ra_VTgeno[j]);
	}

	/*
	tt <- 0;covs[1,1]<-aa[1,1]
	for(i in 2:nSNP) covs[i:nSNP,1] <- covs[1,i:nSNP] <- covs[1,i-1]+aa[1,i]
	for(i in 2:(nSNP-1) ) {
	covs[i,i] <- covs[i,i-1]*2-covs[i-1,i-1]+aa[i,i]
	for(j in (i+1):nSNP){
	covs[i,j]<-covs[j,i]<-covs[i,j-1]+aa[i,j]
	}
	}
	covs[nSNP,nSNP] <- covs[nSNP,nSNP-1]*2 - covs[nSNP-1,nSNP-1]+aa[nSNP,nSNP]
	*/
	wsSym Ra_cov = sseSymMat(N_gs);
	Ra_cov[0][0] = Ra_aa[0][0];
	for (wsUint i=1 ; i<N_gs ; i++) {
		wsReal R_ssum = W0;
		for (wsUint j=0 ; j<i ; j++) {
			wsReal Q = Ra_aa[i][j];
			Ra_cov[i][j] = Ra_cov[i-1][j] + Q + R_ssum;
			R_ssum += Q;
		}
		Ra_cov[i][i] = Ra_cov[i][i-1] + Ra_aa[i][i] + R_ssum;
		/* Ra_cov[i-1] + aa[i] */
		//sseVaV(Ra_cov[i-1], Ra_aa[i], Ra_cov[i], i);
	}
// 	Ra_cov[N_gs-1][N_gs-1] = -Ra_cov[N_gs-1][N_gs-1] +
// 		W2 * Ra_cov[N_gs][N_gs-1] + Ra_aa[i][i];

	/*
	stats     <- rep(W[1]*sum(genes[,1]*YY),nSNP)
	for(i in 2:nSNP){
	stats[i:nSNP] <- stats[i-1] + W[i]*sum(genes[,i]*YY)
	}
	*/
	wsReal R_init = Ra_wGSvt[0] * sseVV(Ra_VTgeno[0], N_anaSamp, Ra_ySKATo);
	cVector V_stat(R_init, N_gs);
	wsReal*	Ra_stat = V_stat.get();
	for (wsUint i=1 ; i<N_gs ; i++) {
		wsReal R_mul = Ra_wGSvt[i] * sseVV(Ra_VTgeno[i], N_anaSamp, Ra_ySKATo);
		Ra_stat[i] = Ra_stat[i-1] + R_mul;
	}

	/* get sd of varcov matrix */
	double *Ra_invSD = NULL;
	wsAlloc(Ra_invSD, double, N_gs);
	for (wsUint i=0 ; i<N_gs ; i++)
		Ra_invSD[i] = W1 / sqrt(Ra_cov[i][i]);
	/*
	tstat <- abs(stats/sds)
	mstat <- max(tstat)	
	*/
	wsReal R_mstat = Ra_stat[0] * Ra_invSD[0];
	for (wsUint i=1 ; i<N_gs ; i++) { 
		wsReal R_tmp = fabs(Ra_stat[i] * Ra_invSD[i]);
		if (R_tmp > R_mstat)
			R_mstat = R_tmp;
	}
				
	wsSym Ra_cor = sseSymMat(N_gs);
	for (wsUint i=0 ; i<N_gs ; i++)
		for (wsUint j=0 ; j<=i ; j++)
			Ra_cor[i][j] = Ra_cov[i][j] * Ra_invSD[i] * Ra_invSD[j];
	DEALLOC(Ra_invSD);

	sseUnmat(Ra_aa, N_gs);
	sseUnmat(Ra_cov, N_gs);
	/* pvals <- 1-pmnorm(rep(mstat,nSNP),mean(0,nSNP),varcov=cors) */
	wsReal R_pvalVTfam = 1.0 - pmnorm(R_mstat, N_gs, Ra_cor);
	sseUnmat(Ra_cor, N_gs);
	return R_pvalVTfam;
}

double m_checkAdaptivePvalue(wsUint N_permCnt1, wsUint N_permCnt2,
	wsUint N_curIdx, wsUint checkPoint, wsUint alternative, wsReal m_alpha)
{
	if (N_curIdx % checkPoint == 0 && checkPoint > 5) {
		//!- adaptive p-value calculation, at an interval of #checkPoint permutations 
		// apply the "six-sigma" rule

		double R_pValAdp = 1.0;
		if (alternative == 1 || alternative == 0) { 
			R_pValAdp = (1.0 * N_permCnt1 + 1.0) / (1.0 * N_curIdx + 1.0);
		}
		else {
			double permcount = (N_permCnt1 < N_permCnt2) ? N_permCnt1 : N_permCnt2;
			R_pValAdp = (2.0 * permcount + 1.0) / (1.0 * N_curIdx + 1.0);
		}

		double R_sd = sqrt(R_pValAdp * (1.0 - R_pValAdp) / (1.0 * N_curIdx));
		double R_6sig = R_pValAdp - 6.0 * R_sd;

		if (R_6sig > m_alpha) 
			return R_pValAdp;
		else
			return 9.0;
	}
	else 
		return 9.0;
}

void _getGGemma(cStdMatrix &M_Xt, wsUintCst N_anaSamp,
	wsUintCst N_covar,
	wsReal R_l0, wsMat Ra_Pt,
	cVector& V_D, cVector *Vp_yp_t, cVector &V_y2p_t,
	wsReal R_pvals[3], cVector &V_clumped)
{
	cDiagMatrix M_D(V_D);
	M_D.setDontDealloc();

	/* Replace last vector */
	wsUint N_var0	= M_Xt.row()-1;
	wsReal *Rp_tmp = M_Xt.get()[N_var0];
	M_Xt.get()[N_var0] = V_clumped.get();

	/* Do GEMMA */
	xFemmaRes X_r;
	cFemmaAnalysis::_doTest(M_Xt, N_anaSamp, N_covar,
		Ra_Pt, &M_D, R_l0, 1, Vp_yp_t, V_y2p_t, &X_r);
	if (X_r.B_failed)
		R_pvals[0] = R_pvals[1] = R_pvals[2] = WISARD_NAN;
	else {
		R_pvals[0] = X_r.R_Plrt;
		R_pvals[1] = X_r.R_Pscore;
		R_pvals[2] = X_r.R_Pwald;
	}

	/* REcover */
	M_Xt.get()[N_var0] = Rp_tmp;
}

wsReal _getKBAC(cStdMatrix& M_genes, char* Ba_misPheno,
	cVector* Vp_yGSunAdj,
	wsUint N_gs, wsUint N_sample,
	char *Ba_filtAll, wsUint N_case, wsUint N_ctrl, wsUint N_actualSamp,
	wsUint N_perm, wsReal R_alphaKBAC)
{
	wsMat	Ra_genes	= M_genes.get();
	wsReal	R_pvalKBAC	= W2;
	char	B_adaptive	= R_alphaKBAC < W1;

	/* Genotype pattern making */
	vDbl Xv_genoID(N_actualSamp);
	wsUint N_side = OPT_ENABLED(kbac2side) ? 2 : 1;

	//!-Compute unique genotype patterns (string) as ID scores (double) 
	bool hasWt = false;
	for (wsUint i=0 ; i<N_actualSamp ; ++i) {
		double vntIdL = 0.0; 
		double vntIdR = 0.0;
		const double ixiix= pow(9.0, 10.0);
		unsigned int lastCnt = 0;
		unsigned int tmpCnt = 0;

		for (wsUint j=0 ; j<N_gs ; ++j) { 
			if (!isMissingReal(Ra_genes[j][i]) &&
				Ra_genes[j][i] >= W1)
				vntIdR += pow(3.0, 1.0 * (j - lastCnt)) * Ra_genes[j][i];
			else 
				continue;
			if (vntIdR >= ixiix) {
				vntIdL = vntIdL + 1.0;
				vntIdR = vntIdR - ixiix;
				lastCnt = lastCnt + tmpCnt + 1;
				tmpCnt = 0;
				continue;
			}
			else { 
				++tmpCnt; 
				continue; 
			}
		}

		// one-to-one "ID number" for a genotype pattern
		Xv_genoID[i] = vntIdL + vntIdR * 1e-10;
		if (!(Xv_genoID[i] != W0) && !hasWt)
			hasWt = true;
	}
	_rank(Xv_genoID, Xv_genoID, RM_DEFAULT);
	pverbose("Genotype pattern loading ranks: ");
	if (OPT_ENABLED(verbose)) FOREACH (vDbl_it, Xv_genoID, iii) {
		verbosenf("%g ", *iii);
	}
	verbosenf("\n");

	vDbl Xv_uniqPatn = Xv_genoID;
	sort(Xv_uniqPatn.begin(), Xv_uniqPatn.end());
	vDbl_it	it = unique(Xv_uniqPatn.begin(), Xv_uniqPatn.end());
	Xv_uniqPatn.resize(it - Xv_uniqPatn.begin()); 
	if (hasWt) Xv_uniqPatn.erase(Xv_uniqPatn.begin());

	/* non-wildtype genotype data is empty. KBAC has nothing to work on. Return p-value 1.0 */
	if (Xv_uniqPatn.size() == 0)
		//std::cout << "**Warning** non-wildtype genotype data is empty. KBAC has nothing to work on. Return p-value 1.0" << std::endl;
			R_pvalKBAC = W1;
	else {

	}

	// count number of sample individuals for each genotype pattern
	wsUint *Na_uniqPatnCnt = NULL;
	wsAlloc(Na_uniqPatnCnt, wsUint, Xv_uniqPatn.size());
	for (size_t u=0 ; u<Xv_uniqPatn.size() ; ++u) 
		Na_uniqPatnCnt[u] = 0;

	for (wsUint i=0 ; i<N_actualSamp ; ++i) {
		// for each sample, identify/count its genotype pattern

		for (size_t u=0 ; u<Xv_uniqPatn.size() ; ++u) {

			if (Xv_genoID[i] == Xv_uniqPatn[u]) {
				// genotype pattern identified
				++Na_uniqPatnCnt[u];
				// count this genotype pattern
				break;
			}
			else;
			// genotype pattern not found -- move on to next pattern
		}
	}
	if (OPT_ENABLED(verbose)) {
		LOG("Number of each unique individual genotype patterns (totaling %d patterns excluding wildtype): ",
                        Xv_uniqPatn.size());
		LOOP(k, Xv_uniqPatn.size())
			LOGnf("%d, ", Na_uniqPatnCnt[k]);
		LOGnf("\n");
        }

	cVector	V_yGSunadjFin(N_actualSamp);
	wsReal*	Ra_yGSunadjFin	= V_yGSunadjFin.get();
	wsReal*	Ra_yGSunadj	= Vp_yGSunAdj->get();

	/* Get Ra_yGSunadjFin */
	wsUint _I = 0;
	for (wsUint i=0,I=0 ; i<N_sample ; ++i) {
		if (Ba_misPheno[i]) continue;
		_I++;
		if (Ba_filtAll[i]) continue;
		Ra_yGSunadjFin[I++] = Ra_yGSunadj[_I - 1];
	}
	if (_I != N_actualSamp) halt("ERROR on counting");
if (OPT_ENABLED(verbose)) {
                LOG("Phenotype pattern : ");
                LOOP(i, N_actualSamp)
                        LOGnf("%g ", Ra_yGSunadjFin[i]);
                LOGnf("\n");

        }
	unsigned int iPermutation = 0;
	unsigned int permcount1 = 0, permcount2 = 0;
	double observedStatistic = 0.0;
	while (iPermutation <= N_perm) {
		// the KBAC statistic. Will be of length 1 or 2
		vDbl Xv_statKBAC(0);

		// two models
		for (wsUint s=0 ; s<N_side ; ++s) {

			//!- count number of sample cases (for the 1st model, or ctrls for the 2nd model) for each genotype pattern
			wsUint *uniquePatternCountsSub = NULL;
			wsCalloc(uniquePatternCountsSub, wsUint, Xv_uniqPatn.size());
			// genotype pattern counts in cases (for the 1st model, or ctrls for the 2nd model) 

			for (wsUint i=0,I=0 ; i<N_sample ; ++i) {
				if (Ba_filtAll[i]) continue;

				if ((s==0 && Ra_yGSunadjFin[I]==WISARD_AFFECTED) ||
					(s==1 && Ra_yGSunadjFin[I]!=WISARD_AFFECTED)) {
					// for each "case (for the 1st model, or ctrls for 2nd model)", identify/count its genotype pattern
					for (unsigned int u = 0; u != Xv_uniqPatn.size(); ++u) {
						if (Xv_genoID[I] != Xv_uniqPatn[u]) continue;

						// genotype pattern identified in cases (for the 1st model, or ctrls for 2nd model)
						++uniquePatternCountsSub[u];
						// count this genotype pattern
						break;
					}
				}
				I++;
			}

			//!- KBAC weights
			vDbl Ra_wgtUniqPatn(Xv_uniqPatn.size());
			// genotype pattern weights, the hypergeometric distribution cmf
			for (size_t u=0 ; u<Xv_uniqPatn.size() ; ++u)
				Ra_wgtUniqPatn[u] = 0.0;

			for (size_t u=0 ; u<Xv_uniqPatn.size() ; ++u) {
				wsUint Q = s==0 ? N_case : N_ctrl;
				/* Hypergeometric kernel */
				Ra_wgtUniqPatn[u] = phyper(uniquePatternCountsSub[u], Na_uniqPatnCnt[u], N_case+N_ctrl-Na_uniqPatnCnt[u], Q);
//	Q - uniquePatternCountsSub[u],
//	Na_uniqPatnCnt[u] - uniquePatternCountsSub[u],
//	N_case+N_ctrl - Na_uniqPatnCnt[u] - Q + uniquePatternCountsSub[u]);

//fisher(Na_uniqPatnCnt[u] - uniquePatternCountsSub[u],
//					uniquePatternCountsSub[u],
//					N_case+N_ctrl - Na_uniqPatnCnt[u] + uniquePatternCountsSub[u] - Q,
//					Q - uniquePatternCountsSub[u]);
			}
				
			/* Binomial kernel */
			//  					Ra_wgtUniqPatn[u] = pbinom(uniquePatternCountsSub[u],
			// 						Na_uniqPatnCnt[u], N_case/N_actualSamp);
			// 						

			// 					if (m_quiet == false && iPermutation == 0) {
			// 						std::cout << "\nUnique genotype patterns weights (model " << s+1 << "):" << std::endl;
			// 						std::cout << uniquePatternWeights << std::endl;
			// 					}

			//!- KBAC statistic: sum of genotype pattern frequencies differences in cases vs. controls, weighted by the hypergeometric distribution kernel
			double R_subKBAC = 0.0;
			double R_case = (double)N_case;
			double R_ctrl = (double)N_ctrl;
			for (size_t u=0 ; u<Xv_uniqPatn.size() ; ++u) { 
				if (s == 0)
					R_subKBAC += ((wsReal)
					uniquePatternCountsSub[u] / R_case -
					((wsReal)(Na_uniqPatnCnt[u] -
					uniquePatternCountsSub[u])) / R_ctrl ) *
					Ra_wgtUniqPatn[u];
				else
					R_subKBAC += ((1.0 *
					uniquePatternCountsSub[u]) / R_ctrl -
					((wsReal)(Na_uniqPatnCnt[u] -
					uniquePatternCountsSub[u])) / R_case ) *
					Ra_wgtUniqPatn[u];
			}

			//std::cout << kbac << std::endl;

			//FIXME
			//gw_round(kbac, 0.0001);
			Xv_statKBAC.push_back(R_subKBAC);
			if (iPermutation == 0 && OPT_ENABLED(verbose)) {
                               	LOG("Unique genotype patterns weights (model %d) : ", s);
	                        LOOP (jj, Xv_uniqPatn.size())
        	                        LOGnf("%g ", Ra_wgtUniqPatn[jj]);
                	        LOGnf("\n");
				LOG("Unique genotype patterns weights uniq (model %d) : ", s);
                                LOOP (jj, Xv_uniqPatn.size())
                                        LOGnf("%d ", uniquePatternCountsSub[jj]);
                                LOGnf("\n");
	                }
		}

		/* Get current permutation's stat */
		double R_statKBAC = 0.0;
		if (Xv_statKBAC.size() == 1)
			//!- one model statistic
			R_statKBAC = Xv_statKBAC[0];
		else if (Xv_statKBAC.size() == 2)
			//!- two model statistic
			R_statKBAC = max(Xv_statKBAC[0], Xv_statKBAC[1]);
		else
			halt("**Error KBAC statistic (Error code -5)");

		if (iPermutation == 0) 
			observedStatistic = R_statKBAC;
		else {
			if (R_statKBAC >= observedStatistic) 
				++permcount1;
			if (R_statKBAC <= observedStatistic)
				++permcount2;
			if (B_adaptive)
				R_pvalKBAC = m_checkAdaptivePvalue(permcount1,
					permcount2, iPermutation, 5000, 0,
					R_alphaKBAC);
		}
		if (R_pvalKBAC <= 1.0)
			break;

		//!- Permutation
		memset(Ra_yGSunadjFin, 0x00, sizeof(wsReal)*N_actualSamp);
		for (wsUint i=0 ; i<N_case ; i++) {
			wsUint N_idx = 0;
			do {
				N_idx = wsRand() % N_actualSamp;
			} while (Ra_yGSunadjFin[N_idx] != WISARD_UNAFFECTED);
			Ra_yGSunadjFin[N_idx] = WISARD_AFFECTED;
		}
		//random_shuffle(m_ydat.begin(), m_ydat.end());
		++iPermutation;
	}

	if (R_pvalKBAC > W1)
		R_pvalKBAC = ((wsReal)permcount1 + 1.0) / ((wsReal)N_perm + 1.0);

	return R_pvalKBAC;
}

template <typename T>
inline void _setFilter(vInt& Xa_set, vVariant& Xv_vrt, wsUint N_sample,
	vInt& Xv_setOri, char* Ba_misPheno, T** T_data, int& N_start, int& N_end)
{
	wsUint		j;

	/* At first, set every variants to exclude */
	vInt Xa_inc;
	Xa_inc.resize(Xv_setOri.size(), 0);

	/* See the genotype of every samples, and marking it to include
	 * if it have genotype that not 0 */
	LOOP (i, N_sample) {
		if (Ba_misPheno && Ba_misPheno[i]) continue;

		j = 0;
		FOREACHDO (vInt_it, Xv_setOri, iit, j++) {
			if (isAvailable(T_data[i][*iit]) && T_data[i][*iit]>0)
				Xa_inc[j] = 1;
		}
	}
	/*
	 * Chromosome filtering
	 */
	j = 0;
	if (OPT_ENABLED(farvatx)) FOREACHDO (vInt_it, Xv_setOri, iit, j++) {
		/* Include only X-chromosome variants w/ --farvatx */
		if (isXChromosome(Xv_vrt[*iit])) Xa_inc[j] = 1;
	} else  FOREACHDO (vInt_it, Xv_setOri, iit, j++) {
		/* Include only autosome w/o --farvatx */
		if (isAutosome(Xv_vrt[*iit])) Xa_inc[j] = 1;
	}

	/* Make final gene-set with variants that mark to be included */
	j = 0;
	FOREACHDO (vInt_it, Xv_setOri, iit, j++) {
		int pos = (int)Xv_vrt[*iit].pos;
		if (Xa_inc[j]) {
			if (j == 0) {
				N_start = N_end = pos;
			} else {
				if (pos < N_start) N_start = pos;
				if (pos > N_end)   N_end   = pos;
			}
			Xa_set.push_back(*iit);
		}
	}
}

/*
 * FARVATX - robust statistics
 */
inline wsReal getRobust(wsSym Ra_cov, wsRealCst _SS, wsUintCst rL, wsRealCst _ES2)
{
	wsReal Sumcov1	= sseSsum(Ra_cov, rL<<1);
	wsReal varS1	= REAL_CONST(8.0)*rL + W2*Sumcov1;
	wsReal ff1		= W2*SQR(_ES2) / varS1;
	wsReal cc1		= varS1 / (W2*_ES2);

	return PVchisq(_SS/cc1, ff1);
}

void do_geneTest(wsUintCst I_gene, wsUintCst N_gene,
	cIO* Cp_IO, char* Ba_misPheno, wsUintCst N_anaSamp, wsUintCst N_perm,
	wsUintCst N_pheno, wsUintCst N_case, wsUintCst N_ctrl,
	wsStrCst S_gs, vInt& Nv_vrtIdx, xMaf* Xp_maf,
	cVector* Vp_yGSunAdj, wsMat Ra_pheVT, wsReal R_pheVTmean,
	cFemmaAnalysis*	Cp_gemma, cVector& V_xx, xGeneTestParam& X)
{
	wsUint		N_szVec = N_anaSamp * N_pheno;
	cVector		V_pldCMC(N_szVec);
	cVector		V_pldClp(N_szVec);

	int			N_start	= -1;					///< Gene start position (based on the dataset)
	int			N_end	= -1;					///< Gene end position (based on the dataset)
	wsUintCst	N_samp	= Cp_IO->sizeSample();	///< Number of samples in final dataset
	wsFmat		Ra_dsg	= Cp_IO->getDosage();
	char**		Na_geno	= Cp_IO->getGenotype();
	vVariant&	Xv_vrt	= Cp_IO->getVariant();	///< Variant info vector
#ifdef GSvalidate
	char fn[256];
#endif

	/*
	 * 
	 * Preparing header
	 *
	 */
	char	S_bufCMC[128]		= { 0, };
	char	S_bufTime[128]		= { 0, };
	char	S_bufClp[128]		= { 0, };
	char	S_bufVT[128]		= { 0, };
	char	S_bufVTfam[128]		= { 0, };
	char	S_bufGGemma[128]	= { 0, };
	char	S_bufAsum[128]		= { 0, };
	char	S_bufPgene[128]		= { 0, };
	char	S_bufSKAT[128]		= { 0, };
	char	S_bufFARVATx[4096]	= { 0, };
	char	S_bufSKATO[1024]	= { 0, };
	char	S_bufWsum[128]		= { 0, };
	char	S_bufPEDCMC[256]	= { 0, };
	char	S_bufKBAC[128]		= { 0, };

	/*
	 *
	 * Main definitions
	 *
	 */
	wsReal	R_statSKAT	= WISARD_NAN,
		R_statCMC, R_statClp;
	wsReal	R_pvSKAT	= WISARD_NAN,
		R_pvCMC, R_pvClp, R_pvSKATO, R_rhoSKATO;
	wsUint	*Na_clp = NULL;
	cTimer c;

	/* Filtering set */
	vInt Xv_set;
	if (Ra_dsg)
		_setFilter(Xv_set, Xv_vrt, N_samp, Nv_vrtIdx, Ba_misPheno,
			Ra_dsg, N_start, N_end);
	else
		_setFilter(Xv_set, Xv_vrt, N_samp, Nv_vrtIdx, Ba_misPheno,
			Na_geno, N_start, N_end);
	wsUint	N_gs		= (wsUint)Xv_set.size();
	char*	S_bufSKATO_sub2 = new char[20*(N_gs*N_gs)];
	S_bufSKATO_sub2[0] = '\0';

	/* Print out stat */
	notice("%d/%d gene processing...(%d->%d variants)\r", I_gene,
		N_gene, Nv_vrtIdx.size(), Xv_set.size());
	if (Xv_set.size() == 0) return;

	/*
		*
		* Sub-test definitions
		* 
		*/

	/* PEDCMC */
	wsReal	R_statPEDCMC, R_pvPEDCMC;
	wsUint	N_miCase, N_miCtrl;
	wsUint	N_maCase, N_maCtrl;
	wsUint	N_PEDCMCdf;
	wsMat	Ra_PEDCMCca = NULL;	///< Buffer contains case genotypes
	wsMat	Ra_PEDCMCct = NULL;	///< Buffer contains control genotypes

	/* Variable threshold */
	wsVec	Ra_mafBins	= NULL;
	wsUint	N_mafBins	= 0;
	wsMat	Ra_VTsum1	= NULL;
	wsMat	Ra_VTgeno	= NULL;
	wsVec	Ra_VTsum2	= NULL;
	if (OPT_ENABLED(vt) || OPT_ENABLED(famvt)) {
		/* Assign memory space */
		Ra_VTsum1  = sseEmptyMatrix(N_gs, OPT_NUMBER(nperm));
		Ra_VTgeno = sseMatrix(N_gs, N_anaSamp);
		sseCalloc(Ra_VTsum2, wsReal, N_gs);

		/* If VT is on, make a new set that sorted by its MAF */
		wsAlloc(Ra_mafBins, wsReal, Xv_set.size());

		/* Insert sets */
		wsUint j = 0;
 		FOREACHDO (vInt_it, Xv_set, iit, j++) {
			wsUint k = 0;
			for (k=0 ; k<N_mafBins ; k++)
				if (Xp_maf[*iit].R_maf == Ra_mafBins[k]) break;
			if (k >= N_mafBins)
				Ra_mafBins[N_mafBins++] = Xp_maf[*iit].R_maf;
 		}
// 
 		/* Sort set */
 		qsort(Ra_mafBins, N_mafBins, sizeof(wsReal), sort_real_unorder);
	}

	/* Weighted-sum */
	wsMat	Ra_Wsum		= NULL;
	wsReal	R_statWsum	= WISARD_NAN;
	wsReal	R_pvalWsum	= WISARD_NAN;
	if (OPT_ENABLED(wsum))
		Ra_Wsum = sseMatrix(N_gs, N_perm);

		/* Is this sample filtered or not in this gene-set?
		 * This value can vary on geneset because the genotyping rate of
		 * each geneset can be differ along with the sample */
/**/	char	*Ba_filtAll	= NULL; ///< [#samp]	Ba_misPheno + exclusion in current gene-set
//		wsReal	*Ra_maskTest	= NULL;	///< [#tsama]	Exclusion in current gene-set
	cMask	V_mask(N_anaSamp, 1);

	/* Initialize */
	if (Ba_misPheno) {
		sseMalloc(Ba_filtAll, char, N_samp);
		memcpy(Ba_filtAll, Ba_misPheno, N_samp);
	} else {
		sseCalloc(Ba_filtAll, char, N_samp);
	}

	c.start();
	/* Individual-wise filtering by genotyping rate */
	wsUintCst N_filt = 0;//_filterByGenotypingRate(S_gs, Xa_set, Ba_filtAll,
					//	V_mask.get());//, &N_case);
	wsUint N_curSamp = N_anaSamp - N_filt;
	wsUint N_imp = 0;

	/* If --pedcmc is on,
	 * Get the number of common/rare variants */
	wsUint N_rare = 0, N_common = 0, N_col = 0;
	if (OPT_ENABLED(pedcmc)) {
		FOREACH (vInt_it, Xv_set, git)
			Xp_maf[*git].R_maf < OPT_REAL(raremaf) ? N_rare++ : N_common++;

		/* Allocate matrix for counting CMC data */
		N_col = N_common+(N_rare > 0);
		Ra_PEDCMCca = sseMatrix(N_col, N_case);
		Ra_PEDCMCct = sseMatrix(N_col, N_ctrl);

		/* Fill to 0 if RV exists */
		if (N_rare) {
			memset(Ra_PEDCMCca[N_col-1], 0x00, sizeof(wsReal)*N_case);
			memset(Ra_PEDCMCct[N_col-1], 0x00, sizeof(wsReal)*N_ctrl);
		}
	}

#ifdef GSvalidate
 	sprintf(fn, "%s.Yt", S_gs);
	Vp_yGS->file(fn);
// 		exportMatrix(fn, &Ra_yGS, 1, N_sampGS);
#endif
	wsUint	N_curCase=0, N_curCtrl=0;
	char	B_SKAToFailed	= 0;
	/* SKAT-o */
	wsVec	Ra_xxG		= NULL;	///< [N_gs] Dichotomous
	cStdMatrix	M_genes;		///< [] Dichotomous
	wsMat	Ra_Wg_t		= NULL;	///< [] Continuous
	wsVec	Ra_wGScurr	= NULL; ///< [N_gs] Both
	if (OPT_ENABLED(skato) || OPT_ENABLED(kbac) || OPT_ENABLED(famvt)) {
		if (X.B_isCont) {
			wsAlloc(Ra_Wg_t, wsReal*, N_curSamp*N_pheno);
			for (wsUint i=0,k=0 ; i<N_curSamp ; i++) {
				wsReal *Ra_t = sseEmptyVec(N_gs);
				for (wsUint l=0 ; l<N_pheno ; l++,k++)
					Ra_Wg_t[k] = Ra_t;
			}
			//Ra_Wg_t		= sseEmptyMatrix(N_actualSamp*N_pheno, N_gs);
		} else {
			/* For dichotomous phenotype */
			sseCalloc(Ra_xxG, wsReal, N_gs);
			M_genes.init(N_gs, N_curSamp*N_pheno);
		}
		sseMalloc(Ra_wGScurr, wsReal, N_gs);
	}
	if ((OPT_ENABLED(kbac) || OPT_ENABLED(makegeno) || OPT_ENABLED(mfhet)
		|| OPT_ENABLED(mfhom) || OPT_ENABLED(pedgene)) &&
		(!OPT_ENABLED(skato) || X.B_isCont)) {
		M_genes.init(N_gs, N_curSamp*N_pheno);
		if (!Ra_wGScurr) sseMalloc(Ra_wGScurr, wsReal, N_gs);
	}

	/* Get the pooled data */
	sseCalloc(Na_clp, wsUint, Xv_set.size());

	// Ra_adjG [vrt * obs]
	// * Ra_pldCMC [obs]
	// * Ra_pldClp [obs]3
	// * Ra_xxG	[vrt]
	// * Ra_genes [vrt * obs]
	// * Ra_Wg_t  [obs * vrt]
	// * Na_clp [vrt]
		
	cStdMatrix	M_adjGeno;
	if (Ra_dsg)
		M_adjGeno = _getGeneWiseStatV2(X, S_gs,
			Ra_dsg,
			Xv_set, Ba_filtAll,
			N_anaSamp, N_curSamp, V_pldCMC, V_pldClp, Na_clp, &R_statSKAT,
			N_common, Ra_PEDCMCca, Ra_PEDCMCct,
			&N_curCase, &N_curCtrl, &N_imp, V_xx.get(), Ra_xxG, Ra_Wg_t,
			M_genes, Ra_wGScurr, N_case, N_ctrl,
			&N_miCase, &N_miCtrl, &N_maCase, &N_maCtrl,
			N_col, Ra_VTsum1, Ra_VTsum2, Ra_VTgeno, Ra_Wsum, X.Na_permPhe);
	else
		M_adjGeno = _getGeneWiseStatV2(X, S_gs,
			Na_geno,
			Xv_set, Ba_filtAll,
			N_anaSamp, N_curSamp, V_pldCMC, V_pldClp, Na_clp, &R_statSKAT,
			N_common, Ra_PEDCMCca, Ra_PEDCMCct,
			&N_curCase, &N_curCtrl, &N_imp, V_xx.get(), Ra_xxG, Ra_Wg_t,
			M_genes, Ra_wGScurr, N_case, N_ctrl,
			&N_miCase, &N_miCtrl, &N_maCase, &N_maCtrl,
			N_col, Ra_VTsum1, Ra_VTsum2, Ra_VTgeno, Ra_Wsum, X.Na_permPhe);

	/* Export if --makegeno */
	if (OPT_ENABLED(makegeno)) {
		char S_fn[512];
		sprintf(S_fn, "%s.geno", S_gs);
		M_genes.file(S_fn, &(X.Xv_smpbuf));
	}

	/* Check --gsmacthr and # of singleton variants */
	wsUint N_singleton = 0;
	int N_mac = 0;
	for (wsUint i=0 ; i<Xv_set.size() ; i++) {
		if (Na_clp[i] == 1) N_singleton++;
		N_mac += Na_clp[i];
	}
	sseFree(Na_clp);

	if (IS_ASSIGNED(gsmacthr) && N_mac < OPT_NUMBER(gsmacthr)) {
		sseFree(Ba_filtAll);
		//			if (OPT_ENABLED(skat)) sseUnmat(M_adjGeno, N_gs);
		if (OPT_ENABLED(skato)) {
			sseFree(Ra_xxG);	
			sseFree(Ra_wGScurr);
			sseUnmat(Ra_VTgeno, N_gs);

			M_genes.rem();
			if (Ra_Wg_t != NULL) {
				for (wsUint l=0 ; l<N_curSamp ; l++)
					sseFree(Ra_Wg_t[l*N_pheno]);
				DEALLOC(Ra_Wg_t);
			}
		}
		if (OPT_ENABLED(pedcmc)) {
			sseUnmat(Ra_PEDCMCca, N_col);
			sseUnmat(Ra_PEDCMCct, N_col);
		}

		// FIXME: Memory dealloc OK?
		return;
	}

	/*
	 *
	 * Sub-test computations
	 *
	 */

	/* KBAC */
	wsReal R_pvalKBAC = WISARD_NAN;
	if (OPT_ENABLED(kbac)) {
		R_pvalKBAC = _getKBAC(M_genes, Ba_misPheno, Vp_yGSunAdj,
			N_gs, N_samp, Ba_filtAll,
			N_case, N_ctrl, N_curSamp, N_perm, X.R_alphaKB);
		sprintf(S_bufKBAC, "	%g", R_pvalKBAC);
	}

	/* Gene-level GEMMA */
	wsReal R_pvalGGemma[3] = { WISARD_NAN, };
	if (OPT_ENABLED(ggemma)) {
		_getGGemma(X.M_XtGG, N_curSamp, Cp_IO->sizeCovar(),
			X.R_L0GG, Cp_gemma->getPt(),
			Cp_gemma->getD(), Cp_gemma->getYp(), Cp_gemma->getY2p(),
			R_pvalGGemma, V_pldClp);
		char *p = S_bufGGemma;
		if (NA(R_pvalGGemma[0])) p += sprintf(S_bufGGemma, "	NA");
		else p += sprintf(S_bufGGemma, "	%g", R_pvalGGemma[0]);
		if (NA(R_pvalGGemma[1])) p += sprintf(p, "	NA");
		else p += sprintf(p, "	%g", R_pvalGGemma[1]);
		if (NA(R_pvalGGemma[2])) sprintf(p, "	NA");
		else sprintf(p, "	%g", R_pvalGGemma[2]);
	}

	/* FIXME: Adaptive sum test */
// 		if (OPT_ENABLED(asum))
// 			_getAsum();

	/* Variable threshold test */
	if (OPT_ENABLED(vt)) {
		wsUint NP = IS_ASSIGNED(nperm) ? min(MIN_PERM_VT, OPT_NUMBER(nperm)) : MIN_PERM_VT;
		wsReal R_statVT = WISARD_NAN;
		wsReal R_thrVT = WISARD_NAN;
		wsReal R_pvalVT = WISARD_NAN;

		/* Do permutation */
		getVT(NP, Xp_maf, Xv_set, Ra_mafBins, N_mafBins, Ra_VTsum1, Ra_VTsum2, R_pvalVT, R_thrVT, R_statVT);
		if (R_pvalVT < REAL_CONST(.1) && N_perm>MIN_PERM_VT) {
			wsUint k = 0;
			FOREACH (vInt_it, Xv_set, it) {
				wsUint		N_idx	= (wsUint)*it;
				xMaf&		X_maf	= Xp_maf[N_idx];

				/* For all sample */
				wsUint		i, _i;
				for (i=_i=0 ; i<N_samp ; i++) {
					if (Ba_misPheno && Ba_misPheno[i]) continue;
					if (Ba_filtAll[i]) continue;
					wsReal	T_geno	= Ra_dsg ? Ra_dsg[i][N_idx] : (wsReal)Na_geno[i][N_idx];
					bool	B_miss	= isMissing(T_geno);
					/* If missing, impute its genotype into MAF*2 */
					wsReal	R_geno	= B_miss ? X_maf.R_maf*W2 : (wsReal)T_geno;

					/* For each variant, sum(geno * (PHENO - mean(PHENO))) */
					/* For each variant, geno^2 */
					for (wsUint l=NP+1 ; l<N_perm ; l++)
						Ra_VTsum1[k][l] += R_geno * (Ra_pheVT[_i][l] - R_pheVTmean);
					_i += N_pheno;
				} /* END OF sample loop */
				k++;
			}
			getVT(N_perm, Xp_maf, Xv_set, Ra_mafBins, N_mafBins, Ra_VTsum1, Ra_VTsum2, R_pvalVT, R_thrVT, R_statVT);
		}

		sprintf(S_bufVT, "	%g	%g	%g", R_statVT, R_thrVT, R_pvalVT);
	} else
		S_bufVT[0] = '\0';

	/* Family-level variable threshold test */
	if (OPT_ENABLED(famvt)) {
		if (N_gs == 0)
			strcpy(S_bufVTfam, "	NA");
		else {
			wsReal R_pvalVTfam = (wsReal)_getFamVT(Xp_maf, N_anaSamp, X.Mp_P,
				Ra_VTgeno, Ra_wGScurr,
				Xv_vrt, Xv_set, N_gs, X.Ra_ySO);
			sprintf(S_bufVTfam, "	%g", R_pvalVTfam);
		}
	}
	sseUnmat(Ra_VTgeno, N_gs);

	/* Weighted-sum test */
	if (OPT_ENABLED(wsum)) {
		wsMat	Ra_gamma	= sseEmptyMatrix(N_perm+1, N_curSamp);
		//wsReal*	Ra_yGSunadjFin	= sseVector(N_actualSamp);
		wsReal*	Ra_yGSunadj	= Vp_yGSunAdj->get();

		/* Compute gamma */
		LOOP (i, N_gs) {
			wsUint	N_midx	= Xv_set[i];
			wsReal	R_curws	= Ra_Wsum[i][0];

			for (wsUint j=0,J=0 ; j<N_samp ; j++) {
				if (Ba_filtAll[j]) continue;
				/* Ensures genotype in here is NOT NA */
				wsReal R_geno = (wsReal)Na_geno[j][N_midx];
				if (R_geno == W0) continue;

				Ra_gamma[0][J] += R_geno/R_curws;
				J++;
			}

			for (wsUint k=0 ; k<N_samp ; k++) {
				if (Ba_filtAll[k]) continue;
				wsReal R_geno = (wsReal)Na_geno[k][N_midx];

				for (wsUint j=1 ; j<=N_perm ; j++)
					Ra_gamma[j][k] += R_geno/Ra_Wsum[i][j];
			}
		}

		/* Sort gammas */
		xRealSort *Xa_sort = NULL;
		wsAlloc(Xa_sort, xRealSort, N_curSamp);
		for (wsUint i=0,I=0,_I=0 ; i<N_samp ; i++) {
			if (Ba_misPheno && Ba_misPheno[i]) continue;
			_I++;
			if (Ba_filtAll[i]) continue;
			Xa_sort[I].i = _I;
			Xa_sort[I].V = Ra_gamma[0][I];
			I++;
		}
		qsort(Xa_sort, N_curSamp, sizeof(xRealSort), sort_real);

		/* get rank-sum */
		wsReal R_x = W0;
		for (wsUint i=0 ; i<N_curSamp ; i++)
			if (Ra_yGSunadj[Xa_sort[i].i] == WISARD_AFFECTED)
				R_x += i+1; /* 1-base rank */
		DEALLOC(Xa_sort);
		//sseFree(Ra_yGSunadjFin);

		/* Get rank-sum of nulls */
		wsReal *Ra_xs = sseVector(N_perm);
		for (wsUint j=1 ; j<=N_perm ; j++) {
			xRealSort *Xa_sort = buildRealSort(Ra_gamma[j], N_curSamp);
			qsort(Xa_sort, N_curSamp, sizeof(xRealSort), sort_real);

			/* get rank-sum */
			wsReal N_x = W0;
			for (wsUint i=0 ; i<N_curSamp ; i++)
				if (X.Na_permPhe[Xa_sort[i].i][j-1] == WISARD_AFFECTED)
					N_x += i+1; /* 1-base rank */
			DEALLOC(Xa_sort);
			Ra_xs[j-1] = N_x;
		}
		sseUnmat(Ra_gamma, N_perm+1);

		/* Do test */
		wsReal	R_mean	= WISARD_NAN;
		wsReal	R_var	= sseVvar(Ra_xs, N_perm, &R_mean);
		sseFree(Ra_xs);
		R_statWsum		= (R_x - R_mean) / sqrt(R_var);
		R_pvalWsum		= PVchisq(SQR(R_statWsum), 1.0);

		if (NA(R_statWsum) || NA(R_pvalWsum))
			sprintf(S_bufWsum, "	NA	NA");
		else
			sprintf(S_bufWsum, "	%g	%g", R_statWsum, R_pvalWsum);
		sseUnmat(Ra_Wsum, N_gs);
	}

	/* actualSamp-safe : Perform SKAT */
	wsUint B_isLiuApplied = 0;
	if (OPT_ENABLED(skat))
		B_isLiuApplied = _getSKAT(X.R_pppMult,
			S_gs, Xv_set, M_adjGeno, R_statSKAT,
			N_curSamp, &R_pvSKAT);
	//notice("[S]");

	/* pedgene */
	if (OPT_ENABLED(pedgene)) {
		/* Statistics */
		wsReal	R_Tbdn	= WISARD_NAN, R_Pbdn = WISARD_NAN;
		wsReal	R_Tkrn	= WISARD_NAN, R_Pkrn = WISARD_NAN;
		cVector V_wgt(Ra_wGScurr, N_gs, 1);

		_getPedgene(Xp_maf, S_gs, Xv_set, V_wgt, X.V_resPG,
			M_genes, N_curSamp, X.R_cFtrPG,
			R_Tbdn, R_Pbdn, R_Tkrn, R_Pkrn);

		/* Print out */
		char S_tb[64], S_pb[64];
		char S_tk[64], S_pk[64];
		if (NA(R_Tbdn)) strcpy(S_tb, "NA"); else sprintf(S_tb, "%g", R_Tbdn);
		if (NA(R_Pbdn)) strcpy(S_pb, "NA"); else sprintf(S_pb, "%g", R_Pbdn);
		if (NA(R_Tkrn)) strcpy(S_tk, "NA"); else sprintf(S_tk, "%g", R_Tkrn);
		if (NA(R_Pkrn)) strcpy(S_pk, "NA"); else sprintf(S_pk, "%g", R_Pkrn);
		sprintf(S_bufPgene, "	%s	%s	%s	%s", S_tb, S_pb,
			S_tk, S_pk);
	}
		
	/* Perform SKAT-O */
	if (OPT_ENABLED(skato)) {
		if (OPT_ENABLED(farvatx)) {
			vReal	Q;
			cVector	II(W1, N_anaSamp);
			wsVec	Ra_W_X = sseVector(N_gs);
			FOREACH (vInt_it, X.Nv_idxFemale, ii) II.get()[*ii] = W2;

			/* Recalc W */
			cVector		V_maf;
			wsStrCst		S_freq = IS_ASSIGNED(freq) ? OPT_STRING(freq) : OPTION().getDefVal("freq");
			if (!stricmp(S_freq, "founder")) {
				cStdMatrix M_subGenes = M_genes.subsetCol(X.Nv_idxFnder);
				V_maf = M_subGenes * X.V_phiInv1_FX;
			} else
				V_maf = M_genes * X.V_phiInv1_FX;
			V_maf *= X.R_1phiInv1_inv_FX;
			wsVec		Ra_maf = V_maf.get();
			LOOP(i, N_gs) Ra_W_X[i] = dbeta(Ra_maf[i], 1.0, 25.0, 0);
			V_maf.rem();					

			/* Xphi           <- t(XX)%*%akinshipInv */
			cStdMatrix	M_Xphi = M_genes * *((cSymMatrix *)X.Mp_phiInv);
			/* Xphi           <- Xphi %*% II */
			cVector		V_Xphi = M_Xphi * II;

			/*
			 * sig.mat
			 */
			/* XphiX          <- t(XX)%*%akinshipInv%*%XX */
			cSymMatrix	M_sigmat(sym_sseMpMt(M_Xphi.get(), M_Xphi.row(), M_Xphi.col(),
				M_genes.get(), M_genes.row(), M_genes.col()), M_Xphi.row());
			/* subXphiX       <- as.numeric(IIphiInvII_Inv)*Xphi%*%t(Xphi) */
			cSymMatrix	M_sigmat_minus = V_Xphi.tV(X.R_1phiInv1_inv_Xall);
			/* Sig.mat        <- (XphiX - subXphiX)/(nrow(XX)-1) */
			M_sigmat -= M_sigmat_minus;
			M_sigmat /= (wsReal)(N_anaSamp - 1);

			/* WsW            <- W %*% Sig.mat %*% W */
			cDiagMatrix	M_W(N_gs, Ra_W_X, 1);
			cSymMatrix	M_WsW		= M_W.MMt(M_sigmat);
			/* WsW2           <- WsW %*% WsW */
			cSymMatrix	M_WsW2		= M_WsW.Mt();
			/* V_eigC         <- eigen(WsW)$values */
			cVector		V_eigC		= M_WsW.eigen();
			wsRealCst		R_sumWsW	= M_WsW.sum();
			wsRealCst		R_diagWsW2	= M_WsW2.tr();
			wsRealCst		R_sumWsW2	= M_WsW2.sum();


			/* Under --farvatx, should consider dd */
			cVector		V_wGScurr(Ra_W_X, N_gs, 1);
			wsUint		N_dd = OPT_NUMBER(farvatxndiv);
			wsReal		R_div = W1 / (wsReal)(N_dd - 1);
			wsReal		R_dd = W0;
/**/			wsMat		Ra_nom = sseMatrix(N_dd, N_anaSamp);
/**/			wsVec		Ra_denom = sseVector(N_dd);
			char*		Sp = (char *)S_bufFARVATx;
			for (wsUint i=0 ; i<N_dd ; i++,R_dd+=R_div) {
				/* R_dd = current dd */

				/* Generate TT %*% Dd matrix */
				wsVec		Ra_TTDd = sseVectorP(N_anaSamp, X.Vp_yGS->get());
				/* Set male to R_dd and compute t(TT) %*% Dd */
				FOREACH (vInt_it, X.Nv_idxFemale, j) Ra_TTDd[*j] *= R_dd;
				cVector		V_TTDd(Ra_TTDd, N_anaSamp);
				memcpy(Ra_nom[i], Ra_TTDd, sizeof(wsReal)*N_anaSamp);
				wsRealCst		R_denomC = sqrt(W1 / V_TTDd.qf(*(cSymMatrix *)(&X.M_HH_SK)));
				X.R_C = W1 / SQR(R_denomC);

				/* I - akinshipInv*IIphiInvII_Inv */
				cStdMatrix	M_UU	= M_W * M_genes;

				/* Ad          <- W%*%t(XX)%*% M_I_pi1pi1i %*% TD */
				cVector		V_Ad0	= V_TTDd * X.M_I_pi1pi1i;
				cVector		V_Ad	= V_Ad0.Mt(M_UU);
				V_Ad *= R_denomC;
				cVector		V_xxG	= V_Ad0.Mt(M_genes);
				// assign(paste("Q",j,sep=""), t(Ad)%*%Ad)
				wsReal		R_Tcal	= V_Ad.ss();
				Q.push_back(R_Tcal);
				// assign(paste("Q",j,sep=""), t(Ad)%*%IM%*%t(IM)%*%Ad)
				wsReal		R_Tbdn	= V_Ad.sum();
				R_Tbdn = SQR(R_Tbdn);
				Q.push_back(R_Tbdn);

				/* Compute SKATo */
				cStdMatrix	M_U;
				if (Ra_Wg_t) {
					M_U.init(N_curSamp*N_pheno, N_gs, Ra_Wg_t);
					M_U.setClean();
				}
				else
					M_U.init(N_curSamp*N_pheno, N_gs);
/**/			xLiuCov**	Xa_lius	= _getSKATo(N_pheno,
					S_gs, Xv_set, M_U,
					V_wGScurr, M_genes, V_xxG.get(), N_curSamp,
					V_TTDd, &R_rhoSKATO, &R_pvSKATO,
					&M_sigmat, X);
				/* Compute burden */
				// (denom.C)^2*(nom.C%*%(IM%*%t(IM))%*%t(nom.C)
				wsReal		R_pBdnX	= PVchisq(R_Tbdn/R_sumWsW, 1.0);
				/* Compute C-alpha */
				wsReal		R_pCalX	= davies(R_Tcal, V_eigC.get(), N_gs);

				/* Export */
				Sp += sprintf(Sp, "	%g	%g	%g", R_pBdnX, R_pCalX, R_pvSKATO);

				/* Clean */
				LOOP (j, X.N_rho) {
					DEALLOC(Xa_lius[j]);
				}
				DEALLOC(Xa_lius);

				Ra_denom[i] = R_denomC;
			} /* END of N_dd */

			wsUint	N_2dd		= N_dd * 2;
/**/			wsSym	Ra_sigQQ	= sseSymMat(N_2dd);
/**/			wsSym	Ra_rhoQQ	= sseSymMat(N_2dd);
			LOOP(i, N_2dd) {
				wsUint N = i >> 1;
				LOOP(j, i) {
					wsUint NN = j >> 1;
					wsReal R_nom = sseVpSpV(Ra_nom[N], X.M_HH_SK.get(),
						N_anaSamp, Ra_nom[NN]);
					wsReal R_denom = Ra_denom[N] * Ra_denom[NN];
					if ((i % 2 == 0) && (j % 2 == 0))
						Ra_sigQQ[i][j] = R_nom*R_denom*R_diagWsW2;
					if ((i % 2 != 0) && (j % 2 != 0)) {
						Ra_sigQQ[i][j] = R_nom*R_denom*SQR(R_sumWsW);
					}
					else {
						Ra_sigQQ[i][j] = R_nom*R_denom*R_sumWsW2;
					}
				}

				/* Diagonal */
				Ra_sigQQ[i][i] = i % 2 ? SQR(R_sumWsW) : R_diagWsW2;
			}
			sseUnmat(Ra_nom, N_dd);
			sseFree(Ra_denom);
			// Calculating rhoQQ
			LOOP(i, N_2dd)
				for (wsUint j = 0; j <= i; j++)
					Ra_rhoQQ[i][j] = Ra_sigQQ[i][j] / sqrt(Ra_sigQQ[i][i] * Ra_sigQQ[j][j]);
			sseUnmat(Ra_sigQQ, N_2dd);
			// Calculating p-values
			/**/			wsVec Ra_pval = sseVector(N_2dd);
			LOOP(i, N_2dd) {
				if (i % 2 == 0)
					Ra_pval[i] = davies(Q[i], V_eigC.get(), N_gs);
				else
					Ra_pval[i] = PVchisq(Q[i] / R_sumWsW, 1.0);
			}
			wsRealCst	_SS = -2 * sseVlogsum(Ra_pval, N_2dd);
			sseFree(Ra_pval);
			wsRealCst	ES2 = 4 * N_dd;

/**/			wsSym	Ra_cov1 = sseSymMat(N_2dd);
/**/			wsSym	Ra_cov2 = sseSymMat(N_2dd);
			LOOP(i, N_2dd) for (wsUint j = 0; j <= i; j++) {
				wsReal R_cc = Ra_rhoQQ[i][j];
				Ra_cov1[i][j] = 3.279*R_cc + 0.711*SQR(R_cc);
				Ra_cov2[i][j] = 3.263*R_cc + 0.710*SQR(R_cc) + 0.027*CUBE(R_cc);
			}
			sseUnmat(Ra_rhoQQ, N_2dd);
			wsReal R_pRob1 = getRobust(Ra_cov1, _SS, N_dd, ES2);
			wsReal R_pRob2 = getRobust(Ra_cov2, _SS, N_dd, ES2);
			sseUnmat(Ra_cov1, N_2dd);
			sseUnmat(Ra_cov2, N_2dd);
			Sp += sprintf(Sp, "	%g	%g", R_pRob1, R_pRob2);

			/* END OF FARVATX */
			sseFree(Ra_W_X);
		} else {
			char	S_bufSKATO_sub[1024];

			if (N_gs > 1) {
				cStdMatrix	M_U;
				if (Ra_Wg_t) {
					M_U.init(N_curSamp*N_pheno, N_gs, Ra_Wg_t);
					M_U.setClean();
				}
				else
					M_U.init(N_curSamp*N_pheno, N_gs);
				cVector		V_ySKATo(X.Ra_ySO, N_szVec, 1);
				cVector		V_wGScurr(Ra_wGScurr, N_gs, 1);
				/* actualSamp-safe : Perform SKAT-o */
/**/			xLiuCov **Xa_lius = _getSKATo(N_pheno,
					S_gs, Xv_set, M_U,
					V_wGScurr, M_genes, Ra_xxG, N_curSamp,
					V_ySKATo, &R_rhoSKATO, &R_pvSKATO, NULL, X);

				sprintf(S_bufSKATO, "	%g	%g	%g", Xa_lius[0]->R_pVal,
					Xa_lius[0]->R_stat, Xa_lius[0]->R_l);
// 					for (i=1 ; i<(N_rho-1) ; i++) {
// 						sprintf(S_bufSKATO_sub, "	%g	%g	%g	%g", Xa_lius[i]->R_pVal,
// 							Xa_lius[i]->R_stat, Xa_lius[i]->R_l, Xa_lius[i]->R_delta);
// 						strcat(S_bufSKATO, S_bufSKATO_sub);
// 					}
				if (isMissingReal(R_pvSKATO) || NA(R_pvSKATO))
					sprintf(S_bufSKATO_sub, "	%g	%g	%g	%g	NA",
						Xa_lius[X.N_rho - 1]->R_pVal, Xa_lius[X.N_rho - 1]->R_stat,
						Xa_lius[X.N_rho - 1]->R_l, R_rhoSKATO);
				else
					sprintf(S_bufSKATO_sub, "	%g	%g	%g	%g	%g",
						Xa_lius[X.N_rho - 1]->R_pVal, Xa_lius[X.N_rho - 1]->R_stat,
						Xa_lius[X.N_rho - 1]->R_l, R_rhoSKATO, R_pvSKATO);
				strcat(S_bufSKATO, S_bufSKATO_sub);

				/* --makefarvat */
				if (OPT_ENABLED(makefarvat)) {
					cVector	V_stats	= M_genes * X.V_metaFarvat;
					wsVec	Ra_stts	= V_stats.get();
					vInt_it	X		= Xv_set.begin();
					wsUint	pos = 0;

					/* Print variant names */
					pos = sprintf(S_bufSKATO_sub2, "	%s", Xv_vrt[*X].name);
					for (X++; X != Xv_set.end(); X++)
						pos += sprintf(S_bufSKATO_sub2 + pos, ",%s", Xv_vrt[*X].name);

					/* Print out MAC and nind */
					X = Xv_set.begin();
					pos += sprintf(S_bufSKATO_sub2 + pos, "	%d",
						Xp_maf[*X].N_allMac);
					for (X++; X != Xv_set.end(); X++)
						pos += sprintf(S_bufSKATO_sub2 + pos, ",%d",
							Xp_maf[*X].N_allMac);

					/* Print out MAC and nind */
/* COMMENTED 160429
					X = Xv_set.begin();
					pos += sprintf(S_bufSKATO_sub2 + pos, "	%d",
						Xa_maf[*X].N_allele);
					for (X++; X != Xv_set.end(); X++)
						pos += sprintf(S_bufSKATO_sub2 + pos, ",%d",
						Xa_maf[*X].N_allele);
					S_bufSKATO_sub2[pos++] = '\t';
*/

					/* Print out stats */
					pos += sprintf(S_bufSKATO_sub2 + pos, "	%g", Ra_stts[0]);
					for (wsUint i = 1; i < N_gs; i++)
						pos += sprintf(S_bufSKATO_sub2 + pos, ",%g", Ra_stts[i]);

					/* Print out matrix */
					cSymMatrix	M_cov = M_genes.covR();
					wsSym		Ra_cov = M_cov.get();
					pos += sprintf(S_bufSKATO_sub2 + pos, "	%g", Ra_cov[0][0]);
					for (wsUint i = 1; i < N_gs; i++) for (wsUint j = 0; j <= i; j++)
						pos += sprintf(S_bufSKATO_sub2 + pos, ",%g", Ra_cov[i][j]);
				}

				// Clear
				LOOP (i, X.N_rho) { DEALLOC(Xa_lius[i]); }
				DEALLOC(Xa_lius);
			}
			else {
				if (X.Cp_naTest == NULL) {
					LOGoutput("List of genes including only one variant "
						"is exported to [%s.gene.fail.lst], since SKAT-o "
						"test cannot be performed on the genes\n",
						OPT_STRING(out));
/**/				X.Cp_naTest = cExporter::summon("gene.fail.lst");
//					Cp_naTest->put("GENE	SKATo	PEDCMC\n");
				}
				B_SKAToFailed = 1;

				/* Do not perform SKAT-o when gs == 1 */
				strcpy(S_bufSKATO, "	NA	NA	NA");
				// 				for (i=1 ; i<(N_rho-1) ; i++)
				// 					strcat(S_bufSKATO, "	NA	NA	NA");
				strcat(S_bufSKATO, "	NA	NA	NA	NA	NA");
			}
		}
	} else S_bufSKATO[0] = '\0';

	sseFree(Ra_xxG);
	sseFree(Ra_wGScurr);
	M_genes.rem();
	if (Ra_Wg_t != NULL) {
		for (wsUint l=0 ; l<N_curSamp ; l++)
			sseFree(Ra_Wg_t[l*N_pheno]);
		DEALLOC(Ra_Wg_t);
	}
#ifdef GSvalidate
	sprintf(fn, "gs.CMC.%s", S_gs);
	V_pldCMC.file(fn);
//	exportMatrix(fn, &Ra_pldCMC, 1, N_sampGS);
	sprintf(fn, "gs.collapsing.%s", S_gs);
//	exportMatrix(fn, &Ra_pldClp, 1, N_sampGS);
	V_pldClp.file(fn);
#endif

	/* Perform CMC */
 	cMask	V_maskRep(V_mask, N_szVec==N_anaSamp ? 1 : N_pheno);
 	R_pvCMC = _getCMC(V_maskRep, X.V_hat, X.R_vars, X.V_aV, V_pldCMC, &R_statCMC, N_szVec);
	if (isMissingReal(R_pvCMC) || R_statCMC == numeric_limits<wsReal>::infinity())
		sprintf(S_bufCMC, "	NA	NA");
	else
		sprintf(S_bufCMC, "	%g	%g", R_statCMC, R_pvCMC);

	/* Perform Collapsing */
	wsReal R_aVyClp;
	R_pvClp = _getCollapsing(V_maskRep, X.R_vars, X.V_aV, V_pldClp, &R_statClp,
		N_szVec, &R_aVyClp);
	R_pvClp = _getCollapsing(V_maskRep, X.R_vars, X.V_aV, V_pldClp, &R_statClp,
		N_szVec, &R_aVyClp);
	if (isMissingReal(R_pvClp))
		sprintf(S_bufClp, "	NA	NA");
	else
		sprintf(S_bufClp, "	%g	%g", R_statClp, R_pvClp);

	/* Perform PEDCMC */
	if (OPT_ENABLED(pedcmc)) {
		/* Counting */
//			if ((N_common+N_rare) > 1) {
			wsUint N_typeTest = 0;

			/* V_mask in here should be union of V_mask and V_maskPheno */
			cMask V_uniMask(N_anaSamp, 1);
			wsReal *Ra_mAll = V_mask.get();
			wsReal *Ra_mPhe = X.V_availPhenoMask.get();
			wsReal *Ra_uni = V_uniMask.get();
			LOOP (i, N_anaSamp) {
				if (Ra_mAll[i] == W0 ||
					Ra_mPhe[i] == W0)
					Ra_uni[i] = W0;
			}

			R_statPEDCMC	= W0;
			R_pvPEDCMC		= _getPEDCMC(Cp_IO, N_samp, X.Mp_phi, Vp_yGSunAdj,
				S_gs, Xv_set, Ra_PEDCMCca,
				Ra_PEDCMCct, &R_statPEDCMC, N_curCtrl, N_curCase,
				N_common, N_rare, V_uniMask, Ba_filtAll, &N_typeTest,
				&N_PEDCMCdf);

			if (isMissingReal(R_statPEDCMC))
				sprintf(S_bufPEDCMC, "	%d	%d	%d	%d	NA	%d	%d	%d	NA	NA",
					N_maCase, N_maCtrl, N_PEDCMCdf,
					N_miCase, N_miCtrl, N_common, N_rare);
			else
				sprintf(S_bufPEDCMC, "	%d	%d	%d	%d	%s	%d	%d	%d	%g	%g",
					N_maCase, N_maCtrl,
					N_miCase, N_miCtrl, N_typeTest==2?"non2x2":(N_typeTest?"2x2fisher":"2x2chisq"),
					N_PEDCMCdf,
					N_common, N_rare,
					R_statPEDCMC, R_pvPEDCMC);
// 		} else {
// 			if (Cp_naTest == NULL) {
// 				LOG("PEDCMC failed due to there is single variant within"
// 					" given gene-set, list of those gene-set will reported."
// 					" However, their result will be equivalent to p-value"
// 					" of collapsing test, please refer that.\n");
// 				Cp_naTest	= new cExporter("gs.fail");
// 				Cp_naTest->put("GENE	SKATo	PEDCMC\n");
// 			}
// 			B_PEDCMCFailed = 1;
// 			sprintf(S_bufPEDCMC, "	NA	NA	%d	%d	NA	NA", N_common, N_rare);
// 		}

		sseUnmat(Ra_PEDCMCca, N_col);
		sseUnmat(Ra_PEDCMCct, N_col);
	} else S_bufPEDCMC[0] = '\0';

	if (B_SKAToFailed/* || B_PEDCMCFailed*/) {
// 		Cp_naTest->fmt("%s	%s	%s\n", S_gs,
// 			B_SKAToFailed?"FAILED":"SUCCESS",
// 			B_PEDCMCFailed?"FAILED":"SUCCESS");
		X.Cp_naTest->fmt("%s\n", S_gs);
	}

	if (OPT_ENABLED(skat))
		sprintf(S_bufSKAT, "	%g	%g	%d", R_statSKAT, R_pvSKAT,
			B_isLiuApplied);

	/* Print out statistics */
	char	*Sp_chr = NULL;
	char	*S_na = (char *)"NA";
	if (Nv_vrtIdx.size())
		Sp_chr = (char *)getChrName2(Xv_vrt[Nv_vrtIdx.at(0)].chr);
	else
		Sp_chr = S_na;

	if (OPT_ENABLED(time))
		strcpy(S_bufTime, c.getReadable());
			
	if (OPT_ENABLED(farvatx))
		X.Cp_exporter->fmt("%s	%s	"
			"%d	%d	%d	%d	%d	%d	%d	%d"
			"%s%s\n",
			Sp_chr, S_gs,
			N_curSamp, (int)Nv_vrtIdx.size(), (int)Xv_set.size(),
			N_singleton, N_mac, N_imp, N_start, N_end,
			S_bufTime, S_bufFARVATx);
	else
		X.Cp_exporter->fmt("%s	%s	"
			"%d	%d	%d	%d	%d	%d	%d	%d"
			"%s%s%s	%g%s%s%s%s%s%s%s%s%s%s%s%s\n",
			Sp_chr, S_gs,
			N_curSamp, (int)Nv_vrtIdx.size(), (int)Xv_set.size(),
			N_singleton, N_mac, N_imp, N_start, N_end,
			S_bufTime,
//			c.getReadable(),
			S_bufCMC, S_bufClp, R_aVyClp, S_bufKBAC, S_bufWsum, S_bufAsum,
			S_bufVTfam, S_bufVT, S_bufGGemma, S_bufPgene, S_bufSKAT,
			S_bufPEDCMC, S_bufSKATO, S_bufSKATO_sub2, S_bufFARVATx);
	delete [] S_bufSKATO_sub2;

	sseFree(Ba_filtAll);
	DEALLOC(Ra_mafBins);
	if (Ra_VTsum1)
		sseUnmat(Ra_VTsum1, N_gs);
	sseFree(Ra_VTsum2);
//		sseFree(Ra_maskTest);
}

void cGSTestAnalysis::run()
{
	/*
	 *
	 * Checking conditions
	 *
	 */
	if (OPT_ENABLED(pedcmc)) {
		if (!IS_ASSIGNED(raremaf)) {
			OPTION().assign("raremaf", OPTION().getDefVal("raremaf"));
			OPTION().FORCE_OPT_REAL(raremaf);
		}

		if (OPT_ENABLED(imputepheno))
			LOG("--imputepheno is enabled, but --pedcmc cannot include "
				"missing individual, they will be excluded\n");
	}

	if (N_anaSamp == 0) {
		LOG("No available sample left, can't do this analysis\n");
		return;
	}

	/*
	 *
	 * Commonly used variables
	 * 
	 */
	wsUintCst	N_sample	= Cp_IO->sizeSample();	///< Number of samples in final dataset
	vVariant&	Xv_vrt		= Cp_IO->getVariant();		///< Variant info vector
	xMaf*		Xp_maf		= Cp_IO->getMAF();
	wsRealCst*	Ra_pheno	= Cp_IO->getPheno();
	wsReal**	Ra_phenos	= Cp_IO->getPhenos();
	mGeneDef&	Xa_gs		= Cp_anaGsm->getGeneDef();	///< Given gene-sets
	wsUint		N_pheno		= IS_ASSIGNED(longitudinal) ? Cp_IO->sizePheno() : 1;

	wsUint		i, _i, j;
/**/cExporter*	Cp_exporter = cExporter::summon("gene.res");						///< Output exporter
	LOGoutput("Result of gene-level test is exported to [%s.gene.res]\n",
		OPT_STRING(out)); /* CONFIRMED 140109 */

	wsReal		R_IIa		= W0;	///< Sum of adjusted phenotype
//	wsReal		**Ra_HH		= NULL;

	/* PEDCMC part */
	wsUint		N_ctrl		= 0;	///< Number of cases in dataset
	wsUint		N_case		= 0;	///< Number of controls in dataset
//	char		OPT_ENABLED(pedcmc)	= 0;	///< (==1) perform PEDCMC or not

	if (OPT_ENABLED(pedcmc) || OPT_ENABLED(wsum) || OPT_ENABLED(kbac)) {
		/* Count the number of included individual */
		for (i=_i=0 ; i<N_sample ; i++) {
			if (Ba_misPheno[i]) continue;

			if (Ra_pheno[i]==WISARD_AFFECTED) N_case++;
			else if (Ra_pheno[i]==WISARD_UNAFFECTED) N_ctrl++;
			_i++;
		}

		/* Sum of case and control should be equal to the number of samples
		 * should be included to analysis */
		if ((N_case+N_ctrl) != N_availPhenoSamp) {
			LOGwarn("Counted sample size [case %d+control %d=%d] "
				"does not match with expected number of samples [%d]\n",
				N_case, N_ctrl, N_availPhenoSamp);
			LOG("Thus, PEDCMC will not performed\n");
		}// else if (OPT_ENABLED(pedcmc))
//			OPT_ENABLED(pedcmc) = 1;

		if (N_case == 0 && N_ctrl == 0)
			halt("No control nor case sample is available, can't test!");
	}

	/* QLSuniv part */
	// [1]
	// IIa     <- sum(YY)
	R_IIa = Vp_yGS->sum();

	// [#pheno*#pheno] = [#pheno*#samp]%*%[#samp*#samp]%*%[#samp*#pheno] - [???]
	// vars    <- (t(YY)%*%phi%*%YY)[1,1] - IIa^2*IIphiInvII_Inv
	if (IS_ASSIGNED(longitudinal) || OPT_ENABLED(indep) || OPT_ENABLED(logistic))
		R_vars = Vp_yGS->qf(*(cIdtMatrix *)Mp_phi) - SQR(R_IIa)*R_1phiInv1_inv;
// 	else if (OPT_ENABLED(logistic))
// 		R_vars = Vp_yGS->qf(*(cDiagMatrix *)Mp_phi) - SQR(R_IIa)*R_1phiInv1_inv;
	else
		R_vars = Vp_yGS->qf(*(cSymMatrix *)Mp_phi) - SQR(R_IIa)*R_1phiInv1_inv;
//	R_vars		= sseVpMpV(Ra_yGS, Ra_phiGS, N_sampGS) -
//		R_IIa*R_IIa*R_1phiInv1_inv;

// 	wsReal *Rp_vars = &R_vars;
// 	exportMatrix("gs.vars", &Rp_vars, 1, 1);
	LOG("qlsUinv vars = %g\n", R_vars);

	/* SKAT part */
	// [#samp*#samp] = [#samp*#samp] - [1]
	// HH      <- akinship - matrix(IIphiInvII_Inv,totind,totind)
	cMatrix&	M_HH = *Mp_phi - R_1phiInv1_inv;
	cSymMatrix	M_HH_X;
	if (OPT_ENABLED(farvatx)) {
		cVector II(W1, N_anaSamp);
		FOREACH (vInt_it, Nv_idxFemale, i) II.get()[*i] = W2;

		// HH             <- akinship - II%*%IIphiInvII_Inv%*%t(II)
		wsSym Ra_HH_X = sseSymMatP(N_anaSamp, Mp_phi->get());
		LOOP(i, N_anaSamp) {
			wsReal R_mul1 = II.get()[i];
			LOOP(j, i) Ra_HH_X[i][j] -= R_mul1*II.get()[j] * R_1phiInv1_inv_Xall;
			Ra_HH_X[i][i] -= SQR(R_mul1) * R_1phiInv1_inv_Xall;
		}
		M_HH_X.init(Ra_HH_X, N_anaSamp);
	}

//	Ra_HH		= sseMatrix(N_sampGS, N_sampGS); /* 120917 */
//	sseMsS(Ra_phiGS, R_1phiInv1_inv, Ra_HH, N_sampGS, N_sampGS);
#ifdef GSvalidate
	M_HH.file("gs.HH");
#endif

	// [#pheno*#pheno] = [#pheno*#samp]%*%[#samp*#samp]%*%[#samp*#pheno]
	// pppMult <- (t(YY)%*%HH%*%YY)[1,1]
//	if (IS_ASSIGNED(longitudinal))
//		R_pppMult	= Vp_yGS->qf(*(cBlkMatrix *)(&M_HH));
//	else
// 	if (OPT_ENABLED(indep) || OPT_ENABLED(logistic))
// 		R_pppMult	= Vp_yGS->qf(*(cDiagMatrix *)(&M_HH));
// 	else if (OPT_ENABLED(logistic))
// 		R_pppMult	= Vp_yGS->qf(*(cDiagMatrix *)(&M_HH));
//	else
		R_pppMult	= Vp_yGS->qf(*(cSymMatrix *)(&M_HH));
//	R_pppMult	= sseVpMpV(Ra_yGS, Ra_HH, N_sampGS);
	R_C			= R_pppMult;
//	sseUnmat(Ra_HH, N_sampGS);
#ifdef GSvalidate
	wsReal *Rp_pppMult = &R_pppMult;
	exportMatrix("gs.pppMult", &Rp_pppMult, 1, 1);
#endif

#if 0
	/* --metafarvat */
	if (OPT_ENABLED(makefarvat)) {
#ifdef GSvalidate
		V_metaFarvat.file("gs.mf.before");
#endif
//		V_metaFarvat /= sqrt(R_C);
#ifdef GSvalidate
		V_metaFarvat.file("gs.mf.after");
#endif
	}
#endif

	/* SKAT-o part */
// 	wsReal	**Ra_xxM	= NULL;
// 	wsReal	*Ra_xx		= NULL;
	wsUint	N_rho		= 0;	/* Number of rhos to evaluate */
	wsVec	Ra_rhos		= NULL;	/* Array of rhos */
	wsVec	Ra_ySKATo	= NULL;
	cVector	V_xx(N_anaSamp*N_pheno);

	cSymMatrix	C_Psym;

	if (OPT_ENABLED(skato) || OPT_ENABLED(famvt)) {
		/*
		 *
		 * Set ��
		 *
		 */
		if (OPT_ENABLED(skato))
			Ra_rhos = getOptimalRhos(N_rho);

		/* HH2     <- diag(totind) - matrix(IIphiInvII_Inv, totind,totind)%*%akinshipInv
		 *  - Since akinshipInv is symmetric, matrix(IIphiInvII_Inv, totind,totind)%*%akinshipInv
		 *   is equivalent to i,j = sum(phiInv_j+) */

		/* xxM_(i) is same to y(i) - \Sigma_(k=1)^n [ \Phi^-1_ik * 1phiInv1inv ] * sumY */
		wsReal	*Rp_xx		= V_xx.get();
		wsReal	*Rp_yGS		= Vp_yGS->get();
		wsReal	R_yGSsum	= Vp_yGS->sum();
		wsReal	*Rp_phiInv1	= V_phiInv1.get();
		wsUint	N_szY		= N_anaSamp * N_pheno;
		for (i=0 ; i<N_szY ; i++)
			Rp_xx[i] = Rp_yGS[i] - Rp_phiInv1[i] * R_yGSsum / R_1phiInv1;

		/* For continuous trait */
		if (B_isCont) {
			wsUint		N_inpCol	= Cp_IO->sizeCovar()+1;
			cStdMatrix	C_Xt(N_inpCol, N_anaSamp); ///< [#cov*#tsamp]
			wsMat		Ra_Xt		= C_Xt.get();
			wsMat		Ra_cov		= Cp_IO->getCovariates();

			/* Build X=intercept+covariates */
			for (i=_i=0 ; i<N_sample ; i++) {
				if (Ba_misPheno[i]) continue;

				Ra_Xt[0][_i] = W1;
				for (j=1 ; j<N_inpCol ; j++)
					Ra_Xt[j][_i] = Ra_cov[j-1][i];
				_i++;
			}
#ifdef GSvalidate
			exportMatrix("skato.Xt", Ra_Xt, N_inpCol, N_anaSamp);
//			exportMatrix("skato.sig", &Ra_est, 1, 2);
#endif
			/* (--longitudinal)
			 *     [ ��     ]
			 * V = [  ...  ], �� repeats N times
			 *     [     �� ]
			 *
			 * (otherwise)
			 *
			 * V = sig2*I + sig2g*��
			 *
			 * P = Vi - Vi X (X' Vi X)^-1 X' Vi
			 *
			 * Psq * Psq = P
			 */

			if (IS_ASSIGNED(longitudinal)) {
				cTimer		t;
				t.start();
				cSymMatrix	C_varPheno	= cSymMatrix(Ra_phiGS, N_szPhi, 1);
				//cSymMatrix	C_varPinvSq;
				cSymMatrix	C_varPinv	= C_varPheno.einv();//&C_varPinvSq);
				// Get I - X(X'X)^-1X'
				cSymMatrix	C_proj		= C_Xt.proj();
				cSymMatrix	C_pInvSq;
//				cSymMatrix	&C_pSq		= C_proj.sqrt();
				// Get (I - X(X'X)^-1X')^-1 %x% Vinv
#if 1
				cStdMatrix	M_projEV;
				cVector		V_projEV = C_proj.eigen(M_projEV);
				cStdMatrix	M_piEV;
				cVector		V_piEV = C_varPinv.eigen(M_piEV);
//				C_proj.file("proj");
//				C_varPinv.file("pi");
//				M_projEV.file("EXproj");
//				V_projEV.file("EVproj");
//				M_piEV.file("EXpi");
//				V_piEV.file("EVpi");

				cStdMatrix	M_krEV	= M_projEV.krnk(M_piEV);
//				M_krEV.file("EXkr");
				cVector		V_krEV	= V_projEV.krnk(V_piEV);
//				V_krEV.file("EVkr");
				Mp_P	= &(M_krEV.MMt(V_krEV));
//				Mp_P->file("P");
//				Cp_Psq	= &(C_pSq.krnk(C_varPinvSq));
				Mp_Psq	= &(((cSymMatrix *)Mp_P)->sqrt(M_krEV.get(),
					V_krEV.get(), V_krEV.size()));
//				Cp_Psq->file("Psq");
//				cSymMatrix M_P2 = ((cSymMatrix *)Cp_Psq)->Mt();
//				LOG("Diff %g\n", compSSfrobenius(M_P2.get(), M_P2.row(), Mp_P->get(), Mp_P->row()));
//				delete &C_pSq;
#else
				Mp_P	= &(C_proj.krnk(C_varPinv));
				//				Cp_Psq	= &(C_pSq.krnk(C_varPinvSq));
				Mp_Psq	= &(Mp_P->sqrt());
#endif
				LOG("Longitudinal precomputation [%s] taken", t.getReadable());
			} else {
				if (OPT_ENABLED(logistic)) {
					/* For --logistic, V = diag(pHat) */
//					cDiagMatrix V(sqrt(R_logisticS), N_anaSamp);
					cDiagMatrix V(N_anaSamp, Ra_logisticV);
					Mp_P = &(V.clone());
				} else if (OPT_ENABLED(indep)) {
					/* sig_g^2 is ignored, but sig^2 is remained */
					wsRealCst *Ra_est	= Cp_anaEmai->getEstimates(0);
					cDiagMatrix V(Ra_est[0], N_anaSamp);
					Mp_P = &(V.inv());
				} else {
					/* V = sig2*I + sig2g*�� */
					wsRealCst		*Ra_est	= Cp_anaEmai->getEstimates(0);
					cSymMatrix	&V		= *((cSymMatrix *)Mp_phi)*Ra_est[1];
					wsReal		**Ra_V	= V.get();
					for (i=0 ; i<N_anaSamp ; i++)
						Ra_V[i][i] += Ra_est[0];
					Mp_P = &(V.inv());
					V.rem();
				}

				cMatrix &C_Vinv = *Mp_P;

				// ViX  <- Vinv%*%X
				//      <- Vinv %**% t(X)
				cStdMatrix	C_ViX = OPT_ENABLED(indep) || OPT_ENABLED(logistic) ?
					((cDiagMatrix *)&C_Vinv)->mt(C_Xt) :
					((cSymMatrix *)&C_Vinv)->mt(C_Xt);

				// XVIXinv <- solve(t(X)%*%Vinv%*%X)
				//         <- solve(t(X) %*% Vinv %**% t(X))
				cSymMatrix C_XtViX = OPT_ENABLED(indep) || OPT_ENABLED(logistic) ?
					C_Xt.MMt(*((cDiagMatrix *)&C_Vinv)) :
					C_Xt.MMt(*((cSymMatrix *)&C_Vinv));

				cSymMatrix &C_XtViX_inv = C_XtViX.inv();

				// [NP][PN][NP] = NP(N+P)
				// Vdec <- ViX %*%XVIXinv%*%t(X)%*%Vinv == t(ViX)
				//      <- ViX %*% XVIXinv %**% ViX
				cSymMatrix C_Videc = C_ViX.MMt(C_XtViX_inv);
				delete &C_XtViX_inv;
				C_ViX.rem();

				if (OPT_ENABLED(logistic) || OPT_ENABLED(indep)) {
					C_Psym = *((cDiagMatrix *)Mp_P) - C_Videc;
					Mp_P->rem();
					Mp_P = &C_Psym;
#ifdef GSvalidate
					Mp_P->file("skato.P");
#endif
					// P0.5 <- sqMat(P)
					Mp_Psq = &(*((cDiagMatrix *)Mp_P) * 0.5);
				}
				else {
					C_Psym = *((cSymMatrix *)Mp_P) - C_Videc;
					Mp_P->rem();
					Mp_P = &C_Psym;
#ifdef GSvalidate
					Mp_P->file("skato.P");
#endif
					// P0.5 <- sqMat(P)
					Mp_Psq = &(Mp_P->sqrt());
				}
			}
			
			// YY    <- P%*%pheMat
			wsUint N_szVec = N_anaSamp*N_pheno;
			wsUint N_szPhen = N_pheno;
			if (OPT_ENABLED(logistic)) {
				/* 140111 with --logistic, (y-p) = YY */
				Ra_ySKATo = sseVector(N_szVec);
				memcpy(Ra_ySKATo, Vp_yGS->get(), sizeof(wsReal)*N_szVec);
			} else {
				cVector V_yGSunAdj(N_szVec);

				vSampPtr& Xa_samp = getIO()->getSample();
				wsReal *Ra_yUnadjGS = V_yGSunAdj.get();
				//		sseMalloc(Ra_yUnadjGS, wsReal, N_sampGS);
				wsUint _i = 0;
				i = 0;
				FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
					if (Ba_misPheno[i]) continue;

					for (wsUint j=0 ; j<N_szPhen ; j++)
						Ra_yUnadjGS[_i+j] = Ra_phenos[j][i];
					_i += N_szPhen;
				}

				cVector V_yy = *Mp_P * V_yGSunAdj;
				V_yy.setDontDealloc();
				Ra_ySKATo = V_yy.get();
			}
#ifdef GSvalidate
			exportMatrix("skato.YY", &Ra_ySKATo, 1, Mp_P->row());
#endif
#ifdef GSvalidate
			exportMatrix("skato.Psq", Mp_Psq->get(), Mp_P->row(), Mp_P->row());
#endif
		}
	}

	/*
	 * 
	 * Preparing header
	 *
	 */
	char	S_bufCMC[128]		= { 0, };
	char	S_bufTime[128]		= { 0, };
	char	S_bufClp[128]		= { 0, };
	char	S_bufVT[128]		= { 0, };
	char	S_bufVTfam[128]		= { 0, };
	char	S_bufGGemma[128]	= { 0, };
	char	S_bufAsum[128]		= { 0, };
	char	S_bufPgene[128]		= { 0, };
	char	S_bufSKAT[128]		= { 0, };
	char	S_bufFARVATx[4096]	= { 0, };
	char	S_bufSKATO[1024]	= { 0, };
	char	S_bufWsum[128]		= { 0, };
	char	S_bufPEDCMC[256]	= { 0, };
	char	S_bufKBAC[128]		= { 0, };

	/* FIXME: Bugfix needed Adaptive sum */
//	if (OPT_ENABLED(asum)) {
//		LOG("Adaptive-sum test will performed\n");
//		sprintf(S_bufAsum, "	P_ASUM");
//	}

	/* KBAC */
	wsReal R_alphaKB = WISARD_NAN;
	if (OPT_ENABLED(kbac)) {
		sprintf(S_bufKBAC, "	P_KBAC");
		/* Force assignment of KBAC alpha to default value */
		if (!IS_ASSIGNED(kbacalpha)) {
			OPTION().assign("kbacalpha", OPTION().getDefVal("kbacalpha"));
			OPTION().FORCE_OPT_REAL(kbacalpha);
		}
		/* Set the alpha level */
		R_alphaKB = OPT_REAL(kbacalpha);

		/* Kernel check */
		if (!IS_ASSIGNED(kbackernel)) {
			OPTION().assign("kbackernel", OPTION().getDefVal("kbackernel"));
			OPTION().FORCE_OPT_STRING(kbackernel);
		}
		wsStrCst S_kernel = OPT_STRING(kbackernel);
		if (stricmp(S_kernel, "hypergeometric"))
			halt("Only 'hypergeometric' kernel is currently supported");

		/* Print */
		LOG("KBAC test will performed [%s-sided, alpha=%g, kernel=%s%s]\n", OPT_ENABLED(kbac2side)?
			"two":"one", R_alphaKB, S_kernel, R_alphaKB < W1 ? ", adaptive":"");
	}

	/* Wsum */
	wsUint	N_perm		= IS_ASSIGNED(nperm) ? OPT_NUMBER(nperm) : 1000;
	char**	Na_permPhe	= NULL; ///< [#actSamp * #perm]
	if (OPT_ENABLED(wsum)) {
		LOG("Weighted-sum test will performed\n");

		sprintf(S_bufWsum, "	STAT_WSUM	P_WSUM");
		wsAlloc(Na_permPhe, char*, N_sample);
		for (i=0 ; i<N_sample ; i++)
			wsAlloc(Na_permPhe[i], char, N_perm);
	}
	
	/* VT */
	if (OPT_ENABLED(vt)) {
		LOG("Variable-threshold test will performed\n");
		sprintf(S_bufVT, "	STAT_VT	VTthr	P_VT");
	}

	/* famVT */
	if (OPT_ENABLED(famvt)) {
		LOG("Variable-threshold test for family will performed\n");
		sprintf(S_bufVTfam, "	P_FAMVT");
	}

	/* Gene GEMMA */
	cFemmaAnalysis*	Cp_gemma	= NULL;
	wsReal			R_gemmal0	= WISARD_NAN;
	cStdMatrix		M_gemmaXt;
	if (OPT_ENABLED(ggemma)) {
		/* Get reduced covariates matrix */
		char*	Ba_filt	= Cp_IO->getPheCovMissing();
		wsMat	Ra_DM	= getIO()->getFiltNullDesignMat(Ba_filt, 1, N_anaSamp, 1);
		Cp_gemma	= new cFemmaAnalysis(getIO(), Cp_anaCorr);
		R_gemmal0	= Cp_gemma->get0().R_logL;
		M_gemmaXt.init(getIO()->sizeCovar()+2, N_anaSamp, Ra_DM);
		/* Set 1 */
		M_gemmaXt.setRow(0, 1);

		LOG("Gene-level GEMMA test will performed\n");
		sprintf(S_bufGGemma, "	P_GGEMMA_LRT	P_GGEMMA_SCORE	P_GGEMMA_WALD");
	}

	/* PEDCMC */
	if (OPT_ENABLED(pedcmc))
		sprintf(S_bufPEDCMC, "	PEDcaMa	PEDctMa	PEDcaMi	PEDctMi	PEDtest	"
			"PEDdf	PEDcommon	PEDrare	STAT_PEDCMC	P_PEDCMC");

	/* pedgene */
	wsReal	R_cFtr = WISARD_NAN;
	cVector	V_resid;
	if (OPT_ENABLED(pedgene)) {
		V_resid	= Vp_yGS->clone();

		R_cFtr	= quadfactor(Mp_phi->get(), N_anaSamp, 0, V_resid.get());
		sprintf(S_bufPgene, "	STAT_BDPG	P_BDPG	STAT_KNPG	P_KNPG");
	}

	/* SKAT */
	if (OPT_ENABLED(skat))
		sprintf(S_bufSKAT, "	STAT_SKAT	P_SKAT	SKATliu");

	/* FARVATX */
	if (OPT_ENABLED(farvatx)) {
		wsUint N_dd = OPT_NUMBER(farvatxndiv);
		wsReal R_span = W1 / (wsReal)(N_dd - 1);
		wsReal R_div = W0;
		char*	Sp = (char *)S_bufFARVATx;
		for (wsUint i=0 ; i<N_dd ; i++,R_div+=R_span)
			Sp += sprintf(Sp, "P_FARVATXB_%g	P_FARVATXC_%g	P_FARVATXO_%g	",
				R_div, R_div, R_div);
		Sp += sprintf(Sp, "P_FARVAT_R1	P_FARVAT_R2");
	}

	/* SKAT-o */
	if (OPT_ENABLED(skato)) {
//		char S_bufSKATO_sub[256];
		if (OPT_ENABLED(farvat))
			sprintf(S_bufSKATO, "	P_FARVATC	STAT_FARVATC	cAlphaDf");
		else
			sprintf(S_bufSKATO, "	P_CALPHA	STAT_CALPHA	cAlphaDf");
// 		for (i=1 ; i<(N_rho-1) ; i++) {
// 			sprintf(S_bufSKATO_sub, "	Rho_%g_pv	Rho_%g_stat	Rho_%g_df"
// 				"	Rho_%g_ncp", Ra_rhos[i], Ra_rhos[i], Ra_rhos[i], Ra_rhos[i]);
// 			strcat(S_bufSKATO, S_bufSKATO_sub);
// 		}
		if (OPT_ENABLED(farvat))
			strcat(S_bufSKATO, "	P_FARVATB	STAT_FARVATB	Burden_df	"
				"RHO_FARVATO	P_FARVATO");
		else
			strcat(S_bufSKATO, "	P_BURDEN	STAT_BURDEN	Burden_df	"
				"RHO_SKATO	P_SKATO");
		if (OPT_ENABLED(makefarvat))
//			strcat(S_bufSKATO, "	FV_VARIANTS	FV_MACS	FV_NALLELES	FV_STATS	FV_COVAR");
			strcat(S_bufSKATO, "    FV_VARIANTS     FV_MACS	FV_STATS        FV_COVAR");
	}
	if (OPT_ENABLED(time))
		strcpy(S_bufTime, "TIME	");
	/* Print header */
	if (OPT_ENABLED(farvatx))
		Cp_exporter->fmt("CHR	GENE	NSAMP	SZSET	NVARIANT	NSINGLETON	MAC"
			"	NIMP	START	END	%s%s\n",
			S_bufTime, S_bufFARVATx);
	else
		Cp_exporter->fmt("CHR	GENE	NSAMP	SZSET	NVARIANT	NSINGLETON	MAC"
			"	NIMP	START	END	%sSTAT_CMC	P_CMC	STAT_CLP"
			"	P_CLP	ClpAVY%s%s%s%s%s%s%s%s%s%s\n",
			S_bufTime,
			S_bufKBAC, S_bufWsum, S_bufAsum, S_bufVTfam, S_bufVT,
			S_bufGGemma, S_bufPgene, S_bufSKAT, S_bufPEDCMC, S_bufSKATO);

	wsUint		I_gene		= 1;	/* Gene loop variable */

	/* Get data */
	char**		Na_geno		= Cp_IO->getGenotype();
	/* IID vector */
	vSampPtr&	Xv_smp		= Cp_IO->getSample();
	vStr		Xv_smpbuf; {
		wsUint i = 0;
		FOREACHDO (vSampPtr_it, Xv_smp, I, i++)
			if (Ba_misPheno && !Ba_misPheno[i])
				Xv_smpbuf.push_back((*I)->S_IID);
	}

/**/cExporter	*Cp_naTest	= NULL;
	FOREACHDO (mGeneDef_it, Xa_gs, it, I_gene++) {
		xGeneTestParam X = {
			// Common
			Cp_IO,
			Xv_smpbuf, R_C,
			V_hat, R_vars,
			B_isCont,
			Na_permPhe,
			V_aV,
			Vp_yGS,
			Vp_wGS,
			Ba_misPheno,
			// Output
			Cp_exporter, Cp_naTest,
			// GGEMMA
			M_gemmaXt, R_gemmal0,
			// SKAT
			M_HH_X, R_pppMult,
			// Pedgene
			V_resid, R_cFtr,
			// PEDCMC
			V_availPhenoMask,
			// FARVATx
			V_phiInv1_X, R_1phiInv1_inv_X,
			Nv_idxFemale, Nv_idxFnder,
			R_1phiInv1_inv_Xall,
			M_I_pi1pi1i,
			// MFARVAT
			V_1TPA, R_sumTPAPT,
			M_TPAPT, M_Tt, M_TPA,
			// metaFARVAT
			V_metaFarvat,
			// KBAC
			R_alphaKB,
			// VT
			R_pheVTmean, Ra_pheVT,
			// SKAT-o
			Ra_ySKATo, N_rho, Ra_rhos,
			Mp_P, Mp_Psq, Mp_phi, Mp_phiInv,
		};

		do_geneTest(I_gene, Xa_gs.size(),
			Cp_IO, Ba_misPheno, N_anaSamp, N_perm,
			N_pheno, N_case, N_ctrl,
			it->first.c_str(), it->second, Xp_maf,
			Vp_yGSunAdj, Ra_pheVT, R_pheVTmean,
			Cp_gemma, V_xx, X);
	}

	delete& M_HH;

	/* SKAT-o part */
	if (OPT_ENABLED(skato))
		sseFree(Ra_rhos);

	if (Cp_naTest)
		delete Cp_naTest;

	delete Cp_exporter;
	sseFree(Ra_ySKATo);
// 	sseFree(Ra_pldCMC);
// 	sseFree(Ra_pldClp);
	LOG("%d/%d geneset processed...\n", Xa_gs.size(), Xa_gs.size());

	if (OPT_ENABLED(wsum)) {
		for (i=0 ; i<N_sample ; i++)
			DEALLOC(Na_permPhe[i]);
		DEALLOC(Na_permPhe);
	}
	delete Vp_wGS;
//	sseFree(Ra_wGS);

	if (Ra_phiGS != Ra_origPhi)
		sseUnmat(Ra_phiGS, N_anaSamp);
	if (IS_ASSIGNED(longitudinal)) {
		sseUnmat(Ra_origPhi, N_pheno);
//	if (Vp_yGSunAdj != NULL)
//		delete [] Vp_yGSunAdj;
//	Mp_phi->rem();
//	if (Cp_P) delete Cp_P;
	} else if (OPT_ENABLED(indep)) {
		sseUnmat(Ra_origPhi, N_anaSamp);
	}
	if (Mp_P && Mp_P->get())
		Mp_P->rem();
//	Cp_P = NULL;
}

cStdMatrix rho1(wsReal R_rho, cStdMatrix& M_r1, cStdMatrix& M_r2)
{
	wsReal	R_1rho	= W1 - R_rho;
	wsUint	N		= M_r1.row();
	wsMat	Ra_ret	= sseMatrix(N, N);
	wsMat	Ra_r1	= M_r1.get();
	wsMat	Ra_r2	= M_r2.get();
	wsUint	N_med	= 0;
#ifdef USE_SSE
	N_med = getMed(N);
	sse_t	sse_r	= sseSet(R_rho);
	sse_t	sse_1r	= sseSet(R_1rho);
#endif

	for (wsUint i=0 ; i<N ; i++) {
		wsReal*	Ra_r	= Ra_ret[i];
		wsReal*	Ra_1	= Ra_r1[i];
		wsReal*	Ra_2	= Ra_r2[i];
#ifdef USE_SSE
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t*	sse_1	= (sse_t *)(Ra_1 + j);
			sse_t*	sse_2	= (sse_t *)(Ra_2 + j);
			sse_t*	sse_R	= (sse_t *)(Ra_r + j);
			*sse_R = sseAdd(sseMul(*sse_1, sse_r), sseMul(*sse_2, sse_1r));
		}
#endif
		for (wsUint j=N_med ; j<N ; j++)
			Ra_r[j] = Ra_1[j]*R_rho + Ra_2[j]*R_1rho;
	}
	return cStdMatrix(N, N, Ra_ret);
}

cStdMatrix rho2(wsReal R_rho, cStdMatrix& M_r1, cStdMatrix& M_r2, cStdMatrix& M_r3)
{
	wsReal	R_1rho	= W1 - R_rho;
	wsReal	R_r2	= SQR(R_rho);
	wsReal	R_r1r	= R_rho * R_1rho;
	wsReal	R_1r2	= SQR(R_1rho);

	wsUint	N		= M_r1.row();
	wsMat	Ra_ret	= sseMatrix(N, N);
	wsMat	Ra_r1	= M_r1.get();
	wsMat	Ra_r2	= M_r2.get();
	wsMat	Ra_r3	= M_r3.get();
	wsUint	N_med	= 0;
#ifdef USE_SSE
	N_med = getMed(N);
	sse_t	sse_r2	= sseSet(R_r2);
	sse_t	sse_r1r	= sseSet(R_r1r);
	sse_t	sse_1r2	= sseSet(R_1r2);
#endif

	LOOP(i, N) {
		wsReal*	Ra_r	= Ra_ret[i];
		wsReal*	Ra_1	= Ra_r1[i];
		wsReal*	Ra_2	= Ra_r2[i];
		wsReal*	Ra_3	= Ra_r3[i];
#ifdef USE_SSE
		SSE_LOOP(j, N_med) {
			sse_t*	sse_1	= (sse_t *)(Ra_1 + j);
			sse_t*	sse_2	= (sse_t *)(Ra_2 + j);
			sse_t*	sse_3	= (sse_t *)(Ra_3 + j);
			sse_t*	sse_R	= (sse_t *)(Ra_r + j);
			*sse_R = sseAdd(
				sseAdd(sseMul(*sse_1, sse_r2), sseMul(*sse_2, sse_r1r)),
				sseMul(*sse_3, sse_1r2)
			);
		}
#endif
		LOOPV(j, N_med, N)
			Ra_r[j] = Ra_1[j]*R_r2 + Ra_2[j]*R_r1r + Ra_3[j]*R_1r2;
	}
	return cStdMatrix(N, N, Ra_ret);
}


wsMat _getIMV2_p3(cStdMatrix& M_Z, wsReal R_z2)
{
	// Calculate Zm=Z'Z
	cSymMatrix	M_ZtZ	= M_Z.tM();
///**/MAT_t	Ra_Zp		= sseMtM(Ra_Z, N_szLiu, N_gs, Ra_Z, N_szLiu, N_gs);

	// Compute ZZp
	cStdMatrix	M_ZZtZ	= M_Z * M_ZtZ;
///**/MAT_t	Ra_ZZp		= sseMpMt(Ra_Z, N_szLiu, N_gs, Ra_Zp, N_gs, N_gs);
	// ZZp1 == Rowsum of ZZp
	cVector		V_ZZtZ1	= M_ZZtZ.sumR();
///**/wsReal	*Ra_ZZp1	= sseMsumR(Ra_ZZp, N_szLiu, N_gs);

	// 1'Zp == Zp1 == Rowsum of Zp since Zp is sym.
	cVector		V_ZtZ1	= M_ZtZ.sumR();
///**/wsReal *Ra_Zp1		= sseMsumR(Ra_Zp, N_gs);
	/* 1'ZmZm1 */
	wsReal R_1ZpZp1 = V_ZtZ1.ss();
	M_ZtZ.rem();
//	wsReal R_1ZpZp1 = sseVV(Ra_Zp1, N_gs);
//	sseUnmat(Ra_Zp, N_gs);

	// Z1 == Rowsum of Z
	cVector	V_Z1	= M_Z.sumR();
///**/wsReal *Ra_Z1 = sseMsumR(Ra_Z, N_szLiu, N_gs);

	// 1st minus term : ZZp1(ZZp1)'
	cSymMatrix	M_ZZtZ11tZtZZt = V_ZZtZ1.tV();
///**/MAT_t Ra_ZZp1_1tZpZt = sseVtV(Ra_ZZp1, N_szLiu);

	// 1'ZmZmZ'
	cVector		V_1tZtZZtZZt	= V_ZtZ1.Mt(M_ZZtZ);
	V_ZtZ1.rem();
//**/MAT_t Ra_1tZpZpZt = sseMpMt(&Ra_Zp1, 1, N_gs, Ra_ZZp, N_szLiu, N_gs);
//	sseFree(Ra_Zp1);

	// 2nd minus term : Z1 %*% 1'ZmZmZ'
	cStdMatrix	M_Z11tZtZZtZZt	= V_Z1.tV(V_1tZtZZtZZt);
///**/wsReal **Ra_Z11tZmZmZt = sseVtV(Ra_Z1, N_szLiu, Ra_1tZpZpZt[0]);

	//                  [                            ]
	// Final minus term [ (ZZp1)(ZZp1)' + Z11'ZpZpZ' ] / (p*p)
	//                  [                            ]
	wsUint N_szLiu = M_Z.row();
	wsUint	N_gs	= M_Z.col();
	sseMaM(M_ZZtZ11tZtZZt.get(), M_Z11tZtZZtZZt.get(), M_ZZtZ11tZtZZt.get(),
		N_szLiu, N_szLiu, W1 / (wsReal)(SQR(N_gs)) / R_z2);
// 	sseMaM(Ra_ZZp1_1tZpZt, Ra_Z11tZmZmZt, Ra_ZZp1_1tZpZt,
// 		N_szLiu, N_szLiu, W1 / (wsReal)(SQR(N_gs)) / R_z2);
	M_Z11tZtZZtZZt.rem();
	V_1tZtZZtZZt.rem();
// 	sseUnmat(Ra_Z11tZmZmZt, N_szLiu);
// 	sseUnmat(Ra_1tZpZpZt, 1);

	// Z1 %*% ZZm1
	cStdMatrix	M_Z1ZZtZ1	= V_Z1.tV(V_ZZtZ1);
	V_Z1.rem();
	V_ZZtZ1.rem();
// /**/wsReal **Ra_Z1ZZm1 = sseVtV(Ra_Z1, N_szLiu, Ra_ZZp1);
// 	sseFree(Ra_Z1);
// 	sseFree(Ra_ZZp1);

	// plus term : Z1ZZm1 * 1ZmZm1 / m^4
	sseMpC(M_Z1ZZtZ1.cget(), R_1ZpZp1/(wsReal)(SQR(N_gs)*SQR(N_gs)) / SQR(R_z2),
		M_Z1ZZtZ1.get(), N_szLiu);

	// Compute ZZmZt
	cStdMatrix	M_ZZtZZt	= M_ZZtZ.Mt(M_Z);
//	MAT_t Ra_ZZmZt = sseMpMt(M_ZZtZ, N_szLiu, N_gs, Ra_Z, N_szLiu, N_gs);
	M_ZZtZ.rem();
//	sseUnmat(Ra_ZZp, N_szLiu);
	// ZZmZt - minus term
	M_ZZtZZt -= M_ZZtZ11tZtZZt;
	M_ZZtZ11tZtZZt.rem();
//	sseMsM(Ra_ZZmZt, Ra_ZZp1_1tZpZt, Ra_ZZmZt, N_szLiu);
	// ZZmZt - minus term + plus term
	M_ZZtZZt += M_Z1ZZtZ1;
	M_Z1ZZtZ1.rem();
//	sseMaM(Ra_ZZmZt, Ra_Z1ZZm1, Ra_ZZmZt, N_szLiu);
// 	sseUnmat(Ra_ZZp1_1tZpZt, N_szLiu);
// 	sseUnmat(Ra_Z1ZZm1, N_szLiu);

	M_ZZtZZt.setClean();
	return M_ZZtZZt.get();
}

xLiuCov* _getSKATo_rho(wsRealCst R_rho, wsUintCst N_szLiu,
	wsUintCst N_curSamp, wsUintCst N_gs, wsUintCst N_phe, cVector& V_wGScurr,
	cStdMatrix& M_genes, cSymMatrix& M_covs, cSymMatrix& M_varS,
	wsRealCst R_xx2sum, wsRealCst R_xxsum2, xGeneTestParam& X)
{
	xLiuCov *Xp_ret = NULL;
	pverbose("Loop for rho %g [%d?=%d]\n", R_rho, N_szLiu, N_curSamp);

	/*** MULTI-FARVAT ***/
	wsReal	R_mf	= WISARD_NAN;
	if (OPT_ENABLED(mfhom) || OPT_ENABLED(mfhet)) {
		/*** MULTI-FARVAT ***/
		if (OPT_ENABLED(mfhom)) {
			cSymMatrix	M_K(N_gs);
			wsVec		Ra_wGScurr	= V_wGScurr.get();
			wsSym		Ra_K		= M_K.get();

			LOOP(j, N_gs) {
				cVector V_Kj = M_K.r2v_ptr(j);
				sseVpC(Ra_wGScurr, Ra_wGScurr[j] * R_rho, V_Kj.get(), j);
//				V_wGScurr.mul(Ra_wGScurr[j] * R_rho, V_Kj);
				Ra_K[j][j] = Ra_wGScurr[j]*Ra_wGScurr[j];
			}
			// stats[i] <- t(IIq)%*%S%*%R%*%t(S)%*%IIq/sum(TphiAphiT) L263
			// S = TPAG
			cVector		V_1TPAG = X.V_1TPA.Mt(M_genes);
			/* 1TPA G R G' (1TPA)' */ 
			R_mf = V_1TPAG.qf(M_K);
			// pvals[i] <- pliu(stats[i],lius[[i]]) L264
			Xp_ret = liuCov(M_covs.get(), N_szLiu, Ra_K);

			M_varS.init(M_covs);
		} else if (OPT_ENABLED(mfhet)) {
			/* Expand K to [pq] */
			cSymMatrix M_K(R_rho, N_szLiu);
			M_K.setDiag(W1);

			/* Get varS = psi %*% TPAPT */
			cSymMatrix M_psi(M_covs.get(), N_gs, 1);
			M_varS = M_psi.krnk(X.M_TPAPT);

			/* Expand W */
			cDiagMatrix	M_W(N_gs, V_wGScurr.get(), 1);
			cIdtMatrix	M_q(N_phe);
			cDiagMatrix	M_Wq = M_q.krnk(M_W);

			/* Get R */
			cSymMatrix  M_R = M_Wq.MMt(M_K);

			/* Get S */
			cStdMatrix M_S = X.M_TPA.Mt(M_genes);

			Xp_ret = liuCov(M_varS.get(), N_szLiu, M_R.get());
			cVector V_S = M_S.vectorize();
			R_mf = V_S.qf(M_R);
		}
	} else if (!X.B_isCont) {
		/* Dichotomous trait */

		cSymMatrix	M_R(N_szLiu);
		wsVec		Ra_wGScurr	= V_wGScurr.get();
		wsSym		Ra_R		= M_R.get();

		LOOP(j, N_gs) {
			cVector	V_wc = V_wGScurr.clone(j+1);
			cVector V_Rj = M_R.r2v_ptr(j);
			V_wc.mul(Ra_wGScurr[j] * R_rho, V_Rj);
			Ra_R[j][j] += (W1-R_rho)*Ra_wGScurr[j]*Ra_wGScurr[j];
		}

		Xp_ret = liuCov(M_covs.get(), N_szLiu, Ra_R);
	} else {
		/* Continuous trait */

		/* Calculating K = (1-rho)*I + rho*11' */
#ifdef USE_CONT_SKATO_V2
		cStdMatrix	M_ttt1		= rho1(R_rho, X.M_r, X.M_1r);
		cStdMatrix	M_ttt2		= rho2(R_rho, X.M_r2, X.M_r1r, X.M_1r2);
		M_ttt1.setClean();
		M_ttt2.setClean();
		wsMat	Ra_ttt = M_ttt1.get();
		wsMat	Ra_ttt2 = M_ttt2.get();
#else
		/**/		wsSym	Ra_K = sseSymMat(N_gs);
		if (R_rho == W0) {
			for (j=0 ; j<N_gs ; j++) {
				sseVset0(Ra_K[j], j);
				Ra_K[j][j] = W1;
			}
		} else if (R_rho == W1) {
			for (j=0 ; j<N_gs ; j++)
				sseVinit(Ra_K[j], j+1, R_rho);
		} else for (j=0 ; j<N_gs ; j++) {
			sseVinit(Ra_K[j], j, R_rho);
			Ra_K[j][j] = W1;
		}
		cSymMatrix M_K(Ra_K, N_gs, 1);

		/* Get W_g' K */
		cStdMatrix M_UK = M_U * M_K;
		// /**/		MAT_t	Ra_UK	= sseMpS(M_U, N_curSamp*N_phe, N_gs,
		// 				Ra_K, N_gs);
		/* Get W_g P */
		cStdMatrix M_Ut = M_U.transpose();
		///**/		MAT_t	Ra_Ut	= transpose(M_U, N_curSamp*N_phe, N_gs);
		cSymMatrix M_P(Mp_P->get(), N_curSamp*N_phe, 1);
		cStdMatrix	M_UtP = M_U * M_P;
		///**/		MAT_t	Ra_UtP	= sseMpS(Ra_Ut, N_gs, N_curSamp*N_phe,
		//				Mp_P->get(), N_curSamp*N_phe);
		M_Ut.rem();
		//			sseUnmat(Ra_Ut, N_gs);
		/* Get ttt1 */
		cStdMatrix	M_ttt	= M_UK * M_UtP;
		// /**/		MAT_t	Ra_ttt		= sseMpM(Ra_UK, N_curSamp*N_phe, N_gs,
		// 				Ra_UtP, N_gs, N_curSamp*N_phe);
		/* Get ttt2 */
		wsMat Ra_UtP = M_UtP.get();
		/**/		wsSym	Ra_UtPU		= sym_sseMpM(Ra_UtP, N_gs, N_curSamp*N_phe,
			M_U.get(), N_curSamp*N_phe, N_gs);
		/**/		wsMat	Ra_KUtP		= sseSpM(Ra_K, N_gs,
			Ra_UtP, N_gs, N_curSamp*N_phe);
		sseUnmat(Ra_K, N_gs);
		M_UtP.rem();
		//			sseUnmat(Ra_UtP, N_gs);
		/**/		wsMat	Ra_ttt2		= sseMSMt(M_UK.get(), N_curSamp*N_phe, N_gs,
			Ra_UtPU, N_gs, N_gs, Ra_KUtP, N_gs, N_curSamp*N_phe);
		M_UK.rem();
		//			deallocMatrix(Ra_UK, N_curSamp*N_phe, (void *)1);
		sseUnmat(Ra_KUtP, N_gs);
		sseUnmat(Ra_UtPU, N_gs);
		//#endif
#ifdef GSvalidate
		char S_bufbuf[256];
		sprintf(S_bufbuf, "%s.rho%g.ttt", S_gs, R_rho);
		exportMatrix(S_bufbuf, Ra_ttt, N_szLiu, N_szLiu);
		sprintf(S_bufbuf, "%s.rho%g.ttt2", S_gs, R_rho);
		exportMatrix(S_bufbuf, Ra_ttt2, N_szLiu, N_szLiu);
#endif

		// 			LOG("ttt1 %g ttt2 %g\n",
		// 				compMMfrobenius(Ra_ttt, N_szLiu, N_szLiu, Ra_ttt_1,
		// 					N_szLiu, N_szLiu),
		// 				compMMfrobenius(Ra_ttt2, N_szLiu, N_szLiu, Ra_ttt2_1,
		// 					N_szLiu, N_szLiu));
#endif
		Xp_ret = liuCov(NULL, N_szLiu, NULL, Ra_ttt, Ra_ttt2);
		deallocMatrix(Ra_ttt, N_curSamp*N_phe, (void *)1);
		deallocMatrix(Ra_ttt2, N_curSamp*N_phe, (void *)1);
	}

	/* Derive statistics */
	if (OPT_ENABLED(mfhom) || OPT_ENABLED(mfhet))
		/* MULTI-FARVAT */
		Xp_ret->R_stat = R_mf / X.R_C;
	else
		/* Other case */
		Xp_ret->R_stat = (R_xx2sum*(W1-R_rho) + R_xxsum2*R_rho) / X.R_C;

	/* Derive p-value */
	// pvals[jj] <- pliu(stats[jj],lius[[jj]])	
	double R_nc = (double)(Xp_ret->R_delta);
	Xp_ret->R_pVal = (wsReal)PVchisq((Xp_ret->R_stat -
		Xp_ret->R_muQ) / Xp_ret->R_sigQ * Xp_ret->R_sigX
		+ Xp_ret->R_muX, Xp_ret->R_l, &R_nc);
	pverbose("Rho [%g], stat [%g], p-value [%g]\n", R_rho, Xp_ret->R_stat, Xp_ret->R_pVal);

	return Xp_ret;
}

xLiuCov** _getSKATo(wsUintCst N_phe,
	const char *S_gs, vInt &Xa_set,
	cStdMatrix&	M_U,		/* Checked = GW */
	cVector&	V_wGScurr,
//	VEC_t		Ra_wGScurr,	/* Checked */
	cStdMatrix&	M_genes,
	wsVec		Ra_xxG,			/* Checked */
	wsUintCst	N_curSamp,
	cVector&	V_ySKATo,
	wsReal*		Rp_rhoSKATO,
	wsReal*		Rp_Popt,
	cSymMatrix*	Mp_covs,
	xGeneTestParam& X)
{
	wsUint		N_gs		= (wsUint)Xa_set.size();
	cSymMatrix	M_covs; /* Checked */
	wsReal		R_xx2sum	= W0;
	wsReal		R_xxsum2	= W0;
	wsUint		N_szLiu;

	/* For continuous trait */
	if (OPT_ENABLED(mfhom)) {
		N_szLiu = N_gs;

		/* C is sum(TPAPT) */
		X.R_C = X.R_sumTPAPT;

		// psi <- cov(Genes) L541
		M_covs = M_genes.covR();
//		M_covs = sseMcovRow(M_genes, N_gs, N_curSamp, 1); /* Checked */
	} else if (OPT_ENABLED(mfhet)) {
		N_szLiu = N_gs*N_phe;
		/* C = 1 */
		X.R_C = W1;

		// psi <- cov(genes) L317
		M_covs = M_genes.covR();
//		M_covs = sseMcovRow(M_genes, N_gs, N_curSamp, 1); /* Checked */
	} else if (X.B_isCont) {
		cStdMatrix	M_Ut		= M_U.transpose();
		cVector		V_uuxx		= M_Ut * V_ySKATo;
		wsReal		R_uuxx2sum	= W0, R_uuxxsum2;

#ifdef GSvalidate
		M_U.fmtfile("skato.Wg.%s", S_gs);
#endif

		R_uuxxsum2	= V_uuxx.sum(&R_uuxx2sum);
		R_uuxxsum2	*= R_uuxxsum2;

		M_covs.init(X.Mp_P->get(), WISARD_NAN, N_curSamp*N_phe, 0, 1);
		R_xx2sum	= R_uuxx2sum;
		R_xxsum2	= R_uuxxsum2;

#ifdef USE_CONT_SKATO_V2				/* [#samp * #phe] * [#gs] */
		cVector		V_1Ut	= M_U.sumR();	/* [#samp * #phe] */
		cStdMatrix	M_UtP	= M_Ut * *((cSymMatrix *)(X.Mp_P));
			/* #gs * [#samp * #phe] */
		cVector		V_1UtP	= V_1Ut * *((cSymMatrix *)(X.Mp_P));
			/* [#samp * #phe] */
		M_Ut.rem();

		cVector		V_1UtPU	= V_1UtP * M_U; /* [#gs] */
		cSymMatrix	M_UtPU(sym_sseMpM(M_UtP.get(), N_gs, N_curSamp*N_phe,
			M_U.get(), N_curSamp*N_phe, N_gs), N_gs);

		cVector		V_1UtPUUtP	= V_1UtPU * M_UtP;
		X.M_r1r		= V_1Ut.tV(V_1UtPUUtP); // 1u' 1u P U U' P

		cVector		V_UUtP1U	= V_1UtPU.Mt(M_U); // U U' P 1u
		cStdMatrix	M_1rho_2	= V_UUtP1U.tV(V_1UtP); // U U' P 1u 1u' P

		X.M_r1r += M_1rho_2;
		/*         [n1] [1n][nn][np] [pn][nn] */
		/* (1-p) [ 1u (1u)'PU U'P + U U'P1u (1u)'P */

		cStdMatrix	M_UUtPU		= M_U * M_UtPU; // U U' P U
		X.M_1r2		= M_UUtPU * M_UtP; // U U' P U U' P
		/* (1-p)^2 [ U U'PU U'P ] */

		X.M_r.rem();
		X.M_r		= V_1Ut.tV(V_1UtP);	// 1u 1u' P
		X.M_r2		= X.M_r * V_1UtP.sum(V_1Ut); // 
		/* p^2 [ 1u (1u)'P1u (1u)'P ] */

		X.M_1r.rem();
		X.M_1r		= M_U * M_UtP;
		//M_1r	= V_1Ut.tV(V_1UtP);
#endif

		/* C should be 1 */
		X.R_C = W1;

		N_szLiu = N_curSamp*N_phe;
	} else {
		/* For dichotomous trait */

#ifdef GSvalidate
		V_wGScurr.fmtfile("skato.weighZt.%s", S_gs);
#endif

		// xx <- xxG%*%W
		// xx2sum   <- sum(xx^2)
		// xxsum2   <- sum(xx)^2
		wsVec Ra_wGScurr = V_wGScurr.get();
		LOOP(i, N_gs) {
			Ra_xxG[i]	*= Ra_wGScurr[i];
			R_xx2sum	+= Ra_xxG[i]*Ra_xxG[i];
			R_xxsum2	+= Ra_xxG[i];
		}
		//sseFree(Ra_xxG);
		R_xxsum2 *= R_xxsum2;

#ifdef GSvalidate
		M_genes.fmtfile("skato.genes.%s", S_gs);
#endif

		// covs     <- cov(genes)
		if (Mp_covs) {
			M_covs = *Mp_covs;
			M_covs.setDontDealloc();
		} else
			M_covs = M_genes.covR();
#ifdef GSvalidate
		M_covs.fmtfile("skato.cov.%s", S_gs);
#endif

		N_szLiu = N_gs;
	}

	cTimer		c2;
	wsReal		R_minPv	= W2;
	xLiuCov**	Xa_liu	= NULL; /* Checked */
	cSymMatrix	M_varS;	/* For MULTI-FARVAT */
	wsAlloc(Xa_liu, xLiuCov*, X.N_rho);
	//wsUint N_idxMin = 0;

	/* For each rho, compute statistics and p-value */
	LOOP (i, X.N_rho) {
		/* Get LIU result */
		Xa_liu[i] = _getSKATo_rho(X.Ra_rho[i], N_szLiu, N_curSamp, N_gs,
			N_phe, V_wGScurr, M_genes, M_covs, M_varS, R_xx2sum,
			R_xxsum2, X);

		/* Minimum p-value check */
		if (Xa_liu[i]->R_pVal < R_minPv) {
			R_minPv = Xa_liu[i]->R_pVal;
			*Rp_rhoSKATO = X.Ra_rho[i];
		}
	}
// 	deallocMatrix(Ra_Ir1, N_actualSamp*N_pheno, (void *)1);
// 	deallocMatrix(Ra_Ir2, N_actualSamp*N_pheno, (void *)1);

#ifdef USE_CONT_SKATO_V2
	X.M_r2.rem();
#endif

	wsMat	Ra_IMV	= NULL;
	wsMat	Ra_IMV2	= NULL;
	wsReal*	Ra_rRho	= sseVector(X.N_rho);
	wsReal	R_muQ	= WISARD_NAN;
	wsReal	R_sigQ	= WISARD_NAN;
	wsReal	R_sigXi	= WISARD_NAN;
	if (OPT_ENABLED(mfhom) || OPT_ENABLED(mfhet)) {
		/*
		e <- eigen(psi)
		evectors <- Re(e$vectors)
		evalues <-Re(e$values)
		o <- evalues<0
		evalues[o] <- 0
		Z <- evectors%*%diag(sqrt(evalues))%*%t(evectors)%*%W
		*/
		//cSymMatrix	M_psi(Ra_covs, N_szLiu, 1);
		cDiagMatrix	M_W;
		cStdMatrix	M_Z;
		if (OPT_ENABLED(mfhom)) {
			cSymMatrix M_psi(M_covs.get(), N_gs, 1);
			cSymMatrix&	M_psi05 = M_psi.sqrt();
			M_W.init(N_gs, V_wGScurr.get(), 0, 1);
			M_Z = M_psi05 * M_W;
			delete &M_psi05;
		} else {
			cSymMatrix&	M_psi05 = M_varS.sqrt();
			wsReal *Ra_wGScurrRep = sseVrep(V_wGScurr.get(), N_gs, N_phe);
			M_W.init(N_gs*N_phe, Ra_wGScurrRep);
			M_Z = M_psi05 * M_W;
			delete &M_psi05;
		}

		/* zbar <- matrix(rowMeans(Z), ncol=1) */
		cVector		V_zBar	= M_Z.meanR();
		/* zbar2 <- t(zbar)%*%zbar */
		wsReal		R_zBar2 = V_zBar.ss();
		/* M <- zbar%*%solve(t(zbar)%*%zbar)%*%t(zbar) */
		cSymMatrix	M_M		= V_zBar.tV(W1/R_zBar2);

		/* IMV <-  (Im-M)%*%Z%*%t(Z)%*%(Im-M) */
		cIdtMatrix M_1(N_szLiu);
		cSymMatrix M_1M		= M_M.subI();
		cStdMatrix M_1MZ	= M_1M * M_Z;

		cSymMatrix M_IMV	= M_1MZ.Mt();
		cSymMatrix M_IMV2	= M_IMV.Mt();
		Ra_IMV	= M_IMV.get();
		Ra_IMV2	= M_IMV2.get();
		M_IMV.setDontDealloc();
		M_IMV2.setDontDealloc();
		/* muQ <- sum(diag(IMV)) */
		/*
		 *	IMV2 <- IMV%*%IMV
		 *	muQ2 <- sum(diag(IMV2)) == sum(IMV2)
		 */
		R_muQ = M_IMV.tr();
		wsReal R_muQ2 = M_IMV2.tr();

		/* sigXi <- 2*sqrt( sum( diag(t(Z)%*%M%*%Z%*%t(Z)%*%(Im-M)%*%Z) ) ) */
		cStdMatrix	M_ZtM	= M_Z.tM(M_M);
		cSymMatrix	M_ZZt	= M_Z.Mt();
		cStdMatrix	M_ZZIMZ	= M_ZZt * M_1MZ;
		R_sigXi = diagSum2matrix(M_ZtM.get(), N_szLiu, N_szLiu,
			M_ZZIMZ.get(), N_szLiu, N_szLiu);
		R_sigXi = W2 * sqrt(R_sigXi);

		/* sigQ <- sqrt( 2*muQ2 + sigXi^2 ) */
		R_sigQ = sqrt(W2*R_muQ2 + SQR(R_sigXi));

		/* rRho <- nSNP^2*rhos*zbar2 + (1-rhos)/zbar2*sum((t(zbar)%*%Z)^2) */
		cVector		V_zBz	= V_zBar * M_Z;
		wsReal		R_z2zz	= V_zBz.ss() / R_zBar2;
		wsReal		R_xdim	= (wsReal)N_szLiu * (wsReal)N_szLiu;

		/* Set Ra_rRho */
		LOOP (i, X.N_rho) Ra_rRho[i] = R_xdim*X.Ra_rho[i]*R_zBar2
				+ (W1-X.Ra_rho[i]) * R_z2zz;
	} else {
		/*
		 * Normal univariate FARVAT
		 */

		wsMat	Ra_Z		= NULL;	/* Checked */
		wsReal	*Ra_zBar	= NULL;	/* Checked */
		wsReal	R_z2		= W0;
		c2.start();

		if (X.B_isCont) {
			/*
			 * Continuous case
			 *
			 */

			// Z     <- P0.5%*%t(W)	# [n*n] %*% [n*p] = [n*p]
			// zbar  <- matrix(rowMeans(Z),nrow=1)
			Ra_Z	= sseSpM(X.Mp_Psq->get(), N_curSamp*N_phe,
				M_U.get(), N_curSamp*N_phe, N_gs);
// 			Ra_zBar	= sseMmeanR(Ra_Z, N_actualSamp*N_pheno, N_gs);
// 			// z2    <- sum(zbar^2)
// 			for (i=0 ; i<N_actualSamp*N_pheno ; i++)
// 				R_z2 += SQR(Ra_zBar[i]);
		} else {
			/*
			 * Binary case, comparison w/ FARVAT R code
			 */

			/*
e <- eigen(covs)
evectors <- e$vectors
evalues  <- e$values
o <- evalues <  0
evalues[o]   <- 0
			*/
			char		B_isFail = 0;
			cStdMatrix	M_eVec;
			cVector		V_eVal	= M_covs.eigen(M_eVec);
			if (B_isFail) {
				char fn[256];
				sprintf(fn, "skato.error.covs.%s", S_gs);
				M_covs.file(fn);
				LOGwarn("[%s] : Eigendecomposition of covar. mat. failed",
					S_gs);
			}

			/* Take square root to get ppp^(1/2) */
			V_eVal.selfSqrt(1);

			/* Calculate Z matrix */
			Ra_Z		= sseEmptyMatrix(N_gs, N_gs);
			sseCalloc(Ra_zBar, wsReal, N_gs);

			// Z <- evectors %*% diag(sqrt(evalues)) %*% t(evectors) %*% W
			// zbar <- matrix(rowMeans(Z),nrow=1)
			// z2   <- sum(zbar^2)
			wsMat Ra_eVec = M_eVec.get();
			LOOP(i, N_gs) {
				cVector	V_curRow	= M_eVec.r2v_ptr(i);
				wsVec	Ra_wGScurr	= V_wGScurr.get();
				cVector	V_b			= V_curRow * V_eVal;
				wsVec	Ra_b		= V_b.get();
// 				for (wsUint k=0 ; k<N_gs ; k++)
// 					Ra_b[k] = Ra_eVec[i][k]*Ra_eVal[k];
				LOOP(j, N_gs) {
					for (wsUint k=0 ; k<N_gs ; k++)
						Ra_Z[i][j] += Ra_b[k]*Ra_eVec[j][k];
					Ra_Z[i][j] *= Ra_wGScurr[j];
					Ra_zBar[i] += Ra_Z[i][j];
				}
				Ra_zBar[i]	/= (wsReal)N_gs;
				R_z2		+= SQR(Ra_zBar[i]);
			}
//			sseUnmat(Ra_eVec, N_gs);
		}
		//sseFree(Ra_wGScurr);

		/* Set Z matrix */
		cStdMatrix M_Z(N_szLiu, N_gs, Ra_Z);

		// Z2    <- Z%*%t(Z)
		cSymMatrix M_ZZt = M_Z.Mt();
#ifdef GSvalidate
		M_ZZt.fmtfile("skato.Z2.%s", S_gs);
#endif

		wsReal R_sumZbarZ = W0;
		// ZbZ   <- zbar%*%Z
		if (!X.B_isCont) {

			/*
			 * BINARY SKAT-O
			 */

			// zbar%*%Z
			cVector V_zBar(Ra_zBar, N_szLiu);
			cVector V_ZbarZ = V_zBar * M_Z;

			// sum( (zbar%*%Z)^2 )
			R_sumZbarZ = V_ZbarZ.ss();
			V_ZbarZ.rem();
			pverbose("[%s] Prepare Z\n", c2.getReadable());

		// M    <- t(zbar)%*%zbar / z2
		// IM   <- I-M
		// IMV = (I-M)ZZ'

		/* Define u(j) = \Sigma_(k=1)^n( z\bar_k (ZZ')_jk ) */
//	/**/wsReal	**Ra_u = NULL;
	//	for (i=0 ; i<N_szLiu ; i++)
//			Ra_u = sseMS(&Ra_zBar, 1, N_szLiu, Ra_Z2, N_szLiu);

		/* Now [(I-M)ZZ']_ij = (ZZ')_ij - Z\bar_i * u(j) */
			c2.start();
 			wsSym Ra_IM = sseSymMat(N_szLiu);
			LOOP(i, N_szLiu) {
				wsReal R_mul = -Ra_zBar[i] / R_z2;
				sseVpC(Ra_zBar, R_mul, Ra_IM[i], i);

				wsReal R_ii = -Ra_zBar[i]*R_mul;
				Ra_IM[i][i] = W1 - R_ii;
			}
			V_zBar.rem();
 		//sseFree(Ra_zBar);
			cSymMatrix M_IM(Ra_IM, N_szLiu);
			cStdMatrix M_IMZ;
			/* IMV = (I-M) ZZ' */
			cStdMatrix M_IMV = M_IM * M_ZZt;
			M_IMV.setClean();
			Ra_IMV = M_IMV.get();
	// 
	// 	// IMV  <- IM %*% Z2
	// 	//wsReal **Ra_IMV = sseThrRowMMt(Ra_IM, N_szLiu, N_szLiu, Ra_Z2,
	// 	//	N_szLiu, N_szLiu); /* Checked */
	// 	wsReal **Ra_IMV = sseThrRowSS(Ra_IM, N_szLiu, Ra_Z2, N_szLiu);
	// 	
			cStdMatrix M_Zt;
			cSymMatrix M_ZtIMZ;
			if (X.B_isCont) {
				/* GetIMV2 = non-n^3 version */
				cStdMatrix M_Z(N_szLiu, N_gs, Ra_Z, MATDEL_NONE);
				Ra_IMV2 = _getIMV2_p3(M_Z, R_z2);
			} else {
			// IMV2 <- IMV%*%IMV
				/* GetIMV2 = non-p^3 version */
				//cSymMatrix M_IMV2 = M_IMV.Mt();
				cStdMatrix M_IMV(N_szLiu, N_szLiu, Ra_IMV, MATDEL_NONE);
				cStdMatrix M_IMV2 = M_IMV * M_IMV;
				M_IMV2.setClean();
				Ra_IMV2 = M_IMV2.get();
			}
		//sseUnmat(Ra_Z, N_szLiu);
		} else {
			/* X <- t(Z) %*% Z */
			cSymMatrix	M_X		= M_Z.tM();
			/* iXi <- 1 / sum(X) */
			wsReal		R_Xsum	= M_X.sum();
			wsReal		R_invXs	= W1/R_Xsum;

			/* ZX <- Z %*% X */
			cStdMatrix	M_ZX	= M_Z * M_X;
			/* Z1 <- matrix(rowSums(Z), ncol=1) */
			cVector		V_Z1	= M_Z.sumR();
			/* ZX1 <- matrix(rowSums(ZX), ncol=1) */
			cVector		V_1XZ	= M_ZX.sumR();
			/* IMVs <- Z1 %*% t(ZX1) * iXi */
			cStdMatrix	M_IMVs	= V_Z1.tV(V_1XZ, R_invXs);
			/* IMV.my <- Z%*%t(Z) - IMVs */
			cStdMatrix	M_IMV	= M_ZZt - M_IMVs;
			M_IMV.setClean();
			Ra_IMV = M_IMV.get();

			/* ZXZt <- ZX %*% t(Z) */
			wsMat Ra_ZX = M_ZX.get();
			//SYM_t Ra_X = M_X.get();
			wsSym Ra_ZXZt = sym_sseMpMt(Ra_ZX, N_szLiu, N_gs, Ra_Z, N_szLiu, N_gs);
			cSymMatrix M_ZXZt(Ra_ZXZt, N_szLiu);
			/* IMVs1 <- ZX1 %*% t(ZX1) / iXi */
			cSymMatrix M_IMVs1 = V_1XZ.tV(R_invXs);
			/* ZXZt <- ZXZt - IMVs1 */
			M_ZXZt -= M_IMVs1;

			cVector V_1X = M_X.sumR(&R_z2);
			R_z2 /= ((wsReal)N_gs*(wsReal)N_gs);
			/* iXXZ2 <- t(X1) %*% t(ZX) */
			cVector V_1XXZt = V_1X.Mt(M_ZX);
			/* IMVs2 <- Z1 %*% iXXZ / iXi */
			cStdMatrix M_IMVs2 = V_Z1.tV(V_1XXZt, R_invXs);
			/* iXXi <- sum(X %*% X) */
			wsReal R_1XX1 = V_1X.ss();
			R_sumZbarZ = R_1XX1 / ((wsReal)N_gs*(wsReal)N_gs);
			/* IMVs3 <- Z1 %*% t(ZX1) * iXXi/ iXi^2 */
			cStdMatrix M_IMVs3 = V_Z1.tV(V_1XZ, R_1XX1*SQR(R_invXs));
			/* IMVs2 <- IMVs2 - IMVs3 */
			M_IMVs2 -= M_IMVs3;
			/* IMV2.my <- ZXZt - IMVs2 */
			cStdMatrix M_IMV2 = M_ZXZt - M_IMVs2;
			M_IMV2.setClean();
			Ra_IMV2 = M_IMV2.get();
		}

		// muQ  <- sum(diag(IMV))
		R_muQ = W0;
		wsReal R_sumIMV2 = W0;
		LOOP (i, N_szLiu) {
			R_muQ += Ra_IMV[i][i];
			R_sumIMV2 += Ra_IMV2[i][i];
		}
		pverbose("[%s] Prepare IMV\n", c2.getReadable());

//		if (!B_isCont)
//			R_muQ = R_sumIMV2;

		// rRho <- nSNP^2*rhos*z2 + (1-rhos)/z2 * ssZbZ
		LOOP (i, X.N_rho)
			Ra_rRho[i] = SQR(N_gs)*X.Ra_rho[i]*R_z2
				+ (W1-X.Ra_rho[i])/R_z2 * R_sumZbarZ;

		// t(Z)%*%M%*%Z%*%t(Z)%*%(Imq-M)%*%Z)
		wsReal R_dsum = WISARD_NAN;
		// Z2IMV <- Z2-IMV
		c2.start();
//			cSymMatrix M_PsIMV = M_ZZt - M_IMV;
 		wsMat Ra_PsIMV = sseMatrix(N_szLiu, N_szLiu); /* Checked */
 		sseSsM(M_ZZt.get(), Ra_IMV, Ra_PsIMV, N_szLiu);

		// sigXi <- 2 * sqrt( sum(diag( Z2IMV%*%IMV )) )
		//wsReal R_dsum = M_PsIMV.tr(M_IMV);
		R_dsum = diagSum2matrix(Ra_PsIMV, N_szLiu, N_szLiu,
			Ra_IMV, N_szLiu, N_szLiu);
		sseUnmat(Ra_PsIMV, N_szLiu);
		if (R_dsum < W0) /* Possible to down to 0 from round error */
			R_dsum = W0;
		R_sigXi = W2 * sqrt(R_dsum);
		//sseUnmat(Ra_PsIMV, N_szLiu);
		pverbose("[%s] Prepare sigXi\n", c2.getReadable());

		// sigQ <- sqrt( 2*sum(diag(IMV2)) + sigXi^2 )
		R_sigQ = sqrt(W2*R_sumIMV2 + SQR(R_sigXi));
	}
	if (!(X.Mp_P) || M_covs.get() != X.Mp_P->get())
		M_covs.rem();

	// id   <- order(pvals)[1]
	// rhoMin <- rhos[id]
	// pMin   <- pvals[id]
	// for(jj in seq(along=rhos)) {
	c2.start();
	wsReal *Ra_qLiu = NULL; /* Checked */
	wsAlloc(Ra_qLiu, wsReal, X.N_rho);
	LOOP (i, X.N_rho) {
		// rho <- rhos[jj]
		// qmins[jj] <- qliu(pMin,lius[[jj]])
		Ra_qLiu[i] = qLiu(R_minPv, Xa_liu[i]);
	}

	// lius      <- getliu(covs,I - M)
	//[121210] lius      <- getliu(Z2,I - M)
	xLiuCov *Xp_liuFinal = NULL;
	if (OPT_ENABLED(mfhet) || OPT_ENABLED(mfhom))
		Xp_liuFinal = liuCovSS(Ra_IMV, Ra_IMV2, N_szLiu); /* Checked */
	else
		Xp_liuFinal = liuCov(NULL, N_szLiu, NULL, Ra_IMV, Ra_IMV2);
	sseUnmat(Ra_IMV, N_szLiu);
	sseUnmat(Ra_IMV2, N_szLiu);
	pverbose("[%s] Get liuFinal\n", c2.getReadable());

	// 	pvalSKATO <- 1 - unlist(integrate( integrateSKATo, 0, Inf, qmins=qmins, rRho=rRho,
	// 	rhos=rhos, muQ=muQ, sigQ=sigQ, sigXi=sigXi,lius=lius))$value
	xIntgSKATO X_par;
	X_par.R_muQ			= R_muQ;
	X_par.R_sigQ		= R_sigQ;
	X_par.R_sigXi		= R_sigXi;
	X_par.Ra_qLiu		= Ra_qLiu;
	X_par.Ra_rho		= X.Ra_rho;
	X_par.Ra_rRho		= Ra_rRho;
	X_par.Xp_liuFinal	= Xp_liuFinal;
	X_par.N_rho			= X.N_rho;
	xIntgInp X_inp(integrateOptimal, 0, 40.0, &X_par, 500, pow(numeric_limits<double>::epsilon(), 0.25)); // 140507 update
	xIntgRes *Xp_res = integrate(&X_inp); /* Checked */
//	LOG("%g %d\n", Xp_res->R_value, Xp_res->N_err);
	if (NA(Xp_res->R_value) || Xp_res->R_value == 0 || Xp_res->N_err || Xp_res->R_value > 0.99999) {
		xIntgInp X_inp2(integrateOptimal, 0, 40.0, &X_par, 2000, 1e-25); // 140507 update
//		DEALLOC(Xp_res);
//		LOG("Do compensate...\n");
		xIntgRes *Xp_res2 = integrate(&X_inp2); /* Checked */
		if (Xp_res2->R_value == 0 || NA(Xp_res2->R_value)) {
			DEALLOC(Xp_res2);
			xIntgInp X_inp3(integrateOptimal, 0, numeric_limits<double>::infinity(), &X_par, 2000, 1e-25); // 140507 update
			//		DEALLOC(Xp_res);
			//		LOG("Do compensate...\n");
			xIntgRes *Xp_res3 = integrate(&X_inp3); /* Checked */
			if (Xp_res3->R_value == 0 || NA(Xp_res3->R_value)) {
				DEALLOC(Xp_res3);
			} else {
				DEALLOC(Xp_res);
				Xp_res = Xp_res3;
			}
		} else {
			DEALLOC(Xp_res);
			Xp_res = Xp_res2;
		}
	}
	sseFree(Ra_rRho);
	DEALLOC(Ra_qLiu);
	DEALLOC(Xp_liuFinal);

//	if (Xp_res->N_err)
//		*Rp_pValSKATO = WISARD_NA_REAL;
//	else
		*Rp_Popt = (wsReal)(1.0 - Xp_res->R_value);
	DEALLOC(Xp_res);

	return Xa_liu;
}

wsUint _getSKAT(wsRealCst R_pppMult, const char *S_gs, vInt &Xa_set,
	cStdMatrix& M_adjGeno, wsRealCst R_stat, wsUintCst N_actualSamp, wsReal *Rp_pVal)
{
#ifdef GSvalidate
	char fn[256];
	sprintf(fn, "%s.X", S_gs);
#endif

	/* Review of original code
	 * 
		II   <- matrix(1, nrow=totind, ncol=1)	# [#i,1]
		phiInvII       <- akinshipInv%*%II		# [#i,#i]%*%[#i,1]
		IIphiInvII     <- sum(phiInvII)			# [1]
		IIphiInvII_Inv <- 1/IIphiInvII			# [1]
		hat    <- t(phiInvII)*IIphiInvII_Inv	# [1,#i]%*%[1]
		hatMat <- II%*%hat						# [#i,#i]
		# Each cell of g2 (ith ind, jth snp) is come from (1-2)
		#  (1) Original value of data[i][j]
		#  (2) hatMat[i,]%*%genes[,j]
		g2 <- genes - hatMat%*%genes
		# [#i,#s]
		# set YY to equal across all SNPs for same individual
		YY     <- adjpheno(phens, prev)
		pheMat <- matrix(rep(YY,nSNP),nrow=nInd,ncol=nSNP)
		# Each cell of xx before colSums() (ith ind, jth snp)
		#  is YY[i]*g2[i,j]
		# xx have 1*#s dim. so all info of individual for same SNP
		#  should be concatenated
		# xx[j] = sum(YY[*]*g2[*,j])
		xx    <- colSums(pheMat*g2)
		# Final stats is SCALAR, which is the squared sum of xx
		stats <- sum(xx^2)
	 */

	//exportMatrix("Y.genetest", &Ra_Y, 1, N_availPhenoSamp);
	wsUint	N_gs	= (wsUint)Xa_set.size();		///< Size of current gene-set

#ifdef GSvalidate
	char fs[256];
	sprintf(fs, "gs.SKATg2.%s", S_gs);
	M_adjGeno.file(fs);
#endif

	// 	ppp   <- pppMult*cov(g2)
	
	/* Get the variance-covariance matrix of Ra_covGeno */

	// [#gs*#gs] pppCov <- cov(g2)
//	wsUint		N_pheno	= getIO()->getPhenoCount();
	cSymMatrix	M_ppp	= M_adjGeno.covR();

	// ppp <- pppMult*pppCov
	// 121026 After [ppp <- pppCov]
	//sseMmS(Ra_ppp, R_pppMult, Ra_ppp, N_gs);
	M_adjGeno.rem();
//	sseUnmat(M_adjGeno, N_gs);

#ifdef GSvalidate
	M_ppp.fmtfile("gs.SKATppp.%s", S_gs);
#endif

	//		eig   <- eigen(ppp)$values
	cVector V_eigValue; /* 120917 */
#ifndef USE_DBL
	double *Ra_eigValueD	= NULL; /* 120917 */
	wsAlloc(Ra_eigValueD, double, N_gs);
#endif

// 	sprintf(fn, "%s.ppp", S_gs);
// 	exportMatrix(fn, Ra_ppp, N_gs, N_gs);
 	wsUint N_ev;

	/* Perform triangular decomposition */
#ifdef USE_JCBI
	jcbi(M_ppp, V_eigValue, N_gs, 0.000001, 100000);
	N_ev = N_gs;
#elif USE_POWER
	if (N_gs > 800) {
		N_ev = 50;
		V_eigValue = EVpower(M_ppp, N_gs, N_ev, 1e-10);

		/* Get trace of ppp */
		wsReal R_trace = W0;
		for (wsUint i=0 ; i<N_gs ; i++)
			R_trace += M_ppp[i][i];
		wsReal R_sumEV = W0;
		for (wsUint i=0 ; i<N_ev ; i++)
			R_sumEV += V_eigValue[i];

		wsReal *Ra_newEV = NULL;
		wsAlloc(Ra_newEV, wsReal, ++N_ev);
		memcpy(Ra_newEV, V_eigValue, sizeof(wsReal)*(N_ev-1));
		Ra_newEV[N_ev-1] = R_trace-R_sumEV;
		LOG("Set last EV to %g\n", R_trace-R_sumEV);
		DEALLOC(V_eigValue);
		V_eigValue = Ra_newEV;
	} else {
		char B_ret;
		V_eigValue = eigen(M_ppp, N_gs, &B_ret);
		if (B_ret) {
			LOG("Gene set [%s] have nonnegative definitive eigenvalues (mult %g)\n",
				S_gs, R_pppMult);

			char fn[256];
			sprintf(fn, "%s.ppp", S_gs);
			exportMatrix(fn, M_ppp, N_gs, N_gs);
		}
		N_ev = N_gs;
	}
#elif USE_LANCZOS
	if (N_gs > 800) {
		V_eigValue = EVlanczos(M_ppp, N_gs, N_gs, &N_ev);
		LOG("Lanczos to [%s] : %d EVs, last EV %g\n", S_gs, N_ev,
			V_eigValue[N_ev-1]);

		/* Get trace of ppp */
		wsReal R_trace = W0;
		for (wsUint i=0 ; i<N_gs ; i++)
			R_trace += M_ppp[i][i];
		wsReal R_sumEV = W0;
		for (wsUint i=0 ; i<N_ev ; i++)
			R_sumEV += V_eigValue[i];

		wsReal *Ra_newEV = NULL;
		wsAlloc(Ra_newEV, wsReal, ++N_ev);
		memcpy(Ra_newEV, V_eigValue, sizeof(wsReal)*(N_ev-1));
		Ra_newEV[N_ev-1] = R_trace-R_sumEV;
		LOG("Set last EV to %g\n", R_trace-R_sumEV);
		DEALLOC(V_eigValue);
		V_eigValue = Ra_newEV;
	} else {
		char B_ret;
		V_eigValue = eigen(M_ppp, N_gs, &B_ret);
		if (B_ret) {
			LOG("Gene set [%s] have nonnegative definitive eigenvalues (mult %g)\n",
				S_gs, R_pppMult);

			char fn[256];
			sprintf(fn, "%s.ppp", S_gs);
			exportMatrix(fn, M_ppp, N_gs, N_gs);
		}
		N_ev = N_gs;
	}
#else
	char B_ret = 0;
#ifdef USE_ED2
	V_eigValue = M_ppp.eigen();//EIGENDECOMPOSITION(M_ppp, N_gs);
#else
	V_eigValue = eigen(M_ppp, N_gs, &B_ret);
#endif
	if (B_ret) {
		LOG("Gene [%s] have nonnegative definitive eigenvalues (mult %g)\n",
			S_gs, R_pppMult);

		char fn[256];
		sprintf(fn, "%s.ppp", S_gs);
		M_ppp.file(fn);
	}
	N_ev = N_gs;
#endif

	if (OPT_ENABLED(verbose)) {
		char fs[256];
		sprintf(fs, "gs.EValue.%s", S_gs);
		V_eigValue.file(fs);
	}

#ifndef USE_DBL
	for (wsUint i=0 ; i<N_ev ; i++)
		Ra_eigValueD[i] = V_eigValue[i];
#endif
	M_ppp.rem();
//	sseUnmat(M_ppp, N_gs);

	// 121016 After [stat divided into R_pppMult]
	//		pval  <- davies(231217.64734915376,eig)$Qq
#ifndef USE_DBL
	wsReal R_pVal = (wsReal)davies(R_stat/R_pppMult, Ra_eigValueD, N_ev);
	DEALLOC(Ra_eigValueD);
#else
	wsReal R_pVal = (wsReal)davies(R_stat/R_pppMult, V_eigValue.get(), N_ev);
#endif
	wsUint B_isLiu = 0;
	if (R_pVal <= W0) {
#ifndef USE_DBL
		R_pVal = (wsReal)liuEV(R_stat, Ra_eigValueD, N_ev);
#else
		R_pVal = liuEV(R_stat, V_eigValue.get(), N_ev);
#endif
		B_isLiu = 1;
	}
	V_eigValue.rem();
	//sseFree(V_eigValue);

	//	c(stat=stats, pval=pval)
	*Rp_pVal = R_pVal;
	return B_isLiu;
}

wsRealCst _getCMC(cMask &V_mask, cVector& V_hat, wsRealCst R_vars,
	cVector& V_aV,
	cVector &V_pldCMC, wsReal *Rp_statCMC, wsUintCst N_samp)
{
	// 	ppp   <- (hat%*%pooledGene)[1,1]
	wsReal	R_ppp			= W0;
	wsReal	*Ra_maskTest	= V_mask.get();
	wsReal	*Ra_pldCMC		= V_pldCMC.get();
//	wsReal	*Ra_hat			= V_hat.get();

	/* [N_testSamp]	Ra_hat
	 * [N_testSamp]	Ra_pldCMC */
	R_ppp = V_mask.vv(V_hat, V_pldCMC);
// 	for (wsUint i=0 ; i<N_samp ; i++) {
// 		if (Ra_maskTest[i])
// 			R_ppp += Ra_hat[i]*Ra_pldCMC[i];
// 	}

	if (R_ppp < 1.0) {
		// 	XX    <- pooledGene/sqrt(ppp*(1 - ppp))
		wsReal R_sqrtp1p = sqrt(R_ppp*(W1-R_ppp));
		for (wsUint i=0 ; i<N_samp ; i++)
			if (Ra_maskTest[i])
				Ra_pldCMC[i] /= R_sqrtp1p;

		return _qlsUniv(V_mask, V_aV, V_pldCMC, R_vars, Rp_statCMC);
	} else {
		*Rp_statCMC	= numeric_limits<wsReal>::infinity();

		return W1;
	}
}

wsRealCst _getCollapsing(cMask &V_mask, wsRealCst R_vars, cVector& V_aV,
	cVector &V_pldClp,
	wsReal *Rp_statClp, wsUint N_samp, wsReal *Rp_aVy)
{
	if (Rp_statClp == NULL)
		halt_fmt(WISARD_SYST_NULL_CLPGENO);
	wsReal R_mean = W0;
	wsReal R_var1 = W0, R_var2 = W0;
	wsReal *Ra_pldClp = V_pldClp.get();

	// 	ttt   <- var(pooledGene)[1,1]
	//_getGeneWiseStat(it, GSTEST_COLLAPSING, Ra_w, Ra_pooled, Ba_misPheno);
	//exportMatrix("xg.genetest", Ra_ret, N_availPhenoSamp, 1);
	/* More precise estimation of mean/variance */
	wsUint	N = 0;
	wsReal	*Ra_maskTest = V_mask.get();
	for (wsUint i=0 ; i<N_samp ; i++)
		if (Ra_maskTest[i]) {
			R_mean += Ra_pldClp[i];
			N++;
		}
	R_mean /= (wsReal)N;
	for (wsUint i=0 ; i<N_samp ; i++) {
		if (Ra_maskTest[i]) {
			R_var1 += SQR(Ra_pldClp[i]);
			R_var2 -= Ra_pldClp[i];
		}
	}
	R_var1 += W2*R_mean*R_var2 + N*SQR(R_mean);
	R_var1 /= (wsReal)(N-1);
	R_var1 = sqrt(R_var1);
	if (R_var1 != R_var1)
		return WISARD_NA_REAL;

	// 	XX    <- pooledGene/sqrt(ttt)
	for (wsUint i=0 ; i<N_samp ; i++)
		if (Ra_maskTest[i])
			Ra_pldClp[i] /= R_var1;

	return _qlsUniv(V_mask, V_aV, V_pldClp, R_vars, Rp_statClp, Rp_aVy);
}

wsRealCst _getPEDCMC(cIO* Cp_IO, wsUintCst N_samp, cMatrix* Mp_phi, cVector* Vp_yGSunAdj,
	const char *S_gs, vInt &Xa_set,
	wsReal **Ra_PEDCMCca, wsReal **Ra_PEDCMCct, wsReal *Rp_stat,
	wsUint N_ctrl, wsUint N_case, wsUint N_common, wsUint N_rare,
	cMask &V_mask, char *Ba_filt, wsUint *Np_typeTest, wsUint *Np_rank)
{
//	vVariant&	Xa_vrt		= Cp_IO->getVariants();
	wsReal*		Ra_Yp		= NULL;
	wsReal		R_pV;
	wsUint		i, Ica, Ico, I;

	/* Process special case */
	if (N_ctrl <= 1 || N_case <= 1) {
		if (Rp_stat) *Rp_stat = numeric_limits<wsReal>::infinity();
		return W1;
	}

// 	/* Count the number of included individual */
// 	for (i=I=0 ; i<N_samp ; i++) {
// 		if (Ba_misPheno[i]) continue;
// 
// 		if (Ra_Y[I] == 0)		N_ctrl++;
// 		else if (Ra_Y[I] == 1)	N_case++;
// 		else halt("Cannot perform --pedcmc due to this data is not dichotomous phenotype");
// 		I++;
// 	}
// 	if ((N_case+N_ctrl) != N_testSamp)
// 		halt("SYSERR: Counted sample size [case %d+control %d=%d] does not match"
// 			"with expected number of samples [%d]", N_case, N_ctrl, N_testSamp);
// 
// 	FOREACH (vector<int>::iterator, gs->second, it) {
// 		wsUint	N_idx	= *it;
// 		xSNP	&X_SNP	= Xa_SNP[N_idx];
// 
// 		if (X_SNP.maf < OPT_REAL(raremaf))
// 			N_rare++;
// 		else
// 			N_common++;
// 	}
// 
// 	/* Allocate matrix */
	wsUint N_col = N_common+(N_rare > 0);
// 	Ra_mat = sseMatrix(N_testSamp, N_col);
// 
// 	// m <- apply( as.matrix(x[,ic1]), 1, max )
// 	// reduced.x <- cbind(x[,ic2], m)
// 	j = 0;
// 	Ico = 0;
// 	Ica = N_ctrl;
// 	FOREACHDO (vector<int>::iterator, Xa_set, it, j++) {
// 		wsUint	N_idx	= *it;
// 		xSNP	&X_SNP	= Xa_SNP[N_idx];
// 
// 		for (i=0 ; i<N_samp ; i++) {
// 			if (Ba_misPheno[i]) continue;
// 
// 			wsUint X = Ra_Y[I]==1?Ica:Ico;
// 
// 			if (X_SNP.maf < OPT_REAL(raremaf)) {
// 				if (isAvailable(Na_data[i][N_idx]) && Na_data[i][N_idx])
// 					Ra_mat[X][N_common] = 1;
// 			} else
// 				Ra_mat[X][j] = (wsReal)Na_data[i][N_idx];
// 
// 			Ra_Y[I]==1?Ica++:Ico++;
// 			I++;
// 		}
// 	}

	/* If some geneset have multiple variants */
	if (N_col > 1) {
		wsReal **Ra_covCase = sseMcovRow(Ra_PEDCMCca, N_col, N_case);	/* 120927 */
		wsReal **Ra_covCtrl = sseMcovRow(Ra_PEDCMCct, N_col, N_ctrl);	/* 120927 */
		wsReal **Ra_s0 = Ra_covCtrl;
		wsReal *Ra_diffMean = NULL;	/* 120927 */
		sseMalloc(Ra_diffMean, wsReal, N_col);

#ifdef GSvalidate
		char S_buf[256];
		sprintf(S_buf, "%s.PEDCMC.covCase", S_gs);
		exportMatrix(S_buf, Ra_covCase, N_col, N_col);
		sprintf(S_buf, "%s.PEDCMC.covCtrl", S_gs);
		exportMatrix(S_buf, Ra_covCtrl, N_col, N_col);
		sprintf(S_buf, "%s.PEDCMC.case", S_gs);
		exportMatrix(S_buf, Ra_PEDCMCca, N_col, N_case);
		sprintf(S_buf, "%s.PEDCMC.ctrl", S_gs);
		exportMatrix(S_buf, Ra_PEDCMCct, N_col, N_ctrl);
#endif

		// cG <- (nG-1)*cov(co)
		sseMpC((wsMatCst)Ra_covCtrl, (wsReal)(N_ctrl-1), Ra_covCtrl, N_col);
		// cA <- (nA-1)*cov(ca)
		sseMpC((wsMatCst)Ra_covCase, (wsReal)(N_case-1), Ra_covCase, N_col);
		// cAG <- cA+cG
		sseMaM(Ra_covCtrl, Ra_covCase, Ra_covCtrl, N_col);
		sseUnmat(Ra_covCase, N_col);
		// s0 <- 1/(nA+nG-2) * cAG
		sseMpC((wsMatCst)Ra_covCtrl, W1/(wsReal)(N_ctrl+N_case-2),
			Ra_s0, N_col);

		// diff <- m1-m2	
		wsUint j;
		for (i=0 ; i<N_col ; i++) {
			wsReal R_meanCtrl = W0;
			wsReal R_meanCase = W0;

			for (j=0 ; j<N_ctrl ; j++)
				R_meanCtrl += Ra_PEDCMCct[i][j];
			for (j=0 ; j<N_case ; j++)
				R_meanCase += Ra_PEDCMCca[i][j];

			R_meanCtrl /= (wsReal)N_ctrl;
			R_meanCase /= (wsReal)N_case;
			Ra_diffMean[i] = R_meanCase-R_meanCtrl;
		}
		//	Ra_diffMean[i] = sseVsum(Ra_PEDCMCct[i], N_ctrl)/(wsReal)N_ctrl
		//		- sseVsum(Ra_PEDCMCca[i], N_case)/(wsReal)N_case;

#ifdef GSvalidate
		sprintf(S_buf, "%s.PEDCMC.diffMean", S_gs);
		exportMatrix(S_buf, &Ra_diffMean, 1, N_col);
#endif
		// Get the pseudoinverse of s0 and rank concurrently */
		wsUint	N_s0rank;
		wsReal	**Ra_s0inv = NULL;

		// rnk <- qr(s0)$rank
		// s0inv <- ginv(s0)
		SVDinverse(Ra_s0, N_col, &Ra_s0inv, &N_s0rank);
#ifdef GSvalidate
		sprintf(S_buf, "%s.PEDCMC.s0inv", S_gs);
		exportMatrix(S_buf, Ra_s0inv, N_col, N_col);
		wsReal R_s0rank = (wsReal)N_s0rank;
		wsReal *Rp_s0rank = &R_s0rank;
		sprintf(S_buf, "%s.PEDCMC.s0rank", S_gs);
		exportMatrix(S_buf, &Rp_s0rank, 1, 1);
#endif

		// t2 <- (nA*nG)/(nA+nG) * t(diff)%*%s0inv%*%diff
		wsReal R_t2 = (wsReal)(N_case*N_ctrl)/(wsReal)(N_case+N_ctrl) *
			sseVpMpV(Ra_diffMean, Ra_s0inv, N_col);
 		sseUnmat(Ra_s0, N_col);
		sseUnmat(Ra_s0inv, N_col);
		sseFree(Ra_diffMean);

//		sseMalloc(Ra_Yp, wsReal, N_sampGS);	/* 120927 */
		// Yp <- cc-nA/(nA+nG)
		cVector V_yP = *Vp_yGSunAdj - (wsReal)N_case/(wsReal)(N_case+N_ctrl);
//		sseVaC(Ra_yUnadjGS, -(wsReal)N_case/(wsReal)(N_case+N_ctrl), Ra_Yp,
//			N_sampGS);
		wsReal R_denom = V_mask.qf(V_yP, *((cSymMatrix *)Mp_phi));
#ifdef GSvalidate
		wsReal *Rp_denom = &R_denom;
		sprintf(S_buf, "%s.PEDCMC.denom", S_gs);
		exportMatrix(S_buf, &Rp_denom, 1, 1);
		sprintf(S_buf, "%s.PEDCMC.mask", S_gs);
		V_mask.file(S_buf);

		sprintf(S_buf, "%s.PEDCMC.yP", S_gs);
		V_yP.file(S_buf);
//		exportMatrix(S_buf, &Ra_Yp, 1, N_sampGS);
#endif
	// t2 <- t2*nA*nG/(nA+nG)/( t(Yp) %*% phi %*% Yp )
		R_t2 *= (wsReal)(N_case*N_ctrl)/(wsReal)(N_case+N_ctrl)/
			R_denom;
//			sseVpMpV(Ra_Yp, Ra_phiGS, N_sampGS, Ra_Yp, Ra_mask);
//		sseFree(Ra_Yp);

		if (Rp_stat) *Rp_stat = R_t2;
		R_pV = PVchisq(R_t2, N_s0rank);

		*Np_rank = N_s0rank;
		*Np_typeTest = 2;
		V_yP.rem();
	} else {
// 		if ((N_common+N_rare) == 1) {
// 			/* Currently, below is not reliable, do not calculate */
//  			*Rp_stat = WISARD_NA_REAL;
//  			R_pV = WISARD_NAN;
// 		} else {
			wsRealCst *Ra_yGS = Cp_IO->getPheno();

			/* Otherwise, using chi^2 */
			Ico = 0;
			Ica = 0;
			wsUint N_AA = 0, N_aA = 0;
			wsUint N_Aa, N_aa;
			for (i=I=0 ; i<N_samp ; i++) {
				if (Ba_filt[i]) continue;

				if (Ra_yGS[I] == WISARD_AFFECTED) {
					if (Ra_PEDCMCca[0][Ica++]) N_AA++;
				} else if (Ra_yGS[I] == WISARD_UNAFFECTED) {
					if (Ra_PEDCMCct[0][Ico++]) N_aA++;
				}
				I++;
			}
			N_Aa = N_case - N_AA;
			N_aa = N_ctrl - N_aA;

			R_pV = test2x2_PEDCMC(N_AA, N_aA, N_Aa, N_aa, TRUE, Np_typeTest);

			// Yp <- cc-nA/(nA+nG)
			cVector V_yP = *Vp_yGSunAdj - (wsReal)N_case/(wsReal)(N_case+N_ctrl);

			if (*Np_typeTest == 1) {
				/* In fisher case
				 * rlt <- qchisq(rlt, 2, lower.tail=F)
				 * x <- dd-nA/(nA+nG)
				 * rlt = rlt * nA*nG/(nA+nG)/( t(x)%*%phi%*%x )
				 */
				int N_which = 2, N_stat;
				double R_pValInp = 1.0 - (double)R_pV;
				double R_X, R_bound;
				double R_df = 2.0;
				double R_nc = 0.0;
				if (R_pValInp == 1.0) /* In case of invalid range */
					R_X = numeric_limits<wsReal>::infinity();
				else
					cdfchn(&N_which, &R_pValInp, NULL, &R_X, &R_df,
					&R_nc, &N_stat, &R_bound);
				R_X *= N_case*N_ctrl/(wsReal)(N_case+N_ctrl) /
					V_mask.qf(V_yP, *((cSymMatrix *)Mp_phi));
				R_pV = PVchisq(R_X, 2.0);
			} else {
				/* In chisq case
				 * x <- dd-nA/(nA+nG)
				 * rlt = rlt * nA*nG/(nA+nG)/( t(x)%*%phi%*%x )
				 */
				R_pV *= N_case*N_ctrl/(wsReal)(N_case+N_ctrl) /
					V_mask.qf(V_yP, *((cSymMatrix *)Mp_phi));
				R_pV = PVchisq(R_pV, 2.0);
			}
			*Np_rank = 2;

//			sseMalloc(Ra_Yp, wsReal, N_sampGS);	/* 120927 */
//			sseVaC(Ra_yUnadjGS, -(wsReal)N_case/(wsReal)N_sampGS, Ra_Yp,
//				N_sampGS);
// 			R_pV *= N_case*N_ctrl/(wsReal)(N_case+N_ctrl) /
// 				V_mask.qf(V_yP, *((cSymMatrix *)Mp_phi));
//				sseVpMpV(Ra_Yp, Ra_phiGS, N_sampGS, Ra_Yp, Ra_mask);
			sseFree(Ra_Yp);
//		}
	}

	return R_pV;
}

mGeneDef& cGSTestAnalysis::getGeneDef()
{ return Cp_anaGsm->getGeneDef(); }

wsMat cGSTestAnalysis::_makeCorr(wsUint *Np_szPhi)
{
	wsMat Ra_phi = NULL;

	if (IS_ASSIGNED(longitudinal)) {
		wsUint N_r, N_c;
		wsUintCst N_pheno = Cp_IO->sizePheno();

		/* Load matrix */
		Ra_phi = makeMatrixSSE(OPT_STRING(longitudinal), &N_r, &N_c);
		/* Dimension check */
		if (N_r!=N_c || N_r<=0 || N_pheno!=N_r)
			halt_fmt(WISARD_INVL_LONGIVARDIM, N_r, N_c, N_pheno, N_pheno);
		*Np_szPhi = N_r;
	} else if (OPT_ENABLED(indep)) {
		wsUint N_sz = getIO()->sizeSample();
		Ra_phi		= sseEmptyMatrix(N_sz, N_sz);
		for (wsUint i=0 ; i<N_sz ; i++)
			Ra_phi[i][i] = W1;
		*Np_szPhi	= N_sz;
	} else {
		Ra_phi		= getFullCorMat(Cp_anaCorr);
		*Np_szPhi	= getIO()->sizeSample();
	}

	return Ra_phi;
}

#endif

} // End namespace ONETOOL
