#pragma once
#ifndef __WISARD_REGR_H__
#define __WISARD_REGR_H__

#include "global/analysis.h"

#define USE_LINEAR2

namespace ONETOOL {

typedef enum _xRegrT {
	REGR_LINEAR,
//#ifdef USE_LINEAR2
	REGR_LINEAR2, // Fast version
//#endif
	REGR_LOGISTIC
} xRegrT;

typedef struct _xAnaRegr {
	xRegrT	X_type;
	wsUint	N_pheno;
	wsUint	N_samp;
	wsUint	N_col;
	wsUint	N_gxe;
	char	B_gInv;
	short	N_iter;
	wsMat	Ra_Xt;
	wsMat	Ra_Ys;
	wsMat	Ra_t;
	wsMat	Ra_b;
	wsVec	Ra_swgt;

#ifdef USE_LINEAR2
	wsMat	Ra_Xy;		// _linear2
	wsMat	Ra_XX;		// _linear2
	wsUint	N_idxSNP;	// _linear2
#endif
} xAnaRegr;

class cSetManagerAnalysis;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cRegrAnalysis : public cAnalysis
{
	wsUint		N_pheno;
	char		B_isCont;	///< Whether the phenotype to be analyzed is continuous(1) or binary(0)
//	wsReal		*Ra_anaY;	///< Array[#tsamp] phenotype value
	wsMat		Ra_anaYs;	///< Array[#tsamp] phenotype value
	wsMat		Ra_anaXt;	///< Matrix[#cov+2,#tsamp] design matrix
	wsUint		N_anaSamp;	///< [=#tsamp] Sample to be included
	wsUint		N_col;
	xRegrT		X_type;
	const char	*Ba_isExcl;	///< Array[#samp] Excluded from the analysis(1) or not(0)
	cSetManagerAnalysis *Cp_anaSM;
	wsUint		_makeY(const char *Ba_filter, wsReal **Ra_resY);
	wsUint		_makeYs(const char *Ba_filter, wsMat *Ra_resY);
	wsUint		_makeX(const char *Ba_filter, wsMat *Ra_resXt);

	/* Single-thread linear regression */
	void		_regrLinear(cExporter &C_regr);
	/* Single-thread logistic regression */
	void		_regressionMain(cExporter &C_regr, wsMat Ra_pvals, int** Na_chrss,
		wsUint *Na_pval, int N_chr=-1);
	void		_doNull(wsMat Ra_curXt);
public:
	cRegrAnalysis(cIO* Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaMgr);
	~cRegrAnalysis();
	void		run();
	static int	fitLogisticByNR(wsUint N_curSamp, wsUint N_col,
		wsReal **Ra_curXt, wsReal *Ra_curY, wsReal *Ra_coef, wsReal *Ra_V,
		wsUint N_maxIter, bool &B_conv, wsReal **Rp_phiHat=NULL,
		wsReal *Rp_sigma=NULL);
//	wsReal*		getAnaY() { return Ra_anaY; }
	wsMat		getAnaYs() { return Ra_anaYs; }
	xRegrT		getType() { return X_type; }
	void		_getX(wsReal ***Ra_retAnaXt);
	wsUintCst		_setX(xVariant &X_snp, wsUint N_idxSNP, wsMat Ra_anaXt,
		wsReal *Ra_anaY, wsMat *Rp_curXt, wsReal **Rp_curY, wsUint *Np_imp);
	wsUintCst		_setXs(xVariant &X_snp, wsUint N_idxSNP, wsMat Ra_anaXt,
		wsMat Ra_anaYs, wsMat *Rp_curXt, wsMat *Rp_curYs, wsUint *Np_imp,
		wsUint N_idxGxG=0xffffffff);
	wsUintCst		_setXSs(vInt& X_vIdx, wsMat Ra_anaXt,
		wsMat Ra_anaYs, wsMat *Rp_curXt, wsMat *Rp_curYs, wsUint *Np_imp);
	wsUintCst		_setXXs(xVariant &X_s1, xVariant &X_s2, wsUint N_x1,
		wsUint N_x2, wsMat Ra_anaXt, wsMat Ra_anaYs, wsMat *Rp_curXt,
		wsMat *Rp_curYs, wsUint *Np_imp);
	wsUintCst		getBetaSize() { return N_col; }
	wsUintCst		getAnaSampSize() { return N_anaSamp; }
};

#else

/* DUMMY DEFINITION */
class cRegrAnalysis : public cAnalysis
{
public:
	cRegrAnalysis(cIO* Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaMgr) : cAnalysis(Cp_inpIO) {}
	~cRegrAnalysis() {}
	void		run() {}
	static int	fitLogisticByNR(wsUint N_curSamp, wsUint N_col,
		wsReal **Ra_curXt, wsReal *Ra_curY, wsReal *Ra_coef, wsReal *Ra_V,
		wsUint N_maxIter, bool &B_conv, wsReal **Rp_phiHat=NULL,
		wsReal *Rp_sigma=NULL) { return 0; }
	wsMat		getAnaYs() { return NULL; }
	void		_getX(wsReal ***Ra_retAnaXt) { }
	wsUintCst		_setX(xVariant &X_snp, wsUint N_idxSNP, wsMat Ra_anaXt,
		wsReal *Ra_anaY, wsMat *Rp_curXt, wsReal **Rp_curY, wsUint *Np_imp) { return 0; }
	wsUintCst		_setXs(xVariant &X_snp, wsUint N_idxSNP, wsMat Ra_anaXt,
		wsMat Ra_anaYs, wsMat *Rp_curXt, wsMat *Rp_curYs, wsUint *Np_imp) { return 0; }
	wsUintCst		_setXSs(vInt& X_vIdx, wsMat Ra_anaXt,
		wsMat Ra_anaYs, wsMat *Rp_curXt, wsMat *Rp_curYs, wsUint *Np_imp);
	wsUintCst		_setXXs(xVariant &X_snp, wsUint N_x1, wsUint N_x2, wsMat Ra_anaXt,
		wsMat Ra_anaYs, wsMat *Rp_curXt, wsMat *Rp_curYs, wsUint *Np_imp) { return 0; }
	wsUintCst		getBetaSize() { return 0; }
	wsUintCst		getAnaSampSize() { return 0; }
};

#endif

} // End namespace ONETOOL

#endif
