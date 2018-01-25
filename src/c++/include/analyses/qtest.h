#pragma once
#ifndef __WISARD_QTEST_H__
#define __WISARD_QTEST_H__
#include "global/analysis.h"
#include "global/io.h"

namespace ONETOOL {

typedef struct _xRegression
{
	wsReal	*Ra_residual;
	wsReal	**Ra_covBeta;
	wsReal	*Ra_hatBeta;
} xRegression;

class cSetManagerAnalysis;
class cPPPAnalysisV2;
class cFemmaAnalysis;
class cSymMatrix;
class cVector;

/* !!! WISARD-specific analysis !!! */
#if ((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE & FUNCTION_QTEST)

class cQtestAnalysis : public cAnalysis
{
	cFemmaAnalysis
				*Cp_anaFemma;
	cPPPAnalysisV2
				*Cp_anaPPP;
	cSetManagerAnalysis
				*Cp_anaGsm;
	wsUint		N_anaSamp, N_oriSamp;
	xRegression *Xp_0;
	wsReal		*Ra_anaY;
	wsReal		R_stt;
	char		*Ba_misPheno;
	wsUint		*Na_null;

	wsReal**		_getClumpedGeno(wsStrCst S_gname, vInt &Xa_snp, char **Na_data,
		wsUint *Np_finMarker, wsReal **Rp_W, char *Ba_filt=NULL);
	void			_doTest(wsStrCst Sp_gname, xRegression *Xp_0,
		vInt &Xa_set, wsRealCst *Ra_origY, wsUint N_samp, char **Na_data,
		wsReal *Ra_pVals, wsUint *Np_szFin, wsUint *Np_sampFin,
		wsReal *Rp_zQ1, wsReal *Rp_zQ2);
	wsReal*			_makeWeight();
	wsReal			_getWeight(wsUint N_idxSnp);
	wsReal			_getWeight(vInt& Na_idxSnps);
	wsReal			_decomp(cSymMatrix &M_v, cVector &V_beta1, wsReal *Rp_Q2);
	xRegression*	regression(wsUint N_obs, wsReal *Ra_Y, wsUint N_var,
		wsReal **Ra_X);
public:
	cQtestAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP,
		cSetManagerAnalysis *Cp_inpGsmAna, cFemmaAnalysis *Cp_inpAnaFemma=NULL);
	~cQtestAnalysis();
	void run();
};

#else

/* DUMMY DEFINITION */
class cQtestAnalysis : public cAnalysis
{
public:
	cQtestAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP,
		cSetManagerAnalysis *Cp_inpGsmAna, cFemmaAnalysis *Cp_inpAnaFemma=NULL) : cAnalysis(Cp_inpIO) {}
	~cQtestAnalysis() {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
