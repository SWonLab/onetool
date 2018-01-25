#pragma once
#ifndef __WISARD_SCOREGEN_H__
#define __WISARD_SCOREGEN_H__
#include "global/analysis.h"

namespace ONETOOL {

class cEmAiAnalysisV2;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cTestGenScoreAnalysis : public cAnalysis
{
	cEmAiAnalysisV2
				*Cp_anaEmAi;
	wsReal		R_Tgs;
	wsReal		R_pValue;

//	wsReal		*Ra_snp;
	wsReal		**Ra_phi;
//	wsReal		*Ra_residual;
	cVector		*Vp_RtP;
	cStdMatrix	*Mp_WtP;

	wsUint		*Na_szFam;
	wsUint		N_samp, Ra_colX;
//	wsReal		*Ra_pheno;
	cVector		*Vp_YtP;

//	wsReal		**Ra_invV;
	cVector		*Vp_XtP;
	cDiagMatrix	*Mp_invV;
	cStdMatrix	*Mp_EVX;
public:
	cTestGenScoreAnalysis(cIO *Cp_IO, wsUintCst N_inpSamp,
		cVector	&V_YtP,
		wsReal **Ra_inpPhi, wsUint *Na_inpSzFam, cEmAiAnalysisV2 *Cp_inpAnaEmAi,
		cStdMatrix *Mp_inpWtP, wsUintCst N_colDesignMat,
		cVector &V_inpRtP,
//		wsReal *Ra_inpResidual,
		cDiagMatrix &M_inpInvV, cStdMatrix *Mp_inpEVX, cVector *Vp_inpXtP);
	~cTestGenScoreAnalysis();
	void run();
	wsReal getPvalue() { return R_pValue; }
	wsReal getStatistics() { return R_Tgs; }
};

#else

/* DUMMY DEFINITION */
class cTestGenScoreAnalysis : public cAnalysis
{
public:
	cTestGenScoreAnalysis(cIO *Cp_IO, wsUintCst N_inpSamp,
		cVector	&V_YtP,
		wsReal **Ra_inpPhi, wsUint *Na_inpSzFam, cEmAiAnalysisV2 *Cp_inpAnaEmAi,
		cStdMatrix *Mp_inpWtP, wsUintCst N_colDesignMat,
		cVector &V_inpRtP,
		cDiagMatrix &M_inpInvV, cStdMatrix *Mp_inpEVX, cVector *Vp_inpXtP)
		: cAnalysis(Cp_IO) {}
	~cTestGenScoreAnalysis() {}
	void run() {}
	wsReal getPvalue() { return WISARD_NAN; }
	wsReal getStatistics() { return WISARD_NAN; }
};

#endif

} // End namespace ONETOOL

#endif
