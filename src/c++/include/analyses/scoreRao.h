#pragma once
#ifndef __WISARD_SCORERAO_H__
#define __WISARD_SCORERAO_H__
#include "global/analysis.h"

namespace ONETOOL {

class cEmAiAnalysisV2;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cTestRaoAnalysis : public cAnalysis
{
	cEmAiAnalysisV2	*Cp_anaEmAi;
	wsReal			R_Trs;
	wsReal			R_pValue;

//	wsReal		*Ra_snp;
	cVector		*Vp_XtP;
	wsReal		**Ra_phi;
//	wsReal		*Ra_residual;
	cVector		*Vp_RtP;
//	wsReal		**Ra_designMat;
//	wsReal		**Ra_Xt;
	cStdMatrix	*Mp_WtP;
	wsUint		N_samp, N_colX;
	cVector		*Vp_YtP;
//	wsReal		*Ra_pheno;
	cDiagMatrix	*Mp_Vinv;
	cStdMatrix	*Mp_EVX;
public:
//	static wsReal	**Ra_XVIXinv;
	static cSymMatrix	M_WViWi;

	cTestRaoAnalysis(cIO *Cp_IO, wsUintCst N_inpSamp,
		cVector &V_inpYtP, wsReal **Ra_inpPhi, cEmAiAnalysisV2 *Cp_inpAnaEmAi,
		cStdMatrix *Mp_inpWtP,
		wsUintCst N_inpColDesignMat, cVector &V_inpRtP,
		cDiagMatrix &M_inpInvV, cStdMatrix *Mp_inpEVX, cVector *Vp_inpXtP);
	~cTestRaoAnalysis();
	void run();
	wsReal getPvalue() { return R_pValue; }
	wsReal getStatistics() { return R_Trs; }
};

#else

/* DUMMY DEFINITION */
class cTestRaoAnalysis : public cAnalysis
{
public:
	//	static wsReal	**Ra_XVIXinv;
	static cSymMatrix	M_WViWi;

	cTestRaoAnalysis(cIO *Cp_IO, wsUintCst N_inpSamp,
		cVector &V_inpYtP, wsReal **Ra_inpPhi, cEmAiAnalysisV2 *Cp_inpAnaEmAi,
		cStdMatrix *Mp_inpWtP,
		wsUintCst N_inpColDesignMat, cVector &V_inpRtP,
		cDiagMatrix &M_inpInvV, cStdMatrix *Mp_inpEVX, cVector *Vp_inpXtP)
		: cAnalysis(Cp_IO) {}
	~cTestRaoAnalysis() {}
	void run() {}
	wsReal getPvalue() { return WISARD_NAN; }
	wsReal getStatistics() { return WISARD_NAN; }
};

#endif

} // End namespace ONETOOL

#endif
