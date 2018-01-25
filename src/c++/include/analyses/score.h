#pragma once
#ifndef __WISARD_SCORE_H__
#define __WISARD_SCORE_H__
#include "global/analysis.h"
#include "utils/matrix.h"
#include "utils/vector.h"

namespace ONETOOL {

class cEmAiAnalysisV2;
class cFamStrAnalysis;

/* !!! FARVAT-specific analysis !!! */
#if ((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE & FUNCTION_FARVAT) || \
	(TOOLSET_TYPE & FUNCTION_MFARVAT) || (TOOLSET_TYPE & FUNCTION_FARVATX)

class cTestScoreAnalysis : public cAnalysis
{
	cEmAiAnalysisV2	*Cp_res0EmAi;
	cAnalysis		*Cp_anaPhi;
	cFamStrAnalysis *Cp_anaFS;
	wsReal			**Ra_phi;
//	wsReal			**Ra_phiInv;
	wsUint			N_testSamp;
	char			B_isREML;
	char			*Ba_excl;
//	wsUint			*Na_szFam;
//	wsUint			*Na_sampIdx;
//	wsReal			**Ra_Y;
	wsReal			*Ra_stat;
	cStdMatrix		*Mp_Yt;
	cStdMatrix		M_YtP;
//	wsReal			**Ra_Vinv;
	cDiagMatrix		*Mp_Vinv;
//	wsReal			**Ra_X, **Ra_Xt;
	cStdMatrix		M_Xt;
	cStdMatrix		M_WtP;
	cVector			V_residP;
	cSymMatrix		M_WViWi;
//	wsReal			*Ra_residual;
	cStdMatrix		*Mp_EVX;
	cDiagMatrix		*Mp_EV;

	cStdMatrix		**Ma_famWtP;
	cDiagMatrix		*Ma_famVinv;
	wsUint*			Na_seq;			/* Family-wise sorted sequence */

	void _test0();
public:
	cTestScoreAnalysis(cIO *Cp_IO, cAnalysis *Cp_inpAnaPhi,
		cFamStrAnalysis *Cp_inpAnaFamStr, char B_inpIsREML=1);
	~cTestScoreAnalysis();

	void			run();
	wsUintCst			getTestSampSize() { return N_testSamp; }
//	wsUint*			getMapSample() { return Na_sampIdx; }
	const char*		getExcl() { return Ba_excl; }
//	wsUint*			getFamSize() { return Na_szFam; }
	cEmAiAnalysisV2*
					getNullEmAi() { return Cp_res0EmAi; }
//	wsReal**		getDesignMatT() { return Ra_Xt; }
	cStdMatrix&		getWtP() { return M_WtP; }
	cStdMatrix&		getYtP() { return M_YtP; }
//	wsReal**		getY() { return Ra_Y; }
//	wsReal**		getVinv() { return Ra_Vinv; }
	cDiagMatrix&	getVinv() { return *Mp_Vinv; }
	cStdMatrix*		getEVX() { return Mp_EVX; }
	cSymMatrix&		getWViWi() { return M_WViWi; }
	wsReal*			getTs() { return Ra_stat; }
//	wsReal**		getPhi() { return Ra_phi; }
//	wsReal*			getResidual() { return Ra_residual; }
	cVector&		getResidP() { return V_residP; }
	cStdMatrix**	getFamWtP() { return Ma_famWtP; }
	cDiagMatrix*	getFamVinv() { return Ma_famVinv; }
	wsUint*			getSampSeqByFam() { return Na_seq; }
};

#else

/* DUMMY DEFINITION */

class cTestScoreAnalysis : public cAnalysis
{
public:
	cTestScoreAnalysis(cIO *Cp_IO, cAnalysis *Cp_inpAnaPhi,
		cFamStrAnalysis *Cp_inpAnaFamStr, char B_inpIsREML=1) : cAnalysis(Cp_IO) {}
	~cTestScoreAnalysis() {}

	void			run() {}
	wsUintCst			getTestSampSize() { return 0; }
	wsUint*			getMapSample() { return NULL; }
	wsUint*			getFamSize() { return NULL; }
	cEmAiAnalysisV2*
		getNullEmAi() { return NULL; }
	//	wsReal**		getDesignMatT() { return Ra_Xt; }
	cStdMatrix&		getWtP() { return *(new cStdMatrix()); }
	cStdMatrix&		getYtP() { return *(new cStdMatrix()); }
	//	wsReal**		getY() { return Ra_Y; }
	//	wsReal**		getVinv() { return Ra_Vinv; }
	cDiagMatrix&	getVinv() { return *(new cDiagMatrix()); }
	cStdMatrix*		getEVX() { return NULL; }
	cSymMatrix&		getWViWi() { return *(new cSymMatrix()); }
	//	wsReal**		getPhi() { return Ra_phi; }
	//	wsReal*			getResidual() { return Ra_residual; }
	cVector&		getResidP() { return *(new cVector()); }
	cStdMatrix**	getFamWtP() { return NULL; }
	cDiagMatrix*	getFamVinv() { return NULL; }
};

#endif

} // End namespace ONETOOL

#endif
