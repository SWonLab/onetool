#pragma once
#include "global/analysis.h"
#include "global/io.h"
#include "emai.h"

namespace ONETOOL {

wsVec _buildFamI(cIO* Cp_IO, xFamily &X_fam, wsUint N_szFam,
	wsUint *Na_idxMem, wsVec Ra_adjY, wsSym Ra_cor, wsVec Ra_newI_a,
	wsVecCst Ra_fullY, wsReal R_h2);
wsVec _buildFamII(cIO* Cp_IO, xFamily &X_fam, wsUint N_szFam,
	wsUint* Na_idxMem, wsVec Ra_adjY, wsSym Ra_cor, wsVec Ra_newII_a,
	wsVecCst Ra_fullY, wsReal R_h2);

class cEmAiAnalysisV2;
wsVec	genAdjPheno(cIO *Cp_IO, cEmAiAnalysisV2 *Cp_anaEmAi, cVector& V_prev,
	wsReal **Ra_adjFullY_mqls);
void	genFqlsPheno(cIO* Cp_IO, wsVec Ra_adjFullY_fqls, wsSym Ra_offsetCor,
	wsVec* Ra_newI_a, wsVec* Ra_newII_a);

typedef struct _resPAformula
{
	wsReal *Ra_muHat;
	wsReal **Ra_varHat;
} resPAfmlar;

class cPPPAnalysisV2;
class cFamStrAnalysis;

int fetchProband(xSample *Xp_s, void *Vp_data);
int fetchFamMem(xSample *Xp_s, void *Vp_data);

/* !!! FQLS-specific analysis !!! */
#if (TOOLSET_TYPE & FUNCTION_MFQLS) || ((TOOLSET_TYPE & 0x100) == 0x100)

class cVector;
class cStdMatrix;
class cEmAiAnalysisV2;

class cFamQlsAnalysis : public cAnalysis
{
	cAnalysis		*Cp_anaCorr;
	cFamStrAnalysis	*Cp_anaFam;
	cEmAiAnalysisV2	*Cp_anaEmAi;
	cPPPAnalysisV2	*Cp_anaPPP;
	const char		*Ba_incSamp;
	wsUint			N_nonMissSamp;
	wsReal			R_thr, R_prev;
	char			B_doFqls;

	cStdMatrix		_getAdjPheno(cVector& V_prev, cStdMatrix **Mp_adjFYmqls);
	int				_doTest(wsUint i, xVariant &X_snp, wsReal *Rp_TmQLS,
		wsReal *Rp_TfqlsI, wsReal *Rp_TfqlsII, char *Ba_misGeno,
		char **Na_data, wsReal *Ra_adjFullY, wsReal *Ra_adjFullY_mqls,
		wsReal **Ra_fullCor, wsReal *Ra_full_Hrow, wsReal *Ra_newI_a, wsReal *Ra_newII_a,
		wsReal **Ra_full_phiInv, wsReal *Ra_full_1phiInv,
		wsReal R_full_1phiInv1_inv, wsReal **Ra_full_Vx, char B_availOnly);
public:
	cFamQlsAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP,
		cFamStrAnalysis *Cp_inpAnaFam,
		cAnalysis *Cp_inpAnaCorr, cEmAiAnalysisV2* Cp_inpAnaNullEmAi);
	~cFamQlsAnalysis();

	void		run();
};

#else

/* DUMMY DEFINITION */
class cFamQlsAnalysis : public cAnalysis
{
public:
	cFamQlsAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP,
		cFamStrAnalysis *Cp_inpAnaFam,
		cAnalysis *Cp_inpAnaCorr, cEmAiAnalysisV2* Cp_inpAnaNullEmAi) : cAnalysis(Cp_inpIO) {}
	~cFamQlsAnalysis() {}
	void		run() {}
};

#endif

} // End namespace ONETOOL
