#pragma once
#ifndef __WISARD_DESIGN_H__
#define __WISARD_DESIGN_H__

#include "global/analysis.h"
#include "utils/vector.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cStudyAnalysis : public cAnalysis
{
	cAnalysis	*Cp_anaCor;
	cVector		V_varSamp;
public:
	cStudyAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaCor);
	~cStudyAnalysis();
	wsReal _exhaustive(wsUint **Na_Lcomb, wsUint *Np_comb, wsUint N_combReq,
		wsUint N_samp, wsMat Ra_cor);
	wsReal _exp(wsUint N_szComb, wsUint N_samp,
		wsSym Ra_cor, xRealSort *Xa_sort, wsReal *Rp_expSS);
	wsReal _stepwise(wsUint **Np_comb, wsUint *Np_retSize, wsUint N_samp2rem,
		wsUint N_samp, wsSym Ra_cor, xRealSort *Xa_sort);
	void run();
};

#else

/* DUMMY DEFINITION */
class cStudyAnalysis : public cAnalysis
{
public:
	cStudyAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaCor) : cAnalysis(Cp_inpIO) {}
	~cStudyAnalysis() {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
