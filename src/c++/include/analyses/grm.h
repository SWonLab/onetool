#pragma once
#include "global/analysis.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cPDDTAnalysis;

int corrTest(double R_corr, double R_pddt2, wsUint N_sz, double R_alpha, FILE *H_detOut);

class cGRMAnalysis : public cAnalysis
{
	cPDDTAnalysis	*Cp_anaPDDT;

public:
	cGRMAnalysis(cIO *Cp_inpIO, cPDDTAnalysis *Cp_inpAnaPDDT);
	~cGRMAnalysis();
	static int	sort(const void *p, const void *q);
	void		findSampleRemovals(wsReal **Rp_haveError, int N_threshold);
	
	void		run();
};

#else

/* DUMMY DEFINITION */
class cPDDTAnalysis;

int corrTest(double R_corr, double R_pddt2, wsUint N_sz, double R_alpha, FILE *H_detOut);

class cGRMAnalysis : public cAnalysis
{
	cPDDTAnalysis	*Cp_anaPDDT;

public:
	cGRMAnalysis(cIO *Cp_inpIO, cPDDTAnalysis *Cp_inpAnaPDDT) : cAnalysis(Cp_inpIO) {}
	~cGRMAnalysis() {}
	static int	sort(const void *p, const void *q) { return 0; }
	void		findSampleRemovals(wsReal **Rp_haveError, int N_threshold) {}

	void		run() {}
};

#endif

} // End namespace ONETOOL
