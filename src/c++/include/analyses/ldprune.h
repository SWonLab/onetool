#ifndef __WISARD_LDPRUNE_H__
#define __WISARD_LDPRUNE_H__

#pragma once
#include "global/analysis.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cLdPruneAnalysis : public cAnalysis
{
	vector<xUintSort> Xm_mmap[MAX_NCHR + 1];
	wsUint N_pruneLdWin;
	wsUint N_pruneStep;
	wsReal R_pruneThr;
public:
	cLdPruneAnalysis(cIO *Cp_inpIO);
	~cLdPruneAnalysis();
	char*		vif_prune(cSymMatrix& m, double threshold, vInt& Nv_idxVrt);
	cSymMatrix	calcSetCovarianceMatrix(vInt Nv_idx);

	void run();
};

#else

/* DUMMY DEFINITION */
class cCnvAnalysis : public cTestScoreAnalysis
{
	vector<xUintSort> Xm_mmap[MAX_NCHR + 1];
	wsUint N_pruneLdWin;
	wsUint N_pruneStep;
	wsReal R_pruneThr;
public:
	cLdPruneAnalysis(cIO *Cp_inpIO) {}
	~cLdPruneAnalysis() {}
	vBool		vif_prune(cSymMatrix& m, double threshold, vInt& Nv_idxVrt) { return vBool(); }
	cSymMatrix	calcSetCovarianceMatrix(vInt Nv_idx) { return cSymMatrix(); }
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
