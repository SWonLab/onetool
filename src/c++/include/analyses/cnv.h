#ifndef __WISARD_CNV_H__
#define __WISARD_CNV_H__

#pragma once
#include "score.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cCnvAnalysis : public cTestScoreAnalysis
{
	wsUint	N_cnvClus;
	void _doCnv(wsReal **Ra_CNV, wsUint N_szCnv);
public:
	cCnvAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaPhi);
	~cCnvAnalysis();

	void run();
};

#else

/* DUMMY DEFINITION */
class cCnvAnalysis : public cTestScoreAnalysis
{
	wsUint	N_cnvClus;
	void _doCnv(wsReal **Ra_CNV, wsUint N_szCnv);
public:
	cCnvAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaPhi) : cTestScoreAnalysis(Cp_inpIO, NULL, 0) {}
	~cCnvAnalysis() {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
