#ifndef __WISARD_LD_H__
#define __WISARD_LD_H__
#pragma once
#include "global/analysis.h"

namespace ONETOOL {

class cVariantMap;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cLdAnalysis : public cAnalysis
{
	vector<vInt>	Na_SE;
	wsUint			*Na_sset, N_set;
	cVariantMap		*Cp_anaMM;
public:
	cLdAnalysis(cIO *Cp_inpIO, wsUint *Np_sset=NULL, wsUint N_size=0);
	cLdAnalysis(cIO *Cp_inpIO, cVariantMap *Cp_inpAnaMM);
	~cLdAnalysis();
	void run();
};

#else

/* DUMMY DEFINITION */
class cLdAnalysis : public cAnalysis
{
public:
	cLdAnalysis(cIO *Cp_inpIO, wsUint *Np_sset=NULL, wsUint N_size=0) :
		cAnalysis(Cp_inpIO) {}
	cLdAnalysis(cIO *Cp_inpIO, cVariantMap *Cp_inpAnaMM) :
		cAnalysis(Cp_inpIO) {}
	~cLdAnalysis() {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
