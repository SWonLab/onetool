#pragma once
#ifndef __WISARD_WINDOW_H__
#define __WISARD_WINDOW_H__
#include "global/analysis.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cWindowAnalysis : public cAnalysis
{
	wsUint **Na_bin;
	wsUint **Na_trs;
	wsUint *Na_szChr;
	wsUint *Na_szBin;
	wsUint N_szWindow;
public:
	cWindowAnalysis(cIO *Cp_inpIO);
	~cWindowAnalysis();
	void run();
};

#else

/* DUMMY DEFINITION */
class cWindowAnalysis : public cAnalysis
{
public:
	cWindowAnalysis(cIO *Cp_inpIO) : cAnalysis(Cp_inpIO) {}
	~cWindowAnalysis() {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
