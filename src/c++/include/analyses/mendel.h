#pragma once
#ifndef __WISARD_MENDEL_H__
#define __WISARD_MENDEL_H__
#include "global/analysis.h"

namespace ONETOOL {

typedef struct _xMendelStat {
	map<string,__int64>	Xm_fam;
	vInt				Xv_samp;
	vInt				Xv_marker;
	map<string,__int64>	Xm_ffam;
	vInt				Xv_fsamp;
	vInt				Xv_fmarker;
} xMendelStat;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cFamStrAnalysis;
class cMendelAnalysis : public cAnalysis
{
	xMendelStat				X_t;
	cFamStrAnalysis*	Cp_anaFS;
public:
	cMendelAnalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS);
	~cMendelAnalysis();

	void run();
	xMendelStat& getMendelStat() { return X_t; }
};

#else

/* DUMMY DEFINITION */
class cMendelAnalysis : public cAnalysis
{
public:
	cMendelAnalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS) : cAnalysis(Cp_inpIO) {}
	~cMendelAnalysis() {}
	xMendelStat& getMendelStat() { xMendelStat u; return u; }
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
