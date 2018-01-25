#pragma once
#ifndef __WISARD_MARKER_H__
#define __WISARD_MARKER_H__

#include "global/analysis.h"

namespace ONETOOL {

class cFamStrAnalysis;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cVariantAnalysis : public cAnalysis
{
	cFamStrAnalysis *Cp_anaFS;
public:
	cVariantAnalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS);
	~cVariantAnalysis();
	void run();
};

#else

/* DUMMY DEFINITION */

class cVariantAnalysis : public cAnalysis
{
public:
	cVariantAnalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS) : cAnalysis(Cp_inpIO) {}
	~cVariantAnalysis() {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif