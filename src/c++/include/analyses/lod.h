#ifndef __WISARD_LOD_H__
#define __WISARD_LOD_H__
#pragma once
#include "global/analysis.h"

namespace ONETOOL {

class cFamStrAnalysis;
class cNucFam;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cLodAnalysis : public cAnalysis
{
	cSetManagerAnalysis*	Cp_anaSM;
	mNucFam					Xm_curNF;
	cFamStrAnalysis*		Cp_anaFS;
public:
	cLodAnalysis(cIO* Cp_inpIO, cFamStrAnalysis* Cp_inpAnaFS,
		cSetManagerAnalysis* Cp_inpAnaSM);
	void run();
};

#else

/* DUMMY DEFINITION */
class cLodAnalysis : public cAnalysis
{
public:
	cLodAnalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS,
		cSetManagerAnalysis* Cp_inpAnaSM) :
		cAnalysis(Cp_inpIO) {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif