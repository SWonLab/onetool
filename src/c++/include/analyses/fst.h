#ifndef __WISARD_FST_H__
#define __WISARD_FST_H__
#pragma once

#include "global/analysis.h"
#include "utils/vector.h"

namespace ONETOOL {

class cSetManagerAnalysis;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cFstAnalysis : public cAnalysis
{
	cSetManagerAnalysis*	Cp_anaSM;
//	cAnalysis*				Cp_anaCor;
	cVector					V_varSamp;
public:
	cFstAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSM);
	~cFstAnalysis();
	void run();
};

#else

/* DUMMY DEFINITION */
class cFstAnalysis : public cAnalysis
{
public:
	cFstAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSM) : cAnalysis(Cp_inpIO) {}
	~cFstAnalysis() {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif