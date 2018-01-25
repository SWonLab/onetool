#pragma once
#ifndef __WISARD_TRIDGE_H__
#define __WISARD_TRIDGE_H__
#include "score.h"

namespace ONETOOL {

typedef enum _xSeqType {
	SEQ_ASIS,
	SEQ_LOG,
	SEQ_EXP
} xSeqType;

cVector seq(wsReal R_s, wsReal R_e, wsUint N_len, xSeqType X_in, xSeqType X_out);

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cTridgeAnalysis : public cAnalysis
{
public:
	cTridgeAnalysis(cIO *Cp_inpIO);
	~cTridgeAnalysis();

	void run();
};

#else

/* DUMMY DEFINITION */
class cTridgeAnalysis : public cAnalysis
{
public:
	cTridgeAnalysis(cIO *Cp_inpIO) : cAnalysis(Cp_inpIO) {}
	~cTridgeAnalysis() {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
