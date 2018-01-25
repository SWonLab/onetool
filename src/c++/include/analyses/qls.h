#pragma once
#ifndef __WISARD_QLS_H__
#define __WISARD_QLS_H__
#include "global/analysis.h"
#include "global/io.h"

namespace ONETOOL {

class cPPPAnalysisV2;
class cMatrix;
class cVector;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cQlsAnalysis : public cAnalysis
{
	cAnalysis		*Cp_anaCorr;
	cPPPAnalysisV2	*Cp_anaPPP;
	const char		*Ba_isInc;
	wsReal			*Ra_availY;
	wsUint			N_incSamp;
	void			_doTest(wsUint i, char **Na_data,
		cMatrix *Mp_phi, wsRealCst *Ra_ppp, cVector& V_yV, wsReal R_yVy,
		double &R_Tqls, double &R_Pqls);
public:
	cQlsAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP,
		cAnalysis *Cp_inpAnaCorr);
	~cQlsAnalysis();

	void		run();
};

#else

/* DUMMY DEFINITION */
class cQlsAnalysis : public cAnalysis
{
public:
	cQlsAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP,
		cAnalysis *Cp_inpAnaCorr) : cAnalysis(Cp_inpIO) {}
	~cQlsAnalysis() {}
	void		run() {}
};

#endif

} // End namespace ONETOOL

#endif
