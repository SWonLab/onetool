#ifndef __WISARD_POWERCALC_H__
#define __WISARD_POWERCALC_H__
#pragma once

#include "global/analysis.h"
#include "global/io.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cPowerCalcAnalysis : public cAnalysis
{
	cAnalysis *Cp_anaCorr;
	void	_getVinv(wsUint N_samp, wsUint N_t, wsReal R_rho,
		wsReal R_sigg, wsReal R_sigc, wsReal *Rp_1, wsReal *Rp_G,
		wsReal *Rp_C, wsReal *Rp_0);
	void	_getSinglePheno(wsUint N_perm);
	void	_getMultiPheno(wsUint N_t, wsUint N_s, wsUint N_e=0xffffffff);
public:
	cPowerCalcAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaCorr);
	~cPowerCalcAnalysis();
	void run();
};

class cPowerCalcAnalysisV2 : public cAnalysis
{
	cAnalysis *Cp_anaCorr;
	void	_getSinglePheno(wsUint N_samp);
	void	_getMultiPheno(wsUint N_t, wsUint N_s, wsUint N_e=0xffffffff);
public:
	cPowerCalcAnalysisV2(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaCorr);
	~cPowerCalcAnalysisV2();
	void run();
};

#else

/* DUMMY DEFINITION */
class cPowerCalcAnalysis : public cAnalysis
{
public:
	cPowerCalcAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaCorr) : cAnalysis(Cp_inpIO) {}
	~cPowerCalcAnalysis() {}
	void run() {}
};

/* DUMMY DEFINITION */
class cPowerCalcAnalysisV2 : public cAnalysis
{
public:
	cPowerCalcAnalysisV2(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaCorr) : cAnalysis(Cp_inpIO) {}
	~cPowerCalcAnalysisV2() {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
