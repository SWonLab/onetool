#pragma once
#include "../analysis.h"

class cPPPAnalysisV2;

/* !!! WISARD-specific analysis !!! */
#if TOOLSET_TYPE == TOOLSET_WISARD

class cFinalStatAnalysisV2 : public cAnalysis
{
	cAnalysis		*Cp_anaCorr;	///< Pointer to correlation analysis instance
	cPPPAnalysisV2	*Cp_anaPPP;		///< Pointer to PPP analysis instance
	REAL_t			*Ra_stat;		///< Array[size #SNP], Final statistics
public:
	REAL_t			*Ra_stdizeY;	///< Array[size #samp], Standardized by the formula Y-E(Y)

	cFinalStatAnalysisV2(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP, cAnalysis *Cp_inpAnaCorr);
	~cFinalStatAnalysisV2();
	void			run();

	/**
	 * cFinalStatAnalysis::getCorr Get correlation matrix from Cp_anaCorr
	 *
	 * @return    (REAL_t**)
	 */
	cAnalysis*		getCorrAnalysis() { return Cp_anaCorr; }

	/**
	 * cFinalStatAnalysis::getPPPAnalysis Get the pointer of PPP analysis instance
	 *
	 * @return    (cPPPAnalysis*)
	 */
	cPPPAnalysisV2*	getPPPAnalysis();

	static int		calcFinalStat(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result);
};

#else

/* DUMMY DEFINITION */
class cFinalStatAnalysisV2 : public cAnalysis
{
public:
	cFinalStatAnalysisV2(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP, cAnalysis *Cp_inpAnaCorr) :
		cAnalysis(Cp_inpIO) {}
	~cFinalStatAnalysisV2() {}
	void			run() {}
	cAnalysis*		getCorrAnalysis() { return NULL; }
	cPPPAnalysisV2*	getPPPAnalysis() { return NULL; }
	static int		calcFinalStat(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
	{
		return 0;
	}
};

#endif