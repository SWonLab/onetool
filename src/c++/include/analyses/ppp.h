#ifndef __WISARD_PPP_H__
#define __WISARD_PPP_H__
#pragma once
#include "global/analysis.h"
#define USE_CORR00

namespace ONETOOL {

struct xNucFamily;
class cFamStrAnalysis;
class cPPPAnalysisV2 : public cAnalysis
{
	wsUint	N_noIncCor, N_sexMarker;
	wsFloat*	Ra_corrMask;
	char*	Na_corrMask;
	wsReal	*Ra_meanX;				///< Array[size #SNP], Mean of founders' genotype
	char	B_doNorm;				///< Do normalization of genotype if 1 based on ppp
	mSampPtr
			Xm_misFounder;
	strMap
			Xm_misFounderDtIdx;
	wsReal**
			Ra_misFounderGeno;
	cFamStrAnalysis
			*Cp_anaFamStr;

	void	_doImpute(xSample *Xp_misFounder, wsUint N_dtIdxMisFounderGeno); /* Dosage-safe */
	void	_doCalcPPP(); /* Dosage-safe */
	void	_calcPPPcontribution(wsUint N_sIdx, xNucFamily &X_nf,
		wsReal *Rp_misDen, wsReal *Rp_misDiv); /* Dosage-safe */
public:
	wsReal	*Ra_ppp;				///< Array[size #SNP], PPP statistics for all SNPs
	wsReal	*Ra_ppp_sqrt2pq;		///< Array[size #SNP], Values sqrt(2*ppp*(1-ppp))
	/* ONLY FOR THREAD */
	wsUint	*Na_ind;
	wsReal	*Ra_sum;

	cPPPAnalysisV2(cIO *Cp_inpIO,
		cFamStrAnalysis *Cp_inpAnaFamStr, char B_inpDoNorm=1);
	~cPPPAnalysisV2();

	void		normalize(wsFloat **Ra_data);
	void		run();
	wsRealCst*		getFounderMean() { return Ra_meanX; }
	wsRealCst*		getPPPsq() { return Ra_ppp_sqrt2pq; }
	wsRealCst*		getPPP() {
		getIO()->getMAF();
		return Ra_ppp;
	}
	wsReal**	getMisFounderGeno() { return Ra_misFounderGeno; }
	wsUint		getCorMask(wsFloat** Rp_corMask, char **Np_corMask,
		wsUint *Np_nSexMarker) {
		*Rp_corMask = Ra_corrMask;
		*Np_corMask = Na_corrMask;
		if (Np_nSexMarker) *Np_nSexMarker = N_sexMarker;
		return N_noIncCor;
	}
};

typedef struct _xNormThr
{
	cPPPAnalysisV2*	Cp_anaPPP;
	wsFloat**		Ra_data;
	char**			Na_data;
	wsFloat**		Ra_dsg;
} xNormThr;

} // End namespace ONETOOL

#endif
