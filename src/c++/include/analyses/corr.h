#pragma once
#include "global/analysis.h"
#include "global/io.h"

namespace ONETOOL {

struct xCorrThread
{
	wsUint		N_variant;
	wsUint		N_sample;
	wsFloat**	Ra_data;
	char		**Na_data;
	wsReal		**Ra_corr;
	int			**Na_corrN;
	wsFloatCst	*Ra_corrMask;
	char		*Na_corrMask;
	WISARD_pc	**Na_gt;	/* For Kendall's tau */
	WISARD_pc	**Na_neq;	/* For Kendall's tau */
	__int64		N_szComb;	/* For Kendall's tau */
	wsUint		N_szBlk;	/* For Kendall's tau */
	wsUint		N_szSNP;	/* For Kendall's tau */
	wsUint		**Na_same;
	wsUint		N_last;
	int			B_isComplete;
	vVariant*	Xp_SNP;

	/* For countN */
	wsUint		N_mBlk;
	WISARD_pc**	Na_mMissing;
	wsUint*		Na_cachedN;
	wsUint		N_noIncCor;
};

wsReal corPearsonComplete(wsReal *Na_v1, wsReal *Na_v2, wsUint N_SNP,
   wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsRealCst *Ra_cm=NULL);
wsReal corPearsonIncomplete(wsReal *Na_v1, wsReal *Na_v2, wsUint N_SNP,
   wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsRealCst *Ra_cm=NULL);

wsReal _corMedianComplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corMeanComplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corPearsonComplete(wsFloat *Na_v1, wsFloat *Na_v2, wsUint N_SNP,
   wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corPearsonIncomplete(wsFloat *Na_v1, wsFloat *Na_v2, wsUint N_SNP,
   wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corPearsonIncomplete(char *Na_v1, char *Na_v2, wsUint N_obs,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, char *Na_cm=NULL);
wsReal _corMeanComplete(char *Na_v1, char *Na_v2, wsUint N_SNP,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corMedianIncomplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corMeanIncomplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corMeanIncomplete(char *Na_v1, char *Na_v2, wsUint N_SNP,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corTauIncomplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corTauIncomplete(char *Na_v1, char *Na_v2, wsUint N_SNP,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corTauComplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corTauComplete(char *Na_v1, char *Na_v2, wsUint N_SNP,
	wsUint *Np_inc=NULL, wsReal *Ra_cors=NULL, wsFloatCst *Ra_cm=NULL);
wsReal _corTauComplete(char *Na_v1, char *Na_v2, wsUint N_SNP,
	wsUint *Np_inc, wsReal *Ra_cors, char *Na_cm);

class cMatrix;
class cCorrAnalysis : public cAnalysis
{
	cMatrix		*Cp_corr;
//	wsReal		**Ra_corr;	///< Matrix[size #samp*#samp], Correlation between sample (i) and (j) */
	int			**Na_corrN;	///< Matrix[size #samp*#samp], # of included samples(i.e. non-missing samples) of correlation between sample (i) and (j)
	wsReal		*Ra_ysum;	///< Array[size #SNP], Values sum(x[1:#SNP])
	wsReal		*Ra_y2sum;	///< Array[size #SNP], Values sum(x[1:#SNP]^2)
	wsUint		N_corSamp;
	cPPPAnalysisV2
				*Cp_anaPPP;

	void		_export();
	void		_pairedAltCorLoad(cStream &C_altcor, char *Sp_buf,
		wsMat Ra_corr, wsStrCst Sp_corPath);
	void		_emmaxAltCorLoad(cStream &C_altcor, char *Sp_buf,
		wsMat Ra_corr, wsStrCst Sp_corPath);
	template <typename T>
	void		_getGRM(wsUint N_corSamp, wsUint N_vrt, wsVecCst Ra_ppp, T** Na_data, wsMat Ra_corr);
public:
	cCorrAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP);
	cCorrAnalysis(cIO *Cp_inpIO, const char *Sp_corPath);
	~cCorrAnalysis();

	void		run();
	cMatrix&	getCorrClass() { return *Cp_corr; }
	wsReal**	getCorr();
	static int	calcCorr(int N_idx, void *Vp_shareData, void *Vp_data,
		void *Vp_result);
	static int	calcCorrV2(int N_idx, void *Vp_shareData, void *Vp_data,
		void *Vp_result);
	static int	calcCorrV2pearson(int N_idx, void *Vp_shareData, void *Vp_data,
		void *Vp_result);
	static int	calcCorrV2tauS1(int N_idx, void *Vp_shareData, void *Vp_data,
		void *Vp_result);
	static int	calcCorrV2tauS2(int N_idx, void *Vp_shareData, void *Vp_data,
		void *Vp_result);
	static int	calcCorrV2tau(int N_idx, void *Vp_shareData, void *Vp_data,
		void *Vp_result);
	static int	calcCorrV2_gpu(int N_idx, void *Vp_shareData,
		void *Vp_data, void *Vp_result);
	static int	calcBN1(int N_idx, void *Vp_shareData, void *Vp_data,
		void *Vp_result);
	static int	calcIBS(int N_idx, void *Vp_shareData, void *Vp_data,
		void *Vp_result);
	static int calcHamming(int N_idx, void *Vp_shareData, void *Vp_data,
		void *Vp_result);
	static int	calcTau(int N_idx, void *Vp_shareData, void *Vp_data,
		void *Vp_result);
};

int forAllSample_corr(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);

} // End namespace ONETOOL
