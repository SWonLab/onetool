#ifndef __WISARD_PCA_H__
#define __WISARD_PCA_H__
#pragma once

#include "global/analysis.h"
#include "fam.h"

namespace ONETOOL {

class cPPPAnalysis;
class cPPPAnalysisV2;
class cFamStrAnalysis;

typedef map<string,xNucFamily>	mNucFam;
typedef mNucFam::iterator		mNucFam_it;

/* #SNP		Number of markers included in the dataset analyzed
 * #samp	Number of samples included in the dataset analyzed
 * #founder	Number of founders included in the dataset analyzed
 */

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cPcaAnalysisV2 : public cAnalysis
{
	wsReal		**Ra_cmds;
	wsReal		**Ra_corr;			///< Matrix[#samp*#samp], Correlation between sample (i) and (j)
	int			**Na_corrN;			///< Matrix[#samp*#samp], # of included samples(i.e. non-missing samples) of correlation between sample (i) and (j)
	wsReal		*Ra_ysum;			///< Array[#SNP], Values sum(x[1:#SNP])
	wsReal		*Ra_y2sum;			///< Array[#SNP], Values sum(x[1:#SNP]^2)
public:
	strMap		X_iid2fouIdx;		///< Map[#samp], Convert IID into index on founder array
	wsFloat		**Ra_pcaInputData;	///< Matrix[#samp*#SNP], Genotypes of (i)th founder, and its sequence is same with in the dataset
	char		B_isPcaDataComplete;
	wsReal		*Ra_expectedX;		///< Array[#founder], E(Xi) i=1:#SNP 
	wsReal		**Ra_cov;			///< Matrix[#founder*#founder], Variance-covariance matrix of Ra_corr
	wsUint		N_pcaSample;		///< Number of samples included in PC analysis
	cFamStrAnalysis
				*Cp_anaFamStr;
	cPPPAnalysisV2
				*Cp_anaPPP;
	vStr		Sa_fidFounder, Sa_iidFounder;
	vBool		Xa_isBothParentsMissed;

	cPcaAnalysisV2(cIO *Cp_IO, cPPPAnalysisV2 *Cp_inpAnaPPP, cFamStrAnalysis *Cp_inpAnaFamStr);
	~cPcaAnalysisV2();
	void		run();
	/**
	 * cPcaAnalysisV2::calcEx Calculate E(X) from correlation matrix from founder dataset
	 *
	 * @param     N_idx			An index of thread, starts from 0
	 * @param     Vp_shareData	A pointer to referenced data, which is shared for all threads
	 * @param     Vp_data		A pointer to referenced data, which is thread-specific
	 * @param     Vp_result		An array of result, NULL if it is not necessary
	 * @return    (int)			0 if normally executed
	 */
	static int	calcEx(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result);
	/**
	 * cPcaAnalysisV2::calcXeX Normalizes correlation matrix from founder dataset with the formula X-E(X)
	 *
	 * @param     N_idx			An index of thread, starts from 0
	 * @param     Vp_shareData	A pointer to referenced data, which is shared for all threads
	 * @param     Vp_data		A pointer to referenced data, which is thread-specific
	 * @param     Vp_result		An array of result, NULL if it is not necessary
	 * @return    (int)			0 if normally executed
	 */
	static int	calcXeX(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result);
	/**
	 * cPcaAnalysisV2::calcCov Calculate variance-covariance matrix of correlation matrix from founder dataset
	 *
	 * @param     N_idx			An index of thread, starts from 0
	 * @param     Vp_shareData	A pointer to referenced data, which is shared for all threads
	 * @param     Vp_data		A pointer to referenced data, which is thread-specific
	 * @param     Vp_result		An array of result, NULL if it is not necessary
	 * @return    (int)			0 if normally executed
	 */
	static int	calcCov(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result);

	void		_stepPolyEffectMat(wsReal **Rp_poly);
	int			_stepFullPCmatrix_MF(vStr_it &iit, mNucFam_it &it,
		wsMat Ra_pc, wsMat Ra_mds, wsUint *Np_fillCount, wsUint N_srcIdx,
		wsUint N_origSample, wsUint N_tempIdxForMissingSample, wsUint N_pc);
	/**
	 * cPcaAnalysisV2::_stepFullPCmatrix Calculates PC scores for all samples from
	 *   the samples that included in PCA
	 *
	 * @param     N_pc			Number of PCs to be calculated to all samples
	 * @return    (MAT_t)		Computed full PC matrix (row=sample, col=#pc)
	 */
	wsMat		_stepFullPCmatrix(wsUint N_pc);
	/**
	 * cPcaAnalysisV2::getPcaSampleSize Returns the number of samples that
	 *   included in PCA
	 *
	 * @return    (uint)		the number of samples included in PCA
	 */
	wsUint		getPcaSampleSize();
};

#else

/* DUMMY DEFINITION */
class cPcaAnalysisV2 : public cAnalysis
{
public:
	strMap		X_iid2fouIdx;		///< Map[#samp], Convert IID into index on founder array
	wsReal		**Ra_pcaInputData;	///< Matrix[#samp*#SNP], Genotypes of (i)th founder, and its sequence is same with in the dataset
	char		B_isPcaDataComplete;
	wsReal		*Ra_expectedX;		///< Array[#founder], E(Xi) i=1:#SNP 
	wsReal		**Ra_cov;			///< Matrix[#founder*#founder], Variance-covariance matrix of Ra_corr
	wsUint		N_pcaSample;		///< Number of samples included in PC analysis
	cFamStrAnalysis
				*Cp_anaFamStr;
	cPPPAnalysisV2
				*Cp_anaPPP;
	vStr		Sa_fidFounder, Sa_iidFounder;
	vBool		Xa_isBothParentsMissed;

	cPcaAnalysisV2(cIO *Cp_IO, cPPPAnalysisV2 *Cp_inpAnaPPP, cFamStrAnalysis *Cp_inpAnaFamStr)
	 : cAnalysis(Cp_IO) {}
	~cPcaAnalysisV2() {}
	void		run() {}
	static int	calcEx(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result) { return 0; }
	static int	calcXeX(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result) { return 0; }
	static int	calcCov(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result) { return 0; }

	void		_stepPolyEffectMat(wsReal **Rp_poly) {}
	int			_stepFullPCmatrix_MF(vStr_it &iit, mNucFam_it &it,
		wsReal **Ra_pc, wsUint *Np_fillCount, wsUint N_srcIdx,
		wsUint N_origSample, wsUint N_tempIdxForMissingSample, wsUint N_pc) { return 0; }
	wsMat		_stepFullPCmatrix(wsUint N_pc) { return NULL; }
	wsUint		getPcaSampleSize() { return 0; }
};

#endif

} // End namespace ONETOOL

#endif
