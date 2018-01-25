#pragma once
#include "global/analysis.h"
#include "utils/matrix.h"
#include "utils/vector.h"
#include "utils/integrate.h"

#define USE_CONT_SKATO_V2

namespace ONETOOL {

int forGene_equal(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);

class cPDDTAnalysis;
class cEmAiAnalysisV2;
class cMatrix;
class cSetManagerAnalysis;

typedef map<string,vInt>	mGeneDef;
typedef mGeneDef::iterator	mGeneDef_it;

/* #tsamp	Number of samples included in gene-set analysis, samples having
 *			missing phenotype is excluded in this number
 * #samp	Number of samples included in the dataset analyzed
 * #SNP		Number of markers included in the dataset analyzed
 * 
 */

typedef struct _xLiuCov
{
	wsReal	R_muX, R_sigX, R_muQ, R_sigQ;
	wsReal	R_l, R_delta, R_c2, R_c3, R_c4;
	wsReal	R_stat, R_pVal;
} xLiuCov;

typedef struct _xGeneTestParam {
	// Common
	cIO*		Cp_IO;
	vStr&		Xv_smpbuf;
	wsReal		R_C;
	cVector&	V_hat;
	wsRealCst	R_vars;
	char		B_isCont;
	char**		Na_permPhe;
	cVector&	V_aV;
	cVector*	Vp_yGS;
	cVector*	Vp_wGS;
	char*		Ba_misPheno;
	// Output
	cExporter*	Cp_exporter;
	cExporter*	Cp_naTest;
	// GG - GGEMMA
	cStdMatrix&	M_XtGG;
	wsRealCst	R_L0GG; // LogL_0
	// SK - SKAT
	cSymMatrix&	M_HH_SK;
	wsRealCst	R_pppMult;
	// PG - Pedgene
	cVector&	V_resPG;
	wsRealCst	R_cFtrPG;
	// PC - PEDCMC
	cVector&	V_availPhenoMask;
	// FX - FARVATx
	cVector&	V_phiInv1_FX;
	wsRealCst	R_1phiInv1_inv_FX;
	vInt&		Nv_idxFemale;
	vInt&		Nv_idxFnder;
	wsRealCst	R_1phiInv1_inv_Xall;
	cStdMatrix&	M_I_pi1pi1i;
	// FM - MFARVAT
	cVector&	V_1TPA;
	wsRealCst	R_sumTPAPT;
	cSymMatrix&	M_TPAPT;
	cStdMatrix&	M_Tt;
	cStdMatrix&	M_TPA;
	// Fm - metaFARVAT
	cVector&	V_metaFarvat;
	// KB - KBAC
	wsRealCst	R_alphaKB;
	// VT - wVT
	wsRealCst	R_pheVTmean;
	wsMat		Ra_pheVT;
	// SO - SKAT-o
	wsVec		Ra_ySO;
	wsUintCst	N_rho;
	wsVec		Ra_rho;
	cMatrix*	Mp_P;
	cMatrix*	Mp_Psq;
	cMatrix*	Mp_phi;
	cMatrix*	Mp_phiInv;
#ifdef USE_CONT_SKATO_V2
	cStdMatrix	M_r2, M_r1r, M_1r2;
	cStdMatrix	M_r, M_1r;
#endif
} xGeneTestParam;

typedef struct _xIntgSKATO
{
	wsReal	*Ra_qLiu;		// qmins	[nRho]
	wsReal	*Ra_rRho;		// rRho		[nRho]
	wsVecCst	Ra_rho;		// rhos		[nRho]
	wsReal	R_muQ;			// muQ
	wsReal	R_sigQ;			// sigQ
	wsReal	R_sigXi;		// sigXi
	xLiuCov	*Xp_liuFinal;	// lius		[nRho]
	wsUint	N_rho;
	wsReal	R_absTol;
} xIntgSKATO;

xLiuCov*	liuCov(wsSym Ra_cov, wsUint N_sz, wsSym Ra_R,
	wsReal **Ra_ttt=NULL, wsReal **Ra_ttt2=NULL,
	wsReal R_muQ=WISARD_NAN);
wsVec		getOptimalRhos(wsUint& N_rho);
wsReal		qLiu(wsReal R_pVal, xLiuCov *Xp_liuCov);
void		integrateOptimal(double *Ra_X, int N, xIntgInp *Vp_i);

/* !!! FARVAT-specific analysis !!! */
#if ((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE & FUNCTION_FARVAT) || \
	(TOOLSET_TYPE & FUNCTION_MFARVAT) || (TOOLSET_TYPE & FUNCTION_FARVATX)

class cGSTestAnalysis : public cAnalysis
{
	cSetManagerAnalysis
				*Cp_anaGsm;
	cPPPAnalysisV2
				*Cp_anaPPP;
	cEmAiAnalysisV2
				*Cp_anaEmai;
	cAnalysis	*Cp_anaCorr;			///< A pointer to class instance having phi matrix

	wsUint		N_szPhi;
	wsSym       Ra_phiGS;			///< Matrix[#pheno,#pheno], [--longitudinal] Variance-covariance matrix of phenotypes
									///< Matrix[#tsamp,#tsamp], [otherwise]      Sample relatedness matrix
	cSymMatrix	M_A;				///< Matrix[#tsamp,#tsamp], Projection of phiInv, only available for multi FARVAT
	wsReal		**Ra_origPhi;
	wsReal		R_1phiInv1_inv;		///< (1^T %*% phi %*% 1)^-1
	vInt		Nv_idxFemale;			///< [#tsamp] Indices of male samples
	vInt		Nv_idxFnder;
	wsReal		R_1phiInv1_inv_X;	///< (1^T %*% phi %*% 1)^-1 for --farvatx
	wsReal		R_1phiInv1_inv_Xall;	///< (1^T %*% phi %*% 1)^-1 for --farvatx
	cStdMatrix	M_I_pi1pi1i;		///< Matrix[#tsamp,#tsamp], I - akinshipInv*IIphiInvII_Inv, for --farvatx
	cMatrix*    Mp_phi;
//	wsReal		**Ra_phiGS;			///< Matrix[#tsamp,#tsamp], resized phi matrix according to Ba_misPheno
//	wsReal		*Ra_hat;			///< Array[#tsamp], Equivalent to h
	cVector		V_hat;				///< Array[#tsamp*#pheno]
	char*		Ba_misPheno;		///< Array[#samp], 1 when the phenotype of ith sample is missing, 0 otherwise
	wsUint		N_anaSamp;			///< Number of samples included in gene-set analysis
	wsUint		N_availPhenoSamp;
	cMask		V_availPhenoMask;	///< Array[#samp], 1 (should be included), 0 (not to be included)
//	wsReal		*Ra_yGS;			///< Array[#tsamp], resized phenotype array according to Ba_misPheno
	cVector		*Vp_yGS;			///< Array[#tsamp*#pheno], resized phenotype array according to Ba_misPheno
	cVector		*Vp_yGSunAdj;
//	wsReal		*Ra_yUnadjGS;		///< Array[#tsamp], resized phenotype array according to Ba_misPheno
//	wsReal		**Ra_aV;
	cVector		V_aV;				///< Array[#tsamp*#pheno]
//	wsReal		*Ra_wGS;			///< Array[#SNP], weights of SNPs determined by founder's MAF
	cVector*	Vp_wGS;			
	cVector		V_phiInv1;			///< Array[#tsamp], rowSum of Mp_phi
	cVector		V_phiInv1_X;		///< Array[#tsamp], rowSum of Mp_phi for --farvatx
	cVector		V_phiInv1_Xall;		///< Array[#tsamp], rowSum of Mp_phi for --farvatx
	wsReal		R_1phiInv1;

	cMatrix*	Mp_phiInv;

	/*** MULTI-FARVAT ***/
	cVector		V_1TPA;
	wsReal		R_sumTPAPT;
	cSymMatrix	M_TPAPT;
	cStdMatrix	M_Tt;				///< Array[#pheno*#tsamp], adjusted phenotype for --mfhom
	cStdMatrix	M_TPA;				///< Array[#pheno*#tsamp], adjusted phenotype for --mfhet

	/* Note that below two parameters will be determined at run() */
	wsReal		R_pppMult, R_vars;

	/* VT */
	wsMat		Ra_pheVT;
	wsReal		R_pheVTmean;

	/* SKAT-o */
	cMatrix*	Mp_P;
	cMatrix*	Mp_Psq;
	wsReal		R_C;

	/* Meta-FARVAT */
	cVector		V_metaFarvat;

	/* Logistic indep. */
	char		B_isCont;
	wsReal		*Ra_logisticV;
	wsReal		R_logisticS;
	wsReal		*Ra_phiHat;

	wsReal*		_makePheno(char &B_isContinuous);
//	wsReal*		_makeWeight();
	wsMat		_makeCorr(wsUint *Np_szPhi);
public:
	cGSTestAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP,
		cAnalysis *Cp_inpAnaPhi,
		cSetManagerAnalysis *Cp_inpGsmAna=NULL,
		cEmAiAnalysisV2 *Cp_inpEmaiAna=NULL);
	~cGSTestAnalysis();

	wsUintCst		_filterByGenotypingRate(const char *S_gs, vInt &X_set,
					char *Ba_filt, wsReal *Ra_mask);//, wsUint *Np_case);

	wsUintCst	getAvailableSampleSize() { return N_anaSamp; }
	mGeneDef&	getGeneDef();

	void		run();
};

#else

/* DUMMY DEFINITION */
class cGSTestAnalysis : public cAnalysis
{
public:
	cGSTestAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP,
		cAnalysis *Cp_inpAnaPhi,
		cSetManagerAnalysis *Cp_inpGsmAna=NULL,
		cEmAiAnalysisV2 *Cp_inpEmaiAna=NULL) : cAnalysis(Cp_inpIO) {}
	~cGSTestAnalysis() {}

	wsUintCst		_filterByGenotypingRate(const char *S_gs, vInt &X_set,
		char *Ba_filt, wsReal *Ra_mask) { return 0; }

	wsUintCst	getAvailableSampleSize() { return 0; }
	mGeneDef&	getGeneDef() { mGeneDef *r = new mGeneDef; return *r; }

	void		run() {}
};

#endif

} // End namespace ONETOOL
