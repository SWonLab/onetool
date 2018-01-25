#ifndef __WISARD_GSCA_H__
#define __WISARD_GSCA_H__
#pragma once
#include "global/analysis.h"
#include "setmgr.h"

namespace ONETOOL {

typedef struct _xGeSCAres {
	_xGeSCAres(wsUint Nw, wsUint Nc, wsUint Nb,
 		wsUint No, wsUint Nl) {
// 		wsAlloc(vec_FIT, wsFloat, N_perm);
// 		wsAlloc(vec_AFIT, wsFloat, N_perm);
// 
// 		wsAlloc(MatW, wsFvec, N_perm);
// 		wsAlloc(MatB, wsFvec, N_perm);
// 		wsAlloc(MatC, wsFvec, N_perm);
// 		wsAlloc(MatR2M, wsFvec, N_perm);
// 		wsAlloc(MatR2S, wsFvec, N_perm);
// 		LOOP (i, N_perm) {
// 			wsAlloc(MatW[i], wsFloat, Nw);
// 			wsAlloc(MatB[i], wsFloat, Nb);
// 			wsAlloc(MatC[i], wsFloat, Nc);
// 			wsAlloc(MatR2M[i], wsFloat, No);
// 			wsAlloc(MatR2S[i], wsFloat, No+Nl);
// 		}
	};
	_xGeSCAres() {
	};
	~_xGeSCAres() {
		FOREACH(vector<wsFvec>::iterator, MatW, i) DEALLOC(*i);
		FOREACH(vector<wsFvec>::iterator, MatB, i) DEALLOC(*i);
		FOREACH(vector<wsFvec>::iterator, MatC, i) DEALLOC(*i);
		FOREACH(vector<wsFvec>::iterator, MatR2M, i) DEALLOC(*i);
		FOREACH(vector<wsFvec>::iterator, MatR2S, i) DEALLOC(*i);
	};
	vector<wsFloat>	vec_FIT;
	vector<wsFloat>	vec_AFIT;
	vector<wsFvec>	MatW, MatB, MatC;
	vector<wsFvec>	MatR2M, MatR2S;
} xGeSCAres;

typedef struct _xCoord {
	_xCoord(USHORT_t _x, USHORT_t _y) {
		x = _x;
		y = _y;
	}
	_xCoord(wsUint _x, wsUint _y) {
		x = (USHORT_t)_x;
		y = (USHORT_t)_y;
	}
	USHORT_t x, y;
} xCoord;

typedef vector<xCoord>		vCoord;
typedef vCoord::iterator	vCoord_it;

typedef vector<vReal>		vLambda;
typedef vLambda::iterator	vLambda_it;

class cPermuter;

typedef struct _xGescaThread {
	/* For gesca */
	char			B_formative;
	cStdMatrix*		Mp_W0_t;
	cStdMatrix*		Mp_B0t;
	cStdMatrix*		Mp_C0t;
	vInt&			bindex, cindex;
	xGeSCAres		X_res;
	/* For fixed permutation sequence */
	cPermuter*		Cp_perm;
} xGescaThread;

typedef struct _xGscaAnaThread {
	char			B_overlap;
	vInt&			Nv_dist;
	wsUint			N_mnfs;
	cStdMatrix*		Mp_Xt;
	cStdMatrix*		Mp_W0;
	cStdMatrix&		M_origYt;
	mDataIdx&		Xm_mnfs2idx;
	vStr&			Sv_idx2mnfs;
	vStr&			Sv_latent;
	mvStr&			Sv_mnfsInLatent;
	vvInt&			windex;
	/* Shared between doubleRidge & optLambda */
	vLambda&		Rv_lambdas;
	wsUint*			Na_samp;
	/* For optLambda */
	wsUint			N_cv;
	wsUint**		Na_indices;
	wsVec			Ra_devs;
	wsUint**		Na_iters;
	/* For doubleRidge */
	wsVec			Ra_trueA;
	wsVec			Ra_trueW;
	vector<wsVec>&	Rv_permA;
	vector<wsVec>&	Rv_permW;

	/* For fixed permutation sequence */
//	vector<wsUint*>&
//				Nv_seqPerm;
//	wsMat		Ra_valPerm;
	cPermuter*		Cp_perm;
	vvInt&			M_phenoRel;

	/* For GeSCA */
	cStdMatrix*		Mp_B0t;
	cStdMatrix*		Mp_C0t;
	vvInt&			cindex;
	vInt&			bindex;
	xGeSCAres&		X_res;
	wsUint			N_obs;
	wsUint			N_ltt;
} xGscaAnaThread;

typedef struct _xGescaCvThread {
	char		B_formative;
	cStdMatrix&	M_Yt;
	cStdMatrix*	Mp_Xt;
	cStdMatrix&	M_W0;
	cStdMatrix& M_B0_t;
	cStdMatrix& M_C0_t;
	vLambda		Rv_lambdas;
	vvInt&		windex;
	vInt&		bindex;
	vvInt&		cindex;
	wsUintCst	N_obs;
	wsUintCst	N_ltt;
	wsUintCst	N_cv;
	vStr&		Sv_idx2mnfs;
	vStr&		Sv_latent;
	mvStr&		Sv_mnfsInLatent;
	wsUint*		Na_samp;
	wsUint**	Na_idxCV;
	wsUint**	Na_iters;
	wsVec		Ra_devs;
} xGescaCvThread;

/* !!! GSCA-specific analysis !!! */
#if ((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE == TOOLSET_PHARAOH)

class cPharaohAnalysis : public cAnalysis
{
	cSetManagerAnalysis*	Cp_anaSM;
	cPPPAnalysisV2*			Cp_anaPPP;
	mGeneSet				Xm_gs;
	mGeneDef				Xm_gd;
	char					_initGeneDef(mGeneDef& Xm_oriGD, string j, xMaf* Xp_maf, wsUintCst NFILTV);
	void					_makeGeneSummary(mGeneDef& Xm_gd, const char* Ba_miss,
		map<string, wsVec>& Xm_collapsed, mDataIdx& Xm_nvrtGene, wsUint N_anaSamp);
	void					_makePathwaySummary(char B_noOverlap,
		vector<wsVec>& Xv_X, map<string, wsVec>& Xm_collapsed,
		mDataIdx& Xm_idxGene, mDataIdx& Xm_nvrtGene,
		vStr& Sv_gnames, vStr& Sv_latent, mvStr& Sv_mnfsInLatent,
		vInt& Nv_nvrtInLatent);
	void					run_GESCA(cPermuter* Cp_perm, char B_formative,
		wsUint N_perm, cStdMatrix& M_Yt, cStdMatrix& M_Xt,
		cStdMatrix& M_Wr, vStr& Sv_idx2ltt, vvInt& windex,
		vInt& Nv_pathCoords, wsRealCst R_lambdaG, wsRealCst R_lambdaP,
		vStr& Sv_idx2mnfs);
	void					run2(wsUintCst N_pheno, wsUint N_anaSamp, char* Ba_miss, char B_formative, vInt Nv_comb);
public:
	cPharaohAnalysis(cIO* Cp_inpIO, cSetManagerAnalysis* Cp_inpAnaSM,
		cPPPAnalysisV2* Cp_inpAnaPPP);
	~cPharaohAnalysis();

	void run();
};

#else

/* DUMMY DEFINITION */
class cPharaohAnalysis : public cAnalysis
{
public:
	cPharaohAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis* Cp_inpAnaSM,
		cPPPAnalysisV2* Cp_inpAnaPPP) : cAnalysis(Cp_inpIO) {}
	~cPharaohAnalysis(){}

	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
