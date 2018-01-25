#pragma once
#include "utils/data.h"
#include "global/analysis.h"
#include "utils/vector.h"

#define EMAI_EST_SIZE	8

namespace ONETOOL {

extern const char	*Sa_EmaiEst_headers[EMAI_EST_SIZE];

typedef wsReal**	(*invFunc)(wsReal**, wsUint);

invFunc getInvMethod();
wsReal diagSum3matrix(wsReal **Ra_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Ra_mat2, wsUint N_mat2Row, wsUint N_mat2Col,
	wsReal **Ra_mat3, wsUint N_mat3Row, wsUint N_mat3Col);
wsReal diagSum3matrixSSE(wsReal **Ra_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Ra_mat2, wsUint N_mat2Row, wsUint N_mat2Col,
	wsReal **Ra_mat3, wsUint N_mat3Row, wsUint N_mat3Col, char B_mat3tp=0);
wsReal** GJ_invMatrix(wsReal **Ra_mat, wsUint N_row);
wsReal** subsetMatrix(wsReal **Rp_inpMatrix,
	wsUint N_rowStart, wsUint N_rowCount,
	wsUint N_colStart, wsUint N_colCount);

class cEmAiAnalysisV2 : public cAnalysis
{
	char		B_quiet;
	char*		Ba_isFailed;
	wsUint		N_anaSamp;
	cStdMatrix	*Mp_EVX;
	cDiagMatrix	*Mp_EV;
	char		B_dealloc;
	wsReal		R_trPhiInv;
	cMatrix		*Mp_phi;
	cSymMatrix	M_XtX, M_XtXi;
	cStdMatrix	M_XtP, *Mp_Xt;
	cDiagMatrix	M_dInv;
	cStdMatrix	M_Yt;
	cStdMatrix	M_XtY_t, M_YtP;
	cVector		V_YtY;
	wsReal*		Ra_sampwgt;

	char		*Ba_excl;		///< [#samp] Indicates that ith sample in final
	char		B_doSpecDcmp;
	char		B_doEM;			///< [1] Whether applying EM algorithm(1) or not(0)
	// dataset included in this model(1) or not(0)
//	wsUint		*Na_sampIdx;
	wsReal		**Ra_estimates;
	wsReal		**Ra_calcBLUP;
	wsReal		**Ra_calcPred;
	wsReal		**Ra_alpha;
	void EM_ML(wsReal *Rp_sig2, wsReal *Rp_sig2g, cDiagMatrix &M_d,
		cDiagMatrix &M_dInv, cStdMatrix &M_XtP, cVector &V_YtP,
		cSymMatrix &M_XtX, cVector &V_XtY, cStdMatrix &M_Xt, wsReal R_YtY,
		cVector &V_Y, wsUint N_idxPheno);
	void NR_ML(wsReal *Rp_sig2, wsReal *Rp_sig2g, cDiagMatrix &M_d,
		cDiagMatrix &M_dInv, cStdMatrix &M_XtP, cVector &V_YtP,
		cSymMatrix &M_XtX, cVector &V_XtY, cStdMatrix &M_Xt, wsReal R_YtY,
		cVector &V_Y, wsUint N_idxPheno);
	void EM_REML(wsReal *Rp_sig2, wsReal *Rp_sig2g, cDiagMatrix &M_d,
		cDiagMatrix &M_dInv, cStdMatrix &M_XtP, cVector &V_YtP,
		cSymMatrix &M_XtX, cVector &V_XtY, cStdMatrix &M_Xt, wsReal R_YtY,
		cVector &V_Y, wsUint N_idxPheno);
	void AI_REML(wsReal *Rp_curTheta, wsReal *Rp_curK, cDiagMatrix &M_d,
		cDiagMatrix &M_dInv, cStdMatrix &M_XtP, cVector &V_YtP,
		cSymMatrix &M_XtX, cVector &V_XtY, cStdMatrix &M_Xt, wsReal R_YtY,
		cVector &V_Y, wsUint N_idxPheno);
	void AI_REMLv2(wsReal *Rp_curTheta, wsReal *Rp_curK, cDiagMatrix &M_d,
		cDiagMatrix &M_dInv, cStdMatrix &M_XtP, cVector &V_YtP,
		cSymMatrix &M_XtX, cVector &V_XtY, cStdMatrix &M_Xt, wsReal R_YtY,
		cVector &V_Y, wsUint N_idxPheno);
	void spectralDecomposition(cVector &V_sig2, cVector &V_sig2g);
	void _init(cIO *Cp_IO, cStdMatrix &M_inpXt, cMatrix &M_phi,
		char B_isREML, char *Ba_iExcl,
		cStdMatrix *Mp_eVec, cDiagMatrix *Mp_eVal);
	void _exportBLUP();
public:
	cEmAiAnalysisV2(cIO *Cp_IO, cStdMatrix &M_inpXt, cVector &V_inpY,
		cMatrix &M_phi, char B_isREML, char *Ba_iExcl,
		cStdMatrix *Mp_eVec=NULL, cDiagMatrix *Mp_eVal=NULL);
	cEmAiAnalysisV2(cIO *Cp_IO, cStdMatrix &M_inpXt, cStdMatrix &M_inpYt,
		cMatrix &M_phi, char B_isREML, char *Ba_iExcl,
		cStdMatrix *Mp_eVec=NULL, cDiagMatrix *Mp_eVal=NULL);
	~cEmAiAnalysisV2();
	void		quiet() { B_quiet = 1; }
	char		getIsFailed() { for (wsUint i=0 ; i<Cp_IO->sizePheno() ; i++)
		if (Ba_isFailed[i]) return 1; return 0; }
	void		run();
	void		doSpecDcmp() { B_doSpecDcmp = 1; }
	wsReal**	getBeta() { return Ra_alpha; }
	int			getEst();
	wsRealCst		getML(wsUint i=0) { return Ra_estimates[i][2]; }
	wsRealCst*		getBLUP(wsUint i=0) { return Ra_calcBLUP[i]; }
	wsRealCst*		getPred(wsUint i=0) { return Ra_calcPred[i]; }
	wsRealCst*		getEstimates(wsUint i) { return Ra_estimates[i]; }
	wsReal**	getEstimates() { return Ra_estimates; }
};

} // End namespace ONETOOL
