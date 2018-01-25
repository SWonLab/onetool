#pragma once
#include "global/analysis.h"
#include "utils/vector.h"
#include "utils/matrix.h"

namespace ONETOOL {

typedef struct _xFemmaRes {
	wsReal R_lambda;
	wsReal R_logL;
	wsReal R_Twald, R_Pwald;
	wsReal R_Tscore, R_Pscore;
	wsReal R_Tlrt, R_Plrt;
	wsUint B_failed;
	wsReal	R_time;
} xFemmaRes;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cStdMatrix;
class cFemmaAnalysis : public cAnalysis
{
	cStdMatrix	UtW;
	cVector		Uty;
	cStdMatrix	Uab_t;
	double		l_mle_null;

	const char	*Ba_isInc;
	cAnalysis	*Cp_anaPhi;
	wsUint		N_anaRow;	///< Number of rows of design matrix in alternative test [Value : #cov+2]
	wsReal		**Ra_phi;
	wsReal		**Ra_Pt;
	wsReal		*Ra_yGS;
//	wsReal		R_LdetWtW; /* |W'W| */
	wsUint		N_anaSamp;
//	MAT_t		*Ra_fDMs;
	wsMat		Ra_DM; ///< Design matrix [#cov*#tsamp]

	cStdMatrix	*Mp_Wp_t, *Mp_Wp;
	cVector		*Vp_yp_t;
	cVector		V_y2p_t;
	cVector		V_D;
	xFemmaRes	X_0;
public:
	cFemmaAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaPhi);
	~cFemmaAnalysis();

	xFemmaRes&	get0() { return X_0; }
	wsReal		NR(wsReal R_param, wsUint N_idxSNP, cStdMatrix &M_XtP,
		wsReal *Rp_tau, wsReal ***Rp_ab, cSymMatrix *Mp_WHWinv,
		wsReal *Rp_LdetH, wsReal *Rp_yPy, wsReal *Rp_LdetXHX=NULL);
	void		run();
	static void	_doTest(cStdMatrix &M_Xt, wsUint N_anaSamp,
		wsUint N_cova, wsMat Ra_Pt, cDiagMatrix *Mp_D, wsReal R_l0,
		int N_idxSNP, cVector *Vp_yp_t, cVector &V_y2p_t, xFemmaRes *Xp_ret);
	wsMat		getPt() { return Ra_Pt; }
	cVector&	getD() { return V_D; }
	cVector*	getYp() { return Vp_yp_t; }
	cVector&	getY2p() { return V_y2p_t; }
	cStdMatrix&	getUtW() { return UtW; }
	cVector&	getUty() { return Uty; }
	wsRealCst	getLmleNull() { return l_mle_null; }
	wsUintCst	getAnaSampSize() { return N_anaSamp; }
//	MAT_t		getDM(wsUint N_idx=0) { return Ra_fDMs[N_idx]; }
	wsMat		get0DM() { return Ra_DM; }
	const char*	getIsInc() { return Ba_isInc; }
};

#else

/* DUMMY DEFINITION */

class cStdMatrix;
class cFemmaAnalysis : public cAnalysis
{
	const char	*Ba_isInc;
	cAnalysis	*Cp_anaPhi;
	wsUint		N_anaRow;	///< Number of rows of design matrix in alternative test [Value : #cov+2]
	wsReal		**Ra_phi;
	wsReal		**Ra_Pt;
	wsReal		*Ra_yGS;
	//	wsReal		R_LdetWtW; /* |W'W| */
	wsUint		N_anaSamp;
	//	MAT_t		*Ra_fDMs;
	wsMat		Ra_DM; ///< Design matrix [#cov*#tsamp]

	cStdMatrix	*Mp_Wp_t, *Mp_Wp;
	cVector		*Vp_yp_t;
	cVector		V_y2p_t;
	cDiagMatrix	*Mp_D;
	xFemmaRes	X_0;
public:
	cFemmaAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaPhi) : cAnalysis(Cp_inpIO) {}
	~cFemmaAnalysis(){}

	xFemmaRes&	get0() { return X_0; }
	wsReal		NR(wsReal R_param, wsUint N_idxSNP, cStdMatrix &M_XtP,
		wsReal *Rp_tau, wsReal ***Rp_ab, cSymMatrix *Mp_WHWinv,
		wsReal *Rp_LdetH, wsReal *Rp_yPy, wsReal *Rp_LdetXHX=NULL) {
		return WISARD_NAN;
	}
	void		run() {}
	static void	_doTest(cStdMatrix &M_Xt, wsUint N_anaSamp,
		wsUint N_cova, wsMat Ra_Pt, cDiagMatrix *Mp_D, wsReal R_l0,
		int N_idxSNP, cVector *Vp_yp_t, cVector &V_y2p_t, xFemmaRes *Xp_ret) {}
	wsMat		getPt() { return NULL; }
	cDiagMatrix*
		getD() { return NULL; }
	cVector*	getYp() { return NULL; }
	cVector&	getY2p() { cVector u; return u; }
	wsUintCst		getAnaSampSize() { return 0; }
	//	MAT_t		getDM(wsUint N_idx=0) { return Ra_fDMs[N_idx]; }
	wsMat		get0DM() { return NULL; }
	const char*	getIsInc() { return NULL; }

};

#endif

} // End namespace ONETOOL
