#pragma once
#ifndef __WISARD_MFQLS_H__
#define __WISARD_MFQLS_H__
#include "global/analysis.h"
#include "fqls.h"
#include <vector>
using namespace std;

namespace ONETOOL {

class cSymMatrix;

#define MQLS_NTEST	3

typedef struct _xMqlsRes {
	wsReal	R_t[MQLS_NTEST];
	wsReal	R_p[MQLS_NTEST];
	char	B_ginv;
} xMqlsRes;

typedef struct _xMqlsTest
{
	wsStrCst		Sp_name;
	vInt&		Xa_snps;
	xMqlsRes	X_res;

	cIO			*Cp_IO;
	char		*Ba_misPheno;
	wsUint		N_pheno;
	wsUint		N_anaSamp;

	vector<cSymMatrix*>&	Ma_TtAphiAtT_inv;
	vector<wsMat>&			Ma_TtA;
} xMqlsTest;

class cSetManagerAnalysis;
class cVariantMap;
class cStdMatrix;

/* !!! FQLS-specific analysis !!! */
#if (TOOLSET_TYPE == TOOLSET_MFQLS) || ((TOOLSET_TYPE & 0x100) == 0x100)

class cMFqlsAnalysis : public cAnalysis
{
	//MAT_t	Ra_A;	/* [N*N] I - 1(1'��^-11)^-11'��^-1 */
	vector<wsMat>		Ma_TtA; /* [Q*N] T'A */
	vector<cSymMatrix*>	Ma_TtAphiAtT_inv;
	//SYM_t	Ra_TtAphiAtT_inv; /* [Q*Q] (T'A��A'T)^-1 */
	wsUintCst	N_pheno;/* Q */
	wsUint	N_anaSamp;
	const char	*Ba_misPhCv;
	const char	*Ba_misCov;
	cSetManagerAnalysis *Cp_anaGS;
	cVariantMap			*Cp_mm;
	cStdMatrix			M_Y, M_Yf1, M_Yf2;
//	MAT_t	Ra_Y; /* [Q*N] Y-E(Y) */
	static int		_doTest(xMqlsTest *Xp_t);
//	int				_doTest(str_c Sp_name, vInt& Xa_snps, xMqlsRes &X_res);
	wsReal*			_buildFamI(xFamily &X_fam, wsUint N_szFam,
		wsUint *Na_idxMem, wsReal *Ra_adjY, wsReal **Ra_cor,
		wsReal *Ra_newI_a, wsRealCst *Ra_fullY, wsReal R_h2, wsReal R_p1,
		wsReal R_p2);
	wsReal*			_buildFamII(xFamily &X_fam, wsUint N_szFam,
		wsUint *Na_idxMem, wsReal *Ra_adjY, wsReal **Ra_cor,
		wsReal *Ra_newII_a, wsRealCst *Ra_fullY, wsReal R_h2, wsReal R_p1,
		wsReal R_p2);
	cStdMatrix		_getFQLSphe(cVector& V_prev);
	resPAfmlar*	_getProbandStat(wsUint *N_idxProband,
		wsUint N_pb, wsRealCst *Ra_fullY, wsReal R_h2, wsReal **Ra_cor,
		wsReal R_p1, wsReal R_p2);
public:
	cMFqlsAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpCorr,
		cSetManagerAnalysis *Cp_inpAnaGS, cVariantMap *Cp_inpMM);
	~cMFqlsAnalysis();
	void run();
	vector<wsMat>& getTtA() { return Ma_TtA; }
	vector<cSymMatrix*>&	getTtAphiAtT() { return Ma_TtAphiAtT_inv; }
	cSetManagerAnalysis*	getSet() { return Cp_anaGS; }
};

#else

/* DUMMY DEFINITION */
class cMFqlsAnalysis : public cAnalysis
{
public:
	cMFqlsAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpCorr,
		cSetManagerAnalysis *Cp_inpAnaGS, cVariantMap *Cp_inpMM) : cAnalysis(Cp_inpIO) {}
	~cMFqlsAnalysis() {}
	void run() {}
};

#endif

} // End namespace ONETOOL

#endif
