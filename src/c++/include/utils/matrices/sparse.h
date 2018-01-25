#pragma once
#include "utils/matrix.h"

namespace ONETOOL {

class cSpsMatrix : public cMatrix
{
	char	B_dtdealloc;
	char	B_sym;
	wsUint*	Na_rEnd;
	wsUint*	Na_cIdx;
	wsUint	N_size;
public:
	cSpsMatrix(wsUint N_r, wsUint N_c, wsUint N_sz=0);
	cSpsMatrix(cStdMatrix& M);
	cSpsMatrix(cSymMatrix& S);
	cSpsMatrix(cDiagMatrix& D);
	cSpsMatrix(cSpsMatrix& P);
	cSpsMatrix(cSpsMatrix& P, wsVec Ra_new);
	cSpsMatrix(cSpsMatrix&& P);
        cSpsMatrix(cSpsMatrix&& P, wsVec Ra_new);
	void init(wsUint N_r, wsUint N_c, wsUint N_sz=0);

	/* Special operations */
	wsUintCst			size();
	wsUint*			getRowEnds() const;
	wsUint*			getColIdxs() const;
	cStdMatrix		toStd();
	cStdMatrix		operator*(cStdMatrix& M);
	cSpsMatrix		operator*(cDiagMatrix& D);
	cSpsMatrix		operator*(cIdtMatrix& I);
	cSpsMatrix		operator*(cSymMatrix& S);
	cSpsMatrix&		addDiag(wsRealCst R);
	void			setDontDealloc();

	/* Predefined operations */
	wsRealCst			sum();
	bool			sym();
	cMatrix&		inv();
	void			rem();
	cVector			sumR(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	cVector			sumC(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	void			file(wsStrCst S_ext, wsUint N_prec/*=0*/);
	void			file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/);
	void			file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/);
	cMatrix&		operator-(wsReal R);/*NS*/
	cMatrix&		operator+(wsReal R);/*NS*/
	cSpsMatrix&		operator*(wsReal R);/*NS*/
	cVector			operator*(cVector &V);/*NS*/
	void			setR(wsStrCst S_varName);/*NS*/
	wsReal			tr();
	cVector			diag();
	cSpsMatrix&		clone();
	cMatrix&		sqrt();/*NS*/
	wsRealCst			mean();
	cVector			eigen(cStdMatrix &M_eVec);/*NS*/
	cVector			eigen();/*NS*/
	cSpsMatrix&		subset(char *Ba_YorN, wsUint N_Yval);
};

} // End namespace ONETOOL
