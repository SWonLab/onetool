#pragma once
#include "utils/matrix.h"

namespace ONETOOL {

class cSubMatrix : public cMatrix
{
	vMatrix		V_mats;
public:
	cSubMatrix(vMatrix &V_inp);
	cSubMatrix(wsUint N_mat, ...);
	~cSubMatrix();

	vMatrix&		get();
	
	/* Predefined operations */
	wsRealCst			sum();
	bool			sym();
	cMatrix&		inv();
	void			rem();
	cVector			sumR(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	cVector			sumC(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	void			file(wsStrCst S_ext, wsUint N_prec);
	void			file(wsStrCst S_ext, vStr *Sa_rowNames, wsUint N_prec=0);
	void			file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/);
	cMatrix&		operator-(wsReal R);
	cMatrix&		operator+(wsReal R);
	cMatrix&		operator*(wsReal R);
	cVector			operator*(cVector &V);
	cMatrix&		addDiag(wsRealCst R);
	void			setR(wsStrCst S_varName);
	wsReal			tr();
	cVector			diag();
	cSubMatrix&		clone();
	cSubMatrix&		sqrt();
	wsRealCst			mean();
	cVector			eigen(cStdMatrix &Mp_eVec);
	cVector			eigen();
	cSubMatrix&		subset(char *Ba_YorN, wsUint N_Yval);

	/* Special operations */
	cStdMatrix&		mt(cStdMatrix &C);
	void			sub(cSymMatrix &S);
	cSubMatrix&		operator=(cSubMatrix &B);
};

} // End namespace ONETOOL
