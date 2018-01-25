#pragma once
#include "utils/matrix.h"

namespace ONETOOL {

class cDiagMatrix : public cMatrix
{
	char	B_dDealloc;
public:
	cDiagMatrix();
	cDiagMatrix(wsUint N_sz, char B_rand=0);
	cDiagMatrix(wsUint N_sz, wsReal *Ra_V, char B_inpDDealloc=0);
	cDiagMatrix(wsReal R_V, wsUint N_sz);
	cDiagMatrix(cDiagMatrix &&D);
	cDiagMatrix(cDiagMatrix &D);
	cDiagMatrix(cVector &V);
	~cDiagMatrix();
	void			init(wsUint N_sz, wsReal *Ra_V=NULL, char B_rand=0,
		char B_inpDDealloc=0);
	void			init(cDiagMatrix &D);

	/* Multiplication */
	cStdMatrix		operator*(const cStdMatrix &C);
	cSymMatrix		operator*(const cSymMatrix &C);
	cStdMatrix		operator*(const cBlkMatrix &C);
	cDiagMatrix&	operator*=(cDiagMatrix &D);
	cDiagMatrix&	operator*=(wsRealCst C);
	cDiagMatrix		operator*(const cDiagMatrix &D);

	/* Special operations */
	void			setDontDealloc() { B_dDealloc = 1; }
	cStdMatrix		mt(const cStdMatrix &M);
	cSymMatrix		MMt(const cSymMatrix &S);
	cDiagMatrix		MMt(const cDiagMatrix &D);
	cDiagMatrix&	operator=(cDiagMatrix &D);
	cDiagMatrix&	operator=(cVector &V);
	cDiagMatrix&	operator=(cDiagMatrix &&D);
	cSymMatrix		operator-(const cSymMatrix &S);
	cDiagMatrix&	operator+=(wsReal N);
	cDiagMatrix&	operator+=(cDiagMatrix &D);
	wsReal			trMMt(cSymMatrix &S);
	wsReal			detL();
	wsReal			tr2();
	cVector			MV(cDiagMatrix &D, cVector &V);
	cDiagMatrix		subset(wsUint N_start, wsUint N_sz);
	cDiagMatrix		inv2();
	wsRealCst		Lsum();
	wsRealCst		prod();
	cDiagMatrix		krnk(const cIdtMatrix &I);
	wsUintCst		nn0();
	cDiagMatrix		v_kv(wsRealCst R_k);

	/* Predefined operations */
	wsRealCst		sum();
	bool			sym();
	cDiagMatrix&	inv();
	void			rem();
	cVector			sumR(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	cVector			sumC(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	void			file(wsStrCst S_ext, wsUint N_prec/*=0*/);
	void			file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/);
	void			file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/);
	cDiagMatrix&	operator-(wsReal R);
	cDiagMatrix&	operator+(wsReal R);
	cDiagMatrix		operator+(cDiagMatrix &D);
	cDiagMatrix&	operator*(wsReal R);
	cVector			operator*(cVector &V);
	cDiagMatrix&	addDiag(wsRealCst R);
	void			setR(wsStrCst S_varName);
	wsReal			tr();
	cVector			diag();
	cDiagMatrix&	clone();
	cDiagMatrix&	sqrt();
	wsRealCst		mean();
	cVector			eigen(cStdMatrix &M_eVec);
	cVector			eigen();
	cDiagMatrix&	subset(char *Ba_YorN, wsUint N_Yval);
};

} // End namespace ONETOOL
