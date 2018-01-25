#pragma once
#include "utils/matrix.h"

namespace ONETOOL {

class cBlkMatrix : public cMatrix
{
	char	B_sym, B_ddealloc;
	wsUint	N_szBlk, N_rep;
public:
	cBlkMatrix(wsUint N_szBlk, wsUint N_rep, char B_sym=0, char B_rand=0);
	cBlkMatrix(cMatrix& C_mat, wsUint N_rep);
	cBlkMatrix(wsMat Ra_blkMat, wsUint N_szBlk, wsUint N_rep, char B_sym);
	cBlkMatrix(cBlkMatrix &&B);
	~cBlkMatrix();

	void		init(wsMat Ra_blkMat, wsUint N_szBlk, wsUint N_rep, char B_sym, char B_rand=0);
	void		setDontDealloc() { B_ddealloc = 1; }
	wsUintCst		blk() const;
	wsUintCst		rep() const;

	cMatrix&	operator-(cSymMatrix &S);

	cStdMatrix&	operator*(cStdMatrix &C);
	cMatrix&	operator*=(cStdMatrix &C);

	/* Special operations */
	bool		sym();
	cStdMatrix	mt(const cStdMatrix &C);
	void		sub(cSymMatrix &S);
	cBlkMatrix&	operator=(cBlkMatrix &B);
	cBlkMatrix&	operator=(cBlkMatrix &&B);

	/* Predefined operations */
	wsRealCst		sum();
	cMatrix&	inv();
	void		rem();
	cVector		sumR(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	cVector		sumC(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	void		file(wsStrCst S_ext, wsUint N_prec/*=0*/);
	void		file(wsStrCst S_ext, vStr *Sa_rowNames, wsUint N_prec/*=0*/);
	void		file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/);
	cBlkMatrix&	operator-(wsReal R);
	cBlkMatrix&	operator+(wsReal R);
	cBlkMatrix&	operator*(wsReal R);
	cVector		operator*(cVector &V);
	cBlkMatrix& addDiag(wsRealCst R);
	void		setR(wsStrCst S_varName);
	wsReal		tr();
	cVector		diag();
	cBlkMatrix&	clone();
	cBlkMatrix&	sqrt();
	wsRealCst		mean();
	cVector		eigen(cStdMatrix &Mp_eVec);
	cVector		eigen();
	cMatrix&	subset(char *Ba_YorN, wsUint N_Yval);
};

} // End namespace ONETOOL
