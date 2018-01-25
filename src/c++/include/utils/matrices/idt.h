#pragma once
#include "utils/matrix.h"

namespace ONETOOL {

class cZroMatrix;
class cIdtMatrix : public cMatrix
{
	void			_init(wsUint N_sz);
public:
	cIdtMatrix(wsUint N_sz);
	cIdtMatrix(cIdtMatrix &&V_noname);
	~cIdtMatrix();

	virtual wsMat	get();

	/* Special operations */
	wsRealCst		normF(cStdMatrix &M);
	wsRealCst		normF(cSymMatrix &S);
	wsRealCst		normF();
	cStdMatrix&		mt(cStdMatrix &M);
	void			sub(cSymMatrix &S);
	cSymMatrix		operator-(const cSymMatrix &S);
	cZroMatrix		operator-(const cIdtMatrix &I);
	cIdtMatrix&		operator=(cIdtMatrix &I);
	cIdtMatrix&		operator=(cIdtMatrix &&I);
	cSymMatrix&		operator*(cSymMatrix &S);
	cDiagMatrix		krnk(cDiagMatrix &D);
	cStdMatrix		rmerge(cStdMatrix& M);
	cStdMatrix		bmerge(cStdMatrix& M);

	/* Pre-defined operations */
	wsRealCst			sum();
	bool			sym();
	cIdtMatrix&		inv();
	void			rem();
	cVector			sumR(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	cVector			sumC(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	void			file(wsStrCst S_ext, wsUint N_prec=0);
	void			file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/);
	void			file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/);
	cMatrix&		operator-(wsReal R); /* cSymMatrix */
	cMatrix&		operator+(wsReal R); /* cSymMatrix */
	cMatrix&		operator*(wsReal R); /* cDiagMatrix */
	cVector			operator*(cVector &V);
	cMatrix&		addDiag(wsRealCst R);

	void			setR(wsStrCst S_varName);
	wsReal			tr();
	cVector			diag();
	cIdtMatrix&		clone();
	cIdtMatrix&		sqrt();
	wsRealCst			mean();
	cVector			eigen(cStdMatrix &M_eVec);
	cVector			eigen();
	cIdtMatrix&		subset(char *Ba_YorN, wsUint N_Yval);
};

} // End namespace ONETOOL
