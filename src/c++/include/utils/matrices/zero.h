#pragma once
#include "utils/matrix.h"
#include "global/analysis.h"

namespace ONETOOL {

class cZroMatrix : public cMatrix
{
	void			_init(wsUint N_sz);
public:
	cZroMatrix();
	cZroMatrix(wsUint N_sz);
	cZroMatrix(cZroMatrix &&V_noname);
	~cZroMatrix();

	virtual wsMat	get();

	/* Pre-defined operations */
	wsRealCst		sum();
	bool		sym();
	cZroMatrix&	inv();
	void		rem();
	cVector		sumR(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	cVector		sumC(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	void		file(wsStrCst S_ext, wsUint N_prec=0);
	void		file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/);
	void		file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/);
	cMatrix&	operator-(wsReal R); /* cSymMatrix */
	cMatrix&	operator+(wsReal R); /* cSymMatrix */
	cMatrix&	operator*(wsReal R); /* cDiagMatrix */
	cVector		operator*(cVector &V);
	cMatrix&	addDiag(wsRealCst R);

	void		setR(wsStrCst S_varName);
	wsReal		tr();
	cVector		diag();
	cZroMatrix&	clone();
	cZroMatrix&	sqrt();
	wsRealCst		mean();
	cVector		eigen(cStdMatrix &M_eVec);
	cVector		eigen();
	cZroMatrix&	subset(char *Ba_YorN, wsUint N_Yval);
};

} // End namespace ONETOOL
