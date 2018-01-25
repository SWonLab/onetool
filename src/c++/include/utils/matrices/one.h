#pragma once
#include "utils/matrix.h"

namespace ONETOOL {

class cOneMatrix : public cMatrix
{
	void			_init(wsUint N_inpRow, wsUint N_inpCol);
public:
	cOneMatrix();
	cOneMatrix(wsUint N_inpRow, wsUint N_inpCol);
	cOneMatrix(cOneMatrix &&V_noname);
	~cOneMatrix();

	virtual wsMat	get();

	/* Pre-defined operations */
	wsRealCst		sum();
	bool		sym();
	cOneMatrix&	inv();
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
	cOneMatrix&	clone();
	cOneMatrix&	sqrt();
	wsRealCst		mean();
	cVector		eigen(cStdMatrix &M_eVec);
	cVector		eigen();
	cOneMatrix&	subset(char *Ba_YorN, wsUint N_Yval);
};

} // End namespace ONETOOL
