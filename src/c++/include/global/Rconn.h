#pragma once
#ifndef __WISARD_RCONN_H__
#define __WISARD_RCONN_H__

#ifdef USE_R
#	include <Rinternals.h>
#else
#	define SEXP void*
#endif

namespace ONETOOL {

class cMatrix;
class cRvar
{
	char	*S_varName;
	SEXP	X_var;
public:
	cRvar() {
		S_varName = NULL;
	}
	cRvar(const char *S_newVarName, SEXP X_exp) {
		set(S_newVarName, X_exp);
	}
	~cRvar();
	void set(const char *S_newVarName, SEXP X_exp);
	operator SEXP() {
		return X_var;
	}
	void operator=(cMatrix &C_mat);
};

class cOption;
class cRconn
{
	char	B_init;
	void	_setRpath(char &B_init);
public:
	cRconn();
	~cRconn();
	cRvar&	operator[](wsStrCst S_varName);
	void	init();
	char	isInit() { return B_init; }
	SEXP	Rparse(wsStrCst S_Rcmd, ...);
	SEXP	Reval(wsStrCst S_buf);
	void	sendVector(wsStrCst S_varName, wsReal *Rp_data, wsUint N_sz);
	void	sendVector(const char *S_varName, int* Np_data, wsUint N_sz);
	void	sendVector(wsStrCst S_varName, char *Np_data, wsUint N_sz, char N_varNA=WISARD_NA);
	void	sendVector(wsStrCst S_varName, char **Sp_strs, wsUint N_sz);
	void	sendVector(const char *S_varName, vStr &Xv_strs);
	void	sendMatrix(wsStrCst S_varName, wsReal **Rp_data, wsUint N_row, wsUint N_col);
	void	recvVector(wsStrCst S_varName, wsReal *Rp_data, wsUint N_sz);
	void	recvMatrix(wsStrCst S_varName, wsReal **Rp_data, wsUint N_row, wsUint N_col);
	void	sendOption(cOption &X_opt);
	void	dim(wsStrCst S_varName, wsUint *R, wsUint *C);
};

cRconn& R();

} // End namespace ONETOOL

#endif
