#pragma once
#include "utils/matrix.h"

namespace ONETOOL {

wsSym	sseSsubsetRect(wsSym Ra_m, wsUintCst N_sz, const char *Ba_filt,
	char N_incVal, wsUintCst N_szAfter);

/**
 * sseSS	Performs SS multiplication, return matrix will be standard
 *
 * @param	Ra_sym1	Former symmetric matrix (S) [n*n] to be multiplicated
 * @param	N_s1	size of S
 * @param	Ra_sym2	Latter matrix (T) [n*n] to be multiplicated
 * @param	N_s2	size of T
 * @return	(MAT_t)	Result of SS [n*n] = [n*n] %*% [n*n]
 */
wsMat	sseSpS(wsSym Ra_sym1, wsUint N_s1, wsSym Ra_sym2, wsUint N_s2);

/**
 * sseSS	Performs SS multiplication, return matrix will be symmetric
 *
 * @param	Ra_sym1	Symmetric matrix (S) [n*n] to be multiplicated
 * @param	N_s1	size of S
 * @return	(SYM_t)	Result of SS [n*n] = [n*n] %*% [n*n]
 */
wsSym	sseSpS(wsSym Ra_sym1, wsUint N_s1);

/**
 * sseMpS	Performs MS multiplication, return matrix will be standard
 *
 * @param	Ra_mat	Former matrix (M) [n*p] to be multiplicated
 * @param	N_row	# row of M
 * @param	N_col	# column of M
 * @param	Ra_s	Latter symmetric matrix (S) [p*p] to be multiplicated
 * @param	N_sz	size of S
 * @return	(MAT_t)	Return of MS [n*p] = [n*p] %*% [p*p]
 */
wsMat	sseMpS(wsMat Ra_mat, wsUint N_row, wsUint N_col, wsSym Ra_s,
	wsUint N_sz);

/**
 * sseMSMt	Performs MSM' multiplication, return matrix will be symmetric
 *
 * @param	Ra_mat		Former matrix (M) [n*p] to be multiplicated
 * @param	N_row		# row of M
 * @param	N_col		# column of M
 * @param	Ra_symMat	Latter symmetric matrix (S) [p*p] to be multiplicated
 * @param	N_sz		size of S
 * @return	(SYM_t)		Return of MSM' [p*p] = [n*p] %*% [p*p] %*% [p*n]
 */
wsSym	sseMSMt(wsMat Ra_mat, wsUint N_row, wsUint N_col, wsSym Ra_symMat,
	wsUint N_sz);

/**
 * sseMSMt	Performs SSS multiplication, return matrix will be symmetric
 *          note that the first and the third matrices are same matrix
 *
 * @param	Ra_s1		Former sym matrix (S) [p*p] to be multiplicated
 * @param	N_s1		size of the first & third matrix
 * @param	Ra_s2		Latter symmetric matrix (S) [p*p] to be multiplicated
 * @param	N_s2		size of the second matrix
 * @return	(SYM_t)		Return of SSS [p*p] = [p*p] %*% [p*p] %*% [p*p]
 */
wsSym	sseSSS(wsSym Ra_s1, wsUint N_s1, wsSym Ra_s2, wsUint N_s2);

/**
 * sseSpM	Performs SM multiplication, return matrix will be standard
 *
 * @param	Ra_symMat	Former symmetric matrix (S) [n*n] to be multiplicated
 * @param	N_sz		size of S
 * @param	Ra_mat		Latter matrix (M) [n*p] to be multiplicated
 * @param	N_row		# row of M
 * @param	N_col		# column of M
 * @return	(MAT_t)		Return of SM [n*p] = [n*n] %*% [n*p]
 */
wsMat	sseSpM(wsSym Ra_symMat, wsUint N_sz, wsReal **Ra_mat, wsUint N_row,
	wsUint N_col);

/**
 * sseSpMt	Performs SM' multiplication, return matrix will be standard
 *
 * @param	Ra_m1	Former symmetric matrix (S) [n*n] to be multiplicated
 * @param	N_sz	size of S
 * @param	Ra_m2	Latter matrix (M) [p*n] to be multiplicated
 * @param	N_r2	# row of M
 * @param	N_c2	# column of M
 * @return	(SYM_t)	Return of SM' [n*p] = [n*n] %*% [n*p]
 */
wsSym	sseSpMt(wsSym Ra_m1, wsUint N_sz, wsMat Ra_m2, wsUint N_r2,
	wsUint N_c2);

/**
 * sym_sseMS	Performs MS multiplication, but assumes the return matrix
 *				will be symmetric
 *
 * @param	Ra_mat	Former matrix (M) [n*p] to be multiplicated
 * @param	N_row	# row of M
 * @param	N_col	# column of M
 * @param	Ra_s	Latter symmetric matrix (S) [p*p] to be multiplicated
 * @param	N_sz	size of S
 * @return	(SYM_t)	Return of MS [n*p] = [n*p] %*% [p*p]
 */
wsSym	sym_sseMS(wsMat Ra_mat, wsUint N_row, wsUint N_col, wsSym Ra_s,
	wsUint N_sz);

/* Calculate Ax
 *  - A is symmetric matrix
 *  - x is transposed vector (=row vector)
 */
wsVec	sseSpV(wsSymCst Ra_symMat, wsUintCst N_sz, wsVecCst Ra_vec);
wsReal*	sseVpS(wsReal *Ra_vec, wsUint N_sz, wsSym Ra_symMat);
wsReal*	sseVpS(wsReal *Ra_v, wsUint N_vLen, wsSym Ra_symMat, wsUint N_sz,
	wsReal *Ra_mask=NULL, wsReal *Rp_r=NULL);

/**
 * sseVpS	Performs VS multiplication
 *
 * @param	Ra_v		Former vector (V) [n*1] to be multiplicated
 * @param	N_vLen		size of V
 * @param	Ra_symMat	Latter symmetric matrix (S) [n*n] to be multiplicated
 * @param	Ra_mask		(optional) Computation mask [size N_sz],
 *						0 then skip, otherwise compute
 * @return	(VEC_t)		Return of VS [n*1] = [1*n] %*% [n*n]
 */
wsVec	sseVpS(wsVec Ra_v, wsUint N_vLen, wsSym Ra_symMat,
	wsReal *Ra_mask=NULL);

class cSymMatrix : public cMatrix
{
	char B_ddealloc;
public:
	cSymMatrix();
	cSymMatrix(wsUint N_sz, char B_rand=0);
	cSymMatrix(wsSym Ra_mat, wsUintCst N_sz, char B_inpDDealloc=0);
	cSymMatrix(wsRealCst R_Vdiag, wsRealCst R_VoffDiag, wsUintCst N_sz);
	cSymMatrix(wsRealCst R_val, wsUintCst N_sz);
	cSymMatrix(wsVec Ra_mat, wsUintCst N_sz);
	cSymMatrix(cSymMatrix &&S);
	~cSymMatrix();

	void		init(wsUint N_sz);
	void		init(wsVec Ra_mat, wsUint N_sz);
	void		init(wsSym Ra_mat, wsUint N_sz);
	void		init(wsSym Ra_mat, wsReal R_val, wsUint N_sz, char B_rand=0,
		char B_inpDDealloc=0);
	void		init(cSymMatrix &S);
	cStdMatrix	operator*(const cStdMatrix& M);
	cStdMatrix	operator*(const cSymMatrix& S);
	cStdMatrix	operator*(const cDiagMatrix& D);
	cStdMatrix	operator*(const cBlkMatrix& B);
	cSymMatrix&	operator*=(wsReal R_val);
	cSymMatrix&	operator+=(wsReal R_val);

	cStdMatrix	operator+(const cStdMatrix& M);
	cSymMatrix	operator+(const cSymMatrix& S);
	cSymMatrix	operator+(const cDiagMatrix& D);
	cStdMatrix	operator+(const cBlkMatrix& B);

	cSymMatrix&	operator-=(cSymMatrix &S);
	cSymMatrix&	operator+=(cSymMatrix &S);
	cSymMatrix&	operator+=(cIdtMatrix &I);
	cSymMatrix&	operator-=(wsReal R_val);
	cSymMatrix	operator-(const cSymMatrix &S);
	cSymMatrix	operator-(const cIdtMatrix &I);
	cStdMatrix	operator-(const cStdMatrix &M);
	cSymMatrix&	operator/=(wsRealCst R_val);

	/* Special operations */
	wsRealCst	det();
	wsRealCst	detL();
	wsRealCst	normF(cStdMatrix &M);
	wsRealCst	normF(cSymMatrix &S);
	wsRealCst	normF();
	cStdMatrix	mt(const cStdMatrix &M);
	void		sub(cSymMatrix &S);
	cSymMatrix&	operator=(cSymMatrix &S);
	cSymMatrix&	operator=(cSymMatrix &&S);
	cVector		diagInv2vec();
	wsRealCst	tr2();
	wsReal		tr(cSymMatrix &S);
	void		setDontDealloc() { B_ddealloc = 1; }
	cSymMatrix&	ginv(wsReal *Rp_Ldet=NULL);
	cSymMatrix&	krnk(const cSymMatrix &S);
	cSymMatrix&	krnk(const cIdtMatrix &I);
	cSymMatrix	einv(cSymMatrix *Mp_sqInv=NULL);
	cSymMatrix	MMt(cSymMatrix &S);
	cSymMatrix	Mt();
	void		setDiag(wsReal R_val);
	cSymMatrix&	neg();
	cSymMatrix&	addDiag(wsRealCst R_val);
	/**
	 * cSymMatrix::subIself Performs I-S operation to itself
	 *
	 * @return    (void)
	 */
	void		subIself();
	/**
	 * cStdMatrix::subI     Performs I-S / I-S*mul operation and return it
	 *
	 * @return    (void)
	 */
	cSymMatrix	subI(wsRealCst R_mul=W1);

	/**
	 * cSymMatrix::r2v_ptr  Pointing a row in matrix and return it as cVector.
	 *                      Note that the vector returned here is non-deallocatable.
	 *                      And the size of returned vector reflects the
	 *                      symmetry of matrix.
	 *
	 * @param     N_idx     An index of row to be a vector, should between 0 and (#row-1)
	 * @return    (cVector) Row vector of matrix, data is referenced
	 */
	cVector		r2v_ptr(int N_idx);
	cSymMatrix	tM();

	/* Predefined operations */
	wsRealCst	sum();
	bool		sym();
	void		rem();
	cSymMatrix&	inv();
	cSymMatrix&	inv(char &B_isGinv);
	cSymMatrix&	inv(wsReal **Ra_eVec, wsReal *Ra_eVal, wsUint N_sz);
	cSymMatrix&	inv(wsReal *Rp_Ldet, char *B_isGinv=NULL, char B_halt=0);
	cSymMatrix&	invOnly();
	cVector		sumR(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	cVector		sumC(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	cVector		meanR();
	void		file(wsStrCst S_ext, wsUint N_prec=0);
	void		file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec=0);
	void		file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec=0);
	cSymMatrix&	operator-(wsReal R);
	cSymMatrix&	operator+(wsReal R);
	cSymMatrix&	operator*(wsReal R);
	cVector		operator*(cVector& V);
	void		setR(wsStrCst S_varName);
	wsReal		tr();
	cVector		diag();
	cSymMatrix&	clone();
	cSymMatrix&	sqrt();
	cSymMatrix&	sqrt(wsReal **Ra_eVec, wsReal *Ra_eVal, wsUint N_sz);
	wsRealCst	mean();
	cVector		eigen(cStdMatrix &Mp_eVec);
	cVector		eigen();
	cSymMatrix&	subset(char *Ba_YorN, wsUint N_Yval);
};

} // End namespace ONETOOL
