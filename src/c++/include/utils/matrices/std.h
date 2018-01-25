#pragma once
#include "utils/matrix.h"

namespace ONETOOL {

wsVec	sseMpV(wsMat Ra_mat, wsUintCst N_r, wsUintCst N_c,
	wsVecCst Ra_v, wsRealCst R_mul=WISARD_NAN);

/* Return Sum(i=1:N_sz) Ra_v[i] */

wsMat	sseMcovCol(wsReal **Ra_mat, wsUint N_row, wsUint N_col);
wsMat	sym_sseCovRow(wsReal **Ra_mat, wsUint N_row, wsUint N_col,
	char B_isComplete=0);
wsMat	sseMcovRow(wsMat Ra_mat, wsUintCst N_row, wsUintCst N_col,
	char B_isComplete=0, char B_retSym=0);
wsSym	sseMcorRow(wsMat Ra_mat, wsUintCst N_r, wsUintCst N_c);
wsVec	sseMsumCol(wsMat Ra_mat, wsUintCst N_row, wsUint N_col=0xffffffff);

/**
 * sseMtV	Performs M'V multiplication, 
 *
 * @param	Ra_m		Former matrix (M) [n*p] to be multiplicated
 * @param	N_r1		# of row M
 * @param	N_c1		# of column M
 * @param	Ra_v		Latter vector (V) [n*1] to be multiplicated
 * @param	N_sz		size of vector
 * @return	(wsReal*)	Result of M'V, vector [p*1] = [p*n] %*% [n*1]
 */
wsReal*	sseMtV(wsMatCst Ra_m, wsUintCst N_r1, wsUintCst N_c1, wsVecCst Ra_v, wsUintCst N_sz);

wsMat	sseMsubsetRect(wsMat Ra_m, wsUintCst N_sz, const char *Ba_filt,
	char N_incVal, wsUintCst N_szAfter);
wsMat	sseMsubsetRect(wsMat Ra_m, wsUintCst N_sz, vBool &Ba_filt,
	char N_incVal, wsUintCst N_szAfter);

/* stdMatrix flag */
typedef enum _xMatOpt {
	MATOPT_DEFAULT		= 0x00,
	MATOPT_DELROW		= 0x01,
	MATOPT_DELNONE		= 0x02,
	MATOPT_TRANSPOSE	= 0x04,
	MATOPT_RANDOMIZE	= 0x08,
	MATOPT_UNIFINIT		= 0x10,
} xMatOpt;

typedef enum _xMatDealloc {
	MATDEL_ALL			= 0x00,
	MATDEL_ROW			= 0x01,
	MATDEL_NONE			= 0x02
} xMatCleanAct;

class cStdMatrix : public cMatrix
{
	xMatCleanAct X_actClean;
public:
	void init(xMatOpt X_flag, wsUint N_R, wsUint N_C, wsReal R_val,
		wsReal **Rp_data=NULL);
	void init(wsUint N_R, wsUint N_C, wsReal **Rp_data=NULL);
	cStdMatrix();
	cStdMatrix(cStdMatrix &&M);
	cStdMatrix(cSpsMatrix &P);
	/* Make matrix and copy data from Rp_data */
	cStdMatrix(wsUint N_R, wsUint N_C, wsReal **Rp_data,
		xMatCleanAct X_cleanAct=MATDEL_ALL);
	cStdMatrix(xMatOpt X_opt, wsUint N_R, wsUint N_C, wsReal **Rp_data);
	cStdMatrix(wsUint N_R, wsUint N_C, char B_rand=0);
	cStdMatrix(wsUint N_R, wsUint N_C, wsReal R_val);
	~cStdMatrix();

	/* Summation */
	cStdMatrix&	operator+=(cStdMatrix &C);
	cStdMatrix&	operator+=(cSymMatrix &C);
	cStdMatrix&	operator+=(cDiagMatrix &C);
	cStdMatrix&	operator+=(cBlkMatrix &C);
	cStdMatrix&	operator+=(wsReal R);
	/* Subtraction */
	cStdMatrix&	operator-=(cStdMatrix &C);
	cStdMatrix&	operator-=(cSymMatrix &C);
	cStdMatrix&	operator-=(cDiagMatrix &C);
	cStdMatrix&	operator-=(cBlkMatrix &C);
	cStdMatrix&	operator-=(wsReal R);
	cStdMatrix	operator-(const cStdMatrix &C);
	/* Multiplication */
	cStdMatrix&	operator*=(cStdMatrix &C);	///< [N*Q] %*% [Q*R] = [N*R]
	cStdMatrix&	operator*=(cSymMatrix &C);	///< [N*Q] %*% [Q*Q] = [N*Q]
	cStdMatrix&	operator*=(cDiagMatrix &C);	///< [N*Q] %*% [Q*Q] = [N*Q]
	cStdMatrix&	operator*=(cBlkMatrix &C);	///< [N*Q] %*% [Q*Q] = [N*Q]
	cStdMatrix	operator*(const cStdMatrix &C);	///< [N*Q] %*% [Q*R] = [N*R]
	cStdMatrix	operator*(const cSymMatrix &C);		///< [N*Q] %*% [Q*Q] = [N*Q]
	cStdMatrix	operator*(const cDiagMatrix &C);	///< [N*Q] %*% [Q*Q] = [N*Q]
	cStdMatrix	operator*(const cBlkMatrix &C);	///< [N*Q] %*% [Q*Q] = [N*Q]
	cStdMatrix&	operator*=(wsReal R);		///< [N*Q] *   [1]   = [N*Q]
	/* Division */
	cStdMatrix& operator/=(wsReal R);

	/* Special operation */
	void		setClean(xMatCleanAct X_opt=MATDEL_NONE) {
		X_actClean = X_opt;
	}
	xMatCleanAct	getClean() { return X_actClean; }
	cStdMatrix	transpose();
	cStdMatrix*	transposeNew();
	cSymMatrix	MMt(const cStdMatrix &M);
	cStdMatrix	MMt(const cStdMatrix &M, const cStdMatrix &N);
	cStdMatrix	MMt(const cSymMatrix &S, const cStdMatrix &M);
	cStdMatrix	MMt(const cSymMatrix &S, cSpsMatrix &P);
	cStdMatrix	MV(const cStdMatrix &M, const cVector &V);
	cSymMatrix	MMt(cMatrix &M);
	cSymMatrix	MMt(const cBlkMatrix &M);
	cSymMatrix	MMt(const cSymMatrix &M);
	cSymMatrix	MMt(const cIdtMatrix &M);
	cSymMatrix	MMt(const cDiagMatrix &M);
	cStdMatrix	MMt(const cDiagMatrix &D, const cStdMatrix &M);
	cSymMatrix&	MMt(const cVector &V);
	cVector		MV(cDiagMatrix &M, cVector &V);
	cVector		MV(cSymMatrix &M, cVector &V);
	cStdMatrix	Mt(const cStdMatrix &M);
	cSymMatrix	Mts(const cStdMatrix &M); /* Do Mt, but assume return is sym */
	cSymMatrix	Ms(const cStdMatrix &M); /* Do *, but assume return is sym */
	cSymMatrix	Mt();
	cVector		tV(wsReal *V, wsUint L);
	cVector		tV(cVector &V);
	cStdMatrix	tM(cStdMatrix &M);
	cStdMatrix	tM(cSymMatrix &S);
	cStdMatrix	tM(cDiagMatrix &D);
	cSymMatrix	tM(wsRealCst R_mul=W0);
	void		getR(wsStrCst S_varName);
	wsRealCst	normF(cStdMatrix &M);
	wsRealCst	normF();
	cSymMatrix	proj();	/* Get projection I-X (X'X)^-1 X' : Assume it is X' */
	cStdMatrix&	normalize(wsRealCst R_mul=W1);
	cStdMatrix&	operator=(cStdMatrix &M);
	cStdMatrix&	operator=(cStdMatrix &&M);
	cStdMatrix	subsetPtr(wsUint N_s, wsUint N_e);
	cStdMatrix	subsetPtr(vInt& Nv_idx);
	cStdMatrix	subsetCol(vInt& Nv_idx);
	wsVec		operator[](wsUint N_rIdx);

	/**
	 * cStdMatrix::subIself Performs I - M operation to itself
	 *
	 * @return    (void)
	 */
	void		subIself();
	/**
	 * cStdMatrix::subI     Performs I - M operation and return it
	 *
	 * @return    (void)
	 */
	cStdMatrix	subI();
	cStdMatrix	divide(wsUint* Na_idx, wsUint N_sz, cStdMatrix* Mp_rest=NULL);
	cStdMatrix	divideCol(wsUint* Na_idx, wsUint N_sz, cStdMatrix* Mp_rest=NULL);
	void		trcLess(wsRealCst R_val, wsUint N_from=0);
	cStdMatrix	log1exp();
	/**
	 * cStdMatrix::r2v      Extracts a row in matrix and return it as cVector
	 *
	 * @param     N_idx     An index of row to be a vector, should between 0 and (#row-1)
	 * @return    (cVector) Row vector of matrix, data is duplicated
	 */
	cVector		r2v(int N_idx);
	/**
	 * cStdMatrix::r2v_ptr  Pointing a row in matrix and return it as cVector.
	 *                      Note that the vector returned here is non-deallocatable.
	 *
	 * @param     N_idx     An index of row to be a vector, should between 0 and (#row-1)
	 * @return    (cVector) Row vector of matrix, data is referenced
	 */
	cVector		r2v_ptr(int N_idx);
	/**
	 * cStdMatrix::tr2     Computes tr(M %*% M), thus the matrix should be square matrix
	 *
	 * @return    (REAL_c) Result of tr(M %*% M)
	 */
	wsRealCst	tr2();
	cVector		rSS();	/* Row-wise sum of squared */
	cStdMatrix	rDiv(const cVector& V_divider);	/* Row-wise division */
	/**
	 * cStdMatrix::setRow	Set specific row to have same specific values
	 *
	 * @param     N_row			An index of row, should be placed between 0 and (#row-1)
	 * @param     R_val			Value to be set to the row, any real value
	 * @return    (void)
	 */
	void		setRow(wsUint N_row, wsReal R_val);
	/**
	 * cStdMatrix::setCol	Set specific column to have same specific values
	 *
	 * @param     N_col			An index of column, should be placed between 0 and (#col-1)
	 * @param     R_val			Value to be set to the column, any real value
	 * @return    (void)
	 */
	void		setCol(wsUint N_col, wsReal R_val);	
	/**
	 * cStdMatrix::addCol	Add a value to specific column
	 *
	 * @param     N_col			An index of column, should be placed between 0 and (#col-1)
	 * @param     R_val			Value to be added to the column, any real value
	 * @return    (void)
	 */
	void		addCol(wsUint N_col, wsReal R_val);	
	/**
	 * cStdMatrix::vectorize Do vectorize the matrix, just same as vec() operator
	 *
	 * @return    (cVector)  Vectorized result, data is duplicated 
	 */
	cVector		vectorize(); /* vec() operator */
	/**
	 * cStdMatrix::trMtSM  Compute the trace of t(M) %*% S %*% M,
	 *                     M is this matrix and S is given matrix.
	 *                     Note that M is firstly transposed, thus
	 *                     the dimension of S should be #row * #row
	 *
	 * @param     S        A symmetric matrix to be included to the computation
	 * @return    (REAL_c) Computed trace of t(M) %*% S %*% M
	 */
	wsRealCst		trMtSM(cSymMatrix &S);
	/**
	 * cStdMatrix::trMSMt  Compute the trace of M %*% S %*% t(M),
	 *                     M is this matrix and S is given matrix.
	 *                     Note that M is not transposed, thus
	 *                     the dimension of S should be #col * #col
	 *
	 * @param     S        A symmetric matrix to be included to the computation
	 * @return    (REAL_c) Computed trace of M %*% S %*% t(M)
	 */
	wsRealCst		trMSMt(cSymMatrix &S);
	/**
	 * cStdMatrix::diagMtSM Compute the diagonals of t(M) %*% S %*% M,
	 *                      M is this matrix and S is given matrix.
	 *                      Note that M is firstly transposed, thus
	 *                      te dimension of S should be #row * #row
	 *
	 * @param     S         A symmetric matrix to be included to the computation
	 * @return    (cVector) Computed diagonals with cVector from t(M) %*% S %*% M
	 */
	cVector		diag(cSymMatrix &S);
	cVector		diag(cDiagMatrix &D);
	cStdMatrix	krnk(cStdMatrix &M);		/* X %x% M */
	cStdMatrix	mldivide(cStdMatrix& M, char B_trans=1);	/* (X'X)^-1X'y */
	cVector		mldivide(cVector& V, char B_trans=1);		/* (X'X)^-1X'y */

	/* Predefined operation */
	wsRealCst	sum();
	wsRealCst	ssum();
	bool		sym();
	cStdMatrix&	inv();
	void		rem();
	/**
	 * cStdMatrix::sumR     Performs row-wise sum of this matrix
	 *
	 * @param     Rp_sumAll If this assigned, sum of matrix will be returned.
	 *                      R_mul will also affects to this
	 * @param     R_mul     If this assigned, result*R_mul will be returned
	 * @return    (cVector) Row-wise sum vector of matrix
	 */
	cVector		sumR(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	/**
	 * cStdMatrix::sumC     Performs column-wise sum of this matrix
	 *
	 * @param     Rp_sumAll If this assigned, sum of matrix will be returned.
	 *                      R_mul will also affects to this
	 * @param     R_mul     If this assigned, result*R_mul will be returned
	 * @return    (cVector) Column-wise sum vector of matrix
	 */
	cVector		sumC(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	/**
	 * cStdMatrix::asumC     Performs column-wise absoulte sum of this matrix
	 *
	 * @param     Rp_sumAll If this assigned, absolute sum of matrix will be returned.
	 *                      R_mul will also affects to this
	 * @param     R_mul     If this assigned, result*R_mul will be returned
	 * @return    (cVector) Column-wise absoulte sum vector of matrix
	 */
	cVector		asumC(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN);
	cStdMatrix	acDiv(cVector &V_divider);
	/**
	 * cStdMatrix::rmerge      Merges input matrix to right-side and return new cStdMatrix
	 *
	 * @param     M_merged     cStdMatrix to be merged to right-side, #row must be same
	 * @return    (cStdMatrix) A merged matrix, data is duplicated
	 */
	cStdMatrix	rmerge(cStdMatrix& M_merged);
	/**
	 * cStdMatrix::rmerge      Merges input matrix to bottom-side and return new cStdMatrix
	 *
	 * @param     M_merged     cStdMatrix to be merged to right-side, #col must be same
	 * @return    (cStdMatrix) A merged matrix, data is duplicated
	 */
	cStdMatrix	bmerge(cStdMatrix& M_merged);
	/**
	 * cStdMatrix::rRepeat     Repeats entire matrix k times for row-direction
	 *
	 * @param     N_rep        Number of repetition
	 * @return    (cStdMatrix) A merged matrix, data is duplicated
	 */
	cStdMatrix	rRepeat(wsUint N_rep);
	/**
	 * cStdMatrix::ssRow       Extracts multiple rows in matrix and return it as cStdMatrix
	 *
	 * @param     Xv_idx       Indices of row to be a new cStdMatrix, should between 0 and (#row-1)
	 * @return    (cStdMatrix) A new matrix, data is duplicated
	 */
	cStdMatrix	ssRow(vInt Xv_idx);
	/**
	 * cStdMatrix::meanR    Computes row-wise mean of this matrix
	 *
	 * @return    (cVector) Row-wise mean vector of matrix
	 */
	cVector		meanR();
	/**
	 * cStdMatrix::meanC    Computes column-wise mean of this matrix
	 *
	 * @return    (cVector) Column-wise mean vector of matrix
	 */
	cVector		meanC();
	cSymMatrix	covR();
	void		file(wsStrCst S_ext, wsUint N_prec=0);
	void		file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec=0);
	void		file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec=0);
	cStdMatrix&	operator-(wsReal R);
	cStdMatrix&	operator+(wsReal R);
	cStdMatrix&	operator*(wsReal R);
	cVector		operator*(cVector &V);		///< [N*Q] %*% [Q*1] = [N*1]
	cStdMatrix&	addDiag(wsRealCst R);
	void		setR(wsStrCst S_varName);
	wsReal		tr();
	cVector		diag();
	cVector		diag2();
	cStdMatrix&	clone();
	cStdMatrix&	sqrt();
	wsRealCst	mean();
	cVector		eigen(cStdMatrix &M_eVec);
	cVector		eigen();
	cStdMatrix&	subset(char *Ba_YorN, wsUint N_Yval);
	cVector		vec();	/* vec() operator */
	// min, max
	wsRealCst	getMin(wsUint* Np_idxR=NULL, wsUint* Np_idxC=NULL);
	wsRealCst	getMax(wsUint* Np_idxR=NULL, wsUint* Np_idxC=NULL);
	// Subtract ith vector's value from ith row's values
	cStdMatrix	subR(cVector& V_val);
};

} // End namespace ONETOOL
