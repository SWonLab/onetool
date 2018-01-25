#pragma once
#ifndef __WISARD_MATRIX_H__
#define __WISARD_MATRIX_H__

#include <vector>
#include "global/common.h"
#include "utils/util.h"
#include "global/io.h"

namespace ONETOOL {

void exportMatrix(const char *S_ext, char **Na_mat, wsUint N_row, wsUint N_col);
void exportMatrix(const char *S_ext, vector<wsFvec>& Na_mat, wsUint N_col,
	vStr *Xp_header=NULL);
void exportMatrix(const char *S_ext, wsMat Rp_mat, wsUint N_row,
	wsUint N_col, vStr *Xp_header, vStr *Xp_rowNames, wsUint N_prec,
	wsReal *Rp_miss = NULL);
void exportMatrix(const char *S_ext, wsMat Rp_mat, wsUint N_row,
	wsUint N_col, vStr *Xp_header);
#ifdef USE_DBL
void exportMatrix(const char *S_ext, wsFmat Rp_mat, wsUint N_row,
	wsUint N_col);
void exportMatrix(const char *S_ext, wsFmat Rp_mat, wsUint N_row,
	wsUint N_col, vStr *Xp_header, vStr *Xp_rowNames, wsUint N_prec,
	wsReal *Rp_miss = NULL);
#endif
void exportMatrix(const char *S_ext, wsMat Rp_mat, wsUint N_row,
	wsUint N_col, vStr *Xp_header, vStr *Xp_rowNames, wsReal *Rp_miss = NULL);
void exportMatrix(const char *S_ext, vector<wsVec>& Rv_mat, wsUint N_col,
	vStr *Xp_header, vStr *Xp_rowNames, wsUint N_prec, wsReal *Rp_miss = NULL);

void exportMatrix(const char *S_ext, wsMat Rp_mat, wsUint N_row,
	wsUint N_col);
void exportMatrix(const char *S_ext, vector<wsVec>& Rv_mat, wsUint N_col);
void exportSpsMat(const char *S_ext, wsReal *Rp_mat, wsUint N_r, wsUint N_c,
	wsUint* Na_rEnd, wsUint* Na_cIdx, vStr *Xp_header = NULL, wsUint N_prec = 0);
void exportSymMat(const char *S_ext, wsMat Rp_mat, wsUint N_sz,
	vStr *Xp_header = NULL, vStr *Xp_rowNames = NULL, wsUint N_prec = 0);
void exportBlkMat(wsStrCst S_ext, wsMat Rp_mat, wsUint N_szBlk,
	wsUint N_rep, wsUint N_prec = 0);
void exportBlkMat(wsStrCst S_ext, wsMat Rp_mat, wsUint N_szBlk,
	wsUint N_rep, vStr *Sa_rowNames = NULL, wsUint N_prec = 0);
void exportDiagMat(wsStrCst S_ext, wsMat Rp_mat, wsUint N_sz,
	wsUint N_prec = 0);
void exportDiagMat(wsStrCst S_ext, wsMat Rp_mat, wsUint N_sz,
	vStr *Sa_rowNames = NULL, wsUint N_prec = 0);
void exportIdtMat(wsStrCst S_ext, wsUint N_sz, vStr *Sa_colNames = NULL, vStr *Sa_rowNames = NULL);
void exportValMat(wsReal R_val, wsStrCst S_ext, wsUint N_sz, vStr *Sa_rowNames = NULL);

class cVector;
class cStdMatrix;
class cSymMatrix;
class cBlkMatrix;
class cDiagMatrix;
class cIdtMatrix;
class cSubMatrix;
class cSpsMatrix;

enum xMatErr
{
	MATERR_DIM_INCORRECT,
	MATERR_TYPE_INCORRECT,
	MATERR_PARAM_INCORRECT,
};

enum xMatType
{
	MATTYPE_STD,
	MATTYPE_SYM,
	MATTYPE_BLK,
	MATTYPE_DIAG,
	MATTYPE_IDT,
	MATTYPE_SUB,
};

class cMatrix
{
protected:
	xMatType			X_type;
	/* Data pointer, actual dimension may be differ with N_row & N_col */
	wsReal				**Ra_data;
	/* VIRTUAL size of matrix */
	wsUint				N_row, N_col;
	/* Raise matrix-specific error */
	void				raise(xMatErr X_err);
	void				__init(wsMat Ra_data, wsUint N_row, wsUint N_col,
		char B_passNull=0);
public:
	cMatrix();
	virtual ~cMatrix();
	xMatType			getType() { return X_type; }
	char				ok() { return Ra_data != NULL; }
	virtual wsMat		get() const;
	wsMatCst			cget() { return (wsMatCst)Ra_data; }
	wsUintCst			row() const;
	wsUintCst			col() const;
	const char			sqr() const;
	const char			sane();
	void				fmtfile(wsStrCst S_fmtExt, ...);

	/* Predefined operations */
	virtual wsRealCst		sum()=0;
	virtual bool		sym()=0;
	virtual void		rem()=0;
	virtual cMatrix&	inv()=0;
	virtual cVector		sumR(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN)=0;
	virtual cVector		sumC(wsReal *Rp_sumAll=NULL, wsReal R_mul=WISARD_NAN)=0;
	virtual void		file(wsStrCst S_ext, wsUint N_prec=0)=0;
	virtual void		file(wsStrCst S_ext, vStr *Sa_rowNames, wsUint N_prec=0)=0;
	virtual void		file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec=0)=0;
	virtual cMatrix&	operator-(wsReal R)=0;
	virtual cMatrix&	operator+(wsReal R)=0;
	virtual cMatrix&	operator*(wsReal R)=0;
	virtual cVector		operator*(cVector &V)=0;
	virtual cMatrix&	addDiag(wsRealCst R) = 0;
	virtual void		setR(wsStrCst S_varName)=0;
	virtual	wsReal		tr()=0;
	virtual	cVector		diag()=0;
	virtual cMatrix&	clone()=0;
	virtual cMatrix&	sqrt()=0;
	virtual wsRealCst	mean()=0;
	virtual	cVector		eigen(cStdMatrix &Mp_eVec)=0;
	virtual	cVector		eigen()=0;
	virtual	cMatrix&	subset(char *Ba_YorN, wsUint N_Yval)=0;
};

typedef vector<cMatrix*>				vMatrix;
typedef vMatrix::iterator				vMatrix_it;
typedef std::map<std::string, vMatrix*>	mStrMat;
typedef mStrMat::iterator				mStrMat_it;

} // End namespace ONETOOL

#include "utils/matrices/blk.h"
#include "utils/matrices/diag.h"
#include "utils/matrices/std.h"
#include "utils/matrices/sym.h"
#include "utils/matrices/idt.h"
#include "utils/matrices/zero.h"
#include "utils/matrices/sparse.h"

namespace ONETOOL {

/*
 *
 * Matrix initialization
 * 
 */

/**
 * sseMinit	Initiates standard matrix with specific value
 *
 * @param	Ra_m	Matrix (M) [n*p] to be initated
 * @param	N_row	# row of M
 * @param	N_col	# column of M
 * @param	R_val	Initial value
 * @return	(void)
 */
void	sseMinit(wsMat Ra_m, wsUintCst N_row, wsUintCst N_col, wsRealCst R_val);

wsMat	sseMrand(wsUint N_r, wsUint N_c);
wsMat	sseMrand01(wsUint N_r, wsUint N_c, wsReal R_prob);

/**
 * sseMset0	Zero-filling matrix, a bit faster than sseMinit with 0
 *
 * @param	Ra_m	Matrix (M) [n*p] to be zero-filled
 * @param	N_row	# row of M
 * @param	N_col	# column of M
 * @return	(void)
 */
void	sseMset0(wsMat Ra_m, wsUintCst N_row, wsUintCst N_col);

wsRealCst sseMmin(wsMat Ra_m, wsUintCst N_row, wsUintCst N_col,
	wsUint* Np_minIdxR=NULL, wsUint* Np_minIdxC=NULL);
wsRealCst sseMmax(wsMat Ra_m, wsUintCst N_row, wsUintCst N_col,
	wsUint* Np_maxIdxR=NULL, wsUint* Np_maxIdxC=NULL);

/**
 * sseSinit	Initiates symmetric matrix with specific value
 *
 * @param	Ra_m	Matrix (M) [n*p] to be initated
 * @param	N_sz	Size of sym. matrix
 * @param	R_val	Initial value
 * @return	(void)
 */
void	sseSinit(wsMat Ra_m, wsUintCst N_sz, wsRealCst R_val);

wsMat	sseMsubset(wsMat Rp_inpMatrix,
	wsUintCst N_rowStart, wsUintCst N_rowCount,
	wsUintCst N_colStart, wsUintCst N_colCount);
/* Make a subset of given matrix with given sequence of indices */
wsMat	sseMsubset(wsMat Rp_inpMatrix, wsUint N_sz, vInt& Nv_newSeq);
wsMat	sseMsubset(wsMat Rp_inpMatrix, wsUint N_sz, wsUint *Na_newSeq,
	wsUint N_newSz=0xffffffff);
wsMat	sseShuffleMatrix(wsMat Ra_mat, wsUint N_msz, wsUint N_sz,
	char B_noReplace=0);
/* Get (ret)11, 12, 22 */
wsSym	sseSplitMatrix(wsMat Rp_inpMatrix, wsUint N_sz, char *Ba_isInc,
	wsUint N_tsz, wsMat *Ra_12, wsSym *Ra_22);
/* Get (ret)22, 12 */
wsSym	sseSplitMatrixV2(wsMat Rp_inpMatrix, wsUint N_sz, char *Ba_isInc,
	wsUint N_szX1, wsMat *Ra_12);
wsSym	sseSplitMatrix2(wsMat Rp_inpMatrix, wsUint N_sz, char *Ba_isExc,
	wsUint N_tsz, wsMat *Ra_12, wsSym *Ra_22);
wsMat	sseSubremMatrix(wsReal **Ra_m, wsUint N_sz, wsUint N_st,
	wsUint N_end=0xffffffff);

wsSym	sseSsubset(wsSym Rp_inpMatrix, wsUint N_start, wsUint N_count);
wsSym	sseSsubset(wsMat Rp_inpMatrix, wsUint N_sz, vInt& Nv_newSeq);
wsSym	sseSrand(wsUint N_sz);

/*
 *
 * Matrix addition
 * 
 */

/* Std. matrix + scalar */
void sseMaC(wsMat Ra_m, wsUintCst N_r, wsUintCst N_c, wsRealCst R_v,
	wsMat Ra_t, wsUintCst N_y=0, wsUintCst N_x=0);
/* Std. matrix + Std. matrix */
void	sseMaM(wsMat Ra_m1, wsMat Ra_m2, wsMat Ra_mt, wsUint N_row,
	wsUint N_col=0xffffffff, wsReal R_multiplier=REAL_CONST(.0));
/* Do mat + mat with coordinate */
void	sseMaM(wsMat Ra_m1, wsMat Ra_m2, wsMat Ra_mt, wsUint N_y,
	wsUint N_x, wsUint N_row, wsUint N_col=0xffffffff,
	wsReal R_multiplier=REAL_CONST(.0));
/* Std. matrix - matrix from (Vector * Vector) */
void	sseMaVtV(wsMat Ra_m1, wsVec Ra_v1, wsVec Ra_v2, wsUint N_row, wsUint N_col);

/* Do sym^2 */
wsSym	sseSsq(wsSym Ra_s, wsUint N_sz);
/* Do mat^2 */
wsMat	sseMsq(wsMat Ra_m, wsUint N_row, wsUint N_col);
/* Std. matrix + Sym. matrix */
void sseMaS(wsMatCst Ra_m, wsSym Ra_s, wsMat Ra_t, wsUintCst N_sz);
/* Std. matrix - Sym. matrix */
void sseMsS(wsMatCst Ra_m, wsSymCst Ra_s, wsMat Ra_t, wsUintCst N_sz);

/* Std. matrix - Std. matrix
 * (Std. matrix - Std. matrix)*const */
void	sseMsM(wsMat Ra_m1, wsMat Ra_m2, wsMat Ra_mt,
	wsUint N_row, wsUint N_col=0xffffffff, wsReal R_multiplier=0.0);
/* Std. matrix - matrix from (Vector * Vector) */
void	sseMsVtV(wsMat Ra_m1, wsVec Ra_v1, wsVec Ra_v2, wsUint N_row, wsUint N_col);
void sseMsMsq(wsMat Ra_m1, wsMat Ra_m2, wsMat Ra_mt,
	wsUint N_row, wsUint N_col=0xffffffff);
void sseMsVsq(wsMat Ra_m1, wsReal *Ra_m2, wsMat Ra_mt,
	wsUint N_row, wsUint N_col=0xffffffff);
/* Do mat * (1-mat) */
void	sseMp1p(wsMat Ra_m1, wsMat Ra_mt, wsUint N_row, wsUint N_col=0xffffffff);
/* Do sym - mat */
void	sseSsM(wsMat Ra_m1, wsMat Ra_m2, wsMat Ra_mt,
	wsUint N_row, wsUint N_col=0xffffffff, wsReal R_multiplier=0.0);
/* Do mat - mat with coordinate */
void	sseMsM(wsMat Ra_m1, wsMat Ra_m2, wsMat Ra_mt, wsUint N_y,
	wsUint N_x, wsUint N_row, wsUint N_col, wsReal R_multiplier=0.0);
void	sseMsC(wsMatCst Ra_m1, wsRealCst R_s, wsMat Ra_mt, wsUintCst N_row,
	wsUint N_col=0xffffffff);

/* [SSE] mat(A) %*% mat(B) */
wsMat	sseMpM(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsReal **Rp_rMean=NULL);

/**
 * sseMtM	Performs M'M multiplication, return matrix will be symmetric
 *
 * @param	Ra_m	Matrix (M) [n*p] to be multiplicated
 * @param	N_r		# row of matrix
 * @param	N_c		# column of matrix
 * @param	R_mul	Value to be multiplied, 0 if no multiplier
 * @return	(SYM_t)	Result of M'M [p*p] = [p*n] %*% [n*p]
 */
wsSym	sseMtM(wsMat Ra_m, wsUint N_r, wsUint N_c, wsRealCst R_mul=W0);

/**
 * sseMtM	Performs M'N multiplication, return matrix will be standard
 *
 * @param	Ra_m1	Former matrix (M) [n*p] to be multiplicated
 * @param	N_r1	# row of M
 * @param	N_c1	# column of M
 * @param	Ra_m2	Latter matrix (N) [n*m] to be multiplicated
 * @param	N_r2	# row of N
 * @param	N_c2	# column of N
 * @return	(MAT_t)	Result of M'N [p*m] = [p*n] %*% [n*m]
 */
wsMat	sseMtM(wsMat Ra_m1, wsUint N_r1, wsUint N_c1, wsMat Ra_m2,
	wsUint N_r2, wsUint N_c2);

/**
 * sym_sseMtM	Performs M'N multiplication, but assumes the result will be
 *				symmetric
 *
 * @param	Ra_m1	Former matrix (M) [n*p] to be multiplicated
 * @param	N_r1	# row of M
 * @param	N_c1	# column of M
 * @param	Ra_m2	Latter matrix (N) [n*p] to be multiplicated
 * @param	N_r2	# row of N
 * @param	N_c2	# column of N
 * @return	(SYM_t)	Result of M'N, symmetric matrix [p*p] = [p*n] %*% [n*p]
 */
wsSym	sym_sseMtM(wsMat Ra_m1, wsUint N_r1, wsUint N_c1, wsMat Ra_m2,
	wsUint N_r2, wsUint N_c2);

/**
 * sseMtS	Performs M'S multiplication, return matrix will be standard
 *
 * @param	Ra_m1	Matrix (M) [n*p] to be multiplicated
 * @param	N_r1	# row of matrix
 * @param	N_c1	# column of matrix
 * @param	Ra_m2	Symmetric matrix (S) [n*n] to be multiplicated
 * @param	N_s2	size of S
 * @return	(MAT_t)	Result M'S [p*n] = [p*n] %*% [n*n]
 */
wsMat	sseMtS(wsMat Ra_m1, wsUint N_r1, wsUint N_c1, wsSym Ra_m2,
	wsUint N_s2);

wsVec sseVpM(wsVec Rp_m1, wsUint N_c1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsReal **Rp_rMean=NULL);
/* [SSE] mat(A) %*% mat(B), assum symmetry of result */
wsSym	sym_sseMpM(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2);
/* [TRD] mat(A) %*% mat(B) */
void	multMpM(wsMat Rp_m1, wsUint N_r1, wsUint N_c1, wsUint N_y1, wsUint N_x1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsUint N_y2, wsUint N_x2,
	wsMat Ra_T, wsUint N_rT, wsUint N_cT, wsUint N_yT, wsUint N_xT);

/* [SSE] mat(A) %*% mat(B)' */
wsMat	sseMpMt(wsMat Rp_mat1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_mat2, wsUint N_r2, wsUint N_c2, wsReal **Rp_rMean=NULL);
/* [SSE] mat(A) %*% mat(A)' */
wsSym	sseMpMt(wsMat Rp_m, wsUint N_r1, wsUint N_c1);
/* [SSE] mat(A) %*% mat(B)', assum symmetry of result */
wsSym	sym_sseMpMt(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_m2=NULL, wsUint N_r2=0, wsUint N_c2=0,
	wsReal **Rp_rMean=NULL);


void	multMpMt(wsMat Rp_m1, wsUint N_r1, wsUint N_c1, wsUint N_y1, wsUint N_x1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsUint N_y2, wsUint N_x2,
	wsMat Ra_T, wsUint N_rT, wsUint N_cT, wsUint N_yT, wsUint N_xT,
	wsReal **Rp_rMean=NULL);
wsMat	multMpMt(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsReal **Rp_rMean=NULL);
wsMat	multMS(wsMat Ra_mat, wsUint N_row, wsUint N_col, wsSym Ra_symMat,
	wsUint N_sz);
/**
 * multMP	Performs MP multiplication, return matrix will be standard
 *
 * @param	Ra_mat	Former matrix (M) [n*p] to be multiplicated
 * @param	N_row	# row of M
 * @param	N_col	# column of M
 * @param	Rp_data	Value vector of sparse matrix (P) [p*q] to be multiplicated
 * @param	N_pr	# row of P
 * @param	N_pc	# column of P
 * @param	Na_rEnd	Offset vector of row [size N_pr] end
 * @param	Na_cIdx	Index vector of column [size Na_rEnd[END]] for items
 * @return	(MAT_t)	Return of MP [n*q] = [n*p] %*% [p*q]
 */
wsMat	multMP(wsMat Ra_mat, wsUint N_row, wsUint N_col, wsReal *Rp_data,
	wsUint N_pr, wsUint N_pc, wsUint *Na_rEnd, wsUint *Na_cIdx);

wsMat	multSMt(wsSym Ra_m1, wsUint N_sz, wsMat Ra_m2, wsUint N_r2,
	wsUint N_c2);

wsMat	sseMpD(wsMat Ra_m1, wsUint N_r, wsUint N_c, wsReal *Ra_D);
void	sseMpC(wsMatCst Ra_m, wsRealCst R_s, wsMat Ra_mt, wsUintCst N_row,
	wsUint N_col=0xffffffff);
void	sseSpC(wsSymCst Ra_m, wsRealCst R_s, wsSym Ra_mt, wsUintCst N_sz);

wsSym	sseStS(wsSym Ra_s, wsUint N_sz);
wsSym	sseSkS(wsSym Ra_s1, wsUint N_s1, wsSym Ra_s2, wsUint N_s2);
void	sseSsS(wsSym Ra_m1, wsSym Ra_m2, wsSym Ra_mt, wsUint N_sz,
			wsReal R_multiplier=REAL_CONST(.0));


/* Get the trace of matrix s1*s1, where s1 is symmetric matrix */
wsReal sseTrSS(wsSym Ra_s1, wsUint N_s1);
wsReal sseTrSS(wsSym Ra_s1, wsUint N_s1, wsSym Ra_s2);

/* Do vt = -v1 or v1 = -v1 */
void sseVneg(wsReal *Ra_v1, wsUint N_sz, wsReal *Ra_vt=NULL, wsReal R_add=WISARD_NAN);
/* Do vt = -v1*c or v1 = -v1*c */
void sseVnegMul(wsReal R_mul, wsReal *Ra_v1, wsUint N_sz, wsReal *Ra_vt=NULL);
/**

 * sseMBt
 */
// Can calculate X[p*m]  B[mk*mk]
//               X[p*mk] B[mk*mk]
wsMat sseMpB(wsMat Rp_mat1, wsUint N_r1, wsUint N_c1,
			 wsMat Rp_bMat, wsUint N_s2, wsUint N_blkRep);
wsMat sseMpBt(wsMat Rp_mat1, wsUint N_r1, wsUint N_c1,
			  wsMat Rp_bMat, wsUint N_s2, wsUint N_blkRep);
wsMat sseSpD(wsSym Ra_s, wsUint N_sz, wsDiag Ra_d);
wsSym sseSpDpSt(wsSym Ra_s, wsUint N_sz, wsDiag Ra_d);
wsSym sseSpB(wsSym Ra_s, wsUint N_sz, wsMat Ra_B, wsUint N_szBlk, wsUint N_rep);

/**
 * sseBM
 */
wsMat sseBpM(wsMat Rp_bMat, wsUint N_szBlk, wsUint N_rep,
			 wsMat Rp_mat, wsUint N_r, wsUint N_c);
wsMat sseBpMt(wsMat Rp_bMat, wsUint N_szBlk, wsUint N_rep,
			  wsMat Rp_mat, wsUint N_r, wsUint N_c);

wsMat		sseSym2Mat(wsSym Ra_s, wsUint N_sz);

/* Row-mean */
wsReal*		sseMmeanR(wsMat Ra_mat, wsUint N_row, wsUint N_col=0xffffffff);
wsReal*		sseMmeanRavail(wsMat Ra_mat, wsUint N_row, wsUint N_col=0xffffffff);
wsReal*		sseMmeanC(wsMat Ra_mat, wsUint N_row, wsUint N_col=0xffffffff);
wsReal*		sseMmeanCavail(wsMat Ra_mat, wsUint N_row, wsUint N_col=0xffffffff);
wsReal*		sseSmeanR(wsMat Ra_mat, wsUint N_row);

/* Row-sum */
wsReal*		sseMsumR(wsMat Ra_mat, wsUint N_row, wsUint N_col=0xffffffff,
	wsReal *Rp_sumAll=NULL);
wsReal*		sseSsumR(wsMat Ra_mat, wsUint N_sz, wsReal *Rp_sumAll=NULL);
wsReal*		sseBsumR(wsMat Ra_mat, wsUint N_szBlk, wsUint N_sz,
	wsReal *Rp_sumAll=NULL);

/* Row-sqsum */
wsReal*		sseMssR(wsMat Ra_mat, wsUint N_row, wsUint N_col=0xffffffff,
	wsReal *Rp_sumAll=NULL);

/* Col-wise sd & mean */
void		sseMsdmeanC(wsMat Ra_m, wsUintCst N_r, wsUintCst N_c, wsVec* Rp_mean, wsVec* Rp_sd);

/* Col-sum */
wsVec		sseMsumC(wsMat Ra_mat, wsUint N_row, wsUint N_col,
				wsReal *Rp_sumAll=NULL);
wsReal*		sseMasumC(wsMat Ra_mat, wsUint N_row, wsUint N_col,
				wsReal *Rp_sumAll=NULL);
wsReal*		sseBsumC(wsMat Ra_mat, wsUint N_szBlk, wsUint N_sz,
				wsReal *Rp_sumAll=NULL);

/* Entire-sum */
wsRealCst	sseMsum(wsReal **Ra_m, wsUint N_w, wsUint N_h=0xffffffff);
wsRealCst	sseMsqsum(wsReal **Ra_m, wsUint N_w, wsUint N_h=0xffffffff);
wsRealCst	sseBsum(wsReal **Ra_m, wsUint N_szBlk, wsUint N_rep);
wsRealCst	sseSsum(wsReal **Ra_m, wsUint N_sz);

/* Family-wise eigenvalue decomposition */
wsReal*		famEIGENDECOMPOSITION(wsSym Ra_phi, wsUint N_sz,
	wsUint *Na_fam, wsMat *Rp_eVec=NULL, char B_noTranspose=0);

/* Matrix walking */
typedef wsSym (*SYMMATWALKPROC)(wsSym S, wsUint N);
typedef wsMat (*STDMATWALKPROC)(wsMat M, wsUint R, wsUint C);

/* Diagonal of matrix */
wsVec sseSdiag(wsSym Ra_mat, wsUint N_sz);
wsVec sseMdiag(wsMat Ra_mat, wsUint N_r, wsUint N_c);

/**
 * sseMMM is SSE-enabled form of multMMM()
 */
wsMat sseMMM(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_m2,  wsUint N_r2, wsUint N_c2,
	wsMat Rp_m3,  wsUint N_r3, wsUint N_c3);
wsMat sseMSM(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsSym Rp_m2,  wsUint N_r2, wsUint N_c2,
	wsMat Rp_m3,  wsUint N_r3, wsUint N_c3);
wsMat sseMDM(wsMat Ra_m1, wsUint N_r1, wsUint N_c1,
	wsDiag Ra_diag, wsMat Ra_m3,  wsUint N_r3, wsUint N_c3);

/**
 * sseMMMt performs multiplication of three matrices, but more memory-efficient
 *			than two invocation of multMMt(), and third matrix should be transposed form
 *			when the arguments after N_mat2Col is omitted, third matrix will be assumed
 *			to equal with Rp_mat1
 *
 * @param     Rp_mat1		A pointer to first matrix
 * @param     N_mat1Row		Row size of first matrix
 * @param     N_mat1Col		Column size of first matrix, this value should equal to N_mat2Row
 * @param     Rp_mat2		A pointer to second matrix
 * @param     N_mat2Row		Row size of second matrix, this value should equal to N_mat1Col
 * @param     N_mat2Col		Column size of second matrix, this value should equal to N_mat3Col
 * @param     Rp_mat3		A pointer to third matrix, if omitted, Rp_mat1 will be assumed
 * @param     N_mat3Row		Row size of third matrix
 * @param     N_mat3Col		Column size of third matrix, this value should equal to N_mat2Col
 * @return    (wsReal**)	A pointer to matrix have its size to (N_mat1Row*N_mat3Row)
 *							Contains result of multiplication
 */
wsMat sseMMMt(wsMat Ra_m1, wsUint N_r1, wsUint N_c1,
	wsMat Ra_m2, wsUint N_r2, wsUint N_c2,
	wsMat Ra_m3, wsUint N_r3, wsUint N_c3);
wsSym sseMMMt(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_m2,  wsUint N_r2, wsUint N_c2);
wsMat sseMDMt(wsMat Ra_m1, wsUint N_r1, wsUint N_c1,
	wsDiag Ra_d, wsMat Ra_m3, wsUint N_r3, wsUint N_c3);
wsSym sseMDMt(wsReal **Rp_m1, wsUint N_r1, wsUint N_c1,
	wsReal *Rp_d);
wsVec	multMV(wsMat Rp_mat, wsUintCst N_row, wsUintCst N_col, wsVecCst Rp_v);
wsReal	multVMV(wsVecCst Ra_v, wsMatCst Ra_m, wsUintCst N_len, wsVecCst Ra_v2=NULL);
wsMat	multMtM(wsMatCst Rp_mat, wsUint N_row, wsUint N_col);
wsMat	multMMt(wsMatCst Rp_mat, wsUint N_row, wsUint N_col);
wsReal	sseVpMpV(wsVec Ra_v, wsMat Ra_m, wsUint N_len,
	wsVecCst Ra_v2=NULL, wsReal *Ra_mask=NULL);
wsReal	sseVpSpV(wsVec Ra_V, wsSym Ra_S, wsUintCst N_sz, wsVec Ra_W=NULL,
	wsReal *Ra_mask=NULL);

wsReal sseVVmiss(wsRealCst *Ra_v1, wsUint N_sz, wsRealCst *Ra_v2=NULL);
wsReal sseVVsel(wsReal *Ra_mask, wsRealCst *Ra_v1, wsUint N_sz,
	wsRealCst *Ra_v2=NULL);

wsRealCst	diagSum(wsReal **Ra_m, wsUint N_sz);
wsSym	diagMatrix(wsReal R_val, wsUint N_size);

wsReal diagSum2matrix(wsMat Ra_m1, wsUint N_r1, wsUint N_c1,
	wsMat Ra_m2, wsUint N_r2, wsUint N_c2);

wsReal** invSymMat(wsReal **Ra_mat, wsUint N_row, wsUint B_noParallel=0);

wsReal** SO_invMatrix(wsReal **Ra_mat, wsUint N_row,
					  wsUint N_start, wsUint N_end=0xffffffff, wsReal *Rp_det=NULL, char B_isDebug=0);
wsReal** SO_invMatrix(wsReal **Ra_mat, wsUint N_row, wsReal *Rp_det);
wsReal** dbgSO_invMatrix(wsReal **Ra_matrix, wsUint N_row);
wsReal** SO_invMatrix(wsReal **Ra_mat, wsUint N_row);
wsReal** SinvSymMat(wsReal **Ra_mat, wsUint N_row, wsUint B_noParallel=0);

/* Do Kronecker product */
wsMat	krnkMM(wsMat Rp_m1, wsUint N_r1, wsUint N_c1, wsMat Rp_m2,
	wsUint N_r2, wsUint N_c2);
wsSym	krnkSS(wsSym Rp_m1, wsUint N_s1, wsSym Rp_m2, wsUint N_s2);
wsSym	krnkSI(wsSym Rp_m1, wsUint N_s1, wsUint N_s2);

} // End namespace ONETOOL

#endif
