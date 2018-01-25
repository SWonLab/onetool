#pragma once
#ifndef __WISARD_VECTOR_H__
#define __WISARD_VECTOR_H__

#include "global/common.h"
#include "utils/util.h"
#define ASSERT_SIZE(a,b,c) if ((a) != (b)) \
	halt_fmt(WISARD_SYST_INVL_DIM_MSG, (c), (a), (b))

namespace ONETOOL {

void exportVector(const char *S_ext, vInt& Nv_vec);
void exportVector(const char *S_ext, vector<wsFloat>& Nv_vec);
void exportVector(const char *S_ext, vStr& Sv_vec, wsUint N_sz);
void exportVector(const char *S_ext, wsVecCst Rp_vec, wsUint N_sz, wsUint N_prec=0);
void exportVector(const char *S_ext, const wsFvec Rp_vec, wsUint N_sz, wsUint N_prec=0);
void exportUnitVector(const char *S_ext, wsUint N_sz);
void exportValVector(const char *S_ext, wsReal R_val, wsUint N_sz, wsUint N_prec = 0);

class cStdMatrix;
class cBlkMatrix;
class cIdtMatrix;
class cDiagMatrix;
class cSymMatrix;
class cSpsMatrix;

wsVec	vectorRepeat(wsVecCst Ra_vec, vInt& V_cntRep, wsUint N_sz);
wsVec	vectorRepeat(wsVecCst Ra_vec, mvStr& V_cntRep, wsUint N_sz);
wsVec	vectorCummax(wsVecCst Ra_vec, wsUint N_sz, xRealSort* Xa_rs=NULL);
wsVec	vectorCummax2(wsVecCst Ra_vec, wsUint N_sz, xRealSort* Xa_rs=NULL);
wsReal  sseVse(wsVec Ra_t, wsUint N_sz, wsReal *Rp_mean=NULL);
void    sseVstdize(wsVec Ra_t, wsUint N_sz);

class cVector
{
	static char B_report;
	char	B_dontDealloc;
protected:
	wsUint N_sz;
	wsReal *Ra_buf;
public:
	static void report(char V) { B_report = V; }
	cVector();
	cVector(cVector &&V_noname);
	//cVector(cVector &C_vec);
	cVector(cVector &C_vec, wsUint N_rep);
	cVector(wsReal *Ra_val, wsUint N_inpSz, char B_inpDDealloc=0);
	cVector(wsReal R_val, wsUint N_inpSz);
	cVector(wsUint N_inpSz, char B_rand=0);
	cVector(xRealSort *Xa_rs, wsUint N_inpSz);
	~cVector();

	void		init(wsUint N_inpSz, wsReal *Ra_val=NULL, wsReal *Rp_val=NULL,
		char B_rand=0, char B_inpDDealloc=0);
	wsReal*		get() const;
	/**
	 * cVector::size	Returns the size of vector, if the vector is not
	 *					initialized, returns 0
	 *
	 * @return    (wsUintCst)	The size of vector
	 */
	wsUintCst		size() const;
	/**
	 * cVector::sum	Computes sum of the vector and squared sum optionally.
	 *
	 * @param     Rp_sqSum	(optional) If this parameter given, squared sum
	 *						will be stored to here
	 * @return    (REAL_c)	Sum of the vector
	 */
	wsRealCst		sum(wsReal *Rp_sqSum=NULL);
	wsRealCst		asum(wsReal *Rp_sqSum=NULL);
	wsRealCst		aasum(cVector& Vsub);
	wsRealCst		asum(cVector& Vsub);
	wsRealCst		ssum(cVector& Vsub);
	/**
	 * cVector::sum	Computes sum with the other vector == t(this) %*% V
	 *
	 * @param     Ra_vec	A pointer to the vector to be included to sum
	 * @param     N_sz		# elements of Ra_vec
	 * @return    (REAL_c)	Result of t(this) %*% V
	 */
	wsRealCst		sum(wsVec Ra_vec, wsUint N_sz);
	/**
	 * cVector::sum	Computes sum with the other vector == t(this) %*% V
	 *
	 * @param     cVector&	Vector instance to be included to sum
	 * @return    (REAL_c)	Result of t(this) %*% V
	 */
	wsRealCst		sum(cVector &C_vec);
	/**
	 * cVector::mean	Computes the mean(average) of the vector
	 *
	 * @return    (REAL_c)	Mean of the vector
	 */
	wsRealCst		mean();
	/**
	 * cVector::rem	Truncates the contents of vector and empty it
	 *
	 * @return    (void)
	 */
	void		rem();
	/**
	 * cVector::set0	Initiates the contents of vector to 0, note that
	 *					this method assumes the buffer of vector is already
	 *					allocated.
	 *
	 * @return    (void)
	 */
	void		set0();
	/**
	 * cVector::setDontDealloc	Protects the buffer from uninitialization
	 *
	 * @return    (void)
	 */
	void		setDontDealloc() { B_dontDealloc = 1; }
	cVector&	inv();		/* Compute self-inversion */
	cVector		p1p();		/* Compute v(1-v) */
	cVector		p21p2();	/* Compute v^2 * (1-v)^2 */
	cStdMatrix	toMt();
	char		isDdealloc() { return B_dontDealloc; }

	/**
	 * cVector::colP  Performs column-wise product, vector size must be
	 *                identical to the row size of input matrix
	 *
	 *
	 */
	cStdMatrix	colP(cStdMatrix& M);

	/* Calculate quadratic form
	 * this %*% mat %*% t(vec) */
	wsRealCst		qf(wsMat Ra_mat, wsUint N_sz, wsVec Ra_vec);
	/* Calculate quadratic form
	 * this %*% mat %*% t(this) */
	wsReal		qf(wsMat Ra_mat, wsUint N_sz);
	wsReal		qf(cBlkMatrix &B);
	wsReal		qf(cBlkMatrix &B, cVector &V);
	wsReal		qf(cStdMatrix &M);
	wsReal		qf(cStdMatrix &M, cVector &V);
	wsReal		qf(cIdtMatrix &I);
	wsReal		qf(cIdtMatrix &I, cVector &V);
	wsReal		qf(cDiagMatrix &D);
	wsReal		qf(cDiagMatrix &D, cVector &V);
	wsReal		qf(cSymMatrix &S);
	wsReal		qf(cSymMatrix &S, cVector &V);
	wsReal		qf(cSpsMatrix &P);
	wsReal		qf(cSpsMatrix &P, cVector &V);

	cVector&	operator=(const cVector &V);
	cVector&	operator=(cVector &&V);

	/* Self subtraction */
	void		operator-=(wsReal R_val);
	void		operator-=(cVector &C_vec);
	/* Self addition */
	cVector&	operator+=(wsReal R_val);
	cVector&	operator+=(wsUint N_val);
	cVector&	operator+=(cVector &C_vec);
	/* Self production */
	void		operator*=(wsRealCst R_val);
	void		operator*=(const cVector &C_vec);
	/* Self division */
	void		operator/=(wsReal R_val);
	void		operator/=(cVector &C_vec);
	/* Subtract and return */
	cVector		operator-(wsRealCst R_val);
	cVector		operator-(const cVector &C_vec);
	/* Addition and return */
	cVector		operator+(wsRealCst R_val);
	cVector		operator+(const cVector &C_vec);
	/* Production and return */
	cVector		operator*(wsRealCst R_val);
	cVector		operator*(const cVector &C_vec);
	cVector		operator*(const cStdMatrix &C_mat);
	cVector		operator*(const cSymMatrix &S);
	cVector		Mt(const cStdMatrix &M);
	cVector&	operator*=(cDiagMatrix &D);
	cVector&	operator*=(cIdtMatrix &I);
	cVector&	operator*=(cSpsMatrix &P);
	cVector		operator*(const cDiagMatrix &D);
	cVector		operator*(const cIdtMatrix &I);
	cVector		operator*(const cSpsMatrix &P);
	cVector		sq();
	cVector		sq(cVector& V);
	cVector		sqrt();
	void		selfSqrt(char B_abs=0);
	cVector		cb();
	void		mul(wsRealCst R_multipier, cVector &V_dest);
	void		mul(wsRealCst R_multipier, wsRealCst R_adder, cVector &V_dest);
	void		mul(wsRealCst R_multipier, wsRealCst R_adder);

	/* Division and return */
	cVector		operator/(wsRealCst R_val);
	cVector		operator/(const cVector &C_vec);

	/* Special operations */
	void		file(wsStrCst S_ext);
	void		fmtfile(wsStrCst S_fmtExt, ...);
	cVector&	sub(cVector &V_s, cVector &V_t);
	// Squared sum of the matrix
	wsRealCst	ss();
	cVector		clone(wsUint N_len=0);
	cVector		subset(wsUint N_start, wsUint N_sz);
	cSymMatrix	tV(wsReal R_prod=WISARD_NAN);
	cStdMatrix	tV(cVector &V, wsReal R_prod=WISARD_NAN);
	cVector&	neg(wsReal R_add=WISARD_NAN);
	cVector		reorder(wsUint *Na_seq);
	cVector		krnk(cVector &V);
	cStdMatrix	krnkV(cStdMatrix &M);
	// min(scalar, v)
	void		minSelf(wsRealCst R_val);
	// max(scalar, v)
	void		maxSelf(wsRealCst R_val);
	// log(1+exp(v))
	cVector		log1exp();
	// logit = exp(v)/(exp(v)+1)
	cVector		logit();
	// divide
	cVector		divide(wsUint* Na_idx, wsUint N_sz, cVector* Vp_rest=NULL);
	// replace
	void		replace(wsUint* Na_idx, wsUint N_sz, wsVec Ra_rep, char B_inv=0);
	// exp(v)
	cVector		exp();
	// V * v1 * v2
	wsRealCst	VV(wsReal *Ra_1, wsReal *Ra_2, wsUint N_sz);
	// V^2 * v
	wsRealCst	V2(wsReal *Ra_v, wsUint N_sz);
	// log(1+
	// return 1 if at least one is negative
	char		isAnyNeg();
	void		set(wsRealCst R_v);
	// weighted mean vec * Ra_w * N_wsz / sum(Ra_w)
	wsRealCst		wMean(wsReal *Ra_w, wsUint N_wsz);
	// Do self pinv
	void		pinv();
	// Do v - this = this
	cVector&	subFrom(cVector &V);
	// return permuted this
	cVector		permute();
	cVector		permute(wsUint* Na_seq);
};

class cSpsVector : public cVector
{
	wsUint		N_elem;
	wsUint*		Na_idx;
public:
	cSpsVector(wsUint N_inpSz, wsReal R_inpVal);
	cSpsVector(wsVec Ra_val, wsUint N_inpSz, wsUint N_inpElem, wsUint* Na_inpIdx);
	cSpsVector(cSpsVector &&V_noname);
	void		init(wsUint N_inpSz, wsReal *Ra_val=NULL, wsReal *Rp_val=NULL,
		char B_rand=0, char B_inpDDealloc=0);

	wsReal*		get();

	wsRealCst		sum(wsReal *Rp_sqSum=NULL); //impld
	wsRealCst		sum(wsReal *Ra_vec, wsUint N_sz); //impld
	wsRealCst		sum(cVector &V); //impld
	void		set0(); //impld
	cVector&	inv();
	cSpsVector	p1p(); /* Compute v(1-v) */ //impld
	cStdMatrix	toMt();

	/* Calculate quadratic form
	 * this %*% mat %*% t(vec) */
	wsReal		qf(wsMat Ra_mat, wsUint N_sz, wsReal *Ra_vec);
	/* Calculate quadratic form
	 * this %*% mat %*% t(this) */
	wsReal		qf(wsMat Ra_mat, wsUint N_sz);
	wsReal		qf(cBlkMatrix &B);
	wsReal		qf(cBlkMatrix &B, cVector &V);
	wsReal		qf(cStdMatrix &M);
	wsReal		qf(cStdMatrix &M, cVector &V);
	wsReal		qf(cIdtMatrix &I);
	wsReal		qf(cIdtMatrix &I, cVector &V);
	wsReal		qf(cDiagMatrix &D);
	wsReal		qf(cDiagMatrix &D, cVector &V);
	wsReal		qf(cSymMatrix &S);
	wsReal		qf(cSymMatrix &S, cVector &V);
	wsReal		qf(cSpsMatrix &S);
	wsReal		qf(cSpsMatrix &S, cVector &V);
	cVector&	operator=(const cVector &V);
	cVector&	operator=(cVector &&V);

	/* Self subtraction */
	void		operator-=(wsReal R_val);
	void		operator-=(cVector &C_vec);
	/* Self addition */
	cVector&	operator+=(wsReal R_val);
	cVector&	operator+=(wsUint N_val);
	cVector&	operator+=(cVector &C_vec);
	/* Self production */
	void		operator*=(wsRealCst R_val);
	void		operator*=(const cVector &C_vec);
	/* Self division */
	void		operator/=(wsReal R_val);
	void		operator/=(cVector &C_vec);
	/* Subtract and return */
	cSpsVector	operator-(wsRealCst R_inpVal);
	cVector		operator-(const cVector &V);
	/* Addition and return */
	cVector		operator+(wsRealCst R_inpVal);
	cVector		operator+(const cVector &V);
	/* Production and return */
	cVector		operator*(wsRealCst R_inpVal);
	cVector		operator*(cVector &V);
	cVector		operator*(cStdMatrix &M);
	cVector		operator*(cSymMatrix &S);
	cVector		Mt(cStdMatrix &M);
	cVector&	operator*=(cDiagMatrix &D);
	cVector&	operator*=(cIdtMatrix &I);
	cVector		operator*(cDiagMatrix &D);
	cVector		operator*(const cIdtMatrix &I);
	cVector		sq();
	cVector		sqrt();
	cVector		cb();
	/* Division and return */
	cSpsVector	operator/(wsRealCst R_inpVal);
	cVector		operator/(cVector &V);

	/* Special operations */
	void		file(wsStrCst S_ext);
	cVector&	sub(cVector &V_s, cVector &V_t);
	wsRealCst		ss();
	cSpsVector	clone();
	cSpsVector	subset(wsUint N_start, wsUint N_sz);
	cSymMatrix	tV(wsReal R_prod=WISARD_NAN);
	cStdMatrix	tV(cVector &V, wsReal R_prod=WISARD_NAN);
	cVector&	neg(wsReal R_add=WISARD_NAN);
};

class cValVector : public cVector
{
	wsReal		R_val;
public:
	cValVector(wsUint N_inpSz, wsReal R_inpVal);
	cValVector(cValVector &&V_noname);
	~cValVector();
	void		init(wsUint N_inpSz, wsReal *Ra_val=NULL, wsReal *Rp_val=NULL,
		char B_rand=0, char B_inpDDealloc=0);

	wsReal*		get();
	wsReal		getV() { return R_val; }

	wsRealCst		sum(wsReal *Rp_sqSum=NULL);
	wsRealCst		sum(wsReal *Ra_vec, wsUint N_sz);
	wsRealCst		sum(cVector &V);
	void		set0();
	cVector&	inv();
	cValVector	p1p(); /* Compute v(1-v) */
	cStdMatrix	toMt();

	/* Calculate quadratic form
	 * this %*% mat %*% t(vec) */
	wsReal		qf(wsMat Ra_mat, wsUint N_sz, wsReal *Ra_vec);
	/* Calculate quadratic form
	 * this %*% mat %*% t(this) */
	wsReal		qf(wsMat Ra_mat, wsUint N_sz);
	wsReal		qf(cBlkMatrix &B);
	wsReal		qf(cBlkMatrix &B, cVector &V);
	wsReal		qf(cStdMatrix &M);
	wsReal		qf(cStdMatrix &M, cVector &V);
	wsReal		qf(cIdtMatrix &I);
	wsReal		qf(cIdtMatrix &I, cVector &V);
	wsReal		qf(cDiagMatrix &D);
	wsReal		qf(cDiagMatrix &D, cVector &V);
	wsReal		qf(cSymMatrix &S);
	wsReal		qf(cSymMatrix &S, cVector &V);
	wsReal		qf(cSpsMatrix &S);
	wsReal		qf(cSpsMatrix &S, cVector &V);
	cVector&	operator=(const cVector &V);
	cVector&	operator=(cVector &&V);

	/* Self subtraction */
	void		operator-=(wsReal R_val);
	void		operator-=(cVector &C_vec);
	/* Self addition */
	cVector&	operator+=(wsReal R_val);
	cVector&	operator+=(wsUint N_val);
	cVector&	operator+=(cVector &C_vec);
	/* Self production */
	void		operator*=(wsRealCst R_val);
	void		operator*=(const cVector &C_vec);
	/* Self division */
	void		operator/=(wsReal R_val);
	void		operator/=(cVector &C_vec);
	/* Subtract and return */
	cValVector	operator-(wsRealCst R_inpVal);
	cVector		operator-(const cVector &V);
	/* Addition and return */
	cVector		operator+(wsRealCst R_inpVal);
	cVector		operator+(const cVector &V);
	/* Production and return */
	cVector		operator*(wsRealCst R_inpVal);
	cVector		operator*(cVector &V);
	cVector		operator*(cStdMatrix &M);
	cVector		operator*(cSymMatrix &S);
	cVector		Mt(cStdMatrix &M);
	cVector&	operator*=(cDiagMatrix &D);
	cVector&	operator*=(cIdtMatrix &I);
	cVector		operator*(cDiagMatrix &D);
	cVector		operator*(const cIdtMatrix &I);
	cVector		sq();
	cVector		sqrt();
	cVector		cb();
	/* Division and return */
	cValVector	operator/(wsRealCst R_inpVal);
	cVector		operator/(cVector &V);

	/* Special operations */
	void		file(wsStrCst S_ext);
	cVector&	sub(cVector &V_s, cVector &V_t);
	wsRealCst		ss();
	cValVector	clone();
	cValVector	subset(wsUint N_start, wsUint N_sz);
	cSymMatrix	tV(wsReal R_prod=WISARD_NAN);
	cStdMatrix	tV(cVector &V, wsReal R_prod=WISARD_NAN);
	cVector&	neg(wsReal R_add=WISARD_NAN);
};

class cUnitVector : public cVector
{
public:
	cUnitVector(wsUint N_inpSz);
	cUnitVector(cUnitVector &&V_noname);
	void		init(wsUint N_inpSz, wsReal *Ra_val=NULL, wsReal *Rp_val=NULL,
		char B_rand=0, char B_inpDDealloc=0);
	wsReal*		get() const;
	wsRealCst		sum(wsReal *Rp_sqSum=NULL);
	wsRealCst		sum(wsReal *Ra_vec, wsUint N_sz);
	wsRealCst		sum(cVector &V);
	void		set0();
	cVector&	inv();
	cVector		p1p(); /* Compute v(1-v) */
	cStdMatrix	toMt();

	/* Calculate quadratic form
	 * this %*% mat %*% t(vec) */
	wsReal		qf(wsMat Ra_mat, wsUint N_sz, wsReal *Ra_vec);
	/* Calculate quadratic form
	 * this %*% mat %*% t(this) */
	wsReal		qf(wsMat Ra_mat, wsUint N_sz);
	wsReal		qf(cBlkMatrix &B);
	wsReal		qf(cBlkMatrix &B, cVector &V);
	wsReal		qf(cStdMatrix &M);
	wsReal		qf(cStdMatrix &M, cVector &V);
	wsReal		qf(cIdtMatrix &I);
	wsReal		qf(cIdtMatrix &I, cVector &V);
	wsReal		qf(cDiagMatrix &D);
	wsReal		qf(cDiagMatrix &D, cVector &V);
	wsReal		qf(cSymMatrix &S);
	wsReal		qf(cSymMatrix &S, cVector &V);
	wsReal		qf(cSpsMatrix &S);
	wsReal		qf(cSpsMatrix &S, cVector &V);
	cVector&	operator=(const cVector &V);
	cVector&	operator=(cVector &&V);

	/* Self subtraction */
	void		operator-=(wsReal R_val);
	void		operator-=(cVector &C_vec);
	/* Self addition */
	cVector&	operator+=(wsReal R_val);
	cVector&	operator+=(wsUint N_val);
	cVector&	operator+=(cVector &C_vec);
	/* Self production */
	void		operator*=(wsRealCst R_val);
	void		operator*=(const cVector &C_vec);
	/* Self division */
	void		operator/=(wsReal R_val);
	void		operator/=(cVector &C_vec);
	/* Subtract and return */
	cValVector	operator-(wsRealCst R_val);
	cVector		operator-(const cVector &C_vec);
	/* Addition and return */
	cVector		operator+(wsRealCst R_val);
	cVector		operator+(const cVector &C_vec);
	/* Production and return */
	cVector		operator*(wsRealCst R_val);
	cVector		operator*(cVector &C_vec);
	cVector		operator*(cStdMatrix &M);
	cVector		operator*(cSymMatrix &S);
	cVector		Mt(cStdMatrix &M);
	cVector&	operator*=(cDiagMatrix &D);
	cVector&	operator*=(cIdtMatrix &I);
	cVector		operator*(cDiagMatrix &D);
	cVector		operator*(const cIdtMatrix &I);
	cVector		sq();
	cVector		sqrt();
	cVector		cb();
	/* Division and return */
	cValVector	operator/(wsRealCst R_val);
	cVector		operator/(cVector &C_vec);

	/* Special operations */
	void		file(wsStrCst S_ext);
	cVector&	sub(cVector &V_s, cVector &V_t);
	wsRealCst		ss();
	cVector		clone();
	cVector		subset(wsUint N_start, wsUint N_sz);
	cSymMatrix	tV(wsReal R_prod=WISARD_NAN);
	cStdMatrix	tV(cVector &V, wsReal R_prod=WISARD_NAN);
	cVector&	neg(wsReal R_add=WISARD_NAN);
};

class cMask : public cVector
{
public:
	cMask();
	cMask(cMask &C_ref, wsUint N_rep=1);
	cMask(wsUint N_sz, char B_on);
	~cMask();
	void	init(wsUint N_sz, char B_on);
	wsRealCst	vv(cVector &V, cVector &W);
	wsRealCst	vv(cVector &V);
	wsRealCst	qf(cVector &V, cStdMatrix &M);
	wsRealCst	qf(cVector &V, cStdMatrix &M, cVector &W);
	wsRealCst	qf(cVector &V, cSymMatrix &M);
	wsRealCst	qf(cVector &V, cSymMatrix &M, cVector &W);
	wsRealCst	qf(cVector &V, cIdtMatrix &M);
	wsRealCst	qf(cVector &V, cIdtMatrix &M, cVector &W);

	/* Special operations */
	void	file(wsStrCst S_ext);
};

/*
 *
 * Vector initialization
 * 
 */

/**
 * sseVinit	Initiates vector with specific value
 *
 * @param	Ra_v	Vector (V) [n*1] to be initiated
 * @param	N_sz	size of V
 * @param	R_val	Initial value
 * @return	(void)
 */
void	sseVinit(wsVec Ra_v, wsUintCst N_sz, wsRealCst R_val);

/**
 * sseVset0	Zero-filling vector
 *
 * @param	Ra_v	Vector (V) [n*1] to be zero-filled
 * @param	N_sz	size of V
 * @return	(void)
 */
void	sseVset0(wsVec Ra_v, wsUintCst N_sz);
wsVec	sseVrep(wsVecCst Ra_v, wsUintCst N_sz, wsUintCst N_times);

/*
 *
 * Vector scalar return
 *
 */

/* Vector variance, var(v) and optionally mean(v) */
wsReal	sseVvar(wsVecCst V, wsUintCst N, wsRealPtr Rp_mean=NULL);
/* Vector sum, sum(v) and optionally sum(v^2) */
wsReal	sseVsum(wsVecCst Ra_v, wsUintCst N_sz, wsRealPtr Rp_sqSum=NULL);
/* Vector squared sum, sum(v^2) or sum((v+c)^2) */
wsReal	sseVsqsum(wsVecCst Ra_v, wsUintCst N_sz, wsRealCst R_adder=WISARD_NAN);
/* Vector sum of absolute, sum(|v|) and optionally sum(v^2) */
wsReal	sseVasum(wsVecCst Ra_v, wsUintCst N_sz, wsRealPtr Rp_sqSum=NULL);
/* Vector sum of absolute log, sum(log(|v|)) */
wsReal	sseVlogasum(wsVecCst Ra_v, wsUintCst N_sz);
/* Vector sum of log, sum(log(v)) */
wsReal	sseVlogsum(wsVecCst Ra_v, wsUintCst N_sz);
/* Vector production, prod(v) */
wsRealCst	sseVprod(wsVecCst Ra_v, wsUintCst N_sz);
/* Vector sum of addition, sum(v1+v2) */
wsReal	sseVaVsum(wsVecCst Ra_v1, wsVecCst Ra_v2, wsUintCst N_sz,
	wsRealCst R_mul=WISARD_NAN);
/* Vector sum of subtraction, sum(v1-v2) */
wsReal sseVsVsum(wsVecCst Ra_v1, wsVecCst Ra_v2, wsUintCst N_sz, char B_sq=0,
	wsRealCst R_mul=WISARD_NAN);
/* Vector sum of scalar subtraction, sum(v1-c) */
wsReal	sseVsCsum(wsVecCst Ra_v1, wsRealCst R_v, wsUintCst N_sz, char B_sq=0);
/* Return Sum(i=1:N_sz) Ra_v1[i]*Ra_v2[i] or Ra_v1[i]*Ra_v1[i] */
wsReal	sseVV(wsVecCst Ra_v1, wsUintCst N_sz, wsVec Ra_v2=NULL,
	wsRealCst R_mul=WISARD_NAN);
/* Return Sum(i=1:N_sz) Ra_v1[i]*Ra_v2[i]*Ra_v3[i] */
wsReal	sseVVV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVecCst Ra_v3, wsUintCst N_sz);

/*
 *
 * Vector power
 *
 */

/* Vector square root, sqrt(v) */
void	sseVsqrt(wsVec Ra_V, wsUintCst N_sz, wsVec Ra_tgt=NULL);
/* Vector square root of absolute, sqrt(|v|) */
void	sseVabsqrt(wsVec Ra_V, wsUintCst N_sz, wsVec Ra_tgt=NULL);
void	sseVsquare(wsVecCst Ra_V, wsUintCst N_sz, wsVec Ra_vt, wsRealCst R_mul=WISARD_NAN);
void	sseVsquare(wsVecCst Ra_V, wsUintCst N_sz, wsVec Ra_vt, wsVec Ra_v2);
void	sseVcube(wsVecCst Ra_V, wsUintCst N_sz, wsVec Ra_v3);

/*
 *
 * Vector addition
 * 
 */

/* Addition of vector and scalar, v+c */
void	sseVaC(wsVecCst Ra_v, wsRealCst R_const, wsVec Ra_vt, wsUintCst N_sz);
/* Addition of vector and vector, v1+v2 or (v1+v2)*m */
void	sseVaV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
	wsRealCst R_mul=WISARD_NAN);
/* Addition of vector abs. and vector, |v1|+v2 or (|v1|+v2)*m */
void	sseAVaV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
	wsRealCst R_mul=WISARD_NAN);
/* Subtraction of vector and vector, v1-v2 or (v1-v2)*m */
void	sseVsV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
	wsRealCst R_mul=WISARD_NAN);
/* Square of "subtraction of vector and vector", (v1-v2)^2 */
void	sseVsVsq(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz);
/* Subtraction of vector and vector by scalar, v1-(v2*m) */
void	sseVsVC(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
	wsRealCst R_mul);
/* Subtraction of vector and scalar, v1-c or (v1-c)*m */
void	sseVsC(wsVecCst Ra_v1, wsRealCst R_s, wsVec Ra_vt, wsUintCst N_sz,
	wsRealCst R_mul=WISARD_NAN);

/*
 *
 * Vector production/division
 *
 */
void	sseVpV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
	wsRealCst R_mul=WISARD_NAN);
void	sseVpVsqrt(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz);
void	sseVdV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz,
	wsRealCst R_mul=WISARD_NAN);
void	sseVpC(wsVecCst Ra_v, wsRealCst R_s, wsVec Ra_vt, wsUintCst N_sz,
	wsRealPtr Rp_sqsum=NULL);
/* Do v*s + add = vt */
void	sseVpC(wsVecCst Ra_v, wsRealCst R_s, wsVec Ra_vt, wsUintCst N_sz,
	wsRealCst R_add);
/* Do v*s + vt + add = vt */
void	sseVpCadd(wsVecCst Ra_v, wsRealCst R_s, wsVec Ra_vt, wsUintCst N_sz,
	wsRealCst R_add=WISARD_NAN);
void	sseAVdV(wsVecCst Ra_v1, wsVecCst Ra_v2, wsVec Ra_vt, wsUintCst N_sz);
/**
 * sseVp1p	Performs V(1-V) or V^2 * (1-V)^2 multiplication, return will be vector
 *          although V should be [0,1], do not check it explicitly
 *
 * @param	Ra_V	Vector (V) [n*1] to be multiplicated
 * @param	N_sz	size of V
 * @param	Ra_ret	Return vector buffer, should be allocated one
 * @param   B_sq    (optional) if true, then perform V^2 * (1-V)^2
 * @return	(VEC_t)	Return of V(1-V) or V^2 * (1-V)^2
 */
wsVec	sseVp1p(wsVecCst Ra_v, wsUintCst N_sz, wsVec Ra_ret, char B_sq=0);

/**
 * sseVtV	Performs VV' or (VV')p multiplication, return will be matrix
 *
 * @param	Ra_V	Vector (V) [n*1] to be multiplicated
 * @param	N_sz	size of V
 * @param	R_mul	(optional) multiplier to be applied
 * @return	(MAT_t)	Return of VV' [n*n] = [n*1] %*% [1*n]
 */
wsMat	sseVtV(wsVecCst Ra_V, wsUintCst N_sz, wsRealCst R_mul=WISARD_NAN);

/**
 * sseVtV	Performs VW' or (VW')p multiplication, return will be matrix
 *
 * @param	Ra_v1	Former vector (V) [n*1] to be multiplicated
 * @param	N_sz	size of V and W
 * @param	Ra_v2	Former vector (W) [p*1] to be multiplicated
 * @param	R_mul	(optional) multiplier to be applied
 * @return	(MAT_t)	Return of VW' [n*n] = [n*1] %*% [1*n]
 */
wsMat	sseVtV(wsVecCst Ra_v1, wsUintCst N_sz, wsVecCst Ra_v2,
	wsRealCst R_mul=WISARD_NAN);

/**
 * sseVtV	Performs VW' or (VW')p multiplication, return will be matrix
 *
 * @param	Ra_v1	Former vector (V) [n*1] to be multiplicated
 * @param	N_s1	size of V
 * @param	Ra_v2	Former vector (W) [p*1] to be multiplicated
 * @param	N_s2	size of W
 * @param	R_mul	(optional) multiplier to be applied
 * @return	(MAT_t)	Return of VW' [n*p] = [n*1] %*% [1*p]
 */
wsMat	sseVtV(wsVecCst Ra_v1, wsUintCst N_s1, wsVecCst Ra_v2, wsUintCst N_s2,
	wsRealCst R_mul=WISARD_NAN);

void	sseVinv(wsVec Ra_v, wsUintCst N_sz, wsVec Ra_vt=NULL);

/**
 * sseVpermSelf	Permutes V with given order within itself
 *
 * @param	Ra_v	A vector [n*1] to be multiplicated
 * @param	N_sz	size of V (n)
 * @param	Na_seq	Length [n] UINT vector of offset, to be permuted
 */
void	sseVpermSelf(wsVec Ra_v, wsUintCst N_sz, wsUintCst* Na_seq);

} // End namespace ONETOOL

#endif
