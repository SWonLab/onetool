#include <math.h>
#include "utils/matrices/diag.h"
#include "utils/vector.h"
#include "global/Rconn.h"

namespace ONETOOL {

cDiagMatrix::cDiagMatrix()
{
	Ra_data	= NULL;
	N_row = N_col = 0;
}

cDiagMatrix::cDiagMatrix(cDiagMatrix &&D)
{
	init(D.row(), D.get()[0], 0, 0);
	D.setDontDealloc();
}

cDiagMatrix::cDiagMatrix(cDiagMatrix &D)
{
	init(D);
}

cDiagMatrix::cDiagMatrix(wsUint N_sz, char B_rand/*=0*/)
{
	Ra_data = NULL;
	init(N_sz, NULL, B_rand);
}

cDiagMatrix::cDiagMatrix(cVector& V)
{
	Ra_data = NULL;
	wsUint N_s = V.size();
	init(N_s, sseVectorP(N_s, V.get()));
}

cDiagMatrix::cDiagMatrix(wsUint N_sz, wsReal *Ra_V, char B_inpDDealloc/*=0*/)
{
	Ra_data = NULL;
	init(N_sz, Ra_V, 0, B_inpDDealloc);
}

cDiagMatrix::cDiagMatrix(wsReal R_V, wsUint N_sz)
{
	Ra_data = NULL;
	init(N_sz);
	sseVinit(Ra_data[0], N_sz, R_V);
// 	for (wsUint i=0 ; i<N_sz ; i++)
// 		Ra_data[0][i] = R_V;
}

cDiagMatrix::~cDiagMatrix()
{
	rem();
}

void cDiagMatrix::init(wsUint N_sz, wsReal *Ra_V/*=NULL*/, char B_rand/*=0*/,
	char B_inpDDealloc/*=0*/)
{
	if (Ra_data) sseUnmat(Ra_data, 1);

	wsMat Ra_ptr = NULL;
	if (Ra_V) {
		wsAlloc(Ra_ptr, wsReal*, 1);
		Ra_ptr[0] = Ra_V;
	} else {
		Ra_ptr = sseMatrix(1, N_sz);
	}
	__init(Ra_ptr, N_sz, N_sz);

	B_dDealloc = B_inpDDealloc;
	if (B_rand)
		for (wsUint i=0 ; i<N_sz ; i++)
			Ra_data[0][i] = (wsReal)(rand()%1000) - REAL_CONST(500.0);
}

void cDiagMatrix::init(cDiagMatrix &D)
{
	init(D.row(), sseVectorP(D.row(), D.get()[0]), 0, 0);
}

wsRealCst cDiagMatrix::sum()
{
	return sseVsum(get()[0], row());
}

wsRealCst cDiagMatrix::Lsum()
{
	return sseVlogasum(get()[0], row());
}

wsRealCst cDiagMatrix::prod()
{
	return sseVprod(get()[0], row());
}

cDiagMatrix cDiagMatrix::krnk(const cIdtMatrix &I)
{
	wsUint Ni = I.row();
	wsUint Nd = row();
	wsUint N	= Nd*Ni;
	wsReal *Ra_ret = sseVector(N);
	wsReal *S	= get()[0];

	for (wsUint i=0,j=0 ; i<Nd ; i++,j+=Ni)
		sseVinit(Ra_ret+j, Ni, S[i]);

	return cDiagMatrix(N, Ra_ret);
}

wsUintCst cDiagMatrix::nn0()
{
	/* Count the number of nonzero values */
	wsReal *Ra_S = get()[0];
	wsUint N_med = 0;
	wsUint N = row();
	wsUint R = 0;
#ifdef USE_SSE
	N_med = getMed(N);
	sse_t sse_N = sseSet(0);
	sse_t sse_0 = sseSet(0);
	sse_t sse_1 = sseSet(1);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_p = (sse_t *)(Ra_S + i);
		sse_N = sseAdd(sse_N, sseAnd(sse_1, sseNeq(*sse_p, sse_0)));
	}
	sseUsum(sse_N, R);
#endif
	for (wsUint i=N_med ; i<N ; i++)
		if (Ra_S[i] != W0)
			R++;
	return R;
}

void sseV_kV(wsReal *S, wsUint N, wsReal k, wsReal *D=NULL)
{
	wsUint N_med = 0;
	if (D) {
#ifdef USE_SSE
		N_med = getMed(N);
		sse_t sse_K = sseSet(k);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_V = (sse_t *	)(S + i);
			sse_t *sse_D = (sse_t *)(D + i);
			*sse_D = sseDiv(*sse_V, sseAdd(sse_K, *sse_V));
		}
#endif
		for (wsUint i=N_med ; i<N ; i++)
			D[i] = S[i] / (S[i] + k);
	} else {
#ifdef USE_SSE
		N_med = getMed(N);
		sse_t sse_K = sseSet(k);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_V = (sse_t *	)(S + i);
			*sse_V = sseDiv(*sse_V, sseAdd(sse_K, *sse_V));
		}
#endif
		for (wsUint i=N_med ; i<N ; i++)
			S[i] /= (S[i] + k);
	}
}

cDiagMatrix cDiagMatrix::v_kv(wsRealCst R_k)
{
	wsUint N = row();
	wsReal *R = sseVector(N);
	wsReal *S = get()[0];

	sseV_kV(S, N, R_k, R);

	return cDiagMatrix(N, R);
}

bool cDiagMatrix::sym()
{
	/* Since any diagonal matrix is symmetric, return true */
	return true;
}

cDiagMatrix& cDiagMatrix::inv()
{
	wsReal	*V		= get()[0];
	/* Prepare buffer */
	wsReal	*Ra_res	= NULL;
	sseMalloc(Ra_res, wsReal, row());
	/* Calculate inverse */
	for (wsUint i=0,N=row() ; i<N ; i++)
		Ra_res[i] = W1 / V[i];
	/* Return */
	return *(new cDiagMatrix(row(), Ra_res));
}

void cDiagMatrix::rem()
{
	if (Ra_data) {
		if (B_dDealloc) { DEALLOC(Ra_data); }
		else sseUnmat(Ra_data, 1);
	}
	Ra_data	= NULL;
	N_row	= 0;
	N_col	= 0;
}

cVector cDiagMatrix::sumR(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	/* Make return buffer */
	wsReal *Ra_ret = sseVector(row());
	/* rowSum of diagonal matrix is matrix itself, so just copy */
	memcpy(Ra_ret, get()[0], sizeof(wsReal)*row());
	*Rp_sumAll = sseVsum(Ra_ret, row());
	if (!NA(R_mul)) sseVpC(Ra_ret, R_mul, Ra_ret, row());
	/* Return */
	return cVector(Ra_ret, row());
}

cVector cDiagMatrix::sumC(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	/* Make return buffer */
	wsReal *Ra_ret = sseVector(row());
	/* rowSum of diagonal matrix is matrix itself, so just copy */
	memcpy(Ra_ret, get()[0], sizeof(wsReal)*row());
	*Rp_sumAll = sseVsum(Ra_ret, row());
	if (!NA(R_mul)) sseVpC(Ra_ret, R_mul, Ra_ret, row());
	/* Return */
	return cVector(Ra_ret, row());
}

/*
# d <- vector, length L
# m <- matrix, dim L*X
DM <- function(d, m) {
L <- length(d); X <- dim(m)[2]
R <- matrix(0, nrow=L, ncol=X)
for (i in 1:L) for (j in 1:X) {//
R[i,j] <- d[i] * m[i,j]
}
}
*/
cStdMatrix cDiagMatrix::operator*(const cStdMatrix &C)
{
	/* 150127 bugfix
	 * Check */
	if (C.row()!=row()) raise(MATERR_DIM_INCORRECT);
	/* Alloc */
	wsReal **Ra_ret = sseMatrix(row(), C.col());
	for (wsUint i=0,N=row() ; i<N ; i++)
		sseVpC(C.get()[i], get()[0][i], Ra_ret[i], C.col());
	/* Return */
	return cStdMatrix(row(), C.col(), Ra_ret);
}

/*
# d <- vector, length L
# s <- symmat, dim L*L
DS <- function(d, s) {
L <- length(d)
R <- matrix(0, nrow=L, ncol=L)
for (i in 1:L) for (j in 1:i) {
R[i,j] <- d[i] * s[i,j]
}
}
*/
cSymMatrix cDiagMatrix::operator*(const cSymMatrix &C)
{
	/* Check */
	if (C.row()!=row()) raise(MATERR_DIM_INCORRECT);
	/* Allocate */
	wsReal **Ra_ret	= sseSymMat(row());
	wsReal **Ra_M	= C.get();
	for (wsUint i=0,N=row() ; i<N ; i++) {
		wsReal R_di = get()[0][i];
		sseVpC(Ra_M[i], R_di, Ra_ret[i], i+1);
	}

	return cSymMatrix(Ra_ret, row());
}

/*
# d <- vector, length L
# e <- vector, length L
DD <- function(d, e) {
L <- length(d)
R <- diag(d*e)
R
}
*/
cDiagMatrix& cDiagMatrix::operator*=(cDiagMatrix &D)
{
	/* Check */
	if (D.row()!=row()) raise(MATERR_DIM_INCORRECT);
	/* Calc */
	sseVpV(get()[0], D.get()[0], get()[0], row());
	/* Does not requires substitution */
	return *this;
}

cDiagMatrix& cDiagMatrix::operator*=(wsRealCst C)
{
	/* Calc */
	sseVpC(get()[0], C, get()[0], row());
	/* Does not requires substitution */
	return *this;
}

cDiagMatrix cDiagMatrix::operator*(const cDiagMatrix &D)
{
	/* Check */
	if (D.row()!=row()) raise(MATERR_DIM_INCORRECT);
	/* Make */
	wsReal		*R		= sseVector(row());
	/* Calc */
	sseVpV(get()[0], D.get()[0], R, row());
	/* Return */
	return cDiagMatrix(row(), R);
}

/*
 *
 * Special operations definition
 *
 */


cStdMatrix cDiagMatrix::mt(const cStdMatrix &M)
{
	if (row()!=M.col()) raise(MATERR_DIM_INCORRECT);
	wsReal *S = get()[0];
	wsReal **R = sseMatrix(row(), M.row());
	wsReal **T = M.get();

	/* R_ij = S_ii * T_ji */
	for (wsUint i=0,N=row() ; i<N ; i++) {
		wsReal R_s = S[i];
		wsReal *Rp = R[i];
		for (wsUint j=0,X=M.row() ; j<X ; j++)
			Rp[j] = R_s * T[j][i]; 
	}

	return cStdMatrix(row(), M.row(), R);
}

wsSym sseDSDt(wsReal *Ra_D, wsUint N_szD, wsSym Ra_S, wsUint N_szS)
{
	if (N_szD != N_szS) halt("Dimension does not match");

	wsSym Ra_R = sseSymMat(N_szD);
	for (wsUint i=0 ; i<N_szS ; i++)
		for (wsUint j=0 ; j<=i ; j++)
			Ra_R[i][j] = Ra_D[j]*Ra_D[i]*Ra_S[i][j];

	return Ra_R;
}

cSymMatrix cDiagMatrix::MMt(const cSymMatrix &S)
{
	if (row()!=S.col()) raise(MATERR_DIM_INCORRECT);
	wsReal **R = sseDSDt(get()[0], row(), S.get(), S.row()); 
	return cSymMatrix(R, row());
}

cDiagMatrix cDiagMatrix::MMt(const cDiagMatrix &D)
{
	if (row()!=D.row()) raise(MATERR_DIM_INCORRECT);
	wsUint	N_med	= 0;
	wsReal	*v1		= get()[0];
	wsReal	*v2		= D.get()[0];
	wsReal	*vt		= sseVector(row());
#ifdef USE_SSE
	N_med = getMed(row());
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t	*sse_v1	= (sse_t *)(v1 + i);
		sse_t	*sse_v2	= (sse_t *)(v2 + i);
		sse_t	*sse_vt	= (sse_t *)(vt + i);
		*sse_vt = sseMul(sseMul(*sse_v1, *sse_v1), *sse_v2);
	}
#endif
	for (wsUint i=N_med ; i<row() ; i++)
		vt[i] = SQR(v1[i]) * v2[i];
	return cDiagMatrix(row(), vt);
}

cDiagMatrix& cDiagMatrix::operator=(cDiagMatrix &D)
{
	init(D.row(), D.get()[0]);
	return *this;
}

cDiagMatrix& cDiagMatrix::operator=(cDiagMatrix &&D)
{
	init(D.row(), D.get()[0]);
	D.setDontDealloc();
	return *this;
}

cSymMatrix cDiagMatrix::operator-(const cSymMatrix &S)
{
	if (row()!=S.col()) raise(MATERR_DIM_INCORRECT);
	wsReal **R = sseSymMat(S.row());
	wsReal **U = S.get();
	wsReal *V = get()[0];

	for (wsUint i=0,N=row() ; i<N ; i++) {
		sseVneg(U[i], i+1, R[i]);
		R[i][i] += V[i];
	}

	return cSymMatrix(R, row());
}

/*
 *
 * Pre-defined operations defition
 *
 */

void cDiagMatrix::file(wsStrCst S_ext, wsUint N_prec/*=0*/)
{
	exportDiagMat(S_ext, get(), row(), N_prec);
}

void cDiagMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/)
{
	exportDiagMat(S_ext, get(), row(), Sa_colNames, N_prec);
}

void cDiagMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/)
{
	/*FIXME*/
	halt("This should not be called");
}

/* CONFIRMED 130314 */
cDiagMatrix& cDiagMatrix::operator-(wsReal R)
{
	/* Prepare */
	wsUint N = row();
	cDiagMatrix *T = new cDiagMatrix(N);
	/* Calc */
	sseVaC(get()[0], -R, T->get()[0], N);
	/* Return */
	return *T;
}

cDiagMatrix& cDiagMatrix::operator+=(wsReal N)
{
	/* Prepare */
	wsUint sz = row();
	/* Calc */
	sseVaC(get()[0], N, get()[0], sz);
	/* Return */
	return *this;
}

cDiagMatrix& cDiagMatrix::operator+=(cDiagMatrix &D)
{
	/* Prepare */
	wsUint sz = row();
	/* Calc */
	sseVaV(get()[0], D.get()[0], get()[0], sz);
	/* Return */
	return *this;
}

wsReal cDiagMatrix::trMMt(cSymMatrix &S)
{
	if (row()!=S.col()) raise(MATERR_DIM_INCORRECT);
	wsReal	R_ret	= W0;
	wsUint	N		= S.col();
	wsReal	**Sp	= S.get();
	wsReal	*Dp		= get()[0];

	for (wsUint i=0 ; i<N ; i++)
		R_ret += SQR(Dp[i]) * Sp[i][i];

	return R_ret;
}

wsReal cDiagMatrix::detL()
{
	wsReal	R_ret	= W0;
	wsReal	*Dp		= get()[0];
	for (wsUint i=0,N=row() ; i<N ; i++)
		R_ret += log(Dp[i]);
	return R_ret;
}

wsReal cDiagMatrix::tr2()
{
	return sseVV(get()[0], row());
}

cVector cDiagMatrix::MV(cDiagMatrix &D, cVector &V)
{
	if (row()!=D.col() || D.col()!=V.size()) raise(MATERR_DIM_INCORRECT);
	wsReal	*Ra_p	= sseVector(row());
	//wsReal	R_ret	= W0;
	wsUint	N		= D.col();
	wsReal	*Vp		= V.get();
	wsReal	*Dp		= D.get()[0];
	wsReal	*Pp		= get()[0];
	wsUint	N_med	= 0;

#ifdef USE_SSE
	N_med = getMed(N);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *R_	= (sse_t *)(Ra_p + i);
		sse_t *S_	= (sse_t *)(Pp + i);
		sse_t *D_	= (sse_t *)(Dp + i);
		sse_t *V_	= (sse_t *)(Vp + i);

		*R_	= sseMul(sseMul(*S_, *D_), *V_);
	}
#endif
	for (wsUint i=N_med ; i<N ; i++)
		Ra_p[i] = Dp[i] * Vp[i] * Pp[i];

	return cVector(Ra_p, N);
}

cDiagMatrix cDiagMatrix::subset(wsUint N_start, wsUint N_sz)
{
	wsReal	*V	= get()[0];
	wsReal	*R	= sseVector(N_sz);

	/* Overflow check */
	if ((N_start+N_sz) > row())
		halt("SYSERR : Subset range (%d~%d) overflow (>%d)", N_start,
			N_start+N_sz-1, row());
	/* Copy */
	memcpy(R, V+N_start, sizeof(wsReal)*N_sz);

	/* Return */
	return cDiagMatrix(N_sz, R);
}

cDiagMatrix cDiagMatrix::inv2()
{
	wsReal	*V	= get()[0];
	wsUint	N	= row();
	wsReal	*R	= sseVector(N);

	/* Copy */
	for (wsUint i=0 ; i<N ; i++)
		if (V[i] == W0)
			R[i] = W0;
		else
			R[i] = W1/V[i];
	/* Return */
	return cDiagMatrix(N, R);
}

/* CONFIRMED 130316 */
cDiagMatrix& cDiagMatrix::operator+(wsReal R)
{
	/* Prepare */
	wsUint	N	= row();
	wsVec	T	= sseVector(N);
	/* Calc */
	sseVaC(get()[0], R, T, N);
	/* Return */
	return *(new cDiagMatrix(N, T));
}

cDiagMatrix cDiagMatrix::operator+(cDiagMatrix &D)
{
	/* Prepare */
	wsUint	N	= row();
	wsReal*	T	= sseVector(N);
	/* Calc */
	sseVaV(get()[0], D.get()[0], T, N);
	/* Return */
	return cDiagMatrix(N, T);
}

/* CONFIRMED 130318 */
cDiagMatrix& cDiagMatrix::operator*(wsReal R)
{
	/* Prepare */
	wsUint N = row();
	cDiagMatrix *T = new cDiagMatrix(N);
	/* Calc */
	sseVpC(get()[0], R, T->get()[0], N);
	/* Return */
	return *T;
}

cVector cDiagMatrix::operator*(cVector &V)
{
	if (row() != V.size())
		halt_fmt(WISARD_SYST_INVL_DIM, "matrix-vector multiplication",
			row(), V.size());
	/* Prepare */
	wsReal *R = NULL;
	sseMalloc(R, wsReal, V.size());
	/* Calc */
	sseVpV(get()[0], V.get(), R, row());
	/* Return */
	return cVector(R, row());
}

cDiagMatrix& cDiagMatrix::addDiag(wsRealCst R)
{
	sseVaC(get()[0], R, get()[0], row());
	return *this;
}

void cDiagMatrix::setR(wsStrCst S_varName)
{
#ifdef USE_R
	if (row() == 1) {
		R().sendMatrix(S_varName, get(), 1, 1);
	} else {
		char S[1024];
		R().sendVector(S_varName, get()[0], row());
		sprintf(S, "%s<-diag(as.numeric(%s))", S_varName, S_varName);
		R().Rparse(S);
	}
#else
	halt_fmt(WISARD_NULL_FEATURE, "R support");
	//halt("USE_R is required");
#endif
}

wsReal cDiagMatrix::tr()
{
	return sseVsum(get()[0], row());
}

cVector cDiagMatrix::diag()
{
	return cVector(sseVectorP(row(), get()[0]), row());
}

cDiagMatrix& cDiagMatrix::clone()
{
	wsReal *Ra_buf = NULL;
	sseMalloc(Ra_buf, wsReal, row());
	memcpy(Ra_buf, get()[0], sizeof(wsReal)*row());
	return *(new cDiagMatrix(row(), Ra_buf));
}

cDiagMatrix& cDiagMatrix::sqrt()
{
	/* Square root of diagonal matrix is just a square root of diagonal terms */
	wsReal *Ra_ret = sseVector(row());
	sseVsqrt(get()[0], row(), Ra_ret);
	return *(new cDiagMatrix(row(), Ra_ret));
}

wsRealCst cDiagMatrix::mean()
{
	wsUint N = row();
	return sum() / (wsReal)(N*N);
}

cVector cDiagMatrix::eigen(cStdMatrix &M_eVec)
{
	wsUint N_sz = row();

	/* Sort diagonal elements */
	xRealSort *Xa_rs = buildRealSort(get()[0], N_sz);
	qsort(Xa_rs, N_sz, sizeof(xRealSort), sort_real_desc);

	/* Matrix = identity permutation */
	M_eVec.init(N_sz, N_sz);
	wsMat Ra_eVec = M_eVec.get();
	for (wsUint i=0 ; i<N_sz ; i++)
		Ra_eVec[i][Xa_rs[i].i] = W1;

	/* Vector = 0 */
	return cVector(Xa_rs, N_sz);
}

cVector cDiagMatrix::eigen()
{
	wsUint N_sz = row();

	/* Sort diagonal elements */
	xRealSort *Xa_rs = buildRealSort(get()[0], N_sz);
	qsort(Xa_rs, N_sz, sizeof(xRealSort), sort_real_desc);

	/* Vector = 0 */
	return cVector(Xa_rs, N_sz);
}

cDiagMatrix& cDiagMatrix::subset(char *Ba_YorN, wsUint N_Yval)
{
	wsUint N_sz = row();
	wsReal* S = get()[0];

	/* Count how much Y's are */
	wsUint N_szAfter = 0;
	for (wsUint i=0 ; i<N_sz ; i++)
		if (Ba_YorN[i] == (char)N_Yval) N_szAfter++;

	/* Make return matrix */
	wsReal *Ra_vec = sseVector(N_szAfter);
	for (wsUint i=0,I=0 ; i<N_sz ; i++)
		if (Ba_YorN[i] == (char)N_Yval)
			Ra_vec[I++] = S[i];

	/* Vector = 0 */
	return *(new cDiagMatrix(N_szAfter, Ra_vec));
}

} // End namespace ONETOOL
