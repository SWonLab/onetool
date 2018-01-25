#include "utils/matrices/blk.h"
#include "utils/matrices/std.h"
#include "global/Rconn.h"

namespace ONETOOL {

cBlkMatrix::cBlkMatrix(cBlkMatrix &&B)
{
	init(B.get(), B.blk(), B.rep(), B.sym(), 0);
	B.setDontDealloc();
}

cBlkMatrix::cBlkMatrix(wsUint N_szBlk, wsUint N_rep, char B_sym,
	char B_rand/*=0*/)
{
	Ra_data = NULL;
	init(NULL, N_szBlk, N_rep, B_sym, B_rand);
}

cBlkMatrix::cBlkMatrix(wsMat Ra_blkMat, wsUint N_szBlk, wsUint N_rep,
	char B_sym)
{
	Ra_data = NULL;
	init(Ra_blkMat, N_szBlk, N_rep, B_sym, 0);
}

cBlkMatrix::~cBlkMatrix()
{
	rem();
}

void cBlkMatrix::init(wsMat Ra_blkMat, wsUint N_inpSzBlk, wsUint N_inpRep,
	char B_inpSym, char B_rand/*=0*/)
{
	if (B_rand && Ra_blkMat) raise(MATERR_PARAM_INCORRECT);
	if (Ra_data) deallocMatrix(Ra_data, blk(), (void *)1);

	/* Initialize variables */
	this->N_szBlk		= N_inpSzBlk;
	this->N_rep			= N_inpRep;
	this->B_sym			= B_inpSym;
	this->B_ddealloc	= 0;

	/* Make matrix */
//	if (B_sym)
//		Ra_data = sseSymMat(N_szBlk, Ra_blkMat);
//	else

// 	this->N_row			= N_inpSzBlk*N_inpRep;
// 	this->N_col			= N_inpSzBlk*N_inpRep;
// 	Ra_data = ;
	wsUint N_sz = N_inpSzBlk*N_inpRep;
	__init(sseMatrixP(N_szBlk, N_szBlk, Ra_blkMat), N_sz, N_sz);
	if (B_rand) {
		if (B_sym) {
			for (wsUint i=0 ; i<N_szBlk ; i++)
				for (wsUint j=0 ; j<=i ; j++)
					Ra_data[i][j] = (wsReal)(rand()%1000) - REAL_CONST(500.0);
		} else {
			for (wsUint i=0 ; i<N_szBlk ; i++)
				for (wsUint j=0 ; j<N_szBlk ; j++)
					Ra_data[i][j] = (wsReal)(rand()%1000) - REAL_CONST(500.0);
		}
	}
}

wsUintCst cBlkMatrix::blk() const
{
	return N_szBlk;
}

wsUintCst cBlkMatrix::rep() const
{
	return N_rep;
}

cMatrix& cBlkMatrix::operator-(cSymMatrix &S)
{
	/* Check */
	if (row()!=S.row()) raise(MATERR_DIM_INCORRECT);
	/* Calc */
	if (B_sym == 1) {
		/* Return matrix will be sym. */
		cSymMatrix *Cp_ret = new cSymMatrix(row());
		wsReal		**R		= Cp_ret->get();
		wsReal		**Sp	= S.get();
		wsReal		**Bp	= get();
		for (wsUint i=0,E=rep(),B=blk() ; i<E ; i++) {
			for (wsUint j=0 ; j<i ; j++) {
				/* Do negation to off-block-diag part */
				for (wsUint k=0 ; k<B ; k++)
					for (wsUint l=0 ; l<B ; l++)
						R[i*B+k][j*B+l] = -Sp[i*B+k][j*B+l];
			}
			for (wsUint k=0 ; k<B ; k++)
				for (wsUint l=0 ; l<=k ; l++)
					R[i*B+k][i*B+l] = Bp[k][l]-Sp[i*B+k][i*B+l];
		}

		return *Cp_ret;
	} else {
		halt("Still not supported");
	}
	halt("Still not supported");

	return *(new cStdMatrix());
}

cStdMatrix& cBlkMatrix::operator*(cStdMatrix &C)
{
	/* Dimension check */
	if (col() != C.row()) raise(MATERR_DIM_INCORRECT);
	wsReal	**O = C.get();

	/* Make return matrix */
	cStdMatrix	*R	= new cStdMatrix(row(), C.col());
	wsReal	**T	= R->get();

	/* DO calculation */
	for (wsUint i=0 ; i<N_rep ; i++) {
		multMpM(T, row(), C.col(), i*N_szBlk, 0,
			get(), N_szBlk, N_szBlk, 0, 0,
			O, C.row(), C.col(), i*N_szBlk, 0); 
	}

	return *R;
}

cMatrix& cBlkMatrix::operator*=(cStdMatrix &C)
{
	/* Dimension check */
	if (col() != C.row()) raise(MATERR_DIM_INCORRECT);
	wsReal	**O = C.get();

	/* Make return matrix */
	wsReal	**T	= get();

	/* DO calculation */
	for (wsUint i=0 ; i<N_rep ; i++) {
		multMpM(T, row(), C.col(), i*N_szBlk, 0,
			get(), N_szBlk, N_szBlk, 0, 0,
			O, C.row(), C.col(), i*N_szBlk, 0); 
	}
	
	return *this;
}

bool cBlkMatrix::sym()
{
	return B_sym ? true : false;
}

cStdMatrix cBlkMatrix::mt(const cStdMatrix &M)
{
	/* Dimension check */
	if (this->col() != M.col() && this->rep() != M.col()) raise(MATERR_DIM_INCORRECT);
	//wsReal	**O = M.get();

	/* Make return matrix */
// 	cMatrix	*R	= new cStdMatrix(this->row(), C.row());
// 	wsReal	**T	= R->get();
// 
// 	/* DO calculation */
// 	for (wsUint i=0 ; i<N_rep ; i++) {
// 		multMpMt(T, row(), C.row(), i*N_szBlk, 0,
// 			get(), N_szBlk, N_szBlk, 0, 0,
// 			O, C.row(), C.col(), 0, i*N_szBlk);
// 	}
	return cStdMatrix(this->row(), M.row(), sseBpMt(get(), blk(), rep(),
		M.get(), M.row(), M.col()));
}

void cBlkMatrix::sub(cSymMatrix &S)
{

}

cBlkMatrix& cBlkMatrix::operator=(cBlkMatrix &B)
{
	init(B.get(), B.blk(), B.rep(), B.sym());
	return *this;
}

cBlkMatrix& cBlkMatrix::operator=(cBlkMatrix &&B)
{
	init(B.get(), B.blk(), B.rep(), B.sym(), 0);
	B.setDontDealloc();
	return *this;
}

/*
 *
 * Predefined operations definition
 *
 */

wsRealCst cBlkMatrix::sum()
{
	return sseBsum(get(), blk(), rep());
}

cMatrix& cBlkMatrix::inv()
{
	wsMat Ra_ret = NULL;

	if (B_sym)
		Ra_ret = invSymMat(get(), blk());
	else
		Ra_ret = SO_invMatrix(get(), blk());

	return *(new cBlkMatrix(Ra_ret, blk(), rep(), B_sym));
}

void cBlkMatrix::rem()
{
	if (Ra_data) deallocMatrix(Ra_data, blk(), (void *)1);
	Ra_data	= NULL;
	N_row	= 0;
	N_col	= 0;
	N_rep	= 0;
	N_szBlk	= 0;
}

cVector cBlkMatrix::sumR(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsReal *Ra_ret = sseBsumR(get(), blk(), row(), Rp_sumAll);
	if (!NA(R_mul)) sseVpC(Ra_ret, R_mul, Ra_ret, row());
	return cVector(Ra_ret, row());
}

cVector cBlkMatrix::sumC(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsReal *Ra_ret = sseBsumC(get(), blk(), row(), Rp_sumAll);
	if (!NA(R_mul)) sseVpC(Ra_ret, R_mul, Ra_ret, col());
	return cVector(Ra_ret, col());
}

void cBlkMatrix::file(wsStrCst S_ext, wsUint N_prec/*=0*/)
{
	exportBlkMat(S_ext, get(), blk(), rep(), N_prec);
}

void cBlkMatrix::file(wsStrCst S_ext, vStr *Sa_colNames/*=NULL*/, wsUint N_prec/*=0*/)
{
	exportBlkMat(S_ext, get(), blk(), rep(), Sa_colNames, N_prec);
}

void cBlkMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/)
{
	/*FIXME*/
	halt("This should not be called");
}

/* CONFIRMED 130316 */
cBlkMatrix& cBlkMatrix::operator-(wsReal R)
{
	/* Make */
	wsUint B = blk();
	cBlkMatrix *T = new cBlkMatrix(B, rep(), B_sym);
	/* Calc */
	sseMaC(get(), B, B, -R, T->get());
	/* Return */
	return *T;
}

/* CONFIRMED 130316 */
cBlkMatrix& cBlkMatrix::operator+(wsReal R)
{
	/* Make */
	wsUint B = blk();
	cBlkMatrix *T = new cBlkMatrix(B, rep(), B_sym);
	/* Calc */
	sseMaC(get(), B, B, R, T->get());
	/* Return */
	return *T;
}

/* CONFIRMED 130318 */
cBlkMatrix& cBlkMatrix::operator*(wsReal R)
{
	/* Make */
	wsUint B = blk();
	cBlkMatrix *T = new cBlkMatrix(B, rep(), B_sym);
	/* Calc */
	sseMpC(cget(), R, T->get(), B);
	/* Return */
	return *T;
}

cVector cBlkMatrix::operator*(cVector &V)
{
	/* If V is a multiply of this matrix, repeat each element of V
	 * EXACTLY N_szBlk times */
	wsReal R_tRep = row() / (wsReal)(V.size());
	wsUint N_tRep = (wsUint )R_tRep;
	/* Check */
	if (N_tRep!=R_tRep || N_tRep!=blk()) raise(MATERR_DIM_INCORRECT);

	/* Make return vector */
	wsReal	*P	= sseVector(row());

	/* Get rowSums of block matrix */
	wsReal	*Ra_sumR = sseMsumR(get(), row());
	wsReal	*X	= V.get();

	/* i == vector idx
	 * q == result idx */
	for (wsUint q=0,i=0,N=V.size() ; i<N ; i++,q+=blk()) {
		wsUint N_med = 0;
#ifdef USE_SSE
		N_med = getMed(blk());
		sse_t sse_V = sseSet(X[i]);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_T = (sse_t *)(P + q+j);
			sse_t *sse_B = (sse_t *)(Ra_sumR + j);
			*sse_T = sseMul(sse_V, *sse_B);
		}
#endif
		for (wsUint j=N_med ; j<blk() ; j++)
			P[q+j] = X[i]*Ra_sumR[j];
	}

	/* Deallocate */
	sseFree(Ra_sumR);

	return cVector(P, row());
}

cBlkMatrix& cBlkMatrix::addDiag(wsRealCst R)
{
	halt("This function is still not supported!");
	return *this;
}

void cBlkMatrix::setR(wsStrCst S_varName)
{
#ifdef USE_R
	char S[1024];
	if (rep() == 1) {
		R().sendMatrix(S_varName, get(), blk(), blk());
	} else {
		R().sendMatrix("tmp", get(), blk(), blk());
		sprintf(S, "%s <- matrix(0, nrow=%d, ncol=%d)", S_varName, row(), row());
		R().Rparse(S);
		sprintf(S, "eval(parse(text=\"for(i in 1:%d)%s[((i-1)*%d+1):(i*%d),((i-1)*%d+1):(i*%d)]<-tmp\"))",
			rep(), S_varName, blk(), blk(), blk(), blk());
		R().Rparse(S);
		R().Rparse("rm(tmp)");
	}
#else
	halt("USE_R is required");
#endif
}

wsReal cBlkMatrix::tr()
{
	wsUint B = blk();
	wsReal **R = get();
	wsReal R_ret = W0;
	for (wsUint i=0 ; i<B ; i++)
		R_ret += R[i][i];
	return R_ret*rep();
}

cVector cBlkMatrix::diag()
{
	wsUint	B		= blk();
	wsMat	R		= get();
	wsVec	Ra_ret	= sseVector(row());

	LOOP (i, B) {
		wsReal R_v = R[i][i];
		for (wsUint j=i ; j<rep() ; j+=B)
			Ra_ret[j] = R_v;
	}

	return cVector(Ra_ret, row());
}

cBlkMatrix& cBlkMatrix::clone()
{
	wsReal **Ra_ret = NULL;
	if (sym())
		Ra_ret = sseSymMatP(blk(), get());
	else
		Ra_ret = sseMatrixP(blk(), blk(), get());

	return *(new cBlkMatrix(Ra_ret, blk(), rep(), sym()));
}

cBlkMatrix& cBlkMatrix::sqrt()
{
	wsReal **Ra_ret = NULL;
	if (sym())
		Ra_ret = sqrtMatrix(get(), blk(), 1);
	else
		halt("SYSERR : Unsupported feature [Square root of block matrix with nonsymmetric matrix] called");

	return *(new cBlkMatrix(Ra_ret, blk(), rep(), sym()));
}

wsRealCst cBlkMatrix::mean()
{
	wsUint N = row();
	return sum() / (wsReal)(N*N);
}

cVector cBlkMatrix::eigen(cStdMatrix &Mp_eVec)
{
	/* FIXME : Implement this */
	halt_fmt(WISARD_SYST_NULL_IMPLEMENT, "Eigendecomposition of block-diagonal matrix");
	return cVector();
}

cVector cBlkMatrix::eigen()
{
	/* FIXME : Implement this */
	halt_fmt(WISARD_SYST_NULL_IMPLEMENT, "Eigendecomposition of block-diagonal matrix");
	return cVector();
}

cMatrix& cBlkMatrix::subset(char *Ba_YorN, wsUint N_Yval)
{
	/* FIXME : Implement this */
	halt_fmt(WISARD_SYST_NULL_IMPLEMENT, "Subsetting of block-diagonal matrix");
	return *(new cSymMatrix());
}

} // End namespace ONETOOL
