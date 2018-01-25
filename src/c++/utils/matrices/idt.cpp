#include "utils/matrices/idt.h"
#include "utils/matrices/sym.h"
#include "utils/vector.h"
#include "global/Rconn.h"

namespace ONETOOL {

void cIdtMatrix::_init(wsUint N_sz)
{
	__init(NULL, N_sz, N_sz, 1);
}

cIdtMatrix::cIdtMatrix(wsUint N_sz)
{
	_init(N_sz);
}

cIdtMatrix::cIdtMatrix(cIdtMatrix &&I)
{
	_init(I.row());
}

cIdtMatrix::~cIdtMatrix()
{
	rem();
}

wsMat cIdtMatrix::get()
{
	if (Ra_data == NULL) {
		Ra_data = sseEmptyMatrix(N_row, N_row);
		for (wsUint i=0 ; i<N_row ; i++)
			Ra_data[i][i] = W1;
	}

	return Ra_data;
}

wsRealCst cIdtMatrix::sum()
{
	return (wsReal)(N_row);
}

bool cIdtMatrix::sym()
{
	/* Since any identity matrix is symmetric, return true */
	return true;
}

cIdtMatrix& cIdtMatrix::inv()
{
	return *(new cIdtMatrix(N_row));
}

void cIdtMatrix::rem()
{
	N_row = N_col = 0;
}

cVector cIdtMatrix::sumR(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	if (Rp_sumAll) *Rp_sumAll = (wsReal)N_row;
	if (NA(R_mul))
		return cVector(W1, N_row);
	else
		return cVector(R_mul, N_row);
}

cVector cIdtMatrix::sumC(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	if (Rp_sumAll) *Rp_sumAll = (wsReal)N_row;
	if (NA(R_mul))
		return cVector(W1, N_row);
	else
		return cVector(R_mul, N_row);
}

void cIdtMatrix::file(wsStrCst S_ext, wsUint N_prec/*=0*/)
{
	exportIdtMat(S_ext, N_row);
}

void cIdtMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/)
{
	exportIdtMat(S_ext, N_row, Sa_colNames);
}

void cIdtMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/)
{
	exportIdtMat(S_ext, N_row, Sa_colNames, Sa_rowNames);
}

cMatrix& cIdtMatrix::operator-(wsReal R)
{
	cSymMatrix *Cp_R = new cSymMatrix(N_row);
	wsReal		**Rp = Cp_R->get();
	for (wsUint i=0 ; i<N_row ; i++) {
		for (wsUint j=0 ; j<i ; j++)
			Rp[i][j] = -R;
		Rp[i][i] = W1-R;
//		for (wsUint j=i+1 ; j<N_row ; j++)
//			Rp[i][j] = -R;
	}

	return *Cp_R;
}

cMatrix& cIdtMatrix::operator+(wsReal R)
{
	cSymMatrix *Cp_R = new cSymMatrix(N_row);
	wsReal		**Rp = Cp_R->get();
	for (wsUint i=0 ; i<N_row ; i++) {
		for (wsUint j=0 ; j<i ; j++)
			Rp[i][j] = R;
		Rp[i][i] = W1+R;
		for (wsUint j=i+1 ; j<N_row ; j++)
			Rp[i][j] = R;
	}

	return *Cp_R;
}

cMatrix& cIdtMatrix::operator*(wsReal R)
{
	cDiagMatrix *Cp_R = new cDiagMatrix(R, N_row);
	return *Cp_R;
}

void cIdtMatrix::setR(wsStrCst S_varName)
{
#ifdef USE_R
	char S[64];
	sprintf(S, "%s<-diag(%d)", S_varName, N_row);
	R().Rparse(S);
#else
	halt("R functionality required");
#endif
}

/*
 *
 * Special operations definition
 *
 */

wsRealCst cIdtMatrix::normF(cStdMatrix &M)
{
	wsUint i, j, N;
	if (!M.sqr() || row()!=M.row()) raise(MATERR_DIM_INCORRECT);

	/* Off-diagonal -Mij^2 = Mij^2
	 *     Diagonal (1-Mii)^2 = Mii^2-2Mii+1
	 */
	wsReal R_ret = W0;
	for (i=0,N=row() ; i<N ; i++) {
		wsReal *R = M.get()[i];
		R_ret += sseVV(R, i);
		/* Diagonal */
		R_ret += SQR(R[i]) - W2*R[i] + W1;
		
		for (j=i+1 ; (size_t)(R+j)%16 ; j++)
			R_ret += SQR(R[j]);
		R_ret += sseVV(R+j, N-j);
	}

	return R_ret;
}

wsRealCst cIdtMatrix::normF(cSymMatrix &M)
{
	wsUint i, N;
	if (row()!=M.row()) raise(MATERR_DIM_INCORRECT);

	/* Off-diagonal -Mij^2 = Mij^2
	 *     Diagonal (1-Mii)^2 = Mii^2-2Mii+1
	 */
	wsReal R_ret = W0;
	for (i=0,N=row() ; i<N ; i++) {
		wsReal *R = M.get()[i];
		R_ret += W2 * sseVV(R, i);
		/* Diagonal */
		R_ret += SQR(R[i]) - W2*R[i] + W1;
	}

	return R_ret;
}

/* IM = M itself */
cStdMatrix& cIdtMatrix::mt(cStdMatrix &M)
{
	if (row()!=M.row()) raise(MATERR_DIM_INCORRECT);

	wsReal **Ra_ret = sseMatrixP(M.row(), M.col(), M.get());
	return *(new cStdMatrix(M.row(), M.col(), Ra_ret));
}

/* I-S */
void cIdtMatrix::sub(cSymMatrix &S)
{
}

cSymMatrix cIdtMatrix::operator-(const cSymMatrix &S)
{
	if (row()!=S.row()) raise(MATERR_DIM_INCORRECT);
	wsUint N	= row();
	wsReal **R	= sseSymMat(N);

	for (wsUint i=0 ; i<N ; i++) {
		wsReal *T = S.get()[i];
		sseVneg(T, i, R[i]);
		R[i][i] = W1 - T[i];
	}

	return cSymMatrix(R, N);
}

cZroMatrix cIdtMatrix::operator-(const cIdtMatrix &S)
{
	if (row()!=S.row()) raise(MATERR_DIM_INCORRECT);

	return cZroMatrix(row());
}

cIdtMatrix& cIdtMatrix::operator=(cIdtMatrix &I)
{
	_init(I.row());
	return *this;
}

cIdtMatrix& cIdtMatrix::operator=(cIdtMatrix &&I)
{
	_init(I.row());
	return *this;
}

cSymMatrix& cIdtMatrix::operator*(cSymMatrix &S)
{
	if (S.row()!=row()) raise(MATERR_DIM_INCORRECT);
	return S.clone();
}

cDiagMatrix cIdtMatrix::krnk(cDiagMatrix &D)
{
	wsUint	N = D.row() * row();
	wsVec	R = sseVrep(D.get()[0], D.row(), row());
	return cDiagMatrix(N, R);
}

cStdMatrix cIdtMatrix::rmerge(cStdMatrix& M)
{
	if (row() != M.row()) raise(MATERR_DIM_INCORRECT);
	wsUintCst N = row();
	wsUintCst C = col();
	wsMat Ra_src = M.get();
	wsMat Ra_ret = sseMatrix(N, C+M.col());
	LOOP (i, N) {
		memset(Ra_ret[i], 0x00, sizeof(wsReal)*C);
		Ra_ret[i][i] = W1;
		memcpy(Ra_ret[i]+C, Ra_src[i], sizeof(wsReal)*M.col());
	}

	return cStdMatrix(N, C+M.col(), Ra_ret);
}

cStdMatrix cIdtMatrix::bmerge(cStdMatrix& M)
{
	wsUintCst C = col();
	if (C != M.col()) raise(MATERR_DIM_INCORRECT);

	wsUintCst N = row();
	wsUintCst Ns = M.row();
	wsMat Ra_src = M.get();
	wsMat Ra_ret = sseMatrix(N+M.row(), C);
	LOOP(i, N) {
		memset(Ra_ret[i], 0x00, sizeof(wsReal)*C);
		Ra_ret[i][i] = W1;
	}
	LOOP(i, Ns)
		memcpy(Ra_ret[i+N], Ra_src[i], sizeof(wsReal)*C);

	return cStdMatrix(N+M.row(), C, Ra_ret);
}

cVector cIdtMatrix::operator*(cVector &V)
{
	if (V.size()!=row()) raise(MATERR_DIM_INCORRECT);
	return V.clone();
}

cMatrix& cIdtMatrix::addDiag(wsRealCst R)
{
	halt("This function is still not supported!");
	return *this;
}

wsReal cIdtMatrix::tr()
{
	/* Trace of identity matrix is # row itself */
	return (wsReal)row();
}

cVector cIdtMatrix::diag()
{
	return cUnitVector(row());
}

cIdtMatrix& cIdtMatrix::clone()
{
	return *(new cIdtMatrix(row()));
}

cIdtMatrix& cIdtMatrix::sqrt()
{
	/* Square root of identity matrix is itself */
	return *(new cIdtMatrix(row()));
}

wsRealCst cIdtMatrix::mean()
{
	return W1 / (wsReal)row();
}

cVector cIdtMatrix::eigen(cStdMatrix &M_eVec)
{
	wsUint N_sz = row();

	/* Matrix = identity */
	M_eVec.init(N_sz, N_sz);
	wsMat Ra_eVec = M_eVec.get();
	for (wsUint i=0 ; i<N_sz ; i++)
		Ra_eVec[i][i] = W1;

	/* Vector = 1 */
	return cVector(W1, N_sz);
}

cVector cIdtMatrix::eigen()
{
	wsUint N_sz = row();

	/* Vector = 1 */
	return cVector(W1, N_sz);
}

cIdtMatrix& cIdtMatrix::subset(char *Ba_YorN, wsUint N_Yval)
{
	wsUint N_sz = row();

	/* Count how much Y's are */
	wsUint N_szAfter = 0;
	for (wsUint i=0 ; i<N_sz ; i++)
		if (Ba_YorN[i] == (char)N_Yval) N_szAfter++;

	/* Vector = 0 */
	return *(new cIdtMatrix(N_szAfter));
}

} // End namespace ONETOOL
