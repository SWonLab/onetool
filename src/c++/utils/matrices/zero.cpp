#include "utils/matrices/zero.h"
#include "utils/matrices/sym.h"
#include "utils/vector.h"
#include "global/Rconn.h"

namespace ONETOOL {

void cZroMatrix::_init(wsUint N_sz)
{
	__init(NULL, N_sz, N_sz, 1);
}

cZroMatrix::cZroMatrix()
{
	__init(NULL, 0, 0, 1);
}

cZroMatrix::cZroMatrix(wsUint N_sz)
{
	_init(N_sz);
}

cZroMatrix::cZroMatrix(cZroMatrix &&I)
{
	_init(I.row());
}

cZroMatrix::~cZroMatrix()
{
	rem();
}

wsMat cZroMatrix::get()
{
	if (Ra_data == NULL)
		Ra_data = sseEmptyMatrix(N_row, N_row);

	return Ra_data;
}

wsRealCst cZroMatrix::sum()
{
	return (wsReal)(N_row);
}

bool cZroMatrix::sym()
{
	/* Since any zero matrix is symmetric, return true */
	return true;
}

cZroMatrix& cZroMatrix::inv()
{
	/* There is NO inversion of ZERO matrix */
	halt("An inversion of null matrix was attempted");
	return *(new cZroMatrix());
}

void cZroMatrix::rem()
{
	N_row = N_col = 0;
}

cVector cZroMatrix::sumR(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	if (Rp_sumAll) *Rp_sumAll = W0;
	return cVector(W0, N_row);
}

cVector cZroMatrix::sumC(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	if (Rp_sumAll) *Rp_sumAll = W0;
	return cVector(W0, N_row);
}

void cZroMatrix::file(wsStrCst S_ext, wsUint N_prec/*=0*/)
{
	exportValMat(W0, S_ext, N_row);
}

void cZroMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/)
{
	exportValMat(W0, S_ext, N_row, Sa_colNames);
}

void cZroMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/)
{
	/*FIXME*/
	halt("This should not be called");
}

cMatrix& cZroMatrix::operator-(wsReal R)
{
	cSymMatrix *Cp_R = new cSymMatrix(N_row);
	wsReal		**Rp = Cp_R->get();
	for (wsUint i=0 ; i<N_row ; i++) {
		for (wsUint j=0 ; j<=i ; j++)
			Rp[i][j] = -R;
//		Rp[i][i] = W1-R;
//		for (wsUint j=i+1 ; j<N_row ; j++)
//			Rp[i][j] = -R;
	}

	return *Cp_R;
}

cMatrix& cZroMatrix::operator+(wsReal R)
{
	cSymMatrix *Cp_R = new cSymMatrix(N_row);
	wsReal		**Rp = Cp_R->get();
	for (wsUint i=0 ; i<N_row ; i++)
		sseVinit(Rp[i], i+1, R);
	return *Cp_R;
}

cMatrix& cZroMatrix::operator*(wsReal R)
{
	return *(new cZroMatrix(row()));
}

void cZroMatrix::setR(wsStrCst S_varName)
{
#ifdef USE_R
	char S[64];
	sprintf(S, "%s<-matrix(0, nrow=%d, ncol=%d)", S_varName, N_row, N_row);
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

wsReal cZroMatrix::tr()
{
	/* Trace of identity matrix is # row itself */
	return (wsReal)row();
}

cVector cZroMatrix::diag()
{
	return cValVector(row(), W0);
}

cZroMatrix& cZroMatrix::clone()
{
	return *(new cZroMatrix(row()));
}

cZroMatrix& cZroMatrix::sqrt()
{
	/* Square root of zero matrix is itself */
	return *(new cZroMatrix(row()));
}

cVector cZroMatrix::operator*(cVector &V)
{
	if (V.size()!=row()) raise(MATERR_DIM_INCORRECT);
	return cVector(row());
//	return V.clone();
}

cMatrix& cZroMatrix::addDiag(wsRealCst R)
{
	halt("This function is still not supported!");
	return *this;
}

wsRealCst cZroMatrix::mean()
{
	return W0;
}

cVector cZroMatrix::eigen(cStdMatrix &M_eVec)
{
	wsUint N_sz = row();

	/* Matrix = identity */
	M_eVec.init(N_sz, N_sz);
	wsMat Ra_eVec = M_eVec.get();
	for (wsUint i=0 ; i<N_sz ; i++)
		Ra_eVec[i][i] = W1;

	/* Vector = 0 */
	return cVector(W0, N_sz);
}

cVector cZroMatrix::eigen()
{
	wsUint N_sz = row();

	/* Vector = 0 */
	return cVector(W0, N_sz);
}

cZroMatrix& cZroMatrix::subset(char *Ba_YorN, wsUint N_Yval)
{
	wsUint N_sz = row();

	/* Count how much Y's are */
	wsUint N_szAfter = 0;
	for (wsUint i=0 ; i<N_sz ; i++)
		if (Ba_YorN[i] == (char)N_Yval) N_szAfter++;

	/* Vector = 0 */
	return *(new cZroMatrix(N_szAfter));
}

} // End namespace ONETOOL
