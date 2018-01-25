#include "utils/matrices/one.h"
#include "utils/matrices/zero.h"
#include "utils/matrices/sym.h"
#include "utils/vector.h"
#include "global/Rconn.h"

namespace ONETOOL {

void cOneMatrix::_init(wsUint N_inpRow, wsUint N_inpCol)
{
	__init(NULL, N_inpRow, N_inpCol, 1);
}

cOneMatrix::cOneMatrix()
{
	__init(NULL, 0, 0, 1);
}

cOneMatrix::cOneMatrix(wsUint N_inpRow, wsUint N_inpCol)
{
	_init(N_inpRow, N_inpCol);
}

cOneMatrix::cOneMatrix(cOneMatrix &&O)
{
	_init(O.row(), O.col());
}

cOneMatrix::~cOneMatrix()
{
	rem();
}

wsReal** cOneMatrix::get()
{
	if (Ra_data == NULL) {
		Ra_data = sseMatrix(N_row, N_col);
		for (wsUint i=0 ; i<N_row ; i++)
			sseVinit(Ra_data[i], N_col, W1);
	}

	return Ra_data;
}

wsRealCst cOneMatrix::sum()
{
	return (wsReal)(N_row*N_col);
}

bool cOneMatrix::sym()
{
	return N_row == N_col;
}

cOneMatrix& cOneMatrix::inv()
{
	/* There is NO inversion of ZERO matrix */
	halt("An inversion of 1 matrix was attempted");
	return *(new cOneMatrix());
}

void cOneMatrix::rem()
{
	N_row = N_col = 0;
}

cVector cOneMatrix::sumR(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	if (Rp_sumAll) *Rp_sumAll = (wsReal)N_row * (wsReal)N_col;
	if (NA(R_mul))
		return cVector((wsReal)N_col, N_row);
	else
		return cVector((wsReal)N_col * R_mul, N_row);
}

cVector cOneMatrix::sumC(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	if (Rp_sumAll) *Rp_sumAll = (wsReal)N_row * (wsReal)N_col;
	if (NA(R_mul))
		return cVector((wsReal)N_col, N_row);
	else
		return cVector((wsReal)N_col * R_mul, N_row);
}

void cOneMatrix::file(wsStrCst S_ext, wsUint N_prec/*=0*/)
{
	exportValMat(W1, S_ext, N_row);
}

void cOneMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/)
{
	exportValMat(W1, S_ext, N_row, Sa_colNames);
}

void cOneMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/)
{
	/*FIXME*/
	halt("This should not be called");
}

cMatrix& cOneMatrix::operator-(wsReal R)
{
	cSymMatrix *Cp_R = new cSymMatrix(N_row);
	wsReal		**Rp = Cp_R->get();
	for (wsUint i=0 ; i<N_row ; i++)
		sseVinit(Rp[i], i+1, W1-R);
	return *Cp_R;
}

cMatrix& cOneMatrix::operator+(wsReal R)
{
	cSymMatrix *Cp_R = new cSymMatrix(N_row);
	wsReal		**Rp = Cp_R->get();
	for (wsUint i=0 ; i<N_row ; i++)
		sseVinit(Rp[i], i+1, W1+R);
	return *Cp_R;
}

/* FIXME : Implement below */

cMatrix& cOneMatrix::operator*(wsReal R)
{
	return *(new cOneMatrix(row(), col()));
}

cMatrix& cOneMatrix::addDiag(wsRealCst R)
{
	halt("This function is still not supported!");
	return *this;
}

void cOneMatrix::setR(wsStrCst S_varName)
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

wsReal cOneMatrix::tr()
{
	/* Trace of identity matrix is # row itself */
	return (wsReal)row();
}

cVector cOneMatrix::diag()
{
	return cUnitVector(row());
}

cOneMatrix& cOneMatrix::clone()
{
	return *(new cOneMatrix(row(), col()));
}

cOneMatrix& cOneMatrix::sqrt()
{
	/* Square root of zero matrix is itself */
	return *(new cOneMatrix(row(), col()));
}

cVector cOneMatrix::operator*(cVector &V)
{
	if (V.size()!=row()) raise(MATERR_DIM_INCORRECT);
	return cVector(row());
//	return V.clone();
}

wsRealCst cOneMatrix::mean()
{
	return W0;
}

cVector cOneMatrix::eigen(cStdMatrix &M_eVec)
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

cVector cOneMatrix::eigen()
{
	wsUint N_sz = row();

	/* Vector = 0 */
	return cVector(W0, N_sz);
}

cOneMatrix& cOneMatrix::subset(char *Ba_YorN, wsUint N_Yval)
{
	wsUint N_sz = row();

	/* Count how much Y's are */
	wsUint N_szAfter = 0;
	for (wsUint i=0 ; i<N_sz ; i++)
		if (Ba_YorN[i] == (char)N_Yval) N_szAfter++;

	/* Vector = 0 */
	return *(new cOneMatrix(N_szAfter, N_szAfter));
}

} // End namespace ONETOOL
