#include "utils/vector.h"
#include "utils/matrices/sparse.h"

namespace ONETOOL {

cSpsMatrix::cSpsMatrix(wsUint N_r, wsUint N_c, wsUint N_sz/*=0*/)
{
	Na_rEnd = NULL;
	Na_cIdx = NULL;

	init(N_r, N_c, N_sz);
}

cSpsMatrix::cSpsMatrix(cStdMatrix &M)
{
	wsUint	N_r		= M.row();
	wsUint	N_c		= M.col();
	wsMat	Rp_mat	= M.get();

	Na_rEnd = NULL;
	Na_cIdx = NULL;

	init(N_r, N_c);

	vInt	Nv_cIdx;
	vInt	Nv_rEnd;
	vReal	Rv_val;

	/* Identifying data */
	for (wsUint i=0 ; i<N_r ; i++) {
		wsReal* Rp_m = Rp_mat[i];
		for (wsUint j=0 ; j<N_c ; j++)
			if (Rp_m[j] != W0) {
				Nv_cIdx.push_back(j);
				Rv_val.push_back(Rp_m[j]);
			}
		Nv_rEnd.push_back((int)Rv_val.size());
	}

	/* Size check */
	N_size = (wsUint)Rv_val.size();
	if ((wsReal)N_size / ((wsReal)row() * (wsReal)col()) > REAL_CONST(0.1))
		LOGwarn("A sparse matrix over 10%% filled was made\n");

	/* Filling data */
	wsAlloc(Na_rEnd, wsUint, N_r);
	wsAlloc(Na_cIdx, wsUint, N_size);
	Ra_data = sseMatrix(1, N_size);
	wsUint i = 0;
	FOREACHDO (vInt_it, Nv_cIdx, it, i++) {
		Na_cIdx[i] = *it;
		Ra_data[0][i] = Rv_val[i];
	}
	i = 0;
	FOREACHDO (vInt_it, Nv_rEnd, it, i++)
		Na_rEnd[i] = *it;
}

cSpsMatrix::cSpsMatrix(cSymMatrix &S)
{
	wsUint	N_r		= S.row();
	wsSym	Rp_mat	= S.get();

	Na_rEnd = NULL;
	Na_cIdx = NULL;

	init(N_r, N_r);
	/* It is symmetric */
	B_sym = 1;

	vInt	Nv_cIdx;
	vInt	Nv_rEnd;
	vReal	Rv_val;

	/* Identifying data */
	for (wsUint i=0 ; i<N_r ; i++) {
		wsReal* Rp_m = Rp_mat[i];
		for (wsUint j=0 ; j<=i ; j++)
			if (Rp_m[j] != W0) {
				Nv_cIdx.push_back(j);
				Rv_val.push_back(Rp_m[j]);
			}
		Nv_rEnd.push_back((int)Rv_val.size());
	}

	/* Size check */
	N_size = (wsUint)Rv_val.size();
	if ((wsReal)N_size / ((wsReal)N_r * (wsReal)(N_r+1) / W2) > REAL_CONST(0.1))
		LOGwarn("A symmetric sparse matrix over 10% filled was made\n");

	/* Filling data */
	wsAlloc(Na_rEnd, wsUint, N_r);
	wsAlloc(Na_cIdx, wsUint, N_size);
	Ra_data = sseMatrix(1, N_size);
	wsUint i = 0;
	FOREACHDO (vInt_it, Nv_cIdx, it, i++) {
		Na_cIdx[i] = *it;
		Ra_data[0][i] = Rv_val[i];
	}
	i = 0;
	FOREACHDO (vInt_it, Nv_rEnd, it, i++)
		Na_rEnd[i] = *it;
}

cSpsMatrix::cSpsMatrix(cDiagMatrix &D)
{
	wsUint i, _i, N=D.row();

	Na_rEnd = NULL;
	Na_cIdx = NULL;

	/* Initialize buffer */
	init(N, N, N);
	/* it is symmetric */
	B_sym = 1;

	/* Count # of elements */
	wsReal *Ra_V = D.get()[0];
	wsReal *Ra_T = get()[0];
	for (i=_i=0 ; i<N ; i++) {
		if (Ra_V[i] != W0) {
			Ra_T[_i]	= Ra_V[i];
			Na_cIdx[_i++] = i;
		}
		Na_rEnd[i] = _i;
	}
}

cSpsMatrix::cSpsMatrix(cSpsMatrix& P)
{
	wsUint N = P.size();

	Na_rEnd = NULL;
	Na_cIdx = NULL;

	/* Init class and data */
	init(P.row(), P.col(), N);
	memcpy(get()[0], P.get()[0], sizeof(wsReal)*N);

	wsUint* Na_orEnd = P.getRowEnds();
	wsUint* Na_ocIdx = P.getColIdxs();

	memcpy(Na_rEnd, Na_orEnd, sizeof(wsUint)*P.row());
	memcpy(Na_cIdx, Na_ocIdx, sizeof(wsUint)*N);
}

cSpsMatrix::cSpsMatrix(cSpsMatrix& P, wsVec Ra_new)
{
	wsUint N = P.size();

	Na_rEnd = NULL;
	Na_cIdx = NULL;

	/* Init class and data */
	init(P.row(), P.col(), N);
	memcpy(get()[0], Ra_new, sizeof(wsReal)*N);

	wsUint* Na_orEnd = P.getRowEnds();
	wsUint* Na_ocIdx = P.getColIdxs();

	memcpy(Na_rEnd, Na_orEnd, sizeof(wsUint)*P.row());
	memcpy(Na_cIdx, Na_ocIdx, sizeof(wsUint)*N);
}

cSpsMatrix::cSpsMatrix(cSpsMatrix&& P)
{
	P.setDontDealloc();
        wsUint N = P.size();

        Na_rEnd = NULL;
        Na_cIdx = NULL;

        /* Init class and data */
        init(P.row(), P.col(), N);
        memcpy(get()[0], P.get()[0], sizeof(wsReal)*N);

        wsUint* Na_orEnd = P.getRowEnds();
        wsUint* Na_ocIdx = P.getColIdxs();

        memcpy(Na_rEnd, Na_orEnd, sizeof(wsUint)*P.row());
        memcpy(Na_cIdx, Na_ocIdx, sizeof(wsUint)*N);
}

cSpsMatrix::cSpsMatrix(cSpsMatrix&& P, wsVec Ra_new)
{
	P.setDontDealloc();
        wsUint N = P.size();

        Na_rEnd = NULL;
        Na_cIdx = NULL;

        /* Init class and data */
        init(P.row(), P.col(), N);
        memcpy(get()[0], Ra_new, sizeof(wsReal)*N);

        wsUint* Na_orEnd = P.getRowEnds();
        wsUint* Na_ocIdx = P.getColIdxs();

        memcpy(Na_rEnd, Na_orEnd, sizeof(wsUint)*P.row());
        memcpy(Na_cIdx, Na_ocIdx, sizeof(wsUint)*N);
}

void cSpsMatrix::init(wsUint N_r, wsUint N_c, wsUint N_sz/*=0*/)
{
	sseFree(Na_rEnd);
	sseFree(Na_cIdx);
	sseUnmat(Ra_data, 1);

	B_sym	= 0;
	N_row	= N_r;
	N_col	= N_c;
	Ra_data	= NULL;
	N_size	= N_sz;

	sseMalloc(Na_rEnd, wsUint, N_r);
	sseMalloc(Na_cIdx, wsUint, N_size);
	if (N_size)
		Ra_data = sseMatrix(1, N_size);
}

wsUintCst cSpsMatrix::size()
{
	return N_size;
}

wsUint* cSpsMatrix::getRowEnds() const
{
	return Na_rEnd;
}

wsUint* cSpsMatrix::getColIdxs() const
{
	return Na_cIdx;
}

cStdMatrix cSpsMatrix::toStd()
{
	wsUint	N_r		= row();
	wsUint	N_c		= col();
	wsMat	Ra_ret	= sseEmptyMatrix(N_r, N_c);

	wsUint	i		= 0;
	wsReal*	Rp_data	= get()[0];
	wsUint	N_s		= 0;
	wsUint	N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e)
			for (wsUint j=N_s ; j<N_e ; j++)
				Ra_ret[i][Na_cIdx[j]] += Rp_data[j];

		/* Update start point */
		i++;
		if (i < N_r)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	return cStdMatrix(N_r, N_c, Ra_ret);
}

cStdMatrix cSpsMatrix::operator*(cStdMatrix& M)
{
	wsUint	N_r		= row();
	wsUint	N_c		= col();
	wsMat	Ra_ret	= sseEmptyMatrix(N_r, N_c);

	wsUint	i		= 0;
	wsReal*	Rp_data	= get()[0];
	wsUint	N_s		= 0;
	wsUint	N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e)
			for (wsUint j=N_s ; j<N_e ; j++)
				LOOP (k, M.row())
					sseVpCadd(M.get()[k], Rp_data[j], Ra_ret[k], M.col());
				//Ra_ret[i][Na_cIdx[j]] += Rp_data[j];
	
		/* Update start point */
		i++;
		if (i < N_r)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	return cStdMatrix(N_r, N_c, Ra_ret);
}

cSpsMatrix cSpsMatrix::operator*(cDiagMatrix& D)
{
	wsUint	N_r		= row();
	wsVec	Ra_ret	= sseVector(size());
	wsVec	Dp		= D.get()[0];

	wsUint	i		= 0;
	wsReal*	Rp_data	= get()[0];
	wsUint	N_s		= 0;
	wsUint	N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e)
			for (wsUint j=N_s ; j<N_e ; j++)
				Ra_ret[j] += Rp_data[j] * Dp[Na_cIdx[j]];

		/* Update start point */
		i++;
		if (i < N_r)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	return cSpsMatrix(*this, Ra_ret);
}

cSpsMatrix cSpsMatrix::operator*(cIdtMatrix& I)
{
	return cSpsMatrix(*this);
}

wsRealCst cSpsMatrix::sum()
{
	return sseVsum(Ra_data[0], size());
}

bool cSpsMatrix::sym()
{
	/* Need to check all elements are symmetric....
	 * FIXME : Should be, currently returns false */
	return false;
}

cMatrix& cSpsMatrix::inv()
{
	halt("This operation is not supported!");
	return *(new cStdMatrix());
}

void cSpsMatrix::rem()
{
	__init(NULL, 0, 0, 1);
}

cVector cSpsMatrix::sumR(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsUint	N_r		= row();
	wsReal*	Rp_ret	= sseVector(N_r);
	wsUint	i		= 0;
	wsReal*	Rp_data	= get()[0];
	wsUint	N_s		= 0;
	wsUint	N_e;
	do {
		wsReal	R_ret	= REAL_CONST(0.0);
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e)  {
			/* Sum up them */
			for (wsUint j=N_s ; j<N_e ; j++)
				R_ret += Rp_data[j];
		}
		Rp_ret[i] = R_ret;

		/* Update start point */
		i++;
		if (i < N_r)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	if (!NA(R_mul)) sseVpC(Rp_ret, R_mul, Rp_ret, N_r);
	if (Rp_sumAll) *Rp_sumAll = sseVsum(Rp_ret, N_r);

	return cVector(Rp_ret, N_r);
}

cVector cSpsMatrix::sumC(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsUint	N_c		= col();
	wsUint	N_r		= row();
	wsReal*	Rp_ret	= sseVector(N_c);
	wsUint	i		= 0;
	wsReal*	Rp_data	= get()[0];
	wsUint	N_s		= 0;
	wsUint	N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e)
			/* Sum up them */
			for (wsUint j=N_s ; j<N_e ; j++)
				Rp_ret[Na_cIdx[j]] += Rp_data[j];

		/* Update start point */
		i++;
		if (i < N_r)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	if (!NA(R_mul)) sseVpC(Rp_ret, R_mul, Rp_ret, N_c);
	if (Rp_sumAll) *Rp_sumAll = sseVsum(Rp_ret, N_c);

	return cVector(Rp_ret, N_c);
}

void cSpsMatrix::file(wsStrCst S_ext, wsUint N_prec/*=0*/)
{
	exportSpsMat(S_ext, Ra_data[0], row(), col(), Na_rEnd, Na_cIdx, NULL, N_prec);
}

void cSpsMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/)
{
	exportSpsMat(S_ext, Ra_data[0], row(), col(), Na_rEnd, Na_cIdx, Sa_colNames,
		N_prec);
}

void cSpsMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/)
{
	/*FIXME*/
	halt("This should not be called");
}

cMatrix& cSpsMatrix::operator-(wsReal R)
{
	halt("This operation is not supported!");
	return *(new cStdMatrix());
}

cMatrix& cSpsMatrix::operator+(wsReal R)
{
	halt("This operation is not supported!");
	return *(new cStdMatrix());
}

cSpsMatrix& cSpsMatrix::operator*(wsReal R)
{
	cSpsMatrix*	Cp_ret	= new cSpsMatrix(*this);
	wsReal*		Rp_dat	= Cp_ret->get()[0];
	wsUint		N_sz	= Cp_ret->size();

	sseVpC(Rp_dat, R, Rp_dat, N_sz);

	return *Cp_ret;
}

cVector cSpsMatrix::operator*(cVector &V)
{
	if (V.size() != col()) raise(MATERR_DIM_INCORRECT);
	wsUint		N		= row();
	wsReal*		Rp_dat	= sseVector(N);
	wsReal*		T		= V.get();
	wsUint		N_sz	= size();

	for (wsUint i=0 ; i<N_sz ; i++)
		Rp_dat[i] *= T[Na_cIdx[i]];

	return cVector(Rp_dat, N);
}

cSpsMatrix& cSpsMatrix::addDiag(wsRealCst R)
{
	halt("This function is still not supported!");
	return *this;
}

void cSpsMatrix::setDontDealloc()
{
	B_dtdealloc = 1;
}

void cSpsMatrix::setR(wsStrCst S_varName)
{
	halt("This operation is not supported!");
}

wsReal cSpsMatrix::tr()
{
	/* Print values along with R start point */
	wsUint	i		= 0;
	wsReal	R_ret	= REAL_CONST(0.0);
	wsReal*	Rp_data	= get()[0];
	wsUint	N_r		= row();
	wsUint	N_s		= 0;
	wsUint	N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e)  {
			/* From first column to last column */
			wsUint k = 0;
			for (wsUint j=N_s ; j<N_e ; j++) {
				/* Front fill */
				for ( ; k<Na_cIdx[j] ; k++) {
					if (Na_cIdx[j] > i) break;
					else if (Na_cIdx[j] == i) {
						R_ret += Rp_data[j];
						break;
					}
				}
			}
		}

		/* Update start point */
		i++;
		if (i < N_r)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	return R_ret;
}

cVector cSpsMatrix::diag()
{
	/* Print values along with R start point */
	wsUint	i		= 0;
	wsReal*	Rp_data	= get()[0];
	wsUint	N_r		= row();
	wsVec	Ra_ret	= sseVector(N_r);
	wsUint	N_s		= 0;
	wsUint	N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e)  {
			/* From first column to last column */
			wsUint k = 0;
			for (wsUint j=N_s ; j<N_e ; j++) {
				/* Front fill */
				for ( ; k<Na_cIdx[j] ; k++) {
					if (Na_cIdx[j] > i) break;
					else if (Na_cIdx[j] == i) {
						Ra_ret[i] = Rp_data[j];
						break;
					}
				}
			}
		} else
			Ra_ret[i] = W0;

		/* Update start point */
		i++;
		if (i < N_r)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	return cVector(Ra_ret, N_r);
}

cSpsMatrix& cSpsMatrix::clone()
{
	return *(new cSpsMatrix(*this));
}

cMatrix& cSpsMatrix::sqrt()
{
	halt("This operation is not supported!");
	return *(new cStdMatrix());
}

wsRealCst cSpsMatrix::mean()
{
	return sseVsum(get()[0], size()) / (wsReal)size();
}

cVector cSpsMatrix::eigen(cStdMatrix &Mp_eVec)
{
	halt("This operation is not supported!");
	return cVector();
}

cVector cSpsMatrix::eigen()
{
	halt("This operation is not supported!");
	return cVector();
}

cSpsMatrix& cSpsMatrix::subset(char *Ba_YorN, wsUint N_Yval)
{
	/* Print values along with R start point */
	wsUint	i		= 0;
	wsReal*	Rp_data	= get()[0];
	wsUint	N_r		= row();
	wsUint	N_c		= col();
	wsUint	N_s		= 0;
	wsUint	N_e;
	vReal	Xv_val;
	vInt	Xv_rEnd;
	vInt	Xv_cIdx;
	char	Y		= (char)N_Yval;

	if (N_r != N_c) halt("subset() requires square matrix");

	/* Count up the number of Ys */
	wsUint	N_nr = 0;
	for (wsUint i=0 ; i<N_r ; i++)
		if (Ba_YorN[i] == Y) N_nr++;

	do {
		/* Set end point */
		N_e = Na_rEnd[i];
		/* No element in the row */
		if (N_s != N_e && 
		/* If this row is not choice then skip */
			Ba_YorN[i] == Y) {
			/* From first column to last column */
//			wsUint k = 0;
			for (wsUint j=N_s ; j<N_e ; j++) {
				wsUint N_c = Na_cIdx[j];
				if (Ba_YorN[N_c] != Y) continue;

				Xv_val.push_back(Rp_data[j]);
				Xv_cIdx.push_back(N_c);
			}
		}
		/* Update row end */
		Xv_rEnd.push_back((int)Xv_cIdx.size());

		/* Update start point */
		i++;
		if (i < N_r)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	wsUint N = (wsUint)Xv_val.size();

	/* Init class and data */
	cSpsMatrix *Cp_ret = new cSpsMatrix(N_nr, N_nr, N);
	memcpy(Cp_ret->get()[0], get()[0], sizeof(wsReal)*N);

	memcpy(Cp_ret->getRowEnds(), getRowEnds(), sizeof(wsUint)*N_nr);
	memcpy(Cp_ret->getColIdxs(), getColIdxs(), sizeof(wsUint)*N);

	return *Cp_ret;
}

} // End namespace ONETOOL
