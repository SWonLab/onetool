#include <math.h>
#include "utils/matrices/std.h"
#include "utils/matrices/sym.h"
#include "utils/matrices/diag.h"
#include "utils/matrices/blk.h"
#include "global/Rconn.h"

namespace ONETOOL {

wsVec sseMpV(wsMat Ra_mat, wsUintCst N_r, wsUintCst N_c, wsVecCst Ra_v,
	wsRealCst R_mul/*=WISARD_NAN*/)
{
	/* Make return matrix */
	wsReal*	Ra_ret	= NULL;
#ifdef USE_SSE
	wsUint	N_med	= getMed(N_c);
#else
	wsUint	N_med	= 0;
#endif
	sseCalloc(Ra_ret, wsReal, N_r);

	if (NA(R_mul)) LOOP (i, N_r) {
		wsReal	R_sum	= W0;
		sse_t	sse_sum	= sseSet(0.0f);

#ifdef USE_SSE
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m = (sse_t *)(Ra_mat[i]+j);
			sse_t *sse_v = (sse_t *)(Ra_v+j);
			sse_sum = sseAdd(sse_sum, sseMul(*sse_m, *sse_v));
		}
		sseSum(sse_sum, R_sum);
#endif

		for (wsUint j=N_med ; j<N_c ; j++)
			R_sum += Ra_mat[i][j] * Ra_v[j];
		Ra_ret[i] = R_sum;
	} else LOOP (i, N_r) {
		wsReal	R_sum	= W0;
		sse_t	sse_sum	= sseSet(0.0f);

#ifdef USE_SSE
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m = (sse_t *)(Ra_mat[i]+j);
			sse_t *sse_v = (sse_t *)(Ra_v+j);
			sse_sum = sseAdd(sse_sum, sseMul(*sse_m, *sse_v));
		}
		sseSum(sse_sum, R_sum);
#endif

		for (wsUint j=N_med ; j<N_c ; j++)
			R_sum += Ra_mat[i][j] * Ra_v[j];
		Ra_ret[i] = R_sum * R_mul;
	}

	return Ra_ret;
}

wsReal** sseMcovCol(wsReal **Ra_mat, wsUint N_row, wsUint N_col)
{
	wsReal **Ra_R = sseMatrix(N_col, N_col);

	for (wsUint i=0 ; i<N_col ; i++) {
		for (wsUint j=i ; j<N_col ; j++) {
			Ra_R[i][j]	= W0;

			wsReal R_sumI = W0, R_sumJ = W0;
			// S <- S + (m[k,i]*m[k,j])
			for (wsUint k=0 ; k<N_row ; k++) {
				R_sumI		+= Ra_mat[k][i];
				R_sumJ		+= Ra_mat[k][j];
				Ra_R[i][j]	+= Ra_mat[k][j]*Ra_mat[k][i];
			}
			Ra_R[i][j]	= Ra_R[i][j]*(wsReal)N_row - R_sumI*R_sumJ;
			// ret[i,j] <- S/(d[1]-1)
			Ra_R[i][j]	/= (wsReal)(N_row*(N_row-1));
			Ra_R[j][i]	= Ra_R[i][j];
		}
	}

	return Ra_R;
}

wsReal** sym_sseCovRow(wsReal **Ra_mat, wsUint N_row, wsUint N_col,
					   char B_isComplete/*=0*/)
{
	return sseMcovRow(Ra_mat, N_row, N_col, B_isComplete, 1);
}

wsMat sseMcovRow(wsMat Ra_mat, wsUintCst N_row, wsUintCst N_col,
	char B_isComplete/*=0*/, char B_retSym/*=0*/)
{
	wsReal **Ra_R	= B_retSym ? sseSymMat(N_row) : sseMatrix(N_row, N_row);
	wsUint N_vMed = 0;
#ifdef USE_SSE
	N_vMed = getMed(N_col);
#endif
	wsUint i, j, k;

	if (B_isComplete) {
		for (i=0 ; i<N_row ; i++) {
			for (j=0 ; j<=i ; j++) {
				Ra_R[i][j]	= W0;

				wsReal N = (wsReal)N_col;
				wsReal sI = W0, sJ = W0;
				wsReal sIJ = W0;

				/* Get a_ijk */
				wsUint N_med = N_vMed;
#ifdef USE_SSE
				sse_t sse_sI = sseSet(0.0);
				sse_t sse_sJ = sseSet(0.0);
				sse_t sse_sIJ = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_mI = (sse_t *)(Ra_mat[i]+k);
					sse_t *sse_mJ = (sse_t *)(Ra_mat[j]+k);

					sse_sIJ = sseAdd(sse_sIJ, sseMul(*sse_mI, *sse_mJ));
					sse_sI = sseAdd(sse_sI, *sse_mI);
					sse_sJ = sseAdd(sse_sJ, *sse_mJ);
				}
				sseSum(sse_sI, sI);
				sseSum(sse_sJ, sJ);
				sseSum(sse_sIJ, sIJ);
#endif
				for (k=N_med ; k<N_col ; k++) {
					sI += Ra_mat[i][k];
					sJ += Ra_mat[j][k];
					sIJ += Ra_mat[i][k]*Ra_mat[j][k];
				}

				Ra_R[i][j] = (sIJ-sI/N*sJ-sJ/N*sI+sI*sJ/N) / (N-W1);
				if (!B_retSym)
					Ra_R[j][i] = Ra_R[i][j];
			}
		}
	} else {
		for (i=0 ; i<N_row ; i++) {
			for (j=0 ; j<=i ; j++) {
				Ra_R[i][j]	= W0;

				wsReal N = W0;
				wsReal sI = W0, sJ = W0;
				wsReal sIJ = W0;

				wsUint N_med = N_vMed;
				/* Get a_ijk */
#ifdef USE_SSE
				sse_t sse_sA = sseSet(0.0);
				sse_t sse_sI = sseSet(0.0);
				sse_t sse_sJ = sseSet(0.0);
				sse_t sse_sIJ = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t sse_A;
					sse_t *sse_mI = (sse_t *)(Ra_mat[i]+k);
					sse_t *sse_mJ = (sse_t *)(Ra_mat[j]+k);

					sse_t sse_tmp = sseSet(WISARD_NA_REAL);	/* tmp = MISSING */
					sse_A = sseAnd(sseNeq(*sse_mI, sse_tmp), sseNeq(*sse_mJ, sse_tmp));
					sse_tmp = sseSet(1.0);			/* tmp = 1 */
					sse_sA = sseAdd(sse_sA, sseAnd(sse_tmp, sse_A));

					sse_t sse_tI = sseAnd(*sse_mI, sse_A); // ti=mI if both 1
					sse_t sse_tJ = sseAnd(*sse_mJ, sse_A);
					sse_sIJ = sseAdd(sse_sIJ, sseMul(sse_tI, sse_tJ));
					sse_sI = sseAdd(sse_sI, sse_tI);
					sse_sJ = sseAdd(sse_sJ, sse_tJ);
				}
				sseSum(sse_sA, N);
				sseSum(sse_sI, sI);
				sseSum(sse_sJ, sJ);
				sseSum(sse_sIJ, sIJ);
#endif
				for (k=N_med ; k<N_col ; k++) {
					if (isAvailableReal(Ra_mat[i][k]) && isAvailableReal(Ra_mat[j][k])) {
						sI += Ra_mat[i][k];
						sJ += Ra_mat[j][k];
						sIJ += Ra_mat[i][k]*Ra_mat[j][k];
						N += W1;
					}
				}

				Ra_R[i][j] = (sIJ-sI/N*sJ-sJ/N*sI+sI*sJ/N) / (N-W1);
				if (!B_retSym)
					Ra_R[j][i] = Ra_R[i][j];
			}
		}
	}

	return Ra_R;
}

wsSym sseMcorRow(wsMat Ra_mat, wsUintCst N_r, wsUintCst N_c)
{
	wsVec	Ra_s	= sseVector(N_r);
	wsVec	Ra_ss	= sseVector(N_r);
	wsSym	Ra_r	= sseSymMat(N_r);

	/* Compute element-wise sum */
	LOOP (i, N_r)
		Ra_s[i] = sseVsum(Ra_mat[i], N_c, &(Ra_ss[i]));

	/* Compute corr coef */
	LOOP (i, N_r) {
		wsReal R_x	= Ra_s[i];
		wsReal R_x2	= Ra_ss[i];

		for (wsUint j=0 ; j<i ; j++) {
			wsReal	R_y		= Ra_s[j];
			wsReal	R_y2	= Ra_ss[j];

			wsReal	R_xy	= sseVV(Ra_mat[i], N_c, Ra_mat[j]);
			wsReal	num		= N_c*R_xy - R_x*R_y;
			wsReal	deno	= (N_c*R_x2 - SQR(R_x)) * (N_c*R_y2 - SQR(R_y));

			/* calculate correlation coefficient */
			Ra_r[i][j] = num / sqrt(deno);
		}
		Ra_r[i][i] = W1;
	}

	sseFree(Ra_s);
	sseFree(Ra_ss);
	return Ra_r;
}

wsVec sseMsumCol(wsMat Ra_mat, wsUintCst N_row, wsUint N_col/*=0xffffffff*/)
{
	wsVec Ra_ret = sseEmptyVec(N_row);

	/* Allocate default value for N_col as equal to N_row */
	if (N_col == 0xffffffff)
		N_col = N_row;

#ifdef USE_SSE
	wsUint N_med = getMed(N_col);
#endif

	for (wsUint i=0 ; i<N_row ; i++) {
		wsReal *Rp_m = Ra_mat[i];

#ifdef USE_SSE
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_r = (sse_t *)(Ra_ret + j);
			sse_t *sse_m = (sse_t *)(Rp_m + j);
			*sse_r = sseAdd(*sse_r, *sse_m);
		}
		for (wsUint j=N_med ; j<N_col ; j++)
			Ra_ret[j] += Rp_m[j];
#else
		for (wsUint j=0 ; j<N_col ; j++)
			Ra_ret[j] += Rp_m[j];
#endif
	}

	return Ra_ret;
}

wsReal* sseMtV(wsMatCst Ra_m, wsUintCst N_r1, wsUintCst N_c1, wsVecCst Ra_v, wsUintCst N_sz)
{
	if (N_r1!=N_sz && N_sz%N_r1 && N_r1%N_sz) 
		halt_fmt(WISARD_SYST_INVL_DIM, "matrix-vector multiplication", N_r1, N_sz);

	wsUint i, j;
	wsReal *Ra_r = NULL;
	sseCalloc(Ra_r, wsReal, N_c1);

	if (N_r1 < N_sz) {
		wsUint N_div = N_sz / N_r1;
		for (i=0 ; i<N_sz ; i++)  {
			wsReal	R_v		= Ra_v[i];
			wsUint	N_med	= 0;
			wsVecCst	Ra_mi	= Ra_m[(wsUint)(i/N_div)];
#ifdef USE_SSE
			N_med = getMed(N_c1);
			sse_t sse_v = sseSet(R_v);
			for (j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_r = (sse_t *)(Ra_r + j);
				sse_t *sse_m = (sse_t *)(Ra_mi + j);
				*sse_r = sseAdd(*sse_r, sseMul(*sse_m, sse_v));
			}
#endif
			for (j=N_med ; j<N_c1 ; j++)
				Ra_r[j] += Ra_mi[j] * R_v;
		}
	} else if (N_r1 > N_sz) {
		wsUint N_div = N_r1 / N_sz;
		for (i=0 ; i<N_r1 ; i++)  {
			wsReal	R_v		= Ra_v[(wsUint)(i/N_div)];
			wsUint	N_med	= 0;
			wsVecCst	Ra_mi	= Ra_m[i];
#ifdef USE_SSE
			N_med = getMed(N_c1);
			sse_t sse_v = sseSet(R_v);
			for (j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_r = (sse_t *)(Ra_r + j);
				sse_t *sse_m = (sse_t *)(Ra_mi + j);
				*sse_r = sseAdd(*sse_r, sseMul(*sse_m, sse_v));
			}
#endif
			for (j=N_med ; j<N_c1 ; j++)
				Ra_r[j] += Ra_mi[j] * R_v;
		}
	} else for (i=0 ; i<N_sz ; i++)  {
		wsReal	R_v		= Ra_v[i];
		wsUint	N_med	= 0;
		wsVecCst	Ra_mi	= Ra_m[i];
#ifdef USE_SSE
		N_med = getMed(N_c1);
		sse_t sse_v = sseSet(R_v);
		for (j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_r = (sse_t *)(Ra_r + j);
			sse_t *sse_m = (sse_t *)(Ra_mi + j);
			*sse_r = sseAdd(*sse_r, sseMul(*sse_m, sse_v));
		}
#endif
		for (j=N_med ; j<N_c1 ; j++)
			Ra_r[j] += Ra_mi[j] * R_v;
	}

	return Ra_r;
}

wsMat sseMsubsetRect(wsMat Ra_m, wsUintCst N_sz, const char *Ba_filt,
	char N_incVal, wsUintCst N_szAfter)
{
	wsReal **Ra_ret = sseMatrix(N_szAfter, N_szAfter);

	wsUint i, j, _i, _j;
	for (i=_i=0 ; i<N_sz ; i++) {
		if (Ba_filt[i] != N_incVal) continue;

		for (j=_j=0 ; j<N_sz ; j++) {
			if (Ba_filt[j] != N_incVal) continue;

			Ra_ret[_i][_j] = Ra_m[i][j];
			_j++;
		}
		_i++;
	}
	if (_i != N_szAfter)
		halt("Invalid subset making, desired [%d] but subsetted [%d]",
		N_szAfter, _i);

	return Ra_ret;
}

wsMat sseMsubsetRect(wsMat Ra_m, wsUintCst N_sz, vBool &Ba_filt,
	char N_incVal, wsUintCst N_szAfter)
{
	wsMat Ra_ret = sseMatrix(N_szAfter, N_szAfter);

	wsUint i, j, _i, _j;
	for (i=_i=0 ; i<N_sz ; i++) {
		if (Ba_filt[i] != N_incVal) continue;

		for (j=_j=0 ; j<N_sz ; j++) {
			if (Ba_filt[j] != N_incVal) continue;

			Ra_ret[_i][_j] = Ra_m[i][j];
			_j++;
		}
		_i++;
	}
	if (_i != N_szAfter)
		halt("Invalid subset making, desired [%d] but subsetted [%d]",
		N_szAfter, _i);

	return Ra_ret;
}

cStdMatrix::cStdMatrix()
{
	__init(NULL, 0, 0, 1);
}

cStdMatrix::cStdMatrix(cStdMatrix &&M)
{
	rem();
	X_actClean = M.getClean();
	__init(M.get(), M.row(), M.col());

	M.setClean();
}

cStdMatrix::cStdMatrix(cSpsMatrix &P)
{
	wsUint	N_r		= P.row();
	wsUint	N_c		= P.col();
	wsUint*	Na_rEnd	= P.getRowEnds();
	wsUint*	Na_cIdx	= P.getColIdxs();
	__init(NULL, N_r, N_c, 1);

	/* Allocate */
	Ra_data = sseEmptyMatrix(N_r, N_c);

	/* Move values from sparse matrix */
	wsUint	i		= 0;
	wsReal*	Rp_data	= P.get()[0];
	wsUint	N_s		= 0;
	wsUint	N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e)
			for (wsUint j=N_s ; j<N_e ; j++)
				Ra_data[i][Na_cIdx[j]] = Rp_data[j];

		/* Update start point */
		i++;
		N_s = N_e;
	} while (i < N_r);
}

cStdMatrix& cStdMatrix::operator=(cStdMatrix &&M)
{
	rem();
	X_actClean	= M.getClean();
	Ra_data		= M.get();
	N_row		= M.row();
	N_col		= M.col();

	M.setClean();
	return *this;
}

void cStdMatrix::init(wsUint N_R, wsUint N_C, wsReal **Rp_data/*=NULL*/)
{
	init(MATOPT_DEFAULT, N_R, N_C, 0, Rp_data);
}

void cStdMatrix::init(xMatOpt X_flag, wsUint N_R, wsUint N_C, wsReal R_val,
	wsReal **Rp_data/*=NULL*/)
{
	/* Remove previously allocated space */
	if (Ra_data && !X_actClean) deallocMatrix(get(), row(), (void *)1);

	/* If transposed, exchange row/col */
	if (X_flag & MATOPT_TRANSPOSE) {
		/* Data must exists */
		if (Rp_data == NULL)
			halt_fmt(WISARD_SYST_NULL_DATA, "matrix to be transposed");
		wsUint N_tmp = N_R;
		N_R = N_C;
		N_C = N_tmp;
	}
	/* Is this matrix should be deallocated after the end of this class? */
	if (X_flag&MATOPT_DELROW)
		X_actClean	= MATDEL_ROW;
	else if (X_flag&MATOPT_DELNONE)
		X_actClean	= MATDEL_NONE;
	else
		X_actClean	= MATDEL_ALL;

	/* With given data it cannot be randomized */
	if ((X_flag&MATOPT_RANDOMIZE) && Rp_data) raise(MATERR_PARAM_INCORRECT);

	/* Allocate memory or set the pointer */
	Ra_data	= Rp_data==NULL || X_flag&MATOPT_TRANSPOSE ?
		sseEmptyMatrix(N_R, N_C) : Rp_data;
	N_row	= N_R;
	N_col	= N_C;

	/* If the transposition is activated */
	if (X_flag&MATOPT_TRANSPOSE)
		for (wsUint i=0 ; i<N_row ; i++)
			for (wsUint j=0 ; j<N_col ; j++)
				Ra_data[i][j] = Rp_data[j][i];

	/* If the randomization is activated */
	if (X_flag&MATOPT_RANDOMIZE)
		for (wsUint i=0 ; i<N_row ; i++)
			for (wsUint j=0 ; j<N_col ; j++)
				Ra_data[i][j] = (wsReal)(rand()%1000)-REAL_CONST(500.0);
	/* If the initialization is activated */
	else if (X_flag&MATOPT_UNIFINIT)
		for (wsUint i=0 ; i<N_row ; i++)
			sseVinit(Ra_data[i], N_col, R_val);
}

cStdMatrix::cStdMatrix(wsUint N_R, wsUint N_C, char B_rand/*=0*/)
{
	Ra_data = NULL;
	init(B_rand?MATOPT_RANDOMIZE:MATOPT_DEFAULT, N_R, N_C, 0);
}

cStdMatrix::cStdMatrix(wsUint N_R, wsUint N_C, wsReal R_val)
{
	Ra_data = NULL;
	init(MATOPT_UNIFINIT, N_R, N_C, R_val);
}

cStdMatrix::cStdMatrix(wsUint N_R, wsUint N_C, wsReal **Rp_data/*=NULL*/,
	xMatCleanAct X_cleanAct/*=0*/)
{
	Ra_data = NULL;
	init(X_cleanAct==MATDEL_NONE?MATOPT_DELNONE:
		(X_cleanAct==MATDEL_ROW?MATOPT_DELROW:MATOPT_DEFAULT), N_R, N_C, 0, Rp_data);
}

cStdMatrix::cStdMatrix(xMatOpt X_opt, wsUint N_R, wsUint N_C,
	wsReal **Rp_data/*=NULL*/)
{
	Ra_data = NULL;
	init(X_opt, N_R, N_C, 0, Rp_data);
}

cStdMatrix::~cStdMatrix()
{
	rem();
}

/*
 *
 * Addition
 * 
 */

cStdMatrix& cStdMatrix::operator+=(cStdMatrix &C)
{
	/* Check */
	if (col()!=C.col() || row()!=C.row()) raise(MATERR_DIM_INCORRECT);
	/* Do */
	sseMaM(get(), C.get(), get(), row(), col());
	return *this;
}

cStdMatrix& cStdMatrix::operator+=(cSymMatrix &C)
{
	/* Check */
	if (!sqr() || !C.sqr() || col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	/* Do */
	sseMaS(cget(), C.get(), get(), row());
	return *this;
}

cStdMatrix& cStdMatrix::operator+=(cDiagMatrix &C)
{
	/* Check */
	if (!sqr() || !C.sqr() || col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	/* Do */
	wsReal **T	= get();
	wsReal *O	= C.get()[0];
	for (wsUint i=0,N=col() ; i<N ; i++) T[i][i] += O[i];
	return *this;
}

cStdMatrix& cStdMatrix::operator+=(cBlkMatrix &C)
{
	/* Check */
	if (!sqr() || col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	for (wsUint i=0,N=col() ; i<N ; i+=C.blk())
		sseMaM(get(), C.get(), get(), i, i, C.blk());
	return *this;
}

cStdMatrix& cStdMatrix::operator+=(wsReal R)
{
	sseMaC(get(), row(), col(), R, get());
	return *this;
}

/*
 *
 * Subtraction
 * 
 */

cStdMatrix& cStdMatrix::operator-=(cStdMatrix &C)
{
	/* Check */
	if (col()!=C.col() || row()!=C.row()) raise(MATERR_DIM_INCORRECT);
	/* Do */
	sseMsM(get(), C.get(), get(), row(), col());
	return *this;
}

cStdMatrix& cStdMatrix::operator-=(cSymMatrix &C)
{
	/* Check */
	if (!sqr() || !C.sqr() || col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	/* Do */
	sseMsS(cget(), C.cget(), get(), row());
	return *this;
}

cStdMatrix& cStdMatrix::operator-=(cDiagMatrix &C)
{
	/* Check */
	if (!sqr() || !C.sqr() || col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	/* Do */
	wsReal **T	= get();
	wsReal *O	= C.get()[0];
	for (wsUint i=0,N=col() ; i<N ; i++) T[i][i] -= O[i];
	return *this;
}

cStdMatrix& cStdMatrix::operator-=(cBlkMatrix &C)
{
	/* Check */
	if (!sqr() || col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	for (wsUint i=0,N=col() ; i<N ; i+=C.blk())
		sseMsM(get(), C.get(), get(), i, i, C.blk());
	return *this;
}

cStdMatrix& cStdMatrix::operator-=(wsReal R)
{
	sseMsC(cget(), R, get(), row(), col());
	return *this;
}

cStdMatrix cStdMatrix::operator-(const cStdMatrix &C)
{
	/* Check */
	if (col()!=C.col() || row()!=C.row()) raise(MATERR_DIM_INCORRECT);
	wsMat	Ra_ret = sseMatrix(row(), col());
	/* Do */
	sseMsM(get(), C.get(), Ra_ret, row(), col());
	return cStdMatrix(row(), col(), Ra_ret);
}

/*
 *
 * Multiplication
 *
 */

cStdMatrix& cStdMatrix::operator*=(cStdMatrix &C)
{
	/* Check */
	if (col()!=C.row()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseMpM(get(), row(), col(), C.get(), C.row(), C.col());
	/* Substitute matrix */
	init(row(), C.col(), Ra_ret);
	return *this;
}

cStdMatrix& cStdMatrix::operator*=(cSymMatrix &C)
{
	/* Check */
	if (col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseMpS(get(), row(), col(), C.get(), col());
	/* Substitute matrix */
	init(row(), C.col(), Ra_ret);
	return *this;
}

cStdMatrix& cStdMatrix::operator*=(cDiagMatrix &C)
{
	/* Check */
	if (!sqr() || col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseMpD(get(), row(), col(), C.get()[0]);
	/* Substitute matrix */
	init(row(), C.col(), Ra_ret);
	return *this;
}

cStdMatrix& cStdMatrix::operator*=(cBlkMatrix &C)
{
	/* Check */
	if (row()!=C.row() || col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseMpB(get(), row(), col(), C.get(), C.blk(),
		C.rep());
	/* Substitute matrix */
	init(row(), C.col(), Ra_ret);
	return *this;
}

cStdMatrix cStdMatrix::operator*(const cStdMatrix &C)
{
	/* Check */
	if (col()!=C.row()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseMpM(get(), row(), col(), C.get(), C.row(), C.col());
	/* Return */
	return cStdMatrix(row(), C.col(), Ra_ret);
}

cStdMatrix cStdMatrix::operator*(const cSymMatrix &C)
{
	/* Check */
	if (col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseMpS(get(), row(), col(), C.get(), col());
	/* Return */
	return cStdMatrix(row(), C.col(), Ra_ret);
}

cStdMatrix cStdMatrix::operator*(const cDiagMatrix &C)
{
	/* Check */
	if (col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseMpD(get(), row(), col(), C.get()[0]);
	/* Return */
	return cStdMatrix(row(), C.col(), Ra_ret);
}

cStdMatrix cStdMatrix::operator*(const cBlkMatrix &C)
{
	/* Check */
	if (row()!=C.row() || col()!=C.col()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseMpB(get(), row(), col(), C.get(), C.blk(),
		C.rep());
	/* Return */
	return cStdMatrix(row(), C.col(), Ra_ret);
}

cStdMatrix& cStdMatrix::operator*=(wsReal R)
{
	sseMpC(cget(), R, get(), row(), col());
	return *this;
}

cStdMatrix cStdMatrix::transpose()
{
	wsReal **Ra_t = ONETOOL::transpose(get(), row(), col());
	return cStdMatrix(col(), row(), Ra_t);
}

cStdMatrix* cStdMatrix::transposeNew()
{
	wsReal **Ra_t = ONETOOL::transpose(get(), row(), col());
	return new cStdMatrix(col(), row(), Ra_t);
}

cSymMatrix cStdMatrix::MMt(const cStdMatrix &M)
{
	/* [N*Q] %*% [Q*Q] %*% [Q*N] so
	 *  M should be square matrix
	 *  row(M) should be equal to col(this) */
	if (col()!=M.row() || !M.sqr()) raise(MATERR_DIM_INCORRECT);
	wsSym Ra_ret = sseMMMt(get(), row(), col(), M.get(), M.row(), M.col());

	return cSymMatrix(Ra_ret, row());
}

cStdMatrix cStdMatrix::MMt(const cStdMatrix &M, const cStdMatrix &N)
{
	/* [N*Q] %*% [Q*R] %*% [R*S] so
	 *  row(M) should be equal to col(this)
	 *  row(N) should be equal to col(M) */
	if (col()!=M.row() || M.col()!=N.row()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseMMMt(get(), row(), col(), M.get(), M.row(), M.col());

	return cStdMatrix(row(), N.col(), Ra_ret);
}

cStdMatrix cStdMatrix::MMt(const cSymMatrix &S, const cStdMatrix &M)
{
	if (col()!=S.row() || S.col()!=M.row()) raise(MATERR_DIM_INCORRECT);

	/* [N*Q] %*% [Q*R] %*% [R*S] so
	 *  row(M) should be equal to col(this)
	 *  row(N) should be equal to col(M) */
	wsMat	Ra_r1	= sseMpS(get(), row(), col(), S.get(), S.row());
	wsMat	Ra_ret	= sseMpMt(Ra_r1, row(), col(), M.get(), M.row(), M.col());

	return cStdMatrix(row(), M.col(), Ra_ret);
}

cStdMatrix cStdMatrix::MMt(const cSymMatrix &S, cSpsMatrix &P)
{
	if (col()!=S.row() || S.col()!=P.row()) raise(MATERR_DIM_INCORRECT);

	/* [N*Q] %*% [Q*R] %*% [R*S] so
	 *  row(M) should be equal to col(this)
	 *  row(N) should be equal to col(M) */
	wsMat	Ra_r1	= sseMpS(get(), row(), col(), S.get(), S.row());
	wsMat	Ra_ret	= multMP(Ra_r1, row(), col(), P.get()[0], P.row(),
		P.col(), P.getRowEnds(), P.getColIdxs());

	return cStdMatrix(row(), P.col(), Ra_ret);
}

cStdMatrix cStdMatrix::MV(const cStdMatrix &M, const cVector &V)
{
	/* [N*Q] %*% [Q*R] %*% [R*S] so
	 *  row(M) should be equal to col(this)
	 *  row(N) should be equal to col(M) */
	if (col()!=M.row() || M.col()!=V.size()) raise(MATERR_DIM_INCORRECT);
	wsReal *Ra_v = V.get();
	wsReal **Ra_ret = sseMMMt(get(), row(), col(), M.get(), M.row(), M.col(),
		&Ra_v, 1, V.size());

	return cStdMatrix(row(), V.size(), Ra_ret);
}

cSymMatrix cStdMatrix::MMt(cMatrix &M)
{
	switch (M.getType()) {
	case MATTYPE_SYM:
		return MMt(*(cSymMatrix *)&M);
	case MATTYPE_IDT:
		return Mt();
	case MATTYPE_STD:
		return MMt(*(cStdMatrix *)&M);
	case MATTYPE_BLK:
		return MMt(*(cBlkMatrix *)&M);
	default:
		raise(MATERR_TYPE_INCORRECT);
	}
	/* Should not reach here */
	return cSymMatrix();
}

cSymMatrix cStdMatrix::MMt(const cBlkMatrix &M)
{
	//wsReal	**Ra_int	= sseMpB(get(), row(), col(), M.get(), M.blk(), M.rep());
	wsReal **Ra_ret	= NULL;
	if (col() == M.rep()) {
		Ra_ret = sseSymMat(row());
		wsReal	R_bSum	= sseMsum(M.get(), M.blk());
		for (wsUint i=0,N=row() ; i<N ; i++)
			for (wsUint j=0 ; j<=i ; j++)
				Ra_ret[i][j] = R_bSum * sseVV(get()[i], col(), get()[j]);
	} else {
		wsReal **Ra_int = sseMpB(get(), row(), col(),
			M.get(), M.blk(), M.rep());
		Ra_ret = sym_sseMpMt(Ra_int, row(), col(), get(), row(), col());
		//		wsReal	R_bSum	= M.sum();
		//		for (wsUint i=0,N=row() ; i<N ; i++)
		//			for (wsUint j=0 ; j<=i ; j++)
		//				Ra_ret[i][j] = R_bSum * sseVV(get()[i], col(), get()[j]);
		deallocMatrix(Ra_int, row(), (void *)1);
	}

	return cSymMatrix(Ra_ret, row());
}

cSymMatrix cStdMatrix::MMt(const cSymMatrix &M)
{
	wsReal **Ra_int	= sseMpS(get(), row(), col(), M.get(), M.row());
	wsReal **Ra_sym = sym_sseMpMt(Ra_int, row(), M.col(), get(), row(),
		col());
	deallocMatrix(Ra_int, row(), (void *)1);

	return cSymMatrix(Ra_sym, row());
}

/*
mdmt <- function(m, d) {
	L <- dim(m)[1]
	ret <- list()
	ret$expected <- m %*% diag(d) %*% t(m)
	ret$my <- matrix(0, nrow=L, ncol=L)
	for (i in 1:L) {
		intm <- m[i,] * d
		for (j in 1:L)
			ret$my[i,j] <- sum(intm * m[j,])
	}
	return(ret)
}
m <- matrix(rnorm(80), nrow=4)
d <- rnorm(20)
r <- mdmt(m, d)
cat("Diff : ", sum((r$expected-r$my)^2), "\n")
 */
cSymMatrix cStdMatrix::MMt(const cDiagMatrix &M)
{
	/* R_ij = X_i+ * D_+j = X_ij * D_jj */
	/* RR_ij = R_i+ * X_j+ = X_i+ * D_++ */
	if (col()!=M.row()) raise(MATERR_DIM_INCORRECT);

	wsUint Ro = row();
	wsReal *TE = NULL;
	wsReal **S = get();
	wsReal *D = M.get()[0];
	wsUint N = M.row();
	wsReal **R = sseSymMat(Ro);
	sseMalloc(TE, wsReal, N);

	for (wsUint i=0 ; i<Ro ; i++) {
		/* Calc X_i+*D_++ */
		sseVpV(S[i], D, TE, N);
		/* Calc R_i+ * X_j+ */
		for (wsUint j=0 ; j<=i ; j++)
			R[i][j] = sseVV(TE, N, S[j]);
	}

	sseFree(TE);

	return cSymMatrix(R, Ro);
}

cStdMatrix cStdMatrix::MMt(const cDiagMatrix &D, const cStdMatrix &M)
{
	wsMat Ra_ret = sseMDMt(get(), row(), col(), D.get()[0], M.get(), M.row(),
		M.col());
	return cStdMatrix(row(), M.row(), Ra_ret);
}

cSymMatrix& cStdMatrix::MMt(const cVector &V)
{
	if (col() != V.size()) raise(MATERR_DIM_INCORRECT);
	wsSym Ra_ret = sseMDMt(get(), row(), col(), V.get());
	return *(new cSymMatrix(Ra_ret, row()));
}

cVector cStdMatrix::MV(cDiagMatrix &M, cVector &V)
{
	/* R_ij = X_i+ * D_+j = X_ij * D_jj */
	/* RR_ij = R_i+ * X_j+ = X_i+ * D_++ */
	if (col()!=M.row()) raise(MATERR_DIM_INCORRECT);

	wsUint Ro = row();
	wsReal *TE = NULL;
	wsReal **S = get();
	wsReal *D = M.get()[0];
	wsUint N = M.row();
	wsReal *R = sseVector(Ro);
	sseMalloc(TE, wsReal, N);

	for (wsUint i=0 ; i<Ro ; i++) {
		/* Calc X_i+*D_++ */
		sseVpV(S[i], D, TE, N);
		/* Calc R_i+ * X_j+ */
		R[i] = sseVV(TE, N, V.get());
	}

	sseFree(TE);

	return cVector(R, Ro);
}

cVector cStdMatrix::MV(cSymMatrix &M, cVector &V)
{
	if (col()!=M.row() || M.row()!=V.size()) raise(MATERR_DIM_INCORRECT);

	wsMat	intm	= sseMpS(get(), row(), col(), M.get(), M.row());
	wsReal	*res	= sseMpV(intm, row(), col(), V.get());
	deallocMatrix(intm, row(), (void *)1);

	return cVector(res, row());
}

cStdMatrix cStdMatrix::Mt(const cStdMatrix &M)
{
	/* [N*Q] %*% [P*Q] = [N*P] so
	 *  col(M) should be equal to col(this) */
	if (col()!=M.col()) raise(MATERR_DIM_INCORRECT);
	return cStdMatrix(row(), M.row(),
		sseMpMt(get(), row(), col(), M.get(), M.row(), M.col()));
}

cSymMatrix cStdMatrix::Mts(const cStdMatrix &M)
{
	/* [N*Q] %*% [P*Q] = [N*P] so
	 *  col(M) should be equal to col(this) */
	if (col()!=M.col()) raise(MATERR_DIM_INCORRECT);
	return cSymMatrix(
		sym_sseMpMt(get(), row(), col(), M.get(), M.row(), M.col()),
		row());
}

cSymMatrix cStdMatrix::Mt()
{
	/* Self-multiplication, MM' type */
	/* [N*Q] %*% [Q*N] = [N*N] */
	wsReal **Ra_ret = sym_sseMpMt(get(), row(), col());

	return cSymMatrix(Ra_ret, row());
}

cSymMatrix cStdMatrix::Ms(const cStdMatrix &M)
{
	wsReal **Ra_ret = sym_sseMpM(get(), row(), col(), M.get(), M.row(),
		M.col());

	return cSymMatrix(Ra_ret, row());
}

cVector cStdMatrix::tV(wsReal *V, wsUint L)
{
	if (row()!=L) raise(MATERR_DIM_INCORRECT);
	wsReal *Ra_ret = sseMtV(cget(), row(), col(), V, L);
	return cVector(Ra_ret, col());
}

cVector cStdMatrix::tV(cVector &V)
{
	if (row()!=V.size()) raise(MATERR_DIM_INCORRECT);
	wsReal *Ra_ret = sseMtV(cget(), row(), col(), V.get(), V.size());
	return cVector(Ra_ret, col());
}

cStdMatrix cStdMatrix::tM(cStdMatrix &M)
{
	if (row()!=M.row()) raise(MATERR_DIM_INCORRECT);
	wsMat Ra_ret = sseMtM(get(), row(), col(), M.get(), M.row(), M.col());
	return cStdMatrix(col(), M.col(), Ra_ret);
}

cStdMatrix cStdMatrix::tM(cSymMatrix &S)
{
	if (row()!=S.row()) raise(MATERR_DIM_INCORRECT);
	wsMat Ra_ret = sseMtS(get(), row(), col(), S.get(), S.row());
	return cStdMatrix(col(), S.col(), Ra_ret);
}

cStdMatrix cStdMatrix::tM(cDiagMatrix &D)
{
	if (row()!=D.row()) raise(MATERR_DIM_INCORRECT);
	wsUint R = row();
	wsUint C = col();
	wsMat S = get();
	wsReal *V = D.get()[0];
	wsMat Ra_ret = sseMatrix(C, R);

	for (wsUint i=0 ; i<C ; i++) {
		wsReal *T = Ra_ret[i];
		for (wsUint j=0 ; j<R ; j++)
			T[j] = V[j] * S[j][i];
	}
	return cStdMatrix(C, R, Ra_ret);
}

cSymMatrix cStdMatrix::tM(wsRealCst R_mul/*=W0*/)
{
	wsSym Ra_ret = sseMtM(get(), row(), col(), R_mul);

	return cSymMatrix(Ra_ret, col());
}

/*
 *
 * Special operations definition
 *
 */

wsRealCst cStdMatrix::normF(cStdMatrix &M)
{
	return compMMfrobenius(get(), row(), col(), M.get(), M.row(), M.col());
}

wsRealCst cStdMatrix::normF()
{
	return normFrobenius(get(), row(), col());
}

cSymMatrix cStdMatrix::proj()
{
	/* Do I - X (X'X)^-1 X' */
	cSymMatrix	C_XtX	= Mt();
	cSymMatrix	&C_XtXi	= C_XtX.inv();
	C_XtX.rem();
	cStdMatrix	C_XXtXi	= tM(C_XtXi);
	C_XtXi.rem();
	cSymMatrix	C_XXXiX	= C_XXtXi.Ms(*this);
	C_XXtXi.rem();
	cIdtMatrix	C_1n(C_XXXiX.row());

	return C_1n - C_XXXiX;
}

cStdMatrix& cStdMatrix::normalize(wsRealCst R_mul/*=W1*/)
{
	wsUint	N = row();
	wsUint	C = col();
	wsMat	R = get();
	LOOP (i, N) {
		wsVec	V		= R[i];
		wsRealCst	R_mean	= sseVsum(V, C) / C;
		wsRealCst	R_sd	= ::sqrt(sseVsCsum(V, R_mean, C, 1) / (wsReal)(C-1));
		sseVsC(V, R_mean, V, C, R_mul/R_sd);
	}
	return *this;
}

cStdMatrix&	cStdMatrix::operator=(cStdMatrix &M)
{
	init(M.row(), M.col(), M.get());
	return *this;	
}

void cStdMatrix::subIself()
{
	wsUint	N = row();
	wsUint	C = col();
	wsMat	R = get();
	LOOP (i, N) {
		sseVneg(R[i], C, R[i]);
		R[i][i] += W1;
	}
}

cStdMatrix cStdMatrix::subI()
{
	wsUint	N = row();
	wsUint	C = col();
	wsMat	V = get();
	wsMat	R = sseMatrix(N, C);	
	LOOP (i, N) {
		sseVneg(V[i], C, R[i]);
		R[i][i] += W1;
	}
	return cStdMatrix(N, C, R);
}

cStdMatrix cStdMatrix::divide(wsUint* Na_idx, wsUint N_szIdx, cStdMatrix* Mp_rest/*=NULL*/)
{
	wsUint		N		= row();
	wsUint		C		= col();
	wsMat		Ra_ret	= sseMatrix(N_szIdx, C);
	char*		Ba_map	= NULL;
	wsCalloc(Ba_map, char, N);
	LOOP (i, N_szIdx)
		Ba_map[Na_idx[i]] = 1;
#ifdef _DEBUG
	wsUint		Nz		= 0;
	LOOP (i, N) if (Ba_map[i]) Nz++;
	if (Nz != N_szIdx) halt("SYSERR: Divide error, duplication found");
#endif
	wsMat		V		= get();
	wsUint		j		= 0;
	wsUint		k		= 0;

	if (Mp_rest) {
		Mp_rest->init(N - N_szIdx, C);
		wsMat	Ra_rest = Mp_rest->get();
		LOOP (i, N) {
			if (Ba_map[i] == 1)
				memcpy(Ra_ret[j++], V[i], sizeof(wsReal)*C);
			else
				memcpy(Ra_rest[k++], V[i], sizeof(wsReal)*C);
		}
	} else LOOP(i, N)
		if (Ba_map[i] == 1) memcpy(Ra_ret[j++], V[i], sizeof(wsReal)*C);

	return cStdMatrix(N_szIdx, C, Ra_ret);
}

cStdMatrix cStdMatrix::divideCol(wsUint* Na_idx, wsUint N_szIdx, cStdMatrix* Mp_rest/*=NULL*/)
{
	wsUint		N		= row();
	wsUint		C		= col();
	wsMat		Ra_ret	= sseMatrix(N, N_szIdx);
	char*		Ba_map	= NULL;
	wsCalloc(Ba_map, char, C);
	LOOP(i, N_szIdx)
		Ba_map[Na_idx[i]] = 1;
	wsMat		V		= get();

	if (Mp_rest) {
		Mp_rest->init(N, C - N_szIdx);
		wsMat	Ra_rest = Mp_rest->get();
		LOOP (i, N) {
			wsUint j = 0, k = 0;
			LOOP (l, C)
				if (Ba_map[l] == 1)
					Ra_ret[i][j++] = V[i][l];
				else
					Ra_rest[i][k++] = V[i][l];
		}
	} else LOOP (i, N) {
		wsUint j = 0;
		LOOP (l, C) if (Ba_map[l] == 1)
			Ra_ret[i][j++] = V[i][l];
	}

	return cStdMatrix(N, N_szIdx, Ra_ret);
}

void cStdMatrix::trcLess(wsRealCst R_val, wsUint N_from/*=0*/)
{
	wsUint		N		= row();
	wsUint		C		= col();
	wsMat		V		= get();

	LOOP (i, N) for (wsUint j=N_from ; j<C ; j++)
		if (V[i][j] < R_val) V[i][j] = W0;
}

cStdMatrix cStdMatrix::log1exp()
{
	wsUint		N		= row();
	wsUint		C		= col();
	wsMat		V		= get();
	wsMat		Ra_ret	= sseMatrix(N, C);

	LOOP(i, N) LOOP(j, C)
		Ra_ret[i][j] = log(W1 + exp(V[i][j]));

	return cStdMatrix(N, C, Ra_ret);
}

cVector cStdMatrix::r2v(int N_idx)
{
	wsUint N_finIdx = 0xffffffff;

	/* Do row-N_idx if negative */
	if (N_idx < 0)
		N_finIdx = row() + N_idx;
	else
		N_finIdx = N_idx;

	/* Range check */
	if (N_finIdx == 0xffffffff || N_finIdx >= row())
		halt("Too large/small row index[%d] requested", N_idx);

	/* Make ret vec */
	wsReal	*Rp_ret = sseVector(col());
	memcpy(Rp_ret, get()[N_finIdx], sizeof(wsReal)*col());

	return cVector(Rp_ret, col());
}

cVector cStdMatrix::r2v_ptr(int N_idx)
{
	wsUint N_finIdx = 0xffffffff;

	/* Do row-N_idx if negative */
	if (N_idx < 0)
		N_finIdx = row() + N_idx;
	else
		N_finIdx = N_idx;

	/* Range check */
	if (N_finIdx == 0xffffffff || N_finIdx >= row())
		halt("Too large/small row index[%d] requested", N_idx);

	/* Make ret vec */
	return cVector(get()[N_finIdx], col(), 1);
}

wsRealCst cStdMatrix::tr2()
{
	wsReal R_ret = W0;
	LOOP (i, row())
		R_ret += sseVV(Ra_data[i], col());
	return R_ret;
//	return diagSum2matrix(get(), row(), col(), get(), row(), col());
}

cVector cStdMatrix::rSS()
{
	wsReal *Ra_ret = sseMssR(get(), row(), col());
	return cVector(Ra_ret, row());
}

cStdMatrix cStdMatrix::rDiv(const cVector &V_divider)
{
	if (row()!=V_divider.size()) raise(MATERR_DIM_INCORRECT);
	wsReal	*Ra_div		= V_divider.get();
	wsReal	**Ra_ret	= sseMatrix(row(), col());
	for (wsUint i=0 ; i<row() ; i++)
		sseVpC(get()[i], W1/Ra_div[i], Ra_ret[i], col());
	return cStdMatrix(row(), col(), Ra_ret);
}

void cStdMatrix::setRow(wsUint N_row, wsReal R_val)
{
	if (row() <= N_row) halt("Row index[%d] larger than its size[%d] called", N_row+1, row());
	sseVinit(get()[N_row], col(), R_val);
}

void cStdMatrix::setCol(wsUint N_col, wsReal R_val)
{
	wsUint	N_row	= row();
	wsMat	Ra_data	= get();
	if (col() <= N_col) halt("Column index[%d] larger than its size[%d] called", N_col+1, col());
	/* Since there is no way to set the value of columns in SIMD way, just to it */
	for (wsUint i=0 ; i<N_row ; i++)
		Ra_data[i][N_col] = R_val;
}

void cStdMatrix::addCol(wsUint N_col, wsReal R_val)
{
	wsUint	N_row	= row();
	wsMat	Ra_data	= get();
	if (col() <= N_col) halt("Column index[%d] larger than its size[%d] called", N_col+1, col());
	/* Since there is no way to set the value of columns in SIMD way, just to it */
	for (wsUint i=0 ; i<N_row ; i++)
		Ra_data[i][N_col] += R_val;
}

cVector cStdMatrix::vectorize()
{
	wsUint R = row();
	wsUint C = col();
	wsReal *Ra_ret = sseVector(R*C);

	for (wsUint i=0 ; i<R ; i++) {
		wsUint K=i;
		wsReal *S = get()[i];
		for (wsUint j=0; j<C ; j++,K+=R)
			Ra_ret[K] = S[j];
	}
	return cVector(Ra_ret, R*C);
}

wsRealCst cStdMatrix::trMSMt(cSymMatrix &S)
{
	wsReal	R	= W0;
	wsMat	A	= get();
//	SYM_t	B	= S.get();
	wsUint	p	= S.row();
	wsUint	n	= col();
	wsUint	m	= row();

	wsMat U = sseMpS(A, m, n, S.get(), p);
	for (wsUint i=0 ; i<m ; i++)
		R += sseVV(A[i], n, U[i]);
	return R;
}

wsRealCst cStdMatrix::trMtSM(cSymMatrix &S)
{
	wsReal	R	= W0;
	wsMat	A	= get();
	wsSym	B	= S.get();
	wsUint	p	= S.row();
	wsUint	n	= col();
/*

a <- Xt
b <- S
r <- 0
m <- c()
for (i in 1:n) { #x */
	for (wsUint i=0 ; i<n ; i++) {
		wsReal *Ra_Ak_ = sseVector(p);
		wsReal *Ra_ij = sseEmptyVec(p);
		for (wsUint k=0 ; k<p ; k++)
			Ra_Ak_[k] = A[k][i];
		for (wsUint k=0 ; k<p ; k++) {
/*
r.ij <- rep(0, p)
a.k. <- a[,i]
for (k in 1:p) {
a.ki <- a[k,i]
*/
			Ra_ij[k] += sseVV(Ra_Ak_, k+1, B[k]);
			sseVaV(Ra_ij, B[k], Ra_ij, k, Ra_Ak_[k]);
		}
		R += sseVV(Ra_ij, p, Ra_Ak_);
		sseFree(Ra_Ak_);
		sseFree(Ra_ij);
	}
/*
 *# j
if (k>1) for (j in 1:(k-1)) {
r.ij[k] <- r.ij[k] + a.k.[j] * b[k,j]
r.ij[j] <- r.ij[j] + a.ki * b[k,j]
#     for (j in 1:p)
#      r.ij[j] <- r.ij[j] + a.ki * b[k,j]
}
r.ij[k] <- r.ij[k] + a.ki * b[k,k]
}
cat(r.ij, "\n")
m[i] <- sum(r.ij * a.k.)
r <- r + m[i]
#  for (j in 1:p)
#    for (k in 1:p)
#      t <- t + a[k,i]*b[k,j]
#    r <- t * a[j,i]
}

*/
	return R;
}

cVector cStdMatrix::diag(cSymMatrix &S)
{
	wsMat	A	= get();
	wsSym	B	= S.get();
	wsUint	p	= S.row();
	wsUint	n	= col();
	wsReal*	R	= sseVector(n);
/*

a <- Xt
b <- S
r <- 0
m <- c()
for (i in 1:n) { #x */
	for (wsUint i=0 ; i<n ; i++) {
		wsReal *Ra_Ak_ = sseVector(p);
		wsReal *Ra_ij = sseEmptyVec(p);
		for (wsUint k=0 ; k<p ; k++)
			Ra_Ak_[k] = A[k][i];
		for (wsUint k=0 ; k<p ; k++) {
/*
r.ij <- rep(0, p)
a.k. <- a[,i]
for (k in 1:p) {
a.ki <- a[k,i]
*/
			Ra_ij[k] += sseVV(Ra_Ak_, k+1, B[k]);
			sseVaV(Ra_ij, B[k], Ra_ij, k, Ra_Ak_[k]);
		}
		R[i] = sseVV(Ra_ij, p, Ra_Ak_);
		sseFree(Ra_Ak_);
		sseFree(Ra_ij);
	}
/*
 *# j
if (k>1) for (j in 1:(k-1)) {
r.ij[k] <- r.ij[k] + a.k.[j] * b[k,j]
r.ij[j] <- r.ij[j] + a.ki * b[k,j]
#     for (j in 1:p)
#      r.ij[j] <- r.ij[j] + a.ki * b[k,j]
}
r.ij[k] <- r.ij[k] + a.ki * b[k,k]
}
cat(r.ij, "\n")
m[i] <- sum(r.ij * a.k.)
r <- r + m[i]
#  for (j in 1:p)
#    for (k in 1:p)
#      t <- t + a[k,i]*b[k,j]
#    r <- t * a[j,i]
}

*/
	return cVector(R, n);
}

cVector cStdMatrix::diag(cDiagMatrix &D)
{
	cStdMatrix RR = *this * D;
	wsVec R = sseVector(row());
	for (wsUint i=0 ; i<row() ; i++)
		R[i] = sseVV(RR.get()[i], col(), get()[i]);
	return cVector(R, row());
}

cStdMatrix cStdMatrix::krnk(cStdMatrix &M)
{
	wsUint	N_r1	= row();
	wsUint	N_c1	= col();
	wsUint	N_r2	= M.row();
	wsUint	N_c2	= M.col();
	wsMat	Ra_ret	= krnkMM(get(), N_r1, N_c1, M.get(), N_r2, N_c2);
	return cStdMatrix(N_r1*N_r2, N_c1*N_c2, Ra_ret);
}

cStdMatrix cStdMatrix::mldivide(cStdMatrix& M, char B_trans/*=1*/)
{
	wsMat		Ra_mat	= get();
	wsUint		N_row	= row();
	wsUint		N_col	= col();
	wsUint		N_sz	= B_trans == 1 ? N_col : N_row;
	if (M.row() != N_sz) raise(MATERR_DIM_INCORRECT);
	wsSym		Ra_XX	= B_trans == 1 ? sseMtM(Ra_mat, N_row, N_col) :
			sseMpMt(Ra_mat, N_row, N_col);
	cSymMatrix	M_XX(Ra_XX, N_sz);
	cSymMatrix&	M_XXi	= M_XX.inv();
	cStdMatrix	M_Xy	= B_trans == 1 ? this->tM(M) : *this * M;
	cStdMatrix	M_ret	= M_XXi * M_Xy;
	delete &M_XXi;
	M_ret.setClean();

	return cStdMatrix(N_sz, M.col(), M_ret.get());
}

cVector cStdMatrix::mldivide(cVector& V, char B_trans/*=1*/)
{
	wsMat		Ra_mat	= get();
	wsUint		N_row	= row();
	wsUint		N_col	= col();
	wsUint		N_sz	= B_trans == 1 ? N_col : N_row;
	if (V.size() != N_sz) raise(MATERR_DIM_INCORRECT);
	wsSym		Ra_XX	= B_trans == 1 ? sseMtM(Ra_mat, N_row, N_col) :
		sseMpMt(Ra_mat, N_row, N_col);
	cSymMatrix	M_XX(Ra_XX, N_sz);
	cSymMatrix&	M_XXi	= M_XX.inv();
	cVector		V_Xy	= B_trans == 1 ? this->tV(V) : *this * V;
	cVector		V_ret	= M_XXi * V_Xy;
	delete &M_XXi;
	V_ret.setDontDealloc();

	return cVector(V_ret.get(), N_sz);
}

/*
 *
 * Predefined operations definition
 *
 */

wsRealCst cStdMatrix::sum()
{
	return sseMsum(get(), col(), row());
}

wsRealCst cStdMatrix::ssum()
{
	return sseMsqsum(get(), col(), row());
}

bool cStdMatrix::sym()
{
	/* Standard matrix is basically assumes NOT symmetric */
	return false;
}

cVector cStdMatrix::sumR(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsReal *Ra_ret = sseMsumR(get(), row(), col(), Rp_sumAll);

	if (!NA(R_mul)) {
		sseVpC(Ra_ret, R_mul, Ra_ret, row());
		if (Rp_sumAll) *Rp_sumAll *= R_mul;
	}
	return cVector(Ra_ret, row());
}

cVector cStdMatrix::sumC(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsReal *Ra_ret = sseMsumC(get(), row(), col(), Rp_sumAll);

	if (!NA(R_mul)) {
		sseVpC(Ra_ret, R_mul, Ra_ret, col());
		if (Rp_sumAll) *Rp_sumAll *= R_mul;
	}
	return cVector(Ra_ret, col());
}

cVector cStdMatrix::asumC(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsReal *Ra_ret = sseMasumC(get(), row(), col(), Rp_sumAll);

	if (!NA(R_mul)) {
		sseVpC(Ra_ret, R_mul, Ra_ret, col());
		if (Rp_sumAll) *Rp_sumAll *= R_mul;
	}
	return cVector(Ra_ret, col());
}

cStdMatrix cStdMatrix::acDiv(cVector &V_divider)
{
	wsUint	N_col	= col();
	if (N_col != V_divider.size()) raise(MATERR_DIM_INCORRECT);
	wsMat	Ra_s	= get();
	wsReal*	Ra_v	= V_divider.get();
	wsUint	N_row	= row();
	wsMat	Ra_r	= sseMatrix(N_row, N_col);

	for (wsUint i=0 ; i<N_row ; i++)
		sseAVdV(Ra_s[i], Ra_v, Ra_r[i], N_col);

	return cStdMatrix(N_row, N_col, Ra_r);
}

cStdMatrix cStdMatrix::rmerge(cStdMatrix& M_merged)
{
	/* Dim check */
	if (row() != M_merged.row()) raise(MATERR_DIM_INCORRECT);

	/* Create space */
	wsUint N_c1	= col(), N_c2 = M_merged.col();
	wsMat Ra_ret = sseMatrix(row(), N_c1 + N_c2);
	LOOP(i, row()) {
		memcpy(Ra_ret[i], get()[i], sizeof(wsReal)*N_c1);
		memcpy(Ra_ret[i]+N_c1, M_merged.get()[i], sizeof(wsReal)*N_c2);
	}

	return cStdMatrix(row(), N_c1+N_c2, Ra_ret);
}

cStdMatrix cStdMatrix::bmerge(cStdMatrix& M_merged)
{
	wsUint C = col();
	/* Dim check */
	if (C != M_merged.col()) raise(MATERR_DIM_INCORRECT);

	/* Create space */
	wsUint I = 0;
	wsUint N_r1	= row(), N_r2 = M_merged.row();
	wsMat Ra_ret = sseMatrix(N_r1 + N_r2, C);
	for (wsUint i=0 ; i<N_r1 ; i++,I++)
		memcpy(Ra_ret[I], get()[i], sizeof(wsReal)*C);
	for (wsUint i=0 ; i<N_r2 ; i++,I++)
		memcpy(Ra_ret[I], M_merged[i], sizeof(wsReal)*C);

	return cStdMatrix(N_r1+N_r2, col(), Ra_ret);
}

cStdMatrix cStdMatrix::rRepeat(wsUint N_rep)
{
	wsUint C = col();

	/* Create space */
	wsUint N_r1	= row(), N_r2 = N_r1 * N_rep;
	wsMat Ra_ret = sseMatrix(N_r2, C);
	for (wsUint i=0 ; i<N_r1 ; i++)
		for (wsUint j=0, l=i ; j<N_rep ; j++,l+=N_r1)
			memcpy(Ra_ret[l], get()[i], sizeof(wsReal)*C);

	return cStdMatrix(N_r2, col(), Ra_ret);
}

cStdMatrix cStdMatrix::ssRow(vInt Xv_idx)
{
	wsMat	R	= get();
	wsUint	N	= (wsUint)Xv_idx.size();
	/* Create space */
	wsMat Ra_ret = sseMatrix(N, col());
	wsUint I = 0;
	FOREACHDO (vInt_it, Xv_idx, i, I++) {
		if ((int)row() <= *i) halt("Index out of bounds from [%d]th element [%d >= %d]",
			I, *i, row());
		memcpy(Ra_ret[I], R[*i], sizeof(wsReal)*col());
	}

	return cStdMatrix(N, col(), Ra_ret);
}

cVector cStdMatrix::meanR()
{
	wsVec Ra_ret = sseMmeanR(get(), row(), col());
	return cVector(Ra_ret, row());
}

cVector cStdMatrix::meanC()
{
	wsVec Ra_ret = sseMmeanC(get(), row(), col());
	return cVector(Ra_ret, col());
}

cSymMatrix cStdMatrix::covR()
{
	return cSymMatrix(sseMcovRow(get(), row(), col()), row());
}

cStdMatrix& cStdMatrix::inv()
{
	/* Check */
	if (!sqr()) raise(MATERR_DIM_INCORRECT);
	wsMat	Ra_ret = SO_invMatrix(get(), row());
	return *(new cStdMatrix(row(), row(), Ra_ret));
}

void cStdMatrix::rem()
{
	if (Ra_data) switch (X_actClean) {
		case MATDEL_ALL:
			deallocMatrix(get(), row(), (void *)1);
			break;
		case MATDEL_ROW: {
			wsMat Ra_v = get();
			DEALLOC(Ra_v);
			} break;
		case MATDEL_NONE:
			/* Do nothing when don't delete any */
			break;
	}
	Ra_data	= NULL;
	N_row	= 0;
	N_col	= 0;
}

void cStdMatrix::file(wsStrCst S_ext, wsUint N_prec/*=0*/)
{
	exportMatrix(S_ext, get(), row(), col(), NULL, NULL, N_prec);
}

void cStdMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/)
{
	exportMatrix(S_ext, get(), row(), col(), Sa_colNames, NULL, N_prec);
}

void cStdMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/)
{
	exportMatrix(S_ext, get(), row(), col(), Sa_colNames, Sa_rowNames, N_prec);
}

/* CONFIRMED 130316 */
cStdMatrix& cStdMatrix::operator-(wsReal R)
{
	cStdMatrix *T = new cStdMatrix(row(), col());
	sseMaC(get(), row(), col(), -R, T->get());
	return *T;
}

/* CONFIRMED 130316 */
cStdMatrix& cStdMatrix::operator+(wsReal R)
{
	cStdMatrix *T = new cStdMatrix(row(), col());
	sseMaC(get(), row(), col(), R, T->get());
	return *T;
}

/* CONFIRMED 130318 */
cStdMatrix& cStdMatrix::operator*(wsReal R)
{
	cStdMatrix *T = new cStdMatrix(row(), col());
	sseMpC(cget(), R, T->get(), row(), col());
	return *T;
}

cVector cStdMatrix::operator*(cVector &V)
{
	/* Check */
	if (col()!=V.size()) raise(MATERR_DIM_INCORRECT);
	wsReal *Ra_pV = sseMpV(get(), row(), col(), V.get());
	return cVector(Ra_pV, row());
}

cStdMatrix& cStdMatrix::addDiag(wsRealCst R)
{
	wsUint	N = row() > col() ? col() : row();
	wsMat	P = get();
	LOOP (i, N) P[i][i] += R;
	return *this;
}

void cStdMatrix::setR(wsStrCst S_varName)
{
#ifdef USE_R
	R().sendMatrix(S_varName, get(), row(), col());
#else
	halt("USE_R is required");
#endif
}

void cStdMatrix::getR(wsStrCst S_varName)
{
#ifdef USE_R
	wsUint N_r, N_c;
	R().dim(S_varName, &N_r, &N_c);
	init(N_r, N_c);
	R().recvMatrix(S_varName, get(), N_r, N_c);
#else
	halt("USE_R is required");
#endif
}

wsReal cStdMatrix::tr()
{
	if (!sqr()) halt("Cannot get trace of non-square matrix");
	wsUint N_r		= row();
	wsReal R_ret	= W0;
	wsReal **Ra_m	= get();

	LOOP (i, N_r) R_ret += Ra_m[i][i];

	return R_ret;
}

cVector cStdMatrix::diag()
{
	if (!sqr()) halt("Cannot get trace of non-square matrix");
	wsUintCst	N_r		= row();
	wsVec		Ra_ret	= sseVector(N_r);
	wsMat		Ra_m	= get();

	LOOP (i, N_r) Ra_ret[i] = Ra_m[i][i];

	return cVector(Ra_ret, N_r);
}

cVector cStdMatrix::diag2()
{
	wsUintCst	N_r		= row();
	wsUintCst	N_c		= col();
	wsVec		Ra_ret	= sseVector(N_r);
	wsMat		Ra_m	= get();

	LOOP(i, N_r) Ra_ret[i] = sseVV(Ra_m[i], N_c);

	return cVector(Ra_ret, N_r);
}

cStdMatrix& cStdMatrix::clone()
{
	wsReal **Ra_ret = sseMatrixP(row(), col(), get());
	return *(new cStdMatrix(row(), col(), Ra_ret));
}

cStdMatrix& cStdMatrix::sqrt()
{
	halt("SYSERR : Unsupported feature [Square root of normal matrix] called");
	return *(new cStdMatrix());
}

wsRealCst cStdMatrix::mean()
{
	return sum() / (wsReal)(row()*col());
}

cVector cStdMatrix::eigen(cStdMatrix &M_eVecT)
{
	char	B_fail	= 1;
	wsUint	N_sz	= row();
	if (N_sz != col())
		halt("SYSERR: Non-square matrix eigen() called");
	wsMat	Ra_full = sseMatrixP(N_sz, N_sz, get());
	wsMat	Ra_E	= NULL;
	wsReal*	Ra_eI	= EIGENDECOMPOSITION(Ra_full, N_sz, &Ra_E);
	if (B_fail) halt("Eigendecomposition failed");
	M_eVecT.init(N_sz, N_sz, Ra_E);
	sseUnmat(Ra_full, N_sz);

	return cVector(Ra_eI, N_sz);
}

cVector cStdMatrix::eigen()
{
	char	B_fail	= 1;
	wsUint	N_sz	= row();
	if (N_sz != col())
		halt("SYSERR: Non-square matrix eigen() called");
	wsMat	Ra_full = sseMatrixP(N_sz, N_sz, get());
	wsReal*	Ra_eI	= EIGENDECOMPOSITION(Ra_full, N_sz);
	if (B_fail) halt("Eigendecomposition failed");
	sseUnmat(Ra_full, N_sz);

	return cVector(Ra_eI, N_sz);
}

cStdMatrix cStdMatrix::subsetPtr(wsUint N_s, wsUint N_e)
{
	wsUint	N_sz	= N_e - N_s + 1;
	wsMat	Ra_r	= NULL;
	wsMat	Ra_p	= get();
	wsUint	N_row	= row();

	/* Range check */
	if (N_s>=N_row || N_e>=N_row || N_s > N_e)
		halt("Invalid range requested");

	/* Allocate space */
	wsAlloc(Ra_r, wsVec, N_sz);

	/* Set the pointer */
	for (wsUint i=N_s,j=0 ; i<=N_e ; i++,j++)
		Ra_r[j] = Ra_p[i];

	return cStdMatrix(N_sz, col(), Ra_r, MATDEL_ROW);
}

cStdMatrix cStdMatrix::subsetPtr(vInt& Nv_idx)
{
	wsUint	N_sz	= (wsUint)Nv_idx.size();
	wsMat	Ra_r	= NULL;
	wsMat	Ra_p	= get();
	wsUint	N_row	= row();

	/* Allocate space */
	wsAlloc(Ra_r, wsVec, N_sz);

	/* Set the pointer */
	wsUint j = 0;
	FOREACHDO (vInt_it, Nv_idx, i, j++) {
		/* Range check */
		if ((wsUint)*i>=N_row) halt("Invalid range requested");
		Ra_r[j] = Ra_p[*i];
	}

	return cStdMatrix(N_sz, col(), Ra_r, MATDEL_ROW);
}

cStdMatrix cStdMatrix::subsetCol(vInt& Nv_idx)
{
	wsMat Ra_t = get();
	wsMat Ra_r = sseMatrix(row(), (wsUint)Nv_idx.size());
	LOOP (i, row()) {
		wsUint J = 0;
		FOREACHDO (vInt_it, Nv_idx, j, J++)
			Ra_r[i][J] = Ra_t[i][*j];
	}
	return cStdMatrix(row(), (wsUint)Nv_idx.size(), Ra_r);
}

wsVec cStdMatrix::operator[](wsUint N_rIdx)
{
	if (N_rIdx >= row())
		raise(MATERR_DIM_INCORRECT);
	return Ra_data[N_rIdx];
}

cStdMatrix& cStdMatrix::subset(char *Ba_YorN, wsUint N_Yval)
{
	if (row() != col()) raise(MATERR_DIM_INCORRECT);

	wsUint N_nLen = 0;

	if (N_Yval == 1) {
		/* 1 then include, 0 exclude */

		/* Get length */
		for (wsUint i=0,N=row() ; i<N ; i++)
			if (Ba_YorN[i]) N_nLen++;
	} else {
		/* 1 then exclude, 0 include */

		/* Get length */
		for (wsUint i=0,N=row() ; i<N ; i++)
			if (!Ba_YorN[i]) N_nLen++;
	}

	/* Make matrix */
	wsMat	Ra_ret	= sseMatrix(N_nLen, N_nLen);
	wsMat	Ra_src	= get();

	/* Make subset */
	if (N_Yval == 1) {
		/* 1 then include, 0 exclude */

		/* Get length */
		for (wsUint i=0,I=0,N=row() ; i<N ; i++) {
			if (!Ba_YorN[i]) continue;

			for (wsUint j=0,J=0 ; j<N ; j++) {
				if (!Ba_YorN[j]) continue;

				Ra_ret[I][J] = Ra_src[i][j];
				J++;
			}
			I++;
		}
	} else {
		/* 1 then exclude, 0 include */

		/* Get length */
		for (wsUint i=0,I=0,N=row() ; i<N ; i++) {
			if (Ba_YorN[i]) continue;

			for (wsUint j=0,J=0 ; j<N ; j++) {
				if (Ba_YorN[j]) continue;

				Ra_ret[I][J] = Ra_src[i][j];
				J++;
			}
			I++;
		}
	}

	return *(new cStdMatrix(N_nLen, N_nLen, Ra_ret));
}

/* Vectorize */
cVector cStdMatrix::vec()
{
	wsUint	r = row();
	wsUint	c = col();
	wsUint	N = r * c;
	wsMat	S = get();
	wsReal*	R = sseVector(N);
	for (wsUint i=0 ; i<r ; i++) {
		wsUint k = i;
		for (wsUint j=0 ; j<c ; j++,k+=r)
			R[k] = S[i][j];
	}

	return cVector(R, N);
}

wsRealCst cStdMatrix::getMin(wsUint* Np_idxR/*=NULL*/, wsUint* Np_idxC/*=NULL*/)
{
	return sseMmin(get(), row(), col(), Np_idxR, Np_idxC);
}

wsRealCst cStdMatrix::getMax(wsUint* Np_idxR/*=NULL*/, wsUint* Np_idxC/*=NULL*/)
{
	return sseMmax(get(), row(), col(), Np_idxR, Np_idxC);
}

cStdMatrix cStdMatrix::subR(cVector& V_val)
{
	if (row() != V_val.size())
		raise(MATERR_DIM_INCORRECT);
	wsMat	Ra_full = sseMatrix(row(), col());
	wsMat	Ra_R = get();
	wsVec	Ra_v = V_val.get();
	LOOP (i, row())
		sseVsC(Ra_R[i], Ra_v[i], Ra_full[i], col());
	return cStdMatrix(row(), col(), Ra_full);
}

} // End namespace ONETOOL
