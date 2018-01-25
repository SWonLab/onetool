#include <limits>
#include "analyses/emai.h"
#include "input/stream.h"
#include "utils/util.h"
#include "utils/matrix.h"

namespace ONETOOL {

const char	*Sa_EmaiEst_headers[EMAI_EST_SIZE]	=
	{	"sig2",			"sig2g",	"logL",		"Var(sig2)",
		"Var(sig2g)",	"h^2",		"Var(h^2)",	"Var(sig2+sig2g)" };

typedef struct _xSpecDecomp
{
	wsReal **Ra_phi;
	wsReal **Ra_X, **Ra_XXinv;
	wsUint N_sample, N_inpCol;
} xSpecDecomp;

/* Spectral decomposition method */
//	loglikelihood2 <- function(ratio,dat,K,method="ML"){
// double logL_spectralDecomp(double R_ratio, void *Vp_data)
// {
// 	xSpecDecomp *Xp = (xSpecDecomp *)Vp_data;
// 	//	delta <- ratio*ratio
// 	wsReal R_det = SQR(R_ratio);
// 
// 	/* Do eigendecomposition of phi matrix, but not need */
// 	//	k <- eigen(K)
// 	//	xi <- k$values
// 
// 	//	N <- nrow(dat)					# # of inds
// 	//	X <- matrix(rep(1,N),ncol=1)	# Intercept matrix
// 	//	S <- diag(N)-X%*%solve(t(X)%*%X)%*%t(X)	# matrix I-H
// 	wsReal **Ra_S = sseMMMt(Xp->Ra_X, Xp->N_sample, Xp->N_inpCol,
// 		Xp->Ra_XXinv, Xp->N_inpCol, Xp->N_inpCol,
// 		Xp->Ra_X, Xp->N_sample, Xp->N_inpCol);
// 	//	s <- eigen(S);
// 
// 	//	nminusp <- N-dim(X)[2]
// 	wsUint N_NminusP = Xp->N_sample - Xp->N_inpCol;
// 	//	A <- s$vectors[,1:nminusp]
// 	wsReal **Ra_At = NULL;
// 
// 	//	AZKZA <- t(A)%*%K%*%A # -> t(A) %*% K %*% t(A)
// 	wsReal **Ra_AtKA = sseMMMt(Ra_At, N_NminusP, Xp->N_sample,
// 		Xp->Ra_phi, Xp->N_sample, Xp->N_sample,
// 		Ra_At, Xp->N_sample, N_NminusP);
// 	//	azkza <- eigen(AZKZA)
// 	//	lambda <- azkza$values
// 	//	UR <- A%*%azkza$vectors
// 
// 	//	eta <- t(UR)%*%dat$y
// 	//	logR <- log(sum(eta^2/(lambda+delta)))
// 	wsReal R_logR = W0;
// 	wsReal R_log2pi = log(W2*WISARD_PI);
// 	if (OPT_ENABLED(ml)) {
// 		wsReal R_sumLogXiDet = W0;
// 
// 		return REAL_CONST(0.5)*(Xp->N_sample*
// 			(log(Xp->N_sample) - R_log2pi - W1)
// 			- Xp->N_sample*R_logR - R_sumLogXiDet);
// 		//		l <- 0.5*( N*(log(N)-log(2*pi)-1) - N*logR - sum(log(xi+delta)) ) else
// 	} else {
// 		wsReal R_sumLogLiDet = W0;
// 
// 		return REAL_CONST(0.5)*(N_NminusP*
// 			(log(N_NminusP) - R_log2pi - W1)
// 			- N_NminusP*R_logR - R_sumLogLiDet);
// 		//		l <- 0.5*( nminusp*(log(nminusp)-log(2*pi)-1) - nminusp*logR - sum(log(lambda+delta)))
// 	}
// }
// 

wsReal** subsetMatrix(wsReal **Rp_inpMatrix,
	wsUint N_rowStart, wsUint N_rowCount,
	wsUint N_colStart, wsUint N_colCount)
{
	wsReal **Ra_ret = NULL;
	wsAlloc(Ra_ret, wsReal*, N_rowCount);
	for (wsUint i=0 ; i<N_rowCount ; i++) {
		wsAlloc(Ra_ret[i], wsReal, N_colCount);

		memcpy(Ra_ret[i], Rp_inpMatrix[i+N_rowStart]+N_colStart, sizeof(wsReal)*N_colCount);
	}

	return Ra_ret;
}

wsReal** GJ_invMatrix(wsReal **Ra_matrix, wsUint N_row)
{
	wsUint i, j, k;

	wsReal **Ra_inv = NULL, **Ra_matCopy = NULL;
	wsAlloc(Ra_inv, wsReal*, N_row);
	wsAlloc(Ra_matCopy, wsReal*, N_row);
	for (i=0 ; i<N_row ; i++) {
		wsCalloc(Ra_inv[i], wsReal, N_row);
		wsAlloc(Ra_matCopy[i], wsReal, N_row);
		memcpy(Ra_matCopy[i], Ra_matrix[i], sizeof(wsReal)*N_row);
	}
	for (i=0 ; i<N_row ; i++)
		Ra_inv[i][i] = 1;
	for (k=0 ; k<N_row ; k++) {
		for (i=k; i<N_row ; i++) {
			wsReal val = 1.0f / Ra_matCopy[i][k];
			for (j=k ; j<N_row ; j++)
				Ra_matCopy[i][j] *= val;
			for (j=0 ; j<N_row ; j++)
				Ra_inv[i][j] *= val;
		}

		for (i=k+1 ; i<N_row ; i++) {
			for (j=k ; j<N_row ; j++)
				Ra_matCopy[i][j] -= Ra_matCopy[k][j];
			for (j=0 ; j<N_row ; j++)
				Ra_inv[i][j] -= Ra_inv[k][j];
		}
	}

	for (int i=N_row-2 ; i>=0; i--) {
 		for (int j=N_row-1; j>i ; j--) {
			for (k=0 ; k<N_row ; k++) //this part can be can canceled (the two lines)
				Ra_inv[i][k] -= Ra_matCopy[i][j] * Ra_inv[j][k];
			for (k=0 ; k<N_row ; k++) //this part can be can canceled (the two lines)
				Ra_matCopy[i][k] -= Ra_matCopy[i][j] * Ra_matCopy[j][k];
		}
		if (i == 0)
			break;
	}

	return Ra_inv;
}

// wsReal** LU_invMatrix(wsReal **Ra_matrix, wsUint N_row)
// {
// 	wsUint i, j, k;
// 
// 	wsReal **Ra_inv = NULL;
// 	MULTI_MALLOC(Ra_inv, wsReal*, N_row);
// 	for (i=0 ; i<N_row ; i++) {
// 		MULTI_MALLOC(Ra_inv[i], wsReal, N_row);
// 		memcpy(Ra_inv[i], Ra_matrix[i], sizeof(wsReal)*N_row);
// 	}
// 
// 	for (k=0 ; k<(N_row-1) ; k++) {
// 		for (j=k+1 ; j<N_row ; j++) {
// 			wsReal R_x = Ra_inv[j][k]/Ra_inv[k][k];
// 
// 			for(i=k ; i<N_row ; i++)
// 				Ra_inv[j][i] = Ra_inv[j][i]-R_x*Ra_inv[k][i];
// 
// 			Ra_inv[j][k] = R_x;
// 		}
// 	}
// 
// 	return Ra_inv;
// }
// 

wsReal** LUdecomp(wsReal **A, int n)
{
	int i, j, k, p;

	if (A[0][0] == 0.0) return NULL;
	for (i=1 ; i<n ; i++)
		A[i][0] /= A[0][0];

	for (k=1 ; k<n ; k++) {
		for (j=k ; j<n ; j++) {
			wsReal s = 0.0;
			for (p=0 ; p<k ; p++)
				s += A[k][p] * A[p][j];
			A[k][j] -= s;
		}
		if (A[k][k] == 0.0) return NULL;

		for (i=k+1 ; i<n ; i++) {
			wsReal s = 0.0;
			for (p=0 ; p<k ; p++)
				s += A[i][p] * A[p][k];
			A[i][k] -= s;
			A[i][k] /= A[k][k];
		}
	}
	return A;
}

void Unit_Lower_Triangular_Solve(wsReal **L, wsReal B[], wsReal x[], wsUint n)
{
	wsUint i, k;

	//         Solve the linear equation Lx = B for x, where L is a unit lower
	//         triangular matrix.                                      

	x[0] = B[0];
	for (k = 1 ; k < n; k++) {
#ifdef USE_SSE
		sse_t s = sseSet(0.0);
		wsReal ss;
		wsUint N_med = getMed(k);
		for (i=0 ; i<N_med ; i+=sseJmp)
			s = sseAdd(s, sseMul(*((sse_t *)(x+i)), *((sse_t *)(L[k]+i))));
		sseSum(s, ss);
		for (i=N_med ; i < k; i++)
			ss += x[i] * L[k][i];
		x[k] = B[k] - ss;
#else
		for (i=0, x[k]=B[k] ; i<k ; i++)
			x[k] -= x[i] * L[k][i];
#endif
	}
}

void Unit_Lower_Triangular_Solve2(wsReal **L, wsUint p, wsReal x[], wsUint n)
{
	wsUint i, k;

	//         Solve the linear equation Lx = B for x, where L is a unit lower
	//         triangular matrix.                                      

	x[0] = p==0?W1:W0;
	for (k = 1 ; k < n; k++) {
		wsReal xx = k==p?W1:W0;
#ifdef USE_SSE
		wsReal	ss;
		sse_t	s		= sseSet(0.0);
		wsUint	N_med	= getMed(k);
		for (i=0 ; i<N_med ; i+=sseJmp)
			s = sseAdd(s, sseMul(*((sse_t *)(x+i)), *((sse_t *)(L[k]+i))));
		sseSum(s, ss);
		for (i=N_med ; i<k ; i++)
			ss += x[i] * L[k][i];
		x[k] = xx - ss;
#else
		for (i=0 ; i < k; i++)
			x[k] -= x[i] * L[k][i];
#endif
	}
}

int Upper_Triangular_Solve(wsReal **U, wsUint p, wsReal B[], wsReal x[], wsUint n,
	wsReal *Ra_ret)
{
	wsUint	i;
	int		k;

	//         Solve the linear equation Ux = B for x, where U is an upper
	//         triangular matrix.                                      
	for (k=n-1 ; k>=0 ; k--) {
		if (U[k][k] == W0) return -1;
//		if (*(U + k) == 0.0) return -1;           // The matrix U is singular
#ifdef USE_SSE
		wsUint N_s = getMed(k+sseJmp);
		if (N_s > n) N_s = n;
		wsUint N_e = getMed(n);
		wsReal xs = W0;
		if (N_s > n)
			N_s = n;
		for (i=k+1 ; i<N_s ; i++)
			xs += x[i] * U[k][i];
		sse_t s = sseSet(0.0);
		for (i=N_s ; i<N_e ; i+=sseJmp)
			s = sseAdd(s, sseMul(*((sse_t*)(x+i)), *((sse_t*)(U[k]+i))));
		sseAppendSum(s, xs);
		for (i=N_e ; i<n ; i++)
			xs += x[i] * U[k][i];
		x[k] = (B[k]-xs)/U[k][k];
//		if (k<=p)
			Ra_ret[k] = x[k];
#else
		x[k] = B[k];
		for (i = k + 1; i < n; i++) x[k] -= x[i] * U[k][i];//*(U + i);
		x[k] /= U[k][k];
//		if (k<=p)
			Ra_ret[k] = x[k];//*(U + k);
#endif
	}

	return 0;
}

void LTsolve(wsReal **L, wsReal x[], wsUint n, wsUint p)
{
	wsUint i, k;

	//         Solve the linear equation Lx = B for x, where L is a unit lower
	//         triangular matrix.                                      

	x[0] = !p ? W1 : W0;
	for (k=1 ; k<n ; k++) {
		x[k] = k==p?W1 : W0;

		wsReal s = W0;
		wsUint N_med = 0;
#ifdef USE_SSE
		N_med = getMed(k);
		sse_t ss = sseSet(0.0);
		for (i=0 ; i<N_med ; i+=sseJmp)
			ss = sseAdd(ss, sseMul(*((sse_t *)(x+i)), *((sse_t *)(L[k]+i))));
		sseSum(ss, s);
#endif
		for (i=N_med ; i<k ; i++)
			s += x[i] * L[k][i];

		x[k] -= s;
	}
}
// int LTsolve(wsReal **L, double x[], int n, int p)
// {
// 	int i, k;
// 
// 	//         Solve the linear equation Lx = B for x, where L is a lower
// 	//         triangular matrix.                                      
// 	// for some p=k, B is zero vector but kth element is 1
// 
// 	for (k=0 ; k<n ; k++) {
// 		if (L[k][k] == 0.0) return -1;           // The matrix L is singular
// 		x[k] = k==p ? 1.0 : 0.0;
// 		wsReal s = 0.0;
// 		for (i = 0; i < k; i++) s += x[i] * L[k][i];
// 		x[k] -= s;
// 		x[k] /= L[k][k];
// 	}
// 
// 	return 0;
// }

int UTsolve(wsReal **U, wsReal x[], wsUint n, wsReal r[])
{
	wsUint i;
	int k;

	//         Solve the linear equation Ux = B for x, where U is an upper
	//         triangular matrix.                                      

	// k=n-1
	if (U[n-1][n-1] == W0) return -1;           // The matrix U is singular
	x[n-1] /= U[n-1][n-1];
	r[n-1] = x[n-1];

	for (k=n-2 ; k>=0 ; k--) {
		if (U[k][k] == W0) return -1;           // The matrix U is singular

		wsReal s = W0;

#ifdef USE_SSE
		wsUint N_s = getMed(k+2);
		wsUint N_e = getMed(n);
		if (N_s > n) N_s = n;
		for (i=k+1 ; i<N_s ; i++) s += x[i] * U[k][i];
		sse_t ss = sseSet(0.0);
		for (i=N_s ; i<N_e ; i+=sseJmp)
			ss = sseAdd(ss, sseMul(*((sse_t *)(x+i)), *((sse_t *)(U[k]+i))));
		sseAppendSum(ss, s);
		for (i=N_e ; i<n ; i++) s += x[i] * U[k][i];
#else
		for (i=k+1 ; i<n ; i++) s += x[i] * U[k][i];
#endif
		x[k] -= s;
		x[k] /= U[k][k];
		r[k] = x[k];
	}

	return 0;
}

int symLUsolve(wsReal **LU, wsUint p, wsReal x[], int n)
{
	wsReal *Ra_interm = NULL;
	sseMalloc(Ra_interm, wsReal, n);
	//         Solve the linear equation Lx = B for x, where L is a lower
	//         triangular matrix with an implied 1 along the diagonal.

//	Unit_Lower_Triangular_Solve(LU, B, x, n);
//	Unit_Lower_Triangular_Solve2(LU, p, Ra_interm, n);
	LTsolve(LU, Ra_interm, n, p);

	//         Solve the linear equation Ux = y, where y is the solution
	//         obtained above of Lx = B and U is an upper triangular matrix.

//	int q = Upper_Triangular_Solve(LU, p, Ra_interm, Ra_interm, n, x);
	int q = UTsolve(LU, Ra_interm, n, x);
	sseFree(Ra_interm);
	return q;
}

#define USE_SSE

void LU(wsReal **Ra_mat, wsUint N_row)
{
	wsUint	i,j,k;
	wsReal	x;

	for(k=0;k<N_row-1;k++) {
		for(j=k+1;j<N_row;j++) {
			x = Ra_mat[j][k]/Ra_mat[k][k];
#if 0
			UINT_t N_s = getMed(k+sseJmp-1);
			UINT_t N_e = getMed(N_row+1);
			sse_t sse_x = sseSet(x);
			for(i=k;i<N_s;i++)
				Ra_mat[j][i] -= x*Ra_mat[k][i];

			for(i=N_s;i<N_e;i+=sseJmp) {
				sse_t *sse_m1 = (sse_t *)(Ra_mat[j] + i);
				sse_t *sse_m2 = (sse_t *)(Ra_mat[k] + i);

				*sse_m1 = sseSub(*sse_m1, sseMul(sse_x, *sse_m2));
//				Ra_mat[j][i] -= x*Ra_mat[k][i];
			}

			for(i=N_e;i<=N_row;i++)
				Ra_mat[j][i] -= x*Ra_mat[k][i];
#else
			for(i=k;i<N_row;i++)
				Ra_mat[j][i] -= x*Ra_mat[k][i];
#endif
			Ra_mat[j][k] = x;
		}
	}
}

wsReal diagSum3matrix(wsReal **Ra_m1, wsUint N_row1, wsUint N_col1,
	wsReal **Ra_m2, wsUint N_row2, wsUint N_col2,
	wsReal **Ra_m3, wsUint N_row3, wsUint N_col3) {

	wsUint N_i = N_row1 > N_col3 ? N_col3 : N_row1;

	/* p = N_mat2Col
	 * m = N_mat1Col
	 * (ABC)_ii=¥Ò_(k=1)^p (¥Ò_(j=1)^m A_ij B_jk ) C_ki
	 */

	wsReal R_ret = W0;

	for (wsUint i=0 ; i<N_i ; i++) {
		for (wsUint k=0 ; k<N_col2 ; k++) {
			wsReal R_subSum = W0;

			for (wsUint j=0 ; j<N_col1 ; j++)
				R_subSum += Ra_m1[i][j] * Ra_m2[j][k];

			R_ret += R_subSum * Ra_m3[k][i];
		}
	}

	return R_ret;
}

wsReal diagSum3matrixSSE(wsReal **Ra_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Ra_mat2, wsUint N_mat2Row, wsUint N_mat2Col,
	wsReal **Ra_mat3, wsUint N_mat3Row, wsUint N_mat3Col, char B_mat3tp/*=0*/) {

	wsUint N_i = N_mat1Row > N_mat3Col ? N_mat3Col : N_mat1Row;

	/* p = N_mat2Col
	 * m = N_mat1Col
	 * (ABC)_ii=¥Ò_(k=1)^p (¥Ò_(j=1)^m A_ij B_jk ) C_ki
	 */

	wsReal R_ret = W0;

#ifdef USE_SSE
	wsReal *Ra_interm = NULL;
	sseMalloc(Ra_interm, wsReal, N_mat2Col);
	wsUint N_med = getMed(N_mat2Col);

	for (wsUint i=0 ; i<N_i ; i++) {
		wsReal *Rp_m1 = Ra_mat1[i];
		memset(Ra_interm, 0x00, sizeof(wsReal)*N_mat2Col);

		for (wsUint j=0 ; j<N_mat1Col ; j++) {
			sse_t sse_m1_il = sseSet(Rp_m1[j]); // m1[i,j]

			/* SSE-enabled part */
			for (wsUint k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_m2_lj	= (sse_t *)(&(Ra_mat2[j][k])); // m2[j,k:(k+3)]
				sse_t *sse_mR		= (sse_t *)(&(Ra_interm[k])); // m3[i,k:(k+3)] (j)
				*sse_mR = sseAdd(*sse_mR, sseMul(sse_m1_il, *sse_m2_lj));
			}
			/* Rest part */
			for (wsUint k=N_med ; k<N_mat2Col ; k++)
				Ra_interm[k] += Rp_m1[j] * Ra_mat2[j][k];
		}
		if (B_mat3tp == 0) {
			for (wsUint j=0 ; j<N_mat2Col ; j++)
				R_ret += Ra_interm[j]*Ra_mat3[j][i];
		} else {
#ifdef USE_SSE
			sse_t sse_sum = sseSet(0.0);
			for (wsUint j=0 ; j<N_mat2Col ; j+=sseJmp) {
				sse_t *sse_int	= (sse_t *)(Ra_interm + j);
				sse_t *sse_m3	= (sse_t *)(Ra_mat3[i] + j);
				sse_sum = sseAdd(sse_sum, sseMul(*sse_int, *sse_m3));
			}
			sseAppendSum(sse_sum, R_ret);
			for (wsUint j=N_med ; j<N_mat2Col ; j++)
				R_ret += Ra_interm[j]*Ra_mat3[i][j];
#else
			for (wsUint j=0 ; j<N_mat2Col ; j++)
				R_ret += Ra_interm[j]*Ra_mat3[i][j];
#endif
		}
	}
	sseFree(Ra_interm);
#else
	for (wsUint i=0 ; i<N_i ; i++) {
		for (wsUint k=0 ; k<N_mat2Col ; k++) {
			wsReal R_subSum = 0.0f;

			for (wsUint j=0 ; j<N_mat1Col ; j++)
				R_subSum += Ra_mat1[i][j] * Ra_mat2[j][k];

			if (B_mat3tp == 0)
				R_ret += R_subSum * Ra_mat3[k][i];
			else
				R_ret += R_subSum * Ra_mat3[i][k];
		}
	}
#endif

	return R_ret;
}

wsReal** (*getInvMethod())(wsReal**, wsUint)
{
	return SO_invMatrix;
}

void cEmAiAnalysisV2::_init(cIO *Cp_IO, cStdMatrix &M_inpXt,
	cMatrix &M_phi, char B_isREML, char *Ba_iExcl,
	cStdMatrix *Mp_eVec, cDiagMatrix *Mp_eVal)
{
	B_quiet			= 0;
	Mp_Xt			= &M_inpXt;
	Mp_phi			= &M_phi;
	Ba_excl			= Ba_iExcl;
	B_doSpecDcmp	= OPT_ENABLED(specdcmp);
	B_doEM			= 1;
	Ra_sampwgt		= NULL;

	/* Load sample weights if assigned */
	wsUint		N_samp = getIO()->sizeSample();
	if (IS_ASSIGNED(sampleweight)) {
		cStrFile	C_sw(OPT_STRING(sampleweight), "Sample weighting file");
		mSamp&		Xm_samp = getIO()->getSampleData();
		vSampPtr&	Xv_samp = getIO()->getSample();
		char*		Sp_buf = NULL;
		wsAlloc(Sp_buf, char, 4096);

		Ra_sampwgt = sseVector(N_samp);
		sseVinit(Ra_sampwgt, N_samp, WISARD_NAN);

		for (wsUint L=1 ; C_sw.gets(Sp_buf, 4096) ; L++) {
			char *a, *b, *c, *d;
			getString(&Sp_buf, &a);
			getString(&a, &b);
			getString(&b, &c);
			mSamp_it X_find = Xm_samp.find(a);
			if (X_find == Xm_samp.end()) continue;
			xSample&	X_samp = X_find->second;

			/* Get weight */
			wsReal R_wgt = (wsReal)str2dbl(b, &d);
			if (d && d[0])
				halt_fmt(WISARD_INVL_FILE_INVALID, "Sample weight file",
					OPT_STRING(sampleweight), c, L);

			/* If the sample is not available, skip */
			if (X_samp.N_idx == SAMP_NODATA) continue;
			Ra_sampwgt[X_samp.N_idx] = R_wgt;
		}
		DEALLOC(Sp_buf);
		/* Check all samples have their weight */
		for (wsUint i=0 ; i<N_samp ; i++)
			if (NA(Ra_sampwgt[i])) halt("Sample [%s::%s] does not have "
				"sample weight!", Xv_samp[i]->S_FID.c_str(),
				Xv_samp[i]->S_IID.c_str());

		LOGnote("Sample weights successfully loaded\n");
	}

	/* Do eigendecomposition if required */
	if (Mp_eVec == NULL) {
		wsReal	**Ra_eVec	= NULL; /* P^t */
		char	B_failed	= 0;
		/* Eigendecomposition with NO negative correction */
#ifdef USE_ED2
		wsReal	*Ra_eVal	= EIGENDECOMPOSITION(M_phi.get(), M_phi.row(),
			&Ra_eVec);
#else
		wsReal	*Ra_eVal	= eigenDecomp(M_phi.get(), M_phi.row(),
			&B_failed, &Ra_eVec, 0, 0);
#endif
		if (B_failed != 0) halt("Eigendecomposition failed");
		Mp_EV		= new cDiagMatrix(N_anaSamp, Ra_eVal);
		Mp_EVX		= new cStdMatrix(N_anaSamp, N_anaSamp, Ra_eVec);
		B_dealloc	= 1;
	} else {
		Mp_EVX		= Mp_eVec;
		Mp_EV		= Mp_eVal;
		B_dealloc	= 0;
	}

	/* Get the number of nonzero eigenvalues */
	if (Mp_EV->nn0() != Mp_EV->row())
		/* Turn off EM algorithm when the matrix is not full rank */
		B_doEM = 0;

	// sig2, sig2g, logL, Var(sig2), Var(sig2g), h^2, Var(h^2), Var(s2+s2g)
	// ^b_i Var(^b_i) <- Starts from index EMAI_EST_SIZE, paired
	wsUint N_est = EMAI_EST_SIZE+Mp_Xt->row()*2;

	/*
	 *
	 * Preparation
	 *
	 */
	wsUint N_pheno	= M_Yt.row();
	Ra_estimates	= sseMatrix(N_pheno, N_est);
	Ra_calcBLUP		= sseMatrix(N_pheno, Mp_Xt->col());
	Ra_calcPred		= sseMatrix(N_pheno, Mp_Xt->col());
	Ra_alpha		= sseMatrix(N_pheno, Mp_Xt->row());

	V_YtY	= M_Yt.rSS();			/* yy <- t(y)%*%y for each y */
	M_XtY_t	= M_Yt.Mt(*Mp_Xt);		/* Xy <- t(X)%*%y */
	M_XtX	= Mp_Xt->Mt();			/* XX <- t(X)%*%X */
	/* 130714 */
#ifdef USE_ED2
	M_XtP	= Mp_Xt->Mt(*Mp_EVX);	/* XtP <- t(X)%*%P */
	M_YtP	= M_Yt.Mt(*Mp_EVX);		/* YtP <- t(y)%*%P */
#else
	M_XtP	= *Mp_Xt * *Mp_EVX;		/* XtP <- t(X)%*%P */
	M_YtP	= M_Yt * *Mp_EVX;		/* YtP <- t(y)%*%P */
#endif
	M_dInv	= Mp_EV->inv2();		/* dInv <- solve(D) */

	/* Under --est, load estimates from --est */
	if (IS_ASSIGNED(est)) {
		loadEstimates(OPT_STRING(est), Ra_estimates, N_pheno, Mp_Xt->row());

		/* Set alpha */
		for (wsUint j=0 ; j<N_pheno ; j++) for (wsUint i=0 ; i<Mp_Xt->row() ; i++)
			Ra_alpha[j][i] = Ra_estimates[j][EMAI_EST_SIZE + i*2];

		/* Set to not failed */
		memset(Ba_isFailed, 0x00, sizeof(char)*N_pheno);

		/* If --est is on, just getting BLUP */
		if (!IS_ASSIGNED(blup)) {
			for (wsUint i=0 ; i<N_pheno ; i++) {
				wsReal	R_curK = Ra_estimates[i][1] / Ra_estimates[i][0];

				wsReal		R_1delta	= W1 / R_curK;
				cDiagMatrix	M_dzz;
				if (IS_ASSIGNED(sampleweight)) {
					cDiagMatrix M_dX(N_samp);
					sseVpC(Ra_sampwgt, R_curK, M_dX.get()[0], N_samp);
					cDiagMatrix& M_dXi	= M_dX.inv();
					M_dzz	= M_dInv + M_dXi;
				} else
					M_dzz		= M_dInv + R_1delta;
//				M_dzz		= M_dInv + R_1delta;
				cDiagMatrix&	M_dzzInv	= M_dzz.inv();
				M_dzz.rem();

				/* XViX and XViy */
				cSymMatrix	M_XtDX		= M_XtP.MMt(M_dzzInv);
				M_XtDX *= R_1delta;
				cSymMatrix	M_XtVX		= M_XtX - M_XtDX;
				M_XtVX *= R_1delta;
				wsReal		R_XtVXdetL;
				cSymMatrix	&M_XtVXinv	= M_XtVX.inv(&R_XtVXdetL);

				cVector		V_curXtY	= M_XtY_t.r2v(i);
				cVector		V_curY		= M_Yt.r2v(i);
				cVector		V_XtDY		= M_XtP.MV(M_dzzInv, V_curY);
				V_XtDY *= R_1delta;
				cVector		V_XtVY		= V_curXtY - V_XtDY;
				V_XtVY *= R_1delta;

				/* BLUE */
				cVector		V_blue		= M_XtVXinv * V_XtVY;
				/* yHat
				 * predicted<- X%*%blue */
				cVector		V_yHat		= Mp_Xt->tV(V_blue);
				/* residual
				 * resid    <- (y - predicted) */
				cVector		V_resid		= V_curY - V_yHat;
				/* BLUP
				 * blup     <- czzInv %*% resid /curK */
#ifdef USE_ED2
				cVector		V_Presid	= V_resid.Mt(*Mp_EVX);
				cStdMatrix	M_PdzzT		= M_dzzInv * *Mp_EVX;
				cVector		V_blup		= V_Presid * M_PdzzT;
#else
				cVector		V_Presid	= V_resid * *Mp_EVX;
				cStdMatrix	M_Pdzz		= *Mp_EVX * M_dzzInv;
				cVector		V_blup		= M_Pdzz * V_Presid;
#endif
				V_blup *= R_1delta;

				memcpy(Ra_calcBLUP[i], V_blup.get(), sizeof(wsReal)*Mp_Xt->col());
				memcpy(Ra_calcPred[i], V_yHat.get(), sizeof(wsReal)*Mp_Xt->col());
			}
		}

		return;
	}

	/* trPhiInv */
	R_trPhiInv = M_dInv.sum();
	if (R_trPhiInv != R_trPhiInv)
		halt("Invalid sample correlation, cannot proceed...");
}

void cEmAiAnalysisV2::_exportBLUP()
{
	if (!OPT_ENABLED(makeblup)) return;
	wsUint N_pheno = M_Yt.row();
	vPheno& V_pheno = Cp_IO->getPhenoInfo();

	char **S_names = NULL;
	wsAlloc(S_names, char *, N_pheno * 2);
	for (wsUint i=0 ; i<N_pheno ; i++) {
		S_names[i*2] = new char[256];
		sprintf(S_names[i*2], "BLUP_%s", V_pheno[i].S_name.c_str());
		S_names[i*2+1] = new char[256];
		sprintf(S_names[i*2+1], "PRED_%s", V_pheno[i].S_name.c_str());
	//	S_names[i*2+1] = new char[256];
	//	sprintf(S_names[i*2+1], "%s.resid", V_pheno[i].S_name.c_str());
	}
	wsReal **Ra_expMat = NULL;
	wsAlloc(Ra_expMat, wsReal*, N_pheno*2);
	for (wsUint i=0 ; i<N_pheno ; i++)
		Ra_expMat[i*2] = Ra_calcBLUP[i];
	for (wsUint i=0 ; i<N_pheno ; i++)
		Ra_expMat[i*2+1] = Ra_calcPred[i];
	wsStrCst S_ext = B_doSpecDcmp?"SD.blup":(OPT_ENABLED(ml)?"NR.blup":
		"AI.blup");
	LOGoutput("BLUP is exported to [%s.%s]\n", OPT_STRING(out), S_ext);	/* 140109 CONFIRMED */
	exportSampleMatrix(S_ext, getIO(), (const char *)Ba_excl, (char)1, N_pheno*2,
		S_names, Ra_expMat);
	LOOP (i, N_pheno*2) delete [] S_names[i];
	DEALLOC(S_names);
	DEALLOC(Ra_expMat);
}

cEmAiAnalysisV2::cEmAiAnalysisV2(cIO *Cp_IO, cStdMatrix &M_inpXt,
	cVector &V_inpY, cMatrix &M_phi, char B_isREML, char *Ba_iExcl,
	cStdMatrix *Mp_eVec/*=NULL*/,
	cDiagMatrix *Mp_eVal/*=NULL*/) : cAnalysis(Cp_IO)
{
	Ra_sampwgt	= NULL;

	/* col(M_Xt) == col(M_Yt) == col(M_phi) */
	if (M_inpXt.col() != V_inpY.size() || V_inpY.size() != M_phi.row())
		halt("Dimension of given matrices is unmatched");
	N_anaSamp	= M_inpXt.col();

	/* Allocate initial failure status */
	wsAlloc(Ba_isFailed, char, Cp_IO->sizePheno());
	for (wsUint i=0 ; i<Cp_IO->sizePheno() ; i++)
		Ba_isFailed[i] = 1;
	M_Yt		= V_inpY.toMt();

	_init(Cp_IO, M_inpXt, M_phi, B_isREML, Ba_iExcl,
		Mp_eVec, Mp_eVal);
}

cEmAiAnalysisV2::cEmAiAnalysisV2(cIO *Cp_IO, cStdMatrix &M_inpXt,
	cStdMatrix &M_inpY, cMatrix &M_phi, char B_isREML, char *Ba_iExcl,
	cStdMatrix *Mp_eVec/*=NULL*/,
	cDiagMatrix *Mp_eVal/*=NULL*/) : cAnalysis(Cp_IO)
{
	/*	 col(M_Xt) == col(M_Yt) == col(M_phi) */
	if (M_inpXt.col() != M_inpY.col())
		halt_fmt(WISARD_SYST_INVL_DIM, M_inpXt.col(), M_inpY.col());
	if (M_inpY.col() != M_phi.row())
		halt_fmt(WISARD_SYST_INVL_DIM, M_inpY.col(), M_phi.row());
	N_anaSamp	= M_inpXt.col();
	/* Allocate initial failure status */
	wsAlloc(Ba_isFailed, char, Cp_IO->sizePheno());
	for (wsUint i=0 ; i<Cp_IO->sizePheno() ; i++)
		Ba_isFailed[i] = 1;
	M_Yt.init(M_inpY.row(), M_inpY.col(), sseMatrixP(M_inpY.row(), M_inpY.col(), M_inpY.get()));

	_init(Cp_IO, M_inpXt, M_phi, B_isREML, Ba_iExcl,
		Mp_eVec, Mp_eVal);
}

cEmAiAnalysisV2::~cEmAiAnalysisV2()
{
	deallocMatrix(Ra_estimates, M_Yt.row(), (void *)1);
	deallocMatrix(Ra_calcBLUP, M_Yt.row(), (void *)1);
	deallocMatrix(Ra_calcPred, M_Yt.row(), (void *)1);
	deallocMatrix(Ra_alpha, M_Yt.row(), (void *)1);
	if (B_dealloc) {
		delete Mp_EV;
		delete Mp_EVX;
		Mp_EV = NULL;
		Mp_EVX = NULL;
	}
//  	if (M_XtXi.row())
//  		delete &M_XtXi;
	M_XtXi.rem();
	DEALLOC(Ba_isFailed);
	DEALLOC(Ra_sampwgt);
	M_Yt.rem();
//	delete &M_Yt;
}

void cEmAiAnalysisV2::run()
{
	wsReal	R_nk		= (wsReal)(N_anaSamp-Mp_Xt->row());
	wsUint	N_pheno		= M_Yt.row();
	cVector V_sig2(N_pheno), V_sig2g(N_pheno);


	/* Do not run if all of B_isFailed == 0 */ {
		wsUint i = 0;
		for ( ; i<N_pheno ; i++)
			if (Ba_isFailed[i] != 0)
				break;
		if (i == N_pheno) return;
	}

	if (!IS_ASSIGNED(emcount)) {
		OPTION().assign("emcount");
		OPTION().FORCE_OPT_NUMBER(emcount);
	}
	if (!IS_ASSIGNED(aithr)) {
		OPTION().assign("aithr");
		OPTION().FORCE_OPT_REAL(aithr);
	}

	/*
	 *
	 * Initial guess
	 *
	 */
	cSymMatrix& M_tmpXtXi = M_XtX.inv();
	M_XtXi.init(M_tmpXtXi);
	delete& M_tmpXtXi;
	cStdMatrix	M_bHat_t= M_XtY_t * M_XtXi;
	//M_XtXi.rem();

	/* Spectral decomposition */
	if (B_doSpecDcmp) {
		spectralDecomposition(V_sig2, V_sig2g);
		_exportBLUP();
		return;
	}

	cStdMatrix	M_yHat	= M_bHat_t * *Mp_Xt;	/* yHat <- X%*%matrix(solve(XX)%*%Xy,ncol=1) */
	cVector		V_btXtY(N_pheno);
	wsReal		*Ra_btXtY = V_btXtY.get();
	for (wsUint i=0 ; i<N_pheno ; i++) {
		cVector V_curBhat	= M_bHat_t.r2v(i);
		cVector V_curXtY	= M_XtY_t.r2v(i);
		Ra_btXtY[i] = V_curBhat.sum(V_curXtY);
	}
 	//wsReal		R_btXtY	= V_bhat.sum(M_XtY_t);
	cVector		V_var	= (V_YtY-V_btXtY) / (R_nk - W1);
	cStdMatrix	M_resid	= M_Yt - M_yHat;		/* residual <- y - yHat */

#ifdef EMAIvalidate
	if (M_Xt.row() == 1)
		M_resid.file("ea.residual0"); // EMAIvalidate
	else
		M_resid.file("ea.residual1"); // EMAIvalidate
#endif

	/* Get the log-likelihood under the assumption of g=0 */
	/* Equation is equivalent to logL of simple linreg.
		* -n/2*log(2*PI)-n/2*log(s^2)-1/(2*s^2)*SSR */
	cVector	V_ssResid	= M_resid.rSS();
// 	wsReal	R_ssResid	= V_resid.ss();
// 	wsReal	R_logL_g0	= -R_ssResid;
// //	wsReal R_sig2_g0 = numeric_limits<wsReal>::infinity();
// 	wsReal R_sVar;
// 	R_sVar		= -R_logL_g0/(wsReal)N_anaSamp;
// 	R_logL_g0	/= W2*SQR(R_sVar);
// 	R_logL_g0	-= N_anaSamp/W2*log(W2*WISARD_PI)
// 		+ N_anaSamp/W2*log(SQR(R_sVar));

	/*
	 *
	 * getInitial part
	 *
	 */
	wsReal		R_ttXmean		= W0;
	cVector		V_sig2_sig2g;
	cStdMatrix	V_stdRes;
	cSymMatrix	M_ttX;
	
	if (!OPT_ENABLED(indep)) {
		cVector		V_sig2_sig2g_sq;
		V_sig2_sig2g	= V_ssResid / (wsReal)(N_anaSamp - 1);	/* sig2_sig2g <- sum(residual^2)/(S-1) */
		V_sig2_sig2g_sq = V_sig2_sig2g.sqrt();
		V_stdRes		= M_resid.rDiv(V_sig2_sig2g_sq);		/* standres   <- residual/sqrt(sig2_sig2g) */
		cIdtMatrix	M_1n(N_anaSamp);
		M_ttX			= *((cSymMatrix *)Mp_phi) - M_1n;
		// ttxMean   ` <- mean(ttX)             # [1]
		R_ttXmean		= M_ttX.sum() / SQR((wsReal)N_anaSamp);
		// ttX        <- ttX - ttxMean         # [S*S]
		M_ttX -= R_ttXmean;
	}

	/* In case of phi-I == 0 */
	if (R_ttXmean == W0) {
		V_stdRes.rem();
		M_ttX.rem();
		LOG("Found that X is independent, sig2g is set to 0\n");

		memcpy(V_sig2.get(), V_var.get(), sizeof(wsReal)*N_pheno);
		V_sig2g.set0();
	} else {
//		MAT_t		Ra_M		= sseEmptyMatrix(N_anaSamp, N_anaSamp);
		wsSym		Ra_ttX		= M_ttX.get();
		// sxy        <- sxx <- 0              # [1]
		cVector	V_sxy(N_pheno);
		wsReal	R_sxx = W0;
		wsReal	**Ra_stdRes = V_stdRes.get();
		V_sxy.set0();
		wsReal	*Ra_sxy = V_sxy.get();

		// for(i in seq(S-1) ) for(j in 2:S){
		//   sxy <- sxy + standres[i]*standres[j]*ttX[i,j]
		//   sxx <- sxx + ttX[i,j]^2
		// }
		wsUint i, j;
		for (i=1 ; i<(N_anaSamp-1) ; i++) {
			/* For 1 */
			for (wsUint k=0 ; k<N_pheno ; k++) {
				Ra_sxy[k] += Ra_stdRes[k][i]*Ra_stdRes[k][0]*Ra_ttX[i][0];
//				Ra_M[i][0] += W1;
			}
			R_sxx += SQR(Ra_ttX[i][0]);

			for (j=1 ; j<i ; j++) {
				for (wsUint k=0 ; k<N_pheno ; k++) {
					Ra_sxy[k] += Ra_stdRes[k][i]*Ra_stdRes[k][j]*Ra_ttX[i][j]*W2;
				}
// 				Ra_M[i][j] += W1;
// 				Ra_M[j][i] += W1;
				R_sxx += SQR(Ra_ttX[i][j])*W2;
			}

			for (wsUint k=0 ; k<N_pheno ; k++) {
				Ra_sxy[k] += SQR(Ra_stdRes[k][i])*Ra_ttX[i][i];
				//Ra_M[i][i] += W1;
			}
			R_sxx += SQR(Ra_ttX[i][i]);
		}
		for (j=0 ; j<(N_anaSamp-1) ; j++) {
			for (wsUint k=0 ; k<N_pheno ; k++)
				Ra_sxy[k] += Ra_stdRes[k][N_anaSamp-1]*Ra_stdRes[k][j]*Ra_ttX[N_anaSamp-1][j];
			R_sxx += SQR(Ra_ttX[N_anaSamp-1][j]);
//			Ra_M[N_anaSamp-1][j] += W1;
		}
		M_ttX.rem();

		// prop       <- sxy/sxx               # [1]
		cVector	V_prop	= V_sxy / R_sxx;

		// if(prop<0.01) prop <- 0.01
		wsReal	*Ra_prop = V_prop.get();
		for (wsUint k=0 ; k<N_pheno ; k++) {
			if (Ra_prop[k] <= REAL_CONST(0.01)) {
				LOGwarn("Inappropriate prop[%g] was adjusted\n",
					Ra_prop[k]);
				Ra_prop[k] = REAL_CONST(0.01);
			}
			// if(prop<0.99) prop <- 0.99
			if (Ra_prop[k] >= REAL_CONST(0.99)) {
				LOGwarn("Inappropriate prop[%g] was adjusted\n",
					Ra_prop[k]);
				Ra_prop[k] = REAL_CONST(0.99);
			}
		}

		// sig2g      <- prop*sig2_sig2g       # [1]
		V_sig2g	= V_prop * V_sig2_sig2g;

		// sig2       <- sig2_sig2g - sig2g    # [1]
		V_sig2	= V_sig2_sig2g - V_sig2g;
	}

	if (!B_doEM)
		LOG("[NOTICE] Skip EM part...\n");

	/* Print initially estimated value */
	for (wsUint k=0 ; k<N_pheno ; k++) {
		wsReal	R_sig2		= V_sig2.get()[k];
		wsReal	R_sig2g		= V_sig2g.get()[k];
		wsReal	R_YtY		= V_YtY.get()[k];
		cVector	V_Y			= M_Yt.r2v(k);
		cVector	V_curXtY	= M_XtY_t.r2v(k);
		cVector	V_curYtP	= M_YtP.r2v(k);
		if (!B_quiet) LOG("sig2 [%g] sig2g [%g]\n", R_sig2, R_sig2g);

		/*
		 *
		 * EM part
		 *
		 */
		if (R_sig2g) {
			if (OPT_ENABLED(ml)) {
				/* EM, ML case */
				if (B_doEM)
					EM_ML(&R_sig2, &R_sig2g, *Mp_EV, M_dInv, M_XtP, V_curYtP,
						M_XtX, V_curXtY, *Mp_Xt, R_YtY, V_Y, k);
				NR_ML(&R_sig2, &R_sig2g, *Mp_EV, M_dInv, M_XtP, V_curYtP,
					M_XtX, V_curXtY, *Mp_Xt, R_YtY, V_Y, k);
			} else {
				/* EM, REML case */
				if (B_doEM)
					EM_REML(&R_sig2, &R_sig2g, *Mp_EV, M_dInv, M_XtP, V_curYtP,
						M_XtX, V_curXtY, *Mp_Xt, R_YtY, V_Y, k);
				else {
					/* Convert */
					wsReal R_tmp = R_sig2;
					R_sig2 = R_sig2g;
					R_sig2g = R_tmp / R_sig2g;
				}
				/* Now R_sig2  = Ra_curTheta */
				/*     R_sig2g = R_curK */
				AI_REML(&R_sig2, &R_sig2g, *Mp_EV, M_dInv, M_XtP, V_curYtP,
					M_XtX, V_curXtY, *Mp_Xt, R_YtY, V_Y, k);
				if (NA(R_sig2) || NA(R_sig2g)) { 
					Ba_isFailed[k] = 1;
					continue;
				}
			}
		} else {
			Ra_estimates[k][0] = R_sig2;
			Ra_estimates[k][1] = R_sig2g;
			Ra_estimates[k][2] = W0; /* logL... */

			wsUint N_inpCol = Mp_Xt->row();
			memcpy(Ra_alpha[k], M_bHat_t.get()[k], sizeof(wsReal)*N_inpCol);
			memset(Ra_calcBLUP[k], 0x00, sizeof(wsReal)*N_anaSamp);
			memcpy(Ra_calcPred[k], M_yHat.get()[k], sizeof(wsReal)*N_anaSamp);

			Ba_isFailed[k] = 0;
		}
	}
	/* Do spectral decomposition on failure */ {
		wsUint i = 0;
		for ( ; i<N_pheno ; i++)
			if (Ba_isFailed[i])
				break;
		if (i != N_pheno) {
			B_doSpecDcmp = 1;
			spectralDecomposition(V_sig2, V_sig2g);
		}
	}
	_exportBLUP();
}

void cEmAiAnalysisV2::EM_ML(wsReal *Rp_sig2, wsReal *Rp_sig2g, cDiagMatrix &M_d,
	cDiagMatrix &M_dInv, cStdMatrix &M_XtP, cVector &V_YtP,
	cSymMatrix &M_XtX, cVector &V_XtY, cStdMatrix &M_Xt, wsReal R_YtY,
	cVector &V_Y, wsUint N_idxPheno)
{
	wsReal R_cur2	= *Rp_sig2;
	wsReal R_cur2g	= *Rp_sig2g;

	for (int N_loop=0 ; N_loop<OPT_NUMBER(emcount) ; N_loop++) {
		cDiagMatrix&	M_dL	= M_d * R_cur2g;
		if (IS_ASSIGNED(sampleweight)) {
			cDiagMatrix M_dX(M_d.row());
			sseVpC(Ra_sampwgt, R_cur2, M_dX.get()[0], M_d.row());
			M_dL += M_dX;
		} else
			M_dL += R_cur2;
		cDiagMatrix&	M_dLinv	= M_dL.inv();
		M_dL.rem();
		
		/* curK    <- cursig2/cursig2g */
		wsReal		R_curK		= R_cur2 / R_cur2g;
		/* czz     <- 1/curK*I + phiInv */
		wsReal		R_invK		= W1 / R_curK;
		cDiagMatrix	M_dzz;
		if (IS_ASSIGNED(sampleweight)) {
			cDiagMatrix M_dX(M_d.row());
			sseVpC(Ra_sampwgt, R_curK, M_dX.get()[0], M_d.row());
			cDiagMatrix& M_dXi	= M_dX.inv();
			M_dzz	= M_dInv + M_dXi;
		} else
			M_dzz		= M_dInv + R_invK;
		cDiagMatrix&	M_dzzInv	= M_dzz.inv();
		M_dzz.rem();

		/* XViX and XViy */
		cSymMatrix	M_XtDX		= M_XtP.MMt(M_dzzInv);
		cVector		V_XtDY		= M_XtP.MV(M_dzzInv, V_YtP);
		M_XtDX *= R_invK;
		V_XtDY *= R_invK;
		cSymMatrix	M_XtVX		= M_XtX - M_XtDX;
		M_XtDX.rem();
		cVector		V_XtVY		= V_XtY - V_XtDY;
		V_XtDY.rem();
		M_XtVX *= R_invK;
		V_XtVY *= R_invK;
		wsReal		R_XtVXdetL;
		cSymMatrix	&M_XtVXinv	= M_XtVX.inv(&R_XtVXdetL);
		M_XtVX.rem();

		/* BLUE */
		cVector		V_blue		= M_XtVXinv * V_XtVY;
		M_XtVXinv.rem();
		/* yHat
			* predicted<- X%*%blue */
		cVector		V_yHat		= M_Xt.tV(V_blue);
		/* residual
			* resid    <- (y - predicted) */
		cVector		V_resid		= V_Y - V_yHat;
		/* BLUP
			* blup     <- czzInv %*% resid * invK */
#ifdef USE_ED2
		cVector		V_Presid	= V_resid.Mt(*Mp_EVX);
		cStdMatrix	M_PdzzT		= M_dzzInv * *Mp_EVX;
		cVector		V_blup		= V_Presid * M_PdzzT;
		M_PdzzT.rem();
#else
		cVector		V_Presid	= V_resid * *Mp_EVX;
		cStdMatrix	M_Pdzz		= *Mp_EVX * M_dzzInv;
		cVector		V_blup		= M_Pdzz * V_Presid;
		M_Pdzz.rem();
#endif
		V_Presid.rem();
		V_blup *= R_invK;

		cVector		V_rmb		= V_resid - V_blup;
#ifdef USE_ED2
		cVector		V_residP	= V_resid.Mt(*Mp_EVX);
#else
		cVector		V_residP	= V_resid * *Mp_EVX;
#endif
		cDiagMatrix	M_DzDDz		= M_dzzInv.MMt(M_dInv);
		M_dzzInv.rem();
		wsReal		R_bPb		= V_residP.qf(M_DzDDz) * SQR(R_invK);
		V_resid.rem();

		/* sig2  <- cursig2 - cursig2^2*mean(diag(W)) + mean(err^2) */
		wsReal		R_Nsig2		= R_cur2 - SQR(R_cur2)*M_dLinv.tr()/(wsReal)N_anaSamp
			+ V_rmb.ss() / (wsReal)N_anaSamp;
		/* sig2g <- cursig2g- cursig2g^2*mean(diag(W%*%phi))+as.numeric(t(blup)%*%phiInv%*%blup/S) */
		M_dLinv *= *Mp_EV;
		wsReal		R_Nsig2g	= R_cur2g - SQR(R_cur2g)*M_dLinv.tr()/(wsReal)N_anaSamp
			+ R_bPb / (wsReal)N_anaSamp;
		M_dLinv.rem();

		wsReal R_cvgRes = sqrt((R_cur2g-R_Nsig2g)*(R_cur2g-R_Nsig2g) +
			(R_cur2-R_Nsig2)*(R_cur2-R_Nsig2));
		if (!B_quiet) LOG("(sig2) %g (sig2g) %g (tV) %g\n", R_Nsig2, R_Nsig2g, R_cvgRes);
		/* if( sqrt((cursig2g-sig2g)^2 + (cursig2-sig2)^2) < threshold) return(c(sig2g,sig2/sig2g)) */
		R_cur2	= R_Nsig2;
		R_cur2g	= R_Nsig2g;
	}

	*Rp_sig2 = R_cur2;
	*Rp_sig2g = R_cur2g;
}

void cEmAiAnalysisV2::NR_ML(wsReal *Rp_sig2, wsReal *Rp_sig2g, cDiagMatrix &M_d,
	cDiagMatrix &M_dInv, cStdMatrix &M_XtP, cVector &V_YtP,
	cSymMatrix &M_XtX, cVector &V_XtY, cStdMatrix &M_Xt, wsReal R_YtY,
	cVector &V_Y, wsUint N_idxPheno)
{
	wsReal	*Ra_est		= Ra_estimates[N_idxPheno];
	wsReal	R_thrNR		= REAL_CONST(1e-5);
	wsReal	R_curTheta	= *Rp_sig2g;
	wsReal	R_curK		= *Rp_sig2 / *Rp_sig2g;

	char	B_converged	= 0;
	for (wsUint N_loop=0 ; N_loop<30 ; N_loop++) {
		//wsReal		R_invK	= W1 / R_curK;
		cDiagMatrix	M_dL;
		if (IS_ASSIGNED(sampleweight)) {
			cDiagMatrix M_dX(M_d.row());
			sseVpC(Ra_sampwgt, R_curK, M_dX.get()[0], M_d.row());
			M_dL = *Mp_EV + M_dX;
		} else
			M_dL = *Mp_EV + R_curK;
		cDiagMatrix	&M_dLinv= M_dL.inv();

		cStdMatrix	M_dLXp	= M_dLinv.mt(M_XtP);
		cSymMatrix	M_XDX	= M_XtP.MMt(M_dLinv);
		cSymMatrix	&M_XDXi	= M_XDX.inv();
		wsReal		M_dLdL	= M_dL.detL();
		M_dL.rem();
		if (B_converged == 1) {
			cVector	V_XDy	= M_XtP.MV(M_dLinv, V_YtP);
			cVector	V_beta	= M_XDXi * V_XDy;
			*Rp_sig2		= R_curTheta * R_curK;
			cVector	V_pred	= M_Xt.tV(V_beta);
			cVector V_resid	= V_Y - V_pred;
#ifdef USE_ED2
			cVector	V_resP	= V_resid.Mt(*Mp_EVX);
#else
			cVector	V_resP	= V_resid * *Mp_EVX;
#endif
			wsReal	R_rVIr	= V_resP.qf(M_dLinv);
			Ra_est[2]		= -REAL_CONST(0.5)
				* ( M_dLdL + (wsReal)N_anaSamp*log(R_rVIr) );
			cVector	V_DdLPr	= Mp_EV->MV(M_dLinv, V_resP);
#ifdef USE_ED2
			cVector	V_blup	= Mp_EVX->tV(V_DdLPr);
#else
			cVector	V_blup	= *Mp_EVX * V_DdLPr;
#endif

			/* Export BLUP */
// 			if (OPT_ENABLED(makeblup)) {
// 				LOG("Export BLUP to [%s.NR.blup]\n", OPT_STRING(out));
// 				exportSampleMatrix("NR.blup", getIO(), Ba_incl, Na_sampIdx, 2,
// 					"NR.blup", V_blup.get(), "NR.resid", V_resid.get());
// 			}

			memcpy(Ra_calcBLUP[N_idxPheno], V_blup.get(), sizeof(wsReal)*Mp_Xt->col());
			memcpy(Ra_calcPred[N_idxPheno], V_pred.get(), sizeof(wsReal)*Mp_Xt->col());
			memcpy(Ra_alpha[N_idxPheno], V_beta.get(), sizeof(wsReal)*Mp_Xt->row());

			Ra_est[1]	= R_curTheta;
			Ra_est[0]	= R_curTheta * R_curK;
			Ba_isFailed[N_idxPheno]	= 0;
			break;
			//Ra_estimates[2] = R_logL;
		}
		M_XDX.rem();
		cSymMatrix	M_dRdc	= M_dLXp.MMt(M_XDXi);
		M_XDXi.rem();
		cSymMatrix	M_Dst	= M_dLinv - M_dRdc;

		cVector		V_yP	= V_YtP * M_Dst;
		wsReal		R_yPy	= V_yP.sum(V_YtP);
		wsReal		R_yP2y	= V_yP.sum(V_yP);
		wsReal		R_yP3y	= V_yP.qf(M_Dst);

		wsReal		R_trV	= M_dLinv.tr();
		wsReal		R_trVV	= M_dLinv.tr2();
		M_dLinv.rem();

		/* uTheta  <- -0.5*S/curTheta + 0.5*yPy/curTheta^2 */
		wsReal		R_uThta	= -REAL_CONST(0.5)*(wsReal)N_anaSamp/R_curTheta
			+ REAL_CONST(0.5)*R_yPy/SQR(R_curTheta);
		/* uK      <- -0.5*trVI + 0.5/curTheta*yPPy */
		wsReal		R_uK	= -REAL_CONST(0.5)*R_trV
			+ REAL_CONST(0.5)/R_curTheta*R_yP2y;

		/*
		I11     <- S/(2*curTheta^2) - yPy/curTheta^3
		I12     <- -1/(2*curTheta^2)*yPPy
		I22     <- 0.5*sum(diag(VI%*%VI)) - 1/curTheta * (t(Py)%*%P%*%Py)[1,1]
		*/
		wsReal		R_I11	= (wsReal)N_anaSamp/(W2 * SQR(R_curTheta))
			- R_yPy/CUBE(R_curTheta);
		wsReal		R_I12	= -R_yP2y/(W2*SQR(R_curTheta));
		wsReal		R_I22	= REAL_CONST(0.5)*R_trVV - R_yP3y/R_curTheta;

		// 	aa        <- try(IIInv     <- solve(II),silent=T)
		// 	if(class(aa)!='try-error') {
		wsReal R_det = R_I11*R_I22 - R_I12*R_I12;
		wsReal R_updTheta=W0, R_updK=W0;
		if (R_det) {
			// updates   <- matrix(c(curTheta,curK),ncol=1) - IIInv%*%matrix(c(uTheta,uK),ncol=1)
			// 
			// [ curTheta ] - [ I22/R_det  -I12/R_det ] [ R_newTheta ]
			// [ curK     ]   [ -I12/R_det I11/R_det  ] [ R_newK     ]
			// 
			// Thus,
			// updates[1] <- curTheta - (I22*R_newTheta)/R_det + (I12*R_newK)/R_det
			// updates[2] <- curK + (I12*R_newTheat)/R_det - (I11*R_newK)/R_det
			R_updTheta	= R_curTheta
				- (R_I22*R_uThta)/R_det + (R_I12*R_uK)/R_det;
			R_updK		= R_curK
				+ (R_I12*R_uThta)/R_det - (R_I11*R_uK)/R_det;	
		} else {
			// print("Errors!!")
			// break
			LOG("NR-ML error, not invertiable II\n");
			break;
		}

		// cvgRes   <- sqrt( sum( (updates-c(curTheta,curK))^2 ) )
		// ( [ updates[1] ] - [ curTheta ] )^2
		// ( [ updates[2] ]   [ curK     ] )^2
		wsReal R_cvgRes = sqrt( SQR(R_updTheta-R_curTheta) + SQR(R_updK-R_curK) );
		if (R_cvgRes != R_cvgRes)
			break;
		/* Set failed if NA */
		if (NA(R_updK) || NA(R_updTheta)) {
			Ba_isFailed[N_idxPheno] = 1;
			break;
		}
		if (!B_quiet) LOG("(sig2) %g (sig2g) %g (tV) %g\n",
			R_updTheta*R_updK, R_updTheta, R_cvgRes);
		// 				if( sqrt( sum( (updates-c(curTheta,curK))^2 ) ) < threshold) {
		if (R_cvgRes < R_thrNR)
			B_converged = 1;

		// 				curTheta <-updates[1]
		// 				curK     <-updates[2]
		R_curTheta	= R_updTheta;
		R_curK		= R_updK;
	}

	if (Ba_isFailed[N_idxPheno])
		LOG("[%d]th phenotype failed to converge...\n", N_idxPheno+1);
}

void cEmAiAnalysisV2::EM_REML(wsReal *Rp_sig2, wsReal *Rp_sig2g, cDiagMatrix &M_d,
	cDiagMatrix &M_dInv, cStdMatrix &M_XtP, cVector &V_YtP,
	cSymMatrix &M_XtX, cVector &V_XtY, cStdMatrix &M_Xt, wsReal R_YtY,
	cVector &V_Y, wsUint N_idxPheno)
{
	/*
	czz    <- 1/curK*diag(S) + phiInv
	czzInv <- solve(czz)

	dzzInv <- 1 / (diag(S)/curK + dInv)
	*/
	wsReal	R_nk		= (wsReal)(N_anaSamp-M_Xt.row());
	wsReal	R_curK		= *Rp_sig2 / *Rp_sig2g;
	wsReal	R_curTheta	= *Rp_sig2g;
//	wsReal	R_prevLogL;
	wsUint	N_totLoop	= OPT_NUMBER(emcount);

	wsReal	R_newTheta	= R_curTheta;
	wsReal	R_newK		= R_curK;
	for (wsUint N_loop=1 ; N_loop<N_totLoop ; N_loop++) {
		wsReal		R_invK		= W1 / R_curK;
		cDiagMatrix	M_dzz;

		if (IS_ASSIGNED(sampleweight)) {
			cDiagMatrix M_dX(M_d.row());
			sseVpC(Ra_sampwgt, R_curK, M_dX.get()[0], M_d.row());
			cDiagMatrix& M_dXi	= M_dX.inv();
			M_dzz	= M_dInv + M_dXi;
			delete& M_dXi;
		} else {
			cDiagMatrix& M_dzzT	= M_dInv + R_invK;
			M_dzz.init(M_dzzT);
			delete& M_dzzT;
		}
		cDiagMatrix	&M_dzzInv	= M_dzz.inv();
//		delete &M_dzz;

		/* XViX and XViy */
		cSymMatrix	M_XtDX		= M_XtP.MMt(M_dzzInv);
		cVector		V_XtDY		= M_XtP.MV(M_dzzInv, V_YtP);
		M_XtDX *= R_invK;
		V_XtDY *= R_invK;
		cSymMatrix	M_XtVX		= M_XtX - M_XtDX;
		M_XtDX.rem();
		cVector		V_XtVY		= V_XtY - V_XtDY;
		V_XtDY.rem();
		M_XtVX *= R_invK;
		V_XtVY *= R_invK;
		wsReal		R_XtVXdetL;
		cSymMatrix	&M_XtVXinv	= M_XtVX.inv(&R_XtVXdetL);
		M_XtVX.rem();
		wsReal		R_YtVY		= (R_YtY - V_YtP.qf(M_dzzInv) * R_invK) * R_invK;

		/* BLUE */
		cVector		V_blue		= M_XtVXinv * V_XtVY;
		/* yHat
			* predicted<- X%*%blue */
		cVector		V_yHat		= M_Xt.tV(V_blue);
		/* residual
			* resid    <- (y - predicted) */
		cVector		V_resid		= V_Y - V_yHat;
		/* BLUP
		 * blup     <- czzInv %*% resid * invK */
#ifdef USE_ED2
		cVector		V_Presid	= V_resid.Mt(*Mp_EVX);
		cStdMatrix	M_PdzzT		= M_dzzInv * *Mp_EVX;
		cVector		V_blup		= V_Presid * M_PdzzT;
		M_PdzzT.rem();
#else
		cVector		V_Presid	= V_resid * *Mp_EVX;
		cStdMatrix	M_Pdzz		= *Mp_EVX * M_dzzInv;
		cVector		V_blup		= M_Pdzz * V_Presid;
		M_Pdzz.rem();
#endif
		V_Presid.rem();
		V_blup *= R_invK;

		cVector		V_rmb		= V_resid - V_blup;
		V_resid.rem();

		/* yPy      <- yVIy - t(XVIy)%*%XVIXInv%*%XVIy */
		wsReal		R_yPy		= R_YtVY - V_XtVY.qf(M_XtVXinv);
		V_XtVY.rem();

		/* newTheta <- yPy /(S-q) */
		R_newTheta	= R_yPy / R_nk;

		/* c22   <- czzInv + czzInv%*%X%*%XVIXInv%*%t(X)%*%czzInv/curK^2 */
		cStdMatrix	M_dXP		= M_dzzInv.mt(M_XtP);
		cSymMatrix	M_DXXVXXD	= M_dXP.MMt(M_XtVXinv);
		delete &M_XtVXinv;
		M_dXP.rem();
		M_DXXVXXD *= SQR(R_invK);
		cSymMatrix	M_c22		= M_DXXVXXD + M_dzzInv;
		delete &M_dzzInv;
		/* pcp      <- sum(diag(phiInv%*%c22%*%phiInv)) */
		wsReal		R_pcp		= M_dInv.trMMt(M_c22);
		M_c22.rem();


		/* lDet_pK  <- determinant(phi+diag(S)*curK)$modulus[1] */
		cDiagMatrix	&M_dK		= M_d + R_curK;
		wsReal		R_dKdetL	= M_dK.detL();
		delete &M_dK;

		/* newK     <- curK  - (trphiInv - pcp)*curK^2/S + sum( (resid-blup)^2 )/curTheta/S */
		R_newK		= R_curK
			- (R_trPhiInv - R_pcp)*SQR(R_curK) / (wsReal)N_anaSamp
			+ V_rmb.ss() / R_curTheta / (wsReal)N_anaSamp;
			

		/* newlogL  <- -0.5*( lDet_XVIX + lDet_pK + (S-q)*log(curTheta) + yPy/curTheta ) */
		wsReal		R_logL		= -REAL_CONST(0.5) * (fabs(R_XtVXdetL) + R_dKdetL
			+ R_nk*log(R_curTheta) + R_yPy/R_curTheta);
		wsReal R_tV = sqrt(SQR(R_newTheta-R_curTheta) + SQR(R_newK-R_curK));
		if (N_loop && R_tV < OPT_REAL(aithr)) {
			*Rp_sig2 = R_newTheta;
			*Rp_sig2g = R_newK;
		}
//		R_prevLogL	= R_logL;

		R_curTheta	= R_newTheta;
		R_curK		= R_newK;
		if (!B_quiet) LOG("(sig2) %g (sig2g) %g (logL) %g (tV) %g\n",
			R_curTheta*R_curK, R_curTheta, R_logL, R_tV);
	}
	*Rp_sig2	= R_newTheta;
	*Rp_sig2g	= R_newK;
}

void cEmAiAnalysisV2::AI_REML(wsReal *Rp_curTheta, wsReal *Rp_curK,
	cDiagMatrix &M_d, cDiagMatrix &M_dInv, cStdMatrix &M_XtP, cVector &V_YtP,
	cSymMatrix &M_XtX, cVector &V_XtY, cStdMatrix &M_Xt, wsReal R_YtY,
	cVector &V_Y, wsUint N_idxPheno)
{
	wsReal	R_curK		= *Rp_curK;
	wsReal	R_curTheta	= *Rp_curTheta;
	wsReal	R_nk		= (wsReal)(N_anaSamp-M_Xt.row());
	char	B_converged	= 0;
	wsReal	R_lastQQ	= WISARD_NAN;
	wsReal	*Ra_est		= Ra_estimates[N_idxPheno];

	for (wsUint N_loop=0 ; N_loop<30 || OPT_ENABLED(nostop) ; N_loop++) {
		wsReal		R_1delta	= W1 / R_curK;
		cDiagMatrix	M_dzz;

		if (IS_ASSIGNED(sampleweight)) {
			cDiagMatrix M_dX(M_d.row());
			sseVpC(Ra_sampwgt, R_curK, M_dX.get()[0], M_d.row());
			cDiagMatrix& M_dXi	= M_dX.inv();
			M_dzz	= M_dInv + M_dXi;
			delete &M_dXi;
		} else {
			cDiagMatrix& M_dzzT	= M_dInv + R_1delta;
			M_dzz.init(M_dzzT);
			delete& M_dzzT;
		}
		cDiagMatrix	&M_dzzInv	= M_dzz.inv();
		M_dzz.rem();
		//delete &M_dzz;

// 		cDiagMatrix	&M_dzz		= M_dInv + R_1delta;
// 		cDiagMatrix	&M_dzzInv	= M_dzz.inv();
// 		delete &M_dzz;

		/* XViX and XViy */
		cSymMatrix	M_XtDX		= M_XtP.MMt(M_dzzInv);
		cVector		V_XtDY		= M_XtP.MV(M_dzzInv, V_YtP);
		M_XtDX *= R_1delta;
		V_XtDY *= R_1delta;
		cSymMatrix	M_XtVX		= M_XtX - M_XtDX;
		cVector		V_XtVY		= V_XtY - V_XtDY;
		M_XtVX *= R_1delta;
		V_XtVY *= R_1delta;
		wsReal		R_XtVXdetL;
		cSymMatrix	&M_XtVXinv	= M_XtVX.inv(&R_XtVXdetL);

		/* BLUE */
		cVector		V_blue		= M_XtVXinv * V_XtVY;
		/* yHat
			* predicted<- X%*%blue */
		cVector		V_yHat		= M_Xt.tV(V_blue);
		/* residual
			* resid    <- (y - predicted) */
		cVector		V_resid		= V_Y - V_yHat;
		/* BLUP
			* blup     <- czzInv %*% resid /curK */
#ifdef USE_ED2
		cVector		V_Presid	= V_resid.Mt(*Mp_EVX);
		cStdMatrix	M_PdzzT		= M_dzzInv * *Mp_EVX;
		cVector		V_blup		= V_Presid * M_PdzzT;
#else
		cVector		V_Presid	= V_resid * *Mp_EVX;
		cStdMatrix	M_Pdzz		= *Mp_EVX * M_dzzInv;
		cVector		V_blup		= M_Pdzz * V_Presid;
#endif
		V_blup *= R_1delta;

		if (B_converged || ((N_loop+1)==30 && OPT_ENABLED(forceconv))) {
			/* Warning if forcely converged */
			if ((N_loop+1)==30 && OPT_ENABLED(forceconv) && !B_converged)
				LOGwarn("[%d]th phenotype was forcely converged\n", N_idxPheno);
			Ba_isFailed[N_idxPheno]	= 0;

			/* lDet_pK  <- determinant(phi+diag(S)*curK)$modulus[1] */
			cDiagMatrix	&M_dK		= M_d + R_curK;
			wsReal		R_dKdetL	= M_dK.detL();
			delete &M_dK;

			/* Export BLUP */
// 			if (OPT_ENABLED(makeblup)) {
// 				LOG("Export BLUP to [%s.AI.blup]\n", OPT_STRING(out));
// 				exportSampleMatrix("AI.blup", getIO(), Ba_incl, Na_sampIdx, 2,
// 					"AI.blup", V_blup.get(), "AI.resid", V_resid.get());
// 			}

			memcpy(Ra_calcBLUP[N_idxPheno], V_blup.get(), sizeof(wsReal)*Mp_Xt->col());
			memcpy(Ra_calcPred[N_idxPheno], V_yHat.get(), sizeof(wsReal)*Mp_Xt->col());
			memcpy(Ra_alpha[N_idxPheno], V_blue.get(), sizeof(wsReal)*Mp_Xt->row());

			//wsReal		R_YtVY		= (R_YtY - V_YtP.qf(M_dzzInv) * R_1delta) * R_1delta;
			/* yPy      <- yVIy - t(XVIy)%*%XVIXInv%*%XVIy */
			//wsReal		R_yPy		= R_YtVY - V_XtVY.qf(M_XtVXinv);
			/* newlogL  <- -0.5*( lDet_XVIX + lDet_pK + (S-q)*log(curTheta) + QPQ[1,1]*curTheta ) */
			wsReal		R_logL		= -REAL_CONST(0.5) * (R_XtVXdetL + R_dKdetL
				+ R_nk * log(R_curTheta) + R_lastQQ*R_curTheta);

			Ra_est[1] = R_curTheta;
			Ra_est[0] = R_curTheta * R_curK;
			Ra_est[2] = R_logL;
			delete &M_dzzInv;
			delete& M_XtVXinv;
			break;
		}

		/* Py <- (resid-blup)/curK */
		cVector		V_Py		= V_resid - V_blup;
		V_Py *= R_1delta;

		/* Q  <- matrix(c( resid/curTheta, Py), nrow=S) */
		cStdMatrix	M_Q(2, N_anaSamp);
		wsMat		Ra_Q		= M_Q.get();
		cVector		V_Q0		= V_resid / R_curTheta;
		memcpy(Ra_Q[0], V_Q0.get(), sizeof(wsReal)*N_anaSamp);
		memcpy(Ra_Q[1], V_Py.get(), sizeof(wsReal)*N_anaSamp);

		/* QQ <- t(Q)%*%Q */
		cSymMatrix	M_QQ		= M_Q.Mt();
		/* QX <- t(Q)%*%X */
		cStdMatrix	M_QX		= M_Q.Mt(M_Xt);

		/* Qp */
#ifdef USE_ED2
		cStdMatrix	M_Qp		= M_Q.Mt(*Mp_EVX);
#else
		cStdMatrix	M_Qp		= M_Q * *Mp_EVX;
#endif

		/* QVIQ <- QQ/curK - t(Q)%*%czzInv%*%Q/curK^2 */
		cSymMatrix	M_QDQ		= M_Qp.MMt(M_dzzInv);
		M_QDQ *= R_1delta;
		M_QQ -= M_QDQ;
		M_QQ *= R_1delta;
		/* QX/curK - t(Q)%*%czzInv%*%X/curK^2 */
		cStdMatrix	M_QD		= M_Qp * M_dzzInv;
		cStdMatrix	M_QDX		= M_QD.Mt(M_XtP);
		M_QDX *= R_1delta;
		M_QX -= M_QDX;
		M_QX *= R_1delta;
		/* QPQ  <- QVIQ - QVIX%*%XVIXInv%*%t(QVIX) */
		cSymMatrix	M_QXQ		= M_QX.MMt(M_XtVXinv);
		M_QQ -= M_QXQ;
		/* Utheta<- -0.5*( (S-q)/curTheta  - QPQ[1,1]) */
		wsReal		R_Utheta	= -REAL_CONST(0.5)
			* (R_nk/R_curTheta - M_QQ.get()[0][0]);

		/* c22   <- czzInv + czzInv%*%X%*%XVIXInv%*%t(X)%*%czzInv/curK^2 */
		cStdMatrix	M_dXP		= M_dzzInv.mt(M_XtP);
		cSymMatrix	M_DXXVXXD	= M_dXP.MMt(M_XtVXinv);
		delete &M_XtVXinv;
		M_DXXVXXD *= SQR(R_1delta);
		cSymMatrix	M_c22		= M_DXXVXXD + M_dzzInv;
		delete &M_dzzInv;
		/* pcp      <- sum(diag(phiInv%*%c22%*%phiInv)) */
		wsReal		R_pcp		= M_dInv.trMMt(M_c22);

		/* UK      <- -0.5*( trphiInv - pcp - t(Py)%*% Py /(curTheta) ) */
		wsReal		R_Py2		= V_Py.ss();
		wsReal		R_UK		= -REAL_CONST(0.5)
			* (R_trPhiInv - R_pcp - R_Py2/R_curTheta);

		/* updates <- matrix(c(curTheta,curK),ncol=1)
			+ 2*curTheta*solve(QPQ)%*%matrix(c(Utheta,UK),ncol=1) */
// 		cSymMatrix	M_QPQinv	= M_QQ.inv();
// 		MAT_t		Ra_QPQinv	= M_QPQinv.get();
// 		wsReal		R_upd1		= R_curTheta
// 			+ W2*R_curTheta
// 			* (Ra_QPQinv[0][0]*R_Utheta + Ra_QPQinv[1][0]*R_UK);
// 		wsReal		R_upd2		= R_curK
// 			+ W2*R_curTheta
// 			* (Ra_QPQinv[1][0]*R_Utheta + Ra_QPQinv[1][1]*R_UK);
		wsReal **Ra_QVIQ = M_QQ.get();
		wsReal R_det = 2*R_curTheta / (Ra_QVIQ[0][0]*Ra_QVIQ[1][1]-Ra_QVIQ[1][0]*Ra_QVIQ[1][0]);
		wsReal R_upd1 = R_curTheta + R_det*(Ra_QVIQ[1][1]*R_Utheta-Ra_QVIQ[1][0]*R_UK);
		wsReal R_upd2 = R_curK + R_det*(-Ra_QVIQ[1][0]*R_Utheta+Ra_QVIQ[0][0]*R_UK);

		/* if(any(updates<=0)) { */
		if (R_upd1 <= W0 || R_upd2 <= W0) {
				/* updates[1] <- QPQ[1,1]*curTheta^2/(S-q) */
			R_upd1 = M_QQ.get()[0][0] * SQR(R_curTheta) / R_nk;
			/* updates[2] <- curK  - (trphiInv - pcp)*curK^2/S
				+ sum( (resid-blup)^2 )/curTheta/S */
			cVector		V_rmb		= V_resid - V_blup;
			R_upd2 = R_curK
				- (R_trPhiInv - R_pcp)*SQR(R_curK)/(wsReal)N_anaSamp
				+ V_rmb.ss() / R_curTheta / (wsReal)N_anaSamp;
		}

		/* if( sqrt( sum( (updates-c(curTheta,curK))^2 ) ) < threshold) { */
//		wsReal R_adet = sqrt(SQR(R_upd1-R_curTheta) + SQR(R_upd2-R_curK));
		wsReal R_adet = sqrt(
			SQR(R_upd1-R_curTheta)
			+
			SQR(R_upd2-R_curK)
			) / sqrt(SQR(R_upd1)+SQR(R_upd2));
		if (R_adet < OPT_REAL(aithr)) {
			B_converged = 1;
			R_lastQQ = M_QQ.get()[0][0];
		}
		R_curTheta	= R_upd1;
		R_curK		= R_upd2;

		/* Stop if NA */
		if (NA(R_curK) || NA(R_curTheta)) {
			Ba_isFailed[N_idxPheno] = 1;
			break;
		}

		if (!B_quiet) LOG("(sig2) %g (sig2g) %g (tV) %g\n",
			R_curTheta*R_curK, R_curTheta, R_adet);
	}

	if (Ba_isFailed[N_idxPheno])
		LOG("[%d]th phenotype failed to converge...\n", N_idxPheno+1);
}

void cEmAiAnalysisV2::AI_REMLv2(wsReal *Rp_curTheta, wsReal *Rp_curK,
	cDiagMatrix &M_d, cDiagMatrix &M_dInv, cStdMatrix &M_XtP, cVector &V_YtP,
	cSymMatrix &M_XtX, cVector &V_XtY, cStdMatrix &M_Xt, wsReal R_YtY,
	cVector &V_Y, wsUint N_idxPheno)
{
	wsReal	R_curK		= *Rp_curK;
	wsReal	R_curTheta	= *Rp_curTheta;
	wsReal	R_nk		= (wsReal)(N_anaSamp-M_Xt.row());
	char	B_converged	= 0;
	wsReal	R_lastQQ	= WISARD_NAN;
	wsReal	*Ra_est		= Ra_estimates[N_idxPheno];

	for (wsUint N_loop=0 ; N_loop<30 || OPT_ENABLED(nostop) ; N_loop++) {
		cDiagMatrix	M_di	= M_d.v_kv(R_curK);
		cVector		V_di(M_di.get()[0], M_di.row(), 1);

		/* XVIX but no /curK */
		cSymMatrix	M_XdX	= M_XtP.MMt(M_di);
		cSymMatrix	M_XdiX	= M_XtX - M_XdX;
		wsReal		R_XtVXdetL;
		cSymMatrix&	M_XVXi	= M_XdiX.inv(&R_XtVXdetL);
		R_XtVXdetL -= log(R_curK);
		/* XVIy but no /curK */
		cVector		V_Xdy	= M_XtP.MV(M_di, V_YtP);
		cVector		V_Xdiy	= V_XtY - V_Xdy;
		/* (XVIX)^-1XVIy cancels curK */
		cVector		V_blue	= M_XVXi * V_Xdiy;
		cVector		V_yHatP	= V_blue * M_XtP;
		cVector		Y_yP	= V_YtP - V_yHatP;
#ifdef USE_ED2
		cStdMatrix	M_dP	= M_di * *Mp_EVX;
		cVector		V_blup	= Y_yP * M_dP;
#else
		cStdMatrix	M_dP	= *Mp_EVX * M_di;
		cVector		V_blup	= M_dP * Y_yP;
#endif
		cVector		V_yHat	= M_Xt.tV(V_blue);
		cVector		V_resid	= V_Y - V_yHat;

		if (B_converged) {
			Ba_isFailed[N_idxPheno]	= 0;

			/* lDet_pK  <- determinant(phi+diag(S)*curK)$modulus[1] */
			cDiagMatrix	&M_dK		= M_d + R_curK;
			wsReal		R_dKdetL	= M_dK.detL();
			delete &M_dK;

			/* Export BLUP */
// 			if (OPT_ENABLED(makeblup)) {
// 				LOG("Export BLUP to [%s.AI.blup]\n", OPT_STRING(out));
// 				exportSampleMatrix("AI.blup", getIO(), Ba_incl, Na_sampIdx, 2,
// 					"AI.blup", V_blup.get(), "AI.resid", V_resid.get());
// 			}

			memcpy(Ra_calcBLUP[N_idxPheno], V_blup.get(), sizeof(wsReal)*Mp_Xt->col());
			memcpy(Ra_calcPred[N_idxPheno], V_yHat.get(), sizeof(wsReal)*Mp_Xt->col());
			memcpy(Ra_alpha[N_idxPheno], V_blue.get(), sizeof(wsReal)*Mp_Xt->row());

			//wsReal		R_YtVY		= (R_YtY - V_YtP.qf(M_dzzInv) * R_1delta) * R_1delta;
			/* yPy      <- yVIy - t(XVIy)%*%XVIXInv%*%XVIy */
			//wsReal		R_yPy		= R_YtVY - V_XtVY.qf(M_XtVXinv);
			/* newlogL  <- -0.5*( lDet_XVIX + lDet_pK + (S-q)*log(curTheta) + QPQ[1,1]*curTheta ) */
			wsReal		R_logL		= -REAL_CONST(0.5) * (R_XtVXdetL + R_dKdetL
				+ R_nk * log(R_curTheta) + R_lastQQ*R_curTheta);

			Ra_est[1] = R_curTheta;
			Ra_est[0] = R_curTheta * R_curK;
			Ra_est[2] = R_logL;
			delete &M_XVXi;

			break;
		}

		cStdMatrix	M_XtD	= M_XtP * M_di;
		cVector		V_c22	= M_XtD.diag(M_XVXi);
		V_c22 += V_di;
		V_c22 *= R_curK;
		cVector		V_dInv(M_dInv.get()[0], M_dInv.row(), 1);
		cVector		V_dInv2 = V_dInv.sq();
		wsReal		R_pcp	= V_c22.sum(V_dInv2);

		/* Py <- (resid-blup)/curK */
		cVector		V_Py		= V_resid - V_blup;
		V_Py /= R_curK;

		/* Q  <- matrix(c( resid/curTheta, Py), nrow=S) */
		cStdMatrix	M_Q(2, N_anaSamp);
		wsMat		Ra_Q		= M_Q.get();
		cVector		V_Q0		= V_resid / R_curTheta;
		memcpy(Ra_Q[0], V_Q0.get(), sizeof(wsReal)*N_anaSamp);
		memcpy(Ra_Q[1], V_Py.get(), sizeof(wsReal)*N_anaSamp);

		/* QQ <- t(Q)%*%Q */
		cSymMatrix	M_QQ		= M_Q.Mt();
		/* QX <- t(Q)%*%X */
		cStdMatrix	M_QX		= M_Q.Mt(M_Xt);

		/* Qp */
#ifdef USE_ED2
		cStdMatrix	M_Qp		= M_Q.Mt(*Mp_EVX);
#else
		cStdMatrix	M_Qp		= M_Q * *Mp_EVX;
#endif

		/* QVIQ <- QQ/curK - t(Q)%*%czzInv%*%Q/curK^2 */
		cSymMatrix	M_QDQ		= M_Qp.MMt(M_di);
		M_QQ -= M_QDQ;

		/* QX/curK - t(Q)%*%czzInv%*%X/curK^2 */
		cStdMatrix	M_QDX		= M_Qp.MMt(M_di, M_XtP);
		M_QX -= M_QDX;
		cSymMatrix	M_QXQ		= M_QX.MMt(M_XVXi);
		delete &M_XVXi;
		M_QQ -= M_QXQ;
		M_QQ *= W1/R_curK;

		/* Utheta<- -0.5*( (S-q)/curTheta  - QPQ[1,1]) */
		wsReal		R_Utheta	= -REAL_CONST(0.5)
			* (R_nk/R_curTheta - M_QQ.get()[0][0]);

		/* c22   <- czzInv + czzInv%*%X%*%XVIXInv%*%t(X)%*%czzInv/curK^2 */
// 		cStdMatrix	M_dXP		= M_dzzInv.mt(M_XtP);
// 		cSymMatrix	M_DXXVXXD	= M_dXP.MMt(M_XtVXinv);
// 		delete &M_XtVXinv;
// 		M_DXXVXXD *= SQR(R_1delta);
// 		cSymMatrix	M_c22		= M_DXXVXXD + M_dzzInv;
// 		delete &M_dzzInv;
		/* pcp      <- sum(diag(phiInv%*%c22%*%phiInv)) */
//		wsReal		R_pcp		= M_dInv.trMMt(M_c22);

		/* UK      <- -0.5*( trphiInv - pcp - t(Py)%*% Py /(curTheta) ) */
		wsReal		R_Py2		= V_Py.ss();
		wsReal		R_UK		= -REAL_CONST(0.5)
			* (R_trPhiInv - R_pcp - R_Py2/R_curTheta);

		/* updates <- matrix(c(curTheta,curK),ncol=1)
			+ 2*curTheta*solve(QPQ)%*%matrix(c(Utheta,UK),ncol=1) */
// 		cSymMatrix	M_QPQinv	= M_QQ.inv();
// 		MAT_t		Ra_QPQinv	= M_QPQinv.get();
// 		wsReal		R_upd1		= R_curTheta
// 			+ W2*R_curTheta
// 			* (Ra_QPQinv[0][0]*R_Utheta + Ra_QPQinv[1][0]*R_UK);
// 		wsReal		R_upd2		= R_curK
// 			+ W2*R_curTheta
// 			* (Ra_QPQinv[1][0]*R_Utheta + Ra_QPQinv[1][1]*R_UK);
		wsReal **Ra_QVIQ = M_QQ.get();
		wsReal R_det = 2*R_curTheta / (Ra_QVIQ[0][0]*Ra_QVIQ[1][1]-Ra_QVIQ[1][0]*Ra_QVIQ[1][0]);
		wsReal R_upd1 = R_curTheta + R_det*(Ra_QVIQ[1][1]*R_Utheta-Ra_QVIQ[1][0]*R_UK);
		wsReal R_upd2 = R_curK + R_det*(-Ra_QVIQ[1][0]*R_Utheta+Ra_QVIQ[0][0]*R_UK);

		/* if(any(updates<=0)) { */
		if (R_upd1 <= W0 || R_upd2 <= W0) {
				/* updates[1] <- QPQ[1,1]*curTheta^2/(S-q) */
			R_upd1 = M_QQ.get()[0][0] * SQR(R_curTheta) / R_nk;
			/* updates[2] <- curK  - (trphiInv - pcp)*curK^2/S
				+ sum( (resid-blup)^2 )/curTheta/S */
			cVector		V_rmb		= V_resid - V_blup;
			R_upd2 = R_curK
				- (R_trPhiInv - R_pcp)*SQR(R_curK)/(wsReal)N_anaSamp
				+ V_rmb.ss() / R_curTheta / (wsReal)N_anaSamp;
		}

		/* if( sqrt( sum( (updates-c(curTheta,curK))^2 ) ) < threshold) { */
//		wsReal R_adet = sqrt(SQR(R_upd1-R_curTheta) + SQR(R_upd2-R_curK));
		wsReal R_adet = sqrt(
			SQR(R_upd1-R_curTheta)
			+
			SQR(R_upd2-R_curK)
			) / sqrt(SQR(R_upd1)+SQR(R_upd2));
		if (R_adet < OPT_REAL(aithr)) {
			B_converged = 1;
			R_lastQQ = M_QQ.get()[0][0];
		}
		R_curTheta	= R_upd1;
		R_curK		= R_upd2;
		if (!B_quiet) LOG("(sig2) %g (sig2g) %g (tV) %g\n",
			R_curTheta*R_curK, R_curTheta, R_adet);
	}

	if (Ba_isFailed[N_idxPheno])
		LOG("[%d]th phenotype failed to converge...\n", N_idxPheno+1);
}

int cEmAiAnalysisV2::getEst()
{
	wsUint i;
	cTimer t;

	pverbose("Get estimates from model...\n");
	wsUint	N_inpCol	= Mp_Xt->row();
	wsUint	N_pheno		= M_Yt.row();

	/* lambda */
	for (wsUint k=0 ; k<N_pheno ; k++) {
		wsReal *Ra_est	= Ra_estimates[k];
		wsReal R_l = Ra_est[1] / Ra_est[0];

		cDiagMatrix	&M_dL		= *Mp_EV * R_l;
		if (IS_ASSIGNED(sampleweight)) {
			cDiagMatrix M_dX(M_dL.row(), Ra_sampwgt, 1);
			M_dL += M_dX;
		} else
			M_dL += W1;

		cDiagMatrix &M_dInvL	= M_dL.inv();
		delete &M_dL;
		cStdMatrix	M_dLXtP		= M_dInvL.mt(M_XtP);
		cSymMatrix	M_XDX		= M_XtP.MMt(M_dInvL);
		cSymMatrix	&M_XDXinv	= M_XDX.inv();
		cSymMatrix	M_DXXDXXD	= M_dLXtP.MMt(M_XDXinv);
		cSymMatrix	M_Ppart		= M_dInvL - M_DXXDXXD;
		delete &M_dInvL;
		cStdMatrix	M_Pphi		= M_Ppart * *Mp_EV;

		wsReal R_2sig2 = W2 * Ra_est[0];

		wsReal R_11 = (wsReal)(N_anaSamp-N_inpCol) / (W2*SQR(Ra_est[0]));
		wsReal R_12 = M_Pphi.tr() / R_2sig2;
		wsReal R_22 = REAL_CONST(0.5) * M_Pphi.tr2();

		/* Get Var(bHat_i) = sig^2 * (Xt %*% V^-1 %*% X)^-1_ii */
		wsMat	Ra_XtViX_inv = M_XDXinv.get();
		for (i=0 ; i<N_inpCol ; i++) {
			Ra_est[EMAI_EST_SIZE+(i*2)]		= Ra_alpha[k][i];
			Ra_est[EMAI_EST_SIZE+(i*2)+1]	= Ra_est[0] * Ra_XtViX_inv[i][i];
		}
		delete &M_XDXinv;

		/* Take inverse */
		wsReal R_det = R_11*R_22 - R_12*R_12;
		if (R_det == W0) {
			wsMat M = sseMatrix(2, 2);
			wsMat Mi = NULL;
			M[0][0] = R_11;
			M[0][1] = M[1][0] = R_12;
			M[1][1] = R_22;
			if (SVDinverse(M, 2, &Mi) == false)
				halt("Can't get determinant of varcov mat of sigma, maybe error?");
			R_11 = Mi[0][0];
			R_22 = Mi[1][1];
			R_12 = Mi[1][0];
		} else {
			wsReal tmp = R_11;
			R_11 = R_22 / R_det;
			R_22 = tmp / R_det;
			R_12 = -R_12 / R_det;
		}

		Ra_est[3] = R_11;	///< Var(sig2)
		Ra_est[4] = R_11*SQR(Ra_est[1]/Ra_est[0]) +
			W2*R_12*Ra_est[1] + R_22*Ra_est[0]*Ra_est[0];	///< Var(sig2g)
		Ra_est[5] = Ra_est[1] / (Ra_est[0]+Ra_est[1]); ///< h^2

		/* Get Var(h^2)
		 *     [    Var(s2)      Cov(s2, sg2) ]      [ s2 / (sg2+sg)^2   ]
		 * X = [                              ], A = [                   ]
		 *     [ Cov(s2, sg2)      Var(sg2)   ]      [ -sg2 / (sg2+sg)^2 ]
		 *
		 * Var(h^2) = A' %*% X %*% A
		 */
	//	wsReal R_sqS2Sg2	= pow(Ra_est[0]+Ra_est[1], W2);
	// 	wsReal R_a1			= Ra_est[0] / R_sqS2Sg2;
	// 	wsReal R_a2			= -Ra_est[1] / R_sqS2Sg2;
	// 	Ra_est[6] = R_11*R_a1*R_a1 + W2*R_12*R_a1*R_a2 +
	// 		R_22*R_a2*R_a2;
		Ra_est[6] = R_22*pow(Ra_est[0]/(Ra_est[1]+Ra_est[0]), REAL_CONST(4.0));
		// Var(sig2+sig2g) = Var(sig2) + Var(sig2g) + 2*(a* sig2g/sig2 + b*sig2)
		Ra_est[7] = Ra_est[3]+Ra_est[4] + W2*
			(R_11*Ra_est[1]/Ra_est[0] + R_12*Ra_est[0]);
	}

	return 0;
}

void cEmAiAnalysisV2::spectralDecomposition(cVector &V_sig2, cVector &V_sig2g)
{
	wsUint	N_inpCol	= Mp_Xt->row();
	wsUint	N_pheno		= M_Yt.row();

	if (OPT_ENABLED(nospecdcmp)) {
		LOGwarn("--nospecdcmp defined, use initial guess\n");
		for (wsUint i=0 ; i<N_pheno ; i++) {
			wsReal *Ra_est = Ra_estimates[i];
			Ra_est[0] = V_sig2.get()[i];
			Ra_est[1] = V_sig2g.get()[i];
			Ra_est[2] = WISARD_NAN;
		}
	} else {
		// Look for singular values
		xRealSort *Xa_Xi = NULL; /* checked */
		wsAlloc(Xa_Xi, xRealSort, N_anaSamp);

		wsReal *Ra_eval = Mp_EV->get()[0];
		for (wsUint i=0; i<N_anaSamp; i++) {
			Xa_Xi[i].i = i;
			if (Ra_eval[i] != W0)
				Xa_Xi[i].V = W1 / Ra_eval[i];
			else
				Xa_Xi[i].V = W0;
		}
		qsort(Xa_Xi, N_anaSamp, sizeof(xRealSort), sort_real);
		notice("Step 1/7...\r");

		// [nxn].[t(v)] 
		//	printf("O matrix[%x]\n", O);60
#ifdef USE_ED2
		cStdMatrix	M_EVXt = Mp_EVX->transpose();
		cSymMatrix	M_phiInvSD = M_EVXt.MMt(M_dInv);
		M_EVXt.rem();
#else
		cSymMatrix	M_phiInvSD = Mp_EVX->MMt(M_dInv);
#endif
		notice("Step 2/7...\r");

		/* Get optimum using fmin method of R() */
	//  wsReal R_ratio = (wsReal)Brent_fmin(0.01, 2.0, 0,
	//  	logL_spectralDecomp, &X_sd, pow(numeric_limits<double>::epsilon(), 0.25));
	//  wsReal R_delta = SQR(R_ratio);
	
		/* Sp <- X%*%solve(t(X)%*%X)%*%t(X) */
		cStdMatrix	M_1		= M_XtXi * *Mp_Xt;
		cStdMatrix	M_S		= Mp_Xt->tM(M_1);
		wsMat		Ra_S	= M_S.get();

		/* S <- I - Sp # diag(N) - Sp */
		for (wsUint i=0 ; i<N_anaSamp ; i++) {
			for (wsUint j=0 ; j<N_anaSamp ; j++) {
				if (i == j)
					Ra_S[i][j] = W1-Ra_S[i][j];
				else
					Ra_S[i][j] = -Ra_S[i][j];
			}
		}

		wsMat Ra_At = NULL; /* Checked, P^t */
		wsUint N_NminusP = N_anaSamp - N_inpCol;
		char B_res = 0;
		// At <- t(s$vectors[,1:nminusp])
	#ifdef USE_ED2
		wsReal *Ra_Eval_S = EIGENDECOMPOSITION(Ra_S, N_anaSamp, &Ra_At);
	#else
		wsReal *Ra_Eval_S = eigenDecomp(Ra_S, N_anaSamp, &B_res, &Ra_At, 1);
	#endif
		if (B_res)
			halt("Eigendecomposition failed");
		sseFree(Ra_Eval_S); /* Checked */
		M_S.rem();
		cStdMatrix M_At(N_anaSamp, N_anaSamp, Ra_At, MATDEL_NONE);
		//sseUnmat(Ra_S, N_anaSamp);

		/* Extract top N-p eigenvectors */

		notice("Step 3/7...\r");

		// A <- s$vectors[,1:nminusp]
		// AZKZA <- At%*%K%*%t(At)
		cSymMatrix M_AtKA = M_At.MMt(*Mp_phi);
	// 	wsReal **Ra_AtKA = sseMMMt(Ra_At, N_NminusP, N_anaSamp,
	// 		Mp_phi->get(), N_anaSamp, N_anaSamp, Ra_At, N_NminusP, N_anaSamp); /* Checked */
		// azkza <- eigen(AZKZA)
		// lambda <- azkza$values
		wsReal **Ra_Evec_AtKA_t = NULL; /* Checked, P^t */
	#ifdef USE_ED2
		wsReal *Ra_lambda = EIGENDECOMPOSITION(M_AtKA.get(), N_NminusP,
			&Ra_Evec_AtKA_t); /* Checked */
	#else
		wsReal *Ra_lambda = eigenDecomp(M_AtKA.get(), N_NminusP, &B_res,
			&Ra_Evec_AtKA_t, 1); /* Checked */
	#endif
		if (B_res)
			halt("Eigendecomposition failed");
		M_AtKA.rem();

		// UR_t <- t(azkza$vectors)%*%At
		wsReal **Ra_URt = sseMpM(Ra_Evec_AtKA_t, N_NminusP, N_NminusP,
			Ra_At, N_NminusP, N_anaSamp); /* Checked */
		sseUnmat(Ra_Evec_AtKA_t, N_NminusP);
		sseUnmat(Ra_At, N_anaSamp);

		// eta <- UR_t%*%dat$y		# [(n-p) %*% 1]
		cSymMatrix	M_XphiX		= Mp_Xt->MMt(M_phiInvSD);
		// XphiXinv <- solve(XphiX)
		cSymMatrix	&M_XphiXi	= M_XphiX.inv();
		for (wsUint k=0 ; k<N_pheno ; k++) {
			wsReal	R_delta		= W1;
			wsReal	*Ra_Y		= M_Yt.get()[k];
			wsReal	*Ra_est		= Ra_estimates[k];
			cVector	V_Y			= M_Yt.r2v(k);
			wsReal	**Ra_eta	= sseMpMt(Ra_URt, N_NminusP, N_anaSamp,
				&Ra_Y, 1, N_anaSamp); /* Checked */

			notice("Step 4/7...\r");
			/* Calculate common part #2 = \eta_s^2
			 * 
			 * \Sigma_s(\eta_s^2 / (\lambda_s+\delta)^2)
			 * -----------------------------------------
			 *  \Sigma_s(\eta_s^2 / (\lambda_s+\delta)) */
			wsReal R_ssEta		= W0;
			wsReal R_ssLambda	= W0;
			wsReal R_sLambda	= W0;
			for (wsUint i=0 ; i<N_NminusP ; i++) {
				R_ssEta		+= SQR(Ra_eta[i][0]);
				R_sLambda	+= Ra_lambda[i];
				R_ssLambda	+= SQR(Ra_lambda[i]);
			}

			wsUint	N_iter = 100;
			wsReal	R_logL = WISARD_NAN;
			char	B_forceStop = 0;
			while (--N_iter) {
				wsReal R_denom	= W0;
				wsReal R_div	= W0;
				wsReal R_denom2	= W0;
				wsReal R_denom3	= W0;

				wsReal R_E2LD	= W0;
				wsReal R_1LD	= W0;
				wsReal R_1LD2	= W0;
				wsReal R_1XD	= W0;
				wsReal R_1XD2	= W0;
				wsReal R_E2LD2	= W0;
				wsReal R_E2LD3	= W0;

				for (wsUint i=0 ; i<N_NminusP ; i++) {
					wsReal R_E2	= SQR(Ra_eta[i][0]);
					wsReal R_LD	= Ra_lambda[i] + R_delta;
					wsReal R_LD2= SQR(R_LD);
					wsReal R_XD	= Xa_Xi[i].V + R_delta;
					wsReal R_XD2= SQR(R_XD);

					R_1LD		+= W1 / R_LD;
					R_1LD2		+= W1 / R_LD2;
					R_1XD		+= W1 / R_XD;
					R_1XD2		+= W1 / R_XD2;
					R_E2LD		+= R_E2 / R_LD;
					R_E2LD2		+= R_E2 / R_LD2;
					R_E2LD3		+= R_E2 / (R_LD2*R_LD);

					R_denom += SQR(Ra_eta[i][0]) / SQR(Ra_lambda[i]+R_delta);
					R_div += SQR(Ra_eta[i][0]) / (Ra_lambda[i]+R_delta);
					R_denom2 += Ra_lambda[i]+R_delta;
				}
				for (wsUint i=N_NminusP ; i<N_anaSamp ; i++) {
					wsReal R_XD	= Xa_Xi[i].V + R_delta;
					wsReal R_XD2= SQR(R_XD);
					R_1XD		+= W1 / R_XD;
					R_1XD2		+= W1 / R_XD2;
				}

				/* Get f(x)/f'(x) */
				wsReal R_fp, R_fpp;
		//		wsReal R_trueDenom;
				wsReal R_newDelta;
				if (OPT_ENABLED(ml)) {
					for (wsUint i=0 ; i<N_NminusP ; i++) {
						wsReal R_s = Ra_eval[i]+R_delta;
						R_denom2 += R_s;
						R_denom3 += log(R_s);
					}
					R_logL = REAL_CONST(0.5) * ((wsReal)N_anaSamp *
						(log((wsReal)N_anaSamp/((W2*WISARD_PI)))
						- W1
						- log(R_denom)) - R_denom3);
		// 			R_trueDenom = (wsReal)(N_anaSamp)/W2
		// 				* R_div / R_denom - REAL_CONST(0.5)/R_denom2;

					wsReal R_EL2EL = R_E2LD2/R_E2LD;
					R_fp = REAL_CONST(0.5) * ((wsReal)N_anaSamp*R_EL2EL - R_1XD);
					R_fpp = REAL_CONST(0.5)*R_1XD2 - (wsReal)N_NminusP*(
						R_E2LD3/R_E2LD - REAL_CONST(0.5)*SQR(R_EL2EL)
						);
				} else {
					wsReal R_denom4 = W0;
					for (wsUint i=0 ; i<N_NminusP ; i++) {
						wsReal R_s = Ra_lambda[i]+R_delta;
						R_denom4 += W1/R_s;
						R_denom2 += R_s;
						R_denom3 += log(R_s);
					}
					R_logL = REAL_CONST(0.5) * ((wsReal)N_NminusP *
						(log((wsReal)N_NminusP/((W2*WISARD_PI)))
						- W1
						- log(R_div)) - R_denom3);
		// 			R_trueDenom = (wsReal)(N_NminusP)/W2
		// 				* R_denom / R_div - REAL_CONST(0.5)*R_denom4;
					wsReal R_EL2EL = R_E2LD2/R_E2LD;
					R_fp	= REAL_CONST(0.5) * ((wsReal)N_NminusP*R_EL2EL - R_1LD);
					R_fpp	= REAL_CONST(0.5)*R_1LD2 - (wsReal)N_NminusP*(
						R_E2LD3/R_E2LD - REAL_CONST(0.5)*SQR(R_EL2EL)
						);
				}
				R_newDelta = R_delta - R_fp/R_fpp;

				if (fabs(R_newDelta - R_delta) < 1e-5 || B_forceStop) {
					if (B_forceStop == 0)
						R_newDelta = R_delta;
					break;
				}
				if (R_newDelta < REAL_CONST(0.01) || R_newDelta > REAL_CONST(100.0)) {
					if (R_newDelta < REAL_CONST(0.01)) R_newDelta = REAL_CONST(0.01);
					if (R_newDelta > REAL_CONST(100.0)) R_newDelta = REAL_CONST(100.0);
					LOGwarn("Too distant(%g) delta value, it will be clamped\n", R_newDelta);
					B_forceStop = 1;
				}
				R_delta = R_newDelta;
			}
			if (N_iter == 0)
				halt("Cannot converged");

			notice("Step 5/7...\r");
			/* Now we have accurate delta */

			// R <- sum(eta^2/(lambda+delta))		[1]
			wsReal R_R = W0;
			for (wsUint i=0 ; i<N_NminusP ; i++)
				R_R += SQR(Ra_eta[i][0]) / (Ra_lambda[i]+R_delta);

			// temp <- t(X)%*%solve(K)				[#inp*#samp]
			cVector		V_XphiY		= Mp_Xt->MV(M_phiInvSD, V_Y);
			// betahat <- XphiXinv%*%temp%*%dat$y	[#inp*#pheno]
			cVector		V_betaHat	= M_XphiXi * V_XphiY;
			wsReal		*Ra_betaHat	= V_betaHat.get();
			for (wsUint i=0 ; i<N_inpCol ; i++)
				Ra_alpha[k][i] = Ra_betaHat[i];
			sseUnmat(Ra_eta, N_NminusP);

			if (OPT_ENABLED(ml))
				Ra_est[1] = R_R / (wsReal)N_anaSamp; /* sig2g */
			else
				Ra_est[1] = R_R / (wsReal)N_NminusP; /* sig2g */
			Ra_est[0] = Ra_est[1] * R_delta; /* sig2 */
			Ra_est[2] = R_logL;		/* logL */

			LOG("Spectral decomposition succeeded : sig2 (%g) sig2g (%g) logL (%g)\n",
				Ra_est[0], Ra_est[1], R_logL);//EMAIvalidate
		}
		delete &M_XphiXi;
		sseUnmat(Ra_URt, N_NminusP);
		DEALLOC(Xa_Xi);
		sseFree(Ra_lambda);

		notice("Step 6/7...\r");
	}

	/* Now get BLUP */
	for (wsUint k=0 ; k<N_pheno ; k++) {
		if (Ba_isFailed[k] == 0) continue;
		wsReal	*Ra_est		= Ra_estimates[k];
		wsReal	R_curK		= Ra_est[0]/Ra_est[1];

		cVector	V_Y			= M_Yt.r2v(k);
		cVector	V_curXtY	= M_XtY_t.r2v(k);
		cVector	V_curYtP	= M_YtP.r2v(k);

		wsReal		R_1delta	= W1 / R_curK;

		cDiagMatrix	M_dzz;
		if (IS_ASSIGNED(sampleweight)) {
			cDiagMatrix M_dX(V_Y.size());
			sseVpC(Ra_sampwgt, R_curK, M_dX.get()[0], V_Y.size());
			cDiagMatrix& M_dXi	= M_dX.inv();
			M_dzz	= M_dInv + M_dXi;
		} else {
			cDiagMatrix& M_dzzT	= M_dInv + R_1delta;
			M_dzzT.setDontDealloc();
			M_dzz.init(M_dzzT.row(), M_dzzT.get()[0]);
			delete& M_dzzT;
		}

//		cDiagMatrix	&M_dzz		= M_dInv + R_1delta;
		cDiagMatrix	&M_dzzInv	= M_dzz.inv();
		//delete &M_dzz;

		/* XViX and XViy */
		cSymMatrix	M_XtDX		= M_XtP.MMt(M_dzzInv);
		cVector		V_XtDY		= M_XtP.MV(M_dzzInv, V_curYtP);
		M_XtDX *= R_1delta;
		V_XtDY *= R_1delta;
		cSymMatrix	M_XtVX		= M_XtX - M_XtDX;
		cVector		V_XtVY		= V_curXtY - V_XtDY;
		M_XtVX *= R_1delta;
		V_XtVY *= R_1delta;
		wsReal		R_XtVXdetL;
		cSymMatrix	&M_XtVXinv	= M_XtVX.inv(&R_XtVXdetL);

		/* BLUE */
		cVector		V_blue		= M_XtVXinv * V_XtVY;
		delete &M_XtVXinv;
		/* yHat
		 * predicted<- X%*%blue */
		cVector		V_yHat		= Mp_Xt->tV(V_blue);
		/* residual
		 * resid    <- (y - predicted) */
		cVector		V_resid		= V_Y - V_yHat;
		/* BLUP
		 * blup     <- czzInv %*% resid /curK */
#ifdef USE_ED2
		cVector		V_Presid	= V_resid.Mt(*Mp_EVX);
		cStdMatrix	M_PdzzT		= M_dzzInv * *Mp_EVX;
		delete &M_dzzInv;
		cVector		V_blup		= V_Presid * M_PdzzT;
#else
		cVector		V_Presid	= V_resid * *Mp_EVX;
		cStdMatrix	M_Pdzz		= *Mp_EVX * M_dzzInv;
		delete &M_dzzInv;
		cVector		V_blup		= M_Pdzz * V_Presid;
#endif
		V_blup *= R_1delta;

		memcpy(Ra_calcBLUP[k], V_blup.get(), sizeof(wsReal)*Mp_Xt->col());
		memcpy(Ra_calcPred[k], V_yHat.get(), sizeof(wsReal)*Mp_Xt->col());

// 		if (OPT_ENABLED(makeblup)) {
// 			LOG("Export BLUP to [%s.SD.blup]\n", OPT_STRING(out));
// 			exportSampleMatrix("SD.blup", getIO(), Ba_incl, Na_sampIdx,
// 			2,
// 			"SD.blup", Ra_calcBLUP,
// 			"SD.resid", V_resid.get());
// 		}
	}
	memset(Ba_isFailed, 0x00, sizeof(char)*N_pheno);


	LOG("Spectral decomposition complete\n");

	//	halt("phiInv should be given");
//	deallocMatrix(Ra_Xt, N_inpCol);
}

} // End namespace ONETOOL
