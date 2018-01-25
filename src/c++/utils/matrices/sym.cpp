#include <math.h>
#include "utils/matrices/sym.h"
#include "global/Rconn.h"

namespace ONETOOL {

wsSym sseSsubsetRect(wsSym Ra_m, wsUintCst N_sz, const char *Ba_filt,
	char N_incVal, wsUintCst N_szAfter)
{
	wsSym Ra_ret = sseSymMat(N_szAfter);

	wsUint i, j, _i, _j;
	for (i=_i=0 ; i<N_sz ; i++) {
		if (Ba_filt[i] != N_incVal) continue;

		for (j=_j=0 ; j<=i ; j++) {
			if (Ba_filt[j] != N_incVal) continue;

			Ra_ret[_i][_j++] = Ra_m[i][j];
		}
		_i++;
	}
	if (_i != N_szAfter)
		halt("Invalid subset making, desired [%d] but subsetted [%d]",
		N_szAfter, _i);

	return Ra_ret;
}

/*
SS <- function(s, t)
{
  if (dim(s)[2] != dim(t)[1])
    stop("Can't do SS")

  sz <- dim(s)[1]
  ret <- matrix(0, nrow=sz, ncol=sz)
  tmp <- c()
  for (i in 1:sz) {
   if (i!=sz) {
     b <- (i+1):sz
     tmp[b] <- s[b,i]
	}

    for (a in 1:i) {
      b <- 1:(a-1)
      ret[i,a] <- ret[i,a]+sum(s[i,b]*t[a,b])
      ret[i,b] <- ret[i,b]+s[i,a]*t[a,b]
      ret[i,a] <- ret[i,a]+s[i,a]*t[a,a]
    }
    if (i!=sz)
    for (a in (i+1):sz) {
      ai <- tmp[a]

      b <- 1:i
      B <- 1:(i-1)
      ret[i,a] <- ret[i,a]+sum(s[i,B]*t[a,B])
      ret[i,b] <- ret[i,b]+ai*t[a,b]
      b <- (i+1):sz
      ret[i,a] <- ret[i,a]+sum(tmp[b]*t[a,b])
      ret[i,b] <- ret[i,b]+ai*t[a,b]
    }
  }
  ret
}
 */
wsMat sseSpS(wsSym A, wsUint N_s1, wsSym B, wsUint N_s2)
{
	if (N_s1 != N_s2)
		halt_fmt(WISARD_SYST_INVL_DIM, "matrix multiplication", N_s1, N_s2);
#ifdef USE_SYM
	wsUint	i, a, b;
	wsMat	C = sseEmptyMatrix(N_s1, N_s1);
	wsReal*	Ra_iCol = NULL;
	sseMalloc(Ra_iCol, wsReal, N_s1);
	
	for (i=0 ; i<N_s1 ; i++) {
		/* Set i */
		for (a=i ; a<N_s1 ; a++)
			Ra_iCol[a] = A[a][i];

		for (a=0 ; a<=i ; a++) {
			wsUint N_med = 0;
#ifdef USE_SSE
			sse_t sse_Aia = sseSet(A[i][a]);
			sse_t sse_C = sseSet(0.0);
			N_med = getMed(a);
			for(b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_Bab = (sse_t *)(B[a] + b);
				sse_t *sse_Ai = (sse_t *)(A[i] + b);

				sse_C = sseAdd(sse_C, sseMul(*sse_Ai, *sse_Bab));

				sse_t *sse_Ci = (sse_t *)(C[i] + b);
				*sse_Ci = sseAdd(*sse_Ci, sseMul(sse_Aia, *sse_Bab));
			}
			sseSum(sse_C, C[i][a]);
#endif
			for(b=N_med ; b<a ; b++) {
				C[i][a] += A[i][b] * B[a][b];
				C[i][b] += A[i][a] * B[a][b];
			}
			C[i][a] += A[i][a] * B[a][a];
		}

		/* When a == i+1 */
		if (i < (N_s1-1)) {
			a = i+1;
			wsUint N_med = 0;
#ifdef USE_SSE
			sse_t sse_Aia = sseSet(Ra_iCol[a]);
			sse_t sse_C = sseSet(0.0);
			N_med = getMed(i+1);
			for(b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_Bab = (sse_t *)(B[a] + b);
				sse_t *sse_Ai = (sse_t *)(A[i] + b);

				sse_C = sseAdd(sse_C, sseMul(*sse_Ai, *sse_Bab));

				sse_t *sse_Ci = (sse_t *)(C[i] + b);
				*sse_Ci = sseAdd(*sse_Ci, sseMul(sse_Aia, *sse_Bab));
			}
			sseSum(sse_C, C[i][a]);
#endif
			for(b=N_med ; b<=i ; b++) {
				C[i][a] += A[i][b] * B[a][b];
				C[i][b] += Ra_iCol[a] * B[a][b];
			}
			C[i][a] += A[a][i] * B[a][a];
		}

		for(a=i+2 ; a<N_s1 ; a++) {
			wsUint N_med = 0;
#ifdef USE_SSE
			sse_t sse_Aia = sseSet(Ra_iCol[a]);
			sse_t sse_C = sseSet(0.0);
			N_med = getMed(i+1);
			for(b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_Bab = (sse_t *)(B[a] + b);
				sse_t *sse_Ai = (sse_t *)(A[i] + b);

				sse_C = sseAdd(sse_C, sseMul(*sse_Ai, *sse_Bab));

				sse_t *sse_Ci = (sse_t *)(C[i] + b);
				*sse_Ci = sseAdd(*sse_Ci, sseMul(sse_Aia, *sse_Bab));
			}
			sseSum(sse_C, C[i][a]);
#endif
			for(b=N_med ; b<=i ; b++) {
				C[i][a] += A[i][b] * B[a][b];
				C[i][b] += Ra_iCol[a] * B[a][b];
			}
#ifdef USE_SSE
			wsUint N_s = getMed(i+2);
			N_med = getMed(a);
			if (N_s > N_med) N_s = N_med;
			for(b=i+1 ; b<N_s ; b++) {
				C[i][a] += Ra_iCol[b] * B[a][b];
				C[i][b] += Ra_iCol[a] * B[a][b];
			}
			sse_C = sseSet(0.0);
			for(b=N_s ; b<N_med ; b+=sseJmp) {
				sse_t *sse_Bab = (sse_t *)(B[a] + b);
				sse_t *sse_Ai = (sse_t *)(Ra_iCol + b);

				sse_C = sseAdd(sse_C, sseMul(*sse_Ai, *sse_Bab));

				sse_t *sse_Ci = (sse_t *)(C[i] + b);
				*sse_Ci = sseAdd(*sse_Ci, sseMul(sse_Aia, *sse_Bab));
			}
			sseAppendSum(sse_C, C[i][a]);

			for(b=N_med ; b<a ; b++) {
				C[i][a] += Ra_iCol[b] * B[a][b];
				C[i][b] += Ra_iCol[a] * B[a][b];
			}
#else
			for(b=i+1 ; b<a ; b++) {
				C[i][a] += Ra_iCol[b] * B[a][b];
				C[i][b] += Ra_iCol[a] * B[a][b];
			}
#endif
			C[i][a] += A[a][i] * B[a][a];
		}
	}

	sseFree(Ra_iCol);
	return C;
#else
	return sseMpMt(A, N_s1, N_s1, B, N_s2, N_s2);
#endif
}

wsSym sseSpS(wsSym A, wsUint N_s1)
{
#ifdef USE_SYM
	wsUint	i, a, b;
	wsSym	C = sseEmptySymMat(N_s1);
	wsReal*	Ra_iCol = NULL;
	sseMalloc(Ra_iCol, wsReal, N_s1);

	for (i=0 ; i<N_s1 ; i++) {
		/* Set i */
		for (a=i ; a<N_s1 ; a++)
			Ra_iCol[a] = A[a][i];

		for (a=0 ; a<=i ; a++) {
			wsUint N_med = 0;
#ifdef USE_SSE
			sse_t sse_Aia = sseSet(A[i][a]);
			sse_t sse_C = sseSet(0.0);
			N_med = getMed(a);
			for(b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_Bab = (sse_t *)(A[a] + b);
				sse_t *sse_Ai = (sse_t *)(A[i] + b);

				sse_C = sseAdd(sse_C, sseMul(*sse_Ai, *sse_Bab));

				sse_t *sse_Ci = (sse_t *)(C[i] + b);
				*sse_Ci = sseAdd(*sse_Ci, sseMul(sse_Aia, *sse_Bab));
			}
			sseSum(sse_C, C[i][a]);
#endif
			for(b=N_med ; b<a ; b++) {
				C[i][a] += A[i][b] * A[a][b];
				C[i][b] += A[i][a] * A[a][b];
			}
			C[i][a] += A[i][a] * A[a][a];
		}

		/* When a == i+1 */
		if (i < (N_s1-1)) {
			a = i+1;
			wsUint N_med = 0;
#ifdef USE_SSE
			sse_t sse_Aia = sseSet(Ra_iCol[a]);
			N_med = getMed(i+1);
			for(b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_Bab = (sse_t *)(A[a] + b);

				sse_t *sse_Ci = (sse_t *)(C[i] + b);
				*sse_Ci = sseAdd(*sse_Ci, sseMul(sse_Aia, *sse_Bab));
			}
#endif
			for(b=N_med ; b<=i ; b++)
				C[i][b] += Ra_iCol[a] * A[a][b];
		}

		for(a=i+2 ; a<N_s1 ; a++) {
			wsUint N_med = 0;
#ifdef USE_SSE
			sse_t sse_Aia = sseSet(Ra_iCol[a]);
			N_med = getMed(i+1);
			for(b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_Bab = (sse_t *)(A[a] + b);

				sse_t *sse_Ci = (sse_t *)(C[i] + b);
				*sse_Ci = sseAdd(*sse_Ci, sseMul(sse_Aia, *sse_Bab));
			}
#endif
			for(b=N_med ; b<=i ; b++)
				C[i][b] += Ra_iCol[a] * A[a][b];
		}
	}

	sseFree(Ra_iCol);
	return C;
#else
	return sseMpMt(A, N_s1, N_s1, B, N_s2, N_s2);
#endif
}

wsMat sseMpS(wsMat Ra_mat, wsUint N_row, wsUint N_col, wsSym Ra_s, wsUint N_sz)
{
	wsUint i;
	if (N_col != N_sz)
		halt("Column size of first matrix [%d] is not match with dimension "
		"of second symmetric matrix [%d]", N_col, N_sz);
	wsMat	Ra_ret = sseMatrix(N_row, N_sz);
#ifdef USE_SYM
	wsUint a, b;

	for (i=0; i < N_row; i++) {
		memset(Ra_ret[i], 0x00, sizeof(wsReal)*N_sz);
		for(a = 0; a < N_sz; a++) {
			wsReal *Ra_m = Ra_mat[i];
			register wsReal R_a = Ra_m[a];
			wsReal *Ra_x = Ra_s[a];
			wsReal *Ra_r = Ra_ret[i];

#ifdef USE_SSE
			wsUint	N_med = getMed(a);
			sse_t	sse_a = sseSet(R_a);
			sse_t	sse_ra = sseSet(0.0);

			for (b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_r = (sse_t *)(Ra_r + b);
				sse_t *sse_s = (sse_t *)(Ra_x + b);
				sse_t *sse_m = (sse_t *)(Ra_m + b);
				*sse_r = sseAdd(*sse_r, sseMul(sse_a, *sse_s));
				sse_ra = sseAdd(sse_ra, sseMul(*sse_m, *sse_s));
			}
			sseAppendSum(sse_ra, Ra_r[a]);
			for(b=N_med ; b < a; b++) {
				Ra_r[a] += Ra_m[b] * Ra_x[b];
				Ra_r[b] += R_a * Ra_x[b];
			}
#else
			for(b = 0; b < a; b++) {
				Ra_r[a] += Ra_m[b] * Ra_x[b];
				Ra_r[b] += R_a * Ra_x[b];
			}
#endif
			Ra_r[a] += R_a * Ra_s[a][a];
		}
	}

	return Ra_ret;
#else
	wsAlloc(Ra_ret, wsReal*, N_row);

	for (i=0 ; i<N_row ; i++)
		Ra_ret[i] = sseVpS(Ra_mat[i], N_col, Ra_s, N_sz);
	return Ra_ret;
#endif
}

/* Perform ABA^T */
wsSym sseMSMt(wsReal **Ra_mat, wsUint N_row, wsUint N_col,
			  wsReal **Ra_symMat, wsUint N_sz)
{
#ifdef USE_SYM
	wsSym	Ra_ret = sseSymMat(N_row);
	wsUint	i;

	for (i=0 ; i<N_row ; i++) {
		/* Calculate A[i,]*B */
		wsReal *Ra_interm = sseVpS(Ra_mat[i], N_col, Ra_symMat, N_sz);

		for (wsUint j=0 ; j<=i ; j++)
			Ra_ret[i][j] = sseVV(Ra_interm, N_sz, Ra_mat[j]);
		sseFree(Ra_interm);
	}

	return Ra_ret;
#else
	return sseMMMt(Ra_mat, N_row, N_col, Ra_symMat, N_sz, N_sz);	
#endif
}

wsSym sseSSS(wsSym Ra_s1, wsUint N_s1, wsSym Ra_s2, wsUint N_s2)
{
	wsSym Ra_sMed = sseSpS(Ra_s1, N_s1, Ra_s2, N_s2);
	wsSym Ra_ret = sseSpS(Ra_sMed, N_s1, Ra_s1, N_s1);
	sseUnmat(Ra_sMed, N_s1);
	return Ra_ret;
}

/*
sm <- function(s, m) {
  if (dim(s)[2] != dim(m)[1])
    stop("Can't do sm")

  ret <- matrix(0, dim(s)[1], dim(m)[2])
  b <- 1:dim(s)[2]

  ret[1,b] <- ret[1,b]+s[1,1]*m[1,b]

  for (i in 2:dim(s)[1]) {
    for (a in 1:(i-1)) {
      ret[i,b] <- ret[i,b]+s[i,a]*m[a,b]
      ret[a,b] <- ret[a,b]+s[i,a]*m[i,b]
    }
    ret[i,b] <- ret[i,b]+s[i,i]*m[i,b]
  }
  ret
}
 */
wsReal** sseSpM(wsReal **Ra_symMat, wsUint N_sz, wsReal **Ra_mat,
	wsUint N_row, wsUint N_col)
{
#ifdef USE_SYM
	wsUint i, a, b;
	if (N_row != N_sz)
		halt("Row size of second matrix [%d] is not match with dimension of first symmetric matrix [%d]",
		N_row, N_sz);
	wsReal **Ra_ret = sseEmptyMatrix(N_sz, N_col);
	wsUint N_med;

	for (b=0 ; b<N_col ; b++)
		Ra_ret[0][b] += Ra_symMat[0][0] * Ra_mat[0][b];

	for (i=1 ; i<N_sz ; i++) {
		wsReal *Ra_si = Ra_symMat[i];
		wsReal *Ra_ri = Ra_ret[i];
		wsReal *Ra_mi = Ra_mat[i];

		for (a=0 ; a<i ; a++) {
			wsReal *Ra_ra = Ra_ret[a];
			wsReal *Ra_ma = Ra_mat[a];
			wsReal R_s_ia = Ra_si[a];
#ifdef USE_SSE
			N_med = getMed(N_col);
			sse_t sse_s_ia = sseSet(R_s_ia);

			for (b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_ri = (sse_t *)(Ra_ri + b);
				sse_t *sse_ma = (sse_t *)(Ra_ma + b);
				//Ra_ri[b] += R_as * Ra_ma[b];
				*sse_ri = sseAdd(*sse_ri, sseMul(sse_s_ia, *sse_ma));
				sse_t *sse_ra = (sse_t *)(Ra_ra + b);
				sse_t *sse_mi = (sse_t *)(Ra_mi + b);
				//Ra_ra[b] += R_as * Ra_mi[b];
				*sse_ra = sseAdd(*sse_ra, sseMul(sse_s_ia, *sse_mi));
			}
#else
			N_med = 0;
#endif
			for (b=N_med ; b<N_col ; b++) {
				Ra_ri[b] += R_s_ia * Ra_ma[b];
				Ra_ra[b] += R_s_ia * Ra_mi[b];
			}
		}
#ifdef USE_SSE
		N_med = getMed(N_col);
		sse_t sse_s_ii = sseSet(Ra_si[i]);
		for (b=0 ; b<N_med ; b+=sseJmp) {
			sse_t *sse_ri = (sse_t *)(Ra_ri + b);
			sse_t *sse_mi = (sse_t *)(Ra_mi + b);
			*sse_ri = sseAdd(*sse_ri, sseMul(sse_s_ii, *sse_mi));
			// Ra_ri[b] += Ra_si[i] * Ra_mi[b];
		}
#else
		N_med = 0;
#endif
		for (b=N_med ; b<N_col ; b++)
			Ra_ri[b] += Ra_si[i] * Ra_mi[b];
	}
	return Ra_ret;
#else
	return sseMpM(Ra_symMat, N_sz, N_sz, Ra_mat, N_row, N_col);
#endif
// #ifdef USE_SSE
// 	sse_t *sse_m = NULL;
// 	sse_t *sse_r = NULL;
// #endif
// 	wsReal *Ra_tmp = NULL;
// 	sseMalloc(Ra_tmp, wsReal, N_sz);
// 	for(int i = N_sz-1; i >= 0; i--) {
// 		for(int a = 0; a <=i ; a++) {
// 			wsReal R_ia = Ra_symMat[i][a];
// 			wsUint N_med = 0;
// // #ifdef USE_SSE
// // 			N_med = getMed(N_col);
// // 			sse_t sse_ia = sseSet(R_ia);
// // 			for(int b = 0; b < N_med; b+=sseJmp) {
// // 				sse_m = (sse_t *)(Ra_mat[a] + b);
// // 				sse_r = (sse_t *)(Ra_ret[i] + b);
// // 				*sse_r = sseAdd(*sse_r, sseMul(sse_ia, *sse_m));
// // //				Ra_ret[i][b] += R_ia * Ra_mat[a][b];
// // 			}
// // #endif
// 			for(int b = N_med; b < N_col; b++) {
// 				Ra_ret[i][b] += R_ia * Ra_mat[a][b];
// 			}
// 		}
// 		for(int a = i+1; a < N_sz; a++)
// 			Ra_tmp[a] = Ra_symMat[a][i];
// 
// 		for(int a = i+1; a < N_sz; a++) {
// 			wsReal R_ai = Ra_tmp[a];
// 			wsUint N_med = 0;
// // #ifdef USE_SSE
// // 			N_med = getMed(N_col);
// // 			sse_t sse_ia = sseSet(R_ai);
// // 			for(int b = 0; b < N_med; b+=sseJmp) {
// // 				sse_m = (sse_t *)(Ra_mat[a] + b);
// // 				sse_r = (sse_t *)(Ra_ret[i] + b);
// // 				*sse_r = sseAdd(*sse_r, sseMul(sse_ia, *sse_m));
// // 				//				Ra_ret[i][b] += R_ia * Ra_mat[a][b];
// // 			}
// // #endif
// 			for(int b = N_med; b < N_col; b++) {
// //			for(int b = 0; b < N_col; b++) {
// 				Ra_ret[i][b] += R_ai * Ra_mat[a][b];
// 			}
// 		}
// 	}
}

wsSym sseSpMt(wsSym Ra_m1, wsUint N_sz, wsReal **Ra_m2, wsUint N_r2,
	wsUint N_c2)
{
#ifdef USE_SYM
	wsUint i, j;
	/* Check sanity */
	if (N_sz != N_c2)
		halt_fmt(WISARD_SYST_INVL_DIM, "matrix multiplication", N_sz, N_c2);
	/* Make return matrix */
	wsReal **Ra_ret = sseMatrix(N_sz, N_r2);
	/* Make temporary buffer */
	wsReal *Ra_iCol = NULL;
	sseMalloc(Ra_iCol, wsReal, N_sz);

	for (i=0 ; i<N_sz ; i++) {
		/* Fill the ith column to array */
		for (j=i ; j<N_sz ; j++)
			Ra_iCol[j] = Ra_m1[j][i];

//		wsUint L=0;
		for (j=0; j < N_r2 ; j++)  {
			wsUint N_med = 0;
			wsReal R_sum = W0;

			/* Horizontal part */
#ifdef USE_SSE
			sse_t sse_s = sseSet(0.0);
			N_med = getMed(i+1);
			for (wsUint a=0 ; a<N_med ; a+=sseJmp) {
				sse_t *sse_m1 = (sse_t *)(Ra_m1[i] + a);
				sse_t *sse_m2 = (sse_t *)(Ra_m2[j] + a);
				sse_s = sseAdd(sse_s, sseMul(*sse_m1, *sse_m2));
			}
			sseSum(sse_s, R_sum);
//			/* SANITY */ L = N_med;
#endif
//			if (L != N_med) halt("ERROR1");
			for(wsUint a=N_med ; a <=i ; a++)
				R_sum += Ra_m1[i][a] * Ra_m2[j][a];
//			/* SANITY */ L = i+1;

			/* Vertical part */
			N_med = i+1;
			if ((i+1) < N_sz) {
				wsUint N_s = getMed(i+2);
#ifdef USE_SSE
				N_med = getMed(N_sz);
				if (N_s > N_med) N_s = N_med;
				sse_s = sseSet(0.0);
	//			if (L != (i+1)) halt("ERROR2");
				for (wsUint a=i+1 ; a<N_s ; a++)
					R_sum += Ra_iCol[a] * Ra_m2[j][a];
	//			/* SANITY */ L = N_s;

	//			if (L != N_s) halt("ERROR3");
				for (wsUint a=N_s ; a<N_med ; a+=sseJmp) {
					sse_t *sse_m1 = (sse_t *)(Ra_iCol + a);
					sse_t *sse_m2 = (sse_t *)(Ra_m2[j] + a);
					sse_s = sseAdd(sse_s, sseMul(*sse_m1, *sse_m2));
				}
	//			L = N_med;
				sseAppendSum(sse_s, R_sum);
#endif
			}
//			if (L != N_med) halt("ERROR4");z
			for(wsUint a=N_med ; a<N_sz ; a++)
				R_sum += Ra_iCol[a] * Ra_m2[j][a];

			Ra_ret[i][j] = R_sum;
		}
	}
	/* Free temporary buffer */
	sseFree(Ra_iCol);

	return Ra_ret;
#else
	return sseMpMt(Ra_m1, N_sz, N_sz, Ra_m2, N_r2, N_c2);
#endif
}

wsSym sym_sseMS(wsMat Ra_mat, wsUint N_row, wsUint N_col, wsSym Ra_s, wsUint N_sz)
{
	wsUint i;
	if (N_row != N_col)
		halt("First matrix must be square, but dimension is [%d*%d]",
		N_row, N_col);
	if (N_col != N_sz)
		halt("Column size of first matrix [%d] is not match with dimension of second symmetric matrix [%d]",
		N_col, N_sz);
	wsSym	Ra_ret = sseSymMat(N_row);
#ifdef USE_SYM
	wsUint a, b;

	for (i=0 ; i<N_row ; i++) {
		wsReal *Ra_m = Ra_mat[i];
		wsReal *Ra_r = Ra_ret[i];
		memset(Ra_r, 0x00, sizeof(wsReal)*(i+1));

		for(a=0 ; a<=i ; a++) {
			register wsReal R_a = Ra_m[a];
			wsReal *Ra_x = Ra_s[a];
			wsUint	N_med	= 0;

#ifdef USE_SSE
			N_med = getMed(a);
			sse_t	sse_a = sseSet(R_a);
			sse_t	sse_ra = sseSet(0.0);

			for (b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_r = (sse_t *)(Ra_r + b);
				sse_t *sse_s = (sse_t *)(Ra_x + b);
				sse_t *sse_m = (sse_t *)(Ra_m + b);
				*sse_r = sseAdd(*sse_r, sseMul(sse_a, *sse_s));
				sse_ra = sseAdd(sse_ra, sseMul(*sse_m, *sse_s));
			}
			sseAppendSum(sse_ra, Ra_r[a]);
			for (b=N_med ; b<a ; b++) {
				Ra_r[a] += Ra_m[b] * Ra_x[b];
				Ra_r[b] += R_a * Ra_x[b];
			}
#endif
			for (b=N_med ; b<a ; b++) {
				Ra_r[a] += Ra_m[b] * Ra_x[b];
				Ra_r[b] += R_a * Ra_x[b];
			}
			Ra_r[a] += R_a * Ra_s[a][a];
		}
		for(a=i+1 ; a<N_sz ; a++) {
			register wsReal R_a = Ra_m[a];
			wsReal *Ra_x	= Ra_s[a];
			wsUint	N_med	= 0;

#ifdef USE_SSE
			N_med = getMed(i);
			sse_t	sse_a = sseSet(R_a);
			//sse_t	sse_ra = sseSet(0.0);

			for (b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_r = (sse_t *)(Ra_r + b);
				sse_t *sse_s = (sse_t *)(Ra_x + b);
				//sse_t *sse_m = (sse_t *)(Ra_m + b);
				*sse_r = sseAdd(*sse_r, sseMul(sse_a, *sse_s));
				//sse_ra = sseAdd(sse_ra, sseMul(*sse_m, *sse_s));
			}
			//sseAppendSum(sse_ra, Ra_r[a]);
			for(b=N_med ; b<i ; b++)
				Ra_r[b] += R_a * Ra_x[b];
#endif
			for(b=N_med ; b<i ; b++)
				Ra_r[b] += R_a * Ra_x[b];
		}
	}

	return Ra_ret;
#else
	wsAlloc(Ra_ret, wsReal*, N_row);

	for (i=0 ; i<N_row ; i++)
		Ra_ret[i] = sseVpS(Ra_mat[i], N_col, Ra_s, N_sz);
	return Ra_ret;
#endif
}

wsVec sseSpV(wsSymCst Ra_symMat, wsUintCst N_sz, wsVecCst Ra_vec)
{
#ifdef USE_SYM
	wsVec Ra_ret = sseEmptyVec(N_sz);

	LOOP (i, N_sz) {
		wsReal R_v = Ra_vec[i];
		LOOP (j, i) {
			Ra_ret[j] += Ra_symMat[i][j] * R_v;
			Ra_ret[i] += Ra_symMat[i][j] * Ra_vec[j];
		}
		/* When i == j */
		Ra_ret[i] += Ra_symMat[i][i]*Ra_vec[i];
	}

	return Ra_ret;
#else
	return sseMpV(Ra_symMat, N_sz, N_sz, Ra_vec);
#endif
}

wsReal* sseVpS(wsReal *Ra_vec, wsUint N_sz, wsSym Ra_symMat,
	wsReal *Ra_mask/*=NULL*/)
{
	return sseVpS(Ra_vec, N_sz, Ra_symMat, N_sz, Ra_mask);
}

wsReal* sseVpS(wsReal *Ra_v, wsUint N_vLen, wsSym Ra_symMat, wsUint N_sz,
			   wsReal *Ra_mask/*=NULL*/, wsReal *Rp_r/*=NULL*/)
{
	wsUint a, b;
	wsReal *Ra_r = NULL;

	if (N_vLen != N_sz)
		halt("Vector size[%d] is not match to the size of symmetric matrix[%d]",
		N_vLen, N_sz);
	if (Rp_r) {
		Ra_r = Rp_r;
		memset(Ra_r, 0x00, sizeof(wsReal)*N_sz);
	} else
		sseCalloc(Ra_r, wsReal, N_sz);

	if (Ra_mask == NULL) {
		for (a=0 ; a<N_sz ; a++) {
			register wsReal R_a = Ra_v[a];
			wsReal *Ra_s = Ra_symMat[a];
			wsUint	N_med = 0;
#ifdef USE_SSE
			N_med = getMed(a);
			sse_t	sse_a = sseSet(R_a);
			sse_t	sse_ra = sseSet(0.0);

			for (b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_r = (sse_t *)(Ra_r + b);
				sse_t *sse_s = (sse_t *)(Ra_s + b);
				sse_t *sse_m = (sse_t *)(Ra_v + b);

				*sse_r = sseAdd(*sse_r, sseMul(sse_a, *sse_s));
				sse_ra = sseAdd(sse_ra, sseMul(*sse_m, *sse_s));
			}
			sseAppendSum(sse_ra, Ra_r[a]);
#endif
			for(b=N_med ; b<a ; b++) {
				Ra_r[a] += Ra_v[b] * Ra_s[b];
				Ra_r[b] += R_a * Ra_s[b];
			}
			/* When b == a */
			Ra_r[a] += R_a * Ra_s[a];
		}
	} else {
		for (a=0 ; a<N_sz ; a++) {
			if (Ra_mask[a] == W0) continue;

			register wsReal R_a = Ra_v[a];
			wsReal *Ra_s = Ra_symMat[a];
			wsUint	N_med = 0;
#ifdef USE_SSE
			N_med = getMed(a);
			sse_t	sse_a = sseSet(R_a);
			sse_t	sse_ra = sseSet(0.0);

			for (b=0 ; b<N_med ; b+=sseJmp) {
				sse_t *sse_r = (sse_t *)(Ra_r + b);
				sse_t *sse_s = (sse_t *)(Ra_s + b);
				sse_t *sse_m = (sse_t *)(Ra_v + b);
				sse_t *sse_M = (sse_t *)(Ra_mask + b);

				*sse_r = sseAdd(*sse_r, sseAnd(*sse_M, sseMul(sse_a, *sse_s)));
				sse_ra = sseAdd(sse_ra, sseAnd(*sse_M, sseMul(*sse_m, *sse_s)));
			}
			sseAppendSum(sse_ra, Ra_r[a]);
#endif
			for(b=N_med ; b<a ; b++) {
				if (Ra_mask[b]) {
					Ra_r[a] += Ra_v[b] * Ra_s[b];
					Ra_r[b] += R_a * Ra_s[b];
				}
			}
			/* When b == a */
			if (Ra_mask[a])
				Ra_r[a] += R_a * Ra_s[a];
		}
	}

	return Ra_r;
}

void sseVpS(wsReal *Ra_v, wsUint N_vLen, wsSym Ra_symMat, wsUint N_sz,
			wsUint N_calc, wsReal *Rp_r)
{
	wsUint a, b;

	if (N_vLen != N_sz)
		halt("Vector size[%d] is not match to the size of symmetric matrix[%d]",
		N_vLen, N_sz);
	memset(Rp_r, 0x00, sizeof(wsReal)*N_calc);

	for (a=0 ; a<N_sz ; a++) {
		register wsReal R_a = Ra_v[a];
		wsReal *Ra_s = Ra_symMat[a];
		wsUint	N_med	= 0;
		wsUint	N_to	= a>N_calc ? N_calc : a;
#ifdef USE_SSE
		N_med = getMed(N_to);
		sse_t	sse_a = sseSet(R_a);
		sse_t	sse_ra = sseSet(0.0);

		for (b=0 ; b<N_med ; b+=sseJmp) {
			sse_t *sse_r = (sse_t *)(Rp_r + b);
			sse_t *sse_s = (sse_t *)(Ra_s + b);
			sse_t *sse_m = (sse_t *)(Ra_v + b);

			*sse_r = sseAdd(*sse_r, sseMul(sse_a, *sse_s));
			sse_ra = sseAdd(sse_ra, sseMul(*sse_m, *sse_s));
		}
		sseAppendSum(sse_ra, Rp_r[a]);
#endif
		for(b=N_med ; b<N_to ; b++) {
			Rp_r[a] += Ra_v[b] * Ra_s[b];
			Rp_r[b] += R_a * Ra_s[b];
		}
		/* When b == a */
		if (a < N_calc)
			Rp_r[a] += R_a * Ra_s[a];
	}
}

bool symSVDinverse(wsReal **Ra_matOrig, wsUint N_sz, wsSym *Rp_res,
				wsUint *Np_rank=NULL)
{
	wsUint			i;//, j, k;
//	const wsReal	R_eps		= numeric_limits<wsReal>::epsilon();
	const wsReal	R_eps		= 3.66685286250103598289229456241855586995370686054229736328125e-11;
	wsReal			**Ra_v		= sseEmptyMatrix(N_sz, N_sz); /* checked */
	wsReal			*Ra_wgt		= NULL; /* checked */
	sseCalloc(Ra_wgt, wsReal, N_sz);
	wsReal			**Ra_mat	= sseMatrixP(N_sz, N_sz, Ra_matOrig); /* checked */
	bool			B_cvged		= SingularValueDecomposition(Ra_mat, N_sz, Ra_wgt, Ra_v);

#ifdef _DEBUG
// 	exportMatrix("m1", Ra_mat, N_sz, N_sz);
// 	exportVector("v2", Ra_wgt, N_sz);
// 	exportMatrix("m3", Ra_v, N_sz, N_sz);
#endif

	// Look for singular values
	wsReal wmax = W0;
	for (i=0; i<N_sz; i++)
		wmax = Ra_wgt[i] > wmax ? Ra_wgt[i] : wmax;
	wsReal wmin = wmax * R_eps;
	for (i=0; i<N_sz; i++)
		Ra_wgt[i] = Ra_wgt[i] < wmin ? 0 : 1/Ra_wgt[i];

	// Get rank of Ra_mat if assigned
	if (Np_rank) {
		*Np_rank = 0;
		for (i=0; i<N_sz; i++)
			if (Ra_wgt[i]) (*Np_rank)++;
	}

	// u w t(v)
	// row U * 1/w
	// results matrix
	for (i=0 ; i<N_sz ; i++)
		sseVpV(Ra_mat[i], Ra_wgt, Ra_mat[i], N_sz);
	sseFree(Ra_wgt);
// 	for (i=0; i<N_sz; i++)
// 		for (j=0; j<N_sz; j++)
// 			Ra_mat[i][j] *= Ra_wgt[j];

	// [nxn].[t(v)] 
	//	printf("O matrix[%x]\n", O);
	wsSym Ra_ret = sym_sseMpMt(Ra_mat, N_sz, N_sz, Ra_v, N_sz, N_sz);
// 	for (i=0 ; i<N_sz ; i++)
// 		for (j=0 ; j<N_sz ; j++)
// 			for (k=0 ; k<N_sz ; k++)
// 				Ra_ret[i][j] += Ra_mat[i][k] * Ra_v[j][k];

	sseUnmat(Ra_mat, N_sz);
	sseUnmat(Ra_v, N_sz);
	*Rp_res = Ra_ret;
	return B_cvged;
}

bool symSVDinverse(wsReal **Ra_matOrig, wsUint N_sz, wsSym *Rp_res,
				   wsReal *Rp_logDet=NULL)
{
	wsUint			i;//, j, k;
//	const wsReal	R_eps		= numeric_limits<wsReal>::epsilon();
	const wsReal	R_eps		= 3.66685286250103598289229456241855586995370686054229736328125e-11;
	wsReal			**Ra_v		= sseEmptyMatrix(N_sz, N_sz); /* checked */
	wsReal			*Ra_wgt		= NULL; /* checked */
	sseCalloc(Ra_wgt, wsReal, N_sz);
	wsReal			**Ra_mat	= sseMatrixP(N_sz, N_sz, Ra_matOrig); /* checked */
	bool			B_cvged		= SingularValueDecomposition(Ra_mat, N_sz, Ra_wgt, Ra_v);

#ifdef _DEBUG
// 	exportMatrix("m1", Ra_mat, N_sz, N_sz);
// 	exportVector("v2", Ra_wgt, N_sz);
// 	exportMatrix("m3", Ra_v, N_sz, N_sz);
#endif

	// Look for singular values
	wsReal wmax = W0;
	for (i=0; i<N_sz; i++)
		wmax = Ra_wgt[i] > wmax ? Ra_wgt[i] : wmax;
	wsReal wmin = wmax * R_eps;
	for (i=0; i<N_sz; i++)
		Ra_wgt[i] = Ra_wgt[i] < wmin ? 0 : 1/Ra_wgt[i];

	// Get log(determinant) of Ra_mat if assigned
	if (Rp_logDet) {
		*Rp_logDet = 0;
		for (i=0; i<N_sz; i++)
			if (Ra_wgt[i] > W0) *Rp_logDet += log(Ra_wgt[i]);
	}

	// u w t(v)
	// row U * 1/w
	// results matrix
	for (i=0 ; i<N_sz ; i++)
		sseVpV(Ra_mat[i], Ra_wgt, Ra_mat[i], N_sz);
	sseFree(Ra_wgt);
	// 	for (i=0; i<N_sz; i++)
	// 		for (j=0; j<N_sz; j++)
	// 			Ra_mat[i][j] *= Ra_wgt[j];

	// [nxn].[t(v)] 
	//	printf("O matrix[%x]\n", O);
	wsSym Ra_ret = sym_sseMpMt(Ra_mat, N_sz, N_sz, Ra_v, N_sz, N_sz);
	// 	for (i=0 ; i<N_sz ; i++)
	// 		for (j=0 ; j<N_sz ; j++)
	// 			for (k=0 ; k<N_sz ; k++)
	// 				Ra_ret[i][j] += Ra_mat[i][k] * Ra_v[j][k];

	sseUnmat(Ra_mat, N_sz);
	sseUnmat(Ra_v, N_sz);
	*Rp_res = Ra_ret;
	return B_cvged;
}

cSymMatrix::cSymMatrix()
{
	X_type = MATTYPE_SYM;
	__init(NULL, 0, 0, 1);
}

cSymMatrix::cSymMatrix(cSymMatrix &&S)
{
	X_type = MATTYPE_SYM;
	B_ddealloc = 0;
	Ra_data = S.get();
	N_row = S.row();
	N_col = S.col();

	S.setDontDealloc();
}

cSymMatrix::cSymMatrix(wsUintCst N_sz, char B_rand/*=0*/)
{
	Ra_data = NULL;
	init(NULL, WISARD_NAN, N_sz, B_rand, 0);
}

cSymMatrix::cSymMatrix(wsSym Ra_mat, wsUintCst N_sz, char B_inpDDealloc/*=0*/)
{
	Ra_data = NULL;
	init(Ra_mat, WISARD_NAN, N_sz, 0, B_inpDDealloc);
}

cSymMatrix::cSymMatrix(wsRealCst R_val, wsUintCst N_sz)
{
	Ra_data = NULL;
	init(NULL, R_val, N_sz);
}

cSymMatrix::cSymMatrix(wsRealCst R_Vdiag, wsRealCst R_VoffDiag, wsUintCst N_sz)
{
	Ra_data = NULL;
	init(NULL, R_VoffDiag, N_sz);
	LOOP (i, N_sz) Ra_data[i][i] = R_Vdiag;
}

cSymMatrix::cSymMatrix(wsVec Ra_mat, wsUintCst N_sz)
{
	Ra_data = NULL;
	init(Ra_mat, N_sz);
}

cSymMatrix::~cSymMatrix()
{
//	LOG("~cSymMatrix() [%x] called\n", this);
	rem();
}

void cSymMatrix::init(wsUint N_sz)
{
	init(NULL, WISARD_NAN, N_sz);
}

void cSymMatrix::init(wsSym Ra_mat, wsUint N_sz)
{
	__init(Ra_mat, N_sz, N_sz);
	B_ddealloc = 0;
}

void cSymMatrix::init(wsVec Ra_mat, wsUint N_sz)
{
	/* Covert matrix into two-dim */
	wsSym Ra_dst = sseSymMat(N_sz);
	wsUint N_idx = 0;
	LOOP (i, N_sz) {
		memcpy(Ra_dst[i], Ra_mat+N_idx, sizeof(wsReal)*(i+1));
		N_idx += i+1;
	}
	__init(Ra_dst, N_sz, N_sz);
	B_ddealloc = 0;
}


void cSymMatrix::init(wsSym Ra_mat, wsReal R_val, wsUint N_sz, char B_rand/*=0*/,
	char B_inpDDealloc/*=0*/)
{
	if (Ra_data && !B_inpDDealloc) sseUnmat(Ra_data, N_row);

	B_ddealloc	= B_inpDDealloc;
	X_type		= MATTYPE_SYM;

	/* Allocate matrix */
	if (Ra_mat)
		__init(Ra_mat, N_sz, N_sz);
//		Ra_data = Ra_mat;
	else if (!NA(R_val)) {
		__init(sseSymMat(N_sz), N_sz, N_sz);
		sseSinit(Ra_data, N_sz, R_val);
	} else
		__init(sseEmptySymMat(N_sz), N_sz, N_sz);

//		Ra_data = sseEmptySymMat(N_sz);
	/* Definitely square matrix */
//	N_row = N_col = N_sz;
	if (B_rand)
		for (wsUint i=0 ; i<row() ; i++)
			for (wsUint j=0 ; j<=i ; j++)
				Ra_data[i][j] = (wsReal)(rand()%1000) - REAL_CONST(500.0);
}

void cSymMatrix::init(cSymMatrix &S)
{
	/* Do hard-copy given symmetric matrix */
	wsSym Ra_buf = sseSymMatP(S.row(), S.get());
	init(Ra_buf, WISARD_NAN, S.row());
}

cStdMatrix cSymMatrix::operator*(const cStdMatrix& M)
{
	/* col(this) == S.row() */
	if (col()!=M.row()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseSpM(get(), row(), M.get(), M.row(), M.col());
	return cStdMatrix(row(), M.col(), Ra_ret);
}

cStdMatrix cSymMatrix::operator*(const cSymMatrix& S)
{
	/* col(this) == S.row() */
	if (col()!=S.row()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseSpS(get(), row(), S.get(), row());
	return cStdMatrix(row(), S.col(), Ra_ret);
}

cStdMatrix cSymMatrix::operator*(const cDiagMatrix& D)
{
	/* col(this) == S.row() */
	if (col()!=D.row()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseSpD(get(), row(), D.get()[0]);
	return cStdMatrix(row(), D.col(), Ra_ret);
}

cStdMatrix cSymMatrix::operator*(const cBlkMatrix& B)
{
	/* col(this) == S.row() */
	if (col()!=B.row()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseSpB(get(), row(), B.get(), B.blk(), B.rep());
	return cStdMatrix(row(), B.col(), Ra_ret);
}

cSymMatrix& cSymMatrix::operator*=(wsReal R)
{
	wsUint N_r = row();
	/* Calc */
	for (wsUint i=0 ; i<N_r ; i++)
		sseVpC(get()[i], R, get()[i], i+1);
	/* Return */
	return *this;
}

cSymMatrix& cSymMatrix::operator+=(wsReal R)
{
	wsUint N_r = row();
	/* Calc */
	for (wsUint i=0 ; i<N_r ; i++)
		sseVaC(get()[i], R, get()[i], i+1);
	/* Return */
	return *this;
}

cSymMatrix	cSymMatrix::operator+(const cSymMatrix& S)
{
	/* col(this) == S.row() */
	if (col()!=S.row()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_R = sseSymMat(col());
	/* Calc */
	for (wsUint i=0,N=col() ; i<N ; i++)
		sseVaV(get()[i], S.get()[i], Ra_R[i], i+1);
	/* Return */
	return cSymMatrix(Ra_R, col());
}

cSymMatrix	cSymMatrix::operator+(const cDiagMatrix& D)
{
	/* col(this) == S.row() */
	if (col()!=D.row()) raise(MATERR_DIM_INCORRECT);
	/* Prepare */
	wsReal **Ra_R = sseSymMatP(col(), get());
	wsReal *Dp = D.get()[0];
	/* Calc */
	for (wsUint i=0,N=col() ; i<N ; i++) Ra_R[i][i] += Dp[i];
	/* Return */
	return cSymMatrix(Ra_R, col());
}

cSymMatrix& cSymMatrix::operator-=(cSymMatrix &S)
{
	if (row()!=S.row()) raise(MATERR_DIM_INCORRECT);
	/* Prepare */
	wsReal **Sp	= S.get();
	wsReal **T	= get();
	for (wsUint i=0,N=row() ; i<N ; i++)
		sseVsV(T[i], Sp[i], T[i], i+1);
	return *this;
}

cSymMatrix& cSymMatrix::operator+=(cSymMatrix &S)
{
	if (row()!=S.row()) raise(MATERR_DIM_INCORRECT);
	/* Prepare */
	wsReal **Sp	= S.get();
	wsReal **T	= get();
	for (wsUint i=0,N=row() ; i<N ; i++)
		sseVaV(T[i], Sp[i], T[i], i+1);
	return *this;
}

cSymMatrix& cSymMatrix::operator+=(cIdtMatrix &I)
{
	if (row()!=I.row()) raise(MATERR_DIM_INCORRECT);
	/* Prepare */
	wsReal **T	= get();
	for (wsUint i=0,N=row() ; i<N ; i++)
		T[i][i] += W1;
	return *this;
}

cSymMatrix& cSymMatrix::operator-=(wsReal R_val)
{
	/* Prepare */
	wsReal **T	= get();
	for (wsUint i=0,N=row() ; i<N ; i++)
		sseVsC(T[i], R_val, T[i], i+1);
	return *this;
}

cSymMatrix cSymMatrix::operator-(const cSymMatrix &S)
{
	if (row()!=S.row()) raise(MATERR_DIM_INCORRECT);
	/* Prepare */
	wsReal **Sp	= S.get();
	wsReal **R	= sseSymMat(row());
	wsReal **T	= get();
	for (wsUint i=0,N=row() ; i<N ; i++)
		sseVsV(T[i], Sp[i], R[i], i+1);
	return cSymMatrix(R, row());
}

cSymMatrix cSymMatrix::operator-(const cIdtMatrix &I)
{
	if (row()!=I.row()) raise(MATERR_DIM_INCORRECT);
	/* Prepare */
	wsReal **R	= sseSymMatP(row(), get());
	for (wsUint i=0,N=row() ; i<N ; i++)
		R[i][i] -= W1;
	return cSymMatrix(R, row());
}

cStdMatrix cSymMatrix::operator-(const cStdMatrix &M)
{
	if (row()!=M.row() || !M.sqr()) raise(MATERR_DIM_INCORRECT);
	/* Prepare */
	wsUint N = row();
	wsReal **S = get();
	wsReal **X = M.get();
	wsReal **R = sseMatrix(N, N);
	for (wsUint i=0 ; i<N ; i++) {
		wsReal *Si = S[i];
		/* Fill lower-triangular part */
		sseVsV(Si, X[i], R[i], i+1);
		/* Fill upper-triangular part */
		for (wsUint j=0 ; j<i ; j++)
			R[j][i] = Si[j] - X[j][i];
	}
	return cStdMatrix(N, N, R);
}

cSymMatrix& cSymMatrix::operator/=(wsRealCst R_val)
{
	sseSpC((wsSymCst)get(), W1/R_val, get(), row());
	return *this;
}

/*
 *
 * Special operators definition
 *
 */

wsRealCst cSymMatrix::det()
{
	/* FIXME : Need implementation */
	return W0;
}

wsRealCst cSymMatrix::detL()
{
	/* FIXME : Need implementation */
	return W0;
}

wsRealCst cSymMatrix::normF(cStdMatrix &M)
{
	return compSMfrobenius(get(), row(), M.get(), M.row(), M.col());
}

wsRealCst cSymMatrix::normF(cSymMatrix &S)
{
	return compSSfrobenius(get(), row(), S.get(), S.row());
}

wsRealCst cSymMatrix::normF()
{
	return normFrobenius(get(), row());
}

cStdMatrix cSymMatrix::mt(const cStdMatrix &M)
{
	if (col()!=M.col()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret = sseSpMt(get(), row(), M.get(), M.row(), M.col());
	return cStdMatrix(row(), M.row(), Ra_ret);
}

cSymMatrix&	cSymMatrix::operator=(cSymMatrix &S)
{
	init(S.get(), WISARD_NAN, S.row());
	return *this;
}

cSymMatrix& cSymMatrix::operator=(cSymMatrix &&S)
{
	B_ddealloc = 0;
	Ra_data = S.get();
	N_row = S.row();
	N_col = S.col();

	S.setDontDealloc();
	return *this;
}

cVector cSymMatrix::diagInv2vec()
{
	wsMat	M	= get();

	wsReal	*pR	= sseVector(row());
	for (wsUint i=0,N=row() ; i<N ; i++)
		pR[i] = W1 / M[i][i];

	return cVector(pR, row());
}

wsRealCst cSymMatrix::tr2()
{
	return sseTrSS(get(), row());
}

wsReal cSymMatrix::tr(cSymMatrix &S)
{
	if (row() != S.row()) raise(MATERR_DIM_INCORRECT);
	return sseTrSS(get(), row(), S.get());
}

/*
 *
 * Predefined operators definition
 *
 */

wsRealCst cSymMatrix::sum()
{
	return sseSsum(get(), row());
}

bool cSymMatrix::sym()
{
	/* Since cSymMatrix is ALWAYS symmetric, return true */
	return true;
}

/*
mychol <- function(M) {
L <- dim(M)[1]

p <- rep(0, L)
S <- M[1,1]
X <- M
if (S <= 0) stop("Chol failed\n")
p[1] <- sqrt(S)

M[1:L,1] <- M[1:L,1] / p[1]

if(L>=2) for (i in 2:L) {
R <- 1:(i-1)
S <- M[i,i] - sum(M[i,R]^2)
if (S <= 0) stop("Chol failed (", i, ")\n")
p[i] <- sqrt(S)

if ((i+1)<=L) for (j in (i+1):L) {
S <- X[i,j] - sum(M[i,R] * M[j,R])
M[j,i] <- S/p[i]
}
}
diag(M) <- p

return(M * (lower.tri(M)+diag(L)))
}
}*/

/*
LtL <- function(L) {
  d <- dim(L)[1]
  R <- matrix(0, nrow=d, ncol=d)
  for (i in 1:d) {
    for (j in 1:i) {
      K <- 1:j
      R[j,K] <- R[j,K] + L[i,j]*L[i*K]
    }
  }
  R
}
*/
cSymMatrix& cSymMatrix::inv()
{
	return inv(NULL);
}

cSymMatrix& cSymMatrix::inv(char &B_isGinv)
{
	return inv(NULL, &B_isGinv);
}

cSymMatrix& cSymMatrix::inv(wsReal **Ra_eVec, wsReal *Ra_eVal, wsUint N_sz)
{
	const wsReal	R_eps		= numeric_limits<wsReal>::epsilon();
// 	for (wsUint i=0 ; i<N_sz ; i++)
// 		if (Ra_eVal[i] != W0)
// 			Ra_eVal[i] = W1 / Ra_eVal[i];

	// Look for singular values

	wsReal wmax = W0;
	for (wsUint i=0; i<N_sz; i++)
		wmax = Ra_eVal[i] > wmax ? Ra_eVal[i] : wmax;
	wsReal wmin = wmax * R_eps;
	for (wsUint i=0; i<N_sz; i++)
		Ra_eVal[i] = Ra_eVal[i] < wmin ? 0 : 1/Ra_eVal[i];

	//sseVinv(Ra_eI, N_sz);

	/* Calc inverse */
	wsSym	Ra_ret	= sseMDMt(Ra_eVec, N_sz, N_sz, Ra_eVal);
	return *(new cSymMatrix(Ra_ret, N_sz));

}

//#define OPTZ

cSymMatrix& cSymMatrix::inv(wsReal *Rp_logDet, char *Bp_isGinv/*=NULL*/,
	char B_halt/*=1*/)
{
	if (Bp_isGinv) *Bp_isGinv = 0;
	wsUint Z = row(), i, j, k;
	/* Allocate matrix for L */
	wsReal **R = sseEmptySymMat(Z);
	wsReal **Q = get(), **X;
	wsReal **L = NULL;
	cVector D(Z);
	wsReal *Dp = D.get();
	wsReal R_logDet = W0;

	if (Z == 1) {
		if (Rp_logDet)
			*Rp_logDet = log(fabs(Q[0][0]));
		if (Q[0][0] == W0) {
			if (B_halt)
				halt("Matrix is irreversible");
			else return *(new cSymMatrix());
		}
		R[0][0] = W1 / Q[0][0];

		return *(new cSymMatrix(R, 1));
	} else if (Z == 2) {
		wsReal R_det = Q[0][0]*Q[1][1] - SQR(Q[1][0]);
		if (Rp_logDet)
			*Rp_logDet = log(fabs(R_det));
		if (R_det == W0)
			goto _ginv;
		R[0][0] = Q[1][1] / R_det;
		R[1][1] = Q[0][0] / R_det;
		R[1][0] = -Q[1][0] / R_det;

		return *(new cSymMatrix(R, 2));
	}
	L = sseSymMatP(Z, Q);

	/* Do Cholesky decomposition against lower triangle */
	if (L[0][0] <= W0) {
		sseUnmat(L, Z);
		sseUnmat(R, Z);
//		LOG("Cholesky decomposition failed at 0, "
//			"may be the matrix is not positive definite\n");
		goto _ginv;
	}
	Dp[0] = ::sqrt(L[0][0]);
	for (i=0 ; i<Z ; i++) L[i][0] /= Dp[0];
	R[0][0] = W1 / Dp[0];

	for (i=1 ; i<Z ; i++) {
		Dp[i] = L[i][i] - sseVV(L[i], i);
		if (Dp[i] <= W0) {
			sseUnmat(L, Z);
			sseUnmat(R, Z);
//			LOG("Cholesky decomposition failed at %d, "
//				"may be the matrix is not positive definite\n", i);
			goto _ginv;
		}
		R_logDet += log(fabs(Dp[i]));
		Dp[i] = ::sqrt(Dp[i]);
//		printf("Dp[%d] %g\n", i, Dp[i]);

//S <- X[i,j] - sum(M[i,R] * M[j,R])
//M[j,i] <- S/p[i]

		for (j=i+1 ; j<Z ; j++)
			L[j][i] = (Q[j][i] - sseVV(L[i], i, L[j])) / Dp[i];

#ifdef OPTZ
		wsReal R_max = W0;
#endif
		for (j=0 ; j<i ; j++) {
			wsReal R_s = W0;
			for (k=j ; k<i ; k++)
				R_s -= L[i][k]*R[k][j];
			R[i][j] = R_s / Dp[i];
#ifdef OPTZ
			R_max = max(fabs(R[i][j]), R_max);
#endif
		}
		R[i][i] = W1 / Dp[i];
	}
	sseUnmat(L, Z);
//	for (i=0 ; i<Z ; i++) L[i][i] = Dp[i];
#ifdef OPTZ
	if (R_logDet < REAL_CONST(-20.0))
		goto _ginv;
#endif
	if (Rp_logDet)
		*Rp_logDet = R_logDet;

	/* Do t(L) %*% L */
	X = sseStS(R, Z);
	sseUnmat(R, Z);

	return *(new cSymMatrix(X, Z));
_ginv:
	if (Bp_isGinv) *Bp_isGinv = 1;
#ifdef _DEBUG
//	LOG("GINV applied\n");
#endif
	string s = "1";
	OPTION().assign("ginv", s.c_str(), 1);
	OPTION().FORCE_OPT_NUMBER(ginv);

	wsSym Ra_ret = NULL;
	wsReal **Ra_mat = sseSym2Mat(Q, Z);
#ifdef _DEBUG
//	exportMatrix("m", Ra_mat, Z, Z);
#endif
	bool ret = symSVDinverse(Ra_mat, N_row, &Ra_ret, Rp_logDet);
// 	if (Rp_logDet)
// 		halt("Failed to get determinant");
	sseUnmat(Ra_mat, Z);
	if (ret == false) {
		if (B_halt)
			halt("Generalized inverse with SVD failed");
		else
			return *(new cSymMatrix());
//		halt("Generalized inverse with SVD failed");
	}
	return *(new cSymMatrix(Ra_ret, Z));	
}

cSymMatrix& cSymMatrix::invOnly()
{
	wsUint Z = row(), i, j, k;
	/* Allocate matrix for L */
	wsReal **R = sseEmptySymMat(Z);
	wsReal **Q = get(), **X;
	wsReal **L = NULL;
	cVector D(Z);
	wsReal *Dp = D.get();
	wsReal R_logDet = W0;

	if (Z == 1) {
		if (Q[0][0] == W0)
			goto _ginv;
		R[0][0] = W1 / Q[0][0];

		return *(new cSymMatrix(R, 1));
	} else if (Z == 2) {
		wsReal R_det = Q[0][0]*Q[1][1] - SQR(Q[1][0]);
		if (R_det == W0)
			goto _ginv;
		R[0][0] = Q[1][1] / R_det;
		R[1][1] = Q[0][0] / R_det;
		R[1][0] = -Q[1][0] / R_det;

		return *(new cSymMatrix(R, 2));
	}
	L = sseSymMatP(Z, Q);

	/* Do Cholesky decomposition against lower triangle */
	if (L[0][0] <= W0) {
		sseUnmat(L, Z);
		sseUnmat(R, Z);
		//		LOG("Cholesky decomposition failed at 0, "
		//			"may be the matrix is not positive definite\n");
		goto _ginv;
	}
	Dp[0] = ::sqrt(L[0][0]);
	for (i=0 ; i<Z ; i++) L[i][0] /= Dp[0];
	R[0][0] = W1 / Dp[0];

	for (i=1 ; i<Z ; i++) {
		Dp[i] = L[i][i] - sseVV(L[i], i);
		if (Dp[i] <= W0) {
			sseUnmat(L, Z);
			sseUnmat(R, Z);
			//			LOG("Cholesky decomposition failed at %d, "
			//				"may be the matrix is not positive definite\n", i);
			goto _ginv;
		}
		R_logDet += log(fabs(Dp[i]));
		Dp[i] = ::sqrt(Dp[i]);
		//		printf("Dp[%d] %g\n", i, Dp[i]);

		//S <- X[i,j] - sum(M[i,R] * M[j,R])
		//M[j,i] <- S/p[i]

		for (j=i+1 ; j<Z ; j++)
			L[j][i] = (Q[j][i] - sseVV(L[i], i, L[j])) / Dp[i];

		for (j=0 ; j<i ; j++) {
			wsReal R_s = W0;
			for (k=j ; k<i ; k++)
				R_s -= L[i][k]*R[k][j];
			R[i][j] = R_s / Dp[i];
		}
		R[i][i] = W1 / Dp[i];
	}
	sseUnmat(L, Z);

	/* Do t(L) %*% L */
	X = sseStS(R, Z);
	sseUnmat(R, Z);

	return *(new cSymMatrix(X, Z));
_ginv:
	return *(new cSymMatrix());
}

cSymMatrix& cSymMatrix::ginv(wsReal *Rp_logDet/*=NULL*/)
{
	wsUint	Z		= row();
	/* Allocate matrix for L */
	wsMat	Q		= get();
	wsSym	Ra_ret	= NULL;
	wsMat	Ra_mat	= sseSym2Mat(Q, Z);
	bool	ret		= symSVDinverse(Ra_mat, N_row, &Ra_ret, Rp_logDet);
	//	if (Rp_logDet)
	//		halt("Failed to get determinant");
	sseUnmat(Ra_mat, Z);

	if (ret == false) {
		return *(new cSymMatrix());
		//		halt("Generalized inverse with SVD failed");
	}
	return *(new cSymMatrix(Ra_ret, Z));
}

cSymMatrix& cSymMatrix::subset(char *Ba_YorN, wsUint N_Yval)
{
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
	wsSym	Ra_ret	= sseSymMat(N_nLen);
	wsSym	Ra_src	= get();

	/* Make subset */
	if (N_Yval == 1) {
		/* 1 then include, 0 exclude */

		/* Get length */
		for (wsUint i=0,I=0,N=row() ; i<N ; i++) {
			if (!Ba_YorN[i]) continue;

			for (wsUint j=0,J=0 ; j<=i ; j++) {
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

			for (wsUint j=0,J=0 ; j<=i ; j++) {
				if (Ba_YorN[j]) continue;

				Ra_ret[I][J] = Ra_src[i][j];
				J++;
			}
			I++;
		}
	}

	return *(new cSymMatrix(Ra_ret, N_nLen));
}

cSymMatrix& cSymMatrix::krnk(const cSymMatrix &S)
{
	wsSym Ra_ret = krnkSS(get(), row(), S.get(), S.row());

	return *(new cSymMatrix(Ra_ret, row()*S.row()));
}

cSymMatrix& cSymMatrix::krnk(const cIdtMatrix &I)
{
	wsSym Ra_ret = krnkSI(get(), row(), I.row());

	return *(new cSymMatrix(Ra_ret, row()*I.row()));
}

cSymMatrix cSymMatrix::einv(cSymMatrix *Mp_sqInv/*=NULL*/)
{
	char	B_fail	= 0;
	wsUint	N_sz	= row();
	wsMat	Ra_full = sseSym2Mat(get(), N_sz);
	wsMat	Ra_E	= NULL; /* P */
	wsReal*	Ra_eI	= EIGENDECOMPOSITION(Ra_full, N_sz, &Ra_E, 1);
	if (B_fail) halt("Eigendecomposition failed");
	wsReal*	Ra_esqI	= sseVector(N_sz);
	for (wsUint i=0 ; i<N_sz ; i++)
		if (Ra_eI[i] != W0)
			Ra_eI[i] = W1 / Ra_eI[i];
	//sseVinv(Ra_eI, N_sz);
	if (Mp_sqInv) sseVsqrt(Ra_eI, N_sz, Ra_esqI);

	/* Calc inverse */
	wsSym	Ra_ret	= sseMDMt(Ra_E, N_sz, N_sz, Ra_eI);
	if (Mp_sqInv) {
		wsSym	Ra_sqRet = sseMDMt(Ra_E, N_sz, N_sz, Ra_esqI);
		*Mp_sqInv = *(new cSymMatrix(Ra_sqRet, N_sz));
	}
	return cSymMatrix(Ra_ret, N_sz);
}

cSymMatrix cSymMatrix::MMt(cSymMatrix &S)
{
	wsUint	N_sz	= row();
	if (N_sz!=S.row()) raise(MATERR_DIM_INCORRECT);
	/* Full SS */
	wsMat	Ra_SS	= sseSpS(get(), N_sz, S.get(), S.row());
	/* Do MS, sym form */
	wsSym	Ra_ret	= sym_sseMS(Ra_SS, N_sz, N_sz, get(), N_sz);
	sseUnmat(Ra_SS, N_sz);

	return cSymMatrix(Ra_ret, N_sz);
}

cSymMatrix cSymMatrix::Mt()
{
	wsSym Ra_ret = sseSpS(get(), row());
	return cSymMatrix(Ra_ret, row());
}

void cSymMatrix::setDiag(wsReal R_val)
{
	wsSym	R	= get();
	wsUint	N	= row();
	LOOP (i, N) R[i][i] = R_val;
}

cSymMatrix& cSymMatrix::neg()
{
	wsSym	R	= get();
	wsUint	N	= row();
	LOOP (i, N) sseVneg(R[i], i+1);
	return *this;
}

cSymMatrix& cSymMatrix::addDiag(wsRealCst R_val)
{
	wsSym	R	= get();
	wsUint	N	= row();
	LOOP (i, N) R[i][i] += R_val;
	return *this;
}

void cSymMatrix::subIself()
{
	wsSym	R	= get();
	wsUint	N	= row();
	LOOP (i, N) {
		sseVneg(R[i], i);
		R[i][i] = W1 - R[i][i];
	}
}

cSymMatrix cSymMatrix::subI(wsRealCst R_mul/*=W1*/)
{
	wsSym	V	= get();
	wsUint	N	= row();
	wsSym	R	= sseSymMat(N);
	if (R_mul == W1) LOOP (i, N) {
		sseVneg(V[i], i, R[i]);
		R[i][i] = W1 - V[i][i];
	} else LOOP (i, N) {
		sseVnegMul(R_mul, V[i], i, R[i]);
		R[i][i] = W1 - (V[i][i] * R_mul);
	}
	return cSymMatrix(R, N);
}

cVector cSymMatrix::r2v_ptr(int N_idx)
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
	return cVector(get()[N_finIdx], N_finIdx+1, 1);
}

cSymMatrix cSymMatrix::tM()
{
	return cSymMatrix(sseStS(get(), row()), row());
}

void cSymMatrix::rem()
{
//	LOG("rem() [%x] called\n", this);
	if (Ra_data && !B_ddealloc) sseUnmat(Ra_data, N_row);
	Ra_data	= NULL;
	N_row	= 0;
	N_col	= 0;
}

cVector cSymMatrix::sumR(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsReal *Ra_ret = sseSsumR(get(), row(), Rp_sumAll);
	if (!NA(R_mul)) sseVpC(Ra_ret, R_mul, Ra_ret, row());
	return cVector(Ra_ret, row());
}

cVector cSymMatrix::meanR()
{
	wsReal *Ra_ret = sseSmeanR(get(), row());
	return cVector(Ra_ret, row());
}

cVector cSymMatrix::sumC(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsReal *Ra_ret = sseSsumR(get(), row(), Rp_sumAll);
	if (!NA(R_mul)) sseVpC(Ra_ret, R_mul, Ra_ret, col());
	return cVector(Ra_ret, col());
}

void cSymMatrix::file(wsStrCst S_ext, wsUint N_prec/*=0*/)
{
	exportSymMat(S_ext, get(), row(), NULL, NULL, N_prec);
}

void cSymMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec/*=0*/)
{
	exportSymMat(S_ext, get(), row(), Sa_colNames, NULL, N_prec);
}

void cSymMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/)
{
	exportSymMat(S_ext, get(), row(), Sa_colNames, Sa_rowNames, N_prec);
}

/* CONFIRMED 130316 */
cSymMatrix& cSymMatrix::operator-(wsReal R)
{
	wsUint N_r = row();
	cSymMatrix *T = new cSymMatrix(N_r);
	/* Calc */
	for (wsUint i=0 ; i<N_r ; i++)
		sseVaC(get()[i], -R, T->get()[i], i+1);
	/* Return */
	return *T;
}

/* CONFIRMED 130316 */
cSymMatrix& cSymMatrix::operator+(wsReal R)
{
	wsUint N_r = row();
	cSymMatrix *T = new cSymMatrix(N_r);
	/* Calc */
	for (wsUint i=0 ; i<N_r ; i++)
		sseVaC(get()[i], R, T->get()[i], i+1);
	/* Return */
	return *T;
}

/* CONFIRMED 130318 */
cSymMatrix& cSymMatrix::operator*(wsReal R)
{
	wsUint N_r = row();
	cSymMatrix *T = new cSymMatrix(N_r);
	/* Calc */
	for (wsUint i=0 ; i<N_r ; i++)
		sseVpC(get()[i], R, T->get()[i], i+1);
	/* Return */
	return *T;
}

cVector cSymMatrix::operator*(cVector &V)
{
	/* col(this) == V.size() */
	if (col()!=V.size()) raise(MATERR_DIM_INCORRECT);
	wsReal *Ra_ret = sseSpV(cget(), col(), V.get());
	return cVector(Ra_ret, col());
}

void cSymMatrix::setR(wsStrCst S_varName)
{
#ifdef USE_R
	wsReal **Ra_full = sseSym2Mat(get(), row());
	R().sendMatrix(S_varName, Ra_full, row(), row());
#else
	halt("USE_R is required");
#endif
}

wsReal cSymMatrix::tr()
{
	wsUint N_r		= row();
	wsReal R_ret	= W0;
	wsReal **Ra_m	= get();

	for (wsUint i=0 ; i<N_r ; i++)
		R_ret += Ra_m[i][i];

	return R_ret;
}

cVector cSymMatrix::diag()
{
	wsUint	N_r		= row();
	wsVec	Ra_ret	= sseVector(N_r);
	wsMat	Ra_m	= get();

	LOOP (i, N_r) Ra_ret[i] = Ra_m[i][i];

	return cVector(Ra_ret, N_r);
}

cSymMatrix& cSymMatrix::clone()
{
	wsReal **Ra_ret	= sseSymMatP(row(), get());
	return *(new cSymMatrix(Ra_ret, row()));
}

cSymMatrix& cSymMatrix::sqrt()
{
	wsReal **Ra_ret	= sqrtMatrix(get(), row(), 1);
	return *(new cSymMatrix(Ra_ret, row()));
}

cSymMatrix& cSymMatrix::sqrt(wsReal **Ra_eVec, wsReal *Ra_eVal, wsUint N_sz)
{
	if (N_sz != row()) raise(MATERR_DIM_INCORRECT);
	wsReal **Ra_ret	= sqrtMatrixEV(Ra_eVec, Ra_eVal, row(), 1);
	return *(new cSymMatrix(Ra_ret, row()));
}

wsRealCst cSymMatrix::mean()
{
	wsUint N = row();
	return sum() / (wsReal)(N*N);
}

cVector cSymMatrix::eigen(cStdMatrix &M_eVec)
{
	char	B_fail	= 0;
	wsUint	N_sz	= row();
	wsMat	Ra_full = sseSym2Mat(get(), N_sz);
	wsMat	Ra_E	= NULL; /* P */
	wsReal*	Ra_eI	= EIGENDECOMPOSITION(Ra_full, N_sz, &Ra_E, 1);
	if (B_fail) halt("Eigendecomposition failed");
	M_eVec.init(N_sz, N_sz, Ra_E);
	sseUnmat(Ra_full, N_sz);

	return cVector(Ra_eI, N_sz);
}

cVector cSymMatrix::eigen()
{
	char	B_fail	= 0;
	wsUint	N_sz	= row();
	wsMat	Ra_full = sseSym2Mat(get(), N_sz);
	wsReal*	Ra_eI	= EIGENDECOMPOSITION(Ra_full, N_sz, NULL, 1);
	if (B_fail) halt("Eigendecomposition failed");
	sseUnmat(Ra_full, N_sz);

	return cVector(Ra_eI, N_sz);
}

} // End namespace ONETOOL
