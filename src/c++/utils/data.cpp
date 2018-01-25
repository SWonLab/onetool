#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "global/common.h"
#include "global/option.h"
#include "utils/util.h"
#include "utils/data.h"
#include "global/worker.h"
#include "global/io.h"
#include "utils/matrix.h"
#include "utils/ftp.h"
#include "input/stream.h"

#ifndef _WIN32
#	include	<sys/time.h>
#	include <sys/stat.h>
#	define	wisard_mkdir mkdir
#else
#	include <direct.h>
/* Ignores prm, do not needed in Windows because of lack of perm */
#	define	wisard_mkdir(pth,prm) _mkdir(pth)
#endif

namespace ONETOOL {

cSharedMatrix Y;
cSharedMatrix& S() { return Y; }

cDataStorage Q;
cDataStorage& D() { return Q; }

void cDataStorage::init(wsStrCst S_inpDP)
{
	if (!S_inpDP) halt("SYSERR: Data storage path is NULL");
	/* Check the existence of data storage */
	/* FIXME : impl */
	S_dataPath = strdup(S_inpDP);
}

cDataStorage::cDataStorage()
{
	S_dataPath = NULL;
}

cDataStorage::~cDataStorage()
{
	if (S_dataPath) free(S_dataPath);
}

static wsStrCst wsMkdir(wsStrCst dir, char B_noLast=0)
{
	static char tmp[MAX_PATH];
	char *ptmp = NULL;
	char *p = NULL;
	size_t len;

	sprintf(tmp, "%s", dir);
	len = strlen(tmp);
	if (tmp[len-1] == WISARD_DIR_CHAR)
		tmp[len-1] = 0;
	for (p=tmp+1 ; *p ; p++)
		if (*p == WISARD_DIR_CHAR) {
			*p = 0;
			wisard_mkdir(tmp, 0755);
			*p = WISARD_DIR_CHAR;
			ptmp = p;
		}
	if (B_noLast == 0)
		wisard_mkdir(tmp, 0755);
	else
		*ptmp = '\0';
	return tmp;
}

char* unixPath2system(wsStrCst S_path)
{
	char* S_ret = strdup(S_path);
	for (char* p=S_ret ; *p ; p++)
		if (*p == '/') *p = WISARD_DIR_CHAR;
	return S_ret;
}

cStrFile* cDataStorage::query(wsStrCst S_id)
{
	char	S_path[MAX_PATH];
	char*	S_sysPath = unixPath2system(S_id);

	/* Create directory if need */
	wsStrCst Sp_realDir = wsMkdir(S_sysPath, 1);

	sprintf(S_path, "%s" WISARD_DIR_LETTER "%s", S_dataPath, S_sysPath);
	cStrFile* Cp_ret = new cStrFile(S_sysPath, "", 1);
	if (Cp_ret->isFailed()) {
		/* Find file part */
		char* Sp_id = strdup(S_id);
		char* S_file = Sp_id+strlen(Sp_id)-1;
		for ( ; *S_file != '/' && Sp_id<S_file ; S_file--);
		if (*S_file == '/') {
			*S_file = '\0';
			S_file++;
		}

#ifdef USE_FTP
		cFtp C_ftp("ftp.ncbi.nlm.nih.gov", Sp_id);
		/* Move to target directory */
		if (S_dataPath && S_dataPath[0])
			sprintf(S_path, "CD %s" WISARD_DIR_LETTER "%s", S_dataPath,
				Sp_realDir);
		else
			sprintf(S_path, "CD %s", Sp_realDir);
		C_ftp.DoLCD(S_path);
		/* Download dataset */
		C_ftp.DoBinary();
		C_ftp.DoGet(S_file);
#endif

		delete Cp_ret;
		Cp_ret = new cStrFile(S_sysPath, "", 1);
		free(Sp_id);
	}

	free(S_sysPath);
	return Cp_ret;
}

/*
 *
 * Other functions
 * 
 */

bool svdcmp(wsReal **a, int m, int n,
	vector<double> & w, 
	vector<vector<double> > &v)
{
	bool flag;
	int i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z;
	double volatile temp;

	if (m == 0) halt("Internal problem in SVD function (no observations left?)");

	vector<double> rv1(n);
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+2;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i < m) {
			for (k=i;k<m;k++) scale += fabs(a[k][i]);
			if (scale != 0.0) {
				for (k=i;k<m;k++) {
					a[k][i] = (wsReal)(a[k][i]/scale);
					s += a[k][i]*a[k][i];
				}
				f		= a[i][i];
				g		= -SIGN(sqrt(s),f);
				h		= f*g-s;
				a[i][i]	= (wsReal)(f-g);
				for (j=l-1;j<n;j++) {
					for (s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<m;k++) a[k][j] += (wsReal)(f*a[k][i]);
				}
				for (k=i;k<m;k++) a[k][i] = (wsReal)(a[k][i]*scale);
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i+1 <= m && i+1 != n) {
			for (k=l-1;k<n;k++) scale += fabs(a[i][k]);
			if (scale != 0.0) {
				for (k=l-1;k<n;k++) {
					a[i][k] = (wsReal)(a[i][k]/scale);
					s += a[i][k]*a[i][k];
				}
				f			= a[i][l-1];
				g			= -SIGN(sqrt(s),f);
				h			= f*g-s;
				a[i][l-1]	= (wsReal)(f-g);
				for (k=l-1;k<n;k++) rv1[k]=a[i][k]/h;
				for (j=l-1;j<m;j++) {
					for (s=0.0,k=l-1;k<n;k++) s += a[j][k]*a[i][k];
					for (k=l-1;k<n;k++) a[j][k] += (wsReal)(s*rv1[k]);
				}
				for (k=l-1;k<n;k++) a[i][k] = (wsReal)(a[i][k]*scale);
			}
		}
		anorm=MAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g != 0.0) {
				for (j=l;j<n;j++)
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
					for (k=l;k<n;k++) v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=MIN(m,n)-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i][j]=0.0;
		if (g != 0.0) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<m;k++) s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i;k<m;k++) a[k][j] += (wsReal)(f*a[k][i]);
			}
			for (j=i;j<m;j++) a[j][i] = (wsReal)(a[j][i]*g);
		} else for (j=i;j<m;j++) a[j][i]=0.0;
		++a[i][i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=0;its<30;its++) {
			flag=true;
			for (l=k;l>=0;l--) {
				nm=l-1;
				temp=fabs(rv1[l])+anorm;
				if (temp == anorm) {
					flag=false;
					break;
				}
				temp=fabs(w[nm])+anorm;
				if (temp == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<k+1;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					temp = fabs(f)+anorm;
					if (temp == anorm) break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<m;j++) {
						y			= a[j][nm];
						z			= a[j][i];
						a[j][nm]	= (wsReal)(y*c+z*s);
						a[j][i]		= (wsReal)(z*c-y*s);
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j][k] = -v[j][k];
				}
				break;
			}
			if (its == 29) 
				return false; // cannot converge: multi-collinearity?
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g=g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<m;jj++) {
					y			= a[jj][j];
					z			= a[jj][i];
					a[jj][j]	= (wsReal)(y*c+z*s);
					a[jj][i]	= (wsReal)(z*c-y*s);
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	return true;
}

/* Matrix multiplication */

int thr_mult2Matrix_array(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	wsUint	i = *((wsUint *)Vp_data), j, k;
	xMatMul	*Rp_mm		= (xMatMul *)Vp_shareData;
	wsReal	*Rp_matR	= Rp_mm->Ra_ret[i];
	wsReal	*Rp_m1		= Rp_mm->Ra_mat[0][i];
	wsReal	**Rp_m2		= Rp_mm->Ra_mat[1];
	wsUint	N_mat1Col	= Rp_mm->Na_col[0];
	wsUint	N_mat2Col	= Rp_mm->Na_col[1];

	for (j=0 ; j<N_mat2Col ; j++) {
		wsReal R_sum = W0;

		for (k=0 ; k<N_mat1Col ; k++)
			R_sum += Rp_m1[k] * Rp_m2[k][j];

		Rp_matR[j] = R_sum;
	}

	return 0;
}

int thr_mult2Matrix_arraySSE(int N_idx, void *Vp_shareData,
	void *Vp_data, void *Vp_result)
{
	wsUint	i = *((wsUint *)Vp_data), j, k;
	xMatMul	*Rp_mm	= (xMatMul *)Vp_shareData;
	wsReal	*Rp_matR	= Rp_mm->Ra_ret[i];
	wsReal	*Rp_m1	= Rp_mm->Ra_mat[0][i];
	wsReal	**Rp_m2	= Rp_mm->Ra_mat[1];
	wsUint	N_col1	= Rp_mm->Na_col[0];
	wsUint	N_col2	= Rp_mm->Na_col[1];
	wsUint	N_med	= getMed(N_col2);

	for (j=0 ; j<N_col1 ; j++) {
		sse_t sse_m1_il = sseSet(Rp_m1[j]); // m1[i,j]

		/* SSE-enabled part */
		for (k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *sse_m2_lj	= (sse_t *)(&(Rp_m2[j][k])); // m2[j,k:(k+3)]
			sse_t *sse_mR		= (sse_t *)(&(Rp_matR[k])); // m3[i,k:(k+3)] (j)
			*sse_mR = sseAdd(*sse_mR, sseMul(sse_m1_il, *sse_m2_lj));
		}
		/* Rest part */
		for (k=N_med ; k<N_col2 ; k++)
			Rp_matR[k] += Rp_m1[j] * Rp_m2[j][k];
	}

	return 0;
}

int thr_multMM_arraySSE_equal(int N_idx, void *Vp_shareData,
	void *Vp_data, void *Vp_result)
{
	wsUint	i, j, k;
	wsUint	*pData	= (wsUint *)Vp_data;
	xMatMul	*Rp_mm	= (xMatMul *)Vp_shareData;
	wsReal	**Rp_m2	= Rp_mm->Ra_mat[1];
	wsUint	N_col1	= Rp_mm->Na_col[0];
	wsUint	N_col2	= Rp_mm->Na_col[1];

//	fprintf(stderr, "Thread %d : [%d-%d)\n", N_idx, pData[0], pData[1]);
	for (i=pData[0] ; i<pData[1] ; i++,pData[2]++) {
		wsReal	*Rp_matR	= Rp_mm->Ra_ret[i];
		wsReal	*Rp_m1		= Rp_mm->Ra_mat[0][i];
		for (j=0 ; j<N_col1 ; j++) {
			wsUint N_med	= 0;
#ifdef USE_SSE
			N_med			= getMed(N_col2);
			sse_t sse_m1_il	= sseSet(Rp_m1[j]); // m1[i,j]

			/* SSE-enabled part */
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_m2_lj	= (sse_t *)(&(Rp_m2[j][k])); // m2[j,k:(k+3)]
				sse_t *sse_mR		= (sse_t *)(&(Rp_matR[k])); // m3[i,k:(k+3)] (j)
				*sse_mR = sseAdd(*sse_mR, sseMul(sse_m1_il, *sse_m2_lj));
			}
#endif
			/* Rest part */
			for (k=N_med ; k<N_col2 ; k++)
				Rp_matR[k] += Rp_m1[j] * Rp_m2[j][k];
		}

// 		if (N_idx == 0 && (i%100) == 0) {
// 			wsUint s = 0;
// 			wsUint *xData;
// 			for (j=0,xData=(wsUint *)Vp_data ; j<OPT_NUMBER(thread) ; j++,xData+=3)
// 				s += xData[2];
// 			pverbose("%d rows processed...\r", s);
// 		}
	}
// 	if (N_idx == 0) {
// 		wsUint s = 0;
// 		wsUint *xData;
// 		for (j=0,xData=(wsUint *)Vp_data ; j<OPT_NUMBER(thread) ; j++,xData+=3)
// 			s += xData[2];
// 		pverbose("%d rows processed...\n", s);
// 	}

	return 0;
}

int thr_multSS_arraySSE_equal(int N_idx, void *Vp_shareData,
							  void *Vp_data, void *Vp_result)
{
	wsUint	i, a, b;
	wsUint	*pData	= (wsUint *)Vp_data;
	xMatMul	*Rp_mm	= (xMatMul *)Vp_shareData;
//	wsReal	**Rp_m2	= Rp_mm->Ra_mat[1];
	wsUint	N_col1	= Rp_mm->Na_col[0];
//	wsUint	N_col2	= Rp_mm->Na_col[1];

	wsReal **A = Rp_mm->Ra_mat[0];
	wsReal **B = Rp_mm->Ra_mat[1];
	wsReal **C = Rp_mm->Ra_ret;

	wsReal *Ra_iCol = NULL;
	sseMalloc(Ra_iCol, wsReal, N_col1);

	//	fprintf(stderr, "Thread %d : [%d-%d)\n", N_idx, pData[0], pData[1]);
	for (i=pData[0] ; i<pData[1] ; i++,pData[2]++) {
		/* Set i */
		for (a=i ; a<N_col1 ; a++)
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
		if (i < (N_col1-1)) {
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

		for(a=i+2 ; a<N_col1 ; a++) {
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
	// 	if (N_idx == 0) {
	// 		wsUint s = 0;
	// 		wsUint *xData;
	// 		for (j=0,xData=(wsUint *)Vp_data ; j<OPT_NUMBER(thread) ; j++,xData+=3)
	// 			s += xData[2];
	// 		pverbose("%d rows processed...\n", s);
	// 	}
	// 	
	sseFree(Ra_iCol);

	return 0;
}

int thr_multMMt_arraySSE_equal(int N_idx, void *Vp_shareData,
	void *Vp_data, void *Vp_result)
{
	wsUint	i, j, k;
	wsUint	*pData	= (wsUint *)Vp_data;
	xMatMul	*Rp_mm	= (xMatMul *)Vp_shareData;
	wsReal	**Rp_m2	= Rp_mm->Ra_mat[1];
	wsUint	N_col1	= Rp_mm->Na_col[0];
	wsUint	N_row2	= Rp_mm->Na_row[1];
	wsUint	N_med	= getMed(N_col1);

	if (Rp_mm->Ra_mat[0] == Rp_mm->Ra_mat[1]) {
		/* It can be symmetrized */
		for (i=pData[0] ; i<pData[1] ; i++,pData[2]++) {
			wsReal	*Rp_matR	= Rp_mm->Ra_ret[i];
			wsReal	*Rp_m1		= Rp_mm->Ra_mat[0][i];

			for (j=0 ; j<=i ; j++) {
				sse_t sse_sum = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_m1 = (sse_t *)(Rp_m1 + k);
					sse_t *sse_m2 = (sse_t *)(Rp_m2[j] + k);

					sse_sum = sseAdd(sse_sum, sseMul(*sse_m1, *sse_m2));
				}
				sseSum(sse_sum, Rp_matR[j]);
				for (k=N_med ; k<N_col1 ; k++)
					Rp_matR[j] += Rp_m1[k]*Rp_m2[j][k];
				/* Copy symm. data */
				Rp_mm->Ra_ret[j][i] = Rp_matR[j];
			}
		}
	} else {
		for (i=pData[0] ; i<pData[1] ; i++,pData[2]++) {
			wsReal	*Rp_matR	= Rp_mm->Ra_ret[i];
			wsReal	*Rp_m1		= Rp_mm->Ra_mat[0][i];

			for (j=0 ; j<N_row2 ; j++) {
				sse_t sse_sum = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_m1 = (sse_t *)(Rp_m1 + k);
					sse_t *sse_m2 = (sse_t *)(Rp_m2[j] + k);

					sse_sum = sseAdd(sse_sum, sseMul(*sse_m1, *sse_m2));
				}
				sseSum(sse_sum, Rp_matR[j]);
				for (k=N_med ; k<N_col1 ; k++)
					Rp_matR[j] += Rp_m1[k]*Rp_m2[j][k];
			}
		}
	}

	return 0;
}

int thr_mult3Matrix_array(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	wsUint	i = *((wsUint *)Vp_data);
	xMatMul	*Rp_mm		= (xMatMul *)Vp_shareData;
	wsReal	*Rp_matR	= Rp_mm->Ra_ret[i];
	wsReal	*Rp_m1		= Rp_mm->Ra_mat[0][i];
	wsReal	**Rp_m2		= Rp_mm->Ra_mat[1];
	wsReal	**Rp_m3		= Rp_mm->Ra_mat[2];
	wsUint	N_mat1Col	= Rp_mm->Na_col[0];
	wsUint	N_mat2Col	= Rp_mm->Na_col[1];
	wsUint	N_mat3Col	= Rp_mm->Na_col[2];

	/* Create intermediate buffer */
	wsReal	*Ra_interm = NULL;
	wsAlloc(Ra_interm, wsReal, N_mat2Col);

	for (wsUint l=0 ; l<N_mat2Col ; l++) {
		wsReal R_sum = 0.0f;

		for (wsUint j=0 ; j<N_mat1Col ; j++)
			R_sum += Rp_m1[j] * Rp_m2[j][l];

		Ra_interm[l] = R_sum;
	}

	for (wsUint n=0 ; n<N_mat3Col ; n++) {
		wsReal R_sum = 0.0f;

		for (wsUint l=0 ; l<N_mat2Col ; l++)
			R_sum += Ra_interm[l] * Rp_m3[l][n];

		Rp_matR[n] = R_sum;
	}

	DEALLOC(Ra_interm);
	return 0;
}

int thr_mult3Matrix_arraySSE_AB_C(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	wsUint	i = *((wsUint *)Vp_data), j, k;
	xMatMul	*Rp_mm		= (xMatMul *)Vp_shareData;
	wsReal	*Rp_matR	= Rp_mm->Ra_ret[i];
	wsReal	*Rp_m1		= Rp_mm->Ra_mat[0][i];
	wsReal	**Rp_m2		= Rp_mm->Ra_mat[1];
	wsReal	**Rp_m3		= Rp_mm->Ra_mat[2];
	wsUint	N_mat1Col	= Rp_mm->Na_col[0];
	wsUint	N_mat2Col	= Rp_mm->Na_col[1];
	wsUint	N_mat3Col	= Rp_mm->Na_col[2];
	wsUint	N_med		= getMed(N_mat2Col);

	/* Create intermediate buffer */
	wsReal	*Ra_interm = NULL;
	sseCalloc(Ra_interm, wsReal, N_mat2Col);

	/* Calculate A[i,]*B */
	for (j=0 ; j<N_mat1Col ; j++) {
		sse_t sse_m1_il = sseSet(Rp_m1[j]);

		/* SSE-enabled part */
		for (k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *sse_m2_lj	= (sse_t *)(&(Rp_m2[j][k]));
			sse_t *sse_mR		= (sse_t *)(&(Ra_interm[k]));
			*sse_mR = sseAdd(*sse_mR, sseMul(sse_m1_il, *sse_m2_lj));
		}
		/* Rest part */
		for (k=N_med ; k<N_mat2Col ; k++)
			Ra_interm[k] += Rp_m1[j] * Rp_m2[j][k];
	}

	if (N_mat3Col < sseJmp) {
		for (j=0 ; j<N_mat3Col ; j++) {
			wsReal R_sum = 0.0f;

			for (k=0 ; k<N_mat2Col ; k++)
				R_sum += Ra_interm[k] * Rp_m3[k][j];

			Rp_matR[j] = R_sum;
		}
	} else {
		N_med = getMed(N_mat3Col);
		for (j=0 ; j<N_mat2Col ; j++) {
			wsReal	R_m		= Ra_interm[j];
			sse_t	sse_m	= sseSet(R_m);

			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_m2_lj	= (sse_t *)(&(Rp_m3[j][k]));
				sse_t *sse_mR		= (sse_t *)(&(Rp_matR[k]));
				*sse_mR = sseAdd(*sse_mR, sseMul(sse_m, *sse_m2_lj));
			}
			for (k=N_med ; k<N_mat3Col ; k++)
				Rp_matR[k] += R_m * Rp_m3[j][k];
		}
	}

	sseFree(Ra_interm);
	return 0;
}

int forAllRow_array(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx = (int *)Vp_thrData;
	xMatMul	*Rp_mm	= (xMatMul *)Vp_shareData;
	wsUint	N_r1	= Rp_mm->Na_row[0];
	static wsUint N_idxRow	= 0;//Cp_Matrices[1]->N_row;
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		N_idxRow = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++) {
			if ((N_idxRow+i) >= N_r1)
				break;
			Np_idx[i] = N_idxRow+i;
			//printf("[%d] %d row allocated\n", i, N_idxRow+i);
			if (N_idxRow && ((N_idxRow+i)%10) == 0)
				notice("%d/%d Rows processed\r", N_idxRow+i, N_r1);
		}
		N_idxRow += N_thread;
		break;
	case DISTFUN_UNINIT:
		N_idxRow = 0;
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt_fmt(WISARD_SYST_INVL_DISTFUNC_CMD, N_mode);
		break;
	}
	return i;
}

int forAllRow_array_equal(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx = (int *)Vp_thrData;
	xMatMul	*Rp_mm	= (xMatMul *)Vp_shareData;
	wsUint	N_r1	= Rp_mm->Na_row[0];
	if ((int)N_r1 < N_thread)
		N_thread = (int)N_r1;
	int		i=0;
	wsReal	j=W0;
	static wsUint
		N_proc	= 0;
	wsReal	R_szDiv	= (wsReal)N_r1/(wsReal)N_thread;
	int		N_ret = N_proc != 0 ? 0 : 1;

	switch (N_mode) {
	case DISTFUN_INIT:
		i = 0;
		j = W0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		if (N_ret == 0)
			return 0;
		for ( ; i<N_thread ; i++,j+=R_szDiv) {
			wsUint E = (wsUint)(j+R_szDiv+REAL_CONST(0.5));
			if (E > N_r1)
				E = N_r1;
			Np_idx[3*i] = (wsUint)(j+REAL_CONST(0.5));
			Np_idx[3*i+1] = E;
			Np_idx[3*i+2] = 0;
		}
		N_proc = N_ret = N_thread;
		break;
	case DISTFUN_UNINIT:
		N_proc = 0;
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt_fmt(WISARD_SYST_INVL_DISTFUNC_CMD, N_mode);
		break;
	}

	return N_ret;
}

wsReal** multThrRowMM(wsReal **Rp_mat1, wsUint N_r1, wsUint N_c1,
	wsReal **Rp_mat2,  wsUint N_r2, wsUint N_c2)
{
	if (N_c1 != N_r2)
		halt_fmt(WISARD_SYST_INVL_DIM, N_c1, N_r2);

	/* Make return matrix */
	wsReal **Ra_ret = NULL;
	wsCalloc(Ra_ret, wsReal*, N_r1);
	for (wsUint i=0 ; i<N_r1 ; i++) {
		wsCalloc(Ra_ret[i], wsReal, N_c2);
	}

	//	printf("\nMult2 START\n");
	wsReal	**Rp_mat[]	= { Rp_mat1, Rp_mat2 };
	wsUint	Np_row[]	= { N_r1, N_r2 };
	wsUint	Np_col[]	= { N_c1, N_c2 };
	xMatMul	X_mul		= { Ra_ret, 2, Rp_mat, Np_row, Np_col };

	WORKER().run(thr_mult2Matrix_array, forAllRow_array, &X_mul, NULL);
	//	printf("\nMult2 END\n");

	return Ra_ret;
}

wsReal** multMM(wsReal **Rp_mat1, wsUint N_r1, wsUint N_c1,
	wsReal **Rp_mat2,  wsUint N_r2, wsUint N_c2)
{
	if (N_c1 != N_r2)
		halt_fmt(WISARD_SYST_INVL_DIM, N_c1, N_r2);

	wsUint	i, j, k;
	wsReal	**Ra_ret = NULL;
	wsAlloc(Ra_ret, wsReal*, N_r1);

	for (i=0 ; i<N_r1 ; i++) {
		wsAlloc(Ra_ret[i], wsReal, N_c2);
		wsReal	*Rp_matR	= Ra_ret[i];
		wsReal	*Rp_m1		= Rp_mat1[i];

		for (j=0 ; j<N_c2 ; j++) {
			wsReal R_sum = W0;

			for (k=0 ; k<N_c1 ; k++)
				R_sum += Rp_m1[k] * Rp_mat2[k][j];

			Rp_matR[j] = R_sum;
		}
	}

	return Ra_ret;
}

wsReal** multMMt(wsReal **Rp_mat1, wsUint N_r1, wsUint N_c1,
				wsReal **Rp_mat2,  wsUint N_r2, wsUint N_c2)
{
	if (N_c1 != N_c2)
		halt_fmt(WISARD_SYST_INVL_DIM, N_c1, N_c2);

	wsUint	i, j, k;
	wsReal	**Ra_ret = NULL;
	wsAlloc(Ra_ret, wsReal*, N_r1);

	for (i=0 ; i<N_r1 ; i++) {
		wsAlloc(Ra_ret[i], wsReal, N_r2);
		wsReal	*Rp_matR	= Ra_ret[i];
		wsReal	*Rp_m1		= Rp_mat1[i];

		for (j=0 ; j<N_r2 ; j++) {
			wsReal	R_sum	= W0;
			wsReal	*Rp_m2	= Rp_mat2[j];

			for (k=0 ; k<N_c1 ; k++)
				R_sum += Rp_m1[k] * Rp_m2[k];

			Rp_matR[j] = R_sum;
		}
	}

	return Ra_ret;
}

wsReal** sseThrRowMM(wsReal **Rp_m1, wsUint N_r1, wsUint N_c1,
	wsReal **Rp_m2, wsUint N_r2, wsUint N_c2)
{
    if (OPT_NUMBER(thread) == 1)
        return sseMpM(Rp_m1, N_r1, N_c1, Rp_m2, N_r2, N_c2);
	if (N_c1 != N_r2)
		halt_fmt(WISARD_SYST_INVL_DIM, N_c1, N_r2);

	/* Make return matrix */
	wsReal **Ra_ret = sseEmptyMatrix(N_r1, N_c2);

	wsReal	**Rp_mat[]	= { Rp_m1, Rp_m2 };
	wsUint	Np_row[]	= { N_r1, N_r2 };
	wsUint	Np_col[]	= { N_c1, N_c2 };
	xMatMul	X_mul		= { Ra_ret, 2, Rp_mat, Np_row, Np_col };

#ifdef USE_NEQ
	WORKER().run(thr_mult2Matrix_arraySSE, forAllRow_array,
		&X_mul, NULL, sizeof(int)*3);
#else
	WORKER().run(thr_multMM_arraySSE_equal, forAllRow_array_equal,
		&X_mul, NULL, sizeof(int)*3);
#endif
	//	printf("\nMult2 END\n");

	return Ra_ret;
}

wsReal** sseThrRowSS(wsSym Rp_mat1, wsUint N_s1, wsSym Rp_mat2, wsUint N_s2)
{
    if (OPT_NUMBER(thread) == 1)
        return sseSpS(Rp_mat1, N_s1, Rp_mat2, N_s2);
	if (N_s1 != N_s2)
		halt_fmt(WISARD_SYST_INVL_DIM, N_s1, N_s2);
	wsUint N_sz = N_s1;

	/* Make return matrix */
	wsReal **Ra_ret = sseEmptyMatrix(N_sz, N_sz);

	wsReal	**Rp_mat[]	= { Rp_mat1, Rp_mat2 };
	wsUint	Np_row[]	= { N_sz, N_sz };
	wsUint	Np_col[]	= { N_sz, N_sz };
	xMatMul	X_mul		= { Ra_ret, 2, Rp_mat, Np_row, Np_col };

#ifdef USE_NEQ
	WORKER().run(thr_mult2Matrix_arraySSE, forAllRow_array,
		&X_mul, NULL, sizeof(int)*3);
#else
	WORKER().run(thr_multSS_arraySSE_equal, forAllRow_array_equal,
		&X_mul, NULL, sizeof(int)*3);
#endif
	//	printf("\nMult2 END\n");

	return Ra_ret;
}

wsReal** sseThrRowMMt(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2/*=NULL*/,  wsUint N_mat2Row/*=0xffffffff*/,
	wsUint N_mat2Col/*=0xffffffff*/)
{
	if (Rp_mat2 == NULL) {
		Rp_mat2 = Rp_mat1;
		N_mat2Row = N_mat1Row;
		N_mat2Col = N_mat1Col;
	}
	if (N_mat1Col != N_mat2Col)
		halt_fmt(WISARD_SYST_INVL_DIM, N_mat1Col, N_mat2Col);

	/* Make return matrix */
	wsReal **Ra_ret = sseEmptyMatrix(N_mat1Row, N_mat2Row);

	wsReal	**Rp_mat[]	= { Rp_mat1, Rp_mat2 };
	wsUint	Np_row[]	= { N_mat1Row, N_mat2Row };
	wsUint	Np_col[]	= { N_mat1Col, N_mat2Col };
	xMatMul	X_mul		= { Ra_ret, 2, Rp_mat, Np_row, Np_col };

	WORKER().run(thr_multMMt_arraySSE_equal, forAllRow_array_equal,
		&X_mul, NULL, sizeof(int)*3);
	//	printf("\nMult2 END\n");

	return Ra_ret;
}

wsReal** multMMM(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2,  wsUint N_mat2Row, wsUint N_mat2Col,
	wsReal **Rp_mat3,  wsUint N_mat3Row, wsUint N_mat3Col, int B_export)
{
	if (N_mat1Col != N_mat2Row)
		halt_fmt(WISARD_SYST_INVL_DIM, N_mat1Col, N_mat2Row);
	if (N_mat2Col != N_mat3Row)
		halt_fmt(WISARD_SYST_INVL_DIM, N_mat2Col, N_mat3Row);

	/* Make return matrix */
	wsReal **Ra_ret = NULL;
	wsCalloc(Ra_ret, wsReal*, N_mat1Row);
	for (wsUint i=0 ; i<N_mat1Row ; i++) {
		wsCalloc(Ra_ret[i], wsReal, N_mat3Col);
	}

	//	printf("\nMult3 START\n");
	wsReal	**Rp_mat[]	= { Rp_mat1, Rp_mat2, Rp_mat3 };
	wsUint	Np_row[]	= { N_mat1Row, N_mat2Row, N_mat3Row };
	wsUint	Np_col[]	= { N_mat1Col, N_mat2Col, N_mat3Col };
	xMatMul	X_mul		= { Ra_ret, 3, Rp_mat, Np_row, Np_col };
	if (B_export == 1) {
		exportMatrix("ma1", Rp_mat1, N_mat1Row, N_mat1Col);
		exportMatrix("ma2", Rp_mat2, N_mat2Row, N_mat2Col);
		exportMatrix("ma3", Rp_mat3, N_mat3Row, N_mat3Col);
	}
	WORKER().run(thr_mult3Matrix_array, forAllRow_array, &X_mul, NULL);
	if (B_export == 1)
		exportMatrix("maR", Ra_ret, N_mat1Row, N_mat3Col);
	//	printf("\nMult3 END\n");

	return Ra_ret;
}

wsReal** mult3MatrixSSE2(wsReal **Ra_m1, wsUint N_r1, wsUint N_c1,
	wsReal **Ra_m2,  wsUint N_r2, wsUint N_c2,
	wsReal **Ra_m3,  wsUint N_r3, wsUint N_c3)
{
	/* Dimension check */
	if (N_c1 != N_r2)
		halt_fmt(WISARD_SYST_INVL_DIM, N_c1, N_r2);
	if (N_c2 != N_r3)
		halt_fmt(WISARD_SYST_INVL_DIM, N_c2, N_r3);

	double R_s1 = (double)N_r1*(double)N_c2*(double)(N_c1 + N_c3);
	double R_s2 = (double)N_c1*(double)N_c3*(double)(N_r1 + N_c2);

	/* Make return matrix */
	wsReal **Ra_ret = sseEmptyMatrix(N_r1, N_c3);

	wsReal	**Rp_mat[]	= { Ra_m1, Ra_m2, Ra_m3 };
	wsUint	Np_row[]	= { N_r1, N_r2, N_r3 };
	wsUint	Np_col[]	= { N_c1, N_c2, N_c3 };
	xMatMul	X_mul		= { Ra_ret, 3, Rp_mat, Np_row, Np_col };

	/* (AB)C */
	if (R_s1 < R_s2)
		WORKER().run(thr_mult3Matrix_arraySSE_AB_C, forAllRow_array, &X_mul, NULL);
	else {
		notice("Warning, %.2g times inefficient\n", R_s1/R_s2);
		WORKER().run(thr_mult3Matrix_arraySSE_AB_C, forAllRow_array, &X_mul, NULL);
	}

	return Ra_ret;
}

wsReal** mult3MatrixSSE(wsReal **Rp_mat1, wsUint N_r1, wsUint N_c1,
	wsReal **Rp_mat2,  wsUint N_r2, wsUint N_c2,
	wsReal **Rp_mat3,  wsUint N_r3, wsUint N_c3, int B_export)
{
	/* DImension check */
	if (N_c1 != N_r2) halt_fmt(WISARD_SYST_INVL_DIM, N_c1, N_r2);
	if (N_c2 != N_r3) halt_fmt(WISARD_SYST_INVL_DIM, N_c2, N_r3);

	/* Make return matrix */
	wsReal **Ra_ret = NULL;
	wsCalloc(Ra_ret, wsReal*, N_r1);
	for (wsUint i=0 ; i<N_r1 ; i++) {
		sseCalloc(Ra_ret[i], wsReal, N_c3);
	}

	//	printf("\nMult3 START\n");
	wsReal	**Rp_mat[]	= { Rp_mat1, Rp_mat2, Rp_mat3 };
	wsUint	Np_row[]	= { N_r1, N_r2, N_r3 };
	wsUint	Np_col[]	= { N_c1, N_c2, N_c3 };
	xMatMul	X_mul		= { Ra_ret, 3, Rp_mat, Np_row, Np_col };

	if (B_export == 1) {
		exportMatrix("ma1", Rp_mat1, N_r1, N_c1);
		exportMatrix("ma2", Rp_mat2, N_r2, N_c2);
		exportMatrix("ma3", Rp_mat3, N_r3, N_c3);
	}
	WORKER().run(thr_mult3Matrix_array, forAllRow_array, &X_mul, NULL);
	if (B_export == 1)
		exportMatrix("maR", Ra_ret, N_r1, N_c3);
	//	printf("\nMult3 END\n");

	return Ra_ret;
}

wsReal** _mult2Matrix(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2, wsUint N_mat2Row, wsUint N_mat2Col)
{
	/* Dimension checking */
	if (N_mat1Col != N_mat2Row)
		halt_fmt(WISARD_SYST_INVL_DIM, N_mat1Col, N_mat2Row);

	wsReal **Rp_resMatrix = NULL;
	wsCalloc(Rp_resMatrix, wsReal*, N_mat1Row);
	for (wsUint i=0 ; i<N_mat1Row ; i++) {
		wsCalloc(Rp_resMatrix[i], wsReal, N_mat2Col);
	}

	for (wsUint i=0 ; i<N_mat1Row ; i++) {
		wsReal	*Rp_matR = Rp_resMatrix[i];

		for (wsUint l=0 ; l<N_mat2Col ; l++) {
			wsReal R_sum = 0.0f;

			for (wsUint j=0 ; j<N_mat1Col ; j++)
				R_sum += Rp_mat1[i][j] * Rp_mat2[j][l];

			Rp_matR[l] = R_sum;
			//	printf("%d,%d -> %g\n", i, l, R_sum);
		}
	}

	return Rp_resMatrix;
}

wsReal** _mult3Matrix(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2, wsUint N_mat2Row, wsUint N_mat2Col,
	wsReal **Rp_mat3, wsUint N_mat3Row, wsUint N_mat3Col)
{
	if (N_mat1Col != N_mat2Row)
		halt_fmt(WISARD_SYST_INVL_DIM, N_mat1Col, N_mat2Row);
	if (N_mat2Col != N_mat3Row)
		halt_fmt(WISARD_SYST_INVL_DIM, N_mat2Col, N_mat3Row);

	wsReal **Rp_resMatrix = NULL;
	wsCalloc(Rp_resMatrix, wsReal*, N_mat1Row);
	for (wsUint i=0 ; i<N_mat1Row ; i++) {
		wsCalloc(Rp_resMatrix[i], wsReal, N_mat3Col);
	}

	wsReal *Ra_inter = NULL;
	wsCalloc(Ra_inter, wsReal, N_mat2Col);

	for (wsUint i=0 ; i<N_mat1Row ; i++) {
		wsReal	*Rp_matR = Rp_resMatrix[i];

		/* Calculate intermediate */
		for (wsUint l=0 ; l<N_mat2Col ; l++) {
			wsReal R_sum = 0.0f;

			for (wsUint j=0 ; j<N_mat1Col ; j++)
				R_sum += Rp_mat1[i][j] * Rp_mat2[j][l];

			Ra_inter[l] = R_sum;
		}

		for (wsUint l=0 ; l<N_mat3Col ; l++) {
			wsReal R_sum = 0.0f;

			for (wsUint j=0 ; j<N_mat2Col ; j++)
				R_sum += Ra_inter[j] * Rp_mat3[j][l];

			Rp_matR[l] = R_sum;
		}
	}
	DEALLOC(Ra_inter);

	return Rp_resMatrix;
}

void getIdentityError(wsReal **Ra_matrix, wsUint N_sz)
{
	wsReal R_max=Ra_matrix[0][0]-1.0f, R_min=Ra_matrix[0][0]-1.0f;

	for (wsUint i=0 ; i<N_sz ; i++) {
		for (wsUint j=0 ; j<N_sz ; j++) {
			wsReal R_val = (i == j) ? Ra_matrix[i][j]-1.0f : Ra_matrix[i][j];

			if (R_val > R_max)
				R_max = R_val;
			if (R_val < R_min)
				R_min = R_val;
		}
	}

	LOG("Matrix error is (%g) %g~%g\n", fabs(R_min)+fabs(R_max), R_min, R_max);
}

wsReal multVV(wsRealCst *Ra_v1, wsUint N_sz, wsRealCst *Ra_v2/*=NULL*/)
{
	wsReal R_ret = W0;
	if (Ra_v2 == NULL)
		Ra_v2 = Ra_v1;

	for (wsUint i=0 ; i<N_sz ; i++)
		R_ret += Ra_v1[i]*Ra_v2[i];

	return R_ret;
}

cVariantMap::cVariantMap()
{
	Cp_IO = NULL;
	/* Initialize to NULL */
	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) {
		Na_chr[i] = 0;
		Xa_chr[i] = NULL;
	}
}

cVariantMap::cVariantMap(cIO *Cp_inpIO)
{
	init(Cp_inpIO);
}

cVariantMap::~cVariantMap()
{
	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) {
		if (Xa_chr[i]) DEALLOC(Xa_chr[i]);
		if (Va_chr[i]) DEALLOC(Va_chr[i]);
	}
}

void cVariantMap::init(cIO *Cp_inpIO)
{
	Cp_IO = Cp_inpIO;

	/* Get markers and counting */
	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++)
		Na_chr[i] = 0;
	wsUint Na_ins[MAX_NCHR+1] = { 0, };
	vVariant&	Xa_snp = Cp_IO->getVariant();
	FOREACH (vVariant_it, Xa_snp, i) {
		int chr = i->chr < 0 ? 0 : i->chr;
		Na_chr[chr]++;
	}

	/* Now allocating the memory */
	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) {
		wsAlloc(Xa_chr[i], xVariant*, Na_chr[i]);
		wsAlloc(Va_chr[i], wsUint, Na_chr[i]);
	}

	/* Set the point */
	cTimer t;
	t.start();
	wsUint j = 0;
	FOREACHDO (vVariant_it, Xa_snp, i, j++) {
		int chr = i->chr < 0 ? 0 : i->chr;

		xVariant**	X_cc	= Xa_chr[chr];
		wsUint*		V_cc	= Va_chr[chr];

		/* Search */
		wsUint	N_s		= 0, N_e = Na_ins[chr];
		wsUint	N_pT	= i->pos;
		wsUint	N_ip	= Na_ins[chr]>>1;
		if (Na_ins[chr]) {
			if (N_pT > X_cc[N_ip]->pos)
				N_s = N_ip+1;
			else if (N_pT < X_cc[N_ip]->pos)
				N_e = N_ip;
			else
				N_s = N_e = N_ip;

			/* Until the insertion point is made */
			while (N_s < N_e) {
				N_ip = (N_s + N_e) >> 1;
				 
				if (N_pT > X_cc[N_ip]->pos)
					N_s = N_ip+1;
				else if (N_pT < X_cc[N_ip]->pos)
					N_e = N_ip;
				else
					N_s = N_e = N_ip;
			}

			/* Insert to N_s */
			int N_mov = Na_ins[chr] - N_s;
			if (N_mov > 0) {
				memmove(X_cc + N_s + 1, X_cc + N_s, sizeof(xVariant*)*N_mov);
				memmove(V_cc + N_s + 1, V_cc + N_s, sizeof(wsUint)*N_mov);
			}
		}
		X_cc[N_s] = &(*i);
		V_cc[N_s] = j;
		Na_ins[chr]++;
	}
	if (OPT_ENABLED(verbose)) {
		LOGoutput("Mapping status is exported to [%s.variant.map.res]\n",
			OPT_STRING(out));
		cExporter *C_ms = cExporter::summon("variant.map.res");
		for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) {
			if (!Na_ins[i]) continue;

			C_ms->put(getChrName2(i));
			for (wsUint j=0 ; j<Na_ins[i] ; j++) {
				if (Xa_chr[i][j]->name != Xa_snp[Va_chr[i][j]].name)
					halt("Variant [%s] have different index [%d]\n", Xa_chr[i][j]->name,
						Va_chr[i][j]);
				C_ms->fmt("	%s[%d]", Xa_chr[i][j]->name, Xa_chr[i][j]->pos);
			}
			C_ms->put("\n");
		}
		delete C_ms;
	}
	lverbose("Variant categorization took [%s]\n", t.getReadable());
}

xVariant*** cVariantMap::getMarkerMap(wsUint **Np_chr)
{
	if (Np_chr) *Np_chr = Na_chr;
	return Xa_chr;
}

wsUint** cVariantMap::getVariantPosMap(wsUint **Np_chr)
{
	if (Np_chr) *Np_chr = Na_chr;
	return Va_chr;
}

} // End namespace ONETOOL
