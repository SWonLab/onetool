#include "global/worker.h"
#include "utils/util.h"
#include "utils/matrix.h"
#include <stdarg.h>

namespace ONETOOL {

void exportMatrix(const char *S_ext, char **Na_mat, wsUint N_row, wsUint N_col)
{
	wsUint		i;
//	char		S_fmt[32];
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	LOGoutput("Export matrix [%s.%s] (%d * %d)\n", OPT_STRING(out), S_ext, N_row, N_col);
	// 	if (Xp_header) {
	// 		for (i=0 ; i<N_col ; i++)
	// 			Cp_mat->fmt("%s	", (*Xp_header)[i].c_str());
	// 		Cp_mat->put("\n");
	// 	}

	for (i = 0; i < N_row; i++) {
		// 		if (Xp_rowNames)
		// 			Cp_mat->fmt("%s	", (*Xp_rowNames)[i].c_str());
		for (wsUint j = 0; j < N_col; j++)
			isMissing(Na_mat[i][j]) ? Cp_mat->put("NA	") :
			Cp_mat->fmt("%d	", Na_mat[i][j]);
		Cp_mat->put("\n");
	}
	delete Cp_mat;
}

void exportMatrix(const char *S_ext, vector<wsFvec>& Na_mat, wsUint N_col,
	vStr *Xp_header/*=NULL*/)
{
	wsUint		N_row = (wsUint)Na_mat.size();
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	LOGoutput("Export matrix [%s.%s] (%d * %d)\n", OPT_STRING(out), S_ext, N_row, N_col);
	// 	if (Xp_header) {
	// 		for (i=0 ; i<N_col ; i++)
	// 			Cp_mat->fmt("%s	", (*Xp_header)[i].c_str());
	// 		Cp_mat->put("\n");
	// 	}
	
	if (Xp_header) {
		FOREACH (vStr_it, *Xp_header, i) Cp_mat->fmt("%s	", i->c_str());
		Cp_mat->put("\n");
	}

	FOREACH (vector<wsFvec>::iterator, Na_mat, i) {
		// 		if (Xp_rowNames)
		// 			Cp_mat->fmt("%s	", (*Xp_rowNames)[i].c_str());
		for (wsUint j = 0; j < N_col; j++)
			isMissingReal((*i)[j]) ? Cp_mat->put("NA	") :
			Cp_mat->fmt("%g	", (*i)[j]);
		Cp_mat->put("\n");
	}
	delete Cp_mat;
}

void exportMatrix(const char *S_ext, wsMat Rp_mat, wsUint N_row, wsUint N_col)
{
	exportMatrix(S_ext, Rp_mat, N_row, N_col, NULL, NULL, 0, NULL);
}

void exportMatrix(const char *S_ext, vector<wsVec>& Rv_mat, wsUint N_row, wsUint N_col)
{
	exportMatrix(S_ext, Rv_mat, N_col, NULL, NULL, 0, NULL);
}

void exportMatrix(const char *S_ext, wsMat Rp_mat, wsUint N_row, wsUint N_col,
	vStr *Xp_header)
{
	exportMatrix(S_ext, Rp_mat, N_row, N_col, Xp_header, NULL, 0, NULL);
}

void exportMatrix(const char *S_ext, wsMat Rp_mat, wsUint N_row, wsUint N_col,
	vStr *Xp_header, vStr *Xp_rowNames, wsReal *Rp_miss/*=NULL*/)
{
	exportMatrix(S_ext, Rp_mat, N_row, N_col, Xp_header, Xp_rowNames, 0, Rp_miss);
}

void exportMatrix(const char *S_ext, vector<wsVec>& Rv_mat, wsUint N_col)
{
	exportMatrix(S_ext, Rv_mat, N_col, NULL, NULL, 0, NULL);
}

#ifdef USE_DBL
void exportMatrix(wsStrCst S_ext, wsFmat Rp_mat, wsUint N_row, wsUint N_col)
{
	exportMatrix(S_ext, Rp_mat, N_row, N_col, NULL, NULL, 0, NULL);
}

void exportMatrix(wsStrCst S_ext, wsFmat Rp_mat, wsUint N_row, wsUint N_col,
	vStr *Xp_header, vStr *Xp_rowNames, wsUint N_prec, wsReal *Rp_miss/*=NULL*/)
{
	wsUint		i;
	char		S_fmt[32];
	/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	if (N_prec == 0)
		sprintf(S_fmt, "%%g	");
	else
		sprintf(S_fmt, "%%.%de	", N_prec);

	LOGoutput("Export matrix [%s.%s] [%d * %d]\n", OPT_STRING(out), S_ext, N_row, N_col);
	if (Xp_header) {
		for (i = 0; i < N_col; i++)
			Cp_mat->fmt("%s	", (*Xp_header)[i].c_str());
		Cp_mat->put("\n");
	}

	if (Rp_miss) for (i = 0; i < N_row; i++) {
		if (Xp_rowNames)
			Cp_mat->fmt("%s	", (*Xp_rowNames)[i].c_str());
		for (wsUint j = 0; j < N_col; j++)
			Rp_mat[i][j] == *Rp_miss ? Cp_mat->put("NA	") : Cp_mat->fmt(S_fmt, Rp_mat[i][j]);
		Cp_mat->put("\n");
	}
	else for (i = 0; i < N_row; i++) {
		if (Xp_rowNames)
			Cp_mat->fmt("%s	", (*Xp_rowNames)[i].c_str());
		for (wsUint j = 0; j < N_col; j++)
			Cp_mat->fmt(S_fmt, Rp_mat[i][j]);
		Cp_mat->put("\n");
	}

	//	pverbose("Trying to delte..");
	delete Cp_mat;
	//	pverbose("Trying to delte..OK\n");
}
#endif

void exportMatrix(wsStrCst S_ext, wsMat Rp_mat, wsUint N_row, wsUint N_col,
	vStr *Xp_header, vStr *Xp_rowNames, wsUint N_prec, wsReal *Rp_miss/*=NULL*/)
{
	wsUint		i;
	char		S_fmt[32];
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	if (N_prec == 0)
		sprintf(S_fmt, "%%g	");
	else
		sprintf(S_fmt, "%%.%de	", N_prec);

	LOGoutput("Export matrix [%s.%s] [%d * %d]\n", OPT_STRING(out), S_ext, N_row, N_col);
	if (Xp_header) {
		for (i = 0; i < N_col; i++)
			Cp_mat->fmt("%s	", (*Xp_header)[i].c_str());
		Cp_mat->put("\n");
	}

	if (Rp_miss) for (i = 0; i < N_row; i++) {
		if (Xp_rowNames)
			Cp_mat->fmt("%s	", (*Xp_rowNames)[i].c_str());
		for (wsUint j = 0; j < N_col; j++)
			Rp_mat[i][j] == *Rp_miss ? Cp_mat->put("NA	") : Cp_mat->fmt(S_fmt, Rp_mat[i][j]);
		Cp_mat->put("\n");
	}
	else for (i = 0; i < N_row; i++) {
		if (Xp_rowNames)
			Cp_mat->fmt("%s	", (*Xp_rowNames)[i].c_str());
		for (wsUint j = 0; j < N_col; j++)
			Cp_mat->fmt(S_fmt, Rp_mat[i][j]);
		Cp_mat->put("\n");
	}

	//	pverbose("Trying to delte..");
	delete Cp_mat;
	//	pverbose("Trying to delte..OK\n");
}

void exportMatrix(wsStrCst S_ext, vector<wsVec>& Rv_mat, wsUint N_col,
	vStr *Xp_header, vStr *Xp_rowNames, wsUint N_prec, wsReal *Rp_miss/*=NULL*/)
{
	wsUintCst	N_row = (wsUint)Rv_mat.size();
	wsUint		i;
	char		S_fmt[32];
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	if (N_prec == 0)
		sprintf(S_fmt, "%%g	");
	else
		sprintf(S_fmt, "%%.%de	", N_prec);

	LOGoutput("Export matrix [%s.%s] [%d * %d]\n", OPT_STRING(out), S_ext, N_row, N_col);
	if (Xp_header) {
		LOOP (i, N_col)
			Cp_mat->fmt("%s	", (*Xp_header)[i].c_str());
		Cp_mat->put("\n");
	}

	i=0;
	if (Rp_miss) FOREACHDO (vector<wsVec>::iterator, Rv_mat, I, i++) {
		wsVec Ra_v = *I;
		if (Xp_rowNames)
			Cp_mat->fmt("%s	", (*Xp_rowNames)[i].c_str());
		LOOP (j, N_col)
			Ra_v[j] == *Rp_miss ? Cp_mat->put("NA	") : Cp_mat->fmt(S_fmt, Ra_v[j]);
		Cp_mat->put("\n");
	} else FOREACHDO (vector<wsVec>::iterator, Rv_mat, I, i++) {
		wsVec Ra_v = *I;
		if (Xp_rowNames)
			Cp_mat->fmt("%s	", (*Xp_rowNames)[i].c_str());
		LOOP (j, N_col) Cp_mat->fmt(S_fmt, Ra_v[j]);
		Cp_mat->put("\n");
	}

	delete Cp_mat;
}

void exportSpsMat(const char *S_ext, wsReal *Rp_mat, wsUint N_r, wsUint N_c,
	wsUint* Na_rEnd, wsUint* Na_cIdx, vStr *Xp_header/*=NULL*/,
	wsUint N_prec/*=0*/)
{
	wsUint		i;
	char		S_fmt[32];
	/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	if (N_prec == 0)
		sprintf(S_fmt, "%%g	");
	else
		sprintf(S_fmt, "%%.%de	", N_prec);

	LOGoutput("Export sparse matrix [%s.%s] [%d * %d]\n", OPT_STRING(out), S_ext, N_r, N_c);

	/* Print header */
	if (Xp_header) {
		for (i = 0; i < N_c; i++)
			Cp_mat->fmt("%s	", (*Xp_header)[i].c_str());
		Cp_mat->put("\n");
	}

	/* Print values along with R start point */
	i = 0;
	wsUint	N_s = 0;
	wsUint	N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];
		/* No element in the row */
		if (N_s == N_e) for (wsUint j = 0; j < N_c; j++)
			Cp_mat->put("0	");

		/* From first column to last column */
		wsUint k = 0;
		for (wsUint j = N_s; j < N_e; j++) {
			/* Front fill */
			for (; k < Na_cIdx[j]; k++) Cp_mat->put("0	");
			Cp_mat->fmt(S_fmt, Rp_mat[j]);
			k++;
		}
		/* Final fill */
		for (; k < N_c; k++) Cp_mat->put("0	");

		/* Do line feed */
		Cp_mat->put("\n");

		/* Update start point */
		i++;
		N_s = N_e;
	} while (i < N_r);

	delete Cp_mat;
}

void exportSymMat(const char *S_ext, wsReal **Rp_mat, wsUint N_sz,
	vStr *Xp_header/*=NULL*/, vStr *Xp_rowNames/*=NULL*/, wsUint N_prec/*=0*/)
{
	wsUint		i;
	char		S_fmt[32];
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	if (N_prec == 0)
		sprintf(S_fmt, "%%g	");
	else
		sprintf(S_fmt, "%%.%de	", N_prec);

	LOGoutput("Export symmetric matrix [%s.%s] [%d]\n", OPT_STRING(out), S_ext, N_sz);
	if (Xp_header) {
		for (i = 0; i < N_sz; i++)
			Cp_mat->fmt("%s	", (*Xp_header)[i].c_str());
		Cp_mat->put("\n");
	}

	for (i = 0; i < N_sz; i++) {
		if (Xp_rowNames) Cp_mat->fmt("%s ", (*Xp_rowNames)[i].c_str());
		for (wsUint j = 0; j <= i; j++)
			Cp_mat->fmt(S_fmt, Rp_mat[i][j]);
		for (wsUint j = i + 1; j < N_sz; j++)
			Cp_mat->fmt(S_fmt, Rp_mat[j][i]);
		Cp_mat->put("\n");
	}
	delete Cp_mat;
}

void exportBlkMat(const char *S_ext, wsReal **Rp_mat, wsUint N_szBlk,
	wsUint N_rep, vStr *Xa_colNames, wsUint N_prec/*=0*/)
{
	LOG("colNames will be ignored\n"); /* FIXME */
	exportBlkMat(S_ext, Rp_mat, N_szBlk, N_rep, N_prec);
}

void exportBlkMat(const char *S_ext, wsReal **Rp_mat, wsUint N_szBlk,
	wsUint N_rep, wsUint N_prec/*=0*/)
{
	wsUint		i, j, k, l;
	char		S_fmt[32];
	/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	if (N_prec == 0)
		sprintf(S_fmt, "%%g	");
	else
		sprintf(S_fmt, "%%.%de	", N_prec);

	LOGoutput("Export block matrix [%s.%s] [%d*%d * %d*%d]\n", OPT_STRING(out),
		S_ext, N_szBlk, N_rep, N_szBlk, N_rep);

	for (i = 0; i < N_rep; i++) {
		for (l = 0; l < N_szBlk; l++) {
			for (j = 0; j < i; j++) for (k = 0; k < N_szBlk; k++)
				Cp_mat->put("0 ");
			for (k = 0; k < N_szBlk; k++)
				Cp_mat->fmt(S_fmt, Rp_mat[l][k]);
			for (j = i + 1; j < N_rep; j++) for (k = 0; k < N_szBlk; k++)
				Cp_mat->put("0 ");
			Cp_mat->put("\n");
		}
	}
	delete Cp_mat;
}

void exportDiagMat(const char *S_ext, wsReal **Rp_mat, wsUint N_sz,
	wsUint N_prec/*=0*/)
{
	wsUint		i;
	char		S_fmt[32];
	/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	if (N_prec == 0)
		sprintf(S_fmt, "%%g	");
	else
		sprintf(S_fmt, "%%.%de	", N_prec);

	for (i = 0; i < N_sz; i++) {
		for (wsUint j = 0; j < N_sz; j++) {
			if (i == j)
				Cp_mat->fmt(S_fmt, Rp_mat[0][i]);
			else
				Cp_mat->put("0 ");
		}
		Cp_mat->put("\n");
	}
	delete Cp_mat;
}

void exportDiagMat(const char *S_ext, wsReal **Rp_mat, wsUint N_sz,
	vStr *Sa_colNames, wsUint N_prec/*=0*/)
{
	wsUint		i;
	char		S_fmt[32];
	/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */
	if (Sa_colNames && Sa_colNames->size() != N_sz)
		halt("Column names given but length[%d] is not expected[%d]", Sa_colNames->size(), N_sz);

	if (N_prec == 0)
		sprintf(S_fmt, "%%g	");
	else
		sprintf(S_fmt, "%%.%de	", N_prec);

	/* Print header if exists */
	if (Sa_colNames) {
		FOREACH(vStr_it, (*Sa_colNames), it)
			Cp_mat->fmt("%s	", it->c_str());
		Cp_mat->put("\n");
	}

	for (i = 0; i < N_sz; i++) {
		for (wsUint j = 0; j < N_sz; j++) {
			if (i == j)
				Cp_mat->fmt(S_fmt, Rp_mat[0][i]);
			else
				Cp_mat->put("0 ");
		}
		Cp_mat->put("\n");
	}
	delete Cp_mat;
}

void exportIdtMat(const char *S_ext, wsUint N_sz, vStr *Sa_colNames/*=NULL*/, vStr *Sa_rowNames/*=NULL*/)
{
	wsUint		i;
	/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */
	if (Sa_colNames && Sa_colNames->size() != N_sz)
		halt("Column names given but length[%d] is not expected[%d]", Sa_colNames->size(), N_sz);

	/* Print header if exists */
	if (Sa_colNames) {
		FOREACH(vStr_it, (*Sa_colNames), it)
			Cp_mat->fmt("%s	", it->c_str());
		Cp_mat->put("\n");
	}

	for (i = 0; i < N_sz; i++) {
		if (Sa_rowNames) Cp_mat->fmt("%s ", (*Sa_rowNames)[i].c_str());
		for (wsUint j = 0; j < N_sz; j++) {
			if (i == j)
				Cp_mat->put("1 ");
			else
				Cp_mat->put("0 ");
		}
		Cp_mat->put("\n");
	}
	delete Cp_mat;
}

void exportValMat(wsReal R_val, const char *S_ext, wsUint N_sz, vStr *Sa_colNames/*=NULL*/)
{
	wsUint		i;
	/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */
	if (Sa_colNames && Sa_colNames->size() != N_sz)
		halt("Column names given but length[%d] is not expected[%d]", Sa_colNames->size(), N_sz);

	/* Print header if exists */
	if (Sa_colNames) {
		FOREACH(vStr_it, (*Sa_colNames), it)
			Cp_mat->fmt("%s	", it->c_str());
		Cp_mat->put("\n");
	}

	for (i = 0; i < N_sz; i++) {
		for (wsUint j = 0; j < N_sz; j++)
			Cp_mat->fmt("%g ", R_val);
		Cp_mat->put("\n");
	}
	delete Cp_mat;
}

/*
 *
 * cBaseMatrix definition
 * 
 */

cMatrix::cMatrix()
{
	Ra_data = NULL;
	N_row = N_col = 0;
}

cMatrix::~cMatrix()
{
}

wsReal** cMatrix::get() const
{
	return Ra_data;
}

wsUintCst cMatrix::row() const
{
	return N_row;
}

wsUintCst cMatrix::col() const
{
	return N_col;
}

const char cMatrix::sqr() const
{
	return row()==col();
}

const char cMatrix::sane()
{
	return (row() > 0) && (col() > 0);
}

void cMatrix::fmtfile(wsStrCst S_fmtExt, ...)
{
	va_list	H_varList;
	char	S_buf[256];
	va_start(H_varList, S_fmtExt);
	vsprintf(S_buf, S_fmtExt, H_varList);
	file(S_buf);
}

void cMatrix::raise(xMatErr X_err)
{
	switch (X_err) {
	case MATERR_DIM_INCORRECT:
		halt("Incorrect dimension");
		break;
	case MATERR_TYPE_INCORRECT:
		halt("Incorrect matrix type");
		break;
	default:
		halt("Invalid error code");
		break;
	}
}

void cMatrix::__init(wsMat Ra_data, wsUint N_row, wsUint N_col,
	char B_passNull/*=0*/)
{
	if (!Ra_data && !B_passNull) halt_fmt(WISARD_SYST_NULL_DATA,
		"Matrix data");
	this->Ra_data	= Ra_data;
	this->N_row		= N_row;
	this->N_col		= N_col;
}

#if 0
cMatrix& cMatrix::crop(UINT_t N_rowStart, UINT_t N_colStart, int N_rowCount, int N_colCount)
{
	REAL_t **Ra_ret = NULL;
	UINT_t N_realRowCount;

	/* Allocate default value */
	if (N_rowCount == -1)
		N_realRowCount = N_row - N_rowStart;
	else
		N_realRowCount = N_rowCount;

	if (N_colCount == -1)
		N_colCount = N_col - N_colStart;
	if (N_rowCount <= 0 || N_colCount <= 0 ||
		N_rowStart+N_rowCount >= N_row || N_colStart+N_colCount >= N_col)
		halt("Invalid range given");

	/* Fill matrix */
	MULTI_MALLOC(Ra_ret, REAL_t*, N_rowCount);
	for (UINT_t i=0,j=N_rowStart ; i<N_realRowCount ; i++,j++) {
		sseMalloc(Ra_ret[i], REAL_t, N_colCount);
		memcpy(Ra_ret[i], Ra_data[j]+N_colStart, sizeof(REAL_t)*N_colCount);
	}

	return *new cMatrix(N_rowCount, N_colCount, Ra_ret);
}

REAL_t cMatrix::sumDiag()
{
	UINT_t	N_i = N_row > N_col ? N_col : N_row;
	REAL_t			R_sum = 0.0f;

	for (UINT_t i=0 ; i<N_i ; i++)
		R_sum += Ra_data[i][i];

	return R_sum;
}

REAL_t cMatrix::sumDiag(cMatrix& C_mat)
{
	if (C_mat.row() != N_col)
		halt("Can't product two matrix mat1col[%d] != mat2row[%d]", N_col, C_mat.row());

	UINT_t	N_col2		=  C_mat.col();
	register REAL_t	R_ret		= 0.0f;
	REAL_t			**Rp_data	=  C_mat.getData();
	UINT_t	N_loop		=  N_row > N_col2 ? N_col2 : N_row;

	/* Blocked manner, non-SSE */
	for (UINT_t i=0 ; i<N_loop ; i++)
		for (UINT_t j=0 ; j<N_col ; j++)
			R_ret += Ra_data[i][j] * Rp_data[j][i];

	return R_ret;
}

REAL_t cMatrix::sumDiag(cMatrix& C_mat1, cMatrix& C_mat2)
{
	if (C_mat1.row() != N_col)
		halt("Can't product two matrix mat1col[%d] != mat2row[%d]", N_col, C_mat1.row());
	if (C_mat2.row() != C_mat1.col())
		halt("Can't product two matrix mat2col[%d] != mat3row[%d]", C_mat1.col(), C_mat2.row());

	REAL_t	**Rp_data1 = C_mat1.getData();
	REAL_t	**Rp_data2 = C_mat2.getData();
	UINT_t	N_col2 = C_mat1.col();
	UINT_t	N_col3 = C_mat2.col();
	UINT_t	N_i = N_row > N_col3 ? N_col3 : N_row;

	/* p = N_mat2Col
	 * m = N_mat1Col
	 * (ABC)_ii=Σ_(k=1)^p (Σ_(j=1)^m A_ij B_jk ) C_ki
	 */
	REAL_t R_ret = 0.0f;
	
	for (UINT_t i=0 ; i<N_i ; i++) {
		for (UINT_t k=0 ; k<N_col2 ; k++) {
			REAL_t R_subSum = 0.0f;

			for (UINT_t j=0 ; j<N_col ; j++)
				R_subSum += Ra_data[i][j] * Rp_data1[j][k];

			R_ret += R_subSum * Rp_data2[k][i];
		}
	}


	return R_ret;
}

cMatrix& cMatrix::genInv(bool& flag)
{
	const double eps = 1e-24;
	REAL_t **u = getData();

	if (N_row == 0) 
		halt("Internal problem: matrix with no rows (inverse function)");
	if (N_row != N_col) 
		halt("Internal problem: Cannot invert non-square matrix");
	int n = N_row;

	vector<double> w(n,0);

	vector<vector<double> > v(n);
	for (int i=0; i<n; i++) 
		v[i].resize(n,0);

	flag = svdcmp(u,N_row,N_col,w,v); 

	// Look for singular values
	double wmax = 0;
	for (int i=0; i<n; i++)
		wmax = w[i] > wmax ? w[i] : wmax;
	double wmin = wmax * eps;
	for (int i=0; i<n; i++)
		w[i] = w[i] < wmin ? 0 : 1/w[i];

	// u w t(v)

	// row U * 1/w

	// results matrix
	cMatrix *Cp_ret;
	MULTI_MALLOC(Cp_ret, cMatrix, 1);
	Cp_ret->init(n, n);
	REAL_t **r = Cp_ret->getData();
	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			u[i][j] = (REAL_t)(u[i][j] * w[j]);

	// [nxn].[t(v)] 
	for (int i=0; i<n; i++)
		for (int j=0; j<n; j++)
			for (int k=0; k<n; k++)
				r[i][j] += (REAL_t)(u[i][k] * v[j][k]);

	return *Cp_ret;
}

cMatrix& cMatrix::product(cMatrix& C_mat1, cMatrix& C_mat2)
{
	if (C_mat1.row() != N_col)
		halt("Can't product two matrix mat1col[%d] != mat2row[%d]", N_col, C_mat1.row());
	if (C_mat2.row() != C_mat1.col())
		halt("Can't product two matrix mat2col[%d] != mat3row[%d]", C_mat1.col(), C_mat2.row());

	UINT_t	N_col2		= C_mat1.col();
	UINT_t	N_col3		= C_mat2.col();
	REAL_t	**Ra_ret	= sseEmptyMatrix(N_row, N_col3);
	REAL_t	**Rp_data1	= C_mat1.getData();
	REAL_t	**Rp_data2	= C_mat2.getData();

	REAL_t	*Ra_inter;
	sseMalloc(Ra_inter, REAL_t, N_col2);

#ifdef USE_SSE
	for (UINT_t i=0 ; i<N_row ; i++) {
		REAL_t	*Rp_matR = Ra_ret[i];
		memset(Ra_inter, 0x00, sizeof(REAL_t)*N_col2);

		/* Calculate intermediate */
		for (UINT_t j=0 ; j<N_col ; j++) {
			sse_t A = sseSet(Ra_data[i][j]);

			for (UINT_t k=0 ; k<N_col2 ; k+=sseJmp) {
				sse_t *C = (sse_t *)(&Ra_inter[k]);
				sse_t *B = (sse_t *)(&Rp_data1[j][k]);

				*C = sseAdd(*C, sseMul(A, *B));
			}
		}

		for (UINT_t j=0 ; j<N_col2 ; j++) {
			sse_t A = sseSet(Ra_inter[j]);

			for (UINT_t l=0 ; l<N_col3 ; l+=sseJmp) {
				sse_t *B = (sse_t *)(&Rp_data2[j][l]);
				sse_t *C = (sse_t *)(&Rp_matR[l]);

				*C = sseAdd(*C, sseMul(A, *B));
			}
		}
	}
#else
	for (UINT_t i=0 ; i<N_row ; i++) {
		REAL_t	*Rp_matR = Ra_ret[i];
		memset(Ra_inter, 0x00, sizeof(REAL_t)*N_col2);

		/* Calculate intermediate */
		for (UINT_t j=0 ; j<N_col ; j++) {
			REAL_t R_elem = Ra_data[i][j];

			for (UINT_t k=0 ; k<N_col2 ; k++)
				Ra_inter[k] += R_elem * Rp_data1[j][k];
		}

		for (UINT_t l=0 ; l<N_col3 ; l++) {
			REAL_t R_sum = 0.0f;

			for (UINT_t j=0 ; j<N_col2 ; j++)
				R_sum += Ra_inter[j] * Rp_data2[j][l];

			Rp_matR[l] = R_sum;
		}
	}
#endif

	sseFree(Ra_inter);

	return *new cMatrix(N_row, N_col3, Ra_ret);
}

void cMatrix::print()
{
	for (wsUint i=0 ; i<N_row ; i++) {
		for (wsUint j=0 ; j<N_col ; j++) {
			printf("%-10g	", Ra_data[i][j]);
		}
		printf("\n");
	}
}

void cMatrix::file(const char *S_ext)
{
	char S_fn[MAX_PATH];
	sprintf(S_fn, "%s.%s", OPT_STRING(out), S_ext);
	FILE *H_fp = fopen(S_fn, "w+");
	if (H_fp == NULL)
		halt("Can't open output file for matrix `%s`", S_fn);

	for (wsUint i=0 ; i<N_row ; i++) {
		for (wsUint j=0 ; j<N_col ; j++) {
			fprintf(H_fp, "%-10g	", Ra_data[i][j]);
		}
		fprintf(H_fp, "\n");
	}
	fclose(H_fp);
}

cMatrix& cMatrix::sweepOp(wsUint N_start, int N_end)
{
	wsUint i, j, k;
	wsUint N_actualEnd;
	wsReal D, sD, sK;
	wsReal *Ra_row;

	if (N_end == -1)
		N_actualEnd = N_row-1;
	else
		N_actualEnd = N_end;

	/* Make duplication of input matrix */
	wsReal **Ra_inv = NULL;
	wsAlloc(Ra_inv, wsReal*, N_row);
	sseMalloc(Ra_row, wsReal, N_row);

	for (i=0 ; i<N_row ; i++) {
		sseMalloc(Ra_inv[i], wsReal, N_row);
		memcpy(Ra_inv[i], Ra_data[i], sizeof(wsReal)*N_row);
	}

	/* When K == 1 */
	if (N_start == 0) {
		// # i : 1
		D = Ra_inv[0][0];
		if (D == 0.0f)
			LOG("This matrix is singular\n");

#ifdef USE_SSE
		sse_t sse_D = sseSet(D);
#endif

		Ra_inv[0][0] = 1.0f/D; // case 1
		memcpy(Ra_row, Ra_inv[0], sizeof(wsReal)*N_row);

#ifdef USE_SSE
		//		printf("[1 ~ 4) ");
		for (j=1 ; j<4 ; j++)
			Ra_inv[0][j] /= D; // case 2
		//		printf("[4 ~ %d)\n", N_row);
		for (j=4 ; j<N_row ; j+=sseJmp) {
			sse_t *sse_pInv = (sse_t *)(Ra_inv[0] + j);
			*sse_pInv = sseDiv(*sse_pInv, sse_D);
		}
#else
		for (j=1 ; j<N_row ; j++)
			Ra_inv[0][j] /= D; // case 2
#endif

		// # i : (1,N]
		for (i=1 ; i<N_row ; i++) {
			sD = Ra_inv[i][0];
#ifdef USE_SSE
			sse_t sse_sD = sseSet(sD);
#endif
			Ra_inv[i][0] = -sD / D; // case 3
#ifdef USE_SSE
			//			printf("[1 ~ 4) ");
			for (j=1 ; j<4 ; j++)
				Ra_inv[i][j] -= sD*Ra_row[j]/D; // case 4
			//			printf("[4 ~ %d)\n", N_row);
			for (j=4 ; j<N_row ; j+=sseJmp) {
				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *vT = (sse_t *)(Ra_inv[i] + j);
				*vT = sseSub(*vT, sseDiv(sseMul(sse_sD, *sse_pRow), sse_D));
			}
#else
			for (j=1 ; j<N_row ; j++)
				Ra_inv[i][j] -= sD*Ra_row[j]/D; // case 4

			Ra_inv[i][j] = Ra_inv[i][j] - sD*Ra_row[j]/D; // case 4
#endif
		}
	}

	/* When K == (1,N) */
	wsUint N_kEnd = N_actualEnd<(N_row-1) ? N_actualEnd : N_row-1;
	for (k=N_start<2?1:N_start ; k<N_kEnd ; k++) {
		memcpy(Ra_row, Ra_inv[k], sizeof(wsReal)*N_row);
		D = Ra_inv[k][k];
		if (D == 0.0f)
			LOG("This matrix is singular\n");
		//		LOG("DO k = %d\n", k+1);

#ifdef USE_SSE
		wsUint adjK_front = getMed(k);
		wsUint adjK_back = getMed(k+1+3);
		sse_t sse_D = sseSet(D);
#endif

		// # i : [1, k)
		for (i=0 ; i<k ; i++) {
			sK = Ra_inv[i][k];
#ifdef USE_SSE
			sse_t sse_sK = sseSet(sK);
#endif

#ifdef USE_SSE
			//			printf("[0 ~ %d) ", adjK);
			for (j=0 ; j<adjK_front ; j+=sseJmp) {
				//				wsReal res[4];
				//				for (MYuint x=0 ; x<4 ; x++)
				//					res[x] = Ra_inv[i][j+x] - sK*Ra_row[j+x]/D;

				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

				*sse_pInv = sseSub(*sse_pInv, sseDiv(sseMul(sse_sK, *sse_pRow), sse_D));
				//				for (MYuint x=0 ; x<4 ; x++)
				//					if ((j+x) < N_row && (Ra_inv[i][j+x]-res[x])>0.0001)
				//						halt("%g", Ra_inv[i][j+x]-res[x]);
			}
			//			printf("[%d ~ %d) ", adjK, k);
			for (j=adjK_front ; j<k ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#else
			for (j=0 ; j<k ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#endif

			//			printf("[%d] ", k);
			Ra_inv[i][k] = -sK / D; // case 3
#ifdef USE_SSE
			//			printf("[%d~%d) ", k+1, adjK);
			for (j=k+1 ; j<adjK_back ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
			//			printf("[%d~%d)\n", adjK, N_row);
			for (j=adjK_back ; j<N_row ; j+=sseJmp) {
				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

				*sse_pInv = sseSub(*sse_pInv, sseDiv(sseMul(sse_sK, *sse_pRow), sse_D));
			}
#else
			for (j=k+1 ; j<N_row ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#endif
		}

		// # i == k
#ifdef USE_SSE
		//			printf("[0 ~ %d) ", adjK);
		for (j=0 ; j<adjK_front ; j+=sseJmp) {
			//				wsReal res[4];
			//				for (MYuint x=0 ; x<4 ; x++)
			//					res[x] = Ra_inv[i][j+x] - sK*Ra_row[j+x]/D;

			sse_t *sse_pRow = (sse_t *)(Ra_row + j);
			sse_t *sse_pInv = (sse_t *)(Ra_inv[k] + j);

			*sse_pInv = sseDiv(*sse_pRow, sse_D);
			//				for (MYuint x=0 ; x<4 ; x++)
			//					if ((j+x) < N_row && (Ra_inv[i][j+x]-res[x])>0.0001)
			//						halt("%g", Ra_inv[i][j+x]-res[x]);
		}
		//			printf("[%d ~ %d) ", adjK, k);
		for (j=adjK_front ; j<k ; j++)
			Ra_inv[k][j] = Ra_row[j] / D; // case 2
#else
		for (j=0 ; j<k ; j++)
			Ra_inv[k][j] = Ra_row[j] / D; // case 2
#endif
		Ra_inv[k][k] = 1.0f/D; // case 1
#ifdef USE_SSE
		//			printf("[%d~%d) ", k+1, adjK);
		for (j=k+1 ; j<adjK_back ; j++)
			Ra_inv[i][j] = Ra_row[j] / D; // case 2
		//			printf("[%d~%d)\n", adjK, N_row);
		for (j=adjK_back ; j<N_row ; j+=sseJmp) {
			sse_t *sse_pRow = (sse_t *)(Ra_row + j);
			sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

			*sse_pInv = sseDiv(*sse_pRow, sse_D);
		}
#else
		for (j=k+1 ; j<N_row ; j++)
			Ra_inv[k][j] = Ra_row[j] / D; // case 2
#endif

		// # i == (k, N)
		for (i=k+1 ; i<N_row ; i++) {
			sK = Ra_inv[i][k];
#ifdef USE_SSE
			sse_t sse_sK = sseSet(sK);
#endif

#ifdef USE_SSE
			//			printf("[0 ~ %d) ", adjK);
			for (j=0 ; j<adjK_front ; j+=sseJmp) {
				//				wsReal res[4];
				//				for (MYuint x=0 ; x<4 ; x++)
				//					res[x] = Ra_inv[i][j+x] - sK*Ra_row[j+x]/D;

				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

				*sse_pInv = sseSub(*sse_pInv, sseDiv(sseMul(sse_sK, *sse_pRow), sse_D));
				//				for (MYuint x=0 ; x<4 ; x++)
				//					if ((j+x) < N_row && (Ra_inv[i][j+x]-res[x])>0.0001)
				//						halt("%g", Ra_inv[i][j+x]-res[x]);
			}
			//			printf("[%d ~ %d) ", adjK, k);
			for (j=adjK_front ; j<k ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#else
			for (j=0 ; j<k ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#endif
			Ra_inv[i][k] = -sK / D; // case 3
#ifdef USE_SSE
			//			printf("[%d~%d) ", k+1, adjK);
			for (j=k+1 ; j<adjK_back ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
			//			printf("[%d~%d)\n", adjK, N_row);
			for (j=adjK_back ; j<N_row ; j+=sseJmp) {
				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

				*sse_pInv = sseSub(*sse_pInv, sseDiv(sseMul(sse_sK, *sse_pRow), sse_D));
			}
#else
			for (j=k+1 ; j<N_row ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#endif
		}
	}

	/* When K == N */
	if (N_actualEnd == (N_row-1)) {
		// # i == [1,N)
		wsUint N = N_row-1;
		//		LOG("DO k = E\n");
		memcpy(Ra_row, Ra_inv[N], sizeof(wsReal)*N_row);
		D = Ra_inv[N][N];
		if (D == 0.0f)
			LOG("This matrix is singular\n");

#ifdef USE_SSE
		sse_t sse_D = sseSet(D);
#endif
		for (i=0 ; i<N ; i++) {
#ifdef USE_SSE
			sse_t sse_N = sseSet(Ra_inv[i][N]);
			wsUint adjN = getMed(N);
			for (j=0 ; j<adjN ; j+=sseJmp) {
				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

				*sse_pInv = sseSub(*sse_pInv, sseDiv(sseMul(sse_N, *sse_pRow), sse_D));
			}
			for (j=adjN ; j<N ; j++)
				Ra_inv[i][j] -= Ra_inv[i][N]*Ra_row[j]/D; // case 4
#else
			for (j=0 ; j<N ; j++)
				Ra_inv[i][j] -= Ra_inv[i][N]*Ra_row[j]/D; // case 4
#endif
			Ra_inv[i][N] = -Ra_inv[i][N] / D; // case 3
		}
		// # i == N
		for (j=0 ; j<N ; j++)
			Ra_inv[N][j] = Ra_row[j] / D; // case 2
		Ra_inv[N][N] = 1/D; // case 1
	}

	return *new cMatrix(N_row, N_col, Ra_inv);
}

#endif

// Do mt = m+s
void sseMaC(wsMat Ra_m, wsUintCst N_r, wsUintCst N_c, wsRealCst R_v,
			wsMat Ra_t, wsUintCst N_y/*=0*/, wsUintCst N_x/*=0*/)
{
	wsUint	i, j;
	sse_t	sse_s = sseSet(R_v);

	if (Ra_t == NULL) {
		if (N_y || N_x) halt("ERR");

		for (i=0 ; i<N_r ; i++) {
			wsVec	Ra_1	= Ra_m[i];
			wsUint	N_med	= 0;
#ifdef USE_SSE
			N_med = getMed(N_c);
			for (j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_1 = (sse_t *)(Ra_1 + j);
				*sse_1 = sseAdd(*sse_1, sse_s);
			}
#endif
			/* Do rest part */
			for (j=N_med ; j<N_c ; j++)
				Ra_1[j] += R_v;
		}
	} else {
		for (i=0 ; i<N_r ; i++) {
			wsVecCst	Ra_1	= Ra_m[i];
			wsVec	Ra_T	= Ra_t[i+N_y];
			wsUint	N_med	= 0;
#ifdef USE_SSE
			N_med = getMed(N_c);
			for (j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_1 = (sse_t *)(Ra_1 + j);
				sse_t *sse_T = (sse_t *)(Ra_T + j+N_x);
				*sse_T = sseAdd(*sse_1, sse_s);
			}
#endif
			/* Do rest part */
			for (j=N_med ; j<N_c ; j++)
				Ra_T[j] = Ra_1[j] + R_v;
		}
	}
}

void sseMaM(wsReal **Ra_m1, wsReal **Ra_m2, wsReal **Ra_mt,
			wsUint N_row, wsUint N_col/*=0xfffffff*/,
			wsReal R_multiplier/*=REAL_CONST(.0)*/)
{
	sseMaM(Ra_m1, Ra_m2, Ra_mt, 0, 0, N_row, N_col, R_multiplier);
}

void sseMsM(wsReal **Ra_m1, wsReal **Ra_m2, wsReal **Ra_mt,
			wsUint N_row, wsUint N_col/*=0xfffffff*/)
{
	sseMsM(Ra_m1, Ra_m2, Ra_mt, 0, 0, N_row, N_col);
}

void sseMsVtV(wsMat Ra_m1, wsVec Ra_v1, wsVec Ra_v2, wsUint N_row, wsUint N_col)
{
#ifdef USE_SSE
	wsUint N_med = getMed(N_col);
	LOOP (i, N_row) {
		wsVec Ra_1 = Ra_m1[i];
		wsUint j=0;
		sse_t sse_v1i = sseSet(Ra_v1[i]);
		for ( ; j<N_med ; j+=sseJmp) {
			sse_t *sse_mij = (sse_t *)(Ra_1 + j);
			sse_t *sse_v2j = (sse_t *)(Ra_v2 + j);

			*sse_mij = sseSub(*sse_mij, sseMul(sse_v1i, *sse_v2j));
		}
		for (; j<N_col ; j++)
			Ra_1[j] -= Ra_v1[i] * Ra_v2[j];
	}
#else
	LOOP (i, N_row) LOOP (j, N_col)
		Ra_m1[i][j] -= Ra_v1[i] * Ra_v2[j];
#endif
}

//void sseMsMsq(wsReal **Ra_m1, wsReal **Ra_m2, wsReal **Ra_mt,
//			wsUint N_row, wsUint N_col/*=0xfffffff*/)
//{
//	sseMsMsq(Ra_m1, Ra_m2, Ra_mt, 0, 0, N_row, N_col);
//}

void sseMaM(wsReal **Ra_m1, wsReal **Ra_m2, wsReal **Ra_mt, wsUint N_y,
			wsUint N_x, wsUint N_row, wsUint N_col/*=0xfffffff*/,
			wsReal R_multiplier/*=REAL_CONST(.0)*/)
{
	/* Default value correction */
	if (N_col == 0xffffffff) N_col = N_row;

	/* If case m1 -= m2 */
	if (R_multiplier == W0) {
		if (Ra_m1 == Ra_mt) {
			for (wsUint i=0 ; i<N_row ; i++) {
				wsUint	N_med	= 0;
				wsReal	*Ra_1	= Ra_m1[i + N_y];
				wsReal	*Ra_2	= Ra_m2[i];

#ifdef USE_SSE
				N_med = getMed(N_col);
				for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
					sse_t *sse_m1 = (sse_t *)(Ra_1 + j+N_x);
					sse_t *sse_m2 = (sse_t *)(Ra_2 + j);

					*sse_m1 = sseAdd(*sse_m1, *sse_m2);
				}
#endif
				/* Do rest part */
				for (wsUint j=N_med ; j<N_col ; j++)
					Ra_1[j+N_x] += Ra_m2[i][j];
			}
		} else for (wsUint i=0 ; i<N_row ; i++) {
			wsUint	N_med	= 0;
			wsReal	*Ra_1	= Ra_m1[i];
			wsReal	*Ra_2	= Ra_m2[i];
			wsReal	*Ra_T	= Ra_mt[i + N_y];

#ifdef USE_SSE
			N_med = getMed(N_col);
			for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_m1 = (sse_t *)(Ra_1 + j);
				sse_t *sse_m2 = (sse_t *)(Ra_2 + j);
				sse_t *sse_mt = (sse_t *)(Ra_T + j+N_x);

				*sse_mt = sseAdd(*sse_m1, *sse_m2);
			}
#endif
			/* Do rest part */
			for (wsUint j=N_med ; j<N_col ; j++)
				Ra_T[j+N_x] = Ra_1[j] + Ra_2[j];
		}
	} else {
		/* In case of have multiplier */
#ifdef USE_SSE
		sse_t sse_MUL = sseSet(R_multiplier);
#endif
		if (Ra_m1 == Ra_mt) {
			for (wsUint i=0 ; i<N_row ; i++) {
				wsUint	N_med	= 0;
				wsReal	*Ra_1	= Ra_m1[i + N_y];
				wsReal	*Ra_2	= Ra_m2[i];

#ifdef USE_SSE
				N_med = getMed(N_col);
				for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
					sse_t *sse_m1 = (sse_t *)(Ra_1 + j+N_x);
					sse_t *sse_m2 = (sse_t *)(Ra_2 + j);

					*sse_m1 = sseMul(sseAdd(*sse_m1, *sse_m2), sse_MUL);
				}
#endif
				/* Do rest part */
				for (wsUint j=N_med ; j<N_col ; j++)
					Ra_1[j+N_x] = (Ra_1[j+N_x]+Ra_m2[i][j]) * R_multiplier;
			}
		} else for (wsUint i=0 ; i<N_row ; i++) {
			wsUint	N_med	= 0;
			wsReal	*Ra_1	= Ra_m1[i];
			wsReal	*Ra_2	= Ra_m2[i];
			wsReal	*Ra_T	= Ra_mt[i + N_y];

#ifdef USE_SSE
			N_med = getMed(N_col);
			for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_m1 = (sse_t *)(Ra_1 + j);
				sse_t *sse_m2 = (sse_t *)(Ra_2 + j);
				sse_t *sse_mt = (sse_t *)(Ra_T + j+N_x);

				*sse_mt = sseMul(sseAdd(*sse_m1, *sse_m2), sse_MUL);
			}
#endif
			/* Do rest part */
			for (wsUint j=N_med ; j<N_col ; j++)
				Ra_T[j+N_x] = (Ra_1[j] + Ra_2[j]) * R_multiplier;
		}
	}
}

void sseMaVtV(wsMat Ra_m1, wsVec Ra_v1, wsVec Ra_v2, wsUint N_row, wsUint N_col)
{
#ifdef USE_SSE
	wsUint N_med = getMed(N_col);
	LOOP(i, N_row) {
		wsVec Ra_1 = Ra_m1[i];
		wsUint j=0;
		sse_t sse_v1i = sseSet(Ra_v1[i]);
		for (; j < N_med ; j+=sseJmp) {
			sse_t *sse_mij = (sse_t *)(Ra_1 + j);
			sse_t *sse_v2j = (sse_t *)(Ra_v2 + j);

			*sse_mij = sseAdd(*sse_mij, sseMul(sse_v1i, *sse_v2j));
		}
		for (; j < N_col ; j++)
			Ra_1[j] += Ra_v1[i] * Ra_v2[j];
	}
#else
	LOOP(i, N_row) LOOP(j, N_col)
		Ra_m1[i][j] += Ra_v1[i] * Ra_v2[j];
#endif
}

void sseMsM(wsReal **Ra_m1, wsReal **Ra_m2, wsReal **Ra_mt, wsUint N_y,
			wsUint N_x, wsUint N_row, wsUint N_col/*=0xfffffff*/)
{
	/* Default value correction */
	if (N_col == 0xffffffff) N_col = N_row;

	/* If case m1 -= m2 */
	if (Ra_m1 == Ra_mt) for (wsUint i=0 ; i<N_row ; i++) {
		wsUint	N_med	= 0;
		wsReal	*Ra_1	= Ra_m1[i + N_y];
		wsReal	*Ra_2	= Ra_m2[i];

#ifdef USE_SSE
		N_med = getMed(N_col);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m1 = (sse_t *)(Ra_1 + j+N_x);
			sse_t *sse_m2 = (sse_t *)(Ra_2 + j);

			*sse_m1 = sseSub(*sse_m1, *sse_m2);
		}
#endif
		/* Do rest part */
		for (wsUint j=N_med ; j<N_col ; j++)
			Ra_1[j+N_x] -= Ra_m2[i][j];
	} else for (wsUint i=0 ; i<N_row ; i++) {
		wsUint	N_med	= 0;
		wsReal	*Ra_1	= Ra_m1[i];
		wsReal	*Ra_2	= Ra_m2[i];
		wsReal	*Ra_T	= Ra_mt[i + N_y];

#ifdef USE_SSE
		N_med = getMed(N_col);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m1 = (sse_t *)(Ra_1 + j);
			sse_t *sse_m2 = (sse_t *)(Ra_2 + j);
			sse_t *sse_mt = (sse_t *)(Ra_T + j+N_x);

			*sse_mt = sseSub(*sse_m1, *sse_m2);
		}
#endif
		/* Do rest part */
		for (wsUint j=N_med ; j<N_col ; j++)
			Ra_T[j+N_x] = Ra_1[j] - Ra_2[j];
	}
}

void sseMsMsq(wsReal **Ra_m1, wsReal **Ra_m2, wsReal **Ra_mt, wsUint N_y,
	wsUint N_x, wsUint N_row, wsUint N_col/*=0xfffffff*/)
{
	/* Default value correction */
	if (N_col == 0xffffffff) N_col = N_row;

	/* If case m1 -= m2 */
	if (Ra_m1 == Ra_mt) for (wsUint i=0 ; i<N_row ; i++) {
		wsUint	N_med	= 0;
		wsReal	*Ra_1	= Ra_m1[i + N_y];
		wsReal	*Ra_2	= Ra_m2[i];

#ifdef USE_SSE
		N_med = getMed(N_col);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m1 = (sse_t *)(Ra_1 + j+N_x);
			sse_t *sse_m2 = (sse_t *)(Ra_2 + j);

			sse_t sse_tmp = sseSub(*sse_m1, *sse_m2);
			*sse_m1 = sseMul(sse_tmp, sse_tmp);
		}
#endif
		/* Do rest part */
		for (wsUint j=N_med ; j<N_col ; j++) {
			Ra_1[j+N_x] -= Ra_m2[i][j];
			Ra_1[j+N_x] = SQR(Ra_1[j+N_x]);
		}
	} else for (wsUint i=0 ; i<N_row ; i++) {
		wsUint	N_med	= 0;
		wsReal	*Ra_1	= Ra_m1[i];
		wsReal	*Ra_2	= Ra_m2[i];
		wsReal	*Ra_T	= Ra_mt[i + N_y];

#ifdef USE_SSE
		N_med = getMed(N_col);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m1 = (sse_t *)(Ra_1 + j);
			sse_t *sse_m2 = (sse_t *)(Ra_2 + j);
			sse_t *sse_mt = (sse_t *)(Ra_T + j+N_x);

			sse_t sse_tmp = sseSub(*sse_m1, *sse_m2);
			*sse_mt = sseMul(sse_tmp, sse_tmp);
		}
#endif
		/* Do rest part */
		for (wsUint j=N_med ; j<N_col ; j++) {
			Ra_T[j+N_x] = Ra_1[j] - Ra_2[j];
			Ra_T[j+N_x] = SQR(Ra_T[j+N_x]);
		}
	}

}

wsSym sseSsq(wsSym Ra_s, wsUint N_sz)
{
	wsSym Ra_ret = sseSymMat(N_sz);
	LOOP (i, N_sz) sseVsquare(Ra_s[i], i+1, Ra_ret[i]);
	return Ra_ret;
}

wsMat sseMsq(wsMat Ra_m, wsUint N_row, wsUint N_col)
{
	wsMat Ra_ret = sseMatrix(N_row, N_col);
	LOOP (i, N_row) sseVsquare(Ra_m[i], N_col, Ra_ret[i]);
	return Ra_ret;
}

void sseMaS(wsMatCst Ra_m, wsSym Ra_s, wsMat Ra_t, wsUintCst N_sz)
{
	for (wsUint i=0 ; i<N_sz ; i++) {
		for (wsUint j=0 ; j<i ; j++) {
			Ra_t[i][j] = Ra_m[i][j] + Ra_s[i][j];
			Ra_t[j][i] = Ra_m[j][i] + Ra_s[i][j];
		}
	}
}

void sseMsS(wsMatCst Ra_m, wsSymCst Ra_s, wsMat Ra_t, wsUintCst N_sz)
{
	for (wsUint i=0 ; i<N_sz ; i++) {
		for (wsUint j=0 ; j<i ; j++) {
			Ra_t[i][j] = Ra_m[i][j] - Ra_s[i][j];
			Ra_t[j][i] = Ra_m[j][i] - Ra_s[i][j];
		}
	}
}

// Do mt = m1-m2 or (m1-m2)*s
void sseMsM(wsReal **Ra_m1, wsReal **Ra_m2, wsReal **Ra_mt, wsUint N_y,
			wsUint N_x, wsUint N_row, wsUint N_col/*=0xfffffff*/,
			wsReal R_multiplier/*=W0*/)
{
	if (N_col == 0xffffffff)
		N_col = N_row;

	if (R_multiplier == W0) {
		for (wsUint i=0 ; i<N_row ; i++) {
			wsReal	*Ra_1	= Ra_m1[i];
			wsReal	*Ra_2	= Ra_m2[i];
			wsReal	*Ra_T	= Ra_mt[i + N_y];
			wsUint	N_med	= 0;
#ifdef USE_SSE
			N_med = getMed(N_col);
			for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_1 = (sse_t *)(Ra_1 + j);
				sse_t *sse_2 = (sse_t *)(Ra_2 + j);
				sse_t *sse_T = (sse_t *)(Ra_T + j+N_x);

				*sse_T = sseSub(*sse_1, *sse_2);
			}
#endif
			/* Do rest part */
			for (wsUint j=N_med ; j<N_col ; j++)
				Ra_mt[i][j] = Ra_m1[i][j]-Ra_m2[i][j];
		}
	} else {
		sse_t	sse_M = sseSet(R_multiplier);
		for (wsUint i=0 ; i<N_row ; i++) {
			wsReal	*Ra_1	= Ra_m1[i];
			wsReal	*Ra_2	= Ra_m2[i];
			wsReal	*Ra_T	= Ra_mt[i + N_y];
			wsUint	N_med	= 0;
#ifdef USE_SSE
			N_med = getMed(N_col);
			for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_m1 = (sse_t *)(Ra_1 + j);
				sse_t *sse_m2 = (sse_t *)(Ra_2 + j);
				sse_t *sse_mt = (sse_t *)(Ra_T + j);

				*sse_mt = sseMul(sseSub(*sse_m1, *sse_m2), sse_M);
			}
#endif
			/* Do rest part */
			for (wsUint j=N_med ; j<N_col ; j++)
				Ra_mt[i][j] = (Ra_m1[i][j]-Ra_m2[i][j])*R_multiplier;
		}
	}
}

void sseSsM(wsSym Ra_m1, wsMat Ra_m2, wsMat Ra_mt, wsUint N_row,
	wsUint N_col/*=0xffffffff*/, wsReal R_multiplier/*=0.0*/)
{
	/* Default value correction */
	if (N_col == 0xffffffff) N_col = N_row;

	if (N_row != N_col) halt("SsM must be square");

	if (R_multiplier != W0) {
		for (wsUint i=0 ; i<N_row ; i++) {
			for (wsUint j=0 ; j<i ; j++) {
				Ra_mt[i][j] = (Ra_m1[i][j] - Ra_m2[i][j])*R_multiplier;
				Ra_mt[j][i] = (Ra_m1[i][j] - Ra_m2[j][i])*R_multiplier;
			}
			Ra_mt[i][i] = (Ra_m1[i][i] - Ra_m2[i][i])*R_multiplier;
		}
	} else {
		for (wsUint i=0 ; i<N_row ; i++) {
			for (wsUint j=0 ; j<i ; j++) {
				Ra_mt[i][j] = Ra_m1[i][j] - Ra_m2[i][j];
				Ra_mt[j][i] = Ra_m1[i][j] - Ra_m2[j][i];
			}
			Ra_mt[i][i] = Ra_m1[i][i] - Ra_m2[i][i];
		}
	}
}

// Do mt = m1-m2 or (m1-m2)*s
void sseMsM(wsReal **Ra_m1, wsReal **Ra_m2, wsReal **Ra_mt, wsUint N_row,
			wsUint N_col, wsReal R_multiplier/*=W0*/)
{
	/* Default value correction */
	if (N_col == 0xffffffff) N_col = N_row;

	if (R_multiplier == W0) {
		for (wsUint i=0 ; i<N_row ; i++) {
			wsUint N_med = 0;
#ifdef USE_SSE
			N_med = getMed(N_col);
			for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_m1 = (sse_t *)(Ra_m1[i] + j);
				sse_t *sse_m2 = (sse_t *)(Ra_m2[i] + j);
				sse_t *sse_mt = (sse_t *)(Ra_mt[i] + j);

				*sse_mt = sseSub(*sse_m1, *sse_m2);
			}
#endif
			/* Do rest part */
			for (wsUint j=N_med ; j<N_col ; j++)
				Ra_mt[i][j] = Ra_m1[i][j]-Ra_m2[i][j];
		}
	} else {
		sse_t	sse_M = sseSet(R_multiplier);
		for (wsUint i=0 ; i<N_row ; i++) {
			wsUint N_med = 0;
#ifdef USE_SSE
			N_med = getMed(N_col);
			for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_m1 = (sse_t *)(Ra_m1[i] + j);
				sse_t *sse_m2 = (sse_t *)(Ra_m2[i] + j);
				sse_t *sse_mt = (sse_t *)(Ra_mt[i] + j);

				*sse_mt = sseMul(sseSub(*sse_m1, *sse_m2), sse_M);
			}
#endif
			/* Do rest part */
			for (wsUint j=N_med ; j<N_col ; j++)
				Ra_mt[i][j] = (Ra_m1[i][j]-Ra_m2[i][j])*R_multiplier;
		}
	}
}

// Do mt = m-s
void sseMsC(wsMatCst Ra_m1, wsRealCst R_s, wsMat Ra_mt, wsUintCst N_row,
			wsUint N_col/*=0xffffffff*/)
{
	if (N_col == 0xffffffff)
		N_col = N_row;
#ifdef USE_SSE
	sse_t	sse_S = sseSet(R_s);
	wsUint	N_med = getMed(N_col);
	for (wsUint i=0 ; i<N_row ; i++) {
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m1 = (sse_t *)(Ra_m1[i] + j);
			sse_t *sse_mt = (sse_t *)(Ra_mt[i] + j);

			*sse_mt = sseSub(*sse_m1, sse_S);
		}
		/* Do rest part */
		for (wsUint j=N_med ; j<N_col ; j++)
			Ra_mt[i][j] = Ra_m1[i][j]-R_s;
	}
#else
	for (wsUint i=0 ; i<N_row ; i++)
		for (wsUint j=0 ; j<N_col ; j++)
			Ra_mt[i][j] = Ra_m1[i][j]-R_s;
#endif
}

// Do mt = (m1-m2)^2
void sseMsMsq(wsReal **Ra_m1, wsReal **Ra_m2, wsReal **Ra_mt, wsUint N_row,
			wsUint N_col)
{
	/* Default value correction */
	if (N_col == 0xffffffff) N_col = N_row;

	for (wsUint i=0 ; i<N_row ; i++) {
		wsUint N_med = 0;
#ifdef USE_SSE
		N_med = getMed(N_col);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m1 = (sse_t *)(Ra_m1[i] + j);
			sse_t *sse_m2 = (sse_t *)(Ra_m2[i] + j);
			sse_t *sse_mt = (sse_t *)(Ra_mt[i] + j);

			sse_t sse_tmp = sseSub(*sse_m1, *sse_m2);
			*sse_mt = sseMul(sse_tmp, sse_tmp);
		}
#endif
		/* Do rest part */
		for (wsUint j=N_med ; j<N_col ; j++) {
			Ra_mt[i][j] = Ra_m1[i][j]-Ra_m2[i][j];
			Ra_mt[i][j] = SQR(Ra_mt[i][j]);
		}
	}
}

// Do mt = (m1-v)^2, len of v must be eq. w/ row of m1
void sseMsVsq(wsReal **Ra_m1, wsReal *Ra_V, wsReal **Ra_mt, wsUint N_row,
			  wsUint N_col)
{
	/* Default value correction */
	if (N_col == 0xffffffff) N_col = N_row;

	for (wsUint i=0 ; i<N_row ; i++) {
		wsUint N_med = 0;
#ifdef USE_SSE
		N_med = getMed(N_col);
		sse_t sse_m2 = sseSet(Ra_V[i]);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m1 = (sse_t *)(Ra_m1[i] + j);
			sse_t *sse_mt = (sse_t *)(Ra_mt[i] + j);

			sse_t sse_tmp = sseSub(*sse_m1, sse_m2);
			*sse_mt = sseMul(sse_tmp, sse_tmp);
		}
#endif
		/* Do rest part */
		for (wsUint j=N_med ; j<N_col ; j++) {
			Ra_mt[i][j] = Ra_m1[i][j]-Ra_V[i];
			Ra_mt[i][j] = SQR(Ra_mt[i][j]);
		}
	}
}

// Do mt = (m1*(1-m1))
void sseMp1p(wsReal **Ra_m1, wsReal **Ra_mt, wsUint N_row, wsUint N_col)
{
	/* Default value correction */
	if (N_col == 0xffffffff) N_col = N_row;

	for (wsUint i=0 ; i<N_row ; i++) {
		wsUint N_med = 0;
#ifdef USE_SSE
		N_med = getMed(N_col);
		sse_t sse_1 = sseSet(1.0);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m1 = (sse_t *)(Ra_m1[i] + j);
			sse_t *sse_mt = (sse_t *)(Ra_mt[i] + j);

			*sse_mt = sseMul(*sse_m1, sseSub(sse_1, *sse_m1));
		}
#endif
		/* Do rest part */
		for (wsUint j=N_med ; j<N_col ; j++)
			Ra_mt[i][j] = Ra_m1[i][j] * (W1 - Ra_m1[i][j]);
	}
}

wsMat sseMpM(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsReal **Rp_rMean/*=NULL*/)
{
	/* Check */
	if (N_c1 != N_r2)
		halt("Can't do matrix multplication, C(mat1)=[%u] != [%u]=R(mat2)",
			N_c1, N_r2);

	/* Make return matrix */
	wsMat Ra_ret	= sseEmptyMatrix(N_r1, N_c2);
	wsUint	i, j, k;
	if (Rp_rMean != NULL) {
		sseCalloc(*Rp_rMean, wsReal, N_r1);
	}

	/* Calc */
	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_R		= Ra_ret[i];
		wsReal	*Rp_1		= Rp_m1[i];

		for (j=0 ; j<N_c1 ; j++) {
			wsReal *Rp_2	= Rp_m2[j];
			wsReal R_1j		= Rp_1[j];
			wsUint N_med	= 0;
#ifdef USE_SSE
			N_med = getMed(N_c2);
			sse_t sse_1ij = sseSet(R_1j);
			/* SSE-enabled part */
			for (k=0 ; k<N_med ; k+=sseJmp) {
				// m2[j,k:(k+3)]
				sse_t *sse_2jk	= (sse_t *)(Rp_2 + k);
				// m3[i,k:(k+3)] (j)
				sse_t *sse_mR	= (sse_t *)(Rp_R + k);

				*sse_mR = sseAdd(*sse_mR, sseMul(sse_1ij, *sse_2jk));
			}
#endif
			/* Rest part */
			for (k=N_med ; k<N_c2 ; k++)
				Rp_R[k] += R_1j * Rp_2[k];
		}
		/* Mean process */
		if (Rp_rMean)
			(*Rp_rMean)[i] = sseVsum(Rp_R, N_c2) / (wsReal)N_c2;
	}
	/* Return */
	return Ra_ret;
}

/* Calculate M' %*% M and returns symmetric matrix
 *
 *
 Essentially equivalent to below R code

 sym.mtm <- function(m) {
 r <- dim(m)[1]
 c <- dim(m)[2]

 ret <- matrix(0, nrow=c, ncol=c)
 for (i in 1:r) {
 for (j in 1:c) {
 ret[j,1:j] <- ret[j,1:j] + m[i,j]*m[i,1:j]
 }
 }
 return(ret)
 }
 m <- matrix(rnorm(10), nrow=2)
 sym.mtm(m)

 */
wsSym sseMtM(wsReal **Ra_m1, wsUint N_r1, wsUint N_c1, wsRealCst R_mul/*=W0*/)
{
	wsUint i, j, k;

	wsSym Ra_ret = sseEmptySymMat(N_c1);

	for (i=0 ; i<N_r1 ; i++) {
		wsReal *Ra_1 = Ra_m1[i];

		for (j=0 ; j<N_c1 ; j++) {
			wsUint N_med	= 0;
			wsReal R_1j		= Ra_1[j];

#ifdef USE_SSE
			sse_t sse_A = sseSet(R_1j);
			N_med = getMed(j+1);
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_R = (sse_t *)(Ra_ret[j] + k);
				sse_t *sse_T = (sse_t *)(Ra_1 + k);
				*sse_R = sseAdd(*sse_R, sseMul(sse_A, *sse_T));
			}
#endif
			for (k=N_med ; k<=j ; k++)
				Ra_ret[j][k] += R_1j*Ra_1[k];
		}
	}

	/* R_mul process */
	if (R_mul != W0) LOOP (i, N_r1)
		sseVpC(Ra_ret[i], R_mul, Ra_ret[i], i+1);

	return Ra_ret;
}

/*
mtm <- function(a, b) {
d1 <- dim(a)
d2 <- dim(b)

ret <- matrix(0, nrow=d1[2], ncol=d2[2])
for (i in 1:d1[1]) {
for (j in 1:d1[2]) {
Sa <- a[i,j]
k <- 1:d2[2]
ret[j,k] <- ret[j,k] + Sa*b[i,k]
}
}
return(ret)
}
*/
wsMat sseMtM(wsReal **Ra_m1, wsUint N_r1, wsUint N_c1, wsReal **Ra_m2,
	wsUint N_r2, wsUint N_c2)
{
	wsUint k;
	/* N_r1 and N_r2 should be equal */
	if (N_r1 != N_r2)
		halt_fmt(WISARD_SYST_INVL_DIM, "matrix multiplication", N_r1, N_r2);

	wsMat Ra_ret = sseEmptyMatrix(N_c1, N_c2);

	LOOP (i, N_r1) LOOP (j, N_c1) {
		wsUint N_med = 0;
#ifdef USE_SSE
		sse_t sse_A = sseSet(Ra_m1[i][j]);
		N_med = getMed(N_c2);
		for (k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *sse_R = (sse_t *)(Ra_ret[j] + k);
			sse_t *sse_T = (sse_t *)(Ra_m2[i] + k);
			*sse_R = sseAdd(*sse_R, sseMul(sse_A, *sse_T));
		}
#endif
		for (k=N_med ; k<N_c2 ; k++)
			Ra_ret[j][k] += Ra_m1[i][j]*Ra_m2[i][k];
	}

	return Ra_ret;
}

wsSym sym_sseMtM(wsReal **Ra_m1, wsUint N_r1, wsUint N_c1, wsReal **Ra_m2,
	wsUint N_r2, wsUint N_c2)
{
	wsUint i, j, k;
	/* N_r1 and N_r2 should be equal */
	if (N_r1 != N_r2)
		halt_fmt(WISARD_SYST_INVL_DIM, "matrix multiplication", N_r1, N_r2);
	if (N_c1 != N_c2)
		halt_fmt(WISARD_SYST_INVL_DIM, "matrix multiplication", N_c1, N_c2);

	wsSym Ra_ret = sseEmptySymMat(N_c1);
	for (i=0 ; i<N_r1 ; i++) {
		for (j=0 ; j<N_c1 ; j++) {
			wsUint N_med = 0;
#ifdef USE_SSE
			sse_t sse_A = sseSet(Ra_m1[i][j]);
			N_med = getMed(j+1);
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_R = (sse_t *)(Ra_ret[j] + k);
				sse_t *sse_T = (sse_t *)(Ra_m2[i] + k);
				*sse_R = sseAdd(*sse_R, sseMul(sse_A, *sse_T));
			}
#endif
			for (k=N_med ; k<=j ; k++)
				Ra_ret[j][k] += Ra_m1[i][j]*Ra_m2[i][k];
		}
	}

	return Ra_ret;
}

wsMat sseMtS(wsMat Ra_m1, wsUint N_r1, wsUint N_c1, wsSym Ra_m2,
			 wsUint N_s2)
{
	wsUint i, j, k;
	/* N_r1 and N_s2 should be equal */
	if (N_r1 != N_s2)
		halt_fmt(WISARD_SYST_INVL_DIM, "matrix multiplication", N_r1, N_s2);

	wsMat	Ra_ret	= sseEmptyMatrix(N_c1, N_s2);
	wsReal	*Ra_itm	= sseVector(N_s2);
	for (i=0 ; i<N_r1 ; i++) {
		/* Build itm for ith row of SYMmat */
		memcpy(Ra_itm, Ra_m2[i], sizeof(wsReal)*(i+1));
		for (wsUint q=i+1 ; q<N_s2 ; q++)
			Ra_itm[q] = Ra_m2[q][i];
		for (j=0 ; j<N_c1 ; j++) {
			wsUint N_med = 0;
#ifdef USE_SSE
			sse_t sse_A = sseSet(Ra_m1[i][j]);
			N_med = getMed(N_s2);
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_R = (sse_t *)(Ra_ret[j] + k);
				sse_t *sse_T = (sse_t *)(Ra_itm + k);
				*sse_R = sseAdd(*sse_R, sseMul(sse_A, *sse_T));
			}
#endif
			for (k=N_med ; k<N_s2 ; k++)
				Ra_ret[j][k] += Ra_m1[i][j]*Ra_itm[k];
		}
	}
	sseFree(Ra_itm);

	return Ra_ret;
}

wsVec sseVpM(wsVec Rp_m1, wsUint N_c1,
	 wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsReal **Rp_rMean/*=NULL*/)
{
	/* Check */
	if (N_c1 != N_r2)
		halt("Can't do matrix multplication, C(mat1)=[%u] != [%u]=R(mat2)",
		N_c1, N_r2);

	/* Make return matrix */
	wsVec	Ra_ret	= sseEmptyVec(N_c2);
	if (Rp_rMean != NULL) {
		sseCalloc(*Rp_rMean, wsReal, 1);
	}

	/* Calc */

	LOOP(j, N_c1) {
		wsReal *Rp_2	= Rp_m2[j];
		wsReal R_1j		= Rp_m1[j];
		wsUint N_med	= 0;
#ifdef USE_SSE
		N_med = getMed(N_c2);
		sse_t sse_1ij = sseSet(R_1j);
		/* SSE-enabled part */
		SSE_LOOP(k, N_med) {
			// m2[j,k:(k+3)]
			sse_t *sse_2jk	= (sse_t *)(Rp_2 + k);
			// m3[i,k:(k+3)] (j)
			sse_t *sse_mR	= (sse_t *)(Ra_ret + k);

			*sse_mR = sseAdd(*sse_mR, sseMul(sse_1ij, *sse_2jk));
		}
#endif
		/* Rest part */
		LOOPV(k, N_med, N_c2) Ra_ret[k] += R_1j * Rp_2[k];
	}
	/* Mean process */
	if (Rp_rMean)
		(*Rp_rMean)[0] = sseVsum(Ra_ret, N_c2) / (wsReal)N_c2;

	/* Return */
	return Ra_ret;
}

wsSym sym_sseMpM(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2)
{
	/* Check */
	if (N_c1 != N_r2)
		halt("Can't do matrix multplication, C(mat1)=[%u] != [%u]=R(mat2)",
			N_c1, N_r2);

	/* Make return matrix */
	wsMat Ra_ret	= sseEmptyMatrix(N_r1, N_c2);
	wsUint	i, j, k;

	/* Calc */
	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_R		= Ra_ret[i];
		wsReal	*Rp_1		= Rp_m1[i];

		for (j=0 ; j<N_c1 ; j++) {
			wsReal *Rp_2	= Rp_m2[j];
			wsReal R_1j		= Rp_1[j];
			wsUint N_med	= 0;
#ifdef USE_SSE
			N_med = getMed(i+1);
			sse_t sse_1ij = sseSet(R_1j);
			/* SSE-enabled part */
			for (k=0 ; k<N_med ; k+=sseJmp) {
				// m2[j,k:(k+3)]
				sse_t *sse_2jk	= (sse_t *)(Rp_2 + k);
				// m3[i,k:(k+3)] (j)
				sse_t *sse_mR	= (sse_t *)(Rp_R + k);

				*sse_mR = sseAdd(*sse_mR, sseMul(sse_1ij, *sse_2jk));
			}
#endif
			/* Rest part */
			for (k=N_med ; k<=i ; k++)
				Rp_R[k] += R_1j * Rp_2[k];
		}
	}

	/* Return */
	return Ra_ret;
}

void multMpM(wsMat Rp_m1, wsUint N_r1, wsUint N_c1, wsUint N_y1, wsUint N_x1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsUint N_y2, wsUint N_x2,
	wsMat Ra_T, wsUint N_rT, wsUint N_cT, wsUint N_yT, wsUint N_xT)
{
	wsUint	i, j, k;

	N_c1 -= N_x1;
	N_r1 -= N_y1;
	N_c2 -= N_x2;
	N_r2 -= N_y2;

	/* Check */
	if (N_c1 != N_r2)
		halt("Can't do matrix multplication, C(mat1)=[%u] != [%u]=R(mat2)",
			N_c1, N_r2);

	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_matR	= Ra_T[i+N_yT] + N_xT;
		wsReal	*Rp_1		= Rp_m1[i+N_y1] + N_x1;

		for (j=0 ; j<N_c1 ; j++) {
			wsReal R_1j		= Rp_1[j];
			wsReal *Rp_2	= Rp_m2[j+N_y2] + N_x2;

			/* Rest part */
			for (k=0 ; k<N_c2 ; k++)
				Rp_matR[k] += R_1j * Rp_2[k];
		}
	}
}

wsMat sseMpMt(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsReal **Rp_rMean/*=NULL*/)
{
	wsUint	i, j, k;
	/* Make return matrix */
	wsMat	Ra_T	= sseEmptyMatrix(N_r1, N_r2);
//	wsUint	N_med	= getMed(N_c2);

	if (Rp_rMean != NULL) {
		wsCalloc(*Rp_rMean, wsReal, N_r1);
	}
	/* Check */
	if (N_c1 != N_c2)
		halt("Can't do matrix multplication, C(mat1)=[%u] != [%u]=C(mat2)",
		N_c1, N_c2);

	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_matR	= Ra_T[i];
		wsReal	*Rp_1		= Rp_m1[i];

		for (j=0 ; j<N_r2 ; j++) {
			wsReal *Rp_2	= Rp_m2[j];
			wsUint N_med	= 0;

#ifdef USE_SSE
			N_med = getMed(N_c2);
			if (N_med != 0) {
				sse_t sse_R = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_m1 = (sse_t *)(Rp_1 + k);
					sse_t *sse_m2 = (sse_t *)(Rp_2 + k);

					sse_R = sseAdd(sse_R, sseMul(*sse_m1, *sse_m2));
				}
				sseSum(sse_R, Rp_matR[j]);
			}
#endif
			for (k=N_med ; k<N_c2 ; k++)
				Rp_matR[j] += Rp_1[k] * Rp_2[k];
			/* Mean process */
			if (Rp_rMean) (*Rp_rMean)[i] += Rp_matR[j];
		}
		/* Mean process */
		if (Rp_rMean) (*Rp_rMean)[i] /= (wsReal)N_r2;
	}
	return Ra_T;
}

void multMpMt(wsMat Rp_m1, wsUint N_r1, wsUint N_c1, wsUint N_y1, wsUint N_x1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsUint N_y2, wsUint N_x2,
	wsMat Ra_T, wsUint N_rT, wsUint N_cT, wsUint N_yT, wsUint N_xT,
	wsReal **Rp_rMean/*=NULL*/)
{
	wsUint	i, j, k;

	N_c1 -= N_x1;
	N_r1 -= N_y1;
	N_c2 -= N_x2;
	N_r2 -= N_y2;

	if (Rp_rMean != NULL) {
		wsCalloc(*Rp_rMean, wsReal, N_r1);
	}
	/* Check */
	if (N_c1 != N_c2)
		halt("Can't do matrix multplication, C(mat1)=[%u] != [%u]=C(mat2)",
		N_c1, N_c2);

	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_matR	= Ra_T[i+N_yT] + N_xT;
		wsReal	*Rp_1		= Rp_m1[i+N_y1] + N_x1;

		for (j=0 ; j<N_r2 ; j++) {
			wsReal *Rp_2	= Rp_m2[j+N_y2] + N_x2;
			wsReal R_R		= W0;
			for (k=0 ; k<N_c2 ; k++)
				R_R += Rp_1[k] * Rp_2[k];
			Rp_matR[j] = R_R;
			/* Mean process */
			if (Rp_rMean) (*Rp_rMean)[i] += Rp_matR[j];
		}
		/* Mean process */
		if (Rp_rMean) (*Rp_rMean)[i] /= (wsReal)N_r2;
	}
}

/* sseMMt, return sym. */
wsSym sseMpMt(wsReal **Rp_m, wsUint N_r1, wsUint N_c1)
{
	/* Make return matrix */
	wsReal **Ra_ret	= sseEmptyMatrix(N_r1, N_r1);
	wsUint	N_med	= getMed(N_c1);
	wsUint	i, j, k;

	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_matR	= Ra_ret[i];
		wsReal	*Rp_m1		= Rp_m[i];

//#ifdef USE_SYM
//		for (j=0 ; j<=i ; j++) {
//#else
		for (j=0 ; j<N_r1 ; j++) {
//#endif
			wsReal *Rp_m2	= Rp_m[j];

#ifdef USE_SSE
			if (N_med != 0) {
				sse_t sse_R = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_m1 = (sse_t *)(Rp_m1 + k);
					sse_t *sse_m2 = (sse_t *)(Rp_m2 + k);

					sse_R = sseAdd(sse_R, sseMul(*sse_m1, *sse_m2));
				}
				sseSum(sse_R, Rp_matR[j]);
			}
			for (k=N_med ; k<N_c1 ; k++)
				Rp_matR[j] += Rp_m1[k] * Rp_m2[k];
#else
			for (k=0 ; k<N_c1 ; k++)
				Rp_matR[j] += Rp_m1[k] * Rp_m2[k];
#endif
		}
	}

	return Ra_ret;
}

wsSym sym_sseMpMt(wsReal **Rp_m1, wsUint N_r1, wsUint N_c1,
	wsReal **Rp_m2/*=NULL*/, wsUint N_r2/*=0*/, wsUint N_c2/*=0*/,
	wsReal **Rp_rMean/*=NULL*/)
{
	wsUint	i, j, k;
	/* Make return matrix */
	wsSym Ra_T	= sseEmptySymMat(N_r1);
//	wsUint	N_med	= getMed(N_c2);
	if (Rp_m2 == NULL) {
		Rp_m2 = Rp_m1;
		N_r2 = N_r1;
		N_c2 = N_c1;
	}

	if (Rp_rMean != NULL) {
		wsCalloc(*Rp_rMean, wsReal, N_r1);
	}
	/* Check */
	if (N_c1 != N_c2)
		halt("Can't do matrix multplication, C(mat1)=[%u] != [%u]=C(mat2)",
			N_c1, N_c2);
	if (N_r1 != N_r2)
		halt("Can't do matrix multplication, R(mat1)=[%u] != [%u]=R(mat2)",
			N_r1, N_r2);

	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_matR	= Ra_T[i];
		wsReal	*Rp_1		= Rp_m1[i];

		for (j=0 ; j<=i ; j++) {
			wsReal *Rp_2	= Rp_m2[j];
			wsUint N_med	= 0;

#ifdef USE_SSE
			N_med = getMed(N_c1);
			if (N_med != 0) {
				sse_t sse_R = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_m1 = (sse_t *)(Rp_1 + k);
					sse_t *sse_m2 = (sse_t *)(Rp_2 + k);

					sse_R = sseAdd(sse_R, sseMul(*sse_m1, *sse_m2));
				}
				sseSum(sse_R, Rp_matR[j]);
			}
#endif
			for (k=N_med ; k<N_c1 ; k++)
				Rp_matR[j] += Rp_1[k] * Rp_2[k];
			/* Mean process */
			if (Rp_rMean) (*Rp_rMean)[i] += Rp_matR[j];
		}
		/* Mean process */
		if (Rp_rMean) (*Rp_rMean)[i] /= (wsReal)N_c1;
	}
	return Ra_T;
}

wsMat sseMpD(wsMat Ra_m1, wsUint N_r, wsUint N_c, wsReal *Ra_D)
{
	wsReal **Ra_ret = sseMatrix(N_r, N_c);
	/* MD[i,j] = M[i,j]*D[j] */
	for (wsUint i=0 ; i<N_r ; i++) {
		wsReal *Ra_T = Ra_ret[i];
		wsReal *Ra_1 = Ra_m1[i];
		wsUint j, N_med=0;
#ifdef USE_SSE
		for (j=N_med ; j<N_c ; j+=sseJmp) {
			sse_t *sse_T = (sse_t *)(Ra_T + j);
			sse_t *sse_1 = (sse_t *)(Ra_1 + j);
			sse_t *sse_2 = (sse_t *)(Ra_D + j);
			*sse_T = sseMul(*sse_1, *sse_2);
			//			Ra_T[j] = Ra_1[j] * Ra_2[j];
		}
#endif
		for (j=N_med ; j<N_c ; j++)
			Ra_T[j] = Ra_1[j] * Ra_D[j];
	}

	return Ra_ret;
}

// Do mt = m*s
void sseMpC(wsMatCst Ra_m, wsRealCst R_s, wsMat Ra_mt, wsUintCst N_row,
			wsUint N_col/*=0xffffffff*/)
{
	wsUint i, j;
	if (N_col == 0xffffffff)
		N_col = N_row;

	sse_t	sse_S	= sseSet(R_s);

	for (i=0 ; i<N_row ; i++) {
		wsUint	N_med	= 0;
#ifdef USE_SSE
		N_med	= getMed(N_col);
		for (j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m = (sse_t *)(Ra_m[i] + j);
			sse_t *sse_mt = (sse_t *)(Ra_mt[i] + j);

			*sse_mt = sseMul(*sse_m, sse_S);
		}
#endif
		/* Do rest part */
		for (j=N_med ; j<N_col ; j++)
			Ra_mt[i][j] = Ra_m[i][j]*R_s;
	}
}

// Do mt = m*s
void sseSpC(wsSymCst Ra_m, wsRealCst R_s, wsSym Ra_mt, wsUintCst N_sz)
{
	wsUint i, j;

	sse_t	sse_S	= sseSet(R_s);

	for (i=0 ; i<N_sz ; i++) {
		wsUint	N_med	= 0;
#ifdef USE_SSE
		N_med	= getMed(i+1);
		for (j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m	= (sse_t *)(Ra_m[i] + j);
			sse_t *sse_mt	= (sse_t *)(Ra_mt[i] + j);

			*sse_mt = sseMul(*sse_m, sse_S);
		}
#endif
		/* Do rest part */
		for (j=N_med ; j<=i ; j++)
			Ra_mt[i][j] = Ra_m[i][j]*R_s;
	}
}

wsSym sseStS(wsSym Ra_s, wsUint N_sz)
{
	wsSym Ra_ret = sseEmptySymMat(N_sz);
	for (wsUint i=0 ; i<N_sz ; i++) {
		wsReal *S = Ra_s[i];

		for (wsUint j=0 ; j<=i ; j++) {
			wsUint	N_med	= 0;
			wsReal	R_Sij	= S[j];
			wsReal	*R		= Ra_ret[j];
#ifdef USE_SSE
			N_med = getMed(j+1);
			sse_t sse_Sij = sseSet(R_Sij);
			for (wsUint k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_Rjk = (sse_t *)(R + k);
				sse_t *sse_Sik = (sse_t *)(S + k);
				*sse_Rjk = sseAdd(*sse_Rjk, sseMul(sse_Sij, *sse_Sik));
			}
#endif
			for (wsUint k=N_med ; k<=j ; k++)
				R[k] += R_Sij * S[k];
		}
	}

	return Ra_ret;
}

wsSym sseSkS(wsSym Ra_s1, wsUint N_s1, wsSym Ra_s2, wsUint N_s2)
{
	wsSym Ra_ret = sseSymMat(N_s1*N_s2);

	for (wsUint i=0 ; i<N_s1 ; i++) {
		wsReal R_s1ij;
		for (wsUint j=0 ; j<i ; j++) {
			R_s1ij = Ra_s1[i][j];
			for (wsUint k=0 ; k<N_s2 ; k++) {
				for (wsUint l=0 ; l<=k ; l++)
					Ra_ret[i*N_s2+k][j*N_s2+l] =
					Ra_ret[i*N_s2+l][j*N_s2+k] = R_s1ij*Ra_s2[k][l];
// 				for (wsUint l=k+1 ; l<N_s2 ; l++)
// 					Ra_ret[i*N_s2+k][j*N_s2+l] = R_s1ij*Ra_s2[l][k];
			}
		}
		R_s1ij = Ra_s1[i][i];
		for (wsUint k=0 ; k<N_s2 ; k++)
			for (wsUint l=0 ; l<=k ; l++)
				Ra_ret[i*N_s2+k][i*N_s2+l] = R_s1ij * Ra_s2[k][l];
	}

	return Ra_ret;
}

void sseSsS(wsSym Ra_m1, wsSym Ra_m2, wsSym Ra_mt, wsUint N_sz,
	wsReal R_multiplier/*=W1*/)
{
	if (R_multiplier == W0) {
		for (wsUint i=0 ; i<N_sz ; i++) {
			wsUint N_med = 0;
#ifdef USE_SSE
			N_med = getMed(i+1);
			for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_m1 = (sse_t *)(Ra_m1[i] + j);
				sse_t *sse_m2 = (sse_t *)(Ra_m2[i] + j);
				sse_t *sse_mt = (sse_t *)(Ra_mt[i] + j);

				*sse_mt = sseSub(*sse_m1, *sse_m2);
			}
#endif
			/* Do rest part */
			for (wsUint j=N_med ; j<=i ; j++)
				Ra_mt[i][j] = Ra_m1[i][j]-Ra_m2[i][j];
		}
	} else {
		sse_t	sse_M = sseSet(R_multiplier);
		for (wsUint i=0 ; i<N_sz ; i++) {
			wsUint N_med = 0;
#ifdef USE_SSE
			N_med = getMed(i+1);
			for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
				sse_t *sse_m1 = (sse_t *)(Ra_m1[i] + j);
				sse_t *sse_m2 = (sse_t *)(Ra_m2[i] + j);
				sse_t *sse_mt = (sse_t *)(Ra_mt[i] + j);

				*sse_mt = sseMul(sseSub(*sse_m1, *sse_m2), sse_M);
			}
#endif
			/* Do rest part */
			for (wsUint j=N_med ; j<=i ; j++)
				Ra_mt[i][j] = (Ra_m1[i][j]-Ra_m2[i][j])*R_multiplier;
		}
	}
}

wsMat multMpMt(wsMat Rp_m1, wsUint N_r1, wsUint N_c1,
	wsMat Rp_m2, wsUint N_r2, wsUint N_c2, wsReal **Rp_rMean/*=NULL*/)
{
	if (N_c1 != N_c2)
		halt("DIMERR");
	wsMat Ra_ret = sseMatrix(N_r1, N_r2);
	multMpMt(Rp_m1, N_r1, N_c1, 0, 0, Rp_m2, N_r2, N_c2, 0, 0,
		Ra_ret, N_r1, N_r2, 0, 0, Rp_rMean);
	return Ra_ret;
}

wsReal** multMS(wsReal **Ra_mat, wsUint N_row, wsUint N_col, wsReal **Ra_symMat, wsUint N_sz)
{
	wsUint i, j, k;
	if (N_col != N_sz)
		halt("Column size of first matrix [%d] is not match with dimension of second symmetric matrix [%d]",
		N_col, N_sz);
	wsReal **Ra_ret = sseEmptyMatrix(N_row, N_sz);

	for (i=0; i<N_row ; i++) {			//MS Multiplication 수행
		for (j=0; j<N_sz ; j++) {
			Ra_ret[i][j] = 0;
			for (k=0; k<N_sz; k++) {
				if (j >= k)
					Ra_ret[i][j] +=Ra_mat[i][k]*Ra_symMat[k][j];
				else
					Ra_ret[i][j] +=Ra_mat[i][k]*Ra_symMat[j][k];
			}
		}
	}

	return Ra_ret;
}

wsMat multMP(wsMat Ra_mat, wsUint N_row, wsUint N_col, wsReal *Rp_data,
			 wsUint N_pr, wsUint N_pc, wsUint *Na_rEnd, wsUint *Na_cIdx)
{
	if (N_col != N_pr)
		halt("Column size of first matrix [%d] is not match with row size "
		"of second sparase matrix [%d]", N_col, N_pr);

	wsMat		Ra_ret	= sseEmptyMatrix(N_row, N_pc);
	wsUint		i		= 0;
	wsUint		N_s		= 0;
	wsUint		N_e;
	do {
		/* Set end point */
		N_e = Na_rEnd[i];

		/* No element in the row */
		if (N_s != N_e) {
			/* Sum up them */
			for (wsUint r=0 ; r<N_row ; r++)
				for (wsUint j=N_s ; j<N_e ; j++) {
					wsUint TC = Na_cIdx[j];
					Ra_ret[r][TC] += Ra_mat[r][i] * Rp_data[j];
				}
		}

		/* Update start point */
		i++;
		if (i < N_pr)
			N_s = Na_rEnd[i];
		else break;
	} while (1);

	return Ra_ret;
}

wsReal** multSMt(wsSym Ra_m1, wsUint N_sz, wsReal **Ra_m2, wsUint N_r2,
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

			//			if (L != N_med) halt("ERROR1");
			for(wsUint a=N_med ; a <=i ; a++)
				R_sum += Ra_m1[i][a] * Ra_m2[j][a];
			//			/* SANITY */ L = i+1;

			/* Vertical part */
			N_med = i+1;
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
	return multMMt(Ra_m1, N_sz, N_sz, Ra_m2, N_r2, N_c2);
#endif
}

/*
sseTrSS <- function(s1) {
  s1.lower <- lower.tri(s1)
  return (2*sum(s1[s1.lower]^2) + sum(diag(s1)^2))
}
sseTrSS2 <- function(s1) {
  r <- 0
  for (i in 1:nrow(s1)) {
    if (i>1) r <- r + 2*sum(s1[i,1:(i-1)]^2)
	r <- r + s1[i,i]^2
  }
  r
}
*/
wsReal sseTrSS(wsSym Ra_s1, wsUint N_s1)
{
	wsReal R_ret = W0;

	/* Sqsum of each cell and return
		* diag term is *1, otherwise *2 */
	for (wsUint i=0 ; i<N_s1 ; i++)
		R_ret += W2 * sseVsqsum(Ra_s1[i], i) + SQR(Ra_s1[i][i]);

	return R_ret;
}

wsReal sseTrSS(wsSym Ra_s1, wsUint N_s1, wsSym Ra_s2)
{
	wsReal R_ret = W0;

	/* Sqsum of each cell and return
		* diag term is *1, otherwise *2 */
	for (wsUint i=0 ; i<N_s1 ; i++) {
		R_ret += W2 * sseVV(Ra_s1[i], i, Ra_s2[i])
			+ Ra_s1[i][i]*Ra_s2[i][i];
	}

	return R_ret;
}

inline sse_t sseNeg(sse_t &x) {
	static const sse_t sign_mask = sseSet(REAL_CONST(-0.)); // -0.f = 1 << 31
	return sseXor(sign_mask, x);
}

void sseVneg(wsReal *Ra_v1, wsUint N_sz, wsReal *Ra_vt/*=NULL*/, wsReal R_add/*=WISARD_NAN*/)
{
	wsUint N_med = 0, i;

#ifdef USE_SSE
	N_med = getMed(N_sz);
#endif

	if (Ra_vt == NULL) {
		if (NA(R_add)) {
#ifdef USE_SSE
			for (i=0 ; i<N_med ; i+=sseJmp) {
				sse_t *sse_S = (sse_t *)(Ra_v1 + i);
				*sse_S = sseNeg(*sse_S);
			}
#endif
			for (i=N_med ; i<N_sz ; i++)
				Ra_v1[i] = -Ra_v1[i];
		} else {
#ifdef USE_SSE
			sse_t sse_A = sseSet(R_add);
			for (i=0 ; i<N_med ; i+=sseJmp) {
				sse_t *sse_S = (sse_t *)(Ra_v1 + i);
				*sse_S = sseNeg(*sse_S);
				*sse_S = sseAdd(*sse_S, sse_A);
			}
#endif
			for (i=N_med ; i<N_sz ; i++)
				Ra_v1[i] = -Ra_v1[i] + R_add;
		}
	} else {
		if (NA(R_add)) {
#ifdef USE_SSE
			for (i=0 ; i<N_med ; i+=sseJmp) {
				sse_t *sse_T = (sse_t *)(Ra_vt + i);
				sse_t *sse_S = (sse_t *)(Ra_v1 + i);

				*sse_T = sseNeg(*sse_S);
			}
#endif
			for (i=N_med ; i<N_sz ; i++)
				Ra_vt[i] = -Ra_v1[i];
		} else {
#ifdef USE_SSE
			sse_t sse_A = sseSet(R_add);
			for (i=0 ; i<N_med ; i+=sseJmp) {
				sse_t *sse_T = (sse_t *)(Ra_vt + i);
				sse_t *sse_S = (sse_t *)(Ra_v1 + i);

				*sse_T = sseNeg(*sse_S);
				*sse_T = sseAdd(*sse_T, sse_A);
			}
#endif
			for (i=N_med ; i<N_sz ; i++)
				Ra_vt[i] = -Ra_v1[i] + R_add;
		}
	}
}

void sseVnegMul(wsReal R_mul, wsReal *Ra_v1, wsUint N_sz, wsReal *Ra_vt/*=NULL*/)
{
	wsUint N_med = 0, i;

#ifdef USE_SSE
	N_med = getMed(N_sz);
#endif

	if (Ra_vt == NULL) {
#ifdef USE_SSE
		sse_t sse_M = sseSet(R_mul);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_S = (sse_t *)(Ra_v1 + i);
			*sse_S = sseMul(sseNeg(*sse_S), sse_M);
		}
#endif
		for (i=N_med ; i<N_sz ; i++)
			Ra_v1[i] = -Ra_v1[i] * R_mul;
	} else {
#ifdef USE_SSE
		sse_t sse_M = sseSet(R_mul);
		for (i=0 ; i<N_med ; i+=sseJmp) {
			sse_t *sse_T = (sse_t *)(Ra_vt + i);
			sse_t *sse_S = (sse_t *)(Ra_v1 + i);

			*sse_T = sseMul(sseNeg(*sse_S), sse_M);
		}
#endif
		for (i=N_med ; i<N_sz ; i++)
			Ra_vt[i] = -Ra_v1[i] * R_mul;
	}
}

wsReal** sseMpB(wsReal **Rp_mat1, wsUint N_r1, wsUint N_c1,
	wsReal **Rp_bMat, wsUint N_s2, wsUint N_blkRep)
{
	/* N_c1*N_blkRep should be N_s2 */
	if (N_c1 != /*N_s2**/N_blkRep)
		halt("Column of first matrix[%d] should be dimension of second "
		"block matrix[%d] times repeatition[%d]", N_c1, N_s2, N_blkRep);

	/*
	rmat <- matrix(0, nrow=dm[1], ncol=dm[2])
	for (i in 1:dm[1]) {
		for (j in 1:dm[2]) {
			x <- ceiling(j/db[1])
			rmat[i,j] <- sum(m[i,((x-1)*db[1]+1):(x*db[1])]*b[,((j-1)%%db[1])+1])
		}
	}
	*/
	wsReal **Ra_ret = NULL;
	if (N_blkRep == 1) {
		Ra_ret = sseEmptyMatrix(N_r1, N_c1);
		for (wsUint j=0 ; j<N_blkRep ; j++) {
			multMpM(Rp_mat1, N_r1, N_c1-(N_blkRep-(j+1))*N_s2, 0, j*N_s2,
				Rp_bMat, N_s2, N_s2, 0, 0,
				Ra_ret, N_r1, N_c1, 0, j*N_s2);
//		wsReal **Ra_1 = sseSubsetMatrix(Rp_mat1, 0, N_r1, j*N_s2, N_s2);
//		wsReal **Ra_R = sseMpM(Ra_1, N_r1, N_s2, Rp_bMat, N_s2, N_s2);
//		for (wsUint k=0 ; k<N_c1 ; k++)
//			memcpy(Ra_ret[k]+j*N_s2, Ra_R[k], sizeof(wsReal)*N_s2);
		}
	} else if (N_blkRep == N_c1) {
/*
		cb <- colSums(b)
		for (i in 1:cx) {
			s <- (i-1)*K + 1
			e <- i*K
			# x[1,i]*col(b)[1] ... x[1,i]*col(b)[k]
			#   ...
			# x[p,i]*col(b)[1] ... x[p,i]*col(b)[k]
			R[,s:e] <- x[,i,drop=F] %x% t(cb)
		}
*/
		wsUint N_sz = N_s2*N_blkRep;
		Ra_ret = sseEmptyMatrix(N_r1, N_sz);
		wsReal *colS = sseMsumCol(Rp_bMat, N_s2);
		for (wsUint j=0 ; j<N_r1 ; j++) for (wsUint i=0,k=0 ; i<N_c1 ; i++) {
			for (wsUint l=0 ; l<N_s2 ; l++,k++)
				Ra_ret[j][k] = Rp_mat1[j][i] * colS[l];
		}
	} else halt("Noncompatible calculation tried");

	return Ra_ret; 
}

/* Do multiply M %*% t(B), and assume that B is block-diagonal matrix of B'
 * which B' repeats	 N_blkRep times
 * 
 * So return matrix become row(M) * (N_s2*N_blkRep) dimension */
wsReal** sseMpBt(wsReal **Rp_mat1, wsUint N_r1, wsUint N_c1,
	wsReal **Rp_bMat, wsUint N_s2, wsUint N_blkRep)
{
#if 0
	/* N_c1*N_blkRep should be N_s2 */
	if (N_c1*N_blkRep != N_s2)
		halt("Column of first matrix[%d] should be dimension of second "
			"block matrix[%d] times repeatition[%d]", N_c1, N_s2, N_blkRep);

	REAL_t **Ra_ret = sseMatrix(N_r1, N_c1);
	for (UINT_t i=0 ; i<N_r1 ; i++) {
		for (UINT_t j=0 ; j<N_blkRep ; j++) {
			REAL_t *Rp_curPos = Rp_mat1[i] + (j*N_s2);
			/* Loop N_blkRep times, multiply A[i,[j~(j+1)*N_sz)] %*% B = C[i,[j~(j+1)*N_sz)] */
			REAL_t **Ra_ret =
				sseMpMt(&Rp_curPos, 1, N_s2, Rp_bMat, N_s2, N_s2);
			memcpy(Ra_ret[i] + (j*N_s2), Ra_ret[0], sizeof(REAL_t)*N_s2);
			sseUnmat(Ra_ret, 1);
		}
	}

	return Ra_ret;
#else
	/* N_c1*N_blkRep should be N_s2 */
	if (N_c1*N_blkRep != N_s2)
		halt("Column of first matrix[%d] should be dimension of second "
		"block matrix[%d] times repeatition[%d]", N_c1, N_s2, N_blkRep);

	/*
	for (i in 1:dm[1]) {
		for (j in 1:dm[2]) {
			x <- ceiling(j/db[1])
			rmat[i,j] <- sum(m[i,((x-1)*db[1]+1):(x*db[1])]*b[,((j-1)%%db[1])+1])
		}
	}
	*/
	wsReal **Ra_ret = sseMatrix(N_r1, N_c1);
	for (wsUint j=0 ; j<N_blkRep ; j++) {
		wsReal **Ra_1 = sseMsubset(Rp_mat1, 0, N_r1, j*N_s2, N_s2);
		wsReal **Ra_R = sseMpMt(Ra_1, N_r1, N_s2, Rp_bMat, N_s2, N_s2);
		for (wsUint k=0 ; k<N_c1 ; k++)
			memcpy(Ra_ret[k]+j*N_s2, Ra_R[k], sizeof(wsReal)*N_s2);
	}
	return Ra_ret; 
#endif
}

wsMat sseSpD(wsSym Ra_s, wsUint N_sz, wsDiag Ra_d)
{
	wsReal **Ra_ret = sseMatrix(N_sz, N_sz);
	for (wsUint i=0 ; i<N_sz ; i++) {
		wsReal *Ra_R = Ra_ret[i];
		wsReal *Ra_S = Ra_s[i];
		wsUint N_med = 0;
#ifdef USE_SSE
		N_med = getMed(i+1);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_D = (sse_t *)(Ra_d + j);
			sse_t *sse_S = (sse_t *)(Ra_S + j);
			sse_t *sse_R = (sse_t *)(Ra_R + j);

			*sse_R = sseMul(*sse_D, *sse_S);
		}
#endif
		for (wsUint j=N_med ; j<=i ; j++)
			Ra_R[j] = Ra_d[j] * Ra_S[j];
		for (wsUint j=i+1 ; j<N_sz ; j++)
			Ra_R[j] = Ra_d[j] * Ra_s[j][i];
	}
	return Ra_ret;
}

wsSym sseSpDpSt(wsSym Ra_s, wsUint N_sz, wsDiag Ra_d)
{
	wsSym	Ra_ret	= sseSymMat(N_sz);
	wsReal	*Ra_v	= sseVector(N_sz);

	for (wsUint i=0 ; i<N_sz ; i++) {
		//wsReal *Ra_R = Ra_ret[i];

		wsReal *Ra_S = Ra_s[i];
		wsUint N_med = 0;
#ifdef USE_SSE
		N_med = getMed(i+1);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_D = (sse_t *)(Ra_d + j);
			sse_t *sse_S = (sse_t *)(Ra_S + j);
			sse_t *sse_R = (sse_t *)(Ra_v + j);

			*sse_R = sseMul(*sse_D, *sse_S);
		}
#endif
		for (wsUint j=N_med ; j<=i ; j++)
			Ra_v[j] = Ra_d[j] * Ra_S[j];
		for (wsUint j=i+1 ; j<N_sz ; j++)
			Ra_v[j] = Ra_d[j] * Ra_s[j][i];

		sseVpS(Ra_v, N_sz, Ra_s, i+1, Ra_ret[i]);
	}
	return Ra_ret;
}

wsSym sseSpB(wsSym Ra_s, wsUint N_sz, wsMat Ra_B, wsUint N_szBlk, wsUint N_rep)
{
	/* N_c1*N_blkRep should be N_s2 */
	if (N_sz != N_szBlk*N_rep)
		halt("Size of first matrix[%d] should be dimension of second "
			"block matrix[%d] times repeatition[%d]", N_sz, N_szBlk, N_rep);

	/*
	for (i in 1:dm[1]) {
		for (j in 1:dm[2]) {
			x <- ceiling(j/db[1])
			rmat[i,j] <- sum(m[i,((x-1)*db[1]+1):(x*db[1])]*b[,((j-1)%%db[1])+1])
		}
	}
	*/
	wsReal **Ra_ret = sseEmptyMatrix(N_sz, N_sz);
	for (wsUint i=0 ; i<N_sz ; i++) {
		for (wsUint j=0 ; j<N_szBlk ; j++) {
			for (wsUint k=0 ; k<N_sz ; k++) {
				wsUint K=k%N_szBlk;
				wsUint J=(wsUint)(k/N_szBlk);
				if (k<=i)
					Ra_ret[i][J*N_szBlk+j] += Ra_s[i][k]*Ra_B[K][j];
				else
					Ra_ret[i][J*N_szBlk+j] += Ra_s[k][i]*Ra_B[K][j];
			}
		}
	}

//		multMpM(Rp_mat1, N_r1, N_c1-(N_blkRep-(j+1))*N_s2, 0, j*N_s2,
//			Rp_bMat, N_s2, N_s2, 0, 0,
//			Ra_ret, N_r1, N_c1, 0, j*N_s2);
//		wsReal **Ra_1 = sseSubsetMatrix(Rp_mat1, 0, N_r1, j*N_s2, N_s2);
//		wsReal **Ra_R = sseMpM(Ra_1, N_r1, N_s2, Rp_bMat, N_s2, N_s2);
//		for (wsUint k=0 ; k<N_c1 ; k++)
//			memcpy(Ra_ret[k]+j*N_s2, Ra_R[k], sizeof(wsReal)*N_s2);
//	}
	return Ra_ret; 
}

wsReal** sseSym2Mat(wsSym Ra_s, wsUint N_sz)
{
	wsReal **Ra_ret = sseMatrix(N_sz, N_sz);

	for (wsUint i=0 ; i<N_sz ; i++) {
		for (wsUint j=0 ; j<i ; j++)
			Ra_ret[i][j] = Ra_ret[j][i] = Ra_s[i][j];
		Ra_ret[i][i] = Ra_s[i][i];
	}

	return Ra_ret;
}

wsReal* sseMmeanR(wsReal **Ra_mat, wsUint N_row, wsUint N_col/*=0xffffffff*/)
{
	wsReal *Ra_ret	= NULL;
	sseCalloc(Ra_ret, wsReal, N_row);

	/* Allocate default value for N_col as equal to N_row */
	if (N_col == 0xffffffff)
		N_col = N_row;

	for (wsUint i=0 ; i<N_row ; i++) {
		wsReal *Rp_m	= Ra_mat[i];
		wsUint N_med	= 0;
#ifdef USE_SSE
		N_med = getMed(N_col);
		sse_t sse_sum = sseSet(0.0);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m = (sse_t *)(Rp_m + j);
			sse_sum = sseAdd(sse_sum, *sse_m);
		}
		sseSum(sse_sum, Ra_ret[i]);
#endif
		for (wsUint j=N_med ; j<N_col ; j++)
			Ra_ret[i] += Rp_m[j];
		Ra_ret[i] /= (wsReal)N_col;
	}

	return Ra_ret;
}

wsReal* sseSmeanR(wsSym Ra_mat, wsUint N_sz)
{
	//wsReal R_sum	= W0;
	wsReal *Ra_ret	= NULL;
	sseCalloc(Ra_ret, wsReal, N_sz);

	for (wsUint i=0 ; i<N_sz ; i++) {
		wsUint N_med = 0;
#ifdef USE_SSE
		N_med = getMed(i);
		sse_t sse_sumI = sseSet(0.0);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_M = (sse_t *)(Ra_mat[i] + j);
			sse_t *sse_T = (sse_t *)(Ra_ret + j);
			*sse_T		= sseAdd(*sse_T, *sse_M);
			sse_sumI	= sseAdd(sse_sumI, *sse_M);
		}
		sseAppendSum(sse_sumI, Ra_ret[i]);
#endif
		for (wsUint j=N_med ; j<i ; j++) {
			/* Add sum[j] and sum[i] */
			Ra_ret[j] += Ra_mat[i][j];
			Ra_ret[i] += Ra_mat[i][j];
		}
		Ra_ret[i] += Ra_mat[i][i];
	}
	sseVpC(Ra_ret, W1/(wsReal)N_sz, Ra_ret, N_sz);

	return Ra_ret;
}

wsReal* sseMmeanRavail(wsReal **Ra_mat, wsUint N_row, wsUint N_col/*=0xffffffff*/)
{
	wsReal *Ra_ret	= NULL;
	sseCalloc(Ra_ret, wsReal, N_row);

	/* Allocate default value for N_col as equal to N_row */
	if (N_col == 0xffffffff)
		N_col = N_row;

	for (wsUint i=0 ; i<N_row ; i++) {
		wsReal *Rp_m	= Ra_mat[i];
		wsUint N_med	= 0;
#ifdef USE_SSE
		N_med = getMed(N_col);
		sse_t sse_sum = sseSet(0.0);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m = (sse_t *)(Rp_m + j);
			sse_t sse_tmp = sseSet(WISARD_NA_REAL);	/* tmp = MISSING */
			sse_tmp = sseAnd(sseNeq(*sse_m, sse_tmp), *sse_m);

			sse_sum = sseAdd(sse_sum, sse_tmp);
		}
		sseSum(sse_sum, Ra_ret[i]);
#endif
		for (wsUint j=N_med ; j<N_col ; j++)
			if (!isMissingReal(Rp_m[j]))
				Ra_ret[i] += Rp_m[j];
		Ra_ret[i] /= (wsReal)N_col;
	}

	return Ra_ret;
}

wsReal* sseMmeanC(wsMat Ra_mat, wsUint N_row, wsUint N_col/*=0xffffffff*/)
{
	wsReal*	Ra_ret	= NULL;

	/* Allocate default value for N_col as equal to N_row */
	if (N_col == 0xffffffff)
		N_col = N_row;
	Ra_ret = sseMsumC(Ra_mat, N_row, N_col);
	sseVpC(Ra_ret, W1/N_row, Ra_ret, N_col);

	return Ra_ret;
}

wsReal* sseMmeanCavail(wsReal **Ra_mat, wsUint N_row, wsUint N_col/*=0xffffffff*/)
{
	wsReal*	Ra_ret	= NULL;
	wsReal*	Ra_cnt	= NULL;

	/* Allocate default value for N_col as equal to N_row */
	if (N_col == 0xffffffff)
		N_col = N_row;
	sseCalloc(Ra_ret, wsReal, N_col);
	sseCalloc(Ra_cnt, wsReal, N_col);

	for (wsUint i=0 ; i<N_row ; i++) {
		wsReal *Rp_m	= Ra_mat[i];
		wsUint N_med	= 0;
#ifdef USE_SSE
		N_med = getMed(N_col);
		//sse_t sse_sum = sseSet(0.0);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m = (sse_t *)(Rp_m + j);
			sse_t *sse_r = (sse_t *)(Ra_ret + j);
			sse_t *sse_c = (sse_t *)(Ra_cnt + j);
			sse_t sse_tmp = sseSet(WISARD_NA_REAL);	/* tmp = MISSING */
			sse_t	sse_1 = sseSet(1.0);
			sse_1 = sseAnd(sseNeq(*sse_m, sse_tmp), sse_1);
			sse_tmp = sseAnd(sseNeq(*sse_m, sse_tmp), *sse_m);

			*sse_r = sseAdd(*sse_r, sse_tmp);
			*sse_c = sseAdd(*sse_c, sse_1);
		}
#endif
		for (wsUint j=N_med ; j<N_col ; j++)
			if (!isMissingReal(Rp_m[j])) {
				Ra_ret[j] += Rp_m[j];
				Ra_cnt[j] += W1;
			}
	}
	sseVdV(Ra_ret, Ra_cnt, Ra_ret, N_col);

	return Ra_ret;
}

/*
 *
 * Row-sum
 * 
 */

wsReal* sseMsumR(wsReal **Ra_mat, wsUint N_row, wsUint N_col/*=0xffffffff*/,
				 wsReal *Rp_sumAll/*=NULL*/)
{
	wsReal R_sum	= W0;
	wsReal *Ra_ret	= NULL;
	sseCalloc(Ra_ret, wsReal, N_row);

	/* Allocate default value for N_col as equal to N_row */
	if (N_col == 0xffffffff)
		N_col = N_row;

	for (wsUint i=0 ; i<N_row ; i++) {
		wsReal *Rp_m	= Ra_mat[i];
		wsUint N_med	= 0;
#ifdef USE_SSE
		N_med = getMed(N_col);
		sse_t sse_sum = sseSet(0.0);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m = (sse_t *)(Rp_m + j);
			sse_sum = sseAdd(sse_sum, *sse_m);
		}
		sseSum(sse_sum, Ra_ret[i]);
#endif
		for (wsUint j=N_med ; j<N_col ; j++)
			Ra_ret[i] += Rp_m[j];
		R_sum += Ra_ret[i];
	}
	if (Rp_sumAll) *Rp_sumAll = R_sum;

	return Ra_ret;
}

wsReal* sseSsumR(wsSym Ra_mat, wsUint N_sz, wsReal *Rp_sumAll/*=NULL*/)
{
	//wsReal R_sum	= W0;
	wsReal *Ra_ret	= NULL;
	sseCalloc(Ra_ret, wsReal, N_sz);

	for (wsUint i=0 ; i<N_sz ; i++) {
		wsUint N_med = 0;
#ifdef USE_SSE
		N_med = getMed(i);
		sse_t sse_sumI = sseSet(0.0);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_M = (sse_t *)(Ra_mat[i] + j);
			sse_t *sse_T = (sse_t *)(Ra_ret + j);
			*sse_T		= sseAdd(*sse_T, *sse_M);
			sse_sumI	= sseAdd(sse_sumI, *sse_M);
		}
		sseAppendSum(sse_sumI, Ra_ret[i]);
#endif
		for (wsUint j=N_med ; j<i ; j++) {
			/* Add sum[j] and sum[i] */
			Ra_ret[j] += Ra_mat[i][j];
			Ra_ret[i] += Ra_mat[i][j];
		}
		Ra_ret[i] += Ra_mat[i][i];
	}
	if (Rp_sumAll)
		*Rp_sumAll = sseVsum(Ra_ret, N_sz);

	return Ra_ret;
}

wsReal* sseBsumR(wsMat Ra_mat, wsUint N_szBlk, wsUint N_sz,
				 wsReal *Rp_sumAll/*=NULL*/)
{
	wsUint N_cntRep = N_sz / N_szBlk;

	/* Allocate buffer */
	wsReal	R_sumA	= W0;
	wsReal	*Ra_ret	= sseVector(N_sz);

	/* Get block-sum */
	for (wsUint i=0 ; i<N_szBlk ; i++) {
		Ra_ret[i]	=  sseVsum(Ra_mat[i], N_szBlk);
		R_sumA		+= Ra_ret[i];
	}
	/* Get all-sum */
	if (Rp_sumAll) *Rp_sumAll = R_sumA * N_cntRep;
	/* Copy block-sum across entire */
	for (wsUint i=1,j=N_szBlk ; i<N_cntRep ; i++,j+=N_szBlk)
		memcpy(Ra_ret+j, Ra_ret, sizeof(wsReal)*N_szBlk);

	return Ra_ret;
}

/*
 *
 * Row-sqsum
 *
 */

wsReal* sseMssR(wsReal **Ra_mat, wsUint N_row, wsUint N_col/*=0xffffffff*/,
				 wsReal *Rp_sumAll/*=NULL*/)
{
	wsReal R_sum	= W0;
	wsReal *Ra_ret	= NULL;
	sseCalloc(Ra_ret, wsReal, N_row);

	/* Allocate default value for N_col as equal to N_row */
	if (N_col == 0xffffffff)
		N_col = N_row;

	for (wsUint i=0 ; i<N_row ; i++) {
		wsReal *Rp_m	= Ra_mat[i];
		wsUint N_med	= 0;
#ifdef USE_SSE
		N_med = getMed(N_col);
		sse_t sse_sum = sseSet(0.0);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_m = (sse_t *)(Rp_m + j);
			sse_sum = sseAdd(sse_sum, sseMul(*sse_m, *sse_m));
		}
		sseSum(sse_sum, Ra_ret[i]);
#endif
		for (wsUint j=N_med ; j<N_col ; j++)
			Ra_ret[i] += SQR(Rp_m[j]);
		R_sum += Ra_ret[i];
	}
	if (Rp_sumAll) *Rp_sumAll = R_sum;

	return Ra_ret;
}

/*
 *
 * Col-sum
 * 
 */

wsVec sseMsumC(wsReal **Ra_mat, wsUint N_row, wsUint N_col,
	wsReal *Rp_sumAll/*=NULL*/)
{
	wsVec Ra_ret	= NULL;
	if (N_row) {
		sseMalloc(Ra_ret, wsReal, N_col);
	} else {
		sseCalloc(Ra_ret, wsReal, N_col);
		if (Rp_sumAll) *Rp_sumAll = W0;
		return Ra_ret;
	}

	/* Copy first row */
	memcpy(Ra_ret, Ra_mat[0], sizeof(wsReal)*N_col);

	/* Across each row, sum up */
	for (wsUint i=1 ; i<N_row ; i++)
		sseVaV(Ra_ret, Ra_mat[i], Ra_ret, N_col);
	if (Rp_sumAll) *Rp_sumAll = sseVsum(Ra_ret, N_col);

	return Ra_ret;
}

wsReal* sseMasumC(wsReal **Ra_mat, wsUint N_row, wsUint N_col,
	wsReal *Rp_sumAll/*=NULL*/)
{
	wsReal *Ra_ret	= NULL;
	if (N_row) {
		sseMalloc(Ra_ret, wsReal, N_col);
	} else {
		sseCalloc(Ra_ret, wsReal, N_col);
		if (Rp_sumAll) *Rp_sumAll = W0;
		return Ra_ret;
	}

	/* Copy first row */
	memcpy(Ra_ret, Ra_mat[0], sizeof(wsReal)*N_col);

	/* Across each row, sum up */
	for (wsUint i=1 ; i<N_row ; i++)
		sseAVaV(Ra_ret, Ra_mat[i], Ra_ret, N_col);
	if (Rp_sumAll) *Rp_sumAll = sseVsum(Ra_ret, N_col);

	return Ra_ret;
}

wsReal* sseBsumC(wsMat Ra_mat, wsUint N_szBlk, wsUint N_sz,
				 wsReal *Rp_sumAll/*=NULL*/)
{
	wsUint N_cntRep = N_sz / N_szBlk;

	/* Allocate buffer */
	wsReal	*Ra_ret	= NULL;
	sseMalloc(Ra_ret, wsReal, N_sz);
	/* Copy first line */
	memcpy(Ra_ret, Ra_mat[0], sizeof(wsReal)*N_szBlk);

	/* Sum line by line */
	for (wsUint i=1 ; i<N_szBlk ; i++)
		sseVaV(Ra_ret, Ra_mat[i], Ra_ret, N_szBlk);
	/* Get all-sum */
	if (Rp_sumAll) *Rp_sumAll = sseVsum(Ra_ret, N_szBlk) * N_cntRep;
	/* Copy block-sum across entire */
	for (wsUint i=1,j=N_szBlk ; i<N_cntRep ; i++,j+=N_szBlk)
		memcpy(Ra_ret+j, Ra_ret, sizeof(wsReal)*N_szBlk);

	return Ra_ret;
}

/*
 *
 * Entire-sum
 * 
 */

wsRealCst sseMsum(wsReal **Ra_m, wsUint N_w, wsUint N_h/*=0xffffffff*/)
{
	wsReal R_ret = W0;
	if (N_h == 0xffffffff)
		N_h = N_w;

	for (wsUint i=0 ; i<N_h ; i++)
		R_ret += sseVsum(Ra_m[i], N_w);

	return R_ret;
}

wsRealCst sseMsqsum(wsReal **Ra_m, wsUint N_w, wsUint N_h/*=0xffffffff*/)
{
	wsReal R_ret = W0;
	if (N_h == 0xffffffff)
		N_h = N_w;

	for (wsUint i=0 ; i < N_h ; i++)
		R_ret += sseVsqsum(Ra_m[i], N_w);

	return R_ret;
}

wsRealCst sseBsum(wsReal **Ra_m, wsUint N_szBlk, wsUint N_rep)
{
	wsReal R_ret = sseMsum(Ra_m, N_szBlk);
	return R_ret*(wsReal)N_rep;
}

wsRealCst sseSsum(wsReal **Ra_m, wsUint N_sz)
{
	wsReal R_ret = W0;

	for (wsUint i=0 ; i<N_sz ; i++) {
		wsUint N_med = 0;
#ifdef USE_SSE
		N_med = getMed(i);
		sse_t sse_S = sseSet(0.0);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_M = (sse_t *)(Ra_m[i] + j);
			sse_S = sseAdd(sse_S, sseAdd(*sse_M, *sse_M));
		}
		sseAppendSum(sse_S, R_ret);
#endif
		for (wsUint j=N_med ; j<i ; j++)
			R_ret += Ra_m[i][j] + Ra_m[i][j];
		R_ret += Ra_m[i][i];
	}

	return R_ret;
}
/*
bm <- function(b, m, rep) {
	szblk <- dim(b)[1]
	szm <- dim(m)
	szb <- szblk * rep
	# b      m11 m12
	#   b    m21 m22
      #     b  m31 m32
	rb <- rowSums(b)
	R <- matrix(0, nrow=szb, ncol=szm[2])
	if (szm[1] == rep) {
		for (i in 1:rep) {
			for (j in 1:szblk) {
				for (k in 1:szm[2]) {
					R[(i-1)*szblk+j,k] <- rb[j]*m[i,k]
				}
			}
		}
		R
	} else {
		for (i in 1:rep) {
			s <- (i-1)*szblk + 1
			e <- i*szblk
			cat(s, "~", e,"\n")
			for (j in 1:szm[2]) {
				R[s:e,j] <- b%*%m[s:e,j]
			}
		}
		R
	}
}

*/
wsMat sseBpM(wsMat Rp_bMat, wsUint N_szBlk, wsUint N_rep,
			 wsMat Rp_mat, wsUint N_r, wsUint N_c)
{
	wsMat Ra_ret = sseMatrix(N_szBlk*N_rep, N_c);
	wsUint N_s = N_szBlk*N_rep;

	if (N_s == N_r) {
		/* Have same dimension, just multiply it */
/*
		for (i in 1:rep) {
			s <- (i-1)*szblk + 1
			e <- i*szblk
			cat(s, "~", e,"\n")
			for (j in 1:szm[2]) {
				R[s:e,j] <- b%*%m[s:e,j]
			}
		}
*/
		for (wsUint i=0 ; i<N_rep ; i++)
			multMpM(Rp_bMat, N_s, N_s, 0, 0,
				Rp_mat, N_r-(N_r-i-1)*N_szBlk, N_c, i*N_szBlk, 0,
				Ra_ret, N_r, N_c, i*N_szBlk, 0);
	} else {
		/* Have different dimension, check N_rep == N_r */
		if (N_rep != N_r)
			halt("Incompatible dimension");
/*
		for (i in 1:rep) {
			for (j in 1:szblk) {
				for (k in 1:szm[2]) {
					R[(i-1)*szblk+j,k] <- rb[j]*m[i,k]
				}
			}
		}
*/
		wsReal *Ra_sr = sseBsumR(Rp_bMat, N_szBlk, N_szBlk*N_rep);
		for (wsUint i=0,l=0 ; i<N_rep ; i++)
			for (wsUint j=0 ; j<N_szBlk ; j++,l++)
				sseVpC(Rp_mat[i], Ra_sr[j], Ra_ret[l], N_c);
		sseFree(Ra_sr);
	}

	return Ra_ret;
}

wsMat sseBpMt(wsMat Rp_bMat, wsUint N_szBlk, wsUint N_rep,
			 wsMat Rp_mat, wsUint N_r, wsUint N_c)
{
	wsMat Ra_ret = sseMatrix(N_szBlk*N_rep, N_r);
	wsUint N_s = N_szBlk*N_rep;

	if (N_s == N_c) {
		/* Have same dimension, just multiply it */
/*
		for (i in 1:rep) {
			s <- (i-1)*szblk + 1
			e <- i*szblk
			cat(s, "~", e,"\n")
			for (j in 1:szm[2]) {
				R[s:e,j] <- b%*%m[j,s:e]
			}
		}
*/
		for (wsUint i=0 ; i<N_rep ; i++)
			multMpMt(Rp_bMat, N_s, N_s, 0, 0,
				Rp_mat, N_r, N_c-(N_c-i-1)*N_szBlk, 0, i*N_szBlk,
				Ra_ret, N_r, N_c, i*N_szBlk, 0);
	} else {
		/* Have different dimension, check N_rep == N_c */
		if (N_rep != N_c)
			halt("Incompatible dimension");
/*
		for (i in 1:rep) {
			for (j in 1:szblk) {
				for (k in 1:szm[2]) {
					R[(i-1)*szblk+j,k] <- rb[j]*m[k,i]
				}
			}
		}
*/
		wsReal *Ra_sr = sseBsumR(Rp_bMat, N_szBlk, N_szBlk*N_rep);
		for (wsUint i=0,l=0 ; i<N_rep ; i++)
			for (wsUint j=0 ; j<N_szBlk ; j++,l++)
				for (wsUint k=0 ; k<N_r ; k++)
					Ra_ret[l][k] = Ra_sr[j] * Rp_mat[k][i];
		sseFree(Ra_sr);
	}

	return Ra_ret;
}

/*
 *
 * cMask definition
 *
 */

cMask::cMask() : cVector()
{
}

cMask::cMask(wsUint N_sz, char B_on) : cVector(N_sz)
{
	init(N_sz, B_on);
}

cMask::cMask(cMask &C_ref, wsUint N_rep/*=1*/) : cVector(C_ref.size()*N_rep)
{
	wsReal *R = C_ref.get();
	for (wsUint i=0,j=0,N=C_ref.size() ; i<N ; i++,j+=N_rep)
		memset(Ra_buf+j, R[i]?0xff:0x00, sizeof(wsReal)*N_rep);
}

cMask::~cMask()
{

}

void cMask::init(wsUint N_inpSz, char B_on)
{
	if (Ra_buf == NULL) {
		N_sz = N_inpSz;
		sseMalloc(Ra_buf, wsReal, N_sz);
	}
	memset(Ra_buf, B_on?0xff:0x00, sizeof(wsReal)*N_sz);
}

void cMask::file(wsStrCst S_ext)
{
/**/cExporter* E = cExporter::summon(S_ext);

	for (wsUint i=0,N=size() ; i<N ; i++)
		E->fmt("%g	", Ra_buf[i]?1.0:0.0);
	delete E;
}

void cVector::rem()
{
	sseFree(Ra_buf);
	Ra_buf = NULL;
}

void cVector::set0()
{
	if (Ra_buf)
		memset(Ra_buf, 0x00, sizeof(wsReal)*size());
}

cVector& cVector::inv()
{
	wsReal *Ra_d = get();
	wsUint N_med = 0;
#ifdef USE_SSE
	N_med = getMed(size());
	sse_t sse_1 = sseSet(1.0);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp)
		*((sse_t *)(Ra_d + i)) = sseDiv(sse_1, *((sse_t *)(Ra_d + i)));
#endif
	/* Rest part */
	for (wsUint i=N_med ; i<size() ; i++)
		Ra_d[i] = W1 / Ra_d[i];
		
	return *this;
}

cVector cVector::p1p()
{
	wsUint	N	= size();
	wsVec	R	= sseVector(N);
	sseVp1p(get(), N, R);
	return cVector(R, N);
}

cVector cVector::p21p2()
{
	wsUint	N	= size();
	wsVec	R	= sseVector(N);
	sseVp1p(get(), N, R, 1);
	return cVector(R, N);
}

/* Converts vector to matrix */
cStdMatrix cVector::toMt()
{
	wsReal **Ra_ret = sseMatrix(1, size());
	memcpy(Ra_ret[0], get(), sizeof(wsReal)*size());
	return cStdMatrix(1, size(), Ra_ret);
}

wsReal* famEIGENDECOMPOSITION(wsSym Ra_S, wsUint N_sz, wsUint *Na_mask,
	wsMat *Rp_eVec/*=NULL*/, char B_noTranspose/*=0*/)
{
/**/char *Ba_mask = NULL;
	wsAlloc(Ba_mask, char, N_sz);

	/* Copy mask */
/**/wsUint *Na_c = NULL;
	wsAlloc(Na_c, wsUint, N_sz);
	memcpy(Na_c, Na_mask, sizeof(wsUint)*N_sz);

	/* Eigenvalues vector */
/**/wsReal *Ra_allEval = sseVector(N_sz);
/**/wsUint *Na_indices = NULL;
	wsAlloc(Na_indices, wsUint, N_sz);

	/* Walking matrix */
	wsUint N_walkedLoop = 0;
	wsUint N_idx = 0;
	wsUint N_insIdx = 0;
	vector<wsMat>	Rv_mats;
	vInt			Nv_sizes;
	do {
		vInt Nv_indices;
		/* Clean mask */
		memset(Ba_mask, 0x00, sizeof(char)*N_sz);

		/* Skip -1 elements */
		for ( ; N_idx<N_sz ; N_idx++) if (Na_c[N_idx] != 0xffffffff) break;
		/* Stop if no element left */
		if (N_idx == N_sz) break;

		wsUint N_v = Na_c[N_idx];
		wsUint J = 0;
		/* Record it and find sames */
		for (wsUint i=N_idx ; i<N_sz ; i++)
			if (N_v == Na_c[i]) {
				Nv_indices.push_back(Na_c[i]);
				Na_c[i] = 0xffffffff;
				Ba_mask[i] = 1;
				Na_indices[i] = (N_walkedLoop << 16) | (J++);
			}
		/* Resize matrix */
		wsUint N_sub = (wsUint)Nv_indices.size();
/**/	wsSym Ra_sub = sseSsubsetRect(Ra_S, N_sz, Ba_mask, 1, N_sub);

		/* Do eigendecomposition with NO sort */
/**/	wsMat	Ra_subEvec	= NULL;
		wsReal*	Ra_subEval	= NULL;
		if (Rp_eVec)
/**/		Ra_subEval = EIGENDECOMPOSITION(Ra_sub, N_sub, &Ra_subEvec,
				B_noTranspose, 1);
		else
			Ra_subEval = EIGENDECOMPOSITION(Ra_sub, N_sub, NULL,
				B_noTranspose, 1);
		sseUnmat(Ra_sub, N_sub);
		/* Insert eigenvectors */
		if (Rp_eVec)
			Rv_mats.push_back(Ra_subEvec);
		Nv_sizes.push_back(N_sub);
		/* Copy eigenvalues */
		memcpy(Ra_allEval+N_insIdx, Ra_subEval, sizeof(wsReal)*N_sub);
		sseFree(Ra_subEval);
		N_insIdx += N_sub;
		N_walkedLoop++;
	} while (N_idx < N_sz);
	DEALLOC(Ba_mask);

	/* Sort eigenvalues */
/**/xRealSort *Xa_r = buildRealSort(Ra_allEval, N_sz);
	qsort(Xa_r, N_sz, sizeof(xRealSort), sort_real);
	/* Re-arrange indices */
	wsUint *Na_newIndices = NULL;
	if (Rp_eVec) {
/**/	wsAlloc(Na_newIndices, wsUint, N_sz);
		for (wsUint i=0 ; i<N_sz ; i++)
			Na_newIndices[i] = Na_indices[Xa_r[i].i];
	}
	/* Re-arrange eigenvalues */
	for (wsUint i=0 ; i<N_sz ; i++)
		Ra_allEval[i] = Xa_r[i].V;
	DEALLOC(Xa_r);
	DEALLOC(Na_indices);

	/* Make eVec */
	if (Rp_eVec) {
/**/	wsMat Ra_allEvec = NULL;
		Ra_allEvec = sseEmptyMatrix(N_sz, N_sz);

		/* Re-arrange eigenvectors by the order of eigenvalues */
		memcpy(Na_c, Na_mask, sizeof(wsUint)*N_sz);
		wsUint N_insStart = 0;
		/* Walking matrix */
		for (wsUint K=0 ; K<(wsUint)Nv_sizes.size() ; K++) {
			wsUint *Na_insIndices = NULL;
/**/		wsAlloc(Na_insIndices, wsUint, Nv_sizes[K]);
			/* Clean mask */

			/* Find i(th)s */
			for (wsUint i=0 ; i<N_sz ; i++)
				if ((Na_newIndices[i]>>16) == K) {
					if ((wsUint)Nv_sizes[K] <= (Na_newIndices[i]&0xffff))
						halt("ERR");
					Na_insIndices[Na_newIndices[i]&0xffff] = i;
				}

			/* Mapping */
			wsMat Ra_cur = Rv_mats[K];
			for (int i=0 ; i<Nv_sizes[K] ; i++) {
				wsUint Ni = Na_insIndices[i];
				for (int j=0 ; j<Nv_sizes[K] ; j++) {
//					wsUint Nj = Na_insIndices[j];
					Ra_allEvec[Ni][N_insStart+j] = Ra_cur[i][j];
				}
			}
			sseUnmat(Ra_cur, Nv_sizes[K]);

			DEALLOC(Na_insIndices);
			N_insStart += Nv_sizes[K];
		}
		DEALLOC(Na_newIndices);
		*Rp_eVec = Ra_allEvec;
	}
	DEALLOC(Na_c);

	return Ra_allEval;
}

wsSym sseSwalk(wsSym Ra_S, wsUint N_sz, short *Na_mask, SYMMATWALKPROC H_proc)
{
	char *Ba_mask = NULL;
	wsAlloc(Ba_mask, char, N_sz);

	/* Copy mask */
	short *Na_c = NULL;
	wsAlloc(Na_c, short, N_sz);
	memcpy(Na_c, Na_mask, sizeof(short)*N_sz);

	/* Making return matrix */
	wsSym Ra_ret = sseEmptySymMat(N_sz);

	/* Walking matrix */
//	wsUint N_walked = 0;
	wsUint N_idx = 0;
	do {
		vInt Nv_indices;
		/* Clean mask */
		memset(Ba_mask, 0x00, sizeof(char)*N_sz);

		/* Skip -1 elements */
		for ( ; N_idx<N_sz ; N_idx++) if (Na_c[N_idx] != -1) break;
		/* Stop if no element left */
		if (N_idx == N_sz) break;

		short N_v = Na_c[N_idx];
		/* Record it and find sames */
		for (wsUint i=N_idx ; i<N_sz ; i++)
			if (N_v == Na_c[i]) {
				Nv_indices.push_back(Na_c[N_idx]);
				Na_c[N_idx] = -1;
				Ba_mask[N_idx] = 1;
			}
		/* Resize matrix */
		wsUint N_sub = (wsUint)Nv_indices.size();
		wsSym Ra_sub = sseSsubsetRect(Ra_S, N_sz, Ba_mask, 1, N_sub);

		/* Do transform */
		wsSym Ra_sret = H_proc(Ra_sub, N_sub);

		/* Unmapping */
		for (wsUint i=0 ; i<N_sub ; i++) {
			wsUint Ni = Nv_indices[i];
			for (wsUint j=0 ; j<=i ; j++)
				Ra_ret[Ni][Nv_indices[j]] = Ra_sret[i][j];
		}
	} while (N_idx < N_sz);

	return Ra_ret;
}

wsVec sseSdiag(wsSym Ra_mat, wsUint N_sz)
{
	return sseMdiag(Ra_mat, N_sz, N_sz);
}

wsVec sseMdiag(wsMat Ra_mat, wsUint N_r, wsUint N_c)
{
	wsUint	N_l		= N_r > N_c ? N_c : N_r;
	wsVec	Ra_r	= sseVector(N_l);
	LOOP (i, N_l)
		Ra_r[i] = Ra_mat[i][i];
	return Ra_r;
}

void SP1(wsSym Lmat, wsUint nr, wsReal edgecut=0.3, wsReal mcut=0.3)
{
	//SP1=function(Lmat,edgecut=0.3,mcut=0.3){
	//nr=dim(Lmat)[1]  #modified
	//rv=NULL

	//SYM_t WO = sseSsq(Lmat, nr);
	//WO<-Lmat^2
// 		I=diag(nr)
// 		WO<- WO-I
// 		W=WO
// 		W[W<edgecut]<-0
// 		cmpl=split(Lmat,edgecut)
// 
// 		done<-0
// 		binvector<-rep(0,dim(Lmat)[1])
// 		for (n in 1:max(cmpl)){
// 			cmp<-which(cmpl==n)
// 				if(length(cmp)==1){
// 					done<-max(binvector)+1
// 						binvector[cmp]<-done
// 				}else{
// 					rW<-W[cmp,cmp]
// 					rs<-rowSums(rW)
// 						D<-diag(rs)
// 						iD=solve(D)
// 						P=iD%*%rW
// 						Qlist=NULL
// 						taken=NULL
// 						l=length(cmp)-1
// 						ev<-eigen(P)$values
// 						sev<-sort(ev,decreasing=TRUE)
// 						m<-match(sev,ev)
// 						peV<-eigen(P)$vectors
// 						eV<-peV[,m]
// 					modieV<-eV[,1:length(cmp)]
// 					for (k in 1:l){
// 						if(k==1){Qith=1
// 						}else{
// 							Uk<-modieV[,2:k,drop=F]
// 							cl <- kmeans(Uk,k,iter.max=200)
// 								cv<-cl$cluster
// 								Qlist=NULL
// 								Qith=0
// 								for (j in 1:k){
// 									A=c(1:length(cv))[cv==j]
// 									inA<-sum(rW[A,A])/sum(rW)
// 										itAV<-sum(rW[A,])/sum(rW)
// 										Qjth=inA-(itAV*itAV)
// 										Qith=Qith+Qjth
// 								}
// 						}	
// 						Qlist=c(taken,Qith)
// 							taken=Qlist
// 					}
// 					Qvl=min(Qlist[Qlist>mcut])
// 						fcz<-which(Qlist==Qvl)
// 						if(Qvl!=1){
// 							Uk<-modieV[,2:fcz,drop=F]
// 							FinalCL <- kmeans(Uk, fcz,iter.max=200)
// 								pbinvt<-FinalCL$cluster
// 								binvt<-pbinvt+max(binvector)
// 						}else{
// 							binvt<-max(binvector)+1
// 						}
// 						binvector[cmp]<-binvt
// 							done<-max(binvector)
// 				}
// 		}
// 
// 		return(binvector)
}

wsMat sseMMM(wsMat Rp_mat1, wsUint N_r1, wsUint N_c1,
			 wsMat Rp_m2,  wsUint N_r2, wsUint N_c2,
			 wsMat Rp_m3,  wsUint N_r3, wsUint N_c3)
{
	wsUint i, j, k;
	wsMat	Ra_ret		= sseEmptyMatrix(N_r1, N_c3);
	/* Create intermediate buffer */
	wsVec	Ra_interm	= sseVector(N_c2*2);

	for (i=0 ; i<N_r1 ; i++) {
		wsVec	Rp_matR	= Ra_ret[i];
		wsVec	Rp_m1	= Rp_mat1[i];
		wsUint	N_med	= getMed(N_c2);

		memset(Ra_interm, 0x00, sizeof(wsReal)*N_c2);

		/* Calculate A[i,]*B */
		for (j=0 ; j<N_c1 ; j++) {
			sse_t sse_m1_il = sseSet(Rp_m1[j]);

			/* SSE-enabled part */
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_m2_lj	= (sse_t *)(&(Rp_m2[j][k]));
				sse_t *sse_mR		= (sse_t *)(&(Ra_interm[k]));
				*sse_mR = sseAdd(*sse_mR, sseMul(sse_m1_il, *sse_m2_lj));
			}
			/* Rest part */
			for (k=N_med ; k<N_c2 ; k++)
				Ra_interm[k] += Rp_m1[j] * Rp_m2[j][k];
		}

		if (N_c3 < sseJmp) {
			for (j=0 ; j<N_c3 ; j++) {
				wsReal R_sum = 0.0f;

				for (k=0 ; k<N_c2 ; k++)
					R_sum += Ra_interm[k] * Rp_m3[k][j];

				Rp_matR[j] = R_sum;
			}
		} else {
			N_med = getMed(N_c3);
			for (j=0 ; j<N_c2 ; j++) {
				wsRealCst	R_m		= Ra_interm[j];
				sse_t	sse_m	= sseSet(R_m);

				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_m2_lj	= (sse_t *)(&(Rp_m3[j][k]));
					sse_t *sse_mR		= (sse_t *)(&(Rp_matR[k]));
					*sse_mR = sseAdd(*sse_mR, sseMul(sse_m, *sse_m2_lj));
				}
				for (k=N_med ; k<N_c3 ; k++)
					Rp_matR[k] += R_m * Rp_m3[j][k];
			}
		}
	}

	sseFree(Ra_interm);
	return Ra_ret;
}

wsMat sseMSM(wsMat Rp_mat1, wsUint N_r1, wsUint N_c1,
			 wsMat Ra_s,  wsUint N_sz,
			 wsMat Rp_m3,  wsUint N_r3, wsUint N_c3)
{
	/* Size check */
	if (N_c1 != N_sz) halt_fmt(WISARD_SYST_INVL_DIM, N_c1, N_sz);

	/* Memory allocation */
	wsVec	Ra_interm	= sseVector(N_c1*2);
	wsMat	Ra_ret		= sseMatrix(N_r1, N_sz);

	/* Computation */
#ifdef USE_SYM
	LOOP(i, N_r1) {
		memset(Ra_interm, 0x00, sizeof(wsReal)*N_sz);

		LOOP(a, N_sz) {
			wsReal *Ra_m = Rp_mat1[i];
			register wsReal R_a = Ra_m[a];
			wsReal *Ra_x = Ra_s[a];

#ifdef USE_SSE
			wsUint	N_med = getMed(a);
			sse_t	sse_a = sseSet(R_a);
			sse_t	sse_ra = sseSet(0.0);

			SSE_LOOP(b, N_med) {
				sse_t *sse_r = (sse_t *)(Ra_interm + b);
				sse_t *sse_s = (sse_t *)(Ra_x + b);
				sse_t *sse_m = (sse_t *)(Ra_m + b);
				*sse_r = sseAdd(*sse_r, sseMul(sse_a, *sse_s));
				sse_ra = sseAdd(sse_ra, sseMul(*sse_m, *sse_s));
			}
			sseAppendSum(sse_ra, Ra_interm[a]);

			LOOPV(b, N_med, a) {
				Ra_interm[a] += Ra_m[b] * Ra_x[b];
				Ra_interm[b] += R_a * Ra_x[b];
			}
#else
			LOOP(b, a) {
				Ra_interm[a] += Ra_m[b] * Ra_x[b];
				Ra_interm[b] += R_a * Ra_x[b];
			}
#endif
			Ra_interm[a] += R_a * Ra_s[a][a];
		}

		Ra_ret[i] = sseVpM(Ra_interm, N_sz, Rp_m3, N_r3, N_c3);
	}

	sseFree(Ra_interm);
	return Ra_ret;
#else
#	error "ERROR"
#endif
}

wsMat sseMDM(wsMat Ra_m1, wsUint N_r1, wsUint N_c1,
	wsDiag Ra_diag, wsMat Ra_m3,  wsUint N_r3, wsUint N_c3)
{
	/* Size check */
	if (N_c1 != N_r3) halt_fmt(WISARD_SYST_INVL_DIM, "MDM", N_c1, N_r3);

	/* Memory allocation */
	wsMat	Ra_ret		= sseEmptyMatrix(N_r1, N_c3);
	wsVec	Ra_interm	= sseVector(N_c1*2);

	/* Compute */
	LOOP(i, N_r1) {
		wsReal	*Rp_matR	= Ra_ret[i];
		wsReal	*Rp_m1		= Ra_m1[i];
		wsUint	N_med		= getMed(N_c1);

		memset(Ra_interm, 0x00, sizeof(wsReal)*N_c1);

		/* Calculate A[i,]*B */
		sseVpV(Rp_m1, Ra_diag, Ra_interm, N_c1);

		if (N_c3 < sseJmp) {
			LOOP(j, N_c3) {
				wsReal R_sum = W0;

				LOOP(k, N_c1)
					R_sum += Ra_interm[k] * Ra_m3[k][j];

				Rp_matR[j] = R_sum;
			}
		} else {
#ifdef USE_SSE
			N_med = getMed(N_c3);
#endif
			LOOP(j, N_c1) {
				wsRealCst	R_m		= Ra_interm[j];
#ifdef USE_SSE
				sse_t	sse_m	= sseSet(R_m);

				SSE_LOOP(k, N_med) {
					sse_t *sse_m2_lj	= (sse_t *)(&(Ra_m3[j][k]));
					sse_t *sse_mR		= (sse_t *)(&(Rp_matR[k]));
					*sse_mR = sseAdd(*sse_mR, sseMul(sse_m, *sse_m2_lj));
				}
#endif
				LOOPV(k, N_med, N_c3)
					Rp_matR[k] += R_m * Ra_m3[j][k];
			}
		}
	}

	sseFree(Ra_interm);
	return Ra_ret;
}

wsMat sseMMMt(wsMat Ra_m1, wsUint N_r1, wsUint N_c1, wsMat Ra_m2,
	wsUint N_r2, wsUint N_c2, wsMat Ra_m3, wsUint N_r3, wsUint N_c3)
{
	wsUint i, j, k;

	/* Memory allocation */
	wsMat	Ra_ret		= sseEmptyMatrix(N_r1, N_r3);
	wsVec	Ra_interm	= sseVector(N_c2*2);

	/* In case of third matrix is omitted */
	if (Ra_m3 == NULL) {
		Ra_m3 = Ra_m1;
		N_r3 = N_r1;
		N_c3 = N_c1;
	}

	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_matR	= Ra_ret[i];
		wsReal	*Rp_m1		= Ra_m1[i];
		wsUint	N_med		= getMed(N_c2);

		memset(Ra_interm, 0x00, sizeof(wsReal)*N_c2);

		/* Calculate A[i,]*B */
		for (j=0 ; j<N_c1 ; j++) {
			sse_t sse_m1_il = sseSet(Rp_m1[j]);

			/* SSE-enabled part */
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_m2_lj	= (sse_t *)(&(Ra_m2[j][k]));
				sse_t *sse_mR		= (sse_t *)(&(Ra_interm[k]));
				*sse_mR = sseAdd(*sse_mR, sseMul(sse_m1_il, *sse_m2_lj));
			}
			/* Rest part */
			for (k=N_med ; k<N_c2 ; k++)
				Ra_interm[k] += Rp_m1[j] * Ra_m2[j][k];
		}

		N_med = getMed(N_c3);
		for (j=0 ; j<N_r3 ; j++) {
			if (N_med) {
				sse_t sse_R = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_interm	= (sse_t *)(Ra_interm + k);
					sse_t *sse_m2		= (sse_t *)(Ra_m3[j] + k);
					sse_R = sseAdd(sse_R, sseMul(*sse_interm, *sse_m2));
				}
				sseSum(sse_R, Rp_matR[j]);
			} else
				Rp_matR[j] = W0;
			for (k=N_med ; k<N_c3 ; k++)
				Rp_matR[j] += Ra_interm[k] * Ra_m3[j][k];
		}
	}

	sseFree(Ra_interm);
	return Ra_ret;
}

wsSym sseMMMt(wsReal **Rp_mat1, wsUint N_r1, wsUint N_c1,
			  wsReal **Rp_m2,  wsUint N_r2, wsUint N_c2)
{
	wsUint i, j, k;
	wsSym Ra_ret	= sseEmptySymMat(N_r1);
	/* Create intermediate buffer */
	wsReal	*Ra_interm = NULL;
	sseMalloc(Ra_interm, wsReal, N_c2*2);

	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_matR	= Ra_ret[i];
		wsReal	*Rp_m1		= Rp_mat1[i];
		wsUint	N_med		= getMed(N_c2);

		memset(Ra_interm, 0x00, sizeof(wsReal)*N_c2);

		/* Calculate A[i,]*B */
		for (j=0 ; j<N_c1 ; j++) {
			sse_t sse_m1_il = sseSet(Rp_m1[j]);

			/* SSE-enabled part */
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_m2_lj	= (sse_t *)(&(Rp_m2[j][k]));
				sse_t *sse_mR		= (sse_t *)(&(Ra_interm[k]));
				*sse_mR = sseAdd(*sse_mR, sseMul(sse_m1_il, *sse_m2_lj));
			}
			/* Rest part */
			for (k=N_med ; k<N_c2 ; k++)
				Ra_interm[k] += Rp_m1[j] * Rp_m2[j][k];
		}

		N_med = getMed(N_c1);
		for (j=0 ; j<=i ; j++) {
			wsReal *Rp_m1 = Rp_mat1[j];
			if (N_med) {
				sse_t sse_R = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_interm	= (sse_t *)(Ra_interm + k);
					sse_t *sse_m2		= (sse_t *)(Rp_m1 + k);
					sse_R = sseAdd(sse_R, sseMul(*sse_interm, *sse_m2));
				}
				sseSum(sse_R, Rp_matR[j]);
			} else
				Rp_matR[j] = W0;
			for (k=N_med ; k<N_c1 ; k++)
				Rp_matR[j] += Ra_interm[k] * Rp_m1[k];
		}
	}

	sseFree(Ra_interm);
	return Ra_ret;
}

wsMat sseMDMt(wsReal **Rp_mat1, wsUint N_r1, wsUint N_c1,
			  wsReal *Rp_diag, wsReal **Rp_m3, wsUint N_r3, wsUint N_c3)
{
	wsUint i, j, k;

	/* Create intermediate buffer */
	wsReal	*Ra_interm = NULL;
	sseMalloc(Ra_interm, wsReal, N_c1*2);

	/* In case of third matrix is omitted */
	if (Rp_m3 == NULL) {
		Rp_m3 = Rp_mat1;
		N_r3 = N_r1;
		N_c3 = N_c1;
	}
	if (N_c1 != N_c3)
		halt_fmt(WISARD_SYST_INVL_DIM, N_c1, N_c3);
	wsMat Ra_ret	= sseEmptyMatrix(N_r1, N_r3);

	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_matR	= Ra_ret[i];
		wsReal	*Rp_m1		= Rp_mat1[i];
		wsUint	N_med		= getMed(N_c1);

		memset(Ra_interm, 0x00, sizeof(wsReal)*N_c1);

		/* Calculate A[i,]*B */
		sseVpV(Rp_m1, Rp_diag, Ra_interm, N_c1);

		N_med = getMed(N_c3);
		for (j=0 ; j<N_r3 ; j++) {
			if (N_med) {
				sse_t sse_R = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_interm	= (sse_t *)(Ra_interm + k);
					sse_t *sse_m2		= (sse_t *)(Rp_m3[j] + k);
					sse_R = sseAdd(sse_R, sseMul(*sse_interm, *sse_m2));
				}
				sseSum(sse_R, Rp_matR[j]);
			} else
				Rp_matR[j] = W0;
			for (k=N_med ; k<N_c3 ; k++)
				Rp_matR[j] += Ra_interm[k] * Rp_m3[j][k];
		}
	}

	sseFree(Ra_interm);
	return Ra_ret;
}

wsSym sseMDMt(wsReal **Rp_mat1, wsUint N_r1, wsUint N_c1, wsReal *Rp_diag)
{
	wsUint i, j, k;

	/* Create intermediate buffer */
	wsReal	*Ra_interm = NULL;
	sseMalloc(Ra_interm, wsReal, N_c1*2);

	/* In case of third matrix is omitted */
	wsSym	Ra_ret	= sseEmptySymMat(N_r1);

	for (i=0 ; i<N_r1 ; i++) {
		wsReal	*Rp_matR	= Ra_ret[i];
		wsReal	*Rp_m1		= Rp_mat1[i];
		wsUint	N_med		= getMed(N_c1);

		memset(Ra_interm, 0x00, sizeof(wsReal)*N_c1);

		/* Calculate A[i,]*B */
		sseVpV(Rp_m1, Rp_diag, Ra_interm, N_c1);

		N_med = getMed(N_c1);
		for (j=0 ; j<=i ; j++) {
			if (N_med) {
				sse_t sse_R = sseSet(0.0);
				for (k=0 ; k<N_med ; k+=sseJmp) {
					sse_t *sse_interm	= (sse_t *)(Ra_interm + k);
					sse_t *sse_m2		= (sse_t *)(Rp_mat1[j] + k);
					sse_R = sseAdd(sse_R, sseMul(*sse_interm, *sse_m2));
				}
				sseSum(sse_R, Rp_matR[j]);
			} else
				Rp_matR[j] = W0;
			for (k=N_med ; k<N_c1 ; k++)
				Rp_matR[j] += Ra_interm[k] * Rp_mat1[j][k];
		}
	}

	sseFree(Ra_interm);
	return Ra_ret;
}

wsVec multMV(wsMat Rp_mat, wsUintCst N_row, wsUintCst N_col, wsVecCst Rp_v)
{
	/* Make return matrix */
	wsVec Ra_ret = sseVector(N_row);

	for (wsUint i=0 ; i<N_row ; i++) {
		wsReal R_sum = W0;
		for (wsUint j=0 ; j<N_col ; j++)
			R_sum += Rp_mat[i][j] * Rp_v[j];
		Ra_ret[i] = R_sum;
	}

	return Ra_ret;
}

wsReal multVMV(wsVecCst Ra_v, wsMatCst Ra_m, wsUintCst N_len, wsVecCst Ra_v2)
{
	wsReal R_res = W0;
	if (Ra_v2 == NULL) LOOP(i, N_len) {
		wsReal R_sub = W0;
		for (wsUint j=0 ; j<N_len ; j++)
			R_sub += Ra_v[j] * Ra_m[j][i];
		//printf("%g + %g * %g\n", R_res, R_sub, Ra_v2[i]);
		R_res += R_sub*Ra_v2[i];
	} else  LOOP (i, N_len) {
		wsReal R_sub = W0;
		for (wsUint j=0 ; j<N_len ; j++)
			R_sub += Ra_v[j] * Ra_m[j][i];
		//printf("%g + %g * %g\n", R_res, R_sub, Ra_v2[i]);
		R_res += R_sub*Ra_v2[i];
	}

	return R_res;
}

wsMat multMtM(wsMatCst Rp_mat, wsUint N_row, wsUint N_col)
{
	wsMat Ra_ret = NULL;
	wsAlloc(Ra_ret, wsReal*, N_col);

	LOOP (i, N_col) {
		sseCalloc(Ra_ret[i], wsReal, N_col);

		LOOP (j, N_col) LOOP (k, N_row)
			Ra_ret[i][j] += Rp_mat[k][i] * Rp_mat[k][j];
	}

	return Ra_ret;
}

wsMat multMMt(wsMatCst Rp_mat, wsUint N_row, wsUint N_col)
{
	wsReal **Ra_ret = NULL;
	wsAlloc(Ra_ret, wsReal*, N_row);

	LOOP (i, N_row) {
		wsCalloc(Ra_ret[i], wsReal, N_row);
		LOOP (j, N_row) LOOP (k, N_col)
			Ra_ret[i][j] += Rp_mat[i][k] * Rp_mat[j][k];
	}

	return Ra_ret;
}

wsReal sseVpMpV(wsVec Ra_v, wsMat Ra_m, wsUint N_len, wsVecCst Ra_v2/*=NULL*/,
				wsReal *Ra_mask/*=NULL*/)
{
	if (Ra_v2 == NULL)
		Ra_v2 = Ra_v;
	wsReal	R_ret = W0;
	wsUint	j, k;

	if (Ra_mask == NULL) {
#ifdef USE_SSE
		wsUint	N_med		= getMed(N_len);
		wsVec	Ra_interm	= NULL;
		sseCalloc(Ra_interm, wsReal, N_len);

		for (j=0 ; j<N_len ; j++) {
			sse_t sse_v = sseSet(Ra_v[j]);

			/* SSE-enabled part */
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_m_jk	= (sse_t *)(&(Ra_m[j][k]));
				sse_t *sse_mR	= (sse_t *)(&(Ra_interm[k]));
				*sse_mR = sseAdd(*sse_mR, sseMul(sse_v, *sse_m_jk));
			}
			/* Rest part */
			for (k=N_med ; k<N_len ; k++)
				Ra_interm[k] += Ra_v[j] * Ra_m[j][k];
		}
		sse_t sse_sum = sseSet(0.0f);
		for (j=0 ; j<N_med ; j+=sseJmp) {
			sse_t	*sse_mR	= (sse_t *)(&(Ra_interm[j]));
			sse_t	*sse_mK	= (sse_t *)(&(Ra_v2[j]));
			sse_sum = sseAdd(sse_sum, sseMul(*sse_mR, *sse_mK));
		}
		sseSum(sse_sum, R_ret);
		for (j=N_med ; j<N_len ; j++)
			R_ret += Ra_interm[j]*Ra_v2[j];

		sseFree(Ra_interm);
#else
		for (wsUint i=0 ; i<N_len ; i++) {
			wsReal R_sub = W0;
			for (wsUint j=0 ; j<N_len ; j++)
				R_sub += Ra_v[j] * Ra_m[j][i];
			//printf("%g + %g * %g\n", R_res, R_sub, Ra_v2[i]);
			R_ret += R_sub*Ra_v2[i];
		}
#endif
	} else {
#ifdef USE_SSE
		wsUint	N_med		= getMed(N_len);
		wsVec	Ra_interm	= NULL;
		sseCalloc(Ra_interm, wsReal, N_len);

		for (j=0 ; j<N_len ; j++) {
			sse_t sse_v = sseSet(Ra_v[j]);

			/* SSE-enabled part */
			for (k=0 ; k<N_med ; k+=sseJmp) {
				sse_t *sse_m_jk	= (sse_t *)(&(Ra_m[j][k]));
				sse_t *sse_mR	= (sse_t *)(&(Ra_interm[k]));
				sse_t *sse_M	= (sse_t *)(Ra_mask + k);
				*sse_mR = sseAdd(*sse_mR,
					sseAnd(*sse_M, sseMul(sse_v, *sse_m_jk))
					);
			}
			/* Rest part */
			for (k=N_med ; k<N_len ; k++)
				if (Ra_mask[k])
					Ra_interm[k] += Ra_v[j] * Ra_m[j][k];
		}
		sse_t sse_sum = sseSet(0.0f);
		for (j=0 ; j<N_med ; j+=sseJmp) {
			sse_t*	sse_mR	= (sse_t *)(&(Ra_interm[j]));
			sse_t*	sse_mK	= (sse_t *)(&(Ra_v2[j]));
			sse_t*	sse_M	= (sse_t *)(Ra_mask + j);
			sse_sum = sseAdd(sse_sum, 
				sseAnd(*sse_M, sseMul(*sse_mR, *sse_mK))
				);
		}
		sseSum(sse_sum, R_ret);
		for (j=N_med ; j<N_len ; j++)
			if (Ra_mask[j])
				R_ret += Ra_interm[j]*Ra_v2[j];

		sseFree(Ra_interm);
#else
		for (wsUint i=0 ; i<N_len ; i++) {
			wsReal R_sub = W0;
			for (wsUint j=0 ; j<N_len ; j++)
				if (Ra_mask[j])
					R_sub += Ra_v[j] * Ra_m[j][i];
			//printf("%g + %g * %g\n", R_res, R_sub, Ra_v2[i]);
			if (Ra_mask[i])
				R_ret += R_sub*Ra_v2[i];
		}
#endif
	}

	return R_ret;
}

wsReal sseVpSpV(wsVec Ra_V, wsSym Ra_S, wsUintCst N_sz, wsVec Ra_W/*=NULL*/,
	wsReal *Ra_mask/*=NULL*/)
{
	if (Ra_W == NULL) Ra_W = Ra_V;

	wsReal *Ra_vM = sseVpS(Ra_V, N_sz, Ra_S, Ra_mask);
	wsReal R;
	if (Ra_mask != NULL)
		R = sseVVsel(Ra_mask, Ra_vM, N_sz, Ra_W);
	else
		R = sseVV(Ra_vM, N_sz, Ra_W);
	sseFree(Ra_vM);

	return R;
}

wsReal sseVVmiss(wsRealCst *Ra_v1, wsUint N_sz, wsRealCst *Ra_v2/*=NULL*/)
{
	wsReal R_ret = W0;
	if (Ra_v2 == NULL)
		Ra_v2 = Ra_v1;

#ifdef USE_SSE
	wsUint N_med = getMed(N_sz);
	sse_t sse_ret = sseSet(0.0);
	sse_t sse_M = sseSet(WISARD_NA_REAL);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
		sse_t sse_m1 = sseAnd(*sse_v1, sseNeq(*sse_v1, sse_M));
		sse_t sse_m2 = sseAnd(*sse_v2, sseNeq(*sse_v2, sse_M));

		sse_ret = sseAdd(sse_ret, sseMul(sse_m1, sse_m2));
	}
	sseSum(sse_ret, R_ret);
	for (wsUint i=N_med ; i<N_sz ; i++)
		if (isAvailableReal(Ra_v1[i]) && isAvailableReal(Ra_v2[i]))
			R_ret += Ra_v1[i]*Ra_v2[i];
#else
	for (wsUint i=0 ; i<N_sz ; i++) {
		if (isAvailableReal(Ra_v1[i]) && isAvailableReal(Ra_v2[i]))
			R_ret += Ra_v1[i]*Ra_v2[i];
	}
#endif

	return R_ret;
}

wsReal sseVVsel(wsReal *Ra_mask,
				wsRealCst *Ra_v1, wsUint N_sz, wsRealCst *Ra_v2/*=NULL*/)
{
	wsReal R_ret = W0;
	if (Ra_v2 == NULL)
		Ra_v2 = Ra_v1;

#ifdef USE_SSE
	wsUint N_med = getMed(N_sz);
	sse_t sse_ret = sseSet(0.0);
	for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
		sse_t *sse_v1 = (sse_t *)(Ra_v1 + i);
		sse_t *sse_v2 = (sse_t *)(Ra_v2 + i);
		sse_t *sse_mask = (sse_t *)(Ra_mask + i);
		sse_t sse_m1 = sseAnd(*sse_v1, *sse_mask);
		sse_t sse_m2 = sseAnd(*sse_v2, *sse_mask);

		sse_ret = sseAdd(sse_ret, sseMul(sse_m1, sse_m2));
	}
	sseSum(sse_ret, R_ret);
	for (wsUint i=N_med ; i<N_sz ; i++)
		if (Ra_mask[i])
			R_ret += Ra_v1[i]*Ra_v2[i];
#else
	for (wsUint i=0 ; i<N_sz ; i++) {
		if (Ra_mask[i])
			R_ret += Ra_v1[i]*Ra_v2[i];
	}
#endif

	return R_ret;
}

wsRealCst diagSum(wsReal **Ra_m, wsUint N_sz)
{
	wsReal R_ret = W0;

	for (wsUint i=0 ; i<N_sz ; i++)
		R_ret += Ra_m[i][i];

	return R_ret;
}

wsSym diagMatrix(wsReal R_val, wsUint N_size)
{
	wsReal **Ra_mat = NULL;
	wsAlloc(Ra_mat, wsReal*, N_size);
	for (wsUint i=0 ; i<N_size ; i++) {
		sseMalloc(Ra_mat[i], wsReal, i+1);

		memset(Ra_mat[i], 0x00, sizeof(wsReal)*i);
		Ra_mat[i][i] = R_val;
	}
	return Ra_mat;
}

wsReal diagSum2matrix(wsMat Ra_m1, wsUint N_r1, wsUint N_c1,
	wsMat Ra_m2, wsUint N_r2, wsUint N_c2)
{
	wsUint N_i = N_r1 > N_c2 ? N_c2 : N_r1;

	/* p = N_mat2Col
	 * m = N_mat1Col
	 * (ABC)_ii=Σ_(k=1)^p (Σ_(j=1)^m A_ij B_jk ) C_ki
	 */
	wsReal R_ret = W0;
	
	for (wsUint i=0 ; i<N_i ; i++)
		for (wsUint k=0 ; k<N_c1 ; k++)
			R_ret += Ra_m1[i][k] * Ra_m2[k][i];

	return R_ret;
}

typedef struct _xLUsolve
{
	wsUint N_row;
	wsReal **Ra_ret;
	wsReal **Ra_LU;
} xLUsolve;

int LUdist(int N_mode, int N_thread, void *Vp_thrData,
		   wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	int			*Np_idx	= (int *)Vp_thrData;
	xLUsolve	*Xp_lu	= (xLUsolve *)Vp_shareData;
	static wsUint
		N_proc	= 0;
	wsReal	R_szDiv	= (wsReal)(Xp_lu->N_row)/(wsReal)N_thread;
	int		i=0;
	wsReal	j=W0;
	int		N_ret = N_proc != 0 ? 0 : 1;

	switch (N_mode) {
	case DISTFUN_INIT:
		i = 0;
		j=W0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		if (N_ret == 0)
			return 0;
		for ( ; i<N_thread ; i++,j+=R_szDiv) {
			wsUint E = (wsUint)(j+R_szDiv+REAL_CONST(0.5));
			if (E > Xp_lu->N_row)
				E = Xp_lu->N_row;
			Np_idx[2*i] = (wsUint)(j+REAL_CONST(0.5));
			Np_idx[2*i+1] = E;
		}
		N_proc = N_ret = N_thread;
		break;
	case DISTFUN_UNINIT:
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}

	return N_ret;
}

int symLUsolve(wsReal **LU, wsUint p, wsReal x[], int n);

int LUsolve(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	xLUsolve	*Xp_lu		= (xLUsolve *)Vp_shareData;
	int			*Np_data	= (int *)Vp_data;
	wsUint		N_s			= (wsUint)Np_data[0];
	wsUint		N_e			= (wsUint)Np_data[1];

	//fprintf(stdout, "Thread %d traverse [%d, %d)\n", N_idx, N_s, N_e);
	for (wsUint i=N_s ; i<N_e ; i++)
		symLUsolve(Xp_lu->Ra_LU, i, Xp_lu->Ra_ret[i], Xp_lu->N_row);

	return 0;
}

wsReal** symLUdecomp(wsReal **Ra_mat, wsUint N_sz)
{
	wsUint i, j, k;

	//         For each row and column, k = 0, ..., n-1,
	//            find the upper triangular matrix elements for row k
	//            and if the matrix is non-singular (nonzero diagonal element).
	//            find the lower triangular matrix elements for column k. 
	//wsReal **s = sseEmptyMatrix(N_sz, N_sz);
//	wsReal *q = NULL;
//	sseMalloc(q, wsReal, N_sz);
	//int u=0;
	//cTimer t;

	/* Set H and V */
	wsReal **Ra_H = sseEmptyMatrix(N_sz, N_sz);
	wsReal **Ra_V = sseEmptyMatrix(N_sz, N_sz);

	/* Set H0 to be m_00, ... , m_0n */
	memcpy(Ra_H[0], Ra_mat[0], sizeof(wsReal)*N_sz);

	/* Set V0 to be m_00/m_00, m_10/m_00, m_20/m_00, ... m_n0/m_00 */
	wsReal D = Ra_mat[0][0];
	for (i=1 ; i<N_sz ; i++)
		Ra_V[0][i] = Ra_mat[0][i] / D; /* If symm */

	/* Calculate H1,V1 ... Hn,Vn */
	//t.start();
#ifdef USE_SSE
	wsUint N_e = getMed(N_sz);
#endif
	for (i=1 ; i<N_sz ; i++) {
		/* Calculate Hi */
		wsUint N_s = getMed(i+sseJmp-1);
		if (N_s > N_sz)
			N_s = N_sz;

		/* For V0,H0 to V(i-1),H(i-1) */
		for (j=0 ; j<i ; j++) {
			/* Calc Vj_i * Hj_(i:n)<-k */
#ifdef USE_SSE
			sse_t sse_ji = sseSet(Ra_V[j][i]);

			for (k=i ; k<N_s ; k++)
				Ra_H[i][k] -= Ra_V[j][i]*Ra_H[j][k];
			for (k=N_s ; k<N_e ; k+=sseJmp) {
				sse_t *sse_ik = (sse_t *)(Ra_H[i] + k);
				sse_t *sse_jk = (sse_t *)(Ra_H[j] + k);

				*sse_ik = sseSub(*sse_ik, sseMul(sse_ji, *sse_jk));
				//Ra_H[i][k] -= Ra_V[j][i]*Ra_H[j][k];
			}
			for (k=N_e ; k<N_sz ; k++)
				Ra_H[i][k] -= Ra_V[j][i]*Ra_H[j][k];
#else
			for (k=i ; k<N_sz ; k++)
				Ra_H[i][k] -= Ra_V[j][i]*Ra_H[j][k];
#endif
		}
		/* Settle Hi */
#ifdef USE_SSE
		for (j=i ; j<N_s ; j++)
			Ra_H[i][j] += Ra_mat[i][j];
		for (j=N_s ; j<N_e ; j+=sseJmp) {
			sse_t *sse_Hij = (sse_t *)(Ra_H[i] + j);
			sse_t *sse_Mij = (sse_t *)(Ra_mat[i] + j);
			*sse_Hij = sseAdd(*sse_Hij, *sse_Mij);
			//Ra_H[i][j] += Ra_mat[i][j];
		}
		for (j=N_e ; j<N_sz ; j++)
			Ra_H[i][j] += Ra_mat[i][j];
#else
		for (j=i ; j<N_sz ; j++)
			Ra_H[i][j] += Ra_mat[i][j];
#endif

		if (Ra_H[i][i] == W0) {
			sseUnmat(Ra_H, N_sz);
			sseUnmat(Ra_V, N_sz);
			LOG("This matrix is singular\n");	
			return NULL;
		}

		/* Calculate Vi */

		N_s = getMed(i+sseJmp);
		/* Clamping */
		if (N_s > N_sz)
			N_s = N_sz;
		/* For V0,H0 to V(i-1),H(i-1) */
		for (j=0 ; j<i ; j++) {
			/* Calc Hj_i * Vj_(i+1:n)<-k */
#ifdef USE_SSE
			sse_t sse_ji = sseSet(Ra_H[j][i]);
			for (k=i+1 ; k<N_s ; k++)
				Ra_V[i][k] -= Ra_H[j][i]*Ra_V[j][k];
			for (k=N_s ; k<N_e ; k+=sseJmp) {
				sse_t *sse_ik = (sse_t *)(Ra_V[i] + k);
				sse_t *sse_jk = (sse_t *)(Ra_V[j] + k);

				*sse_ik = sseSub(*sse_ik, sseMul(sse_ji, *sse_jk));
				//Ra_V[i][k] -= Ra_H[j][i]*Ra_V[j][k];
			}
			for (k=N_e ; k<N_sz ; k++)
				Ra_V[i][k] -= Ra_H[j][i]*Ra_V[j][k];
#else
			for (k=i+1 ; k<N_sz ; k++)
				Ra_V[i][k] -= Ra_H[j][i]*Ra_V[j][k];
#endif
		}
		/* Settle Hi */
#ifdef USE_SSE
		sse_t sse_ii = sseSet(Ra_H[i][i]);
		for (j=i ; j<N_s ; j++) {
			Ra_V[i][j] += Ra_mat[i][j];
			Ra_V[i][j] /= Ra_H[i][i];
		}
		for (j=N_s ; j<N_e ; j+=sseJmp) {
			sse_t *sse_Vij = (sse_t *)(Ra_V[i] + j);
			sse_t *sse_Mij = (sse_t *)(Ra_mat[i] + j);
			*sse_Vij = sseDiv(sseAdd(*sse_Vij, *sse_Mij), sse_ii);
			//Ra_H[i][j] += Ra_mat[i][j];
		}
		for (j=N_e ; j<N_sz ; j++) {
			Ra_V[i][j] += Ra_mat[i][j];
			Ra_V[i][j] /= Ra_H[i][i];
		}
#else
		for (j=i ; j<N_sz ; j++) {
			Ra_V[i][j] += Ra_mat[i][j];
			Ra_V[i][j] /= Ra_H[i][i];
		}
#endif
	}
	//LOG("HV process %s\n", t.getReadable());

	//t.start();
	/* Set V to H */
	for (i=0 ; i<N_sz ; i++)
		for (j=i+1 ; j<N_sz ; j++) {
			if (Ra_V[i][j] != Ra_V[i][j]) {
				sseUnmat(Ra_V, N_sz);
				sseUnmat(Ra_H, N_sz);
				return NULL;
//				return (wsReal **)1;
			}
			Ra_H[j][i] = Ra_V[i][j];
		}
	//LOG("VH process %s\n", t.getReadable());

	//exportMatrix("LU", Ra_H, N_sz, N_sz);

	sseUnmat(Ra_V, N_sz);
	return Ra_H;

// 	for (k = 0 ; k < N_sz; k++) {
// 		memset(q, 0x0, sizeof(wsReal)*N_sz);
// 		for (j=0 ; j<k ; j++) {
// 			wsReal m = Ra_mat[k][j];
// 			for (p=k ; p<N_sz ; p++)
// 				q[p] += m*Ra_mat[j][p];
// 		}
// 		for (j=k ; j<N_sz ; j++)
// 			Ra_mat[k][j] -= q[j];
// // 		for (j = k; j < N_sz; j++) {
// // 			printf("Fetch k,j [%d,%d]\n", k, j);
// // 			for (p = 0 ; p < k; p++) {
// // 				Ra_mat[k][j] -= Ra_mat[k][p] * Ra_mat[p][j];
// // 				printf("\tAccess k,p [%d,%d]\n", k, p);
// // 				printf("\tAccess p,j [%d,%d]\n\n", p, j);
// // 			}
// // 			s[k][j] = ++u;
// //		}
// 		if ( Ra_mat[k][k] == 0.0 ) return -1;
// 		for (i = k+1 ; i < N_sz; i++) {
// 			printf("Fetch i,k [%d,%d]\n", i, k);
// 			for (p = 0 ; p < k; p++) {
// 				Ra_mat[i][k] -= Ra_mat[i][p] * Ra_mat[p][k];
// 				//*(p_row + k) -= *(p_row + p) * *(p_col + k);
// 				printf("\tAccess i,p [%d,%d]\n", i, p);
// 				printf("\tAccess p,k [%d,%d]\n\n", p, k);
// 			}
// 			Ra_mat[i][k] /= Ra_mat[k][k];
// 			//*(p_row + k) /= *(p_k + k);
// 			//s[i][k] = ++u;
// 		}  
// 	}
	//exportMatrix("ss", s, N_sz, N_sz);
	/*return 0;*/
}

wsReal** invSymMat(wsReal **Ra_mat, wsUint N_row, wsUint B_noParallel/*=0*/)
{
	//	wsReal **Ra_ret = sseSymMat(N_row);
	wsReal **Ra_ret = sseMatrix(N_row, N_row);

	if (N_row == 1) {
		if (Ra_mat[0][0] == W0) {
			LOG("Matrix is singular\n");
			sseUnmat(Ra_ret, 1);

			return NULL;
		}
		Ra_ret[0][0] = W1 / Ra_mat[0][0];

		return Ra_ret;
	} else if (N_row == 2) {
		wsReal R_det = Ra_mat[0][0]*Ra_mat[1][1] - SQR(Ra_mat[1][0]);
		if (R_det == W0) {
			LOG("Matrix is singular\n");
			sseUnmat(Ra_ret, 2);

			return NULL;
		}
		Ra_ret[0][0] = Ra_mat[1][1]/R_det;
		Ra_ret[1][0] = Ra_ret[0][1] = -Ra_mat[1][0]/R_det;
		Ra_ret[1][1] = Ra_mat[0][0]/R_det;

		return Ra_ret;
	}
	wsReal **Ra_LU1 = symLUdecomp(Ra_mat, N_row);
	//	exportMatrix("LU", Ra_LU1, N_row, N_row, 0, 20);
	cTimer t;
	t.start();
	if (Ra_LU1 == NULL) {
		LOG("Matrix is not invertible so --ginv is used instead\n");
		OPTION().assign("ginv", (char *)"1", 1);
		OPTION().FORCE_OPT_NUMBER(ginv);

		sseUnmat(Ra_ret, N_row);
		bool ret = SVDinverse(Ra_mat, N_row, &Ra_ret);
		if (ret == false)
			halt("Generalized inverse with SVD failed");
		exportMatrix("SVD.inv.cor", Ra_ret, N_row, N_row);
	} else if (Ra_LU1 == (wsReal**)1) {
		return NULL;
	} else {
		xLUsolve X_lu = { N_row, Ra_ret, Ra_LU1 };

		if (OPT_NUMBER(thread) > 1 && B_noParallel) {
			WORKER().run(LUsolve, LUdist, &X_lu, NULL, sizeof(wsUint)*2);
		} else for (wsUint i=0 ; i<N_row ; i++)
			symLUsolve(Ra_LU1, i, Ra_ret[i], N_row);
	}
	sseUnmat(Ra_LU1, N_row);
	if (t.get() >= REAL_CONST(1000.0))
		pverbose("solve %s\n", t.getReadable());

	return Ra_ret;
}

wsReal** SsymLUdecomp(wsSym Ra_mat, wsUint N_sz)
{
	wsUint i, j, k;

	//         For each row and column, k = 0, ..., n-1,
	//            find the upper triangular matrix elements for row k
	//            and if the matrix is non-singular (nonzero diagonal element).
	//            find the lower triangular matrix elements for column k. 
	//wsReal **s = sseEmptyMatrix(N_sz, N_sz);
	wsReal *q = NULL;
	sseMalloc(q, wsReal, N_sz);
	//int u=0;
	//cTimer t;

	/* Set H and V */
	wsReal **Ra_H = sseEmptyMatrix(N_sz, N_sz);
	wsReal **Ra_V = sseEmptyMatrix(N_sz, N_sz);

	/* Set H0 to be m_00, ... , m_0n */
	/* Set V0 to be m_00/m_00, m_10/m_00, m_20/m_00, ... m_n0/m_00 */
	wsReal D = Ra_mat[0][0];
	for (i=0 ; i<N_sz ; i++) {
		Ra_H[0][i] = Ra_mat[i][0]; /* If symm */
		Ra_V[0][i] = Ra_mat[i][0] / D; /* If symm */
	}

	/* Calculate H1,V1 ... Hn,Vn */
	//t.start();
	for (i=1 ; i<N_sz ; i++) {
		/* Calculate Hi */

		/* For V0,H0 to V(i-1),H(i-1) */
		for (j=0 ; j<i ; j++) {
			/* Calc Vj_i * Hj_(i:n)<-k */
			for (k=i ; k<N_sz ; k++)
				Ra_H[i][k] -= Ra_V[j][i]*Ra_H[j][k];
		}
		/* Settle Hi */
		for (j=i ; j<N_sz ; j++)
			Ra_H[i][j] += Ra_mat[j][i];

		if (Ra_H[i][i] == W0) {
			sseUnmat(Ra_H, N_sz);
			sseUnmat(Ra_V, N_sz);
			LOG("This matrix is singular\n");	
			return NULL;
		}

		/* Calculate Vi */

		/* For V0,H0 to V(i-1),H(i-1) */
		for (j=0 ; j<i ; j++) {
			/* Calc Hj_i * Vj_(i+1:n)<-k */
			for (k=i+1 ; k<N_sz ; k++)
				Ra_V[i][k] -= Ra_H[j][i]*Ra_V[j][k];
		}
		/* Settle Hi */
		for (j=i ; j<N_sz ; j++) {
			Ra_V[i][j] += Ra_mat[j][i];
			Ra_V[i][j] /= Ra_H[i][i];
		}
	}
	//LOG("HV process %s\n", t.getReadable());

	//t.start();
	/* Set V to H */
	for (i=0 ; i<N_sz ; i++)
		for (j=i+1 ; j<N_sz ; j++)
			Ra_H[j][i] = Ra_V[i][j];
	//LOG("VH process %s\n", t.getReadable());

	//exportMatrix("LU", Ra_H, N_sz, N_sz);

	sseUnmat(Ra_V, N_sz);
	return Ra_H;

// 	for (k = 0 ; k < N_sz; k++) {
// 		memset(q, 0x0, sizeof(wsReal)*N_sz);
// 		for (j=0 ; j<k ; j++) {
// 			wsReal m = Ra_mat[k][j];
// 			for (p=k ; p<N_sz ; p++)
// 				q[p] += m*Ra_mat[j][p];
// 		}
// 		for (j=k ; j<N_sz ; j++)
// 			Ra_mat[k][j] -= q[j];
// // 		for (j = k; j < N_sz; j++) {
// // 			printf("Fetch k,j [%d,%d]\n", k, j);
// // 			for (p = 0 ; p < k; p++) {
// // 				Ra_mat[k][j] -= Ra_mat[k][p] * Ra_mat[p][j];
// // 				printf("\tAccess k,p [%d,%d]\n", k, p);
// // 				printf("\tAccess p,j [%d,%d]\n\n", p, j);
// // 			}
// // 			s[k][j] = ++u;
// //		}
// 		if ( Ra_mat[k][k] == 0.0 ) return -1;
// 		for (i = k+1 ; i < N_sz; i++) {
// 			printf("Fetch i,k [%d,%d]\n", i, k);
// 			for (p = 0 ; p < k; p++) {
// 				Ra_mat[i][k] -= Ra_mat[i][p] * Ra_mat[p][k];
// 				//*(p_row + k) -= *(p_row + p) * *(p_col + k);
// 				printf("\tAccess i,p [%d,%d]\n", i, p);
// 				printf("\tAccess p,k [%d,%d]\n\n", p, k);
// 			}
// 			Ra_mat[i][k] /= Ra_mat[k][k];
// 			//*(p_row + k) /= *(p_k + k);
// 			//s[i][k] = ++u;
// 		}  
// 	}
	//exportMatrix("ss", s, N_sz, N_sz);
	/*return 0;*/
}

wsReal** SinvSymMat(wsReal **Ra_mat, wsUint N_row, wsUint B_noParallel/*=0*/)
{
	//	wsReal **Ra_ret = sseSymMat(N_row);
	wsReal **Ra_ret = sseMatrix(N_row, N_row);

	if (N_row == 1) {
		if (Ra_mat[0][0] == W0) {
			LOG("Matrix is singular\n");
			sseUnmat(Ra_ret, 1);

			return NULL;
		}
		Ra_ret[0][0] = W1 / Ra_mat[0][0];

		return Ra_ret;
	} else if (N_row == 2) {
		wsReal R_det = Ra_mat[0][0]*Ra_mat[1][1] - SQR(Ra_mat[1][0]);
		if (R_det == W0) {
			LOG("Matrix is singular\n");
			sseUnmat(Ra_ret, 2);

			return NULL;
		}
		Ra_ret[0][0] = Ra_mat[1][1]/R_det;
		Ra_ret[1][0] = Ra_ret[0][1] = -Ra_mat[1][0]/R_det;
		Ra_ret[1][1] = Ra_mat[0][0]/R_det;

		return Ra_ret;
	}
	wsReal **Ra_LU1 = SsymLUdecomp(Ra_mat, N_row);
	//	exportMatrix("LU", Ra_LU1, N_row, N_row, 0, 20);
	cTimer t;
	t.start();
	if (Ra_LU1 == NULL) {
		LOG("Matrix is not invertible so --ginv is used instead\n");
		OPTION().assign("ginv", (char *)"1", 1);
		OPTION().FORCE_OPT_NUMBER(ginv);

		sseUnmat(Ra_ret, N_row);
		bool ret = SVDinverse(Ra_mat, N_row, &Ra_ret);
		if (ret == false)
			halt("Generalized inverse with SVD failed");
	} else {
		xLUsolve X_lu = { N_row, Ra_ret, Ra_LU1 };

		if (OPT_NUMBER(thread) > 1 && B_noParallel) {
			WORKER().run(LUsolve, LUdist, &X_lu, NULL, sizeof(wsUint)*2);
		} else for (wsUint i=0 ; i<N_row ; i++)
			symLUsolve(Ra_LU1, i, Ra_ret[i], N_row);
	}
	sseUnmat(Ra_LU1, N_row);
	pverbose("solve %s\n", t.getReadable());

	return Ra_ret;
}

wsMat krnkMM(wsMat Rp_m1, wsUint N_r1, wsUint N_c1, wsMat Rp_m2,
			 wsUint N_r2, wsUint N_c2)
{
	/* Build kronecker matrix */
	wsReal **Ra_kMat = NULL;
	wsUint	N_rkMat = N_r1*N_r2;
	wsUint	N_ckMat = N_c1*N_c2;
	wsAlloc(Ra_kMat, wsReal*, N_rkMat);
	for (wsUint i=0 ; i<N_rkMat ; i++) {
		sseMalloc(Ra_kMat[i], wsReal, N_ckMat);
	}

	for (wsUint i=0 ; i<N_r1 ; i++) {
		wsUint	N_offI	= i*N_r2;
		wsReal	*Ra_1I	= Rp_m1[i];

		for (wsUint j=0 ; j<N_c1 ; j++) {
			wsUint	N_offJ	= j*N_c2;
			register wsReal
				R_1ij	= Ra_1I[j];

			for (wsUint k=0 ; k<N_r2 ; k++) {
				for (wsUint l=0 ; l<N_c2 ; l++) {
					Ra_kMat[N_offI + k][N_offJ + l] =
						R_1ij * Rp_m2[k][l];
				}
			}
		}
	}

	return Ra_kMat;
}

wsMat krnkSS(wsMat Rp_m1, wsUint N_s1, wsMat Rp_m2, wsUint N_s2)
{
	/* Build kronecker matrix */
	wsUint	N_sR	= N_s1*N_s2;
	wsSym	Ra_kMat	= sseSymMat(N_sR);

	/* Do kronecker product */
	for (wsUint i=0 ; i<N_s1 ; i++) {
		wsUint	N_offI	= i*N_s2;
		wsReal*	Ra_1I	= Rp_m1[i];

		/* Full part */
		for (wsUint j=0 ; j<i ; j++) {
			wsUint	N_offJ	= j*N_s2;
			wsReal	R_1ij	= Ra_1I[j];

			for (wsUint k=0 ; k<N_s2 ; k++) {
				for (wsUint l=0 ; l<=k ; l++) {
					Ra_kMat[N_offI+k][N_offJ+l] =
						Ra_kMat[N_offI+l][N_offJ+k] = 
						R_1ij * Rp_m2[k][l];
				}
			}
		}

		/* Sym part */
		wsUint	N_offJ	= i*N_s2;
		wsReal	R_1ii	= Ra_1I[i];
		for (wsUint k=0 ; k<N_s2 ; k++) {
			for (wsUint l=0 ; l<=k ; l++) {
				Ra_kMat[N_offI+k][N_offJ+l] =
					R_1ii * Rp_m2[k][l];
			}
		}
	}

	return Ra_kMat;
}

wsMat krnkSI(wsMat Rp_m1, wsUint N_s1, wsUint N_s2)
{
	/* Build kronecker matrix */
	wsUint	N_sR	= N_s1*N_s2;
	wsSym	Ra_kMat	= sseSymMat(N_sR);

	/* Do kronecker product */
	for (wsUint i=0 ; i<N_s1 ; i++) {
		wsUint	N_offI	= i*N_s2;
		wsReal*	Ra_1I	= Rp_m1[i];

		/* Full part */
		for (wsUint j=0 ; j<i ; j++) {
			wsUint	N_offJ	= j*N_s2;
			wsReal	R_1ij	= Ra_1I[j];

			for (wsUint k=0 ; k<N_s2 ; k++)
				Ra_kMat[N_offI+k][N_offJ+k] = R_1ij;
		}

		/* Sym part */
		wsUint	N_offJ	= i*N_s2;
		wsReal	R_1ii	= Ra_1I[i];
		for (wsUint k=0 ; k<N_s2 ; k++)
			Ra_kMat[N_offI+k][N_offJ+k] = R_1ii;
	}

	return Ra_kMat;
}

wsReal** SO_invMatrix(wsReal **Ra_matrix, wsUint N_row, wsReal *Rp_det)
{
	return SO_invMatrix(Ra_matrix, N_row, 0, 0xffffffff, Rp_det, 0);
}

wsReal** SO_invMatrix(wsReal **Ra_matrix, wsUint N_row)
{
	return SO_invMatrix(Ra_matrix, N_row, 0, 0xffffffff, NULL, 0);
}

wsReal** dbgSO_invMatrix(wsReal **Ra_matrix, wsUint N_row)
{
	return SO_invMatrix(Ra_matrix, N_row, 0, 0xffffffff, NULL, 1);
}

wsReal** SO_invMatrix(wsReal **Ra_matrix, wsUint N_row,
					  wsUint N_start, wsUint N_end, wsReal *Rp_det, char B_isDebug)
{
	wsUint	i, j, k;
	wsReal	D, sD, sK;
	wsReal	*Ra_row;
	wsReal	R_det = W0;

	if (N_end == 0xffffffff)
		N_end = N_row-1;

#ifdef SODBG
	if (B_isDebug)
		printf("SO_invMatrix[%x] %d~%d\n", Ra_matrix, N_start, N_end);
#endif

	/* Make duplication of input matrix */
	wsReal **Ra_inv = NULL;
	wsAlloc(Ra_inv, wsReal*, N_row);
	sseMalloc(Ra_row, wsReal, N_row);
	for (i=0 ; i<N_row ; i++) {
		sseMalloc(Ra_inv[i], wsReal, N_row);
		memcpy(Ra_inv[i], Ra_matrix[i], sizeof(wsReal)*N_row);
	}

	/* When K == 1 */
	if (N_start == 0) {
#ifdef SODBG
		if (B_isDebug)
			printf("SO_invMatrix[%x] Start part\n", Ra_matrix);
#endif

		// # i : 1
		D = Ra_inv[0][0];
		//		if (D == 0.0f)
		//			LOG("This matrix is singular\n");
		R_det = D;
#ifdef USE_SSE
		sse_t sse_D = sseSet(D);
#endif
		wsUint N_e = getMed(N_row);

		Ra_inv[0][0] = W1/D; // case 1
		memcpy(Ra_row, Ra_inv[0], sizeof(wsReal)*N_row);

#ifdef USE_SSE
		//		printf("[1 ~ sseJmp) ");
		if (N_row < sseJmp) {
#ifdef SODBG
			if (B_isDebug)
				printf("SO_invMatrix[%x]::0, Non-sse part\n", Ra_matrix);
#endif

			for (j=1 ; j<N_row ; j++)
				Ra_inv[0][j] /= D; // case 2
		} else {
#ifdef SODBG
			if (B_isDebug)
				printf("SO_invMatrix[%x]::0, sse part\n", Ra_matrix);
#endif

			for (j=1 ; j<sseJmp ; j++)
				Ra_inv[0][j] /= D; // case 2
			//		printf("[sseJmp ~ %d)\n", N_row);
			for (j=sseJmp ; j<N_e ; j+=sseJmp) {
#ifdef SODBG
				if (B_isDebug) {
					if (j+(sseJmp-1) >= N_row)
						halt("SO_invMatrix[%x]::0sse, invalid ragne detected", Ra_matrix);
				}
#endif
				sse_t *sse_pInv = (sse_t *)(Ra_inv[0] + j);
				*sse_pInv = sseDiv(*sse_pInv, sse_D);
			}
			for (j=N_e ; j<N_row ; j++)
				Ra_inv[0][j] /= D; // case 2
		}
#else
		for (j=1 ; j<N_row ; j++)
			Ra_inv[0][j] /= D; // case 2
#endif

		// # i : (1,N]
		for (i=1 ; i<N_row ; i++) {
			sD = Ra_inv[i][0];
#ifdef USE_SSE
			sse_t sse_sD = sseSet(sD);
#endif
			Ra_inv[i][0] = -sD / D; // case 3
#ifdef USE_SSE
			//			printf("[1 ~ sseJmp) ");
			for (j=1 ; j<sseJmp ; j++)
				Ra_inv[i][j] -= sD*Ra_row[j]/D; // case 4
			//			printf("[sseJmp ~ %d)\n", N_row);
			for (j=sseJmp ; j<N_e ; j+=sseJmp) {
				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *vT = (sse_t *)(Ra_inv[i] + j);
				*vT = sseSub(*vT, sseDiv(sseMul(sse_sD, *sse_pRow), sse_D));
			}
			for (j=N_e ; j<N_row ; j++)
				Ra_inv[i][j] -= sD*Ra_row[j]/D; // case 4
#else
			for (j=1 ; j<N_row ; j++)
				Ra_inv[i][j] -= sD*Ra_row[j]/D; // case 4

			Ra_inv[i][j] = Ra_inv[i][j] - sD*Ra_row[j]/D; // case 4
#endif
		}
	}

	/* When K == (1,N) */
	wsUint N_kEnd = N_end<(N_row-1) ? N_end : N_row-1;
	for (k=N_start<2?1:N_start ; k<N_kEnd ; k++) {
#ifdef SODBG
		if (B_isDebug)
			printf("SO_invMatrix[%x]::%d, middle part\n", Ra_matrix, k);
#endif

		memcpy(Ra_row, Ra_inv[k], sizeof(wsReal)*N_row);
		D = Ra_inv[k][k];
		if (D == W0) {
			for (i=0 ; i<N_row ; i++)
				sseFree(Ra_inv[i]);
			DEALLOC(Ra_inv);
			LOG("This matrix is singular\n");
			return NULL;
		}
		R_det *= D;
		wsUint N_e = getMed(N_row);
		//		LOG("DO k = %d\n", k+1);
#ifdef USE_SSE
		wsUint adjK_front = getMed(k);
		wsUint adjK_back = getMed(k+sseJmp);
		if (adjK_back > N_row)
			adjK_back = N_row;
		sse_t sse_D = sseSet(D);
#ifdef SODBG
		if (B_isDebug)
			printf("SO_invMatrix[%x]::0~%d, %d~%d~%d, %d~%d, %d~e[%d], middle part\n", Ra_matrix,
			adjK_front, adjK_front, k, adjK_back, adjK_back, N_e, N_e, N_row);
#endif
#endif

		// # i : [1, k)
		for (i=0 ; i<k ; i++) {
			sK = Ra_inv[i][k];
#ifdef USE_SSE
			sse_t sse_sK = sseSet(sK);
#endif

#ifdef USE_SSE
			//			printf("[0 ~ %d) ", adjK);
			for (j=0 ; j<adjK_front ; j+=sseJmp) {
				//				wsReal res[4];
				//				for (MYuint x=0 ; x<4 ; x++)
				//					res[x] = Ra_inv[i][j+x] - sK*Ra_row[j+x]/D;

				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

				*sse_pInv = sseSub(*sse_pInv, sseDiv(sseMul(sse_sK, *sse_pRow), sse_D));
				//				for (MYuint x=0 ; x<4 ; x++)
				//					if ((j+x) < N_row && (Ra_inv[i][j+x]-res[x])>0.0001)
				//						halt("%g", Ra_inv[i][j+x]-res[x]);
			}
			//			printf("[%d ~ %d) ", adjK, k);
			for (j=adjK_front ; j<k ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#else
			for (j=0 ; j<k ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#endif

			//			printf("[%d] ", k);
			Ra_inv[i][k] = -sK / D; // case 3
#ifdef USE_SSE
			//			printf("[%d~%d) ", k+1, adjK);
			for (j=k+1 ; j<adjK_back ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
			//			printf("[%d~%d)\n", adjK, N_row);
			for (j=adjK_back ; j<N_e ; j+=sseJmp) {
				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

				*sse_pInv = sseSub(*sse_pInv, sseDiv(sseMul(sse_sK, *sse_pRow), sse_D));
			}
			for (j=N_e ; j<N_row ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#else
			for (j=k+1 ; j<N_row ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#endif
		}

		// # i == k
#ifdef USE_SSE
		//			printf("[0 ~ %d) ", adjK);
		for (j=0 ; j<adjK_front ; j+=sseJmp) {
			//				wsReal res[4];
			//				for (MYuint x=0 ; x<4 ; x++)
			//					res[x] = Ra_inv[i][j+x] - sK*Ra_row[j+x]/D;

			sse_t *sse_pRow = (sse_t *)(Ra_row + j);
			sse_t *sse_pInv = (sse_t *)(Ra_inv[k] + j);

			*sse_pInv = sseDiv(*sse_pRow, sse_D);
			//				for (MYuint x=0 ; x<4 ; x++)
			//					if ((j+x) < N_row && (Ra_inv[i][j+x]-res[x])>0.0001)
			//						halt("%g", Ra_inv[i][j+x]-res[x]);
		}
		//			printf("[%d ~ %d) ", adjK, k);
		for (j=adjK_front ; j<k ; j++)
			Ra_inv[k][j] = Ra_row[j] / D; // case 2
#else
		for (j=0 ; j<k ; j++)
			Ra_inv[k][j] = Ra_row[j] / D; // case 2
#endif
		Ra_inv[k][k] = 1.0f/D; // case 1
#ifdef USE_SSE
		//			printf("[%d~%d) ", k+1, adjK);
		for (j=k+1 ; j<adjK_back ; j++)
			Ra_inv[i][j] = Ra_row[j] / D; // case 2
		//			printf("[%d~%d)\n", adjK, N_row);
		for (j=adjK_back ; j<N_e ; j+=sseJmp) {
			sse_t *sse_pRow = (sse_t *)(Ra_row + j);
			sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

			*sse_pInv = sseDiv(*sse_pRow, sse_D);
		}
		for (j=N_e ; j<N_row ; j++)
			Ra_inv[i][j] = Ra_row[j] / D; // case 2
#else
		for (j=k+1 ; j<N_row ; j++)
			Ra_inv[k][j] = Ra_row[j] / D; // case 2
#endif

		// # i == (k, N)
		for (i=k+1 ; i<N_row ; i++) {
			sK = Ra_inv[i][k];
#ifdef USE_SSE
			sse_t sse_sK = sseSet(sK);
#endif

#ifdef USE_SSE
			//			printf("[0 ~ %d) ", adjK);
			for (j=0 ; j<adjK_front ; j+=sseJmp) {
				//				wsReal res[4];
				//				for (MYuint x=0 ; x<4 ; x++)
				//					res[x] = Ra_inv[i][j+x] - sK*Ra_row[j+x]/D;

				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

				*sse_pInv = sseSub(*sse_pInv, sseDiv(sseMul(sse_sK, *sse_pRow), sse_D));
				//				for (MYuint x=0 ; x<4 ; x++)
				//					if ((j+x) < N_row && (Ra_inv[i][j+x]-res[x])>0.0001)
				//						halt("%g", Ra_inv[i][j+x]-res[x]);
			}
			//			printf("[%d ~ %d) ", adjK, k);
			for (j=adjK_front ; j<k ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#else
			for (j=0 ; j<k ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#endif
			Ra_inv[i][k] = -sK / D; // case 3
#ifdef USE_SSE
			//			printf("[%d~%d) ", k+1, adjK);
			for (j=k+1 ; j<adjK_back ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
			//			printf("[%d~%d)\n", adjK, N_row);
			for (j=adjK_back ; j<N_e ; j+=sseJmp) {
				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

				*sse_pInv = sseSub(*sse_pInv, sseDiv(sseMul(sse_sK, *sse_pRow), sse_D));
			}
			for (j=N_e ; j<N_row ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#else
			for (j=k+1 ; j<N_row ; j++)
				Ra_inv[i][j] -= sK*Ra_row[j]/D; // case 4
#endif
		}
	}

	/* When K == N */
	if (N_row > 1 && N_end == (N_row-1)) {
		// # i == [1,N)
		wsUint N = N_row-1;
#ifdef SODBG
		if (B_isDebug)
			printf("SO_invMatrix[%x]::e, last part N=%d\n", Ra_matrix, N);
#endif
		//		LOG("DO k = E\n");
		memcpy(Ra_row, Ra_inv[N], sizeof(wsReal)*N_row);
		D = Ra_inv[N][N];
		if (D == W0) {
			for (i=0 ; i<N_row ; i++)
				sseFree(Ra_inv[i]);
			DEALLOC(Ra_inv);
			LOG("This matrix is singular\n");
			return NULL;
		}
		R_det *= D;

#ifdef USE_SSE
		sse_t sse_D = sseSet(D);
#endif
		for (i=0 ; i<N ; i++) {
#ifdef USE_SSE
			sse_t sse_N = sseSet(Ra_inv[i][N]);
			wsUint adjN = getMed(N);
			for (j=0 ; j<adjN ; j+=sseJmp) {
#ifdef SODBG
				if (B_isDebug) {
					if ((j+(sseJmp-1))>=N)
						halt("SO_invMatrix[%x]::e, Invalid adjN [%d] j [%d~%d]\n", Ra_matrix, adjN, j, j+3);
				}
#endif
				sse_t *sse_pRow = (sse_t *)(Ra_row + j);
				sse_t *sse_pInv = (sse_t *)(Ra_inv[i] + j);

				*sse_pInv = sseSub(*sse_pInv, sseDiv(sseMul(sse_N, *sse_pRow), sse_D));
			}
			for (j=adjN ; j<N ; j++)
				Ra_inv[i][j] -= Ra_inv[i][N]*Ra_row[j]/D; // case 4
#else
			for (j=0 ; j<N ; j++)
				Ra_inv[i][j] -= Ra_inv[i][N]*Ra_row[j]/D; // case 4
#endif
			Ra_inv[i][N] = -Ra_inv[i][N] / D; // case 3
		}
		// # i == N
		for (j=0 ; j<N ; j++)
			Ra_inv[N][j] = Ra_row[j] / D; // case 2
		Ra_inv[N][N] = 1/D; // case 1
	}
	sseFree(Ra_row);

	if (Rp_det)
		*Rp_det = R_det;
	return Ra_inv;
}

void sseMset0(wsMat Ra_m, wsUintCst N_row, wsUintCst N_col)
{
	wsUint N_med = 0;
#ifdef USE_SSE
	N_med = getMed(N_col);
#endif
	LOOP (i, N_row) {
		wsVec	Ra_v = Ra_m[i];
#ifdef USE_SSE
		for (wsUint j=0 ; j<N_med ; j+=sseJmp) {
			sse_t *sse_p = (sse_t *)(Ra_v + j);
			*sse_p = sseXor(*sse_p, *sse_p);
		}
#endif
		for (wsUint j=N_med ; j<N_col ; j++)
			Ra_v[j] = REAL_CONST(0.0);
	}
}

wsRealCst sseMmin(wsMat Ra_m, wsUintCst N_row, wsUintCst N_col,
	wsUint* Np_minIdxR/*=NULL*/, wsUint* Np_minIdxC/*=NULL*/)
{
	wsReal R_val = Ra_m[0][0];
	if (Np_minIdxR) *Np_minIdxR = 0;
	if (Np_minIdxC) *Np_minIdxC = 0;
	/* Find min value */
	LOOP(i, N_row) LOOP(j, N_col) {
		if (Ra_m[i][j] < R_val) {
			R_val = Ra_m[i][j];
			if (Np_minIdxR) *Np_minIdxR = i;
			if (Np_minIdxC) *Np_minIdxC = j;
		}
	}
	return R_val;
}

wsRealCst sseMmax(wsMat Ra_m, wsUintCst N_row, wsUintCst N_col,
	wsUint* Np_maxIdxR/*=NULL*/, wsUint* Np_maxIdxC/*=NULL*/)
{
	wsReal R_val = Ra_m[0][0];
	if (Np_maxIdxR) *Np_maxIdxR = 0;
	if (Np_maxIdxC) *Np_maxIdxC = 0;
	/* Find min value */
	LOOP(i, N_row) LOOP(j, N_col) {
		if (Ra_m[i][j] > R_val) {
			R_val = Ra_m[i][j];
			if (Np_maxIdxR) *Np_maxIdxR = i;
			if (Np_maxIdxC) *Np_maxIdxC = j;
		}
	}
	return R_val;
}

void sseSinit(wsMat Ra_m, wsUintCst N_sz, wsRealCst R_val)
{
	wsUint N_med = 0;
	LOOP (i, N_sz) {
		wsVec	Ra_v = Ra_m[i];
#ifdef USE_SSE
		N_med = getMed(i+1);
		sse_t sse_v = sseSet(R_val);
		for (wsUint j=0 ; j<N_med ; j+=sseJmp)
			*((sse_t *)(Ra_v + j)) = sse_v;
#endif
		for (wsUint j=N_med ; j<=i ; j++)
			Ra_v[j] = R_val;
	}
}

wsMat sseMsubset(wsMat Rp_inpMatrix,
	wsUintCst N_rowStart, wsUintCst N_rowCount,
	wsUintCst N_colStart, wsUintCst N_colCount)
{
	wsReal **Ra_ret = NULL;
	wsAlloc(Ra_ret, wsReal*, N_rowCount);
	for (wsUint i=0 ; i<N_rowCount ; i++) {
		sseMalloc(Ra_ret[i], wsReal, N_colCount);
		memcpy(Ra_ret[i], Rp_inpMatrix[i+N_rowStart]+N_colStart,
			sizeof(wsReal)*N_colCount);
	}

	return Ra_ret;
}

wsMat sseMsubset(wsMat Rp_inpMatrix, wsUint N_sz, vInt& Nv_newSeq)
{
	wsUint	N_newSz	= (wsUint)Nv_newSeq.size();
	wsMat	Ra_ret	= sseMatrix(N_newSz, N_newSz);
	wsUint	I		= 0;

	FOREACHDO (vInt_it, Nv_newSeq, i, I++) {
		wsVec Ra_curInp = Rp_inpMatrix[*i];
		wsVec Ra_curRet	= Ra_ret[I];
		wsUint J = 0;
		FOREACHDO (vInt_it, Nv_newSeq, j, J++)
			Ra_curRet[J] = Ra_curInp[*j];
	}		

	return Ra_ret;
}

wsMat sseMsubset(wsMat Rp_inpMatrix, wsUint N_sz,
	wsUint *Na_newSeq, wsUint N_newSz/*=0xffffffff*/)
{
	if (N_newSz == 0xffffffff)
		N_newSz = N_sz;
	wsReal **Ra_ret = sseMatrix(N_newSz, N_newSz);

	for (wsUint i=0 ; i<N_newSz ; i++)
		for (wsUint j=0 ; j<N_newSz ; j++)
			Ra_ret[i][j] = Rp_inpMatrix[Na_newSeq[i]][Na_newSeq[j]];

	return Ra_ret;
}

wsReal** sseShuffleMatrix(wsReal **Ra_mat, wsUint N_msz, wsUint N_sz, char B_noReplace/*=0*/)
{
	/**/wsUint *Na_idxShuf = NULL;
	wsAlloc(Na_idxShuf, wsUint, N_sz);

	/* Create a vector 0 ... N_samp-1 */
	/**/wsUint *Na_indices = NULL;
	wsAlloc(Na_indices, wsUint, N_msz);
	for (wsUint i=0 ; i<N_msz ; i++)
		Na_indices[i] = i;

	/* Randomly assign (allow oversampling) */
	if (B_noReplace) for (wsUint i=0 ; i<N_sz ; i++) {
		wsUint N_pick = wsRand()%(N_msz-i);
		Na_idxShuf[i] = Na_indices[N_pick];
		/* Delete picked elements */
		memmove(Na_indices+N_pick, Na_indices+N_pick+1,
			sizeof(wsUint)*(N_sz-N_pick-1));
	} else for (wsUint i=0 ; i<N_sz ; i++) {
		wsUint N_pick = wsRand()%N_msz;
		Na_idxShuf[i] = Na_indices[N_pick];
	}

	DEALLOC(Na_indices);

	wsReal **Ra_ret = sseMsubset(Ra_mat, N_sz, Na_idxShuf);
	DEALLOC(Na_idxShuf);

	return Ra_ret;
}

wsSym sseSplitMatrix(wsMat Rp_inpMatrix, wsUint N_sz, char *Ba_isInc,
					 wsUint N_tsz, wsMat *Rp_12, wsSym *Rp_22)
{
	wsMat	Ra_r12	= sseMatrix(N_tsz, N_sz-N_tsz);
	wsSym	Ra_r11	= sseSymMat(N_tsz);
	wsSym	Ra_r22	= sseSymMat(N_sz-N_tsz);

	for (wsUint i=0,I=0,L=0 ; i<N_sz ; i++) {
		if (Ba_isInc[i]) {
			wsUint K=0;
			for (wsUint j=0,J=0 ; j<=i ; j++)
				Ba_isInc[j] ? Ra_r11[I][J++] = Rp_inpMatrix[i][j] :
				Ra_r12[I][K++] = Rp_inpMatrix[i][j];
			for (wsUint j=i+1 ; j<N_sz ; j++)
				if (!Ba_isInc[j]) Ra_r12[I][K++] = Rp_inpMatrix[i][j];
			I++;
		} else {
			for (wsUint j=0,J=0 ; j<=i ; j++)
				if (!Ba_isInc[j]) Ra_r22[L][J++] = Rp_inpMatrix[i][j];
			L++;
		}
	}

	if (Rp_12) *Rp_12 = Ra_r12;
	if (Rp_22) *Rp_22 = Ra_r22;
	return Ra_r11;
}

wsSym sseSplitMatrixV2(wsMat Rp_inpMatrix, wsUint N_sz, char *Ba_isInc,
					   wsUint N_szX1, wsMat *Rp_12)
{
	wsUint	N_szX2	= N_sz - N_szX1;
	wsMat	Ra_r12	= sseMatrix(N_szX1, N_szX2);
	wsSym	Ra_r22	= sseSymMat(N_szX2);

	for (wsUint i=0,I=0,L=0 ; i<N_sz ; i++) {
		if (Ba_isInc[i]) {
			wsUint K=0;
			for (wsUint j=0 ; j<=i ; j++)
				if (!Ba_isInc[j])
					Ra_r12[I][K++] = Rp_inpMatrix[i][j];
			for (wsUint j=i+1 ; j<N_sz ; j++)
				if (!Ba_isInc[j]) Ra_r12[I][K++] = Rp_inpMatrix[i][j];
			I++;
		} else {
			for (wsUint j=0,J=0 ; j<=i ; j++)
				if (!Ba_isInc[j]) Ra_r22[L][J++] = Rp_inpMatrix[i][j];
			L++;
		}
	}

	if (Rp_12) *Rp_12 = Ra_r12;
	return Ra_r22;
}

wsSym sseSplitMatrix2(wsMat Rp_inpMatrix, wsUint N_sz, char *Ba_isExc,
					  wsUint N_tsz, wsMat *Rp_12, wsSym *Rp_22)
{
	// 	MAT_t	Ra_r12	= sseMatrix(N_tsz, N_sz-N_tsz);
	// 	SYM_t	Ra_r11	= sseSymMat(N_tsz);
	// 	SYM_t	Ra_r22	= sseSymMat(N_sz-N_tsz);
	wsMat	Ra_r12	= sseMatrix(N_sz-N_tsz, N_tsz);
	wsSym	Ra_r11	= sseSymMat(N_sz-N_tsz);
	wsSym	Ra_r22	= sseSymMat(N_tsz);

	for (wsUint i=0,I=0,L=0 ; i<N_sz ; i++) {
		if (!Ba_isExc[i]) {
			wsUint K=0;
			for (wsUint j=0,J=0 ; j<=i ; j++)
				!Ba_isExc[j] ? Ra_r11[I][J++] = Rp_inpMatrix[i][j] :
				Ra_r12[I][K++] = Rp_inpMatrix[i][j];
			for (wsUint j=i+1 ; j<N_sz ; j++)
				if (Ba_isExc[j]) Ra_r12[I][K++] = Rp_inpMatrix[i][j];
			I++;
		} else {
			for (wsUint j=0,J=0 ; j<=i ; j++)
				if (Ba_isExc[j]) Ra_r22[L][J++] = Rp_inpMatrix[i][j];
			L++;
		}
	}

	if (Rp_12) *Rp_12 = Ra_r12;
	if (Rp_22) *Rp_22 = Ra_r22;
	return Ra_r11;
}

wsMat sseSubremMatrix(wsReal **Ra_m, wsUint N_sz, wsUint N_st,
	wsUint N_end/*=0xffffffff*/)
{
	/* [st,END] */
	if (N_end == 0xffffffff) {
		wsMat Ra_ret = sseMatrix(N_sz-1, N_sz-1);
		for (wsUint i=0 ; i<N_st ; i++) {
			memcpy(Ra_ret[i], Ra_m[i], sizeof(wsReal)*N_st);
			memcpy(Ra_ret[i]+N_st, Ra_m[i]+N_st+1, sizeof(wsReal)*(N_sz-N_st-1));
		}
		for (wsUint i=N_st+1 ; i<N_sz ; i++) {
			memcpy(Ra_ret[i-1], Ra_m[i], sizeof(wsReal)*N_st);
			memcpy(Ra_ret[i-1]+N_st, Ra_m[i]+N_st+1, sizeof(wsReal)*(N_sz-N_st-1));
		}
		return Ra_ret;
	}

	wsMat Ra_ret = sseMatrix(N_sz-(N_end-N_st+1), N_sz-(N_end-N_st+1));
	for (wsUint i=0 ; i<N_st ; i++) {
		memcpy(Ra_ret[i], Ra_m[i], sizeof(wsReal)*N_st);
		memcpy(Ra_ret[i]+N_st, Ra_m[i]+N_st+1, sizeof(wsReal)*(N_sz-N_end-1));
	}
	for (wsUint i=N_end+1 ; i<N_sz ; i++) {
		memcpy(Ra_ret[i-(N_end-N_st+1)], Ra_m[i], sizeof(wsReal)*N_st);
		memcpy(Ra_ret[i-(N_end-N_st+1)]+N_st, Ra_m[i]+N_st+1, sizeof(wsReal)*(N_sz-N_end-1));
	}
	return Ra_ret;
}

wsSym sseSsubset(wsSym Rp_inpMatrix, wsUint N_start, wsUint N_count)
{
	wsSym Ra_ret = sseSymMat(N_count);

	for (wsUint i=0 ; i<N_count ; i++)
		memcpy(Ra_ret[i], Rp_inpMatrix[i+N_start]+N_start, sizeof(wsReal)*(i+1));

	return Ra_ret;
}

wsSym sseSsubset(wsSym Rp_inpMatrix, wsUint N_sz, vInt& Nv_newSeq)
{
	wsUint	N_newSz	= (wsUint)Nv_newSeq.size();
	wsMat	Ra_ret	= sseMatrix(N_newSz, N_newSz);
	wsUint	I		= 0;

	FOREACHDO (vInt_it, Nv_newSeq, i, I++) {
		wsVec Ra_curInp = Rp_inpMatrix[*i];
		wsVec Ra_curRet	= Ra_ret[I];
		for (wsUint j=0 ; j<=I ; j++)
			Ra_curRet[j] = Ra_curInp[Nv_newSeq[j]];
	}		

	return Ra_ret;
}

wsSym sseSrand(wsUint N_sz)
{
	wsSym Ra_ret = sseSymMat(N_sz);

	for (wsUint i=0 ; i<N_sz ; i++)
#ifdef USE_SYM
		for (wsUint j=0 ; j<=i ; j++)
#else
		for (wsUint j=0 ; j<N_sz ; j++)
#endif
			Ra_ret[i][j] = (wsReal)(wsRand()%100000) / REAL_CONST(1000.0);

	return Ra_ret;
}

void sseMsdmeanC(wsMat Ra_m, wsUintCst N_r, wsUintCst N_c, wsVec* Rp_mean, wsVec* Rp_sd)
{
	/* Compute mean */
	wsVec	Ra_cSum		= sseMsumC(Ra_m, N_r, N_c);
	sseVpC(Ra_cSum, W1 / N_r, Ra_cSum, N_c);

	/* Compute sd */
	wsVec	Ra_cSdsum	= sseEmptyVec(N_c);
	wsVec	Ra_cTmp		= sseVector(N_c);
	LOOP(i, N_r) {
		sseVsVsq(Ra_m[i], Ra_cSum, Ra_cTmp, N_c);
		sseVaV(Ra_cSdsum, Ra_cTmp, Ra_cSdsum, N_c);
	}
	sseVpC(Ra_cSdsum, W1 / (N_r - W1), Ra_cSdsum, N_c);
	sseVsqrt(Ra_cSdsum, N_c);

	*Rp_mean	= Ra_cSum;
	*Rp_sd		= Ra_cSdsum;
	sseFree(Ra_cTmp);
}

} // End namespace ONETOOL
