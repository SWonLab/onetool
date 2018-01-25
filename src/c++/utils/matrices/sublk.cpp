#include "utils/matrices/sublk.h"
#include "utils/vector.h"
#include <stdarg.h>
#include "global/Rconn.h"

namespace ONETOOL {

cSubMatrix::cSubMatrix(vMatrix &V_inp)
{
	for (vMatrix_it i=V_inp.begin() ; i!=V_inp.end() ; i++)
		V_mats.push_back(*i);
}

cSubMatrix::cSubMatrix(wsUint N_mat, ...)
{
	va_list	H_varList;
	va_start(H_varList, N_mat);

	N_row = 0;
	for (wsUint i=0 ; i<N_mat ; i++) {
		cMatrix *Cp_mat = va_arg(H_varList, cMatrix *);
		if (!Cp_mat->sqr())
			halt("SYSERR: Only square matrix can be the part of sub-matrix");
		V_mats.push_back(Cp_mat);
		N_row += Cp_mat->row();
	}
	N_col = N_row;

	va_end(H_varList);
}

cSubMatrix::~cSubMatrix()
{
	rem();
}

vMatrix& cSubMatrix::get()
{
	return V_mats;
}

/*
 *
 * Essential function definition
 *
 */

wsRealCst cSubMatrix::sum()
{
	/* Sum of submatrix is the sum of sub-matrices of this matrix */
	wsReal R_sum = W0;

	/* For all sub-matrices in this matrix */
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++)
		R_sum += (*i)->sum();

	return R_sum;
}

cMatrix& cSubMatrix::inv()
{
	vMatrix V_ret;

	N_row = 0;
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++) {
		V_ret.push_back(&((*i)->inv()));
		N_row += (*i)->row();
	}
	N_col = N_row;

	return *(new cSubMatrix(V_ret));
}

void cSubMatrix::rem()
{
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++)
		delete *i;
	V_mats.clear();
	N_row = N_col = 0;
}

cVector cSubMatrix::sumR(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsReal	*Ra_ret	= sseVector(N_row);

	wsUint j=0;
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++) {
		cVector V_part = (*i)->sumR(NULL, R_mul);
		memcpy(Ra_ret+j, V_part.get(), sizeof(wsReal)*V_part.size());
		j += V_part.size();
	}

	return cVector(Ra_ret, N_row);
}

cVector cSubMatrix::sumC(wsReal *Rp_sumAll/*=NULL*/, wsReal R_mul/*=WISARD_NAN*/)
{
	wsReal	*Ra_ret	= sseVector(N_row);

	wsUint j=0;
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++) {
		cVector V_part = (*i)->sumC(NULL, R_mul);
		memcpy(Ra_ret+j, V_part.get(), sizeof(wsReal)*V_part.size());
		j += V_part.size();
	}

	return cVector(Ra_ret, N_row);
}

void cSubMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, wsUint N_prec)
{
	LOG("colNames will be ignored\n"); /* FIXME */
	file(S_ext, N_prec);
}

void cSubMatrix::file(wsStrCst S_ext, wsUint N_prec)
{
/**/cExporter* Cp_exporter = cExporter::summon(S_ext);

	wsUint j=0;
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++) {
		wsReal **Rp_cur = (*i)->get();
		wsUint N_cur = (*i)->row();
		wsUint q = row() - j - N_cur;

		for (wsUint k=0 ; k<N_cur ; k++) {
			for (wsUint l=0 ; l<j ; l++)
				Cp_exporter->put("0	");
			for (wsUint l=0 ; l<N_cur ; l++)
				Cp_exporter->fmt("%g	", Rp_cur[k][l]);
			for (wsUint l=0 ; l<q ; l++)
				Cp_exporter->put("0	");
		}

		j += N_cur;
	}

	delete Cp_exporter;
}

void cSubMatrix::file(wsStrCst S_ext, vStr *Sa_colNames, vStr *Sa_rowNames, wsUint N_prec/*=0*/)
{
	/*FIXME*/
	halt("This should not be called");
}

/* For a given constant,
 * Result of a subtraction between the constant and sub-block matrix is
 * a sub-block matrix from the subtraction of each submatrices and constant */
cMatrix& cSubMatrix::operator-(wsReal R)
{
	cSubMatrix	*Cp_ret	= new cSubMatrix(get());
	vMatrix&	V_ret	= Cp_ret->get();

	for (vMatrix_it i=V_mats.begin(),j=V_ret.begin() ; i!=V_mats.end() ;
		i++,j++)
		*(*j) = *(*i) - R;

	return *Cp_ret;
}

/* For a given constant,
 * Result of a addition between the constant and sub-block matrix is
 * a sub-block matrix from the addition of each submatrices and constant */
cMatrix& cSubMatrix::operator+(wsReal R)
{
	cSubMatrix	*Cp_ret	= new cSubMatrix(get());
	vMatrix&	V_ret	= Cp_ret->get();

	for (vMatrix_it i=V_mats.begin(),j=V_ret.begin() ; i!=V_mats.end() ;
		i++,j++)
		*(*j) = *(*i) + R;

	return *Cp_ret;
}

/* For a given constant,
 * Result of a multiplication between the constant and sub-block matrix is
 * a sub-block matrix from the multiplication of each submatrices and constant */
cMatrix& cSubMatrix::operator*(wsReal R)
{
	cSubMatrix	*Cp_ret	= new cSubMatrix(get());
	vMatrix&	V_ret	= Cp_ret->get();

	for (vMatrix_it i=V_mats.begin(),j=V_ret.begin() ; i!=V_mats.end() ;
		i++,j++)
		*(*j) = *(*i) * R;

	return *Cp_ret;
}

cVector cSubMatrix::operator*(cVector &V)
{
	if (row() != V.size()) raise(MATERR_DIM_INCORRECT);
	wsReal	*pR	=sseVector(row());
	wsReal	*pV	= V.get();

	FOREACH (vMatrix_it, V_mats, i) {
		wsMat	Ra_cur	= (*i)->get();
		wsMat	Ra_ret	= NULL;

		/* Calc */
		if (sym())
			Ra_ret = multSMt(Ra_cur, (*i)->row(), &pV, 1, (*i)->row());
		else
			Ra_ret = multMpMt(Ra_cur, (*i)->row(), (*i)->col(), &pV, 1,
				(*i)->row());
		/* Store */
		for (wsUint j=0 ; j<(*i)->row() ; j++)
			pR[j] = Ra_ret[j][0];

		deallocMatrix(Ra_ret, (*i)->row(), (void *)1);
		
		pR += (*i)->row();
		pV += (*i)->row();
	}

	return cVector(pR, row());
}

cMatrix& cSubMatrix::addDiag(wsRealCst R)
{
	halt("This function is still not supported!");
	return *this;
}

void cSubMatrix::setR(wsStrCst S_varName)
{
#ifdef USE_R
	wsReal **Ra_mat = sseEmptyMatrix(row(), col());
	wsUint j=0;
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++) {
		wsReal **Rp_cur = (*i)->get();
		wsUint N_cur = (*i)->row();
		wsUint q = row() - j - N_cur;

		for (wsUint k=0 ; k<N_cur ; k++)
			memcpy(Ra_mat[j+k] + j, Rp_cur[k], sizeof(wsReal)*N_cur);

		j += N_cur;
	}
	R().sendMatrix(S_varName, Ra_mat, row(), row());
	deallocMatrix(Ra_mat, row(), (void *)1);
#else
	halt("USE_R is required");
#endif
}

bool cSubMatrix::sym()
{
	/* Sub-block matrix is virtually symmetric, until before ANY of
	 * matrix is NOT symmetric */
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++)
		if (!(*i)->sym()) return false;
	return true;
}

/* Since sub-block matrix is a sequence of submatrices,
 * the trace of sub-block matrix is a sum of trace from submatrices */
wsReal cSubMatrix::tr()
{
	wsReal R_ret = W0;

	/* For each submatrix, calculate trace and sum up them */
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++)
		R_ret += (*i)->tr();

	return R_ret;
}

cVector cSubMatrix::diag()
{
	wsVec Ra_ret = sseVector(row());

	/* For each submatrix, calculate trace and sum up them */
	wsUint x = 0;
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++) {
		cVector V = (*i)->diag();
		memcpy(Ra_ret+x, V.get(), sizeof(wsReal)*V.size());
		x += V.size();
	}

	return cVector(Ra_ret, row());
}

/*
 *
 * Special operations definition
 *
 */

cStdMatrix& cSubMatrix::mt(cStdMatrix &M)
{
	/* FIXME : Implementation is needed */
	return *(new cStdMatrix());
}

cSubMatrix& cSubMatrix::operator=(cSubMatrix &X)
{
	/* Clear current contents */
	for (vMatrix_it i=V_mats.begin() ; i!=V_mats.end() ; i++)
		delete *i;
	V_mats.clear();

	/* Fill again */
	vMatrix &V_src = X.get();
	for (vMatrix_it i=V_src.begin() ; i!=V_src.end() ; i++)
		V_mats.push_back(*i);
	N_row = N_col = X.row();

	return *this;
}

cSubMatrix& cSubMatrix::clone()
{
	vMatrix V_clone;
	wsUint i=0;
	FOREACHDO (vMatrix_it, V_mats, it, i++)
		V_clone.push_back(&((*it)->clone()));
	return *(new cSubMatrix(V_clone));
}

cSubMatrix& cSubMatrix::sqrt()
{
	vMatrix V_clone;
	wsUint i=0;
	FOREACHDO (vMatrix_it, V_mats, it, i++)
		V_clone.push_back(&((*it)->sqrt()));
	return *(new cSubMatrix(V_clone));
}

wsRealCst cSubMatrix::mean()
{
	wsUint N = row();
	return sum() / (wsReal)(N*N);
}

cVector cSubMatrix::eigen(cStdMatrix &Mp_eVec)
{
	/* FIXME : Implement this */
	halt_fmt(WISARD_SYST_NULL_IMPLEMENT, "Eigendecomposition of sub-block matrix");
	return cVector();
}

cVector cSubMatrix::eigen()
{
	/* FIXME : Implement this */
	halt_fmt(WISARD_SYST_NULL_IMPLEMENT, "Eigendecomposition of sub-block matrix");
	return cVector();
}

cSubMatrix& cSubMatrix::subset(char *Ba_YorN, wsUint N_Yval)
{
	vMatrix V_clone;
	wsUint i=0, j=0;
	FOREACHDO (vMatrix_it, V_mats, it, i++) {
		V_clone.push_back(&((*it)->subset(Ba_YorN+j, N_Yval)));
		j += (*it)->row();
	}
	return *(new cSubMatrix(V_clone));
}

} // End namespace ONETOOL
