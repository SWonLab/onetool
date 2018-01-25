#pragma once
#ifndef __WISARD_DATA_H__
#define __WISARD_DATA_H__

#include "global/common.h"
#include "util.h"
#include "utils/matrix.h"
//#include "io/stream.h"

namespace ONETOOL {

typedef map<string,cMatrix*>	mMatrix;
typedef mMatrix::iterator		mMatrix_it;

class cStrFile;
class cDataStorage
{
	char* S_dataPath;
public:
	cDataStorage();
	~cDataStorage();
	void		init(wsStrCst S_inpDP);
	cStrFile*	query(wsStrCst S_path);
};

cDataStorage& D();

class cSharedMatrix
{
	mMatrix		Xa_mats;
public:
	void		set(cMatrix &X_mat, wsStrCst name) {
		mMatrix_it Xa_find = Xa_mats.find(name);
		if (Xa_mats.end() != Xa_find)
			LOGwarn("Shared matrix [%s] replaced\n", name);
		Xa_mats[name] = &X_mat;
	}
	bool		find(wsStrCst S_name) {
		mMatrix_it Xa_find = Xa_mats.find(S_name);
		return Xa_find != Xa_mats.end();
	}
	cMatrix&	get(wsStrCst S_name) {
		mMatrix_it Xa_find = Xa_mats.find(S_name);
		if (Xa_mats.end() == Xa_find)
			halt("Cannot find shared matrix [%s]", S_name);
		return *(Xa_find->second);
	}
};
cSharedMatrix&	S();

class cVariantMap
{
	cIO*	Cp_IO;
	xVariant**	Xa_chr[MAX_NCHR+1];
	wsUint*	Va_chr[MAX_NCHR+1];
	wsUint	Na_chr[MAX_NCHR+1];
public:
	cVariantMap();
	cVariantMap(cIO *Cp_inpIO);
	~cVariantMap();
	void		init(cIO *Cp_inpIO);
	xVariant***	getMarkerMap(wsUint **Np_chr);
	wsUint**	getVariantPosMap(wsUint **Np_chr);
};

struct xMatMul
{
	wsReal	**Ra_ret;
	wsUint	N_mat;
	wsReal	***Ra_mat;
	wsUint	*Na_row;
	wsUint	*Na_col;
};

//cMatrix* multThrRowMM(cMatrix *Cp_mat1, cMatrix *Cp_mat2);
wsReal** multThrRowMM(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2,  wsUint N_mat2Row, wsUint N_mat2Col);

/**
 * multMM performs multiplication of two matrices
 *
 * @param     Rp_mat1		A pointer to first matrix
 * @param     N_mat1Row		Row size of first matrix
 * @param     N_mat1Col		Column size of first matrix, this value should equal to N_mat2Row
 * @param     Rp_mat2		A pointer to second matrix
 * @param     N_mat2Row		Row size of second matrix, this value should equal to N_mat1Col
 * @param     N_mat2Col		Column size of second matrix
 * @return    (wsReal**)	A pointer to matrix have its size to (N_mat1Row*N_mat2Col)
 *							Contains result of multiplication
 */
wsReal** multMM(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2,  wsUint N_mat2Row, wsUint N_mat2Col);

/**
 * multMMt performs multiplication of two matrices, but second matrix should be transposed form
 *
 * @param     Rp_mat1		A pointer to first matrix
 * @param     N_mat1Row		Row size of first matrix
 * @param     N_mat1Col		Column size of first matrix, this value should equal to N_mat2Col
 * @param     Rp_mat2		A pointer to second matrix
 * @param     N_mat2Row		Row size of second matrix
 * @param     N_mat2Col		Column size of second matrix, this value should equal to N_mat1Col
 * @return    (wsReal**)	A pointer to matrix have its size to (N_mat1Row*N_mat2Row)
 *							Contains result of multiplication
 */
wsReal** multMMt(wsReal **Rp_mat1, wsUint N_r1, wsUint N_c1,
	wsReal **Rp_mat2,  wsUint N_r2, wsUint N_c2);

wsReal** sseThrRowMM(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2,  wsUint N_mat2Row, wsUint N_mat2Col);
wsReal** sseThrRowSS(wsReal **Rp_mat1, wsUint N_mat1sz,
	wsReal **Rp_mat2,  wsUint N_mat2sz);
wsReal** sseThrRowMMt(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2=NULL,  wsUint N_mat2Row=0xffffffff, wsUint N_mat2Col=0xffffffff);

/**
 * multMMM performs multiplication of three matrices, but more memory-efficient
 *			than two invocation of multMM()
 *
 * @param     Rp_mat1		A pointer to first matrix
 * @param     N_mat1Row		Row size of first matrix
 * @param     N_mat1Col		Column size of first matrix, this value should equal to N_mat2Row
 * @param     Rp_mat2		A pointer to second matrix
 * @param     N_mat2Row		Row size of second matrix, this value should equal to N_mat1Col
 * @param     N_mat2Col		Column size of second matrix, this value should equal to N_mat3Row
 * @param     Rp_mat3		A pointer to third matrix
 * @param     N_mat3Row		Row size of third matrix, this value should equal to N_mat2Col
 * @param     N_mat3Col		Column size of third matrix
 * @return    (wsReal**)	A pointer to matrix have its size to (N_mat1Row*N_mat3Col)
 *							Contains result of multiplication
 */
wsReal** multMMM(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2,  wsUint N_mat2Row, wsUint N_mat2Col,
	wsReal **Rp_mat3,  wsUint N_mat3Row, wsUint N_mat3Col, int B_export=0);

wsReal multVV(wsRealCst *Ra_v1, wsUint N_sz, wsRealCst *Ra_v2=NULL);

void getIdentityError(wsReal **Ra_matrix, wsUint N_sz);
wsReal** _mult3Matrix(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2, wsUint N_mat2Row, wsUint N_mat2Col,
	wsReal **Rp_mat3, wsUint N_mat3Row, wsUint N_mat3Col);
wsReal** mult3MatrixSSE(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2,  wsUint N_mat2Row, wsUint N_mat2Col,
	wsReal **Rp_mat3,  wsUint N_mat3Row, wsUint N_mat3Col, int B_export=0);
wsReal** mult3MatrixSSE2(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2,  wsUint N_mat2Row, wsUint N_mat2Col,
	wsReal **Rp_mat3,  wsUint N_mat3Row, wsUint N_mat3Col);
wsReal** _mult2Matrix(wsReal **Rp_mat1, wsUint N_mat1Row, wsUint N_mat1Col,
	wsReal **Rp_mat2, wsUint N_mat2Row, wsUint N_mat2Col);

} // End namespace ONETOOL

#endif
