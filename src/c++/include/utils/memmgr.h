#pragma once
#ifndef __WISARD_MEMMGR_H__
#define __WISARD_MEMMGR_H__

#include <vector>
#include "utils/util.h"
using namespace std;

namespace ONETOOL {

typedef enum _xMemType
{
	MEMTYPE_NORMAL,
	MEMTYPE_SSE,
} xMemType;

typedef struct _xMemAllocInfo
{
	xMemType	X_type;
	string		S_file, S_varName;
	wsUint		N_size, N_line;
	void		*Vp_ptr;
} xMemAllocInfo;

class cMemManager
{
	vector<xMemAllocInfo> Xa_allocInfo;
	void _report(wsStrCst S_ext=NULL);
public:
	~cMemManager();
	void* alloc(const char *S_varName, const char *S_file, wsUint N_line, wsUint N_size, char B_nullify=0);
	void* alloc16(const char *S_varName, const char *S_file, wsUint N_line, wsUint N_size, char B_nullify=0);
	void memfree(const char *S_proc, int N_line, const void *Vp_ptr);
};

#ifdef USE_MEMDBG
wsMat	_sseMatrix(wsStrCst S_file, wsUint N_line, wsUint N_row, wsUint N_col, wsMat Ra_src = NULL);
wsUmat	_sseUmatrix(wsStrCst S_file, wsUint N_line, wsUint N_row, wsUint N_col, wsMat Ra_src = NULL);
DMAT_t	_dblMatrix(wsStrCst S_file, wsUint N_line, wsUint N_row, wsUint N_col, DMAT_t Ra_src = NULL);
wsMat	_sseEmptyMatrix(wsStrCst S_file, wsUint N_line, wsUint N_row, wsUint N_col);
wsMat	_sseSymMat(wsStrCst S_file, wsUint N_line, wsUint N_sz, wsSym Ra_src = NULL);
wsMat	_sseEmptySymMat(wsStrCst S_file, wsUint N_line, wsUint N_sz);
wsVec	_sseVector(wsStrCst S_file, wsUint N_line, wsUint N_sz);
wsVec	_sseVector(wsStrCst S_file, wsUint N_line, wsUint N_sz, wsVecCst Ra_src);
wsVec	_sseEmptyVec(wsStrCst S_file, wsUint N_line, wsUint N_sz);
#	define sseVector(sz)		_sseVector(__FILE__, __LINE__, sz)
#	define sseVectorP(sz, p)	_sseVector(__FILE__, __LINE__, sz, p)
#	define dblVector(sz)		_dblVector(__FILE__, __LINE__, sz)
#	define sseEmptyVec(sz)		_sseEmptyVec(__FILE__, __LINE__, sz)
#	define sseSymMat(sz)		_sseSymMat(__FILE__, __LINE__, sz)
#	define sseMatrix(r, c)		_sseMatrix(__FILE__, __LINE__, r, c)
#	define sseUmatrix(r, c)		_sseUmatrix(__FILE__, __LINE__, r, c)
#	define dblMatrix(r, c)		_dblMatrix(__FILE__, __LINE__, r, c)
#	define sseSymMatP(sz, p)	_sseSymMat(__FILE__, __LINE__, sz, p)
#	define sseMatrixP(r, c, p)	_sseMatrix(__FILE__, __LINE__, r, c, p)
#	define sseEmptyMatrix(r, c)	_sseEmptyMatrix(__FILE__, __LINE__, r, c)
#	define sseEmptySymMat(sz)	_sseEmptySymMat(__FILE__, __LINE__, sz)
#else
wsMat	_sseMatrix(wsUint N_row, wsUint N_col, wsMat Ra_src = NULL);
wsUmat	_sseUmatrix(wsUint N_row, wsUint N_col, wsMat Ra_src = NULL);
DMAT_t	_dblMatrix(wsUint N_row, wsUint N_col, DMAT_t Ra_src = NULL);
wsMat	_sseEmptyMatrix(wsUint N_row, wsUint N_col);
wsMat	_sseSymMat(wsUint N_sz, wsSym Ra_src = NULL);
wsMat	_sseEmptySymMat(wsUint N_sz);
wsVec	_sseVector(wsUint N_sz);
wsVec	_sseVector(wsUint N_sz, wsVecCst Ra_src);
wsVec	_sseEmptyVec(wsUint N_sz);
#	define sseVector(sz)		_sseVector(sz)
#	define sseVectorP(sz, p)	_sseVector(sz, p)
#	define dblVector(sz)		_dblVector(sz)
#	define sseEmptyVec(sz)		_sseEmptyVec(sz)
#	define sseSymMat(sz)		_sseSymMat(sz)
#	define sseMatrix(r, c)		_sseMatrix(r, c)
#	define sseUmatrix(r, c)		_sseUmatrix(r, c)
#	define dblMatrix(r, c)		_dblMatrix(r, c)
#	define sseSymMatP(sz, p)	_sseSymMat(sz, p)
#	define sseMatrixP(r, c, p)	_sseMatrix(r, c, p)
#	define sseEmptyMatrix(r, c)	_sseEmptyMatrix(r, c)
#	define sseEmptySymMat(sz)	_sseEmptySymMat(sz)
#endif
#define sseUnmat(p,s)	{ deallocMatrix((p), (s), (void *)1); p = NULL; }

} // End namespace ONETOOL

#endif
