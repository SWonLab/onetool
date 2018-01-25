#include "global/common.h"
#include "utils/util.h"
#include "utils/memmgr.h"

namespace ONETOOL {

#ifdef USE_MEMDBG
	cMemManager C_memMgr;
	cMemManager& MEM() { return C_memMgr; }
#endif

cMemManager::~cMemManager()
{
	_report();
}

void* cMemManager::alloc(wsStrCst S_varName, wsStrCst S_file, wsUint N_line, wsUint N_size, char B_nullify/*=0*/)
{
	void *Vp_ret = NULL;
	xMemAllocInfo X_mai;

	/* Allocate */
	Vp_ret = B_nullify ? calloc(1, N_size) : malloc(N_size);
	if (Vp_ret == NULL)
		halt_fmt(WISARD_FAIL_MEMALLOC_W_PROC, N_size, S_file);

	lverbose("0x%x [%s] allocated %d bytes to %s\n", Vp_ret, S_file, N_size,
		S_varName);

	X_mai.X_type = MEMTYPE_NORMAL;
	X_mai.S_file = (string)S_file;
	X_mai.S_varName = (string)S_varName;
	X_mai.Vp_ptr = Vp_ret;
	X_mai.N_size = N_size;
	X_mai.N_line = N_line;
	Xa_allocInfo.push_back(X_mai);

	return Vp_ret;
}

void* cMemManager::alloc16(wsStrCst S_varName, wsStrCst S_file, wsUint N_line, wsUint N_size, char B_nullify/*=0*/)
{
	xMemAllocInfo X_mai;

	/* Adjust size of memory will be allocated, set to multiply of 16 */
	N_size = ((N_size+15)>>4)<<4;
#ifdef __FreeBSD__
	void *Vp_ret = NULL;
	posix_memalign(&Vp_ret, 16, N_size);
#else
	void *Vp_ret = _aligned_malloc(N_size, 16);
#endif

	if (Vp_ret == NULL)
		halt_fmt(WISARD_FAIL_MEMALLOC_W_PROC, N_size, S_file);
	if (B_nullify) memset(Vp_ret, 0x00, N_size);

	lverbose("0x%x [%s] allocated %d SSE bytes to %s\n", Vp_ret, S_file, N_size,
		S_varName);

	X_mai.X_type = MEMTYPE_SSE;
	X_mai.S_file = (string)S_file;
	X_mai.S_varName = (string)S_varName;
	X_mai.Vp_ptr = Vp_ret;
	X_mai.N_size = N_size;
	X_mai.N_line = N_line;
	Xa_allocInfo.push_back(X_mai);	

	return Vp_ret;
}

void cMemManager::memfree(wsStrCst S_proc, int N_line, const void *Vp_ptr)
{
	/* Find this pointer */
	for (vector<xMemAllocInfo>::iterator i=Xa_allocInfo.begin() ;
		i<Xa_allocInfo.end() ;) {
		/* It is! */
		if (i->Vp_ptr == Vp_ptr) {
			lverbose("0x%x freed\n", i->Vp_ptr);

			/* Call free function, according to its type */
			switch (i->X_type) {
			case MEMTYPE_NORMAL:
				if (i->Vp_ptr) ::free(i->Vp_ptr); i->Vp_ptr = NULL; break;
			case MEMTYPE_SSE: MFREE_ALIGNED(i->Vp_ptr); break;
			default:
				halt_fmt(WISARD_SYST_INVL_MEMTYPE, i->X_type);
				break;
			}
			Xa_allocInfo.erase(i);

			return;
		} else i++;
	}

	/* Not found the allocated address having given address */
	LOG("[%s::%d] Address `0x%x` is not registered\n", S_proc, N_line, Vp_ptr);
}

void cMemManager::_report(wsStrCst S_ext/*=NULL*/)
{
	if (S_ext == NULL) {
		LOG("---- A LIST OF NON-FREED MEMORIES ---\n");
		/* Report everything */
		for (vector<xMemAllocInfo>::iterator i=Xa_allocInfo.begin() ;
			i<Xa_allocInfo.end() ; i++) {
				LOG("%s(line %d)::%s : %d bytes\n", i->S_file.c_str(),
					i->N_line, i->S_varName.c_str(), i->N_size);
		}
	} else {
/**/	cExporter *Cp_mem = cExporter::summon(S_ext);
		Cp_mem->put("FILE	LINE	VARNAME	SZ\n");
		for (vector<xMemAllocInfo>::iterator i=Xa_allocInfo.begin() ;
			i<Xa_allocInfo.end() ; i++) {
				Cp_mem->fmt("%s	%d	%s	%d\n", i->S_file.c_str(), i->N_line,
					i->S_varName.c_str(), i->N_size);
		}

		delete Cp_mem;
	}
}

void deallocMatrix(wsUmat Np_mat, wsUint N_row, void *Vp_supp)
{
	if (Np_mat == NULL) return;
#ifdef USE_SSE
	if (Vp_supp != NULL) {
		for (wsUint i = 0; i < N_row; i++) {
			sseFree(Np_mat[i]);
		}
	}
	else for (wsUint i = 0; i < N_row; i++) {
		DEALLOC(Np_mat[i]);
	}
#else
	for (wsUint i = 0; i < N_row; i++) {
		DEALLOC(Np_mat[i]);
	}
#endif
	DEALLOC(Np_mat);
}

#ifdef USE_DBL
void deallocMatrix(wsFloat **Rp_mat, wsUint N_row, void *Vp_supp)
{
	if (Rp_mat == NULL) return;
#ifdef USE_SSE
	if (Vp_supp != NULL) {
		for (wsUint i = 0; i < N_row; i++) {
			sseFree(Rp_mat[i]);
		}
	}
	else for (wsUint i = 0; i < N_row; i++) {
		DEALLOC(Rp_mat[i]);
	}
#else
	for (wsUint i = 0; i < N_row; i++) {
		DEALLOC(Rp_mat[i]);
	}
#endif
	DEALLOC(Rp_mat);
}
#endif

void deallocMatrix(wsReal **Rp_mat, wsUint N_row, void *Vp_supp)
{
	if (Rp_mat == NULL) return;
#ifdef USE_SSE
	if (Vp_supp != NULL) {
		for (wsUint i = 0; i < N_row; i++) {
#ifdef _DEBUG
			if (!Rp_mat[i])
				halt("Tried to deallocate already deallocated part");
#endif
				sseFree(Rp_mat[i]);
		}
	}
	else for (wsUint i = 0; i < N_row; i++) {
		DEALLOC(Rp_mat[i]);
	}
#else
	for (wsUint i = 0; i < N_row; i++) {
		DEALLOC(Rp_mat[i]);
	}
#endif
	DEALLOC(Rp_mat);
}

void deallocMatrix(USHORT_t **Np_mat, wsUint N_row, void *Vp_supp)
{
	if (Np_mat == NULL) return;
#ifdef USE_SSE
	if (Vp_supp != NULL) {
		for (wsUint i = 0; i < N_row; i++) {
			sseFree(Np_mat[i]);
		}
	}
	else for (wsUint i = 0; i < N_row; i++) {
		DEALLOC(Np_mat[i]);
	}
#else
	for (wsUint i = 0; i < N_row; i++) {
		DEALLOC(Np_mat[i]);
	}
#endif
	DEALLOC(Np_mat);
}

wsReal** allocateMatrix(wsUint N_row, wsUint N_col)
{
	wsReal **Ra_ret = NULL;
	wsAlloc(Ra_ret, wsReal*, N_row);
	for (wsUint i = 0; i < N_row; i++) {
		wsAlloc(Ra_ret[i], wsReal, N_col);
	}

	return Ra_ret;
}

#ifdef USE_MEMDBG
wsReal* _sseVector(wsStrCst S_file, wsUint N_line, wsUint N_sz)
#else
wsReal* _sseVector(wsUint N_sz)
#endif
{
	wsReal *Ra_ret = NULL;
	sseMalloc(Ra_ret, wsReal, N_sz);
	return Ra_ret;
}

#ifdef USE_MEMDBG
wsReal* _sseVector(wsStrCst S_file, wsUint N_line, wsUint N_sz, wsVecCst Ra_src)
#else
wsReal* _sseVector(wsUint N_sz, wsVecCst Ra_src)
#endif
{
	wsReal *Ra_ret = NULL;
	sseMalloc(Ra_ret, wsReal, N_sz);
	memcpy(Ra_ret, Ra_src, sizeof(wsReal)*N_sz);
	return Ra_ret;
}

#ifdef USE_MEMDBG
double* _dblVector(wsStrCst S_file, wsUint N_line, wsUint N_sz)
#else
double* _dblVector(wsUint N_sz)
#endif
{
	double *Ra_ret = NULL;
	sseMalloc(Ra_ret, double, N_sz);
	return Ra_ret;
}

#ifdef USE_MEMDBG
wsReal* _sseEmptyVec(wsStrCst S_file, wsUint N_line, wsUint N_sz)
#else
wsReal* _sseEmptyVec(wsUint N_sz)
#endif
{
	wsReal *Ra_ret = NULL;
	sseCalloc(Ra_ret, wsReal, N_sz);
	return Ra_ret;
}

#ifdef USE_MEMDBG
wsMat _sseMatrix(wsStrCst S_file, wsUint N_line, wsUint N_row, wsUint N_col, wsMat Ra_src/*=NULL*/)
#else
wsMat _sseMatrix(wsUint N_row, wsUint N_col, wsMat Ra_src/*=NULL*/)
#endif
{
	wsReal **Ra_ret = NULL;
	wsAlloc(Ra_ret, wsReal*, N_row);

	if (Ra_src) {
		for (wsUint i = 0; i < N_row; i++) {
			sseMallocSym(Ra_ret[i], wsReal, N_col);
			memcpy(Ra_ret[i], Ra_src[i], sizeof(wsReal)*N_col);
		}
	}
	else {
		for (wsUint i = 0; i < N_row; i++) {
			sseMallocSym(Ra_ret[i], wsReal, N_col);
		}
	}

	return Ra_ret;
}

#ifdef USE_MEMDBG
wsUmat _sseUmatrix(wsStrCst S_file, wsUint N_line, wsUint N_row, wsUint N_col, wsMat Ra_src/*=NULL*/)
#else
wsUmat _sseUmatrix(wsUint N_row, wsUint N_col, wsMat Ra_src/*=NULL*/)
#endif
{
	wsUmat Ra_ret = NULL;
	wsAlloc(Ra_ret, wsUint*, N_row);

	if (Ra_src) {
		for (wsUint i = 0; i < N_row; i++) {
			sseMallocSym(Ra_ret[i], wsUint, N_col);
			memcpy(Ra_ret[i], Ra_src[i], sizeof(wsUint)*N_col);
		}
	}
	else {
		for (wsUint i = 0; i < N_row; i++) {
			sseMallocSym(Ra_ret[i], wsUint, N_col);
		}
	}

	return Ra_ret;
}

#ifdef USE_MEMDBG
DMAT_t _dblMatrix(wsStrCst S_file, wsUint N_line, wsUint N_row, wsUint N_col, DMAT_t Ra_src/*=NULL*/)
#else
DMAT_t _dblMatrix(wsUint N_row, wsUint N_col, DMAT_t Ra_src/*=NULL*/)
#endif
{
	wsReal **Ra_ret = NULL;
	wsAlloc(Ra_ret, wsReal*, N_row);

	if (Ra_src) {
		for (wsUint i = 0; i < N_row; i++) {
			sseMallocSym(Ra_ret[i], wsReal, N_col);
			memcpy(Ra_ret[i], Ra_src[i], sizeof(wsReal)*N_col);
		}
	}
	else {
		for (wsUint i = 0; i < N_row; i++) {
			sseMallocSym(Ra_ret[i], wsReal, N_col);
		}
	}

	return Ra_ret;
}

#ifdef USE_MEMDBG
wsReal** _sseSymMat(wsStrCst S_file, wsUint N_line, wsUint N_sz, wsReal **Ra_src/*=NULL*/)
#else
wsReal** _sseSymMat(wsUint N_sz, wsReal **Ra_src/*=NULL*/)
#endif
{
	wsUint	i;
	wsReal	**Ra_ret = NULL;
	MULTI_MALLOC_SYM(Ra_ret, wsReal*, N_sz);

	if (Ra_src) {
		for (i = 0; i < N_sz; i++) {
#ifdef USE_SYM
			sseMallocSym(Ra_ret[i], wsReal, i + 1);
			memcpy(Ra_ret[i], Ra_src[i], sizeof(wsReal)*(i + 1));
#else
			sseMalloc(Ra_ret[i], wsReal, N_sz);
			memcpy(Ra_ret[i], Ra_src[i], sizeof(wsReal)*N_sz);
#endif
		}
	}
	else {
		for (i = 0; i < N_sz; i++) {
#ifdef USE_SYM
			sseMallocSym(Ra_ret[i], wsReal, i + 1);
#else
			sseMalloc(Ra_ret[i], wsReal, N_sz);
#endif
		}
	}

	return Ra_ret;
}

} // End namespace ONETOOL
