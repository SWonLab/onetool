#include "utils/util.h"
#include "output/result.h"

#define RESULT_BUFFER_INCSIZE	2048

namespace ONETOOL {

cResult::cResult()
{
	Ra_result = NULL;
	N_sz = N_pos = N_szIdx = N_szRes = 0;
}

cResult::cResult(wsUint N_szBuf, wsUint N_inpSzIdx, wsUint N_inpSzRes)
{
	N_sz	= N_szBuf;
	N_szIdx	= N_inpSzIdx;
	N_szRes	= N_inpSzRes;

	/* Allocate memory space */
	if (N_szBuf == 0) {
		sseMalloc(Ra_result, float, RESULT_BUFFER_INCSIZE*(N_szRes+N_szIdx)+2);
	} else {
		sseMalloc(Ra_result, float, N_sz*(N_szRes+N_szIdx)+2);
	}
	N_pos	= 0;
}

cResult::~cResult()
{
	sseFree(Ra_result);
}

float* cResult::getResult()
{
	return Ra_result;
}

wsUintCst cResult::getResultSize()
{
	return N_pos;
}

cTopResult::cTopResult(wsUint N_top, wsUint N_szIdx, wsUint N_szRes,
	wsUint N_inpIdxCmp/*=0*/, char B_inpDesc/*=1*/) : cResult(N_top, N_szIdx, N_szRes)
{
	B_desc		= B_inpDesc;
	N_idxCmp	= N_inpIdxCmp;
}

void cTopResult::update(float *Ra_upd, wsUint N_szUpd) {
	float *Rp_upd = Ra_upd, *Rp_res;
	wsUint i, j;
	wsUint N_szEach = N_szIdx+N_szRes;

	/* For all updates, do comparison and insertion */
	if (B_desc == 1) for (i=0 ; i<N_szUpd ; i++,Rp_upd+=N_szEach) {
		for (j=0,Rp_res=Ra_result ; j<N_pos ; j++,Rp_res+=N_szEach) {
			/* If it is bigger than jth */
			if (Rp_upd[N_szIdx+N_idxCmp] > Rp_res[N_szIdx+N_idxCmp]) {
				/* Determine move size */
				wsUint N_szMove = N_pos-j;
				if (N_szMove >= N_sz) {
					N_szMove = N_sz-1;
					/* If there is no room for last entry, cut it off */
					if (N_pos == N_sz)
						N_szMove--;
				} else if (N_sz > N_pos)
					N_pos++;
				/* Move */
				memmove(Rp_res+N_szEach, Rp_res,
					sizeof(float)*N_szEach*N_szMove);
				/* Insert */
				memcpy(Rp_res, Rp_upd, sizeof(float)*N_szEach);
				break;
			}
		}
		/* Should placed to last */
		if (j == N_pos) {
			/* UNLESS it can be stored, otherwise drop */
			if (N_pos < N_sz) {
				memcpy(Ra_result+N_pos*N_szEach, Rp_upd,
					sizeof(float)*N_szEach);
				N_pos++;
			}
		}
	} else for (i=0 ; i<N_szUpd ; i++,Rp_upd+=N_szEach) {
		for (j=0,Rp_res=Ra_result ; j<N_pos ; j++,Rp_res+=N_szEach) {
			/* If it is smaller than jth */
			if (Rp_upd[N_szIdx+N_idxCmp] < Rp_res[N_szIdx+N_idxCmp]) {
				/* Determine move size */
				wsUint N_szMove = N_pos-j;
				if (N_szMove >= N_sz) {
					N_szMove = N_sz-1;
					/* If there is no room for last entry, cut it off */
					if (N_pos == N_sz)
						N_szMove--;
				} else if (N_sz > N_pos)
					N_pos++;
				/* Move */
				memmove(Rp_res+N_szEach, Rp_res,
					sizeof(float)*N_szEach*N_szMove);
				/* Insert */
				memcpy(Rp_res, Rp_upd, sizeof(float)*N_szEach);
				break;
			}
		}
		/* Should placed to last */
		if (j == N_pos) {
			/* UNLESS it can be stored, otherwise drop */
			if (N_pos < N_sz) {
				memcpy(Ra_result+N_pos*N_szEach, Rp_upd,
					sizeof(float)*N_szEach);
				N_pos++;
			}
		}
	}
}

cAllResult::cAllResult(
	wsUint N_sz,
	wsUint N_szIdx,	///< Size of indices
	wsUint N_szRes	///< Size of result
	) : cResult(N_sz, N_szIdx, N_szRes)
{
}

void cAllResult::update(float *Ra_upd, wsUint N_szUpd)
{
	wsUint i;
	wsUint N_szEach = N_szIdx+N_szRes;
	float *Rp_upd = Ra_upd, *Rp_res = Ra_result+N_szEach*N_pos;

	/* For all updates, do comparison and insertion */
	for (i=0 ; i<N_szUpd ; i++,Rp_upd+=N_szEach,Rp_res+=N_szEach) {
		if (N_pos == N_sz)
			halt("SYSERR: AllResult overflow [%d]", N_sz);
		memcpy(Rp_res, Rp_upd, sizeof(float)*N_szEach);
		N_pos++;
	}
}

cCondResult::cCondResult(wsUint N_szIdx, wsUint N_szRes,
	wsReal R_inpCondVal, wsUint N_inpIdxCmp/*=0*/, char B_inpDesc/*=1*/) :
		cResult(0, N_szIdx, N_szRes)
{
	B_desc		= B_inpDesc;
	N_idxCmp	= N_inpIdxCmp;
	R_condVal	= R_inpCondVal;
}

void cCondResult::update(float *Ra_upd, wsUint N_szUpd)
{
	/* FIXME : Implementation is required */
	halt("This option should be forbidden");
}

} // End namespace ONETOOL
