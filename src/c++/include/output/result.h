#pragma once
#include "global/common.h"

namespace ONETOOL {

class cResult
{
protected:
	float	*Ra_result;
	char	B_dyn;		///< 1 if memory space is dynamically allocated
	wsUint	N_sz;		///< Allocated size of buffer
	wsUint	N_pos;		///< Filled size of buffer
	int		N_szIdx;	///< Size of indices for each result
	int		N_szRes;	///< Size of results for each result
public:
	cResult();
	cResult(wsUint N_szBuf, wsUint N_inpSzIdx, wsUint N_inpSzRes);
	~cResult();

	virtual void	update(float *Ra_upd, wsUint N_sz)=0;
	float*			getResult();
	wsUintCst			getResultSize();
};

class cCondResult : public cResult
{
	char	B_desc;
	wsUint	N_idxCmp;
	wsReal	R_condVal;
public:
	cCondResult(
		wsUint N_szIdx,			///< Size of indices
		wsUint N_szRes,			///< Size of result
		wsReal R_inpCondVal,	///< Value to be compared
		wsUint N_inpIdxCmp=0,	///< Index of result to be compared
		char B_inpDesc=1		///< Will it compared descending order?
		);
	void update(float *Ra_upd, wsUint N_szUpd);
};

class cTopResult : public cResult
{
	char	B_desc;
	wsUint	N_idxCmp;
public:
	cTopResult(
		wsUint N_top,			///< Number of results will stored
		wsUint N_szIdx,			///< Size of indices
		wsUint N_szRes,			///< Size of result
		wsUint N_inpIdxCmp=0,	///< Index of result to be compared
		char B_inpDesc=1		///< Will it compared descending order?
		);
	void update(float *Ra_upd, wsUint N_szUpd);
};

class cAllResult : public cResult
{
public:
	cAllResult(
		wsUint N_sz,
		wsUint N_szIdx,	///< Size of indices
		wsUint N_szRes	///< Size of result
		);
	void update(float *Ra_upd, wsUint N_szUpd);
};

} // End namespace ONETOOL
