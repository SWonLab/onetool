#pragma once
#ifndef __WISARD_PERMUTE_H__
#define __WISARD_PERMUTE_H__

#include <stdio.h>
#include "global/common.h"
#include "utils/util.h"
#include "global/datatype.h"

namespace ONETOOL {

class cStream;
class cIO;
class cVector;
class cStdMatrix;
class cExporter;

class cPermuter {
	pthread_mutex_t	H_permSync;
#ifdef _DEBUG
	cExporter*		Cp_perm;
#endif
	wsUint			N_sz;
protected:
	wsUint			N_pheno;
	wsUint			N_cov;
	wsUint			N_cnt;
	void			setSize(wsUint N_inpSz, wsUintCst N_inpCov=0);
	wsUintCst		getSize();
	wsUintCst		getSizeCov();
	virtual bool	isEnd()=0;
	virtual void	get(wsVec Ra_bufP, wsMat Ra_bufC)=0;
	virtual void	get(wsMat Ra_bufP, wsMat Ra_bufC)=0;
	virtual void	getAgain(wsVec Ra_bufP, wsMat Ra_bufC)=0;
	virtual void	getAgain(wsMat Ra_bufP, wsMat Ra_bufC)=0;
public:
	cPermuter();
	virtual wsUint	getPermCnt() = 0;
	virtual ~cPermuter();
	bool			fetch(cVector& V_permY, cStdMatrix& M_permCov);
	bool			fetch(cStdMatrix& M_permYt, cStdMatrix& M_permCov);
	void			fetchAgain(cVector& V_permY, cStdMatrix& M_permCov);
	void			fetchAgain(cStdMatrix& M_permYt, cStdMatrix& M_permCov);
	static cPermuter* summon(cIO* Cp_IO, wsVec Ra_Y, wsUintCst N_anaSamp,
		wsMat Ra_inpCov, wsUintCst N_cov, char* Ba_miss=NULL);
	static cPermuter* summon(cIO* Cp_IO, wsMat Ra_Y, wsUintCst N_pheno, wsUintCst N_anaSamp,
		wsMat Ra_inpCov, wsUintCst N_cov, char* Ba_miss=NULL);
};

class cSimplePermuter : public cPermuter {
#ifdef _DEBUG
	cExporter*		Cp_seq;
#endif
	wsMat	Ra_ref;
	wsMat	Ra_refCov;
	wsMat	Ra_tmp;		/*S*/
	wsMat	Ra_tmpCov;
	wsUint	N_perm;
protected:
	bool	isEnd();
	void	get(wsVec Ra_bufP, wsMat Ra_bufC);
	void	get(wsMat Ra_bufP, wsMat Ra_bufC);
	void	getAgain(wsVec Ra_bufP, wsMat Ra_bufC);
	void	getAgain(wsMat Ra_bufP, wsMat Ra_bufC);
public:
	cSimplePermuter(wsVec Ra_inpRef, wsUintCst N_inpSz, wsMat Ra_inpCov,
		wsUintCst N_inpCov, wsUintCst N_perm=0);
	cSimplePermuter(wsMat Ra_inpRef, wsUintCst N_row, wsUintCst N_inpSz, wsMat Ra_inpCov,
		wsUintCst N_inpCov, wsUintCst N_perm=0);
	~cSimplePermuter();
	wsUint	getPermCnt();
};

class cSeqPermuter : public cPermuter {
	static wsUint	N_szBuf;
	wsMat			Ra_ref;
	wsMat			Ra_refCov;
	char*			S_buf;
	cStream*		C_src;
	wsUint			N_line;
protected:
	bool			isEnd();
	void			get(wsVec Ra_bufP, wsMat Ra_bufC);
	void			get(wsMat Ra_bufP, wsMat Ra_bufC);
	void			getAgain(wsVec Ra_bufP, wsMat Ra_bufC);
	void			getAgain(wsMat Ra_bufP, wsMat Ra_bufC);
public:
	cSeqPermuter(wsVec Ra_inpRef, wsUintCst N_inpSz, wsMat Ra_inpCov,
		wsUintCst N_inpCov, wsStrCst S_path);
	cSeqPermuter(wsMat Ra_inpRef, wsUintCst N_row, wsUintCst N_inpSz, wsMat Ra_inpCov,
		wsUintCst N_inpCov, wsStrCst S_path);
	~cSeqPermuter();
	wsUint			getPermCnt();
};

class cFilePermuter : public cPermuter {
	static wsUint	N_szBuf;
	wsUint*			Na_permIdx;
	wsMat			Ra_refCov;
	char*			S_buf;
	cStream*		C_src;
	wsUint			N_line;
protected:
	bool			isEnd();
	void			get(wsVec Ra_bufP, wsMat Ra_bufC);
	void			get(wsMat Ra_bufP, wsMat Ra_bufC);
	void			getAgain(wsVec Ra_bufP, wsMat Ra_bufC);
	void			getAgain(wsMat Ra_bufP, wsMat Ra_bufC);
public:
	cFilePermuter(cIO* Cp_IO, wsMat Ra_inpCov, wsUintCst N_inpCov,
		wsStrCst S_path, const char* Ba_miss=NULL);
	~cFilePermuter();
	wsUint			getPermCnt();
};

} // End namespace ONETOOL

#endif
