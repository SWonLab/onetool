#include "utils/permute.h"
#include "global/io.h"
#include "global/option.h"
#include "input/stream.h"
#include "utils/util.h"
#include "utils/matrix.h"

namespace ONETOOL {

cPermuter::cPermuter()
{
	/* No permutation generated */
	N_cnt	= 0;
	/* Size unknown yet */
	N_sz	= 0;
#ifdef _DEBUG
	Cp_perm = cExporter::summon("perm.res");
#endif // _DEBUG
	pthread_mutex_init(H_permSync);
}

cPermuter::~cPermuter()
{
	pthread_mutex_uninit(H_permSync);
#ifdef _DEBUG
	delete Cp_perm;
#endif
}

void cPermuter::setSize(wsUintCst N_inpSz, wsUintCst N_inpCov/*=0*/)
{
	N_sz	= N_inpSz;
	N_cov	= N_inpCov;
}

wsUintCst cPermuter::getSize()
{
	return N_sz;
}

wsUintCst cPermuter::getSizeCov()
{
	return N_cov;
}

bool cPermuter::fetch(cVector& V_permY, cStdMatrix& M_permCov)
{
#ifdef _DEBUG
	/* Size check */
	if (V_permY.size() != N_sz)
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "Permuted Y storage",
		N_sz, V_permY.size());
#endif
	/* Return NULL if ended */
	pthread_mutex_lock(&H_permSync);
	if (isEnd()) {
		pthread_mutex_unlock(&H_permSync);
		return false;
	}

	wsVec Ra_bufP = V_permY.get();
	wsMat Ra_bufC = M_permCov.get();

	/* Generate permutation */
	get(Ra_bufP, Ra_bufC);
	/* Raise count */
	N_cnt++;
	//	LOG("%d\n", N_cnt);
#ifdef _DEBUG
	LOOP(i, N_sz)
		Cp_perm->fmt("%g	", Ra_bufP[i]);
	Cp_perm->put("\n");
#endif // _DEBUG
	/* Copy */
	pthread_mutex_unlock(&H_permSync);

	return true;
}

bool cPermuter::fetch(cStdMatrix& M_permYt, cStdMatrix& M_permCov)
{
#ifdef _DEBUG
	/* Size check */
	if (M_permYt.col() != N_sz)
		halt_fmt(WISARD_SYST_INVL_DIM_MSG, "Permuted Y storage",
		N_sz, M_permYt.col());
#endif
	/* Return NULL if ended */
	pthread_mutex_lock(&H_permSync);
	if (isEnd()) {
		pthread_mutex_unlock(&H_permSync);
		return false;
	}

	wsMat Ra_bufP = M_permYt.get();
	wsMat Ra_bufC = M_permCov.get();

	/* Generate permutation */
	get(Ra_bufP, Ra_bufC);
	/* Raise count */
	N_cnt++;
//	LOG("%d\n", N_cnt);
#ifdef _DEBUG
	LOOP (j, N_pheno) {
		Cp_perm->fmt("PHENOTYPE%d	", j);
		LOOP (i, N_sz)
			Cp_perm->fmt("%g	", Ra_bufP[j][i]);
		Cp_perm->put("\n");
	}
#endif // _DEBUG
	/* Copy */
	pthread_mutex_unlock(&H_permSync);

	return true;
}

void cPermuter::fetchAgain(cVector& V_permY, cStdMatrix& M_permCov)
{
	pthread_mutex_lock(&H_permSync);
	getAgain(V_permY.get(), M_permCov.get());
	pthread_mutex_unlock(&H_permSync);
}

void cPermuter::fetchAgain(cStdMatrix& M_permYt, cStdMatrix& M_permCov)
{
	pthread_mutex_lock(&H_permSync);
	getAgain(M_permYt.get(), M_permCov.get());
	pthread_mutex_unlock(&H_permSync);
}

cPermuter* cPermuter::summon(cIO* Cp_IO, wsMat Ra_Yt, wsUintCst N_pheno, wsUintCst N_anaSamp,
	wsMat Ra_cov, wsUintCst N_inpCov, char* Ba_miss/*=NULL*/)
{
	/* Check options */
	wsUint	N_excl		= IS_ASSIGNED(nperm) + IS_ASSIGNED(seqperm) + IS_ASSIGNED(permfile);
	if (N_excl == 0)
		halt("One of --nperm, --seqperm, --permfile is required!");
	else if (N_excl > 1)
		halt("--nperm, --seqperm, --permfile are mutually exclusive!");

	wsUint		N_perm		= OPT_NUMBER(nperm);// Number of permutations
	wsMat		Ra_valPerm	= NULL;				// values from --permfile [#perm*#anasamp]

	/* --permfile */
	if (IS_ASSIGNED(permfile)) {
		cStrFile	C_permfile(OPT_STRING(permfile), "Permutation values file");
		char*		S_ret		= C_permfile.ugets();
		vSampPtr&	Xv_samp		= Cp_IO->getSample();
		wsUint*		Na_permIdx	= NULL;
		if (S_ret == NULL)
			halt("No data is found in permutation values file!");

		/* Now find out all of the samples in the memory in here */
		mDataIdx Xm_exist;
		FOREACH(vSampPtr_it, Xv_samp, i)
			if (!Ba_miss[(*i)->N_idx]) Xm_exist[(*i)->S_IID] = -1;

		wsUint N_curSamp;
		char** Sa_ret = loadStringValuesByWS(S_ret, &N_curSamp);

		LOOP(i, N_curSamp) {
			if (Xm_exist.find(Sa_ret[i]) == Xm_exist.end()) continue;
			Xm_exist[Sa_ret[i]] = i;
		}

		/* Integrity check */
		wsUint N_noExist = 0;
		FOREACH(mDataIdx_it, Xm_exist, i)
			if (i->second == -1) {
			N_noExist++;
			LOGwarn("Sample IID [%s] does not exist in permutation values file!\n",
				i->first.c_str());
			}
		if (N_noExist)
			halt("[%d] samples in the dataset does not exist in permutation values file, can't proceed!",
			N_noExist);

		/* Get matrix */
		vector<wsVec> Rv_valPerm;
		//		Ra_valPerm = sseMatrix(N_perm, N_anaSamp);
		wsAlloc(Na_permIdx, wsUint, N_anaSamp);

		/* Find out its final indices */
		wsUint I = 0;
		FOREACH(vSampPtr_it, Xv_samp, i) {
			if (Ba_miss[(*i)->N_idx]) continue;
			Na_permIdx[I++] = Xm_exist[(*i)->S_IID];
		}

		/* Retrieve all */
		DEALLOC(S_ret);
		for (wsUint i=0 ; (S_ret=C_permfile.ugets()) ; i++) {
			wsUint N_ret = 0;
			wsVec Ra_tmp = loadRealValuesByWS(S_ret, &N_ret);
			if (N_ret != N_anaSamp)
				halt("Line [%d] of permutation values file contains [%d] "
				"samples (expected [%d] samples)", i+2, N_ret, N_anaSamp);
			wsVec Ra_ins = sseVector(N_anaSamp);
			LOOP(j, N_anaSamp)
				Ra_ins[j] = Ra_tmp[Na_permIdx[j]];
			Rv_valPerm.push_back(Ra_ins);
			sseFree(Ra_tmp);
		}
		N_perm = (wsUint)Rv_valPerm.size();

		/* Cvt it into matrix */
		wsAlloc(Ra_valPerm, wsVec, N_perm);
		I = 0;
		FOREACHDO(vector<wsVec>::iterator, Rv_valPerm, i, I++)
			Ra_valPerm[I] = *i;

		/* Notify and change # of permutations */
		LOGnote("External permutation values file [%s] loaded, [%d] permutations\n",
			OPT_STRING(permfile), Rv_valPerm.size());

		LOOP(i, N_curSamp)
			DEALLOC(Sa_ret[i]);
		DEALLOC(Sa_ret);
	}

	/* --seqperm */
	vector<wsUint*> Nv_seqPerm;
	if (IS_ASSIGNED(seqperm)) {
		/* Load sequence information from file */
		cStrFile	C_seq(OPT_STRING(seqperm), "Permutation sequence file");
		char*		S_ret = NULL;
		for (wsUint L=1 ; (S_ret = C_seq.ugets()) ; L++) {
			wsUint N_curSamp = 0;
			if (S_ret[0]) {
				wsUint* Na_ret = loadIntValuesByWS(S_ret, &N_curSamp);

				/* Check # of samples (should match w/ after-filtering) */
				if (N_anaSamp != N_curSamp)
					halt("Number of sample should be [%d], but [%d] in line [%d]",
						N_anaSamp, N_curSamp, L);

				/* Insert the sequence to the buffer */
				Nv_seqPerm.push_back(Na_ret);

				sseFree(Na_ret);
			}
			DEALLOC(S_ret);
		}

		/* Notify and change # of permutations */
		LOGnote("External permutation sequence file [%s] loaded, [%d] permutations\n",
			OPT_STRING(seqperm), Nv_seqPerm.size());
		N_perm = (wsUint)Nv_seqPerm.size();
	}

	/* Generate instance */
	if (IS_ASSIGNED(nperm))
		return new cSimplePermuter(Ra_Yt, N_pheno, N_anaSamp, Ra_cov, N_inpCov, N_perm);
	else if (IS_ASSIGNED(seqperm))
		return new cSeqPermuter(Ra_Yt, N_pheno, N_anaSamp, Ra_cov, N_inpCov, OPT_STRING(seqperm));
	else if (IS_ASSIGNED(permfile))
		return new cFilePermuter(Cp_IO, Ra_cov, N_inpCov, OPT_STRING(permfile), Ba_miss);

	return NULL;
}

cSimplePermuter::cSimplePermuter(wsVec Ra_inpRef, wsUintCst N_inpSz,
	wsMat Ra_inpCov, wsUintCst N_inpCov, wsUintCst N_inpPerm/*=0*/)
{
	N_pheno = 1;
#ifdef _DEBUG
	Cp_seq = cExporter::summon("perm.seq");
#endif
	/* Default parameter setting */
	if (N_inpPerm == 0) N_perm = OPT_NUMBER(nperm);
	else N_perm = N_inpPerm;

	/* Bind the reference */
	wsAlloc(Ra_ref, wsVec, N_pheno);
	Ra_ref[0]	= Ra_inpRef;
	Ra_refCov	= Ra_inpCov;
	setSize(N_inpSz, N_inpCov);
	Ra_tmp		= sseMatrix(N_pheno, N_inpSz);
	Ra_tmpCov	= sseMatrix(N_inpCov, N_inpSz);
}

cSimplePermuter::cSimplePermuter(wsMat Ra_inpRef, wsUintCst N_row, wsUintCst N_inpSz,
	wsMat Ra_inpCov, wsUintCst N_inpCov, wsUintCst N_inpPerm/*=0*/)
{
	N_pheno = N_row;
#ifdef _DEBUG
	Cp_seq = cExporter::summon("perm.seq");
#endif
	/* Default parameter setting */
	if (N_inpPerm == 0) N_perm = OPT_NUMBER(nperm);
	else N_perm = N_inpPerm;

	/* Bind the reference */
	wsAlloc(Ra_ref, wsVec, N_pheno);
	LOOP (i, N_pheno)
		Ra_ref[i]	= Ra_inpRef[i];
	Ra_refCov	= Ra_inpCov;
	setSize(N_inpSz, N_inpCov);
	Ra_tmp		= sseMatrix(N_pheno, N_inpSz);
	Ra_tmpCov	= sseMatrix(N_inpCov, N_inpSz);
}

cSimplePermuter::~cSimplePermuter()
{
	DEALLOC(Ra_ref);
	sseUnmat(Ra_tmp, N_pheno);
	sseUnmat(Ra_tmpCov, getSizeCov());
#ifdef _DEBUG
	delete Cp_seq;
#endif
}

bool cSimplePermuter::isEnd()
{
	return N_perm <= N_cnt;
}


wsUint cSimplePermuter::getPermCnt()
{
	return N_cnt;
}

void cSimplePermuter::get(wsVec Ra_bufP, wsMat Ra_bufC)
{
	if (N_pheno != 1) halt("Have multiple phenotypes but single permutation asked");
	/* Copy ref -> tmp */
	wsUintCst	N_sz		= getSize();
	wsUintCst	N_cov		= getSizeCov();
	wsUint		N_remained	= N_sz;

	/* Permute it */
#ifdef _DEBUG
	/* Generate sequence */
	wsUint*		Na_seq		= NULL;
	wsAlloc(Na_seq, wsUint, N_sz);
	LOOP (i, N_sz) Na_seq[i] = i;

	LOOP (j, N_sz) {
		wsUintCst	N_idxPick	= wsRand()%N_remained;
		wsUintCst	N_pos		= Na_seq[N_idxPick];

		/* Set covar and pheno */
		Ra_bufP[j] = Ra_ref[0][N_pos];
		LOOP (k, N_cov)
			Ra_bufC[k][j] = Ra_refCov[k][N_pos];

		Cp_seq->fmt("%d	", Na_seq[N_idxPick]);
		memmove(Na_seq+N_idxPick, Na_seq+N_idxPick+1,
			sizeof(wsUint)*(N_remained-N_idxPick-1));
		N_remained--;
	}
	Cp_seq->put("\n");
#else
	memcpy(Ra_tmp[0], Ra_ref[0], sizeof(wsReal)*N_remained);
	LOOP (j, N_cov) memcpy(Ra_tmpCov[j], Ra_refCov[j], sizeof(wsReal)*N_remained);
	LOOP (j, N_sz) {
		wsUint N_idxPick = wsRand()%N_remained;

		/* Set covar and pheno */
		Ra_bufP[j] = Ra_tmp[0][N_idxPick];
		memmove(Ra_tmp[0] + N_idxPick, Ra_tmp[0] + N_idxPick + 1,
			sizeof(wsReal)*(N_remained-N_idxPick-1));
		LOOP (k, N_cov) {
			Ra_bufC[k][j] = Ra_tmpCov[k][N_idxPick];
			memmove(Ra_tmpCov[k]+N_idxPick, Ra_tmpCov[k]+N_idxPick+1,
				sizeof(wsReal)*(N_remained-N_idxPick-1));
		}

		N_remained--;
	}
#endif
}

void cSimplePermuter::get(wsMat Ra_bufP, wsMat Ra_bufC)
{
	/* Copy ref -> tmp */
	wsUintCst	N_sz		= getSize();
	wsUintCst	N_cov		= getSizeCov();
	wsUint		N_remained	= N_sz;

	/* Permute it */
#ifdef _DEBUG
	/* Generate sequence */
	wsUint*		Na_seq		= NULL;
	wsAlloc(Na_seq, wsUint, N_sz);
	LOOP (i, N_sz) Na_seq[i] = i;

	LOOP (j, N_sz) {
		wsUintCst	N_idxPick	= wsRand()%N_remained;
		wsUintCst	N_pos		= Na_seq[N_idxPick];

		/* Set covar and pheno */
		LOOP (k, N_pheno)
			Ra_bufP[k][j] = Ra_ref[k][N_pos];
		LOOP (k, N_cov)
			Ra_bufC[k][j] = Ra_refCov[k][N_pos];

		Cp_seq->fmt("%d	", Na_seq[N_idxPick]);
		memmove(Na_seq+N_idxPick, Na_seq+N_idxPick+1,
			sizeof(wsUint)*(N_remained-N_idxPick-1));
		N_remained--;
	}
	Cp_seq->put("\n");
	DEALLOC(Na_seq);
#else
	LOOP (i, N_pheno)	memcpy(Ra_tmp[i], Ra_ref[i], sizeof(wsReal)*N_remained);
	LOOP (j, N_cov)		memcpy(Ra_tmpCov[j], Ra_refCov[j], sizeof(wsReal)*N_remained);
	LOOP (j, N_sz) {
		wsUint N_idxPick = wsRand()%N_remained;

		/* Set covar and pheno */
		LOOP (k, N_pheno) {
			Ra_bufP[k][j] = Ra_tmp[k][N_idxPick];
			memmove(Ra_tmp[k]+N_idxPick, Ra_tmp[k]+N_idxPick+1,
				sizeof(wsReal)*(N_remained-N_idxPick-1));
		}
		LOOP (k, N_cov) {
			Ra_bufC[k][j] = Ra_tmpCov[k][N_idxPick];
			memmove(Ra_tmpCov[k]+N_idxPick, Ra_tmpCov[k]+N_idxPick+1,
				sizeof(wsReal)*(N_remained-N_idxPick-1));
		}

		N_remained--;
	}
#endif
}

/* getAgain() of simplepermuter is exactly identical to normal get() */
void cSimplePermuter::getAgain(wsVec Ra_bufP, wsMat Ra_bufC)
{
	/* Copy ref -> tmp */
	wsUintCst	N_sz		= getSize();
	wsUintCst	N_cov		= getSizeCov();
	wsUint		N_remained	= N_sz;
	memcpy(Ra_tmp[0], Ra_ref[0], sizeof(wsReal)*N_remained);
	LOOP (j, N_cov) memcpy(Ra_tmpCov[j], Ra_refCov[j], sizeof(wsReal)*N_remained);

	/* Permute it */
	LOOP (j, N_sz) {
		wsUint N_idxPick = wsRand()%N_remained;

		/* Set covar and pheno */
		Ra_bufP[j] = Ra_tmp[0][N_idxPick];
		memmove(Ra_tmp+N_idxPick, Ra_tmp+N_idxPick+1,
			sizeof(wsReal)*(N_remained-N_idxPick-1));
		LOOP (k, N_cov) {
			Ra_bufC[k][j] = Ra_tmpCov[k][N_idxPick];
			memmove(Ra_tmpCov[k]+N_idxPick, Ra_tmpCov[k]+N_idxPick+1,
				sizeof(wsReal)*(N_remained-N_idxPick-1));
		}

		N_remained--;
	}
}

/* getAgain() of simplepermuter is exactly identical to normal get() */
void cSimplePermuter::getAgain(wsMat Ra_bufP, wsMat Ra_bufC)
{
	/* Copy ref -> tmp */
	wsUintCst	N_sz		= getSize();
	wsUintCst	N_cov		= getSizeCov();
	wsUint		N_remained	= N_sz;
	LOOP (j, N_pheno) memcpy(Ra_tmp[j], Ra_ref[j], sizeof(wsReal)*N_remained);
	LOOP (j, N_cov) memcpy(Ra_tmpCov[j], Ra_refCov[j], sizeof(wsReal)*N_remained);

	/* Permute it */
	LOOP(j, N_sz) {
		wsUint N_idxPick = wsRand() % N_remained;

		/* Set covar and pheno */
		LOOP(k, N_pheno) {
			Ra_bufP[k][j] = Ra_tmp[k][N_idxPick];
			memmove(Ra_tmp[k] + N_idxPick, Ra_tmp[k] + N_idxPick + 1,
				sizeof(wsReal)*(N_remained - N_idxPick - 1));
		}
		LOOP(k, N_cov) {
			Ra_bufC[k][j] = Ra_tmpCov[k][N_idxPick];
			memmove(Ra_tmpCov[k] + N_idxPick, Ra_tmpCov[k] + N_idxPick + 1,
				sizeof(wsReal)*(N_remained - N_idxPick - 1));
		}

		N_remained--;
	}
}

wsUint cSeqPermuter::N_szBuf = 1024 * 1024;

cSeqPermuter::cSeqPermuter(wsVec Ra_inpRef, wsUintCst N_inpSz,
	wsMat Ra_inpCov, wsUintCst N_inpCov, wsStrCst S_path)
{
	N_pheno = 1;
	/* Parameter check */
	if (S_path == NULL)
		halt("SYSERR: Pre-sequence permutation requires permutation sequence file!");

	/* Bind the reference */
	wsAlloc(Ra_ref, wsVec, N_pheno);
	Ra_ref[0]	= Ra_inpRef;
	Ra_refCov	= Ra_inpCov;
	setSize(N_inpSz, N_inpCov);

	/* Initiate stream & buffer */
	C_src	= new cStrFile(S_path, "Sequence permutation file");
	N_line	= 0;
	wsAlloc(S_buf, char, N_szBuf);
}

cSeqPermuter::cSeqPermuter(wsMat Ra_inpRef, wsUintCst N_row, wsUintCst N_inpSz,
	wsMat Ra_inpCov, wsUintCst N_inpCov, wsStrCst S_path)
{
	N_pheno = N_row;
	/* Parameter check */
	if (S_path == NULL)
		halt("SYSERR: Pre-sequence permutation requires permutation sequence file!");

	/* Bind the reference */
	Ra_ref		= Ra_inpRef;
	Ra_refCov	= Ra_inpCov;
	setSize(N_inpSz, N_inpCov);

	/* Initiate stream & buffer */
	C_src	= new cStrFile(S_path, "Sequence permutation file");
	N_line	= 0;
	wsAlloc(S_buf, char, N_szBuf);
}

bool cSeqPermuter::isEnd()
{
	N_line++;
	return C_src->end() || (C_src->gets(S_buf, N_szBuf) == NULL);
}


wsUint cSeqPermuter::getPermCnt()
{
	return N_cnt;
}

void cSeqPermuter::get(wsVec Ra_bufP, wsMat Ra_bufC)
{
	if (N_pheno != 1) halt("Have multiple phenotypes but single permutation asked");
	wsUintCst	N_sz		= getSize();
	wsUintCst	N_cov		= getSizeCov();
	wsUint		N_curSamp	= 0;

	if (S_buf[0]) {
		wsUint* Na_ret = loadIntValuesByWS(S_buf, &N_curSamp);

		/* Check # of samples (should match w/ after-filtering) */
		if (N_sz != N_curSamp)
			halt("Number of sample should be [%d], but [%d] in line [%d]",
				N_sz, N_curSamp, N_line);

		/* Fill buffer */
		LOOP (i, N_sz) {
			wsUintCst N_curIdx = Na_ret[i];
			/* Boundary check */
			if (N_curIdx >= N_sz)
				halt("[%d]th sample index [%d] in line [%d] exceeds # of samples [%d]",
					i+1, N_curIdx, N_line, N_sz);

			/* Set cov and pheno */
			Ra_bufP[i] = Ra_ref[0][N_curIdx];
			LOOP (j, N_cov)
				Ra_bufC[j][i] = Ra_refCov[j][N_curIdx];
		}

		sseFree(Na_ret);
	} else
		halt("SYSERR: EOF check of permutation sequence file incomplete");
}

void cSeqPermuter::get(wsMat Ra_bufP, wsMat Ra_bufC)
{
	wsUintCst	N_sz		= getSize();
	wsUintCst	N_cov		= getSizeCov();
	wsUint		N_curSamp	= 0;

	if (S_buf[0]) {
		wsUint* Na_ret = loadIntValuesByWS(S_buf, &N_curSamp);

		/* Check # of samples (should match w/ after-filtering) */
		if (N_sz != N_curSamp)
			halt("Number of sample should be [%d], but [%d] in line [%d]",
			N_sz, N_curSamp, N_line);

		/* Fill buffer */
		LOOP(i, N_sz) {
			wsUintCst N_curIdx = Na_ret[i];
			/* Boundary check */
			if (N_curIdx >= N_sz)
				halt("[%d]th sample index [%d] in line [%d] exceeds # of samples [%d]",
				i + 1, N_curIdx, N_line, N_sz);

			/* Set cov and pheno */
			LOOP(j, N_pheno)
				Ra_bufP[j][i] = Ra_ref[j][N_curIdx];
			LOOP(j, N_cov)
				Ra_bufC[j][i] = Ra_refCov[j][N_curIdx];
		}

		sseFree(Na_ret);
	} else
		halt("SYSERR: EOF check of permutation sequence file incomplete");
}

void cSeqPermuter::getAgain(wsVec Ra_bufC, wsMat Ra_bufP)
{
	/* This permuter does not support repermutation, do nothing */
}

void cSeqPermuter::getAgain(wsMat Ra_bufC, wsMat Ra_bufP)
{
	/* This permuter does not support repermutation, do nothing */
}

cSeqPermuter::~cSeqPermuter()
{
	DEALLOC(S_buf);
	delete C_src;
}

wsUint cFilePermuter::N_szBuf = 1024*1024;

cFilePermuter::cFilePermuter(cIO* Cp_IO, wsMat Ra_inpCov,
	wsUintCst N_inpCov, wsStrCst S_path, const char* Ba_miss/*=NULL*/)
{
	/* Parameter check */
	if (S_path == NULL)
		halt("SYSERR: Pre-sequence permutation requires permutation sequence file!");
	if (N_inpCov > 0)
		LOGwarn("File-based permutation does not support covariate permutation!\n");

	/* Initiate stream & buffer */
	C_src	= new cStrFile(S_path, "Permutation values file");
	N_line	= 0;
	wsAlloc(S_buf, char, N_szBuf);

	/* Set N_sz */
	if (Ba_miss == NULL) setSize(Cp_IO->sizeSample(), N_inpCov);
	else {
		wsUint N_sz = 0;
		wsUint N_samp = Cp_IO->sizeSample();
		LOOP (i, N_samp)
			if (!Ba_miss[i]) N_sz++;
		setSize(N_sz, N_inpCov);
	}

	/* Read first line */
	wsStrCst	S_ret		= C_src->gets(S_buf, N_szBuf);
	vSampPtr&	Xv_samp		= Cp_IO->getSample();

	if (S_ret == NULL)
		halt("No data is found in permutation values file!");

	/* Now find out all of the samples in the memory in here */
	mDataIdx Xm_exist;
	FOREACH (vSampPtr_it, Xv_samp, i)
		if (!Ba_miss[(*i)->N_idx]) Xm_exist[(*i)->S_IID] = -1;

	wsUint N_curSamp;
	char** Sa_ret = loadStringValuesByWS(S_ret, &N_curSamp);

	LOOP(i, N_curSamp) {
		if (Xm_exist.find(Sa_ret[i]) == Xm_exist.end()) continue;
		Xm_exist[Sa_ret[i]] = i;
	}

	/* Integrity check */
	wsUint N_noExist = 0;
	FOREACH (mDataIdx_it, Xm_exist, i)
		if (i->second == -1) {
			N_noExist++;
			LOGwarn("Sample IID [%s] does not exist in permutation values file!\n",
				i->first.c_str());
		}
	if (N_noExist)
		halt("[%d] samples in the dataset does not exist in permutation values file, can't proceed!",
		N_noExist);

	/* Find out its final indices */
	wsAlloc(Na_permIdx, wsUint, getSize());
	wsUint I = 0;
	FOREACH (vSampPtr_it, Xv_samp, i) {
		if (Ba_miss[(*i)->N_idx]) continue;
		Na_permIdx[I++] = Xm_exist[(*i)->S_IID];
	}

	LOOP(i, N_curSamp)
		DEALLOC(Sa_ret[i]);
	DEALLOC(Sa_ret);
}

bool cFilePermuter::isEnd()
{
	N_line++;
	return C_src->end() || (C_src->gets(S_buf, N_szBuf) == NULL);
}


wsUint cFilePermuter::getPermCnt()
{
	return N_cnt;
}

void cFilePermuter::get(wsVec Ra_bufP, wsMat Ra_bufC)
{
	wsUint	N_sz	= getSize();
	wsUint	N_ret	= 0;
	wsVec	Ra_tmp	= loadRealValuesByWS(S_buf, &N_ret);
	if (N_ret != N_sz)
		halt("Line [%d] of permutation values file contains [%d] "
			"samples (expected [%d] samples)", N_line+2, N_ret, N_sz);
	LOOP (j, N_sz)
		Ra_bufP[j] = Ra_tmp[Na_permIdx[j]];
//	LOOP (j, N_cov)
//		memcpy(Ra_bufC[j], Ra_refCov[j], sizeof(wsReal)*N_sz);
	sseFree(Ra_tmp);
}

void cFilePermuter::get(wsMat Ra_bufP, wsMat Ra_bufC)
{
	wsUint	N_sz	= getSize();
	wsUint	N_ret	= 0;
	wsVec	Ra_tmp	= loadRealValuesByWS(S_buf, &N_ret);
	if (N_ret != (N_sz*N_pheno))
		halt("Line [%d] of permutation values file contains [%d] "
			"samples (expected [%d] samples)", N_line+2, N_ret, N_sz*N_pheno);
	LOOP (i, N_pheno) LOOP (j, N_sz)
		Ra_bufP[i][j] = Ra_tmp[N_pheno*i + Na_permIdx[j]];
//	LOOP (j, N_cov)
//		memcpy(Ra_bufC[j], Ra_refCov[j], sizeof(wsReal)*N_sz);
	sseFree(Ra_tmp);
}

void cFilePermuter::getAgain(wsVec Ra_bufP, wsMat Ra_bufC)
{
	/* This permuter does not support repermutation, do nothing */
}

void cFilePermuter::getAgain(wsMat Ra_bufP, wsMat Ra_bufC)
{
	/* This permuter does not support repermutation, do nothing */
}

cFilePermuter::~cFilePermuter()
{
	DEALLOC(S_buf);
	delete C_src;
}

} // End namespace ONETOOL
