#include <math.h>
#include "global/io.h"
#include "analyses/pharaoh.h"
#include "analyses/fam.h"
#include "analyses/setmgr.h"
#include "utils/matrix.h"
#include "input/stream.h"
#include "utils/util.h"
#include "utils/permute.h"
#include "utils/stat.h"

#define GESCA_Winit 99

namespace ONETOOL {

static pthread_mutex_t	H_permInsertSync;
static wsUintCst		N_lambda = 5;

//#define PHARAOHvalidate
#define PHARAOHopti

int do_PHARAOH(char B_overlap,
	wsReal R_lambdaG, wsReal R_lambdaP, cStdMatrix& M_Xt, cStdMatrix& M_y,
	vvInt& M_phenoRel, wsUint N_samp, vInt& Nv_dist, wsUint N_mnfs,
	vStr Sv_latent, mvStr& Sv_mnfsInLatent,
	cStdMatrix& M_W0, vvInt& windex,
	wsVec Ra_retA, wsVec Ra_retW, wsReal* Rp_dev=NULL, wsUint N_maxIter=100, wsVec Ra_seA=NULL);
wsUint do_GeSCA(char B_formative, wsUint b, wsRealCst lambda_b, wsRealCst lambda_w,
	cStdMatrix& M_Xt, cStdMatrix& M_W0_t, cStdMatrix& M_B0_t, cStdMatrix& M_C0_t,
	vvInt& windex, vInt& bindex, vvInt& cindex, xGeSCAres* Xp_res,
	wsUintCst N_obs, wsUintCst N_ltt);

/* !!! PHARAOH-specific analysis !!! */
#if ((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE & FUNCTION_PHARAOH)

wsRealCst doCV_PHARAOH(char B_overlap,
	cStdMatrix& M_Xt, cStdMatrix& M_origYt, vvInt& M_phenoRel,
	wsUint N_curSamp, vInt& Nv_dist, vStr& Sv_idx2mnfs,
	vStr& Sv_latent, mvStr& Sv_mnfsInLatent, cStdMatrix& M_W0,
	vvInt& windex,
	wsUint* Na_curIdx, vReal& Rv_lambdas, wsUint* Np_iter)
{
	wsUint	N_pheno		= M_origYt.row();
	wsUint	N_latent	= (wsUint)Sv_latent.size();
	wsMat	Ra_curXt	= sseMatrix(M_Xt.row(), N_curSamp);
	wsMat	Ra_Xt		= M_Xt.get();
	wsMat	Ra_Yt		= M_origYt.get();
	wsMat	Ra_curYt		= sseMatrix(N_pheno, N_curSamp);

	wsMat	Ra_curW0	= NULL;
	wsAlloc(Ra_curW0, wsVec, M_Xt.row());
	vvInt	curWindex(N_latent);
	mvStr	Sv_curMnfsInLatent;
	FOREACH (mvStr_it, Sv_mnfsInLatent, i)
		Sv_curMnfsInLatent[i->first] = vStr();

	wsReal	R_lambdaG	= Rv_lambdas[0];
	wsReal	R_lambdaP	= Rv_lambdas.size() > 1 ? Rv_lambdas[1] : R_lambdaG;

	/* Subset X */
	wsUint L = 0;
	wsUint N_idxP = 0, N_idxG = 0;
	LOOP (l, M_Xt.row()) {
		char	B_ok = false;

		/* For each row, find abnormality */
		LOOP(o, N_curSamp) {
			wsReal R_v = Ra_Xt[l][Na_curIdx[o]];
			if (R_v) B_ok = true;
			Ra_curXt[L][o] = R_v;
		}

		if (1) {
			if (B_ok) {
				/* Reconstruct W0 */
				Ra_curW0[L] = M_W0.get()[l];

				/* Reconstruct windex */
				LOOP (k, N_latent)
					FOREACH (vInt_it, windex[k], kk)
						if (*kk == (int)l) {
							curWindex[k].push_back(L);//*kk);
							Sv_curMnfsInLatent[Sv_latent[k]].push_back(Sv_idx2mnfs[l]);
						}

				/* Reconstruct Nv_genes */
				L++;
			}
		} else L++;
		
		N_idxG++;
		if (N_idxG == (wsUint)(Sv_mnfsInLatent[Sv_latent[N_idxP]].size())) {
			N_idxG = 0;
			N_idxP++;
		}
	}
	if (0) {//N_idxP != N_latent) {
		LOG("Summary of Sv_mnfsInLatent\n");
		LOG("------------------------------\n");
		FOREACH (mvStr_it, Sv_mnfsInLatent, i)
			LOG("%s %d\n", i->first.c_str(), i->second.size());
		LOG("------------------------------\n");
		LOG("Summary of Sv_latent\n");
		LOG("------------------------------\n");
		FOREACH (vStr_it, Sv_latent, i)
			LOG("%s\n", i->c_str());
		LOG("------------------------------\n");
		halt("SYSERR: Reconstruction error, expect [%d] get [%d]", N_latent, N_idxP);
	}
	wsUint N_curGene = L;
	/* Subset Y */
	LOOP (i, N_pheno) LOOP (o, N_curSamp) Ra_curYt[i][o] = Ra_Yt[i][Na_curIdx[o]];
	cStdMatrix	M_curXt(N_curGene, N_curSamp, Ra_curXt, MATDEL_NONE);
	cStdMatrix	M_curW0(N_curGene, M_W0.col(), Ra_curW0);
	M_curW0.setClean();

	/* Normalize X */
	// X = sqrt(ncase) * zscore(X)/sqrt(ncase-1)
	M_curXt.normalize(sqrt(N_curSamp/(wsReal)(N_curSamp-1)));

	cStdMatrix	M_curYt(N_pheno, N_curSamp, Ra_curYt, MATDEL_NONE);
	wsReal	R_dev	= W0;
	int	N_iter	= 0;
	if (1) N_iter = do_PHARAOH(B_overlap, R_lambdaG, R_lambdaP, M_curXt, M_curYt,
		M_phenoRel, N_curSamp, Nv_dist, N_curGene, Sv_latent, Sv_curMnfsInLatent, M_curW0,
		curWindex, NULL, NULL, &R_dev, OPT_NUMBER(promaxiter));
	else N_iter = do_PHARAOH(B_overlap, R_lambdaG, R_lambdaP, M_curXt, M_curYt,
		M_phenoRel, N_curSamp, Nv_dist, N_curGene, Sv_latent, Sv_mnfsInLatent, M_W0,
		windex, NULL, NULL, &R_dev, OPT_NUMBER(promaxiter));
	if (N_iter <= 0) {
		N_iter = -1;
		LOGwarn("Estimation failed with lambda_G [%g] and lambda_P [%g]\n",
			R_lambdaG, R_lambdaP);
		R_dev = WISARD_NAN;
	}
	if (Np_iter) *Np_iter = N_iter;

	sseUnmat(Ra_curXt, M_Xt.row());
	DEALLOC(Ra_curW0);
	sseUnmat(Ra_curYt, N_pheno);
	return R_dev;
}

wsRealCst doCV_GESCA(char B_formative, cStdMatrix& M_Yt,
	cStdMatrix& M_Xt, cStdMatrix& M_W0, cStdMatrix& M_B0_t, cStdMatrix& M_C0_t,
	vvInt& windex, vInt& bindex, vvInt& cindex,
	wsUintCst N_obs, wsUintCst N_ltt, wsUintCst N_curSamp,
	vStr& Sv_idx2mnfs, vStr& Sv_latent, mvStr& Sv_mnfsInLatent,
	wsUint* Na_curIdx, vReal& Rv_lambdas, wsUint* Np_iter)
{
	wsReal	R_lambdaG	= Rv_lambdas[0];
	wsReal	R_lambdaP	= Rv_lambdas.size() > 1 ? Rv_lambdas[1] : R_lambdaG;
	wsUint	N_anaSamp	= M_Xt.col();
/**/char*	Ba_train	= NULL;
	wsCalloc(Ba_train, char, N_anaSamp);
	LOOP(i, N_curSamp) Ba_train[Na_curIdx[i]] = 1;

	wsUint	N_pheno		= M_Yt.row();
	wsUint	N_latent	= (wsUint)Sv_latent.size();
	wsMat	Ra_curXtTR	= sseMatrix(M_Xt.row(), N_curSamp);
	wsMat	Ra_curXtTE	= sseMatrix(M_Xt.row(), N_anaSamp - N_curSamp);
	wsMat	Ra_Xt		= M_Xt.get();
	wsMat	Ra_Yt		= M_Yt.get();
	wsMat	Ra_curYtTR	= sseMatrix(N_pheno, N_curSamp);
	wsMat	Ra_curYtTE	= sseMatrix(N_pheno, N_anaSamp - N_curSamp);

	wsMat	Ra_curW0	= NULL;
	wsAlloc(Ra_curW0, wsVec, M_Xt.row()+N_pheno);
	vvInt	curWindex(N_latent);
	mvStr	Sv_curMnfsInLatent;
	FOREACH(mvStr_it, Sv_mnfsInLatent, i)
		Sv_curMnfsInLatent[i->first] = vStr();

	/* Subset X */
	wsUint L = 0;
	wsUint N_idxP = 0, N_idxG = 0;
	LOOP(l, M_Xt.row()) {
		char	B_ok = false;

		/* For each row, find abnormality */
		wsUint N_idxTR = 0, N_idxTE = 0;
		LOOP(o, N_anaSamp) {
			wsReal R_v = Ra_Xt[l][o];
			if (Ba_train[o]) {
				if (R_v) B_ok = true;
				Ra_curXtTR[L][N_idxTR++] = R_v;
			} else Ra_curXtTE[L][N_idxTE++] = R_v;
		}

		if (1) {
			if (B_ok) {
				/* Reconstruct W0 */
				Ra_curW0[L] = M_W0.get()[l];

				/* Reconstruct windex */
				LOOP(k, N_latent)
					FOREACH(vInt_it, windex[k], kk)
					if (*kk == (int)l) {
						curWindex[k].push_back(L);//*kk);
						Sv_curMnfsInLatent[Sv_latent[k]].push_back(Sv_idx2mnfs[l]);
					}

				/* Reconstruct Nv_genes */
				L++;
			}
		}
		else L++;

		N_idxG++;
		if (N_idxG == (wsUint)(Sv_mnfsInLatent[Sv_latent[N_idxP]].size())) {
			N_idxG = 0;
			N_idxP++;
		}
	}
	if (0) {//N_idxP != N_latent) {
		LOG("Summary of Sv_mnfsInLatent\n");
		LOG("------------------------------\n");
		FOREACH(mvStr_it, Sv_mnfsInLatent, i)
			LOG("%s %d\n", i->first.c_str(), i->second.size());
		LOG("------------------------------\n");
		LOG("Summary of Sv_latent\n");
		LOG("------------------------------\n");
		FOREACH(vStr_it, Sv_latent, i)
			LOG("%s\n", i->c_str());
		LOG("------------------------------\n");
		halt("SYSERR: Reconstruction error, expect [%d] get [%d]", N_latent, N_idxP);
	}
	wsUint N_curGene = L;
	/* Subset Y */
	LOOP(i, N_pheno) {
		wsUint N_idxTR = 0, N_idxTE = 0;
		LOOP(o, N_anaSamp)
			if (Ba_train[o]) Ra_curYtTR[i][N_idxTR++] = Ra_Yt[i][o];
			else Ra_curYtTE[i][N_idxTE++] = Ra_Yt[i][o];
	}
	DEALLOC(Ba_train);
	LOOP(i, N_pheno)
		Ra_curW0[N_curGene+i] = M_W0.get()[N_obs-N_pheno+i];
	cStdMatrix	M_curW0(N_curGene+N_pheno, M_W0.col(), Ra_curW0);
	M_curW0.setClean(MATDEL_ROW);
	cStdMatrix	M_curW0t = M_curW0.transpose();

	// Add y to both latent and observations
	wsUint	N_newObs	= N_curGene+N_pheno;
	wsMat	__XtTR		= NULL;
	wsAlloc(__XtTR, wsVec, N_newObs);
	LOOP(i, N_curGene)	__XtTR[i] = Ra_curXtTR[i];
	LOOP(i, N_pheno)	__XtTR[N_curGene+i] = Ra_curYtTR[i];
	wsMat	__XtTE		= NULL;
	wsAlloc(__XtTE, wsVec, N_newObs);
	LOOP(i, N_curGene)	__XtTE[i] = Ra_curXtTE[i];
	LOOP(i, N_pheno)	__XtTE[N_curGene+i] = Ra_curYtTE[i];

	xGeSCAres       X_res;
	if (1) {
		cStdMatrix	M_curXtTR(N_newObs, N_curSamp, __XtTR);
		M_curXtTR.setClean(MATDEL_ROW);
		cStdMatrix	M_curYtTR(N_pheno, N_curSamp, Ra_curYtTR);
		wsReal		R_dev	= W0;
		int	N_iter	= do_GeSCA(B_formative, 1,
			R_lambdaP, R_lambdaG,
			M_curXtTR,
			M_curW0t,
			M_B0_t,
			M_C0_t,
			curWindex,
			bindex,
			cindex,
			&X_res,
			N_newObs, N_ltt);
		if (N_iter <= 0) {
			N_iter = -1;
			LOGwarn("Estimation failed with lambda_G [%g] and lambda_P [%g]\n",
				R_lambdaG, R_lambdaP);
			R_dev = WISARD_NAN;
		}
		if (Np_iter) *Np_iter = N_iter;

		sseUnmat(Ra_curXtTR, M_Xt.row());
	}

	wsReal R_PE = WISARD_NAN;
	if (1) {
		char		B_path		= M_curW0t.row() > N_ltt;
		wsUintCst	N_tot		= N_newObs + N_ltt;
		wsUintCst	N_dim		= B_path ? N_tot : N_ltt;
		cStdMatrix	M_curXtTE(N_newObs, N_anaSamp-N_curSamp, __XtTE);
		M_curXtTE.setClean(MATDEL_ROW);
		/* Compute prediction error */ {
			wsUintCst	N_testSamp	= M_curXtTE.col();
			wsUintCst	N_origin	= B_path ? N_newObs : 0;

			// Fill W
			wsUint j = 0;
			LOOP(i, (wsUint)curWindex.size())
				FOREACH(vInt_it, curWindex[i], xx)
					M_curW0t[i+N_origin][*xx] = X_res.MatW[0][j++];
			// Fill B
			j = 0;
			wsMat _Bt = M_B0_t.get();
			FOREACH(vInt_it, bindex, i) {
				wsUint x = (wsUint)(*i/N_dim);
				wsUint y = *i%N_dim;
				_Bt[x][y] = X_res.MatB[0][j++];
			}
			// Fill C
			j = 0;
			wsMat _Ct = M_C0_t.get();
			wsUint x = 0;
			FOREACHDO(vvInt_it, cindex, i, x++) {
				FOREACH(vInt_it, *i, Y) {
					wsUint y = *Y;
					_Ct[x][y] = X_res.MatC[0][j++];
				}
			}
			cStdMatrix& Zt = M_curXtTE.normalize(sqrt(N_testSamp/(wsReal)(N_testSamp-1)));
			wsMat Ra_Zt = Zt.get();
			LOOP (i, Zt.row()) {
				char ok = 1;
				LOOP (j, Zt.col()) if (NA(Ra_Zt[i][j])) ok = 0;
				if (!ok) memset(Ra_Zt[i], 0x00, sizeof(wsReal)*Zt.col());
			}

			// Gamma   <- Z_test %*% W
			cStdMatrix M_Gamma = M_curW0t * Zt;

			// V <- cbind(diag(J),W)
			cIdtMatrix	I(N_newObs);
			cStdMatrix	Vt;
			if (M_C0_t.row()) {
				Vt = I.bmerge(M_curW0t);
			}
			else {
				Vt = M_curW0t;
				Vt.setClean(MATDEL_NONE);
			}
			// Psi     <- Z_test %*% V
			cStdMatrix	PsiT	= Vt * Zt;

			// A
			cStdMatrix At;
			wsUint N_dimA = 0;
			if (M_C0_t.row()) {
				At = M_C0_t.bmerge(M_B0_t);
				N_dimA = N_newObs + N_dim;
			}
			else {
				At = M_B0_t;
				At.setClean(MATDEL_NONE);
				N_dimA = N_dim;
			}

			// dif     <-Psi - Gamma %*% A
			cStdMatrix M_AG = At * M_Gamma;
			cStdMatrix dif = PsiT - M_AG;
			R_PE = dif.tr2();
		}

		// Recover B
		wsUint j = 0;
		wsMat _Bt = M_B0_t.get();
		FOREACH(vInt_it, bindex, i) {
			wsUint x = (wsUint)(*i/N_dim);
			wsUint y = *i%N_dim;
			_Bt[x][y] = GESCA_Winit;
		}
		// Fill C
		j = 0;
		wsMat _Ct = M_C0_t.get();
		wsUint x = 0;
		FOREACHDO(vvInt_it, cindex, i, x++) {
			FOREACH(vInt_it, *i, Y) {
				wsUint y = *Y;
				_Ct[x][y] = GESCA_Winit;
			}
		}

		sseUnmat(Ra_curXtTE, M_Xt.row());
		sseUnmat(Ra_curYtTE, N_pheno);
	}
//	DEALLOC(Ra_curW0);
	return R_PE;
}

/* Running only once
 * Divide SNPs ALMOST equally to each thread */
int forLambda(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
//	int			*Np_idx	= (int *)Vp_thrData;
//	xAnaThread	*Cp_ana	= (xAnaThread *)Vp_shareData;
//	cIO			*Cp_IO	= Cp_ana->Cp_IO;
	static wsUint
				N_proc	= 0;
	int			N_ret	= N_proc != 0 ? 0 : 1;
//	int			i=0;
//	wsReal		j=W0;

	switch (N_mode) {
	case DISTFUN_INIT:
// 		i = 0;
// 		j = W0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		if (N_ret == 0)
			return 0;
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

int forAllPerm_equal(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
//	int		*Np_idx	= (int *)Vp_thrData;
//	xAnaPermThread*
//			Xp_dt	= (xAnaPermThread *)Vp_shareData;
	static wsUint
			N_proc	= 0;
//	int		i=0;
//	wsReal	j=W0;
	int		N_ret = N_proc != 0 ? 0 : 1;

	switch (N_mode) {
	case DISTFUN_INIT:
//		i = 0;
//		j = W0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		if (!N_proc)
			N_ret = N_proc = N_thread;
		break;
	case DISTFUN_UNINIT:
//		LOG("[%d/%d] permutations done\n", N_vrt, N_vrt);
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}

	return N_ret;
}

int thr_PHARAOH(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	xAnaPermThread*
				Xp_dt		= (xAnaPermThread *)Vp_shareData;
	xGscaAnaThread*
				Xp_gs		= (xGscaAnaThread *)Xp_dt->Vp_data;
//	cIO*		Cp_IO		= Xp_dt->Cp_IO;
	char		B_overlap	= Xp_gs->B_overlap;
//	int*		Np_data		= (int *)Vp_data;
//	wsUintCst	N_s			= (wsUint)Np_data[0];
//	wsUintCst	N_e			= (wsUint)Np_data[1];
	wsUintCst	N_samp		= Xp_gs->Na_samp[0];
	wsUintCst	N_maxIter	= OPT_NUMBER(promaxiter) * 2;
	vvInt&		M_phenoRel	= Xp_gs->M_phenoRel;
	cStdMatrix&	M_Yt		= Xp_gs->M_origYt;
	wsUintCst	N_anaSamp	= M_Yt.col();
	wsUintCst	N_pheno		= M_Yt.row();
	wsUintCst	N_cov		= Xp_dt->Cp_IO->sizeCovar();
	cPermuter*	Cp_perm		= Xp_gs->Cp_perm;
	cStdMatrix	M_permYt(N_pheno, N_anaSamp);
	cStdMatrix	M_permCov(N_cov, N_anaSamp);
	wsUint		N_sumWindex = 0;
	FOREACH(vvInt_it, Xp_gs->windex, i) N_sumWindex += (wsUint)(i->size());

	for (wsUint j=0 ; ; j++) {
		if (j == 0 && N_idx == 0) {
			int it = do_PHARAOH(B_overlap, Xp_gs->Rv_lambdas[0][0], Xp_gs->Rv_lambdas[0][1],
				*(Xp_gs->Mp_Xt), Xp_gs->M_origYt, M_phenoRel, N_samp, Xp_gs->Nv_dist,
				Xp_gs->N_mnfs,
				Xp_gs->Sv_latent, Xp_gs->Sv_mnfsInLatent,
				*(Xp_gs->Mp_W0), Xp_gs->windex, Xp_gs->Ra_trueA,
				Xp_gs->Ra_trueW, NULL, N_maxIter);
			if (it <= 0)
				halt("True estimation has failed. Maybe the current lambda"
					" [%g,%g] is too low?", Xp_gs->Rv_lambdas[0][0], Xp_gs->Rv_lambdas[0][1]);
			LOG("True estimation over after [%d] iterations\n", it);
		} else {
			int it = 0;

			/* Fetch anyway */
			if (!Cp_perm->fetch(M_permYt, M_permCov)) {
				LOG("[%d] times of permutation has finished\n", j-1);
				break;
			}

			/* FIXME: Should handle M_permCov */
			char B_success = 0;
			for (wsUint q=0 ; q<3 ; q++) {
				cStdMatrix&	M_curXt = Xp_gs->Mp_Xt->clone();

				/* 150901 Permute cov if required */
//				if (OPT_ENABLED(propermcov)) {
//					wsMat Ra_curXt = M_curXt.get();
//					LOOP(l, N_cov)
//						sseVpermSelf(Ra_curXt[N_gene+l], N_anaSamp,
//							Nv_seqPerm[j]);
//				}

				wsVec Ra_curPermA = sseVector((wsUint)Xp_gs->Sv_mnfsInLatent.size()*N_pheno);
				wsVec Ra_curPermW = sseVector(N_sumWindex);
				it = do_PHARAOH(B_overlap, Xp_gs->Rv_lambdas[0][0], Xp_gs->Rv_lambdas[0][1],
					M_curXt, M_permYt, Xp_gs->M_phenoRel, N_samp, Xp_gs->Nv_dist, Xp_gs->N_mnfs,
					Xp_gs->Sv_latent, Xp_gs->Sv_mnfsInLatent, *(Xp_gs->Mp_W0),
					Xp_gs->windex, Ra_curPermA, Ra_curPermW,
					NULL, N_maxIter);

				delete& M_curXt;
				if (it > 0) {
					B_success = 1;
pthread_mutex_lock(&H_permInsertSync);
					Xp_gs->Rv_permA.push_back(Ra_curPermA);
					Xp_gs->Rv_permW.push_back(Ra_curPermW);
pthread_mutex_unlock(&H_permInsertSync);
					break;
				} else {
					LOG("Permutation [%d] failed, retry...\n", j);
					sseFree(Ra_curPermA);
					sseFree(Ra_curPermW);
					/* Fetch again */
					Cp_perm->fetchAgain(M_permYt, M_permCov);

					/* FIXME: Should handle M_permCov */
				}
			}
			if (!B_success)
				halt("Estimation using doubleRidge failed");
			LOG("Permutation [%d] on thread [%d] over after [%d] iterations\n", j, N_idx, it);
		}
	}

	return 0;
}

int thr_gesca(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	xAnaPermThread*
		Xp_dt		= (xAnaPermThread *)Vp_shareData;
	xGscaAnaThread*
		Xp_gs		= (xGscaAnaThread *)Xp_dt->Vp_data;
	//	cIO*		Cp_IO		= Xp_dt->Cp_IO;
	char		B_fmtv		= true; /* FIXME : Temp set */
	cStdMatrix&	M_Yt		= Xp_gs->M_origYt;
	wsUintCst	N_anaSamp	= M_Yt.col();
	wsUintCst	N_pheno		= M_Yt.row();
	wsUintCst	N_obs		= Xp_gs->N_obs;
	wsUintCst	N_cov		= Xp_dt->Cp_IO->sizeCovar();
	wsRealCst	R_lambdaG	= Xp_gs->Rv_lambdas[0][0];
	wsRealCst	R_lambdaP	= Xp_gs->Rv_lambdas[0][1];
	cPermuter*	Cp_perm		= Xp_gs->Cp_perm;
	cStdMatrix	M_permYt(N_pheno, N_anaSamp);
	cStdMatrix	M_permCov(N_cov, N_anaSamp);
	wsUint		N_sumWindex = 0;
	char		B_pairwise	= OPT_ENABLED(gxg) && !OPT_ENABLED(gxgall);
	FOREACH(vvInt_it, Xp_gs->windex, i) N_sumWindex += (wsUint)(i->size());

	for (wsUint j=0 ; ; j++) {
		if (j == 0 && N_idx == 0) {
			wsUint it = do_GeSCA(B_fmtv, 0, R_lambdaP, R_lambdaG,
				*(Xp_gs->Mp_Xt), *(Xp_gs->Mp_W0), *(Xp_gs->Mp_B0t),
				*(Xp_gs->Mp_C0t), Xp_gs->windex, Xp_gs->bindex,
				Xp_gs->cindex, &Xp_gs->X_res, Xp_gs->N_obs, Xp_gs->N_ltt);
			if (it == 0)
				halt("True estimation has failed. Maybe the current lambda"
				" [%g,%g] is too low?", Xp_gs->Rv_lambdas[0][0], Xp_gs->Rv_lambdas[0][1]);
			LOG("True estimation over after [%d] iterations\n", it);
		} else {
			wsUint it = 0;

			/* Fetch anyway */
			if (!Cp_perm->fetch(M_permYt, M_permCov)) {
				LOG("[%d] times of permutation has finished\n", j-1);
				break;
			}

			/* FIXME: Should handle M_permCov */

			for (wsUint q=0 ; q<3 ; q++) {
				cStdMatrix&	M_curXt = Xp_gs->Mp_Xt->clone();
				wsMat __Xt = M_curXt.get();

				/* Replace phenotype */
				LOOP(j, N_pheno)
					memcpy(__Xt[N_obs - 1 - j], M_permYt.get()[N_pheno - 1 - j], sizeof(wsReal)*N_anaSamp);
				/* Replace covariates */
				if (!IS_ASSIGNED(permfile))
					LOOP(j, N_cov) memcpy(__Xt[N_obs - N_pheno - N_cov + j], M_permCov.get()[j], sizeof(wsReal)*N_anaSamp);

				/* 150901 Permute cov if required */
				//				if (OPT_ENABLED(propermcov)) {
				//					wsMat Ra_curXt = M_curXt.get();
				//					LOOP(l, N_cov)
				//						sseVpermSelf(Ra_curXt[N_gene+l], N_anaSamp,
				//							Nv_seqPerm[j]);
				//				}

				wsVec Ra_curPermA = sseVector((wsUint)Xp_gs->Sv_mnfsInLatent.size()+1);
				wsVec Ra_curPermW = sseVector(N_sumWindex);
				it = do_GeSCA(B_fmtv, 1, R_lambdaP, R_lambdaG,
					M_curXt, *(Xp_gs->Mp_W0), *(Xp_gs->Mp_B0t),
					*(Xp_gs->Mp_C0t), Xp_gs->windex, Xp_gs->bindex,
					Xp_gs->cindex, &Xp_gs->X_res, N_obs, Xp_gs->N_ltt);

				delete& M_curXt;
				if (it != 0) {
					pthread_mutex_lock(&H_permInsertSync);
					Xp_gs->Rv_permA.push_back(Ra_curPermA);
					Xp_gs->Rv_permW.push_back(Ra_curPermW);
					pthread_mutex_unlock(&H_permInsertSync);
					break;
				} else {
					sseFree(Ra_curPermA);
					sseFree(Ra_curPermW);
					/* Fetch again */
					Cp_perm->fetchAgain(M_permYt, M_permCov);

					/* FIXME: Should handle M_permCov */
				}
			}
			if (!B_pairwise || OPT_ENABLED(verbose))
				LOG("Permutation [%d] on thread [%d] over after [%d] iterations\n", j, N_idx, it);
			else
				notice("Permutation [%d] on thread [%d] over after [%d] iterations\r", j, N_idx, it);
		}
	}

	return 0;
}

/* Note : Na_indices will be modified */
wsUint* getPermutedSeq(wsUint N_sz, wsUint* Na_indices=NULL)
{
	wsUintPtr	Na_tmpl		= NULL;
	wsUintPtr	Na_ret		= NULL;
	wsUint		N_remained	= N_sz;
	wsAlloc(Na_ret, wsUint, N_sz);
//	map<int,int> Xa_map1;

	/* Set the template */
	if (Na_indices) {
		Na_tmpl = Na_indices;
//		LOOP (i, N_sz)
//			Xa_map1[Na_indices[i]] = 1;
	} else {
		wsAlloc(Na_tmpl, wsUint, N_sz*10);
		/* Create template */
		LOOP (i, N_sz) {
			Na_tmpl[i] = i;
//			Xa_map1[i] = 1;
		}
	}

	/* Permute it */
	LOOP (j, N_sz) {
		wsUint N_idxPick = wsRand()%N_remained;
		Na_ret[j] = Na_tmpl[N_idxPick];
//		Xa_map1[Na_ret[j]]--;
		memmove(Na_tmpl+N_idxPick, Na_tmpl+N_idxPick+1,
			sizeof(wsUint)*(N_remained-N_idxPick-1));
		N_remained--;
	}

//	for (map<int,int>::iterator i=Xa_map1.begin() ; i!=Xa_map1.end() ; i++)
//		if (i->second) halt("SYSERR: This is an error!");

	/* Dealloc template if it was generated from here */
//	if (!Na_indices)
//		DEALLOC(Na_tmpl);

	return Na_ret;
}

void getCactIndex(wsVecCst Ra_origY, wsUint N_samp, wsUint* Np_case, wsUint* Np_ctrl,
	wsUint** Na_idxPermCa, wsUint** Na_idxPermCt, bool B_perm)
{
	wsUint	N_case		= 0;
	wsUint	N_ctrl		= 0;
	LOOP (i, N_samp)
		if (Ra_origY[i] == WISARD_AFFECTED) N_case++;
		else if (Ra_origY[i] == WISARD_UNAFFECTED) N_ctrl++;
		else halt("SYSERR: Missing phenotype found in the array of non-missing samples!");

	/* Make arrays of indices for cases/controls */
	wsUintPtr	Na_idxCases	= NULL;
	wsUintPtr	Na_idxCtrls	= NULL;
	wsAlloc(Na_idxCases, wsUint, N_case);
	wsAlloc(Na_idxCtrls, wsUint, N_ctrl);
	N_case = N_ctrl = 0;
	LOOP (i, N_samp)
		if (Ra_origY[i] == WISARD_AFFECTED) Na_idxCases[N_case++] = i;
		else Na_idxCtrls[N_ctrl++] = i;
	/* Make permuted */
	if (B_perm) {
		*Na_idxPermCa	= getPermutedSeq(N_case, Na_idxCases);
		*Na_idxPermCt	= getPermutedSeq(N_ctrl, Na_idxCtrls);
		DEALLOC(Na_idxCases);
		DEALLOC(Na_idxCtrls);
	} else {
		*Na_idxPermCa	= Na_idxCases;
		*Na_idxPermCt	= Na_idxCtrls;
	}
	if (Np_case) *Np_case = N_case;
	if (Np_ctrl) *Np_ctrl = N_ctrl;
}

//#define PHARAOHvalidate

wsReal getqi(wsVec yi, wsVec mui, wsUint N_sz, wsUint datatype)
{
	wsReal	R_ret	= W0;

	if (datatype == 1) /* normal - identity */
		R_ret = -sseVsVsum(yi, mui, N_sz, 1) / W2;
	else if (datatype == 2) LOOP (i, N_sz) { /* binary - logit */
		wsReal R_1m = W1 - mui[i];
		R_ret += yi[i]*log(mui[i] / R_1m) + log(R_1m);
		//QI = sum(yi*log(mui/(1-mui)) + log(1 - mui))
	} else if (datatype == 3) LOOP (i, N_sz) /* poisson - log */
		R_ret += yi[i]*log(mui[i]) - mui[i];
		//QI = sum(yi*log(mui) - mui)
	else if (datatype == 4) LOOP (i, N_sz) /* gamma - inverse */
		R_ret += -yi[i]/mui[i] - log(mui[i]);
		//QI = sum((-yi/mui) - log(mui))
	return R_ret;
}

/* Calculates mean vector of current GLM 
 * lp = linear predictor */
cVector getmu(cVector& V_lp, wsUint datatype)
{
	wsUint	N_sz	= V_lp.size();
	wsVec	lp		= V_lp.get();
	wsVec	Ra_ret	= sseVector(N_sz);
	if (datatype == 1) /* normal - identity */
		memcpy(Ra_ret, lp, sizeof(wsReal)*N_sz);
	else if (datatype == 2) LOOP (i, N_sz) { /* binary - logit */
		wsReal R_tmp = exp(lp[i]) / (W1 + exp(lp[i]));
		if (NA(R_tmp) || R_tmp > REAL_CONST(1000.0))
			R_tmp = REAL_CONST(1000.0);
		Ra_ret[i] = R_tmp;
	} else if (datatype == 3) LOOP (i, N_sz)/* poisson - log */
		Ra_ret[i] = exp(lp[i]);
	else if (datatype == 4) LOOP (i, N_sz) { /* gamma - inverse  */
		wsReal R_tmp = N_sz / lp[i];
		Ra_ret[i] = R_tmp <= W0 ?  W1 : R_tmp;
// 		mu	= ones(size(lp)) / lp
// 		mu	= mu * (mu>0) + (mu<=0) #makes sure that mui>0
	}
	return cVector(Ra_ret, N_sz);
}

cVector getvmu(cVector& V_mu, wsUint dist)
{
	wsUintCst	N_sz	= V_mu.size();
	wsVecCst	mu		= V_mu.get();
	wsVec	Ra_ret	= sseVector(N_sz);
	if (dist == 1) /*Normal */
		sseVinit(Ra_ret, N_sz, W1);
		//vmu = ones(size(mu));
	else if (dist == 2) LOOP (i, N_sz) { /* Binomial/m -1 */
		wsReal R_mu = mu[i]<REAL_CONST(0.0001) ? mu[i]+0.0001 : mu[i];
		if (R_mu == W1) R_mu -= 0.0001;
		// mu = mu + 0.0001*(mu<0.0001);
		Ra_ret[i] = R_mu * (W1 - R_mu);
		// vmu = mu.*(1-mu);
	} else if (dist == 3) LOOP (i, N_sz) /* Poisson */
		Ra_ret[i] = mu[i]==W0 ? 0.00001 : mu[i];
		// vmu = mu + 0.00001*(mu==0);
	else if (dist == 4) { /* Gamma */
		sseVsquare(mu, N_sz, Ra_ret);
		LOOP (i, N_sz) if (Ra_ret[i] == W0) Ra_ret[i] = 0.00000001;
		// vmu = (mu.^2)+0.00000001*(mu==0);      %in case  mu=0
	} else if (dist == 5) /* Binomial/m -2 */
		sseVp1p(mu, N_sz, Ra_ret, 1);
		// vmu = (mu.^2).*(1-mu).^2;
	else if (dist == 6) { /* Inverse Gaussian */
		sseVcube(mu, N_sz, Ra_ret);
		LOOP (i, N_sz) if (Ra_ret[i] == W0) Ra_ret[i] = 0.00000001;
		//vmu = (mu.^3)+0.00000001*(mu==0);
	} else if (dist == 7) { /* Negative Binomial */
		halt("SYSERR: k is not yet defined");
		wsReal k = W1;
		sseVsquare(mu, N_sz, Ra_ret, W1/k);
		sseVaV(Ra_ret, mu, Ra_ret, N_sz);
		// vmu = mu + (mu.^2)/k;
	}
	return cVector(Ra_ret, N_sz);
}

/* Calculates the Deviance Score (DS) for each distribution
 * based on McCullagh and Nelder (1989) p. 34
 * y = N by 1 vector of responses
 * mu = N ny 1 vector of means */

wsReal getdeviance(wsVec y, wsVec mu, wsUint N_sz, wsUint datatype)
{
	wsReal R_ret = W0;

	if (datatype == 1) /* normal - identity  */
		R_ret = sseVsVsum(y, mu, N_sz, 1);
		// DS = sum((y - mu).^2);
	else if (datatype == 2) {
		LOOP (i, N_sz) { /* binary - logit  */
			if (y[i] == WISARD_UNAFFECTED)
				R_ret += log(W1/((W1-mu[i]))+ 1e-5);
			else if (y[i] == WISARD_AFFECTED)
				R_ret += log(W1/(mu[i] + 1e-5));
		}
		R_ret *= W2;
	} else if (datatype == 3) { /* poisson - log */
		LOOP (i, N_sz)
			R_ret += y[i]*log(y[i]/mu[i]) - (y[i] - mu[i]);
		R_ret *= W2;
		// DS = 2*sum(y.*log(y./mu) - (y - mu));
	} else if (datatype == 4) { /* gamma - inverse  */
		LOOP (i, N_sz)
			R_ret += -log(y[i]/mu[i]) + (y[i]-mu[i])/mu[i];
		R_ret *= W2;
		// DS = 2*sum(-log(y./mu) + (y-mu)./mu);
	} else if (datatype == 5) LOOP (i, N_sz) { /* inverse Gaussian */
		wsReal R_tmp = y[i] - mu[i];
		R_ret += SQR(R_tmp)/(SQR(mu[i]) * y[i]);
		// DS = sum((y - mu).^2/(mu.^2.*y)); 
	}
	return R_ret;
}

/* !!! PHARAOH-specific analysis !!! */
#if ((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE & FUNCTION_PHARAOH)

/*
 * M_Xt : varsum * N_sample
 */
cVector mldivide(cStdMatrix& M_1, cVector& V_2)
{
	cSymMatrix	M_1t1	= M_1.tM();
	cVector		V_1t2	= M_1 * V_2;
	cSymMatrix&	M_1t1i	= M_1t1.inv();
	if (M_1t1i.row() == 0) {
		delete& M_1t1i;
		cSymMatrix	M_1t1e	= M_1t1.einv();
		if (M_1t1e.row() == 0)
			return cVector();
		cVector		V_ret	= M_1t1e * V_1t2;
		V_ret.setDontDealloc();
		return cVector(V_ret.get(), V_ret.size());
	} else {
	//	M_1t1.rem();
		cVector		V_ret	= M_1t1i * V_1t2;
		V_ret.setDontDealloc();
		delete& M_1t1i;
		return cVector(V_ret.get(), V_ret.size());
	}
}

cVector mldivide(cSymMatrix& M_1, cVector& V_2)
{
	cSymMatrix	M_1t1	= M_1.Mt();
	cVector		V_1t2	= M_1 * V_2;
	cSymMatrix&	M_1t1i	= M_1t1.inv();
	if (M_1t1i.row() == 0) return cVector();

//	M_1t1.rem();
	cVector		V_ret	= M_1t1i * V_1t2;
//	for (wsUint i=0 ; i<V_ret.size() ; i++)
//		if (fabs(V_ret.get()[i]) > 10000.0)
//			printf("Too big\n");
	V_ret.setDontDealloc();
	delete& M_1t1i;
	return cVector(V_ret.get(), V_ret.size());
}

cStdMatrix mldivide(cSymMatrix& M_1, cStdMatrix& M_2)
{
	cSymMatrix	M_1t1	= M_1.Mt();
	cSymMatrix&	M_1t1i	= M_1t1.inv();
	if (M_1t1i.row() == 0) return cStdMatrix();

	//	M_1t1.rem();
	cStdMatrix	M_1t2	= M_1 * M_2;
	cStdMatrix	M_ret	= M_1t1i * M_1t2;
	M_ret.setClean(MATDEL_NONE);
	delete& M_1t1i;
	return cStdMatrix(M_ret.row(), M_ret.col(), M_ret.get());
}

int do_PHARAOH(char B_overlap, wsReal R_lambdaG, wsReal R_lambdaP,
	cStdMatrix& M_Xt, cStdMatrix& M_Yt, vvInt& M_phenoRel,
	wsUint N_anaSamp, vInt& Nv_dist, wsUint N_mnfs,
	vStr Sv_latent, mvStr& Sv_mnfsInLatent, cStdMatrix& M_W0,
	vvInt& windex, wsVec Ra_retA, wsVec Ra_retW, wsReal* Rp_dev/*=NULL*/,
	wsUint N_maxIter/*=100*/, wsVec Ra_seA/*=NULL*/)
{
	static int KKKK = 0;
	KKKK++;
	cTimer t4;
	if (M_Xt.col() != N_anaSamp) halt("SYSERR: Dimension error M_Xt");
	if (M_Xt.row() != N_mnfs) halt("SYSERR: Dimension error M_Xt");
	wsMat		Ra_W0		= M_W0.get();
	wsUint		N_latent	= (wsUint)Sv_mnfsInLatent.size();
	wsUint		N_pheno		= M_Yt.row();

//	windex		= find(W0 == 99)
	cStdMatrix	M_W(N_mnfs, N_latent);			/* [Ng*np] */
	wsMat		Ra_W		= M_W.get();			/* [Ng*np] */
	LOOP (i, N_latent) FOREACH (vInt_it, windex[i], j) {
#ifdef _DEBUG
		if (Ra_W0[*j][i] != GESCA_Winit)
			halt("SYSERR: Indexing error");
		if (*j >= (int)N_mnfs || *j < 0 || i >= N_latent) {
			halt("Index overflow error");
		}
#endif
		Ra_W[*j][i]	= wsUnifrand();
	}
#ifdef PHARAOHvalidate
	char S_fbuf[512];
	sprintf(S_fbuf, "gsca.W0.%d", KKKK);
	M_W.file(S_fbuf);
#endif
	// vecW	= W[windex]
	wsUint N_sumWindex = 0;
	FOREACH (vvInt_it, windex, i) N_sumWindex += (wsUint)(i->size());
	cVector		vecW(N_sumWindex);
//	VEC_t		Ra_vecW		= vecW.get();
	// W[windex]	= rand(num_windex,1)
	// 
	/* [Np*Ns] = [Ng*Np]' %*% [Ng*Ns]  */
	cStdMatrix	M_Ft	= M_W.tM(M_Xt);
	wsMat		Ra_Ft	= M_Ft.get();		/* [Np*Ns] */
// F		= X %*% W

#ifdef PHARAOHvalidate
	sprintf(S_fbuf, "gsca.F0.%d", KKKK);
	M_Ft.file(S_fbuf);
#endif

	/* Normalize F matrix */
	// F[,j] = sqrt(ncase) * F[,j]/norm(F[,j])
	LOOP (j, N_latent)
		sseVpC(Ra_Ft[j], sqrt(N_anaSamp/sseVV(Ra_Ft[j], N_anaSamp)),
			Ra_Ft[j], N_anaSamp);
	
	// A = (t(F) %*% F) %l% t(F) %*% y
 	cSymMatrix		M_FF		= M_Ft.Mt();	/* [Np*Np] t(F) %*% F */
	cSymMatrix&		M_FFi		= M_FF.inv();
	cStdMatrix		M_FFF		= M_FF * M_Ft;
	delete &M_FFi;
	cVector			V_alpha(N_pheno);
	wsVec			Ra_alpha	= V_alpha.get();
	vector<wsVec>	Vv_A;
	vector<wsVec>	Vv_z;
	vector<wsVec>	Mv_V;

	LOOP (i, N_pheno) {
		cVector V_yi	= M_Yt.r2v_ptr(i);
		cVector	V_Fyi	= M_Ft * V_yi;	/* [Np*1] = [Np*Ns] %*% [Ns*1] */
		// cVector		V_A		= mldivide(M_FF, V_Fy);
		cVector V_Ai	= M_FFF * V_yi;
		if (V_Ai.size() == 0) {
			if (Rp_dev) *Rp_dev = WISARD_NAN;
			return -1;
		}
		V_Ai.setDontDealloc();
		// alpha  = mean(y)                  # intercept
		Ra_alpha[i]		= V_yi.mean();			/* [1] */
#ifdef PHARAOHvalidate
		LOGnote("Initial intercept for [%d]th phenotype [%g]\n", i, Ra_alpha[i]);
#endif

#ifdef PHARAOHvalidate
		/* Copy variable */
		wsVec		Ra_Ai	= V_Ai.get();		/* [Np] */
		sprintf(S_fbuf, "gsca.A.phen%d.%d", i, KKKK);
		V_Ai.file(S_fbuf);
#endif

//		cDiagMatrix	V;
		cVector		V_z;
		// lp     = alpha + F %*% A          # lp = linear predictor (eta)
		cVector		V_lp	= V_Ai * M_Ft;	/* [Ns] */
		V_lp += Ra_alpha[i];
		if (V_lp.size() != N_anaSamp) halt("SYSERR: Unexpected V_lp size");

#ifdef PHARAOHvalidate
		sprintf(S_fbuf, "gsca.lp0.phen%d.%d", i, KKKK);
		V_lp.file(S_fbuf);
#endif
		// mu     = getmu(lp,dist)
		cVector		V_mu	= getmu(V_lp, Nv_dist[i]);	/* [Ns] */
#ifdef PHARAOHvalidate
		sprintf(S_fbuf, "gsca.mu0.phen%d.%d", i, KKKK);
		V_mu.file(S_fbuf);
#endif

		// z      = lp + (y - mu)/v;             # adjusted DV
		cVector		vo		= getvmu(V_mu, Nv_dist[i]);	/* [Ns] */
#ifdef PHARAOHvalidate
		sprintf(S_fbuf, "gsca.v0.phen%d.%d", i, KKKK);
		vo.file(S_fbuf);
#endif

		// V      = diag(as.vector(v))
//		V.init(N_anaSamp, vo.get());					/* [Ns*Ns] */
		if (fabs(vo.get()[0]) > 1e5 && OPT_ENABLED(verbose))
			printf("Too large V!\n");
//		V.setDontDealloc();
		vo.setDontDealloc();

		// z      = lp + (y - mu)/v;             # adjusted DV
		V_mu.subFrom(V_yi);
		V_mu /= vo;
		V_mu += V_lp;
		V_mu.setDontDealloc();
		V_z.init(N_anaSamp, V_mu.get());				/* [Ns] */
		V_z.setDontDealloc();

#ifdef PHARAOHvalidate
		sprintf(S_fbuf, "gsca.z0.%d", KKKK);
		V_z.file(S_fbuf);
#endif

		Vv_A.push_back(V_Ai.get());
		Vv_z.push_back(V_z.get());
		Mv_V.push_back(vo.get());
	} /* for(N_pheno) */

	/*
	 *
	 * IRLS Algorithm
	 *
	 */
	wsUint		it			= 0;
	wsReal		ceps		= OPT_REAL(prothr);
	if (ceps < 1e-8)
		halt("Too strict or invalid PHARAOH threshold [%g] given", ceps);
// 		est_new	= matrix(c(vecW, AA),ncol=1)
// 		est_old	= zeros(nrow(est_new),ncol(est_new))
	wsReal		R_crit		= W0;
	char		B_lasso		= OPT_ENABLED(lasso);

	/* Size of estimators = N_mnfs + #pheno(intercepts) + N_latent*N_pheno */
	wsUint		N_sumLatent	= N_latent*N_pheno;
	cVector V_estOld, V_estNew;
	V_estOld.init((wsUint)(N_sumWindex + N_pheno + N_sumLatent));
	V_estNew.init((wsUint)(N_sumWindex + N_pheno + N_sumLatent));
	// AA      = rbind(alpha, A)
	wsVec	Ra_estOld = V_estOld.get();

	/* Copy initial values */
	memcpy(Ra_estOld, vecW.get(), sizeof(wsReal)*N_sumWindex);
	LOOP (i, N_pheno)
		Ra_estOld[N_sumWindex+i] = W1;
	for (wsUint i=0,j=0 ; i<N_pheno ; i++,j+=N_pheno)
		memcpy(Ra_estOld+N_sumWindex+N_pheno+j, Vv_A[i], sizeof(wsReal)*N_latent);

	do {
//		LOG("%g\n", R_alpha);
//#ifdef _DEBUG
		if (it) {
			if (R_lambdaG == R_lambdaP)
				lverbose("lambda [%g] crit [%g] [%s]\n", R_lambdaG,
					R_crit, t4.getReadable());
			else
				lverbose("lambda [%g,%g] crit [%g] [%s]\n", R_lambdaG,
					R_lambdaP, R_crit, t4.getReadable());
		}
//#endif
		t4.start();

		/* Check iteration limit */
		if (it++ == N_maxIter) {
			LOOP (i, N_pheno) Ra_alpha[i] = WISARD_NAN;
			it = 0;
			break;
		}

		/* Step 1: Update W */
//		int	kk		= -1;
		LOOP (j, N_latent) {
			vInt&	wIdxj	= windex[j];
			wsUint	Nj		= (wsUint)Sv_mnfsInLatent[Sv_latent[j]].size();
			if (Nj == 0)
				continue;

			// k               = kk + 1
//			wsUint k		= kk + 1;
			// kk              = kk + Nj
//			kk				+= Nj;

			// X_j             = X[,k:kk]
			cStdMatrix	M_Xjt;
//			if (B_overlap)
				M_Xjt = M_Xt.subsetPtr(wIdxj);
//			else
//				M_Xjt = M_Xt.subsetPtr(k, kk);

			// 			MAT_t	Ra_Xj	= NULL; /* [Nj*Ns] */
			// 			wsAlloc(Ra_Xj, VEC_t, Nj);
			// 			for (int K=k,L=0 ; K<=kk ; K++,L++)
			// 				Ra_Xj[L] = Ra_Xt[K];
			// 			cStdMatrix	M_Xjt(Nj, N_samp, Ra_Xj);		/* [Nj*Ns] */
			// 			M_Xjt.setDontDealloc();

			cVector V_Wj;
			LOOP (i, N_pheno) {
				if (!M_phenoRel[i][j])
					continue;
				wsVec	Ra_A	= Vv_A[i];

				// za = z - alpha
				cVector V_z(Vv_z[i], N_anaSamp, 1);
				cVector V_za	= V_z - Ra_alpha[i];	/* [Ns] */

				// FHA             = F %*% H %*% A
// 				VEC_t	Ra_old	= Ra_Ft[j];		/* [Ns] */
				wsReal	R_v		= Ra_A[j];
				Ra_A[j]			= W0;
//				Ra_Ft[j]		= Ra_0;
//				cVector	V_FjA	= V_A * M_Ft;	/* [Ns] */
				cVector V_A(Vv_A[i], N_latent, 1);
				cVector	V_FjA	= V_A * M_Ft;	/* [Ns] */
				// z1              = za - FHA
				cVector	V_z1	= V_za - V_FjA;	/* [Ns] */
				Ra_A[j]			= R_v;
//				Ra_Ft[j]		= Ra_old;
#ifdef PHARAOHvalidate
//				lverbose("[%d] range [%d - %d]\n", j, k, kk);
#endif

				// A1             = A[j]^2 * (t(X_j) %*% V %*% X_j) + lambda * eye(Nj)
				cDiagMatrix M_V(N_anaSamp, Mv_V[i], 1);
				cStdMatrix	M_XjtV	= M_Xjt * M_V;		/* [Nj*Ns] = [Nj*Ns]*[Ns*Ns] */
				cSymMatrix	M_A1	= M_XjtV.Mts(M_Xjt);	/* [Nj*Nj] */
				M_A1 *= SQR(R_v);
				M_A1.addDiag(R_lambdaG);

				// A2             = A[j] * t(X_j) %*% V %*% z1
				cVector	V_A2 = M_XjtV * V_z1;				/* [Nj] */
				V_A2 *= R_v;
				// w_j       = A1 %l% A2
				cVector		V_Wjk	= mldivide(M_A1, V_A2);	/* [Nj] */
				if (V_Wjk.size() == 0) {
					if (Rp_dev) *Rp_dev = WISARD_NAN;
					return -1;
				}
//				w_j       = sqrt(ncase) * w_j/norm(X_j %*% w_j)
				cVector V_XjWj = V_Wjk * M_Xjt;
				V_Wjk *= sqrt(N_anaSamp / V_XjWj.ss());

				if (V_Wj.size())
					V_Wj += V_Wjk;
				else
					V_Wj = V_Wjk;
			}
			//				W[k:kk,j]       = w_j
//			if (B_overlap) {
				wsUint L = 0;
				FOREACHDO(vInt_it, windex[j], K, L++) {
#ifdef _DEBUG
					if (*K >= (int)N_mnfs || *K < 0 || j >= N_latent) halt("Index overflow error");
#endif
					Ra_W[*K][j] = V_Wj.get()[L];
				}
//			} else for (int K=k, L=0 ; K <= kk ; K++, L++)
//				Ra_W[K][j] = V_Wj.get()[L];

			//				F[,j]           = X_j %*% w_j
			cVector	V_Xwj2 = V_Wj * M_Xjt;
			memcpy(Ra_Ft[j], V_Xwj2.get(), sizeof(wsReal)*N_anaSamp);
		}

//			sseFree(Ra_0);
#ifdef PHARAOHvalidate
		sprintf(S_fbuf, "gsca.F%d.%d", it, KKKK);
		M_Ft.file(S_fbuf);
#endif
		/* Step 2: Update A */
		// Q       = cbind(ones(ncase, 1), F)
		wsMat Ra_Qt = NULL;
		wsAlloc(Ra_Qt, wsVec, N_latent+1);
		Ra_Qt[0] = sseVector(N_anaSamp);
		sseVinit(Ra_Qt[0], N_anaSamp, W1);
		for (wsUint i=1 ; i<=N_latent ; i++)
			Ra_Qt[i] = Ra_Ft[i-1];
		cStdMatrix	M_Qt(N_latent+1, N_anaSamp, Ra_Qt);	/* [Np1*Ns] */
		M_Qt.setClean();
#ifdef PHARAOHvalidate
		sprintf(S_fbuf, "gsca.Q%d.%d", it, KKKK);
		M_Qt.file(S_fbuf);
#endif
		LOOP(i, N_pheno) {
			wsVec		V = Mv_V[i];
			// tQV     = t(Q) %*% V
			cDiagMatrix M_V(N_anaSamp, V, 1);
			cStdMatrix	M_tQV = M_Qt * M_V;					/* [Np1*Ns] */
			// AA1     = tQV %*% Q + lambda*eye(ndset + 1)		/* [Np1*Np1] */
			cSymMatrix	M_AA1(sym_sseMpMt(M_tQV.get(), M_tQV.row(), M_tQV.col(),
				M_Qt.get(), M_Qt.row(), M_Qt.col()), M_tQV.row());
			M_AA1.addDiag(R_lambdaP);
#ifdef PHARAOHvalidate
			sprintf(S_fbuf, "gsca.AA1.phen%d.%d.%d", i, it, KKKK);
			M_AA1.file(S_fbuf);
#endif
			// AA2     = tQV %*% z
			cVector		V_z(Vv_z[i], N_anaSamp, 1);
			cVector		V_AA2		= V_z.Mt(M_tQV);		/* [Np1] */
#ifdef PHARAOHvalidate
			sprintf(S_fbuf, "gsca.AA2.phen%d.%d.%d", i, it, KKKK);
			V_AA2.file(S_fbuf);
#endif
			// AA      = AA1 %l% AA2
			cVector		V_AA		= mldivide(M_AA1, V_AA2); /* [Np1] */
#ifdef PHARAOHvalidate
			sprintf(S_fbuf, "gsca.AA.phen%d.%d.%d", i, it, KKKK);
			V_AA.file(S_fbuf);
#endif
			if (V_AA.size() != (N_latent + 1)) {
				if (Rp_dev) *Rp_dev = WISARD_NAN;
				sseFree(Ra_Qt[0]);
				DEALLOC(Ra_Qt);
				return -1;
			}
			// alpha   = AA[1]  # New intercept
			// wsReal R_alpDif = fabs(Ra_alpha[i] - V_AA.get()[0]);
			Ra_alpha[i] = V_AA.get()[0];
			/* Too large? */
//			if (Ra_alpha[i] > 10e8) {
//				printf("Too large crit?\n");
//			}
			// A       = AA[2:nrow(AA),]
			memcpy(Vv_A[i], V_AA.get() + 1, sizeof(wsReal)*N_latent);

			/* Step 3: Update V and z */
			// lp      = Q %*% AA
			cVector V_lp	= V_AA * M_Qt;							/* [Ns] */

			// mu      = getmu(lp, dist)
			cVector V_newZ	= getmu(V_lp, Nv_dist[i]);

#ifdef PHARAOHvalidate
			sprintf(S_fbuf, "gsca.mu%d.%d", it, KKKK);
			V_newZ.file(S_fbuf);
#endif
			// v       = getvmu(mu, dist)
			// V       = diag(as.vector(v))
			cVector v		= getvmu(V_newZ, Nv_dist[i]);
			memcpy(Mv_V[i], v.get(), sizeof(wsReal)*N_anaSamp);
#ifdef PHARAOHvalidate
			sprintf(S_fbuf, "gsca.v%d.%d", it, KKKK);
			v.file(S_fbuf);
#endif

			// z       = (y - mu)/v + lp             # adjusted DV
			cVector V_y = M_Yt.r2v_ptr(i);
			V_newZ.subFrom(V_y);
			V_newZ /= v;
			V_newZ += V_lp;
			memcpy(Vv_z[i], V_newZ.get(), sizeof(wsReal)*N_anaSamp);

#ifdef PHARAOHvalidate
			sprintf(S_fbuf, "gsca.z%d.%d", it, KKKK);
			V_newZ.file(S_fbuf);
#endif
		}

#ifdef PHARAOHvalidate
		char S_fbuf[512];
		sprintf(S_fbuf, "gsca.W%d.%d", it, KKKK);
		M_W.file(S_fbuf);
#endif
		wsVec Ra_vecW = vecW.get();
		wsUint I = 0;
		LOOP (i, N_latent) FOREACHDO (vInt_it, windex[i], j, I++) {
#ifdef _DEBUG
			if (I == vecW.size()) halt("SYSERR: Debug fail");
#endif
#ifdef _DEBUG
			if (*j >= (int)N_mnfs || i >= N_latent) {
				halt("Index overflow error");
			}
//			halt("Index overflow error");
#endif
			Ra_vecW[I] = Ra_W[*j][i];
		}
		// vecW    = W[windex]
		wsVec Ra_estNew = V_estNew.get();
		memcpy(Ra_estNew, Ra_vecW, sizeof(wsReal)*vecW.size());
		memcpy(Ra_estNew+vecW.size(), Ra_alpha, N_pheno);
		LOOP (i, N_pheno)
			memcpy(Ra_estNew+vecW.size()+N_pheno+i*N_latent, Vv_A[i], sizeof(wsReal)*N_latent);
		//est_new = matrix(c(vecW, AA),ncol=1)

#ifdef PHARAOHvalidate
		char S_fmt[32];
		sprintf(S_fmt, "gsca.e%d.%d", it, KKKK);
		V_estOld.file(S_fmt);
#endif
//			cVector V_diff = V_estNew - V_estOld;
		R_crit = V_estNew.asum(V_estOld) ;
		if (NA(R_crit)) {
			if (Rp_dev) *Rp_dev = WISARD_NAN;
			return -1;
		}
		memcpy(V_estOld.get(), V_estNew.get(), sizeof(wsReal)*V_estNew.size());
		/* Store AA and W */
		if (R_crit <= ceps) {
#if 0
			if (Ra_seA) {
				cStdMatrix	M_repQt = M_Qt.rRepeat(N_pheno);
				cDiagMatrix	M_repV(M_)
				wsVec		V = Mv_V[i];


				LOOP (i, N_pheno) {
					wsVec		V = Mv_V[i];
					// tQV     = t(Q) %*% V
					cDiagMatrix M_V(N_anaSamp, V, 1);
					cStdMatrix	M_tQV = M_Qt * M_V;					/* [Np1*Ns] */
				}
				// -t(Q) %*% V %*% Q - lambda2 * G2
				M_Qt.MMt(MvV)
			}
#endif
			if (Ra_retA) LOOP (i, N_pheno)
				memcpy(Ra_retA + i*N_latent, Vv_A[i], sizeof(wsReal)*N_latent);
			if (Ra_retW) memcpy(Ra_retW, Ra_vecW, sizeof(wsReal)*vecW.size());
		}

		sseFree(Ra_Qt[0]);
		DEALLOC(Ra_Qt);
	} while (R_crit > ceps);

	/* Compute deviance if required */
	if (Rp_dev) {
		cStdMatrix	M_newFt	= M_W.tM(M_Xt);
		*Rp_dev = W0;
		LOOP (i, N_pheno) {
			cVector V_y = M_Yt.r2v_ptr(i);
		/* Deviance computation relies on the condition we are going to satisfy */
//		alpha       = ret$A[1]
//		A           = ret$A[-1]
//		Wn[IwN]     = ret$W
			cVector		V_a(Vv_A[i], N_latent, 1);
			cVector		V_lp	= V_a * M_newFt;
			V_lp += Ra_alpha[i];
			cVector		V_mu	= getmu(V_lp, Nv_dist[i]);
			*Rp_dev		+= -getqi(V_y.get(), V_mu.get(), V_mu.size(), Nv_dist[i]) * W2 + 1e-5;
		}

//		if (NA(*Rp_dev) && it)
//			halt("ERROR");

// 		Ft			= X %*% Wn
// 		lpt			= alpha + Ft %*% A
// 		mut			= getmu(lpt, dist)
// 		dev			= -2*getqi(y, mut, dist) + 1E-5
// 		Deviance	= Deviance + dev
// 		It          = It + ret$it
//		Ok          = Ok + 1

	}

	FOREACH (vector<wsVec>::iterator, Vv_A, i) sseFree(*i);
	FOREACH (vector<wsVec>::iterator, Vv_z, i) sseFree(*i);
	FOREACH (vector<wsVec>::iterator, Mv_V, i) sseFree(*i);

	return it;
}

cPharaohAnalysis::cPharaohAnalysis(cIO* Cp_inpIO, cSetManagerAnalysis* Cp_inpAnaSM,
	cPPPAnalysisV2* Cp_inpAnaPPP)
	: cAnalysis(Cp_inpIO)
{
	wsUintCst	NFILTV = 0;
	wsUintCst	NFILTG = 1;
	/* Option check */
	if (OPT_ENABLED(pharaoh) && OPT_ENABLED(proopt))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--pharaoh", "--proopt");
	/* Variable check */
	if (!Cp_inpAnaSM)
		halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--pharaoh", "no set definition");

	Cp_anaSM	= Cp_inpAnaSM;
	Cp_anaPPP	= Cp_inpAnaPPP;

	/* Choose only for #genes > 1 */
	mGeneSet&	Xm_oriGS	= Cp_anaSM->getGeneSetDef();
	mGeneDef&	Xm_oriGD	= Cp_anaSM->getGeneDef();
	if (Xm_oriGS.size() == 0 || Xm_oriGD.size() == 0) {
		if (!OPT_ENABLED(gesca) || !OPT_ENABLED(gxg))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, OPT_ENABLED(gesca) ? "--gesca" : "--pharaoh", "no available gene or gene-set");
		else if (OPT_ENABLED(gesca) && OPT_ENABLED(gxg) && Xm_oriGD.size() == 0)
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--gesca with --gxg", "no available gene");
	}
//	vVariant&	Xv_vrt		= Cp_IO->getVariants();
	xMaf*		Xp_maf		= Cp_IO->getMAF();
//	mDataIdx	Xm_gene2storePos;
//	wsUint I=0;
// 	FOREACHDO (mGeneDef_it, Xm_gdef, it, I++) {
// 		vInt	Xa_set;
// 		wsUint	N_mac = 0;
//
// 		Xm_gene2storePos[it->first] = I;
// 	}

	if (OPT_ENABLED(gesca) && OPT_ENABLED(gxg)) {
		FOREACH(mGeneDef_it, Xm_oriGD, i)
			_initGeneDef(Xm_oriGD, i->first, Xp_maf, NFILTV);
	} else {
		cTimer t3;
		LOGnote("Collapsing genotypes for [%d] pathways\n", Xm_oriGS.size());
		t3.start();
		FOREACH (mGeneSet_it, Xm_oriGS, i) {
			vStr	X_newGS;
			vStr&	Xv_genes	= i->second;
			/* Do not check if #gene <= NFILT */
			if (Xv_genes.size() <= NFILTG && !IS_ASSIGNED(expression)) continue;

			/* For each gene, check integrity */
			FOREACH (vStr_it, Xv_genes, j) {
				if (!_initGeneDef(Xm_oriGD, *j, Xp_maf, NFILTV)) continue;

				/* Insert new entry */
				X_newGS.push_back(*j);
			}
			/* If # of genes are not enough pass */
			if (X_newGS.size() <= NFILTG) continue;

			/* Insert new entry */
			Xm_gs.insert(make_pair(i->first, X_newGS));
		}
		LOG("[%d] analyzable gene-sets found with [%d] genes, took [%s]\n",
			Xm_gs.size(), Xm_gd.size(), t3.getReadable());
	}
}

cPharaohAnalysis::~cPharaohAnalysis()
{

}

int thr_optLambda_PHARAOH(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	xAnaThread*	Xp_dt		= (xAnaThread *)Vp_shareData;
//	cIO*		Cp_IO		= Xp_dt->Cp_IO;
	xGscaAnaThread*
				Xp_gd		= (xGscaAnaThread *)(Xp_dt->Vp_data);
	cStdMatrix&	M_Xt		= *(Xp_gd->Mp_Xt);
	cStdMatrix&	M_W0		= *(Xp_gd->Mp_W0);
	wsUint		N_loop		= (wsUint)Xp_gd->Rv_lambdas.size();

	for (wsUint j=N_idx ; j<N_loop ; j+=OPT_NUMBER(thread)) {
//		wsReal	R_curDev	= W0;
		wsVec	Ra_curDev	= sseVector(Xp_gd->N_cv);
		LOOP (k, Xp_gd->N_cv) {
			Ra_curDev[k] = doCV_PHARAOH(Xp_gd->B_overlap,
				M_Xt, Xp_gd->M_origYt, Xp_gd->M_phenoRel,
				Xp_gd->Na_samp[k], Xp_gd->Nv_dist, Xp_gd->Sv_idx2mnfs,
				Xp_gd->Sv_latent, Xp_gd->Sv_mnfsInLatent, M_W0,
				Xp_gd->windex,
				Xp_gd->Na_indices[k], Xp_gd->Rv_lambdas[j], Xp_gd->Na_iters[j]+k);
			if (OPT_ENABLED(verbose)) {
				if (NA(Ra_curDev[k])) {
					if (Xp_gd->Rv_lambdas[j].size() == 1)
						LOGwarn("Lambda [%g] CV [%d] failed\n", Xp_gd->Rv_lambdas[j][0],
							k+1);
					else
						LOGwarn("Lambda [%g,%g] CV [%d] failed\n", Xp_gd->Rv_lambdas[j][0],
							Xp_gd->Rv_lambdas[j][1], k+1);
				} else {
					if (Xp_gd->Rv_lambdas[j].size() == 1)
						LOGwarn("Lambda [%g] CV [%d] dev [%g] iter [%d]\n",
							Xp_gd->Rv_lambdas[j][0], k+1,
							Ra_curDev[k], Xp_gd->Na_iters[j][k]);
					else
						LOGwarn("Lambda [%g,%g] CV [%d] dev [%g] iter [%d]\n",
							Xp_gd->Rv_lambdas[j][0], Xp_gd->Rv_lambdas[j][1], k+1,
							Ra_curDev[k], Xp_gd->Na_iters[j][k]);
				}
			}
		}

		/* Take the mean of deviances from K cvs */
		wsUint N_avail = 0;
		wsReal R_sum = W0;
		LOOP (k, Xp_gd->N_cv) {
			if (NA(Ra_curDev[k])) continue;
			R_sum += Ra_curDev[k];
			N_avail++;
		}
		if (N_avail == Xp_gd->N_cv || (Xp_gd->N_cv>4 && ((N_avail+1) == Xp_gd->N_cv)))
			Xp_gd->Ra_devs[j] = R_sum / (wsReal)N_avail;
		else
			Xp_gd->Ra_devs[j] = WISARD_NAN;

		/* Take the mean of deviances from K cvs */
//		Xp_gd->Ra_devs[j] = sseVsum(Ra_curDev, Xp_gd->N_cv) / (wsReal)Xp_gd->N_cv;
	}

	return 0;
}

int thr_optLambda_GESCA(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	xAnaThread*	Xp_dt		= (xAnaThread *)Vp_shareData;
	//	cIO*		Cp_IO		= Xp_dt->Cp_IO;
	xGescaCvThread*
				Xp_gc		= (xGescaCvThread *)(Xp_dt->Vp_data);
	cStdMatrix&	M_Xt		= *(Xp_gc->Mp_Xt);
	cStdMatrix&	M_W0		= Xp_gc->M_W0;
	wsUint		N_loop		= (wsUint)Xp_gc->Rv_lambdas.size();

	for (wsUint j=N_idx ; j<N_loop ; j+=OPT_NUMBER(thread)) {
		//		wsReal	R_curDev	= W0;
		wsVec	Ra_curDev	= sseVector(Xp_gc->N_cv);
		LOOP(k, Xp_gc->N_cv) {
			Ra_curDev[k] = doCV_GESCA(Xp_gc->B_formative, Xp_gc->M_Yt,
				M_Xt, M_W0, Xp_gc->M_B0_t, Xp_gc->M_C0_t, Xp_gc->windex, Xp_gc->bindex, Xp_gc->cindex,
				Xp_gc->N_obs, Xp_gc->N_ltt, Xp_gc->Na_samp[k], Xp_gc->Sv_idx2mnfs, Xp_gc->Sv_latent, Xp_gc->Sv_mnfsInLatent,
				Xp_gc->Na_idxCV[k], Xp_gc->Rv_lambdas[j], Xp_gc->Na_iters[j] + k);
			if (OPT_ENABLED(verbose)) {
				if (NA(Ra_curDev[k])) {
					if (Xp_gc->Rv_lambdas[j].size() == 1)
						LOGwarn("Lambda [%g] CV [%d] failed\n", Xp_gc->Rv_lambdas[j][0],
						k+1);
					else
						LOGwarn("Lambda [%g,%g] CV [%d] failed\n", Xp_gc->Rv_lambdas[j][0],
						Xp_gc->Rv_lambdas[j][1], k+1);
				}
				else {
					if (Xp_gc->Rv_lambdas[j].size() == 1)
						LOGwarn("Lambda [%g] CV [%d] dev [%g] iter [%d]\n",
						Xp_gc->Rv_lambdas[j][0], k+1,
						Ra_curDev[k], Xp_gc->Na_iters[j][k]);
					else
						LOGwarn("Lambda [%g,%g] CV [%d] dev [%g] iter [%d]\n",
						Xp_gc->Rv_lambdas[j][0], Xp_gc->Rv_lambdas[j][1], k+1,
						Ra_curDev[k], Xp_gc->Na_iters[j][k]);
				}
			}
		}

		/* Take the mean of deviances from K cvs */
		wsUint N_avail = 0;
		wsReal R_sum = W0;
		LOOP(k, Xp_gc->N_cv) {
			if (NA(Ra_curDev[k])) continue;
			R_sum += Ra_curDev[k];
			N_avail++;
		}
		if (N_avail == Xp_gc->N_cv || (Xp_gc->N_cv>4 && ((N_avail+1) == Xp_gc->N_cv)))
			Xp_gc->Ra_devs[j] = R_sum / (wsReal)N_avail;
		else
			Xp_gc->Ra_devs[j] = WISARD_NAN;

		/* Take the mean of deviances from K cvs */
		//		Xp_gd->Ra_devs[j] = sseVsum(Ra_curDev, Xp_gd->N_cv) / (wsReal)Xp_gd->N_cv;
	}

	return 0;
}

inline void sub_optLambda_PHARAOH(char B_overlap, cStdMatrix& M_origYt, cStdMatrix& M_Xt, cStdMatrix& M_W0,
	vvInt& windex, mDataIdx& Xm_mnfs2idx, vStr& Sv_idx2mnfs,
	vStr& Sv_latent, mvStr& Sv_mnfsInLatent, vvInt& M_phenoRel, wsUint* Na_samp, wsUint** Na_indices,
	vInt& Nv_dist, wsUint N_cv,
	vector<vReal>& Rv_lambdas,
	wsUint N_sameLambda,
	wsVec Ra_devs, wsUint** Na_iters)
{
	wsUint L = 0;
	FOREACHDO (vLambda_it, Rv_lambdas, l, L++) {
		wsVec Ra_curDev = sseVector(N_cv);
		LOOP (k, N_cv) {
			Ra_curDev[k] = doCV_PHARAOH(B_overlap,
				M_Xt, M_origYt, M_phenoRel,
				Na_samp[k], Nv_dist, Sv_idx2mnfs,
				Sv_latent, Sv_mnfsInLatent, M_W0,
				windex,
				Na_indices[k], *l, Na_iters[L] + k);
		}

		/* Take the mean of deviances from K cvs */
		wsUint N_avail = 0;
		wsReal R_sum = W0;
		LOOP (i, N_cv) {
			if (NA(Ra_curDev[i])) continue;
			R_sum += Ra_curDev[i];
			N_avail++;
		}
		sseFree(Ra_curDev);
		if (N_avail != N_cv)
			Ra_devs[L] = WISARD_NAN;
		else
			Ra_devs[L] = R_sum / (wsReal)N_avail;
	}
}

inline void sub_optLambda_GESCA(char B_formative,
	cStdMatrix& M_Yt,
	cStdMatrix* Mp_Xt, cStdMatrix& M_W0, cStdMatrix& M_B0_t,
	cStdMatrix& M_C0_t, vvInt& windex, vInt& bindex,
	vvInt& cindex, wsUintCst N_obs, wsUintCst N_ltt,
	vStr& Sv_idx2mnfs, vStr& Sv_latent, mvStr& Sv_mnfsInLatent,
	wsUint N_cv,
	wsUint* Na_samp, wsUint** Na_idxCV,
	vector<vReal>& Rv_lambdas,
	wsUint N_sameLambda,
	wsVec Ra_devs, wsUint** Na_iters)
{
	wsUint	L			= 0;
	char	B_pairwise	= OPT_ENABLED(gxg) && !OPT_ENABLED(gxgall);

	FOREACHDO(vLambda_it, Rv_lambdas, l, L++) {
		wsVec Ra_curDev = sseVector(N_cv);
		LOOP(k, N_cv) {
			Ra_curDev[k] = doCV_GESCA(B_formative, M_Yt,
				*Mp_Xt, M_W0, M_B0_t, M_C0_t, windex, bindex, cindex,
				N_obs, N_ltt, Na_samp[k], Sv_idx2mnfs, Sv_latent, Sv_mnfsInLatent,
				Na_idxCV[k], *l, Na_iters[L] + k);
		}

		/* Take the mean of deviances from K cvs */
		wsUint N_avail = 0;
		wsReal R_sum = W0;
		LOOP(i, N_cv) {
			if (NA(Ra_curDev[i])) continue;
			R_sum += Ra_curDev[i];
			N_avail++;
		}
		sseFree(Ra_curDev);
		if (N_avail != N_cv && !B_pairwise)
			Ra_devs[L] = WISARD_NAN;
		else
			Ra_devs[L] = R_sum / (wsReal)N_avail;
	}
}

void setLambdaDefault(wsUint N_tryG, wsUint N_tryP, char B_separate,
	vReal& Rv_intvLambda, vLambda& Rv_lambdas)
{
	Rv_intvLambda.push_back(1);
	if (B_separate) {
		Rv_intvLambda.push_back(1);
		for (int i=0, j=0 ; j<N_lambda ; i++, j++) {
			wsReal R_v1 = (wsReal)pow(10.0, i);
			for (int k=2, l=0 ; l<N_lambda ; k++, l++) {
				vReal V_cur;
				V_cur.push_back(R_v1);
				V_cur.push_back((wsReal)pow(10.0, k));
				Rv_lambdas.push_back(V_cur);
			}
		}
	} else for (int i=1, j=0 ; j<N_lambda ; i++, j++) {
		vReal V_cur;
		V_cur.push_back((wsReal)pow(10.0, i));
		Rv_lambdas.push_back(V_cur);
	}
}

void setLambdaRange(wsReal R_s, wsReal R_e, wsUint N_try, char B_separate,
	vReal& Rv_intvLambda, vLambda& Rv_lambdas)
{
	wsReal R_s10 = log10(R_s);
	wsReal R_e10 = log10(R_e);
	wsReal R_intv = (R_e10 - R_s10) / (N_lambda - 1);

	Rv_intvLambda.push_back(R_intv);
	if (B_separate) {
		Rv_intvLambda.push_back(R_intv);
		LOOP(i, N_lambda) {
			wsReal R_v1 = (wsReal)pow(10.0, R_s10);
			wsReal R_s10_ = log10(R_s);
			LOOP(j, N_lambda) {
				vReal	V_cur;
				V_cur.push_back(R_v1);
				V_cur.push_back((wsReal)pow(10.0, R_s10_));
				Rv_lambdas.push_back(V_cur);
				R_s10_ += R_intv;
			}
			R_s10 += R_intv;
		}
	} else LOOP(i, N_lambda) {
		wsReal R_v1 = (wsReal)pow(10.0, R_s10);
		vReal	V_cur;
		V_cur.push_back(R_v1);
		Rv_lambdas.push_back(V_cur);
		R_s10 += R_intv;
	}
}

void setCVseqBinary(cVector& V_origY, wsUintCst N_samp, wsUintCst N_cv,
	wsUint** Na_idxCV, wsUint* Na_samp)
{
	lverbose("Distribution : binary\n");
	/* Find cases and controls */
	wsUint	N_case, N_ctrl;
	wsUint*	Na_idxPermCa	= NULL;
	wsUint*	Na_idxPermCt	= NULL;

	/* Count the number of cases and controls */
	getCactIndex(V_origY.get(), N_samp, &N_case, &N_ctrl, &Na_idxPermCa,
		&Na_idxPermCt, true);

	wsReal	R_caseS = W0;
	wsReal	R_ctrlS = W0;
	wsReal	R_caseItv	= N_case / (wsReal)N_cv;
	wsReal	R_ctrlItv	= N_ctrl / (wsReal)N_cv;
	LOOP(i, N_cv) {
		wsReal	R_caseE	= R_caseS + R_caseItv;
		wsReal	R_ctrlE	= R_ctrlS + R_ctrlItv;
		wsUint	N_caseS	= (wsUint)(R_caseS + REAL_CONST(0.5));
		wsUint	N_caseE	= (wsUint)(R_caseE + REAL_CONST(0.5));
		wsUint	N_ctrlS	= (wsUint)(R_ctrlS + REAL_CONST(0.5));
		wsUint	N_ctrlE	= (wsUint)(R_ctrlE + REAL_CONST(0.5));

		/* Copy indices */
		memcpy(Na_idxCV[i], Na_idxPermCa+N_caseS, sizeof(wsUint)*(N_caseE-N_caseS));
		memcpy(Na_idxCV[i]+(N_caseE-N_caseS), Na_idxPermCt+N_ctrlS, sizeof(wsUint)*(N_ctrlE-N_ctrlS));
		Na_samp[i] = (N_caseE-N_caseS) + (N_ctrlE-N_ctrlS);

		if (OPT_ENABLED(verbose)) {
			LOG("CV [%d] [%d] samples :", i, Na_samp[i]);
			LOOP(j, Na_samp[i])
				LOGnf(" %d[%g]", Na_idxCV[i][j], V_origY.get()[Na_idxCV[i][j]]);
			LOGnf("\n");
		}

		/* Set the start point to the end point */
		R_caseS = R_caseE;
		R_ctrlS = R_ctrlE;
	}
	DEALLOC(Na_idxPermCa);
	DEALLOC(Na_idxPermCt);
}

void setCVdefault(wsReal R_interval, wsUintCst N_samp, wsUintCst N_cv, wsUint** Na_idxCV, wsUint* Na_samp)
{
	lverbose("Distribution : normal\n");
	/* Mix sample and divide by #cv */
	wsUint* Na_permSeq = getPermutedSeq(N_samp);
	wsReal R_s = W0;
	LOOP(i, N_cv) {
		wsReal R_e = R_s + R_interval;
		wsUint N_s = (wsUint)(R_s + REAL_CONST(0.5));
		wsUint N_e = (wsUint)(R_e + REAL_CONST(0.5));

		/* Copy indices */
		memcpy(Na_idxCV[i], Na_permSeq+N_s, sizeof(wsUint)*(N_e-N_s));
		Na_samp[i] = N_e-N_s;
		if (OPT_ENABLED(verbose)) {
			LOG("CV [%d] [%d] samples :", i, Na_samp[i]);
			LOOP(j, Na_samp[i])
				LOGnf(" %d", Na_idxCV[i][j]);
			LOGnf("\n");
		}

		/* Set the start point to the end point */
		R_s = R_e;
	}
	DEALLOC(Na_permSeq);
}

void updateLambda(vReal& Rv_intvLambda, vLambda& Rv_lambdas,
	wsUint N_idxMin, char B_separate)
{
	vReal Rv_curopt(Rv_intvLambda.size());
	copy(Rv_lambdas[N_idxMin].begin(), Rv_lambdas[N_idxMin].end(), Rv_curopt.begin());
	Rv_lambdas.clear();
	if (B_separate) {
		wsReal R_oldIntv1 = Rv_intvLambda[0];
		wsReal R_oldIntv2 = Rv_intvLambda[1];
		wsReal R_s1 = N_idxMin == 0 ? log10(Rv_curopt[0]) : log10(Rv_curopt[0]) - R_oldIntv1 / 2;
		wsReal R_e1 = N_idxMin == (N_lambda - 1) ? log10(Rv_curopt[0]) : log10(Rv_curopt[0]) + R_oldIntv1 / 2;
		wsReal R_intv1 = (R_e1 - R_s1) / (N_lambda - 1);
		Rv_intvLambda[0] = R_intv1;
		LOOP(i, N_lambda) {
			wsReal R_s2 = log10(Rv_curopt[1]) - R_oldIntv2 / 2;
			wsReal R_e2 = log10(Rv_curopt[1]) + R_oldIntv2 / 2;
			wsReal R_intv2 = (R_e2 - R_s2) / (N_lambda - 1);
			Rv_intvLambda[1] = R_intv2;
			LOOP(j, N_lambda) {
				vReal V_cur;
				//LOG("%g,%g\n", R_s1, R_s2);
				V_cur.push_back(pow(10.0, R_s1));
				V_cur.push_back(pow(10.0, R_s2));
				Rv_lambdas.push_back(V_cur);

				/* Do single if already converged */
				if (R_intv2 < 0.01) break;
				R_s2 += R_intv2;
			}

			/* Do single if already converged */
			if (R_intv1 < 0.01) break;
			R_s1 += R_intv1;
		}
	} else {
		wsReal R_s = log10(Rv_curopt[0]) - Rv_intvLambda[0];
		wsReal R_e = log10(Rv_curopt[0]) + Rv_intvLambda[0];
		wsReal R_intv = (R_e - R_s) / (N_lambda - 1);
		LOOP(i, N_lambda) {
			vReal V_cur;
			V_cur.push_back(pow(10.0, R_s));
			Rv_lambdas.push_back(V_cur);
			R_s += R_intv;
		}
		Rv_intvLambda[0] = R_intv;
	}
}

void optLambda_PHARAOH(cIO* Cp_IO, char B_overlap, cStdMatrix& M_Xt, cStdMatrix& M_origYt,
	vvInt& M_phenoRel, vInt& Nv_dist, wsUint N_mnfs,
	mDataIdx& Xm_mnfs2idx, vStr& Sv_idx2mnfs,
	vStr& Sv_latent, mvStr& Sv_mnfsInLatent,
	cStdMatrix& M_W0, vvInt& windex,
	wsReal* Rp_lambdaG, wsReal* Rp_lambdaP)
{
	ASSERT_OPTION(cv);

	wsUintCst	N_samp		= M_origYt.col();
	wsUintCst	N_pheno		= M_origYt.row();
	vReal		Rv_intvLambda;
	vLambda		Rv_lambdas;
	wsVec		Ra_devs		= sseVector(N_lambda*N_lambda);
	wsUint		N_cv		= OPT_NUMBER(cv);
	wsReal		R_interval	= N_samp / (wsReal)N_cv;
	wsUint		N_maxSamp	= (wsUint)(R_interval * W2);

	wsUint*		Na_samp		= NULL;
	wsUint**	Na_idxCV	= NULL;
	wsAlloc(Na_idxCV, wsUint*, N_cv);
	LOOP (i, N_cv) wsAlloc(Na_idxCV[i], wsUint, N_maxSamp*2);
	wsAlloc(Na_samp, wsUint, N_cv);

	/*
	 * Create CV set
	 */
	char B_bin = 0;
	char B_cont = 0;
	FOREACH(vInt_it, Nv_dist, i) {
		if (*i == 2) B_bin = 1;
		else B_cont = 1;
	}
	if (B_cont & B_bin) {
		LOGwarn("Both binary and continuous phenotype found, CV will not consider binary phenotype's distribution");
		B_bin = 0;
	} else if (N_pheno > 1 && B_bin) {
		LOGwarn("Binary phenotypes found but their distribution will not be considered due to complexity");
		B_bin = 0;
	}

	if (!B_bin)
		setCVdefault(R_interval, N_samp, N_cv, Na_idxCV, Na_samp);
	else {
		cVector V_y = M_origYt.r2v_ptr(0);
		setCVseqBinary(V_y, N_samp, N_cv, Na_idxCV, Na_samp);
	}// else
	//	halt("SYSERR: Permutation scheme for distribution [%d] is not defined!", Nv_dist);

	/*
	 * Set initial lambda values to try
	 */
	if (IS_ASSIGNED(prorange))
		setLambdaRange(OPT_RANGE(prorange).R_s, OPT_RANGE(prorange).R_e,
			N_lambda, Rp_lambdaP?1:0, Rv_intvLambda, Rv_lambdas);
	else
		setLambdaDefault(N_lambda, N_lambda, Rp_lambdaP?1:0,
			Rv_intvLambda, Rv_lambdas);

//	wsUint		N_idxMin	= 0;
	bool		B_ok		= false;
#if defined(_DEBUG) || defined(PHARAOHopti)
	cExporter*	Cp_opti	= cExporter::summon("pharaoh.opti.res");
#endif
	do {
		/* Print current range of lambda search */
		LOG("Lambda search range [%g", Rv_lambdas[0][0]);
		if (Rv_lambdas.size() > 1) LOGnf(" ~ %g", (*(Rv_lambdas.rbegin()))[0]);
		LOGnf("]");
		if (Rp_lambdaP) {
			LOGnf(" for gene, [%g", Rv_lambdas[0][1]);
			if (Rv_lambdas.size() > 1) LOGnf(" ~ %g", (*(Rv_lambdas.rbegin()))[1]);
			LOGnf("] for pathway");
		}
		LOGnf("\n");

		vector<wsUint*> Nv_tmp;
		wsUint** Na_iters = NULL;
		wsAlloc(Na_iters, wsUint*, N_lambda*N_lambda);
		LOOP (i, N_lambda*N_lambda) wsAlloc(Na_iters[i], wsUint, N_cv);

		if (IS_MULTITHREAD) {
			vector<wsVec> X_tmp;
			vInt tmp;
			vvInt tmp2;
			xGeSCAres tmp3;
			xGscaAnaThread X_gd = {
				B_overlap,
				Nv_dist,
				N_mnfs,
				&M_Xt, &M_W0,
				M_origYt,
				Xm_mnfs2idx,
				Sv_idx2mnfs,
				Sv_latent, Sv_mnfsInLatent,
				windex,

				Rv_lambdas,
				Na_samp,

				N_cv,
				Na_idxCV,
				Ra_devs,
				Na_iters,

				NULL, NULL,
				X_tmp, X_tmp,

				NULL,
				M_phenoRel,

				NULL, NULL, tmp2, tmp, tmp3, 0, 0,
			};
			xAnaThread X_dt = { Cp_IO, &X_gd };
			WORKER().run(thr_optLambda_PHARAOH, forLambda, &X_dt, NULL, sizeof(int)*3);
		} else
			sub_optLambda_PHARAOH(B_overlap, M_origYt, M_Xt, M_W0, windex,
			Xm_mnfs2idx, Sv_idx2mnfs,
				Sv_latent, Sv_mnfsInLatent, M_phenoRel,
				Na_samp, Na_idxCV, Nv_dist, N_cv, Rv_lambdas,
				Rp_lambdaP?1:0, Ra_devs, Na_iters);

		/* Find the minimum across CV */
		wsReal	R_minDev	= Ra_devs[0];
		wsUint	N_idxMin	= 0;
		wsUint	I			= 0;
		FOREACHDO (vLambda_it, Rv_lambdas, i, I++) {
#if defined(_DEBUG) || defined(PHARAOHopti)
			Cp_opti->fmt("%g", (*i)[0]);
#endif
			LOG("Lambda [%g", (*i)[0]);
			if (Rp_lambdaP) {
#if defined(_DEBUG) || defined(PHARAOHopti)
				Cp_opti->fmt("	%g", (*i)[1]);
#endif
				LOGnf(",%g", (*i)[1]);
			}
#if defined(_DEBUG) || defined(PHARAOHopti)
			Cp_opti->fmt("	%g\n", Ra_devs[I]);
#endif
			LOGnf("] deviance [%g] iters [", Ra_devs[I]);
			LOOP (j, N_cv)
				LOGnf((j+1)==N_cv?"%d":"%d,", Na_iters[I][j]);
			LOGnf("]\n");

			if (NA(Ra_devs[I])) continue;
			if (R_minDev > Ra_devs[I] || NA(R_minDev)) {
				N_idxMin	= I;
				R_minDev	= Ra_devs[I];
			}
		}
		LOOP (i, N_lambda*N_lambda) DEALLOC(Na_iters[i]);
		DEALLOC(Na_iters);

		/* Check failure */
		if (NA(R_minDev))
			halt("All optimal candidates were failed!");

		/* Print current status */
		if (Rp_lambdaP)
			LOG("Minimum lambda [%g,%g] with deviance [%g]\n",
				Rv_lambdas[N_idxMin][0], Rv_lambdas[N_idxMin][1], R_minDev);
		else
			LOG("Minimum lambda [%g] with deviance [%g]\n", Rv_lambdas[N_idxMin][0],
				R_minDev);

		/* Set the optimal lambda */
		*Rp_lambdaG = Rv_lambdas[N_idxMin][0];
		if (Rp_lambdaP) *Rp_lambdaP = Rv_lambdas[N_idxMin][1];

		/* Set the new lambdas */
		updateLambda(Rv_intvLambda, Rv_lambdas, N_idxMin, Rp_lambdaP?1:0);

		/* Met the criterion? */
		B_ok = true;
		FOREACH(vReal_it, Rv_intvLambda, ld) {
			//LOG("Intv %g\n", *ld);
			if (*ld >= 0.01) {
				B_ok = false;
				break;
			}
		}
	} while (!B_ok);

#if defined(_DEBUG) || defined(PHARAOHopti)
	delete Cp_opti;
#endif

	sseFree(Ra_devs);
	LOOP (i, N_cv) DEALLOC(Na_idxCV[i]);
	DEALLOC(Na_idxCV);
	DEALLOC(Na_samp);
}

void optLambda_GESCA(cIO* Cp_IO, cPermuter* Cp_perm, char B_formative,
	wsUint N_perm, cStdMatrix& M_Yt, cStdMatrix& M_Xt0,
	cStdMatrix& M_Wr0, vStr& Sv_idx2ltt, vvInt& windex,
	vStr& Sv_idx2mnfs, vStr& Sv_latent, mvStr& Sv_mnfsInLatent,
	vInt& Nv_pathCoords, wsRealPtr Rp_lambdaG, wsRealPtr Rp_lambdaP)
{
	ASSERT_OPTION(cv);

	/* parameter setting */
	wsUint N_obs = (wsUint)M_Wr0.row();
	wsUint N_ltt = (wsUint)M_Wr0.col();

	wsUint		N_cov		= Cp_IO->sizeCovar();
	wsUintCst	N_anaSamp	= M_Yt.col();
	wsUintCst	N_pheno		= M_Yt.row();
	vReal		Rv_intvLambda;
	vLambda		Rv_lambdas;
	wsVec		Ra_devs		= sseVector(N_lambda*N_lambda);
	wsUint		N_cv		= OPT_NUMBER(cv);
	wsReal		R_interval	= N_anaSamp / (wsReal)N_cv;
	wsUint		N_maxSamp	= (wsUint)(R_interval * W2);

	wsUint*		Na_samp = NULL;
	wsUint**	Na_idxCV = NULL;
	wsAlloc(Na_idxCV, wsUint*, N_cv);
	LOOP(i, N_cv) wsAlloc(Na_idxCV[i], wsUint, N_maxSamp * 2);
	wsAlloc(Na_samp, wsUint, N_cv);

	/*
	* Create CV set
	*/
	setCVdefault(R_interval, N_anaSamp, N_cv, Na_idxCV, Na_samp);

	/*
	* Set initial lambda values to try
	*/
	if (IS_ASSIGNED(gxglambda)) {
		wsUint N_lambdaGrp = 0;
		wsUint N_lambda1 = 0, N_lambda2 = 0;
		wsVec Ra_lambda1 = NULL, Ra_lambda2 = NULL;
		char** Sa_lambda = loadStringValues2(OPT_STRING(gxglambda), &N_lambdaGrp, '|');

		switch (N_lambdaGrp) {
		case 1:
			Ra_lambda1 = loadRealValues(Sa_lambda[0], &N_lambda1);
			Ra_lambda2 = Ra_lambda1;
			N_lambda2 = N_lambda1;
			break;
		case 2:
			Ra_lambda1 = loadRealValues(Sa_lambda[0], &N_lambda1);
			Ra_lambda2 = loadRealValues(Sa_lambda[1], &N_lambda2);
			break;
		default:
			halt("Invalid value of --gxglambda found");
			break;
		}
		LOOP(i, N_lambda1) LOOP(j, N_lambda2) {
			vReal tmp;
			tmp.push_back(Ra_lambda1[i]);
			tmp.push_back(Ra_lambda2[j]);
			Rv_lambdas.push_back(tmp);
		}
	} else if (IS_ASSIGNED(prorange))
		setLambdaRange(OPT_RANGE(prorange).R_s, OPT_RANGE(prorange).R_e,
			N_lambda, Rp_lambdaP ? 1 : 0, Rv_intvLambda, Rv_lambdas);
	else
		setLambdaDefault(N_lambda, N_lambda, Rp_lambdaP ? 1 : 0,
			Rv_intvLambda, Rv_lambdas);

	//wsUint		N_idxMin = 0;
	bool		B_ok = false;
#if defined(_DEBUG) || defined(PHARAOHopti)
	cExporter*	Cp_opti = cExporter::summon("gesca.opti.res");
#endif
	/* Add phenotype latent name */
	Sv_idx2ltt.push_back("PHENOTYPE");

	wsMat Wr = sseMatrix(N_obs+N_pheno, N_ltt+1);
	LOOP(i, N_obs) {
		memcpy(Wr[i], M_Wr0[i], sizeof(wsReal)*N_ltt);
		Wr[i][N_ltt] = W0;
	}
	/* W matrix, covariate part */
	wsUint N_ataObs = N_obs - N_cov;
	LOOP(i, N_cov) {
		Wr[N_obs-i-1][N_ltt-i-1] = W1;
		windex.pop_back();
	}
	LOOP(i, N_cov)
		windex.push_back(vInt());
	/* W matrix, phenotype part */
	vInt tmp;
	LOOP(i, N_pheno) {
		memset(Wr[N_obs+i], 0x00, sizeof(wsReal)*N_ltt);
		if (N_pheno == 1)
			Wr[N_obs+i][N_ltt] = 1;
		else {
			Wr[N_obs+i][N_ltt] = GESCA_Winit;
			tmp.push_back(N_obs+i);
		}
	}
	/* Add index */
	windex.push_back(tmp);

	/* Increase count */
	N_ltt++;			// #ltt += 1 (pheno latent)
	N_obs += N_pheno;	// #obs += #pheno
	//N_ataObs++;			// FIXME: What's this?

	/* Set matrix */
	cStdMatrix M_Wr(N_obs, N_ltt, Wr);

	// #latent - #cov - phenolatent
	//wsUint	N_pwy = N_ltt - N_cov - 1;

	// input Z0 (observed matrix)
	wsUint	N_tot = N_obs + N_ltt;

	vInt bindex;

	// Build W0(weight matrix)
	//         obs          ltt (pwy + cov + phe)
	//     [ 1             | v       ]
	//     [   1           | v       ]
	// obs [     1         |   v     ]
	//     [       1       |   v     ]
	//     [         1     |   v     ]
	//     [           1   |     c   ]
	//     [             1 |       p ]
	cStdMatrix M_W0;
	if (Nv_pathCoords.size() == 0) {
		M_W0.init(M_Wr.row(), M_Wr.col(), M_Wr.get());
		M_W0.setClean(MATDEL_NONE);
	}
	else {
		cIdtMatrix M_Io(N_obs);
		M_W0 = M_Io.rmerge(M_Wr);
	}
	wsMat Ra_W0 = M_W0.get();

	wsUint N_dim = Nv_pathCoords.size() ? N_tot : N_ltt;
	wsUint N_origin = Nv_pathCoords.size() ? N_obs : 0;
	vvInt cindex(N_obs);

	// Generate B0 (beta coefficients vector)
	cStdMatrix	M_B0t(N_dim, N_dim);
	wsMat		B0t = M_B0t.get();
	// 151125 Path coefficients added
	FOREACH(vInt_it, Nv_pathCoords, i) {
		wsUint x = *i >> 16;
		wsUint y = *i & 0xffff;
		B0t[y][x] = GESCA_Winit;
		bindex.push_back(y*N_dim + x);
	}
	for (wsUint i=N_origin ; i<N_dim-1 ; i++) {
		B0t[N_dim-1][i] = GESCA_Winit;
		bindex.push_back((N_dim-1)*N_dim+i);
	}

	// Build C0 (loading matrix)
	cStdMatrix	M_C0t, M_A0t;
	wsMat		C0t		= NULL;
	wsUint		N_parC	= 0;
	if (!B_formative || N_pheno > 1) {
		M_C0t.init(N_obs, N_dim);
		C0t = M_C0t.get();

		wsUint p = 0;
		// Fill the identity part
		if (Nv_pathCoords.size())
			for (; p<N_obs ; p++) C0t[p][p] = W1;

		// Fill the actual W part
		for (; p<N_dim ; p++)
			LOOP(i, N_obs) {
			/* 170115 multiple phenotype */
			if (i < N_ataObs) continue;
			if (Ra_W0[i][p] == GESCA_Winit) {
				cindex[i].push_back(p);
				N_parC++;
				C0t[i][p] = GESCA_Winit;
			}
		}
		M_A0t = M_C0t.bmerge(M_B0t);
	}
	else {
		LOG("No reflective component found, C matrix will be omitted\n");
		M_A0t = M_B0t;
		M_A0t.setClean(MATDEL_NONE);
	}

	wsUint		N_sumWidx = 0;
	FOREACH(vvInt_it, windex, I) N_sumWidx += (wsUint)I->size();

	do {
		/* Print current range of lambda search */
		LOG("Lambda search range [%g", Rv_lambdas[0][0]);
		if (Rv_lambdas.size() > 1) LOGnf(" ~ %g", (*(Rv_lambdas.rbegin()))[0]);
		LOGnf("]");
		if (Rp_lambdaP) {
			LOGnf(" for gene, [%g", Rv_lambdas[0][1]);
			if (Rv_lambdas.size() > 1) LOGnf(" ~ %g", (*(Rv_lambdas.rbegin()))[1]);
			LOGnf("] for pathway");
		}
		LOGnf("\n");

		vector<wsUint*> Nv_tmp;
		wsUint** Na_iters = NULL;
		wsAlloc(Na_iters, wsUint*, N_lambda*N_lambda);
		LOOP(i, N_lambda*N_lambda) wsAlloc(Na_iters[i], wsUint, N_cv);

		if (IS_MULTITHREAD) {
			vector<wsVec> X_tmp;
			vInt tmp;
			vvInt tmp2;
			xGeSCAres tmp3;
			xGescaCvThread X_gc = {
				B_formative,
				M_Yt,
				&M_Xt0,
				M_W0,
				M_B0t,
				M_C0t,
				Rv_lambdas,
				windex,
				bindex,
				cindex,
				N_obs,
				N_ltt,
				N_cv,
				Sv_idx2mnfs,
				Sv_latent,
				Sv_mnfsInLatent,
				Na_samp,
				Na_idxCV,
				Na_iters,
				Ra_devs,
			};
			xAnaThread X_dt = { Cp_IO, &X_gc };
			WORKER().run(thr_optLambda_GESCA, forLambda, &X_dt, NULL, sizeof(int) * 3);
		} else
			sub_optLambda_GESCA(B_formative, M_Yt, &M_Xt0,
				M_W0, M_B0t, M_C0t, windex, bindex, cindex,
				N_obs, N_ltt, Sv_idx2mnfs, Sv_latent, Sv_mnfsInLatent,
				N_cv, Na_samp, Na_idxCV, Rv_lambdas, Rp_lambdaP ? 1 : 0,
				Ra_devs, Na_iters);

		/* Find the minimum across CV */
		wsReal	R_minDev = Ra_devs[0];
		wsUint	N_idxMin = 0;
		wsUint	I = 0;
		char	B_pairwise	= OPT_ENABLED(gxg) && !OPT_ENABLED(gxgall);
		FOREACHDO(vLambda_it, Rv_lambdas, i, I++) {
			if (!B_pairwise || OPT_ENABLED(verbose)) {
#if defined(_DEBUG) || defined(PHARAOHopti)
				Cp_opti->fmt("%g", (*i)[0]);
#endif
				LOG("Lambda [%g", (*i)[0]);
				if (Rp_lambdaP) {
#if defined(_DEBUG) || defined(PHARAOHopti)
					Cp_opti->fmt("	%g", (*i)[1]);
#endif
					LOGnf(",%g", (*i)[1]);
				}
#if defined(_DEBUG) || defined(PHARAOHopti)
				Cp_opti->fmt("	%g\n", Ra_devs[I]);
#endif
				LOGnf("] deviance [%g] iters [", Ra_devs[I]);
				LOOP(j, N_cv)
					LOGnf((j + 1) == N_cv ? "%d" : "%d,", Na_iters[I][j]);
				LOGnf("]\n");
			}

			if (NA(Ra_devs[I])) continue;
			if (R_minDev > Ra_devs[I] || NA(R_minDev)) {
				N_idxMin = I;
				R_minDev = Ra_devs[I];
			}
		}
		LOOP(i, N_lambda*N_lambda) DEALLOC(Na_iters[i]);
		DEALLOC(Na_iters);

		/* Check failure */
		if (NA(R_minDev))
			halt("All optimal candidates were failed!");

		/* Print current status */
		if (Rp_lambdaP)
			LOG("Minimum lambda [%g,%g] with deviance [%g]\n",
			Rv_lambdas[N_idxMin][0], Rv_lambdas[N_idxMin][1], R_minDev);
		else
			LOG("Minimum lambda [%g] with deviance [%g]\n", Rv_lambdas[N_idxMin][0],
			R_minDev);

		/* Set the optimal lambda */
		*Rp_lambdaG = Rv_lambdas[N_idxMin][0];
		if (Rp_lambdaP) *Rp_lambdaP = Rv_lambdas[N_idxMin][1];

		/* Set the new lambdas */
		B_ok = true;
		if (!IS_ASSIGNED(gxglambda)) {
			updateLambda(Rv_intvLambda, Rv_lambdas, N_idxMin, Rp_lambdaP ? 1 : 0);

			/* Met the criterion? */
			FOREACH(vReal_it, Rv_intvLambda, ld) {
				//LOG("Intv %g\n", *ld);
				if (*ld >= 0.01) {
					B_ok = false;
					break;
				}
			}
		}
	} while (!B_ok);

#if defined(_DEBUG) || defined(PHARAOHopti)
	delete Cp_opti;
#endif

	sseFree(Ra_devs);
	LOOP(i, N_cv) DEALLOC(Na_idxCV[i]);
	DEALLOC(Na_idxCV);
	DEALLOC(Na_samp);
}

wsUint do_GeSCA(char B_formative, wsUint b, wsRealCst lambda_b, wsRealCst lambda_w,
	cStdMatrix& M_Xt, cStdMatrix& M_W0_t, cStdMatrix& M_B0_t, cStdMatrix& M_C0_t,
	vvInt& windex, vInt& bindex, vvInt& cindex, xGeSCAres* Xp_res,
	wsUintCst N_obs, wsUintCst N_ltt)
{
	wsUintCst	N_tot		= N_obs + N_ltt;
	char		B_path		= M_W0_t.row() > N_ltt;
	wsUintCst	N_anaSamp	= M_Xt.col();
	wsUintCst	N_origin	= B_path ? N_obs : 0;
	wsUintCst	N_dim		= B_path ? N_tot : N_ltt;

	if (!b) {
		LOG("# of latent variables [%d]\n", N_ltt);
		LOG("# of manifest variables [%d]\n", N_obs);
	}
	cStdMatrix	PsiT;

	// Assign random values to W, B, and C
	cStdMatrix&	Xt	= M_Xt.clone();
	cStdMatrix&	Wt	= M_W0_t.clone();
	cStdMatrix&	Bt	= M_B0_t.clone();
	cStdMatrix&	Ct	= M_C0_t.clone();
//	wsMat		W0t	= M_W0_t.get();
	wsMat		B0t	= M_B0_t.get();
	wsMat		C0t	= M_C0_t.get();
	wsMat		_Wt	= Wt.get();
	wsMat		_Bt	= Bt.get();
	wsMat		_Ct	= Ct.get();

	// W[is.na(W)] <- runif(length(windex))
	LOOP(i, (wsUint)windex.size())
		FOREACH(vInt_it, windex[i], xx) {
#if defined(_DEBUG) || defined(GESCAvalidate)
			if (*xx >=(int) Wt.col()) halt("ERR");
//			if (_Wt[i+N_origin][*xx] != GESCA_Winit) halt("ERR2");
#endif
			_Wt[i+N_origin][*xx] = (wsReal)wsUnifrand();
		}

	// A[is.na(A)] <- runif(length(aindex))
	FOREACH (vInt_it, bindex, i) {
		wsUint x = (wsUint)(*i/N_dim);
		wsUint y = *i%N_dim;
#if defined(_DEBUG) || defined(GESCAvalidate)
		if (x >= Bt.row() || y >= Bt.col()) halt("Dimension error on B");
		if (_Bt[x][y] != GESCA_Winit) halt("Access error on B");
#endif
		_Bt[x][y] = (wsReal)wsUnifrand();
	}
	wsUint x = 0;
	FOREACHDO (vvInt_it, cindex, i, x++) {
		FOREACH (vInt_it, *i, Y) {
			wsUint y = *Y;
			//wsUint x = (wsUint)(*i/N_dim);
			//wsUint y = *i%N_dim;
#if defined(_DEBUG) || defined(GESCAvalidate)
			if (x >= Ct.row() || y >= Ct.col()) halt("Dimension error on C");
			if (_Ct[x][y] != GESCA_Winit) halt("Access error on C");
#endif
			_Ct[x][y] = (wsReal)wsUnifrand();
		}
	}
	//B[is.na(B)] <-runif(length(bindex))

	cStdMatrix At;
	wsUint N_dimA = 0;
	if (Ct.row()) {
		At = Ct.bmerge(Bt);
		N_dimA = N_obs + N_dim;
	} else {
		At = Bt;
		At.setClean(MATDEL_NONE);
		N_dimA = N_dim;
	}

	// data standardization
	// Z    = zscore(bz0)*sqrt(N)/sqrt(N-1)
	cStdMatrix& Zt = Xt.normalize(sqrt(N_anaSamp/(wsReal)(N_anaSamp-1)));
	//Zt.file("Z", 3);
	//Zt.setClean(MATDEL_NONE);

	// Compute matrices
	cIdtMatrix	I(N_obs);
	cStdMatrix	Vt;
	if (Ct.row()) {
		Vt = I.bmerge(Wt);
	} else {
		Vt = Wt;
		Vt.setClean(MATDEL_NONE);
	}
	wsMat		_Vt		= Vt.get();
	// Gamma   <- Z %*% W;
	cStdMatrix	Gamma_t		= Wt * Zt;
	// Psi     <- Z %*% V;
	PsiT	= Vt * Zt;
	wsVec		V0tot	= sseEmptyVec(N_tot);

	wsUint		it		= 0;
	wsReal		f0		= 1000000000000000;
	wsReal		imp		= 1000000000000000;
	wsReal		f		= W0;
	vInt		windex_0;
#define OPTZ
	wsReal		tr_w	= W0;
#ifdef _DEBUG
	if (!b) {
		Wt.file("Winit");
		Bt.file("Binit");
		Ct.file("Cinit");
	}
#endif
	while (it<5000 && imp>1e-3) {
		it++;

#ifdef OPTZ
		cStdMatrix WAt = At * Wt;
		cStdMatrix WAVt = WAt - Vt;
		wsMat Ra_WAVt = WAVt.get();
#endif

		// for (p in 1:N_tot)
		// Step 1: Update W
		LOOP (p, N_dim) {
//			if (N_obs > 5000 && N_ltt > 500)
				pverbose("W update, latent [%d/%d]...\r", p, N_dim);
			// t <- N_obs + p
			wsUint t = Ct.row() ? N_obs + p : p;

			// windex_p <-which(is.na(W0_t[p, ]))
			if (p < N_origin) continue;
			vInt& windex_p = p < N_origin ? windex_0 : windex[p - N_origin];
			if (windex_p.size() == 0) continue;

			// e = zeros(1,N_tot+N_obs)
			// e[1, t] <-1
			// beta <-e - At[,p]
			wsVec _Ap = sseVector(N_dimA);
			LOOP (i, N_dimA)
				_Ap[i] = -At[i][p];

#ifdef OPTZ
			// 170117 optimize (1)
			// WAVij - WipApj + Vij*I(j!=t)
			wsMat curWAVt = sseMatrixP(N_dimA, N_obs, Ra_WAVt);
			sseMaVtV(curWAVt, _Ap, Wt[p], N_dimA, N_obs);
			sseVaV(curWAVt[t], Vt[t], curWAVt[t], N_obs);
			cStdMatrix Delta(N_dimA, N_obs, curWAVt);
#endif
			_Ap[t] += W1;
			cVector beta(_Ap, N_dimA);

#ifdef OPTZ
			// 170117 optimize
			// Store current Wt[p-N_origin]
			wsVec prevW = sseVector((wsUint)windex_p.size());
			wsUint II = 0;
			FOREACHDO(vInt_it, windex_p, i, II++) prevW[II] = Wt[p][*i];
#else
			// 170117 optimize
			// Substitute below to above (1)

			// Delta = W%*%H2%*%	At - V%*%H1;
			wsVec _Wtp = Wt[p];
			_Wt[p] = V0tot;
			cStdMatrix Delta = At * Wt;

			wsVec _Vtt = Vt[t];
			_Vt[t] = V0tot;
			Delta -= Vt;

			_Vt[t] = _Vtt;
			_Wt[p] = _Wtp;
#endif

			// Zp = as.matrix(Z[, windex_p]);
			cStdMatrix Zpt = Zt.subsetPtr(windex_p);

			// theta = try(solve(t(Zp)%*%Zp) %*% (t(Zp)%*%(Z%*%Delta))%*%t(beta) %*% solve(beta%*%t(beta)), T)

			cSymMatrix	ZpZp	= Zpt.Mt();
			ZpZp *= beta.ss();
			if (!cindex[p - N_origin].size()) // If formative
				ZpZp.addDiag(lambda_w);
//			cStdMatrix	ZD		= Delta * Zt;
//			cStdMatrix	ZDZp	= ZD.Mt(Zpt);
//			cVector		ZDZpB	= beta * ZDZp;
//halt("?");
			cSymMatrix&	iZpZp	= ZpZp.inv();
			if (iZpZp.get()[0][0] != iZpZp.get()[0][0])
				halt("X");
			cVector		DB		= beta * Delta;
			cVector		DBZ		= DB * Zt;
			cVector		DBZZp	= Zpt * DBZ;
			cVector		theta	= DBZZp * iZpZp;
//			cVector		theta	= ZDZpB * iZpZp;
			delete& iZpZp;
//			theta /= beta.ss();
//			theta.pinv();

			// zw = Zp%*%theta
			cVector zw = theta * Zpt;
			//				cVector		zw		= Zp * theta;
//			if (B_formative) {
//			}

			// theta = sqrt(N_anaSamp)*theta/norm(zw)
			theta *= sqrt((wsReal)N_anaSamp / zw.ss());
			tr_w += theta.ss();

			wsVec _theta = theta.get();
			wsUint I = 0;
			FOREACHDO (vInt_it, windex_p, i, I++) {
				_Wt[p][*i] = _theta[I];
#ifdef OPTZ
				// 170117 optimize
				// Piecewise update WA
				LOOP (j, N_dimA) Ra_WAVt[j][*i] += (_theta[I] - prevW[I])*At[j][p];
				wsReal Vtti = _Vt[t][*i];
#endif
				_Vt[t][*i] = _theta[I];
#ifdef OPTZ
				// 170117 optimize
				// Piecewise update WAV
				Ra_WAVt[t][*i] += _theta[I] - Vtti;
#endif
			}
#ifdef OPTZ
			// 170117 optimize
			// Free previous W
			sseFree(prevW);
#endif
		}
		Gamma_t = Wt * Zt;

		/* Update C if reflective */
		if (Ct.row()) {
			LOOP (j, N_obs) {
				wsVec	C0j = C0t[j];
				wsVec	Ctj = Ct[j];
				vInt	cindex_j;
				LOOP(i, N_dim) /* FIXME: Check index */
					if (C0j[i] == GESCA_Winit) cindex_j.push_back(i);
				if (cindex_j.size() == 0) continue;

				// Need to compute jth col of Y
				cStdMatrix	Gamma_c	= Gamma_t.subsetPtr(cindex_j);
				cSymMatrix	FjtFj	= Gamma_c.Mt();

				//  mldivide(X%*%t(X), X)
				cVector		Zt_j	= Zt.r2v_ptr(j);
				cVector		V_xxx	= Gamma_c * Zt_j;
				cVector		V_res	= mldivide(FjtFj, V_xxx);

				// C[cindex_j,j] <- res
				wsVec		Ra_res = V_res.get();
				wsUint	I = 0;
				FOREACHDO(vInt_it, cindex_j, i, I++)
					Ctj[*i] = Ra_res[I];
			}
		}

		wsReal tr_b = W0;
		LOOP(p, N_dim) {
			wsVec	B0p = B0t[p];
			wsVec	Btp = Bt[p];
			vInt	bindex_p;
			LOOP(i, N_dim) /* FIXME: Check index */
				if (B0p[i] == GESCA_Winit) bindex_p.push_back(i);
			if (bindex_p.size() == 0) continue;

			// Need to compute jth col of Y
			cStdMatrix	Gamma_b		= Gamma_t.subsetPtr(bindex_p);
			cSymMatrix	FjtFj		= Gamma_b.Mt();
			FjtFj.addDiag(lambda_b);

			//  mldivide(X%*%t(X), X)
			cVector		Gamma_p	= Gamma_t.r2v_ptr(p);
			cVector		V_xxx	= Gamma_b * Gamma_p;
			cVector		V_res	= mldivide(FjtFj, V_xxx);

			// C[cindex_j,j] <- res
			wsVec		Ra_res	= V_res.get();
			wsUint		I		= 0;
			FOREACHDO(vInt_it, bindex_p, i, I++)
				Btp[*i] = Ra_res[I];

			// tr_b           <- tr_b + Bt[p,bindex_p] %*% t(Bt[p,bindex_p])
			tr_b += V_res.ss();
		}
		if (Ct.row())
			At = Ct.bmerge(Bt);
		else {
			At = Bt;
			At.setClean(MATDEL_NONE);
		}

		PsiT = Vt * Xt;
		cStdMatrix AF = At * Gamma_t;
		cStdMatrix dif = PsiT - AF;
		//		f   = tr(t(dif) %*% dif + lambda_b %*% tr_b);
		f = dif.tr2() +tr_b*lambda_b + tr_w*lambda_w;
		imp = f0 - f;
		pverbose("[%d] th iteration, f [%g], imp [%g]\n", it, f, imp);

		f0 = f;
	} /* end(it & imp) */
	sseFree(V0tot);

	wsUint	N_df = N_anaSamp * N_obs;
	wsUint N_sumWindex = 0;
	FOREACH(vvInt_it, windex, i) N_sumWindex += (wsUint)(i->size());
	wsUint	N_npar = N_sumWindex + (wsUint)cindex.size() + (wsUint)bindex.size();

	/*
	*
	* 5. Model fit measures
	*
	*/
	// Fit   <- 1-f/tr(t(Psi) %*% Psi)
	wsReal R_fit = W1 - f / PsiT.tr2();
	// Afit  <-1-(((1-Fit)*DF)/(DF-NPAR))
	wsReal R_Afit = W1 - (((W1 - R_fit)*(wsReal)N_df)/(wsReal)(N_df-N_npar));

	// r2_m  = t(diag(t(CR) %*% t(F) %*% F %*% CR)) / t(diag(t(Z) %*% Z))
	cVector r2_m;
	/* r2_m */
	if (Ct.row()) {
		cStdMatrix CRt = At.subsetPtr(0, N_obs-1);
		cStdMatrix CRtFt = CRt * Gamma_t;
		r2_m = CRtFt.diag2();
		cVector r2_m_denom = Xt.diag2();
		sseVdV(r2_m.get(), r2_m_denom.get(), r2_m.get(), N_obs);
	}

	// r2_s  = t(diag(t(BR) %*% t(F) %*% F %*% BR)) / t(diag(t(F) %*% F))
	cVector r2_s;
	/* r2_s */ {
		cStdMatrix BRt = At.subsetPtr(Ct.row()?N_obs:0, N_dimA-1);
		cStdMatrix BRtFt = BRt * Gamma_t;
		r2_s = BRtFt.diag2();
		cVector r2_s_denom = Gamma_t.diag2();
		sseVdV(r2_s.get(), r2_s_denom.get(), r2_s.get(), N_dim);
	}

	Xp_res->vec_FIT.push_back((wsFloat)R_fit);
	Xp_res->vec_AFIT.push_back((wsFloat)R_Afit);

	wsUint WI = 0, WJ=0, CJ=0, BJ=0;
	wsFvec MatW = NULL;
	wsFvec MatB = NULL;
	wsFvec MatC = NULL;
	wsAlloc(MatW, wsFloat, N_sumWindex);
	wsAlloc(MatB, wsFloat, bindex.size());
	wsCalloc(MatC, wsFloat, cindex.size());
	FOREACHDO (vvInt_it, windex, i, WI++) FOREACH (vInt_it, *i, j) {
		MatW[WJ++] = (wsFloat)_Wt[WI+N_origin][*j];
	}
	x = 0;
	FOREACHDO (vvInt_it, cindex, i, x++) FOREACH (vInt_it, *i, Y) {
		wsUint y = *Y;
		MatC[CJ++] = (wsFloat)_Ct[x][y];
//		MatC[CJ++] = (wsFloat)_Ct[(wsUint)(*i/N_dim)][*i%N_dim];
	}
	FOREACH (vInt_it, bindex, i) {
		MatB[BJ++] = (wsFloat)_Bt[(wsUint)(*i/N_dim)][*i%N_dim];
	}
	wsFvec Ra_r2mDst = NULL;
	wsFvec Ra_r2sDst = NULL;
	if (OPT_ENABLED(verbose)) {
		wsVec Ra_r2mSrc = r2_m.get();
		wsVec Ra_r2sSrc = r2_s.get();
		wsAlloc(Ra_r2mDst, wsFloat, N_obs);
		wsAlloc(Ra_r2sDst, wsFloat, N_tot);
		if (B_path) {
			LOOP (i, N_obs) {
				Ra_r2mDst[i] = (wsFloat)Ra_r2mSrc[i];
				Ra_r2sDst[i] = (wsFloat)Ra_r2sSrc[i];
			}
			for (wsUint i=N_obs ; i<N_tot ; i++)
				Ra_r2sDst[i] = (wsFloat)Ra_r2sSrc[i];
		} else {
			if (Ct.row()) LOOP (i, N_obs)
				Ra_r2mDst[i] = (wsFloat)Ra_r2mSrc[i];
			LOOP (i, N_ltt)
				Ra_r2sDst[i] = (wsFloat)Ra_r2sSrc[i];
		}
	}
	pthread_mutex_lock(&H_permInsertSync);
	Xp_res->MatW.push_back(MatW);
	Xp_res->MatB.push_back(MatB);
	Xp_res->MatC.push_back(MatC);
	if (OPT_ENABLED(verbose)) {
		Xp_res->MatR2M.push_back(Ra_r2mDst);
		Xp_res->MatR2S.push_back(Ra_r2sDst);
	}
	pthread_mutex_unlock(&H_permInsertSync);

	delete& Xt;
	delete& Wt;
	delete& Bt;
	delete& Ct;

	return it;
}

void cPharaohAnalysis::run_GESCA(cPermuter* Cp_perm, char B_formative,
	wsUint N_perm, cStdMatrix& M_Yt, cStdMatrix& M_Xt0,
	cStdMatrix& M_Wr0, vStr& Sv_idx2ltt, vvInt& windex,
	vInt& Nv_pathCoords, wsRealCst R_lambdaG, wsRealCst R_lambdaP,
	vStr& Sv_idx2mnfs)
{
	char		B_pairwise	= OPT_ENABLED(gxg) && !OPT_ENABLED(gxgall);
	wsUintCst	N_pheno		= M_Yt.row();
	wsUint		N_cov		= Cp_IO->sizeCovar();
	wsUint		N_anaSamp	= M_Yt.col();

	/* parameter setting */
	wsUint N_obs = (wsUint)M_Wr0.row();
	wsUint N_ltt = (wsUint)M_Wr0.col();

	/* Add phenotype latent name */
	Sv_idx2ltt.push_back("PHENOTYPE");

	wsMat Wr = sseMatrix(N_obs+N_pheno, N_ltt+1);
	LOOP (i, N_obs) {
		memcpy(Wr[i], M_Wr0[i], sizeof(wsReal)*N_ltt);
		Wr[i][N_ltt] = W0;
	}
	/* W matrix, covariate part */
	wsUint N_ataObs = N_obs - N_cov;
	LOOP (i, N_cov) {
		Wr[N_obs-i-1][N_ltt-i-1] = W1;
		windex.pop_back();
	}
	LOOP (i, N_cov)
		windex.push_back(vInt());
	/* W matrix, phenotype part */
	vInt tmp;
	vPheno& V_phen = Cp_IO->getPhenoInfo();
	LOOP(i, N_pheno) {
		memset(Wr[N_obs+i], 0x00, sizeof(wsReal)*N_ltt);
		if (N_pheno == 1)
			Wr[N_obs+i][N_ltt] = 1;
		else {
			Wr[N_obs+i][N_ltt] = GESCA_Winit;
			tmp.push_back(N_obs+i);
			Sv_idx2mnfs.push_back(V_phen[i].S_name);
		}
	}
	/* Add index */
	windex.push_back(tmp);

	/* Increase count */
	N_ltt++;			// #ltt += 1 (pheno latent)
	N_obs += N_pheno;	// #obs += #pheno
	//N_ataObs++;			// FIXME: What's this?

	/* Set matrix */
	cStdMatrix M_Wr(N_obs, N_ltt, Wr);

	// #latent - #cov - phenolatent
	//wsUint	N_pwy = N_ltt - N_cov - 1;

	// input Z0 (observed matrix)
	wsUint	N_tot = N_obs + N_ltt;

	vInt bindex;

	// Build W0(weight matrix)
	//         obs          ltt (pwy + cov + phe)
	//     [ 1             | v       ]
	//     [   1           | v       ]
	// obs [     1         |   v     ]
	//     [       1       |   v     ]
	//     [         1     |   v     ]
	//     [           1   |     c   ]
	//     [             1 |       p ]
	cStdMatrix M_W0;
	if (Nv_pathCoords.size() == 0) {
		M_W0.init(M_Wr.row(), M_Wr.col(), M_Wr.get());
		M_W0.setClean(MATDEL_NONE);
	} else {
		cIdtMatrix M_Io(N_obs);
		M_W0 = M_Io.rmerge(M_Wr);
	}
	cStdMatrix M_W0_t = M_W0.transpose();
	wsMat W0_t = M_W0_t.get();

	wsUint N_dim = Nv_pathCoords.size() ? N_tot : N_ltt;
	wsUint N_origin = Nv_pathCoords.size() ? N_obs : 0;
	vvInt cindex(N_obs);

	// Generate B0 (beta coefficients vector)
	cStdMatrix	M_B0t(N_dim, N_dim);
	wsMat		B0t = M_B0t.get();
	// 151125 Path coefficients added
	FOREACH (vInt_it, Nv_pathCoords, i) {
		wsUint x = *i >> 16;
		wsUint y = *i & 0xffff;
		B0t[y][x] = GESCA_Winit;
		bindex.push_back(y*N_dim + x);
	}
	for (wsUint i=N_origin ; i<N_dim-1 ; i++) {
		B0t[N_dim-1][i] = GESCA_Winit;
		bindex.push_back((N_dim-1)*N_dim+i);
	}

	// Build C0 (loading matrix)
	cStdMatrix	M_C0t, M_A0t;
	wsMat		C0t		= NULL;
	wsUint		N_parC	= 0;
	if (!B_formative || N_pheno > 1) {
		M_C0t.init(N_obs, N_dim);
		C0t = M_C0t.get();

		wsUint p = 0;
		// Fill the identity part
		if (Nv_pathCoords.size())
			for ( ; p<N_obs ; p++) C0t[p][p] = W1;

		// Fill the actual W part
		for ( ; p<N_dim ; p++)
			LOOP (i, N_obs) {
				/* 170115 multiple phenotype */
				if (i < N_ataObs) continue;
				if (W0_t[p][i] == GESCA_Winit) {
					cindex[i].push_back(p);
					N_parC++;
					C0t[i][p] = GESCA_Winit;
				}
			}
		M_A0t = M_C0t.bmerge(M_B0t);
	} else {
		LOG("No reflective component found, C matrix will be omitted\n");
		M_A0t = M_B0t;
		M_A0t.setClean(MATDEL_NONE);
	}

	wsUint		N_sumWidx = 0;
	FOREACH(vvInt_it, windex, I) N_sumWidx += (wsUint)I->size();
	xGeSCAres	X_res(N_sumWidx, (wsUint)cindex.size(),
		(wsUint)bindex.size(), N_obs, N_ltt);
#ifdef GESCAvalidate
	LOGwarn("#w [%d] #c [%d] #b [%d] #obs [%d] #ltt [%d] #cov [%d]\n",
		N_sumWidx, (wsUint)cindex.size(), (wsUint)bindex.size(),
		N_obs, N_ltt, N_cov);
	LOGwarn("X0 [%d  * %d] A0 [%d * %d] B0 [%d * %d] W0 [%d * %d]\n",
		M_Xt.row(), M_Xt.col(),
		M_A0t.row(), M_A0t.col(), M_B0t.row(), M_B0t.col(),
		M_W0_t.row(), M_W0_t.col());
#endif

	// Add y to both latent and observations
	wsMat __Xt = sseMatrix(N_obs, N_anaSamp);
	LOOP(i, (N_obs-N_pheno))
		memcpy(__Xt[i], M_Xt0.get()[i], sizeof(wsReal)*N_anaSamp);
	LOOP(i, N_pheno)
		memcpy(__Xt[N_obs+i-N_pheno], M_Yt.get()[i], sizeof(wsReal)*N_anaSamp);
	cStdMatrix M_Xt(N_obs, N_anaSamp, __Xt);

#if defined(_DEBUG) || defined(GESCAvalidate)
	M_Xt.file("Xt");
	M_W0_t.file("Wt");
	M_B0t.file("Bt");
	M_C0t.file("Ct");
#endif

	// permutation loop
	if (IS_MULTITHREAD) {
		vLambda Rv_lambdas;
		vReal Rv_lambda;
		Rv_lambda.push_back(R_lambdaG);
		Rv_lambda.push_back(R_lambdaP);
		Rv_lambdas.push_back(Rv_lambda);
		vector<wsVec> tmp;
		vvInt tmp2;
		mDataIdx tmp3;
		vInt tmp4;
		vStr tmp5;
		mvStr tmp6;
		xGscaAnaThread X_gd ={
			1,
			tmp4,
			N_obs,
			&M_Xt, &M_W0_t,
			M_Yt,
			tmp3,
			Sv_idx2mnfs,
			tmp5, tmp6,
			windex,

			Rv_lambdas,
			&N_anaSamp,

			0, NULL, NULL, NULL,

			NULL, NULL,
			tmp, tmp,

			Cp_perm, tmp2,

			&M_B0t, &M_C0t,
			cindex, bindex,
			X_res, N_obs, N_ltt,
		};
		xAnaPermThread X_dt ={ Cp_IO, &X_gd };
		WORKER().run(thr_gesca, forAllPerm_equal, &X_dt, NULL, sizeof(int) * 3);
		FOREACH(vector<wsVec>::iterator, tmp, i) sseFree(*i);
	} else {
		wsUint	N_it = do_GeSCA(B_formative, 0, R_lambdaP, R_lambdaG,
			M_Xt, M_W0_t, M_B0t, M_C0t, windex, bindex, cindex, &X_res,
			N_obs, N_ltt);
		LOG("True estimation of GESCA over with [%d] iterations\n", N_it);
		cStdMatrix	M_permYt(N_pheno, N_anaSamp);
		cStdMatrix	M_permCov(N_cov, N_anaSamp);

		/* Do permutation */ {
			wsUint i = 1;
			for (; ; i++) {
				/* Get permuted phenotype & covariates */
				if (!Cp_perm->fetch(M_permYt, M_permCov)) {
					LOG("[%d] times of permutation has finished\n", i - 1);
					break;
				}
				/* Replace phenotype */
				LOOP(j, N_pheno)
					memcpy(__Xt[N_obs - 1 - j], M_permYt.get()[N_pheno - 1 - j], sizeof(wsReal)*N_anaSamp);
				/* Replace covariates */
				if (!IS_ASSIGNED(permfile))
					LOOP(j, N_cov) memcpy(__Xt[N_obs - N_pheno - N_cov + j], M_permCov.get()[j], sizeof(wsReal)*N_anaSamp);

				/* Do GeSCA */
				N_it = do_GeSCA(B_formative, i, R_lambdaP, R_lambdaG,
					M_Xt, M_W0_t, M_B0t, M_C0t, windex, bindex, cindex,
					&X_res, N_obs, N_ltt);

				if (!B_pairwise || OPT_ENABLED(verbose))
					LOG("Permutation [%d] over after [%d] iterations\n", i + 1, N_it);
				else
					notice("Permutation [%d] over after [%d] iterations\r", i + 1, N_it);
			}
			N_perm = i;
		}
	}

	vStr Sv_hdrB, Sv_hdrW;
	wsUint II = 0;
	FOREACHDO (vvInt_it, windex, i, II++) {
		char S_bname[128];
		FOREACH (vInt_it, *i, j) {
			sprintf(S_bname, "%s->%s", Sv_idx2mnfs[*j].c_str(), Sv_idx2ltt[II].c_str());
			Sv_hdrW.push_back(S_bname);
		}
	}
	FOREACH (vInt_it, bindex, i) {
		char S_bname[128];
		wsUint x = (wsUint)(*i/N_dim);
		wsUint y = *i%N_dim;
//		if (y < N_gene)
//			sprintf(S_bname, "%s->%s", Sv_idx2mnfs[x].c_str(), Sv_idx2mnfs[y].c_str());
//		else if (y >= N_gene)
//			sprintf(S_bname, "%s", Sv_idx2ltt[y-Sv_idx2mnfs.size()].c_str());
		sprintf(S_bname, "%s->%s", Sv_idx2ltt[y].c_str(), Sv_idx2ltt[x].c_str());
		Sv_hdrB.push_back(S_bname);
	}

	exportVector("gesca.fit", X_res.vec_FIT);
	exportVector("gesca.afit", X_res.vec_AFIT);
	exportMatrix("gesca.W", X_res.MatW, N_sumWidx, &Sv_hdrW);
	exportMatrix("gesca.B", X_res.MatB, (wsUint)bindex.size(), &Sv_hdrB);
	exportMatrix("gesca.C", X_res.MatC, (wsUint)cindex.size());
	if (OPT_ENABLED(verbose)) {
		exportMatrix("gesca.R2M", X_res.MatR2M, N_obs);
		exportMatrix("gesca.R2S", X_res.MatR2S, N_tot);
	}
	vReal Rv_Bp(bindex.size());
	vReal Rv_Wp(N_sumWidx);
	cTableExporter	C_exP("gesca.latent.res", "siir",
		"GeSCA latent-level result summary", 0,
		4, "LATENT", "NPERM", "NMANIFEST",
		"P_GESCA");
	/* Get P_latent */ {
		vInt Nv_cntP(bindex.size());
		wsFvec Ra_trueB = X_res.MatB[0];
		vector<wsFvec>::iterator it=X_res.MatB.begin();
		for (it++ ; it != X_res.MatB.end() ; it++)
			LOOP (i, (wsUint)bindex.size())
				Nv_cntP[i] += (*it)[i] > Ra_trueB[i];
		LOOP (i, (wsUint)bindex.size()) {
			Rv_Bp[i] = (Nv_cntP[i] + 1) / ((wsReal)N_perm + 1);
			C_exP.write(4, Sv_hdrB[i].c_str(), N_perm, windex[i].size(), Rv_Bp[i]);
		}
	}

	cTableExporter	C_exW("gesca.manifest.res", "sir",
		"GeSCA manifest-level result summary", 0,
		3, "MANIFEST_TO_LATENT", "NPERM", "P_GESCA");
	/* Construct W*B */
	wsFmat Ra_compW = new wsFvec[X_res.MatW.size()];
	LOOP (i, X_res.MatW.size()) {
		Ra_compW[i] = new wsFloat[N_sumWidx];
		memcpy(Ra_compW[i], X_res.MatW[i], sizeof(wsFloat)*N_sumWidx);
		/* Repeat MatB by each windex times */
		wsUint curS = 0, k = 0;
		FOREACHDO (vvInt_it, windex, it, k++)
			LOOP (j, it->size()) Ra_compW[i][curS++] *= X_res.MatB[i][k];
	}
	/* Get P_manifest */ {
		vInt Nv_cntP(N_sumWidx);
		wsFvec Ra_trueW = Ra_compW[0];
		for (wsUint i=1 ; i< X_res.MatW.size() ; i++)
			LOOP (j, N_sumWidx)
				Nv_cntP[j] += Ra_compW[i][j] > Ra_trueW[j];
		LOOP (i, N_sumWidx) {
			Rv_Wp[i] = (Nv_cntP[i] + 1) / ((wsReal)N_perm + 1);
			C_exP.write(3, Sv_hdrW[i].c_str(), N_perm, Rv_Wp[i]);
		}
	}
	LOOP (i, X_res.MatW.size()) delete[] Ra_compW[i];
	delete[] Ra_compW;

	FOREACH (vector<wsFvec>::iterator, X_res.MatW, i) DEALLOC(*i);
	FOREACH (vector<wsFvec>::iterator, X_res.MatB, i) DEALLOC(*i);
	FOREACH (vector<wsFvec>::iterator, X_res.MatC, i) DEALLOC(*i);
	if (OPT_ENABLED(verbose)) {
		FOREACH(vector<wsFvec>::iterator, X_res.MatR2M, i) DEALLOC(*i);
		FOREACH (vector<wsFvec>::iterator, X_res.MatR2S, i) DEALLOC(*i);
	}
}

#endif

#define NS 1
typedef map<int,vInt> mNetGset;
typedef mNetGset::iterator mNetGset_it;
void _doTest(cIO* Cp_IO, vStr& Xv_genes, mGeneDef& Xm_gd, wsUint N_anaSamp,
	char** Na_data, char* Ba_formative, wsUint N_bst)
{
	char B_bst = N_bst > 0;
	mNetGset Xm_netIdx;
	wsUint I = 0;

	/* Sum out net # variants */
	I = 0;
	FOREACHDO (vStr_it, Xv_genes, i, I++) {
		vInt& Xv_idx = Xm_gd[*i];
		FOREACH (vInt_it, Xv_idx, j)
			Xm_netIdx[*j].push_back(I);
	}
	wsUint N_netVrt	= (wsUint)Xm_netIdx.size();
	wsUint N_ss		= (wsUint)Xv_genes.size()+1;
	vCoord	windex, cindex, bindex;

	/* Construct W
	 * Indicating (first)th variant should be mapped to (seconds) pos */
	cStdMatrix	W(N_netVrt, N_ss);
	wsMat		Wp = W.get();
	I = 0;
	FOREACHDO (mNetGset_it, Xm_netIdx, i, I++) {
		FOREACH (vInt_it, i->second, j) {
			Wp[I][*j] = GESCA_Winit;
			windex.push_back(xCoord(I, *j));
		}
	}

	/* Construct C */
	cStdMatrix C(N_ss, N_netVrt);
	wsMat		Cp = C.get();
	LOOP (p, N_netVrt)
		if (Ba_formative[p] != 1)
			LOOP (q, N_netVrt) {
				Cp[p][q] = Wp[q][p];
				cindex.push_back(xCoord(p, q));
			}

	/* Construct B */
	cStdMatrix	B(N_ss, N_ss);
	wsMat		Bp = B.get();
	LOOP (i, (N_ss-1)) {
		Bp[i][N_ss-1] = GESCA_Winit; // in case of the num of latent pheno = 1
		bindex.push_back(xCoord(i, N_ss-1));
	}
	wsUint	L_w		= (wsUint)windex.size();
	wsUint	L_c		= (wsUint)cindex.size();
	wsUint	L_b		= (wsUint)bindex.size();
	//wsUint	NPAR	= L_w + L_c + L_b;

// 	vec_FIT   = zeros(nbt,1)
// 	vec_AFIT  = zeros(nbt,1)
// 	MatW      = zeros(nbt,length(windex))
// 	MatC      = zeros(nbt,length(cindex))
// 	MatB      = zeros(nbt,length(bindex))
// 	FW        = NULL
	cVector		vec_FIT(N_bst), vec_AFIT(N_bst);
	cStdMatrix	MatW(N_bst, L_w);
	cStdMatrix	MatC(N_bst, L_c);
	cStdMatrix	MatB(N_bst, L_b);

//     MatR2M    = zeros(nbt,J)
//     MatR2S    = zeros(nbt,P)
//     lambda_w  = 0
//     lambda_b  = 0
//     vec_FIT_m = zeros(nbt,1)
//     vec_FIT_s = zeros(nbt,1)
//     vec_GFI   = zeros(nbt,1)
//     vec_SRMR  = zeros(nbt,1)
	wsReal		lambda_w = W0;
	//wsReal		lambda_b	= W0;
	cVector		vec_FIT_m(N_bst), vec_FIT_s(N_bst);
	cVector		vec_GFI(N_bst), vec_SRMR(N_bst);
	cStdMatrix	MatR2M(N_bst, N_netVrt);
	cStdMatrix	MatR2S(N_bst, N_ss);

	/* Construct Z matrix */
	wsMat	Ra_Z	= sseMatrix(N_netVrt, N_anaSamp);
	I = 0;
 	FOREACHDO (mNetGset_it, Xm_netIdx, i, I++) {
		wsUint J = 0;
		FOREACHDO (vInt_it, i->second, j, J++) {
			setGeno(Cp_IO, *j, Ra_Z[J], NULL, 1);
			sseVstdize(Ra_Z[J], N_anaSamp);
			sseVpC(Ra_Z[J], sqrt(N_anaSamp)/sqrt(N_anaSamp-1), Ra_Z[J],
				N_anaSamp);
		}
	}
	cStdMatrix M_Z(N_netVrt, N_anaSamp, Ra_Z);

//	if rank(Z'*Z) == J
//		Z = chol(Z'*Z);
//	end
// 	cStdMatrix M_Z(N_netVrt, N_anaSamp, Ra_Z);
// 	cSymMatrix M_ZtZ = M_Z.tM();
//	wsReal f1 = numeric_limits<wsReal>::max();
	wsMat Ra_bZ0 = Ra_Z;
	if (B_bst != 0) {
		/* For bootstrap */
		wsUint d = 0;
		while (d == 0) {
			/* bootstrap sample */
			wsAlloc(Ra_bZ0, wsVec, N_netVrt);
			LOOP (i, N_anaSamp)
				Ra_bZ0[i] = Ra_Z[wsRand()%N_anaSamp];
//			rb  = rand(N,1)*N
//			rrb = ceil(rb)
//			bz0 = Z0[rrb,]
//			bb  = svd(t(bz0) %*% bz0)$d
			/* Should be transposed due C bz0 = t(bz0) */
			wsSym	Ra_bZZt		= sseMpMt(Ra_bZ0, N_netVrt, N_anaSamp);
			char	B_isFailed	= 0;
			wsVec	Ra_evZZt	= EIGENDECOMPOSITION(Ra_bZZt, N_netVrt, &B_isFailed, NULL);
			d = sseVprod(Ra_evZZt, N_netVrt) > numeric_limits<wsReal>::epsilon();
			//bb[bb<1e-8] = NA
			//d   = abs(prod(bb,na.rm=T)) > eps()
#           //cat("d", d, "\n")
		}

	}
	LOOP (s, NS) {
		wsUint it        = 0;
		//wsReal f0        = 100000000;
		wsReal imp       = 100000000;
		cStdMatrix	A	= C.rmerge(B);
		cStdMatrix	V(N_netVrt, N_netVrt+N_ss);
//		V         = cbind(eye(J), W)
		wsMat		Vp	= V.get();
		LOOP (i, N_netVrt) {
			Vp[i][i] = W1;
			memcpy(Vp[i]+N_netVrt, Wp[i], sizeof(wsReal)*N_ss);
		}
		while (imp>1e-5){
			pverbose("imp %d\r", imp);
//			wsReal tr_w = 0;
			it = it+1;

			LOOP (p, N_ss) {
				wsUint t        = N_netVrt + p;

				vInt windex_p;
				LOOP (q, N_netVrt)
					if (Wp[q][p] == GESCA_Winit)
						windex_p.push_back(q);

// 				m        = zeros(1, T)
// 				m[t]     = 1
//
// 				a        = A[p,]
// 				beta     = m - a
				cVector beta = A.r2v(p);
				beta.neg();
				beta.get()[t] += W1;

// 				WW       = W
// 				WW[,p]   = 0
				cStdMatrix WW;//FIXME: = W.clone();
				WW.setCol(p, W0);
// 				VV       = V
// 				VV[,t]   = 0
// 				Delta    = WW %*% A - VV
				cStdMatrix Delta = WW * A;
				Delta -= V;
				LOOP (i, N_netVrt)
					Delta.get()[i][p] += V.get()[i][p];
//				Zp       = Z[,windex_p,drop=F]
				cStdMatrix Zp = M_Z.ssRow(windex_p);
				// a = (beta %*% t(beta))[1,1] * t(Zp)%*%Zp
				cSymMatrix Zx = Zp.tM(beta.ss());
				// b = lambda_w * eye(length(windex_p))
				Zx.addDiag(lambda_w);

//				theta    = mldivide(a + b, t(Zp))  %*% (Z %*% Delta) %*% t(beta)
				cSymMatrix	ZxZx	= Zx.tM();
				cSymMatrix&	ZxZxi	= ZxZx.inv();
				cStdMatrix	ZZZ		= ZxZxi * Zx;
				delete &ZxZxi;
				cStdMatrix	theta0 = ZZZ.Mt(Zp);

				cVector		DB		= Delta * beta;
				cVector		ZDB		= M_Z * DB;
				cVector		theta	= theta0 * ZDB;

				// zw       = Zp %*% theta
				cVector zw = Zp * theta;
				// theta    = sqrt(N) * theta/norm(zw)
				theta *= sqrt(N_anaSamp) / sqrt(zw.ss());

				wsVecCst theta_p = theta.get();
				wsUint I = 0;
				FOREACHDO (vInt_it, windex_p, i, I++) {
					Wp[*i][p] = theta_p[I];
					Vp[*i][t] = theta_p[I];
				}

				// tr_w     = tr_w + t(theta) %*% theta
//				tr_w = theta.ss();
			}
			cStdMatrix Gamma = M_Z * W;
		}
	}
}

char cPharaohAnalysis::_initGeneDef(mGeneDef& Xm_oriGD, string j, xMaf* Xp_maf, wsUintCst NFILTV)
{
	/* Check is this entry exists already */
	if (Xm_gd.find(j) == Xm_gd.end()) {
		vInt		X_newGD;
		//				mDataIdx_it	X_find	= Xm_gene2storePos.find(*j);
		//	 			/* Can't found this gene */
		//	 			if (X_find == Xm_gene2storePos.end()) continue;
		vInt&		Xv_idxG	= Xm_oriGD[j];

		// 150824 Do not founded for this when --expression
		if (!IS_ASSIGNED(expression)) {

			if (Xv_idxG.size() <= NFILTV) return 0;

			FOREACH(vInt_it, Xv_idxG, k)
				if (Xp_maf[*k].R_allMaf > 0) X_newGD.push_back(*k);

			if (X_newGD.size() <= NFILTV) return 0;
		} else
			X_newGD.push_back(Xv_idxG[0]);

		/* Insert new entry */
		Xm_gd.insert(make_pair(j, X_newGD));
	}

	return 1;
}

void cPharaohAnalysis::_makeGeneSummary(mGeneDef& Xm_gd, const char* Ba_miss,
	map<string, wsVec>& Xm_collapsed, mDataIdx& Xm_nvrtGene, wsUint N_anaSamp)
{
	wsUint	N_sample	= Cp_IO->sizeSample();
	char**	Na_data		= Cp_IO->getGenotype();
	wsVec	Ra_wgt		= Cp_IO->getWeight(Cp_anaPPP);

	/* For each gene */
	FOREACH (mGeneDef_it, Xm_gd, I) {
		vInt&	Xv_indices = I->second;

		wsVec Ra_sumGeno = sseEmptyVec(N_anaSamp);
		if (Xv_indices.size() && (!IS_ASSIGNED(progenesize) ||
			isInRange(OPT_RANGE(progenesize), (wsUint)Xv_indices.size()))) {
			wsUint J = 0;

			/* --expression */
			if (IS_ASSIGNED(expression)) LOOP(i, N_sample) {
				wsFloat** Ra_expr = getIO()->getDosage();

				/* Skip if missing */
				if (Ba_miss[i]) continue;
				wsFloat*	Ra_curExpr	= Ra_expr[i];
				FOREACH(vInt_it, Xv_indices, j) {
					wsFloat R_expr = Ra_curExpr[*j];
					if (isMissingReal(R_expr)) continue;
					Ra_sumGeno[J] += R_expr;
				}
				J++;
			} else LOOP(i, N_sample) {
				/* Skip if missing */
				if (Ba_miss[i]) continue;
				char*	Na_curGeno	= Na_data[i];
				FOREACH(vInt_it, Xv_indices, j) {
					char N_geno = Na_curGeno[*j];
					if (isMissing(N_geno) || !N_geno) continue;
					if (N_geno < 0 || Ra_wgt[*j] < W0)
						halt("SRA");
					Ra_sumGeno[J] += N_geno * Ra_wgt[*j];
				}
				J++;
			}
		}
		char B_sane = 0;
		char B_same = -1;
		wsReal R_same = 0;
		LOOP(i, N_anaSamp) {
			wsReal v = Ra_sumGeno[i];
			if (B_same == -1) {
				R_same = v;
				B_same = 1;
			} else if (B_same && R_same != v)
				B_same = 0;
			if (v)
				B_sane = 1;
		}
		// All-same-value is unwanted
		if (B_same == 1)
			B_sane = 0;

		/* Can use this gene */
		if (B_sane) {
			Xm_collapsed.insert(make_pair(I->first, Ra_sumGeno));
			Xm_nvrtGene[I->first] = (int)Xv_indices.size();
		} else
			sseFree(Ra_sumGeno);
	}
}

void cPharaohAnalysis::_makePathwaySummary(char B_overlap,
	vector<wsVec>& Xv_X, map<string, wsVec>& Xm_collapsed,
	mDataIdx& Xm_mnfs2idx, mDataIdx& Xm_nvrtGene,
	vStr& Sv_idx2mnfs, vStr& Sv_latent, mvStr& Sv_mnfsInLatent,
	vInt& Nv_nvrtInLatent)
{
	eStr			Xm_colGeneUsed;
	cExporter*		Cp_det = NULL;
	if (OPT_ENABLED(verbose))
		 Cp_det = cExporter::summon("gsca.stat");

	/* For each pathway */
	FOREACH (mGeneSet_it, Xm_gs, i) {
		Sv_latent.push_back(i->first);
		vStr&	Xv_curG = i->second;
		/* List of genes in the dataset */
		vStr	Sv_okGenes;

		/* For genes in the pathway */
		FOREACH (vStr_it, Xv_curG, j) {
			/* Is this gene in the dataset ? */
			map<string,wsVec>::iterator X_find = Xm_collapsed.find(*j);

			if (Cp_det) Cp_det->fmt("%s	%s	%s\n", i->first.c_str(),
				j->c_str(), Xm_collapsed.end() == X_find ? "FAILED" : "PASSED");

			/* Skip if there is no such gene or the gene was filtered out
			 * due to its insaneness */
			if (Xm_collapsed.end() == X_find) continue;
			Sv_okGenes.push_back(*j);
		}

		if (Sv_okGenes.size() && (!IS_ASSIGNED(progsetsize) ||
			isInRange(OPT_RANGE(progsetsize), (double)Sv_okGenes.size()))) {

			wsUint N_vrtInPathway = 0;
			if (B_overlap) {
				FOREACH(vStr_it, Sv_okGenes, j) { // If overlap
					Xm_colGeneUsed.insert(*j);
					mDataIdx_it i = Xm_mnfs2idx.find(*j);
					wsUint		I = (wsUint)Xm_mnfs2idx.size();
					if (i == Xm_mnfs2idx.end())
						Xm_mnfs2idx[*j] = I;
					N_vrtInPathway += Xm_nvrtGene[*j];
				}
			} else {
				FOREACH(vStr_it, Sv_okGenes, j) {
					Xm_colGeneUsed.insert(*j); // Check this gene has appeared
					N_vrtInPathway += Xm_nvrtGene[*j];
					Xv_X.push_back(Xm_collapsed[*j]);
					Sv_idx2mnfs.push_back(*j);
				}
			}
			Sv_mnfsInLatent.insert(make_pair(i->first, Sv_okGenes));
			Nv_nvrtInLatent.push_back(N_vrtInPathway);
		}
	}
	wsUint N_sum = 0;
	FOREACH (mvStr_it, Sv_mnfsInLatent, i) N_sum += (wsUint)i->second.size();
	if (!B_overlap && N_sum != (wsUint)Sv_idx2mnfs.size())
		halt("SYSERR: Sv_idx2mnfs building bug");

	/* Make Xv_X if overlap is permitted, using index information */
	if (B_overlap) {
		Xv_X.resize(Xm_mnfs2idx.size());
		Sv_idx2mnfs.resize(Xm_mnfs2idx.size());

		FOREACH (mDataIdx_it, Xm_mnfs2idx, i) {
			wsVec Ra_v = Xm_collapsed[i->first];
#ifdef _DEBUG
			if (!Ra_v) halt("Collapsed gene for [%s] NULL", i->first.c_str());
#endif
			Xv_X[i->second] = Ra_v;
			Sv_idx2mnfs[i->second] = i->first;
		}
	}

	if (Cp_det)
		delete Cp_det;
	/* Report summary */
	if (OPT_ENABLED(gesca))
		LOGnote("[%d] genes mapped onto [%d/%d] pathways\n",
			Xm_colGeneUsed.size(), Sv_mnfsInLatent.size(), Xm_gs.size());
	else {
		wsUint N_var = 0;
		FOREACH(eStr_it, Xm_colGeneUsed, i)
			N_var += (wsUint)(Xm_gd[*i].size());
		LOGnote("[%d] genes with [%d] variants mapped onto [%d/%d] pathways\n",
			Xm_colGeneUsed.size(), N_var, Sv_mnfsInLatent.size(), Xm_gs.size());
		lverbose("# of manifest variable [%d]\n", N_sum);
	}

	if (OPT_ENABLED(verbose)) {
//		exportVector("gsca.ngene", Nv_genes, N_anaSamp);
		exportVector("gsca.gname", Sv_idx2mnfs, (wsUint)Sv_idx2mnfs.size());
		exportVector("gsca.pname", Sv_latent, (wsUint)Xm_gs.size());
	}
}

vvInt _buildMatrixW(char B_overlap, wsUintCst N_pheno, wsUintCst N_pathway,
	cStdMatrix& M_W0, mvStr& Smv_mnfsInLatent, vStr& Sv_latent, mDataIdx& Xm_mnfs2idx,
	wsUint* Np_sumWidx)
{
	int			kk			= -1;
	wsUintCst	N_latent	= M_W0.col(); // #col = N_latent
	wsMat		Ra_W0		= M_W0.get();
	vvInt		windex(N_latent);

#ifdef GESCAvalidate
	LOG("Latents: ");
#endif
	if (B_overlap) LOOP(j, N_latent) {
#ifdef GESCAvalidate
		LOGnf("%s ", Sv_latent[j].c_str());
#endif
		vStr& Sv_curMnfs = Smv_mnfsInLatent[Sv_latent[j]];
		FOREACH(vStr_it, Sv_curMnfs, k) {
			wsUint idx = Xm_mnfs2idx[*k];
			Ra_W0[idx][j] = GESCA_Winit;
			windex[j].push_back(idx);
		}
	}
	else LOOP(j, N_latent) {
#ifdef GESCAvalidate
		LOGnf("%s ", Sv_latent[j].c_str());
#endif
		if (j < N_pathway || N_pheno == 1) {
			wsUint Nj = (wsUint)Smv_mnfsInLatent[Sv_latent[j]].size();
			wsUint k = kk + 1;
			kk += Nj;
			for (int K = k; K <= kk; K++) {
				Ra_W0[K][j] = GESCA_Winit;
				windex[j].push_back(K);
			}
		} else {
			kk++;
			LOOP(k, N_pheno) {
				Ra_W0[kk][j + k] = GESCA_Winit;
				windex[j + k].push_back(kk);
			}
			j += N_pheno - 1;
		}
	}

#ifdef GESCAvalidate
	LOGnf("\n");
#endif
	wsUint N_sumWidx = 0;
	LOOP (i, N_latent) N_sumWidx += (wsUint)windex[i].size();
	if (N_sumWidx == 0)
		windex.clear();
	if (Np_sumWidx) *Np_sumWidx = N_sumWidx;

	return windex;
}

void _buildLambda(wsReal& R_lambdaG, wsReal& R_lambdaP)
{
	if (!IS_ASSIGNED(prolambda)) return;
	wsUint	N_lambda = 0;
	wsVec	Ra_lambda = loadRealValues(OPT_STRING(prolambda), &N_lambda);
	switch (N_lambda) {
	case 1:
		R_lambdaG = R_lambdaP = Ra_lambda[0];
		LOGnote("Lambda is fixed to [%g]\n", R_lambdaG);
		break;
	case 2:
		R_lambdaG = Ra_lambda[0];
		R_lambdaP = Ra_lambda[1];
		LOGnote("Lambda is fixed to [%g] for gene and [%g] for pathway\n",
			R_lambdaG, R_lambdaP);
		break;
	default:
		halt("Number of given lambda [%d] is invalid!", N_lambda);
		break;
	}
	sseFree(Ra_lambda);
}

typedef map<string, wsVec> mColGeno;
typedef mColGeno::iterator mColGeno_it;

void cPharaohAnalysis::run()
{
	pthread_mutex_init(H_permInsertSync);

	/* Set the default value for promaxiter */
	if (!IS_ASSIGNED(promaxiter)) {
		OPTION().assign("promaxiter");
		OPTION().FORCE_OPT_NUMBER(promaxiter);
	}
	/* Set the default value for prothr */
	if (!IS_ASSIGNED(prothr)) {
		OPTION().assign("prothr");
		OPTION().FORCE_OPT_REAL(prothr);
	}

	/* Check missingness */
	char	B_formative	= 1;
	wsUint	N_pheno		= Cp_IO->sizePheno();
	wsUint	N_anaSamp	= 0;
	char*	Ba_miss	= Cp_IO->getPheCovMissing(&N_anaSamp);
	if (N_anaSamp == 0) {
		LOGwarn("No sample left after missingness filtering, PHARAOH cannot be performed!\n");
		return;
	}

	if (OPT_ENABLED(gesca)) {
		// Check model type
		if (!IS_ASSIGNED(modeltype)) {
			OPTION().assign("modeltype", OPTION().getDefVal("modeltype"), 1);
			OPTION().FORCE_OPT_STRING(modeltype);
		}
		if (!stricmp(OPT_STRING(modeltype), "formative")) {
			if (N_pheno > 1)
				halt("--modeltype should be 'multivariate' for multiple phenotypes");
		} else if (!stricmp(OPT_STRING(modeltype), "reflective")) {
			if (N_pheno > 1)
				halt("--modeltype should be 'multivariate' for multiple phenotypes");
			B_formative = 0;
		} else if (!stricmp(OPT_STRING(modeltype), "multivariate")) {
			if (N_pheno < 2)
				halt("--modeltype multivariate cannot be used for single phenotype");
		} else
			halt("--modeltype must be either of formative / reflective for single phenotype or multivariate for multiple phenotypes!");
	}

	if (OPT_ENABLED(gesca) && OPT_ENABLED(gxg) && !OPT_ENABLED(gxgall)) {
		if (IS_ASSIGNED(gxglist)) {
			/* Selective pairwise */
			cStrFile C_gxglist(OPT_STRING(gxglist), "File contains gene pairs to test GxG");
			mDataIdx Xm_gene2idx;
			wsUint X = 0;
			FOREACHDO(mGeneDef_it, Xm_gd, I, X++) Xm_gene2idx.insert(make_pair(I->first, X));

			char *S_buf = NULL;
			wsAlloc(S_buf, char, 65536);

			wsUint L = 0;
			for (char *a=NULL, *b=NULL ; C_gxglist.gets(S_buf, 2048) ; L++) {
				getString(&S_buf, &a); // S_buf == gene1
				getString(&a, &b); // a = gene2
				mDataIdx_it X_find1 = Xm_gene2idx.find(S_buf);
				mDataIdx_it X_find2 = Xm_gene2idx.find(a);
				if (X_find1 == Xm_gene2idx.end() || X_find2 == Xm_gene2idx.end()) continue;
				vInt tmp;
				if (X_find1->second > X_find2->second) {
					tmp.push_back(X_find2->second);
					tmp.push_back(X_find1->second);
				} else {
					tmp.push_back(X_find1->second);
					tmp.push_back(X_find2->second);
				}
				run2(N_pheno, N_anaSamp, Ba_miss, B_formative, tmp);
			}

			DEALLOC(S_buf);
		} else {
			/* Exhaustive pairwise */
			wsUint X = 0;
			FOREACHDO(mGeneDef_it, Xm_gd, I, X++) {
				wsUint Y = 0;
				FOREACHDO(mGeneDef_it, Xm_gd, J, Y++) {
					if (Y<=X) continue;
					vInt tmp;
					tmp.push_back(X);
					tmp.push_back(Y);
					run2(N_pheno, N_anaSamp, Ba_miss, B_formative, tmp);
				}
			}
		}
	} else
		run2(N_pheno, N_anaSamp, Ba_miss, B_formative, vInt());
	DEALLOC(Ba_miss);
}

void cPharaohAnalysis::run2(wsUintCst N_pheno, wsUint N_anaSamp, char* Ba_miss,
	char B_formative, vInt Nv_comb)
{
	char		B_overlap	= OPT_ENABLED(gesca) ? 1 : 0;

	mvStr		Smv_mnfsInLatent;	///< [#latent * #mnfs_i]	i(th) latent variable's manifest variables' names
	vStr		Sv_latent;			///< [#latent],				i(th) latent variable's name
	vStr		Sv_idx2mnfs;		///< [#mnfs],				i(th) manifest variable's name
	vInt		Nv_variants;
	mDataIdx	Xm_mnfs2idx;

	/* Make burden data */
	mColGeno	Xm_collapsed;
	mDataIdx	Xm_nvrtGene;

	cStdMatrix	M_Xt;
	wsUint		N_gene		= 0;
	wsUint		N_mnfs		= 0;
	wsUint		N_pathway	= 0;
	wsMat		Ra_Xt		= NULL;					///< [Ng*Ns]
	vInt		Nv_dist;

	char*		Sp_oldPrefix= NULL;
	string		S_oldOut	= OPT_STRING(out);

	/* Initiate necessary values */
	wsUint		N_sample	= Cp_IO->sizeSample();
	vPheno&		Xv_phe		= Cp_IO->getPhenoInfo();
	wsMat		Ra_cov		= Cp_IO->getCovariates();
	wsUint		N_cov		= getIO()->sizeCovar();
	vCovar&		Xv_cov		= getIO()->getCovInfo();

	/* Set the lambda */
	wsReal		R_lambdaG = W0;
	wsReal		R_lambdaP = W0;
	_buildLambda(R_lambdaG, R_lambdaP);

	/* Get available y and covs */
	vector<cVector> Xv_cov2rem(N_cov);

	if (OPT_ENABLED(gesca) && OPT_ENABLED(gxg)) {
		wsUint			B_all	= OPT_ENABLED(gxgall) ? 1 : 0;
		vVariant&		Xv_vrt	= getIO()->getVariant();
		char**			Na_data	= getIO()->getGenotype();
		vector<wsVec>	Xv_X;
		if (Xv_vrt.size() > 300 && B_all)
			halt("GeSCA with GxG can be done with equal/less than 300 variants");
		// LOGnote("GeSCA with GxG for [%d] variants will be done pair-wise manner, [%d] pairs", Xv_vrt.size(),
		// (wsReal)(Xv_vrt.size()*(Xv_vrt.size() - 1)) / W2);

		/* Make burden data */
		/* For each gene */
		int		X = 0, Y = 0;
		vStr	Sv_gene;
		FOREACHDO (mGeneDef_it, Xm_gd, I, X++) {
			/* In case of 'selected subset', see only all of them */
			if (!B_all) {
				if (Y == (int)Nv_comb.size()) break;
				if (Nv_comb[Y] != X) continue;
				Y++;
			}
			Sv_gene.push_back(I->first);
			vInt&	Xv_indices = I->second;
			vStr	Sv_curVrt;

			if (Xv_indices.size() && (!IS_ASSIGNED(progenesize) ||
				isInRange(OPT_RANGE(progenesize), (wsUint)Xv_indices.size()))) {

				Sv_latent.push_back(I->first);

				FOREACH (vInt_it, Xv_indices, K) {
					xVariant&	X_vrt	= Xv_vrt[*K];
					string		S_name	= X_vrt.name;

					// Insert if not inserted
					if (Xm_mnfs2idx.find(S_name) == Xm_mnfs2idx.end()) {
						// Generate vector
						wsVec		Ra_geno	= sseVector(N_anaSamp);
						for (wsUint j=0, J=0 ; j < N_sample ; j++) {
							if (Ba_miss[j]) continue;
							Ra_geno[J++] = (wsReal)Na_data[j][*K];
						}
						// Insert index
						wsUint N_idxMnfs = (wsUint)Xm_mnfs2idx.size();
						Xm_mnfs2idx.insert(make_pair(S_name, N_idxMnfs));
						Nv_variants.push_back(1);
						Xm_nvrtGene[S_name]++;
						Xv_X.push_back(Ra_geno);

						// Although it has already inserted, name should be inserted again
						Sv_idx2mnfs.push_back(S_name);
						Sv_curVrt.push_back(S_name);
					}
				}

				//
				Smv_mnfsInLatent.insert(make_pair(I->first, Sv_curVrt));
			}
		}
		Sp_oldPrefix = OPT_STRING(out);
		string S_out = "";
		FOREACH(vStr_it, Sv_gene, it) {
			S_out += S_out.length() ? "_"+*it : *it;
		}
		S_oldOut += ".";
		S_oldOut += S_out;
		OPT_STRING(out) = (char *)S_oldOut.c_str();

		/* For each gene*gene part */
		X = Y = 0;
		FOREACHDO (mGeneDef_it, Xm_gd, I, X++) {
			vInt&	Xv_indices1 = I->second;
			/* In case of 'selected subset', see only all of them */
			if (!B_all) {
				if (Y == (int)Nv_comb.size()) break;
				if (Nv_comb[Y] != X) continue;
				Y++;
			}

			if (Xv_indices1.size() && (!IS_ASSIGNED(progenesize) ||
				isInRange(OPT_RANGE(progenesize), (wsUint)Xv_indices1.size()))) {

				char B_start = 0;
				int XX = 0, YY = 0;
				FOREACHDO (mGeneDef_it, Xm_gd, II, XX++) {
					/* In case of 'selected subset', see only all of them */
					if (!B_all) {
						if (YY == (int)Nv_comb.size()) break;
						if (Nv_comb[YY] != XX) continue;
						YY++;
					}
					if (!B_start) {
						if (I->first == II->first) B_start = 1;
						continue;
					}
					vInt&	Xv_indices2 = II->second;
					vStr	Sv_curVrt;

					if (Xv_indices2.size() && (!IS_ASSIGNED(progenesize) ||
						isInRange(OPT_RANGE(progenesize), (wsUint)Xv_indices2.size()))) {

						// GxG term name
						string S_nameLatent = I->first + "*" + II->first;
						Sv_latent.push_back(S_nameLatent);

						FOREACH (vInt_it, Xv_indices1, K) FOREACH (vInt_it, Xv_indices2, KK) {
							xVariant&	X_vrt1	= Xv_vrt[*K];
							xVariant&	X_vrt2	= Xv_vrt[*KK];

							// G*G variant name
							string		S_name1	= X_vrt1.name;
							S_name1 += "*";
							S_name1 += X_vrt2.name;
							string		S_name2	= X_vrt2.name;
							S_name2 += "*";
							S_name2 += X_vrt1.name;

							// Insert if not inserted
							if (Xm_mnfs2idx.find(S_name1) == Xm_mnfs2idx.end() &&
								Xm_mnfs2idx.find(S_name2) == Xm_mnfs2idx.end()) {
								// Generate vector
								wsReal		Xv = 0;
								wsVec		Ra_geno	= sseVector(N_anaSamp);
								for (wsUint j=0, J=0 ; j < N_sample ; j++) {
									if (Ba_miss[j]) continue;
									wsReal v = (wsReal)Na_data[j][*K] * Na_data[j][*KK];
									if (!J) Xv = v;
									else if (Xv != v) Xv = WISARD_NAN;
									Ra_geno[J++] = (wsReal)Na_data[j][*K] * Na_data[j][*KK];
								}
								if (Xv != Xv) {
									// Insert index
									wsUint N_idxMnfs = (wsUint)Xm_mnfs2idx.size();
									Xm_mnfs2idx.insert(make_pair(S_name1, N_idxMnfs));
									Nv_variants.push_back(1);
									Xm_nvrtGene[S_name1]++;
									Xv_X.push_back(Ra_geno);

									// Although it has already inserted, name should be inserted again
									Sv_idx2mnfs.push_back(S_name1);
									Sv_curVrt.push_back(S_name1);
								}
							}
						}

						//
						Smv_mnfsInLatent.insert(make_pair(S_nameLatent, Sv_curVrt));
					}
				}
			}
		}

		// Now we have gene-level genotypes, generate gxg latents
		N_pathway = (wsUint)Sv_latent.size();

		/* Now we can get #gene and #mnfs */
		N_gene	= (wsUint)Xv_X.size();
		N_mnfs	= N_gene + N_cov;

		/* Convert design matrix into "real" matrix */
		wsCalloc(Ra_Xt, wsVec, N_mnfs);
		wsUint I = 0;
		FOREACHDO(vector<wsVec>::iterator, Xv_X, i, I++)
			Ra_Xt[I] = *i;
	} else {
		ASSERT_OPTS_AND(geneset, gsetconsec);

		/* Check distribution */
		LOOP(i, N_pheno)
			Nv_dist.push_back(Cp_IO->isContinuous() ? 1 : 2);
	//	wsUint		N_dist		= (!) + 1;
		/* Make burden data */
		_makeGeneSummary(Xm_gd, Ba_miss, Xm_collapsed, Xm_nvrtGene, N_anaSamp);

#if TOOLSET_TYPE != TOOLSET_ONETOOL
		if (OPT_ENABLED(makecolgeno)) {
			cExporter* Cp_cg = cExporter::summon("collapsed.geno");
			/* Export header */
			Cp_cg->put("FID	IID	PAT	MAT	SEX	PHENOTYPE");
			FOREACH (mColGeno_it, Xm_collapsed, i)
				Cp_cg->fmt("	%s", i->first.c_str());
			Cp_cg->put("\n");

			/* Export collapsed genotypes */
			vSampPtr&	Xv_samp		= getIO()->getSample();
			wsMat		Ra_pheno	= getIO()->getPhenos();
			wsUint		I			= 0;
			FOREACHDO (vSampPtr_it, Xv_samp, i, I++) {
				if (Ba_miss[I]) continue;
				Cp_cg->fmt("%s	%s	%s	%s	%d	%g", (*i)->S_FID.c_str(),
					(*i)->S_IID.c_str(), (*i)->Xp_pat?(*i)->Xp_pat->S_IID.c_str():"0",
					(*i)->Xp_mat?(*i)->Xp_mat->S_IID.c_str():"0",
					(*i)->N_sex, Ra_pheno[0][I]);
	//			for (wsUint j=0 ; j<N_pheno ; j++)
	//				Cp_cg->fmt("	%g", Ra_pheno[j][I]);
				FOREACH (mColGeno_it, Xm_collapsed, j)
					Cp_cg->fmt("	%g", j->second[I]);
				Cp_cg->put("\n");
			}
			delete Cp_cg;
		}
#endif
		/* Make pathway-wide thing */ {
			vector<wsVec>	Xv_X;				///< [#gene * #latent]		design matrix consists of collapsed genes
			_makePathwaySummary(B_overlap, Xv_X, Xm_collapsed, Xm_mnfs2idx,
				Xm_nvrtGene, Sv_idx2mnfs, Sv_latent, Smv_mnfsInLatent,
				Nv_variants);
			N_pathway = (wsUint)Sv_latent.size();

			/* Now we can get #gene and #mnfs */
			N_gene	= (wsUint)Xv_X.size();
			N_mnfs	= N_gene+N_cov;

			/* Convert design matrix into "real" matrix */
			wsCalloc(Ra_Xt, wsVec, N_mnfs);
			wsUint I = 0;
			FOREACHDO (vector<wsVec>::iterator, Xv_X, i, I++)
				Ra_Xt[I] = *i;
		}
	}

	/* Add Nv_genes and Ra_Xt by # of covariates */
	LOOP (i, N_cov) {
		string	S_covName = Xv_cov[i].Sp_varName;

		// If overlap, insert covariate name to mnfs2index
		if (B_overlap) Xm_mnfs2idx[S_covName] = (wsUint)Sv_idx2mnfs.size();
		// also create reverse index
		Sv_idx2mnfs.push_back(S_covName);

		/* 170115 If multiple phenotype */
		if (N_pheno > 1 && OPT_ENABLED(pharaoh)) {
			vStr	Sv_covName;
			Sv_covName.push_back(S_covName);
			string	S_latentName;
			LOOP (j, N_pheno) {
				S_latentName = Xv_phe[j].S_name + "->" + Xv_cov[i].Sp_varName;

				Smv_mnfsInLatent.insert(make_pair(S_latentName, Sv_covName));

				/* Add latent */
				Sv_latent.push_back(S_latentName);
			}
		} else {
			vStr	Sv_covName;
			Sv_covName.push_back(S_covName);
			Smv_mnfsInLatent.insert(make_pair(S_covName, Sv_covName));

			/* And latent */
			Sv_latent.push_back(S_covName);
		}
	}

	/* Adjust matrix if samples are adjusted */
	if (N_sample == N_anaSamp) LOOP(i, N_cov) {
		Ra_Xt[N_gene+i] = Ra_cov[i];
	} else LOOP(i, N_cov) {
		wsVec	Ra_ncov		= sseVector(N_anaSamp);
		wsVec	Ra_curCov	= Ra_cov[i];
		wsUint	J			= 0;
		LOOP(j, N_sample) {
			if (Ba_miss[j]) continue;
			Ra_ncov[J++] = Ra_curCov[j];
		}
		Ra_Xt[N_gene+i] = Ra_ncov;
		Xv_cov2rem[i].init(N_anaSamp, Ra_ncov, NULL, 0, 0);
	}

	/* Now it is complete design matrix */
	M_Xt.init(N_mnfs, N_anaSamp, Ra_Xt);
	M_Xt.setClean(MATDEL_ROW);
	if (OPT_ENABLED(verbose)) {
		cStdMatrix M_X = M_Xt.transpose();
		exportMatrix("gsca.mat", M_X.get(), N_anaSamp, N_gene, &Sv_idx2mnfs, NULL);
	}
#ifdef PHARAOHvalidate
	M_Xt.file("gsca.X");
#endif

	cStdMatrix M_cov;
	if (N_cov) M_cov = M_Xt.subsetPtr(N_gene, N_gene + N_cov - 1);

	/* phenoRel */
	vvInt M_phenoRel;
	LOOP(i, N_pheno) {
		vInt V_ins;
		if (OPT_ENABLED(pharaoh)) {
			LOOP(j, N_pathway) V_ins.push_back(1);
			LOOP(j, N_cov) LOOP(k, N_pheno) V_ins.push_back(i == k ? 1 : 0);
		} else {
			// FIXME : Impl required
		}
		M_phenoRel.push_back(V_ins);
	}

	/* Load path data */
	vInt		Nv_pathCoords;
	if (IS_ASSIGNED(ggpath)) {
		cStrFile	C_path(OPT_STRING(ggpath), "Gene-gene path file");
		/* idx->str to str->idx */
		mDataIdx Xm_latent2idx;
		wsUint j = 0;
		FOREACHDO (vStr_it, Sv_latent, i, j++) Xm_latent2idx[*i] = j;

		char*	S_bufGG = new char[512];
		wsUint	N_recv = 0;
		while (C_path.gets(S_bufGG, 512)) {
			char *Sp_g1 = NULL, *Sp_g2 = NULL, *Sp_e;
			getString(&S_bufGG, &Sp_g1);
			if (!Sp_g1) continue;
			getString(&Sp_g1, &Sp_g2);
			if (!Sp_g2) continue;
			getString(&Sp_g2, &Sp_e);
			// S_bufGG = pathway name
			mDataIdx_it X_find = Xm_latent2idx.find(S_bufGG);
			if (X_find == Xm_latent2idx.end()) continue;

			// Sp_g1 = gene 1
			mDataIdx_it X_findg1 = Xm_mnfs2idx.find(Sp_g1);
			if (X_findg1 == Xm_mnfs2idx.end()) continue;

			// Sp_g2 = gene 2
			mDataIdx_it X_findg2 = Xm_mnfs2idx.find(Sp_g2);
			if (X_findg2 == Xm_mnfs2idx.end()) continue;

			// All available
			Nv_pathCoords.push_back(((X_findg1->second)<<16) | X_findg2->second);

			N_recv++;
		}
		LOG("Gene-gene path [%s] assigned with [%d] assignable paths\n",
			OPT_STRING(ggpath), N_recv);
	}
	// Reconstruct Sv_latent
	vStr Sv_latentNew;
	FOREACH(vStr_it, Sv_latent, i)
		if (Smv_mnfsInLatent[*i].size() != 0)
			Sv_latentNew.push_back(*i);
	Sv_latent.clear();
	FOREACH(vStr_it, Sv_latentNew, i)
		Sv_latent.push_back(*i);
	for (auto it=Smv_mnfsInLatent.cbegin() ; it != Smv_mnfsInLatent.cend() ;) {
		if (!it->second.size()) {
			LOGnote("Relation [%s] has removed due to missingness\n",
				it->first.c_str());
			N_pathway--;
			it = Smv_mnfsInLatent.erase(it);
		}  else
			it++;
	}
	wsUint		N_latent	= (wsUint)Smv_mnfsInLatent.size();

	// W0 = zeros(sum_nvar,ndset)
	cStdMatrix	M_W0(N_mnfs, N_latent);
	wsUint		N_sumWidx	= 0;
	vvInt		windex		= _buildMatrixW(B_overlap, N_pheno, N_pathway, M_W0,
		Smv_mnfsInLatent, Sv_latent, Xm_mnfs2idx, &N_sumWidx);
	if (!windex.size()) {
		LOGwarn("PHARAOH cannot be proceed, no available gene left");
		return;
	}

	/* Get available y and covs */
	wsMat		Ra_Yt	= NULL;
	if (N_sample == N_anaSamp)
		Ra_Yt = sseMatrixP(N_pheno, N_anaSamp, Cp_IO->getPhenos());
	else {
		Ra_Yt = sseMatrix(N_pheno, N_anaSamp);
		wsMat	Ra_oPhe = Cp_IO->getPhenos();
		LOOP (j, N_pheno) {
			wsUint I = 0;
			LOOP (i, N_sample) {
				if (Ba_miss[i]) continue;
				Ra_Yt[j][I++] = Ra_oPhe[j][i];
			}
		}
	}
	cStdMatrix M_Yt(N_pheno, N_anaSamp, Ra_Yt);

	/* Initiate permuter */
	cPermuter*	Cp_perm		= cPermuter::summon(getIO(), Ra_Yt, N_pheno, N_anaSamp, M_cov.get(), N_cov, Ba_miss);
	if (Cp_perm == NULL) halt_fmt(WISARD_SYST_NULL_DATA, "Permutation class");

	if (OPT_ENABLED(gesca)) {
		// Check model type
		if (!IS_ASSIGNED(modeltype)) {
			OPTION().assign("modeltype", OPTION().getDefVal("modeltype"), 1);
			OPTION().FORCE_OPT_STRING(modeltype);
		}
		char B_formative = 1;
		if (!stricmp(OPT_STRING(modeltype), "formative")) {
			if (N_pheno > 1)
				halt("--modeltype should be 'multivariate' for multiple phenotypes");
		} else if (!stricmp(OPT_STRING(modeltype), "reflective")) {
			if (N_pheno > 1)
				halt("--modeltype should be 'multivariate' for multiple phenotypes");
			B_formative = 0;
		} else if (!stricmp(OPT_STRING(modeltype), "multivariate")) {
			if (N_pheno < 2)
				halt("--modeltype multivariate cannot be used for single phenotype");
		} else
			halt("--modeltype must be either of formative / reflective for single phenotype or multivariate for multiple phenotypes!");

		// Run GeSCA
		LOG("Run GeSCA with [%s] %s gene-gene path\n",
			N_pheno>1 ? "multiple phenotypes" : (B_formative ? "formative" : "reflective"), Nv_pathCoords.size() ? "with" : "without");
		/* Set the lambda */
		if (!IS_ASSIGNED(prolambda)) {
			vvInt windex_(windex.size());
			copy(windex.begin(), windex.end(), windex_.begin());
			LOGnote("Lambda will be chosen optimally\n");
			if (OPT_ENABLED(prosingle)) {
				optLambda_GESCA(getIO(), Cp_perm, B_formative,
					OPT_NUMBER(nperm), M_Yt, M_Xt, M_W0,
					Sv_latent, windex_, Sv_idx2mnfs, Sv_latent,
					Smv_mnfsInLatent, Nv_pathCoords,
					&R_lambdaG, NULL);
				R_lambdaP = R_lambdaG;
			} else
				optLambda_GESCA(getIO(), Cp_perm, B_formative,
					OPT_NUMBER(nperm), M_Yt, M_Xt, M_W0,
					Sv_latent, windex_, Sv_idx2mnfs, Sv_latent,
					Smv_mnfsInLatent, Nv_pathCoords,
					&R_lambdaG, &R_lambdaP);
//			R_lambdaG *= (wsReal)OPT_NUMBER(cv);
//			R_lambdaP *= (wsReal)OPT_NUMBER(cv);
			if (R_lambdaG == R_lambdaP)
				LOGnote("Optimal lambda is chosen to [%g]\n", R_lambdaG);
			else
				LOGnote("Optimal lambda is chosen to [%g,%g]\n", R_lambdaG, R_lambdaP);
		}

		/* Stop if --proopt */
		if (OPT_ENABLED(proopt))
			return;

		run_GESCA(Cp_perm, B_formative, OPT_NUMBER(nperm), M_Yt, M_Xt,
			M_W0, Sv_latent, windex, Nv_pathCoords, R_lambdaG, R_lambdaP, Sv_idx2mnfs);
	} else {
		/* Set the lambda */
		if (!IS_ASSIGNED(prolambda)) {
			LOGnote("Lambda will be chosen optimally\n");
			if (OPT_ENABLED(prosingle)) {
				optLambda_PHARAOH(getIO(), B_overlap, M_Xt, M_Yt, M_phenoRel, Nv_dist, N_mnfs,
					Xm_mnfs2idx, Sv_idx2mnfs,
					Sv_latent, Smv_mnfsInLatent,
					M_W0, windex, &R_lambdaG, NULL);
				R_lambdaP = R_lambdaG;
			} else
				optLambda_PHARAOH(getIO(), B_overlap, M_Xt, M_Yt, M_phenoRel, Nv_dist, N_mnfs,
					Xm_mnfs2idx, Sv_idx2mnfs,
					Sv_latent, Smv_mnfsInLatent,
					M_W0, windex, &R_lambdaG, &R_lambdaP);
			R_lambdaG *= (wsReal)OPT_NUMBER(cv);
			R_lambdaP *= (wsReal)OPT_NUMBER(cv);
			if (R_lambdaG == R_lambdaP)
				LOGnote("Optimal lambda is chosen to [%g]\n", R_lambdaG);
			else
				LOGnote("Optimal lambda is chosen to [%g,%g]\n", R_lambdaG, R_lambdaP);
		}

		/* Stop if --proopt */
		if (OPT_ENABLED(proopt))
			return;

		/* Normalize X */
		// X = sqrt(ncase) * zscore(X)/sqrt(ncase-1)
		M_Xt.normalize(sqrt(N_anaSamp/(wsReal)(N_anaSamp-1)));
	#ifdef PHARAOHvalidate
		M_Xt.file("gsca.normX");
	#endif

		/* Get true */
		wsVec	Ra_truePwCoef	= sseVector(N_latent*N_pheno);
		wsVec	Ra_trueGnCoef	= sseVector(N_sumWidx);
		wsMat	Ra_trueGnPwCoef	= sseMatrix(N_pheno, N_sumWidx);
		/* Get perm */
		vector<wsVec>	Rv_permPwCoef;//	= sseMatrix(N_perm, N_latent);
		vector<wsVec>	Rv_permGnCoef;//	= sseMatrix(N_perm, N_sumWindex);
	//	LOG("PHARAOH with [%d] thread(s) and [%d] permutations\n", OPT_NUMBER(thread),
	//		N_perm);
		if (IS_MULTITHREAD) {
			vLambda Rv_lambdas;
			vReal Rv_lambda;
			Rv_lambda.push_back(R_lambdaG);
			Rv_lambda.push_back(R_lambdaP);
			Rv_lambdas.push_back(Rv_lambda);
			vInt tmp;
			vvInt tmp2;
			xGeSCAres tmp3;
			xGscaAnaThread X_gd = {
				B_overlap,
				Nv_dist,
				N_mnfs,
				&M_Xt, &M_W0,
				M_Yt,
				Xm_mnfs2idx,
				Sv_idx2mnfs,
				Sv_latent, Smv_mnfsInLatent,
				windex,

				Rv_lambdas,
				&N_anaSamp,

				0, NULL, NULL, NULL,

				Ra_truePwCoef, Ra_trueGnCoef,
				Rv_permPwCoef, Rv_permGnCoef,

				Cp_perm,
				M_phenoRel,

				NULL, NULL, tmp2, tmp, tmp3, 0, 0,
			};
			xAnaPermThread X_dt = { Cp_IO, &X_gd };
			WORKER().run(thr_PHARAOH, forAllPerm_equal, &X_dt, NULL, sizeof(int)*3);
		} else {
			cStdMatrix M_permCov(N_cov, N_anaSamp);

			/* Perform PHARAOH for true */
			int it = do_PHARAOH(B_overlap, R_lambdaG, R_lambdaP, M_Xt, M_Yt,
				M_phenoRel, N_anaSamp, Nv_dist, N_mnfs,
				Sv_latent, Smv_mnfsInLatent, M_W0,
				windex, Ra_truePwCoef, Ra_trueGnCoef, NULL, OPT_NUMBER(promaxiter)*2);
//			V_y.file("y");
//			LOG("sum %g\n", V_y.sum());
//			LOG("sumX %g\n", M_Xt.sum());
			if (it <= 0) {
				if (R_lambdaG == R_lambdaP)
					halt("True estimation has failed. Maybe the current lambda [%g]"
						" is too low?", R_lambdaG);
				else
					halt("True estimation has failed. Maybe the current lambda [%g,%g]"
					" is too low?", R_lambdaG, R_lambdaP);
			}
			LOG("True estimation over after [%d] iterations\n", it);

			/* Perform PHARAOH for permutations */
			cStdMatrix M_permYt(N_pheno, N_anaSamp);
			//cVector V_permY(N_anaSamp);
			wsUint	i			= 0;
			char B_success = 0;
			for (Cp_perm->fetch(M_permYt, M_permCov) ; Cp_perm->fetch(M_permYt, M_permCov) == true ;
				i++) {

				/* FIXME: Should handle M_permCov */

//				LOG("sum %g\n", V_permY.sum());
				for (wsUint j=0 ; j<3 ; j++) {
					cStdMatrix&	M_curXt = M_Xt.clone();
//					LOG("sumX %g\n", M_curXt.sum());

					/* 150901 Permute cov if required */
					if (OPT_ENABLED(propermcov)) {
	//					wsMat Ra_curXt = M_curXt.get();
	//					LOOP (j, N_cov)
	//						sseVpermSelf(Ra_curXt[N_gene+j], N_anaSamp,
	//							Nv_seqPerm[i]);
					}

					wsVec Ra_curPermPwCoef = sseVector((wsUint)Smv_mnfsInLatent.size() * N_pheno);
					wsVec Ra_curPermGnCoef = sseVector(N_sumWidx);
//					V_permY.file("yp");
					it = do_PHARAOH(B_overlap, R_lambdaG, R_lambdaP, M_curXt, M_permYt,
						M_phenoRel, N_anaSamp, Nv_dist, N_mnfs,
						Sv_latent, Smv_mnfsInLatent,
						M_W0,
						windex, Ra_curPermPwCoef, Ra_curPermGnCoef, NULL,
						OPT_NUMBER(promaxiter)*2);
					delete& M_curXt;
					if (it != 0) {
						Rv_permPwCoef.push_back(Ra_curPermPwCoef);
						Rv_permGnCoef.push_back(Ra_curPermGnCoef);
						B_success = 1;
						break;
					} else {
						LOG("Permutation [%d] failed, retry...\n", i);
						sseFree(Ra_curPermPwCoef);
						sseFree(Ra_curPermGnCoef);
						Cp_perm->fetchAgain(M_permYt, M_permCov);

						/* FIXME: Should handle M_permCov */
					}
				}
				if (!B_success)
					LOGwarn("Estimation failed with lambda_G [%g] and lambda_P [%g]\n",
						R_lambdaG, R_lambdaP);
				LOG("Permutation [%d] over after [%d] iterations\n", i, it);
			}
		}

		wsUint N_totLatent = N_latent * N_pheno;
		/* Compute trueZ for each phenotype */ LOOP (i, N_pheno) {
			wsVec Ra_trueAZ = vectorRepeat(Ra_truePwCoef + i*N_latent, Smv_mnfsInLatent, N_sumWidx);
			sseVpV(Ra_trueGnCoef, Ra_trueAZ, Ra_trueGnPwCoef[i], N_sumWidx);
			sseFree(Ra_trueAZ);
		}
		exportVector("pharaoh.trueA", Ra_truePwCoef, N_totLatent);
		exportMatrix("pharaoh.trueZ", Ra_trueGnPwCoef, N_pheno, N_sumWidx);
		exportMatrix("pharaoh.permA", Rv_permPwCoef, N_totLatent);
		if (OPT_ENABLED(verbose)) {
			exportVector("pharaoh.trueW", Ra_trueGnCoef, N_sumWidx);
			exportMatrix("pharaoh.permW", Rv_permGnCoef, N_sumWidx);
		}
		wsUintCst	N_perm		= (wsUint)Rv_permPwCoef.size();

		/* 170115 get correlation matrix of beta */
		vector<wsSym>	Rv_pwayCor;
		vector<wsSym>	Rv_geneCor;
		if (N_pheno > 1) {
			LOOP (i, N_pathway) {
				wsMat Ra_beta = sseMatrix(N_pheno, N_perm);
				LOOP (k, N_pheno) LOOP(j, N_perm) Ra_beta[k][j] = Rv_permPwCoef[j][k*N_latent + i];
				wsSym Ra_cor = sseMcorRow(Ra_beta, N_pheno, N_perm);
				Rv_pwayCor.push_back(Ra_cor);
				sseUnmat(Ra_beta, N_pheno);
			}
			wsUint I = 0, J = 0;
			FOREACHDO (vStr_it, Sv_latent, i, I++) {
				vInt& Nv_curWidx = windex[I];

				FOREACH (vInt_it, Nv_curWidx, L) {
					wsMat Ra_w = sseMatrix(N_pheno, N_perm);
					LOOP (k, N_pheno) LOOP(j, N_perm) Ra_w[k][j] = Rv_permPwCoef[j][k*N_latent + I] * Rv_permGnCoef[j][J];
					J++;
					wsSym Ra_cor = sseMcorRow(Ra_w, N_pheno, N_perm);
					Rv_geneCor.push_back(Ra_cor);
					sseUnmat(Ra_w, N_pheno);
				}
		}
		}
		/* 170129 get correlation matrix of w */

		wsVec	Ra_cpTruePwCoef	= NULL;
		wsVec	Ra_sumPwNeg		= sseEmptyVec(N_totLatent);
		wsVec	Ra_sumPwPos		= sseEmptyVec(N_totLatent);
		wsVec	Ra_sePwDist		= sseVector(N_totLatent);
		/* Theoretical A */ {
			wsMat	Ra_pwDist	= sseMatrix(N_totLatent, N_perm);
			wsUint*	Na_negPw	= NULL;
			wsUint*	Na_posPw	= NULL;
			wsUint*	Na_pwDist	= NULL;
			wsCalloc(Na_negPw, wsUint, N_totLatent);
			wsCalloc(Na_posPw, wsUint, N_totLatent);
			wsCalloc(Na_pwDist, wsUint, N_totLatent);
			LOOP (i, N_perm) LOOP (j, N_totLatent) {
				wsRealCst R_v = Rv_permPwCoef[i][j];
				if (R_v < 1e-5 && R_v > -1e-5) continue;
				Ra_pwDist[j][Na_pwDist[j]++] = R_v;
				if (R_v < W0) {
					Ra_sumPwNeg[j] += R_v;
					Na_negPw[j]++;
				} else {
					Ra_sumPwPos[j] += R_v;
					Na_posPw[j]++;
				}
			}
			/* Centering */
			LOOP (j, N_totLatent) {
				Ra_sumPwNeg[j] /= (wsReal)Na_negPw[j];
				Ra_sumPwPos[j] /= (wsReal)Na_posPw[j];
			}
	//		LOG("Pathway coefficient +/- : %g/%g\n", R_sumPwPos, R_sumPwNeg);
			LOOP (j, N_totLatent) {
				wsRealCst R_sumPwNeg = Ra_sumPwNeg[j];
				wsRealCst R_sumPwPos = Ra_sumPwPos[j];
				LOOP(i, Na_pwDist[j])
					Ra_pwDist[j][i] -= Ra_pwDist[j][i] < W0 ? R_sumPwNeg : R_sumPwPos;
				Ra_sePwDist[j] = sseVse(Ra_pwDist[j], Na_pwDist[j]);
			}
			Ra_cpTruePwCoef = sseVectorP(N_totLatent, Ra_truePwCoef);
			sseUnmat(Ra_pwDist, N_totLatent);

			DEALLOC(Na_pwDist);
			DEALLOC(Na_negPw);
			DEALLOC(Na_posPw);
		}

#ifdef _DEBUG
		wsVec	Ra_cpTrueGnCoef = NULL;
		wsReal	R_sumGnNeg = W0, R_sumGnPos = W0;
		wsReal	R_seGnDist		= W0;
#endif
		wsMat	Ra_cpTrueGnPwCoef = NULL;
		wsAlloc(Ra_cpTrueGnPwCoef, wsVec, N_pheno);
		wsVec	Ra_sumGnPwNeg		= sseVector(N_pheno);
		wsVec	Ra_sumGnPwPos		= sseVector(N_pheno);
		wsVec	Ra_seGnPwDist		= sseVector(N_pheno);

#ifdef _DEBUG
		/* Theoretical W */ {
			wsVec	Ra_gnDist = sseVector(N_perm * N_sumWidx);
			wsUint	N_negGn = 0, N_posGn = 0, N_gnDist = 0;
			LOOP(i, N_perm) LOOP(j, N_sumWidx) {
				wsRealCst R_vg = Rv_permGnCoef[i][j];
				if (NA(R_vg)) continue;
				if (R_vg < 1e-5 && R_vg > -1e-5) continue;
				Ra_gnDist[N_gnDist++] = R_vg;
				if (R_vg < W0) {
					R_sumGnNeg += R_vg;
					N_negGn++;
				} else {
					R_sumGnPos += R_vg;
					N_posGn++;
				}
			}
			/* Centering */
			R_sumGnNeg /= (wsReal)N_negGn;
			R_sumGnPos /= (wsReal)N_posGn;
			LOG("Gene coefficient +/- : %g/%g\n", R_sumGnPos, R_sumGnNeg);
			LOOP(i, N_gnDist) Ra_gnDist[i] -= Ra_gnDist[i] < W0 ? R_sumGnNeg : R_sumGnPos;
			R_seGnDist = sseVse(Ra_gnDist, N_gnDist);
			Ra_cpTrueGnCoef = sseVectorP(N_sumWidx, Ra_trueGnCoef);
			sseFree(Ra_gnDist);
		}
#endif

		wsUint	N_rankLatent	= N_pathway;							// Effective number of latent for rank calculation
		wsMat	Ra_pvPwCoef		= sseMatrix(N_pheno, N_latent);			// p-value, pathway, coefficient, permutation-based
		wsMat	Ra_pvGnPwCoef	= sseMatrix(N_pheno, N_sumWidx);		// p-value, gene-global, coefficient, permutation-based
		wsMat	Ra_pvPwRank		= sseMatrix(N_pheno, N_rankLatent);		// p-value, pathway, rank, permutation-based
		wsMat	Ra_pvPwZrank	= sseMatrix(N_pheno, N_rankLatent);		// p-value, pathway, rank, z-transform
		wsUmat	Na_loQpwCoef	= sseUmatrix(N_pheno, N_latent);
		wsUmat	Na_hiQpwCoef	= sseUmatrix(N_pheno, N_latent);
		wsUmat	Na_loQgnpwCoef	= sseUmatrix(N_pheno, N_sumWidx);
		wsUmat	Na_hiQgnpwCoef	= sseUmatrix(N_pheno, N_sumWidx);
#ifdef _DEBUG
		wsVec	Ra_pvGnCoef		= sseVector(N_sumWidx);					// p-value, gene-local, coefficient, permutation-based
#endif
#ifdef _DEBUG
		wsVec	Ra_pvGnRank		= sseVector(N_sumWidx);					// p-value, gene-local, rank, permutation-based
		wsVec	Ra_pvGnZrank	= sseVector(N_sumWidx);					// p-value, gene-local, rank, z-transform
#endif
		wsMat	Ra_pvGnPwRank	= sseMatrix(N_pheno, N_sumWidx);		// p-value, gene-global, rank, permutation-based
		wsMat	Ra_pvGnPwZrank	= sseMatrix(N_pheno, N_sumWidx);		// p-value, gene-global, rank, z-transform

		wsVec	Ra_trueGnRank	= sseVarank(Ra_trueGnCoef, N_sumWidx);
		LOOP (k, N_pheno) {
			wsUint	N_idxCur		= k*N_latent;
			wsVec	Ra_curTruePwCoef= Ra_truePwCoef + N_idxCur;

			memset(Na_loQpwCoef[k], 0x00, sizeof(wsUint)*N_latent);
			memset(Na_hiQpwCoef[k], 0x00, sizeof(wsUint)*N_latent);
			memset(Na_loQgnpwCoef[k], 0x00, sizeof(wsUint)*N_sumWidx);
			memset(Na_hiQgnpwCoef[k], 0x00, sizeof(wsUint)*N_sumWidx);

			/* Ranks for competitive method */
			wsVec	Ra_truePwRank	= sseVarank(Ra_curTruePwCoef, N_rankLatent);
			wsVec	Ra_trueGnPwRank	= sseVarank(Ra_trueGnPwCoef[k], N_sumWidx);
			wsMat	Ra_permPwRank	= NULL;
			wsMat	Ra_permGnRank	= NULL;
			wsMat	Ra_permGnPwRank	= NULL;
			wsCalloc(Ra_permPwRank, wsVec, N_perm);
			wsCalloc(Ra_permGnRank, wsVec, N_perm);
			wsCalloc(Ra_permGnPwRank, wsVec, N_perm);

			/* Get the p-values */
			wsUint	N_mnfst	= N_sumWidx;
#ifdef _DEBUG
			wsUint*	Na_loQgnCoef	= NULL;
			wsUint*	Na_hiQgnCoef	= NULL;
			wsCalloc(Na_loQgnCoef, wsUint, N_mnfst);
			wsCalloc(Na_hiQgnCoef, wsUint, N_mnfst);
#endif

			LOOP(i, N_perm) {
				wsVecCst	Ra_curPwCoef	= Rv_permPwCoef[i] + N_idxCur;
				wsVecCst	Ra_curGnCoef	= Rv_permGnCoef[i];
				wsVec		Ra_curGnPwCoef	= sseVector(N_mnfst);

				/* Pathway/gene-level rank */
				Ra_permPwRank[i] = sseVarank(Ra_curPwCoef, N_rankLatent);
				Ra_permGnRank[i] = sseVarank(Ra_curGnCoef, N_mnfst);

				/* Pathway-level coefficient */
				LOOP(j, N_latent) {
					Na_loQpwCoef[k][j] += Ra_curTruePwCoef[j] <= Ra_curPwCoef[j];
					Na_hiQpwCoef[k][j] += Ra_curTruePwCoef[j] >= Ra_curPwCoef[j];
				}

#ifdef _DEBUG
				/* Gene-level coefficient */
				LOOP(j, N_mnfst) {
					Na_loQgnCoef[j] += Ra_trueGnCoef[j] <= Ra_curGnCoef[j];
					Na_hiQgnCoef[j] += Ra_trueGnCoef[j] >= Ra_curGnCoef[j];
				}
#endif

				/* Build expanded A */
				wsVec Ra_permWZ = vectorRepeat(Rv_permPwCoef[i] + N_idxCur, Smv_mnfsInLatent, N_sumWidx);
				sseVpV(Ra_curGnCoef, Ra_permWZ, Ra_curGnPwCoef, N_sumWidx);
				sseFree(Ra_permWZ);

				/* Gene-level coefficient */
				LOOP(j, N_mnfst) {
					Na_loQgnpwCoef[k][j] += Ra_trueGnPwCoef[k][j] <= Ra_curGnPwCoef[j];
					Na_hiQgnpwCoef[k][j] += Ra_trueGnPwCoef[k][j] >= Ra_curGnPwCoef[j];
				}

				/* Expanded A coeff rank */
				Ra_permGnPwRank[i] = sseVarank(Ra_curGnPwCoef, N_mnfst);

				sseFree(Ra_curGnPwCoef);
			}
			char S_fn[256];
			sprintf(S_fn, "pharaoh.phen%d.permZ", k);
			exportMatrix(S_fn, Rv_permGnCoef, N_sumWidx);

			/* Theoretical Z/W */ {
				wsVec	Ra_distGnPw = sseVector(N_perm * N_sumWidx);
				wsUint	N_negGnPw = 0, N_posGnPw = 0, N_distGnPw = 0;
				LOOP(i, N_perm) LOOP(j, N_sumWidx) {
					wsRealCst R_vgp = Rv_permGnCoef[i][j];
					if (NA(R_vgp) || (R_vgp < 1e-5 && R_vgp > -1e-5)) continue;
					Ra_distGnPw[N_distGnPw++] = R_vgp;

					if (R_vgp < W0) {
						Ra_sumGnPwNeg[k] += R_vgp;
						N_negGnPw++;
					} else {
						Ra_sumGnPwPos[k] += R_vgp;
						N_posGnPw++;
					}
				}
				/* Centering */
				Ra_sumGnPwPos[k] /= (wsReal)N_posGnPw;
				Ra_sumGnPwNeg[k] /= (wsReal)N_negGnPw;
				LOG("Gene-pathway coefficient +/- : %g/%g\n", Ra_sumGnPwPos, Ra_sumGnPwNeg);
				LOOP(i, N_distGnPw) Ra_distGnPw[i] -= Ra_distGnPw[i] < W0 ? Ra_sumGnPwNeg[k] : Ra_sumGnPwPos[k];
				Ra_seGnPwDist[k] = sseVse(Ra_distGnPw, N_distGnPw);
				Ra_cpTrueGnPwCoef[k] = sseVectorP(N_sumWidx, Ra_trueGnPwCoef[k]);
				sseFree(Ra_distGnPw);
			}

			LOOP(i, N_latent) {
				wsReal R_Ap1 = (Na_loQpwCoef[k][i] + 1) / (wsReal)(N_perm + 1);
				wsReal R_Ap2 = (Na_hiQpwCoef[k][i] + 1) / (wsReal)(N_perm + 1);
				Ra_pvPwCoef[k][i] = min(min(R_Ap1, R_Ap2) * W2, W1);
			}
			LOOP(i, N_mnfst) {
#ifdef _DEBUG
				if (k == 0) { // Since it is indep. w/ k
					wsReal R_Wp1 = (Na_loQgnCoef[i] + 1) / (wsReal)(N_perm + 1);
					wsReal R_Wp2 = (Na_hiQgnCoef[i] + 1) / (wsReal)(N_perm + 1);
					Ra_pvGnCoef[i] = min(min(R_Wp1, R_Wp2) * W2, W1);
				}
#endif
				wsReal R_Zp1 = (Na_loQgnpwCoef[k][i] + 1) / (wsReal)(N_perm + 1);
				wsReal R_Zp2 = (Na_hiQgnpwCoef[k][i] + 1) / (wsReal)(N_perm + 1);
				Ra_pvGnPwCoef[k][i] = min(min(R_Zp1, R_Zp2) * W2, W1);
			}
#ifdef _DEBUG
			DEALLOC(Na_loQgnCoef);
			DEALLOC(Na_hiQgnCoef);
#endif

			/* Get the competitive p-values */
			wsUint*	Na_qPwRank		= NULL;
			wsCalloc(Na_qPwRank, wsUint, N_rankLatent);
			// 	wsUint*	Na_Az	= NULL;
			// 	wsCalloc(Na_Az, wsUint, N_latent);
#ifdef _DEBUG
			wsUint*	Na_qGnRank		= NULL;
			wsCalloc(Na_qGnRank, wsUint, N_mnfst);
#endif
			wsUint*	Na_qGnPwRank	= NULL;
			wsCalloc(Na_qGnPwRank, wsUint, N_mnfst);

			LOOP(i, N_perm) {
				wsVecCst	Ra_curA = Ra_permPwRank[i];
				wsVecCst	Ra_curW = Ra_permGnRank[i];
				wsVec		Ra_curZ = Ra_permGnPwRank[i];

				/* Pathway-level coefficient */
				LOOP(j, N_rankLatent)	Na_qPwRank[j] += Ra_truePwRank[j] <= Ra_curA[j];
#ifdef _DEBUG
				/* Gene-local coefficient */
				LOOP(j, N_mnfst)	Na_qGnRank[j] += Ra_trueGnRank[j] <= Ra_curW[j];
#endif
				/* Gene-marginal coefficient */
				LOOP(j, N_mnfst)	Na_qGnPwRank[j] += Ra_trueGnPwRank[j] <= Ra_curZ[j];
			}
			if (OPT_ENABLED(verbose)) {
				exportMatrix("pharaoh.permAr", Ra_permPwRank, N_perm, N_rankLatent);
				exportMatrix("pharaoh.permWr", Ra_permGnRank, N_perm, N_sumWidx);
				exportMatrix("pharaoh.permZr", Ra_permGnPwRank, N_perm, N_sumWidx);
			}

			/* Standardize rank matrix */ {
				wsVec Ra_meanPwRank = NULL, Ra_sdPwRank = NULL;
				sseMsdmeanC(Ra_permPwRank, N_perm, N_rankLatent, &Ra_meanPwRank, &Ra_sdPwRank);
				/* Stdize true */ {
					sseVsV(Ra_truePwRank, Ra_meanPwRank, Ra_truePwRank, N_rankLatent);
					sseVdV(Ra_truePwRank, Ra_sdPwRank, Ra_truePwRank, N_rankLatent);
				}
				// 		LOOP(i, N_perm) {
				// 			VEC_t Ra_curA = Ra_Ar[i];
				// 			sseVsV(Ra_curA, Ra_mAr, Ra_curA, N_rankLatent);
				// 			sseVdV(Ra_curA, Ra_sdAr, Ra_curA, N_rankLatent);
				//
				// 			/* Pathway-level coefficient */
				// 			LOOP(j, N_rankLatent) Na_Az[j] += Ra_trueAr[j] <= Ra_curA[j];
				// 		}
				sseFree(Ra_meanPwRank);
				sseFree(Ra_sdPwRank);

#ifdef _DEBUG
				wsVec Ra_meanGnRank = NULL, Ra_sdGnRank = NULL;
				sseMsdmeanC(Ra_permGnRank, N_perm, N_sumWidx,
					&Ra_meanGnRank, &Ra_sdGnRank);
				/* Stdize true */ {
					sseVsV(Ra_trueGnRank, Ra_meanGnRank, Ra_trueGnRank, N_sumWidx);
					sseVdV(Ra_trueGnRank, Ra_sdGnRank, Ra_trueGnRank, N_sumWidx);
				}
				sseFree(Ra_meanGnRank);
				sseFree(Ra_sdGnRank);
#endif

				wsVec Ra_meanGnPwRank = NULL, Ra_sdGnPwRank = NULL;
				sseMsdmeanC(Ra_permGnPwRank, N_perm, N_sumWidx,
					&Ra_meanGnPwRank, &Ra_sdGnPwRank);
				/* Stdize true */ {
					sseVsV(Ra_trueGnPwRank, Ra_meanGnPwRank, Ra_trueGnPwRank, N_sumWidx);
					sseVdV(Ra_trueGnPwRank, Ra_sdGnPwRank, Ra_trueGnPwRank, N_sumWidx);
				}
				sseFree(Ra_meanGnPwRank);
				sseFree(Ra_sdGnPwRank);
			}

			sseUnmat(Ra_permPwRank, N_perm);
			sseUnmat(Ra_permGnRank, N_perm);
			sseUnmat(Ra_permGnPwRank, N_perm);

			LOOP(i, N_rankLatent) {
				Ra_pvPwRank[k][i] = (Na_qPwRank[i] + 1) / (wsReal)(N_perm + 1);
				if (Ra_pvPwRank[k][i] > W1) Ra_pvPwRank[k][i] = W1;
				Ra_pvPwZrank[k][i] = W1 - pnorm(Ra_truePwRank[i]);
				//			(Na_Az[i] + 1) / (wsReal)(N_perm + 1);
				//		if (Ra_zApv[i] > W1) Ra_zApv[i] = W1;
			}
			LOOP(i, N_mnfst) {
#ifdef _DEBUG
				if (k == 0) { // Since it is indep. w/ k
					Ra_pvGnRank[i] =
						(Na_qGnRank[i] + 1) / (wsReal)(N_perm + 1);
					Ra_pvGnZrank[i] = W1 - pnorm(Ra_trueGnRank[i]);
				}
#endif
				Ra_pvGnPwRank[k][i] =
					(Na_qGnPwRank[i] + 1) / (wsReal)(N_perm + 1);
				Ra_pvGnPwZrank[k][i] = W1 - pnorm(Ra_trueGnPwRank[i]);
				//		if (Ra_rWpv[i] > W1) Ra_rWpv[i] = W1;
				//		if (Ra_rZpv[i] > W1) Ra_rZpv[i] = W1;
			}
			sseFree(Ra_truePwRank);
			sseFree(Ra_trueGnPwRank);
			DEALLOC(Na_qPwRank);
#ifdef _DEBUG
			DEALLOC(Na_qGnRank);
#endif
			DEALLOC(Na_qGnPwRank);
		}
		sseFree(Ra_trueGnRank);
		sseUnmat(Na_loQpwCoef, N_pheno);
		sseUnmat(Na_hiQpwCoef, N_pheno);
		sseUnmat(Na_loQgnpwCoef, N_pheno);
		sseUnmat(Na_hiQgnpwCoef, N_pheno);
		FOREACH(vector<wsVec>::iterator, Rv_permPwCoef, i)
			sseFree(*i);
		FOREACH(vector<wsVec>::iterator, Rv_permGnCoef, i)
			sseFree(*i);

		/* 170115 multiple phenotype */
		if (N_pheno > 1) {
			cTableExporter C_exMulPway("pharaoh.multi.res", "siiirrr",
				"PHARAOH multiple phenotype result summary", 0,
				7, "PATHWAY", "NPERM", "NGENE", "NVARIANT",
				"P_FISHER", "P_STOUFFER", "P_KOST");
			cTableExporter C_exMulGene("pharaoh.gene.multi.res", "ssiirr",
				"PHARAOH multiple phenotype gene-level result summary", 0,
				6, "PATHWAY", "GENE", "NPERM", "NVARIANT",
				"P_FISHER", "P_KOST");

				wsVec Ra_pval = sseVector(N_pheno);
				wsVec Ra_coef = sseVector(N_pheno);
			LOOP(i, N_pathway) {
				// Gather p-values
				LOOP(k, N_pheno) {
					Ra_pval[k] = Ra_pvPwCoef[k][i];
					Ra_coef[k] = Ra_truePwCoef[k*N_latent + i];
				}
				// Get Fisher
				wsReal R_fisher = metaFisher(N_pheno, Ra_pval);
				wsReal R_stouffer = metaStouffer(N_pheno, Ra_pval, Ra_coef);
				wsReal R_kost = metaKost(N_pheno, Ra_pval, Rv_pwayCor[i]);

				C_exMulPway.write(7, Sv_latent[i].c_str(), N_perm, N_gene,
					Nv_variants[i], R_fisher, R_stouffer, R_kost);
			}

			wsUint I = 0, J = 0;
			FOREACHDO (vStr_it, Sv_latent, i, I++) {
				vStr& Sv_mnfsIdx = Smv_mnfsInLatent[*i];
				vInt& Nv_curWidx = windex[I];

				wsUint l = 0;
				FOREACHDO (vInt_it, Nv_curWidx, L, l++) {
					string& S_gene = Sv_mnfsIdx[l];

					LOG("%s -> %s :", i->c_str(), S_gene.c_str());

					// Gather p-values
					LOOP (k, N_pheno) {
						LOGnf(" %g", Ra_pvGnPwCoef[k][J]);
						Ra_pval[k] = Ra_pvGnPwCoef[k][J];
					}
					LOGnf("\n");

					// Get Fisher
					wsReal R_fisher = metaFisher(N_pheno, Ra_pval);
					wsReal R_kost = metaKost(N_pheno, Ra_pval, Rv_geneCor[J++]);

					C_exMulGene.write(6, i->c_str(), S_gene.c_str(), N_perm, Xm_nvrtGene[S_gene],
						R_fisher, R_kost);
				}
			}
				sseFree(Ra_pval);
				sseFree(Ra_coef);
		}

		FOREACH(vector<wsSym>::iterator, Rv_pwayCor, i) sseUnmat(*i, N_pheno);
		FOREACH(vector<wsSym>::iterator, Rv_geneCor, i) sseUnmat(*i, N_pheno);

		sseFree(Ra_truePwCoef);
		sseFree(Ra_trueGnCoef);
		sseUnmat(Ra_trueGnPwCoef, N_pheno);

		/* Print out result */
		//LOGoutput("PHARAOH result is exported to [%s.gsca.res]\n", OPT_STRING(out));
		wsUint k = 0;
		vPheno& Xv_pheno = getIO()->getPhenoInfo();
		FOREACHDO(vPheno_it, Xv_pheno, P, k++) {
			char S_fn[512];

			/* Pathway-level result file */
			if (N_pheno == 1) sprintf(S_fn, "pharaoh.pathway.res");
			else sprintf(S_fn, "pharaoh.pathway.%s.res", P->S_name.c_str());

			cTableExporter	C_exP(S_fn, "siiirrrr",
				"PHARAOH pathway-level result summary", 0,
				8, "PATHWAY", "NPERM", "NGENE", "NVARIANT",
				"P_PHARAOH_SC", "P_PHARAOH_SCZ",
				"P_PHARAOH_CP", "P_PHARAOH_CPZ");

			/* Gene-level result file */
			if (N_pheno == 1) sprintf(S_fn, "pharaoh.gene.res");
			else sprintf(S_fn, "pharaoh.gene.%s.res", P->S_name.c_str());
#ifdef _DEBUG
			cTableExporter	C_exG(S_fn, "ssiirrrrrrrr",
				"PHARAOH gene-level result summary", 0,
				12, "PATHWAY", "GENE", "NPERM", "NVARIANT",
				"P_PHARAOH_SC_GENE_LOCAL", "P_PHARAOH_SC_GENE_MARGINAL",
				"P_PHARAOH_SCZ_GENE_LOCAL", "P_PHARAOH_SCZ_GENE_MARGINAL",
				"P_PHARAOH_CP_GENE_LOCAL", "P_PHARAOH_CP_GENE_MARGINAL",
				"P_PHARAOH_CPZ_GENE_LOCAL", "P_PHARAOH_CPZ_GENE_MARGINAL");
#else
			cTableExporter	C_exG(S_fn, "ssiirrrr",
				"PHARAOH gene-level result summary", 0,
				8, "PATHWAY", "GENE", "NPERM", "NVARIANT",
				"P_PHARAOH_SC_GENE_MARGINAL", "P_PHARAOH_SCZ_GENE_MARGINAL",
				"P_PHARAOH_CP_GENE_MARGINAL", "P_PHARAOH_CPZ_GENE_MARGINAL");
#endif

			wsUint I = 0;
			/* Insert covariate names */
			FOREACH(vCovar_it, Xv_cov, i) {
				char S_covName[256];
				sprintf(S_covName, "COV_%s", i->Sp_varName);
				Sv_latent.push_back(S_covName);
			}
			LOOP(i, N_latent) {
				if (i >= (Smv_mnfsInLatent.size() - N_cov*N_pheno)) {
					wsReal R_theoPv = fabs(Ra_cpTruePwCoef[i]) < 1e-5 ? WISARD_NAN :
						(Ra_cpTruePwCoef[i] > W0 ?
						pnorm((Ra_cpTruePwCoef[i] - Ra_sumPwPos[i]) / Ra_sePwDist[i], W0, W1, 0) :
						pnorm((Ra_cpTruePwCoef[i] - Ra_sumPwNeg[i]) / Ra_sePwDist[i]));

					/* For covariates */
					C_exP.write(8, Sv_latent[i].c_str(), N_perm, WIS_I32NA,
						WIS_I32NA, Ra_pvPwCoef[k][i], R_theoPv,
						WISARD_NAN, WISARD_NAN);
				} else {
					wsReal R_theoPv = fabs(Ra_cpTruePwCoef[i]) < 1e-5 ? WISARD_NAN :
						(Ra_cpTruePwCoef[i] > W0 ?
						pnorm((Ra_cpTruePwCoef[i] - Ra_sumPwPos[i]) / Ra_sePwDist[i], W0, W1, 0) :
						pnorm((Ra_cpTruePwCoef[i] - Ra_sumPwNeg[i]) / Ra_sePwDist[i]));

					/* For genes */
					string	S_latent	= Sv_latent[i];
					wsUint	N_gene		= (wsUint)Smv_mnfsInLatent[S_latent].size();
					if (i < N_rankLatent)
						C_exP.write(8, S_latent.c_str(), N_perm, N_gene,
							Nv_variants[i], Ra_pvPwCoef[k][i], R_theoPv,
							Ra_pvPwRank[k][i], Ra_pvPwZrank[k][i]);
					else
						C_exP.write(8, S_latent.c_str(), N_perm, N_gene,
						Nv_variants.size() <= i ? 0 : Nv_variants[i], Ra_pvPwCoef[k][i], R_theoPv,
							WISARD_NAN, WISARD_NAN);
					LOOP(j, N_gene) {
						//char S_buf[1024];
						wsReal R_theoGnPwPv = fabs(Ra_cpTrueGnPwCoef[k][j]) < 1e-5 ? WISARD_NAN :
							(Ra_cpTrueGnPwCoef[k][j] > W0 ?
							pnorm((Ra_cpTrueGnPwCoef[k][j] - Ra_sumGnPwPos[k]) / Ra_seGnPwDist[k], W0, W1, 0) :
							pnorm((Ra_cpTrueGnPwCoef[k][j] - Ra_sumGnPwNeg[k]) / Ra_seGnPwDist[k]));
#ifdef _DEBUG
						wsReal R_theoGnPv = fabs(Ra_cpTrueGnCoef[j]) < 1e-5 ? WISARD_NAN :
							(Ra_cpTrueGnCoef[j] > W0 ?
							pnorm((Ra_cpTrueGnCoef[j] - R_sumGnPos) / R_seGnDist, W0, W1, 0) :
							pnorm((Ra_cpTrueGnCoef[j] - R_sumGnNeg) / R_seGnDist));

						C_exG.write(12, S_latent.c_str(),
							Sv_idx2mnfs[I].c_str(), N_perm, Xm_nvrtGene[Sv_idx2mnfs[I]],
							Ra_pvGnCoef[I], Ra_pvGnPwCoef[k][I], R_theoGnPv, R_theoGnPwPv,
							Ra_pvGnRank[I], Ra_pvGnPwRank[k][I], Ra_pvGnZrank[I], Ra_pvGnPwZrank[k][I]);
#else
							C_exG.write(8, S_latent.c_str(),
								Sv_idx2mnfs[I].c_str(), N_perm, Xm_nvrtGene[Sv_idx2mnfs[I]],
								Ra_pvGnPwCoef[k][I], R_theoGnPwPv,
								Ra_pvGnPwRank[k][I], Ra_pvGnPwZrank[k][I]);
#endif
							I++;
					}
				}
			}
		}
		sseFree(Ra_sumPwNeg);
		sseFree(Ra_sumPwPos);
		sseFree(Ra_sePwDist);
		sseFree(Ra_cpTruePwCoef);
		sseUnmat(Ra_cpTrueGnPwCoef, N_pheno);

		sseFree(Ra_sumGnPwNeg);
		sseFree(Ra_sumGnPwPos);
		sseFree(Ra_seGnPwDist);


		sseUnmat(Ra_pvPwCoef, N_pheno);
		sseUnmat(Ra_pvPwRank, N_pheno);
		sseUnmat(Ra_pvPwZrank, N_pheno);
		sseUnmat(Ra_pvGnPwZrank, N_pheno);
	#ifdef _DEBUG
		sseFree(Ra_cpTrueGnCoef);
		sseFree(Ra_pvGnCoef);
		sseFree(Ra_pvGnRank);
		sseFree(Ra_pvGnZrank);
	#endif
		sseUnmat(Ra_pvGnPwCoef, N_pheno);
		sseUnmat(Ra_pvGnPwRank, N_pheno);
		/* For each gene-set */
	// 	FOREACH (mGeneSet_it, Xm_gs, i)
	// 		_doTest(getIO(), i->second, Xm_gd, getIO()->getSampleSize(),
	// 			getIO()->getGenotype(), REAL_CONST(9.0), NULL, 0);
	}
	/* Cleanup */
	for (map<string, wsVec>::iterator i=Xm_collapsed.begin() ; i != Xm_collapsed.end() ; i++)
		sseFree(i->second);
	if (Sp_oldPrefix)
		OPT_STRING(out) = Sp_oldPrefix;
	delete Cp_perm;
/*	map<wsVec, bool> isRemoved;
	wsMat _Xt = M_Xt.get();
	LOOP(i, M_Xt.row()) {
		if (isRemoved.find(_Xt[i]) != isRemoved.end()) continue;
		isRemoved[_Xt[i]] = true;
		sseFree(_Xt[i]);
	}*/
}

#endif

} // End namespace ONETOOL
