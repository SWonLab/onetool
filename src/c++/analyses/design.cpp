#include "analyses/design.h"
#include "analyses/emai.h"
#include "utils/util.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

#define CONDVAR	MVNcondvarV2

inline wsReal getEffSmpSz(wsSym Ra_cor, wsUint N_sz)
{
	wsReal*	Ra_eval		= EIGENDECOMPOSITION(Ra_cor, N_sz);
	wsReal	R_nEval		= sseVvar(Ra_eval, N_sz);
	wsReal	R_effSmpSz	= W1 + (wsReal)(N_sz - 1)
		* (W1 - R_nEval/(wsReal)N_sz);
	sseFree(Ra_eval);

	return R_effSmpSz;
}

inline char* dupSel(char *Ba_sel, wsUint N_sz)
{
	char*	Ba_nsel	= NULL;
	wsAlloc(Ba_nsel, char, N_sz);
	memcpy(Ba_nsel, Ba_sel, sizeof(char)*N_sz);
	return Ba_nsel;
}

int forAllcomb_anaEqual(int N_mode, int N_thread, void *Vp_data,
	void *Vp_shareData, wsUint *Na_waitQueue)
{
	int			*Np_idx	= (int *)Vp_data;
	xAnaThread	*Cp_ana	= (xAnaThread *)Vp_shareData;
	wsUint		N_comb	= *((wsUint *)Cp_ana->Vp_data);
//	cIO			*Cp_IO	= Cp_ana->Cp_IO;
	static wsUint
		N_proc	= 0;
	wsReal		R_szDiv	= (wsReal)N_comb/(wsReal)N_thread;
	int			N_ret	= N_proc != 0 ? 0 : 1;
	int			i=0;
	wsReal		j=W0;

	switch (N_mode) {
	case DISTFUN_INIT:
		i = 0;
		j = W0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		if (N_ret == 0)
			return 0;
		for ( ; i<N_thread ; i++,j+=R_szDiv) {
			wsUint E = (wsUint)(j+R_szDiv+REAL_CONST(0.5));
			if (E > N_comb)
				E = N_comb;
			Np_idx[3*i] = (wsUint)(j+REAL_CONST(0.5));
			Np_idx[3*i+1] = E;
			Np_idx[3*i+2] = 0;
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

int thrMVNcondvar(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	if (N_idx < 0 || N_idx >= OPT_NUMBER(thread))
		halt("Error index [%d]", N_idx);
	xAnaThread		*Xp_at		= (xAnaThread *)Vp_shareData;
	cStudyAnalysis	*Cp_fa		= (cStudyAnalysis *)(Xp_at->Vp_data);
	int				*Np_data	= (int *)Vp_data;
	cIO				*Cp_IO		= Cp_fa->getIO();
	vVariant			&Xa_snp		= Cp_IO->getVariant();
//	wsUint			N_anaRow	= Cp_IO->getCovariateSize()+2;

	wsUint			N_s			= (wsUint)Np_data[0];
	wsUint			N_e			= (wsUint)Np_data[1];

	for (wsUint Q=N_s ; Q<N_e ; Q++) {
		Np_data[2]++;

		/* Print progress */
		if (N_idx==0 && (Q%100)==0) {
			wsUint N_sum = 0;
			for (int x=2 ; x<OPT_NUMBER(thread) ; x+=3)
				N_sum += Np_data[x];
			LOG("%d/%d SNPs tested...\r", N_sum, Xa_snp.size());
		}
	}

	return 0;
}

cStudyAnalysis::cStudyAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaCor) :
	cAnalysis(Cp_inpIO)
{
	Cp_anaCor = Cp_inpAnaCor;
}

cStudyAnalysis::~cStudyAnalysis()
{

}

wsReal MVNcondvar(wsMat Ra_cor, wsUint N_samp, char *Ba_isInc, wsUint N_tsz,
	wsReal *Rp_effSmpSz=NULL)
{
	wsSym	Ra_22	= NULL;
	wsMat	Ra_12	= NULL;
	wsSym	Ra_11	= sseSplitMatrix(Ra_cor, N_samp, Ba_isInc, N_tsz,
		&Ra_12, &Ra_22);
	cSymMatrix	M_22(Ra_22, N_samp-N_tsz);
	cSymMatrix	M_11(Ra_11, N_tsz);
	cStdMatrix	M_12(N_tsz, N_samp-N_tsz, Ra_12);

	if (Rp_effSmpSz) {
		wsReal*	Ra_eval		= EIGENDECOMPOSITION(Ra_22, N_tsz);
		wsReal	R_vEval		= sseVvar(Ra_eval, N_tsz);
		wsReal	R_effSmpSz	= W1 + (wsReal)(N_tsz - 1)
			* (W1 - R_vEval/(wsReal)N_tsz);
		*Rp_effSmpSz = R_effSmpSz;
	}

	cSymMatrix&	M_22i = M_22.inv();
	M_22.rem();

	cSymMatrix M_va = M_12.MMt(M_22i);
	delete &M_22i;
	M_11 -= M_va;
	wsReal	R_cur = M_11.sum() / (wsReal)SQR(N_tsz);

	return R_cur;
}

/* Ba_isInc == Samples to be removed
 * N_szX1 == # of 1's in Ba_isInc
 */
wsReal MVNcondvarV2(cIO* Cp_IO, wsMat Ra_cor, wsUint N_samp, char *Ba_isInc, wsUint N_szX2,
	wsReal *Rp_effSmpSz=NULL)
{
	wsUint	N_szX1	= N_samp - N_szX2;
	wsMat	Ra_21	= NULL;
	wsSym	Ra_11	= sseSplitMatrixV2(Ra_cor, N_samp, Ba_isInc, N_szX2,
		&Ra_21);
	cSymMatrix	M_11(Ra_11, N_szX1);
	cStdMatrix	M_21(N_szX2, N_szX1, Ra_21);

	if (Rp_effSmpSz) {
		wsReal*	Ra_eval		= EIGENDECOMPOSITION(Ra_11, N_szX1);
		wsReal	R_vEval		= sseVvar(Ra_eval, N_szX1);
		wsReal	R_effSmpSz	= W1 + (wsReal)(N_szX1 - 1)
			* (W1 - R_vEval/(wsReal)N_szX1);
		*Rp_effSmpSz = R_effSmpSz;
	}

	wsReal R_cur = W0;
	if (OPT_ENABLED(kinship)) {
		wsMat Ra_famInv = Cp_IO->getFamInv(Ba_isInc, Ra_cor, 0);
		cSymMatrix M_famInv(Ra_famInv, N_szX1);
		M_11.rem();
		R_cur = M_21.trMSMt(M_famInv);
	} else {
		cSymMatrix&	M_11i = M_11.inv();
		M_11.rem();
		R_cur = M_21.trMSMt(M_11i);
	}
	return R_cur;
}

wsReal MVNcondvar2(wsMat Ra_cor, wsUint N_samp, char *Ba_isExc, wsUint N_tsz)
{
	wsSym	Ra_22	= NULL;
	wsMat	Ra_12	= NULL;
	wsSym	Ra_11	= sseSplitMatrix2(Ra_cor, N_samp, Ba_isExc, N_tsz,
		&Ra_12, &Ra_22);
// 	cSymMatrix	M_22(Ra_22, N_samp-N_tsz);
// 	cSymMatrix	M_11(Ra_11, N_tsz);
// 	cStdMatrix	M_12(N_tsz, N_samp-N_tsz, Ra_12);
	cSymMatrix	M_22(Ra_22, N_tsz);
	cSymMatrix	M_11(Ra_11, N_samp-N_tsz);
	cStdMatrix	M_12(N_samp-N_tsz, N_tsz, Ra_12);

	cSymMatrix&	M_22i = M_22.inv();
	M_22.rem();

	cSymMatrix M_va = M_12.MMt(M_22i);
	delete &M_22i;
	M_11 -= M_va;
	wsReal	R_cur = M_11.sum() / (wsReal)SQR(N_samp-N_tsz);

	return R_cur;
}

#define SATISFY(orig, nw) ((orig) > (nw))

wsReal cStudyAnalysis::_exhaustive(wsUint **Na_Lcomb, wsUint *Np_comb, wsUint N_combReq,
	wsUint N_samp, wsMat Ra_cor)
{
	cComb X;
	wsUint *Na_comb = NULL;
	wsAlloc(Na_comb, wsUint, 1000*N_combReq);
	wsAlloc(*Na_Lcomb, wsUint, N_combReq);
	wsReal R_exLeastVar = WISARD_NAN;
	X.init(2, N_samp);
	wsUint N_rcomb = X.get(1000, Na_comb);
	while (N_rcomb) {
		wsUint *Na_pcomb = Na_comb;
		for (wsUint i=0 ; i<1000 ; i++,Na_pcomb+=N_combReq) {
			char *Ba_cinc = NULL;
			wsCalloc(Ba_cinc, char, N_samp);
			for (wsUint j=0 ; j<N_combReq ; i++)
				Ba_cinc[Na_pcomb[j]] = 1;

			wsReal R_curV = CONDVAR(getIO(), Ra_cor, N_samp, Ba_cinc, N_combReq);

			if (NA(R_exLeastVar) || SATISFY(R_exLeastVar, R_curV)) {
				memcpy(Na_Lcomb, Na_pcomb, sizeof(wsUint)*N_combReq);
				R_exLeastVar = R_curV;
// 				LOG("Comb %d,%d var [%g]\n", Na_pcomb[0], Na_pcomb[1],
// 					R_exLeastVar);
			}
		}
		N_rcomb = X.get(1000, Na_comb);
	}

	/* Copy combination */
	return R_exLeastVar;
}

wsReal cStudyAnalysis::_exp(wsUint N_szComb, wsUint N_samp,
	wsSym Ra_cor, xRealSort *Xa_sort, wsReal *Rp_expSS)
{
	vInt	Xv_finSels;
//	wsReal	R_finVar	= REAL_CONST(999999999);
	//wsUint	i			= 0;
	// 	Xv_sels.push_back(Xa_sort[0].i);
	// 	Ba_isInc[Xa_sort[0].i] = 1;
	// 	for (wsUint i=1 ; i<N_samp && Xa_sort[i].V == R_v ; i++) {
	// 		Xv_sels.push_back(Xa_sort[i].i);
	// 		Ba_isInc[Xa_sort[i].i] = 1;
	// 		N_tsz++;
	// 	}
	wsReal *Ra_v = sseVector(1000);
	wsReal R_sumSS = W0;
	for (wsUint z=0 ; z<1000 ; z++) {
		vInt Xv_sels;
		char	*Ba_isInc	= NULL;
		char	*Ba_prhb	= NULL;
		wsCalloc(Ba_prhb, char, N_samp);
		wsCalloc(Ba_isInc, char, N_samp);
		for (wsUint i=0,I=0 ; I<N_szComb ; i++) {
			wsUint N_idx = Xa_sort[rand()%N_samp].i;
			if (Ba_prhb[N_idx] || Ba_isInc[N_idx]) continue;

			Xv_sels.push_back(N_idx);
			Ba_isInc[N_idx] = 1;
			I++;
		}
		wsReal R_curSS = W0;
		Ra_v[z] = CONDVAR(getIO(), Ra_cor, N_samp, Ba_isInc, N_szComb, &R_curSS);
		R_sumSS += R_curSS;
		notice("Loop %d/%d...\r", z, 1000);
		DEALLOC(Ba_prhb);
		DEALLOC(Ba_isInc);
	}

	if (OPT_ENABLED(verbose)) {
		exportVector("exp.var.res", Ra_v, 1000);
	}

	/* Determine final set! */
	wsReal R_avg = sseVsum(Ra_v, 1000)/REAL_CONST(1000.0);
	*Rp_expSS = R_sumSS / REAL_CONST(1000.0);
	return R_avg;
}

wsUint* true2idx(char *Ba_set, wsUint N_samp, wsUint N_expect)
{
	wsUint	i=0;
	wsUint*	Na_nset = NULL;
	wsAlloc(Na_nset, wsUint, N_expect);

	for (wsUint r=0 ; r<N_samp ; r++) {
		if (Ba_set[r]) {
			if (i == N_expect) halt("Expectation[%d] overflow!", N_expect);
			Na_nset[i++] = r;
		}
	}
	if (i != N_expect)
		halt("Expectation[%d] underflow[%d]!", N_expect, i);

	return Na_nset;
}

wsReal cStudyAnalysis::_stepwise(wsUint **Np_comb, wsUint *Np_retSize, wsUint N_samp2rem,
	wsUint N_samp, wsSym Ra_cor, xRealSort *Xa_sort)
{
	vInt			Xv_finSels;
	wsReal			R_finLSV	= REAL_CONST(999999999.0);
	int				N_pm		= 2;
	char			B_findExact	= 0;
	vector<vInt>	Xv_candSels;
	vReal			Xv_candLSV;

	/* N_wantSize == # of samples to be REMOVED */
	LOG("[%d] samples required to chosen in step-wise to removed\n", N_samp2rem);

	/* For each sample, do step-wise */
	for (wsUint z=0 ; z<N_samp ; z++) {
		vInt Xv_sels;
		char	*Ba_rem		= NULL; /* 1=Samp. to be removed, 0 otherwise */
		char	*Ba_prhb	= NULL;
		wsUint	N_sel2rem	= 1;
		wsCalloc(Ba_prhb, char, N_samp);
		wsCalloc(Ba_rem, char, N_samp);

		lverbose("Perform step-wise for sample [%d]\n", z);

		/* Init selection buffer */
		Xv_sels.push_back(Xa_sort[z].i);

		/* Anyway select this sample at first */
		Ba_rem[Xa_sort[z].i] = 1;
		wsReal	R_leastVar	= CONDVAR(getIO(), Ra_cor, N_samp, Ba_rem, N_sel2rem);
		wsReal	R_initVar	= R_leastVar;

		/* Add other babies and get conditional variance again */
		for (wsUint N_step=0 ; N_samp > Xv_sels.size() ; N_step++) {
			lverbose("Step %d : Set size [%d], variance [%g]\n", N_step, Xv_sels.size(), R_leastVar);
			if (OPT_ENABLED(verbose)) {
				LOG("Set :");
				FOREACH (vInt_it, Xv_sels, x)
					LOGnf(" %d", *x);
				LOGnf("\n");
			}
			wsReal	R_curLSV	= WISARD_NAN;
			int		N_idx		= -1;

			/* Update selected sample size */
			N_sel2rem = (wsUint)Xv_sels.size();
			for (wsUint i=0 ; i<N_samp ; i++) {
				/* Pass if this sample is already chosen or should be ignored */
				if (Ba_rem[i] || Ba_prhb[i]) continue;

				/* Copy and insert new one */
				char *Ba_newIncl = dupSel(Ba_rem, N_samp);
				Ba_newIncl[i] = 1;

				/* Get conditional variance */
				wsReal R_curVar = CONDVAR(getIO(), Ra_cor, N_samp, Ba_newIncl, N_sel2rem + 1);
				if (N_idx == -1 || SATISFY(R_curLSV, R_curVar)) {
					N_idx = i;
					R_curLSV = R_curVar;
				}
			}

			/* Failed to increase set size */
			if (!SATISFY(R_leastVar, R_curLSV)) {
				/* Backward */
				N_idx		= -1;
				R_curLSV	= WISARD_NAN;

				/* Remove one of them to reduce */
				if (N_sel2rem > 1) for (wsUint i=0 ; i<N_sel2rem ; i++) {
					/* Copy and remove one */
					char *Ba_isNewInc = dupSel(Ba_rem, N_samp);
					Ba_isNewInc[Xv_sels[i]] = 0;

					/* Get conditional variance */
					wsReal R_curVar = CONDVAR(getIO(), Ra_cor, N_samp, Ba_isNewInc, N_sel2rem - 1);
					if (N_idx == -1 || SATISFY(R_curLSV, R_curVar)) {
						N_idx		= Xv_sels[i];
						R_curLSV	= R_curVar;
					}
				}

				/* Also, failed to reduce, then stop */
				if (N_idx == -1 || !SATISFY(R_leastVar, R_curLSV))
					break;

				/* Else, remove it from set and mark as it does not need to be added */
				FOREACH (vInt_it, Xv_sels, j) if (*j == N_idx) {
					Xv_sels.erase(j);
					Ba_rem[N_idx]	= 0;
					R_leastVar		= R_curLSV;
					Ba_prhb[N_idx]	= 1;
					break;
				}
				continue;
			}

			/* Add item if success to increase */
			R_leastVar		= R_curLSV;
			Ba_rem[N_idx]	= 1;
			Xv_sels.push_back(N_idx);
		}
		/* Uninitialize */
		DEALLOC(Ba_prhb);
		DEALLOC(Ba_rem);

		/* Current size */
		wsUint N_cur = (wsUint)Xv_sels.size();

		/* Insert 'candidate set' if the result set within WANT +- N_pm */
		int N_diff = abs((int)N_samp2rem - (int)N_cur);
		if (N_diff <= N_pm) {
			pverbose("Requested sample size [%d], selected with [%d] samples\n",
				N_samp2rem, N_cur);
			Xv_candSels.push_back(Xv_sels);
			Xv_candLSV.push_back(R_leastVar);
		}

		/* Replace final set if required */
		if (SATISFY(R_finLSV, R_leastVar) && (wsUint)Xv_sels.size() <= N_samp2rem) {
			Xv_finSels.resize(Xv_sels.size());
			copy(Xv_sels.begin(), Xv_sels.end(), Xv_finSels.begin());
			R_finLSV = R_leastVar;
			LOG("Update : Set size [%d], variance [%g->%g]\n", Xv_sels.size(),
				R_initVar, R_leastVar);
		} else
			lverbose("Loop : Set size [%d], variance [%g->%g]\n", Xv_sels.size(),
				R_initVar, R_leastVar);

		if (OPT_ENABLED(verbose)) {
			LOG("Set :");
			FOREACH (vInt_it, Xv_sels, x)
				LOGnf(" %d", *x);
			LOGnf("\n");
		}
	}
	if ((wsUint)Xv_finSels.size() == N_samp2rem) {
		B_findExact = 1;
		LOG("Exact size[%d] with variance [%g] found\n", N_samp2rem,
			R_finLSV);
	} else LOG("Inexact size[%d] with variance [%g] found\n", Xv_finSels.size(),
		R_finLSV);

	/* Post-finding */
	wsReal*	Ra_finVars = sseVector((wsUint)Xv_candSels.size());
	wsUint	m		= 0;
	if (Xv_candSels.size())
		LOG("Post-finding process is now applied to [%d] sets\n", Xv_candSels.size());
	else
		LOG("No post-finding process is possible\n");
	/* For all candidate sets */
	FOREACHDO (vector<vInt>::iterator, Xv_candSels, i, m++) {
		vInt&	Xv_set		= *i;
		int		N_ssz		= (int)Xv_set.size();
		wsReal	R_mVar		= Xv_candLSV[m];
		wsUint	N_upd		= 0;

		/* New set */
		wsUint*	Na_nset		= NULL;
		
		/* Added or removed? */
		int		N_supp		= (int)N_samp2rem - N_ssz;
		wsUint	N_absSupp	= N_supp<0?-N_supp:N_supp;

		/* Mark selected samples */
		char*	Ba_sel		= NULL;
		wsCalloc(Ba_sel, char, N_samp);
		FOREACH (vInt_it, Xv_set, j)
			Ba_sel[*j] = 1;
		if (OPT_ENABLED(verbose)) {
			LOG("Post-upd set :");
			FOREACH (vInt_it, Xv_set, x)
				LOGnf(" %d", *x);
			LOGnf("\n");
		}

		/* Should be added */
		if (N_supp > 0) {
			/* Make unselected set */
			wsUint *Na_usel	= NULL, I=0;
			wsUint N_usel	= N_samp - N_ssz;
			wsAlloc(Na_usel, wsUint, N_usel);
			for (wsUint i=0 ; i<N_samp ; i++)
				if (!Ba_sel[i]) Na_usel[I++] = i;
			if (I != N_usel) halt("STOP");

			cComb X(N_absSupp, N_usel);
			wsUint N_gen, *Na_comb = NULL;
			wsAlloc(Na_comb, wsUint, 1000*N_absSupp);

			/* Add and test */
			for ( ; (N_gen=X.get(1000, Na_comb)) ; ) {
				for (wsUint q=0,*Np_comb=Na_comb ; q<N_gen ; q++,Np_comb+=N_absSupp) {
					char*	Ba_nsel	= dupSel(Ba_sel, N_samp);

					/* Add */
					for (wsUint x=0 ; x<N_absSupp ; x++)
						Ba_nsel[Na_usel[Np_comb[x]]] = 1;

					/* Add variable and get variance */
					wsReal R_vv = CONDVAR(getIO(), Ra_cor, N_samp, Ba_nsel, N_samp2rem);
					//pverbose("Set [%d], comb [%g]\n", m, R_vv);
					if (!SATISFY(R_mVar, R_vv)) continue;

					/* Copy */
					DEALLOC(Na_nset);
					Na_nset = true2idx(Ba_nsel, N_samp, N_samp2rem);
					DEALLOC(Ba_nsel);
					R_mVar = R_vv;
					N_upd++;
				}
			}

			DEALLOC(Na_usel);
		} else if (N_supp < 0) {
			/* Make selected set */
			wsUint *Na_sel	= NULL, I=0;
			wsUint N_sel	= N_ssz;
			wsAlloc(Na_sel, wsUint, N_sel);
			FOREACHDO (vInt_it, Xv_set, j, I++)
				Na_sel[I] = *j;

			cComb X(N_absSupp, N_sel);
			wsUint N_gen, *Na_comb = NULL;
			wsAlloc(Na_comb, wsUint, 1000*N_absSupp);

			/* Remove and test */
			for ( ; (N_gen=X.get(1000, Na_comb)) ; ) {
				for (wsUint q=0,*Np_comb=Na_comb ; q<N_gen ; q++,Np_comb+=N_absSupp) {
					char*	Ba_nsel	= dupSel(Ba_sel, N_samp);

					/* Remove */
					for (wsUint x=0 ; x<N_absSupp ; x++)
						Ba_nsel[Na_sel[Np_comb[x]]] = 0;

					/* Add variable and get variance */
					wsReal R_vv = CONDVAR(getIO(), Ra_cor, N_samp, Ba_nsel, N_samp2rem);
					//pverbose("Set [%d], comb [%g]\n", m, R_vv);
					if (!SATISFY(R_mVar, R_vv)) continue;

					/* Copy */
					DEALLOC(Na_nset);
					Na_nset = true2idx(Ba_nsel, N_samp, N_samp2rem);
					DEALLOC(Ba_nsel);
					R_mVar = R_vv;
					N_upd++;
				}
			}

			DEALLOC(Na_sel);
		}

		DEALLOC(Ba_sel);

		/* Get complete set */
		Ra_finVars[m] = R_mVar;
		pverbose("Post-go : Set size [%d], variance [%g->%g]\n", N_samp2rem,
			Xv_candLSV[m], R_mVar);
		/* Skip if no update */
		if (N_upd == 0) continue;

		/* Replace final set if required */
		if (SATISFY(R_finLSV, R_mVar) || B_findExact == 0) {
			Xv_finSels.clear();
			for (wsUint x=0 ; x<N_samp2rem ; x++)
				Xv_finSels.push_back(Na_nset[x]);
			LOG("Post-upd : Set size [%d], variance [%g->%g]\n", N_samp2rem,
				R_finLSV, R_mVar);
			R_finLSV = R_mVar;
			B_findExact = 1;
		}
		DEALLOC(Na_nset);
	}
	if (Xv_candSels.size()) {
		exportVector("var.res", Ra_finVars, (wsUint)Xv_candSels.size());
	}
	/* If there is no exact, do Monte-carlo */
	if (B_findExact == 0) {
		LOG("Cannot found exact size of exclusion, Monte-carlo start\n");

		/* Build permutation buffer with exact length */
		wsUint *Na_fin = NULL;
		wsAlloc(Na_fin, wsUint, N_samp2rem);
		/* Mark selected samples */
		char*	Ba_sel		= NULL;
		wsCalloc(Ba_sel, char, N_samp);

		/* Make base-set */
		wsUint K=0;
		FOREACHDO (vInt_it, Xv_finSels, i, K++) {
			Ba_sel[*i] = 1;
			Na_fin[K] = *i;
		}

		/* Do permutation */
		wsReal R_lgst = WISARD_NAN;

		for (wsUint i=0 ; i<1000 ; i++) {
			char *Ba_nsel = dupSel(Ba_sel, N_samp);

			for (wsUint j=K ; j<N_samp2rem ; j++) {
				wsUint Q = 0;
				do {
					Q = rand()%N_samp;
				} while (Ba_nsel[Q]);
				Ba_nsel[Q] = 1;
				Na_fin[j] = Q;
			}

			/* Add variable and get variance */
			wsReal R_vv = CONDVAR(getIO(), Ra_cor, N_samp, Ba_nsel, N_samp2rem);
			if (NA(R_lgst) || SATISFY(R_lgst, R_vv)) {
				Xv_finSels.clear();
				for (wsUint x=0 ; x<N_samp2rem ; x++)
					Xv_finSels.push_back(Na_fin[x]);
				LOG("Monte-upd : Set size [%d], variance [%g->%g]\n", N_samp2rem,
					R_lgst, R_vv);
				R_lgst = R_vv;
			}

			DEALLOC(Ba_nsel);
		}

		DEALLOC(Ba_sel);
//		LOG("Exact same numbers of samples couldn't chosen\n");
		R_finLSV = R_lgst;
	}

	/* Allocate and copy to return */
	wsAlloc(*Np_comb, wsUint, Xv_finSels.size());
	for (size_t i=0 ; i<Xv_finSels.size() ; i++)
		(*Np_comb)[i] = Xv_finSels[i];
	*Np_retSize	= (wsUint)Xv_finSels.size();

	return R_finLSV;
}

void cStudyAnalysis::run()
{
	wsMat		Ra_cor		= getFullCorMat(Cp_anaCor);
	vSampPtr	&Xa_samp	= Cp_IO->getSample();
	wsUint		N_samp		= Cp_IO->sizeSample();
	xRealSort*	Xa_sort		= NULL;

	/*
	 *
	 * Preparing
	 * 
	 */

	/* Init result vector */
	V_varSamp.init(N_samp);
	wsReal*		Ra_vs = V_varSamp.get();
	cSymMatrix	M_cor(Ra_cor, N_samp, 1);

	/* Compute conditional variance for each sample */
	cTimer	H_timeCondVar;
	H_timeCondVar.start();
	Ra_vs	= getCondVar(getIO(), Xa_samp, M_cor);
	LOG("Conditional variance computation took [%s]\n", H_timeCondVar.getReadable());

	/* Export result *as is* */ {
		cTableExporter	C_vs("sampvar.res", "ssr",
			"Estimation of conditional variance per sample", 0,
			3, "FID", "IID", "CONDVAR");
		C_vs.sampleWise(getIO(), Ra_vs);
	}

	/* Sort */
	Xa_sort	= buildRealSort(Ra_vs, N_samp);
	qsort(Xa_sort, N_samp, sizeof(xRealSort), sort_real);

	/* Print top K */
	wsUint N_prt = IS_ASSIGNED(nsamp) ? N_samp-(wsUint)(OPT_RANGE(nsamp).R_s) :
		N_samp;
	LOGoutput("Top [%d] conditional variance per sample is exported " /* CONFIRMED 140109 */
		"to [%s.top.sampvar.res]\n", N_prt, OPT_STRING(out));
	cExporter	*C_pt	= cExporter::summon("top.sampvar.res");
	cExporter	*C_pt2	= cExporter::summon("top.sampvar.lst");
	C_pt->put("RANK	FID	IID	CONDVAR\n");
	for (wsUint i=0 ; i<N_prt ; i++) {
		wsUint j = Xa_sort[i].i;
		C_pt->fmt("%d	%s	%s	%g\n", i+1, Xa_samp[j]->S_FID.c_str(),
			Xa_samp[j]->S_IID.c_str(), Xa_sort[i].V);
		C_pt2->fmt("%s	%s\n", Xa_samp[j]->S_FID.c_str(),
			Xa_samp[j]->S_IID.c_str());
	}
	for (wsUint i=N_prt ; i<N_samp ; i++) {
		wsUint j = Xa_sort[i].i;
		C_pt->fmt("%d	%s	%s	%g\n", i+1, Xa_samp[j]->S_FID.c_str(),
			Xa_samp[j]->S_IID.c_str(), Xa_sort[i].V);
	}
	delete C_pt;
	delete C_pt2;

	wsReal	R_leastVar	= WISARD_NAN;	///< Least conditional variance (LSV)
	wsUint*	Na_comb		= NULL;			///< Sample indices giving LSV
	wsUint	N_szComb	= 0;			///< Size of Na_comb
	wsUint	N_samp2sel	= (wsUint)(OPT_RANGE(nsamp).R_s);

	/* Exhaustive method */
//	R_leastVar = _exhaustive(&Ra_comb, &N_szComb, N_combReq, N_samp, Ra_cor);

	/* Step-wise method
	 * 1 : Choose top-ranked babies */
	LOG("[STEP 2] Selecting samples to be removed with step-wise method...\n");
	wsUint	N_samp2rem	= N_samp - N_samp2sel;
	/* Via step-wise, should select N_samp2rem samples */
	R_leastVar = _stepwise(&Na_comb, &N_szComb, N_samp2rem, N_samp,
		Ra_cor, Xa_sort);
	wsReal R_lXLeast = W1 - R_leastVar/(wsReal)N_samp2rem;
	/* Print top K */
	LOGoutput("Top [%d] conditional variance by step-wise is exported "
		"to [%s.top.sampvar2.lst]\n", N_szComb, OPT_STRING(out));
	C_pt2	= cExporter::summon("top.sampvar2.lst");
	for (wsUint i=0 ; i<N_szComb ; i++)
		C_pt2->fmt("%s	%s\n", Xa_samp[Na_comb[i]]->S_FID.c_str(),
			Xa_samp[Na_comb[i]]->S_IID.c_str());
	delete C_pt2;

	/* Getting expectation */
	LOG("[STEP 3] Computing expectation...\n");
	wsReal R_expExclSS	= WISARD_NAN, R_expInclSS = WISARD_NAN;
	/* Expected cond. var. with # samples to be removed */
	wsReal R_expExcl	= _exp(N_szComb, N_samp, Ra_cor, Xa_sort, &R_expExclSS);
	wsReal R_lXExcl		= W1 - R_expExcl/(wsReal)N_szComb;
	/* Expected cond. var. with # samples to be chosen */
	wsReal R_expIncl	= _exp(N_samp2sel, N_samp, Ra_cor, Xa_sort, &R_expInclSS);
	wsReal R_lXIncl		= W1 - R_expIncl/(wsReal)N_samp2sel;
	LOG("Expected-exclude : Set size [%d], variance [%g], Lambda(X) [%g]\n", N_szComb,
		R_expExcl, R_lXExcl);
	LOG("Expected-include : Set size [%d], variance [%g], Lambda(X) [%g]\n", N_szComb,
		R_expIncl, R_lXIncl);

	/* Get the effective sample size */
	LOG("[STEP 4] Getting effective sample size...\n");
	char *Na_excl = NULL;
	wsCalloc(Na_excl, char, N_samp);
	for (wsUint i=0 ; i<N_szComb ; i++)
		Na_excl[Na_comb[i]] = 1;
	/* Actual sample size = entire - selected */
	wsUint	N_actualSel	= N_samp - N_szComb;
	wsMat	Ra_corSel	= sseMsubsetRect(Ra_cor, N_samp, Na_excl, 0, N_actualSel);
	wsReal	R_effSmpSz	= getEffSmpSz(Ra_corSel, N_actualSel);
	sseUnmat(Ra_corSel, N_actualSel);
	wsMat	Ra_corNsel	= sseMsubsetRect(Ra_cor, N_samp, Na_excl, 1, N_szComb);
	wsReal	R_NeffSmpSz	= getEffSmpSz(Ra_corNsel, N_szComb);
	sseUnmat(Ra_corNsel, N_szComb);
	/* Entire variance */
	wsReal	R_AeffSmpSz	= getEffSmpSz(Ra_cor, N_samp);
	
	/* Export effective sample size */
	cExporter *C_pt3 = cExporter::summon("sampsel.res");
	C_pt3->fmt("N_SELSAMP	%d\n", N_actualSel);
	C_pt3->put("IID_SELSAMP	");
	for (wsUint i=0,j=0 ; i<N_samp ; i++) {
		if (Na_excl[i]) continue;
		if (!j) C_pt3->fmt("%s", Xa_samp[i]->S_IID.c_str());
		else C_pt3->fmt(",%s", Xa_samp[i]->S_IID.c_str());
		j++;
	}
	C_pt3->put("\n");
	C_pt3->fmt("VAR_REMSAMP	%g\n", R_leastVar);
	C_pt3->fmt("LAMBDAX_REMSAMP	%g\n", R_lXLeast);
	C_pt3->fmt("ALL_EFFSAMP	%g\n", R_AeffSmpSz);
	C_pt3->fmt("SEL_EFFSAMP	%g\n", R_effSmpSz);
	C_pt3->fmt("NSEL_EFFSAMP	%g\n", R_NeffSmpSz);
	C_pt3->fmt("EXP_EXCL	%g\n", R_expExcl);
	C_pt3->fmt("EFFSAMP_EXCL	%g\n", R_expExclSS);
	C_pt3->fmt("LAMBDAX_EXCL	%g\n", R_lXExcl);
	C_pt3->fmt("EXP_INCL	%g\n", R_expIncl);
	C_pt3->fmt("EFFSAMP_INCL	%g\n", R_expInclSS);
	C_pt3->fmt("LAMBDAX_INCL	%g\n", R_lXIncl);

	/* Determine final set! */
	LOG("Final : Removal set size [%d], variance [%g], Lambda(X) [%g]\n",
		N_szComb, R_leastVar, R_lXLeast);
	LOG("        Eff. samp. size of entire   [%g]\n", R_AeffSmpSz);
	LOG("        Eff. samp. size of selected [%g]\n", R_effSmpSz);
}

#endif

} // End namespace ONETOOL
