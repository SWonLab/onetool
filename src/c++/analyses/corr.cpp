#include "global/common.h"
#include "analyses/corr.h"
#include "analyses/ppp.h"
#include "utils/matrix.h"
#include "input/stream.h"
#include "analyses/pddt.h"
#include "app/app_wisard.h"

namespace ONETOOL {

#if TOOLSET_TYPE & FUNCTION_ONETOOL

extern cFamStrAnalysis*	Cp_anaFamStr;
extern cPPPAnalysisV2*		Cp_PPP;
extern cAnalysis*			Cp_corr;

cAnalysis*	getCorrelation(cIO *io)
{
	cTimer t2;

	if (Cp_corr == NULL) {
		if (IS_ASSIGNED(cor)) {
			if (OPT_ENABLED(founderonly))
				halt("--cor option cannot be used when --founderonly is assigned");
			if (OPT_ENABLED(kinship))
				LOG("--kinship detected, input correlation will be considered as kinship coefficient\n");
			if (OPT_ENABLED(ibs))
				halt("--cor option cannot be used when --ibs is assigned");
			if (OPT_ENABLED(indep))
				halt("--cor option cannot be used when --indep is assigned");
			if (OPT_ENABLED(ktau))
				halt("--cor option cannot be used when --ktau is assigned");
			if (OPT_ENABLED(empktau))
				halt("--cor option cannot be used when --empktau is assigned");

			getPPP(io);

			t2.start();
			Cp_corr = new cCorrAnalysis(io, OPT_STRING(cor));
			LOG("[%s] Pre-calculated correlation loaded\n", t2.getReadable());
		} else if (OPT_ENABLED(kinship) || OPT_ENABLED(hybrid)) {
			if (IS_ASSIGNED(cor))		halt("--cor option cannot be used when --kinship is assigned");
			if (IS_ASSIGNED(ibs))		halt("--ibs option cannot be used when --kinship is assigned");
			if (IS_ASSIGNED(hamming))	halt("--hamming option cannot be used when --kinship is assigned");
			if (IS_ASSIGNED(bn))		halt("--bn option cannot be used when --kinship is assigned");
			if (IS_ASSIGNED(indep))		halt("--indep option cannot be used when --kinship is assigned");

			/* If --hybrid is on */
			/* Calc PPP */
			if (OPT_ENABLED(hybrid))
				getPPP(io);
			t2.start();
			Cp_corr = new cPDDTAnalysis(io);
			Cp_corr->run();
			LOG("[%s] Kinship coefficient calculated\n", t2.getReadable());
		} else {
			cPPPAnalysisV2 *Cp_anaPPP = getPPP(io);

			t2.start();
			Cp_corr = new cCorrAnalysis(io, Cp_anaPPP);
			Cp_corr->run();
			LOG("[%s] Correlation calculated\n", t2.getReadable());
		}
	}

	return Cp_corr;
}

#endif

void getX2coeff(map<string,xPDDTnode> &X_map, vector<wsReal> &Xa_scores,
				xSample *Xp_src, xSample *Xp_cur, xSample *Xp_tgt, int B_visitingTwin,
				int B_visitingFemTopNode, wsReal R_score,
				char B_isX);

void sym_sseMpMt(wsFloat **Rp_m1, wsUint N_r1, wsUint N_c1,
	wsReal **Rp_save)
{
	wsUint		i, j, k;

	/* For all row */
	for (i=0 ; i<N_r1 ; i++) {
		if ((i%10) == 0)
			notice("[%d/%d] samples processed...\r", i, N_r1);
		wsReal	*Rp_matR	= Rp_save[i];
		wsFloat	*Rp_1		= Rp_m1[i];

		for (j=0 ; j<=i ; j++) {
			wsFloat *Rp_2	= Rp_m1[j];
			wsUint N_med	= 0;

#ifdef USE_SSE
			N_med = corrGetMed(N_c1);
			if (N_med != 0) {
				corrsse_t sse_R = corrsseSet(0.0);
				for (k=0 ; k<N_med ; k+=corrsseJmp) {
					corrsse_t *sse_m1 = (corrsse_t *)(Rp_1 + k);
					corrsse_t *sse_m2 = (corrsse_t *)(Rp_2 + k);

					sse_R = corrsseAdd(sse_R, corrsseMul(*sse_m1, *sse_m2));
				}
				corrsseSum(sse_R, Rp_matR[j]);
			}
#endif
			for (k=N_med ; k<N_c1 ; k++)
				Rp_matR[j] += Rp_1[k] * Rp_2[k];
		}
	}
}

int thr_count(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	if (N_idx < 0 || N_idx >= OPT_NUMBER(thread))
		halt("Error index [%d]", N_idx);
	xCorrThread	*Xp_at		= (xCorrThread *)Vp_shareData;
	int*		Np_data		= (int *)Vp_data;
	wsMat		Ra_corr		= Xp_at->Ra_corr;
//	wsFloat**	Ra_data		= Xp_at->Ra_data;
	wsUint		N_var		= Xp_at->N_variant;

	wsUint		N_mBlk		= Xp_at->N_mBlk;
	WISARD_pc**	Na_mMissing	= Xp_at->Na_mMissing;
	wsUint*		Na_cachedN	= Xp_at->Na_cachedN;
	wsUint		N_noIncCor	= Xp_at->N_noIncCor;
	int**		Na_corrN	= Xp_at->Na_corrN;

	wsUint		N_s			= (wsUint)Np_data[0];
	wsUint		N_e			= (wsUint)Np_data[1];
	wsUint		N_jmp		= OPT_NUMBER(thread);
//	wsUint		N_med		= 0;
	wsUint		i, j, k;
	pverbose("Thread [%d] computes from [%d] to [%d] with interval [%d]\n",
		N_idx, N_s, N_e, N_jmp);

//#ifdef USE_SSE
//	N_med = corrGetMed(N_var);
//#endif

	for (i=N_s ; i<=N_e ; i+=N_jmp) {
		if ((i%10) == 0) notice("[%d/%d] samples processed...\r", i, Xp_at->N_sample);
		for (j=0 ; j<=i ; j++) {
			WISARD_pc N = 0;

			if (Na_mMissing[i]) {
				if (Na_mMissing[j]) for (k=0 ; k<N_mBlk ; k++)
					/* Noncomplete - Noncomplete case */
						N += WISARD_POPCNT(Na_mMissing[i][k] & Na_mMissing[j][k]);
				else {
					/* Noncomplete - Complete case */
					N = Na_cachedN[i];
				}
			} else {
				if (Na_mMissing[j]) {
					/* Complete - Noncomplete case */
					N = Na_cachedN[j];
				} else
					/* Complete - Complete case */
					N = N_var-N_noIncCor;
			}

			Ra_corr[i][j] /= (wsReal)N;
			Ra_corr[j][i] = Ra_corr[i][j];
			Na_corrN[i][j] = (int)N;
		}
	}

	return 0;
}

int thr_sym_sseMpMt(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	if (N_idx < 0 || N_idx >= OPT_NUMBER(thread))
		halt("Error index [%d]", N_idx);
	xCorrThread		*Xp_at		= (xCorrThread *)Vp_shareData;
	int				*Np_data	= (int *)Vp_data;
	wsMat			Ra_corr		= Xp_at->Ra_corr;
	wsFloat**		Ra_data		= Xp_at->Ra_data;
	wsUint			N_var		= Xp_at->N_variant;

	wsUint			N_s			= (wsUint)Np_data[0];
	wsUint			N_e			= (wsUint)Np_data[1];
	wsUint			N_jmp		= OPT_NUMBER(thread);
	wsUint			N_med		= 0;
	wsUint			i, j, k;
	pverbose("Thread [%d] computes from [%d] to [%d] with interval [%d]\n",
		N_idx, N_s, N_e, N_jmp);

#ifdef USE_SSE
	N_med = corrGetMed(N_var);
#endif

	/* For all row */
	for (i=N_s ; i<=N_e ; i+=N_jmp) {
		if ((i%50) == 0) notice("[%d/%d] samples processed...\r", i, Xp_at->N_sample);
		wsReal	*Rp_matR	= Ra_corr[i];
		wsFloat	*Rp_1		= Ra_data[i];

		for (j=0 ; j<=i ; j++) {
			wsFloat *Rp_2	= Ra_data[j];

#ifdef USE_SSE
			if (N_med != 0) {
				corrsse_t sse_R = corrsseSet(0.0);
				for (k=0 ; k<N_med ; k+=corrsseJmp) {
					corrsse_t *sse_m1 = (corrsse_t *)(Rp_1 + k);
					corrsse_t *sse_m2 = (corrsse_t *)(Rp_2 + k);

					sse_R = corrsseAdd(sse_R, corrsseMul(*sse_m1, *sse_m2));
				}
				corrsseSum(sse_R, Rp_matR[j]);
			}
#endif
			for (k=N_med ; k<N_var ; k++)
				Rp_matR[j] += Rp_1[k] * Rp_2[k];
		}
	}

	return 0;
}

#define ELEM_SWAP(a,b,x) { register x t=(a);(a)=(b);(b)=t; }  

wsFloat quickselect(wsFloat *arr, int nelems, int select)  
{  
	int low, high, middle, ll, hh;  

	low = 0;  
	high = nelems - 1;  

	for (;;) {  
		if (high <= low) /* One element only */  
			return arr[select];  

		if (high == low + 1) {  /* Two elements only */  
			if (arr[low] > arr[high])  
				ELEM_SWAP(arr[low], arr[high], wsFloat);
			return arr[select];  
		}  

		/* Find median of low, middle and high items; swap into position low */  
		middle = (low + high) / 2;  
		if (arr[middle] > arr[high]) ELEM_SWAP(arr[middle], arr[high], wsFloat);
		if (arr[low]	> arr[high]) ELEM_SWAP(arr[low],	arr[high], wsFloat);
		if (arr[middle] > arr[low])	 ELEM_SWAP(arr[middle], arr[low] , wsFloat);

		/* Swap low item (now in position middle) into position (low+1) */  
		ELEM_SWAP(arr[middle], arr[low+1], wsFloat);

		/* Nibble from each end towards middle, swapping items when stuck */  
		ll = low + 1;  
		hh = high;  
		for (;;) {  
			do ll++; while (arr[low] > arr[ll]);  
			do hh--; while (arr[hh]  > arr[low]);  

			if (hh < ll)  
				break;  

			ELEM_SWAP(arr[ll], arr[hh], wsFloat);
		}  

		/* Swap middle item (in position low) back into correct position */  
		ELEM_SWAP(arr[low], arr[hh], wsFloat);

		/* Re-set active partition */  
		if (hh <= select)  
			low = ll;  
		if (hh >= select)  
			high = hh - 1;  
	}  
}  

wsFloat quickmedian(wsFloat *arr, int nelems)  
{  
	return quickselect(arr, nelems, (nelems - 1) / 2);  
}  

wsReal quickselect(wsReal *arr, int nelems, int select)  
{  
	int low, high, middle, ll, hh;  

	low = 0;  
	high = nelems - 1;  

	for (;;) {  
		if (high <= low) /* One element only */  
			return arr[select];  

		if (high == low + 1) {  /* Two elements only */  
			if (arr[low] > arr[high])  
				ELEM_SWAP(arr[low], arr[high], wsReal);
			return arr[select];  
		}  

		/* Find median of low, middle and high items; swap into position low */  
		middle = (low + high) / 2;  
		if (arr[middle] > arr[high]) ELEM_SWAP(arr[middle], arr[high], wsReal);
		if (arr[low]	> arr[high]) ELEM_SWAP(arr[low],	arr[high], wsReal);
		if (arr[middle] > arr[low])	 ELEM_SWAP(arr[middle], arr[low] , wsReal);

		/* Swap low item (now in position middle) into position (low+1) */  
		ELEM_SWAP(arr[middle], arr[low+1], wsReal);

		/* Nibble from each end towards middle, swapping items when stuck */  
		ll = low + 1;  
		hh = high;  
		for (;;) {  
			do ll++; while (arr[low] > arr[ll]);  
			do hh--; while (arr[hh]  > arr[low]);  

			if (hh < ll)  
				break;  

			ELEM_SWAP(arr[ll], arr[hh], wsReal);
		}  

		/* Swap middle item (in position low) back into correct position */  
		ELEM_SWAP(arr[low], arr[hh], wsReal);

		/* Re-set active partition */  
		if (hh <= select)  
			low = ll;  
		if (hh >= select)  
			high = hh - 1;  
	}  
}  

wsReal quickmedian(wsReal *arr, int nelems)  
{  
	return quickselect(arr, nelems, (nelems - 1) / 2);  
}  

#undef ELEM_SWAP

int forAllSample_corr(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx = (int *)Vp_thrData;
	static wsUint N_idxSample;
	wsUint N_sample = ((xCorrThread *)Vp_shareData)->N_sample;
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		N_idxSample = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++) {
			if ((N_idxSample+i) >= N_sample)
				break;
			Np_idx[i] = N_idxSample+i;
			if (((N_idxSample+i)%10) == 0)
				notice("%d/%d samples processed\r", N_idxSample+i, N_sample);
		}
		N_idxSample += N_thread;
		break;
	case DISTFUN_UNINIT:
		LOG("%d/%d samples processed\n", N_sample, N_sample);
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt_fmt(WISARD_SYST_INVL_DISTFUNC_CMD, N_mode);
//		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}
	return i;
}

int forAllSample_corr_equal(int N_mode, int N_thread, void *Vp_data,
					  void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx = (int *)Vp_data;
	wsUint	N_sample = ((xCorrThread *)Vp_shareData)->N_sample;
	static wsUint
			N_proc	= 0;
	wsReal	R_szDiv	= (wsReal)N_sample/(wsReal)N_thread;
	int		i=0;
	wsReal	j=W0;
	int		N_ret = N_proc != 0 ? 0 : 1;

	switch (N_mode) {
	case DISTFUN_INIT:
		i = 0;
		j = W0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++,j+=R_szDiv) {
			wsUint E = (wsUint)(j+R_szDiv+REAL_CONST(0.5));
			if (E > N_sample)
				E = N_sample;
			Np_idx[i*3] = (wsUint)(j+REAL_CONST(0.5));
			Np_idx[i*3+1] = E;
			Np_idx[i*3+2] = 0;
		}
		N_proc = N_ret = N_thread;
		break;
	case DISTFUN_UNINIT:
		LOG("%d/%d samples processed\n", N_sample, N_sample);
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt_fmt(WISARD_SYST_INVL_DISTFUNC_CMD, N_mode);
//		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}
	return i;
}

int forAllSample_corrsym_equal(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx = (int *)Vp_thrData;
	wsUint	N_sample = ((xCorrThread *)Vp_shareData)->N_sample - 1;
	static wsUint
		N_proc	= 0;
//	wsReal	R_szDiv	= (wsReal)N_sample/(wsReal)N_thread;
	int		i=0;
	int		N_ret = N_proc != 0 ? 0 : N_thread;

	pverbose("mode [%d] thread [%d]\n", N_mode, N_thread);

	switch (N_mode) {
	case DISTFUN_INIT:
		i = 0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++) {
			Np_idx[i*3] = i;
			Np_idx[i*3+1] = i + (wsUint)((wsReal)(N_sample - i) / (wsReal)N_thread) * N_thread;
			Np_idx[i*3+2] = 0;
		}
		N_proc = N_thread;
		break;
	case DISTFUN_UNINIT:
//		LOG("[%d/%d] samples processed\n", N_sample, N_sample);
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt_fmt(WISARD_SYST_INVL_DISTFUNC_CMD, N_mode);
		//		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}
	return N_ret;
}

#ifdef USE_CUDA
int forAllSample_corrGpu(int N_mode, int N_thread, void *Vp_data,
					  void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx		= (int *)Vp_data;
	static wsUint
			N_idxSample;
	wsUint	N_sample	= ((xCorrThread *)Vp_shareData)->N_sample;
	wsUint	N_SNP		= ((xCorrThread *)Vp_shareData)->N_variant;
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		N_idxSample = 0;
		break;
	case DISTFUN_DISTRIBUTE: {
			wsUint N_completedSample = N_idxSample;
			wsUint nK;

			for ( ; i<N_thread ; i++) {
				if ((N_idxSample+i) >= N_sample)
					break;
				/* Calculate K
				 * Set K to use up to 90% of GPU memory */
				wsReal K = (wsReal)(WORKER().getThrData(i)->X_propGPU.totalGlobalMem) /
					(W2*N_SNP*sizeof(wsReal)) * 0.9;

				/* If K is greater than the expected number for each thread,
				 * adjust K */
				if (K > ((N_sample-N_completedSample)/(wsUint)N_thread))
					nK = (int)((N_sample-N_completedSample)/(wsUint)N_thread)+1;
				else
					nK = (int)K;

				/* if K is greater than remained sample,
				 * adjust K */
				if (nK > (N_sample-N_idxSample))
					nK = N_sample-N_idxSample;

				Np_idx[i*2] = N_idxSample;
				Np_idx[i*2+1] = nK;

				N_idxSample += nK;
			}
		} break;
	case DISTFUN_UNINIT:
		LOG("%d/%d samples processed\n", N_sample, N_sample);
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt_fmt(WISARD_SYST_INVL_DISTFUNC_CMD, N_mode);
//		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}
	return i;
}
#endif

cCorrAnalysis::cCorrAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP) : cAnalysis(Cp_inpIO)
{
	N_corSamp	= OPT_ENABLED(founderonly) ? Cp_IO->sizeFounder() :
		Cp_IO->sizeSample();
	wsUint		N_SNP	= Cp_IO->sizeVariant();
	wsFloat**	Ra_data	= Cp_IO->getData();
	Cp_anaPPP = Cp_inpAnaPPP;

	/* If the data is dosage data, force --corpearson */
	if (IS_ASSIGNED(dosage)) {
		if (OPT_ENABLED(corpearson))
			LOG("Sample relatedness matrix with dosage will be computed "
				"using Pearson's correlation coefficient\n");
		else
			LOG("Empirical sample relatedness matrix with dosage will be computed\n");
//		OPTION().assign("corpearson", "1", 1);
//		OPTION().FORCE_OPT_NUMBER(corpearson);
	}

	if (Cp_IO->sizeVariant() == 0 || OPT_ENABLED(indep))
		Cp_corr = new cIdtMatrix(N_corSamp);
	else
		Cp_corr = new cStdMatrix(N_corSamp, N_corSamp);

	wsAlloc(Na_corrN, int*, N_corSamp);
	/* Create corr matrix */
	for (wsUint i=0 ; i<N_corSamp ; i++)
		wsAlloc(Na_corrN[i], int, N_corSamp);
	/* Allocate extra memory if data are complete */
	if (Cp_IO->isDataComplete() == 1) {
		wsAlloc(Ra_ysum , wsReal, N_corSamp);
		wsAlloc(Ra_y2sum, wsReal, N_corSamp);

		/* Calculate ysum and y2sum if data are complete */
		for (wsUint i=0 ; i<N_corSamp ; i++) {
			wsFloat	*Rp_row = Ra_data[i];

			Ra_ysum[i]	= W0;
			Ra_y2sum[i]	= W0;

			for (wsUint j=0 ; j<N_SNP ; j++) {
				Ra_ysum[i]	+= (wsReal)Rp_row[j];
				Ra_y2sum[i]	+= (wsReal)(Rp_row[j]*Rp_row[j]);
			}
		}
	} else {
		Ra_ysum		= NULL;
		Ra_y2sum	= NULL;
	}
}

void cCorrAnalysis::_pairedAltCorLoad(cStream &C_altcor, char *Sp_buf,
	wsMat Ra_corr, wsStrCst Sp_corPath)
{
	wsUint	N_sample	= Cp_IO->sizeSample();
	wsUint i;
	LOG("AltCor format : String header, paired format\n\n");
	//B_is3col = 1;

	wsUint N_s1, N_s2;
	char S_s1[64], S_s2[64];
	wsReal R_val;

	sseMinit(Ra_corr, N_sample, N_sample, WISARD_NAN);

	mDataIdx Xm_retrieved;
	mSamp &Xm_samp = Cp_IO->getSampleData();
	for (i=0 ; C_altcor.gets(Sp_buf, 1024*1024) ; i++) {
		wsUint N_get = sscanf(Sp_buf, "%s\t%s\t" REAL_FMT, S_s1,
			S_s2, &R_val);
		if (N_get != 3)
			halt_fmt(WISARD_INVL_CORFORMAT, i+2, N_get, 3);

		string Sx_s1=S_s1, Sx_s2=S_s2;
		mSamp_it X_r1 = Xm_samp.find(Sx_s1);
		mSamp_it X_r2 = Xm_samp.find(Sx_s2);

		/* Cannot found one of sample in the pair from dataset */
		if (X_r1 == Xm_samp.end() || X_r2 == Xm_samp.end())
			halt_fmt(WISARD_NULL_IID_IN_DATASET, S_s1, S_s2, i+2);

		/* Sample found, but it does not have its data */
		if (X_r1->second.B_isMissing) {
			if (mapIsExist(Xm_retrieved, S_s1)) continue;
			Xm_retrieved[S_s1] = 2;
			LOG("IID [%s] have not its data, will be skipped\n", S_s1);
			continue;
		}
		if (X_r2->second.B_isMissing) {
			if (mapIsExist(Xm_retrieved, S_s2)) continue;
			Xm_retrieved[S_s2] = 2;
			LOG("IID [%s] have not its data, will be skipped\n", S_s2);
			continue;
		}
		/* Mark to be retrived successfully */
		Xm_retrieved[S_s1] = 1;
		Xm_retrieved[S_s2] = 1;

		N_s1 = X_r1->second.N_idx;
		N_s2 = X_r2->second.N_idx;
		Ra_corr[N_s1][N_s2] = Ra_corr[N_s2][N_s1] = R_val;
	}
	vSampPtr& Xv_samp = getIO()->getSample();
	LOOP (i, N_sample) LOOP(j, N_sample)
		if (NA(Ra_corr[i][j]))
			halt("Correlation for sampie pair [%s::%s] and [%s::%s] is lacked!",
				Xv_samp[i]->S_FID.c_str(), Xv_samp[i]->S_IID.c_str(),
				Xv_samp[j]->S_FID.c_str(), Xv_samp[j]->S_IID.c_str());
// 	if (i != (N_sample*(N_sample+1))>>1)
// 		halt_fmt(WISARD_INVL_ALTCORSIZE, Sp_corPath, i,
// 			(N_sample*(N_sample+1))>>1);

	i = 0;
	FOREACH (mDataIdx_it, Xm_retrieved, it)
		if (it->second == 1)
			i++;
}

void cCorrAnalysis::_emmaxAltCorLoad(cStream &C_altcor, char *Sp_buf,
	wsMat Ra_corr, wsStrCst Sp_corPath)
{
	mSamp&	Xm_iid2dt	= Cp_IO->getSampleData();
	LOG("AltCor format : String header, EMMAX format\n\n");

	/* At this point, 8 bytes were already read */
	wsUint N_altSamp[2];
	C_altcor.read(N_altSamp, 2);
	/* Position (N_altSamp * (N_altSamp+1)) << 2 + 16 + 8 bytes */
	size_t N_pos = N_altSamp[0];
	N_pos = ((N_pos * (N_pos+1)) << 2) + 16 + 8;
	C_altcor.seek(N_pos, SEEK_SET);

	/* Read sample headers */
	int* Na_mapIdx = NULL;
	wsAlloc(Na_mapIdx, int, N_altSamp[0]);
	LOOP (i, N_altSamp[0]) {
		char *Sp_label = NULL;

		/* Read size */
		wsUint L_label;
		C_altcor.read(&L_label);
		wsAlloc(Sp_label, char, L_label);
		C_altcor.read(Sp_label, L_label);

		/* Check whether this sample is existing or not */
		mSamp_it X_find = Xm_iid2dt.find(Sp_label);	
		/* Mark skip if this sample is not available */
		if (X_find == Xm_iid2dt.end() || X_find->second.N_idx == SAMP_NODATA)
			Na_mapIdx[i] = SAMP_NODATA;
		else
			/* Have data */
			Na_mapIdx[i] = X_find->second.N_idx;
	}
	/* Move to pos 16 - data */
	C_altcor.seek(16, SEEK_SET);
	LOOP (i, N_altSamp[0]) {
		int I = Na_mapIdx[i];
		/* Skip this row */
		if (I == -1)
			C_altcor.seek((i+1)*sizeof(double), SEEK_CUR);

		for (int j=(int)i ; j>=0 ; j--) {
			double v;
			C_altcor.read(&v);
			int J = Na_mapIdx[j];
			if (J == -1) continue;
			Ra_corr[I][J] = Ra_corr[J][I] = (wsReal)v;
		}
	}
	DEALLOC(Na_mapIdx);
}

cCorrAnalysis::cCorrAnalysis(cIO *Cp_inpIO, const char *Sp_corPath)
	: cAnalysis(Cp_inpIO)
{
	wsUint	N_sample	= Cp_IO->sizeSample();
	wsUint i, j;
	Cp_anaPPP = NULL;

	/* Option checking */
	if (OPT_ENABLED(corpearson) || OPT_ENABLED(ibs) || OPT_ENABLED(hybrid) ||
		OPT_ENABLED(medcor) || OPT_ENABLED(ktau) || OPT_ENABLED(empktau))
		halt_fmt(WISARD_CANT_DO_W_SOMEOPT, "Computation of sample relatedness",
			"--cor");
		//halt_fmt(WISARD_CANT_CALCORR_W_ALTCOR);

	N_corSamp = OPT_ENABLED(founderonly) ? Cp_IO->sizeFounder() :
		Cp_IO->sizeSample();

	Cp_corr = new cStdMatrix(N_sample, N_sample);
	wsReal **Ra_corr = Cp_corr->get();
//	MULTI_MALLOC(Ra_corr , wsReal*, N_sample);
	Na_corrN = NULL;
	/* Create corr matrix */
// 	for (i=0 ; i<N_sample ; i++) {
// 		sseMalloc(Ra_corr[i] , wsReal, N_sample);
// 	}
	Ra_ysum		= NULL;
	Ra_y2sum	= NULL;

	cStrFile	C_altcor(Sp_corPath, "Alternative correlation file", 1);
	char		*Sp_buf	= NULL;
	char		*a, *b;
	vStr		Xv_header;
	mDataIdx	Xm_header;
	char		B_haveHeader = 0;
	vector<vector<wsReal> >
				M_inCorr;

	// GRM?
	wsAlloc(Sp_buf, char, 1024*1024);
	if (!stricmp(Sp_corPath + strlen(Sp_corPath) - 4, ".grm")) {
		vInt		Xv_found(Cp_IO->sizeSample());
		vInt		Xv_map;
		vSampPtr&	Xv_samp	= Cp_IO->getSample();
		mDataIdx	Xm_sampkey2idx;
		LOGnote("Given sample relationship matrix is GRM format\n");

		i = 0;
		FOREACHDO (vSampPtr_it, Xv_samp, I, i++)
			Xm_sampkey2idx.insert(make_pair((*I)->S_FID + "|" + (*I)->S_IID, i));

		/* GRM always have header */
		B_haveHeader = 1;
		FOREACH (vInt_it, Xv_found, I) *I = -1;

		/* Read ids first */
		char S_newPath[MAX_PATH];
		sprintf(S_newPath, "%s.id", Sp_corPath);
		cStrFile C_altid(S_newPath, "GRM id file");
		for (i=0 ; C_altid.gets(Sp_buf, 1024*1024) ; i++) {
			char *a = NULL, *b = NULL;
			/* Get FID & IID */
			getString(&Sp_buf, &a); //Sp_buf = FID
			if (!a) halt("Line [%d] in sample relatedness matrix file [%s] "
				"is invalid", i+1, S_newPath);
			getString(&a, &b); // a = IID


			/* Exists in the dataset? */
			string S_sampkey = Sp_buf;
			S_sampkey += "|";
			S_sampkey += a;
			mDataIdx_it X_r1 = Xm_sampkey2idx.find(S_sampkey);
			if (X_r1 != Xm_sampkey2idx.end()) {
				Xv_found[X_r1->second] = i;
				Xv_map.push_back(X_r1->second);
			} else
				Xv_map.push_back(-1);

			Xv_header.push_back(Sp_buf);
			Xm_header.insert(make_pair(Sp_buf, i));
		}
		C_altid.close();

		/* Check all samples have their idx */
		i = 0;
		FOREACHDO (vInt_it, Xv_found, I, i++)
			if (*I == -1) {
				xSample *X = Cp_IO->getSample()[i];
				halt("Sample [%s::%s] not exists in the altcor file",
					X->S_FID.c_str(), X->S_IID.c_str());
			}

		wsUint N_sampGRM = (wsUint)Xv_header.size();
		LOG("[%d] samples are available in GRM file\n", N_sampGRM);
		/* Read cells */
		sprintf(S_newPath, "%s.bin", Sp_corPath);
		C_altcor.open(S_newPath, "GRM file", 0, 1);
		wsFloat *Ra_vec = NULL;
		wsAlloc(Ra_vec, wsFloat, N_sampGRM);
		for (i=0 ; i<N_sampGRM ; i++) {
			/* In default, 4 */
			if (C_altcor.read(Ra_vec, i+1) != (i+1))
				halt("Failed to read GRM file");
			wsUint N_i1 = Xv_map[i];
			wsReal *Ra_tgt = Ra_corr[N_i1];
			for (wsUint j=0 ; j<=i ; j++) {
				wsUint N_i2 = Xv_map[j];
				Ra_corr[N_i2][N_i1] = Ra_tgt[N_i2] = Ra_vec[j];
			}
			/* Make sym to tt */
			for (wsUint j=0 ; j<=i ; j++)
				Ra_corr[j][i] = Ra_corr[i][j];
		}
		C_altcor.close();
	} else {
		/* Check the handle is sane */
		if (C_altcor.isFailed())
			halt_fmt(WISARD_FAIL_OPEN_FILEWITHDESC, "Sample relatedness matrix file",
				Sp_corPath);
		char B_loaded = 0;

		/* Read the first 8 byte to check is it EMMA_KIN */
		char S_checkEmmaxKin[9] = { 0, };
		size_t L_read = C_altcor.read(S_checkEmmaxKin, 8);
		if (L_read != 8 || strcmp("EMMA_KIN", S_checkEmmaxKin)) {
			C_altcor.rewind();
			/* Read data and process */
			for (i=0 ; C_altcor.gets(Sp_buf, 1024*1024) ; i++) {
				vector<wsReal> V_curLine;

				/* Non-matrix style loading */
				if (i == 0 && !memcmp("SAMP1", Sp_buf, 4) &&
					(Sp_buf[5] == ' ' || Sp_buf[5] == '\t')) {
					_pairedAltCorLoad(C_altcor, Sp_buf, Ra_corr, Sp_corPath);
					B_loaded = 1;
					break;
				}		
			
				/* Check the number of lines and matching the expected size of samples */
		// 		if (i > N_sample)
		// 			halt_fmt(WISARD_SYST_INVL_ALTCORSIZE, i, Sp_corPath, N_sample);

				if (i == 0) for (j=0,a=Sp_buf ; a ; j++,a=b) {
					getString(&a, &b);
		// 			if (j >= N_sample) {
		// 				if (!a[0])
		// 					break;
		// 				halt_fmt(WISARD_INVL_CORFILE_MORESZSAMP, i+1, /*Sp_corPath, */N_sample);
		// 				// 				halt("Correlation alternation file [%s] have %d th line "
		// 				//					"have more than %d samples", Sp_corPath, i+1, N_sample);
		// 
		// 			}
					/* For the first line, check is it header or not */
					char	*E		= NULL;
					wsReal	R_vals	= (wsReal)str2dbl(a, &E);
					if (E[0] && !B_haveHeader) {
						LOG("Non-numerical element[%s->%g] found in 1st line of alternative "
							"correlation file, it will be regarded as the header\n",
							a, R_vals);
						B_haveHeader = 1;
					}
					Xv_header.push_back((string)a);
					Xm_header.insert(make_pair(a, j));
					//Ra_corr[i][j] = (wsReal)atof(a);
				} else for (j=0,a=Sp_buf ; a ; j++,a=b) {
					getString(&a, &b);
					if (j >= Xv_header.size()) {
						if (!a[0])
							break;
						halt_fmt(WISARD_INVL_CORFILE_MORESZSAMP, i+1, /*Sp_corPath, */Xv_header.size());
		// 				halt("Correlation alternation file [%s] have %d th line "
		//					"have more than %d samples", Sp_corPath, i+1, N_sample);
					}
					V_curLine.push_back((wsReal)atof(a));
					//Ra_corr[i-1][j] = (wsReal)atof(a);
				}

 				if (j != Xv_header.size())
 					halt_fmt(WISARD_INVL_CORFILE_LESSSZSAMP, i+1,/* Sp_corPath,*/ Xv_header.size());
				if (i != 0)
					M_inCorr.push_back(V_curLine);
		//			halt("Correlation alternation file [%s] have %d th line have "
		//				"less than %d samples (%d retrieved)", Sp_corPath, i+1,
		//				N_sample, j);
			}
			/* Possibility
			 * 1: hh=0 and (i+1) == N_sample : (MAYBE) HAVE HEADER
			 * 2: hh=0 and i == N_sample     : NO HEADER
			 * 3: hh=1 and (i+1) == N_sample : HAVE HEADER
			 * 4: hh=1 and i == N_sample     : ERROR */
			if (i == N_sample) {
				if (B_haveHeader)
					halt("This alternative correlation file have header, but # of samples are less");
				/* Shift all of them +1 */
				for (int X=N_sample-2 ; X>=0 ; X--)
					memcpy(Ra_corr[X+1], Ra_corr[X], sizeof(wsReal)*N_sample);
				/* Nullify header, it is ACTUALLY not header */
				wsUint x=0, y;
				FOREACHDO (vStr_it, Xv_header, X, x++)
					Ra_corr[0][x] = (wsReal)atof(X->c_str());
				x = 1;
				FOREACHDO (vector<vector<wsReal> >::iterator, M_inCorr, X, x++) {
					y = 0;
					FOREACHDO (vector<wsReal>::iterator, (*X), XX, y++)
						Ra_corr[x][y] = *XX;
					X->clear();
				}
				M_inCorr.clear();

				Xv_header.clear();

				LOG("AltCor format : No header, matrix format\n");
			} else if (!B_loaded) { /* Different line? */
				if (B_haveHeader)
					LOG("AltCor format : String header, matrix format\n");
				else
					LOG("AltCor format : Non-string header, matrix format\n");

				/* Check map and vector have same length */
				if (Xv_header.size() != Xm_header.size())
					halt("Duplicated key found in header of alternative correlation matrix");

				/* Check correctness of header */
				vSampPtr&	Xa_sample	= getIO()->getSample();
				wsUint		N_exist		= 0;
				FOREACH (vSampPtr_it, Xa_sample, it) {
					mDataIdx_it X_find = Xm_header.find((*it)->S_IID);
					if (X_find != Xm_header.end())
						N_exist++;
				}
				LOG("[%d/%d] samples matched from alternative correlation matrix\n",
					N_exist, Xm_header.size());
				if (N_exist != N_sample)
					halt("Insufficient number of samples found[%d] in alternative "
						"correlation matrix", N_exist);
				wsUint		*Na_idx		= NULL;
				sseMalloc(Na_idx, wsUint, N_sample);
				wsUint x = 0;
				FOREACH (vSampPtr_it, Xa_sample, it) {
					mDataIdx_it X_find = Xm_header.find((*it)->S_IID);
					if (X_find != Xm_header.end())
						Na_idx[x++] = X_find->second;
					X_find->second = -1;
				}
				/* Check which sample were not included */
				cExporter	*Cp_E = NULL;
				wsUint		N_missAltCor = 0;
				FOREACH (mDataIdx_it, Xm_header, it)
					if (it->second != -1) {
						if (!Cp_E) {
							Cp_E = cExporter::summon("altcor.miss.lst");
							Cp_E->put("IID\n");
							LOGoutput("Samples in alternative correlation but not "
								"in dataset are exported to [%s.altcor.miss.lst]\n",
								OPT_STRING(out));
						}
						if (N_missAltCor < 5)
							LOG("Sample [%s] exists in alternative correlation file but "
								"not exists in dataset\n", it->first.c_str());
						Cp_E->fmt("%s\n", it->first.c_str());
						N_missAltCor++;
					}
				if (Cp_E) {
					delete Cp_E;
					LOG("[%d] samples in alternative correlation are missing in dataset\n", N_missAltCor);
				}

				/* Rebuild */
				for (wsUint I=0 ; I<N_sample ; I++) {
					vector<wsReal> &X = M_inCorr[Na_idx[I]];
					for (wsUint J=0 ; J<N_sample ; J++)
						Ra_corr[I][J] = X[Na_idx[J]];
				}
				sseFree(Na_idx);
			}
		} else
			_emmaxAltCorLoad(C_altcor, Sp_buf, Ra_corr, Sp_corPath);
	}

	/* Process --cordiag1 */
	if (OPT_ENABLED(cordiag1))
		for (i=0 ; i<N_sample ; i++)
			if (Ra_corr[i][i] < W1)
				Ra_corr[i][i] = W1;

	/* Adjusting if --x2 */
	/* --x2 */
	wsUint *Na_corIdx = NULL;
	if (OPT_ENABLED(x2)) {
		mFam&		Xm_fam	= getIO()->getFamilyData();
		vSampPtr&	Xv_samp	= getIO()->getSample();

		/* For all sample, do compensation */
		FOREACH (mFam_it, Xm_fam, i) {
			vInt&	Xv_mem	= i->second.Xv_members;

			for (wsUint j=0 ; j<(wsUint)Xv_mem.size() ; j++) {
				char*	Ba_fill	= NULL;
				wsCalloc(Ba_fill, char, N_corSamp);

				wsUint J = Xv_mem[j];
				wsUint _J = Na_corIdx ? Na_corIdx[J] : J+1;
				if (!Na_corIdx || (Na_corIdx && _J)) {
					for (wsUint k=0 ; k<=j ; k++) {
						wsUint K = Xv_mem[k];
						wsUint _K = Na_corIdx ? Na_corIdx[K] : K+1;

						/* Do not fill if not target */
						if (Na_corIdx && _K == 0) continue;

						map<string,xPDDTnode> X_map;
						buildGraph(X_map, Xv_samp[K], 0, buildGraph(X_map, Xv_samp[J], 1));
						/* Xp_s1���� Xp_s2�� ã�´� */
						vector<wsReal> Xa_scores;
						wsReal R_mul = W0;
						getX2coeff(X_map, Xa_scores, Xv_samp[J], Xv_samp[J], Xv_samp[K], 0, -1, W1, 1);

						if (Xa_scores.size() == 1)
							R_mul = Xa_scores[0];
						else if (Xa_scores.size() == 2)
							R_mul = (Xa_scores[0] + Xa_scores[1]) / W2;

						pverbose("[%s] <-> [%s] : %g\n", Xv_samp[J]->S_IID.c_str(), Xv_samp[K]->S_IID.c_str(), R_mul);

						/* Mul */
						Ra_corr[_J-1][_K-1] *= R_mul;
						Ra_corr[_K-1][_J-1] = Ra_corr[_J-1][_K-1];
						Ba_fill[_K-1] = 1;
					}
				}
				/* Set 0 */
				for (wsUint k=0 ; k<N_corSamp ; k++)
					if (!Ba_fill[k]) Ra_corr[_J-1][k] = W0;

				DEALLOC(Ba_fill);
			}
		}
	}

	/* Export if need */
	/* Name as vStr */
	vStr		Sa_names;
	vSampPtr	&Xa_samples = getIO()->getSample();
	FOREACH (vSampPtr_it, Xa_samples, i)
		Sa_names.push_back((*i)->S_IID);
	/* Export as matrix form */

	if (IS_ASSIGNED(makecor))
		exportMatrix("alt.cor", Ra_corr, N_sample, N_sample, &Sa_names);

	DEALLOC(Sp_buf);
}

void getX2coeff(map<string,xPDDTnode> &X_map, vector<wsReal> &Xa_scores,
	xSample *Xp_src, xSample *Xp_cur, xSample *Xp_tgt, int B_visitingTwin,
	int B_visitingFemTopNode, wsReal R_score,
	char B_isX)
{
	xPDDTnode &X_currNode = X_map[Xp_cur->S_IID];

	/* Target found */
	if (Xp_cur == Xp_tgt) {
		/* twin�� ����Ǿ����ٸ�?*2�Ǿ��?�Ѵ� */
		Xa_scores.push_back(R_score * (1+B_visitingTwin) * (1+(Xp_tgt->N_sex==1)));
		return;
	}

	/* For all connected node, search and expand */
	FOREACH (vSampPtr_it, X_currNode.Xa_next, it) {
		/* Twin�� �ְ�, �׷����� �����Ѵٸ� ����ؾ���?twin�̹Ƿ�
			* �÷��׸� �����ش� */
		if (B_visitingTwin == 0) {
			FOREACH (vSampPtr_it, (*it)->Xp_twins, iit)
				if (mapIsExist(X_map, (*iit)->S_IID)) {
					B_visitingTwin = 1;
					break;
				}
		}

		/* If this is top level parent and FEMALE */
		if ((X_currNode.N_depth < X_map[(*it)->S_IID].N_depth)) {
			if (Xp_cur->N_sex == 2)
				B_visitingFemTopNode = 1;
			else
				B_visitingFemTopNode = 0;
		}

		wsReal R_mult = W1;

		/* Special for X chromosome */
		if (Xp_cur->N_sex != 1 && Xp_cur->N_sex != 2)
			/* Current have no sex then set to 0 */
			R_mult = W0;
		else {
			/* Target is mother of current or
				* Current is mother of target */
// 			if (*it == Xp_cur->Xp_mat || (*it)->Xp_mat == Xp_cur)
// 				R_mult = W2;
// 			/* Target is father of current or
// 				* Current is father of target */
// 			else if (*it == Xp_cur->Xp_pat ||
// 				(*it)->Xp_pat == Xp_cur) {
				/* Both male */
				if ((Xp_cur->N_sex == 1) && ((*it)->N_sex == 1))
					R_mult = W0;
				/* Me and target have diff. sex (NOT MISSING) */
				else if (Xp_cur->N_sex==2 && (*it)->N_sex==1) 
					R_mult = W1;
				else if (Xp_cur->N_sex==1 && (*it)->N_sex==2)
					R_mult = W2;
				else
					R_mult = W1;
//			}
		}
// 		pverbose("	[%s]->[%s], %g*%g=%g, twin %d, fem %d\n", Xp_cur->S_IID.c_str(),
// 			(*it)->S_IID.c_str(), R_score, R_mult, R_score*R_mult,
// 			B_visitingTwin, B_visitingFemTopNode);
		getX2coeff(X_map, Xa_scores, Xp_src, *it, Xp_tgt, B_visitingTwin,
			B_visitingFemTopNode, R_score*R_mult, B_isX);
	}
}

cCorrAnalysis::~cCorrAnalysis()
{
//	LOG("~cCorrAnalysis() called\n");

	DEALLOC(Ra_ysum);
	DEALLOC(Ra_y2sum);
	if (Na_corrN) for (wsUint i=0 ; i<N_corSamp ; i++) {
//		sseFree(Ra_corr[i]);
		DEALLOC(Na_corrN[i]);
	}
	delete Cp_corr;
//	DEALLOC(Ra_corr);
	DEALLOC(Na_corrN);
}

template <typename T>
void cCorrAnalysis::_getGRM(wsUint N_corSamp, wsUint N_vrt, wsVecCst Ra_ppp, T** Na_data, wsMat Ra_corr)
{
	cIO*		Cp_IO			= getIO();
	xMaf*		Xp_maf			= Cp_IO->getMAF();
	wsFloat**	Rp_data			= NULL;
	T**			Np_data			= NULL;
	wsUint*		Na_idx			= NULL;
	wsUint*		Na_corIdx		= NULL; /* Index in cor matrix, starting from 1 */
	/* Draw mask */
	wsFloat*		Ra_corrMask	= NULL;
	char*		Na_corrMask	= NULL;
	wsUint		N_sexMarker	= 0;
	wsUint		N_noIncCor		= Cp_anaPPP->getCorMask(&Ra_corrMask,
		&Na_corrMask, &N_sexMarker);

	/* Run calculation */
	if (OPT_ENABLED(founderonly)) {
		wsCalloc(Na_corIdx, wsUint, Cp_IO->sizeSample());

		const char	*Bp_isFounder	= Cp_IO->getIsFounder();
		wsFloat**	Ra_oriData		= Cp_IO->getData();
		T**			Na_oriData		= Na_data;

		/* In this case, reconstruct Rp_data to having founders only */
		wsAlloc(Rp_data, wsFloat*, N_corSamp);
		wsAlloc(Np_data, T*, N_corSamp);
		wsAlloc(Na_idx, wsUint, N_corSamp);
		for (wsUint i=0,_i=0 ; i<N_corSamp ; i++) {
			if (Bp_isFounder[i] == 0)
				continue;

			Rp_data[_i]	= Ra_oriData[i];
			Np_data[_i]	= Na_oriData[i];
			Na_idx[_i]	= i;
			Na_corIdx[i]	= _i+1;
			_i++;
		}
	} else {
		if (!OPT_ENABLED(ktau))
			Rp_data = Cp_IO->getData();
		Np_data = Na_data;
	}

#ifdef USE_CORR00
	/* Generate bitwise missing map */
	WISARD_pc**	Na_mMissing = NULL;
	wsCalloc(Na_mMissing, WISARD_pc*, N_corSamp);
	wsUint		N_mBlk		= (N_vrt+(WISARD_szpc-1))>>WISARD_sftpc;
	vSampPtr&	Xa_samp		= Cp_IO->getSample();

	if (!Np_data) {
		if (!IS_ASSIGNED(dosage)) halt("SYSERR: Non-dosage data have NULL genotype!");
	} else for (wsUint i=0 ; i<N_corSamp ; i++) {
		/* Do not make mask if the sample is complete */
		if ((Na_idx && Xa_samp[Na_idx[i]]->B_isComplete) ||
			Xa_samp[i]->B_isComplete) continue;

		sseCalloc(Na_mMissing[i], WISARD_pc, N_mBlk);
		for (wsUint j=0 ; j<N_vrt ; j++)
			if (!isMissing(Np_data[i][j]) && Na_corrMask[j])
				Na_mMissing[i][j>>WISARD_sftpc] |= (WISARD_pc)1 << (j&(WISARD_szpc-1));
	}
#endif

// 		if (OPT_ENABLED(indep)) {
// 			for (wsUint i=0 ; i<N_corSamp ; i++) {
// 				memset(Ra_corr[i], 0x00, sizeof(wsReal)*N_corSamp);
// 				Ra_corr[i][i] = W1;
// 				for (wsUint j=0 ; j<N_corSamp ; j++)
// 					Na_corrN[i][j] = N_corSamp;
// 			}
// 		} else {
		xCorrThread	X_ct;
		X_ct.Na_gt			= NULL;
		X_ct.Xp_SNP			= &Cp_IO->getVariant();
		X_ct.N_sample		= N_corSamp;
		X_ct.N_variant		= Cp_IO->sizeVariant();
		X_ct.Na_corrN		= Na_corrN;
		X_ct.Ra_corr		= getCorr();
		X_ct.B_isComplete	= Cp_IO->isDataComplete();
		X_ct.Ra_corrMask	= Ra_corrMask;//Cp_IO->getCorrMask();
		X_ct.Na_corrMask	= Na_corrMask;
		vVariant&	Xv_vrt		= Cp_IO->getVariant();
		if (OPT_ENABLED(verbose)) {
			cExporter*	Cp_incStat	= cExporter::summon("cor.inclusion.variant");
			for (wsUint i=0 ; i<X_ct.N_variant ; i++)
				Cp_incStat->fmt("%s	%g	%g	%d\n", Xv_vrt[i].name,
					Xp_maf[i].R_maf, Ra_ppp[i], X_ct.Ra_corrMask[i]==W0?0:1);
			delete Cp_incStat;
		}
		if (OPT_ENABLED(ibs)) {
			LOG("IBS will be calculated instead of correlation\n");
			X_ct.Ra_data	= NULL;
			X_ct.Na_data	= (char **)Np_data;
			WORKER().run(cCorrAnalysis::calcIBS, forAllSample_corr, &X_ct, NULL);
		} else if (OPT_ENABLED(hamming)) {
			LOG("Hamming distance will be calculated instead of correlation\n");
			X_ct.Ra_data	= NULL;
			X_ct.Na_data	= (char **)Np_data;
			WORKER().run(cCorrAnalysis::calcHamming, forAllSample_corr, &X_ct, NULL);
		} else if (OPT_ENABLED(bn)) {
			LOG("Balding-Nichols model will be used\n");
			X_ct.Ra_data	= NULL;
			X_ct.Na_data	= (char **)Np_data;

			// For each SNP
			wsVec X = sseVector(X_ct.N_sample);
			xMaf* Xp_maf = Cp_IO->getMAF();
			for (wsUint i=0 ; i<X_ct.N_variant ; i++) {
				if (X_ct.Ra_corrMask[i] == W0) continue;
				wsReal mean = Xp_maf[i].R_allMaf * W2;
				wsReal scale = sqrt(W1 / (W1 - mean / W2) / mean);

				/* Set genotype */
				wsUint _j=0;
				for (wsUint j=0 ; j<X_ct.N_sample ; j++) {
					char G = (char)Np_data[j][i];
					if (G != WISARD_NA) X[_j++] = (wsReal)G;
				}
				sseVsC(X, mean, X, _j, scale);
				LOOP(j, X_ct.N_sample)
					sseVpCadd(X, X[j], Ra_corr[j], j+1);
			}
			//WORKER().run(cCorrAnalysis::calcBN1, forAllSample_corr, &X_ct, NULL);

		} else if (OPT_ENABLED(corpearson)) {
			LOG("Pearson correlation will be calculated\n");
			X_ct.Ra_data	= Rp_data;
			X_ct.Na_data	= (char **)Np_data;
			WORKER().run(cCorrAnalysis::calcCorrV2pearson, forAllSample_corr, &X_ct, NULL);
		} else if (OPT_ENABLED(ktau) || OPT_ENABLED(empktau)) {
			LOG("Kendall's tau rank correlation will be calculated\n");
			X_ct.Ra_data	= Rp_data;
			X_ct.Na_data	= (char **)Np_data;

			if (Cp_IO->isDataComplete()) {
				LOG("Confirmed that the data have complete genotype, use faster manner\n");

				/* For all sample, calculate combination comparison result */
				wsAlloc(X_ct.Na_gt, WISARD_pc *, N_corSamp);
				wsAlloc(X_ct.Na_neq, WISARD_pc *, N_corSamp);
				wsAlloc(X_ct.Na_same, wsUint *, N_corSamp);

				/* Count how much of SNPs will be included */
				__int64 N_sz = 0;
				if (X_ct.Ra_corrMask) {
					for (wsUint k=0 ; k<X_ct.N_variant ; k++)
						if (X_ct.Ra_corrMask[k] != W0) N_sz++;
				} else
					N_sz = X_ct.N_variant;
				X_ct.N_szComb = (N_sz*(N_sz-1))>>1;
				LOG("Use [%d] SNPs to calculate Kendall's tau\n", N_sz);
				/* Calculate the number of combinations */
				__int64 N_comb		= X_ct.N_szComb;
				__int64 N_szByte64	= (N_comb+(WISARD_szpc-1))>>WISARD_sftpc;
				X_ct.N_last		= (wsUint)(N_szByte64&(WISARD_szpc-1));
				if (N_szByte64 >= 0xffffffff) halt("Too large number of combinations are required, cannot proceed!");
				wsUint	N_szBlk	= (wsUint)N_szByte64;
				X_ct.N_szBlk	= N_szBlk;
				X_ct.N_szSNP	= (wsUint)N_sz;
				LOG("Use [%d] bytes for each sample\n", N_szBlk);


				/* Do phase 1 */
				WORKER().run(cCorrAnalysis::calcCorrV2tauS1, forAllSample_corr, &X_ct, NULL);

				/* Do phase 2 */
				WORKER().run(cCorrAnalysis::calcCorrV2tauS2, forAllSample_corr, &X_ct, NULL);

				/* Cleanup */
				for (wsUint i=0 ; i<N_corSamp ; i++) {
					DEALLOC(X_ct.Na_gt[i]);
					DEALLOC(X_ct.Na_neq[i]);
					DEALLOC(X_ct.Na_same[i]);
				}
				DEALLOC(X_ct.Na_gt);
				DEALLOC(X_ct.Na_neq);
				DEALLOC(X_ct.Na_same);
			} else {
				/* Otherwise, use slower manner */
				WORKER().run(cCorrAnalysis::calcCorrV2tau, forAllSample_corr, &X_ct, NULL);
			}
		} else {
			X_ct.Ra_data	= Rp_data;
			X_ct.Na_data	= (char **)Np_data;
#ifdef USE_CORR00
			if (OPT_ENABLED(ktau) || OPT_ENABLED(medcor)) {
				WORKER().run(cCorrAnalysis::calcCorrV2, forAllSample_corr, &X_ct, NULL,
					sizeof(int));
			} else {
				if (OPT_NUMBER(thread) == 1)
					sym_sseMpMt(Rp_data, N_corSamp, N_vrt, Ra_corr);
				else
					WORKER().run(thr_sym_sseMpMt, forAllSample_corrsym_equal,
					&X_ct, NULL, sizeof(int) * 3);
				wsUint *Na_cachedN = NULL;
				wsCalloc(Na_cachedN, wsUint, N_corSamp);

				if (OPT_NUMBER(thread) == 1) for (wsUint i=0 ; i<N_corSamp ; i++) {
					for (wsUint j=0 ; j<=i ; j++) {
						WISARD_pc N = 0;

						if (Na_mMissing[i]) {
							if (Na_mMissing[j]) for (wsUint k=0 ; k<N_mBlk ; k++)
								/* Noncomplete - Noncomplete case */
								N += WISARD_POPCNT(Na_mMissing[i][k] & Na_mMissing[j][k]);
							else {
								/* Noncomplete - Complete case */
								if (Na_cachedN[i] == 0) for (wsUint k=0 ; k<N_mBlk ; k++)
									Na_cachedN[i] += (wsUint)WISARD_POPCNT(Na_mMissing[i][k]);
								N = Na_cachedN[i];
							}
						} else {
							if (Na_mMissing[j]) {
								/* Complete - Noncomplete case */
								if (Na_cachedN[j] == 0) for (wsUint k=0 ; k<N_mBlk ; k++)
									Na_cachedN[j] += (wsUint)WISARD_POPCNT(Na_mMissing[j][k]);
								N = Na_cachedN[j];
							} else
								/* Complete - Complete case */
								N = N_vrt-N_noIncCor;
						}

						Ra_corr[i][j]	/= (wsReal)N;
						Ra_corr[j][i]	= Ra_corr[i][j];
						Na_corrN[i][j]	= (int)N;
					}
				} else {
					for (wsUint i=0 ; i<N_corSamp ; i++) {
						if (!Na_mMissing[i]) continue;

						/* Noncomplete - Complete case */
						for (wsUint k=0 ; k<N_mBlk ; k++)
							Na_cachedN[i] += (wsUint)WISARD_POPCNT(Na_mMissing[i][k]);
					}

					X_ct.Na_corrN		= Na_corrN;
					X_ct.N_mBlk			= N_mBlk;
					X_ct.N_noIncCor		= N_noIncCor;
					X_ct.Na_cachedN		= Na_cachedN;
					X_ct.Na_mMissing	= Na_mMissing;

					WORKER().run(thr_count, forAllSample_corrsym_equal,
						&X_ct, NULL, sizeof(int) * 3);
				}
				DEALLOC(Na_cachedN);

				if (OPT_ENABLED(verbose)) {
					char S_fn[256];
					sprintf(S_fn, "%s.empi.corN", OPT_STRING(out));
					FILE *fp = fopen(S_fn, "w+");
					for (wsUint i=0 ; i<N_corSamp ; i++) {
						for (wsUint j=0 ; j<=i ; j++)
							fprintf(fp, "%d", Na_corrN[i][j]);
						fprintf(fp, "\n");
					}
					fclose(fp);
				}
			}

			/* Adjust diagonal */
			if (!OPT_ENABLED(ktau)) {
				cTimer q;
				q.start();

				if (OPT_ENABLED(medcor)) {
					wsReal *Ra_vars = sseVector(N_vrt);

					for (wsUint i=0 ; i<N_corSamp ; i++) {
						wsUint _j = 0;
						for (wsUint j=0 ; j<N_vrt ; j++) {
							T N_geno = Na_data[i][j];
							if (isMissing(N_geno)) continue;
							wsReal R_p = Ra_ppp[j];
							if (N_geno) {
								wsReal R_2p = R_p*W2;
								Ra_vars[_j++] = ((wsReal)(N_geno*N_geno) - (W1 + R_2p)*(wsReal)N_geno + SQR(R_p) * W2)
									/ (R_2p * (W1 - R_p));
							} else
								Ra_vars[_j++] = R_p / (W1 - R_p);
						}

						if (_j)
							Ra_corr[i][i] = quickmedian(Ra_vars, _j) + REAL_CONST(1.0);
						else
							Ra_corr[i][i] = REAL_CONST(1.0);
					}
				} else for (wsUint i=0 ; i<N_corSamp ; i++) {
					wsUint _j = 0;
					wsReal R_s = W0;
					for (wsUint j=0 ; j<N_vrt ; j++) {
						T N_geno = Na_data[i][j];
						if (isMissing(N_geno)) continue;
						wsReal R_p = Ra_ppp[j];
						if (N_geno) {
							wsReal R_2p = R_p*W2;
							R_s += ((wsReal)(N_geno*N_geno) - (W1 + R_2p)*(wsReal)N_geno + SQR(R_p) * W2)
								/ (R_2p * (W1 - R_p));
						} else
							R_s += R_p / (W1 - R_p);
						_j++;
					}
					if (_j)
						Ra_corr[i][i] = W1 + R_s / (wsReal)_j;
				}
				LOG("Adjustment taken [%s]\n", q.getReadable());
			}
#else
#ifdef USE_CUDA
			/* Each thread requires two int element:
				*  starting point, and calc length */
			WORKER().run(calcCorrV2_gpu, forAllSample_corrGpu, &X_ct, NULL,
				sizeof(int)*2);
#else
#if 0
			WORKER().run(calcCorrV2_equal, forAllSample_corr_equal, &X_ct, NULL,
				sizeof(int)*3);
#else
			WORKER().run(calcCorrV2, forAllSample_corr, &X_ct, NULL,
				sizeof(int));
#endif
#endif
#endif
		}

		/* Force diagonal terms to 1 if --cordiag1 */
		if (OPT_ENABLED(cordiag1))
			for (wsUint i=0 ; i<N_corSamp ; i++)
				if (Ra_corr[i][i] < W1)
					Ra_corr[i][i] = W1;
//		}

	/* --x2 */
	if (OPT_ENABLED(x2)) {
		mFam&		Xm_fam	= Cp_IO->getFamilyData();
		vSampPtr&	Xv_samp	= Cp_IO->getSample();

		/* For all sample, do compensation */
		FOREACH (mFam_it, Xm_fam, i) {
			vInt&	Xv_mem	= i->second.Xv_members;

			for (wsUint j=0 ; j<(wsUint)Xv_mem.size() ; j++) {
				char*	Ba_fill	= NULL;
				wsCalloc(Ba_fill, char, N_corSamp);

				wsUint J = Xv_mem[j];
				wsUint _J = Na_corIdx ? Na_corIdx[J] : J+1;
				if (!Na_corIdx || (Na_corIdx && _J)) {
					for (wsUint k=0 ; k<=j ; k++) {
						wsUint K = Xv_mem[k];
						wsUint _K = Na_corIdx ? Na_corIdx[K] : K+1;

						/* Do not fill if not target */
						if (Na_corIdx && _K == 0) continue;

						map<string,xPDDTnode> X_map;
						buildGraph(X_map, Xv_samp[K], 0, buildGraph(X_map, Xv_samp[J], 1));
						/* Xp_s1���� Xp_s2�� ã�´� */
						vector<wsReal> Xa_scores;
						wsReal R_mul = W0;
						getX2coeff(X_map, Xa_scores, Xv_samp[J], Xv_samp[J], Xv_samp[K], 0, -1, W1, 1);

						if (Xa_scores.size() == 1)
							R_mul = Xa_scores[0];
						else if (Xa_scores.size() == 2)
							R_mul = (Xa_scores[0] + Xa_scores[1]) / W2;

						pverbose("[%s] <-> [%s] : %g\n", Xv_samp[J]->S_IID.c_str(), Xv_samp[K]->S_IID.c_str(), R_mul);

						/* Mul */
						Ra_corr[_J-1][_K-1] *= R_mul;
						Ra_corr[_K-1][_J-1] = Ra_corr[_J-1][_K-1];
						Ba_fill[_K-1] = 1;
					}
				}
				/* Set 0 */
				for (wsUint k=0 ; k<N_corSamp ; k++)
					if (!Ba_fill[k]) Ra_corr[_J-1][k] = W0;

				DEALLOC(Ba_fill);
			}
		}
	}
		
#ifdef USE_CORR00
	for (wsUint i=0 ; i<N_corSamp ; i++)
		sseFree(Na_mMissing[i]);
	DEALLOC(Na_mMissing);
#endif	
	/* Deallocate pointer array */
	if (OPT_ENABLED(founderonly)) {
		DEALLOC(Rp_data);
		DEALLOC(Np_data);
	}
// 	} else {
// 		for (wsUint i=0  ; i<N_corSamp ; i++) {
// 			memset(Ra_corr[i], 0x00, sizeof(wsReal)*N_corSamp);
// 			Ra_corr[i][i] = W1;
// 			memset(Na_corrN[i], 0x00, sizeof(int)*N_corSamp);
// 		}
}

void cCorrAnalysis::run()
{
	// 150520 REMOVED
//	if (IS_ASSIGNED(dosage)) return;
	wsUint		N_vrt		= Cp_IO->sizeVariant();
//	xMaf*		Xp_maf		= Cp_IO->getMAF();
	/* Draw mask */
	wsFloat*	Ra_corrMask	= NULL;
	char*	Na_corrMask	= NULL;
	wsUint	N_sexMarker	= 0;
	wsUint	N_noIncCor	= Cp_anaPPP->getCorMask(&Ra_corrMask, &Na_corrMask,
		&N_sexMarker);

	//vMarker	&Xv_mk		= Cp_IO->getMarkers();
	wsRealCst*	Ra_ppp		= Cp_anaPPP->getPPP();

	if (N_sexMarker != 0)
		LOGnote("[%d] sex-chromosome variants are excluded from correlation "
			"analysis\n", N_sexMarker);
			

	/* If no SNPs are incorporated */
	if (N_noIncCor == N_vrt) {
		/* Halt if it is impossible to calculate empirical correlation
		 * or IBS matrix */
		if (!OPT_ENABLED(kinship) && !IS_ASSIGNED(cor)) {
 			LOGwarn("Warning : There is no variant to available to calculate empirical "
 				"correlation, the range of --cormaf should be extended!\n");
			/* Force --indep */
			OPTION().assign("indep", "1", 1);
			OPTION().FORCE_OPT_NUMBER(indep);
			delete Cp_corr;
			Cp_corr = new cIdtMatrix(N_corSamp);
		}
		LOGwarn("No variants available for correlation computation\n");
	} else if ((N_vrt-N_noIncCor) < WINSARD_WARNSIZE_EMPCORCALC) {
		if (OPT_ENABLED(ibs))
			LOGwarn("Too small number of variants will be incorporated "
				"to IBS estimation, it is recommended to extend the range "
				"of --cormaf (%d variants included)\n", N_vrt-N_noIncCor);
		else if (!OPT_ENABLED(kinship) && !IS_ASSIGNED(cor))
			LOGwarn("Too small number of variants will be incorporated "
				"to empirical correlation estimation, it is recommended to "
				"extend the range of --cormaf (%d variants included)\n",
			N_vrt-N_noIncCor);
	} else
		LOG("Correlation computation with [%d] variants\n", N_vrt-N_noIncCor);

	/* When correlation matrix is not need, calc is also not need */
	if (!OPT_ENABLED(indep) && N_vrt) {
		if (IS_ASSIGNED(dosage)) _getGRM(N_corSamp, N_vrt, Ra_ppp, getIO()->getDosage(), getCorr());
		else _getGRM(N_corSamp, N_vrt, Ra_ppp, getIO()->getGenotype(), getCorr());
	}  else {
		/* Force --indep */
		if (!IS_ASSIGNED(indep)) {
			OPTION().assign("indep", "1", 1);
			OPTION().FORCE_OPT_NUMBER(indep);

			wsMat Ra_corr = getCorr();
			for (wsUint i=0 ; i<N_corSamp ; i++) {
				memset(Ra_corr[i], 0x00, sizeof(wsReal)*N_corSamp);
				Ra_corr[i][i] = W1;
			}
		}
	}

//	sseFree(Ra_corrMask);
//	sseFree(Na_corrMask);

	/* Export correlation matrix */
	_export();
	if (OPT_ENABLED(founderonly))
		halt("Cannot perform subsequent analysis due to --founderonly");
}

wsReal** cCorrAnalysis::getCorr() { return Cp_corr->get(); }

wsReal _corMedianComplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
//	wsReal	sumXY	= W0;
	wsUint	N		= 0;
	char	B_dalc	= 0;

	if (Ra_cors == NULL) {
		sseMalloc(Ra_cors, wsReal, N_SNP);
		B_dalc = 1;
	}

	wsReal R_ret	= W0; /* 0 in default */

	/* Do NOT consider missing, because the data is complete */
	if (Ra_cm) for (wsUint k=0 ; k<N_SNP ; k++) {
		/* Do not consider this marker if it should be filtered */
		if (Ra_cm[k] == W0) continue;

		Ra_cors[N++] = (Ra_v1[k] * Ra_v2[k]);
	} else for (N=0 ; N<N_SNP ; N++)
		Ra_cors[N] = (Ra_v1[N] * Ra_v2[N]);

	if (N) {
// 		qsort(Ra_cors, N, sizeof(wsReal), sort_real_unorder);
// 		/* Get the median */
// 		if ((N%2) == 0)
// 			sumXY = (Ra_cors[N>>1] + Ra_cors[(N>>1)+1]) / W2; 
// 		else
// 			sumXY = Ra_cors[N>>1];
// 		R_ret = sumXY;

		R_ret = quickmedian(Ra_cors, N);
	}

	/* Store included number of samples */
	if (Np_inc) *Np_inc = N; 

	/* FIXME : This part should SSE'd */
	if (B_dalc)
		sseFree(Ra_cors);

	return R_ret;
}

wsReal _corMeanComplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	wsReal sumXY	= W0;
	wsUint N		= 0;
	wsUint N_med	= getMed(N_SNP);
	wsReal R_ret	= W0;

	if (Ra_cm) {
#ifdef USE_SSE
	corrsse_t sse_sA = corrsseSet(0.0);
	corrsse_t sse_sumXY = corrsseSet(0.0);
	for (wsUint k=0 ; k<N_med ; k+=corrsseJmp) {
		corrsse_t *sse_mI = (corrsse_t *)(Ra_v1+k);
		corrsse_t *sse_mJ = (corrsse_t *)(Ra_v2+k);
		corrsse_t *sse_cm = (corrsse_t *)(Ra_cm+k);

		/* Modify sse_A to apply the filter */
		corrsse_t sse_tmp = corrsseSet(1.0);			/* tmp = 1 */
		sse_sA = corrsseAdd(sse_sA, corrsseAnd(sse_tmp, *sse_cm));

		corrsse_t sse_tI = corrsseAnd(*sse_mI, *sse_cm); // ti=mI if both 1
		corrsse_t sse_tJ = corrsseAnd(*sse_mJ, *sse_cm);

		sse_sumXY = corrsseAdd(sse_sumXY, corrsseMul(sse_tI, sse_tJ));
	}
	corrsseUsum(sse_sA, N);
	corrsseSum(sse_sumXY, sumXY);
	for (wsUint k=N_med ; k<N_SNP ; k++) {
		/* Do not consider this marker if it should be filtered */
		if (Ra_cm[k] == W0) continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		sumXY += Ra_v1[k] * Ra_v2[k];
		N++;
	}
#else
	for (wsUint k=0 ; k<N_SNP ; k++) {
		/* Do not consider this marker if it should be filtered */
		if (Ra_cm[k] == W0) continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
		sumXY += Ra_v1[k] * Ra_v2[k];
		N++;
	}
#endif
	} else {
#ifdef USE_SSE
		corrsse_t sse_sumXY = corrsseSet(0.0);
		for (wsUint k=0 ; k<N_med ; k+=corrsseJmp) {
			corrsse_t *sse_mI = (corrsse_t *)(Ra_v1+k);
			corrsse_t *sse_mJ = (corrsse_t *)(Ra_v2+k);

			sse_sumXY = corrsseAdd(sse_sumXY, corrsseMul(*sse_mI, *sse_mJ));
		}
		N = N_med;
		corrsseSum(sse_sumXY, sumXY);
		for (wsUint k=N_med ; k<N_SNP ; k++) {
			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
			sumXY += Ra_v1[k] * Ra_v2[k];
			N++;
		}
#else
		for (wsUint k=0 ; k<N_SNP ; k++) {
			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			sumXY += Ra_v1[k] * Ra_v2[k];
			N++;
		}
#endif
	}
	if (N)
		R_ret = sumXY/(wsReal)N;
	if (Np_inc) *Np_inc = N;

	return R_ret;
} 

void _corTest()
{
#define N_CORTEST 100
	/* Generate small matrix to calculate correlation */
}



wsReal _corPearsonIncomplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_obs,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	wsReal sumXY	= 0;
	wsReal sumX		= 0;
	wsReal sumX2	= 0;
	wsReal sumY		= 0;
	wsReal sumY2	= 0;
	wsUint N		= 0;
	//wsUint N_med	= getMed(N_SNP);
	wsReal R_ret	= W0;

	if (Ra_cm) {
		for (wsUint k=0 ; k<N_obs ; k++) {
			/* Do not consider this marker if it should be filtered */
			if (Ra_cm[k] == W0) continue;

			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			if (isAvailableReal(Ra_v1[k]) && isAvailableReal(Ra_v2[k])) {
				sumX += Ra_v1[k];
				sumX2 += Ra_v1[k] * Ra_v1[k];
				sumY += Ra_v2[k];
				sumY2 += Ra_v2[k] * Ra_v2[k];
				sumXY += Ra_v1[k] * Ra_v2[k];
				N++;
			}
		}
	} else {
		wsUint N_med = 0;

#ifdef USE_SSE
		corrsse_t x_sumx = corrsseSet(0.0);
		corrsse_t x_sumy = corrsseSet(0.0);
		corrsse_t x_sumx2 = corrsseSet(0.0);
		corrsse_t x_sumy2 = corrsseSet(0.0);
		corrsse_t x_sumxy = corrsseSet(0.0);
		corrsse_t x_1c = corrsseSet(1.0);
		N_med = getMed(N_obs);
		corrsse_t x_na = corrsseSet(WISARD_NA_REAL);
		corrsse_t x_n = corrsseSet(0.0);
		for (wsUint k=0 ; k<N_med ; k+=sseJmp) {
			corrsse_t *x_1 = (corrsse_t *)(Ra_v1 + k);
			corrsse_t *x_2 = (corrsse_t *)(Ra_v2 + k);
			corrsse_t x_mask = corrsseAnd(corrsseNeq(*x_1, x_na), corrsseNeq(*x_2, x_na));
			x_n = corrsseAdd(x_n, corrsseAnd(x_1c, x_mask));
			corrsse_t x_v1 = corrsseAnd(*x_1, x_mask);
			corrsse_t x_v2 = corrsseAnd(*x_2, x_mask);
			x_sumx = corrsseAdd(x_sumx, x_v1);
			x_sumy = corrsseAdd(x_sumy, x_v2);
			x_sumxy = corrsseAdd(x_sumxy, corrsseMul(x_v1, x_v2));
			x_sumx2 = corrsseAdd(x_sumx2, corrsseMul(x_v1, x_v1));
			x_sumy2 = corrsseAdd(x_sumy2, corrsseMul(x_v2, x_v2));
		}
		corrsseUsum(x_n, N);
		corrsseSum(x_sumx, sumX);
		corrsseSum(x_sumy, sumY);
		corrsseSum(x_sumxy, sumXY);
		corrsseSum(x_sumx2, sumX2);
		corrsseSum(x_sumy2, sumY2);
#endif
		for (wsUint k=N_med ; k<N_obs ; k++) {
			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			if (isAvailableReal(Ra_v1[k]) && isAvailableReal(Ra_v2[k])) {
				sumX += Ra_v1[k];
				sumX2 += Ra_v1[k] * Ra_v1[k];
				sumY += Ra_v2[k];
				sumY2 += Ra_v2[k] * Ra_v2[k];
				sumXY += Ra_v1[k] * Ra_v2[k];
				N++;
			}
		}
	}
	if (N)
		R_ret = ((wsReal)N*sumXY - sumX*sumY)
		/ (sqrt((wsReal)N*sumX2 - sumX*sumX) *
		sqrt((wsReal)N*sumY2 - sumY*sumY));
	if (Np_inc) *Np_inc = N;

	return R_ret;
}

wsReal corPearsonIncomplete(wsReal *Ra_v1, wsReal *Ra_v2, wsUint N_obs,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsRealCst *Ra_cm/*=NULL*/)
{
	wsReal sumXY	= 0;
	wsReal sumX		= 0;
	wsReal sumX2	= 0;
	wsReal sumY		= 0;
	wsReal sumY2	= 0;
	wsUint N		= 0;
	//wsUint N_med	= getMed(N_SNP);
	wsReal R_ret	= W0;

	if (Ra_cm) {
		for (wsUint k=0 ; k<N_obs ; k++) {
			/* Do not consider this marker if it should be filtered */
			if (Ra_cm[k] == W0) continue;

			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			if (isAvailableReal(Ra_v1[k]) && isAvailableReal(Ra_v2[k])) {
				sumX += Ra_v1[k];
				sumX2 += Ra_v1[k] * Ra_v1[k];
				sumY += Ra_v2[k];
				sumY2 += Ra_v2[k] * Ra_v2[k];
				sumXY += Ra_v1[k] * Ra_v2[k];
				N++;
			}
		}
	} else {
		wsUint N_med = 0;

#ifdef USE_SSE
		sse_t x_sumx = sseSet(0.0);
		sse_t x_sumy = sseSet(0.0);
		sse_t x_sumx2 = sseSet(0.0);
		sse_t x_sumy2 = sseSet(0.0);
		sse_t x_sumxy = sseSet(0.0);
		sse_t x_1c = sseSet(1.0);
		N_med = getMed(N_obs);
		sse_t x_na = sseSet(WISARD_NA_REAL);
		sse_t x_n = sseSet(0.0);
		for (wsUint k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *x_1 = (sse_t *)(Ra_v1 + k);
			sse_t *x_2 = (sse_t *)(Ra_v2 + k);
			sse_t x_mask = sseAnd(sseNeq(*x_1, x_na), sseNeq(*x_2, x_na));
			x_n = sseAdd(x_n, sseAnd(x_1c, x_mask));
			sse_t x_v1 = sseAnd(*x_1, x_mask);
			sse_t x_v2 = sseAnd(*x_2, x_mask);
			x_sumx = sseAdd(x_sumx, x_v1);
			x_sumy = sseAdd(x_sumy, x_v2);
			x_sumxy = sseAdd(x_sumxy, sseMul(x_v1, x_v2));
			x_sumx2 = sseAdd(x_sumx2, sseMul(x_v1, x_v1));
			x_sumy2 = sseAdd(x_sumy2, sseMul(x_v2, x_v2));
		}
		sseUsum(x_n, N);
		sseSum(x_sumx, sumX);
		sseSum(x_sumy, sumY);
		sseSum(x_sumxy, sumXY);
		sseSum(x_sumx2, sumX2);
		sseSum(x_sumy2, sumY2);
#endif
		for (wsUint k=N_med ; k<N_obs ; k++) {
			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			if (isAvailableReal(Ra_v1[k]) && isAvailableReal(Ra_v2[k])) {
				sumX += Ra_v1[k];
				sumX2 += Ra_v1[k] * Ra_v1[k];
				sumY += Ra_v2[k];
				sumY2 += Ra_v2[k] * Ra_v2[k];
				sumXY += Ra_v1[k] * Ra_v2[k];
				N++;
			}
		}
	}
	if (N)
		R_ret = ((wsReal)N*sumXY - sumX*sumY)
		/ (sqrt((wsReal)N*sumX2 - sumX*sumX) *
		sqrt((wsReal)N*sumY2 - sumY*sumY));
	if (Np_inc) *Np_inc = N;

	return R_ret;
}

wsReal _corPearsonIncomplete(char *Na_v1, char *Na_v2, wsUint N_obs,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, char *Na_cm/*=NULL*/)
{
	wsReal sumXY	= 0;
	wsReal sumX		= 0;
	wsReal sumX2	= 0;
	wsReal sumY		= 0;
	wsReal sumY2	= 0;
	wsUint N		= 0;
	//wsUint N_med	= getMed(N_SNP);
	wsReal R_ret	= W0;

	if (Na_cm) {
		wsUint	N_med = 0;

		if (OK(SSSE3)) {
			N_med = get16(N_obs);

			/* Valuemask to find out which is missing */
			__m128i sse_m1		= _mm_set1_epi8(-1);
			__m128i sse_1		= _mm_set1_epi8(1);

			__m128i sse_XY		= sse128zero; /* 8 x 16bit */
			__m128i sse_X2		= sse128zero; /* 8 x 16bit */
			__m128i sse_Y2		= sse128zero; /* 8 x 16bit */
			__m128i sse_X		= sse128zero; /* 16 x 8bit */
			__m128i sse_Y		= sse128zero; /* 16 x 8bit */
			__m128i sse_N		= sse128zero; /* 16 x 8bit */

			/* Process 16 genotypes in same time */
			for (wsUint i=0,j=0 ; i<N_med ; i+=16) {
				__m128i *p1 = (__m128i *)(Na_v1 + i);
				__m128i *p2 = (__m128i *)(Na_v2 + i);
				__m128i *ma = (__m128i *)(Na_cm + i);
				/* Make unified mask */
				__m128i fi = sse128and(*ma, sse128and(
					_mm_cmpgt_epi8(*p1, sse_m1),
					_mm_cmpgt_epi8(*p2, sse_m1)));

				/* Make filtered value */
				__m128i v1 = sse128and(*p1, fi);
				__m128i v2 = sse128and(*p2, fi);
				/* Set N */
				sse_N = _mm_add_epi8(sse_N, sse128and(fi, sse_1));

				/* Compute */
				sse_XY = _mm_add_epi16(sse_XY, _mm_maddubs_epi16(v1, v2)); //SSSE3
				sse_X2 = _mm_add_epi16(sse_X2, _mm_maddubs_epi16(v1, v1)); //SSSE3
				sse_Y2 = _mm_add_epi16(sse_Y2, _mm_maddubs_epi16(v2, v2)); //SSSE3
				sse_X = _mm_add_epi8(sse_X, v1); //SSSE3
				sse_Y = _mm_add_epi8(sse_Y, v2); //SSSE3

				// 				wsUint Q = 0;
				// 				for (wsUint k=i,l=0 ; l<16 ; l++,k++) {
				// 					if (isMissing(Na_1[k]) || isMissing(Na_2[k])) continue;
				// 					Q += Na_1[k] * Na_2[k];
				// 				}
				/* Avoid overflow */
				if (++j == 32) {
					sse8SHORTasum(sse_XY, sumXY);
					sse8SHORTasum(sse_X2, sumX2);
					sse8SHORTasum(sse_Y2, sumY2);
					sse16CHARasum(sse_X, sumX);
					sse16CHARasum(sse_Y, sumY);
					sse16CHARasum(sse_N, N);

					sse_XY	= sse128zero; //SSE2
					sse_X2	= sse128zero; //SSE2
					sse_Y2	= sse128zero; //SSE2
					sse_X	= sse128zero; //SSE2
					sse_Y	= sse128zero; //SSE2
					sse_N	= sse128zero; //SSE2
					j		= 0;
				}
			}
			sse8SHORTasum(sse_XY, sumXY);
			sse8SHORTasum(sse_X2, sumX2);
			sse8SHORTasum(sse_Y2, sumY2);
			sse16CHARasum(sse_X, sumX);
			sse16CHARasum(sse_Y, sumY);
			sse16CHARasum(sse_N, N);
		}

		/* Summation using CPU */
		for (wsUint i=N_med ; i<N_obs ; i++) {
			if (isMissing(Na_v1[i]) || isMissing(Na_v2[i])) continue;

			sumXY	+= Na_v1[i] * Na_v2[i];
			sumX2	+= Na_v1[i] * Na_v1[i];
			sumY2	+= Na_v2[i] * Na_v2[i];
			sumX	+= Na_v1[i];
			sumY	+= Na_v2[i];
			N++;
		}
	} else {
		wsUint	N_med = 0;

		if (OK(SSSE3)) {
			N_med = get16(N_obs);
			/* Valuemask to find out which is missing */
			__m128i sse_m1		= _mm_set_epi8(-1, -1, -1, -1, -1, -1, -1, -1,
				-1, -1, -1, -1, -1, -1, -1, -1);
			__m128i sse_1		= _mm_set_epi8(1, 1, 1, 1, 1, 1, 1, 1,
				1, 1, 1, 1, 1, 1, 1, 1);
			__m128i sse_XY		= _mm_setzero_si128(); /* 8 x 16bit */
			__m128i sse_X2		= _mm_setzero_si128(); /* 8 x 16bit */
			__m128i sse_Y2		= _mm_setzero_si128(); /* 8 x 16bit */
			__m128i sse_X		= _mm_setzero_si128(); /* 16 x 8bit */
			__m128i sse_Y		= _mm_setzero_si128(); /* 16 x 8bit */
			__m128i sse_N		= _mm_setzero_si128(); /* 16 x 8bit */

			/* Process 16 genotypes in same time */
			for (wsUint i=0,j=0 ; i<N_med ; i+=16) {
				__m128i *p1 = (__m128i *)(Na_v1 + i);
				__m128i *p2 = (__m128i *)(Na_v2 + i);
				/* Make unified mask */
				__m128i fi = _mm_and_si128(
					_mm_cmpgt_epi8(*p1, sse_m1),
					_mm_cmpgt_epi8(*p2, sse_m1));

				/* Make filtered value */
				__m128i v1 = _mm_and_si128(*p1, fi);
				__m128i v2 = _mm_and_si128(*p2, fi);
				/* Set N */
				sse_N = _mm_add_epi8(sse_N, _mm_and_si128(fi, sse_1));

				/* Compute */
				sse_XY = _mm_add_epi16(sse_XY, _mm_maddubs_epi16(v1, v2)); //SSSE3
				sse_X2 = _mm_add_epi16(sse_X2, _mm_maddubs_epi16(v1, v1)); //SSSE3
				sse_Y2 = _mm_add_epi16(sse_Y2, _mm_maddubs_epi16(v2, v2)); //SSSE3
				sse_X = _mm_add_epi8(sse_X, v1); //SSSE3
				sse_Y = _mm_add_epi8(sse_Y, v2); //SSSE3

// 				wsUint Q = 0;
// 				for (wsUint k=i,l=0 ; l<16 ; l++,k++) {
// 					if (isMissing(Na_1[k]) || isMissing(Na_2[k])) continue;
// 					Q += Na_1[k] * Na_2[k];
// 				}
				/* Avoid overflow */
				if (++j == 32) {
					sse8SHORTasum(sse_XY, sumXY);
					sse8SHORTasum(sse_X2, sumX2);
					sse8SHORTasum(sse_Y2, sumY2);
					sse16CHARasum(sse_X, sumX);
					sse16CHARasum(sse_Y, sumY);
					sse16CHARasum(sse_N, N);

					sse_XY	= _mm_setzero_si128(); //SSE2
					sse_X2	= _mm_setzero_si128(); //SSE2
					sse_Y2	= _mm_setzero_si128(); //SSE2
					sse_X	= _mm_setzero_si128(); //SSE2
					sse_Y	= _mm_setzero_si128(); //SSE2
					sse_N	= _mm_setzero_si128(); //SSE2
					j		= 0;
				}
			}
			sse8SHORTasum(sse_XY, sumXY);
			sse8SHORTasum(sse_X2, sumX2);
			sse8SHORTasum(sse_Y2, sumY2);
			sse16CHARasum(sse_X, sumX);
			sse16CHARasum(sse_Y, sumY);
			sse16CHARasum(sse_N, N);
		}

		/* Summation using CPU */
		for (wsUint i=N_med ; i<N_obs ; i++) {
			if (isMissing(Na_v1[i]) || isMissing(Na_v2[i])) continue;

			sumXY	+= Na_v1[i] * Na_v2[i];
			sumX2	+= Na_v1[i] * Na_v1[i];
			sumY2	+= Na_v2[i] * Na_v2[i];
			sumX	+= Na_v1[i];
			sumY	+= Na_v2[i];
			N++;
		}
	}
	if (N)
		R_ret = ((wsReal)N*sumXY - sumX*sumY)
			/ (sqrt((wsReal)N*sumX2 - sumX*sumX) *
			sqrt((wsReal)N*sumY2 - sumY*sumY));
	if (Np_inc) *Np_inc = N;

	return R_ret;
}

inline wsReal _computeCorr(wsUint N_col, wsReal *Ra_1, wsReal *Ra_2,
	wsReal R_Xi, wsReal R_X2i, wsReal R_Xj, wsReal R_X2j)
{
	wsReal R_col = (wsReal)N_col;
	return (R_col*(sseVV(Ra_1, N_col, Ra_2) - R_Xi*R_Xj) /
		sqrt((R_col*R_X2i - SQR(R_Xi)) * (R_col*R_X2j - SQR(R_Xj))));
}

#ifdef USE_SSE4

typedef union _w128 {
	__m128i s;
	unsigned int i[4];
} w128;

/* Reduce sixteen 8-bit integers into eight 16-bit integers by their adjacent sum */
#define _mm_rdc_8to16(s,Z)	_mm_add_epi16(_mm_cvtepu8_epi16(s), _mm_unpackhi_epi8((s), (Z)))
/* Reduce eight 16-bit integers into four 32-bit integers by their adjacent sum */
#define _mm_rdc_16to32(s,Z)	_mm_add_epi32(_mm_cvtepu16_epi32(s), _mm_unpackhi_epi16((s), (Z)))
#define _mm_part_ing_sum(part,zero,dest) { __m128i s16 = _mm_rdc_8to16(part, zero); \
	dest = _mm_add_epi32(dest, _mm_rdc_16to32(s16, zero)); \
	part = _mm_setzero_si128(); }
#define _mm_part_last_sum(part,zero,accu,dest) { __m128i s16 = _mm_rdc_8to16(part, zero); \
	dest = _mm_add_epi32(accu, _mm_rdc_16to32(s16, zero)); }

/*
 * Pearson's correlation coefficient
 *  - INPUT : char
 *  - DATA  : complete
 *  - RANGE : Pair-wise level
 */

inline __m128i sseGenoMul(__m128i aNew, __m128i bNew) {
#define MASK (char)0x80
	/* Separate LO part */
	__m128i Z = _mm_setzero_si128();
	__m128i Rlo = _mm_mullo_epi16(_mm_cvtepu8_epi16(aNew), _mm_cvtepu8_epi16(bNew));
	/* Separate HI part */
	__m128i Rhi = _mm_mullo_epi16(_mm_unpackhi_epi8(aNew,Z), _mm_unpackhi_epi8(bNew,Z));
	/* Combin & ret */
	__m128i maskLo = _mm_set_epi8(MASK, MASK, MASK, MASK, MASK, MASK, MASK, MASK, 14, 12, 10, 8, 6, 4, 2, 0);
	__m128i maskHi = _mm_set_epi8(14, 12, 10, 8, 6, 4, 2, 0, MASK, MASK, MASK, MASK, MASK, MASK, MASK, MASK);
	return _mm_or_si128(_mm_shuffle_epi8(Rlo, maskLo), _mm_shuffle_epi8(Rhi, maskHi));
#undef MASK
}

union char32
{
	WISARD_pc i[2];
	unsigned char v[16];
	__m128i s;
};

inline __m128i sseGenoMulIncomplete(__m128i a, __m128i b, int *Np_avail=NULL)
{
	__m128i mN = _mm_set_epi8(-9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9, -9);
	/* Apply mask to both value */
	__m128i aMask = _mm_cmpgt_epi8(a, mN);
	__m128i bMask = _mm_cmpgt_epi8(b, mN);
	__m128i aNew = _mm_and_si128(a, aMask);
	__m128i bNew = _mm_and_si128(b, bMask);
	char32 uMask;
	uMask.s = _mm_and_si128(aMask, bMask);
	if (Np_avail) {
		*Np_avail = 0;
#if !defined(_WIN32) || defined(_WIN32) && defined(_M_X64)
		for (wsUint i=0 ; i<2 ; i++)
#else
		for (wsUint i=0 ; i<4 ; i++)
#endif
			*Np_avail += (int)(WISARD_POPCNT(uMask.i[i]))>>3;
	}
//	printf("%d available\n", N_avail);

	return sseGenoMul(aNew, bNew);
}

inline __m128i sseGenoMulCompl(__m128i a, __m128i b) {
	return sseGenoMul(a, b);
}

inline wsReal sseVVchar(char *Na_1, wsUint N_col, char *Na_2)
{
	wsUint	N_med = 0;
	wsUint	N_ret = 0;
#ifdef USE_SSE
	w128	W_sse;
	memset(&W_sse, 0x00, sizeof(w128));
	N_med = get16(N_col);
	__m128i sse_part	= _mm_setzero_si128();
	__m128i sse_XY		= _mm_setzero_si128();
	for (wsUint i=0 ; i<N_med ; i+=16) {
		__m128i *a = (__m128i *)(Na_1 + i);
		__m128i *b = (__m128i *)(Na_2 + i);
		sse_part = _mm_add_epi8(sse_part, sseGenoMulCompl(*a, *b));

		/* Avoid an overflow of partial sum */
		if ((i%512) == (512-16)) {
			/* Lower part */
			__m128i Z	= _mm_setzero_si128();
			/* As eight 16-bit */
			__m128i s16	= _mm_rdc_8to16(sse_part, Z);
			sse_XY		= _mm_add_epi32(sse_XY, _mm_rdc_16to32(s16, Z));
			sse_part	= _mm_setzero_si128();
		}
	}
	/* Lower part */
	__m128i Z	= _mm_setzero_si128();
	/* As eight 16-bit */
	__m128i s16	= _mm_rdc_8to16(sse_part, Z);
	W_sse.s		= _mm_add_epi32(sse_XY, _mm_rdc_16to32(s16, Z));
	N_ret += W_sse.i[0] + W_sse.i[1] + W_sse.i[2] + W_sse.i[3];
#endif
	/* Summation using CPU */
	for (wsUint i=N_med ; i<N_col ; i++)
		N_ret += Na_1[i] * Na_2[i];

	return (wsReal)N_ret;
}

inline wsReal sseVsumChar(char *Na_1, wsUint N_col, wsReal *Rp_sqSum)
{
	wsUint	N_med	= 0;
	wsUint	N_retX	= 0;
	wsUint	N_retXX	= 0;
#ifdef USE_SSE
	w128	W_sse;
	memset(&W_sse, 0x00, sizeof(w128));
	N_med = get16(N_col);
	__m128i sse_pXX	= _mm_setzero_si128();
	__m128i sse_pX	= _mm_setzero_si128();
	__m128i sse_XX	= _mm_setzero_si128();
	__m128i sse_X	= _mm_setzero_si128();
	for (wsUint i=0 ; i<N_med ; i+=16) {
		__m128i *a	= (__m128i *)(Na_1 + i);
		sse_pXX		= _mm_add_epi8(sse_pXX, sseGenoMulCompl(*a, *a));
		sse_pX		= _mm_add_epi8(sse_pX, *a);

		/* Avoid an overflow of partial sum */
		if ((i%512) == (512-16)) {
			__m128i Z	= _mm_setzero_si128();
			/* pXX -> XX */
			__m128i s16	= _mm_rdc_8to16(sse_pXX, Z);
			sse_XX		= _mm_add_epi32(sse_XX, _mm_rdc_16to32(s16, Z));
			sse_pXX		= _mm_setzero_si128();
			/* pX -> X */
			s16			= _mm_rdc_8to16(sse_pX, Z);
			sse_X		= _mm_add_epi32(sse_X, _mm_rdc_16to32(s16, Z));
			sse_pX		= _mm_setzero_si128();
		}
	}
	/* Lower part */
	__m128i Z	= _mm_setzero_si128();
	/* pXX -> XX */
	__m128i s16	= _mm_rdc_8to16(sse_pXX, Z);
	W_sse.s		= _mm_add_epi32(sse_XX, _mm_rdc_16to32(s16, Z));
	N_retXX		+= W_sse.i[0] + W_sse.i[1] + W_sse.i[2] + W_sse.i[3];
	/* pX -> X */
	s16			= _mm_rdc_8to16(sse_pX, Z);
	W_sse.s		= _mm_add_epi32(sse_X, _mm_rdc_16to32(s16, Z));
	N_retX		+= W_sse.i[0] + W_sse.i[1] + W_sse.i[2] + W_sse.i[3];
#endif
	/* Summation using CPU */
	for (wsUint i=N_med ; i<N_col ; i++) {
		N_retX += Na_1[i];
		N_retX += Na_1[i]*Na_1[i];
	}

	*Rp_sqSum = (wsReal)N_retXX;
	return (wsReal)N_retX;
}

inline wsReal _sseCorCharComplete(wsUint N_col, char *Na_1, char *Na_2,
	wsReal R_Xi, wsReal R_X2i, wsReal R_Xj, wsReal R_X2j)
{
	wsReal R_col = (wsReal)N_col;
	return (R_col*sseVVchar(Na_1, N_col, Na_2) - R_Xi*R_Xj) /
		sqrt((R_col*R_X2i - SQR(R_Xi)) * (R_col*R_X2j - SQR(R_Xj)));
}

/**************************** END OF SECTION ****************************/

/*
 * Pearson's correlation coefficient
 *  - INPUT : char
 *  - DATA  : incomplete
 *  - RANGE : Pair-wise level
 */

inline wsReal sseVVcharIncomplete(char *Na_1, wsUint N_col, char *Na_2,
	wsUint *Np_ret)
{
	wsUint	N_med = 0;
	wsUint	N_ret = 0;
	*Np_ret = 0;
#ifdef USE_SSE
	w128	W_sse;
	memset(&W_sse, 0x00, sizeof(w128));
	N_med = get16(N_col);
	__m128i sse_part	= _mm_setzero_si128();
	__m128i sse_XY		= _mm_setzero_si128();
	for (wsUint i=0 ; i<N_med ; i+=16) {
		int N_partialRet = 0;
		__m128i *a = (__m128i *)(Na_1 + i);
		__m128i *b = (__m128i *)(Na_2 + i);
		sse_part = _mm_add_epi8(sse_part, sseGenoMulIncomplete(*a, *b,
			&N_partialRet));
		/* Add the number of incorporated genotype */
		*Np_ret += N_partialRet;

		/* Avoid an overflow of partial sum */
		if ((i%512) == (512-16)) {
			/* Lower part */
			__m128i Z	= _mm_setzero_si128();
			/* As eight 16-bit */
			__m128i s16	= _mm_add_epi16(_mm_cvtepu8_epi16(sse_part),
				_mm_unpackhi_epi8(sse_part, Z));
			sse_XY		= _mm_add_epi32(sse_XY,
				_mm_add_epi32(
					_mm_cvtepu16_epi32(s16), _mm_unpackhi_epi16(s16, Z)
				));
			sse_part	= _mm_setzero_si128();
		}
	}
	/* Lower part */
	__m128i Z	= _mm_setzero_si128();
	/* As eight 16-bit */
	__m128i s16	= _mm_add_epi16(_mm_cvtepu8_epi16(sse_part),
		_mm_unpackhi_epi8(sse_part, Z));
	W_sse.s		= _mm_add_epi32(sse_XY,
		_mm_add_epi32(
			_mm_cvtepu16_epi32(s16), _mm_unpackhi_epi16(s16, Z)
		));
	N_ret += W_sse.i[0] + W_sse.i[1] + W_sse.i[2] + W_sse.i[3];
#endif
	/* Summation using CPU */
	for (wsUint i=N_med ; i<N_col ; i++)
		if (isAvailable(Na_1[i]) && isAvailable(Na_2[i])) {
			N_ret += Na_1[i] * Na_2[i];
			(*Np_ret)++;
		}

	return (wsReal)N_ret;
}

inline wsReal _sseCorCharIncomplete(wsUint N_col, char *Na_1, char *Na_2)
{
	wsUint	N_med	= 0;
	wsUint	N_retX	= 0;
	wsUint	N_retXX	= 0;
	wsUint	N_retY	= 0;
	wsUint	N_retYY	= 0;
	wsUint	N_retXY	= 0;
#ifdef USE_SSE
	w128	W_sse;
	memset(&W_sse, 0x00, sizeof(w128));
	N_med = get16(N_col);
	__m128i sse_pXX	= _mm_setzero_si128();
	__m128i sse_pX	= _mm_setzero_si128();
	__m128i sse_pYY	= _mm_setzero_si128();
	__m128i sse_pY	= _mm_setzero_si128();
	__m128i sse_pXY	= _mm_setzero_si128();
	__m128i sse_XX	= _mm_setzero_si128();
	__m128i sse_X	= _mm_setzero_si128();
	__m128i sse_YY	= _mm_setzero_si128();
	__m128i sse_Y	= _mm_setzero_si128();
	__m128i sse_XY	= _mm_setzero_si128();
	wsUint	N_ent	= 0;
	for (wsUint i=0 ; i<N_med ; i+=16) {
		__m128i *a	= (__m128i *)(Na_1 + i);
		__m128i *b	= (__m128i *)(Na_2 + i);
		int		N	= 0;
		sse_pXX		= _mm_add_epi8(sse_pXX, sseGenoMulIncomplete(*a, *a));
		sse_pX		= _mm_add_epi8(sse_pX, *a);
		sse_pYY		= _mm_add_epi8(sse_pYY, sseGenoMulIncomplete(*b, *b));
		sse_pY		= _mm_add_epi8(sse_pY, *b);
		sse_pXY		= _mm_add_epi8(sse_pXX, sseGenoMulIncomplete(*a, *b, &N));
		N_ent		+= N;

		/* Avoid an overflow of partial sum */
		if ((i%512) == (512-16)) {
			__m128i Z	= _mm_setzero_si128();
				
			_mm_part_ing_sum(sse_pXX, Z, sse_XX);	/* pXX -> XX */
			_mm_part_ing_sum(sse_pX , Z, sse_X );	/* pX  -> X  */
			_mm_part_ing_sum(sse_pYY, Z, sse_YY);	/* pYY -> YY */
			_mm_part_ing_sum(sse_pY , Z, sse_Y );	/* pY  -> Y  */
			_mm_part_ing_sum(sse_pXY, Z, sse_XY);	/* pXY -> XY */
		}
	}
	/* Lower part */
	__m128i Z	= _mm_setzero_si128();

	_mm_part_last_sum(sse_pXX, Z, sse_XX, W_sse.s);	/* pXX -> XX */
	N_retXX = W_sse.i[0]+W_sse.i[1]+W_sse.i[2]+W_sse.i[3];
	_mm_part_last_sum(sse_pX , Z, sse_X , W_sse.s);	/* pX  -> X  */
	N_retX  = W_sse.i[0]+W_sse.i[1]+W_sse.i[2]+W_sse.i[3];
	_mm_part_last_sum(sse_pYY, Z, sse_YY, W_sse.s);	/* pYY -> YY */
	N_retYY = W_sse.i[0]+W_sse.i[1]+W_sse.i[2]+W_sse.i[3];
	_mm_part_last_sum(sse_pY , Z, sse_Y , W_sse.s);	/* pY  -> Y  */
	N_retY  = W_sse.i[0]+W_sse.i[1]+W_sse.i[2]+W_sse.i[3];
	_mm_part_last_sum(sse_pXY, Z, sse_XY, W_sse.s);	/* pXY -> XY */
	N_retXY = W_sse.i[0]+W_sse.i[1]+W_sse.i[2]+W_sse.i[3];
#endif
	/* Summation using CPU */
	for (wsUint i=N_med ; i<N_col ; i++) {
		N_retX	+= Na_1[i];
		N_retXX	+= Na_1[i]*Na_1[i];
		N_retY	+= Na_2[i];
		N_retYY	+= Na_2[i]*Na_2[i];
		N_retXY	+= Na_1[i]*Na_2[i];
	}

	wsReal R_avail = (wsReal)N_ent;
	return (R_avail*(wsReal)N_retXY - (wsReal)N_retX*(wsReal)N_retY) /
		sqrt(
			(R_avail*(wsReal)N_retXX - SQR((wsReal)N_retX)) *
			(R_avail*(wsReal)N_retYY - SQR((wsReal)N_retY))
		);
}

/**************************** END OF SECTION ****************************/

/*
 * Pearson's correlation coefficient
 *  - INPUT : char
 *  - DATA  : incomplete
 *  - RANGE : Dataset level
 */

wsMat corPearson(char *Ba_isComplete, char **Na_data, wsUint N_samp,
	wsUint N_marker, char B_compl)
{
	/* N in the computation of correlation == # of snps == # of columns */
	//wsReal	R_mk	= (wsReal)N_marker;
	/* Temporal buffer for cachable results */
	wsReal	*Ra_X	= NULL;
	wsReal	*Ra_X2	= NULL;
	/* Allocate */
	wsMat	Ra_ret	= sseMatrix(N_samp, N_samp);
	sseMalloc(Ra_X, wsReal, N_samp);
	sseMalloc(Ra_X2, wsReal, N_samp);

	/* Calculate X and X^2 */
	for (wsUint i=0 ; i<N_samp ; i++)
		Ra_X[i] = sseVsumChar(Na_data[i], N_marker, Ra_X2 + i);

	if (B_compl) {
		/* For all samples */
		for (wsUint i=0 ; i<N_samp ; i++) {
			for (wsUint j=0 ; j<i ; j++)
				Ra_ret[j][i] = Ra_ret[i][j] = _sseCorCharComplete(N_marker,
					Na_data[i], Na_data[j], Ra_X[i], Ra_X2[i], Ra_X[j],
					Ra_X2[j]);

			/* Self-self correlation in Pearson's correlation coefficient is always 1 */
			Ra_ret[i][i] = W1;
		}
	} else {
		for (wsUint i=0 ; i<N_samp ; i++) {
			for (wsUint j=0 ; j<i ; j++)
				/* Use complete method if both sample have complete genotype */
				if (Ba_isComplete[i] && Ba_isComplete[j])
					Ra_ret[j][i] = Ra_ret[i][j] = _sseCorCharComplete(
						N_marker, Na_data[i], Na_data[j], Ra_X[i],
						Ra_X2[i], Ra_X[j], Ra_X2[j]);
				else
					/* Otherwise, using slow method */
					Ra_ret[j][i] = Ra_ret[i][j] = _sseCorCharIncomplete(N_marker,
						Na_data[i], Na_data[j]);

			/* Self-self correlation in Pearson's correlation coefficient is always 1 */
			Ra_ret[i][i] = W1;
		}
	}

	return Ra_ret;
}
#endif

/**************************** END OF SECTION ****************************/

wsSym pearson(char *Ba_isComplete, wsMat Ra_mat, wsUint N_row,
	wsUint N_col, char B_isComplete=1)
{
	/* N in the computation of correlation == # of snps == # of columns */
	//wsReal	R_col	= (wsReal)N_col;
	/* Temporal buffer for cachable results */
	wsReal	*Ra_X	= NULL;
	wsReal	*Ra_X2	= NULL;
	/* Allocate */
	wsSym	Ra_ret	= sseSymMat(N_row);
	sseMalloc(Ra_X, wsReal, N_row);
	sseMalloc(Ra_X2, wsReal, N_row);

	/* Calculate X and X^2 */
	for (wsUint i=0 ; i<N_row ; i++)
		Ra_X[i] = sseVsum(Ra_mat[i], N_col, Ra_X2 + i);

	if (B_isComplete) {
		/* For all samples */
		for (wsUint i=0 ; i<N_row ; i++) {
			for (wsUint j=0 ; j<i ; j++)
				Ra_ret[i][j] = _computeCorr(N_col, Ra_mat[i], Ra_mat[j],
					Ra_X[i], Ra_X2[i], Ra_X[j], Ra_X2[j]);

			/* Self-self correlation in Pearson's correlation coefficient is always 1 */
			Ra_ret[i][i] = W1;
		}
	} else {
		for (wsUint i=0 ; i<N_row ; i++) {
			for (wsUint j=0 ; j<i ; j++)
				/* Use complete method if both sample have complete genotype */
				if (Ba_isComplete[i] && Ba_isComplete[j])
					Ra_ret[i][j] = _computeCorr(N_col, Ra_mat[i], Ra_mat[j],
					Ra_X[i], Ra_X2[i], Ra_X[j], Ra_X2[j]);
				else
					/* Otherwise, using slow method */
					Ra_ret[i][j] = corPearsonIncomplete(Ra_mat[i], Ra_mat[j], N_col);

			/* Self-self correlation in Pearson's correlation coefficient is always 1 */
			Ra_ret[i][i] = W1;
		}
	}

	sseFree(Ra_X);
	sseFree(Ra_X2);
	return Ra_ret;
}

wsReal _corPearsonComplete(wsFloat *Na_v1, wsFloat *Na_v2, wsUint N_SNP,
   wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	wsFloat sumXY	= 0;
	wsFloat sumX	= 0;
	wsFloat sumX2	= 0;
	wsFloat	sumY	= 0;
	wsFloat sumY2	= 0;
	wsUint	N		= 0;
	//wsUint N_med	= getMed(N_SNP);
	wsFloat	R_ret	= 0;

	if (Ra_cm) {
		for (wsUint k=0 ; k<N_SNP ; k++) {
			/* Do not consider this marker if it should be filtered */
			if (Ra_cm[k] == W0) continue;

			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			sumX += Na_v1[k];
			sumX2 += Na_v1[k] * Na_v1[k];
			sumY += Na_v2[k];
			sumY2 += Na_v2[k] * Na_v2[k];
			sumXY += Na_v1[k] * Na_v2[k];
			N++;
		}
	} else {
		wsUint N_med = 0;

#ifdef USE_SSE
		corrsse_t x_sumx	= corrsseSet(0.0);
		corrsse_t x_sumy	= corrsseSet(0.0);
		corrsse_t x_sumx2	= corrsseSet(0.0);
		corrsse_t x_sumy2	= corrsseSet(0.0);
		corrsse_t x_sumxy	= corrsseSet(0.0);
		corrsse_t x_1c		= corrsseSet(1.0);
		N_med = corrGetMed(N_SNP);
		corrsse_t x_n = corrsseSet(0.0);
		for (wsUint k=0 ; k<N_med ; k+=corrsseJmp) {
			corrsse_t *x_1 = (corrsse_t *)(Na_v1 + k);
			corrsse_t *x_2 = (corrsse_t *)(Na_v2 + k);
			x_n = corrsseAdd(x_n, x_1c);
			x_sumx = corrsseAdd(x_sumx, *x_1);
			x_sumy = corrsseAdd(x_sumy, *x_2);
			x_sumxy = corrsseAdd(x_sumxy, corrsseMul(*x_1, *x_2));
			x_sumx2 = corrsseAdd(x_sumx2, corrsseMul(*x_1, *x_1));
			x_sumy2 = corrsseAdd(x_sumy2, corrsseMul(*x_2, *x_2));
		}
		corrsseUsum(x_n, N);
		corrsseSum(x_sumx, sumX);
		corrsseSum(x_sumy, sumY);
		corrsseSum(x_sumxy, sumXY);
		corrsseSum(x_sumx2, sumX2);
		corrsseSum(x_sumy2, sumY2);
#endif
		for (wsUint k=N_med ; k<N_SNP ; k++) {
			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			sumX += Na_v1[k];
			sumX2 += Na_v1[k] * Na_v1[k];
			sumY += Na_v2[k];
			sumY2 += Na_v2[k] * Na_v2[k];
			sumXY += Na_v1[k] * Na_v2[k];
			N++;
		}
	}
	if (N)
		R_ret = ((wsFloat)N*sumXY - sumX*sumY)
			/ (sqrt((wsFloat)N*sumX2 - sumX*sumX) *
			sqrt((wsFloat)N*sumY2 - sumY*sumY));
	if (Np_inc) *Np_inc = N;

	return R_ret;
}

wsReal corPearsonComplete(wsReal *Na_v1, wsReal *Na_v2, wsUint N_SNP,
   wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsRealCst *Ra_cm/*=NULL*/)
{
	wsReal sumXY	= 0;
	wsReal sumX		= 0;
	wsReal sumX2	= 0;
	wsReal sumY		= 0;
	wsReal sumY2	= 0;
	wsUint N		= 0;
	//wsUint N_med	= getMed(N_SNP);
	wsReal R_ret	= W0;

	if (Ra_cm) {
		for (wsUint k=0 ; k<N_SNP ; k++) {
			/* Do not consider this marker if it should be filtered */
			if (Ra_cm[k] == W0) continue;

			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			sumX += Na_v1[k];
			sumX2 += Na_v1[k] * Na_v1[k];
			sumY += Na_v2[k];
			sumY2 += Na_v2[k] * Na_v2[k];
			sumXY += Na_v1[k] * Na_v2[k];
			N++;
		}
	} else {
		wsUint N_med = 0;

#ifdef USE_SSE
		sse_t x_sumx = sseSet(0.0);
		sse_t x_sumy = sseSet(0.0);
		sse_t x_sumx2 = sseSet(0.0);
		sse_t x_sumy2 = sseSet(0.0);
		sse_t x_sumxy = sseSet(0.0);
		sse_t x_1c = sseSet(1.0);
		N_med = getMed(N_SNP);
		sse_t x_n = sseSet(0.0);
		for (wsUint k=0 ; k<N_med ; k+=sseJmp) {
			sse_t *x_1 = (sse_t *)(Na_v1 + k);
			sse_t *x_2 = (sse_t *)(Na_v2 + k);
			x_n = sseAdd(x_n, x_1c);
			x_sumx = sseAdd(x_sumx, *x_1);
			x_sumy = sseAdd(x_sumy, *x_2);
			x_sumxy = sseAdd(x_sumxy, sseMul(*x_1, *x_2));
			x_sumx2 = sseAdd(x_sumx2, sseMul(*x_1, *x_1));
			x_sumy2 = sseAdd(x_sumy2, sseMul(*x_2, *x_2));
		}
		sseUsum(x_n, N);
		sseSum(x_sumx, sumX);
		sseSum(x_sumy, sumY);
		sseSum(x_sumxy, sumXY);
		sseSum(x_sumx2, sumX2);
		sseSum(x_sumy2, sumY2);
#endif
		for (wsUint k=N_med ; k<N_SNP ; k++) {
			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			sumX += Na_v1[k];
			sumX2 += Na_v1[k] * Na_v1[k];
			sumY += Na_v2[k];
			sumY2 += Na_v2[k] * Na_v2[k];
			sumXY += Na_v1[k] * Na_v2[k];
			N++;
		}
	}
	if (N)
		R_ret = ((wsReal)N*sumXY - sumX*sumY)
		/ (sqrt((wsReal)N*sumX2 - sumX*sumX) *
		sqrt((wsReal)N*sumY2 - sumY*sumY));
	if (Np_inc) *Np_inc = N;

	return R_ret;
}

wsReal _corMeanComplete(char *Na_v1, char *Na_v2, wsUint N_SNP,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	wsUint sumXY	= 0;
	wsUint sumX		= 0;
	wsUint sumX2	= 0;
	wsUint sumY		= 0;
	wsUint sumY2	= 0;
	wsUint N		= 0;
	//wsUint N_med	= getMed(N_SNP);
	wsReal R_ret	= W0;

	if (Ra_cm) {
		for (wsUint k=0 ; k<N_SNP ; k++) {
			/* Do not consider this marker if it should be filtered */
			if (Ra_cm[k] == W0) continue;

			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			sumX += Na_v1[k];
			sumX2 += Na_v1[k] * Na_v1[k];
			sumY += Na_v2[k];
			sumY2 += Na_v2[k] * Na_v2[k];
			sumXY += Na_v1[k] * Na_v2[k];
			N++;
		}
	} else {
		for (wsUint k=0 ; k<N_SNP ; k++) {
			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			sumX += Na_v1[k];
			sumX2 += Na_v1[k] * Na_v1[k];
			sumY += Na_v2[k];
			sumY2 += Na_v2[k] * Na_v2[k];
			sumXY += Na_v1[k] * Na_v2[k];
			N++;
		}
	}
	if (N)
		R_ret = ((wsReal)N*(wsReal)sumXY - (wsReal)sumX*(wsReal)sumY)
			/ (sqrt((wsReal)N*(wsReal)sumX2 - (wsReal)sumX*(wsReal)sumX) *
			sqrt((wsReal)N*(wsReal)sumY2 - (wsReal)sumY*(wsReal)sumY));
	if (Np_inc) *Np_inc = N;

	return R_ret;
}

wsReal _corMedianIncomplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	wsUint N	= 0;

	if (Ra_cm) for (wsUint k=0 ; k<N_SNP ; k++) {
		/* Do not consider this marker if it should be filtered */
		if (Ra_cm[k] == W0) continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		if (isAvailableReal(Ra_v1[k]) && isAvailableReal(Ra_v2[k])) {
			Ra_cors[N] = (Ra_v1[k] * Ra_v2[k]);
			N++;
		}
	} else for (wsUint k=0 ; k<N_SNP ; k++) {
		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		if (isAvailableReal(Ra_v1[k]) && isAvailableReal(Ra_v2[k])) {
			Ra_cors[N] = (Ra_v1[k] * Ra_v2[k]);
			N++;
		}
	}

	wsReal R_ret = W0; /* 0 in default */
	if (N) {
		R_ret = quickmedian(Ra_cors, N);
// 		qsort(Ra_cors, N, sizeof(wsReal), sort_real_unorder);
// 		/* Get the median */
// 		if ((N%2) == 0)
// 			R_ret = (Ra_cors[N>>1] + Ra_cors[(N>>1)+1]) / W2; 
// 		else
// 			R_ret = Ra_cors[N>>1];
	}
	if (Np_inc) *Np_inc = N;

	return R_ret;
}

wsReal _corMeanIncomplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	wsReal sumXY	= W0;
	wsUint N		= 0;
	wsUint N_med	= 0;

#ifdef USE_SSE
	N_med	= getMed(N_SNP);
	corrsse_t sse_sA = corrsseSet(0.0);
	corrsse_t sse_sumXY = corrsseSet(0.0);

	if (Ra_cm) for (wsUint k=0 ; k<N_med ; k+=corrsseJmp) {
		corrsse_t sse_A;
		corrsse_t *sse_mI = (corrsse_t *)(Ra_v1+k);
		corrsse_t *sse_mJ = (corrsse_t *)(Ra_v2+k);
		corrsse_t *sse_cm = (corrsse_t *)(Ra_cm+k);

		corrsse_t sse_tmp = corrsseSet(WISARD_NA_REAL);	/* tmp = MISSING */
		sse_A = corrsseAnd(corrsseNeq(*sse_mI, sse_tmp), corrsseNeq(*sse_mJ, sse_tmp));
		/* Modify sse_A to apply the filter */
		sse_A = corrsseAnd(sse_A, *sse_cm);
		sse_tmp = corrsseSet(1.0);			/* tmp = 1 */
		sse_sA = corrsseAdd(sse_sA, corrsseAnd(sse_tmp, sse_A));

		corrsse_t sse_tI = corrsseAnd(*sse_mI, sse_A); // ti=mI if both 1
		corrsse_t sse_tJ = corrsseAnd(*sse_mJ, sse_A);

		sse_sumXY = corrsseAdd(sse_sumXY, corrsseMul(sse_tI, sse_tJ));
	} else for (wsUint k=0 ; k<N_med ; k+=corrsseJmp) {
		corrsse_t sse_A;
		corrsse_t *sse_mI = (corrsse_t *)(Ra_v1+k);
		corrsse_t *sse_mJ = (corrsse_t *)(Ra_v2+k);

		corrsse_t sse_tmp = corrsseSet(WISARD_NA_REAL);	/* tmp = MISSING */
		sse_A = corrsseAnd(corrsseNeq(*sse_mI, sse_tmp), corrsseNeq(*sse_mJ, sse_tmp));
		sse_tmp = corrsseSet(1.0);			/* tmp = 1 */
		sse_sA = corrsseAdd(sse_sA, corrsseAnd(sse_tmp, sse_A));

		corrsse_t sse_tI = corrsseAnd(*sse_mI, sse_A); // ti=mI if both 1
		corrsse_t sse_tJ = corrsseAnd(*sse_mJ, sse_A);

		sse_sumXY = corrsseAdd(sse_sumXY, corrsseMul(sse_tI, sse_tJ));
	}
	corrsseUsum(sse_sA, N);
	corrsseSum(sse_sumXY, sumXY);
#endif
	if (Ra_cm) for (wsUint k=N_med ; k<N_SNP ; k++) {
		/* Do not consider this marker if it should be filtered */
		if (Ra_cm[k] == W0) continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		if (isAvailableReal(Ra_v1[k]) && isAvailableReal(Ra_v2[k])) {
			sumXY += Ra_v1[k] * Ra_v2[k];
			N++;
		}
	} else for (wsUint k=N_med ; k<N_SNP ; k++) {
		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		if (isAvailableReal(Ra_v1[k]) && isAvailableReal(Ra_v2[k])) {
			sumXY += Ra_v1[k] * Ra_v2[k];
			N++;
		}
	}

	wsReal R_ret = W0; /* 0 in default */
	if (N)
		R_ret = sumXY/(wsReal)N;
	if (Np_inc) *Np_inc = N;

	return R_ret;
}

/* Assume all NAs are 0 so does not affect to the result */
wsReal _corMeanIncomplete2(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	char *Na_v1, char *Na_v2, wsUint *Np_inc/*=NULL*/, char *Na_cm/*=NULL*/)
{
	wsReal	sumXY	= W0;
	wsUint	N		= 0;
	wsUint	N_med	= 0;
	wsUint	N_med16	= 0;
	char*	Na_ccm	= new char[N_SNP];
	memcpy(Na_ccm, Na_cm, sizeof(char)*N_SNP);
	
	/* Update ccm */
	for (wsUint i=0 ; i<N_SNP ; i++) {
		unsigned char U = isMissing(Na_v1[i]) | isMissing(Na_v2[i]);
		unsigned char UU = U-1;
		Na_ccm[i] &= UU;
	}

#ifdef USE_SSE
	N_med	= corrGetMed(N_SNP);
	N_med16	= get16(N_SNP);
	corrsse_t sse_sumXY = corrsseSet(0.0);

	for (wsUint k=0 ; k<N_med ; k+=corrsseJmp) {
		corrsse_t *sse_mI = (corrsse_t *)(Ra_v1+k);
		corrsse_t *sse_mJ = (corrsse_t *)(Ra_v2+k);

		sse_sumXY = corrsseAdd(sse_sumXY, corrsseMul(*sse_mI, *sse_mJ));
	}
	corrsseSum(sse_sumXY, sumXY);

	__m128i sse_A = _mm_set1_epi8(0);
	__m128i sse_1 = _mm_set1_epi8(1);
	for (wsUint k=0 ; k<N_med16 ; k+=16) {
		__m128i *sse_cm = (__m128i *)(Na_cm + k);
		sse_A = _mm_add_epi8(sse128and(*sse_cm, sse_1), sse_A);
		if (k%(16<<6)) {
			sse16CHARasum(sse_A, N);
			sse_A = _mm_set1_epi8(0);
		}
	}
	sse16CHARasum(sse_A, N);
#endif
	for (wsUint k=N_med ; k<N_SNP ; k++) {
		/* Do not consider this marker if it should be filtered */
		if (!Na_cm[k]) continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		sumXY += Ra_v1[k] * Ra_v2[k];
		N++;
	}
	for (wsUint k=N_med16 ; k<N_SNP ; k++)
		if (Na_ccm[k]) N++;

	wsReal R_ret = W0; /* 0 in default */
	if (N)
		R_ret = sumXY/(wsReal)N;
	if (Np_inc) *Np_inc = N;

	return R_ret;
}

wsReal _corMeanIncomplete(char *Na_v1, char *Na_v2, wsUint N_SNP,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	wsUint sumXY	= 0;
	wsUint sumX		= 0;
	wsUint sumX2	= 0;
	wsUint sumY		= 0;
	wsUint sumY2	= 0;
	wsUint N		= 0;
	//wsUint N_med	= getMed(N_SNP);
	wsReal R_ret	= W0;

/*
  0 1 2 0
0 0 0 0 0
1 0 1 2 0
2 0 2 4 0
0 0 0 0 0

v1<<v2

00 11 11
01 11 10
10 11 01
11 11 00
         00000000 00000001 00000010 00000000
00000000
00000001
00000010
11110111

V^(V&0xf)

01 01 001
01 10 010
10 01 010
10 10 100
int v1 = i[j];
int v2 = i[k];
int m1 = m(v1);
int m2 = m(v2);
int x1 = (v2<<(m1-1));


 */
#define m(v) (v)&(((v)^(((v)>>4)&0xf))&0xf)

	if (Ra_cm) {
		for (wsUint k=0 ; k<N_SNP ; k++) {
			/* Do not consider this marker if it should be filtered */
			if (Ra_cm[k] == W0) continue;

			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			if (isAvailable(Na_v1[k]) && isAvailable(Na_v2[k])) {
				sumX	+= Na_v1[k];
				sumX2	+= Na_v1[k] * Na_v1[k];
				sumY	+= Na_v2[k];
				sumY2	+= Na_v2[k] * Na_v2[k];
				sumXY	+= Na_v1[k] * Na_v2[k];
				N++;
			}
		}
	} else {
		for (wsUint k=0 ; k<N_SNP ; k++) {
			/* use="pair", i�� j ���� �� �� k��° SNP�� ���� gxenotype�� �־�߸� ��� */
			if (isAvailable(Na_v1[k]) && isAvailable(Na_v2[k])) {
				sumX += Na_v1[k];
				sumX2 += Na_v1[k] * Na_v1[k];
				sumY += Na_v2[k];
				sumY2 += Na_v2[k] * Na_v2[k];
				sumXY += Na_v1[k] * Na_v2[k];
				N++;
			}
		}
	}
	if (N)
		R_ret = ((wsReal)N*(wsReal)sumXY - (wsReal)sumX*(wsReal)sumY)
		/ (sqrt((wsReal)N*(wsReal)sumX2 - (wsReal)sumX*(wsReal)sumX) *
		sqrt((wsReal)N*(wsReal)sumY2 - (wsReal)sumY*(wsReal)sumY));
	if (Np_inc) *Np_inc = N;

	return R_ret;
}

wsReal _corTauIncomplete(char *Na_v1, char *Na_v2, wsUint N_SNP,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	wsUint N_dis = 0;
	wsUint N_con = 0;
	wsUint N = 0;
	char *Ba_t1 = NULL;
	char *Ba_t2 = NULL;
	wsCalloc(Ba_t1, char, N_SNP);
	wsCalloc(Ba_t2, char, N_SNP);

	if (Ra_cm) for (wsUint i=0 ; i<N_SNP ; i++) {
		if (isMissing(Na_v1[i]) || isMissing(Na_v2[i]) ||
			Ra_cm[i] == W0) continue;

		for (wsUint j=i+1 ; j<N_SNP ; j++) {
			if (isMissing(Na_v1[j]) || isMissing(Na_v2[j]) ||
				Ra_cm[i] == W0) continue;

			bool B_eq1 = Na_v1[i] == Na_v1[j];
			bool B_eq2 = Na_v2[i] == Na_v2[j];
			if (B_eq1 || B_eq2) {
				if (B_eq1) Ba_t1[i] = Ba_t1[j] = 1;
				if (B_eq2) Ba_t2[i] = Ba_t2[j] = 1;
				continue;
			}

			bool B_gt1 = Na_v1[i] > Na_v1[j];
			bool B_gt2 = Na_v2[i] > Na_v2[j];

			(B_gt1^B_gt2) ? N_dis++ : N_con++;
//			N++;
		}
	} else for (wsUint i=0 ; i<N_SNP ; i++) {
		if (isMissing(Na_v1[i]) || isMissing(Na_v2[i])) continue;

		for (wsUint j=i+1 ; j<N_SNP ; j++) {
			if (isMissing(Na_v1[j]) || isMissing(Na_v2[j])) continue;

			bool B_eq1 = Na_v1[i] == Na_v1[j];
			bool B_eq2 = Na_v2[i] == Na_v2[j];
			if (B_eq1 || B_eq2) {
				if (B_eq1) Ba_t1[i] = Ba_t1[j] = 1;
				if (B_eq2) Ba_t2[i] = Ba_t2[j] = 1;
				continue;
			}

			bool B_gt1 = Na_v1[i] > Na_v1[j];
			bool B_gt2 = Na_v2[i] > Na_v2[j];

			(B_gt1^B_gt2) ? N_dis++ : N_con++;
//			N++;
		}
	}

	wsUint N_t1[3] = { 0, };
	wsUint N_t2[3] = { 0, };
	for (wsUint i=0 ; i<N_SNP ; i++) {
		if (Ra_cm && Ra_cm[i] == W0) continue;
		if (Ba_t1[i]) N_t1[(int)Na_v1[i]]++;
		if (Ba_t2[i]) N_t2[(int)Na_v2[i]]++;
		N++;
	}
	N = N&1 ? N * ((N-1)>>1) : (N>>1) * (N-1);
	/* Calculate tie term */
	wsReal R_tm1 = (wsReal)(N_t1[0]*(N_t1[0]-1) + N_t1[1]*(N_t1[1]-1) +
		N_t1[2]*(N_t1[2]-1)) / W2;
	wsReal R_tm2 = (wsReal)(N_t2[0]*(N_t2[0]-1) + N_t2[1]*(N_t2[1]-1) +
		N_t2[2]*(N_t2[2]-1)) / W2;
	wsReal R_denom = ((wsReal)N-R_tm1)*((wsReal)N-R_tm2);
	if (Np_inc) *Np_inc = N;

	DEALLOC(Ba_t1);
	DEALLOC(Ba_t2);

	return R_denom<=W0 ? WISARD_NAN :
		((wsReal)N_con - (wsReal)N_dis) / sqrt(((wsReal)N-R_tm1)*((wsReal)N-R_tm2));
}

// wsReal _corTauComplete(char *Na_v1, char *Na_v2, wsUint N_SNP,
// 	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
// {
// 	wsUint N_dis = 0;
// 	wsUint N_con = 0;
// 	wsUint N = 0;
// 	char *Ba_t1 = NULL;
// 	char *Ba_t2 = NULL;
// 	MULTI_CALLOC(Ba_t1, char, N_SNP);
// 	MULTI_CALLOC(Ba_t2, char, N_SNP);
// 
// 	if (Ra_cm) for (wsUint i=0 ; i<N_SNP ; i++) {
// 		if (isMissing(Na_v1[i]) || isMissing(Na_v2[i]) ||
// 			Ra_cm[i] == W0) continue;
// 
// 		for (wsUint j=i+1 ; j<N_SNP ; j++) {
// 			if (isMissing(Na_v1[j]) || isMissing(Na_v2[j]) ||
// 				Ra_cm[i] == W0) continue;
// 
// 			bool B_eq1 = Na_v1[i] == Na_v1[j];
// 			bool B_eq2 = Na_v2[i] == Na_v2[j];
// 			if (B_eq1 || B_eq2) {
// 				if (B_eq1) Ba_t1[i] = Ba_t1[j] = 1;
// 				if (B_eq2) Ba_t2[i] = Ba_t2[j] = 1;
// 				continue;
// 			}
// 
// 			bool B_gt1 = Na_v1[i] > Na_v1[j];
// 			bool B_gt2 = Na_v2[i] > Na_v2[j];
// 
// 			(B_gt1^B_gt2) ? N_dis++ : N_con++;
// 			//			N++;
// 		}
// 	} else for (wsUint i=0 ; i<N_SNP ; i++) {
// 		if (isMissing(Na_v1[i]) || isMissing(Na_v2[i])) continue;
// 
// 		for (wsUint j=i+1 ; j<N_SNP ; j++) {
// 			if (isMissing(Na_v1[j]) || isMissing(Na_v2[j])) continue;
// 
// 			bool B_eq1 = Na_v1[i] == Na_v1[j];
// 			bool B_eq2 = Na_v2[i] == Na_v2[j];
// 			if (B_eq1 || B_eq2) {
// 				if (B_eq1) Ba_t1[i] = Ba_t1[j] = 1;
// 				if (B_eq2) Ba_t2[i] = Ba_t2[j] = 1;
// 				continue;
// 			}
// 
// 			bool B_gt1 = Na_v1[i] > Na_v1[j];
// 			bool B_gt2 = Na_v2[i] > Na_v2[j];
// 
// 			(B_gt1^B_gt2) ? N_dis++ : N_con++;
// 			//			N++;
// 		}
// 	}
// 
// 	wsUint N_t1[3] = { 0, };
// 	wsUint N_t2[3] = { 0, };
// 	for (wsUint i=0 ; i<N_SNP ; i++) {
// 		if (Ra_cm && Ra_cm[i] == W0) continue;
// 		if (Ba_t1[i]) N_t1[(int)Na_v1[i]]++;
// 		if (Ba_t2[i]) N_t2[(int)Na_v2[i]]++;
// 		N++;
// 	}
// 	N = N&1 ? N * ((N-1)>>1) : (N>>1) * (N-1);
// 	/* Calculate tie term */
// 	wsReal R_tm1 = (wsReal)(N_t1[0]*(N_t1[0]-1) + N_t1[1]*(N_t1[1]-1) +
// 		N_t1[2]*(N_t1[2]-1)) / W2;
// 	wsReal R_tm2 = (wsReal)(N_t2[0]*(N_t2[0]-1) + N_t2[1]*(N_t2[1]-1) +
// 		N_t2[2]*(N_t2[2]-1)) / W2;
// 	wsReal R_denom = ((wsReal)N-R_tm1)*((wsReal)N-R_tm2);
// 	if (Np_inc) *Np_inc = N;
// 
// 	DEALLOC(Ba_t1);
// 	DEALLOC(Ba_t2);
// 
// 	return R_denom<=W0 ? WISARD_NAN :
// 		((wsReal)N_con - (wsReal)N_dis) / sqrt(((wsReal)N-R_tm1)*((wsReal)N-R_tm2));
// }

wsReal _corTauIncomplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	wsUint N_dis = 0;
	wsUint N_con = 0;
	wsUint N = 0;
	char *Ba_t1 = NULL;
	char *Ba_t2 = NULL;
	wsCalloc(Ba_t1, char, N_SNP);
	wsCalloc(Ba_t2, char, N_SNP);

	if (Ra_cm) for (wsUint i=0 ; i<N_SNP ; i++) {
		if (isMissingReal(Ra_v1[i]) || isMissingReal(Ra_v2[i]) ||
			Ra_cm[i] == W0) continue;

		for (wsUint j=i+1 ; j<N_SNP ; j++) {
			if (isMissingReal(Ra_v1[j]) || isMissingReal(Ra_v2[j]) ||
				Ra_cm[i] == W0) continue;

			bool B_eq1 = Ra_v1[i] == Ra_v1[j];
			bool B_eq2 = Ra_v2[i] == Ra_v2[j];
			if (B_eq1 || B_eq2) {
				if (B_eq1) Ba_t1[i] = Ba_t1[j] = 1;
				if (B_eq2) Ba_t2[i] = Ba_t2[j] = 1;
				continue;
			}

			bool B_gt1 = Ra_v1[i] > Ra_v1[j];
			bool B_gt2 = Ra_v2[i] > Ra_v2[j];

			(B_gt1^B_gt2) ? N_dis++ : N_con++;
			N++;
		}
	} else for (wsUint i=0 ; i<N_SNP ; i++) {
		if (isMissingReal(Ra_v1[i]) || isMissingReal(Ra_v2[i])) continue;

		for (wsUint j=i+1 ; j<N_SNP ; j++) {
			if (isMissingReal(Ra_v1[j]) || isMissingReal(Ra_v2[j])) continue;

			bool B_eq1 = Ra_v1[i] == Ra_v1[j];
			bool B_eq2 = Ra_v2[i] == Ra_v2[j];
			if (B_eq1 || B_eq2) {
				if (B_eq1) Ba_t1[i] = Ba_t1[j] = 1;
				if (B_eq2) Ba_t2[i] = Ba_t2[j] = 1;
				continue;
			}

			bool B_gt1 = Ra_v1[i] > Ra_v1[j];
			bool B_gt2 = Ra_v2[i] > Ra_v2[j];

			(B_gt1^B_gt2) ? N_dis++ : N_con++;
			N++;
		}
	}
	wsUint N_t1[3] = { 0, };
	wsUint N_t2[3] = { 0, };
	for (wsUint i=0 ; i<N_SNP ; i++) {
		if (Ba_t1[i]) N_t1[(int)Ra_v1[i]]++;
		if (Ba_t2[i]) N_t2[(int)Ra_v2[i]]++;
	}
	/* Calculate tie term */
	wsReal R_tm1 = (wsReal)(N_t1[0]*(N_t1[0]-1) + N_t1[1]*(N_t1[1]-1) +
		N_t1[2]*(N_t1[2]-1)) / W2;
	wsReal R_tm2 = (wsReal)(N_t2[0]*(N_t2[0]-1) + N_t2[1]*(N_t2[1]-1) +
		N_t2[2]*(N_t2[2]-1)) / W2;
	wsReal R_denom = ((wsReal)N-R_tm1)*((wsReal)N-R_tm2);
	if (Np_inc) *Np_inc = N;

	DEALLOC(Ba_t1);
	DEALLOC(Ba_t2);

	return R_denom<=W0 ? WISARD_NAN :
		((wsReal)N_con - (wsReal)N_dis) / sqrt(((wsReal)N-R_tm1)*((wsReal)N-R_tm2));
}

inline void _ssss(const char *Na_v1, const char *Na_v2, wsUint i, wsUint j,
	char *Ba_t1, char *Ba_t2, wsUint &N_dis, wsUint &N_con)
{
	bool B_eq1 = Na_v1[i] == Na_v1[j];
	bool B_eq2 = Na_v2[i] == Na_v2[j];
	if (B_eq1 || B_eq2) {
		if (B_eq1) Ba_t1[i] = Ba_t1[j] = 1;
		if (B_eq2) Ba_t2[i] = Ba_t2[j] = 1;
	} else {
		bool B_gt1 = Na_v1[i] > Na_v1[j];
		bool B_gt2 = Na_v2[i] > Na_v2[j];

		(B_gt1^B_gt2) ? N_dis++ : N_con++;
	}
}

wsReal _corTauComplete(char *Na_v1, char *Na_v2, wsUint N_SNP,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, char *Na_cm/*=NULL*/)
{
	wsUint	N_dis		= 0;
	wsUint	N_con		= 0;
	wsUint	N_avail		= 0;

	/* Recording tie status */
	char *Ba_t1 = NULL;
	char *Ba_t2 = NULL;

	/* Rebuild Na_t1 & Na_t2 to have nonmissing */
	wsUint	N_med	= 0;
	unsigned char*	Ba_mis	= NULL;
	MULTI_ALIGNED_MALLOC(Ba_mis, unsigned char, N_SNP);
	if (OK(SSSE3)) {
		/* Check missingness */
		N_med	= get16(N_SNP);

		__m128i mi = _mm_set1_epi8(-1);
		__m128i x1 = _mm_set1_epi8(1);
		__m128i n1 = sse128zero;

		if (Na_cm) for (wsUint i=0,j=0 ; i<N_med ; i+=16) {
			__m128i *v1 = (__m128i *)(Na_v1 + i);
			__m128i *v2 = (__m128i *)(Na_v2 + i);
			__m128i *cc	= (__m128i *)(Na_cm + i);
			__m128i *cm = (__m128i *)(Ba_mis + i);

			/* ff if both available */
			*cm = sse128and(sse128and(_mm_cmpgt_epi8(*v1, mi),
				_mm_cmpgt_epi8(*v2, mi)), *cc);
			n1 = _mm_add_epi8(n1, sse128and(*cm, x1));

			/* Sum up for every 128 step to avoid overflow */
			if (++j == 128) {
				sse16CHARasum(n1, N_avail);
				n1 = sse128zero;
				j = 0;
			}
		} else for (wsUint i=0,j=0 ; i<N_med ; i+=16) {
			__m128i *v1 = (__m128i *)(Na_v1 + i);
			__m128i *v2 = (__m128i *)(Na_v2 + i);
			__m128i *cm = (__m128i *)(Ba_mis + i);

			/* ff if both available */
			*cm = sse128and(_mm_cmpgt_epi8(*v1, mi),
				_mm_cmpgt_epi8(*v2, mi));
			n1 = _mm_add_epi8(n1, sse128and(*cm, x1));

			/* Sum up for every 128 step to avoid overflow */
			if (++j == 128) {
				sse16CHARasum(n1, N_avail);
				n1 = sse128zero;
				j = 0;
			}
		}
		sse16CHARasum(n1, N_avail);
	}
	/* For non-SSE */
	if (Na_cm) for (wsUint i=N_med ; i<N_SNP ; i++) {
		if (!isMissing(Na_v1[i]) && !isMissing(Na_v2[i]) && Na_cm[i]) {
			Ba_mis[i] = 0xff;
			N_avail++;
		} else Ba_mis[i] = 0x00;
	} else for (wsUint i=N_med ; i<N_SNP ; i++)
		if (!isMissing(Na_v1[i]) && !isMissing(Na_v2[i])) {
			Ba_mis[i] = 0xff;
			N_avail++;
		} else Ba_mis[i] = 0x00;

	/* Build new vector */
	if (N_avail != N_SNP) {
		char*	Na_nv1	= NULL;
		char*	Na_nv2	= NULL;
		MULTI_ALIGNED_MALLOC(Na_nv1, char, N_avail);
		MULTI_ALIGNED_MALLOC(Na_nv2, char, N_avail);

		for (wsUint i=0,j=0 ; i<N_SNP ; i++)
			if (Ba_mis[i]) {
				Na_nv1[j] = Na_v1[i];
				Na_nv2[j] = Na_v2[i];
				j++;
			}
		Na_v1 = Na_nv1;
		Na_v2 = Na_nv2;
	}
	MULTI_ALIGNED_CALLOC(Ba_t1, char, N_avail);
	MULTI_ALIGNED_CALLOC(Ba_t2, char, N_avail);

	wsUint b2 = get16(N_avail);
	for (wsUint i=0 ; i<N_avail ; i++) {
		if (OK(SSSE3)) {
			__m128i v1 = _mm_set1_epi8(Na_v1[i]);
			__m128i v2 = _mm_set1_epi8(Na_v2[i]);

			/* Find b1 such that i+1 <= b1 < N_avail and multiply of 16 */
			wsUint b1 = ((i+16)>>4)<<4;

			__m128i di = sse128zero;
			__m128i co = sse128zero;
			__m128i t1 = sse128zero;
			__m128i t2 = sse128zero;
			__m128i x1 = _mm_set1_epi8(1);
			__m128i m1 = _mm_set1_epi8(-1);
			/* Only possible when b1 < b2 because SSE can be b1 ~ b2 */
			if (b1 < b2) {
				for (wsUint j=i+1 ; j<b1 ; j++)
					_ssss(Na_v1, Na_v2, i, j, Ba_t1, Ba_t2, N_dis, N_con);
				for (wsUint j=b1,k=0 ; j<b2 ; j+=16,k++) {
					__m128i *w1 = (__m128i *)(Na_v1 + j);
					__m128i *w2 = (__m128i *)(Na_v2 + j);
					__m128i *b1 = (__m128i *)(Ba_t1 + j);
					__m128i *b2 = (__m128i *)(Ba_t2 + j);

					__m128i e1 = _mm_cmpeq_epi8(v1, *w1);
					__m128i e2 = _mm_cmpeq_epi8(v2, *w2);
					*b1 = sse128or(*b1, sse128and(x1, e1));
					*b2 = sse128or(*b2, sse128and(x1, e2));
					__m128i ce = sse128or(e1, e2);
					t1 = sse128or(t1, e1);
					t2 = sse128or(t2, e2);
					ce = sse128xor(ce, m1);
					__m128i g1 = _mm_cmpgt_epi8(v1, *w1);
					__m128i g2 = _mm_cmpgt_epi8(v2, *w2);

					/* (g1^g2) & 0x01 because res of xor is ff / 00 */
					__m128i cdi = sse128and(x1,
						sse128and(ce, sse128xor(g1, g2)));
					g2 = sse128xor(g2, m1);
					__m128i cco = sse128and(x1,
						sse128and(ce, sse128xor(g1, g2)));

					/* Sum up disconcordance and concordance */
					di = _mm_add_epi8(di, cdi);
					co = _mm_add_epi8(co, cco);

					if (++k ==128) {
						sse16CHARasum(di, N_dis);
						sse16CHARasum(co, N_con);
						di = sse128zero;
						co = sse128zero;
						k = 0;
					}
				}
				sse16CHARasum(di, N_dis);
				sse16CHARasum(co, N_con);
				for (wsUint j=b2 ; j<N_avail ; j++)
					_ssss(Na_v1, Na_v2, i, j, Ba_t1, Ba_t2, N_dis, N_con);
					
				/* check tc == 1 then set Ba_t1[i] = 1 */
				wsUint N_d = 0;
				sse16CHARasum(t1, N_d);
				if (N_d) Ba_t1[i] = 1;
				sse16CHARasum(t2, N_d);
				if (N_d) Ba_t2[i] = 1;
			} else for (wsUint j=i+1 ; j<N_avail ; j++)
				_ssss(Na_v1, Na_v2, i, j, Ba_t1, Ba_t2, N_dis, N_con);
			/* END OF b1<b2 */
		} else for (wsUint j=i+1 ; j<N_avail ; j++)
			_ssss(Na_v1, Na_v2, i, j, Ba_t1, Ba_t2, N_dis, N_con);
	} /* END OF for(i) */
	wsReal N = ((wsReal)N_avail*(wsReal)(N_avail-1)) / W2;

	/* Count up ties */
	wsUint N_t1[3] = { 0, };
	wsUint N_t2[3] = { 0, };
	for (wsUint i=0 ; i<N_avail ; i++) {
		if (Ba_t1[i]) N_t1[(int)Na_v1[i]]++;
		if (Ba_t2[i]) N_t2[(int)Na_v2[i]]++;
	}

	/* Calculate tie term */
	wsReal R_tm1 = (wsReal)(N_t1[0]*(N_t1[0]-1) + N_t1[1]*(N_t1[1]-1) +
		N_t1[2]*(N_t1[2]-1)) / W2;
	wsReal R_tm2 = (wsReal)(N_t2[0]*(N_t2[0]-1) + N_t2[1]*(N_t2[1]-1) +
		N_t2[2]*(N_t2[2]-1)) / W2;
	wsReal R_denom = (N-R_tm1)*(N-R_tm2);
	if (Np_inc) *Np_inc = (wsUint)N;

	DEALLOC_ALIGNED(Ba_t1);
	DEALLOC_ALIGNED(Ba_t2);

	return R_denom<=W0 ? WISARD_NAN :
		((wsReal)N_con - (wsReal)N_dis) / sqrt((N-R_tm1)*(N-R_tm2));
}

inline void _ssss(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint i, wsUint j,
	char *Ba_t1, char *Ba_t2, wsUint &N_dis, wsUint &N_con)
{
	bool B_eq1 = Ra_v1[i] == Ra_v1[j];
	bool B_eq2 = Ra_v2[i] == Ra_v2[j];
	if (B_eq1 || B_eq2) {
		if (B_eq1) Ba_t1[i] = Ba_t1[j] = 1;
		if (B_eq2) Ba_t2[i] = Ba_t2[j] = 1;
	} else {
		bool B_gt1 = Ra_v1[i] > Ra_v1[j];
		bool B_gt2 = Ra_v2[i] > Ra_v2[j];

		(B_gt1^B_gt2) ? N_dis++ : N_con++;
	}
}

wsReal _corTauComplete(wsFloat *Ra_v1, wsFloat *Ra_v2, wsUint N_SNP,
	wsUint *Np_inc/*=NULL*/, wsReal *Ra_cors/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	wsUint N_dis = 0;
	wsUint N_con = 0;
	wsUint N = 0;
	char *Ba_t1 = NULL;
	char *Ba_t2 = NULL;
	wsCalloc(Ba_t1, char, N_SNP);
	wsCalloc(Ba_t2, char, N_SNP);

	if (Ra_cm) {
		for (wsUint i=0 ; i<N_SNP ; i++) {
			if (Ra_cm[i] == W0) continue;

			for (wsUint j=i+1 ; j<N_SNP ; j++) {
				if (Ra_cm[i] == W0) continue;

				bool B_eq1 = Ra_v1[i] == Ra_v1[j];
				bool B_eq2 = Ra_v2[i] == Ra_v2[j];
				if (B_eq1 || B_eq2) {
					if (B_eq1) Ba_t1[i] = Ba_t1[j] = 1;
					if (B_eq2) Ba_t2[i] = Ba_t2[j] = 1;
					continue;
				}

				bool B_gt1 = Ra_v1[i] > Ra_v1[j];
				bool B_gt2 = Ra_v2[i] > Ra_v2[j];

				(B_gt1^B_gt2) ? N_dis++ : N_con++;
				N++;
			}
		}
	} else {
		//wsUint b2 = get16(N_SNP);

		for (wsUint i=0 ; i<N_SNP ; i++)
			for (wsUint j=i+1 ; j<N_SNP ; j++)
				_ssss(Ra_v1, Ra_v2, i, j, Ba_t1, Ba_t2, N_dis, N_con);

		N = (N_SNP*(N_SNP-1))>>1;
	}
	wsUint N_t1[3] = { 0, };
	wsUint N_t2[3] = { 0, };
	for (wsUint i=0 ; i<N_SNP ; i++) {
		if (Ba_t1[i]) N_t1[(int)Ra_v1[i]]++;
		if (Ba_t2[i]) N_t2[(int)Ra_v2[i]]++;
	}
	/* Calculate tie term */
	wsReal R_tm1 = (wsReal)(N_t1[0]*(N_t1[0]-1) + N_t1[1]*(N_t1[1]-1) +
		N_t1[2]*(N_t1[2]-1)) / W2;
	wsReal R_tm2 = (wsReal)(N_t2[0]*(N_t2[0]-1) + N_t2[1]*(N_t2[1]-1) +
		N_t2[2]*(N_t2[2]-1)) / W2;
	wsReal R_denom = ((wsReal)N-R_tm1)*((wsReal)N-R_tm2);
	if (Np_inc) *Np_inc = N;

	DEALLOC(Ba_t1);
	DEALLOC(Ba_t2);

	return R_denom<=W0 ? WISARD_NAN :
		((wsReal)N_con - (wsReal)N_dis) / sqrt(((wsReal)N-R_tm1)*((wsReal)N-R_tm2));
}

wsReal _ibsComplete(char *Na_v1, char *Na_v2, wsUint N_SNP, vVariant &Xa_snp,
	wsUint *Np_inc/*=NULL*/, char *Na_cm/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	int N[3] = { 0, 0, 0 };
	static int Na_IBS[3][3] = {
		{ 2, 1, 0 },
		{ 1, 2, 1 },
		{ 0, 1, 2 }
	};

#ifdef _M_ARM
#	if 1
	if (Ra_cm) for (wsUint k=0 ; k<N_SNP ; k++) {
		if (Ra_cm[k] == W0) continue;

		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		N[Na_IBS[(int)Na_v1[k]][(int)Na_v2[k]]]++;
	} else for (wsUint k=0 ; k<N_SNP ; k++) {
		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		N[Na_IBS[(int)Na_v1[k]][(int)Na_v2[k]]]++;
	}
#	else
	if (Ra_cm) for (wsUint k=0 ; k<N_SNP ; k++) {
		if (Ra_cm[k] == W0) continue;

		if (isSexChromosome(Xa_snp[k]))
			continue;

		int X = ~(( (Na_v1[k]^Na_v2[k]) + 3) | 0x02) & 0x07;
		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		N[(X | (X >> 1))&0x03]++;
	} else for (wsUint k=0 ; k<N_SNP ; k++) {
		if (isSexChromosome(Xa_snp[k]))
			continue;

		int X = ~(( (Na_v1[k]^Na_v2[k]) + 3) | 0x02) & 0x07;
		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		N[(X | (X >> 1))&0x03]++;
	}
#	endif
#else
	/*
	(sse2) _mm_cmpeq_epi8
	(sse2) _mm_set1_epi8
	(sse2) _mm_srai_epi16
	(
	*/
	wsUint N_med = get16(N_SNP);
	__m128i n0 = _mm_set1_epi8(0);
	__m128i n1 = _mm_set1_epi8(0);
	__m128i n2 = _mm_set1_epi8(0);
	for (wsUint K=0 ; K<N_med ; K+=16) {
		__m128i *v1 = (__m128i *)(Na_v1 + K);
		__m128i *v2 = (__m128i *)(Na_v2 + K);
		__m128i *nm = (__m128i *)(Na_cm + K);
		__m128i vC = _mm_set1_epi8(3);
		__m128i vI = _mm_add_epi8(sse128xor(*v1, *v2), vC);
		vI = sse128or(vI, _mm_set1_epi8(2));
		vI = sse128and(sse128xor(_mm_cmpeq_epi8(vI, vI), vI), _mm_set1_epi8(7));
		
		vI = sse128and(sse128or(vI, _mm_srai_epi16(vI, 1)), vC);
		vC = sse128and(_mm_set1_epi8(1), *nm);
		n0 = _mm_add_epi8(sse128and(_mm_cmpeq_epi8(vI, _mm_set1_epi8(0)), vC), n0);
		n1 = _mm_add_epi8(sse128and(_mm_cmpeq_epi8(vI, vC), vC), n1);
		n2 = _mm_add_epi8(sse128and(_mm_cmpeq_epi8(vI, _mm_set1_epi8(2)), vC), n2);

		if ((K%(16<<6)) == 0) {
			sse16CHARasum(n0, N[0]);
			sse16CHARasum(n1, N[1]);
			sse16CHARasum(n2, N[2]);
			n0 = _mm_set1_epi8(0);
			n1 = _mm_set1_epi8(0);
			n2 = _mm_set1_epi8(0);
		}
	}
	sse16CHARasum(n0, N[0]);
	sse16CHARasum(n1, N[1]);
	sse16CHARasum(n2, N[2]);
	for (wsUint k=N_med ; k<N_SNP ; k++) {
		if (Ra_cm[k] == W0) continue;

		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		N[Na_IBS[(int)Na_v1[k]][(int)Na_v2[k]]]++;
	}
#endif

#if 1
	wsReal	R_ret	= W0;
	wsUint	Nsum	= N[0] + N[1] + N[2];
	if (Nsum)
		R_ret	= ((wsReal)N[1]*REAL_CONST(0.5) + (wsReal)N[2])
			/(wsReal)Nsum;
#else
	wsUint Nsum = N[0] + N[1];
	Ra_corr[i][j]	= sqrt(
		((wsReal)N[1]*REAL_CONST(0.5) + (wsReal)N[2]*W2)
		/
		(wsReal)(Nsum + N[2]*W2)
		);
#endif
	if (Np_inc) *Np_inc = Nsum;
	return R_ret;
}

wsReal _hammingComplete(char *Na_v1, char *Na_v2, wsUint N_SNP, vVariant &Xa_snp,
	wsUint *Np_inc/*=NULL*/, char *Na_cm/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	int N[2] ={ 0, 0 };

#ifdef _M_ARM
	if (Ra_cm) for (wsUint k=0 ; k<N_SNP ; k++) {
		if (Ra_cm[k] == W0) continue;

		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		N[Na_v1[k] == Na_v2[k]]++;
	} else for (wsUint k=0 ; k<N_SNP ; k++) {
		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		N[Na_v1[k] == Na_v2[k]]++;
	}
#else
	/*
	(sse2) _mm_cmpeq_epi8
	(sse2) _mm_set1_epi8
	(sse2) _mm_srai_epi16
	(
	*/
	wsUint N_med = get16(N_SNP);
	__m128i n0 = _mm_set1_epi8(0);
	__m128i n1 = _mm_set1_epi8(0);
	for (wsUint K=0 ; K<N_med ; K+=16) {
		__m128i *v1 = (__m128i *)(Na_v1 + K);
		__m128i *v2 = (__m128i *)(Na_v2 + K);
		__m128i *nm = (__m128i *)(Na_cm + K);
		__m128i vC = sse128and(_mm_set1_epi8(1), *nm);
		n0 = _mm_add_epi8(sse128and(_mm_cmpeq_epi8(*v1, *v2), vC), n0);
		n1 = _mm_add_epi8(sse128and(sse128xor(*v1, *v2), vC), n0);

		if ((K%256) == 0) {
			sse16CHARasum(n0, N[0]);
			sse16CHARasum(n1, N[1]);
			n0 = _mm_set1_epi8(0);
			n1 = _mm_set1_epi8(0);
		}
	}
	sse16CHARasum(n0, N[0]);
	sse16CHARasum(n1, N[1]);
	for (wsUint k=N_med ; k<N_SNP ; k++) {
		if (Ra_cm[k] == W0) continue;

		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		N[Na_v1[k] == Na_v2[k]]++;
	}
#endif

	wsReal	R_ret	= W0;
	wsUint	Nsum	= N[0] + N[1];
	if (Nsum)
		R_ret	= (wsReal)N[0]/(wsReal)Nsum;

	if (Np_inc) *Np_inc = Nsum;
	return R_ret;
}

wsReal _ibsIncomplete(char *Na_v1, char *Na_v2, wsUint N_SNP, vVariant &Xa_snp,
	wsUint *Np_inc/*=NULL*/, char *Na_cm/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	static int Na_IBS[3][3] = {
		{ 2, 1, 0 },
		{ 1, 2, 1 },
		{ 0, 1, 2 }
	};
	int N[3] = { 0, 0, 0 };

#ifdef _M_ARM
	if (Ra_cm) for (wsUint k=0 ; k<N_SNP ; k++) {
		if (Ra_cm[k] == W0) continue;

		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		if (isAvailable(Na_v1[k]) && isAvailable(Na_v2[k]))
			N[Na_IBS[(int)Na_v1[k]][(int)Na_v2[k]]]++;
	} else for (wsUint k=0 ; k<N_SNP ; k++) {
		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		if (isAvailable(Na_v1[k]) && isAvailable(Na_v2[k]))
			N[Na_IBS[(int)Na_v1[k]][(int)Na_v2[k]]]++;
	}
#else
	/*
	(sse2) _mm_cmpeq_epi8
	(sse2) _mm_set1_epi8
	(sse2) _mm_srai_epi16
	(
	*/
	wsUint N_med = get16(N_SNP);
	__m128i n0 = _mm_set1_epi8(0);
	__m128i n1 = _mm_set1_epi8(0);
	__m128i n2 = _mm_set1_epi8(0);
	for (wsUint K=0 ; K<N_med ; K+=16) {
		__m128i *v1 = (__m128i *)(Na_v1 + K);
		__m128i *v2 = (__m128i *)(Na_v2 + K);
		__m128i *nm = (__m128i *)(Na_cm + K);
		__m128i vC = _mm_set1_epi8(-9);
		__m128i vX = sse128or(_mm_cmpgt_epi8(*v1, vC), _mm_cmpgt_epi8(*v2, vC));
		vC = _mm_set1_epi8(3);
		__m128i vI = _mm_add_epi8(sse128xor(*v1, *v2), vC);
		vI = sse128or(vI, _mm_set1_epi8(2));
		vI = sse128and(sse128xor(_mm_cmpeq_epi8(vI, vI), vI), _mm_set1_epi8(7));

		vI = sse128and(sse128or(vI, _mm_srai_epi16(vI, 1)), vC);
		vC = sse128and(_mm_set1_epi8(1), sse128and(vX, *nm));
		n0 = _mm_add_epi8(sse128and(_mm_cmpeq_epi8(vI, _mm_set1_epi8(0)), vC), n0);
		n1 = _mm_add_epi8(sse128and(_mm_cmpeq_epi8(vI, vC), vC), n1);
		n2 = _mm_add_epi8(sse128and(_mm_cmpeq_epi8(vI, _mm_set1_epi8(2)), vC), n2);

		if ((K%(16<<6)) == 0) {
			sse16CHARasum(n0, N[0]);
			sse16CHARasum(n1, N[1]);
			sse16CHARasum(n2, N[2]);
			n0 = _mm_set1_epi8(0);
			n1 = _mm_set1_epi8(0);
			n2 = _mm_set1_epi8(0);
		}
	}
	sse16CHARasum(n0, N[0]);
	sse16CHARasum(n1, N[1]);
	sse16CHARasum(n2, N[2]);
	for (wsUint k=N_med ; k<N_SNP ; k++) {
		if (Ra_cm[k] == W0) continue;

		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸� ��� */
		if (isAvailable(Na_v1[k]) && isAvailable(Na_v2[k]))
			N[Na_IBS[(int)Na_v1[k]][(int)Na_v2[k]]]++;
	}
#endif

#if 1
	wsReal	R_ret	= W0;
	wsUint	Nsum	= N[0] + N[1] + N[2];
	if (Nsum)
		R_ret	= ((wsReal)N[1]*REAL_CONST(0.5) + (wsReal)N[2])
		/(wsReal)Nsum;
#else
	wsUint Nsum = N[0] + N[1];
	Ra_corr[i][j]	= sqrt(
		((wsReal)N[1]*REAL_CONST(0.5) + (wsReal)N[2]*W2)
		/
		(wsReal)(Nsum + N[2]*W2)
		);
#endif
	if (Np_inc) *Np_inc = Nsum;
	return R_ret;
}
#if 0
int cCorrAnalysis::calcCorrV2(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	xCorrThread
			*Xp_ct		= (xCorrThread *)Vp_shareData;
	wsFloat	**Ra_data	= Xp_ct->Ra_data;
	wsReal	**Ra_corr	= Xp_ct->Ra_corr;
	char	**Na_data	= Xp_ct->Na_data;
	int		**Na_corrN	= Xp_ct->Na_corrN;
	wsUint	N_sample	= Xp_ct->N_sample;
	wsUint	N_SNP		= Xp_ct->N_variant;
	wsFloatCst	*Ra_cm		= Xp_ct->Ra_corrMask;
	char	*Na_cm		= Xp_ct->Na_corrMask;
	int		i			= *((int *)Vp_data);
	wsReal	*Ra_cors	= NULL;
	wsUint	N_inc;

	/* Allocate buffer for calculating correlation by --medcor */
	if (OPT_ENABLED(medcor))
		sseMalloc(Ra_cors, wsReal, N_SNP);

 	if (0) {//Xp_ct->B_isComplete) {
		if (OPT_ENABLED(medcor)) {
			/* Complete genotype, median */
			for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corMedianComplete(Ra_data[i], Ra_data[j],
					N_SNP, &N_inc, Ra_cors, Ra_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
		} else if (OPT_ENABLED(ktau)) {
			halt("Can't be ktau in here!");
		} else {
			/* Complete genotype, mean */
			for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corMeanComplete(Ra_data[i], Ra_data[j],
					N_SNP, &N_inc, Ra_cors, Ra_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
		}
 	} else {
		if (OPT_ENABLED(medcor)) {
			/* Incomplete genotype, median */
			for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corMedianIncomplete(Ra_data[i], Ra_data[j],
					N_SNP, &N_inc, Ra_cors, Ra_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
		} else if (OPT_ENABLED(ktau)) {
			/* Incomplete genotype, Kendall's tau */
			if (IS_ASSIGNED(dosage)) for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corTauComplete(Ra_data[i], Ra_data[j],
					N_SNP, &N_inc, Ra_cors, Ra_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			} else for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corTauComplete(Na_data[i], Na_data[j],
					N_SNP, &N_inc, Ra_cors, Na_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
		} else {
			/* Incomplete genotype, mean */
			for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corMeanIncomplete(Ra_data[i], Ra_data[j],
					N_SNP, &N_inc, Ra_cors, Ra_cm);

				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
		}
	}
	if (OPT_ENABLED(medcor))
		sseFree(Ra_cors);

	return 0;
}
#endif // endif
wsReal _hammingIncomplete(char *Na_v1, char *Na_v2, wsUint N_SNP, vVariant &Xa_snp,
	wsUint *Np_inc/*=NULL*/, char *Na_cm/*=NULL*/, wsFloatCst *Ra_cm/*=NULL*/)
{
	int N[2] ={ 0, 0 };

#ifdef _M_ARM
	if (Ra_cm) for (wsUint k=0 ; k<N_SNP ; k++) {
		if (Ra_cm[k] == W0) continue;

		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸�?���?*/
		if (isAvailable(Na_v1[k]) && isAvailable(Na_v2[k]))
			N[Na_v1[k] == Na_v2[k]]++;
	}
	else for (wsUint k=0 ; k<N_SNP ; k++) {
		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸�?���?*/
		if (isAvailable(Na_v1[k]) && isAvailable(Na_v2[k]))
			N[Na_v1[k] == Na_v2[k]]++;
	}
#else
	/*
	(sse2) _mm_cmpeq_epi8
	(sse2) _mm_set1_epi8
	(sse2) _mm_srai_epi16
	(
	*/
	wsUint N_med = get16(N_SNP);
	__m128i n0 = _mm_set1_epi8(0);
	__m128i n1 = _mm_set1_epi8(0);
	for (wsUint K=0 ; K<N_med ; K+=16) {
		__m128i *v1 = (__m128i *)(Na_v1 + K);
		__m128i *v2 = (__m128i *)(Na_v2 + K);
		__m128i *nm = (__m128i *)(Na_cm + K);
		__m128i vC = _mm_set1_epi8(-9);
		__m128i vX = sse128or(_mm_cmpgt_epi8(*v1, vC), _mm_cmpgt_epi8(*v2, vC));
		vC = sse128and(_mm_set1_epi8(1), sse128and(vX, *nm));
		n0 = _mm_add_epi8(sse128and(_mm_cmpeq_epi8(*v1, *v2), vC), n0);
		n1 = _mm_add_epi8(sse128and(sse128xor(*v1, *v2), vC), n1);

		if ((K%(16<<6)) == 0) {
			sse16CHARasum(n0, N[0]);
			sse16CHARasum(n1, N[1]);
			n0 = _mm_set1_epi8(0);
			n1 = _mm_set1_epi8(0);
		}
	}
	sse16CHARasum(n0, N[0]);
	sse16CHARasum(n1, N[1]);
	for (wsUint k=N_med ; k<N_SNP ; k++) {
		if (Ra_cm[k] == W0) continue;

		if (isSexChromosome(Xa_snp[k]))
			continue;

		/* use="pair", i�� j ���� �� �� k��° SNP�� ���� genotype�� �־�߸�?���?*/
		if (isAvailable(Na_v1[k]) && isAvailable(Na_v2[k]))
			N[Na_v1[k] == Na_v2[k]]++;
	}
#endif

	wsReal	R_ret	= W0;
	wsUint	Nsum	= N[0] + N[1];
	if (Nsum)
		R_ret	= (wsReal)N[0]/(wsReal)Nsum;

	if (Np_inc) *Np_inc = Nsum;
	return R_ret;
}

int cCorrAnalysis::calcCorrV2(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	xCorrThread
		*Xp_ct		= (xCorrThread *)Vp_shareData;
	wsFloat	**Ra_data	= Xp_ct->Ra_data;
	wsReal	**Ra_corr	= Xp_ct->Ra_corr;
	char	**Na_data	= Xp_ct->Na_data;
	int		**Na_corrN	= Xp_ct->Na_corrN;
	wsUint	N_sample	= Xp_ct->N_sample;
	wsUint	N_SNP		= Xp_ct->N_variant;
	wsFloatCst	*Ra_cm		= Xp_ct->Ra_corrMask;
	char	*Na_cm		= Xp_ct->Na_corrMask;
	int		i			= *((int *)Vp_data);
	wsReal	*Ra_cors	= NULL;
	wsUint	N_inc;

	/* Allocate buffer for calculating correlation by --medcor */
	if (OPT_ENABLED(medcor)) {
		sseMalloc(Ra_cors, wsReal, N_SNP);
	}

	if (0) {//Xp_ct->B_isComplete) {
		if (OPT_ENABLED(medcor)) {
			/* Complete genotype, median */
			for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corMedianComplete(Ra_data[i], Ra_data[j],
					N_SNP, &N_inc, Ra_cors, Ra_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
		}
		else if (OPT_ENABLED(ktau)) {
			halt("Can't be ktau in here!");
		}
		else {
			/* Complete genotype, mean */
			for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corMeanComplete(Ra_data[i], Ra_data[j],
					N_SNP, &N_inc, Ra_cors, Ra_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
		}
	}
	else {
		if (OPT_ENABLED(medcor)) {
			/* Incomplete genotype, median */
			for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corMedianIncomplete(Ra_data[i], Ra_data[j],
					N_SNP, &N_inc, Ra_cors, Ra_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
		}
		else if (OPT_ENABLED(ktau)) {
			/* Incomplete genotype, Kendall's tau */
			if (IS_ASSIGNED(dosage)) for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corTauComplete(Ra_data[i], Ra_data[j],
					N_SNP, &N_inc, Ra_cors, Ra_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
			else for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corTauComplete(Na_data[i], Na_data[j],
					N_SNP, &N_inc, Ra_cors, Na_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
		}
		else {
			/* Incomplete genotype, mean */
			for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corMeanIncomplete(Ra_data[i], Ra_data[j],
					N_SNP, &N_inc, Ra_cors, Ra_cm);

				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				Na_corrN[i][j]	= N_inc;
				Na_corrN[j][i]	= N_inc;
			}
		}
	}
	if (OPT_ENABLED(medcor))
		sseFree(Ra_cors);

	return 0;
}

int cCorrAnalysis::calcCorrV2pearson(int N_idx, void *Vp_shareData,
	void *Vp_data, void *Vp_result)
{
	xCorrThread
		*Xp_ct		= (xCorrThread *)Vp_shareData;
	char	**Na_data	= Xp_ct->Na_data;
	wsFloat	**Ra_data	= Xp_ct->Ra_data;
	wsReal	**Ra_corr	= Xp_ct->Ra_corr;
	int		**Na_corrN	= Xp_ct->Na_corrN;
	wsUint	N_sample	= Xp_ct->N_sample;
	wsUint	N_SNP		= Xp_ct->N_variant;
	wsFloatCst	*Ra_cm		= Xp_ct->Ra_corrMask;
	char	*Na_cm		= Xp_ct->Na_corrMask;
	int		i			= *((int *)Vp_data);
	wsReal	*Ra_cors	= NULL;
	wsUint	N_inc;

	if (Na_data == NULL) {
		if (!IS_ASSIGNED(dosage))
			halt("SYSERR : Genotype data is NULL w/o --dosage");

		if (Xp_ct->B_isComplete) {
			if (OPT_ENABLED(medcor)) {
				halt_fmt(WISARD_CANT_MEDCOR_W_CORPEARSON);
				//			halt("--medcor cannot be used with --corpearson");
			} else {
				/* Complete genotype, mean */
				for (wsUint j=i ; j<N_sample ; j++) {
					wsReal R_res = _corPearsonComplete(Ra_data[i], Ra_data[j],
						N_SNP, &N_inc, Ra_cors, Ra_cm);
					Ra_corr[i][j]	= R_res;
					Ra_corr[j][i]	= R_res;
					if (Na_corrN) {
						Na_corrN[i][j]	= N_inc;
						Na_corrN[j][i]	= N_inc;
					}
				}
			}
		} else {
			if (OPT_ENABLED(medcor)) {
				halt_fmt(WISARD_CANT_MEDCOR_W_CORPEARSON);
				//			halt("--medcor cannot be used with --corpearson");
			} else {
				/* Incomplete genotype, mean */
				for (wsUint j=i ; j<N_sample ; j++) {
					wsReal R_res = _corPearsonIncomplete(Ra_data[i], Ra_data[j],
						N_SNP, &N_inc, Ra_cors, Ra_cm);
					Ra_corr[i][j]	= R_res;
					Ra_corr[j][i]	= R_res;
					if (Na_corrN) {
						Na_corrN[i][j]	= N_inc;
						Na_corrN[j][i]	= N_inc;
					}
				}
			}
		}
	} else {
		/* Allocate buffer for calculating correlation by --medcor */
		if (OPT_ENABLED(medcor))
			sseMalloc(Ra_cors, wsReal, N_SNP);

		if (Xp_ct->B_isComplete) {
			if (OPT_ENABLED(medcor)) {
				halt_fmt(WISARD_CANT_MEDCOR_W_CORPEARSON);
	//			halt("--medcor cannot be used with --corpearson");
			} else {
				/* Complete genotype, mean */
				for (wsUint j=i ; j<N_sample ; j++) {
					wsReal R_res = _corMeanComplete(Na_data[i], Na_data[j],
						N_SNP, &N_inc, Ra_cors, Ra_cm);
					Ra_corr[i][j]	= R_res;
					Ra_corr[j][i]	= R_res;
					if (Na_corrN) {
						Na_corrN[i][j]	= N_inc;
						Na_corrN[j][i]	= N_inc;
					}
				}
			}
		} else {
			if (OPT_ENABLED(medcor)) {
				halt_fmt(WISARD_CANT_MEDCOR_W_CORPEARSON);
	//			halt("--medcor cannot be used with --corpearson");
			} else {
				/* Incomplete genotype, mean */
				for (wsUint j=i ; j<N_sample ; j++) {
// 					wsReal R_res = _corMeanIncomplete(Na_data[i], Na_data[j],
// 						N_SNP, &N_inc, Ra_cors, Ra_cm);
					wsReal R_res = _corPearsonIncomplete(Na_data[i], Na_data[j],
						N_SNP, &N_inc, Ra_cors, Na_cm);
					Ra_corr[i][j]	= R_res;
					Ra_corr[j][i]	= R_res;
					if (Na_corrN) {
						Na_corrN[i][j]	= N_inc;
						Na_corrN[j][i]	= N_inc;
					}
				}
			}
		}
		if (OPT_ENABLED(medcor))
			sseFree(Ra_cors);
	}

	return 0;
}

int cCorrAnalysis::calcCorrV2tauS1(int N_idx, void *Vp_shareData,
	void *Vp_data, void *Vp_result)
{
	xCorrThread	*Xp_ct		= (xCorrThread *)Vp_shareData;
	WISARD_pc	**Na_gt		= Xp_ct->Na_gt;
	WISARD_pc	**Na_neq	= Xp_ct->Na_neq;
	wsUint		**Na_same	= Xp_ct->Na_same;
	wsUint		N_SNP		= Xp_ct->N_variant;
	wsFloatCst		*Ra_cm		= Xp_ct->Ra_corrMask;
	int			i			= *((int *)Vp_data);
	wsUint		N_szBlk		= Xp_ct->N_szBlk;

	/* Allocate memory and set size */
	wsAlloc(Na_gt[i], WISARD_pc, N_szBlk);
	wsAlloc(Na_neq[i], WISARD_pc, N_szBlk);
	WISARD_pc	*Np_gt		= Na_gt[i];
	WISARD_pc	*Np_eq		= Na_neq[i];

	/* For the given sample, calculate comparison result */
	WISARD_pc N_gt = 0, N_eq = 0;
	wsUint X = 0, Y = 0;

	/* If --ktau */
	if (OPT_ENABLED(ktau)) {
		char		**Na_data	= Xp_ct->Na_data;	/* Assumed to be PERFECTLY genotyped */

		/* Set same-counting buffer */
		wsAlloc(Na_same[i], wsUint, 3);
		wsUint		*Np_same	= Na_same[i];
		Np_same[0] = Np_same[1] = Np_same[2] = 0;

		/* For the given sample, calculate comparison result */
		char *Np_data = Na_data[i];

		if (Ra_cm) for (wsUint j=0 ; j<N_SNP ; j++) {
			if (Ra_cm[j] == W0) continue;
			Np_same[(wsUint)Np_data[j]]++;

			char N_geno = Np_data[j];
			for (wsUint k=j+1 ; k<N_SNP ; k++) {
				if (Ra_cm[k] == W0) continue;

				N_gt |= N_geno > Np_data[k];
				N_eq |= N_geno != Np_data[k];

				if (++X == WISARD_szpc) {
					Np_gt[Y] = N_gt;
					Np_eq[Y] = N_eq;
					Y++;
					X = 0;
					N_gt = N_eq = 0;
				} else {
					N_gt <<= 1;
					N_eq <<= 1;
				}
			}
		} else  for (wsUint j=0 ; j<N_SNP ; j++) {
			Np_same[(wsUint)Np_data[j]]++;

			char N_geno = Np_data[j];
			for (wsUint k=j+1 ; k<N_SNP ; k++) {
				N_gt |= N_geno > Np_data[k];
				N_eq |= N_geno != Np_data[k];

				if (++X == WISARD_szpc) {
					Np_gt[Y] = N_gt;
					Np_eq[Y] = N_eq;
					Y++;
					X = 0;
					N_gt = N_eq = 0;
				} else {
					N_gt <<= 1;
					N_eq <<= 1;
				}
			}
		}
		if (X) {
			N_gt >>= 1;
			N_eq >>= 1;
			Np_gt[Y] = N_gt;
			Np_eq[Y] = N_eq;
			Y++;
		}
		/* Adjust same */
		if (Np_same[0]==1) Np_same[0] = 0;
		if (Np_same[1]==1) Np_same[1] = 0;
		if (Np_same[2]==1) Np_same[2] = 0;
	} else if (OPT_ENABLED(empktau)) {
		wsFloat **Ra_data	= Xp_ct->Ra_data;	/* Assumed to be PERFECTLY genotyped */

		/* Set same-counting MAP */
		realMap	Xm_same;
		wsFloat *Rp_data = Ra_data[i];

		if (Ra_cm) for (wsUint j=0 ; j<N_SNP ; j++) {
			if (Ra_cm[j] == W0) continue;
			wsReal R_geno = Rp_data[j];

			realMap_it X_find = Xm_same.find((float)R_geno);
			if (X_find == Xm_same.end())
				Xm_same[(float)R_geno] = 1;
			else
				Xm_same[(float)R_geno]++;

			for (wsUint k=j+1 ; k<N_SNP ; k++) {
				if (Ra_cm[k] == W0) continue;

				N_gt |= R_geno > Rp_data[k];
				N_eq |= R_geno != Rp_data[k];

				if (++X == WISARD_szpc) {
					Np_gt[Y] = N_gt;
					Np_eq[Y] = N_eq;
					Y++;
					X = 0;
					N_gt = N_eq = 0;
				} else {
					N_gt <<= 1;
					N_eq <<= 1;
				}
			}
		} else for (wsUint j=0 ; j<N_SNP ; j++) {
			wsReal R_geno = Rp_data[j];

			realMap_it X_find = Xm_same.find((float)R_geno);
			if (X_find == Xm_same.end())
				Xm_same[(float)R_geno] = 1;
			else
				Xm_same[(float)R_geno]++;

			for (wsUint k=j+1 ; k<N_SNP ; k++) {
				N_gt |= R_geno > Rp_data[k];
				N_eq |= R_geno != Rp_data[k];

				if (++X == WISARD_szpc) {
					Np_gt[Y] = N_gt;
					Np_eq[Y] = N_eq;
					Y++;
					X = 0;
					N_gt = N_eq = 0;
				} else {
					N_gt <<= 1;
					N_eq <<= 1;
				}
			}
		}
		if (X) {
			N_gt >>= 1;
			N_eq >>= 1;
			Np_gt[Y] = N_gt;
			Np_eq[Y] = N_eq;
			Y++;
		}
		/* Make same buffer */
		size_t N_szBuf = Xm_same.size();
		wsAlloc(Na_same[i], wsUint, N_szBuf+1);

		/* Insert and adjust same */
		wsUint *Np_same = Na_same[i];
		wsUint k=0;
		FOREACHDO (realMap_it, Xm_same, it, k++)
			Np_same[k] = it->second> 1 ? it->second : 0;
		/* End mark */
		Np_same[k] = 0xffffffff;
	}
	/* Sanity check */
	if (Y != N_szBlk) halt("COUNTING ERR");

	return 0;
}

		
int cCorrAnalysis::calcCorrV2tauS2(int N_idx, void *Vp_shareData,
	void *Vp_data, void *Vp_result)
{
	xCorrThread	*Xp_ct		= (xCorrThread *)Vp_shareData;
	wsReal		**Ra_corr	= Xp_ct->Ra_corr;
	int			**Na_corrN	= Xp_ct->Na_corrN;
	wsUint		N_sample	= Xp_ct->N_sample;
	wsUint		N_SNP		= Xp_ct->N_szSNP;
	int			i			= *((int *)Vp_data);
	wsReal		*Ra_cors	= NULL;
	WISARD_pc	*Na_g1		= Xp_ct->Na_gt[i];
	WISARD_pc	*Na_n1		= Xp_ct->Na_neq[i];
	wsUint		*Na_s1		= Xp_ct->Na_same[i];
	wsUint		N_szBlk		= Xp_ct->N_szBlk;

	/* Allocate buffer for calculating correlation by --medcor */
	if (OPT_ENABLED(medcor))
		halt_fmt(WISARD_CANT_MEDCOR_W_CORPEARSON);

	wsReal R_N = ((wsReal)N_SNP * (wsReal)(N_SNP-1)) / W2;

	if (Xp_ct->B_isComplete) {
		/* Complete genotype, mean */
		for (wsUint j=i ; j<N_sample ; j++) {
			WISARD_pc	*Na_g2		= Xp_ct->Na_gt[j];
			WISARD_pc	*Na_n2		= Xp_ct->Na_neq[j];
			wsUint		*Na_s2		= Xp_ct->Na_same[j];
			wsUint		N_sumDisc	= 0;
			wsUint		N_sumConc	= 0;

			for (wsUint k=0 ; k<N_szBlk ; k++) {
				WISARD_pc N_norSame	= Na_n1[k] & Na_n2[k];
				WISARD_pc N_disc	= N_norSame & (Na_g1[k] ^ Na_g2[k]);
				WISARD_pc N_conc	= N_norSame & (Na_g1[k] ^ ~Na_g2[k]);

				N_sumDisc += (wsUint)WISARD_POPCNT(N_disc);
				N_sumConc += (wsUint)WISARD_POPCNT(N_conc);
			}

			wsReal R_t = W0;
			wsReal R_u = W0;
			if (OPT_ENABLED(ktau)) {
				R_t = (wsReal)(Na_s1[0]*(Na_s1[0]-1) + Na_s1[1]*(Na_s1[1]-1)
					+ Na_s1[2]*(Na_s1[2]-1)) / W2;
				R_u = (wsReal)(Na_s2[0]*(Na_s2[0]-1) + Na_s2[1]*(Na_s2[1]-1)
					+ Na_s2[2]*(Na_s2[2]-1)) / W2;
			} else if (OPT_ENABLED(empktau)) {
				for (wsUint k=0 ; Na_s1[k]!=0xffffffff ; k++)
					R_t += (wsReal)Na_s1[k] * (wsReal)(Na_s1[k]-1);
				for (wsUint k=0 ; Na_s2[k]!=0xffffffff ; k++)
					R_u += (wsReal)Na_s2[k] * (wsReal)(Na_s2[k]-1);
				R_t /= W2;
				R_u /= W2;
			}
			wsReal R_denom = (R_N - R_t) * (R_N - R_u);

			wsReal R_res;
			if (R_denom > W0)
				R_res = ((wsReal)N_sumConc - (wsReal)N_sumDisc) / sqrt(R_denom);
			else
				R_res = WISARD_NAN;

			Ra_corr[i][j]	= R_res;
			Ra_corr[j][i]	= R_res;
			if (Na_corrN) {
				Na_corrN[i][j]	= (wsUint)R_N;
				Na_corrN[j][i]	= (wsUint)R_N;
			}
		}
	} else {
		halt("Still not support!");
	}
	if (OPT_ENABLED(medcor))
		sseFree(Ra_cors);

	return 0;
}

int cCorrAnalysis::calcCorrV2tau(int N_idx, void *Vp_shareData,
	void *Vp_data, void *Vp_result)
{
	xCorrThread
		*Xp_ct		= (xCorrThread *)Vp_shareData;
	char	**Na_data	= Xp_ct->Na_data;
	wsReal	**Ra_corr	= Xp_ct->Ra_corr;
	int		**Na_corrN	= Xp_ct->Na_corrN;
	wsUint	N_sample	= Xp_ct->N_sample;
	wsUint	N_SNP		= Xp_ct->N_variant;
	//REAL_c	*Ra_cm		= Xp_ct->Ra_corrMask;
	char	*Na_cm		= Xp_ct->Na_corrMask;
	int		i			= *((int *)Vp_data);
	wsReal	*Ra_cors	= NULL;
	wsUint	N_inc;

	/* Allocate buffer for calculating correlation by --medcor */
	if (OPT_ENABLED(medcor))
		sseMalloc(Ra_cors, wsReal, N_SNP);

	if (Xp_ct->B_isComplete) {
		if (OPT_ENABLED(medcor)) {
			halt_fmt(WISARD_CANT_MEDCOR_W_CORPEARSON);
			//			halt("--medcor cannot be used with --corpearson");
		} else {
			/* Complete genotype, mean */
			for (wsUint j=i ; j<N_sample ; j++) {
				wsReal R_res = _corTauComplete(Na_data[i], Na_data[j],
					N_SNP, &N_inc, Ra_cors, Na_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				if (Na_corrN) {
					Na_corrN[i][j]	= N_inc;
					Na_corrN[j][i]	= N_inc;
				}
			}
		}
	} else {
		if (OPT_ENABLED(medcor)) {
			halt_fmt(WISARD_CANT_MEDCOR_W_CORPEARSON);
			//			halt("--medcor cannot be used with --corpearson");
		} else {
			/* Incomplete genotype, mean */
			for (wsUint j=i ; j<N_sample ; j++) {
// 				wsReal R_res = _corTauIncomplete(Na_data[i], Na_data[j],
// 					N_SNP, &N_inc, Ra_cors, Ra_cm);
				wsReal R_res = _corTauComplete(Na_data[i], Na_data[j],
				N_SNP, &N_inc, Ra_cors, Na_cm);
				Ra_corr[i][j]	= R_res;
				Ra_corr[j][i]	= R_res;
				if (Na_corrN) {
					Na_corrN[i][j]	= N_inc;
					Na_corrN[j][i]	= N_inc;
				}
			}
		}
	}
	if (OPT_ENABLED(medcor))
		sseFree(Ra_cors);

	return 0;
}

#ifdef USE_CUDA

#define CORR_THREADS 32
int cCorrAnalysis::calcCorrV2_gpu(int N_idx, void *Vp_shareData,
	void *Vp_data, void *Vp_result)
{
	xCorrThread
			*Xp_ct		= (xCorrThread *)Vp_shareData;
	wsReal	**Ra_data	= Xp_ct->Ra_data;
	wsReal	**Ra_corr	= Xp_ct->Ra_corr;
	int		**Na_corrN	= Xp_ct->Na_corrN;
	wsUint	N_sample	= Xp_ct->N_sample;
	wsUint	N_SNP		= Xp_ct->N_variant;
	int		i			= *((int *)Vp_data);
	wsUint	sz			= *(((wsUint *)Vp_data)+1);
	wsReal	*Ra_cors	= NULL;
	wsUint	N_med		= getMed(N_SNP);

	if (OPT_ENABLED(medcor))
		sseMalloc(Ra_cors, wsReal, N_SNP);
	//	printf("%d,*th column\n", i);

	wsReal *g_Ra_template	= NULL;
	wsReal *g_Ra_target	= NULL;
	wsReal *g_Ra_ret		= NULL;

	if (Xp_ct->B_isComplete) {
		if (OPT_ENABLED(medcor)) {
			halt("Still not support this");
		} else {
			/* Initially, copy 'sz' samples into memory */
			cuAlloc(g_Ra_template, sizeof(wsReal)*sz*N_SNP);
			cuAlloc(g_Ra_target, sizeof(wsReal)*sz*N_SNP);
			cuAlloc(g_Ra_ret, sizeof(wsReal)*sz*sz);

			/* Copy template */
			wsUint N_offset = 0;
			for (wsUint j=0 ; j<sz ; j++) {
				cuSend(g_Ra_template+N_offset, Ra_data[i+j],
					sizeof(wsReal)*N_SNP);
				N_offset += N_SNP;
			}

			/* Launch first kernel that g_Ra_template == g_Ra_target */
			notice("Calculating %d-%d...\r", 0, sz-1);
			wisard_corrMain(g_Ra_template, g_Ra_template, N_SNP, g_Ra_ret,
				sz, sz, CORR_THREADS);

			/* Fill the square part of diagonal direction
			 * Result size is sz*sz */
			N_offset = 0;
			for (wsUint j=0 ; j<sz ; j++) {
				cuRecv(Ra_corr[i+j]+i, g_Ra_ret+N_offset,
					sizeof(wsReal)*sz);
				N_offset += sz;

				for (wsUint k=0 ; k<sz ; k++)
					Na_corrN[i+j][i+k] = N_SNP;
			}

			/* Successively launch kernels */
			for (wsUint j=i+sz ; j<N_sample ; j+=sz) {
				wsUint N_curSz = j+sz > N_sample ? N_sample-j : sz;
				if ((j+sz-1) >= N_sample)
					notice("Calculating %d-%d...\r", j, N_sample);
				else
					notice("Calculating %d-%d...\r", j, j+sz-1);

				/* Make target */
				N_offset = 0;
				for (wsUint k=0 ; k<N_curSz ; k++) {
					cuSend(g_Ra_target+N_offset, Ra_data[j],
						sizeof(wsReal)*N_SNP);
					N_offset += N_SNP;
				}

				wisard_corrMain(g_Ra_template, g_Ra_target, N_SNP, g_Ra_ret,
					sz, N_curSz, CORR_THREADS);

				/* Result size will be sz * N_curSz */
				N_offset = 0;
				for (wsUint k=0 ; k<sz ; k++) {
					/* Fill horizontal part */
					cuRecv(Ra_corr[i+k]+j, g_Ra_ret+N_offset,
						sizeof(wsReal)*N_curSz);

					/* Fill vertical part */
					for (wsUint l=j ; l<N_curSz ; l++)
						Ra_corr[l][i+k] = Ra_corr[i+k][l];

					/* One line = N_curSz */
					N_offset += N_curSz;
				}

				for (wsUint k=0 ; k<sz ; k++) {
					Na_corrN[j][i+k] = N_SNP;
					Na_corrN[i+k][j] = N_SNP;
				}
			}

			cuFree(g_Ra_template);
			cuFree(g_Ra_target);
			cuFree(g_Ra_ret);
		}
	} else {
		if (OPT_ENABLED(medcor)) {
			halt("Still not support this");
		} else {
			wsUint *g_Na_cntInc = NULL;

			/* Initially, copy 'sz' samples into memory */
			cuAlloc(g_Ra_template, sizeof(wsReal)*sz*N_SNP);
			cuAlloc(g_Ra_target, sizeof(wsReal)*sz*N_SNP);
			cuAlloc(g_Ra_ret, sizeof(wsReal)*sz*sz);
			cuAlloc(g_Na_cntInc, sizeof(wsReal)*sz*sz);

			/* Copy template */
			wsUint N_offset = 0;
			for (wsUint j=0 ; j<sz ; j++) {
				cuSend(g_Ra_template+N_offset, Ra_data[i+j],
					sizeof(wsReal)*N_SNP);
				N_offset += N_SNP;
			}

			/* Launch first kernel that g_Ra_template == g_Ra_target */
			wisard_corrMain(g_Ra_template, g_Ra_template, N_SNP, g_Ra_ret,
				sz, sz, 32, g_Na_cntInc);
			N_offset = 0;
			for (wsUint j=0 ; j<sz ; j++) {
				cuRecv(Ra_corr[i+j]+i, g_Ra_ret+N_offset,
					sizeof(wsReal)*sz);
				N_offset += sz;

				for (wsUint k=0 ; k<sz ; k++)
					Na_corrN[i+j][i+k] = N_SNP;
			}

			/* Successively launch kernels */
			for (wsUint j=i+sz ; j<N_sample ; j+=sz) {
				wsUint N_curSz = j+sz > N_sample ? N_sample-j : sz;

				/* Make target */
				N_offset = 0;
				for (wsUint k=0 ; k<N_curSz ; k++) {
					cuSend(g_Ra_target+N_offset, Ra_data[j],
						sizeof(wsReal)*N_SNP);
					N_offset += N_SNP;
				}

				wisard_corrMain(g_Ra_template, g_Ra_target, N_SNP, g_Ra_ret,
					sz, N_curSz, 32, g_Na_cntInc);
				N_offset = 0;
				for (wsUint k=0 ; k<sz ; k++) {
					cuRecv(Ra_corr[i+k]+j, g_Ra_ret+N_offset,
						sizeof(wsReal)*N_curSz);
					for (wsUint l=j ; l<N_curSz ; l++)
						Ra_corr[l][i+k] = Ra_corr[i+k][l];
//					memcpy(Ra_corr[j]+i+k, Ra_corr[i+k]+j,
//						sizeof(wsReal)*N_curSz);
					N_offset += N_curSz;
				}

				for (wsUint k=0 ; k<sz ; k++) {
					Na_corrN[j][i+k] = N_SNP;
					Na_corrN[i+k][j] = N_SNP;
				}
			}

			cuFree(g_Ra_template);
			cuFree(g_Ra_target);
			cuFree(g_Ra_ret);
		}
		//		Ra_corr[i][i]	= W1;
		//		Na_corrN[i][i]	= N_SNP;
	}
	sseFree(Ra_cors);

	return 0;
}

#endif

int cCorrAnalysis::calcIBS(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	xCorrThread
		*Xp_ct		= (xCorrThread *)Vp_shareData;
	char	**Na_data	= Xp_ct->Na_data;
	wsReal	**Ra_corr	= Xp_ct->Ra_corr;
	int		**Na_corrN	= Xp_ct->Na_corrN;
	wsUint	N_sample	= Xp_ct->N_sample;
	wsUint	N_SNP		= Xp_ct->N_variant;
	char*	Na_cm		= Xp_ct->Na_corrMask;
	wsFloatCst*	Ra_cm		= Xp_ct->Ra_corrMask;
	vVariant	&Xa_SNP		= *(Xp_ct->Xp_SNP);
	int		i			= *((int *)Vp_data);

	if (Xp_ct->B_isComplete) {
		for (wsUint j=i ; j<N_sample ; j++) {
			wsUint N_inc;
			wsReal R_res = _ibsComplete(Na_data[i], Na_data[j], N_SNP,
				Xa_SNP, &N_inc, Na_cm, Ra_cm);
			Ra_corr[i][j]	= R_res;
			Ra_corr[j][i]	= R_res;
			Na_corrN[i][j]	= N_inc;
			Na_corrN[j][i]	= N_inc;

		}
	} else {
		for (wsUint j=i ; j<N_sample ; j++) {
			wsUint N_inc;
			wsReal R_res = _ibsIncomplete(Na_data[i], Na_data[j], N_SNP,
				Xa_SNP, &N_inc, Na_cm, Ra_cm);
			Ra_corr[i][j]	= R_res;
			Ra_corr[j][i]	= R_res;
			Na_corrN[i][j]	= N_inc;
			Na_corrN[j][i]	= N_inc;

		}
	}

	return 0;
}

int cCorrAnalysis::calcHamming(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	xCorrThread
		*Xp_ct		= (xCorrThread *)Vp_shareData;
	char	**Na_data	= Xp_ct->Na_data;
	wsReal	**Ra_corr	= Xp_ct->Ra_corr;
	int		**Na_corrN	= Xp_ct->Na_corrN;
	wsUint	N_sample	= Xp_ct->N_sample;
	wsUint	N_SNP		= Xp_ct->N_variant;
	char*	Na_cm		= Xp_ct->Na_corrMask;
	wsFloatCst*	Ra_cm		= Xp_ct->Ra_corrMask;
	vVariant	&Xa_SNP		= *(Xp_ct->Xp_SNP);
	int		i			= *((int *)Vp_data);

	if (Xp_ct->B_isComplete) {
		for (wsUint j=i ; j<N_sample ; j++) {
			wsUint N_inc;
			wsReal R_res = _ibsComplete(Na_data[i], Na_data[j], N_SNP,
				Xa_SNP, &N_inc, Na_cm, Ra_cm);
			Ra_corr[i][j]	= R_res;
			Ra_corr[j][i]	= R_res;
			Na_corrN[i][j]	= N_inc;
			Na_corrN[j][i]	= N_inc;

		}
	}
	else {
		for (wsUint j=i ; j<N_sample ; j++) {
			wsUint N_inc;
			wsReal R_res = _ibsIncomplete(Na_data[i], Na_data[j], N_SNP,
				Xa_SNP, &N_inc, Na_cm, Ra_cm);
			Ra_corr[i][j]	= R_res;
			Ra_corr[j][i]	= R_res;
			Na_corrN[i][j]	= N_inc;
			Na_corrN[j][i]	= N_inc;

		}
	}

	return 0;
}

int cCorrAnalysis::calcBN1(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	xCorrThread	*Xp_ct		= (xCorrThread *)Vp_shareData;
	char**		Na_data		= Xp_ct->Na_data;
	wsReal**	Ra_corr		= Xp_ct->Ra_corr;
	int**		Na_corrN	= Xp_ct->Na_corrN;
	wsUint		N_sample	= Xp_ct->N_sample;
	wsUint		N_SNP		= Xp_ct->N_variant;
	char*		Na_cm		= Xp_ct->Na_corrMask;
	wsFloatCst*	Ra_cm		= Xp_ct->Ra_corrMask;
	vVariant	&Xa_SNP		= *(Xp_ct->Xp_SNP);
	int			i			= *((int *)Vp_data);

	if (Xp_ct->B_isComplete) {
		for (wsUint j=i ; j<N_sample ; j++) {
			wsUint N_inc;
			wsReal R_res = _ibsComplete(Na_data[i], Na_data[j], N_SNP,
				Xa_SNP, &N_inc, Na_cm, Ra_cm);
			Ra_corr[i][j]	= R_res;
			Ra_corr[j][i]	= R_res;
			Na_corrN[i][j]	= N_inc;
			Na_corrN[j][i]	= N_inc;

		}
	}
	else {
		for (wsUint j=i ; j<N_sample ; j++) {
			wsUint N_inc;
			wsReal R_res = _ibsIncomplete(Na_data[i], Na_data[j], N_SNP,
				Xa_SNP, &N_inc, Na_cm, Ra_cm);
			Ra_corr[i][j]	= R_res;
			Ra_corr[j][i]	= R_res;
			Na_corrN[i][j]	= N_inc;
			Na_corrN[j][i]	= N_inc;

		}
	}

	return 0;
}

void cCorrAnalysis::_export()
{
	wsUint		i, j;
	vSampPtr	&Xa_samples = Cp_IO->getSample();
	wsUint		N_sample = Cp_IO->sizeSample();
	char		S_fn[32];

	if (!IS_ASSIGNED(makecor)) return;

	/* --corpair & --corgrm/--corepacts is m.e. */
	if (OPT_ENABLED(corpair) && (OPT_ENABLED(corgrm) || OPT_ENABLED(corepacts)))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--corpair", "--corgrm/--corepacts");

	if (OPT_ENABLED(ibs))
		strcpy(S_fn, "ibs.cor");
	else if (OPT_ENABLED(medcor))
		strcpy(S_fn, "empi.med.cor");
	else if (OPT_ENABLED(ktau))
		strcpy(S_fn, "ktau.cor");
	else if (OPT_ENABLED(empktau))
		strcpy(S_fn, "empi.ktau.cor");
	else if (OPT_ENABLED(corpearson))
		strcpy(S_fn, "pearson.cor");
	else
		strcpy(S_fn, "empi.cor");

	if (OPT_NUMBER(corpair) || !stricmp(OPT_STRING(makecor), "pair")) {
/**/	cExporter*	Cp_exporter		= cExporter::summon(S_fn); /* CHECKED */
		wsMat		Ra_corr			= getCorr();

		Cp_exporter->put("IID1	IID2	CORR	N\n");

		if (OPT_ENABLED(founderonly)) {
			const char *Ba_isFounder = Cp_IO->getIsFounder();

			wsUint _i, _j;
			if (OPT_ENABLED(indep)) for (i=_i=0 ; i<N_sample ; i++) {
				if (Ba_isFounder[i] == 0) continue;

				for (j=_j=0 ; j<i ; j++) {
					if (Ba_isFounder[j] == 0) continue;

					Cp_exporter->fmt("%s	%s	0.0	%d\n",
						Xa_samples[_i]->S_IID.c_str(),
						Xa_samples[_j]->S_IID.c_str(),
						N_sample);
					_j++;
				}
				Cp_exporter->fmt("%s	%s	1.0	%d\n",
					Xa_samples[_i]->S_IID.c_str(),
					Xa_samples[_i]->S_IID.c_str(), N_sample);
				_i++;
			} else for (i=_i=0 ; i<N_sample ; i++) {
				if (Ba_isFounder[i] == 0) continue;

				for (j=_j=0 ; j<=i ; j++) {
					if (Ba_isFounder[j] == 0) continue;

					if (Ra_corr[_i][_j] != Ra_corr[_i][_j])
						Cp_exporter->fmt("%s	%s	NA	%d\n",
							Xa_samples[_i]->S_IID.c_str(),
							Xa_samples[_j]->S_IID.c_str(),
							Na_corrN[_i][_j]);
					else
						Cp_exporter->fmt("%s	%s	%g	%d\n",
							Xa_samples[_i]->S_IID.c_str(),
							Xa_samples[_j]->S_IID.c_str(),
							Ra_corr[_i][_j], Na_corrN[_i][_j]);
					_j++;
				}
				_i++;
			}
		} else {
			if (!OPT_ENABLED(indep)) for (i=0 ; i<N_sample ; i++)
				for (j=0 ; j<=i ; j++) {
					if (Ra_corr[i][j] != Ra_corr[i][j])
						Cp_exporter->fmt("%s	%s	NA	%d\n",
							Xa_samples[i]->S_IID.c_str(),
							Xa_samples[j]->S_IID.c_str(),
							Na_corrN[i][j]);
					else
						Cp_exporter->fmt("%s	%s	%g	%d\n",
							Xa_samples[i]->S_IID.c_str(),
							Xa_samples[j]->S_IID.c_str(),
							Ra_corr[i][j], Na_corrN[i][j]);
				}
			else for (i=0 ; i<N_sample ; i++) {
				for (j=0 ; j<i ; j++)
					Cp_exporter->fmt("%s	%s	0.0	%d\n",
						Xa_samples[i]->S_IID.c_str(),
						Xa_samples[j]->S_IID.c_str(),
						N_sample);
				Cp_exporter->fmt("%s	%s	1.0	%d\n",
					Xa_samples[i]->S_IID.c_str(),
					Xa_samples[i]->S_IID.c_str(),
					N_sample);
			}
		}

		delete Cp_exporter;
	} else if (OPT_ENABLED(corepacts) || !stricmp(OPT_STRING(makecor), "epacts")) {
		cExporter*	Cp_es		= cExporter::summon("kin", 0, ET_BIN);
		vSampPtr&	Xv_smp		= Cp_IO->getSample();
		wsMat		Ra_corr		= getCorr();

		/* Write out header */
		wsUint N_empty = 0;
		Cp_es->write("EMMA_KIN", 8);
		Cp_es->write(&N_sample, 4);
		Cp_es->write(&N_empty, 4);

		/* Write out values */
		LOOP (i, N_sample) {
			for (int j=(int)i ; j>=0 ; j--) {
				double r = (double)Ra_corr[i][j];
				Cp_es->write(&r, sizeof(double));
			}
		}

		/* Write labels */
		double v = 1.0;
		Cp_es->write(&v, sizeof(double));
		FOREACH (vSampPtr_it, Xv_smp, i) {
			wsUint L_samp = (wsUint)(*i)->S_IID.size();
			Cp_es->write(&L_samp, 4);
			Cp_es->write((*i)->S_IID.c_str(), L_samp);
		}

		delete Cp_es;
	} else if (OPT_ENABLED(corgrm) || !stricmp(OPT_STRING(makecor), "gcta")) {
		cExporter*	Cp_grmbin	= cExporter::summon("grm.bin", 0, ET_BIN);
		wsMat		Ra_corr		= getCorr();

		/* Write out corr val */
#ifdef USE_DBL
		wsFloat*		Ra_tmp		= NULL;
		sseMalloc(Ra_tmp, wsFloat, N_sample);
#endif
		for (wsUint i=0 ; i<N_sample ; i++) {
#ifdef USE_DBL
			wsFloat*	Ra_ocor		= Ra_tmp;
			for (wsUint j=0 ; j<=i ; j++)
				Ra_ocor[j]		= (wsFloat)Ra_corr[i][j];
#else
			wsFloat*	Ra_ocor		= Ra_corr[i];
#endif
			Cp_grmbin->write(Ra_ocor, sizeof(wsFloat)*(i+1));
		}
#ifdef USE_DBL
		sseFree(Ra_tmp);
#endif
		delete Cp_grmbin;

		/* Write out N */
		cExporter *Cp_grmN = cExporter::summon("grm.N");
		for (wsUint i=0 ; i<N_sample ; i++)
			Cp_grmN->write(Na_corrN[i], sizeof(wsUint)*(i+1));
		delete Cp_grmN;

		/* Write out IDs */
		cExporter*	Cp_grmID	= cExporter::summon("grm.id");
		vSampPtr&	Xv_samp		= getIO()->getSample();
		FOREACH (vSampPtr_it, Xv_samp, I)
			Cp_grmID->fmt("%s	%s\n", (*I)->S_FID.c_str(), (*I)->S_IID.c_str());
		delete Cp_grmID;
	} else {
		/* Name as vStr */
		vStr Sa_rowNames, Sa_colNames;
		if (!stricmp(OPT_STRING(makecor), "rvtests")) {
			Sa_colNames.push_back("FID");
			Sa_colNames.push_back("IID");
			FOREACH(vSampPtr_it, Xa_samples, i) {
				Sa_rowNames.push_back((*i)->S_FID + "\t" + (*i)->S_IID);
				Sa_colNames.push_back((*i)->S_IID);
			}
		} else FOREACH(vSampPtr_it, Xa_samples, i)
			Sa_rowNames.push_back((*i)->S_IID);

		/* Export as matrix form */
		Cp_corr->file(S_fn, &Sa_rowNames, 10);
	}
}

} // End namespace ONETOOL
