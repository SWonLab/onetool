#include <limits>
#include "analyses/score.h"
#include "analyses/scoreGen.h"
#include "analyses/scoreRao.h"
#include "analyses/emai.h"
#include "analyses/fam.h"
#include "analyses/annot.h"
#include "utils/matrix.h"
#include "utils/stat.h"
#define SCORE_THREAD_INDEPXT
using namespace std;

namespace ONETOOL {

/* !!! FARVAT-specific analysis !!! */
#if ((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE & FUNCTION_FARVAT) || \
	(TOOLSET_TYPE & FUNCTION_MFARVAT) || (TOOLSET_TYPE & FUNCTION_FARVATX)

typedef struct _xResScore {
	wsUint	N_miss;
	wsReal	R_gscStat, R_gscPval;
	wsReal	R_raoStat, R_raoPval;
	wsReal	R_lrtStat, R_lrtPval;
} xResScore;

typedef struct _xInpScore
{
	wsUint		N_idx, N_miss;
	cVector		*Vp_snp;
	cStdMatrix	*Mp_EVX;
	cStdMatrix	*Mp_WtP;
	cStdMatrix	**Ma_famWtP;
	cDiagMatrix	*Ma_famVinv;
	char		B_test;
} xInpScore;

void _gen(wsUint N_colX, cVector &V_XtP, cVector &V_RtP, size_t N_fam,
	wsUint *Na_seq, cStdMatrix **Ma_famWtP, cDiagMatrix *Ma_famVinv,
	wsReal *Rp_stat, wsReal *Rp_pVal)
{
	cVector		V_XVW(N_colX);
	cStdMatrix	M_WVW(N_colX, N_colX);
// 		d11 <- 0
// 		s1 <- 0
// 		d12tran <- 0
// 		d22 <- 0
	wsReal		R_d11		= W0;
	wsReal		R_s1		= W0;
	cVector		V_d12t(N_colX);
	cSymMatrix	M_d22(N_colX);

	/* Split V_XtP */
	wsUint		N_testSamp	= V_XtP.size();
	cVector		V_famXtP(N_testSamp);
	wsReal*		Ra_famXtP	= V_famXtP.get();
	cVector		V_famRtP(N_testSamp);
	wsReal*		Ra_famRtP	= V_famRtP.get();
	wsReal*		Ra_XtP		= V_XtP.get();
	wsReal*		Ra_RtP		= V_RtP.get();
	for (wsUint i=0 ; i<N_testSamp ; i++) {
		Ra_famXtP[i] = Ra_XtP[Na_seq[i]];
		Ra_famRtP[i] = Ra_RtP[Na_seq[i]];
	}

	// [1]	F <- length(fam_num)
	// for(i in 1:F){
	wsUint I=0, J;
	for (size_t iii=0 ; iii<N_fam ; iii++) {
		wsUint N_sz = Ma_famVinv[iii].row();

		// phi_list <- subset of correlation structure
		// #mem_i	<- # of members of ith family
// 		cStdMatrix	M_pWtP(N_colX, N_sz);
// 		wsReal	**Ra_famXt	= M_pWtP.get();

		// [#mem_i*#snp]
		// fam_n[[i]] <- snp[(cut[i]+1):cut[i+1]]
		// Equivalent to Ra_snp

		// [#mem_i*#inp]	fam_x[[i]] <- x[(cut[i]+1):cut[i+1],]
		J = I;
// 		for (i=0 ; i<N_colX ; i++)
// 			memcpy(Ra_famXt[i], M_WtP.get()[i]+J, sizeof(wsReal)*N_sz);
		I += N_sz;

		// [#mem_i*#mem_i]	cov <- sig2g*phi_list[[i]]+sig2*diag(fam_num[i])
		// [#mem_i*#mem_i]	Invcov <- solve(cov)
		//cDiagMatrix	M_pVinv		= M_Vinv.subset(J, N_sz);

		// [#snp*#mem_i]%*%[#mem_i*#mem_i]
		// [#snp*#mem_i]	ntInvcov <- t(fam_n[[i]])%*%Invcov
		cVector		V_pXtP		= V_famXtP.subset(J, N_sz);
		cVector		V_pXpVi		= V_pXtP * Ma_famVinv[iii];

		// [#snp*#mem_i]%*%[#mem_i*#pheno]
		// [#snp*#pheno]	nve[[i]] <- ntInvcov%*%fam_e[[i]]
		cVector		V_pRtP		= V_famRtP.subset(J, N_sz);
		wsReal		R_curXVe	= V_pXpVi.sum(V_pRtP);

		// [#snp*#mem_i]%*%[#mem_i*#inp]
		// [#snp*#inp]		cur_nvx <- ntInvcov%*%fam_x[[i]]
		cVector		V_curXpViWp	= V_pXpVi.Mt(*(Ma_famWtP[iii]));

		// [#snp*#inp]		nvx <- nvx + cur_nvx
		V_XVW					+= V_curXpViWp;

		// [#inp*#mem_i]%*%[#mem_i*#mem_i]
		// [#inp*#mem_i]	xtInvcov <- t(fam_x[[i]])%*%Invcov
		cStdMatrix	M_WpVi		= *(Ma_famWtP[iii]) * Ma_famVinv[iii];

		// [#inp*#mem_i]%*%[#mem_i*#inp]
		// [#inp*#inp]		cur_xvx <- xtInvcov%*%fam_x[[i]]
		cStdMatrix	M_curWpViWp	= M_WpVi.Mt(*(Ma_famWtP[iii]));
		// 		wsReal **Ra_curXtVx = sseMpMt(Ra_xtInvcov, N_colX, N_sz,
		// 			Ra_famXt, N_colX, N_sz); /* 120917 */
		// 		sseUnmat(Ra_famXt, N_colX);

		// [#inp*#inp]		xvx <- xvx + cur_xvx
		M_WVW					+= M_curWpViWp;

		// [#inp*#mem_i]%*%[#mem_i*#pheno]
		// [#inp*#pheno] xve[[i]] <- xtInvcov%*%fam_e[[i]]
		cVector	V_curWVe		= M_WpVi * V_pRtP;

		// + [#snp*#pheno]%*%[#pheno*#snp]
		// [#snp*#snp]	d11 <- d11+nve[[i]]%*%t(nve[[i]])
		R_d11					+= SQR(R_curXVe);

		// [#snp*#snp]	s1 <- s1+nve[[i]]
		R_s1					+= R_curXVe;

		// FIXME : [#inp*#pheno]%*%[#snp*#pheno] in code of Meiling
		// Maybe should be t(nve[[i]])
		// 
		// [#inp*#pheno]%*%[#snp*#pheno]
		// [#inp*#snp]	d12tran <- d12tran+xve[[i]]%*%nve[[i]]
		cVector	V_WVeXVe		= V_curWVe * R_curXVe;
		V_d12t					+= V_WVeXVe;

		// [#inp*#pheno]%*%[#pheno*#inp]
		// [#inp*#inp]	d22 <- d22+xve[[i]]%*%t(xve[[i]])
		cSymMatrix M_WVe2		= V_curWVe.tV();
		M_d22					+= M_WVe2;
	}
	if (I != N_testSamp)
		halt("SYSERR: Sum of available family member expected %d, counted %d",
		N_testSamp, I);
	//	sseFree(Ra_snpSub);
#ifdef EMAIvalidate
	exportMatrix("ea.gen.xtVx", Ra_xtVx, N_colX, N_colX);//EMAIvalidate
	exportMatrix("ea.gen.nvx", Ra_ntVx, 1, N_colX);//EMAIvalidate
	exportMatrix("ea.gen.d22", Ra_d22, N_colX, N_colX);//EMAIvalidate
#endif

	// [#inp*#inp]		Invxvx	<- solve(xvx)
	cStdMatrix&	M_WVWi		= M_WVW.inv();

	// [#snp*#inp]%*%[#inp*#inp]
	// [#snp*#inp]		nvxInvxvx <- nvx%*%Invxvx
	cVector		V_XVWWVWi	= V_XVW * M_WVWi;

	// [#snp*#inp]%*%[#inp*#inp]%*%[#inp*#snp]
	// [#snp*#snp]		IID		<- nvxInvxvx%*%d12tran
	wsReal		R_IID		= V_XVWWVWi.sum(V_d12t);

	// [#snp*#inp] %*% [#inp*#inp]
	// [#snp*#inp]		nvxInvxvxd22 <- nvxInvxvx %*% d22
	cVector		V_XVWWVWiD	= V_XVWWVWi * M_d22;

	// [#snp*#inp] %*% [#inp*#inp] %*% [#inp*#snp]
	// [#snp*#snp]		IIDII <- nvxInvxvxd22 %*% Invxvx %*% nvx_t
	wsReal		R_IIDII		= V_XVWWVWiD.qf(M_WVWi, V_XVW);
	delete &M_WVWi;

	// [#snp*#snp] - [#snp*#snp] - [#snp*#snp] + [#snp*#snp]
	// In case of #snp=1, sx <- d11-2*IID+IIDII cause IID=t(IID)
	// [#snp*#snp]		sx <- d11-IID-t(IID)+IIDII
	R_d11					-= R_IID + R_IID - R_IIDII;

	// [#snp*#snp] %*% [#snp*#snp] %*% [#snp*#snp]
	// [#snp*#snp]		Tgs <- t(s1) %*% solve(sx) %*% s1
	wsReal		R_Tgs		= R_s1 * (W1/R_d11) * R_s1;

	// FIXME : How to calculate the p-value under chisq dist. with Tgs?
	// In case of #snp>1, Tgs is MATRIX with [#snp*#snp] dim.
	// Spvalue <- 1-pchisq(Tgs, #snp)
	wsReal		R_pValue	= (wsReal)PVchisq(R_Tgs, 1.0);

	*Rp_stat = R_Tgs;
	*Rp_pVal = R_pValue;
}

void _rao(cVector &V_XtP, cVector &V_RtP, cStdMatrix &M_WtP, cSymMatrix &M_WViWi,
	cDiagMatrix &M_Vinv, wsReal *Rp_stat, wsReal *Rp_pVal)
{
	//	wsUint	i;

	// [#snp*#samp]%*%[#samp*#samp]
	// [#snp*#samp]		snpInvcov <- t(snp)%*%Invcov
	cVector	V_XVi	= V_XtP * M_Vinv;

#ifdef EMAIvalidate
	exportMatrix("em.rao.snpInvcov", Ra_snpVinv, 1, N_samp);//EMAIvalidate
#endif

	// [#snp*#samp]%*%[#samp*#pheno]
	// [#snp*#pheno]	nve <- snpInvcov%*%e 
	wsReal	R_xVe	= V_XVi.sum(V_RtP);

	// [#snp*#samp]%*%[#samp*#snp]
	// [#snp*#snp]		nvn <- snpInvcov%*%snp
	wsReal	R_xVx	= V_XVi.sum(V_XtP);

	// [#snp*#inp]		nvx <- snpInvcov%*%x
	cVector	V_xVW	= V_XVi.Mt(M_WtP);

#ifdef EMAIvalidate
	exportMatrix("em.rao.nvX", Ra_nvX, 1, N_colX);//EMAIvalidate
	//exportMatrix("em.rao.XVX", Ra_XVX, N_colX, N_colX);//EMAIvalidate
#endif

	// Trs <- nve%*%solve(nvn-nvx%*%solve(xvx)%*%t(nvx))%*%t(nve)
	wsReal R_Trs = SQR(R_xVe) / (R_xVx - V_xVW.qf(M_WViWi));

	// RSpvalue <- 1-pchisq(Trs,1)
	wsReal R_pValue = (wsReal)PVchisq(R_Trs, 1.0);
	
	*Rp_stat	= R_Trs;
	*Rp_pVal	= R_pValue;
}

int thr_scoreTest(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	cTestScoreAnalysis	*Cp_anaSC	= (cTestScoreAnalysis *)Vp_shareData;		///< Shared data is its instance
//	cIO					*Cp_IO		= Cp_anaSC->getIO();
	cDiagMatrix			&M_Vinv		= Cp_anaSC->getVinv();
//	wsUint				*Na_szFam	= Cp_anaSC->getFamSize();
	cVector&			V_residP	= Cp_anaSC->getResidP();
	cSymMatrix&			M_WViWi		= Cp_anaSC->getWViWi();
	wsReal*				Ra_Ts		= Cp_anaSC->getTs();
	xInpScore			*Xp_inp		= (xInpScore *)Vp_data;
	xResScore			*Xp_res		= ((xResScore *)Vp_result) + Xp_inp->N_idx;
	int					N_idxVrt	= *((int *)Vp_data);
	if (Xp_inp->B_test == 0) {
		Xp_res->R_raoStat = WISARD_NAN;
		Xp_res->R_raoPval = WISARD_NAN;
		if (Ra_Ts) Ra_Ts[N_idxVrt] = Xp_res->R_raoStat;
		return 0;
	}
#ifdef USE_ED2
	cVector				V_XtP		= Xp_inp->Vp_snp->Mt(*(Xp_inp->Mp_EVX));
#else
	cVector				V_XtP		= *(Xp_inp->Vp_snp) *(Xp_inp->Mp_EVX);
#endif

	Xp_res->R_gscStat	= numeric_limits<wsReal>::infinity();
	Xp_res->R_gscPval	= numeric_limits<wsReal>::infinity();
	Xp_res->N_miss		= Xp_inp->N_miss;
	if (OPT_ENABLED(kinship)) {
		wsUint*			Na_seq		= Cp_anaSC->getSampSeqByFam();
		/* FIXME : Reimpl */
		_gen(M_WViWi.row(), V_XtP, V_residP, Cp_anaSC->getIO()->getFamilyData().size(),
			Na_seq, Xp_inp->Ma_famWtP, Xp_inp->Ma_famVinv,
			&Xp_res->R_gscStat, &Xp_res->R_gscPval);
	}

	/* Do Rao's Score test */
	_rao(V_XtP, V_residP, *(Xp_inp->Mp_WtP), M_WViWi, M_Vinv,
		&Xp_res->R_raoStat, &Xp_res->R_raoPval);
	if (Ra_Ts) Ra_Ts[N_idxVrt] = Xp_res->R_raoStat;

	return 0;
}

int forAllSNP_score(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	static wsUint		N_idxVrt;
	static cVector**	Va_snps		= NULL;
	cTestScoreAnalysis*	Cp_anaScr	= (cTestScoreAnalysis *)Vp_shareData;
	xInpScore*			Xp_inp		= (xInpScore *)Vp_thrData;
	cIO*				Cp_IO		= Cp_anaScr->getIO();
	xMaf*				Xp_maf		= Cp_IO->getMAF();
	wsUintCst				N_testSamp	= Cp_anaScr->getTestSampSize();
	const char*			Ba_excl		= Cp_anaScr->getExcl();
//	wsUint				*Na_mSamp	= Cp_anaScr->getMapSample();
	wsUintCst				N_vrt		= Cp_IO->sizeVariant();
//	char				**Na_data	= Cp_IO->getGenotype();
	static cStdMatrix	*Ma_WtP		= NULL;
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		wsAlloc(Ma_WtP, cStdMatrix, N_thread);
		for (int ii=0 ; ii<N_thread ; ii++)
			Ma_WtP[ii] = Cp_anaScr->getWtP().clone();
		wsAlloc(Va_snps, cVector *, N_thread);
		for (int ii=0 ; ii<N_thread ; ii++)
			Va_snps[ii] = new cVector(N_testSamp);
		N_idxVrt = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for (i=0 ; i<N_thread ; i++) {
			/* Checking border */
			if ((N_idxVrt+i) >= N_vrt) break;

			Xp_inp[i].N_idx			= N_idxVrt+i;
			Xp_inp[i].B_test		= 0;
			if (Xp_maf[N_idxVrt+i].R_allMaf == W0)
				continue;

			/* Copy SNP data into Ra_snp */
			wsReal*	Ra_snp	= Va_snps[i]->get();
			wsUint	N_miss	= setGeno(Cp_IO, N_idxVrt+i, Ra_snp,
				Ba_excl, 1);

			Xp_inp[i].B_test		= 1;
			Xp_inp[i].Vp_snp		= Va_snps[i];
			Xp_inp[i].N_miss		= N_miss;
			Xp_inp[i].Mp_EVX		= Cp_anaScr->getEVX();
			Xp_inp[i].Mp_WtP		= &(Ma_WtP[i]);
			Xp_inp[i].Ma_famVinv	= Cp_anaScr->getFamVinv();
			Xp_inp[i].Ma_famWtP		= Cp_anaScr->getFamWtP();
			if (((N_idxVrt+i)%1000) == 0)
				notice("%d/%d SNPs processed\r", N_idxVrt+i, N_vrt);
		}
		N_idxVrt += N_thread;
		break;
	case DISTFUN_UNINIT:
		for (int ii=0 ; ii<N_thread ; ii++)
			delete Va_snps[ii];
		delete [] Ma_WtP;
		DEALLOC(Va_snps);
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

// int getFamMember(xSample *Xp_samp, void *Vp_data)
// {
// 	map<string,int> *Xp_famMem = (map<string,int> *)Vp_data;
// 	if (Xp_samp->N_idx != SAMP_NODATA)
// 		Xp_famMem->insert(make_pair(Xp_samp->S_IID, Xp_samp->N_idx));
// 	return 0;
// }

struct xCntSamp
{
	char *Ba_incl;
	wsUint N_cnt;
};

int countActualFamSize(xSample *Xp_samp, void *Vp_data)
{
	xCntSamp *Xp_cntSamp = (xCntSamp *)Vp_data;

	if (Xp_samp->N_idx != SAMP_NODATA) {
		if (Xp_cntSamp->Ba_incl[Xp_samp->N_idx] == 1)
			Xp_cntSamp->N_cnt++;
	}

	return 0;
}

cTestScoreAnalysis::cTestScoreAnalysis(cIO *Cp_IO, cAnalysis *Cp_inpAnaPhi,
	cFamStrAnalysis *Cp_inpAnaFS, char B_inpIsREML) : cAnalysis(Cp_IO)
{
	Cp_res0EmAi	= NULL;
	Cp_anaPhi	= Cp_inpAnaPhi;
	Cp_anaFS	= Cp_inpAnaFS;
	B_isREML	= B_inpIsREML;
	Mp_EV		= NULL;
	Mp_EVX		= NULL;
	Mp_Yt		= NULL;
	Ra_phi		= NULL;
	Ra_stat		= NULL;
	Mp_Vinv		= NULL;
	Na_seq		= NULL;
//	Ra_Vinv		= NULL;
// 	Ra_X		= NULL;
// 	Ra_Xt		= NULL;
//	Ra_residual	= NULL;

	wsUint N_sample		= Cp_IO->sizeSample();
	wsUint N_pheno		= Cp_IO->sizePheno();
// 	wsUint N_cov		= Cp_IO->getCovariateSize();	
// 	wsReal **Ra_cov		= Cp_IO->getCovariates();
// 	wsReal **Ra_pheno	= Cp_IO->getPhenos();
	Ba_excl				= NULL;

	/* Which sample should be included? */
	wsUint i, j;

	/* Check its phenotype and covariates */
	Ba_excl = Cp_IO->getPheCovMissing(&N_testSamp);
	/* If --kinship, get famIdx and resize famIdx */
	wsUint*	Na_testFamIdx = NULL;
	if (OPT_ENABLED(kinship)) {
		wsUint*	Na_famIdx		= Cp_anaFS->getFamIndices();
		wsAlloc(Na_testFamIdx, wsUint, N_testSamp);
		for (wsUint k=0,K=0 ; k<N_sample ; k++)
			if (!Ba_excl[k])
				Na_testFamIdx[K++] = Na_famIdx[k];
	}

	/* 0 check */
	if (N_testSamp == 0) {
		LOG("No available sample left, can't do this analysis\n");
		return;
	}
		
	if (N_sample == N_testSamp)
		LOG("No missing on both of covariates and phenotype, not resized\n");
	else
		LOG("Target is resized to %d->%d samples\n", N_sample, N_testSamp);

	/* Get phi matrix */
	wsReal **Ra_origPhi = getFullCorMat(Cp_anaPhi);

	/* Rearrange phi matrix into family-wise */
// 	MULTI_MALLOC(Na_sampIdx, wsUint, N_testSamp); /* 120917 */
// 	Ra_phi = sseEmptyMatrix(N_testSamp, N_testSamp); /* 120917 */

// 	vSampPtr	&Xa_samp	= Cp_IO->getSamplesV2();

#ifdef EMAIvalidate
	mFam		&Xa_fam		= Cp_IO->getFamilyData();
/**/cExporter* E = cExporter::summon("ea.scr.famsz");//EMAIvalidate
	i = 0;
	FOREACHDO (mFam_it, Xa_fam, fit, i++) {
		E->fmt("%d ", Na_szFam[i]);
	}
	delete E;
#endif

	/* Resize phi matrix */
	if (Ra_origPhi != NULL)
		Ra_phi = sseMsubsetRect(Ra_origPhi, N_sample, Ba_excl, 0, N_testSamp);
// 		for (i=0 ; i<N_testSamp ; i++)
// 			for (j=0 ; j<N_testSamp ; j++)
// 				Ra_phi[i][j] = Ra_origPhi[Na_sampIdx[i]][Na_sampIdx[j]];
	else {
		Ra_phi = sseEmptyMatrix(N_testSamp, N_testSamp);
		for (i=0 ; i<N_testSamp ; i++)
			Ra_phi[i][i] = W1;
	}
#ifdef EMAIvalidate
	exportMatrix("ea.phi", Ra_phi, N_testSamp, N_testSamp);
#endif
// 	if (!OPT_ENABLED(specdcmp)) {
// 		if (Cp_IO->isTwinAssigned() || OPT_ENABLED(ginv))
// 			SVDinverse(Ra_phi, N_testSamp, &Ra_phiInv); /* 120917 */
// 		else {
// 			if (OPT_ENABLED(kinship))
// 				Ra_phiInv = Cp_IO->getFamInv(Na_sampIdx, N_testSamp, Ra_phi); /* 120917 */
// 			else
// 				Ra_phiInv = invSymMat(Ra_phi, N_testSamp); /* 120917 */
// 			if (Ra_phiInv == NULL)
// 				LOG("Can't get inverse of phi matrix, spectral decomposition will be used\n");
// 		}
// 		LOG("phiInv calculated\n");
// 	} else {
// 		LOG("Spectral decomposition will be used\n");
// 		Ra_phiInv = NULL;
// 	}

	/* Do eigendecomposition if required */
	cTimer tt;
	tt.start();

	wsReal	**Ra_eVec	= NULL;
	wsReal	*Ra_eVal	= NULL;
	getEigenPhi(Ra_phi, N_testSamp, &Ra_eVec, &Ra_eVal, Na_testFamIdx);
	DEALLOC(Na_testFamIdx);
	Mp_EV		= new cDiagMatrix(N_testSamp, Ra_eVal);
	/* Mp_EVX = P^t, not P */
	Mp_EVX		= new cStdMatrix(N_testSamp, N_testSamp, Ra_eVec);

	LOG("Eigendecomposition [%s]\n", tt.getReadable());

#ifdef EMAIvalidate
	exportMatrix("ea.scr.phi", Ra_phi, N_testSamp, N_testSamp);//EMAIvalidate
// 	if (Ra_phiInv != NULL)
// 		exportMatrix("ea.scr.phiInv", Ra_phiInv, N_testSamp, N_testSamp);//EMAIvalidate
#endif

	/* Rearrange y matrix */
	Mp_Yt = new cStdMatrix(N_pheno, N_testSamp);
	wsMat		Ra_Y		= Mp_Yt->get();//sseMatrix(N_pheno, N_testSamp); /* 120917 */
	wsMat		Ra_origY	= Cp_IO->getPhenos();
	for (j=0 ; j<N_pheno ; j++) {
		wsUint		I = 0;
		for (i=0 ; i<N_sample ; i++)
			if (!Ba_excl[i])
				Ra_Y[j][I++] = Ra_origY[j][i];
	}
	//}
#ifdef USE_ED2
	M_YtP	= Mp_Yt->Mt(*Mp_EVX);
#else
	M_YtP	= *Mp_Yt * *Mp_EVX;
#endif
//	LOG("phenotype rearranged\n");

#ifdef EMAIvalidate
	exportMatrix("ea.scr.y.t", Ra_Y, N_pheno, N_testSamp);//EMAIvalidate
#endif

	/* Perform NULL test */
	_test0();
}

cTestScoreAnalysis::~cTestScoreAnalysis()
{
	delete Mp_EV;
	delete Mp_EVX;
	delete Mp_Yt;
	sseUnmat(Ra_phi, N_testSamp);
	M_YtP.rem();
	DEALLOC(Ba_excl);
	sseFree(Ra_stat);
	DEALLOC(Na_seq);

	if (Cp_res0EmAi)	delete Cp_res0EmAi;
	if (Mp_Vinv)		delete Mp_Vinv;

	/* Deallocate static element of Rao's score test */
	cTestRaoAnalysis::M_WViWi.rem();
}

void cTestScoreAnalysis::run()
{
	if (N_testSamp == 0) return;

	vVariant	Xv_vrt			= Cp_IO->getVariant();
	xMaf*		Xp_maf			= Cp_IO->getMAF();
	wsUint		N_vrt			= Cp_IO->sizeVariant();
//	MAT_t		Ra_cov			= Cp_IO->getCovariates();
	wsUint		N_cov			= Cp_IO->sizeCovar();
	char		S_bufLRT[256]	= { 0, };

	/* Do not perform score test w/o --scoretest */
	if (!OPT_ENABLED(scoretest)) {
//		LOG("Score test will not performed due to --heritability or --makeblup were only assigned\n");
		return;
	}

	if (!OPT_ENABLED(kinship))
		LOGwarn("--kinship is disabled, generalized score test will not performed\n");

	/* Set 1st ~ N_cov+1th column */
	wsMat	Ra_Xt	= getIO()->getFiltNullDesignMat(Ba_excl, 1, N_testSamp, 1);
	wsUint	N_col	= N_cov+1;
	wsUint	N_samp	= Cp_IO->sizeSample();
	M_Xt.init(MATOPT_DELNONE, N_col, N_testSamp, WISARD_NAN, Ra_Xt);

	/* Print header */
/**/cExporter*	Cp_score = cExporter::summon("scoretest.res"); /* 140109 CONFIRMED */
	LOGoutput("Result of score test is exported to [%s.scoretest.res]\n",
		OPT_STRING(out));
	headerVariant(Cp_score);
	Cp_score->put("	NMISS");
	if (OPT_ENABLED(kinship))
		Cp_score->put("	STAT_GEN	P_GEN");
	if (OPT_ENABLED(lrt)) sprintf(S_bufLRT, "	STAT_LRT	P_LRT");
	Cp_score->fmt("	STAT_RAO	P_RAO%s\n", S_bufLRT);

	/* Calculate invV using sig2, sig2g */
	wsReal R_sig2	= Cp_res0EmAi->getEstimates(0)[0];
	wsReal R_sig2g	= Cp_res0EmAi->getEstimates(0)[1];

	// [#samp*#samp]	Invcov <- solve(V)
	cDiagMatrix &M_VEV = *Mp_EV * R_sig2g;
	M_VEV += R_sig2;
	cDiagMatrix &M_Vinv = M_VEV.inv();
	Mp_Vinv = &M_Vinv;
	delete& M_VEV;

	// [#samp*#pheno] - [#samp*#inp]%*%[#inp*#pheno]
	// [#samp*#pheno]	e <- pheno-x%*%alpha
	/* Get x%*%alpha */
	wsReal*	Ra_beta	= Cp_res0EmAi->getBeta()[0];
	wsRealCst* Ra_blupP = Cp_res0EmAi->getBLUP();
	wsReal* Ra_blup = sseVector(N_testSamp);
	memcpy(Ra_blup, Ra_blupP, sizeof(wsReal)*N_testSamp);
	cVector V_blup(Ra_blup, N_testSamp);
	cVector	V_beta(Ra_beta, N_col, 1);
	cVector	V_pred	= M_Xt.tV(V_beta);
#ifdef USE_ED2
	M_WtP			= M_Xt.Mt(*Mp_EVX);
#else
	M_WtP			= M_Xt * *Mp_EVX;
#endif

	/* Get e */
	cVector V_res	= Mp_Yt->r2v(0) - V_pred;
	V_res -= V_blup;
#ifdef USE_ED2
	V_residP		= V_res.Mt(*Mp_EVX);
#else
	V_residP		= V_res * *Mp_EVX;
#endif
	//sseVsV(Ra_Y[0], Ra_predicted, Ra_residual, N_testSamp);

#ifdef EMAIvalidate
//	exportMatrix("scr.residual", &Ra_residual, 1, N_testSamp);//EMAIvalidate
	V_res.file("scr.residual");
#endif
	V_res.rem();

	/* Calculate XVIXinv */
	cTimer t;
//	LOG("Calculating XVIXinv\n");
	/* This term requires only one-time calculation
	 * Because all X and Vi will be equal on all tests
	 * t(X) %*% Vinv %*% X */
	cSymMatrix	M_XViX = M_WtP.MMt(*Mp_Vinv);
// 	wsReal **Ra_XVX = sseMMMt(Ra_Xt, N_cov+1, N_testSamp,
// 		M_Vinv, Mp_EVX, N_testSamp, N_testSamp,
// 		Ra_Xt, N_cov+1, N_testSamp);
// 	cTestRaoAnalysis::Ra_XVIXinv = SO_invMatrix(Ra_XVX, N_cov+1);
// 	if (cTestRaoAnalysis::Ra_XVIXinv == NULL)
// 		halt("Can't get XVIXinv matrix");
	cSymMatrix& M_tmpWViWi = M_XViX.inv();
	M_WViWi.init(M_tmpWViWi);
	delete& M_tmpWViWi;
	M_XViX.rem();
//	deallocMatrix(Ra_XVX, N_cov+1, (void *)1);

	/* Divide if required */
	mFam&	Xm_fam	= getIO()->getFamilyData();
	wsUint	N_fam	= (wsUint)getIO()->getFamilyData().size();
	Ma_famWtP		= NULL;
	Ma_famVinv		= NULL;
	/* Family sequence */
	if (OPT_ENABLED(kinship)) {
		wsAlloc(Na_seq, wsUint, N_testSamp);
		wsAlloc(Ma_famWtP, cStdMatrix *, N_fam);
		wsAlloc(Ma_famVinv, cDiagMatrix, N_fam);
		wsUint	N_colX	= M_WtP.row();
		wsUint	I=0, J=0;

		/* Get sample's re-ordered sequence */
		wsUint*	Na_newSeq = NULL;
		wsAlloc(Na_newSeq, wsUint, N_samp);
		for (wsUint j=0 ; j<N_samp ; j++)
			Na_newSeq[j] = Ba_excl[j] ? 0xffffffff : J++;

		/* Make family-wise WtP and Vinv */
		FOREACHDO (mFam_it, Xm_fam, fit, I++) {
			wsUint		N_sz	= 0;
			xFamily&	X_fam	= fit->second;
			FOREACH (vInt_it, X_fam.Xv_members, j) {
				if (Ba_excl[*j]) continue;
				N_sz++;
			}

			Ma_famWtP[I] = new cStdMatrix(N_colX, N_sz);
			Ma_famVinv[I].init(N_sz);
		}

		/* For each family, split WtP and Vinv */
		I = 0;
		J = 0;
		FOREACHDO (mFam_it, Xm_fam, fit, I++) {
			wsMat	Ra_famWtP	= Ma_famWtP[I]->get();
			wsReal*	Ra_Vinv		= Ma_famVinv[I].get()[0];
			xFamily &X_fam		= fit->second;
			wsUint	l			= 0;

			/* Family-wise */
			FOREACHDO (vInt_it, X_fam.Xv_members, j, l++) {
				if (Ba_excl[*j]) continue;
	
				for (wsUint k=0 ; k<N_col ; k++)
					Ra_famWtP[k][l] = M_WtP.get()[k][*j];
				Ra_Vinv[l]	= Mp_Vinv->get()[0][*j];
				/* ith element of Na_seq = value'th row should be mapped */
				Na_seq[J++]	= Na_newSeq[*j];
			}
		}
		DEALLOC(Na_newSeq);
	}

	/* --genoctrl or --mutlticomp */
	if (OPT_ENABLED(genoctrl) || OPT_ENABLED(adjust)) {
		/* Allocate buffer to storing statistics */
		Ra_stat = sseVector(N_vrt);
		sseVinit(Ra_stat, N_vrt, WISARD_NAN);
	}
		

	t.start();
	wsUint N_prt = 0xffffffff;
	if (OPT_NUMBER(thread) > 1 && !OPT_ENABLED(lrt)) {
		/* If multithreaded */
		xResScore *Xa_resScore = NULL;
		wsAlloc(Xa_resScore, xResScore, N_vrt);
		WORKER().run(thr_scoreTest, forAllSNP_score, this, Xa_resScore, sizeof(xInpScore));

		/* Print output */
		for (wsUint i=0 ; i<N_vrt ; i++) {
			xResScore *Xp_r = Xa_resScore + i;
			/* --remna */
			if (OPT_ENABLED(remna) && NA(Xp_r->R_gscPval) && NA(Xp_r->R_raoPval))
				continue;
			/* --pvalrange */
			if (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), Xp_r->R_gscPval) &&
				!isInRange(OPT_RANGE(pvalrange), Xp_r->R_raoPval))
				continue;

			entryVariant(Cp_score, Xv_vrt[i]);
			Cp_score->fmt("	%d", Xp_r->N_miss);
			if (OPT_ENABLED(kinship))
				Cp_score->fmt("	%g	%g", Xp_r->R_gscStat, Xp_r->R_gscPval);
			Cp_score->fmt("	%g	%g\n", Xp_r->R_raoStat, Xp_r->R_raoPval);
			N_prt++;
		}
		DEALLOC(Xa_resScore);
	} else {
		/* Resize Xt matrix */
		M_Xt.init(MATOPT_DELNONE, N_col+1, N_testSamp, WISARD_NAN, Ra_Xt);

		N_prt = 0;
		for (wsUint i=0 ; i<N_vrt ; i++) {
			wsReal	R_gscStat, R_gscPval;
			wsReal	R_raoStat, R_raoPval;
			wsReal	R_lrtStat=WISARD_NAN, R_lrtPval=WISARD_NAN;
			wsUint	N_miss = 0;

			R_gscStat = WISARD_NAN;
			R_gscPval = WISARD_NAN;

			if (Xp_maf[i].R_allMaf != W0) {
				/* Set genotype */
				N_miss = setGeno(Cp_IO, i, Ra_Xt[N_cov+1],
					Ba_excl, 1);
#ifdef EMAIvalidate
				exportMatrix("ea.scr.x1", Ra_Xt, N_cov+2, N_testSamp);
#endif
				cVector	V_snp(Ra_Xt[N_cov+1], N_testSamp, 1);
#ifdef USE_ED2
				cVector V_XtP = V_snp.Mt(*Mp_EVX);
#else
				cVector V_XtP = V_snp * *Mp_EVX;
#endif

				/* Do EM-AI only --lrt */
				cEmAiAnalysisV2 *Cp_currEmAi = NULL;
				if (OPT_ENABLED(lrt)) {
					/* Make Full XT */
					cVector		V_Y(Mp_Yt->get()[0], N_testSamp, 1);
					cSymMatrix	M_phi(Ra_phi, N_testSamp, 1);
					cStdMatrix	M_Xt(N_col+1, N_testSamp, Ra_Xt, MATDEL_NONE);
					Cp_currEmAi = new cEmAiAnalysisV2(Cp_IO, M_Xt, V_Y, M_phi,
						B_isREML, Ba_excl, Mp_EVX, Mp_EV);
					Cp_currEmAi->quiet();
					Cp_currEmAi->run();
				}


				/* What entities are required to perform test?
				 * Adjusted size of phi matrix
				 * Family-wise adjusted phenotype vector
				 * Family-wise adjusted phi matrix
				 * Which samples were included to analysis compared to its original
				 * EM-AI estimates (sig2 and sig2g)
				 * Design matrix
				 * The column size of design matrix */

				/* Do generalized score test */
				if (OPT_ENABLED(kinship)) {
					/* FIXME : Reimpl */
					_gen(N_col, V_XtP, V_residP, N_fam, Na_seq, Ma_famWtP,
						Ma_famVinv, &R_gscStat, &R_gscPval);
				}

				/* Do Rao's Score test */
				_rao(V_XtP, V_residP, M_WtP, M_WViWi, *Mp_Vinv,
					&R_raoStat, &R_raoPval);

				/* --genoctrl */
				if (Ra_stat) Ra_stat[i] = R_raoStat;

				/* If --lrt, perform LRT */
				if (OPT_ENABLED(lrt) && !Cp_res0EmAi->getIsFailed())
					LRtest(Cp_currEmAi->getML(), Cp_res0EmAi->getML(),
						&R_lrtStat, &R_lrtPval, W2);
				if (Cp_currEmAi)
					delete Cp_currEmAi;
			} else {
				R_raoPval = R_raoStat = WISARD_NAN;
				R_gscPval = R_gscStat = WISARD_NAN;
				R_lrtPval = R_lrtStat = WISARD_NAN;
			}
			char B_na = NA(R_lrtPval) || NA(R_gscPval) || NA(R_raoPval);
			/* --remna */
			if (OPT_ENABLED(remna) && B_na) goto _next;

			/* --pvalrange */
			if (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), R_gscPval) &&
				!isInRange(OPT_RANGE(pvalrange), R_raoPval))
				goto _next;
			
			/* Make output */
			if (OPT_ENABLED(lrt)) {
				if (Cp_res0EmAi->getIsFailed())
					strcpy(S_bufLRT, "	NA	NA");
				else
					sprintf(S_bufLRT, "	%g	%g", R_lrtStat, R_lrtPval);
			} else S_bufLRT[0] = '\0';

			entryVariant(Cp_score, Xv_vrt[i]);
			Cp_score->fmt("	%d", N_miss);
			if (OPT_ENABLED(kinship)) {
				if (NA(R_gscPval)) Cp_score->put("	NA	NA");
				else Cp_score->fmt("	%g	%g", R_gscStat, R_gscPval);
			}
			if (NA(R_raoPval)) {
				Cp_score->fmt("	NA	NA%s\n", S_bufLRT);
			}
			else Cp_score->fmt("	%g	%g%s\n", R_raoStat, R_raoPval, S_bufLRT);
			N_prt++;
			/* Filtering log */
//			if (N_prt != N_vrt)
//				LOGnote("[%d] results were omitted by p-value filtering\n", N_vrt-N_prt);
_next:

			if ((i%1000) == 0)
				notice("%d/%d variants tested...\r", i, Xv_vrt.size());

#ifdef EMAIvalidate
			halt("Planned exit");//EMAIvalidate
#endif
		}
//		sseFree(Ra_snp);
	}
	/* Filtering log */
	if (N_prt != 0xffffffff && N_prt != N_vrt)
		LOGnote("[%d] results were omitted by p-value filtering\n", N_vrt-N_prt);
	LOG("[%d] varinats tested, [%s] taken\n", N_vrt, t.getReadable());
	delete Cp_score;

	/* --genoctrl, --adjust */
	if (Ra_stat)
		multcomp(getIO(), Ra_stat, N_vrt, "scoretest.adjust.res");

	if (OPT_ENABLED(kinship)) {
		for (wsUint i=0 ; i<N_fam ; i++)
			delete Ma_famWtP[i];
		DEALLOC(Ma_famWtP);
		DEALLOC(Ma_famVinv);
	}
	sseUnmat(Ra_Xt, N_col+1);
}

void cTestScoreAnalysis::_test0()
{
//	wsReal	**Ra_X;
	wsReal	**Ra_cov	= Cp_IO->getCovariates();
	wsUint	N_cov		= Cp_IO->sizeCovar();
	wsUint	N_pheno		= Cp_IO->sizePheno();
	wsUint	N_samp		= Cp_IO->sizeSample();
	wsUint	i, I, j;

	/* Allocate matrix */
// 	Ra_X	= sseMatrix(N_testSamp, N_cov+2);
// 	if (Ra_Xt)
// 		deallocMatrix(Ra_Xt, N_cov+2, (void *)1);
	wsMat	Ra_Xt		= sseMatrix(N_cov+1, N_testSamp);

	/* Make design matrix */
	for (i=0,I=0 ; i<N_samp ; i++) {
		if (Ba_excl[i] == 1) continue;

		//Ra_X[I][0] = W1;
		Ra_Xt[0][I] = W1;
		for (j=0 ; j<N_cov ; j++)
			Ra_Xt[j+1][I] = /*Ra_X[I][j+1] = */Ra_cov[j][i];
		I++;
	}
	//exportMatrix("X0", Ra_X, N_testSamp, N_cov+1);
	//exportMatrix("Yt", Ra_y, 1, N_testSamp);
// 	cEmAiAnalysis(cIO *Cp_IO,
// 		wsUint N_samp,
// 		wsReal **Ra_inpY, 
// 		wsUint N_inpPheno,
// 		wsReal **Ra_inpData,
// 		wsUint N_inpCol, wsUint N_idxSNP,
// 		wsReal **Ra_inpPhi,
// 		wsReal **Ra_inpPhiInv=NULL,
// 		char B_inpIsREML=1);
//	exportMatrix("test0.scr", Ra_onlyCov, N_sample, N_cov+1);

	/* Perform EM-AI */

	cStdMatrix	M_Xt(N_cov+1, N_testSamp, Ra_Xt, MATDEL_NONE);
//	cVector		V_Y(Ra_Y[0], N_testSamp, 1);
//	cStdMatrix	V_Y(N_pheno, N_testSamp, Ra_Y, 1);
	cSymMatrix	M_phi(Ra_phi, N_testSamp, 1);
	if (IS_ASSIGNED(est))
		LOG("Pre-assigned estimations loading from [%s]\n", OPT_STRING(est));
	else if (IS_ASSIGNED(blup))
		halt("BLUP loaded but no estimation was loaded, please assing --est with --blup");
	cTimer t;
	t.start();
	Cp_res0EmAi = new cEmAiAnalysisV2(Cp_IO, M_Xt, *Mp_Yt, M_phi, B_isREML, Ba_excl,
		Mp_EVX, Mp_EV);
 	Cp_res0EmAi->run();
	wsStrCst S = t.getReadable();

	if (Cp_res0EmAi->getIsFailed()) {
		if (!OPT_ENABLED(specdcmp)) {
			LOG("Null test failed using EM-AI algorithm has failed\n");
			LOG("Spectral decomposition will be used\n");

			delete Cp_res0EmAi;
			Cp_res0EmAi = new cEmAiAnalysisV2(Cp_IO, M_Xt, *Mp_Yt, M_phi,
				B_isREML, Ba_excl);
			Cp_res0EmAi->doSpecDcmp();
			Cp_res0EmAi->run();
		}

		if (Cp_res0EmAi->getIsFailed())
			halt_fmt(WISARD_FAIL_FIT_NULLMODEL);
//			halt("Failed to fitting NULL model");
	}
//	sseUnmat(Ra_X, N_testSamp);
	LOG("[%s] Null model fitted\n", S);

	/* Get the other estimates */
	if (!IS_ASSIGNED(est)) {
		Cp_res0EmAi->getEst();

		/* Print out estimates */
/**/	cExporter*	Cp_est			= cExporter::summon("poly.est.res");
		LOGoutput("Polygenic estimations are exported to [%s.poly.est.res]\n",
			OPT_STRING(out));
		wsReal		**Ra_est0		= Cp_res0EmAi->getEstimates();
		vPheno&		X_phenos		= Cp_IO->getPhenoInfo();
		const char	*Sa_headers[]	= { "sig2", "sig2g", "logL", "Var(sig2)",
			"Var(sig2g)", "h^2", "Var(h^2)", "Var(sig2+sig2g)" };

		/* Print header */
		Cp_est->put("ESTNAME");
		for (i=0 ; i<N_pheno ; i++)
			Cp_est->fmt("	%s", X_phenos[i].S_name.c_str());
		Cp_est->put("\n");

		Cp_est->put("ESTMETHOD");
		for (i=0 ; i<N_pheno ; i++)
			Cp_est->put(OPT_ENABLED(ml)?"	ML":"	REML");
		Cp_est->put("\n");


		/* Print estimates */
		for (i=0 ; i<EMAI_EST_SIZE ; i++) {
			Cp_est->fmt("%s", Sa_headers[i]);
			for (wsUint k=0 ; k<N_pheno ; k++)
				Cp_est->fmt("	%g", Ra_est0[k][i]);
				/* FIXME ::: */
			Cp_est->put("\n");
		}
		for (i=0 ; i<=N_cov ; i++) {
			Cp_est->fmt("beta_%d", i);
			for (wsUint k=0 ; k<N_pheno ; k++)
				Cp_est->fmt("	%g", Ra_est0[k][EMAI_EST_SIZE+i*2]);
			Cp_est->put("\n");
			Cp_est->fmt("Var(beta_%d)", i);
			for (wsUint k=0 ; k<N_pheno ; k++)
				Cp_est->fmt("	%g", Ra_est0[k][EMAI_EST_SIZE+i*2+1]);
			Cp_est->put("\n");			
		}
		delete Cp_est;
	} /* If --est is on, DO NOT print out estimations because it is all the same */

	deallocMatrix(Ra_Xt, N_cov+1, (void *)1);
	Ra_Xt = NULL;
}

#endif

} // End namespace ONETOOL
