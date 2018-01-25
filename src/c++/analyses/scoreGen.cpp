#include "analyses/emai.h"
#include "analyses/scoreGen.h"
#include "utils/matrix.h"
#include "utils/stat.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cTestGenScoreAnalysis::cTestGenScoreAnalysis(cIO *Cp_IO, wsUintCst N_inpSamp,
	cVector	&V_YtP,
//	wsReal *Ra_inpPheno,
	wsReal **Ra_inpPhi, wsUint *Na_inpSzFam, cEmAiAnalysisV2 *Cp_inpAnaEmAi,
	cStdMatrix *Mp_inpWtP, wsUintCst N_inpColX,
	cVector	&V_inpRtP,
//	wsReal *Ra_inpResidual,
	cDiagMatrix &M_inpInvV, cStdMatrix *Mp_inpEVX, cVector *Vp_inpXtP) : cAnalysis(Cp_IO)
{
	/* Copy SNP indices */
	Vp_XtP		= Vp_inpXtP;
//	Ra_snp		= Ra_inpSnp;
	N_samp		= N_inpSamp;
	Na_szFam	= Na_inpSzFam;
	Vp_YtP		= &V_YtP;
//	Ra_pheno	= Ra_inpPheno;
	Ra_phi		= Ra_inpPhi;
//	Ra_Xt		= Ra_inpXt;
	Mp_WtP		= Mp_inpWtP;
	Ra_colX		= N_inpColX;
	Cp_anaEmAi	= Cp_inpAnaEmAi;
//	Ra_invV		= Ra_inpInvV;
	Mp_invV		= &M_inpInvV;
	Mp_EVX		= Mp_inpEVX;
//	Ra_residual	= Ra_inpResidual;
	Vp_RtP		= &V_inpRtP;
}

cTestGenScoreAnalysis::~cTestGenScoreAnalysis()
{
}

/* X <- covariate matrix
 * Y <- phenotype
 * sig2g, sig2 <- estimate
 */
void cTestGenScoreAnalysis::run()
{
	wsUint	i, j;
	mFam	&Xa_fam		= Cp_IO->getFamilyData();
	wsUint	N_colX		= Ra_colX;

	cVector	V_ntVx(N_colX);
	//wsReal	**Ra_ntVx	= sseEmptyMatrix(1, N_colX);		 /* 120917 */
	cStdMatrix
			M_xtVx(N_colX, N_colX);
//	wsReal	**Ra_xtVx	= sseEmptyMatrix(N_colX, N_colX);	 /* 120917 */
	// 		d11 <- 0
	// 		s1 <- 0
	// 		d12tran <- 0
	// 		d22 <- 0
	wsReal	R_d11		= W0;					 /* 120917 */
	wsReal	R_s1		= W0;					 /* 120917 */
	wsReal	*Ra_d12t	= NULL;								 /* 120917 */
	sseCalloc(Ra_d12t, wsReal, N_colX);
	cStdMatrix	M_d22(N_colX, N_colX);
//	wsReal	**Ra_d22	= sseEmptyMatrix(N_colX, N_colX);	 /* 120917 */
	wsMat		Ra_d22	= M_d22.get();

	/* Which is biggest */
	wsUint Nsz = Na_szFam[0];
	for (i=1 ; i<Xa_fam.size() ; i++)
		if (Nsz < Na_szFam[i]) Nsz = Na_szFam[i];
// 	wsReal	*Ra_snpSub	= NULL;
// 	sseMalloc(Ra_snpSub, wsReal, Nsz);

	// [1]	F <- length(fam_num)
	// for(i in 1:F){
	wsUint I=0, J;
	for (wsUint iii=0 ; iii<Xa_fam.size() ; iii++) {
		wsUint N_sz = Na_szFam[iii];

		// phi_list <- subset of correlation structure
		// #mem_i	<- # of members of ith family
		cStdMatrix	M_famXt(N_colX, N_sz);
		wsReal	**Ra_famXt	= M_famXt.get();

		// [#mem_i*#snp]
		// fam_n[[i]] <- snp[(cut[i]+1):cut[i+1]]
		// Equivalent to Ra_snp

		// [#mem_i*#inp]	fam_x[[i]] <- x[(cut[i]+1):cut[i+1],]
		J = I;
		for (i=0 ; i<N_colX ; i++)
			memcpy(Ra_famXt[i], Mp_WtP->get()[i]+J, sizeof(wsReal)*N_sz);
		I += N_sz;
		// [#mem_i*#mem_i]	cov <- sig2g*phi_list[[i]]+sig2*diag(fam_num[i])
		// [#mem_i*#mem_i]	Invcov <- solve(cov)
		cDiagMatrix	M_invCov	= Mp_invV->subset(J, N_sz);
//		wsReal **Ra_invCov = sseSubsetMatrix(Ra_invV, J, N_sz, J, N_sz);

		cVector		V_snpSub	= Vp_XtP->subset(J, N_sz);
//		memcpy(Ra_snpSub, Ra_snp+J, sizeof(wsReal)*N_sz);

		// [#snp*#mem_i]%*%[#mem_i*#mem_i]
		// [#snp*#mem_i]	ntInvcov <- t(fam_n[[i]])%*%Invcov
		cVector		V_ntInvcov	= V_snpSub * M_invCov;
//		wsReal **Ra_ntInvcov = sseMpM(&Ra_snpSub, 1, N_sz,
//			Ra_invCov, N_sz, N_sz); /* 120917 */

		// [#snp*#mem_i]%*%[#mem_i*#pheno]
		// [#snp*#pheno]	nve[[i]] <- ntInvcov%*%fam_e[[i]]
//		wsReal R_nvE = multVV(Ra_ntInvcov[0], N_sz, Ra_residual+J);
		wsReal R_nvE = multVV(V_ntInvcov.get(), N_sz, Vp_RtP->get()+J);

		// [#snp*#mem_i]%*%[#mem_i*#inp]
		// [#snp*#inp]		cur_nvx <- ntInvcov%*%fam_x[[i]]
		cVector	V_curNtVx	= V_ntInvcov * M_famXt;
// 		wsReal **Ra_curNtVx = sseMpMt(V_ntInvcov.get(), 1, N_sz,
// 			Ra_famXt, N_colX, N_sz); /* 120917 */
// 		sseUnmat(Ra_ntInvcov, 1);

		// [#snp*#inp]		nvx <- nvx + cur_nvx
		V_ntVx += V_curNtVx;
// 		for (j=0 ; j<N_colX ; j++)
// 			Ra_ntVx[0][j] += Ra_curNtVx[0][j];
// 		sseUnmat(Ra_curNtVx, 1);

		// [#inp*#mem_i]%*%[#mem_i*#mem_i]
		// [#inp*#mem_i]	xtInvcov <- t(fam_x[[i]])%*%Invcov
		cStdMatrix	M_xtInvcov = M_famXt * M_invCov;
// 		wsReal **Ra_xtInvcov = sseMpM(Ra_famXt, N_colX, N_sz,
// 			Ra_invCov, N_sz, N_sz); /* 120917 */
//		sseUnmat(Ra_invCov, N_sz);

		// [#inp*#mem_i]%*%[#mem_i*#inp]
		// [#inp*#inp]		cur_xvx <- xtInvcov%*%fam_x[[i]]
		cStdMatrix	M_curXtVx	= M_xtInvcov.Mt(M_famXt);
// 		wsReal **Ra_curXtVx = sseMpMt(Ra_xtInvcov, N_colX, N_sz,
// 			Ra_famXt, N_colX, N_sz); /* 120917 */
// 		sseUnmat(Ra_famXt, N_colX);

		// [#inp*#inp]		xvx <- xvx + cur_xvx
		M_curXtVx += M_curXtVx;
// 		for (i=0 ; i<N_colX ; i++)
// 			for (j=0 ; j<N_colX ; j++)
// 				Ra_xtVx[i][j] += Ra_curXtVx[i][j];
// 		sseUnmat(Ra_curXtVx, N_colX);

		// [#inp*#mem_i]%*%[#mem_i*#pheno]
		// [#inp*#pheno] xve[[i]] <- xtInvcov%*%fam_e[[i]]
		wsReal *Ra_xvE = multMV(M_xtInvcov.get(), N_colX, N_sz, Vp_RtP->get()+J); /* 120917 */
//		sseUnmat(Ra_xtInvcov, N_colX);

		// + [#snp*#pheno]%*%[#pheno*#snp]
		// [#snp*#snp]	d11 <- d11+nve[[i]]%*%t(nve[[i]])
		R_d11	+= SQR(R_nvE);

		// [#snp*#snp]	s1 <- s1+nve[[i]]
		R_s1	+= R_nvE;

		// FIXME : [#inp*#pheno]%*%[#snp*#pheno] in code of Meiling
		// Maybe should be t(nve[[i]])
		// 
		// [#inp*#pheno]%*%[#snp*#pheno]
		// [#inp*#snp]	d12tran <- d12tran+xve[[i]]%*%nve[[i]]
		for (i=0 ; i<N_colX ; i++)
			Ra_d12t[i] += Ra_xvE[i] * R_nvE;

		// [#inp*#pheno]%*%[#pheno*#inp]
		// [#inp*#inp]	d22 <- d22+xve[[i]]%*%t(xve[[i]])
		for (i=0 ; i<N_colX ; i++)
			for (j=0 ; j<N_colX ; j++)
				Ra_d22[i][j] += Ra_xvE[i] * Ra_xvE[j];
		DEALLOC(Ra_xvE);
	}
	if (I != N_samp)
		halt("SYSERR: Sum of available family member expected %d, counted %d",
			N_samp, I);
//	sseFree(Ra_snpSub);
#ifdef EMAIvalidate
	exportMatrix("ea.gen.xtVx", Ra_xtVx, N_colX, N_colX);//EMAIvalidate
	exportMatrix("ea.gen.nvx", Ra_ntVx, 1, N_colX);//EMAIvalidate
	exportMatrix("ea.gen.d22", Ra_d22, N_colX, N_colX);//EMAIvalidate
#endif

	// [#inp*#inp]		Invxvx	<- solve(xvx)
	cStdMatrix &M_xtVx_i = M_xtVx.inv();
//	wsReal **Ra_xtVx_i = dbgSO_invMatrix(Ra_xtVx, N_colX); /* 120917 */
// 	if (Ra_xtVx_i == NULL)
// 		halt("Can't get xtVx inverse matrix");
//	printf("xtVx_i[0] %x\nxtVx_i[1] %x\n", Ra_xtVx_i[0], Ra_xtVx_i[1]);
//	sseUnmat(Ra_xtVx, N_colX);

	// [#snp*#inp]%*%[#inp*#inp]
	// [#snp*#inp]		nvxInvxvx <- nvx%*%Invxvx
	cVector V_ntVx_xtVx_i = V_ntVx * M_xtVx_i;
// 	wsReal **Ra_ntVx_xtVx_i
// 		= multMM(Ra_ntVx, 1, N_colX, Ra_xtVx_i, N_colX, N_colX); /* 120917 */

	// [#snp*#inp]%*%[#inp*#inp]%*%[#inp*#snp]
	// [#snp*#snp]		IID		<- nvxInvxvx%*%d12tran
	wsReal R_IID = V_ntVx_xtVx_i.sum(Ra_d12t, N_colX);
//	wsReal R_IID = multVV(Ra_ntVx_xtVx_i[0], N_colX, Ra_d12t);
	sseFree(Ra_d12t);

	// [#snp*#inp] %*% [#inp*#inp] %*% [#inp*#inp] %*% [#inp*#snp]
	// In case of #snp=1, nvxInvxvx %*% d22
	// and multVMV(nvxInvxvx %*% d22, Invxvx, t(nvx))

	// [#snp*#inp] %*% [#inp*#inp]
	// [#snp*#inp]		nvxInvxvxd22 <- nvxInvxvx %*% d22
// 	wsReal **Ra_nvxInvxvx_d22 = multMM(Ra_ntVx_xtVx_i, 1, N_colX,
// 		Ra_d22, N_colX, N_colX); /* 120917 */
	cVector	V_nvxInvxvx_d22 = V_ntVx_xtVx_i * M_d22;
// 	wsReal **Ra_nvxInvxvx_d22 = multMM(Ra_ntVx_xtVx_i, 1, N_colX,
// 		Ra_d22, N_colX, N_colX); /* 120917 */
// 	deallocMatrix(Ra_ntVx_xtVx_i, 1);
// 	sseUnmat(Ra_d22, N_colX);

	// [#snp*#inp] %*% [#inp*#inp] %*% [#inp*#snp]
	// [#snp*#snp]		IIDII <- nvxInvxvxd22 %*% Invxvx %*% nvx_t
	wsReal R_IIDII = V_nvxInvxvx_d22.qf(M_xtVx_i, V_ntVx);
// 	wsReal R_IIDII = multVMV(Ra_nvxInvxvx_d22[0], Ra_xtVx_i, N_colX,
// 		Ra_ntVx[0]);
// 	sseUnmat(Ra_xtVx_i, N_colX);
// 	deallocMatrix(Ra_nvxInvxvx_d22, 1);
// 	sseUnmat(Ra_ntVx, 1);

	// [#snp*#snp] - [#snp*#snp] - [#snp*#snp] + [#snp*#snp]
	// In case of #snp=1, sx <- d11-2*IID+IIDII cause IID=t(IID)
	// [#snp*#snp]		sx <- d11-IID-t(IID)+IIDII
	R_d11 -= R_IID + R_IID - R_IIDII;

	// [#snp*#snp] %*% [#snp*#snp] %*% [#snp*#snp]
	// [#snp*#snp]		Tgs <- t(s1) %*% solve(sx) %*% s1
	R_Tgs = R_s1 * (W1/R_d11) * R_s1;

	// FIXME : How to calculate the p-value under chisq dist. with Tgs?
	// In case of #snp>1, Tgs is MATRIX with [#snp*#snp] dim.
	// Spvalue <- 1-pchisq(Tgs, #snp)
	R_pValue = (wsReal)PVchisq(R_Tgs, 1.0);
}

#endif

} // End namespace ONETOOL
