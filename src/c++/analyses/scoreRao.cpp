#include "analyses/emai.h"
#include "analyses/scoreRao.h"
#include "utils/matrix.h"
#include "utils/stat.h"

namespace ONETOOL {

//wsReal** cTestRaoAnalysis::Ra_XVIXinv = NULL;
cSymMatrix cTestRaoAnalysis::M_WViWi;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cTestRaoAnalysis::cTestRaoAnalysis(cIO *Cp_IO, wsUint N_inpSamp,
	cVector &V_inpYtP, wsReal **Ra_inpPhi, cEmAiAnalysisV2 *Cp_inpAnaEmAi,
	cStdMatrix *Mp_inpWtP,
	wsUint N_inpColX, cVector &V_inpRtP,
	cDiagMatrix &M_inpInvV, cStdMatrix *Mp_inpEVX, cVector *Vp_inpXtP) : cAnalysis(Cp_IO)
{
//	Ra_snp		= Ra_inpSnp;
	Vp_XtP		= Vp_inpXtP;
	Cp_anaEmAi	= Cp_inpAnaEmAi;
	N_samp		= N_inpSamp;
//	Ra_residual	= Ra_inpResidual;
	Vp_RtP		= &V_inpRtP;
//	Ra_pheno	= Ra_inpPheno;
	Vp_YtP		= &V_inpYtP;
	Ra_phi		= Ra_inpPhi;
//	Ra_Xt		= Ra_inpXt;	/* [#inp*#samp] */
	Mp_WtP		= Mp_inpWtP;
	N_colX		= N_inpColX;
	Mp_Vinv		= &M_inpInvV;
	Mp_EVX		= Mp_inpEVX;
}

cTestRaoAnalysis::~cTestRaoAnalysis()
{

}

/* X <- covariate matrix
 * sig2g, sig2 <- estimate
 */
void cTestRaoAnalysis::run()
{
//	wsUint	i;

	// [#snp*#samp]%*%[#samp*#samp]
	// [#snp*#samp]		snpInvcov <- t(snp)%*%Invcov
//	wsReal *Ra_snpVinv = sseMpV(&Ra_snp, 1, N_samp, Mp_Vinv->get()[0]);
	cVector	V_XVi	= *Vp_XtP * *Mp_Vinv;
//	wsReal *Ra_snpVinv = sseMpV(&Ra_snp, 1, N_samp, Mp_Vinv->get()[0]);
#ifdef EMAIvalidate
	exportMatrix("em.rao.snpInvcov", Ra_snpVinv, 1, N_samp);//EMAIvalidate
#endif

	// [#snp*#samp]%*%[#samp*#pheno]
	// [#snp*#pheno]	nve <- snpInvcov%*%e 
	wsReal	R_xVe	= V_XVi.sum(*Vp_RtP);
//	wsReal	R_nvE	= sseVV(Ra_snpVinv, N_samp, Ra_residual);

	// [#snp*#samp]%*%[#samp*#snp]
	// [#snp*#snp]		nvn <- snpInvcov%*%snp
	wsReal	R_xVx	= V_XVi.sum(*Vp_XtP);
// 	wsReal R_nvN = W0;
// 	for (i=0 ; i<N_samp ; i++)
// 		R_nvN += Ra_snpVinv[i] * Ra_snp[i];

	// [#snp*#samp]%*%[#samp*#inp]
	// [#snp*#inp]		nvx <- snpInvcov%*%x
	cVector	V_xVW	= V_XVi * *Mp_WtP;
// 	wsReal **Ra_nvX = multMMt(&Ra_snpVinv, 1, N_samp,
// 		Ra_Xt, N_colX, N_samp);
// 	sseFree(Ra_snpVinv);
// 	sseUnmat(Ra_snpVinv, 1);

	// [#inp*#samp]%*%[#samp*#samp]%*%[#samp*#inp]
	// [#inp*#inp]		xvx <- t(x)%¤·*%Invcov%*%x
	// Calculated from cTestScoreAnalysis
// 	if (cTestRaoAnalysis::Ra_XVIXinv == NULL) {
// 		LOG("Calculating XVIXinv\n");
// 		/* This term requires only one-time calculation
// 		 * Because all X and Vi will be equal on all tests */
// 		wsReal **Ra_XVX = sseMMMt(Ra_Xt, N_colX, N_samp,
// 			Ra_Vinv, N_samp, N_samp,
// 			Ra_Xt, N_colX, N_samp);
// 		cTestRaoAnalysis::Ra_XVIXinv = SO_invMatrix(Ra_XVX, N_colX);
// 		sseUnmat(Ra_XVX, N_colX);
// 	}

#ifdef EMAIvalidate
	exportMatrix("em.rao.nvX", Ra_nvX, 1, N_colX);//EMAIvalidate
	//exportMatrix("em.rao.XVX", Ra_XVX, N_colX, N_colX);//EMAIvalidate
#endif

	// FIXME : Very weird, it cannot be calculated...
	// 
	// [#snp*#pheno] %*%
	//		([#snp*#snp] - [#snp*#inp]%*%[#inp*#inp]%*%[#inp*#snp]) %*%
	// [#pheno*#snp]
	// 
	// [?*?] Trs <- nve%*%solve(nvn-nvx%*%solve(xvx)%*%t(nvx))%*%t(nve)
	// Anyway, currently it is 1*1
// 	R_Trs = R_nvE *
// 		(W1 / (R_nvN - multVMV(Ra_nvX[0], cTestRaoAnalysis::Ra_XVIXinv,
// 			N_colX))) *
// 		R_nvE;
	R_Trs = SQR(R_xVe) / (R_xVx - V_xVW.qf(cTestRaoAnalysis::M_WViWi));

//	deallocMatrix(Ra_nvX, 1);
	// RSpvalue <- 1-pchisq(Trs,1)
	// FIXME : What governs DF?
	R_pValue = (wsReal)PVchisq(R_Trs, 1.0);
// 		detach(dat)
// 		return(list(RScoreT=Trs,RSpvalue=RSpvalue))
// 		
}

#endif

} // End namespace ONETOOL
