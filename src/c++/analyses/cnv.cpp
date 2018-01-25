#include "analyses/cnv.h"
#include "analyses/emai.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cCnvAnalysis::cCnvAnalysis(cIO *Cp_inpIO, cAnalysis *Cp_inpAnaPhi)
	: cTestScoreAnalysis(Cp_inpIO, Cp_inpAnaPhi, NULL, !OPT_NUMBER(ml))
{
	N_cnvClus = 3;
}

cCnvAnalysis::~cCnvAnalysis()
{
}

void cCnvAnalysis::run()
{
//	wsReal R_varExplProp	= REAL_CONST(0.9);
//	wsReal R_cvgThre		= REAL_CONST(1e-5);
	wsReal **Ra_curCNV		= NULL;
	wsUint N_sample			= Cp_IO->sizeSample();
	wsReal **Ra_cnv			= Cp_IO->getCnvData();
	wsUint N_szCnv			= Cp_IO->getCnvSize();
	
	/* For each CNV, perform indep. */
	for (wsUint i=0 ; i<N_szCnv ; i++) {
		Ra_curCNV = sseMatrix(N_sample, 1);

		/* Make transposed ith CNV */
		for (wsUint j=0 ; j<N_sample ; j++)
			Ra_curCNV[j][0] = Ra_cnv[j][i];
		_doCnv(Ra_curCNV, 1);

		sseUnmat(Ra_curCNV, N_sample);
	}
}

void cCnvAnalysis::_doCnv(wsReal **Ra_CNV, wsUint N_szCnv)
{
	wsUint	i, j, k;
	xRealSort
			*Xa_sort	= NULL;
	wsUint	N_sample	= Cp_IO->sizeSample();
	wsUint	*Na_clus	= kmeans(Ra_CNV, N_szCnv, N_sample, N_cnvClus);


	/* Get the mean of each cluster */
	wsAlloc(Xa_sort, xRealSort, N_cnvClus);
	for (i=0 ; i<N_cnvClus ; i++) {
		wsReal R_mean = W0;

		for (j=0 ; j<N_sample ; j++) {
			if (Na_clus[j] != i) continue;
			for (k=0 ; k<N_szCnv ; k++)
				R_mean += Ra_CNV[k][j];
		}
		Xa_sort[i].i = i;
		Xa_sort[i].V = R_mean / (wsReal)N_sample;
	}
	qsort(Xa_sort, N_cnvClus, sizeof(xRealSort), sort_real);

	/* In this point, the value of Na_clus is irrelevant to the mean of
	 * each cluster. Thus, we need to arrange the value of Na_clus to
	 * increase as the mean of corresponding cluster grows */
	wsUint *Na_afterClus = NULL;
	wsAlloc(Na_afterClus, wsUint, N_cnvClus);
	for (i=0 ; i<N_cnvClus ; i++) {
		for (j=0 ; j<N_cnvClus ; j++)
			if (Xa_sort[j].i == i) {
				Na_afterClus[i] = j;
				break;
			}
	}
	DEALLOC(Xa_sort);
	/* In this point, Na_afterClus[i] -> Arranged index of cluster #i
	 * Thus, Na_afterClus[Na_clus[i]] refers the arranged index of cluster
	 * for sample (i) */
	for (i=0 ; i<N_sample ; i++)
		Na_clus[i] = Na_afterClus[Na_clus[i]];

	/* Account for the number of samples included in (i)th cluster,
	 * and (i) is proportional to the mean of cluster */
	wsUint *Na_szClus = NULL;
	wsAlloc(Na_szClus, wsUint, N_cnvClus);

	/* Rearrange samples into clusters */
	wsReal **Ra_ordCnv = NULL;
	wsUint N_idxOrdCnv = 0;
	wsAlloc(Ra_ordCnv, wsReal*, N_sample);
	for (i=0 ; i<N_cnvClus ; i++) {
		wsUint N_szClus = 0;
		/* Find samples having its cluster == i */
		for (j=0 ; j<N_sample ; j++)
			if (Na_clus[j] == i) {
				Ra_ordCnv[N_idxOrdCnv] = Ra_CNV[j];
				N_idxOrdCnv++;
				N_szClus++;
			}
		Na_szClus[i] = N_szClus;
	}
	if (N_idxOrdCnv != N_sample)
		halt("System rror : Re-ordering process have bug");

	/* Ra_ordCnv now have below form :
	 * 
	 * #idx | #cnv1 ... #cnvN
	 * -----+----------------
	 *  ... |
	 * clus0| (minimum mean)
	 *  ... |
	 * -----+----------------
	 *  ... |
	 * clus1| (second-minimum mean...)
	 *  ... */

	/* Get initial value */
	wsReal *Ra_mu		= NULL;
	wsReal ***Ra_cov	= NULL;
	wsReal *Ra_detCov	= NULL;
	wsReal *Ra_lDetCov	= NULL;
	wsReal ***Ra_covInv	= NULL;
	wsAlloc(Ra_mu, wsReal, N_cnvClus);
	wsAlloc(Ra_cov, wsReal**, N_cnvClus);
	wsAlloc(Ra_detCov, wsReal, N_cnvClus);
	wsAlloc(Ra_lDetCov, wsReal, N_cnvClus);
	wsAlloc(Ra_covInv, wsReal**, N_cnvClus);
	for (i=j=0 ; i<N_cnvClus ; i++) {
		Ra_mu[Na_afterClus[i]] = Xa_sort[i].V;
		Ra_cov[i] = sseMcovCol(Ra_ordCnv+j, Na_szClus[i], N_szCnv);
		Ra_detCov[i] = detCholesky(Ra_cov[i], N_szCnv);
		if (Ra_detCov[i] != Ra_detCov[i])
			halt("Failed to get determinant of [%d]th CNV cluster", i);
		Ra_lDetCov[i] = log(Ra_detCov[i]);
		Ra_covInv[i] = SO_invMatrix(Ra_cov[i], N_szCnv);

		j += Na_szClus[i];
	}

	//wsUint N_iter = 0;
// 	while (1) {
// 		wsReal **Ra_pin = sseMatrix(N_sample, 3);
// 
// 		/* V <- sig2*diag(3)+sig2g*phi1 */
// 		for (i=0 ; i<N_sample ; i++) {
// 			for (j=0 ; j<N_cnvClus ; j++){}
// 				/* pin[i,n]  <- Ddel[n]^(-1/2) * exp(-1/2*(pX[i,]-mu[[n]])
// 				 *		%*% Invdel[[n]] %*% as.matrix(pX[i,]-mu[[n]])) * P[n]*/
// 		}
// 
// 		N_iter++;
// 	}

	sseFree(Na_clus);
}

#endif

} // End namespace ONETOOL
