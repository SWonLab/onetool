#include "finalstat.h"
#include "ppp.h"

/* !!! WISARD-specific analysis !!! */
#if TOOLSET_TYPE == TOOLSET_WISARD

cFinalStatAnalysisV2::cFinalStatAnalysisV2(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP, cAnalysis *Cp_inpAnaCorr) : cAnalysis(Cp_IO)
{
	UINT_t	i;
	UINT_t	N_SNP		= Cp_IO->getMarkerSize();
	UINT_t	N_sample	= Cp_IO->getSampleSize();
	REAL_c	*Ra_pheno	= Cp_IO->getPheno();

	MULTI_MALLOC(Ra_stat, REAL_t, N_SNP);
	MULTI_MALLOC(Ra_stdizeY, REAL_t, N_sample);

	Cp_anaCorr = Cp_inpAnaCorr;
	Cp_anaPPP = Cp_inpAnaPPP;

	/* Y-E(Y) */
	REAL_t eY = 0.0f;
	for (i=0 ; i<N_sample ; i++) /* Missing 고려 안 되어 있음 */
		eY += Ra_pheno[i];
	eY /= (REAL_t)N_sample;
	for (i=0 ; i<N_sample ; i++) /* Missing 고려 안 되어 있음 */
		Ra_stdizeY[i] = Ra_pheno[i] - eY;

	/* Export Y-E(Y) values if verbose */
	if (OPT_NUMBER(verbose)) {
/**/	cExporter* Cp_exporter = cExporter::summon("yey"); /* CHECKED */

		for (i=0 ; i<N_sample ; i++) /* Missing 고려 안 되어 있음 */
			Cp_exporter->fmt("SAMP%d	%g\n", i+1, Ra_stdizeY[i]);

		delete Cp_exporter;
	}
}

cFinalStatAnalysisV2::~cFinalStatAnalysisV2()
{
	LOG("~cFinalStatAnalysisV2() called\n");
	DEALLOC(Ra_stat);
	DEALLOC(Ra_stdizeY);
}

void cFinalStatAnalysisV2::run()
{
	UINT_t	N_SNP	= Cp_IO->getMarkerSize();

	/* Calculate final statistics */
	WORKER().run(calcFinalStat, forAllSNP_ana, this, Ra_stat);

	/* Export result */
/**/cExporter*	Cp_exporter = cExporter::summon("sta"); /* CHECKED */

	Cp_exporter->put("SNP	CHI^2STAT\n");
	for (UINT_t i=0 ; i<N_SNP ; i++) /* Missing 고려 안 되어 있음 */
		Cp_exporter->fmt("SNP%d	%g\n", i+1, Ra_stat[i]);
	delete Cp_exporter;
}

cPPPAnalysisV2* cFinalStatAnalysisV2::getPPPAnalysis() { return Cp_anaPPP; }

int cFinalStatAnalysisV2::calcFinalStat(int N_idx, void *Vp_shareData, void *Vp_data, void *VP_result)
{
	cFinalStatAnalysisV2
			*Cp_anaFS			= (cFinalStatAnalysisV2 *)Vp_shareData;
	cIO		*Cp_IO				= Cp_anaFS->getIO();
	UINT_t	N_sample			= Cp_IO->getSampleSize();
	cPPPAnalysisV2
			*Cp_anaPPP			= Cp_anaFS->getPPPAnalysis();
	REAL_t	**Ra_corr			= getFullCorMat(Cp_anaFS->getCorrAnalysis());
	REAL_t	**Ra_data			= Cp_IO->getData();
	REAL_c	*Ra_meanX			= Cp_anaPPP->getFounderMean();
	REAL_c	*Ra_ppp_sqrt2pq	= Cp_anaPPP->getPPPsq();
	REAL_t	*Ra_stdizeX		= NULL;
	UINT_t	N_idxSNP			= *((UINT_t *)Vp_data);

	MULTI_MALLOC(Ra_stdizeX, REAL_t, N_sample);

	/* X-E(X) */
	int N_availX = 0;
	for (UINT_t i=0 ; i<N_sample ; i++) {
		if (isAvailableReal(Ra_data[i][N_idxSNP])) {
			Ra_stdizeX[i] = Ra_data[i][N_idxSNP] - Ra_meanX[N_idxSNP];
			N_availX++;
		}
	}

	/* (Y-E(Y)) %*% cov(Ii,Ik) %*% (Y-E(Y))' = corr(Ii,Ik)*ppp(Sj) */
	REAL_t sum1 = 0.0f, sum2 = 0.0f;
	for (UINT_t i=0 ; i<N_sample ; i++) {
		if (isAvailableReal(Ra_data[i][N_idxSNP])) {
			sum1 += Ra_stdizeX[i] * Cp_anaFS->Ra_stdizeY[i];

			for (UINT_t k=0 ; k<N_sample ; k++) {
				if (isAvailableReal(Ra_data[i][k])) {
					sum2 += Cp_anaFS->Ra_stdizeY[i]
					*
						(Ra_corr[i][k] * Ra_ppp_sqrt2pq[N_idxSNP])
						*
						Cp_anaFS->Ra_stdizeY[k];
					//cov[i*szSamp + k] = cor[i][k] * ppp[j];
				}
			}
		}
	}
	Cp_anaFS->Ra_stat[N_idxSNP] = sum1 / sqrt(sum2);

	DEALLOC(Ra_stdizeX);

	return 0;
}

#endif
