#include <string.h>
#include <stdlib.h>
#include <map>
#include <iterator>
#include "analyses/pca.h"
#include "analyses/fam.h"
#include "analyses/ppp.h"
#include "analyses/corr.h"
#include "global/Rconn.h"
#include "utils/matrix.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

#define SIGN(Ra_mat, b) ( (b) < 0 ? -fabs(Ra_mat) : fabs(Ra_mat) )

int forAllSample_pca(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData)
{
	static wsUint N_idxSample;
	int		*Np_idx = (int *)Vp_thrData;
	wsUint N_pcaSample = ((cPcaAnalysisV2 *)Vp_shareData)->getPcaSampleSize();
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		N_idxSample = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++) {
			if ((N_idxSample+i) >= N_pcaSample)
				break;
			Np_idx[i] = N_idxSample+i;
			if (((N_idxSample+i)%10) == 0)
				notice("%d/%d samples processed\r", N_idxSample+i, N_pcaSample);
		}
		N_idxSample += N_thread;
		break;
	case DISTFUN_UNINIT:
		LOG("%d/%d samples processed\n", N_pcaSample, N_pcaSample);
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}
	return i;
}

cPcaAnalysisV2::cPcaAnalysisV2(cIO *Cp_IO, cPPPAnalysisV2 *Cp_inpAnaPPP,
	cFamStrAnalysis *Cp_inpAnaFamStr) : cAnalysis(Cp_IO)
{
	wsUint		N_SNP			= Cp_IO->sizeVariant();	///< # of SNPs in dataset
	wsUint		N_founder		= Cp_IO->sizeFounder();	///< # of founders in dataset
//	MYuint		N_sample		= Cp_IO->getSampleSize();	///< # of samples in dataset
//	const char	*Ba_isFounder	= Cp_IO->getIsFounder();	///< Array[size #sample], 1 if (i)th sample is founder, otherwise 0
	wsFloat**	Ra_data			= NULL;//Cp_IO->getData();	///< Matrix[size #sample*#SNP], Normalized by PPP in this point
	vSampPtr	Xa_sample		= Cp_IO->getSample();	///< Array[size #sample], having the description about each sample
	wsUint		i, j;
	wsUint		N_finalFounder	= N_founder;
	mSampPtr	&Xm_mfData		= Cp_inpAnaFamStr->getMissingFounderData();
	mNucFam		&Xm_nfData		= Cp_inpAnaFamStr->getNucFamData();
	wsReal		**Ra_misFounderGeno
								= Cp_inpAnaPPP->getMisFounderGeno();

	Ra_cov = NULL;
	wsUint	N_vtSample = 0, N_unavailVtSamp = 0;
	Cp_anaFamStr	= Cp_inpAnaFamStr;
	Cp_anaPPP		= Cp_inpAnaPPP;
	/* if usemf == 1, N_finalFounder = N_founder + N_psuFounder */
	if (OPT_ENABLED(usemf)) {
		N_finalFounder += (wsUint)Xm_mfData.size();
		N_vtSample = (wsUint)Xm_mfData.size();
		LOG("%d missing parents by --usemf added to PC analysis due to missing parents\n", N_vtSample);
	} else if (!OPT_ENABLED(founderonly)) {
		/* �ٰ��� �� �θ� founder�� ��� virtual sample�� �� �� �ִ� */
		for (mNucFam_it it=Xm_nfData.begin() ; it!=Xm_nfData.end() ; it++) {
			if ((it->second.Xp_mat->Xp_mat == NULL && it->second.Xp_mat->Xp_pat == NULL) ||
				(it->second.Xp_pat->Xp_mat == NULL && it->second.Xp_pat->Xp_pat == NULL)) {

				int B_haveNonmissingChild = 0;
				FOREACH (vSampPtr_it, it->second.Xp_childs, cit)
					if ((*cit)->B_isMissing == 0) {
						B_haveNonmissingChild = 1;
						break;
					}

				if (B_haveNonmissingChild)
					N_vtSample++;
				else
					N_unavailVtSamp++;
			}
		} 
		N_finalFounder += N_vtSample;

		LOG("%d virtual samples added to PC analysis due to missing parents\n", N_vtSample);
		if (N_unavailVtSamp)
			LOG("%d unusable virtual samples found\n", N_unavailVtSamp);
	} else {
		/* When using --founderonly */
	}
	printf("%d final founders found\n", N_finalFounder);


	/* Make founder-only dataset */
	vSampPtr Xa_samples = Cp_IO->getSample();
	j = 0;
	if (!IS_ASSIGNED(cor)) {
		Ra_pcaInputData = NULL;
		wsAlloc(Ra_pcaInputData, wsFloat*, N_finalFounder);
		for (i=0 ; i<N_finalFounder ; i++) {
			sseCalloc(Ra_pcaInputData[i], wsFloat, N_SNP);
		}
		Ra_data		= Cp_IO->getData();			///< Matrix[size #sample*#SNP], Normalized by PPP in this point

		FOREACH(vSampPtr_it, Xa_samples, it) {
			xSample	*Xp_currSample = *it;
			/* It is founder */
			if (Xp_currSample->Xp_mat == NULL && Xp_currSample->Xp_pat == NULL) {
				memcpy(Ra_pcaInputData[j], Ra_data[Xp_currSample->N_idx], sizeof(wsFloat)*N_SNP);
				X_iid2fouIdx.insert(make_pair(Xp_currSample->S_IID, j));

				Sa_fidFounder.push_back(Xp_currSample->S_FID);
				Sa_iidFounder.push_back(Xp_currSample->S_IID);
				j++;
			}
		}
		if (N_founder != j)
			halt("N_founder[%d] in cIO and re-counted number of parents[%d] is not match, maybe system bug?", N_founder, j);
	} else if (0) {
		/* Otherwise use virtual samples only */
		//Ra_pcaInputData = sseMatrix(N_vtSample, N_SNP);
		Ra_pcaInputData = NULL;
		wsAlloc(Ra_pcaInputData, wsFloat*, N_vtSample);
		for (i=0 ; i<N_vtSample ; i++)
			sseCalloc(Ra_pcaInputData[i], wsFloat, N_SNP);	} else {
		Ra_pcaInputData = NULL;
		FOREACH(vSampPtr_it, Xa_samples, it) {
			xSample	*Xp_currSample = *it;
			/* It is founder */
			if (Xp_currSample->Xp_mat == NULL && Xp_currSample->Xp_pat == NULL) {
				X_iid2fouIdx.insert(make_pair(Xp_currSample->S_IID, j));

				Sa_fidFounder.push_back(Xp_currSample->S_FID);
				Sa_iidFounder.push_back(Xp_currSample->S_IID);
			}
		}
	}

	/* Copy psuedo-founder's data if usemf == 1 */
	wsUint	N_missingGeno = 0;
	if (OPT_ENABLED(usemf)) {
		mSampPtr_it it=Xm_mfData.begin();
		for (i=0 ; i<Xm_mfData.size() ; i++,it++) {
			X_iid2fouIdx.insert(make_pair(it->first, j));
#ifdef USE_CORRv2
			LOOP(k, N_SNP)
				Ra_pcaInputData[j][k] = (wsFloat)Ra_misFounderGeno[i][k];
			j++;
#else
			memcpy(Ra_pcaInputData[j++], Ra_misFounderGeno[i], sizeof(wsFloat)*N_SNP);
#endif

			Sa_fidFounder.push_back("MF");
			Sa_iidFounder.push_back(it->first);
		}
	} else if (!OPT_ENABLED(founderonly)) {
		/* Missing founder�� �����ϴ� �ٰ����� �ڽ� ������ �̿��� virtual sample�� �����Ѵ� */
		for (mNucFam_it it=Xm_nfData.begin() ; it!=Xm_nfData.end() ; it++) {
			/* If one of parents in this nuclear family is not founder,
			 * do not consider this family */
			if ((it->second.Xp_mat->Xp_mat || it->second.Xp_mat->Xp_pat) &&
				(it->second.Xp_pat->Xp_mat || it->second.Xp_pat->Xp_pat))
				continue;

			/* Due to this case considers both parents are founder and missing,
			 * each nuclear family will considered as ONE individual set */
			vector<int> Na_childDataIdx;
			Sa_fidFounder.push_back("MF");
			Sa_iidFounder.push_back(it->second.Xp_pat->S_IID+","+it->second.Xp_mat->S_IID);

			/* Insert all non-missing childs included in this nuclear family */
// 			LOG("Nuclear family %s,%s :\n",
// 				it->second.Xp_pat->S_IID.c_str(), it->second.Xp_mat->S_IID.c_str());
			FOREACH (vSampPtr_it, it->second.Xp_childs, cit) {
				if ((*cit)->B_isMissing == 0) {
					X_iid2fouIdx.insert(make_pair((*cit)->S_IID, j));

					/* Insert this child's data index into the vecotr of
					 * nuclear family what this child belonged */
					Na_childDataIdx.push_back((*cit)->N_idx);
//					LOG("	Child idx %d\n", (*cit)->N_idx);
				}
			}

			/* If there is no one who non-missing its data,
			 * skip this nuclear family, they will excluded */
			if (Na_childDataIdx.size() == 0) {
				LOG("Nuclear family from [%s,%s] from family [%s] does not have nonmissing child\n",
					it->second.Xp_pat->S_IID.c_str(), it->second.Xp_mat->S_IID.c_str(), it->second.Xp_mat->S_FID.c_str());
				continue;
			}

			/* Aggregate non-missing child's genotype into
			 * one representative */
			wsUint N_sz = (wsUint)Na_childDataIdx.size();
			if (Cp_IO->isDataComplete()) {
				for (wsUint k=0 ; k<N_SNP ; k++) {
					wsReal R_sum = W0;
					for (i=0 ; i<N_sz ; i++)
						R_sum += Ra_data[Na_childDataIdx[i]][k];
					Ra_pcaInputData[j][k] = (wsFloat)(R_sum / N_sz);
				}
			} else {
				for (wsUint k=0 ; k<N_SNP ; k++) {
					wsUint N = 0;
					wsReal R_sum = W0;
					for (i=0 ; i<N_sz ; i++) {
						if (isAvailableReal(Ra_data[Na_childDataIdx[i]][k])) {
							R_sum += Ra_data[Na_childDataIdx[i]][k];
							N++;
						}
					}
					if (N)
						Ra_pcaInputData[j][k] = (wsFloat)(R_sum / (wsReal)N);
					else {
						Ra_pcaInputData[j][k] = (wsFloat)WISARD_NA_REAL;
						N_missingGeno++;
					}
				}
			}
			j++;
		}
	}
	B_isPcaDataComplete = N_missingGeno == 0;
	/* Sanity check */
	if ((IS_ASSIGNED(cor) && N_vtSample != j) ||
		(!IS_ASSIGNED(cor) && N_finalFounder != j))
		halt("N_finalFounder[%d] in cIO and re-counted number of founder-array-included[%d] is not match, maybe system bug?", N_finalFounder, j);

	/* Export PCA input dataset if verbose */
// 	if (OPT_NUMBER(verbose)) {
// 		cExporter* Cp_exporter = cExporter::summon("pds"); /* CHECKED */
// 		LOG("Export standardized PCA input dataset by --verbose, virtual sample only\n");
// 
// 		FILE *fp = Cp_exporter->getHandle();
// 		for (i=0 ; i<N_finalFounder ; i++) {
// 			for (j=0 ; j<N_SNP ; j++)
// 				fprintf(fp, "%g ", Ra_pcaInputData[i][j]);
// 			fprintf(fp, "\n");
// 
// 			if ((i%10) == 9)
// 				notice("%d samples were exported\r", i);
// 		}
// 		LOG("%d[%d+%d] samples were exported\n", i, N_founder, N_vtSample);
// 	}

	/* Create correlation-related matrix */
	if (!IS_ASSIGNED(cor)) {
		wsAlloc(Ra_corr , wsReal*, N_finalFounder);
		wsAlloc(Na_corrN, int*  , N_finalFounder);
		for (wsUint i=0 ; i<N_finalFounder ; i++) {
			sseMalloc(Ra_corr[i] , wsReal, N_finalFounder);
			wsAlloc(Na_corrN[i], int  , N_finalFounder);
		}
		/* Allocate extra memory if data are complete */
		if (B_isPcaDataComplete) {
			wsAlloc(Ra_ysum , wsReal, N_finalFounder);
			wsAlloc(Ra_y2sum, wsReal, N_finalFounder);

			/* Calculate ysum and y2sum if data are complete */
			for (wsUint i=0 ; i<N_finalFounder ; i++) {
				wsFloat	*Rp_row = Ra_pcaInputData[i];

				Ra_ysum[i]	= W0;
				Ra_y2sum[i]	= W0;

				for (wsUint j=0 ; j<N_SNP ; j++) {
					Ra_ysum[i]	+= Rp_row[j];
					Ra_y2sum[i]	+= Rp_row[j]*Rp_row[j];
				}
			}
		} else {
			/* Ra_ysum and Ra_y2sum should be marked as NULL if the dataset is NOT complete */
			Ra_ysum	 = NULL;
			Ra_y2sum = NULL;
		}
	} else { 
		LOG("Use previously computed correlation 'as is'\n");
		cCorrAnalysis *Cp_cor = (cCorrAnalysis *)getCorrelation(getIO());
		Ra_corr		= Cp_cor->getCorr();
		Na_corrN	= NULL;
		Ra_ysum		= NULL;
		Ra_y2sum	= NULL;
	}

	N_pcaSample = N_finalFounder;
}

cPcaAnalysisV2::~cPcaAnalysisV2()
{
	/* Deallocate memories */
	if (Ra_pcaInputData) {
		for (wsUint i=0 ; i<N_pcaSample ; i++) {
			sseFree(Ra_pcaInputData[i]);
			sseFree(Ra_corr[i]);
			DEALLOC(Na_corrN[i]);
		}
		DEALLOC(Ra_pcaInputData);
		DEALLOC(Ra_corr);
		DEALLOC(Na_corrN);
		DEALLOC(Ra_ysum);
		DEALLOC(Ra_y2sum);
	}
}

wsReal** powerMethod(double **Ra_mat, double *Ra_eval, wsUint N_sz,
	int &N_pc, double R_thr, double R_varAll, double R_prop)
{
	wsUint i, j, k;
	wsReal **Ra_ev	= allocateMatrix(N_pc==-1?N_sz:N_pc, N_sz);
	double *ttt		= NULL;

	/* Matrix x (initial guess) */
	double *xk;	///< Matrix xk
	double m0;	///< x * x variable
	double m1;	///< x * xk variable
	double m2;	///< xk * xk variable
	double q = 0;	///< Calculated Raleigh quotient
	double e = 0;	///< Calculated error bound
	double max;	///< Maximum magnitude of xk

	sseMalloc(xk, double, N_sz);
	sseMalloc(ttt, double, N_sz);

	LOG("Initial : trace(A) = %g\n", R_varAll);
	wsReal R_explained = W0;

	int N_iter = 0;
	/* Calculate xk vector, Raleigh Quotient, and Error. */
	for ( ; N_iter<N_pc || !NA(R_prop) ; N_iter++) {
		// x0<-matrix(rep(1,ncol(dat)),ncol=1)
#ifdef USE_DBL
		sseVinit(ttt, N_sz, W1);
#else
		for (i=0 ; i<N_sz ; i++)
			ttt[i] = 1.0;
#endif
		double pe = -999;

		for (k=0 ; k<1000 ; k++) {
			/* Load matrix 'x' w/ zeros. */
			memset(xk, 0x00, sizeof(double)*N_sz);

			/* Set coeffs to zero */
			m0 = m1 = m2 = 0;
			max = 0;

			/* Get Ax */
			for(i=0; i<N_sz ; i++) {
#ifdef USE_DBL
				xk[i] = sseVV(Ra_mat[i], N_sz, ttt);
#else
 				for(j=0; j<N_sz; j++) /* Calculate A * x. */
 					xk[i] += Ra_mat[i][j]*ttt[j];
#endif
				m0 += ttt[i]*ttt[i];   /* x * x vector   */
				m1 += ttt[i]*xk[i];  /* x * xk vector  */
				m2 += xk[i]*xk[i]; /* xk * xk vector */

				/* Find max magnitude of xk vector */
				if (xk[i]>max) max=xk[i];
				else if (0-xk[i]>max) max=0-xk[i];
			}

			/* Display x(k) vector */
#ifdef USE_DBL
			sseVpC(xk, W1/max, ttt, N_sz);
#else
 			for(i=0; i<N_sz ; i++)
 				ttt[i] = xk[i]/max;
#endif

			q = m1/m0; /* Calculate Raleigh quotient */
			e = sqrt((m2/m0)-(q*q)); /* Calculate error bound */

			if (pe != -999) {
//				LOG("Loop %d e %g...\r", k, e);
				if (e < R_thr) {
// 					LOG("\nx%d=[",k);
// 					for(i=0; i<N ; i++)
// 						LOG("%lf ",ttt[i]);
// 					if (N_sz!=N)
// 						LOG("...");
// 					LOG("] q%d=%g e=%g\n",q,N_iter,e);
					break;
				}
			} else
				pe = e;
		}
		wsReal R_varExp = q/R_varAll*100.0;
		R_explained += R_varExp;
		LOG("Iteration %d : (%s:%d) e.value=%g (%g%%), error=%g\n",
			N_iter+1, k==100?"STOP":"OK", k, q, R_varExp, e);
		/* Fill eigenvalue */
		if (Ra_eval) Ra_eval[N_iter] = q;
		/* If prop mode */
		if (!NA(R_prop) && R_explained >= R_prop)
			R_prop = WISARD_NAN;

		/* Substract */
#ifdef USE_DBL
		double xxx = sseVsqsum(ttt, N_sz);//0;
#else
		double xxx = 0.0;
 		for (j=0 ; j<N_sz ; j++)
 			xxx += ttt[j]*ttt[j];
#endif
		xxx = sqrt(xxx);
#ifdef USE_DBL
		sseVpC(ttt, W1/xxx, Ra_ev[N_iter], N_sz);
#else
 		for(i=0; i<N_sz ; i++)
 			Ra_ev[N_iter][i] = (wsReal)(ttt[i]/xxx);
#endif

		for (j=0 ; j<N_sz ; j++)
			ttt[j] /= xxx;
//		LOG("Extract matrix\n");
		for (i=0 ; i<N_sz ; i++) {
#ifdef USE_DBL
			double xx = sseVV(Ra_mat[i], N_sz, ttt);//0;
#else
			double xx = 0.0;
 			for (j=0 ; j<N_sz ; j++)
 				xx += Ra_mat[i][j]*ttt[j];
#endif

			/* Get vec */
			for (j=0 ; j<N_sz ; j++) {
				Ra_mat[i][j] -= xx*ttt[j];
				//				printf("%g ", A[i][j]);
			}
			//			printf("\n");
		}
	}
	LOG("[%d] PCs explained [%g%%] of entire variance\n", N_iter, R_explained);
	/* If N_pc == -1, adjust return matrix */
	if (N_pc == -1) {
		wsReal **Ra_newEV = NULL;
		wsAlloc(Ra_newEV, wsReal*, N_iter);
		/* Only retains N_iter PCs */
		for (int i=0 ; i<N_iter ; i++)
			Ra_newEV[i] = Ra_ev[i];
		/* Remove all others */
		for (int i=N_iter ; i<(int)N_sz ; i++)
			DEALLOC(Ra_ev[i]);
		DEALLOC(Ra_ev);
		/* Substitute to new */
		Ra_ev = Ra_newEV;
		/* Update # of PCs */
		N_pc = N_iter;
	}
	sseFree(ttt);
	sseFree(xk);

	return Ra_ev;
}

void cPcaAnalysisV2::run()
{
	cTimer	X_timer;
	int		N_spc = 0;

	if (IS_ASSIGNED(proppc)) {
		N_spc = -1;
	} if (!IS_ASSIGNED(npc)) {
		OPTION().assign("npc");
		OPTION().FORCE_OPT_NUMBER(npc);
		N_spc = OPT_NUMBER(npc);
	} else
		N_spc = OPT_NUMBER(npc);

	/* Allocate required memory */
	wsAlloc(Ra_expectedX, wsReal, N_pcaSample);

	/* Get correlation */
	X_timer.start();
	LOG("[STEP1] Correlation calculation...\n");

	if (!IS_ASSIGNED(cor)) {
		wsUint		N_vrt	= Cp_IO->sizeVariant();
		vVariant&	Xa_vrt	= Cp_IO->getVariant();
		
		/* Draw mask */
		wsFloat *Ra_corrMask = NULL;
		sseMalloc(Ra_corrMask, wsFloat, N_vrt);
		memset(Ra_corrMask, 0xff, sizeof(wsFloat)*N_vrt);
		wsRealCst *Ra_ppp = Cp_anaPPP->getPPP();
		wsUint N_noIncCor = 0;
		for (wsUint i=0 ; i<N_vrt ; i++)
			if (!isInRange(OPT_RANGE(cormaf), Ra_ppp[i]) ||
				isSexChromosome(Xa_vrt[i])) {
				Ra_corrMask[i] = CORR_CONST(0.0);
				N_noIncCor++;
			}

		/* If no SNPs are incorporated */
		if (N_noIncCor == N_vrt) {
			/* Halt if it is impossible to calculate empirical correlation
			 * or IBS matrix */
			if (!OPT_ENABLED(kinship) && !IS_ASSIGNED(cor)) {
	// 			LOG("Warning : There is no variant to available to calculate empirical "
	// 				"correlation, the range of --cormaf should be extended!\n");
				/* Force --indep */
				OPTION().assign("indep", "1", 1);
				OPTION().FORCE_OPT_NUMBER(indep);
			}
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
		}

#ifdef USE_CORRv2
		xCorrThread	X_ct;
		X_ct.Xp_SNP			= NULL;
		X_ct.Na_data		= NULL;
		X_ct.N_sample		= N_pcaSample;
		X_ct.N_variant			= Cp_IO->sizeVariant();
		X_ct.Na_corrN		= this->Na_corrN;
		X_ct.Ra_corr		= this->Ra_corr;
		X_ct.Ra_data		= this->Ra_pcaInputData;
		X_ct.B_isComplete	= B_isPcaDataComplete;
		X_ct.Ra_corrMask	= Ra_corrMask;//Cp_IO->getCorrMask();

		WORKER().run(cCorrAnalysis::calcCorrV2, forAllSample_corr, &X_ct, NULL);
#else
		WORKER().run(cCorrAnalysis::calcCorr, forAllSample_pca, this, NULL);
#endif
		LOG("%s taken\n", X_timer.getReadable());

		/* Export correlation matrix if verbose */
		if (OPT_ENABLED(verbose))
			exportMatrix("fcv", Ra_corr, N_pcaSample, N_pcaSample);
		sseFree(Ra_corrMask);
	}

	double		R_varAll	= 0.0;
	double**	Ra_pwmat	= NULL;
	wsAlloc(Ra_pwmat, double*, N_pcaSample);
	/* Do not perform standardize & covariance matrix induction */
	for (wsUint i=0 ; i<N_pcaSample ; i++) {
		wsAlloc(Ra_pwmat[i], double, N_pcaSample);
		for (wsUint j=0 ; j<N_pcaSample ; j++)
			Ra_pwmat[i][j] = Ra_corr[i][j];
		R_varAll += Ra_pwmat[i][i];
	}

	if (R_varAll == W0)
		halt("Cannot proceed anymore : Sample relatedness is complete wrong!");

	int N_pc = N_spc;
	if (N_pc > (int)N_pcaSample) {
		LOG("# of PC adjusted from %d to %d by # of pca-included sample size\n",
			N_pc, N_pcaSample);
		N_pc = (int)N_pcaSample;
	}

	/* Perform eigenreduction for variance-covariance matrix */
 	X_timer.start();
	//LOG("[STEP5] Eigenreduction...");
	wsReal *Ra_temp=NULL;

	/* For MDS */
	Ra_cmds = NULL;

	wsAlloc(Ra_temp, wsReal, N_pcaSample+2);
	if (OPT_ENABLED(fullpca)) {
//		Ra_cov = allocateMatrix(N_pcaSample, N_pcaSample);
		wsReal *Ra_eigValue=NULL;
		wsAlloc(Ra_eigValue, wsReal, N_pcaSample);

		/* Perform triangular decomposition */
		tred2_dbl(Ra_pwmat, N_pcaSample, Ra_eigValue, Ra_temp);
		// 	if (OPT_NUMBER(verbose)) {
		// 		cExporter* Cp_exporter = cExporter::summon("pp1"); /* CHECKED */
		// 		FILE *fp = Cp_exporter->getHandle();
		// 		fprintf(fp, "EVAL	INTERM\n");
		// 		for (i=0 ; i<N_founder ; i++)
		// 			fprintf(fp, "%g	%g\n", Ra_evals[i], Ra_interm[i]);
		// 	}
		/* Reduce matrix into symmetric & tridiagonal matrix */
		tqli_dbl(Ra_eigValue, Ra_temp, N_pcaSample, Ra_pwmat);
		if (OPT_NUMBER(verbose)) {
/**/		cExporter* Cp_pp2 = cExporter::summon("pp2"); /* CHECKED */

			Cp_pp2->put("IDX	EIGENVALUE\n");
			for (wsUint i=0 ; i<N_pcaSample ; i++)
				Cp_pp2->fmt("%d	%g\n", i+1, Ra_eigValue[i]);

			delete Cp_pp2;
		}
		LOG("%s taken\n", X_timer.getReadable());
		/* Transpose */
		transposeSelf(Ra_pwmat, N_pcaSample, N_pcaSample);

		/* Sort eigenvalues by descending order */
		xRealSort *Xa_eValSorted = buildRealSort(Ra_eigValue, N_pcaSample);

// 		if (OPT_NUMBER(verbose)) {
// 			exportMatrix("fev", Ra_pwmat, N_pcaSample, N_pcaSample);
// 		}
// 		LOG("%d eigenvectors : \n", N_pcaSample);
// 		for (i=0 ; i<N_pcaSample ; i++) {
// 			printf(" EV %d : ", i+1);
// 			for (j=0 ; j<N_pcaSample ; j++) {
// 				printf("%g ", Ra_cov[j][i]);
// 			}
// 			printf("\n");
// 		}
// 		printf("\n");

		/* ���밪���� �����ؾ� �Ѵ� */
		qsort(Xa_eValSorted, N_pcaSample, sizeof(xRealSort), sortabs_real);
		if (N_pcaSample > 1 && Xa_eValSorted[0].V < Xa_eValSorted[1].V)
			halt("Sort is incrorrect!");

		/* Re-order matrix */
		double **Ra_newMat = NULL;
		wsAlloc(Ra_newMat, double*, N_pcaSample);
		for (wsUint i=0 ; i<N_pcaSample ; i++)
			Ra_newMat[i] = Ra_pwmat[Xa_eValSorted[i].i];
		DEALLOC(Ra_pwmat);
		Ra_pwmat = Ra_newMat;

		/* For MDS */
		if (OPT_ENABLED(mds))
			Ra_cmds = sseMatrix(N_pcaSample, N_pcaSample);

// 		LOG("After sort : ");
// 		for (i=0 ; i<N_pcaSample ; i++)
// 			printf("%g(%d) ", Xa_eValSorted[i].V, Xa_eValSorted[i].i);
// 		printf("\n");

// 		LOG("%d eigenvectors : \n", N_pcaSample);
// 		for (i=0 ; i<N_pcaSample ; i++) {
// 			printf(" EV %d : ", i+1);
// 			for (j=0 ; j<N_pcaSample ; j++)
// 				printf("%g ", Ra_cov[j][Xa_eValSorted[i].i]);
// 			printf("\n");
// 		}
// 		printf("\n");

		/* Extract projections of N_pc principle components, it overwrites Ra_corr */
		X_timer.start();
		LOG("[STEP6] PC score calculation...");
		double *Ra_Csum = sseMsumC(Ra_pwmat, N_pc, N_pcaSample);;
		for (wsUint i=0 ; i<N_pcaSample ; i++) {
			/* Ra_interm�� i��° individual�� ���� ��� corr �� */
			memcpy(Ra_temp, Ra_corr[i], sizeof(wsReal)*N_pcaSample);

			/* PC ������ŭ ������ ���鼭 */
			for (int j=0 ; j<N_pc ; j++) {
				/* i��° individual�� corr�� * ���̰պ��� j��° -> i,j�� �� */
//				wsUint _j = Xa_eValSorted[j].i; /* _j��° �����͵���� ó���϶�� �ǹ��� */
				Ra_corr[i][j] = W0; /* i��° individual�� j��° pc�� ���� �� */

				for (wsUint k=0 ; k<N_pcaSample ; k++)
//					Ra_corr[i][j] += Ra_temp[k] * Ra_pwmat[k][_j];
					Ra_corr[i][j] += Ra_temp[k] * Ra_pwmat[j][k];

				/* For MDS */
				if (Ra_cmds)
					sseVpC(Ra_Csum, Xa_eValSorted[j].V, Ra_cmds[j], N_pcaSample);
			}
		}
		DEALLOC(Ra_Csum);
		DEALLOC(Xa_eValSorted);
		DEALLOC(Ra_eigValue);
//		deallocMatrix(Ra_cov, N_pcaSample);
		for (wsUint i=0 ; i<N_pcaSample ; i++) {
			DEALLOC(Ra_pwmat[i]);
		}
		DEALLOC(Ra_pwmat);
	} else {
		/* Space for eigenvalues */
		double*		Ra_eVal	= NULL;
		sseMalloc(Ra_eVal, double, N_pcaSample);
		double**	Ra_eVec	= powerMethod(Ra_pwmat, Ra_eVal, N_pcaSample, N_pc,
			10e-5, R_varAll, IS_ASSIGNED(proppc)?OPT_REAL(proppc):WISARD_NAN);
		if (OPT_ENABLED(verbose))
			exportMatrix("fev", Ra_eVec, N_pc, N_pcaSample);

		/* For MDS */
		if (OPT_ENABLED(mds))
			Ra_cmds = sseMatrix(N_pc, N_pcaSample);

		X_timer.start();
		LOG("[STEP6] PC score calculation...");
		wsVec Ra_Csum = sseMsumC(Ra_eVec, N_pc, N_pcaSample);;
		for (wsUint i=0 ; i<N_pcaSample ; i++) {
			/* Ra_interm�� i��° individual�� ���� ��� corr �� */
			memcpy(Ra_temp, Ra_corr[i], sizeof(wsReal)*N_pcaSample);

			/* PC ������ŭ ������ ���鼭 */
			for (int j=0 ; j<N_pc ; j++) {
//				Ra_corr[i][j] = 0.0; /* i��° individual�� j��° pc�� ���� �� */

#ifdef USE_DBL
				Ra_corr[i][j] = sseVV(Ra_temp, N_pcaSample, Ra_eVec[j]);

				/* For MDS */
				if (Ra_cmds) {
					sseVpC(Ra_Csum, Ra_eVal[j], Ra_cmds[j], N_pcaSample);
				}
#else
 				for (wsUint k=0 ; k<N_pcaSample ; k++)
 					Ra_corr[i][j] += Ra_temp[k] * Ra_eigValue[j][k];
#endif
			}
			DEALLOC(Ra_pwmat[i]);
		}
		sseFree(Ra_Csum);
		DEALLOC(Ra_pwmat);
		deallocMatrix(Ra_eVec, N_pc);
	}
	DEALLOC(Ra_temp);
	LOG("%s taken\n", X_timer.getReadable());

	/* Export PC score if verbose */
	if (OPT_NUMBER(verbose)) {
/**/	cExporter* Cp_pca = cExporter::summon("pca"); /* CHECKED */

		LOG("Export PC score by --verbose\n");

#if 0
		Cp_pca->put("#PC	FID	IID	PC_SCORE\n");
		for (i=0 ; i<N_pc ; i++)
			for (j=0 ; j<N_pcaSample; j++)
				Cp_pca->fmt("%d	%s	%s	%g\n", i+1,
					Sa_fidFounder[j].c_str(), Sa_iidFounder[j].c_str(),
					Ra_corr[j][i]);
#else
		Cp_pca->put("FID	IID");
		for (int i=0 ; i<N_pc ; i++)
			Cp_pca->fmt("	PC%d", i+1);
		Cp_pca->put("\n");
		for (wsUint j=0 ; j<N_pcaSample; j++) {
			Cp_pca->fmt("%d	%s	%s", j+1,
				Sa_fidFounder[j].c_str(), Sa_iidFounder[j].c_str());
			for (int i=0 ; i<N_pc ; i++)
				Cp_pca->fmt("	%g", Ra_corr[j][i]);
			Cp_pca->put("\n");
		}
#endif

		delete Cp_pca;
	}
	
	/* Build full PC matrix using pedigree structure */
	wsMat Ra_pc = _stepFullPCmatrix(N_pc);

	/* Export MDS plot points */
	if (Ra_cmds) {
// 		cExporter *Cp_mds = cExporter::summon("mds.res");
// 		delete Cp_mds;
	}

	/* Add to covariates if required */
	if (OPT_ENABLED(pc2cov)) {
		vCovar&	Xa_cov	= getIO()->getCovInfo();
		wsMat	Ra_cov	= getIO()->getCovariates();
		wsUint	N_cov	= getIO()->sizeCovar();
		wsUint	N_samp	= getIO()->sizeSample();

		wsUint	N_nCov	= N_cov + N_pc;
		wsMat	Ra_nCov	= NULL;
		wsAlloc(Ra_nCov, wsReal*, N_nCov);

		/* Copy original */
		for (wsUint i=0 ; i<N_cov ; i++)
			Ra_nCov[i] = Ra_cov[i];
		DEALLOC(Ra_cov); /* Deallocate original */

		/* Copy newly computed PCs into covariates */
		for (int i=0 ; i<N_pc ; i++) {
			/* Add to covar */
			xCovar X_new;
			X_new.N_idx			= N_cov+i;
			X_new.N_idxFile		= -1; /* Not originated from file */
			X_new.N_szFactor	= 0; /* Set to 0 since it is NOT factor */
			X_new.Sp_bl			= NULL;
			char S_pcName[64];
			sprintf(S_pcName, "PC%d", i+1);
			X_new.Sp_varName	= strdup(S_pcName);
			X_new.X_type		= WISARD_VAR_REAL;
			Xa_cov.push_back(X_new);

			/* Copy to covar */
			sseMalloc(Ra_nCov[N_cov+i], wsReal, N_samp);
			for (wsUint j=0; j<N_samp ; j++)
				Ra_nCov[N_cov+i][j] = Ra_pc[j][i];
		}

		getIO()->chgCovInfo(N_nCov, Ra_nCov);
	}
	/* Deallocates expired memories */
	wsUint	N_origSample	= Cp_IO->getSampleSizeBeforeFiltering();
	for (wsUint i=0 ; i<N_origSample ; i++) {
		DEALLOC(Ra_pc[i]);
	}
	DEALLOC(Ra_pc);

	DEALLOC(Ra_expectedX);
}

void cPcaAnalysisV2::_stepPolyEffectMat(wsReal **Rp_poly)
{
	halt("Unimplemented");
}

int cPcaAnalysisV2::_stepFullPCmatrix_MF(vStr_it &iit, mNucFam_it &it,
	wsMat Ra_pc, wsMat Ra_mds, wsUint *Np_fillCount, wsUint N_srcIdx,
	wsUint N_origSample, wsUint N_tmpIdxForMissingSamp, wsUint N_pc)
{
	/* Make concatenated identifier */
	string	S_keyCouple = it->second.Xp_pat->S_IID + "," +
		it->second.Xp_mat->S_IID;
	/* Set the copy size */
	wsUint	N_szCopy = sizeof(wsReal)*N_pc;

	if (N_srcIdx >= N_origSample)
		halt("Invalid index given #0 (%d>%d)\n", N_srcIdx, N_origSample);

	/* If both parents are not founder, skip */
	if ((it->second.Xp_mat->Xp_mat || it->second.Xp_mat->Xp_pat)
		&& (it->second.Xp_pat->Xp_mat || it->second.Xp_pat->Xp_pat))
		return 0;
	/* If cannot found, also skip */
	if (iit->compare(S_keyCouple))
		return 0;

	/* Set to 'both missing' if they do not have any data */
	char B_bothMissing = (it->second.Xp_mat->N_idx == SAMP_NODATA &&
		it->second.Xp_pat->N_idx == SAMP_NODATA);

	/* DO log */
	//LOG("Group [%s] filling : %d childs%s\n", iit->c_str(),
	//	it->second.Xp_childs.size(), B_bothMissing?" (bothMiss)":"");

	/* Propagate calculated PC score to the offsprings have its parent
	 * as missing (both or single */
	if (B_bothMissing) {
		/* Both founders are missing */
		FOREACH (vSampPtr_it, it->second.Xp_childs, cit) {
			/* Missing�� ��� �ӽ� ��ġ�� ������ �д� */
			if ((*cit)->B_isMissing) {
				if (N_tmpIdxForMissingSamp >= N_origSample)
					halt("Invalid index given #1\n");
				memcpy(Ra_pc[N_tmpIdxForMissingSamp],
					Ra_corr[N_srcIdx], N_szCopy);
				if (Ra_mds) for (wsUint j=0; j<N_pc ; j++)
					Ra_mds[N_tmpIdxForMissingSamp][j]
						= Ra_cmds[j][N_srcIdx];

				(*cit)->N_isVisited = short(N_tmpIdxForMissingSamp + 1);
				N_tmpIdxForMissingSamp++;
			} else {
				if ((*cit)->N_idx >= (int)N_origSample || (*cit)->N_idx < 0)
					halt("Invalid index given #2\n");
				memcpy(Ra_pc[(*cit)->N_idx], Ra_corr[N_srcIdx],
					N_szCopy);
				/* MDS filling */
				if (Ra_mds) for (wsUint j=0; j<N_pc ; j++)
					Ra_mds[(*cit)->N_idx][j] = Ra_cmds[j][N_srcIdx];

				(*cit)->N_isVisited = short((*cit)->N_idx + 1);
				//LOG("%s (bothmissing) filled\n", (*cit)->S_IID.c_str());
				(*Np_fillCount)++;
			}
		}
	} else {
		/* Single missing�� ��� */
		wsReal	*Ra_newPC			= NULL;
		wsReal	*Ra_newMds			= NULL;
		xSample	*Xp_availableSrc	= it->second.Xp_mat->B_isMissing ?
			it->second.Xp_pat : it->second.Xp_mat;
		wsAlloc(Ra_newPC, wsReal, N_pc);

		/* For MDS */
		if (Ra_mds)
			wsAlloc(Ra_newMds, wsReal, N_origSample);

		/* Missing�� �ƴ� ������ Ra_corr�� ��� ����ִ��� Ȯ�� */
		wsUint N_aSideSrcIdx = X_iid2fouIdx[Xp_availableSrc->S_IID];

		/* ���ھ� ����� */
		for (wsUint j=0 ; j<N_pc ; j++) {
			Ra_newPC[j] = W2*Ra_corr[N_srcIdx][j] -
				Ra_corr[N_aSideSrcIdx][j];
			if (Ra_mds) Ra_newMds[j] = W2*Ra_cmds[j][N_srcIdx] -
				Ra_cmds[j][N_aSideSrcIdx];
		}

		FOREACH (vSampPtr_it, it->second.Xp_childs, cit) {
			/* Missing�� ��� �ӽ� ��ġ�� ������ �д� */
			if ((*cit)->B_isMissing) {
				if (N_tmpIdxForMissingSamp >= N_origSample)
					halt("Invalid index given #3\n");
				memcpy(Ra_pc[N_tmpIdxForMissingSamp],
					Ra_newPC, N_szCopy);
				if (Ra_mds) memcpy(Ra_mds[N_tmpIdxForMissingSamp],
					Ra_newMds, N_szCopy);

				(*cit)->N_isVisited = short(N_tmpIdxForMissingSamp + 1);
				N_tmpIdxForMissingSamp++;
			} else {
				if ((*cit)->N_idx >= (int)N_origSample || (*cit)->N_idx < 0)
					halt("Invalid index given #4\n");
				memcpy(Ra_pc[(*cit)->N_idx], Ra_newPC, N_szCopy);

				(*cit)->N_isVisited = short((*cit)->N_idx + 1);
				//LOG("%s (singmissing) filled\n", (*cit)->S_IID.c_str());
				(*Np_fillCount)++;
			}
		}

		DEALLOC(Ra_newPC);
		DEALLOC(Ra_newMds);
	}

	return 1;
}

wsMat cPcaAnalysisV2::_stepFullPCmatrix(wsUint N_pc)
{
	wsUint	i = 0;
	wsUint	N_sample		= Cp_IO->sizeSample();
	wsUint	N_origSample	= Cp_IO->getSampleSizeBeforeFiltering();
	mNucFam	&Xm_nfData		= Cp_anaFamStr->getNucFamData();
	cTimer	X_timer;

	/* Allocate full PC matrix */
	wsMat	Ra_pc			= allocateMatrix(N_origSample, N_pc);
	wsMat	Ra_mds			= NULL;
	/* For MDS */
	if (OPT_ENABLED(mds))
		Ra_mds = sseMatrix(N_origSample, N_pc);

	/* Fill PC matrix with NaN */
	for (i=0 ; i<N_origSample ; i++) {
		for (wsUint j=0 ; j<N_pc ; j++)
			Ra_pc[i][j] = WISARD_NAN;
	}

	/* Ra_corr [#pcaSamp * #pc] ������ ������ ����
	 *  - i��° ���� Sa_iidFounder[i] �� PC score�� ������ �ִ�
	 *  - j��° ���� j��° PC score�� ������ �ִ�
	 */

	/* Make PC matrix using family structure */
	/* PCA �м��� ���� ��� sample ������ Ȯ���ϰ�, �̵鿡�� üũ �ο� */
	/* B_isVisited ���
	 * ���� �ǹ� : 0�� �ƴϸ� �ش� ��-1�� PCA matrix������ �ش� ���ÿ� �����Ǵ� PC score */
	mSamp	&Xa_sampleData
						= Cp_IO->getSampleData();
	wsUint N_tmpIdxForMissingSamp
						= N_sample;
	wsUint N_srcIdx		= 0;
	wsUint N_fillCount	= 0;

	/* For all founders included in PCA */
	i = 0;
	for (vStr_it iit=Sa_iidFounder.begin(),fit=Sa_fidFounder.begin() ;
		iit!=Sa_iidFounder.end() ; iit++,fit++,i++) {

		/* �� ���� Missing founder ���̳� �ƴϳĿ� ���� ó���� �޶����� �Ѵ� */
		if (fit->compare("MF") == 0) {
			if (OPT_ENABLED(usemf)) {
				xSample	&X_sample = Xa_sampleData[*iit];

				/* usemf Ȱ��ȭ ���¿����� missing founder���� �ϴ� üũ�� �� ���� */
				X_sample.N_isVisited = short(N_tmpIdxForMissingSamp+1);
			} else if (!OPT_ENABLED(founderonly)) {
				/* usemf ��Ȱ��ȭ ���¿����� ���� �ٰ��� ����Ʈ���� ��� */
				FOREACH (mNucFam_it, Xm_nfData, it) {
					if (_stepFullPCmatrix_MF(iit, it, Ra_pc, Ra_mds,
						&N_fillCount, N_srcIdx, N_origSample,
						N_tmpIdxForMissingSamp, N_pc))
						N_srcIdx++;
				}
			}
			N_tmpIdxForMissingSamp++;
			continue;
		} else {
			xSample	&X_sample = Xa_sampleData[*iit];

			wsUint N_insIdx = X_sample.N_idx;
			wsUint N_srcIdx = i;

			if (N_insIdx >= N_origSample)
				halt("Invalid index [%d] given, exceeds the total [%d]",
					N_insIdx, N_origSample);

			X_sample.N_isVisited = short(N_insIdx+1);
			memcpy(Ra_pc[N_insIdx], Ra_corr[N_srcIdx], sizeof(wsReal)*N_pc);
			/* For MDS */
			if (Ra_cmds) for (wsUint j=0 ; j<N_pc ; j++)
				Ra_mds[N_insIdx][j] = Ra_cmds[j][N_srcIdx];
			pverbose("%s [%d->%d] filled\n", X_sample.S_IID.c_str(),
				N_srcIdx, N_insIdx);
			//LOG("%s (have) filled\n", X_sample.S_IID.c_str());
			N_fillCount++;
		}
	}
	LOG("Initial filling status [%d filled]\n", N_fillCount);
	N_tmpIdxForMissingSamp = N_sample;

	/* Fill Ra_pc matrix with its PC score, in iterative manner.
	 * Maximum iteration count will be the maximum generation count in dataset */
	/* PCA�� ���� ��� sample���� �̹� �Ҵ��� �� �Ǿ����� ���� */
	/* ���� �� ����Ʈ�� ���Ե� �ڽĵ��� ������� �����ؾ� �Ѵ� */
	mSampPtr	Xa_currTarget;
	mSamp		&Xm_sampleData = Cp_IO->getSampleData();

	/* X_iid2fouIdx�� ���Ե� sample���� �ڽ����� �ʱ� Ÿ�� ���� */
	FOREACH (strMap_it, X_iid2fouIdx, it) {
		/* xSample ����ü Ž�� */
		xSample &X_currSample = Xm_sampleData[it->first];

		/* ���� sample�� �ڽĵ��� Xa_currTarget�� ���� */
		FOREACH (vSampPtr_it, X_currSample.Xp_childs, cit)
			Xa_currTarget.insert(make_pair((*cit)->S_IID, *cit));
	}

	/* Xa_currTarget�� ������� ���� �� ��ŷ�� ����� �����ִ� ���̹Ƿ� ���� */
	mSampPtr	Xa_newTarget;
	wsUint		N_szTarget = 0;
	wsUint		N_stop = 0;
	for (int i=1 ; Xa_currTarget.size() ; i++) {
		wsUint	N_filledThisLoop = 0;
		LOG("Filling loop %d\n", i);

		FOREACH (mSampPtr_it, Xa_currTarget, it) {
			/* xSample ����ü Ž�� */
			xSample		*Xp_samp = it->second;
			wsUint		N_insertSelf = 1;
			
//			LOG("Check %s[%d] mot %d fat %d\n", Xp_samp->S_IID.c_str(), Xp_samp->N_isVisited,
//				Xp_samp->Xp_mat->N_isVisited, Xp_samp->Xp_pat->N_isVisited);
			/* �湮 ���� ���߰� �θ� �� ä������� ä��⸦ ���� */
			if (Xp_samp->N_isVisited == 0) {
				/* Make same PC score IF only one parents is available */
				if ((Xp_samp->Xp_mat && !Xp_samp->Xp_pat) ||
					(!Xp_samp->Xp_mat && Xp_samp->Xp_pat)) {
					wsUint N_idx = Xp_samp->Xp_pat ?
						Xp_samp->Xp_pat->N_isVisited-1 :
						Xp_samp->Xp_mat->N_isVisited-1;

					/* Missing�̸� �ӽ� ������ ä��� */
					if (Xp_samp->B_isMissing) {
						/* Make PC score of this node as the average of
						 * PC score of two parents */
						for (wsUint k=0 ; k<N_pc ; k++) {
							Ra_pc[N_tmpIdxForMissingSamp][k]
								= Ra_pc[N_idx][k];
							if (Ra_cmds)
								Ra_mds[N_tmpIdxForMissingSamp][k]
									= Ra_mds[N_idx][k];
						}

						it->second->N_isVisited = short(N_tmpIdxForMissingSamp+1);
						N_tmpIdxForMissingSamp++;
					} else {
						/* Make PC score of this node as the average of
						 * PC score of two parents */
						for (wsUint k=0 ; k<N_pc ; k++) {
							Ra_pc[Xp_samp->N_idx][k]
								= Ra_pc[N_idx][k];
							if (Ra_cmds)
								Ra_mds[Xp_samp->N_idx][k]
									= Ra_mds[N_idx][k];
						}
//						LOG("%s filled\n", it->second->S_IID.c_str());
						it->second->N_isVisited = short(Xp_samp->N_idx+1);
						N_fillCount++;
						N_filledThisLoop++;
					}
					N_insertSelf = 0;

					pverbose(" > %s filled [at %d]\n", Xp_samp->S_IID.c_str(), Xp_samp->N_idx);
				}
				/* Fill PC score of this node IF both parents are filled */
				else if (Xp_samp->Xp_mat->N_isVisited && 
					Xp_samp->Xp_pat->N_isVisited) {
					wsUint N_pi = Xp_samp->Xp_pat->N_isVisited-1;
					wsUint N_mi = Xp_samp->Xp_mat->N_isVisited-1;

					/* Missing�̸� �ӽ� ������ ä��� */
					if (Xp_samp->B_isMissing) {
						/* Make PC score of this node as the average of
						 * PC score of two parents */
						for (wsUint k=0 ; k<N_pc ; k++) {
							Ra_pc[N_tmpIdxForMissingSamp][k]
								= (Ra_pc[N_pi][k] + Ra_pc[N_mi][k]) /
									W2;
								if (Ra_cmds)
									Ra_mds[N_tmpIdxForMissingSamp][k]
										= (Ra_mds[N_pi][k] + Ra_mds[N_mi][k]) /
											W2;
						}

						it->second->N_isVisited = short(N_tmpIdxForMissingSamp+1);
						N_tmpIdxForMissingSamp++;
					} else {
						/* Make PC score of this node as the average of
						 * PC score of two parents */
						for (wsUint k=0 ; k<N_pc ; k++) {
							Ra_pc[Xp_samp->N_idx][k]
								= (Ra_pc[N_pi][k] + Ra_pc[N_mi][k]) /
									W2;
							if (Ra_cmds)
								Ra_mds[Xp_samp->N_idx][k]
									= (Ra_mds[N_pi][k] + Ra_mds[N_mi][k]) /
										W2;
						}

//						LOG("%s filled\n", it->second->S_IID.c_str());
						it->second->N_isVisited = short(Xp_samp->N_idx+1);
						N_fillCount++;
						N_filledThisLoop++;
					}
					N_insertSelf = 0;
					pverbose(" > %s filled [at %d]\n", Xp_samp->S_IID.c_str(), Xp_samp->N_idx);
				}
			} else
				N_insertSelf = 0;

			if (N_insertSelf == 0) {
				FOREACH (vSampPtr_it, it->second->Xp_childs, cit) {
//					LOG("%s newly inserted\n", (*cit)->S_IID.c_str());
					Xa_newTarget[(*cit)->S_IID] = *cit;
				}
			} else
				Xa_newTarget[it->second->S_IID] = it->second;
		}
		/* Clear current target */
		Xa_currTarget.clear();

		if (N_filledThisLoop == 0) {
			LOG("Nothing were filled\n");

			/* Check if nothing were filled and target size is also same */
			if (Xa_currTarget.size() == N_szTarget) {
				N_stop++;
				if (N_stop == 5) {
					LOG("Seems to there is not enough information to fill\n"
						"Remained samples will not have their PC score\n");
					break;
				}
			}
		} else
			LOG("%d filled\n", N_filledThisLoop);

		/* Re-make current target */
		N_szTarget = (wsUint)Xa_currTarget.size();
		copy(Xa_newTarget.begin(), Xa_newTarget.end(),
			inserter(Xa_currTarget, Xa_currTarget.begin()));
		Xa_newTarget.clear();
	}

	/* Export induced PCs if verbose */
	if (1) {
/**/	cExporter*	Cp_fpc = cExporter::summon("pca.res"); /* CHECKED */
		vSampPtr	&Xa_sample	= Cp_IO->getSample();

		Cp_fpc->put("FID	IID");
		for (i=0 ; i<N_pc ; i++)
			Cp_fpc->fmt("	PC%d", i+1);
		Cp_fpc->put("\n");
		for (i=0 ; i<N_sample ; i++) {
			Cp_fpc->fmt("%s	%s", Xa_sample[i]->S_FID.c_str(), Xa_sample[i]->S_IID.c_str());
			for (wsUint j=0 ; j<N_pc ; j++) {
				if (NA(Ra_pc[i][j]))
					Cp_fpc->put("	NA");
				else
					Cp_fpc->fmt("	%g", Ra_pc[i][j]);
			}
			Cp_fpc->put("\n");
		}
		delete Cp_fpc;

		/* For MDS */
		if (OPT_ENABLED(mds)) {
			cExporter*	Cp_mds = cExporter::summon("mds.res"); /* CHECKED */

			Cp_mds->put("FID	IID");
			for (i=0 ; i<N_pc ; i++)
				Cp_mds->fmt("	MDS%d", i+1);
			Cp_mds->put("\n");
			for (i=0 ; i<N_sample ; i++) {
				Cp_mds->fmt("%s	%s", Xa_sample[i]->S_FID.c_str(), Xa_sample[i]->S_IID.c_str());
				for (wsUint j=0 ; j<N_pc ; j++) {
					if (NA(Ra_mds[i][j]))
						Cp_mds->put("	NA");
					else
						Cp_mds->fmt("	%g", Ra_mds[i][j]);
				}
				Cp_mds->put("\n");
			}

			delete Cp_mds;
		}
	}

	/* Build polygenic effects */
// 	X_timer.start();
// 	LOG("[STEP8] Make polygenic effect matrix..."); 
// 	_stepPolyEffectMat(Ra_poly);
// 	LOG("%s taken\n", X_timer.getReadable());

	X_timer.start();

	/* Remove visitation markers */
	mFam &Xa_fam = Cp_IO->getFamilyData();
	FOREACH (mFam_it, Xa_fam, it) {
		it->second.setUnvisited();
	}

	return Ra_pc;
}

wsUint cPcaAnalysisV2::getPcaSampleSize()
{
	return N_pcaSample;
}

int cPcaAnalysisV2::calcEx(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	cPcaAnalysisV2
			*Cp_anaVC	= (cPcaAnalysisV2 *)Vp_shareData;		///< Shared data is its instance
	cIO		*Cp_IO		= Cp_anaVC->getIO();
	wsUint	N_founder	= Cp_IO->sizeFounder();				///< # of founders in dataset
	wsUint	N_idxSamp	= *((int *)Vp_data);					///< Thread-specific data is index of founder for this thread
	wsReal	R_sum		= 0.0f;									///< Sum of genotype
	int		N_incSamp	= 0;									///< # of founders that included to summation
	wsReal	*Ra_fData	= Cp_anaVC->Ra_corr[N_idxSamp];
	wsReal	*Rp_eX		= ((wsReal *)Vp_result) + N_idxSamp;	///< E(X) storing address

	if (Cp_IO->isDataComplete() == 1) {
		/* N_incSample == N_founder if data is perfect */
		for (wsUint i=0 ; i<N_founder ; i++)
			R_sum += Ra_fData[i];
		N_incSamp = N_founder;
	} else {
		/* Otherwise, it should be counted to get accurate E(X) */
		for (wsUint i=0 ; i<N_founder ; i++)
			if (isAvailableReal(Ra_fData[i])) {
				R_sum += Ra_fData[i];
				N_incSamp++;
			}
	}
	*Rp_eX = R_sum / (wsReal)N_founder;

	return 0;
}

int cPcaAnalysisV2::calcXeX(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	cPcaAnalysisV2
			*Cp_anaVC	= (cPcaAnalysisV2 *)Vp_shareData;		///< Shared data is its instance
	cIO		*Cp_IO		= Cp_anaVC->getIO();
	wsUint	N_founder	= Cp_IO->sizeFounder();				///< # of founders in dataset
	int		N_idxSamp	= *((int *)Vp_data);					///< Thread-specific data is index of founder for this thread
	wsReal	*Ra_fData	= Cp_anaVC->Ra_corr[N_idxSamp];
	wsReal	*Rp_eX		= Cp_anaVC->Ra_expectedX;

	if (Cp_IO->isDataComplete() == 1) {
		for (wsUint i=0 ; i<N_founder ; i++)
			Ra_fData[i] -= Rp_eX[i];
	} else {
		/* Missing-considered */
		for (wsUint i=0 ; i<N_founder ; i++)
			if (isAvailableReal(Ra_fData[i]))
				Ra_fData[i] -= Rp_eX[i];
	}

	return 0;
}

int cPcaAnalysisV2::calcCov(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	cPcaAnalysisV2
			*Cp_anaVC	= (cPcaAnalysisV2 *)Vp_shareData;
	cIO		*Cp_IO		= Cp_anaVC->getIO();
	wsReal	**Ra_fData	= Cp_anaVC->Ra_corr;
	wsUint	N_idxSamp	= *((int *)Vp_data);
	wsUint	N_founder	= Cp_IO->sizeFounder();

	for (wsUint j=N_idxSamp ; j<N_founder ; j++) {
		Cp_anaVC->Ra_cov[N_idxSamp][j] = 0.0f;
		// S <- S + (m[k,i]*m[k,j])
		for (wsUint i=0 ; i<N_founder ; i++)
			Cp_anaVC->Ra_cov[N_idxSamp][j] += Ra_fData[i][N_idxSamp]*Ra_fData[i][j];
		// ret[i,j] <- S/(d[1]-1)
		Cp_anaVC->Ra_cov[N_idxSamp][j] /= (N_founder-1);
		Cp_anaVC->Ra_cov[j][N_idxSamp] = Cp_anaVC->Ra_cov[N_idxSamp][j];
	}

	return 0;
}

#endif

} // End namespace ONETOOL
