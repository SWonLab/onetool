#include <stdlib.h>
#include "analyses/ppp.h"
#include "analyses/fam.h"
#include "analyses/annot.h"
#include "utils/matrix.h"

namespace ONETOOL {

template <typename T>
void _nrmIncomplete(wsFloat **Ra_data, T** Ta_data, wsUint j,
	wsUint N_SNP, wsUint N_sex, wsReal *Ra_sum, wsUint *Na_ind, vVariant &Xa_SNP,
	wsReal *Ra_ppp, wsReal *Ra_ppp_sqrt2pq, char *Na_cm)
{
	wsFloat	*Rp_data	= Ra_data[j];
	T*Np_data	= Ta_data[j];
	wsReal	R_mulX		= N_sex==1 ? W1 : W2;
	wsReal	R_mulY		= N_sex==1 ? W1 : W0;

#ifdef USE_CORR00
#	define NORM_MISVAL	0
#else
#	define NORM_MISVAL	WISARD_NA_REAL
#endif

	for (wsUint i=0 ; i<N_SNP ; i++) {
		if (!Na_cm[i]) continue;

		if (isXChromosome(Xa_SNP[i])) {
			if (isAvailable(Np_data[i])) {
				Rp_data[i]	= (wsFloat)(((wsReal)Np_data[i] - R_mulX*Ra_ppp[i])
					/ Ra_ppp_sqrt2pq[i]);
				Ra_sum[i]	+= (wsReal)Rp_data[i];
				Na_ind[i]++;
			} else
				Rp_data[i]	= NORM_MISVAL;
		} else if (isYChromosome(Xa_SNP[i])) {
			if (isAvailable(Np_data[i])) {
				Rp_data[i]	= (wsFloat)(((wsReal)Np_data[i] - R_mulY*Ra_ppp[i])
					/ Ra_ppp_sqrt2pq[i]);
				Ra_sum[i]	+= (wsReal)Rp_data[i];
				Na_ind[i]++;
			} else
				Rp_data[i]	= NORM_MISVAL;
		} else {
			if (isAvailable(Ta_data[j][i])) {
				Rp_data[i]	= (wsFloat)(((wsReal)Np_data[i] - W2*Ra_ppp[i])
					/ Ra_ppp_sqrt2pq[i]);
				Ra_sum[i]	+= (wsReal)Rp_data[i];
				Na_ind[i]++;
			} else
				Rp_data[i]	= NORM_MISVAL;
		}
	}
}

template <typename T>
void _nrmComplete(wsFloat **Ra_data, T** Ta_data, wsUint j,
	wsUint N_SNP, wsUint N_sex, wsReal *Ra_sum, wsUint *Na_ind, vVariant &Xa_SNP,
	wsReal *Ra_ppp, wsReal *Ra_ppp_sqrt2pq, char *Na_cm)
{
	wsFloat*	Rp_data	= Ra_data[j];
	T*		Np_data	= Ta_data[j];
	wsReal	R_mulX	= N_sex==1 ? W1 : W2;
	wsReal	R_mulY	= N_sex==1 ? W1 : W0;

	for (wsUint i=0 ; i<N_SNP ; i++) {
		if (!Na_cm[i]) continue;

		if (isXChromosome(Xa_SNP[i])) {
			Rp_data[i]	= (wsFloat)(((wsReal)Np_data[i] - R_mulX*Ra_ppp[i])
				/ Ra_ppp_sqrt2pq[i]);
			Ra_sum[i]	+= (wsReal)Rp_data[i];
		} else if (isYChromosome(Xa_SNP[i])) {
			Rp_data[i]	= (wsFloat)(((wsReal)Np_data[i] - R_mulY*Ra_ppp[i])
				/ Ra_ppp_sqrt2pq[i]);
			Ra_sum[i]	+= (wsReal)Rp_data[i];
		} else {
			Rp_data[i]	= (wsFloat)(((wsReal)Np_data[i] - W2*Ra_ppp[i])
				/ Ra_ppp_sqrt2pq[i]);
			Ra_sum[i]	+= (wsReal)Rp_data[i];
		}
	}
}

cPPPAnalysisV2::cPPPAnalysisV2(cIO *Cp_inpIO,
	cFamStrAnalysis *Cp_inpAnaFamStr, char B_inpDoNorm) : cAnalysis(Cp_inpIO)
{
	wsUint	N_SNP	= Cp_IO->sizeVariant();
	Cp_anaFamStr	= Cp_inpAnaFamStr;
	B_doNorm		= B_inpDoNorm;

	/* Check option validity */
	if (OPT_ENABLED(founderonly) && OPT_ENABLED(usemf))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--founderonly", "--usemf");

	Ra_misFounderGeno = NULL;
	if (IS_ASSIGNED(dosage) && OPT_ENABLED(corpearson)) {
		LOG("--dosage with --corpearson, PPP does not required so pass\n");
		Ra_ppp = Ra_ppp_sqrt2pq = Ra_meanX = NULL;
	} else {
		wsAlloc(Ra_ppp		   , wsReal, N_SNP);
		sseMalloc(Ra_ppp_sqrt2pq, wsReal, N_SNP);
		wsAlloc(Ra_meanX	   , wsReal, N_SNP);
	}
	Ra_corrMask = NULL;
	Na_corrMask = NULL;
}

cPPPAnalysisV2::~cPPPAnalysisV2()
{
	deallocMatrix(Ra_misFounderGeno, (wsUint)Xm_misFounder.size());
	DEALLOC(Ra_ppp		  );
	sseFree(Ra_ppp_sqrt2pq);
	DEALLOC(Ra_meanX	  );
	sseFree(Ra_corrMask);
	sseFree(Na_corrMask);
}

void cPPPAnalysisV2::run()
{
	/* Can't do when expression */
	if (IS_ASSIGNED(expression)) {
		LOGwarn("Skip PPP calculation for expression dataset\n");
		return;
	}
	/* --cormaf if not set */
	if (!IS_ASSIGNED(cormaf)) {
		OPTION().assign("cormaf", OPTION().getDefVal("cormaf"), 1);
		FORCE_OPT_RANGE(cormaf);
	}

	mFam	&Xa_family	= Cp_IO->getFamilyData();
	wsUint	N_SNP		= Cp_IO->sizeVariant();

	/* Step 1 : Find missing founders from family structure */
	FOREACH (mFam_it, Xa_family, it) {
		vSampPtr &Xp_founders = it->second.Xp_founders;
		FOREACH (vSampPtr_it, Xp_founders, fit) {
			if ((*fit)->B_isMissing)
				Xm_misFounder.insert(make_pair((*fit)->S_IID, *fit));
		}
	}
	LOG("%d missing founders were found\n", Xm_misFounder.size());

	/* Step 2 : Impute missing founder's genotype data from its nuclear families */
	if (OPT_ENABLED(usemf)) {
		wsUint i=0;
		wsAlloc(Ra_misFounderGeno, wsReal*, Xm_misFounder.size());
		LOG("--usemf assigned, performing imputation of missing founders...\n");

		/* Impute genotype for each missing founder */
		FOREACHDO (mSampPtr_it, Xm_misFounder, it, i++) {
			wsCalloc(Ra_misFounderGeno[i], wsReal, N_SNP);

			/* In this point, using Na_data */
			_doImpute(it->second, i);

			Xm_misFounderDtIdx.insert(make_pair(it->first, i));
			notice("%d/%d missing founders imputed\r", i+1, Xm_misFounder.size());
		}

		LOGnote("[%d] missing founders imputed\n", i);
		if (i != Xm_misFounder.size())
			halt_fmt(WISARD_SYST_INVL_EXPCOUNT, "# of imputed samples", i,
				"# of missing founders", Xm_misFounder.size());
	}

	/* Step 3 : Calculate PPP */
	/* In this point, using Na_data, Ra_data will made in this function */
	_doCalcPPP();

	/* Step 4 : Export if verbose */
	/* Export ppp values if VERBOSE */
	if (OPT_ENABLED(verbose)) {
		cTableExporter	C_ppp("ppp", "r", "PPP coefficient", 1, 1, "PPP");
		vVariant&		Xv_vrt	= Cp_IO->getVariant();
		wsUint			I		= 0;

		FOREACHDO (vVariant_it, Xv_vrt, i, I++)
			C_ppp.writeVariant(&(*i), Ra_ppp[I]);
	}

	/* Step 5 : Build mask */
	sseMalloc(Ra_corrMask, wsFloat, N_SNP);
	sseMalloc(Na_corrMask, char, N_SNP);
	memset(Ra_corrMask, 0xff, sizeof(wsFloat)*N_SNP);
	memset(Na_corrMask, 0xff, sizeof(char)*N_SNP);
	wsRealCst*	Ra_ppp		= getPPP();
	N_noIncCor	= 0;
	vVariant	&Xv_mk		= Cp_IO->getVariant();

	/* UPDATE 140114 : Do not include sex-chromosome */
	N_sexMarker = 0;
	for (wsUint i=0 ; i<N_SNP ; i++)
		if (!isInRange(OPT_RANGE(cormaf), Ra_ppp[i]) ||
			isSexChromosome(Xv_mk[i])) {
			if (isSexChromosome(Xv_mk[i]))
				N_sexMarker++;
			Ra_corrMask[i] = CORR_CONST(0.0);
			Na_corrMask[i] = 0x00;
			N_noIncCor++;
		}
}

wsReal _inferGeno(xVariant &X_snp, xSample *Xp_misFounder, vNucFamPtr& Xa_relNucFam,
	char **Na_data, wsUint &N, wsUint i)
{
	N = 0;
	wsReal R_inferGeno = W0;

	if (isXChromosome(X_snp)) {
		if (Xp_misFounder->N_sex == 1) {
			FOREACH (vNucFamPtr_it, Xa_relNucFam, nit) {
				xNucFamily &X_nf = *(*nit);

				/* Missing founder is male */
				char *Na_availData;
				if ((*nit)->Xp_mat->N_idx == SAMP_NODATA &&
					(*nit)->Xp_pat->N_idx == SAMP_NODATA)
					Na_availData = NULL;
				else
					Na_availData = X_nf.Xp_mat->N_idx != SAMP_NODATA ?
					Na_data[X_nf.Xp_mat->N_idx] :
				Na_data[X_nf.Xp_pat->N_idx];

				FOREACH (vSampPtr_it, (*nit)->Xp_childs, cit) {
					if ((*cit)->B_isMissing) continue;
					if ((*cit)->N_idx == SAMP_NODATA)
						halt_fmt(WISARD_SYST_NULL_CHILDGENO,
							(*cit)->S_FID.c_str(), (*cit)->S_IID.c_str(),
							X_nf.Xp_pat->S_IID.c_str(), X_nf.Xp_mat->S_IID.c_str());
					char *Na_offData = Na_data[(*cit)->N_idx];
					if (!isAvailable(Na_offData[i]))
						continue;
					/* 아빠의 X genotype을 추정할 때 아들은 무시한다 */
					if ((*cit)->N_sex == 1)
						continue;

					/* 엄마가 있으면 엄마의 효과를 제거 */
					if (Na_availData)
						R_inferGeno += W2*(wsReal)Na_offData[i]
					- (wsReal)Na_availData[i];
					else /* 아니면 1/2만 */
						R_inferGeno += (wsReal)Na_offData[i]/W2;
					N++;
				}
			}
		} else {
			FOREACH(vNucFamPtr_it, Xa_relNucFam, nit) {
				xNucFamily &X_nf = *(*nit);

				/* Missing founder is female */
				char *Na_availData;
				if ((*nit)->Xp_mat->N_idx == SAMP_NODATA &&
					(*nit)->Xp_pat->N_idx == SAMP_NODATA)
					Na_availData = NULL;
				else
					Na_availData = X_nf.Xp_mat->N_idx != SAMP_NODATA ?
					Na_data[X_nf.Xp_mat->N_idx] :
				Na_data[X_nf.Xp_pat->N_idx];

				FOREACH (vSampPtr_it, X_nf.Xp_childs, cit) {
					if ((*cit)->B_isMissing) continue;
					if ((*cit)->N_idx == SAMP_NODATA)
						halt_fmt(WISARD_SYST_NULL_CHILDGENO,
							(*cit)->S_FID.c_str(), (*cit)->S_IID.c_str(),
							X_nf.Xp_pat->S_IID.c_str(), X_nf.Xp_mat->S_IID.c_str());
					char *Na_offData = Na_data[(*cit)->N_idx];
					if (!isAvailable(Na_offData[i]))
						continue;

					/* 아빠가 있으면 아빠의 효과를 제거 */
					if (Na_availData)
						R_inferGeno += W2*(wsReal)Na_offData[i]
					- (wsReal)Na_availData[i];
					else /* 아니면 그냥 가져온다 */
						R_inferGeno += (wsReal)Na_offData[i];
					N++;
				}
			}
		}
	} else if (isYChromosome(X_snp)) {
		FOREACH (vNucFamPtr_it, Xa_relNucFam, nit) {
			xNucFamily &X_nf = *(*nit);

			/* In case of chromosome Y, only paternal side will considered */
			if (Xp_misFounder->N_sex == 1) {
				FOREACH (vSampPtr_it, (*nit)->Xp_childs, cit) {
					if ((*cit)->B_isMissing) continue;
					if ((*cit)->N_idx == SAMP_NODATA)
						halt_fmt(WISARD_SYST_NULL_CHILDGENO,
							(*cit)->S_FID.c_str(), (*cit)->S_IID.c_str(),
							X_nf.Xp_pat->S_IID.c_str(), X_nf.Xp_mat->S_IID.c_str());
					char *Na_offData = Na_data[(*cit)->N_idx];
					if (!isAvailable(Na_offData[i]))
						continue;

					/* If this child is female, skip */
					if ((*cit)->N_sex != 1)
						continue;

					/* Y는 그냥 더하면 됨 */
					R_inferGeno += (wsReal)Na_offData[i];
					N++;
				}
			}
		}
	} else {
		/* Autosome */
		FOREACH (vNucFamPtr_it, Xa_relNucFam, nit) {
			char *Na_availData;
			if ((*nit)->Xp_mat->N_idx == SAMP_NODATA &&
				(*nit)->Xp_pat->N_idx == SAMP_NODATA)
				Na_availData = NULL;
			else
				Na_availData = (*nit)->Xp_mat->N_idx != SAMP_NODATA ?
				Na_data[(*nit)->Xp_mat->N_idx] :
				Na_data[(*nit)->Xp_pat->N_idx];
			xNucFamily &X_nf = *(*nit);

			FOREACH (vSampPtr_it, (*nit)->Xp_childs, cit) {
				if ((*cit)->B_isMissing) continue;
				if ((*cit)->N_idx == SAMP_NODATA)
					halt_fmt(WISARD_SYST_NULL_CHILDGENO,
						(*cit)->S_FID.c_str(), (*cit)->S_IID.c_str(),
						X_nf.Xp_pat->S_IID.c_str(), X_nf.Xp_mat->S_IID.c_str());
				char *Na_offData = Na_data[(*cit)->N_idx];
				if (!isAvailable(Na_offData[i]))
					continue;


				if (Na_availData)
					R_inferGeno += W2*(wsReal)Na_offData[i]
						- (wsReal)Na_availData[i];
				else
					R_inferGeno += (wsReal)Na_offData[i]; 
				N++;
			}
		}
	}

	return R_inferGeno;
}

wsReal _inferGeno(xVariant &X_snp, xSample *Xp_misFounder, vNucFamPtr& Xa_relNucFam,
	wsFloat **Ra_data, wsUint &N, wsUint i)
{
	N = 0;
	wsReal R_inferGeno = W0;

	if (isXChromosome(X_snp)) {
		if (Xp_misFounder->N_sex == 1) {
			FOREACH (vNucFamPtr_it, Xa_relNucFam, nit) {
				xNucFamily &X_nf = *(*nit);

				/* Missing founder is male */
				wsFloat *Ra_availData = NULL;
				if ((*nit)->Xp_mat->N_idx != SAMP_NODATA ||
					(*nit)->Xp_pat->N_idx != SAMP_NODATA)
					Ra_availData = X_nf.Xp_mat->N_idx != SAMP_NODATA ?
						Ra_data[X_nf.Xp_mat->N_idx] :
						Ra_data[X_nf.Xp_pat->N_idx];

				FOREACH (vSampPtr_it, (*nit)->Xp_childs, cit) {
					if ((*cit)->B_isMissing) continue;
					if ((*cit)->N_idx == SAMP_NODATA)
						halt_fmt(WISARD_SYST_NULL_CHILDGENO,
							(*cit)->S_FID.c_str(), (*cit)->S_IID.c_str(),
							X_nf.Xp_pat->S_IID.c_str(), X_nf.Xp_mat->S_IID.c_str());

					wsFloat *Ra_offData = Ra_data[(*cit)->N_idx];
					if (!isAvailableReal(Ra_offData[i]))
						continue;
					/* 아빠의 X genotype을 추정할 때 아들은 무시한다 */
					if ((*cit)->N_sex == 1)
						continue;

					/* 엄마가 있으면 엄마의 효과를 제거 */
					if (Ra_availData)
						R_inferGeno += CORR_CONST(2.0)*Ra_offData[i]
							- Ra_availData[i];
					else /* 아니면 1/2만 */
						R_inferGeno += Ra_offData[i]/CORR_CONST(2.0);
					N++;
				}
			}
		} else {
			FOREACH(vNucFamPtr_it, Xa_relNucFam, nit) {
				xNucFamily &X_nf = *(*nit);

				/* Missing founder is female */
				wsFloat *Ra_availData = NULL;
				if ((*nit)->Xp_mat->N_idx != SAMP_NODATA ||
					(*nit)->Xp_pat->N_idx != SAMP_NODATA)
					Ra_availData = X_nf.Xp_mat->N_idx != SAMP_NODATA ?
						Ra_data[X_nf.Xp_mat->N_idx] :
						Ra_data[X_nf.Xp_pat->N_idx];

				FOREACH (vSampPtr_it, X_nf.Xp_childs, cit) {
					if ((*cit)->B_isMissing) continue;
					if ((*cit)->N_idx == SAMP_NODATA)
						halt_fmt(WISARD_SYST_NULL_CHILDGENO,
							(*cit)->S_FID.c_str(), (*cit)->S_IID.c_str(),
							X_nf.Xp_pat->S_IID.c_str(), X_nf.Xp_mat->S_IID.c_str());
					wsFvec Ra_offData = Ra_data[(*cit)->N_idx];
					if (!isAvailableReal(Ra_offData[i]))
						continue;

					/* 아빠가 있으면 아빠의 효과를 제거 */
					if (Ra_availData)
						R_inferGeno += CORR_CONST(2.0)*Ra_offData[i]
							- Ra_availData[i];
					else /* 아니면 그냥 가져온다 */
						R_inferGeno += Ra_offData[i];
					N++;
				}
			}
		}
	} else if (isYChromosome(X_snp)) {
		FOREACH (vNucFamPtr_it, Xa_relNucFam, nit) {
			xNucFamily &X_nf = *(*nit);

			/* In case of chromosome Y, only paternal side will considered */
			if (Xp_misFounder->N_sex == 1) {
				FOREACH (vSampPtr_it, (*nit)->Xp_childs, cit) {
					if ((*cit)->B_isMissing) continue;
					if ((*cit)->N_idx == SAMP_NODATA)
						halt_fmt(WISARD_SYST_NULL_CHILDGENO,
							(*cit)->S_FID.c_str(), (*cit)->S_IID.c_str(),
							X_nf.Xp_pat->S_IID.c_str(), X_nf.Xp_mat->S_IID.c_str());
					wsFvec Ra_offData = Ra_data[(*cit)->N_idx];
					if (!isAvailableReal(Ra_offData[i]))
						continue;

					/* If this child is female, skip */
					if ((*cit)->N_sex != 1)
						continue;

					/* Y는 그냥 더하면 됨 */
					R_inferGeno += Ra_offData[i];
					N++;
				}
			}
		}
	} else {
		/* Autosome */
		FOREACH (vNucFamPtr_it, Xa_relNucFam, nit) {
			wsFvec Ra_availData = NULL;
			if ((*nit)->Xp_mat->N_idx != SAMP_NODATA ||
				(*nit)->Xp_pat->N_idx != SAMP_NODATA)
				Ra_availData = (*nit)->Xp_mat->N_idx != SAMP_NODATA ?
					Ra_data[(*nit)->Xp_mat->N_idx] :
					Ra_data[(*nit)->Xp_pat->N_idx];
			xNucFamily &X_nf = *(*nit);

			FOREACH (vSampPtr_it, (*nit)->Xp_childs, cit) {
				if ((*cit)->B_isMissing) continue;
				if ((*cit)->N_idx == SAMP_NODATA)
					halt_fmt(WISARD_SYST_NULL_CHILDGENO,
						(*cit)->S_FID.c_str(), (*cit)->S_IID.c_str(),
						X_nf.Xp_pat->S_IID.c_str(), X_nf.Xp_mat->S_IID.c_str());
				wsFloat *Ra_offData = Ra_data[(*cit)->N_idx];
				if (!isAvailableReal(Ra_offData[i]))
					continue;

				if (Ra_availData)
					R_inferGeno += CORR_CONST(2.0)*Ra_offData[i]
						- Ra_availData[i];
				else
					R_inferGeno += Ra_offData[i]; 
				N++;
			}
		}
	}

	return R_inferGeno;
}

/**
 * cPPPAnalysisV2::_doImpute Perform imputation of genotypes of given missing founder Xp_misFounder, from the related individuals such as its children or spouse(s)
 *
 * @param     Xp_misFounder			A pointer to specific missing founder to impute its genotype
 * @param     N_dtIdxMisFounderGeno	An index of Ra_misFounderGeno that contains specific missing founder
 * @return    (void)
 */
void cPPPAnalysisV2::_doImpute(xSample *Xp_misFounder, wsUint N_dtIdxMisFounderGeno)
{
	mNucFam		&Xa_nucFam	= Cp_anaFamStr->getNucFamData();
	char		**Na_data	= Cp_IO->getGenotype();
	wsUint		N_SNP		= Cp_IO->sizeVariant();
	vNucFamPtr	Xa_relNucFam;

	if (!Xp_misFounder->B_isMissing)
		halt_fmt(WISARD_SYST_INVL_MISSFOUNDER,
			Xp_misFounder->S_FID.c_str(), Xp_misFounder->S_IID.c_str());

	/* Find matching nuclear family with given missing founder */
	/* 현재 missing founder와 관련이 있는 모든 nuclear family를 찾아낸다 */
	wsUint N_childs = 0;
	FOREACH (mNucFam_it, Xa_nucFam, it) {
		xNucFamily &X_nf = it->second;

		if (X_nf.Xp_mat == Xp_misFounder || X_nf.Xp_pat == Xp_misFounder) {
			Xa_relNucFam.push_back(&X_nf);
			/* # of nonmissing child = # of child - # of missing child */
			N_childs += (wsUint)X_nf.Xp_childs.size() - X_nf.getMissingChildSize();
		}
	}

	/* For all SNP, perform imputation */
	vVariant	Xa_SNP = Cp_IO->getVariant();
	for (wsUint i=0 ; i<N_SNP ; i++) {
		wsReal	R_inferGeno = W0;
		wsUint	N=0;
		
		if (IS_ASSIGNED(dosage))
			R_inferGeno = _inferGeno(Xa_SNP[i], Xp_misFounder, Xa_relNucFam,
				getIO()->getData(), N, i);
		else
			R_inferGeno = _inferGeno(Xa_SNP[i], Xp_misFounder, Xa_relNucFam,
			Na_data, N, i);

		Ra_misFounderGeno[N_dtIdxMisFounderGeno][i] =
			N ? R_inferGeno/(wsReal)N : WISARD_NA_REAL;
	}
	return;
	
}

void cPPPAnalysisV2::_calcPPPcontribution(wsUint N_sIdx, xNucFamily &X_nf,
	wsReal *Rp_misDen, wsReal *Rp_misDiv)
{
//	const char	*Ba_isFounder	= Cp_IO->getIsFounder();
	char		B_misPat		= X_nf.Xp_pat->B_isMissing;
	char		B_misMat		= X_nf.Xp_mat->B_isMissing;
	char		**Na_data		= Cp_IO->getGenotype();

	/* If both parents are available, do nothing */
	if (!B_misPat && !B_misMat)
		return;
	/* If both parents are not founder, do nothing */
	if (X_nf.Xp_pat->Xp_pat || X_nf.Xp_pat->Xp_mat ||
		X_nf.Xp_mat->Xp_pat || X_nf.Xp_mat->Xp_mat)
//	if (!Ba_isFounder[X_nf.Xp_pat->N_idx] || Ba_isFounder[X_nf.Xp_mat->N_idx])
		return;

	/* Get the number of available samples in NF */
	wsUint N_avail = 2 - B_misPat - B_misMat;
	xSample *Xp_avail = !B_misPat ? X_nf.Xp_pat :
		(!B_misMat ? X_nf.Xp_mat : NULL);

	FOREACH (vSampPtr_it, X_nf.Xp_childs, cit) {
		if (!(*cit)->B_isMissing)
			N_avail++;
	}

	/* If there is no avail samples in this NF, do nothing */
	if (N_avail == 0)
		return;
	
	wsReal R_d, R_de;

	if (N_avail%2) {
		R_de = (wsReal)(N_avail>>1);
		R_d = R_de / (R_de+W1);
	} else {
		R_de = (wsReal)N_avail;
		R_d = (R_de-W1)/(R_de+W1);
	}

	/* Add to missingDenom and missingDivid */
	*Rp_misDen += R_d + W1;
	R_d = W1 - R_d;

	/* Add existing parents' genotype if available */
	if (IS_ASSIGNED(dosage)) {
		wsFmat Ra_data = getIO()->getData();
		if (Xp_avail != NULL)
			*Rp_misDiv += R_d * Ra_data[Xp_avail->N_idx][N_sIdx];
		/* Add existing childs' genotype */
		FOREACH (vSampPtr_it, X_nf.Xp_childs, cit) {
			if (!(*cit)->B_isMissing)
				*Rp_misDiv += R_d * Ra_data[(*cit)->N_idx][N_sIdx];
		}
	} else {
		if (Xp_avail != NULL)
			*Rp_misDiv += R_d * Na_data[Xp_avail->N_idx][N_sIdx];
		/* Add existing childs' genotype */
		FOREACH (vSampPtr_it, X_nf.Xp_childs, cit) {
			if (!(*cit)->B_isMissing)
				*Rp_misDiv += R_d * Na_data[(*cit)->N_idx][N_sIdx];
		}
	}
}

int thr_normalize(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	xAnaThread*		X			= (xAnaThread *)Vp_shareData;
	xNormThr*		Xp			= (xNormThr *)X->Vp_data;
	cIO*			Cp_IO		= X->Cp_IO;
	cPPPAnalysisV2*	Cp_anaPPP	= Xp->Cp_anaPPP;
	wsFloat*			Ra_cm		= NULL;
	char*			Na_cm		= NULL;
	Cp_anaPPP->getCorMask(&Ra_cm, &Na_cm, NULL);
	wsUint			N_SNP		= Cp_IO->sizeVariant();
	vSampPtr&		Xa_samp		= Cp_IO->getSample();
	vVariant&		Xa_SNP		= Cp_IO->getVariant();
	int				j			= *((int *)Vp_data);

	if (IS_ASSIGNED(dosage)) {
		/* Dosage case */
		if (!Cp_IO->isDataComplete()) {
			/* Dosage, incomplete */
			_nrmIncomplete(Xp->Ra_data, Xp->Ra_dsg, j, N_SNP, Xa_samp[j]->N_sex,
				Cp_anaPPP->Ra_sum, Cp_anaPPP->Na_ind, Xa_SNP, Cp_anaPPP->Ra_ppp,
				Cp_anaPPP->Ra_ppp_sqrt2pq, Na_cm);
		} else {
			/* Dosage, complete */
			_nrmComplete(Xp->Ra_data, Xp->Ra_dsg, j, N_SNP, Xa_samp[j]->N_sex,
				Cp_anaPPP->Ra_sum, Cp_anaPPP->Na_ind, Xa_SNP, Cp_anaPPP->Ra_ppp,
				Cp_anaPPP->Ra_ppp_sqrt2pq, Na_cm);
		}
	} else {
		if (!Cp_IO->isDataComplete()) {
			/* Genotype, incomplete */
			_nrmIncomplete(Xp->Ra_data, Xp->Na_data, j, N_SNP, Xa_samp[j]->N_sex,
				Cp_anaPPP->Ra_sum, Cp_anaPPP->Na_ind, Xa_SNP, Cp_anaPPP->Ra_ppp,
				Cp_anaPPP->Ra_ppp_sqrt2pq, Na_cm);
		} else {
			/* Genotype, complete */
			_nrmComplete(Xp->Ra_data, Xp->Na_data, j, N_SNP, Xa_samp[j]->N_sex,
				Cp_anaPPP->Ra_sum, Cp_anaPPP->Na_ind, Xa_SNP, Cp_anaPPP->Ra_ppp,
				Cp_anaPPP->Ra_ppp_sqrt2pq, Na_cm);
		}
	}

	return 0;
}

void cPPPAnalysisV2::normalize(wsFloat **Ra_data)
{
	vVariant	Xa_SNP		= Cp_IO->getVariant();
	vSampPtr&	Xa_samp		= Cp_IO->getSample();
	char**		Na_data		= Cp_IO->getGenotype();
	wsUint		N_SNP		= Cp_IO->sizeVariant();
	wsUint		N_sample	= Cp_IO->sizeSample();
	wsUint		i, j;
	sseCalloc(Na_ind, wsUint, N_SNP*3);
	sseCalloc(Ra_sum, wsReal, N_SNP*3);

	if (OPT_NUMBER(thread) > 1) {
		xNormThr	Xp	= { this, Ra_data, Na_data };
		xAnaThread	X	= { getIO(), &Xp };
		WORKER().run(thr_normalize, forAllSample_anaXthr, &X, NULL);
	} else {
		if (IS_ASSIGNED(dosage)) {
			wsFloat**	Ra_dsg = Cp_IO->getDosage();
			/* Dosage case */
			if (!Cp_IO->isDataComplete()) {
				/* Dosage case, incomplete data */
				for (j=0 ; j<N_sample ; j++) {
					_nrmIncomplete(Ra_data, Ra_dsg, j, N_SNP, Xa_samp[j]->N_sex, Ra_sum,
						Na_ind, Xa_SNP, Ra_ppp, Ra_ppp_sqrt2pq, Na_corrMask);
					if ((j%100) == 0)
						notice("%d/%d samples were normalized\r", j, N_sample);
				}
				for (i=0 ; i<N_SNP ; i++)
					Ra_meanX[i] = Ra_sum[i] / (wsReal)Na_ind[i];
			} else {
				/* Dosage case, complete data */
				for (j=0 ; j<N_sample ; j++) {
					_nrmComplete(Ra_data, Ra_dsg, j, N_SNP, Xa_samp[j]->N_sex, Ra_sum,
						Na_ind, Xa_SNP, Ra_ppp, Ra_ppp_sqrt2pq, Na_corrMask);
					if ((j%100) == 0)
						notice("%d/%d samples were normalized\r", j, N_sample);
				}
				for (i=0 ; i<N_SNP ; i++)
					Ra_meanX[i] = Ra_sum[i] / (wsReal)N_SNP;
			}
		} else {
			/* Hard-genotype case */
			if (!Cp_IO->isDataComplete()) {
				/* Hard-genotype case, incomplete data */
				for (j=0 ; j<N_sample ; j++) {
					_nrmIncomplete(Ra_data, Na_data, j, N_SNP,
						Xa_samp[j]->N_sex, Ra_sum, Na_ind, Xa_SNP,
						Ra_ppp, Ra_ppp_sqrt2pq, Na_corrMask);
					if ((j%100) == 0)
						notice("%d/%d samples were normalized\r", j, N_sample);
				}
				for (i=0 ; i<N_SNP ; i++)
					Ra_meanX[i] = Ra_sum[i] / (wsReal)Na_ind[i];
			} else {
				/* Hard-genotype case, complete data */
				for (j=0 ; j<N_sample ; j++) {
					_nrmComplete(Ra_data, Na_data, j, N_SNP,
						Xa_samp[j]->N_sex, Ra_sum, Na_ind, Xa_SNP,
						Ra_ppp, Ra_ppp_sqrt2pq, Na_corrMask);
					if ((j%100) == 0)
						notice("%d/%d samples were normalized\r", j, N_sample);
				}
				for (i=0 ; i<N_SNP ; i++)
					Ra_meanX[i] = Ra_sum[i] / (wsReal)N_SNP;
			}
		}
	}
	sseFree(Na_ind);
	sseFree(Ra_sum);
	LOG("%d/%d samples were normalized\n", N_sample, N_sample);

	/* For dosage, do normalize */
// 	if (IS_ASSIGNED(dosage))
// 		normalize(Ra_data);

	/* Export normalized genotype matrix --verbose */
	if (OPT_ENABLED(makenrm)) {
		wsFloat**	Ra_data		= Cp_IO->getData();
		wsUint		N_sample	= Cp_IO->sizeSample();

		exportMatrix("nrm", Ra_data, N_sample, N_SNP, NULL, NULL, 20);
	}
}

#define USE_SAMPLEWISE_PPP

void _getSampWiseSum(wsReal *Ra_sum, wsUint *Na_ind, char B_isComplete,
	wsUint N_sample, wsUint N_SNP, const char *Ba_isFounder, vSampPtr &Xa_samp,
	char **Na_data)
{
	wsUint i, j;

	if (B_isComplete) {
		wsUint N_founder = 0;
		if (Ba_isFounder) for (i=0 ; i<N_sample ; i++) {
			wsUint N_pIdx = 0;
			if (Ba_isFounder[i] == 0) continue;
			N_founder++;

			switch (Xa_samp[i]->N_sex) {
			case 1: /* male */ N_pIdx = 0; break;
			case 2: /* female */ N_pIdx = 1; break;
			default: /* undet. */ N_pIdx = 2; break;
			}

			for (j=0 ; j<N_SNP ; j++) {
				Ra_sum[j*3+N_pIdx] += (wsReal)Na_data[i][j];
				Na_ind[j*3+N_pIdx]++;
			}
		} else for (i=0 ; i<N_sample ; i++) {
			wsUint N_pIdx = 0;
			N_founder++;

			switch (Xa_samp[i]->N_sex) {
			case 1: /* male */ N_pIdx = 0; break;
			case 2: /* female */ N_pIdx = 1; break;
			default: /* undet. */ N_pIdx = 2; break;
			}

			for (j=0 ; j<N_SNP ; j++) {
				Ra_sum[j*3+N_pIdx] += (wsReal)Na_data[i][j];
				Na_ind[j*3+N_pIdx]++;
			}
		}

		return;
	}

	/* Non-complete case */
	if (Ba_isFounder) for (i=0 ; i<N_sample ; i++) {
		wsUint N_pIdx = 0;
		if (Ba_isFounder[i] == 0) continue;

		switch (Xa_samp[i]->N_sex) {
		case 1: /* male */ N_pIdx = 0; break;
		case 2: /* female */ N_pIdx = 1; break;
		default: /* undet. */ N_pIdx = 2; break;
		}

		for (j=0 ; j<N_SNP ; j++) {
			if (isAvailable(Na_data[i][j])) {
				Ra_sum[j*3+N_pIdx] += (wsReal)Na_data[i][j];
				Na_ind[j*3+N_pIdx]++;
			}
		}
	} else for (i=0 ; i<N_sample ; i++) {
		wsUint N_pIdx = 0;

		switch (Xa_samp[i]->N_sex) {
		case 1: /* male */ N_pIdx = 0; break;
		case 2: /* female */ N_pIdx = 1; break;
		default: /* undet. */ N_pIdx = 2; break;
		}

		for (j=0 ; j<N_SNP ; j++) {
			if (isAvailable(Na_data[i][j])) {
				Ra_sum[j*3+N_pIdx] += (wsReal)Na_data[i][j];
				Na_ind[j*3+N_pIdx]++;
			}
		}
	}
}

void _getSampWiseSum(wsReal *Ra_sum, wsUint *Na_ind, char B_isComplete,
	wsUint N_sample, wsUint N_SNP, const char *Ba_isFounder, vSampPtr &Xa_samp,
	wsFloat **Ra_data)
{
	wsUint i, j;

	if (B_isComplete) {
		wsUint N_founder = 0;
		if (Ba_isFounder) for (i=0 ; i<N_sample ; i++) {
			wsUint N_pIdx = 0;
			if (Ba_isFounder[i] == 0) continue;
			N_founder++;

			switch (Xa_samp[i]->N_sex) {
			case 1: /* male */ N_pIdx = 0;
			case 2: /* female */ N_pIdx = 1;
			default: /* undet. */ N_pIdx = 2;
			}

			for (j=0 ; j<N_SNP ; j++) {
				Ra_sum[j*3+N_pIdx] += Ra_data[i][j];
				Na_ind[j*3+N_pIdx]++;
			}
		} else for (i=0 ; i<N_sample ; i++) {
			wsUint N_pIdx = 0;
			N_founder++;

			switch (Xa_samp[i]->N_sex) {
			case 1: /* male */ N_pIdx = 0;
			case 2: /* female */ N_pIdx = 1;
			default: /* undet. */ N_pIdx = 2;
			}

			for (j=0 ; j<N_SNP ; j++) {
				Ra_sum[j*3+N_pIdx] += Ra_data[i][j];
				Na_ind[j*3+N_pIdx]++;
			}
		}
	}

	/* Non-complete case */
	if (Ba_isFounder) for (i=0 ; i<N_sample ; i++) {
		wsUint N_pIdx = 0;
		if (Ba_isFounder[i] == 0) continue;

		switch (Xa_samp[i]->N_sex) {
		case 1: /* male */ N_pIdx = 0;
		case 2: /* female */ N_pIdx = 1;
		default: /* undet. */ N_pIdx = 2;
		}

		for (j=0 ; j<N_SNP ; j++) {
			if (isAvailable(Ra_data[i][j])) {
				Ra_sum[j*3+N_pIdx] += Ra_data[i][j];
				Na_ind[j*3+N_pIdx]++;
			}
		}
	} else for (i=0 ; i<N_sample ; i++) {
		wsUint N_pIdx = 0;

		switch (Xa_samp[i]->N_sex) {
		case 1: /* male */ N_pIdx = 0;
		case 2: /* female */ N_pIdx = 1;
		default: /* undet. */ N_pIdx = 2;
		}

		for (j=0 ; j<N_SNP ; j++) {
			if (isAvailable(Ra_data[i][j])) {
				Ra_sum[j*3+N_pIdx] += Ra_data[i][j];
				Na_ind[j*3+N_pIdx]++;
			}
		}
	}
}

void cPPPAnalysisV2::_doCalcPPP()
{
	const char
			*Ba_isFounder	= Cp_IO->getIsFounder();
	wsUint	N_SNP			= Cp_IO->sizeVariant();
	wsUint	N_sample		= Cp_IO->sizeSample();
	char	**Na_data		= Cp_IO->getGenotype();
	wsUint	i;
	vVariant	Xa_SNP			= Cp_IO->getVariant();
	vSampPtr&
			Xa_samp			= Cp_IO->getSample();

#ifdef USE_SAMPLEWISE_PPP
	wsUint	*Na_ind			= NULL;
	wsReal	*Ra_sum			= NULL;
	sseCalloc(Na_ind, wsUint, N_SNP*3);
	sseCalloc(Ra_sum, wsReal, N_SNP*3);

	const char *Bp_isFounder = OPT_ENABLED(empiall) ? NULL : Ba_isFounder;
	if (IS_ASSIGNED(dosage))
		_getSampWiseSum(Ra_sum, Na_ind, Cp_IO->isDataComplete(), N_sample,
			N_SNP, Bp_isFounder, Xa_samp, getIO()->getDosage());
	else
		_getSampWiseSum(Ra_sum, Na_ind, Cp_IO->isDataComplete(), N_sample,
			N_SNP, Bp_isFounder, Xa_samp, Na_data);
	LOG("Sample-wise sum calculated\n");
#endif

	for (i=0 ; i<N_SNP ; i++) {
#ifndef USE_SAMPLEWISE_PPP
		wsUint	N[3] = {0, };
		wsReal	Ra_sumSNP[3]	= {W0, };
		wsReal	R_normSum		= W0;
		wsReal	R_mult			= isSexChromosome(Xa_SNP[i]) ?
			W1 : W2;

		/* Sum up count */
		for (j=0 ; j<N_sample ; j++) {
			wsUint N_pIdx = 0;

			switch (Xa_samp[j]->N_sex) {
			case 1: /* male */ N_pIdx = 0;
			case 2: /* female */ N_pIdx = 1;
			default: /* undet. */ N_pIdx = 2;
			}

			if (Ba_isFounder[j] == 1 && isAvailable(Na_data[j][i])) {
				Ra_sumSNP[N_pIdx] += (wsReal)Na_data[j][i];
				N[N_pIdx]++;
			}
		}
#endif
		/* For missing parents */
		wsReal R_missingDenom = W0;
		wsReal R_missingDivid = W0;
		if (OPT_ENABLED(usemf)) {
			mNucFam Xa_nf = Cp_anaFamStr->getNucFamData();
			FOREACH (mNucFam_it, Xa_nf, it) {
				xNucFamily &X_nf = it->second;

				_calcPPPcontribution(i, X_nf, &R_missingDenom, &R_missingDivid);
			}
		}

		/* Sum up pseudo-founder's genotype */
// 		if (OPT_ENABLED(usemf)) {
// 			for (j=0 ; j<Xm_misFounder.size() ; j++) {
// 				if (isAvailableReal(Ra_misFounderGeno[j][i])) {
// 					R_sumSNP += Ra_misFounderGeno[j][i];
// 					N++;
// 				}
// 			}
// 		}

		/* Calculate PPP */
		wsReal R_sum	= W0;
		wsUint N_ind	= 0;
		wsReal R_mul	= 0;
#ifdef USE_SAMPLEWISE_PPP
		wsUint I		= i*3;

		if (isXChromosome(Xa_SNP[i])) {
			/* In case of X-chr
			 * Considering genotypes having its sex */
			R_sum = Ra_sum[I]+Ra_sum[I+1];
			N_ind = Na_ind[I]+(Na_ind[I+1]<<1);
			R_mul = W1;
		} else if (isYChromosome(Xa_SNP[i])) {
			/* In case of Y-chr
			 * Considering genotypes for male */
			R_sum = Ra_sum[I];
			N_ind = Na_ind[I];
			R_mul = W1;
		} else {
			/* In case of autosome */
			R_sum = Ra_sum[I]+Ra_sum[I+1]+Ra_sum[I+2];
			N_ind = (Na_ind[I]+Na_ind[I+1]+Na_ind[I+2])<<1;
			R_mul = W2;
		}

		Ra_ppp[i] = (R_sum + REAL_CONST(0.5) + R_missingDivid) /
			(wsReal)(N_ind + W1 + R_missingDenom);
		Ra_ppp_sqrt2pq[i] = sqrt(R_mul*Ra_ppp[i]*(1-Ra_ppp[i]));
#else
		if (isXChromosome(Xa_SNP[i])) {
			/* In case of X-chr */
			R_sum = Ra_sumSNP[0]+Ra_sumSNP[1];
			N_ind = Na_ind[0]+W2*Na_ind[1];
			R_mul = W1;
		} else if (isYChromosome(Xa_SNP[i])) {
		} else {
			/* In case of autosome */
			R_sum = Ra_sumSNP[0]+Ra_sumSNP[1]+Ra_sumSNP[2];
			N_ind = W2 * (N[0]+N[1]+N[2]);
			R_mul = W2;
		}

		Ra_ppp[i] = (R_sum + REAL_CONST(0.5) + R_missingDivid) /
			(wsReal)(N_ind + W1 + R_missingDenom);
		Ra_ppp_sqrt2pq[i] = sqrt(R_mul*Ra_ppp[i]*(1-Ra_ppp[i]));

		/* Normalize Ra_data */
		if (B_doNorm) {
			if (isXChromosome(Xa_SNP[i])) for (j=N=0 ; j<N_sample ; j++) {
				R_mult = Xa_samp[j]->N_sex==1 ? W1 : W2;
				if (isAvailable(Na_data[j][i])) {
					Ra_data[j][i]	= ((wsReal)Na_data[j][i] - R_mult*Ra_ppp[i])
						/ Ra_ppp_sqrt2pq[i];
					R_normSum		+= (wsReal)Ra_data[j][i];
					N++;
				} else
					Ra_data[j][i]	= WISARD_NA_REAL;
			} else if (isYChromosome(Xa_SNP[i])) for (j=N=0 ; j<N_sample ; j++) {
				R_mult = Xa_samp[j]->N_sex==1 ? W1 : W0;
				if (isAvailable(Na_data[j][i])) {
					Ra_data[j][i]	= ((wsReal)Na_data[j][i] - R_mult*Ra_ppp[i])
						/ Ra_ppp_sqrt2pq[i];
					R_normSum		+= (wsReal)Ra_data[j][i];
					N++;
				} else
					Ra_data[j][i]	= WISARD_NA_REAL;
			} else for (j=N=0 ; j<N_sample ; j++) {
				if (isAvailable(Na_data[j][i])) {
					Ra_data[j][i]	= ((wsReal)Na_data[j][i] - W2*Ra_ppp[i])
						/ Ra_ppp_sqrt2pq[i];
					R_normSum		+= (wsReal)Ra_data[j][i];
					N++;
				} else
					Ra_data[j][i]	= WISARD_NA_REAL;
			}
		}

		/* Get "Normalized E(Xj)", j is idx of SNP */
		Ra_meanX[i] = R_normSum / (wsReal)N;
#endif

		if ((i%1000) == 0)
			notice("%d/%d PPPs were processed\r", i, N_SNP);
	}
	LOG("%d/%d PPPs were processed\n", i, N_SNP);

#ifdef USE_SAMPLEWISE_PPP
//	if (B_doNorm)
//		normalize(Ra_data);

	sseFree(Na_ind);
	sseFree(Ra_sum);
#endif
}

} // End namespace ONETOOL
