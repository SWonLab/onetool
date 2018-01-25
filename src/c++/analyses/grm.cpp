#include "global/common.h"
#ifdef USE_GZ
	#include <zlib.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <limits>
#include "analyses/grm.h"
#include "analyses/pddt.h"
#include "utils/stat.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

inline int corrTest(double R_corr, double R_pddt2, wsUint N_sz, double R_alpha, FILE *H_detOut)
{
	/*1/2 * log (1+p)/log(1-p)  -> N(0.5*log(1+r)/log(1-r),1/N) N <- # of snp, p<-data로부터의 corr
	r = 2*kinship coeff, */
	double tt1 = 0.5f * log((1+R_corr) / (1-R_corr));
	double tt2 = qnorm(R_alpha, 0, 1.0/sqrt((double)(N_sz-3)), 0, 0);
	double LB = 1.0 - 2.0 / (exp(2.0 * (tt1 - tt2)) + 1.0);
	double UB = 1.0 - 2.0 / (exp(2.0 * (tt1 + tt2)) + 1.0);

//	printf("Corr	: %g\n", R_corr);
//	printf("PDDT*2	: %g\n", R_pddt2);
//	printf("Size	: %d\n", N_sz);
//	printf("tt1		: %g\n", tt1);
//	printf("tt2		: %g\n", tt2);
//	printf("LB		: %g\n", LB);
//	printf("UB		: %g\n", UB);

	return (R_pddt2 > UB || R_pddt2 < LB);
}

cGRMAnalysis::cGRMAnalysis(cIO *Cp_inpIO, cPDDTAnalysis *Cp_inpAnaPDDT) : cAnalysis(Cp_inpIO)
{
	if (Cp_inpAnaPDDT == NULL)
		halt("SYSERR: Kinship coefficient instance should exist!");
	Cp_anaPDDT = Cp_inpAnaPDDT;
}

cGRMAnalysis::~cGRMAnalysis()
{

}

void cGRMAnalysis::run()
{
	vSampPtr	&Xa_samples = Cp_IO->getSample();
	wsUint		N_sample	 = Cp_IO->sizeSample();	
	wsUint		N_SNP		 = Cp_IO->sizeVariant();
	wsReal		**Rp_pddt	 = Cp_anaPDDT->getPDDT();
	FILE		*H_grm;

	/* Open input file */
	H_grm = fopen(OPT_STRING(grm), "rb");
	if (H_grm == NULL)
		halt("Can't open GRM matrix input");
	if (N_sample != Xa_samples.size())
		halt("# of sample [%d] is not match with sample information length [%d]", N_sample, Xa_samples.size());

	wsUint N_sig;
	if (fread(&N_sig, sizeof(wsUint), 1, H_grm) != 1)
		halt("Can't read signature of GRM matrix from file!");

	/* Open output */
/**/cExporter* Cp_exporter	= cExporter::summon("comp.grm"); /* CHECKED */
/**/cExporter* Cp_exporter2	= NULL;

	if (OPT_NUMBER(verbose)) {
		Cp_exporter2 = cExporter::summon("detail.grm");/* CHECKED */
		Cp_exporter2->put("IND1	IDX1	IND2	IDX2	CORR	tt1	tt2	"
			"LB	UB	KINSHIP\n");
	}

	wsReal **Ra_haveError = NULL;
	wsAlloc(Ra_haveError, wsReal*, N_sample);
	for (wsUint i=0 ; i<N_sample ; i++) {
		wsAlloc(Ra_haveError[i], wsReal, N_SNP);
		for (wsUint j=0 ; j<N_SNP ; j++)
			Ra_haveError[i][j] = WISARD_NA_REAL;
	}

	/* Not gz file? */
	LOG("GRM alpha %g\n", OPT_REAL(grmalpha));
	char *S_buf = NULL;
	wsAlloc(S_buf, char, 1024);
	LOG("GRM file signature %08X\n", N_sig);
	if (N_sig != 0x00088b1f) {
		rewind(H_grm);
		int		N_idx1=0, N_idx2;

		for (wsUint L=1 ; fgets(S_buf, 1024, H_grm) ; L++) {
			char	*a, *b;

			/* S_buf == Pair index #1 */
			getString(&S_buf, &a);
			N_idx1 = atoi(S_buf);
			/* a == Pair index #2 */
			getString(&a, &b);
			N_idx2 = atoi(a);

			/* Data error */
			if (N_idx1 <= 0 || N_idx2 <= 0 || N_idx1 < N_idx2)
				halt("Data error in line %d on GRM matrix", L);
			if (Xa_samples[N_idx1-1]->S_FID.compare(Xa_samples[N_idx2-1]->S_FID) != 0)
				continue;

			/* Considered N's */
			getString(&b, &a);
			/* a == corr */
			getString(&a, &b);
			wsReal R_grmCorr = (wsReal)atof(a);

			int N_result = 0;
			if (N_idx1 != N_idx2) {
				double	R_corr = R_grmCorr;
				double	R_pddt2 = 2.0 * Rp_pddt[N_idx1-1][N_idx2-1];
				wsUint	N_sz = Cp_IO->sizeVariant();
				double	R_alpha = OPT_REAL(grmalpha);

				/*1/2 * log (1+p)/log(1-p)  -> N(0.5*log(1+r)/log(1-r),1/N) N <- # of snp, p<-data로부터의 corr
				r = 2*kinship coeff, */
				double	tt1	= 0.5 * log((1.0+R_corr) / (1.0-R_corr));
				double	tt2	= qnorm(R_alpha, 0, 1.0/sqrt((double)(N_sz-3)), 0, 0);
				double	LB	= 1.0 - 2.0 / (exp(2.0 * (tt1 - tt2)) + 1.0);
				double	UB	= 1.0 - 2.0 / (exp(2.0 * (tt1 + tt2)) + 1.0);

			//	printf("Corr	: %g\n", R_corr);
			//	printf("PDDT*2	: %g\n", R_pddt2);
			//	printf("Size	: %d\n", N_sz);
			//	printf("tt1		: %g\n", tt1);
			//	printf("tt2		: %g\n", tt2);
			//	printf("LB		: %g\n", LB);
			//	printf("UB		: %g\n", UB);

				N_result = (R_pddt2 > UB || R_pddt2 < LB);

				if (OPT_NUMBER(verbose) && N_result == 1)
// C_exporter2.fmt("IND1 IDX1    IND2    IDX2    CORR    tt1     tt2     LB      UB      KINSHIP\n");
					Cp_exporter2->fmt("%s	%d	%s	%d	%g	%g	%g	%g	%g	%g\n",
						Xa_samples[N_idx1-1]->S_IID.c_str(), N_idx1-1,
						Xa_samples[N_idx2-1]->S_IID.c_str(), N_idx2-1,
						R_corr, tt1, tt2, LB, UB, R_pddt2);
			}

			if (N_result == 1)
				Ra_haveError[N_idx1-1][N_idx2-1] = R_grmCorr;
//				fprintf(H_grm, "%s	%s	%g	%g\n", Xa_samples[N_idx1-1].iid.c_str(), Xa_samples[N_idx2-1].iid.c_str(),
//					R_grmCorr, 2.0 * Rp_pddt[N_idx1-1] [N_idx2-1]);
		}
		/* Data size is not match */
		if ((wsUint)N_idx1 != Cp_IO->sizeSample())
			halt("Input matrix is maybe invalid, size not match");

		fclose(H_grm);
	} else {
		fclose(H_grm);
#ifdef USE_GZ
		wsUint N_idx1 = 0, N_idx2;
		gzFile H_gzGrm = gzopen(OPT_STRING(grm), "r");
		if (H_gzGrm == NULL)
			halt("Can't open GRM matrix input, maybe this is invalid gzip file");

		for (wsUint L=1 ; gzgets(H_gzGrm, S_buf, 1024) ; L++) {
//			LOG("buffer [%s]\n", S_buf);
			char	*a, *b;

			/* S_buf == Pair index #1 */
			getString(&S_buf, &a);
			N_idx1 = atoi(S_buf);
			/* a == Pair index #2 */
			getString(&a, &b);
			N_idx2 = atoi(a);

			/* Data error */
			if (N_idx1 <= 0 || N_idx2 <= 0 || N_idx1 < N_idx2)
				halt("Data error in line %d on GRM matrix", L);
			if (Xa_samples[N_idx1-1]->S_FID.compare(Xa_samples[N_idx2-1]->S_FID) != 0)
				continue;

			/* Considered N's */
			getString(&b, &a);
			/* a == corr */
			getString(&a, &b);
			wsReal R_grmCorr = (wsReal)atof(a);

			int N_result = 0;
			if (N_idx1 != N_idx2) {
				double R_corr = R_grmCorr;
				double R_pddt2 = 2.0 * Rp_pddt[N_idx1-1][N_idx2-1];
				wsUint N_sz = Cp_IO->sizeVariant();
				double R_alpha = OPT_REAL(grmalpha);

				/*1/2 * log (1+p)/log(1-p)  -> N(0.5*log(1+r)/log(1-r),1/N) N <- # of snp, p<-data로부터의 corr
				r = 2*kinship coeff, */
				double tt1 = 0.5 * log((1.0+R_corr) / (1.0-R_corr));
				double tt2 = qnorm(R_alpha, 0, 1.0/sqrt((double)(N_sz-3)), 0, 0);
				double LB = 1.0 - 2.0 / (exp(2.0 * (tt1 - tt2)) + 1.0);
				double UB = 1.0 - 2.0 / (exp(2.0 * (tt1 + tt2)) + 1.0);

			//	printf("Corr	: %g\n", R_corr);
			//	printf("PDDT*2	: %g\n", R_pddt2);
			//	printf("Size	: %d\n", N_sz);
			//	printf("tt1		: %g\n", tt1);
			//	printf("tt2		: %g\n", tt2);
			//	printf("LB		: %g\n", LB);
			//	printf("UB		: %g\n", UB);

				N_result = (R_pddt2 > UB || R_pddt2 < LB);

				if (OPT_NUMBER(verbose) && N_result == 1)
// C_exporter2.fmt("IND1 IDX1    IND2    IDX2    CORR    tt1     tt2     LB      UB      KINSHIP\n");
					Cp_exporter2->fmt("%s   %d      %s      %d      %g      %g      %g      %g      %g      %g\n", Xa_samples[N_idx1-1]->S_IID.c_str(), N_idx1-1, Xa_samples[N_idx2-1]->S_IID.c_str(), N_idx2-1,
                        R_corr, tt1, tt2, LB, UB, R_pddt2);
			}

			if (N_result == 1)
				Ra_haveError[N_idx1-1][N_idx2-1] = R_grmCorr;
//				fprintf(H_grm, "%s	%s	%g	%g\n", Xa_samples[N_idx1-1].iid.c_str(), Xa_samples[N_idx2-1].iid.c_str(),
//					R_grmCorr, 2.0 * Rp_pddt[N_idx1-1][N_idx2-1]);
		}
		/* Data size is not match */
		if (N_idx1 != Cp_IO->sizeSample())
			halt("Input matrix is maybe invalid, size not match");

		gzclose(H_gzGrm);
#else
		halt("GRM matrix is gzipped but this compilation does not support GZ");
#endif
	}
	DEALLOC(S_buf);
	if (Cp_exporter2)
		delete Cp_exporter2;

	Cp_exporter->put("FID	IID	IDX	CONFLICT_WITH	CONF_IDX	RELATION	CORR	KINCOEF\n");
	/* For all family, print */
// 	famMap &X_famMap = Cp_IO->getFamData();
// 	for (famMap_it i=X_famMap.begin() ; i!=X_famMap.end() ; i++) {
// 		sampMap	&X_fam = *i;
// 		char	S_tmp[1024];	
// 
// 		/* For all members of each family */
// 		for (sampMap_it j=X_fam.begin() ; j!=X_fam.end() ; j++) {
// 			SAMPLE X_samp = (*j).second;
// 
// 			sprintf(S_tmp, "%s	%s	%d	", X_samp.fid.c_str(), X_samp.iid.c_str(), X_samp.idx);
// 
// 			/* Show inner confliction */
// 			int K=0;
// 			for (MYuint k=0 ; k<N_sample ; k++) {
// 				if (isAvailable(Ra_haveError[X_samp.idx][k]) && X_samp.idx != (int)k) {
// 					if (K != 0)
// 						fprintf(H_out, "			");
// 					else
// 						fprintf(H_out, S_tmp);
// 					fprintf(H_out, "%s	%d	%s	%g	%g\n", Xa_samples[k].iid.c_str(),
// 						k,
// 						Xa_samples[k].fid.compare(X_samp.fid) == 0 ? "FAMILY" : "UNRELATED",
// 						Ra_haveError[X_samp.idx][k],
// 						Rp_pddt[X_samp.idx][k]);
// 					K++;
// 				}
// 			}
// 
// 			if (K)
// 				fprintf(H_out, "\n");
// 		}
// 	} /* END OF X_famMap*/

	findSampleRemovals(Ra_haveError, 2);

	/* Memory deallocation */
	for (wsUint i=0 ; i<N_sample ; i++) {
		DEALLOC(Ra_haveError[i]);
	}
	DEALLOC(Ra_haveError);
	delete Cp_exporter;
}


int cGRMAnalysis:: sort(const void *p, const void *q)
{
	return (*(int **)q)[1]-(*(int **)p)[1];
}

void cGRMAnalysis::findSampleRemovals(wsReal **Rp_haveError, int N_threshold)
{
	wsUint		N_sample = Cp_IO->sizeSample();
/**/cExporter*	Cp_exporter = cExporter::summon("samples.grm"); /* CHECKED */

	int **Np_result=NULL;
	wsAlloc(Np_result, int*, N_sample);
	
	for (wsUint i=0 ; i<N_sample ; i++) {
		wsAlloc(Np_result[i], int, 2);
		Np_result[i][0] = i;
		Np_result[i][1] = 0;
	}
	// 1. N_sample 개수의 array를 만들고 0으로 초기화

	wsUint	i,j; // for iteration
	int		temp;
	// Rp_haveError : N_sample * N_sample matrix, 각 cell은 MISSING면 edge 없음이며 나머지 경우는 edge 있음임

	for(i=0 ; i<N_sample ; i++) {
		temp=0;
		for(j=0; j<N_sample; j++) {
			if (isAvailableReal(Rp_haveError[i][j]))
				temp++;
		}
		Np_result[i][1]=temp;
	}
	
	// 2. Rp_haveError 모든 cell 탐색하여 array에 사람마다 다른 사람과 연결된 edge 수를 세서 채움
	
	qsort(Np_result,N_sample, sizeof(int **), sort);
	// 3. 내림차순으로 정렬함 (정렬한 후 index가 보존되어 있어야 함)

	vSampPtr	&Xa_samples = Cp_IO->getSample();
	
	for(i=0; i<N_sample; i++) {
		if(Np_result[i][1]>=N_threshold)
			Cp_exporter->fmt("%s	%s	%d\n",
				Xa_samples[Np_result[i][0]]->S_FID.c_str(),
				Xa_samples[Np_result[i][0]]->S_IID.c_str(), Np_result[i][1]);
	}
	// 4. N_threshold 이상의 edge 수를 가지는 모든 individual에 대해 다음 포맷으로 출력 (H_out에)
	//	- 각 줄은 세 개의 항목을 가지며 탭으로 구분됨, 각 줄은 하나의 individual을 나타냄
	//	- 세 개의 항목은 순서대로 FID, IID, edge 수가 되어야 함
	//	
	// *** 여기에 구현하세요 ***
	// 
	delete Cp_exporter;
}

#endif

} // End namespace ONETOOL
