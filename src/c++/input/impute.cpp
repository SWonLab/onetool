#include "utils/matrix.h"
#include "input/impute.h"
#include "global/io.h"
#include "utils/vector.h"
#include "global/worker.h"
#include "global/analysis.h"
#include <math.h>

namespace ONETOOL {

cImpute::cImpute(cIO *Cp_inpIO, char B_useNorm/*=0*/)
{
	this->B_useNorm = B_useNorm;
	wsUint	N_samp	= Cp_inpIO->sizeSample();
	wsUint	N_snp	= Cp_inpIO->sizeVariant();
	char**	Na_data	= Cp_inpIO->getGenotype();
	Cp_IO	= Cp_inpIO;

	/* Prepare data and copy snapshot */
	Ra_data	= sseMatrix(N_snp, N_samp);
	if (B_useNorm) {
		wsFloat **Ra_xdata = Cp_inpIO->getData();
		for (wsUint i=0 ; i<N_samp ; i++) for (wsUint j=0 ; j<N_snp ; j++) {
			wsReal R_geno = Ra_xdata[i][j];
			Ra_data[j][i] = R_geno;
		}
	} else for (wsUint i=0 ; i<N_samp ; i++)
		for (wsUint j=0 ; j<N_snp ; j++) {
			char N_geno = Na_data[i][j];
			Ra_data[j][i] = isMissing(N_geno) ? WISARD_NA_REAL : (wsReal)N_geno;
		}
}

cImpute::~cImpute()
{

}

cMAFimpute::cMAFimpute(cIO *Cp_inpIO) : cImpute(Cp_inpIO)
{
}

void cMAFimpute::impute()
{
	vVariant&	Xa_snp	= Cp_IO->getVariant();
	vSampPtr&	Xa_samp	= Cp_IO->getSample();
	xMaf*		Xp_maf	= Cp_IO->getMAF();
	wsUint		N_snp	= Cp_IO->sizeVariant();

	/* For each SNPs */
	wsUint I=0, J=0;
	FOREACHDO (vVariant_it, Xa_snp, i, I++) {
		xVariant& X_snp = *i;

		FOREACHDO (vSampPtr_it, Xa_samp, j, J++) {
			xSample *Xp_samp = *j;

			if (!isMissingReal(Ra_data[I][J]))
				continue;

			wsReal R_mult = W2;
			if (isXChromosome(X_snp)) {
				switch (Xp_samp->N_sex) {
				case 2: break;
				case 1: R_mult = W1; break;
				default: R_mult = W0; break;
				}
			} else if (isYChromosome(X_snp)) {
				switch (Xp_samp->N_sex) {
				case 1: R_mult = W1; break;
				default: R_mult = W0; break;
				}
			}
			Ra_data[I][J] = R_mult * Xp_maf[I].R_maf;
		}

		notice("%d/%d variants imputed\r", I, N_snp);
	}
	LOG("%d/%d variants imputed\n", N_snp, N_snp);
}

cCorImpute::cCorImpute(cIO *Cp_inpIO, wsSym Ra_cor)
	: cImpute(Cp_inpIO, 0)
{
	M_cor.init(Ra_cor, WISARD_NAN, Cp_IO->sizeSample(), 0, 1);
}

int doCorImpute(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	xAnaThread*	Xp_at		= (xAnaThread *)Vp_shareData;
	cCorImpute*	Cp_ci		= (cCorImpute *)Xp_at->Vp_data;
	cIO*		Cp_IO		= Xp_at->Cp_IO;
	wsMat		Ra_data		= Cp_ci->getData();
	wsSym		Ra_cor		= Cp_ci->getCor();
	wsUint		N_samp		= Cp_IO->sizeSample();
	wsUintCst		N_vrt		= (wsUintCst)Cp_IO->getVariant().size();
	xMaf*		Xp_maf		= Cp_IO->getMAF();
	int*		Np_data		= (int *)Vp_data;
	wsUint		N_s			= (wsUint)Np_data[0];
	wsUint		N_e			= (wsUint)Np_data[1];

	//LOG("Thread %d : Range %d~%d\n", N_idx, N_s, N_e);
	char*	Ba_m	= NULL;
	wsAlloc(Ba_m, char, N_samp);
	wsReal*	Ra_G	= sseVector(N_samp);
	for (wsUint j=N_s ; j<N_e ; j++) {
		xMaf&		X_maf	= Xp_maf[j];
		memset(Ba_m, 0x00, sizeof(char)*N_samp);

		wsUint I=0, J=0;
		LOOP (i, N_samp)
			if (!isMissingReal(Ra_data[j][i])) {
				//				Ra_m2[I]	= Ra_spMean[i];
				Ra_G[I++]	= Ra_data[j][i];
			} else {
				J++;
				//				Ra_m1[J++]	= Ra_spMean[i];
				Ba_m[i]		= 1;
			}
		if (J == 0) continue;

		cVector a(Ra_G, I, 1);
		//		cVector m2(Ra_m2, I, 1);
		//		a += m2;
		a -= W2 * X_maf.R_maf;

		/* Make 12 and 22 */
		cStdMatrix	M_12(J, I);
		cSymMatrix	M_22(I);
		wsMat	Ra_12 = M_12.get();
		wsSym	Ra_22 = M_22.get();
		I = J = 0;
		for (wsUint i=0 ; i<N_samp ; i++) {
			if (Ba_m[i]) /* Go 12 */ {
				if (J >= M_12.row()) halt("STOP2");
				wsUint K=0;
				for (wsUint j=0 ; j<N_samp ; j++)
					if (!Ba_m[j]) Ra_12[J][K++] = Ra_cor[i][j];
				if (K != M_12.col())
					halt("STOP3");
				J++;
				//LOG("J to %d\n", J);
			} else /* Go 22 */ {
				wsUint K=0;
				for (wsUint j=0; K<=I ; j++)
					if (!Ba_m[j]) Ra_22[I][K++] = Ra_cor[i][j];
				if (I != (K-1)) halt("STOP");
				I++;
			}
		}

		if (OPT_ENABLED(kinship)) {
			wsSym Ra_22i = Cp_IO->getFamInv(Ba_m, Ra_cor);
			cSymMatrix M_22i(Ra_22i, N_samp-J);
			M_12 *= M_22i;
		} else {
			cSymMatrix &M_22i = M_22.inv();
			M_12 *= M_22i;
			delete &M_22i;
		}

		cVector Q = M_12*a;
		Q += W2 * X_maf.R_maf;
		//		sseVsV(Ra_m1, Q.get(), Ra_m1, J);
		//		wsReal R_mMean = sseVsum(Ra_m1, J) / (wsReal)J;

		I = 0;
		//		wsReal R_p = W1 - X_mkr.maf;
		// 		wsReal R_t0 = qnorm(SQR(R_p));
		// 		wsReal R_t1 = qnorm(SQR(R_p) + W2*R_p*X_mkr.maf);
		wsReal *Ra_m1 = Q.get();
		for (wsUint i=0 ; i<N_samp ; i++) {
			if (!Ba_m[i]) continue;

			if (fabs(Ra_m1[I]) < 1e-10) Ra_m1[I] = 0;
			Ra_data[j][i] = Ra_m1[I];// / sqrt(Ra_cv[i]);
			//			Ra_m1[I] /= Ra_cv[i];
			// 			if (Ra_m1[I] <= R_t0)
			// 				Ra_data[j][i] = 0;
			// 			else if (Ra_m1[I] <= R_t1)
			// 				Ra_data[j][i] = 1;
			// 			else
			// 				Ra_data[j][i] = 2;
			//			Ra_data[j][i] = (wsUint)(Ra_m1[I]/W2);

			I++;
		}

		if (N_idx == 0 && (j%10) == 0) {
			int S = 0;
			for (int i=0 ; i<OPT_NUMBER(thread) ; i++)
				S += Np_data[i*3+2];
			notice("%d/%d variants processed...\r", S, N_vrt);
		}
		Np_data[2]++;
	}
	sseFree(Ra_G);
	DEALLOC(Ba_m);

	/* Finish normally */
	return 0;
}

void cCorImpute::impute()
{
	wsMat		Ra_cor	= M_cor.get();
	wsUint		N_samp	= Cp_IO->sizeSample();
	wsUint		N_vrt	= Cp_IO->sizeVariant();
	xMaf*		Xp_maf	= Cp_IO->getMAF();

	/* Getting sample-wise mean */
//	wsReal*	Ra_spMean	= sseMmeanCavail(Ra_data, N_snp, N_samp);

	/* Get the cond. dist. of each sample */
//	wsReal* Ra_cv = getCondVar(Cp_IO->getSample(), M_cor);
//	sseVsqrt(Ra_cv, N_samp);

	/* Now we have mean, for all marker, do impute */
// 	wsReal*	Ra_m1	= sseVector(N_samp);
// 	wsReal*	Ra_m2	= sseVector(N_samp);
	if (IS_MULTITHREAD) {
		xAnaThread X_ana = {
			Cp_IO,
			this
		};
		WORKER().run(doCorImpute, forVariant_equal, &X_ana, NULL, sizeof(int)*3);
	} else {
		char*	Ba_m	= NULL;
		wsAlloc(Ba_m, char, N_samp);
		wsReal*	Ra_G	= sseVector(N_samp);
		LOOP (j, N_vrt) {
			xMaf&		X_maf	= Xp_maf[j];
			memset(Ba_m, 0x00, sizeof(char)*N_samp);

			wsUint I=0, J=0;
			LOOP (i, N_samp)
				if (!isMissingReal(Ra_data[j][i])) {
	//				Ra_m2[I]	= Ra_spMean[i];
					Ra_G[I++]	= Ra_data[j][i];
				} else {
					J++;
	//				Ra_m1[J++]	= Ra_spMean[i];
					Ba_m[i]		= 1;
				}
			if (J == 0) continue;

			cVector a(Ra_G, I, 1);
	//		cVector m2(Ra_m2, I, 1);
	//		a += m2;
			a -= W2 * X_maf.R_maf;

			/* Make 12 and 22 */
			cStdMatrix	M_12(J, I);
			cSymMatrix	M_22(I);
			wsMat	Ra_12 = M_12.get();
			wsSym	Ra_22 = M_22.get();
			I = J = 0;
			for (wsUint i=0 ; i<N_samp ; i++) {
				if (Ba_m[i]) /* Go 12 */ {
					if (J >= M_12.row()) halt("STOP2");
					wsUint K=0;
					for (wsUint j=0 ; j<N_samp ; j++)
						if (!Ba_m[j]) Ra_12[J][K++] = Ra_cor[i][j];
					if (K != M_12.col())
						halt("STOP3");
					J++;
					//LOG("J to %d\n", J);
				} else /* Go 22 */ {
					wsUint K=0;
					for (wsUint j=0; K<=I ; j++)
						if (!Ba_m[j]) Ra_22[I][K++] = Ra_cor[i][j];
					if (I != (K-1)) halt("STOP");
					I++;
				}
			}

			if (OPT_ENABLED(kinship)) {
				wsSym Ra_22i = Cp_IO->getFamInv(Ba_m, Ra_cor);
				cSymMatrix M_22i(Ra_22i, N_samp-J);
				M_12 *= M_22i;
			} else {
				cSymMatrix &M_22i = M_22.inv();
				M_12 *= M_22i;
				delete &M_22i;
			}

			cVector Q = M_12*a;
			Q += W2 * X_maf.R_maf;
	//		sseVsV(Ra_m1, Q.get(), Ra_m1, J);
	//		wsReal R_mMean = sseVsum(Ra_m1, J) / (wsReal)J;

			I = 0;
	//		wsReal R_p = W1 - X_mkr.maf;
	// 		wsReal R_t0 = qnorm(SQR(R_p));
	// 		wsReal R_t1 = qnorm(SQR(R_p) + W2*R_p*X_mkr.maf);
			wsReal *Ra_m1 = Q.get();
			for (wsUint i=0 ; i<N_samp ; i++) {
				if (!Ba_m[i]) continue;

				if (fabs(Ra_m1[I]) < 1e-10) Ra_m1[I] = 0;
				Ra_data[j][i] = Ra_m1[I];// / sqrt(Ra_cv[i]);
	//			Ra_m1[I] /= Ra_cv[i];
	// 			if (Ra_m1[I] <= R_t0)
	// 				Ra_data[j][i] = 0;
	// 			else if (Ra_m1[I] <= R_t1)
	// 				Ra_data[j][i] = 1;
	// 			else
	// 				Ra_data[j][i] = 2;
	//			Ra_data[j][i] = (wsUint)(Ra_m1[I]/W2);

				I++;
			}

			notice("%d/%d variants imputed\r", j, N_vrt);
		}
		DEALLOC(Ba_m);
		sseFree(Ra_G);
	}
	/* if B_useNorm */
	if (B_useNorm) {
		/* Substitute to raw genotype if nonmissing */
		char **Na_data = Cp_IO->getGenotype();
		for (wsUint i=0 ; i<N_samp ; i++) for (wsUint j=0 ; j<N_vrt ; j++) {
			char N_geno = Na_data[i][j];
			if (isMissing(N_geno)) continue;
			Ra_data[j][i] = (wsReal)N_geno;
		}
	}
	LOG("%d/%d variants imputed\n", N_vrt, N_vrt);
}

} // End namespace ONETOOL
