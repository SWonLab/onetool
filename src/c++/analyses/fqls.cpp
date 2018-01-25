#include "analyses/fqls.h"
#include "analyses/emai.h"
#include "analyses/pddt.h"
#include "utils/util.h"
#include "utils/dcdflib.h"
#include "utils/matrix.h"
#include "analyses/annot.h"
#include "utils/stat.h"

#define FQLS_V2

namespace ONETOOL {

typedef struct _xFqls
{
	char		B_availOnly;
	vVariant&	Xa_snp;
	xMaf*		Xp_maf;
	wsVec		Ra_adjFullY;
	wsVec		Ra_adjFullY_mqls;
	wsSym		Ra_fullCor;
	wsReal		R_full_1phiInv1_inv;
	wsSym		Ra_full_phiInv;
	wsReal**	Ra_full_Vx;
	wsVec		Ra_full_1phiInv;
	wsVec		Ra_full_Hrow;
	wsReal	*Ra_newI_a, *Ra_newII_a;
	char	B_doFqls;
} xFqls;

typedef struct _xFqlsMulti
{
	char		B_slow;
	vVariant&	Xa_snp;
	xMaf*		Xp_maf;
	cStdMatrix	&M_adjFYfqls, &M_adjFYmqls;
	wsReal	**Ra_fullCor;
	wsReal	R_full_1phiInv1_inv;
	wsReal	**Ra_full_phiInv;
	cSymMatrix	&M_full_Vx;
	wsReal	*Ra_full_1phiInv;
//	wsReal	*Ra_full_Hrow;
	cVector	&V_full_Hrow;
	cStdMatrix	&M_a1, &M_a2;
	//wsReal	*Ra_newI_a, *Ra_newII_a;
} xFqlsMulti;

typedef struct _xFqlsResult
{
	wsUint	N_imp;
	wsReal	R_TfqlsI, R_TfqlsII;
	wsReal	R_Tmqls;
	wsReal	R_time;
} xFqlsResult;

int fetchProband(xSample *Xp_s, void *Vp_data);
resPAfmlar* _getProbandStat(cAnalysis* Cp_anaCorr, wsUint* N_idxProband,
	wsUint N_pb, wsVecCst Ra_fullY, wsReal R_h2, wsSym Ra_cor/*=NULL*/);

/* !!! FQLS-specific analysis !!! */
#if (TOOLSET_TYPE & FUNCTION_MFQLS) || ((TOOLSET_TYPE & 0x100) == 0x100)

wsReal getStat(wsReal *Ra_Hrow, wsUint N_sz, wsReal *Ra_X, wsReal *Ra_a,
	wsReal **Ra_phiInv, wsReal R_1phiInv1_inv, wsReal **Ra_Vx)
{
	/* Get \Sigma_(i=1)^n(a_i) */
	wsReal R_aSum = sseVsum(Ra_a, N_sz);

	/* Get p = (1'��^-11)^-11'��^-1X */
	wsReal R_p = sseVV(Ra_Hrow, N_sz, Ra_X);
	if (R_p <= W0)
		return WISARD_NAN;

	/* Calculate X^T %*% V^T %*% a
	 * 1*N %*% N*N %*% N*1
	 * Due to V^T = I^T-H^T and all columns of H^T are same,
	 * V[j,i] = -H_j + I(j==i)
	 *                       [ x1(-H1+1)+ x1(-H1)+   ... x1(-H1)+  ]
	 * [x1 ... xn] %*% V^T = [ x2(-H2)+   x2(-H2+1)+     x2(-H2)+  ]
	 *                       [ ....       ...                      ]
	 *                       [ xn(-Hn)    xn(-Hn)        xn(-Hn+1) ]
	 *                       
	 *                     = [ x*(-H*)+x1 x*(-H*)+x2 ... x*(-H*)+xn ]
	 * Let Xh = x*H* = x1H1 + x2H2 + ... + xnHn
	 * Thus X^T %*% V^T %*% a = (-Xh+x1)a1 + (-Xh+x2)a2 + ... + (-Xh+xn)an
	 *                        = \Sigma_(i=1)^n(x_i*a_i) - Xh\Sigma_(i=1)^n(a_i)
	 */

	/* Get \Sigma_(i=1)^n(x_i*a_i) */
	wsReal R_TmX = sseVV(Ra_X, N_sz, Ra_a);

	/* Get final XVa */
	wsReal R_XVa = R_TmX - R_p*R_aSum;

	/* Calculate a'��^-1a == (a'��^-1a)/2p(1-p) for autosome, X(fem. only)
	 *                       (a'��^-1a)/p(1-p)  for Y
	 * How it should be changed if there are mixture of male/female for testing variant of X chr.? */
//	wsReal R_asa = sseVpMpV(Ra_a, Ra_phiInv, N_sz) / R_p;

	/* Calculate (a' %*% Vx %*% a) * pmu * (1-pmu/2) */
	wsReal R_denum = sqrt(sseVpMpV(Ra_a, Ra_Vx, N_sz)* R_p*(W1-R_p/W2));
	if (NA(R_denum) || R_denum == 0)
		return WISARD_NAN;
// 		halt("ERROR");
// 	}

//	LOG("Bunja : %g\nBunmo : %g\n", R_XVa, R_denum);

	/* Get 1s1inv == (1'��^-11)^-1 == 2p(1-p)/(1'��^-11) */
//	wsReal R_1s1_inv = R_p * R_1phiInv1_inv;

	/* Get a'1(1s1inv)1'a == sum(a)(1s1inv)sum(a) */
//	wsReal R_asaDec = R_aSum * R_1s1_inv * R_aSum;

	/* Final statistics */
//	return SQR(R_XVa) / (R_asa - R_asaDec);
	return R_XVa / R_denum;
}

wsReal getStatMulti(cVector &V_Hrow, wsUint N_pheno,
	wsReal *Ra_X, cStdMatrix &M_a, wsReal **Ra_phiInv,
	wsReal R_1phiInv1_inv, cSymMatrix &M_Vx)
{
	wsUint	N_sz		= V_Hrow.size();
	/* Get \Sigma_(i=1)^n(a_i) */
	cVector	V_aSum		= M_a.sumR();
	wsReal	*Ra_aSum	= V_aSum.get();

	/* Get p = (1'��^-11)^-11'��^-1X */
	//wsReal R_p = sseVV(V_Hrow.get(), N_sz, Ra_X);

	/* Calculate X^T %*% V^T %*% a
	 * 1*N %*% N*N %*% N*1
	 * Due to V^T = I^T-H^T and all columns of H^T are same,
	 * V[j,i] = -H_j + I(j==i)
	 *                       [ x1(-H1+1)+ x1(-H1)+   ... x1(-H1)+  ]
	 * [x1 ... xn] %*% V^T = [ x2(-H2)+   x2(-H2+1)+     x2(-H2)+  ]
	 *                       [ ....       ...                      ]
	 *                       [ xn(-Hn)    xn(-Hn)        xn(-Hn+1) ]
	 *                       
	 *                     = [ x*(-H*)+x1 x*(-H*)+x2 ... x*(-H*)+xn ]
	 * Let Xh = x*H* = x1H1 + x2H2 + ... + xnHn
	 * Thus X^T %*% V^T %*% a = (-Xh+x1)a1 + (-Xh+x2)a2 + ... + (-Xh+xn)an
	 *                        = \Sigma_(i=1)^n(x_i*a_i) - Xh\Sigma_(i=1)^n(a_i)
	 */
	wsReal *Ra_XVa = sseVector(N_pheno);
	/* Get Xh */
	wsReal R_Xh = sseVV(Ra_X, N_sz, V_Hrow.get());
	for (wsUint i=0 ; i<N_pheno ; i++) {

		/* Get \Sigma_(i=1)^n(x_i*a_i) */
		wsReal R_xa = sseVV(Ra_X, N_sz, M_a.get()[i]);

		/* Get final XVa */
		Ra_XVa[i] = R_xa - R_Xh*Ra_aSum[i];
	}
	cVector V_XVa(Ra_XVa, N_pheno);

	/* Calculate a'��^-1a == (a'��^-1a)/2p(1-p) for autosome, X(fem. only)
	 *                       (a'��^-1a)/p(1-p)  for Y
	 * How it should be changed if there are mixture of male/female for testing variant of X chr.? */
//	wsReal R_asa = sseVpMpV(Ra_a, Ra_phiInv, N_sz) / R_p;

	/* Calculate (a' %*% Vx %*% a) * pmu * (1-pmu/2) */
	cSymMatrix	Ra_aVa	= M_a.MMt(M_Vx);
	wsReal		Ra_qInv	= W1/M_Vx.sum();
	wsSym		Ra_a11a	= sseMtM(&Ra_aSum, 1, N_pheno);
	sseSpC((wsSymCst)Ra_a11a, Ra_qInv, Ra_a11a, N_pheno);
	cSymMatrix	M_a11a(Ra_a11a, N_pheno);
	Ra_aVa -= M_a11a;
// 	sseSsS(Ra_aVa, Ra_a11a, Ra_aVa, N_sz);
	cSymMatrix&	M_aVai = Ra_aVa.inv();
//	MAT_t	Ra_aVai	= invSymMat(Ra_aVa, N_pheno);
	wsReal R_stat = V_XVa.qf(M_aVai);
	delete &M_aVai;

	return R_stat;//sseVpSpV(Ra_XVa, Ra_aVai, N_pheno);
}

int testFQLS(int N_idx, void *Vp_shareData, void *Vp_data,
   void *Vp_result)
{
	xAnaThread	*Xp_at		= (xAnaThread *)Vp_shareData;
	cIO			*Cp_IO		= Xp_at->Cp_IO;
	wsUint	i, j;
	cTimer		t;
	xFqlsResult	*Xp_resFqls	= (xFqlsResult *)Vp_result;
	wsUint		N_samp		= Cp_IO->sizeSample();
	xFqls*		Xp_F		= (xFqls *)(Xp_at->Vp_data);
	char		B_doFqls	= Xp_F->B_doFqls;
	vSampPtr&	Xa_samp		= Cp_IO->getSample();
	char**		Na_data		= Cp_IO->getGenotype();

	char*		Ba_misGeno	= NULL;
	wsAlloc(Ba_misGeno, char, N_samp);

	int		*Np_data	= (int *)Vp_data;
	wsUint	N_s			= (wsUint)Np_data[0];
	wsUint	N_e			= (wsUint)Np_data[1];
	Np_data[2] = 0;

	/* MQLS and B_availOnly==0 should not be activated concurrently */
// 	if (OPT_ENABLED(mqls) && Xp_F->B_availOnly == 0)
// 		halt_fmt(WISARD_CANT_DO_W_SOMEOPT, "Testing with genotype imputation",
// 			"--mqls");

	wsReal *Ra_aM		= NULL; /* Checked */
	wsReal *Ra_aF1		= NULL;
	wsReal *Ra_aF2		= NULL;
	wsReal *Ra_fullX	= sseVector(N_samp); /* Checked */
	wsReal *Ra_availX	= sseVector(N_samp); /* Checked */	

	for (i=N_s ; i<N_e ; i++) {
		wsReal **Ra_currPhiInv	= NULL; /* Checked */
		wsReal **Ra_currPhi	= NULL;
		wsReal **Ra_currVx		= NULL;

		xFqlsResult&	X_res	= Xp_resFqls[i];
		xVariant&		X_snp	= Xp_F->Xa_snp[i];
		xMaf&			X_maf	= Xp_F->Xp_maf[i];
		wsUint			N_miss	= 0, N_nonMiss = 0;

		t.start();
		/* Set to all non-missing */
		memset(Ba_misGeno, 0x00, sizeof(char)*N_samp);

		/* Due to Mqls or Fqls with available genotype
		 * requires checking which SNP is missing, find it */
		if (OPT_ENABLED(impute)) {
			setGeno(Cp_IO, i, Ra_fullX, NULL);
			memcpy(Ra_availX, Ra_fullX, sizeof(wsReal)*N_samp);
		} else for (j=0 ; j<N_samp ; j++) {
			xSample	*Xp_samp	= Xa_samp[j];
			char N_geno = Na_data[j][i];

			if (isMissing(N_geno)) {
				N_miss++;
				Ba_misGeno[j] = 1;

				/* Impute */
				Ra_fullX[j] = _imputeGeno(X_snp, Xp_samp, X_maf.R_maf);
			} else
				/* Insert to both vector */
				Ra_availX[N_nonMiss++] = Ra_fullX[j] = N_geno;
		}
		/* Ra_availX	: Length N_nonMiss */
		/* Ra_fullX		: Length N_sample  */

		// pmu_num <- t(rep(1,nfam*l))%*%ivxphi%*%X
		// pmu <- as.numeric(pmu_num/hat_phi)
		//wsReal R_Wqls_pmu = sseVV(Ra_rSum, N_nonMissSamp, Ra_fullX) / R_IIphiInvII;

		// res <- Y%*%ivxphi / (pmu*(1-pmu/2))
		/*** CALCULATION OF a VECTOR for W_QLS ***/
	// 		wsReal **Ra_res = sseMMt(&Ra_availY, 1, N_nonMissSamp, Ra_Rinv,
	// 			N_nonMissSamp, N_nonMissSamp);
	// 		Ra_Wqls_a = Ra_res[0];
	// 		DEALLOC(Ra_res);
	// 		sseVpS(Ra_Wqls_a, W1/(R_Wqls_pmu*(W1-R_Wqls_pmu/W2)),
	// 			Ra_Wqls_a, N_nonMissSamp);
		X_res.N_imp = N_miss;

		if (N_miss == 0 || !Xp_F->B_availOnly) {
			/* Set Ra_Mqls_a to equal with Ra_adjFullY */
			sseMalloc(Ra_aM, wsReal, N_samp);
			memcpy(Ra_aM, Xp_F->Ra_adjFullY_mqls, sizeof(wsReal)*N_samp);
		} else if (Xp_F->B_availOnly) {
			/* Resizing full phi matrix */
			wsReal **Ra_availPhi	= sseMatrix(N_nonMiss, N_nonMiss); /* checked */
			wsReal **Ra_missPhi	= sseMatrix(N_miss, N_nonMiss); /* checked */
			wsReal *Ra_curAy		= NULL; /* checked */
			wsReal *Ra_curMy		= NULL; /* checked */

			sseMalloc(Ra_curAy, wsReal, N_nonMiss);
			sseMalloc(Ra_curMy, wsReal, N_miss);

			/* Divide Ra_fullPhi into Ra_availPhi and Ra_missPhi */
			for (wsUint j=0,jA=0,jM=0 ; j<N_samp ; j++) {
				wsUint jI;
				wsReal **Ra_tPhi, *Ra_tY;

				if (Ba_misGeno[j]) {
					/* When this sample's genotype is missing,
					 * this sample should be included to phi_n,m and A_m */
					jI		= jM++;
					Ra_tPhi	= Ra_missPhi;
					Ra_tY	= Ra_curMy;	
				} else {
					/* Otherwise,
						* this sample should be included to phi_n,n and A_n */
					jI		= jA++;
					Ra_tPhi	= Ra_availPhi;
					Ra_tY	= Ra_curAy;
				}

				/* Build Y vector */
				Ra_tY[jI]	= Xp_F->Ra_adjFullY_mqls[j];

				/* Build phi matrix */
				for (wsUint k=0,_k=0 ; k<N_samp ; k++) {
					if (Ba_misGeno[k]) continue;

					Ra_tPhi[jI][_k] = Xp_F->Ra_fullCor[j][k];
					_k++;
				}
			}

			/* Ra_availPhi^-1 %*% Ra_missPhi %*% A_m */
			Ra_currPhiInv	= invSymMat(Ra_availPhi, N_nonMiss); /* checked */
			Ra_currPhi		= Ra_availPhi;

			wsReal **Ra_pp = sseMpMt(Ra_currPhiInv, N_nonMiss, N_nonMiss,
				Ra_missPhi, N_miss, N_nonMiss); /* checked */
			sseUnmat(Ra_missPhi, N_miss);

			wsReal *Ra_Av = sseMpV(Ra_pp, N_nonMiss, N_miss, Ra_curMy); /* checked */
			sseFree(Ra_curMy);
			sseUnmat(Ra_pp, N_nonMiss);
			sseVaV(Ra_curAy, Ra_Av, Ra_curAy, N_nonMiss);
			sseFree(Ra_Av);

			/* Set Ra_Mqls_a */
			/*** CALCULATION OF a VECTOR for M_QLS ***/
			Ra_aM = Ra_curAy;
		}

		/* Get statistics for newI, newII when --impute */
		X_res.R_TfqlsI		= WISARD_NAN;
		X_res.R_TfqlsII	= WISARD_NAN;
		if (B_doFqls && Xp_F->B_availOnly) {
			if (N_nonMiss == N_samp) {
				/* There is no missing */
				Ra_aF1	= Xp_F->Ra_newI_a;
				Ra_aF2	= Xp_F->Ra_newII_a;
			} else {
				/* Some genotype missing in there */
				sseMalloc(Ra_aF1, wsReal, N_nonMiss);
				sseMalloc(Ra_aF2, wsReal, N_nonMiss);

				/* Build reduced FqlsI and FqlsII */
				for (wsUint j=0,J=0 ; j<N_samp ; j++) {
					if (Ba_misGeno[j]) continue;

					Ra_aF1[J]	= Xp_F->Ra_newI_a[j];
					Ra_aF2[J]	= Xp_F->Ra_newII_a[j];
					J++;
				}
			}
		} else if (B_doFqls) {
			/* Get statistics and p-value directly */
			X_res.R_TfqlsI		= getStat(Xp_F->Ra_full_Hrow, N_samp,
				Ra_fullX, Xp_F->Ra_newI_a, Xp_F->Ra_full_phiInv,
				Xp_F->R_full_1phiInv1_inv, Xp_F->Ra_full_Vx);
			if (Cp_IO->isProbandAssigned())
				X_res.R_TfqlsII	= getStat(Xp_F->Ra_full_Hrow,
					N_samp, Ra_fullX, Xp_F->Ra_newII_a,
					Xp_F->Ra_full_phiInv, Xp_F->R_full_1phiInv1_inv,
					Xp_F->Ra_full_Vx);
		}
	
		/* Calculate for Mqls */

		/* Get (1^Tphi^-1) */
		wsReal *Ra_Mqls_1phiInv = NULL;
		if (Ra_currPhiInv != NULL)
			Ra_Mqls_1phiInv	= sseMsumR(Ra_currPhiInv, N_nonMiss); /* Checked */
		else {
			Ra_Mqls_1phiInv	= Xp_F->Ra_full_1phiInv;
			Ra_currPhiInv	= Xp_F->Ra_full_phiInv;
			Ra_currPhi		= Xp_F->Ra_fullCor;
		}

		/* Calculate (1^Tphi^-11)^-1 == 1 / sum of 1^Tphi^-1 */
		wsReal R_Mqls_1p1_inv = W1/sseVsum(Ra_Mqls_1phiInv, N_nonMiss);
		wsReal *Ra_Mqls_Hrow = NULL;
		if (!Xp_F->B_availOnly || N_nonMiss==N_samp) {
			Ra_currVx = Xp_F->Ra_full_Vx;
			Ra_Mqls_Hrow = Xp_F->Ra_full_Hrow;
		} else {
			/* Vx = phi * (1'phi^-11)^-1 (MODIFIED into) */
			/* Vx = phi - 1(1^Tphi^-11)^-11^T */
			Ra_currVx = sseMatrix(N_nonMiss, N_nonMiss);
			sseMsC((wsMatCst)Ra_currPhi, R_Mqls_1p1_inv, Ra_currVx, N_nonMiss);
			/* Calculate (1^Tphi^-11)^-1 %*% 1^Tphi^-1 */
			sseVpC(Ra_Mqls_1phiInv, R_Mqls_1p1_inv, Ra_Mqls_1phiInv, N_nonMiss);
			// 	if (N_nonMissSamp != N_sample)
			// 		sseVpS(Ra_avail_1phiInv, R_avail_1phiInv1_inv, Ra_avail_1phiInv, N_nonMissSamp);
			Ra_Mqls_Hrow = Ra_Mqls_1phiInv; /* Checked */
		}


		/* Get statistics for Mqls */
		if (OPT_ENABLED(mqls) || OPT_ENABLED(fastmqls))
			X_res.R_Tmqls = getStat(Ra_Mqls_Hrow, N_nonMiss, Ra_availX,
				Ra_aM, Ra_currPhiInv, R_Mqls_1p1_inv, Ra_currVx);
		sseFree(Ra_aM);

		if (B_doFqls && Xp_F->B_availOnly) {
			/* If missing genotype is not ignored, do calculate now */
			X_res.R_TfqlsI =   getStat(Ra_Mqls_Hrow, Ra_availX?N_nonMiss:N_samp, Ra_availX?Ra_availX:Ra_fullX,
				Ra_aF1, Ra_currPhiInv, R_Mqls_1p1_inv, Ra_currVx);
			if (Cp_IO->isProbandAssigned())
				X_res.R_TfqlsII = getStat(Ra_Mqls_Hrow, N_nonMiss, Ra_availX,
					Ra_aF2, Ra_currPhiInv, R_Mqls_1p1_inv, Ra_currVx);
		}

		X_res.R_time = t.get();

		if (N_samp != N_nonMiss) {
			if (Xp_F->B_availOnly) {
				sseFree(Ra_aF1);
				sseFree(Ra_aF2);
				sseUnmat(Ra_currPhiInv, N_nonMiss);
				sseUnmat(Ra_currPhi, N_nonMiss);
				sseUnmat(Ra_currVx, N_nonMiss);
				sseFree(Ra_Mqls_1phiInv);
			}
		}

		/* Count done */
		Np_data[2]++;

		/* Print progress */
		if (N_idx==0 && (i%10)==0) {
			wsUint N_sum = 0;
			for (int x=0 ; x<OPT_NUMBER(thread) ; x++)
				N_sum += Np_data[x*3 + 2];
			LOG("%d/%d markers tested...\r", N_sum, Xp_F->Xa_snp.size());
		}
	}


	sseFree(Ra_fullX);
	sseFree(Ra_availX);
	return 0;
}

int testFQLSmulti(int N_idx, void *Vp_shareData, void *Vp_data,
   void *Vp_result)
{
	xAnaThread	*Xp_at		= (xAnaThread *)Vp_shareData;
	cIO			*Cp_IO		= Xp_at->Cp_IO;
	wsUint		N_pheno		= Cp_IO->sizePheno();
	wsUint	i, j;
	cTimer		t;
	xFqlsResult	*Xp_resFqls	= (xFqlsResult *)Vp_result;
	wsUint		N_samp		= Cp_IO->sizeSample();
	xFqlsMulti	*Xp_fqls	= (xFqlsMulti *)(Xp_at->Vp_data);
	char		**Na_data	= Cp_IO->getGenotype();

	char		*Ba_misGeno = NULL;
	wsAlloc(Ba_misGeno, char, N_samp);

	int		*Np_data	= (int *)Vp_data;
	wsUint	N_s			= (wsUint)Np_data[0];
	wsUint	N_e			= (wsUint)Np_data[1];
	Np_data[2] = 0;

	/* MQLS and B_availOnly==0 should not be activated concurrently */
	if ((OPT_ENABLED(mqls) || OPT_ENABLED(fastmqls)) && Xp_fqls->B_slow == 0)
		halt_fmt(WISARD_CANT_DO_W_SOMEOPT, "Testing with genotype imputation",
			"--mqls");

	wsReal **Ra_Mqls_a		= NULL; /* Checked */
	wsReal **Ra_FqlsI_a	= NULL;
	wsReal **Ra_FqlsII_a	= NULL;

	wsReal *Ra_fullX	= NULL; /* Checked */
	wsReal *Ra_availX	= NULL; /* Checked */
	sseMalloc(Ra_fullX, wsReal, N_samp);
	sseMalloc(Ra_availX, wsReal, N_samp);

	for (i=N_s ; i<N_e ; i++) {
		wsReal **Ra_currPhiInv	= NULL; /* Checked */
		wsReal **Ra_currPhi	= NULL;
		wsReal **Ra_currVx		= NULL;

		xFqlsResult&	X_res	= Xp_resFqls[i];
		xMaf&			X_maf	= Xp_fqls->Xp_maf[i];
		wsUint			N_miss	= 0, N_nonMiss = 0;
		t.start();
		/* Set to all non-missing */
		memset(Ba_misGeno, 0x00, sizeof(char)*N_samp);

		/* Due to Mqls or Fqls with available genotype
		 * requires checking which SNP is missing, find it */
		if (OPT_ENABLED(impute)) {
			setGeno(Cp_IO, i, Ra_fullX, NULL);
			memcpy(Ra_availX, Ra_fullX, sizeof(wsReal)*N_samp);
		} else for (j=0 ; j<N_samp ; j++) {
			char N_geno = Na_data[j][i];

			if (isMissing(N_geno)) {
				N_miss++;
				Ba_misGeno[j] = 1;

				/* Impute */
				Ra_fullX[j] = X_maf.R_maf * W2;
			} else
				/* Insert to both vector */
				Ra_availX[N_nonMiss++] = Ra_fullX[j] = N_geno;
		}
		/* Ra_availX	: Length N_nonMiss */
		/* Ra_fullX		: Length N_sample  */

		// pmu_num <- t(rep(1,nfam*l))%*%ivxphi%*%X
		// pmu <- as.numeric(pmu_num/hat_phi)
		//wsReal R_Wqls_pmu = sseVV(Ra_rSum, N_nonMissSamp, Ra_fullX) / R_IIphiInvII;

		// res <- Y%*%ivxphi / (pmu*(1-pmu/2))
		/*** CALCULATION OF a VECTOR for W_QLS ***/
	// 		wsReal **Ra_res = sseMMt(&Ra_availY, 1, N_nonMissSamp, Ra_Rinv,
	// 			N_nonMissSamp, N_nonMissSamp);
	// 		Ra_Wqls_a = Ra_res[0];
	// 		DEALLOC(Ra_res);
	// 		sseVpS(Ra_Wqls_a, W1/(R_Wqls_pmu*(W1-R_Wqls_pmu/W2)),
	// 			Ra_Wqls_a, N_nonMissSamp);
		X_res.N_imp = N_miss;

		if (N_miss == 0)
			/* Set Ra_Mqls_a to equal with Ra_adjFullY */
			Ra_Mqls_a = Xp_fqls->M_adjFYmqls.get();
		else if (Xp_fqls->B_slow) {
			/* Resizing full phi matrix */
			wsMat	Ra_availPhi	= sseMatrix(N_nonMiss, N_nonMiss); /* checked */
			wsMat	Ra_missPhi		= sseMatrix(N_miss, N_nonMiss); /* checked */
			wsMat	Ra_curAy		= sseMatrix(N_pheno, N_nonMiss); /* checked */
			wsMat	Ra_curMy		= sseMatrix(N_pheno, N_miss); /* checked */

// 			sseMalloc(Ra_curAy, wsReal, N_nonMiss);
// 			sseMalloc(Ra_curMy, wsReal, N_miss);

			/* Divide Ra_fullPhi into Ra_availPhi and Ra_missPhi */
			for (wsUint j=0,jA=0,jM=0 ; j<N_samp ; j++) {
				wsUint jI;
				wsReal **Ra_tPhi, **Ra_tY;

				if (Ba_misGeno[j]) {
					/* When this sample's genotype is missing,
					 * this sample should be included to phi_n,m and A_m */
					jI		= jM++;
					Ra_tPhi	= Ra_missPhi;
					Ra_tY	= Ra_curMy;	
				} else {
					/* Otherwise,
						* this sample should be included to phi_n,n and A_n */
					jI		= jA++;
					Ra_tPhi	= Ra_availPhi;
					Ra_tY	= Ra_curAy;
				}

				/* Build Y vector */
				for (wsUint k=0 ; k<N_pheno ; k++)
					Ra_tY[k][jI]	= Xp_fqls->M_adjFYmqls.get()[k][j];

				/* Build phi matrix */
				for (wsUint k=0,_k=0 ; k<N_samp ; k++) {
					if (Ba_misGeno[k]) continue;

					Ra_tPhi[jI][_k] = Xp_fqls->Ra_fullCor[j][k];
					_k++;
				}
			}

			/* Ra_availPhi^-1 %*% Ra_missPhi %*% A_m */
			Ra_currPhiInv	= invSymMat(Ra_availPhi, N_nonMiss); /* checked */
			Ra_currPhi		= Ra_availPhi;

			wsReal **Ra_pp = sseMpMt(Ra_currPhiInv, N_nonMiss, N_nonMiss,
				Ra_missPhi, N_miss, N_nonMiss); /* checked */
			sseUnmat(Ra_missPhi, N_miss);

			wsReal **Ra_Av = sseMpMt(Ra_pp, N_nonMiss, N_miss, Ra_curMy,
				N_pheno, N_miss); /* checked */
			sseFree(Ra_curMy);
			sseUnmat(Ra_pp, N_nonMiss);
			sseMaM(Ra_curAy, Ra_Av, Ra_curAy, N_nonMiss, N_samp);
			sseFree(Ra_Av);

			/* Set Ra_Mqls_a */
			/*** CALCULATION OF a VECTOR for M_QLS ***/
			Ra_Mqls_a = Ra_curAy;
		}

		/* Get statistics for newI, newII when --impute */
		X_res.R_TfqlsI		= WISARD_NAN;
		X_res.R_TfqlsII	= WISARD_NAN;
		if (Xp_fqls->B_slow) {
			if (N_nonMiss == N_samp) {
				/* There is no missing */
				Ra_FqlsI_a		= Xp_fqls->M_a1.get();
				Ra_FqlsII_a	= Xp_fqls->M_a2.get();
			} else {
				/* Some genotype missing in there */
				Ra_FqlsI_a		= sseMatrix(N_nonMiss, N_samp);
				Ra_FqlsII_a	= sseMatrix(N_nonMiss, N_samp);

				/* Build reduced FqlsI and FqlsII */
				for (wsUint j=0,J=0 ; j<N_samp ; j++) {
					if (Ba_misGeno[j]) continue;

					for (wsUint k=0 ; k<N_pheno ; k++) {
						Ra_FqlsI_a[k][J]	= Xp_fqls->M_a1.get()[k][j];
						Ra_FqlsII_a[k][J]	= Xp_fqls->M_a2.get()[k][j];
					}
					J++;
				}
			}
		} else {
			/* Get statistics and p-value directly */
			X_res.R_TfqlsI		= getStatMulti(Xp_fqls->V_full_Hrow, N_pheno,
				Ra_fullX, Xp_fqls->M_a1,
				Xp_fqls->Ra_full_phiInv, Xp_fqls->R_full_1phiInv1_inv,
				Xp_fqls->M_full_Vx);
			if (Cp_IO->isProbandAssigned())
				X_res.R_TfqlsII	= getStatMulti(Xp_fqls->V_full_Hrow,
					N_pheno, Ra_fullX, Xp_fqls->M_a2,
					Xp_fqls->Ra_full_phiInv, Xp_fqls->R_full_1phiInv1_inv,
					Xp_fqls->M_full_Vx);
		}
	
		/* Calculate for Mqls */

		/* Get (1^Tphi^-1) */
		wsReal *Ra_Mqls_1phiInv = NULL;
		if (Ra_currPhiInv != NULL)
			Ra_Mqls_1phiInv	= sseMsumR(Ra_currPhiInv, N_nonMiss); /* Checked */
		else {
			Ra_Mqls_1phiInv	= Xp_fqls->Ra_full_1phiInv;
			Ra_currPhiInv	= Xp_fqls->Ra_full_phiInv;
			Ra_currPhi		= Xp_fqls->Ra_fullCor;
		}

		/* Calculate (1^Tphi^-11)^-1 == 1 / sum of 1^Tphi^-1 */
		wsReal R_Mqls_1p1_inv = W1/sseVsum(Ra_Mqls_1phiInv, N_nonMiss);

		cSymMatrix *Mp_currVx = NULL;
		if (!Xp_fqls->B_slow || N_nonMiss==N_samp)
			Mp_currVx = new cSymMatrix(Xp_fqls->M_full_Vx.get(), N_samp, 1);
		else {
			/* Vx = phi * (1'phi^-11)^-1 (MODIFIED into) */
			/* Vx = phi - 1(1^Tphi^-11)^-11^T */
			Ra_currVx = sseMatrix(N_nonMiss, N_nonMiss);
			sseMsC((wsMatCst)Ra_currPhi, R_Mqls_1p1_inv, Ra_currVx, N_nonMiss);
			Mp_currVx = new cSymMatrix(Ra_currVx, N_nonMiss);
		}

		/* Calculate (1^Tphi^-11)^-1 %*% 1^Tphi^-1 */
		sseVpC(Ra_Mqls_1phiInv, R_Mqls_1p1_inv, Ra_Mqls_1phiInv, N_nonMiss);
		// 	if (N_nonMissSamp != N_sample)
		// 		sseVpS(Ra_avail_1phiInv, R_avail_1phiInv1_inv, Ra_avail_1phiInv, N_nonMissSamp);
		wsReal	*Ra_Mqls_Hrow = Ra_Mqls_1phiInv; /* Checked */
		cVector	V_Mqls_Hrow(Ra_Mqls_Hrow, N_nonMiss);

		/* Get statistics for Mqls */
		cStdMatrix M_Mqls_a(N_pheno, N_nonMiss, Ra_Mqls_a);
		cStdMatrix M_FqlsI_a(N_pheno, N_nonMiss, Ra_FqlsI_a);
		cStdMatrix M_FqlsII_a(N_pheno, N_nonMiss, Ra_FqlsII_a);
		if (OPT_ENABLED(mqls) || OPT_ENABLED(fastmqls))
			X_res.R_Tmqls = getStatMulti(V_Mqls_Hrow, N_pheno, Ra_availX,
				M_Mqls_a, Ra_currPhiInv, R_Mqls_1p1_inv, *Mp_currVx);

		if (Xp_fqls->B_slow) {
			/* If missing genotype is not ignored, do calculate now */
			X_res.R_TfqlsI =   getStatMulti(V_Mqls_Hrow, N_pheno, Ra_availX,
				M_FqlsI_a, Ra_currPhiInv, R_Mqls_1p1_inv, *Mp_currVx);
			if (Cp_IO->isProbandAssigned())
				X_res.R_TfqlsII = getStatMulti(V_Mqls_Hrow, N_pheno, Ra_availX,
					M_FqlsII_a, Ra_currPhiInv, R_Mqls_1p1_inv, *Mp_currVx);
		}

		X_res.R_time = t.get();

		if (N_samp != N_nonMiss) {
			if (Xp_fqls->B_slow) {
				sseUnmat(Ra_currPhiInv, N_nonMiss);
				sseUnmat(Ra_currPhi, N_nonMiss);
				sseFree(Ra_Mqls_1phiInv);
			}
		} else {
			M_FqlsI_a.setClean();
			M_FqlsII_a.setClean();
			Mp_currVx->setDontDealloc();
		}

		/* Count done */
		Np_data[2]++;

		/* Print progress */
		if (N_idx==0 && (i%100)==0) {
			wsUint N_sum = 0;
			for (int x=2 ; x<OPT_NUMBER(thread) ; x+=3)
				N_sum += Np_data[x];
			LOG("%d/%d SNPs tested...\r", N_sum, Xp_fqls->Xa_snp.size());
		}
	}


	sseFree(Ra_fullX);
	sseFree(Ra_availX);
	return 0;
}

cStdMatrix cFamQlsAnalysis::_getAdjPheno(cVector& V_prev, cStdMatrix **Mp_adjFYmqls)
{
	wsUintCst		N_pheno		= Cp_IO->sizePheno();
	wsUintCst		N_sample	= Cp_IO->sizeSample();
	wsReal*		Rp_prev		= V_prev.get();
	cStdMatrix	M_FY(N_pheno, N_sample, Cp_IO->getPhenos());
	cStdMatrix	M_adjFYfqls(N_pheno, N_sample);
	*Mp_adjFYmqls = new cStdMatrix(N_pheno, N_sample);

	/* Check # of prevalences */
	if (V_prev.size()!=N_pheno && V_prev.size()!=(N_pheno<<1))
		 halt_fmt(WISARD_INVL_NPREV_NPHENO, V_prev.size(), N_pheno);
	wsUint N_sz = V_prev.size() / N_pheno;

	/* Phenotype adjustment by --prevalence
	 * 
	 * 1. In case of one prevalence
	 *  CTRL -> -P/(1-P)
	 *  CASE -> 1
	 *  MISS -> 0
	 * 
	 * 2. In case of two prevalence
	 *  CTRL -> -P(SEX)/(1-P(SEX)) if SEX is available
	 *       -> 0                  if SEX is missing
	 *  CASE -> 1
	 *  MISS -> 0
	 */
	vSampPtr	&Xa_samp	= Cp_IO->getSample();

	/* Adjust phenotype */
	wsReal *Ra_prev = Rp_prev;
	for (wsUint k=0 ; k<N_pheno ; k++,Ra_prev+=N_sz) {
		wsReal*	Ra_fullY	= M_FY.get()[k];
		wsReal*	Ra_ret		= M_adjFYfqls.get()[k];
		wsReal*	Ra_adjFullY_mqls = (*Mp_adjFYmqls)->get()[k];

		wsUint		i = 0;
		if (Ra_prev) {
			switch (N_sz) {
			case 1: /* Equal to male/female */
				FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
					if (Ra_fullY[i] == WISARD_UNAFFECTED) {
						Ra_ret[i] = W0;
						Ra_adjFullY_mqls[i] = -Ra_prev[0]/(W1-Ra_prev[0]);
					} else if (Ra_fullY[i] == WISARD_AFFECTED)
						Ra_adjFullY_mqls[i] = Ra_ret[i] = W1;
					else {
						Ra_ret[i] = Ra_prev[0];
						Ra_adjFullY_mqls[i] = W0;
					}
				}
				break;
			case 2: /* Non-equal to male/female */
				FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
					if (Ra_fullY[i] == WISARD_UNAFFECTED) {
						Ra_ret[i] = W0;
						switch ((*it)->N_sex) {
						case 1: case 2:
							Ra_adjFullY_mqls[i] = -Ra_prev[(*it)->N_sex-1]/
								(W1-Ra_prev[(*it)->N_sex-1]);
							break;
						default: /* FIXME : How to treat in case of missing sex? */
							halt("Missing sex found, two-prevalence assignment");
						}
					} else if (Ra_fullY[i] == WISARD_AFFECTED) {
						Ra_adjFullY_mqls[i] = Ra_ret[i] = W1;
					} else {
						Ra_adjFullY_mqls[i] = W0;
						switch ((*it)->N_sex) {
						case 1: case 2:
							Ra_ret[i] = Ra_prev[(*it)->N_sex-1];
							break;
						default: /* FIXME : How to treat in case of missing sex? */
							halt("Missing sex found, two-prevalence assignment");
						}
					}
				}
				break;
			}
		} else if (Cp_anaEmAi) {
			wsRealCst *Ra_BLUP = Cp_anaEmAi->getBLUP(i);
			wsRealCst *Ra_pred = Cp_anaEmAi->getPred(i);
			FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
				if (isMissingReal(Ra_fullY[i]))
					Ra_adjFullY_mqls[i] = Ra_ret[i] = W0;
				else
					Ra_adjFullY_mqls[i] = Ra_ret[i] = Ra_fullY[i] - (Ra_BLUP[i] + Ra_pred[i]);
			}
		} else
			halt("SYSERR: No appropriate object found for adjusting phenotype");
	}

	return M_adjFYfqls;
}

cFamQlsAnalysis::cFamQlsAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP,
	cFamStrAnalysis *Cp_inpAnaFam,
	cAnalysis *Cp_inpAnaCorr, cEmAiAnalysisV2* Cp_inpAnaNullEmAi) : cAnalysis(Cp_inpIO)
{
	Cp_anaEmAi = Cp_inpAnaNullEmAi;
	//wsUintCst	N_sample	= Cp_IO->getSampleSize();
	//REAL_c	*Ra_Y		= Cp_IO->getPheno();
	B_doFqls = OPT_ENABLED(fqls) || OPT_ENABLED(fastfqls);

	// Fast version and normal version are m.e.
	if (OPT_ENABLED(mqls) && OPT_ENABLED(fastmqls)) halt_fmt(WISARD_CANT_EXCL_OPT, "--mqls", "--fastmqls");
	if (OPT_ENABLED(fqls) && OPT_ENABLED(fastfqls)) halt_fmt(WISARD_CANT_EXCL_OPT, "--fqls", "--fastfqls");
	if (OPT_ENABLED(mqls) && OPT_ENABLED(fastfqls)) halt_fmt(WISARD_CANT_EXCL_OPT, "--mqls", "--fastfqls");
	if (OPT_ENABLED(fqls) && OPT_ENABLED(fastmqls)) halt_fmt(WISARD_CANT_EXCL_OPT, "--fqls", "--fastmqls");

	// With FQLS, binary phenotype requires both prevalence and heritability
	if (B_doFqls && !Cp_IO->isContinuous() && (!IS_ASSIGNED(prevalence) || !IS_ASSIGNED(heri)))
		halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--fqls with binary phenotype", "requires both --prevalence & --heri");

	/* Option check */
//		halt("Family QLS analysis does not requires BLUP");
	if (Cp_anaEmAi) {
//		if (!Cp_anaEmAi)
	//		halt("SYSERR: EM-AI object is NULL");
		if (Cp_IO->sizePheno() > 1)
			halt_fmt(WISARD_CANT_DO_W_SOMEOPT, "Multiple phenotype",
				"--fqls with continuous phenotype");
		// MQLS cannot be executed with EM-AI since a mixed choice of EM-AI and prevalence is not possible
		if (OPT_ENABLED(mqls) || OPT_ENABLED(fastmqls))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--mqls", "without prevalence");
	} else {
		if (IS_ASSIGNED(blup))
			halt_fmt(WISARD_CANT_FQLS_W_BLUP);
		ASSERT_OPTION(prevalence);

		/* Set R_thr via `prevalence` */
		wsUint N_prev;
		wsReal *Ra_prev = loadRealValues(OPT_STRING(prevalence), &N_prev);
		R_prev = Ra_prev[0];
		/* lower=T */
		R_thr = (wsReal)qnorm(R_prev, 0, 1, 0, 0);
		sseFree(Ra_prev);
	}
	if ((OPT_ENABLED(mqls) || OPT_ENABLED(fastmqls)) && Cp_IO->sizePheno() > 1)
		halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--mqls", "multiple phenotype");
//		halt("Family QLS analysis is only applicable to dichotomous phenotype!");
	// 170315 removed, FQLS with continuous phenotype not require --heri
//	if (B_doFqls)
//		ASSERT_OPTION(heri);
//	ASSERT_OPTION(kinship);

	/* Marking which Y is missing */
	Ba_incSamp = Cp_IO->getPheCovMissing(&N_nonMissSamp);

	/* Set classes */
	Cp_anaFam	= Cp_inpAnaFam;
	Cp_anaCorr	= Cp_inpAnaCorr;
	Cp_anaPPP	= Cp_inpAnaPPP;
}

cFamQlsAnalysis::~cFamQlsAnalysis()
{
	DEALLOC(Ba_incSamp);
}

int cFamQlsAnalysis::_doTest(wsUint i, xVariant &X_snp, wsReal *Rp_TmQLS,
	wsReal *Rp_TfqlsI, wsReal *Rp_TfqlsII, char *Ba_misGeno,
	char **Na_data, wsReal *Ra_adjFullY, wsReal *Ra_adjFullY_mqls,
	wsReal **Ra_fullCor,
	wsReal *Ra_full_Hrow, wsReal *Ra_newI_a, wsReal *Ra_newII_a,
	wsReal **Ra_full_phiInv, wsReal *Ra_full_1phiInv,
	wsReal R_full_1phiInv1_inv, wsReal **Ra_full_Vx,
	char B_slow)
{
#ifdef FQLSvalidate
	char	S_fn[256];
#endif
	wsUint		N_sample	= Cp_IO->sizeSample();
	wsUint		N_miss	= 0, N_nonMiss = 0;
	vSampPtr	&Xa_samp	= Cp_IO->getSample();
	xMaf*		Xp_maf		= Cp_IO->getMAF();

	/* MQLS and B_availOnly==0 should not be activated concurrently */
// 	if (OPT_ENABLED(mqls) && B_availOnly == 0)
// 		halt_fmt(WISARD_CANT_DO_W_SOMEOPT, "Testing with genotype imputation",
// 			"--mqls");
// 		halt("Both of MQLS and genotype imputation cannot performed"
// 			"concurrently");

	wsReal *Ra_aM		= NULL; /* Checked */
	wsReal *Ra_aF1		= NULL;
	wsReal *Ra_aF2		= NULL;
	wsReal *Ra_fullX	= sseVector(N_sample); /* Checked */
	wsReal *Ra_availX	= sseVector(N_sample); /* Checked */


	wsReal **Ra_currPhiInv	= NULL; /* Checked */
	wsReal **Ra_currPhi	= NULL;
	wsReal **Ra_currVx		= NULL;
//		wsReal *Ra_Wqls_a = NULL;

	/* Set to all non-missing */
	memset(Ba_misGeno, 0x00, sizeof(char)*N_sample);

	/* Due to Mqls or Fqls with available genotype
		* requires checking which SNP is missing, find it */
	if (OPT_ENABLED(impute)) {
		setGeno(Cp_IO, i, Ra_fullX, NULL);
		memcpy(Ra_availX, Ra_fullX, sizeof(wsReal)*N_sample);
	} else for (wsUint j=0 ; j<N_sample ; j++) {
		xSample	*Xp_samp	= Xa_samp[j];
		char N_geno = Na_data[j][i];

		if (isMissing(N_geno)) {
			N_miss++;
			Ba_misGeno[j] = 1;

			/* Impute */
			Ra_fullX[j] = _imputeGeno(X_snp, Xp_samp, Xp_maf[i].R_maf);
		} else
			/* Insert to both vector */
			Ra_availX[N_nonMiss++] = Ra_fullX[j] = N_geno;
	}
	/* Ra_availX	: Length N_nonMiss */
	/* Ra_fullX		: Length N_sample  */

	// pmu_num <- t(rep(1,nfam*l))%*%ivxphi%*%X
	// pmu <- as.numeric(pmu_num/hat_phi)
	//wsReal R_Wqls_pmu = sseVV(Ra_rSum, N_nonMissSamp, Ra_fullX) / R_IIphiInvII;

	// res <- Y%*%ivxphi / (pmu*(1-pmu/2))
	/*** CALCULATION OF a VECTOR for W_QLS ***/
// 		wsReal **Ra_res = sseMMt(&Ra_availY, 1, N_nonMissSamp, Ra_Rinv,
// 			N_nonMissSamp, N_nonMissSamp);
// 		Ra_Wqls_a = Ra_res[0];
// 		DEALLOC(Ra_res);
// 		sseVpS(Ra_Wqls_a, W1/(R_Wqls_pmu*(W1-R_Wqls_pmu/W2)),
// 			Ra_Wqls_a, N_nonMissSamp);

	if (N_miss == 0 || !B_slow) {
		/* Set Ra_Mqls_a to equal with Ra_adjFullY */
		sseMalloc(Ra_aM, wsReal, N_sample);
		memcpy(Ra_aM, Ra_adjFullY_mqls, sizeof(wsReal)*N_sample);
	} else if (B_slow) {
		/* Resizing full phi matrix */
		wsReal **Ra_availPhi	= sseMatrix(N_nonMiss, N_nonMiss);	/* checked */
		wsReal **Ra_missPhi		= sseMatrix(N_miss, N_nonMiss);	/* checked */
		wsReal *Ra_curAy		= sseVector(N_nonMiss);				/* checked */
		wsReal *Ra_curMy		= sseVector(N_miss);				/* checked */

		/* Divide Ra_fullPhi into Ra_availPhi and Ra_missPhi */
		for (wsUint j=0,jA=0,jM=0 ; j<N_sample ; j++) {
			wsUint jI;
			wsReal **Ra_tPhi, *Ra_tY;

			if (Ba_misGeno[j]) {
				/* When this sample's genotype is missing,
					* this sample should be included to phi_n,m and A_m */
				jI = jM++;
				Ra_tPhi	= Ra_missPhi;
				Ra_tY	= Ra_curMy;	
			} else {
				/* Otherwise,
					* this sample should be included to phi_n,n and A_n */
				jI = jA++;
				Ra_tPhi	= Ra_availPhi;
				Ra_tY	= Ra_curAy;
			}

			/* Build Y vector */
			Ra_tY[jI] = Ra_adjFullY_mqls[j];

			/* Build phi matrix */
			for (wsUint k=0,_k=0 ; k<N_sample ; k++) {
				if (Ba_misGeno[k]) continue;
				Ra_tPhi[jI][_k] = Ra_fullCor[j][k];
				_k++;
			}
		}

#ifdef FQLSvalidate
		sprintf(S_fn, "%s.fqls.availPhi", X_snp.name);
		exportMatrix(S_fn, Ra_availPhi, N_nonMiss, N_nonMiss);
		sprintf(S_fn, "%s.fqls.missPhi", X_snp.name);
		exportMatrix(S_fn, Ra_missPhi, N_miss, N_nonMiss);
#endif

		/* Ra_availPhi^-1 %*% Ra_missPhi %*% A_m */
		Ra_currPhiInv = invSymMat(Ra_availPhi, N_nonMiss); /* checked */
		Ra_currPhi = Ra_availPhi;

		wsReal **Ra_pp = sseMpMt(Ra_currPhiInv, N_nonMiss, N_nonMiss,
			Ra_missPhi, N_miss, N_nonMiss); /* checked */
		sseUnmat(Ra_missPhi, N_miss);

		wsReal *Ra_Av = sseMpV(Ra_pp, N_nonMiss, N_miss, Ra_curMy); /* checked */
		sseFree(Ra_curMy);
		sseUnmat(Ra_pp, N_nonMiss);
		sseVaV(Ra_curAy, Ra_Av, Ra_curAy, N_nonMiss);
		sseFree(Ra_Av);

		/* Set Ra_Mqls_a */
		/*** CALCULATION OF a VECTOR for M_QLS ***/
		Ra_aM = Ra_curAy;
	}

#ifdef FQLSvalidate
	sprintf(S_fn, "%s.fqls.MqlsA", X_snp.name);
	exportMatrix(S_fn, &Ra_aM, 1, N_nonMiss);
#endif

	/* Get statistics for newI, newII when --impute */
	*Rp_TfqlsI		= WISARD_NAN;
	*Rp_TfqlsII		= WISARD_NAN;
	if (B_doFqls) {
		if (B_slow) {
			if (N_nonMiss == N_sample) {
				Ra_aF1	= Ra_newI_a;
				Ra_aF2	= Ra_newII_a;
			} else {
				sseMalloc(Ra_aF1, wsReal, N_nonMiss);
				sseMalloc(Ra_aF2, wsReal, N_nonMiss);

				/* Build reduced FqlsI and FqlsII */
				for (wsUint j=0,J=0 ; j<N_sample ; j++) {
					if (Ba_misGeno[j]) continue;

					Ra_aF1[J] = Ra_newI_a[j];
					Ra_aF2[J] = Ra_newII_a[j];
					J++;
				}
			}
		} else {
			*Rp_TfqlsI		= getStat(Ra_full_Hrow, N_sample, Ra_fullX,
				Ra_newI_a, Ra_full_phiInv, R_full_1phiInv1_inv, Ra_full_Vx);
			if (Cp_IO->isProbandAssigned())
				*Rp_TfqlsII	= getStat(Ra_full_Hrow, N_sample, Ra_fullX,
					Ra_newII_a, Ra_full_phiInv, R_full_1phiInv1_inv, Ra_full_Vx);
		}
	}
	sseFree(Ra_fullX);
#ifdef FQLSvalidate
	sprintf(S_fn, "%s.fqls.FqlsI", X_snp.name);
	exportMatrix(S_fn, &Ra_aF1, 1, N_nonMiss);
	sprintf(S_fn, "%s.fqls.FqlsII", X_snp.name);
	exportMatrix(S_fn, &Ra_aF2, 1, N_nonMiss);
#endif	
	/* Calculate for Mqls */

	/* Get (1^Tphi^-1) */
	wsReal *Ra_Mqls_1phiInv = NULL;
	if (Ra_currPhiInv != NULL)
		Ra_Mqls_1phiInv	= sseMsumR(Ra_currPhiInv, N_nonMiss); /* Checked */
	else {
		Ra_Mqls_1phiInv	= Ra_full_1phiInv;
		Ra_currPhiInv		= Ra_full_phiInv;
		Ra_currPhi			= Ra_fullCor;
	}

	/* Calculate (1^Tphi^-11)^-1 == 1 / sum of 1^Tphi^-1 */
	wsReal R_Mqls_1p1_inv = W1/sseVsum(Ra_Mqls_1phiInv, N_nonMiss);
	wsReal *Ra_Mqls_Hrow = NULL;
	if (!B_slow || N_nonMiss==N_sample) {
		Ra_currVx = Ra_full_Vx;
		Ra_Mqls_Hrow = Ra_full_Hrow;
	} else {
		/* Vx = phi * (1'phi^-11)^-1 (MODIFIED into) */
		/* Vx = phi - 1(1^Tphi^-11)^-11^T */
		Ra_currVx = sseMatrix(N_nonMiss, N_nonMiss);
		sseMsC((wsMatCst)Ra_currPhi, R_Mqls_1p1_inv, Ra_currVx, N_nonMiss);

		/* Calculate (1^Tphi^-11)^-1 %*% 1^Tphi^-1 */
		sseVpC(Ra_Mqls_1phiInv, R_Mqls_1p1_inv, Ra_Mqls_1phiInv, N_nonMiss);
		// 	if (N_nonMissSamp != N_sample)
		// 		sseVpS(Ra_avail_1phiInv, R_avail_1phiInv1_inv, Ra_avail_1phiInv, N_nonMissSamp);
		Ra_Mqls_Hrow = Ra_Mqls_1phiInv; /* Checked */
	}


	/* Get statistics for Mqls */
	if (OPT_ENABLED(mqls) || OPT_ENABLED(fastmqls))
		*Rp_TmQLS = getStat(Ra_Mqls_Hrow, Ra_availX?N_nonMiss:N_sample,
			Ra_availX?Ra_availX:Ra_fullX,
			Ra_aM, Ra_currPhiInv, R_Mqls_1p1_inv, Ra_currVx);
	sseFree(Ra_aM);

	if (B_doFqls && B_slow) {
		*Rp_TfqlsI = getStat(Ra_Mqls_Hrow, N_nonMiss, Ra_availX,
			Ra_aF1, Ra_currPhiInv, R_Mqls_1p1_inv, Ra_currVx);
		if (Cp_IO->isProbandAssigned())
			*Rp_TfqlsII = getStat(Ra_Mqls_Hrow, N_nonMiss, Ra_availX,
				Ra_aF2, Ra_currPhiInv, R_Mqls_1p1_inv, Ra_currVx);
	}
	sseFree(Ra_availX);
		
	if (N_sample != N_nonMiss) {
		if (B_slow) {
			sseFree(Ra_aF1);
			sseFree(Ra_aF2);
			sseUnmat(Ra_currPhiInv, N_nonMiss);
			sseUnmat(Ra_currPhi, N_nonMiss);
			sseUnmat(Ra_currVx, N_nonMiss);
			sseFree(Ra_Mqls_1phiInv);
		}
	}

	return N_miss;
}

void cFamQlsAnalysis::run()
{
	wsUint		i;
	wsUintCst	N_sample		= Cp_IO->sizeSample();
	/* Perform two-step in default */
	/* Set to 0 in default, do not impute if 1 */
	char	B_slow	= OPT_ENABLED(mqls) && OPT_ENABLED(fqls);
//	char	B_twoStep	= 1;

	/* FIXME: --fastmqls, --fastfqls, --fqls, --mqls mix-up fix */
	if (OPT_ENABLED(avail))
		LOGwarn("--avail option no more effective with MQLS/FQLS\n");

	/* This analysis requires samples */
	if (N_sample == 0) {
		LOG("No sample in analysis, cannot perform\n");
		return;
	}

	if (!B_slow)
		LOGnote("MQLS/FQLS runs with fast mode, variants will be imputed\n");

	/*
	 * Export preparation
	 * 
	 */

	/* Prepare output */
/**/cExporter*	Cp_fqls			= NULL;
	if ((OPT_ENABLED(fqls) && OPT_ENABLED(mqls)) || (OPT_ENABLED(fastfqls) && OPT_ENABLED(fastmqls)))
		Cp_fqls = cExporter::summon(B_slow ? "fqls.mqls.res" : "fast.fqls.mqls.res");
	else if (OPT_ENABLED(fastmqls)) Cp_fqls = cExporter::summon("fast.mqls.res");
	else if (OPT_ENABLED(fastfqls)) Cp_fqls = cExporter::summon("fast.fqls.res");
	else if (OPT_ENABLED(fqls)) Cp_fqls = cExporter::summon("fqls.res");
	else Cp_fqls = cExporter::summon("mqls.res");

	char		S_newIIbuf[256]	= { 0, };
	char		S_mqlsBuf[256]	= { 0, };
	char		S_fqlsBuf[256]	= { 0, };
	char		S_impBuf[256]	= { 0, };
	char		S_time[256]		= { 0, };
	/* Print header */ {
		if (OPT_ENABLED(mqls) || OPT_ENABLED(fastmqls))
			strcpy(S_mqlsBuf, "	MqlsStat	MqlsPval");
		if (OPT_ENABLED(time))
			strcpy(S_time, "	TIME");
		if (OPT_ENABLED(fqls) || OPT_ENABLED(fastfqls)) {
			strcpy(S_fqlsBuf, "	FqlsStat	FqlsPval");
			if (Cp_IO->isProbandAssigned())
				strcpy(S_newIIbuf, "	FqlsStatProband	FqlsPvalProband");
		}
		sprintf(S_impBuf, "	%s", !B_slow?"NIMP":"NMISSING");

		headerVariant(Cp_fqls);
		Cp_fqls->fmt("%s%s%s%s%s\n", S_time,
			S_impBuf, S_mqlsBuf, S_fqlsBuf, S_newIIbuf);
	}

	/*
	 * Data preparation
	 *
	 */
	vVariant		&Xa_snp			= Cp_IO->getVariant();

	/* Prevalence */
	wsUint		N_prev;
	wsReal		*Ra_prev		= loadRealValues(OPT_STRING(prevalence), &N_prev);
	cVector		V_prev(Ra_prev, N_prev);

	/* X */
	char		**Na_data		= Cp_IO->getGenotype();
	/* Correlation matrix */
	wsReal		**Ra_fullCor	= NULL; /* CHECKED */
	wsReal		**Ra_offsetCor	= NULL;
	cPDDTAnalysis
		*Cp_pddt		= NULL;

	/* Get full correlation matrix */
	Ra_offsetCor = Ra_fullCor = getFullCorMat(Cp_anaCorr);

	if (!OPT_ENABLED(kinship) || IS_ASSIGNED(cor)) {
		/* If --fqlsnopddt, force use to PDDT */
		if (!OPT_ENABLED(fqlsnopddt)) {
			/* Calculate PDDT */
			Cp_pddt = new cPDDTAnalysis(getIO());
			Cp_pddt->run();
			Ra_offsetCor = Cp_pddt->getPDDT();
		}
	}

	/* Retest threshold */
	if (!IS_ASSIGNED(retestthr)) {
		OPTION().assign("retestthr", OPTION().getDefVal("retestthr"));
		OPTION().FORCE_OPT_REAL(retestthr);
	}
	wsReal		R_thrRetest		= OPT_REAL(retestthr);

	/*
	 * Data calculation
	 * 
	 */

	/* Inverse of correlation matrix */
	wsReal **Ra_full_pInv	= NULL;
	if (OPT_ENABLED(kinship))
		/* Get family inverse if it is possible to divide corr. mat. into
		 * family-wise form */
		Ra_full_pInv		= Cp_IO->getFamInv(NULL, Ra_fullCor, 0);
	else
		Ra_full_pInv		= invSymMat(Ra_fullCor, N_sample);
#ifdef FQLSvalidate
	exportMatrix("fqls.phiInv", Ra_full_pInv, N_sample, N_sample);
#endif

	/* Calculate 1^Tphi^-1 == rowsum of phi^-1 */
	wsReal *Ra_full_1pInv	= sseMsumR(Ra_full_pInv, N_sample);
	/* Calculate 1^Tphi^-11 == sum of 1^Tphi^-1 */
	wsReal R_full_1pInv1_inv = W1/sseVsum(Ra_full_1pInv, N_sample);

	/* Calculate (1^Tphi^-11)^-1 %*% 1^Tphi^-1 */
	/* H matrix is repeated Ra_*_Hrow n times
	 * Used to get V = I-H */
	sseVpC(Ra_full_1pInv, R_full_1pInv1_inv, Ra_full_1pInv, N_sample);
	wsReal *Ra_full_Hrow = Ra_full_1pInv;
	
	/* Get Vx = phi - 1(1^Tphi^-11)^-11^T */
	wsReal **Ra_full_Vx = sseMatrix(N_sample, N_sample);
	sseMsC((wsMatCst)Ra_fullCor, R_full_1pInv1_inv, Ra_full_Vx, N_sample);

	/* For all SNPs */
	char *Ba_misGeno = NULL;
	sseMalloc(Ba_misGeno, char, N_sample);



	wsUint		N_pheno		= getIO()->sizePheno();
	/* Get --heri */
	wsUint		N_heri		= 0;
	wsReal		*Ra_heri	= NULL;
	if (B_doFqls) {
		Ra_heri = loadRealValues(OPT_STRING(heri), &N_heri);
		if (N_heri != N_pheno)
			halt("Given number of --heri[%d] is not match with phenotypes[%d]",
				N_heri, N_pheno);
//		sseFree(Ra_vec);
	}

	wsReal		*Ra_adjFullY_mqls	= NULL;
	wsReal		*Ra_adjFullY_fqls	= NULL;
	wsReal		*Ra_newI_a			= NULL; /* Checked */
	wsReal		*Ra_newII_a			= NULL; /* Checked */
	cStdMatrix	*Mp_newI_a			= NULL;
	cStdMatrix	*Mp_newII_a			= NULL;
	cStdMatrix	*Mp_adjFYmqls		= NULL;
	cStdMatrix	*Mp_adjFYfqls		= NULL;
	if (N_pheno == 1) {
		/* Unadjusted Y */
		//REAL_c		*Ra_origY		= Cp_IO->getPheno();
		/* Adjusted Y */
		Ra_adjFullY_mqls	= NULL;
		Ra_adjFullY_fqls	= genAdjPheno(getIO(), Cp_anaEmAi, V_prev,
			&Ra_adjFullY_mqls);
		if (B_doFqls) {
			sseCalloc(Ra_newI_a, wsReal, N_sample); /* Default set to 0 */
			sseCalloc(Ra_newII_a, wsReal, N_sample); /* Default set to 0 */

			genFqlsPheno(getIO(), Ra_adjFullY_fqls, Ra_offsetCor,
				&Ra_newI_a, &Ra_newII_a);
		} /* end of if - for fqls */
	} else {
//		*Mp_adjFYmqls = NULL;
		*Mp_adjFYfqls = _getAdjPheno(V_prev, &Mp_adjFYmqls);

		Mp_newI_a	= new cStdMatrix(N_pheno, N_sample);
		Mp_newII_a	= new cStdMatrix(N_pheno, N_sample);
		wsMat Ra_a1 = Mp_newI_a->get();
		wsMat Ra_a2 = Mp_newII_a->get();

		/* For all sample we can calculate each getStat for newI and newII */
		/* Family-wise, newI */
		if (Cp_IO->isContinuous()) {
			memcpy(Ra_a1, Ra_adjFullY_fqls, sizeof(wsReal)*N_sample);
			memcpy(Ra_a2, Ra_adjFullY_fqls, sizeof(wsReal)*N_sample);
		} else {
			mFam &Xa_fam = Cp_IO->getFamilyData();
			FOREACH (mFam_it, Xa_fam, ft) {
				pverbose("Checking family %s\n", ft->first.c_str());
				xFamily &X_fam = ft->second;

				/* Fetch indices of family members */
				vInt X_idxMem;
				X_fam.visit(fetchFamMem, &X_idxMem);
				wsUint *Na_idxMem = NULL; /* Checked */
				wsAlloc(Na_idxMem, wsUint, X_idxMem.size());
				wsUint i=0;
				FOREACHDO (vInt_it, X_idxMem, xt, i++)
					Na_idxMem[i] = *xt;

				/* If data have proband information, newII can be performed */
				for (wsUint P=0 ; P<N_pheno ; P++) {
					if (Cp_IO->isProbandAssigned()) {
						wsReal *R = _buildFamII(getIO(), X_fam, (wsUint)X_idxMem.size(),
							Na_idxMem, Ra_a2[P], Ra_offsetCor, Ra_a1[P],
							getIO()->getPheno(), Ra_heri[P]);
						if (!R)
							halt_fmt(WISARD_SYST_FAIL_GET_FQLS, ft->first.c_str(), "type 2");
					}
					wsReal *R = _buildFamI(getIO(), X_fam, (wsUint)X_idxMem.size(),
						Na_idxMem, Ra_a1[P], Ra_offsetCor, Ra_a1[P],
						getIO()->getPheno(), Ra_heri[P]);
					if (!R)
						halt_fmt(WISARD_SYST_FAIL_GET_FQLS, ft->first.c_str(), "type 1");
				}
				DEALLOC(Na_idxMem);
			}
		} /* END OF if binary */
	}
	sseFree(Ra_heri);
	pverbose("Family traversing end\n");
	if (Cp_pddt) delete Cp_pddt;

#ifdef FQLSvalidate
	exportMatrix("fqls.newI", &Ra_newI_a, 1, N_sample);
	if (Cp_IO->isProbandAssigned())
		exportMatrix("fqls.newII", &Ra_newII_a, 1, N_sample);
#endif

	vInt Xa_idxRetest;
	i = 0;

	char	S_bufAnno[512] = { 0, };
	cTimer	t;
	t.start();
	if (OPT_NUMBER(thread) > 1) {
		xFqlsResult *Xa_res = NULL;
		wsAlloc(Xa_res, xFqlsResult, Xa_snp.size());

		if (N_pheno == 1) {
			xFqls X_runReqs = {
				(char)(B_slow ? 1 : 0),
				Xa_snp,
				getIO()->getMAF(),
				Ra_adjFullY_fqls, Ra_adjFullY_mqls,
				Ra_fullCor, R_full_1pInv1_inv,
				Ra_full_pInv, Ra_full_Vx, Ra_full_1pInv,
				Ra_full_Hrow, Ra_newI_a, Ra_newII_a,
				B_doFqls
			};
			xAnaThread X_at = { getIO(), &X_runReqs };
			WORKER().run(testFQLS, forVariant_equal, &X_at, Xa_res,
				sizeof(int)*3);
		} else {
			cVector V_full_Hrow(Ra_full_Hrow, N_sample);
			cSymMatrix M_full_Vx(Ra_full_Vx, N_sample);
			xFqlsMulti X_runReqs = {
				(char)(B_slow ? 1 : 0),
				Xa_snp,
				getIO()->getMAF(),
				*Mp_adjFYfqls, *Mp_adjFYmqls,
				Ra_fullCor, R_full_1pInv1_inv,
				Ra_full_pInv, M_full_Vx, Ra_full_1pInv,
				V_full_Hrow, *Mp_newI_a, *Mp_newII_a
			};
			xAnaThread X_at = { getIO(), &X_runReqs };
			WORKER().run(testFQLSmulti, forVariant_equal, &X_at, Xa_res,
				sizeof(int)*3);
		}
		/* Export results */
		for (i=0 ; i<Xa_snp.size() ; i++) {
			xFqlsResult&	X_res	= Xa_res[i];
			xVariant		X_snp	= Xa_snp[i];

			/* Export result */
			double	R_PfqlsI	= WISARD_NAN;
			double	R_PfqlsII	= numeric_limits<double>::quiet_NaN();
			char	B_na1		= 0, B_na2 = 1, B_na3 = 1;

			if (B_doFqls) {
				R_PfqlsI	= pnorm(fabs(X_res.R_TfqlsI), 0, 1, 0) * 2.0;
				B_na1		= NA(R_PfqlsI);
				if (Cp_IO->isProbandAssigned()) {
					if (X_res.R_TfqlsII != X_res.R_TfqlsII) {
						sprintf(S_newIIbuf, "	NA	NA");
						R_PfqlsII	= W1;
					} else {
						B_na2 = 0;
						R_PfqlsII	= pnorm(fabs(X_res.R_TfqlsII), 0, 1, 0) * 2.0;
						sprintf(S_newIIbuf, "	%g	%g", X_res.R_TfqlsII, R_PfqlsII);
					}
				} else S_newIIbuf[0] = '\0';

				if (B_na1) strcpy(S_fqlsBuf, "	NA	NA");
				else sprintf(S_fqlsBuf, "	%g	%g", X_res.R_TfqlsI, R_PfqlsI);
			}

			/* Set mQLS result */
			wsReal R_Pmqls = WISARD_NAN;
			if (OPT_ENABLED(mqls) || OPT_ENABLED(fastmqls)) {
				if (X_res.R_Tmqls != X_res.R_Tmqls) {
					sprintf(S_mqlsBuf, "	NA	NA");
				} else {
					B_na3 = 0;
					R_Pmqls = pnorm(fabs(X_res.R_Tmqls), 0, 1, 0) * 2.0;
					sprintf(S_mqlsBuf, "	%g	%g", X_res.R_Tmqls, R_Pmqls);
				}
			}
			else S_mqlsBuf[0] = '\0';

//			if (!B_availOnly)
				sprintf(S_impBuf, "	%d", X_res.N_imp);
//			else S_impBuf[0] = '\0';

			if (IS_ASSIGNED(annogene))
				sprintf(S_bufAnno, "	%s", X_snp.anno);
			if (OPT_ENABLED(time))
				sprintf(S_time, "	%s", cTimer::fmt(X_res.R_time));

			/* Print result */
			if (OPT_ENABLED(remna) && B_na1 && B_na2 && B_na3)
				continue;
			entryVariant(Cp_fqls, X_snp);
			Cp_fqls->fmt("%s%s%s%s%s\n", S_time, S_impBuf,
				S_mqlsBuf, S_fqlsBuf, S_newIIbuf);

			/* If imputation testing was performed and one of FQLS measure
			 * satisfies p-value threshold, put in this index to re-testing
			 * vector */
			if (B_slow==0) if (R_PfqlsI<R_thrRetest ||
				(Cp_IO->isProbandAssigned() && R_PfqlsII<R_thrRetest) ||
				(!B_doFqls && R_Pmqls<R_thrRetest))
				Xa_idxRetest.push_back(i);
		}

		DEALLOC(Xa_res);
	} else FOREACHDO (vVariant_it, Xa_snp, it, i++) {
		wsReal R_TmQLS, R_TfqlsI, R_TfqlsII;
		cTimer t2;

		/* Do family QLS and mQLS test */
		t2.start();
		wsUint N_miss = _doTest(i, *it, &R_TmQLS, &R_TfqlsI, &R_TfqlsII, Ba_misGeno,
			Na_data, Ra_adjFullY_fqls, Ra_adjFullY_mqls,
			Ra_fullCor, Ra_full_Hrow, Ra_newI_a,
			Ra_newII_a, Ra_full_pInv, Ra_full_1pInv,
			R_full_1pInv1_inv, Ra_full_Vx, B_slow ? 1 : 0);

		/* Export result */
		double	R_PfqlsI	= WISARD_NAN;
		char	B_na1		= 0, B_na2 = 1, B_na3 = 1;
		double	R_PfqlsII	= WISARD_NAN;
		if (B_doFqls) {
			R_PfqlsI	= pnorm(fabs(R_TfqlsI), 0, 1, 0) * 2.0;
			B_na1		= NA(R_PfqlsI);
			if (Cp_IO->isProbandAssigned()) {
				if (R_TfqlsII != R_TfqlsII) {
					sprintf(S_newIIbuf, "	NA	NA");
					R_PfqlsII	= W1;
				} else {
					B_na2 = 0;
					R_PfqlsII	= pnorm(fabs(R_TfqlsII), 0, 1, 0) * 2.0;
						sprintf(S_newIIbuf, "	%g	%g", R_TfqlsII, R_PfqlsII);
				}
			} else S_newIIbuf[0] = '\0';

			if (B_na1) strcpy(S_fqlsBuf, "	NA	NA");
			else sprintf(S_fqlsBuf, "	%g	%g", R_TfqlsI, R_PfqlsI);
		}

		/* Set mQLS result */
		wsReal R_Pmqls = WISARD_NAN;
		if (OPT_ENABLED(mqls) || OPT_ENABLED(fastmqls)) {
			if (R_TmQLS != R_TmQLS)
				sprintf(S_mqlsBuf, "	NA	NA");
			else {
				B_na3 = 0;
				R_Pmqls = pnorm(fabs(R_TmQLS), 0, 1, 0) * 2.0;
				sprintf(S_mqlsBuf, "	%g	%g", R_TmQLS, R_Pmqls);
			}
		}
		else S_mqlsBuf[0] = '\0';

//		if (!B_availOnly)
			sprintf(S_impBuf, "	%d", N_miss);
//		else S_impBuf[0] = '\0';

		if (IS_ASSIGNED(annogene))
			sprintf(S_bufAnno, "	%s", it->anno);
		if (OPT_ENABLED(time))
			sprintf(S_time, "	%s", t2.getReadable());

		/* Print result */
		if (OPT_ENABLED(remna) && B_na1 && B_na2 && B_na3)
			continue;
		entryVariant(Cp_fqls, *it);
		Cp_fqls->fmt("%s%s%s%s%s\n", S_time, S_impBuf,
			S_mqlsBuf, S_fqlsBuf, S_newIIbuf);

		/* If imputation testing was performed and one of FQLS measure
		 * satisfies p-value threshold, put in this index to re-testing
		 * vector */
		if (B_slow==0) if (R_PfqlsI<R_thrRetest ||
			(Cp_IO->isProbandAssigned() && R_PfqlsII<R_thrRetest) ||
			(!B_doFqls && R_Pmqls<R_thrRetest))
			Xa_idxRetest.push_back(i);

		if ((i%100) == 0)
			notice("%d variants processed...\r", i);
	}
	LOG("[%s] %d variants successfully tested\n", t.getReadable(),
		Xa_snp.size());

	delete Cp_fqls;
	sseFree(Ba_misGeno);
	sseFree(Ra_newI_a);
	sseFree(Ra_newII_a);
	sseFree(Ra_adjFullY_fqls);
	sseFree(Ra_adjFullY_mqls);
	sseFree(Ra_full_1pInv);
	sseUnmat(Ra_full_pInv, N_sample);
	sseUnmat(Ra_full_Vx, N_sample);
}

#endif

wsVec genAdjPheno(cIO *Cp_IO, cEmAiAnalysisV2* Cp_anaEmAi, cVector& V_prev,
	wsReal **Ra_adjFullY_mqls)
{
	wsUintCst	N_sample	= Cp_IO->sizeSample();
	wsRealCst	*Ra_fullY	= Cp_IO->getPheno();
	wsReal	*Ra_prev	= V_prev.get();
	wsReal	*Ra_ret		= NULL;
	sseMalloc(Ra_ret, wsReal, N_sample);
	if (Ra_adjFullY_mqls)
		sseMalloc(*Ra_adjFullY_mqls, wsReal, N_sample);
	

	/* Phenotype adjustment by --prevalence
	 * 
	 * 1. In case of one prevalence
	 *  CTRL -> -P/(1-P)
	 *  CASE -> 1
	 *  MISS -> 0
	 * 
	 * 2. In case of two prevalence
	 *  CTRL -> -P(SEX)/(1-P(SEX)) if SEX is available
	 *       -> 0                  if SEX is missing
	 *  CASE -> 1
	 *  MISS -> 0
	 */
	vSampPtr	&Xa_samp	= Cp_IO->getSample();
	wsUint		N_prev		= V_prev.size();

	wsUint		i = 0;
	if (N_prev) {
		switch (N_prev) {
		case 1: /* Equal to male/female */
			FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
				if (Ra_fullY[i] == WISARD_UNAFFECTED) {
					Ra_ret[i] = W0;
					if (Ra_adjFullY_mqls)
						(*Ra_adjFullY_mqls)[i] = -Ra_prev[0]/(W1-Ra_prev[0]);
				} else if (Ra_fullY[i] == WISARD_AFFECTED) {
					if (Ra_adjFullY_mqls)
						(*Ra_adjFullY_mqls)[i] = Ra_ret[i] = W1;
					else
						Ra_ret[i] = W1;
				} else {
					Ra_ret[i] = Ra_prev[0];
					if (Ra_adjFullY_mqls)
						(*Ra_adjFullY_mqls)[i] = W0;
				}
			}
			break;
		case 2: /* Non-equal to male/female */
			FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
				if (Ra_fullY[i] == WISARD_UNAFFECTED) {
					Ra_ret[i] = W0;
					if (Ra_adjFullY_mqls) switch ((*it)->N_sex) {
					case 1: case 2:
						(*Ra_adjFullY_mqls)[i] = -Ra_prev[(*it)->N_sex-1]/
							(W1-Ra_prev[(*it)->N_sex-1]);
						break;
					default: /* FIXME : How to treat in case of missing sex? */
						halt("Missing sex found, two-prevalence assignment");
					}
				} else if (Ra_fullY[i] == WISARD_AFFECTED) {
					if (Ra_adjFullY_mqls)
						(*Ra_adjFullY_mqls)[i] = Ra_ret[i] = W1;
					else
						Ra_ret[i] = W1;
				} else {
					if (Ra_adjFullY_mqls)
						(*Ra_adjFullY_mqls)[i] = W0;
					switch ((*it)->N_sex) {
					case 1: case 2:
						Ra_ret[i] = Ra_prev[(*it)->N_sex-1];
						break;
					default: /* FIXME : How to treat in case of missing sex? */
						halt("Missing sex found, two-prevalence assignment");
					}
				}
			}
			break;
		}
	} else if (Cp_anaEmAi) {
		wsRealCst *Ra_BLUP = Cp_anaEmAi->getBLUP(i);
		wsRealCst *Ra_pred = Cp_anaEmAi->getPred(i);
		FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
			if (Ra_adjFullY_mqls) {
				if (isMissingReal(Ra_fullY[i]))
					(*Ra_adjFullY_mqls)[i] = Ra_ret[i] = W0;
				else
					(*Ra_adjFullY_mqls)[i] = Ra_ret[i] = Ra_fullY[i] - (Ra_BLUP[i] + Ra_pred[i]);
			} else {
				if (isMissingReal(Ra_fullY[i]))
					Ra_ret[i] = W0;
				else
					Ra_ret[i] = Ra_fullY[i] - (Ra_BLUP[i] + Ra_pred[i]);
			}
		}
// 		exportVector("blup", Ra_BLUP, N_sample, 5);
// 		exportVector("pred", Ra_pred, N_sample, 5);
	} else
		halt("SYSERR: No appropriate object found for adjusting phenotype");
// 	if (Ra_adjFullY_mqls)
// 		exportVector("newIa", *Ra_adjFullY_mqls, N_sample, 5);
	return Ra_ret;
}

void genFqlsPheno(cIO* Cp_IO, wsVec Ra_adjFullY_fqls, wsSym Ra_offsetCor,
				  wsVec* Ra_newI_a, wsVec* Ra_newII_a)
{
	/* For all sample we can calculate each getStat for newI and newII */
	/* Family-wise, newI */

	if (Cp_IO->isContinuous()) {
		if (Ra_newI_a)
			memcpy(Ra_newI_a[0], Ra_adjFullY_fqls ,sizeof(wsReal)*Cp_IO->sizeSample());
		if (Cp_IO->isProbandAssigned())
			memcpy(Ra_newII_a[0], Ra_adjFullY_fqls ,sizeof(wsReal)*Cp_IO->sizeSample());
	} else {
		mFam &Xa_fam = Cp_IO->getFamilyData();
		FOREACH (mFam_it, Xa_fam, ft) {
			pverbose("Checking family %s\n", ft->first.c_str());
			xFamily &X_fam = ft->second;

			/* Fetch indices of family members */
			vInt X_idxMem;
			X_fam.visit(fetchFamMem, &X_idxMem);
			wsUint *Na_idxMem = NULL; /* Checked */
			wsAlloc(Na_idxMem, wsUint, X_idxMem.size());
			wsUint i=0;
			FOREACHDO (vInt_it, X_idxMem, xt, i++)
				Na_idxMem[i] = *xt;

			/* If data have proband information, newII can be performed */
			if (Cp_IO->isProbandAssigned() && OPT_ENABLED(adjf2)) {
				*Ra_newII_a = _buildFamII(Cp_IO, X_fam, (wsUint)X_idxMem.size(),
					Na_idxMem, Ra_adjFullY_fqls, Ra_offsetCor,
					*Ra_newII_a, Cp_IO->getPheno(), (wsReal)atof(OPT_STRING(heri)));
				if (Ra_newII_a == NULL)
					halt_fmt(WISARD_SYST_FAIL_GET_FQLS, ft->first.c_str(), "type 2");
				//				halt("Family %s error in newII", ft->first.c_str());
			}
			/* Update newI for this family */
			if (Ra_newI_a) {
				*Ra_newI_a = _buildFamI(Cp_IO, X_fam, (wsUint)X_idxMem.size(),
					Na_idxMem, Ra_adjFullY_fqls, Ra_offsetCor, *Ra_newI_a,
					Cp_IO->getPheno(), (wsReal)atof(OPT_STRING(heri)));
				if (Ra_newI_a == NULL)
					halt_fmt(WISARD_SYST_FAIL_GET_FQLS, ft->first.c_str(), "type 1");
			}

			DEALLOC(Na_idxMem);
		}
	}
}

int fetchFamMem(xSample *Xp_s, void *Vp_data)
{
	vInt *Xp_idxFamMem = (vInt *)Vp_data;

	if (Xp_s->N_idx == SAMP_NODATA)
		return 0;
	Xp_idxFamMem->push_back(Xp_s->N_idx);

	/* No meaning what it returns */
	return 0;
}

wsVec _buildFamII(cIO* Cp_IO, xFamily &X_fam, wsUint N_szFam,
	wsUint* Na_idxMem, wsVec Ra_adjY, wsSym Ra_cor, wsVec Ra_newII_a,
	wsVecCst Ra_fullY, wsReal R_h2)
{
	wsUint i;
	/* Set R_thr via `prevalence` */
	wsUint N_prev;
	wsReal *Ra_prev = loadRealValues(OPT_STRING(prevalence), &N_prev);
	wsReal	R_prev = Ra_prev[0];
	/* lower=T */
	wsReal	R_thr = (wsReal)qnorm(R_prev, 0, 1, 0, 0);
	sseFree(Ra_prev);

	/* Fetch the information about proband */
	vInt Xa_idxProband;
	X_fam.visit(fetchProband, &Xa_idxProband);
	wsUint N_proband = (wsUint)Xa_idxProband.size();
	if (N_proband == 0)
		halt("Family [%s] does not have proband information!",
			X_fam.S_FID.c_str());

	wsUint *Na_idxProband = NULL; /* Checked */
	wsAlloc(Na_idxProband, wsUint, N_proband);
	i=0;
	FOREACHDO (vInt_it, Xa_idxProband, xt, i++)
		Na_idxProband[i] = *xt;

	/* Consist non-probands */
	wsUint *Na_idxNonPb = NULL; /* Checked */
	wsUint N_nonPb = N_szFam-N_proband;
	sseMalloc(Na_idxNonPb, wsUint, N_nonPb);

	wsUint N_realNonPb = 0;
	for (i=0 ; i<N_szFam ; i++) {
		char B_found = 0;
		FOREACH (vInt_it, Xa_idxProband, xt)
			if (Na_idxMem[i] == (wsUint)(*xt)) {
				B_found = 1;
				break;
			}
			if (B_found == (char)0)
				Na_idxNonPb[N_realNonPb++] = Na_idxMem[i];
	}
	if (N_realNonPb != N_nonPb)
		return NULL;

	/* Get Vyx */
	// Vxy      <- V[-j,j,drop=F]
	wsReal **Ra_Vyx	= sseMatrix(N_nonPb, N_proband); /* checked */
	wsReal **Ra_Vy	= sseMatrix(N_nonPb, N_nonPb); /* Checked */
	wsReal **Ra_Vx	= sseMatrix(N_proband, N_proband); /* checked */
	if (Ra_cor) for (i=0 ; i<N_proband ; i++) {
		for (wsUint j=0 ; j<N_nonPb ; j++)
			Ra_Vyx[j][i] = Ra_cor[Na_idxNonPb[j]][Na_idxProband[i]]*R_h2;
		for (wsUint j=0 ; j<N_proband ; j++)
			Ra_Vx[i][j] = Ra_cor[Na_idxProband[i]][Na_idxProband[j]]*R_h2;
		Ra_Vx[i][i] += W1;
	} else for (i=0 ; i<N_proband ; i++) {
		for (wsUint j=0 ; j<N_nonPb ; j++) Ra_Vyx[j][i] = W0;
		for (wsUint j=0 ; j<N_proband ; j++)
			Ra_Vx[i][j] = i==j ? R_h2 : W0;
		Ra_Vx[i][i] += W1;
	}
	if (Ra_cor) for (i=0 ; i<N_nonPb ; i++) {
		for (wsUint j=0 ; j<N_nonPb ; j++)
			Ra_Vy[i][j] = Ra_cor[Na_idxNonPb[i]][Na_idxNonPb[j]]*R_h2;
		Ra_Vy[i][i] += W1;
	} else for (i=0 ; i<N_nonPb ; i++) {
		for (wsUint j=0 ; j<N_nonPb ; j++)
			Ra_Vy[i][j] = i==j ? R_h2 : W0;
		Ra_Vy[i][i] += W1;
	}
	wsReal **Ra_invVx = invSymMat(Ra_Vx, N_proband); /* checked */
	sseUnmat(Ra_Vx, N_proband);

	/* Get proband stat */
	// probSum  <- get_mu_var(eProb,T,l-1)
	resPAfmlar *Xp_r = _getProbandStat(getCorrelation(Cp_IO),
		Na_idxProband, N_proband, Ra_fullY, R_h2, Ra_cor); /* checked */

	// mu=Vyx%*%invVx%*%probSum$mu [#nonpb*1]
	wsReal **R_mu = sseMMMt(Ra_Vyx, N_nonPb, N_proband,
		Ra_invVx, N_proband, N_proband,
		&(Xp_r->Ra_muHat), 1, N_proband); /* checked */

	// invVxDec <- invVx%*%probSum$va%*%invVx
	wsReal **Ra_invVx_dec = sseMMMt(Ra_invVx, N_proband,
		N_proband, Xp_r->Ra_varHat, N_proband, N_proband,
		Ra_invVx, N_proband, N_proband); /* checked */
	deallocMatrix(Xp_r->Ra_varHat, N_proband, (void *)1);
	sseFree(Xp_r->Ra_muHat);
	DEALLOC(Xp_r);
	sseMsM(Ra_invVx, Ra_invVx_dec, Ra_invVx, N_proband);
	sseUnmat(Ra_invVx_dec, N_proband);

	// var=eNprob$vars+Vyx%*%(invVx-invVxDec)%*%Vxy
	wsReal **Ra_Vy_P = sseMMMt(Ra_Vyx, N_nonPb, N_proband,
		Ra_invVx, N_proband, N_proband,
		Ra_Vyx, N_nonPb, N_proband); /* checked */
	sseMaM(Ra_Vy, Ra_Vy_P, Ra_Vy, N_nonPb);
	sseUnmat(Ra_Vy_P, N_nonPb);
	sseUnmat(Ra_invVx, N_proband);
	sseUnmat(Ra_Vyx, N_nonPb);

	/* adjY for non-proband */
	for (i=0 ; i<N_nonPb ; i++) {
		// adjY[-o]-pnorm( (T-nprobSum$mu)/sqrt(nprobSum$va[obs]),lower=F )
		if (isMissingReal(Ra_fullY[Na_idxNonPb[i]]))
			Ra_newII_a[Na_idxNonPb[i]] = W0;
		else
			Ra_newII_a[Na_idxNonPb[i]] = Ra_adjY[Na_idxNonPb[i]]
		- (wsReal)(pnorm((R_thr - R_mu[i][0])/sqrt(Ra_Vy[i][i]), 0.0, 1.0, 0));
	}
	sseUnmat(Ra_Vy, N_nonPb);
	sseUnmat(R_mu, N_nonPb);

	/* adjY for proband */
	for (i=0 ; i<N_proband ; i++) {
		// adjY[o]-prev
		if (isMissingReal(Ra_fullY[Na_idxProband[i]]))
			Ra_newII_a[Na_idxProband[i]] = W0;
		else
			Ra_newII_a[Na_idxProband[i]] = Ra_adjY[Na_idxProband[i]]
		- R_prev;
	}

	sseFree(Na_idxNonPb);
	DEALLOC(Na_idxProband);

	return Ra_newII_a;
}

wsVec _buildFamI(cIO* Cp_IO, xFamily &X_fam, wsUint N_szFam,
	wsUint *Na_idxMem, wsVec Ra_adjY, wsSym Ra_cor, wsVec Ra_newI_a,
	wsVecCst Ra_fullY, wsReal R_h2)
{
	wsUint i, j, k;
	
	if (!IS_ASSIGNED(prevalence))
		halt("--prevalence is required");
	
	/* Set R_thr via `prevalence` */
	wsUint N_prev;
	wsReal *Ra_prev = loadRealValues(OPT_STRING(prevalence), &N_prev);
	wsReal	R_prev = Ra_prev[0];
	/* lower=T */
	wsReal	R_thr = (wsReal)qnorm(R_prev, 0, 1, 0, 0);
	sseFree(Ra_prev);
	wsUint N_proband	= N_szFam-1;

	if (N_proband) for (k=0 ; k<N_szFam ; k++) {
		wsUint xt = Na_idxMem[k];

		/* Members except current member will be proband */
		wsUint *Na_idxPb = NULL; /* checked */
		wsUint *Na_currIdxPb = NULL;
		wsAlloc(Na_idxPb, wsUint, N_proband);
		wsAlloc(Na_currIdxPb, wsUint, N_proband);
		memcpy(Na_idxPb, Na_idxMem, sizeof(wsUint)*k);
		if ((N_szFam-k-1) > 0)
			memcpy(Na_idxPb+k, Na_idxMem+k+1, sizeof(wsUint)*(N_szFam-k-1));

		wsUint I, J;
		for (i=I=0 ; i<N_proband ; i++) {
			if (isMissingReal(Ra_fullY[Na_idxPb[i]])) continue;
			Na_currIdxPb[I++] = Na_idxPb[i];
		}
		wsUint N_currProband = I;
		
		if (N_currProband == 0) {
			DEALLOC(Na_idxPb);
			DEALLOC(Na_currIdxPb);

			/* In case of single sample family */
			wsUint		N_prev;
			wsReal*		Ra_prev	= loadRealValues(OPT_STRING(prevalence), &N_prev);
			vSampPtr&	Xa_samp	= Cp_IO->getSample();

			for (k=0 ; k<N_szFam ; k++) {
				wsUint xt = Na_idxMem[k];

				if (N_prev == 1) {
					Ra_newI_a[xt] = Ra_adjY[xt] - Ra_prev[0];
				} else {
					char N_sex = Xa_samp[xt]->N_sex;

					/* If there one individual */
					switch (N_sex) {
					case 1: case 2:
						Ra_newI_a[xt] = Ra_adjY[xt] - Ra_prev[N_sex-1];
						break;
					default:
						halt("Invalid sex value[%d] for two-prevalence",
							N_sex);
					}
				}
			}
			sseFree(Ra_prev);
		} else {
			/* Get Vyx */
			// Vxy      <- V[-j,j,drop=F]
			wsReal *Ra_Vyx	= NULL; /* Checked */
#ifndef FQLS_V2
			wsReal **Ra_Vx	= sseMatrix(N_currProband, N_currProband); /* checked */
			sseMalloc(Ra_Vyx, wsReal, N_currProband);
			for (i=I=0 ; i<N_proband ; i++) {
				if (isMissingReal(Ra_fullY[Na_idxPb[i]])) continue;

				Ra_Vyx[I] = Ra_cor[xt][Na_idxPb[i]]*R_h2;
				for (j=J=0 ; j<N_proband ; j++) {
					if (isMissingReal(Ra_fullY[Na_idxPb[j]])) continue;

					Ra_Vx[I][J] = Ra_cor[Na_idxPb[i]][Na_idxPb[j]] * R_h2;
					J++;
				}
				Ra_Vx[I][I] += W1; //*= R_h2;
				I++;
			}
			wsMat Ra_invVx = invSymMat(Ra_Vx, N_currProband); /* checked */
#else
			cSymMatrix	M_Vx(N_currProband);
			wsSym		Ra_Vx	= M_Vx.get(); /* checked */
			sseMalloc(Ra_Vyx, wsReal, N_currProband);
			if (Ra_cor) for (i=I=0 ; i<N_proband ; i++) {
				if (isMissingReal(Ra_fullY[Na_idxPb[i]])) continue;

				Ra_Vyx[I] = Ra_cor[xt][Na_idxPb[i]]*R_h2;
				wsVec Ra_corI = Ra_cor[Na_idxPb[i]];
				for (j=J=0 ; J<=I ; j++) {
					if (isMissingReal(Ra_fullY[Na_idxPb[j]])) continue;

					Ra_Vx[I][J] = Ra_corI[Na_idxPb[j]] * R_h2;
					J++;
				}
				Ra_Vx[I][I] += W1; //*= R_h2;
				I++;
			} else for (i=I=0 ; i<N_proband ; i++) {
				if (isMissingReal(Ra_fullY[Na_idxPb[i]])) continue;

				Ra_Vyx[I] = W0;
				for (j=J=0 ; J<=I ; j++) {
					if (isMissingReal(Ra_fullY[Na_idxPb[j]])) continue;

					Ra_Vx[I][J] = I==J ? R_h2 : W0;
					J++;
				}
				Ra_Vx[I][I] += W1; //*= R_h2;
				I++;
			}
			cSymMatrix&	M_invVx		= M_Vx.inv();
			wsSym		Ra_invVx	= M_invVx.get();
			M_Vx.rem();
#endif
//			sseUnmat(Ra_Vx, N_currProband);

			/* Get proband stat */
			// probSum  <- get_mu_var(eProb,T,l-1)
			resPAfmlar *Xp_r = _getProbandStat(getCorrelation(Cp_IO),
				Na_currIdxPb, N_currProband, Cp_IO->getPheno(), R_h2, Ra_cor);
			DEALLOC(Na_idxPb);
			DEALLOC(Na_currIdxPb);

			// mu=Vyx%*%invVx%*%probSum$mu
#ifndef FQLS_V2
			wsReal R_mu = sseVpMpV(Ra_Vyx, Ra_invVx, N_currProband,
				Xp_r->Ra_muHat);
			// invVxDec <- invVx%*%probSum$va%*%invVx
			wsSym Ra_invVx_dec = sseMMMt(Ra_invVx, N_currProband,
				N_currProband, Xp_r->Ra_varHat, N_currProband, N_currProband); /* checked */
#else
			wsReal R_mu = sseVpSpV(Ra_Vyx, Ra_invVx, N_currProband,
				Xp_r->Ra_muHat);
			// invVxDec <- invVx%*%probSum$va%*%invVx
			wsSym Ra_invVx_dec = sseSSS(Ra_invVx, N_currProband,
				Xp_r->Ra_varHat, N_currProband); /* checked */
#endif

			deallocMatrix(Xp_r->Ra_varHat, N_currProband, (void *)1);
			sseFree(Xp_r->Ra_muHat);
			DEALLOC(Xp_r);
#ifndef FQLS_V2
			sseMsM(Ra_invVx, Ra_invVx_dec, Ra_invVx, N_currProband);
#else
			sseSsS(Ra_invVx, Ra_invVx_dec, Ra_invVx, N_currProband);
#endif
			sseUnmat(Ra_invVx_dec, N_currProband);

#ifndef FQLS_V2
			// var=eNprob$vars+Vyx%*%(invVx-invVxDec)%*%Vxy
			wsReal R_var = (Ra_cor ? (Ra_cor[xt][xt]*R_h2+W1) : (R_h2+W1)) +
				sseVpMpV(Ra_Vyx, Ra_invVx, N_currProband);
			sseUnmat(Ra_invVx, N_currProband);
#else
			// var=eNprob$vars+Vyx%*%(invVx-invVxDec)%*%Vxy
			wsReal R_var = (Ra_cor ? (Ra_cor[xt][xt]*R_h2+W1) : (R_h2+W1)) +
				sseVpSpV(Ra_Vyx, Ra_invVx, N_currProband);
			delete& M_invVx;
#endif
			sseFree(Ra_Vyx);

			// adjY[j]<- adjY[j]-pnorm( (T-nprobSum$mu[1])/sqrt(nprobSum$va[1]),lower=F )
			Ra_newI_a[xt] = Ra_adjY[xt] -
				(wsReal)(pnorm((R_thr - R_mu)/sqrt(R_var), 0.0, 1.0, 0));
		}
	} else {
		/* In case of single sample family */
		wsUint		N_prev;
		wsReal*		Ra_prev	= loadRealValues(OPT_STRING(prevalence), &N_prev);
		vSampPtr&	Xa_samp	= Cp_IO->getSample();

		for (k=0 ; k<N_szFam ; k++) {
			wsUint xt = Na_idxMem[k];

			if (N_prev == 1) {
				Ra_newI_a[xt] = Ra_adjY[xt] - Ra_prev[0];
			} else {
				char N_sex = Xa_samp[xt]->N_sex;

				/* If there one individual */
				switch (N_sex) {
				case 1: case 2:
					Ra_newI_a[xt] = Ra_adjY[xt] - Ra_prev[N_sex-1];
					break;
				default:
					halt("Invalid sex value[%d] for two-prevalence",
						N_sex);
				}
			}
		}
		sseFree(Ra_prev);
	}

	return Ra_newI_a;
}

resPAfmlar* _getProbandStat(cAnalysis* Cp_anaCorr, wsUint* N_idxProband,
							wsUint N_pb, wsVecCst Ra_fullY, wsReal R_h2, wsSym Ra_cor/*=NULL*/)
{
	/* Set R_thr via `prevalence` */
	wsUint N_prev;
	wsReal *Ra_prev = loadRealValues(OPT_STRING(prevalence), &N_prev);
	wsReal	R_prev = Ra_prev[0];
	/* lower=T */
	wsReal	R_thr = (wsReal)qnorm(R_prev, 0, 1, 0, 0);
	sseFree(Ra_prev);

	wsReal	**Ra_fullCor	= NULL; /* CHECKED */
	// Ra_ret[0] = mu_hat (length N_pb)
	// Ra_ret[1] = var_hat (length N_pb)

	resPAfmlar
		*Xa_ret			= new resPAfmlar;
	//	MULTI_MALLOC(Xa_ret, resPAformula, 1);

	sseMalloc(Xa_ret->Ra_muHat, wsReal, N_pb);
	Xa_ret->Ra_varHat = sseMatrix(N_pb, N_pb);

	/* Get correlation matrix */
	if (Ra_cor == NULL)
		Ra_fullCor = getFullCorMat(Cp_anaCorr);
	else
		Ra_fullCor = Ra_cor;

	if (N_pb == 1) {
		/* In case of single proband */
		int N_idx = N_idxProband[0];

		// if(phen==0) {
		if (Ra_fullY[N_idx] == WISARD_UNAFFECTED)
			// mu_hat <- -dnorm(T)/pnorm(T,lower.tail=T)
				Xa_ret->Ra_muHat[0] = (wsReal)(-dnorm(R_thr)/(1.-pnorm(R_thr)));
		else if (Ra_fullY[N_idx] == WISARD_AFFECTED)
			// mu_hat <-  dnorm(T)/pnorm(T,lower.tail=F)
			Xa_ret->Ra_muHat[0] = (wsReal)(dnorm(R_thr)/pnorm(R_thr, 0.0, 1.0, 0));
		else {
			halt("SYSERR: Continuous phenotype cannot be handled with this method");
			//halt("Unhandled here, implement is needed");
		}

		// var_hat <- 1-mu_hat^2+mu_hat*T
		Xa_ret->Ra_varHat[0][0] = (wsReal)(1.0 - SQR(Xa_ret->Ra_muHat[0]) +
			Xa_ret->Ra_muHat[0]*R_thr);
	} else {
		/* In case of multiple probands */

		/* Make subset of phi matrix with probands */
		wsReal **Ra_subPhi = sseMatrix(N_pb, N_pb);
		if (Ra_fullCor) for (wsUint i=0 ; i<N_pb ; i++) {
			wsVec Ra_fi = Ra_fullCor[N_idxProband[i]];
			for (wsUint j=0 ; j<N_pb ; j++)
				Ra_subPhi[i][j] = Ra_fi[N_idxProband[j]]*R_h2;
			Ra_subPhi[i][i] += W1;
		} else for (wsUint i=0 ; i<N_pb ; i++) {
			for (wsUint j=0 ; j<N_pb ; j++)
				Ra_subPhi[i][j] = i==j ? R_h2 : W0;
			Ra_subPhi[i][i] += W1;
		}

		/* Do cholesky decomposition */
		// G  <- t(chol(vars))
		chol(Ra_subPhi, N_pb);
		wsReal **Ra_Ginv = SO_invMatrix(Ra_subPhi, N_pb);

		// TG <- c(solve(G)%*%rep(T,nSize))
		wsReal *Ra_TG = NULL;
		sseCalloc(Ra_TG, wsReal, N_pb);
		for (wsUint i=0 ; i<N_pb ; i++) for (wsUint j=0 ; j<=i ; j++)
			Ra_TG[i] += Ra_Ginv[i][j]*R_thr;
		sseUnmat(Ra_Ginv, N_pb);

		/* For each proband */
		wsReal *Ra_A = NULL;
		sseMalloc(Ra_A, wsReal, N_pb);
		for (wsUint i=0 ; i<N_pb ; i++) {
			int N_idx = N_idxProband[i];

			if (Ra_fullY[N_idx] == WISARD_UNAFFECTED)
				// A[oo1,1]  <- -dnorm(TG[oo1])/pnorm(TG[oo1],lower.tail=T)
					Ra_A[i] = (wsReal)(-dnorm(Ra_TG[i])/pnorm(Ra_TG[i]));
			else if (Ra_fullY[N_idx] == WISARD_AFFECTED)
				Ra_A[i] = (wsReal)(dnorm(Ra_TG[i])/pnorm(Ra_TG[i], 0.0, 1.0, 0));
		}
		// mu_hat <- G%*%A
		sseFree(Xa_ret->Ra_muHat);
		Xa_ret->Ra_muHat = sseMpV(Ra_subPhi, N_pb, N_pb, Ra_A);

		// B <- 1 - A^2+A*TG
		for (wsUint i=0 ; i<N_pb ; i++)
			Ra_A[i] = W1 - SQR(Ra_A[i]) + Ra_A[i]*Ra_TG[i];

		// var_hat <- G%*%diag(as.numeric(B))%*%t(G)
		//         [ g11*b1*g11
		// GBG^T = [ g21*b1*g11 g21*b1*g21 + g22*b2*g22 ...
		//         [     ...              ...
		//         [ gn1*b1*g11 gn1*b1*g21 + gn2*b2*g22
		// 
		// Column X_1 1st, b1*g11
		//        X_2 2nd, b1*g21, b2*g22
		//        X_3 ith, b1*gi1, bk*gik, bi*gii
		//           /* Depends on j=column X_j */
		//           
		// GBG^T_i1 = g_i1*X_1                           (i >= 1)
		//      _i2 = [ g_i1 g_i2 ]^T * X_2              (i >= 2)
		//      _ij = [ g_i1 ... g_ik ... g_ij ]^T * X_j (i >= j)
		sseLDLt(Ra_subPhi, N_pb, Ra_A, Xa_ret->Ra_varHat);
		sseUnmat(Ra_subPhi, N_pb);
		sseFree(Ra_A);
		sseFree(Ra_TG);
	}

	return Xa_ret;
}

int fetchProband(xSample *Xp_s, void *Vp_data)
{
	vInt *Xp_idxPb = (vInt *)Vp_data;

	if (Xp_s->N_idx == SAMP_NODATA)
		return 0;
	if (Xp_s->B_isProband)
		Xp_idxPb->push_back(Xp_s->N_idx);

	/* No meaning what it returns */
	return 0;
}

} // End namespace ONETOOL
