#include "utils/matrix.h"
#include "analyses/mfqls.h"
#include "analyses/setmgr.h"
#include "analyses/emai.h"
#include "analyses/fqls.h"
#include "analyses/pddt.h"
#include "utils/dcdflib.h"
#include "analyses/gstest.h"
#include "utils/stat.h"

namespace ONETOOL {

/* !!! FQLS-specific analysis !!! */
#if (TOOLSET_TYPE & FUNCTION_MFQLS) || ((TOOLSET_TYPE & 0x100) == 0x100)

int _doTest(xMqlsTest *Xp_t);
int _doTest(cIO *Cp_IO, wsStrCst Sp_name, vInt& Xa_snps, xMqlsRes &X_res,
	vector<wsMat>& Ma_TtA, vector<cSymMatrix*>& Ma_TtAphiAtT_inv,
	const char* Ba_misPhCv, wsUint N_anaSamp);

int thrMFQLSvariant(int N_idx, void *Vp_shareData, void *Vp_data,
	  void *Vp_result)
{
	xAnaThread*		Xp_at			= (xAnaThread *)Vp_shareData;
	cMFqlsAnalysis*	Cp_mfqls		= (cMFqlsAnalysis *)(Xp_at->Vp_data);
	cIO*			Cp_IO			= Xp_at->Cp_IO;
	vVariant&		Xa_SNP			= Cp_IO->getVariant();
	int*			Np_data			= (int *)Vp_data;
	wsUint			N_s				= (wsUint)Np_data[0];
	wsUint			N_e				= (wsUint)Np_data[1];

	/* Return pointer */
	wsMat			Ra_rets			= (wsMat)Vp_result;
	/* Check missing status */
	wsUint			N_anaSamp		= 0;
/**/char*			Ba_misPhCv		= Cp_IO->getPheCovMissing(&N_anaSamp);
	vector<wsMat>&	Ma_TtA			= Cp_mfqls->getTtA();
	vector<cSymMatrix*>&
					Ma_TtAphiAtT_inv = Cp_mfqls->getTtAphiAtT();

	//LOG("Thread %d : Range %d~%d\n", N_idx, N_s, N_e);
	for (wsUint j=N_s ; j<N_e ; j++) {
		xVariant&	X_var = Xa_SNP[j];
		cTimer		t;
		vInt		X_snp;
		xMqlsRes	X_res = { { W0, }, { W0, }, 0 };
		X_snp.push_back(j);

		t.start();
		_doTest(Cp_IO, X_var.name, X_snp, X_res, Ma_TtA, Ma_TtAphiAtT_inv,
			Ba_misPhCv, N_anaSamp);

		Ra_rets[1][j] = X_res.R_t[0];
		Ra_rets[4][j] = X_res.R_p[0];

		if (OPT_ENABLED(time))
			Ra_rets[0][j] = t.get();

		if (OPT_ENABLED(fqls)) {
			Ra_rets[2][j] = X_res.R_t[1];
			Ra_rets[3][j] = X_res.R_t[2];
			Ra_rets[5][j] = X_res.R_p[1];
			Ra_rets[6][j] = X_res.R_p[2];
		}

		/* Count done */
		Np_data[2]++;

		/* Print progress */
		if (N_idx==0 && (j%100)==0) {
			wsUint N_sum = 0;
			for (int x=2,y=0 ; y<OPT_NUMBER(thread) ; x+=3,y++)
				N_sum += Np_data[x];
			notice("[%d/%d] variants tested...\r", N_sum, Xa_SNP.size());
		}
	}

	return 0;
}

int thrMFQLSgene(int N_idx, void *Vp_shareData, void *Vp_data,
					void *Vp_result)
{
	xAnaSetThread*	Xp_at			= (xAnaSetThread *)Vp_shareData;
	cMFqlsAnalysis*	Cp_mfqls		= (cMFqlsAnalysis *)(Xp_at->Vp_data);
	cSetManagerAnalysis*
					Cp_anaGS		= Cp_mfqls->getSet();
	cIO*			Cp_IO			= Xp_at->Cp_IO;
	int*			Np_data			= (int *)Vp_data;
	wsUint			N_s				= (wsUint)Np_data[0];
	wsUint			N_e				= (wsUint)Np_data[1];

	/* Return pointer */
	wsMat			Ra_rets			= (wsMat)Vp_result;
	/* Check missing status */
	wsUint			N_anaSamp		= 0;
/**/char*			Ba_misPhCv		= Cp_IO->getPheCovMissing(&N_anaSamp);
	vector<wsMat>&	Ma_TtA			= Cp_mfqls->getTtA();
	vector<cSymMatrix*>&
		Ma_TtAphiAtT_inv = Cp_mfqls->getTtAphiAtT();

	/* If it is gene-based */
	mGeneDef&	Xa_set	= Cp_anaGS->getGeneDef();

	//LOG("Thread %d : Range %d~%d\n", N_idx, N_s, N_e);
	/* Set position */
	mGeneDef_it it = Xa_set.begin();
	advance(it, N_s);
	for (wsUint j=N_s ; j<N_e ; j++,it++) {
		if (it == Xa_set.end()) halt("SYSERR: Gene definition range out");
		cTimer		t;
		xMqlsRes	X_res = { { W0, }, { W0, }, 0 };

		t.start();
#ifdef ONETOOL
		wsUint N_snp = ONETOOL::_doTest(Cp_IO, it->first.c_str(), it->second,
			X_res, Ma_TtA, Ma_TtAphiAtT_inv, Ba_misPhCv, N_anaSamp);
#else
		wsUint N_snp = _doTest(Cp_IO, it->first.c_str(), it->second,
			X_res, Ma_TtA, Ma_TtAphiAtT_inv, Ba_misPhCv, N_anaSamp);
#endif

		Ra_rets[1][j] = X_res.R_t[0];
		Ra_rets[4][j] = X_res.R_p[0];

		if (OPT_ENABLED(time))
			Ra_rets[0][j] = t.get();

		if (OPT_ENABLED(fqls)) {
			Ra_rets[2][j] = X_res.R_t[1];
			Ra_rets[3][j] = X_res.R_t[2];
			Ra_rets[5][j] = X_res.R_p[1];
			Ra_rets[6][j] = X_res.R_p[2];
		}
		Ra_rets[7][j] = (wsReal)N_snp;

		/* Count done */
		Np_data[2]++;

		/* Print progress */
		if (N_idx==0 && (j%50)==0) {
			wsUint N_sum = 0;
			for (int x=2,y=0 ; y<OPT_NUMBER(thread) ; x+=3,y++)
				N_sum += Np_data[x];
			notice("[%d/%d] genes tested...\r", N_sum, Xa_set.size());
		}
	}
	DEALLOC(Ba_misPhCv);

	return 0;
}

resPAfmlar* cMFqlsAnalysis::_getProbandStat(wsUint *N_idxProband,
	wsUint N_pb, wsRealCst *Ra_fullY, wsReal R_h2, wsReal **Ra_cor,
	wsReal R_p1, wsReal R_p2)
{
	vSampPtr&		Xa_samp		= Cp_IO->getSample();
	wsMat			Ra_fullCor	= NULL; /* CHECKED */
	// Ra_ret[0] = mu_hat (length N_pb)
	// Ra_ret[1] = var_hat (length N_pb)
	
	resPAfmlar*	Xa_ret		= new resPAfmlar;
//	MULTI_MALLOC(Xa_ret, resPAformula, 1);

	sseMalloc(Xa_ret->Ra_muHat, wsReal, N_pb);
	Xa_ret->Ra_varHat = sseMatrix(N_pb, N_pb);

	/* Get correlation matrix */
	Ra_fullCor = Ra_cor;

	if (N_pb == 1) {
		/* In case of single proband */
		int		N_idx	= N_idxProband[0];
		wsReal	R_thr	= (wsReal)qnorm(
			R_p2 != R_p2 ? R_p1 : (Xa_samp[N_idx]->N_sex == 1 ? R_p1 : R_p2),
			0, 1, 0, 0);

		// if(phen==0) {
		if (Ra_fullY[N_idx] == WISARD_UNAFFECTED)
			// mu_hat <- -dnorm(T)/pnorm(T,lower.tail=T)
			Xa_ret->Ra_muHat[0] = (wsReal)(-dnorm(R_thr)/pnorm(R_thr));
		else if (Ra_fullY[N_idx] == WISARD_AFFECTED)
			// mu_hat <-  dnorm(T)/pnorm(T,lower.tail=F)
			Xa_ret->Ra_muHat[0] = (wsReal)(dnorm(R_thr)/pnorm(R_thr, 0.0, 1.0, 0));
		else
			halt("Unhandled here, implement is needed");

		// var_hat <- 1-mu_hat^2+mu_hat*T
		Xa_ret->Ra_varHat[0][0] = (wsReal)(1.0 - SQR(Xa_ret->Ra_muHat[0]) +
			Xa_ret->Ra_muHat[0]*R_thr);
	} else {
		/* In case of multiple probands */

		/* Make subset of phi matrix with probands */
		wsReal **Ra_subPhi = sseMatrix(N_pb, N_pb);
		for (wsUint i=0 ; i<N_pb ; i++) {
			wsReal *Ra_ii = Ra_fullCor[N_idxProband[i]];

			for (wsUint j=0 ; j<N_pb ; j++)
				Ra_subPhi[i][j] = Ra_ii[N_idxProband[j]] * R_h2;
			Ra_subPhi[i][i] += W1;
		}

		/* Do cholesky decomposition */
		// G  <- t(chol(vars))
		chol(Ra_subPhi, N_pb);
		wsReal **Ra_Ginv = SO_invMatrix(Ra_subPhi, N_pb);

		// TG <- c(solve(G)%*%rep(T,nSize))
		wsReal *Ra_TG = NULL;
		sseCalloc(Ra_TG, wsReal, N_pb);
		if (R_p2 != R_p2) for (wsUint i=0 ; i<N_pb ; i++) { /* Single prevalence */
			wsReal R_thr = (wsReal)qnorm(R_p1, 0, 1, 0, 0);
			for (wsUint j=0 ; j<=i ; j++)
				Ra_TG[i] += Ra_Ginv[i][j]*R_thr;
		} else for (wsUint i=0 ; i<N_pb ; i++) { /* Double prevalence */
			wsReal R_thr = (wsReal)qnorm(
				Xa_samp[N_idxProband[i]]->N_sex == 1 ? R_p1 : R_p2,
				0, 1, 0, 0);
			
			//Xa_samp[N_idxProband[i]]->N_sex ? R_p1 : R_p2;
			for (wsUint j=0 ; j<=i ; j++)
				Ra_TG[i] += Ra_Ginv[i][j]*R_thr;
		}
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

cStdMatrix cMFqlsAnalysis::_getFQLSphe(cVector& V_prev)
{
	wsUintCst		N_pheno		= Cp_IO->sizePheno();
	wsUintCst		N_sample	= Cp_IO->sizeSample();
	wsReal*		Rp_prev		= V_prev.get();
	cStdMatrix	M_FY(N_pheno, N_sample, Cp_IO->getPhenos(), MATDEL_NONE);
	cStdMatrix	M_adjFYfqls(N_pheno, N_sample);

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

		wsUint		i = 0;
		switch (N_sz) {
		case 1: /* Equal to male/female */
			FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
				if (Ra_fullY[i] == WISARD_UNAFFECTED)
					Ra_ret[i] = W0;
				else if (Ra_fullY[i] == WISARD_AFFECTED)
					Ra_ret[i] = W1;
				else
					Ra_ret[i] = Ra_prev[0];
			}
			break;
		case 2: /* Non-equal to male/female */
			FOREACHDO (vSampPtr_it, Xa_samp, it, i++) {
				if (Ra_fullY[i] == WISARD_UNAFFECTED) {
					Ra_ret[i] = W0;
				} else if (Ra_fullY[i] == WISARD_AFFECTED) {
					Ra_ret[i] = W1;
				} else {
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
	}

	return M_adjFYfqls;
}

wsReal* cMFqlsAnalysis::_buildFamII(xFamily &X_fam, wsUint N_szFam,
	wsUint *Na_idxMem, wsReal *Ra_adjY, wsReal **Ra_cor,
	wsReal *Ra_newII_a, wsRealCst *Ra_fullY, wsReal R_h2, wsReal R_p1,
	wsReal R_p2)
{
	wsUint		i;
	vSampPtr&	Xa_samp = Cp_IO->getSample();

	/* Fetch the information about proband */
	vInt Xa_idxProband;
	X_fam.visit(fetchProband, &Xa_idxProband);
	wsUint N_proband = (wsUint)Xa_idxProband.size();

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
	cSymMatrix	M_Vx(N_proband);
	cSymMatrix	M_Vy(N_nonPb);
	cStdMatrix	M_Vyx(N_nonPb, N_proband);
	wsMat		Ra_Vyx	= M_Vyx.get(); /* checked */
	wsSym		Ra_Vy	= M_Vy.get(); /* Checked */
	wsSym		Ra_Vx	= M_Vx.get(); /* checked */

	LOOP (i, N_proband) {
		for (wsUint j=0 ; j<N_nonPb ; j++)
			Ra_Vyx[j][i] = Ra_cor[Na_idxNonPb[j]][Na_idxProband[i]]*R_h2;
		for (wsUint j=0 ; j<i ; j++)
			Ra_Vx[i][j] = Ra_cor[Na_idxProband[i]][Na_idxProband[j]]*R_h2;
		Ra_Vx[i][i] = Ra_cor[Na_idxProband[i]][Na_idxProband[i]]*R_h2 + W1;
	}
	LOOP (i, N_nonPb) {
		for (wsUint j=0 ; j<i ; j++)
			Ra_Vy[i][j] = Ra_cor[Na_idxNonPb[i]][Na_idxNonPb[j]]*R_h2;
		Ra_Vy[i][i] = Ra_cor[Na_idxNonPb[i]][Na_idxNonPb[i]]*R_h2 + W1;
	}

	cSymMatrix	&M_Vxi = M_Vx.inv();
	M_Vx.rem();

	/* Get proband stat */
	// probSum  <- get_mu_var(eProb,T,l-1)
	resPAfmlar*	Xp_r = _getProbandStat(Na_idxProband, N_proband,
		Ra_fullY, R_h2, Ra_cor, R_p1, R_p2); /* checked */
	cVector		V_muHat(Xp_r->Ra_muHat, N_proband, 1);
	cSymMatrix	M_varHat(Xp_r->Ra_varHat, N_proband, 1);

	// mu=Vyx%*%invVx%*%probSum$mu [#nonpb*1]
	cVector	V_mu	= M_Vyx.MV(M_Vxi, V_muHat);
	wsReal*	Ra_mu	= V_mu.get();
// 	wsReal **R_mu = sseMMMt(Ra_Vyx, N_nonPb, N_proband,
// 		Ra_invVx, N_proband, N_proband,
// 		&(Xp_r->Ra_muHat), 1, N_proband); /* checked */

	// invVxDec <- invVx%*%probSum$va%*%invVx
	cSymMatrix M_Vxi_dec = M_Vxi.MMt(M_varHat);
//	wsReal **Ra_invVx_dec = sseMMMt(Ra_invVx, N_proband,
//		N_proband, Xp_r->Ra_varHat, N_proband, N_proband,
//			Ra_invVx, N_proband, N_proband); /* checked */
	deallocMatrix(Xp_r->Ra_varHat, N_proband, (void *)1);
	sseFree(Xp_r->Ra_muHat);
	DEALLOC(Xp_r);

	M_Vxi -= M_Vxi_dec;
	M_Vxi_dec.rem();
// 	sseMsM(Ra_invVx, Ra_invVx_dec, Ra_invVx, N_proband);
//	sseUnmat(Ra_invVx_dec, N_proband);

	// var=eNprob$vars+Vyx%*%(invVx-invVxDec)%*%Vxy
	cSymMatrix M_VyP = M_Vyx.MMt(M_Vxi);
// 	wsReal **Ra_Vy_P = sseMMMt(Ra_Vyx, N_nonPb, N_proband,
// 		Ra_invVx, N_proband, N_proband,
// 		Ra_Vyx, N_nonPb, N_proband); /* checked */
	M_Vy += M_VyP;
//	sseMaM(Ra_Vy, Ra_Vy_P, Ra_Vy, N_nonPb);
	M_VyP.rem();
	M_Vxi.rem();
	M_Vyx.rem();
// 	sseUnmat(Ra_Vy_P, N_nonPb);
// 	sseUnmat(Ra_invVx, N_proband);
// 	sseUnmat(Ra_Vyx, N_nonPb);

	/* adjY for non-proband */
	if (R_p2 != R_p2) LOOP (i, N_nonPb) { /* Single prevalence */
		wsReal R_thr = (wsReal)qnorm(R_p1, 0, 1, 0, 0);

		// adjY[-o]-pnorm( (T-nprobSum$mu)/sqrt(nprobSum$va[obs]),lower=F )
		if (isMissingReal(Ra_fullY[Na_idxNonPb[i]]))
			Ra_newII_a[Na_idxNonPb[i]] = W0;
		else
			Ra_newII_a[Na_idxNonPb[i]] = Ra_adjY[Na_idxNonPb[i]]
				- (wsReal)(pnorm((R_thr - Ra_mu[i])/sqrt(Ra_Vy[i][i]), 0.0, 1.0, 0));
//		- (wsReal)(pnorm((R_thr - R_mu[i][0])/sqrt(Ra_Vy[i][i]), 0.0, 1.0, 0));
	} else LOOP (i, N_nonPb) { /* Double prevalence */
		wsReal R_thr = (wsReal)qnorm(
			Xa_samp[Na_idxNonPb[i]]->N_sex == 1 ? R_p1 : R_p2, 0, 1, 0, 0);

		// adjY[-o]-pnorm( (T-nprobSum$mu)/sqrt(nprobSum$va[obs]),lower=F )
		if (isMissingReal(Ra_fullY[Na_idxNonPb[i]]))
			Ra_newII_a[Na_idxNonPb[i]] = W0;
		else
			Ra_newII_a[Na_idxNonPb[i]] = Ra_adjY[Na_idxNonPb[i]]
				- (wsReal)(pnorm((R_thr - Ra_mu[i])/sqrt(Ra_Vy[i][i]), 0.0, 1.0, 0));
//		- (wsReal)(pnorm((R_thr - R_mu[i][0])/sqrt(Ra_Vy[i][i]), 0.0, 1.0, 0));
	}
	M_Vy.rem();
//	sseUnmat(Ra_Vy, N_nonPb);
	V_mu.rem();
//	sseUnmat(R_mu, N_nonPb);

	/* adjY for proband */
	if (R_p2 != R_p2) LOOP (i, N_proband) {
		wsReal R_prev = R_p1;

		// adjY[o]-prev
		if (isMissingReal(Ra_fullY[Na_idxProband[i]]))
			Ra_newII_a[Na_idxProband[i]] = W0;
		else
			Ra_newII_a[Na_idxProband[i]] = Ra_adjY[Na_idxProband[i]]
				- R_prev;
	} else LOOP (i, N_proband) {
		wsReal R_prev = Xa_samp[Na_idxNonPb[i]]->N_sex == 1 ? R_p1 : R_p2;

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

wsReal* cMFqlsAnalysis::_buildFamI(xFamily &X_fam, wsUint N_szFam,
	wsUint *Na_idxMem, wsReal *Ra_adjY, wsReal **Ra_cor,
	wsReal *Ra_newI_a, wsRealCst *Ra_fullY, wsReal R_h2, wsReal R_p1,
	wsReal R_p2)
{
	wsUint		i, j, k;
	vSampPtr&	Xa_samp		= Cp_IO->getSample();
	wsUint		N_proband	= N_szFam-1;

	if (N_proband) for (k=0 ; k<N_szFam ; k++) {
		wsUint	xt				= Na_idxMem[k];
		wsReal	R_thr			= (wsReal)qnorm(
			R_p2 != R_p2 ? R_p1 : (Xa_samp[xt]->N_sex == 1 ? R_p1 : R_p2),
			0, 1, 0, 0);

		/* Members except current member will be proband */
		wsUint*	Na_idxPb		= NULL; /* checked */
		wsUint*	Na_currIdxPb	= NULL;
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

		/* Get Vyx */
		// Vxy      <- V[-j,j,drop=F]
		cSymMatrix	M_Vx(N_currProband);
		cVector		V_Vyx(N_currProband);
// 		wsReal *Ra_Vyx	= V_Vyx.get(); //NULL; /* Checked */
// 		wsReal **Ra_Vx	= M_VxsseMatrix(N_currProband, N_currProband); /* checked */
		wsReal *Ra_Vyx	= V_Vyx.get(); //NULL; /* Checked */
		wsSym	Ra_Vx	= M_Vx.get();
//		sseMalloc(Ra_Vyx, wsReal, N_currProband);
		for (i=I=0 ; i<N_proband ; i++) {
			wsUint N_ii = Na_idxPb[i];
			if (isMissingReal(Ra_fullY[N_ii])) continue;

			Ra_Vyx[I] = Ra_cor[xt][N_ii]*R_h2;
			for (j=J=0 ; J<I ; j++) {
				if (isMissingReal(Ra_fullY[Na_idxPb[j]])) continue;

				Ra_Vx[I][J] = Ra_cor[N_ii][Na_idxPb[j]] * R_h2;
				J++;
			}
			Ra_Vx[I][I] = Ra_cor[N_ii][N_ii] * R_h2 + W1;
			I++;
		}
		cSymMatrix	&M_Vxi = M_Vx.inv();
		M_Vx.rem();
//		wsReal **Ra_invVx = invSymMat(Ra_Vx, N_currProband); /* checked */
//		sseUnmat(Ra_Vx, N_currProband);

		/* Get proband stat */
		// probSum  <- get_mu_var(eProb,T,l-1)
		resPAfmlar *Xp_r = _getProbandStat(Na_currIdxPb, N_currProband,
			Ra_fullY, R_h2, Ra_cor, R_p1, R_p2);
		cVector		V_muHat(Xp_r->Ra_muHat, N_currProband);
		cSymMatrix	M_varHat(Xp_r->Ra_varHat, N_currProband);
		DEALLOC(Na_idxPb);
		DEALLOC(Na_currIdxPb);

		// mu=Vyx%*%invVx%*%probSum$mu
// 		wsReal R_mu = sseVpMpV(Ra_Vyx, Ra_invVx, N_currProband,
// 			Xp_r->Ra_muHat);
		wsReal R_mu = V_Vyx.qf(M_Vxi, V_muHat);

		// invVxDec <- invVx%*%probSum$va%*%invVx
		cSymMatrix	M_Vxi_dec = M_Vxi.MMt(M_varHat);
// 		wsReal **Ra_invVx_dec = sseMMMt(Ra_invVx, N_currProband,
// 			N_currProband, Xp_r->Ra_varHat, N_currProband, N_currProband,
// 			Ra_invVx, N_currProband, N_currProband); /* checked */
//		deallocMatrix(Xp_r->Ra_varHat, N_currProband, (void *)1);
//		sseFree(Xp_r->Ra_muHat);
		M_varHat.rem();
		V_muHat.rem();
		DEALLOC(Xp_r);
		M_Vxi -= M_Vxi_dec;
//		sseMsM(Ra_invVx, Ra_invVx_dec, Ra_invVx, N_currProband);
		M_Vxi_dec.rem();
//		sseUnmat(Ra_invVx_dec, N_currProband);

		// var=eNprob$vars+Vyx%*%(invVx-invVxDec)%*%Vxy
		wsReal R_var = (Ra_cor[xt][xt]*R_h2+W1) +
//			sseVpMpV(Ra_Vyx, Ra_invVx, N_currProband, Ra_Vyx);
			V_Vyx.qf(M_Vxi);
		delete &M_Vxi;
		V_Vyx.rem();
//		sseUnmat(Ra_invVx, N_currProband);
//		sseFree(Ra_Vyx);

		// adjY[j]<- adjY[j]-pnorm( (T-nprobSum$mu[1])/sqrt(nprobSum$va[1]),lower=F )
		Ra_newI_a[xt] = Ra_adjY[xt] -
			(wsReal)(pnorm((R_thr - R_mu)/sqrt(R_var), 0.0, 1.0, 0));
	} else {
		/* In case of single family */
		vSampPtr	&Xa_samp	= Cp_IO->getSample();

		for (k=0 ; k<N_szFam ; k++) {
			wsUint xt = Na_idxMem[k];

			if (NA(R_p2)) {
				Ra_newI_a[xt] = Ra_adjY[xt] - R_p1;
			} else {
				char N_sex = Xa_samp[xt]->N_sex;

				/* If there one individual */
				switch (N_sex) {
				case 1:
					Ra_newI_a[xt] = Ra_adjY[xt] - R_p1;
					break;
				case 2:
					Ra_newI_a[xt] = Ra_adjY[xt] - R_p2;
					break;
				default:
					halt("Invalid sex value[%d] for two-prevalence",
						N_sex);
				}
			}
		}
	}

	return Ra_newI_a;
}

cMFqlsAnalysis::cMFqlsAnalysis(cIO *Cp_IO, cAnalysis *Cp_inpCorr,
	cSetManagerAnalysis *Cp_inpAnaGS, cVariantMap *Cp_inpMM) :
	cAnalysis(Cp_IO), N_pheno(Cp_IO->sizePheno())
{
	wsUintCst	N_samp			= Cp_IO->sizeSample();
	wsUint	N_nonPhCvMiss	= 0;
	wsUint	N_nonCvMiss	= 0;
	wsMat	Ra_pheno		= Cp_IO->getPhenos();
	wsUint	N_pheno			= Cp_IO->sizePheno();
	char	B_indep			= 0;

	/* Assign variables */
	Cp_anaGS	= Cp_inpAnaGS;
	Cp_mm		= Cp_inpMM;

	/* One of --blup/--prevalence is required if dichotomous */
	char	B_isCont	= Cp_IO->isContinuous(0);
	for (wsUint i=1 ; i<N_pheno ; i++)
		if (B_isCont != Cp_IO->isContinuous(i))
			halt("Phenotypes must be all dichotomous/all continuous!");

	/* --fqls checking */
	if (OPT_ENABLED(fqls)) {
		if (B_isCont)
			halt_fmt(WISARD_CANT_DO_W_SOMEOPT, "Continuous phenotype", "--fqls");
		ASSERT_OPTION(prevalence);
		ASSERT_OPTION(heri);
	}
// 	if (!B_isCont) {
// 		if (!IS_ASSIGNED(blup) && !IS_ASSIGNED(prevalence))
// 			halt_fmt(WISARD_CANT_MQLS_WO_BLUPORPREV);
// 	}
	if (IS_ASSIGNED(blup) && IS_ASSIGNED(prevalence))
		halt_fmt(WISARD_CANT_MQLS_W_BOTHBLUPPREV);
	if (SET_AVAIL() && IS_ASSIGNED(mqlsconsec))
		halt_fmt(WISARD_CANT_OPT_W_OTHEROPT, "--set", "--mqlsconsec");

	/* Check missing status */
/**/Ba_misPhCv	= Cp_IO->getPheCovMissing(&N_nonPhCvMiss);
	if (N_nonPhCvMiss == 0) halt("There is no sample to perform MQLS");
	N_anaSamp	= N_nonPhCvMiss;
	LOG("[%d] samples used in multi FQLS\n", N_anaSamp);

	/* Resize phenotype */
	M_Y.init(N_pheno, N_anaSamp);
	wsMat Ra_Y = M_Y.get();
	for (wsUint j=0 ; j<N_pheno ; j++)
		for (wsUint i=0,I=0 ; i<N_samp ; i++)
			if (!Ba_misPhCv[i])
				Ra_Y[j][I++] = Ra_pheno[j][i];
#ifdef MQLSvalidate
	exportMatrix("Y", Ra_Y, N_pheno, N_anaSamp);
#endif

	/* Get subset of phi matrix */
/**/wsReal **Ra_phi		= getFullCorMat(Cp_inpCorr);
	cMatrix	*Cp_phiMQ	= NULL;
	if (Ra_phi == NULL) {
		Cp_phiMQ	= new cIdtMatrix(N_anaSamp);
		B_indep		= 1;
	} else {
		Cp_phiMQ	= new cSymMatrix(N_anaSamp);
		wsMat Ra_phiMQ = Cp_phiMQ->get();
		for (wsUint i=0,I=0 ; i<N_samp ; i++) {
			if (Ba_misPhCv[i]) continue;
			for (wsUint j=0,J=0 ; j<=i ; j++)
				if (!Ba_misPhCv[j])
					Ra_phiMQ[I][J++] = Ra_phi[i][j];
			I++;
		}
	}
#ifdef MQLSvalidate
	exportSymMat("phi", Ra_phiMQ, N_anaSamp);
#endif

	/* For --fqls */
	wsMat Ra_a1 = NULL;
	wsMat Ra_a2 = NULL;
	if (OPT_ENABLED(fqls)) {
		/* Prevalence */
		wsUint		N_prev;
		wsReal		*Ra_prev	= loadRealValues(OPT_STRING(prevalence), &N_prev);
		if ((N_prev != N_pheno && N_prev != (N_pheno<<1)) ||
			(N_pheno == 1 && N_prev != 1 && N_prev != 2))
			halt("Given number of --prev[%d] is not match with phenotypes[%d]",
				N_prev, N_pheno);
		char		B_prevMF	= N_prev == (N_pheno<<1);
		cVector		V_prev(Ra_prev, N_prev);

		cPDDTAnalysis	*Cp_pddt	= NULL;
		wsMat			Ra_ofsCor	= NULL;

		cStdMatrix	M_Yfqls = _getFQLSphe(V_prev);
		M_Yf1.init(N_pheno, N_samp);
		Ra_a1 = M_Yf1.get();
		if (Cp_IO->isProbandAssigned()) {
			M_Yf2.init(N_pheno, N_samp);
			Ra_a2 = M_Yf2.get();
		}

		wsUint		N_heri		= 0;
		wsReal		*Ra_heri	= loadRealValues(OPT_STRING(heri), &N_heri);
		if (N_heri != N_pheno)
			halt("Given number of --heri[%d] is not match with phenotypes[%d]",
			N_heri, N_pheno);

		if (!OPT_ENABLED(kinship) || IS_ASSIGNED(cor)) {
			/* If --fqlsnopddt, force use to PDDT */
			if (!OPT_ENABLED(fqlsnopddt)) {
				/* Calculate PDDT */
				Cp_pddt = new cPDDTAnalysis(getIO());
				Cp_pddt->run();
				Ra_ofsCor = Cp_pddt->getPDDT();
			}
		}
		if (Ra_ofsCor == NULL)
			Ra_ofsCor = Ra_phi;

		/* For all sample we can calculate each getStat for newI and newII */
		/* Family-wise, newI */
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
				wsReal R_p1 = B_prevMF ? Ra_prev[P<<1] : Ra_prev[P];
				wsReal R_p2 = B_prevMF ? Ra_prev[(P<<1)+1] : WISARD_NAN;

				if (Cp_IO->isProbandAssigned()) {
					wsReal *R = _buildFamII(X_fam, (wsUint)X_idxMem.size(),
						Na_idxMem, M_Yfqls.get()[P], Ra_ofsCor, Ra_a2[P],
						getIO()->getPheno(P), Ra_heri[P], R_p1, R_p2);
					if (!R)
						halt_fmt(WISARD_SYST_FAIL_GET_FQLS, ft->first.c_str(), "type 2");
				}
				wsReal *R = _buildFamI(X_fam, (wsUint)X_idxMem.size(),
					Na_idxMem, M_Yfqls.get()[P], Ra_ofsCor, Ra_a1[P],
					getIO()->getPheno(P), Ra_heri[P], R_p1, R_p2);
				if (!R)
					halt_fmt(WISARD_SYST_FAIL_GET_FQLS, ft->first.c_str(), "type 1");
			}
			DEALLOC(Na_idxMem);
		}

	}

	/* Get phiInv */
	cMatrix *Cp_phiInv = NULL;
#ifdef MQLSvalidate
	Cp_phiInv->file("phiInv");
#endif

	if (!IS_ASSIGNED(blup)) {
		/* For each phenotype, fitting EM-AI algorithm and get BLUP */
		wsUint	N_cov	= getIO()->sizeCovar();
		wsUint	N_prev	= 0;
/**/	wsReal*	Ra_prev	= loadRealValues(OPT_STRING(prevalence), &N_prev);

		if (N_prev == 0 || getIO()->isContinuous()) {
			if (!getIO()->isContinuous())
				LOG("Dichotomous phenotype assigned, but no prevalence assigned so LMM will be applied\n");
//			cEmAiAnalysisV2	**Cp_EMAIs	= NULL;

// 			wsReal	**Ra_cov	= sseMatrix(N_sampMQ, N_cov+1);
// 			wsReal	**Ra_covT	= sseMatrix(N_cov+1, N_sampMQ);
//			MAT_t	Ra_cov		= NULL;
			wsMat	Ra_covT		= NULL;

			/* Build covariate matrix */
			Ra_covT = getIO()->getFiltNullDesignMat(Ba_misPhCv, 1, N_anaSamp, 0);
// 			for (wsUint i=0,I=0 ; i<N_samp ; i++) {
// 				if (Ba_misPheno[i]) continue;
// 
// 				Ra_covT[0][I] = Ra_cov[I][0] = W1;
// 				for (wsUint j=0 ; j<N_cov ; j++)
// 					Ra_covT[j+1][I] = Ra_cov[I][j+1] = Ra_cov[j][I];
// 				I++;
// 			}
// 			
			/* ED */
			wsMat	Ra_eVec	= NULL;
			wsReal*	Ra_eVal	= NULL;
			getEigenPhi(Cp_phiMQ->get(), N_anaSamp, &Ra_eVec, &Ra_eVal);
			/* eVec = P^t, not P */
			cStdMatrix	M_eVec(N_anaSamp, N_anaSamp, Ra_eVec);
			cDiagMatrix	M_eVal(N_anaSamp, Ra_eVal);

//			MULTI_MALLOC(Cp_EMAIs, cEmAiAnalysisV2 *, N_pheno);
//			for (wsUint i=0 ; i<N_pheno ; i++) {
				cStdMatrix	M_Y(N_pheno, N_anaSamp, Ra_Y, MATDEL_NONE);
//				cSymMatrix	M_phiMQ(Cp_phiMQ->get(), N_sampMQ, 1);
				cStdMatrix	M_Xt(N_cov+1, N_anaSamp, Ra_covT);
// 				Cp_EMAIs[i] = new cEmAiAnalysisV2(Cp_IO, M_Xt, V_Y, M_phiMQ,
// 					!OPT_ENABLED(ml), NULL, NULL);
// 				Cp_EMAIs[i]->run();
// 				if (Cp_EMAIs[i]->getIsFailed())
// 					halt("LMM fitting with EM-AI algorithm was failed for [%d]th phenotype", i+1);
// 				REAL_c *Ra_BLUP = Cp_EMAIs[i]->getBLUP();

				char *Ba_isExcl = NULL;
				wsAlloc(Ba_isExcl, char, N_samp);
				for (wsUint x=0 ; x<N_samp ; x++)
					Ba_isExcl[x] = Ba_misPhCv[x];
// 				char *Ba_isIncl = NULL;
// 				MULTI_MALLOC(Ba_isIncl, char, N_samp);
// 				for (wsUint x=0 ; x<N_samp ; x++)
// 					Ba_isIncl[x] = !Ba_misPhCv[x];

				cEmAiAnalysisV2	C_EMAI(Cp_IO, M_Xt, M_Y, *Cp_phiMQ,
					!OPT_ENABLED(ml), Ba_isExcl, &M_eVec, &M_eVal);
				C_EMAI.run();
				if (C_EMAI.getIsFailed())
					halt("LMM fitting with EM-AI algorithm was failed");// for [%d]th phenotype", i+1);

				DEALLOC(Ba_isExcl);

				/* Substitute BLUP value */
				for (wsUint i=0 ; i<N_pheno ; i++) {
					wsRealCst *Ra_BLUP = C_EMAI.getBLUP(i);
					wsRealCst *Ra_pred = C_EMAI.getPred(i);
					for (wsUint j=0 ; j<N_anaSamp ; j++)
						Ra_Y[i][j] -= Ra_BLUP[j] + Ra_pred[j];
				}
				//cSymMatrix *Cp_sym = 
				transposeSelf(M_eVec.get(), N_anaSamp, N_anaSamp);
				wsSym Ra_pi = eigenInvSymMatrix(M_eVec.get(), M_eVal.get()[0],
					N_anaSamp);
				Cp_phiInv = new cSymMatrix(Ra_pi, N_anaSamp);
//			}
		} else {
			/* Otherwise, use prevalence */

			/* Load prevalence if possible */
			/* Check the number of prevalences is match with the number of phenotypes */
			if (Ra_prev && ((N_pheno > 1 && N_prev != N_pheno) ||
				(N_pheno == 1 && N_prev > 2)))
				halt_fmt(WISARD_INVL_NPREV_NPHENO, N_prev, N_pheno);

			/* Substitute prevalence to phenotype
			 * CASE = 1
			 * CTRL = -P/(1-P) */
			if (N_prev == 1 && N_pheno == 2) {
				wsUint I = 0, J = 0;
				vSampPtr&	Xa_samp = getIO()->getSample();
				FOREACHDO(vSampPtr_it, Xa_samp, it, J++) {
					if (Ba_misPhCv[J]) continue;

					if (Ra_Y[0][I] == WISARD_UNAFFECTED) {
						switch ((*it)->N_sex) {
						case 1: case 2:
							Ra_Y[0][I] = -Ra_prev[(*it)->N_sex - 1] /
								(W1 - Ra_prev[(*it)->N_sex - 1]);
							break;
						default: /* FIXME : How to treat in case of missing sex? */
							halt("Missing sex found, two-prevalence assignment");
						}
					}
					I++;
				}
			} else for (wsUint i=0 ; i<N_pheno ; i++)
				for (wsUint j=0 ; j<N_anaSamp ; j++)
					if (Ra_Y[i][j] == WISARD_UNAFFECTED)
						Ra_Y[i][j] = -Ra_prev[i] / (W1 - Ra_prev[i]);
			Cp_phiInv = &(Cp_phiMQ->inv());
		}
		sseFree(Ra_prev);
	} else
		Cp_phiInv = &(Cp_phiMQ->inv());

	/* Expand Ra_Y if requires */
	if (OPT_ENABLED(imputepheno)) {
/**/	Ba_misCov	= Cp_IO->getCovMissing(&N_nonCvMiss);
		wsMat Ra_newY = sseMatrix(N_pheno, N_nonCvMiss);

		LOG("--imputepheno found, sample size is adjusted from [%d] to [%d]\n",
			N_anaSamp, N_nonCvMiss);

		/* Impute phenotype missing to 0 */
		for (wsUint i=0 ; i<N_pheno ; i++)
			for (wsUint I=0,j=0,J=0 ; j<N_samp ; j++) {
				if (Ba_misCov[j]) continue;
				if (Ba_misPhCv[j])
					Ra_newY[i][I++] = W0;
				else
					Ra_newY[i][I++] = Ra_Y[i][J++];
			}

/**/		wsMat Ra_phi		= getFullCorMat(Cp_inpCorr);
			cMatrix	*Cp_phiMQ	= NULL;
			if (Ra_phi == NULL) {
				Cp_phiMQ	= new cIdtMatrix(N_nonCvMiss);
				B_indep		= 1;
			} else {
				Cp_phiMQ = new cStdMatrix(N_nonCvMiss, N_nonCvMiss);
				wsReal **Ra_phiMQ = Cp_phiMQ->get();
				for (wsUint i=0,I=0 ; i<N_samp ; i++) {
					if (Ba_misCov[i]) continue;
					for (wsUint j=0,J=0 ; j<=i ; j++) {
						if (Ba_misCov[j]) continue;

						Ra_phiMQ[J][I] = Ra_phiMQ[I][J] = Ra_phi[i][j];
						J++;
					}
					I++;
				}
			}
#ifdef MQLSvalidate
		exportSymMat("phi", Ra_phiMQ, N_nonCvMiss);
#endif

		/* Get phiInv */
		Cp_phiInv = &(Cp_phiMQ->inv());

		/* Substitute */
		N_anaSamp = N_nonCvMiss;

		/* Adjust M_Yf1, M_Yf2 */
		wsMat Ra_newYf1 = sseMatrix(N_pheno, N_anaSamp);
		for (wsUint x=0 ; x<N_pheno ; x++)
		for (wsUint i=0,I=0 ; i<N_samp ; i++) {
			if (Ba_misCov[i]) continue;
			Ra_newYf1[x][I] = Ba_misPhCv[i] ? W0 : Ra_a1[x][i];
			I++;
		}
		M_Yf1.init(N_pheno, N_anaSamp, Ra_newYf1);
		if (Ra_a2) {
			wsMat Ra_newYf2 = sseMatrix(N_pheno, N_anaSamp);
			for (wsUint x=0 ; x<N_pheno ; x++)
				for (wsUint i=0,I=0 ; i<N_samp ; i++) {
					if (Ba_misCov[i]) continue;
					Ra_newYf2[x][I] = Ba_misPhCv[i] ? W0 : Ra_a2[x][i];
					I++;
				}
			M_Yf2.init(N_pheno, N_anaSamp, Ra_newYf2);
		}
	} else {
		Ba_misCov = NULL;

		if (Ra_a1) {
			wsMat Ra_newYf1 = sseMatrix(N_pheno, N_anaSamp);
			for (wsUint x=0 ; x<N_pheno ; x++)
				for (wsUint i=0,I=0 ; i<N_samp ; i++) {
					if (Ba_misPhCv[i]) continue;
					Ra_newYf1[x][I] = Ra_a1[x][i];
					I++;
				}
			M_Yf1.init(N_pheno, N_anaSamp, Ra_newYf1);
			Ra_a1 = M_Yf1.get();
		}
		if (Ra_a2) {
			wsMat Ra_newYf2 = sseMatrix(N_pheno, N_anaSamp);
			for (wsUint x=0 ; x<N_pheno ; x++)
				for (wsUint i=0,I=0 ; i<N_samp ; i++) {
					if (Ba_misPhCv[i]) continue;
					Ra_newYf2[x][I] = Ra_a2[x][i];
					I++;
				}
			M_Yf2.init(N_pheno, N_anaSamp, Ra_newYf2);
			Ra_a2 = M_Yf1.get();
		}
	}


	/* (1'��^-11)^-1 */
	wsReal	R_1p1Sum	= W0;
	cVector	V_Apart		= Cp_phiInv->sumR(&R_1p1Sum);
	wsReal	*Ra_Apart	= V_Apart.get();
	delete Cp_phiInv;
	wsReal	R_1p1inv	= W1 / R_1p1Sum;
	V_Apart	*= R_1p1inv;

	/* Get Y'A */
	wsMat	Ra_TtA		= sseMatrix(N_pheno, N_anaSamp);
	wsMat	Ra_TtAf1	= Ra_a1 ? sseMatrix(N_pheno, N_anaSamp) : NULL;
	wsMat	Ra_TtAf2	= Ra_a2 ? sseMatrix(N_pheno, N_anaSamp) : NULL;
	for (wsUint i=0 ; i<N_pheno ; i++) {
		/* Get sum of ith phenotype */
		wsReal R_ySum = -sseVsum(Ra_Y[i], N_anaSamp);
		/* T'A_ij = (1p1)^-1 * (-T'_(i+) *  phi^-1_(+j) + y_ij) */
		for (wsUint j=0 ; j<N_anaSamp ; j++) {
			wsReal R_interm = R_ySum * Ra_Apart[j];

			Ra_TtA[i][j] = R_interm + Ra_Y[i][j];
			if (Ra_TtAf1) Ra_TtAf1[i][j] = R_interm + Ra_a1[i][j];
			if (Ra_TtAf2) Ra_TtAf2[i][j] = R_interm + Ra_a2[i][j];
		}
	}
	/* Insert */
	Ma_TtA.push_back(Ra_TtA);
	if (Ra_a1) Ma_TtA.push_back(Ra_TtAf1);
	if (Ra_a2) Ma_TtA.push_back(Ra_TtAf2);
#ifdef MQLSvalidate
	exportMatrix("TtA", Ra_TtA, N_pheno, N_anaSamp);
#endif

	// Set Ap = -Ap
//	sseVneg(Ra_Apart, N_sampMQ);

	/* I - 1(1'��^-11)^-11'��^-1
	 * Let z=(1'��^-11)^-1 then
	 *                    [ z .. z ]
	 * 1(1'��^-11)^-11' = [ ..     ] so (1z1')_i+ %*% ��^-1 = rowSum(��^-1) * z
	 *                    [ z .. z ] */
// /**/MAT_t Ra_A = sseMatrix(N_sampMQ, N_sampMQ);
// 	for (wsUint i=0 ; i<N_sampMQ ; i++) {
// 		memcpy(Ra_A[i], Ra_Apart, sizeof(wsReal)*N_sampMQ);
// 		Ra_A[i][i] += W1;
// 	}
	V_Apart.rem();

	/* Get T^T %*% A [Q*N][N*N] */
///**/Ra_TtA = sseMpM(Ra_Y, N_pheno, N_sampMQ, Ra_A, N_sampMQ, N_sampMQ);
//	sseUnmat(Ra_A, N_sampMQ);

	/* Get (T'A��A'T)^-1 for all A's */
	FOREACH (vector<wsMat>::iterator, Ma_TtA, i) {
		wsMat	Ra_TtA		= *i;
/**/	wsSym	Ra_TAPAT	= B_indep ?
			sym_sseMpMt(Ra_TtA, N_pheno, N_anaSamp)
				: sseMSMt(Ra_TtA, N_pheno, N_anaSamp, Cp_phiMQ->get(),
					N_anaSamp);
		cSymMatrix	M_TAPAT(Ra_TAPAT, N_pheno);
		cSymMatrix*	Mp_TAPTAi = &(M_TAPAT.inv());
		Ma_TtAphiAtT_inv.push_back(Mp_TAPTAi);
	}

	/* Remove reduced phi matrix */
	delete Cp_phiMQ;
}

cMFqlsAnalysis::~cMFqlsAnalysis()
{
	FOREACH (vector<cSymMatrix*>::iterator, Ma_TtAphiAtT_inv, i)
		delete *i;
	DEALLOC(Ba_misPhCv);
	DEALLOC(Ba_misCov);
	FOREACH (vector<wsMat>::iterator, Ma_TtA, i)
		deallocMatrix(*i, N_pheno, (void *)1);
}

int _doTest(xMqlsTest *Xp_t)
{
	cIO*		Cp_IO		= Xp_t->Cp_IO;
	vVariant&	Xa_mkr		= Cp_IO->getVariant();
	wsUintCst		N_sample	= Cp_IO->sizeSample();
	char**		Na_data		= Cp_IO->getGenotype();
	/* Start and end of gene-set */
	wsUint		N_start		= 0;
	wsUint		N_end		= 0;
	vInt	Xa_set;
	wsUint	i, j;
#ifdef MQLSvalidate
	char	S_buf[128];
#endif

	/* Remove MAF=0 */

	/* At first, set every variants to exclude */
	vInt Xa_inc;
	Xa_inc.resize(Xp_t->Xa_snps.size(), 0);
	//FOREACH (vInt_it, it->second, iit) Xa_inc.push_back(0);

	/* See the genotype of every samples, and marking it to include
	 * if it have genotype that not 0 */
	for (i=0 ; i<N_sample ; i++) {
		if (Xp_t->Ba_misPheno[i]) continue;

		j = 0;
		FOREACHDO (vInt_it, Xp_t->Xa_snps, iit, j++) {
			if (isAvailable(Na_data[i][*iit]) && Na_data[i][*iit]>0)
				Xa_inc[j] = 1;
		}
	}

	/* Make final gene-set with variants that mark to be included */
	i = j = 0;
	FOREACHDO (vInt_it, Xp_t->Xa_snps, iit, j++) {
		if (Xa_inc[j]) {
			if (i == 0) {
				N_start = N_end = Xa_mkr[*iit].pos;
			} else {
				if (Xa_mkr[*iit].pos < N_start)
					N_start = Xa_mkr[*iit].pos;
				if (Xa_mkr[*iit].pos > N_end)
					N_end = Xa_mkr[*iit].pos;
			}
			Xa_set.push_back(*iit);
			i++;
		}
	}
	wsUint	N_gs		= (wsUint)Xa_set.size();
	/* No available variants */
	if (N_gs == 0) {
		for (wsUint i=0 ; i<MQLS_NTEST ; i++) {
			Xp_t->X_res.R_p[i] = WISARD_NAN;
			Xp_t->X_res.R_t[i] = WISARD_NAN;
		}
		return 0;
	}

	/* Make subset of genotypes */
/**/wsReal **Ra_Xt = sseMatrix(N_gs, Xp_t->N_anaSamp); /* [M*N] since it is T */
	j = 0;
	FOREACHDO (vInt_it, Xa_set, it, j++) {
		//wsUint I = 0;

		setGeno(Cp_IO, *it, Ra_Xt[j], Xp_t->Ba_misPheno, 1);
// 		for (i=0 ; i<N_sample ; i++) {
// 			if (Xp_t->Ba_misPheno[i]) continue;
// 
// 			if (isMissing(Na_data[i][*it]))
// 				Ra_Xt[j][I] = Xa_mkr[*it].maf * W2;
// 			else
// 				Ra_Xt[j][I] = (wsReal)Na_data[i][*it];
// 
// 			I++;
// 		}
	}
#ifdef MQLSvalidate
	sprintf(S_buf, "%s.Xt", Xp_t->Sp_name);
	exportMatrix(S_buf, Ra_Xt, N_gs, Xp_t->N_anaSamp);
#endif

	/* Get AX [N*P] <= Calcelled
	 *               [  A' ]                     1
	 * Since A = I - [ ... ] and A'_ij = - ------------- * sum(��^-1_j+)
	 *               [  A' ]               sum(��^-1_++)
	 *
	 * AX_ij = -phi * sum(��^-1_j+ * X'_j+) + X_ji             1
	 *                                        where phi = - -------------
	 *                                                      sum(��^-1_++)
	 */

	/* Precompute M_psiInv if required */
	cSymMatrix	*Mp_psiInv = NULL;
	if (N_gs != 1) {
		wsSym		Ra_psi = sym_sseCovRow(Ra_Xt, N_gs, Xp_t->N_anaSamp);
		cSymMatrix	M_psi(Ra_psi, N_gs);
		Mp_psiInv = &(M_psi.inv(Xp_t->X_res.B_ginv));
#ifdef MQLSvalidate
		sprintf(S_buf, "%s.psi", Xp_t->Sp_name);
		exportSymMat(S_buf, Ra_psi, N_gs);
#endif
		if (Mp_psiInv->get() == NULL)
			LOG("Failed to test [%s]\n", Xp_t->Sp_name);
	}

	/* Get cov(X) */
	for (wsUint X=0 ; X<Xp_t->Ma_TtA.size() ; X++) {
/**/	wsSym	Ra_varSinv	= NULL;
/**/	wsReal*	Ra_vecTtAX = NULL;
		wsUint	N_pheno		= Xp_t->N_pheno;
		wsMat	Ra_TtA		= Xp_t->Ma_TtA[X];
		wsSym	Ra_TAPATi	= Xp_t->Ma_TtAphiAtT_inv[X]->get();

		if (N_gs != 1) {
			/* Calc TtA %*% X -> [Q*M] */
/**/		wsMat Ra_TtAX		= NULL;
			Ra_TtAX = sseMpMt(Ra_TtA, N_pheno, Xp_t->N_anaSamp,
				Ra_Xt, N_gs, Xp_t->N_anaSamp);
			sseUnmat(Ra_Xt, N_gs);

			/* Do vec(T'AX)' (�� KRONECKER (T'A��A'T)^-1) vec(T'AX) */
			Ra_varSinv = sseSkS(Mp_psiInv->get(), N_gs, Ra_TAPATi, N_pheno);

			sseMalloc(Ra_vecTtAX, wsReal, N_pheno*N_gs);
			/* vec(T'AX) */
			for (wsUint i=0 ; i<N_pheno ; i++) {
				wsUint k = i;
				for (wsUint j=0 ; j<N_gs ; j++,k+=N_pheno)
					Ra_vecTtAX[k] = Ra_TtAX[i][j];
			}
			sseUnmat(Ra_TtAX, N_pheno);
		} else {
			Ra_vecTtAX = sseMpV(Ra_TtA, N_pheno, Xp_t->N_anaSamp, Ra_Xt[0]);
			Ra_varSinv = sseSymMat(N_pheno);
			sseSpC((wsSymCst)Ra_TAPATi, W1/ sseVvar(Ra_Xt[0], Xp_t->N_anaSamp),
				Ra_varSinv, N_pheno);
		}

		/* Get statistics */
		if (Ra_vecTtAX) {
			Xp_t->X_res.R_t[X] = sseVpSpV(Ra_vecTtAX, Ra_varSinv, N_pheno*N_gs);
			sseFree(Ra_vecTtAX);
			deallocMatrix(Ra_varSinv, N_pheno*N_gs, (void *)1);

			Xp_t->X_res.R_p[X] = (wsReal)PVchisq(Xp_t->X_res.R_t[X],
				(double)(N_pheno*N_gs));
		} else {
			Xp_t->X_res.R_t[X] = Xp_t->X_res.R_p[X] = WISARD_NAN;
		}
	}
	/* Remove psiInv */
	if (Mp_psiInv) delete Mp_psiInv;
	
	return N_gs;	
}

int _doTest(cIO *Cp_IO, wsStrCst Sp_name, vInt& Xa_snps, xMqlsRes &X_res,
	vector<wsMat>& Ma_TtA, vector<cSymMatrix*>& Ma_TtAphiAtT_inv,
	const char* Ba_misPhCv, wsUint N_anaSamp)
{
	vVariant&	Xa_mkr		= Cp_IO->getVariant();
	wsUint		N_pheno		= Cp_IO->sizePheno();
	wsUintCst		N_sample	= Cp_IO->sizeSample();
	char**		Na_data		= Cp_IO->getGenotype();
	/* Start and end of gene-set */
	wsUint		N_start		= 0;
	wsUint		N_end		= 0;
	vInt		Xa_set;
	wsUint		i, j;
#ifdef MQLSvalidate
	char	S_buf[256];
#endif

	/* Remove MAF=0 */

	/* At first, set every variants to exclude */
	vInt Xa_inc;
	Xa_inc.resize(Xa_snps.size(), 0);
	//FOREACH (vInt_it, it->second, iit) Xa_inc.push_back(0);

	/* See the genotype of every samples, and marking it to include
	 * if it have genotype that not 0 */
	for (i=0 ; i<N_sample ; i++) {
		if (Ba_misPhCv[i]) continue;

		j = 0;
		FOREACHDO (vInt_it, Xa_snps, iit, j++) {
			if (isAvailable(Na_data[i][*iit]) && Na_data[i][*iit]>0)
				Xa_inc[j] = 1;
		}
	}

	/* Make final gene-set with variants that mark to be included */
	i = j = 0;
	FOREACHDO (vInt_it, Xa_snps, iit, j++) {
		if (Xa_inc[j]) {
			if (i == 0) {
				N_start = N_end = Xa_mkr[*iit].pos;
			} else {
				if (Xa_mkr[*iit].pos < N_start)
					N_start = Xa_mkr[*iit].pos;
				if (Xa_mkr[*iit].pos > N_end)
					N_end = Xa_mkr[*iit].pos;
			}
			Xa_set.push_back(*iit);
			i++;
		}
	}
	wsUint	N_gs	= (wsUint)Xa_set.size();
	if (N_gs == 0) {
		for (wsUint i=0 ; i<MQLS_NTEST ; i++) {
			X_res.R_t[i] = WISARD_NAN;
			X_res.R_p[i] = WISARD_NAN;
		}
		return 0;
	}

	/* Make subset of genotypes */
/**/wsReal **Ra_Xt = sseMatrix(N_gs, N_anaSamp); /* [M*N] since it is T */
	j = 0;
	FOREACHDO (vInt_it, Xa_set, it, j++) {
//		wsUint I = 0;

		setGeno(Cp_IO, *it, Ra_Xt[j], Ba_misPhCv, 1);

// 		for (i=0 ; i<N_sample ; i++) {
// 			if (Ba_misPhCv[i]) continue;
// 
// 			if (isMissing(Na_data[i][*it]))
// 				Ra_Xt[j][I] = Xa_mkr[*it].maf * W2;
// 			else
// 				Ra_Xt[j][I] = (wsReal)Na_data[i][*it];
// 
// 			I++;
// 		}
	}
	cStdMatrix M_Xt(N_gs, N_anaSamp, Ra_Xt, MATDEL_NONE);
#ifdef MQLSvalidate
	sprintf(S_buf, "%s.Xt", Sp_name);
	exportMatrix(S_buf, Ra_Xt, N_gs, N_anaSamp);
#endif

	/* Get AX [N*P] <= Calcelled
	 *               [  A' ]                     1
	 * Since A = I - [ ... ] and A'_ij = - ------------- * sum(��^-1_j+)
	 *               [  A' ]               sum(��^-1_++)
	 *
	 * AX_ij = -phi * sum(��^-1_j+ * X'_j+) + X_ji             1
	 *                                        where phi = - -------------
	 *                                                      sum(��^-1_++)
	 */

	/* Precompute M_psiInv if required */
	cSymMatrix	*Mp_psiInv = NULL;
	if (N_gs != 1) {
		wsSym		Ra_psi = sym_sseCovRow(Ra_Xt, N_gs, N_anaSamp);
		cSymMatrix	M_psi(Ra_psi, N_gs);
		Mp_psiInv = &(M_psi.inv(X_res.B_ginv));
#ifdef MQLSvalidate
		sprintf(S_buf, "%s.psi", Xp_t->Sp_name);
		exportSymMat(S_buf, Ra_psi, N_gs);
#endif
		if (Mp_psiInv->get() == NULL)
			LOG("Failed to test [%s]\n", Sp_name);
	}

	/* Get cov(X) */
	if (N_gs == 1 || Mp_psiInv->get()) for (wsUint X=0 ; X<Ma_TtA.size() ; X++) {
///**/	SYM_t	Ra_varSinv	= NULL;
///**/	wsReal*	Ra_vecTtAX = NULL;
		cVector	V_vecTtAX;
		wsMat	Ra_TtA		= Ma_TtA[X];
		cSymMatrix&	M_TAPATi	= *(Ma_TtAphiAtT_inv[X]);
		cSymMatrix*	Mp_varSinv	= NULL;
//		SYM_t	Ra_TAPATi	= Ma_TtAphiAtT_inv[X]->get();
		cStdMatrix	M_TtA(N_pheno, N_anaSamp, Ra_TtA, MATDEL_NONE);

		if (N_gs != 1) {
			/* Calc TtA %*% X -> [Q*M] */
///**/		MAT_t Ra_TtAX		= NULL;
			cStdMatrix M_TtAX	= M_TtA.Mt(M_Xt);
// 			Ra_TtAX = sseMpMt(Ra_TtA, N_pheno, N_anaSamp,
// 				Ra_Xt, N_gs, N_anaSamp);
			//sseUnmat(Ra_Xt, N_gs);

			/* Do vec(T'AX)' (�� KRONECKER (T'A��A'T)^-1) vec(T'AX) */
			Mp_varSinv = &(Mp_psiInv->krnk(M_TAPATi));
//			Ra_varSinv = sseSkS(Mp_psiInv->get(), N_gs, Ra_TAPATi, N_pheno);

//			sseMalloc(Ra_vecTtAX, wsReal, N_pheno*N_gs);
			/* vec(T'AX) */
// 			for (wsUint i=0 ; i<N_pheno ; i++) {
// 				wsUint k = i;
// 				for (wsUint j=0 ; j<N_gs ; j++,k+=N_pheno)
// 					Ra_vecTtAX[k] = Ra_TtAX[i][j];
// 			}
			V_vecTtAX = M_TtAX.vec();
//			sseUnmat(Ra_TtAX, N_pheno);
		} else {
			cVector X = M_Xt.r2v_ptr(0);
			X.setDontDealloc();
			V_vecTtAX = M_TtA * X;
//			Ra_vecTtAX = sseMpV(Ra_TtA, N_pheno, N_anaSamp, Ra_Xt[0]);
//			Ra_varSinv = sseSymMat(N_pheno);
			Mp_varSinv = &(M_TAPATi * (W1/ sseVvar(Ra_Xt[0], N_anaSamp)));
// 			sseSpC(Ra_TAPATi, W1/ sseVvar(Ra_Xt[0], N_anaSamp),
// 				Ra_varSinv, N_pheno);
		}

		/* Get statistics */
//		if (Ra_vecTtAX) {
			//X_res.R_t[X] = sseVpSpV(Ra_vecTtAX, Ra_varSinv, N_pheno*N_gs);
			X_res.R_t[X] = V_vecTtAX.qf(*Mp_varSinv);
			//sseFree(Ra_vecTtAX);
			//deallocMatrix(Ra_varSinv, N_pheno*N_gs, (void *)1);

			X_res.R_p[X] = (wsReal)PVchisq(X_res.R_t[X],
				(double)(N_pheno*N_gs));
// 		} else {
// 			X_res.R_t[X] = X_res.R_p[X] = WISARD_NAN;
// 		}
		delete Mp_varSinv;
	} else for (wsUint X=0 ; X<Ma_TtA.size() ; X++) {
		X_res.R_t[X] = WISARD_NAN;
		X_res.R_p[X] = WISARD_NAN;
	}
	sseUnmat(Ra_Xt, N_gs);
	/* Remove psiInv */
	if (Mp_psiInv) delete Mp_psiInv;

	return N_gs;
}

void prepMfqlsOutput(cIO* Cp_IO, cTimer& t, xMqlsRes& X_res, char *S_time,
	char* S_fqls)
{
	if (OPT_ENABLED(time))
		sprintf(S_time, "%s	", t.getReadable());

	if (OPT_ENABLED(fqls)) {
		char S_buf1[128] = { 0, };
		char S_buf2[128] = { 0, };

		if (NA(X_res.R_t[1]) || NA(X_res.R_p[1]))
			strcpy(S_buf1, "	NA	NA");
		else
			sprintf(S_buf1, "	%g	%g", X_res.R_t[1], X_res.R_p[1]);
		if (Cp_IO->isProbandAssigned()) {
			if (NA(X_res.R_t[2]) || NA(X_res.R_p[2]))
				strcpy(S_buf2, "	NA	NA");
			else
				sprintf(S_buf2, "	%g	%g", X_res.R_t[2],
				X_res.R_p[2]);
		}
		sprintf(S_fqls, "%s%s", S_buf1, S_buf2);
	}
}

void prepMfqlsOutput(cIO* Cp_IO, wsMat Ra_res, wsUint N_idx, char *S_time,
	char* S_fqls)
{
	if (OPT_ENABLED(time))
		sprintf(S_time, "%s	", cTimer::fmt(Ra_res[0][N_idx]));

	if (OPT_ENABLED(fqls)) {
		char S_buf1[128] = { 0, };
		char S_buf2[128] = { 0, };

		if (NA(Ra_res[2][N_idx]) || NA(Ra_res[2+MQLS_NTEST][N_idx]))
			strcpy(S_buf1, "	NA	NA");
		else
			sprintf(S_buf1, "	%g	%g", Ra_res[2][N_idx], Ra_res[2+MQLS_NTEST][N_idx]);
		if (Cp_IO->isProbandAssigned()) {
			if (NA(Ra_res[3][N_idx]) || NA(Ra_res[3+MQLS_NTEST][N_idx]))
				strcpy(S_buf2, "	NA	NA");
			else
				sprintf(S_buf2, "	%g	%g", Ra_res[3][N_idx],
					Ra_res[3+MQLS_NTEST][N_idx]);
		}
		sprintf(S_fqls, "%s%s", S_buf1, S_buf2);
	}
}

void cMFqlsAnalysis::run()
{
/**/cExporter*	Cp_mqls = cExporter::summon("extended.qls.res");
	LOGoutput("Result of extended QLS test is exported to [%s.extended.qls.res]\n", /* CONFIRMED 140109 */
		OPT_STRING(out));
	char		S_time[256] = { 0, };
	char		S_fqls[256] = { 0, };

	if (OPT_ENABLED(time))
		sprintf(S_time, "TIME	");
	if (OPT_ENABLED(fqls)) {
		if (Cp_IO->isProbandAssigned())
			sprintf(S_fqls, "	STAT_FQLS1	P_FQLS1	STAT_FQLS2	P_FQLS2");
		else
			sprintf(S_fqls, "	STAT_FQLS1	P_FQLS1");
	}

	if (IS_ASSIGNED(mqlsconsec)) {
		wsUint N_cs = OPT_NUMBER(mqlsconsec);
		LOG("Consecutive [%d] variants will be tested...\n", N_cs);

		/* Print header */
		Cp_mqls->fmt("PAIR	%sGINV	STAT_MFQLS	P_MFQLS%s\n", S_time, S_fqls);

		/* Concatenate 2 snps into 1 */
		cTimer		t;
		vVariant	&Xa_snp	= getIO()->getVariant();
		wsUint		N_vrt	= getIO()->sizeVariant();
		if (Cp_mm) {
			wsUint *Na_chr = NULL;
			/* Mode 1 : Use 'actually' consecutive region */
			wsUint **Na_pos = Cp_mm->getVariantPosMap(&Na_chr);

			for (wsUint chr=0 ; chr<=NCHR_SPECIES ; chr++) {
				wsUint *Np_pos = Na_pos[chr];
				/* Skip if there is not enough number of variants in the chromosome
				 * to perform N_cs consecutive variants */
				if (Na_chr[chr] < N_cs) continue;
				wsUint N_loop = Na_chr[chr] - N_cs + 1;

				for (wsUint i=0 ; i<N_loop ; i++) {
					xMqlsRes	X_res = { { W0, }, { W0, }, 0 };

					/* Make pairing name */
					wsUint	N_pos = Np_pos[i];
					vInt	Xa_snps;
					string	S_pName(Xa_snp[N_pos].name);
					int		N_chr = Xa_snp[N_pos].chr;
					char	B_test = 1;
					Xa_snps.push_back(N_pos);
					for (wsUint j=1 ; j<N_cs ; j++) {
						xVariant &X_s = Xa_snp[Np_pos[i+j]];

						/* Make do not test if chr are disagree */
						if (X_s.chr != N_chr) {
							B_test = 0;
							break;
						}

						S_pName += ",";
						S_pName += X_s.name;
						Xa_snps.push_back(Np_pos[i+j]);
					}
					/* Do not test if B_test is disabled */
					if (B_test == 0) continue;

					/* Do test */
					t.start();
					ONETOOL::_doTest(getIO(), S_pName.c_str(), Xa_snps, X_res,
						Ma_TtA, Ma_TtAphiAtT_inv, Ba_misPhCv, N_anaSamp);
					prepMfqlsOutput(getIO(), t, X_res, S_time, S_fqls);

					/* Make output */
					if (NA(X_res.R_t[0]) || NA(X_res.R_p[0])) {
						/* --remna */
						if (OPT_ENABLED(remna)) continue;
						Cp_mqls->fmt("%s	%s%d	NA	NA%s\n", S_pName.c_str(),
							S_time, X_res.B_ginv, S_fqls);
					} else {
						/* --pvalrange */
						if (IS_ASSIGNED(pvalrange) &&
							!isInRange(OPT_RANGE(pvalrange), X_res.R_p[0]))
							continue;
						Cp_mqls->fmt("%s	%s%d	%g	%g%s\n", S_pName.c_str(),
							S_time, X_res.B_ginv, X_res.R_t[0], X_res.R_p[0],
							S_fqls);
					}

					if ((i%1000) == 0)
						notice("[Chr %d] %d/%d pairs tested...\r", chr, i, Na_chr[chr]);
				}
			}
		} else {
			wsUint	N_pair	= N_vrt - N_cs + 1; /* Since it is "pairing" */
			for (wsUint i=0 ; i<N_pair ; i++) {
				xMqlsRes	X_res = { { W0, }, { W0, }, 0 };
				/* Make pairing name */
				vInt		Xa_snps;
				string		S_pName(Xa_snp[i].name);
				int			N_chr	= Xa_snp[i].chr;
				char		B_test	= 1;
				Xa_snps.push_back(i);
				for (wsUint j=1 ; j<N_cs ; j++) {
					xVariant &X_s = Xa_snp[i+j];

					/* Make do not test if chr are disagree */
					if (X_s.chr != N_chr) {
						B_test = 0;
						break;
					}

					S_pName += ",";
					S_pName += Xa_snp[i+j].name;
					Xa_snps.push_back(i+j);
				}
				/* Do not test if B_test is disabled */
				if (B_test == 0) continue;

				/* Do test */
				t.start();
				wsUint N_snp = ONETOOL::_doTest(getIO(), S_pName.c_str(), Xa_snps,
					X_res, Ma_TtA, Ma_TtAphiAtT_inv, Ba_misPhCv, N_anaSamp);
				prepMfqlsOutput(getIO(), t, X_res, S_time, S_fqls);

				if (NA(X_res.R_t[0]) || NA(X_res.R_p[0])) {
					/* --remna */
					if (OPT_ENABLED(remna)) continue;
					Cp_mqls->fmt("%s	%s%d	NA	NA%s\n", S_pName.c_str(),
						S_time, X_res.B_ginv, S_fqls);
				} else {
					/* --pvalrange */
					if (IS_ASSIGNED(pvalrange) &&
						!isInRange(OPT_RANGE(pvalrange), X_res.R_p[0])) continue;
					Cp_mqls->fmt("%s	%s%d	%g	%g%s\n", S_pName.c_str(),
						S_time, X_res.B_ginv, X_res.R_t[0], X_res.R_p[0],
						S_fqls);
				}

				if ((N_pair%1000) == 0 && N_snp)
					notice("%d/%d pairs tested...\r", i, N_pair);
			}
			LOG("%d/%d pairs tested...\n", N_pair, N_pair);
		}
	} else if (SET_AVAIL()) {
		/* Print header */
		Cp_mqls->fmt("GENE	SZSET	NVAR	%sGINV	STAT_MFQLS	P_MFQLS%s\n",
			S_time, S_fqls);

		/* If it is gene-based */
		mGeneDef&	Xa_set	= Cp_anaGS->getGeneDef();
		wsUint		i		= 0;
		wsUint		N_gene	= (wsUint)Xa_set.size();
		cTimer		t, tAll;

		tAll.start();
		if (OPT_NUMBER(thread) == 1) FOREACHDO (mGeneDef_it, Xa_set, it, i++) {
			xMqlsRes	X_res = { { W0, }, { W0, }, 0 };

			t.start();
			wsUint N_snp = ONETOOL::_doTest(getIO(), it->first.c_str(), it->second,
				X_res, Ma_TtA, Ma_TtAphiAtT_inv, Ba_misPhCv, N_anaSamp);
			prepMfqlsOutput(getIO(), t, X_res, S_time, S_fqls);

			if (NA(X_res.R_t[0]) || NA(X_res.R_p[0])) {
				if (OPT_ENABLED(remna)) continue;
				Cp_mqls->fmt("%s	%d	%d	%s%d	NA	NA%s\n",
					it->first.c_str(), it->second.size(), N_snp, S_time,
					X_res.B_ginv, S_fqls);
			} else {
				/* --pvalrange */
				if (IS_ASSIGNED(pvalrange) &&
					!isInRange(OPT_RANGE(pvalrange), X_res.R_p[0])) continue;
				Cp_mqls->fmt("%s	%d	%d	%s%d	%g	%g%s\n",
					it->first.c_str(), it->second.size(), N_snp, S_time,
					X_res.B_ginv, X_res.R_t[0], X_res.R_p[0], S_fqls);
			}

			notice("[%d/%d] gene-sets tested...\r", i, Xa_set.size());
		} else {
			/* Prepare buffer */
			wsMat	Ra_rets = sseMatrix(1+MQLS_NTEST*2+1, N_gene);

			/* Do regression */
			xAnaSetThread X_at = { getIO(), Cp_anaGS, this };
			WORKER().run(thrMFQLSgene, forGene_equal, &X_at, Ra_rets,
				sizeof(int)*3);

			FOREACHDO (mGeneDef_it, Xa_set, it, i++) {
				xMqlsRes	X_res = { { W0, }, { W0, }, 0 };
				prepMfqlsOutput(getIO(), Ra_rets, i, S_time, S_fqls);

				wsReal R_t = Ra_rets[1][i];
				wsReal R_p = Ra_rets[1+MQLS_NTEST][i];
				if (NA(R_t) || NA(R_p)) {
					if (OPT_ENABLED(remna)) continue;
					Cp_mqls->fmt("%s	%d	%d	%s%d	NA	NA%s\n",
						it->first.c_str(), it->second.size(),
						(wsUint)Ra_rets[MQLS_NTEST*2+1][i], S_time,
						X_res.B_ginv, S_fqls);
				} else {
					/* --pvalrange */
					if (IS_ASSIGNED(pvalrange) &&
						!isInRange(OPT_RANGE(pvalrange), X_res.R_p[0])) continue;
					Cp_mqls->fmt("%s	%d	%d	%s%d	%g	%g%s\n",
						it->first.c_str(), it->second.size(),
						(wsUint)Ra_rets[MQLS_NTEST*2+1][i], S_time,
						X_res.B_ginv, R_t, R_p, S_fqls);
				}
			}

			sseUnmat(Ra_rets, 1+MQLS_NTEST*2+1);
		}
		LOG("[%d/%d] genes tested in [%s]\n", i, Xa_set.size(),
			tAll.getReadable());
	} else {
		/* Print header regarding --annogene,--time */
		headerVariant(Cp_mqls);
		Cp_mqls->fmt("%s	STAT	P_MFQLS%s\n", S_time, S_fqls);

		/* If it is variant-based */
		cTimer		t, tAll;
		vVariant	Xv_vrt			= Cp_IO->getVariant();
		wsUint		N_vrt			= (wsUint)Xv_vrt.size();
		wsUint		i				= 0;
		char		S_bufStat[256]	= { 0, };

		tAll.start();
		if (OPT_NUMBER(thread) == 1) FOREACHDO (vVariant_it, Xv_vrt, it, i++) {
			vInt		X_snp;
			xMqlsRes	X_res = { { W0, }, { W0, }, 0 };
			X_snp.push_back(i);

			t.start();
			ONETOOL::_doTest(getIO(), it->name, X_snp, X_res, Ma_TtA,
				Ma_TtAphiAtT_inv, Ba_misPhCv, N_anaSamp);
			prepMfqlsOutput(getIO(), t, X_res, S_time, S_fqls);

			if (NA(X_res.R_t[0]) || NA(X_res.R_p[0])) {
				if (OPT_ENABLED(remna)) continue;
				strcpy(S_bufStat, "	NA	NA");
			} else {
				/* --pvalrange */
				if (IS_ASSIGNED(pvalrange) &&
					!isInRange(OPT_RANGE(pvalrange), X_res.R_p[0])) continue;
				sprintf(S_bufStat, "	%g	%g", X_res.R_t[0], X_res.R_p[0]);
			}

			entryVariant(Cp_mqls, *it);
			Cp_mqls->fmt("%s%s%s\n", S_time, S_bufStat, S_fqls);

			if ((i%100) == 0)
				notice("[%d/%d] variants tested...\r", i, Xv_vrt.size());
		} else {
			/* Prepare buffer */
			wsMat Ra_rets = sseMatrix(1+MQLS_NTEST*2, N_vrt);

			/* Do regression */
			xAnaThread X_at = { getIO(), this };
			WORKER().run(thrMFQLSvariant, forVariant_equal, &X_at, Ra_rets,
				sizeof(int)*3);

			LOOP (i, N_vrt) {
				prepMfqlsOutput(getIO(), Ra_rets, i, S_time, S_fqls);

				wsReal R_t = Ra_rets[1][i];
				wsReal R_p = Ra_rets[1+MQLS_NTEST][i];
				if (NA(R_t) || NA(R_p)) {
					if (OPT_ENABLED(remna)) continue;
					strcpy(S_bufStat, "	NA	NA");
				} else {
					/* --pvalrange */
					if (IS_ASSIGNED(pvalrange) &&
						!isInRange(OPT_RANGE(pvalrange), R_p)) continue;
					sprintf(S_bufStat, "	%g	%g", R_t, R_p);
				}

				entryVariant(Cp_mqls, Xv_vrt[i]);
				Cp_mqls->fmt("%s%s%s\n", S_time, S_bufStat, S_fqls);
			}

			sseUnmat(Ra_rets, 1+MQLS_NTEST*2);
		}
		LOG("[%d/%d] variants tested in [%s]\n", Xv_vrt.size(), Xv_vrt.size(),
			tAll.getReadable());
	}

	delete Cp_mqls;
}

#endif

} // End namespace ONETOOL
