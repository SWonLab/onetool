#include "analyses/qls.h"
#include "analyses/emai.h"
#include "utils/util.h"
#include "utils/matrix.h"
#include "analyses/annot.h"
#include "analyses/ppp.h"
#include "utils/stat.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cQlsAnalysis::cQlsAnalysis(cIO *Cp_inpIO, cPPPAnalysisV2 *Cp_inpAnaPPP,
	cAnalysis *Cp_inpAnaCorr) : cAnalysis(Cp_inpIO)
{
	/* Initialize member */
	Cp_anaCorr	= Cp_inpAnaCorr;
	Cp_anaPPP	= Cp_inpAnaPPP;

	/* Option check */
	if (!IS_ASSIGNED(blup))
		LOG("--blup is not assigned, assume it is integrated into phenotype\n");

	if (!IS_ASSIGNED(genemiss)) {
		OPTION().assign("genemiss", OPTION().getDefVal("genemiss"));
		OPTION().FORCE_OPT_REAL(genemiss);
	}

	/* Make phenotype-available index */
	Ba_isInc = Cp_IO->getAvailPheno(&N_incSamp, &Ra_availY);
	LOG("%d samples will include to QLS analysis\n", N_incSamp);
}

cQlsAnalysis::~cQlsAnalysis()
{
	DEALLOC(Ba_isInc);
	sseFree(Ra_availY);
}

void cQlsAnalysis::_doTest(wsUint i, char **Na_data,
	cMatrix *Mp_phi, wsRealCst *Ra_ppp, cVector& V_yV, wsReal R_yVy,
	double &R_Tqls, double &R_Pqls)
{
	xMaf&		X_maf	= getIO()->getMAF()[i];
	wsUint		N_samp	= getIO()->sizeSample();
	wsReal*		Ra_curX	= sseVector(N_incSamp);

	if (OPT_ENABLED(avail)) {
		cVector	V_curY(N_incSamp);
		wsReal*	Ra_curY		= V_curY.get();

		/* Copy Ba_isInc */
		char*	Ba_curIsInc	= NULL;
		wsAlloc(Ba_curIsInc, char, N_samp);
		memcpy(Ba_curIsInc, Ba_isInc, sizeof(char)*N_samp);

		/* Which are missing? */
		wsUint N_curSamp = 0;
		for (wsUint j=0,K=0 ; j<N_samp ; j++) {
			if (!Ba_isInc[j]) continue;

			char N_geno = Na_data[j][i];
			if (isMissing(N_geno))
				Ba_curIsInc[j] = 1;
			else {
				Ra_curY[N_curSamp] = Ra_availY[K];
				Ra_curX[N_curSamp++] = N_geno;
			}
			K++;
		}

		/* Get phi */
		cMatrix&	M_cPhi	= Mp_phi->subset(Ba_curIsInc, 1);
		cMatrix&	M_cPinv	= M_cPhi.inv();
			
		/* STEP 1 : Get (1' %*% Rinv %*% 1)^-1 */
		wsReal		R_c1p1	= M_cPinv.sum();
		wsMat		Ra_cPinv= M_cPinv.get();

		/* STEP 2 : Get rowsum */
		cVector		V_cR	= M_cPinv.sumR();
		wsReal		*Ra_cR	= V_cR.get();

		/* STEP 3 : Get R^-1 - V */
		for (i=0 ; i<N_incSamp ; i++)
			sseVsV(Ra_cPinv[i], Ra_cR, Ra_cPinv[i], i+1, Ra_cR[i]/R_c1p1);
		V_cR.rem();

		/* STEP 4 : Get divider part : yV <- (Y-mu)'(R^-1-V) */
		cVector		V_yV;
		if (OPT_ENABLED(indep))
			V_yV = V_curY * (cIdtMatrix &)M_cPinv;
		else
			V_yV = V_curY * (cSymMatrix &)M_cPinv;

		/* STEP 5 : Get denominator part : (Y-mu)'(R^-1-V)(Y-mu)
			*  - can be optimized to yV %*% (Y-mu) */
		wsReal		R_cYVY	= V_yV.sum(V_curY);

		/* Get the full denominator = sqrt(sig2*yVy) */
		wsReal R_denom = sqrt(W2*Ra_ppp[i]*(W1 - Ra_ppp[i])*R_cYVY);

		/* Y already BLUP enabled */
		R_Tqls = sseVV(V_yV.get(), N_curSamp, Ra_curX) / R_denom;
		R_Pqls = norprobP(R_Tqls);
	} else if (!IS_ASSIGNED(genemiss) || (N_samp-(X_maf.N_allele>>1))/(wsReal)N_samp >= OPT_REAL(genemiss)) {
		/* Build ith data */
 		setGeno(getIO(), i, Ra_curX, Ba_isInc);

		/* Get the full denominator = sqrt(sig2*yVy) */
		wsReal R_denom = sqrt(W2*Ra_ppp[i]*(W1 - Ra_ppp[i])*R_yVy);

		/* Y already BLUP enabled */
		cVector V_curX(Ra_curX, N_incSamp, 1);
		R_Tqls = V_yV.sum(V_curX) / R_denom;
		R_Pqls = norprobP(R_Tqls);
	}

	sseFree(Ra_curX);
}

void cQlsAnalysis::run()
{
	if (N_incSamp == 0) {
		LOG("No sample in analysis, cannot perform\n");
		return;
	}

	vVariant	&Xa_snp		= Cp_IO->getVariant();
	wsUintCst	N_sample	= Cp_IO->sizeSample();
	wsRealCst	*Ra_ppp		= Cp_anaPPP->getPPP();
	wsReal	*Ra_curX	= NULL; /* CHECKED */
	char	**Na_data	= Cp_IO->getGenotype();
	wsReal	**Ra_oCorr	= NULL; /* CHECKED */
//	wsReal	**Ra_Rinv	= NULL; /* CHECKED */
	wsReal	**Ra_corr	= NULL; /* CHECKED */
	wsUint	i;

	/* Prepare exporter */
	vStr		Sv_headers;
	char		S_fmts[16];
	wsUint		N_elem = 0;
	if (IS_ASSIGNED(time)) {
		Sv_headers.push_back("TIME");
		S_fmts[N_elem++] = 's';
	}
	Sv_headers.push_back("STAT");
	Sv_headers.push_back("P_QLS");
	S_fmts[N_elem++] = 'r';
	S_fmts[N_elem++] = 'r';
	S_fmts[N_elem] = '\0';
	cTableExporter	C_qls("simple.qls.res", Sv_headers, S_fmts,
		"Simple QLS result");

	/* Get correlation matrix */
	Ra_oCorr = getFullCorMat(Cp_anaCorr);

	/* Make correlation of samples included in analysis */
	cMatrix *Mp_phi = NULL;
	if (OPT_ENABLED(indep)) {
		Mp_phi = new cIdtMatrix(N_incSamp);
	} else {
		if (N_incSamp == N_sample)
			Ra_corr = Ra_oCorr;
		else
			Ra_corr = sseMsubsetRect(Ra_oCorr, N_sample, Ba_isInc, 1,
			N_incSamp);
		Mp_phi = new cSymMatrix(Ra_corr, N_incSamp, 1);
	}

	/* Get the inverse of correlation matrix */
	cMatrix& M_phiInv = Mp_phi->inv();
// 	if (Cp_IO->getTwinSize() || OPT_ENABLED(ginv))
// 		SVDinverse(Ra_corr, N_incSamp, &Ra_Rinv);
// 	else {
// 		Ra_Rinv = invSymMat(Ra_corr, N_incSamp);
// 		if (Ra_Rinv == NULL)
// 			halt("Failed to get an inverse of correlation matrix");
// 	}

	/* Deallocate correlation matrix if required */
	if (N_incSamp != N_sample)
		sseUnmat(Ra_corr, N_incSamp);

	/* Calculate V = Rinv %*% 1 %*% (1' %*% Rinv %*% 1)^-1 %*% 1' %*% Rinv
	 *  equal to rowsum's wise-product
	 *     [ r1*r1 r1*r2 r1*r3 ... r1*rN ]
	 * V = [       r2*r2    ....     .   ] / (1' %*% Rinv %*% 1)^-1 %*% 1'
	 *   = [  sym.                 rN*rN ]
	 *   
	 * And get R^-1 - V */
	sseMalloc(Ra_curX, wsReal, N_incSamp);

	/* STEP 1 : Get (1' %*% Rinv %*% 1)^-1 */
// 	wsReal R_1Rinv1 = W0;
// // 	for (i=0 ; i<N_incSamp ; i++)
// // 		for (wsUint j=0 ; j<N_incSamp ; j++)
// // 			R_1Rinv1 += Ra_Rinv[i][j];
// 	R_1Rinv1 = M_phiInv.sum();
// //	wsReal **Ra_Rinv = M_phiInv.get();
// 
// 	/* STEP 2 : Get rowsum */
// //	wsReal *Ra_r = sseMsumR(Ra_Rinv, N_incSamp); /* CHECKED */
// 	cVector	V_r		= M_phiInv.sumR();
// 	//wsReal	*Ra_r	= V_r.get();
// 
// 	/* STEP 3 : Get R^-1 - V */
// // 	for (i=0 ; i<N_incSamp ; i++) {
// // 		sseVsV(Ra_Rinv[i], Ra_r, Ra_Rinv[i], i+1, Ra_r[i]/R_1Rinv1);
// // // 		for (wsUint j=0 ; j<N_incSamp ; j++)
// // // 			Ra_Rinv[i][j] -= Ra_r[i]*Ra_r[j] / R_1Rinv1;
// // 	}
// 	V_r.rem();

	/* STEP 4 : Get divider part : yV <- (Y-mu)'(R^-1-V) */
	cVector V_availY(Ra_availY, N_incSamp, 1);
	cVector V_yV;
	if (OPT_ENABLED(indep))
		V_yV = V_availY * (cIdtMatrix &)M_phiInv;
	else
		V_yV = V_availY * (cSymMatrix &)M_phiInv;
//	wsReal *Ra_yV = sseVS(Ra_availY, N_incSamp, Ra_Rinv, N_incSamp); /* CHECKED */

	/* STEP 5 : Get denominator part : (Y-mu)'(R^-1-V)(Y-mu)
	 *  - can be optimized to yV %*% (Y-mu) */
	wsReal R_yVy = V_yV.sum(V_availY);
	delete &M_phiInv;

	/* For all SNPs */
	i = 0;
	FOREACHDO (vVariant_it, Xa_snp, it, i++) {
		cTimer	t;
		xMaf&	X_maf	= getIO()->getMAF()[i];
		wsReal	R_Tqls	= WISARD_NAN;
		wsReal	R_pVal	= WISARD_NAN;

		_doTest(i, Na_data, Mp_phi, Ra_ppp, V_yV, R_yVy, R_Tqls, R_pVal);

		t.start();
		if (OPT_ENABLED(avail)) {
			/* Copy Ba_isInc */
			char *Ba_curIsInc = NULL;
			wsAlloc(Ba_curIsInc, char, N_sample);
			memcpy(Ba_curIsInc, Ba_isInc, sizeof(char)*N_sample);
			cVector	V_curY(N_incSamp);
			wsReal	*Ra_curY	= V_curY.get();

			/* Which are missing? */
			wsUint N_curSamp = 0;
			for (wsUint j=0,K=0 ; j<N_sample ; j++) {
				if (!Ba_isInc[j]) continue;

				char N_geno = Na_data[j][i];
				if (isMissing(N_geno))
					Ba_curIsInc[j] = 1;
				else {
					Ra_curY[N_curSamp] = Ra_availY[K];
					Ra_curX[N_curSamp++] = N_geno;
				}
				K++;
			}

			/* Get phi */
			cMatrix&	M_cPhi	= Mp_phi->subset(Ba_curIsInc, 1);
			cMatrix&	M_cPinv	= M_cPhi.inv();
			
			/* STEP 1 : Get (1' %*% Rinv %*% 1)^-1 */
			wsReal		R_c1p1	= M_cPinv.sum();
			wsMat		Ra_cPinv= M_cPinv.get();

			/* STEP 2 : Get rowsum */
			cVector		V_cR	= M_cPinv.sumR();
			wsReal		*Ra_cR	= V_cR.get();

			/* STEP 3 : Get R^-1 - V */
			for (i=0 ; i<N_incSamp ; i++)
				sseVsV(Ra_cPinv[i], Ra_cR, Ra_cPinv[i], i+1, Ra_cR[i]/R_c1p1);
			V_cR.rem();

			/* STEP 4 : Get divider part : yV <- (Y-mu)'(R^-1-V) */
			cVector		V_yV;
			if (OPT_ENABLED(indep))
				V_yV = V_curY * (cIdtMatrix &)M_cPinv;
			else
				V_yV = V_curY * (cSymMatrix &)M_cPinv;

			/* STEP 5 : Get denominator part : (Y-mu)'(R^-1-V)(Y-mu)
			 *  - can be optimized to yV %*% (Y-mu) */
			wsReal		R_cYVY	= V_yV.sum(V_curY);

			/* Get the full denominator = sqrt(sig2*yVy) */
			wsReal R_denom = sqrt(W2*Ra_ppp[i]*(W1 - Ra_ppp[i])*R_cYVY);

			/* Y already BLUP enabled */
			R_Tqls = sseVV(V_yV.get(), N_curSamp, Ra_curX) / R_denom;
			R_pVal = norprobP(R_Tqls);
		} else if (!IS_ASSIGNED(genemiss) || X_maf.R_maf >= OPT_REAL(genemiss)) {
			/* Build ith data */
 			setGeno(getIO(), i, Ra_curX, Ba_isInc);

			/* Get the full denominator = sqrt(sig2*yVy) */
			wsReal R_denom = sqrt(W2*Ra_ppp[i]*(W1 - Ra_ppp[i])*R_yVy);

			/* Y already BLUP enabled */
			cVector V_curX(Ra_curX, N_incSamp, 1);
			R_Tqls = V_yV.sum(V_curX) / R_denom;
			R_pVal = norprobP(R_Tqls);
		}

		/* --remna */
		if (OPT_ENABLED(remna) && NA(R_pVal))
			continue;
		/* --pvalrange */
		if (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), R_pVal))
			continue;

		if (OPT_ENABLED(time))
			C_qls.writeVariant(&(*it), t.getReadable(), R_Tqls, R_pVal);
		else
			C_qls.writeVariant(&(*it), R_Tqls, R_pVal);
	}

	/* Deallocate */
	sseFree(Ra_availY);
	//sseFree(Ra_yV);
	sseFree(Ra_curX);
}

#endif

} // End namespace ONETOOL
