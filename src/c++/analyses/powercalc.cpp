#include "utils/dcdflib.h"
#include "utils/matrix.h"
#include "analyses/powercalc.h"
#include "analyses/emai.h"
#include "input/sim.h"
#include "utils/stat.h"

//#define VALIDATE_POWER

namespace ONETOOL {

wsReal* _get_pStarVinv(wsUint N_samp, wsUint N_t, wsRealCst *Ra_pStar,
	wsReal R_1, wsReal R_G, wsReal R_C, wsReal R_0)
{
	wsUint j, k, l;
	/**/wsReal *Ra_pStarVinv = NULL;

	sseCalloc(Ra_pStarVinv, wsReal, N_samp*N_t);
	for (wsUint T=0,i=0 ; T<N_t ; T++) {
		for (wsUint N=0 ; N<N_samp ; N++,i++) {
			j = 0;

			/* Former C part */
			for (l=0 ; l<T ; l++) {
				for (k=0 ; k<N ; k++,j++) {
					Ra_pStarVinv[i] += R_0*Ra_pStar[j];
					//printf("0");
				}
				Ra_pStarVinv[i] += R_C*Ra_pStar[j];
				j++;
				//printf("C");
				for (k=N+1 ; k<N_samp ; k++,j++) {
					Ra_pStarVinv[i] += R_0*Ra_pStar[j];
					//printf("0");
				}
			}

			/* D part */
			for (k=0 ; k<N ; k++,j++) {
				Ra_pStarVinv[i] += R_G*Ra_pStar[j];
				//printf("G");
			}
			//printf("1");
			Ra_pStarVinv[i] += R_1*Ra_pStar[j];
			j++;
			for (k=N+1 ; k<N_samp ; k++,j++) {
				Ra_pStarVinv[i] += R_G*Ra_pStar[j];
				//printf("G");
			}

			/* Latter C part */
			for (l=T+1 ; l<N_t ; l++) {
				for (k=0 ; k<N ; k++,j++) {
					Ra_pStarVinv[i] += R_0*Ra_pStar[j];
					//printf("0");
				}
				Ra_pStarVinv[i] += R_C*Ra_pStar[j];
				j++;
				//printf("C");
				for (k=N+1 ; k<N_samp ; k++,j++) {
					Ra_pStarVinv[i] += R_0*Ra_pStar[j];
					//printf("0");
				}
			}
			//printf("\n");

			if (j != (N_samp*N_t))
				halt("SYSERR: Implemntation wrong");
		} /* END OF N loop */
	} /* END OF T loop */

	return Ra_pStarVinv;
}

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cPowerCalcAnalysis::cPowerCalcAnalysis(cIO *Cp_IO,
	cAnalysis *Cp_inpAnaCorr) : cAnalysis(Cp_IO)
{
	Cp_anaCorr = Cp_inpAnaCorr;
}

cPowerCalcAnalysis::~cPowerCalcAnalysis()
{

}

inline char doMendelianTransmission(char v1, char v2)
{
	int V = v1*3 + v2;

	switch (V) {
		/* 0 0 */ case 0: return 0;
		/* 0 1, 1 0 */ case 1: case 3:
			return (char)genBinomial(1.0, 0.5);
		/* 0 2, 2 0 */ case 2: case 6: return 1;
		/* 1 1 */ case 4:
			return (char)genBinomial(2.0, 0.5);
		/* 1 2, 2 1 */ case 5: case 7:
			return (char)genBinomial(1.0, 0.5) + 1;
		/* 2 2 */ case 8: return 2;
	}

	/* Halt because it cannot be occured */
	halt("SYSERR : Errornous transmission value [%d]", V);
	return 0;
}

void genGenotype(cIO *Cp_IO, wsReal *Ra_X, wsReal R_p)
{
	wsUint	Na_s[3]	= { 0, };
	//wsUint	SS		= 0;

	wsUint*			Na_t1	= NULL;
	wsUint*			Na_t2	= NULL;
	mDataIdx		Xm_insIdx;
	vector<vInt>	Xv_prop	= _buildPropagatePath(Cp_IO, &Na_t1, &Na_t2,
		&Xm_insIdx, 0);
	wsUint			N_tSamp	= Cp_IO->sizeSample();
	wsUint			N_aSamp	= (wsUint)Cp_IO->getSampleData().size();
	vSampPtr&		Xa_samp	= Cp_IO->getSample();
	mSamp&			Xm_samp	= Cp_IO->getSampleData();

	/* minimum allele must over */
	if (W1/(wsReal)(N_aSamp<<1) > R_p)
		halt("Too low MAF requested [%g]", R_p);

	do {
		char *Na_ovals = NULL;
		wsAlloc(Na_ovals, char, N_aSamp);

		/* Generate by pedigree */
		wsUint	N_fnd	= 0;
		FOREACH (mSamp_it, Xm_samp, i)
			if (!i->second.Xp_mat && !i->second.Xp_pat) N_fnd++;
		wsUint	S		= 0;
		wsUint	N_try	= 0;
		/* Visiting init */
		// 	int*	Na_vst	= NULL;
		// 	MULTI_MALLOC(Na_vst, int, N_aSamp);
		// 	for (wsUint i=0 ; i<N_aSamp ; i++)
		// 		Na_vst[i] = -1;

		do {
			if (++N_try == 1000)
				halt("SYSERR : Too much tries in generating genotype, maybe error?");
			S = 0;
			/* 1st pass */
			FOREACH (vInt_it, Xv_prop[0], j) {
				int x = *j;
				char V = (char)genBinomial(2, R_p);
				pverbose("Initial [%d]\n", x);
				S += Na_ovals[x] = V;
			}			
			wsReal maf = S / (wsReal)(N_fnd<<1);
			if (maf > 0.5) { S = 0; continue; }

			/* impute descendants */
			wsUint N_lv = 0;
			FOREACHDO (vector<vInt>::iterator, Xv_prop, j, N_lv++) {
				if (N_lv == 0) continue;

				FOREACH (vInt_it, *j, k) {
					int	x		= *k;
					int	N_i1	= Na_t1[x];
					int	N_i2	= Na_t2[x];

// 					if (Xa_samp[x]->N_isVisited != (N_lv+1))
// 						halt("SYSERR: Sample [%s] leveled as [%d] but visited at level [%d]",
// 							Xa_samp[x]->S_IID.c_str(), Xa_samp[x]->N_isVisited,
// 							N_lv+1);

					if (OPT_ENABLED(verbose)) {
						string S_f, S_m, S_s;
						if (N_i1 >= (int)N_tSamp) {
							FOREACH (mDataIdx_it, Xm_insIdx, q) {
								if (q->second == N_i1) {
									S_f = q->first;
									break;
								}
							}
						} else S_f = Xa_samp[N_i1]->S_IID;
						if (N_i2 >= (int)N_tSamp) {
							FOREACH (mDataIdx_it, Xm_insIdx, q) {
								if (q->second == N_i2) {
									S_m = q->first;
									break;
								}
							}
						} else S_m = Xa_samp[N_i2]->S_IID;
						if (x >= (int)N_tSamp) {
							FOREACH (mDataIdx_it, Xm_insIdx, q) {
								if (q->second == x) {
									S_s = q->first;
									break;
								}
							}
						} else S_s = Xa_samp[x]->S_IID;
						pverbose("[%d,%d] -> %d [%s,%s] -> %s\n", N_i1, N_i2, x,
							S_f.c_str(), S_m.c_str(), S_s.c_str());
					}
					char V = doMendelianTransmission(Na_ovals[N_i1], Na_ovals[N_i2]);
					S += Na_ovals[x] = V;
				}
			}
			// 		X_snv.mac = S;
			// 		X_snv.allmaf = S / (wsReal)(N_samp<<1);
		} while (S == 0); /* Avoid monomorphic SNP */
		Xm_insIdx.clear();

		/* Copy it to X and check integrity */
		Na_s[0] = Na_s[1] = Na_s[2] = 0;
		for (wsUint i=0 ; i<N_tSamp ; i++) {
			int v =  Na_ovals[i];
			Ra_X[i] = (wsReal)v;
			Na_s[v]++;
		}
	} while (((Na_s[0]!=0)+(Na_s[1]!=0)+(Na_s[2]!=0)) <= 1);
}

void cPowerCalcAnalysis::_getSinglePheno(wsUint N_perm)
{
	ASSERT_OPTION(prevalence);
	ASSERT_OPTION(heri);
//	ASSERT_OPTION(beta);

	//p<-0.2
	wsUint	N_prev	= 0;
	wsReal*	Ra_p	= loadRealValues(OPT_STRING(prevalence), &N_prev);
	wsReal	R_p		= Ra_p[0];
	wsReal	R_p1p	= R_p * (W1 - R_p);
	wsReal	R_h2	= atof(OPT_STRING(heri));
	wsReal	R_beta	= OPT_REAL(beta);
	wsUint	N_samp	= getIO()->sizeSample();
//	wsReal	R_beta	= sqrt(0.005 / (REAL_CONST(1.99)*R_p1p));

	LOG("Perform power analysis of sample size [%d] with [%d] permutations\n",
		N_samp, N_perm);
	if (N_prev == 2)
		LOGwarn("Sex-wise prevalence is given, only the first prevalence "
			"will used for sample size calculation");

// 	p0<-dbinom(0,2,p)
// 	wsReal R_p0 = dbinom(0, 2, R_p);
// // 	p1<-dbinom(1,2,p)
// 	wsReal R_p1 = dbinom(1, 2, R_p);
// 	p2<-dbinom(2,2,p)
	//wsReal R_p2 = dbinom(2, 2, R_p);

	wsReal R_2bbp1p = W2 * SQR(R_beta) * R_p1p;
	// sig.g<-h2-2*beta*beta*p*(1-p)
//	wsReal R_sigg = R_h2*(W1+R_2bbp1p) - R_2bbp1p;
	wsReal R_sigg = R_h2 - R_2bbp1p;
	// sig.e<-1-sig.g-2*beta*beta*p*(1-p)
//	wsReal R_sige = W1 - R_sigg;// - R_2bbp1p;
	wsReal R_sige = W1 - R_sigg - R_2bbp1p;

	if (R_sigg < W0)
		halt("Negative sig^2(g) induced, too low h^2[%g] or too high beta[%g]",
			R_h2, R_beta);
		
	LOG("sig^2(g) : %g\n", R_sigg);
	LOG("sig^2(e) : %g\n", R_sige);

	/* Correlation matrix */
	wsSym		Ra_corOrig	= getFullCorMat(Cp_anaCorr);
	wsSym		Ra_cor		= sseSymMatP(N_samp, Ra_corOrig);
	cSymMatrix	M_cor(Ra_cor, N_samp);
#ifdef VALIDATE_POWER
	M_cor.file("cor.g");
#endif

	/* Correlation-related */
	//var.v<-sig.g*cor.g+sig.e*i.n */
	M_cor *= R_sigg;
	for (wsUint i=0 ; i<N_samp ; i++)
		Ra_cor[i][i] += R_sige;
#ifdef VALIDATE_POWER
	M_cor.file("var.v");
#endif

	// sol.v<-solve(var.v)
	cSymMatrix&	M_Vinv		= M_cor.inv();
#ifdef VALIDATE_POWER
	M_Vinv.file("sol.v");
#endif
	wsReal		R_sumVinv	= W0;
///**/wsReal **Ra_Vinv = invSymMat(Ra_cor, N_samp);
//	sseUnmat(Ra_cor, N_samp);
//	wsReal R_sumVinv = sseMsum(Ra_Vinv, N_samp);

	/* sol.v%*%colv.n%*%solve(sumVinv)%*%t(colv.n)
	 * n*n matrix but elements in same row are all same */
	cVector		V_sumVrow	= M_Vinv.sumR(&R_sumVinv);
	wsReal		R_invSVi	= W1 / R_sumVinv;
	cVector		V_svv		= V_sumVrow * R_invSVi;

	/* Genotype buffer */
	wsReal *Ra_X = sseVector(N_samp);

	/* Generate storing buffer */
	wsReal *Ra_pow1	= sseVector(N_perm);
	wsReal *Ra_pow2	= sseVector(N_perm);
	wsReal *Ra_nc	= sseVector(N_perm);
	for (wsUint i=0 ; i<N_perm ; i++) {
		/* Generate genotype */
		genGenotype(getIO(), Ra_X, R_p);
		cVector		V_x(Ra_X, N_samp, 1);
#ifdef VALIDATE_POWER
		V_x.file("x");
#endif

		wsReal		R_xV1	= V_x.sum(V_sumVrow);
		wsReal		R_xSub	= V_x.sum(V_svv);
	
		/* pstar<-t(x)-t(x)%*%sol.v%*%colv.n%*%solve(sumVinv)%*%t(colv.n)
		 * pstar<-t(x)-xSub %*% t(colv.n) */
		cVector		V_pStar		= V_x - R_xSub;
#ifdef VALIDATE_POWER
		V_pStar.file("pstar");
#endif
		/* xi<-pstar%*%sol.v%*%(x*beta-colv.n%*%solve(sumVinv)%*%t(t(x)%*%sol.v%*%colv.n)*beta)
		 * xi<-pstar%*%sol.v%*%(x*beta-colv.n%*%xV1*beta/sumVinv) */
		wsReal		R_xbSub		= R_xV1*R_invSVi * R_beta;
		cVector		V_xBeta		= V_x*R_beta;
		V_xBeta -= R_xbSub;
#ifdef VALIDATE_POWER
		V_xBeta.file("xisub");
#endif
		wsReal	R_xi		= V_pStar.qf(M_Vinv, V_xBeta);
		wsReal	R_sigma1	= V_pStar.qf(M_Vinv);
	
		// non_cen_par<-(t(xi)/sigma1)%*%xi
		wsReal	R_nc		= (R_xi/R_sigma1) * R_xi;
		Ra_nc[i]			= R_nc;
		double	R_dnc		= (double)R_nc;

		// power<-pchisq(qchisq(0.999999,df=1),df=1,ncp=non_cen_par,lower.tail=F)
		//wsReal R_pow = chiprobP(23.92812697687946865699, 1.0, &R_nc);
		Ra_pow1[i]			= PVchisq(3.841458820694124032258, 1.0, &R_dnc);
		Ra_pow2[i]			= PVchisq(28.373987362798125389, 1.0, &R_dnc);
	}

	LOGoutput("Power estimation on single phenotype is exported to [%s.power.single.res]\n",
		OPT_STRING(out));
/**/cExporter* Cp_exp = cExporter::summon("power.single.res");
	Cp_exp->fmt("sample_size	%d\n", N_samp);
	Cp_exp->fmt("pheno_size	1\n");
//	Cp_exp->fmt("noncen_param	%g\n", R_nc);
	Cp_exp->fmt("sig2g	%g\n", R_sigg);
	Cp_exp->fmt("sig2e	%g\n", R_sige);
	Cp_exp->fmt("power	%g\n", sseVsum(Ra_pow1, N_perm)/(wsReal)N_perm);
	Cp_exp->fmt("power_gwas	%g\n", sseVsum(Ra_pow2, N_perm)/(wsReal)N_perm);

	exportVector("power.single.nc.res", Ra_nc, N_perm);
	sseFree(Ra_pow1);
	sseFree(Ra_pow2);
	sseFree(Ra_nc);

	delete Cp_exp;
	sseFree(Ra_p);
}

void cPowerCalcAnalysis::_getVinv(wsUint N_samp, wsUint N_t, wsReal R_l,
	wsReal R_sigg, wsReal R_sigc, wsReal *Rp_1, wsReal *Rp_G, wsReal *Rp_C,
	wsReal *Rp_0)
{
	wsReal R_vg	= R_l * R_sigg;
	wsReal R_vc	= R_sigc;

	// g2 <- vg*vg
	wsReal R_g2 = SQR(R_vg);

	if (N_t > 1) {
		// c2 <- vc*vc
		wsReal R_c2 = SQR(R_vc);
		// n1g2 <- n1*g2
		wsReal R_n1g2 = (wsReal)(N_samp-1)*R_g2;
		// n2vg <- n2*vg
		wsReal R_n2vg = (wsReal)(N_samp-2)*R_vg;
		// t2vc <- t2*vc
		wsReal R_t2vc = (wsReal)(N_t - 2)*R_vc;
		// t1c2 <- t1*c2
		wsReal R_t1c2 = (wsReal)(N_t - 1)*R_c2;
		// t1vc <- t1*vc
		wsReal R_t1vc = (wsReal)(N_t - 1)*R_vc;
		// R <- (t2vc+1)*(n1g2-n2vg-1)+(n2vg*vc+vc)*t1vc
		wsReal R_vR = (R_t2vc+W1)*(R_n1g2-R_n2vg-W1)
			+(R_n2vg*R_vc + R_vc)*R_t1vc;
		// Q <- (t1c2-n1g2 -t2vc-1)*t1vc
		wsReal R_vQ = (R_t1c2 - R_n1g2 - R_t2vc - W1) * R_t1vc;
		// G <- n1g2 - n2vg - 1 - t1c2
		wsReal R_vG = R_n1g2 - R_n2vg - W1 - R_t1c2;
		// H <- (2 + n2vg + t2vc)*t1vc
		wsReal R_vH = (W2 + R_n2vg + R_t2vc) * R_t1vc;
		// VgDet <- R*H + Q*G
		// Vg <- (vg*(t2vc+1)*H+Q*vg) / VgDet
		*Rp_G = (R_vg * (R_t2vc + W1) * R_vH + R_vQ*R_vg)
			/ (R_vR*R_vH + R_vQ*R_vG);
		// V0Det <- Q*G + R*H
		// V0 <- (G*vg*(t2vc + 1) - vg*R) / V0Det
		*Rp_0 = (R_vG*R_vg*(R_t2vc+W1) - R_vg*R_vR) /
			(R_vQ*R_vG + R_vR*R_vH);
		// Vc = (-vc*Vg - (1+n2vg + t2vc)*V0) / vg
		*Rp_C = (-R_vc**Rp_G - (W1+R_n2vg + R_t2vc)**Rp_0) / R_vg;
		// V1 <- 1 - n1*vg*Vg - t1vc*Vc
		*Rp_1 = W1 - (wsReal)(N_samp-1)*R_vg**Rp_G - R_t1vc**Rp_C;
	} else {
		*Rp_G = R_vg / ((wsReal)(N_samp-1)*R_g2 - (wsReal)(N_samp-2)*R_vg -
			W1);
		*Rp_1 = W1 - (wsReal)(N_samp-1)*R_vg**Rp_G;
		*Rp_C = W0;
		*Rp_0 = W0;
	}
}

void cPowerCalcAnalysis::_getMultiPheno(wsUint N_t, wsUint N_s,
	wsUint N_e/*=0xffffffff*/)
{
	ASSERT_OPTION(prevalence);
	ASSERT_OPTION(rhopheno);
	ASSERT_OPTION(rho);
	ASSERT_OPTION(heri);

	/* Cannot concurrently define --noshuffle and --nperm */
	if (OPT_ENABLED(noshuffle) && IS_ASSIGNED(nperm))
		halt("--noshuffle and --nperm cannot concurrently assigned");

/**/cExporter* Cp_exp = cExporter::summon("power.multi.res");
	wsReal R_rho	= OPT_REAL(rho);
	wsUint N_prev	= 0;
	/*	rho=0.7###*/
	wsReal *Ra_p	= loadRealValues(OPT_STRING(prevalence), &N_prev);
	wsReal R_p		= Ra_p[0];
	wsReal R_h2		= atof(OPT_STRING(heri));
	wsReal R_1p		= W1-R_p;
	if (N_t < 2)
		halt("Cannot perform this analysis : #pheno should be > 1");

	if (N_prev == 2)
		LOGwarn("Sex-wise prevalence is given, only the first prevalence "
			"will used for sample size calculation\n");

	// 	p0<-dbinom(0,2,p)
	wsReal R_p0 = dbinom(0, 2, R_p);
	// 	p1<-dbinom(1,2,p)
	wsReal R_p1 = dbinom(1, 2, R_p);
	// 	p2<-dbinom(2,2,p)
	//wsReal R_p2 = dbinom(2, 2, R_p);

	// beta<-sqrt(0.005/(2*(1-0.005)*p*(1-p)))
	wsReal R_beta = sqrt(REAL_CONST(0.005) /
		(W2 * (wsReal)(1.0-0.005) * R_p * R_1p));
	/*
	sig<-1
	sig.g<-h2*(1+2*beta*beta*p*(1-p))-2*beta*beta*p*(1-p)*/
	wsReal R_2b2p1p = W2*SQR(R_beta)*R_p*R_1p;
	wsReal R_sigg = R_h2*(W1+R_2b2p1p)-R_2b2p1p;
	// sig.c<-rho-sig.g
	wsReal R_sigc = R_rho - R_sigg;
	// sig.e<-1-sig.g-sig.c
	//wsReal R_sige = W1 - R_sigg - R_sigc;

	// n1<-round(n*p0)
	if (N_e == 0xffffffff) N_e = N_s;

	/* Set the number of permutation */
	wsUint N_perm = IS_ASSIGNED(nperm) ? OPT_NUMBER(nperm) : 1;

	/* Export header if have range */
	if (N_s != N_e) {
		N_perm==1 ?
			Cp_exp->put("SAMPSIZE	PHENOSIZE	NONCEN_PARAM	POWER\n") :
			Cp_exp->put("SAMPSIZE	PHENOSIZE	NUM_PERM	POWER\n");
	}

	// vc <- sig.c
	// vg <- l*sig.g
	for (wsUint N_samp=N_s ; N_samp<=N_e ; N_samp++) {
		notice("Calculating power for sample size=%d\r", N_samp);

		double R_nc = W0;
		wsReal *Ra_pow = NULL;
		sseMalloc(Ra_pow, wsReal, N_perm);

		for (wsUint N_loop=0 ; N_loop<N_perm ; N_loop++) {
			wsUint N_n1 = (wsUint)(R_p0 * (wsReal)N_samp);
			// n2<-round(n*p1)
			wsUint N_n2 = (wsUint)(R_p1 * (wsReal)N_samp);
			// n3<-n-n1-n2
			wsUint N_n3 = N_samp - N_n1 - N_n2;

			/* x<-matrix(rep(sample(xx),t),ncol=1) */
		/**/wsReal *Ra_X = NULL;
			sseMalloc(Ra_X, wsReal, N_samp*N_t);
			if (OPT_ENABLED(noshuffle)) {
				/* Do not sampling samples */
				wsUint j=0;
				for (wsUint k=0 ; k<N_t ; k++) {
					for (wsUint i=0 ; i<N_n1 ; i++,j++) Ra_X[j] = 0;
					for (wsUint i=0 ; i<N_n2 ; i++,j++) Ra_X[j] = 1;
					for (wsUint i=0 ; i<N_n3 ; i++,j++) Ra_X[j] = 2;
				}
			} else {
				/*
				x1<-matrix(rep(0,n1),ncol=1)
				x2<-matrix(rep(1,n2),ncol=1)
				x3<-matrix(rep(2,n3),ncol=1)
				xx<-rbind(x1,x2,x3)*/
		/**/	char *Na_vals = NULL;
				wsAlloc(Na_vals, char, N_samp*N_t);
				memset(Na_vals, 0x00, N_n1*N_t);
				memset(Na_vals+N_n1*N_t, 0x01, N_n2*N_t);
				memset(Na_vals+(N_n1+N_n2)*N_t, 0x02, N_n3*N_t);
				for (wsUint i=0,j=N_samp*N_t ; i<(N_samp*N_t) ; i++,j--) {
					wsUint N_selIdx = rand()%j;
	
					Ra_X[i] = Na_vals[N_selIdx];
					/* Do erase */
					memmove(Na_vals+N_selIdx, Na_vals+N_selIdx+1,
						sizeof(char)*(j-N_selIdx-1));
				}
				DEALLOC(Na_vals);
			}

			wsReal R_1, R_G;
			wsReal R_C, R_0;

			/* Get V^-1 */
			_getVinv(N_samp, N_t, OPT_REAL(rhopheno), R_sigg, R_sigc,
				&R_1, &R_G, &R_C, &R_0);

			// sumVinv <- n*t*V1 + n*t*n1*Vg + n*t*t1*Vc + n*t*t1*n1*V0
			wsReal R_sumVinv = (wsReal)(N_samp*N_t)*R_1 +
				(wsReal)(N_samp*N_t*(N_samp-1))*R_G +
				(wsReal)(N_samp*N_t*(N_t-1))*R_C +
				(wsReal)(N_samp*N_t*(N_t-1)*(N_samp-1))*R_0;

			/*
			l=0.1####
			cor.g<-matrix(rep(l,n*n),ncol=n)
			diag(cor.g)<-1


			i.t<-diag(t)
			i.n<-diag(n)
			i.nt<-diag(n*t)
			colv.t<-matrix(rep(1,t),ncol=1)
			colv.nt<-matrix(rep(1,n*t),ncol=1)
			// Randomly assigned column vector, having n1 0s, n2 1s, n3 2s
			x<-matrix(rep(sample(xx),t),ncol=1) */


			// var.v<-sig.g*kronecker(i.t,cor.g)+
			//		sig.c*kronecker(colv.t%*%t(colv.t),i.n)+
			//		sig.e*kronecker(i.t,i.n) == i_tn
			// 
		// 	sol.v<-solve(var.v)
			/* Now we have sol.v as the form of [ [R_1   R_C  ] [R_G   R_0  ]
			 *                                  [ [    ...    ] [    ...    ] ...
			 *                                  [ [  R_C   R_1] [  R_0   R_G]
			 *                                  [               ...
			 */
			/* t(colv.nt)%*%sol.v%*%colv.nt have
			 *  - nt           R_1
			 *  - nt(n-1)      R_G
			 *  - nt(t-1)      R_C
			 *  - nt(t-1)(n-1) R_0 */

			// 	Vinv1 <- sol.v%*%colv.nt # Repeatition of r1 * (n-1)*rg + (k-1)*rc + (n-1)*(k-1)r0
			wsReal R_rowSumVinv = R_1 + (wsReal)(N_samp-1)*R_G + (wsReal)(N_t-1)*R_C + 
				(wsReal)(N_samp-1) * (wsReal)(N_t-1) * R_0;
			// 	Vq <- Vinv1%*%solve(sumVinv)%*%t(colv.nt)
			// 	    is repeatition of 1/sumVinv * Vinv[1]
			wsReal R_Vq = W1/R_sumVinv * R_rowSumVinv;

			wsReal R_xSum = sseVsum(Ra_X, N_samp*N_t);

			// 	pstar<-t(x)-t(x)%*%Vq
			//   Get sum(x)*Vq
			wsReal R_xVq = R_xSum * R_Vq;
		/**/wsReal *Ra_pStar = NULL;
			sseMalloc(Ra_pStar, wsReal, N_samp*N_t);
			sseVsC(Ra_X, R_xVq, Ra_pStar, N_samp*N_t);

			// xVb <- t(t(x)%*%Vinv1)*beta)
			wsReal R_xVb = R_xSum*R_rowSumVinv*R_beta;

			// xBeta <- x*beta-colv.nt%*%solve(sumVinv)%*%xVb
		/**/wsReal *Ra_xBeta = NULL;
			sseMalloc(Ra_xBeta, wsReal, N_samp*N_t);
			sseVpC(Ra_X, R_beta, Ra_xBeta, N_samp*N_t);
			sseFree(Ra_X);
			sseVsC(Ra_xBeta, R_xVb/R_sumVinv, Ra_xBeta, N_samp*N_t);

			// xi<-pstar%*%sol.v%*%xBeta
			//                [ D C .. C ] p1*1 + p2...pn*g + p(n+1)*c + ... + p((k-1)n+1)*c
			// [p1 p2 ... pn] [ C D .. C ] p1*g + p2*1 + p3...pn*g + p(n+2)*c + ... + p((k-1)n+2)*c
			//                [ C C .. D ]          ...
			//                             p1*c + p(n+1) + p(n+2)...p(2n)*g + p(2n+1)*c + ... + p((k-1)n+2)*c
/**/		wsReal *Ra_pStarVinv = _get_pStarVinv(N_samp, N_t, Ra_pStar,
				R_1, R_G, R_C, R_0);

			wsReal R_xi = sseVV(Ra_pStarVinv, N_samp*N_t, Ra_xBeta);
			sseFree(Ra_xBeta);
			// 	sigma1<-pstar%*%sol.v%*%t(pstar)
			wsReal R_sigma1 = sseVV(Ra_pStarVinv, N_samp*N_t, Ra_pStar);
			sseFree(Ra_pStarVinv);
			sseFree(Ra_pStar);

			// 	non_cen_par<-(t(xi)/sigma1)%*%xi
			R_nc = (double)((R_xi/R_sigma1)*R_xi);
		// 	power<-pchisq(qchisq(0.95,df=1),df=1,ncp=non_cen_par,lower.tail=F)
			Ra_pow[N_loop] = PVchisq(3.841458820694124032258, 1.0, &R_nc);
		}

		/* Get the mean of power */
		wsReal R_meanPow = sseVsum(Ra_pow, N_perm) / (wsReal)N_perm;

		if (N_s == N_e) {
			/* For the output of sample range is 1 */
			Cp_exp->fmt("sample_size	%d\n", N_samp);
			Cp_exp->fmt("pheno_size	%d\n", N_t);
			N_perm == 1 ?
				Cp_exp->fmt("noncen_param	%g\n", R_nc) :
				Cp_exp->fmt("num_perm	%d\n", N_perm);
			Cp_exp->fmt("Power	%g\n", R_meanPow);
		} else {
			N_perm == 1 ?
				Cp_exp->fmt("%d	%d	%g	%g\n", N_samp, N_t, R_nc, R_meanPow) :
				Cp_exp->fmt("%d	%d	%d	%g\n", N_samp, N_t, N_perm, R_meanPow);
		}
	}
	notice("Calculating powers for %d sample sizes performed\n", N_e-N_s+1);

	sseFree(Ra_p);
	delete Cp_exp;
// 	res1<-power
// 	res<-cbind(h2,t,rho,res1,n,l)
// 	return(res)
// 	*/
}

void cPowerCalcAnalysis::run()
{
	wsUint N_t = 0, N_samp = 0xffffffff;

	if (IS_ASSIGNED(npheno)) {
		N_t = OPT_NUMBER(npheno);
	} else {
		if (Cp_IO == NULL)
			halt("Power calculation requires appropriate option");
		N_t = Cp_IO->sizePheno();
	}

	if (!IS_ASSIGNED(nsamp)) {
		if (Cp_IO == NULL)
			halt("Power calculation requires appropriate option");
		N_samp = Cp_IO->sizeSample();
	}

	if (N_t > 1) {
		if (N_samp != 0xffffffff)
			_getMultiPheno(N_t, N_samp);
		else {
			xOptRange &X_rng = OPT_RANGE(nsamp);

			if (X_rng.R_s == numeric_limits<wsReal>::infinity() ||
				X_rng.R_e == numeric_limits<wsReal>::infinity())
				halt("Power calculation requires exact value or range");
			if (X_rng.R_s != X_rng.R_s || X_rng.R_e != X_rng.R_e)
				halt("SYSERR: nsamp range is invalid");

			wsUint N_s = X_rng.R_sEQ ? (wsUint)X_rng.R_s : ((wsUint)X_rng.R_s)+1;
			wsUint N_e = X_rng.R_eEQ ? (wsUint)X_rng.R_e : ((wsUint)X_rng.R_e)-1;

			if (N_s > N_e) halt("Invalid range given, range check is required");

			if (N_s == N_e)
				_getMultiPheno(N_t, N_s);
			else
				_getMultiPheno(N_t, N_s, N_e);
		}
	} else {
		wsUint N_perm = IS_ASSIGNED(nperm) ? OPT_NUMBER(nperm) : 1;
		if (Cp_anaCorr == NULL)
			halt("Power calculation for single phenotype requires dataset");
		else if (IS_ASSIGNED(npheno))// || IS_ASSIGNED(nsamp))
			halt("--npheno"/* or --nsamp*/" cannot be used when correlation matrix "
				"is given (via input or --cor)");
		_getSinglePheno(N_perm);
	}
}

cPowerCalcAnalysisV2::cPowerCalcAnalysisV2(cIO *Cp_IO,
	cAnalysis *Cp_inpAnaCorr) : cAnalysis(Cp_IO)
{
	Cp_anaCorr = Cp_inpAnaCorr;
}

cPowerCalcAnalysisV2::~cPowerCalcAnalysisV2()
{

}

void cPowerCalcAnalysisV2::_getSinglePheno(wsUint N_samp)
{
	ASSERT_OPTION(prevalence);
	ASSERT_OPTION(heri);
//	ASSERT_OPTION(beta);

	//p<-0.2
	wsUint	N_prev	= 0;
	wsReal*	Ra_p	= loadRealValues(OPT_STRING(prevalence), &N_prev);
	wsReal	R_p		= Ra_p[0];
	wsReal	R_p1p	= R_p * (W1 - R_p);
	wsReal	R_h2	= atof(OPT_STRING(heri));
//	wsReal	R_beta	= OPT_REAL(beta);
	wsReal	R_beta	= sqrt(0.005 / (REAL_CONST(1.99)*R_p1p));

	if (N_prev == 2)
		LOGwarn("Sex-wise prevalence is given, only the first prevalence "
			"will used for sample size calculation");

// 	p0<-dbinom(0,2,p)
// 	wsReal R_p0 = dbinom(0, 2, R_p);
// // 	p1<-dbinom(1,2,p)
// 	wsReal R_p1 = dbinom(1, 2, R_p);
// 	p2<-dbinom(2,2,p)
	//wsReal R_p2 = dbinom(2, 2, R_p);

	wsReal R_2bbp1p = W2 * SQR(R_beta) * R_p1p;
	// sig.g<-h2-2*beta*beta*p*(1-p)
	wsReal R_sigg = R_h2 - R_2bbp1p;
	// sig.e<-1-sig.g-2*beta*beta*p*(1-p)
	wsReal R_sige = W1 - R_sigg - R_2bbp1p;

	LOG("sig^2(g) : %g\n", R_sigg);
	LOG("sig^2(e) : %g\n", R_sige);
		
	/*
	x1<-matrix(rep(0,n1),ncol=1)
	x2<-matrix(rep(1,n2),ncol=1)
	x3<-matrix(rep(2,n3),ncol=1)
	xx<-rbind(x1,x2,x3) */
	// n<-i+780####

	/*
	x1<-matrix(rep(0,n1),ncol=1)
	x2<-matrix(rep(1,n2),ncol=1)
	x3<-matrix(rep(2,n3),ncol=1)
	xx<-rbind(x1,x2,x3)*/

	// co<-read.table('kare_fam.cor')
	wsReal **Ra_corOrig = getFullCorMat(Cp_anaCorr);
/**/wsReal **Ra_cor = NULL;
	/* Permute sample */
	// l<-sample(1:nrow(co),n,replace=T)
	// cor.g<-as.matrix(co[l,l])
	Ra_cor = OPT_ENABLED(noshuffle) ?
		/* Do not shuffling correlation matrix */
		sseMatrixP(N_samp, N_samp, Ra_corOrig) :
		/* Do shuffling correlation matrix */
		sseShuffleMatrix(Ra_corOrig, Cp_IO->sizeSample(), N_samp);
	cSymMatrix	M_cor(Ra_cor, N_samp);

	wsSym Ra_cor2 = sseSymMatP(N_samp, Ra_cor);
	cSymMatrix	M_cor2(Ra_cor2, N_samp);
	M_cor2 *= W2 * R_p1p;
	M_cor2 += REAL_CONST(4.0) * R_p * R_p;

	//M_cor.file("xx");
	//LOG("Matrixsum %g\n", sseMsum(Ra_cor, N_samp));

	// I_n
	/*i.n<-diag(n)
	// 1_n
	colv.n<-matrix(rep(1,n),ncol=1)
	// x<-matrix(sample(xx),ncol=1)
	
	//var.v<-sig.g*cor.g+sig.e*i.n */
	M_cor *= R_sigg;
//	sseMpC(Ra_cor, R_sigg, Ra_cor, N_samp);
	for (wsUint i=0 ; i<N_samp ; i++)
		Ra_cor[i][i] += R_sige;

	// sol.v<-solve(var.v)
	cSymMatrix&	M_Vinv		= M_cor.inv();
	wsReal		R_trXVX		= sseTrSS(M_Vinv.get(), N_samp, Ra_cor2);
	double		R_nc		= R_beta * sqrt((double)N_samp) / sqrt(R_trXVX);

	// power<-pchisq(qchisq(0.999999,df=1),df=1,ncp=non_cen_par,lower.tail=F)
	//wsReal R_pow = chiprobP(23.92812697687946865699, 1.0, &R_nc);
	wsReal R_pow = PVchisq(3.841458820694124032258, 1.0, &R_nc);

	LOGoutput("Power estimation on single phenotype is exported to [%s.power.single.res]\n",
		OPT_STRING(out));
/**/cExporter* Cp_exp = cExporter::summon("power2.single.res");
	Cp_exp->fmt("sample_size	%d\n", N_samp);
	Cp_exp->fmt("pheno_size	1\n");
	Cp_exp->fmt("noncen_param	%g\n", R_nc);
	Cp_exp->fmt("sig2g	%g\n", R_sigg);
	Cp_exp->fmt("sig2e	%g\n", R_sige);
	Cp_exp->fmt("Power	%g\n", R_pow);
	delete Cp_exp;
	sseFree(Ra_p);
}

void cPowerCalcAnalysisV2::_getMultiPheno(wsUint N_t, wsUint N_s,
	wsUint N_e/*=0xffffffff*/)
{
	halt("Not implemented");
}

void cPowerCalcAnalysisV2::run()
{
	wsUint N_t = 0, N_samp = 0xffffffff;

	if (IS_ASSIGNED(npheno)) {
		N_t = OPT_NUMBER(npheno);
	} else {
		if (Cp_IO == NULL)
			halt("Power calculation requires appropriate option");
		N_t = Cp_IO->sizePheno();
	}

	if (!IS_ASSIGNED(nsamp)) {
		if (Cp_IO == NULL)
			halt("Power calculation requires appropriate option");
		N_samp = Cp_IO->sizeSample();
	}

	if (N_t > 1) {
		if (N_samp != 0xffffffff)
			_getMultiPheno(N_t, N_samp);
		else {
			xOptRange &X_rng = OPT_RANGE(nsamp);

			if (X_rng.R_s == numeric_limits<wsReal>::infinity() ||
				X_rng.R_e == numeric_limits<wsReal>::infinity())
				halt("Power calculation requires exact value or range");
			if (X_rng.R_s != X_rng.R_s || X_rng.R_e != X_rng.R_e)
				halt("SYSERR: nsamp range is invalid");

			wsUint N_s = X_rng.R_sEQ ? (wsUint)X_rng.R_s : ((wsUint)X_rng.R_s)+1;
			wsUint N_e = X_rng.R_eEQ ? (wsUint)X_rng.R_e : ((wsUint)X_rng.R_e)-1;

			if (N_s > N_e) halt("Invalid range given, range check is required");

			if (N_s == N_e)
				_getMultiPheno(N_t, N_s);
			else
				_getMultiPheno(N_t, N_s, N_e);
		}
	} else {
		if (N_samp == 0xffffffff)
			N_samp = (wsUint)(OPT_RANGE(nsamp).R_s);
		if (Cp_anaCorr == NULL)
			halt("Power calculation for single phenotype requires dataset");
		else if (IS_ASSIGNED(npheno))// || IS_ASSIGNED(nsamp))
			halt("--npheno"/* or --nsamp*/" cannot be used when correlation matrix "
				"is given (via input or --cor)");
		_getSinglePheno(N_samp);
	}
}

#endif

} // End namespace ONETOOL
