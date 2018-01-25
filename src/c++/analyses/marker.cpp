#include "analyses/marker.h"
#include "utils/stat.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cVariantAnalysis::cVariantAnalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS) : cAnalysis(Cp_inpIO)
{
	Cp_anaFS = Cp_inpAnaFS;
}

cVariantAnalysis::~cVariantAnalysis()
{

}

typedef struct _xVariantSumm {
	wsUint	N_snv, N_ins, N_del, N_sub, N_mono, N_single, N_double;
	int		N_minPos, N_maxPos;
	wsUint	N_lenIns, N_minIns, N_maxIns;
	wsUint	N_lenDel, N_minDel, N_maxDel;
	wsUint	N_ti, N_tv;
} xVariantSumm;

/* Fisher's exact test table */
typedef struct _xFisher {
	wsUint Na_freq[2][2];
} xFisher;

/* Cochran-Armitage trend test table */
typedef struct _xTrend {
	wsUint N_freq[2][3];
} xTrend;

void cVariantAnalysis::run()
{
	xVariantSumm	Xa_msum[MAX_NCHR+1] = { { 0, }, };
	wsUint		N_mkr	= getIO()->sizeVariant();
	wsUint		N_samp	= getIO()->sizeSample();
	char**		Na_data	= getIO()->getGenotype();
	vSampPtr&	Xv_samp	= getIO()->getSample();
	vVariant&	Xv_mkr	= getIO()->getVariant();
	xMaf*		Xp_maf	= getIO()->getMAF();
	wsUintCst		N_pheno	= getIO()->sizePheno();
	vPheno&		Xv_pi	= getIO()->getPhenoInfo();

	/* Fisher's exact test */
	xFisher*	Xa_fisher	= NULL;
	cExporter*	Cp_fe		= NULL;
	wsMat		Ra_pheno	= getIO()->getPhenos();
	if (OPT_ENABLED(fisher)) {
		LOG("Fisher's Exact Test for case-control association testing "
			"will be performed\n");
		for (wsUint i=0 ; i<N_pheno ; i++)
			if (getIO()->isContinuous(i))
				halt("[%d]th phenotype [%s] is not binary, Fisher's exact "
					"test cannot be performed!", i+1, Xv_pi[i].S_name.c_str());
		Cp_fe = cExporter::summon("fisher.res");
		LOGoutput("Result of Fisher's exact test is exported to [%s.fisher.res]\n", OPT_STRING(out));
		headerVariant(Cp_fe);
		Cp_fe->put("	PHENO	P_FISHER\n");
		wsAlloc(Xa_fisher, xFisher, N_pheno);
	}

	/* Cochran-Armitage trend test */
	xTrend*		Xa_trend	= NULL;
	cExporter*	Cp_ca		= NULL;
	if (OPT_ENABLED(trend)) {
		LOG("Cochran-Armitage trend test and allelic based test for "
			"case-control association testing will be performed\n");
		for (wsUint i=0 ; i<N_pheno ; i++)
			if (getIO()->isContinuous(i))
				halt("[%d]th phenotype [%s] is not binary, Cochran-Armitage "
					"trend test cannot be performed!", i+1, Xv_pi[i].S_name.c_str());
		Cp_ca = cExporter::summon("trend.res");
		LOGoutput("Result of Cochran-Armitage trend test and allelic based "
			"test are exported to [%s.trend.res]\n", OPT_STRING(out));
		headerVariant(Cp_ca);
		Cp_ca->put("	PHENO	P_TREND	P_ABT	P_GENO2DF\n");
		wsAlloc(Xa_trend, xTrend, N_pheno);
	}

	/* Family-specific variant discovery */
	cExporter *Cp_fs = NULL;
	if (OPT_ENABLED(famuniq)) {
		LOGoutput("List of family-specific variants is exported to [%s.variant.famuniq.res]\n",
			OPT_STRING(out));
		Cp_fs = cExporter::summon("variant.famuniq.res");
		headerVariant(Cp_fs);
		Cp_fs->put("	FID	COUNT\n");
	}

	/* Population-specific variant discovery */
	cExporter *Cp_ps = NULL;
	if (IS_ASSIGNED(popuniq)) {
		LOGoutput("List of population-specific variants is exported to [%s.variant.popuniq.res]\n",
			OPT_STRING(out));
		Cp_ps = cExporter::summon("variant.popuniq.res");
		headerVariant(Cp_ps);
		Cp_ps->put("	POP_ID	COUNT\n");
	}

	cExporter *Cp_st = NULL;
	if (OPT_ENABLED(singleton)) {
		LOGoutput("List of singleton variants is exported to [%s.variant.singleton.res]\n",
			OPT_STRING(out));
		Cp_st = cExporter::summon("variant.singleton.res");
		headerVariant(Cp_st);
	}

	cExporter *Cp_mt = NULL;
	if (OPT_ENABLED(monotone)) {
		LOGoutput("List of monotone variants is exported to [%s.variant.monotone.res]\n",
			OPT_STRING(out));
		Cp_mt = cExporter::summon("variant.monotone.res");
		headerVariant(Cp_mt);
	}

	cExporter *Cp_dt = NULL;
	if (OPT_ENABLED(doubleton)) {
		LOGoutput("List of doubleton variants is exported to [%s.variant.doubleton.res]\n",
			OPT_STRING(out));
		Cp_dt = cExporter::summon("variant.doubleton.res");
		headerVariant(Cp_dt);
	}
	wsUint N_x;
	char *Ba_misPheno = getIO()->getPheMissing(&N_x);

	if (OPT_ENABLED(variantsummary)) for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) {
		Xa_msum[i].N_minPos = -1;
		Xa_msum[i].N_maxPos = -1;
	}

	for (wsUint i=0 ; i<N_mkr ; i++) {
		wsStrCst		Sp_FID	= NULL;
		int			N_pop	= -1;
		bool		B_fs	= true;
		bool		B_ps	= true;
		wsUint		N_cnt	= 0;
		xVariant&	X_m		= Xv_mkr[i];
		xMaf&		X_maf	= Xp_maf[i];

		/* Fisher */
		if (OPT_ENABLED(fisher)) {
			for (wsUint j=0 ; j<N_pheno ; j++)
				memset(Xa_fisher+j, 0x00, sizeof(xFisher));
			if (X_maf.N_allMac != 0) {
				for (wsUint j=0 ; j<N_samp ; j++) {
					if (Ba_misPheno[j]) continue;
					char N_geno = Na_data[j][i];
					if (isMissing(N_geno)) continue;

					wsUint N_ma = N_geno > 0;
					for (wsUint k=0 ; k<N_pheno ; k++)
						Xa_fisher[k].Na_freq[Ra_pheno[k][j]==WISARD_AFFECTED][N_ma]++;
				}
			}
			for (wsUint j=0 ; j<N_pheno ; j++) {
#define n Xa_fisher[j].Na_freq
				if (n[0][0] == 0 && n[1][0] == 0) {
					if (!OPT_ENABLED(remna)) {
						entryVariant(Cp_fe, X_m);
						Cp_fe->fmt("	%s	NA\n", Xv_pi[j].S_name.c_str());
					}
				} else {
					wsReal R_pval = fisher(n[0][0], n[0][1], n[1][0], n[1][1]);
					entryVariant(Cp_fe, X_m);

					/* p-value check */
					if (PVAL_FAIL(R_pval)) continue;
					Cp_fe->fmt("	%s	%g\n", Xv_pi[j].S_name.c_str(),
						R_pval);
				}
#undef n
			}
		}

		/* Cochran-Armitage trend test */
		if (OPT_ENABLED(trend)) {
			for (wsUint j=0 ; j<N_pheno ; j++)
				memset(Xa_trend+j, 0x00, sizeof(xTrend));
			wsReal Ra_t[3] = { 0, 1, 2 };
			for (wsUint j=0 ; j<N_samp ; j++) {
				if (Ba_misPheno[j]) continue;
				char N_geno = Na_data[j][i];
				if (isMissing(N_geno)) continue;

				for (wsUint k=0 ; k<N_pheno ; k++)
					Xa_trend[k].N_freq[Ra_pheno[k][j]==WISARD_AFFECTED][(wsUint)N_geno]++;
			}
			for (wsUint j=0 ; j<N_pheno ; j++) {
#define n Xa_trend[j].N_freq
#if 1
				///////////////////////////////////
				// Calculate association statistics
// 				double obs_A = n[1][0] + n[1][1] + n[1][2];
// 				double obs_U = n[0][0] + n[0][1] + n[0][2];
// 				double obs_T = obs_A + obs_U;
// 
// 				wsUint A11 = n[1][0];
// 				wsUint A12 = n[1][1];
// 				wsUint A22 = n[1][2];
// 				wsUint U11 = n[0][0];
// 				wsUint U12 = n[0][1];
// 				wsUint U22 = n[0][2];
// 
// 				double obs_1 = 2*(A11+U11) + A12 + U12;
// 				double obs_2 = 2*(A22+U22) + A12 + U12;
// 
// 				double obs_11 = n[0][0] + n[1][0];
// 				double obs_12 = n[0][1] + n[1][1];
// 				double obs_22 = n[0][2] + n[1][2];
// 
// 				///////////////////////
// 				// Cochram-Armitage Trend test
// 
// 				double CA = ( ( obs_U / obs_T * A12 ) - ( obs_A / obs_T * U12 ) )
// 					+ 2*( ( obs_U / obs_T * A22 ) - ( obs_A / obs_T * U22 ) ) ;
// 
// 				double varCA = obs_A * obs_U 
// 					* ( ( obs_T * ( obs_12 + 4*obs_22 ) 
// 					- ( obs_12+2*obs_22 ) * ( obs_12+2*obs_22 )  ) 
// 					/ (obs_T * obs_T * obs_T  )) ;
// 
// 				double CA_chisq = (CA*CA) / varCA;
// 				double CA_p = chiprobP(CA_chisq,1);
// 
// 				double mult_p, mult_chisq;
// 
// 				///////////////////////
// 				// Multiplicative model
// 
// 				double obs_A1 = 2*A11 + A12;
// 				double obs_A2 = 2*A22 + A12;
// 				double obs_U1 = 2*U11 + U12;
// 				double obs_U2 = 2*U22 + U12;
// 
// 				if (OPT_ENABLED(fisher))  {
// 					mult_p = fisher((int)obs_A1, (int)obs_U1, (int)obs_A2, (int)obs_U2);
// 				} else {
// 					double exp_A1 = (obs_A * obs_1 ) / obs_T;    // note 2's cancelled for obs_A and obs_T 
// 					double exp_A2 = (obs_A * obs_2 ) / obs_T;    // which are counts of individuals, not 
// 					double exp_U1 = (obs_U * obs_1 ) / obs_T;    // alleles
// 					double exp_U2 = (obs_U * obs_2 ) / obs_T;
// 
// 					mult_chisq =  ( ( obs_A1 - exp_A1 ) * ( obs_A1 - exp_A1 ) ) / exp_A1
// 						+ ( ( obs_A2 - exp_A2 ) * ( obs_A2 - exp_A2 ) ) / exp_A2
// 						+ ( ( obs_U1 - exp_U1 ) * ( obs_U1 - exp_U1 ) ) / exp_U1
// 						+ ( ( obs_U2 - exp_U2 ) * ( obs_U2 - exp_U2 ) ) / exp_U2;
// 
// 
// 					///////////////////////
// 					// Multiplicative model
// 
// 					mult_p = chiprobP(mult_chisq,1);	
// 
// 				}
// 
// 
// 				double gen_p, dom_p, rec_p;
// 				gen_p = dom_p = rec_p = -9;
// 				double dom_chisq, rec_chisq, gen_chisq;
// 
// 				//////////////////////////////////////////////////////////////
// 				// Standard chi-square test, or Fisher's exact
// 
// 				if (OPT_ENABLED(fisher)) {
// 					////////////
// 					// General
// 
// // 					table_t t;
// // 					sizeTable(t,3,2);
// // 					t[0][0] = A11;
// // 					t[1][0] = A12;
// // 					t[2][0] = A22;
// // 					t[0][1] = U11;
// // 					t[1][1] = U12;
// // 					t[2][1] = U22;
// 					gen_p = fisher(U11, U12, U22, A11, A12, A22);
// 
// 					////////////
// 					// Dominant
// 					dom_p = fisher(A11+A12, U11+U12, A22, U22);
// 
// 					/////////////
// 					// Recessive
// 					rec_p = fisher(A11, U11, A12+A22, U12+U22);
// 				} else {
// 					///////////////////////
// 					// General model
// 
// 					double exp_A11 = (obs_A * obs_11 ) / obs_T;
// 					double exp_A12 = (obs_A * obs_12 ) / obs_T;
// 					double exp_A22 = (obs_A * obs_22 ) / obs_T;
// 					double exp_U11 = (obs_U * obs_11 ) / obs_T;
// 					double exp_U12 = (obs_U * obs_12 ) / obs_T;
// 					double exp_U22 = (obs_U * obs_22 ) / obs_T;
// 
// 					gen_chisq =  ( ( A11 - exp_A11 ) * ( A11 - exp_A11 ) ) / exp_A11
// 						+ ( ( A12 - exp_A12 ) * ( A12 - exp_A12 ) ) / exp_A12
// 						+ ( ( A22 - exp_A22 ) * ( A22 - exp_A22 ) ) / exp_A22
// 						+ ( ( U11 - exp_U11 ) * ( U11 - exp_U11 ) ) / exp_U11
// 						+ ( ( U12 - exp_U12 ) * ( U12 - exp_U12 ) ) / exp_U12
// 						+ ( ( U22 - exp_U22 ) * ( U22 - exp_U22 ) ) / exp_U22;
// 
// 
// 					///////////////////////
// 					// Dominant (minor allele) (1) model 
// 
// 					dom_chisq =  ( ( (A11+A12) - (exp_A11+exp_A12) ) * ( (A11+A12) - (exp_A11+exp_A12) ) ) / (exp_A11+exp_A12) 
// 						+ ( ( A22 - exp_A22 ) * ( A22 - exp_A22 ) ) / exp_A22
// 						+ ( ( (U11+U12) - (exp_U11+exp_U12) ) * ( (U11+U12) - (exp_U11+exp_U12) ) ) / (exp_U11+exp_U12) 
// 						+ ( ( U22 - exp_U22 ) * ( U22 - exp_U22 ) ) / exp_U22;
// 
// 
// 					//////////////////////////////////////
// 					// Recessive (minor allele) (1) model 
// 
// 					rec_chisq =  ( ( (A22+A12) - (exp_A22+exp_A12) ) * ( (A22+A12) - (exp_A22+exp_A12) ) ) / (exp_A22+exp_A12) 
// 						+ ( ( A11 - exp_A11 ) * ( A11 - exp_A11 ) ) / exp_A11
// 						+ ( ( (U22+U12) - (exp_U22+exp_U12) ) * ( (U22+U12) - (exp_U22+exp_U12) ) ) / (exp_U22+exp_U12) 
// 						+ ( ( U11 - exp_U11 ) * ( U11 - exp_U11 ) ) / exp_U11;
// 
// 
// 					//////////////////////////////////
// 					// p-values and model comparisons 
// 
// 					gen_p = chiprobP(gen_chisq,2);
// 					dom_p = chiprobP(dom_chisq,1);
// 					rec_p = chiprobP(rec_chisq,1);
// 				}
#endif
				wsReal R_r1 = n[0][0] + n[0][1] + n[0][2];
				wsReal R_c1 = n[0][0] + n[1][0];
				wsReal R_c2 = n[0][1] + n[1][1];
				wsReal R_c3 = n[0][2] + n[1][2];

				if (R_c2 == 0 && R_c3 == 0) {
					if (OPT_ENABLED(remna)) continue;

					entryVariant(Cp_ca, X_m);
					Cp_ca->fmt("	%s	NA	NA\n", Xv_pi[j].S_name.c_str());
				} else {
					wsReal R_N = R_c1 + R_c2 + R_c3;
					wsReal R_r2 = R_N - R_r1;

					wsReal exp_A11 = (R_r2 * R_c1) / R_N;
					wsReal exp_A12 = (R_r2 * R_c2) / R_N;
					wsReal exp_A22 = (R_r2 * R_c3) / R_N;
					wsReal exp_U11 = (R_r1 * R_c1) / R_N;
					wsReal exp_U12 = (R_r1 * R_c2) / R_N;
					wsReal exp_U22 = (R_r1 * R_c3) / R_N;

					wsReal gen_chisq =  ( ( n[1][0] - exp_A11 ) * ( n[1][0] - exp_A11 ) ) / exp_A11
						+ ( ( n[1][1] - exp_A12 ) * ( n[1][1] - exp_A12 ) ) / exp_A12
						+ ( ( n[1][2] - exp_A22 ) * ( n[1][2] - exp_A22 ) ) / exp_A22
						+ ( ( n[0][0] - exp_U11 ) * ( n[0][0] - exp_U11 ) ) / exp_U11
						+ ( ( n[0][1] - exp_U12 ) * ( n[0][1] - exp_U12 ) ) / exp_U12
						+ ( ( n[0][2] - exp_U22 ) * ( n[0][2] - exp_U22 ) ) / exp_U22;
					wsReal R_pGeno = (R_r1 * R_r2 * R_c1 * R_c2 * R_c3 == W0) ? WISARD_NAN : PVchisq(gen_chisq, 2.0);

					wsReal R_bar = (R_c1*Ra_t[0] + R_c2*Ra_t[1] + R_c3*Ra_t[2]) / R_N;
					wsReal R_s1 = (Ra_t[0]-R_bar);
					wsReal R_s2 = (Ra_t[1]-R_bar);
					wsReal R_s3 = (Ra_t[2]-R_bar);
					wsReal R_ss = SQR(R_s1)*R_c1 + SQR(R_s2)*R_c2 + SQR(R_s3)*R_c3;

					wsReal R_phi = R_r1 / R_N;
					wsReal R_tTrend = (n[0][0]*R_s1 + n[0][1]*R_s2 + n[0][2]*R_s3);
					R_tTrend = SQR(R_tTrend) / (R_phi*(W1-R_phi)*R_ss);
					wsReal R_pTrend = PVchisq(R_tTrend, 1.0);

					//## pd is the estimated allele frequency in cases. ##
					wsReal pd = (2*n[0][2]+n[0][1])/(2*R_r1);
					//## ph is the estimated allele frequency in controls. ##
					wsReal ph = (2*n[1][2]+n[1][1])/(2*(R_N-R_r1));
					//## p is the estimated allele frequency under null hypothesis. ##
					wsReal p = (2*R_c3+R_c2)/(2*R_N);
					//## Calculate the test statistic. ##
					wsReal u = ph-pd;
					wsReal v = p*(1-p)*(1/(2*R_r1)+1/(2*(R_N-R_r1)));
					wsReal R_Tabt = SQR(u)/v;
					//## Calculate the p-value. ##
					wsReal R_pABT = PVchisq(R_Tabt, 1.0);

					/* If both p-values are fail */
					if (PVAL_FAIL3(R_pTrend, R_pABT, R_pGeno)) continue;
					entryVariant(Cp_ca, X_m);
					char S_pTrend[128];
					char S_pABT[128];
					char S_pGeno[128];
					if (NA(R_pTrend)) strcpy(S_pTrend, "NA"); else sprintf(S_pTrend, "%g", R_pTrend);
					if (NA(R_pABT)) strcpy(S_pABT, "NA"); else sprintf(S_pABT, "%g", R_pABT);
					if (NA(R_pGeno)) strcpy(S_pGeno, "NA"); else sprintf(S_pGeno, "%g", R_pGeno);
					Cp_ca->fmt("	%s	%s	%s	%s\n", Xv_pi[j].S_name.c_str(),
						S_pTrend, S_pABT, S_pGeno);
				}
#undef n
			}
		}

		/* --variantsummary */
		if (OPT_ENABLED(variantsummary)) {
			wsUint	N_cIdx	= X_m.chr > 0 && (wsUint)X_m.chr <= NCHR_SPECIES ? X_m.chr : 0;
			xVariantSumm
					*Xp_m	= Xa_msum + N_cIdx;

			if (OPT_ENABLED(indel)) {
				if (!X_m.indel2) {
				} else if (X_m.indel1[0] == '-') {
					/* Insertion? */
					wsUint N_lenIns = (wsUint)strlen(X_m.indel2);
					Xp_m->N_ins++;
					if (Xp_m->N_maxIns == 0)
						Xp_m->N_maxIns = Xp_m->N_minIns = N_lenIns;
					else if (Xp_m->N_maxIns < N_lenIns)
						Xp_m->N_maxIns = N_lenIns;
					else if (Xp_m->N_minIns > N_lenIns)
						Xp_m->N_minIns = N_lenIns;
					Xp_m->N_lenIns += N_lenIns;
				} else if (X_m.indel2[0] == '-') {
					/* Deletion? */
					wsUint N_lenDel = (wsUint)strlen(X_m.indel1);
					Xp_m->N_del++;
					if (Xp_m->N_maxDel == 0)
						Xp_m->N_maxDel = Xp_m->N_minDel = N_lenDel;
					else if (Xp_m->N_maxDel < N_lenDel)
						Xp_m->N_maxDel = N_lenDel;
					else if (Xp_m->N_minDel > N_lenDel)
						Xp_m->N_minDel = N_lenDel;
					Xp_m->N_lenDel += N_lenDel;
				} else if (!X_m.indel1[1] && !X_m.indel2[1]) {
					/* SNV */
					Xp_m->N_snv++;
				} else {
					/* Substitution */
					Xp_m->N_sub++;
				}
			} else
				/* Only SNVs... */
				Xp_m->N_snv++;

			/* pos */
			if (Xp_m->N_minPos == -1)
				Xp_m->N_minPos = Xp_m->N_maxPos = X_m.pos;
			else if (Xp_m->N_minPos > (int)X_m.pos)
				Xp_m->N_minPos = X_m.pos;
			else if (Xp_m->N_maxPos < (int)X_m.pos)
				Xp_m->N_maxPos = X_m.pos;
				
			/* w/ --variantsummary */
			switch (X_maf.N_allMac) {
			case 0:
				if (Cp_mt) {
					entryVariant(Cp_mt, X_m);
					Cp_mt->put("\n");
				}
				Xp_m->N_mono++;
				break;
			case 1:
				if (Cp_st) {
					entryVariant(Cp_st, X_m);
					Cp_st->put("\n");
				}
				Xp_m->N_single++;
				break;
			case 2:
				if (Cp_dt) {
					entryVariant(Cp_dt, X_m);
					Cp_dt->put("\n");
				}
				Xp_m->N_double++;
				break;
			}
		} else {
			/* w/o --variantsummary */
			if (X_maf.N_allMac == 0 && Cp_mt) {
				entryVariant(Cp_mt, X_m);
				Cp_mt->put("\n");
			}
			else if (X_maf.N_allMac == 1 && Cp_st) {
				entryVariant(Cp_st, X_m);
				Cp_st->put("\n");
			}
			else if (X_maf.N_allMac == 2 && Cp_dt) {
				entryVariant(Cp_dt, X_m);
				Cp_dt->put("\n");
			}
		}

		if (Cp_ps) {
			for (wsUint j=0 ; j<N_samp ; j++) {
				char N_geno = Na_data[j][i];
				if (isMissing(N_geno) || !N_geno) continue;

				/* Record its FID if NULL */
				if (N_pop == -1)
					N_pop = Xv_samp[j]->N_grpFst;
				/* Stop if both FID is DIFFERENT */
				else if (N_pop != Xv_samp[j]->N_grpFst) {
					B_ps = false;
					break;
				}
				N_cnt++;
			}
			/* Print result if this variant is population specific */
			if (B_ps) {
				entryVariant(Cp_ps, X_m);
				Cp_ps->fmt("	%s	%d\n", Cp_IO->getSampGrpIdx2grp()[N_pop].c_str(), N_cnt);
			}
		}

		if (Cp_fs) {
			for (wsUint j=0 ; j<N_samp ; j++) {
				char N_geno = Na_data[j][i];
				if (isMissing(N_geno) || !N_geno) continue;

				/* Record its FID if NULL */
				if (Sp_FID == NULL)
					Sp_FID = Xv_samp[j]->S_FID.c_str();
				/* Stop if both FID is DIFFERENT */
				else if (!stricmp(Sp_FID, Xv_samp[j]->S_FID.c_str())) {
					B_fs = false;
					break;
				}
				N_cnt++;
			}
			/* Print result if this variant is family specific */
			if (B_fs) {
				entryVariant(Cp_fs, X_m);
				Cp_fs->fmt("	%s	%d\n", Sp_FID, N_cnt);
			}
		}
	}

	/* --variantsummary */
	if (OPT_ENABLED(variantsummary)) {
		cExporter *Cp_ms = cExporter::summon("summary.variant.res");
		LOGoutput("Summary of variants is exported to [%s.summary.variant.res]\n",
			OPT_STRING(out));

		if (OPT_ENABLED(indel)) {
			Cp_ms->put("CHR	MINPOS	MAXPOS	NSNV	NMONO	NSINGLETON	"
				"NDOUBLETON	NSUB	NINS	MININS	MAXINS	AVGINS	NDEL"
				"	MINDEL	MAXDEL	AVGDEL\n");
			xVariantSumm X_overall = { 0, };
			for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) {
				xVariantSumm& X_m = Xa_msum[i];
				/* Print chromosome-wise */
				Cp_ms->fmt("%s	%d	%d	%d	%d	%d	%d	%d	%d	%d	%d	%g"
					"	%d	%d	%d	%g\n", getChrName2(i), X_m.N_minPos,
					X_m.N_maxPos, X_m.N_snv, X_m.N_mono, X_m.N_single,
					X_m.N_double, X_m.N_sub, X_m.N_ins, X_m.N_minIns,
					X_m.N_maxIns, X_m.N_lenIns/(wsReal)X_m.N_ins,
					X_m.N_del, X_m.N_minDel, X_m.N_maxDel,
					X_m.N_lenDel/(wsReal)X_m.N_del);
				X_overall.N_snv += X_m.N_snv;
				X_overall.N_mono += X_m.N_mono;
				X_overall.N_single += X_m.N_single;
				X_overall.N_double += X_m.N_double;
				X_overall.N_sub += X_m.N_sub;
				X_overall.N_ins += X_m.N_ins;
				X_overall.N_lenIns += X_m.N_lenIns;
				X_overall.N_del += X_m.N_del;
				X_overall.N_lenDel += X_m.N_lenDel;
				/* min/max update */
				if (i == 0) {
					X_overall.N_minIns = X_m.N_minIns;
					X_overall.N_maxIns = X_m.N_maxIns;
					X_overall.N_minDel = X_m.N_minDel;
					X_overall.N_maxDel = X_m.N_maxDel;
				} else {
					if (X_overall.N_minIns > X_m.N_minIns)
						X_overall.N_minIns = X_m.N_minIns;
					if (X_overall.N_maxIns > X_m.N_maxIns)
						X_overall.N_maxIns = X_m.N_maxIns;
					if (X_overall.N_minDel > X_m.N_minDel)
						X_overall.N_minDel = X_m.N_minDel;
					if (X_overall.N_maxDel > X_m.N_maxDel)
						X_overall.N_maxDel = X_m.N_maxDel;
				}
			}
			/* Print overall */
			Cp_ms->fmt("OVERALL	-	-	%d	%d	%d	%d	%d	%d	%d	%d	%g"
				"	%d	%d	%d	%g\n", X_overall.N_snv, X_overall.N_mono,
				X_overall.N_single, X_overall.N_double, X_overall.N_sub,
				X_overall.N_ins, X_overall.N_minIns, X_overall.N_maxIns,
				X_overall.N_lenIns/(wsReal)X_overall.N_ins,
				X_overall.N_del, X_overall.N_minDel, X_overall.N_maxDel,
				X_overall.N_lenDel/(wsReal)X_overall.N_del);
		} else {
			Cp_ms->put("CHR	MINPOS	MAXPOS	NSNV	NMONO	NSINGLETON	"
				"NDOUBLETON\n");
			/* Print chromosome-wise */
			xVariantSumm X_overall = { 0, };
			for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) {
				xVariantSumm& X_m = Xa_msum[i];
				Cp_ms->fmt("%s	%d	%d	%d	%d	%d	%d\n",
					getChrName2(i), X_m.N_minPos, X_m.N_maxPos, X_m.N_snv,
					X_m.N_single, X_m.N_double);
				X_overall.N_snv += X_m.N_snv;
				X_overall.N_mono += X_m.N_mono;
				X_overall.N_single += X_m.N_single;
				X_overall.N_double += X_m.N_double;
			}
			/* Print overall */
			Cp_ms->fmt("OVERALL	-	-	%d	%d	%d	%d\n", X_overall.N_snv,
				X_overall.N_mono, X_overall.N_single, X_overall.N_double);
		}
		delete Cp_ms;
	}
	delete Cp_fs;

	if (Cp_ps) delete Cp_ps;
	if (Cp_fe) delete Cp_fe;
	if (Cp_ca) delete Cp_ca;
}

#endif

} // End namespace ONETOOL
