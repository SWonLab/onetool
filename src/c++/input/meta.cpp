#include <math.h>
#include "input/meta.h"
#include "input/stream.h"
#include "utils/integrate.h"
#include "utils/matrix.h"
#include "analyses/gstest.h"
#include "utils/stat.h"
#define NMAX_METAFILE	16

namespace ONETOOL {

typedef map<string,vector<wsReal*> >	mMetaRes;
typedef mMetaRes::iterator				mMetaRes_it;
typedef map<string,vReal>				mMetaWgt;
typedef mMetaWgt::iterator				mMetaWgt_it;

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cMetaIO::cMetaIO(wsStrCst S_fn, char B_inDry)
{
	/* Divide S_fn into multiples */
	wsUint	N_file	= 0;
	char**	S_fns	= loadStringValues(S_fn, &N_file);
	if (N_file == 0) halt("No meta-analysis result was assigned!");
	else if (N_file > NMAX_METAFILE) halt("Too many meta-analysis results "
		"[%d found, %d maximum] were assigned!", N_file, NMAX_METAFILE);
	//char	B_comp	= !OPT_ENABLED(avail);

	/* Check each file path */
	cStrFile	Ha_fns[NMAX_METAFILE];
	char*		Sp_tmp = new char[65536];
	LOOP (i, N_file) {
		sprintf(Sp_tmp, "[%d]th meta-analysis result", i+1);
		Ha_fns[i].open(S_fns[i], Sp_tmp);
	}

	/* Now read the results */
	mDataIdx			Xm_cNames;
	vector<mDataIdx>	Xm_cIndices(N_file);
	LOOP (i, N_file) {
		mDataIdx Xm_curCnames;
		char *a = NULL, *b = Sp_tmp;
		/* Read the header */
		if (Ha_fns[i].gets(Sp_tmp, 65536) == NULL)
			halt("[%d]th meta-analysis result is empty!", i+1);
		wsUint N_idx = 0;
		/* Retrieve header columns */
		getString(&b, &a);
		if (!b[0]) halt("[%d]th meta-analysis result is empty!", i+1);
		do {
			/* Store index */
			Xm_cIndices[i].insert(make_pair(b, N_idx++));

			/* Count the column name */
			Xm_cNames[b]++;
			if (Xm_curCnames[b] == 1)
				halt("Column [%s] appeared twice in file [%s]", b, S_fns[i]);
			Xm_curCnames[b]++;
			/* Retrieve next */
			if (!a) break;
			b = a;
			getString(&b, &a);
		} while (b && b[0]);
	}
	wsStrCst	S_tgtCol	= NULL;
	/* FIXME : Should address the function of this variable */
//	char	B_wgtInv	= 0;	/* Weights should be inverted(1) or not(0) */
	wsStrCst	S_wgtCol	= NULL;
	vStr	Sa_pvCol;

	char	B_searchP_	= 1;
	char*	S_defPcol	= NULL;
	char	N_metaMode	= -1;
	
	/* --farvat */
	wsStrCst	S_fvCol		= NULL;	/* FV_VARIANTS */
	wsStrCst	S_fcCol		= NULL;	/* FV_COVAR */
	wsStrCst	S_fsCol		= NULL;	/* FV_STATS */

	/* Under --farvat,
	 * (1) P_ columns are not be searched
	 * (2) Only P_CALPHA will be regarded
	 */
	if (OPT_ENABLED(farvat)) {
		B_searchP_ = 0;
		S_defPcol = strdup("P");
	} else
		S_defPcol = strdup("P");

	/* Which columns are appeared N_file times? */
	FOREACH (mDataIdx_it, Xm_cNames, i) {
		if (i->second != (int)N_file) continue;

		wsStrCst S_col = i->first.c_str();
// 		LOG("Column [%s] are appeared commonly in all meta-analysis files",
// 			S_col);

		/* Automatically finds required columns */
		if (!stricmp(S_col, "SNP") || !stricmp(S_col, "VARIANT")) {
			if (S_tgtCol)
				halt("Target column [%s] already exists but [%s] "
					"additionally found", S_tgtCol, S_col);
			LOG("Variant column [%s] found\n", S_col);
			S_tgtCol	= S_col;
			N_metaMode	= 0;
		} else if (!stricmp(S_col, "GENE")) {
			if (S_tgtCol)
				halt("Target column [%s] already exists but [%s] "
					"additionally found", S_tgtCol, S_col);
			LOG("Gene column [%s] found\n", S_col);
			S_tgtCol	= S_col;
			N_metaMode	= 1;
		} else if (!stricmp(S_col, "SE") || !stricmp(S_col, "WEIGHT")) {
			if (S_wgtCol)
				halt("Weight column [%s] already exists but [%s] "
					"additionally found", S_wgtCol, S_col);
//			B_wgtInv = 1;
			S_wgtCol = S_col;
		} else if (!memcmp(S_col, "P_", 2) && B_searchP_) {
			LOG("p-value column of [%s] found\n", S_col);
			Sa_pvCol.push_back(S_col);
		} else if (!stricmp(S_col, S_defPcol)) { /* Default p-value column */
			LOG("p-value column [%s] found\n", S_col);
			Sa_pvCol.push_back(S_col);
		}
		/* --farvat */
		else if (!stricmp(S_col, "FV_VARIANTS")) {
			S_fvCol = S_col;
		} else if (!stricmp(S_col, "FV_COVAR")) {
			S_fcCol = S_col;
		} else if (!stricmp(S_col, "FV_STATS")) {
			S_fsCol = S_col;
		}
	}
	/* Free unnecessary S_defPcol */
	free(S_defPcol);

	/* Check target column */
	if (!S_tgtCol)
		halt("Meta-analysis target column not found automatically!");
	/* Note for current analysis */
	LOGnote("Meta-analysis for [%d] tests against [%d] studies will be performed\n",
		Sa_pvCol.size(), N_file);
	/* If --farvat but either S_fvCol/S_fcCol/S_fsCol empty, do error */
	if (OPT_ENABLED(farvat)) {
		if (!S_fcCol || !S_fvCol || !S_fsCol)
			halt("Meta-FARVAT assigned but either FV_VARIANTS/FV_COVAR/FV_STATS is missing!");
		if (!N_metaMode)
			halt("Meta-FARVAT requires GENE field!");
	} else {
		if (Sa_pvCol.size() == 0)
			halt("Meta-analysis p-value column not found automatically!");
		/* If weight column, do Stouffer's Z-score method additionally */
		if (S_wgtCol)
			LOGnote("Stouffer's Z-score method is also performed\n");

	}
	/* B_isGene fixed? */
	if (N_metaMode == -1)
		halt("There is no target column!");

	/* Retrieve */
	mDataIdx Xm_vNames;
	LOOP (i, N_file) {
		wsUint N_tIdx = Xm_cIndices[i][S_tgtCol];
		wsUint L = 2;

		for (L=2 ; Ha_fns[i].gets(Sp_tmp, 65536) ; L++) {
			char *a = Sp_tmp, *b = NULL;
			/* Skip */
			LOOP (j, N_tIdx) {
				if (!a) halt("Line [%d] in file [%s] have incomplete data!",
					L, S_fns[i]);
				getString(&a, &b);
				a = b;
			}
			getString(&a, &b);
			Xm_vNames[a]++;
		}
		LOG("[%d] entries found in file [%s]\n", L-2, S_fns[i]);
	}

	/* Check out which are complete? */
	LOG("[%d] unique entries found\n", Xm_vNames.size());
	for (mDataIdx_it i=Xm_vNames.begin() ; i!=Xm_vNames.end() ; ) {
		if (i->second != (int)N_file)
			Xm_vNames.erase(i++);
		else
			++i;
	}

	if (Xm_vNames.size() == 0)
		halt("No complete entry found, cannot perform meta-analysis!");
	LOG("[%d] complete entries found\n", Xm_vNames.size());

	/* Retrieve again */
	LOOP (i, N_file) {
		Ha_fns[i].rewind();
		Ha_fns[i].gets(Sp_tmp, 65536);
	}

	/* Now perform meta analysis */
	LOGoutput("Meta-analysis result is exported to [%s.meta.res]\n",
		OPT_STRING(out));

	switch (N_metaMode) {
	case 0: /* For SNV */ {
		mMetaRes	Xm_pvals;
		mMetaWgt	Xm_wgts;
		LOOP (i, N_file) {
			int	N_tIdx = Xm_cIndices[i][S_tgtCol];
			int	N_wIdx = -1;
			if (S_wgtCol)
				Xm_cIndices[i][S_wgtCol];
			vInt	Na_pIdx;
			FOREACH (vStr_it, Sa_pvCol, j)
				Na_pIdx.push_back(Xm_cIndices[i][j->c_str()]);

			for (wsUint L=2 ; Ha_fns[i].gets(Sp_tmp, 65536) ; L++) {
				vStr	Sa_vals;
				for (char *a=Sp_tmp,*b=NULL ; a ; a=b) {
					getString(&a, &b);
					Sa_vals.push_back(a);
				}
				wsStrCst	S_tName	= Sa_vals[N_tIdx].c_str();
				/* Skip if not a target */
				mDataIdx_it X_find = Xm_vNames.find(S_tName);
				if (X_find == Xm_vNames.end()) continue;

				/* Find meta result */
				mMetaRes_it X_find2 = Xm_pvals.find(S_tName);
				mMetaWgt_it X_find3 = Xm_wgts.find(S_tName);
				if (X_find2 == Xm_pvals.end()) {
					/* Create new one */
					vector<wsReal*> V1;
					LOOP (j, Na_pIdx.size()) {
						V1.push_back(sseVector(N_file));
					}
					Xm_pvals.insert(make_pair(S_tName, V1));
					if (S_wgtCol)
						Xm_wgts.insert(make_pair(S_tName, vReal()));
					X_find2 = Xm_pvals.find(S_tName);
					X_find3 = Xm_wgts.find(S_tName);
				}

				/* Insert value */
				wsUint J = 0;
				FOREACHDO (vInt_it, Na_pIdx, j, J++) {
					wsReal R_pv = atof(Sa_vals[*j].c_str());
					X_find2->second[J][i] = R_pv;
				}

				/* Insert weight */
				if (S_wgtCol)
					X_find3->second.push_back(atof(Sa_vals[N_wIdx].c_str()));
			}
		}

		/* Export header */
		vStr	Sv_cols;
		string	Sv_params;
		Sv_cols.push_back("TARGET");
		Sv_params += "s";
		wsUint J = 0;
		char	S_colBuf[512];
		FOREACHDO (vStr_it, Sa_pvCol, j, J++) {
			sprintf(S_colBuf, "P_FISHER_%s", j->c_str());
			Sv_cols.push_back(S_colBuf);
			Sv_params += "r";
			if (S_wgtCol) {
				sprintf(S_colBuf, "P_STOUFFER_%s", j->c_str());
				Sv_cols.push_back(S_colBuf);
				Sv_params += "r";
			}
		}
		/* Initialize exporter */
		cTableExporter C_meta("meta.res", Sv_cols, Sv_params.c_str(),
			"Meta-analysis result");

		/* Now perform */
		double R_df = N_file * 2.0;
		mMetaWgt_it ii = Xm_wgts.begin();
		FOREACH (mMetaRes_it, Xm_pvals, i) {
			vReal	Xv_0;
			vReal&	Xv_wgts	= ii==Xm_wgts.end()?Xv_0:ii->second;
			wsStrCst	S_name	= i->first.c_str();
			vector<wsReal*>&
					Ra_pvs = i->second;
		
			C_meta.put(1, S_name);
			wsUint J = 0;
			FOREACHDO (vector<wsReal*>::iterator, Ra_pvs, j, J++) {
				/* Perform Fisher's method */
				wsReal R_sum = W0;
				LOOP (k, N_file) R_sum += log((*j)[k]);
				R_sum *= REAL_CONST(-2.0);
				C_meta.put(1, PVchisq(R_sum, R_df));

				/* FIXME : Perform Stouffer's method */
				if (S_wgtCol) {
					wsReal R_denom = W0;
					wsReal R_numer = W0;
					LOOP (k, N_file) {
						wsReal R_invNorm = qnorm(W1 - (*j)[k]);
						R_numer += Xv_wgts[k] * R_invNorm;
						R_denom += SQR(Xv_wgts[k]);
					}
					C_meta.put(1, PVchisq(R_numer/R_denom, 1.0));
				}
				C_meta.next();

				sseFree(*j);
			}

			if (S_wgtCol) ii++;
		}
	} break;
	case 1: { /* For FARVAT */
typedef map<string,mDataIdx>	mVarCnt;
typedef mVarCnt::iterator		mVarCnt_it;

		/* Init exporter */
		cTableExporter C_meta("meta.res", "srrr", "FARVAT meta analysis result",
			0, 4, "TARGET", "P_METAFARVAT_BURDEN", "P_METAFARVAT_VARCOMP",
			"P_METAFARVAT_OPTIMAL");

		mVarCnt								Xm_varCount;
		map<string,vector<mDataIdx> >		Xm_varPos;
		map<string,vector<wsVec> >			Xm_varStat;
		map<string,vector<cSymMatrix*> >	Xm_varCovar;

		LOOP (i, N_file) {
			int	N_varIdx	= Xm_cIndices[i][S_fvCol];
			int	N_sttIdx	= Xm_cIndices[i][S_fsCol];
			int	N_covIdx	= Xm_cIndices[i][S_fcCol];
			int	N_tIdx		= Xm_cIndices[i][S_tgtCol];
			for (wsUint L=2 ; Ha_fns[i].gets(Sp_tmp, 65536) ; L++) {
				vStr	Sa_vals;
				for (char *a=Sp_tmp,*b=NULL ; a ; a=b) {
					getString(&a, &b);
					Sa_vals.push_back(a);
				}
				string		S_gene		= Sa_vals[N_tIdx];
				mDataIdx&	Xm_curVC	= Xm_varCount[S_gene];
				vector<mDataIdx>&
							Xv_curVP	= Xm_varPos[S_gene];

				/* Get the number of variants in the gene at current file */
				wsUint		N_curVar = 0, N_curStat = 0, N_curCov = 0;
				char**		S_vals = loadStringValues(Sa_vals[N_varIdx].c_str(), &N_curVar);
				mDataIdx	Xm_curVP;
				LOOP (j, N_curVar) {
					char* S_cur = S_vals[j];

					Xm_curVP[S_cur] = (char)j;
					Xm_curVC[S_cur]++;
					free(S_cur);
				}
				DEALLOC(S_vals);
				Xv_curVP.push_back(Xm_curVP);

				/* Get the statistics and matrices */
				Xm_varStat[S_gene].push_back(
					loadRealValues(Sa_vals[N_sttIdx].c_str(), &N_curStat));
				wsVec Ra_covs = loadRealValues(Sa_vals[N_covIdx].c_str(), &N_curCov);
				vector<cSymMatrix*>& Xv_curCovar = Xm_varCovar[S_gene];
				if (Xv_curCovar.size() == 0)
					Xv_curCovar.resize(N_file);
				Xv_curCovar[i] = new cSymMatrix(Ra_covs, N_curStat);
				sseFree(Ra_covs);
			}
		}

		/* Removal threshold */
		char N_thres = (char)N_file;

		/* For each gene */
		FOREACH (mVarCnt_it, Xm_varCount, j) {
			/* Gene name */
			string S_gene = j->first;

			/* Position info, stats and covars for this gene */
			vector<mDataIdx>&		Xv_curVP	= Xm_varPos[S_gene];
			vector<wsVec>&			Xv_curVS	= Xm_varStat[S_gene];
			vector<cSymMatrix*>&	Xv_curVC	= Xm_varCovar[S_gene];

			/* Filtering out the variants */
			vStr Xv_passedVars;
			FOREACH (mDataIdx_it, j->second, i)
				if (i->second >= N_thres)
					Xv_passedVars.push_back(i->first);
			wsUint N_passedVars = (wsUint)Xv_passedVars.size();

			/* Init fin data */
			wsVec	Ra_U	= sseEmptyVec(N_passedVars);
			wsSym	Ra_S	= sseEmptySymMat(N_passedVars);
			wsVec	Ra_W	= sseEmptyVec(N_passedVars);
			sseVinit(Ra_W, N_passedVars, W1);
			cVector		V_U(Ra_U, N_passedVars);
			cDiagMatrix	M_W(N_passedVars, Ra_W);
			cVector		V_W(Ra_W, N_passedVars, 1);
			cSymMatrix	M_S(Ra_S, N_passedVars);

			LOOP (i, N_file) {
				vInt Xv_passedVarIdx;

				/* Check its position */
				FOREACH (vStr_it, Xv_passedVars, k)
					Xv_passedVarIdx.push_back(Xv_curVP[i][*k]);

				/* Resize Xm_varStat and Xm_varCovar and sum up them */
				wsUint	K = 0;
				wsSym	Ra_VC = Xv_curVC[i]->get();
				FOREACHDO (vInt_it, Xv_passedVarIdx, k, K++) {
					Ra_U[K] += Xv_curVS[i][*k];
					wsUint M = 0;
					wsVec Ra_co = Ra_S[K];

					FOREACHDO (vInt_it, Xv_passedVarIdx, l, M++)
						if (M <= K) {
							if (*k > *l)
								Ra_co[M] += Ra_VC[*k][*l];
							else
								Ra_co[M] += Ra_VC[*l][*k];
						}
				}
			}

			/* Get S0.5 */
			cSymMatrix& M_Ssq = M_S.sqrt();

			/* Burden */
			cVector V_UW	= V_U * V_W;
			wsReal	R_Tbdn	= SQR(V_UW.sum());

			cVector V_SsqW	= M_Ssq * V_W;
			wsRealCst	R_df	= V_SsqW.ss();

			double	R_Pbdn	= PVchisq(R_Tbdn, R_df);

			/* Calpha */
			cVector		V_UW2	= V_U * M_W;
			cStdMatrix	M_SsqW	= M_Ssq * M_W;
			delete &M_Ssq;
			cSymMatrix	M_SWWS	= M_SsqW.Mt();
			cVector		V_ev	= M_SWWS.eigen();
			wsReal		R_Tvc	= V_UW2.ss();
			wsReal		R_Pvc	= (wsReal)davies(R_Tvc, V_ev.get(), V_ev.size());
//			wsUint		B_isLiu	= 0;
			if (R_Pvc <= W0) {
				R_Pvc = liuEV(R_Tvc, V_ev.get(), V_ev.size());
//				B_isLiu = 1;
			}
			V_ev.rem();

			/* SKAT-o */
			wsUint N_rho = 0;
			wsVec Ra_rhos = getOptimalRhos(N_rho);

			wsReal		R_minPv		= W2;
			xLiuCov**	Xa_liu		= NULL; /* Checked */
//			wsReal		R_rhoSKATO	= WISARD_NAN;
			wsAlloc(Xa_liu, xLiuCov*, N_rho);
			//wsUint N_idxMin = 0;

			/* For each rho, compute statistics and p-value */
			LOOP(i, N_rho) {
				cSymMatrix M_R(W1, Ra_rhos[i], N_passedVars);
				//REAL_c R_Tstat = V_UW2.qf(M_R);
				Xa_liu[i] = liuCov(M_S.get(), N_passedVars, M_R.get());

				/* Other case */
				Xa_liu[i]->R_stat = (R_Tvc*(W1-Ra_rhos[i]) + R_Tbdn*Ra_rhos[i]);

				/* Derive p-value */
				// pvals[jj] <- pliu(stats[jj],lius[[jj]])	
				double R_nc = (double)(Xa_liu[i]->R_delta);
				Xa_liu[i]->R_pVal = (wsReal)PVchisq((Xa_liu[i]->R_stat -
					Xa_liu[i]->R_muQ) / Xa_liu[i]->R_sigQ * Xa_liu[i]->R_sigX
					+ Xa_liu[i]->R_muX, Xa_liu[i]->R_l, &R_nc);

				/* Minimum p-value check */
				if (Xa_liu[i]->R_pVal < R_minPv) {
					R_minPv		= Xa_liu[i]->R_pVal;
//					R_rhoSKATO	= Ra_rhos[i];
				}
			}

			/* Z = M_SsqW */
			// Z2    <- Z%*%t(Z)
			// 
			wsUint N_szLiu = N_passedVars;
			cSymMatrix M_ZZt = M_SsqW.Mt();
			wsReal R_sumZbarZ = W0;
			cVector V_zBar = M_SsqW.meanR();
			wsVec  Ra_zBar = V_zBar.get();
			wsRealCst R_z2 = V_zBar.ss();
			cVector V_ZbarZ = V_zBar * M_SsqW;

			// ssZbZ <- sum(ZbZ^2)
			R_sumZbarZ = V_ZbarZ.ss();
			V_ZbarZ.rem();

			wsSym Ra_IM = sseSymMat(N_szLiu);
			LOOP(i, N_szLiu) {
				wsReal R_mul = -Ra_zBar[i] / R_z2;
				sseVpC(Ra_zBar, R_mul, Ra_IM[i], i);

				wsReal R_ii = -Ra_zBar[i]*R_mul;
				Ra_IM[i][i] = W1 - R_ii;
			}
			V_zBar.rem();
			cSymMatrix M_IM(Ra_IM, N_szLiu);
			cStdMatrix	M_IMZ;
			wsMat		Ra_IMV = NULL;
			wsMat		Ra_IMV2 = NULL;
			/* IMV = (I-M) ZZ' */
			cStdMatrix M_IMV = M_IM * M_ZZt;
			M_IMV.setClean();
			Ra_IMV = M_IMV.get();

			cStdMatrix M_Zt;
			cSymMatrix M_ZtIMZ;
			// IMV2 <- IMV%*%IMV
			/* GetIMV2 = non-p^3 version */
			//cSymMatrix M_IMV2 = M_IMV.Mt();
			cStdMatrix M_IMV2 = M_IMV * M_IMV;
			M_IMV2.setClean();
			Ra_IMV2 = M_IMV2.get();

			// muQ  <- sum(diag(IMV))
			wsReal R_muQ = W0;
			wsReal R_sumIMV2 = W0;
			LOOP (i, N_szLiu) {
				R_muQ += Ra_IMV[i][i];
				R_sumIMV2 += Ra_IMV2[i][i];
			}

			//		if (!B_isCont)
			//			R_muQ = R_sumIMV2;

			// rRho <- nSNP^2*rhos*z2 + (1-rhos)/z2 * ssZbZ
			wsVec Ra_rRho = sseVector(N_rho);
			LOOP(i, N_rho)
				Ra_rRho[i] = SQR(N_passedVars)*Ra_rhos[i]*R_z2
					+ (W1-Ra_rhos[i])/R_z2 * R_sumZbarZ;

			// t(Z)%*%M%*%Z%*%t(Z)%*%(Imq-M)%*%Z)
			wsReal R_dsum = WISARD_NAN;
			// Z2IMV <- Z2-IMV
			//			cSymMatrix M_PsIMV = M_ZZt - M_IMV;
			wsMat Ra_PsIMV = sseMatrix(N_szLiu, N_szLiu); /* Checked */
			sseSsM(M_ZZt.get(), Ra_IMV, Ra_PsIMV, N_szLiu);

			// sigXi <- 2 * sqrt( sum(diag( Z2IMV%*%IMV )) )
			//wsReal R_dsum = M_PsIMV.tr(M_IMV);
			R_dsum = diagSum2matrix(Ra_PsIMV, N_szLiu, N_szLiu,
				Ra_IMV, N_szLiu, N_szLiu);
			sseUnmat(Ra_PsIMV, N_szLiu);
			if (R_dsum < W0) /* Possible to down to 0 from round error */
				R_dsum = W0;
			wsReal R_sigXi = W2 * sqrt(R_dsum);
			//sseUnmat(Ra_PsIMV, N_szLiu);

			// sigQ <- sqrt( 2*sum(diag(IMV2)) + sigXi^2 )
			wsReal R_sigQ = sqrt(W2*R_sumIMV2 + SQR(R_sigXi));
			// id   <- order(pvals)[1]
			// rhoMin <- rhos[id]
			// pMin   <- pvals[id]
			// for(jj in seq(along=rhos)) {
			wsReal *Ra_qLiu = NULL; /* Checked */
			wsAlloc(Ra_qLiu, wsReal, N_rho);
			LOOP(i, N_rho) {
				// rho <- rhos[jj]
				// qmins[jj] <- qliu(pMin,lius[[jj]])
				Ra_qLiu[i] = qLiu(R_minPv, Xa_liu[i]);
				DEALLOC(Xa_liu[i]);
			}
			DEALLOC(Xa_liu);

			// lius      <- getliu(covs,I - M)
			//[121210] lius      <- getliu(Z2,I - M)
			xLiuCov *Xp_liuFinal = NULL;
			Xp_liuFinal = liuCov(NULL, N_szLiu, NULL, Ra_IMV, Ra_IMV2);
			sseUnmat(Ra_IMV, N_szLiu);
			sseUnmat(Ra_IMV2, N_szLiu);

			// 	pvalSKATO <- 1 - unlist(integrate( integrateSKATo, 0, Inf, qmins=qmins, rRho=rRho,
			// 	rhos=rhos, muQ=muQ, sigQ=sigQ, sigXi=sigXi,lius=lius))$value
			xIntgSKATO X_par;
			X_par.R_muQ			= R_muQ;
			X_par.R_sigQ		= R_sigQ;
			X_par.R_sigXi		= R_sigXi;
			X_par.Ra_qLiu		= Ra_qLiu;
			X_par.Ra_rho		= Ra_rhos;
			X_par.Ra_rRho		= Ra_rRho;
			X_par.Xp_liuFinal	= Xp_liuFinal;
			X_par.N_rho			= N_rho;
			xIntgInp X_inp(integrateOptimal, 0, 40.0, &X_par, 500, pow(numeric_limits<double>::epsilon(), 0.25)); // 140507 update
			xIntgRes *Xp_res = integrate(&X_inp); /* Checked */
			//	LOG("%g %d\n", Xp_res->R_value, Xp_res->N_err);
			if (NA(Xp_res->R_value) || Xp_res->R_value == 0 || Xp_res->N_err || Xp_res->R_value > 0.99999) {
				xIntgInp X_inp2(integrateOptimal, 0, 40.0, &X_par, 2000, 1e-25); // 140507 update
				//		DEALLOC(Xp_res);
				//		LOG("Do compensate...\n");
				xIntgRes *Xp_res2 = integrate(&X_inp2); /* Checked */
				if (Xp_res2->R_value == 0 || NA(Xp_res2->R_value)) {
					xIntgInp X_inp3(integrateOptimal, 0, numeric_limits<double>::infinity(), &X_par, 2000, 1e-25); // 140507 update
					//		DEALLOC(Xp_res);
					//		LOG("Do compensate...\n");
					xIntgRes *Xp_res3 = integrate(&X_inp3); /* Checked */
					if (Xp_res3->R_value == 0 || NA(Xp_res3->R_value)) {
						DEALLOC(Xp_res3);
					} else {
						DEALLOC(Xp_res);
						Xp_res = Xp_res3;
					}
					DEALLOC(Xp_res2);
				} else {
					DEALLOC(Xp_res);
					Xp_res = Xp_res2;
				}
			}
			sseFree(Ra_rRho);
			DEALLOC(Ra_qLiu);
			DEALLOC(Xp_liuFinal);
			sseFree(Ra_rhos);

			//	if (Xp_res->N_err)
			//		*Rp_pValSKATO = WISARD_NA_REAL;
			//	else
			wsReal R_Popt = (wsReal)(1.0 - Xp_res->R_value);
			DEALLOC(Xp_res);

			/* Write out */
			C_meta.write(4, S_gene.c_str(), R_Pbdn, R_Pvc, R_Popt);

			FOREACH (vector<wsVec>::iterator, Xv_curVS, i)
				sseFree(*i);
			FOREACH (vector<cSymMatrix*>::iterator, Xv_curVC, i)
				delete *i;
		}




	} break;
	} /* END OF switch */

	delete [] Sp_tmp;
	LOOP (i, N_file) free(S_fns[i]);
	DEALLOC(S_fns);
}
 
#endif

} // End namespace ONETOOL
