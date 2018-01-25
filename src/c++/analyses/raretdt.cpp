#include <algorithm>
#include "utils/util.h"
#include "analyses/raretdt.h"
#include "analyses/fam.h"
#include "analyses/setmgr.h"
#include "utils/dcdflib.h"
//#include "unphasedTDT.h"

#define PR(a) std::cout << #a " = " << a << std::endl;

namespace ONETOOL {

inline float   abs1(float x)    {return (x>0.0 ? x :-x);}
inline int      abs1(int x)       {return (x>0 ? x :-x);}

/* 140603 From FB-SKAT */
wsReal expected_geno_offspring(int ch, int p1, int p2)
{
	if ((p1 == p2) && (abs1(ch - p1) == 2)) return(0);
	else {
		if((p1 == 1) && (p2 == 1))
			return(pow(0.5,(abs1(ch - 1) + 1)));
		if(abs1(p1 - p2) == 2) {
			if(ch == 1)
				return W1;
			else return W0;
		}
		if((p1 == p2) && (p1 != 1)) {
			if(ch == p1)
				return W1;
			else return W0;
		}
		if(abs1(p1 - p2) == 1) {
			if(ch == p1)
				return(0.5);
			if(ch == p2)
				return(0.5);
			return W0;
		}
	}
	halt("SYSERR : This should not be happen!");
	return WISARD_NAN;
}

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

void cRvtdtAnalysis::m_preprocess(vInt& Xv_set)
{
	mNucFam&	Xm_nf	= Cp_anaFS->getNucFamData();
	char**		Na_geno	= getIO()->getGenotype();
	wsUint		N_var	= (wsUint)Xv_set.size();
// 	if (__genotypes.size()%__rowsPerTrio != 0) {
// 		__informGene=false;
// 		return;
// 	}

	/* Get the number of *trios* */
	__nInformTrios = 0;//__genotypes.size()/__rowsPerTrio;
	FOREACH (mNucFam_it, Xm_nf, i) {
		xNucFamily &X_fam = i->second;
		/* For each child, we can construct trio */
		FOREACH (vSampPtr_it, X_fam.Xp_childs, j) {
			/* FB-SKAT part */
			if (OPT_ENABLED(fbskat)) {
				for (wsUint p1=0 ; p1<=2 ; p1++)
					for (wsUint p2=0 ; p2<=2 ; p2++)
						Echp1p2[p1][p2] = W0;

				for (wsUint p1=0;p1<=2;p1++)
					for (wsUint p2=0;p2<=2;p2++)
						for (wsUint m=0;m<=2;m++)
							Echp1p2[p1][p2] += m*expected_geno_offspring(m,p1,p2);
			}

			__nInformTrios++;
		}
	}

	__missingRatio.assign(N_var, 0);
	__samMafs.assign(N_var, 0);
	vector<wsUint> parentNotMiss(N_var, 0); //need this to calculate sample maf
//	double maxGenos = W2;
	vBool varInParent(N_var, false), varInKid(N_var, false);

	// For each trio
	FOREACH (mNucFam_it, Xm_nf, i) {
		xNucFamily&	X_fam	= i->second;
		vInt		Xv_idx;

		/* Construct map */
		Xv_idx.push_back(X_fam.Xp_pat->N_idx);
		Xv_idx.push_back(X_fam.Xp_mat->N_idx);
		FOREACH (vSampPtr_it, X_fam.Xp_childs, j)
			Xv_idx.push_back((*j)->N_idx);
		/* Do not consider if the family does not have valid child */
		if (Xv_idx.size() == 2) continue;

		wsUint J = 0;
		FOREACHDO (vInt_it, Xv_idx, j, J++) {
			wsUint K = 0;
			FOREACHDO (vInt_it, Xv_set, k, K++) {
				char N_geno = Na_geno[*j][*k];
				if (isMissing(N_geno)) {
					__missingRatio[K]++;
					continue;
				}

				if (J < 2) {
					__samMafs[K] += N_geno; //only consider parents when calculating sample maf
					parentNotMiss[K]++;
					varInParent[K] = true;
				} else {
					varInKid[K] = true;
				}
			}
		}
	}

	LOOP (i, N_var) {
		__samMafs[i] = (__samMafs[i]==0)?1:__samMafs[i]; // in case zero 
		__samMafs[i] = __samMafs[i]/float(parentNotMiss[i]*2);
		__missingRatio[i] = __missingRatio[i] / (wsReal)(Xv_set.size());
	}

	// flip genotype if samMaf< 0.5
// 	bool flip(false);
// 	vectorL needFlip(N_var, false);
// 	LOOP (i, N_var) {
// 		if (__samMafs[i]>0.5) {flip=true; needFlip[i]=true;}
// 	}
// 	if (flip){
// 		LOOP (i, N_var) {
// 			if (!needFlip[i]) continue;
// 
// 			for (wsUint x=0; x<Na_geno.size(); x++) {
// 				if (Na_geno[x][i] >= 0 && Na_geno[x][i] <= maxGenos) Na_geno[x][i]=maxGenos-Na_geno[x][i];
// 			}
// 		}
// 	}

	// which variant sites to be analyzed
	__varToBeAnalyzed.assign(N_var, true);
	if (__missCutoff<W1 && __missCutoff>W0)
		for (wsUint i = 0; i < __varToBeAnalyzed.size(); i++)
			if (__missingRatio[i] > __missCutoff)
				__varToBeAnalyzed[i]=false;

	// write log Json
// 	Json::Value vStatic;
// 	vStatic["populationMafs"] = buildJsonArray(__popMafs);
// 	vStatic["sampleMafs"] = buildJsonArray(__samMafs);
// 	vStatic["missingRatio"] = buildJsonArray(__missingRatio);
// 	vStatic["varFoundInParent"] = buildJsonArray(varInParent);
// 	vStatic["varFoundInKid"] = buildJsonArray(varInKid);
// 	logRoot["variantStatic"] = vStatic;
}


cRvtdtAnalysis::cRvtdtAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSM,
	cFamStrAnalysis *Cp_inpAnaFS, bool B_inpPhased, const vReal& popMafs, bool skipMiss, double missCutoff)
	: cAnalysis(Cp_inpIO)
{
	Cp_anaFS = Cp_inpAnaFS;
	Cp_anaSM = Cp_inpAnaSM;
	__skipMiss = skipMiss;
	__missCutoff = missCutoff;

	/* Set default option value */
	if (!IS_ASSIGNED(nperm)) {
		OPTION().assign("nperm",OPTION().getDefVal("nperm"), 1);
		OPTION().FORCE_OPT_NUMBER(nperm);
	}
	// 	if (__popMafs.size() != __genotypes[0].size()) {
// 		__informGene=false;
// 		return;
// 	}
	/* Type checking */
	if (getIO()->sizePheno() != 1 ||
		getIO()->isContinuous())
		halt("--rarefam/--fbskat only applicable to single dichotomous "
			"phenotype");
}

cRvtdtAnalysis::cRvtdtAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSM,
	cFamStrAnalysis *Cp_inpAnaFS, bool phased, bool skipMiss, double missCutoff) :
	cAnalysis(Cp_inpIO)
{
	Cp_anaFS = Cp_inpAnaFS;
	Cp_anaSM = Cp_inpAnaSM;
	__skipMiss = skipMiss;
	__missCutoff = missCutoff;

	/* Set default option value */
	if (!IS_ASSIGNED(nperm)) {
		OPTION().assign("nperm",OPTION().getDefVal("nperm"), 1);
		OPTION().FORCE_OPT_NUMBER(nperm);
	}
	// 	if (__popMafs.size() != __genotypes[0].size()) {
	// 		__informGene=false;
	// 		return;
	// 	}
	/* Type checking */
	if (getIO()->sizePheno() != 1 ||
		getIO()->isContinuous())
		halt("--rarefam/--fbskat only applicable to single dichotomous "
		"phenotype");
}

inline double skatPval(double v1, double v2){
	// y1=1-pchisq(x[,4],x[,5]);
	double res = (NA(v1) || NA(v2) || v1<0 || v2<0) ? WISARD_NAN : PVchisq(v1, v2);
	return res;
}

void cRvtdtAnalysis::run()
{
	/* Check set */
	if (!Cp_anaSM)
		halt("rvTDT requires either --set or --setconsec!");

	/* Prepare exporter */
	cExporter *Cp_rt = NULL;
	cExporter *Cp_fs = NULL;
	if (OPT_ENABLED(fbskat)) {
		Cp_fs = cExporter::summon("fbskat.res");
		Cp_fs->put("GENE	SZSET	NMARKER	SKATstat	SKATdf	BurdenStat	BurdenDf	P_FBSKAT	P_FBBURDEN\n");
	}
	if (OPT_ENABLED(rvtdt)) {
		Cp_rt = cExporter::summon("rvtdt.res");
		Cp_rt->put("GENE	SZSET	NMARKER	BRV01	BRV05	BRVG01	BRVG05	VTBRV	WSSG\n");
	}

	mGeneDef& Xm_gene = Cp_anaSM->getGeneDef();
	wsUint Q = 0;
	FOREACHDO (mGeneDef_it, Xm_gene, it, Q++) {
		notice("Testing [%d/%d] genes...\r", Q, Xm_gene.size());
		vector<vReal>	__tdtBtable;
		vector<vReal>	__tdtCtable;
		vInt&		Xv_set = it->second;
		m_preprocess(Xv_set);
		/* --genesize option */
		if (IS_ASSIGNED(genesize) && !isInRange(OPT_RANGE(genesize), __getInformVarNum()))
			continue;
		map<string, double> pvals;
		if (!__load(pvals, Xv_set))
			continue;
		m_tdtTableLoad(__tdtBtable, __tdtCtable);

		/* Print out results */
		if (OPT_ENABLED(fbskat)) {
			Cp_fs->fmt("%s	%d	%d	", it->first.c_str(), Xv_set.size(),
				__varToBeAnalyzed.size());
			wsReal R_pSKAT		= skatPval(pvals["FBSKAT1"], pvals["FBSKAT2"]);
			wsReal R_pBurden	= skatPval(pvals["FBSKAT3"], pvals["FBSKAT4"]);
			wsStrCst Sa_cols[] = { "FBSKAT1", "FBSKAT2", "FBSKAT3", "FBSKAT4" };

			wsUint N_test = len(Sa_cols, char*);
			LOOP(i, N_test) {
				map<string, double>::iterator X_find = pvals.find(Sa_cols[i]);
				if (X_find == pvals.end() || NA(X_find->second) || X_find->second == WISARD_INF)
					Cp_fs->put("	NA");
				else
					Cp_fs->fmt("	%g", X_find->second);
			}
			if (NA(R_pSKAT)) Cp_fs->put("	NA"); else Cp_fs->fmt("	%g", R_pSKAT);
			if (NA(R_pBurden)) Cp_fs->put("	NA"); else Cp_fs->fmt("	%g", R_pBurden);
			Cp_fs->put("\n");
		}
		if (OPT_ENABLED(rvtdt)) {
			Cp_rt->fmt("%s	%d	%d	", it->first.c_str(), Xv_set.size(),
				__varToBeAnalyzed.size());
			vStr tests;
			tdtTest("All", 0, 1, OPT_NUMBER(nperm)*2, OPT_NUMBER(nperm), 0.05, pvals,
				tests, __tdtBtable, __tdtCtable, Xv_set);

			wsStrCst Sa_cols[] = { "BRV-T01", "BRV-T05", "BRV-Geno-T01", "BRV-Geno-T05",
				"VT-BRV-Geno", "WSS-Geno" };
			wsUint N_test = len(Sa_cols, char*);
			LOOP(i, N_test) {
				map<string, double>::iterator X_find = pvals.find(Sa_cols[i]);
				if (X_find == pvals.end() || NA(X_find->second) || X_find->second == WISARD_INF)
					Cp_rt->put("	NA");
				else
					Cp_rt->fmt("	%g", X_find->second);
			}
			Cp_rt->put("\n");
		}
	}
	LOG("Testing [%d/%d] genes...\n", Xm_gene.size(), Xm_gene.size());
	/* Close exporter */
	if (OPT_ENABLED(fbskat)) delete Cp_fs;
	if (OPT_ENABLED(rvtdt)) delete Cp_rt;
}

wsUint cRvtdtAnalysis::__getInformVarNum()
{
	wsUint varNum(0);

	FOREACH (vBool_it, __varToBeAnalyzed, yes) {
		if (*yes) varNum++;
	}
	return varNum;
}

#define ro 1
inline wsReal get_skat_stat_fast(wsMat geno, wsReal *Y, wsReal *weight,
	wsUint nbrmarkers, wsReal offset, wsUint N_trio)
{
	wsReal Qstat=0, s;
	if (ro==1) {
		LOOP (i, N_trio) {
			s=0;
			LOOP (j, nbrmarkers)
				s += weight[j]*geno[i][j];
			Qstat += (Y[i]-offset)*s;
		}
		return (Qstat*Qstat);
	}
	if (ro==0)
	{
		LOOP (j, nbrmarkers) {
			s=0;
			LOOP (i, N_trio)
				s += (Y[i]-offset)*geno[i][j];
			Qstat += SQR(weight[j])*SQR(s);
		}
		return (Qstat);
	}
}

char cRvtdtAnalysis::__load(map<string, double> &pvals, vInt& Xv_set)
{
	wsUint			N_var	= (wsUint)Xv_set.size();
	mNucFam&		Xm_nf	= Cp_anaFS->getNucFamData();
	vector<wsMat>	Xv_mats;
	wsUint			N_trio	= 0;
	__trios.clear();
	FOREACH (mNucFam_it, Xm_nf, i) {
		xNucFamily &X_fam = i->second;
		/* For each child, we can construct trio */
		FOREACH (vSampPtr_it, X_fam.Xp_childs, j) {
// 			if (__phased){
// 				phasedTrioPtr oneTrio = phasedTrioPtr(new phasedTrio);
// 				oneTrio->load(getIO(), Xv_set, __varToBeAnalyzed, __skipMiss);
// 				if (oneTrio->checkInformTrio()){
// 					trioPtr baseTrio = oneTrio;
// 					__trios.push_back(baseTrio);
// 				}   
// 				else __nInformTrios--;
// 			} else {
				unPhasedTrioPtr oneTrio = unPhasedTrioPtr(new unPhasedTrio);
				oneTrio->load(getIO(), Xv_set,
					__varToBeAnalyzed, __skipMiss,
					X_fam.Xp_pat->N_idx,
					X_fam.Xp_mat->N_idx, (*j)->N_idx);
				Xv_mats.push_back(oneTrio->getGeno());
				if (oneTrio->checkInformTrio()){
					trioPtr baseTrio = oneTrio;
					__trios.push_back(baseTrio);
				}
				else __nInformTrios--;
//			}
				N_trio++;
		}
	}

	if (OPT_ENABLED(fbskat)) {
		wsVec	Y		= sseVector(N_trio);
		wsVec	zz		= sseEmptyVec(N_var);
		wsVec	pf		= sseVector(N_var);
		wsVec	pf_temp	= sseEmptyVec(N_var);
		wsUint	N_perm	= OPT_NUMBER(nperm);
		wsMat	Ra_tmpG	= sseMatrix(N_trio, N_var);

		wsUint k = 0;
		FOREACH (mNucFam_it, Xm_nf, i) {
			xNucFamily &X_fam = i->second;
			/* For each child, we can construct trio */
			FOREACHDO (vSampPtr_it, X_fam.Xp_childs, j, k++) {
				wsVec p1_genotype = Xv_mats[k][0];
				wsVec p2_genotype = Xv_mats[k][1];
				wsVec x_genotype = Xv_mats[k][2];

				Y[k] = getIO()->getPheno()[(*j)->N_idx];

				for (wsUint J=0 ; J<N_var ; J++) {
					Ra_tmpG[k][J]=x_genotype[J]- Echp1p2[(wsUint)p1_genotype[J]][(wsUint)p2_genotype[J]];
					pf_temp[J] += p1_genotype[J]+p2_genotype[J];
					zz[J] += x_genotype[J];
				}
			}
		}

		for (wsUint J=0 ; J<N_var ; J++)
			pf_temp[J] = (pf_temp[J]+1.0)/(4.0*k+1.0);

		wsMat geno= sseMatrix(N_trio, N_var);
//		MAT_t geno_p = sseMatrix(N_trio, N_var);

		wsUint curr_mark=0;
//		float stat2=0;
		//std::cout << start << " " << nbrmarkers_temp << std::endl;
		LOOP (j, N_var) {
			/* std::cout << zz[j-start+1]+pf_temp[j-start+1] << " " << pf_temp[j-start+1] << " " << mend[j] << " " << pass[j] << " " << weight_ns[j]; */
			// weight_ns[j]>0 OR weight_ns[j]>=0 ???
			if ((zz[j]+pf_temp[j]>0) && (pf_temp[j]<=1)) {
				//if ((zz[j-start+1]+pf_temp[j-start+1]>0) && (pf_temp[j-start+1]<=MAF) && (mend[j]<=mend_th) && (pass[j]>0) && (weight_ns[j]>=0)) {
				/* std::cout << "\tanalyzed"; */
				pf[curr_mark]=pf_temp[j];
				LOOP (i, N_trio)
					geno[i][curr_mark]=Ra_tmpG[i][j];
				curr_mark=curr_mark+1;
			}
			else {
				//std::cout << zz[j-start+1]+pf_temp[j-start+1] << " " << pf_temp[j-start+1] << " " << mend[j] << " " << pass[j] << " " << weight_ns[j] << std::endl;
			}
			/* std::cout << std::endl; */
		}
		wsUint nbrmarkers=curr_mark;
		wsUint min_mark = 1;
		//std::cout << "nbrmarkers: " << nbrmarkers << std::endl;
		if (nbrmarkers >= min_mark) {
			//can also have a weight that is data-dependent
			wsVec weight = sseVector(nbrmarkers);

			if (!IS_ASSIGNED(betaweight)) {
				OPTION().assign("betaweight", OPTION().getDefVal("betaweight"), 1);
				OPTION().FORCE_OPT_STRING(betaweight);
			}
			wsUint N_beta = 0;
			wsVec	Ra_bw		= loadRealValues(OPT_STRING(betaweight), &N_beta);
			if (N_beta != 2) halt("--betaweight should take 2 real values!");

			LOOP (k, nbrmarkers)
				weight[k] = dbeta(pf[k], Ra_bw[0], Ra_bw[1], 0);
				//weight[k]=gsl_ran_beta_pdf(pf[k+1],1,25);
			sseFree(Ra_bw);

			/* for (i=0; i<nbrmarkers; i++) */
			/* 	std::cout << weight[i] << " "; */
			/* std::cout << std::endl; */

			wsReal Ymean = sseVsum(Y, N_trio);
			//if all affected
			if ((Ymean>=N_trio) && (Ymean<=N_trio))  Ymean=0.05;
			else Ymean=Ymean/(wsReal)N_trio;

			/* std::cout << geno[1][1] << " " << Ymean << std::endl; */

			wsReal Qstat=get_skat_stat_fast(geno, Y, weight, nbrmarkers, Ymean, N_trio);

			/* std::cout << Qstat << std::endl; */
			wsVec	Qstat_p = sseVector(OPT_NUMBER(nperm));
			wsVec	permut	= sseVector(N_trio);
			wsMat	geno_p	= sseMatrix(N_trio, nbrmarkers);
			LOOP (b1, N_perm) {
				LOOP (i, N_trio)
					permut[i] = (wsReal)(wsRand()%2);

				LOOP (i, N_trio) LOOP (k, nbrmarkers)
					geno_p[i][k] = permut[i]*geno[i][k];
				Qstat_p[b1]=get_skat_stat_fast(geno_p, Y, weight,nbrmarkers,Ymean, N_trio);
			}

			wsReal mu_Q = sseVsum(Qstat_p, N_perm) / (wsReal)N_perm;
			wsReal var_Q=0;
			LOOP (b, N_perm)
				var_Q += (Qstat_p[b]-mu_Q)*(Qstat_p[b]-mu_Q);
			var_Q=var_Q/(wsReal)N_perm;

			wsReal mu4_Q=0;
			LOOP (b, N_perm)
				mu4_Q=mu4_Q+(Qstat_p[b]-mu_Q)*(Qstat_p[b]-mu_Q)*(Qstat_p[b]-mu_Q)*(Qstat_p[b]-mu_Q);
			mu4_Q=mu4_Q/(wsReal)N_perm;
			wsReal gamma=mu4_Q/(var_Q*var_Q)-3;
			wsReal newdf=12.0/gamma;

			/* std::cout << "===\n" << gene[ig] << " " << start_gene[ig] << " " << end_gene[ig] << " " << Qstat*2*mu_Q/var_Q << " " << 2*mu_Q*mu_Q/var_Q << " " << (Qstat-mu_Q)*sqrt(2*newdf)/sqrt(var_Q)+newdf << " " << newdf << std::endl; */
			if (fabs(var_Q) < 1.0e-6) {
				pvals["FBSKAT1"] = WISARD_NAN;
				pvals["FBSKAT2"] = WISARD_NAN;
				pvals["FBSKAT3"] = WISARD_NAN;
				pvals["FBSKAT4"] = WISARD_NAN;
			} else {
				pvals["FBSKAT1"] = Qstat*2*mu_Q/var_Q;
				pvals["FBSKAT2"] = 2*mu_Q*mu_Q/var_Q;
				pvals["FBSKAT3"] = (Qstat-mu_Q)*sqrt(2*newdf)/sqrt(var_Q)+newdf;
				pvals["FBSKAT4"] = newdf;
			}

			sseFree(Qstat_p);
			sseFree(permut);
			sseFree(weight);
			sseUnmat(geno_p, N_trio);
		}

		sseUnmat(Ra_tmpG, N_trio);
		sseUnmat(geno, N_trio);
		sseFree(zz);
		sseFree(pf);
		sseFree(pf_temp);
		sseFree(Y);
	}

	if (__trios.size()==0)
		return false;
	else
		return true;
}

void cRvtdtAnalysis::m_tdtTableLoad(vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable){
	__denovoCount.resize(0); __effTriosBySite.resize(0); 
	vector<vReal> unTrans;
	for (trioIterator iter=__trios.begin(); iter!=__trios.end(); iter++){
		if (!(*iter)->checkInformTrio()) continue;
		vector<vReal> tdtTable = (*iter)->getTdtTable();
		for (wsUint i=0; i<tdtTable.size(); i++){
			if (i%2==0) __tdtBtable.push_back(tdtTable[i]);
			else __tdtCtable.push_back(tdtTable[i]);
		}
		vector<vReal> unTransObj = (*iter)->getNuTransPatChrs();
		for (wsUint i=0; i<unTransObj.size(); i++) unTrans.push_back(unTransObj[i]);
		m_ifTrueThenAddOne((*iter)->getDenovo(), __denovoCount);
		m_ifTrueThenAddOne((*iter)->getAnalyzed(), __effTriosBySite);
	}
	__wssWeight = m_wssWeight(unTrans); // need __effTriosBySite
	//for (UINT i=0; i<__tdtBtable.size(); i++) std::cout << __tdtBtable[i] << std::endl;
	// write log Json
// 	Json::Value tdtStat;
// 	tdtStat["effectiveTrioNum"] = __nInformTrios;
// 	tdtStat["siteBeAnalyzed"] = buildJsonArray(__varToBeAnalyzed);
// 	tdtStat["Transmitted"] = buildJsonArray(m_colTotal(__tdtCtable));
// 	tdtStat["unTransmitted"] = buildJsonArray(m_colTotal(__tdtBtable));
// 	tdtStat["denovoCount"] = buildJsonArray(__denovoCount);
// 	tdtStat["effectiveTrioNum"] = buildJsonArray(__effTriosBySite);
// 	tdtStat["wssWeight"] = buildJsonArray(__wssWeight);
// 	tdtStat["singleSitePval"] = buildJsonArray(m_singleSiteP());
// 	logRoot["tdtStatic"] = tdtStat;
	return;
}

vReal cRvtdtAnalysis::m_singleSiteP(const vector<vReal> &__tdtBtable, const vector<vReal> &__tdtCtable) const {
	vReal singleSiteP;
	vReal Bcount(m_colTotal(__tdtBtable)), Ccount(m_colTotal(__tdtCtable));

	for (wsUint i=0; i<Bcount.size(); i++){
		if ((Ccount[i] + Bcount[i]) == 0) singleSiteP.push_back(WISARD_INF);
		else singleSiteP.push_back(W1 - pnorm((Ccount[i]-Bcount[i])/sqrt(Ccount[i]+Bcount[i])));
			//gsl_cdf_gaussian_Q((Ccount[i]-Bcount[i])/sqrt(Ccount[i]+Bcount[i]), 1));
	}
	return singleSiteP;

}

vReal cRvtdtAnalysis::m_wssWeight(const vector<vReal>& unTrans)
{
	wsUint N_sz = 0;
	LOOP (i, unTrans.size())
		if (N_sz < (wsUint)unTrans[i].size())
			N_sz = (wsUint)unTrans[i].size();
	vector<wsUint> NotZero(N_sz, 0), NotMiss(N_sz, 0);

	for (wsUint i = 0; i < unTrans.size(); i++){
		for (wsUint j = 0; j < unTrans[i].size(); j++){
			if (j>=NotMiss.size() || j>=NotZero.size())
				halt("Range error");
			if (unTrans[i][j] >= 0.0) {
				NotMiss[j]++;
				NotZero[j] += (wsUint)unTrans[i][j];
			}
		}
	}
	//std::cout << NotMiss << std::endl << NotZero << std::endl;
	vReal weight(N_sz, 0);
	// w=sqrt(n*q*(1-q))
	LOOP (i, N_sz) {
		if (NotMiss[i] == 0) weight[i] = W0;
		else {
			double q = (NotZero[i]+1)/float(NotMiss[i]+1);
			if (q==1) q = (NotZero[i]+1)/float(NotMiss[i]+2);
			else if (q>1) {weight[i] = W0; continue;}
			else {;}
			if (i >= __effTriosBySite.size())
				halt("Range error 2");
			weight[i] = 1/sqrt(__effTriosBySite[i]*4*q*(1-q)); // only consider parents
		}
	}
	return weight;
}

/* set up for rvTDT test */
// shuffle trios, get tdt table and untransmitted haplotypes 
rvTDTshuffleOutput cRvtdtAnalysis::m_shuffle(std::string shuffleMethod, gsl_rng* gslr)
{
	vector<vReal> unTrans, tdtBtable, tdtCtable;
	for (trioIterator iter=__trios.begin(); iter!=__trios.end(); iter++){
		if (!(*iter)->checkInformTrio()) continue;
		shuffleOutput shuffled = (*iter)->shuffle(shuffleMethod,gslr);

		if (shuffled.tdtTable.size()%2 != 0)
			halt("Error: shuffle return a invalid structure");
		for (wsUint i=0; i<shuffled.tdtTable.size(); i++){
			if (i%2 == 0) tdtBtable.push_back(shuffled.tdtTable[i]);
			else tdtCtable.push_back(shuffled.tdtTable[i]);
		}
		unTrans.push_back(shuffled.untransmitted[0]);
		unTrans.push_back(shuffled.untransmitted[1]);
	}
	return rvTDTshuffleOutput(unTrans, tdtBtable, tdtCtable);
}

// aggreate events count
// TODO: how to deal with missing infer in CMC 
tdtResult cRvtdtAnalysis::m_rvAggreate(const vector<vReal>& tdtBtable, const vector<vReal>& tdtCtable, const vBool& tobeAnalyzed, const vReal& weight, bool allowDenovo) const
{
	//for (UINT i=0; i<__tdtBtable.size(); i++) std::cout << __tdtBtable[i] << std::endl;
	//for (UINT i=0; i<__tdtCtable.size(); i++) std::cout << __tdtCtable[i] << std::endl;
	//std::cout << tobeAnalyzed << " | " << weight << " | " << allowDenovo << std::endl;
	tdtResult errorResult(0, 0, WISARD_INF, WISARD_INF);
	if (tdtBtable.size() != tdtCtable.size() ||
		tdtBtable[0].size() != tobeAnalyzed.size() ||
		tdtBtable[0].size() != weight.size()) return errorResult;

	double MZ_B(0), MZ_C(0), CMC_B(0), CMC_C(0);
	for (wsUint i=0; i<tdtBtable.size(); i++){
		if (tdtBtable[i].size() != tdtCtable[i].size()) return errorResult;
		double oneB(0), oneC(0), floorB(0), floorC(0);
		for (wsUint j=0; j<tdtBtable[i].size(); j++){
			if (tobeAnalyzed[j]) {
				oneB += tdtBtable[i][j]*weight[j]; floorB += floor(tdtBtable[i][j])*weight[j];
				// denovo marked as 101
				if (tdtCtable[i][j]>=100){
					double offset=(allowDenovo)?100:101;
					oneC += (tdtCtable[i][j]-offset)*weight[j]; 
					floorC += floor(tdtCtable[i][j]-offset)*weight[j];
				}else {
					oneC += tdtCtable[i][j]*weight[j];
					floorC += floor(tdtCtable[i][j])*weight[j];
				}			}
		}
		MZ_B += oneB; MZ_C += oneC;
		if (__skipMiss && (oneB + oneC !=0)) {CMC_B += oneB/(oneB+oneC); CMC_C += oneC/(oneB+oneC);} // skip miss
		if ((!__skipMiss) && (floorB + floorC !=0)) {CMC_B += floorB/(floorB+floorC); CMC_C += floorC/(floorB+floorC);} // infer miss
	}
	double MZ_chi2(0), CMC_chi2(0);
	double MZ_P(WISARD_INF), CMC_P(WISARD_INF);
	if (MZ_B+MZ_C !=0 ){
		MZ_chi2 = (MZ_C-MZ_B)/sqrt(MZ_B+MZ_C);
		MZ_P = W1 - pnorm((MZ_C-MZ_B)/sqrt(MZ_B+MZ_C));
			//gsl_cdf_gaussian_Q((MZ_C-MZ_B)/sqrt(MZ_B+MZ_C), 1);
	}
	if (CMC_B+CMC_C !=0 ){
		CMC_chi2 = (CMC_C-CMC_B)/sqrt(CMC_B+CMC_C);
		CMC_P = W1 - pnorm((CMC_C-CMC_B)/sqrt(CMC_B+CMC_C));
			//gsl_cdf_gaussian_Q((CMC_C-CMC_B)/sqrt(CMC_B+CMC_C), 1);
	}
	//std::cout << MZ_chi2 <<" " <<  CMC_chi2 <<" " <<  MZ_P <<" " <<  CMC_P << std::endl;
	return tdtResult(MZ_chi2, CMC_chi2, MZ_P, CMC_P);
}

// variable threshold
tdtResult cRvtdtAnalysis::m_VT(const vector<vReal>& tdtBtable, const vector<vReal>& tdtCtable, const vBool& tobeAnalyzed, const vReal& weight, vInt& Xv_set) const
{
	xMaf*		Xp_maf	= getIO()->getMAF();

	tdtResult errorResult(0, 0, WISARD_INF, WISARD_INF);
	if (tdtBtable.size() != tdtCtable.size() ||  tdtBtable[0].size() != tobeAnalyzed.size() || tdtBtable[0].size() != weight.size()) return errorResult;

	vReal validMafs;
	for (wsUint i=0; i<tobeAnalyzed.size(); i++){
		if (tobeAnalyzed[i])
			validMafs.push_back(Xp_maf[Xv_set[i]].R_allMaf);
	}
	std::sort(validMafs.begin(), validMafs.end());
	validMafs.erase(std::unique(validMafs.begin(), validMafs.end()), validMafs.end());

	double maxMZ(-999.0), maxCMC(-999.0);
	for (wsUint i=0; i < validMafs.size(); i++){
		vBool oneAnalyze(tobeAnalyzed);
		for (wsUint j=0; j<tobeAnalyzed.size(); j++){
			if (Xp_maf[Xv_set[j]].R_allMaf > validMafs[i]) oneAnalyze[j] = false;}
		tdtResult res = m_rvAggreate(tdtBtable, tdtCtable, oneAnalyze, weight, false);
		if (res.mz_chi2 > maxMZ) maxMZ = res.mz_chi2;
		if (res.cmc_chi2 > maxCMC) maxCMC = res.cmc_chi2;
	}
	return tdtResult(maxMZ, maxCMC, WISARD_INF, WISARD_INF); // will never use MZ_P, CMC_P
}

/* rvTDT test */ 
// return pvalue for MZ-T1 CMC-T1 MZ-T5 CMC-T5
vReal cRvtdtAnalysis::m_noPermut(const vBool& tobeAnalyzed,
	vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable, vInt& Xv_set)
{
	xMaf*		Xp_maf	= getIO()->getMAF();
	vReal		Rv_fakeWgt(tobeAnalyzed.size(), 1);
	vBool		BL01(__varToBeAnalyzed), BL05(__varToBeAnalyzed);
	
	wsUint i = 0;
	FOREACHDO (vInt_it, Xv_set, I, i++) {
		wsReal R_maf = Xp_maf[*I].R_allMaf;
		if (R_maf > REAL_CONST(0.05)) BL05[i] = false;
		if (R_maf > REAL_CONST(0.01)) BL01[i] = false;
	}
// 	for (wsUint i=0 ; i < __popMafs.size(); i++) {
// 		if (__popMafs[i] > 0.05) BL05[i] = false;
// 		if (__popMafs[i] > 0.01) BL01[i] = false;
// 	}

	tdtResult oriT1 = m_rvAggreate(__tdtBtable, __tdtCtable, BL01, Rv_fakeWgt, true); // MZ CMC 01
	tdtResult oriT5 = m_rvAggreate(__tdtBtable, __tdtCtable, BL05, Rv_fakeWgt, true); // MZ CMC 05

	if (oriT1.mz_pval == WISARD_INF ||
		oriT5.mz_pval == WISARD_INF){
		vReal res(4, WISARD_INF);
		return res;
	}

	vReal nopermut;
	nopermut.push_back(oriT1.mz_pval); //MZ01
	nopermut.push_back(oriT1.cmc_pval); //CMC01
	nopermut.push_back(oriT5.mz_pval); //MZ05
	nopermut.push_back(oriT5.cmc_pval); //CMC05
	return nopermut;
}

/* rvTDT test */ 
// return pvalue for MZ-T1-xxx CMC-T1-xxx MZ-T5-xxx CMC-T5-xxx VT-MZ-xxx, VT-CMC-xxx* and WSS-xxx
vReal cRvtdtAnalysis::m_Permut(std::string shuffleMethod,
	const vBool& tobeAnalyzed, wsUint adaptive, wsUint N_perm,
	double R_alpha, vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable,
	gsl_rng* gslr, vInt& Xv_set)
{
	xMaf*		Xp_maf	= getIO()->getMAF();
	vReal		Rv_fakeWgt(tobeAnalyzed.size(), 1);
	vBool		BL01(__varToBeAnalyzed), BL05(__varToBeAnalyzed);

	wsUint i = 0;
	FOREACHDO (vInt_it, Xv_set, I, i++) {
		wsReal R_maf = Xp_maf[*I].R_allMaf;
		if (R_maf > REAL_CONST(0.50)) BL05[i] = 1;
		if (R_maf > REAL_CONST(0.10)) BL01[i] = 1;
	}
// 	for (wsUint i = 0; i < __popMafs.size(); i++) {
// 		if (__popMafs[i] > 0.05) BL05[i] = false;
// 		if (__popMafs[i] > 0.01) BL01[i] = false;
// 	}
	tdtResult oriT1 = m_rvAggreate(__tdtBtable, __tdtCtable, BL01, Rv_fakeWgt, false); // MZ CMC 01
	tdtResult oriT5 = m_rvAggreate(__tdtBtable, __tdtCtable, BL05, Rv_fakeWgt, false); // MZ CMC 05
	tdtResult oriVT = m_VT(__tdtBtable, __tdtCtable, tobeAnalyzed, Rv_fakeWgt, Xv_set); // VT-MZ VT-CMC
	tdtResult oriWSS = m_rvAggreate(__tdtBtable, __tdtCtable, tobeAnalyzed, __wssWeight, false); // WSS
	//for (UINT x=0; x<__tdtBtable.size(); x++) std::cout << "B:\t" << __tdtBtable[x] << std::endl << "C:\t" << __tdtCtable[x] << std::endl;
	//std::cout << "-------------\n";
	if (oriT1.mz_pval == WISARD_INF || oriT5.mz_pval == WISARD_INF || oriWSS.mz_pval == WISARD_INF || oriVT.mz_chi2 == -99 || oriVT.cmc_chi2 == -99){
		vReal res(7, WISARD_INF);
		return res;
	}
	vReal oriRes;
	oriRes.push_back(oriT1.mz_chi2); oriRes.push_back(oriT1.cmc_chi2); oriRes.push_back(oriT5.mz_chi2); oriRes.push_back(oriT5.cmc_chi2); // MZ CMC 
	oriRes.push_back(oriVT.mz_chi2); oriRes.push_back(oriVT.cmc_chi2); oriRes.push_back(oriWSS.mz_chi2); // VT-MZ; VT-CMC; WSS
	//std::cout << oriRes << std::endl;
	if (shuffleMethod == "NoPermut"){
		vReal nopermut;
		nopermut.push_back(oriT1.mz_pval); //MZ01
		nopermut.push_back(oriT1.cmc_pval); //CMC01
		nopermut.push_back(oriT5.mz_pval); //MZ05
		nopermut.push_back(oriT5.cmc_pval); //CMC05
		return nopermut;
	}

	vector<wsUint> permCount(7,0);
	vReal permPval(7,9.0);
	for (wsUint i=1; i<=N_perm; i++) {
		rvTDTshuffleOutput Obj=m_shuffle(shuffleMethod, gslr);
		//for (UINT x=0; x<Obj.tdtBtable.size(); x++) std::cout << "B:\t" << Obj.tdtBtable[x] << std::endl << "C:\t" << Obj.tdtCtable[x] << std::endl;        
		tdtResult permT1 = m_rvAggreate(Obj.tdtBtable, Obj.tdtCtable, BL01, Rv_fakeWgt, false); // MZ CMC 01
		tdtResult permT5 = m_rvAggreate(Obj.tdtBtable, Obj.tdtCtable, BL05, Rv_fakeWgt, false); // MZ CMC 05
		tdtResult permVT = m_VT(Obj.tdtBtable, Obj.tdtCtable, tobeAnalyzed, Rv_fakeWgt, Xv_set); // VT-MZ, VT-CMC
		tdtResult permWSS = m_rvAggreate(Obj.tdtBtable, Obj.tdtCtable, tobeAnalyzed, m_wssWeight(Obj.unTrans), false);

		vReal permRes; 
		permRes.push_back((permT1.mz_pval==WISARD_INF)?0:permT1.mz_chi2); //MZ01
		permRes.push_back((permT1.cmc_pval==WISARD_INF)?0:permT1.cmc_chi2); //CMC01
		permRes.push_back((permT5.mz_pval==WISARD_INF)?0:permT5.mz_chi2); //MZ05
		permRes.push_back((permT5.cmc_pval==WISARD_INF)?0:permT5.cmc_chi2); //CMC05
		permRes.push_back((permVT.mz_chi2==-99)?0:permVT.mz_chi2); // VTMZ
		permRes.push_back((permVT.cmc_chi2==-99)?0:permVT.cmc_chi2); //VTCMC
		permRes.push_back((permWSS.mz_pval==WISARD_INF)?0:permWSS.mz_chi2); //WSS
		//std::cout << i << "  " << permRes << std::endl;
		m_updateCount(permCount, oriRes, permRes, gslr);
		if (i % adaptive != 0 || i == 0) {;}
		else{
			m_updatePvalue(permPval, permCount, i, R_alpha);
			bool MZok(true), CMCok(true);
			for (wsUint t=0; t<permPval.size(); t++){
				if (permPval[t] > 1) {
					if (t%2) CMCok=false; //CMC series: 1 3 5
					else MZok=false;  //MZ series: 0 2 4 6
				}
			}
// 			if (__phased){
// 				if (MZok && CMCok) break; // CMC works for phased only
// 			} else {
				if (MZok) break;
//			}
		}
	}
	for (wsUint i=0; i<permCount.size(); i++){
		permPval[i] = (permPval[i]<=1.0)?permPval[i]:(1.0*permCount[i]+1)/(1.0*N_perm+1);
	}
	return permPval;
}

void cRvtdtAnalysis::tdtTest(std::string GenoHapo, wsReal R_mafL, wsReal R_mafU,
	wsUint adaptive, wsUint N_perm, wsReal R_alpha, std::map<std::string, double>& pvalues,
	std::vector<std::string>& tests, vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable, vInt& Xv_set)
{
	// TODO: check the parameters
//	RNG rng;
	gsl_rng gslr = { 0, 1 };
	pvalues.clear(); tests.clear();
	vReal fakeWeight(Xv_set.size(), 1); // for MZ CMC, no weight

	vBool selectedVar(__varToBeAnalyzed); // first use user-defined criteria 
	xMaf*		Xp_maf	= getIO()->getMAF();
	wsUint i = 0;
	FOREACHDO (vInt_it, Xv_set, I, i++) {
		wsReal R_maf = Xp_maf[*I].R_allMaf;
		if (R_maf<R_mafL || R_maf>R_mafU) selectedVar[i] = false;
	}
// 	for (wsUint i=0; i<__popMafs.size(); i++) {
// 		if (__popMafs[i]<R_mafL || __popMafs[i]>R_mafU) selectedVar[i] = false;
// 	}
	vReal noPermut = m_noPermut(selectedVar, __tdtBtable, __tdtCtable,
		Xv_set); //"NoPermut", return MZ-T1 CMC-T1 MZ-T5 and CMC-T5

	if (GenoHapo == "NoPermut"){ 
		//std::cout << noPermut << std::endl;
		pvalues["BRV-T01"] = noPermut[0]; tests.push_back("BRV-T01");
//		if (__phased) {pvalues["CMC-T01"] = noPermut[1]; tests.push_back("CMC-T01");}
		pvalues["BRV-T05"] = noPermut[2]; tests.push_back("BRV-T05");
//		if (__phased) {pvalues["CMC-T05"] = noPermut[3]; tests.push_back("CMC-T05");}
	}
	else if (GenoHapo == "GenoPermut" || GenoHapo == "HapoPermut"){
		if (GenoHapo == "HapoPermut")// && __phased == false)
			LOGwarn("Warning: using hapotype permutation on unphased data\n");
		std::string shuffleMethod = (GenoHapo == "GenoPermut")?"GenoShuffle":"HapoShuffle";
		vReal permPval = m_Permut(shuffleMethod, selectedVar, adaptive,
			N_perm, R_alpha, __tdtBtable, __tdtCtable, &gslr, Xv_set);

		if (GenoHapo == "GenoPermut") {
			pvalues["BRV-T01"] = noPermut[0]; tests.push_back("BRV-T01");
//			if (__phased) {pvalues["CMC-T01"] = noPermut[1]; tests.push_back("CMC-T01");}
			pvalues["BRV-T05"] = noPermut[2]; tests.push_back("BRV-T05");
//			if (__phased) {pvalues["CMC-T05"] = noPermut[3]; tests.push_back("CMC-T05");}

			pvalues["BRV-Geno-T01"] = permPval[0]; tests.push_back("BRV-Geno-T01");
//			if (__phased) {pvalues["CMC-Geno-T01"] = permPval[1]; tests.push_back("CMC-Geno-T01");}
			pvalues["BRV-Geno-T05"] = permPval[2]; tests.push_back("BRV-Geno-T05");
//			if (__phased) {pvalues["CMC-Geno-T05"] = permPval[3]; tests.push_back("CMC-Geno-T05");}
			pvalues["VT-BRV-Geno"] = permPval[4]; tests.push_back("VT-BRV-Geno");
//			if (__phased) {pvalues["VT-CMC-Geno"] = permPval[5];  tests.push_back("VT-CMC-Geno");}
			pvalues["WSS-Geno"] = permPval[6];  tests.push_back("WSS-Geno");
		} else {
			pvalues["BRV-T01"] = noPermut[0]; tests.push_back("BRV-T01");
//			if (__phased) {pvalues["CMC-T01"] = noPermut[1]; tests.push_back("CMC-T01");}
			pvalues["BRV-T05"] = noPermut[2]; tests.push_back("BRV-T05");
//			if (__phased) {pvalues["CMC-T05"] = noPermut[3]; tests.push_back("CMC-T05");}

			pvalues["BRV-Hapo-T01"] = permPval[0]; tests.push_back("BRV-Hapo-T01");
//			if (__phased) {pvalues["CMC-Hapo-T01"] = permPval[1]; tests.push_back("CMC-Hapo-T01");}
			pvalues["BRV-Hapo-T05"] = permPval[2]; tests.push_back("BRV-Hapo-T05");
//			if (__phased) {pvalues["CMC-Hapo-T05"] = permPval[3]; tests.push_back("CMC-Hapo-T05");}
			pvalues["VT-BRV-Hapo"] = permPval[4]; tests.push_back("VT-BRV-Hapo");
//			if (__phased) {pvalues["VT-CMC-Hapo"] = permPval[5];  tests.push_back("VT-CMC-Hapo");}
			pvalues["WSS-Hapo"] = permPval[6];  tests.push_back("WSS-Hapo");
		}
	}
	else if (GenoHapo == "All"){
		//vectorF genoPval = m_Permut("GenoShuffle", selectedVar, adaptive, PermutateTimes, alpha, gslr);

// 		if (__phased) {
// 			vectorF hapoPval = m_Permut("HapoShuffle", selectedVar, adaptive,
// 				N_perm, R_alpha, &gslr, Xv_set);
// 
// 			pvalues["BRV-T01"] = noPermut[0]; tests.push_back("BRV-T01");
// 			pvalues["CMC-T01"] = noPermut[1]; tests.push_back("CMC-T01");
// 			pvalues["BRV-T05"] = noPermut[2]; tests.push_back("BRV-T05");
// 			pvalues["CMC-T05"] = noPermut[3]; tests.push_back("CMC-T05");
// 
// 			//pvalues["BRV-Geno-T01"] = genoPval[0]; tests.push_back("BRV-Geno-T01");
// 			//pvalues["CMC-Geno-T01"] = genoPval[1]; tests.push_back("CMC-Geno-T01");
// 			//pvalues["BRV-Geno-T05"] = genoPval[2]; tests.push_back("BRV-Geno-T05");
// 			//pvalues["CMC-Geno-T05"] = genoPval[3]; tests.push_back("CMC-Geno-T05");
// 			//pvalues["VT-BRV-Geno"] = genoPval[4]; tests.push_back("VT-BRV-Geno");
// 			//pvalues["VT-CMC-Geno"] = genoPval[5];  tests.push_back("VT-CMC-Geno");
// 			//pvalues["WSS-Geno"] = genoPval[6];  tests.push_back("WSS-Geno");
// 
// 			pvalues["BRV-Hapo-T01"] = hapoPval[0]; tests.push_back("BRV-Hapo-T01");
// 			pvalues["CMC-Hapo-T01"] = hapoPval[1]; tests.push_back("CMC-Hapo-T01");
// 			pvalues["BRV-Hapo-T05"] = hapoPval[2]; tests.push_back("BRV-Hapo-T05");
// 			pvalues["CMC-Hapo-T05"] = hapoPval[3]; tests.push_back("CMC-Hapo-T05");
// 			pvalues["VT-BRV-Hapo"] = hapoPval[4]; tests.push_back("VT-BRV-Hapo");
// 			pvalues["VT-CMC-Hapo"] = hapoPval[5];  tests.push_back("VT-CMC-Hapo");
// 			pvalues["WSS-Hapo"] = hapoPval[6];  tests.push_back("WSS-Hapo");
// 		} else {
			vReal genoPval = m_Permut("GenoShuffle", selectedVar, adaptive,
				N_perm, R_alpha, __tdtBtable, __tdtCtable, &gslr, Xv_set);
			pvalues["BRV-T01"] = noPermut[0]; tests.push_back("BRV-T01");
			pvalues["BRV-T05"] = noPermut[2]; tests.push_back("BRV-T05");
			pvalues["BRV-Geno-T01"] = genoPval[0]; tests.push_back("BRV-Geno-T01");
			pvalues["BRV-Geno-T05"] = genoPval[2]; tests.push_back("BRV-Geno-T05");
			pvalues["VT-BRV-Geno"] = genoPval[4]; tests.push_back("VT-BRV-Geno");
			pvalues["WSS-Geno"] = genoPval[6];  tests.push_back("WSS-Geno");
//		}
	}
	else {;}
	return;
}

//// for debug: try to see difference between pvals from old and new methods
//vectorF rareTdt::m_CmpPermut(std::string GenoHapo, const vectorL& tobeAnalyzed, UINT adaptive, UINT PermutateTimes, double alpha, gsl_rng* gslr, double lowerMaf, double upperMaf)
//{
//vectorF fakeWeight;
//for (UINT i = 0; i < tobeAnalyzed.size(); i++) fakeWeight.push_back(1);

//vectorF oriRes;
////for (UINT i=0; i<__tdtBtable.size(); i++) std::cout << __tdtBtable[i] << std::endl;
//tdtResult res = m_VT(__tdtBtable, __tdtCtable, tobeAnalyzed, fakeWeight); // MZ CMC

////std::cout << "ori: " << res2.mz_chi2 << "\t" << res.mz_chi2 << "\t" << res.cmc_chi2 << std::endl;
//if (res.mz_pval == WISARD_INF){
//vectorF res(1, WISARD_INF);
//return res;
//}
//vectorF samMafs;
//for (UINT i=0; i<__samMafs.size(); i++) {
//if (__samMafs[i]>0.5) samMafs.push_back(1-__samMafs[i]);
//else samMafs.push_back(__samMafs[i]);
//}
//std::string NotX = "NotX";      
//unphasedTDT tdtTest(__genotypes, __popMafs, samMafs, NotX, __missCutoff*3, __skipMiss);
//vectorF oriVT = tdtTest.m_tdtVT(0, 1);
//double oldVTMZ = oriVT[0];
//oriRes.push_back(res.mz_chi2); oriRes.push_back(oldVTMZ);
////std::cout << res.mz_chi2 << "\t" << oldVTMZ << std::endl;
////exit(1);
//vectorUI permCount(2,0);
//vectorF permPval(2,9.0);
//UINT __nPermut(0);
//for (UINT i=1; i <= PermutateTimes; i++)
//{
//vector2F tdtBtable, tdtCtable;
//vector2F genos;
//for (trioIterator iter=__trios.begin(); iter!=__trios.end(); iter++){
//if (!(*iter)->checkInformTrio()) continue;
//vector2F geno;
//shuffleOutput shuffled = (*iter)->shuffleWithGenos(GenoHapo,gslr,geno);
//for (UINT j=0; j<geno.size(); j++) genos.push_back(geno[j]);

//if (shuffled.tdtTable.size()%2 != 0 ) {std::cerr << "Error: shuffle return a invalid structure" << std::endl; exit(-1);}
//for (UINT i=0; i<shuffled.tdtTable.size(); i++){
//if (i%2 == 0) tdtBtable.push_back(shuffled.tdtTable[i]);
//else tdtCtable.push_back(shuffled.tdtTable[i]);
//}
//}
//// may be need to change a little bit
//tdtResult perm = m_VT(tdtBtable, tdtCtable, tobeAnalyzed, fakeWeight);
//// oldmethod: 
//std::string NotX = "NotX";    
//unphasedTDT tdtTest(genos, __popMafs, samMafs, NotX, __missCutoff*3, __skipMiss);
//vectorF permVT = tdtTest.m_tdtVT(0, 1);
//double pOldVTMZ = permVT[0];

////if (perm.mz_chi2 - pOldVTMZ > 0.0001) {std::cout << perm.mz_chi2 << "\t" << pOldVTMZ << std::endl;}
//vectorF permRes; // 
//permRes.push_back((perm.mz_chi2 == -99)?0:perm.mz_chi2);
//permRes.push_back((oldVTMZ == WISARD_NA)?0:pOldVTMZ);
//__nPermut++;
//m_updateCount(permCount, oriRes, permRes, gslr);
//if (i % adaptive != 0 || i == 0) {;}
//else{
//m_updatePvalue(permPval, permCount, i, alpha);
//if (permPval[0] <= 1.0 && permPval[2] <= 1.0){
//if (__phased){
//if (permPval[1] <= 1.0) break;
//} else break;
//}
//}
//}
//for (UINT i=0; i<permCount.size(); i++){
//permPval[i] = (permPval[i] <= 1.0)? permPval[i] : (1.0 * permCount[i] + 1) / (1.0 * PermutateTimes + 1);
//}
//return permPval;
//}

// for debug: to see if the permuted MZ statistic is normal distribution
vReal cRvtdtAnalysis::m_MzCmcPermut(vBool& tobeAnalyzed, wsUint adaptive, wsUint PermutateTimes, double alpha, gsl_rng* gslr, double mafUpper,
	vector<vReal> &__tdtBtable, vector<vReal> &__tdtCtable, vInt& Xv_set)
{
	xMaf*	Xp_maf	= getIO()->getMAF();
	vReal	Rv_fakeWeight(tobeAnalyzed.size(), 1);
	vBool	Bv_BLcutoff(__varToBeAnalyzed); //, BL05(__varToBeAnalyzed);

	for (wsUint i = 0; i < Xv_set.size(); i++) {
		if (Xp_maf[Xv_set[i]].R_allMaf > mafUpper) Bv_BLcutoff[i] = false;
	}
	vReal oriRes;
	tdtResult res = m_rvAggreate(__tdtBtable, __tdtCtable, Bv_BLcutoff, Rv_fakeWeight, false); // MZ CMC BLcutoff
	//tdtResult res = m_rvAggreate(__tdtBtable, __tdtCtable, BL05, fakeWeight); // MZ BL05
	if (res.mz_pval == WISARD_INF){
		vReal res(6, WISARD_INF);
		return res;
	}
	// original statistic 
	oriRes.push_back(res.mz_chi2), oriRes.push_back(res.cmc_chi2); oriRes.push_back(res.mz_chi2), oriRes.push_back(res.cmc_chi2); // MZ CMC MZ CMC
	vector<wsUint> permCount(4,0);
	vReal permPval(4,9.0);

	for (wsUint i=1; i <= PermutateTimes; i++) {
		vector<vReal> tdtBtableGeno, tdtCtableGeno, tdtBtableHapo, tdtCtableHapo;
		for (trioIterator iter=__trios.begin(); iter!=__trios.end(); iter++){
			if (!(*iter)->checkInformTrio()) continue;
			shuffleOutput shuffled = (*iter)->shuffle("GenoShuffle",gslr);
			for (wsUint i=0; i<shuffled.tdtTable.size(); i++){
				if (i%2 == 0) tdtBtableGeno.push_back(shuffled.tdtTable[i]);
				else tdtCtableGeno.push_back(shuffled.tdtTable[i]);
			}
			shuffleOutput shuffledHapo = (*iter)->shuffle("HapoShuffle",gslr);
			for (wsUint i=0; i<shuffledHapo.tdtTable.size(); i++){
				if (i%2 == 0) tdtBtableHapo.push_back(shuffledHapo.tdtTable[i]);
				else tdtCtableHapo.push_back(shuffledHapo.tdtTable[i]);
			}
		}
		// TODO: may be need to change a little bit
		tdtResult permGeno = m_rvAggreate(tdtBtableGeno, tdtCtableGeno, Bv_BLcutoff, Rv_fakeWeight, false);
		tdtResult permHapo = m_rvAggreate(tdtBtableHapo, tdtCtableHapo, Bv_BLcutoff, Rv_fakeWeight, false);
		//std::cout << permGeno.mz_chi2 << "\t" << permHapo.mz_chi2  << std::endl;
		vReal permRes; // GenoMZ GenoCMC HapoMZ HapoCMC
		permRes.push_back((permGeno.mz_pval == WISARD_INF)?0:permGeno.mz_chi2);
		permRes.push_back((permGeno.cmc_pval == WISARD_INF)?0:permGeno.cmc_chi2);
		permRes.push_back((permHapo.mz_pval == WISARD_INF)?0:permHapo.mz_chi2);
		permRes.push_back((permHapo.cmc_pval == WISARD_INF)?0:permHapo.cmc_chi2);

		m_updateCount(permCount, oriRes, permRes, gslr);
		if (i % adaptive != 0 || i == 0) {;}
		else{
			m_updatePvalue(permPval, permCount, i, alpha);
			if (permPval[0] <= 1.0 && permPval[2] <= 1.0){
// 				if (__phased){
// 					if (permPval[1] <= 1.0 && permPval[3] <= 1.0 ) break;
// 				} else
					break;
			}
		}
	}
	for (wsUint i=0; i<permCount.size(); i++){
		permPval[i] = (permPval[i] <= 1.0)? permPval[i] : (1.0 * permCount[i] + 1) / (1.0 * PermutateTimes + 1);
	}
	permPval.push_back(res.mz_pval); permPval.push_back(res.cmc_pval);
	//std::cout << permPval << std::endl;
	return permPval;
}








/** Base Trio **/
shuffleOutput trio::shuffle(std::string GenoHapo, gsl_rng* gslr){
	shuffleOutput fakeOne;
	return fakeOne;
}

shuffleOutput trio::shuffleWithGenos(std::string GenoHapo, gsl_rng* gslr, vector<vReal>& genos){
	shuffleOutput fakeOne;
	return fakeOne;
}

/** Phased Trio **/
bool ChromEqual(const vReal& ChromOne, const vReal& ChromTwo, bool allowMissing, bool allowDenovo)
{
	//We allow at most 1 de novo mutation
	if(ChromOne.size() != ChromTwo.size()) return false;
	double AllowedErrorNum = 1;
	int ErrorFound(0);
	bool StrictEqual(true), MissEqual(true), DenovoEqual(true);

	for(wsUint i=0; i<ChromOne.size(); i++){
		if(ChromOne[i] != ChromTwo[i]){
			StrictEqual = false;
			if(ChromOne[i] != WISARD_NA && ChromTwo[i] != WISARD_NA) { MissEqual = false; ErrorFound++;}
		}
	}
	if(ErrorFound > AllowedErrorNum) DenovoEqual = false;

	if((!allowMissing) && (!allowDenovo)) {return StrictEqual;} 
	else if(allowMissing && (!allowDenovo)) {return MissEqual;}
	else if(allowMissing && allowDenovo) {return DenovoEqual;}
	else {return false;}
}

/** Can be used by both **/
vector<vReal> trio::m_GenoShuffle(gsl_rng* gslr)
{
	xMaf*		Xp_maf	= Cp_IO->getMAF();

	vector<wsUint> recordMissing(__tobeAnalyzed.size(), 0); // only record informative sites: 0 no missing, 1 father missing, 2 mother missing
	vector<vReal> shuffledGeno(__pesudoGeno);
	if(shuffledGeno.size() != 4)
		halt("Error: pesudo genotype for shuffle is not set up properly.");

	wsUint j = 0;
	FOREACHDO (vInt_it, *Xv_set, it, j++) {
		wsReal R_maf = Xp_maf[*it].R_allMaf;
		if(!__tobeAnalyzed[j]) continue;
		vReal fakeParents;
		for(wsUint i=0; i<4; i++){
			if(shuffledGeno[i][j] < 0){
				if(BINV(gslr, R_maf)) fakeParents.push_back(1);
				else fakeParents.push_back(0);
				recordMissing[j] = (i<=1)?1:2; // record missing
			}
			else fakeParents.push_back(shuffledGeno[i][j]);
		}
		random_shuffle(fakeParents.begin(), fakeParents.begin()+4);
		for(wsUint i=0; i<4; i++) shuffledGeno[i][j] = fakeParents[i];
	}

	if(BIN(gslr)){
		shuffledGeno.push_back(shuffledGeno[wsRand()%2]);
		shuffledGeno.push_back(shuffledGeno[wsRand()%2 + 2]);
	}
	else{
		shuffledGeno.push_back(shuffledGeno[wsRand()%2 + 2]);
		shuffledGeno.push_back(shuffledGeno[wsRand()%2]);
	}

	// recode missing
	for(wsUint i=0; i<recordMissing.size(); i++){
		if (recordMissing[i] == 0) {;}
		else if(recordMissing[i] == 1) {shuffledGeno[0][i] = WISARD_NA; shuffledGeno[1][i] = WISARD_NA;}
		else if(recordMissing[i] == 2) {shuffledGeno[2][i] = WISARD_NA; shuffledGeno[3][i] = WISARD_NA;}
		else {;}    
	}

	if(shuffledGeno.size() != 6)
		halt("Error: genotype shuffle goes wrong.");
	return shuffledGeno;
}

// no hapotype permutation for unphased data
vector<vReal> trio::m_HapoShuffle(gsl_rng* gslr)
{
    vector<wsUint> recordMissing(__tobeAnalyzed.size(), 0); // only record informative sites: 0 no missing, 1 father missing, 2 mother missing
    vector<vReal> shuffledGeno(__pesudoGeno);
    if(shuffledGeno.size() != 4)
		halt("Error: pesudo genotype for shuffle is not set up properly.");
    
    // refill missing site
	xMaf*		Xp_maf	= Cp_IO->getMAF();
    for(wsUint j=0; j<shuffledGeno[0].size(); j++){
        if(!__tobeAnalyzed[j]) continue;
		int N_idx = (*Xv_set)[j];
        
        if(shuffledGeno[0][j] < 0 || shuffledGeno[1][j] < 0 || shuffledGeno[2][j] < 0 || shuffledGeno[3][j] < 0){
            if((shuffledGeno[0][j] < 0 && shuffledGeno[1][j] < 0) || (shuffledGeno[2][j] < 0 && shuffledGeno[3][j] < 0))
				halt("Error: pesudo genotype wrong, no guessed genotype in missing sites.");
            // both parent missing
            if(shuffledGeno[0][j] + shuffledGeno[1][j] + shuffledGeno[2][j] + shuffledGeno[3][j] < WISARD_NA)
				halt("Error: pesudo genotype wrong, both missing are analyzed here");
            
            int base = (shuffledGeno[0][j] + shuffledGeno[1][j] < 0)?0:2; // which parent is missing
            recordMissing[j] = (base == 0)?1:2;
            
            double infered = (shuffledGeno[base+0][j] < 0)?shuffledGeno[base+1][j]:shuffledGeno[base+0][j];
            double guessed = (BINV(gslr, Xp_maf[N_idx].R_allMaf))?1:0;
            if(BIN(gslr)) {shuffledGeno[base+0][j] = infered; shuffledGeno[base+1][j] = guessed;}
            else {shuffledGeno[base+0][j] = guessed; shuffledGeno[base+1][j] = infered;}
        }
    }
    
    random_shuffle(shuffledGeno.begin(), shuffledGeno.begin()+4);

    if(BIN(gslr)){
        shuffledGeno.push_back(shuffledGeno[wsRand()%2]);
        shuffledGeno.push_back(shuffledGeno[wsRand()%2 + 2]);
    }
    else{
        shuffledGeno.push_back(shuffledGeno[wsRand()%2 + 2]);
        shuffledGeno.push_back(shuffledGeno[wsRand()%2]);
    }
    // recode missing
    for(wsUint i=0; i<recordMissing.size(); i++){
        if (recordMissing[i] == 0) {;}
        else if(recordMissing[i] == 1) {shuffledGeno[0][i] = -9; shuffledGeno[1][i] = -9;}
        else if(recordMissing[i] == 2) {shuffledGeno[2][i] = -9; shuffledGeno[3][i] = -9;}
        else {;}    
    }
    
    if(shuffledGeno.size() != 6)
		halt("Error: hapotype shuffle goes wrong.");
    return shuffledGeno;
}


/** Unphased Trio **/
void unPhasedTrio::load(cIO *Cp_IO, vInt& Xv_set, const vBool& tobeAnalyzed,
	bool skipMiss, wsUint N_idxPat, wsUint N_idxMat, wsUint N_idxChi)
{
	__tobeAnalyzed	= tobeAnalyzed;
	__skipMiss = skipMiss;
	__PesudoGenoReady = false;
	__isInformative=true;
	this->Xv_set = &Xv_set;
	this->Cp_IO = Cp_IO;
	N_pat = N_idxPat;
	N_mat = N_idxMat;
	N_chi = N_idxChi;
// 	//data check
// 	if(__genotype.size() != 3 || __genotype[0].size() != __popMafs.size() || __tobeAnalyzed.size() != __popMafs.size()){
// 		__isInformative=false; return;
// 	}
// 	else {
// 		for(UINT i=0; i<__genotype.size(); i++) {
// 			if(__genotype[i].size() != __popMafs.size()) {
// 				__isInformative=false; return;
// 			}
// 		}
// 	}
	tdtTableCount();
	return;
}

void unPhasedTrio::tdtTableCount()
{
	wsUint	N_sz	= (wsUint)Xv_set->size();
	char**	Na_data	= Cp_IO->getGenotype();
	__tdtTable.resize(2); 
	__tdtTable[0].assign(N_sz, 0); __tdtTable[1].assign(N_sz, 0); //__tdtTable[2].assign(__tobeAnalyzed.size(), 0); __tdtTable[3].assign(__tobeAnalyzed.size(), 0); 
	__unTransParentalChrs.resize(2);
	__unTransParentalChrs[0].assign(N_sz, WISARD_NA);
	__unTransParentalChrs[1].assign(N_sz, WISARD_NA);

	if (OPT_ENABLED(fbskat))
		Ra_geno = sseMatrix(3, N_sz);

	// label out uninformative sites
	wsUint i = 0;
	FOREACHDO (vInt_it, *Xv_set, it, i++) {
		wsUint N_gPat = Na_data[N_pat][*it];
		wsUint N_gMat = Na_data[N_mat][*it];
		wsUint N_gChi = Na_data[N_chi][*it];

//	for(wsUint i=0; i<__genotype[0].size(); i++){
		if(isMissing(N_gPat) && isMissing(N_gMat)) __tobeAnalyzed[i]=false; // both parents missing
		else if(isMissing(N_gChi)) __tobeAnalyzed[i]=false; // kid missing
		else if(__skipMiss){
			if(isMissing(N_gPat) & isMissing(N_gMat) & isMissing(N_gChi))
				__tobeAnalyzed[i]=false; // individual missing
		}
		if (OPT_ENABLED(fbskat)) {
			Ra_geno[0][i] = N_gPat;
			Ra_geno[1][i] = N_gMat;
			Ra_geno[2][i] = N_gChi;
		}
// 		else {;}
// 
// 		for(wsUint j =0; j<__genotype.size(); j++){
// 			if(__genotype[j][i] != 0 && __genotype[j][i] != 1 && __genotype[j][i] != 2 && __genotype[j][i] != -9) __tobeAnalyzed[i]=false;
// 		}
	}
	// uninformative sites as missing
//	trimNonSenseSites(__genotype, __tobeAnalyzed);

	__denovoSite.assign(N_sz, 0);
	i = 0;
	FOREACHDO (vInt_it, *Xv_set, it, i++) {
//	for(wsUint i=0; i<__genotype[0].size(); i++){
		if(!__tobeAnalyzed[i]) {continue;}
		wsUint N_gPat = Na_data[N_pat][*it];
		wsUint N_gMat = Na_data[N_mat][*it];
		wsUint N_gChi = Na_data[N_chi][*it];
		m_unPhasedTdtCount(N_gPat, N_gMat, N_gChi, i);
	}
	// in case reverse mutation, double mutations
	bool allSiteWrong(true);
	LOOP (i, N_sz)
		if (__tobeAnalyzed[i]) {
			allSiteWrong=false; break;
		}
	if (allSiteWrong)
		__isInformative=false;
}

// modify __tdtTable, __unTransParentalChrs, __denovoSite and __tobeAnalyzed; through index
// only analyze valid site, that is no double missing, no kid missing, no coding error
void unPhasedTrio::m_unPhasedTdtCount(double fat, double mot, double kid, wsUint index)
{
	wsReal R_maf = Cp_IO->getMAF()[index].R_allMaf;

    if(!__tobeAnalyzed[index]) {return;}
    // no missing
    if(fat >= 0 && mot >= 0){
        if(fat == 0){
            if(mot == 0){
                if(kid == 0) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 101; __denovoSite[index]=true;} // denovo
                else {__tobeAnalyzed[index]=false;} // no double mutation
            }
            else if(mot == 1){
                if(kid == 0) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += 1;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 102; __denovoSite[index]=true;} //transmitted and denovo
                else {;}
            }
            else if(mot == 2){
                if(kid == 0) {__tobeAnalyzed[index]=false;} // no reverse mutation
                else if(kid == 1) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 1; __tdtTable[1][index] += 101;__denovoSite[index]=true;} // denovo 
                else {;}
            }
            else {;}
        }
        else if(fat == 2){
            if(mot == 0){
                if(kid == 0) {__tobeAnalyzed[index]=false;} // no reverse mutation
                else if(kid == 1) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 0;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 101;__denovoSite[index]=true;} // denovo
                else {;}
            }
            else if(mot == 1){
                if(kid == 0) {__tobeAnalyzed[index]=false;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 1;}
                else {;}
            }
            else if(mot == 2){
                if(kid == 2) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 1;}
                else {__tobeAnalyzed[index]=false;}
            }
            else {;}
        }
        else if(fat == 1){
            if(mot == 0){
                if(kid == 0) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 0; __tdtTable[0][index] += 1;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 102; __denovoSite[index]=true;} // denovo
                else {;}
            }
            else if(mot == 1){
                if(kid == 0) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += 2;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 0; __tdtTable[0][index] += 1; __tdtTable[1][index] += 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 2;}
                else {;}
            }
            else if(mot == 2){
                if(kid == 0) {__tobeAnalyzed[index]=false;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 1; __tdtTable[1][index] += 1;}
                else {;}
            }
            else {;}
	}
        else {;}
    }
    // one parent missing: no denovo mutation when missing
    else if(fat == -9 || mot == -9){
        if(fat == -9 && mot == -9) {__tobeAnalyzed[index]=false; return;}
        double missing = (fat==WISARD_NA)?fat:mot;
        mot = (fat==WISARD_NA)?mot:fat; fat = missing; // assume missing father
        
		wsReal	R_mafInv = W1 - R_maf;
        double heteroProb = R_maf*R_mafInv*2;
        double homoProb = SQR(R_maf);
        double wildProb = SQR(R_mafInv);
        
        if(mot == 0){
            if(kid == 0) {__unTransParentalChrs[1][index] = 0; __tdtTable[0][index] += heteroProb/(wildProb+heteroProb);} // -9X0 -> 0 ==> 1X0 -> 0
            else if(kid == 1) {__unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += heteroProb/(homoProb+heteroProb);} // -9X0 -> 1 ==> 1X0 -> 1
            else if(kid == 2) {__tobeAnalyzed[index]=false;}
            else {;}
        }
        else if(mot == 1){
            if(kid == 0) {__unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += (1+heteroProb/(wildProb+heteroProb));} // -9X1 -> 0 ==> 1X1 -> 0 or 0X1 -> 0
            else if(kid == 1) {__unTransParentalChrs[1][index] = 0; __tdtTable[0][index] += homoProb+heteroProb; __tdtTable[1][index] += wildProb+heteroProb;} // -9X1 -> 1 ==> 0X1 -> 1 or 1X1 -> 1 or 2X1 -> 1
            else if(kid == 2) {__unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += (1+heteroProb/(homoProb+heteroProb));} // -9X1 -> 2 ==> 1X1 -> 2 or 2X1 -> 2
            else {;}
        }
        else if(mot == 2){
            if(kid == 0) {__tobeAnalyzed[index]=false;} 
            else if(kid == 1) {__unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += heteroProb/(wildProb+heteroProb);} // -9X2 -> 1 ==> 1X2 -> 1
            else if(kid == 2) {__unTransParentalChrs[1][index] = 1; __tdtTable[1][index] += heteroProb/(homoProb+heteroProb);} // -9X2 -> 2 ==> 1X2 -> 2
            else {;}
        }
        else {;}
    }
    else{;}
}


shuffleOutput unPhasedTrio::shuffle(std::string GenoHapo, gsl_rng* gslr)
{
	if(!__isInformative)
		halt("Error: can't shuffle on non-informative trio.");
	if(__PesudoGenoReady) {;}
	else { m_setUpPesudoGeno(gslr); __PesudoGenoReady=true;}
	vector<vReal> shuffledGeno;
	if(GenoHapo == "GenoShuffle") shuffledGeno = m_GenoShuffle(gslr);
	else if(GenoHapo == "HapoShuffle") shuffledGeno = m_HapoShuffle(gslr);
	else
		halt("Error: invalid shuffle method.");
	// phased to unphased conversion
	vector<vReal> geno(shuffledGeno.begin(), shuffledGeno.begin()+3);
	for(wsUint j=0; j<shuffledGeno[0].size(); j++){
		for(wsUint i=0; i<3; i++){
			if(shuffledGeno[2*i][j] < 0 || shuffledGeno[2*i+1][j] < 0) geno[i][j] = WISARD_NA;
			else geno[i][j] = shuffledGeno[2*i][j] + shuffledGeno[2*i+1][j];
		}
	}
	unPhasedTrio shuffleObj;
	shuffleObj.load(Cp_IO, *Xv_set, __tobeAnalyzed, __skipMiss,
		N_pat, N_mat, N_chi);
	return shuffleOutput(shuffleObj.getTdtTable(), shuffleObj.getNuTransPatChrs());
}

// for debug: shuffle and return shuffled genotypes
shuffleOutput unPhasedTrio::shuffleWithGenos(std::string GenoHapo, gsl_rng* gslr, vector<vReal>& genos)
{
	if(!__isInformative)
		halt("Error: can't shuffle on non-informative trio.");
	if(__PesudoGenoReady) {;}
	else { m_setUpPesudoGeno(gslr); __PesudoGenoReady=true;}
	vector<vReal> shuffledGeno;
	if(GenoHapo == "GenoShuffle") shuffledGeno = m_GenoShuffle(gslr);
	else if(GenoHapo == "HapoShuffle") shuffledGeno = m_HapoShuffle(gslr);
	else
		halt("Error: invalid shuffle method.");
	// phased to unphased conversion
	vector<vReal> geno(shuffledGeno.begin(), shuffledGeno.begin()+3);
	for(wsUint j=0; j<shuffledGeno[0].size(); j++){
		for(wsUint i=0; i<3; i++){
			if(shuffledGeno[2*i][j] < 0 || shuffledGeno[2*i+1][j] < 0) geno[i][j] = WISARD_NA;
			else geno[i][j] = shuffledGeno[2*i][j] + shuffledGeno[2*i+1][j];
		}
	}
	genos = geno;
	unPhasedTrio shuffleObj;
	shuffleObj.load(Cp_IO, *Xv_set, __tobeAnalyzed, __skipMiss,
		N_pat, N_mat, N_chi);
	return shuffleOutput(shuffleObj.getTdtTable(), shuffleObj.getNuTransPatChrs());
}


void unPhasedTrio::m_setUpPesudoGeno(gsl_rng* gslr)
{
	char** Na_geno = Cp_IO->getGenotype();
    //for(UINT i=0; i<2; i++) std::cout << __genotype[i] << std::endl;
	__pesudoGeno.resize(4);
	for(wsUint i=0; i<4; i++)
		__pesudoGeno[i].resize(Xv_set->size());

	wsUint j = 0;
	FOREACHDO (vInt_it, *Xv_set, J, j++) {
		char	N_gPat = Na_geno[N_pat][*J];
		char	N_gMat = Na_geno[N_mat][*J];
		char	N_gChi = Na_geno[N_chi][*J];
		if(!__tobeAnalyzed[j]){
			for(wsUint i=0; i<4; i++)
				__pesudoGeno[i][j]=WISARD_NA; // mark all not analyzed sites as -9
		}
		else if(N_gPat < 0 && N_gMat < 0) { for(wsUint i=0; i<4; i++) __pesudoGeno[i][j]=WISARD_NA;}    
		else if(N_gPat < 0 && N_gMat >= 0)
		{
			double guessed = guessMissing(N_gMat, N_gChi, gslr);
			if(guessed < 0) {for(wsUint i=0; i<4; i++) __pesudoGeno[i][j]=WISARD_NA; __tobeAnalyzed[j] = false;}
			else{
				if(BIN(gslr)) { __pesudoGeno[0][j] = guessed;  __pesudoGeno[1][j] = WISARD_NA;}
				else { __pesudoGeno[1][j] = guessed;  __pesudoGeno[0][j] = WISARD_NA;}
			}

			if(N_gMat == 0) {__pesudoGeno[2][j] = 0; __pesudoGeno[3][j] = 0;}
			else if(N_gMat == 1) {__pesudoGeno[2][j] = 0; __pesudoGeno[3][j] = 1;}
			else if(N_gMat == 2) {__pesudoGeno[2][j] = 1; __pesudoGeno[3][j] = 1;}
			else {for(wsUint i=0; i<4; i++) __pesudoGeno[i][j]=WISARD_NA; __tobeAnalyzed[j] = false;}
		}    
		else if(N_gPat >= 0 && N_gMat < 0)
		{
			double guessed = guessMissing(N_gPat, N_gChi, gslr);
			if(guessed < 0) {for(wsUint i=0; i<4; i++) __pesudoGeno[i][j]=WISARD_NA; __tobeAnalyzed[j] = false;}
			else{
				if(BIN(gslr)) {__pesudoGeno[2][j] = guessed; __pesudoGeno[3][j] = WISARD_NA;}
				else {__pesudoGeno[3][j] = guessed; __pesudoGeno[2][j] = WISARD_NA;}
			}

			if(N_gPat == 0) {__pesudoGeno[0][j] = 0; __pesudoGeno[1][j] = 0;}
			else if(N_gPat == 1) {__pesudoGeno[0][j] = 0; __pesudoGeno[1][j] = 1;}
			else if(N_gPat == 2) {__pesudoGeno[0][j] = 1; __pesudoGeno[1][j] = 1;}
			else {for(wsUint i=0; i<4; i++) __pesudoGeno[i][j]=WISARD_NA; __tobeAnalyzed[j] = false;}
		} 
		else if(N_gPat >= 0 && N_gMat >= 0)
		{
			if(N_gPat == 0) {__pesudoGeno[0][j] = 0; __pesudoGeno[1][j] = 0;}
			else if(N_gPat == 1) {__pesudoGeno[0][j] = 0; __pesudoGeno[1][j] = 1;}
			else if(N_gPat == 2) {__pesudoGeno[0][j] = 1; __pesudoGeno[1][j] = 1;}
			else {for(wsUint i=0; i<4; i++) __pesudoGeno[i][j]=WISARD_NA; __tobeAnalyzed[j] = false;}

			if(N_gMat == 0) {__pesudoGeno[2][j] = 0; __pesudoGeno[3][j] = 0;}
			else if(N_gMat == 1) {__pesudoGeno[2][j] = 0; __pesudoGeno[3][j] = 1;}
			else if(N_gMat == 2) {__pesudoGeno[2][j] = 1; __pesudoGeno[3][j] = 1;}
			else {for(wsUint i=0; i<4; i++) __pesudoGeno[i][j]=WISARD_NA; __tobeAnalyzed[j] = false;}
		} // no missing
		else {;}
	}
    //std::cout << "-----------\n";
    //for(UINT i=0; i<4; i++) std::cout << __pesudoGeno[i] << std::endl;
    //std::cout << "===========\n";
	return;
}

double unPhasedTrio::guessMissing(double knownParent, double kid, gsl_rng* gslr)
{
    double guessedHapo(-1);
    if(knownParent == 0 && kid == 0) guessedHapo = 0;
    else if(knownParent == 0 && kid == 1) guessedHapo = 1;
    else if(knownParent == 0 && kid == 2) guessedHapo = -1; // 0X-9 -> 2
    else if(knownParent == 1 && kid == 0) guessedHapo = 0;
    else if(knownParent == 1 && kid == 1) guessedHapo = (BIN(gslr))?1:0; // TODO: should be more explicit
    else if(knownParent == 1 && kid == 2) guessedHapo = 1;
    else if(knownParent == 2 && kid == 0) guessedHapo = -1; // 2X-9 ->0
    else if(knownParent == 2 && kid == 1) guessedHapo = 0;
    else if(knownParent == 2 && kid == 2) guessedHapo = 1;
    else {;}
    return guessedHapo;
}

#endif

} // End namespace ONETOOL
