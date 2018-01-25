//==========================================================================
// File:      wisard.cpp
//
// Author:    Sungyoung Lee
//
// History:   Initial implementation                              syl Mar 22
//
// Notes:     This source file implements wisard analyses.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include <math.h>
#include "app/app_wisard.h"
#include "global/Rconn.h"
#include "analyses/regr.h"
#include "analyses/corr.h"
#include "analyses/powercalc.h"
#include "analyses/window.h"
#include "analyses/qtest.h"
#include "analyses/emai.h"
#include "analyses/fqls.h"
#include "analyses/qls.h"
#include "analyses/fst.h"
#include "analyses/setmgr.h"
#include "analyses/tridge.h"
#include "analyses/ld.h"
#include "analyses/pca.h"
#include "analyses/mdr.h"
#include "analyses/corr.h"
#include "analyses/pddt.h"
#include "analyses/ppp.h"
#include "analyses/pharaoh.h"
#include "analyses/mfqls.h"
#include "analyses/raretdt.h"
#include "analyses/marker.h"
#include "analyses/design.h"
#include "analyses/gemma.h"
#include "analyses/gstest.h"
#include "analyses/grm.h"
#include "analyses/lod.h"
#include "analyses/ldprune.h"
#include "input/impute.h"
#include "utils/ext.h"

using namespace std;

#if TOOLSET_TYPE == TOOLSET_ONETOOL
char				S_exePath[MAX_PATH] ={ 0, };
#endif

namespace ONETOOL {

cPPPAnalysisV2*		Cp_PPP			= NULL;
cFamStrAnalysis*	Cp_anaFamStr	= NULL;
cAnalysis*			Cp_corr			= NULL;

void getFreqByCorr(cIO *Cp_IO, cAnalysis *Cp_ana)
{
//	wsUint	i;
	wsUint	N_sample	= Cp_IO->sizeSample();
//	wsReal	**Rp_corInv	= NULL;
	// Checked 
// 	wsReal	**Ra_1t		= NULL;
// 
// 	Ra_1t	= sseMatrix(1, N_sample);
// 	for (i=0 ; i<N_sample ; i++)
// 		Ra_1t[0][i] = W1;

	// Get correlation structure 
	cMatrix		*Mp_phiInv;
	if (OPT_ENABLED(indep)) 
		Mp_phiInv = new cIdtMatrix(N_sample);
	else if (IS_ASSIGNED(cor) || (!OPT_ENABLED(kinship) && !OPT_ENABLED(hybrid))) {
		cSymMatrix	M_phi(((cCorrAnalysis *)Cp_ana)->getCorr(), N_sample, 1);
		Mp_phiInv = &(M_phi.inv());
		// CORR 
// 		Rp_corInv = SO_invMatrix(,
// 			N_sample);
// 		if (Rp_corInv == NULL)
// 			halt("Can't get corInv matrix");
	} else {
		// PDDT 
		//wsReal **Ra_pddt = ((cPDDTAnalysis *)Cp_ana)->getPDDT();
		cSymMatrix	M_phi(((cPDDTAnalysis *)Cp_ana)->getPDDT(), N_sample, 1);
//		M_phiInv = M_phi.inv();
 		if (Cp_IO->getTwinSize() || OPT_ENABLED(ginv))
			Mp_phiInv = &(M_phi.ginv());
// 			SVDinverse(Ra_pddt, N_sample, &Rp_corInv);
 		else {
			Mp_phiInv = &(M_phi.inv());
// 			Rp_corInv = SO_invMatrix(Ra_pddt, N_sample);
// 			if (Rp_corInv == NULL)
// 				halt("Can't get corInv matrix");
 		}
	}

	//       T -1  
	//  Get (1 Φ  )  = 1*n matrix 
	// Checked 
// 	wsReal **Ra_1tPhi = sseMpM(Ra_1t, 1, N_sample, Rp_corInv, N_sample,
// 		N_sample);
// 	sseUnmat(Rp_corInv, N_sample);
// 	sseUnmat(Ra_1t, 1);
// 
// 	//       T -1  -1
// 	// Get (1 Φ  1) = Scalar  == Sum of phiInv
// 	// Checked 
// 	wsReal **Ra_corDenom = sseMpMt(Ra_1tPhi, 1, N_sample, Ra_1t, 1, N_sample);
// 	Ra_corDenom[0][0] = W1 / Ra_corDenom[0][0];
	wsReal		R_corDenom	= W0;
	cVector		V_pIsum		= Mp_phiInv->sumR(&R_corDenom);
	delete Mp_phiInv;
	wsReal		*Ra_pIsum	= V_pIsum.get();

	//       T -1  -1  T -1
	// Get (1 Φ  1)  (1 Φ  X) = 1*p matrix 
	wsUint		N_vrt	= Cp_IO->sizeVariant();
	xMaf*		Xp_maf	= Cp_IO->getMAF();
// 	sseMpC(Ra_1tPhi, Ra_corDenom[0][0], Ra_1tPhi, 1, N_sample);
// 	deallocMatrix(Ra_corDenom, 1);

	// Checked 
	cExporter*	Cp_frq	= cExporter::summon("estimated.maf");
	LOGoutput("Estimated MAF from sample relatedness is exported to "
		"[%s.estimated.maf]\n", OPT_STRING(out));
	vVariant		Xa_SNP = Cp_IO->getVariant();
	// Export header 
	if (IS_ASSIGNED(annogene))
		Cp_frq->put("CHR	SNP	ANNOT	MAJOR	MINOR	pHat	pObs\n");
	else
		Cp_frq->put("CHR	SNP	MAJOR	MINOR	pHat	pObs\n");

	wsReal *Ra_buf = sseVector(N_sample);
	LOOP (i, N_vrt) {
		xVariant	X_SNP	= Xa_SNP[i];
		xMaf&		X_maf	= Xp_maf[i];

		setGeno(Cp_IO, i, Ra_buf, NULL);
		wsReal R_est = sseVV(Ra_pIsum, N_sample, Ra_buf) / R_corDenom /
			W2;

		// Export calculated frequency 
		wsStrCst S_chr = getChrName2(X_SNP.chr);
		char S_bufAnno[512] = { 0, };
		if (IS_ASSIGNED(annogene))
			sprintf(S_bufAnno, "	%s", X_SNP.anno);
		if (OPT_ENABLED(indel)) {
			Cp_frq->fmt("%s	%s%s	%s	%s	%g	%g\n", S_chr, X_SNP.name, S_bufAnno,
				X_SNP.indel1, X_SNP.indel2 ? X_SNP.indel2 : "<NA>",
				R_est, X_maf.R_maf);
		} else {
			if (X_SNP.al2)
				Cp_frq->fmt("%s	%s%s	%c	%c	%g	%g\n", S_chr, X_SNP.name, S_bufAnno,
					X_SNP.al1, X_SNP.al2, R_est, X_maf.R_maf);
			else
				Cp_frq->fmt("%s	%s%s	%c	<NA>	%g	%g\n", S_chr, X_SNP.name, S_bufAnno,
					X_SNP.al1, R_est, X_maf.R_maf);
		}
	}

	delete Cp_frq;
}

void _exportWeight(cIO *io, cPPPAnalysisV2* Cp_anaPPP)
{
	wsUint		N_SNP	= io->sizeVariant();
	vVariant&	Xv_vrt	= io->getVariant();
	cExporter*	Cp_wgt	= cExporter::summon("weight.lst");
	wsVecCst		Ra_wgt	= io->getWeight(Cp_anaPPP);
	LOGoutput("Loaded/computed variant weights are exported to "
		"[%s.weight.lst]\n", OPT_STRING(out));

	// Print header 
	headerVariant(Cp_wgt);
	Cp_wgt->put("	WEIGHT\n");

	// For identical weight, just print all 1.0s 
	if (OPT_ENABLED(noweight)) for (wsUint i=0 ; i<N_SNP ; i++) {
		entryVariant(Cp_wgt, Xv_vrt[i]);
		Cp_wgt->put("	1.0\n");
	}
	// If there is customweight buffer 
	else if (Ra_wgt) for (wsUint i=0 ; i<N_SNP ; i++) {
		entryVariant(Cp_wgt, Xv_vrt[i]);
		Cp_wgt->fmt("	%g\n", Ra_wgt[i]);
	}
	// Otherwise, use sqrt2pq 
	else {
		wsRealCst *Ra_sqrt2pq = getPPP(io)->getPPPsq();
		for (wsUint i=0 ; i<N_SNP ; i++) {
			entryVariant(Cp_wgt, Xv_vrt[i]);
			Cp_wgt->fmt("	%g\n", W1/Ra_sqrt2pq[i]);
		}
	}

	delete Cp_wgt;
}

cPPPAnalysisV2* getPPP(cIO *io)
{
	cTimer t2;
	
	if (Cp_anaFamStr == NULL)
		halt_fmt(WISARD_SYST_NULL_FAMSTRUCT);

	if (Cp_PPP == NULL) {
		t2.start();
		Cp_PPP = new cPPPAnalysisV2(io, Cp_anaFamStr);

		Cp_PPP->run();
		if (io->sizeVariant())
			LOG("[%s] PPP calculated\n", t2.getReadable());
		io->setPPP(Cp_PPP);
	}

	return Cp_PPP;
}

int onetool_analysis(cIO *io)
{
	cTimer t2;

	if (io == NULL) {
		cCorrAnalysis *Cp_anaCorr = NULL;

		// Calculated correlation load 
		if (IS_ASSIGNED(cor)) {
			t2.start();
			Cp_anaCorr = new cCorrAnalysis(io, OPT_STRING(cor));
			LOG("[%s] Pre-calculated correlation loaded\n", t2.getReadable());
		}

		// Do sample-size analysis 
		if (OPT_ENABLED(powercalc)) {
			t2.start();
			cPowerCalcAnalysis C_anaSS(NULL, Cp_anaCorr);
			C_anaSS.run();
			LOG("[%s] Sample size estimation performed\n", t2.getReadable());
		}
		if (OPT_ENABLED(powercalc2)) {
			t2.start();
			cPowerCalcAnalysisV2 C_anaSS(NULL, Cp_anaCorr);
			C_anaSS.run();
			LOG("[%s] Sample size estimation performed\n", t2.getReadable());
		}

		LOG("No dataset were found, stop analyses\n");
		return 1;
	}
	
	// Pedigree plot, if --plot
#ifdef USE_R
	if (OPT_ENABLED(plot)) {
		R().init();
		if (io->getFamilyData().size() > 1) {
			if (R().isInit()) {
				vSampPtr& Xa_samp = io->getSample();
				wsVecCst Ra_phe = io->getPheno();
				LOGnote("Try to draw pedigree plot with [%d] families\n", io->getFamilyData().size());
				string out = "pedigree.data <- matrix(c(";
				int j = 1;
				FOREACHDO (vSampPtr_it, Xa_samp, i, j++) {
					xSample& I = **i;
					out += "\"";
					out += I.S_FID;
					out += "\",\"";
					out += I.S_IID;
					if (I.Xp_pat) {
						out += "\",\"";
						out += I.Xp_pat->S_IID;
						out += "\",";
					} else out += "\",NA,";
					if (I.Xp_mat) {
						out += "\"";
						out += I.Xp_mat->S_IID;
						out += "\",";
					} else out += "NA,";
					out += I.N_sex == 1 ? "1" : (I.N_sex == 2 ? "2" : "NA");
					out += ",";
					out += !Ra_phe || WISARD_NA_REAL == Ra_phe[I.N_idx] || io->isContinuous() ? "NA" : (
						Ra_phe[I.N_idx] == WISARD_AFFECTED ? "2" : "1"
					);
					if (j != Xa_samp.size())
						out += ",";
				}
				out += "), ncol=6, byrow=T);if(!exists(\"onetool.plot.pedigree\")){cat(\"No function for pedigree plot found!\\n\");}else onetool.plot.pedigree(pedigree.data, \"";
				out += OPT_STRING(out);
				out += ".pedigree.pdf\")";

				R().Reval(out.c_str());
			} else LOGwarn("Cannot draw plot because R is not properly initialized!\n");
		}
	}
#endif

	// Perform variant addressing 
	cVariantMap	*Cp_anaMM = new cVariantMap(io);

	// Init geneset 
	cSetManagerAnalysis	*Cp_anaSetMgr = NULL;
	t2.start();
	// --setconsec & --set & --setrandom are mut. excl. 
	if (IS_ASSIGNED(setconsec) && IS_ASSIGNED(set))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--setconsec", "--set");
	else if (IS_ASSIGNED(setconsec) && IS_ASSIGNED(setrandom))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--setconsec", "--setrandom");
	else if (IS_ASSIGNED(set) && IS_ASSIGNED(setrandom))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--set", "--setrandom");
	if (IS_ASSIGNED(setconsec) || IS_ASSIGNED(setrandom)) {
		// No --genesize is allowed 
		if (IS_ASSIGNED(genesize) && IS_ASSIGNED(setconsec))
			halt_fmt(WISARD_CANT_EXCL_OPT, "--setconsec", "--genesize");
		if (IS_ASSIGNED(genesize) && IS_ASSIGNED(setrandom))
			halt_fmt(WISARD_CANT_EXCL_OPT, "--setrandom", "--genesize");

		Cp_anaSetMgr = new cSetManagerAnalysis(io, Cp_anaMM);
		Cp_anaSetMgr->run();
	} else if (IS_ASSIGNED(set) || IS_ASSIGNED(geneset)) {
		Cp_anaSetMgr = new cSetManagerAnalysis(io);
		Cp_anaSetMgr->run();
	}

	// If --genesplit 
	if (OPT_ENABLED(genesplit)) {
		if (IS_ASSIGNED(genoprob)) LOGwarn("Genotype probability format does not support conversion, skips --genesplit");
		else io->exportBEDgene(1, Cp_anaSetMgr);
	}

	// Print gene summary 
	if (Cp_anaSetMgr) {
		mGeneDef &Xm_gdef = Cp_anaSetMgr->getGeneDef();
		// Filtering gene-sets if needed 
		for (mGeneDef_it i=Xm_gdef.begin() ; i!=Xm_gdef.end() ; ) {
			if (IS_ASSIGNED(genesize) && !isInRange(OPT_RANGE(genesize),
				(wsReal)i->second.size())) {
					LOG("Gene set `%s`(%d SNPs) excluded\n",
						i->first.c_str(), i->second.size());
					Xm_gdef.erase(i++);
			} else
				i++;
		}
		if (OPT_RANGE(genesize).R_s == W0) {
			if (OPT_RANGE(genesize).R_e == W0)
				LOG("[%s] %d gene set loaded\n", t2.getReadable(), (int)Xm_gdef.size());
			else
				LOG("[%s] %d gene set loaded (Permitted range ~%g)\n",
					t2.getReadable(), (int)Xm_gdef.size(), OPT_RANGE(genesize).R_e);
		} else if (OPT_RANGE(genesize).R_e == W0)
			LOG("[%s] %d gene set loaded (Permitted range %g~)\n",
				t2.getReadable(), (int)Xm_gdef.size(), OPT_RANGE(genesize).R_s);
		else
			LOG("[%s] %d gene set loaded (Permitted range %g~%g)\n",
				t2.getReadable(), (int)Xm_gdef.size(), OPT_RANGE(genesize).R_s,
				OPT_RANGE(genesize).R_e);

		// Check whether there is remained gene set 
		// 150824 --expression do not require gdef
		if (Xm_gdef.size() == 0 && !IS_ASSIGNED(expression))
			halt_fmt(WISARD_NULL_GENEDEF);

		if (OPT_ENABLED(genesummary) || OPT_ENABLED(gmapsummary))
			Cp_anaSetMgr->summary();
		if (OPT_ENABLED(makeset))
			Cp_anaSetMgr->exportSet();
	}

	// Perform family structure-related analysis 
	Cp_anaFamStr = new cFamStrAnalysis(io);
	Cp_anaFamStr->run();

	/* Perform imputation */
#if TOOLSET_TYPE != TOOLSET_ONETOOL
	if (OPT_ENABLED(impute)) {
		getPPP(io);
		io->impute(getCorrelation(io));
		// Export if required 
		if (OPT_ENABLED(makeimpute)) {
			LOGoutput("Imputed dataset is exported to [%s.imputed.dosage]", OPT_STRING(out));
			cImpute *Cp_imp = io->getImpute();
			exportMatrix("imputed.dosage", Cp_imp->getData(),
				io->sizeVariant(), io->sizeSample());
		}
	}
#endif

	//	vector<SAMPLE>	&samp = io->getSamples();

	// 	FILE *temp = fopen("temp.txt", "w+");
	// 	fprintf(temp, "Sample size %d\n", (int)(samp.size()));
	// 	for(vector<SAMPLE>::iterator i=samp.begin(); i!=samp.end(); i++)
	// 	{
	// //		notice("%s	", (*i).iid.c_str());
	// 		fprintf(temp, "%s	", (*i).iid.c_str());
	// 	}
	// 	fclose(temp);
	LOG("%d founders found\n", io->sizeFounder());
//	R().init();

	/* Perform LD prune analysis */
	if (IS_ASSIGNED(prunevif) || IS_ASSIGNED(prunepw)) {
		t2.start();
		cLdPruneAnalysis C_anaLP(io);
		C_anaLP.run();
		LOG("[%s] LD pruning analysis performed\n", t2.getReadable());
	}

	// Perform sample-wise analysis 
	if (OPT_ENABLED(het) || OPT_ENABLED(hethom)) {
		cExporter*	Cp_hh	= NULL;
		cExporter*	Cp_het	= NULL;
		xMaf*		Xp_maf	= io->getMAF();
		vVariant&	Xv_vrt	= io->getVariant();
		vSampPtr&	Xv_samp	= io->getSample();
		char**		Na_data	= io->getGenotype();

		// --hethom 
		if (OPT_ENABLED(hethom)) {   
			Cp_hh = cExporter::summon("hethom.res");
			LOGoutput("Sample-site het/hom ratios is exported to [%s.hethom.res]\n",
				OPT_STRING(out));
			// Export header 
			Cp_hh->put("FID	IID	NOBS	NHET	HETHOM\n");
		}
		// --het 
		if (OPT_ENABLED(het)) {
			Cp_het = cExporter::summon("het.res");
			LOGoutput("Sample heterozygosity is exported to [%s.het.res]\n",
				OPT_STRING(out));
			// Export header 
			Cp_het->put("FID	IID	O(HOM)	E(HOM)	N(NM)	F\n");
		}

		wsUint i = 0;
		FOREACHDO (vSampPtr_it, Xv_samp, I, i++) {
			wsUint	N_hhDenom	= 0;
			wsUint	N_hhNumer	= 0;
			wsUint	j			= 0;
			wsReal	R_hE		= W0;

			FOREACHDO (vVariant_it, Xv_vrt, J, j++)  {
				xMaf& X_maf = Xp_maf[j];

				char N_geno = Na_data[i][j];
				if (isMissing(N_geno)) continue;

				// --het 
				if (isAutosome(*J) && X_maf.N_allMac>1) {
					R_hE += W1 - (W2*X_maf.R_allMaf*(W1-X_maf.R_allMaf))
						* (wsReal)(X_maf.N_allMac/(X_maf.N_allMac-1));

					// --hethom 
					if (N_geno == 1) N_hhNumer++;
					N_hhDenom++;
				}
			}

			// Export statistics 
			if (Cp_hh) Cp_hh->fmt("%s	%s	%d	%d	%g\n", (*I)->S_FID.c_str(),
				(*I)->S_IID.c_str(), N_hhDenom, N_hhNumer,
				N_hhNumer/(wsReal)(N_hhDenom-N_hhNumer));
			if (Cp_het) Cp_het->fmt("%s	%s	%d	%g	%d	%g\n", (*I)->S_FID.c_str(),
				(*I)->S_IID.c_str(), N_hhDenom-N_hhNumer, R_hE,
				N_hhDenom, ((wsReal)(N_hhDenom-N_hhNumer) - R_hE) / (N_hhDenom-R_hE));
		}

		if (Cp_hh) delete Cp_hh;
		if (Cp_het) delete Cp_het;
	}
	// Perform window-based analysis 
	if (IS_ASSIGNED(density) || IS_ASSIGNED(tstv)) {
		t2.start();
		cWindowAnalysis C_anaWnd(io);
		C_anaWnd.run();
		LOG("[%s] Window-based analysis performed\n", t2.getReadable());
	}
	// Perform MDR analysis 
	if (OPT_ENABLED(mdr) || OPT_ENABLED(gmdr)) {
		t2.start();
		cMdrAnalysis C_anaMDR(io, Cp_anaSetMgr);
		C_anaMDR.run();
		LOG("[%s] MDR analysis performed\n", t2.getReadable());
	}
	// Perform LD contrast method 
	if (OPT_ENABLED(ldcontrast)) {
		t2.start();
		// FIXME : Impl 
		halt("This function is under development!");
		LOG("[%s] LD contrast analysis performed\n", t2.getReadable());
	}
	// Perform epistasis analysis (BOOST and PLINK --fast-epistasis)
	if (OPT_ENABLED(boost) || OPT_ENABLED(quickepi)) {
		t2.start();
		cBoostAnalysis C_anaBoost(io);
		C_anaBoost.run();

		if (OPT_ENABLED(boost))
			LOG("[%s] BOOST analysis performed\n", t2.getReadable());
		else if (OPT_ENABLED(quickepi))
			LOG("[%s] Quick epistasis analysis performed\n", t2.getReadable());
		else
			LOG("[%s] Epistasis analyses performed\n", t2.getReadable());
	}
#if TOOLSET_TYPE != TOOLSET_ONETOOL
	// Perform LD analysis 
	if (OPT_ENABLED(ld)) {
		t2.start();
		cLdAnalysis C_LD(io, Cp_anaMM);
		C_LD.run();
		LOG("[%s] LD analysis performed\n", t2.getReadable());
	}
	// Perform LoD analysis 
	if (OPT_ENABLED(lod)) {
		t2.start();
		cLodAnalysis C_LoD(io, Cp_anaFamStr, Cp_anaSetMgr);
		C_LoD.run();
		LOG("[%s] LoD analysis performed\n", t2.getReadable());
	}
	// Perform truncated ridge
	if (OPT_ENABLED(tridge)) {
		t2.start();
		cTridgeAnalysis C_tridge(io);
		C_tridge.run();
		LOG("[%s] Truncated ridge analysis performed\n", t2.getReadable());
	}
#endif
	// Perform Fst 
	if (IS_ASSIGNED(fst)) {
		t2.start();
		cFstAnalysis C_Fst(io, Cp_anaSetMgr);
		C_Fst.run();
		LOG("[%s] Fst analysis performed\n", t2.getReadable());
	}
	// Perform TDT
	if (OPT_ENABLED(tdt)) {
		t2.start();
		cTDTanalysis C_tdt(io, Cp_anaFamStr);
		C_tdt.run();
		LOG("[%s] TDT analysis performed\n", t2.getReadable());
	}
	// Perform SDT 
	if (OPT_ENABLED(sdt)) {
		t2.start();
		cSDTanalysis C_sdt(io, Cp_anaFamStr);
		C_sdt.run();
		LOG("[%s] Sibship disequilibrium test performed\n", t2.getReadable());
	}
	// Perform R analysis 
	if (IS_ASSIGNED(R)) {
#ifdef USE_R
		cStdMatrix U;
#else
#endif
	}
	// Perform variant-level analyses 
	if (OPT_ENABLED(famuniq) || OPT_ENABLED(monotone) || OPT_ENABLED(singleton) ||
		OPT_ENABLED(doubleton) || OPT_ENABLED(fisher) || OPT_ENABLED(trend) ||
		IS_ASSIGNED(popuniq)) {
		t2.start();
		cVariantAnalysis C_mk(io, Cp_anaFamStr);
		C_mk.run();
		LOG("[%s] Variant-target analysis performed\n", t2.getReadable());
	}

	if (IS_ASSIGNED(grm)) {
		// Get GRM matrix 
		ASSERT_OPTION(kinship);

		t2.start();
		cGRMAnalysis	C_anaGRM(io, (cPDDTAnalysis *)getCorrelation(io));
		C_anaGRM.run();
		LOG("[%s] GRM comparison complete\n", t2.getReadable());
	} else {
		// Calc PPP 
		if (io->sizeVariant()) {
			// Calc PCA 
			if (OPT_ENABLED(pca) || OPT_ENABLED(mds)) {
				t2.start();
				cPcaAnalysisV2 C_anaPCA(io, getPPP(io), Cp_anaFamStr);

				C_anaPCA.run();
				LOG("[%s] Principal Component-related analysis performed\n", t2.getReadable());
			}
		}

		// Compute SNP-blup 
		if (IS_ASSIGNED(variantblup)) {
			cMatrix*	Mp_phi		= NULL;
			wsUint		N_samp		= io->sizeSample();
			wsUint		N_vrt		= io->sizeVariant();
			wsUint		N_anaSamp	= 0;
			char*		Ba_miss		= io->getPheCovMissing(&N_anaSamp);
			wsReal		R_rho		= OPT_REAL(variantblup);

			/* Get phi matrix */ {
				cAnalysis*	Cp_ana		= getCorrelation(io);
				if (OPT_ENABLED(indep))
					Mp_phi = new cIdtMatrix(N_anaSamp);
				else if (IS_ASSIGNED(cor) || (!OPT_ENABLED(kinship) && !OPT_ENABLED(hybrid)))
					Mp_phi = new cSymMatrix(((cCorrAnalysis *)Cp_ana)->getCorr(), N_anaSamp, 1);
				else
					Mp_phi = new cSymMatrix(((cPDDTAnalysis *)Cp_ana)->getPDDT(), N_anaSamp, 1);
			}

			cMatrix& M_V = Mp_phi->addDiag(R_rho);
			cMatrix& M_Vinv = M_V.inv();
			delete Mp_phi;

			// Construct X' 
			wsUint N_cov = io->sizeCovar();
			LOGnote("Compute variant-level BLUP with rho [%g] and [%d] covariates\n",
				R_rho, N_cov);

			wsUint N_x = 1 + N_cov;
			cStdMatrix M_Xt(N_x, N_anaSamp);
			wsMat Ra_Xt = M_Xt.get();
			M_Xt.setRow(0, W1);
			wsMat Ra_cov = io->getCovariates();
			LOOP (i, N_cov) {
				wsUint J = 0;
				LOOP (j, N_samp) {
					if (Ba_miss[j]) continue;
					Ra_Xt[i+1][J++] = Ra_cov[i][j];
				}
			}
			cVector V_y(N_anaSamp);
			/* Construct Y */ {
				wsRealCst* Ra_oriY = io->getPheno();
				wsVec Ra_y = V_y.get();
				wsUint J = 0;
				LOOP (j, N_samp) {
					if (Ba_miss[j]) continue;
					Ra_y[J++] = Ra_oriY[j];
				}
			}

			// Compute X'V^-1X 
			cSymMatrix M_XVX = M_Xt.MMt(M_Vinv);
			cSymMatrix& M_XVXinv = M_XVX.inv();
			cVector V_Vy = M_Vinv * V_y;
			cVector V_XVy = M_Xt * V_Vy;
			cVector V_beta = M_XVXinv * V_XVy;
			delete& M_XVXinv;
			cVector V_pred = V_beta * M_Xt;
			cVector V_res = V_y - V_pred;

			// Compute V^-1(Y-Xb^hat) 
			cVector V_Vr = M_Vinv * V_res;
			delete& M_Vinv;

			// Get MAF 
			xMaf* Xa_maf = io->getMAF();

			// Compute BLUP 
			cVector V_blup(N_vrt);
			wsVec Ra_blup = V_blup.get();
			char** Na_data = io->getGenotype();
			cVector V_g(N_anaSamp);
			wsVec Ra_g = V_g.get();
			LOOP (i, N_vrt) {
				// Denominator
				wsReal R_maf = Xa_maf[i].R_allMaf;
				wsReal R_denom = W1 / sqrt(W2*R_maf * (W1 - R_maf));
				wsReal R_2maf = W2 * R_maf;

				wsUint J = 0;
				LOOP (j, N_samp) {
					if (Ba_miss[j]) continue;
					Ra_g[J++] = (wsReal)Na_data[j][i];
				}
				cVector V_gs = V_g - R_2maf;
				Ra_blup[i] = V_gs.sum(V_Vr) * R_denom / sqrt((wsReal)N_vrt);
			}
			cExporter* Cp_blup = cExporter::summon("variant.blup");
			LOGoutput("Variant-level BLUP is exported to [%s.variant.blup]\n",
				OPT_STRING(out));
			headerVariant(Cp_blup);
			Cp_blup->put("	BLUP\n");
			vVariant& Xv_vrt = io->getVariant();
			LOOP (i, N_vrt) {
				entryVariant(Cp_blup, Xv_vrt[i]);
				Cp_blup->fmt("	%g\n", Ra_blup[i]);
			}
			delete Cp_blup;
		}

		// Covariate export moved to here since PCA adds some covariates 
		if (OPT_ENABLED(makecov))
			io->exportCovariates();

		// Perform regression analysis 
		if (OPT_ENABLED(regression)) {
			t2.start();
			cRegrAnalysis C_anaRegr(io, Cp_anaSetMgr);
			C_anaRegr.run();
			LOG("[%s] Regression analysis performed\n", t2.getReadable());
		}

		// Correlation calc 
		//		cFinalStatAnalysisV2	*Cp_anaFS;
		cTestScoreAnalysis	*Cp_inpScore = NULL;

		// --freq by correlation 
		if (IS_ASSIGNED(freq) && !stricmp(OPT_STRING(freq), "blue")) {
			t2.start();
			getFreqByCorr(io, getCorrelation(io));
			LOG("[%s] MAF calculation by correlation\n", t2.getReadable());
		}

		// Do sample-size analysis 
		if (OPT_ENABLED(powercalc)) {
			t2.start();
			cPowerCalcAnalysis C_anaSS(io, getCorrelation(io));
			C_anaSS.run();
			LOG("[%s] Sample size estimation performed\n", t2.getReadable());
		}
		if (OPT_ENABLED(powercalc2)) {
			t2.start();
			cPowerCalcAnalysisV2 C_anaSS(io, getCorrelation(io));
			C_anaSS.run();
			LOG("[%s] Sample size estimation performed\n", t2.getReadable());
		}
		// Do study-related analyses 
		if (OPT_ENABLED(sxa)) {
			t2.start();
			cStudyAnalysis C_anaSA(io, getCorrelation(io));
			C_anaSA.run();
			LOG("[%s] Sample selection procedure performed\n", t2.getReadable());
		}

		if (OPT_ENABLED(qls)) {
			// Perform QLS test 
			t2.start();
			cQlsAnalysis C_anaQls(io, getPPP(io), getCorrelation(io));
			C_anaQls.run();

			LOG("[%s] QLS test performed\n", t2.getReadable());
		}
		if (OPT_ENABLED(multifqls)) {
			getPPP(io);

			// Solo-MQLS 
			t2.start();
			cMFqlsAnalysis C_anaMqls(io, getCorrelation(io), Cp_anaSetMgr,
				Cp_anaMM);
			C_anaMqls.run();

			LOG("[%s] Multi family QLS (MFQLS) test performed\n", t2.getReadable());
		} else if (OPT_ENABLED(fqls) || OPT_ENABLED(mqls)) {
			t2.start();
			cEmAiAnalysisV2 *Cp_nullEmAi = NULL;

			if (OPT_ENABLED(fqls))
				ASSERT_OPTION(heri);

			if (OPT_ENABLED(fqls) && !(io->isContinuous()))
				ASSERT_OPTION(prevalence);

			// --prevalence & cont. pheno are m.e. 
			if (IS_ASSIGNED(prevalence) && io->isContinuous())
				halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--prevalence",
				"continuous phenotype");

			// If the prevalence is not assigned, using BLUP 
			// --logistic then dont' need 
			if (!IS_ASSIGNED(prevalence)) {
				if (Cp_inpScore == NULL) {
					if (!IS_ASSIGNED(blup))
						LOG("--prevalence is not given, BLUP will be calculated"
						" in order to substitute prevalence information\n");
					Cp_inpScore = new cTestScoreAnalysis(io, getCorrelation(io),
						Cp_anaFamStr, !OPT_ENABLED(ml));
				}
				Cp_nullEmAi = Cp_inpScore->getNullEmAi();
			}

			// Perform family QLS test - 
			t2.start();
			cFamQlsAnalysis C_anaFamQls(io, getPPP(io), Cp_anaFamStr,
				getCorrelation(io), Cp_nullEmAi);
			C_anaFamQls.run();

			if (!OPT_ENABLED(fqls) && !OPT_ENABLED(fastfqls))
				LOG("[%s] MQLS test performed\n", t2.getReadable());
			else
				LOG("[%s] Family QLS test performed\n", t2.getReadable());
		}
		cFemmaAnalysis *Cp_anaFemma = NULL;
		if (OPT_ENABLED(gemma)) {
			// Perform FEMMA analysis 
			t2.start();
			Cp_anaFemma = new cFemmaAnalysis(io, getCorrelation(io));
			Cp_anaFemma->run();

			LOG("[%s] GEMMA performed\n", t2.getReadable());
		}
		if (OPT_ENABLED(qtest)) {
			// Perform Q-test analysis 
			t2.start();
			cQtestAnalysis C_anaQtest(io, getPPP(io), Cp_anaSetMgr, NULL);
			C_anaQtest.run();

			LOG("[%s] Q-test performed\n", t2.getReadable());
		}
		if (OPT_ENABLED(pharaoh) || OPT_ENABLED(proopt) || OPT_ENABLED(gesca)) {
			// Perform PHARAOH analysis 
			t2.start();
			cPharaohAnalysis C_anaGsca(io, Cp_anaSetMgr, getPPP(io));
			C_anaGsca.run();

			LOG("[%s] PHARAOH performed\n", t2.getReadable());
		}
		if (OPT_ENABLED(rvtdt) || OPT_ENABLED(fbskat)) {
			// Perform rvTDT analysis 
			t2.start();
			cRvtdtAnalysis C_anaRvtdt(io, Cp_anaSetMgr, Cp_anaFamStr, false, false, 0.05);
			C_anaRvtdt.run();

			LOG("[%s] rvTDT performed\n", t2.getReadable());
		}
		if (OPT_ENABLED(scoretest) || OPT_ENABLED(heritability) ||
			OPT_ENABLED(makeblup)) {
			// Perform score test 
			t2.start();
			Cp_inpScore = new cTestScoreAnalysis(io, getCorrelation(io),
				Cp_anaFamStr, !OPT_ENABLED(ml));
			Cp_inpScore->run();

			if (OPT_ENABLED(scoretest))
				LOG("[%s] Score test performed\n", t2.getReadable());
			if (OPT_ENABLED(heritability) || OPT_ENABLED(makeblup))
				LOG("[%s] Heritability-related analyses performed\n", t2.getReadable());
		}
		if (OPT_ENABLED(genetest) || OPT_ENABLED(skato)) {
			ASSERT_OPTS_AND(set, setconsec);

			t2.start();
			cEmAiAnalysisV2 *Cp_nullEmAi = NULL;

			// --prevalence & cont. pheno are m.e. 
			if (IS_ASSIGNED(prevalence) && io->isContinuous())
				halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--prevalence",
					"continuous phenotype");

			// If the prevalence is not assigned, using BLUP 
			// --logistic then dont' need 
			if (!IS_ASSIGNED(prevalence) && !IS_ASSIGNED(longitudinal) &&
				!OPT_ENABLED(logistic)) {
				if (Cp_inpScore == NULL) {
					if (!IS_ASSIGNED(blup))
						LOG("--prevalence is not given, BLUP will be calculated"
							" in order to substitute prevalence information\n");
					Cp_inpScore = new cTestScoreAnalysis(io, getCorrelation(io),
						Cp_anaFamStr, !OPT_ENABLED(ml));
				}
				Cp_nullEmAi = Cp_inpScore->getNullEmAi();
			}

			if (OPT_ENABLED(skato)) {
				// FIXME : Only against for single phenotype 
				if (io->isContinuous(0) && io->sizePheno() == 1) {
					if (Cp_inpScore == NULL) {
						LOG("--scoretest disabled, null model will be evaluated"
							" in order to perform gene-level test for continuous"
							" trait\n");
						Cp_inpScore = new cTestScoreAnalysis(io, getCorrelation(io),
							Cp_anaFamStr, !OPT_ENABLED(ml));
					}

					Cp_nullEmAi = Cp_inpScore->getNullEmAi();
				}
			}
			cGSTestAnalysis C_anaGS(io, getPPP(io), getCorrelation(io),
				Cp_anaSetMgr, Cp_nullEmAi);
			C_anaGS.run();
			LOG("[%s] Gene-level test performed\n", t2.getReadable());
		}

		/* If explicit correlation export is defined */
		if (OPT_ENABLED(corpair) || IS_ASSIGNED(makecor)) {
			getPPP(io);
			getCorrelation(io);
		}

		/* Perform external analyses */
#if TOOLSET_TYPE != TOOLSET_ONETOOL
		if (IS_ASSIGNED(ext)) {
			wsUint	N_ext	= 0;
			char**	S_exts	= loadStringValues(OPT_STRING(ext), &N_ext);

			LOG("[%d] external scripts found, try to execute...\n", N_ext);

			/* Perform external analyses */
			for (wsUint i=0 ; i<N_ext ; i++) {
				/* Find [ */
				char *Sp_optStart = strchr(S_exts[i], '[');
				if (Sp_optStart) {
					char *Sp_optEnd = strchr(Sp_optStart+1, ']');
					if (Sp_optEnd == NULL)
						halt("Option delimiter '[' assigned, but not closed with ']'");
				}

 				doWISARDext(S_exts[i]);
			}
		}
#endif
		if (Cp_inpScore)
			delete Cp_inpScore;
		if (Cp_anaFemma)
			delete Cp_anaFemma;
	}

#if TOOLSET_TYPE != TOOLSET_ONETOOL
	/* At this point PPP now can be calculated */
	if (OPT_ENABLED(makeweight))
		_exportWeight(io, getPPP(io));
#endif
	if (Cp_anaFamStr) {
		delete Cp_anaFamStr;
		Cp_anaFamStr = NULL;
	}
	if (Cp_anaSetMgr) {
		delete Cp_anaSetMgr;
		Cp_anaSetMgr = NULL;
	}
	if (Cp_anaMM) {
		delete Cp_anaMM;
		Cp_anaMM = NULL;
	}
	if (Cp_PPP) {
		delete Cp_PPP;
		Cp_PPP = NULL;
	}

	// no error
	return 0;
}

} // end of namespace ONETOOL
