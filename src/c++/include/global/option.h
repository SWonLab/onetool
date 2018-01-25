// 120119
// Sungyoung Lee, BIBS
//
// cOption.h
#pragma once
#ifndef __WISARD_OPTION_H__
#define __WISARD_OPTION_H__

#include "global/common.h"
#include "utils/logic.h"

namespace ONETOOL {

struct xOptRange
{
	char	R_sEQ, R_eEQ;
	wsReal	R_s, R_e;
};

enum eOptType {
	OT_ERROR,
	OT_ONOFF,
	OT_NUMBER,
	OT_RANGE,
	OT_STRING,
	OT_REAL,
	OT_EXPR,
};

typedef union _xParam
{
	int			N_intVal;
	xOptRange	X_rngVal;
	wsReal		R_realVal;
	char*		S_strVal;
	xOperation	X_expVal;
} xParam;

typedef enum _xPossibleRange {
	RT_POS,	/* (0,inf) */
	RT_PROPNUM,	/* (0,1) or [1,2,...] */
	RT_NONNEG,	/* [0,inf) */
	RT_M01,		/* (0,1] */
	RT_01,		/* [0,1] */
	RT_ONOFF,	/* 0 or 1 */
	RT_NUM,		/* (-inf,inf) */
	RT_ABS1,	/* [-1,1] */
	RT_OTHER,
} xPossibleRange;

struct xOption {
	bool			B_isApplyDef;
	bool			B_isAssigned;
	char			*S_longName;
	char			*S_longSynonym;
	char			*S_shortSynonym;
	char			*S_shortName;
	char			*S_defValue;
	eOptType		E_type;
	xPossibleRange	X_rng;
	union		{
		int			N_intVal;
		xOptRange	X_rngVal;
		wsReal		R_realVal;
		char		*S_strVal;
		xOperation	X_expVal;
	};
};

#define Z(v) (char *)v

#define LOG_BUFFER ""

#define ASSERT_OPTION(a)		OPTION().assertOption(#a)			///< Assert if false
#define ASSERT_OPTS_AND(a1, a2)	OPTION().assertOptionAnd(#a1, #a2)	///< Assert if both false

#define DEF_OPT_RANGE(a)	xOptRange X_##a
#define OPT_RANGE(a)		(OPTION().X_##a)
#define ASG_OPT_RANGE(a)	memcpy(&X_##a, &(Xa_defOption[findOption(#a)].X_rngVal), sizeof(xOptRange))
#define	FORCE_OPT_RANGE(a)	memcpy(&(OPTION().X_##a), &(Xa_defOption[OPTION().findOption(#a)].X_rngVal), sizeof(xOptRange))

#define DEF_OPT_NUMBER(a)	int	N_##a
#define OPT_NUMBER(a)		(OPTION().N_##a)
#define ASG_OPT_NUMBER(a)	N_##a = Xa_defOption[findOption(#a)].N_intVal
#define FORCE_OPT_NUMBER(a)	N_##a = Xa_defOption[OPTION().findOption(#a)].N_intVal

#define OPT_ENABLED(a)		(OPTION().N_##a == 1)

#define DEF_OPT_STRING(a)	char *S_##a
#define OPT_STRING(a)		(OPTION().S_##a)
#define ASG_OPT_STRING(a)	S_##a = Xa_defOption[findOption(#a)].S_strVal
#define FORCE_OPT_STRING(a)	S_##a = Xa_defOption[OPTION().findOption(#a)].S_strVal

#define DEF_OPT_REAL(a)		wsReal R_##a
#define OPT_REAL(a)			(OPTION().R_##a)
#define ASG_OPT_REAL(a)		R_##a = Xa_defOption[findOption(#a)].R_realVal
#define FORCE_OPT_REAL(a)	R_##a = Xa_defOption[OPTION().findOption(#a)].R_realVal

#define DEF_OPT_EXPR(a)		xOperation	O_##a
#define OPT_EXPR(a)			(OPTION().O_##a)
#define ASG_OPT_EXPR(a)		memcpy(&O_##a, &(Xa_defOption[findOption(#a)].X_expVal), sizeof(xOperation))

#define IS_ASSIGNED(a)		OPTION().isAssigned(#a)
#define IS_GIVEN(a)			(OPTION().S_##a != NULL)

typedef	vector<xOption*>	vOptPtr;
typedef	vOptPtr::iterator	vOptPtr_it;

class cOption
{
	vOptPtr	Xv_asgnOpts;
	wsUint	L_defOpt;

	int		__procOption(char *S_opt, char *Sp_prev);
	int		__procValue(int N_idxOpt, const char *S_val, char B_over=0);
	void	_checkOptionIntegrity();
	void	_runAs(const char *S_fn);
	void	_loadConf();
public:
	/* OPTIONS FOR PROGRAM FLOW */
	DEF_OPT_NUMBER(verbose);	// Whether provides more detailed output or not
	DEF_OPT_NUMBER(seed);
	DEF_OPT_STRING(script);
	DEF_OPT_NUMBER(thread);		// Number of threads utilizing
	DEF_OPT_STRING(out);		// Input dataset file name

	/* OPTIONS FOR DATA INPUT */
	DEF_OPT_STRING(lgen);		//	
	DEF_OPT_STRING(tped);
	DEF_OPT_STRING(data);		// Shared data path
	DEF_OPT_STRING(bed);		// Input dataset file name
	DEF_OPT_STRING(fam);		// Alternative family definition file

	/* OPTIONS FOR ADDITIONAL INPUT */
	DEF_OPT_STRING(bim);
	DEF_OPT_STRING(map);
	DEF_OPT_STRING(pheno);	// Alternative phenotype file name
	DEF_OPT_STRING(sampvarflag);
	DEF_OPT_STRING(geneset);
	DEF_OPT_STRING(set);

	/* OPTIONS FOR HANDLING ADDITIONAL INPUT */
	DEF_OPT_NUMBER(ignorefid);
	DEF_OPT_NUMBER(nofid);
	DEF_OPT_NUMBER(noparent);
	DEF_OPT_NUMBER(nosex);
	DEF_OPT_NUMBER(nopos);
	DEF_OPT_NUMBER(nogdist);
	DEF_OPT_STRING(acgt);
	DEF_OPT_NUMBER(nomap);
	DEF_OPT_NUMBER(1234);
	DEF_OPT_STRING(sepallele);
	DEF_OPT_NUMBER(consecallele);
	DEF_OPT_NUMBER(1case);
	DEF_OPT_STRING(mispheno);
	DEF_OPT_STRING(misgeno);

	/* ELEMENTARY ANALYSIS OPTIONS */
	DEF_OPT_STRING(freq);
	DEF_OPT_STRING(hwe);
	DEF_OPT_NUMBER(heritability);

	/* GENE-BASED TESTS */
	DEF_OPT_NUMBER(noweight);
	DEF_OPT_STRING(betaweight);
	DEF_OPT_STRING(weight);
	DEF_OPT_NUMBER(mafweight);
	DEF_OPT_NUMBER(gsmacthr);

	/* VARIANT-LEVEL TESTS */
	DEF_OPT_STRING(prevalence);
	DEF_OPT_NUMBER(scoretest);

	/* SAMPLE RELATEDNESS OPTIONS */
	DEF_OPT_NUMBER(kinship);		// Use PDDT instead of correlation calculation?
	DEF_OPT_NUMBER(ibs);
	DEF_OPT_NUMBER(ktau);
	DEF_OPT_NUMBER(empktau);
	DEF_OPT_NUMBER(corpearson);
	DEF_OPT_NUMBER(cordiag1);
	DEF_OPT_NUMBER(medcor);
	DEF_OPT_NUMBER(indep);
	DEF_OPT_NUMBER(hybrid);
	DEF_OPT_STRING(cor);
	DEF_OPT_NUMBER(corpair);
	DEF_OPT_NUMBER(corgrm);
	DEF_OPT_NUMBER(corepacts);
	DEF_OPT_NUMBER(x);
	DEF_OPT_NUMBER(x2);
	DEF_OPT_RANGE(cormaf);

	/* EXPORT-RELATED OPTIONS */

	DEF_OPT_NUMBER(quiet);
	DEF_OPT_NUMBER(version);
	DEF_OPT_STRING(ped);		// Input dataset file name
	DEF_OPT_STRING(variantvar);	// Alternative phenotype file name
	DEF_OPT_EXPR(filvariant);
	DEF_OPT_EXPR(incvariant);
	DEF_OPT_EXPR(filgeno);
	DEF_OPT_EXPR(incgeno);
	DEF_OPT_STRING(pname);
	DEF_OPT_STRING(remsamp);
	DEF_OPT_STRING(selsamp);
	DEF_OPT_STRING(nasamp);
	DEF_OPT_REAL(randnasamp);
	DEF_OPT_STRING(remvariant);
	DEF_OPT_STRING(selvariant);
	DEF_OPT_STRING(filrange);
	DEF_OPT_STRING(incrange);
	DEF_OPT_REAL(varresize);
	DEF_OPT_NUMBER(varwindow);
	DEF_OPT_REAL(sampresize);
	DEF_OPT_NUMBER(emcount);
	DEF_OPT_REAL(aithr);
	DEF_OPT_NUMBER(autoonly);
	DEF_OPT_NUMBER(sexonly);
	DEF_OPT_STRING(chr);
	DEF_OPT_NUMBER(makeblup);
	DEF_OPT_RANGE(filmaf);
	DEF_OPT_RANGE(filmac);
	DEF_OPT_RANGE(filhwe);
	DEF_OPT_RANGE(filgind);
	DEF_OPT_RANGE(filgvar);
	DEF_OPT_STRING(remfam);
	DEF_OPT_NUMBER(filnf);
	DEF_OPT_NUMBER(filmf);
	DEF_OPT_NUMBER(filcase);
	DEF_OPT_NUMBER(filcontrol);
	DEF_OPT_NUMBER(filmispheno);
// 	DEF_OPT_STRING(filpheno);
// 	DEF_OPT_STRING(filcov);
	DEF_OPT_EXPR(filsample);
	DEF_OPT_RANGE(filgdist);
	DEF_OPT_NUMBER(filnosex);
	DEF_OPT_NUMBER(filmale);
	DEF_OPT_NUMBER(filfemale);
	DEF_OPT_NUMBER(snvonly);
	DEF_OPT_NUMBER(indelonly);
	DEF_OPT_RANGE(incmaf);
	DEF_OPT_RANGE(incmac);
	DEF_OPT_RANGE(inchwe);
	DEF_OPT_RANGE(incgind);
	DEF_OPT_RANGE(incgvar);
	DEF_OPT_STRING(selfam);
	DEF_OPT_RANGE(incgdist);
// 	DEF_OPT_STRING(incpheno);
// 	DEF_OPT_STRING(inccov);
	DEF_OPT_EXPR(incsample);
	DEF_OPT_NUMBER(ml);
	DEF_OPT_STRING(blup);
	DEF_OPT_NUMBER(imputepheno);
	DEF_OPT_STRING(cname);
	DEF_OPT_STRING(makecor);
	DEF_OPT_NUMBER(empiall);
	DEF_OPT_STRING(cact);
	DEF_OPT_NUMBER(1sex);
	DEF_OPT_STRING(mafe);
	DEF_OPT_NUMBER(ginv);
	DEF_OPT_NUMBER(makecov);
	DEF_OPT_NUMBER(makepheno);
	DEF_OPT_STRING(makeflag);
	DEF_OPT_NUMBER(founderonly);
	DEF_OPT_NUMBER(famsummary);
	DEF_OPT_NUMBER(ncsummary);
	DEF_OPT_NUMBER(specdcmp);
	DEF_OPT_NUMBER(nostop);
	DEF_OPT_NUMBER(ignoreparent);
	DEF_OPT_STRING(sepid);
	DEF_OPT_NUMBER(nopheno);
	DEF_OPT_STRING(probandcol);
	DEF_OPT_STRING(twincol);
	DEF_OPT_NUMBER(chrwise);
	DEF_OPT_STRING(est);
	DEF_OPT_NUMBER(makeclgeno);
	DEF_OPT_NUMBER(indel);
	DEF_OPT_NUMBER(logistic);
	DEF_OPT_STRING(sortvariant);
	DEF_OPT_STRING(sortpos);
	DEF_OPT_STRING(sortsample);
	DEF_OPT_STRING(sortiid);
	DEF_OPT_STRING(baseline);
	DEF_OPT_NUMBER(gz);
	DEF_OPT_NUMBER(makenrm);
	DEF_OPT_NUMBER(makeweight);
	DEF_OPT_NUMBER(nophenohdr);
	DEF_OPT_RANGE(pvalrange);
	DEF_OPT_NUMBER(explore);
	DEF_OPT_NUMBER(time);
	DEF_OPT_NUMBER(phenostdize);
	DEF_OPT_STRING(outmispheno);
	DEF_OPT_STRING(outmisgeno);
	DEF_OPT_STRING(gxecovs);
	DEF_OPT_NUMBER(passemptyline);
	DEF_OPT_NUMBER(listvariant);
	DEF_OPT_NUMBER(listsample);
	DEF_OPT_NUMBER(listfounder);
	DEF_OPT_NUMBER(nolmm);
	DEF_OPT_NUMBER(gmdr);
	DEF_OPT_RANGE(gsetconsec);
	DEF_OPT_NUMBER(nosysfeat);
	DEF_OPT_STRING(fname);
	DEF_OPT_NUMBER(makeev);
	DEF_OPT_STRING(ev);
	DEF_OPT_NUMBER(natural);
	DEF_OPT_NUMBER(dupnaming);
	DEF_OPT_NUMBER(citation);
	DEF_OPT_NUMBER(miss);
	DEF_OPT_STRING(misparent);
	DEF_OPT_NUMBER(sampmajor);
	DEF_OPT_NUMBER(nospecdcmp);
	DEF_OPT_STRING(species);
	DEF_OPT_NUMBER(nodata);
	DEF_OPT_NUMBER(regex);
	DEF_OPT_STRING(model);

//#if (TOOLSET_TYPE == TOOLSET_MFQLS) || ((TOOLSET_TYPE & 0x100) == 0x100)
	DEF_OPT_NUMBER(mqls);
	DEF_OPT_NUMBER(fqls);
	DEF_OPT_STRING(heri);
	DEF_OPT_NUMBER(multifqls);
	DEF_OPT_NUMBER(mqlsconsec);

	DEF_OPT_NUMBER(fqlsnopddt);
	DEF_OPT_REAL(retestthr);
	DEF_OPT_NUMBER(avail);
	DEF_OPT_NUMBER(fastfqls);
	DEF_OPT_NUMBER(fastmqls);
	/* #endif
 
#if (TOOLSET_TYPE == TOOLSET_FARVAT) || (TOOLSET_TYPE == TOOLSET_MFARVAT) || \
	((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE == TOOLSET_FARVATX)*/
	DEF_OPT_NUMBER(genesummary);
	DEF_OPT_NUMBER(genetest);
	DEF_OPT_REAL(raremaf);
	DEF_OPT_RANGE(genesize);
	DEF_OPT_NUMBER(pedcmc);
	DEF_OPT_NUMBER(wsum);
	DEF_OPT_NUMBER(kbac);
	DEF_OPT_NUMBER(asum);
	DEF_OPT_NUMBER(farvat);
	DEF_OPT_NUMBER(pedgene);
	DEF_OPT_NUMBER(skato);
	DEF_OPT_NUMBER(skat);
	DEF_OPT_NUMBER(farvatx);
	DEF_OPT_REAL(farvatxd);

	DEF_OPT_NUMBER(gmapsummary);
	DEF_OPT_REAL(genemiss);
	DEF_OPT_NUMBER(skatondiv);
	DEF_OPT_STRING(skatodivs);
	DEF_OPT_REAL(kbacalpha);
	DEF_OPT_NUMBER(kbac2side);
	DEF_OPT_STRING(kbackernel);
	DEF_OPT_NUMBER(mfhom);
	DEF_OPT_NUMBER(mfhet);
	DEF_OPT_NUMBER(makegeno);
	DEF_OPT_NUMBER(adjf1);
	DEF_OPT_NUMBER(adjf2);
	DEF_OPT_NUMBER(makefarvat);
	DEF_OPT_NUMBER(farvatxndiv);
// #endif
// 
// #if (TOOLSET_TYPE == TOOLSET_MFQLS) || ((TOOLSET_TYPE & 0x100) == 0x100)
	DEF_OPT_NUMBER(qtest);
	DEF_OPT_RANGE(qtestrange);
	DEF_OPT_NUMBER(qtestclump);
	DEF_OPT_REAL(qteststt);
	DEF_OPT_NUMBER(qtestbetacov);
	DEF_OPT_NUMBER(makebeta);
// #endif
// 
// #if (TOOLSET_TYPE == TOOLSET_PHARAOH) || ((TOOLSET_TYPE & 0x100) == 0x100)
	DEF_OPT_NUMBER(pharaoh);
	DEF_OPT_NUMBER(proopt);
	DEF_OPT_STRING(prolambda);
	DEF_OPT_RANGE(prorange);
	DEF_OPT_REAL(prothr);
	DEF_OPT_NUMBER(promaxiter);
	DEF_OPT_RANGE(progenesize);
	DEF_OPT_RANGE(progsetsize);
	DEF_OPT_NUMBER(prosingle);
	DEF_OPT_NUMBER(propermcov);
	DEF_OPT_NUMBER(nperm);
	DEF_OPT_STRING(seqperm);
	DEF_OPT_STRING(permfile);
	DEF_OPT_NUMBER(cv);
	DEF_OPT_NUMBER(gesca);
	DEF_OPT_STRING(modeltype);
	DEF_OPT_STRING(ggpath);
// #endif

//#if TOOLSET_TYPE == TOOLSET_ONETOOL
	DEF_OPT_NUMBER(debug);
	DEF_OPT_NUMBER(help);
//#endif

// #if ((TOOLSET_TYPE & 0x100) == 0x100)
	DEF_OPT_STRING(dosage);
	DEF_OPT_STRING(expression);
	DEF_OPT_STRING(genoprob);
	DEF_OPT_STRING(vcf);		// VCF dataset file name

	/* ELEMENTARY ANALYSIS OPTIONS */
	DEF_OPT_STRING(fst);
	DEF_OPT_NUMBER(mendel);
	DEF_OPT_NUMBER(famuniq);
	DEF_OPT_NUMBER(pca);		// Use PCA correction?
	DEF_OPT_NUMBER(npc);		// Number of PCs to extract

	/* GENE-BASED TESTS */
	DEF_OPT_NUMBER(genesplit);
	DEF_OPT_NUMBER(rvtdt);
	DEF_OPT_NUMBER(fbskat);
	DEF_OPT_NUMBER(famvt);
	DEF_OPT_NUMBER(ggemma);
	DEF_OPT_STRING(longitudinal);
	DEF_OPT_NUMBER(vt);

	/* VARIANT-LEVEL TESTS */
	DEF_OPT_NUMBER(fisher);
	DEF_OPT_NUMBER(trend);
	DEF_OPT_NUMBER(regression);
	DEF_OPT_NUMBER(qls);
	DEF_OPT_NUMBER(tdt);
	DEF_OPT_NUMBER(sdt);
	DEF_OPT_NUMBER(emmax);
	DEF_OPT_NUMBER(gemma);

	/* SAMPLE RELATEDNESS OPTIONS */

	/* EXPORT-RELATED OPTIONS */
	DEF_OPT_NUMBER(makeped);
	DEF_OPT_NUMBER(maketped);
	DEF_OPT_NUMBER(makebed);
	DEF_OPT_NUMBER(makeraw);
	DEF_OPT_NUMBER(makedom);
	DEF_OPT_NUMBER(makerec);
	DEF_OPT_NUMBER(makevcf);
	DEF_OPT_NUMBER(makebcf);
	DEF_OPT_NUMBER(makelgen);
	DEF_OPT_NUMBER(makegen);
	DEF_OPT_NUMBER(makebgen);
	DEF_OPT_NUMBER(makebeagle);

	/* S.A.G.E. OPTIONS */
	DEF_OPT_NUMBER(relpair);
	DEF_OPT_NUMBER(fcor);
	DEF_OPT_NUMBER(segreg);
	DEF_OPT_STRING(par);
	DEF_OPT_NUMBER(lodlink);
	DEF_OPT_STRING(typ);
	DEF_OPT_NUMBER(fcorStdErrOff);
	DEF_OPT_NUMBER(lodlinkLinkageTestOff);
	DEF_OPT_NUMBER(lodlinkLinkageHomogOff);
	DEF_OPT_NUMBER(lodlinkLinkageSexSpecific);
	DEF_OPT_NUMBER(lodlinkSmithHomogTest);
	DEF_OPT_NUMBER(lodlinkGenotypes);

	DEF_OPT_NUMBER(merlin);
	DEF_OPT_NUMBER(plot);

	DEF_OPT_NUMBER(simtrio);	// Use trio simulation instead of data input?
	DEF_OPT_NUMBER(szfam);		// # of families in the simulation
	DEF_OPT_NUMBER(szvar);		// Use trio simulation instead of data input?
	DEF_OPT_NUMBER(simfam);
	DEF_OPT_NUMBER(trio);
	DEF_OPT_NUMBER(extfam);
	DEF_OPT_NUMBER(nsig);
	DEF_OPT_REAL(sigmaf);
	DEF_OPT_STRING(simfreq);
	DEF_OPT_STRING(mafvar);
	DEF_OPT_REAL(proppc);
	DEF_OPT_NUMBER(usemf);		// Use pseudoparent's information?
	DEF_OPT_NUMBER(fullpca);	// Use full eigenreduction in PCA?
	DEF_OPT_STRING(rpath);		// A path to R
	DEF_OPT_STRING(grm);
	DEF_OPT_REAL(grmalpha);
	DEF_OPT_NUMBER(makemdr);
	DEF_OPT_NUMBER(zipbgen);
	DEF_OPT_RANGE(filqual); 
	DEF_OPT_RANGE(incqual);
	DEF_OPT_RANGE(filmendelfam);
	DEF_OPT_RANGE(incmendelfam);
	DEF_OPT_RANGE(filmendelsamp);
	DEF_OPT_RANGE(incmendelsamp);
	DEF_OPT_RANGE(filmendelvar);
	DEF_OPT_RANGE(incmendelvar);
	DEF_OPT_NUMBER(lrt);
	DEF_OPT_NUMBER(vcfqc);
	DEF_OPT_NUMBER(phasedonly);
	DEF_OPT_NUMBER(unphasedonly);
	DEF_OPT_NUMBER(interactive);
	DEF_OPT_NUMBER(mdr);
	DEF_OPT_NUMBER(order);
	DEF_OPT_NUMBER(top);
	DEF_OPT_NUMBER(hmdr);
	DEF_OPT_STRING(hmdrall);
	DEF_OPT_STRING(hmdrprior);
	DEF_OPT_NUMBER(check);
	DEF_OPT_NUMBER(gxe);
	DEF_OPT_REAL(beta);
	DEF_OPT_REAL(rho);
	DEF_OPT_REAL(rhopheno);
	DEF_OPT_RANGE(nsamp);
	DEF_OPT_NUMBER(npheno);
	DEF_OPT_NUMBER(noshuffle);
	DEF_OPT_NUMBER(shuffle);
	DEF_OPT_NUMBER(powercalc);
	DEF_OPT_NUMBER(powercalc2);
	DEF_OPT_NUMBER(split);
	DEF_OPT_STRING(merge);
	DEF_OPT_NUMBER(mergemode);
	DEF_OPT_NUMBER(testmatrix);
	DEF_OPT_NUMBER(testmatfunc);
	DEF_OPT_NUMBER(testmatclass);
	DEF_OPT_NUMBER(testfunc);
	DEF_OPT_NUMBER(ld);
	DEF_OPT_NUMBER(ldcor);
	DEF_OPT_STRING(annogene);
	DEF_OPT_NUMBER(annorange);
	DEF_OPT_STRING(annovar);
	DEF_OPT_NUMBER(density);
	DEF_OPT_NUMBER(tstv);
	DEF_OPT_STRING(sep);
	DEF_OPT_NUMBER(dsgdist);
	DEF_OPT_NUMBER(pc2cov);
	DEF_OPT_NUMBER(donull);
	DEF_OPT_STRING(updvariant);
	DEF_OPT_STRING(updchr);
	DEF_OPT_STRING(updname);
	DEF_OPT_STRING(updgdist);
	DEF_OPT_STRING(updpos);
	DEF_OPT_STRING(updgeno);
	DEF_OPT_STRING(ref);
	DEF_OPT_NUMBER(ldbin);
	DEF_OPT_NUMBER(ldsize);
	DEF_OPT_STRING(ldvar);
	DEF_OPT_NUMBER(remna);
	DEF_OPT_NUMBER(boost);
	DEF_OPT_REAL(thrboost);
	DEF_OPT_NUMBER(quickepi);
	DEF_OPT_STRING(ext);
	DEF_OPT_NUMBER(impute);
	DEF_OPT_NUMBER(lod);
	DEF_OPT_NUMBER(makeimpute);
	DEF_OPT_NUMBER(sxa);
	DEF_OPT_STRING(R);
	DEF_OPT_NUMBER(randbinpheno);
	DEF_OPT_NUMBER(randpheno);
	DEF_OPT_NUMBER(genoctrl);
	DEF_OPT_REAL(usergc);
	DEF_OPT_NUMBER(adjust);
	DEF_OPT_STRING(popuniq);
	DEF_OPT_NUMBER(monotone);
	DEF_OPT_NUMBER(singleton);
	DEF_OPT_NUMBER(doubleton);
	DEF_OPT_NUMBER(genemdr);
	DEF_OPT_STRING(variant2cov);
	DEF_OPT_NUMBER(inbreed);
	DEF_OPT_STRING(group);
	DEF_OPT_NUMBER(filgenic);
	DEF_OPT_NUMBER(filintergenic);
	DEF_OPT_NUMBER(variantsummary);
	DEF_OPT_STRING(sampleorder);
	DEF_OPT_STRING(variantorder);
	DEF_OPT_STRING(cosi);
	DEF_OPT_NUMBER(setconsec);
	DEF_OPT_NUMBER(setoverlap);
	DEF_OPT_RANGE(setrandom);
	DEF_OPT_NUMBER(makeset);
	DEF_OPT_STRING(outcact);
	DEF_OPT_NUMBER(out1case);
	DEF_OPT_NUMBER(settype);
	DEF_OPT_NUMBER(mistest);
	DEF_OPT_RANGE(incmistest);
	DEF_OPT_RANGE(filmistest);
	DEF_OPT_REAL(nageno);
	DEF_OPT_REAL(napheno);
	DEF_OPT_NUMBER(filtreport);
	DEF_OPT_STRING(genofield);	/* FIXME : Impl needed */
	DEF_OPT_NUMBER(outphenoonly);
	DEF_OPT_NUMBER(outnoheader);
	DEF_OPT_NUMBER(makemerlin);
	DEF_OPT_NUMBER(mds);
	DEF_OPT_NUMBER(gxg);
	DEF_OPT_NUMBER(invnorm);
	DEF_OPT_NUMBER(forceconv);
	DEF_OPT_NUMBER(out1234);
	DEF_OPT_STRING(outacgt);
	DEF_OPT_NUMBER(dfam);	/* FIXME : Impl needed */
	DEF_OPT_STRING(updallele);
	DEF_OPT_RANGE(window);
	DEF_OPT_REAL(ci);
	DEF_OPT_NUMBER(lasso);
	DEF_OPT_REAL(lassolambda);
	DEF_OPT_NUMBER(lassoall);
	DEF_OPT_NUMBER(pls);
	DEF_OPT_STRING(sampleweight);
	DEF_OPT_NUMBER(nskip);
	DEF_OPT_NUMBER(singleparent);
	DEF_OPT_STRING(bcf);
	DEF_OPT_NUMBER(loocv);
	DEF_OPT_RANGE(mdrthr);
	DEF_OPT_NUMBER(ldcontrast);
	DEF_OPT_STRING(flip);
	DEF_OPT_STRING(varsubset);
	DEF_OPT_NUMBER(hethom);
	DEF_OPT_NUMBER(markercheck);
	DEF_OPT_STRING(meta);
	DEF_OPT_STRING(fid);
	DEF_OPT_STRING(outformat);
	DEF_OPT_STRING(maf);
	DEF_OPT_NUMBER(het);
	DEF_OPT_REAL(variantblup);
	DEF_OPT_NUMBER(famsplit);
	DEF_OPT_NUMBER(setspan);

	DEF_OPT_NUMBER(tridge);
	DEF_OPT_NUMBER(hamming);
	DEF_OPT_NUMBER(bn);

	DEF_OPT_NUMBER(fuzzymdr);
	DEF_OPT_NUMBER(gxgall);
	DEF_OPT_STRING(gxglist);
	DEF_OPT_STRING(gxglambda);

	DEF_OPT_STRING(prunevif);
	DEF_OPT_STRING(prunepw);
//#endif
	//DEF_OPT_NUMBER(corbn);

	/*** STEP 2 : Add the definition of new option here */

	DEF_OPT_NUMBER(phenoCase);		// Not assignable
	DEF_OPT_NUMBER(phenoCtrl);		// Not assignable
	DEF_OPT_NUMBER(sexMale);		// Not assignable
	DEF_OPT_NUMBER(sexFema);		// Not assignable
	DEF_OPT_NUMBER(maxNumChr);		// Not assignable
	DEF_OPT_NUMBER(maxNumAutoChr);	// Not assignable
	DEF_OPT_STRING(nonAutoChrSeq);	// Not assignable
	DEF_OPT_STRING(cmdLine);		// Not assignable

	cOption();
	cOption(const int N_argc, char *Sa_argv[]);
	~cOption();
	void		clear();
	int			init(const int N_argc, char *Sa_argv[]);
	wsUint		findOption(const char *S_longName, char B_noError=0);
	char		isAssigned(const char *S_longName);
	void		setUnassigned(const char *S_longName);
	void		assertOption(const char *S_longName);
	void		assertOptionAnd(const char *S_ln1, const char *S_ln2);
	void		assign(const char *S_name, const char *S_val=NULL,
		char B_over=0);
	const char*	getDefVal(const char *S_longName);
	int			procRange(xOption *Xp_opt, const char *S_val);
	vOptPtr&	getAsgnOpts() { return Xv_asgnOpts; }
};

cOption& OPTION();

extern xOption Xa_defOption[];

char isInRange(xOptRange& X_range, wsReal R_val);
char isInRangeOR(xOptRange& X_range, wsReal* Ra_val, wsUint N_val);
char isInRangeAND(xOptRange& X_range, wsReal* Ra_val, wsUint N_val);

} // End namespace ONETOOL

#endif
