#pragma once
#ifndef __WISARD_IO_H__
#define __WISARD_IO_H__

#include "global/common.h"
#include "utils/util.h"
#include "utils/vector.h"
#include "global/option.h"
#include "utils/regex.h"
#include "utils/genet.h"
#include "input/sampvar.h"

namespace ONETOOL {

class cVector;
int forAllVariant_equal(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);

#define	MFR_MENDELIAN	10	/* Filtered by Mendelian error */
#define	MFR_MISCACT		11	/* Filtered by missingness test across ca/ct */
#define MFR_RGNGENE		12	/* Filtered by gene-region */
#define	MFR_HWE			13	/* Filtered by Hardy-Weinberg Equilibrium */
#define	MFR_MAC			14	/* Filtered by Minor Allele Count */
#define MFR_GTRATE		15	/* Filtered by genotyping rate */
#define	MFR_MAF			16	/* Filtered by Minor Allele Frequency */
#define MFR_TRIALLELIC	17	/* Filtered by Tri-allelcity */
#define MFR_WINDOW		18	/* Filtered by --window */

typedef enum _xFileType
{
	FT_NONE,
	FT_PED,
	FT_BED,
	FT_RAW,
	FT_LGEN,
	FT_TPED,
	FT_VCF,
	FT_BCF,
	FT_DOSAGE,
	FT_EXPR,
	FT_SIM,
	FT_MACH,
	FT_GEN,
	FT_BGEN,
	FT_BEAGLE,
	FT_META,
} xFileType;

typedef struct _xUpdAlleleInfo {
	char o1, o2;
	char u1, u2;
	char *O1, *O2;
	char *U1, *U2;
} xUpdAlleleInfo;

typedef map<string,xUpdAlleleInfo>	mUAinfo;
typedef mUAinfo::iterator			mUAinfo_it;

class cAnalysis;
class cSetManagerAnalysis;
wsMat getFullCorMat(cAnalysis *Cp_anaPhi);

struct xVariant
{
	int		chr;
	char	*name;
	char	al1, al2;
	char	filter;
	char	*indel1, *indel2;
	char	*anno;
	wsUint	pos, nmis;
	wsReal	hwe, gdist;
};

typedef enum _xDataFiltStrategy {
	FILT_NOTHING,
	FILT_REMOVE,
	FILT_SELECT
} xDataFiltStrategy;

extern const wsUint MAX_TWINS;						// Maximum number of twins acceptable in wisard
extern const wsUint WINSARD_WARNSIZE_EMPCORCALC;	// Minimum required number of markers to calculate empirical cor. correctly
extern const wsUint PED_MAXLEN;

#define isSexChromosome(d)	((d).chr > (int)NAUTO_SPECIES)
#define isXChromosome(d)	((d).chr == (int)NAUTO_SPECIES+1)
#define isYChromosome(d)	(getChrName2((d).chr)[0] == 'Y')
#define isMtChromosome(d)	(getChrName2((d).chr)[0] == 'M')
#define isAutosome(d)		((d).chr <= OPT_NUMBER(maxNumAutoChr) && (d).chr > 0)
//#define PED_MAXLEN			1024*1024*64
#define SAMP_NODATA			-1

#define VARIANT_SSE_PART(vrt, v, u) wsUint N_byte = (vrt+3)>>2; \
	wsUint N_med = (N_byte>>4)<<4; \
	for (wsUint X=0 ; X<N_med ; X+=sseJmp*sizeof(wsReal)) { \
		__m128i* u = (__m128i *)(v + X);
#define VARIANT_NRM_PART(vrt, v, u) } for (wsUint X=N_med ; X<(N_byte-1) ; X++) { \
		unsigned char u = v[X];
#define VARIANT_END_PART(vrt, v, u) } unsigned char u = v[N_byte-1]; \
		unsigned int N_pureByte = (vrt>>2)<<2; \
		if (N_pureByte != vrt) \
			u &= ~((0xff)<<((vrt-N_pureByte)<<1)); \

struct xSample
{
	string				S_FID, S_IID;		///< Sample family ID and individual ID
	xSample				*Xp_pat, *Xp_mat;	///< Pointer to parents, both NULL if founder, one-side NULL is not permitted
	vector<xSample*>	Xp_childs;			///< Pointer vector of childs
	vector<xSample*>	Xp_spouses;			///< Pointer vector of "spouses"
	char				N_sex;				///< Gender, 1=male, 2=female, -9=undetermined, other values are not permitted
	char				B_isAssigned;		///< Is the record for this sample found
	char				B_isMissing;		///< Is this sample "missing its genotypes"
	char				B_isProband;		///< Is this sample "proband"
	char				N_idTwin;			///< Index of "twin"
	char				B_isComplete;		///< Is this sample have no missing
	short				N_isVisited;		///< For 
	int					N_missGeno;			///< # missing genotype
	int					N_idx, N_oriIdx;	///< Index in getSamples() and index of input file
	int					N_grpFst;			///< Unused
	vector<xSample*>	Xp_twins;			///< Pointer vector of twins
};

typedef map<string,xSample*>	mSampPtr;
typedef mSampPtr::iterator		mSampPtr_it;
typedef map<string,xSample>		mSamp;
typedef mSamp::iterator			mSamp_it;
typedef vector<xSample*>		vSampPtr;
typedef vSampPtr::iterator		vSampPtr_it;
typedef vector<xVariant>		vVariant;
typedef vVariant::iterator		vVariant_it;
typedef	vector<xVariant*>		vVrtPtr;
typedef	vVrtPtr::iterator		vVrtPtr_it;

struct xFamily
{
	wsUint		N_size;			///< # member
	string		S_FID;			///< Family identifier
	vSampPtr	Xp_founders;	///< xSample vector of founders
	void setUnvisited()
	{
		mSampPtr Xp_currVisitor;

		/* Set initial visitors to founders */
		FOREACH (vSampPtr_it, Xp_founders, it)
			Xp_currVisitor.insert(make_pair((*it)->S_IID, *it));

		/* While there are remained visitors */
		while (Xp_currVisitor.size()) {
			mSampPtr Xp_newVisitors;

			/* Visiting all visitors */
			FOREACH (mSampPtr_it, Xp_currVisitor, it) {
				xSample *Xp_currSamp = it->second;
				Xp_currSamp->N_isVisited = 0;

				/* Add all childs to new visitors */
				FOREACH (vSampPtr_it, Xp_currSamp->Xp_childs, sit) {
					/* Insert sample if B_isVisited is not zero */
					if ((*sit)->N_isVisited != 0)
						Xp_newVisitors.insert(make_pair((*sit)->S_IID, *sit));
				}
			}

			Xp_currVisitor.clear();
			Xp_currVisitor = Xp_newVisitors;
		}
	}
	void visit(int (*H_func)(xSample*, void *), void *Vp_data)
	{
		mSampPtr Xp_currVisitor;

		setUnvisited();

		/* Set initial visitors to founders */
		FOREACH (vSampPtr_it, Xp_founders, it) {
			Xp_currVisitor.insert(make_pair((*it)->S_IID, *it));
//			H_func(*it, Vp_data);
		}

		/* While there are remained visitors */
		while (Xp_currVisitor.size()) {
			mSampPtr Xp_newVisitors;

			/* Visiting all visitors */
			FOREACH (mSampPtr_it, Xp_currVisitor, it) {
				xSample *Xp_currSamp = it->second;
				if (Xp_currSamp->N_isVisited == 1)
					continue;

				Xp_currSamp->N_isVisited = 1;
				H_func(Xp_currSamp, Vp_data);

				/* Add all childs to new visitors */
				FOREACH (vSampPtr_it, Xp_currSamp->Xp_childs, sit) {
					/* Insert sample if B_isVisited is not zero */
					if ((*sit)->N_isVisited != 1)
						Xp_newVisitors.insert(make_pair((*sit)->S_IID, *sit));
				}
			}

			Xp_currVisitor.clear();
			Xp_currVisitor = Xp_newVisitors;
		}
	}
	void visitWithLevel(int (*H_func)(xSample*, int, void *), void *Vp_data)
	{
		mSampPtr Xp_currVisitor;

		setUnvisited();

		/* Set initial visitors to founders */
		FOREACH (vSampPtr_it, Xp_founders, it) {
			Xp_currVisitor.insert(make_pair((*it)->S_IID, *it));
			//			H_func(*it, Vp_data);
		}

		/* While there are remained visitors */
		for (wsUint L=0 ; Xp_currVisitor.size() ; L++) {
			mSampPtr Xp_newVisitors;

			/* Visiting all visitors */
			FOREACH (mSampPtr_it, Xp_currVisitor, it) {
				xSample *Xp_currSamp = it->second;
				if (Xp_currSamp->N_isVisited == 1)
					continue;

				Xp_currSamp->N_isVisited = 1;
				H_func(Xp_currSamp, L, Vp_data);

				/* Add all childs to new visitors */
				FOREACH (vSampPtr_it, Xp_currSamp->Xp_childs, sit) {
					/* Insert sample if B_isVisited is not zero */
					if ((*sit)->N_isVisited != 1)
						Xp_newVisitors.insert(make_pair((*sit)->S_IID, *sit));
				}
			}

			Xp_currVisitor.clear();
			Xp_currVisitor = Xp_newVisitors;
		}
	}
	vInt	Xv_members;			///< Indices of members in getSamples()
};

typedef map<string,xFamily>		mFam;
typedef mFam::iterator			mFam_it;

typedef enum _xModel {
	MD_ADDITIVE,
	MD_DOMINANT,
	MD_RECESSIVE,
	MD_MULTIPLICATIVE,
	MD_ERROR
} xModel;

typedef struct _xMaf {
	wsReal	R_maf, R_allMaf;
	wsUint	N_allele, N_allMac, N_mac;
	bool	B_flipped;
} xMaf;

int forAllVariant(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);
int forAllSample(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue);

/* Null definition */
class cStream;
class cPPPAnalysisV2;
class cImpute;

class cIO
{
	cPPPAnalysisV2 *Cp_anaPPP;
protected:
/* N_sampOrig : Number of samples ACTUALLY exists in the input
 * Initialization : MANUAL
 * Setting        : MANUAL */
	wsUint		N_sampOrig;
/* N_vrtOrig : Number of variants ACTUALLY exists in the input
 * Initialization : MANUAL
 * Setting
 *		BED  (CONFIRMED)
 */
/* N_sample : Number of samples LOADED to memory
 * Initialization : MANUAL
 * Setting        : _loadRealTimeFamData */
	wsUint		N_sample;
	wsUint		N_vrtOrig;
	wsUint		N_variant;

/* N_missGeno : Indicates the number of missing genotypes across ENTIRE file
 * Initialization : MANUAL
 * Setting        :
 *		PED  (CONFIRMED)
 *		BED  (CONFIRMED)
 *      RAW */
	wsUint		N_missGeno;

/* Ba_chrAllowed : Indicates that (i)th chromosomes should be excluded[0] or not[1] (RESIZE_SAFE)
 * Type : Array[size MAX_CHR]
 * Initialization : cIO()
 * Setting        : cIO() */
	char		Ba_chrAllowed[MAX_NCHR];

	eStr		Xe_varSubset;

/* Ba_chrExists : Indicates that (i)th chromosomes exist[0] or not[1] in dataset (RESIZE_SAFE)
 * Type : Array[size MAX_CHR + 1] (for Un)
 * Initialization : cIO()
 * Setting        : _initFilterings() */
	char		Ba_chrExists[MAX_NCHR+1];
/* Whether alternative marker def. is available[1] or not[0]
 * Type : Scalar
 * Initialization : cIO()
 * Setting        : _initFilterings() */
	char		B_altVariantDef;

/* Whether this instance have actual data[0] or not[1]
 * Type : Scalar
 * Initialization : cIO()
 * Setting        : cIO() with initial assignment */
	char		B_dry;


/* Ba_isFounder : Indicates that (i)th samples is founder[1](mat=0 && pat=0) or not[0] (RESIZE_SAFE)
 * Type : Array[size #samp]
 * Initialization : _initMemories
 * Setting        :
 *		BED  (CONFIRMED)
 *		PED  (CONFIRMED)
 *		RAW  (IMPOSSIBLE)
 *		TPED (CONFIRMED)
 *		VCF  (CONFIRMED)
 *		DSG
 * */
	char		*Ba_isFounder;
/* RESIZE_SAFE
 * Indicates the number of parents in dataset
 * Automatically determined
 * Initialization :
 * Setting        :
 *		BED  (CONFIRMED)
 *		PED  (CONFIRMED)
 *		RAW  (IMPOSSIBLE)
 *		TPED (CONFIRMED)
 *		VCF  (CONFIRMED)
 *		DSG  (CONFIRMED)
 * */
	wsUint		N_founder;

/* RESIZE_SAFE
 * Matrix[size #pheno*#samp]
 */
	wsMat		Ra_pheno;
	wsUint		N_pheno;

/* Array[size #pheno]
 * Indicates that (i)th phenotype is continuous[1] or dichotomous[0]
 * Initialization :
 * Setting        :
 *		BED  (CONFIRMED)
 *		PED  (CONFIRMED)
 *		RAW  (IMPOSSIBLE)
 *		TPED (CONFIRMED)
 *		VCF
 *      DSG
 * */
	char		*Ba_typePheno;
	wsUint		N_cov, N_covGxE;
	wsUint		N_twin;
	char**		Na_geno;
	unsigned char**
				Na_raw;
	wsFloat**	Ra_data;
	wsFloat**	Ra_dosage;
	wsUint**	Na_gProb;
	unsigned char**
				Na_phased1;
	unsigned char**
				Na_phased2;

/* Matrix[size #samp*#vrt]
 * Used only for loading PED file,
 * Contains original PED genotypes into coding */
	USHORT_t	**Na_charData;

/* Matrix[size #cov*#samp]
 * Contains covariates used in analysis */
	wsMat		Ra_cov;
/* Matrix[size #gxecov*#samp]
 * Contains covariates for GxE used in analysis
 * Note that it is just a subset of pointers from Ra_cov */
	wsMat		Ra_covGxE;

/* Indicates whether the given dataset have
 * no missing(1) or not(0)
 * Setting : cIO::_finalize() */
	char		B_isDataComplete;
	vVariant	Xv_variant;
	vVariant*	Xv_altVariantDef;
/* Array[size #vrt]
 * Indicates whether (i)th variant should be incorporated into correlation est.
 * incorporate(filled with 0xff) or not (filled with 0x00) */
//	wsReal		*Ra_corrMask;
/* Array[size #vrt]
 * Includes the weights of variants that customized from file input,
 * should be NULL unless --weight is given */
	wsVec		Ra_customWeight;
/* Array[size #samp]
 * Have pointers to corresponding xSample pointer
 * by its index */
	vSampPtr	Xa_sampleV2;
	mFam		Xm_fid2data;

/* Indicate that data have proband information or not */
	char		B_haveProband;

/* Sample group index to actual group name */
	vStr		Sv_sampGrpIdx2grp;

/* Indicates the mode of filtering sample
 * 1 - Specified sample will be removed
 * 2 - Specified sample will be selected */
	char		N_filterSample;
/* Indicates the mode of filtering variant
 * 1 - Specified variant will be removed
 * 2 - Specified variant will be selected */
	char		N_filterVrt;
/* Indicates the mode of filtering variants by range
 * 1 - variants in specified range will be removed
 * 2 - variants in specified range will be selected */
	char		N_filterRange;
	vBool		Bv_filtSamp;

/* Alternative A/C/G/T coding definition,
 * could be ascii character ranged from 0x00 ~ 0x7f
 * otherwise should be error */
	char*		Na_ACGT;

	wsUint		N_szFilteredSamp, N_szFilteredVrt;
/* MAF information */
	xMaf*		Xp_maf;
	char		B_mafComputed;

/* Array[size #varFile]
 * Indicates that (i)th variant in the file should be filtered out[1] or not[0] */
	vBool		Bv_filtVrt;
	vBool		Bv_filtVrt1pass;
	
	vRange*		Xv_listRanges;
	/* w/o --regex */
	eStr		Xe_listVrt;
	eStr		Xe_listSample;
	eStr		Xe_listNAsamp;
	eStr		Xe_listFam;
	/* w/ --regex */
	vPattern	Xv_listVrt;
	vPattern	Xv_listSample;
	vPattern	Xv_listFam;

	wsUint		N_altSample;

	vPheno		V_phenos;
	vVarRng		V_filtPhenos;
	//vStr		Sa_phenoName;
	vCovar		Xa_covInfo;
	vVarRng		V_filtCovs;

	cImpute*	Cp_impute; ///< An imputation instance

	vInt		Xa_gxeIdx;
	char		B_blupApplied;

	mStrReal	Xm_newMAF;

/* CNV */
	wsMat		Ra_cnv;
	wsUint		N_szCnv;

	cVector		V_sampWgt;	/* [#samp] */

	void		_loadSampleList(char B_select=0);
	void		_loadSampleSelectionList();
	void		_loadSampleNAList();
	void		_loadFamilyList(char B_select=0);
	void		_loadMarkerAltInfo(wsStrCst S_fn, int N_mode);
	vVariant*	_loadVariantDef(wsStrCst S_fn);
	void		_loadVariantRemovalList();
	void		_loadVariantSelectionList();
	void		_loadRangeRemovalList();
	void		_loadRangeSelectionList();
	void		_loadVariantVar(vStr &Xv_reqVars);
	void		_filtVar(wsStrCst Sp_filExpr, char B_filt, char B_forCov);
	void		_loadBlup();
	void		_loadVariantWeight();
	void		_md();
	void		_freq();
	void		_hwe();	/* Do Hardy-Weinberg Equilibrium test */
	void		_mistest();
	void		_buildTwin();
	void		_resizeData();
	xModel		_getModel();
protected:
	mSamp		Xm_iid2data;
	void		_initFilterings();
	void		_initMemories(char B_noGeno=0);
	xSample*	_loadRealTimeFamData(char *Sp_buf, wsUint N_oriIdx);
	void		_registMissingFounder();

	wsUint*		_setSampOrder();
	wsUint*		_setVrtOrder();

	int			getDiploid(xVariant &X_vrt, char *Sp_buf, char **Sp_next);
public:
	/* Imported allele coding to ACGT coding if possible */
	void		code2acgt(xVariant &X_var, char B_indel) {
		if (!Na_ACGT) return;
		/* --varsubset */
		if (Xe_varSubset.size() && Xe_varSubset.find(X_var.name) == Xe_varSubset.end())
			return;

		if (B_indel) {
			for (char *x=X_var.indel1 ; *x ; x++) {
				wsUint X = *x;
				if (!Na_ACGT[X]) halt("Character [%c] in 1st indel on variant [%s] is not in --acgt",
					*x, X_var.name);
				*x = Na_ACGT[X];
			}
			for (char *x=X_var.indel2 ; *x ; x++) {
				wsUint X = *x;
				if (!Na_ACGT[X]) halt("Character [%c] in 2nd indel on variant [%s] is not in --acgt",
					*x, X_var.name);
				*x = Na_ACGT[X];
			}
		} else {
			if (!Na_ACGT[(wsUint)X_var.al1]) halt("Character [%c] in variant [%s] is not in --acgt",
				X_var.al1, X_var.name);
			if (!Na_ACGT[(wsUint)X_var.al2]) halt("Character [%c] in variant [%s] is not in --acgt",
				X_var.al2, X_var.name);
			X_var.al1 = Na_ACGT[(wsUint)X_var.al1];
			X_var.al2 = Na_ACGT[(wsUint)X_var.al2];
		}
	}

	cIO();
	virtual ~cIO();
	static cIO*	summon(xFileType X_type, char *Sp_fn, char B_dry=0);
	void		_finalize();
	void		clear();
	void		setPPP(cPPPAnalysisV2 *Cp_inpAnaPPP) {
		Cp_anaPPP = Cp_inpAnaPPP;
	}
	char*		getACGT();

	mStrReal&	getNewMAF() { return Xm_newMAF; }
	virtual
	wsStrCst	getFormat()=0;
	wsUintCst	getSampleSizeBeforeFiltering() { return N_sampOrig; }

	wsFmat		getData(char B_retOnly = 0);
	wsUint**	getGenoProb();
	void		setGenotype(char **Na_newData) {
		if (Na_geno) {
			for (wsUint i=0 ; i<N_sample ; i++)
				sseFree(Na_geno[i]);
			DEALLOC(Na_geno);
		}
		Na_geno = Na_newData;
	}
	void		setPhenotype(wsVec Ra_newPheno) {
		sseUnmat(Ra_pheno, N_pheno);
		N_pheno = 1;
		wsAlloc(Ra_pheno, wsVec, 1);
		Ra_pheno[0] = Ra_newPheno;
	}
	void		setReservedCol(xReservedCol X_col) {
		switch (X_col) {
		case RC_PROBAND:
			B_haveProband = 1; break;
		default: break;
		}
	}

	wsVec		getWeight(cPPPAnalysisV2 *Cp_anaPPP);
	USHORT_t**	getCharData() { return Na_charData; }
	const char*	getIsFounder() { return Ba_isFounder; }
	vCovar&		getCovInfo() { return Xa_covInfo; }
	void		chgCovInfo(wsUint N_nCov, wsMat Ra_nCov) {
		/* Assume that Xa_covInfo have already change
		 * and previously covars have linked to new covar
		 * and old pointer to covar already deallocated */
		N_cov	= N_nCov;
		Ra_cov	= Ra_nCov;
		if (Xa_covInfo.size() != N_cov)
			halt("SYSERR : New covariates size[%d] should be equal to "
				"requested covariates size[%d]", Xa_covInfo.size(), N_cov);
	}
	wsUintCst	setPhenoBuf(wsUintCst N_phe) {
		sseUnmat(Ra_pheno, N_pheno);
		DEALLOC(Ba_typePheno);
		Ra_pheno	= sseMatrix(N_phe, N_sample);
		wsCalloc(Ba_typePheno, char, N_phe);
		N_pheno		= N_phe;
		return N_pheno;
	}
	char*		getTypePheno() { return Ba_typePheno; }
	wsUint&		getFiltSampSize() { return N_szFilteredSamp; }
	wsUintCst	setCovBuf(wsUintCst N_sz, wsMat Ra_newCov=NULL) {
		sseUnmat(Ra_cov, N_cov);
		N_cov	= N_sz;
		if (Ra_newCov == NULL)
			Ra_cov	= sseMatrix(N_cov, N_sample);
		else
			Ra_cov	= Ra_newCov;
		pverbose("Covariate buffer size [%d * %d]\n", N_cov, N_sample);
		return N_cov;
	}
	void		setGxeBuf(vInt& Xa_idx) {
		Xa_gxeIdx.resize(Xa_idx.size());
		copy(Xa_idx.begin(), Xa_idx.end(), Xa_gxeIdx.begin());

		/* Make actual buffer */
		N_covGxE = (wsUint)Xa_gxeIdx.size();
		wsAlloc(Ra_covGxE, wsReal*, N_covGxE);
		for (wsUint i=0 ; i<N_covGxE ; i++)
			Ra_covGxE[i] = Ra_cov[Xa_gxeIdx[i]];
	}
	const char*	getChrAllowance() { return Ba_chrAllowed; }
	wsMat		getFiltNullDesignMat(const char *Ba_filt, char N_valFilt,
		wsUint N_szEst, wsUint N_moreRow=0, wsMat *Ra_dm=NULL);
	wsUintCst	getTwinSize() { return N_twin; }
	wsUintCst	getMissGenoCount() { return N_missGeno; }
	mFam&		getFamilyData() { return Xm_fid2data; }
	mSamp&		getSampleData() { return Xm_iid2data; }
	xMaf*		getMAF();
	xMaf*		getMAFbuffer();
	vPheno&		getPhenoInfo() { return V_phenos; }
	vStr&		getSampGrpIdx2grp() { return Sv_sampGrpIdx2grp; }

	/* Get the samples in the final dataset */
	vSampPtr&	getSample() { return Xa_sampleV2; }
	/* Get the variants in the final dataset */
	vVariant&	getVariant() { return Xv_variant; }
	/* Get the one phenotype */
	wsRealCst*	getPheno(int N_idxPheno=0) { return Ra_pheno[N_idxPheno]; }
	/* Get all phenotypes */
	wsReal**	getPhenos() { return Ra_pheno; }
	/* Get the covariates */
	wsMat		getCovariates();
	/* Get the covariates assigned to GxE */
	wsMat		getGxEcovariates();
	/* Get the genotypes */
	char**		getGenotype() { return Na_geno; }
	/* Get the dosages */
	wsFloat**	getDosage() { return Ra_dosage; }
	/* Get the raw genotypes */
	unsigned char**
				getRawGeno() { return Na_raw; }

	/* Get the number of samples in the final dataset */
	wsUintCst	sizeSample() { return N_sample; }
	/* Get the number of founders in the final dataset */
	wsUintCst	sizeFounder() { return N_founder; }
	/* Get the number of covariates in the final dataset */
	wsUintCst	sizeCovar();
	/* Get the number of GxE covariates */
	wsUintCst	sizeGxEcovar();
	/* Get the number of phenotypes in the final dataset */
	wsUintCst	sizePheno() { return N_pheno; }
	/* Get the number of variants in the final dataset */
	wsUintCst	sizeVariant() { return N_variant; }
	/* Set the number of twins in the final dataset */
	void		setTwinSize(wsUintCst N_sz) { N_twin = N_sz; }
	/* Set the number of samples in the alternative phenotype */
	void		setAltSampSize(wsUintCst N_sz) { N_altSample = N_sz; }

	vBool&		getIsSampleFiltered() { return Bv_filtSamp; }
//	REAL_c*		getCorrMask() { return Ra_corrMask; }

	char		isDataComplete() { return B_isDataComplete; }
	int			isSampleFiltered(string& S_FID, string &S_IID, char B_silent=0);
	int			isTwinAssigned() { return N_twin?1:0; }
	char		isProbandAssigned() { return B_haveProband?1:0; }
	char		isContinuous(int N_idx=0) { return Ba_typePheno[N_idx]; }
	char		isBLUPapplied() { return B_blupApplied; }
	const char*	isChrExists() { return Ba_chrExists; }

	void		exportPED();
	void		exportLGEN();
	void		exportGEN();
	void		exportBGEN(char B_compress=0);
	void		exportBeagle();
	void		exportMerlin();
	void		getChrwiseVariantCount(wsUint *Na_ret, char *Ba_vrtIncl);
	void		exportBED(wsUint B_vrtMajor, char *Ba_vrtIncl=NULL,
		char *Ba_sampIncl=NULL);
	void		exportBEDgene(wsUint B_vrtMajor, cSetManagerAnalysis *Cp_anaSM,
		char *Ba_sampIncl=NULL);

	wsVec		getSampleWeight() { return V_sampWgt.get(); }
	wsVec		setSampleWeightBuf() {
		wsReal R_val = WISARD_NAN;
		V_sampWgt.init(N_sample, NULL, &R_val);
		return V_sampWgt.get();
	}

	/* mode = 0 additive
	 * mode = 1 dominant
	 * mode = 2 recessive
	 */
	void		exportRAW(int mode=0);
	void		exportVCF(char B_bin=0);
	void		exportCovariates();
	void		exportPhenotypes();
	void		exportFamily();
	void		exportFamDiagram();
	void		exportListVariant();
	void		exportListSample();
	void		exportListFounder();

	wsReal**	getCnvData() { return Ra_cnv; }
	wsUintCst		getCnvSize() { return N_szCnv; }

	wsReal**	getFamInv(const char *Ba_filt, wsReal **Ra_corr,
		char B_filtExcl=1);
	char*		getPheCovMissing(wsUint *Np_nonMiss=NULL) {
		char*		Ba_ret		= NULL;
		wsMat		Ra_origCov	= getCovariates();
		wsMat		Ra_origY	= getPhenos();
		wsUint		N_origSamp	= sizeSample();
		wsUintCst	N_phe		= sizePheno();
		wsCalloc(Ba_ret, char, N_origSamp);

		/* Remove individuals having missing phenotypes or covariates */
		LOOP (j, N_phe) LOOP (i, N_origSamp)
			if (isMissingReal(Ra_origY[j][i])) Ba_ret[i] = 1;
		if (Ra_origCov) LOOP (i, N_cov) LOOP(j, N_origSamp)
			if (isMissingReal(Ra_origCov[i][j])) Ba_ret[j] = 1;

		/* If Np_nonMiss is required, set Np_nonMiss to # of nonmissing phenotype */
		if (Np_nonMiss) {
			*Np_nonMiss = 0;
			LOOP (i, N_origSamp)
				if (Ba_ret[i] == 0) (*Np_nonMiss)++;
		}

		return Ba_ret;
	}
	char*		getPheMissing(wsUint *Np_nonMiss=NULL) {
		wsUint		i, j;
		char	*	Ba_ret		= NULL;
		wsReal**	Ra_origY	= getPhenos();
		wsUint		N_origSamp	= sizeSample();
		wsUintCst	N_phe		= sizePheno();
		wsCalloc(Ba_ret, char, N_origSamp);

		/* Remove individuals having missing phenotypes or covariates */
		for (j=0 ; j<N_phe ; j++)
			for (i=0 ; i<N_origSamp ; i++)
				if (isMissingReal(Ra_origY[j][i])) Ba_ret[i] = 1;

		/* If Np_nonMiss is required, set Np_nonMiss to # of nonmissing phenotype */
		if (Np_nonMiss) {
			*Np_nonMiss = 0;
			for (i=0 ; i<N_origSamp ; i++)
				if (Ba_ret[i] == 0) (*Np_nonMiss)++;
		}

		return Ba_ret;
	}
	const char*	getCovMissing(wsUint *Np_nonMiss=NULL) {
		wsUint	i, j;
		char	*Ba_ret			= NULL;
		wsReal	**Ra_origCov	= getCovariates();
//		wsReal	**Ra_origY		= getPhenos();
		wsUint	N_origSamp		= sizeSample();
//		wsUintCst	N_pheno			= getPhenoCount();
		wsCalloc(Ba_ret, char, N_origSamp);

		/* Remove individuals having missing phenotypes or covariates */
		if (Ra_origCov)
			for (i=0 ; i<N_cov ; i++) for (j=0 ; j<N_origSamp ; j++)
				if (isMissingReal(Ra_origCov[i][j])) Ba_ret[j] = 1;

		/* If Np_nonMiss is required, set Np_nonMiss to # of nonmissing phenotype */
		if (Np_nonMiss) {
			*Np_nonMiss = 0;
			for (i=0 ; i<N_origSamp ; i++)
				if (Ba_ret[i] == 0) (*Np_nonMiss)++;
		}

		return Ba_ret;
	}
	/* Get indices of phenotype-available sample and their number
	 * (optionally, their reduced array (sse'd) */
	const char*	getAvailPheno(wsUint *Np_availSamp, wsReal **Rp_availY=NULL);

	/* Get the names of covariates as string vector
	 * if factor is included in covariates, each level will be exported */
	void	getCovNames(vStr &Sa_ret) {
		FOREACH (vCovar_it, Xa_covInfo, it) {
			if (it->X_type == WISARD_VAR_UNDET)
				halt_fmt(WISARD_SYST_INVL_COVTYPE, it->Sp_varName);

			Sa_ret.push_back(it->Sp_varName);
		}
	}
	char		_filterVariantPosition(xVariant &X_vrt);
	void		exportSampleMap();
	wsReal*		subsetPheCov(char *Ba_filt, char B_incVal, wsUint N_sz,
		wsMat *Rp_cov);
	wsReal*		subsetPheCovItct(char *Ba_filt, char B_incVal, wsUint N_sz,
		wsMat *Rp_cov, wsUint *Np_cov=NULL);
	void		impute(cAnalysis *Cp_anaPhi);
	cImpute*	getImpute() { return Cp_impute; }
	wsMat		getImputeGeno();
	void		setDataSize(wsUint N_inSamp, wsUint N_inMarker) {
		N_sample = N_inSamp;
		N_variant = N_inMarker;
	}
	void		getPhenoCode(wsUint N_idx, char *Sp_buf, wsStrCst Sp_case,
		wsStrCst Sp_ctrl, int N_idxPheno=0)
	{
		if (N_pheno == 0 || !Ra_pheno || !Ra_pheno[N_idxPheno]) {
			strcpy(Sp_buf, "-9");
			return;
		}
		wsReal R_pheno = Ra_pheno[N_idxPheno][N_idx];
		if (isMissingReal(R_pheno))
			sprintf(Sp_buf, "%s", IS_ASSIGNED(outmispheno)?OPT_STRING(outmispheno):"-9");
		else if (isContinuous())
			sprintf(Sp_buf, "%g", R_pheno);
		else
			sprintf(Sp_buf, "%s", R_pheno==WISARD_AFFECTED ? Sp_case : Sp_ctrl);

	}
	vInt&		getGxEcovariatesIndices() { return Xa_gxeIdx; }
};

void exportSampleMatrix(wsStrCst S_ext, cIO *Cp_IO, char *B_incl,
	char B_filtVal, wsUint N_mat, ...);
void exportSampleMatrix(wsStrCst S_ext, cIO *Cp_IO,
	const char *Ba_filter, char B_filtVal,
	wsUint N_mat, char **Sp_name, wsReal **Ra_mats);
void otherExt(wsStrCst S_in, wsStrCst S_ext, char *Sp_dupFn, wsUint N_lvRem=1);
void loadEstimates(wsStrCst Sp_strFn, wsMat Ra_estimates, wsUint N_pheno,
	wsUint N_cov);

inline wsReal _imputeGeno(xVariant &X_vrt, xSample *Xp_samp, wsReal R_maf)
{
	wsReal R_mult = W2;
	if (isXChromosome(X_vrt)) {
		switch (Xp_samp->N_sex) {
		case 2: break;
		case 1: R_mult = W1; break;
		default: R_mult = W0; break;
		}
	} else if (isYChromosome(X_vrt)) {
		switch (Xp_samp->N_sex) {
		case 1: R_mult = W1; break;
		default: R_mult = W0; break;
		}
	}
	return R_maf*R_mult;
}

inline wsUint setGeno(cIO *Cp_IO, wsUint N_mIdx, wsReal *Ra_t, const char *Ba_mask,
	char B_isExcl=0)
{
	wsUint		N_samp	= Cp_IO->sizeSample();
	xVariant&	X_mkr	= Cp_IO->getVariant()[N_mIdx];
	xMaf*		Xp_maf	= Cp_IO->getMAF();
	vSampPtr&	Xa_samp	= Cp_IO->getSample();
	wsMat		Ra_data	= Cp_IO->getImputeGeno();
	wsFloat**	Ra_dos	= !Ra_data && IS_ASSIGNED(dosage) ? Cp_IO->getData() : NULL;
	wsUint**	Na_gprob= !Ra_data && IS_ASSIGNED(genoprob) ? Cp_IO->getGenoProb() : NULL;
	char**		Na_data	= Cp_IO->getGenotype();

	/* Have imputed data */
	if (Ra_data) {
		if (Ba_mask == NULL)
			/* Just fill, because it is already FULLY imputed */
			memcpy(Ra_t, Ra_data[N_mIdx], sizeof(wsReal)*N_samp);
		else if (!B_isExcl) for (wsUint i=0,I=0 ; i<N_samp ; i++) {
			if (Ba_mask[i]) Ra_t[I++] = Ra_data[N_mIdx][i];
		} else for (wsUint i=0,I=0 ; i<N_samp ; i++)
			if (!Ba_mask[i]) Ra_t[I++] = Ra_data[N_mIdx][i];

		/* Now have full genotype */
		return 0;
	} else if (Ra_dos) {
		wsUint N_miss = 0;
		if (Ba_mask == NULL) for (wsUint i=0 ; i<N_samp ; i++) {
                        if (isMissingReal(Ra_dos[i][N_mIdx])) {
                                N_miss++;
                                Ra_t[i]         = W0;
                        } else Ra_t[i]  = (wsReal)Ra_dos[i][N_mIdx];
                } else if (!B_isExcl) for (wsUint i=0,I=0 ; i<N_samp ; i++) {
                        if (Ba_mask[i] == 0) continue;

                        if (isMissingReal(Ra_dos[i][N_mIdx])) {
                                N_miss++;
                                Ra_t[I++]       = W0;
                        } else Ra_t[I++]        = (wsReal)Ra_dos[i][N_mIdx];
                } else for (wsUint i=0,I=0 ; i<N_samp ; i++) {
                        if (Ba_mask[i]) continue;
                        if (isMissingReal(Ra_dos[i][N_mIdx])) {
                                N_miss++;
                                Ra_t[I++]       = W0; //_imputeGeno(X_mkr, Xa_samp[i], Xp_maf[i].R_maf);
                        } else Ra_t[I++]        = (wsReal)Ra_dos[i][N_mIdx];
                }
                /* Now have full genotype */
                return N_miss;
	}

	wsUint N_miss = 0;
	/* Just impute */
//	if (!X_mkr.complete) {
		if (Ba_mask == NULL) for (wsUint i=0 ; i<N_samp ; i++) {
			if (isMissing(Na_data[i][N_mIdx])) {
				N_miss++;
				Ra_t[i]		= _imputeGeno(X_mkr, Xa_samp[i], Xp_maf[i].R_maf);
			} else Ra_t[i]	= (wsReal)Na_data[i][N_mIdx];
		} else if (!B_isExcl) for (wsUint i=0,I=0 ; i<N_samp ; i++) {
			if (Ba_mask[i] == 0) continue;
			
			if (isMissing(Na_data[i][N_mIdx])) {
				N_miss++;
				Ra_t[I++]		= _imputeGeno(X_mkr, Xa_samp[i], Xp_maf[i].R_maf);
			} else Ra_t[I++]	= (wsReal)Na_data[i][N_mIdx];
		} else for (wsUint i=0,I=0 ; i<N_samp ; i++) {
			if (Ba_mask[i]) continue;
			if (isMissing(Na_data[i][N_mIdx])) {
				N_miss++;
				Ra_t[I++]		= _imputeGeno(X_mkr, Xa_samp[i], Xp_maf[i].R_maf);
			} else Ra_t[I++]	= (wsReal)Na_data[i][N_mIdx];
		}
// 	} else {
// 		if (Ba_isInc == NULL) for (wsUint i=0 ; i<N_samp ; i++)
// 			Ra_t[i]		= Na_data[i][N_mIdx];
// 		else if (!B_isExcl) for (wsUint i=0,I=0 ; i<N_samp ; i++) {
// 			if (Ba_isInc[i]) Ra_t[I++]	= Na_data[i][N_mIdx];
// 		} else for (wsUint i=0,I=0 ; i<N_samp ; i++)
// 			if (!Ba_isInc[i]) Ra_t[I++]= Na_data[i][N_mIdx];
// 	}
	return N_miss;
}

inline wsUint setGxGeno(cIO *Cp_IO, wsUint N_m1, wsUint N_m2,
	wsReal *Ra_t1, wsReal *Ra_t2, wsReal *Ra_tx,
	const char *Ba_mask, char B_isExcl=0)
{
	wsUint		N_samp	= Cp_IO->sizeSample();
	vVariant&	Xv_var	= Cp_IO->getVariant();
	xMaf*		Xp_maf	= Cp_IO->getMAF();
	xVariant&	X_m1	= Xv_var[N_m1];
	xVariant&	X_m2	= Xv_var[N_m2];
	vSampPtr&	Xa_samp	= Cp_IO->getSample();
	wsMat		Ra_data	= Cp_IO->getImputeGeno();
	char**		Na_data	= Cp_IO->getGenotype();

	/* Have imputed data */
	if (Ra_data) {
		if (Ba_mask == NULL) {
			/* Just fill, because it is already FULLY imputed */
			memcpy(Ra_t1, Ra_data[N_m1], sizeof(wsReal)*N_samp);
			memcpy(Ra_t2, Ra_data[N_m2], sizeof(wsReal)*N_samp);
			sseVpV(Ra_t1, Ra_t2, Ra_tx, N_samp);
		} else if (!B_isExcl) for (wsUint i=0,I=0 ; i<N_samp ; i++) {
			if (!Ba_mask[i]) continue;
			Ra_t1[I] = Ra_data[N_m1][i];
			Ra_t2[I] = Ra_data[N_m2][i];
			Ra_tx[I] = Ra_t1[I] * Ra_t2[I];
			I++;
		} else for (wsUint i=0,I=0 ; i<N_samp ; i++) {
			if (Ba_mask[i]) continue;
			Ra_t1[I] = Ra_data[N_m1][i];
			Ra_t2[I] = Ra_data[N_m2][i];
			Ra_tx[I] = Ra_t1[I] * Ra_t2[I];
			I++;
		}

		/* Now have full genotype */
		return 0;
	}

	wsUint N_miss = 0;
	/* Just impute */
	//	if (!X_mkr.complete) {
	if (Ba_mask == NULL) for (wsUint i=0 ; i<N_samp ; i++) {
		if (isMissing(Na_data[i][N_m1])) {
			N_miss++;
			Ra_t1[i] = _imputeGeno(X_m1, Xa_samp[i], Xp_maf[i].R_maf);
		} else Ra_t1[i]	= (wsReal)Na_data[i][N_m1];
		
		if (isMissing(Na_data[i][N_m2])) {
			N_miss++;
			Ra_t2[i] = _imputeGeno(X_m2, Xa_samp[i], Xp_maf[i].R_maf);
		} else Ra_t2[i]	= (wsReal)Na_data[i][N_m2];

		Ra_tx[i] = Ra_t1[i] * Ra_t2[i];
	} else if (!B_isExcl) for (wsUint i=0,I=0 ; i<N_samp ; i++) {
		if (Ba_mask[i] == 0) continue;

		if (isMissing(Na_data[i][N_m1])) {
			N_miss++;
			Ra_t1[I]	= _imputeGeno(X_m1, Xa_samp[i], Xp_maf[i].R_maf);
		} else Ra_t1[I]	= (wsReal)Na_data[i][N_m1];

		if (isMissing(Na_data[i][N_m2])) {
			N_miss++;
			Ra_t2[I]	= _imputeGeno(X_m2, Xa_samp[i], Xp_maf[i].R_maf);
		} else Ra_t2[I]	= (wsReal)Na_data[i][N_m2];

		Ra_tx[I] = Ra_t1[I] * Ra_t2[I];
		I++;
	} else for (wsUint i=0,I=0 ; i<N_samp ; i++) {
		if (Ba_mask[i]) continue;

		if (isMissing(Na_data[i][N_m1])) {
			N_miss++;
			Ra_t1[I]	= _imputeGeno(X_m1, Xa_samp[i], Xp_maf[i].R_maf);
		} else Ra_t1[I]	= (wsReal)Na_data[i][N_m1];

		if (isMissing(Na_data[i][N_m1])) {
			N_miss++;
			Ra_t2[I]	= _imputeGeno(X_m1, Xa_samp[i], Xp_maf[i].R_maf);
		} else Ra_t2[I]	= (wsReal)Na_data[i][N_m2];

		Ra_tx[I] = Ra_t1[I] * Ra_t2[I];
		I++;
	}
	// 	} else {
	// 		if (Ba_isInc == NULL) for (wsUint i=0 ; i<N_samp ; i++)
	// 			Ra_t[i]		= Na_data[i][N_mIdx];
	// 		else if (!B_isExcl) for (wsUint i=0,I=0 ; i<N_samp ; i++) {
	// 			if (Ba_isInc[i]) Ra_t[I++]	= Na_data[i][N_mIdx];
	// 		} else for (wsUint i=0,I=0 ; i<N_samp ; i++)
	// 			if (!Ba_isInc[i]) Ra_t[I++]= Na_data[i][N_mIdx];
	// 	}
	return N_miss;
}

#define HEADER_VRT_1	"CHR"
#define HEADER_VRT_2	"VARIANT"
#define HEADER_VRT_3	"POS"
#define HEADER_VRT_4	"ALT"
#define HEADER_VRT_5	"ANNOT"

inline void headerVariant(cExporter *Cp)
{
	/* Note : If this function is changed, fetchMarkerList() must be changed */
	if (IS_ASSIGNED(annogene))
		Cp->put(HEADER_VRT_1 "	" HEADER_VRT_2 "	" HEADER_VRT_3 "	" HEADER_VRT_4 "	" HEADER_VRT_5);
	else
		Cp->put(HEADER_VRT_1 "	" HEADER_VRT_2 "	" HEADER_VRT_3 "	" HEADER_VRT_4);
}

inline wsUint headerVariant(vStr& Sv_hdrs, char* S_fmt=NULL)
{
	wsUint N_elem = 4; /* In default, # of elements is 4 */
	Sv_hdrs.push_back(HEADER_VRT_1);
	Sv_hdrs.push_back(HEADER_VRT_2);
	Sv_hdrs.push_back(HEADER_VRT_3);
	Sv_hdrs.push_back(HEADER_VRT_4);
	if (S_fmt) strcpy(S_fmt, "ssis");

	/* Note : If this function is changed, fetchMarkerList() must be changed */
	if (IS_ASSIGNED(annogene)) {
		Sv_hdrs.push_back(HEADER_VRT_5);
		S_fmt[N_elem++] = 's';
		S_fmt[N_elem] = '\0';
	}

	return N_elem;
}

inline wsUint entryVariant(cExporter *Cp, xVariant &Xs)
{
	if (IS_ASSIGNED(annogene)) {
		if (OPT_ENABLED(indel))
			Cp->fmt("%s	%s	%d	%s	%s", getChrName2(Xs.chr), Xs.name,
				Xs.pos, Xs.indel2?Xs.indel2:"<NA>", Xs.anno);
		else if (Xs.al2)
			Cp->fmt("%s	%s	%d	%c	%s", getChrName2(Xs.chr), Xs.name,
				Xs.pos, Xs.al2, Xs.anno);
		else
			Cp->fmt("%s	%s	%d	<NA>	%s", getChrName2(Xs.chr), Xs.name,
				Xs.pos, Xs.anno);
		return 5;
	} else {
		if (OPT_ENABLED(indel))
			Cp->fmt("%s	%s	%d	%s", getChrName2(Xs.chr), Xs.name,
				Xs.pos, Xs.indel2?Xs.indel2:"<NA>");
		else if (Xs.al2)
			Cp->fmt("%s	%s	%d	%c", getChrName2(Xs.chr), Xs.name,
				Xs.pos, Xs.al2);
		else
			Cp->fmt("%s	%s	%d	<NA>", getChrName2(Xs.chr), Xs.name,
				Xs.pos);
		return 4;
	}
}

inline void getOutCaCt(char **Sp_case, char **Sp_ctrl)
{
	if (IS_ASSIGNED(outcact)) {
		wsUint N_cact = 0;
		char **Sa_cact = loadStringValues(OPT_STRING(outcact), &N_cact);
		if (N_cact != 2)
			halt("Parameter of --outcact must be two strings divided by comma without whitespace!");
		*Sp_case = Sa_cact[0];
		*Sp_ctrl = Sa_cact[1];
		DEALLOC(Sa_cact);
	} else {
		/* Set case/control output character */
		*Sp_case = OPT_ENABLED(out1case) ? Z("1") : Z("2");
		*Sp_ctrl = OPT_ENABLED(out1case) ? Z("0") : Z("1");
	}
}

inline char* loadFIDIID(char* Sp_buf, string& S_curFID, string& S_curIID)
{
	char	*Sp_tmp1 = NULL;
	char	*Sp_tmp2 = NULL;

	if (OPT_ENABLED(nofid) || IS_ASSIGNED(sepid)) Sp_tmp1 = Sp_buf;
	else {
		getString(&Sp_buf, &Sp_tmp1);
		S_curFID = (string)Sp_buf;
	}
	getString(&Sp_tmp1, &Sp_tmp2);

	/* In case of --sepid, separate IID by --sepid and [0]=FID, [1]=IID */
	if (IS_ASSIGNED(sepid)) {
		char* Sp_sep = strstr(Sp_tmp1, OPT_STRING(sepid));
		if (!Sp_sep)
			halt("IID [%s] does not have valid separator [%s]", S_curIID.c_str(), OPT_STRING(sepid));
		*Sp_sep = '\0';
		S_curFID = Sp_tmp1;
		S_curIID = Sp_tmp1+strlen(OPT_STRING(sepid));
	}
	else S_curIID = (string)Sp_tmp1;

	/* In case of --nofid/--ignorefid, assume FID is equal to IID */
	if (OPT_ENABLED(nofid) || OPT_ENABLED(ignorefid)) S_curFID = S_curIID;

	return Sp_tmp2;
}

class cSymMatrix;
wsReal* getCondVar(cIO* Cp_IO, vSampPtr& Xa_samp, cSymMatrix& M_cor, wsMat *Rp_vt=NULL);

} // End namespace ONETOOL

#endif
