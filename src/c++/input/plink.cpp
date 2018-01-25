#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <list>
#include <vector>
#include "input/plink.h"
#include "global/option.h"
#include "utils/util.h"
#include "global/worker.h"
#include "global/io.h"
#include "input/stream.h"

namespace ONETOOL {

typedef struct _xRangedThreadBedLoad : public xRangedThread {
	size_t	N_missGeno;
	wsUint*	Na_sampMissGeno;
	wsUint	N_szFilteredSamp;
} xRangedThreadBedLoad;

typedef struct _xThreadBedLoad  {
	xThread		t;
	wsUint		N_fSamp;
	wsUint		N_fileVrt;
	cIO*		Cp_IO;
	char		B_vrtMajor;
	unsigned char *
				Na_data;
	vBool&		Bv_filtVrtOrig;
	vBool&		Bv_filtSampOrig;
	vBool&		Bv_filtSamp;
} xThreadBedLoad;

inline char setINDEL(xVariant &X_vrt, char *Sp_gen1)
{
	/* If it is missing */
	if (!stricmp(Sp_gen1, OPT_STRING(misgeno)))
		return WISARD_NA;

	if (X_vrt.indel1 == NULL)
		X_vrt.indel1 = strdup(Sp_gen1);
	else {
		/* Have indel1 */

		/* Same with indel1? */
		if (stricmp(X_vrt.indel1, Sp_gen1)) {
			/* indel1 != this */

			if (X_vrt.indel2 == NULL) {
				/* indel2 == NULL */
				X_vrt.indel2 = strdup(Sp_gen1);
			} else if (stricmp(X_vrt.indel2, Sp_gen1))
				/* Trialleleic */
				return 0;

			/* Match with second allele */
			return '2';
		}
	}

	return '1';
}

inline void _indelSub(char *Sp_gen1, char *Sp_gen2, xVariant &X_vrt,
	wsUint N_fileIdx, USHORT_t &N_geno, char *Ba_chrFound)
{
	char B_idx1, B_idx2;
	if (!stricmp(OPT_STRING(misgeno), Sp_gen1)) {
		if (stricmp(OPT_STRING(misgeno), Sp_gen2))
			halt_fmt(WISARD_CANT_ONEALLELE_MISSING, X_vrt.name, N_fileIdx);
		// 				halt("%dth sample, %dth variant, PED missing should be both!",
		// 					fiIdx, fsIdx);
		B_idx1 = '0';
		B_idx2 = '0';
	} else {
		B_idx1 = setINDEL(X_vrt, Sp_gen1);
		B_idx2 = setINDEL(X_vrt, Sp_gen2);

		/* Triallelic */
		if (!B_idx1 || !B_idx2) {
			X_vrt.filter = 1;
			//LOG("Variant [%s] will be excluded\n", X_vrt.name);
		}

		/* At this point, at least one indel1 */
		if (stricmp(X_vrt.indel1, Sp_gen2)) {
			if (X_vrt.indel2 == NULL)
				X_vrt.indel2 = strdup(Sp_gen2);
		}

		/* --snvonly but indel */
		if (X_vrt.indel1[1] || (X_vrt.indel2 && X_vrt.indel2[1]) ||
			X_vrt.indel1[0] == '-' || (X_vrt.indel2 && X_vrt.indel2[0] == '-'))
			X_vrt.filter = OPT_ENABLED(snvonly);
		/* --indelonly & it is SNV then filter */
		if (!X_vrt.indel1[1] && !X_vrt.indel2[1] &&
			X_vrt.indel1[0] != '-' && X_vrt.indel2[0] != '-')
			X_vrt.filter = OPT_ENABLED(indelonly);
	}

	Ba_chrFound[(int)Sp_gen1[0]] = 1;
	Ba_chrFound[(int)Sp_gen2[0]] = 1;

	N_geno = (B_idx1<<8)|B_idx2;
}

cPlinkIO::cPlinkIO(char *S_fn, xFileType X_type, char B_inDry) : cIO()
{
	B_dry = B_inDry;
	init(S_fn, X_type);
}

void cPlinkIO::init(char *S_fn, xFileType X_type)
{
	X_fileType = X_type;

	switch (X_fileType) {
	case FT_TPED:	_loadTPED(S_fn); break;
	case FT_LGEN:	_loadLGEN(S_fn); break;
	case FT_PED:	_loadPED(S_fn); break;
	case FT_BED:	_loadBED(S_fn); break;
	default:
		halt("SYSERR: Unexpected filetype in PLINK I/O");
		break;
	}
}

void cPlinkIO::_loadBED_BIM(char *S_bedFn)
{
	char	S_fnBIM[MAX_PATH]		= { 0, };

	/*
	 * Set BIM load information
	 */
	if (IS_ASSIGNED(bim)) {
		strcpy(S_fnBIM, OPT_STRING(bim));
		LOG("Using alternative BIM file for BED loading [%s]\n",
			OPT_STRING(bim));
	} else if (S_bedFn[0]=='.'&&!S_bedFn[1])
		halt_fmt(WISARD_CANT_DO_WO_SOMEOPT, "BED loading from STDIN",
			"--bim and --fam");
	else
		otherExt(S_bedFn, "bim", S_fnBIM);

	/*
	 * Initiate BIM file stream
	 */
	cStrFile C_bim(S_fnBIM, "BIM file");

	char *Ba_isSelResize	= NULL;
	if (IS_ASSIGNED(varresize)) {
		wsReal R_valSS = OPT_REAL(varresize);
		wsUint N_sel = R_valSS >= W1 ?
			(wsUint)R_valSS : (wsUint)round(R_valSS*N_vrtOrig);

		/* Retrieve the number of variants in BIM file first */
		char *Sp_buf = NULL;
		for (N_vrtOrig=0 ; (Sp_buf=C_bim.ugets(1)) ; N_vrtOrig++) {
			DEALLOC(Sp_buf);
		}
		C_bim.rewind();

		/* Check it should be resized */
		if (N_vrtOrig <= N_sel) {
			LOG("# of variants in dataset[%d] is larger than the number of "
				"variants to be resized[%d] so the dataset will not be resized\n",
				N_vrtOrig, N_sel);
		} else {
			/* Alloc memory */
			wsCalloc(Ba_isSelResize, char, N_vrtOrig);

			/* Randomly select variants */
			for (wsUint N_chosen=0 ; N_chosen<N_sel ; ) {
				/* Select one */
				wsUint N_curr = rand()%N_vrtOrig;

				if (Ba_isSelResize[N_curr] == 0) {
					Ba_isSelResize[N_curr] = 1;
					N_chosen++;
				}
			}
		}		
	}

	/* chr	name	pos	bp	al1	al2 */
	char *a, *b;
	N_variant = 0;
	char *Sp_buf = NULL;
	for (N_vrtOrig=0 ; (Sp_buf=C_bim.ugets()) ; N_vrtOrig++) {
		xVariant X_ent = {0, };
		char B_filter = 0;

		// S_buf = chr
		getString(&Sp_buf, &a);
		X_ent.chr = getChr(Sp_buf);

		/* [FILTERING]
		 * Variant filtering with --autoonly/--sexonly */
		if (OPT_ENABLED(autoonly) || OPT_ENABLED(sexonly)) {
			/* This variant should be filtered */
			if (X_ent.chr == -1 ||
				(OPT_ENABLED(autoonly) && X_ent.chr > (int)NAUTO_SPECIES) ||
				(OPT_ENABLED(sexonly) && X_ent.chr <= (int)NAUTO_SPECIES))
				B_filter = 1;
		}
	
		// a = name
		getString(&a, &b);
		wsAlloc(X_ent.name, char, strlen(a)+1);
		strcpy(X_ent.name, a);

		/* [FILTERING]
		 * Variant filtering with --varresize */
		if (Ba_isSelResize && Ba_isSelResize[N_vrtOrig]==0)
			B_filter = 1; /* Filter if not selected */

		/* b = genetic distance */
		if (!OPT_ENABLED(nogdist)) {
			getString(&b, &a);
			X_ent.gdist = (wsReal)atof(b);
		} else {
			X_ent.gdist = W0;
			a = b;
		}

		/* a = pos */
		if (!OPT_ENABLED(nopos)) {
			getString(&a, &b);
			for (char *p=a ; *p ; p++)
				if (*p < '0' || *p > '9')
					halt_fmt(WISARD_INVL_VARIANTPOS, X_ent.name, a);
			int N_ret = atoi(a);
			if (N_ret < 0) B_filter = 1;
			X_ent.pos = (wsUint)N_ret;
		} else {
			X_ent.pos = 0;
			b = a;
		}

		/* [FILTERING]
		 * Check variant filter */
		B_filter |= _filterVariantPosition(X_ent);

		X_ent.indel1	= NULL;
		X_ent.indel2	= NULL;
		if (b && a) {
			getString(&b, &a);
			/* Do trimming */
			for (char *q=a+strlen(a)-1 ; a<=q && (*q == '\n' || *q == '\t' || *q == ' ' || *q == '\r') ; q--) {
				if (*q == '\0')
					halt_fmt(WISARD_INVL_BIMFORMAT, N_vrtOrig+1);
				*q = '\0';
			}

			if (OPT_ENABLED(indel)) {
				X_ent.indel1 = strdup(a);
				X_ent.indel2 = strdup(b);

				/* --acgt w/ --indel */
				code2acgt(X_ent, 1);

				/* --snvonly & it is indel then filter */
				if (X_ent.indel1[1] || X_ent.indel2[1] ||
					X_ent.indel1[0] == '-' || X_ent.indel2[0] == '-')
					B_filter |= OPT_ENABLED(snvonly);
				/* --indelonly & it is SNV then filter */
				if (!X_ent.indel1[1] && !X_ent.indel2[1] &&
					X_ent.indel1[0] != '-' && X_ent.indel2[0] != '-')
					B_filter |= OPT_ENABLED(indelonly);
			} else {
				if (b[1] || a[1])
					halt_fmt(WISARD_CANT_GETINDEL_WO_INDEL, N_vrtOrig+1);

				// b = al1, a = al2
				if (b[0] == '0') {
					X_ent.al1 = a[0] == '0' ? '\0' : a[0];
					X_ent.al2 = '\0';
				} else if (a[0] == '0') {
					X_ent.al1 = b[0] == '0' ? '\0' : b[0];
					X_ent.al2 = '\0';
				} else {
					X_ent.al1	= a[0];
					X_ent.al2	= b[0];
				}

				/* --acgt w/o --indel */
				code2acgt(X_ent, 0);
			}
		} else {
			if (Na_ACGT)
				halt_fmt(WISARD_INVL_FILE_INVALID, "BIM file", S_fnBIM,
					N_vrtOrig+1);
//				halt("Variant [%s] does not have allele info, --acgt makes halt", X_ent.name);

			/* Just assigns as 1/2 */
			if (OPT_ENABLED(indel)) {
				X_ent.indel1 = strdup("1");
				X_ent.indel2 = strdup("2");
			} else {
				X_ent.al1	= '1';
				X_ent.al2	= '2';
			}
		}
		X_ent.filter	= 0;

		if (B_filter == 1)
			Bv_filtVrt.push_back(1);
		else {
			N_variant++;
			Bv_filtVrt.push_back(0);
		}

		if (B_filter == 0)
			Xv_variant.push_back(X_ent);
		else {
			DEALLOC(X_ent.name);
			if (X_ent.indel1) free(X_ent.indel1);
			if (X_ent.indel2) free(X_ent.indel2);
		}
		DEALLOC(Sp_buf);
	}
	LOG("Variant loaded [BIM format] [%d/%d] chosen\n", N_variant, N_vrtOrig);

	DEALLOC(Ba_isSelResize);
}

void cPlinkIO::_loadBED_FAM(char *S_bedFn)
{
	char	S_famFileName[MAX_PATH];
	char*	S_buf = NULL;
	wsAlloc(S_buf, char, 4096);

	/* Reading FAM file next */
	if (IS_ASSIGNED(fam)) {
		/* Use alternative fam file if --fam given */
		LOG("Use alternative FAM file [%s] instead of expected fam file [%s.fam]\n",
			OPT_STRING(fam), S_bedFn);
		sprintf(S_famFileName, "%s", OPT_STRING(fam));
	} else if (S_bedFn[0]=='.'&&!S_bedFn[1])
		halt_fmt(WISARD_CANT_DO_WO_SOMEOPT, "BED loading from STDIN",
			"--bim and --fam");
	else
		otherExt(S_bedFn, "fam", S_famFileName);

	/* Initiate FAM file stream */
	cStrFile C_fam(S_famFileName, "FAM file");

	/* Reading... */
	N_sample = 0;
	for (N_sampOrig=0 ; C_fam.gets(S_buf, 4096) ; N_sampOrig++)
		_loadRealTimeFamData(S_buf, N_sampOrig);
	C_fam.rewind();

	/* Make complete family data */
	_registMissingFounder();

	/* Init memories */
	_initMemories();

	/* Read FAM file : main */
	N_founder = 0;
	char *Sp_tmp1 = NULL;
	char *Sp_tmp2 = NULL;
	/* should be --1 but forget check */
	wsUint N_0 = 0;
	char B_noPheno = 0;
	for (wsUint i=0,j=0 ; i<N_sampOrig ; i++) {
		if (C_fam.gets(S_buf, 1024) == 0)
			halt_fmt(WISARD_INVL_EOF_LINE, "FAM file", S_famFileName,
				N_sampOrig, i+1);
//			halt("Unexpected end of file found at line %d", i+1);
		string	S_currFID, S_currIID;

		/* Get FID & IID */
		Sp_tmp2 = loadFIDIID(S_buf, S_currIID, S_currIID);

		xSample	&X_sample = Xm_iid2data[S_currIID];
		if (!OPT_ENABLED(noparent)) {
			getString(&Sp_tmp2, &Sp_tmp1);	// Sp_tmp2 = PAT
			getString(&Sp_tmp1, &Sp_tmp2);	// Sp_tmp1 = MAT
		}
		if (OPT_ENABLED(nosex)) Sp_tmp1 = Sp_tmp2;
		else getString(&Sp_tmp2, &Sp_tmp1);	// Sp_tmp2 = SEX

		/* If current samples should be filtered, do not retrieve anything */
		if (Bv_filtSamp[i] == 1) {
			X_sample.N_idx = SAMP_NODATA;
			continue;
		}

		/* Is this sample founder? */
		if ((Ba_isFounder[j] = !X_sample.Xp_mat && !X_sample.Xp_pat) != 0)
			N_founder++;
		/* Allocate initial data index */
		X_sample.N_idx = j;
		Xa_sampleV2.push_back(&X_sample);

		/* a = phenotype, retrieve anyway, because alternative phenotype should
		 * processed AFTER this procedure */
		char *Sp_shouldNull = NULL;
		if (!OPT_ENABLED(nopheno)) {
			if (!Sp_tmp1) {
				if (!B_noPheno) {
					LOGwarn("No phenotype column found, the phenotype will be NA unless alternative phenotype is provided!\n");
					B_noPheno = 1;
				}
			} else {
				getString(&Sp_tmp1, &Sp_tmp2);
				Ra_pheno[0][j] = (wsReal)str2dbl(Sp_tmp1, &Sp_shouldNull);
			}
		} else {
			Ra_pheno[0][j] = WISARD_NA_REAL;
			Sp_tmp2 = Sp_tmp1;
		}

		/* Phenotype value is not control nor case */
		if (!B_noPheno) {
			if (Ra_pheno[0][j] != (wsReal)OPT_NUMBER(phenoCtrl) &&
				Ra_pheno[0][j] != (wsReal)OPT_NUMBER(phenoCase)) {

				if (!strcmp(Sp_tmp1, OPT_STRING(mispheno))) {
					/* It is missing phenotype */
					Ra_pheno[0][j] = WISARD_NA_REAL;
	//				LOG("Missing phenotype sample [%s:%s]\n",
	//					X_sample.S_FID.c_str(), X_sample.S_IID.c_str());
				} else {
					if (Sp_shouldNull && Sp_shouldNull[0])
						halt_fmt(WISARD_INVL_FILE_INVALID_DESC, "FAM file",
							S_famFileName, "phenotype", Sp_tmp1, i+1);
	// 					halt("Phenotype [%s] for sample [%s] looks invalid",
	// 					Sp_tmp1, X_sample.S_FID.c_str());

					/* Special for 0 */
					if (!strcmp(Sp_tmp1, "0")) N_0++;
					/* Otherwise, it is continuous phenotype */
					Ba_typePheno[0] = 1;
				}
			}
		}
		/* Filtering option */
// 		else if ((OPT_ENABLED(filcase) && Ra_pheno[0][j] == OPT_NUMBER(phenoCase)) ||
// 			(OPT_ENABLED(filcontrol) && Ra_pheno[0][j] == OPT_NUMBER(phenoCtrl)))
// 			Ba_isSampFilt[j] = 1;

		/* Status */
		if (((++j)%100) == 0) notice("%d samples retrieved\r", j);
	}
	/* Possible chance to miss --1 */
	if (N_0) {
		LOGwarn("This file could be coded 0/1 but not assumed that. Please check it!\n");
		LOG("        ( 0 found %d times )\n", N_0);
	}

	//fclose(H_famFile);
	LOG("Sample loaded [FAM format] [%d/%d] chosen\n", N_sample, N_sampOrig);

	/* Assign phenotype count to 1 
	 * If alternative phenotype is given, this number will be adjusted */
	if (!OPT_ENABLED(nopheno))
		N_pheno = 1;
	else
		N_pheno = 0;

	DEALLOC(S_buf);
}

char cPlinkIO::_loadBED_header(cStream &C_bed)
{
	// 1) Check for magic number
	// 2) else check for 0.99 variant/Ind coding
	// 3) else print warning that file is too old
	unsigned char ch;
	if (C_bed.read(&ch) != 1)
		halt_fmt(WISARD_INVL_BEDSIGNATURE);
	bool B_vrtMajor = false;
	bool v1_bfile = true;

	// If v1.00 file format
	// Magic numbers for .bed file: 00110110 11011000 = v1.00 bed file
	//	LOG("BED file sig #1 : %02X\n", ch);
	if (ch == 0x6C) {	
		if (C_bed.read(&ch) != 1)
			halt_fmt(WISARD_FAIL_READ_BEDSIG);
//			halt("Failed to read signature of BED file!");
		//		LOG("BED file sig #2 : %02X\n", ch);
		if (ch == 0x1B) {
			// Read variant/Ind major coding
			if (C_bed.read(&ch) != 1)
				halt_fmt(WISARD_FAIL_READ_BEDSIG);
//				halt("Failed to read signature of BED file!");
			if (ch&1) B_vrtMajor = true;
			else B_vrtMajor = false;

			if (B_vrtMajor)
				LOG("BED file v1.00, variant-major mode\n");
			else
				LOG("BED file v1.00, individual-major mode\n");
		} else v1_bfile = false;
	} else v1_bfile = false;

	/* Incompatible/corrupted BED format */
	if (!v1_bfile)
		halt_fmt(WISARD_INVL_BEDVERSION);

	return B_vrtMajor;
}

size_t cPlinkIO::loadBED_single(cStrFile& C_bed, char B_vrtMajor,
	wsUint* Na_sampMissGeno, vBool& Bv_filtSampOrig, vBool& Bv_filtVrtOrig)
{
	/* If it is variant-major mode, should be transposed when reading*/
	// 00     01   10  11
	// MAJHOM MISS HET MINHOM
	wsUint	i		= 0;
	char	tbl[4]	= { 0, };
// 	xModel	X_model	= _getModel();
// 	switch (X_model) {
// 	case MD_DOMINANT:
// 		tbl[0] = tbl[2] = 1; break;
// 	case MD_RECESSIVE:
// 		tbl[0] = 1; break;
// 	case MD_MULTIPLICATIVE:
// 		tbl[0] = 4;
// 		tbl[2] = 1; break;
// 	default:
		tbl[0] = 2;
		tbl[2] = 1;
// 		break;
// 	}
	char*	Sp_skip	= NULL;
	N_missGeno = 0;
	size_t N_filtGeno = 0;

	if (B_vrtMajor) {
		wsUint I_vrtOrig;
		/* Each byte holds 4 genotypes and variant-major
		 * hence one 'line' should be as below */
		wsUint N_szLine = (N_sampOrig+3)>>2;

		/* Initial data index is already given in _loadBED_FAM() */
		wsAlloc(Sp_skip, char, N_szLine);

		// _i == Index of variants in FILE
		// i  == Index of variants in data, filtered variants is not counted
		for (I_vrtOrig=i=0 ; I_vrtOrig<N_vrtOrig ; I_vrtOrig++) {
			unsigned char C_buf;

			/* Skip all samples if this variant filtered */
			if (Bv_filtVrtOrig[I_vrtOrig] == 1) {
				if (C_bed.read(Sp_skip, (N_sampOrig+3)>>2) != ((N_sampOrig+3)>>2))
					halt_fmt(WISARD_INVL_BEDSIZE, "smaller", N_vrtOrig, N_sampOrig);
				N_filtGeno += N_sampOrig;
				continue;
			}

			// Inner loop for variants
			wsUint I_samp=0, I_sampOrig=0;
			for ( ; I_sampOrig<N_sampOrig ; ) { // for all variants
				//				int tbl[4] = { 0, 0, 1, 2 };
				if (C_bed.read(&C_buf) == 0)
					halt_fmt(WISARD_INVL_BEDSIZE, "smaller", N_vrtOrig, N_sampOrig);

				/* This sample is included */
				if (Bv_filtSampOrig[I_sampOrig] == 0) {
					if ((C_buf&3)==1) {
						N_missGeno++;
						Na_sampMissGeno[I_samp]++;
						Xa_sampleV2[I_samp]->B_isComplete = -1;
						Xv_variant[i].nmis++;
						Na_geno[I_samp++][i] = WISARD_NA;
					} else
						Na_geno[I_samp++][i] = tbl[C_buf&3];
				} else
					N_filtGeno++;
				if (++I_sampOrig == N_sampOrig) break;
				C_buf>>=2;

				/* This sample is included */
				if (Bv_filtSampOrig[I_sampOrig] == 0) { 
					if ((C_buf&3)==1) {
						N_missGeno++;
						Na_sampMissGeno[I_samp]++;
						Xa_sampleV2[I_samp]->B_isComplete = -1;
						Xv_variant[i].nmis++;
						Na_geno[I_samp++][i] = WISARD_NA;
					} else
						Na_geno[I_samp++][i] = tbl[C_buf&3];
				} else
					N_filtGeno++;
				if (++I_sampOrig == N_sampOrig) break;
				C_buf>>=2;

				/* This sample is included */
				if (Bv_filtSampOrig[I_sampOrig] == 0) {
					if ((C_buf&3)==1) {
						N_missGeno++;
						Na_sampMissGeno[I_samp]++;
						Xa_sampleV2[I_samp]->B_isComplete = -1;
						Xv_variant[i].nmis++;
						Na_geno[I_samp++][i] = WISARD_NA;
					} else
						Na_geno[I_samp++][i] = tbl[C_buf&3];  
				} else
					N_filtGeno++;
				if (++I_sampOrig == N_sampOrig) break;
				C_buf>>=2;

				/* This sample is included */
				if (Bv_filtSampOrig[I_sampOrig] == 0) { 
					if ((C_buf&3)==1) {
						N_missGeno++;
						Na_sampMissGeno[I_samp]++;
						Xa_sampleV2[I_samp]->B_isComplete = -1;
						Xv_variant[i].nmis++;
						Na_geno[I_samp++][i] = WISARD_NA;
					} else
						Na_geno[I_samp++][i] = tbl[C_buf&3];
				} else
					N_filtGeno++;
				++I_sampOrig;	
			} /* END OF variants loop */
			i++;

			if ((I_vrtOrig%1000) == 0)
				notice("%d/%d variants retrieved (%u filtered)\r", I_vrtOrig, N_vrtOrig, N_filtGeno);
		}
	} else {
		wsUint I_samp;
		wsUint N_passByte = (wsUint)((Xv_variant.size()-1) >> 2) + 1;
		wsAlloc(Sp_skip, char, N_passByte);

		for (i=I_samp=0 ; i<N_sampOrig ; I_samp++) {
			unsigned char	C_buf;
			wsUint	N_sampMissGeno = 0;

			if (Bv_filtSampOrig[I_samp] == 1) {
				// Size check
				if (C_bed.read(Sp_skip, N_passByte) != N_passByte)
					halt_fmt(WISARD_INVL_BEDSIZE, "smaller", N_vrtOrig, N_sampOrig);
				N_filtGeno += Xv_variant.size();
				continue;
			}

			// Inner loop for variants
			wsUint I_vrt=0, I_vrtOrig=0;
			for ( ; I_vrtOrig<N_vrtOrig ; ) { // for all variants
				if (C_bed.read(&C_buf) == 0)
					halt_fmt(WISARD_INVL_BEDSIZE, "smaller", N_vrtOrig, N_sampOrig);

				if (Bv_filtVrtOrig[I_vrtOrig] == 0) {
					if ((C_buf&3)==1) {
						N_missGeno++;
						N_sampMissGeno++;
						Xa_sampleV2[i]->B_isComplete = -1;
						Xv_variant[I_vrt].nmis++;
						Na_geno[i][I_vrt++] = WISARD_NA;
					} else
						Na_geno[i][I_vrt++] = tbl[C_buf&3];
				}
				if (++I_vrtOrig == N_vrtOrig) break;
				C_buf>>=2;

				if (Bv_filtVrtOrig[I_vrtOrig] == 0) {
					if ((C_buf&3)==1) {
						N_missGeno++;
						N_sampMissGeno++;
						Xa_sampleV2[i]->B_isComplete = -1;
						Xv_variant[I_vrt].nmis++;
						Na_geno[i][I_vrt++] = WISARD_NA;
					} else
						Na_geno[i][I_vrt++] = tbl[C_buf&3];
				}
				if (++I_vrtOrig == N_vrtOrig) break;
				C_buf>>=2;

				if (Bv_filtVrtOrig[I_vrtOrig] == 0) {
					if ((C_buf&3)==1) {
						N_missGeno++;
						N_sampMissGeno++;
						Xa_sampleV2[i]->B_isComplete = -1;
						Xv_variant[I_vrt].nmis++;
						Na_geno[i][I_vrt++] = WISARD_NA;
					} else
						Na_geno[i][I_vrt++] = tbl[C_buf&3];
				}
				if (++I_vrtOrig == N_vrtOrig) break;
				C_buf>>=2;

				if (Bv_filtVrtOrig[I_vrtOrig] == 0) {
					if ((C_buf&3)==1) {
						N_missGeno++;
						N_sampMissGeno++;
						Xa_sampleV2[i]->B_isComplete = -1;
						Xv_variant[I_vrt].nmis++;
						Na_geno[i][I_vrt++] = WISARD_NA;
					} else
						Na_geno[i][I_vrt++] = tbl[C_buf&3];
				}
				I_vrtOrig++;
			}

			/* In default, someone having NO genotype will be excluded */
			/* Filtering sample by --geno */
			wsReal R_gr = W1 - N_sampMissGeno / (wsReal)N_variant;
			if ((N_sampMissGeno == N_variant && !OPT_ENABLED(nodata)) ||
				(IS_ASSIGNED(filgind) &&
					isInRange(OPT_RANGE(filgind), R_gr)) ||
				(IS_ASSIGNED(incgind) &&
					!isInRange(OPT_RANGE(incgind), R_gr))) {
					Bv_filtSamp[i] = 1;
					N_szFilteredSamp++;
			}
			Xa_sampleV2[i]->N_missGeno = N_sampMissGeno;

			/* Status */
			if (((i++)%10) == 0) notice("[%d/%d] samples retrieved\r",
				i, N_sample);
		} /* END OF if-else of Xa_isFiltered */
	}
	DEALLOC(Sp_skip);

	return N_filtGeno;
}

static const unsigned char popcount_tab[] =
{
	0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
	3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
}; 

size_t cPlinkIO::loadBED_single2(cStrFile& C_bed, char B_vrtMajor,
	wsUint* Na_sampMissGeno, vBool& Bv_filtSampOrig, vBool& Bv_filtVrtOrig)
{
	wsAlloc(Na_raw, unsigned char*, N_vrtOrig);

	/* If it is variant-major mode, should be transposed when reading*/
	// 00     01   10  11
	// MAJHOM MISS HET MINHOM
	wsUint	i = 0;
//	char	tbl[4] = { 0, };
	// 	xModel	X_model	= _getModel();
	// 	switch (X_model) {
	// 	case MD_DOMINANT:
	// 		tbl[0] = tbl[2] = 1; break;
	// 	case MD_RECESSIVE:
	// 		tbl[0] = 1; break;
	// 	case MD_MULTIPLICATIVE:
	// 		tbl[0] = 4;
	// 		tbl[2] = 1; break;
	// 	default:
//	tbl[0] = 2;
//	tbl[2] = 1;
	// 		break;
	// 	}
	unsigned char*	Sp_skip = NULL;
	N_missGeno = 0;
	size_t N_filtGeno = 0;
	wsUint N_geno = 0;
	LOG("BED file expected to have %d samples and %d variants\n", N_sampOrig, N_vrtOrig);
	if (B_vrtMajor) {
		/* Initial data index is already given in _loadBED_FAM() */
		wsUint N_byte2read = (N_sampOrig + 3) >> 2;
		wsUint fsIdx;
		wsAlloc(Sp_skip, unsigned char, N_byte2read);

		// _i == Index of variants in FILE
		// i  == Index of variants in data, filtered variants is not counted
		for (fsIdx = i = 0; fsIdx < N_vrtOrig; fsIdx++) {
			// Read #samp data
//			size_t N_read = 
			C_bed.read(Sp_skip, N_byte2read);

			// Inner loop for variants
//			wsUint diIdx = 0, fiIdx = 0;//,rByte=0;

			VARIANT_SSE_PART(N_sampOrig, Sp_skip, s)
//				__m128i tt = *s;
				__m128i s2 = _mm_set_epi32(0xaaaaaaaa, 0xaaaaaaaa,
					0xaaaaaaaa, 0xaaaaaaaa);
				__m128i s1 = _mm_set_epi32(0x55555555, 0x55555555,
					0x55555555, 0x55555555);

				// t1 = (a&2)>>1
				__m128i t1 = _mm_and_si128(*s, s2);
				__m128i u1 = _mm_set_epi32(0, 0, 0, 1);
				t1 = _mm_srl_epi32(t1, u1);

				// s = (a^2) ^ ((a&2)>>1)
				*s = _mm_xor_si128(_mm_xor_si128(*s, s2), t1);
				// m = ((a&2)>>1) _& (a&1)
//				tt = _mm_andnot_si128(t1,
//					_mm_and_si128(*s, s1) );
				
				N_geno += 64;
				// FIXME : Linux compatibility
//				N_missGeno += _mm_popcnt_u64(tt.m128i_i64[0]) + _mm_popcnt_u64(tt.m128i_i64[1]);
			VARIANT_NRM_PART(N_sampOrig, Sp_skip, N)
				unsigned char t1 = (N&0xaa)>>1;
				unsigned char S = (N^2)^t1;
				N_missGeno += popcount_tab[~t1 & (S&0x55)];
				N_geno += 4;
			VARIANT_END_PART(N_sampOrig, Sp_skip, N)
				N_geno += N_pureByte != N_sampOrig ?
					N_sampOrig-N_pureByte : 4;
				// Do same thing
				unsigned char t1 = (N & 0xaa) >> 1;
				unsigned char S = (N ^ 2) ^ t1;
				N_missGeno += popcount_tab[~t1 & (S & 0x55)];

			Na_raw[fsIdx] = Sp_skip;
		}
	}
	LOG("[%d] missing genotypes found, [%d] genotypes found\n", N_missGeno, N_geno);
	DEALLOC(Sp_skip);

	return N_filtGeno;
}

int thr_loadBED(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	if (N_idx < 0 || N_idx >= OPT_NUMBER(thread))
		halt("Error index [%d]", N_idx);
	xThreadBedLoad*	Xp_at		= (xThreadBedLoad *)Vp_shareData;
	wsUint			N_fSamp		= Xp_at->N_fSamp;
	wsUint			N_fileVrt	= Xp_at->N_fileVrt;
	char			B_vrtMajor	= Xp_at->B_vrtMajor;
	xRangedThreadBedLoad*
					Xp_data		= (xRangedThreadBedLoad *)Vp_data;

	vSampPtr&		Xa_sampleV2		= Xp_at->Cp_IO->getSample();
	vVariant&		Xv_variant		= Xp_at->Cp_IO->getVariant();
	char**			Na_data			= Xp_at->Cp_IO->getGenotype();
	wsUint			N_variant		= Xp_at->Cp_IO->sizeVariant();

	wsUint			N_s			= Xp_data->N_start;
	wsUint			N_e			= Xp_data->N_end;
	wsUint			N_jmp		= OPT_NUMBER(thread);
//	wsUint			N_med		= 0;-		tt	{m128i_i8=0x000000d8e496e6c0  <���ڿ��� �߸��� ���ڰ� �ֽ��ϴ�.> m128i_i16=0x000000d8e496e6c0 {0x7df7, 0x7fd7, 0xffff, ...} ...}	__m128i

	wsUint			I;
	pverbose("Thread [%d] computes from [%d] to [%d] with interval [%d]\n",
		N_idx, N_s, N_e, N_jmp);

	/* Init buffer */
	if (!Xp_data->Na_sampMissGeno)
		wsCalloc(Xp_data->Na_sampMissGeno, wsUint, B_vrtMajor ? N_fSamp : N_fileVrt);

	/* For all row */
	wsUint N_szBlock = B_vrtMajor ? (N_fSamp+3)>>2 : (N_fileVrt+3)>>2;
	unsigned char *Sp_ptr = Xp_at->Na_data + N_s * N_szBlock;
	for (I=N_s ; I<=N_e ; I+=N_jmp,Sp_ptr+=N_szBlock) {
		if ((I%10) == 0) notice("[%d/%d] samples processed...\r", I, Xp_at->t.N_size);

		/* If it is variant-major mode, should be transposed when reading*/
		// 00     01   10  11
		// MAJHOM MISS HET MINHOM
		wsUint	i		= 0;
		char	tbl[4]	= { 2, 0, 1, 0 };
//		char*	Sp_skip	= NULL;
//		N_missGeno = 0;
		size_t N_filtGeno = 0;
		if (B_vrtMajor) {
			/* Initial data index is already given in _loadBED_FAM() */
			wsUint fsIdx;

			// _i == Index of variants in FILE
			// i  == Index of variants in data, filtered variants is not counted
			for (fsIdx=i=0 ; fsIdx<N_fileVrt ; fsIdx++) {
				unsigned char C_buf;

				/* Skip all samples if this variant set to be filtered */
				if (Xp_at->Bv_filtVrtOrig[fsIdx] == 1) {
					N_filtGeno += N_fSamp;
					continue;
				}

				// Inner loop for variants
				wsUint diIdx=0,fiIdx=0;//,rByte=0;
				for (;fiIdx<N_fSamp;) { // for all variants
					//				int tbl[4] = { 0, 0, 1, 2 };
					C_buf = Sp_ptr[fiIdx];
// 					if (Xp_at->C_bed.read(&C_buf) == 0)
// 						halt_fmt(WISARD_INVL_BEDSIZE, "smaller", N_fileVrt, N_fSamp);
					/* This sample is included */
					if (Xp_at->Bv_filtSampOrig[fiIdx] == 0) {
						if ((C_buf&3)==1) {
							Xp_data->N_missGeno++;
							Xp_data->Na_sampMissGeno[diIdx]++;
							Xa_sampleV2[diIdx]->B_isComplete = -1;
							Xv_variant[i].nmis++;
							Na_data[diIdx++][i] = WISARD_NA;
						} else
							Na_data[diIdx++][i] = tbl[C_buf&3];
					} else
						N_filtGeno++;
					if (++fiIdx == N_fSamp) break;
					C_buf>>=2;

					/* This sample is included */
					if (Xp_at->Bv_filtSampOrig[fiIdx] == 0) { 
						if ((C_buf&3)==1) {
							Xp_data->N_missGeno++;
							Xp_data->Na_sampMissGeno[diIdx]++;
							Xa_sampleV2[diIdx]->B_isComplete = -1;
							Xv_variant[i].nmis++;
							Na_data[diIdx++][i] = WISARD_NA;
						} else
							Na_data[diIdx++][i] = tbl[C_buf&3];
					} else
						N_filtGeno++;
					if (++fiIdx == N_fSamp) break;
					C_buf>>=2;

					/* This sample is included */
					if (Xp_at->Bv_filtSampOrig[fiIdx] == 0) {
						if ((C_buf&3)==1) {
							Xp_data->N_missGeno++;
							Xp_data->Na_sampMissGeno[diIdx]++;
							Xa_sampleV2[diIdx]->B_isComplete = -1;
							Xv_variant[i].nmis++;
							Na_data[diIdx++][i] = WISARD_NA;
						} else
							Na_data[diIdx++][i] = tbl[C_buf&3];  
					} else
						N_filtGeno++;
					if (++fiIdx == N_fSamp) break;
					C_buf>>=2;

					/* This sample is included */
					if (Xp_at->Bv_filtSampOrig[fiIdx] == 0) { 
						if ((C_buf&3)==1) {
							Xp_data->N_missGeno++;
							Xp_data->Na_sampMissGeno[diIdx]++;
							Xa_sampleV2[diIdx]->B_isComplete = -1;
							Xv_variant[i].nmis++;
							Na_data[diIdx++][i] = WISARD_NA;
						} else
							Na_data[diIdx++][i] = tbl[C_buf&3];
					} else
						N_filtGeno++;
					++fiIdx;	
				}
				// 				if (rByte != ((N_fSamp+3)>>2))
				// 					halt("SYSERR: Insufficient read [%d read, expected to read %d]", rByte, ((N_fSamp+3)>>2));

				i++;

				if ((fsIdx%1000) == 0)
					notice("%d/%d variants retrieved (%u filtered)\r", fsIdx, N_fileVrt, N_filtGeno);
			}
		} else {
			wsUint diIdx;
// 			wsUint N_passByte = (wsUint)((Xv_variant.size()-1) >> 2) + 1;
// 			MULTI_MALLOC(Sp_skip, char, N_passByte);

			for (i=diIdx=0 ; i<N_fSamp ; diIdx++) {
				unsigned char	C_buf;
				wsUint	N_sampMissGeno = 0;

				// Inner loop for variants
				if (Xp_at->Bv_filtSampOrig[diIdx] == 1) {
// 					if (Xp_at->C_bed.read(Sp_skip, N_passByte) != N_passByte)
// 						halt_fmt(WISARD_INVL_BEDSIZE, "smaller", N_fileVrt, N_fSamp);
					//halt("%d byte required to pass but failed", N_passByte);
					N_filtGeno += Xv_variant.size();
				} else {
					wsUint dsIdx=0, fsIdx=0;
					while (fsIdx<N_fileVrt) { // for all variants
						C_buf = Sp_ptr[fsIdx];
// 						if (Xp_at->C_bed.read(&C_buf) == 0)
// 							halt_fmt(WISARD_INVL_BEDSIZE, "smaller", N_fileVrt, N_fSamp);
						//						halt("BED size is smaller than expected, possibly error?");

						if (Xp_at->Bv_filtVrtOrig[fsIdx] == 0) {
							if ((C_buf&3)==1) {
								Xp_data->N_missGeno++;
								N_sampMissGeno++;
								Xa_sampleV2[i]->B_isComplete = -1;
								Xv_variant[dsIdx].nmis++;
								Na_data[i][dsIdx++] = WISARD_NA;
							} else
								Na_data[i][dsIdx++] = tbl[C_buf&3];
						}
						if (++fsIdx == N_fileVrt) break;
						C_buf>>=2;

						if (Xp_at->Bv_filtVrtOrig[fsIdx] == 0) {
							if ((C_buf&3)==1) {
								Xp_data->N_missGeno++;
								N_sampMissGeno++;
								Xa_sampleV2[i]->B_isComplete = -1;
								Xv_variant[dsIdx].nmis++;
								Na_data[i][dsIdx++] = WISARD_NA;
							} else
								Na_data[i][dsIdx++] = tbl[C_buf&3];
						}
						if (++fsIdx == N_fileVrt) break;
						C_buf>>=2;

						if (Xp_at->Bv_filtVrtOrig[fsIdx] == 0) {
							if ((C_buf&3)==1) {
								Xp_data->N_missGeno++;
								N_sampMissGeno++;
								Xa_sampleV2[i]->B_isComplete = -1;
								Xv_variant[dsIdx].nmis++;
								Na_data[i][dsIdx++] = WISARD_NA;
							} else
								Na_data[i][dsIdx++] = tbl[C_buf&3];
						}
						if (++fsIdx == N_fileVrt) break;
						C_buf>>=2;

						if (Xp_at->Bv_filtVrtOrig[fsIdx] == 0) {
							if ((C_buf&3)==1) {
								Xp_data->N_missGeno++;
								N_sampMissGeno++;
								Xa_sampleV2[i]->B_isComplete = -1;
								Xv_variant[dsIdx].nmis++;
								Na_data[i][dsIdx++] = WISARD_NA;
							} else
								Na_data[i][dsIdx++] = tbl[C_buf&3];
						}
						fsIdx++;
					}

					/* In default, someone having NO genotype will be excluded */
					/* Filtering sample by --geno */
					wsReal R_gr = W1 - N_sampMissGeno / (wsReal)N_variant;
					if (//N_sampMissGeno == N_variant ||
						(IS_ASSIGNED(filgind) &&
							isInRange(OPT_RANGE(filgind), R_gr)) ||
						(IS_ASSIGNED(incgind) &&
							!isInRange(OPT_RANGE(incgind), R_gr))) {
							Xp_at->Bv_filtSamp[i] = 1;
							Xp_data->N_szFilteredSamp++;
					}

					/* Status */
					if (((i++)%10) == 0) notice("%d/%d samples retrieved\r", i, Xp_at->t.N_size);
				}
			} /* END OF if-else of Xa_isFiltered */
		}
	}

	return 0;
}

void cPlinkIO::_loadBED(char *S_bedIn)
{
	char	S_bedFn[MAX_PATH];
	char	*S_buf = NULL;
	wsAlloc(S_buf, char, 4096);

	strcpy(S_bedFn, S_bedIn);
	/* Prepare handle */
	cStrFile C_bed(S_bedFn, "Binary PED file", 1, 1);
	if (C_bed.isFailed()) {
		/* try to add .bed and again */
		strcat(S_bedFn, ".bed");
		C_bed.open(S_bedFn, "Binary PED file", 0, 1);
	}

	/* Load BIM file */
	_loadBED_BIM(S_bedFn);

	/* Load FAM file */
	_loadBED_FAM(S_bedFn);

	/* Copy original filtering buffer */
	vBool Bv_filtSampOrig;
	Bv_filtSampOrig.resize(Bv_filtSamp.size());
	copy(Bv_filtSamp.begin(), Bv_filtSamp.end(),
		Bv_filtSampOrig.begin());
	vBool Bv_filtVrtOrig;
	Bv_filtVrtOrig.resize(Bv_filtVrt.size());
	copy(Bv_filtVrt.begin(), Bv_filtVrt.end(),
		Bv_filtVrtOrig.begin());

	/* Resize current buffer and clear */
	Bv_filtSamp.clear();
	Bv_filtVrt.clear();
	Bv_filtSamp.resize(N_sample, 0);
	Bv_filtVrt.resize(N_variant, 0);

	/* Find out the type of BED file */
	size_t	N_filtGeno			= 0;
	wsUint	*Na_sampMissGeno	= NULL;
	char	B_vrtMajor			= _loadBED_header(C_bed);
	if (B_vrtMajor) {
		/* If this BED file is variant-major mode, additional space is required for counting available variants for all samples */
		wsCalloc(Na_sampMissGeno, wsUint, N_sample);
	}

	/* Try to allocate file covers entire BED file */
#ifdef USE_MTLOAD
	size_t N_szBed = B_vrtMajor ? ((N_sampOrig+3)>>2) * N_vrtOrig :
		((N_vrtOrig+3)>>2) * N_sampOrig;
	unsigned char *Sa_bed = new unsigned char[N_szBed];
	if (Sa_bed && IS_MULTITHREAD) {
		C_bed.read(Sa_bed, N_szBed);
		LOGnote("Sufficient memory found, data will be loaded multithreaded manner...\n");
		xThreadBedLoad X = {
			{ B_vrtMajor ? N_vrtOrig : N_sampOrig },
			N_sampOrig,
			N_vrtOrig,
			this,
			B_vrtMajor,
			Sa_bed,
			Bv_filtVrtOrig,
			Bv_filtSampOrig,
			Bv_filtSamp
		};
		/* Remove buffer */
		WORKER().run(thr_loadBED, distEqual, &X, NULL, sizeof(xRangedThreadBedLoad));
		delete [] Sa_bed;
	} else
#endif
		N_filtGeno = loadBED_single(C_bed, B_vrtMajor, Na_sampMissGeno,
			Bv_filtSampOrig, Bv_filtVrtOrig);

	/* Check EOF */ {
		unsigned char C_buf;
		if (C_bed.read(&C_buf))
			halt_fmt(WISARD_INVL_BEDSIZE, "larger", N_vrtOrig, N_sampOrig);
	}

	/* Close dataset file */
	C_bed.close();
	LOG("Genotype loaded [BED format] [%u filtered]\n", N_filtGeno);

	/* If this BED file is variant-major, check Na_sampMissGeno
	 * to determine which sample should be excluded */
	if (B_vrtMajor) {
		for (wsUint i=0 ; i<N_sample ; i++) {
			wsReal R_gr = W1 - Na_sampMissGeno[i] / (wsReal)N_variant;
			/* In default, someone having NO genotype will be excluded
			 * Otherwise, it is also possible to exclude samples with
			 *	--filgind option */
			if ((Na_sampMissGeno[i] == N_variant && !OPT_ENABLED(nodata)) ||
				(IS_ASSIGNED(filgind) && isInRange(OPT_RANGE(filgind), R_gr)) ||
				(IS_ASSIGNED(incgind) && !isInRange(OPT_RANGE(incgind), R_gr))) {
				Bv_filtSamp[i] = 1;
				N_szFilteredSamp++;
			}
			Xa_sampleV2[i]->N_missGeno = Na_sampMissGeno[i];
		}
		LOG("%d samples filtered by genotyping rate filtering\n",
			N_szFilteredSamp);
	}
	DEALLOC(Na_sampMissGeno);

	DEALLOC(S_buf);
}

void cPlinkIO::_loadPED_countVariant(char *S_fn, char B_isRawFile, char *Sp_tmp1)
{
	char	B_wereBuffer = Sp_tmp1 != NULL;
	char	*Sp_tmp2 = NULL;

	/* Counting variant */
	N_variant = N_vrtOrig = 0;
	if (B_wereBuffer) {
		if (B_isRawFile) for (Sp_tmp2=Sp_tmp1 ; ; N_vrtOrig++) {
			/* Get next token */
			getString(&Sp_tmp2, &Sp_tmp1);

			/* RAW file contains its variant information in the first line of RAW file */
			if (N_vrtOrig > 4) {
				xVariant	X_currVrt = {0, };
				char	B_filtered = 0;

				/* For convenience, major allele=1 and minor allele=0 */
				X_currVrt.al1		= '0';
				X_currVrt.al2		= '1';
				X_currVrt.indel1	= NULL;
				X_currVrt.indel2	= NULL;
				X_currVrt.filter	= 0;
				X_currVrt.gdist	= W0;
				/* Chromosome number is undeterministic, because RAW file does not contain such information */
				X_currVrt.chr		= 0;	
				X_currVrt.pos		= 0;

				/* Retrieve variant name */
				wsAlloc(X_currVrt.name, char, strlen(Sp_tmp2)+1);
				strcpy(X_currVrt.name, Sp_tmp2);

				/* [FILTERING]
				 * Variant filtering by other options */
				B_filtered = _filterVariantPosition(X_currVrt);

				if (B_filtered == 1)
					Bv_filtVrt.push_back(1);
				else {
					Bv_filtVrt.push_back(0);
					Xv_variant.push_back(X_currVrt);
					N_variant++;
				}
			}

			if (Sp_tmp1 == NULL)
				break;
			Sp_tmp2 = Sp_tmp1;
		} else {
			/* If --sepallele */
			char S_delim = IS_ASSIGNED(sepallele) ? OPT_STRING(sepallele)[0] : ' ';
			for (Sp_tmp2=Sp_tmp1 ; ; N_vrtOrig++) {
				/* Get next token */
				getStringDelim(&Sp_tmp2, &Sp_tmp1, S_delim);
				if (Sp_tmp1 == NULL)
					break;
				Sp_tmp2 = Sp_tmp1;
			}
		}
		lverbose("Number of columns found [%d], sex [%d] fid [%d] pheno [%d] parent [%d]\n",
			N_vrtOrig, OPT_ENABLED(nosex), (OPT_ENABLED(nofid) || IS_ASSIGNED(sepid)), OPT_ENABLED(nopheno), (OPT_ENABLED(noparent)<<1));
		N_vrtOrig -= 4 - OPT_ENABLED(nosex) - (OPT_ENABLED(nofid) || IS_ASSIGNED(sepid)) -
			OPT_ENABLED(nopheno) - (OPT_ENABLED(noparent)<<1);

		/* Checking variant count if it is PED file, it should be multiply of two */
		if (B_isRawFile == false) {
			if ((N_vrtOrig%2) == 1)
				halt("Allele count is not the even!");

			/* 2 allele consists 1 variant */
			N_vrtOrig >>= 1;
		}
	}

	/* PED file have MAP file, which contains variant information */
	if ((!B_isRawFile && N_vrtOrig) || (!B_isRawFile && !B_wereBuffer)) {

		/* --varresize filter */
		char* Ba_isVrtSelected = NULL;
		if (IS_ASSIGNED(varresize)) {
			wsReal R_valSS = OPT_REAL(varresize);
			wsUint N_sel = R_valSS >= W1 ?
				(wsUint)R_valSS : (wsUint)round(R_valSS*N_vrtOrig);
			wsCalloc(Ba_isVrtSelected, char, N_vrtOrig);

			/* Check it should be resized */
			if (N_vrtOrig <= N_sel) {
				LOG("# of variants in dataset[%d] is larger than the number of "
					"variants to be resized[%d] so the dataset will not be resized\n",
					N_vrtOrig, N_sel);
			} else {
				/* Randomly select variants */
				for (wsUint N_chosen = 0; N_chosen < N_sel;) {
					/* Select one */
					wsUint N_curr = rand() % N_vrtOrig;

					if (Ba_isVrtSelected[N_curr] == 0) {
						Ba_isVrtSelected[N_curr] = 1;
						N_chosen++;
					}
				}
			}
		}

		/* Find . from backward */
		char Sp_dupFn[1024];

		/* If --nomap assigned, just make */
		if (OPT_ENABLED(nomap)) {
			/* --nomap & --lgen exclusive mutually */
			if (IS_ASSIGNED(lgen))
				halt_fmt(WISARD_CANT_OPT_W_OTHEROPT, "--nomap", "--lgen");

			LOGwarn("--nomap found, [%d] variants will named sequentially, "
				"and other will be set to default\n", N_vrtOrig);

			if (IS_ASSIGNED(map))
				halt_fmt(WISARD_CANT_OPT_W_OTHEROPT, "--nomap", "--map");
//				halt("--nomap and --map cannot assigned concurrently!");
			for (wsUint i=0 ; i<N_vrtOrig ; i++) {
				xVariant X_currVrt = {0, };

				X_currVrt.chr = 0;
				wsAlloc(X_currVrt.name, char, 64);
				sprintf(X_currVrt.name, "VARIANT_%d", i+1);
				X_currVrt.al1		= '1';
				X_currVrt.al2		= '2';
				X_currVrt.indel1	= NULL;
				X_currVrt.indel2	= NULL;
				X_currVrt.pos		= 0;
				X_currVrt.filter	= 0;
				X_currVrt.gdist		= W0;

				/* Insert anyway */
				Xv_variant.push_back(X_currVrt);
				Bv_filtVrt.push_back(0);
				N_variant++;
			}
		} else {
			/* Make MAP filename */
			if (IS_ASSIGNED(map))
				strcpy(Sp_dupFn, OPT_STRING(map));
			else
				otherExt(S_fn, "map", Sp_dupFn);

			/* Open MAP file */
			cStrFile C_map(Sp_dupFn, "MAP file");
// 			FILE *H_map = fopen(Sp_dupFn, "r");
// 			if (H_map == NULL) {
// 				halt("Can't find MAP file of [%s]", S_fn);
// 			}

			/* Read */
			wsUint	N_vrtInMAP	= 0;
			char	*Sp_buf		= NULL;
			wsAlloc(Sp_buf, char, 1024);

			for ( ; C_map.gets(Sp_buf, 1024) ; N_vrtInMAP++) {
				if (Sp_buf[0] == '\r' || Sp_buf[0] == '\n' || Sp_buf[0] == '\0')
					break;
				xVariant X_currVrt = {0, };
				char B_filtered = 0;

				/* S_buf == chromosome */
				Sp_tmp1 = NULL;
				getString(&Sp_buf, &Sp_tmp1);
				if (OPT_ENABLED(autoonly) || OPT_ENABLED(sexonly)) {
					X_currVrt.chr = getChr(Sp_buf);

					/* If this variant should be skipped */
					if (X_currVrt.chr == -1 ||
						(OPT_ENABLED(autoonly) && X_currVrt.chr > (int)NAUTO_SPECIES) ||
						(OPT_ENABLED(sexonly) && X_currVrt.chr <= (int)NAUTO_SPECIES))
						B_filtered = 1;
				} else
					X_currVrt.chr = getChr(Sp_buf);

				/* a == variant name */
				if (!Sp_tmp1) halt_fmt(WISARD_INVL_FILE_INCMPL, "MAP", C_map.getName(), "Unexpected end of line");
				getString(&Sp_tmp1, &Sp_tmp2);
				wsAlloc(X_currVrt.name, char, strlen(Sp_tmp1)+1);
				strcpy(X_currVrt.name, Sp_tmp1);

				/* Sp_tmp2 == gdist */
				if (!OPT_ENABLED(nogdist)) {
					getString(&Sp_tmp2, &Sp_tmp1);
					X_currVrt.gdist		= (wsReal)atof(Sp_tmp2);
				} else
					X_currVrt.gdist		= W0;

				/* Sp_tmp1 == pos */
				if (!OPT_ENABLED(nopos)) {
					getString(&Sp_tmp1, &Sp_tmp2);

					/* Position validity check */
					for (char *p=Sp_tmp1 ; *p ; p++)
						if (*p < '0' || *p > '9')
							halt_fmt(WISARD_INVL_VARIANTPOS, X_currVrt.name, Sp_tmp1);
							//halt("variant %s have invalid position %s", X_currVrt.name, Sp_tmp1);
					int N_ret = atoi(Sp_tmp1);
					if (N_ret < 0) B_filtered = 1;
					X_currVrt.pos		= (wsUint)N_ret;
				} else
					X_currVrt.pos		= 0;

				/* Filter if not chosen */
				if (Ba_isVrtSelected && !Ba_isVrtSelected[N_vrtInMAP])
					B_filtered = 1;

				/* Check variant filtering conditions */
				B_filtered |= _filterVariantPosition(X_currVrt);

				X_currVrt.indel1	= NULL;
				X_currVrt.indel2	= NULL;
				X_currVrt.filter	= 0;

				if (B_filtered == 1) {
					Bv_filtVrt.push_back(1);
					DEALLOC(X_currVrt.name);
				} else {
					Bv_filtVrt.push_back(0);

					/* al1 and al2 should be determined */
					Xv_variant.push_back(X_currVrt);
					N_variant++;
				}
			}
			//fclose(H_map);
			/* Check variant count between MAP and PED */
			if (B_wereBuffer) {
				if (N_vrtInMAP != N_variant) {
					DEALLOC(Sp_buf);
					halt_fmt(WISARD_INVL_NVRTMAP_W_NVRTPED, N_vrtInMAP, N_variant);
					//halt("Variant count in MAP file [%d] is different in variant count in PED file [%d]");
				}
				N_vrtOrig = N_variant;
			} else
				N_vrtOrig = N_vrtInMAP;

			/* Free up duplicate of filename */
			DEALLOC(Sp_buf);
		}
		DEALLOC(Ba_isVrtSelected);
	} else
		/* If RAW file, MAP file is not exist, so pass */
		N_variant = N_vrtOrig;

	/* --varresize filter */
	if (IS_ASSIGNED(varresize)) {
		wsReal R_valSS = OPT_REAL(varresize);
		wsUint N_sel = R_valSS >= W1 ?
			(wsUint)R_valSS : (wsUint)round(R_valSS*N_vrtOrig);

		/* Check it should be resized */
		if (N_variant <= N_sel) {
			LOG("# of loaded variants in dataset[%d] is larger than the number of "
				"variants to be resized[%d] so the dataset will not be resized\n",
				N_variant, N_sel);
		}
		else {
			/* Randomly select variants */
			for (wsUint N_chosen = 0; N_chosen < N_sel;) {
				/* Select one */
				wsUint N_curr = rand() % N_vrtOrig;

				if (Bv_filtVrt[N_curr] == 0) {
					Bv_filtVrt[N_curr] = 1;
					N_chosen++;
				}
			}
		}
	}

	LOG("Variant loaded [MAP format] [%d/%d] chosen\n", N_variant, N_vrtOrig);
}

char cPlinkIO::_loadPED_init(char *S_fn, char *Sp_buf, cStrFile &C_file,
	wsUint *Np_dataPos, char **Sp_tmp1)
{
	char	B_isRawFile	= false;

	/* Read first row to [Sp_buf]
	 * To confirm the type of file : PED or RAW
	 */
	C_file.open(S_fn, "PED/RAW file", 1);
	if (C_file.isFailed()) {
		/* try to add .ped and again */
		strcat(S_fn, ".ped");
		C_file.open(S_fn, "PED/RAW file", 0);
	}

	/* --nskip */
	wsUint N_actualDtPos = 0;
	wsUint N_skip = OPT_NUMBER(nskip);
	if (IS_ASSIGNED(nskip)) {
		for (wsUint i=0 ; i<N_skip ; i++)
			if (C_file.gets(Sp_buf, PED_MAXLEN) == 0)
				halt_fmt(WISARD_NULL_PEDCONTENTS, S_fn);
		N_actualDtPos = (wsUint)C_file.tell();
	}
	/* Check there is contents */
	if (C_file.gets(Sp_buf, PED_MAXLEN) == 0)
		halt_fmt(WISARD_NULL_PEDCONTENTS, S_fn);
	*Np_dataPos = (wsUint)C_file.tell();

	/* Check file type */
	getString(&Sp_buf, Sp_tmp1);
	if (!strcmp(Sp_buf, "FID") ||
		((OPT_ENABLED(nofid) || IS_ASSIGNED(sepid)) && !strcmp(Sp_buf, "IID"))) { /* Raw file have its first column to 'FID' */
		LOG("Input file is RAW file of PLINK\n");
		LOG("All variants will treated as autosomal variants\n");

		/* RAW file constraint check */
		if (IS_ASSIGNED(incrange) || IS_ASSIGNED(filrange))
			LOG("--incrange and --filrange will be ignored since RAW file "
				"does not have range info\n");

		B_isRawFile = true;
	} else
		*Np_dataPos = N_actualDtPos; /* Data starts from 0 when it is PED file */

	return B_isRawFile;
}

void cPlinkIO::_loadPED_asAdditive(char *Ba_charFound)
{
	/* 1 : Have 1 */
	if (Ba_charFound[(int)'1']) {
		char mask[128] = {
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
		};

		for (int i=0 ; i<128 ; i++)
			if (Ba_charFound[i]&mask[i])
				halt_fmt(WISARD_INVL_PEDCHAR, OPT_STRING(ped), "1/2/3/4", "non-permitted characters", i);
//				halt("This PED file seems 1/2/3/4 coding, but non-permitted character were found\n");

		/* If '2' or '3' exists AT LEAST ONCE */
		if (Ba_charFound[(int)'2'] || Ba_charFound[(int)'3']) {
			/* 0/N cannot exists simultaneously */
			if (Ba_charFound[(int)'0'] && Ba_charFound[(int)'N'])
				halt_fmt(WISARD_INVL_PEDMISCHAR, OPT_STRING(ped), "1/2/3/4", "multiple missing characters", "3 and N");
//				halt("This PED file seems 1/2/3/4 coding, but multiple missing character[3 and N] were found\n");

			if (Ba_charFound[(int)'0'])
				C_pedMissingChar = '0';
			else if (Ba_charFound[(int)'N'])
				C_pedMissingChar = 'N';
			else
				C_pedMissingChar = '-';

			LOG(" ::: PED file, 1234 coding, Missing char %c :::\n", C_pedMissingChar);
		} else {
			/* Otherwise, it is 0/1 coding */
			if (Ba_charFound[(int)'N'])
				C_pedMissingChar = 'N';
			else
				C_pedMissingChar = '-';

			LOG(" ::: PED file, 0/1 coding, Missing char %c :::\n", C_pedMissingChar);
		}
	} else if (Ba_charFound[(int)'D'] || Ba_charFound[(int)'I']) {
		/* [3] : A/C/G/T D/I coding */
		char mask[128] = {
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,0,1,0,0,1,1,0,1,0,1,1,1,1,0,1,
			1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
		};

		for (int i=0 ; i<128 ; i++)
			if (Ba_charFound[i]&mask[i])
				halt_fmt(WISARD_INVL_PEDCHAR, OPT_STRING(ped), "A/C/G/T & D/I", "non-permitted characters", i);
		//				halt("This PED file seems A/C/G/T coding, but non-permitted character (%c[%d]) were found", (char)i, i);

		/* 0/3, 3/N, 0/N cannot exists simultaneously */
		if ((Ba_charFound[(int)'0'] && Ba_charFound[(int)'3']) ||
			(Ba_charFound[(int)'0'] && Ba_charFound[(int)'N']) ||
			(Ba_charFound[(int)'N'] && Ba_charFound[(int)'3']))
			halt_fmt(WISARD_INVL_PEDMISCHAR, OPT_STRING(ped), "A/C/G/T & D/I", "multiple missing characters", "at least 2 of 0/3/N");

		if (Ba_charFound[(int)'0'])
			C_pedMissingChar = '0';
		else if (Ba_charFound[(int)'3'])
			C_pedMissingChar = '3';
		else if (Ba_charFound[(int)'N'])
			C_pedMissingChar = 'N';
		else {
			if (OPT_ENABLED(indel))
				C_pedMissingChar = '0';
			else
				C_pedMissingChar = '-';
		}

		LOG(" ::: PED file, ACGT & DI coding, Missing char %c :::\n", C_pedMissingChar);
	} else {
		/* [3] : A/C/G/T coding */
		char mask[128] = {
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,0,1,0,1,1,1,0,1,1,1,1,1,1,0,1,
			1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
			1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
		};

		for (int i=0 ; i<128 ; i++)
			if (Ba_charFound[i]&mask[i])
				halt_fmt(WISARD_INVL_PEDCHAR, OPT_STRING(ped), "A/C/G/T", "non-permitted characters", i);
//				halt("This PED file seems A/C/G/T coding, but non-permitted character (%c[%d]) were found", (char)i, i);

		/* 0/3, 3/N, 0/N cannot exists simultaneously */
		if ((Ba_charFound[(int)'0'] && Ba_charFound[(int)'3']) ||
			(Ba_charFound[(int)'0'] && Ba_charFound[(int)'N']) ||
			(Ba_charFound[(int)'N'] && Ba_charFound[(int)'3']))
			halt_fmt(WISARD_INVL_PEDMISCHAR, OPT_STRING(ped), "A/C/G/T", "multiple missing characters", "at least 2 of 0/3/N");

		if (Ba_charFound[(int)'0'])
			C_pedMissingChar = '0';
		else if (Ba_charFound[(int)'3'])
			C_pedMissingChar = '3';
		else if (Ba_charFound[(int)'N'])
			C_pedMissingChar = 'N';
		else {
			if (OPT_ENABLED(indel))
				C_pedMissingChar = '0';
			else
				C_pedMissingChar = '-';
		}

		LOG(" ::: PED file, ACGT coding, Missing char %c :::\n", C_pedMissingChar);
	}
	/* Check that this PED file does not have missing genotype */
	if (C_pedMissingChar == '-')
		B_isDataComplete = 1;

	/* Additive re-coding */
	/* For debug */
	int N_thread = OPT_NUMBER(thread);
	int *Na_missGeno = NULL;
	wsCalloc(Na_missGeno, int, N_thread);
	WORKER().run(calcMAFv2, forAllVariant_equal, this, Na_missGeno, sizeof(int)*3);
	notice("%d variants converted\r", N_variant);

	/* Sum up calc */
	N_missGeno = 0;
	for (int i=0 ; i<N_thread ; i++)
		N_missGeno += Na_missGeno[i];
	delete [] Na_missGeno;

	/* Deallocate expired memories */
	deallocMatrix(Na_charData, N_sample);
	Na_charData = NULL;
}

void cPlinkIO::_loadPED(char *S_fn)
{
	wsUint		N_dataStPos	= 0;
	char		B_isRawFile	= 0;
	char		*Sp_buf			= NULL;
	char		*Sp_tmp1		= NULL;
	cStrFile	C_file;
	char*		S_finalFn		= NULL;
	wsAlloc(S_finalFn, char, strlen(S_fn) + 16);
	strcpy(S_finalFn, S_fn);
	wsReal		*Rp_pheno;					///< Phenotype in PED/RAW file

	/* Allocate buffer for read PED file */
	wsAlloc(Sp_buf, char, PED_MAXLEN);

	/* Initialize file stream */
	B_isRawFile = _loadPED_init(S_finalFn, Sp_buf, C_file, &N_dataStPos,
		&Sp_tmp1);
	if (Sp_tmp1 == NULL)
		halt_fmt(WISARD_INVL_PLINK_INPUT, S_finalFn);

	/* --indel,--acgt cannot be used when input is raw file */
	if (B_isRawFile) {
		if (OPT_ENABLED(indel))
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--indel", "input is PLINK RAW format");
		if (Na_ACGT)
			halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--acgt", "input is PLINK RAW format");
	}
		//halt_fmt(WISARD_CANT_INDEL_WITH_RAW);

	/* Get variant count */
	_loadPED_countVariant(S_finalFn, B_isRawFile, Sp_tmp1);
	DEALLOC(S_finalFn);

	/* Go the to genotype data start position */
	C_file.seek(N_dataStPos, SEEK_SET);

	/* Get sample count */
	for (N_sample=N_sampOrig=0 ; C_file.gets(Sp_buf, PED_MAXLEN) ; N_sampOrig++) {
		pverbose("[%d/%d] samples chosen/found...\r", N_sample, N_sampOrig);
		_loadRealTimeFamData(Sp_buf, N_sampOrig);
	}
	LOG("Sample loaded [PED format] [%d/%d] chosen\n", N_sample, N_sampOrig);

	/* Make complete family data */
	_registMissingFounder();

	/* Copy original filtering buffer */
	vBool Ba_isSampFiltOrig;
	Ba_isSampFiltOrig.resize(Bv_filtSamp.size());
	copy(Bv_filtSamp.begin(), Bv_filtSamp.end(),
		Ba_isSampFiltOrig.begin());
	vBool Ba_isVrtFiltOrig;
	Ba_isVrtFiltOrig.resize(Bv_filtVrt.size());
	copy(Bv_filtVrt.begin(), Bv_filtVrt.end(),
		Ba_isVrtFiltOrig.begin());

	/* Resize current buffer and clear */
	Bv_filtSamp.clear();
	Bv_filtVrt.clear();
	Bv_filtSamp.resize(N_sample, 0);
	Bv_filtVrt.resize(N_variant, 0);

	/* Allocate required memory */
	_initMemories();
	Rp_pheno	= Ra_pheno[0];

	/* If it is PED file, extra space is required to calculate MA */
	if (!B_isRawFile) {
		wsAlloc(Na_charData, USHORT_t*, N_sample);
		for (wsUint i=0 ; i<N_sample ; i++) {
			wsAlloc(Na_charData[i], USHORT_t, N_variant);
		}
	}

	/* Read */
	wsUint diIdx=0;
	char Ba_chrFound[128] = {0, };
	map<string, int> Xa_family;

	/* Seek to the start of genotype data */
	C_file.seek(N_dataStPos, SEEK_SET);

	/* Load PED file main */
	cExporter *Cp_misInds = NULL;
	for (wsUint fiIdx=0 ; C_file.gets(Sp_buf, PED_MAXLEN) ; fiIdx++) {
		if (OPT_ENABLED(indel)) {
			if (_loadPED_INDELmain(Sp_buf, fiIdx, diIdx, B_isRawFile,
				Rp_pheno, Ba_isSampFiltOrig, Ba_isVrtFiltOrig, Ba_chrFound,
				&Cp_misInds))
				continue;
		} else {
			if (_loadPED_SNVmain(Sp_buf, fiIdx, diIdx, B_isRawFile, Rp_pheno,
				Ba_isSampFiltOrig, Ba_isVrtFiltOrig, Ba_chrFound,
				&Cp_misInds))
				continue;
		}

		/* Increase data index */
		diIdx++;
		if ((diIdx%10) == 0)
			notice("%d samples retrieved\r", diIdx);
	}
	if (Cp_misInds)
		delete Cp_misInds;
	C_file.close();
	DEALLOC(Sp_buf);
	LOG("%d samples retrieved\n", diIdx);

	/* PED file requires minor-allele determination */
	if (!B_isRawFile) {
		_loadPED_asAdditive(Ba_chrFound);

		wsUint i;
/**/	cExporter *Cp_exp = NULL;
		for (i=0 ; i<N_variant ; i++)
			if (Xv_variant[i].filter) {
/**/			Cp_exp = cExporter::summon("allelic.excl.lst");
				break;
			}
		for ( ; i<N_variant ; i++)
			if (Xv_variant[i].filter) {
				Cp_exp->fmt("%s\n", Xv_variant[i].name);
				Bv_filtVrt[i] = 1;
				N_szFilteredVrt++;
			}
		if (Cp_exp) delete Cp_exp;
		if (N_szFilteredVrt)
			LOGoutput("[%d] variants will be excluded because they have "
				"more than two alleles, the list is exported to "
				"[%s.allelic.excl.lst]\n", N_szFilteredVrt, OPT_STRING(out));

		/* In default, someone having NO genotype will be excluded */
		/* Filtering by --geno */
		for (wsUint i=0 ; i<N_sample ; i++) {
			if (Bv_filtSamp[i]) continue;

			wsReal R_gr = W1 - Xa_sampleV2[i]->N_missGeno / (wsReal)N_variant;
			if ((Xa_sampleV2[i]->N_missGeno == (int)N_variant && !OPT_ENABLED(nodata))
				|| (IS_ASSIGNED(filgind) &&
					isInRange(OPT_RANGE(filgind), R_gr))
				|| (IS_ASSIGNED(incgind) &&
					!isInRange(OPT_RANGE(incgind), R_gr))
				) {
				Bv_filtSamp[i] = 1;
				N_szFilteredSamp++;
			}
		}
	}
}

void cPlinkIO::_loadTPED(char *S_fn)
{
	/* W/o --indel, --misgeno should be 1 byte */
	if (!OPT_ENABLED(indel) && OPT_STRING(misgeno)[1])
		halt_fmt(WISARD_CANT_DO_WO_SOMEOPT, "--misgeno with multiple characters", "--indel");
//		halt_fmt(WISARD_CANT_LONGMISGENO_WO_INDEL);

	char	*Sp_buf		= NULL;
//	char	*Sp_tmp1	= NULL;
	FILE*	H_fpFam		= NULL;
//	wsReal	*Rp_pheno	= NULL;					///< Phenotype in PED/RAW file

	/* --nomap is not allowed in TPED loading */
	if (OPT_ENABLED(nomap))
		halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--nomap", "input is TPED format");
//		halt_fmt(WISARD_CANT_NOMAP_WITH_TPED);

	/* Try to open TPED file */
	cStrFile C_file(S_fn, "TPED dataset");

	/* Allocate buffer for read TPED file */
	wsAlloc(Sp_buf, char, PED_MAXLEN);

	char S_fnTFAM[512] = { 0, };
	/* Load TFAM file first */
	if (IS_ASSIGNED(fam)) {
		strcpy(S_fnTFAM, OPT_STRING(fam));
		H_fpFam = fopen(S_fnTFAM, "r");
		if (H_fpFam == NULL)
			halt_fmt(WISARD_FAIL_OPEN_FILEWITHDESC, "TFAM", S_fnTFAM);
	} else {
		/* Find extension part */
		char *Sp_ext = S_fn + strlen(S_fn) - 1;
		while (*Sp_ext!='.' && Sp_ext>S_fn) Sp_ext--;

		/* If there is NO extension part, TFAM file cannot find */
		if (S_fn == Sp_ext)
			halt_fmt(WISARD_FAIL_AUTOOPEN_TFAMEXTMISS, S_fn);

		/* Make TFAM extension */
		*Sp_ext = '\0';
		sprintf(S_fnTFAM, "%s.tfam", S_fn);
		
		/* Try to open */
		H_fpFam = fopen(S_fnTFAM, "r");
		if (H_fpFam == NULL) {
			/* Try to open .fam */
			sprintf(S_fnTFAM, "%s.fam", S_fn);
			H_fpFam = fopen(S_fnTFAM, "r");

			/* Failed to open both of .fam / .tfam file */
			if (H_fpFam == NULL) {
				*Sp_ext = '.';
				halt_fmt(WISARD_FAIL_AUTOOPEN_TFAMNOTEXIST, S_fn);
			}
		}
		LOG("Corresponding TFAM file [%s] found\n", S_fnTFAM);
	}

	/* OK, read this */
	N_sample = 0;
	for (N_sampOrig=0 ; fgets(Sp_buf, PED_MAXLEN, H_fpFam) ; N_sampOrig++)
		_loadRealTimeFamData(Sp_buf, N_sampOrig);
	LOG("%d samples detected before filtering\n", N_sampOrig);
	LOG("%d samples detected after filtering\n", N_sample);

	/* --nskip */
	wsUint N_skip = OPT_NUMBER(nskip);
	if (IS_ASSIGNED(nskip)) for (wsUint i=0 ; i<N_skip ; i++)
		if (C_file.gets(Sp_buf, PED_MAXLEN) == 0)
			halt_fmt(WISARD_NULL_PEDCONTENTS, S_fn);

	/* Get the number of variants */
	char *a, *b;
	for (N_vrtOrig=N_variant=0 ; C_file.gets(Sp_buf, PED_MAXLEN) ; N_vrtOrig++) {
		xVariant X_vrt	= {0, };
		char B_filt	= 0;

		/* S_buf == chromosome */
		getString(&Sp_buf, &a);
		X_vrt.chr = OPT_ENABLED(autoonly)?getChrAuto(Sp_buf):getChr(Sp_buf);

		/* a == variant name */
		getString(&a, &b);
		wsAlloc(X_vrt.name, char, strlen(a)+1);
		strcpy(X_vrt.name, a);

		/* b = genet. dist */
		if (!OPT_ENABLED(nogdist)) {
			getString(&b, &a);
			X_vrt.gdist		= (wsReal)atof(b);
		} else
			X_vrt.gdist		= W0;

		/* a == pos */
		if (!OPT_ENABLED(nopos)) {
			getString(&a, &b);

			/* Position validity check */
			for (char *p=a ; *p ; p++)
				if (*p < '0' || *p > '9')
					halt_fmt(WISARD_INVL_VARIANTPOS, X_vrt.name, a);
			X_vrt.pos		= atoi(a);
		} else
			X_vrt.pos		= 0;

		/* Check variant filtering conditions */
		B_filt |= _filterVariantPosition(X_vrt);

		/* Insert information */
		if (B_filt) {
			Bv_filtVrt.push_back(1);
			DEALLOC(X_vrt.name);
		} else {
			Xv_variant.push_back(X_vrt);
			Bv_filtVrt.push_back(0);
			N_variant++;
		}
	}
	C_file.rewind();
//	rewind(H_fp);

	/* --nskip */
	if (IS_ASSIGNED(3)) for (wsUint i=0 ; i<N_skip ; i++)
		if (C_file.gets(Sp_buf, PED_MAXLEN) == 0)
			halt_fmt(WISARD_NULL_PEDCONTENTS, S_fn);

	/* Substitute Xa_vrt if --updvariant */
	if (B_altVariantDef) {
		/* --updgdist, --updpos, --updchr, --updname? */
		wsUint N_sub = IS_ASSIGNED(updgdist) + IS_ASSIGNED(updpos) +
			IS_ASSIGNED(updchr) + IS_ASSIGNED(updname);

		if (Bv_filtVrt.size() != Xv_altVariantDef->size())
			halt("Alternative variant definition size[%d] does not match "
				"with retrieved variant count[%d]", Xv_altVariantDef->size(),
				Bv_filtVrt.size());

		LOG("Variant information is altered by [%s]\n", OPT_STRING(updvariant));
		if (N_sub > 0) for (vVariant_it s=Xv_altVariantDef->begin(),t=Xv_variant.begin() ;
			t != Xv_variant.end() ; t++,s++) {
			if (IS_ASSIGNED(updchr)) t->chr		= s->chr;
			if (IS_ASSIGNED(updgdist)) t->gdist	= s->gdist;
			if (IS_ASSIGNED(updpos)) t->pos		= s->pos;
			if (IS_ASSIGNED(updname)) {
				free(t->name);
				t->name		= s->name;
			}
		} else for (vVariant_it s=Xv_altVariantDef->begin(),t=Xv_variant.begin() ;
			t != Xv_variant.end() ; t++,s++) {
			t->chr		= s->chr;
			t->gdist	= s->gdist;
			t->pos		= s->pos;
			free(t->name);
			t->name		= s->name;
		}
	} else {
		/* Nothing can be * */
		if ((IS_ASSIGNED(updchr) && OPT_STRING(updchr)[0] == '*' && !OPT_STRING(updchr)[0]) ||
			(IS_ASSIGNED(updpos) && OPT_STRING(updpos)[0] == '*' && !OPT_STRING(updpos)[0]) ||
			(IS_ASSIGNED(updgdist) && OPT_STRING(updgdist)[0] == '*' && !OPT_STRING(updgdist)[0]) ||
			(IS_ASSIGNED(updname) && OPT_STRING(updname)[0] == '*' && !OPT_STRING(updname)[0]))
			halt("Without --updvariant, nothing of --updchr/--updpos/--updgdist/--updname can be default value[*]");

		/* FIXME : pair-style update */
	}

	/* Allocate required memory */
	_initMemories();

	/* Make complete family data */
	_registMissingFounder();

	/* Re-processing FAM */
	
	/* Read FAM file : main */
	rewind(H_fpFam);
	N_founder = 0;
	char *Sp_tmp1 = NULL;
	char *Sp_tmp2 = NULL;
	/* should be --1 but forget check */
	wsUint N_0 = 0;
	for (wsUint i=0,j=0 ; i<N_sampOrig ; i++) {
		if (fgets(Sp_buf, 1024, H_fpFam) == 0)
			halt_fmt(WISARD_INVL_EOF_LINE, "TFAM file", S_fnTFAM, N_sampOrig,
				i+1);
//			halt("Unexpected end of file found at line %d", i+1);
		string	S_currIID;

		/* Pass FID/IID/PAT/MAT/SEX */
		getString(&Sp_buf, &Sp_tmp1);	// Sp_buf  = FID
		getString(&Sp_tmp1, &Sp_tmp2);	// Sp_tmp1 = IID
		S_currIID = (string)Sp_tmp1;
		xSample	&X_sample = Xm_iid2data[S_currIID];
		getString(&Sp_tmp2, &Sp_tmp1);	// Sp_tmp2 = PAT
		getString(&Sp_tmp1, &Sp_tmp2);	// Sp_tmp1 = MAT
		getString(&Sp_tmp2, &Sp_tmp1);	// Sp_tmp2 = SEX

		/* If current samples should be filtered, do not retrieve anything */
		if (Bv_filtSamp[i] == 1) {
			X_sample.N_idx = SAMP_NODATA;
			continue;
		}

		/* Is this sample founder? */
		if ((Ba_isFounder[j] = !X_sample.Xp_mat && !X_sample.Xp_pat) != 0)
			N_founder++;
		/* Allocate initial data index */
		X_sample.N_idx = j;
		Xa_sampleV2.push_back(&X_sample);

		/* a = phenotype, retrieve anyway, because alternative phenotype should
		 * processed AFTER this procedure */
		getString(&Sp_tmp1, &Sp_tmp2);
		char *Sp_shouldNull = NULL;
		Ra_pheno[0][j] = (wsReal)str2dbl(Sp_tmp1, &Sp_shouldNull);

		/* Phenotype value is not control nor case */
		if (Ra_pheno[0][j] != (wsReal)OPT_NUMBER(phenoCtrl) &&
			Ra_pheno[0][j] != (wsReal)OPT_NUMBER(phenoCase)) {
			if (!strcmp(Sp_tmp1, OPT_STRING(mispheno))) {
				/* It is missing phenotype */
				Ra_pheno[0][j] = WISARD_NA_REAL;
//				LOG("Missing phenotype sample [%s:%s]\n",
//					X_sample.S_FID.c_str(), X_sample.S_IID.c_str());
			} else {
				if (Sp_shouldNull && Sp_shouldNull[0])
					halt_fmt(WISARD_INVL_FILE_INVALID_DESC, "TFAM file",
						S_fnTFAM, "phenotype", Sp_tmp1, i+1);
// 					halt("Phenotype [%s] for sample [%s] looks invalid",
// 					Sp_tmp1, X_sample.S_FID.c_str());

				/* Special for 0 */
				if (!strcmp(Sp_tmp1, "0")) N_0++;
				/* Otherwise, it is continuous phenotype */
				Ba_typePheno[0] = 1;
			}
		}
		/* Filtering option */
// 		else if ((OPT_ENABLED(filcase) && Ra_pheno[0][j] == OPT_NUMBER(phenoCase)) ||
// 			(OPT_ENABLED(filcontrol) && Ra_pheno[0][j] == OPT_NUMBER(phenoCtrl)))
// 			Ba_isSampFilt[j] = 1;

		/* Status */
		if (((++j)%100) == 0) notice("%d samples retrieved\r", j);
	}
	/* Possible chance to miss --1 */
	if (N_0) {
		LOGwarn("This file could be coded 0/1 but not assumed that. Please check it!\n");
		LOG("        ( 0 found %d times )\n", N_0);
	}

	//fclose(H_famFile);
	LOG("%d samples found in FAM file (%d filtered)\n", N_sample, N_sampOrig-N_sample);

	/* close stream */
	fclose(H_fpFam);

	/* Copy original filtering buffer */
	vBool Ba_isSampFiltOrig;
	Ba_isSampFiltOrig.resize(Bv_filtSamp.size());
	copy(Bv_filtSamp.begin(), Bv_filtSamp.end(),
		Ba_isSampFiltOrig.begin());
	vBool Ba_isVrtFiltOrig;
	Ba_isVrtFiltOrig.resize(Bv_filtVrt.size());
	copy(Bv_filtVrt.begin(), Bv_filtVrt.end(),
		Ba_isVrtFiltOrig.begin());

	/* Resize current buffer and clear */
	Bv_filtSamp.clear();
	Bv_filtVrt.clear();
	Bv_filtSamp.resize(N_sample, 0);
	Bv_filtVrt.resize(N_variant, 0);

	//Rp_pheno	= Ra_pheno[0];

	/* --varresize */
	if (IS_ASSIGNED(varresize)) {
		wsReal R_valSS = OPT_REAL(varresize);
		wsUint N_sel = R_valSS >= W1 ?
			(wsUint)R_valSS : (wsUint)round(R_valSS*N_vrtOrig);

		/* Check it should be resized */
		if (N_vrtOrig <= N_sel) {
			LOG("# of variants in dataset[%d] is larger than the number of variants to be resized[%d]"
				" so the dataset will not be resized\n", N_vrtOrig, N_sel);
		} else {
			/* Randomly select variants */
			for (wsUint N_chosen=0 ; N_chosen<N_sel ; ) {
				/* Select one */
				wsUint N_curr = rand()%N_vrtOrig;

				if (Ba_isVrtFiltOrig[N_curr] == 0) {
					Ba_isVrtFiltOrig[N_curr] = 1;
					N_chosen++;
				}
			}
		}		
	}

	char *Na_vrtData = NULL;
	wsAlloc(Na_vrtData, char, N_sample);

	/* Read TPED file */
	//vSampPtr &Xa_samp = getSample();
	wsUint	fsIdx=0, dsIdx=0;
	//xModel	N_gmodel = _getModel();
	for ( ; C_file.gets(Sp_buf, PED_MAXLEN) ; fsIdx++) {
		if (Ba_isVrtFiltOrig[fsIdx]) continue;

		/* S_buf == chromosome */
		getString(&Sp_buf, &a);
		/* a == variant name */
		getString(&a, &b);
		/* b = genet. dist */
		getString(&b, &a);
		/* a == pos */
		getString(&a, &b);
		a = b;

		/* Process variant */
		xVariant&	X_vrt	= Xv_variant[dsIdx];
		wsUint		N_gtypd	= 0;
		wsUint		N_gsum	= 0;
		for (wsUint fiIdx=0,diIdx=0 ; fiIdx<N_sampOrig ; fiIdx++,a=b) {
			if (Ba_isSampFiltOrig[fiIdx]) continue;
			if (diIdx>=N_sampOrig)
				halt("Variant [%s] have more than [%d] "
					"samples, than expected!", X_vrt.name, N_sampOrig);

			/* Retrieve genotype */
			int N_ret = getDiploid(X_vrt, a, &b);

			/* Process genotype */
			switch (N_ret) {
			case -1: halt_fmt(WISARD_INVL_PEDCONTENTS, fsIdx+1, S_fn);
			case -2: halt_fmt(WISARD_CANT_MULTILOCUS_WO_INDEL, fsIdx+1, S_fn);
			case -3: halt_fmt(WISARD_CANT_ONEALLELE_MISSING, X_vrt.name, fsIdx+1);
			/* Successfully interpreted */
			case 0: case 1: case 2:
				/* Model-dependent processing */
// 				switch (N_gmodel) {
// 				case MD_DOMINANT:
// 					N_ret = (N_ret+2)>>2; break;
// 				case MD_RECESSIVE:
// 					N_ret = (N_ret+1)>>1; break;
// 				case MD_MULTIPLICATIVE:
// 					N_ret = (N_ret+2)&0x05; break;
// 				default: break;
// 				}
				N_gtypd++;
				N_gsum += N_ret;
				Na_vrtData[diIdx++] = (char)N_ret;
				break;
			case -9:
				N_missGeno++;
				Xa_sampleV2[diIdx]->B_isMissing = -1;
				Na_vrtData[diIdx++] = (char)N_ret;
				break;
			/* Parsed value is invalid */
			default: halt_fmt(WISARD_SYST_INVL_PEDGENO, N_ret, X_vrt.name, fsIdx+1);
			}
		}
		/* Insert to Na_vrtData */
		if ((wsReal)N_gsum / (wsReal)(N_gtypd<<1) > REAL_CONST(0.5)) {
			/* Flip */
			for (wsUint i=0 ; i<N_sample ; i++)
				if (!isMissing(Na_vrtData[i])) Na_geno[i][dsIdx] = (~Na_vrtData[i])&0x03;
		} else {
			/* No flip */
			for (wsUint i=0 ; i<N_sample ; i++)
				Na_geno[i][dsIdx] = Na_vrtData[i];
		}

		/* Increase the number of variant */
		dsIdx++;
	}

	/* If it is PED file, extra space is required to calculate MA */
// 	MULTI_MALLOC(Na_pedData, USHORT_t*, N_sample);
// 	for (wsUint i=0 ; i<N_sample ; i++) {
// 		MULTI_MALLOC(Na_pedData[i], USHORT_t, N_variant);
// 	}

	/* Read */
// 	wsUint diIdx=0;
// 	char Ba_chrFound[128] = {0, };
// 	map<string, int> Xa_family;
// 
// 	/* Load PED file main */
// 	for (wsUint fiIdx=0 ; fgets(Sp_buf, PED_MAXLEN, H_fp) ; fiIdx++) {
// 		/* Do not process if this sample filtered */
// 		if (Ba_isSampFiltOrig[fiIdx]) continue;
// 
// 		/* Loop expected # of variants times */
// 		char *a=Sp_buf, *b;
// 		for (wsUint fsIdx=0 ; a ; fsIdx++) {
// 			char *Sp_1, *Sp_2;
// 
// 			/* Retrieve */
// 			getString(&a, &b);
// 			Sp_1 = a;
// 			if (b == NULL) halt("Unexpected end");
// 			getString(&b, &a);
// 			Sp_2 = b;
// 
// 			/* Check */
// 			if (OPT_ENABLED(indel)) {
// 				/* Single missing? */
// 				
// 			} else {
// 				/* More than 2 bytes? */
// 				if (Sp_1[1]) halt_fmt(WISARD_INVL_ALLELEFORMAT, fiIdx+1,
// 					fsIdx+1, Sp_1);
// 				if (Sp_2[1]) halt_fmt(WISARD_INVL_ALLELEFORMAT, fiIdx+1,
// 					fsIdx+1, Sp_2);
// 
// 			}
// 		}
// 		
// 		if (OPT_ENABLED(indel)) {
// 			if (_loadPED_INDELmain(Sp_buf, fiIdx, diIdx, B_isRawFile,
// 				Rp_pheno, Ba_isSampFiltOrig, Ba_isVrtFiltOrig, Ba_chrFound))
// 				continue;
// 		} else {
// 			if (_loadPED_SNVmain(Sp_buf, fiIdx, diIdx, B_isRawFile, Rp_pheno,
// 				Ba_isSampFiltOrig, Ba_isVrtFiltOrig, Ba_chrFound))
// 				continue;
// 		}

// 		/* Increase data index */
// 		diIdx++;
// 		if ((diIdx%10) == 0)
// 			notice("%d samples retrieved\r", diIdx);
// 	}
// 	DEALLOC(Sp_buf);
// 	LOG("%d samples retrieved\n", diIdx);
// 
// 	/* PED file requires minor-allele determination */
// 	if (!B_isRawFile) {
// 		_loadPED_asAdditive(Ba_chrFound);
// 
// 		wsUint i;
// 		cExporter *Cp_exp = NULL;
// 		for (i=0 ; i<N_variant ; i++)
// 			if (Xa_vrt[i].filter) {
// 				Cp_exp = new cExporter("alleleic.excl");
// 				break;
// 			}
// 			for ( ; i<N_variant ; i++)
// 				if (Xa_vrt[i].filter) {
// 					Cp_exp->fmt("%s\n", Xa_vrt[i].name);
// 					Ba_isVrtFilt[i] = 1;
// 					N_szFilteredVrt++;
// 				}
// 				if (N_szFilteredVrt)
// 					LOG("[%d] variants will be excluded because they have more than "
// 					"two alleles\n", N_szFilteredVrt);
// 
// 				/* In default, someone having NO genotype will be excluded */
// 				/* Filtering by --geno */
// 				for (wsUint i=0 ; i<N_sample ; i++) {
// 					if (Ba_isSampFilt[i]) continue;
// 
// 					if (Xa_sampleV2[i]->N_missGeno == (int)N_variant || (IS_ASSIGNED(filgind) &&
// 						isInRange(
// 						OPT_RANGE(filgind),
// 						(Xa_sampleV2[i]->N_missGeno/(wsReal)N_variant))
// 						) || (IS_ASSIGNED(incgind) &&
// 						!isInRange(
// 						OPT_RANGE(incgind),
// 						(Xa_sampleV2[i]->N_missGeno/(wsReal)N_variant))
// 						)) {
// 							Ba_isSampFilt[i] = 1;
// 							N_szFilteredSamp++;
// 					}
// 				}
// 	}
// 
// 	LOG("%d samples were retrieved\n", N_sample);
}	

void cPlinkIO::_loadLGEN(char *S_fn)
{
	char		*Sp_buf		= NULL;
	char		*Sp_tmp1	= NULL;
	cStrFile	C_file(S_fn, "LGEN file");
//	wsReal		*Rp_pheno;					///< Phenotype in PED/RAW file

	/* Allocate buffer for read PED file */
	wsAlloc(Sp_buf, char, 2048);

	/* Get variant count */
	char S_fnMap[512];
	otherExt(S_fn, "map", S_fnMap);
	_loadPED_countVariant(S_fn, 0, Sp_tmp1);

	/* Get sample count */
	_loadBED_FAM(S_fn);
	LOG("%d samples detected before filtering\n", N_sampOrig);
	LOG("%d samples detected after filtering\n", N_sample);

	/* Make complete family data */
	_registMissingFounder();

	/* Build map for variants */
	mDataIdx M_vrt;
	wsUint I=0;
	FOREACHDO (vVariant_it, Xv_variant, it, I++)
		M_vrt.insert(make_pair(it->name, I));

	/* Copy original filtering buffer */
	vBool Ba_isSampFiltOrig;
	Ba_isSampFiltOrig.resize(Bv_filtSamp.size());
	copy(Bv_filtSamp.begin(), Bv_filtSamp.end(),
		Ba_isSampFiltOrig.begin());
	vBool Ba_isSVrtFiltOrig;
	Ba_isSVrtFiltOrig.resize(Bv_filtVrt.size());
	copy(Bv_filtVrt.begin(), Bv_filtVrt.end(),
		Ba_isSVrtFiltOrig.begin());

	/* Resize current buffer and clear */
	Bv_filtSamp.clear();
	Bv_filtVrt.clear();
	Bv_filtSamp.resize(N_sample, 0);
	Bv_filtVrt.resize(N_variant, 0);

	/* Allocate required memory */
//	_initMemories();
//	Rp_pheno	= Ra_pheno[0];

	/* If it is PED file, extra space is required to calculate MA */
	wsAlloc(Na_charData, USHORT_t*, N_sample);
	for (wsUint i=0 ; i<N_sample ; i++) {
		wsAlloc(Na_charData[i], USHORT_t, N_variant);
	}

	/* Read */
	wsUint diIdx=0;
	char Ba_chrFound[128] = {0, };
	map<string, int> Xa_family;

	/* Load PED file main */
	char B_consec = OPT_ENABLED(consecallele);
	char S_delim = IS_ASSIGNED(sepallele) ? OPT_STRING(sepallele)[0] : 0xff;
	cExporter *Cp_misInds = NULL;
	for (wsUint fiIdx=0 ; C_file.gets(Sp_buf, PED_MAXLEN) ; fiIdx++) {
		/* Sp_buf = FID */
		char *a, *b;
		getString(&Sp_buf, &a);
		if (a == NULL) halt_fmt(WISARD_INVL_FILE_INCMPL, "LGEN file", S_fn, fiIdx+1);
		/* a = IID */
		getString(&a, &b);
		if (b == NULL) halt_fmt(WISARD_INVL_FILE_INCMPL, "LGEN file", S_fn, fiIdx+1);
		/* Check the sample exists */
		mSamp_it X_res = Xm_iid2data.find(a);
		if (X_res == Xm_iid2data.end())
			halt_fmt(WISARD_INVL_SAMPNAME, Sp_buf, a, fiIdx+1);
// 			halt("[%s::%s] have genotype data, but not exists in FAM file, at line [%d]\n",
// 				Sp_buf, a, fiIdx);
		else if (X_res->second.S_FID.compare(Sp_buf))
			halt("Sample [%s] registered at family [%s] but marked family [%s] in genotype data, at line [%d]\n",
				a, X_res->second.S_FID.c_str(), Sp_buf, fiIdx);

		/* b = vrtid */
		getString(&b, &a);
		if (a == NULL) halt_fmt(WISARD_INVL_FILE_INCMPL, "LGEN file", S_fn, fiIdx+1);
		/* Check variant ID */
		mDataIdx_it X_rVrt = M_vrt.find((string)b);
		if (X_rVrt == M_vrt.end()) continue; /* Pass if no exist */

		/* a = al1 */
		char *c = NULL;
		if (B_consec) {
			getString(&a, &c);
			b = a+1;
			if (a[2]) halt_fmt(WISARD_INVL_FILE_INVALID_DESC, "LGEN file",
				S_fn, "variant data", a, fiIdx+1);
//				halt("Variant at line [%d] is invalid", fiIdx);
		} else {
			getStringDelim(&a, &b, S_delim);
			if (b == NULL) halt_fmt(WISARD_INVL_FILE_INCMPL, "LGEN file", S_fn, fiIdx+1);
			/* b = al2 */
			getString(&b, &c);
		}

		if (OPT_ENABLED(indel))
			_indelSub(a, b, Xv_variant[X_rVrt->second], fiIdx,
				Na_charData[X_res->second.N_idx][X_rVrt->second],
				Ba_chrFound);
		else {
			char gen1, gen2;

			if (a[1] && !B_consec)
				halt("Variant at line [%d] is invalid", fiIdx);
//				halt_fmt(WISARD_INVL_ALLELEFORMAT, dsIdx<<1, diIdx, Sp_tmp1);
			gen1 = a[0];
			if (b[1])
				halt("Variant at line [%d] is invalid", fiIdx);
//				halt_fmt(WISARD_INVL_ALLELEFORMAT, dsIdx<<1, diIdx, Sp_tmp2);
			gen2 = b[0];

			/* Capitalize */
			if (gen1 >= 'a' && gen1 <= 'z') gen1 -= 0x20;
			if (gen2 >= 'a' && gen2 <= 'z') gen2 -= 0x20;

			/* Mark it exist in PED file */
			Ba_chrFound[(int)gen1] = 1;
			Ba_chrFound[(int)gen2] = 1;

			Na_charData[X_res->second.N_idx][X_rVrt->second] = (gen1<<8)|gen2;
		}

		/* Increase data index */
		diIdx++;
		if ((diIdx%10000) == 0)
			notice("%d genotypes retrieved\r", diIdx);
	}
	if (Cp_misInds)
		delete Cp_misInds;
	C_file.close();
	DEALLOC(Sp_buf);
	LOG("%d samples retrieved\n", diIdx);

	/* PED file requires minor-allele determination */
	_loadPED_asAdditive(Ba_chrFound);

	wsUint i;
/**/cExporter *Cp_exp = NULL;
	for (i=0 ; i<N_variant ; i++)
		if (Xv_variant[i].filter) {
/**/		Cp_exp = cExporter::summon("alleleic.excl");
			break;
		}
		for ( ; i<N_variant ; i++)
			if (Xv_variant[i].filter) {
				Cp_exp->fmt("%s\n", Xv_variant[i].name);
				Bv_filtVrt[i] = 1;
				N_szFilteredVrt++;
			}
			if (Cp_exp) delete Cp_exp;
			if (N_szFilteredVrt)
				LOG("[%d] variants will be excluded because they have more than "
				"two alleles\n", N_szFilteredVrt);

			/* In default, someone having NO genotype will be excluded */
			/* Filtering by --geno */
			for (wsUint i=0 ; i<N_sample ; i++) {
				if (Bv_filtSamp[i]) continue;

				wsReal R_gr = W1 - Xa_sampleV2[i]->N_missGeno / (wsReal)N_variant;
				if ((Xa_sampleV2[i]->N_missGeno == (int)N_variant && !OPT_ENABLED(nodata))
					|| (IS_ASSIGNED(filgind) &&
						isInRange(OPT_RANGE(filgind), R_gr)
					) || (IS_ASSIGNED(incgind) &&
						!isInRange(OPT_RANGE(incgind), R_gr)
					)) {
						Bv_filtSamp[i] = 1;
						N_szFilteredSamp++;
				}
			}

	LOG("%d samples were retrieved\n", N_sample);
}

char cPlinkIO::_loadPED_SNVmain(char *Sp_buf, wsUint fiIdx, wsUint diIdx,
	char B_isRawFile, wsReal *Rp_pheno, vBool& Ba_isSampFiltOrig,
	vBool& Ba_isVrtFiltOrig, char *Ba_chrFound, cExporter **Cp_misInds)
{
	char	*Sp_tmp1 = NULL;
	char	*Sp_tmp2 = NULL;
	string	S_currFID, S_currIID;
	//		int		ns = 0;
	int		N_missGenoCurrSamp = 0;

	/* If current samples should be filtered, do not retrieve anything */
	if (Ba_isSampFiltOrig[fiIdx] == 1)
		return 1;

	/* Get FID & IID */
	Sp_tmp2 = loadFIDIID(Sp_buf, S_currFID, S_currIID);

	xSample* Xp_sample = &(Xm_iid2data[S_currIID]);
	if (Xp_sample->N_idx != SAMP_NODATA) {
		if (!OPT_ENABLED(dupnaming))
			halt("ERROR : Duplicated IID [%s]", S_currIID.c_str());
		else {
			for (wsUint NN=1 ; ; NN++) {
				char S[32];
				sprintf(S, "%d", NN);
				string S_altIID = S_currIID + "_" + S;
				Xp_sample = &(Xm_iid2data[S_altIID]);
				if (Xp_sample->N_idx == SAMP_NODATA)
					break;
			}
		}
	}

	if (!OPT_ENABLED(noparent)) {
		getString(&Sp_tmp2, &Sp_tmp1);	// Sp_tmp2 = PAT
		getString(&Sp_tmp1, &Sp_tmp2);	// Sp_tmp1 = MAT
	}
	if (OPT_ENABLED(nosex)) Sp_tmp1 = Sp_tmp2;
	else getString(&Sp_tmp2, &Sp_tmp1);	// Sp_tmp2 = SEX

	/* Is this sample founder? */
	Ba_isFounder[diIdx] = !Xp_sample->Xp_mat && !Xp_sample->Xp_pat;
	if (Ba_isFounder[diIdx]) N_founder++;

	/* Allocate initial data index */
	Xp_sample->N_idx = diIdx;
	Xa_sampleV2.push_back(Xp_sample);

	/* Sp_tmp1 = phenotype */
	char *Sp_shouldNull = NULL;
	if (OPT_ENABLED(nopheno)) {
		Sp_tmp2 = Sp_tmp1;
		Ra_pheno[0][diIdx] = WISARD_NA_REAL;
	} else {
		getString(&Sp_tmp1, &Sp_tmp2);
		Ra_pheno[0][diIdx] = (wsReal)str2dbl(Sp_tmp1, &Sp_shouldNull);
	}
	//	Rp_pheno[diIdx] = (wsReal)atof(Sp_tmp1);

	/* If phenotype is missing */
	if (!strcmp(Sp_tmp1, OPT_STRING(mispheno))) {
		Rp_pheno[diIdx] = WISARD_NA_REAL;
		LOG("Missing phenotype [%s] sample [%s:%s]\n", Sp_tmp1,
			Xp_sample->S_FID.c_str(), Xp_sample->S_IID.c_str());
	} else if (Rp_pheno[diIdx] != OPT_NUMBER(phenoCtrl) &&
		Rp_pheno[diIdx] != OPT_NUMBER(phenoCase)) {
		if (Sp_shouldNull && Sp_shouldNull[0])
			halt("Phenotype [%s] for sample [%s] looks invalid",
				Sp_tmp1, Xp_sample->S_FID.c_str());

		Ba_typePheno[0] = 1; /* Set to continuous */
	}
// 	else if ((OPT_ENABLED(filcase) && Rp_pheno[diIdx] == OPT_NUMBER(phenoCase)) ||
// 		(OPT_ENABLED(filcontrol) && Rp_pheno[diIdx] == OPT_NUMBER(phenoCtrl)))
// 		Ba_isSampFilt[diIdx] = 1;

	Sp_tmp1 = Sp_tmp2;
	char S_delim = IS_ASSIGNED(sepallele) ? OPT_STRING(sepallele)[0] : 0xff;
	wsUint fsIdx = 0, dsIdx = 1;
	wsStrCst S_misGeno = IS_ASSIGNED(misgeno) ? OPT_STRING(misgeno) : NULL;
	for ( ; fsIdx<N_vrtOrig ; fsIdx++) {
		if (Sp_tmp1 == NULL)
			break;

		if (B_isRawFile) {
			getString(&Sp_tmp1, &Sp_tmp2);

			/* Retrieve information only if it is not filtered */
			if (Ba_isVrtFiltOrig[fsIdx] == 0) {
				char *Sp_ptr = NULL;
				char tmp = (char)strtoul(Sp_tmp1, &Sp_ptr, 10);//atoi(Sp_tmp1);
				if (Sp_ptr && Sp_ptr[0]) {
					/* Check is it NA */
					if (S_misGeno && !stricmp(Sp_ptr, S_misGeno)) {
						tmp = -1;
					} else halt("Input [%s] looks like RAW file, but non-numeric "
						"character [%s] found in line [%d]!",
						OPT_STRING(ped), Sp_tmp1, fiIdx);
				}

				/* Only 0/1/2 are acceptable, others are missing */
				if (tmp < 0 || tmp > 2 || Sp_ptr[0]) {
					Xa_sampleV2[diIdx]->B_isComplete = -1;
					Xv_variant[dsIdx-1].nmis++;
					tmp = WISARD_NA;
					N_missGeno++;
					N_missGenoCurrSamp++;
				}

				Na_geno[diIdx][dsIdx-1] = tmp;
				Sp_tmp1 = Sp_tmp2;
				dsIdx++;
			}
		} else {
			char gen1, gen2;

			/* Get biallele from buffer */
			if (OPT_ENABLED(consecallele)) {
				getString(&Sp_tmp1, &Sp_tmp2);
				if (Sp_tmp1[2])
					halt_fmt(WISARD_INVL_ALLELEFORMAT, dsIdx<<1, diIdx, Sp_tmp1);
				gen1 = Sp_tmp1[0];
				gen2 = Sp_tmp1[1];
				Sp_tmp1 = Sp_tmp2;
			} else {
				getStringDelim(&Sp_tmp1, &Sp_tmp2, S_delim);
				if (Sp_tmp1[1])
					halt_fmt(WISARD_INVL_ALLELEFORMAT, dsIdx<<1, diIdx, Sp_tmp1);
				gen1 = Sp_tmp1[0];
				getString(&Sp_tmp2, &Sp_tmp1);
				if (Sp_tmp2[1])
					halt_fmt(WISARD_INVL_ALLELEFORMAT, dsIdx<<1, diIdx, Sp_tmp2);
				gen2 = Sp_tmp2[0];
			}

			/* Retrieve information only if it is not filtered */
			if (Ba_isVrtFiltOrig[fsIdx]) continue;

			/* Capitalize */
			if (gen1 >= 'a' && gen1 <= 'z') gen1 -= 0x20;
			if (gen2 >= 'a' && gen2 <= 'z') gen2 -= 0x20;

			/* Mark it exist in PED file */
			Ba_chrFound[(int)gen1] = 1;
			Ba_chrFound[(int)gen2] = 1;

			Na_charData[diIdx][dsIdx-1] = (gen1<<8)|gen2;
			dsIdx++;
		}
	}
	/* Check sanity of PED file */
	if ((dsIdx-1) != N_variant || fsIdx != N_vrtOrig)
		halt_fmt(WISARD_INVL_NVRTOFSAMPLE, Xp_sample->S_FID.c_str(),
			Xp_sample->S_IID.c_str(), N_vrtOrig, N_variant, fsIdx, dsIdx-1);

	if (N_vrtOrig && B_isRawFile) {
		/* In default, someone having NO genotype will be excluded */
		/* Filtering by --geno when RAW file */
		wsReal R_gr = W1 - N_missGenoCurrSamp / (wsReal)N_variant;
		if ((N_missGenoCurrSamp == (int)N_variant && !OPT_ENABLED(nodata)) ||
				(IS_ASSIGNED(filgind) &&
					isInRange(OPT_RANGE(filgind), R_gr)
				) ||
				(IS_ASSIGNED(incgind) &&
					!isInRange(OPT_RANGE(incgind), R_gr)
				)) {
			if (N_missGenoCurrSamp == (int)N_variant && Cp_misInds) {
				if (*Cp_misInds == NULL) {
					*Cp_misInds = cExporter::summon("allmiss.sample.lst");
					LOGoutput("Sample with all missing genotypes found, it "
						"will exported to [%s.allmiss.sample.lst]\n",
						OPT_STRING(out));
				}
				(*Cp_misInds)->fmt("%s	%s\n", Sp_buf, S_currIID.c_str());
			}
			Bv_filtSamp[diIdx] = 1;
			N_szFilteredSamp++;
		}
	}
	
	return 0;
}

char cPlinkIO::_loadPED_INDELmain(char *Sp_buf, wsUint fiIdx, wsUint diIdx,
	char B_isRawFile, wsReal *Rp_pheno, vBool& Ba_isSampFiltOrig,
	vBool& Ba_isVrtFiltOrig, char *Ba_chrFound, cExporter **Cp_misInds)
{
	char	*Sp_tmp1 = NULL;
	char	*Sp_tmp2 = NULL;
	string	S_currIID;
	//		int		ns = 0;
	int		N_missGenoCurrSamp = 0;

	/* If current samples should be filtered, do not retrieve anything */
	if (Ba_isSampFiltOrig[fiIdx] == 1)
		return 1;

	/* Get FID & IID */
	Sp_tmp2 = loadFIDIID(Sp_buf, S_currIID, S_currIID);

	xSample& X_sample = Xm_iid2data[S_currIID];
	if (X_sample.N_idx != SAMP_NODATA) {
		if (!OPT_ENABLED(dupnaming))
			halt("ERROR : Duplicated IID [%s]", S_currIID.c_str());
		else {
			for (wsUint NN=1 ; ; NN++) {
				char S[32];
				sprintf(S, "%d", NN);
				string S_altIID = S_currIID + "_" + S;
				X_sample = Xm_iid2data[S_currIID];
				if (X_sample.N_idx == SAMP_NODATA)
					break;
			}
		}
	}

	if (!OPT_ENABLED(noparent)) {
		getString(&Sp_tmp2, &Sp_tmp1);	// Sp_tmp2 = PAT
		getString(&Sp_tmp1, &Sp_tmp2);	// Sp_tmp1 = MAT
	}
	if (OPT_ENABLED(nosex)) Sp_tmp1 = Sp_tmp2;
	else getString(&Sp_tmp2, &Sp_tmp1);	// Sp_tmp2 = SEX

	/* Is this sample founder? */
	Ba_isFounder[diIdx] = !X_sample.Xp_mat && !X_sample.Xp_pat;
	if (Ba_isFounder[diIdx]) N_founder++;

	/* Allocate initial data index */
	X_sample.N_idx = diIdx;
	Xa_sampleV2.push_back(&X_sample);

	/* Sp_tmp1 = phenotype */
	char *Sp_shouldNull = NULL;
	if (OPT_ENABLED(nopheno)) {
		Ra_pheno[0][diIdx] = WISARD_NA_REAL;
		Sp_tmp2 = Sp_tmp1;
	} else {
		getString(&Sp_tmp1, &Sp_tmp2);
		Ra_pheno[0][diIdx] = (wsReal)str2dbl(Sp_tmp1, &Sp_shouldNull);
	}
//	Rp_pheno[diIdx] = (wsReal)atof(Sp_tmp1);

	/* If phenotype is missing */
	if (!strcmp(Sp_tmp1, OPT_STRING(mispheno))) {
		Rp_pheno[diIdx] = WISARD_NA_REAL;
		LOG("Missing phenotype [%s] sample [%s:%s]\n", Sp_tmp1,
			X_sample.S_FID.c_str(), X_sample.S_IID.c_str());
	} else if (Rp_pheno[diIdx] != OPT_NUMBER(phenoCtrl) &&
		Rp_pheno[diIdx] != OPT_NUMBER(phenoCase)) {
		if (Sp_shouldNull && Sp_shouldNull[0])
			halt("Phenotype [%s] for sample [%s] looks invalid",
			Sp_tmp1, X_sample.S_FID.c_str());

		Ba_typePheno[0] = 1; /* Set to continuous */
	}
// 	else if ((OPT_ENABLED(filcase) && Rp_pheno[diIdx] == OPT_NUMBER(phenoCase)) ||
// 		(OPT_ENABLED(filcontrol) && Rp_pheno[diIdx] == OPT_NUMBER(phenoCtrl)))
// 		Ba_isSampFilt[diIdx] = 1;

	char S_delim = IS_ASSIGNED(sepallele) ? OPT_STRING(sepallele)[0] : 0xff;
	Sp_tmp1 = Sp_tmp2;
	wsUint fsIdx = 0, dsIdx = 0;
	for ( ; fsIdx<N_vrtOrig ; fsIdx++) {
		//xVariant &X_vrt = getVariants()[fsIdx];
		if (Sp_tmp1 == NULL)
			break;
		if (Ba_isVrtFiltOrig[fsIdx] == 1)
			continue;

		/* Buffer, fileSnpI, dataSnpI, fileIndI, dataIndI */
		//proc2allele(Sp_tmp1, fsIdx, dsIdx, fiIdx, diIdx);

		/* NOTE : --indel won't care --consecallele */

		/* Get biallele from buffer */
		getStringDelim(&Sp_tmp1, &Sp_tmp2, S_delim);
		char *Sp_gen1 = Sp_tmp1;
		getString(&Sp_tmp2, &Sp_tmp1);
		char *Sp_gen2 = Sp_tmp2;

		/* Retrieve information only if it is not filtered */
		if (Bv_filtVrt[fsIdx]) continue;
		_indelSub(Sp_gen1, Sp_gen2, Xv_variant[dsIdx], fiIdx,
			Na_charData[diIdx][dsIdx], Ba_chrFound);
		dsIdx++;
	}
	/* Check sanity of PED file */
	if (dsIdx != N_variant || fsIdx != N_vrtOrig)
		halt_fmt(WISARD_INVL_NVRTOFSAMPLE, X_sample.S_FID.c_str(),
			X_sample.S_IID.c_str(), N_vrtOrig, N_variant, fsIdx, dsIdx);

	if (N_vrtOrig && B_isRawFile) {
		/* In default, someone having NO genotype will be excluded */
		/* Filtering by --geno when RAW file */
		wsReal R_gr = W1 - N_missGenoCurrSamp / (wsReal)N_variant;
		if ((N_missGenoCurrSamp == (int)N_variant && !OPT_ENABLED(nodata)) ||
			(IS_ASSIGNED(filgind) &&
				isInRange(OPT_RANGE(filgind), R_gr)
			) ||
			(IS_ASSIGNED(incgind) &&
				!isInRange(OPT_RANGE(incgind), R_gr)
			)) {
				if (N_missGenoCurrSamp == (int)N_variant) {
					if (Cp_misInds == NULL) {
						*Cp_misInds = cExporter::summon("allmiss.sample.lst");
						LOGoutput("Sample with no avail genotype found, it "
							"is exported to [%s.allmiss.sample.lst]\n", OPT_STRING(out));
					}
					(*Cp_misInds)->fmt("%s	%s\n", Sp_buf, S_currIID.c_str());
				}
				Bv_filtSamp[diIdx] = 1;
				N_szFilteredSamp++;
		}
	}

	return 0;
}

int calcMAF(int tid, void *shareData, void *data, void *res)
{
	int			i = *((int *)data);
	cIO			*Cp_IO = (cIO *)shareData;
	wsFmat		Rp_data = Cp_IO->getData();
	USHORT_t	**Np_pedData = Cp_IO->getCharData();
	int			Na_allele[128] = { 0, };
	char		Ca_allele[2] = { 127, 127 }, N_allele = 0;
	int			*Np_missGeno = (int *)res;
	wsUint		N_sample = Cp_IO->sizeSample();
	char		C_missing = ((cPlinkIO *)Cp_IO)->getPedMissingChar();
	xVariant&	X_vrt = Cp_IO->getVariant()[i];

	for (wsUint j = 0; j < N_sample; j++) {
		char C_al1 = Np_pedData[j][i] >> 8;
		char C_al2 = Np_pedData[j][i] & 0xff;

		/* Is this missing? */
		if (C_al1 == C_missing || C_al2 == C_missing) {
			Rp_data[j][i] = (wsFloat)WISARD_NA_REAL;
			Np_missGeno[tid]++;
		}
		else {
			/* Allele marking */
			if (Na_allele[(int)C_al1] == 0) {
				/* Non-biallele checking */
				if (N_allele >= 2) {
					LOG("Variant [%s] have more than two [%d] alleles! [ ", X_vrt.name, N_allele + 1);
					for (int i = 0; i < 128; i++) {
						if (Na_allele[i] > 0)
							LOG("%c ", i);
					}
					LOG("%c ]\n", C_al1);
					exit(1);
				}

				Ca_allele[(int)N_allele++] = C_al1;
			}
			Na_allele[(int)C_al1]++;

			if (Na_allele[(int)C_al2] == 0) {
				/* Non-biallele checking */
				if (N_allele >= 2)
					halt("Variant [%d] have more than two [%d] alleles! [ %c(%d) %c(%d) %c(%d) ]",
					i + 1, N_allele, Ca_allele[0], Ca_allele[0],
					Ca_allele[1], Ca_allele[1], C_al2, C_al2);

				Ca_allele[(int)N_allele++] = C_al2;
			}
			Na_allele[(int)C_al2]++;
		}
	}

	/* Determine minor and major */
	char C_minor;//, C_major;
	//	int N_minor, N_major;
	if (Na_allele[(int)Ca_allele[0]] > Na_allele[(int)Ca_allele[1]]) {
		//		C_major = Ca_allele[0];
		//		N_major = Na_allele[(int)Ca_allele[0]];
		if (Ca_allele[1] == 127)
			C_minor = 'N';
		else
			C_minor = Ca_allele[1];
		//		N_minor = Na_allele[(int)Ca_allele[1]];
	}
	else {
		//		C_major = Ca_allele[1];
		//		N_major = Na_allele[(int)Ca_allele[1]];
		if (Ca_allele[0] == 127)
			C_minor = 'N';
		else
			C_minor = Ca_allele[0];
		//		N_minor = Na_allele[(int)Ca_allele[0]];
	}
	//			fprintf(fp2, "SNP%d	%c(%d)	%c(%d)	%.2f%%\n", i+1, C_major, N_major, C_minor, N_minor, 100.0f*(N_minor/(wsReal)(N_major+N_minor)));
	if (C_minor == Ca_allele[0]) {
		X_vrt.al1 = Ca_allele[1];
		X_vrt.al2 = Ca_allele[0];
	}
	else {
		X_vrt.al1 = Ca_allele[0];
		X_vrt.al2 = Ca_allele[1];
	}

	for (wsUint j = 0; j < N_sample; j++) {
		char C_al1 = Np_pedData[j][i] >> 8;
		char C_al2 = Np_pedData[j][i] & 0xff;

		int val = 0;
		/* Additive coding */
		if (C_al1 != C_missing && C_al2 != C_missing) {
			if ((Ca_allele[0] == C_al1 && C_minor == Ca_allele[0]) ||
				(Ca_allele[1] == C_al1 && C_minor == Ca_allele[1]))
				val++;
			if ((Ca_allele[0] == C_al2 && C_minor == Ca_allele[0]) ||
				(Ca_allele[1] == C_al2 && C_minor == Ca_allele[1]))
				val++;

			Rp_data[j][i] = (wsFloat)val;
		}
	}

	return 0;
}

/* Calculate both of MAF and genotyping rate PER sample */
int calcMAFv2(int tid, void *shareData, void *data, void *res)
{
	cIO			*Cp_IO = (cIO *)shareData;
	char		**Np_geno = Cp_IO->getGenotype();
	USHORT_t	**Np_pedData = Cp_IO->getCharData();
	int			*Np_missGeno = (int *)res;
	wsUint		N_sample = Cp_IO->sizeSample();
	char		C_missing = ((cPlinkIO *)Cp_IO)->getPedMissingChar();
	//	vBool&		Xa_isSampFilt	= Cp_IO->getIsSampleFiltered();
	//	mSamp&		Xm_iid2data		= Cp_IO->getSampleData();
	int*		Np_data = (int *)data;
	wsUint		N_s = (wsUint)Np_data[0];
	wsUint		N_e = (wsUint)Np_data[1];
	vSampPtr&	Xa_sampleV2 = Cp_IO->getSample();

	if (Xa_sampleV2.size() < N_sample)
		halt("SYSERR: Sample size should be %d, but insufficient",
		N_sample);

	for (wsUint i = N_s; i < N_e; i++) {
		int			Na_allele[128] = { 0, };
		char		Ca_allele[2] = { 127, 127 }, N_allele = 0;
		xVariant&	X_vrt = Cp_IO->getVariant()[i];

		/* Do not process if it has filtered */
		if (X_vrt.filter == 1) continue;

		for (wsUint j = 0; j < N_sample; j++) {
			xSample	*Xp_sample = Xa_sampleV2[j];
			char C_al1 = Np_pedData[j][i] >> 8;
			char C_al2 = Np_pedData[j][i] & 0xff;

			if (C_al2 == -1)
				halt("Invalid allele found at variant [%s] of sample [%s]",
				X_vrt.name, Xp_sample->S_IID.c_str());

			/* Is this missing? */
			if (C_al1 == C_missing || C_al2 == C_missing) {
				Xp_sample->B_isComplete = -1;
				X_vrt.nmis++;
				Np_geno[j][i] = WISARD_NA;
				Np_missGeno[tid]++;
				Xp_sample->N_missGeno++;
			}
			else {
				/* Allele marking */
				if (Na_allele[(int)C_al1] == 0) {
					/* Non-biallele checking */
					if (N_allele >= 2) {
						LOG("Variant [%s] have more than two [%d] alleles! [ ", X_vrt.name, N_allele + 1);
						for (int i = 0; i < 128; i++) {
							if (Na_allele[i] > 0)
								LOGnf("%c ", i);
						}
						LOGnf("%c ]\n", C_al1);
						exit(1);
					}

					Ca_allele[(int)N_allele++] = C_al1;
				}
				Na_allele[(int)C_al1]++;

				if (Na_allele[(int)C_al2] == 0) {
					/* Non-biallele checking */
					if (N_allele >= 2)
						halt("Variant [%s] have more than two [%d] alleles at variant "
						"of sample [%s]! [ %c(%d) %c(%d) %c(%d) ]",
						X_vrt.name, N_allele, X_vrt.name, Xp_sample->S_IID.c_str(),
						Ca_allele[0], Ca_allele[0], Ca_allele[1],
						Ca_allele[1], C_al2, C_al2);

					Ca_allele[(int)N_allele++] = C_al2;
				}
				Na_allele[(int)C_al2]++;
			}
		}

		/* Determine minor and major */
		char C_minor;//, C_major;
		//	int N_minor, N_major;
		if (Na_allele[(int)Ca_allele[0]] > Na_allele[(int)Ca_allele[1]]) {
			//		C_major = Ca_allele[0];
			//		N_major = Na_allele[(int)Ca_allele[0]];
			if (Ca_allele[1] == 127)
				C_minor = 'N';
			else
				C_minor = Ca_allele[1];
			//		N_minor = Na_allele[(int)Ca_allele[1]];
		}
		else {
			//		C_major = Ca_allele[1];
			//		N_major = Na_allele[(int)Ca_allele[1]];
			if (Ca_allele[0] == 127)
				C_minor = 'N';
			else
				C_minor = Ca_allele[0];
			//		N_minor = Na_allele[(int)Ca_allele[0]];
		}
		//			fprintf(fp2, "SNP%d	%c(%d)	%c(%d)	%.2f%%\n", i+1, C_major, N_major, C_minor, N_minor, 100.0f*(N_minor/(wsReal)(N_major+N_minor)));
		if (Ca_allele[1] == 127) Ca_allele[1] = 0;

		/* Determine minor & major */
		if (C_minor == Ca_allele[0]) {
			X_vrt.al1 = Ca_allele[1];
			X_vrt.al2 = Ca_allele[0];
		}
		else {
			X_vrt.al1 = Ca_allele[0];
			X_vrt.al2 = Ca_allele[1];
		}

		/* If --acgt */
		Cp_IO->code2acgt(X_vrt, 0);
		// 	char *Na_ACGT = Cp_IO->getACGT();
		// 	if (Na_ACGT) {
		// 		if (!Na_ACGT[(wsUint)X_vrt.al1]) halt("Character [%c] in variant [%s] is not in --acgt",
		// 			X_vrt.al1, X_vrt.name);
		// 		if (X_vrt.al2 && !Na_ACGT[(wsUint)X_vrt.al2]) halt("Character [%c] in variant [%s] is not in --acgt",
		// 			X_vrt.al2, X_vrt.name);
		// 		X_vrt.al1 = Na_ACGT[(wsUint)X_vrt.al1];
		// 		if (X_vrt.al2) X_vrt.al2 = Na_ACGT[(wsUint)X_vrt.al2];
		// 	}

		//	xModel N_gmodel = Cp_IO->_getModel();
		for (wsUint j = 0; j < N_sample; j++) {
			char C_al1 = Np_pedData[j][i] >> 8;
			char C_al2 = Np_pedData[j][i] & 0xff;

			int N_ret = 0;
			/* Additive coding */
			if (C_al1 != C_missing && C_al2 != C_missing) {
				if ((Ca_allele[0] == C_al1 && C_minor == Ca_allele[0]) ||
					(Ca_allele[1] == C_al1 && C_minor == Ca_allele[1]))
					N_ret++;
				if ((Ca_allele[0] == C_al2 && C_minor == Ca_allele[0]) ||
					(Ca_allele[1] == C_al2 && C_minor == Ca_allele[1]))
					N_ret++;

				/* Model-dependent processing */
				//			switch (N_gmodel) {
				//			case MD_DOMINANT:
				//				N_ret = (N_ret+2)>>2; break;
				//			case MD_RECESSIVE:
				//				N_ret = (N_ret+1)>>1; break;
				//			case MD_MULTIPLICATIVE:
				//				N_ret = (N_ret+2)&0x05; break;
				//			default: break;
				//			}

				Np_geno[j][i] = (char)N_ret;
			}
		}
	}
	return 0;
}

void exportTPED(cIO* Cp_IO)
{
	if (IS_ASSIGNED(genoprob)) {
		LOGwarn("Conversion of genotype probability format does not supported, skips --maketped");
		return;
	}
	/* --outmisgeno : TPED '0' */
	const char*	S_misgeno		= IS_ASSIGNED(outmisgeno) ?
		OPT_STRING(outmisgeno) : "0";

	/* Export TFAM data */
	vSampPtr&	Xa_samp			= Cp_IO->getSample();
	const char*	Ba_chrExists	= Cp_IO->isChrExists();
	//	REAL_c		*Ra_pheno = getPheno();

	/* Cannot apply --outphenoonly */
	if (OPT_ENABLED(outphenoonly))
		LOGwarn("BED file does not support --outphenoonly, will be ignored\n");

	/* Allocate EXPORTERS */
	cExporter**	Cp_tped = NULL;
	cExporter**	Cp_tfam = NULL;
	wsCalloc(Cp_tped, cExporter*, MAX_NCHR + 1);
	wsCalloc(Cp_tfam, cExporter*, MAX_NCHR + 1);

	if (OPT_ENABLED(split)) {
		/* Make instance of EXPOTER ONLY IF exists */
		char S_ext[64];
		for (wsUint i = 0; i < NCHR_SPECIES; i++) if (Ba_chrExists[i]) {
			wsStrCst S_chr = getChrName2(i + 1);
			sprintf(S_ext, "chr%s.tped", S_chr);
			Cp_tped[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "chr%s.tfam", S_chr);
			Cp_tfam[i] = cExporter::summon(S_ext);
		}

		/* For Undetermined IF exists */
		if (Ba_chrExists[NCHR_SPECIES]) {
			Cp_tped[NCHR_SPECIES] = cExporter::summon("chrUn.tped");
			Cp_tfam[NCHR_SPECIES] = cExporter::summon("chrUn.tfam");
		}

		LOGoutput("Final dataset is exported with [%s.chr***.tped] [%s.chr***.tfam]\n",
			OPT_STRING(out), OPT_STRING(out));
	}
	else {
		/* Otherwise, use just one */
		Cp_tped[0] = cExporter::summon("tped"); /* CHECKED */
		Cp_tfam[0] = cExporter::summon("tfam"); /* CHECKED */
		LOGoutput("Final dataset is exported with [%s.tped] [%s.tfam]\n",
			OPT_STRING(out), OPT_STRING(out));
	}

	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

	wsUint I = 0;
	FOREACHDO(vSampPtr_it, Xa_samp, i, I++) {
		/* Determining phenotype */
		char	S_pheno[128];
		Cp_IO->getPhenoCode(I, S_pheno, Sp_case, Sp_ctrl);

		/* For phenotype, export header */
		for (wsUint j = 0; j <= NCHR_SPECIES; j++) if (Cp_tped[j])
			Cp_tfam[j]->fmt("%s %s %s %s %d %s\n", (*i)->S_FID.c_str(),
			(*i)->S_IID.c_str(),
			(*i)->Xp_pat ? (*i)->Xp_pat->S_IID.c_str() : "0",
			(*i)->Xp_mat ? (*i)->Xp_mat->S_IID.c_str() : "0",
			(*i)->N_sex, S_pheno);
	}

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	// 	LOGoutput("Final dataset is exported with [%s.tped] [%s.tfam]\n",
	// 		OPT_STRING(out), OPT_STRING(out));

	/* Export TPED data */
	char**		Na_data	= Cp_IO->getGenotype();
	vVariant&	Xa_snvs	= Cp_IO->getVariant();
	wsFloat**	Ra_data	= Cp_IO->getData(1);
	I = 0;
	FOREACHDO(vVariant_it, Xa_snvs, i, I++) {
		/* Determine target */
		wsUint N_target = 0;
		/* Target could be not 0 ONLY IF --split */
		if (OPT_ENABLED(split)) {
			if (i->chr <= 0 || (wsUint)i->chr > NCHR_SPECIES)
				N_target = NCHR_SPECIES;
			else N_target = i->chr - 1;
		}
		cExporter *Cp_currPED = Cp_tped[N_target];

		/* Print variant information */
		Cp_currPED->fmt("%s	%s	0	%d", getChrName2(i->chr), i->name, i->pos);

		/* Print variant data */
		if (IS_ASSIGNED(dosage)) FOREACH(vSampPtr_it, Xa_samp, j) {
			wsReal R_dsg = Ra_data[(*j)->N_idx][I];
			if (isMissingReal(R_dsg))
				Cp_currPED->fmt("  %s %s", S_misgeno, S_misgeno);
			else if (R_dsg < W0 || R_dsg > W2) {
				halt("Dosage [%s] of sample [%s] have invalid value [%g]",
					i->name, (*j)->S_IID.c_str(), R_dsg);
			}
			else if (!OPT_ENABLED(indel)) {
				if (R_dsg < REAL_CONST(0.5)) /* Major homo */
					Cp_currPED->fmt("  %s %s", i->indel1, i->indel1);
				else if (R_dsg < REAL_CONST(1.5)) /* Hetero */
					Cp_currPED->fmt("  %s %s", i->indel2, i->indel1);
				else /* Minor homo */
					Cp_currPED->fmt("  %s %s", i->indel2, i->indel2);
			}
			else {
				if (R_dsg < REAL_CONST(0.5)) /* Major homo */
					Cp_currPED->fmt("  %c %c", i->al1, i->al1);
				else if (R_dsg < REAL_CONST(1.5)) /* Hetero */
					Cp_currPED->fmt("  %c %c", i->al2, i->al1);
				else /* Minor homo */
					Cp_currPED->fmt("  %c %c", i->al2, i->al2);
			}
		} else FOREACH(vSampPtr_it, Xa_samp, j) {
			if (!OPT_ENABLED(indel)) switch (Na_data[(*j)->N_idx][I]) {
			case 0: /* Major homo */
				Cp_currPED->fmt("  %c %c", i->al1, i->al1); break;
			case 1: /* Hetero */
				Cp_currPED->fmt("  %c %c", i->al1, i->al2); break;
			case 2: /* Minor homo */
				Cp_currPED->fmt("  %c %c", i->al2, i->al2); break;
			case WISARD_NA:
				Cp_currPED->fmt("  %s %s", S_misgeno, S_misgeno); break;
			default:
				halt_fmt(WISARD_SYST_INVL_SNPDATA, i->name, (*j)->S_FID.c_str(),
					(*j)->S_IID.c_str(), Na_data[(*j)->N_idx][I]);
			}
			else switch (Na_data[(*j)->N_idx][I]) {
			case 0: /* Major homo */
				Cp_currPED->fmt("  %s %s", i->indel1, i->indel1); break;
			case 1: /* Hetero */
				Cp_currPED->fmt("  %s %s", i->indel1, i->indel2); break;
			case 2: /* Minor homo */
				Cp_currPED->fmt("  %s %s", i->indel1, i->indel2); break;
			case WISARD_NA:
				Cp_currPED->fmt("  %s %s", S_misgeno, S_misgeno); break;
			default:
				halt_fmt(WISARD_SYST_INVL_SNPDATA, i->name, (*j)->S_FID.c_str(),
					(*j)->S_IID.c_str(), Na_data[(*j)->N_idx][I]);
			}
		}

		Cp_currPED->put("\n");
	}
	for (wsUint i = 0; i <= NCHR_SPECIES; i++) if (Cp_tped[i]) {
		delete Cp_tped[i];
		delete Cp_tfam[i];
	}
	DEALLOC(Cp_tped);
	DEALLOC(Cp_tfam);
}

void exportPED(cIO* Cp_IO)
{
	if (IS_ASSIGNED(genoprob)) {
		LOGwarn("Conversion of genotype probability format does not supported, skips --makeped");
		return;
	}
	vVariant&	Xv_vrt			= Cp_IO->getVariant();
	vSampPtr&	Xv_samp			= Cp_IO->getSample();
	char**		Na_data			= Cp_IO->getGenotype();
	wsFloat**	Ra_data			= Cp_IO->getData(1);
	wsUint		N_variant		= Cp_IO->sizeVariant();
	const char*	Ba_chrExists	= Cp_IO->isChrExists();
	/* Allocate EXPORTERS */
	cExporter**	Cp_ped			= NULL;
	cExporter**	Cp_map			= NULL;
	/* --outmisgeno : PED '0' */
	const char*	S_misgeno		= IS_ASSIGNED(outmisgeno) ?
		OPT_STRING(outmisgeno) : (char *)"0";
	wsUint		N_exporter		= 0;
	/* --famsplit */
	mDataIdx	Xm_fid2idx;

	/* Cannot apply --outphenoonly */
	if (OPT_ENABLED(outphenoonly))
		LOGwarn("BED file does not support --outphenoonly, will be ignored\n");

	if (OPT_ENABLED(split)) {
		N_exporter = MAX_NCHR + 1;
		wsCalloc(Cp_ped, cExporter*, MAX_NCHR + 1);
		wsCalloc(Cp_map, cExporter*, MAX_NCHR + 1);

		/* Make instance of EXPOTER ONLY IF exists */
		char S_ext[64];
		for (wsUint i = 0; i < NCHR_SPECIES; i++) if (Ba_chrExists[i]) {
			wsStrCst S_chr = getChrName2(i + 1);
			sprintf(S_ext, "chr%s.ped", S_chr);
			Cp_ped[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "chr%s.map", S_chr);
			Cp_map[i] = cExporter::summon(S_ext);
		}

		/* For Undetermined IF exists */
		if (Ba_chrExists[NCHR_SPECIES]) {
			Cp_ped[NCHR_SPECIES] = cExporter::summon("chrUn.ped");
			Cp_map[NCHR_SPECIES] = cExporter::summon("chrUn.map");
		}

		LOGoutput("Final dataset is exported with [%s.chr***.ped] [%s.chr***.map]\n",
			OPT_STRING(out), OPT_STRING(out));
	} else if (OPT_ENABLED(famsplit)) {
		mFam& Xm_fam = Cp_IO->getFamilyData();
		N_exporter = (wsUint)Xm_fam.size();
		wsCalloc(Cp_ped, cExporter*, N_exporter);
		wsCalloc(Cp_map, cExporter*, N_exporter);

		/* For each fid, create index mapping */
		wsUint i = 0;
		FOREACH (mFam_it, Xm_fam, X_fam) {
			string S_fid = X_fam->first;
			Xm_fid2idx[S_fid] = i++;
		}

		/* Make instances */
		char S_ext[64];
		i = 0;
		FOREACH (mFam_it, Xm_fam, X_fam) {
			string S_fid = X_fam->first;
			sprintf(S_ext, "fam.%s.ped", S_fid.c_str());
			Cp_ped[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "fam.%s.map", S_fid.c_str());
			Cp_map[i] = cExporter::summon(S_ext);

			i++;
		}

		LOGoutput("Final dataset is exported with [%s.fam.***.ped] [%s.fam.***.map]\n",
			OPT_STRING(out), OPT_STRING(out));
	} else {
		N_exporter = 1;
		wsCalloc(Cp_ped, cExporter*, MAX_NCHR + 1);
		wsCalloc(Cp_map, cExporter*, MAX_NCHR + 1);

		/* Otherwise, use just one */
		Cp_ped[0] = cExporter::summon("ped"); /* CHECKED */
		Cp_map[0] = cExporter::summon("map"); /* CHECKED */
		LOGoutput("Final dataset is exported with [%s.ped] [%s.map]\n",
			OPT_STRING(out), OPT_STRING(out));
	}

	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

	for (vSampPtr_it i = Xv_samp.begin(); i != Xv_samp.end();) {
		xSample *Xp_currSamp = *i;

		/* Determining phenotype */
		char	S_pheno[128];
		Cp_IO->getPhenoCode(Xp_currSamp->N_idx, S_pheno, Sp_case, Sp_ctrl);

		/* For phenotype, export header */
		if (OPT_ENABLED(famsplit))
			Cp_ped[Xm_fid2idx[Xp_currSamp->S_FID]]->fmt("%s %s %s %s %d %s", Xp_currSamp->S_FID.c_str(),
				Xp_currSamp->S_IID.c_str(),
				Xp_currSamp->Xp_pat ? Xp_currSamp->Xp_pat->S_IID.c_str() : "0",
				Xp_currSamp->Xp_mat ? Xp_currSamp->Xp_mat->S_IID.c_str() : "0",
				Xp_currSamp->N_sex, S_pheno);
		else for (wsUint j = 0; j <= NCHR_SPECIES; j++) if (Cp_ped[j])
			Cp_ped[j]->fmt("%s %s %s %s %d %s", Xp_currSamp->S_FID.c_str(),
				Xp_currSamp->S_IID.c_str(),
				Xp_currSamp->Xp_pat ? Xp_currSamp->Xp_pat->S_IID.c_str() : "0",
				Xp_currSamp->Xp_mat ? Xp_currSamp->Xp_mat->S_IID.c_str() : "0",
				Xp_currSamp->N_sex, S_pheno);

		/* For all variants */
		if (!OPT_ENABLED(indel)) for (wsUint j = 0; j < N_variant; j++) {
			/* Determine target */
			wsUint N_target = 0;

			/* For --famsplit, allocate target by its FID */
			if (OPT_ENABLED(famsplit)) {
				N_target = Xm_fid2idx[Xp_currSamp->S_FID];
			} else if (OPT_ENABLED(split)) { /* Target could be not 0 ONLY IF --split */
				if (Xv_vrt[j].chr <= 0 || (wsUint)Xv_vrt[j].chr > NCHR_SPECIES)
					N_target = NCHR_SPECIES;
				else N_target = Xv_vrt[j].chr - 1;
			}
			cExporter*	Cp_currPED = Cp_ped[N_target];
			xVariant&	X_vrt = Xv_vrt[j];

			char al1 = X_vrt.al1;
			char al2 = X_vrt.al2;

			/* --out1234 */
			if (OPT_ENABLED(out1234)) {
				switch (X_vrt.al1) {
				case 'A': case 'a': al1 = '1'; break;
				case 'C': case 'c': al1 = '2'; break;
				case 'G': case 'g': al1 = '3'; break;
				case 'T': case 't': al1 = '4'; break;
				case '\0': al2 = '0'; break;
				default: halt("Invalid character [%c](%d) found, this dataset can't be converted with --out1234!", X_vrt.al1, X_vrt.al1);
					break;
				}
				switch (X_vrt.al2) {
				case 'A': case 'a': al2 = '1'; break;
				case 'C': case 'c': al2 = '2'; break;
				case 'G': case 'g': al2 = '3'; break;
				case 'T': case 't': al2 = '4'; break;
				case '\0': al2 = '0'; break;
				default: halt("Invalid character [%c](%d) found, this dataset can't be converted with --out1234!", X_vrt.al2, X_vrt.al2);
					break;
				}
			}

			if (IS_ASSIGNED(dosage)) {
				wsReal R_dsg = Ra_data[Xp_currSamp->N_idx][j];
				if (isMissingReal(R_dsg))
					Cp_currPED->fmt("  %s %s", S_misgeno, S_misgeno);
				else if (R_dsg < REAL_CONST(0.0)) {
					halt("Dosage [%s] of sample [%s] have invalid value [%g]",
						X_vrt.name, Xp_currSamp->S_IID.c_str(), R_dsg);
				} else if (R_dsg < REAL_CONST(0.5)) /* Major homo */
					Cp_currPED->fmt("  %c %c", al1, al1);
				else if (R_dsg < REAL_CONST(1.5)) /* Hetero */
					Cp_currPED->fmt("  %c %c", al2, al1);
				else if (R_dsg <= REAL_CONST(2.0)) /* Minor homo */
					Cp_currPED->fmt("  %c %c", al2, al2);
				else
					halt("Dosage [%s] of sample [%s] have invalid value [%g]",
						X_vrt.name, Xp_currSamp->S_IID.c_str(), R_dsg);
			} else switch (Na_data[Xp_currSamp->N_idx][j]) {
			case 0: /* Major homo */
				Cp_currPED->fmt("  %c %c", al1, al1); break;
			case 1: /* Hetero */
				Cp_currPED->fmt("  %c %c", al2, al1); break;
			case 2: /* Minor homo */
				Cp_currPED->fmt("  %c %c", al2, al2); break;
			case WISARD_NA:
				Cp_currPED->fmt("  %s %s", S_misgeno, S_misgeno); break;
			default:
				halt_fmt(WISARD_SYST_INVL_SNPDATA, Xv_vrt[j].name, Xp_currSamp->S_FID.c_str(),
					Xp_currSamp->S_IID.c_str(), Na_data[Xp_currSamp->N_idx][j]);
				break;
			}
		} else for (wsUint j = 0; j < N_variant; j++) {
			/* Determine target */
			wsUint N_target = 0;

			/* For --famsplit, allocate target by its FID */
			if (OPT_ENABLED(famsplit)) {
				N_target = Xm_fid2idx[Xp_currSamp->S_FID];
			} else if (OPT_ENABLED(split)) { /* Target could be not 0 ONLY IF --split */
				if (Xv_vrt[j].chr <= 0 || (wsUint)Xv_vrt[j].chr > NCHR_SPECIES)
					N_target = NCHR_SPECIES;
				else N_target = Xv_vrt[j].chr - 1;
			}
			cExporter *Cp_currPED = Cp_ped[N_target];

			switch (Na_data[Xp_currSamp->N_idx][j]) {
			case 0: /* Major homo */
				Cp_currPED->fmt("  %s %s", Xv_vrt[j].indel1, Xv_vrt[j].indel1); break;
			case 1: /* Hetero */
				Cp_currPED->fmt("  %s %s", Xv_vrt[j].indel1, Xv_vrt[j].indel2); break;
			case 2: /* Minor homo */
				Cp_currPED->fmt("  %s %s", Xv_vrt[j].indel1, Xv_vrt[j].indel2); break;
			case WISARD_NA:
				Cp_currPED->fmt("  %s %s", S_misgeno, S_misgeno); break;
			default:
				halt_fmt(WISARD_SYST_INVL_SNPDATA, Xv_vrt[j].name, Xp_currSamp->S_FID.c_str(),
					Xp_currSamp->S_IID.c_str(), Na_data[Xp_currSamp->N_idx][j]);
				break;
			}
		}
		i++;
		if (i != Xv_samp.end()) {
			if (OPT_ENABLED(famsplit))
				Cp_ped[Xm_fid2idx[Xp_currSamp->S_FID]]->put("\n");
			else for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) if (Cp_ped[i])
				Cp_ped[i]->put("\n");
		}
	}

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	/* Export MAP file */

	/* For --famsplit, export all variant info to all files */
	if (OPT_ENABLED(famsplit)) for (vVariant_it i=Xv_vrt.begin() ; i!=Xv_vrt.end() ; i++) {
		for (wsUint j=0 ; j<N_exporter ; j++)
			Cp_map[j]->fmt("%d %s %g %d\n", i->chr, i->name, i->gdist, i->pos);//  break;
	} else for (vVariant_it i=Xv_vrt.begin() ; i!=Xv_vrt.end() ;) {
		/* Determine target */
		wsUint N_target = 0;
		/* Target could be not 0 ONLY IF --split */
		if (OPT_ENABLED(split)) {
			if (i->chr <= 0 || (wsUint)i->chr > NCHR_SPECIES)
				N_target = NCHR_SPECIES;
			else N_target = i->chr - 1;
		}

		Cp_map[N_target]->fmt("%d %s %g %d\n", i->chr, i->name, i->gdist, i->pos);//  break;
		// 		switch ((*i).chr) {
		// 		case 23: /* X */
		// 		case 24: /* Y */
		// 			Cp_map->fmt("Y %s 0 0\n", (*i).name);  break;
		// 		case 25: /* Mt */
		// 			Cp_map->fmt("Mt %s 0 0\n", (*i).name);  break;
		// 		default:
		// 			Cp_map->fmt("%d %s 0 0\n", (*i).chr, (*i).name);  break;
		// 		}
		i++;
	}

	for (wsUint i=0 ; i<N_exporter ; i++) if (Cp_ped[i]) {
		delete Cp_map[i];
		delete Cp_ped[i];
	}
	DEALLOC(Cp_map);
	DEALLOC(Cp_ped);
}

} // End namespace ONETOOL
