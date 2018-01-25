#include <stdio.h>
#include <math.h>
#include "global/common.h"
#include "utils/util.h"
#include "global/option.h"
#include "input/stream.h"
#include "input/vcf.h"
#include "input/typed.h"

#define VCF_LINE_LENGTH	1024*1024

namespace ONETOOL {

typedef struct _xVCFfield {
	wsUint N_idxGT;
	wsUint N_idxDS;
	wsUint N_idxGQ;
	wsUint N_idxDP;
	wsUint N_idxGP;
} xVCFfield;

inline char* getVcfString(char *S_start, char *S_cur, char **Sp_next, wsUint N_len)
{
	/* If over the end of string */
	if ((S_cur - S_start) >= (int)N_len) return NULL;
	/* Find '\0' from S_cur */
	size_t L_str = strlen(S_cur);
	*Sp_next = S_cur + L_str + 1;

	return S_cur;
}

inline void _parseFormatField(char *Sp_fmt, xVCFfield &X_fdef)
{
	wsUint	N_idx=1;
	char	*a, *b;

	for (a=Sp_fmt-1 ; a ; N_idx++) {
		b = strchr(a+1, ':');

		/* Check out what of a+1 have */
		if (a[1] == 'G' && a[2] == 'T')
			X_fdef.N_idxGT = N_idx;
		/* Genotype dosage */
		else if (a[1] == 'D' && a[2] == 'S')
			X_fdef.N_idxDS = N_idx;
		/* Genotype quality */
		else if (a[1] == 'G' && a[2] == 'Q')
			X_fdef.N_idxGQ = N_idx;
		/* Read depth */
		else if (a[1] == 'D' && a[2] == 'P')
			X_fdef.N_idxDP = N_idx;
		/* Genotype posterior probability */
		else if (a[1] == 'G' && a[2] == 'P')
			X_fdef.N_idxGP = N_idx;
		a = b;
	}

	/* Cannot found GT field */
//	return 0xffffffff;
}

inline char* _getSubField(char *Sp_str, wsUint N_idxTarget)
{
	wsUint	N_idx=0;
	char	*a, *b;

	for (a=Sp_str-1 ; a ; N_idx++) {
		b = strchr(a+1, ':');
		if (b)
			*b = '\0';
		if (N_idx == N_idxTarget)
			return a+1;
		a = b;
	}

	/* Cannot found desired index */
	return NULL;
}

void typedValue(char* S_buf)
{
	char N_tbyte = S_buf[0];
	// Read atomic type (lower 4 bits)
	char N_atype = N_tbyte & 0x0f;
	char N_h4 = N_tbyte >> 4;
	if (N_h4 == 15) {
		int N_sz = 0;
		// array size, following 'integer' required
		switch (N_atype) {
		case 1:
			N_sz = (char)S_buf[1]; break;
		case 2:
			N_sz = (short)(*((short *)S_buf + 1)); break;
		case 3:
			N_sz = (int)(*((int *)S_buf + 1)); break;
		default:
			// ERR
			break;
		}
	} else if (N_h4 == 0) {
		// 
	} else {
		// size of following vector
	}

}

cVcfIO::cVcfIO(char *S_fn, xFileType X_type, char B_inDry)
{
	X_fType	= X_type;
	B_dry	= B_inDry;
	init(S_fn);
}

void cVcfIO::init(char *S_fn)
{
	if (X_fType == FT_DOSAGE && (OPT_ENABLED(phasedonly) || OPT_ENABLED(unphasedonly)))
		LOGwarn("--phasedonly/--unphasedonly will be disabled due to dosage input");
	if (S_fn == NULL) halt("SYSERR: NULL file name given");

	/* Init properties */
	cStrFile C_vcf(S_fn, "VCF file");

	/* Check whether it is BCF or not */
	char S_magic[5];
	C_vcf.read(S_magic, 5);
	if (memcmp(S_magic, "BCF", 3) == 0) {
		C_vcf.close();
		/* Re-open it as binary */
		C_vcf.open(S_fn, "Binary VCF file", 0, 1);
		/* Set as binary */
		B_bin = 1;
		/* Set the version */
		N_version = ((S_magic[3] - '0') << 4) | (S_magic[4] - '0');

		/* It seems to BCF file */
		LOGnote("[%s] identified as BCF version %d.%d\n", S_fn,
			N_version>>4, N_version&0xf);

		/* FIXME : It should be implemented */
		halt_fmt(WISARD_SYST_NULL_IMPLEMENT, "BCF retrieval");
	} else {
		B_bin = 0;
		C_vcf.close();
		C_vcf.open(S_fn, "VCF file");
	}

	/* Read alternative FAM file if exists */
	N_sample = 0;
	wsUint N_fSamp = 0;
	if (IS_ASSIGNED(fam)) {
		cStrFile C_fam(OPT_STRING(fam), "Alternative FAM file");
		char S_buf[1024];

		/* Construct sample data */
		for (N_fSamp=0 ; C_fam.gets(S_buf, 1024) ; N_fSamp++)
			_loadRealTimeFamData(S_buf, N_fSamp);
		LOG("Binstatus [%d], [%d] samples found from alternative FAM file\n", B_bin, N_fSamp);
	} else {
		/* Option checking related with sex */
		LOGwarn("Alternative FAM file is not assigned, sex and family-related option will be ignored\n");
		if (OPT_ENABLED(filnosex))	LOG("    --nosex ignored\n");
		if (OPT_ENABLED(filmale))	LOG("    --filmale ignored\n");
		if (OPT_ENABLED(filfemale))	LOG("    --filfemale ignored\n");
		if (IS_ASSIGNED(remfam))	LOG("    --remfam ignored\n");
		if (IS_ASSIGNED(selfam))	LOG("    --selfam ignored\n");
	}

	/* Initialize */
	char *Sp_buf = NULL;
	wsAlloc(Sp_buf, char, VCF_LINE_LENGTH);
	N_loadState = VCF_LOAD_META;
	N_sample	= 0xffffffff; // Undetermined
	N_idxGeno	= 0xffffffff; // Undetermined
	N_dataPos	= 0xffffffff; // Undetermined

	/* Find out whether it is gzipped or not */
	N_vrtOrig = N_variant = 0;
	if (B_bin) {
		/* Get the size of header */
		wsUint N_szHeader;
		C_vcf.read(&N_szHeader, 1);
		/* Read header */
		C_vcf.read(Sp_buf, N_szHeader);

		char *Sp_ptr = Sp_buf, *Sp_end = NULL;
		for (wsUint L=0 ; Sp_ptr && (Sp_end = strchr(Sp_ptr, '\n')) ; Sp_ptr=Sp_end?Sp_end+1:NULL) {
			int N_ret = procVcfLine(Sp_ptr, L+1);

			if (N_ret == 2) {}				
			else if (N_ret)
				halt("Load failed at line %d", L+1);
		}
		N_dataPos = C_vcf.tell();

		/* Load variants */
		while (!C_vcf.end()) {
			xVcfRecord X_rec;

			/* Read record header */
			C_vcf.read(&X_rec, sizeof(xVcfRecord));

			/* Compute remained length to the end of INFO */
			wsUint L_remained = X_rec.L_shared - sizeof(xVcfRecord) + (sizeof(wsUint)<<1);
			
			/* Read ID and FILTER */
			char *S_buf = NULL;
			wsAlloc(S_buf, char, L_remained);

			/* Read data from file */
			C_vcf.read(S_buf, L_remained);

			/* Get ID */
			char *Sp_ptr = NULL;
			char *Sp_ID = getVcfString(S_buf, S_buf, &Sp_ptr, L_remained);

			/* Get the list of alleles */
			vStr Sv_allele;
			LOOP (i, X_rec.N_allele) {
				char *Sp_allele = getVcfString(S_buf, Sp_ptr, &Sp_ptr, L_remained);
				Sv_allele.push_back(Sp_allele);
			}

			/* Read FILTER .*/
			char* Sp_filter = getVcfString(S_buf, Sp_ptr, &Sp_ptr, L_remained);

			/* Read INFO */
			vStr Sv_info;
			LOOP (i, X_rec.N_info) {
				char *Sp_allele = getVcfString(S_buf, Sp_ptr, &Sp_ptr, L_remained);
				Sv_info.push_back(Sp_allele);
			}

			/* Read genotype block (BGZF-compressed?) */
			LOOP(i, X_rec.N_fmt) {
				/* Read triplet */

				/* 1: INFO ID */
				char* Sp_idInfo = getVcfString(S_buf, Sp_ptr, &Sp_ptr, L_remained);
				/* Do lookup on info */
				/* Read typing byte */
				/* Read data */
			}

			/* Free memory */
			DEALLOC(S_buf);
		}
	} else for (wsUint L=0 ; C_vcf.gets(Sp_buf, VCF_LINE_LENGTH) ; L++) {
		int N_ret = procVcfLine(Sp_buf, L+1);

		if (N_ret == 2) {
			N_dataPos = L+1;
			/* Before data load */

			N_fSamp = (wsUint)Bv_filtSamp.size();
			if (X_fType == FT_VCF)
				Xa_geno = new vChar[N_sample];
			else if (X_fType == FT_DOSAGE)
				Xa_dosage = new vShort[N_sample];
			else halt("SYSERR: Non-permitted filetype for VCF handler [%d]", X_fType);

			Bv_filtVrt.clear();
			N_vrtOrig = 0;

			/* If no --fam, all samples are assumed to founders */
			if (!IS_ASSIGNED(fam)) {
				N_founder = N_sample;

				/* All set to founders */
				wsAlloc(Ba_isFounder, char, N_sample);
				memset(Ba_isFounder, 0x01, sizeof(char)*N_sample);
			}

			/* Allocate data space */
			_initMemories(1); // Without variant

			if (IS_ASSIGNED(fam)) {
				cStrFile C_fam(OPT_STRING(fam), "Alternative FAM file");
				char* S_buf = NULL;
				wsAlloc(S_buf, char, 1024);

				/* Read FAM file : main */
				N_founder = 0;
				char *Sp_tmp1 = NULL;
				char *Sp_tmp2 = NULL;
				/* should be --1 but forget check */
				wsUint N_0 = 0;
				char B_noPheno = 0;
				for (wsUint i=0, j=0 ; i<N_sampOrig ; i++) {
					if (C_fam.gets(S_buf, 1024) == 0)
						halt_fmt(WISARD_INVL_EOF_LINE, "FAM file", OPT_STRING(fam),
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
					//Xa_sampleV2.push_back(&X_sample);

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
									OPT_STRING(fam), "phenotype", Sp_tmp1, i+1);
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
				//		wsVec Ra_p = getPheno();
				//		for (int i=0 ; i<sizeSample() ; i++)
				//			LOG("%g\n", Ra_p[i]);

				DEALLOC(S_buf);
			}
		} else if (N_ret)
			halt("Load failed at line %d", L+1);

		if ((N_vrtOrig%1000) == 0)
			notice("[%d] variants scanned\r", N_vrtOrig);
	}
	if (X_fType == FT_VCF) {
		// If --vcf
		wsAlloc(Na_geno, char*, N_sample);
		wsUint sz = (wsUint)Xv_variant.size();
		LOOP(i, N_sample) {
			vChar& X_vgeno = Xa_geno[i];
			if (sz != X_vgeno.size())
				halt("SYSERR: Invalid Vgeno size");
			sseMalloc(Na_geno[i], char, sz);
			memcpy(Na_geno[i], X_vgeno.data(), sizeof(char)*sz);
			X_vgeno.clear();
		}
		DEALLOC(Xa_geno);
	} else {
		// If --dosage
		wsUint sz = (wsUint)Xv_variant.size();
		wsAlloc(Ra_dosage, wsFloat*, N_sample);

		LOOP (i, N_sample) {
			vShort& X_vgeno = Xa_dosage[i];
			if (sz != X_vgeno.size())
				halt("SYSERR: Invalid Vgeno size");
			sseMalloc(Ra_dosage[i], wsFloat, sz);
			LOOP (j, sz)
				if (X_vgeno[j] == (short)0xffff)
					Ra_dosage[i][j] = (wsFloat)WISARD_NA_REAL;
				else
					Ra_dosage[i][j] = X_vgeno[j] / CORR_CONST(10000.0);
			X_vgeno.clear();
		}
		DEALLOC(Xa_dosage);
	}
    Bv_filtSamp.clear();
    Bv_filtSamp.resize(N_sample, 0);
    Bv_filtVrt.clear();
    Bv_filtVrt.resize(Xv_variant.size(), 0);
    N_szFilteredVrt = 0;
    N_szFilteredSamp = 0;

    LOG("[%d/%d] variants and [%d/%d] samples loaded\n", Xv_variant.size(), Bv_filtVrt.size(),
        N_sample, N_fSamp);

	C_vcf.close();

	/* Register missing founder */
	_registMissingFounder();

	/* Reset missing status */
	vSampPtr &Xa_samps = getSample();
	FOREACH (vSampPtr_it, Xa_samps, i)
		if ((*i)->B_isComplete != -1) (*i)->B_isComplete = 0;

	/* Deallocate buffer */
	DEALLOC(Sp_buf);
}

inline char** loadVCFelements(char *Sp_val, wsUint *Np_var)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev		= 0;
	char		**Sp_ret	= NULL;
	if (Sp_val == NULL) {
		*Np_var = 0;
		return NULL;
	}
	char		*Sp_prev	= new char[strlen(Sp_val)+1];
	strcpy(Sp_prev, Sp_val);
	//strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') {
		*Np_var = 0;
		return NULL;
	}

	/* Get the number of given prevalences */
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ':');
		N_prev++;
	}
	// Since --mqls now requires the number of prevalence>2, it removed (130423)
	// 	if (N_prev > 2)
	// 		halt("Too much of prevalences (>2)");

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
	wsAlloc(Sp_ret, char*, N_prev);
	N_prev = 0;
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ':');
		if (b) *b = '\0';
		Sp_ret[N_prev++] = a+1;
	}

	*Np_var = N_prev;
	return Sp_ret;
}

inline void loadVCFinfo(char *Sp_val, mVar &Xa_info)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	//wsUint		N_prev		= 0;
	//char		**Sp_ret	= NULL;
	if (Sp_val == NULL) return;
	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') return;

	/* Get the number of given prevalences */
// 	for (a=Sp_prev-1 ; a ; a=b) {
// 		b = strchr(a+1, ':');
// 		N_prev++;
// 	}
	// Since --mqls now requires the number of prevalence>2, it removed (130423)
	// 	if (N_prev > 2)
	// 		halt("Too much of prevalences (>2)");

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
// 	MULTI_MALLOC(Sp_ret, char*, N_prev);
// 	N_prev = 0;
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ';');
		if (b) *b = '\0';
		/* Find = */
		char *c = strchr(a+1, '=');
		if (c) *(c++) = '\0';
		string S_name = a+1;
		xOperand X;
		X.X_type = OTP_NA;
		if (c) str2op(c, X);
		Xa_info.insert(make_pair(S_name, X));
	}
}

inline char** loadVCFalleles(char *Sp_val, wsUint *Np_var, char* Sp_ref)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev		= 0;
	char		**Sp_ret	= NULL;
	if (Sp_val == NULL) {
		*Np_var = 0;
		return NULL;
	}
	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') {
		*Np_var = 0;
		return NULL;
	}

	/* Get the number of given prevalences */
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		N_prev++;
	}
	// Since --mqls now requires the number of prevalence>2, it removed (130423)
	// 	if (N_prev > 2)
	// 		halt("Too much of prevalences (>2)");

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
	wsAlloc(Sp_ret, char*, N_prev);
	N_prev = 0;
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		if (b) *b = '\0';
		/* If dup. w/ Sp_ref do not count */
		if (Sp_ref &&
			((Sp_ref[1] && !stricmp(Sp_ref, a + 1)) ||
			(!Sp_ref[1] && !a[2] && Sp_ref[0] == a[1]))) {
			pverbose("ALT allele [%s] is duplicated with REF allele [%s]\n", a+1, Sp_ref);
		} else
			Sp_ret[N_prev++] = a+1;
	}

	*Np_var = N_prev;
	return Sp_ret;
}

inline char** loadVCFalleles(char *Sp_val, wsUint *Np_var, char S_ref=0)
{
	/* Get the structure of --prevalence */
	char		*a, *b;
	wsUint		N_prev		= 0;
	char		**Sp_ret	= NULL;
	if (Sp_val == NULL) {
		*Np_var = 0;
		return NULL;
	}
	char		*Sp_prev	= strdup(Sp_val);

	/* If --prevalence is not given, just return */
	if (Sp_prev[0] == '\0') {
		*Np_var = 0;
		return NULL;
	}

	/* Get the number of given prevalences */
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		N_prev++;
	}
	// Since --mqls now requires the number of prevalence>2, it removed (130423)
	// 	if (N_prev > 2)
	// 		halt("Too much of prevalences (>2)");

	/* Process prevalence string into real */
	//LOG("%d prevalences detected\n", N_prev);
	wsAlloc(Sp_ret, char*, N_prev);
	N_prev = 0;
	for (a=Sp_prev-1 ; a ; a=b) {
		b = strchr(a+1, ',');
		if (b) *b = '\0';
		/* If dup. w/ Sp_ref do not count */
		if (S_ref && (!a[2] && (S_ref == a[1]))) {
			pverbose("ALT allele [%s] is duplicated with REF allele [%c]\n", a+1, S_ref);
		} else
			Sp_ret[N_prev++] = a+1;
	}

	*Np_var = N_prev;
	return Sp_ret;
}

int cVcfIO::procVcfLine(char *Sp_buf, wsUint N_line)
{
	char	*a, *b;
//	xModel	X_model = _getModel();

	switch (N_loadState) {
	case VCF_LOAD_META:
		if (Sp_buf[0] == '#') {
			for (char *Sp = Sp_buf + strlen(Sp_buf) - 1 ; *Sp == '\r' || *Sp == '\n' || *Sp == ' ' || *Sp == '\t' ; *(Sp--) = '\0');
			
			if (Sp_buf[1] == 'C') { /* VCF_LOAD_HEADER section */
				/* == END OF meta */
				if (OPT_ENABLED(verbose)) {
					FOREACH(mDataIdx_it, Xm_seqINFO, i)
						LOG("INFO index %d : %s\n", i->second, i->first.c_str());
					FOREACH(mDataIdx_it, Xm_seqCONTIG, i)
						LOG("CONTIG index %d : %s\n", i->second, i->first.c_str());
				}
				N_loadState = VCF_LOAD_HEADER;
				return procVcfLine(Sp_buf, N_line);
			} else if (Sp_buf[1] != '#')
				halt_fmt(WISARD_INVL_VCFMETA, N_line, "Does not begin with double sharp(##)");
			_procVcfLine_meta(Sp_buf+2, N_line);

			//LOG("Metadata %s\n", Sp_buf);
		} else {
			halt_fmt(WISARD_INVL_VCFMETA, N_line, "Does not begin with sharp(#)");
			//LOG("Meta [%s]\n", Sp_buf);
		}
		break;
	case VCF_LOAD_HEADER:
		_procVcfLine_header(Sp_buf, N_line);
		return 2;
	case VCF_LOAD_SNP: {
			//wsUint	N_insIdx	= (wsUint)Xv_variant.size();
			xVariant	X_snp		= { 0, };
			char		B_filter	= 0;

			/* Load 8 mandatory columns
			 *
			 * WISARD will not use latter 3 columns, QUAL / FILTER / INFO
			 */
			getString(&Sp_buf, &a);	// Sp_buf = CHROM
			X_snp.chr = getChr(Sp_buf);

			getString(&a, &b);		// a = POS
			X_snp.pos = atoi(a);

			getString(&b, &a);		// b = ID
			if (b[0] == '.' && b[1] == '\0' || b[0] == '\0') {
				/* Allocate chr:pos style name if there is no name */
				wsAlloc(X_snp.name, char, 256);
				sprintf(X_snp.name, "%s:%u", getChrName2(X_snp.chr), X_snp.pos);
			} else {
				wsAlloc(X_snp.name, char, strlen(b)+1);
				strcpy(X_snp.name, b);
			}
			B_filter = _filterVariantPosition(X_snp);
			string S_snpID = (string)b;

			getString(&a, &b);		// a = REF
			char *Sp_ref = a;
			char B_aSNV = 0;
			char B_bSNV = 0;
			if (a[1]) {
				/* Not this variant is SNP, pass */
				if (!OPT_ENABLED(indel)) {
					if (OPT_ENABLED(snvonly)) B_filter |= 1;
					else halt_fmt(WISARD_CANT_GETINDEL_WO_INDEL, N_line);
				} else
					X_snp.indel1 = strdup(a);
			} else if (a[1] != '-') {
				B_aSNV = 1;
				if (OPT_ENABLED(indel)) X_snp.indel1 = strdup(a);
				else X_snp.al1 = a[0];
			} else B_aSNV = 1;

			getString(&b, &a);		// b = ALT
			/* Split them into multiple */
			wsUint N_alt = 0;
			char **Na_alts = loadVCFalleles(b, &N_alt, Sp_ref);
			/* No definite triallelic allowed */
			if (N_alt > 2)
				B_filter |= 1;
			else if (Na_alts[0][1] || (N_alt == 2 && Na_alts[1][1])) {
				/* At least one indel */
				if (!OPT_ENABLED(indel)) {
					/* Not this variant is SNP or biallelic, pass */
					if (OPT_ENABLED(snvonly)) B_filter |= 1;
					else halt_fmt(WISARD_CANT_GETINDEL_WO_INDEL, N_line);
				}
			} else if (a[2] != '-') {
				B_bSNV = 1;
				/* SNP */
				if (OPT_ENABLED(indel)) {
					if (N_alt == 2) {
						free(X_snp.indel1);
						X_snp.indel2 = strdup(Na_alts[0]);
						X_snp.indel2 = strdup(Na_alts[1]);
					} else
						X_snp.indel2 = strdup(Na_alts[0]);
				} else {
					if (N_alt == 2) {
						X_snp.al1 = Na_alts[0][0];
						X_snp.al2 = Na_alts[1][0];
					} else X_snp.al2 = Na_alts[0][0];
				}
			} else B_bSNV = 1;

			if (N_alt > 0)
				free(Na_alts[0]);
			DEALLOC(Na_alts);

			/* --acgt */
			code2acgt(X_snp, OPT_ENABLED(indel));

			/* One-byte indel filtering */
			if (!(B_aSNV & B_bSNV) && OPT_ENABLED(snvonly))
				B_filter |= 1;
			if ((B_aSNV & B_bSNV) && OPT_ENABLED(indelonly))
				B_filter |= 1;

			/* Skip latter 3 columns */
			getString(&a, &b);		// a = QUAL
			wsReal R_qual = (wsReal)atof(a);
			if ((IS_ASSIGNED(filqual) && isInRange(OPT_RANGE(filqual), R_qual)) ||
				(IS_ASSIGNED(incqual) && isInRange(OPT_RANGE(incqual), R_qual)))
				/* Not this variant passed QUAL filter, pass */
				B_filter |= 1;
			getString(&b, &a);		// b = FILTER
			if (OPT_ENABLED(vcfqc) && stricmp(b, "PASS"))
				/* Not this variant passed QC under --vcfqc, pass */
				B_filter |= 1;

			getString(&a, &b);		// a = INFO
			/* Assign values */
			if (IS_ASSIGNED(filvariant) || IS_ASSIGNED(incvariant) || IS_ASSIGNED(filmaf) || IS_ASSIGNED(incmaf)) {
				/* Fetch info */
				mVar Xa_info;
				loadVCFinfo(a, Xa_info);

				if ((IS_ASSIGNED(filmaf) || IS_ASSIGNED(incmaf)) && mapIsExist(Xa_info, "AF")) {
					xOptRange& X_rng = IS_ASSIGNED(filmaf) ? OPT_RANGE(filmaf) : OPT_RANGE(incmaf);
					xOperand& X_af = Xa_info["AF"];
					float R_af = X_af.R_val > 0.5f ? 1 - X_af.R_val : X_af.R_val;
					if (IS_ASSIGNED(filmaf) && isInRange(X_rng, R_af) || IS_ASSIGNED(incmaf) && !isInRange(X_rng, R_af))
						B_filter |= 1;
				}

				if (IS_ASSIGNED(filvariant) || IS_ASSIGNED(incvariant)) {
					xOperation &X_eorig = IS_ASSIGNED(filvariant) ? OPT_EXPR(filvariant) : OPT_EXPR(incvariant);
					xOperation *Xp_op = X_eorig.clone();

					bool B_satisfy = false;
					/* Assign values */
					FOREACH (mVar_it, Xa_info, i)
						Xp_op->assign(i->first.c_str(), &(i->second));
					/* And evaluate... */
					xOperand *Xp_ret = understand(Xp_op);
					if (Xp_ret->X_type != OTP_LOGICAL)
						halt("Option --%s have non-logical statement!");
					/* Check 'filtering' condition is satisfied */
					B_satisfy = Xp_ret->B_val;
					if (IS_ASSIGNED(incvariant)) B_satisfy = !B_satisfy;
					/* Filter if satisfied */
					if (B_satisfy) B_filter |= 1;

					/* Clean remains */
					remOperand(Xp_ret);
					remOperand(Xp_op->X_op1);
					remOperand(Xp_op->X_op2);
					delete Xp_op;

					FOREACH(mVar_it, Xa_info, q)
						remOperand(&(q->second));
				}
			}

			/* Seems to have FORMAT field, but check */
			xVCFfield X_fdef ={ 0, };
			if (b && b[0]) {
				if (N_sample == 0xffffffff)
					halt("This file indicated to have no data, but data at line %d", N_line);

				getString(&b, &a);
				_parseFormatField(b, X_fdef);
				N_idxGeno = X_fType == FT_VCF ?
					(X_fdef.N_idxGT ? X_fdef.N_idxGT-1 : 0xffffffff) :
					(X_fdef.N_idxDS ? X_fdef.N_idxDS-1 : 0xffffffff);

				if (N_idxGeno == 0xffffffff)
					/* Not this variant have genotype/dosage field, pass */
					B_filter |= 1;
			} else
				B_filter |= 1;

			/* Stop and go next if this variant should be filtered */
			N_vrtOrig++;
			if (B_filter)
				break;

			/* Load genotyping filtering */
			xOperation *Xp_corig = IS_ASSIGNED(filgeno) ?
				&OPT_EXPR(filgeno) : (IS_ASSIGNED(incgeno) ?
				&OPT_EXPR(incgeno) : NULL);
			xOperation *Xp_cond = Xp_corig ? Xp_corig->clone() : NULL;
			char B_skip = 0;

			wsUint N_sampleOrig = (wsUint)Bv_filtSamp.size();
			for (wsUint i=0, _i=0 ; i<N_sampleOrig ; i++, a=b) {
				bool	B_isNA	= false;
				wsUint	N_elem	= 0;
				if (!a)
					halt_fmt(WISARD_INVL_NSAMPLEOFSNP, X_snp.name, N_sample, i);
				getString(&a, &b);

				/* Sample filtering check */
				if (Bv_filtSamp[i]) continue;

				char **Sa_elem = loadVCFelements(a, &N_elem);

				if (Xp_cond) {
					/* GQ found*/
					if (X_fdef.N_idxGQ && N_elem >= X_fdef.N_idxGQ) {
						xOperand *X_op = str2op(Sa_elem[X_fdef.N_idxGQ - 1]);
						Xp_cond->assign("GQ", X_op);
					}
					else Xp_cond->assign("GQ", &X_NA);

					/* DP found */
					if (X_fdef.N_idxDP && Xp_cond && N_elem >= X_fdef.N_idxDP) {
						xOperand *X_op = str2op(Sa_elem[X_fdef.N_idxDP-1]);
						Xp_cond->assign("DP", X_op);
					}
					else Xp_cond->assign("DP", &X_NA);

					/* If Xp_cond exists, evaluate */
					xOperand *X_res = understand(Xp_cond);
					B_isNA = X_res->B_val;
					if (IS_ASSIGNED(incgeno)) B_isNA = !B_isNA;
				}

				/* Get genotype field */
				if (N_elem <= N_idxGeno)
					halt("Data lacked at [%d]th sample on [%d]th variant", i + 1, N_line);
				char *Sp_geno = Sa_elem[N_idxGeno];
				if (Sp_geno == NULL)
					halt("Data lacked at [%d]th sample on [%d]th variant", i + 1, N_line);

				/* If input is VCF genotype (--vcf) */
				if (X_fType == FT_VCF) {
					char N_al1 = Sp_geno[0]-'0';
					if (!Sp_geno[1]) {
						/* Should be Y or MT */
						if (isMtChromosome(X_snp) || isYChromosome(X_snp)) {}
						else halt("Variant [%s] (chr %s, pos %d) is non-haploid but "
							"contains haploid genotype!", X_snp.name,
							getChrName2(X_snp.chr), X_snp.pos);

						/* When N_alt == 2, 0 is not allowed */
						if (N_alt == 2) {
							if (N_al1 == 0) {
								if (N_szFilteredVrt == 100)
									LOGnote("Number of filtered variants is now more than 100, rest are omitted\n");
								else if (N_szFilteredVrt < 100)
									LOG("Variant [%s] (chr %s, pos %d) is identified as tri-allelic\n",
										X_snp.name, getChrName2(X_snp.chr), X_snp.pos);

								/* Pop already inserted genotypes */
								for (wsUint k=0 ; k<_i ; k++) Xa_geno[_i].pop_back();

								B_skip = 1;
								/* Should be decreased */
								Bv_filtVrt.push_back(MFR_TRIALLELIC);
								//N_fileVrt++;
								N_szFilteredVrt++;
								//N_variant--;
								free(Sa_elem[0]);
								DEALLOC(Sa_elem);
								break;
							} else N_al1--;
						}

						/* N_line = N_line th variant */

						/* Retrieve genotype */
						if (Sp_geno[0] == '.' || B_isNA) {
							Xa_geno[_i].push_back(WISARD_NA);
							N_missGeno++;
						} else {
							wsUint N_ret = N_al1;
							if (N_ret > 1) {
								halt("Too large allele value [%d]", N_ret);
							}
							/* Model-dependent processing */
							// 					switch (X_model) {
							// 					case MD_DOMINANT:
							// 						N_ret = (N_ret+2)>>2; break;
							// 					case MD_RECESSIVE:
							// 						N_ret = (N_ret+1)>>1; break;
							// 					case MD_MULTIPLICATIVE:
							// 						N_ret = (N_ret+2)&0x05; break;
							// 					default: break;
							// 					}
							Xa_geno[_i].push_back((char)N_ret);
						}
					} else {
						char B_phased = Sp_geno[1] == '|';
						char N_al2 = Sp_geno[2]-'0';
						if (!B_isNA) {
							if ((OPT_ENABLED(unphasedonly) && B_phased) ||
								(OPT_ENABLED(phasedonly) && !B_phased))
								B_isNA = 1;
						}
						/* When N_alt == 2, 0 is not allowed */
						if (N_alt == 2) {
							if (N_al1 == 0 || N_al2 == 0) {
								if (N_szFilteredVrt == 100)
									LOGnote("Number of filtered variants is now more than 100, rest are omitted\n");
								else if (N_szFilteredVrt < 100)
									LOG("Variant [%s] (chr %s, pos %d) is identified as tri-allelic!\n",
										X_snp.name, getChrName2(X_snp.chr), X_snp.pos);

								/* Pop already inserted genotypes */
								for (wsUint k=0 ; k<_i ; k++) Xa_geno[k].pop_back();

								B_skip = 1;
								/* Should be decreased */
								Bv_filtVrt.push_back(MFR_TRIALLELIC);
								//N_fileVrt++;
								N_szFilteredVrt++;
								//N_variant--;
								free(Sa_elem[0]);
								DEALLOC(Sa_elem);
								break;
							} else {
								N_al1--;
								N_al2--;
							}
						}

						/* N_line = N_line th variant */

						/* Retrieve genotype */
						if (Sp_geno[0] == '.' || Sp_geno[2] == '.' || B_isNA) {
							Xa_geno[_i].push_back(WISARD_NA);
							Xa_sampleV2[_i]->B_isComplete = -1;
							N_missGeno++;
						} else {
							wsUint N_ret = N_al1 + N_al2;
							if (N_ret > 2) {
								halt("Too large allele value [%d]", N_ret);
							}
							/* Model-dependent processing */
							// 					switch (X_model) {
							// 					case MD_DOMINANT:
							// 						N_ret = (N_ret+2)>>2; break;
							// 					case MD_RECESSIVE:
							// 						N_ret = (N_ret+1)>>1; break;
							// 					case MD_MULTIPLICATIVE:
							// 						N_ret = (N_ret+2)&0x05; break;
							// 					default: break;
							// 					}
							Xa_geno[_i].push_back((char)N_ret);
						}
					}
				} else if (X_fType == FT_DOSAGE) {
					char* Sp_tmp = NULL;

					if (Sp_geno[0] == '.') { // NA
						Xa_dosage[_i].push_back((short)0xffff);
						Xa_sampleV2[_i]->B_isComplete = -1;
						N_missGeno++;
					} else {
						double R_ret = strtod(Sp_geno, &Sp_tmp);
						if (!Sp_tmp || Sp_tmp[0])
							halt("Invalid dosage data [%s] in [%d]th sample on [%d]th variant", Sp_geno, i + 1, N_line);

						// 0~2 -> 0~20000
						Xa_dosage[_i].push_back((short)round(R_ret*10000.0));
					}
				}

				free(Sa_elem[0]);
				DEALLOC(Sa_elem);

				_i++;
			}
			if (B_skip == 0) {
				N_variant++;
				/* Insert variant info */
				Xv_variant.push_back(X_snp);
				Bv_filtVrt.push_back(0);
			}
			else
				DEALLOC(X_snp.name);
		} break;
	default:
		halt("SYSERR: Inappropriate state[%d] given", N_loadState);
	}
	/* Successfully terminated */
	return 0;
}

void cVcfIO::_procVcfLine_header(char *Sp_buf, wsUint N_line)
{
	char *a, *b;
	wsUint	i;
	const char *Sa_cols[] = { "#CHROM", "POS", "ID", "REF", "ALT",
		"QUAL", "FILTER", "INFO" };
	/* Check integrity of column of first 8 columns */
	for (i=0,a=Sp_buf ; i<8 ; i++,a=b) {
		getString(&a, &b);

		if (stricmp(a, Sa_cols[i]))
			halt_fmt(WISARD_INVL_VCFMANDTRY_COLNAME, i+1, Sa_cols[i], a);
	}

	/* Check whether there are extra columns */
	wsUint N_filtSamp = 0, N_idxSamp = 0;
	if (IS_ASSIGNED(fam))
		Bv_filtSamp.clear();
	if (a && a[0]) {
		getString(&a, &b);

		if (stricmp(a, "FORMAT"))
			halt("Invalid VCF format, FORMAT field expected, got [%s]", a);

		/* Get the number of samples included in VCF file
			* and make filtering */
		for (i=0,a=b ; a ; a=b,i++) {
			xSample*	Xp_samp		= NULL;
			char		B_inserted	= 0;
			char		B_dyn		= 0;
			getString(&a, &b);

			if (IS_ASSIGNED(fam)) {
				/* Find it from alternative sample */
				mSamp_it X_it = Xm_iid2data.find(a);

				if (X_it == Xm_iid2data.end()) {
					LOG("[%s] exist in data, but not exist in alternative FAM file, will ignored\n", a);
					continue;
				}

				Xp_samp = &(X_it->second);
				pverbose("Sample [%s::%s] found in VCF file\n",
					Xp_samp->S_FID.c_str(), Xp_samp->S_IID.c_str());
			} else {
				Xp_samp					= new xSample;
				//memset(Xp_samp, 0x00, sizeof(xSample));
				Xp_samp->N_oriIdx		= (int)Bv_filtSamp.size();
				Xp_samp->S_FID.assign(a);
				Xp_samp->S_IID.assign(a);
				Xp_samp->B_isAssigned	= 0;
				Xp_samp->B_isComplete	= 0;
				Xp_samp->B_isMissing	= 0;
				Xp_samp->B_isProband	= 0;
				Xp_samp->N_idTwin		= 0;
				Xp_samp->N_idx			= SAMP_NODATA;
				Xp_samp->N_isVisited	= 0;
				Xp_samp->N_missGeno		= 0;
				Xp_samp->N_sex			= 0;
				Xp_samp->Xp_mat			= NULL;
				Xp_samp->Xp_pat			= NULL;
				Xp_samp->N_grpFst		= -1;
				B_dyn					= 1;
			}
			/* Mark it assigned */
			string str2srch = a;
			char B_filtered = N_filterSample && isSampleFiltered(str2srch, str2srch,
				((wsUint)Bv_filtSamp.size()-N_sample)>5);
			if (B_filtered) {
				Bv_filtSamp.push_back(1);
				N_filtSamp++;
				Xp_samp->N_idx = SAMP_NODATA;
			} else {
				Bv_filtSamp.push_back(0);
				Xp_samp->N_idx	= N_idxSamp++;
				B_inserted		= 1;
			}
			if (B_inserted) {
				Xp_samp->B_isAssigned	= 1;
				Xp_samp->B_isMissing	= 0;
			} else
				Xp_samp->B_isMissing	= 1;
			Xm_iid2data.insert(make_pair((string)a, *Xp_samp));
			if (B_inserted)
				Xa_sampleV2.push_back(&Xm_iid2data[a]);
			if (B_dyn)
				delete Xp_samp;
		}
		N_sampOrig = (wsUint)Bv_filtSamp.size();
	}
	LOG("%d samples found in header, %d filtered\n", N_sampOrig,
		N_filtSamp);
	N_sample = (wsUint)Xa_sampleV2.size();
	if ((N_sampOrig - N_filtSamp) != N_sample)
		halt("SYSERR: VCF header loading error N_sampOrig [%d] N_filtSamp [%d] N_sample [%d]",
			N_sampOrig, N_filtSamp, N_sample);

	/* Move state to VCF_LOAD_SNP
		* because this state only requires ONE line */
	N_loadState	= VCF_LOAD_SNP;
	N_variant	= 0;
}

// Input is double-sharp(##) removed character buffer
void cVcfIO::_procVcfLine_meta(char *Sp_buf, wsUint N_line)
{
	string S_raw = Sp_buf;
	char* S_fieldKey	= new char[4096];
	char* S_cont		= new char[8192];
	int N_ret = sscanf(Sp_buf, "%[^=]=<%[^>]>", S_fieldKey, S_cont);
	if (N_ret == 2) {
		map<string,char*> Xm_data;
		char* Sp_tmp = strtok(S_cont, ",");
		do {
			// Sp_tmp = one segment
			// Find =
			char* Sp_div = strchr(Sp_tmp, '=');
			if (Sp_div) {
				char* Sp_key = Sp_tmp;
				char* Sp_val = Sp_div+1;
				*Sp_div = '\0';
				Xm_data[Sp_key] = Sp_val;
			}
		} while (Sp_tmp = strtok(NULL, ","));

		// According to field key
		if (!stricmp(S_fieldKey, "INFO")) {
			// Halt if don't have ID key
			if (Xm_data.find("ID") == Xm_data.end())
				halt("INFO field at line [%d] does not have its ID", N_line);
			// 'INFO' field should be numbered by its sequence of appearance
			Xm_seqINFO[Xm_data["ID"]] = (int)(Xm_seqINFO.size() + 1);
		} else if (!stricmp(S_fieldKey, "contig")) {
			// Halt if don't have ID key
			if (Xm_data.find("ID") == Xm_data.end())
				halt("CONTIG field at line [%d] does not have its ID", N_line);
			// 'contig' field should be numbered by its sequence of appearance
			Xm_seqCONTIG[Xm_data["ID"]] = (int)(Xm_seqCONTIG.size() + 1);
		}
	}
	//
	delete [] S_fieldKey;
	delete [] S_cont;
	Sa_meta.push_back(S_raw);
}


void xSeq::load()
{
	char *Sp_buf = NULL;
	wsAlloc(Sp_buf, char, 8192);
	S_buf.clear();

	cStrFile &X_f = Cp_fa->getHandle();
	X_f.seek(N_pos, SEEK_SET);
	while (X_f.gets(Sp_buf, 8192)) {
		if (Sp_buf[0] == '>') break;
		/* Do trimming */
		for (char *Sp_t=Sp_buf+strlen(Sp_buf)-1 ; Sp_buf<=Sp_t && (*Sp_t == '\r' || *Sp_t == '\n' || *Sp_t == ' ' || *Sp_t == '\t') ; Sp_t--) *Sp_t = '\0';
		/* DO uppercasing */
		for (char *Sp_t=Sp_buf ; *Sp_t ; Sp_t++)
			*Sp_t = (*Sp_t) & (~0x20);
		S_buf += Sp_buf;
	}

	DEALLOC(Sp_buf);
}

typedef enum _xVcfAnnoRes {
	AR_OK,
	AR_INVRANGE,	/* Invalid chromosome position */
	AR_NOREFAL,		/* Non-reference allele */
} xVcfAnnoRes;

void xSeq::unload()
{
	S_buf.clear();
}


cStrFile& cFa::getHandle()
{
	return X_f;
}

vSeq& cFa::getSeqs()
{
	return Xv_seq;
}

bool cFa::init(wsStrCst S_fn, wsStrCst S_desc/*=NULL*/)
{
	char *Sp_buf = NULL;
	wsAlloc(Sp_buf, char, 8192);

	/* Open */
	X_f.open(S_fn, S_desc, 1);
	if (X_f.isFailed()) return false;

	/* Find all >'s */
	for (int L=1 ; X_f.gets(Sp_buf, 8192) ; L++) {
		if (Sp_buf[0] == '>') {
			xSeq	X_ent;
			X_ent.Cp_fa	= this;
			X_ent.S_desc	= (string)(Sp_buf+1);
			char *Sp_k = Sp_buf+1, *tmp;
			getString(&Sp_k, &tmp);
			LOG("Chromosome [%s] found in FASTA\n", Sp_k);
			X_ent.N_pos	= X_f.tell();

			Xv_seq.push_back(X_ent);
		}
	}

	DEALLOC(Sp_buf);
	return true;
}

class cSeqFa
{
	bool	B_init;
	cFa		*Xv_fas[MAX_NCHR+1];
	xSeq	*Xv_pFa[MAX_NCHR+1];
public:
	xSeq* getSeq(int N_chr) {
		if (N_chr < 1 || N_chr > (int)NCHR_SPECIES)
			halt("SYSERR: Invalid chromosome definition [%d]", N_chr);
		return Xv_pFa[N_chr];
	}
	cSeqFa() {
		B_init = false;
	}
	cSeqFa(wsStrCst S_fns) {
		char	*Sp_fns		= strdup(S_fns);
		char	**Sp_path	= NULL;
		wsUint	N_path		= 0;
		bool	B_seq		= false;

		/* Init */
		memset(Xv_pFa, 0x00, sizeof(xSeq*)*(MAX_NCHR+1));

		/* Anyway, try to open this */
		cStrFile X_seq(Sp_fns, "Reference FASTA path file", 1);
		if (!X_seq.isFailed()) {
// 			vector<string>	Xv_fns;
// 			char *Sp_buf = NULL;
// 			MULTI_MALLOC(Sp_buf, char, 8192);
// 
// 			/* Read and replace them */
// 			for (wsUint L=0 ; X_seq.gets(Sp_buf, 8192) ; L++) {
// 				if (Sp_buf[0] && L>=MAX_CHR2)
// 					halt("Too many lines over maximum number of chromosomes "
// 						"found in reference FASTA path file [%s]", Sp_fns);
// 				Xv_fns.push_back((string)Sp_buf);
// 			}
// 			free(Sp_fns);
// 			/* If # of paths is less than autosome size, warning */
// 			if (Xv_fns.size() < MAX_AUTOCHR)
// 				LOGwarn("Reference FASTA path file contains lines[%d] "
// 					"less than the number of autosome[%d]\n", Xv_fns.size(),
// 					MAX_AUTOCHR);
// 
// 			MULTI_MALLOC(Sp_path, char*, Xv_fns.size());
// 			wsUint I=0;
// 			FOREACHDO (vStr_it, Xv_fns, i, I++)
// 				Sp_path[I] = strdup(i->c_str());
// 			N_path	= I;
// 			B_seq	= true;

 			N_path	= 1;
			wsAlloc(Sp_path, char*, 1);
			Sp_path[0] = strdup(Sp_fns);
// 			B_seq	= true;
// 			DEALLOC(Sp_buf);
		} else {
			/* otherwise, get the path by plain and store anyway,
			 * # of chromosome integrity will be confirmed */
			Sp_path = loadStringValues(Sp_fns, &N_path);
		}

		/* Now load FA files */
		if (B_seq) {
			for (wsUint i=1 ; i<=N_path ; i++) {
				Xv_fas[i] = new cFa(Sp_path[i], "Reference FASTA file");
				vSeq &Xv_seqs = Xv_fas[i]->getSeqs();
				if (Xv_seqs.size() != 1)
					halt("Reference FASTA file [%s] for chromosome [%s] "
						"have invalid number of sequences[%d]", Sp_path[i],
						getChrName2(i), Xv_seqs.size());
				Xv_pFa[i] = &(Xv_seqs[0]);
			}
		} else for (wsUint i=1 ; i<=N_path ; i++) {
			Xv_fas[i] = new cFa(Sp_path[i-1], "Reference FASTA file");

			/* Check which seq */
			vSeq &Xv_seqs = Xv_fas[i]->getSeqs();
			FOREACH (vSeq_it, Xv_seqs, it) {
				/* Check chr */
				char *Sp_desc = strdup(it->S_desc.c_str()), *a = NULL;
				getString(&Sp_desc, &a);

				/* (chr) or chr(chr) */
				int chr = !memcmp("chr", Sp_desc, 3) ? getChr(Sp_desc+3)
					: getChr(Sp_desc);
				if (chr >= 1 && chr <= (int)NCHR_SPECIES) {
					if (Xv_pFa[chr] == NULL)
						Xv_pFa[chr] = &(*it);
					else
						halt("Reference FASTA file [%s] found for chromosome [%s] "
						"but it is duplicated", Sp_path[i-1], getChrName2(chr));
				} else {
					char *label = strdup(it->S_desc.c_str()), *tmp = NULL;
					getString(&label, &tmp);
					LOGwarn("Unrecognized chromosome label [%s] in file [%s]\n",
						label, Sp_path[i-1]);
					free(label);
				}
			}
		}

		B_init = true;
	}
};

void cIO::exportVCF(char B_bin/*=0*/)
{
	wsUint		N_sample	= sizeSample();
	//wsUint		N_mk		= getMarkerSize();
	vVariant&	Xa_mk		= getVariant();
	char		B_makeFam = 1;// (N_sample != N_fnder);
	cExporter**	Cp_fam	= NULL;
	cExporter**	Cp_vcf	= NULL;
	if (B_makeFam) {
		wsCalloc(Cp_fam, cExporter*, MAX_NCHR+1);
	}
	wsCalloc(Cp_vcf, cExporter*, MAX_NCHR+1);

	/* Cannot apply --outphenoonly */
	if (OPT_ENABLED(outphenoonly))
		LOGwarn("BED file does not support --outphenoonly, will be ignored\n");

	/* if --ref */
	char*		Ba_ref		= NULL;
	cExporter*	H_anerr		= NULL;
	wsUint		N_annoErr	= 0;
	//cFa*		Cp_fa	= NULL;
	if (IS_ASSIGNED(ref)) {
		cTimer t;
		t.start();

		wsCalloc(Ba_ref, char, N_variant);
		LOG("Reference sequence assigned, VCF annotation now based on this....\n");

		/* Get ref */
		cSeqFa X_seq(OPT_STRING(ref));
		for (int i=1 ; i<=(int)NCHR_SPECIES ; i++) {
			/* Trying to get sequence */
			xSeq *Xp_seq = X_seq.getSeq(i);
			if (Xp_seq == NULL) {
				LOG("Chromosome [%s] does not have sequence, skip...\n",
					getChrName2(i));
				continue;
			} else
				notice("Annotating chromosome [%s] (ref. size %d)...\r",
					getChrName2(i), Xp_seq->S_buf.size());

			/* If available, load sequence */
			Xp_seq->load();
			string	S_buf	= Xp_seq->S_buf;
			wsStrCst	Sp_buf	= S_buf.c_str();
			size_t	N_sz	= S_buf.size();

			/* Do annotation */
			wsUint J=0;
			FOREACHDO (vVariant_it, Xa_mk, j, J++) {
				if (j->chr != i) continue;
				
				/* Check range of pos */
				if (N_sz < (wsUint)j->pos) {
					if (!H_anerr) {
						H_anerr = cExporter::summon("anno.error.lst");
						headerVariant(H_anerr);
						H_anerr->put("	REASON\n");
					}
					entryVariant(H_anerr, *j);
					H_anerr->put("	INVALID_POSITION\n");
					N_annoErr++;
					/* If Ba_ref[J] do not have 0x40, then it means it is
					 * a reference allele that not equivalent to two alleles
					 * in the dataset */
					Ba_ref[J] = 'N' & (~0x40);
// 					halt("Variant [%s] at chr [%s] pos [%d] is invalid, max [%d]",
// 						j->name, getChrName2(j->chr), j->pos, N_sz);
					continue;
				}
				char N_ref = Sp_buf[j->pos - 1]; // pos is 1-base
				char S_ref[2] = { N_ref, 0 };

				/* Check reference allele is exist */
				char B_match1 = OPT_ENABLED(indel) ? !strcmp(j->indel1, S_ref) : j->al1 == N_ref;
				char B_match2 = OPT_ENABLED(indel) ? !strcmp(j->indel2, S_ref) : j->al2 == N_ref;
				if (!B_match1 && !B_match2) {
					if (!H_anerr) {
						H_anerr = cExporter::summon("anno.error.lst");
						headerVariant(H_anerr);
						H_anerr->put("	REASON\n");
					}
					entryVariant(H_anerr, *j);
					H_anerr->put("	REFALLELE_NOTEXIST\n");
// 					halt("Variant [%s] at chr [%s] pos [%d] does not have "
// 						"reference allele [%c]", j->name, getChrName2(j->chr),
// 						j->pos, N_ref);
					N_annoErr++;
					/* If Ba_ref[J] do not have 0x40, then it means it is
					 * a reference allele that not equivalent to two alleles
					 * in the dataset */
					Ba_ref[J] = N_ref & (~0x40);
				} else
					/* Should be flipped? */
					Ba_ref[J] = B_match2 ? N_ref|0x80 : N_ref;
			}

			/* Unload sequence */
			Xp_seq->unload();
		}

		LOG("Annotation took [%s]\n", t.getReadable());
	}
	if (N_annoErr)
		LOGwarn("[%d] errornous variant found during annotation\n", N_annoErr);

	if (IS_ASSIGNED(split)) {
		/* Make instance of EXPOTER ONLY IF exists */
		char S_ext[64];
		for (wsUint i=0 ; i<NCHR_SPECIES ; i++) if (Ba_chrExists[i]) {
			wsStrCst S_chr = getChrName2(i + 1);
			if (B_bin)
				sprintf(S_ext, "chr%s.bcf", S_chr);
			else
				sprintf(S_ext, "chr%s.vcf", S_chr);
			Cp_vcf[i] = cExporter::summon(S_ext, 0,
				B_bin?ET_BGZF:(OPT_ENABLED(gz)?ET_GZIP:ET_PLAIN));
			if (B_makeFam) {
				sprintf(S_ext, "chr%s.fam", S_chr);
				Cp_fam[i] = cExporter::summon(S_ext);
			}
		}

		/* For Undetermined IF exists */
		if (Ba_chrExists[NCHR_SPECIES]) {
			Cp_vcf[NCHR_SPECIES] = cExporter::summon(B_bin?"chrUn.bcf":"chrUn.vcf");
			if (B_makeFam) Cp_fam[NCHR_SPECIES] = cExporter::summon("chrUn.fam");
		}

		if (B_makeFam)
			LOGoutput("Final dataset is exported with [%s.chr***.%s] [%s.chr***.fam]\n",
				OPT_STRING(out), B_bin?"bcf":"vcf", OPT_STRING(out));
		else
			LOGoutput("Final dataset is exported with [%s.chr***.%s]\n",
				OPT_STRING(out), B_bin?"bcf":"vcf");
	} else {
/**/	Cp_vcf[0]	= cExporter::summon(B_bin?"bcf":"vcf"); /* CHECKED */
		
		if (B_makeFam) {
			Cp_fam[0]	= cExporter::summon("fam"); /* CHECKED */
			LOG("Final dataset is exported with [%s.%s] [%s.fam]\n",
				OPT_STRING(out), B_bin?"bcf":"vcf", OPT_STRING(out));
		} else
			LOG("Final dataset is exported with [%s.%s]\n", OPT_STRING(out),
				B_bin?"bcf":"vcf");
	}

	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

	/* If dataset is independent, do not make FAM file */
	if (B_makeFam) {
		wsUint j=0;
		FOREACHDO (vSampPtr_it, Xa_sampleV2, i, j++) {
			xSample *Xp_currSamp = *i;
			wsUint k = Xp_currSamp->N_idx;

			/* Determining phenotype */
			char	S_pheno[128];
			getPhenoCode(k, S_pheno, Sp_case, Sp_ctrl);

			/* For phenotype, export header */
			for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) if (Cp_fam[k])
				Cp_fam[k]->fmt("%s %s %s %s %d %s\n",
					Xp_currSamp->S_FID.c_str(), Xp_currSamp->S_IID.c_str(),
					Xp_currSamp->Xp_pat?Xp_currSamp->Xp_pat->S_IID.c_str():"0",
					Xp_currSamp->Xp_mat?Xp_currSamp->Xp_mat->S_IID.c_str():"0",
					Xp_currSamp->N_sex, S_pheno);
		}
	}

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	/* For BCF, print MAGIC first */
	char magic[5] ={ 'B', 'C', 'F', '\2', '\1' };
	if (B_bin) for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) if (Cp_vcf[k])
		Cp_vcf[k]->write(magic, 5);

	/* Export header if available */
	wsUint N_insert = 1;
	if (!stricmp(getFormat(), "io.vcf") ||
		!stricmp(getFormat(), "io.bcf")) {
		cVcfIO	*Cp_this	= (cVcfIO *)this;
		vStr	&Sa_meta	= Cp_this->getMeta();

		/* Insert 'source' field to metadata unless LAST source filed is
		 * NOT WISARD */
		for (vStr::reverse_iterator it=Sa_meta.rbegin() ; it!=Sa_meta.rend() ;
			it++)
			if (it->compare(0, 7, "source=") == 0) {
				if (it->compare(7, 6, "WISARD") != 0) {
					it--;
					Sa_meta.insert(it.base(), (string)"source=WISARD");
					N_insert = 0;
				}
				break;
			}

		/* if --makebcf, make dictionary */
		if (OPT_ENABLED(makebcf)) {
			char S_tmp[16];
			string S_meta = "dictionary=";
			wsUint I = 0;
			/* INFO */
			FOREACHDO (mDataIdx_it, Cp_this->Xm_seqINFO, i, I++) {
				if (I == 0)
					sprintf(S_tmp, "s%d", I);
				else
					sprintf(S_tmp, ",s%d", I);
				S_meta += S_tmp;
			}
			Sa_meta.push_back(S_meta);
		}

		for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) if (Cp_vcf[k]) {
			FOREACH (vStr_it, Sa_meta, it)
				Cp_vcf[k]->fmt("##%s\n", it->c_str());
		}

		/* FIXME: --makebcf requires implmentation */
	}
	if (N_insert) {
		WISARD_NOW(X_start);
		for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) if (Cp_vcf[k]) {
			Cp_vcf[k]->fmt("##fileformat=VCFv4.1\n");
			Cp_vcf[k]->fmt("##fileDate=%04d%02d%02d\n", X_start.tm_year + 1900,
				X_start.tm_mon+1, X_start.tm_mday);
			Cp_vcf[k]->fmt("##source=WISARD\n");
			Cp_vcf[k]->fmt("##WISARDversion=\"" VERSION_WISARD "\"\n");
			Cp_vcf[k]->fmt("##WISARDcommand=\"%s\"\n", OPT_STRING(cmdLine));
			Cp_vcf[k]->fmt("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
		}
	}

	/* Print out mandatory columns */
	const char *Sa_cols[] = { "#CHROM", "POS", "ID", "REF", "ALT",
		"QUAL", "FILTER", "INFO", "FORMAT" };
	for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) if (Cp_vcf[k]) {
		LOOP (i, 8) Cp_vcf[k]->fmt("%s\t", Sa_cols[i]);
		Cp_vcf[k]->fmt("%s", Sa_cols[8]);

		/* Print out sample IDs */
		vSampPtr &Xa_samp = getSample();
		FOREACH (vSampPtr_it, Xa_samp, it)
			Cp_vcf[k]->fmt("	%s", (*it)->S_IID.c_str());
		Cp_vcf[k]->put("\n");
	}

	/* Now print out more info */
	vVariant&	Xa_snp		= getVariant();
	char**		Na_geno		= getGenotype();
	wsUint		N_idxSNP	= 0;
	if (OPT_ENABLED(indel)) FOREACHDO (vVariant_it, Xa_snp, it, N_idxSNP++) {
		/* Determine target */
		wsUint N_target = 0;
		/* Target could be not 0 ONLY IF --split */
		if (OPT_ENABLED(split)) {
			if (it->chr <=0 || it->chr > (int)NCHR_SPECIES)
				N_target = NCHR_SPECIES;
			else N_target = it->chr - 1;
		}

		wsStrCst	Sp_chr	= getChrName2(it->chr);
		char	N_rInfo	= Ba_ref ? Ba_ref[N_idxSNP] : 0xff;
		char*	Sp_al1	= Ba_ref ? (N_rInfo&0x80 ? it->indel2 : it->indel1) : it->indel1;
		char*	Sp_al2	= Ba_ref ? (N_rInfo&0x80 ? it->indel1 : it->indel2) : it->indel2;
		/* Two non-ref */
		if (!(N_rInfo & 0x40)) {
			Cp_vcf[N_target]->fmt("%s	%d	%s	%c	%s,%s	.	PASS	.	GT",
				Sp_chr, it->pos, it->name, N_rInfo|0x40, Sp_al1, Sp_al2);
		} else
			Cp_vcf[N_target]->fmt("%s	%d	%s	%s	%s	.	PASS	.	GT",
				Sp_chr, it->pos, it->name, Sp_al1, Sp_al2&&Sp_al2[0]?Sp_al2:".");

		if (Ba_ref) LOOP (i, N_sample) {
			char N_geno	= Na_geno[i][N_idxSNP];
			char N_al1	= N_rInfo&0x80 ? N_geno>0 : N_geno>>1;
			char N_al2	= N_rInfo&0x80 ? N_geno>>1 : N_geno>0;

			/* Two non-ref, then +1 */
			if (!(N_rInfo & 0x40)) {
				N_al1++;
				N_al2++;
			}
			if (isMissing(N_geno))
				Cp_vcf[N_target]->put("	./.");
			else
				Cp_vcf[N_target]->fmt("	%d/%d", N_al1, N_al2);
		} else LOOP (i, N_sample) {
			char N_geno = Na_geno[i][N_idxSNP];
			if (isMissing(N_geno))
				Cp_vcf[N_target]->put("	./.");
			else
				Cp_vcf[N_target]->fmt("	%d/%d", N_geno>>1, N_geno>0);
		}
		Cp_vcf[N_target]->put("\n");

		if ((N_idxSNP%100)==0)
			notice("[%d/%d] variants exported...\r", N_idxSNP, N_variant);
	} else  FOREACHDO (vVariant_it, Xa_snp, it, N_idxSNP++) {
		/* Determine target */
		wsUint N_target = 0;
		/* Target could be not 0 ONLY IF --split */
		if (OPT_ENABLED(split)) {
			if (it->chr <=0 || it->chr > (int)NCHR_SPECIES)
				N_target = NCHR_SPECIES;
			else N_target = it->chr - 1;
		}

		wsStrCst	Sp_chr	= getChrName2(it->chr);
		char	N_rInfo	= Ba_ref ? Ba_ref[N_idxSNP] : 0xff;
		char	N_al1	= Ba_ref ? (N_rInfo&0x80 ? it->al2 : it->al1) : it->al1;
		char	N_al2	= Ba_ref ? (N_rInfo&0x80 ? it->al1 : it->al2) : it->al2;
		/* Two non-ref */
		if (!(N_rInfo & 0x40)) {
			Cp_vcf[N_target]->fmt("%s	%d	%s	%c	%c,%c	.	PASS	.	GT",
				Sp_chr, it->pos, it->name, N_rInfo|0x40, N_al1, N_al2);
		} else
			Cp_vcf[N_target]->fmt("%s	%d	%s	%c	%c	.	PASS	.	GT",
				Sp_chr, it->pos, it->name, N_al1, N_al2?N_al2:'.');

		if (Ba_ref) LOOP (i, N_sample) {
			char N_geno	= Na_geno[i][N_idxSNP];
			char N_al1	= Ba_ref[N_idxSNP]&0x80 ? N_geno>0 : N_geno>1;
			char N_al2	= Ba_ref[N_idxSNP]&0x80 ? N_geno>1 : N_geno>0;

			/* Two non-ref, then +1 */
			if (!(N_rInfo & 0x40)) {
				N_al1++;
				N_al2++;
			}
			if (isMissing(N_geno))
				Cp_vcf[N_target]->put("	./.");
			else
				Cp_vcf[N_target]->fmt("	%d/%d", N_al1, N_al2);
		} else LOOP (i, N_sample) {
			char N_geno = Na_geno[i][N_idxSNP];
			if (isMissing(N_geno))
				Cp_vcf[N_target]->put("	./.");
			else
				Cp_vcf[N_target]->fmt("	%d/%d", N_geno>1, N_geno>0);
		}
		Cp_vcf[N_target]->put("\n");

		if ((N_idxSNP%100)==0)
			notice("[%d/%d] variants exported...\r", N_idxSNP, N_variant);
	}
	LOG("[%d/%d] variants exported\n", N_variant, N_variant);


	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) if (Cp_vcf[i]) {
		delete Cp_vcf[i];
		if (B_makeFam) delete Cp_fam[i];
	}
	DEALLOC(Cp_vcf);
	if (B_makeFam) DEALLOC(Cp_fam);
}

} // End namespace ONETOOL
