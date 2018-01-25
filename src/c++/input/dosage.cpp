#include <limits>
#include "global/common.h"
#include "global/option.h"
#include "input/dosage.h"
#include "input/stream.h"
#include "output/exporter.h"

namespace ONETOOL {

typedef vector<wsFloat*>		vArrReal;
typedef vArrReal::iterator	vArrReal_it;
cDosageIO::cDosageIO(wsStrCst S_fn, char B_inDry)
{
	X_type = DSG_UNKNOWN;
	B_dry = B_inDry;
	_init(S_fn);
}

cDosageIO::~cDosageIO()
{

}

void cDosageIO::_exportDosageDist(wsReal R_minDsg, wsReal R_maxDsg)
{
	/* Windows size : [R_minDsg, R_maxDsg] with precision 0.001 */
	int N_maxDsg = (int)(R_maxDsg * REAL_CONST(1000.0) + REAL_CONST(0.5));
	int N_minDsg = (int)(R_minDsg * REAL_CONST(1000.0) + REAL_CONST(0.5));

	LOG("Dosage distribution info will be exported to\n");
	LOG("   [%s.entire.dist.dsg] that contains entire distribution\n", OPT_STRING(out));
	LOG("   [%s.sample.dist.dsg] that contains sample-wise distribution\n", OPT_STRING(out));
//	LOG("	[%s.variant.disg.dsg] that contains variant-wise distribution\n", OPT_STRING(out));

	/* Count'em and allocate */
	wsUint	N_intv	= (N_maxDsg - N_minDsg) + 1;
	/* Entire count */
	__int64	*Na_cnt		= NULL;
	/* Sample-wise count */
	wsUint	*Na_cntSmp	= NULL;
	/* Variant-wise count */
//	wsUint	**Na_cntSnp	= NULL;

	/* Allocate memory */
	sseCalloc(Na_cnt, __int64, N_intv);
	sseCalloc(Na_cntSmp, wsUint, N_intv);
// 	sseMalloc(Na_cntSnp, wsUint*, N_SNP);
// 	for (wsUint i=0 ; i<N_SNP ; i++) {
// 		sseCalloc(Na_cntSnp[i], wsUint, N_intv);
// 	}

	cExporter *Cp_smp = cExporter::summon("sample.dist.dsg");
	LOGoutput("Sample-wise dosage value distribution is exported to "
		"[%s.sample.dist.dsg]\n", OPT_STRING(out));
	/* Print header */
	for (int i=N_minDsg ; i<N_maxDsg ; i++) Cp_smp->fmt("%g	", i/1000.0);
	Cp_smp->fmt("%g\n", N_maxDsg/1000.0);

	/* For dataset */
	if (B_isDataComplete) for (wsUint i=0 ; i<N_sample ; i++) {
		memset(Na_cntSmp, 0x00, sizeof(wsUint)*N_intv);
		for (wsUint j=0 ; j<N_variant ; j++) {
			/* Get index */
			int N_val = (int)(Ra_data[i][j] * REAL_CONST(1000.0) + REAL_CONST(0.5));
			int N_idx = N_val - N_minDsg;

			/* Increase count */
			Na_cnt[N_idx]++;
			Na_cntSmp[N_idx]++;
//			Na_cntSnp[j][N_idx]++;
		}

		/* Print sample-wise result */
		Cp_smp->fmt("%s:%s", Xa_sampleV2[i]->S_FID.c_str(), Xa_sampleV2[i]->S_IID.c_str());
		for (wsUint j=0 ; j<N_intv ; j++)
			Cp_smp->fmt("	%d", Na_cntSmp[j]);
		Cp_smp->put("\n");

		if ((i%10) == 0)
			notice("%d/%d samples computed...\r", i, N_sample);
	} else  for (wsUint i=0 ; i<N_sample ; i++) {
		memset(Na_cntSmp, 0x00, sizeof(wsUint)*N_intv);
		for (wsUint j=0 ; j<N_variant ; j++) {
			if (isMissingReal(Ra_data[i][j])) continue;

			/* Get index */
			int N_val = (int)(Ra_data[i][j] * REAL_CONST(1000.0) + REAL_CONST(0.5));
			int N_idx = N_val - N_minDsg;

			/* Increase count */
			Na_cnt[N_idx]++;
			Na_cntSmp[N_idx]++;
//			Na_cntSnp[j][N_idx]++;
		}

		/* Print sample-wise result */
		Cp_smp->fmt("%s:%s", Xa_sampleV2[i]->S_FID.c_str(), Xa_sampleV2[i]->S_IID.c_str());
		for (wsUint j=0 ; j<N_intv ; j++)
			Cp_smp->fmt("	%d", Na_cntSmp[j]);
		Cp_smp->put("\n");

		if ((i%10) == 0)
			notice("%d/%d samples computed...\r", i, N_sample);
	}
	LOG("Dosage distribution has successfully calculated!\n");

	/* Export data */
	LOGoutput("Entire dosage value distribution is exported to "
		"[%s.entire.dist.dsg]\n", OPT_STRING(out));
	cExporter *Cp_ent = cExporter::summon("entire.dist.dsg");
//	cExporter *Cp_snp = cExporter::summon("variant.dist.dsg");

	/* Print entire */
	for (int i=N_minDsg,j=0 ; i<N_maxDsg ; i++,j++)
		Cp_ent->fmt("%g	%d\n", i/1000.0, Na_cnt[j]);

	/* Print SNP-wise */
// 	for (int i=N_minDsg ; i<N_maxDsg ; i++) Cp_snp->fmt("%g	", i/1000.0);
// 	Cp_snp->fmt("%g\n", N_maxDsg/1000.0);
// 	for (wsUint i=0 ; i<N_SNP ; i++) {
// 		Cp_snp->fmt("%s", Xa_SNP[i].name);
// 		for (wsUint j=0 ; j<N_intv ; j++)
// 			Cp_snp->fmt("	%g", Na_cntSnp[i][j]);
// 		Cp_snp->put("\n");
// 	}

	/* Deallocate */
	sseFree(Na_cnt);
	sseFree(Na_cntSmp);
// 	for (wsUint i=0 ; i<N_SNP ; i++) {
// 		sseFree(Na_cntSnp[i]);
// 	}
// 	sseFree(Na_cntSnp);
}

void cDosageIO::_init(wsStrCst S_fn)
{
	/* Forcely turn on --indep
	 * 150520 removed */
// 	if (!OPT_ENABLED(indep)) {
// 		OPTION().assign("indep", "1");
// 		OPTION().FORCE_OPT_NUMBER(indep);
// 	}

	/* At first, assume it is OTHER type */
	X_type = DSG_OTHER;

	/* Separator helper */
	char		*a			= NULL;
	char		*b			= NULL;
	/* Custom separator */
	char		*Sp_delim	= OPT_STRING(sep);
	if (Sp_delim) LOG("Data field delimiter is set to [%c]\n", *Sp_delim);
	/* Reading buffer */
	char		*S_buf		= NULL;
	wsAlloc(S_buf, char, 1024*1024);
	/* Is the FAM file provided? */
	char		B_isFAM		= 0;
	/* Export dosage distribution? */
	char		B_dsgdist	= OPT_ENABLED(dsgdist);
	wsReal		R_minDsg	= REAL_CONST(999999.0);
	wsReal		R_maxDsg	= REAL_CONST(-999999.0);
	/* Intermediate buffer for dosage */
	vArrWord	Xa_dsg;
	vArrReal	Xa_expr;

	/* Dosage/expression? */
	wsStrCst		S_type		= IS_ASSIGNED(dosage) ? "Dosage" : "Expression";

	/* Init properties */
	cStrFile C_dosage(S_fn, "Dosage file");

	/* Read alternative FAM file if exists */
	N_sample	= 0;
	N_sampOrig	= 0;
	if (IS_ASSIGNED(fam)) {
		cStrFile C_fam(OPT_STRING(fam), "Alternative FAM file");

		/* Construct sample data */
		for (N_sampOrig=0 ; C_fam.gets(S_buf, 1024) ; N_sampOrig++)
			_loadRealTimeFamData(S_buf, N_sampOrig);

		B_isFAM = TRUE;
	} else {
		/* Option checking related with sex */
		LOG("Alternative FAM file is not assigned, sex-related option will be ignored\n");
		if (OPT_ENABLED(filnosex))
			LOG("    --nosex ignored\n");
		if (OPT_ENABLED(filmale))
			LOG("    --filmale ignored\n");
		if (OPT_ENABLED(filfemale))
			LOG("    --filfemale ignored\n");
	}

	/* Read dosage file HEADER */
	if (C_dosage.gets(S_buf, 1024*1024) == NULL)
		halt("Failed to read file, maybe the empty file?");

	/* Read header */
	/* Separate..., assume that the dosage file is variant-wise */
	if (Sp_delim)	getStringDelim(&S_buf, &a, *Sp_delim);
	else			getString(&S_buf, &a);

	char	B_alleleInfo	= -1;
	char	B_vrtWise		= 1;
	char	B_header		= 1;

	/* Check additional file to confirm the filetype */ {
		char* S_newFn = new char[4096];
		sprintf(S_newFn, "%s_info", S_fn);
		FILE* H_test = fopen(S_newFn, "r");
		if (H_test) {
			fclose(H_test);
			B_alleleInfo = 1;
			LOG("%s file [%s] seems IMPUTE2 format\n", S_type, S_fn);
			X_type = DSG_IMPUTE2;
			if (!IS_ASSIGNED(fam))
				halt_fmt(WISARD_CANT_DO_WO_SOMEOPT, "IMPUTE2 format dosage data", "--fam");
			B_header = 0;
		}
		delete [] S_newFn;
	}


	/* S_buf check */
	if (X_type == DSG_UNKNOWN || X_type == DSG_OTHER) {
		if (!stricmp(S_buf, "marker")) { /* BEAGLE */
			B_alleleInfo = 1;
			LOG("%s file [%s] seems BEAGLE format\n", S_type, S_fn);
			X_type = DSG_BEAGLE;

			/* Skip two fields */
			getString(&a, &b);
			getString(&b, &a);
		} else
			B_alleleInfo = 0;
	}

	/*
	 *
	 * Perform FORMAT-wise additional information load
	 *
	 */
	strMap M_qual;

	if (X_type == DSG_BEAGLE) {
		/* Data will be complete */
		B_isDataComplete = 1;

		char	S_r2[512];
		char	*Sp_r2buf = NULL;
		otherExt(S_fn, "r2", S_r2);
		wsAlloc(Sp_r2buf, char, 1024);
		cStrFile	C_r2(S_r2, "r^2 file of BEAGLE", 1);
		if (C_r2.isFailed()) {
			otherExt(S_fn, "r2", S_r2, 2);
			C_r2.open(S_r2, "r^2 file of BEAGLE", 1);
		}
		char		*aa, *bb;

		if (!C_r2.isFailed()) for (wsUint L=1 ; C_r2.gets(Sp_r2buf, 1024) ; L++) {
			getString(&Sp_r2buf, &aa);
			string	S_vrt	= (string)Sp_r2buf;
			/* aa = quality */
			getString(&aa, &bb);
			wsReal	R_qual	= (wsReal)str2dbl(aa, &bb);
			if (bb[0] && stricmp(bb, "NaN")) halt("Invalid BEAGLE r^2 format found at variant [%s] "
				"from line [%d]", Sp_r2buf, L);

			/* Insert this entry UNLESS this entry does have value */
			if (R_qual == R_qual)
				/* 10e-5 precision */
				M_qual.insert(make_pair(S_vrt, (wsUint)(R_qual*REAL_CONST(100000.0))));
		}
		
		DEALLOC(Sp_r2buf);
	}

	wsUint		N_fVrt = 0;
	wsStrCst	Sp_misGeno = "NA";

	/* Pass first field */
	if (Sp_delim) getStringDelim(&a, &b, *Sp_delim);
	else getString(&a, &b);

	N_vrtOrig = 0;
	N_variant = 0;
	/* If a == "DOSE", then it is minimac DOSAGE */
	if (!stricmp(a, "DOSE") || !stricmp(a, "ML_DOSE")) {
		LOG("%s file [%s] seems minimac format\n", S_type, S_fn);
		X_type		= DSG_MINIMAC;
		B_vrtWise	= 0;
		B_header	= 0;
		if (!B_isFAM) N_pheno = 0;

		char S_infoFn[512];
		otherExt(S_fn, "mlinfo", S_infoFn);
		cStrFile C_info(S_infoFn, "Variant information of minimac", 1);
		if (C_info.isFailed()) {
			// Count the variants using the first line
			// b = start
			char *c = NULL;
			while (b) {
				getString(&b, &c);
				N_vrtOrig++;
				N_variant++;
				N_fVrt++;
				Bv_filtVrt.push_back(0);
				b = c;
			}
			LOGwarn("Variant information of given dosage file does not exists, [%d] variants found in the data\n",
				N_vrtOrig);
		} else {
			char* S_buf2 = new char[1024*1024];
			char *aa = NULL, *bb = NULL, *cc = NULL, *dd = NULL;
			C_info.gets(S_buf2, 1024*1024); // Skip header
			while (C_info.gets(S_buf2, 1024*1024)) {
				getString(&S_buf2, &aa); // S_buf2 = name
				getString(&aa, &bb); // aa = al1
				getString(&bb, &cc); // bb = al2
				getString(&cc, &dd); // cc = freq1

				if ((aa[1] || bb[1]) && !OPT_ENABLED(indel))
					halt("--indel is required since the given dosage file contains indel(s)");

				/* variant information it is! */
				xVariant	X_currVrt = { 0, };
				char		B_filtered = 0;

				/* Check freq1 */ {
					char *tmp = NULL;
					double R_freq = strtod(cc, &tmp);
					// Flip
					if (R_freq < 0.5) {
						tmp = aa;
						aa = bb;
						bb = tmp;
					}
				}

				/* For convenience, major allele=1 and minor allele=0 */
				if (OPT_ENABLED(indel)) {
					X_currVrt.indel1 = strdup(aa);
					X_currVrt.indel2 = strdup(bb);
					X_currVrt.al1 = '\0';
					X_currVrt.al2 = '\0';
				} else {
					X_currVrt.indel1 = NULL;
					X_currVrt.indel2 = NULL;
					X_currVrt.al1 = aa[0];
					X_currVrt.al2 = bb[0];
				}
				X_currVrt.filter = 0;
				X_currVrt.gdist = W0;
				/* Chromosome number is undeterministic, because RAW file does not contain such information */
				X_currVrt.chr = 0;
				X_currVrt.pos = 0;

				/* Retrieve variant name */
				wsAlloc(X_currVrt.name, char, strlen(a) + 1);
				strcpy(X_currVrt.name, S_buf2);

				B_filtered = _filterVariantPosition(X_currVrt);

				if (B_filtered == 1)
					Bv_filtVrt.push_back(1);
				else {
					Bv_filtVrt.push_back(0);
					Xv_variant.push_back(X_currVrt);
					N_variant++;
				}
				N_vrtOrig++;
			}
		}
	} else if (!stricmp(S_buf, "FID") || !stricmp(S_buf, "IID")) {
		eStr Xe_key;
		LOG("%s file [%s] seems raw format\n", S_type, S_fn);

		if (!IS_ASSIGNED(misgeno)) {
			OPTION().assign("misgeno", "NA");
			OPTION().FORCE_OPT_STRING(misgeno);
		}

		/* --nofid? */
		if (!stricmp(S_buf, "IID")) {
			OPTION().assign("nofid", "1", 1);
			OPTION().FORCE_OPT_NUMBER(nofid);
		}
		Xe_key.insert(S_buf);
		Xe_key.insert(a);

		/* Move to next */
		for (a=b ; a ; a=b) {
			getString(&a, &b);

			/* Should be FID/IID/PAT/MAT/SEX/PHENOTYPE */
			if (stricmp(a, "FID") && stricmp(a, "IID") &&
				stricmp(a, "PAT") && stricmp(a, "MAT") &&
				stricmp(a, "SEX") && stricmp(a, "PHENOTYPE"))
				break;
			Xe_key.insert(a);
		}

		if (Xe_key.find("FID") == Xe_key.end()) {
			/* --nofid */
		}
		if (Xe_key.find("IID") == Xe_key.end()) {
			/* No IID is not possible */
			halt("No IID column found for raw dosage file!");
		}
		/* --nosex */
		if (Xe_key.find("SEX") == Xe_key.end()) {
			OPTION().assign("nosex", "1", 1);
			OPTION().FORCE_OPT_NUMBER(nosex);
		}
		/* --nopheno */
		if (Xe_key.find("PHENOTYPE") == Xe_key.end()) {
			OPTION().assign("nopheno", "1", 1);
			OPTION().FORCE_OPT_NUMBER(nopheno);
		}
		/* no PAT | MAT */
		eStr_it X_pat = Xe_key.find("PAT");
		eStr_it X_mat = Xe_key.find("MAT");
		if (X_pat == Xe_key.end() || X_mat == Xe_key.end()) {
			if (X_pat == Xe_key.end() && X_mat == Xe_key.end()) {
				/* --noparent */
				OPTION().assign("noparent", "1", 1);
				OPTION().FORCE_OPT_NUMBER(noparent);
			} else
				/* One PAT | MAT exists, then error */
				halt("Column PAT and MAT should exist or not exist simultaneously!");
		}

		X_type		= DSG_RAW;
		B_vrtWise	= 0;
		B_header	= 1;
	}

	/* Set all to missing WHO are marked as NOT filtered
	 * ONLY IF there is header */
	wsUint N_newNsamp = 0;
	if (B_isFAM && B_header) {
		wsUint I=0;
		FOREACHDO (mSamp_it, Xm_iid2data, it, I++) {
			if (!it->second.B_isMissing && it->second.N_oriIdx != SAMP_NODATA) {
				it->second.B_isMissing				= -1;
				Bv_filtSamp[it->second.N_oriIdx]	= 1;
				it->second.N_idx					= -1;
			}
		}
	}	

	/* For all column headers IF B_header */
	wsUint N_xx = 0;
	if (B_header) for ( ; a ; ) {
		if (B_vrtWise) {
			/* Check S_tok exists in the IID map that read from --fam IF --fam assgined */
			if (B_isFAM) {
				mSamp_it it = Xm_iid2data.find(a);

				if (it == Xm_iid2data.end())
					halt("Sample [%s] does not exists in the assigned FAM file", a);
				else if (it->second.B_isMissing == -1) { /* NOT filtered BUT undetermined */

					it->second.B_isMissing = 0;
					it->second.N_idx = N_xx++;
					Bv_filtSamp[it->second.N_oriIdx] = 0;
					N_newNsamp++;
				}
			} else {
				/* Otherwise, insert */
				char S_fam[512];
				/* Make fam buffer by the IID */
				sprintf(S_fam, "%s	%s	0	0	0	-9", a, a);
				/*xSample *Xp_samp = */_loadRealTimeFamData(S_fam, N_sampOrig);
			}

			/* Increase index */
			if (!B_isFAM) N_sampOrig++;

			/* Go to the next ptr */
			a = b;
			if (a) {
				if (Sp_delim) getStringDelim(&a, &b, *Sp_delim);
				else getString(&a, &b);
			}
		} else {
			/* variant information it is! */
			xVariant	X_currVrt = { 0, };
			char		B_filtered = 0;

			/* For convenience, major allele=1 and minor allele=0 */
			X_currVrt.al1 = '0';
			X_currVrt.al2 = '1';
			X_currVrt.indel1 = NULL;
			X_currVrt.indel2 = NULL;
			X_currVrt.filter = 0;
			X_currVrt.gdist = W0;
			/* Chromosome number is undeterministic, because RAW file does not contain such information */
			X_currVrt.chr = 0;
			X_currVrt.pos = 0;

			/* Retrieve variant name */
			wsAlloc(X_currVrt.name, char, strlen(a) + 1);
			strcpy(X_currVrt.name, a);

			B_filtered = _filterVariantPosition(X_currVrt);

			if (B_filtered == 1)
				Bv_filtVrt.push_back(1);
			else {
				Bv_filtVrt.push_back(0);
				Xv_variant.push_back(X_currVrt);
				N_variant++;
			}

			/* Go to the next ptr */
			a = b;
			if (a) {
				if (Sp_delim) getStringDelim(&a, &b, *Sp_delim);
				else getString(&a, &b);
			}
			N_vrtOrig++;
		}
	}

	if (X_type == DSG_RAW && B_vrtWise == 0) {
		size_t N_pos = C_dosage.tell();

		/* Read pedigree info */
		for (char* S_buf=NULL ; (S_buf=C_dosage.ugets()) ; ) {
			_loadRealTimeFamData(S_buf, N_sampOrig++);
			DEALLOC(S_buf);
		}

		C_dosage.seek(N_pos, SEEK_SET);
	}
	
	/* --fam and dosage file have header, it can be compared */
	if (B_isFAM && B_header) {
		LOG("FAM have [%d] samples but actual data have [%d] matched samples,"
			"so # of samples now adjusted\n", N_sample, N_newNsamp);
		N_sample = N_newNsamp;
	}

	/* Make complete family data */
	_registMissingFounder();

	/* Init memories */
	_initMemories(TRUE);
	if (!B_isFAM) {
		N_pheno = 0;
	}

	/* Init filterings */
	_initFilterings();

	/* Get the data-start position */
	//size_t N_dataPos = C_dosage.tell();
	B_isDataComplete = 0;
	N_missGeno = 0;

	/* Get the count of markers and read data
	 * Assume the range of dosage value is -3.2767 ~ 3.2768 */

	/* Read FAM file : main */
	if (B_isFAM) {
		cStrFile	C_fam(OPT_STRING(fam), "Alternative FAM file");

		/* should be --1 but forget check */
		wsUint N_0 = 0;
		N_founder = 0;
		wsUint j = 0;
		for (wsUint i=0 ; i<N_sampOrig ; i++) {
			if (C_fam.gets(S_buf, 1024) == 0)
				halt("Unexpected end of file found at line %d", i+1);
			string	S_currIID;

			/* Pass FID/IID/PAT/MAT/SEX */
			getString(&S_buf, &a);	// Sp_buf  = FID
			getString(&a, &b);	// Sp_tmp1 = IID
			S_currIID = (string)a;
			xSample	&X_sample = Xm_iid2data[S_currIID];
			getString(&b, &a);	// Sp_tmp2 = PAT
			getString(&a, &b);	// Sp_tmp1 = MAT
			getString(&b, &a);	// Sp_tmp2 = SEX

			/* If current samples should be filtered, do not retrieve anything */
			if (Bv_filtSamp[i]) {
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
			getString(&a, &b);
			Ra_pheno[0][j] = (wsReal)atoi(a);

			/* Phenotype value is not control nor case */
			if (Ra_pheno[0][j] != (wsReal)OPT_NUMBER(phenoCtrl) &&
				Ra_pheno[0][j] != (wsReal)OPT_NUMBER(phenoCase)) {

				if (!strcmp(a, OPT_STRING(mispheno))) {
					/* It is missing phenotype */
					Ra_pheno[0][j] = WISARD_NA_REAL;
	//				LOG("Missing phenotype sample `%s:%s`\n",
	//					X_sample.S_FID.c_str(), X_sample.S_IID.c_str());
				} else {
					/* Special for 0 */
					if (!strcmp(a, "0")) N_0++;

					/* Otherwise, it is continuous phenotype */
					Ba_typePheno[0] = 1;
				}
			}
			/* Filtering option */
// 			else if ((OPT_ENABLED(filcase) && Ra_pheno[0][j] == OPT_NUMBER(phenoCase)) ||
// 				(OPT_ENABLED(filcontrol) && Ra_pheno[0][j] == OPT_NUMBER(phenoCtrl)))
//				Ba_isSampFilt[j] = 1;

			/* Status */
			if (((++j)%100) == 0) notice("%d samples retrieved\r", j);
		}
		if (j != N_sample)
			halt("SYSERR : Retrieved sample [%d] is not same with expected [%d]", j, N_sample);
		/* Possible chance to miss --1
		 * 150305 do not print it when there is alternative phenotype */
		if (N_0 && !IS_ASSIGNED(pname)) {
			LOGwarn("This file could be coded 0/1 but not assumed that. Please check it!\n");
			LOG("        ( 0 found %d times )\n", N_0);
		}
		//fclose(H_famFile);
		LOG("%d samples found in FAM file (%d filtered)\n", N_sample, N_sampOrig-N_sample);
	} else if (B_vrtWise && B_header) {
		/* # of founders == # of samples */
		N_founder = N_sample;

		/* Back to the start of file */
		C_dosage.close();
		C_dosage.open(S_fn);

		/* --nskip */
		wsUint N_skip = OPT_NUMBER(nskip);
		if (IS_ASSIGNED(nskip)) for (wsUint i=0 ; i<N_skip ; i++)
			if (C_dosage.gets(S_buf, 1024*1024) == 0)
				halt_fmt(WISARD_NULL_PEDCONTENTS, S_fn);

		C_dosage.gets(S_buf, 1024*1024);

		/* Read header */
		/* Separate..., assume that the dosage file is variant-wise */
		if (Sp_delim)	getStringDelim(&S_buf, &a, *Sp_delim);
		else			getString(&S_buf, &a);

		/* Pass first field */
		if (Sp_delim) getStringDelim(&a, &b, *Sp_delim);
		else getString(&a, &b);

		/* Pass second and third if BEAGLE */
		if (X_type == DSG_BEAGLE) {
			getString(&b, &a);
			getString(&a, &b);
		}

		/* For all column headers */
		for (wsUint i=0,j=0 ; a ; i++) {
			/* Get record */
			/* Otherwise, insert */
			string S_currIID = (string)a;
			xSample	&X_sample = Xm_iid2data[S_currIID];

			/* If current samples should be filtered, do not retrieve anything */
			if (Bv_filtSamp[i] == 1) {
				X_sample.N_idx = SAMP_NODATA;
				continue;
			}

			/* Allocate initial data index */
			X_sample.N_idx = j;
			Xa_sampleV2.push_back(&X_sample);
			j++;

			/* Go to the next ptr */
			a = b;
			if (a) {
				if (Sp_delim) getStringDelim(&a, &b, *Sp_delim);
				else getString(&a, &b);
			}
		}
	}

	/* FIXME : If vrtWise == 0 */
	cExporter *Cp_misInds = NULL;
	wsUint diIdx = 0, fiIdx = 0;
	if (B_vrtWise == 0) for ( ; C_dosage.gets(S_buf, 1024*1024) ; fiIdx++) {
		a = S_buf;
		wsUint N_missGenoCurrSamp = 0;
		string S_currFID, S_currIID;

		/* Get FID & IID */
		if (X_type != DSG_MINIMAC)
			b = loadFIDIID(a, S_currFID, S_currIID);
		else {
			getString(&a, &b);
			// a ==> fid->iid
			char* ptr = strstr(a, "->");
			if (!ptr) halt("[%d] th line of minimac dosage does not have FID/IID information", fiIdx+1);
			*ptr = '\0';
			ptr += 2;
			// a = fid
			// ptr = iid
			S_currFID = a;
			S_currIID = ptr;
		}

		xSample* Xp_sample = &(Xm_iid2data[S_currIID]);
		if (X_type != DSG_MINIMAC) {
			/* Skip PAT & MAT w/o --noparent */
			if (!OPT_ENABLED(noparent)) {
				getString(&b, &a);
				getString(&a, &b);
			}

			/* Skip w/o --nosex */
			getString(&b, &a);

			/* Skip w/o --nopheno */
			char *Sp_shouldNull = NULL;
			wsVec Rp_pheno = Ra_pheno[0];
			getString(&a, &b);
			if (!OPT_ENABLED(nopheno)) {
				Rp_pheno[diIdx] = (wsReal)str2dbl(a, &Sp_shouldNull);

				/* If phenotype is missing */
				if (!strcmp(a, OPT_STRING(mispheno))) {
					Rp_pheno[diIdx] = WISARD_NA_REAL;
					LOG("Missing phenotype [%s] sample [%s:%s]\n", a,
						Xp_sample->S_FID.c_str(), Xp_sample->S_IID.c_str());
				} else if (Rp_pheno[diIdx] != OPT_NUMBER(phenoCtrl) &&
					Rp_pheno[diIdx] != OPT_NUMBER(phenoCase)) {
					if (Sp_shouldNull &&  Sp_shouldNull[0])
						halt("Phenotype [%s] for sample [%s] looks invalid",
							b, Xp_sample->S_FID.c_str());

					Ba_typePheno[0] = 1; /* Set to continuous */
				}
			} else
				Ra_pheno[0][diIdx] = WISARD_NA_REAL;
		} else {
			Xp_sample->S_FID = S_currFID;
			Xp_sample->S_IID = S_currIID;
			// Skip DOSE / ML_DOSE
			char *tmp = NULL;
			getString(&b, &tmp);
			b = tmp;
		}

		/* Is this sample founder? */
		if (B_isFAM) {
			Ba_isFounder[diIdx] = !Xp_sample->Xp_mat && !Xp_sample->Xp_pat;
			if (Ba_isFounder[diIdx]) N_founder++;
		} else N_founder++;

		/* Allocate initial data index */
		Xp_sample->N_idx = diIdx;
		Xa_sampleV2.push_back(Xp_sample);
		if (X_type == DSG_MINIMAC) {
			N_sampOrig++;
			N_sample++;
		}

		wsUint fsIdx = 0, dsIdx = 1;
		wsStrCst S_misGeno = IS_ASSIGNED(misgeno) ? OPT_STRING(misgeno) : NULL;

		/* Now read data */
		short*		Na_tmpDsg	= NULL;
		wsFloat*	Ra_tmpDsg	= NULL;
		if (IS_ASSIGNED(expression)) {
			wsAlloc(Ra_tmpDsg, wsFloat, N_variant);
		} else {
			wsAlloc(Na_tmpDsg, short, N_variant);
		}
		wsUint N_retrived = 0, N_data = 0;

		for (; fsIdx < N_vrtOrig; fsIdx++) {
			if (b == NULL)
				break;

			getString(&b, &a);

			/* Retrieve information only if it is not filtered */
//			if (Ba_isVrtFiltOrig[fsIdx] == 0) {

			// 150824 expression input
				if (IS_ASSIGNED(expression)) {
					char *Sp_ptr = NULL;
//					wsReal tmp =
						(wsReal)strtod(b, &Sp_ptr);//atoi(Sp_tmp1);
					if (Sp_ptr && Sp_ptr[0]) {
						/* Check is it NA */
						if (S_misGeno && !stricmp(Sp_ptr, S_misGeno)) {
//							tmp = -1;
						} else halt("%s file [%s] looks like RAW format, but non-numeric "
							"character [%s] found in line [%d]!",
							S_type, S_fn, b, fiIdx);
					}
				} else {
					char *Sp_ptr = NULL;
					double tmp = strtod(b, &Sp_ptr);//atoi(Sp_tmp1);
					if (Sp_ptr && Sp_ptr[0]) {
						/* Check is it NA */
						if (S_misGeno && !stricmp(Sp_ptr, S_misGeno)) {
							tmp = -1;
						} else halt("%s file [%s] looks like RAW format, but non-numeric "
							"character [%s] found in line [%d]!",
							S_type, S_fn, b, fiIdx);
					}

					/* Only 0/1/2 are acceptable, others are missing */
					if (tmp < 0.0 || tmp > 2.0 || Sp_ptr[0]) {
						Xa_sampleV2[diIdx]->B_isComplete = -1;
						Xv_variant[dsIdx -1].nmis++;
						tmp = WISARD_NA;
						N_missGeno++;
						N_missGenoCurrSamp++;
					}
				}

				/* Pass if filtered */
				if (Bv_filtVrt[fsIdx] == 0) {
					if (stricmp(Sp_misGeno, b)) {
						/* Convert and retrieve */
						wsReal	R_dsgVal = (wsReal)atof(b);
						if (IS_ASSIGNED(expression))
							Ra_tmpDsg[N_data] = (wsFloat)R_dsgVal;
						else {
							/* Check the range of dosage */
							if (R_minDsg > R_dsgVal) R_minDsg = R_dsgVal;
							if (R_maxDsg < R_dsgVal) R_maxDsg = R_dsgVal;

							Na_tmpDsg[N_data] = (short)(R_dsgVal *
								REAL_CONST(10000.0) + REAL_CONST(0.5));
						}
					} else {
						if (IS_ASSIGNED(expression))
							Ra_tmpDsg[N_data] = (wsFloat)WISARD_NA;
						else
							Na_tmpDsg[N_data] = (short)0xffff;

						N_missGeno++;
						Xv_variant[dsIdx-1].nmis++;
						Xa_sampleV2[diIdx]->N_missGeno++;
						Xa_sampleV2[diIdx]->B_isComplete = -1;
					}

					/* Increase data index */
					N_data++;
					dsIdx++;
				}

				N_retrived++;

				b = a;
//			}
		}
		if (b) {
			if (IS_ASSIGNED(dosage))
				halt("Input dosage file [%s] line [%d] have more data than"
					" defined! (expected [%d] variants)", OPT_STRING(dosage),
					fiIdx+2, N_vrtOrig);
			else if (IS_ASSIGNED(expression))
				halt("Input expression file [%s] line [%d] have more data than"
					" defined! (expected [%d] expressions)", OPT_STRING(expression),
					fiIdx+2, N_vrtOrig);
		}
		/* Check sanity of PED file */
		if ((dsIdx - 1) != N_variant || fsIdx != N_vrtOrig)
			halt_fmt(WISARD_INVL_NVRTOFSAMPLE, Xp_sample->S_FID.c_str(),
			Xp_sample->S_IID.c_str(), N_vrtOrig, N_variant, fsIdx, dsIdx - 1);

		if (N_vrtOrig) {
			wsReal R_gr = W1 - N_missGenoCurrSamp / (wsReal)N_variant;
			/* In default, someone having NO genotype will be excluded */
			/* Filtering by --geno when RAW file */
			if ((N_missGenoCurrSamp == (int)N_variant && !OPT_ENABLED(nodata)) ||
				(IS_ASSIGNED(filgind) &&
					isInRange(OPT_RANGE(filgind), R_gr)) ||
				(IS_ASSIGNED(incgind) &&
					!isInRange(OPT_RANGE(incgind), R_gr))) {
				if (N_missGenoCurrSamp == (int)N_variant && Cp_misInds) {
					if (Cp_misInds == NULL) {
						Cp_misInds = cExporter::summon("allmiss.sample.lst");
						LOGoutput("Sample with all missing genotypes found, it "
							"will exported to [%s.allmiss.sample.lst]\n",
							OPT_STRING(out));
					}
					Cp_misInds->fmt("%s	%s\n", S_currFID.c_str(),
						S_currIID.c_str());
				}
				Bv_filtSamp[diIdx] = 1;
				N_szFilteredSamp++;
			}
		}

		diIdx++;
		if (IS_ASSIGNED(expression))
			Xa_expr.push_back(Ra_tmpDsg);
		else
			Xa_dsg.push_back(Na_tmpDsg);
	} else {
		N_variant = 0;

		if (!B_header) C_dosage.rewind();

		for ( ; C_dosage.gets(S_buf, 1024*1024) ; N_fVrt++) {
			if (Sp_delim) getStringDelim(&S_buf, &a, *Sp_delim);
			else getString(&S_buf, &a);

			char* Sp_vrtName = NULL;
			switch (X_type) {
			case DSG_BEAGLE:	Sp_vrtName = S_buf; break;
			case DSG_IMPUTE2:	Sp_vrtName = a; break;
			default:
				halt("SYSERR: Uncontrolled dosage type [%d] in dosage data reading", X_type); break;
			}

			/* Skip one more if IMPUTE2 */
			if (X_type == DSG_IMPUTE2) {
				getString(&a, &b);
				a = b;
			}

			/* Check quality-related */
			xVariant	X_currVrt	= {0, };

			if (B_alleleInfo == 1) {
				getString(&a, &b);
				X_currVrt.al1 = a[0];
				getString(&b, &a);
				X_currVrt.al2 = b[0];
			} else {
				/* For convenience, major allele=1 and minor allele=0 */
				X_currVrt.al1		= '0';
				X_currVrt.al2		= '1';
			}

			char	B_filtered	= 0;
			/* S_tok == Variant name */
			X_currVrt.indel1	= NULL;
			X_currVrt.indel2	= NULL;
			X_currVrt.filter	= 0;
			/* Same for gdist */
			X_currVrt.gdist	= W0;
			/* Chromosome number is undeterministic, because RAW file does not contain such information */
			X_currVrt.chr		= 0;
			/* Same for position */
			X_currVrt.pos		= 0;
			X_currVrt.name		= strdup(Sp_vrtName);

			/* Position section if IMPUTE2 */
			if (X_type == DSG_IMPUTE2) {
				getString(&a, &b);
				X_currVrt.pos = atoi(a);
				a = b;
			}

			/* Check this variant's name if filtering is activated */
			if (N_filterVrt == FILT_REMOVE) {
				/* Remove IF match */
				if (mapIsExist(Xe_listVrt, X_currVrt.name)) B_filtered = 1;
			} else if (N_filterVrt == FILT_SELECT) {
				/* Remove IF NOT match */
				if (mapIsInexist(Xe_listVrt, X_currVrt.name)) B_filtered = 1;
			}

			/* Quality-related */
			if (M_qual.size()) {
				if ((IS_ASSIGNED(incqual) && !isInRange(OPT_RANGE(incqual), (wsReal)M_qual[X_currVrt.name])) ||
					(IS_ASSIGNED(filqual) && isInRange(OPT_RANGE(filqual), (wsReal)M_qual[X_currVrt.name])))
					B_filtered = 1;
			}

			if (B_filtered == 1)
				Bv_filtVrt.push_back(1);
			else {
				Bv_filtVrt.push_back(0);

				/* Now read data */
				short	*Na_vrt	= NULL;
				wsAlloc(Na_vrt, short, N_sample);
				wsUint N_retrived = 0, N_data = 0;

				/* Move to next variant */
				if (Sp_delim) getStringDelim(&a, &b, *Sp_delim);
				else getString(&a, &b);

				while (a) {
					/* Overflow check */
					if (N_retrived >= N_sampOrig)
						halt("Variant [%s] at line [%d] have more than [%d] samples, which does not match with "
							"the # of samples defined in header section", X_currVrt.name, N_fVrt+2, N_retrived);

					/* If IMPUTE2, get two fields more*/
					char *v2 = NULL, *v3 = NULL;
					if (X_type == DSG_IMPUTE2) {
						getString(&b, &v3);
						v2 = b;
						getString(&v3, &b);
						N_fVrt++;
					}

					/* Pass if filtered */
					if (Bv_filtSamp[N_retrived] == 0) {
						if (stricmp(Sp_misGeno, a)) {
							/* Convert and retrieve */
							wsReal	R_dsgVal	= W0;
							if (X_type == DSG_BEAGLE) {
								R_dsgVal	= (wsReal)atof(a);
							} else if (X_type == DSG_IMPUTE2) {
								wsReal	R_d1 = (wsReal)atof(v2);
								wsReal	R_d2 = (wsReal)atof(v3);
								R_dsgVal = R_d1 + R_d2*W2;
							}
							short N_conv = (short)(R_dsgVal *
								REAL_CONST(10000.0) + REAL_CONST(0.5));
							Na_vrt[N_data] = N_conv;

							/* Check the range of dosage */
							if (R_minDsg > R_dsgVal) R_minDsg = R_dsgVal;
							if (R_maxDsg < R_dsgVal) R_maxDsg = R_dsgVal;
						} else {
							Na_vrt[N_data] = (short)0xffff;

							N_missGeno++;
							X_currVrt.nmis++;
							Xa_sampleV2[N_data]->N_missGeno++;
							Xa_sampleV2[N_data]->B_isComplete = -1;
						}

						/* Increase data index */
						N_data++;
					}

					N_retrived++;
					/* Move to next */
					a = b;
					if (a) {
						if (Sp_delim) getStringDelim(&a, &b, *Sp_delim);
						else getString(&a, &b);
					}
				}

				/* Check */
				if (N_retrived != N_sample)
					halt("Number of samples[%d] for variant [%s] at [%d] line "
						"does not match with the # of samples defined in header "
						"section", N_retrived, X_currVrt.name, N_fVrt+2);
				S_buf[0] = '\0';
				N_variant++;

				/* Insert variant info */
				Xv_variant.push_back(X_currVrt);

				/* Insert dosage data */
				Xa_dsg.push_back(Na_vrt);
			}

			if ((N_variant%1000) == 0)
				notice("[%d] variants detected, [%d] variants retrieve...\r", N_fVrt, N_variant);

			/* Range check */
			if (R_minDsg < -3.2767 || R_maxDsg > 3.2767)
				halt("Invalid dosage value found in variant [%s] at line %d",
					X_currVrt.name, N_fVrt+2);
		}
	}
	LOG("[%d] samples and [%d] variants found in dosage data\n", N_sampOrig,
		N_fVrt);

	/* Allocate data */
	wsAlloc(Ra_dosage, wsFloat*, N_sample);
	for (wsUint i=0 ; i<N_sample ; i++) {
		sseMalloc(Ra_dosage[i], wsFloat, N_variant);
	}

	if (IS_ASSIGNED(expression)) {
		/* Convert data into real form */
		for (wsUint i=0 ; i<N_sample ; i++) {
			wsFloat*	Np_data = Xa_expr[i];
			for (wsUint j=0 ; j<N_variant ; j++)
 				Ra_dosage[i][j] = Np_data[j];
			if ((i%1000) == 0)
				notice("%d/%d variants converted...\r", i, N_variant);
		}
		FOREACH(vArrReal_it, Xa_expr, i) {
			DEALLOC(*i);
		}
	} else {
		/* Convert data into real form */
		if (!B_vrtWise) for (wsUint i=0 ; i<N_sample ; i++) {
			short *Np_data = Xa_dsg[i];
			for (wsUint j=0 ; j<N_variant ; j++) {
				if (Np_data[j] == (short)0xffff)
					Ra_dosage[i][j] = (wsFloat)WISARD_NA_REAL;
				else
					Ra_dosage[i][j] = Np_data[j] / CORR_CONST(10000.0);
			}
			if ((i%1000) == 0)
				notice("%d/%d samples converted...\r", i, N_sample);
		} else for (wsUint i = 0; i < N_variant; i++) {
			short *Np_data = Xa_dsg[i];
			for (wsUint j = 0; j < N_sample; j++) {
				if (Np_data[j] == (short)0xffff)
					Ra_dosage[j][i] = (wsFloat)WISARD_NA_REAL;
				else
					Ra_dosage[j][i] = Np_data[j] / CORR_CONST(10000.0);
			}
			if ((i % 1000) == 0)
				notice("%d/%d variants converted...\r", i, N_variant);
		}
		/* Dealloc */
		FOREACH(vArrWord_it, Xa_dsg, i) {
			DEALLOC(*i);
		}
	}
	LOG("[%d] variants converted and loaded into memory\n", N_variant);

	if (N_missGeno == 0)
		B_isDataComplete = 1;
	notice("%d/%d variants retreived...\r", N_fVrt, N_fVrt);
	LOG("Range of dosage value : %g ~ %g\n", R_minDsg, R_maxDsg);

	/* Export dosage distribution */
	if (B_dsgdist)
		_exportDosageDist(R_minDsg, R_maxDsg);


	DEALLOC(S_buf);
}

} // End namespace ONETOOL
