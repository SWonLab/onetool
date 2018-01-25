#include <math.h>
#include "input/mach.h"
#include "input/stream.h"

namespace ONETOOL {

cMachIO::cMachIO(char *S_fn, char B_inDry)
{
	B_dry = B_inDry;
	init(S_fn);
}

char cMachIO::_loadMain(char *Sp_buf, wsUint fiIdx, wsUint diIdx,
	wsReal *Rp_pheno, vBool& Ba_isSampFiltOrig, vBool& Ba_isVrtFiltOrig,
	cExporter **Cp_misInds)
{
	char	*Sp_tmp1 = NULL;
	char	*Sp_tmp2 = NULL;
	string	S_currFID, S_currIID;
	//		int		ns = 0;
	int		N_missGenoCurrSamp = 0;

	/* If current samples should be filtered, do not retrieve anything */
	if (Ba_isSampFiltOrig[fiIdx] == 1)
		return 1;

	/* Pass FID/IID/PAT/MAT/SEX */

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
	/* What is default genotype? */
	wsStrCst Sp_misGeno = IS_ASSIGNED(misgeno) ? OPT_STRING(misgeno) : "N";
// 	else if ((OPT_ENABLED(filcase) && Rp_pheno[diIdx] == OPT_NUMBER(phenoCase)) ||
// 		(OPT_ENABLED(filcontrol) && Rp_pheno[diIdx] == OPT_NUMBER(phenoCtrl)))
// 		Ba_isSampFilt[diIdx] = 1;

	char S_delim = IS_ASSIGNED(sepallele) ? OPT_STRING(sepallele)[0] : 0xff;
	Sp_tmp1 = Sp_tmp2;
	wsUint fsIdx = 0, dsIdx = 0;
	for ( ; fsIdx<N_vrtOrig ; fsIdx++) {
		xVariant &X_vrt = getVariant()[fsIdx];
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

		/* Sp_gen1 & Sp_gen2 makes genotype */
		bool N_m11, N_m12, N_m21, N_m22;
		bool N_mm1, N_mm2;
		if (OPT_ENABLED(indel)) {
			/* Check match */
			N_m11 = !stricmp(Sp_gen1, X_vrt.indel1);
			N_m12 = !stricmp(Sp_gen1, X_vrt.indel2);
			N_m21 = !stricmp(Sp_gen2, X_vrt.indel1);
			N_m22 = !stricmp(Sp_gen2, X_vrt.indel2);
			N_mm1 = !stricmp(Sp_gen1, Sp_misGeno);
			N_mm2 = !stricmp(Sp_gen2, Sp_misGeno);
		} else {
			/* Should not be multiallelic */
			if (Sp_gen1[1] || Sp_gen2[1])
				LOG("Multi-base allele [%s/%s] found in line [%d], use --indel",
					Sp_gen1, Sp_gen2, fiIdx);
			/* Check match */
			N_m11 = Sp_gen1[0] == X_vrt.al1;
			N_m12 = Sp_gen1[0] == X_vrt.al2;
			N_m21 = Sp_gen2[0] == X_vrt.al1;
			N_m22 = Sp_gen2[0] == X_vrt.al2;
			N_mm1 = Sp_gen1[0] == Sp_misGeno[0];
			N_mm2 = Sp_gen2[0] == Sp_misGeno[0];
		}
		/* nor m11 m12 mm1 */
		if (!(N_m11 | N_m12 | N_mm1))
			halt("Unrecognized genotype [%s] found in line [%d]", Sp_gen1, fiIdx);
		/* nor m21 m22 mm2 */
		if (!(N_m21 | N_m22 | N_mm2))
			halt("Unrecognized genotype [%s] found in line [%d]", Sp_gen2, fiIdx);
		/* mm1 but !mm2 or !mm1 but mm2 */
		if ((N_mm1 && !N_mm2) || (!N_mm1 && N_mm2))
			halt_fmt(WISARD_CANT_ONEALLELE_MISSING, X_vrt.name, fiIdx);
		/* Set allele */
		if (N_mm1) {
			N_missGenoCurrSamp++;
			Na_geno[diIdx][dsIdx] = WISARD_NA;
		} else
			Na_geno[diIdx][dsIdx] = N_m12 + N_m22;

		dsIdx++;
	}
	/* Check sanity of PED file */
	if (dsIdx != N_variant || fsIdx != N_vrtOrig)
		halt_fmt(WISARD_INVL_NVRTOFSAMPLE, X_sample.S_FID.c_str(),
			X_sample.S_IID.c_str(), N_vrtOrig, N_variant, fsIdx, dsIdx);

	if (N_vrtOrig) {
		wsReal R_gr = W1 - N_missGenoCurrSamp / (wsReal)N_variant;
		/* In default, someone having NO genotype will be excluded */
		/* Filtering by --geno when RAW file */
		if (N_missGenoCurrSamp == (int)N_variant ||
			(IS_ASSIGNED(filgind) &&
				isInRange(OPT_RANGE(filgind), R_gr)) ||
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

void cMachIO::init(char *S_fn)
{
	if (S_fn == NULL)
		halt("SYSERR: NULL file name given");

	/* Initialize */
	char *Sp_buf = NULL;
	wsAlloc(Sp_buf, char, PED_MAXLEN);

	/* Get .mlinfo */
	char S_mlInfo[512];
	otherExt(S_fn, "mlinfo", S_mlInfo);
	cStrFile C_info(S_mlInfo, "Mach variant info file");
	LOG("Reading variant information from [%s]...\n", S_mlInfo);
	char *Ba_isSelResize	= NULL;
	if (IS_ASSIGNED(varresize)) {
		wsReal R_valSS	= OPT_REAL(varresize);
		wsUint N_sel = R_valSS >= W1 ?
			(wsUint)R_valSS : (wsUint)round(R_valSS*N_vrtOrig);

		/* Ignore the first line : Header */
		C_info.gets(Sp_buf, 1024*1204);

		/* Retrieve the number of variants in BIM file first */
		for (N_vrtOrig=0 ; C_info.gets(Sp_buf, 1024*1024) ; N_vrtOrig++);
		C_info.rewind();

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
	/*
		--autoonly
		--sexonly
		--chr
		--nofid
		--noparent
		--ignoreparent
		--ignorefid
		--filmale
		--filfemale
		--filnosex
	*/
	/* Ignore the first line : Header */
	C_info.gets(Sp_buf, 1024*1204);
	/* VARIANT     Al1     Al2     Freq1   MAF     Quality Rsq */
	char *a, *b;
	N_variant = 0;
	for (N_vrtOrig=0 ; C_info.gets(Sp_buf, 1024*1024) ; N_vrtOrig++) {
		a = Sp_buf;
		/* Do not set chromosome */
		xVariant X_ent = {0, };
		char B_filter = 0;

		getString(&a, &b);
		// a = name
		wsAlloc(X_ent.name, char, strlen(a)+1);
		strcpy(X_ent.name, a);

		/* Check --varresize */
		if (Ba_isSelResize && Ba_isSelResize[N_vrtOrig]==0)
			B_filter = 1; /* Filter if not selected */

		/* Check variant filter */
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

				/* --acgt */
				if (Na_ACGT) {
					for (char *x=X_ent.indel1 ; *x ; x++) {
						wsUint X = *x;
						if (!Na_ACGT[X]) halt("Character [%c] in variant [%s] is not in --acgt",
							X_ent.name, *x);
						*x = Na_ACGT[X];
					}
					if (X_ent.indel2) for (char *x=X_ent.indel2 ; *x ; x++) {
						wsUint X = *x;
						if (!Na_ACGT[X]) halt("Character [%c] in variant [%s] is not in --acgt",
							X_ent.name, *x);
						*x = Na_ACGT[X];
					}
				}

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
				X_ent.al1	= a[0];
				X_ent.al2	= b[0];

				/* --acgt */
				if (Na_ACGT) {
					if (!Na_ACGT[(wsUint)a[0]]) halt("Character [%c] in variant [%s] is not in --acgt",
						X_ent.name, a[0]);
					if (!Na_ACGT[(wsUint)a[1]]) halt("Character [%c] in variant [%s] is not in --acgt",
						X_ent.name, a[1]);
					X_ent.al1 = Na_ACGT[(wsUint)a[0]];
					X_ent.al2 = Na_ACGT[(wsUint)a[1]];
				}
			}
		} else {
			if (Na_ACGT) halt("Variant [%s] does not have allele info, --acgt makes halt", X_ent.name);

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

		/* Get freq for allele 1 */
		getString(&b, &a);
		wsReal R_freq1 = atof(b);
		if (R_freq1 < REAL_CONST(0.5)) {
			/* Switch allele */
			if (X_ent.indel1) {
				char *Sp_tmp = X_ent.indel1;
				X_ent.indel1 = X_ent.indel2;
				X_ent.indel2 = Sp_tmp;
			} else {
				char S_tmp = X_ent.al1;
				X_ent.al1 = X_ent.al2;
				X_ent.al2 = S_tmp;
			}
		}

		/* Skip MAF */
		getString(&a, &b);

		/* Get quality */
		getString(&b, &a);
		wsReal R_qual = atof(b);

		/* --incqual, --filqual */
		if ((IS_ASSIGNED(incqual) && !isInRange(OPT_RANGE(incqual), R_qual)) ||
			(IS_ASSIGNED(filqual) && isInRange(OPT_RANGE(filqual), R_qual)))
			B_filter |= 1;

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
	}
	LOG("%d variants found in BIM file, %d variants were chosen\n", N_vrtOrig, N_variant);
	DEALLOC(Ba_isSelResize);

	/* Read alternative FAM file if exists */
	N_sample = 0;
	wsUint N_fSamp = 0;
	if (IS_ASSIGNED(fam)) {
		cStrFile C_fam(OPT_STRING(fam), "Alternative FAM file");
		char S_buf[1024];

		/* Construct sample data */
		for (N_fSamp=0 ; C_fam.gets(S_buf, 1024) ; N_fSamp++)
			_loadRealTimeFamData(S_buf, N_fSamp);
	} else {
		/* Option checking related with sex */
		LOGwarn("Alternative FAM file is not assigned, sex-related option will be ignored\n");
		if (OPT_ENABLED(filnosex))	LOG("    --nosex ignored\n");
		if (OPT_ENABLED(filmale))	LOG("    --filmale ignored\n");
		if (OPT_ENABLED(filfemale))	LOG("    --filfemale ignored\n");

		/* For each line */
		cStrFile C_mgeno(S_fn, "Mach genotype file");
		for (wsUint L=1 ; C_mgeno.gets(Sp_buf, PED_MAXLEN) ; L++) {
			char *a = NULL;
			getString(&Sp_buf, &a);

			/* Get FID and IID */
			char *Sp_sub = strstr(Sp_buf, "->");
			if (Sp_sub == NULL)
				halt("Invalid FID/IID in Mach genotype file [%s] in line [%d]",
					S_fn, L);
			Sp_sub[0] = '\0';
			Sp_sub += 2;
_domore:
			string S_currFID = (string)Sp_buf;
			string S_currIID = (string)Sp_sub;

			/* Check this sample should be filtered by --samprem/--sampsel */
			char B_filtered = N_filterSample && isSampleFiltered(S_currFID,
				S_currIID, ((wsUint)Bv_filtSamp.size()-N_sample)>5);

			/* Check three samples in this line (sample itself, father and mother) are previously inserted */
			mSamp_it	X_resFind		= Xm_iid2data.find(S_currIID);
			mSamp_it	X_fail			= Xm_iid2data.end();
			xSample *Xp_thisSample = NULL; ///< Pointer of curr line's xSample structure
			xSample X_newMyData;
			string	S_altIID;

			if (X_resFind != X_fail) {
				/* Someone previously inserted this sample itself */
				xSample &X_myData = X_resFind->second;

				/* HALT if isAssigned==1, it means my record previously existed */
				if (X_myData.B_isAssigned == 1) {
					if (!OPT_ENABLED(dupnaming))
						halt("Sample [%s:%s] have duplicated record in line [%d]",
							S_currFID.c_str(), S_currIID.c_str(), L);
					else {
						for (wsUint NN=1 ; ; NN++) {
							char S[32];
							sprintf(S, "%d", NN);
							S_altIID = S_currIID + "_" + S;
							mSamp_it	X_resFind		= Xm_iid2data.find(S_altIID);
							if (X_resFind == Xm_iid2data.end() || X_resFind->second.B_isAssigned == 0)
								break;
						}
						/* Rewrite */
						goto _domore;
					}
				}

				/* Assign */
				X_myData.B_isAssigned	= 1;
				X_myData.N_oriIdx		= L-1;
				X_myData.B_isMissing	= B_filtered;
				X_myData.Xp_mat			= NULL;
				X_myData.Xp_pat			= NULL;
				X_myData.B_isProband	= 0;
				X_myData.B_isComplete	= 0;
				X_myData.N_sex			= 0;
				X_myData.N_grpFst		= -1;

				Xp_thisSample = &X_myData;
			} else {
				/* This sample is currently new in xSample map */

				/* Fill up with basic information */
				X_newMyData.S_FID			= S_currFID;
				X_newMyData.S_IID			= S_currIID;
				X_newMyData.N_sex			= 0;
				X_newMyData.Xp_mat			= NULL;
				X_newMyData.Xp_pat			= NULL;
				X_newMyData.B_isAssigned	= 1;
				X_newMyData.B_isMissing		= B_filtered;
				X_newMyData.N_idx			= SAMP_NODATA;
				X_newMyData.N_oriIdx		= L-1;
				X_newMyData.N_missGeno		= 0;
				X_newMyData.N_idTwin		= 0;
				X_newMyData.B_isProband		= 0;
				X_newMyData.B_isComplete	= 0;
				X_newMyData.N_grpFst		= -1;

				Xp_thisSample = &X_newMyData;
			}
			/* Sample filtering */
			if (B_filtered) {
				Bv_filtSamp.push_back(1);
			} else {
				N_sample++;
				Bv_filtSamp.push_back(0);
			}

			/* Insert current sample */
			Xp_thisSample->N_isVisited	= 0;
			pair<mSamp_it, bool>
						X_insRes		= Xm_iid2data.insert(make_pair(S_currIID, *Xp_thisSample));
			mSamp_it	X_thisSample	= X_insRes.first;

			/* If one of this sample's mother/father does not inserted to xSample map,
				* do insertion
				*/
			/* If founder, insert to family structure */
			mFam_it X_resFamFind = Xm_fid2data.find(S_currFID);

			/* This sample is the founder appeared at first */
			if (X_resFamFind == Xm_fid2data.end()) {
				xFamily X_newFam;

				X_newFam.N_size = 1;
				X_newFam.S_FID = S_currFID;
				X_newFam.Xp_founders.push_back(&X_thisSample->second);

				Xm_fid2data.insert(make_pair(S_currFID, X_newFam));
			} else {
				X_resFamFind->second.N_size++;
				X_resFamFind->second.Xp_founders.push_back(&X_thisSample->second);
			}
		}
	}

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

	/* For each line */
	cExporter *Cp_misInds = NULL;
	cStrFile C_mgeno(S_fn, "Mach genotype file");
	for (wsUint L=1,LL=0 ; C_mgeno.gets(Sp_buf, PED_MAXLEN) ; L++) {
		char *a = NULL, *b = NULL;
		getString(&Sp_buf, &a);

		/* Get FID and IID */
		if (!IS_ASSIGNED(fam)) {
			char *Sp_sub = strstr(Sp_buf, "->");
			if (Sp_sub == NULL)
				halt("Invalid FID/IID in Mach genotype file [%s] in line [%d]",
					S_fn, L);
		}

		getString(&a, &b);
		/* Check whether it is ML_GENO or not */
		if (stricmp("ML_GENO", a))
			halt("Invalid ML_GENO signature [%s] in Mach genotype file [%s] in line [%d]",
				a, S_fn, L);

		/* FIXME : Do filtering */

		/* Now retrieve genotypes by sample-wise */
		_loadMain(b, L-1, LL, NULL, Bv_filtSampOrig, Bv_filtVrtOrig, &Cp_misInds);
	}
}

} // End namespace ONETOOL
