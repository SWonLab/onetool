#include "global/common.h"
#include "global/option.h"
#include "global/io.h"
#include "utils/matrix.h"
#include "input/sampvar.h"
#include "input/stream.h"
#include <math.h>

#define IS_REG_COVAR(it) ((it)->R_times == W1 && (it)->N_idxMul == -1)

namespace ONETOOL {

cStream** _loadSampVar_check(cIO *Cp_IO, vPheno &Xa_asgPhenos, vCovar &Xa_asgCovars,
	mDataIdx &Xm_fileColNames, vStr &Xa_fColNames, vInt& Xv_resvIdx)//wsUint *Np_idxTw,
	//wsUint *Np_idxPb, wsUint *Np_idxFst)
{
	wsUint	N_files		= 0;
	char	**Sa_files	= loadStringValues(OPT_STRING(pheno), &N_files);

	/* Check */
	if (N_files == 0)
		halt("SYSERR : There is no file to load for sample variable!");

	/* Reserved columns */
	vStr	Xv_resvCol;
	/* Insert sequence MUST equal to the sequence of xReservedCol */
	Xv_resvCol.push_back(IS_ASSIGNED(probandcol) ? OPT_STRING(probandcol) : "PROBAND");
	Xv_resvCol.push_back(IS_ASSIGNED(twincol) ? OPT_STRING(twincol) : "TWIN");
	if (IS_ASSIGNED(fst) && IS_ASSIGNED(popuniq))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--fst", "--popuniq");
	Xv_resvCol.push_back(IS_ASSIGNED(fst) ? OPT_STRING(fst) :
		(IS_ASSIGNED(popuniq) ? OPT_STRING(popuniq) : "POP_GROUP"));
	Xv_resvCol.push_back(IS_ASSIGNED(sampleweight) && Cp_IO->getSampleWeight() ?
		OPT_STRING(sampleweight) : "__XXX__");

// 	char	*Sp_colPb	= NULL;
// 	char	*Sp_colTw	= NULL;
// 	char	*Sp_colFt	= NULL;
	char	*S_phenos	= IS_ASSIGNED(pname) ? strdup(OPT_STRING(pname)) : NULL;
	char	*S_covs		= IS_ASSIGNED(cname) ? strdup(OPT_STRING(cname)) : NULL;
	char	*S_fcovs	= IS_ASSIGNED(fname) ? strdup(OPT_STRING(fname)) : NULL;

	/*
	*
	* Initialization of reserved names
	* 
	*/

	/* Set to 0 */
	Xv_resvIdx.push_back(0);
	Xv_resvIdx.push_back(0);
	Xv_resvIdx.push_back(0);
	Xv_resvIdx.push_back(0);
	wsStrCst Sa_resvColID[] = {
		"Proband", "Twin", "Population group", "Sample-wise weight",
	};
// 	*Np_idxTw = 0;
// 	*Np_idxPb = 0;
// 	*Np_idxFst = 0;

	/* Allocate temporal memory */
// 	MULTI_MALLOC(Sp_colPb, char, 512);
// 	MULTI_MALLOC(Sp_colTw, char, 512);
// 	MULTI_MALLOC(Sp_colFt, char ,512);
 
	/* Substitute column name for PROBAND information if assigned */
// 	if (IS_ASSIGNED(probandcol)) {
// 		if (strlen(OPT_STRING(probandcol)) > 63)
// 			halt("Too long name for PROBAND column, should be less than 63 characters");
// 		strcpy(Sp_colPb, OPT_STRING(probandcol));
// 	} else
// 		strcpy(Sp_colPb, "PROBAND");

	/* Substitute column name for TWIN information if assigned */
// 	if (IS_ASSIGNED(twincol)) {
// 		if (strlen(OPT_STRING(twincol)) > 63)
// 			halt("Too long name for TWIN column, should be less than 63 characters");
// 		strcpy(Sp_colTw, OPT_STRING(twincol));
// 	} else
// 		strcpy(Sp_colTw, "TWIN");

	/* Set group information */
//	char *Sp_group = NULL;
// 	if (IS_ASSIGNED(fst) && IS_ASSIGNED(popuniq))
// 		halt_fmt(WISARD_CANT_EXCL_OPT, "--fst", "--popuniq");
// 	else if (IS_ASSIGNED(fst))
// 		Sp_group = OPT_STRING(fst);
// 	else if (IS_ASSIGNED(popuniq))
// 		Sp_group = OPT_STRING(popuniq);

	/* Substitute column name for Fst group information if assigned */
// 	if (Sp_group) {
// 		if (strlen(Sp_group) > 63)
// 			halt("Too long name for group column for population group, should be less than 63 characters");
// 		strcpy(Sp_colFt, Sp_group);
// 	} else
// 		strcpy(Sp_colFt, "POP_GROUP");

	/* Split --pname and --cname if unless they are exist, and not *(all) */
	if (S_phenos && strcmp("*", S_phenos)) {
		wsUint	N_phenos	= 0;
		char	**Sa_phenos	= loadStringValues(S_phenos, &N_phenos);
		for (wsUint i=0 ; i<N_phenos ; i++) {
			xPheno X_pheno = { Sa_phenos[i] };
			Xa_asgPhenos.push_back(X_pheno);
			free(Sa_phenos[i]);
		}
		DEALLOC(Sa_phenos);
	}
	if (S_covs && strcmp("*", S_covs)) {
		wsUint	N_covs		= 0;
		char	**Sa_covs	= loadStringValues(S_covs, &N_covs);
		for (wsUint i=0 ; i<N_covs ; i++) {
			/* Insert as covariates */
			xCovar X_covar = { WISARD_VAR_UNDET, 0, strdup(Sa_covs[i]), (int)0xffffffff,
				(int)0xffffffff, NULL, W1, -1, NULL };
			Xa_asgCovars.push_back(X_covar);
			free(Sa_covs[i]);
		}
		DEALLOC(Sa_covs);		
	}
	if (S_fcovs && strcmp("*", S_fcovs)) {
		wsUint	N_fcovs		= 0;
		char	**Sa_covs	= loadStringValues(S_fcovs, &N_fcovs);
		for (wsUint i=0 ; i<N_fcovs ; i++) {
			/* Insert as covariates */
			xCovar X_covar = { WISARD_VAR_FACTOR, 0, strdup(Sa_covs[i]), (int)0xffffffff,
				(int)0xffffffff, NULL, W1, -1, NULL };
			Xa_asgCovars.push_back(X_covar);
		}
	}

	/* Read sampvar flags if have */
	const char*	S_fid = "FID";
	const char*	S_iid = "IID";
	mStr	Xm_flag;
	if (IS_ASSIGNED(sampvarflag)) {
		loadFlags(OPT_STRING(sampvarflag), Xm_flag);

		if (mapIsExist(Xm_flag, "nofid") && mapIsExist(Xm_flag, "fid"))
			LOGwarn("Flag 'nofid' ignores flag 'fid'\n");

		if (mapIsExist(Xm_flag, "nofid")) S_fid = NULL;
		else if (mapIsExist(Xm_flag, "fid")) S_fid = (char *)Xm_flag["fid"].c_str();

		if (mapIsExist(Xm_flag, "iid")) S_iid = (char *)Xm_flag["iid"].c_str();
	}


	/* Initialize file stream */
	wsUint		N_fileCol	= 0;
	char*		Sp_buf		= NULL;
	cStream**	Cp_pheno	= NULL;
	wsAlloc(Cp_pheno, cStream*, N_files);
	wsAlloc(Sp_buf, char, 8192);
	for (wsUint X=0 ; X<N_files ; X++) {
		char		S_delim		= ' ';
		char		*a, *b;
		char		*Sp_curFn	= Sa_files[X];
		/* Open file */
		Cp_pheno[X]	= new cStrFile(Sp_curFn, "Sample-wise variable file");

		/* Prepare */
		wsUint	N_curFC		= 0;
		wsUint	N_idxVar	= X*10000;

		/* Read first line of input file : Should be header */
		if (Cp_pheno[X]->gets(Sp_buf, 8192) == 0)
			halt_fmt(WISARD_INVL_FILE_GENERAL, "Sample variable file",
				"Failed to read the first line, maybe file is empty?");

		/* Identify separator */ if (0) {
			wsUint N_comma=0, N_tab=0, N_space=0;
			for (char* c=Sp_buf ; *c ; c++) {
				if (*c == ',') N_comma++;
				else if (*c == '\t') N_tab++;
				else if (*c == ' ') N_space++;
			}
			if (N_comma >= N_tab && N_space >= N_tab){}
		}

		if (OPT_ENABLED(nophenohdr)) {
			LOG("--nophenohdr found, assume there is no header...\n");

			getStringDelim(&Sp_buf, &b, S_delim);
			if (b == NULL) {
				S_delim = ','; // 151228 try comma-separated
				getStringDelim(&Sp_buf, &b, S_delim);
				if (b == NULL)
					halt_fmt(WISARD_INVL_FILE_GENERAL, "Sample variable file",
						Sp_curFn, "At least three columns are required at the first row");
			}
			getStringDelim(&b, &a, S_delim);
			if (a == NULL)
				halt_fmt(WISARD_INVL_FILE_GENERAL, "Sample variable file",
					Sp_curFn, "No sample variable found in this file");
	// 		if (stricmp(Sp_buf, "FID") || stricmp(b, "IID"))
	// 			halt("First two columns should be 'FID' and 'IID'");

			/* Read which columns are exist in the file */
			for ( ; ; a=b) {
				getStringDelim(&a, &b, S_delim);

				/* Make arbitrary column name */
				char FX[64];
				sprintf(FX, "V%d", N_fileCol+N_curFC+1);

				/* Insert */
				Xm_fileColNames.insert(make_pair(FX, N_idxVar));
				Xa_fColNames.push_back(FX);
				N_curFC++;
				N_idxVar++;

				if (b == NULL) break;
			}
			Cp_pheno[X]->rewind();

			if (N_curFC) {
				if (N_curFC == 1)
					LOG("In sample variable file [%s], [V%d] found\n",
						Sa_files[X], N_fileCol+1);
				else
					LOG("In sample variable file [%s], [V%d] to [V%d] found\n",
						Sa_files[X],
						N_fileCol+1, N_curFC+N_fileCol);
			}
			N_fileCol += N_curFC;
		} else {
			/* Format check ; first two columns should be 'FID' and 'IID' */
			getStringDelim(&Sp_buf, &b, S_delim);
			if (b == NULL) {
				S_delim = ','; // 151228 try comma-separated
				getStringDelim(&Sp_buf, &b, S_delim);
				if (b == NULL)
					halt_fmt(WISARD_INVL_FILE_GENERAL, "Sample variable file",
						Sp_curFn, "At least three columns are required at the first row");
			}
			if (mapIsExist(Xm_flag, "nofid")) {
				a = b;
				b = Sp_buf;
			} else
				getStringDelim(&b, &a, S_delim);

			if (a == NULL)
				halt_fmt(WISARD_INVL_FILE_GENERAL, "Sample variable file",
				Sp_curFn, "No sample variable found in this file");
			/* Header name check */
			if ((S_fid && stricmp(Sp_buf, S_fid)) ||
				stricmp(b, S_iid))
				halt_fmt(WISARD_INVL_FILE_GENERAL, "Sample variable file",
				Sp_curFn, "No column(s) found for Family ID and Individual ID");

			/* Read which columns are exist in the file */
			for ( ; ; a=b,N_curFC++,N_idxVar++) {
				getStringDelim(&a, &b, S_delim);

				/* Duplication check */
				string		S_fx	= a;
				mDataIdx_it	X_srch	= Xm_fileColNames.find(S_fx);
				if (X_srch != Xm_fileColNames.end())

					halt("Column [%s] is duplicated: Existing in [%s]",
						a, Sa_files[X]);

				Xm_fileColNames.insert(make_pair(a, N_idxVar));
				Xa_fColNames.push_back(a);

				if (b == NULL) break;
			}
			N_fileCol += N_curFC;
		}

		LOG("Sample variable file [%s] loaded\n", Sp_curFn);
	}

	LOOP (i, N_files)
		free(Sa_files[i]);
	DEALLOC(Sa_files);

	/* Check reserved fields exist */ {
		wsUint I = 0;
		FOREACHDO (vStr_it, Xv_resvCol, i, I++)
			FOREACH (mDataIdx_it, Xm_fileColNames, it)
				if (it->first.compare(*i) == 0) {
					Xv_resvIdx[I] = it->second + 1;
					LOG("%s column '%s' found in the alternative phenotype file's"
						" %d th column\n", Sa_resvColID[I], i->c_str(), Xv_resvIdx[I]);
					Cp_IO->setReservedCol((xReservedCol)I);
					break;
				}
	}

	/* Check 'TWIN' exists */
// 	wsUint N_idxTwin = 0;
// 	FOREACH (mDataIdx_it, Xm_fileColNames, it)
// 		if (it->first.compare(Sp_colTw) == 0) {
// 			N_idxTwin = it->second + 1;
// 			LOG("Twin column '%s' found in the alternative phenotype file's"
// 				" %d th column\n", Sp_colTw, N_idxTwin);
// 			break;
// 		}
// 	if ((N_idxTwin%10000) > Xm_fileColNames.size()) N_idxTwin = 0;
// 	*Np_idxTw = N_idxTwin;

	/* Check 'FST_GROUP' exists */
// 	wsUint N_idxFst = 0;
// 	FOREACH (mDataIdx_it, Xm_fileColNames, it)
// 		if (it->first.compare(Sp_colFt) == 0) {
// 			N_idxFst = it->second + 1;
// 			LOG("Sample group column [%s] found in the alternative phenotype file's"
// 				" %d th column\n", Sp_colFt, N_idxFst);
// 			break;
// 		}
// 	if ((N_idxFst%10000) > Xm_fileColNames.size()) N_idxFst = 0;
// 	*Np_idxFst = N_idxFst;

	/* Check 'PROBAND' exists */
// 	wsUint N_idxProband = 0;
// 	FOREACH (mDataIdx_it, Xm_fileColNames, it)
// 		if (it->first.compare(Sp_colPb) == 0) {
// 			N_idxProband = it->second + 1;
// 			LOG("Proband column '%s' found in the alternative phenotype file's"
// 				" %d th column\n", Sp_colPb, N_idxProband);
// 			B_haveProband = 1;
// 			break;
// 		}
// 	if ((N_idxProband%10000) > Xm_fileColNames.size()) N_idxProband = 0;
// 	*Np_idxPb = N_idxProband;

 	N_fileCol++;

	/* Process --baseline option if assigned */
	if (IS_ASSIGNED(baseline)) {
		char *S_bl	= strdup(OPT_STRING(baseline));
		char *Sp_t	= strtok(S_bl, ",");

		do {
			/* stop if Sp_t == NULL */
			if (Sp_t == NULL) break;

			/* Divide again by = */
			char *Sp_div = strchr(Sp_t, '=');
			if (Sp_div == NULL)
				halt("--baseline error, should formatted as column=baseline_value");
			*(Sp_div++) = '\0';

			/* Check Sp_t have valid column names */
			vCovar_it it = Xa_asgCovars.begin();
			for ( ; it!=Xa_asgCovars.end() ; it++) {
				/* Match found, then regist baseline value */
				if (strcmp(it->Sp_varName, Sp_t) == 0) {
					wsAlloc(it->Sp_bl, char, strlen(Sp_div)+1);
					strcpy(it->Sp_bl, Sp_div);
					break;
				}
			}
			/* Not found in assigned covariates */
			if (it == Xa_asgCovars.end())
				halt("--baseline error, column [%s] is not assigned in --cname",
					Sp_t);
		} while ((Sp_t=strtok(NULL, ",")) != NULL);
	}


	DEALLOC(Sp_buf);
// 	DEALLOC(Sp_colPb);
// 	DEALLOC(Sp_colTw);
// 	DEALLOC(Sp_colFt);
	DEALLOC(S_phenos);
	DEALLOC(S_covs);

	return Cp_pheno;
}

void loadSampVar(cIO* Cp_IO)
{
	char	*a, *b;
	char	*S_phenos	= IS_ASSIGNED(pname) ? strdup(OPT_STRING(pname)) : NULL;
	char	*S_covs		= IS_ASSIGNED(cname) ? strdup(OPT_STRING(cname)) : NULL;
	//char	*S_fcovs	= IS_ASSIGNED(fname) ? strdup(OPT_STRING(fname)) : NULL;

	vStr	Xa_phenoNames;

	wsUint	N_sample	= Cp_IO->sizeSample();
	wsUint	N_newPheno;
	wsUint	N_fileCol = 0;
	mSamp&	Xm_sample	= Cp_IO->getSampleData();
	/* An array of indicators that ith column IN FILE will be the part of
	* phenotype or covariates */
	char	*Ba_isIncPheno = NULL, *Ba_isIncCov = NULL;

	char	*Sp_buf		= NULL;
	wsAlloc(Sp_buf, char, 8192);

	/*
	* Processing options
	*
	*/
	/*
	* Process first row
	* 
	*/
	vPheno		Xa_asgPhenos;	///< Column names of phenotypes that assigned from option
	vCovar		Xa_asgCovars;	///< Column names of covariates that assigned from option
	mDataIdx	Xm_fileColNames;	///< Key:Var names, value:10000(file idx),%10000(inner-file idx)
	vStr		Xa_fColNames;
	vInt		Xv_resvIdx;
	vPheno&		V_phenos	= Cp_IO->getPhenoInfo();
//	wsUint		N_idxTwin, N_idxProband, N_idxFst; ///< 1-based, 0 if not exist
	wsUint		N_files		= 0;
	char**		Sa_files	= loadStringValues(OPT_STRING(pheno), &N_files);
	cStream**	Cp_pheno	= _loadSampVar_check(Cp_IO, Xa_asgPhenos, Xa_asgCovars,
		Xm_fileColNames, Xa_fColNames, Xv_resvIdx);//&N_idxTwin, &N_idxProband, &N_idxFst);
	N_fileCol = (wsUint)Xm_fileColNames.size();
	LOG("[%d] columns were found from sample variable file(s)\n", N_fileCol);

	/* Allocate buffer that indicating inclusion for each column */
	wsCalloc(Ba_isIncPheno, char, N_fileCol);
	wsCalloc(Ba_isIncCov, char, N_fileCol);

	if (S_phenos && !strcmp("*", S_phenos)) {
		wsUint N_colSkip = 0, N_stop = 0;
		/* Assign 1 to N_fileCol sequentially to Ba_isIncPheno */
		for (wsUint i=0 ; i<N_fileCol ; i++,N_colSkip+=N_stop) {
			/* Skip */
			N_stop = 0;
			FOREACH (vInt_it, Xv_resvIdx, j)
				if ((*j - 1) == (int)i) {
					N_stop = 1;
					break;
				}
			if (N_stop) break;
// 			if ((N_idxProband-1) == i || (N_idxTwin-1) == i || (N_idxFst-1) == i)
// 				continue;

			Ba_isIncPheno[i] = (char)i+1;
			xPheno X;
			X.S_name = Xa_fColNames[i];
			V_phenos.push_back(X);
			Xa_phenoNames.push_back(Xa_fColNames[i]);
		}
		//memset(Ba_isIncPheno, 1, sizeof(char)*N_fileCol);
		N_newPheno = N_fileCol - N_colSkip;//(N_idxTwin != 0) - (N_idxProband != 0) - (N_idxFst != 0);
	/* Check the all phenotype names in --pname/--cname exists in alternative phenotype file */
	} else {
		N_newPheno = 0;

		/* For each desired phenotypes, check index and DO insert */
		wsUint j=1;
		FOREACH (vPheno_it, Xa_asgPhenos, it) {
			wsUint i=0;

			FOREACHDO (vStr_it, Xa_fColNames, iit, i++) {
				/* Pass if not match */
				if (it->S_name.compare(*iit)) continue;
				if (Ba_isIncPheno[i])
					halt("Column [%s] in [%s] assigned twice in --pname!",
						it->S_name.c_str(), it->S_fn.c_str());

				Ba_isIncPheno[i] = (char)j;
				//Sa_phenoName.push_back(*it);
				Xa_phenoNames.push_back(it->S_name);
				N_newPheno++;
				j++;
				break;
			}

			/* HALT if it is NOT exists on the altpheno */
			if (i == Xm_fileColNames.size()) {
				/* Check -(hyphen) character exists */
				char *S_dup	= strdup(it->S_name.c_str());
				char *S_isH	= strchr(S_dup, '-');
				if (S_isH) {
					/* Split them into two, and check both exists */
					*(S_isH++) = '\0';
					/* Check NULLity */
					if (!S_dup[0] || !S_isH[0])
						halt("Phenotype [%s-%s] was not found on sample "
							"variable file", S_dup, S_isH);
					int N_idx1 = -1;
					int N_idx2 = -1;
					wsUint K = 0;
					FOREACHDO (vStr_it, Xa_fColNames, iit, K++) {
						if (iit->compare(S_dup) == 0) N_idx1 = K;
						if (iit->compare(S_isH) == 0) N_idx2 = K;
					}
					/* Check both found */
					if (N_idx1 == -1)
						halt("Phenotype [%s-%s] recognized ranged type, but first one is does not exists",
							S_dup, S_isH);
					if (N_idx2 == -1)
						halt("Phenotype [%s-%s] recognized ranged type, but second one is does not exists",
							S_dup, S_isH);
					/* Check both are in same file */
					int N_idxS	= Xm_fileColNames[Xa_fColNames[N_idx1]];
					int N_idxE	= Xm_fileColNames[Xa_fColNames[N_idx2]];
					wsUint N_fS	= (wsUint)(N_idxS/10000.0);
					wsUint N_fE	= (wsUint)(N_idxE/10000.0);
					if (N_fS != N_fE)
						halt("Phenotype [%s] in [%s] but [%s] in [%s], "
							"both file must be match!", Sa_files[N_fS],
							Sa_files[N_fE]);

					/* Insert sequentially */
					wsUint N_incre = N_idx1 > N_idx2 ? -1 : 1;
					it++;
					LOG("Phenotype [%s-%s] will include [", S_dup, S_isH);
					for (K=N_idx1 ; ; K+=N_incre) {
						Ba_isIncPheno[K] = (char)j;
						//Sa_phenoName.push_back(*it);
						Xa_phenoNames.push_back(Xa_fColNames[K]);
						N_newPheno++;
						j++;

						xPheno V;
						V.S_name = Xa_fColNames[K];
						it = Xa_asgPhenos.insert(it, V);
						it++;

						/* Now end */
						if (K == (wsUint)N_idx2) {
							LOGnf("%s]\n", V.S_name.c_str());
							break;
						} else
							LOGnf("%s,", V.S_name.c_str());
					}
					/* Go back to original/delete/go back */
					it -= abs(N_idx1-N_idx2) + 2;
					it = Xa_asgPhenos.erase(it);
					it += abs(N_idx1-N_idx2);
				} else {
					char *S_isP = strchr(S_dup, '+');
					if (S_isP) {
						/* Check front exists */
						*(S_isP++) = '\0';

						FOREACHDO (vStr_it, Xa_fColNames, iit, i++) {
							/* Pass if not match */
							if (it->S_name.compare(*iit)) continue;
							if (Ba_isIncPheno[i])
								halt("Column [%s] in [%s] assigned twice in --pname!",
								it->S_name.c_str(), Sa_files);

							Ba_isIncPheno[i] = (char)j;
							//Sa_phenoName.push_back(*it);
							Xa_phenoNames.push_back(it->S_name);
							N_newPheno++;
							j++;
							break;
						}

					} else
						/* Otherwise, it is just an error */
						halt("Phenotype [%s] was not found on altpheno file", it->S_name.c_str());
				}
			}

			xPheno X;
			X.S_name = it->S_name;
			V_phenos.push_back(X);
		}
	}

	/* Build covariate information */
	if (S_covs && !strcmp("*", S_covs)) {
		wsUint i=0;
		/* Insert every columns if ALL variables are chosen to covariates */
		FOREACHDO (vStr_it, Xa_fColNames, it, i++) {
			/* Insert */
			xCovar ins;
			ins.Sp_varName	= strdup(it->c_str());
			ins.X_type		= WISARD_VAR_UNDET;		/* Set to undetermined */
			ins.N_szFactor	= 0;
			ins.N_idxFile	= i;
			ins.N_idx		= i;
			ins.R_times		= W1;
			ins.N_idxMul	= -1;
			ins.Sp_bl		= NULL;
			ins.S_fn		= Sa_files[(int)(Xm_fileColNames[(string)it->c_str()]/10000.0)];
			Xa_asgCovars.push_back(ins);
		}

		/* Assign 1 to N_fileCol sequentially to Ba_isIncCov */
		for (i=0 ; i<N_fileCol ; i++)
			Ba_isIncCov[i] = (char)i+1;
//		memset(Ba_isIncCov, 1, sizeof(char)*N_fileCol);
	} else {
		wsUint j=1,X=1;

		/* For each covariate name that assigned from option */
		FOREACHDO (vCovar_it, Xa_asgCovars, it, X++) {
			wsUint i=0;

			FOREACHDO (vStr_it, Xa_fColNames, iit, i++) {
				/* Pass if does not match */
				if (iit->compare(it->Sp_varName)) continue;

				/* Store the index to be stored (1-based) */
				if (Ba_isIncCov[i]) halt("Column [%s] in alternative phenotype"
					" assigned twice in --cname!", it->Sp_varName);
				Ba_isIncCov[i] = (char)X;
				it->N_idxFile = i;
				it->N_idx = j-1;
				j++;
//				Xa_covInfo.push_back(ins);
				break;
			}

			/* Found then skip */
			if (i != Xa_fColNames.size())
				continue;

			/* Check -(hyphen) character exists */
			char *S_dup	= strdup(it->Sp_varName);
			char *S_isH	= strchr(S_dup, '-');
			if (S_isH) {
				/* Split them into two, and check both exists */
				*(S_isH++) = '\0';
				/* Check NULLity */
				if (!S_dup[0] || !S_isH[0])
					halt("Covariate [%s-%s] was not found on altpheno file", S_dup, S_isH);
				int N_idx1 = -1;
				int N_idx2 = -1;
				wsUint K = 0;
				FOREACHDO (vStr_it, Xa_fColNames, iit, K++) {
					if (iit->compare(S_dup) == 0) N_idx1 = K;
					if (iit->compare(S_isH) == 0) N_idx2 = K;
				}
				/* Check both found */
				if (N_idx1 == -1)
					halt("Covariate [%s-%s] recognized ranged type, but first one does not exists",
						S_dup, S_isH);
				if (N_idx2 == -1)
					halt("Covariate [%s-%s] recognized ranged type, but second one does not exists",
						S_dup, S_isH);
				/* Check both are in same file */
				int N_idxS	= Xm_fileColNames[Xa_fColNames[N_idx1]];
				int N_idxE	= Xm_fileColNames[Xa_fColNames[N_idx2]];
				wsUint N_fS	= (wsUint)(N_idxS/10000.0);
				wsUint N_fE	= (wsUint)(N_idxE/10000.0);
				if (N_fS != N_fE)
					halt("Phenotype [%s] in [%s] but [%s] in [%s], "
					"both file must be match!", S_dup, Sa_files[N_fS],
					S_isH, Sa_files[N_fE]);

				/* Insert sequentially */
				wsUint N_incre = N_idx1 > N_idx2 ? -1 : 1;
				LOG("Covariate [%s-%s] will include [", S_dup, S_isH);
				//char *Sp_vn = it->Sp_varName;
				it++; /* Since insert() inserts element at front of itself */
				for (K=N_idx1 ; ; K+=N_incre) {
					Ba_isIncCov[K] = (char)(X++);
					//Sa_phenoName.push_back(*it);
//						Xa_phenoNames.push_back(it->S_name);
//						N_newPheno++;

					xCovar ins;
					ins.Sp_varName	= strdup(Xa_fColNames[K].c_str());
					ins.X_type		= WISARD_VAR_UNDET;		/* Set to undetermined */
					ins.N_szFactor	= 0;
					ins.N_idxFile	= i;
					ins.N_idx		= j - 1;
					ins.Sp_bl		= NULL;
					ins.N_idxMul	= -1;
					ins.R_times		= W1;
					it = Xa_asgCovars.insert(it, ins);
					it++;
					j++;

					/* Now end */
					if (K == (wsUint)N_idx2) {
						LOGnf("%s]\n", Xa_fColNames[K].c_str());
						break;
					} else
						LOGnf("%s,", Xa_fColNames[K].c_str());
				}
				it -= abs(N_idx1-N_idx2) + 2;
				it = Xa_asgCovars.erase(it);
				it += abs(N_idx1-N_idx2);
				X--;
			} else {
				/* FIXME : Add ^(number) type and *(field) and !(unary) type */
				char *Sp_isV = strchr(S_dup, '^');
				char *Sp_isS = strchr(S_dup, '*');
				char *Sp_isN = strchr(S_dup, '!');
				if (Sp_isV) {
					char *S_e = NULL;
					/* ^(number), number check */
					wsReal R_v = (wsReal)str2dbl(Sp_isV+1, &S_e);
					if (!S_e) goto _halt;
					*Sp_isV = '\0';

					/* Find again */
					wsUint k=0;
					FOREACHDO (vStr_it, Xa_fColNames, iit, k++) {
						/* Pass if does not match */
						if (iit->compare(S_dup)) continue;

						it->X_type		= WISARD_VAR_REAL;		/* Now they are must be REAL */
						it->N_szFactor	= 0;
						it->N_idxFile	= k;
						it->N_idx		= -1;
						it->Sp_bl		= NULL;
						it->N_idxMul	= -1;
						it->R_times		= R_v;
						break;
					}
					free(S_dup);
					continue;
				} else if (Sp_isS) {
					/* ^(number), number check */
					*(Sp_isS++) = '\0';

					/* Find again */
					wsUint k=0;
					FOREACHDO (vStr_it, Xa_fColNames, iit, k++) {
						/* Pass if does not match */
						if (iit->compare(S_dup)) continue;

						wsUint u=0;
						FOREACHDO (vStr_it, Xa_fColNames, iiit, u++) {
							/* Pass if does not match */
							if (iiit->compare(Sp_isS)) continue;

							it->X_type		= WISARD_VAR_REAL;		/* Now they are must be REAL */
							it->N_szFactor	= 0;
							it->N_idxFile	= k;
							it->N_idx		= -1;
							it->Sp_bl		= NULL;
							it->N_idxMul	= u;
							it->R_times		= W1;
							u = 0xffffffff;
							break;
						}
						if (u == 0xffffffff)
							break;
					}
					/* Successfully found A*B type */
					if (k != Xa_fColNames.size()) {
						free(S_dup);
						continue;
					}
				} else if (Sp_isN && Sp_isN == S_dup) {
					/* ! must be place at the first character */
					wsUint u=0;
					FOREACHDO (vStr_it, Xa_fColNames, iiit, u++) {
						/* Pass if does not match */
						if (iiit->compare(Sp_isN+1)) continue;
						
						it->N_idxMul	= -2;
						it->N_idxFile	= u;
						it->X_type		= WISARD_VAR_REAL;		/* Now they are must be REAL */
						it->N_szFactor	= 0;
						it->N_idx		= -1;
						it->Sp_bl		= NULL;
						it->R_times		= W1;
						u = 0xffffffff;
						break;
					}
					if (u == 0xffffffff) {
						free(S_dup);
						continue;
					}
				}
_halt:
				if (stricmp(it->Sp_varName, "SEX"))
					halt("Covariate [%s] was not found on alternative phenotype list", it->Sp_varName);
				else
					/* Index for sex */
					it->N_idx = j++;
			} /* END OF IF, hyphen or not */
			free(S_dup);
		} /* END OF FOREACH, Xa_asgCovars */
	}
	wsUint L=0,K=0;
	FOREACHDO (vCovar_it, Xa_asgCovars, it, L++)
		if (IS_REG_COVAR(it)) {
			LOG("Covariate [%s] stored in [%d]th column\n", it->Sp_varName, K+1);
			K++;
		}
	FOREACH (vCovar_it, Xa_asgCovars, it)
		if (!IS_REG_COVAR(it)) {
			it->N_idx = K; /* Can't set at this point cause don't know the size of factors */
			LOG("Covariate [%s] stored in [%d]th column\n", it->Sp_varName, ++K);
		}

	wsUintCst N_cov = Cp_IO->setCovBuf((wsUint)Xa_asgCovars.size());
	LOG("%d phenotype(s) and %d covariate(s) selected\n", N_newPheno, N_cov);
	if (N_newPheno == 0)
		LOG("No phenotypes were selected, original phenotype will be used\n");

	/* Get the number of samples included in the file,
	* and check each row is exists or not */
	char*		Ba_found		= NULL; ///< Whether (i)th sample found in sampvar file(1) or not(0)
	mDataIdx	Xm_iid2tmpIdx;
	wsUint		N_missSamp		= 0;
	mStr		Xm_flag;
	if (IS_ASSIGNED(sampvarflag)) loadFlags(OPT_STRING(sampvarflag), Xm_flag);
	char		B_nofid			= mapIsExist(Xm_flag, "nofid");
	wsCalloc(Ba_found, char, N_sample);
	for (wsUint Q=0 ; Q<N_files ; Q++) {
		char*		Ba_locFound		= NULL; ///< Whether (i)th sample found in sampvar file(1) or not(0)
		wsCalloc(Ba_locFound, char, N_sample);
		char		S_delim				= ' '; // 151228 delimiter change

		/* Read line... */
		for (wsUint L=1 ; Cp_pheno[Q]->gets(Sp_buf, 4096) ; L++) {
			/* Get FID and IID */
			getStringDelim(&Sp_buf, &b, S_delim);
			if (b == NULL) {
				S_delim = ','; // 151228 Try to change delimiter
				getStringDelim(&Sp_buf, &b, S_delim);
				if (b == NULL)
					halt_fmt(WISARD_INVL_FILE_INCMPL, "Sample variable file",
						Sa_files[Q], L);
			}
			/* Flag 'nofid' */
			if (!B_nofid) {
				getStringDelim(&b, &a, S_delim);
				if (a == NULL)
					halt_fmt(WISARD_INVL_FILE_INCMPL, "Sample variable file",
						Sa_files[Q], L);
			} else {
				a = b;
				b = Sp_buf;
			}

			/* Is this sample in the data? */
			mSamp_it X_resSampleFind = Xm_sample.find(b);
			if (X_resSampleFind == Xm_sample.end()) {
				if (++N_missSamp < 10)
					LOG("Sample [%s::%s] in [%s] is not exists in the dataset\n",
						Sp_buf, b, Sa_files[Q]);
				continue;
			}

			/* Check its FID is same with defined in original dataset */
			xSample	&X_samp = X_resSampleFind->second;
			if (X_samp.S_FID.compare(Sp_buf))
				halt("Sample [%s::%s] registered [%s] family in data",
					Sp_buf, b, X_samp.S_FID.c_str());

			/* Check whether it is filtered or not */
			if (X_samp.B_isMissing || X_samp.N_idx == SAMP_NODATA)
				continue;
			/* Mark as 'found */
			Ba_found[X_samp.N_idx] = 1;
			/* If dup, raise error */
			if (Ba_locFound[X_samp.N_idx])
				halt("Sample [%s::%s] in [%s] is duplicated!", Sp_buf, b,
					Sa_files[Q]);
			Ba_locFound[X_samp.N_idx] = 1;

			/* Assign index to be stored...anyway */
			mDataIdx_it X_res = Xm_iid2tmpIdx.find(X_samp.S_IID);
			if (X_res == Xm_iid2tmpIdx.end()) {
//				wsUint N_idx = Xm_iid2tmpIdx.size();
				Xm_iid2tmpIdx[X_samp.S_IID] = X_samp.N_idx;
			}
		}
		DEALLOC(Ba_locFound);
	}
	/* Incorrect number, some samples are lacked in pheno */
	if (Xm_iid2tmpIdx.size() != N_sample) {
		cExporter*	Cp_missamp	= cExporter::summon("pheno.missing.lst");
		vSampPtr	&Xa_samp	= Cp_IO->getSample();
		wsUint		N_lack		= 0;
		for (wsUint i=0 ; i<N_sample ; i++)
			if (Ba_found[i] == 0) {
				LOG("Sample [%s::%s] does not exists in phenotype file\n",
					Xa_samp[i]->S_FID.c_str(), Xa_samp[i]->S_IID.c_str());
				Cp_missamp->fmt("%s	%s\n", Xa_samp[i]->S_FID.c_str(), Xa_samp[i]->S_IID.c_str());
				N_lack++;
			}
		halt("[%d] samples are lacked in phenotype file", N_lack);
		delete Cp_missamp;
	}
	DEALLOC(Ba_found);

	/* Reallocate phenotype buffer */
	if (N_newPheno)
		Cp_IO->setPhenoBuf(N_newPheno);

	/* Factor-type buffer for Fst group */
	mDataIdx	Xm_grpFst;

	/* Factor type buffer */
	vector<mDataIdx> Xm_factor;
	for (wsUint i=0 ; i<N_cov ; i++)
		Xm_factor.push_back(mDataIdx());

	/* Count # of column in each file */
	vInt Xa_nCol(2);
	FOREACH (mDataIdx_it, Xm_fileColNames, i)
		Xa_nCol[(int)(i->second/10000.0)]++;

	/* Retrieve */
	char ***Sa_rawPheno = NULL;
	wsAlloc(Sa_rawPheno, char**, N_sample);
	for (wsUint i=0 ; i<N_sample ; i++) {
		wsCalloc(Sa_rawPheno[i], char*, Xa_fColNames.size());
	}
	/* For all files */
	vector<char*> Xa_bufs;
	N_missSamp = 0;
	for (wsUint Q=0 ; Q<N_files ; Q++) {
		char S_delim = ' ';
//		wsAlloc(Sp_buf, char, 8192);
		Cp_pheno[Q]->rewind();

		if (!OPT_ENABLED(nophenohdr))
			Cp_pheno[Q]->gets(Sp_buf, 8192);
// 			halt("Can't read phenotype header!");
// 
		/* Set the base */
		wsUint N_offset = 0;
		for (wsUint L=0 ; L<Q ; L++) N_offset += Xa_nCol[L];

		for (wsUint L=1 ; Cp_pheno[Q]->gets(Sp_buf, 8192) ; Sp_buf[0]='\0',L++) {
			/* Get FID and IID */
			getStringDelim(&Sp_buf, &b, S_delim);
			if (b == NULL) {
				S_delim = ',';
				getStringDelim(&Sp_buf, &b, S_delim);
			}
			/* Flag 'nofid' */
			if (!B_nofid)
				getStringDelim(&b, &a, S_delim);
			else {
				a = b;
				b = Sp_buf;
			}

			/* Get sample data */
			string key = (string)b;
			mSamp_it X_resSampleFind = Xm_sample.find(key);
			if (X_resSampleFind == Xm_sample.end()) {
				if (++N_missSamp < 10)
					B_nofid ? LOG("IID [%s] not exists in the dataset\n", Sp_buf) :
						LOG("[%s::%s] not exists in the dataset\n", Sp_buf, b);
				continue;
			}
			xSample	&X_samp = X_resSampleFind->second;

			/* Give the phenotype unless it is not filtered */
			if (X_samp.B_isMissing || X_samp.N_idx == SAMP_NODATA)
				continue;

			/* Get storing index */
			int N_idxStore = Xm_iid2tmpIdx[key];
			/* Retrieve current data */
			int j=0;
			for ( ; j<Xa_nCol[Q] ; j++,a=b) {
				if (a == NULL)
					halt_fmt(WISARD_INVL_FILE_INCMPL, "Sample variable file",
						Sa_files[Q], L);
// 					halt("Sample [%s::%s] have insufficient data",
// 						X_samp.S_FID.c_str(), X_samp.S_IID.c_str());
				getStringDelim(&a, &b, S_delim);

				/* Retrieve anyway */
				Sa_rawPheno[N_idxStore][N_offset+j] = a;
			}
			if (b != NULL) {
				getStringDelim(&a, &b, S_delim);
				if (a[0])
					halt("[%s::%s] have invalid number of phenotypes in [%s] at line [%d]",
						X_samp.S_FID.c_str(), X_samp.S_IID.c_str(),
						Sa_files[Q], L);
			}

			Xa_bufs.push_back(Sp_buf);
			wsAlloc(Sp_buf, char, 8192);
		}

		Cp_pheno[Q]->close();
		delete Cp_pheno[Q];
	}
	DEALLOC(Cp_pheno);

	/* If --filsample, set variable and check */
	vBool&		Bv_filtSamp		= Cp_IO->getIsSampleFiltered();
	vSampPtr&	Xa_samp				= Cp_IO->getSample();
	wsUint&		N_szFilteredSamp	= Cp_IO->getFiltSampSize();
	if (IS_ASSIGNED(filsample) || IS_ASSIGNED(incsample)) {
		xOperation*	Xp_orig	= IS_ASSIGNED(filsample) ? &(OPT_EXPR(filsample)) :
			&(OPT_EXPR(incsample));
		for (wsUint i=0 ; i<N_sample ; i++) {
			xOperation *Xp_o = Xp_orig->clone();
			for (wsUint j=0 ; j<N_fileCol ; j++) {
				xOperand *Xp_op = str2op(Sa_rawPheno[i][j]);
				Xp_o->assign(Xa_fColNames[j].c_str(), Xp_op);
			}
			xOperand *Xp_res = understand(Xp_o, 1);
			if (Xp_res->X_type != OTP_LOGICAL)
				halt("Result of [%s] is not logical for sample [%s::%s]",
					IS_ASSIGNED(filsample)?"--filsample":"--incsample",
					Xa_samp[i]->S_FID.c_str(),
					Xa_samp[i]->S_IID.c_str());
			bool B_filt = IS_ASSIGNED(incsample) ? !Xp_res->B_val : Xp_res->B_val;

			/* Should be filtered */
			if (B_filt && Bv_filtSamp[i] == 0) {
				N_szFilteredSamp++;
				Bv_filtSamp[i] = 1;
			}
		}
	}

	/* Deallocate */
	for (wsUint i=0 ; i<N_files ; i++)
		free(Sa_files[i]);
	DEALLOC(Sa_files);

	/* Flip check */
	for (wsUint j=0 ; j<Xa_asgCovars.size() ; j++) {
		xCovar &X_cov = Xa_asgCovars[j];

		if (X_cov.N_idxMul == -2) {
			/* ! operator */

			/* Check # of value category is 2 */
			map<wsReal,char> Xm_cat;
			for (wsUint i=0 ; i<N_sample ; i++)
				if (stricmp(OPT_STRING(mispheno), Sa_rawPheno[i][X_cov.N_idxFile])) {
					wsReal R_v = (wsReal)atof(Sa_rawPheno[i][X_cov.N_idxFile]);
					Xm_cat[R_v] = 1;

					/* Size check */
					if (Xm_cat.size() > 2)
						halt("More than two value found for column [%s] in "
							"[%s], ! operator cannot be applied",
							Xa_asgCovars[j].Sp_varName,
							Xa_asgCovars[j].S_fn);
				}

			/* Flip value */
			X_cov.R_v0 = Xm_cat.at(0);
			X_cov.R_v1 = Xm_cat.at(1);
		}
	}

	/* if --sampleweight, assign memory */
	wsVec Ra_sampWgt = NULL;
	if (IS_ASSIGNED(sampleweight) && !Cp_IO->getSampleWeight()) {
		/* Init with NAN so the system recognizes lack */
		Ra_sampWgt = Cp_IO->setSampleWeightBuf();
	}

	/* Now process */
	wsMat	Ra_pheno			= Cp_IO->getPhenos();
	wsMat	Ra_cov				= Cp_IO->getCovariates();
	wsUint	N_twin				= 0;
	vStr&	Sv_sampGrpIdx2grp	= Cp_IO->getSampGrpIdx2grp();
	char*	Ba_typePheno		= Cp_IO->getTypePheno();
	for (wsUint i=0 ; i<N_sample ; i++) {
		xSample &X_samp = *(Xa_samp[i]);
		wsUint	Jc = 0;

		for (wsUint j=0 ; j<N_fileCol ; j++) {
			wsUint N_xIdx = Xm_fileColNames[Xa_fColNames[j]];
			a = Sa_rawPheno[i][j];

			if (Xv_resvIdx[RC_TWIN] && (int)(N_xIdx+1) == Xv_resvIdx[RC_TWIN]) {
//			if (N_idxTwin && (N_xIdx+1) == N_idxTwin) {
				if (!a) halt("Sample [%s::%s] does not have twin information!",
					X_samp.S_FID.c_str(), X_samp.S_IID.c_str());

				/* TWIN id should consisted of numbers */
				for (char *p=a ; *p ; p++)
					if (*p < '0' || *p > '9') halt("Sample [%s::%s] have invalid twin id [%s]",
						X_samp.S_FID.c_str(), X_samp.S_IID.c_str(), a);
				wsUint N_idTwin = atoi(a);

				/* Currently twin id is restricted to not over MAX_TWINS-1 */
				if (N_idTwin >= MAX_TWINS)
					halt("Too large twin ID >= %d!", MAX_TWINS);
				if (N_idTwin > N_twin)
					N_twin = N_idTwin;
				/* Registering twin id */
				X_samp.N_idTwin = (char)N_idTwin;

				continue;
			} else if (Xv_resvIdx[RC_PROBAND] && (int)(N_xIdx+1) == Xv_resvIdx[RC_PROBAND]) {
				//if (N_idxProband && (N_xIdx+1) == N_idxProband) {
				if (!a) halt("Sample [%s::%s] does not have proband information!",
					X_samp.S_FID.c_str(), X_samp.S_IID.c_str());

				/* Proband status should be 0 (not proband) or 1 (proband) */
				if (a[1] || (a[0] != '0' && a[0] != '1'))
					halt("Sample [%s::%s] have invalid proband status [%s] :"
						" should be 0 or 1", X_samp.S_FID.c_str(),
						X_samp.S_IID.c_str(), a);
				X_samp.B_isProband = a[0]-'0';

				continue;
			} else if (Xv_resvIdx[RC_POPGROUP] && (int)(N_xIdx+1) == Xv_resvIdx[RC_POPGROUP]) {
				//if (N_idxFst && (N_xIdx+1) == N_idxFst) {
				if (!a) halt("Sample [%s::%s] does not have sample group information!",
					X_samp.S_FID.c_str(), X_samp.S_IID.c_str());

				/* Group status should be factor-type */
				string _a = (string)a;
				mDataIdx_it X_find = Xm_grpFst.find(_a);
				if (X_find == Xm_grpFst.end()) {
					Sv_sampGrpIdx2grp.push_back(_a);
					Xm_grpFst.insert(make_pair(_a, (wsUint)Xm_grpFst.size()));
				}
				pverbose("Group [%s] found\n", a);
				X_samp.N_grpFst = Xm_grpFst[_a];

				continue;
			} else if (Xv_resvIdx[RC_SAMPLEWEIGHT] && (int)(N_xIdx+1) == Xv_resvIdx[RC_SAMPLEWEIGHT]) {
				if (!a) halt("Sample [%s::%s] does not have sample-wise weight!",
					X_samp.S_FID.c_str(), X_samp.S_IID.c_str());
				char *tmp = NULL;
				Ra_sampWgt[i] = str2dbl(a, &tmp);
				if (tmp[0]) halt("Sample [%s::%s] have invalid sample weight [%s] defined in the column [%s], cannot be this value!",
					X_samp.S_FID.c_str(), X_samp.S_IID.c_str(), a, Xa_fColNames[Xv_resvIdx[RC_SAMPLEWEIGHT] - 1].c_str());
				Ra_sampWgt[i] = W1/Ra_sampWgt[i];
				/* FIXME : Implementation needed */
			}

			/* If this column assigned to phenotype */
			if (Ba_isIncPheno[j]) {
				wsUint N_idxIns = Ba_isIncPheno[j] - 1;

				if (!a || !stricmp(OPT_STRING(mispheno), a)) {
					Ra_pheno[N_idxIns][X_samp.N_idx] = WISARD_NA_REAL;
					/* --filmispheno */
					if (OPT_ENABLED(filmispheno)) {
						if (Bv_filtSamp[X_samp.N_idx] == 0)
							N_szFilteredSamp++;
						Bv_filtSamp[X_samp.N_idx] = 1;
					}
				} else {
					char*	p		= NULL;
					wsReal	R_pheno	= (wsReal)str2dbl(a, &p);
					if (p && p[0])
						halt("Phenotype [%s] have invalid value [%s]",
							Xa_fColNames[j].c_str(), a);

					/* Marking continuous if it is continuous */
					if (R_pheno != OPT_NUMBER(phenoCtrl) &&
						R_pheno != OPT_NUMBER(phenoCase))
						Ba_typePheno[N_idxIns] = 1;
					/* Filtering option */
					else if ((OPT_ENABLED(filcase) && R_pheno == OPT_NUMBER(phenoCase)) ||
						(OPT_ENABLED(filcontrol) && R_pheno == OPT_NUMBER(phenoCtrl)))
						Bv_filtSamp[X_samp.N_idx] = 1;

					if (R_pheno != R_pheno)
						Ra_pheno[N_idxIns][X_samp.N_idx] = WISARD_NA_REAL;
					else
						Ra_pheno[N_idxIns][X_samp.N_idx] = R_pheno;
				}
			}

			/* If this column assigned to covariates */
			/* FIXME : Multiple covariates assignment is possible by ^ and * modifier so fix is needed */
			if (Ba_isIncCov[j]) {
				wsUint N_idxIns = Ba_isIncCov[j] - 1;

				/* If missing */
				if (!a || !stricmp(OPT_STRING(mispheno), a)) {
					Ra_cov[N_idxIns][X_samp.N_idx]	= WISARD_NA_REAL;
					Jc++;
					continue;
				}

				/* Try to convert with wsReal */
				char*		Sp_end	= NULL;
				wsReal		R_res	= (wsReal)str2dbl(a, &Sp_end);
				xCovar&		X_curC	= Xa_asgCovars[N_idxIns];
				mDataIdx&	Xm_curF	= Xm_factor[N_idxIns];
				N_idxIns = X_curC.N_idx;

				/* Set the type of covariate if not determined */
				if (X_curC.X_type == WISARD_VAR_UNDET) {
					X_curC.X_type = Sp_end[0] != '\0' ? WISARD_VAR_FACTOR :
						WISARD_VAR_REAL;
					/* NA check */
					if (!stricmp(a, "NA"))
						LOGwarn("Value [NA] found for covariate [%s] in [%s], "
							"it will considered as `FACTOR` but is it true?\n",
							X_curC.Sp_varName, X_curC.S_fn);
					pverbose("Covariate [%s] recognized as [%s]\n",
						X_curC.Sp_varName,
						Sp_end[0] != '\0' ? "FACTOR" : "REAL NUMBER");
				}


				switch (X_curC.X_type) {
				case WISARD_VAR_FACTOR: { /* If this is factor, Sp_end == '\0' is not allowed */
// 						if (Sp_end[0] == '\0')
// 							halt("Covariate [%s] recognized as factor, but "
// 								"real number [%s] found for sample [%s::%s]",
// 								X_curC.Sp_varName, a, X_samp.S_FID.c_str(),
// 								X_samp.S_IID.c_str());

						/* Check this factor value is already inserted */
						string	S_factor = (string)a;
						mDataIdx_it
								X_idxFactor = Xm_curF.find(S_factor);
						if (X_idxFactor == Xm_curF.end()) {
							/* Insert with new index */
							Xm_curF.insert(
								make_pair(S_factor, X_curC.N_szFactor));
							X_curC.N_szFactor++;
						}
						/* Insert factor value */
						wsUint N_idxFactor = Xm_curF[S_factor];
						Ra_cov[N_idxIns][X_samp.N_idx] = (wsReal)N_idxFactor;

						/* FIXME : Value should be stored as factor type */
					} break;
				case WISARD_VAR_REAL: /* If this is real, Sp_end != '\0' is not allowed */
					/* Notice : MUST NOT it have Sp_bl */
					if (X_curC.Sp_bl)
						halt("Covariate [%s] recognized as real number, but "
							"--baseline has assigned to determine baseline value",
							X_curC.Sp_varName);
					/* If it is real number, factor-like value are not permitted 
					* (which have surplus character following by the end of real number) */
					if (Sp_end[0] != '\0')
						halt("Covariate [%s] recognized as real number, but "
							"factor-like [%s] found for sample [%s::%s]",
							X_curC.Sp_varName, a, X_samp.S_FID.c_str(),
							X_samp.S_IID.c_str());

					/* Nan processing */
					if (R_res != R_res)
						Ra_cov[N_idxIns][X_samp.N_idx] = WISARD_NA_REAL;
					else
						Ra_cov[N_idxIns][X_samp.N_idx] = R_res;
					break;
				default: /* Other values are not permitted */
					halt("SYSERR: Covariate [%s] have invalid status %d",
						X_curC.Sp_varName, X_curC.X_type);
					break;
				} /* END OF switch */

				Jc++;
			}
		} /* END OF for */

		for (wsUint j=0 ; j<Xa_asgCovars.size() ; j++) {
			xCovar &X_cov = Xa_asgCovars[j];

			/* Special meaning */
			if (X_cov.R_times != W1) {
				/* Check is it factor */
				if (Ba_isIncCov[X_cov.N_idxFile]) {
					xCovar	&X_curC	= Xa_asgCovars[Ba_isIncCov[X_cov.N_idxFile]-1];
					if (X_curC.X_type == WISARD_VAR_FACTOR)
						halt("Covariate [%s] is factor, so [%s] cannot be applied",
							X_curC.Sp_varName, X_cov.Sp_varName);
				}

				if (!Sa_rawPheno[i][X_cov.N_idxFile]) {
					Ra_cov[Jc++][X_samp.N_idx] = WISARD_NA_REAL;
				} else {
					wsReal R_v = (wsReal)atof(Sa_rawPheno[i][X_cov.N_idxFile]);
					if (!stricmp(OPT_STRING(mispheno), Sa_rawPheno[i][X_cov.N_idxFile]))
						Ra_cov[Jc++][X_samp.N_idx] = WISARD_NA_REAL;
					else
						Ra_cov[Jc++][X_samp.N_idx] = pow(R_v, X_cov.R_times);
				}
			} else if (Xa_asgCovars[j].N_idxMul == -2) {
				/* ! operator */
				
				/* Check is it factor */
				if (Ba_isIncCov[X_cov.N_idxFile]) {
					xCovar	&X_curC	= Xa_asgCovars[Ba_isIncCov[X_cov.N_idxFile]-1];
					if (X_curC.X_type == WISARD_VAR_FACTOR)
						halt("Covariate [%s] is factor, so [%s] cannot be applied",
							X_curC.Sp_varName, X_cov.Sp_varName);
				}

				/* Flip value */
				if (!Sa_rawPheno[i][X_cov.N_idxFile]) {
					Ra_cov[Jc++][X_samp.N_idx] = WISARD_NA_REAL;
				} else {
					wsReal R_v = (wsReal)atof(Sa_rawPheno[i][X_cov.N_idxFile]);
					if (!stricmp(OPT_STRING(mispheno), Sa_rawPheno[i][X_cov.N_idxFile]))
						Ra_cov[Jc++][X_samp.N_idx] = WISARD_NA_REAL;
					else
						Ra_cov[Jc++][X_samp.N_idx] = R_v == X_cov.R_v0 ? X_cov.R_v1 : X_cov.R_v0;
				}
			} else if (Xa_asgCovars[j].N_idxMul != -1) {
				/* Check is it factor */
				if (Ba_isIncCov[X_cov.N_idxFile]) {
					wsUint N_idxIns = Ba_isIncCov[X_cov.N_idxFile] - 1;
					xCovar	&X_curC	= Xa_asgCovars[N_idxIns];
					if (X_curC.X_type == WISARD_VAR_FACTOR)
						halt("Covariate [%s] is factor, so [%s] cannot be applied",
						X_curC.Sp_varName, X_cov.Sp_varName);
				}
				if (Ba_isIncCov[X_cov.N_idxMul]) {
					wsUint N_idxIns = Ba_isIncCov[X_cov.N_idxMul] - 1;
					xCovar	&X_curC	= Xa_asgCovars[N_idxIns];
					if (X_curC.X_type == WISARD_VAR_FACTOR)
						halt("Covariate [%s] is factor, so [%s] cannot be applied",
						X_curC.Sp_varName, X_cov.Sp_varName);
				}

				if (!Sa_rawPheno[i][X_cov.N_idxFile] || !Sa_rawPheno[i][X_cov.N_idxMul]) {
					Ra_cov[Jc++][X_samp.N_idx] = WISARD_NA_REAL;
				} else {
					wsReal R_v1 = (wsReal)atof(Sa_rawPheno[i][X_cov.N_idxFile]);
					wsReal R_v2 = (wsReal)atof(Sa_rawPheno[i][X_cov.N_idxMul]);
					if (!stricmp(OPT_STRING(mispheno), Sa_rawPheno[i][X_cov.N_idxFile]) ||
						!stricmp(OPT_STRING(mispheno), Sa_rawPheno[i][X_cov.N_idxMul]))
						Ra_cov[Jc++][X_samp.N_idx] = WISARD_NA_REAL;
					else
						Ra_cov[Jc++][X_samp.N_idx] = R_v1*R_v2;
				}
			}
		}
	}
	/* Now the # of twins is fixed */
	Cp_IO->setTwinSize(N_twin);

	/* Deallocate */
	wsUint X=0;
	for (wsUint i=0 ; i<N_sample ; i++) {
		DEALLOC(Sa_rawPheno[i]);
	}
	DEALLOC(Sa_rawPheno);
//	FOREACHDO (vector<vector<char*> >::iterator, Xa_rawPheno, it1, X++) {
// 		wsUint Y=0;
// 		FOREACHDO (vector<char*>::iterator, *it1, it2, Y++)
// 			delete *it2;
// 	}
 	FOREACHDO (vector<char*>::iterator, Xa_bufs, it1, X++)
		delete *it1;


	/* Check integrity */
// 	if (N_fAltSample != N_fSamp)
// 		halt("# of samples in alternative phenotype file is not match with"
// 			" original dataset");
	/* Check integrity for sample weight */
	if (IS_ASSIGNED(sampleweight) ) {

	}

	/* Print status if phenotype is altered */
	wsUint N_pheno = Cp_IO->sizePheno();
	if (N_newPheno)
		LOG("%d*%d alternative phenotypes were retrieved\n", N_sample,
			N_pheno);

	/* Print summary */
	wsUint Jc = 0, Jp = 0;
	wsUint N_newCov = 0;
	for (wsUint i=0 ; i<N_fileCol ; i++) {
		if (Ba_isIncPheno[i]) {
			wsUint	N_idxIns	= Ba_isIncPheno[i] - 1;

			LOG("Phenotype [%s] have type `REAL NUMBER`\n",
				Xa_phenoNames[N_idxIns].c_str());
			Jp++;
		}
	}
	wsUint N_fac = 0;
	FOREACH (vCovar_it, Xa_asgCovars, it) {
		xCovar		&X_curC		= *it;
		mDataIdx	&X_curF		= Xm_factor[Jc];

		switch (X_curC.X_type) {
		case WISARD_VAR_FACTOR:
			if (X_curC.Sp_bl) /* If user assigned specific value to be a baseline, print it instead */
				LOG("Covariate [%s] have type [FACTOR] with %d levels "
					"(Baseline %s)\n", X_curC.Sp_varName,
					X_curC.N_szFactor, X_curC.Sp_bl);
			else /* Otherwise, first element is baseline */
				LOG("Covariate [%s] have type [FACTOR] with %d levels "
					"(Baseline %s)\n", X_curC.Sp_varName,
					X_curC.N_szFactor, (X_curF.begin())->first.c_str());

			N_newCov += X_curC.N_szFactor-1;
			N_fac++;
			break;
		case WISARD_VAR_REAL:
			LOG("Covariate [%s] have type [REAL NUMBER]\n", X_curC.Sp_varName);
			N_newCov++;
			break;
		default:
			/* Special treatment for 'SEX' */
			if (!stricmp(X_curC.Sp_varName, "SEX")) {
				LOGnote("Covariate [SEX] will use internal value to [%d]\n",
					X_curC.N_idx);
				X_curC.X_type		= WISARD_VAR_REAL;
				/* Allocate! */
				FOREACH (vSampPtr_it, Xa_samp, iit) {
					wsUint N_sex = (*iit)->N_sex;
					Ra_cov[X_curC.N_idx - 1][(*iit)->N_idx] = N_sex ?
						(N_sex - 1) : WISARD_NA;
				}
				N_newCov++;
			} else
				halt("SYSERR: Covariate type [%d] for [%s] should be "
					"fixed at this point, retrieved data [%d] samples",
					X_curC.X_type, X_curC.Sp_varName, N_sample);
			break;
		}

		Jc++;
	}

	/* Rearrange covariates if factor is assigned as covariate */
	vCovar&	Xa_covInfo		= Cp_IO->getCovInfo();
	wsUint	N_altSample	= N_sample;
	if (N_fac) {
		wsReal	**Ra_newCov = NULL;
		vCovar	Xa_newCI;
		wsAlloc(Ra_newCov, wsReal*, N_newCov);

		Jc = 0;
		wsUint N_idxNew = 0;
		for (wsUint i=0 ; i<Xa_asgCovars.size() ; i++) {
			/* Find corresponding to ith covariates IN file's perspective */
			xCovar	&X_curC		= Xa_asgCovars[i];
			mDataIdx X_curF		= Xm_factor[i];
			wsUint	 N_idxIns	= X_curC.N_idx;

			if (X_curC.X_type == WISARD_VAR_FACTOR) {
				/* Allocate Nfactor-1 rows */
				wsUint N_st = N_idxNew, j;
				for (j=1 ; j<X_curC.N_szFactor ; j++) {
					sseCalloc(Ra_newCov[N_idxNew+j-1], wsReal, N_sample);
				}
				N_idxNew += X_curC.N_szFactor-1;

				size_t N_lenOrigVar = strlen(X_curC.Sp_varName);
				j = N_st;
				wsUint L = j;
				/* Current covariate is factor,
				* will occupy Nfactor-1 slots */
				if (X_curC.Sp_bl) {
					/* Find which index fits to the given baseline */
					mDataIdx_it it=X_curF.begin();
					for ( ; it!=X_curF.end() ; it++,L++) {
						/* Match found */
						if (it->first.compare(X_curC.Sp_bl) == 0)
							break;
					}
					if (it == X_curF.end()) halt("--baseline value [%s] for"
						" covariate column [%s] not found in the possible"
						" values of covariates", X_curC.Sp_bl, X_curC.Sp_varName);
				} /* If there is no --baseline, first element will baseline */

				wsUint m=j;
				FOREACHDO (mDataIdx_it, X_curF, it, j++) {
					/* Skip first row but make CI info */
					if (j == L) continue;
					wsReal R_valTarget = (wsReal)(it->second);

					/* For all sample, marking */
					for (wsUint k=0 ; k<N_sample ; k++) {
						if (isMissingReal(Ra_cov[N_idxIns][k]))
							Ra_newCov[m][k] = WISARD_NA_REAL;
						else if (R_valTarget == Ra_cov[N_idxIns][k])
							Ra_newCov[m][k] = W1;
					}

					/* Making new covariate information */
					Xa_newCI.push_back(X_curC);
					Xa_newCI[m].N_idx = m;
					wsAlloc(Xa_newCI[m].Sp_varName, char,
						N_lenOrigVar+2+it->first.length());
					sprintf(Xa_newCI[m].Sp_varName, "%s=%s",
						X_curC.Sp_varName, it->first.c_str());
					m++;
				}

				/* Delete previous variable name */
				DEALLOC(X_curC.Sp_varName);

				/* Delete original buffer */
//				sseFree(Ra_cov[N_idxIns]);
			} else {
				/* Current covariate is real, copy entire from original */
				Ra_newCov[N_idxNew] = sseVectorP(N_sample, Ra_cov[N_idxIns]);

				Xa_newCI.push_back(X_curC);
				Xa_newCI[N_idxNew].N_idx	= N_idxNew;

				N_idxNew++;
			}

			/* Increase original index */
			Jc++;
		}

		Cp_IO->setCovBuf(N_newCov, Ra_newCov);

		/* Replace Xa_covInfo */
		//Xa_covInfo.resize(N_newCov);
		FOREACH (vCovar_it, Xa_newCI, it)
//			if (it->N_idxMul == -1 && it->R_times == W1)
				Xa_covInfo.push_back(*it);
// 		FOREACH (vCovar_it, Xa_newCI, it)
// 			if (it->N_idxMul != -1 || it->R_times != W1)
// 				Xa_covInfo.push_back(*it);
		//copy(Xa_newCI.begin(), Xa_newCI.end(), Xa_covInfo.begin());

		if (OPT_ENABLED(verbose)) {
			wsReal **Ra_covT = transpose(Ra_cov, N_cov, N_altSample);

			vStr Sa_covNames;
			Cp_IO->getCovNames(Sa_covNames);

			exportMatrix("final.cov", Ra_covT, N_altSample, N_cov, &Sa_covNames);
			sseUnmat(Ra_covT, N_altSample);
			LOG("Final covariates matrix exported\n");
		}
	} else {
//		Xa_covInfo.resize(N_cov);
		FOREACH (vCovar_it, Xa_asgCovars, it)
			if (it->N_idxMul == -1 && it->R_times == W1)
				Xa_covInfo.push_back(*it);
		FOREACH (vCovar_it, Xa_asgCovars, it)
				if (it->N_idxMul != -1 || it->R_times != W1)
					Xa_covInfo.push_back(*it);
//		copy(Xa_asgCovars.begin(), Xa_asgCovars.end(), Xa_covInfo.begin());
	}
	/* Now the number of altsample is fixed */
	Cp_IO->setAltSampSize(N_altSample);

	/* Check --gxecovs */
	vInt	Xa_gxeIdx;
	if (IS_ASSIGNED(gxecovs)) {
		wsUint	N_gc		= 0;
		char**	Sa_gxeCovs	= loadStringValues(OPT_STRING(gxecovs), &N_gc);

		/* Check those variables */
		for (wsUint i=0 ; i<N_gc ; i++) {
			char*	Sp_cur	= Sa_gxeCovs[i];
			wsUint	N_found	= 0;
			FOREACH (vCovar_it, Xa_covInfo, it)
				if (!stricmp(it->Sp_varName, Sp_cur)) {
					Xa_gxeIdx.push_back(it->N_idx);
					N_found = 1;
					break;
				}
			if (N_found) continue;

			/* If not found, check range-type */
			char *Sp_div = strchr(Sp_cur, '-');
			if (!Sp_div) halt("Column name [%s] in --gxecovs does not "
				"exists in the column names from --cname", Sp_cur);
			*(Sp_cur++) = '\0';
			/* Check two variables are existing */
			int N_idxS = -1;
			int N_idxE = -1;
			FOREACH (vCovar_it, Xa_covInfo, it) {
				if (!stricmp(it->Sp_varName, Sp_div)) {
					Xa_gxeIdx.push_back(it->N_idx);
					N_idxS = it->N_idx;
				}
				if (!stricmp(it->Sp_varName, Sp_cur)) {
					Xa_gxeIdx.push_back(it->N_idx);
					N_idxE = it->N_idx;
				}
			}
			/* Determine direction and insert them */
			int N_move = N_idxS <= N_idxE ? 1 : -1;
			for (int j=N_idxS ; ; j+=N_move) {
				Xa_gxeIdx.push_back(j);
				/* Stop can't work anymore */
				if (j == N_idxE) break;
			}
		}
	} else if (OPT_ENABLED(gxe)) {
		/* Use all coviarates as gxe */
		wsUint k = 0;
		FOREACHDO (vCovar_it, Xa_covInfo, it, k++)
			Xa_gxeIdx.push_back(k);
	}
	Cp_IO->setGxeBuf(Xa_gxeIdx);

	DEALLOC(Sp_buf);
	DEALLOC(Ba_isIncPheno);
	DEALLOC(Ba_isIncCov);
	DEALLOC(S_phenos);
	DEALLOC(S_covs);

}

} // End namespace ONETOOL
