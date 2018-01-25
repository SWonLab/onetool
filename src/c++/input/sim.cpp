#include "global/common.h"
#include "utils/vector.h"
#include "input/stream.h"
#include "global/option.h"
#include "utils/stat.h"
#include "input/sim.h"
#include "analyses/pddt.h"
#include "utils/matrix.h"
#ifdef _WIN32
#	include <io.h>
#else
#	include <dirent.h>
#endif

namespace ONETOOL {

int genMarker(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	cSimFamIO	*Cp_IO		= (cSimFamIO *)Vp_shareData;
	wsUint		N_sample	= Cp_IO->sizeSample();
	wsUint		N_SNP		= Cp_IO->sizeVariant();
	wsUint		N_sig		= OPT_NUMBER(nsig);
	vector<vInt>&
				Xv_prop		= Cp_IO->getPropIndices();
	vReal&		Xv_maf		= Cp_IO->getMAFs();

	int*		Np_data		= (int *)Vp_data;
	wsUint		N_s			= (wsUint)Np_data[0];
	wsUint		N_e			= (wsUint)Np_data[1];
	wsReal		R_beta		= WISARD_NAN;
	wsReal		R_sigmaf	= OPT_REAL(sigmaf);
	wsReal		R_p1p		= R_sigmaf * (W1 - R_sigmaf);
	if (IS_ASSIGNED(beta))
		R_beta = OPT_REAL(beta);
	else
		R_beta	= sqrt(0.005 / (REAL_CONST(1.99)*R_p1p));

	//LOG("Thread %d : Range %d~%d\n", N_idx, N_s, N_e);
	for (wsUint j=N_s ; j<N_e ; j++) {
		cVector V_y(Cp_IO->getPhenos()[0], N_sample, 1);
		if (j < N_sig) while (!Cp_IO->_simsig(Xv_prop, j, R_sigmaf, Cp_IO->getV(),
			V_y, R_beta));
		else while (!Cp_IO->_sim(Xv_prop, j, Xv_maf[j]));

		if (N_idx == 0 && (j%1000) == 0) {
			int S = 0;
			for (int i=0 ; i<OPT_NUMBER(thread) ; i++)
				S += Np_data[i*3+2];
			notice("%d/%d SNPs processed...\r", S, N_SNP);
		}
		Np_data[2]++;
	}

	return 0;
}

char* _findfile_sub(char *Sp_fn)
{
	/* Must staring with .hap and .pos */
	char *Sp_h = strstr(Sp_fn, ".hap");
	char *Sp_p = strstr(Sp_fn, ".pos");
	char *Sp_d = Sp_h ? Sp_h : Sp_p;
	if (Sp_d) {
		/* Before it is equivalent to the non-path part of S_prefix? */
		if (memcmp(Sp_d, Sp_fn, Sp_d-Sp_fn)) return NULL;

		return Sp_d + 5;
	}

	return NULL;
}

char** findCOSI(char *S_prefix, wsUint *Np_sz)
{
	mDataIdx Xm_pop;
#ifdef _WIN32
	/* Extract non-path part */
	char *Sp_e = S_prefix + strlen(S_prefix) - 1;
	while (Sp_e >= S_prefix && *Sp_e != '\\') Sp_e--;
	_finddata_t	X_fd;
	intptr_t	H_find;
	int			N_res = 0;
	char		S_buf[MAX_PATH + 3];
	sprintf(S_buf, "%s.*", S_prefix);
	H_find = _findfirst(S_buf, &X_fd);
	if (H_find == -1)
		return NULL;

	while (N_res != -1) {
		/* Pass if directory */
		if (X_fd.attrib & _A_SUBDIR) continue;

		char *Sp_pop = _findfile_sub(X_fd.name);
		if (Sp_pop) Xm_pop[Sp_pop]++;
	}
#else
	DIR *dp;
	struct dirent *dirp;

	/* Get directory point */
	char *ll = strdup(S_prefix);
	char *lp = ll + strlen(ll) - 1;
	while (ll<=lp && *lp != '/') lp--;
	if (ll<=lp) *lp = '\0';

	if ((dp = opendir(ll)) == NULL) {
		printf("Failed to open directory [%s]", ll);
		exit(1);
	}
	while ((dirp = readdir(dp)) != NULL) {
		if(!strcmp(dirp->d_name, lp+1))
			printf("%s\n", dirp->d_name);
	}
	closedir(dp);
#endif
	/* Pick out all 2's */
	for (mDataIdx_it i=Xm_pop.begin() ; i!=Xm_pop.end() ; )
		if (i->second != 2)
			Xm_pop.erase(i++);
		else
			i++;

	/* Report */
	LOG("[%d] COSI populations detected\n", Xm_pop.size());

	/* Summarize */
	char **Sa_ret = NULL;
	wsUint I = 0;
	wsAlloc(Sa_ret, char*, Xm_pop.size());
	FOREACHDO (mDataIdx_it, Xm_pop, i, I++)
		Sa_ret[I] = strdup(i->first.c_str());
	*Np_sz = I;

	return Sa_ret;
}

void loadCOSI(char *S_prefix)
{
	/* Find files starting with S_prefix . * */
	wsUint	N_files		= 0;
	char**	Sa_files	= findCOSI(S_prefix, &N_files);

	/* Check .pos integrity */
	char S_fn[MAX_PATH + 1];
	cStrFile *Ca_files = NULL;
	wsAlloc(Ca_files, cStrFile, N_files);
	for (wsUint i=0 ; i<N_files ; i++) {
		sprintf(S_fn, "%s.pos-%s", S_prefix, Sa_files[i]);
		Ca_files[i].open(S_fn, "COSI result POS");
	}

	/* Now read lines one by one... */
	bool	B_eof = false;
	char*	S_buf = NULL;
	wsAlloc(S_buf, char, 4096);
	for (wsUint L=1 ; !B_eof ; L++) {
		if (Ca_files[0].gets(S_buf, 4096) == NULL) B_eof = true;
		char *a, *b, *c;
		char S_curName[256];

		/* Get CHROM and CHROM_POS */
		getString(&S_buf, &a);
		getString(&a, &b);
		getString(&b, &c);
		strcpy(S_curName, a);
		strcat(S_curName, b);

		for (wsUint i=0 ; i<N_files ; i++) {
			if (Ca_files[i].gets(S_buf, 4096) == NULL) B_eof = true;
			char S_curName2[256];

			/* Get CHROM and CHROM_POS */
			getString(&S_buf, &a);
			getString(&a, &b);
			getString(&b, &c);
			strcpy(S_curName2, a);
			strcat(S_curName2, b);

			/* Same? */
			if (strcmp(S_curName, S_curName2))
				halt("Contents of POS file are different at line [%d]", L);
		}
	}

	/* Every file now should be EOF */
	for (wsUint i=0 ; i<N_files ; i++) {
		if (!Ca_files[i].isEnd())
			halt("File [%s.pos-%s] have different POS size with others!",
				S_prefix, Sa_files[i]);
		Ca_files[i].close();
	}

	/* Now read HAP file... */
	for (wsUint i=0 ; i<N_files ; i++) {
		sprintf(S_fn, "%s.hap-%s", S_prefix, Sa_files[i]);
		Ca_files[i].open(S_fn, "COSI result HAP");

		/* FIXME : Impl req. */
	}
}

vector<vInt> _buildPropagatePath(cIO *Cp_IO, wsUint **Np_t1, wsUint **Np_t2,
	mDataIdx *Xp_insIdx/*=NULL*/, char B_noMissing/*=1*/)
{
	/* Fetch required information from IO class */
	mSamp&	Xa_samp		= Cp_IO->getSampleData();
	mFam&	Xm_fam		= Cp_IO->getFamilyData();
	/* Index for missing sample */
	wsUint	N_misIdx	= Cp_IO->sizeSample();
	/* Number of 'entire' samples */
	int		N_aSamp		= (int)Xa_samp.size();
	/* IID -> propagate index mapping */
	vStr	Xv_iid;
	Xv_iid.resize(N_aSamp);

	/* FID, propagate index pair */
	map<string,vector<int> > Xm_famMems;

/**/char *Na_filled = NULL;
	wsCalloc(Na_filled, char, N_aSamp);
	FOREACH (mSamp_it, Xa_samp, j) {
		int N_idx = j->second.N_idx;
		/* If B_noMissing, skip missings */
		if (N_idx == SAMP_NODATA) {
			if (B_noMissing) continue;
			/* Fill index */
			mDataIdx_it X_pfind = Xp_insIdx->find(j->second.S_IID);
			/* If not, register */
			if (X_pfind == Xp_insIdx->end())
				Xp_insIdx->insert(make_pair(j->second.S_IID, N_misIdx++));
			N_idx = (*Xp_insIdx)[j->second.S_IID];
		} else if (Xp_insIdx)
			Xp_insIdx->insert(make_pair(j->second.S_IID, N_idx));
		/* Duplication / bug check */
		if (N_idx >= N_aSamp || Na_filled[N_idx] || N_idx < 0)
			halt("SYSERR : Duplicated N_idx[%d] found", N_idx);

		Xm_famMems[j->second.S_FID].push_back(N_idx);
		Na_filled[N_idx] = 1;
		Xv_iid[N_idx] = j->second.S_IID;
	}
	DEALLOC(Na_filled);

	vector<vInt> Xv_prop;
	/* Insert level 0 */
	Xv_prop.resize(1);

	wsAlloc(*Np_t1, wsUint, N_aSamp);
	wsAlloc(*Np_t2, wsUint, N_aSamp);
	FOREACH (mFam_it, Xm_fam, i) {
		xFamily &X_fam = i->second;

		X_fam.setUnvisited();

		/* Step 1 : Check OK for all founders */
		FOREACH (vSampPtr_it, X_fam.Xp_founders, ii) {
			int N_idx = (*ii)->N_idx;
			if ((*ii)->N_isVisited == 1) continue;
			if ((*ii)->N_idx == SAMP_NODATA) {
				if (B_noMissing) {
					if (OPT_ENABLED(usemf)) halt("SYSERR : Should not be NULL");
					continue;
				}
				/* Remap index */
				N_idx = (*Xp_insIdx)[(*ii)->S_IID];
			}

			(*ii)->N_isVisited = 1;
			Xv_prop[0].push_back(N_idx);
			pverbose("Step %d : %s [%d]\n", 0, (*ii)->S_IID.c_str(), N_idx);
		}
		 
		/* Step 2 : Mark if both parents are sane */
		vInt&	Xv_mems	= Xm_famMems[i->first];
		/* # of samples filled at this level */
		wsUint	N_fill	= 0;
		/* Current propagation level */
		short	N_level	= 2;
		do {
			N_fill = 0;
			FOREACH (vInt_it, Xv_mems, ii) {
				string	&ID		= Xv_iid[*ii];
				xSample	&X		= Xa_samp[Xv_iid[*ii]];
				char	B_miss	= X.B_isMissing || X.N_idx == SAMP_NODATA;
				xSample	*Xpp	= X.Xp_pat;
				xSample	*Xmm	= X.Xp_mat;

				/* Skip if already visited */
				if (X.N_isVisited) continue;
				if (B_miss && B_noMissing) {
					X.N_isVisited = N_level;
					continue;
				}

				/* No have parents (=founder) should not reach to here */
				if (!Xpp || !Xmm)
					halt("Sample [%s::%s] have non-full parent [%s father, "
						"%s mother], which is not permitted!",
						X.S_FID.c_str(), ID.c_str(),
						Xpp?"HAVE":"NO", Xmm?"HAVE":"NO");
				else if (Xpp->N_idx == SAMP_NODATA || Xmm->N_idx == SAMP_NODATA) {
					/* Halt if B_noMissing == 1 */
					if (B_noMissing)
						halt("At least one parent of sample [%s::%s] is missing, "
							"no missing founder is allowed in --simfam!",
							X.S_FID.c_str(), ID.c_str());
					B_miss = 1;
				/* Skip if its parent are in same level */
				}
				
				/* Skip if its parents are 'NOT' ready */
				if (Xpp->N_isVisited == N_level ||
					Xmm->N_isVisited == N_level ||
					Xpp->N_isVisited == 0 || Xmm->N_isVisited == 0)
					continue;

				/* Use the index 'as is' unless B_noMissing = 1 */
				wsUint N_insIdx = 0;
				if (Xp_insIdx) {
					mDataIdx_it X_find = Xp_insIdx->find(ID);
					if (X_find == Xp_insIdx->end())
						halt("SYSERR: Non-mapped sample [%s]", ID.c_str());
					N_insIdx = X_find->second;
				} else
					N_insIdx = X.N_idx;
				wsUint N_pIdx = Xpp->N_idx == SAMP_NODATA ?
					(*Xp_insIdx)[Xpp->S_IID] : Xpp->N_idx;
				wsUint N_mIdx = Xmm->N_idx == SAMP_NODATA ?
					(*Xp_insIdx)[Xmm->S_IID] : Xmm->N_idx;
				pverbose("Step %d : %s [%d,%d]\n", N_level, ID.c_str(),
					Xpp->N_isVisited, Xmm->N_isVisited);
				if ((short)Xv_prop.size() == (N_level-1))
					Xv_prop.push_back(vInt());
				Xv_prop[N_level-1].push_back(N_insIdx);

				/* Set parents' index */
				(*Np_t1)[N_insIdx] = N_pIdx;
				(*Np_t2)[N_insIdx] = N_mIdx;

				/* Marking */
				X.N_isVisited = N_level;
				N_fill++;
			}
			N_level++;
		} while (N_fill);

		/* Final step : Check orphan */
		FOREACH (vInt_it, Xv_mems, ii) {
			xSample	&X		= Xa_samp[Xv_iid[*ii]];
			if (!X.N_isVisited)
				halt("[%s::%s] is orphan, cannot continue simulation",
					X.S_FID.c_str(), X.S_IID.c_str());
		}
	}

	//if (Xp_insIdx) Xp_insIdx->clear();

	return Xv_prop;
}

cSimTrioIO::cSimTrioIO(int N_simFam, int N_simSNP)
{
	N_variant	= N_simSNP;
	N_sample	= N_simFam*3;
	N_sampOrig	= N_sample;

	_initMemories();

	/* Set to first nFam*2 samples to parent */
	/* Assign their phenotype randomly */
	Ba_typePheno[0] = 0; /* Dichotomous trait */
	for (int i=0 ; i<N_simFam*3 ; i++) {
		Ba_isFounder[i] = i<(N_simFam*2);
		Ra_pheno[0][i] = (wsReal)((rand()%2)+1);
	}
	for (int i=0 ; i<N_simFam*3 ; i++) {
		for (int j=0 ; j<N_simSNP ; j++)
			Na_geno[i][j] = rand()%3;
	}
	N_founder	= N_simFam*2;
}

void cSimFamIO::_procFAM(char *Sp_buf, wsUint i, wsUint &j)
{
	string	S_currIID;
	char	*Sp_tmp1 = NULL, *Sp_tmp2 = NULL;

	/* Get FID & IID */
	Sp_tmp2 = loadFIDIID(Sp_buf, S_currIID, S_currIID);

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
		return;
	}

	/* Is this sample founder? */
	if ((Ba_isFounder[j] = !X_sample.Xp_mat && !X_sample.Xp_pat) != 0)
		N_founder++;
	/* Allocate initial data index */
	X_sample.N_idx = j;
	Xa_sampleV2.push_back(&X_sample);

	/* a = phenotype, retrieve anyway, because alternative phenotype should
		* processed AFTER this procedure */
	//char *Sp_shouldNull = NULL;
	if (!OPT_ENABLED(nopheno)) {
		getString(&Sp_tmp1, &Sp_tmp2);
		char* Sp_shouldNull = NULL;
		Ra_pheno[0][j] = (wsReal)str2dbl(Sp_tmp1, &Sp_shouldNull);
		if (Sp_shouldNull && Sp_shouldNull[0])
			halt("Phenotype [%s] for sample [%s] looks invalid",
				Sp_tmp1, X_sample.S_FID.c_str());
		if (Ra_pheno[0][j] != OPT_NUMBER(phenoCtrl) &&
			Ra_pheno[0][j] != OPT_NUMBER(phenoCase))
			Ba_typePheno[0] = 1;
	} else Sp_tmp2 = Sp_tmp1;

	/* Status */
	if (((++j)%100) == 0) notice("%d samples retrieved\r", j);
}

cSimFamIO::cSimFamIO(wsStrCst S_fn, char B_inDry)
{
	B_dry = B_inDry;
	/* Can't dry */
	if (B_dry) halt("Simulation cannot be [[DRY]]");

	ASSERT_OPTION(szvar);
	N_variant	= OPT_NUMBER(szvar);

	/* Load .fam file */
	char *Sp_buf = NULL;
	wsAlloc(Sp_buf, char, 8192);
	cStrFile *Cp_fam = NULL;
	if (IS_ASSIGNED(fam)) {
		Cp_fam = new cStrFile(S_fn, "Pedigree base file");
		N_sampOrig = 0;
		for (wsUint L=1 ; Cp_fam->gets(Sp_buf, 8192) ; L++,N_sampOrig++)
			_loadRealTimeFamData(Sp_buf, L-1);
		Cp_fam->rewind();
	} else if (OPT_ENABLED(indep)) {
		wsStrCst S_bufs[1] = {
			"FAM	SAMP	0	0	1	-9",
		};
		for (wsUint i=0 ; i<1 ; i++) {
			char *Sp_buf = strdup(S_bufs[i]);
			_loadRealTimeFamData(Sp_buf, 0);
			free(Sp_buf);
		}
	} else if (OPT_ENABLED(trio)) {
		wsStrCst S_bufs[3] = {
			"FAM	SAMP1	0	0	1	-9",
			"FAM	SAMP2	0	0	2	-9",
			"FAM	SAMP3	SAMP1	SAMP2	1	-9",
		};
		for (wsUint i=0 ; i<3 ; i++) {
			char *Sp_buf = strdup(S_bufs[i]);
			_loadRealTimeFamData(Sp_buf, 0);
			free(Sp_buf);
		}
	} else if (OPT_ENABLED(extfam)) {
		wsStrCst S_bufs[10] = {
			"FAM	SAMP1	0	0	1	-9",
			"FAM	SAMP2	0	0	2	-9",
			"FAM	SAMP3	SAMP1	SAMP2	1	-9",
			"FAM	SAMP4	0	0	1	-9",
			"FAM	SAMP5	0	0	2	-9",
			"FAM	SAMP6	SAMP4	SAMP5	2	-9",
			"FAM	SAMP7	SAMP3	SAMP6	1	-9",
			"FAM	SAMP8	SAMP3	SAMP6	2	-9",
			"FAM	SAMP9	SAMP3	SAMP6	1	-9",
			"FAM	SAMP10	SAMP3	SAMP6	2	-9",
		};
		for (wsUint i=0 ; i<10 ; i++) {
			char *Sp_buf = strdup(S_bufs[i]);
			_loadRealTimeFamData(Sp_buf, i);
			free(Sp_buf);
		}
	}

	/* Make complete family data */
	_registMissingFounder();

	/* Init memories */
	_initMemories(N_variant > 0);

	/* Read FAM file : main */
	N_founder = 0;
	/* should be --1 but forget check */
	//wsUint N_0 = 0;
	Ba_typePheno[0] = 0;
	if (IS_ASSIGNED(fam)) for (wsUint i=0,j=0 ; i<N_sampOrig ; i++) {
		if (Cp_fam->gets(Sp_buf, 8192) == 0)
			halt("Unexpected end of file found at line %d", i+1);
		_procFAM(Sp_buf, i, j);
	} else if (OPT_ENABLED(indep)) {
		wsUint j = 0;
		wsStrCst S_bufs[1] = {
			"FAM	SAMP	0	0	1	-9",
		};
		for (wsUint i=0 ; i<1 ; i++) {
			char *Sp_buf = strdup(S_bufs[i]);
			_procFAM(Sp_buf, i, j);
			free(Sp_buf);
		}
	} else if (OPT_ENABLED(trio)) {
		wsUint j = 0;
		wsStrCst S_bufs[3] = {
			"FAM	SAMP1	0	0	1	-9",
			"FAM	SAMP2	0	0	2	-9",
			"FAM	SAMP3	SAMP1	SAMP2	1	-9",
		};
		for (wsUint i=0 ; i<3 ; i++) {
			char *Sp_buf = strdup(S_bufs[i]);
			_procFAM(Sp_buf, i, j);
			free(Sp_buf);
		}
	} else if (OPT_ENABLED(extfam)) {
		wsUint j = 0;
		wsStrCst S_bufs[10] = {
			"FAM	SAMP1	0	0	1	-9",
			"FAM	SAMP2	0	0	2	-9",
			"FAM	SAMP3	SAMP1	SAMP2	1	-9",
			"FAM	SAMP4	0	0	1	-9",
			"FAM	SAMP5	0	0	2	-9",
			"FAM	SAMP6	SAMP4	SAMP5	2	-9",
			"FAM	SAMP7	SAMP3	SAMP6	1	-9",
			"FAM	SAMP8	SAMP3	SAMP6	2	-9",
			"FAM	SAMP9	SAMP3	SAMP6	1	-9",
			"FAM	SAMP10	SAMP3	SAMP6	2	-9",
		};
		for (wsUint i=0 ; i<10 ; i++) {
			char *Sp_buf = strdup(S_bufs[i]);
			_procFAM(Sp_buf, i, j);
			free(Sp_buf);
		}
	}
//	for (int i=0 ; i<N_sample ; i++) LOG("%g\n", Ra_pheno[0][i]);
	if (Cp_fam) delete Cp_fam;
	DEALLOC(Sp_buf);

	/* Usepp */
	if (OPT_ENABLED(usemf)) {
		mFam&	Xm_fam	= getFamilyData();		

		wsUint	N_ins = 0;
		FOREACH (mFam_it, Xm_fam, i) {
			xFamily &X_fam = i->second;

			X_fam.setUnvisited();

			/* Step 1 : Check OK for all founders */
			vSampPtr x = X_fam.Xp_founders;
			FOREACH (vSampPtr_it, x, ii)
				if ((*ii)->N_idx == SAMP_NODATA) {
					(*ii)->N_idx = (int)Xa_sampleV2.size();
					(*ii)->B_isMissing = 0;
					Xa_sampleV2.push_back(*ii);
					N_ins++;
					N_founder++;
//					X_fam.Xp_founders.push_back(*ii);
				}
		}
		LOG("%d missing founders inserted newly\n", N_ins);

		/* Update sample size */
		N_sampOrig = N_sample = (wsUint)Xa_sampleV2.size();
	}

	/* Check szfounder */ {
		mFam&	Xm_fam	= getFamilyData();		
		wsUint N_nsamp = 0;
		FOREACH (mFam_it, Xm_fam, i)
			N_nsamp += (wsUint)(i->second.Xp_founders.size());

		if (N_nsamp != N_founder) {
			halt("Pedigree file [%s] contains missing founders, redo this with --usemf or --filmf");
//			halt("SYSERR : Founder number[%d] does not match with true [%d]",
//			N_founder, N_nsamp);
		}
	}

	/* szfam process */
	if (IS_ASSIGNED(szfam)) {
		wsUint N_szBloat = N_sample;
		N_sample	*= OPT_NUMBER(szfam);
		N_sampOrig	*= OPT_NUMBER(szfam);
		N_founder	*= OPT_NUMBER(szfam);

		/* Copy sample */
		vSampPtr Xa_samp = Xa_sampleV2;

		/* Copy Xm_iid2data */
		mSamp Xm_iidcopy;
// 		copy(Xm_iid2data.begin(), Xm_iid2data.end(), inserter(Xm_iidcopy,
// 			Xm_iidcopy.begin()));

		/* Bloating Xa_sampleV2 and Xm_iid2data */
		//Xa_sampleV2.clear();
		//Xm_iid2data.clear();

		Xa_sampleV2.resize(N_sample);
		FOREACH (mSamp_it, Xm_iid2data, i) {

			wsUint N_idx = i->second.N_idx;
			wsUint N_fam = OPT_NUMBER(szfam);
			for (wsUint j=0 ; j<N_fam ; j++) {
				char Sp_curFam[64];
				sprintf(Sp_curFam, "_%d", j+1);

				xSample X_samp = i->second;
				X_samp.S_FID	+= Sp_curFam;
				X_samp.S_IID	+= Sp_curFam;
				X_samp.N_idx	= i->second.N_idx == SAMP_NODATA ? SAMP_NODATA : N_idx;
				Xm_iidcopy[X_samp.S_IID] = X_samp;

				N_idx += N_szBloat;
			}
		}
		/* Update pat & mat */
		map<string,string> Xp_pat;
		map<string,string> Xp_mat;
		mvStr	Xp_childs;
		mvStr	Xp_spouses;
		wsUint	N_fam = OPT_NUMBER(szfam);
		FOREACH (mSamp_it, Xm_iid2data, i) {
			for (wsUint j=0 ; j<N_fam ; j++) {
				char Sp_curFam[64];
				sprintf(Sp_curFam, "_%d", j+1);

				xSample &X_samp = i->second;
				string X_key = X_samp.S_IID + Sp_curFam;
				xSample &X_smp = Xm_iidcopy[X_key];
				X_smp.Xp_childs.clear(); 
				X_smp.Xp_spouses.clear();

				if (X_smp.Xp_mat) Xp_mat[X_key] = X_samp.Xp_mat->S_IID + Sp_curFam;
				if (X_smp.Xp_pat) Xp_pat[X_key] = X_samp.Xp_pat->S_IID + Sp_curFam;

				FOREACH (vSampPtr_it, X_samp.Xp_spouses, k)
					Xp_spouses[X_key].push_back((*k)->S_IID + Sp_curFam);
				FOREACH (vSampPtr_it, X_samp.Xp_childs, k)
					Xp_childs[X_key].push_back((*k)->S_IID + Sp_curFam);
			}
		}

		/* Update Xm_fam */
		mvStr	Xa_founders;
		mvStr	Xa_members;
		mFam&	Xm_fam = getFamilyData();
		mFam	Xm_fcopy;
		FOREACH (mFam_it, Xm_fam, i) {
			Xm_fcopy.insert(make_pair(i->first, i->second));

			Xa_founders[i->first] = vStr();
			FOREACH (vSampPtr_it, i->second.Xp_founders, k)
				Xa_founders[i->first].push_back((*k)->S_IID);
			Xa_members[i->first] = vStr();
			FOREACH (vInt_it, i->second.Xv_members, k)
				Xa_founders[i->first].push_back(Xa_sampleV2[*k]->S_IID);
		}
		Xm_fam.clear();

		/* Replace Xm_iid2data */
		Xm_iid2data.clear();
		FOREACH (mSamp_it, Xm_iidcopy, i) {
			Xm_iid2data.insert(make_pair(i->first, i->second));

			xSample &V = Xm_iid2data[i->first];
			if (i->second.N_idx != SAMP_NODATA)
				Xa_sampleV2[i->second.N_idx] = &V;
		}
		FOREACH (mSamp_it, Xm_iid2data, i) {
			xSample &V = i->second;

			/* Update pat/mat/spouses/childs */
			V.Xp_pat = Xp_pat[i->first].length() ? &(Xm_iid2data[Xp_pat[i->first]]) : NULL;
			V.Xp_mat = Xp_mat[i->first].length() ? &(Xm_iid2data[Xp_mat[i->first]]) : NULL;
			FOREACH (vStr_it, Xp_spouses[i->first], j) V.Xp_spouses.push_back(&(Xm_iid2data[*j]));
			FOREACH (vStr_it, Xp_childs[i->first], j) V.Xp_childs.push_back(&(Xm_iid2data[*j]));
		}

		FOREACH (mFam_it, Xm_fcopy, i) for (wsUint j=0 ; j<N_fam ; j++) {
			xFamily	&X_ofam = i->second;
			xFamily	X_nfam = i->second;
			char	Sp_curFam[64];
			sprintf(Sp_curFam, "_%d", j+1);

			X_nfam.S_FID += Sp_curFam;

			/* Fill Xp_founders */
			X_nfam.Xp_founders.clear();
			FOREACH (vStr_it, Xa_founders[X_ofam.S_FID], k) {
				string S_nsamp = *k + Sp_curFam;
				X_nfam.Xp_founders.push_back(&(Xm_iid2data[S_nsamp]));
			}
			/* Fill Xp_members */
			X_nfam.Xv_members.clear();
			FOREACH (vStr_it, Xa_members[X_ofam.S_FID], k) {
				string S_nsamp = *k + Sp_curFam;
				X_nfam.Xv_members.push_back(Xm_iid2data[S_nsamp].N_idx);
			}
			Xm_fam.insert(make_pair(X_nfam.S_FID, X_nfam));
		}

//  		copy(Xm_iidcopy.begin(), Xm_iidcopy.end(), inserter(Xm_iid2data,
//  			Xm_iid2data.begin()));

		/* if --verbose */
		if (OPT_ENABLED(verbose)) {
			mFam		&Xa_f	= getFamilyData();
			vSampPtr	&Xa_p	= getSample();

			FOREACH (mFam_it, Xa_f, i) {
				LOG("Family [%s]\n", i->first.c_str());
				FOREACH (vSampPtr_it, i->second.Xp_founders, j)
					LOG("Founder [%s::%s] idx %d\n", (*j)->S_FID.c_str(), (*j)->S_IID.c_str(), (*j)->N_idx);
				FOREACH (vInt_it, i->second.Xv_members, j)
					LOG("Member [%s::%s] idx %d\n", Xa_p[*j]->S_FID.c_str(), Xa_p[*j]->S_IID.c_str(), Xa_p[*j]->N_idx);
			}
		}
	}

	/* Compute propagate path from given pedigree */
	Xv_prop	= _buildPropagatePath(this, &Na_t1, &Na_t2);
	vSampPtr&		Xv_samp	= getSample();
	mFam&			Xm_fam	= getFamilyData();

	/* Make report if required */
	if (OPT_ENABLED(verbose)) {
		cExporter *Cp_o = cExporter::summon("simfam.prop.seq");
		wsUint L=1;
		FOREACHDO (vector<vInt>::iterator, Xv_prop, i, L++) {
			Cp_o->fmt("Level[%d]", L);
			FOREACH (vInt_it, *i, ii)
				Cp_o->fmt("	%s", Xv_samp[*ii]->S_IID.c_str());
			Cp_o->put("\n");
		}
		delete Cp_o;
	}

	char *tmp1 = new char[1];
	tmp1[0] = Ba_typePheno[0];
	wsVec tmp2 = new wsReal[N_sample];
	memcpy(tmp2, Ra_pheno[0], sizeof(wsReal)*N_sample);
	_initMemories();
	memcpy(Ra_pheno[0], tmp2, sizeof(wsReal)*N_sample);
	Ba_typePheno[0] = tmp1[0];
	delete [] tmp1;
	delete [] tmp2;
	/* Ba_isFounder setting */
	FOREACH (mFam_it, Xm_fam, i) {
		xFamily &X_fam = i->second;

		/* Step 1 : Check OK for all founders */
		FOREACH (vSampPtr_it, X_fam.Xp_founders, ii)
			if ((*ii)->N_idx != SAMP_NODATA)
				Ba_isFounder[(*ii)->N_idx] = 1;
	}

	/* Load --mafsnp if required */
	wsUint	N_nAllele	= N_sample - 1;
	wsUint	N_denom		= N_sample << 1;
	/* Load frequency if required */	
	if (IS_ASSIGNED(simfreq)) {
		cStrFile	C_sf(OPT_STRING(simfreq), "Frequency background file for MAF values");
		vReal		Xv_freq;
		wsReal		R_minMAF	= W1 / (wsReal)(N_sample << 1);
		char		*S_buf		= NULL;
		wsAlloc(S_buf, char, 1024);

		wsUint	N_fail = 0;
		for (wsUint L=1 ; C_sf.gets(S_buf, 1024) ; L++) {
			char *a = S_buf, *b = NULL;

			/* Seek final column */
			getString(&a, &b);
			while (b) {
				a = b;
				getString(&a, &b);
			}

			wsReal V = (wsReal)atof(a);
			/* Since it is population allele freq. */
			V = V>REAL_CONST(0.5) ? W1 - V : V;

			/* ...and should be over the minimum freq. */
			if (R_minMAF < V)
				Xv_freq.push_back(V);
			else
				N_fail++;
		}
		/* Randomly assign MAF */
		wsUint N_sz = (wsUint)Xv_freq.size();
		if (N_fail > 0)
			LOGwarn("[%d] MAFs dropped due to the lowest possible MAF [%g]\n",
				N_fail, R_minMAF);
		if (N_sz == 0)
			halt("No available MAFs found in [%s]", OPT_STRING(simfreq));
				
		for (wsUint i=0 ; i<N_variant ; i++)
			Xv_maf.push_back(Xv_freq[rand()%N_sz]);
		LOG("[%d] MAFs chosen from [%d] background MAFs\n", N_variant, Xv_freq.size());
		DEALLOC(S_buf);
	} else if (IS_ASSIGNED(mafvar)) {
		cStrFile	C_ms(OPT_STRING(mafvar), "MAF values file for SNVs simulated", 1);
		if (C_ms.isFailed()) {
			/* Interpret as range */
			xOption X;
			X.S_longName = (char *)"--mafvar";
			OPTION().procRange(&X, OPT_STRING(mafvar));
			xOptRange& R = X.X_rngVal;

			if (NA(R.R_e)) R.R_e = REAL_CONST(0.5);
			else if (!R.R_eEQ) R.R_e -= W1 / (wsReal)N_denom;
			if (NA(R.R_s)) R.R_s = W0;
			else if (!R.R_sEQ) R.R_s += W1 / (wsReal)N_denom;

			if (R.R_e > REAL_CONST(0.5) || R.R_e <= W0 ||
				R.R_s > REAL_CONST(0.5) || R.R_s <= W0)
				halt("Invalid MAF range from --mafvar!");

			/* Initialize MAF */
			for (wsUint i=0 ; i<N_variant ; i++) {
				wsReal R_maf = randInRange(R.R_s, R.R_e);
				Xv_maf.push_back(R_maf);
//				LOG("%g\n", R_maf);
			}
			LOGnote("Simulated MAFs will ranged within [%g,%g]\n",
				R.R_s, R.R_e);
		} else {
			char		*S_buf = NULL;
			wsAlloc(S_buf, char, 1024);
			for (wsUint L=0 ; C_ms.gets(S_buf, 1024) ; L++) {
				if (!S_buf[0]) break;

				wsReal V = (wsReal)atof(S_buf);
				/* Should be placed btwn (0.0,0.5] */
				if (V != V || V <= W0 || V > REAL_CONST(0.5))
					halt("Invalid MAF range [%g] assigned at line [%d] of file [%s]",
						V, L+1, OPT_STRING(mafvar));

				/* Insert */
				Xv_maf.push_back((wsReal)atof(S_buf));
			}
			/* Size check */
			if (Xv_maf.size() != N_variant)
				halt("Number of MAFs [%d] in file [%s] does not match with --szvar [%d]",
					Xv_maf.size(), OPT_STRING(mafvar), OPT_NUMBER(szvar));
			DEALLOC(S_buf);
		}
	} else {
		/* Otherwise initialize MAF */
		for (wsUint i=0 ; i<N_variant ; i++)
			Xv_maf.push_back((wsReal)(rand()%N_nAllele + 1) / (wsReal)N_denom);
	}

	/* Do simulation (0, 0.5] */
	//wsReal R_minMAF = W1 / (wsReal)(N_sample<<1);
	wsUint			N_sig		= IS_ASSIGNED(nsig) ? OPT_NUMBER(nsig) : 0;
	wsReal			R_maf		= WISARD_NAN;
	wsReal			R_heri		= WISARD_NAN;
	wsReal			R_beta		= WISARD_NAN;
	if (N_sig) {
		ASSERT_OPTION(sigmaf);
		ASSERT_OPTION(heri);
		ASSERT_OPTION(randpheno);

		R_maf		= OPT_REAL(sigmaf);
		R_heri		= atof(OPT_STRING(heri));
	}
	cPDDTAnalysis*	Cp_anaPDDT = NULL;
	Mp_V		= NULL;
	cVector			V_y(Ra_pheno[0], N_sample, 1);
	//SYM_t			Ra_cor		= NULL;
	if (N_sig) {
		Cp_anaPDDT = new cPDDTAnalysis(this);
		Cp_anaPDDT->run();

		/* Compute sig^2(g) */
		wsReal	R_p1p	= R_maf * (W1 - R_maf);
		if (IS_ASSIGNED(beta))
			R_beta = OPT_REAL(beta);
		else
			R_beta	= sqrt(0.005 / (REAL_CONST(1.99)*R_p1p));
		LOG("Desired beta is [%g]\n", R_beta);
		wsReal	R_2bp1p	= W2 * SQR(R_beta) * R_p1p;
		wsReal	R_sigg	= R_heri - R_2bp1p;
		wsReal	R_sige	= W1 - R_sigg - R_2bp1p;

		/* Compute V^-1 */
		cSymMatrix M_cor(Cp_anaPDDT->getPDDT(), N_sample);
		M_cor *= R_sigg;
		for (wsUint i=0 ; i<N_sample ; i++)
			M_cor.get()[i][i] += R_sige;
		Mp_V = &(M_cor.inv());
	}
	LOG("[%d] variants will be simulated with [%d] significant variant(s)\n",
		N_variant, N_sig);
	for (wsUint i=0 ; i<N_variant ; i++) {
		xVariant v = { 0, };
		v.al1		= '1';
		v.al2		= '2';
		v.chr		= 1;
		v.nmis		= 0;
		v.filter	= 0;
		v.gdist		= 0;
		v.hwe		= WISARD_NAN;
		v.indel1	= NULL;
		v.indel2	= NULL;
		wsAlloc(v.name, char, 32);
		sprintf(v.name, "SNP%d", i+1);
		v.pos = i+1;
		Xv_variant.push_back(v);
	}
	pverbose("Variant info made\n");

	/* Prepare memory space */
	wsAlloc(Xp_maf, xMaf, N_variant);

	/* Now, simulate markers */
	if (OPT_NUMBER(thread) == 1) for (wsUint i=0 ; i<N_variant ; i++) {
		if (i < N_sig) while (!_simsig(Xv_prop, i, R_maf, Mp_V, V_y, R_beta));
		else while (!_sim(Xv_prop, i, Xv_maf[i]));
		if ((i%1000) == 0)
			notice("[%d/%d] variants simulated...\r", i, N_variant);
	} else WORKER().run(genMarker, forAllVariant_equal, this, NULL, sizeof(int)*3);

	/* MAF is done */
	B_mafComputed = 1;

	Bv_filtVrt.resize(N_variant);
	Bv_filtSamp.resize(N_sample);
	LOG("Simulation data generated, [%d] samples with [%d] variants\n",
		N_sample, N_variant);
	N_sampOrig = N_sample;

	/* Apply missing rate */
	if (IS_ASSIGNED(nageno)) {
		/* Export original */
		LOG("Export complete simulated data to [%s.sim.ori.res]\n",
			OPT_STRING(out));
		exportMatrix("sim.ori.res", Na_geno, N_sample, N_variant);
		
		wsReal	R_geno	= ((wsReal)N_variant * (wsReal)N_sample);
		INT64_t N_mis	= 0;
		wsReal	R_mis	= OPT_REAL(nageno);
		if (R_mis > W0 && R_mis < W1)			
			N_mis = (INT64_t)(R_geno * R_mis);
		else if (R_mis >= W1 && (wsReal)(INT64_t)R_mis == R_mis) {
			/* Range check */
			if (R_mis > R_geno) {
				LOGwarn("Too large number of missing[%g] assigned, adjusted to [%g]\n",
					R_mis, R_geno);
				R_mis = R_geno;
			}
			N_mis = (INT64_t)R_mis;
		} else
			halt("Invalid value of --nageno [%g]", R_mis);
		LOG("[" FMT_INT64 "] genotypes are marked as missing\n", N_mis);
		for (INT64_t i=0 ; i<N_mis ; i++) {
			wsUint N_x = (wsUint)(wsRand()%N_variant);
			wsUint N_y = (wsUint)(wsRand()%N_sample);

			if (isMissing(Na_geno[N_y][N_x])) {
				i--;
				continue;
			}
			N_missGeno++;
			Na_geno[N_y][N_x] = WISARD_NA;
			if ((i%1000) == 0)
				notice("[" FMT_INT64 "] genotypes remained...\r", i);
		}
	}
}

cSimFamIO::~cSimFamIO()
{
	Xv_maf.clear();
	DEALLOC(Na_t1);
	DEALLOC(Na_t2);
}

inline char doMendelianTransmission(char v1, char v2)
{
	int V = v1*3 + v2;
	
	switch (V) {
	/* 0 0 */ case 0: return 0;
	/* 0 1, 1 0 */ case 1: case 3:
		return (char)genBinomial(1.0, 0.5);
	/* 0 2, 2 0 */ case 2: case 6: return 1;
	/* 1 1 */ case 4:
		return (char)genBinomial(2.0, 0.5);
	/* 1 2, 2 1 */ case 5: case 7:
		return (char)genBinomial(1.0, 0.5) + 1;
	/* 2 2 */ case 8: return 2;
	}

	/* Halt because it cannot be occured */
	halt("SYSERR : Errornous transmission value [%d]", V);
	return 0;
}

int cSimFamIO::_sim(vector<vInt> &Xv_prop, wsUint i, wsReal R_maf)
{
	xMaf&		X_maf	= Xp_maf[i];
	char**		Na_data	= getGenotype();
	wsUintCst		N_samp	= sizeSample();
	wsUintCst		N_fnd	= sizeFounder();
	wsUint		S		= 0;
	wsUint		N_try	= 0;

	do {
		if (++N_try == 1000)
			halt("SYSERR : Too much tries in --simfam, maybe error?");
		/* 1st pass */
		FOREACH (vInt_it, Xv_prop[0], j) {
			int x = *j;
			char V = (char)genBinomial(2, R_maf);
			S += Na_data[x][i] = V;
		}
		X_maf.R_maf = S / (wsReal)(N_fnd<<1);
		if (X_maf.R_maf > 0.5) return 0;
		X_maf.N_mac = S;

		/* impute descendants */
		wsUint N_lv = 0;
		FOREACHDO (vector<vInt>::iterator, Xv_prop, j, N_lv++) {
			if (N_lv == 0) continue;

			FOREACH (vInt_it, *j, k) {
				int x = *k;
				wsUint N_i1 = Na_t1[x];
				wsUint N_i2 = Na_t2[x];

				if (Xa_sampleV2[x]->N_isVisited != (short)(N_lv+1))
					halt("SYSERR: Sample [%s] leveled as [%d] but visited at level [%d]",
					Xa_sampleV2[x]->S_IID.c_str(), Xa_sampleV2[x]->N_isVisited,
					N_lv+1);

				char V = doMendelianTransmission(Na_data[N_i1][i], Na_data[N_i2][i]);
				S += Na_data[x][i] = V;
			}
		}
		X_maf.N_allMac	= S;
		X_maf.N_allele	= N_samp << 1;
		X_maf.R_allMaf	= S / (wsReal)(N_samp<<1);
	} while (S == 0); /* Avoid monomorphic SNP */

	return 1;
}

int cSimFamIO::_simsig(vector<vInt> &Xv_prop, wsUint i, wsReal R_maf,
	cSymMatrix *Mp_V, cVector &V_y, wsReal R_beta)
{
	xMaf&	X_maf	= Xp_maf[i];
	char**	Na_data	= getGenotype();
	wsUintCst	N_samp	= sizeSample();
	wsUintCst	N_fnd	= sizeFounder();
	wsUint	S		= 0;
	wsUint	N_try	= 0;
	wsReal	R_bHat	= WISARD_NAN;
	bool	B_ok	= false;

_rerun:
	do {
		S = 0;
		R_bHat	= WISARD_NAN;
		cVector	V_x(WISARD_NAN, N_samp);
		wsReal*	Ra_x	= V_x.get();
		if (++N_try == 1000) {
			if (!OPT_ENABLED(nostop))
				halt("SYSERR : Too much tries in --simfam, maybe too strict beta[%g]?", R_beta);
			else break;
		}
		/* 1st pass */
		FOREACH (vInt_it, Xv_prop[0], j) {
			int x = *j;
			char V = (char)genBinomial(2, R_maf);
			S += Na_data[x][i] = V;
			Ra_x[x] = (wsReal)V;
		}
		X_maf.R_maf	= S / (wsReal)(N_fnd<<1);
		X_maf.N_mac	= S;
		if (X_maf.R_maf > 0.5) return 0;

		/* impute descendants */
		wsUint N_lv = 0;
		FOREACHDO (vector<vInt>::iterator, Xv_prop, j, N_lv++) {
			if (N_lv == 0) continue;

			FOREACH (vInt_it, *j, k) {
				int x = *k;
				wsUint N_i1 = Na_t1[x];
				wsUint N_i2 = Na_t2[x];

				if (Xa_sampleV2[x]->N_isVisited != (short)(N_lv+1))
					halt("SYSERR: Sample [%s] leveled as [%d] but visited at level [%d]",
					Xa_sampleV2[x]->S_IID.c_str(), Xa_sampleV2[x]->N_isVisited,
					N_lv+1);

				char V = doMendelianTransmission(Na_data[N_i1][i], Na_data[N_i2][i]);
				S += Na_data[x][i] = V;
				if (Ra_x[x] == Ra_x[x])
					halt("Duplicated write");
				Ra_x[x] = (wsReal)V;
			}
		}
		X_maf.N_allMac	= S;
		X_maf.N_allele	= N_samp<<1;
		X_maf.R_allMaf	= S / (wsReal)(N_samp<<1);

		/* Get beta estimate */
		if (S) {
			cVector	V_xVinv	= V_x * *Mp_V;
			R_bHat	= V_xVinv.sum(V_y) / V_xVinv.sum(V_x);
			if (R_bHat != R_bHat)
				halt("Error bHat");
		}
		B_ok = fabs(R_bHat) >= fabs(R_beta);
	} while (S == 0 || !B_ok); /* Avoid monomorphic SNP */

	/* Re-generate phenotype */
	if (!B_ok && !IS_ASSIGNED(fam)) {
		wsReal *Ra_y = V_y.get();
		for (wsUint i=0 ; i<V_y.size() ; i++) {
			wsReal R_p = (wsReal)((rand()%(RAND_MAX-2)+1) / (wsReal)(RAND_MAX));
			if (R_p != R_p) halt("Invalid pheno");
			Ra_y[i] = qnorm(R_p);
		}
		pverbose("Retry....\n");
		N_try = 0;
		goto _rerun;
	}

	lverbose("Variant generated with beta[%g > %g] with [%d] tries\n", R_bHat, R_beta, N_try);

	return 1;
}

} // End namespace ONETOOL
