#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <vector>
#include <map>
#include <string>
#include "global/common.h"
#include "global/worker.h"
#include "global/option.h"
#include "utils/util.h"
#include "analyses/annot.h"
#include "global/io.h"
#include "input/stream.h"
#include "input/impute.h"
#include "input/sampvar.h"
#include "utils/dcdflib.h"
#include "analyses/emai.h"
#include "analyses/pddt.h"
#include "analyses/corr.h"
#include "analyses/ppp.h"
#include "analyses/fam.h"
#include "analyses/mendel.h"
#include "analyses/setmgr.h"
#include "utils/vector.h"
#include "utils/ftp.h"
#include "utils/data.h"
#include "utils/stat.h"
#include "utils/vis.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <iterator>
#include <list>
#include <vector>

namespace ONETOOL {

const wsUint MAX_TWINS						= 128;
const wsUint WINSARD_WARNSIZE_EMPCORCALC	= 10000;
const wsUint PED_MAXLEN						= 1024*1024*64;

#define IS_REG_COVAR(it) ((it)->R_times == W1 && (it)->N_idxMul == -1)

double variantHWE(int obs_hets, int obs_hom1, int obs_hom2);

wsReal _inbreed(wsUint N_idx, vVariant &Xv_var, xMaf* Xp_maf, char **Na_data,
	wsReal &O, wsReal &E, wsReal &N)
{
	// P(Homo) = F + (1-F)P(Homo by chance)
	wsUint N_vrt = (wsUint)Xv_var.size();

	// P(Homo by chance) = p^2+q^2 for a biallelic locus.
	// For an individual with N genotyped loci, we
	//   1. count the total observed number of loci which are homozygous (O),
	//   2. calculate the total expected number of loci homozygous by chance (E)
	// Then, using the method of moments, we have
	//    O = NF + (1-F)E
	// Which rearranges to give
	//    F = (O-E)/(N-E)


	// Count of nonmissing loci
	N=0;
	O=0;
	E=0;

	// Consider all variants
	char *Np_data = Na_data[N_idx];
	for (wsUint i=0 ; i<N_vrt ; i++) {
		xVariant&	X_mk	= Xv_var[i];
		xMaf&		X_maf	= Xp_maf[i];

		////////////////////////////////////////////////
		// Skip X and haploid chromosome variants, or not
		if (OPT_ENABLED(x)) { // For sex-checks
			// only consider the X chromosome
			if (!isSexChromosome(X_mk)) continue;
		} else { // normal heterozygosity calculation
			// so skip haploid variants
			if (!isAutosome(X_mk)) continue;
		}

		// Skip monomorphic variants, uninformative variants
// 		if (X_mk.allmaf == W0 || X_mk.mac < 2)
// 			continue;
		if (X_maf.N_allele == 0) continue;

		//////////////////
		// Observed data

		// check not missing:
		if (isMissing(Np_data[i])) continue;

		// homozygous non-missing loci
		if (Np_data[i] == 0 || Np_data[i] == 2) O++;
		// non-missing loci
		N++;


		/////////////////////////
		// Expected homozygousity

		// E = 2pq . 2N/(2N-1)
		// (Using Nei's unbiased estimator)

		E += W1
			- ( 2 * X_maf.R_maf * ( W1 - X_maf.R_maf ) 
			* ( (X_maf.N_allele>>1) / ( (X_maf.N_allele>>1) - 1 ) ) );
	}
	if ((N-E) == 0)
		return WISARD_NAN;
	return (O-E)/(N-E);

	// 	if ( par::check_sex)
	// 	{
	// 		HET << setw(par::pp_maxfid) << p1->fid << " " 
	// 			<< setw(par::pp_maxiid) << p1->iid << " "      
	// 			<< setw(12) << p1->sexcode << " ";
	// 
	// 		if ( F > par::sex_threshold_male ) 
	// 		{
	// 			HET << setw(12) << 1 << " ";
	// 			if (p1->sexcode == "1")
	// 				HET << setw(12) << "OK" << " ";
	// 			else
	// 			{
	// 				HET << setw(12) << "PROBLEM" << " ";
	// 				if (par::impute_sex)
	// 				{
	// 					p1->sexcode = "1";
	// 					p1->sex = true;
	// 				}
	// 			}
	// 		}
	// 		else if ( F < par::sex_threshold_female )
	// 		{
	// 			HET << setw(12) << 2 << " ";
	// 			if (p1->sexcode == "2")
	// 				HET << setw(12) << "OK" << " ";
	// 			else
	// 			{
	// 				HET << setw(12) << "PROBLEM" << " ";		  
	// 				if (par::impute_sex)
	// 				{
	// 					p1->sexcode = "2";
	// 					p1->sex = false;
	// 				}
	// 			}
	// 		}
	// 		else 
	// 		{
	// 			HET << setw(12) << 0 << " "
	// 				<< setw(12) << "PROBLEM" << " ";		  
	// 
	// 			if (par::impute_sex)
	// 			{
	// 				p1->sexcode = "0";
	// 				p1->sex = false;
	// 				if (!par::ignore_missing_sex)
	// 					p1->missing = true;	      
	// 			}
	// 
	// 		}
	// 
	// 
	// 		HET << setw(12) << F << "\n";
	// 	}
	// 	else
	// 		HET << setw(par::pp_maxfid) << p1->fid << " " 
	// 		<< setw(par::pp_maxiid) << p1->iid << " "
	// 		<< setw(12) << (int)O << " " 
	// 		<< setw(12) << E << " "
	// 		<< setw(12) << (int)N << " "
	// 		<< setw(12) << F << "\n";

//	return F;
}

vStr bloating(char *Sp_fn)
{
	vStr X_ret;

	/* Find asterisk(*) */
	char B_bloatable = 0;
	char *Sp_ast = NULL;
	for (char *p=Sp_fn ; *p ; p++) if (*p == '*') {
		if (Sp_ast == NULL) Sp_ast = p;
		B_bloatable++;
	}

	/* # of asterisk must be 1 */
	if (B_bloatable > 1)
		halt("Too many(%d) asterisk(*)s found from an input, please set to 1", B_bloatable);
	else if (B_bloatable == 0) {
		/* If there is no asterisk, just return the filename only */
		X_ret.push_back((string)Sp_fn);
		return X_ret;
	}

	/* Otherwise, try autosome first */
	*(Sp_ast++) = '\0';
	LOG("Multiloading try with [%s <CHROM> %s]\n", Sp_fn, Sp_ast);
	for (wsUint i=1 ; i<=NAUTO_SPECIES ; i++) {
		char	S_fn[512];
		char	B_found = 0;

		/* Make filename */
		sprintf(S_fn, "%s%d%s", Sp_fn, i, Sp_ast);

		/* Try to open */
		FILE *H_fp = fopen(S_fn, "rb");
		if (H_fp != NULL) {
			/* If can open, insert this to vector */
			X_ret.push_back((string)S_fn);
			B_found = 1;
			fclose(H_fp);
		}

		/* If failed, try two-digits with 0-filled */
		if (B_found == 0) {
			sprintf(S_fn, "%s%02d%s", Sp_fn, i, Sp_ast);
			H_fp = fopen(S_fn, "rb");
			if (H_fp != NULL) {
				/* If can open, insert this to vector */
				X_ret.push_back((string)S_fn);
				B_found = 1;
				fclose(H_fp);
			}
		}
		LOG("	Chromosome [%d] %s....\n", i, B_found?"found":"not found");
	}
	LOG("[%d] files found for multiloading\n", X_ret.size());
	return X_ret;
}

wsReal**  getFullCorMat(cAnalysis *Cp_anaPhi)
{
	wsReal** Ra_origPhi = NULL;

	if ((OPT_ENABLED(kinship) && !IS_ASSIGNED(cor)) || OPT_ENABLED(hybrid))
		Ra_origPhi = ((cPDDTAnalysis *)Cp_anaPhi)->getPDDT();
	else
		Ra_origPhi = ((cCorrAnalysis *)Cp_anaPhi)->getCorr();

	return Ra_origPhi;
}

int forAllinput(int N_mode, int N_thread, void *Vp_data,
	void *Vp_shareData, wsUint *Na_waitQueue)
{
	static wsUint N_idxFile;
	int		*Np_idx		= (int *)Vp_data;
	vStr	&Xa_files	= *((vStr *)Vp_shareData);
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		N_idxFile = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++) {
			if ((N_idxFile+i) >= Xa_files.size())
				break;
			Np_idx[i] = N_idxFile+i;
		}
		N_idxFile += N_thread;
		break;
	case DISTFUN_UNINIT:
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}
	return i;
}

int thr_loadFile(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
// 	vStr	&Xa_files	= *((vStr *)Vp_shareData);
// 	int		i			= *((int *)Vp_data);
// 	cIO		*Cp_io		= NULL;
// 	char	S_fn[512];
// 
// 	strcpy(S_fn, Xa_files[i].c_str());
// 
// 	if (IS_ASSIGNED(ped)) {
// 		Cp_io = cIO::summon(FT_PED, S_fn);
// 	} else if (IS_ASSIGNED(bed)) {
// 		Cp_io = cIO::summon(FT_BED, S_fn);
// 	} else if (IS_ASSIGNED(tped)) {
// 		Cp_io = cIO::summon(FT_TPED, S_fn);
// 	} else if (IS_ASSIGNED(vcf)) {
// 		Cp_io = cIO::summon(FT_VCF, S_fn);
// 	} else if (IS_ASSIGNED(dosage)) {
// 		Cp_io = cIO::summon(FT_DOSAGE, S_fn);
// 	}

	/* FIXME : Do more */
	return 0;
}

/* Running only once
 * Divide variants ALMOST equally to each thread */
int forAllVariant_equal(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx	= (int *)Vp_thrData;
	cIO		*Cp_IO	= (cIO *)Vp_shareData;
	wsUintCst	N_vrt	= Cp_IO->sizeVariant();
	static wsUint
			N_proc	= 0;
	wsReal	R_szDiv	= (wsReal)N_vrt/(wsReal)N_thread;
	int		i=0;
	wsReal	j=W0;
	int		N_ret = N_proc != 0 ? 0 : 1;

	switch (N_mode) {
	case DISTFUN_INIT:
		i = 0;
		j=W0;
		N_proc = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		if (N_ret == 0)
			return 0;
		for ( ; i<N_thread ; i++,j+=R_szDiv) {
			wsUint E = (wsUint)(j+R_szDiv+REAL_CONST(0.5));
			if (E > N_vrt)
				E = N_vrt;
			Np_idx[3*i] = (wsUint)(j+REAL_CONST(0.5));
			Np_idx[3*i+1] = E;
			Np_idx[3*i+2] = 0;
		}
		N_proc = N_ret = N_thread;
		break;
	case DISTFUN_UNINIT:
		LOG("%d/%d variants processed\n", N_vrt, N_vrt);
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}

	return N_ret;
}

int mafCalc(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	cIO*		Cp_IO		= (cIO *)Vp_shareData;
	wsUint		N_sample	= Cp_IO->sizeSample();
	char**		Na_data		= Cp_IO->getGenotype();
	wsUint		N_founder	= Cp_IO->sizeFounder();
	xMaf*		Xp_maf		= Cp_IO->getMAFbuffer();
	const char*	Ba_founder	= Cp_IO->getIsFounder();
	mStrReal&	Xm_newMAF	= Cp_IO->getNewMAF();
	vVariant&	Xv_vrt		= Cp_IO->getVariant();
	wsUintCst		N_vrt		= (wsUintCst)Xv_vrt.size();
	int*		Np_data		= (int *)Vp_data;
	wsUint		N_s			= (wsUint)Np_data[0];
	wsUint		N_e			= (wsUint)Np_data[1];

	if (!Na_data) {
		wsFloat**	Ra_dsg		= Cp_IO->getDosage();
		if (!IS_ASSIGNED(dosage) && !IS_ASSIGNED(expression))
			halt("SYSERR: Input is not dosage/expression but genotype data is missing!");

		wsVec	Ra_ma		= sseEmptyVec(N_e-N_s+1);
		wsVec	Ra_maAll	= sseEmptyVec(N_e-N_s+1);

		if (Cp_IO->isDataComplete()) {
			LOOP(i, N_sample) {
				if (Ba_founder[i]) for (wsUint j=N_s, k=0 ; j<N_e ; j++, k++) {
					wsFloat	R_geno	= Ra_dsg[i][j];

					Ra_ma[k]	+= R_geno;
					Ra_maAll[k]	+= R_geno;
				} else for (wsUint j=N_s, k=0 ; j<N_e ; j++, k++) {
					wsFloat	R_geno	= Ra_dsg[i][j];
					Ra_maAll[k]	+= R_geno;
				}
			}
			for (wsUint j=N_s, k=0 ; j<N_e ; j++, k++) {
				xMaf&		X_maf		= Xp_maf[j];
				xVariant&	X_vrt		= Xv_vrt[j];
				wsReal		R_maf		= (wsReal)Ra_ma[k] / (wsReal)(N_founder<<1);
				wsReal		R_mafAll	= (wsReal)Ra_maAll[k] / (wsReal)(N_sample<<1);

				if (R_maf > REAL_CONST(0.5)) {
					/* Flip */
					LOOP(i, N_sample) Ra_dsg[i][j] = (wsFloat)W2 - Ra_dsg[i][j];
					R_maf		= W1 - R_maf;
					R_mafAll	= W1 - R_mafAll;
					X_maf.B_flipped = true;
					if (OPT_ENABLED(indel)) {
						char *Sp_tmp = X_vrt.indel1;
						X_vrt.indel1 = X_vrt.indel2;
						X_vrt.indel2 = Sp_tmp;
					} else {
						char S_tmp = X_vrt.al1;
						X_vrt.al1 = X_vrt.al2;
						X_vrt.al2 = S_tmp;
					}
				} else
					X_maf.B_flipped = false;

				if (Xm_newMAF.size()) {
					/* --maf */
					mStrReal_it X_find = Xm_newMAF.find(X_vrt.name);
					if (X_find == Xm_newMAF.end())
						halt("Pre-defined MAF of variant [%s] is not exists!",
							X_vrt.name);

					X_maf.R_maf = X_maf.R_allMaf = X_find->second;
				} else {
					X_maf.R_maf	= R_maf;
					X_maf.R_allMaf	= R_mafAll;
				}
				X_maf.N_allele	= N_sample<<1;
				/* FIXME: Currently it is rounded */
				X_maf.N_mac		= (wsUint)round(Ra_ma[k]);
				X_maf.N_allMac	= (wsUint)round(Ra_maAll[k]);

				if (N_idx == 0 && (j%1000) == 0) {
					int S = 0;
					for (int i=0 ; i<OPT_NUMBER(thread) ; i++)
						S += Np_data[i*3+2];
					notice("%d/%d variants processed...\r", S, N_vrt);
				}
				Np_data[2]++;
			}
		} else {
			wsUint* Na_N = NULL;
			wsUint* Na_Nall = NULL;
			wsCalloc(Na_N, wsUint, N_e-N_s+1);
			wsCalloc(Na_Nall, wsUint, N_e-N_s+1);

			for (wsUint i=0 ; i<N_sample ; i++) {
				if (Ba_founder[i]) for (wsUint j=N_s, k=0 ; j<N_e ; j++, k++) {
					wsFloat R_geno = Ra_dsg[i][j];
					if (isMissing(R_geno)) continue;

					Ra_ma[k] += R_geno;
					Ra_maAll[k] += R_geno;
					Na_Nall[k]++;
					Na_N[k]++;
				} else for (wsUint j=N_s, k=0 ; j<N_e ; j++, k++) {
					wsFloat R_geno = Ra_dsg[i][j];
					if (isMissing(R_geno)) continue;

					Ra_maAll[k] += R_geno;
					Na_Nall[k]++;
				}
			}

			for (wsUint j=N_s, k=0 ; j<N_e ; j++, k++) {
				xVariant&	X_vrt	= Xv_vrt[j];
				xMaf&		X_maf	= Xp_maf[j];

				wsReal R_maf	= Na_N[k] ? (wsReal)Ra_ma[k] / (wsReal)(Na_N[k]<<1) : W0;
				wsReal R_mafAll	= Na_Nall[k] ? (wsReal)Ra_maAll[k] / (wsReal)(Na_Nall[k]<<1) : W0;
				if (R_maf > REAL_CONST(0.5)) {
					/* Flip */
					LOOP(i, N_sample) {
						wsFloat N_geno = Ra_dsg[i][j];
						if (isMissing(N_geno)) continue;
						Ra_dsg[i][j] = (wsFloat)W2 - Ra_dsg[i][j];
					}
					R_maf		= W1 - R_maf;
					R_mafAll	= W1 - R_mafAll;
					X_maf.B_flipped = true;
					if (OPT_ENABLED(indel)) {
						char *Sp_tmp = X_vrt.indel1;
						X_vrt.indel1 = X_vrt.indel2;
						X_vrt.indel2 = Sp_tmp;
					} else {
						char S_tmp = X_vrt.al1;
						X_vrt.al1 = X_vrt.al2;
						X_vrt.al2 = S_tmp;
					}
				} else
					X_maf.B_flipped = false;

				if (Xm_newMAF.size()) {
					/* --maf */
					mStrReal_it X_find = Xm_newMAF.find(X_vrt.name);
					if (X_find == Xm_newMAF.end())
						halt("Pre-defined MAF of variant [%s] is not exists!",
						X_vrt.name);

					X_maf.R_maf = X_maf.R_allMaf = X_find->second;
				} else {
					X_maf.R_maf	= R_maf;
					X_maf.R_allMaf	= R_mafAll;
				}
				X_maf.N_allele	= Na_Nall[k]<<1;
				/* FIXME: Currently it is rounded */
				X_maf.N_mac		= (wsUint)round(Ra_ma[k]);
				X_maf.N_allMac	= (wsUint)round(Ra_maAll[k]);

				if (N_idx == 0 && (j%1000) == 0) {
					int S = 0;
					for (int i=0 ; i<OPT_NUMBER(thread) ; i++)
						S += Np_data[i*3+2];
					notice("%d/%d variants processed...\r", S, N_vrt);
				}
				Np_data[2]++;
			}

			DEALLOC(Na_N);
			DEALLOC(Na_Nall);
		}

		sseFree(Ra_ma);
		sseFree(Ra_maAll);
		return 0;
	}

	wsUint* Na_ma = NULL;
	wsUint* Na_maAll = NULL;
	wsCalloc(Na_ma, wsUint, N_e-N_s+1);
	wsCalloc(Na_maAll, wsUint, N_e-N_s+1);

	if (Cp_IO->isDataComplete()) {
		LOOP (i, N_sample) {
			if (Ba_founder[i]) for (wsUint j=N_s,k=0 ; j<N_e ; j++,k++) {
				char	N_geno	= Na_data[i][j];

				Na_ma[k]	+= N_geno;
				Na_maAll[k]	+= N_geno;
			} else for (wsUint j=N_s,k=0 ; j<N_e ; j++,k++) {
				char	N_geno	= Na_data[i][j];

				Na_maAll[k]	+= N_geno;
			}
		}
		for (wsUint j=N_s,k=0 ; j<N_e ; j++,k++) {
			xVariant&	X_vrt		= Xv_vrt[j];
			xMaf&		X_maf		= Xp_maf[j];
			wsReal		R_maf		= (wsReal)Na_ma[k] / (wsReal)(N_founder<<1);
			wsReal		R_mafAll	= (wsReal)Na_maAll[k] / (wsReal)(N_sample<<1);

			if (R_maf > REAL_CONST(0.5)) {
				/* Flip */
				LOOP (i, N_sample) Na_data[i][j] = ((Na_data[i][j]-1)^0x1)&0x03;
				R_maf		= W1 - R_maf;
				R_mafAll	= W1 - R_mafAll;
				X_maf.B_flipped = true;
				if (OPT_ENABLED(indel)) {
					char *Sp_tmp = X_vrt.indel1;
					X_vrt.indel1 = X_vrt.indel2;
					X_vrt.indel2 = Sp_tmp;
				} else {
					char S_tmp = X_vrt.al1;
					X_vrt.al1 = X_vrt.al2;
					X_vrt.al2 = S_tmp;
				}
			} else
				X_maf.B_flipped = false;

			if (Xm_newMAF.size()) {
				/* --maf */
				mStrReal_it X_find = Xm_newMAF.find(X_vrt.name);
				if (X_find == Xm_newMAF.end())
					halt("Pre-defined MAF of variant [%s] is not exists!",
					X_vrt.name);

				X_maf.R_maf = X_maf.R_allMaf = X_find->second;
			} else {
				X_maf.R_maf	= R_maf;
				X_maf.R_allMaf	= R_mafAll;
			}
			X_maf.N_allele	= N_sample<<1;
			X_maf.N_mac		= Na_ma[k];
			X_maf.N_allMac	= Na_maAll[k];

			if (N_idx == 0 && (j%1000) == 0) {
				int S = 0;
				for (int i=0 ; i<OPT_NUMBER(thread) ; i++)
					S += Np_data[i*3+2];
				notice("%d/%d variants processed...\r", S, N_vrt);
			}
			Np_data[2]++;
		}
	} else {
		LOG("Sample size %d founder %d\n", N_sample, N_founder);
		wsUint* Na_N = NULL;
		wsUint* Na_Nall = NULL;
		wsCalloc(Na_N, wsUint, N_e-N_s+1);
		wsCalloc(Na_Nall, wsUint, N_e-N_s+1);

		for (wsUint i=0 ; i<N_sample ; i++) {
			if (Ba_founder[i]) for (wsUint j=N_s,k=0 ; j<N_e ; j++,k++) {
				char N_geno = Na_data[i][j];
				if (isMissing(N_geno)) continue;

				Na_ma[k] += N_geno;
				Na_maAll[k] += N_geno;
				Na_Nall[k]++;
				Na_N[k]++;
			} else for (wsUint j=N_s,k=0 ; j<N_e ; j++,k++) {
				char N_geno = Na_data[i][j];
				if (isMissing(N_geno)) continue;

				Na_maAll[k] += N_geno;
				Na_Nall[k]++;
			}
		}

		for (wsUint j=N_s,k=0 ; j<N_e ; j++,k++) {
			xVariant&	X_vrt		= Xv_vrt[j];
			xMaf&		X_maf		= Xp_maf[j];
			wsReal		R_maf		= Na_N[k] ? (wsReal)Na_ma[k] / (wsReal)(Na_N[k]<<1) : W0;
			wsReal		R_mafAll	= Na_Nall[k] ? (wsReal)Na_maAll[k] / (wsReal)(Na_Nall[k]<<1) : W0;

			if (R_maf > REAL_CONST(0.5)) {
				/* Flip */
				LOOP (i, N_sample) {
					char N_geno = Na_data[i][j];
					if (isMissing(N_geno)) continue;
					Na_data[i][j] = ((N_geno-1)^0x1)&0x03;
				}
				R_maf		= W1 - R_maf;
				R_mafAll	= W1 - R_mafAll;
				X_maf.B_flipped = true;
				if (OPT_ENABLED(indel)) {
					char *Sp_tmp = X_vrt.indel1;
					X_vrt.indel1 = X_vrt.indel2;
					X_vrt.indel2 = Sp_tmp;
				} else {
					char S_tmp = X_vrt.al1;
					X_vrt.al1 = X_vrt.al2;
					X_vrt.al2 = S_tmp;
				}
			} else
				X_maf.B_flipped = false;

			if (Xm_newMAF.size()) {
				/* --maf */
				mStrReal_it X_find = Xm_newMAF.find(X_vrt.name);
				if (X_find == Xm_newMAF.end())
					halt("Pre-defined MAF of variant [%s] is not exists!",
					X_vrt.name);

				X_maf.R_maf = X_maf.R_allMaf = X_find->second;
			} else {
				X_maf.R_maf	= R_maf;
				X_maf.R_allMaf	= R_mafAll;
			}
			X_maf.N_allele	= Na_Nall[k]<<1;
			X_maf.N_mac		= Na_ma[k];
			X_maf.N_allMac	= Na_maAll[k];

			if (N_idx == 0 && (j%1000) == 0) {
				int S = 0;
				for (int i=0 ; i<OPT_NUMBER(thread) ; i++)
					S += Np_data[i*3+2];
				notice("%d/%d variants processed...\r", S, N_vrt);
			}
			Np_data[2]++;
		}

		DEALLOC(Na_N);
		DEALLOC(Na_Nall);
	}

	DEALLOC(Na_ma);
	DEALLOC(Na_maAll);

	return 0;
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

int hweCalc(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	cIO*		Cp_IO		= (cIO *)Vp_shareData;
	wsUint		N_sample	= Cp_IO->sizeSample();
	wsUint		N_vrt		= Cp_IO->sizeVariant();
	char**		Na_data		= Cp_IO->getGenotype();
	unsigned char**
				Na_raw		= Cp_IO->getRawGeno();
	vVariant&	Xv_vrt		= Cp_IO->getVariant();
	const char*	Ba_founder	= Cp_IO->getIsFounder();
	int*		Np_data		= (int *)Vp_data;
	wsUint		N_s = (wsUint)Np_data[0];
	wsUint		N_e = (wsUint)Np_data[1];
	wsUint*		Np_res = (wsUint *)Vp_result;
	char		B_complete = Cp_IO->isDataComplete();
	Np_res[0] = Np_res[1] = Np_res[2] = 0;

	//LOG("Thread %d : Range %d~%d\n", N_idx, N_s, N_e);
	for (wsUint i = N_s; i < N_e; i++) {
		wsUint		N = 0;
		wsUint		S = 0;
		xVariant&	X_vrt = Xv_vrt[i];
		//		wsUint	Np_res[3] = { 0, 0, 0 };

		if (Na_raw) {
			// In case of there is raw genotype
			wsUint N_mis = 0;
			VARIANT_SSE_PART(N_sample, Na_raw, s)
				__m128i s2 = _mm_set_epi32(0xaaaaaaaa, 0xaaaaaaaa,
				0xaaaaaaaa, 0xaaaaaaaa);
			__m128i s1 = _mm_set_epi32(0x55555555, 0x55555555,
				0x55555555, 0x55555555);

			// (s&2) >> 1 == # of min. hom. + # of missing
			// FIXME : Linux compatibility
// 			__m128i u1 = _mm_set_epi32(0, 0, 0, 1);
// 			__m128i n2 = _mm_srl_epi32(_mm_and_si128(*s, s2), u1);
// 			Np_res[2] += _mm_popcnt_u64(n2.m128i_i64[0]) + _mm_popcnt_u64(n2.m128i_i64[1]);
// 
// 			// (s&1) == # of maj. hom. + # of missing
// 			__m128i n1 = _mm_and_si128(*s, s1);
// 			Np_res[1] += _mm_popcnt_u64(n1.m128i_i64[0]) + _mm_popcnt_u64(n1.m128i_i64[1]);
// 
// 			// (n1 & n2) == # of missing
// 			__m128i nm = _mm_and_si128(n1, n2);
// 			N_mis += _mm_popcnt_u64(n1.m128i_i64[0]) + _mm_popcnt_u64(n1.m128i_i64[1]);
			VARIANT_NRM_PART(N_sample, Na_raw[i], NN)
				unsigned char n2 = (NN & 0xaa) >> 1;
			unsigned char n1 = NN & 0x55;
			unsigned char nm = n1&n2;
			Np_res[2] += popcount_tab[n2];
			Np_res[1] += popcount_tab[n1];
			N_mis += popcount_tab[nm];
			VARIANT_END_PART(N_sample, Na_raw[i], NN)
				unsigned char n2 = (NN & 0xaa) >> 1;
			unsigned char n1 = NN & 0x55;
			unsigned char nm = n1&n2;
			Np_res[2] += popcount_tab[n2];
			Np_res[1] += popcount_tab[n1];
			N_mis += popcount_tab[nm];

			N = N_sample - N_mis;
		} else {
			/* Get MA count */
			if (!stricmp(OPT_STRING(hwe), "all")) {
				if (B_complete) {
					for (wsUint j = 0; j < N_sample; j++) {
						S += Na_data[j][i];
						Np_res[(wsUint)Na_data[j][i]]++;
					}
					N = N_sample;
				} else {
					for (wsUint j = 0; j < N_sample; j++) {
						int N_genotype = Na_data[j][i];

						if (isAvailable(N_genotype)) {
							S += N_genotype;
							Np_res[N_genotype]++;
							N++;
						}
					}
				}
			} else if (!stricmp(OPT_STRING(hwe), "founder")) {
				if (B_complete) {
					for (wsUint j = 0; j < N_sample; j++) {
						if (Ba_founder[j] == 0) continue;

						S += Na_data[j][i];
						Np_res[(wsUint)Na_data[j][i]]++;
						N++;
					}
				} else {
					for (wsUint j = 0; j < N_sample; j++) {
						if (Ba_founder[j] == 0) continue;
						int N_genotype = Na_data[j][i];

						if (isAvailable(N_genotype)) {
							S += N_genotype;
							Np_res[N_genotype]++;
							N++;
						}
					}
				}
			}
		}

		wsReal	R_maf = S / ((wsReal)N*W2);

		wsReal R_q = N*R_maf*R_maf;
		wsReal R_p = N*(W1 - R_maf)*(W1 - R_maf);
		wsReal R_2pq = N*W2*R_maf*(W1 - R_maf);

		/* Get chi^2 statistics */
		wsReal R_chistat =
			(Np_res[0] - R_p)*(Np_res[0] - R_p) / R_p +
			(Np_res[1] - R_2pq)*(Np_res[1] - R_2pq) / R_2pq +
			(Np_res[2] - R_q)*(Np_res[2] - R_q) / R_q;

#if 0
		REAL_t R_pVal = (REAL_t)chiprobP(R_chistat, W1);
#else
		wsReal R_pVal = variantHWE(Np_res[1], Np_res[0], Np_res[2]);
#endif
		Np_res[3] = N;
		((float *)Np_res)[4] = (float)R_chistat;
		X_vrt.hwe = R_pVal;

		if (N_idx == 0 && (i % 1000) == 0) {
			int S = 0;
			for (int j = 0; j < OPT_NUMBER(thread); j++)
				S += Np_data[j * 3 + 2];
			notice("%d/%d variants processed...\r", S, N_vrt);
		}
		Np_data[2]++;
	}

	return 0;
}

inline char _getIndel(xVariant &X_snv, char *Sp_gen1, char *Na_acgt)
{
	if (X_snv.indel1 == NULL) {
		X_snv.indel1 = strdup(Sp_gen1);
		/* --acgt */
		if (Na_acgt) for (char *x=X_snv.indel1 ; *x ; x++) {
			if (!Na_acgt[(wsUint)*x]) halt("Character [%c] in variant [%s] is not in --acgt",
				X_snv.name, *x);
			*x = Na_acgt[(wsUint)*x];
		}
	} else {
		/* Have indel1 */

		/* Same with indel1? */
		if (stricmp(X_snv.indel1, Sp_gen1)) {
			/* indel1 != this */

			if (X_snv.indel2 == NULL) {
				/* indel2 == NULL */
				X_snv.indel2 = strdup(Sp_gen1);

				/* --acgt */
				if (Na_acgt && X_snv.indel2) for (char *x=X_snv.indel2 ; *x ; x++) {
					if (!Na_acgt[(wsUint)*x]) halt("Character [%c] in variant [%s] is not in --acgt",
						X_snv.name, *x);
					*x = Na_acgt[(wsUint)*x];
				}
			} else if (stricmp(X_snv.indel2, Sp_gen1))
				/* Trialleleic */
				return -1;

			/* Match with second allele */
			return 1;
		}
	}

	return 0;
}

inline char _getVariant(xVariant &X_snv, char C_geno, char *Na_acgt)
{
	if (!X_snv.al1) {
		X_snv.al1 = C_geno;
		/* --acgt */
		if (Na_acgt) {
			if (!Na_acgt[(wsUint)X_snv.al1]) halt("Character [%c] in variant [%s] is not in --acgt",
				X_snv.al1, X_snv.name);
			X_snv.al1 = Na_acgt[(wsUint)X_snv.al1];
		}
	} else {
		/* Have indel1 */

		/* Same with indel1? */
		if (X_snv.al1 != C_geno) {
			/* indel1 != this */

			if (!X_snv.al2) {
				/* indel2 == NULL */
				X_snv.al2 = C_geno;

				/* --acgt */
				if (Na_acgt) {
					if (!Na_acgt[(wsUint)X_snv.al2]) halt("Character [%c] in variant [%s] is not in --acgt",
						X_snv.al2, X_snv.name);
					X_snv.al2 = Na_acgt[(wsUint)X_snv.al2];
				}
			} else if (X_snv.al2 != C_geno)
				/* Trialleleic */
				return -1;

			/* Match with second allele */
			return 1;
		}
	}

	return 0;
}

/* Checklists for PED-like file loading (i.e., having 2 datapoints for
 * each genotype is like below:
 * 
 * 1. Characters consisting one genotype should be two or three types,
 *  including the characters representing missing genotype.
 * 2. All genotypes should be coded as 0/1/2/3 or A/C/G/T, or the combinations
 *  of above characters
 */

using namespace std;

int forAllVariant(int N_mode, int N_thread, void *Vp_thrData, wsUint L_thrData,
	void *Vp_shareData, wsUint *Na_waitQueue)
{
	int		*Np_idx = (int *)Vp_thrData;
	static wsUint N_idxVrt;
	cIO *Cp_IO = (cIO *)Vp_shareData;
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		N_idxVrt = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++) {
			if ((N_idxVrt+i) >= Cp_IO->sizeVariant())
				break;
			Np_idx[i] = N_idxVrt+i;
			if (((N_idxVrt+i)%1000) == 0)
				notice("%d/%d variants processed\r", N_idxVrt+i, Cp_IO->sizeVariant());
		}
		N_idxVrt += N_thread;
		break;
	case DISTFUN_UNINIT:
		LOG("%d/%d variants processed\n", Cp_IO->sizeVariant(), Cp_IO->sizeVariant());
		break;
	case DISTFUN_AFTERLOOP:
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}
	return i;
}

int forAllSample(int N_mode, int N_thread, void *Vp_data, void *Vp_shareData)
{
	int		*Np_idx = (int *)Vp_data;
	static wsUint N_idxSample;
	cIO *Cp_IO = (cIO *)Vp_shareData;
	int i=0;

	switch (N_mode) {
	case DISTFUN_INIT:
		N_idxSample = 0;
		break;
	case DISTFUN_DISTRIBUTE:
		for ( ; i<N_thread ; i++) {
			if ((N_idxSample+i) >= Cp_IO->sizeSample())
				break;
			Np_idx[i] = N_idxSample+i;
			if (((N_idxSample+i)%10) == 0)
				notice("%d/%d samples processed\r", N_idxSample+i, Cp_IO->sizeSample());
		}
		N_idxSample += N_thread;
		break;
	case DISTFUN_UNINIT:
		LOG("%d/%d samples processed\n", Cp_IO->sizeSample(), Cp_IO->sizeSample());
		break;
	}
	return i;
}

void exportSampleMatrix(const char *S_ext, cIO *Cp_IO,
	const char *Ba_filter, char B_filtVal,
	wsUint N_mat, char **Sp_name, wsReal **Ra_mats,
	wsStrCst S_desc="")
{
	wsUint			i, j, I;
	/* Prepare format */
	wsString		S_fmt = "ss";
	S_fmt.append('r', N_mat);
	/* Prepare column names */
	vStr			S_colnames;
	S_colnames.push_back("FID");
	S_colnames.push_back("IID");
	LOOP (ii, N_mat) S_colnames.push_back(Sp_name[ii]);
	/* Init exporter */
	cTableExporter	C_mat(S_ext, S_colnames, S_fmt.get(), S_desc);

	/* Export contents */
	vSampPtr&	Xa_sample	= Cp_IO->getSample();
	wsUint		N_row		= Cp_IO->sizeSample();
	for (i=I=0 ; i<N_row ; i++) {
		if (Ba_filter && Ba_filter[i] == B_filtVal) continue;

		C_mat.put(2, Xa_sample[i]->S_FID.c_str(), Xa_sample[i]->S_IID.c_str());
		/* Export all contents in each corresponding matrix */
		for (j=0 ; j<N_mat ; j++) C_mat.put(1, Ra_mats[j][I]);
		C_mat.next();
		I++;
	}
}

void exportSampleMatrix(const char *S_ext, cIO *Cp_IO,
	const char *Ba_filter, char B_filtVal,
	wsUint N_mat, char **Sp_name, wsReal **Ra_mats)
{
	wsUint		i, j, I, N_sumCol = N_mat;
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */

	//	LOG("Export matrix [%d * %d]\n", N_row, N_sumCol);

	/* Export header */
	Cp_mat->put("FID	IID");
	for (i=0 ; i<N_sumCol ; i++)
		Cp_mat->fmt("	%s", Sp_name[i]);
	Cp_mat->put("\n");

	/* Export contents */
	vSampPtr&	Xa_sample	= Cp_IO->getSample();
	wsUint		N_row		= Cp_IO->sizeSample();
	for (i=I=0 ; i<N_row ; i++) {
		if (Ba_filter && Ba_filter[i] == B_filtVal) continue;

		Cp_mat->fmt("%s	%s", Xa_sample[i]->S_FID.c_str(),
			Xa_sample[i]->S_IID.c_str());

		/* Export all contents in each corresponding matrix */
		for (j=0 ; j<N_mat ; j++) {
			wsReal *Rp_mat = Ra_mats[j];
			Cp_mat->fmt("	%.16g", Rp_mat[I]);
		}

		Cp_mat->put("\n");
		I++;
	}

	delete Cp_mat;
}

void exportSampleMatrix(const char *S_ext, cIO *Cp_IO,
	char *B_filt, char B_filtVal,
	wsUint N_mat, ...)
{
	wsUint		i, j, I, N_sumCol = 0;
/**/cExporter*	Cp_mat = cExporter::summon(S_ext); /* CHECKED */
	wsReal		**Ra_mats	= NULL;
	char		**Sp_name	= NULL;

	va_list		H_vl;
	va_start(H_vl, N_mat);
	wsAlloc(Ra_mats, wsReal*, N_mat);
	wsAlloc(Sp_name, char*, N_mat);

	for (i=0 ; i<N_mat ; i++) {
		Sp_name[i]	= va_arg(H_vl, char*);
		Ra_mats[i]	= va_arg(H_vl, wsReal*);
		N_sumCol++;
	}

//	LOG("Export matrix [%d * %d]\n", N_row, N_sumCol);

	/* Export header */
	Cp_mat->put("FID IID");
	for (i=0 ; i<N_sumCol ; i++)
		Cp_mat->fmt(" %s", Sp_name[i]);
	Cp_mat->put("\n");

	/* Export contents */
	vSampPtr&	Xa_sample	= Cp_IO->getSample();
	wsUint		N_row		= Cp_IO->sizeSample();
	for (i=I=0 ; i<N_row ; i++) {
		if (B_filt && B_filt[i] == B_filtVal) continue;

		Cp_mat->fmt("%s %s", Xa_sample[i]->S_FID.c_str(),
			Xa_sample[i]->S_IID.c_str());

		/* Export all contents in each corresponding matrix */
		for (j=0 ; j<N_mat ; j++) {
			wsReal *Rp_mat = Ra_mats[j];
			Cp_mat->fmt(" %g", Rp_mat[I]);
		}

		Cp_mat->put("\n");
		I++;
	}

	delete Cp_mat;
	/* Terminate */
	va_end(H_vl);
}

cIO::cIO()
{
	N_missGeno			= 0;
	N_sampOrig			= 0;
	N_sample			= 0;
	N_variant			= 0;
	N_founder			= 0;
	N_pheno				= 0;
	N_cov				= 0;
	N_twin				= 0;
	Ba_isFounder		= NULL;
	Ra_pheno			= NULL;
	Ra_data				= NULL;
	Na_gProb			= NULL;
	Na_geno				= NULL;
	Ra_dosage			= NULL;
	Ra_cov				= NULL;
	Ra_covGxE			= NULL;
	Na_charData			= NULL;
	B_isDataComplete	= 0;
	N_filterSample		= FILT_NOTHING;
	N_filterVrt			= FILT_NOTHING;
	Ba_typePheno		= NULL;
	B_haveProband		= 0;
	Ra_customWeight		= NULL;
	Cp_anaPPP			= NULL;
	B_altVariantDef		= 0;
	Xv_altVariantDef	= NULL;
	Cp_impute			= NULL;
	B_dry				= 0;
	N_altSample			= 0;
	Xp_maf				= NULL;
	B_mafComputed		= 0;
	B_blupApplied		= 0;
	N_altSample			= 0;
	N_szFilteredSamp	= 0;
	N_szFilteredVrt		= 0;
	memset(Ba_chrExists, 0x00, sizeof(char)*MAX_NCHR);

	/* Init filterings */
	_initFilterings();
}

cIO::~cIO()
{
	if (N_cov) {
		LOOP (i, N_cov) DEALLOC(Xa_covInfo[i].Sp_varName);
		Xa_covInfo.clear();
	}
	FOREACH(vVariant_it, Xv_variant, i)
		DEALLOC(i->name);
	DEALLOC(Xp_maf);
	Xm_iid2data.clear();
	DEALLOC(Ra_covGxE);
	sseFree(Ra_customWeight);
	clear();
}

wsUint* cIO::_setSampOrder()
{
	wsUint* Na_idxSamp = NULL;

	if (IS_ASSIGNED(sampleorder)) {
		LOG("Samples are ordered by user-defined sequence\n");
		cStrFile C_mo(OPT_STRING(sampleorder), "Sample order file");

		/* Read markers */
		wsAlloc(Na_idxSamp, wsUint, N_sample);
		/* Init with 0xffffffff to prevent multiple */
		for (wsUint i=0 ; i<N_sample ; i++)
			Na_idxSamp[i] = 0xffffffff;

		char *Sp_buf = NULL, *a = NULL, *b = NULL;
		wsAlloc(Sp_buf, char, 2048);
		for (wsUint L=1 ; C_mo.gets(Sp_buf, 2048) ; L++) {
			int		N_insIdx = -1;
			/* FID & IID */
			getString(&Sp_buf, &a);
			getString(&a, &b);
			string IID = a;
			/* find */
			mSamp_it X_find = Xm_iid2data.find(IID);
			/* If find, record index */
			if (X_find != Xm_iid2data.end() && X_find->second.N_idx != SAMP_NODATA) {
				N_insIdx = X_find->second.N_idx;
				break;
			}
			wsStrCst Sp_IID = X_find->second.S_IID.c_str();
			/* Ignore this line if no match found */
			if (N_insIdx == -1 || a == NULL)
				continue;
			/* Check */
			if (Na_idxSamp[L-1] != 0xffffffff)
				halt("Sample [%s] in line [%d] allocated twice!", Sp_IID, L);
			/* Check FID */
			if (Na_idxSamp[L-1] != 0xffffffff)
				halt("Sample [%s] is in family [%s] but requested to family "
					"[%s] in sample order file [%s]", Sp_IID,
					X_find->second.S_FID.c_str(), Sp_buf, OPT_STRING(sampleorder));
			/* Insert */
			Na_idxSamp[L-1] = N_insIdx;
		}

		/* Check lack */
		for (wsUint i=0 ; i<N_sample ; i++)
			if (Na_idxSamp[i] == 0xffffffff)
				halt("Sample [%s] not in variant order file [%s]",
					Xa_sampleV2[i]->S_IID.c_str(), OPT_STRING(sampleorder));
		DEALLOC(Sp_buf);
		C_mo.close();
	} else if (IS_ASSIGNED(sortsample)) {
		LOG("Samples are sorted by its family ID and individual ID\n");
		wsAlloc(Na_idxSamp, wsUint, N_sample);
		wsUint		i = 0;

		/* Build index */
		xStrSort	*Xa_sampids = NULL;
		wsAlloc(Xa_sampids, xStrSort, Xa_sampleV2.size());
		FOREACHDO (vSampPtr_it, Xa_sampleV2, it, i++) {
			Xa_sampids[i].i = i;
			Xa_sampids[i].V = NULL;
			wsAlloc(Xa_sampids[i].V, char, 1024);
			sprintf(Xa_sampids[i].V, "%s~%s", (*it)->S_FID.c_str(),
				(*it)->S_IID.c_str());
		}
		if (OPT_ENABLED(natural))
			qsort(Xa_sampids, Xa_sampleV2.size(), sizeof(xStrSort),
				stricmp(OPT_STRING(sortsample), "asc")?sort_natstr_desc:sort_natstr);
		else
			qsort(Xa_sampids, Xa_sampleV2.size(), sizeof(xStrSort),
				stricmp(OPT_STRING(sortsample), "asc")?sort_str_desc:sort_str);

		/* Set index */
		for (i=0 ; i<N_sample ; i++)
			Na_idxSamp[i] = Xa_sampids[i].i;
	} else if (IS_ASSIGNED(sortiid)) {
		LOG("Samples are sorted by its individual ID\n");
		wsAlloc(Na_idxSamp, wsUint, N_sample);
		wsUint		i = 0;

		/* Build index */
		xStrSort	*Xa_iids = NULL;
		wsAlloc(Xa_iids, xStrSort, Xa_sampleV2.size());
		FOREACHDO (vSampPtr_it, Xa_sampleV2, it, i++) {
			Xa_iids[i].i = i;
			Xa_iids[i].V = strdup((*it)->S_IID.c_str());
		}
		if (OPT_ENABLED(natural))
			qsort(Xa_iids, Xa_sampleV2.size(), sizeof(xStrSort),
				stricmp(OPT_STRING(sortiid), "asc")?sort_natstr_desc:sort_natstr);
		else
			qsort(Xa_iids, Xa_sampleV2.size(), sizeof(xStrSort),
				stricmp(OPT_STRING(sortiid), "asc")?sort_str_desc:sort_str);

		/* Set index */
		for (i=0 ; i<N_sample ; i++)
			Na_idxSamp[i] = Xa_iids[i].i;
	}
	
	return Na_idxSamp;
}

wsUint* cIO::_setVrtOrder()
{
	wsUint* Na_idxMarker = NULL;

	if (IS_ASSIGNED(variantorder)) {
		LOG("Variants are ordered with user-defined sequence\n");

		cStrFile C_mo(OPT_STRING(variantorder), "Variants order file");
		/* Read variants */
		wsAlloc(Na_idxMarker, wsUint, N_variant);
		/* Init with 0xffffffff to prevent multiple */
		for (wsUint i=0 ; i<N_variant ; i++)
			Na_idxMarker[i] = 0xffffffff;

		/* Build marker map */
		mDataIdx	Xm_mmap[MAX_NCHR+1];
		wsUint		j = 0;
		FOREACHDO (vVariant_it, Xv_variant, i, j++) {
			if (i->chr > 0 && (wsUint)i->chr <= NCHR_SPECIES)
				Xm_mmap[i->chr].insert(make_pair((string)i->name, j));
			else
				Xm_mmap[0].insert(make_pair((string)i->name, j));
		}

		char *Sp_buf = NULL, *a = NULL;
		wsAlloc(Sp_buf, char, 2048);
		for (wsUint L=1 ; C_mo.gets(Sp_buf, 2048) ; L++) {
			int		N_insIdx = -1;
			/* name only */
			getString(&Sp_buf, &a);
			/* find */
			for (wsUint N_chr=0 ; N_chr<=NCHR_SPECIES ; N_chr++) {
				mDataIdx_it X_find = Xm_mmap[N_chr].find(Sp_buf);
				/* If find, record index */
				if (X_find != Xm_mmap[N_chr].end()) {
					N_insIdx = X_find->second;
					break;
				}
			}
			/* Ignore this line if no match found */
			if (N_insIdx == -1 || a == NULL)
				continue;
			/* Check */
			if (Na_idxMarker[L-1] != 0xffffffff)
				halt("Variant [%s] in line [%d] allocated twice!",
					Xv_variant[N_insIdx].name, L);
			/* Insert */
			Na_idxMarker[L-1] = N_insIdx;
		}

		/* Check lack */
		for (wsUint i=0 ; i<N_variant ; i++)
			if (Na_idxMarker[i] == 0xffffffff)
				halt("Variant [%s] not in variants order file [%s]",
					Xv_variant[i].name, OPT_STRING(variantorder));
		DEALLOC(Sp_buf);
		C_mo.close();
	} else if (IS_ASSIGNED(sortvariant)) {
		LOG("Variants are sorted by its identifier\n");
		wsAlloc(Na_idxMarker, wsUint, N_variant);
		wsUint i = 0;

		/* Build index */
		xStrSort	*Xa_snvs = NULL;
		wsAlloc(Xa_snvs, xStrSort, N_variant);
		FOREACHDO (vVariant_it, Xv_variant, it, i++) {
			Xa_snvs[i].i = i;
			Xa_snvs[i].V = it->name;
		}
		if (OPT_ENABLED(natural))
			qsort(Xa_snvs, N_variant, sizeof(xStrSort),
				stricmp(OPT_STRING(sortvariant), "asc")?sort_natstr_desc:sort_natstr);
		else
			qsort(Xa_snvs, N_variant, sizeof(xStrSort),
				stricmp(OPT_STRING(sortvariant), "asc")?sort_str_desc:sort_str);

		/* Set index */
		for (i=0 ; i<N_variant ; i++)
			Na_idxMarker[i] = Xa_snvs[i].i;
	} else if (IS_ASSIGNED(sortpos)) {
		LOG("Variants are sorted by its chromosome and position\n");
		wsAlloc(Na_idxMarker, wsUint, N_variant);
		wsUint i = 0;

		/* Split by chromosomes */
		xUintSort**	Xa_vrts	= NULL;
		wsUint*		Na_mkr	= NULL;
		wsCalloc(Na_mkr, wsUint, MAX_NCHR+1);
		wsAlloc(Xa_vrts, xUintSort*, MAX_NCHR+1);
		/* Count each chromosome */
		FOREACHDO (vVariant_it, Xv_variant, it, i++)
			Na_mkr[it->chr > 0 && (wsUint)it->chr <= NCHR_SPECIES ? it->chr : 0]++;
		/* Allocate memory */
		for (wsUint ii=0 ; ii<=NCHR_SPECIES ; ii++)
			wsAlloc(Xa_vrts[ii], xUintSort, Na_mkr[ii]);
		/* Build index... now */
		FOREACHDO (vVariant_it, Xv_variant, it, i++) {
			xUintSort *Xp_vrt = Xa_vrts[it->chr > 0 && (wsUint)it->chr <= NCHR_SPECIES ? it->chr : 0];
			Xp_vrt[i].i = i;
			Xp_vrt[i].V = it->pos;
		}

		/* For each chromosome, do sort */
		for (wsUint i=0 ; i<=NCHR_SPECIES ; i++)
			qsort(Xa_vrts[i], Na_mkr[i], sizeof(xStrSort),
				stricmp(OPT_STRING(sortpos), "asc")?sort_uint_desc:sort_uint);

		/* Set index from mapped chrs */
		wsUint k = 0;
		for (wsUint i=1 ; i<=NCHR_SPECIES ; i++) {
			xUintSort *Xp_vrt = Xa_vrts[i];
			for (wsUint j=0 ; j<Na_mkr[i] ; j++,k++)
				Na_idxMarker[k] = Xp_vrt[j].i;
		}
		/* Then do for Un */ {
			xUintSort *Xp_vrt = Xa_vrts[0];
			for (wsUint j=0 ; j<Na_mkr[i] ; j++,k++)
				Na_idxMarker[k] = Xp_vrt[j].i;
		}
	}

	return Na_idxMarker;
}

#define VARRNG_SZWINDOW	1000
void cIO::_resizeData()
{
	cTimer		t;
	wsUint		N_afterVrt		= N_variant - N_szFilteredVrt;
	wsUint		N_afterSample	= N_sample - N_szFilteredSamp;
	vVariant&	Xv_vrt			= getVariant();

	/* 
	 * Option validity check
	 */

	/* --sortvariant, --sortsample, --sortiid */
	/* --sortiid and --sortsample are mutually exclusive */
	if (IS_ASSIGNED(sortiid) && IS_ASSIGNED(sortsample))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--sortiid", "--sortsample");
	if (IS_ASSIGNED(sortpos) && IS_ASSIGNED(sortvariant))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--sortpos", "--sortvariant");
	/* --sortiid & --sampleorder, --sortsample & --sampleorder are mut. exc. */
	if (IS_ASSIGNED(sortiid) && IS_ASSIGNED(sampleorder))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--sortiid", "--sampleorder");
	if (IS_ASSIGNED(sortsample) && IS_ASSIGNED(sampleorder))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--sortsample", "--sampleorder");
	/* --sortpos or --sortvariant & --markerorder are mut. exc */
	if (IS_ASSIGNED(sortpos) && IS_ASSIGNED(variantorder))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--sortpos", "--variantorder");
	if (IS_ASSIGNED(sortvariant) && IS_ASSIGNED(variantorder))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--sortvariant", "--variantorder");

	/* if --window, set range */
	if (IS_ASSIGNED(window)) {
		wsUint			N_wndFilt	= 0;
		xOptRange&		X_rng		= OPT_RANGE(window);
		vector<vInt>*	Na_filtIdx	= NULL;
		Na_filtIdx = new vector<vInt>[NCHR_SPECIES];

		/* Set range */
		t.start();
		wsUint I = 0;
		FOREACHDO (vBool_it, Bv_filtVrt, i, I++) {
			int N_chr = Xv_vrt[I].chr;
			/* Do not treat when chr == Un or not filtered*/
			if (N_chr == 0 || *i == 0) continue;

			/* Get index and check size */
			wsUint N_1 = (wsUint)(Xv_vrt[I].pos/VARRNG_SZWINDOW);
			wsUint N_2 = (wsUint)((Xv_vrt[I].pos+X_rng.R_s)/VARRNG_SZWINDOW);
			wsUint N_3 = (wsUint)((Xv_vrt[I].pos+X_rng.R_e)/VARRNG_SZWINDOW);
			if (Na_filtIdx[N_chr].size() >= N_3)
				Na_filtIdx[N_chr].resize(N_3+1);

			/* Mark as exist */
			Na_filtIdx[N_chr][N_1].push_back(I);
			if (N_1 != N_2) Na_filtIdx[N_chr][N_2].push_back(I);
			if (N_3 != N_1) Na_filtIdx[N_chr][N_3].push_back(I);
		}
		pverbose("Indexing took [%s]\n", t.getReadable());

		/* For each non-filtered, evaluate */
		FOREACHDO (vBool_it, Bv_filtVrt, i, I++) {
			int N_chr = Xv_vrt[I].chr;
			/* Do not treat when chr == Un or filtered*/
			if (N_chr == 0 || *i != 0) continue;

			/* Check whether this is ranged */
			wsUint N_1 = (wsUint)(Xv_vrt[I].pos/VARRNG_SZWINDOW);
			vInt& Xv_idx = Na_filtIdx[N_chr][N_1];
			if (Xv_idx.size() == 0) continue;

			/* Search the range and do filtering */
			FOREACH (vInt_it, Xv_idx, j) if (isInRange(X_rng, Xv_vrt[I].pos - Xv_vrt[*j].pos)) {
				Bv_filtVrt[I] = MFR_WINDOW;
				N_filterVrt++;
				N_wndFilt++;
				break;
			}
		}

		/* Memory free */
		delete [] Na_filtIdx;

		/* Print log */
		if (N_wndFilt)
			LOGnote("[%d] variants were filtered by --window\n", N_wndFilt);
	}

	/* if --filtreport, set exporter */
	if (OPT_ENABLED(filtreport)) {
		cExporter* Cp_fr = cExporter::summon("filt.variant.lst");
		headerVariant(Cp_fr);
		Cp_fr->put("	REASON\n");

		for (wsUint i=0 ; i<N_variant ; i++)
			if (Bv_filtVrt[i]) {
				entryVariant(Cp_fr, Xv_variant[i]);
				switch (Bv_filtVrt[i]) {
				case MFR_GTRATE:	Cp_fr->put("	GENOTYPE_RATE\n"); break;
				case MFR_MENDELIAN:	Cp_fr->put("	MENDELIAN_ERROR\n"); break;
				case MFR_MISCACT:	Cp_fr->put("	MISSING_TEST\n"); break;
				case MFR_RGNGENE:	Cp_fr->put("	GENE_REGION\n"); break;
				case MFR_HWE:		Cp_fr->put("	HWE_TEST\n"); break;
				case MFR_MAC:		Cp_fr->put("	MINORALLELE_CNT\n"); break;
				case MFR_MAF:		Cp_fr->put("	MINORALLELE_FREQ\n"); break;
				case MFR_WINDOW:	Cp_fr->put("	FILTER_WINDOW\n"); break;
				default:			Cp_fr->put("	UNKNOWN\n"); break;
				}
			}
		delete Cp_fr;
	}
	

	/* If --sampresize, set N_afterSample to be specific value */
	if (IS_ASSIGNED(sampresize)) {
		wsReal R_valSS		= OPT_REAL(sampresize);
		wsUint N_finSamp	= R_valSS >= W1 ?
			(wsUint)R_valSS : (wsUint)round(R_valSS*N_sample);

		/* Set N_afterSample >= OPT_NUMBER(sampresize) */
		if (N_afterSample <= N_finSamp) {
			LOGwarn("Desired reduced sample size[%d] is larger than "
				"the number of remained samples[%d]\n", N_finSamp,
				N_afterSample);
		} else {
			wsUint N_rdcSamp = N_afterSample - N_finSamp;
			while (N_rdcSamp) {
				/* Randomly select one individual */
				wsUint N_idx = wsRand() % N_sample;
				if (Bv_filtSamp[N_idx]) continue;

				/* Mark this sample */
				Bv_filtSamp[N_idx] = 1;
				N_rdcSamp--;
				N_szFilteredSamp++;
			}
			LOG("Samples reduced from [%d] samples to [%d] samples\n", N_afterSample,
				N_finSamp);
			N_afterSample = N_sample - N_szFilteredSamp;
		}
	}

	LOG("Post-filtering step, %d variants and %d samples should be removed from"
		" analysis\n", N_szFilteredVrt, N_szFilteredSamp);
	t.start();

	/* Unset mafComputed */
	B_mafComputed = 0;

	/* Allocate */
	char	**Na_newData = NULL;
	wsAlloc(Na_newData, char*, N_afterSample);

	/* Set ordering-related indices */
	wsUint*	Na_idxSamp		= _setSampOrder();
	wsUint*	Na_idxMarker	= _setVrtOrder();

	/* Remap data
	 * i  = idx of original data
	 * _i = idx of after data
	 */
	vSampPtr Xa_cpSamp;
	Xa_cpSamp.resize(Xa_sampleV2.size());
	copy(Xa_sampleV2.begin(), Xa_sampleV2.end(), Xa_cpSamp.begin());

	/* --randnasamp */
	if (IS_ASSIGNED(randnasamp)) {
		if (IS_ASSIGNED(nasamp))
			halt_fmt(WISARD_CANT_EXCL_OPT, "--nasamp", "--randnasamp");

		wsReal R_naSamp = OPT_REAL(randnasamp);
		wsUint N_naSamp = 0;
		
		if (R_naSamp < W1)	///< proportion
			N_naSamp = (wsUint)((wsReal)N_sample * R_naSamp + REAL_CONST(0.5));
		/* In case of exact number (floored) */
		else				///< exact number (floored)
			N_naSamp = (wsUint)R_naSamp;
		LOG("[%d] samples randomly nullified\n", N_naSamp);
		for (wsUint i=0 ; i<N_naSamp ; i++) {
			xSample* Xp_s = NULL;
			string ID;
			char B_useFull = mapIsInexist(Xe_listSample, "__LIST__");
			do {
				Xp_s = Xa_sampleV2[wsRand()%N_sample];
				ID = B_useFull ? Xp_s->S_FID+","+Xp_s->S_IID : Xp_s->S_IID;
			} while (N_filterSample == FILT_REMOVE && mapIsExist(Xe_listSample, ID));
			Xe_listNAsamp.insert(Xp_s->S_FID+","+Xp_s->S_IID);
		}
	}

	char	B_NAlist = mapIsExist(Xe_listNAsamp, "__LIST__");
	wsUint _i=0;
	N_missGeno = 0;
	if (Na_idxSamp == NULL) {
		vSampPtr_it it=Xa_sampleV2.begin();
		for (wsUint i=0 ; i<N_sample ; i++) {
			wsUint _j=0;

			if (Bv_filtSamp[i]) {
				/* Mark as missing, and delete from data-only sample array */
				(*it)->N_idx = SAMP_NODATA;
				(*it)->B_isMissing = 1;
				it = Xa_sampleV2.erase(it);
				sseFree(Na_geno[i]);
				continue;
			}
			sseMalloc(Na_newData[_i], char, N_afterVrt);

			/* Adjust data index */
			(*it)->N_idx -= i-_i;

			/* if NAsamp, nullifying and continue */
			if ((B_NAlist && mapIsExist(Xe_listNAsamp, (*it)->S_IID)) ||
				mapIsExist(Xe_listNAsamp, (*it)->S_FID+","+(*it)->S_IID)) {
				lverbose("Sample [%s::%s] nullified\n", (*it)->S_FID.c_str(),
					(*it)->S_IID.c_str());
				memset(Na_newData[_i], (char)WISARD_NA, sizeof(char)*N_afterVrt);

				_j = N_afterVrt;
				N_missGeno += N_afterVrt;
			} else if (Na_idxMarker) for (wsUint J=0 ; J<N_variant ; J++) {
				wsUint j = Na_idxMarker[J];
				if (Bv_filtVrt[j])
					continue;
				char N_geno = Na_geno[i][j];

				if (isMissing(N_geno)) N_missGeno++;
				Na_newData[_i][_j] = N_geno;
				_j++;
			} else for (wsUint j=0 ; j<N_variant ; j++) {
				if (Bv_filtVrt[j])
					continue;
				char N_geno = Na_geno[i][j];

				if (isMissing(N_geno)) N_missGeno++;
				Na_newData[_i][_j] = N_geno;
				_j++;
			}
			if (_i == 0 && _j != N_afterVrt)
				halt("SYSERR: Na_data reconstruction have bug, # of variants unmatch ([%d] gathered, [%d] expected)",
					_j, N_afterVrt);
			_i++;
			it++;
			sseFree(Na_geno[i]);
		}
	} else {
		/* Clear all */
		Xa_sampleV2.clear();

		/* For all sample index */
		for (wsUint i=0 ; i<N_sample ; i++) {
			wsUint N_idxCurSamp = Na_idxSamp[i];
			wsUint _j=0;
			xSample *Xp_s = Xa_cpSamp[N_idxCurSamp];

			if (Bv_filtSamp[N_idxCurSamp]) {
				/* Mark as missing, and delete from data-only sample array */
				Xp_s->N_idx = SAMP_NODATA;
				Xp_s->B_isMissing = 1;
				sseFree(Na_geno[N_idxCurSamp]);
				continue;
			}

			sseMalloc(Na_newData[_i], char, N_afterVrt);

			/* Adjust data index */
			Xp_s->N_idx = _i;

			/* if NAsamp, nullifying and continue */
			if ((B_NAlist && mapIsExist(Xe_listNAsamp, Xp_s->S_IID)) ||
				mapIsExist(Xe_listNAsamp, Xp_s->S_FID+","+Xp_s->S_IID)) {
				lverbose("Sample [%s::%s] nullified\n", Xp_s->S_FID.c_str(),
					Xp_s->S_IID.c_str());
				memset(Na_newData[_i], (char)WISARD_NA, sizeof(char)*N_afterVrt);

				_j = N_afterVrt;
				N_missGeno += N_afterVrt;
			} else if (Na_idxMarker) for (wsUint J=0 ; J<N_variant ; J++) {
				/* Realloc */
				wsUint j = Na_idxMarker[J];
				if (Bv_filtVrt[j])
					continue;
				char N_geno = Na_geno[N_idxCurSamp][j];

				if (isMissing(N_geno)) N_missGeno++;
				Na_newData[_i][_j] = N_geno;
				_j++;
			} else for (wsUint j=0 ; j<N_variant ; j++) {
				if (Bv_filtVrt[j])
					continue;
				char N_geno = Na_geno[N_idxCurSamp][j];

				if (isMissing(N_geno)) N_missGeno++;
				Na_newData[_i][_j] = N_geno;
				_j++;
			}
			if (_i == 0 && _j != N_afterVrt)
				halt("SYSERR: Na_data reconstruction have bug, # of variants unmatch");
			_i++;
			sseFree(Na_geno[N_idxCurSamp]);
			Xa_sampleV2.push_back(Xa_cpSamp[N_idxCurSamp]);
		}
	}
	/* CHECKME : Check integrity */
	if (_i != N_afterSample) {
		halt("SYSERR: Na_data reconstruction have bug [Expected %d, retrieved %d]", N_afterSample, _i);
/**/	cExporter* Cp_err = cExporter::summon("filter.sample.lst");
		for (wsUint i=0 ; i<N_sample ; i++) {
			Cp_err->fmt("%s	%s	%s\n", Xa_cpSamp[i]->S_FID.c_str(),
				Xa_cpSamp[i]->S_IID.c_str(), Bv_filtSamp[i]?"FILTERED":"REMAINED");
		}
		delete Cp_err;
	}

	DEALLOC(Na_geno);
	Na_geno = Na_newData;
	LOG("Data reconstructed, %s taken\n", t.getReadable());

	/* E????A???? ??o??o????, validate */ {
		vSampPtr_it it=Xa_sampleV2.begin();
		for (_i=0,it=Xa_sampleV2.begin() ; it!=Xa_sampleV2.end() ; it++,_i++) {
			/* ?IA??C?????? CO */
			if ((*it)->N_idx != (int)_i)
				halt("SYSERR: Sample [%s:%s] have invalid data index %d (Expected %d)",
					(*it)->S_FID.c_str(), (*it)->S_IID.c_str(), (*it)->N_idx, _i);
		}
	}

	/* Allocate */
	t.start();
	wsReal **Ra_newPheno = sseMatrix(N_pheno, N_afterSample);

	/* Remap data
	* i = idx of original data
	* _i = idx of after data
	*/
	/* Determine the type of phenotype */
	DEALLOC(Ba_typePheno);
	wsCalloc(Ba_typePheno, char, N_pheno);
	_i = 0;
	for (wsUint i=0 ; i<N_pheno ; i++) {
		_i = 0;
		if (Na_idxSamp == NULL) for (wsUint j=0 ; j<N_sample ; j++) {
			if (Bv_filtSamp[j])
				continue;
			wsReal R_pheno = Ra_pheno[i][j];

			if (R_pheno != WISARD_NA_REAL && R_pheno != OPT_NUMBER(phenoCase)
				&& R_pheno != OPT_NUMBER(phenoCtrl))
				Ba_typePheno[i] = 1; /* Set to continuous */
			Ra_newPheno[i][_i++] = Ra_pheno[i][j];
		} else for (wsUint _j=0 ; _j<N_sample ; _j++) {
			wsUint j = Na_idxSamp[_j];
			if (Bv_filtSamp[j])
				continue;
			wsReal R_pheno = Ra_pheno[i][j];

			if (R_pheno != WISARD_NA_REAL && R_pheno != OPT_NUMBER(phenoCase)
				&& R_pheno != OPT_NUMBER(phenoCtrl))
				Ba_typePheno[i] = 1; /* Set to continuous */
			Ra_newPheno[i][_i++] = Ra_pheno[i][j];
		}
		if (_i != N_afterSample)
			halt("SYSERR: Ra_pheno reconstruction have bug");
	}
	sseUnmat(Ra_pheno, N_pheno);
	Ra_pheno = Ra_newPheno;
	LOG("Phenotype reconstructed, %s taken\n", t.getReadable());

	/* Reconstruct isFounder information */
	t.start();
	char	*Ba_newFounder = NULL;
	wsAlloc(Ba_newFounder, char, N_afterSample);
	_i = 0;
	N_founder = 0;
	if (Na_idxSamp == NULL) for (wsUint j=0 ; j<N_sample ; j++) {
		if (Bv_filtSamp[j])
			continue;
		if (Ba_isFounder[j]) N_founder++;
		Ba_newFounder[_i++] = Ba_isFounder[j];
	} else for (wsUint _j=0 ; _j<N_sample ; _j++) {
		wsUint j = Na_idxSamp[_j];
		if (Bv_filtSamp[j])
			continue;
		if (Ba_isFounder[j]) N_founder++;
		Ba_newFounder[_i++] = Ba_isFounder[j];
	}
	DEALLOC(Ba_isFounder);
	Ba_isFounder = Ba_newFounder;
	LOG("Founder information reconstructed, %s taken\n", t.getReadable());
	LOG("Final number of founders is %d\n", N_founder);

	if (N_cov) {
		t.start();
		wsReal **Ra_newCov = sseMatrix(N_cov, N_afterSample);
		/* Remap data
		* i = idx of original data
		* _i = idx of after data
		*/
		_i = 0;
		for (wsUint i=0 ; i<N_cov ; i++) {
			_i = 0;
			if (Na_idxSamp == NULL) for (wsUint j=0 ; j<N_sample ; j++) {
				if (Bv_filtSamp[j])
					continue;
				Ra_newCov[i][_i++] = Ra_cov[i][j];
			} else for (wsUint _j=0 ; _j<N_sample ; _j++) {
				wsUint j = Na_idxSamp[_j];
				if (Bv_filtSamp[j])
					continue;
				Ra_newCov[i][_i++] = Ra_cov[i][j];
			}
			if (_i != N_afterSample)
				halt("SYSERR: Ra_cov reconstruction have bug");
		}
		sseUnmat(Ra_cov, N_cov);
		Ra_cov = Ra_newCov;
		LOG("Covariate reconstructed, %s taken\n", t.getReadable());
	}

	/* Xa_vrt */
	t.start();
	vVariant Xa_newVrt;
	if (Na_idxMarker == NULL) for (wsUint i=0 ; i<N_variant ; i++) {
		if (Bv_filtVrt[i] == 0)
			Xa_newVrt.push_back(Xv_variant[i]);
		else
			DEALLOC(Xv_variant[i].name);
	} else for (wsUint _i=0 ; _i<N_variant ; _i++) {
		wsUint i = Na_idxMarker[_i];
		if (Bv_filtVrt[i] == 0)
			Xa_newVrt.push_back(Xv_variant[i]);
		else
			DEALLOC(Xv_variant[i].name);
	}
	DEALLOC(Na_idxMarker);
	DEALLOC(Na_idxSamp);

	if (Xa_newVrt.size() != N_afterVrt)
		halt("SYSERR: Xa_vrt reconstruction have bug");
	Xv_variant.clear();
	Xv_variant = Xa_newVrt;

	/* Finalize */
	N_sample	= N_afterSample;
	N_variant		= N_afterVrt;
	if (N_sample != Xa_sampleV2.size())
		halt("SYSERR: Xa_sample reconstruction have bug");
	LOG("Post-filtering step complete, data have %d variants and %d samples\n", N_variant, N_sample);
}

xModel cIO::_getModel()
{
	if (stricmp(OPT_STRING(model), "additive"))
		return MD_ADDITIVE;
	else if (stricmp(OPT_STRING(model), "recessive"))
		return MD_RECESSIVE;
	else if (stricmp(OPT_STRING(model), "dominant"))
		return MD_DOMINANT;
	else if (stricmp(OPT_STRING(model), "multiplicative"))
		return MD_MULTIPLICATIVE;
	halt("Genetic model [%s] is incorrect!");
	return MD_ERROR;
}

int cIO::isSampleFiltered(string& S_FID, string &S_IID, char B_silent/*=0*/)
{
	char S_buf[256];
	if (N_filterSample == FILT_NOTHING)
		return 0;

	strcpy(S_buf, N_filterSample==FILT_REMOVE ? "filtered" : "selected");

	// ?I ??I??????????????A ?I?U??I ??N??i??A X_sample?C fid & iid??|
	// cPlinkIO::init  ??I???????????? ??O?????O???? fid & iid ???????????I ??AA??CI????
	// ??AA???�T??e 1, ??????I??e 0??? ????AIC?????? CO??I??U
	if (OPT_ENABLED(regex)) {
		FOREACH (vPattern_it, Xv_listSample, i) {
			bool B_res = trex_match(*i, S_IID.c_str());
			if (B_res && N_filterSample == FILT_REMOVE) return 1;
			else if (B_res && N_filterSample == FILT_SELECT) return 0;
		}
		FOREACH (vPattern_it, Xv_listFam, i) {
			bool B_res = trex_match(*i, S_FID.c_str());
			if (B_res && N_filterSample == FILT_REMOVE) return 1;
			else if (B_res && N_filterSample == FILT_SELECT) return 0;
		}
	} else {
		eStr_it	X_resFind	= Xe_listSample.find(S_FID+","+S_IID);
		eStr_it	X_resFind2	= Xe_listSample.find(S_IID);
		eStr_it	X_resFind3	= Xe_listFam.find(S_FID);
		char		B_filtByIID = mapIsExist(Xe_listSample, "__LIST__");
		if (X_resFind != Xe_listSample.end() ||
			(B_filtByIID && X_resFind2 != Xe_listSample.end()) ||
			X_resFind3 != Xe_listFam.end()) {

			if (X_resFind3 != Xe_listFam.end()) {
				if (!B_silent) LOG("Family [%s] member [%s] %s\n", S_FID.c_str(),
					S_IID.c_str(), S_buf);
	//			X_resFind2->second = 2;
			} else if (B_filtByIID) {
				if (!B_silent) LOG("Sample [%s] %s\n", S_IID.c_str(), S_buf);
				// Remove it from the list
				Xe_listSample.erase(X_resFind2);
			} else {
				if (!B_silent) LOG("Sample [%s] from family [%s] %s\n",
					S_IID.c_str(), S_FID.c_str(), S_buf);
				// Remove it from the list
				Xe_listSample.erase(X_resFind);
			}
			if (N_filterSample == FILT_REMOVE)
				return 1;
			else
				return 0;
		}
	}

	if (N_filterSample == FILT_REMOVE)
		return 0;
	else
		return 1;
}

void cIO::exportLGEN()
{
	/* --outmisgeno : LGEN '0' */
	const char*	S_misgeno = IS_ASSIGNED(outmisgeno) ?
		OPT_STRING(outmisgeno) : "0";

	/* Cannot apply --outphenoonly */
	if (OPT_ENABLED(outphenoonly))
		LOGwarn("BED file does not support --outphenoonly, will be ignored\n");

	/* Allocate EXPORTERS */
	cExporter**	Cp_lgen	= NULL;
	cExporter**	Cp_map	= NULL;
	cExporter**	Cp_fam	= NULL;
	wsCalloc(Cp_lgen, cExporter*, MAX_NCHR+1);
	wsCalloc(Cp_map, cExporter*, MAX_NCHR+1);
	wsCalloc(Cp_fam, cExporter*, MAX_NCHR+1);

	if (OPT_ENABLED(split)) {
		/* Make instance of EXPOTER ONLY IF exists */
		char S_ext[64];
		for (wsUint i=0 ; i<NCHR_SPECIES ; i++) if (Ba_chrExists[i]) {
			wsStrCst S_chr = getChrName2(i + 1);
			sprintf(S_ext, "chr%s.lgen", S_chr);
			Cp_lgen[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "chr%s.map", S_chr);
			Cp_map[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "chr%s.fam", S_chr);
			Cp_fam[i] = cExporter::summon(S_ext);
		}

		/* For Undetermined IF exists */
		if (Ba_chrExists[NCHR_SPECIES]) {
			Cp_lgen[NCHR_SPECIES] = cExporter::summon("chrUn.lgen");
			Cp_fam[NCHR_SPECIES] = cExporter::summon("chrUn.fam");
			Cp_map[NCHR_SPECIES] = cExporter::summon("chrUn.map");
		}

		LOGoutput("Final dataset is exported with [%s.chr***.lgen] [%s.chr***.map] [%s.chr***.fam]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
	} else {
/**/	Cp_lgen[0]	= cExporter::summon("lgen"); /* CHECKED */
/**/	Cp_map[0]	= cExporter::summon("map"); /* CHECKED */
/**/	Cp_fam[0]	= cExporter::summon("fam"); /* CHECKED */

		LOGoutput("Final dataset is exported with [%s.lgen] [%s.map] [%s.fam]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
	}

	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

	for (vSampPtr_it i=Xa_sampleV2.begin() ; i!=Xa_sampleV2.end() ; ) {
		xSample *Xp_currSamp = *i;

		/* Determining phenotype */
		char	S_pheno[128];
		getPhenoCode(Xp_currSamp->N_idx, S_pheno, Sp_case, Sp_ctrl);

		/* For phenotype, export header */
		for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_lgen[j])
			Cp_fam[j]->fmt("%s %s %s %s %d %s\n", Xp_currSamp->S_FID.c_str(),
				Xp_currSamp->S_IID.c_str(),
				Xp_currSamp->Xp_pat?Xp_currSamp->Xp_pat->S_IID.c_str():"0",
				Xp_currSamp->Xp_mat?Xp_currSamp->Xp_mat->S_IID.c_str():"0",
				Xp_currSamp->N_sex, S_pheno);

		/* For all variants */
		if (!OPT_ENABLED(indel)) for (wsUint j=0 ; j<N_variant ; j++) {
			/* Determine target */
			wsUint N_target = 0;
			/* Target could be not 0 ONLY IF --split */
			if (OPT_ENABLED(split)) {
				if (Xv_variant[j].chr <=0 || (wsUint)Xv_variant[j].chr > NCHR_SPECIES)
					N_target = NCHR_SPECIES;
				else N_target = Xv_variant[j].chr - 1;
			}
			cExporter *Cp_currPED = Cp_lgen[N_target];

			Cp_currPED->fmt("%s	%s	%s", Xp_currSamp->S_FID.c_str(),
				Xp_currSamp->S_IID.c_str(), Xv_variant[j].name);
			switch (Na_geno[Xp_currSamp->N_idx][j]) {
			case 0: /* Major homo */
				Cp_currPED->fmt("	%c %c\n", Xv_variant[j].al1, Xv_variant[j].al1); break;
			case 1: /* Hetero */
				Cp_currPED->fmt("	%c %c\n", Xv_variant[j].al1, Xv_variant[j].al2); break;
			case 2: /* Minor homo */
				Cp_currPED->fmt("	%c %c\n", Xv_variant[j].al2, Xv_variant[j].al2); break;
			case WISARD_NA:
				Cp_currPED->fmt("	%s %s\n", S_misgeno, S_misgeno); break;
			default:
				halt_fmt(WISARD_SYST_INVL_SNPDATA, Xv_variant[j].name, Xp_currSamp->S_FID.c_str(),
					Xp_currSamp->S_IID.c_str(), Na_geno[Xp_currSamp->N_idx][j]);
				break;
			}
		} else for (wsUint j=0 ; j<N_variant ; j++) {
			/* Determine target */
			wsUint N_target = 0;
			/* Target could be not 0 ONLY IF --split */
			if (OPT_ENABLED(split)) {
				if (Xv_variant[j].chr <=0 || (wsUint)Xv_variant[j].chr > NCHR_SPECIES)
					N_target = NCHR_SPECIES;
				else N_target = Xv_variant[j].chr - 1;
			}
			cExporter *Cp_currPED = Cp_lgen[N_target];

			Cp_currPED->fmt("%s	%s	%s", Xp_currSamp->S_FID.c_str(),
				Xp_currSamp->S_IID.c_str(), Xv_variant[j].name);
			switch (Na_geno[Xp_currSamp->N_idx][j]) {
			case 0: /* Major homo */
				Cp_currPED->fmt("	%s %s\n", Xv_variant[j].indel1, Xv_variant[j].indel1); break;
			case 1: /* Hetero */
				Cp_currPED->fmt("	%s %s\n", Xv_variant[j].indel1, Xv_variant[j].indel2); break;
			case 2: /* Minor homo */
				Cp_currPED->fmt("	%s %s\n", Xv_variant[j].indel1, Xv_variant[j].indel2); break;
			case WISARD_NA:
				Cp_currPED->fmt("	%s %s\n", S_misgeno, S_misgeno); break;
			default:
				halt_fmt(WISARD_SYST_INVL_SNPDATA, Xv_variant[j].name, Xp_currSamp->S_FID.c_str(),
					Xp_currSamp->S_IID.c_str(), Na_geno[Xp_currSamp->N_idx][j]);
				break;
			}
		}
		i++;
	}

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	/* Export MAP file */
	for (vVariant_it i=Xv_variant.begin() ; i!=Xv_variant.end() ; ) {
		/* Determine target */
		wsUint N_target = 0;
		/* Target could be not 0 ONLY IF --split */
		if (OPT_ENABLED(split)) {
			if (i->chr <=0 || (wsUint)i->chr > NCHR_SPECIES)
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
	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) if (Cp_lgen[i]) {
		delete Cp_map[i];
		delete Cp_fam[i];
		delete Cp_lgen[i];
	}
	DEALLOC(Cp_lgen);
	DEALLOC(Cp_map);
	DEALLOC(Cp_fam);
}

void cIO::exportGEN()
{
	/* Allocate EXPORTERS */
	cExporter**	Cp_gen	= NULL;
	cExporter**	Cp_smp	= NULL;
	wsCalloc(Cp_gen, cExporter*, MAX_NCHR+1);
	wsCalloc(Cp_smp, cExporter*, MAX_NCHR+1);

	/* Cannot apply --outmisgeno */
	if (IS_ASSIGNED(outmisgeno))
		LOGwarn("GEN conversion does not affected by --outmisgeno\n");

	/* Cannot apply --outphenoonly */
	if (OPT_ENABLED(outphenoonly))
		LOGwarn("GEN file does not support --outphenoonly, will be ignored\n");

	if (OPT_ENABLED(split)) {
		/* Make instance of EXPOTER ONLY IF exists */
		char S_ext[64];
		for (wsUint i=0 ; i<NCHR_SPECIES ; i++) if (Ba_chrExists[i]) {
			wsStrCst S_chr = getChrName2(i + 1);
			sprintf(S_ext, "chr%s.gen", S_chr);
			Cp_gen[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "chr%s.sample", S_chr);
			Cp_smp[i] = cExporter::summon(S_ext);
		}

		/* For Undetermined IF exists */
		if (Ba_chrExists[NCHR_SPECIES]) {
			Cp_gen[NCHR_SPECIES] = cExporter::summon("chrUn.gen");
			Cp_smp[NCHR_SPECIES] = cExporter::summon("chrUn.sample");
		}

		LOGoutput("Final dataset is exported with [%s.chr***.gen] [%s.chr***.sample]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
	} else {
/**/	Cp_gen[0]	= cExporter::summon("gen"); /* CHECKED */
/**/	Cp_smp[0]	= cExporter::summon("sample"); /* CHECKED */

		LOGoutput("Final dataset is exported with [%s.gen] [%s.sample]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
	}

	for (wsUint i=0 ; i<N_variant ; i++) {
		/* Determine target */
		wsUint N_target = 0;
		/* Target could be not 0 ONLY IF --split */
		if (OPT_ENABLED(split)) {
			if (Xv_variant[i].chr <=0 || (wsUint)Xv_variant[i].chr > NCHR_SPECIES)
				N_target = NCHR_SPECIES;
			else N_target = Xv_variant[i].chr - 1;
		}
		cExporter*	Cp_currGEN	= Cp_gen[N_target];
		xVariant&	X_mk		= Xv_variant[i];

		char S_alBuf[256];
		if (OPT_ENABLED(indel))
			sprintf(S_alBuf, "%s %s", X_mk.indel1, X_mk.indel2 ?
				X_mk.indel2 : "<NA>");
// 		else if (!X_mk.al2)
// 			sprintf(S_alBuf, "%c <NA>", X_mk.al1);
		else if (!X_mk.al2) // Change to equal to GTOOL
			sprintf(S_alBuf, "%c %c", X_mk.al1, X_mk.al1);
		else
			sprintf(S_alBuf, "%c %c", X_mk.al1, X_mk.al2);

		/* Export marker-header */
// 		Cp_currGEN->fmt("SNP%d %s %d %s", i+1, X_mk.name, X_mk.pos,
// 			S_alBuf);
		Cp_currGEN->fmt("0 %s %d %s", X_mk.name, X_mk.pos, S_alBuf);

		/* Export genotype (hard-coded) */
		for (wsUint j=0 ; j<N_sample ; j++) {
			char N_geno = Na_geno[j][i];
			switch (N_geno) {
			case 0: Cp_currGEN->put(" 1 0 0"); break;
			case 1: Cp_currGEN->put(" 0 1 0"); break;
			case 2: Cp_currGEN->put(" 0 0 1"); break;
			default: /* NA */
				Cp_currGEN->put(" 0 0 0"); break;
			}
		}

		Cp_currGEN->put("\n");
	}

	/* Export sample header */
	for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_smp[j]) {
		Cp_smp[j]->put("ID_1 ID_2 missing sex");
		/* Export phenotypes */
		vPheno& Xv_phe = getPhenoInfo();
		FOREACH (vPheno_it, Xv_phe, k)
			Cp_smp[j]->fmt(" %s", k->S_name.c_str());
		/* Export covariates */
		vCovar& Xv_cov = getCovInfo();
		FOREACH (vCovar_it, Xv_cov, k)
			Cp_smp[j]->fmt(" %s", k->Sp_varName);
		Cp_smp[j]->put("\n0 0 0 D");
		/* Export phenotype type (P or B) */
		for (wsUint i=0 ; i<N_pheno ; i++)
			Cp_smp[j]->put(isContinuous(i)?" P":" B");
		/* Export covariate type (all C...at this time) */
		for (wsUint i=0 ; i<N_cov ; i++)
			Cp_smp[j]->put(" C");
		Cp_smp[j]->put("\n");
	}

	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

	/* Export sample file */
	wsUint I = 0;
	FOREACHDO (vSampPtr_it, Xa_sampleV2, i, I++) {
		/* For phenotype, export header */
		for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_smp[j]) {
			Cp_smp[j]->fmt("%s %s %g %d", (*i)->S_FID.c_str(),
				(*i)->S_IID.c_str(), 0.0, (*i)->N_sex);
			/* FIXME : 0.0 should be replaced into missing rate... */

			/* Export phenotype value */
			for (wsUint ii=0 ; ii<N_pheno ; ii++) {
				char	S_pheno[128];
				getPhenoCode(I, S_pheno, Sp_case, Sp_ctrl, ii);

				Cp_smp[j]->fmt(" %s", S_pheno);
			}
			/* Export covariate value */
			LOOP (ii, N_cov) Cp_smp[j]->fmt(" %g", Ra_cov[ii][I]);
			Cp_smp[j]->put("\n");
		}
	}

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	/* Clear */
	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) if (Cp_gen[i]) {
		delete Cp_smp[i];
		delete Cp_gen[i];
	}
	DEALLOC(Cp_gen);
	DEALLOC(Cp_smp);
}

inline wsUint getDigit(wsUint N)
{
	if (N < 10) return 1;
	else if (N < 100) return 2;
	else if (N < 1000) return 3;
	else if (N < 10000) return 4;
	else if (N < 100000) return 5;
	else if (N < 1000000) return 6;
	else if (N < 10000000) return 7;
	else if (N < 100000000) return 8;
	else if (N < 1000000000) return 9;
	else return 10;
}

/* Note : BGEN is coded as 'little-endian' */
void cIO::exportBGEN(char B_compress/*=0*/)
{
#ifndef USE_GZ
	if (B_compress)
		halt("--zipbgen is not supported with this version");
#endif
	/* Allocate EXPORTERS */
	cExporter**	Cp_bgen	= NULL;
	cExporter**	Cp_smp	= NULL;
	wsCalloc(Cp_bgen, cExporter*, MAX_NCHR+1);
	wsCalloc(Cp_smp, cExporter*, MAX_NCHR+1);

	/* Cannot apply --outphenoonly */
	if (OPT_ENABLED(outphenoonly))
		LOGwarn("BED file does not support --outphenoonly, will be ignored\n");

	if (OPT_ENABLED(split)) {
		/* Make instance of EXPOTER ONLY IF exists */
		char S_ext[64];
		for (wsUint i=0 ; i<NCHR_SPECIES ; i++) if (Ba_chrExists[i]) {
			wsStrCst S_chr = getChrName2(i + 1);
			sprintf(S_ext, "chr%s.bgen", S_chr);
			Cp_bgen[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "chr%s.sample", S_chr);
			Cp_smp[i] = cExporter::summon(S_ext);
		}

		/* For Undetermined IF exists */
		if (Ba_chrExists[NCHR_SPECIES]) {
			Cp_bgen[NCHR_SPECIES] = cExporter::summon("chrUn.bgen");
			Cp_smp[NCHR_SPECIES] = cExporter::summon("chrUn.sample");
		}

		LOGoutput("Final dataset is exported with [%s.chr***.bgen] [%s.chr***.sample]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
	} else {
/**/	Cp_bgen[0]	= cExporter::summon("bgen"); /* CHECKED */
/**/	Cp_smp[0]	= cExporter::summon("sample"); /* CHECKED */

		LOGoutput("Final dataset is exported with [%s.bgen] [%s.sample]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
	}

	/* Export BGEN header */
	for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_bgen[j]) {
		wsUint N_header = 20;
		/*
		 * The first four bytes of the file encode an unsigned integer indicating the offset, relative to the 5th byte of the file, of the start of the first snp block (or the end of the file if there are 0 snp blocks). For example, if this offset is 0, the snp blocks start at byte 5.
		 * (Since header is size 20(fixed), it should be 20
		 */
		Cp_bgen[j]->write(&N_header, 4);

		/* type   name  desc
		 * 4_UINT H     the length, in bytes, of the header block
		 * 4_UINT <NA>  the number of snp blocks stored in the file
		 * 4_UINT <NA>  the number of samples represented in the snp blocks in the file
		 * 4_<NA> <NA>  Reserved. (Writers should write 0 here, readers should ignore these bytes.)
		 * H-20_? <NA>  Free data area. This could be used to store, for example, identifying information about the file
		 * 4_<NA> flags bits numbered as for an unsigned integer. See below for flag definitions.
		 */
		wsUint N_zero = 0;
		Cp_bgen[j]->write(&N_header, 4);
		Cp_bgen[j]->write(&N_variant, 4);
		Cp_bgen[j]->write(&N_sample, 4);
		Cp_bgen[j]->write(&N_zero, 4);
		/* 131219 No free data will be written at this time */

		/* Flag
		 * bit 0 : CompressedvariantBlocks(set to 0)
		 * bit 3 : LongIds (--indel=1, otherwise 0) */
		wsUint N_flag = 0;
		if (B_compress) N_flag |= 1;
		if (OPT_ENABLED(indel)) N_flag |= 1<<3;
		Cp_bgen[j]->write(&N_flag, 4);
	}

	/* Write variant block */
	for (wsUint i=0 ; i<N_variant ; i++) {
		/* Determine target */
		wsUint N_target = 0;
		/* Target could be not 0 ONLY IF --split */
		if (OPT_ENABLED(split)) {
			if (Xv_variant[i].chr <=0 || (wsUint)Xv_variant[i].chr > NCHR_SPECIES)
				N_target = NCHR_SPECIES;
			else N_target = Xv_variant[i].chr - 1;
		}
		cExporter*	Cp_currGEN	= Cp_bgen[N_target];
		xVariant&	X_mk		= Xv_variant[i];

		/* Write variant block */

		/* Number of inds (4B) -> Fixed to N_sample */
		Cp_currGEN->write(&N_sample, 4);

		if (OPT_ENABLED(indel)) {
			/* Length of variantID = FIXED as 3(variant) + Ndigit */
			USHORT_t N_digit = (USHORT_t)getDigit(i+1) + 3;
			Cp_currGEN->write(&N_digit, 2);
			/* variantID */
			char S_vrtID[256];
			sprintf(S_vrtID, "VARIANT%d", i+1);
			Cp_currGEN->write(S_vrtID, N_digit);
			/* Length of rsID */
			unsigned short N_rsID = (unsigned short)strlen(X_mk.name);
			Cp_currGEN->write(&N_rsID, 2);
			/* rsID */
			Cp_currGEN->write(X_mk.name, N_rsID);
			/* ChrName length */
			wsStrCst S_chr = getChrName2(X_mk.chr);
			unsigned short N_chr = (unsigned short)strlen(X_mk.name);
			Cp_currGEN->write(&N_chr, 2);
			/* ChrName */
			Cp_currGEN->write(S_chr, N_chr);
			/* Variant position */
			Cp_currGEN->write(&X_mk.pos, 4);
			/* Allele A len */
			wsUint N_al1 = (wsUint)strlen(X_mk.indel1);
			Cp_currGEN->write(&N_al1, 4);
			/* Allele A */
			Cp_currGEN->write(X_mk.indel1, N_al1);
			/* Allele B len */
			wsUint N_al2 = X_mk.indel2 ? (wsUint)strlen(X_mk.indel2) : 0;
			Cp_currGEN->write(&N_al2, 4);
			/* Allele B */
			if (X_mk.indel2) Cp_currGEN->write(X_mk.indel2, N_al2);
		} else {
			/* Length of variantID  */
			unsigned char N_digit = (unsigned char)getDigit(i+1) + 3;
			/* vrtID */
			char S_vrtID[256] = { 0, };
			sprintf(S_vrtID, "VARIANT%d", i+1);
			/* Length of rsID */
			unsigned char N_rsID = (unsigned char)strlen(X_mk.name);

			/* Determine S = which one is greater? */
			unsigned char S = N_digit > N_rsID ? N_digit : N_rsID;

			/* Write down S */
			Cp_currGEN->write(&S, 1);

			/* Write down variantID */
			Cp_currGEN->write(S_vrtID, S);
			/* Write down RSID */
			char S_rsID[256] = { 0, };
			strcpy(S_rsID, X_mk.name);
			Cp_currGEN->write(S_rsID, S);

			/* Write down chromosome */
			unsigned char N_chr = 255; /* Unknown */
			if (isAutosome(X_mk))
				N_chr = (unsigned char)X_mk.chr;
			else if (isXChromosome(X_mk))
				N_chr = 23;
			else if (isYChromosome(X_mk))
				N_chr = 24;
			Cp_currGEN->write(&N_chr, 1);
			/* variant position */
			Cp_currGEN->write(&X_mk.pos, 4);
			/* Allele A */
			Cp_currGEN->write(&X_mk.al1, 1);
			/* Allele B */
			Cp_currGEN->write(&X_mk.al2, 1);
		}

		/* Write down probability... */
		if (B_compress) {
			unsigned short *Na_buf = NULL;
			wsCalloc(Na_buf, unsigned short, N_sample*3);

			for (wsUint j=0 ; j<N_sample ; j++) {
				char N_geno = Na_geno[j][i];
				unsigned short *p = Na_buf + (j*3);

				switch (N_geno) {
				case 0: p[0] = 10000; break;
				case 1: p[1] = 10000; break;
				case 2: p[2] = 10000; break;
				default: break;
				}
			}

#ifdef USE_GZ
			unsigned char *Na_cbuf = NULL;
			wsCalloc(Na_cbuf, unsigned char, N_sample*6);
			uLongf N_len = 0;
			compress((Bytef *)Na_cbuf, &N_len, (Bytef *)Na_buf, N_sample*3*sizeof(unsigned short));
			/* Write down compressed size in 4 byte */
			wsUint N_zipLen = (wsUint)N_len;
			Cp_currGEN->write(&N_zipLen, 4);
			/* Write down zipped contents */
			Cp_currGEN->write(Na_cbuf, N_zipLen);
			DEALLOC(Na_cbuf);
#endif
			DEALLOC(Na_buf);
		} else for (wsUint j=0 ; j<N_sample ; j++) {
			char N_geno = Na_geno[j][i];
			unsigned short v1 = 10000;
			unsigned short v0 = 0;

			switch (N_geno) {
			case 0:
				Cp_currGEN->write(&v1, 2);
				Cp_currGEN->write(&v0, 2);
				Cp_currGEN->write(&v0, 2);
				break;
			case 1:
				Cp_currGEN->write(&v0, 2);
				Cp_currGEN->write(&v1, 2);
				Cp_currGEN->write(&v0, 2);
				break;
			case 2:
				Cp_currGEN->write(&v0, 2);
				Cp_currGEN->write(&v0, 2);
				Cp_currGEN->write(&v1, 2);
				break;
			default:
				Cp_currGEN->write(&v0, 2);
				Cp_currGEN->write(&v0, 2);
				Cp_currGEN->write(&v0, 2);
				break;
			}
		}
	}

	/* Export sample header */
	for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_smp[j]) {
		Cp_smp[j]->put("ID_1 ID_2 missing");
		/* Export phenotypes */
		vPheno& Xv_phe = getPhenoInfo();
		FOREACH (vPheno_it, Xv_phe, k)
			Cp_smp[j]->fmt(" %s", k->S_name.c_str());
		/* Export covariates */
		vCovar& Xv_cov = getCovInfo();
		FOREACH (vCovar_it, Xv_cov, k)
			Cp_smp[j]->fmt(" %s", k->Sp_varName);
		Cp_smp[j]->put("\n0 0 0");
		/* Export phenotype type (P or B) */
		for (wsUint i=0 ; i<N_pheno ; i++)
			Cp_smp[j]->put(isContinuous(i)?" P":" B");
		/* Export covariate type (all C...at this time) */
		for (wsUint i=0 ; i<N_cov ; i++)
			Cp_smp[j]->put(" C");
		Cp_smp[j]->put("\n");
	}

	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

	/* Export sample file */
	wsUint I = 0;
	FOREACHDO (vSampPtr_it, Xa_sampleV2, i, I++) {
		/* For phenotype, export header */
		for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_smp[j]) {
			Cp_smp[j]->fmt("%s %s %g", (*i)->S_FID.c_str(),
				(*i)->S_IID.c_str(), 0.0);
			/* FIXME : 0.0 should be replaced into missing rate... */

			/* Export phenotype value */
			for (wsUint ii=0 ; ii<N_pheno ; ii++) {
				char	S_pheno[128];
				getPhenoCode(I, S_pheno, Sp_case, Sp_ctrl, ii);

				Cp_smp[j]->fmt(" %s", S_pheno);
			}
			/* Export covariate value */
			LOOP (ii, N_cov) Cp_smp[j]->fmt(" %g", Ra_cov[ii][I]);
			Cp_smp[j]->put("\n");
		}
	}

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	/* Clear */
	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) if (Cp_bgen[i]) {
		delete Cp_smp[i];
		delete Cp_bgen[i];
	}
	DEALLOC(Cp_bgen);
	DEALLOC(Cp_smp);
}

void cIO::exportBeagle()
{
	/* --outmisgeno : Beagle '0' */
	const char*	S_misgeno = IS_ASSIGNED(outmisgeno) ?
		OPT_STRING(outmisgeno) : "0";

	/* Export Beagle data */
	vSampPtr	&Xa_samp = getSample();
	//REAL_c		*Ra_pheno = getPheno();

	/* Cannot apply --outphenoonly */
	if (OPT_ENABLED(outphenoonly))
		LOGwarn("Beagle file does not support --outphenoonly, will be ignored\n");

	/* Allocate EXPORTERS */
	cExporter**	Cp_bgl = NULL;
	cExporter**	Cp_mkr = NULL;
	wsCalloc(Cp_bgl, cExporter*, MAX_NCHR+1);
	wsCalloc(Cp_mkr, cExporter*, MAX_NCHR+1);

	if (OPT_ENABLED(split)) {
		/* Make instance of EXPOTER ONLY IF exists */
		char S_ext[64];
		for (wsUint i=0 ; i<NCHR_SPECIES ; i++) if (Ba_chrExists[i]) {
			wsStrCst S_chr = getChrName2(i + 1);
			sprintf(S_ext, "chr%s.bgl", S_chr);
			Cp_bgl[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "chr%s.markers", S_chr);
			Cp_mkr[i] = cExporter::summon(S_ext);
		}

		/* For Undetermined IF exists */
		if (Ba_chrExists[NCHR_SPECIES]) {
			Cp_bgl[NCHR_SPECIES] = cExporter::summon("chrUn.bgl");
			Cp_mkr[NCHR_SPECIES] = cExporter::summon("chrUn.markers");
		}

		LOGoutput("Final dataset is exported with [%s.chr***.bgl] [%s.chr***.markers]\n",
			OPT_STRING(out), OPT_STRING(out));
	} else {
		/* Otherwise, use just one */
		Cp_bgl[0] = cExporter::summon("bgl"); /* CHECKED */
		Cp_mkr[0] = cExporter::summon("markers"); /* CHECKED */
		LOGoutput("Final dataset is exported with [%s.bgl] [%s.markers]\n",
			OPT_STRING(out), OPT_STRING(out));
	}

	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

	/* Export common part */
	for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_bgl[j])
		Cp_bgl[j]->put("I id");

	/* Export IID */
	wsUint I=0;
	FOREACHDO (vSampPtr_it, Xa_samp, i, I++) {
		/* Determining phenotype */
		char	S_pheno[128];
		getPhenoCode(I, S_pheno, Sp_case, Sp_ctrl);

		/* For phenotype, export header */
		for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_bgl[j])
			Cp_bgl[j]->fmt(" %s", (*i)->S_IID.c_str());
	}

	for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_bgl[j])
		Cp_bgl[j]->put("\nA phenotype");

	/* Export phenotype */
	I=0;
	FOREACHDO (vSampPtr_it, Xa_samp, i, I++) {
		/* Determining phenotype */
		char	S_pheno[128];
		getPhenoCode(I, S_pheno, Sp_case, Sp_ctrl);

		/* For phenotype, export header */
		for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_bgl[j])
			Cp_bgl[j]->fmt(" %s", S_pheno);
	}
	for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_bgl[j])
		Cp_bgl[j]->put("\n");

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	// 	LOGoutput("Final dataset is exported with [%s.tped] [%s.tfam]\n",
	// 		OPT_STRING(out), OPT_STRING(out));

	/* Export TPED data */
	char**		Na_data	= getGenotype();
	vVariant&	Xa_snvs	= getVariant();

	I = 0;
	FOREACHDO (vVariant_it, Xa_snvs, i, I++) {
		/* Determine target */
		wsUint N_target = 0;
		/* Target could be not 0 ONLY IF --split */
		if (OPT_ENABLED(split)) {
			if (i->chr <=0 || (wsUint)i->chr > NCHR_SPECIES)
				N_target = NCHR_SPECIES;
			else N_target = i->chr - 1;
		}
		cExporter *Cp_currBGL = Cp_bgl[N_target];
		cExporter *Cp_currMKR = Cp_mkr[N_target];

		/* Print variant information */
		Cp_currBGL->fmt("M %s", i->name);

		/* Export markers */
		if (OPT_ENABLED(indel))
			Cp_currMKR->fmt("%s %d %s %s\n", i->name, i->pos,
				i->indel1, i->indel2);
		else
			Cp_currMKR->fmt("%s %d %c %c\n", i->name, i->pos,
				i->al1, i->al2);

		/* Print variant data */
		FOREACH (vSampPtr_it, Xa_samp, j) {
			if (!OPT_ENABLED(indel)) switch (Na_data[(*j)->N_idx][I]) {
			case 0: /* Major homo */
				Cp_currBGL->fmt(" %c %c", i->al1, i->al1); break;
			case 1: /* Hetero */
				Cp_currBGL->fmt(" %c %c", i->al1, i->al2); break;
			case 2: /* Minor homo */
				Cp_currBGL->fmt(" %c %c", i->al2, i->al2); break;
			case WISARD_NA:
				Cp_currBGL->fmt(" %s %s", S_misgeno, S_misgeno); break;
			default:
				halt_fmt(WISARD_SYST_INVL_SNPDATA, i->name, (*j)->S_FID.c_str(),
					(*j)->S_IID.c_str(), Na_data[(*j)->N_idx][I]);
			} else switch (Na_data[(*j)->N_idx][I]) {
			case 0: /* Major homo */
				Cp_currBGL->fmt(" %s %s", i->indel1, i->indel1); break;
			case 1: /* Hetero */
				Cp_currBGL->fmt(" %s %s", i->indel1, i->indel2); break;
			case 2: /* Minor homo */
				Cp_currBGL->fmt(" %s %s", i->indel1, i->indel2); break;
			case WISARD_NA:
				Cp_currBGL->fmt(" %s %s", S_misgeno, S_misgeno); break;
			default:
				halt_fmt(WISARD_SYST_INVL_SNPDATA, i->name, (*j)->S_FID.c_str(),
					(*j)->S_IID.c_str(), Na_data[(*j)->N_idx][I]);
			}
		}

		Cp_currBGL->put("\n");
	}
	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) if (Cp_bgl[i]) {
		delete Cp_bgl[i];
		delete Cp_mkr[i];
	}
	DEALLOC(Cp_bgl);
	DEALLOC(Cp_mkr);
}

void cIO::exportMerlin()
{
	/* --outmisgeno : Merlin '0' */
	const char*	S_misgeno = IS_ASSIGNED(outmisgeno) ?
		OPT_STRING(outmisgeno) : "0";

	/* Allocate EXPORTERS */
	cExporter**	Cp_ped = NULL;
	cExporter**	Cp_map = NULL;
	cExporter**	Cp_dat = NULL;
	wsCalloc(Cp_ped, cExporter*, MAX_NCHR+1);
	wsCalloc(Cp_map, cExporter*, MAX_NCHR+1);
	wsCalloc(Cp_dat, cExporter*, MAX_NCHR+1);

	/* Cannot apply --outphenoonly */
	if (OPT_ENABLED(outphenoonly))
		LOGwarn("Merlin file does not support --outphenoonly, will be ignored\n");

	if (OPT_ENABLED(split)) {
		/* Make instance of EXPOTER ONLY IF exists */
		char S_ext[64];
		for (wsUint i=0 ; i<NCHR_SPECIES ; i++) if (Ba_chrExists[i]) {
			wsStrCst S_chr = getChrName2(i + 1);
			sprintf(S_ext, "chr%s.ped", S_chr);
			Cp_ped[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "chr%s.map", S_chr);
			Cp_map[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "chr%s.dat", S_chr);
			Cp_dat[i] = cExporter::summon(S_ext);
		}

		/* For Undetermined IF exists */
		if (Ba_chrExists[NCHR_SPECIES]) {
			Cp_ped[NCHR_SPECIES] = cExporter::summon("chrUn.ped");
			Cp_map[NCHR_SPECIES] = cExporter::summon("chrUn.map");
			Cp_dat[NCHR_SPECIES] = cExporter::summon("chrUn.dat");
		}

		LOGoutput("Final dataset is exported with [%s.chr***.ped] [%s.chr***.map] [%s.chr***.dat]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
	} else {
		/* Otherwise, use just one */
		Cp_ped[0] = cExporter::summon("ped"); /* CHECKED */
		Cp_map[0] = cExporter::summon("map"); /* CHECKED */
		Cp_dat[0] = cExporter::summon("dat"); /* CHECKED */
		LOGoutput("Final dataset is exported with [%s.ped] [%s.map] [%s.dat]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
	}

	/* Export common DAT fields (phenotypes) */ {
		vPheno &Xv_phe = getPhenoInfo();
		wsUint j = 0;
		FOREACHDO (vPheno_it, Xv_phe, i, j++)
			for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) if (Cp_dat[k])
				Cp_dat[k]->fmt("%c	%s\n", isContinuous(j)?'T':'A', i->S_name.c_str());
	}
	/* Export common DAT fields (covariates) */ {
		vCovar &Xv_cov = getCovInfo();
		wsUint j = 0;
		FOREACHDO (vCovar_it, Xv_cov, i, j++)
			for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) if (Cp_dat[k])
				Cp_dat[k]->fmt("C	%s\n", i->Sp_varName);
	}

	/* Export map header */
	for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) if (Cp_map[k])
		Cp_map[k]->put("CHROMOSOME	MARKER	POSITION\n");

	/* Export markers */
	for (wsUint j=0 ; j<N_variant ; j++) {
		/* Determine target */
		wsUint N_target = 0;
		/* Target could be not 0 ONLY IF --split */
		if (OPT_ENABLED(split)) {
			if (Xv_variant[j].chr <=0 || (wsUint)Xv_variant[j].chr > NCHR_SPECIES)
				N_target = NCHR_SPECIES;
			else N_target = Xv_variant[j].chr - 1;
		}
		cExporter *Cp_currDAT = Cp_dat[N_target];
		Cp_currDAT->fmt("M	%s\n", Xv_variant[j].name);
	}

	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

//	MAT_t	Ra_pheno	= getPhenos();
	wsMat	Ra_COV		= getCovariates();

	for (vSampPtr_it i=Xa_sampleV2.begin() ; i!=Xa_sampleV2.end() ; ) {
		xSample *Xp_currSamp = *i;

		/* Determining phenotype */
		char	S_pheno[128];

		/* For phenotype, export header */
		for (wsUint j=0 ; j<=NCHR_SPECIES ; j++) if (Cp_ped[j]) {
			/* Print out mandatory */
			Cp_ped[j]->fmt("%s %s %s %s %d", Xp_currSamp->S_FID.c_str(),
				Xp_currSamp->S_IID.c_str(),
				Xp_currSamp->Xp_pat?Xp_currSamp->Xp_pat->S_IID.c_str():"0",
				Xp_currSamp->Xp_mat?Xp_currSamp->Xp_mat->S_IID.c_str():"0",
				Xp_currSamp->N_sex);
			/* Print out phenotypes */
			for (wsUint k=0 ; k<N_pheno ; k++) {
				getPhenoCode(Xp_currSamp->N_idx, S_pheno, Sp_case, Sp_ctrl, k);
					Cp_ped[j]->fmt(" %s", S_pheno);
			}
			/* Print out covariates */
			for (wsUint k=0 ; k<N_cov ; k++)
				Cp_ped[j]->fmt(" %g", Ra_COV[k][Xp_currSamp->N_idx]);
		}

		/* For all variants */
		if (!OPT_ENABLED(indel)) for (wsUint j=0 ; j<N_variant ; j++) {
			/* Determine target */
			wsUint N_target = 0;
			/* Target could be not 0 ONLY IF --split */
			if (OPT_ENABLED(split)) {
				if (Xv_variant[j].chr <=0 || (wsUint)Xv_variant[j].chr > NCHR_SPECIES)
					N_target = NCHR_SPECIES;
				else N_target = Xv_variant[j].chr - 1;
			}
			cExporter *Cp_currPED = Cp_ped[N_target];

			switch (Na_geno[Xp_currSamp->N_idx][j]) {
			case 0: /* Major homo */
				Cp_currPED->fmt("  %c %c", Xv_variant[j].al1, Xv_variant[j].al1); break;
			case 1: /* Hetero */
				Cp_currPED->fmt("  %c %c", Xv_variant[j].al2, Xv_variant[j].al1); break;
			case 2: /* Minor homo */
				Cp_currPED->fmt("  %c %c", Xv_variant[j].al2, Xv_variant[j].al2); break;
			case WISARD_NA:
				Cp_currPED->fmt("  %s %s", S_misgeno, S_misgeno); break;
			default:
				halt_fmt(WISARD_SYST_INVL_SNPDATA, Xv_variant[j].name, Xp_currSamp->S_FID.c_str(),
					Xp_currSamp->S_IID.c_str(), Na_geno[Xp_currSamp->N_idx][j]);
				break;
			}
		} else for (wsUint j=0 ; j<N_variant ; j++) {
			/* Determine target */
			wsUint N_target = 0;
			/* Target could be not 0 ONLY IF --split */
			if (OPT_ENABLED(split)) {
				if (Xv_variant[j].chr <=0 || (wsUint)Xv_variant[j].chr > NCHR_SPECIES)
					N_target = NCHR_SPECIES;
				else N_target = Xv_variant[j].chr - 1;
			}
			cExporter *Cp_currPED = Cp_ped[N_target];

			switch (Na_geno[Xp_currSamp->N_idx][j]) {
			case 0: /* Major homo */
				Cp_currPED->fmt("  %s %s", Xv_variant[j].indel1, Xv_variant[j].indel1); break;
			case 1: /* Hetero */
				Cp_currPED->fmt("  %s %s", Xv_variant[j].indel1, Xv_variant[j].indel2); break;
			case 2: /* Minor homo */
				Cp_currPED->fmt("  %s %s", Xv_variant[j].indel1, Xv_variant[j].indel2); break;
			case WISARD_NA:
				Cp_currPED->fmt("  %s %s", S_misgeno, S_misgeno); break;
			default:
				halt_fmt(WISARD_SYST_INVL_SNPDATA, Xv_variant[j].name, Xp_currSamp->S_FID.c_str(),
					Xp_currSamp->S_IID.c_str(), Na_geno[Xp_currSamp->N_idx][j]);
				break;
			}
		}
		i++;
		if (i != Xa_sampleV2.end())
			for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) if (Cp_ped[i])
				Cp_ped[i]->put("\n");
	}

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	/* Export MAP file */
	for (vVariant_it i=Xv_variant.begin() ; i!=Xv_variant.end() ; ) {
		/* Determine target */
		wsUint N_target = 0;
		/* Target could be not 0 ONLY IF --split */
		if (OPT_ENABLED(split)) {
			if (i->chr <=0 || (wsUint)i->chr > NCHR_SPECIES)
				N_target = NCHR_SPECIES;
			else N_target = i->chr - 1;
		}

		Cp_map[N_target]->fmt("%d	%s	%d\n", i->chr, i->name, i->pos);
		i++;
	}

	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) if (Cp_ped[i]) {
		delete Cp_map[i];
		delete Cp_ped[i];
		delete Cp_dat[i];
	}
	DEALLOC(Cp_map);
	DEALLOC(Cp_ped);
	DEALLOC(Cp_dat);
}

void cIO::getChrwiseVariantCount(wsUint *Na_ret, char *Ba_snvIncl)
{
	for (wsUint i=0 ; i<=NCHR_SPECIES ; i++)
		Na_ret[i] = 0;

	wsUint I=0;
	FOREACHDO (vVariant_it, Xv_variant, i, I++) {
		if (Ba_snvIncl && !Ba_snvIncl[I]) continue;

		if (i->chr <=0 || (wsUint)i->chr > NCHR_SPECIES)
			Na_ret[NCHR_SPECIES]++;
		else
			Na_ret[i->chr-1]++;
	}
}

void cIO::exportBED(wsUint B_vrtMajor, char *Ba_snvIncl/*=NULL*/,
	char *Ba_sampIncl/*=NULL*/)
{
	vVariant&	Xv_vrt		= this->getVariant();
	wsUint		Na_szVrts[MAX_NCHR+1] = {0, };

	/* Cannot apply --outmisgeno */
	if (IS_ASSIGNED(outmisgeno))
		LOGwarn("BED conversion does not affected by --outmisgeno\n");

	/* Cannot apply --expression */
	if (IS_ASSIGNED(expression)) {
		LOGwarn("Expression data cannot be exported with BED format\n");
		return;
	}

	/* Cannot apply --outphenoonly */
	if (OPT_ENABLED(outphenoonly))
		LOGwarn("BED file does not support --outphenoonly, will be ignored\n");

	/* Allocate EXPORTERS */
	cExporter**	Cp_fam		= NULL;
	cExporter**	Cp_bim		= NULL;
	cExporter**	Cp_bed		= NULL;
	wsUint		N_exporter	= 0;
	/* --famsplit */
	mDataIdx	Xm_fid2idx;
	mvInt		Xm_fi2dsampIdx;

	if (OPT_ENABLED(split)) {
		N_exporter = MAX_NCHR+1;
		wsCalloc(Cp_fam, cExporter*, N_exporter);
		wsCalloc(Cp_bim, cExporter*, N_exporter);
		wsCalloc(Cp_bed, cExporter*, N_exporter);

		/* Make instance of EXPOTER ONLY IF exists */
		char S_ext[64];
		for (wsUint i=0 ; i<NCHR_SPECIES ; i++) if (Ba_chrExists[i]) {
			wsStrCst S_chr = getChrName2(i + 1);
			sprintf(S_ext, "chr%s.bed", S_chr);
			Cp_bed[i] = cExporter::summon(S_ext, 0, ET_BIN);
			sprintf(S_ext, "chr%s.fam", S_chr);
			Cp_fam[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "chr%s.bim", S_chr);
			Cp_bim[i] = cExporter::summon(S_ext);
		}

		/* For Undetermined IF exists */
		if (Ba_chrExists[NCHR_SPECIES]) {
			Cp_bed[NCHR_SPECIES] = cExporter::summon("chrUn.bed", 0, ET_BIN);
			Cp_fam[NCHR_SPECIES] = cExporter::summon("chrUn.fam");
			Cp_bim[NCHR_SPECIES] = cExporter::summon("chrUn.bim");
		}

		LOGoutput("Final dataset is exported with [%s.chr***.bed] [%s.chr***.fam] [%s.chr***.bim]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
		getChrwiseVariantCount(Na_szVrts, Ba_snvIncl);
	} else if (OPT_ENABLED(famsplit)) {
		vSampPtr& Xv_samp = getSample();
		mFam& Xm_fam = getFamilyData();
		N_exporter = (wsUint)Xm_fam.size();
		wsCalloc(Cp_fam, cExporter*, N_exporter);
		wsCalloc(Cp_bim, cExporter*, N_exporter);
		wsCalloc(Cp_bed, cExporter*, N_exporter);

		/* For each fid, create index mapping */
		wsUint i = 0;
		FOREACH (mFam_it, Xm_fam, X_fam) {
			string S_fid = X_fam->first;
			Xm_fid2idx[S_fid] = i++;
		}
		/* For each sample, create sequence of access by FID */
		i = 0;
		FOREACHDO (vSampPtr_it, Xv_samp, j, i++) {
			if (Ba_sampIncl && !Ba_sampIncl[i]) continue;
			Xm_fi2dsampIdx[(*j)->S_FID].push_back(i);
		}

		/* Make instances */
		char S_ext[64];
		i = 0;
		FOREACH (mFam_it, Xm_fam, X_fam) {
			string S_fid = X_fam->first;
			sprintf(S_ext, "fam.%s.bed", S_fid.c_str());
			Cp_bed[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "fam.%s.fam", S_fid.c_str());
			Cp_fam[i] = cExporter::summon(S_ext);
			sprintf(S_ext, "fam.%s.bim", S_fid.c_str());
			Cp_bim[i] = cExporter::summon(S_ext);

			i++;
		}

		LOGoutput("Final dataset is exported with [%s.fam.***.bed] [%s.fam.***.fam] [%s.fam.***.bim]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
	} else {
		N_exporter = 1;
		wsCalloc(Cp_fam, cExporter*, MAX_NCHR + 1);
		wsCalloc(Cp_bim, cExporter*, MAX_NCHR + 1);
		wsCalloc(Cp_bed, cExporter*, MAX_NCHR + 1);

/**/	Cp_bed[0]	= cExporter::summon("bed", 0, ET_BIN); /* CHECKED */
/**/	Cp_fam[0]	= cExporter::summon("fam"); /* CHECKED */
/**/	Cp_bim[0]	= cExporter::summon("bim"); /* CHECKED */

		LOGoutput("Final dataset is exported with [%s.bed] [%s.fam] [%s.bim]\n",
			OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));
		Na_szVrts[0] = (wsUint)Xv_vrt.size();
	}

	/* Count samples to be exported */
	wsUint N_expSamp	= 0;
//	wsUint N_expSNP		= 0;
	if (Ba_sampIncl) {
		for (wsUint i=0 ; i<N_sample ; i++)
			if (Ba_sampIncl[i]) N_expSamp++;
	} else
		N_expSamp = N_sample;
// 	if (Ba_snvIncl) {
// 		for (wsUint i=0 ; i<N_SNP ; i++)
// 			if (Ba_snvIncl[i]) N_expSNP++;
// 	} else
// 		N_expSNP = N_SNP;
	if (Xv_vrt.size() != N_variant)
		halt("SYSERR : Variant info array size != # variants");

	/* Set if outcact */
	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

	/* Parse --makeflag */
	mStr Xm_flag;
	if (IS_ASSIGNED(makeflag)) loadFlags(OPT_STRING(makeflag), Xm_flag);

	/* Export FAM file */
	wsUint j=0;
	FOREACHDO (vSampPtr_it, Xa_sampleV2, i, j++) {
		if (Ba_sampIncl && Ba_sampIncl[j] == 0) continue;
		xSample *Xp_currSamp = *i;
		//wsUint k = Xp_currSamp->N_idx;

		/* Determining phenotype */
		char	S_pheno[128];
		getPhenoCode(Xp_currSamp->N_idx, S_pheno, Sp_case, Sp_ctrl);

		/* For --famsplit, assign to their target */
		if (mapIsExist(Xm_flag, "nofid")) {
			if (OPT_ENABLED(famsplit))
				Cp_fam[Xm_fid2idx[Xp_currSamp->S_FID]]->fmt("%s %s %s %d %s\n",
					Xp_currSamp->S_IID.c_str(),
					Xp_currSamp->Xp_pat?Xp_currSamp->Xp_pat->S_IID.c_str():"0",
					Xp_currSamp->Xp_mat?Xp_currSamp->Xp_mat->S_IID.c_str():"0",
					Xp_currSamp->N_sex, S_pheno);
			/* For phenotype, export header */
			else for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) if (Cp_fam[k]) {
				Cp_fam[k]->fmt("%s %s %s %d %s\n", Xp_currSamp->S_IID.c_str(),
					Xp_currSamp->Xp_pat?Xp_currSamp->Xp_pat->S_IID.c_str():"0",
					Xp_currSamp->Xp_mat?Xp_currSamp->Xp_mat->S_IID.c_str():"0",
					Xp_currSamp->N_sex, S_pheno);
			}
		} else {
			if (OPT_ENABLED(famsplit))
				Cp_fam[Xm_fid2idx[Xp_currSamp->S_FID]]->fmt("%s %s %s %s %d %s\n", Xp_currSamp->S_FID.c_str(),
				Xp_currSamp->S_IID.c_str(),
				Xp_currSamp->Xp_pat ? Xp_currSamp->Xp_pat->S_IID.c_str() : "0",
				Xp_currSamp->Xp_mat ? Xp_currSamp->Xp_mat->S_IID.c_str() : "0",
				Xp_currSamp->N_sex, S_pheno);
			/* For phenotype, export header */
			else for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) if (Cp_fam[k]) {
				Cp_fam[k]->fmt("%s %s %s %s %d %s\n", Xp_currSamp->S_FID.c_str(),
					Xp_currSamp->S_IID.c_str(),
					Xp_currSamp->Xp_pat ? Xp_currSamp->Xp_pat->S_IID.c_str() : "0",
					Xp_currSamp->Xp_mat ? Xp_currSamp->Xp_mat->S_IID.c_str() : "0",
					Xp_currSamp->N_sex, S_pheno);
			}
		}
	}

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	/* For --famsplit, export all variant info to all files */
	if (OPT_ENABLED(famsplit)) {
		j = 0;
		FOREACHDO (vVariant_it, Xv_vrt, it, j++) {
			if (Ba_snvIncl && Ba_snvIncl[j] == 0) continue;

			for (wsUint i=0 ; i<N_exporter ; i++) {
				char Sp_chrName[64];// = getChrName(it->chr);
				sprintf(Sp_chrName, "%d", it->chr);
				if (OPT_ENABLED(indel)) {
					if (it->indel1 && it->indel2)
						Cp_bim[i]->fmt("%s %s %g %d %s %s\n", Sp_chrName, it->name, it->gdist, it->pos,
							it->indel2, it->indel1);
					else if (!it->indel2)
						Cp_bim[i]->fmt("%s %s %g %d 0 %s\n", Sp_chrName, it->name, it->gdist, it->pos,
							it->indel1);
					else
						Cp_bim[i]->fmt("%s %s %g %d 0 0\n", Sp_chrName, it->name, it->gdist, it->pos);
				} else {
					if (it->al1 && it->al2)
						Cp_bim[i]->fmt("%s %s %g %d %c %c\n", Sp_chrName, it->name, it->gdist, it->pos,
							it->al2, it->al1);
					else if (!it->al2)
						Cp_bim[i]->fmt("%s %s %g %d 0 %c\n", Sp_chrName, it->name, it->gdist, it->pos,
							it->al1);
					else
						Cp_bim[i]->fmt("%s %s %g %d 0 0\n", Sp_chrName, it->name, it->gdist, it->pos);
				}
			}
		}
	} else {
		/* Export BIM file */
		j = 0;
		FOREACHDO (vVariant_it, Xv_vrt, it, j++) {
			if (Ba_snvIncl && Ba_snvIncl[j] == 0) continue;

			/* Determine target */
			wsUint N_target = 0;
			/* Target could be not 0 ONLY IF --split */
			if (OPT_ENABLED(split)) {
				if (it->chr <= 0 || (wsUint)it->chr > NCHR_SPECIES)
					N_target = NCHR_SPECIES;
				else N_target = it->chr - 1;
			}
			cExporter *Cp_currBIM = Cp_bim[N_target];

			char Sp_chrName[64];// = getChrName(it->chr);
			sprintf(Sp_chrName, "%d", it->chr);
			if (OPT_ENABLED(indel)) {
				if (it->indel1 && it->indel2)
					Cp_currBIM->fmt("%s %s %g %d %s %s\n", Sp_chrName, it->name, it->gdist, it->pos,
						it->indel2, it->indel1);
				else if (!it->indel2)
					Cp_currBIM->fmt("%s %s %g %d 0 %s\n", Sp_chrName, it->name, it->gdist, it->pos,
						it->indel1);
				else
					Cp_currBIM->fmt("%s %s %g %d 0 0\n", Sp_chrName, it->name, it->gdist, it->pos);
			} else {
				if (it->al1 && it->al2)
					Cp_currBIM->fmt("%s %s %g %d %c %c\n", Sp_chrName, it->name, it->gdist, it->pos,
						it->al2, it->al1);
				else if (!it->al2)
					Cp_currBIM->fmt("%s %s %g %d 0 %c\n", Sp_chrName, it->name, it->gdist, it->pos,
						it->al1);
				else
					Cp_currBIM->fmt("%s %s %g %d 0 0\n", Sp_chrName, it->name, it->gdist, it->pos);
			}
			//DEALLOC(Sp_chrName);
		}
	}

	/* Export BED file */
	unsigned char C_bedSig1 = 0x6C;
	unsigned char C_bedSig2 = 0x1B; 
	char C_code[3] = { 0x3, 0x2, 0x0 };
	for (wsUint i=0 ; i<N_exporter ; i++) if (Cp_bed[i]) {
		Cp_bed[i]->write(&C_bedSig1, 1);
		Cp_bed[i]->write(&C_bedSig2, 1);
	}
	if (B_vrtMajor) {
		C_bedSig1 = 0x01;
		for (wsUint i=0 ; i<N_exporter ; i++) if (Cp_bed[i])
			Cp_bed[i]->write(&C_bedSig1, 1);

		/* 1 row = 1 variant for all samples */
// 		wsUint Na_rowInBytes[MAX_CHR+1] = {0, };
// 		for (wsUint i=0 ; i<=MAX_CHR2 ; i++)
// 			Na_rowInBytes[i] = (Na_szSNPs[i]+3)>>2;
// 		LOG("Required row size is %d bytes from %d samples\n",
// 			N_rowInBytes, N_expSamp);

		for (wsUint i=0 ; i<N_variant ; i++) {
			/* Filtering if required */
			if (Ba_snvIncl && Ba_snvIncl[i] == 0) continue;

			/* For --famsplit, export the dataset by its FID sets */
			if (OPT_ENABLED(famsplit)) {
				FOREACH (mvInt_it, Xm_fi2dsampIdx, j) {
					wsUint N_target = Xm_fid2idx[j->first];
					if (IS_ASSIGNED(dosage)) for (wsUint l=0,k=0 ; l<j->second.size() ; ) {
						unsigned char C_write = 0x00;

						for (k=0 ; k<4&&(l+k)<j->second.size() ; k++) {
								/* Filtering if required */
								if (Ba_sampIncl && Ba_sampIncl[j->second[l+k]] == 0) continue;

								wsFloat N_genotype = Ra_data[j->second[l+k]][i];
								if (!isMissingReal(N_genotype)) {
										if (N_genotype < W0 || N_genotype > W2)
												halt("Invalid genotype [%g] found", N_genotype);
										C_write |= C_code[(int)(N_genotype+0.5f)]<<(k<<1);
								} else
										C_write |= (1)<<(k<<1);
						}
						l+=k;
						Cp_bed[N_target]->write(&C_write, 1);
					} else for (wsUint l=0,k=0 ; l<j->second.size() ; ) {
						unsigned char C_write = 0x00;

						for (k=0 ; k<4&&(l+k)<j->second.size() ; k++) {
							/* Filtering if required */
							if (Ba_sampIncl && Ba_sampIncl[j->second[l+k]] == 0) continue;

							int N_genotype = Na_geno[j->second[l+k]][i];
							if (isAvailable(N_genotype)) {
								if (N_genotype < 0 || N_genotype > 2)
									halt("Invalid genotype [%d] found", N_genotype);
								C_write |= (C_code[N_genotype])<<(k<<1);
							} else
								C_write |= (1)<<(k<<1);
						}
						l+=k;
						Cp_bed[N_target]->write(&C_write, 1);
					}
				}
			} else {
				/* Determine target */
				wsUint N_target = 0;
				/* Target could be not 0 ONLY IF --split */
				if (OPT_ENABLED(split)) {
					if (Xv_vrt[i].chr <= 0 || (wsUint)Xv_vrt[i].chr > NCHR_SPECIES)
						N_target = NCHR_SPECIES;
					else N_target = Xv_vrt[i].chr - 1;
				}
				cExporter *Cp_currBED = Cp_bed[N_target];

				if (IS_ASSIGNED(dosage)) for (wsUint j=0,k=0 ; j<N_sample ; ) {
					unsigned char C_write = 0x00;

					for (k=0 ; k<4&&(j+k)<N_sample ; k++) {
							/* Filtering if required */
							if (Ba_sampIncl && Ba_sampIncl[j+k] == 0) continue;

							wsFloat N_genotype = Ra_data[j+k][i];
							if (!isMissingReal(N_genotype)) {
									if (N_genotype < W0 || N_genotype > W2)
											halt("Invalid genotype [%g] found", N_genotype);
									C_write |= C_code[(int)(N_genotype+0.5f)]<<(k<<1);
							} else
									C_write |= (1)<<(k<<1);
					}
					j+=k;
					Cp_currBED->write(&C_write, 1);
				} else for (wsUint j=0,k=0 ; j<N_sample ; ) {
					unsigned char C_write = 0x00;

					for (k=0 ; k<4&&(j+k)<N_sample ; k++) {
						/* Filtering if required */
						if (Ba_sampIncl && Ba_sampIncl[j+k] == 0) continue;

						int N_genotype = Na_geno[j+k][i];
						if (isAvailable(N_genotype)) {
							if (N_genotype < 0 || N_genotype > 2)
								halt("Invalid genotype [%d] found", N_genotype);
							C_write |= (C_code[N_genotype])<<(k<<1);
						} else
							C_write |= (1)<<(k<<1);
					}
					j+=k;
					Cp_currBED->write(&C_write, 1);
				}
			}
		}
	} else {
		/* FIXME : Sample-wise BED export does not support --split, will be ignored */
		/* FIXME : Sample-wise BED export does not support --famsplit, will be ignored */
		LOG("Sample-wise BED export does not support --split and --famsplit, will be ignored\n");

		C_bedSig1 = 0x00;
		Cp_bed[0]->write(&C_bedSig1, 1);

//		wsUint N_rowInBytes = (N_expSamp+3)>>2;
// 		LOG("Required row size is %d bytes from %d SNPs",
// 			N_rowInBytes, N_expSNP);

		for (wsUint i=0 ; i<N_sample ; i++) {
			/* Filtering if required */
			if (Ba_sampIncl && Ba_sampIncl[i] == 0) continue;

			for (wsUint j=0,k=0 ; j<N_variant ; ) {
				char C_code[3] = { 0x0, 0x1, 0x3 };
				unsigned char C_write = 0x00;

				for (k=0 ; k<4&&(j+k)<N_variant ; k++) {
					/* Filtering if required */
					if (Ba_snvIncl && Ba_snvIncl[j+k] == 0) continue;

					int N_genotype = Na_geno[i][j+k];
					C_write <<= 2;
					if (isAvailable(N_genotype))
						C_write |= C_code[N_genotype];
					else
						C_write |= 2;
				}
				j += k;
				Cp_bed[0]->write(&C_write, 1);
			}
		}
	}

	for (wsUint i=0 ; i<N_exporter ; i++) if (Cp_bed[i]) {
		delete Cp_bed[i];
		delete Cp_bim[i];
		delete Cp_fam[i];
	}
	DEALLOC(Cp_bed);
	DEALLOC(Cp_bim);
	DEALLOC(Cp_fam);
	/*** DO SOMETHING ***/
}

void cIO::exportBEDgene(wsUint B_vrtMajor, cSetManagerAnalysis *Cp_anaSM,
					char *Ba_sampIncl/*=NULL*/)
{
	vVariant&	Xv_vrt		= this->getVariant();

	/* Cannot apply --outphenoonly */
	if (OPT_ENABLED(outphenoonly))
		LOGwarn("BED file does not support --outphenoonly, will be ignored\n");

	/* Print log message */
	LOGoutput("Final dataset is exported with [%s.***.bed] [%s.***.fam] [%s.***.bim] by gene definition\n",
		OPT_STRING(out), OPT_STRING(out), OPT_STRING(out));

	/* Count samples to be exported */
	wsUint N_expSamp	= 0;
	//	wsUint N_expSNP		= 0;
	if (Ba_sampIncl) {
		for (wsUint i=0 ; i<N_sample ; i++)
			if (Ba_sampIncl[i]) N_expSamp++;
	} else
		N_expSamp = N_sample;

	// Sanity check
	if (Xv_vrt.size() != N_variant)
		halt("SYSERR : Variant info array size != # variants");

	/* Set if outcact */
	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

	mGeneDef&	Xm_gene	= Cp_anaSM->getGeneDef();
	wsUint		N_file	= (wsUint)Xm_gene.size();
	wsUint		j		= 0;
	FOREACHDO (mGeneDef_it, Xm_gene, i, j++) {
		/* Allocate EXPORTER */
		cExporter*	Cp_fam	= NULL;
		cExporter*	Cp_bim	= NULL;
		cExporter*	Cp_bed	= NULL;

		/* Get the number of available genes */
		unsigned char C_bedSig1 = 0x6C;
		unsigned char C_bedSig2 = 0x1B; 

		/* Make instance of EXPORTER for each genes */
		//LOG("[Step 1] Preparing handles...\n");
		char	S_ext[64];
		wsStrCst	S_gene = i->first.c_str();

		sprintf(S_ext, "%s.bed", S_gene);
		Cp_bed = cExporter::summon(S_ext, 0, ET_BIN);
		sprintf(S_ext, "%s.fam", S_gene);
		Cp_fam = cExporter::summon(S_ext);
		sprintf(S_ext, "%s.bim", S_gene);
		Cp_bim = cExporter::summon(S_ext);

		/* Export FAM file */
		j = 0;
		FOREACHDO (vSampPtr_it, Xa_sampleV2, k, j++) {
			if (Ba_sampIncl && Ba_sampIncl[j] == 0) continue;
			xSample *Xp_currSamp = *k;
			//wsUint k = Xp_currSamp->N_idx;

			/* Determining phenotype */
			char	S_pheno[128];
			getPhenoCode(Xp_currSamp->N_idx, S_pheno, Sp_case, Sp_ctrl);

			/* For phenotype, export header */
			Cp_fam->fmt("%s %s %s %s %d %s\n", Xp_currSamp->S_FID.c_str(),
				Xp_currSamp->S_IID.c_str(),
				Xp_currSamp->Xp_pat?Xp_currSamp->Xp_pat->S_IID.c_str():"0",
				Xp_currSamp->Xp_mat?Xp_currSamp->Xp_mat->S_IID.c_str():"0",
				Xp_currSamp->N_sex, S_pheno);
		}

		/* Export BIM file */
		vInt& Xv_idx = i->second;
		FOREACH (vInt_it, Xv_idx, k) {
			xVariant &X_var = Xv_vrt[*k];

			char Sp_chrName[64];// = getChrName(it->chr);
			sprintf(Sp_chrName, "%d", X_var.chr);
			if (X_var.indel1 && X_var.indel2)
				Cp_bim->fmt("%s %s 0 %d %s %s\n", Sp_chrName, X_var.name,
				X_var.pos, X_var.indel2, X_var.indel1);
			else if (X_var.indel1)
				Cp_bim->fmt("%s %s 0 %d 0 %s\n", Sp_chrName, X_var.name,
				X_var.pos, X_var.indel1);
			else
				Cp_bim->fmt("%s %s 0 %d %c %c\n", Sp_chrName, X_var.name,
				X_var.pos, X_var.al2, X_var.al1);
		}

		/* Export BED file */
		Cp_bed->write(&C_bedSig1, 1);
		Cp_bed->write(&C_bedSig2, 1);

		if (B_vrtMajor) {
			C_bedSig1 = 0x01;
			Cp_bed->write(&C_bedSig1, 1);

			j = 0;
			/* Determine target */
			FOREACH (vInt_it, Xv_idx, l) for (wsUint j=0,k=0 ; j<N_sample ; ) {
				char C_code[3] = { 0x3, 0x2, 0x0 };
				unsigned char C_write = 0x00;

				for (k=0 ; k<4&&(j+k)<N_sample ; k++) {
					/* Filtering if required */
					if (Ba_sampIncl && Ba_sampIncl[j+k] == 0) continue;

					int N_genotype = Na_geno[j+k][*l];
					if (isAvailable(N_genotype)) {
						if (N_genotype < 0 || N_genotype > 2)
							halt("Invalid genotype found");
						C_write |= (C_code[N_genotype])<<(k<<1);
					} else
						C_write |= (1)<<(k<<1);
				}
				j+=k;
				Cp_bed->write(&C_write, 1);
			}
		} else {
			/* FXIME : Sample-wise BED export does not support --split, will be ignored */
			halt("Still not supported");
		}
		notice("[%d/%d] genes written...\r", j, N_file);

		delete Cp_bed;
		delete Cp_bim;
		delete Cp_fam;
	}
	LOG("[%d/%d] genes written...\n", N_file, N_file);

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	/*** DO SOMETHING ***/
}

void cIO::exportRAW(int mode/*=0*/)
{
	/* --outmisgeno : RAW 'NA' */
	const char*	S_misgeno = IS_ASSIGNED(outmisgeno) ?
		OPT_STRING(outmisgeno) : "NA";
	/* --expression will ignore mode */
	if (IS_ASSIGNED(expression) && mode)
		LOGwarn("Genotype mode will be ignored due to --expression\n");

	/* Allocate EXPORTERS */
	cExporter**	Cp_raw		= NULL;
	wsUint		N_exporter	= 0;
	/* --famsplit */
	mDataIdx	Xm_fid2idx;

	if (OPT_ENABLED(split)) {
		N_exporter = MAX_NCHR+1;
		wsCalloc(Cp_raw, cExporter*, N_exporter);

		/* Make instance of EXPOTER ONLY IF exists */
		char S_ext[64];
		for (wsUint i=0 ; i<NCHR_SPECIES ; i++) if (Ba_chrExists[i]) {
			sprintf(S_ext, "chr%s.raw", getChrName2(i + 1));
			Cp_raw[i] = cExporter::summon(S_ext);
		}

		/* For Undetermined IF exists */
		if (Ba_chrExists[NCHR_SPECIES])
			Cp_raw[NCHR_SPECIES] = cExporter::summon("chrUn.raw");

		LOGoutput("Final dataset is exported with [%s.chr***.raw]\n", OPT_STRING(out));
	} else if (OPT_ENABLED(famsplit)) {
		mFam& Xm_fam = getFamilyData();
		N_exporter = (wsUint)Xm_fam.size();
		wsCalloc(Cp_raw, cExporter*, N_exporter);

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
			sprintf(S_ext, "fam.%s.raw", S_fid.c_str());
			Cp_raw[i] = cExporter::summon(S_ext);

			i++;
		}

		LOGoutput("Final dataset is exported with [%s.fam.***.raw]\n",
			OPT_STRING(out));
	} else {
		N_exporter = 1;
		wsCalloc(Cp_raw, cExporter*, NCHR_SPECIES);

/**/	Cp_raw[0]	= cExporter::summon("raw"); /* CHECKED */
		LOGoutput("Final dataset is exported with [%s.raw]\n", OPT_STRING(out));
	}

	if (!OPT_ENABLED(outnoheader)) for (wsUint k=0 ; k<N_exporter ; k++) if (Cp_raw[k]) {
		if (mode != 3) {
			if (!OPT_ENABLED(outphenoonly))
				Cp_raw[k]->put("FID IID PAT MAT SEX PHENOTYPE");
			else
				Cp_raw[k]->put("PHENOTYPE");

			for (wsUint i=0 ; i<N_variant ; i++)
				Cp_raw[k]->fmt(" %s", Xv_variant[i].name);
		} else {
			for (wsUint i=0 ; i<N_variant ; i++)
				Cp_raw[k]->fmt("%s	", Xv_variant[i].name);
			Cp_raw[k]->put("Class");
		}
		Cp_raw[k]->put("\n");
	}

	char *Sp_case = NULL;
	char *Sp_ctrl = NULL;
	getOutCaCt(&Sp_case, &Sp_ctrl);

	wsUint N_exported = 0;
	vSampPtr& Xa_samp = getSample();
	FOREACH (vSampPtr_it, Xa_samp, it) {
		xSample *Xp = *it;
		int i = Xp->N_idx;
		if (i == -1)
			halt("Sample [%s:%s] marked as nonmissing but do not have data",
				Xp->S_FID.c_str(), Xp->S_IID.c_str());
		/* Do not print samples with missing phenotype with --makemdr */
		if (mode == 3 && isMissingReal(Ra_pheno[0][i]))
			continue;

		/* Determining phenotype */
		char	S_pheno[128];
		getPhenoCode(i, S_pheno, Sp_case, Sp_ctrl);

		/* For phenotype, export header */
		if (mode != 3) {
			if (OPT_ENABLED(famsplit)) {
				if (!OPT_ENABLED(outphenoonly)) Cp_raw[Xm_fid2idx[Xp->S_FID]]->fmt("%s %s %s %s %d %s", Xp->S_FID.c_str(),
					Xp->S_IID.c_str(),
					Xp->Xp_pat?Xp->Xp_pat->S_IID.c_str():"0",
					Xp->Xp_mat?Xp->Xp_mat->S_IID.c_str():"0",
					Xp->N_sex, S_pheno);
				else Cp_raw[Xm_fid2idx[Xp->S_FID]]->fmt("%s", S_pheno);
			} else for (wsUint j=0 ; j <= NCHR_SPECIES ; j++) if (Cp_raw[j]) {
				if (!OPT_ENABLED(outphenoonly)) Cp_raw[j]->fmt("%s %s %s %s %d %s", Xp->S_FID.c_str(),
					Xp->S_IID.c_str(),
					Xp->Xp_pat?Xp->Xp_pat->S_IID.c_str():"0",
					Xp->Xp_mat?Xp->Xp_mat->S_IID.c_str():"0",
					Xp->N_sex, S_pheno);
				else Cp_raw[j]->fmt("%s", S_pheno);
			}
		}

		if (Na_geno) for (wsUint j=0 ; j<N_variant ; j++) {
			/* Determine target */
			wsUint N_target = 0;

			/* For --famsplit, allocate target by its FID */
			if (OPT_ENABLED(famsplit)) {
				N_target = Xm_fid2idx[Xp->S_FID];
			} else if (OPT_ENABLED(split)) { /* Target could be not 0 ONLY IF --split */
				if (Xv_variant[j].chr <=0 || (wsUint)Xv_variant[j].chr > NCHR_SPECIES)
					N_target = NCHR_SPECIES;
				else N_target = Xv_variant[j].chr - 1;
			}
			cExporter *Cp_currRAW = Cp_raw[N_target];

			if (isMissing(Na_geno[i][j]))
				Cp_currRAW->fmt(mode==3?"%s	":" %s", S_misgeno);
			else switch (mode) {
			case 0: Cp_currRAW->fmt(" %d", Na_geno[i][j]); break;
			case 1: /* DOM 00 0 01 1 02 1 */ Cp_currRAW->fmt(" %d", (Na_geno[i][j]+1)>>1); break;
			case 2: /* REC 00 0 01 0 02 1 */ Cp_currRAW->fmt(" %d", Na_geno[i][j]>>1); break;
			case 3: Cp_currRAW->fmt("%d	", Na_geno[i][j]); break;
			default: halt("Invalid genotype [%d] found in variant [%s] for sample [%s::%s]",
				Na_geno[i][j], Xv_variant[j].name, Xa_sampleV2[i]->S_FID.c_str(),
				Xa_sampleV2[i]->S_IID.c_str());
			}
		} else if (Ra_data) for (wsUint j=0 ; j<N_variant ; j++) {
			if (mode == 3)
				halt("--makemdr cannot be used to export non-genotype data!");
			/* Determine target */
			wsUint N_target = 0;

			/* For --famsplit, allocate target by its FID */
			if (OPT_ENABLED(famsplit)) {
				N_target = Xm_fid2idx[Xp->S_FID];
			} else if (OPT_ENABLED(split)) { /* Target could be not 0 ONLY IF --split */
				if (Xv_variant[j].chr <=0 || (wsUint)Xv_variant[j].chr > NCHR_SPECIES)
					N_target = NCHR_SPECIES;
				else N_target = Xv_variant[j].chr - 1;
			}
			cExporter *Cp_currRAW = Cp_raw[N_target];

			if (isMissingReal(Ra_data[i][j]))
				Cp_currRAW->fmt(" %s", S_misgeno);
			else switch (mode) {
			case 0: Cp_currRAW->fmt(" %g", Ra_data[i][j]); break;
			case 1: /* DOM 00 0 01 1 02 1 */ Cp_currRAW->fmt(" %d", ((int)Ra_data[i][j]+1)>>1); break;
			case 2: /* REC 00 0 01 0 02 1 */ Cp_currRAW->fmt(" %d", (int)Ra_data[i][j]>>1); break;
			default: halt("Invalid genotype [%d] found in variant [%s] for sample [%s::%s]",
				Na_geno[i][j], Xv_variant[j].name, Xa_sampleV2[i]->S_FID.c_str(),
				Xa_sampleV2[i]->S_IID.c_str());
			}
		} else {
			wsFloat** Ra_dsg = getDosage();
			if (Ra_dsg) for (wsUint j=0 ; j<N_variant ; j++) {
				if (mode == 3)
					halt("--makemdr cannot be used to export non-genotype data!");
				/* Determine target */
				wsUint N_target = 0;

				/* For --famsplit, allocate target by its FID */
				if (OPT_ENABLED(famsplit)) {
					N_target = Xm_fid2idx[Xp->S_FID];
				} else if (OPT_ENABLED(split)) { /* Target could be not 0 ONLY IF --split */
					if (Xv_variant[j].chr <=0 || (wsUint)Xv_variant[j].chr > NCHR_SPECIES)
						N_target = NCHR_SPECIES;
					else N_target = Xv_variant[j].chr - 1;
				}
				cExporter *Cp_currRAW = Cp_raw[N_target];

				if (isMissingReal(Ra_dsg[i][j]))
					Cp_currRAW->fmt(" %s", S_misgeno);
				else switch (mode) {
				case 0: Cp_currRAW->fmt(" %g", Ra_dsg[i][j]); break;
				case 1: /* DOM 00 0 01 1 02 1 */ Cp_currRAW->fmt(" %d", ((int)Ra_dsg[i][j]+1)>>1); break;
				case 2: /* REC 00 0 01 0 02 1 */ Cp_currRAW->fmt(" %d", (int)Ra_dsg[i][j]>>1); break;
				default: halt("Invalid genotype [%d] found in variant [%s] for sample [%s::%s]",
					Na_geno[i][j], Xv_variant[j].name, Xa_sampleV2[i]->S_FID.c_str(),
					Xa_sampleV2[i]->S_IID.c_str());
				}
			}
		}
		/* For --famsplit, allocate target by its FID */
		if (OPT_ENABLED(famsplit)) {
			/* w/ --makemdr, phenotype coding must be 0 or 1 */
			if (mode == 3)
				Cp_raw[Xm_fid2idx[Xp->S_FID]]->fmt("%d\n", Ra_pheno[0][i] == WISARD_AFFECTED ? 1 : 0);
			else
				Cp_raw[Xm_fid2idx[Xp->S_FID]]->put("\n");
		} else for (wsUint k=0 ; k <= NCHR_SPECIES ; k++) if (Cp_raw[k]) {
			/* w/ --makemdr, phenotype coding must be 0 or 1 */
			if (mode == 3)
				Cp_raw[k]->fmt("%d\n", Ra_pheno[0][i]==WISARD_AFFECTED?1:0);
			else
				Cp_raw[k]->put("\n");
		}

		if (((N_exported++)%10) == 9)
			notice(" %d samples exported\r", N_exported+1);
	}
	for (wsUint k=0 ; k<N_exporter ; k++) if (Cp_raw[k])
		delete Cp_raw[k];
	DEALLOC(Cp_raw);

	if (IS_ASSIGNED(outcact)) {
		DEALLOC(Sp_case);
		DEALLOC(Sp_ctrl);
	}

	LOG("[%d] samples exported\n", N_sample);
}

void cIO::exportCovariates()
{
	vSampPtr	&Xa_samples	= getSample();
	wsUintCst		N_cov		= sizeCovar();
	wsUint		i, j;
	if (N_cov == 0) {
		LOGwarn("No covariates found, --makecov have NO effect\n");
		return;
	}
/**/cExporter*	Cp_cov		= cExporter::summon("covar.mat");
	LOGoutput("Covariates with final dataset are exported to [%s.covar.mat]\n",
		OPT_STRING(out));

	Cp_cov->put("FID	IID");
	/* Export label */
	for (i=0 ; i<N_cov ; i++)
		Cp_cov->fmt("	%s", Xa_covInfo[i].Sp_varName);
	Cp_cov->put("\n");

	i=0;
	FOREACHDO (vSampPtr_it, Xa_samples, it, i++) {
		Cp_cov->fmt("%s	%s", (*it)->S_FID.c_str(), (*it)->S_IID.c_str());

		for (j=0 ; j<N_cov ; j++)
			Cp_cov->fmt("	%g", Ra_cov[j][i]);
		Cp_cov->put("\n");
	}
	delete Cp_cov;
}

void cIO::exportPhenotypes()
{
/**/cExporter*	Cp_phe		= cExporter::summon("pheno.mat");
	vSampPtr&	Xa_samples	= getSample();
	vPheno&		X_phenos	= getPhenoInfo();
	wsUintCst	N_phe		= sizePheno();
	wsUint		i, j;
	LOGoutput("Phenotypes with final dataset are exported to [%s.pheno.mat]\n",
		OPT_STRING(out));

	/* Parse --makeflag */
	mStr Xm_flag;
	if (IS_ASSIGNED(makeflag)) loadFlags(OPT_STRING(makeflag), Xm_flag);

	Cp_phe->put(mapIsExist(Xm_flag, "nofid") ? "IID" : "FID	IID");
	/* Export label */
	for (i=0 ; i<N_phe ; i++)
		Cp_phe->fmt("	%s", X_phenos[i].S_name.c_str());
	Cp_phe->put("\n");

	i=0;
	FOREACHDO (vSampPtr_it, Xa_samples, it, i++) {
		if (mapIsExist(Xm_flag, "nofid"))
			Cp_phe->fmt("%s", (*it)->S_IID.c_str());
		else
			Cp_phe->fmt("%s	%s", (*it)->S_FID.c_str(), (*it)->S_IID.c_str());

		for (j=0 ; j<N_phe ; j++)
			Cp_phe->fmt("	%g", Ra_pheno[j][i]);
		Cp_phe->put("\n");
	}
	delete Cp_phe;
}

void cIO::exportFamDiagram()
{
	/* For each family, get the extent */
	mFam &Xm_fam = getFamilyData();

	FOREACH (mFam_it, Xm_fam, i) {
		/* FIXME : For each family.... */
	}
}

void cIO::exportFamily()
{
/**/cExporter* Cp_famstr = cExporter::summon("structure.fam");
	LOGoutput("Family structure summary is exported to [%s.structure.fam]\n",
		OPT_STRING(out));

	/* For all families */
	for (mFam_it it=Xm_fid2data.begin() ; it!=Xm_fid2data.end() ; it++) {
		xFamily&	X_fam = it->second;
		Cp_famstr->fmt("Family [%s] : %lu founders\n", X_fam.S_FID.c_str(),
			X_fam.Xp_founders.size());

		/* Print		*/
		wsUint i=1;
		FOREACHDO (vSampPtr_it, X_fam.Xp_founders, sit, i++)
			Cp_famstr->fmt("	Founder %d : %s%s\n", i, (*sit)->S_IID.c_str(),
				(*sit)->B_isMissing?" (Missing)":"");

		mSampPtr Xp_currVisitor;

		/* Set initial visitors to founders */
		FOREACH(vSampPtr_it, X_fam.Xp_founders, it)
			Xp_currVisitor.insert(make_pair((*it)->S_IID, *it));

		/* While there are remained visitors */
		i = 1;
		while (Xp_currVisitor.size()) {
			mSampPtr Xp_newVisitors;

			/* Visiting all visitors */
			FOREACH(mSampPtr_it, Xp_currVisitor, it) {
				xSample *Xp_currSamp = it->second;
				if (Xp_currSamp->N_isVisited)
					continue;
				Xp_currSamp->N_isVisited = 1;

				Cp_famstr->fmt("	Member %d : %s%s%s\n", i++, Xp_currSamp->S_IID.c_str(),
					Xp_currSamp->Xp_mat||Xp_currSamp->Xp_pat?" (Nonfounder)":" (Founder)",
					Xp_currSamp->B_isMissing?" (Missing)":"");

				/* Add all children to new visitors */
				FOREACH (vSampPtr_it, Xp_currSamp->Xp_childs, sit) {
					if ((*sit)->N_isVisited == 0)
						Xp_newVisitors.insert(make_pair((*sit)->S_IID, *sit));
				}
			}

			Xp_currVisitor.clear();
			Xp_currVisitor = Xp_newVisitors;
		}
	}
	delete Cp_famstr;
}

void cIO::exportListVariant()
{
	vVariant&	Xv_vrt	= getVariant();
	cExporter*	Cp_lm	= cExporter::summon("variant.lst");
	Cp_lm->put("NAME	MAJOR	MINOR\n");

	LOGoutput("List of variants in the final dataset is exported to "
		"[%s.variant.lst]\n", OPT_STRING(out));

	FOREACH (vVariant_it, Xv_vrt, i)
		if (OPT_ENABLED(indel)) {
			if (i->indel2 == NULL)
				Cp_lm->fmt("%s	%s	<NA>\n", i->name, i->indel1);
			else
				Cp_lm->fmt("%s	%s	%s\n", i->name, i->indel1, i->indel2);
		} else {
			if (!i->al2)
				Cp_lm->fmt("%s	%c	<NA>\n", i->name, i->al2);
			else
				Cp_lm->fmt("%s	%c	%c\n", i->name, i->al1, i->al2);
		}

	delete Cp_lm;
}

void cIO::exportListSample()
{
	vSampPtr	&Xa_smp	= getSample();
	cExporter	*Cp_ls	= cExporter::summon("sample.lst");
	wsUint		N_phe	= sizePheno();
	wsRealCst		*Ra_phe	= getPheno();

	LOGoutput("List of samples in the final dataset is exported to "
		"[%s.sample.lst]\n", OPT_STRING(out));

	if (N_phe > 1)
		LOGwarn("Multiple phenotypes assigned, only the first phenotype "
			"[%s] is exported in sample list\n", getPhenoInfo()[0].S_name.c_str());
	else if (Ra_phe == NULL)
		LOGwarn("No phenotype assigned, phenotype column will be miss "
			"in exported in sample list\n");

	/* Print header */
	if (Ra_phe == NULL)
		Cp_ls->put("FID	IID	PAT	MAT	SEX\n");
	else
		Cp_ls->put("FID	IID	PAT	MAT	SEX	PHENOTYPE\n");

	FOREACH (vSampPtr_it, Xa_smp, i)
		if (Ra_phe)
			Cp_ls->fmt("%s %s %s %s %d %g\n", (*i)->S_FID.c_str(),
				(*i)->S_IID.c_str(),
				(*i)->Xp_pat?(*i)->Xp_pat->S_IID.c_str():"0",
				(*i)->Xp_mat?(*i)->Xp_mat->S_IID.c_str():"0",
				(*i)->N_sex, Ra_phe[(*i)->N_idx]);
		else
			Cp_ls->fmt("%s %s %s %s %d\n", (*i)->S_FID.c_str(),
				(*i)->S_IID.c_str(),
				(*i)->Xp_pat?(*i)->Xp_pat->S_IID.c_str():"0",
				(*i)->Xp_mat?(*i)->Xp_mat->S_IID.c_str():"0",
				(*i)->N_sex);

	delete Cp_ls;
}

void cIO::exportListFounder()
{
	vSampPtr	&Xa_smp	= getSample();
	cExporter	*Cp_ls	= cExporter::summon("founder.lst");
	wsUint		N_phe	= sizePheno();
	wsRealCst		*Ra_phe	= getPheno();

	LOGoutput("List of founders in the final dataset is exported to "
		"[%s.founder.lst]\n", OPT_STRING(out));

	if (N_phe > 1)
		LOGwarn("Multiple phenotypes assigned, only the first phenotype "
		"[%s] is exported in founder list\n", getPhenoInfo()[0].S_name.c_str());
	else if (Ra_phe == NULL)
		LOGwarn("No phenotype assigned, phenotype column will be miss "
		"in exported in founder list\n");

	/* Print header */
	if (Ra_phe == NULL)
		Cp_ls->put("FID	IID	SEX\n");
	else
		Cp_ls->put("FID	IID	SEX	PHENOTYPE\n");

	FOREACH (vSampPtr_it, Xa_smp, i) {
		if ((*i)->Xp_mat || (*i)->Xp_pat) continue;

		if (Ra_phe)
			Cp_ls->fmt("%s %s %d %g\n", (*i)->S_FID.c_str(),
			(*i)->S_IID.c_str(),
			(*i)->N_sex, Ra_phe[(*i)->N_idx]);
		else
			Cp_ls->fmt("%s %s %d\n", (*i)->S_FID.c_str(),
			(*i)->S_IID.c_str(),
			(*i)->N_sex);
	}

	delete Cp_ls;
}

wsUint loadSampleList(wsStrCst S_fn, wsStrCst S_desc, eStr& Xm_set)
{
	wsUint		N_nonRegItem = 0;
	cStrFile	C_sampNA(S_fn, S_desc, 1);

	if (C_sampNA.isFailed()) {
		LOG("Failed to open [%s], maybe the IID selection option is given as their FID/IID pairs?\n",
			S_fn);
		char *S_remsamp = NULL;
		wsAlloc(S_remsamp, char, strlen(S_fn) + 1);
		strcpy(S_remsamp, S_fn);

		char *S_tok = strtok(S_remsamp, ",");
		while (S_tok) {
			string S_key = S_tok;
			Xm_set.insert(S_key);
			S_tok = strtok(NULL, ",");
		}
		Xm_set.insert("__LIST__");
		N_nonRegItem++;

		DEALLOC(S_remsamp);
	} else {
		int		L = 1;
		char	*Sp_buf = NULL;
		wsAlloc(Sp_buf, char, 2048);

		for (; C_sampNA.gets(Sp_buf, 2048); L++) {
			char	*a, *b;

			getString(&Sp_buf, &a);
			if (!a)	halt_fmt(WISARD_INVL_SAMPLIST_FORMAT, S_fn, L);
			getString(&a, &b);
			/* 131225 Ignores after second column */
//			if (b)	halt_fmt(WISARD_INVL_SAMPLIST_FORMAT, OPT_STRING(selsamp), L);

			string S_key = Sp_buf;
			S_key += ",";
			S_key += a;
			Xm_set.insert(S_key);
		}
		DEALLOC(Sp_buf);
	}

	return N_nonRegItem;
}

wsUint loadFamilyList(wsStrCst S_fn, wsStrCst S_desc, eStr& Xm_set)
{
	wsUint		N_nonRegItem = 0;
	cStrFile	C_sampNA(S_fn, S_desc, 1);

	if (C_sampNA.isFailed()) {
		LOG("Failed to open [%s], maybe the IID selection option is given as their FID/IID pairs?\n",
			S_fn);
		char *S_remsamp = NULL;
		wsAlloc(S_remsamp, char, strlen(S_fn) + 1);
		strcpy(S_remsamp, S_fn);

		char *S_tok = strtok(S_remsamp, ",");
		while (S_tok) {
			string S_key = S_tok;
			Xm_set.insert(S_key);
			S_tok = strtok(NULL, ",");
		}
		Xm_set.insert("__LIST__");
		N_nonRegItem++;

		DEALLOC(S_remsamp);
	}
	else {
		int		L = 1;
		char	*Sp_buf = NULL;
		wsAlloc(Sp_buf, char, 2048);

		for (; C_sampNA.gets(Sp_buf, 2048); L++) {
			char	*a;

			getString(&Sp_buf, &a);
			if (a && a[0])	halt_fmt(WISARD_INVL_SAMPLIST_FORMAT, S_fn, L);

			string S_key = Sp_buf;
			Xm_set.insert(S_key);
		}
		DEALLOC(Sp_buf);
	}

	return N_nonRegItem;
}

void loadRegexList(wsStrCst S_fn, wsStrCst S_desc, vPattern& Xv_set)
{
	wsUint		N_nonRegItem = 0;
	cStrFile	C_sampNA(S_fn, S_desc, 1);
	const TRexChar*	Sp_err = NULL;

	if (C_sampNA.isFailed()) {
		LOG("Failed to open [%s], maybe it is a single regular expression?\n",
			S_fn);
		TRex* X_regex = trex_compile(S_fn, &Sp_err);
		if (!X_regex) {
			halt("Parameter [%s] is error by [%s]!", S_fn,
				Sp_err ? Sp_err : "???");
			return;
		}
		Xv_set.push_back(X_regex);
	} else {
		int		L = 1;
		char	*Sp_buf = NULL;
		wsAlloc(Sp_buf, char, 2048);

		for (; C_sampNA.gets(Sp_buf, 2048); L++) {
			char	*a, *b;

			getString(&Sp_buf, &a);
			if (a && a[0])	halt_fmt(WISARD_INVL_SAMPLIST_FORMAT, S_fn, L);

			TRex* X_regex = trex_compile(Sp_buf, &Sp_err);
			if (!X_regex) {
				halt("Parameter [%s] on line [%d] is error by [%s]!", S_fn, L,
					Sp_err ? Sp_err : "???");
				return;
			}
			Xv_set.push_back(X_regex);
		}
		DEALLOC(Sp_buf);
	}
}

void cIO::_loadSampleList(char B_select/*=0*/)
{
	wsStrCst	S_msgRemSel		= B_select ? "selection" : "removal";
	if (OPT_ENABLED(regex)) loadRegexList(
			B_select ? OPT_STRING(selsamp) : OPT_STRING(remsamp),
			B_select ? "Sample selection list file" : "Sample removal list file",
			Xv_listSample);
	else {
		wsUint		N_nonRegItem	= loadSampleList(
			B_select?OPT_STRING(selsamp) : OPT_STRING(remsamp),
			B_select ? "Sample selection list file" : "Sample removal list file",
			Xe_listSample);

		/* Logging and clear */
		LOG("[%d] samples were scanned for %s\n",
			(int)(Xe_listSample.size() - N_nonRegItem), S_msgRemSel);
	}
	N_filterSample = char(B_select ? FILT_SELECT : FILT_REMOVE);
}

void cIO::_loadSampleNAList()
{
	wsUint	N_nonRegItem	= loadSampleList(OPT_STRING(nasamp),
		"Sample nullifying list file", Xe_listNAsamp);

	/* Logging and clear */
	LOG("[%d] samples were indicated to nullifying its genotype\n",
		(int)(Xe_listNAsamp.size() - N_nonRegItem));
}

void cIO::_loadFamilyList(char B_select/*=0*/)
{
	wsStrCst	S_msgRemSel	= B_select ? "selection" : "removal";

	if (OPT_ENABLED(regex)) loadRegexList(
		B_select ? OPT_STRING(selfam) : OPT_STRING(remfam),
		B_select ? "FID selection list file" : "FID removal list file",
		Xv_listFam);
	else {
		wsUint	N_nonRegItem	= loadFamilyList(
			B_select ? OPT_STRING(selfam) : OPT_STRING(remfam),
			B_select ? "FID selection list file" : "FID removal list file",
			Xe_listFam);

		/* Logging and clear */
		LOG("[%d] %s of families were found\n",
			(int)(Xe_listFam.size() - N_nonRegItem), S_msgRemSel);
	}
	N_filterSample = char(B_select ? FILT_SELECT : FILT_REMOVE);
}

void cIO::_loadRangeRemovalList()
{
	int			L=1;
	cStrFile	C_rngRem(OPT_STRING(filrange), "Variant removal range list file", 1);
	if (C_rngRem.isFailed()) {
		char	S_chr[32];
		int		s, e, pos;

		/* Init */
		wsAlloc(Xv_listRanges, vector<xRange>, MAX_NCHR+1);

		/* Maybe listed form? */
		char *p = OPT_STRING(filrange);
		while (1) {
			int n = sscanf(p, "%[^[][%d,%d]%n", S_chr, &s, &e, &pos);
			if (*p == '\0') break;

			/* Value check */
			if (n != 3 || s < 0 || e < 0 || s > e)
				halt("Invalid value [%s] for --filrange", OPT_STRING(filrange));
			wsUint N_chr = getChr(S_chr);
			xRange X;
			X.N_start = s;
			X.N_end = e;
			Xv_listRanges[N_chr].push_back(X);

			L++;
			p += pos;
			if (*p != ',' && *p != '\0')
                                halt("Unexpected character [%c](%d) found in [%d]th character in --filrange", *p, *p, p-OPT_STRING(filrange));
                        if (*p == ',') p++;
			if (*p == '\0') break;
		}
		L--;
	} else {
		char		*Sp_buf = NULL, *a, *b, *c;
		wsAlloc(Sp_buf, char, 2048);

		/* Init */
		wsAlloc(Xv_listRanges, vector<xRange>, MAX_NCHR+1);

		for ( ; C_rngRem.gets(Sp_buf, 2048) ; L++) {
			xRange X_rng;

			getString(&Sp_buf, &a);
			wsUint N_chr = getChr(Sp_buf);
			if (!a) halt_fmt(WISARD_INVL_RANGELIST_FORMAT, OPT_STRING(filrange), L);

			getString(&a, &b);
			X_rng.N_start = atoi(a);
			if (!b) halt_fmt(WISARD_INVL_RANGELIST_FORMAT, OPT_STRING(filrange), L);

			getString(&b, &c);
			X_rng.N_end = atoi(b);
			if (c) halt_fmt(WISARD_INVL_RANGELIST_FORMAT, OPT_STRING(filrange), L);

			/* Make data and insert */
			Xv_listRanges[N_chr].push_back(X_rng);
		}
		DEALLOC(Sp_buf);
	}

	/* Logging and clear */
	LOG("%d ranges retrieved to remove\n", L);
	N_filterRange = FILT_REMOVE;
}

void cIO::_loadRangeSelectionList()
{
	cStrFile	C_rngSel(OPT_STRING(incrange), "Variant selection range list file", 1);
	int			L=1;
	if (C_rngSel.isFailed()) {
		char	S_chr[32];
		int		s, e, pos;

		/* Init */
		wsAlloc(Xv_listRanges, vector<xRange>, MAX_NCHR+1);

		/* Maybe listed form? */
		char *p = OPT_STRING(incrange);
		while (1) {
			int n = sscanf(p, "%[^[][%d,%d]%n", S_chr, &s, &e, &pos);
			if (*p == '\0') break;

			/* Value check */
			if (n != 3 || s < 0 || e < 0 || s > e)
				halt("Invalid value [%s] for --incrange, stopped from [%s], # elements [%d]", OPT_STRING(incrange), p, n);
			wsUint N_chr = getChr(S_chr);
			xRange X;
			X.N_start = s;
			X.N_end = e;
			Xv_listRanges[N_chr].push_back(X);

			L++;
			p += pos;
			if (*p != ',' && *p != '\0')
				halt("Unexpected character [%c](%d) found in [%d]th character in --incrange", *p, *p, p-OPT_STRING(incrange));
			if (*p == ',') p++;
			if (*p == '\0') break;
		}
		L--;
	} else {
		char		*Sp_buf = NULL, *a, *b, *c;
		wsAlloc(Sp_buf, char, 2048);

		/* Init */
		wsAlloc(Xv_listRanges, vector<xRange>, MAX_NCHR+1);

		for ( ; C_rngSel.gets(Sp_buf, 2048) ; L++) {
			xRange X_rng;

			getString(&Sp_buf, &a);
			wsUint N_chr = getChr(Sp_buf);
			if (!a) halt_fmt(WISARD_INVL_RANGELIST_FORMAT, OPT_STRING(incrange), L);

			getString(&a, &b);
			X_rng.N_start = atoi(a);
			if (!b) halt_fmt(WISARD_INVL_RANGELIST_FORMAT, OPT_STRING(incrange), L);

			getString(&b, &c);
			X_rng.N_end = atoi(b);
			if (c) halt_fmt(WISARD_INVL_RANGELIST_FORMAT, OPT_STRING(incrange), L);

			/* Make data and insert */
			Xv_listRanges[N_chr].push_back(X_rng);
		}
		DEALLOC(Sp_buf);
	}

	/* Logging and clear */
	LOG("%d ranges retrived to select\n", L);
	N_filterRange = FILT_SELECT;
}

void cIO::_loadVariantSelectionList()
{
	fetchMarkerList(OPT_STRING(selvariant), Xe_listVrt, "Variant inclusion list file");
	LOG("[%d] variants indicated to select\n", (int)(Xe_listVrt.size()));
	N_filterVrt = FILT_SELECT;
}

typedef enum _xMarkerAltMode {
	AM_CHR,
	AM_NAME,
	AM_POS,
	AM_GDIST,
} xMarkerAltMode;

void cIO::_loadMarkerAltInfo(wsStrCst S_fn, int N_mode)
{
	vVariant&	Xv_vrt	= getVariant();
	vVariant*	Xv_ret	= new vVariant;
	memset(Xv_ret, 0x00, sizeof(vVariant));
	char		*Sp_buf = NULL, *a = NULL, *b = NULL;
	wsAlloc(Sp_buf, char, 2048);

	/* Build marker map */
	mDataIdx	Xm_mmap[MAX_NCHR+1];
	wsUint		j = 0;
	FOREACHDO (vVariant_it, Xv_vrt, i, j++) {
		if (i->chr > 0 && (wsUint)i->chr <= NCHR_SPECIES)
			Xm_mmap[i->chr].insert(make_pair((string)i->name, j));
		else
			Xm_mmap[0].insert(make_pair((string)i->name, j));
	}


	cStrFile	C_md(S_fn, "Variant alternative information file");
	wsUint	N_alt = 0;
	for (wsUint L=1 ; C_md.gets(Sp_buf, 2048) ; L++) {
		int		N_insIdx = -1;
		/* name chr */
		getString(&Sp_buf, &a);
		/* find */
		for (wsUint N_chr=0 ; N_chr<=NCHR_SPECIES ; N_chr++) {
			mDataIdx_it X_find = Xm_mmap[N_chr].find(Sp_buf);
			/* If find, record index */
			if (X_find != Xm_mmap[N_chr].end()) {
				N_insIdx = X_find->second;
				break;
			}
		}
		/* Ignore this line if no match found */
		if (N_insIdx == -1 || a == NULL)
			continue;
		N_alt++;
		getString(&a, &b);
		/* Get next */
		switch (N_mode) {
		case AM_CHR: {
			int cur_chr = getChr(a);
			Xv_vrt[N_insIdx].chr = cur_chr;
		} break;
		case AM_GDIST: {
			wsReal cur_chr = (wsReal)atof(a);
			Xv_vrt[N_insIdx].gdist = cur_chr;
		} break;
		case AM_NAME: {
			free(Xv_vrt[N_insIdx].name);
			Xv_vrt[N_insIdx].name = strdup(a);
		} break;
		case AM_POS: {
			Xv_vrt[N_insIdx].pos = atoi(a);
		} break;
		default:
			halt("Invalid alteration mode [%d]", N_mode);
		}
	}
	LOG("[%d] variants altered by [%s]\n", N_alt, S_fn);
}

vVariant* cIO::_loadVariantDef(wsStrCst S_fn)
{
	vVariant *Xv_ret = new vVariant;
	memset(Xv_ret, 0x00, sizeof(vVariant));
	char		*Sp_buf = NULL, *a = NULL, *b = NULL;
	wsAlloc(Sp_buf, char, 2048);

	cStrFile	C_md(S_fn, "Variants information file");
	for (wsUint L=1 ; C_md.gets(Sp_buf, 2048) ; L++) {
		xVariant X_in;
		char *Sp_null = NULL;

		/* Chr */
		getString(&Sp_buf, &a);
		X_in.chr = getChr(Sp_buf);

		/* name */
		if (a == NULL) halt("Unexpected end in line [%d]", L);
		getString(&a, &b);
		X_in.name = strdup(a);

		/* gdist */
		if (b == NULL) halt("Unexpected end in line [%d]", L);
		getString(&b, &a);
		X_in.gdist = (wsReal)str2dbl(b, &Sp_null);
		if (Sp_null) halt("Invalid genetic distance format [%s] in line [%d]", b, L);

		/* pos */
		if (a == NULL) halt("Unexpected end in line [%d]", L);
		getString(&a, &b);
		X_in.pos = (wsUint)strtol(a, &Sp_null, 10);
		if (Sp_null) halt("Invalid physical position format [%s] in line [%d]", a, L);

		Xv_ret->push_back(X_in);
	}

	/* Return */
	return Xv_ret;
}

void cIO::_loadVariantRemovalList()
{
	fetchMarkerList(OPT_STRING(remvariant), Xe_listVrt, "Variants filtering list file");
// 	cStrFile	C_snvRem(OPT_STRING(remmarker), "Variant filtering list file");
// 	int			L=1;
// 	char		*Sp_buf = NULL, *a;
// 
// 	MULTI_MALLOC(Sp_buf, char, 2048);
// 	for ( ; C_snvRem.gets(Sp_buf, 2048) ; L++) {
// 		getString(&Sp_buf, &a);
// 		if (a) halt_fmt(WISARD_INVL_SNPLIST_FORMAT, OPT_STRING(remmarker), L);
// 		Xm_listSNP.insert(make_pair(Sp_buf, 1));
// 	}
// 	DEALLOC(Sp_buf);

	// ?I ??I??????????????A H_sampRem ??A?I ?????IAI??I??IAI
	// FID and IID??| ??a??i CuAA??I ?????????????? ??O?????I???? CI??c,
	// cIO???? ?u?yCN ??a??o ?????o??| ????A??CI??i ??y?????U/??O??e?U???? ?�J?u CO??c?I ??e??i??? ?????i ?I???? ??uCN Ca?�J??? ??i??AC?????? CO??I??U.
	LOG("[%d] variants indicated to remove\n", (int)(Xe_listVrt.size()));
	N_filterVrt = FILT_REMOVE;
}

void cIO::_loadVariantVar(vStr &Xv_reqVars)
{
	if (!IS_ASSIGNED(variantvar)) return;

	char*		Sp_mv	= strdup(OPT_STRING(variantvar));
	wsUint		N_mkrVarFile	= 0;
	char**		Sa_mvs	= loadStringValues(Sp_mv, &N_mkrVarFile);
	cStrFile*	Xa_file	= NULL;

	vStr		Xv_header;		///< [#col] Column names, sequantial by file
	vBool		Xv_read;		///< [#col] Whether read this column or not

	vInt		Xv_ncol;		///< [#file] Number of columns in the file
	vVariant&	Xv_vrt	= getVariant();
	vector<char*>
				Xv_ptrs;
	mDataIdx	Xm_header;		///< [Column name] to [existance]
	mDataIdx	Xm_idxStore;	///< [Required column name] to [stored index]

	wsUint		N_reqCols	= (wsUint)Xv_reqVars.size();
	vBool		Xv_reqVarFound;	///< [#req] Is the requested variable found?
	Xv_reqVarFound.resize(Xv_reqVars.size());

	/*
	 *
	 * Construct marker map
	 * 
	 */
	mDataIdx	Xm_mmap[MAX_NCHR + 1];
	{
		wsUint		j = 0;
		FOREACHDO (vVariant_it, Xv_vrt, i, j++) {
			if (i->chr > 0 && (wsUint)i->chr <= NCHR_SPECIES)
				Xm_mmap[i->chr].insert(make_pair((string)i->name, j));
			else
				Xm_mmap[0].insert(make_pair((string)i->name, j));
		}
	}

	/*
	 *
	 * File check
	 * 
	 */
	wsAlloc(Xa_file, cStrFile, N_mkrVarFile);
	Xv_ncol.resize(N_mkrVarFile);
	char *Sp_buf = NULL;
	wsAlloc(Sp_buf, char, 65536);
	for (wsUint i=0,j=0 ; i<N_mkrVarFile ; i++) {
		Xa_file[i].open(Sa_mvs[i], "Variant variable file");

		/* Try to read one line */
		char *Sp_tmp = Xa_file[i].gets(Sp_buf, 65536);
		if (Sp_tmp == NULL)
			halt("Given variant variable file [%s] is invalid!", Sa_mvs[i]);

		/* Retrieve header from given file */
		char *a = Sp_buf, *b = NULL;
		for ( ; a ; a=b) {
			/* Extract element, and check duplication */
			getString(&a, &b);
			string		S_cur	= (string)a;
			mDataIdx_it	X_find	= Xm_header.find(S_cur);
			if (X_find != Xm_header.end())
				halt("Header [%s] in variant variable file [%s] is duplicated!",
					a, Sa_mvs[i]);

			/* Mark if this column is required */
			int N_reqVarIdx = findVstr(Xv_reqVars, S_cur);
			if (N_reqVarIdx != -1) {
				Xv_reqVarFound[N_reqVarIdx] = 1;
				Xm_idxStore[S_cur] = N_reqVarIdx;
			}

			/* It is sane, insert it to map */
			Xm_header[S_cur] = (char)(j++);
			Xv_header.push_back(S_cur);
			Xv_ncol[i]++;
		}
	}

	/* Halt if there is a required column that cannot found */
	wsUint I = 0;
	FOREACHDO (vBool_it, Xv_reqVarFound, i, I++)
		if (*i == 0) halt("Required variant variable [%s] is not found!",
			Xv_reqVars[I].c_str());

	/* Allocate buffer of retrieved marker variables */
	char***	Sa_vals		= NULL;
	wsAlloc(Sa_vals, char**, N_reqCols);
	for (wsUint i=0 ; i<N_reqCols ; i++)
		wsAlloc(Sa_vals[i], char*, Xv_vrt.size());

	/* Now, read'em all */
	for (wsUint i=0,j=0 ; i<N_mkrVarFile ; i++) {
		/* Try to read one line */
		for (wsUint L=1 ; Xa_file[i].gets(Sp_buf, 1024*1024) ; L++) {
			char *a = Sp_buf, *b = NULL;

			/* First column = MARKER NAME */
			getString(&a, &b);
			int	N_idx = -1;
			string	S_marker = (string)a;
			for (wsUint k=0 ; k<=NCHR_SPECIES ; k++) {
				mDataIdx_it X_find = Xm_mmap[k].find(S_marker);
				if (X_find != Xm_mmap[k].end()) {
					N_idx = X_find->second;
					break;
				}
			}
			/* Skip if not exists */
			if (N_idx == -1) continue;

			/* Store entire buffer */
			Xv_ptrs.push_back(Sp_buf);
			a = b;

			/* Now, get the values */
			for (wsUint Q=j ; a ; getString(&a, &b),Q++) {
				/* IF it is REQUIRED */
				Sa_vals[Q][N_idx] = a;
			}
			wsAlloc(Sp_buf, char, 1024*1024);
		}
		/* Increase index */
		j += Xv_ncol[i];
	}

	DEALLOC(Sp_buf);
}

char findEqualSign(wsStrCst Sp_expr, char *Sp_ineqSign)
{
	char B_1eq = 0;

	/* Determine equal sign */
	if (*(Sp_ineqSign-1) == '=') {
		*(Sp_ineqSign-1) = '\0';
		B_1eq = 1;
	}
	if (*(Sp_ineqSign+1) == '=') {
		/* Check double sign */
		if (B_1eq) {
			*(Sp_ineqSign-1) = '=';
			halt("Double equal sign found at [%s]", Sp_expr);
		}
		B_1eq = 1;
		*(Sp_ineqSign+1) = '\0';
	}

	return B_1eq;
}

void process(xVarRange *Xp_filt, wsStrCst Sp_tok, char *Sp_1, char *Sp_2,
	char N_sign)
{
	char B_1eq = 0, B_2eq = 0;
	wsReal R_1, R_2 = WISARD_NAN;

	/* tok (1) some (2) else */

	/* Should satisfy tok != Sp_1 && *(Sp_2+1) != '\0' */
	if (Sp_tok==Sp_1 || (Sp_2 && *(Sp_2+1)=='\0'))
		halt("Invalid format [%s]", Sp_tok);

	/* Determine inequality sign */
	B_1eq = findEqualSign(Sp_tok, Sp_1);
	if (Sp_2)
		B_2eq = findEqualSign(Sp_tok, Sp_2);

	/* (tok) somevalue (1) somepheno (2) somevalue */
	*(Sp_1++) = '\0';
	if (*Sp_1 == '\0') Sp_1++;
	if (Sp_2) {
		*(Sp_2++) = '\0';
		if (*Sp_2 == '\0') Sp_2++;
	}

	char *Sp_end = NULL;
	char B_phenoPos = -1;
	if (Sp_2) {
		/* value < pheno < value type */
		R_1 = (wsReal)str2dbl(Sp_tok, &Sp_end);
		if (*Sp_end != '\0') {
			/* Is it phenotype name? */
			/* Otherwise, can't proceed = maybe error */
			halt("Invalid format, first element of %s should be real number!", Sp_1);
		}

		R_2 = (wsReal)str2dbl(Sp_2, &Sp_end);
		if (*Sp_end != '\0') halt("Invalid format, second element of %s"
			" should be real number!", Sp_1);
		/* Assign phenotype name */
		strcpy(Xp_filt->S_varName, Sp_1);

		/* Sign check */
		if (R_1==R_2 && (!B_1eq || !B_2eq))
			halt("Invalid range, cannot include anything!");
		if ((N_sign=='<' && R_1>R_2) || (N_sign=='>' && R_1<R_2))
			halt("Invalid sign of %s found", Sp_1);
	} else {
		/* try value < pheno */
		R_1 = (wsReal)str2dbl(Sp_tok, &Sp_end);
		if (*Sp_end != '\0') {
			/* try pheno < value */
			R_1 = (wsReal)str2dbl(Sp_1, &Sp_end);
			if (*Sp_end != '\0') halt("Cannot found any value");
			strcpy(Xp_filt->S_varName, Sp_tok);
			B_phenoPos = 0;
		} else {
			strcpy(Xp_filt->S_varName, Sp_1);
			B_phenoPos = 1;
		}
	}

	/* Assign value */
	if (Sp_2) {
		if (N_sign=='<') {
			Xp_filt->R_s = R_1;
			Xp_filt->R_sEQ = B_1eq;
			Xp_filt->R_e = R_2;
			Xp_filt->R_eEQ = B_2eq;
		} else {
			Xp_filt->R_s = R_2;
			Xp_filt->R_sEQ = B_2eq;
			Xp_filt->R_e = R_1;
			Xp_filt->R_eEQ = B_1eq;
		}
	} else {
		/*
		A<value
		value>A

		A>value
		value<A
		*/

		if ((N_sign=='<'&&B_phenoPos==0) ||
				(N_sign=='>'&&B_phenoPos==1)) {
			Xp_filt->R_s	= numeric_limits<wsReal>::infinity();
			Xp_filt->R_sEQ	= 0;
			Xp_filt->R_e	= R_1;
			Xp_filt->R_eEQ	= B_1eq;
		} else if ((N_sign=='<'&&B_phenoPos==1) ||
				(N_sign=='>'&&B_phenoPos==0)) {
			Xp_filt->R_e	= numeric_limits<wsReal>::infinity();
			Xp_filt->R_eEQ	= 0;
			Xp_filt->R_s	= R_1;
			Xp_filt->R_sEQ	= B_1eq;
		} else
			halt("SYSERR: Should not reach to here");
	}
}

char isInVarRange(xVarRange& X_range, wsMat Ra_data, wsUint N_srcIdx,
	wsUint N_offset)
{
	wsReal S = numeric_limits<wsReal>::infinity();
	wsReal E = numeric_limits<wsReal>::infinity();
	wsReal V = Ra_data[N_srcIdx][N_offset];

	/* Start is INDEX */
	if (X_range.N_sIdx != -1) {
		S = Ra_data[X_range.N_sIdx][N_offset];
		if (isMissingReal(S)) return 1; /* Do NOT filter if missing */
	/* Have specific value */
	} else if (X_range.R_s == X_range.R_e && X_range.N_eIdx == -1)
		return X_range.R_s == V;
	else
		/* Start is VALUE */
		S = X_range.R_s;

	/* End is INDEX/VALUE */
	E = X_range.N_eIdx != -1 ? Ra_data[X_range.N_eIdx][N_offset] :
		X_range.R_e;
	if (isMissingReal(E)) return 1; /* Do NOT filter if missing */

	/* Greater than R_s */
	if (X_range.R_e == numeric_limits<wsReal>::infinity()) {
		return X_range.R_sEQ ? S<=V : S<V;
	} else if (X_range.R_s == numeric_limits<wsReal>::infinity()) {
		return X_range.R_eEQ ? V<=E : V<E;
	} else {
		return X_range.R_sEQ ?
			(X_range.R_eEQ ? S<=V && V<=E : S<=V && V<E)
			: (X_range.R_eEQ ? S<V && V<=E : S<V && V<E);
	}
}

inline char _isThisFiltered(vVarRng &X_pr, wsMat Ra_data, wsUint N_srcIdx,
	wsUint N_offset)
{
	/* For all filter given to this phenotype, check range */
	FOREACH (vVarRng_it, X_pr, i)
		if (!isInVarRange(*i, Ra_data, N_srcIdx, N_offset) ^ i->B_filt) return 1;
	return 0;
}

void cIO::_filtVar(wsStrCst Sp_filExpr, char B_filt, char B_forCov)
{
	vVarRng	*Xp_filts = B_forCov ? &V_filtCovs : &V_filtPhenos;
	if (Sp_filExpr == NULL) return;

	if (!Xa_covInfo.size() && B_forCov)
		halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--filcov/--inccov",
			"no covaritate was assigned");

	/* Copy expressions in order to divide */
	char *Sp_expr = NULL;
	wsAlloc(Sp_expr, char, strlen(Sp_filExpr)+1);
	strcpy(Sp_expr, Sp_filExpr);

	/* Divide expressions into each expression */
	char *Sp_tok = strtok(Sp_expr, ",");
	do {
		xVarRange	X_filt		= { 0, };
		char		B_subFilt	= B_filt;
		
		/* Step 1 : Find < or > */
		char *Sp_2 = NULL;
		char *Sp_1 = strchr(Sp_tok, '<');
		if (Sp_1 == NULL) {
			/* Find > */
			Sp_1 = strchr(Sp_tok, '>');
			if (Sp_1 == NULL) {
				/* Find = */
				Sp_1 = strchr(Sp_tok, '=');
				if (Sp_1 == NULL)
					halt("Invalid phenotype filtering format [%s]", Sp_tok);
				/* Find ! */
				if (*(Sp_1-1) == '!')
					B_subFilt = !B_subFilt;
				*(Sp_1++) = '\0';

				/* Is it 'real number'? */
				char *Sp_end = NULL;
				wsReal R_v = (wsReal)str2dbl(Sp_1, &Sp_end);

				if (Sp_end && Sp_end[0]) {
					/* COVNAME=NONREAL type */
					if (!B_forCov)
						/* Only for covariates */
						halt("[%s] including factor-type, and it can be only applied to the covariates", Sp_filExpr);
					halt("Currently, filters for factor covariates are not supported");
// 					X_filt.Sp_varValue = strdup(Sp_1);
// 					X_filt.R_s		= X_filt.R_e = numeric_limits<wsReal>::infinity();
				} else
					/* PHENONAME=VAL type */
					X_filt.R_s		= X_filt.R_e = R_v;
				strcpy(X_filt.S_varName, Sp_tok);
				X_filt.N_sIdx	= X_filt.N_eIdx = -1;
				X_filt.R_eEQ	= X_filt.R_sEQ = 1;
				X_filt.B_filt	= B_subFilt;
				Xp_filts->push_back(X_filt);
			}
		} else {
			/* Here we have Sp_1, find Sp_2 */
			char S_1 = *Sp_1;
			char S_dame = S_1=='>' ? '<' : '>';

			Sp_2 = strchr(Sp_tok, S_dame);
			if (Sp_2 != NULL)
				/* tok_f < tok_dame > else || tok_f > tok_dame < else */
				halt("Invalid format [%s]", Sp_tok);

			/* tok > else || tok < else VALID
			* Find out else have more sign */
			Sp_2 = strchr(Sp_1+1, S_dame);
			if (Sp_2 != NULL)
				/* tok < else_f > else_dame || tok > else_f < else_dame */
				halt("Invalid format [%s]", Sp_tok);
			Sp_2 = strchr(Sp_1+1, S_1);

			/* tok (1) else */
			process(&X_filt, Sp_tok, Sp_1, Sp_2, S_1);
			X_filt.B_filt = B_subFilt;
		}

		/* Check this phenotype is included */
		if (B_forCov) {
			if (Xa_covInfo.size()) {
				vCovar_it it = Xa_covInfo.begin();
				for ( ; it!=Xa_covInfo.end() ; it++) {
					if (!stricmp(it->Sp_varName, X_filt.S_varName)) {
						it->Xa_filters.push_back(X_filt);
						break;
					}
				}
				if (it == Xa_covInfo.end())
					halt("Covariate [%s] was not found on the available "
						"covariates list", X_filt.S_varName);
			} else
				/* --filcov cannot used if there is no covariates */
				halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--filcov",
					"no covariates assignment");
		} else {
			if (V_phenos.size()) {
				vPheno_it it = V_phenos.begin();
				for ( ; it!=V_phenos.end() ; it++) {
					wsStrCst Sp = it->S_name.c_str();
					if (!stricmp(Sp, X_filt.S_varName)) {
						it->Xa_filters.push_back(X_filt);
						break;
					}
				}
				if (it == V_phenos.end())
					halt("Phenotype [%s] was not found on the available "
						"phenotype list", X_filt.S_varName);
			} else if (stricmp(X_filt.S_varName, "PHENOTYPE"))
				/* If there is no alternative phenotype, the name of phenotype
				* must be PHENOTYPE */
				halt("Default phenotype name must be `PHENOTYPE`");
			else
				/* Insert to DEFAULT phenotype filter */
				V_filtPhenos.push_back(X_filt);
		}
	} while ((Sp_tok = strtok(NULL, ",")));

	/* Print criteria */
	if (B_forCov) FOREACH (vVarRng_it, V_filtCovs, i) {
		LOG("Covariate [%s] filtered with ", i->S_varName);
		if (i->N_sIdx != -1)
			LOGnf("[%s]", Xa_covInfo[i->N_sIdx].Sp_varName);
		else
			LOGnf("[%g]", i->R_s);
		LOGnf("%s VAL %s", i->R_sEQ?"<=":"<", i->R_eEQ?"<=":"<");
		if (i->N_eIdx != -1)
			LOGnf("[%s]\n", Xa_covInfo[i->N_eIdx].Sp_varName);
		else
			LOGnf("[%g]\n", i->R_e);
	} else FOREACH (vVarRng_it, V_filtPhenos, i) {
		LOG("Phenotype [%s] filtered with ", i->S_varName);
		if (i->N_sIdx != -1)
			LOGnf("[%s]", V_phenos[i->N_sIdx].S_name.c_str());
		else
			LOGnf("[%g]", i->R_s);
		LOGnf("%s VAL %s", i->R_sEQ?"<=":"<", i->R_eEQ?"<=":"<");
		if (i->N_eIdx != -1)
			LOGnf("[%s]\n", V_phenos[i->N_eIdx].S_name.c_str());
		else
			LOGnf("[%g]\n", i->R_e);
	}

	/* Do filtering */
	wsUint	I=0;
	wsMat	Ra_allPheno	= getPhenos();
	wsMat	Ra_covs		= getCovariates();
	if (B_forCov) FOREACHDO (vCovar_it, Xa_covInfo, i, I++) {
		if (i->X_type == WISARD_VAR_FACTOR) {
			/* FIXME : Update needed */
// 			for (wsUint j=0 ; j<N_sample ; j++)
// 				FOREACH (vVarRng_it, i->Xa_filters, ii)
// 					if (ii->Sp_varValue) {
// 						if (ii->B_filt && !stricmp(ii->Sp_varValue, i->Sp_bl))
// 					}
		} else for (wsUint j=0 ; j<N_sample ; j++)
			/* Perform filtering */
			if (!Bv_filtSamp[j] && _isThisFiltered(i->Xa_filters, Ra_covs, I, j)) {
				N_szFilteredSamp++;
				Bv_filtSamp[j] = 1;
			}
	} else FOREACHDO (vPheno_it, V_phenos, i, I++) {
		for (wsUint j=0 ; j<N_sample ; j++)
			/* Perform filtering */
			if (!Bv_filtSamp[j] && _isThisFiltered(i->Xa_filters, Ra_allPheno, I, j)) {
				N_szFilteredSamp++;
				Bv_filtSamp[j] = 1;
			}
	}

	/* Cleanup */
	DEALLOC(Sp_expr);
}

void cIO::_loadBlup()
{
	char		*a, *b;
	cStrFile	C_blup(OPT_STRING(blup), "Phenotype BLUP file");
	char		*Sp_buf	= NULL;
	wsUint		N_pheno	= sizePheno();
	vStr		Xa_phenos, Xa_covs;
	strMap	Xa_phenoExists;
	wsUint	N_blupField = 0;
	char	*Ba_isIncPheno = NULL, *Ba_isIncCov = NULL;
	mSamp&	Xm_sample = getSampleData();
	wsUint	N_blupSample = 0;

// 	if (N_pheno != 1)
// 		halt("Can't apply --blup, because it is multi-phenotype (%d phenos)", N_pheno);
	wsAlloc(Sp_buf, char, 4096);

	/* First column should be header */
	if (!C_blup.gets(Sp_buf, 4096))
		halt("Can't read blup file header!");

	/* Format check ; first two columns should be 'FID' and 'IID' */
	getString(&Sp_buf, &b);
	if (b == NULL)
		halt("Given blup file have less two columns in the first row");
	getString(&b, &a);
	if (a == NULL)
		halt("Given blup file have no field");
	if (strcmp(Sp_buf, "FID") || strcmp(b, "IID"))
		halt("First two columns of BLUP file should be 'FID' and 'IID'");

	vStr	X_blupFields;
	vInt	X_blupIdx;
	for ( ; ; a=b,N_blupField++) {
		getString(&a, &b);
		if ((N_blupField%2) == 0) {
			X_blupFields.push_back(a);
			X_blupIdx.push_back(-1);
		}
		if (b == NULL) break;
	}
	N_blupField++;
	if (N_blupField%2)
		halt("# of BLUP fields[%d] cannot be odd", N_blupField);
	N_blupField >>= 1;
	LOG("%d fields were found in blup file\n", N_blupField);
	
	/* # of fields must match with the # of phenotypes */
	if (N_blupField != N_pheno)
		halt("# of fields[%d] does not match with # of phenotypes[%d]",
			N_blupField, N_pheno);
	/* Which idx to be assigned? */
	vPheno	&X_phenos = getPhenoInfo();
	wsUint I=0;
	FOREACHDO (vPheno_it, X_phenos, i, I++) {
		string S = "BLUP_";
		S += i->S_name;
		wsUint J=0;
		FOREACHDO (vStr_it, X_blupFields, j, J++)
			if (!j->compare(S)) {
				X_blupIdx[J] = I;
				break;
			}
	}
	/* Sanity check */
	for (wsUint i=0 ; i<N_pheno ; i++)
		if (X_blupIdx[i] == -1)
			halt("BLUP column [%s] failed to matching, is there proper phenotype?",
				X_blupFields[i].c_str());

	/* N_blupSample the number */
	for (N_blupSample=0 ; C_blup.gets(Sp_buf, 4096) ; N_blupSample++);
	LOG("%d samples were found at BLUP file\n", N_blupSample);

	/* Retrieve */
	char *Ba_blupApplied = NULL;
	wsCalloc(Ba_blupApplied, char, sizeof(char)*N_sample);
	C_blup.rewind();
	if (!C_blup.gets(Sp_buf, 4096))
		halt("Can't read phenotype header!");
	for (N_blupSample=0 ; C_blup.gets(Sp_buf, 4096) ; N_blupSample++) {
		/* Get FID and IID */
		getString(&Sp_buf, &b);
		if (b == NULL) halt("%d th line have insufficient data", N_blupSample+1);
		getString(&b, &a);
		if (a == NULL) halt("%d th line have no phenotype data", N_blupSample+1);

		/* Get sample data */
		mSamp_it X_resSampleFind =
			Xm_sample.find(b);
		if (X_resSampleFind == Xm_sample.end())
			halt("%s::%s is not exist in the data", Sp_buf, b);
		xSample	&X_samp = X_resSampleFind->second;
		if (X_samp.S_FID.compare(Sp_buf))
			halt("Sample [%s::%s] registered [%s] family in data",
			Sp_buf, b, X_samp.S_FID.c_str());

		/* Give the phenotype unless it is not filtered */
		if (X_samp.N_idx != SAMP_NODATA) {
//			MYuint j=0, Jp=0, Jc=0;
			wsUint K=0, j=0;
			for ( ; j<N_blupField ; j++,K++) {
				char *q = NULL;
				getString(&a, &b);
				getString(&b, &q);

				/* Cannot be NA */
				if (isMissingReal(Ra_pheno[X_blupIdx[K]][X_samp.N_idx]))
					halt("BLUP for sample [%s] assigned, but its phenotype [%s]"
						" seems to NA", X_samp.S_IID.c_str(),
						V_phenos[X_blupIdx[K]].S_name.c_str());

				/* Set Y-blup ~ variant */
				Ra_pheno[X_blupIdx[K]][X_samp.N_idx] -= (wsReal)atof(a)
					+ (wsReal)atof(b);
				/* And set it to continuous */
				//Ba_typePheno[0] = 1;
				Ba_blupApplied[X_samp.N_idx] = 1;
				a = q;
			}
// 			if (b != NULL) {
// 				getString(&a, &b);
// 				if (a[0])
// 					halt("%d th line [%s::%s] have invalid number of phenotypes",
// 					N_altSample+1, X_samp.S_FID.c_str(), X_samp.S_IID.c_str());
// 			}
		}
	}
	wsUint N_blupExcl = 0;
	for (wsUint i=0 ; i<N_sample ; i++) {
		if (Ba_blupApplied[i]) continue;

		wsUint j=0;
		for (j=0 ; j<N_pheno ; j++)
			if (isMissingReal(Ra_pheno[j][i])) break;
		if (j == N_pheno) {
// 			halt("Sample [%s:%s] have nonmissing phenotype but BLUP was not found",
// 				Xa_sampleV2[i]->S_FID.c_str(),
// 				Xa_sampleV2[i]->S_IID.c_str());
			LOG("Sample [%s:%s] excluded from the subsequent analyses\n",
				Xa_sampleV2[i]->S_FID.c_str(), Xa_sampleV2[i]->S_IID.c_str());
			for (j=0 ; j<N_pheno ; j++)
				Ra_pheno[j][i] = WISARD_NA_REAL;
			N_blupExcl++;
		}
	}
	if (N_blupExcl)
		LOG("[%d] samples were additionally excluded excluded from subsequent "
			"analyses because those samples does not have BLUP\n", N_blupExcl);
//	if (N_blupSample != N_fSamp)
//		halt("# of samples in BLUP file is not match with original dataset");
	LOG("BLUP for [%d] phenotypes with [%d] samples retrieved\n", N_pheno,
		N_blupSample);

	/* Change phenotype name with '-BLUP' */
	FOREACHDO (vPheno_it, X_phenos, i, I++) {
		string S = i->S_name;
		S += "-BLUP";
		i->S_name = S;
	}
	B_blupApplied = 1;

	DEALLOC(Sp_buf);
	DEALLOC(Ba_isIncPheno);
	DEALLOC(Ba_isIncCov);
}

void cIO::_loadVariantWeight()
{
	/* Execute if only --weight have given */
	if (!IS_ASSIGNED(weight))
		return;
	wsUint		N_idxLine	= 1;
	char		*S_buf		= NULL;
	cStrFile	C_sw(OPT_STRING(weight), "Variant-wise weight file");
	wsAlloc(S_buf, char, 2048);

	/* This weight file should satisfy below conditions
	* 
	* 1. Each row should be consist of two columns
	* 2. The first column should be equal to the name of variant in the map file
	*  that describing corresponding variant
	* 3. The second column should be the value that will be the weight of
	*  corresponding variant, and the range of the weight should be placed between
	*  [0,1]
	*/

	/* Allocate space for customized weight */
	sseMalloc(Ra_customWeight, wsReal, sizeVariant());

	/* Build variant-index map */
	vVariant&	Xa_snv = getVariant();
	multimap<string, int> Xa_snvMap;
	{
		int index;
		vector<xVariant>::iterator iter;
		/* Build variant map */
		for( iter = Xa_snv.begin(), index = 0; iter != Xa_snv.end(); iter++, index++)
			Xa_snvMap.insert( pair<string, int>(string(iter->name), index) );
	}

	/* Retrieving variant-wise weight information */
	wsUint N_inserted = 0;
	for ( ; C_sw.gets(S_buf, 2048) ; N_idxLine++) {
		char *Sp_wgt = NULL;
		if ((N_idxLine%10000) == 0)
			notice("%d lines retrieved\r", N_idxLine);

		/* Get the name of variant */
		getString(&S_buf, &Sp_wgt);
		if (Sp_wgt == NULL)
			halt("Invalid line formatting in line %d", N_idxLine);

		/* Find out the index of corresponding variant */
		multimap<string,int>::iterator X_res = Xa_snvMap.find(string(S_buf));
		if (X_res == Xa_snvMap.end())
			continue;
		Ra_customWeight[X_res->second] = (wsReal)atof(Sp_wgt);
		N_inserted++;
	}
	LOG("%d lines retrieved from weight definiton file\n", N_idxLine);
	/* Number checking */
	if (N_inserted != sizeVariant())
		halt("# of matched weights in the file[%d] is not match with the "
			"number of variants in current dataset[%d]", N_inserted, sizeVariant());

	DEALLOC(S_buf);
}

inline wsUint buildVariantFilterByMAF(xOptRange& X_rng, char B_filtVal,
	vBool &Ba_filter, wsUint N_vrt, vVariant &Xv_vrt, xMaf*	Xp_maf)
{
	wsUint N_szFiltVrt = 0;
	LOOP (i, N_vrt) {
		if (Ba_filter[i]) continue;
		xMaf&		X_maf	= Xp_maf[i];

		/* Should filtered */
		if (isInRange(X_rng, X_maf.R_maf) == B_filtVal) {
			N_szFiltVrt++;
			Ba_filter[i] = MFR_MAF;
		}
	}

	return N_szFiltVrt;
}

inline wsUint buildVariantFilterByMAC(xOptRange& X_rng, char B_filtVal,
	vBool &Ba_filter, wsUint N_vrt, vVariant &Xv_vrt, xMaf*	Xp_maf)
{
	wsUint N_szFiltVrt = 0;
	for (wsUint i=0 ; i<N_vrt ; i++) {
		/* Do not double-check */
		if (Ba_filter[i]) continue;
		xMaf&		X_maf	= Xp_maf[i];

		/* Should filtered */
		if (isInRange(X_rng, (wsReal)X_maf.N_allMac) == B_filtVal) {
			N_szFiltVrt++;
			Ba_filter[i] = MFR_MAC;
		}
	}

	return N_szFiltVrt;
}

inline wsUint buildVariantFilterByGRATE(xOptRange& X_rng, char B_filtVal,
	vBool &Ba_filter, wsUint N_vrt, vVariant &Xa_vrt, wsUint N_samp)
{
	wsUint N_szFiltVrt = 0;
	for (wsUint i=0 ; i<N_vrt ; i++) {
		xVariant &X_vrt = Xa_vrt[i];
		/* Do not double-check */
		if (Ba_filter[i]) continue;
		wsReal R_grate = W1 - X_vrt.nmis / (wsReal)N_samp;

		/* Should filtered */
		if (isInRange(X_rng, R_grate) == B_filtVal) {
			N_szFiltVrt++;
			Ba_filter[i] = MFR_GTRATE;
		}
	}

	return N_szFiltVrt;
}

void cIO::_md()
{
	if (!OPT_ENABLED(mendel) &&
		!IS_ASSIGNED(filmendelfam) && !IS_ASSIGNED(incmendelfam) &&
		!IS_ASSIGNED(filmendelsamp) && !IS_ASSIGNED(incmendelsamp) &&
		!IS_ASSIGNED(filmendelvar) && !IS_ASSIGNED(incmendelvar))
		return;

	/* Do FS */
	cFamStrAnalysis C_anaFS(this);
	C_anaFS.run();

	/* Do ME */
	cMendelAnalysis C_anaMD(this, &C_anaFS);
	C_anaMD.run();
	xMendelStat& X_t = C_anaMD.getMendelStat();

	/* If famrem/inc */
	if (IS_ASSIGNED(filmendelfam) || IS_ASSIGNED(incmendelfam)) {
		char		B_filtVal	= IS_ASSIGNED(filmendelfam) ? 1 : 0;
		xOptRange	X_rng		= IS_ASSIGNED(filmendelfam) ?
			OPT_RANGE(filmendelfam) : OPT_RANGE(incmendelfam);

		for (wsUint i=0 ; i<N_sample ; i++) {
			xSample	*Xp_samp	= Xa_sampleV2[i];
			string	&S_FID		= Xp_samp->S_FID;	
			/* Do not double-check */
			if (Bv_filtSamp[i]) continue;

			/* Should filtered */
			wsReal R_ratio = (wsReal)X_t.Xm_ffam[S_FID] / (wsReal)X_t.Xm_fam[S_FID];
			if (isInRange(X_rng, R_ratio) == B_filtVal) {
				N_szFilteredSamp++;
				Bv_filtSamp[i] = 1;
			}
		}
	}
	/* If samprem/inc */
	if (IS_ASSIGNED(filmendelsamp) || IS_ASSIGNED(incmendelsamp)) {
		char		B_filtVal	= IS_ASSIGNED(filmendelsamp) ? 1 : 0;
		xOptRange	X_rng		= IS_ASSIGNED(filmendelsamp) ?
			OPT_RANGE(filmendelsamp) : OPT_RANGE(incmendelsamp);

		for (wsUint i=0 ; i<N_sample ; i++) {
			xSample *Xp_samp = Xa_sampleV2[i];
			/* Do not double-check */
			if (Bv_filtSamp[i]) continue;

			/* Should filtered */
			wsReal R_ratio = X_t.Xv_fsamp[Xp_samp->N_idx] / X_t.Xv_samp[Xp_samp->N_idx];
			if (isInRange(X_rng, R_ratio) == B_filtVal) {
				N_szFilteredSamp++;
				Bv_filtSamp[i] = 1;
			}
		}
	}
	/* If markerrem/inc */
	if (IS_ASSIGNED(filmendelvar) || IS_ASSIGNED(incmendelvar)) {
		char		B_filtVal	= IS_ASSIGNED(filmendelvar) ? 1 : 0;
		xOptRange	X_rng		= IS_ASSIGNED(filmendelvar) ?
			OPT_RANGE(filmendelvar) : OPT_RANGE(incmendelvar);

		for (wsUint i=0 ; i<N_variant ; i++) {
//			xMarker &X_var = Xa_marker[i];
			/* Do not double-check */
			if (Bv_filtVrt[i]) continue;

			/* Should filtered */
			wsReal R_ratio = X_t.Xv_fmarker[i] / X_t.Xv_marker[i];
			if (isInRange(X_rng, R_ratio) == B_filtVal) {
				N_szFilteredVrt++;
				Bv_filtVrt[i] = MFR_MENDELIAN;
			}
		}
	}
}

void cIO::_freq()
{
	if (!IS_ASSIGNED(freq) &&
		!IS_ASSIGNED(filmaf) && !IS_ASSIGNED(filmac) &&
		!IS_ASSIGNED(incmaf) && !IS_ASSIGNED(incmac) &&
		!IS_ASSIGNED(filgvar) && !IS_ASSIGNED(incgvar))
		return;

	char *Sp_freq = OPT_STRING(freq);
	if (!Sp_freq || Sp_freq[0] == '\0') {
		/* If --freq is off but --filmaf or --filmac assigned, do it with
		* --freq all */
		if (!IS_ASSIGNED(filmaf) && !IS_ASSIGNED(filmac) &&
			!IS_ASSIGNED(incmaf) && !IS_ASSIGNED(incmac) &&
			!IS_ASSIGNED(filgvar) && !IS_ASSIGNED(incgvar))
			return;
		Sp_freq = Z("all");
	}

	/* Option check */
	const char	*Sp_optFreq = Sp_freq;
	if (stricmp(Sp_freq, "all") && stricmp(Sp_freq, "founder") &&
		stricmp(Sp_freq, "blue"))
		halt("Invalid option parameter [%s] for --freq", Sp_optFreq);

/**/cExporter*	Cp_frq		= NULL;
	cExporter*	Cp_sampFrq	= NULL;
	if (stricmp(Sp_optFreq, "all") == 0) {
		Cp_frq = cExporter::summon("all.maf");
		LOGoutput("Variant MAF using all samples is exported to "
			"[%s.all.maf]\n", OPT_STRING(out));
	} else {
		Cp_frq = cExporter::summon("founders.maf");
		LOGoutput("Variant MAF using founder samples is exported to "
			"[%s.founders.maf]\n", OPT_STRING(out));
	}

	/* FIXME : Chromosome?? ???? N?? ?????? ???????? ??? */
	wsUint **Na_stat = NULL;
	if (stricmp(Sp_optFreq, "all")     == 0) {
		LOGoutput("Sample allele info is exported to [%s.sample.maf]\n",
			OPT_STRING(out));
		Cp_sampFrq = cExporter::summon("sample.maf");
		Cp_sampFrq->put("FID	IID	N_AA	N_Aa	N_aa	N_miss\n");
		wsAlloc(Na_stat, wsUint *, N_sample);
		for (wsUint i=0 ; i<N_sample ; i++) {
			wsCalloc(Na_stat[i], wsUint, 4);
		}
	}
	
	if (IS_ASSIGNED(annogene))
		Cp_frq->put("CHR	VARIANT	ANNOT	MAJOR	MINOR	MAF	MAC	FMAC	NIND\n");
	else
		Cp_frq->put("CHR	VARIANT	MAJOR	MINOR	MAF	MAC	FMAC	NIND\n");

	/* Forcely do computation */
	(void)getMAF();
	/* Generate histogram */
#ifdef USE_R
	wsVec Ra_mafs = sseVector(N_variant);
	LOOP (i, N_variant) Ra_mafs[i] = Xp_maf[i].R_maf;
	hist(stricmp(Sp_optFreq, "all") ? "MAF of all samples" : "MAF of founders",
		"maf", "Minor Allele Frequency (MAF)", Ra_mafs, N_variant);
	sseFree(Ra_mafs);
#endif

	if (stricmp(Sp_optFreq, "all") == 0) {
		LOOP (i, N_variant) {
			xVariant&	X_vrt	= Xv_variant[i];
			xMaf&		X_maf	= Xp_maf[i];
			wsStrCst		S_chr	= getChrName2(X_vrt.chr);
			if (IS_ASSIGNED(annogene)) {
				if (OPT_ENABLED(indel))
					Cp_frq->fmt("%s	%s	%s	%s	%s	%g	%d	%d	%d\n", S_chr, X_vrt.name,
						X_vrt.anno, X_vrt.indel1, X_vrt.indel2,
						X_maf.R_allMaf, X_maf.N_allMac, X_maf.N_mac, X_maf.N_allele>>1);
				else
					Cp_frq->fmt("%s	%s	%s	%c	%c	%g	%d	%d	%d\n", S_chr, X_vrt.name,
						X_vrt.anno, X_vrt.al1, X_vrt.al2,
						X_maf.R_allMaf, X_maf.N_allMac, X_maf.N_mac, X_maf.N_allele>>1);
			} else {
				if (OPT_ENABLED(indel))
					Cp_frq->fmt("%s	%s	%s	%s	%g	%d	%d	%d\n", S_chr,
						X_vrt.name, X_vrt.indel1, X_vrt.indel2,
						X_maf.R_allMaf, X_maf.N_allMac, X_maf.N_mac,
						X_maf.N_allele>>1);
				else
					Cp_frq->fmt("%s	%s	%c	%c	%g	%d	%d	%d\n", S_chr,
						X_vrt.name, X_vrt.al1, X_vrt.al2, X_maf.R_allMaf,
						X_maf.N_allMac, X_maf.N_mac, X_maf.N_allele>>1);
			}
			/* Set maf to allmaf */
			X_maf.R_maf = X_maf.R_allMaf;
		}

		for (wsUint i=0 ; i<N_sample ; i++) {
			Cp_sampFrq->fmt("%s	%s	%d	%d	%d	%d\n",
				Xa_sampleV2[i]->S_FID.c_str(), Xa_sampleV2[i]->S_IID.c_str(),
				Na_stat[i][0], Na_stat[i][1], Na_stat[i][2], Na_stat[i][3]);
			DEALLOC(Na_stat[i]);
		}
		delete Cp_sampFrq;
		DEALLOC(Na_stat);
	} else {
		wsUint		PN		= 0;
		/* Check the number of founders */
		for (wsUint j=0 ; j<N_sample ; j++)
			if (Ba_isFounder[j])
				PN++;
		if (PN != N_founder) halt("SYSERR: # of founders from "
			"founder indicator[%d] unmatch! Expected [%d]", PN,
			N_founder);

		LOOP (i, N_variant) {
			xVariant&	X_vrt	= Xv_variant[i];
			xMaf&		X_maf	= Xp_maf[i];
			wsStrCst		S_chr	= getChrName2(X_vrt.chr);

			if (IS_ASSIGNED(annogene)) {
				if (OPT_ENABLED(indel))
					Cp_frq->fmt("%s	%s	%s	%s	%s	%g	%d	%d	%d\n", S_chr,
						X_vrt.name, X_vrt.anno, X_vrt.indel1, X_vrt.indel2,
						X_maf.R_maf, X_maf.N_allMac, X_maf.N_mac, 
						X_maf.N_allele>>1);
				else
					Cp_frq->fmt("%s	%s	%s	%c	%c	%g	%d	%d	%d\n", S_chr,
						X_vrt.name, X_vrt.anno, X_vrt.al1, X_vrt.al2,
						X_maf.R_maf, X_maf.N_allMac, X_maf.N_mac,
						X_maf.N_allele>>1);
			} else {
				if (OPT_ENABLED(indel))
					Cp_frq->fmt("%s	%s	%s	%s	%g	%d	%d	%d\n", S_chr,
						X_vrt.name, X_vrt.indel1, X_vrt.indel2, X_maf.R_maf,
						X_maf.N_allMac, X_maf.N_mac, X_maf.N_allele>>1);
				else
					Cp_frq->fmt("%s	%s	%c	%c	%g	%d	%d	%d\n", S_chr,
						X_vrt.name, X_vrt.al1, X_vrt.al2, X_maf.R_maf,
						X_maf.N_allMac, X_maf.N_mac, X_maf.N_allele>>1);
			}
		} /* END of N_variant LOOP */
	}
	delete Cp_frq;

	/* Do --filmaf/--incmaf if --freq is not corr */
	if (IS_ASSIGNED(filmaf) && stricmp(Sp_optFreq, "blue"))
		N_szFilteredVrt += buildVariantFilterByMAF(OPT_RANGE(filmaf), 1,
			Bv_filtVrt, N_variant, Xv_variant, Xp_maf);
	else if (IS_ASSIGNED(incmaf) && stricmp(Sp_optFreq, "blue"))
		N_szFilteredVrt += buildVariantFilterByMAF(OPT_RANGE(incmaf), 0,
			Bv_filtVrt, N_variant, Xv_variant, Xp_maf);

	/* Do --filmac/--incmac */
	if (IS_ASSIGNED(filmac))
		N_szFilteredVrt += buildVariantFilterByMAC(OPT_RANGE(filmac), 1,
			Bv_filtVrt, N_variant, Xv_variant, Xp_maf);
	else if (IS_ASSIGNED(incmac))
		N_szFilteredVrt += buildVariantFilterByMAC(OPT_RANGE(incmac), 0,
			Bv_filtVrt, N_variant, Xv_variant, Xp_maf);

	/* Do --filgvar/--incgvar */
	if (IS_ASSIGNED(filgvar))
		N_szFilteredVrt += buildVariantFilterByGRATE(OPT_RANGE(filgvar), 1,
			Bv_filtVrt, N_variant, Xv_variant, N_sample);
	else if (IS_ASSIGNED(incgvar))
		N_szFilteredVrt += buildVariantFilterByGRATE(OPT_RANGE(incgvar), 0,
			Bv_filtVrt, N_variant, Xv_variant, N_sample);

	if (N_szFilteredVrt != 0)
		LOG("[%d/%d] variants will be removed from dataset\n", N_szFilteredVrt,
		N_variant);
}

/*
// This function implements an exact variant test of Hardy-Weinberg
// Equilibrium as described in Wigginton, JE, Cutler, DJ, and
// Abecasis, GR (2005) A Note on Exact Tests of Hardy-Weinberg
// Equilibrium. American Journal of Human Genetics. 76: 000 - 000
//
// Written by Jan Wigginton
*/

double variantHWE(int obs_hets, int obs_hom1, int obs_hom2)
{
  
  if (obs_hom1 + obs_hom2 + obs_hets == 0 ) return 1;
  
  if (obs_hom1 < 0 || obs_hom2 < 0 || obs_hets < 0)
	 halt("SYSERR: negative count in HWE test (%d/%d/%d)",
		 obs_hets, obs_hom1, obs_hom2);

  int obs_homc = obs_hom1 < obs_hom2 ? obs_hom2 : obs_hom1;
  int obs_homr = obs_hom1 < obs_hom2 ? obs_hom1 : obs_hom2;

  int rare_copies = 2 * obs_homr + obs_hets;
  int genotypes   = obs_hets + obs_homc + obs_homr;

  double * het_probs = NULL;
  wsAlloc(het_probs, double, (rare_copies + 1));
  
  int i;
  for (i = 0; i <= rare_copies; i++)
    het_probs[i] = 0.0;

  /* start at midpoint */
  int mid = rare_copies * (2 * genotypes - rare_copies) / (2 * genotypes);
  
  /* check to ensure that midpoint and rare alleles have same parity */
  if ((rare_copies & 1) ^ (mid & 1))
    mid++;
  
  int curr_hets = mid;
  int curr_homr = (rare_copies - mid) / 2;
  int curr_homc = genotypes - curr_hets - curr_homr;

  het_probs[mid] = 1.0;
  double sum = het_probs[mid];
  for (curr_hets = mid; curr_hets > 1; curr_hets -= 2)
    {
      het_probs[curr_hets - 2] = het_probs[curr_hets] * curr_hets * (curr_hets - 1.0)
	/ (4.0 * (curr_homr + 1.0) * (curr_homc + 1.0));
      sum += het_probs[curr_hets - 2];

      /* 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote */
      curr_homr++;
      curr_homc++;
    }

  curr_hets = mid;
  curr_homr = (rare_copies - mid) / 2;
  curr_homc = genotypes - curr_hets - curr_homr;
  for (curr_hets = mid; curr_hets <= rare_copies - 2; curr_hets += 2)
    {
      het_probs[curr_hets + 2] = het_probs[curr_hets] * 4.0 * curr_homr * curr_homc
	/((curr_hets + 2.0) * (curr_hets + 1.0));
      sum += het_probs[curr_hets + 2];
      
      /* add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote */
      curr_homr--;
      curr_homc--;
    }
  
  for (i = 0; i <= rare_copies; i++)
    het_probs[i] /= sum;

  /* alternate p-value calculation for p_hi/p_lo
   double p_hi = het_probs[obs_hets];
   for (i = obs_hets + 1; i <= rare_copies; i++)
     p_hi += het_probs[i];
   
   double p_lo = het_probs[obs_hets];
   for (i = obs_hets - 1; i >= 0; i--)
      p_lo += het_probs[i];
   
   double p_hi_lo = p_hi < p_lo ? 2.0 * p_hi : 2.0 * p_lo;
  */

  double p_hwe = 0.0;
  /*  p-value calculation for p_hwe  */
  for (i = 0; i <= rare_copies; i++)
    {
      if (het_probs[i] > het_probs[obs_hets])
	continue;
      p_hwe += het_probs[i];
    }
   
  p_hwe = p_hwe > 1.0 ? 1.0 : p_hwe;

  DEALLOC(het_probs);

  return p_hwe;
}

void cIO::_hwe()
{
	if (!IS_ASSIGNED(hwe) && !IS_ASSIGNED(filhwe) && !IS_ASSIGNED(inchwe))
		return;
	char	*Sp_optHWE = OPT_STRING(hwe);
	if (!Sp_optHWE || Sp_optHWE[0] == '\0') {
		if (IS_ASSIGNED(filhwe) || IS_ASSIGNED(inchwe)) {
			Sp_optHWE = strdup("founder");
			LOGwarn("--hwe is not defined, assumed to be 'founder'\n");
		} else
			return;
	}
/**/cExporter*	Cp_hwe = NULL;

	if (stricmp(Sp_optHWE, "all") == 0) {
		Cp_hwe = cExporter::summon("all.hwe");
		LOGoutput("Computed HWEs from all samples are exported to "
			"[%s.all.hwe]\n", OPT_STRING(out));
	} else if (stricmp(Sp_optHWE, "founder") == 0) {
		Cp_hwe = cExporter::summon("founders.hwe");
		LOGoutput("Computed HWEs from only founder samples are exported to "
			"[%s.founders.hwe]\n", OPT_STRING(out));
	}

	/* FIXME : Chromosome?? ???? N?? ?????? ???????? ??? */
	headerVariant(Cp_hwe);
	Cp_hwe->put("	GENO	O(HET)	E(HET)	STAT_HWE	P_EXACTHWE	P_ASYMHWE	NIND\n");

	if        (stricmp(Sp_optHWE, "all")     == 0) {
//		wsUint *Na_res = NULL;
		// c0 c1 c2 N stat
//		MULTI_MALLOC(Na_res, wsUint, sizeof(wsUint)*5);
//		WORKER().run(hweCalc, forAllSNP_equal, this, Na_res, sizeof(wsUint)*5);
		for (wsUint i=0 ; i<N_variant ; i++) {
			wsUint		N		= 0;
			wsUint		S		= 0;
			xVariant&	X_vrt	= Xv_variant[i];
			wsUint	N_occur[3] = { 0, 0, 0 };

			/* Get MA count */
			if (isDataComplete()) {
				for (wsUint j=0 ; j<N_sample ; j++) {
					S += Na_geno[j][i];
					N_occur[(wsUint)Na_geno[j][i]]++;
				}
				N = N_sample;
			} else {
				for (wsUint j=0 ; j<N_sample ; j++) {
					int N_genotype = Na_geno[j][i];

					if (isAvailable(N_genotype)) {
						S += N_genotype;
						N_occur[N_genotype]++;
						N++;
					}
				}
			}

			wsReal	R_maf = S/((wsReal)N*W2);

			wsReal R_q = N*R_maf*R_maf;
			wsReal R_p = N*(W1-R_maf)*(W1-R_maf);
			wsReal R_2pq = N*W2*R_maf*(W1-R_maf);

			/* Get chi^2 statistics */
			wsReal R_chistat =
				(N_occur[0] - R_p)*(N_occur[0] - R_p)/R_p +
				(N_occur[1] - R_2pq)*(N_occur[1] - R_2pq)/R_2pq +
				(N_occur[2] - R_q)*(N_occur[2] - R_q)/R_q;

			wsReal R_pVala = NA(R_chistat) ? WISARD_NAN :
				(wsReal)PVchisq(R_chistat, W1);
			wsReal R_pVale = variantHWE(N_occur[1], N_occur[0], N_occur[2]);
			/* Export output */
			entryVariant(Cp_hwe, X_vrt);
			if (NA(R_pVale)) {
				if (NA(R_pVala))
					Cp_hwe->fmt("	%d/%d/%d	%g	%g	%g	NA	NA	%d\n",
						N_occur[2], N_occur[1], N_occur[0],
						N_occur[1]/(wsReal)N, W2*R_maf*(W1-R_maf),
						R_chistat, N);
				else
					Cp_hwe->fmt("	%d/%d/%d	%g	%g	%g	NA	%g	%d\n",
						N_occur[2], N_occur[1], N_occur[0],
						N_occur[1]/(wsReal)N, W2*R_maf*(W1-R_maf),
						R_chistat, R_pVala, N);
			} else if (NA(R_pVala)) Cp_hwe->fmt("	%d/%d/%d	%g	%g	%g	%g	NA	%d\n",
				N_occur[2], N_occur[1], N_occur[0],
				N_occur[1]/(wsReal)N, W2*R_maf*(W1-R_maf),
				R_chistat, R_pVale, N);
			else Cp_hwe->fmt("	%d/%d/%d	%g	%g	%g	%g	%g	%d\n",
				N_occur[2], N_occur[1], N_occur[0],
				N_occur[1]/(wsReal)N, W2*R_maf*(W1-R_maf),
				R_chistat, R_pVale, R_pVala, N);

			X_vrt.hwe = R_pVale;
		}
	} else if (stricmp(Sp_optHWE, "founder") == 0) {
		for (wsUint i=0 ; i<N_variant ; i++) {
			wsUint		N		= 0;
			wsUint		S		= 0;
			xVariant&	X_vrt	= Xv_variant[i];
			wsUint	N_occur[3] = { 0, 0, 0 };

			/* Get MA count */
			if (isDataComplete()) {
				for (wsUint j=0 ; j<N_sample ; j++) {
					if (Ba_isFounder[j] == 0) continue;

					S += Na_geno[j][i];
					N_occur[(wsUint)Na_geno[j][i]]++;
					N++;
				}
			} else {
				for (wsUint j=0 ; j<N_sample ; j++) {
					if (Ba_isFounder[j] == 0) continue;
					int N_genotype = Na_geno[j][i];

					if (isAvailable(N_genotype)) {
						S += N_genotype;
						N_occur[N_genotype]++;
						N++;
					}
				}
			}

			wsReal	R_maf = S/((wsReal)N*W2);

			wsReal R_q = N*R_maf*R_maf;
			wsReal R_p = N*(W1-R_maf)*(W1-R_maf);
			wsReal R_2pq = N*W2*R_maf*(W1-R_maf);

			/* Get chi^2 statistics */
			wsReal R_chistat =
				(N_occur[0] - R_p)*(N_occur[0] - R_p)/R_p +
				(N_occur[1] - R_2pq)*(N_occur[1] - R_2pq)/R_2pq +
				(N_occur[2] - R_q)*(N_occur[2] - R_q)/R_q;

			wsReal R_pVala = R_q == W0 ? 1.0 : (wsReal)PVchisq(R_chistat, W1);
			double R_pVale = variantHWE(N_occur[1], N_occur[0], N_occur[2]);

			/* Export output */
			entryVariant(Cp_hwe, X_vrt);
			if (NA(R_pVale)) {
				if (NA(R_pVala))
					Cp_hwe->fmt("	%d/%d/%d	%g	%g	%g	NA	NA	%d\n",
					N_occur[2], N_occur[1], N_occur[0],
					N_occur[1]/(wsReal)N, W2*R_maf*(W1-R_maf),
					R_chistat, N);
				else
					Cp_hwe->fmt("	%d/%d/%d	%g	%g	%g	NA	%g	%d\n",
					N_occur[2], N_occur[1], N_occur[0],
					N_occur[1]/(wsReal)N, W2*R_maf*(W1-R_maf),
					R_chistat, R_pVala, N);
			} else if (NA(R_pVala)) Cp_hwe->fmt("	%d/%d/%d	%g	%g	%g	%g	NA	%d\n",
				N_occur[2], N_occur[1], N_occur[0],
				N_occur[1]/(wsReal)N, W2*R_maf*(W1-R_maf),
				R_chistat, R_pVale, N);
			else Cp_hwe->fmt("	%d/%d/%d	%g	%g	%g	%g	%g	%d\n",
				N_occur[2], N_occur[1], N_occur[0],
				N_occur[1]/(wsReal)N, W2*R_maf*(W1-R_maf),
				R_chistat, R_pVale, R_pVala, N);

			X_vrt.hwe = (wsReal)R_pVale;
		}
//	} else if (stricmp(Sp_optHWE, "blue") == 0) {
		/* Calculate later */
//		LOG("--freq=corr will calculate p^ value after getting correlation structure, such as correlation/kinship coefficient\n");
	} else
		halt("Invalid option parameter [%s] for --freq", Sp_optHWE);
	delete Cp_hwe;

	/* Do --filhwe if --hwe is not corr */
	if (IS_ASSIGNED(filhwe) && stricmp(Sp_optHWE, "blue")) {
		for (wsUint i=0 ; i<N_variant ; i++) {
			xVariant& X_vrt = Xv_variant[i];

			/* Should filtered */
			if (isInRange(OPT_RANGE(filhwe), X_vrt.hwe) && !Bv_filtVrt[i]) {
				N_szFilteredVrt++;
				Bv_filtVrt[i] = MFR_HWE;
			}
		}
	} else if (IS_ASSIGNED(inchwe) && stricmp(Sp_optHWE, "blue")) {
		for (wsUint i=0 ; i<N_variant ; i++) {
			xVariant& X_vrt = Xv_variant[i];

			/* Should filtered */
			if (!isInRange(OPT_RANGE(inchwe), X_vrt.hwe) && !Bv_filtVrt[i]) {
				N_szFilteredVrt++;
				Bv_filtVrt[i] = MFR_HWE;
			}
		}
	}
}

void _mistest_sub(wsUint N_samp, xVariant &X_m, wsUint I, char** Na_data,
	wsUint N_phe, wsMat Ra_phe, wsReal *Rp_ret=NULL)
{
	wsUint N_11 = 0, N_12 = 0, N_21 = 0, N_22 = 0;
	wsUint *Na_11 = NULL;
	wsUint *Na_12 = NULL;
	wsUint *Na_21 = NULL;
	wsUint *Na_22 = NULL;

	if (N_phe == 1) {
		/* Get and return p-value */
		for (wsUint i=0 ; i<N_samp ; i++) {
			if (isMissingReal(Ra_phe[0][i])) continue;
			char N_geno = Na_data[i][I];

			if (isMissing(N_geno))
				Ra_phe[0][i] == OPT_NUMBER(phenoCase) ? N_11++ : N_12++;
			else
				Ra_phe[0][i] == OPT_NUMBER(phenoCase) ? N_21++ : N_22++;
		}

		Rp_ret[0] = fisher(N_11, N_12, N_21, N_22);
	} else {
		wsCalloc(Na_11, wsUint, N_phe);
		wsCalloc(Na_12, wsUint, N_phe);
		wsCalloc(Na_21, wsUint, N_phe);
		wsCalloc(Na_22, wsUint, N_phe);
		for (wsUint i=0 ; i<N_samp ; i++) {
			char N_geno = Na_data[i][I];

			if (isMissing(N_geno)) for (wsUint j=0 ; j<N_phe ; j++) {
				if (isMissingReal(Ra_phe[j][i])) continue;
				Ra_phe[j][i] == OPT_NUMBER(phenoCase) ? Na_11[j]++ : Na_12[j]++;
			} else for (wsUint j=0 ; j<N_phe ; j++) {
				if (isMissingReal(Ra_phe[j][i])) continue;
				Ra_phe[j][i] == OPT_NUMBER(phenoCase) ? Na_21[j]++ : Na_22[j]++;
			}
		}
		/* Get the ith p-value */
		for (wsUint i=0 ; i<N_phe ; i++)
			Rp_ret[i] = fisher(Na_11[i], Na_12[i], Na_21[i], Na_22[i]);
	}
}

void cIO::_mistest()
{
//	wsUint	N_mkr	= getMarkerSize();
	vVariant	&Xa_mkr	= getVariant();
	char**	Na_data	= getGenotype();
	wsUint	N_phe	= sizePheno();
	vPheno&	Xv_phe	= getPhenoInfo();
	wsMat	Ra_phe	= getPhenos();
	if (!OPT_ENABLED(mistest) && !IS_ASSIGNED(filmistest) && !IS_ASSIGNED(incmistest))
		return;
	/* Do not perform if there is no phenotype or have continuous phenotype */
	if (N_phe == 0) {
		LOGwarn("--mistest is not performed because there is no phenotype\n");
		return;
	}
	for (wsUint i=0 ; i<N_phe ; i++)
		if (isContinuous(i)) {
			LOGwarn("[%d]th phenotype [%s] is not dichotomous, --mistest is"
				" not performed\n", i+1, Xv_phe[i].S_name.c_str());
			return;
		}

	/* Init exporter */
	cExporter *Cp_mt = cExporter::summon("mistest.res");
	LOGoutput("Result of variant-level missingness test by dichotomous trait "
		"is exported to [%s.mistest.res]\n", OPT_STRING(out));
	/* Print header */
	headerVariant(Cp_mt);
	FOREACH (vPheno_it, Xv_phe, i)
		Cp_mt->fmt("	%s", i->S_name.c_str());
	Cp_mt->put("\n");

	/* DO test */
	wsUint I=0;
	wsReal *Ra_ret = sseVector(N_phe);
	if (N_phe == 1) FOREACHDO (vVariant_it, Xa_mkr, i, I++) {
		_mistest_sub(N_sample, *i, I, Na_data, N_phe, Ra_phe, Ra_ret);

		/* Print output */
		if (NA(Ra_ret[0])) {
			if (!OPT_ENABLED(remna)) {
				entryVariant(Cp_mt, *i);
				Cp_mt->fmt("	NA\n");
			}
		} else if (!IS_ASSIGNED(pvalrange) ||
			(IS_ASSIGNED(pvalrange) && isInRange(OPT_RANGE(pvalrange), Ra_ret[0]))) {
			entryVariant(Cp_mt, *i);
			Cp_mt->fmt("	%g\n", Ra_ret[0]);
		}

		/* Do filtering */
		if (!Bv_filtVrt[I]) {
			if ((IS_ASSIGNED(filmistest) && isInRange(OPT_RANGE(filmistest), Ra_ret[0])) ||
				(IS_ASSIGNED(incmistest) && !isInRange(OPT_RANGE(incmistest), Ra_ret[0]))) {
				N_szFilteredVrt++;
				Bv_filtVrt[I] = MFR_MISCACT;
			}
		}
	} else {
		FOREACHDO (vVariant_it, Xa_mkr, i, I++) {
			_mistest_sub(N_sample, *i, I, Na_data, N_phe, Ra_phe, Ra_ret);

			/* Print output */
			if (OPT_ENABLED(remna)) {
				/* Check all na */
				wsUint j=0;
				for (j=0 ; j<N_phe ; j++)
					if (!NA(Ra_ret[j])) break;
				if (j == N_phe) continue;
			}
			if (IS_ASSIGNED(pvalrange)) {
				/* Check at least one pass */
				wsUint j=0;
				for (j=0 ; j<N_phe ; j++)
					if (isInRange(OPT_RANGE(pvalrange), Ra_ret[j]))
						break;
				if (j == N_phe) continue;
			}
			entryVariant(Cp_mt, *i);
			for (wsUint j=0 ; j<N_phe ; j++)
				if (NA(Ra_ret[j]))
					Cp_mt->put("	NA");
				else
					Cp_mt->fmt("	%g", Ra_ret[j]);
			Cp_mt->put("\n");

			/* Do filtering */
			if (!Bv_filtVrt[I]) {
				if ((IS_ASSIGNED(filmistest) && isInRangeOR(OPT_RANGE(filmistest), Ra_ret, N_phe)) ||
					(IS_ASSIGNED(incmistest) && !isInRangeOR(OPT_RANGE(incmistest), Ra_ret, N_phe))) {
					N_szFilteredVrt++;
					Bv_filtVrt[I] = MFR_MISCACT;
				}
			}
		}
	}

	delete Cp_mt;
}

wsMat cIO::getCovariates()
{
	return Ra_cov;
}

wsMat cIO::getGxEcovariates()
{
	return Ra_covGxE ? Ra_covGxE : Ra_cov;
}

const wsUint cIO::sizeCovar()
{
	return N_cov;
}

const wsUint cIO::sizeGxEcovar()
{
	return OPT_ENABLED(gxe) ? (Ra_covGxE ? N_covGxE : N_cov) : 0;
}

wsMat cIO::getFiltNullDesignMat(const char *Ba_filt, char N_valFilt,
	wsUint N_szEst, wsUint N_moreRow/*=0*/, wsMat *Ra_dm/*=NULL*/)
{
	wsUint	N_cova	= sizeCovar();
	wsUint	N_samp	= sizeSample();
	wsMat	Ra_ret	= sseMatrix(N_cova+1+N_moreRow, N_szEst);
	wsUint	I		= 0;

	if (Ra_dm) {
		/* Return non-transposed also */
		*Ra_dm = sseMatrix(N_szEst, N_cova+1);

		for (wsUint i=0 ; i<N_samp ; i++) {
			if (Ba_filt[i] == N_valFilt) continue;

			Ra_ret[0][I] = (*Ra_dm)[I][0] = W1;
			for (wsUint j=0 ; j<N_cova ; j++)
				Ra_ret[j+1][I] = (*Ra_dm)[I][j+1] = Ra_cov[j][i];
			I++;
		}
	} else {
		/* Return transposed only */
		for (wsUint i=0 ; i<N_samp ; i++) {
			if (Ba_filt[i] == N_valFilt) continue;

			Ra_ret[0][I] = W1;
			for (wsUint j=0 ; j<N_cova ; j++)
				Ra_ret[j+1][I] = Ra_cov[j][i];
			I++;
		}
	}

	if (I != N_szEst)
		halt("SYSERR : Built size[%d] is not match with expected size[%d]",
			I, N_szEst);

	return Ra_ret;
}

xMaf* cIO::getMAF()
{
	// 150824 Do not require it --expression
	if (IS_ASSIGNED(expression)) {
		LOGwarn("--expression do not require MAF computation\n");
		return NULL;
	}
	if (B_mafComputed) {
		if (!Xp_maf) halt("SYSERR: MAF computed but the pointer is NULL!");
		return Xp_maf;
	}
	if (Xp_maf) DEALLOC(Xp_maf);
	wsCalloc(Xp_maf, xMaf, N_variant);

	/* Update MAF */
	wsUint N_noIncCor = 0;
	if (N_variant) {
		cTimer c;
		wsUint N_noMinor = 0;
		LOG("MAF recalculation start...\n");

		/* --maf */
		c.start();
		if (IS_ASSIGNED(maf)) {
			cStrFile X_maf(OPT_STRING(maf), "Pre-defined MAF");
			char* S_buf = NULL;
			wsAlloc(S_buf, char, 1024);

			/* Pair of variant - MAF */
			char *a = NULL, *b = NULL;
			for (wsUint L=0 ; X_maf.gets(S_buf, 1024) ; L++) {
				getString(&S_buf, &a);
				if (a == NULL || !a[0])
					halt_fmt(WISARD_INVL_FILE_INVALID, "Pre-defined MAF",
						OPT_STRING(maf), S_buf, L+1);
				getString(&a, &b);
				Xm_newMAF.insert(make_pair(S_buf, atof(a)));
			}
			LOG("[%d] pre-defined MAFs loaded", Xm_newMAF.size());

			DEALLOC(S_buf);
		}

		WORKER().run(mafCalc, forAllVariant_equal, this, NULL, sizeof(int)*3);

		/* Check # of flips */
		wsUint N_flip = 0;
		LOOP (j, N_variant)
			if (Xp_maf[j].B_flipped) N_flip++;
		if (N_flip > 0)
			LOGnote("[%d] variants were flipped due to flipped MAF\n", N_flip);

//		150520 removed
// 		if (IS_ASSIGNED(dosage)) {
// 			/* Do not perform MAF-related thing w/ dosage */
// 			LOG("MAF calculation will not performed due to it is dosage data\n");
// 		} else
		LOOP (j, N_variant) {
			xVariant&	X_snv	= Xv_variant[j];
			xMaf&		X_maf	= Xp_maf[j];

			if (X_maf.R_maf == W0)
				N_noMinor++;
			/* Correlation incorporation determination */
			if (OPT_ENABLED(x) && X_snv.chr == (int)CHR_X) {
			} else if (!isInRange(OPT_RANGE(cormaf), X_maf.R_maf)) {
				//				Ra_corrMask[j] = W0;
				N_noIncCor++;
			} else if (!isAutosome(X_snv)) {
				//				if (!OPT_ENABLED(x) || OPT_ENABLED(x) && X_snv.chr != CHR_X)
				//				Ra_corrMask[j] = W0;
				//				N_noIncCor++;
			}
		}

		LOG("[%s] MAF recalculation end, %d variants having no minor allele found\n",
			c.getReadable(), N_noMinor);
	} else
		LOGwarn("No MAF recalculation performed due to the lack of variants\n");


	B_mafComputed = 1;
	return Xp_maf;
}

xMaf* cIO::getMAFbuffer()
{
	if (B_mafComputed || !Xp_maf) halt("SYSERR: Unmatched calling condition!");
	return Xp_maf;
}

void cIO::_buildTwin()
{
	/* STEP1 : twin ID should continuous */
	for (wsUint N_curTwinIdx = 1 ; N_curTwinIdx<=N_twin ; N_curTwinIdx++) {
		string		S_currFID = "";
		xSample		*Xp_p = NULL, *Xp_m = NULL;
		vSampPtr	Xa_currTwins;

		FOREACH (mSamp_it, Xm_iid2data, it) {
			xSample	&X_samp = it->second;

			/* Collect current loop's twins */
			if ((wsUint)(X_samp.N_idTwin) != N_curTwinIdx)
				continue;

			if (S_currFID == "") {
				S_currFID = X_samp.S_FID;
				Xp_m = X_samp.Xp_mat;
				Xp_p = X_samp.Xp_pat;
			}
			/* Sanity check */
			else if (S_currFID.compare(X_samp.S_FID))
				halt("Twin index %d have more than two FID [%s,%s...?]",
					N_curTwinIdx,
					S_currFID.c_str(), X_samp.S_FID.c_str());
			else if (Xp_m != X_samp.Xp_mat || Xp_p != X_samp.Xp_pat)
				halt("Twin index %d have inequal parents", N_curTwinIdx);

			Xa_currTwins.push_back(&X_samp);
		}
		/* Check twin index continues */
		if (Xa_currTwins.size() == 0)
			halt("Twin index %d have any twin, twin index should strictly continuous!",
				N_curTwinIdx);
		/* Register twin to its structure */
		FOREACH (vSampPtr_it, Xa_currTwins, it) {
			xSample *Xp_samp = *it;
			FOREACH (vSampPtr_it, Xa_currTwins, iit) 
				if (*iit != Xp_samp)
					Xp_samp->Xp_twins.push_back(*iit);
		}
	}
}

wsFmat cIO::getData(char B_retOnly/*=0*/)
{
//	if (OPT_ENABLED(indep)) return;
	if (Ra_data)
		return Ra_data;
	else if (B_retOnly)
		return NULL;
	if (!Cp_anaPPP)
		halt("SYSERR: PPP analysis required priorly!");

	Ra_data = NULL;
	wsAlloc(Ra_data, wsFvec, N_sample);
	for (wsUint i=0 ; i<N_sample ; i++) {
		sseCalloc(Ra_data[i], wsFloat, N_variant);
	}
	Cp_anaPPP->normalize(Ra_data);

	return Ra_data;
}

wsUint** cIO::getGenoProb()
{
	if (Na_gProb)
		return Na_gProb;
	return NULL;
}

void cIO::_initFilterings()
{
	N_filterVrt		= FILT_NOTHING; /* Default : Do not filtering */
	N_filterSample	= FILT_NOTHING; /* Default : Do not filtering */
	N_filterRange	= FILT_NOTHING; /* Default : Do not filtering */

	/* Initialize chr_allowance info */
	/* Make filter for chromosome if defined */
	memset(Ba_chrAllowed, 0x00, MAX_NCHR);
	if (IS_ASSIGNED(chr)) {
		char *a, *b = NULL;
		char *Sp_chr = OPT_STRING(chr);

		for (a = Sp_chr; a; a = b) {
			b = strchr(a, ',');
			if (b) {
				*b = '\0';
				b++;
			}
			/* Check is it ranged */
			char *p = strchr(a, '-');
			if (p) {
				*(p++) = '\0';
				int N_chrS = getChr(a);
				int N_chrE = getChr(p);
				if (N_chrS == -1 || N_chrE == -1)
					halt("[%s-%s] is invalid chromosome definition", a, p);
				for (int u = N_chrS; u <= N_chrE; u++)
					Ba_chrAllowed[u - 1] = 1;
			} else {
				int N_chr = getChr(a);
				if (N_chr == -1)
					halt("[%s] is invalid chromosome definition", a);
				Ba_chrAllowed[N_chr - 1] = 1;
			}
		}

		/* --autoonly checking */
		if (OPT_ENABLED(autoonly)) for (wsUint i=NAUTO_SPECIES ; i<NCHR_SPECIES ; i++)
			if (Ba_chrAllowed[i])
				halt_fmt(WISARD_CANT_DO_W_SOMEOPT,
				"Non-autosome selection with --chr", "--autoonly");
		/* --sexonly checking */
		if (OPT_ENABLED(autoonly)) for (wsUint i = 0; i < NAUTO_SPECIES; i++)
			if (Ba_chrAllowed[i])
				halt_fmt(WISARD_CANT_DO_W_SOMEOPT,
				"Autosome selection with --chr", "--sexonly");

		/* Print */
		LOG("Following chromosomes will included :");
		for (wsUint i = 1; i <= NCHR_SPECIES; i++) {
			if (Ba_chrAllowed[i - 1] == 0) continue;
			LOGnf(" chr%s", getChrName2(i));
		}
		LOGnf("\n");
	} else for (wsUint i = 0; i < NCHR_SPECIES; i++)
		/* Permit all chromosomes */
		Ba_chrAllowed[i] = 1;

	/* --varsubset? */
	if (IS_ASSIGNED(varsubset)) {
		/* File open */
		cStrFile H_vs(OPT_STRING(varsubset), "Variants subset file");
		char *Sp_buf = NULL;
		wsAlloc(Sp_buf, char, 1024);
		for (wsUint L=1 ; H_vs.gets(Sp_buf, 1024) ; L++) {
			char *tmp = NULL;
			getString(&Sp_buf, &tmp);
			if (!OPT_ENABLED(passemptyline) && !Sp_buf[0])
				halt("Line [%d] of variant subset file [%s] is invalid!",
					L, OPT_STRING(varsubset));
			Xe_varSubset.insert(Sp_buf);
		}
		DEALLOC(Sp_buf);
	}

	/* if --flip, then convert it to acgt */
	if (IS_ASSIGNED(flip)) {
		wsStrCst	Sa_flip	= OPT_STRING(flip);
		char*	Sp_dup	= strdup(Sa_flip);
		if (!Sp_dup[0] || !Sp_dup[1] || !Sp_dup[2] || !Sp_dup[3] || Sp_dup[4])
			halt("--flip value must be four characters, but the value is [%s]!", Sp_dup);
		OPTION().FORCE_OPT_STRING(acgt) = Sp_dup;
		LOG("--flip assigned [A->%c, C->%c, G->%c, T->%c]\n", Sp_dup[0],
			Sp_dup[1], Sp_dup[2], Sp_dup[3]);
	}

	/* --acgt */
	if (IS_ASSIGNED(acgt)) {
		char *Sa_acgt = OPT_STRING(acgt);
		wsCalloc(Na_ACGT, char, 0x80);
		/* Must be 4 chars */
		if (!Sa_acgt[0] || !Sa_acgt[1] || !Sa_acgt[2] || !Sa_acgt[3] || Sa_acgt[4])
			halt("Invalid --acgt definition, must be 4 ascii characters!");
		/* Should be 0x00~0x7f */
		if (Sa_acgt[0]<0 || Sa_acgt[1]<0 || Sa_acgt[2]<0 || Sa_acgt[3]<0)
			halt("Invalid --acgt definition, must be ascii character!");
		/* Define A */
		Na_ACGT[(wsUint)Sa_acgt[0]] = 'A';
		/* Define C */
		if (Na_ACGT[(wsUint)Sa_acgt[1]]) halt("C[%c] duplicates A in --acgt", Sa_acgt[1]);
		Na_ACGT[(wsUint)Sa_acgt[1]] = 'C';
		/* Define G */
		if (Na_ACGT[(wsUint)Sa_acgt[2]]) halt("G[%c] duplicates %c in --acgt",
			Sa_acgt[2], Na_ACGT[(wsUint)Sa_acgt[2]]);
		Na_ACGT[(wsUint)Sa_acgt[2]] = 'G';
		/* Define T */
		if (Na_ACGT[(wsUint)Sa_acgt[3]]) halt("T[%c] duplicates %c in --acgt",
			Sa_acgt[3], Na_ACGT[(wsUint)Sa_acgt[3]]);
		Na_ACGT[(wsUint)Sa_acgt[3]] = 'T';

		LOG("A/C/G/T assigned to the characters [%c] [%c] [%c] [%c]\n",
			Sa_acgt[0], Sa_acgt[1], Sa_acgt[2], Sa_acgt[3]);
	} else if (OPT_ENABLED(1234)) {
		wsCalloc(Na_ACGT, char, 0x80);
		Na_ACGT[(int)'1'] = 'A';
		Na_ACGT[(int)'2'] = 'C';
		Na_ACGT[(int)'3'] = 'G';
		Na_ACGT[(int)'4'] = 'T';
	} else
		Na_ACGT = NULL;

	/* Load marker info update list if assigned */
	if (IS_ASSIGNED(updvariant)) {
		Xv_altVariantDef = _loadVariantDef(OPT_STRING(updvariant));
		B_altVariantDef = 1;
	}

	/* Load filtering variant list if assigned */
	if (IS_ASSIGNED(remvariant) && IS_ASSIGNED(selvariant))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--remvariant", "--selvariant");
	else if (IS_ASSIGNED(remvariant))
		_loadVariantRemovalList();
	else if (IS_ASSIGNED(selvariant))
		_loadVariantSelectionList();

	/* Load range-filtering variant list if assigned */
	if (IS_ASSIGNED(filrange) && IS_ASSIGNED(incrange))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--filrange", "--incrange");
	else if (IS_ASSIGNED(filrange))
		_loadRangeRemovalList();
	else if (IS_ASSIGNED(incrange))
		_loadRangeSelectionList();
	else
		Xv_listRanges = NULL;

	/* Load filtering sample list if assigned */
	if (IS_ASSIGNED(remsamp) && IS_ASSIGNED(selsamp))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--remsamp", "--selsamp");
	else if (IS_ASSIGNED(remsamp) || IS_ASSIGNED(selsamp))
		_loadSampleList(IS_ASSIGNED(selsamp));
	if (IS_ASSIGNED(nasamp))
		_loadSampleNAList();

	/* Filtering direction should be equal */
	if ((IS_ASSIGNED(remfam) && IS_ASSIGNED(indsel)))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--remfam", "--indsel");
	if ((IS_ASSIGNED(selfam) && IS_ASSIGNED(remsamp)))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--selfam", "--remsamp");

	/* Load filtering family list if assigned */
	if (IS_ASSIGNED(remfam) && IS_ASSIGNED(selfam))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--selfam", "--remfam");
	else if (IS_ASSIGNED(remfam) || IS_ASSIGNED(selfam))
		_loadFamilyList(IS_ASSIGNED(selfam));
}

void def_xSample(xSample& X_dst, string& S_FID, string& S_IID, char N_sex)
{
	X_dst.S_FID = S_FID;
	X_dst.S_IID = S_IID;
	X_dst.N_sex = N_sex;

	/* In default, assume the sample is founderl */
	X_dst.Xp_mat = NULL;
	X_dst.Xp_pat = NULL;

	/* Undeterministic */
	X_dst.B_isComplete = 0;
	X_dst.B_isAssigned = 0;
	X_dst.B_isMissing = 1;
	X_dst.N_isVisited = 0;
	X_dst.N_idx = SAMP_NODATA;
	X_dst.N_oriIdx = SAMP_NODATA;
	X_dst.N_missGeno = 0;
	X_dst.N_idTwin = 0;
	X_dst.B_isProband = 0;
	X_dst.N_grpFst = -1;
}

void cIO::_initMemories(char B_noGeno/*=0*/)
{
	if (!B_noGeno) {
		wsAlloc(Na_geno, char*, N_sample);
		for (wsUint i=0 ; i<N_sample ; i++) {
			sseMalloc(Na_geno[i], char, N_variant);
		}
	}
	DEALLOC(Ba_isFounder);
	wsAlloc(Ba_isFounder, char, N_sample);
	for (wsUint i=0 ; i<N_sample ; i++)
		Ba_isFounder[i] = 0; // In default, all samples are considered as non-founder

	if (Ra_pheno)
		sseUnmat(Ra_pheno, N_pheno);
	if (IS_ASSIGNED(randpheno))
		N_pheno	= OPT_NUMBER(randpheno);
	else
		N_pheno	= 1;	// In default, single phenotype is assumed
	Ra_pheno	= sseMatrix(N_pheno, N_sample);
	if (Ba_typePheno) DEALLOC(Ba_typePheno);
	wsCalloc(Ba_typePheno, char, N_pheno);
	sseInit(Ra_pheno[0], N_sample, WISARD_NA_REAL);
}

xSample* cIO::_loadRealTimeFamData(char *Sp_buf, wsUint N_oriIdx)
{
	wsUint		L			= N_oriIdx+1;
	wsStrCst	S_misParent	= IS_ASSIGNED(misparent) ? OPT_STRING(misparent) : "0";
	string		S_dep		= Sp_buf;
	string		S_altIID;
	string		S_currFID, S_currIID;
	string		S_currPAT, S_currMAT;
	char		N_sex;
	char		*Sp_tmp1, *Sp_tmp2;
	static wsUint N	= 1; /* ID for --singleparent */
	char		B_filtered	= 0; /* In default, do not filter the sample */

_domore:

	/*
	 *
	 * Retrieve basic information from the buffer
	 *
	 */

	/*
	 * S_buf = FID (w/o --nofid)
	 */
	if (!OPT_ENABLED(nofid) && !IS_ASSIGNED(sepid)) {
		getString(&Sp_buf, &Sp_tmp1);
		S_currFID = (string)Sp_buf;
	} else Sp_tmp1 = Sp_buf;

	/*
	 * Sp_tmp1 = IID
	 */
	getString(&Sp_tmp1, &Sp_tmp2);
	/* In case of --sepid, separate IID by --sepid and [0]=FID, [1]=IID */
	if (IS_ASSIGNED(sepid)) {
		char* Sp_sep = strstr(Sp_tmp1, OPT_STRING(sepid));
		if (!Sp_sep)
			halt("IID [%s] does not have valid separator [%s]", S_currIID.c_str(), OPT_STRING(sepid));
		*Sp_sep = '\0';
		S_currFID = Sp_tmp1;
		S_currIID = Sp_sep + strlen(OPT_STRING(sepid));
	} else S_currIID = (string)Sp_tmp1;

	/* Replace IID if there is an alternative */
	if (S_altIID.length()) S_currIID = S_altIID;
	/* In case of --nofid/--ignorefid, assume FID is equal to IID */
	if (OPT_ENABLED(nofid) || OPT_ENABLED(ignorefid)) S_currFID = S_currIID;

	/* [FILTERING]
	 * Check this sample should be filtered by --samprem/--sampsel
	 */
	B_filtered = N_filterSample && isSampleFiltered(S_currFID,
		S_currIID, ((wsUint)Bv_filtSamp.size()-N_sample)>5 );

	/*
	 * Sp_tmp2 = PAT (w/o --noparent)
	 */
	if (!OPT_ENABLED(noparent)) {
		getString(&Sp_tmp2, &Sp_tmp1);
		S_currPAT = (string)Sp_tmp2;
		/* Sp_tmp1 = MAT */
		getString(&Sp_tmp1, &Sp_tmp2);
		S_currMAT = (string)Sp_tmp1;
	} else {
		/* Set to NO PARENT w/ --noparent */
		S_currPAT = S_misParent;
		S_currMAT = S_misParent;
	}
	/* Remove parent information if --ignoreparent is enabled */
	if (OPT_ENABLED(ignoreparent)) S_currMAT = S_currPAT = S_misParent;

	/*
	 * Sp_tmp2 = SEX (w/o --nosex)
	 */
	if (!OPT_ENABLED(nosex)) getString(&Sp_tmp2, &Sp_tmp1);
	else Sp_tmp1 = Sp_tmp2;
	/* Check sex validity */
	char N_chrMale = (char)OPT_NUMBER(sexMale);
	char N_chrFema = (char)OPT_NUMBER(sexFema);
	//LOG("male %d female %d\n", N_chrMale, N_chrFema);
	if (!OPT_ENABLED(nosex) && (Sp_tmp2[0] == N_chrMale || Sp_tmp2[0] == N_chrFema)
		&& Sp_tmp2[1] == '\0') {
		N_sex = Sp_tmp2[0];

		/* [FILTERING]
		 * Sex-related filtering
		 */
		if (N_sex == N_chrMale) {
			if (OPT_ENABLED(filmale)) B_filtered = 1;
			N_sex = 1;
		} else if (N_sex == N_chrFema) {
			if (OPT_ENABLED(filfemale)) B_filtered = 1;
			N_sex = 2;
		} else
			halt("SYSERR: Unhandled case for sex processing, value [%s]\n",
				Sp_tmp2);
	} else {
		/* Mark as MISSING */
		N_sex = WISARD_NA;

		/* Exclude this sample IF sex-related option is activated */
		if (OPT_ENABLED(autoonly)) {
			LOG("Sample [%s:%s] at line [%d] will be excluded because its "
				"sex is unknown[%s] by --autoonly\n",
				S_currFID.c_str(), S_currIID.c_str(), L, Sp_tmp2);
			B_filtered = 1;
		}

		/* Exclude this sample IF --filnosex */
		if (OPT_ENABLED(filnosex))
			B_filtered = 1;
	}

	/* Set the parent-related flags */
	char B_haveFather	= S_currPAT.compare(S_misParent) ? 1 : 0;
	char B_haveMother	= S_currMAT.compare(S_misParent) ? 1 : 0;
	char B_isFounder	= B_haveFather | B_haveMother ? 0 : 1;

	/* Check one-side parent */
	if (B_haveMother^B_haveFather) {
		if (OPT_ENABLED(singleparent)) {
			char S_tmp[512];
			sprintf(S_tmp, "WISARDTEMP_%d", ++N);
			if (B_haveFather) S_currMAT = S_tmp;
			else S_currPAT = S_tmp;
		} else
			halt_fmt(WISARD_CANT_ONESIDE_PARENT, S_currFID.c_str(),
				S_currIID.c_str(), S_currPAT.c_str(), S_currMAT.c_str());
	}

	/* 170711 --filnf added */
	if (OPT_ENABLED(filnf) && !B_isFounder)
		B_filtered = 1;

	/* Check three samples in this line (sample itself, father and mother)
	 * are previously inserted
	 */
	mSamp_it	X_resFind		= Xm_iid2data.find(S_currIID);
	mSamp_it	X_fail			= Xm_iid2data.end();
	mSamp_it	X_resPatFind, X_resMatFind;
	xSample		X_newFather;
	xSample		X_newMother;

	if (!B_isFounder) {
		X_resPatFind = Xm_iid2data.find(S_currPAT);
		X_resMatFind = Xm_iid2data.find(S_currMAT);

		/* If some parent is not exists, insert anyway */
		if (X_resPatFind == X_fail) {
			def_xSample(X_newFather, S_currFID, S_currPAT, 1);
			Xm_iid2data.insert(make_pair(S_currPAT, X_newFather));
			X_resPatFind = Xm_iid2data.find(S_currPAT);
		}
		if (X_resMatFind == X_fail) {
			def_xSample(X_newMother, S_currFID, S_currMAT, 2);
			Xm_iid2data.insert(make_pair(S_currMAT, X_newMother));
			X_resMatFind = Xm_iid2data.find(S_currMAT);
		}
	} else {
		X_resPatFind = X_fail;
		X_resMatFind = X_fail;
	}

	xSample*	Xp_thisSample = NULL; ///< Pointer of curr line's xSample structure
	xSample		X_mySamp;

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
				strcpy(Sp_buf, S_dep.c_str());
				goto _domore;
			}
		}

		/* Check sex and FID wrote on this line, and previously inserted info.
		 * These information should be identical
		 */
		if (N_sex != X_myData.N_sex || S_currFID.compare(X_myData.S_FID))
			halt("Sample [%s:%s (%s)] expected have [%s:%s (%s)], but not matched",
				S_currFID.c_str(), S_currIID.c_str(), N_sex==1?"Male":"Female",
				X_myData.S_FID.c_str(), X_myData.S_IID.c_str(), X_myData.N_sex==1?"Male":"Female");

		/* Assign */
		X_myData.B_isAssigned	= 1;
		X_myData.N_oriIdx		= N_oriIdx;
		X_myData.B_isMissing	= B_filtered;
		X_myData.Xp_mat			= B_isFounder ? NULL : &(X_resMatFind->second);
		X_myData.Xp_pat			= B_isFounder ? NULL : &(X_resPatFind->second);
		X_myData.B_isProband	= 0;
		X_myData.B_isComplete	= 0;
		X_myData.N_grpFst		= -1;

		Xp_thisSample = &X_myData;
	} else {
		/* This sample is currently new in xSample map */

		/* Fill up with basic information */
		X_mySamp.S_FID			= S_currFID;
		X_mySamp.S_IID			= S_currIID;
		X_mySamp.N_sex			= N_sex;
		X_mySamp.Xp_mat			= B_isFounder ? NULL : &(X_resMatFind->second);
		X_mySamp.Xp_pat			= B_isFounder ? NULL : &(X_resPatFind->second);
		X_mySamp.B_isAssigned	= 1;
		X_mySamp.B_isMissing	= B_filtered;
		X_mySamp.N_idx			= SAMP_NODATA;
		X_mySamp.N_oriIdx		= N_oriIdx;
		X_mySamp.N_missGeno		= 0;
		X_mySamp.N_idTwin		= 0;
		X_mySamp.B_isProband	= 0;
		X_mySamp.B_isComplete	= 0;
		X_mySamp.N_grpFst		= -1;

		Xp_thisSample = &X_mySamp;
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
	if (!B_isFounder) {
		/* Check paternal sex/FID */
		xSample& X_pat = X_resPatFind->second;
		if (X_pat.N_sex != 1 || X_pat.S_FID.compare(S_currFID))
			halt("Sample [%s:%s (%s,%d)] registered as the father of [%s:%s] but have unmatched information",
				X_pat.S_FID.c_str(), X_pat.S_IID.c_str(), X_pat.N_sex==1?"Male":"Female", X_pat.N_sex,
				S_currFID.c_str(), S_currIID.c_str());
		/* Check maternal sex/FID */
		xSample& X_mat = X_resMatFind->second;
		if (X_mat.N_sex != 2 || X_mat.S_FID.compare(S_currFID))
			halt("Sample [%s:%s (%s,%d)] registered as the mother of [%s:%s] but have unmatched information",
				X_mat.S_FID.c_str(), X_mat.S_IID.c_str(), X_mat.N_sex==1?"Male":"Female", X_mat.N_sex,
				S_currMAT.c_str(), S_currIID.c_str());

		/* Insert this sample as parent's child */
		X_pat.Xp_childs.push_back(&X_thisSample->second);
		X_mat.Xp_childs.push_back(&X_thisSample->second);

		/* Check this sample's mother has registered as the spouse of this sample's father */
		vSampPtr_it it;
		for (it=X_pat.Xp_spouses.begin() ; it!=X_pat.Xp_spouses.end() ; it++) {
			if ((*it)->S_IID.compare(S_currMAT) == 0)
				break;
		}
		/* And register if unregistered */
		if (it == X_pat.Xp_spouses.end())
			X_pat.Xp_spouses.push_back(&X_resMatFind->second);
		/* Do same to this sample's mother */
		for (it=X_mat.Xp_spouses.begin() ; it!=X_mat.Xp_spouses.end() ; it++) {
			if ((*it)->S_IID.compare(S_currPAT) == 0)
				break;
		}
		if (it == X_mat.Xp_spouses.end()) /* Not inserted */
			X_mat.Xp_spouses.push_back(&X_resPatFind->second);

		/* Insert to the family structure */
		mFam_it X_resFamFind = Xm_fid2data.find(S_currFID);

		/* This sample is the non-founder appeared at first in that family */
		if (X_resFamFind == Xm_fid2data.end()) {
			xFamily X_newFam;

			X_newFam.N_size = 1;
			X_newFam.S_FID = S_currFID;

			Xm_fid2data.insert(make_pair(S_currFID, X_newFam));
		} else {
			X_resFamFind->second.N_size++;
		}
	} else {
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

	return Xp_thisSample;
}

void cIO::_registMissingFounder()
{
	/* Register missing-founder to Xfamily */
	if (!OPT_ENABLED(filmf)) FOREACH (mSamp_it, Xm_iid2data, it) {
		xSample &X_samp = it->second;

		/* If missing founder */
		if (X_samp.Xp_mat == NULL && X_samp.Xp_pat == NULL &&
			X_samp.B_isMissing) {
			/* If founder, insert to family structure */
			mFam_it X_resFamFind = Xm_fid2data.find(X_samp.S_FID);

			/* This sample is the founder appeared at first */
			if (X_resFamFind == Xm_fid2data.end()) {
				xFamily X_newFam;

				X_newFam.N_size = 0;
				X_newFam.S_FID = X_samp.S_FID;
				X_newFam.Xp_founders.push_back(&X_samp);

				Xm_fid2data.insert(make_pair(X_samp.S_FID, X_newFam));
			} else
				X_resFamFind->second.Xp_founders.push_back(&X_samp);
		}
	} else {
		mSampPtr	Xm_ds;

		/* Remove all missing founders from available samples */
		FOREACH (mSamp_it, Xm_iid2data, i) {
			xSample	&X_s	= i->second;
			char	B_missM	= X_s.Xp_mat && X_s.Xp_mat->B_isMissing;
			char	B_missP	= X_s.Xp_pat && X_s.Xp_pat->B_isMissing;

			/* Promote to founder if satisfies */
			if (B_missM && B_missP) {
				xFamily &X_fam = Xm_fid2data[X_s.S_FID];
				X_fam.Xp_founders.push_back(&X_s);
				LOG("Sample [%s] promoted to founder because its founders are all missing\n",
					X_s.S_IID.c_str());
			}
			if (B_missM) {
//				Xm_ds.insert(make_pair(X_s.Xp_mat->S_IID, X_s.Xp_mat));
				X_s.Xp_mat->B_isMissing = -9;
				X_s.Xp_mat = NULL;
			}
			if (B_missP) {
//				Xm_ds.insert(make_pair(X_s.Xp_pat->S_IID, X_s.Xp_pat));
				X_s.Xp_pat->B_isMissing = -9;
				X_s.Xp_pat = NULL;
			}
		}

		for (mSamp_it i=Xm_iid2data.begin() ; i!=Xm_iid2data.end() ; ) {
			if (i->second.B_isMissing == -9) {
				/* Remove spouse info */
				FOREACH (vSampPtr_it, i->second.Xp_spouses, j) {
					for (vSampPtr_it k=(*j)->Xp_spouses.begin() ;
						k!=(*j)->Xp_spouses.end() ; )
						if ((*k) == &(i->second))
							k = (*j)->Xp_spouses.erase(k);
						else
							k++;
				}

				LOG("Missing founder [%s] erased\n", i->second.S_IID.c_str());
				Xm_iid2data.erase(i++);
			} else ++i;
		}

		/* Remove all missing foundrers */
// 		FOREACH (mFam_it, Xm_fid2data, i) {
// 			FOREACH (vSampPtr_it, i->second.Xp_founders, it) {
// 				if ((*it)->B_isMissing == -9) {
// 					LOG("Founder [%s] erased\n", (*it)->S_IID.c_str());
// 					it = i->second.Xp_founders.erase(it);
// 					it--;
// 				}
// 			}
// 		}
	}

	/* Export family structure */
	if (OPT_ENABLED(famsummary))
		exportFamily();

	/* Export family diagram */
	if (0) //OPT_ENABLED(famdiagram))
		exportFamDiagram();
}

char* cIO::getACGT()
{
	return Na_ACGT;
}

void _inverseNormalize(wsReal *Ra_vec, wsUint N_sz)
{
	/* Build */
	wsUint N_szAfter = 0;
	xRealSort *Xp_sort = buildRealSort(Ra_vec, N_sz, &N_szAfter);

	/* Sort */
	qsort(Xp_sort, N_szAfter, sizeof(xRealSort), sort_real);

	/* Now have rank */
	for (wsUint i=0 ; i<N_szAfter ; i++)
		Ra_vec[Xp_sort[i].i] = qnorm((0.5 + i) / N_szAfter);
}

void cIO::_finalize()
{
	/* Update genotype */
	if (IS_ASSIGNED(updgeno)) {
		cStrFile S(OPT_STRING(updgeno), "Genotype alternate file");
		char *B = new char[4096*4096];
		
		for (wsUint L=0 ; S.gets(B, 4096*4096) ; L++) {
			char *a=B, *b=NULL;
			for (wsUint M=0 ; a ; M++,a=b) {
				getString(&a, &b);
				Na_geno[L][M] = (char)atoi(a);
			}
		}

		LOG("Genotype successfully updated to [%s]\n", OPT_STRING(updgeno));
		delete [] B;
	}
// 	N_szFilteredSamp	= 0;
// 	N_szFilteredVrt		= 0;

	/* Update allele info */
	if (IS_ASSIGNED(updallele)) {
		vVariant& Xv_vrt = getVariant();

		/* Load update allele info */
		cStrFile	X_ua(OPT_STRING(updallele), "Allele update information");
		char		*S_buf = NULL;
		wsAlloc(S_buf, char, 1024);
		mUAinfo		Xm_ua;

		if (!OPT_ENABLED(indel)) for (wsUint L=1 ; X_ua.gets(S_buf, 1024) ; L++) {
			xUpdAlleleInfo	X_ua;
			string			S_n;
			char *a = NULL, *b = NULL;
			getString(&S_buf, &a);	// S_buf = rsid
			if (!a) goto _stop;

			S_n = S_buf;
			getString(&a, &b);		// a = o1
			if (a[1]) goto _stop;
			X_ua.o1 = a[0];
			if (!b) goto _stop;

			getString(&b, &a);		// b = o2
			if (b[1]) goto _stop;
			X_ua.o2 = b[0];
			if (!a) goto _stop;

			getString(&a, &b);		// a = u1
			if (a[1]) goto _stop;
			X_ua.u1 = a[0];
			if (!b) goto _stop;

			getString(&b, &a);		// b = u2
			if (b[1]) goto _stop;
			X_ua.u2 = b[0];

			X_ua.O1 = X_ua.O2 = X_ua.U1 = X_ua.U2 = NULL;

			/* Insert data */
			Xm_ua.insert(make_pair(S_n, X_ua));
			continue;
_stop:
			halt_fmt(WISARD_INVL_FILE_INVALID_DESC, "allele update information",
				OPT_STRING(updallele), "number of columns",
				"rsid	orig_al1	orig_al2	upd_al1	upd_al2", L);
		} else for (wsUint L=1 ; X_ua.gets(S_buf, 1024) ; L++) {
			xUpdAlleleInfo	X_ua;
			string			S_n;
			char *a = NULL, *b = NULL;
			getString(&S_buf, &a);	// S_buf = rsid
			if (!a) goto _stop2;

			S_n = S_buf;
			getString(&a, &b);		// a = o1
			X_ua.O1 = strdup(a);
			if (!b) goto _stop2;

			getString(&b, &a);		// b = o2
			X_ua.O2 = strdup(b);
			if (!a) goto _stop2;

			getString(&a, &b);		// a = u1
			X_ua.U1 = strdup(a);
			if (!b) goto _stop2;

			getString(&b, &a);		// b = u2
			X_ua.U2 = strdup(b);

			X_ua.o1 = X_ua.o2 = X_ua.u1 = X_ua.u2 = '\0';

			/* Insert data */
			Xm_ua.insert(make_pair(S_n, X_ua));
			continue;
_stop2:
			halt_fmt(WISARD_INVL_FILE_INVALID_DESC, "allele update information",
				OPT_STRING(updallele), "number of columns",
				"rsid	orig_al1	orig_al2	upd_al1	upd_al2", L);
		}
		/* Free resources */
		Xm_ua.clear();
		DEALLOC(S_buf);
		X_ua.close();

		/* Now check variants */
		wsUint N_update = 0;
		FOREACH (vVariant_it, Xv_vrt, i) {
			/* Get entry */
			mUAinfo_it X_find = Xm_ua.find(i->name);
			if (X_find == Xm_ua.end()) continue;
			xUpdAlleleInfo &X_ua = X_find->second;

			N_update++;
			/* Check allele information */
			if (!OPT_ENABLED(indel)) {
				if (i->al1 == X_ua.o1) i->al1 = X_ua.u1;
				else if (i->al2 == X_ua.o1) i->al2 = X_ua.u1;
				if (i->al1 == X_ua.o2) i->al1 = X_ua.u2;
				else if (i->al2 == X_ua.o2) i->al2 = X_ua.u2;
				continue;
			}

			if (!stricmp(i->indel1, X_ua.O1)) {
				DEALLOC(i->indel1);
				i->indel1 = strdup(X_ua.U1);
			} else if (!stricmp(i->indel2, X_ua.O1)) {
				DEALLOC(i->indel2);
				i->indel2 = strdup(X_ua.U1);
			}
			if (!stricmp(i->indel1, X_ua.O2)) {
				DEALLOC(i->indel1);
				i->indel1 = strdup(X_ua.U2);
			} else if (!stricmp(i->indel2, X_ua.O2)) {
				DEALLOC(i->indel2);
				i->indel2 = strdup(X_ua.U2);
			}
		}
		LOG("--updallele updated [%d] variants\n", N_update);

		/* Deallocate if --indel */
		if (OPT_ENABLED(indel)) FOREACH (mUAinfo_it, Xm_ua, i) {
			if (i->second.O1) DEALLOC(i->second.O1);
			if (i->second.O2) DEALLOC(i->second.O2);
		}
	}

	/* Apply random missing (--simfam already have nageno) */
	if (IS_ASSIGNED(nageno) && stricmp(getFormat(), "sim.fam")) {
		INT64_t	N_cntGeno	=  0;
		if (OPT_REAL(nageno) < 1.0)
			N_cntGeno = (INT64_t)((double)N_sample*(double)N_variant * OPT_REAL(nageno));
		else
			N_cntGeno = (INT64_t)OPT_REAL(nageno);

		LOG("[" FMT_INT64 "] genotypes are masked to NA by --nageno\n", N_cntGeno);
		for (wsUint i=0 ; i<N_cntGeno ; i++) {
			wsUint k = wsRand()%N_sample;
			wsUint j = wsRand()%N_variant;
			if (!isMissing(Na_geno[k][j]))
				N_missGeno++;
			Na_geno[k][j] = WISARD_NA;
		}
	}

	/* varwindow */
	if (IS_ASSIGNED(varwindow)) {
		wsUint N_szWnd = OPT_NUMBER(varwindow);
		vVariant& Xv_variant = getVariant();

		/* Split by chromosomes */
		xUintSort**	Xa_vrts	= NULL;
		wsUint*		Na_mkr	= NULL;
		wsCalloc(Na_mkr, wsUint, MAX_NCHR+1);
		wsAlloc(Xa_vrts, xUintSort*, MAX_NCHR+1);
		/* Count each chromosome */
		wsUint i = 0;
		FOREACHDO(vVariant_it, Xv_variant, it, i++)
			Na_mkr[it->chr > 0 && (wsUint)it->chr <= NCHR_SPECIES ? it->chr : 0]++;
		/* Allocate memory */
		for (wsUint i=0 ; i<=NCHR_SPECIES ; i++)
			wsAlloc(Xa_vrts[i], xUintSort, Na_mkr[i]);
		/* Build index... now */
		i = 0;
		FOREACHDO(vVariant_it, Xv_variant, it, i++) {
			xUintSort *Xp_vrt = Xa_vrts[it->chr > 0 && (wsUint)it->chr <= NCHR_SPECIES ? it->chr : 0];
			Xp_vrt[i].i = i;
			Xp_vrt[i].V = it->pos;
		}

		/* For each chromosome, do sort */
		for (wsUint i=0 ; i<=NCHR_SPECIES ; i++)
			qsort(Xa_vrts[i], Na_mkr[i], sizeof(xStrSort),
				stricmp(OPT_STRING(sortpos), "asc") ? sort_uint_desc : sort_uint);

		/* By sorting, do picking */
		for (wsUint i=0 ; i<=NCHR_SPECIES ; i++) {
			wsUint N_mkr = Na_mkr[i];
			if (N_mkr) continue;
			wsUint N_lastPick = Xa_vrts[i][0].V;

			for (wsUint j=1 ; j<N_mkr ;j++) {
				if ((Xa_vrts[i][j].V - N_lastPick) <= N_szWnd)
					Bv_filtVrt[Xa_vrts[i][j].i] = MFR_WINDOW;
				else
					N_lastPick = Xa_vrts[i][j].V;
			}
		}
	}

	/* Anno */
	if (IS_ASSIGNED(annogene)) wsAnnot(this);

	/* --mendel related option */
	_md();
	/* --freq related option */
	_freq();
	/* --hwe related option */
	_hwe();

	/* --filgenic/--filintergenic */
	if (OPT_ENABLED(filgenic) || OPT_ENABLED(filintergenic)) {
		/* Have annotation? */
		if (wsAnnot(this) == NULL)
			halt_fmt(WISARD_CANT_DO_WO_SOMEOPT, "--filgenic/--filintergenic",
				"--annogene");
		wsUint j=0;
		if (OPT_ENABLED(filgenic)) FOREACHDO (vVariant_it, Xv_variant, i, j++) {
			if  (i->anno && !Bv_filtVrt[j]) {
				Bv_filtVrt[j] = MFR_RGNGENE;
				N_szFilteredVrt++;
			}
		} else if (OPT_ENABLED(filintergenic)) FOREACHDO (vVariant_it, Xv_variant, i, j++) {
			if (!i->anno && !Bv_filtVrt[j]) {
				Bv_filtVrt[j] = MFR_RGNGENE;
				N_szFilteredVrt++;
			}
		}
	}

	/* cVariantMap initialization */

	/* Load extra phenotype information if it assigned */
	if (IS_ASSIGNED(sampleweight)) {
		cStrFile	C_sw(OPT_STRING(sampleweight), "Sample weighting file", 1);
		if (!C_sw.isFailed()) {
			LOGnote("Retrieving sample weight from file [%s]\n", OPT_STRING(sampleweight));
			mSamp&		Xm_samp = getSampleData();
			vSampPtr&	Xv_samp = getSample();
			wsUint		N_samp	= sizeSample();
			char*		Sp_buf	= NULL;
			wsAlloc(Sp_buf, char, 4096);

			wsReal R_val = WISARD_NA;
			V_sampWgt.init(N_samp, NULL, &R_val);
			wsVec Ra_sampwgt = V_sampWgt.get();

			for (wsUint L=1 ; C_sw.gets(Sp_buf, 4096) ; L++) {
				char *a, *b, *c, *d;
				getString(&Sp_buf, &a);
				getString(&a, &b);
				getString(&b, &c);
				mSamp_it X_find = Xm_samp.find(a);
				if (X_find == Xm_samp.end()) continue;
				xSample&	X_samp = X_find->second;

				/* Get weight */
				wsReal R_wgt = (wsReal)str2dbl(b, &d);
				if (d && d[0]) halt("Invalid weight value [%s] in line [%d]",
					c, L);
				/* Weight can't be zero since it takes an inverse */
				if (R_wgt == W0)
					halt("Sample weight of [%s::%s] cannot be zero",
					Sp_buf, a);

				/* If the sample is not available, skip */
				if (X_samp.N_idx == SAMP_NODATA) continue;
				Ra_sampwgt[X_samp.N_idx] = W1/R_wgt;
			}
			DEALLOC(Sp_buf);
			/* Check all samples have their weight */
			for (wsUint i=0 ; i<N_samp ; i++)
				if (NA(Ra_sampwgt[i])) halt("Sample [%s::%s] does not have "
					"sample weight!", Xv_samp[i]->S_FID.c_str(),
					Xv_samp[i]->S_IID.c_str());

			LOGnote("Sample weights successfully loaded\n");
		} else if (IS_ASSIGNED(pheno)) {
			LOG("--pheno [%s] seems not a file, maybe in the sample variable file?",
				OPT_STRING(sampleweight));
		} else
			halt("INvalid sample weighting value [%s]", OPT_STRING(sampleweight));
	}
	if (IS_ASSIGNED(pheno))
		loadSampVar(this);
	if (IS_ASSIGNED(variantvar)) {
		vStr FUCK;
		_loadVariantVar(FUCK);
	}
	char B_haveAltPHENO = IS_ASSIGNED(pheno) && IS_ASSIGNED(pname);
	if (!B_haveAltPHENO && (OPT_ENABLED(filcase) || OPT_ENABLED(filcontrol)
		|| OPT_ENABLED(filmispheno))) {
		/* Determine filtering by --filcase/--filcontrol */
		for (wsUint i=0 ; i<N_sample ; i++) {
			if (Bv_filtSamp[i]) continue;
			
			if ((OPT_ENABLED(filmispheno) && isMissingReal(Ra_pheno[0][i])) ||
				(!isContinuous() && OPT_ENABLED(filcase) && Ra_pheno[0][i] == OPT_NUMBER(phenoCase)) ||
				(!isContinuous() && OPT_ENABLED(filcontrol) && Ra_pheno[0][i] == OPT_NUMBER(phenoCtrl))) {
				Bv_filtSamp[i] = 1;
				N_szFilteredSamp++;
			}
		}
	}
	/* Set default phenotype name */
	if ((IS_ASSIGNED(pheno) && !IS_ASSIGNED(pname)) || !IS_ASSIGNED(pheno)) {
      /* Set default */
		if (N_pheno == 1) {
			xPheno X_phe;
			X_phe.S_name = "PHENOTYPE";
			X_phe.N_idx  = 0;
			V_phenos.push_back(X_phe);
		} else for (wsUint i=0 ; i<N_pheno ; i++) {
			xPheno	X_phe;
			char	S_phen[256];
			sprintf(S_phen, "PHENOTYPE%d", i+1);
			X_phe.S_name = S_phen;
			X_phe.N_idx  = i;
			V_phenos.push_back(X_phe);
		}
	}

	/* --mistest related option */
	_mistest();

	/* Filtering samples by variable conditions */
	if (IS_ASSIGNED(filsample) && IS_ASSIGNED(incsample))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--filsample", "--incsample");
// 	else if (IS_ASSIGNED(filsample)) _filtVar(OPT_STRING(filsample), 1, 0);
// 	else if (IS_ASSIGNED(incsample)) _filtVar(OPT_STRING(incsample), 0, 0);

	/* If there is more than single phenotype, check option
	*  --filcase, --filcontrol */
	if (sizePheno() > 1 && (OPT_ENABLED(filcase) || OPT_ENABLED(filcontrol)))
		halt_fmt(WISARD_CANT_OPT_W_SOMESTATE, "--filcase/--filcontrol",
			"multiple phenotype assignment");
		//halt_fmt(WISARD_CANT_FILCACT_W_MULTIPHENO, getPhenoCount());
	if (N_altSample && N_altSample != N_sample) {
		vSampPtr	&Xa_samp	= getSample();
		//MAT_t		Ra_cov		= getCovariates();
		wsMat		Ra_pheno	= getPhenos();
		for (wsUint j=0 ; j<N_pheno ; j++)
			for (wsUint i=0 ; i<N_sample ; i++)
				if (Ra_pheno[j][i] != Ra_pheno[j][i])
					LOG("Sample [%s::%s] have incomplete phenotype data\n",
					Xa_samp[i]->S_FID.c_str(),
					Xa_samp[i]->S_IID.c_str());

		halt("# of samples in alternative phenotype file[%d] is not same with dataset[%d]!",
			N_altSample, N_sample);
	}

	/* Check sample filtering irf enabled */
	if (N_filterSample != FILT_NOTHING && mapIsExist(Xe_listSample, "__LIST__")) {
		/* If sample list was given list format, all of them should be
		* processed, otherwise error */
		FOREACH (eStr_it, Xe_listSample, it)
			if (it->compare("__LIST__"))
				halt("Sample filtering keyword [%s] looks invalid (stopped by [%s])",
					IS_ASSIGNED(remsamp)?OPT_STRING(remsamp):OPT_STRING(selsamp),
					it->c_str());
	}
	
	/* Check missingness and set as below:
	* Not checked as missing(0) -> 1
	* Checked as have missing(-1) -> 0
	*/
	{
		vSampPtr &Xa_samps = getSample();
		FOREACH (vSampPtr_it, Xa_samps, i) {
			if ((*i)->B_isComplete != 0 && (*i)->B_isComplete != -1)
				halt("Missingness of sample [%s:%s] have invalid value [%d]",
				(*i)->S_FID.c_str(), (*i)->S_IID.c_str(),
				(*i)->B_isComplete);
			(*i)->B_isComplete = !((*i)->B_isComplete);
		}
	}

	/* Register twin information */
	_buildTwin();

	/* if Xm_listSample, check which samples were not filtered even they are listed */
	cExporter *Cp_E = NULL;
	wsUint		N_missSampList = 0;
	FOREACH (eStr_it, Xe_listSample, it)
		if (it->compare("__LIST__")) {
			if (!Cp_E) {
				Cp_E = cExporter::summon("samplist.miss.lst");
				Cp_E->put("IID\n");
				LOGoutput("Samples in sample list but not in dataset are "
					"exported to [%s.samplist.miss.lst]\n", OPT_STRING(out));
			}
			LOG("Sample [%s] listed in sample list, but not exists in dataset\n",
				it->c_str());
			N_missSampList++;
			Cp_E->fmt("%s\n", it->c_str());
		}
	if (Cp_E) {
		LOG("[%d] samples from sample list are missing in dataset\n", N_missSampList);
		delete Cp_E;
	}

	/* Have filtering or sorting*/
	if (N_szFilteredVrt || N_szFilteredSamp ||
		IS_ASSIGNED(sortvariant) || IS_ASSIGNED(sortsample) ||
		IS_ASSIGNED(sampresize) || IS_ASSIGNED(nasamp) || IS_ASSIGNED(randnasamp) ||
		IS_ASSIGNED(sortiid) || IS_ASSIGNED(sortfid)) {
		INT64_t	N_cntGeno	= (INT64_t)N_sample*(INT64_t)N_variant;

		if (N_variant)
			LOG("Pre-resize genotyping rate : %.2f%% [" FMT_INT64 "/" FMT_INT64 "]\n",
				REAL_CONST(100.0)
				*(W1-(N_missGeno/((wsReal)N_sample*(wsReal)N_variant))),
				N_cntGeno-N_missGeno,
				N_cntGeno);
		else
			LOGwarn("No variant found in the dataset\n");

		_resizeData();
	}

	/* Phenotype variance check */
	char N_phenoIntegrity = 1;
	LOOP (i, N_pheno) {
		wsUint j = 1;
		wsReal R_v = Ra_pheno[i][0];
		for (j=1 ; j<N_sample ; j++)
			if (!isMissingReal(Ra_pheno[i][j]) && (R_v != Ra_pheno[i][j])) break;
		if (j == N_sample && !isMissingReal(R_v)) {
			vPheno& Xv_pheno = getPhenoInfo();
			LOGwarn("Phenotype [%s] consists of single value [%g], cannot continue!\n",
				Xv_pheno[i].S_name.c_str(), R_v);
			N_phenoIntegrity = 0;
		}
	}
	if (N_phenoIntegrity == 0)
		LOGwarn("At least one of phenotype have problem, please check it!");

	/* --napheno */
	if (IS_ASSIGNED(napheno)) {
		wsUint	N_cntPheno	=  0;
		if (OPT_REAL(napheno) < 1.0)
			N_cntPheno = (wsUint)((double)N_sample*(double)N_pheno * OPT_REAL(napheno));
		else
			N_cntPheno = (wsUint)OPT_REAL(napheno);

		LOG("[%d] phenotypes are masked to NA by --napheno\n", N_cntPheno);
		for (wsUint i=0 ; i<N_cntPheno ; i++) {
			wsUint k = wsRand()%N_sample;
			wsUint j = wsRand()%N_pheno;
// 			if (!isMissing(Ra_pheno[j][k]))
// 				N_missGeno++;
			Ra_pheno[j][k] = WISARD_NA;
		}
	}

	/* Check variant validity */
	if (OPT_ENABLED(markercheck)) {
		wsUint N_inv = 0;
		cExporter* Cp_ex = cExporter::summon("rs.invalid.lst");
		LOGoutput("rsID validation failed markers are exported to [%s.rs.invalid.lst]\n",
			OPT_STRING(out));
		int N_maxChr = (int)NCHR_SPECIES;
		vector<map<int,int>> Xv_map;
		Xv_map.resize(N_maxChr+1);

		/* Set the species */
		char S_species[64];
		if (!IS_ASSIGNED(species) || !stricmp(OPT_STRING(species), "human")) {
			strcpy(S_species, "human_9606");
		} else {
			if (!stricmp(OPT_STRING(species), "rice")) {
				strcpy(S_species, "rice_4530");
			} if (!stricmp(OPT_STRING(species), "mouse")) {
				strcpy(S_species, "mouse_10090");
			} else if (!stricmp(OPT_STRING(species), "rat")) {
				strcpy(S_species, "rat_10116");
			} else if (!stricmp(OPT_STRING(species), "rabbit")) {
				strcpy(S_species, "rat_9986");
			} else if (!stricmp(OPT_STRING(species), "sheep")) {
				strcpy(S_species, "sheep_9940");
			} else if (!stricmp(OPT_STRING(species), "cow")) {
				strcpy(S_species, "cow_9913");
			} else if (!stricmp(OPT_STRING(species), "horse")) {
				strcpy(S_species, "horse_9796");
			} else if (!stricmp(OPT_STRING(species), "dog")) {
				strcpy(S_species, "dog_9615");
			} else if (!stricmp(OPT_STRING(species), "pig")) {
				strcpy(S_species, "pig_9823");
			}
		}

		/* Fetch reference data and build the map */
		for (wsUint i=1 ; i<=NCHR_SPECIES ; i++) {
			char S_bedPath[512];
			char *Sp_buf = NULL;
			wsAlloc(Sp_buf, char, 1024);

			/* Check file existence */
			wsStrCst Sp_chr = getChrName2(i);
			sprintf(S_bedPath, "snp/organisms/%s/BED/bed_chr_%s.bed.gz",
				S_species, Sp_chr);
			cStrFile* X = D().query(S_bedPath);

			/* Skip the first line */
			X->gets(Sp_buf, 1024);

			char *a = NULL;
			char *b = NULL;
			for (wsUint L=1 ; X->gets(Sp_buf, 1024) ; L++) {
				getString(&Sp_buf, &a);
				/* second column is 0-based */
				getString(&a, &b);
				wsUint N_pos = atoi(a);
				getString(&b, &a);
				/* a = rsname */
				getString(&a, &b);
				/* Should be 'rs' */
				if (a[0]!='r' || a[1]!='s') continue;
				/* rsname */
				wsUint N_rsid = atoi(a+2);

				/* Make entry and insert */
				Xv_map[i].insert(make_pair(N_rsid, N_pos));
			}
			delete X;

			DEALLOC(Sp_buf);
		}

		/* Validate */
		FOREACH (vVariant_it, Xv_variant, i) {
			/* Skip if no chr */
			if (i->chr == 0 || i->chr > N_maxChr) continue;
			/* Skip if not starts with 'rs' */
			if (!i->name[0] || !i->name[1] || !i->name[2] ||
				(i->name[0] != 'r' && i->name[0] != 'R') ||
				(i->name[1] != 's' && i->name[1] != 'S')) continue;
			map<int,int>& X_curMap = Xv_map[i->chr];
			int N_rsn = atoi(i->name + 2);
			string S_vn = i->name;
			map<int,int>::iterator X_find = X_curMap.find(N_rsn);
			/* If RSID not exists */
			if (X_find == X_curMap.end()) {
				entryVariant(Cp_ex, *i);
				Cp_ex->fmt("	RSID_NOT_EXIST\n");
				N_inv++;
			}
			/* If position differs */
			if (X_find->second != (int)i->pos) {
				entryVariant(Cp_ex, *i);
				Cp_ex->fmt("	POSITION_ERR[%d]", X_find->second);
				N_inv++;
			}
		}
		if (N_inv)
			LOGwarn("[%d] variants marked as unmatch with reference information\n",
				N_inv);

		/* Remove */
		delete Cp_ex;
	}

#if TOOLSET_TYPE == TOOLSET_ONETOOL
	/* ONETOOL requires MAF computation */
	(void)getMAF();
#endif

	/* Make random phenotype */
	if (OPT_ENABLED(randbinpheno)) {
		LOG("Generate random binary phenotype...\n");

		if (OPT_ENABLED(1case)) for (wsUint i=0 ; i<N_sample ; i++)
			Ra_pheno[0][i] = (wsReal)(rand()%2);
		else for (wsUint i=0 ; i<N_sample ; i++)
			Ra_pheno[0][i] = (wsReal)(rand()%2)+1;
	} else if (IS_ASSIGNED(randpheno)) {
		LOG("Generate [%d] random continuous phenotype...\n", OPT_NUMBER(randpheno));

		N_pheno = OPT_NUMBER(randpheno);
		for (wsUint k=0 ; k<N_pheno ; k++) {
			for (wsUint i=0 ; i<N_sample ; i++) {
				wsReal R_p = (wsReal)((rand()%(RAND_MAX-2)+1) / (wsReal)(RAND_MAX));
				if (NA(R_p)) halt("Invalid pheno");
				Ra_pheno[k][i] = qnorm(R_p);
			}
			Ba_typePheno[k] = 1;
		}
	}

	/* Convert dichotomous phenotype into AFF/UNAFF */
	for (wsUint i=0 ; i<N_pheno ; i++) {
		if (isContinuous(i)) continue;
		for (wsUint j=0 ; j<N_sample ; j++) {
			if (Ra_pheno[i][j] == OPT_NUMBER(phenoCase))
				Ra_pheno[i][j] = WISARD_AFFECTED;
			else if (Ra_pheno[i][j] == OPT_NUMBER(phenoCtrl))
				Ra_pheno[i][j] = WISARD_UNAFFECTED;
		}
	}

	if (OPT_ENABLED(verbose)) {
		cExporter	*C			= cExporter::summon("stat.sample.lst");
		mSamp		&Xm_sampAll	= getSampleData();
		LOGoutput("Sample imported status is exported to [%s.stat.sample.lst]\n",
			OPT_STRING(out));
		wsUint i = 0;
		C->put("IDX	FILEIDX	FID	IID	PHENO\n");
		FOREACHDO (mSamp_it, Xm_sampAll, it, i++) {
			xSample &X_samp = it->second;
			if (it->second.N_idx == SAMP_NODATA)
				C->fmt("NA	%d	%s	%s	NA\n", X_samp.N_oriIdx,
					X_samp.S_FID.c_str(), X_samp.S_IID.c_str());
			else
				C->fmt("%d	%d	%s	%s	%g\n", X_samp.N_idx, X_samp.N_oriIdx,
					X_samp.S_FID.c_str(), X_samp.S_IID.c_str(),
					Ra_pheno[0][it->second.N_idx]);
		}
		delete C;
	}

	/* Make random phenotype */
// 	if (OPT_ENABLED(randbinpheno)) {
// 		LOG("Generate random binary phenotype...\n");
// 
// 		sseUnmat(Ra_pheno, N_pheno);
// 		N_pheno = 1;
// 		Ra_pheno = sseMatrix(1, N_sample);
// 		for (wsUint i=0 ; i<N_sample ; i++)
// 			Ra_pheno[0][i] = (wsReal)(rand()%2);
// 	} else if (OPT_ENABLED(randpheno)) {
// 		LOG("Generate random continuous phenotype...\n");
// 
// 		sseUnmat(Ra_pheno, N_pheno);
// 		N_pheno = 1;
// 		Ra_pheno = sseMatrix(1, N_sample);
// 		for (wsUint i=0 ; i<N_sample ; i++) {
// 			wsReal R_p = rand() / RAND_MAX;
// 			Ra_pheno[0][i] = qnorm(R_p);
// 		}
// 	}

	/* Check phenotype&sex missingness */
	char	*Ba_misPheno	= NULL;
	char	*Ba_misSex		= NULL;
	char	B_pheMissExist	= 0;
	char	B_sexMissExist	= 0;
	wsUint	*Na_case		= NULL;
	wsUint	*Na_ctrl		= NULL;
	wsUint	*Na_miss		= NULL;
	sseCalloc(Ba_misPheno, char, N_sample);
	sseCalloc(Ba_misSex, char, N_sample);
	sseCalloc(Na_case, wsUint, N_pheno);
	sseCalloc(Na_ctrl, wsUint, N_pheno);
	sseCalloc(Na_miss, wsUint, N_pheno);
	for (wsUint i=0 ; i<N_pheno ; i++) {
		/* If dichotomous, force to 0(control) and 1(case) */
		if (isContinuous(i)) {
			for (wsUint j=0 ; j<N_sample ; j++) {
				if (isMissingReal(Ra_pheno[i][j])) {
					Ba_misPheno[j]	= 1;
					B_pheMissExist	= 1;
					Na_miss[i]++;
				}
				if (isMissing(Xa_sampleV2[j]->N_sex)) {
					Ba_misSex[j]	= 1;
					B_sexMissExist	= 1;
				}
			}
			if (Na_miss[i])
				LOG("Continuous phenotype [%s] have [%d] missings\n",
					V_phenos[i].S_name.c_str(), Na_miss[i]);
		} else {
			for (wsUint j=0 ; j<N_sample ; j++) {
				if (isMissing(Xa_sampleV2[j]->N_sex)) {
					Ba_misSex[j]	= 1;
					B_sexMissExist	= 1;
				}
				if (isMissingReal(Ra_pheno[i][j])) {
					Ba_misPheno[j] = 1;
					B_pheMissExist = 1;
					Na_miss[i]++;
				} else {
					/* Force it */
					if (Ra_pheno[i][j] == WISARD_AFFECTED) {
						Ra_pheno[i][j] = W1;
						Na_case[i]++;
					} else if (Ra_pheno[i][j] == WISARD_UNAFFECTED) {
						Ra_pheno[i][j] = W0;
						Na_ctrl[i]++;
					} else
						halt("Invalid phenotype value [%g] found", Ra_pheno[i][j]);
				}
			}
			LOG("Binary phenotype [%s] have [%d] cases, [%d] controls, and [%d] missings\n",
				V_phenos[i].S_name.c_str(), Na_case[i], Na_ctrl[i], Na_miss[i]);
		}
	}
	/* Stdize phenotype if enabled, only dividing by its s.d. */
	if (OPT_ENABLED(phenostdize)) {
		if (IS_ASSIGNED(est)) {
			int			N_est	= EMAI_EST_SIZE + ((N_cov+1)<<1);
			wsMat		Ra_est = sseMatrix(N_pheno, N_est<<1);

			/* If --est, use sig2 in here instead */
			LOG("--est found, --phenostdize performs standardization using this file\n");
			loadEstimates(OPT_STRING(est), Ra_est, N_pheno, N_cov+1);

			for (wsUint i=0 ; i<N_pheno ; i++) {
				wsReal*	Ra_ph	= Ra_pheno[i];

				/* Get sd */
				wsReal R_sd = sqrt(Ra_est[i][0]);

				LOG("Phenotype [%s] standardized by s.d.[%g] from --est\n",
					V_phenos[i].S_name.c_str(), R_sd);

				/* Divide them */
				for (wsUint j=0 ; j<N_sample ; j++)
					if (!isMissingReal(Ra_ph[j]))
						Ra_ph[j] /= R_sd;
			}

			sseUnmat(Ra_est, N_pheno);
		} else {
			for (wsUint i=0 ; i<N_pheno ; i++) {
				wsReal	R_m		= W0;
				wsReal	R_sd	= W0;
				wsReal*	Ra_ph	= Ra_pheno[i];
				wsUint	N		= 0;

				/* Get mean */
				for (wsUint j=0 ; j<N_sample ; j++)
					if (!isMissingReal(Ra_ph[j])) {
						N++;
						R_m += Ra_ph[j];
					}
				R_m /= (wsReal)N;

				/* Get sd */
				for (wsUint j=0 ; j<N_sample ; j++)
					if (!isMissingReal(Ra_ph[j])) {
						wsReal R_v = Ra_ph[j] - R_m;
						R_sd += SQR(R_v);
					}
				R_sd /= (wsReal)(N-1);
				R_sd = sqrt(R_sd);

				LOG("Phenotype [%s] standardized by s.d.[%g] from [%d] samples\n",
					V_phenos[i].S_name.c_str(), R_sd, N);

				/* Divide them */
				for (wsUint j=0 ; j<N_sample ; j++)
					if (!isMissingReal(Ra_ph[j]))
						Ra_ph[j] /= R_sd;
			}
		}
	}

	/* Remove all missing founders if required */
// 	if (OPT_ENABLED(filmf)) {
// 	}

	/* Export phenotype missingness information */
	if (B_pheMissExist) {
/**/	cExporter*	Cp_misPheno = cExporter::summon("miss.pheno.lst");
		wsUint		N_miss = 0;
		Cp_misPheno->put("FID	IID\n");
		for (wsUint i=0 ; i<N_sample ; i++) {
			if (Ba_misPheno[i]) {
				Cp_misPheno->fmt("%s	%s\n", Xa_sampleV2[i]->S_FID.c_str(),
					Xa_sampleV2[i]->S_IID.c_str());
				N_miss++;
			}
		}
		delete Cp_misPheno;
		LOG("Found that [%d] samples having at least one missing phenotype\n", N_miss);
		LOGoutput("List of them exported to [%s.miss.pheno.lst]\n", OPT_STRING(out));
	}
	/* Export sex missingness information */
	if (B_sexMissExist) {
/**/	cExporter*	Cp_misSex = cExporter::summon("miss.sex.lst");
		wsUint		N_miss = 0;
		Cp_misSex->put("FID	IID\n");
		for (wsUint i=0 ; i<N_sample ; i++) {
			if (Ba_misSex[i]) {
				Cp_misSex->fmt("%s	%s\n", Xa_sampleV2[i]->S_FID.c_str(),
					Xa_sampleV2[i]->S_IID.c_str());
				N_miss++;
			}
		}
		delete Cp_misSex;
		LOG("Found that [%d] samples having sex missing\n", N_miss);
		LOGoutput("List of them exported to [%s.miss.sex.lst]\n", OPT_STRING(out));
	}
	sseFree(Ba_misPheno);
	sseFree(Ba_misSex);
	sseFree(Na_case);
	sseFree(Na_ctrl);
	sseFree(Na_miss);

	/* Process --blup */
	if (IS_ASSIGNED(blup))
		_loadBlup();

	/* Process --weight */
	if (IS_ASSIGNED(weight))
		_loadVariantWeight();

	/* Inbreeding coefficient */
	if (OPT_ENABLED(inbreed)) {
		cExporter *Cp_ib = cExporter::summon("inbreed.coef.res");
		LOGoutput("Sample-wise inbreeding coefficient F is exported to "
			"[%s.inbreed.coef.res]\n", OPT_STRING(out));
		/* Print header */
		Cp_ib->put("FID	IID	OBSHOM	EXPHOM	MAC	F\n");
		for (wsUint i=0 ; i<N_sample ; i++) {
			xSample *Xp_s = Xa_sampleV2[i];
			wsReal O, E, N;
			wsReal F = _inbreed(i, getVariant(), getMAF(), getGenotype(), O, E, N);

			/* NA */
			if (NA(F)) {
				if (!OPT_ENABLED(remna))
					Cp_ib->fmt("%s	%s	%g	%g	%g	NA\n", Xp_s->S_FID.c_str(),
						Xp_s->S_IID.c_str(), O, E, N);
			} 	else
				Cp_ib->fmt("%s	%s	%g	%g	%g	%g\n", Xp_s->S_FID.c_str(),
					Xp_s->S_IID.c_str(), O, E, N, F);
		}
		delete Cp_ib;
	}

	/* --variant2cov */
	/* Add to covariates if required */
	if (IS_ASSIGNED(variant2cov)) {
		vCovar&		Xa_cov	= getCovInfo();
		vVariant&	Xv_vrt	= getVariant();
		wsMat		Ra_cov	= getCovariates();
		wsUint		N_cov	= sizeCovar();
		wsUint		N_samp	= sizeSample();

		/* Try to open */
		wsUint	N_pc		= 0;
		char**	Sa_variants	= NULL;
		wsUint*	Na_indices	= NULL;
		cStrFile C_m2(OPT_STRING(variant2cov), "Covariate variants list", 1);
		if (C_m2.isFailed())
			Sa_variants = loadStringValues(OPT_STRING(variant2cov), &N_pc);
		else {
			char *S_buf = NULL, *a = NULL;
			wsAlloc(S_buf, char, 512);
			for ( ; C_m2.gets(S_buf, 512) ; N_pc++) {
				getString(&S_buf, &a);
				if (!S_buf[0]) break;
			}
			C_m2.rewind();

			wsAlloc(Sa_variants, char*, N_pc);
			for (wsUint i=0 ; C_m2.gets(S_buf, 512) ; i++) {
				getString(&S_buf, &a);
				Sa_variants[i] = strdup(S_buf);
				if (!S_buf[0]) break;
			}
			C_m2.close();

			DEALLOC(S_buf);
		}
		wsAlloc(Na_indices, wsUint, N_pc);

		/* Build variants map */
		mDataIdx	Xm_mmap[MAX_NCHR+1];
		wsUint		j = 0;
		FOREACHDO (vVariant_it, Xv_vrt, i, j++) {
			if (i->chr > 0 && (wsUint)i->chr <= NCHR_SPECIES)
				Xm_mmap[i->chr].insert(make_pair((string)i->name, j));
			else
				Xm_mmap[0].insert(make_pair((string)i->name, j));
		}

		/* Get their index */
		for (wsUint i=0 ; i<N_pc ; i++) {
			string	S_key	= (string)Sa_variants[i];
			wsUint	N_idx	= 0;
			for (j=0 ; j<=NCHR_SPECIES ; j++) {
				mDataIdx_it X_find = Xm_mmap[j].find(S_key);
				if (X_find != Xm_mmap[j].end()) {
					N_idx = X_find->second;
					break;
				}
			}
			/* Failed to found */
			if (j > NCHR_SPECIES)
				halt("Variant [%s] in --variant2cov is not exists in the dataset",
					Sa_variants[i]);
			Na_indices[i] = N_idx;
		}


		wsUint	N_nCov	= N_cov + N_pc;
		wsMat	Ra_nCov	= NULL;
		wsAlloc(Ra_nCov, wsReal*, N_nCov);

		/* Copy original */
		for (wsUint i=0 ; i<N_cov ; i++)
			Ra_nCov[i] = Ra_cov[i];
		DEALLOC(Ra_cov); /* Deallocate original */

		/* Copy newly computed PCs into covariates */
		char **Na_gn = getGenotype();
		for (wsUint i=0 ; i<N_pc ; i++) {
			/* Add to covar */
			xCovar X_new;
			X_new.N_idx			= N_cov+i;
			X_new.N_idxFile		= -1; /* Not originated from file */
			X_new.N_szFactor	= 0; /* Set to 0 since it is NOT factor */
			X_new.Sp_bl			= NULL;
			char S_pcName[64];
			sprintf(S_pcName, "VARIANT_%s", Sa_variants[i]);
			X_new.Sp_varName	= strdup(S_pcName);
			X_new.X_type		= WISARD_VAR_REAL;
			Xa_cov.push_back(X_new);

			/* Copy to covar */
			wsUint N_idx = Na_indices[i];
			sseMalloc(Ra_nCov[N_cov+i], wsReal, N_samp);
			for (wsUint j=0; j<N_samp ; j++)
				Ra_nCov[N_cov+i][j] = Na_gn[j][N_idx];
		}

		/* Free memory */
		for (wsUint i=0 ; i<N_pc ; i++)
			free(Sa_variants[i]);
		DEALLOC(Sa_variants);
		DEALLOC(Na_indices);

		chgCovInfo(N_nCov, Ra_nCov);
	}


	LOG("%d family read\n", Xm_fid2data.size());
	int N_indCount = 0;
	FOREACH (mFam_it, Xm_fid2data, Xp_it)
		//		printf("Family %s have %u members\n", (*Xp_it).first.c_str(), (MYuint)X_famData[N_famIdx].size());
		N_indCount += (*Xp_it).second.N_size;
	LOG("%d[?] individuals read\n", N_indCount);

	/* Perform inverse normalization */
	if (OPT_ENABLED(invnorm)) {
		/* All phenotypes must be continuous */
		for (wsUint i=0 ; i<N_pheno ; i++)
			if (!isContinuous(i)) halt("[%d]th phenotype is dichotomous, "
				"--invnorm cannot be applied!", i+1);
		/* Now apply */
		for (wsUint i=0 ; i<N_pheno ; i++)
			_inverseNormalize(getPhenos()[i], N_sample);
	}

	/* Post-processing */
	wsUintCst	N_missGeno	= getMissGenoCount();
	wsUintCst	N_vrt		= sizeVariant();
	wsUintCst	N_sample	= sizeSample();
	INT64_t	N_cntGeno	= (INT64_t)N_sample*(INT64_t)N_vrt;

	if (N_vrt && N_sample)
		LOG("Total genotyping rate : %.2f%% [" FMT_INT64 "/" FMT_INT64 "]\n",
			REAL_CONST(100.0)
			*(W1-(N_missGeno/((wsReal)N_sample*(wsReal)N_vrt))),
			N_cntGeno-N_missGeno,
			N_cntGeno);
	else
		LOG("No genotype detected\n");

	/* Deallocate space */
	Xe_varSubset.clear();

	/* Integrity check */
	if (V_sampWgt.size() && N_sample != V_sampWgt.size())
		halt("Sample weight #[%d] is not match with final sample size[%d]",
			V_sampWgt.size(), N_sample);

	/* Given dataset is complete(Have no missing) */
	if (N_missGeno == 0) {
		LOG("This dataset have full genotypes, so analysis speed will be boosted!\n");
		B_isDataComplete = 1;
	}
}

int cIO::getDiploid(xVariant &X_snv, char *Sp_buf, char **Sp_next)
{
	char	B_consec = OPT_ENABLED(consecallele);
	char	S_delim	= IS_ASSIGNED(sepallele) ? OPT_STRING(sepallele)[0] : 0xff;
	char*	Np_ACGT	= !Xe_varSubset.size() ||
		(Xe_varSubset.size() && Xe_varSubset.find(X_snv.name) != Xe_varSubset.end()) ?
		Na_ACGT : NULL;
	char*	Sp_1 = NULL, *Sp_2 = NULL;
	if (!Sp_buf) return -1;		/* FATAL - Null buffer */

	/* Get two genotypes */
	if (B_consec) {
		getString(&Sp_buf, Sp_next);
		Sp_1 = Sp_buf;
		Sp_2 = Sp_buf + 1;
		if (Sp_buf[2]) return -2; /* FATAL - NOSINGLEGENOTYPE */
	} else {
		getStringDelim(&Sp_buf, &Sp_2, S_delim);
		Sp_1 = Sp_buf;
		if (!Sp_2) return -1;		/* FATAL - Null buffer */
		getString(&Sp_2, Sp_next);
	} 

	/* If w/o --indel, check they are over than single byte */
	/* Check genotype missingness */
	char B_miss1, B_miss2;
	if (!OPT_ENABLED(indel)) {
		if ((!B_consec && Sp_1[1]) || Sp_2[1]) return -2; /* FATAL - NOSINGLEGENOTYPE */

		B_miss1 = Sp_1[0] == OPT_STRING(misgeno)[0];
		B_miss2 = Sp_2[0] == OPT_STRING(misgeno)[0];
	} else {
		B_miss1 = !stricmp(Sp_1, OPT_STRING(misgeno));
		B_miss2 = !stricmp(Sp_2, OPT_STRING(misgeno));
	}

	/* Check single-missing */
	if (B_miss1^B_miss2)
		return -3; /* FATAL - SINGLEMISSING */
	else if (B_miss1&B_miss2)
		return WISARD_NA;

	/* Check genotype */
	char C_geno1 = OPT_ENABLED(indel)?_getIndel(X_snv, Sp_1, Np_ACGT):
		_getVariant(X_snv, Sp_1[0], Np_ACGT);
	char C_geno2 = OPT_ENABLED(indel)?_getIndel(X_snv, Sp_2, Np_ACGT):
		_getVariant(X_snv, Sp_2[0], Np_ACGT);

	/* Check trialleleic */
	if (C_geno1 == -1 || C_geno2 == -1) {
		X_snv.filter = 1;
		return -4;
	}

	/* Sum of genotype is FINAL */
	return C_geno1 + C_geno2;
}

typedef struct _xFamMembers
{
	const char	*Ba_filt;
	char		B_filtExcl;
	vInt		Xa_members;
} xFamMembers;

int getFamMember(xSample *Xp_s, void *Vp_data)
{
	xFamMembers *Xp_fm = (xFamMembers *)Vp_data;
	if (Xp_s->B_isMissing)
		return 0;

	/* Return if this sample should be filtered */
	if (Xp_fm->Ba_filt)
		if ((Xp_fm->Ba_filt[Xp_s->N_idx] && Xp_fm->B_filtExcl) ||
			(!(Xp_fm->Ba_filt[Xp_s->N_idx]) && !(Xp_fm->B_filtExcl)))
			return 0;

	/* Insert to vector */
	Xp_fm->Xa_members.push_back(Xp_s->N_idx);
	return 0;
}

typedef struct _xFamMembers2
{
	int		*Na_posSamp;
	vInt	Xa_members;
} xFamMembers2;

int getFamMember2(xSample *Xp_s, void *Vp_data)
{
	xFamMembers2 *Xp_fm = (xFamMembers2 *)Vp_data;
	if (Xp_s->B_isMissing)
		return 0;
	/* Return if this sample should be filtered */
	if (Xp_fm->Na_posSamp[Xp_s->N_idx] == -1)
		return 0;

	/* Insert to vector */
	Xp_fm->Xa_members.push_back(Xp_fm->Na_posSamp[Xp_s->N_idx]);
	return 0;
}

wsReal** cIO::getFamInv(const char *Ba_filt, wsReal **Ra_corr,
	char B_filtExcl/*=1*/)
{
	wsUint*		Na_updIdx	= NULL;
	wsUint*		Na_szFam	= NULL;
	mFam&		Xa_fam		= getFamilyData();
	vSampPtr&	Xa_samp		= getSample();
	wsCalloc(Na_szFam, wsUint, Xa_fam.size());
	wsUint*		Na_fromIdx	= NULL;
	wsUint*		Na_toIdx	= NULL;
	wsAlloc(Na_fromIdx, wsUint, N_sample);
	wsAlloc(Na_toIdx, wsUint, N_sample);

	wsUint N_dif = 0;
	wsAlloc(Na_updIdx, wsUint, N_sample);
	LOOP (i, N_sample) {
		if (Ba_filt && B_filtExcl == Ba_filt[i]) N_dif++;
		Na_updIdx[i] = N_dif;
	}

	wsUint i=0, I=0;
	wsUint N_sz = 0;
	FOREACHDO (mFam_it, Xa_fam, fit, i++) {
		xFamily &X_fam = fit->second;
		for (vSampPtr_it it=Xa_samp.begin() ; it !=Xa_samp.end() ; it++) {
			if ((*it)->N_idx == SAMP_NODATA) continue;

			if ((*it)->S_FID.compare(X_fam.S_FID) == 0 && (!Ba_filt ||
				Ba_filt[(*it)->N_idx] != B_filtExcl)) {
				Na_fromIdx[I] = (*it)->N_idx;
				Na_toIdx[I++] = (*it)->N_idx - Na_updIdx[(*it)->N_idx];
				Na_szFam[i]++;
				N_sz++;
			}
		}
		/* Rebuild matrix */
	}
	wsMat	Ra_ret		= sseEmptyMatrix(N_sz, N_sz);
	i = 0;
	wsUint X = 0;
	FOREACHDO (mFam_it, Xa_fam, fit, i++) {
		wsUint		N_fsz	= Na_szFam[i];
		if (N_fsz == 0) continue;
		cSymMatrix	M_sub(N_fsz);
		wsSym		Ra_sub	= M_sub.get();
		wsUint		N_sIdx	= X;

		/* Rebuild corr matrix */
		for (wsUint j=0 ; j<N_fsz ; j++,X++) {
			wsUint V = Na_fromIdx[X];
			for (wsUint k=0,Z=N_sIdx ; k<j ; k++,Z++) {
				wsUint W = Na_fromIdx[Z];
				Ra_sub[j][k] = V>W?Ra_corr[V][W]:Ra_corr[W][V];
			}
			Ra_sub[j][j] = Ra_corr[V][V];
		}
		/* Make inversion */
		cSymMatrix&	M_sInv	= M_sub.inv();
		wsSym		Ra_sInv	= M_sInv.get();

		/* Remap */
		X = N_sIdx;
		for (wsUint j=0 ; j<N_fsz ; j++,X++) {
			wsUint V = Na_toIdx[X];
			for (wsUint k=0,Z=N_sIdx ; k<j ; k++,Z++) {
				wsUint W = Na_toIdx[Z];
				Ra_ret[V][W] = Ra_ret[W][V] = Ra_sInv[j][k];
			}
			Ra_ret[V][V] = Ra_sInv[j][j];
		}

		delete &M_sInv;
	}
	DEALLOC(Na_fromIdx);
	DEALLOC(Na_szFam);

	/* Now we got the inversion of entire data */
	return Ra_ret;
}

const char*	cIO::getAvailPheno(wsUint *Np_availSamp, wsReal **Rp_availY/*=NULL*/)
{
	wsVecCst Ra_1stpheno = getPheno();
	if (Np_availSamp == NULL)
		halt("SYSERR: Pointer to store # of available sample is empty");
	/* FIXME : Does not compatible with multiple phenotypes */

	wsUint	i, j;
	/* 1 if nonmissing, 0 if missing */
	char	*Ba_ret		= NULL;
	wsCalloc(Ba_ret, char, sizeSample());

	/* Count */
	for (i=j=0 ; i<N_sample ; i++)
		if (!isMissingReal(Ra_1stpheno[i])) {
			Ba_ret[i] = 1;
			j++;
		}
	*Np_availSamp = j;

	/* Allocate and fill phenotypes of phenotype-available samples if required */
	if (Np_availSamp) {
		sseMalloc(*Rp_availY, wsReal, j);
		for (i=j=0 ; i<N_sample ; i++)
			if (!isMissingReal(Ra_1stpheno[i]))
				(*Rp_availY)[j++] = Ra_1stpheno[i];

		/* Sanity check */
		if (j != (*Np_availSamp))
			halt("SYSERR: Available phenotype buffer is not properly constructed");
	}

	return Ba_ret;
}

} // End namespace ONETOOL

#include "input/plink.h"
#include "input/vcf.h"
#include "input/dosage.h"
#include "input/sim.h"
#include "input/meta.h"

namespace ONETOOL {
/*
#if TOOLSET_TYPE == TOOLSET_ONETOOL
cIO* cIO::summon(xFileType X_type, char *Sp_fn, app_log& alog, char B_dry)
{
//	cIO	**Cp_IOs = NULL;

	// Get files
//	vStr X_ret = bloating(Sp_fn);
	
//	if (X_ret.size() == 1) {
		switch (X_type) {
		case FT_SIM:
			return new cSimFamIO(Sp_fn, B_dry);
		case FT_PED: case FT_TPED: case FT_BED: case FT_LGEN: case FT_BEAGLE:
			return new cPlinkIO(Sp_fn, X_type, B_dry);
		case FT_VCF:
			return new cVcfIO(Sp_fn, alog, B_dry);
		case FT_DOSAGE:
		case FT_EXPR:
			return new cDosageIO(Sp_fn, B_dry);
		case FT_META:
			return new cMetaIO(Sp_fn, B_dry);
		default:
			halt("SYSERR : Invalid format requested to make I/O object");
		}
//	}

// 	LOG("Below files will be loaded...\n");
// 	FOREACH (vStr_it, X_ret, i)
// 		LOG("	%s\n", i->c_str());

	//WORKER().run(thrLoad, dist_fileLoad, this, Cp_IOs);
	// Now aggregate them
	return NULL;
}
#else
*/
cIO* cIO::summon(xFileType X_type, char *Sp_fn, char B_dry/*=0*/)
{
//	cIO	**Cp_IOs = NULL;

	/* Get files */
//	vStr X_ret = bloating(Sp_fn);
	
//	if (X_ret.size() == 1) {
		switch (X_type) {
		case FT_SIM:
			return new cSimFamIO(Sp_fn, B_dry);
		case FT_PED: case FT_TPED: case FT_BED: case FT_LGEN: case FT_BEAGLE:
			return new cPlinkIO(Sp_fn, X_type, B_dry);
		case FT_VCF:
			return new cVcfIO(Sp_fn, X_type, B_dry);
		case FT_DOSAGE:
		case FT_EXPR:
			/* Check Sp_fn is VCF / VCF.GZ */ {
				char* Sp_end = Sp_fn + strlen(Sp_fn);
				if (!stricmp(Sp_end - 3, "vcf") || !stricmp(Sp_end - 6, "vcf.gz"))
					return new cVcfIO(Sp_fn, X_type, B_dry);
				else
					return new cDosageIO(Sp_fn, B_dry);
			}
		case FT_META:
			return new cMetaIO(Sp_fn, B_dry);
		default:
			halt("SYSERR : Invalid format requested to make I/O object");
		}
//	}

// 	LOG("Below files will be loaded...\n");
// 	FOREACH (vStr_it, X_ret, i)
// 		LOG("	%s\n", i->c_str());

	//WORKER().run(thrLoad, dist_fileLoad, this, Cp_IOs);
	/* Now aggregate them */
	return NULL;
}
//#endif

void cIO::clear()
{
	DEALLOC(Ba_isFounder);
	if (Ra_pheno) sseUnmat(Ra_pheno, N_pheno);
	if (Ra_data) sseUnmat(Ra_data, N_sample);
	if (Ra_dosage) sseUnmat(Ra_dosage, N_sample);
	if (Na_geno) {
		for (wsUint i=0 ; i<N_sample ; i++) {
			sseFree(Na_geno[i]);
		}
		DEALLOC(Na_geno);
	}
	DEALLOC(Na_charData);
	DEALLOC(Ba_typePheno);
	if (Ra_cov) sseUnmat(Ra_cov, N_cov);
	sseFree(Ra_customWeight);
	/* Remove string of xVariant */
	FOREACH (vVariant_it, Xv_variant, it) {
		//LOG(it->name);
		DEALLOC(it->name);
		DEALLOC(it->anno);
	}
	//sseFree(Ra_corrMask);
	if (Xv_altVariantDef)
		delete Xv_altVariantDef;

	/* cIO clear */
	Xa_sampleV2.clear();
	/* xVariant clear */
	FOREACH (vVariant_it, Xv_variant, i) {
		if (OPT_ENABLED(indel)) {
			free(i->indel1);
			free(i->indel2);
		}
		DEALLOC(i->name);
	}
}

char cIO::_filterVariantPosition(xVariant &X_snv)
{
	char B_filter = 0;

	/* Range-based filter */
	if (N_filterRange != FILT_NOTHING) {
		char B_inRange = 0;
		FOREACH (vector<xRange>::iterator, Xv_listRanges[X_snv.chr], i) {
			/* Check if IN RANGE */
			if (i->N_start <= X_snv.pos && X_snv.pos <= i->N_end) {
				B_inRange = 1;
				break;
			}
		}

		if (N_filterRange == FILT_SELECT)
			B_filter = !B_inRange;
		else if (N_filterRange == FILT_REMOVE)
			B_filter = B_inRange;
		else
			halt("SYSERR : Invalid range filtering option [%d]", N_filterRange);
	}

	if (IS_ASSIGNED(incgdist) && !isInRange(OPT_RANGE(incgdist), X_snv.gdist))
		B_filter = 1;
	if (IS_ASSIGNED(filgdist) && isInRange(OPT_RANGE(filgdist), X_snv.gdist))
		B_filter = 1;

#ifdef _DEBUG
// 	if (Xm_listVrt.size()) {
// 		pverbose("Following variants will be selected/filtered\n");
// 		FOREACH (mStrKey_it, Xm_listVrt, i)
// 			pverbose("	[%s]\n", i->first.c_str());
// 	}
#endif

	/* Name-based filter */
	if (N_filterVrt) {
		if (!OPT_ENABLED(regex)) {
			eStr_it X_find = Xe_listVrt.find(X_snv.name);
			if (X_find != Xe_listVrt.end()) {
				if ((N_filterVrt == FILT_REMOVE && X_find != Xe_listVrt.end()) ||
					(N_filterVrt == FILT_SELECT && X_find == Xe_listVrt.end())) /* Remove */
					B_filter = 1;
			} else if (N_filterVrt == FILT_SELECT)
				/* Filter if not listed */
				B_filter = 1;
		} else FOREACH (vPattern_it, Xv_listVrt, i) {
			/* --regex */
			bool B_res = trex_match(*i, X_snv.name);
			if ((B_res && N_vrtOrig != FILT_REMOVE) ||
				(!B_res && N_vrtOrig != FILT_SELECT))
				continue;
			/* Filter if found&remove || notfound&select */
			B_filter = 1;
		}
	}

	/* 150216 fixed to not filter 0
	 * Chromosome-based filter */
	if ((X_snv.chr > 0 && Ba_chrAllowed[X_snv.chr - 1] == 0) ||
		(X_snv.chr == 0 && IS_ASSIGNED(chr)) ||
		(OPT_ENABLED(autoonly) && getChrAuto(X_snv.chr) == -1) ||
		(OPT_ENABLED(sexonly) && getChrAuto(X_snv.chr) != -1))
		B_filter = 1;

	/* Marking it exists if not filtered */
	if (B_filter == 0) {
		if (X_snv.chr != -1) Ba_chrExists[X_snv.chr-1] = 1;
		else Ba_chrExists[MAX_NCHR] = 1;
	}

	/* Return filtering result */
	return B_filter;
}

void cIO::exportSampleMap()
{
	cExporter *Cp_e = cExporter::summon("all.samp.lst");
	LOGoutput("Detailed status of imported samples is exported to [%s.all.samp.lst]\n",
		OPT_STRING(out));
	Cp_e->put("FID	IID	MISSING	ORIIDX	IDX	SEX	MISGENO\n");

	FOREACH (mSamp_it, Xm_iid2data, i) {
		xSample &X = i->second;
		Cp_e->fmt("%s	%s	%d	%d	%d	%d	%d\n", X.S_FID.c_str(),
			X.S_IID.c_str(), X.B_isMissing, X.N_oriIdx, X.N_idx, X.N_sex,
			X.N_missGeno);
	}
	delete Cp_e;
}

wsReal*	cIO::subsetPheCov(char *Ba_filt, char B_incVal, wsUint N_sz,
	wsMat *Rp_cov)
{
	wsRealCst*	Ra_phe	= getPheno();
	wsUint	N_samp	= sizeSample();
	wsReal*	Ra_ret = sseVector(N_sz);
	*Rp_cov = sseMatrix(sizeCovar(), N_sz);

	if (B_incVal == 0) for (wsUint i=0,j=0 ; i<N_samp ; i++) {
		if (!Ba_filt[i]) continue;
		Ra_ret[j] = Ra_phe[i];
		for (wsUint k=0 ; k<N_cov ; k++)
			(*Rp_cov)[k][j] = Ra_cov[k][i];
	} else for (wsUint i=0,j=0 ; i<N_samp ; i++) {
		if (Ba_filt[i]) continue;
		Ra_ret[j] = Ra_phe[i];
		for (wsUint k=0 ; k<N_cov ; k++)
			(*Rp_cov)[k][j] = Ra_cov[k][i];
	}

	return Ra_ret;
}

wsReal*	cIO::subsetPheCovItct(char *Ba_filt, char B_incVal, wsUint N_sz,
	wsMat *Rp_cov, wsUint *Np_cov/*=NULL*/)
{
	wsRealCst*	Ra_phe	= getPheno();
	wsUint	N_samp	= sizeSample();
	wsReal*	Ra_ret = sseVector(N_sz);
	*Rp_cov = sseMatrix(sizeCovar()+1, N_sz);

	/* Fill the fist line of return covariate matrix with 1 */
	sseInit((*Rp_cov)[0], N_sz, 1);

	if (B_incVal == 0) for (wsUint i=0,j=0 ; i<N_samp ; i++) {
		if (!Ba_filt[i]) continue;
		Ra_ret[j] = Ra_phe[i];
		for (wsUint k=0 ; k<N_cov ; k++)
			(*Rp_cov)[k+1][j] = Ra_cov[k][i];
	} else for (wsUint i=0,j=0 ; i<N_samp ; i++) {
		if (Ba_filt[i]) continue;
		Ra_ret[j] = Ra_phe[i];
		for (wsUint k=0 ; k<N_cov ; k++)
			(*Rp_cov)[k+1][j] = Ra_cov[k][i];
	}

	if (Np_cov) *Np_cov = N_cov+1;
	return Ra_ret;
}

void cIO::impute(cAnalysis *Cp_anaPhi)
{
	/* FIXME : Do imputation */
	if (sizeFounder() == sizeSample()) {
		LOG("Impute with minor allele frequency\n");
		Cp_impute = new cMAFimpute(this);
	} else {
		LOG("Impute with family structure\n");
		Cp_impute = new cCorImpute(this, getFullCorMat(Cp_anaPhi));
	}
	Cp_impute->impute();
}

wsMat cIO::getImputeGeno()
{
	return Cp_impute ? Cp_impute->getData() : NULL;
}

void otherExt(wsStrCst S_fn, wsStrCst S_ext, char *Sp_dupFn, wsUint N_lvRem/*=1*/)
{
	strcpy(Sp_dupFn, S_fn);
	char *p = Sp_dupFn + strlen(Sp_dupFn) - 1;

	while (*p != '.' && p>=Sp_dupFn) p--;
	// 170911 gz process
	if (!stricmp(*p=='.'?p+1:p, "gz")) for (p-- ; *p!='.' && p>=Sp_dupFn ; p--);
	if (N_lvRem > 1) {
		for (wsUint i=1 ; p>=Sp_dupFn && i<N_lvRem ; i++) {
			p--;
			while (*p != '.') p--;
		}
	}
	size_t N_l = strlen(S_ext);
	p[0] = '.';
	memcpy(p+1, S_ext, N_l);
	p[N_l+1] = '\0';
}

void loadEstimates(wsStrCst Sp_strFn, wsMat Ra_estimates, wsUint N_pheno,
	wsUint N_cov)
{
	int			N_est		= EMAI_EST_SIZE + (N_cov<<1);
	int			N_idx		= 0;
	cStrFile	C_est(Sp_strFn);
	char		*S_buf		= NULL;
	char		*a			= NULL;
	char		*b			= NULL;
	wsAlloc(S_buf, char, 4096);

	/* Read first line */
	C_est.gets(S_buf, 4096);
	getString(&S_buf, &a);

	/* Is it "EST"? */
	if (stricmp(S_buf, "ESTNAME")) halt("Invalid poly.est, must start with 'ESTNAME'!");
	/* Get the number of columns */
	for (wsUint i=0 ; i<N_pheno ; i++,a=b) {
		if (a == NULL) halt("# of column[%d] less than # of phenotype[%d]",
			i, N_pheno);
		getString(&a, &b);
	}
	if (a) halt("# of column more than # of phenotype[%d]", N_pheno);
	C_est.gets(S_buf, 4096);

	/* Read */
	// 0      1       2      3           4
	// [sig2] [sig2g] [logL] [Var(sig2)] [Var(sig2g)]
	// 5     6          7
	// [h^2] [Var(h^2)] [Var(sig2+sig2g)]
	for (N_idx=0 ; C_est.gets(S_buf, 4096) ; N_idx++) {
		getString(&S_buf, &a);
		if (S_buf[0] && N_idx == N_est)
			halt("Too much number of estimates were found [%d expected], "
				"maybe the covariates are different? [%s]", N_est, S_buf);

		/* Get data */
		for (wsUint i=0 ; i<N_pheno ; i++,a=b) {
			if (a == NULL) halt("# of column[%d] less than # of phenotype[%d]",
				i, N_pheno);
			Ra_estimates[i][N_idx] = (wsReal)atof(a);
			getString(&a, &b);
		}
		if (a) halt("# of column more than # of phenotype[%d]", N_pheno);
	}
	if (N_idx != N_est)
		halt("Too less number of estimates were retrieved [%d expected, %d retrieved]",
			N_est, N_idx);

	DEALLOC(S_buf);
}

typedef map<string,wsMat>	mMatData;
typedef mMatData::iterator	mMatData_it;

inline void build(mMatData& Xm_mats, mvInt& Xm_idx, wsVec Ra_dest,
	wsVec Ra_src, wsUintCst N_sIdx, wsUintCst N_remIdx=0xffffffff)
{
	mvInt_it II = Xm_idx.begin();
	FOREACH (mMatData_it, Xm_mats, JJ) {
//		wsMat	Ra_tFam	= JJ->second;
		vInt&	Xv_tIdx	= (II++)->second;

		FOREACH (vInt_it, Xv_tIdx, J) {
			wsUint	N_tIdx = *J;
			/* Adjust index */
			if (N_tIdx == N_remIdx) continue;
			else if (N_tIdx > N_remIdx) N_tIdx--;
			/* Skip since it is sym. */
			if (N_tIdx < N_sIdx) continue;

			Ra_dest[N_tIdx] = Ra_src[*J];
		}
	}

}

wsSym buildMatrix(mMatData& Xm_mats, mvInt& Xm_idx, wsUintCst N_remIdx=0xffffffff,
	wsVec* Rp_12=NULL)
{
	wsUint N_sz = 0;
	FOREACH (mvInt_it, Xm_idx, i) N_sz += (wsUint)i->second.size();
	/* If N_remIdx, dim will decrease 1 */
	if (N_remIdx != 0xffffffff) N_sz--;

	wsSym	Ra_ret = sseEmptySymMat(N_sz);
//	wsUint	N_fam	= (wsUint)Xm_mats.size();
	if (N_remIdx == 0xffffffff) {
		mvInt_it II = Xm_idx.begin();
		FOREACH (mMatData_it, Xm_mats, JJ) {
			/* When it is just a generation */
//			wsMat	Ra_sFam	= JJ->second;
			vInt&	Xv_sIdx	= (II++)->second;

			FOREACH (vInt_it, Xv_sIdx, I) {
				wsUint	N_sIdx = *I;
				mvInt_it KK = Xm_idx.begin();
				FOREACH (mMatData_it, Xm_mats, LL) {
//					wsMat	Ra_tFam	= LL->second;
					vInt&	Xv_tIdx	= (KK++)->second;
					FOREACH (vInt_it, Xv_tIdx, J) {
						wsUint	N_tIdx = *J;
						if (N_tIdx < N_sIdx) continue;
						Ra_ret[N_sIdx][N_tIdx] = Ra_ret[*I][*J];
					}
				}
			}
		}
	} else {
		/* When there is a sample to exclude */
		mvInt_it II = Xm_idx.begin();
		FOREACH (mMatData_it, Xm_mats, JJ) {
//			wsMat	Ra_sFam	= JJ->second;
			vInt&	Xv_sIdx	= (II++)->second;

			FOREACH (vInt_it, Xv_sIdx, I) {
				wsUint	N_sIdx = *I;
				if (N_sIdx == N_remIdx) {
					*Rp_12 = sseVector(N_sz);
					/* Generate Ra_12 if desired */
					build(Xm_mats, Xm_idx, *Rp_12, Ra_ret[*I], N_sIdx, N_remIdx);
				} else if (N_sIdx > N_remIdx) N_sIdx--;

				/* Generation of Ra_22 */
				build(Xm_mats, Xm_idx, Ra_ret[N_sIdx], Ra_ret[*I],
					N_sIdx, N_remIdx);
			}
		}
	}
	return Ra_ret;
}

typedef struct _xThreadGetCondVar {
	cIO*		Cp_IO;
	vSampPtr&	Xa_samp;
	cSymMatrix&	M_cor;
	wsUint		N_samp;
} xThreadGetCondVar;

int thrGetCondVar(int N_idx, void *Vp_shareData, void *Vp_data,
	void *Vp_result)
{
	xThreadGetCondVar*	Xp_gc	= (xThreadGetCondVar *)Vp_shareData;
	int					i		= *((int *)Vp_data);
	xSample*			I		= Xp_gc->Xa_samp[i];
	wsSym				Ra_cor	= Xp_gc->M_cor.get();
	wsUint				N_samp	= Xp_gc->N_samp;
	wsVec				Ra_vs	= (wsVec)Vp_result;

	wsUint	j		= I->N_idx;
	wsReal	R_11	= Ra_cor[j][j];

	wsMat	Ra_22	= sseSubremMatrix(Ra_cor, N_samp, j);
	wsMat	Ra_22i	= invSymMat(Ra_22, N_samp-1);
	deallocMatrix(Ra_22, N_samp-1, (void *)1);
	wsVec	Ra_12	= sseRowMatrix(Ra_cor, N_samp, j, j);

	wsReal	R_var = R_11 - sseVpSpV(Ra_12, Ra_22i, N_samp - 1);
	deallocMatrix(Ra_22i, N_samp - 1, (void *)1);
	sseFree(Ra_12);

	Ra_vs[j] = R_var;

	return 0;
}

wsReal* getCondVar(cIO* Cp_IO, vSampPtr& Xa_samp, cSymMatrix& M_cor, wsMat *Rp_vt/*=NULL*/)
{
	if (Xa_samp.size() != M_cor.row())
		halt("Sample size[%d] and corr matrix[%d] is unmatch",
		Xa_samp.size(), M_cor.row());
	wsUint	N_samp = M_cor.row();
	wsReal*	Ra_vs = sseVector(N_samp);
	if (Rp_vt)
		*Rp_vt = sseMatrix(N_samp, N_samp); /* cond. var. except this sample */
	wsSym	Ra_cor = M_cor.get();
	wsUint k=0;

	mvInt		Xm_famIdx;
	mMatData	Xm_subFam;
	mMatData	Xm_subFamInv;
	/* if --kinship */ {
		wsUint I = 0;
		FOREACHDO (vSampPtr_it, Xa_samp, i, I++)
			Xm_famIdx[(*i)->S_FID].push_back(I);

		/* Make subset */
		FOREACH (mvInt_it, Xm_famIdx, i) {
			wsSym Ra_curFam = sseSsubset(Ra_cor, N_samp, i->second);
			Xm_subFam.insert(make_pair(i->first, Ra_curFam));
			/* Do inverse */
			cSymMatrix M_mat(Ra_curFam, (wsUintCst)i->second.size(), 1);
			cSymMatrix& M_inv = M_mat.inv();
			if (M_inv.row() == 0)
				halt("Failed to get an inverse of sample relatedness for FID [%s]", i->first.c_str());
			Xm_subFamInv.insert(make_pair(i->first, M_inv.get()));
			M_inv.setDontDealloc();
			delete& M_inv;
		}
	}

	if (Rp_vt) {
		FOREACHDO(vSampPtr_it, Xa_samp, i, k++) {
			wsUint	j = (*i)->N_idx;
			wsReal	R_11 = Ra_cor[j][j];

			wsMat	Ra_22 = sseSubremMatrix(Ra_cor, N_samp, j);
			wsMat	Ra_22i = NULL;
			Ra_22i = invSymMat(Ra_22, N_samp - 1);
			//deallocMatrix(Ra_22, N_samp-1, (void *)1);
			wsVec	Ra_12 = sseRowMatrix(Ra_cor, N_samp, j, j);

			wsSym	Ra_vx = sseVtV(Ra_12, N_samp - 1, W1 / R_11);
			sseSsS(Ra_22, Ra_vx, Ra_22, N_samp - 1);
			for (wsUint l = 0; l < k; l++)
				(*Rp_vt)[k][l] = Ra_22[l][l];
			(*Rp_vt)[k][k] = WISARD_NAN;
			for (wsUint l = k + 1; l < N_samp; l++)
				(*Rp_vt)[k][l] = Ra_22[l - 1][l - 1];
			deallocMatrix(Ra_22, N_samp - 1, (void *)1);

			wsReal	R_var = R_11 - sseVpSpV(Ra_12, Ra_22i, N_samp - 1);
			deallocMatrix(Ra_22i, N_samp - 1, (void *)1);
			sseFree(Ra_12);
			Ra_vs[j] = R_var;
			notice("%d samples processed...\r", k + 1);
		}
	} else {
		if (IS_MULTITHREAD) {
			xThreadGetCondVar X_gc = {
				Cp_IO, Xa_samp, M_cor, N_samp,
			};
			WORKER().run(thrGetCondVar, forAllSample_anaXthr, &X_gc, Ra_vs);
		} else FOREACHDO(vSampPtr_it, Xa_samp, i, k++) {
			wsUint	j		= (*i)->N_idx;
			wsReal	R_11	= Ra_cor[j][j];

			wsMat	Ra_22	= sseSubremMatrix(Ra_cor, N_samp, j);
			wsMat	Ra_22i	= invSymMat(Ra_22, N_samp-1);
			deallocMatrix(Ra_22, N_samp-1, (void *)1);
			wsVec	Ra_12	= sseRowMatrix(Ra_cor, N_samp, j, j);

			wsReal	R_var	= R_11 - sseVpSpV(Ra_12, Ra_22i, N_samp-1);
			deallocMatrix(Ra_22i, N_samp-1, (void *)1);
			sseFree(Ra_12);
			Ra_vs[j] = R_var;
			notice("%d samples processed...\r", k+1);
		}
	}

	return Ra_vs;
}

/* tr(phixy phiyy^-1 phiyx) */
wsReal* getCondVar2(vSampPtr& Xa_samp, cSymMatrix& M_cor)
{
	if (Xa_samp.size() != M_cor.row())
		halt("Sample size[%d] and corr matrix[%d] is unmatch",
		Xa_samp.size(), M_cor.row());
	wsUint	N_samp = M_cor.row();
	wsReal*	Ra_vs = sseVector(N_samp);
	wsSym	Ra_cor = M_cor.get();
	wsUint k=0;

	FOREACHDO (vSampPtr_it, Xa_samp, i, k++) {
		wsUint	j		= (*i)->N_idx;
//		wsReal	R_11	= Ra_cor[j][j];

		wsMat	Ra_22	= sseSubremMatrix(Ra_cor, N_samp, j);
		wsMat	Ra_22i	= invSymMat(Ra_22, N_samp-1);
		deallocMatrix(Ra_22, N_samp-1, (void *)1);
		wsReal	*Ra_12	= sseRowMatrix(Ra_cor, N_samp, j, j);


		wsReal	R_var	= sseVpSpV(Ra_12, Ra_22i, N_samp-1);
		deallocMatrix(Ra_22i, N_samp-1, (void *)1);
		sseFree(Ra_12);
		Ra_vs[j] = R_var;
		notice("%d samples processed...\r", k+1);
	}

	return Ra_vs;
}

/* Set weight */
wsVec cIO::getWeight(cPPPAnalysisV2 *Cp_anaPPP)
{
	// 150824 --noweight when --expression
	if (IS_ASSIGNED(expression)) {
		Ra_customWeight = NULL;
		return Ra_customWeight;
	}
	if (Ra_customWeight)
		return Ra_customWeight;

	/* Weight-related option check : one of --noweight, --mafweight, --betaweight, --weight */
	wsUint N_opt = OPT_ENABLED(noweight) + OPT_ENABLED(mafweight) +
		IS_ASSIGNED(betaweight) + IS_ASSIGNED(weight);
	if (N_opt > 1)
		halt("Only one weight-related option can be applied! [%s%s%s%sassigned]",
			OPT_ENABLED(noweight) ? "--noweight" : "",
			OPT_ENABLED(mafweight) ? "--mafweight" : "",
			IS_ASSIGNED(betaweight) ? "--betaweight" : "",
			IS_ASSIGNED(weight) ? "--weight" : "");

	/* If custom weight */
	if (IS_ASSIGNED(weight))
		_loadVariantWeight();
	else if (IS_ASSIGNED(mafweight)) {
		/* Apply MAF-weight */
		wsUint	N_vrt			= sizeVariant();
		wsVecCst	Ra_ppp_sqrt2pq	= Cp_anaPPP->getPPPsq();
		wsUint	N_med = 0;
		sseMalloc(Ra_customWeight, wsReal, N_variant);

#ifdef USE_SSE
		sse_t	sse_H	= sseSet(1.0);
		N_med	= getMed(N_vrt);
		for (wsUint i=0 ; i<N_med ; i+=sseJmp) {
			sse_t*	sse_pW	= (sse_t *)(Ra_customWeight + i);
			sse_t*	sse_pP	= (sse_t *)(Ra_ppp_sqrt2pq + i);
			*sse_pW = sseDiv(sse_H, *sse_pP);
		}
#endif
		for (wsUint i=N_med ; i<N_vrt ; i++)
			Ra_customWeight[i] = W1/Ra_ppp_sqrt2pq[i];
	} else if (OPT_ENABLED(noweight))
		Ra_customWeight = NULL;
	else {
		if (!IS_ASSIGNED(betaweight)) {
			OPTION().assign("betaweight", OPTION().getDefVal("betaweight"));
			OPTION().FORCE_OPT_STRING(betaweight);
		}
		/* Default = beta-weight */
		vVariant&	Xa_variant	= getVariant();
		wsUint i= 0;

		wsUint	N_beta	= 0;
		wsVec	Ra_bw	= loadRealValues(OPT_STRING(betaweight), &N_beta);
		sseMalloc(Ra_customWeight, wsReal, N_variant);
		if (N_beta != 2) halt("--betaweight should take 2 real values!");
		FOREACHDO (vVariant_it, Xa_variant, I, i++)
			Ra_customWeight[i] = (wsReal)dbeta(Xp_maf[i].R_maf, Ra_bw[0], Ra_bw[1], 0);
		sseFree(Ra_bw);
	}

	return Ra_customWeight;
}

} // End namespace ONETOOL
