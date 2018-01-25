#include <algorithm>
#include <numeric>
#include "analyses/setmgr.h"
#include "analyses/mdr.h"
#include "analyses/regr.h"
#include "output/result.h"
#include "utils/vector.h"
#include "utils/matrix.h"
#include "input/stream.h"
#include "utils/stat.h"

#define COMB_SIZE	2048
#define RES_SIZE	1

namespace ONETOOL {

typedef struct _xMdrCvcResult {
	string S_comb;
	wsReal R_mean, R_min, R_max;
} xMdrCvcResult;
typedef vector<xMdrCvcResult> vMdrCvcResult;
typedef vMdrCvcResult::iterator vMdrCvcResult_it;
typedef vMdrCvcResult::reverse_iterator vMdrCvcResult_rit;
typedef vector<vMdrCvcResult> vvMdrCvcResult;
typedef vvMdrCvcResult::iterator vvMdrCvcResult_it;
typedef vvMdrCvcResult::reverse_iterator vvMdrCvcResult_rit;

typedef struct _xThreadMdr {
	wsUint			N_curCV;
	vInt&			Na_testIdxCase;
	vInt&			Na_testIdxCtrl;
	cMdrAnalysis*	Cp_anaMdr;
} xThreadMdr;

typedef struct _xMdr {
	xMdrType	X_type;
	bool		B_haveTest;
	wsUint		N_order;			/* Order of MDR */
	wsUint		N_sample;			/* Samples included in MDR analysis */
	char**		Na_data;			/* Genotype data */
	wsVecCst		Ra_Y;				/* Phenotype data */
	wsUmat		Na_cntCase;			/* #sample for case in train(0), test(1) */
	wsUmat		Na_cntCtrl;			/* #sample for control in train(0), test(1) */
	wsUint*		Na_mapMarker;		/* Index mapping for marker */
	wsVecCst		Ra_sc;
	wsVec		Ra_scca;
	wsVec		Ra_scct;
	vInt&		Nv_testIdxCase;	/* Ref of indices of case samples in test */
	vInt&		Nv_testIdxCtrl;	/* Ref of indices of case samples in test */
	wsMat		Ra_permY;			/* Permuted Y's, dim is q * p, where q is # of perms */
} xMdr;

void _hmdr(xMdr* Xp_mdr, vHmdrEntry& Xv_he, vReal& Rv_minP, char** Na_data,
	wsUint N_vrt, cExporter *Cp_hmdr, vVariant& Xv_vrt);
void _hmdr2(xMdr* Xp_mdr, vHmdrEntry& Xv_he, vReal& Rv_minP, char** Na_data,
	wsUint N_vrt, cExporter *Cp_hmdr, vVariant& Xv_vrt);

wsUint nCk(wsUint n, wsUint k)
{
	if (k > n) return 0;
	if (k * 2 > n) k = n-k;
	if (k == 0) return 1;

	wsUint N_ret = n;
	for (wsUint i=2 ; i<=k ; i++) {
		N_ret *= (n-i+1);
		N_ret /= i;
	}

	return N_ret;
}

/* !!! HIMINI-specific analysis !!! */
#if ((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE & FUNCTION_HIMINI)

void CalculateMarginalEntropyNdist(char **Na_geno, wsUintCst N_mkr, wsUintCst N_samp,
	wsReal *MarginalEntropySNP, wsReal *MarginalEntropySNP_Y, wsReal *Ra_Y,
	wsUint ***Np_genoMD, wsUint ***Np_genoMDs)
{
	wsReal		ptmp;
	wsUint**	Na_genoMD	= NULL;
	wsUint**	Na_genoMDs	= NULL;
	wsAlloc(Na_genoMD, wsUint*, 6);
	wsAlloc(Na_genoMDs, wsUint*, 3);
	for (wsUint i=0 ; i<3 ; i++) {
		wsCalloc(Na_genoMDs[i], wsUint, N_mkr);
		wsCalloc(Na_genoMD[i], wsUint, N_mkr);
	}
	for (wsUint i=3 ; i<6 ; i++) {
		wsCalloc(Na_genoMD[i], wsUint, N_mkr);
	}
	//int GenoMarginalDistr[3][2];

	for (wsUint i=0 ; i<N_samp ; i++) {
		if (isMissingReal(Ra_Y[i])) continue;
		wsUint N_idx = (wsUint)(W1 - Ra_Y[i]) * 3;

		for (wsUint j=0 ; j<N_mkr ; j++) {
			char G = Na_geno[i][j];
			if (isMissing(G)) continue;
			Na_genoMD[N_idx + G][j]++;
		}
	}

	for (wsUint i=0 ; i<N_mkr ; i++) {
		for (wsUint j=0 ; j<3 ; j++) {
			wsReal tmp0 = (wsReal) Na_genoMD[j][i] + Na_genoMD[3+j][i];
			Na_genoMDs[j][i] = (wsUint)tmp0;
			if (tmp0 > 0) {
				ptmp = tmp0/N_samp;
				MarginalEntropySNP[i] += -ptmp*log(ptmp);
			}
			if (Na_genoMD[j][i]) {
				ptmp = (wsReal)Na_genoMD[j][i]/N_samp;
				MarginalEntropySNP_Y[i] += -ptmp*log(ptmp);
			}
			if (Na_genoMD[3+j][i]) {
				ptmp = (wsReal)Na_genoMD[3+j][i]/N_samp;
				MarginalEntropySNP_Y[i] += -ptmp*log(ptmp);
			}
		}
	}

	*Np_genoMD	= Na_genoMD;
	*Np_genoMDs	= Na_genoMDs;
}

inline void calcGenoJD(xBoostGeno *Xa_cnt, int N_snp, int N_caL, int N_ctL,
	int *Np_genoJD, int I, int J, wsUint **Na_genoMgnDist)
{
	int i1, i2, i3;

	register WISARD_pc count;
	//	uint64 tmp;

	// 00 00
	

	for (i1=0 ; i1<2 ; i1++) {
		for (i2=0 ; i2<2 ; i2++) {
			count = 0;
			for (i3=0 ; i3<N_caL ; i3++) {
				//count += bitCount(pgeno[j1*3 + i1].genocase[i3] & pgeno[j2*3 + i2].genocase[i3]);
				count += WISARD_POPCNT(Xa_cnt[I*2 + i1].ca[i3] & Xa_cnt[J*2 + i2].ca[i3]);
				//GenoDistr[i1*3 + i2] += bitCount(pgeno[j1*3 + i1].genocase[i3] & pgeno[j2*3 + i2].genocase[i3]);
			}
			Np_genoJD[i1*3 + i2] = (int)count;
			count = 0;

			for (i3=0 ; i3<N_ctL ; i3++) {
				//count += bitCount(pgeno[j1*3 + i1].genoctrl[i3] & pgeno[j2*3 + i2].genoctrl[i3]);
				count += WISARD_POPCNT(Xa_cnt[I*2 + i1].ct[i3] & Xa_cnt[J*2 + i2].ct[i3]);
				//GenoDistr[9 + i1*3 + i2] += bitCount(pgeno[j1*3 + i1].genoctrl[i3] & pgeno[j2*3 + i2].genoctrl[i3]);
			}
			Np_genoJD[9 + i1*3 + i2] = (int)count;
		}
	}
	//for case
	
	Np_genoJD[2] = Na_genoMgnDist[0][I] - Np_genoJD[0] - Np_genoJD[1];
 	Np_genoJD[5] = Na_genoMgnDist[1][I] - Np_genoJD[3] - Np_genoJD[4];

	Np_genoJD[6] = Na_genoMgnDist[0][J] - Np_genoJD[0] - Np_genoJD[3];
	Np_genoJD[7] = Na_genoMgnDist[1][J] - Np_genoJD[1] - Np_genoJD[4];

	Np_genoJD[8] = Na_genoMgnDist[2][J] - Np_genoJD[2] - Np_genoJD[5];

	//for ctrl
	Np_genoJD[11] = Na_genoMgnDist[3][I] - Np_genoJD[9] - Np_genoJD[10];
	Np_genoJD[14] = Na_genoMgnDist[4][I] - Np_genoJD[12] - Np_genoJD[13];

	Np_genoJD[15] = Na_genoMgnDist[3][J] - Np_genoJD[9] - Np_genoJD[12];
	Np_genoJD[16] = Na_genoMgnDist[4][J] - Np_genoJD[10] - Np_genoJD[13];

	Np_genoJD[17] = Na_genoMgnDist[5][J] - Np_genoJD[11] - Np_genoJD[14];
}

void _comp2(xBoostGeno *Xa_cnt, wsUintCst N_samp, wsUintCst N_snp, wsUint N_caL,
	wsUint N_ctL, wsUint *InteractionSNPpairs,
	wsUint **Na_genoMD, wsUintCst InteractionCount, wsReal *InteractionMeasureSNPpair,
	wsReal *Ra_zVal)
{
	wsUint	i, j, k;
	static double mu[3][3][2]	= {
		{ { 1,1 }, { 1,1 }, { 1,1 } },
		{ { 1,1 }, { 1,1 }, { 1,1 } },
		{ { 1,1 }, { 1,1 }, { 1,1 } }
	};
	static double mu0[3][3][2]	= {
		{ { 0,0 }, { 0,0 }, { 0,0 } },
		{ { 0,0 }, { 0,0 }, { 0,0 } },
		{ { 0,0 }, { 0,0 }, { 0,0 } }
	};
	static double mutmp[3][3][2];
	static double mu0tmp[3][3][2];

	static double mu_ij[3][3];
	static double mu_ik[3][2];
	static double mu_jk[3][2];
	double muError;
	int Na_genoJD[18];

	for (wsUint ii=0 ; ii<InteractionCount ; ii++) {
		wsUint j1 = InteractionSNPpairs[2*ii];
		wsUint j2 = InteractionSNPpairs[2*ii + 1];
		calcGenoJD(Xa_cnt, N_snp-1, N_caL, N_ctL,
			Na_genoJD, j1, j2, Na_genoMD);

		memcpy(mutmp, mu, 18*sizeof(double));
		memcpy(mu0tmp, mu0, 18*sizeof(double));

		//Iterative Proportional Fitting homogeneous model (see section 8.7.2 of Categorical Data Analysis (2nd Ed.))
		muError = 0.0;
		for (i=0 ; i<3 ; i++)
			for(j=0 ; j<3 ; j++)
				for (k=0 ; k<2 ; k++)
					muError += fabs(mutmp[i][j][k]-mu0tmp[i][j][k]);

		while (muError > 0.001) {
			memcpy(mu0tmp, mutmp, 18*sizeof(double)); //mu0tmp = mutmp;

			// mu_ij
			for (i = 0; i<3; i++)
				for (j = 0; j<3; j++)
					mu_ij[i][j] = mutmp[i][j][0] + mutmp[i][j][1];

			//mu_ijk = mu_ijk*n_ij/mu_ij
			for (i = 0; i<3; i++)
				for (j = 0; j <3; j++)
					for (k = 0; k <2; k++) {
				if (mu_ij[i][j]>0)
				{
					mutmp[i][j][k] = mutmp[i][j][k] *
						(Na_genoJD[3*i + j]+Na_genoJD[9 + 3*i + j])/mu_ij[i][j];
				} else
					mutmp[i][j][k] = 0;

					}

			// mu_ik
			for (i = 0; i<3; i++)
				for (k = 0; k<2; k++)
					mu_ik[i][k] = mutmp[i][0][k] + mutmp[i][1][k] + mutmp[i][2][k];

			//mu_ijk = mu_ijk*n_ik/mu_ik
			for (i = 0; i<3; i++)
				for (j = 0; j <3; j++)
					for (k = 0; k <2; k++) {
				// mu(i,j,k) = mu(i,j,k) * n_ik(i,k)/mu_ik(i,1,k);
				if (mu_ik[i][k] > 0) {
					mutmp[i][j][k] = mutmp[i][j][k] *
						(Na_genoJD[9*k + 3*i]+Na_genoJD[9*k + 3*i + 1] +
						Na_genoJD[9*k + 3*i + 2])/mu_ik[i][k];
				} else
					mutmp[i][j][k] = 0;
					}

			// mu_jk
			for (j = 0; j<3; j++)
				for (k = 0; k<2; k++)
					mu_jk[j][k] = mutmp[0][j][k] + mutmp[1][j][k] + mutmp[2][j][k];

			//mu_ijk = mu_ijk*n_jk/mu_jk
			for (i = 0; i<3; i++)
				for (j = 0; j <3; j++)
					for (k = 0; k <2; k++) {
				// mu(i,j,k) = mu(i,j,k) * n_jk(k,j)/mu_jk(1,j,k);
				if (mu_jk[j][k] > 0)
				{
					mutmp[i][j][k] = mutmp[i][j][k] *
						(Na_genoJD[9*k + j]+Na_genoJD[9*k + j + 3] +
						Na_genoJD[9*k + j + 6])/mu_jk[j][k];
				} else
					mutmp[i][j][k] = 0;
					}

			//calculate Error
			muError = 0.0;
			for (i = 0; i<3; i++)
				for (j = 0; j <3; j++)
					for (k = 0; k <2; k++)
						muError += fabs(mutmp[i][j][k]-mu0tmp[i][j][k]);
		}// end for while

		double tao = 0.0;
		double InteractionMeasure=0.0;
		for (i = 0; i<3; i++) // index for A
			for (j = 0; j<3; j++) //index for B
				for (k = 0; k<2; k++) { //index for C
					double ptmp1 = (double)Na_genoJD[k*9 + i*3 +j] / N_samp;
					if (ptmp1>0)
						InteractionMeasure += ptmp1 * log(ptmp1);
					double ptmp2 = mutmp[i][j][k] / N_samp;
					if (ptmp2>0) {
						InteractionMeasure += -ptmp1*log(ptmp2);
						tao += ptmp2;
					}
				}

		InteractionMeasure = (InteractionMeasure+log(tao))*N_samp*2;
// 		if (InteractionMeasure < 30)
// 		{
// 			//printf("SNPpair (%7d,%7d) : Approx: %f; Exact: %f\n", j1,j2, InteractionMeasureSNPpair[ii], InteractionMeasure);
// 		}

		InteractionMeasureSNPpair[ii] = InteractionMeasure; // update the interactionMeasure;


		//PLINK
		//     BB Bb  bb
		//     AA  a  b  c
		//     Aa  d  e  f
		//     aa  g  h  i
		//          B            b
		//     A  4a+2b+2d+e   4c+2b+2f+e
		//     a  4g+2h+2d+e   4i+2h+2f+e

		int *AlleleJointDistr = new int[8];

		//  case
		AlleleJointDistr[0] = 4*Na_genoJD[0] + 2*Na_genoJD[1] + 2*Na_genoJD[3] + Na_genoJD[4];
		AlleleJointDistr[1] = 4*Na_genoJD[2] + 2*Na_genoJD[1] + 2*Na_genoJD[5] + Na_genoJD[4];
		AlleleJointDistr[2] = 4*Na_genoJD[6] + 2*Na_genoJD[7] + 2*Na_genoJD[3] + Na_genoJD[4];
		AlleleJointDistr[3] = 4*Na_genoJD[8] + 2*Na_genoJD[7] + 2*Na_genoJD[5] + Na_genoJD[4];
		// control
		AlleleJointDistr[4] = 4*Na_genoJD[9] + 2*Na_genoJD[10] + 2*Na_genoJD[12] + Na_genoJD[13];
		AlleleJointDistr[5] = 4*Na_genoJD[11] + 2*Na_genoJD[10] + 2*Na_genoJD[14] + Na_genoJD[13];
		AlleleJointDistr[6] = 4*Na_genoJD[15] + 2*Na_genoJD[16] + 2*Na_genoJD[12] + Na_genoJD[13];
		AlleleJointDistr[7] = 4*Na_genoJD[17] + 2*Na_genoJD[16] + 2*Na_genoJD[14] + Na_genoJD[13];
		//
		double or_aff = log( (double)(AlleleJointDistr[0]*AlleleJointDistr[3])/ (double)(AlleleJointDistr[1]*AlleleJointDistr[2]) );
		double v_aff = 1/(double)AlleleJointDistr[0]
			+ 1/(double)AlleleJointDistr[1]
			+ 1/(double)AlleleJointDistr[2]
			+ 1/(double)AlleleJointDistr[3];

		double or_unf = log( (double)(AlleleJointDistr[4]*AlleleJointDistr[7])/ (double)(AlleleJointDistr[5]*AlleleJointDistr[6]) );
		double v_unf = 1/(double)AlleleJointDistr[4] + 1/(double)AlleleJointDistr[5] + 1/(double)AlleleJointDistr[6] + 1/(double)AlleleJointDistr[7];

 		if (fabs(v_aff) == numeric_limits<double>::infinity() ||
 			fabs(v_unf) == numeric_limits<double>::infinity())
			Ra_zVal[ii] = WISARD_NAN;
		else
// 			halt("ERR");

		//zval[ii] = Abs( (or_aff - or_unf) / sqrt ( v_aff + v_unf ) );
			Ra_zVal[ii] = ( (or_aff - or_unf) / sqrt ( v_aff + v_unf ) );

		delete [] AlleleJointDistr;

		if ((ii+1)%10000==0)
			notice("iteration [%d/%d]\r", ii+1, InteractionCount);
	}
}

double fastEpistasis(wsUint N_samp, wsReal *Ra_phe, char **Na_data, wsUintCst e1, wsUintCst e2)
{
	double z;  // statistic from either method
	wsUint person = 0;

	// Odds ratio test
	// make two 2x2 tables

	int a11, a12, a21, a22;
	int u11, u12, u21, u22;
	a11=a12=a21=a22=0;
	u11=u12=u21=u22=0;

	for (wsUint i=0 ; i<N_samp ; i++) {
		/* If missing */
		if (isMissingReal(Ra_phe[i])) {
			// Next person
			person++;
			continue;
		}

		char Q11[16] = { 4, 2, 0, 0, 2, 1, 0, };
		char Q12[16] = { 0, 0, 0, 0, 2, 1, 0, 0, 4, 2, 0, };
		char Q21[16] = { 0, 2, 4, 0, 0, 1, 2, 0, };
		char Q22[16] = { 0, 0, 0, 0, 0, 1, 2, 0, 0, 2, 4, 0, };
		char G1 = isMissing(Na_data[i][e1]) ? 3 : Na_data[i][e1];
		char G2 = isMissing(Na_data[i][e2]) ? 3 : Na_data[i][e2];
		wsUint G = (G1<<2) + G2;
		
		if (Ra_phe[i] == W1) { // if affected
			a11 += Q11[G];
			a12 += Q12[G];
			a21 += Q21[G];
			a22 += Q22[G];
		} else { // unaffected 
			u11 += Q11[G];
			u12 += Q12[G];
			u21 += Q21[G];
			u22 += Q22[G];
		}

		// Next person
		person++;
	}


	// Calculate log(OR) and SEs
	double or_aff, v_aff, or_unf, v_unf;

	or_aff = log( (double)(a11*a22)/ (double)(a12*a21) );
	v_aff = 1/(double)a11 + 1/(double)a12 + 1/(double)a21 + 1/(double)a22;

	// Case-only z-score (if requested)
// 	if (par::epi_caseonly)
// 		z = fabs( or_aff / sqrt(v_aff) );
//	else // Standard case-control analysis 
	{
		or_unf = log( (double)(u11*u22)/ (double)(u12*u21) );
		v_unf = 1/(double)u11 + 1/(double)u12 + 1/(double)u21 + 1/(double)u22;
		z = fabs( (or_aff - or_unf) / sqrt ( v_aff + v_unf ) );
	}


	//////////////////////////////
	// --nop option in effect 
	// Just output z score, if valid & above threshold

// 	if (par::epi_quickscan)
// 	{
// 		// Is this worth recording?
// 		if ( realnum(z) )
// 		{
// 			nepi++;
// 			if (z >= par::epi_alpha1)
// 				EPI << setw(4) << locus[e1]->chr << " " 
// 				<< setw(par::pp_maxsnp) << locus[e1]->name << " "
// 				<< setw(4) << locus[e2]->chr << " " 
// 				<< setw(par::pp_maxsnp) << locus[e2]->name << " "
// 				<< setw(12) << z*z << "\n";
// 			EPI.flush();
// 			continue;
// 		}
// 	}

	return z;
}

void cBoostAnalysis::_comp(wsUint I, wsUint N_samp, wsUint N_snp,
	/*int *DistrCollection, */wsReal &R_minItrc, wsReal &R_maxItrc,
	wsUint &buffersize, wsReal thresholdRecord, wsUint &N_cntItrc,
	wsUint **Na_itrcPairs, wsReal **Ra_itrcValPairs,
	wsUint *N_ca, wsUint *N_ct, wsUint *N_caL, wsUint *N_ctL,
	wsUint **Na_genoMD, wsUint **Na_genoMDs, wsReal *Ra_phe, char **Na_data,
	xFastEpi *Xp_e)
{
	wsUint ncase = N_ca[0];
	wsUint nctrl = N_ct[0];
	int Na_genoJD[18];

	wsReal *Pca = new wsReal[6];
	//P(C|A)
	/*index
			A=0		A=1		A=2
	C = 0    0		2		4
	C = 1	 1		3		5
	*/
	for (wsUint i=0 ; i<3 ; i++)// i for A
		for (wsUint j=0 ; j<2 ; j++) // j for C
			Pca[2*i + j] = (wsReal)Na_genoMD[j*3+i][I]/ Na_genoMDs[i][I];

	for (wsUint J=I+1 ; J<N_snp ; J++) {
		wsReal *Pab = new wsReal[9];
		wsReal *Pbc = new wsReal[6];

		if (OPT_ENABLED(quickepi)) {
			wsReal z = fastEpistasis(N_samp, Ra_phe, Na_data, I, J);
			/////////////////////////////////	
			// More full parsing of results
			//wsReal zero = 0;

			// Check this is a proper result
			if (z == z) {
				// One more test performed
				Xp_e->nepi++;

				// Count as a good result
				Xp_e->summary_good[I]++;
				Xp_e->summary_good[J]++;

				// Do we want to record this as part of the summary for the first set? 	  
				if (z >= 0.01) {
					// first variable will always be in A set
					Xp_e->summary_sig[I]++;
					// but the second may also be in A set
					Xp_e->summary_sig[J]++;

					wsStrCst Sp_c1 = getChrName2(Xp_e->Xa_snp[I].chr);
					wsStrCst Sp_c2 = getChrName2(Xp_e->Xa_snp[J].chr);
					Xp_e->C_fe->fmt("%s	%s	%d	%s	%s	%d	%g	%g\n",
						Sp_c1, Xp_e->Xa_snp[I].name, Xp_e->Xa_snp[I].pos,
						Sp_c2, Xp_e->Xa_snp[J].name, Xp_e->Xa_snp[J].pos,
						z, pnorm(z, 0, 1, 0) * 2.0);
				}

				// Is this result the best scrore yet for marker in set A?
				if (z > Xp_e->best_score[I] || NA(Xp_e->best_score[I])) {
					Xp_e->best_score[I] = z;
					Xp_e->best_partner[I] = J;
				}

				// The second marker might also be in set A 
				if (z > Xp_e->best_score[J] || NA(Xp_e->best_score[J])) {
					Xp_e->best_score[J] = z;
					Xp_e->best_partner[J] = I;
				}
			}
		}

		calcGenoJD(Xa_cnt, N_snp, N_caL[0], N_ctL[0], Na_genoJD, I, J, Na_genoMD);
/*
	GenoJointDistr: the index is as follows:
			AABB
	Case	0	1	2	3	4	5	6	7	8
	Ctrl	9	10	11	12	13	14	15	16	17
*/
		// P(A|B)
		/* index
				B = 0   B = 1	B = 2
		A = 0	0		3		6
		A = 1	1		4		7
		A = 2	2		5		8
		*/

		// i for B, j for A
		for (wsUint i=0 ; i<3 ; i++)
			for (wsUint j=0 ; j<3 ; j++)
				Pab[3*i+j] = (wsReal) (Na_genoJD[3*j+i] + Na_genoJD[3*j+i+9])
					/ Na_genoMDs[i][J];
					//pMarginalDistr[j2].MarginalDistrSNP[i];

		//P(B|C)
		/* index
				C = 0   C = 1
		B = 0	0		3
		B = 1	1		4
		B = 2	2		5
		*/
		//i for C (Y, class label), j for B

		for (wsUint j=0 ; j<3 ; j++)
			Pbc[j] = (wsReal) Na_genoMD[j][J]/ncase;
		for (wsUint j=0 ; j<3 ; j++)
			Pbc[3 + j] = (wsReal) Na_genoMD[j+3][J]/nctrl;

		// Papprx = Pab*Pbc*Pca
		// tao = sum(Papprx)
		// sum(p.* log(p) - p.*log(Pappr) + p.* log(tao))
		wsReal tao = 0.0;
		wsReal R_itrcVal=0.0;
		for (wsUint i=0 ; i<3 ; i++) { // index for A 
			for (wsUint j=0 ; j<3 ; j++) { //index for B
				for (wsUint k=0 ; k<2 ; k++) { //index for C
					wsReal ptmp1 = (wsReal)Na_genoJD[k*9 + i*3 +j] / N_samp;
					if (ptmp1>0)
						R_itrcVal += ptmp1 * log(ptmp1);

					wsReal ptmp2 = Pab[3*j+i]*Pbc[3*k+j]*Pca[2*i+k];
					if (ptmp2>0) {
						R_itrcVal += -ptmp1*log(ptmp2);
						tao += ptmp2;
					}
				}
			}
		}

		R_itrcVal = (R_itrcVal+log(tao))*N_samp*2;

		/* Range update */
		if (R_itrcVal > R_maxItrc)
			R_maxItrc = R_itrcVal;
		if (R_itrcVal < R_minItrc)
			R_minItrc = R_itrcVal;

		if (R_itrcVal > thresholdRecord) {
			if (N_cntItrc >= buffersize) { // buffersize is not enough, calloc new memory
				buffersize = buffersize * 2;
				wsUint *Na_nItrcPairs = NULL;
				wsAlloc(Na_nItrcPairs, wsUint, 2*buffersize);
				memcpy(Na_nItrcPairs, *Na_itrcPairs, sizeof(wsUint)*buffersize);
				DEALLOC(*Na_itrcPairs);
				*Na_itrcPairs = Na_nItrcPairs;

				wsReal *Ra_nItrcValPairs = NULL;
				wsAlloc(Ra_nItrcValPairs, wsReal, buffersize);
				memcpy(Ra_nItrcValPairs, *Ra_itrcValPairs, sizeof(wsReal)*(buffersize>>1));
				DEALLOC(*Ra_itrcValPairs);
				*Ra_itrcValPairs = Ra_nItrcValPairs;

				wsUint *Np_itrcPairs = *Na_itrcPairs;
				wsReal *Rp_itrcValPairs = *Ra_itrcValPairs;

//				Na_itrcPairs = (wsUint *)realloc(Na_itrcPairs, buffersize*2*sizeof(int));
//				Ra_itrcValPairs = (wsReal *)realloc(Ra_itrcValPairs, buffersize*sizeof(wsReal));

				Np_itrcPairs[2*N_cntItrc] = I;
				Np_itrcPairs[2*N_cntItrc + 1] = J;

				Rp_itrcValPairs[N_cntItrc] = R_itrcVal;
				N_cntItrc ++;
				//printf("InteractionCount: %6d\t(%6d,%6d)\t %f\n", InteractionCount,j1,j2, InteractionMeasure);
			} else {
				wsUint *Np_itrcPairs = *Na_itrcPairs;
				wsReal *Rp_itrcValPairs = *Ra_itrcValPairs;

				Np_itrcPairs[2*N_cntItrc] = I;
				Np_itrcPairs[2*N_cntItrc + 1] = J;

				Rp_itrcValPairs[N_cntItrc] = R_itrcVal;
				N_cntItrc ++;
				//printf("InteractionCount: %6d\t(%6d,%6d)\t %f\n", InteractionCount,j1,j2, InteractionMeasure);
			}
		}

		// collect distribution based on KSA
// 		if (InteractionMeasure > 100)
// 			DistrCollection[1000]++;
// 			//printf("pair:%d,%d\t %f\n", j1,j2, InteractionMeasure);
// 		else if (InteractionMeasure > 0)
// 			DistrCollection[(int)(InteractionMeasure/0.1)]++;
		delete [] Pab;
		delete [] Pbc;
	}
	delete [] Pca;
}

cBoostAnalysis::cBoostAnalysis(cIO *Cp_inpIO) : cAnalysis(Cp_inpIO)
{
	wsUintCst	N_pheno	= Cp_IO->sizePheno();
	wsUintCst	N_samp	= Cp_IO->sizeSample();
	wsUintCst	N_snp	= Cp_IO->sizeVariant();
	wsMat	Ra_phe	= Cp_IO->getPhenos();
	char**	Na_data	= Cp_IO->getGenotype();
	N_marker = N_snp;

	/*
	 *
	 * Build BOOST-like data structure
	 * 
	 */
	
	/* Step 1 : Count ca/ct */
	wsCalloc(N_ca, wsUint, N_pheno);
	wsCalloc(N_ct, wsUint, N_pheno);
	wsCalloc(N_caL, wsUint, N_pheno);
	wsCalloc(N_ctL, wsUint, N_pheno);
	for (wsUint i=0 ; i<N_pheno ; i++) {
		for (wsUint j=0 ; j<N_samp ; j++) {
			if ((wsUint)Ra_phe[i][j] == WISARD_AFFECTED)
				N_ca[i]++;
			else if ((wsUint)Ra_phe[i][j] == WISARD_UNAFFECTED)
				/* Control */ N_ct[i]++;
		}
		if (N_ca[i] == 0 || N_ct[i] == 0)
			halt("[%d]th phenotype have no case or no control!", i);
		N_caL[i] = (N_ca[i]+(WISARD_szpc)) >> WISARD_sftpc;
		N_ctL[i] = (N_ct[i]+(WISARD_szpc)) >> WISARD_sftpc;
	}

	/* Step 2 : Count genotype */
	wsUint N_sx = (N_samp+WISARD_szpc) >> WISARD_sftpc;
	wsAlloc(Xa_cnt, xBoostGeno, N_snp*2);
	/* Generate genotype buffer */
	for (wsUint i=0 ; i<N_snp ; i++) {
		wsCalloc(Xa_cnt[i*2].ca, WISARD_pc, N_sx);
		wsCalloc(Xa_cnt[i*2].ct, WISARD_pc, N_sx);
		wsCalloc(Xa_cnt[i*2+1].ca, WISARD_pc, N_sx);
		wsCalloc(Xa_cnt[i*2+1].ct, WISARD_pc, N_sx);
	}
	wsUint N_x[2] = { 0, };
	for (wsUint i=0 ; i<N_samp ; i++) {
		wsReal		R_v	= Ra_phe[0][i];
		if (isMissingReal(R_v)) continue;
		wsUint&		Q	= R_v==W0 ? N_x[0] : N_x[1];
		wsUint		X	= Q >> WISARD_sftpc;
		uint64_t	V	= (uint64_t)1 << ((Q++)&(WISARD_szpc-1));

		if (R_v == W0) for (wsUint j=0 ; j<N_snp ; j++) {
			switch (Na_data[i][j]) {
			case 0:
				Xa_cnt[j*2].ct[X] |= V;
				break;
			case 1:
				Xa_cnt[j*2+1].ct[X] |= V;
				break;
			}
		} else for (wsUint j=0 ; j<N_snp ; j++) {
			switch (Na_data[i][j]) {
			case 0:
				//LOG("[%d,%d]->0 , [%d,%d],%x->%x\n", i, j, j*2, X, V, Xa_cnt[j*2].ca[X]);
				Xa_cnt[j*2].ca[X] |= V;
				break;
			case 1:
				Xa_cnt[j*2+1].ca[X] |= V;
				break;
			}
		}
	}
}

cBoostAnalysis::~cBoostAnalysis()
{
	DEALLOC(N_ca);
	DEALLOC(N_ct);
	DEALLOC(N_caL);
	DEALLOC(N_ctL);
	 
	if (Xa_cnt) {
		/* Generate genotype buffer */
		for (wsUint i=0 ; i<N_marker ; i++) {
			DEALLOC(Xa_cnt[i*2].ca);
			DEALLOC(Xa_cnt[i*2].ct);
			DEALLOC(Xa_cnt[i*2+1].ca);
			DEALLOC(Xa_cnt[i*2+1].ct);
		}
		DEALLOC(Xa_cnt);		
	}
}

void cBoostAnalysis::run()
{
	wsUintCst		N_samp		= Cp_IO->sizeSample();
	wsUintCst		N_var		= Cp_IO->sizeVariant();
	wsMat		Ra_phe		= Cp_IO->getPhenos();
	char**		Na_data		= Cp_IO->getGenotype();
	vVariant&	Xa_snp		= Cp_IO->getVariant();

	/* Step 3 : Compute marginal */
	wsUint	**Na_genoMD		= NULL;
	wsUint	**Na_genoMDs	= NULL;
	wsReal	ptmp1			= (wsReal)N_ca[0]/N_samp;
	wsReal	R_mgnEtrpY		= -ptmp1 *log(ptmp1) - (1-ptmp1) *log(1-ptmp1);

	wsReal *Ra_mgnAssoc		= sseEmptyVec(N_var);
	wsReal *Ra_mgnEtrp		= sseEmptyVec(N_var);
	wsReal *Ra_mgnEtrpY		= sseEmptyVec(N_var);
	CalculateMarginalEntropyNdist(Na_data, N_var, N_samp, Ra_mgnEtrp,
		Ra_mgnEtrpY, Ra_phe[0], &Na_genoMD, &Na_genoMDs);

	/* Step 4 : Compute marginal association */
	for (wsUint i=0 ; i<N_var ; i++)
		Ra_mgnAssoc[i] = (-Ra_mgnEtrpY[i] + Ra_mgnEtrp[i] + R_mgnEtrpY)*N_samp*2;

	/* Export result */ {
		cExporter	*C_ma	= cExporter::summon("boost.marginal.res");
		vVariant		&Xa_snp	= Cp_IO->getVariant();
		headerVariant(C_ma);
		C_ma->put("	ASSOC\n");

		wsUint j=0;
		FOREACHDO (vVariant_it, Xa_snp, i, j++) {
			entryVariant(C_ma, *i);
			C_ma->fmt("	%g\n", Ra_mgnAssoc[j]);
		}
		delete C_ma;
	}
	sseFree(Ra_mgnEtrp);
	sseFree(Ra_mgnEtrpY);

	/* Set variables */
	if (!IS_ASSIGNED(thrboost)) {
		OPTION().assign("thrboost");
		OPTION().FORCE_OPT_REAL(thrboost);
	}

	/* For each combination */
	wsReal	R_maxItrc			= -9999999;
	wsReal	R_minItrc			= 9999999;
	//int *DistrCollection = NULL;
	wsUint	buffersize			= 50000;
	wsReal	thresholdRecord	= OPT_REAL(thrboost);
	wsUint	N_cntItrc			= 0;
	wsUint*	Na_itrcPairs		= NULL;
	wsCalloc(Na_itrcPairs, wsUint, buffersize*2);
	wsReal*	Ra_itrcValPairs		= NULL;
	wsCalloc(Ra_itrcValPairs, wsReal, buffersize);
	//MULTI_CALLOC(DistrCollection, int, 1001);

	/* Compute combination */
	xFastEpi V = { Cp_IO->getVariant(), 0, };
	if (OPT_ENABLED(quickepi)) {
		wsCalloc(V.best_partner, wsUint, N_var);
		wsCalloc(V.best_score, wsReal, N_var);
		wsCalloc(V.summary_good, wsUint, N_var);
		wsCalloc(V.summary_sig, wsUint, N_var);
		V.C_fe = cExporter::summon("qepi.res");
		V.C_fe->put("CHR1	VARIANT1	POS1	CHR2	VARIANT2	POS2	STAT_QEPI	P_QEPI\n");
		LOGoutput("Quick epistasis result is exported to [%s.qepi.res]\n",
			OPT_STRING(out));
	}

	LOG("[GxG] Combination calculation [%d choose 2] combinations...\n", N_var);
	xMaf*	Xp_maf = getIO()->getMAF();
	for (wsUint i=0 ; i<(N_var-1) ; i++) {
		/* Should have MAC */
		if (Xp_maf[i].N_allMac)
			_comp(i, N_samp, N_var, /*DistrCollection, */R_minItrc, R_maxItrc,
				buffersize, thresholdRecord, N_cntItrc, &Na_itrcPairs,
				&Ra_itrcValPairs, N_ca, N_ct, N_caL, N_ctL, Na_genoMD, Na_genoMDs,
				Ra_phe[0], Na_data, &V);
		else if (V.best_score[i] == W0)
			V.best_score[i] = WISARD_NAN;

		if ((i%10) == 0)
			notice("Combination [%d/%d]...\r", i, N_var);
	}
	LOG("Combination [%d/%d]...\n", N_var, N_var);

	if (OPT_ENABLED(quickepi)) {
		cExporter*	Cp_qe	= cExporter::summon("qepi.summary.res");

		LOGoutput("Summary of quick epistasis test is exported to [%s.qepi.summary.res]",
			OPT_STRING(out));
		headerVariant(Cp_qe);
		Cp_qe->put("	NSIG	NTOTAL	PROP	BESTCHISQ	BESTCHR	BESTSNP\n");

		for (wsUint i=0 ; i<N_var ; i++) {
			xVariant&	X_snp	= Xa_snp[i];

			entryVariant(Cp_qe, X_snp);
			if (V.best_score[i] != W0 && !NA(V.best_score[i])) {
				xVariant&	X_psnp	= Xa_snp[V.best_partner[i]];
				wsStrCst		Sp_c2	= getChrName2(X_psnp.chr);

				Cp_qe->fmt("	%d	%d	%g	%g	%s	%s\n", V.summary_sig[i],
					V.summary_good[i], (wsReal)V.summary_sig[i]/V.summary_good[i],
					V.best_score[i], Sp_c2, X_psnp.name);
			} else {
				Cp_qe->fmt("	%d	%d	%g	NA	NA	NA\n", V.summary_sig[i],
					V.summary_good[i], (wsReal)V.summary_sig[i]/V.summary_good[i]);
			}
		}
		delete Cp_qe;
	}

	if (N_cntItrc) {
		LOG("[GxG] p-value computation for [%d] chosen combinations...\n", N_cntItrc);
		wsReal*	Ra_zVal				= (wsReal *)calloc(N_cntItrc, sizeof(wsReal));	
		_comp2(Xa_cnt, N_samp, N_var, N_caL[0], N_ctL[0], Na_itrcPairs,
			Na_genoMD, N_cntItrc, Ra_itrcValPairs, Ra_zVal);

		/* Print interaction list */ {
			vVariant&		Xa_snp	= Cp_IO->getVariant();
			cExporter	*C_il	= cExporter::summon("boost.interaction.res");
			C_il->put("RANK	CHR1	SNP1	POS1	CHR2	SNP2	POS2	"
				"MARGIN1	MARGIN2	INTERACTION	ZSTAT\n");

			for (wsUint i=0 ; i<N_cntItrc ; i++) {
				if (Ra_itrcValPairs[i] > thresholdRecord && !NA(Ra_zVal[i])) {
					xVariant &s1 = Xa_snp[Na_itrcPairs[2*i]];
					xVariant &s2 = Xa_snp[Na_itrcPairs[2*i+1]];
					wsStrCst c1 = getChrName2(s1.chr);
					wsStrCst c2 = getChrName2(s2.chr);
					C_il->fmt("%d	%s	%s	%d	%s	%s	%d	%g	%g	%g	%g\n",
						i, c1, s1.name, s1.pos, c2, s2.name, s2.pos, 
						Ra_mgnAssoc[Na_itrcPairs[2*i]], Ra_mgnAssoc[Na_itrcPairs[2*i+1]],
						Ra_itrcValPairs[i], Ra_zVal[i]);
				}
			}

			delete C_il;
		}

		free(Ra_zVal);
	}
	sseFree(Ra_mgnAssoc);


	for (wsUint i=0 ; i<3 ; i++) {
		DEALLOC(Na_genoMD[i]);
		DEALLOC(Na_genoMDs[i]);
	}
	for (wsUint i=3 ; i<6 ; i++)
		DEALLOC(Na_genoMD[i]);
	DEALLOC(Na_genoMD);
	DEALLOC(Na_genoMDs);

	free(Ra_itrcValPairs);
	free(Na_itrcPairs);

	if (OPT_ENABLED(quickepi)) {
		DEALLOC(V.best_partner);
		DEALLOC(V.best_score);
		DEALLOC(V.summary_good);
		DEALLOC(V.summary_sig);
		delete V.C_fe;
	}
}

cMdrAnalysis::cMdrAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSetMgr) :
	cAnalysis(Cp_inpIO)
{
	Ra_sc			= NULL;
	Ra_permY		= NULL;
	Cp_anaSetMgr	= Cp_inpAnaSetMgr;
	Na_idxCase		= NULL;
	Na_idxCtrl		= NULL;
	N_anaCase		= 0;
	N_anaCtrl		= 0;

	/* --loocv and --cv */
	if (IS_ASSIGNED(cv) && OPT_ENABLED(loocv))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--cv", "--loocv");
	/* --cv and --nperm */
	if (IS_ASSIGNED(cv) && IS_ASSIGNED(nperm))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--cv", "--nperm");
	/* --hmdr and --gmdr */
	if (OPT_ENABLED(hmdr) && OPT_ENABLED(gmdr))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--hmdr", "--gmdr");
	/* --hmdr and --fuzzymdr */
	if (OPT_ENABLED(hmdr) && OPT_ENABLED(fuzzymdr))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--hmdr", "--fuzzymdr");
	/* --fuzzymdr and --gmdr */
	if (OPT_ENABLED(fuzzymdr) && OPT_ENABLED(gmdr))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--fuzzymdr", "--gmdr");

	/* Determine the order */
	wsUint N_order;
	Ra_sc	= NULL;
	if (!IS_ASSIGNED(order)) {
		OPTION().assign("order");
		OPTION().FORCE_OPT_NUMBER(order);
	}
	if (OPT_ENABLED(fuzzymdr)) {
		if (OPT_ENABLED(mdr)) {
			LOG("MDR analysis cannot concurrently performed with "
				"Fuzzy MDR analysis, MDR analysis will not performed\n");
		}

		/* Allocate MDR type */
		X_mdrType = MT_FUZZYMDR;
	} else if (IS_ASSIGNED(hmdr)) {
// 		if (IS_ASSIGNED(order)) {
// 			LOG("Hierarchical does not requires --order, it will be ignored\n");
// 
// 			OPTION().assign("order", "1");
// 			OPTION().FORCE_OPT_NUMBER(order);
// 		}
		if (OPT_ENABLED(mdr)) {
			LOG("MDR analysis cannot concurrently performed with "
				"Hierarchical MDR analysis, MDR analysis will not performed\n");
		}
		/* HMDR requires top N value. If not assigned, forcely assigns */
		if (!IS_ASSIGNED(top)) {
			OPTION().assign("top", OPTION().getDefVal("top"));
			OPTION().FORCE_OPT_NUMBER(top);
		}

		/* Allocate MDR type */
		X_mdrType = MT_HMDR;
	} else if (OPT_ENABLED(gmdr))
		/* Allocate MDR type */
		X_mdrType = MT_GMDR;
	else
		/* Allocate MDR type */
		X_mdrType = MT_MDR;

	/* If there is cv, top also should be in there */
	if (IS_ASSIGNED(cv)) {
		if (!IS_ASSIGNED(top))
			halt("--top should be defined to use --cv in MDR");
	}
	/* Gene-MDR */
	if (OPT_ENABLED(genemdr)) {
		ASSERT_OPTS_AND(set, setconsec);

		LOG("Gene-level MDR analysis is now performed\n");
	}

	/* If phenotype is continuous form */
	wsVecCst	Ra_Y	= Cp_IO->getPheno();
	wsUint	N_samp	= Cp_IO->sizeSample();

	/* Get # of samples */
	const char *Ba_misPheno = Cp_IO->getPheCovMissing(&N_anaSamp);
	sseMalloc(Ra_anaY, wsReal, N_anaSamp);
	if (Cp_IO->isContinuous() || IS_ASSIGNED(blup) || OPT_ENABLED(phenostdize)) {
		if (!OPT_ENABLED(gmdr))
			halt("MDR analysis for continuous phenotype requires --gmdr");

		/* Set GMDR flag */
		X_mdrType = MT_GMDR;
		/* Allocate score buffer */
		Ra_sc = sseVector(N_anaSamp);

		wsUint i, _i;
		/* If GMDR, fit model and get y-^y */
		wsUint	N_cov		= Cp_IO->sizeCovar();
		wsMat	Ra_cov		= Cp_IO->getCovariates();
		wsReal*	Ra_regY		= sseVector(N_anaSamp);
/**/	wsMat	Ra_anaXt	= sseMatrix(N_cov+1, N_anaSamp);
		for (i=0 ; i<N_anaSamp ; i++)
			Ra_anaXt[0][i] = W1;
		for (wsUint j=1 ; j<=N_cov ; j++) {
			for (i=_i=0 ; i<N_samp ; i++) {
				if (Ba_misPheno[i]) continue;
				Ra_anaXt[j][_i] = Ra_cov[j-1][i];
				_i++;
			}
		}
		for (i=_i=0 ; i<N_samp ; i++) {
			if (Ba_misPheno[i]) continue;
			Ra_regY[_i] = Ra_Y[i];
			_i++;
		}

		/* Initialize */
		wsReal *Ra_logisticV = NULL;
		wsVec	Ra_coef = sseEmptyVec(N_cov + 1);

		wsVec Ra_res = NULL;
		if (Cp_IO->isContinuous()) {
			cVector		V_y(Ra_regY, N_anaSamp);
			cStdMatrix	M_Xt(N_cov+1, N_anaSamp, Ra_anaXt);
			cSymMatrix	M_XtX	= M_Xt.Mt();
			cVector		V_Xty	= M_Xt * V_y;
			cSymMatrix&	M_XtXi	= M_XtX.inv();
			cVector		V_bHat	= M_XtXi * V_Xty;
			delete& M_XtXi;
			//M_XtXi.rem();

			cVector		V_yHat	= V_bHat * M_Xt;	/* yHat <- X%*%matrix(solve(XX)%*%Xy,ncol=1) */
			cVector		V_resid	= V_y - V_yHat;		/* residual <- y - yHat */
			V_resid.setDontDealloc();
			Ra_res = V_resid.get();
		} else {
			wsReal *Ra_phiHat = NULL;
			sseMalloc(Ra_logisticV, wsReal, N_anaSamp);
			bool B_conv = false;
			/*int N_iter = */cRegrAnalysis::fitLogisticByNR(N_anaSamp, N_cov+1, Ra_anaXt,
				Ra_anaY, Ra_coef, Ra_logisticV, 20, B_conv, &Ra_phiHat);
			/* Use raw residual (See Lou XY GMDR impl) */
			sseVsV(Ra_anaY, Ra_phiHat, Ra_sc, N_anaSamp);
		}

		for (wsUint i=0,I=0 ; i<N_anaSamp ; i++) {
			if (isMissingReal(Ra_Y[i])) continue;

			Ra_sc[I] = Ra_res[I];
			if (Ra_sc[I] < W0) {
				Ra_anaY[I++] = WISARD_UNAFFECTED;
				N_anaCtrl++;
			} else {
				Ra_anaY[I++] = WISARD_AFFECTED;
				N_anaCase++;
			}
		}
		sseFree(Ra_res);
		sseFree(Ra_coef);
	} else for (wsUint i=0,I=0 ; i<N_anaSamp ; i++) {
		if (isMissingReal(Ra_Y[i])) continue;

		if (Ra_Y[i] == WISARD_UNAFFECTED) {
			Ra_anaY[I++] = WISARD_UNAFFECTED;
			N_anaCtrl++;
		} else {
			Ra_anaY[I++] = WISARD_AFFECTED;
			N_anaCase++;
		}
	}
	DEALLOC(Ba_misPheno);
	if (N_anaCase == 0 || N_anaCtrl == 0)
		halt("Either of # of affected [%d] or unaffected [%d] is zero, cannot proceed!",
			N_anaCase, N_anaCtrl);
	N_order = OPT_NUMBER(order);
	LOG("%sMDR is now performed (order %d)\n", X_mdrType == MT_GMDR?"Generalized ":
		(X_mdrType == MT_FUZZYMDR?"Fuzzy ":""),
		N_order);

	/* If there is --cv */
	if (IS_ASSIGNED(cv)) {
		wsUint N_cv = (wsUint)OPT_NUMBER(cv);
		/* FInd out # of case/ctrl & record their pos */
		N_anaCase = N_anaCtrl = 0;
		wsAlloc(Na_idxCase, wsUint, N_anaSamp);
		wsAlloc(Na_idxCtrl, wsUint, N_anaSamp);
		LOOP (i, N_anaSamp) {
			if (Ra_anaY[i] == WISARD_AFFECTED)
				Na_idxCase[N_anaCase++] = i;
			else if (Ra_anaY[i] == WISARD_UNAFFECTED)
				Na_idxCtrl[N_anaCtrl++] = i;
			else halt("Phenotype of [%d]th sample is not case nor control!",
				i);
		}
		/* #--cv should not exceed #case/#ctrl */
		if (N_cv>N_anaCase || N_cv>N_anaCtrl)
			halt("# of cross-validation [%d] cannot exceed the number of "
				"case[%d] or control[%d]!", N_cv, N_anaCase, N_anaCtrl);
	}

	wsUint N_res = IS_ASSIGNED(cv) ? OPT_NUMBER(cv) : 1;
	wsAlloc(Cp_res, cResult*, N_res);
	LOOP (i, N_res) {
		if (IS_ASSIGNED(top))
			Cp_res[i] = new cTopResult(OPT_NUMBER(top), N_order, RES_SIZE, 0, !IS_ASSIGNED(nperm));
		else {
			wsUint N_mkr = Cp_IO->sizeVariant();
			wsUint N_comb = nCk(N_mkr, N_order);
			Cp_res[i] = new cAllResult(N_comb, N_order, RES_SIZE);
		}
	}
}

cMdrAnalysis::~cMdrAnalysis()
{
	wsUint N_res = IS_ASSIGNED(cv) ? OPT_NUMBER(cv) : 1;
	LOOP(i, N_res)
		delete Cp_res[i];
	DEALLOC(Cp_res);
	sseFree(Ra_sc);
	DEALLOC(Na_idxCase);
	DEALLOC(Na_idxCtrl);
	sseUnmat(Ra_permY, OPT_NUMBER(nperm));
}

inline void _mdrY(xMdr *Xp_mdr, wsUint N_cell, wsVecCst Ra_Y, char** Na_data,
	wsUint* Na_idxVrt, vInt_it X_curTestCaIdx, vInt_it X_curTestCtIdx,
	vInt& Nv_addIdxVrt, vInt& Na_cellIdx, vInt& Na_newCellIdx)
{
	wsVecCst	Ra_sc	= Xp_mdr->Ra_sc;
	wsVec	Ra_scca	= Xp_mdr->Ra_scca;
	wsVec	Ra_scct	= Xp_mdr->Ra_scct;

	memset(Xp_mdr->Na_cntCase[0], 0x00, sizeof(wsUint)*N_cell);
	memset(Xp_mdr->Na_cntCtrl[0], 0x00, sizeof(wsUint)*N_cell);
	if (Xp_mdr->B_haveTest) {
		memset(Xp_mdr->Na_cntCase[1], 0x00, sizeof(wsUint)*N_cell);
		memset(Xp_mdr->Na_cntCtrl[1], 0x00, sizeof(wsUint)*N_cell);
	}

	/* Initiate score */
	if (Xp_mdr->X_type & MT_GMDR) {
		memset(Ra_scca, 0x00, sizeof(wsReal)*N_cell);
		memset(Ra_scct, 0x00, sizeof(wsReal)*N_cell);
	}

	LOOP (j, Xp_mdr->N_sample) {
		wsUint	B_isCurTest	= 0;
		wsReal	Y			= Ra_Y[j];

		if (isMissingReal(Y) || (Na_cellIdx.size() && Na_cellIdx[j] == -1)) continue;

		/* If Na_testIdxCtrl/Case and its index, then continue */
		if (X_curTestCaIdx != Xp_mdr->Nv_testIdxCase.end() && *X_curTestCaIdx == (int)j) {
			X_curTestCaIdx++;
			B_isCurTest = 1;
		}
		if (X_curTestCtIdx != Xp_mdr->Nv_testIdxCtrl.end() && *X_curTestCtIdx == (int)j) {
			X_curTestCtIdx++;
			B_isCurTest = 1;
		}

		wsUint N_idx = 0;
		LOOP (k, Xp_mdr->N_order) {
			wsUint	N_idxMkr	= Na_idxVrt[k];
			if (Xp_mdr->Na_mapMarker)
				N_idxMkr		= Xp_mdr->Na_mapMarker[N_idxMkr];
			char	N_geno		= Na_data[j][N_idxMkr];

			if (isMissing(N_geno)) {
				if (Na_newCellIdx.size()) Na_newCellIdx[j] = -1;
				goto _end;
			}
			N_idx += N_geno;
			N_idx *= 3; /* 0/1/2 coding */
		}
		N_idx /= 3;
		if (N_idx >= N_cell)
			halt_fmt(WISARD_SYST_INVL_MDR_CELLRANGE, N_idx, N_cell);

		/* Counting */
		Y==WISARD_AFFECTED ?
			Xp_mdr->Na_cntCase[B_isCurTest][N_idx]++ :
			Xp_mdr->Na_cntCtrl[B_isCurTest][N_idx]++;
		/* Store index if needed */
		if (Na_newCellIdx.size()) Na_newCellIdx[j] = N_idx;

		/* Counting score */
		if (Xp_mdr->X_type & MT_GMDR)
			Y==WISARD_AFFECTED ? Ra_scca[N_idx] += Ra_sc[j] :
				Ra_scct[N_idx] += Ra_sc[j];
_end:;
	}
}

inline void _mdrY(xMdr *Xp_mdr, wsUint N_cell, wsVecCst Ra_Y, char** Na_data,
	vInt& Nv_idxVrt, vInt_it X_curTestCaIdx, vInt_it X_curTestCtIdx,
	vInt& Nv_addIdxVrt, vInt& Na_cellIdx, vInt& Na_newCellIdx)
{
	wsVecCst	Ra_sc	= Xp_mdr->Ra_sc;
	wsVec	Ra_scca	= Xp_mdr->Ra_scca;
	wsVec	Ra_scct	= Xp_mdr->Ra_scct;
	char	B_add	= 1;

	/* If Nv_addIdxVrt */
	vInt*	Nvp_idx	= NULL;
	if (Nv_addIdxVrt.size()) {
		/* If there is no cache, add Nv_idxVrt to Nv_addIdxVrt */
		if (!Na_cellIdx.size()) {
			FOREACH (vInt_it, Nv_idxVrt, i)
				Nv_addIdxVrt.push_back(*i);
			B_add = 0;
		}			
		Nvp_idx = &Nv_addIdxVrt;
	} else
		Nvp_idx = &Nv_idxVrt;
	vInt&	Nv_idx	= *Nvp_idx;

	memset(Xp_mdr->Na_cntCase[0], 0x00, sizeof(wsUint)*N_cell);
	memset(Xp_mdr->Na_cntCtrl[0], 0x00, sizeof(wsUint)*N_cell);
	if (Xp_mdr->B_haveTest) {
		memset(Xp_mdr->Na_cntCase[1], 0x00, sizeof(wsUint)*N_cell);
		memset(Xp_mdr->Na_cntCtrl[1], 0x00, sizeof(wsUint)*N_cell);
	}

	/* Initiate score */
	if (Xp_mdr->X_type & MT_GMDR) {
		memset(Ra_scca, 0x00, sizeof(wsReal)*N_cell);
		memset(Ra_scct, 0x00, sizeof(wsReal)*N_cell);
	}

	LOOP (j, Xp_mdr->N_sample) {
		wsUint	B_isCurTest	= 0;
		wsReal	Y			= Ra_Y[j];

		if (isMissingReal(Y) || (Na_cellIdx.size() &&
			Na_cellIdx[j] == -1)) continue;

		/* If Na_testIdxCtrl/Case and its index, then continue */
		if (X_curTestCaIdx != Xp_mdr->Nv_testIdxCase.end() && *X_curTestCaIdx == (int)j) {
			X_curTestCaIdx++;
			B_isCurTest = 1;
		}
		if (X_curTestCtIdx != Xp_mdr->Nv_testIdxCtrl.end() && *X_curTestCtIdx == (int)j) {
			X_curTestCtIdx++;
			B_isCurTest = 1;
		}

		wsUint N_idx = Na_cellIdx.size() ? Na_cellIdx[j] * 3 : 0;
		FOREACH (vInt_it, Nv_idx, k) {
			wsUint	N_idxMkr	= *k;
			if (Xp_mdr->Na_mapMarker)
				N_idxMkr		= Xp_mdr->Na_mapMarker[N_idxMkr];
			char	N_geno		= Na_data[j][N_idxMkr];

			if (isMissing(N_geno)) {
				if (Na_newCellIdx.size()) Na_newCellIdx[j] = -1;
				goto _end;
			}
			N_idx += N_geno;
			N_idx *= 3; /* 0/1/2 coding */
		}
		N_idx /= 3;
		if (N_idx >= N_cell)
			halt_fmt(WISARD_SYST_INVL_MDR_CELLRANGE, N_idx, N_cell);

		/* Counting */
		Y==WISARD_AFFECTED ?
			Xp_mdr->Na_cntCase[B_isCurTest][N_idx]++ :
		Xp_mdr->Na_cntCtrl[B_isCurTest][N_idx]++;
		/* Store index if needed */
		if (Na_newCellIdx.size()) Na_newCellIdx[j] = N_idx;

		/* Counting score */
		if (Xp_mdr->X_type & MT_GMDR)
			Y==WISARD_AFFECTED ? Ra_scca[N_idx] += Ra_sc[j] :
			Ra_scct[N_idx] += Ra_sc[j];
_end:;
	}

	if (Nv_addIdxVrt.size() && B_add) FOREACH (vInt_it, Nv_idxVrt, i)
		Nv_addIdxVrt.push_back(*i);
}


inline float _measure(xMdr* Xp_mdr, wsUint N_cell, char B_haveTest)
{
	wsVec	Ra_scca	= Xp_mdr->Ra_scca;
	wsVec	Ra_scct	= Xp_mdr->Ra_scct;
	wsUint	N_TP=0, N_FP=0, N_FN=0, N_TN=0;
	wsReal	R_TP=W0, R_FP=W0, R_FN=W0, R_TN=W0;

	wsUint*	Na_cca	= Xp_mdr->Na_cntCase[B_haveTest ? 1 : 0];
	wsUint*	Na_cct	= Xp_mdr->Na_cntCtrl[B_haveTest ? 1 : 0];

	switch (Xp_mdr->X_type) {
	case MT_FUZZYMDR: {
			wsVec Ra_1 = sseVector(N_cell);
			wsVec Ra_2 = sseVector(N_cell);
			wsVec Ra_3 = sseVector(N_cell);
			wsVec Ra_4 = sseVector(N_cell);
			wsVec Ra_5 = sseVector(N_cell);
			wsVec Ra_6 = sseVector(N_cell);
			wsRealCst W05 = REAL_CONST(0.5);
			wsRealCst W_05 = REAL_CONST(-0.5);

			LOOP(j, N_cell) {
				Ra_1[j] = Xp_mdr->Na_cntCase[0][j];
				Ra_2[j] = Xp_mdr->Na_cntCtrl[0][j];
			}
			wsReal R_ratio = sseVsum(Ra_1, N_cell) / sseVsum(Ra_2, N_cell);

			LOOP (j, N_cell) {
				// if (counts[c, 2]==0){ counts[c, 3] <-counts[c, 1]/(counts[c, 2]+0.1) }
				if (Xp_mdr->Na_cntCtrl[0][j] == 0)
					Ra_3[j] = Xp_mdr->Na_cntCase[0][j] / ((wsReal)Xp_mdr->Na_cntCtrl[0][j] + REAL_CONST(0.1));
				//  else{ counts[, 3] <-counts[, 1]/counts[, 2] }
				else
					Ra_3[j] = Xp_mdr->Na_cntCase[0][j] / (wsReal)Xp_mdr->Na_cntCtrl[0][j];

				// if (counts[c, 1]==0){ counts[c, 4]<-log((counts[c, 3]+0.1)/ratio, 10) }
				if (Xp_mdr->Na_cntCase[0][j] == 0)
					Ra_4[j] = log10((Ra_3[j]+REAL_CONST(0.1))/R_ratio);
				// else{ counts[c, 4] <-log(counts[c, 3]/ratio, 10) } #hy
				else
					Ra_4[j] = log10(Ra_3[j]/R_ratio);
				// if (counts[c, 4]<=-0.5){
				if (Ra_4[j] <= W_05)
					Ra_5[j] = W1;
					// counts[c, 5] <-1
				// } else{
				else {
					// if (counts[c, 4]>=0.5) {counts[c,5]<-0}
					if (Ra_4[j] >= W05)
						Ra_5[j] = W0;
					// else{ counts[c, 5]<-1/
					//	(
					//		(
					//			(counts[c, 4]-(-0.5)) / (counts[c, 4]-(-0.5))
					//		)^2
					//		+
					//		(
					//			(counts[c, 4]-(-0.5))/(counts[c, 4]-0.5)
					//		)^2
					//	) }
					else
						Ra_5[j] = W1/(
							pow((Ra_4[j]-W_05)/(Ra_4[j]-W_05), 2) +
							pow((Ra_4[j]-W_05)/(Ra_4[j]-W05), 2)
							);
				}

				// if (counts[c, 4]>=0.5){ counts[c, 6] <-1 } else{
				if (Ra_4[j]>=W05)
					Ra_6[j] = W1;
				else {
					// if (counts[c, 4]<=-0.5){ counts[c, 6] <-0 }
					if (Ra_4[j]<=W_05) Ra_6[j] = W0;
					else Ra_6[j] = W1/(
						pow((Ra_4[j]-W05)/(Ra_4[j]-W_05), 2) +
						pow((Ra_4[j]-W05)/(Ra_4[j]-W05), 2)
						);
					// else{ counts[c, 6]<-1/(((counts[c, 4]-0.5)/(counts[c, 4]-(-0.5)))^2+
					// ((counts[c, 4]-0.5)/(counts[c, 4]-0.5))^2) }
				}
			}
			R_TP = sseVV(Ra_1, N_cell, Ra_6);
			R_FP = sseVV(Ra_2, N_cell, Ra_6);
			R_TN = sseVV(Ra_2, N_cell, Ra_5);
			R_FN = sseVV(Ra_1, N_cell, Ra_5);

			sseFree(Ra_1);
			sseFree(Ra_2);
			sseFree(Ra_3);
			sseFree(Ra_4);
			sseFree(Ra_5);
			sseFree(Ra_6);
		}
		return (float)(R_TP+R_TN)/(float)(R_TP+R_TN+R_FP+R_FN);
	default:
		LOOP (j, N_cell) {
			wsReal R_ratio = Xp_mdr->X_type & MT_GMDR ? Ra_scca[j] + Ra_scct[j]
				: (Xp_mdr->Na_cntCase[0][j]>Xp_mdr->Na_cntCtrl[0][j] ? W1 : REAL_CONST(-1.0));
			if (R_ratio > W0) {
				/* High risk cell,
					* case in here should be TP
					* controls in here should be FP */
				N_TP += Na_cca[j];
				N_FP += Na_cct[j];
			} else {
				/* Low risk cell,
					* case in here should be FN
					* controls in here should be TN */
				N_FN += Na_cca[j];
				N_TN += Na_cct[j];
			}
		}
		/* In default, return BA */
		return (float)(N_TP+N_TN)/(float)(N_TP+N_TN+N_FP+N_FN);
	}

	/* Should not reach to here */
	return (float)W0;
}

inline void _mdr(xMdr* Xp_mdr, wsUint N_cell, wsUint N_comb, wsUint *Np_data)
{
	vInt	tmp;
	/* FIXME : MDR should support scca/scct */

// 	VEC_c	Ra_sc	= Xp_mdr->Ra_sc;
// 	VEC_t	Ra_scca	= Xp_mdr->Ra_scca;
// 	VEC_t	Ra_scct	= Xp_mdr->Ra_scct;
	wsUint	N_perm	= OPT_NUMBER(nperm);
	char**	Na_data	= Xp_mdr->Na_data;
	/* Is it train/test scheme? */
	char B_haveTest = Xp_mdr->Nv_testIdxCase.size() || Xp_mdr->Nv_testIdxCtrl.size();

	for (wsUint i=0 ; i<N_comb ; i++,Np_data+=Xp_mdr->N_order+RES_SIZE) {
		// Current combination of markers is stored in Np_data ~ Np_data+order-1
		vInt_it X_curTestCaIdx = Xp_mdr->Nv_testIdxCase.begin();
		vInt_it X_curTestCtIdx = Xp_mdr->Nv_testIdxCtrl.begin();

		/* Using permutation */
		if (Xp_mdr->Ra_permY) {
			/* Observed Y */
			_mdrY(Xp_mdr, N_cell, Xp_mdr->Ra_Y, Na_data, Np_data,
				X_curTestCaIdx, X_curTestCtIdx, tmp, tmp, tmp);
			wsReal R_obsBA = (wsReal)_measure(Xp_mdr, N_cell, B_haveTest);
			wsUint N_over = 1;

			LOOP (j, N_perm) {
				_mdrY(Xp_mdr, N_cell, Xp_mdr->Ra_permY[j], Na_data, Np_data,
					X_curTestCaIdx, X_curTestCtIdx, tmp, tmp, tmp);

				// Classify high/low
				wsReal R_permBA = (wsReal)_measure(Xp_mdr, N_cell, B_haveTest);
				if (R_permBA >= R_obsBA) N_over++;
			}

			/* Here the result will be stored */
			float *Rp_res = ((float *)Np_data)+Xp_mdr->N_order;
			Rp_res[0] = (float)N_over / (float)(N_perm + 1);
		} else {
			_mdrY(Xp_mdr, N_cell, Xp_mdr->Ra_Y, Na_data, Np_data,
				X_curTestCaIdx, X_curTestCtIdx, tmp, tmp, tmp);

			// Classify high/low
			float R_ba = _measure(Xp_mdr, N_cell, B_haveTest);

			/* Here the result will be stored */
			float *Rp_res = ((float *)Np_data)+Xp_mdr->N_order;
			Rp_res[0] = R_ba;
		}
	}
}

int thr_mdr(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	xThreadMdr*		Xp_mdr		= (xThreadMdr *)Vp_shareData;
	cMdrAnalysis*	Cp_anaMdr	= Xp_mdr->Cp_anaMdr;
	cIO*			Cp_IO		= Cp_anaMdr->getIO();
	wsUint			N_sample	= Cp_IO->sizeSample();
	xMdrType		X_type		= Cp_anaMdr->getMdrType();
	char**			Na_data		= Cp_IO->getGenotype();
	wsVecCst		Ra_Y		= Cp_anaMdr->getAnaY();
	wsVecCst		Ra_sc		= Cp_anaMdr->getScore();
	wsMat			Ra_permY	= Cp_anaMdr->getPermY();

	wsUint	*Np_data	= (wsUint *)Vp_data;
	/* First element of thread-specific data is equivalent to
	 * the number of combinations should be processed in this thread */
	wsUint	N_comb		= *(Np_data++);	
	wsUint	N_order		= OPT_NUMBER(order);			///< Order of combination
	wsUint	N_cell		= (wsUint)pow(3.0, N_order);	///< Number of possible cells for given order of combination

	wsUmat	Na_case		= sseUmatrix(2, N_cell); /* Number of cases for each cell */
	wsUmat	Na_ctrl		= sseUmatrix(2, N_cell); /* Number of ctrls for each cell */
	wsReal	*Ra_scca	= NULL;
	wsReal	*Ra_scct	= NULL;
	if (X_type == MT_GMDR) {
		Ra_scct = sseEmptyVec(N_cell);
		Ra_scca = sseEmptyVec(N_cell);
	}

	xMdr X_mdr = {
		X_type,
		(bool)IS_ASSIGNED(cv),
		N_order,
		N_sample,
		Na_data,
		Ra_Y,
		Na_case,
		Na_ctrl,
		NULL,
		Ra_sc, Ra_scca, Ra_scct,
		Xp_mdr->Na_testIdxCase, Xp_mdr->Na_testIdxCtrl,
		Ra_permY,
	};
	_mdr(&X_mdr, N_cell, N_comb, Np_data);

	/* Cleanup */
	sseUnmat(Na_case, 2);
	sseUnmat(Na_ctrl, 2);
	sseFree(Ra_scca);
	sseFree(Ra_scct);

	/* This means nothing */
	return 0;
}

int dist_mdr(int N_mode, int N_thread, void *Vp_thrData,
	wsUint L_thrData, void *Vp_shareData, wsUint *Na_waitQueue)
{
	wsUint	i=0;
	wsUint	*Np_curc;
	wsUint	*Np_comb	= (wsUint *)Vp_thrData;
	cMdrAnalysis*
			Cp_anaMdr	= ((xThreadMdr *)Vp_shareData)->Cp_anaMdr;
	wsUint	N_cv		= ((xThreadMdr *)Vp_shareData)->N_curCV;
	cIO		*Cp_IO		= Cp_anaMdr->getIO();
	cResult	*Cp_res = Cp_anaMdr->getResult(N_cv);
	static cComb
			C_comb;

	switch (N_mode) {
	case DISTFUN_INIT:
		C_comb.init(OPT_NUMBER(order), Cp_IO->sizeVariant());
		break;
	case DISTFUN_DISTRIBUTE:
		Np_curc = Np_comb;
		for (i=0,Np_curc=Np_comb ; i<(wsUint)N_thread ; i++,
			Np_curc+=COMB_SIZE*(OPT_NUMBER(order)+RES_SIZE)+1) {
			/* Generate combination */
			wsUint N_genComb = C_comb.get(COMB_SIZE, Np_curc+1, OPT_NUMBER(order)+RES_SIZE);
			*Np_curc = N_genComb;
			if (N_genComb != COMB_SIZE) {
				i += N_genComb>0;
				break;
			}
		}
		break;
	case DISTFUN_AFTERLOOP:
		Np_curc = Np_comb;
		for (i=0,Np_curc=Np_comb ; i<(wsUint)N_thread ; i++,
			Np_curc+=COMB_SIZE*(OPT_NUMBER(order)+RES_SIZE)+1) {
			/* Generate combination */
			wsUint N_genComb = *Np_curc;
			Cp_res->update((float *)(Np_curc+1), N_genComb);
		}
		break;
	case DISTFUN_UNINIT:
		C_comb.uninit();
		break;
	default:
		halt("Unsupported distFunc command [%d]", N_mode);
		break;
	}

	return i;
}

inline wsUint fround(wsReal v)
{
	return (wsUint)(v + REAL_CONST(0.5));
}
void cMdrAnalysis::run()
{
	vVariant&	Xv_var	= getIO()->getVariant();
	wsUint		N_order	= OPT_NUMBER(order);
	wsUint		N_cxval	= IS_ASSIGNED(cv) ? OPT_NUMBER(cv) : 1;
	wsReal		R_cxval	= (wsReal)N_cxval;
	char**		Na_data	= Cp_IO->getGenotype();

	LOG("Exhaustive MDR analysis with %d-order start\n", N_order);
	if (IS_ASSIGNED(cv))
		LOG("Cross-validation count is [%d]\n", OPT_NUMBER(cv));
	else if (IS_ASSIGNED(nperm)) {
		LOG("Permutation count for MDR is [%d]\n", OPT_NUMBER(nperm));
		Ra_permY = sseMatrix(OPT_NUMBER(nperm), N_anaSamp);
		wsVec Ra_subY = sseVector(N_anaSamp);
		LOOP (j, (wsUint)OPT_NUMBER(nperm)) {
// 			LOOP (k, N_anaSamp) Ra_subY[k] = k;
// 			for (wsUint k=0 ; k<N_anaSamp ; k++) {
// 				wsUint N_idx = wsRand()%(N_anaSamp-k);
// 				Ra_permY[j][k] = Ra_subY[N_idx];
//				memmove(Ra_subY+N_idx, Ra_subY+N_idx+1, sizeof(wsUint)*(N_anaSamp-N_idx-1));
//			}
			memcpy(Ra_subY, Ra_anaY, sizeof(wsReal)*N_anaSamp);
			for (wsUint k=0 ; k<N_anaSamp ; k++) {
				wsUint N_idx = wsRand()%(N_anaSamp-k);
				Ra_permY[j][k] = Ra_subY[N_idx];
				memmove(Ra_subY+N_idx, Ra_subY+N_idx+1, sizeof(wsUint)*(N_anaSamp-N_idx-1));
			}
		}					
	}

	if (IS_ASSIGNED(hmdr))
		LOG("Hierarchical MDR analysis start with initial order [%d], "
		"top [%d] results\n", OPT_NUMBER(order), OPT_NUMBER(top));

	if (!IS_ASSIGNED(hmdrprior)) LOOP(cv, N_cxval) {
		vInt Na_testIdxCase;
		vInt Na_testIdxCtrl;
		/* Make test set indices if cv > 1 */
		if (N_cxval > 1) {
			wsUint	N_sCase = fround((wsReal)N_anaCase / R_cxval * (wsReal)cv);
			wsUint	N_eCase = fround((wsReal)N_anaCase / R_cxval * (wsReal)(cv+1));
			wsUint	N_sCtrl = fround((wsReal)N_anaCtrl / R_cxval * (wsReal)cv);
			wsUint	N_eCtrl = fround((wsReal)N_anaCtrl / R_cxval * (wsReal)(cv+1)); 

			/* Insert [s,e) case/ctrl */
			for (wsUint i=N_sCase ; i<N_eCase ; i++) Na_testIdxCase.push_back(Na_idxCase[i]);
			for (wsUint i=N_sCtrl ; i<N_eCtrl ; i++) Na_testIdxCtrl.push_back(Na_idxCtrl[i]);
		}

		/* Do run */
		if (OPT_ENABLED(genemdr)) {
			wsUint	N_cell	= (wsUint)pow(3.0, N_order);	///< Number of possible cells for given order of combination

			wsUmat	Na_case	= sseUmatrix(2, N_cell); /* Number of cases for each cell */
			wsUmat	Na_ctrl	= sseUmatrix(2, N_cell); /* Number of ctrls for each cell */
			wsVec	Ra_scca	= NULL;
			wsVec	Ra_scct	= NULL;
			if (X_mdrType == MT_GMDR) {
				Ra_scct = sseVector(N_cell);
				Ra_scca = sseVector(N_cell);
			}

			/* If --cv && --shuffle */
			if (OPT_ENABLED(shuffle)) {
				/* FIXME */
			}

			/* Do gene-wise, and choose best */
			mGeneDef& Xm_gd = Cp_anaSetMgr->getGeneDef();

			/* For all genes, do MDR */
			wsUint	*Np_comb = NULL;
			wsAlloc(Np_comb, wsUint, (RES_SIZE+OPT_NUMBER(order))*COMB_SIZE+1);
			FOREACH (mGeneDef_it, Xm_gd, gene) {
				wsUint	N_sz		= (wsUint)gene->second.size();
				vInt&	Xv_idx		= gene->second;
				wsUint*	Na_idxMkr	= NULL;
				cComb	C_comb(N_order, N_sz);
				wsAlloc(Na_idxMkr, wsUint, N_sz);

				/* Build */
				wsUint I=0;
				FOREACHDO (vInt_it, Xv_idx, idx, I++)
					Na_idxMkr[I] = *idx;

				/* Generate combination */
				wsUint N_genComb = 0;
				do {
					N_genComb = C_comb.get(COMB_SIZE, Np_comb+1, OPT_NUMBER(order)+RES_SIZE);

					/* Do MDR */
					xMdr X_mdr = {
						X_mdrType,
						(bool)IS_ASSIGNED(cv),
						N_order,
						N_anaSamp,
						Na_data,
						Ra_anaY,
						Na_case,
						Na_ctrl,
						Na_idxMkr,
						Ra_sc, Ra_scca, Ra_scct,
						Na_testIdxCase, Na_testIdxCtrl,
						Ra_permY
					};
					_mdr(&X_mdr, N_cell, N_genComb, Np_comb);

					/* Update result */
					Cp_res[cv]->update((float *)(Np_comb+1), N_genComb);
				} while (N_genComb == COMB_SIZE);

				/* Get highest...? */
				/* FIXME : Implementation required */
				/* Export */
				float*		Ra_res	= Cp_res[cv]->getResult();
				wsUint		N_res	= Cp_res[cv]->getResultSize();
				vVariant&	Xa_mkr	= Cp_IO->getVariant();
				wsStrCst		Sa_labels[RES_SIZE] = { NULL, };
				if (IS_ASSIGNED(nperm))
					Sa_labels[0] = "P_MDR";
				else
					Sa_labels[0] = "BAL_ACC";

/**/			cExporter* Cp_mdr = NULL;
				if (N_cxval > 1) {
					char S_ext[512];
					sprintf(S_ext, "mdr.cv%d.res", cv);
					Cp_mdr = OPT_ENABLED(verbose) ? cExporter::summon(S_ext) : NULL;
				} else
					Cp_mdr = cExporter::summon("mdr.res");

				if (Cp_mdr) {
					Cp_mdr->put("SNP1");
					for (wsUint i = 1; i < N_order; i++)
						Cp_mdr->fmt("	SNP%d", i + 1);
					for (wsUint i = 0; i < RES_SIZE; i++)
						Cp_mdr->fmt("	%s", Sa_labels[i]);
					Cp_mdr->put("\n");

					for (wsUint i = 0; i < N_res; i++, Ra_res += N_order + RES_SIZE) {
						if (N_cxval == 1 && IS_ASSIGNED(pvalrange) &&
							!isInRange(OPT_RANGE(pvalrange), Ra_res[N_order + 0]))
							continue;
						/* 0~[order-1] th element represents the index of curr. comb. */
						Cp_mdr->fmt("%s", Xa_mkr[*(wsUint *)(Ra_res + 0)].name);
						for (wsUint j = 1; j < N_order; j++)
							Cp_mdr->fmt("	%s", Xa_mkr[*(wsUint *)(Ra_res + j)].name);
						for (wsUint j = 0; j < RES_SIZE; j++)
							Cp_mdr->fmt("	%g", Ra_res[N_order + j]);
						Cp_mdr->put("\n");
					}
					delete Cp_mdr;
				}
			}
		} else {
			/* 1st element of data for each thread : Number of combs. assigned
				* 2nd~(1+(order+1)*(1st elem.))th element of data for each thread
				*  : Combinations[0~(order-1)],result[order~(order+RES_SIZE-1)] */
			/* Example of order 2, RES_SIZE=3, COMB_SIZE=2
			 * Ncomb   Order1 Res1          Order2 Res2
			 * [1]   [ [2][3] [4][5][6] ] [ [7][8] [9][10][11] ]
			 *       comb1                comb2
			 */
			xThreadMdr X_thr = { cv, Na_testIdxCase, Na_testIdxCtrl, this };
			WORKER().run(thr_mdr, dist_mdr, &X_thr, NULL,
				sizeof(wsUint)*((RES_SIZE+OPT_NUMBER(order))*COMB_SIZE+1));

			/* Export */
			float*		Ra_res	= Cp_res[cv]->getResult();
			wsUint		N_res	= Cp_res[cv]->getResultSize();
			vVariant&	Xa_mkr	= Cp_IO->getVariant();
			wsStrCst		Sa_labels[RES_SIZE] = { NULL, };
			Sa_labels[0] = IS_ASSIGNED(nperm) ? "P_MDR" : "BAL_ACC";
/**/		cExporter*	Cp_mdr	= NULL;
			if (N_cxval > 1) {
				char S_ext[512];
				sprintf(S_ext, "mdr.cv%d.res", cv);
				Cp_mdr = cExporter::summon(S_ext);
			} else
				Cp_mdr = cExporter::summon("mdr.res");
			Cp_mdr->put("SNP1");
			for (wsUint i=1 ; i<N_order ; i++)	Cp_mdr->fmt("	SNP%d", i+1);
			for (wsUint i=0 ; i<RES_SIZE ; i++)	Cp_mdr->fmt("	%s", Sa_labels[i]);
			Cp_mdr->put("\n");

			for (wsUint i=0 ; i<N_res ; i++,Ra_res+=N_order+RES_SIZE) {
				wsReal R_mea = Ra_res[N_order + 0];
				if (IS_ASSIGNED(mdrthr) && !isInRange(OPT_RANGE(mdrthr), R_mea))
					continue;
				/* 0~[order-1] th element represents the index of curr. comb. */
				Cp_mdr->fmt("%s", Xa_mkr[*(wsUint *)(Ra_res+0)].name);
				for (wsUint j=1 ; j<N_order ; j++)
					Cp_mdr->fmt("	%s", Xa_mkr[*(wsUint *)(Ra_res+j)].name);
				for (wsUint j=0 ; j<RES_SIZE ; j++)
					Cp_mdr->fmt("	%g", Ra_res[N_order + j]);
				Cp_mdr->put("\n");
			}
			delete Cp_mdr;
		}
	} else {
		cStrFile	C_hp(OPT_STRING(hmdrprior), "HMDR prior interaction pairs");
		mDataIdx	Xm_varmap;

		/* Build variant map */
		wsUint		i		= 0;
		FOREACHDO (vVariant_it, Xv_var, I, i++)
			Xm_varmap.insert(make_pair(I->name, i));

		/* Prepare buffer */
		char*	Sp_buf	= NULL;
		wsAlloc(Sp_buf, char, 4096);

		/* Count the lines */
		wsUint	N_line	= 0;
		while (C_hp.gets(Sp_buf, 4096)) N_line++;
		C_hp.rewind();

		/* Read lines */
		wsUint*		Na_res	= NULL;
		wsUint		j		= 0;
		wsUint		k		= 0;
		mDataIdx	Xm_exist;
		wsAlloc(Na_res, wsUint, N_line * (N_order + 1));
		N_line	= 0;
		for (i=0 ; C_hp.gets(Sp_buf, 4096) ; i++) {
			char *a = NULL, *b = NULL;
			k = N_line * (N_order + 1);

			for (a=Sp_buf,j=0 ; a && j<N_order ; a=b,j++) {
				getString(&a, &b); // a = snp#k
				if (!a || !a[0]) halt("Line [%d] in file [%s] have invalid "
					"data", i+1, OPT_STRING(hmdrprior));
				mDataIdx_it X_find = Xm_varmap.find(a);
				if (X_find == Xm_varmap.end()) {
					if (!Xm_exist[a]) {
						LOGwarn("Failed to find variant [%s] in the final dataset\n",
							a);
						Xm_exist[a] = 1;
					}
					break;
				} else
					Na_res[k+j] = (wsUint)(X_find->second);
			}
			if (j == N_order) {
				wsReal R_mea = atof(b);

				if (IS_ASSIGNED(mdrthr) && !isInRange(OPT_RANGE(mdrthr), R_mea))
					continue;

				/* Fill stat */
				((float *)Na_res)[k+N_order] = (float)R_mea;
				N_line++;
			}
		}
		/* Update result */
		Cp_res[0]->update((float *)Na_res, N_line);
		LOG("[%d] of HMDR prior results loaded and updated\n", N_line);

		/* Uninit */
		C_hp.close();
		DEALLOC(Na_res);
		DEALLOC(Sp_buf);
	}

	/* If --hmdr, create initial set */
	if (IS_ASSIGNED(hmdr)) {
		float*		Ra_res = Cp_res[0]->getResult();
		vHmdrEntry	Xv_network;
		wsStrCst		Sa_labels[RES_SIZE] = { NULL, };
		if (IS_ASSIGNED(nperm))
			Sa_labels[0] = "P_MDR";
		else
			Sa_labels[0] = "BAL_ACC";
/**/	cExporter* Cp_hmdr = cExporter::summon("hmdr.res");
		Cp_hmdr->put("COMBIN");
		for (wsUint i=0 ; i<RES_SIZE ; i++)
			Cp_hmdr->fmt("	%s", Sa_labels[i]);
		Cp_hmdr->put("\n");

		for (wsUint i=0,j=0 ; i<Cp_res[0]->getResultSize() ; i++,j+=N_order+RES_SIZE) {
			xHmdrEntry X;
			/* Set combination */
			LOOP (k, N_order)
				X.Nv_comb.push_back(*((wsUint *)(Ra_res+j+k)));
			/* Insert last result measure as determinant measure */
			X.R_mea = Ra_res[j+N_order+RES_SIZE-1];

			Xv_network.push_back(X);
		}

		/* Based on initial graph, do hierarchical MDR */
		/* Do MDR */
		wsUint	N_maxO = 0;
		vInt	tmp;
		wsUint	N_targetO = OPT_NUMBER(hmdr);
		xMdr	X_mdr = {
			X_mdrType,
			(bool)IS_ASSIGNED(cv),
			N_order,
			N_anaSamp,
			Na_data,
			Ra_anaY,
			NULL, NULL, NULL,
			Ra_sc, NULL, NULL,
			tmp, tmp,
			Ra_permY
		};
		vReal	Rv_minP(20, W0);
		do {
			_hmdr2(&X_mdr, Xv_network, Rv_minP, Na_data, Cp_IO->sizeVariant(),
				Cp_hmdr, getIO()->getVariant());
			if (Xv_network.size() == 0) {
				LOGwarn("Failed to find better interaction scores, cannot perform HMDR\n");
				return;
			}
			FOREACH (vHmdrEntry_it, Xv_network, i) {
				wsUint N_curSz = (wsUint)(i->Nv_comb.size());
				if (N_curSz > N_maxO) N_maxO = N_curSz;
			}
		} while (N_maxO < N_targetO);
		delete Cp_hmdr;
	}
	sseFree(Ra_anaY);

	/* If --cv > 1, do summarization */
	if (N_cxval > 1) {
		mvReal Sa_vx;

		/* Count combinations */
		wsUint N_selMeasue = 0;
		LOOP (cv, N_cxval) {
			float*	Ra_res	= Cp_res[cv]->getResult();
			wsUint	N_res	= Cp_res[cv]->getResultSize();
			for (wsUint i = 0; i < N_res; i++, Ra_res += N_order + RES_SIZE) {
				wsReal R_mea = Ra_res[N_order + 0];
				if (IS_ASSIGNED(mdrthr) && !isInRange(OPT_RANGE(mdrthr), R_mea))
					continue;
				/* Make comb name */
				char Sa_name[512], *Sp=Sa_name;
				LOOP(j, N_order) Sp += sprintf(Sp, "%d ", *(wsUint*)(Ra_res + j));
				Sa_vx[Sa_name].push_back(Ra_res[N_order + N_selMeasue]);
			}
		}
		/* Divide by CVC */
		wsUint N_maxCVC = 1;
		FOREACH(mvReal_it, Sa_vx, it)
			if (N_maxCVC < (wsUint)it->second.size()) N_maxCVC = (wsUint)it->second.size();

		/* Reconstruct result */
		vector<vector<xMdrCvcResult> > Xv_resByCvc;
		Xv_resByCvc.resize(N_maxCVC);
		FOREACH(mvReal_it, Sa_vx, it) {
//			if (it->second.size() < 2) continue;
			/* Compute mean */
			wsReal R_sum = W0;
			wsReal R_min = W1;
			wsReal R_max = W0;
			FOREACH(vReal_it, it->second, iit) {
				wsReal R_vv = *iit;
				if (R_min > R_vv) R_min = R_vv;
				if (R_max < R_vv) R_max = R_vv;
				R_sum += R_vv;
			}
			R_sum /= (wsReal)it->second.size();
			vMdrCvcResult& Xv_cvc = Xv_resByCvc[it->second.size() - 1];
			/* Find stored index */
			vMdrCvcResult_it iit = Xv_cvc.begin();
			for (; iit != Xv_cvc.end(); iit++)
				if (iit->R_mean < R_sum) break;
			xMdrCvcResult X_res = {
				it->first,
				R_sum,
				R_min,
				R_max
			};
			Xv_cvc.insert(iit, X_res);
		}


		/* Export result */
		cExporter* Cp_mdr = cExporter::summon("mdr.cv.res");
		Cp_mdr->put("SNP1");
		for (wsUint i=1 ; i<N_order ; i++)	Cp_mdr->fmt("	SNP%d", i+1);
		Cp_mdr->put(" CVC");
		Cp_mdr->put(" MIN_BA	MAX_BA	MEAN_BA");
		Cp_mdr->put("\n");
		RFOREACH(vvMdrCvcResult_rit, Xv_resByCvc, it) {
			FOREACH(vMdrCvcResult_it, *it, iit) {
				char* S_comb = strdup(iit->S_comb.c_str());
				char* S_tok = strtok(S_comb, " ");
				Cp_mdr->fmt("%s", Xv_var[atoi(S_tok)].name);
				while (S_tok) {
					S_tok = strtok(NULL, " ");
					if (!S_tok) break;
					Cp_mdr->fmt("	%s", Xv_var[atoi(S_tok)].name);
				}
				free(S_comb);
				Cp_mdr->fmt("	%d	%g	%g	%g\n", N_maxCVC,
					iit->R_min, iit->R_max, iit->R_mean);
			}
			N_maxCVC--;
		}
		delete Cp_mdr;
	}
}

bool sort_xHmdrEntry(const xHmdrEntry& a, const xHmdrEntry& b) {
	return a.R_mea > b.R_mea;
}

typedef map<wsUint,int> mInt;
typedef mInt::iterator	mInt_it;

void _hmdr(xMdr* Xp_mdr, vHmdrEntry& Xv_he, vReal& Rv_minP, char** Na_data,
	wsUint N_vrt, cExporter *Cp_hmdr, vVariant& Xv_vrt)
{
	/* FIXME : HMDR supports GMDR */
//	VEC_c		Ra_sc	= Xp_mdr->Ra_sc;
	wsUint		N_entry	= (wsUint)Xv_he.size();
	vHmdrEntry	Xv_newHE;
	if (N_entry == 0)
		halt("SYSERR : No entries were left!");
	wsUint		Na_thrHmdrAll[3]	= { 1, 1, 1 };
	wsUint		N_rr		= 0;
	wsUint*		Na_pp		= loadIntValues(OPT_STRING(hmdrall), &N_rr);
	LOOP (i, N_rr)
		Na_thrHmdrAll[i] = Na_pp[i];

	float		R_minV	= 999.0f;
	float		R_maxV	= -999.0f;
	wsUint		N_maxO	= 0;
	wsUint		N_flmP	= 0;
	vReal		Rv_curMinP(20, W1);
	vInt		tmp;
	LOOP (i, N_entry-1) {
		vHmdrEntry	Xv_curH1, Xv_curH2, Xv_curH3;

		if ((i%10) == 0) notice("[%d/%d] checked, size [%d], "
			"BA [%g ~ %g], max. order [%d]...\r", i, N_entry, Xv_newHE.size(),
			R_minV, R_maxV, N_maxO);

		xHmdrEntry&	I = Xv_he[i];
		xHmdrEntry	X_sing;
		if (Na_thrHmdrAll[0]) for (wsUint j=i+1 ; j<N_entry ; j++) {
			xHmdrEntry&	J = Xv_he[j];
			xHmdrEntry	X_new;

			/* Generate new combination */
			mInt Xm_elem;
			FOREACH (vInt_it, I.Nv_comb, II) Xm_elem[*II] = 1;
			FOREACH (vInt_it, J.Nv_comb, JJ) Xm_elem[*JJ] = 1;
			FOREACH (mInt_it, Xm_elem, KK) X_new.Nv_comb.push_back(KK->first);
			std::sort(X_new.Nv_comb.begin(), X_new.Nv_comb.end());

			wsUint	N_curO	= (wsUint)Xm_elem.size();
			wsUint	N_cell	= (wsUint)pow(3.0, N_curO);	///< Number of possible cells for given order of combination
			// Current combination of markers is stored in Np_data ~ Np_data+order-1
			wsUint* Np_data = NULL;
			wsAlloc(Np_data, wsUint, N_curO+RES_SIZE);
			for (wsUint k=0 ; k<(wsUint)X_new.Nv_comb.size() ; k++)
				Np_data[k] = X_new.Nv_comb[k];

			wsUmat	Na_cntCase	= sseUmatrix(2, N_cell); /* Number of cases for each cell (0)train (1)test */
			wsUmat	Na_cntCtrl	= sseUmatrix(2, N_cell); /* Number of ctrls for each cell (0)train (1)test */
			wsVec	Ra_scca		= NULL;
			wsVec	Ra_scct		= NULL;
			if (OPT_ENABLED(gmdr)) {
				Ra_scct = sseVector(N_cell);
				Ra_scca = sseVector(N_cell);
			}

			Xp_mdr->Na_cntCase	= Na_cntCase;
			Xp_mdr->Na_cntCtrl	= Na_cntCtrl;
			Xp_mdr->N_order		= N_curO;
			_mdrY(Xp_mdr, N_cell, Xp_mdr->Ra_Y, Na_data, Np_data,
				Xp_mdr->Nv_testIdxCase.begin(), Xp_mdr->Nv_testIdxCtrl.begin(),
				tmp, tmp, tmp);

			// Classify high/low
			X_new.R_mea = _measure(Xp_mdr, N_cell, 0);
			if (X_new.R_mea > W1) halt("ERROR");

			sseUnmat(Na_cntCase, 2);
			sseUnmat(Na_cntCtrl, 2);
			sseFree(Ra_scca);
			sseFree(Ra_scct);

			// If higher achieved
			if (X_new.R_mea > I.R_mea) {
				if (IS_ASSIGNED(hmdrall)) {
					if (R_minV > X_new.R_mea) R_minV = X_new.R_mea;
					if (R_maxV < X_new.R_mea) R_maxV = X_new.R_mea;
					if (X_new.Nv_comb.size() > N_maxO)
						N_maxO = (wsUint)X_new.Nv_comb.size();

					/* Store p-values only exceeding curMinP */
					if (Rv_minP[X_new.Nv_comb.size()] < X_new.R_mea) {
						/* Update curMinP */
						if (Rv_curMinP[X_new.Nv_comb.size()] > X_new.R_mea)
							Rv_curMinP[X_new.Nv_comb.size()] = X_new.R_mea;
						Xv_curH1.push_back(X_new);
					} else N_flmP++;
				} else {
					X_sing.Nv_comb.resize(X_new.Nv_comb.size());
					copy(X_new.Nv_comb.begin(), X_new.Nv_comb.end(), X_sing.Nv_comb.begin());
					X_sing.R_mea = X_new.R_mea;
				}
			}
//				Xv_newHE.push_back(X_new);
		}
		if (X_sing.Nv_comb.size() && !IS_ASSIGNED(hmdrall)) {
			if (X_sing.R_mea > W1) halt("ERROR");
			if (R_minV > X_sing.R_mea) R_minV = X_sing.R_mea;
			if (R_maxV < X_sing.R_mea) R_maxV = X_sing.R_mea;
			if (X_sing.Nv_comb.size() > N_maxO)
				N_maxO = (wsUint)X_sing.Nv_comb.size();
			Xv_newHE.push_back(X_sing);
		}
		wsUint	N_order = (wsUint)I.Nv_comb.size() + 1;
		wsUint	N_cell	= (wsUint)pow(3.0, N_order);
		wsUint*	Np_data	= NULL;
		wsAlloc(Np_data, wsUint, N_order+RES_SIZE);
		mInt	Xm_pass;
		LOOP (j, N_order-1) {
			wsUint N_idx = I.Nv_comb[j];
			Xm_pass[N_idx] = 1;
			Np_data[j] = N_idx;
		}
		wsUmat	Na_cntCase	= sseUmatrix(2, N_cell); /* Number of cases for each cell */
		wsUmat	Na_cntCtrl	= sseUmatrix(2, N_cell); /* Number of ctrls for each cell */
		wsVec	Ra_scca		= NULL;
		wsVec	Ra_scct		= NULL;
		if (OPT_ENABLED(gmdr)) {
			Ra_scct = sseVector(N_cell);
			Ra_scca = sseVector(N_cell);
		}
		Xp_mdr->Na_cntCase	= Na_cntCase;
		Xp_mdr->Na_cntCtrl	= Na_cntCtrl;
		Xp_mdr->N_order		= N_order;

		xHmdrEntry X_sing2;
		X_sing2.R_mea = (float)W0;
		if (Na_thrHmdrAll[1]) LOOP (j, N_vrt) {
			xHmdrEntry X_new;
			/* If this variant is already in the set, just skip */
			if (Xm_pass[j]) continue;

			/* Generate combination */
			FOREACH (vInt_it, I.Nv_comb, KK) X_new.Nv_comb.push_back(*KK);
			X_new.Nv_comb.push_back(j);

			/* Fill combination for MDR */
			Np_data[N_order-1] = j;

			/* Perform MDR itself */
			_mdrY(Xp_mdr, N_cell, Xp_mdr->Ra_Y, Na_data, Np_data, 
				Xp_mdr->Nv_testIdxCase.begin(), Xp_mdr->Nv_testIdxCtrl.begin(),
				tmp, tmp, tmp);

			/* Classify high/low */
			X_new.R_mea = _measure(Xp_mdr, N_cell, 0);
			if (X_new.R_mea > W1) halt("ERROR");

			/* If higher achieved */
			if (X_new.R_mea > I.R_mea) {
				if (IS_ASSIGNED(hmdrall)) {
					if (R_minV > X_new.R_mea) R_minV = X_new.R_mea;
					if (R_maxV < X_new.R_mea) R_maxV = X_new.R_mea;
					if (X_new.Nv_comb.size() > N_maxO)
						N_maxO = (wsUint)X_new.Nv_comb.size();

					/* Store p-values only exceeding curMinP */
					if (Rv_minP[X_new.Nv_comb.size()] < X_new.R_mea) {
						/* Update curMinP */
						if (Rv_curMinP[X_new.Nv_comb.size()] > X_new.R_mea)
							Rv_curMinP[X_new.Nv_comb.size()] = X_new.R_mea;
						Xv_curH2.push_back(X_new);
					} else N_flmP++;
				} else {
					X_sing2.Nv_comb.resize(X_new.Nv_comb.size());
					copy(X_new.Nv_comb.begin(), X_new.Nv_comb.end(), X_sing2.Nv_comb.begin());
					X_sing2.R_mea = X_new.R_mea;
				}
			}
		}
		if (X_sing2.Nv_comb.size() && !IS_ASSIGNED(hmdrall)) {
			if (X_sing2.R_mea > W1) halt("ERROR");
			if (R_minV > X_sing2.R_mea) R_minV = X_sing2.R_mea;
			if (R_maxV < X_sing2.R_mea) R_maxV = X_sing2.R_mea;
			if (X_sing2.Nv_comb.size() > N_maxO)
				N_maxO = (wsUint)X_sing2.Nv_comb.size();
			Xv_newHE.push_back(X_sing2);
		}

		xHmdrEntry	X_sing3;
		wsUint		N_szOri	= (wsUint)I.Nv_comb.size();
		X_sing3.R_mea		= (float)W0;
		Xp_mdr->N_order		= N_szOri - 1;
		if (Na_thrHmdrAll[2]) if (OPT_NUMBER(order) < (int)N_szOri) {
			wsUint	N_cell	= (wsUint)pow(3.0, Xp_mdr->N_order);
			LOOP (j, N_szOri) {
				xHmdrEntry X_new;

				/* Generate combination */
				LOOP (k, j)
					X_new.Nv_comb.push_back(I.Nv_comb[k]);
				for (wsUint k=j+1 ; k<N_szOri ; k++)
					X_new.Nv_comb.push_back(I.Nv_comb[k]);

				/* Perform MDR itself */
				_mdrY(Xp_mdr, N_cell, Xp_mdr->Ra_Y, Na_data, Np_data,
					Xp_mdr->Nv_testIdxCase.begin(), Xp_mdr->Nv_testIdxCtrl.begin(),
					tmp, tmp, tmp);

				/* Classify high/low */
				X_new.R_mea = _measure(Xp_mdr, N_cell, 0);
				if (X_new.R_mea > W1) halt("ERROR");

				/* If higher achieved */
				if (X_new.R_mea <= I.R_mea) continue;
				if (IS_ASSIGNED(hmdrall)) {
					if (R_minV > X_new.R_mea) R_minV = X_new.R_mea;
					if (R_maxV < X_new.R_mea) R_maxV = X_new.R_mea;
					if (X_new.Nv_comb.size() > N_maxO)
						N_maxO = (wsUint)X_new.Nv_comb.size();

					/* Store p-values only exceeding curMinP */
					if (Rv_minP[X_new.Nv_comb.size()] < X_new.R_mea) {
						/* Update curMinP */
						if (Rv_curMinP[X_new.Nv_comb.size()] > X_new.R_mea)
							Rv_curMinP[X_new.Nv_comb.size()] = X_new.R_mea;
						Xv_curH3.push_back(X_new);
					} else N_flmP++;
				} else {
					X_sing3.Nv_comb.resize(X_new.Nv_comb.size());
					copy(X_new.Nv_comb.begin(), X_new.Nv_comb.end(), X_sing3.Nv_comb.begin());
					X_sing3.R_mea = X_new.R_mea;
				}
			}
		}
		if (X_sing3.Nv_comb.size() && !IS_ASSIGNED(hmdrall)) {
			if (X_sing3.R_mea > W1) halt("ERROR");
			if (R_minV > X_sing3.R_mea) R_minV = X_sing3.R_mea;
			if (R_maxV < X_sing3.R_mea) R_maxV = X_sing3.R_mea;
			if (X_sing3.Nv_comb.size() > N_maxO)
				N_maxO = (wsUint)X_sing3.Nv_comb.size();
			Xv_newHE.push_back(X_sing3);
		}
		/* Sort and insert */
		if (Xv_curH1.size()) {
			sort(Xv_curH1.begin(), Xv_curH1.end(), sort_xHmdrEntry);
			sort(Xv_curH2.begin(), Xv_curH2.end(), sort_xHmdrEntry);
			sort(Xv_curH3.begin(), Xv_curH3.end(), sort_xHmdrEntry);
			for (wsUint i=0 ; i<Na_thrHmdrAll[0] ; i++)
				if (Xv_curH1.size() > i) Xv_newHE.push_back(Xv_curH1[i]);
			for (wsUint i=0 ; i<Na_thrHmdrAll[1] ; i++)
				if (Xv_curH2.size() > i) Xv_newHE.push_back(Xv_curH2[i]);
			for (wsUint i=0 ; i<Na_thrHmdrAll[2] ; i++)
				if (Xv_curH3.size() > i) Xv_newHE.push_back(Xv_curH3[i]);
		}

		sseUnmat(Na_cntCase, 2);
		sseUnmat(Na_cntCtrl, 2);
		sseFree(Ra_scca);
		sseFree(Ra_scct);
	}
	/* Update minP */
	LOOP (j, 20)
		if (Rv_curMinP[j] != W1)
			Rv_minP[j] = Rv_curMinP[j];

	if (0) { //OPT_ENABLED(hmdrall)) {
		wsUint			N_del	= 0;
		vector<vReal>	Xv_vals;
		vReal			Xv_mean;

		Xv_vals.resize(20);
		Xv_mean.resize(20);

		FOREACH (vHmdrEntry_it, Xv_newHE, i)
			Xv_vals[i->Nv_comb.size()].push_back(i->R_mea);

		LOOP (i, 20) {
			if (Xv_vals[i].size())
				Xv_mean[i] = accumulate(Xv_vals[i].begin(), Xv_vals[i].end(), 0.0) / (wsReal)Xv_vals[i].size();
		}

		for (vHmdrEntry_it i=Xv_newHE.begin() ; i!=Xv_newHE.end() ; ) {
			if (Xv_mean[i->Nv_comb.size()] > i->R_mea) {
				i = Xv_newHE.erase(i);
				N_del++;
				continue;
			}
			i++;
		}
		LOGnote("[%d] entries removed\n", N_del);
	}
	LOG("[%d/%d] checked, size [%d], BA [%g ~ %g], max. order [%d]\n", N_entry,
		N_entry, Xv_newHE.size(), R_minV, R_maxV, N_maxO);

	/* Export if available */
	if (Cp_hmdr) {
		FOREACH (vHmdrEntry_it, Xv_newHE, i) {
			Cp_hmdr->put(Xv_vrt[i->Nv_comb[0]].name);
			for (wsUint j=1 ; j<i->Nv_comb.size() ; j++)
				Cp_hmdr->fmt(",%s", Xv_vrt[i->Nv_comb[j]].name);
			Cp_hmdr->fmt("	%g\n", i->R_mea);
		}
	}


	LOOP (i, 20)
		if (Rv_minP[i] != W0)
			pverbose("minP for order [%d] is [%g]\n", i, Rv_minP[i]);
	pverbose("[%d] entries were filtered by minP\n", N_flmP);

	/* Replace buffer */
	Xv_he.resize(Xv_newHE.size());
	copy(Xv_newHE.begin(), Xv_newHE.end(), Xv_he.begin());
}

void _hmdr2(xMdr* Xp_mdr, vHmdrEntry& Xv_he, vReal& Rv_minP, char** Na_data,
	wsUint N_vrt, cExporter *Cp_hmdr, vVariant& Xv_vrt)
{	
	/* FIXME : HMDR supports GMDR */
//	VEC_c		Ra_sc	= Xp_mdr->Ra_sc;
	wsUint		N_entry	= (wsUint)Xv_he.size();
	vHmdrEntry	Xv_newHE;
	if (N_entry == 0)
		halt("SYSERR : No entries were left!");
	wsUint		Na_thrHmdrAll[3]	= { 1, 1, 1 };
	wsUint		N_rr		= 0;
	wsUint*		Na_pp		= loadIntValues(OPT_STRING(hmdrall), &N_rr);
	LOOP (i, N_rr)
		Na_thrHmdrAll[i] = Na_pp[i];

	/* Draw combination map */
	eStr Xm_isExist;
	FOREACH (vHmdrEntry_it, Xv_he, i) {
		char	S_buf[1024], *p = S_buf;
		FOREACH (vInt_it, i->Nv_comb, j)
			p += sprintf(p, "%d,", *j);
		Xm_isExist.insert(S_buf);
	}

	float		R_minV	= 999.0f;
	float		R_maxV	= -999.0f;
	wsUint		N_maxO	= 0;
	wsUint		N_flmP	= 0;
	vReal		Rv_curMinP(20, W1);
	vInt		tmp;
	LOOP (i, N_entry-1) {
		vHmdrEntry	Xv_curH1, Xv_curH2, Xv_curH3;

		if ((i%10) == 0) notice("[%d/%d] checked, size [%d], "
			"BA [%g ~ %g], max. order [%d]...\r", i, N_entry, Xv_newHE.size(),
			R_minV, R_maxV, N_maxO);

		xHmdrEntry&	I = Xv_he[i];
		xHmdrEntry	X_sing;

		mInt		Xm_elem;
		FOREACH (vInt_it, I.Nv_comb, II) Xm_elem[*II] = 1;

		if (Na_thrHmdrAll[0]) for (wsUint j=i+1 ; j<N_entry ; j++) {
			xHmdrEntry&	J = Xv_he[j];
			xHmdrEntry	X_new;
			X_new.Nv_cellIdx.resize(Xp_mdr->N_sample);

			/* Generate additional combination */
			FOREACH (vInt_it, J.Nv_comb, KK)
				if (Xm_elem.find(*KK) == Xm_elem.end()) X_new.Nv_comb.push_back(*KK);
			/* If no addition is made, just skip */
			if (X_new.Nv_comb.size() == 0) continue;
			//std::sort(X_new.Nv_comb.begin(), X_new.Nv_comb.end());

			wsUint	N_curO	= (wsUint)(Xm_elem.size() + X_new.Nv_comb.size());
			wsUint	N_cell	= (wsUint)pow(3.0, N_curO);	///< Number of possible cells for given order of combination

			wsUmat	Na_cntCase	= sseUmatrix(2, N_cell); /* Number of cases for each cell (0)train (1)test */
			wsUmat	Na_cntCtrl	= sseUmatrix(2, N_cell); /* Number of ctrls for each cell (0)train (1)test */
			wsVec	Ra_scca		= NULL;
			wsVec	Ra_scct		= NULL;
			if (OPT_ENABLED(gmdr)) {
				Ra_scct = sseVector(N_cell);
				Ra_scca = sseVector(N_cell);
			}

			Xp_mdr->Na_cntCase	= Na_cntCase;
			Xp_mdr->Na_cntCtrl	= Na_cntCtrl;
			Xp_mdr->N_order		= N_curO;
			_mdrY(Xp_mdr, N_cell, Xp_mdr->Ra_Y, Na_data, I.Nv_comb,
				Xp_mdr->Nv_testIdxCase.begin(), Xp_mdr->Nv_testIdxCtrl.begin(),
				X_new.Nv_comb, I.Nv_cellIdx, X_new.Nv_cellIdx);

			// Classify high/low
			X_new.R_mea = _measure(Xp_mdr, N_cell, 0);
			if (X_new.R_mea > W1) halt("ERROR");

// 			vInt tmp;
// 			_mdrY(Xp_mdr, N_cell, Xp_mdr->Ra_Y, Na_data, X_new.Nv_comb,
// 				Xp_mdr->Nv_testIdxCase.begin(), Xp_mdr->Nv_testIdxCtrl.begin(),
// 				tmp, tmp, tmp);
// 			float U = _measure(Xp_mdr, N_cell, 0);
// 			if (U != X_new.R_mea) halt("ERROR2");

			sseUnmat(Na_cntCase, 2);
			sseUnmat(Na_cntCtrl, 2);
			sseFree(Ra_scca);
			sseFree(Ra_scct);

			// If higher achieved
			if (X_new.R_mea > I.R_mea) {
				// Fill combination
				std::sort(X_new.Nv_comb.begin(), X_new.Nv_comb.end());

				if (IS_ASSIGNED(hmdrall)) {
					if (R_minV > X_new.R_mea) R_minV = X_new.R_mea;
					if (R_maxV < X_new.R_mea) R_maxV = X_new.R_mea;
					if (X_new.Nv_comb.size() > N_maxO)
						N_maxO = (wsUint)X_new.Nv_comb.size();

					/* Store p-values only exceeding curMinP */
					if (Rv_minP[X_new.Nv_comb.size()] < X_new.R_mea) {
						/* Update curMinP */
						if (Rv_curMinP[X_new.Nv_comb.size()] > X_new.R_mea)
							Rv_curMinP[X_new.Nv_comb.size()] = X_new.R_mea;
						Xv_curH1.push_back(X_new);
					} else N_flmP++;
				} else {
					X_sing.Nv_comb.resize(X_new.Nv_comb.size());
					copy(X_new.Nv_comb.begin(), X_new.Nv_comb.end(), X_sing.Nv_comb.begin());
					X_sing.R_mea = X_new.R_mea;

					X_sing.Nv_cellIdx.resize(Xp_mdr->N_sample);
					copy(X_sing.Nv_cellIdx.begin(), X_sing.Nv_cellIdx.end(),
						X_new.Nv_cellIdx.begin());
				}
			}
			//				Xv_newHE.push_back(X_new);
		}
		if (X_sing.Nv_comb.size() && !IS_ASSIGNED(hmdrall)) {
			if (X_sing.R_mea > W1) halt("ERROR");
			if (R_minV > X_sing.R_mea) R_minV = X_sing.R_mea;
			if (R_maxV < X_sing.R_mea) R_maxV = X_sing.R_mea;
			if (X_sing.Nv_comb.size() > N_maxO)
				N_maxO = (wsUint)X_sing.Nv_comb.size();
			Xv_newHE.push_back(X_sing);
		}
		wsUint	N_order = (wsUint)I.Nv_comb.size() + 1;
		wsUint	N_cell	= (wsUint)pow(3.0, N_order);
		mInt	Xm_pass;
		LOOP (j, N_order-1) {
			wsUint N_idx = I.Nv_comb[j];
			Xm_pass[N_idx] = 1;
		}
		wsUmat	Na_cntCase	= sseUmatrix(2, N_cell); /* Number of cases for each cell */
		wsUmat	Na_cntCtrl	= sseUmatrix(2, N_cell); /* Number of ctrls for each cell */
		wsVec	Ra_scca		= NULL;
		wsVec	Ra_scct		= NULL;
		if (OPT_ENABLED(gmdr)) {
			Ra_scct = sseVector(N_cell);
			Ra_scca = sseVector(N_cell);
		}
		Xp_mdr->Na_cntCase	= Na_cntCase;
		Xp_mdr->Na_cntCtrl	= Na_cntCtrl;
		Xp_mdr->N_order		= N_order;

		xHmdrEntry X_sing2;
		X_sing2.R_mea = (float)W0;
		if (Na_thrHmdrAll[1]) LOOP (j, N_vrt) {
			xHmdrEntry X_new;
			/* If this variant is already in the set, just skip */
			if (Xm_pass[j]) continue;
			X_new.Nv_cellIdx.resize(Xp_mdr->N_sample);

			/* Generate combination */
			//FOREACH (vInt_it, I.Nv_comb, KK) X_new.Nv_comb.push_back(*KK);
			X_new.Nv_comb.push_back(j);

			/* Perform MDR itself */
			_mdrY(Xp_mdr, N_cell, Xp_mdr->Ra_Y, Na_data, I.Nv_comb,
				Xp_mdr->Nv_testIdxCase.begin(), Xp_mdr->Nv_testIdxCtrl.begin(),
				X_new.Nv_comb, I.Nv_cellIdx, X_new.Nv_cellIdx);

			/* Classify high/low */
			X_new.R_mea = _measure(Xp_mdr, N_cell, 0);
			if (X_new.R_mea > W1) halt("ERROR");

			/* If higher achieved */
			if (X_new.R_mea > I.R_mea) {
				// Fill combination
				std::sort(X_new.Nv_comb.begin(), X_new.Nv_comb.end());

				if (IS_ASSIGNED(hmdrall)) {
					if (R_minV > X_new.R_mea) R_minV = X_new.R_mea;
					if (R_maxV < X_new.R_mea) R_maxV = X_new.R_mea;
					if (X_new.Nv_comb.size() > N_maxO)
						N_maxO = (wsUint)X_new.Nv_comb.size();

					/* Store p-values only exceeding curMinP */
					if (Rv_minP[X_new.Nv_comb.size()] < X_new.R_mea) {
						/* Update curMinP */
						if (Rv_curMinP[X_new.Nv_comb.size()] > X_new.R_mea)
							Rv_curMinP[X_new.Nv_comb.size()] = X_new.R_mea;
						Xv_curH2.push_back(X_new);
					} else N_flmP++;
				} else {
					X_sing2.Nv_comb.resize(X_new.Nv_comb.size());
					copy(X_new.Nv_comb.begin(), X_new.Nv_comb.end(), X_sing2.Nv_comb.begin());
					X_sing2.R_mea = X_new.R_mea;

					X_sing2.Nv_cellIdx.resize(Xp_mdr->N_sample);
					copy(X_sing2.Nv_cellIdx.begin(), X_sing2.Nv_cellIdx.end(),
						X_new.Nv_cellIdx.begin());
				}
			}
		}
		if (X_sing2.Nv_comb.size() && !IS_ASSIGNED(hmdrall)) {
			if (X_sing2.R_mea > W1) halt("ERROR");
			if (R_minV > X_sing2.R_mea) R_minV = X_sing2.R_mea;
			if (R_maxV < X_sing2.R_mea) R_maxV = X_sing2.R_mea;
			if (X_sing2.Nv_comb.size() > N_maxO)
				N_maxO = (wsUint)X_sing2.Nv_comb.size();
			Xv_newHE.push_back(X_sing2);
		}

		xHmdrEntry	X_sing3;
		wsUint		N_szOri	= (wsUint)I.Nv_comb.size();
		X_sing3.R_mea		= (float)W0;
		Xp_mdr->N_order		= N_szOri - 1;
		if (0) if (Na_thrHmdrAll[2]) if (OPT_NUMBER(order) < (int)N_szOri) {
			wsUint	N_cell	= (wsUint)pow(3.0, Xp_mdr->N_order);
			LOOP (j, N_szOri) {
				xHmdrEntry X_new;
				X_new.Nv_cellIdx.resize(Xp_mdr->N_sample);

				/* Generate combination */
				LOOP (k, j)
					X_new.Nv_comb.push_back(I.Nv_comb[k]);
				for (wsUint k=j+1 ; k<N_szOri ; k++)
					X_new.Nv_comb.push_back(I.Nv_comb[k]);

				/* Perform MDR itself */
				_mdrY(Xp_mdr, N_cell, Xp_mdr->Ra_Y, Na_data, X_new.Nv_comb,
					Xp_mdr->Nv_testIdxCase.begin(), Xp_mdr->Nv_testIdxCtrl.begin(),
					tmp, tmp, tmp);

				/* Classify high/low */
				X_new.R_mea = _measure(Xp_mdr, N_cell, 0);
				if (X_new.R_mea > W1) halt("ERROR");

				/* If higher achieved */
				if (X_new.R_mea <= I.R_mea) continue;
				if (IS_ASSIGNED(hmdrall)) {
					if (R_minV > X_new.R_mea) R_minV = X_new.R_mea;
					if (R_maxV < X_new.R_mea) R_maxV = X_new.R_mea;
					if (X_new.Nv_comb.size() > N_maxO)
						N_maxO = (wsUint)X_new.Nv_comb.size();

					/* Store p-values only exceeding curMinP */
					if (Rv_minP[X_new.Nv_comb.size()] < X_new.R_mea) {
						/* Update curMinP */
						if (Rv_curMinP[X_new.Nv_comb.size()] > X_new.R_mea)
							Rv_curMinP[X_new.Nv_comb.size()] = X_new.R_mea;
						Xv_curH3.push_back(X_new);
					} else N_flmP++;
				} else {
					X_sing3.Nv_comb.resize(X_new.Nv_comb.size());
					copy(X_new.Nv_comb.begin(), X_new.Nv_comb.end(), X_sing3.Nv_comb.begin());
					X_sing3.R_mea = X_new.R_mea;

					X_sing3.Nv_cellIdx.resize(Xp_mdr->N_sample);
					copy(X_sing3.Nv_cellIdx.begin(), X_sing3.Nv_cellIdx.end(),
						X_new.Nv_cellIdx.begin());
				}
			}
		}
		if (X_sing3.Nv_comb.size() && !IS_ASSIGNED(hmdrall)) {
			if (X_sing3.R_mea > W1) halt("ERROR");
			if (R_minV > X_sing3.R_mea) R_minV = X_sing3.R_mea;
			if (R_maxV < X_sing3.R_mea) R_maxV = X_sing3.R_mea;
			if (X_sing3.Nv_comb.size() > N_maxO)
				N_maxO = (wsUint)X_sing3.Nv_comb.size();
			Xv_newHE.push_back(X_sing3);
		}
		/* Sort and insert */
		if (Xv_curH1.size()) {
			sort(Xv_curH1.begin(), Xv_curH1.end(), sort_xHmdrEntry);
			sort(Xv_curH2.begin(), Xv_curH2.end(), sort_xHmdrEntry);
			sort(Xv_curH3.begin(), Xv_curH3.end(), sort_xHmdrEntry);
			for (wsUint i=0 ; i<Na_thrHmdrAll[0] ; i++)
				if (Xv_curH1.size() > i) {
					char	S_buf[1024], *p = S_buf;
					FOREACH (vInt_it, Xv_curH1[i].Nv_comb, j)
						p += sprintf(p, "%d,", *j);
					if (Xm_isExist.find(S_buf) == Xm_isExist.end()) {
						Xm_isExist.insert(S_buf);
						Xv_newHE.push_back(Xv_curH1[i]);
					}
				}
			for (wsUint i=0 ; i<Na_thrHmdrAll[1] ; i++)
				if (Xv_curH2.size() > i) {
					char	S_buf[1024], *p = S_buf;
					FOREACH (vInt_it, Xv_curH2[i].Nv_comb, j)
						p += sprintf(p, "%d,", *j);
					if (Xm_isExist.find(S_buf) == Xm_isExist.end()) {
						Xm_isExist.insert(S_buf);
						Xv_newHE.push_back(Xv_curH2[i]);
					}
				}
			for (wsUint i=0 ; i<Na_thrHmdrAll[2] ; i++)
				if (Xv_curH3.size() > i) {
					char	S_buf[1024], *p = S_buf;
					FOREACH (vInt_it, Xv_curH3[i].Nv_comb, j)
						p += sprintf(p, "%d,", *j);
					if (Xm_isExist.find(S_buf) == Xm_isExist.end()) {
						Xm_isExist.insert(S_buf);
						Xv_newHE.push_back(Xv_curH3[i]);
					}
				}
		}

		sseUnmat(Na_cntCase, 2);
		sseUnmat(Na_cntCtrl, 2);
		sseFree(Ra_scca);
		sseFree(Ra_scct);
	}
	/* Update minP */
	LOOP (j, 20)
		if (Rv_curMinP[j] != W1)
			Rv_minP[j] = Rv_curMinP[j];

	if (0) { //OPT_ENABLED(hmdrall)) {
		wsUint			N_del	= 0;
		vector<vReal>	Xv_vals;
		vReal			Xv_mean;

		Xv_vals.resize(20);
		Xv_mean.resize(20);

		FOREACH (vHmdrEntry_it, Xv_newHE, i)
			Xv_vals[i->Nv_comb.size()].push_back(i->R_mea);

		LOOP (i, 20) {
			if (Xv_vals[i].size())
				Xv_mean[i] = accumulate(Xv_vals[i].begin(), Xv_vals[i].end(), 0.0) / (wsReal)Xv_vals[i].size();
		}

		for (vHmdrEntry_it i=Xv_newHE.begin() ; i!=Xv_newHE.end() ; ) {
			if (Xv_mean[i->Nv_comb.size()] > i->R_mea) {
				i = Xv_newHE.erase(i);
				N_del++;
				continue;
			}
			i++;
		}
		LOGnote("[%d] entries removed\n", N_del);
	}
	LOG("[%d/%d] checked, size [%d], BA [%g ~ %g], max. order [%d]\n", N_entry,
		N_entry, Xv_newHE.size(), R_minV, R_maxV, N_maxO);

	/* Export if available */
	if (Cp_hmdr) {
		FOREACH (vHmdrEntry_it, Xv_newHE, i) {
			Cp_hmdr->put(Xv_vrt[i->Nv_comb[0]].name);
			for (wsUint j=1 ; j<i->Nv_comb.size() ; j++)
				Cp_hmdr->fmt(",%s", Xv_vrt[i->Nv_comb[j]].name);
			Cp_hmdr->fmt("	%g\n", i->R_mea);
		}
	}


	LOOP (i, 20)
		if (Rv_minP[i] != W0)
			pverbose("minP for order [%d] is [%g]\n", i, Rv_minP[i]);
	pverbose("[%d] entries were filtered by minP\n", N_flmP);

	/* Replace buffer */
	Xv_he.resize(Xv_newHE.size());
	copy(Xv_newHE.begin(), Xv_newHE.end(), Xv_he.begin());
}

cResult* cMdrAnalysis::getResult(wsUint N_idx/*=0*/)
{
	if (Cp_res == NULL)
		halt_fmt(WISARD_SYST_NULL_MDR_RESULT);
		//halt("SYSERR: Null result class");
	return Cp_res[N_idx];
}

#endif

} // End namespace ONETOOL
