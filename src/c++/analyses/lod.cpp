#include "analyses/fam.h"
#include "analyses/lod.h"
#include "analyses/setmgr.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cLodAnalysis::cLodAnalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS,
	cSetManagerAnalysis* Cp_inpAnaSM) : cAnalysis(Cp_inpIO)
{
	wsVecCst Ra_phe = Cp_IO->getPheno();
	Cp_anaFS = Cp_inpAnaFS;
	Cp_anaSM = Cp_inpAnaSM;

	if (Cp_IO->isContinuous())
		halt("Can't do LoD analysis on continuous phenotype");

	/* For all nuclear families, select case/control NFs */
	mNucFam &Xm_nf = Cp_anaFS->getNucFamData();
	FOREACH (mNucFam_it, Xm_nf, i) {
		wsUint N_pi = i->second.Xp_pat->N_idx;
		wsUint N_mi = i->second.Xp_mat->N_idx;

		/* Pass if either is missing */
		if (isMissingReal(Ra_phe[N_pi]) || isMissingReal(Ra_phe[N_mi]) ||
			fabs(Ra_phe[N_pi] - Ra_phe[N_mi]) != W1)
			continue;

		/* Copy */
		xNucFamily	X_nf = i->second;

		/* Reset and fill with nonmissings */
		char B_rebuild = 0;
		FOREACH (vSampPtr_it, X_nf.Xp_childs, j)
			if (isMissingReal(Ra_phe[(*j)->N_idx])) {
				B_rebuild = 1;
				break;
			}
		/* Rebuild if required */
		if (B_rebuild) {
			X_nf.Xp_childs.clear();
			FOREACH (vSampPtr_it, i->second.Xp_childs, j)
				if (!isMissingReal(Ra_phe[(*j)->N_idx]))
					X_nf.Xp_childs.push_back(*j);
		}
		if (X_nf.Xp_childs.size())
			Xm_curNF.insert(make_pair(i->first, X_nf));
	}
	if (Xm_curNF.size() == 0)
		LOGwarn("No available nuclear family in the dataset for LoD analysis!\n");
	else
		LOG("[%d/%d] nuclear families are included to LoD analysis\n",
			Xm_curNF.size(), Xm_nf.size());
}

typedef struct _xLoD {
	mNucFam&	Xm_nf;
	wsRealCst*		Ra_phe;
	char**		Na_data;
	wsUint		j;
	xVariant&	X_snp;
} xLoD;

typedef struct _xLoD2 {
	mNucFam&	Xm_nf;
	wsRealCst*		Ra_phe;
	char**		Na_data;
	vInt&		Xv_snp;
} xLoD2;

double f(double v, void *x)
{
	xLoD *R = (xLoD *)x;
	wsReal	R_LoD	= W0;
	wsUint	j		= R->j;
	FOREACH (mNucFam_it, R->Xm_nf, i) {
		wsUint	N_pi = i->second.Xp_pat->N_idx;
		wsUint	N_mi = i->second.Xp_mat->N_idx;

		char	N_g0 = R->Ra_phe[N_pi] == WISARD_AFFECTED ? R->Na_data[N_mi][j]
			: R->Na_data[N_pi][j];
		char	N_g1 = R->Ra_phe[N_pi] == WISARD_AFFECTED ? R->Na_data[N_pi][j]
			: R->Na_data[N_mi][j];

		/* Cannot test if both genotype is same */
		if (N_g0 == N_g1) continue;

		wsUint N_sum = 0;
		FOREACH (vSampPtr_it, i->second.Xp_childs, k) {
			char N_geno = R->Na_data[(*k)->N_idx][j];

			/* Same genotype but different label, count this */
			if ((R->Ra_phe[(*k)->N_idx] == WISARD_AFFECTED && N_g1 != N_geno) ||
				(R->Ra_phe[(*k)->N_idx] == WISARD_UNAFFECTED && N_g0 != N_geno))
				N_sum++;
		}

		/* Current family's LoD result */
		wsUint N[2] = { N_sum, (wsUint)(i->second.Xp_childs.size()) };

		wsReal q = (W1 - (wsReal)v) / W2;
		R_LoD += log10(pow(q, N[0]) * pow(q, N[1]-N[0]) / pow(0.25, N[1]));
	}

	return -R_LoD;
}

double f2(double v, void *x)
{
	xLoD2 *R = (xLoD2 *)x;

	wsReal	R_LoD	= W0;
	FOREACH (mNucFam_it, R->Xm_nf, i) {
		wsUint	N_pi = i->second.Xp_pat->N_idx;
		wsUint	N_mi = i->second.Xp_mat->N_idx;

		wsUint	N_g0 = 0, N_g1 = 0;
		wsUint	N_i1 = R->Ra_phe[N_pi] == WISARD_AFFECTED ? N_mi : N_pi;
		wsUint	N_i2 = R->Ra_phe[N_pi] == WISARD_AFFECTED ? N_pi : N_mi;
			
		FOREACH (vInt_it, R->Xv_snp, j) {
			N_g0 <<= 2;
			N_g1 <<= 2;

			N_g0 += R->Na_data[N_i1][*j];
			N_g1 += R->Na_data[N_i2][*j];
		}

		/* Cannot test if both genotype is same */
		if (N_g0 == N_g1) continue;

		wsUint N_sum = 0;
		FOREACH (vSampPtr_it, i->second.Xp_childs, k) {
			wsUint N_geno = 0;
			FOREACH (vInt_it, R->Xv_snp, j) {
				N_geno <<= 2;
				N_geno += R->Na_data[(*k)->N_idx][*j];
			}

			/* Same genotype but different label, count this */
			if ((R->Ra_phe[(*k)->N_idx] == WISARD_AFFECTED && N_g1 != N_geno) ||
				(R->Ra_phe[(*k)->N_idx] == WISARD_UNAFFECTED && N_g0 != N_geno))
				N_sum++;
		}

		/* Current family's LoD result */
		wsUint N[2] = { N_sum, (wsUint)(i->second.Xp_childs.size()) };

		wsReal q = (W1 - (wsReal)v) / W2;
		R_LoD += log10(pow(q, N[0]) * pow(q, N[1]-N[0]) / pow(0.25, N[1]));
	}

	return -R_LoD;
}

wsReal _test(xLoD &X, wsReal *Rp)
{
	double R = Brent_fmin(0, 0.5, f, (void *)(&X), 1e-5);
	wsReal q = (W1 - (wsReal)R) / W2;
	wsUint j = X.j;

	wsReal R_LoD = W0;
	FOREACH (mNucFam_it, X.Xm_nf, i) {
		wsUint	N_pi = i->second.Xp_pat->N_idx;
		wsUint	N_mi = i->second.Xp_mat->N_idx;

		char	N_g0 = X.Ra_phe[N_pi] == WISARD_AFFECTED ? X.Na_data[N_mi][j]
		: X.Na_data[N_pi][j];
		char	N_g1 = X.Ra_phe[N_pi] == WISARD_AFFECTED ? X.Na_data[N_pi][j]
		: X.Na_data[N_mi][j];

		/* Cannot test if both genotype is same */
		if (N_g0 == N_g1) continue;

		wsUint N_sum = 0;
		FOREACH (vSampPtr_it, i->second.Xp_childs, k) {
			char N_geno = X.Na_data[(*k)->N_idx][j];

			/* Same genotype but different label, count this */
			if ((X.Ra_phe[(*k)->N_idx] == WISARD_AFFECTED && N_g1 != N_geno) ||
				(X.Ra_phe[(*k)->N_idx] == WISARD_UNAFFECTED && N_g0 != N_geno))
				N_sum++;
		}

		/* Current family's LoD result */
		wsUint N[2] = { N_sum, (wsUint)(i->second.Xp_childs.size()) };
		R_LoD += log10(pow(q, N[0]) * pow(q, N[1]-N[0]) / pow(0.25, N[1]));
	}

	*Rp = R;
	return R_LoD / 2;
}

wsReal _test(xLoD2 &X, wsReal *Rp)
{
	double R = Brent_fmin(0, 0.5, f2, (void *)(&X), 1e-5);
	wsReal q = (W1 - (wsReal)R) / W2;

	wsReal	R_LoD	= W0;
	FOREACH (mNucFam_it, X.Xm_nf, i) {
		wsUint	N_pi = i->second.Xp_pat->N_idx;
		wsUint	N_mi = i->second.Xp_mat->N_idx;

		wsUint	N_g0 = 0, N_g1 = 0;
		wsUint	N_i1 = X.Ra_phe[N_pi] == WISARD_AFFECTED ? N_mi : N_pi;
		wsUint	N_i2 = X.Ra_phe[N_pi] == WISARD_AFFECTED ? N_pi : N_mi;

		FOREACH (vInt_it, X.Xv_snp, j) {
			N_g0 <<= 2;
			N_g1 <<= 2;

			N_g0 += X.Na_data[N_i1][*j];
			N_g1 += X.Na_data[N_i2][*j];
		}

		/* Cannot test if both genotype is same */
		if (N_g0 == N_g1) continue;

		wsUint N_sum = 0;
		FOREACH (vSampPtr_it, i->second.Xp_childs, k) {
			wsUint N_geno = 0;
			FOREACH (vInt_it, X.Xv_snp, j) {
				N_geno <<= 2;
				N_geno += X.Na_data[(*k)->N_idx][*j];
			}

			/* Same genotype but different label, count this */
			if ((X.Ra_phe[(*k)->N_idx] == WISARD_AFFECTED && N_g1 != N_geno) ||
				(X.Ra_phe[(*k)->N_idx] == WISARD_UNAFFECTED && N_g0 != N_geno))
				N_sum++;
		}

		/* Current family's LoD result */
		wsUint N[2] = { N_sum, (wsUint)(i->second.Xp_childs.size()) };
		R_LoD += log10(pow(q, N[0]) * pow(q, N[1]-N[0]) / pow(0.25, N[1]));
	}

	*Rp = R;
	return R_LoD / 2;
}

void cLodAnalysis::run()
{
	wsRealCst*		Ra_phe	= Cp_IO->getPheno();
	vVariant&	Xv_vrt	= Cp_IO->getVariant();
	char**		Na_data	= Cp_IO->getGenotype();

	if (Xm_curNF.size() == 0) return;

	wsUint j=0;
	if (Cp_anaSM) {
		cTableExporter	C_la("lod.res", "srrs", "Gene-level LoD result", 0,
			4, "GENE", "LOD", "R", "DECISION");
		mGeneDef&		Xm_gene = Cp_anaSM->getGeneDef();

		FOREACHDO (mGeneDef_it, Xm_gene, i, j++) {
			wsReal	R;
			xLoD2	X		= { Xm_curNF, Ra_phe, Na_data, i->second };
			wsReal	R_LoD	= _test(X, &R);

			if (NA(R_LoD)) {
				if (!OPT_ENABLED(remna))
					C_la.write(4, i->first.c_str(), WISARD_NAN, WISARD_NAN, "<NA>");
			} else
				C_la.write(4, i->first.c_str(), R_LoD, R,
				R_LoD>REAL_CONST(3.0) ? "LINKED" : (R_LoD>REAL_CONST(1.73) ?
				"POSSIBLE" : (R_LoD<REAL_CONST(-2.0) ? "NOT_LINKED": "AMBIGUOUS")
				));
		}
	} else {
		cTableExporter C_la("lod.res", "rrs", "Variant-wise LoD result", 1,
			3, "LOD", "R", "DECISION");
		FOREACHDO (vVariant_it, Xv_vrt, i, j++) {
			wsReal	R;
			xLoD	X		= { Xm_curNF, Ra_phe, Na_data, j, *i };
			wsReal	R_LoD	= _test(X, &R);

			if (NA(R_LoD)) {
				if (!OPT_ENABLED(remna))
					C_la.writeVariant(&(*i), WISARD_NAN, WISARD_NAN, "<NA>");
			} else
				C_la.writeVariant(&(*i), R_LoD, R,
					R_LoD>REAL_CONST(3.0) ? "LINKED" : (R_LoD>REAL_CONST(1.73) ?
						"POSSIBLE" : (R_LoD<REAL_CONST(-2.0) ? "NOT_LINKED": "AMBIGUOUS")
				));
		}
	}
}

#endif

} // End namespace ONETOOL
