#include "analyses/fst.h"
#include "analyses/setmgr.h"
#include "utils/matrix.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cFstAnalysis::cFstAnalysis(cIO *Cp_inpIO, cSetManagerAnalysis *Cp_inpAnaSM) :
	cAnalysis(Cp_inpIO)
{
	Cp_anaSM = Cp_inpAnaSM;
}

cFstAnalysis::~cFstAnalysis()
{

}

wsReal _getFst(wsUint N_samp, char **Na_data, vInt& Xv_idx_, xMaf* Xp_maf,
	vector<vInt>& Xv_pop)
{
	/* Resize Xv_idx */
	vInt Xv_idx;
	FOREACH (vInt_it, Xv_idx_, i) {
		if (Xp_maf[*i].R_allMaf > W0)
			Xv_idx.push_back(*i);
	}
	if (Xv_idx.size() == 0)
		return WISARD_NAN;

	/* Step 1 : Get population-wise MAF */
	wsUint		N_vrt		= (wsUint)Xv_idx.size();
	wsUint		N_pop		= (wsUint)Xv_pop.size();
/**/wsUint**	Na_nallele	= NULL;
/**/wsUint**	Na_cnt		= NULL;
	wsAlloc(Na_nallele, wsUint*, N_vrt);
	wsAlloc(Na_cnt, wsUint*, N_vrt);
	for (wsUint i=0 ; i<N_vrt ; i++) {
		wsCalloc(Na_nallele[i], wsUint, N_pop);
		wsCalloc(Na_cnt[i], wsUint, N_pop);
	}

	/* For each pop */
	wsUint J = 0;
	FOREACHDO (vector<vInt>::iterator, Xv_pop, i, J++) {
		vInt& Xv_cur = *i;

		FOREACH (vInt_it, Xv_cur, k) {
			wsUint I = 0;
			FOREACHDO (vInt_it, Xv_idx, i, I++) {
				/* *k = smp index */
				/* *i = vrt index */
				char N_geno = Na_data[*k][*i];
				if (isMissing(N_geno)) continue;

				/* J = pop index from 0 */
				/* I = vrt index from 0 */
				Na_nallele[I][J] += N_geno;
				Na_cnt[I][J]++;
			}
		}
	}

	/* Get pop-wise MAF */
/**/wsMat		Ra_maf		= sseMatrix(N_vrt, N_pop);
	for (wsUint i=0 ; i<N_vrt ; i++)
		for (wsUint j=0 ; j<N_pop ; j++)
			Ra_maf[i][j] = (wsReal)Na_nallele[i][j] / (wsReal)(Na_cnt[i][j] << 1);

	for (wsUint i=0 ; i<N_vrt ; i++) {
		DEALLOC(Na_nallele[i]);
		DEALLOC(Na_cnt[i]);
	}
	DEALLOC(Na_nallele);
	DEALLOC(Na_cnt);

	/* n */
/**/wsReal*		Ra_n		= sseVector((wsUint)Xv_pop.size());
	J = 0;
	FOREACH (vector<vInt>::iterator, Xv_pop, i)
		Ra_n[J++] = (wsReal)i->size();

	/* Get variant-wise mean of MAFs */
/**/wsReal*		Ra_varMaf	= sseMpV(Ra_maf, N_vrt, N_pop, Ra_n, W1/(wsReal)N_samp);

	/* msp_1<-1/(npop-1)*(apply(p_1,1,function(x){((x-p_1_avg)^2)})%*%n) */
/**/wsMat		Ra_msp11	= sseMatrix(N_vrt, N_pop);
	sseMsVsq(Ra_maf, Ra_varMaf, Ra_msp11, N_vrt, N_pop);
	sseFree(Ra_varMaf);
/**/wsReal*		Ra_msp1		= sseMpV(Ra_msp11, N_vrt, N_pop, Ra_n,
		W1/(wsReal)(N_pop-1));
	sseUnmat(Ra_msp11, N_vrt);
	wsMat		Ra_msp21	= sseMatrix(N_vrt, N_pop);
	sseMp1p(Ra_maf, Ra_msp21, N_vrt, N_pop);
	sseUnmat(Ra_maf, N_vrt);
/**/wsReal*		Ra_msg1		= sseMpV(Ra_msp21, N_vrt, N_pop, Ra_n,
		W1/(wsReal)(N_samp-N_pop));
	sseUnmat(Ra_msp21, N_vrt);
	sseFree(Ra_n);
	/* sum(msp_1-msg_1+msp_2-msg_2) */
	wsReal R_nume = sseVsVsum(Ra_msp1, Ra_msg1, N_vrt, 0, W2);

	/*  
	msg_1<-1/(N-npop)*(apply(p_1,1,function(x){x*(1-x)})%*%n) */
	/* msp_2<-1/(npop-1)*(apply(p_2,1,function(x){((x-p_2_avg)^2)})%*%n)
	msg_2<-1/(N-npop)*(apply(p_2,1,function(x){x*(1-x)})%*%n) */

	/* n_c<-1/(npop-1)*(N-sum(n^2)/N) */
	wsReal R_sqSum = W0;
	FOREACH (vector<vInt>::iterator, Xv_pop, i) {
		wsReal	R_sz	= (wsReal)i->size();
		R_sqSum += SQR(R_sz);
	}
	wsReal	R_samp	= (wsReal)N_samp;
	wsReal	R_c		= (R_samp - R_sqSum / R_samp) / (wsReal)(N_pop-1);

	/* sum(msp_1+(n_c-1)*msg_1+msp_2+(n_c-1)*msg_2) */
	sseVpC(Ra_msg1, R_c-W1, Ra_msg1, N_vrt);
	wsReal R_deno = sseVaVsum(Ra_msp1, Ra_msg1, N_vrt, W2);
	sseFree(Ra_msg1);
	sseFree(Ra_msp1);

	/* theta<-sum(msp_1-msg_1+msp_2-msg_2)/sum(msp_1+(n_c-1)*msg_1+msp_2+(n_c-1)*msg_2) */
	wsReal R_theta = R_nume / R_deno;

	return R_theta;
}

void cFstAnalysis::run()
{
	vector<vInt>	Xv_pop;
	wsUintCst			N_samp	= getIO()->sizeSample();
	char**			Na_data	= getIO()->getGenotype();
	vVariant&		Xv_vrt	= getIO()->getVariant();
	xMaf*			Xp_maf	= getIO()->getMAF();
	vSampPtr&		Xv_samp	= getIO()->getSample();

	/* Build population info */

	/* Find maximum number */
	int N_maxPop = 0;
	FOREACH (vSampPtr_it, Xv_samp, i) {
		xSample& X_s = *(*i);
		if (X_s.N_grpFst == -1)
			halt("Sample [%s:%s] does not contains Fst group info!",
				X_s.S_FID.c_str(), X_s.S_IID.c_str());
		if (X_s.N_grpFst > N_maxPop)
			N_maxPop = X_s.N_grpFst;
	}
	/* Can't do Fst when # group = 1 */
	if (N_maxPop == 0)
		halt("Only one group exists in Fst group info, can't do anymore!");
	LOG("[%d] groups found for Fst computation\n", N_maxPop+1);
	Xv_pop.resize(N_maxPop+1);
	/* Build info */
	FOREACH (vSampPtr_it, Xv_samp, i)
		Xv_pop[(*i)->N_grpFst].push_back((*i)->N_idx);

	if (Cp_anaSM) {
		cTableExporter C_fst("fst.res", "sir", "Gene-level Fst result", 0,
			3, "GENE", "NVARIANT", "F_ST");
		mGeneDef&	Xm_gd = Cp_anaSM->getGeneDef();
		FOREACH (mGeneDef_it, Xm_gd, i) {
			wsReal R_Fst = _getFst(N_samp, Na_data, i->second, Xp_maf, Xv_pop);
			/* --remna */
			if (NA(R_Fst) && OPT_ENABLED(remna)) continue;
			/* Write out the result */
			C_fst.write(3, i->first.c_str(), i->second.size(), R_Fst);
		}
	} else {
		vStr Sv_hdr;
		Sv_hdr.push_back("F_ST");

		cTableExporter C_fst("fst.res", Sv_hdr, "r",
			"Variant-level Fst result", 1);
		wsUint			I		= 0;
		FOREACHDO (vVariant_it, Xv_vrt, i, I++) {
			vInt X;
			X.push_back(I);
			wsReal R_Fst = _getFst(N_samp, Na_data, X, Xp_maf, Xv_pop);
			/* --remna */
			if (NA(R_Fst) && OPT_ENABLED(remna)) continue;
			/* Write out the result */
			C_fst.writeVariant(&(*i), R_Fst);
		}
	}
}

#endif

} // End namespace ONETOOL
