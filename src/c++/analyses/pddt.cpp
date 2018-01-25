#include <stdlib.h>
#include <string.h>
#include "analyses/pddt.h"
#include "analyses/corr.h"
#include "utils/matrix.h"

namespace ONETOOL {

void findGraph(map<string,xPDDTnode> &X_map, vector<wsReal> &Xa_scores,
	xSample *Xp_src, xSample *Xp_curr, xSample *Xp_tgt, int B_visitingTwin,
	int B_visitingFemTopNode, wsReal R_score, char B_isX);

cPDDTAnalysis::cPDDTAnalysis(cIO *Cp_inpIO) : cAnalysis(Cp_inpIO)
{
	N_sampPDDT = Cp_IO->sizeSample();

	Ra_pddt = sseEmptyMatrix(N_sampPDDT, N_sampPDDT);
}

cPDDTAnalysis::~cPDDTAnalysis()
{
	sseUnmat(Ra_pddt, N_sampPDDT);
}

void cPDDTAnalysis::_getEmpiCorr(wsReal **Ra_pddt, wsFloat **Ra_data,
	wsUint N_1, wsUint N_2, wsUint N_SNP)
{
	if (Cp_IO->isDataComplete()) {
		if (OPT_ENABLED(medcor)) {
			/* Complete genotype, median */
			Ra_pddt[N_1][N_2] = _corMedianComplete(
				Ra_data[N_1], Ra_data[N_2], N_SNP);
		} else if (OPT_ENABLED(ktau)) {
			/* Complete genotype, Kendall's tau */
			Ra_pddt[N_1][N_2] = _corTauComplete(
				Ra_data[N_1], Ra_data[N_2], N_SNP);
		} else {
			/* Complete genotype, mean */
			Ra_pddt[N_1][N_2] = _corMeanComplete(
				Ra_data[N_1], Ra_data[N_2], N_SNP);
		}
	} else {
		if (OPT_ENABLED(medcor)) {
			/* Incomplete genotype, median */
			Ra_pddt[N_1][N_2] = _corMedianComplete(
				Ra_data[N_1], Ra_data[N_2], N_SNP);
		} else if (OPT_ENABLED(ktau)) {
			/* Incomplete genotype, Kendall's tau */
			Ra_pddt[N_1][N_2] = _corTauComplete(
				Ra_data[N_1], Ra_data[N_2], N_SNP);
		} else {
			/* Incomplete genotype, mean */
			Ra_pddt[N_1][N_2] = _corMeanIncomplete(
				Ra_data[N_1], Ra_data[N_2], N_SNP);
		}
	}
}

void cPDDTAnalysis::run()
{
	if (OPT_ENABLED(x) && OPT_ENABLED(x2))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--x", "--x2");

/**/cExporter*	Cp_exporter	= NULL;
	if (IS_ASSIGNED(makecor))
		Cp_exporter = cExporter::summon(
			OPT_ENABLED(hybrid)?"hybrid.cor":"theo.cor"); /* CHECKED */

	wsFloat**	Ra_data		= NULL;
	if (OPT_ENABLED(hybrid))
		Ra_data = Cp_IO->getData();
	wsUintCst		N_SNP		= Cp_IO->sizeVariant();

	/* Print header */
// 	fprintf(H_out, "		");
// 	for (vector<strMap_p>::iterator Xp_itL=X_indExistInSeq.begin() ; Xp_itL!=X_indExistInSeq.end() ; Xp_itL++)
// 		fprintf(H_out, "%s	", X_famIndexToName[(*Xp_itL).second].c_str());
// 	fprintf(H_out, "\n		");
// 	for (vector<strMap_p>::iterator Xp_itL=X_indExistInSeq.begin() ; Xp_itL!=X_indExistInSeq.end() ; Xp_itL++)
// 		fprintf(H_out, "%s	", (*Xp_itL).first.c_str());
// 	fprintf(H_out, "\n");

	/* For sample pair, get PDDT */
	vSampPtr	&Xa_samples = Cp_IO->getSample();
	wsUint		I=0;
	FOREACHDO (vSampPtr_it, Xa_samples, it, I++) {
		xSample *Xp_s1 = *it;
		if (Xp_s1->B_isMissing || Xp_s1->N_idx == SAMP_NODATA)
			halt_fmt(WISARD_SYST_NULL_PHENOTYPE, Xp_s1->S_FID.c_str(),
				Xp_s1->S_IID.c_str());

		/* Set self PDDT to 1.0 */
		if (OPT_ENABLED(x2))
			Ra_pddt[Xp_s1->N_idx][Xp_s1->N_idx] = W1 + (wsReal)(Xp_s1->N_sex==2);
		else
			Ra_pddt[Xp_s1->N_idx][Xp_s1->N_idx] = W1;

		/* (it) against (ait) */
		FOREACH (vSampPtr_it, Xa_samples, ait) {
			xSample	*Xp_s2	= *ait;

			/* PDDT=0 if both FID is different */
			if (Xp_s1->S_FID.compare(Xp_s2->S_FID))	{
				if (OPT_ENABLED(hybrid)) {
					_getEmpiCorr(Ra_pddt, Ra_data, Xp_s1->N_idx, Xp_s2->N_idx, N_SNP);
					Ra_pddt[Xp_s2->N_idx][Xp_s1->N_idx] = Ra_pddt[Xp_s1->N_idx][Xp_s2->N_idx];
				}
				continue;
			}
			/* PDDT=1 if same sample */
			if (Xp_s1 == Xp_s2) {
				if (OPT_ENABLED(x2))
					Ra_pddt[Xp_s1->N_idx][Xp_s1->N_idx] = W1 + (wsReal)(Xp_s1->N_sex==2);
				else
					Ra_pddt[Xp_s1->N_idx][Xp_s1->N_idx] = W1;
				break;
			}
			map<string,xPDDTnode> X_map;

			/* Graph is grown from both of Xp_s1 and Xp_s2,
			 *  bipolar between two samples to get its PDDT
			 *
			 * However, a tree grown from Xp_s1(indicated by 1) should be directed from
			 *   itself to ascendants, and another tree grown from Xp_s2 (indicated by 0)
			 *   should be directed from ascendants to itself */

//			pverbose("[%s] -> [%s]\n", Xp_s1->S_IID.c_str(),
//				Xp_s2->S_IID.c_str());
			buildGraph(X_map, Xp_s2, 0, buildGraph(X_map, Xp_s1, 1));

			/* Xp_s1???? Xp_s2?? a?ве? */
			vector<wsReal> Xa_scores;
			char B_isX = OPT_ENABLED(x);
			if (OPT_ENABLED(x2))
				B_isX = 2;

#ifdef DEBUG_KINSHIP
			pverbose("[%s] <=> [%s]\n", Xp_s1->S_IID.c_str(),
				Xp_s2->S_IID.c_str());
#endif

			/* In case of getting X, it should be 0.5 */
			wsReal R_init = B_isX?(Xp_s1->N_sex==2?W2:W1):W1;
			findGraph(X_map, Xa_scores, Xp_s1, Xp_s1, Xp_s2, 0, -1, R_init, B_isX);

			/* Get PDDT value */
			wsReal R_finVal = W0;
			wsUint N_tPath = 0;
#ifdef DEBUG_KINSHIP
			pverbose("Coefficients :");
#endif
			FOREACH (vector<wsReal>::iterator, Xa_scores, sit) {
#ifdef DEBUG_KINSHIP
				verbosenf(" %g", *sit);
#endif
				R_finVal += *sit;
				if (*sit > W0) N_tPath++;
			}

			/* If both are UNRELATED and --hybrid */
			/* FIXME : Ra_cm should be applied, currently set to NULL */
			if (OPT_ENABLED(hybrid) && Ra_pddt[Xp_s1->N_idx][Xp_s2->N_idx] == W0)
				_getEmpiCorr(Ra_pddt, Ra_data, Xp_s1->N_idx, Xp_s2->N_idx, N_SNP);
#ifdef DEBUG_KINSHIP
			verbosenf(" = %g\n\n", R_finVal);
#endif
			/* Since it is symmetric, copy here to symmetric side */
			Ra_pddt[Xp_s1->N_idx][Xp_s2->N_idx] = R_finVal;
			Ra_pddt[Xp_s2->N_idx][Xp_s1->N_idx] = Ra_pddt[Xp_s1->N_idx][Xp_s2->N_idx];
			if (Cp_exporter && OPT_ENABLED(mqls)) {
				Cp_exporter->fmt("%s	%s	%s	%g\n",
					Xp_s1->S_FID.c_str(), Xp_s1->S_IID.c_str(),
					Xp_s2->S_IID.c_str(),
					Ra_pddt[Xp_s1->N_idx][Xp_s2->N_idx] / W2);
			}
		}
		notice("%d/%d kinship coefficient calculated\r", I, N_sampPDDT);
	}
	LOG("%d/%d kinship coefficient calculated\n", I, N_sampPDDT);

	if (Cp_exporter) {
		if (!OPT_ENABLED(corpair)) {
			cSymMatrix	C_pddt(Ra_pddt, N_sampPDDT, 1);
			vStr		V_iid;
			delete		Cp_exporter;

			FOREACH (vSampPtr_it, Xa_samples, i)
				V_iid.push_back((*i)->S_IID);
			C_pddt.file(OPT_ENABLED(hybrid)?"hybrid.cor":"theo.cor",
				&V_iid);
// 		for (wsUint i=0 ; i<N_sampPDDT ; i++) {
// 			for (wsUint j=0 ; j<N_sampPDDT ; j++) {
// 				Cp_exporter->fmt("%g	", Ra_pddt[i][j]);
// 			}
// 			Cp_exporter->put("\n");
// 		}
		} else {
			if (!OPT_ENABLED(mqls)) {
				Cp_exporter->put("SAMP1	SAMP2	KINSHIP2\n");
				for (wsUint i=0 ; i<N_sampPDDT ; i++) {
					string &S_i1 = Xa_samples[i]->S_IID;

					for (wsUint j=0 ; j<=i ; j++) {
						Cp_exporter->fmt("%s	%s	%g\n", S_i1.c_str(),
							Xa_samples[j]->S_IID.c_str(), Ra_pddt[i][j]);
					}
				}
			}
			delete Cp_exporter;
		}
	}
}

void findGraph(map<string,xPDDTnode> &X_map, vector<wsReal> &Xa_scores,
	xSample *Xp_src, xSample *Xp_cur, xSample *Xp_tgt, int B_visitingTwin,
	int B_visitingFemTopNode, wsReal R_score,
	char B_isX)
{
	xPDDTnode &X_currNode = X_map[Xp_cur->S_IID];

	/* Target found */
	if (Xp_cur == Xp_tgt) {
// 		if (B_isX) {
// 			if (B_visitingFemTopNode == -1) {
// 				/* Direct relation */
// 				if (X_map[Xp_src->S_IID].N_depth > X_map[Xp_tgt->S_IID].N_depth)
// 					B_visitingFemTopNode = Xp_tgt->N_sex-1;
// 				else
// 					B_visitingFemTopNode = Xp_src->N_sex-1;
// 			}
// 		} else
			B_visitingFemTopNode = 0;

//		pverbose("	FOUND [%s] : score %g twin %d fem %d\n",
//			Xp_tgt->S_IID.c_str(), R_score,
//			B_visitingTwin, B_visitingFemTopNode);

		/* twin?? ??????????? *2???? ??? */
		Xa_scores.push_back(R_score * (1+B_visitingTwin) * (1+B_visitingFemTopNode));
		return;
	}

	/* For all connected node, search and expand */
	FOREACH (vSampPtr_it, X_currNode.Xa_next, it) {
		/* Twin?? ???, ??????? ???????? ???????? twin????
			* ?б└???? ??????? */
		if (B_visitingTwin == 0) {
			FOREACH (vSampPtr_it, (*it)->Xp_twins, iit)
				if (mapIsExist(X_map, (*iit)->S_IID)) {
					B_visitingTwin = 1;
					break;
				}
		}

		/* If this is top level parent and FEMALE */
		if ((X_currNode.N_depth < X_map[(*it)->S_IID].N_depth)) {
			if (Xp_cur->N_sex == 2)
				B_visitingFemTopNode = 1;
			else
				B_visitingFemTopNode = 0;
		}

		wsReal R_mult = REAL_CONST(0.5);

		/* Special for X chromosome */
		if (B_isX) {
 			if (Xp_cur->N_sex == 1)
				/* M -> F = 1 */
				R_mult = (*it)->N_sex == 2 ? W1 : W0;
			else if (Xp_cur->N_sex == 2) {
				if ((*it)->N_sex != 1 && (*it)->N_sex != 2)
					R_mult = W0;
			} else
				R_mult = W0;
		}
#ifdef DEBUG_KINSHIP
 		pverbose("	[%s]->[%s], %g*%g=%g, twin %d, fem %d\n", Xp_cur->S_IID.c_str(),
 			(*it)->S_IID.c_str(), R_score, R_mult, R_score*R_mult,
 			B_visitingTwin, B_visitingFemTopNode);
#endif
		findGraph(X_map, Xa_scores, Xp_src, *it, Xp_tgt, B_visitingTwin,
			B_visitingFemTopNode, R_score*R_mult, B_isX);
	}
}

int buildGraph(map<string,xPDDTnode> &X_map, xSample *Xp_t, char B_fromMe,
	int N_depth/*=0*/)
{
	/* Mark myself to N_depth */
	xPDDTnode &X_currNode = X_map[Xp_t->S_IID];
	X_currNode.N_depth = N_depth;

	/* If this is NOT founder */
	if (Xp_t->Xp_pat && Xp_t->Xp_mat) {
		if (B_fromMe) {
			/* If B_fromMe==1, extend from this sample to its parent */
			X_currNode.Xa_next.push_back(Xp_t->Xp_pat);
			X_currNode.Xa_next.push_back(Xp_t->Xp_mat);
		} else {
			/* Otherwise, extend from its parent to this sample */
			char B_insP=1, B_insF=1;

			/* Xp_t?? ??? Xp_pat?? ???? ??? ?? ?? ????????? ?? */
			vSampPtr &Xa_vt = X_currNode.Xa_next;
			FOREACH (vSampPtr_it, Xa_vt, it) {
				if (*it == Xp_t->Xp_pat) {
					B_insP = 0;
					it=Xa_vt.erase(it);
					if (it==Xa_vt.end()) break;
				}
				if (*it == Xp_t->Xp_mat) {
					B_insF = 0;
					it=Xa_vt.erase(it);
					if (it==Xa_vt.end()) break;
				}
			}
			if (B_insP) X_map[Xp_t->Xp_pat->S_IID].Xa_next.push_back(Xp_t);
			if (B_insF) X_map[Xp_t->Xp_mat->S_IID].Xa_next.push_back(Xp_t);
		}
		int s1 = buildGraph(X_map, Xp_t->Xp_pat, B_fromMe, N_depth-1);
		int s2 = buildGraph(X_map, Xp_t->Xp_mat, B_fromMe, N_depth-1);

		return s1>s2?s2:s1;
	} else if (Xp_t->Xp_pat || Xp_t->Xp_mat)
		halt_fmt(WISARD_SYST_NULL_PARENT, Xp_t->Xp_pat?"Mother":"Father",
			Xp_t->S_FID.c_str(), Xp_t->S_IID.c_str(),
			Xp_t->Xp_pat?"Father":"Mother");

	if (Xp_t->Xp_twins.size()) {
		if (!B_fromMe) {
			/* In case of this sample is twin */
			FOREACH (vSampPtr_it, Xp_t->Xp_twins, it)
				X_map[(*it)->S_IID].Xa_next.push_back(Xp_t);
		}
	}

	return N_depth;
}

} // End namespace ONETOOL
