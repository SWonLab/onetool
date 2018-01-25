#include "analyses/fam.h"
#include "utils/matrix.h"
#include "utils/vis.h"
#include "utils/stat.h"

namespace ONETOOL {

cFamStrAnalysis::cFamStrAnalysis(cIO *Cp_inpIO) :
	cAnalysis(Cp_inpIO)
{
	Na_famIndices	= NULL;
	B_sameFamStr	= 0;
}

cFamStrAnalysis::~cFamStrAnalysis()
{
	Xm_couple2data.clear();
	Xm_missingFounderIID2data.clear();
	DEALLOC(Na_famIndices);
}

int getFamMemberAvail(xSample *Xp_s, void *Vp_data)
{
	xFamily *Xp_fm = (xFamily *)Vp_data;
	if (Xp_s->B_isMissing || Xp_s->N_idx == SAMP_NODATA) {
		Xp_s->N_isVisited = 1;
		return 0;
	}

	/* Insert to vector */
	Xp_fm->Xv_members.push_back(Xp_s->N_idx);
	Xp_s->N_isVisited = 1;
	return 0;
}

void cFamStrAnalysis::run()
{
	mFam&		Xm_famData = Cp_IO->getFamilyData();
	vSampPtr	Xa_currTargets;
	/* Is all family in the dataset have same family structure? */
	B_sameFamStr	= 1;
	wsUint		N_sameDeg		= 0;
	cSymMatrix	M_detMat;
	wsUint		N_idxFam		= 0;
	wsAlloc(Na_famIndices, wsUint, Cp_IO->sizeSample());

	FOREACHDO (mFam_it, Xm_famData, it, N_idxFam++) {
		vSampPtr &Xa_samp = Cp_IO->getSample();

		/* Insert family members */
		it->second.setUnvisited();
		it->second.visit(getFamMemberAvail, &(it->second));

		/* Skip if the family have no available member */
		if (it->second.Xv_members.size() == 0)
			continue;

		/* Recording fammem index */
		vInt& Xv_mem = it->second.Xv_members;
		FOREACH (vInt_it, Xv_mem, i)
			Na_famIndices[*i] = N_idxFam;

		/* For avail members, build matrix */
		if (B_sameFamStr) {
			cSymMatrix *Mp_cDetMat = NULL;

			if (M_detMat.row() == 0) {
				M_detMat.init((wsUint)Xv_mem.size());
				Mp_cDetMat = &M_detMat;
			}
			/* Unvalidate if family size are diff. */
			else if (M_detMat.row() != Xv_mem.size()) {
				B_sameFamStr = 0;
				M_detMat.rem();
				continue;
			} else {
				Mp_cDetMat = new cSymMatrix((wsUint)Xv_mem.size());
			}
			//LOG("Matrix %x made\n", Mp_cDetMat);

			/* If same, construct */
			int I = 0;
			wsSym Ra_detMat = Mp_cDetMat->get();
			FOREACHDO (vInt_it, Xv_mem, i, I++) {
				/* Skip if founder */
				if (!Xa_samp[*i]->Xp_mat) continue;

				/* Check index of parents, and set to 1 */
				xSample *Xp_pat = Xa_samp[*i]->Xp_pat;
				xSample *Xp_mat = Xa_samp[*i]->Xp_mat;
				int		B_pfound = -1;
				int		B_mfound = -1;
				int		II = 0;
				FOREACHDO (vInt_it, Xv_mem, ii, II++) {
					if (Xp_pat == Xa_samp[*ii])
						B_pfound = II;
					else if (Xp_mat == Xa_samp[*ii])
						B_mfound = II;
					if (B_pfound!=-1 && B_mfound!=-1) break;
				}
				/* Mark if found */
				if (B_pfound != -1) {
					if (I<B_pfound) Ra_detMat[B_pfound][I] = 1;
					else Ra_detMat[I][B_pfound] = 1;
				}
				if (B_mfound != -1) {
					if (I<B_mfound) Ra_detMat[B_mfound][I] = 1;
					else Ra_detMat[I][B_mfound] = 1;
				}
			}
			/* Get sum */
			wsUint N_sum = (wsUint)M_detMat.sum();
			if (N_sameDeg == 0)
				/* Set initial degree */
				N_sameDeg = N_sum;
			else if (N_sameDeg != N_sum)
				/* Unvalidate if degree is not same */
				B_sameFamStr = 0;
			if (Mp_cDetMat != &M_detMat) {
//				LOG("Matrix %x removed\n", Mp_cDetMat);
				delete Mp_cDetMat;
			}
		}
	}
	M_detMat.rem();

	/* Make notification */
	if (B_sameFamStr)
		LOG("This dataset have *identical* family structures\n");

	/* Construct initial sample list(consisted of founders) for build nuclear family structure */
	FOREACH (mFam_it, Xm_famData, it) {
		vSampPtr Xp_founders = it->second.Xp_founders;
		it->second.setUnvisited();

		FOREACH (vSampPtr_it, Xp_founders, fit) {
			if ((*fit)->B_isMissing)
				Xm_missingFounderIID2data.insert(make_pair((*fit)->S_IID, *fit));
			Xa_currTargets.push_back(*fit);
			(*fit)->N_isVisited = 1;
		}
	}
	//LOG("%d missing founder were found\n", Xm_missingFounderIID2data.size());

	/* Do searching */
	vSampPtr	Xa_newTargets;
	for (;Xa_currTargets.size();) {
		Xa_newTargets.clear();

		/* For all samples in current loop */
		FOREACH (vSampPtr_it, Xa_currTargets, it) {
			xSample &X_curSamp = *(*it);

			/* 현재 sample의 모든 배우자와 pair string을 만든 후 확인 */
			FOREACH (vSampPtr_it, X_curSamp.Xp_spouses, sit) {
				if (X_curSamp.N_sex == 1)	_buildNucFam(X_curSamp, *(*sit));
				else						_buildNucFam(*(*sit), X_curSamp);
			}

			/* 현재 sample의 자식을 Xa_currTargets에 삽입 (이미 있는 대상 빼고 */
			FOREACH (vSampPtr_it, X_curSamp.Xp_childs, cit) {
				/* Pass already visited sample */
				if ((*cit)->N_isVisited) continue;

				vSampPtr_it tit;
				for (tit=Xa_newTargets.begin() ; tit!=Xa_newTargets.end() ; tit++)
					if ((*tit)->S_IID.compare((*cit)->S_IID) == 0)
						break;

				if (tit==Xa_newTargets.end()) {
					(*cit)->N_isVisited = 1;
					Xa_newTargets.push_back(*cit);
				}
			}
		}

		/* 리스트 대치 */
		Xa_currTargets.clear();
		Xa_currTargets = Xa_newTargets;
	}

	/* Clear visiting status */
	FOREACH (mFam_it, Xm_famData, it)
		it->second.setUnvisited();

	/* Export */
	if (OPT_ENABLED(ncsummary)) {
/**/	cExporter* Cp_nuc = cExporter::summon("structure.nucl.fam");

		FOREACH (mNucFam_it, Xm_couple2data, it) {
			Cp_nuc->fmt("Nuclear family [%s]\n", it->first.c_str());
			if (it->second.Xp_mat->Xp_mat == NULL && it->second.Xp_mat->Xp_pat == NULL &&
				it->second.Xp_pat->Xp_mat == NULL && it->second.Xp_pat->Xp_pat == NULL)
				Cp_nuc->put("Both-founder\n");


			FOREACH (vSampPtr_it, it->second.Xp_childs, cit) {
				Cp_nuc->fmt("\tChild %s\n", (*cit)->S_IID.c_str());
			}
		}
		delete Cp_nuc;
	}
}

void cFamStrAnalysis::_buildNucFam(xSample &X_husband, xSample &X_wife)
{
	/* Key = <husband,wife> */
	string X_coupleKey = X_husband.S_IID + "," + X_wife.S_IID;

	/* 맵에 들어가 있는지 확인하고 없으면 넣기 */
	mNucFam_it X_resFind = Xm_couple2data.find(X_coupleKey);
	if (X_resFind != Xm_couple2data.end())
		return;

	/* Set the nuclear family consists of above husband and wife */
	xNucFamily X_newNucFam;
	X_newNucFam.Xp_pat = &X_husband;
	X_newNucFam.Xp_mat = &X_wife;

	/* Insert children having above parents to be their parents */
	FOREACH (vSampPtr_it, X_husband.Xp_childs, hit) {
		FOREACH (vSampPtr_it, X_wife.Xp_childs, wit) {
			if ((*hit)->S_IID.compare((*wit)->S_IID) == 0) {
				X_newNucFam.Xp_childs.push_back(*hit);
				break;
			}
		}
	}
	Xm_couple2data.insert(make_pair(X_coupleKey, X_newNucFam));
}

mNucFam& cFamStrAnalysis::getNucFamData()
{
	return Xm_couple2data;
}

mSampPtr& cFamStrAnalysis::getMissingFounderData()
{
	return Xm_missingFounderIID2data;
}

double chi2_3x3(wsUint t[3][3])
{
	int a = 3;
	int b = 3;

	wsUint rows[3] = { 0, };
	wsUint cols[3] = { 0, };
	double sum = 0;

	for (int i=0; i<a; i++)
		for (int j=0; j<b; j++)
		{
			rows[i] += t[i][j];
			cols[j] += t[i][j];
			sum += t[i][j];
		}

	// Sum (O-E)^2 / E
	wsReal exp[3][3] = { { REAL_CONST(0), }, };
	wsReal chisq = 0;
	for (int i=0; i<a; i++)
		for (int j=0; j<b; j++) {
			exp[i][j] += ( rows[i] * cols[j] ) / sum;
			wsReal tmp = t[i][j] - exp[i][j];
			tmp *= tmp;
			tmp /= exp[i][j];
			chisq += tmp;
		}
	return chisq;
}

double chiPsym_3x3(wsUint t[3][3])
{
	// Test for symmetry in a n x n table
	int a = 3;
	int b = 3;
	if ( a != b )
		return -1;

	int df = a*(a-1)/2;
	double x = 0;
	for (int i=0; i<a; i++)
		for (int j=0; j<i; j++) {
			double tmp = t[i][j] - t[j][i];
			tmp *= tmp;
			tmp /= t[i][j] + t[j][i];
			x += tmp;
		}
	return PVchisq(x, df);
}

void parenTDT(cIO *Cp_IO, wsUint N_idx, wsReal *Ra_Y, cFamStrAnalysis *Cp_anaFS)
{
	char**		Na_data	= Cp_IO->getGenotype();
	mNucFam&	Xv_nfam	= Cp_anaFS->getNucFamData();
	xVariant&	X_var	= Cp_IO->getVariant()[N_idx];
	// Transmission counts
	double t1 = 0;
	double t2 = 0;

	// Count over families

	wsUint f = 0;
	FOREACHDO (mNucFam_it, Xv_nfam, _F, f++) {
		xNucFamily& NF = _F->second;

		int trA = 0;  // transmitted allele from first het parent
		int unA = 0;  // untransmitted allele from first het parent

		int trB = 0;  // transmitted allele from second het parent
		int unB = 0;  // untransmitted allele from second het parent

		xSample *pat = NF.Xp_pat;
		xSample *mat = NF.Xp_mat;
		vSampPtr &kid = NF.Xp_childs;

		bool pat1 = Na_data[pat->N_idx][N_idx] > 0;
		bool pat2 = Na_data[pat->N_idx][N_idx] == 2;

		bool mat1 = Na_data[mat->N_idx][N_idx] > 0;
		bool mat2 = Na_data[mat->N_idx][N_idx] == 2;

		// We need two genotyped parents, with 
		// at least one het
		if (pat1 == pat2 && mat1 == mat2 ) continue;
		if ( ( pat1 && !pat2 ) || ( mat1 && !mat2 ) )  continue;


		// Consider all offspring in nuclear family
		wsUint c = 0;
		FOREACHDO (vSampPtr_it, kid, _K, c++) {
			xSample *K  = *_K;

			// Only consider affected children
			if (Ra_Y[K->N_idx] != WISARD_AFFECTED) continue;
//			if ( ! kid[c]->aff ) continue;

			bool kid1 = Na_data[K->N_idx][N_idx] > 0;
			bool kid2 = Na_data[K->N_idx][N_idx] == 2;
//			bool kid1 = kid[c]->one[l];
	//		bool kid2 = kid[c]->two[l];

			// Skip if offspring has missing genotype
			if ( kid1 && !kid2 ) continue;

			// We've now established: no missing genotypes
			// and at least one heterozygous parent

			// Kid is 00

			if ( (!kid1) && (!kid2) ) 
			{
				if ( ( (!pat1) && pat2 ) && 
					( (!mat1) && mat2 ) )
				{ trA=1; unA=2; trB=1; unB=2; }
				else 
				{ trA=1; unA=2; } 
			}
			else if ( (!kid1) && kid2 )  // Kid is 01
			{
				// het dad
				if (pat1 != pat2 )
				{
					// het mum
					if ( mat1 != mat2 )
					{ trA=1; trB=2; unA=2; unB=1; }
					else if ( !mat1 ) 
					{ trA=2; unA=1; }
					else { trA=1; unA=2; }
				}
				else if ( !pat1 ) 
				{
					trA=2; unA=1; 
				}		    
				else
				{
					trA=1; unA=2;
				}
			}
			else // kid is 1/1
			{

				if ( ( (!pat1) && pat2 ) && 
					( (!mat1) && mat2 ) )
				{ trA=2; unA=1; trB=2; unB=1; }
				else 
				{ 
					trA=2; unA=1;
				}
			}

			// We have now populated trA (first transmission) 
			// and possibly trB also 

			// Increment transmission counts
			if (trA==1) t1++;
			if (trB==1) t1++;
			if (trA==2) t2++;
			if (trB==2) t2++; 
		} // next offspring in family
	}  // next nuclear family



	/////////////////////////////////////////////
	// Consider parental discordance information

	double p1 = 0; 
	double p2 = 0;
	double d1 = 0;
	double d2 = 0;

//	true if there is at least one nuc fam have different phenotypes between pat and mat
//	if (par::discordant_parents) {
		// Count over families
		f = 0;
		FOREACHDO (mNucFam_it, Xv_nfam, _F, f++) {
			xNucFamily &NF = _F->second;

			// Requires parental discordance...( == phenotype is diff. )
			if (Ra_Y[NF.Xp_mat->N_idx] == Ra_Y[NF.Xp_pat->N_idx] ||
				isMissingReal(Ra_Y[NF.Xp_mat->N_idx]) ||
				isMissingReal(Ra_Y[NF.Xp_pat->N_idx])) continue;
//			if ( ! family[f]->discordant_parents ) 
	//			continue;

			xSample *pat = NF.Xp_pat;
			xSample *mat = NF.Xp_mat;

			char N_Gp = Na_data[pat->N_idx][N_idx];
			char N_Gm = Na_data[mat->N_idx][N_idx];

			// ...and that both are genotyped
			if (isMissing(N_Gp) || isMissing(N_Gm)) continue;

			// Get number of 'F' alleles that the affected parent has
			// above the unaffected; this count is p1/d1

			//  excess T alleles in affected -> p1/d1
			//  excess F alleles in unaffected -> p2/d2

			// skip if both hom
			if (N_Gp != 1 && N_Gm != -1)
				continue;  // d = 0;
			else if (Ra_Y[pat->N_idx] == WISARD_AFFECTED) { // affected pat 
				if (N_Gp == 0)
//				if ( (!pat1) && (!pat2) ) // F/F
				{
					// mat will either be T/T or F/T
					//if ( mat1 ) d1++; // two extra T
					if (N_Gm) d1++;
					else p1++; // one extra T
				}
				else if (N_Gp == 1)
//				else if ( (!pat1 ) && pat2 ) // pat F/T
				{
					// mat either T/T or F/F
					if (N_Gm)
//					if ( mat1 ) 
						p1++; // one extra T
					else
						p2++; // one less T
				}
				else // pat must be T/T
				{
					if (N_Gm < 2) d2++;
					// mat will either be F/F or F/T
					//if ( ! mat2 ) d2++; // two less T
					else p2++; // one less T
				}
			} 
			else // affected mat / score other direction
			{
				if (!N_Gp)
//				if ( (!pat1) && (!pat2) ) // F/F
				{
					// mat will either be T/T or F/T
					if (N_Gm) d2++;
//					if ( mat1 ) d2++;
					else p2++;
				}
				else if (N_Gp == 1)
//				else if ( (!pat1 ) && pat2 ) // pat F/T
				{
					// mat either T/T or F/F
					if (N_Gm)
//					if ( mat1 ) 
						p2++;
					else
						p1++;	      
				}
				else // pat must be T/T
				{
					// mat will either be F/F or F/T
					if (N_Gm < 2) d1++;
//					if ( ! mat2 ) d1++;
					else p1++;
				}
			}
		}
//	}


	///////////////////////////////////
	// General family test

// 	if ( par::mating_tests )  {
// 		wsUint parenMT[3][3] = { 0, };
// 
// 		// Count over families
// 		wsUint f = 0;
// 		FOREACHDO (mNucFam_it, Xv_nfam, _F, f++) {
// 			xNucFamily &NF = _F->second;
// 
// 			xSample *pat = NF.Xp_pat;
// 			xSample *mat = NF.Xp_mat;
// 
// 			char N_Gp = Na_data[pat->N_idx][N_idx];
// 			char N_Gm = Na_data[mat->N_idx][N_idx];
// 
// 			// ...and that both are genotyped
// 			if (isMissing(N_Gp) || isMissing(N_Gm)) continue;
// 
// 			// Do count
// 			++parenMT[N_Gp][N_Gm];
// 		}
// 
// 
// 		double mean1 = 0, mean2 = 0; 
// 		int total = 0;
// 		for(int i=0; i<=2; i++)
// 			for (int j=0; j<=2; j++) {
// 				mean1 += parenMT[i][j] * i;
// 				mean2 += parenMT[i][j] * j;
// 				total += parenMT[i][j];
// 			}
// 
// 		mean1 /= (double)total;
// 		mean2 /= (double)total;
// 
// 		double var1 = 0, var2 = 0, covar = 0;
// 		for(int i=0; i<=2; i++)
// 			for (int j=0; j<=2; j++)
// 			{
// 				var1 += ( i - mean1 )*(i-mean1)*parenMT[i][j];
// 				var2 += ( j - mean2 )*(j-mean2)*parenMT[i][j];
// 				covar += ( i - mean1 )*(j-mean2)*parenMT[i][j];
// 			}
// 
// 		var1 /= (double)total - 1.0;
// 		var2 /= (double)total - 1.0;
// 		covar /= (double)total - 1.0;
// 
// 		double r = covar / sqrt( var1 * var2 );
// 		double t1 = chi2_3x3(parenMT);
// 
// 		MT << setw(4) << locus[l]->chr << " " 
// 			<< setw(par::pp_maxsnp) << locus[l]->name << " " 	
// 			<< setw(12) << locus[l]->bp << " "
// 			<< setw(8) << locus[l]->freq << " "
// 			<< setw(12) << t1 << " ";	  
// 
// 		double t2 = chiPsym_3x3(parenMT);
// 		MT << setw(12) << t2 << " ";
// 		MT << setw(12) << r << " ";
// 
// 		for(int i=0; i<=2; i++)
// 			for (int j=0; j<=2; j++)
// 				MT <<  setw(5) << parenMT[i][j] << " ";
// 
// 		MT << "\n";
//	}




	/////////////////////////////
	// Finished counting: now compute
	// the statistics
	double tdt_chisq, par_chisq, com_chisq;
	tdt_chisq = par_chisq = com_chisq = -1;

	// Basic TDT test

	if (t1+t2 > 0)
		tdt_chisq = ((t1-t2)*(t1-t2))/(t1+t2);

//	if (par::discordant_parents) {
		// parenTDT
		if ( p1+p2+d1+d2 > 0 ) 
			par_chisq = (((p1+2*d1)-(p2+2*d2))*((p1+2*d1)-(p2+2*d2)))
			/(p1+p2+4*(d1+d2));


		// Combined test

		if ( t1+p1+4*d1+t2+p2+4*d2 > 0 )
			com_chisq = ( ( (t1+p1+2*d1) - (t2+p2+2*d2) ) 
			* ( (t1+p1+2*d1) - (t2+p2+2*d2) ) ) 
			/ ( t1+p1+4*d1+t2+p2+4*d2 ) ;
//	}


	// Display asymptotic results
	cExporter *Cp_pt = cExporter::summon("parent.tdt.res");

	// Get the final p-value
	double pvalue = PVchisq(tdt_chisq,1);

	// Skip?, if filtering p-values
	if (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), pvalue))
		return;
//		if ( par::pfilter && pvalue > par::pfvalue ) 
//			continue;

	entryVariant(Cp_pt, X_var);
	Cp_pt->fmt("%g	%g", t1, t2);

	// Odds ratio for T:U
	double OR = t1 / t2;

	if (NA(OR))
		Cp_pt->put(IS_ASSIGNED(ci) ?"	NA	NA	NA":"	NA");
	else if (IS_ASSIGNED(ci)) {
		double zt = ltqnorm( 1 - (1 - OPT_REAL(ci)) / 2  ) ; 
		double OR_lower = exp( log(OR) - zt * sqrt(1/t1+1/t2)) ;
		double OR_upper = exp( log(OR) + zt * sqrt(1/t1+1/t2)) ;

		Cp_pt->fmt("	%g	%g	%g", OR, OR_lower, OR_upper);
	} else
		Cp_pt->fmt("	%g", OR);

	if (tdt_chisq>=0)
		Cp_pt->fmt("	%g	%g", tdt_chisq, PVchisq(tdt_chisq, 1));
	else
		Cp_pt->put("	NA	NA");

	Cp_pt->fmt("	%g:%g", p1+2*d1, p2+2*d2);

	if (par_chisq>=0)
		Cp_pt->fmt("	%g	%g", par_chisq, PVchisq(par_chisq,1));
	else
		Cp_pt->put("	NA	NA");

	if (com_chisq>=0)
		Cp_pt->fmt("	%g	%g\n", com_chisq, PVchisq(com_chisq,1));
	else
		Cp_pt->put("	NA	NA\n");

	///////////////////////////////////////////
	// Choose which statistic for permutation

// 	if (par::perm_TDT_basic) res[l] = tdt_chisq;
// 	else if (par::perm_TDT_parent) res[l] = par_chisq;
// 	else res[l] = com_chisq;

}

cTDTanalysis::cTDTanalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS)
	: cAnalysis(Cp_inpIO)
{
	/* --blup should not be here */
	if (IS_ASSIGNED(blup))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--blup", "--tdt");

	/* All phenotypes should be dichotomous */
	wsUint N_pheno = Cp_IO->sizePheno();
	for (wsUint i=0 ; i<N_pheno ; i++)
		if (Cp_IO->isContinuous(i))
			halt("[%d]th phenotype is continuous, TDT cannot be applied\n");

	/* Set default variables */
	Cp_anaFS	= Cp_inpAnaFS;
}

cTDTanalysis::~cTDTanalysis()
{

}

int thr_tdt(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	xAnaThread*		Xp_at	= (xAnaThread *)Vp_shareData;
	cTDTanalysis*	Cp_tdt	= (cTDTanalysis *)(Xp_at->Vp_data);
	cFamStrAnalysis	*Cp_fs	= Cp_tdt->getFamStrAna();
	mNucFam&		Xm_nf	= Cp_fs->getNucFamData();
	cIO*			Cp_IO	= Xp_at->Cp_IO;
	char**			Na_data	= Cp_IO->getGenotype();
	vVariant&		Xa_vrt	= Cp_IO->getVariant();
	int*			Np_data	= (int *)Vp_data;

	wsUint	N_pheno		= Cp_IO->sizePheno();
	wsUint	N_s			= (wsUint)Np_data[0];
	wsUint	N_e			= (wsUint)Np_data[1];

	wsMat	Ra_res		= (wsMat)Vp_result;

	for (wsUint j=N_s ; j<N_e ; j++) {
		wsReal *Ra_rets = sseVector(N_pheno*2);
		cTDTanalysis::_test(N_pheno, Cp_IO->getPhenos(),
			Na_data, Xm_nf, Xa_vrt[j], j, Ra_rets);
		for (wsUint i=0,k=0 ; i<N_pheno ; i++,k+=N_pheno) {
			Ra_res[i][j*2] = Ra_rets[k];
			Ra_res[i][j*2+1] = Ra_rets[k+1];
		}
	}

	return 0;
}

int thr_sdt(int N_idx, void *Vp_shareData, void *Vp_data, void *Vp_result)
{
	xAnaThread*		Xp_at	= (xAnaThread *)Vp_shareData;
	cSDTanalysis*	Cp_sdt	= (cSDTanalysis *)(Xp_at->Vp_data);
	cFamStrAnalysis	*Cp_fs	= Cp_sdt->getFamStrAna();
	mNucFam&		Xm_nf	= Cp_fs->getNucFamData();
	cIO*			Cp_IO	= Xp_at->Cp_IO;
	char**			Na_data	= Cp_IO->getGenotype();
	vVariant&		Xa_vrt	= Cp_IO->getVariant();
	int*			Np_data	= (int *)Vp_data;

	wsUint	N_pheno		= Cp_IO->sizePheno();
	wsUint	N_s			= (wsUint)Np_data[0];
	wsUint	N_e			= (wsUint)Np_data[1];

	wsMat	Ra_res		= (wsMat)Vp_result;

	for (wsUint j=N_s ; j<N_e ; j++) {
		wsReal *Ra_rets = sseVector(N_pheno*2);
		cSDTanalysis::_test(N_pheno, Cp_IO->getPhenos(),
			Na_data, Xm_nf, Xa_vrt[j], j, Ra_rets);
		for (wsUint i=0,k=0 ; i<N_pheno ; i++,k+=N_pheno) {
			Ra_res[i][j*2] = Ra_rets[k];
			Ra_res[i][j*2+1] = Ra_rets[k+1];
		}
	}

	return 0;
}

void cTDTanalysis::run()
{
	vVariant&	Xv_vrt	= Cp_IO->getVariant();
	wsUint		N_vrt	= (wsUint)Xv_vrt.size();
	mNucFam&	Xm_nf	= Cp_anaFS->getNucFamData();
	char**		Na_data	= Cp_IO->getGenotype();
	wsUint		N_pheno	= Cp_IO->sizePheno();
	vPheno&		Xa_pi	= Cp_IO->getPhenoInfo();

	/* Prepare output */
	cExporter	*Cp_tdt	= cExporter::summon("tdt.res");
	LOGoutput("Result of Transmission Disequilibrium Test is exported to "
		"[%s.tdt.res]\n", OPT_STRING(out));
	headerVariant(Cp_tdt);
	Cp_tdt->put("	PHENO	STAT	P_TDT\n");

	wsMat Ra_pTDT = sseMatrix(N_pheno, N_vrt);
	wsUint j=0;
	if (IS_MULTITHREAD) {
		xAnaThread	X_at	= { getIO(), this };
		wsMat		Ra_rets	= sseMatrix(N_pheno, 2*N_vrt);
		WORKER().run(thr_tdt, forVariant_equal, &X_at, Ra_rets,
			sizeof(int)*3);

		/* --adjust, --genoctrl */
		if (OPT_ENABLED(adjust) || OPT_ENABLED(genoctrl)) {
			double *Ra_stat = sseVector(N_vrt);
			for (wsUint j=0 ; j<N_pheno ; j++) {
				for (wsUint i=0 ; i<N_vrt ; i++)
					Ra_stat[i] = Ra_rets[j][i*2];

				vPheno &Xa_phe = getIO()->getPhenoInfo();
				/* Make name */
				char S_buf[512] = { 0, };
				if (N_pheno == 1)
					strcpy(S_buf, "tdt.adjust.res");
				else
					sprintf(S_buf, "tdt.%s.adjust.res", Xa_phe[j].S_name.c_str());
				multcomp(getIO(), Ra_stat, N_vrt, S_buf);
			}
			sseFree(Ra_stat);
		}

		/* Print results */
		FOREACHDO (vVariant_it, Xv_vrt, i, j++) {
			/* Print result */
			for (wsUint k=0,l=0 ; k<N_pheno ; k++,l+=2) {
				wsReal R_pTDT = Ra_rets[k][j*2 + 1];
				Ra_pTDT[k][j] = R_pTDT;

				/* --remna */
				if (OPT_ENABLED(remna) && NA(R_pTDT)) continue;
				/* --pvalrange */
				if (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), R_pTDT))
					continue;

				entryVariant(Cp_tdt, *i);
				Cp_tdt->fmt("	%s	%g	%g\n", Xa_pi[k].S_name.c_str(),
					Ra_rets[l], R_pTDT);
			}
		}
	} else {
		DMAT_t Ra_stat = NULL;
		
		/* --adjust, --genoctrl */
		if (OPT_ENABLED(adjust) || OPT_ENABLED(genoctrl))
			dblMatrix(N_pheno, N_vrt);

		wsReal *Ra_rets = sseVector(N_pheno*2);
		FOREACHDO (vVariant_it, Xv_vrt, i, j++) {
			_test(N_pheno, Cp_IO->getPhenos(),
				Na_data, Xm_nf, *i, j, Ra_rets);

			/* Print result */
			for (wsUint k=0,l=0 ; k<N_pheno ; k++,l+=2) {
				wsReal R_pTDT = Ra_rets[l+1];
				Ra_pTDT[k][j] = R_pTDT;

				/* --remna */
				if (OPT_ENABLED(remna) && NA(R_pTDT)) continue;
				/* --pvalrange */
				if (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), R_pTDT))
					continue;

				entryVariant(Cp_tdt, *i);
				Cp_tdt->fmt("	%s	%g	%g\n", Xa_pi[k].S_name.c_str(),
					Ra_rets[l], R_pTDT);

				/* --adjust, --genoctrl */
				if (Ra_stat) Ra_stat[k][j] = (double)Ra_rets[l];
			}
		}
		sseFree(Ra_rets);

		/* --adjust, --genoctrl */
		if (Ra_stat) {
			for (wsUint j=0 ; j<N_pheno ; j++) {
				vPheno &Xa_phe = getIO()->getPhenoInfo();
				/* Make name */
				char S_buf[512] = { 0, };
				if (N_pheno == 1)
					strcpy(S_buf, "tdt.adjust.res");
				else
					sprintf(S_buf, "tdt.%s.adjust.res", Xa_phe[j].S_name.c_str());
				multcomp(getIO(), Ra_stat[j], N_vrt, S_buf);
			}
			sseFree(Ra_stat);
		}
	}
	/* [[R]] Visualization */
	qqVariant("TDT analysis", "tdt", Ra_pTDT, getIO());
	sseUnmat(Ra_pTDT, N_pheno);
}

void cTDTanalysis::_test(wsUint N_pheno, wsReal **Ra_pheno, char **Na_data,
	mNucFam &Xm_nf, xVariant &X_snp, wsUint idx, wsReal *Rp_rets)
{
	wsUint P1[16] = { 0, };
	wsUint P2[16] = { 0, };

	/* For each family */
	FOREACH (mNucFam_it, Xm_nf, i) {
		xNucFamily &n = i->second;

		/* Pass if cannot calculate */
		if (!n.Xp_pat || !n.Xp_mat || n.Xp_pat->N_idx == SAMP_NODATA ||
			n.Xp_mat->N_idx == SAMP_NODATA) continue;

		/* Two indices */
		wsUint N1 = i->second.Xp_pat->N_idx;
		wsUint N2 = i->second.Xp_mat->N_idx;
		char Gp = Na_data[N1][idx];
		char Gm = Na_data[N2][idx];
		wsUint Nhet = (Gp%2) + (Gm%2);

		/* Should be both genotyped and at least one hetero */
		if (isMissing(Na_data[N1][idx]) || isMissing(Na_data[N2][idx]) ||
			Nhet == 0) continue;

		FOREACH (vSampPtr_it, n.Xp_childs, j) {
			wsUint	C = (*j)->N_idx;
			char	G = Na_data[C][idx];

			/* Skip if missing */
			if (isMissing(G)) continue;

			/* For all phenotypes */
			for (wsUint k=0 ; k<N_pheno ; k++) {
				if (Ra_pheno[k][C] != WISARD_AFFECTED) continue;

				char tA = 0, uA = 0, tB = 0, uB = 0;
				switch (G) {
				case 0:
					if (Nhet == 2) { tA = 1; uA = 2; tB = 1; uB = 2; }
					else { tA = 1; uA = 2; }
					break;
				case 1:
					if (Gp%2) {
						if (Gm%2) { tA = 1; tB = 2; uA = 2; uB =1; }
						else if (!Gm) { tA = 2; uA = 1; }
						else { tA = 1; uA = 2; }
					} else if (!Gp) { tA = 2; uA = 1; }
					else { tA = 1; uA = 2; }
					break;
				case 2:
					if (Nhet) { tA = 2; uA = 1; tB = 2; uB = 1; }
					else { tA = 2; uA = 1; }
					break;
				default: halt("Errornous case");
				}

				if (tA == 1) P1[k]++;
				if (tB == 1) P1[k]++;
				if (tA == 2) P2[k]++;
				if (tB == 2) P2[k]++;
			}
		}
	}

	for (wsUint i=0,j=0 ; i<N_pheno ; i++) {
		wsReal R_tTDT = WISARD_NAN;
		if (P1[i] + P2[i] > 0)
			R_tTDT = SQR(P1[i]-P2[i]) / (P1[i]+P2[i]);
		Rp_rets[j++] = R_tTDT;
		Rp_rets[j++] = PVchisq(R_tTDT, 1.0);
	}
}

cSDTanalysis::cSDTanalysis(cIO *Cp_inpIO, cFamStrAnalysis *Cp_inpAnaFS)
	: cAnalysis(Cp_inpIO)
{
	/* --blup should not be here */
	if (IS_ASSIGNED(blup))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--blup", "--tdt");

	/* All phenotypes should be dichotomous */
	wsUint N_pheno = Cp_IO->sizePheno();
	for (wsUint i=0 ; i<N_pheno ; i++)
		if (Cp_IO->isContinuous(i))
			halt("[%d]th phenotype is continuous, TDT cannot be applied\n");

	/* Set default variables */
	Cp_anaFS = Cp_inpAnaFS;
}

cSDTanalysis::~cSDTanalysis()
{

}

void cSDTanalysis::run()
{
	vVariant&	Xv_vrt	= Cp_IO->getVariant();
	wsUint		N_vrt	= (wsUint)Xv_vrt.size();
	mNucFam&	Xm_nf	= Cp_anaFS->getNucFamData();
	char**		Na_data	= Cp_IO->getGenotype();
	wsUint		N_pheno	= Cp_IO->sizePheno();
	vPheno&		Xa_pi	= Cp_IO->getPhenoInfo();
	cExporter	*Cp_tdt	= cExporter::summon("sdt.res");
	LOGoutput("Result of Sibship Disequilibrium Test is exported to "
		"[%s.sdt.res]\n", OPT_STRING(out));
	headerVariant(Cp_tdt);
	if (N_pheno > 1)
		Cp_tdt->put("	PHENO	STAT	P_SDT\n");
	else
		Cp_tdt->put("	STAT	P_SDT\n");

	wsUint j=0;
	if (IS_MULTITHREAD) {
		xAnaThread	X_at	= { getIO(), this };
		wsMat		Ra_rets	= sseMatrix(N_pheno, 2*N_vrt);
		WORKER().run(thr_sdt, forVariant_equal, &X_at, Ra_rets,
			sizeof(int)*3);

		/* --adjust, --genoctrl */
		if (OPT_ENABLED(adjust) || OPT_ENABLED(genoctrl)) {
			double *Ra_stat = sseVector(N_vrt);
			for (wsUint j=0 ; j<N_pheno ; j++) {
				for (wsUint i=0 ; i<N_vrt ; i++)
					Ra_stat[i] = Ra_rets[j][i*2];

				vPheno &Xa_phe = getIO()->getPhenoInfo();
				/* Make name */
				char S_buf[512] = { 0, };
				if (N_pheno == 1)
					strcpy(S_buf, "sdt.adjust.res");
				else
					sprintf(S_buf, "sdt.%s.adjust.res", Xa_phe[j].S_name.c_str());
				multcomp(getIO(), Ra_stat, N_vrt, S_buf);
			}
			sseFree(Ra_stat);
		}

		/* Print results */
		FOREACHDO (vVariant_it, Xv_vrt, i, j++) {
			/* Print result */
			for (wsUint k=0,l=0 ; k<N_pheno ; k++,l+=2) {
				wsReal R_pSDT = Ra_rets[k][j*2 + 1];

				/* --remna */
				if (OPT_ENABLED(remna) && NA(R_pSDT)) continue;
				/* --pvalrange */
				if (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), R_pSDT))
					continue;

				entryVariant(Cp_tdt, *i);
				if (N_pheno > 1) Cp_tdt->fmt("	%s", Xa_pi[k].S_name.c_str());
				if (NA(R_pSDT)) {
					if (!OPT_ENABLED(remna))
						Cp_tdt->put("	NA	NA\n");
				} else 
					Cp_tdt->fmt("	%g	%g\n", Ra_rets[l], R_pSDT);
			}
		}
	} else {
		DMAT_t Ra_stat = NULL;

		/* --adjust, --genoctrl */
		if (OPT_ENABLED(adjust) || OPT_ENABLED(genoctrl))
			dblMatrix(N_pheno, N_vrt);

		wsReal *Ra_rets = sseVector(N_pheno*2);
		FOREACHDO (vVariant_it, Xv_vrt, i, j++) {
			_test(N_pheno, Cp_IO->getPhenos(),
				Na_data, Xm_nf, *i, j, Ra_rets);

			/* Print result */
			for (wsUint k=0,l=0 ; k<N_pheno ; k++,l+=2) {
				wsReal R_pSDT = Ra_rets[l+1];

				/* --remna */
				if (OPT_ENABLED(remna) && NA(R_pSDT)) continue;
				/* --pvalrange */
				if (IS_ASSIGNED(pvalrange) && !isInRange(OPT_RANGE(pvalrange), R_pSDT))
					continue;

				entryVariant(Cp_tdt, *i);
				if (N_pheno > 1) Cp_tdt->fmt("	%s", Xa_pi[k].S_name.c_str());
				if (NA(R_pSDT)) {
					if (!OPT_ENABLED(remna))
						Cp_tdt->put("	NA	NA\n");
				} else
					Cp_tdt->fmt("	%g	%g\n", Ra_rets[l], R_pSDT);
			}
		}
		sseFree(Ra_rets);

		/* --adjust, --genoctrl */
		if (Ra_stat) {
			for (wsUint j=0 ; j<N_pheno ; j++) {
				vPheno &Xa_phe = getIO()->getPhenoInfo();
				/* Make name */
				char S_buf[512] = { 0, };
				if (N_pheno == 1)
					strcpy(S_buf, "sdt.adjust.res");
				else
					sprintf(S_buf, "sdt.%s.adjust.res", Xa_phe[j].S_name.c_str());
				multcomp(getIO(), Ra_stat[j], N_vrt, S_buf);
			}
			sseFree(Ra_stat);
		}
	}
}

void cSDTanalysis::_test(wsUint N_pheno, wsReal **Ra_pheno, char **Na_data,
						 mNucFam &Xm_nf, xVariant &X_snp, wsUint idx, wsReal *Rp_rets)
{
// 	wsUint P1[16] = { 0, };
// 	wsUint P2[16] = { 0, };

	wsUint	*mA = NULL, *nA = NULL;
	wsUint	*mU = NULL, *nU = NULL;
	wsAlloc(mA, wsUint, N_pheno);
	wsAlloc(mU, wsUint, N_pheno);
	wsAlloc(nA, wsUint, N_pheno);
	wsAlloc(nU, wsUint, N_pheno);

	wsUint *Shit = NULL;
	wsUint *Smiss = NULL;
	wsCalloc(Shit, wsUint, N_pheno);
	wsCalloc(Smiss, wsUint, N_pheno);

	/* For each family */
	FOREACH (mNucFam_it, Xm_nf, i) {
		xNucFamily &n = i->second;

		memset(mA, 0x00, sizeof(wsUint)*N_pheno);
		memset(mU, 0x00, sizeof(wsUint)*N_pheno);
		memset(nA, 0x00, sizeof(wsUint)*N_pheno);
		memset(nU, 0x00, sizeof(wsUint)*N_pheno);

		/* Only sees childs */
		FOREACH (vSampPtr_it, n.Xp_childs, j) {
			wsUint	C = (*j)->N_idx;
			char	G = Na_data[C][idx];

			/* Skip if missing */
			if (isMissing(G)) continue;

			/* For all phenotypes */
			for (wsUint k=0 ; k<N_pheno ; k++) {
				Ra_pheno[k][C] == WISARD_AFFECTED ? nA[k]++ : nU[k]++;
				if (G)
					Ra_pheno[k][C] == WISARD_AFFECTED ? mA[k]++ : mU[k]++;
			}
		}

		/* Get score */
		for (wsUint k=0 ; k<N_pheno ; k++) {
			wsReal rA = (wsReal)mA[k] / (wsReal)nA[k];
			wsReal rU = (wsReal)mU[k] / (wsReal)nU[k];
			if (rA > rU) Shit[k]++;
			else if (rA < rU) Smiss[k]++;
		}
	}

	for (wsUint i=0,j=0 ; i<N_pheno ; i++) {
		wsReal R_tSDT = WISARD_NAN;
		if (Shit[i] + Smiss[i] > 0)
			R_tSDT = SQR(Shit[i]-Smiss[i]) / (Shit[i]+Smiss[i]);
		Rp_rets[j++] = R_tSDT;
		if (R_tSDT != R_tSDT)
			Rp_rets[j++] = WISARD_NAN;
		else
			Rp_rets[j++] = pnorm(R_tSDT, 0.0, 1.0, 0) * 2.0;
	}

	DEALLOC(Shit);
	DEALLOC(Smiss);
	DEALLOC(mA);
	DEALLOC(mU);
	DEALLOC(nA);
	DEALLOC(nU);
}

wsUint* cFamStrAnalysis::getFamIndices()
{
	return Na_famIndices;
}

} // End namespace ONETOOL
