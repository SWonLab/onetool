#include "analyses/ldprune.h"
#include "utils/matrix.h"
#include "utils/memmgr.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cSymMatrix cLdPruneAnalysis::calcSetCovarianceMatrix(vInt Nv_idx)
{
	wsUint		N_var	= (wsUint)Nv_idx.size();
	wsUint		N_samp	= Cp_IO->sizeSample();
	const char*	Ba_fnd	= Cp_IO->getIsFounder();
	wsUint		N_fnd	= Cp_IO->sizeFounder();
	vSampPtr&	Xv_samp	= Cp_IO->getSample();
	char**		Na_geno	= Cp_IO->getGenotype();
	wsSym		Ra_ret	= sseSymMat(N_var);
	char**		Na_curGeno = new char*[N_var];
	LOOP (i, N_var) {
		Na_curGeno[i] = new char[N_fnd];
		wsUint I = Nv_idx[i];
		for (wsUint j=0,J=0 ; j<N_samp ; j++) if (Ba_fnd[j]) Na_curGeno[i][J++] = Na_geno[j][I];
	}
	LOOP (i, N_var) for (wsUint j=0 ; j<=i ; j++) {
		char* Na_i = Na_curGeno[i];
		char* Na_j = Na_curGeno[j];
		wsUint N_xy = 0;
		wsUint N_x = 0;
		wsUint N_y = 0;
		wsUint N_f = 0;
		LOOP (k, N_fnd) {
			if (Na_i[k] == WISARD_NA || Na_j[k] == WISARD_NA) continue;
			N_xy += Na_i[k]*Na_j[k];
			N_x += Na_i[k];
			N_y += Na_j[k];
			N_f++;
		}
		Ra_ret[i][j] = N_f>1 ? ((wsReal)N_xy - (N_x*N_y)/(wsReal)N_f) / wsReal(N_f) : W0;
	}
	LOOP (i, N_var) delete[] Na_curGeno[i];
	delete[] Na_curGeno;

	return cSymMatrix(Ra_ret, N_var);
}

char* cLdPruneAnalysis::vif_prune(cSymMatrix& m, double threshold, vInt& Nv_idxVrt)
{
	xMaf* Xa_maf = Cp_IO->getMAF();
	// Number of variables
	wsUint p = m.row();
	char* cur = new char[p];
	memset(cur, 0x01, sizeof(char)*p);

	// This only is needed if we have 2+ SNPs
	if (p < 2) return cur;
	wsSym mp = m.get();
	wsSym rp = sseSymMatP(p, mp);
	cSymMatrix r(rp, p);

	// Make 'm' a correlation matrix
	LOOP (i, p) for (wsUint j=0 ; j<=i ; j++) rp[i][j] = mp[i][j] / sqrt(mp[i][i] * mp[j][j]);
#if 0
	for (int i = 0; i < p; i++) {
		for (int j = 0; j <=i; j++)
			printf("%g ", rp[i][j]);
		printf("\n");
	}
	printf("\n");
#endif

	// Number of excluded items
	int it = 0;

	// Any SNPs with zero variance should be automatically excluded
	LOOP (i, p)
		if (rp[i][i] == 0 || NA(rp[i][i])) {
			cur[i] = 0;
			it++;
		}


	// For any pair of perfectly correlated SNPs, exclude 1
	wsReal R_pruneThrR2 = IS_ASSIGNED(prunevif) ? W1 - 1e-6 : R_pruneThr;
//	LOG("Pruning threshold [%g]\n", R_pruneThrR2);
	while (1) {
		bool done = true;
		for (wsUint i = 1; i < p ; i++) {
			if (!cur[i]) continue;
			for (wsUint j = 0 ; j < i; j++) {
				if (!cur[j]) continue;
				if (fabs(rp[i][j]) <= R_pruneThrR2) continue;
				if (true) {
					// Remove SNP with lower MAF
					if (Xa_maf[Nv_idxVrt[i]].R_maf < Xa_maf[Nv_idxVrt[j]].R_maf)
						cur[i] = 0;
					else
						cur[j] = 0;
				} else {
					// Just remove first
					cur[i] = 0;
				}
				it++;
				done = false;
				break;
			}
		}
		if (done) break;
	}

	// Skip VIF calculation?
	if (IS_ASSIGNED(prunepw)) return cur;

	// Calculate r^2 for each element versus all others
	// considering only the current non-pruned elements
	while (1) {
		// Build correlation matrix all included items
		cSymMatrix& u = r.subset(cur, 1);

		// Check enough markers left  
		if (u.row() < 2) {
			delete& u;
			break;
		}

		// Get inverse
		cSymMatrix& ui = u.ginv();
		if (!ui.row()) {
		//	LOG("Failed to invert\n");
			delete& ui;
			break;
		}

		// Calculate VIFs
		double maxVIF = 0;
		int maxI;
		int cnt = 0;
		LOOP (i, p)
			if (cur[i]) {
				// r^2 = 1 - 1/x where x is diagonal element of inverted
				// correlation matrix
				// As VIF = 1 / ( 1 - r^2 ) , implies VIF = x
				double vif = ui.get()[cnt][cnt];

				if (maxVIF < vif) {
					maxVIF = vif;
					maxI = i;
				}

				cnt++;
			}
		delete& u;
		delete& ui;

		// How are we doing?
		//LOG("maxvif %g>%g index %d\n", maxVIF, threshold, maxI);
		if (maxVIF > threshold) {
			// exclude this item
			cur[maxI] = 0;
		} else break;

		// Increase count of removed items
		it++;

		// Down to a single item or worse?
		if (it == p - 1) break;
	}

	return cur;
}

cLdPruneAnalysis::cLdPruneAnalysis(cIO *Cp_inpIO) : cAnalysis(Cp_inpIO)
{
	if (IS_ASSIGNED(prunevif) && IS_ASSIGNED(prunepw))
		halt_fmt(WISARD_CANT_EXCL_OPT, "--prunevif", "--prunepw");
	if (IS_ASSIGNED(prunevif)) {
		wsUint N_param = 0;
		char** Sa_ret = loadStringValues(OPT_STRING(prunevif), &N_param);
		if (N_param != 3) halt("Argument [%s] is not valid for --prunevif", OPT_STRING(prunevif));

		N_pruneLdWin = (wsUint)atoi(Sa_ret[0]);
		N_pruneStep = (wsUint)atoi(Sa_ret[1]);
		R_pruneThr = (wsReal)atof(Sa_ret[2]);
	} else {
		wsUint N_param = 0;
		char** Sa_ret = loadStringValues(OPT_STRING(prunepw), &N_param);
		if (N_param != 3) halt("Argument [%s] is not valid for --prunepw", OPT_STRING(prunepw));

		N_pruneLdWin = (wsUint)atoi(Sa_ret[0]);
		N_pruneStep = (wsUint)atoi(Sa_ret[1]);
		R_pruneThr = (wsReal)sqrt(atof(Sa_ret[2]));
	}

	vVariant& Xv_snp = Cp_IO->getVariant();
	wsUint	j = 0;
	FOREACHDO (vVariant_it, Xv_snp, i, j++) {
		if (i->chr <= 0 || (wsUint)i->chr > NCHR_SPECIES) continue;
		Xm_mmap[i->chr].push_back(xUintSort(j, i->pos));
	}
}

cLdPruneAnalysis::~cLdPruneAnalysis()
{

}

void cLdPruneAnalysis::run()
{
	LOG("Performing LD-based pruning...\n");

	vVariant&	Xv_vrt		= Cp_IO->getVariant(); 
	cExporter*	Cp_prIn		= cExporter::summon("prune.in.lst");
	cExporter*	Cp_prOut	= cExporter::summon("prune.out.lst");

	LOGoutput("Writing pruned-in SNPs to [%s.prune.in.lst]\n", OPT_STRING(out));
	LOGoutput("Writing pruned-out SNPs to [%s.prune.out.lst]\n", OPT_STRING(out));

	////////////////////////
	// Scan each chromosome

	// Inclusion or no?
	vBool include(Xv_vrt.size(), true);

	// Only consider founders (set flag)

	// Scan each chromosome
	for (wsUint N_chr=1 ; N_chr<=NCHR_SPECIES ; N_chr++) {
		// Skip chromosome 0

		wsUint Ns = 0;
		wsUint Ne = (wsUint)Xm_mmap[N_chr].size();

		for (wsUint s1=Ns ; s1<Ne ; s1+=N_pruneStep) {
			wsUint s2 = s1 + N_pruneLdWin - 1;

			// MW fix
			if (s2 >= Ne) s2 = Ne-1;

			// calc VIF and set 
			vector<int> nSNP(0);
			for (wsUint l=s1 ; l<=s2 ; l++) if (include[l]) nSNP.push_back(l);

			// Skip if we only have a single SNP left
			if (nSNP.size() < 2) continue;

			cSymMatrix variance;
			if (OPT_ENABLED(verbose))
				pverbose("Pruning SNPs [%d] to [%d] of [%d]\r", s1 - Ns + 1, s2 - Ns + 1, Ne - Ns + 1);

			// Calculate covariance matrices
			variance = calcSetCovarianceMatrix(nSNP);

			// Calculate VIFs
			char* cur = vif_prune(variance, R_pruneThr, nSNP);

			// Update main list
			int k = 0;
			for (wsUint l = s1; l <= s2; l++) {
				// Update main list, but do not get back
				// already excluded SNPs

				if (include[l] && !cur[k++])
					include[l] = false;
			}

			delete[] cur;
		} // next window
	} // next chromosome

	// Record what is in, what is out
	wsUint cnt_in = 0, cnt_out = 0;
	LOOP (i, Xv_vrt.size()) {
		if (include[i]) {
			Cp_prIn->fmt("%s\n", Xv_vrt[i].name);
			cnt_in++;
		} else {
			Cp_prOut->fmt("%s\n", Xv_vrt[i].name);
			cnt_out++;
		}
	}
	LOG("[%d] variants pruned out and [%d] remaining\n", cnt_out, cnt_in);
}

#endif

} // End namespace ONETOOL
