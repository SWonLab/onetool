#include "analyses/annot.h"
#include "input/stream.h"

namespace ONETOOL {

cAnnotAnalysis *Cp_anno = NULL;

cAnnotAnalysis* wsAnnot(cIO *Cp_IO)
{
	if (Cp_anno == NULL) {
		if (Cp_IO == NULL) halt_fmt(WISARD_SYST_NULL_IO);
		Cp_anno = new cAnnotAnalysis(Cp_IO);
		Cp_anno->run();
	}
	return Cp_anno;
}

void ANNOclear()
{
	if (Cp_anno != NULL)
		delete Cp_anno;
}

void cAnnotAnalysis::_loadAnnoGene(wsStrCst Sp_gene)
{
	char *S_buf = NULL;
	wsAlloc(S_buf, char, 65536);

	/* Open stream */
	cStrFile	C_gene(Sp_gene, "Gene annotation file");

	/* Read entries */
	wsUint L = 1;
	for (char *a=NULL,*b=NULL ; C_gene.gets(S_buf, 2048) ; L++) {
		xAnnoGene X_anno = { 0, };

		getString(&S_buf, &a); // S_buf == chr
		if (a == NULL) halt_fmt(WISARD_INVL_ANNOGENE_FORMAT, L, 4, 1);

		/* Retrieve chromosome */
		int N_chr = getChr(S_buf);

		getString(&a, &b); // a == start
		if (b == NULL) halt_fmt(WISARD_INVL_ANNOGENE_FORMAT, L, 4, 2);

		/* Retrieve start position */
		X_anno.N_start = atoi(a);

		getString(&b, &a); // b == end
		if (a == NULL) halt_fmt(WISARD_INVL_ANNOGENE_FORMAT, L, 4, 3);

		/* Retrieve end position */
		X_anno.N_end = atoi(b);

		// a == genename
		// trim a
		char *x = a + strlen(a) - 1;
		while (*x != '\0' &&
			(*x == '\r' || *x == '\n' || *x == '\t' || *x == ' '))
			*(x--) = '\0';
		/* Check name exists */
		if (a[0] == '\0') halt_fmt(WISARD_NULL_ANNOGENENAME, L);

		wsAlloc(X_anno.S_name, char, strlen(a)+1);
		strcpy(X_anno.S_name, a);

		/* Insert */
		Xa_annoGene[N_chr].push_back(X_anno);
	}
	LOGnote("[%d] annotation entries loaded\n", L);

	/* Close */
	DEALLOC(S_buf);
}	

cAnnotAnalysis::cAnnotAnalysis(cIO *Cp_IO) : cAnalysis(Cp_IO)
{
	wsAlloc(Xa_annoGene, vector<xAnnoGene>, MAX_NCHR+1);

	if (IS_ASSIGNED(annogene))
		_loadAnnoGene(OPT_STRING(annogene));
 	if (IS_ASSIGNED(annosnp))
		halt_fmt(WISARD_SYST_NULL_IMPLEMENT, "SNP-level annotation");
}

cAnnotAnalysis::~cAnnotAnalysis()
{
	/* Free up memories */
	for (wsUint i=0 ; i<NCHR_SPECIES ; i++) {
		for (vector<xAnnoGene>::iterator j=Xa_annoGene[i].begin() ;
			j!=Xa_annoGene[i].end() ; j++)
			DEALLOC(j->S_name);
	}
	DEALLOC(Xa_annoGene);
}

void cAnnotAnalysis::run()
{
	//halt_fmt(WISARD_SYST_CANT_RUNANALYSIS);
	vVariant& Xv_var = getIO()->getVariant();
	FOREACH (vVariant_it, Xv_var, i)
		anno(*i);
}

void _printDist(string &S_prt, int N_dist, int N_isPlus)
{
	S_prt.append(N_isPlus?"+":"-");

	char S_dist[256];
	if (N_dist < 1000) {
		sprintf(S_dist, "%d", N_dist);
		S_prt.append(S_dist);
		S_prt.append("bp");
	} else if (N_dist < 1000000) {
		sprintf(S_dist, "%d", (int)(N_dist/REAL_CONST(1000.0) + REAL_CONST(0.5)));
		S_prt.append(S_dist);
		S_prt.append("Kbp");
	} else {
		sprintf(S_dist, "%d", (int)(N_dist/REAL_CONST(1000000.0) + REAL_CONST(0.5)));
		S_prt.append(S_dist);
		S_prt.append("Mbp");
	}
}

void cAnnotAnalysis::anno(xVariant &X_snp)
{
	wsStrCst	S_div	= ",";
	string	S_str;
	S_str = "<NA>";
	wsUint	N_anno	= 0;

	/* Check range of SNP to annotate */
	if (X_snp.chr <= 0 || (wsUint)X_snp.chr > NCHR_SPECIES) {
		X_snp.anno = strdup(S_str.c_str());
		return;
	}

	/* Find appropriate annotation */
	vector<xAnnoGene> &Xa_anno = Xa_annoGene[X_snp.chr];
	if (IS_ASSIGNED(annorange)) {
		wsUint N_annoRange = OPT_NUMBER(annorange);

		FOREACH (vector<xAnnoGene>::iterator, Xa_anno, i) {
			int N_s = i->N_start < N_annoRange ? 0 : i->N_start-N_annoRange;
			/* Within gene */
			if (i->N_start <= X_snp.pos && X_snp.pos <= i->N_end) {
				if (N_anno)
					S_str += S_div;
				else
					S_str.clear();
				S_str += i->S_name;
				N_anno++;
			}

			/* Have proximity */
			else if (N_s <= (int)X_snp.pos &&
				X_snp.pos <= i->N_end+N_annoRange) {
				if (N_anno)
					S_str += S_div;
				else
					S_str.clear();
				S_str += i->S_name;
				N_anno++;

				int N_distFromStart = (int)(i->N_start - X_snp.pos);
				if (N_distFromStart > 0) {
					S_str += "(";
					_printDist(S_str, N_distFromStart, 0);
					S_str += ")";
				} else {
					int N_distFromEnd = X_snp.pos - i->N_end;
					S_str += "(";
					_printDist(S_str, N_distFromEnd, 1);
					S_str += ")";
				}
			}
		}
	} else {
		FOREACH (vector<xAnnoGene>::iterator, Xa_anno, i)
			/* Within gene */
			if (i->N_start <= X_snp.pos && X_snp.pos <= i->N_end) {
				if (N_anno)
					S_str += S_div;
				else
					S_str.clear();
				S_str += i->S_name;
				N_anno++;
			}
	}
	X_snp.anno = strdup(S_str.c_str());
}

} // End namespace ONETOOL
