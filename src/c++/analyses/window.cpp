#include <sstream>
#include "analyses/window.h"

#define WISARD_MAX_SZWINDOW 1000000000

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

cWindowAnalysis::cWindowAnalysis(cIO *Cp_IO) : cAnalysis(Cp_IO)
{
	vVariant &Xa_snp = Cp_IO->getVariant();

	/* Check out the required parameters */
	if (!IS_ASSIGNED(density) && !IS_ASSIGNED(tstv))
		halt("--density or --tstv required but not assigned");

	/* Check the range of the value */
	if (IS_ASSIGNED(density))
		N_szWindow = OPT_NUMBER(density);
	else
		N_szWindow = OPT_NUMBER(tstv);
	if (N_szWindow > WISARD_MAX_SZWINDOW)
		halt("Too large window size[%d] given, should not over [%d]",
			N_szWindow, WISARD_MAX_SZWINDOW);

	/* Make the bins */

	/* Step 1 : Find max value for each chromosome */
	wsUint I = 0;
	sseCalloc(Na_szChr, wsUint, MAX_NCHR);
	FOREACHDO (vVariant_it, Xa_snp, i, I++) {
		if (i->chr <= 0 || (wsUint)i->chr > NCHR_SPECIES) continue;

		if (Na_szChr[i->chr - 1] < i->pos)
			Na_szChr[i->chr - 1] = i->pos;
	}

	/* Step 2 : Make bin */
	if (IS_ASSIGNED(tstv))
		wsAlloc(Na_trs, wsUint*, MAX_NCHR);
	wsAlloc(Na_bin, wsUint*, MAX_NCHR);
	sseMalloc(Na_szBin, wsUint, MAX_NCHR);
	for (wsUint i=0 ; i<NCHR_SPECIES ; i++) {
		Na_szBin[i] = (wsUint)((wsReal)Na_szChr[i]/
			(wsReal)N_szWindow+REAL_CONST(0.5)) + 1;
		if (IS_ASSIGNED(tstv)) {
			sseCalloc(Na_trs[i], wsUint, Na_szBin[i]);
		}
		sseCalloc(Na_bin[i], wsUint, Na_szBin[i]);
	}
}

cWindowAnalysis::~cWindowAnalysis()
{
	for (wsUint i=0 ; i<NCHR_SPECIES ; i++)
		sseFree(Na_bin[i]);
	DEALLOC(Na_bin);
	sseFree(Na_szBin);
	sseFree(Na_szChr);
}

void cWindowAnalysis::run()
{
	vVariant &Xa_snp = Cp_IO->getVariant();
	wsUint*	Na_cbin	= NULL;
	wsUint*	Na_tbin	= NULL;
	wsCalloc(Na_cbin, wsUint, MAX_NCHR);
	wsCalloc(Na_tbin, wsUint, MAX_NCHR);

	/* Step 3 : Fill bin by the position */
	wsUint I = 0;
	FOREACHDO (vVariant_it, Xa_snp, i, I++) {
		/* Only account for 'valid' chromosomes */
		if (i->chr <= 0 || (wsUint)i->chr > NCHR_SPECIES) continue;

		wsUint N_bin = (wsUint)((wsReal)i->pos / (wsReal)N_szWindow);
		if (N_bin >= Na_szBin[i->chr - 1])
			halt("Invalid range");
		Na_bin[i->chr-1][N_bin]++;
		Na_cbin[i->chr-1]++;

		if (i->indel1 && i->indel2) {
			/* In case of indel,
			 * it only possible to calculate ts/tv ratio unless
			 * its length are all 1, and none of them is '-'
			 */
		} else {
			if ((i->al1 == 'A' && i->al2 == 'G') ||
				(i->al1 == 'G' && i->al2 == 'A') ||
				(i->al1 == 'C' && i->al2 == 'T') ||
				(i->al1 == 'T' && i->al2 == 'C'))
			Na_trs[i->chr-1][N_bin]++;
			Na_tbin[i->chr-1]++;
		}
		//		if (Na_szChr[i->chr - 1] < i->pos) 
		//			Na_szChr[i->chr - 1] = i->pos;
	}

	/* Step 4 : Export result */
/**/cTableExporter*	Cp_tstv		= NULL;
/**/cTableExporter*	Cp_ctstv	= NULL;
	cTableExporter	C_dens("density.variant.lst", "siiis",
		"Variant density by physical position", 0, 5,
		"CHR", "SZBIN", "SZCHR", "NBIN", "DENSITY");
	if (IS_ASSIGNED(tstv)) {
		Cp_tstv		= new cTableExporter("tstv.variant.lst", "siiis",
			"Ts/Tv ratio by physical position", 0, 5,
			"CHR", "SZBIN", "SZCHR", "NBIN", "TSTV");
		Cp_ctstv		= new cTableExporter("tstv.chr.lst", "siiir",
			"Ts/Tv ratio for chromosome", 0, 5,
			"CHR", "SZBIN", "SZCHR", "NBIN", "TSTV");
	}
	for (wsUint i=0 ; i<NCHR_SPECIES ; i++) {
		wsStrCst Sp_chr = getChrName2(i+1);

		wsString S_dens;
		wsString S_tstv;
		for (wsUint j=0 ; j<Na_szBin[i] ; j++) {
			if (j) S_dens += ',';
			S_dens += Na_bin[i][j];

			if (IS_ASSIGNED(tstv)) {
				if (j) S_tstv += ',';
				S_tstv += (wsReal)(Na_bin[i][j]-Na_trs[i][j]) /
					(wsReal)(Na_trs[i][j]);
			}
		}
		C_dens.write(5, Sp_chr, N_szWindow, Na_szChr[i], Na_szBin[i],
			S_dens.get());

		if (IS_ASSIGNED(tstv)) {
			Cp_tstv->write(5, Sp_chr, N_szWindow, Na_szChr[i], Na_szBin[i],
				S_tstv.get());
			Cp_ctstv->write(5, Sp_chr, N_szWindow, Na_szChr[i], Na_szBin[i],
				(wsReal)(Na_cbin[i]-Na_tbin[i])/(wsReal)Na_tbin[i]);
		}
	}
	DEALLOC(Na_cbin);
	DEALLOC(Na_tbin);
	if (Cp_tstv) {
		delete Cp_tstv;
		delete Cp_ctstv;
	}
}

#endif

} // End namespace ONETOOL
