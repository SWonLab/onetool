#include "global/common.h"
#include "utils/util.h"
#include "utils/log.h"
#include "global/Rconn.h"
#include "global/option.h"

namespace ONETOOL {

/*
 * s	string
 * u	unsigned int
 * x	structure
 * 
 */
xRepData Xa_rd[] = {
	{ WISARD_REP_NONE,        "REP_NONE",        '\0', 0, { 0, 0 }, 0, },
	{ WISARD_REP_DATESTART,   "REP_DATESTART",   't',  0, { 0, 0 }, 0, },
	{ WISARD_REP_DATEEND,     "REP_DATEEND",     't',  0, { 0, 0 }, 0, },
	{ WISARD_REP_TOTEXECTIME, "REP_TOTEXECTIME", 's',  0, { 0, 0 }, 0, },
	{ WISARD_REP_RANDSEED,    "REP_RANDSEED",    'u',  0, { 0, 0 }, 0, },
	{ WISARD_REP_CACT,        "REP_CACT",        's',  0, { 0, 0 }, 0, },
	{ WISARD_REP_CORSTR,      "REP_CORSTR",      's',  0, { 0, 0 }, 0, },
	{ WISARD_REP_OUTPREFIX,   "REP_OUTPREFIX",   's',  0, { 0, 0 }, 0, },
	{ WISARD_REP_WARNING,     "REP_WARNING",     'x',  0, { 0, 0 }, 0, },
	{ WISARD_REP_EXECRES,     "REP_EXECRES",     's',  0, { 0, 0 }, 0, },
	{ WISARD_REP_TOTNCORE,    "REP_TOTNCORE",    'u',  0, { 0, 0 }, 0, },
	{ WISARD_REP_USEDNCORE,   "REP_USEDNCORE",   'u',  0, { 0, 0 }, 0, },
};

void report(xRepType X_rtype, wsStrCst S_val)
{
	/* Type MUST BE 's', otherwise ERROR */
	xRepData &X_d = Xa_rd[X_rtype];
	if (X_d.S_fmt != 's') halt_fmt(WISARD_SYST_INVL_REPTYPE, X_d.name,
		"string", S_val);

	/* Already assigned */
	if (X_d.Sp_val) {
		/* If the entry is NOT replaceable ERROR */
		if (!X_d.B_replaceable) halt_fmt(WISARD_SYST_CANT_REPREPLACE,
			X_d.name, X_d.Sp_val);
		/* Otherwise deallocate */
		DEALLOC(X_d.Sp_val);
	}

	/* Allocate value */
	wsAlloc(X_d.Sp_val, char, strlen(S_val)+1);
	strcpy(X_d.Sp_val, S_val);

	/* Update time */
	gettimeofday(&(X_d.X_time), NULL);
}

void report(xRepType X_rtype, wsUint N_val)
{
	/* Make print */
	char *Sp_tmp = NULL;
	wsAlloc(Sp_tmp, char, 64);
	sprintf(Sp_tmp, "%u", N_val);

	/* Type MUST BE 'u' or 'i', otherwise ERROR */
	xRepData &X_d = Xa_rd[X_rtype];
	if (X_d.S_fmt != 'u' && X_d.S_fmt != 'i')
		halt_fmt(WISARD_SYST_INVL_REPTYPE, X_d.name, "integer", Sp_tmp);

	/* Already assigned */
	if (X_d.Sp_val) {
		/* If the entry is NOT replaceable ERROR */
		if (!X_d.B_replaceable) halt_fmt(WISARD_SYST_CANT_REPREPLACE,
			X_d.name, X_d.Sp_val);
		/* Otherwise deallocate */
		DEALLOC(X_d.Sp_val);
	}

	/* Allocate value (up to 64 digits) */
	X_d.Sp_val = Sp_tmp;

	/* Update time */
	gettimeofday(&(X_d.X_time), NULL);
}

void report(xRepType X_rtype)
{
	/* Type MUST BE 't', otherwise ERROR */
	xRepData &X_d = Xa_rd[X_rtype];
	if (X_d.S_fmt != 't')
		halt_fmt(WISARD_SYST_INVL_REPTYPE, X_d.name, "time", "");

	/* Already assigned */
	if (X_d.Sp_val) {
		/* If the entry is NOT replaceable ERROR */
		if (!X_d.B_replaceable) halt_fmt(WISARD_SYST_CANT_REPREPLACE,
			X_d.name, X_d.Sp_val);
		/* Otherwise deallocate */
		DEALLOC(X_d.Sp_val);
	}

	/* Allocate value (up to 64 digits) */
	X_d.Sp_val = NULL;

	/* Update time */
	gettimeofday(&(X_d.X_time), NULL);
}

void reportEnd()
{
	for (wsUint i=0 ; i<len(Xa_rd, xRepData) ; i++) {
		if (Xa_rd[i].Sp_val) free(Xa_rd[i].Sp_val);
	}
#ifdef USE_R
	vStr	Xa_subj;
	vStr	Xa_valu;

	/* Build report vector */

	/* Send to R */
	R().sendVector("rnames", Xa_subj);
	R().sendVector("rvals", Xa_valu);

	/* Draw */
//	R().Rparse("report('%s.report.pdf', rnames, rvalus)", OPT_STRING(out));
#endif
}

} // End namespace ONETOOL
