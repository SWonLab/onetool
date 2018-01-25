#include <stdlib.h>
#include <stdarg.h>
#include "global/common.h"
#include "global/option.h"
#include "global/interactive.h"
#include "utils/util.h"
#include "global/io.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

typedef struct _xItrtData_EXPL
{
	cIO			*Cp_IO;
	mDataIdx	&Xm_vrt2idx;
	size_t		N_sampName;		// Maximum length of all sample IID
	size_t		N_vrtName;		// Maximum length of all variant identifier
} xItrtData_EXPL;

union _answer
{
	wsUint	N_sel;
	char	S_sel[256];
	char	B_yes;
} ANS;

typedef enum _askType
{
	ASK_SELECT,
	ASK_YN,
	ASK_STRING,
	ASK_INTEGER,
	ASK_INTEGERNULL,
} askType;

void ask(askType X_type, const char *S_question, wsUint N_numSel, ...)
{
	va_list	H_varList;

	memset(&ANS, 0x00, sizeof(ANS));
__ask:
	va_start(H_varList, N_numSel);
	printf("-------------------------------------------------------------------------------\n"
		"*** %s\n"
		"-------------------------------------------------------------------------------\n",
		S_question);
	for (wsUint i=1 ; i<=N_numSel ; i++) {
		printf(" %d : %s\n", i, va_arg(H_varList, char *));
	}
	switch (X_type) {
	case ASK_YN:
		printf("(y/n) ");
		break;
	default:
		break;
	}
	printf("> ");
	if (!fgets(ANS.S_sel, 256, stdin))
		exit(1);

	switch (X_type) {
	case ASK_SELECT:
	case ASK_INTEGER:
	case ASK_INTEGERNULL: {
			/* Trim */
			for (char *Sp=ANS.S_sel+strlen(ANS.S_sel)-1 ; Sp>=ANS.S_sel &&
				(*Sp == '\t' || *Sp == '\n' || *Sp == '\r') ; Sp--)
				*Sp = '\0';
			/* Allow null value */
			if (ANS.S_sel[0] == '\0' && X_type == ASK_INTEGERNULL) {
				ANS.N_sel = 0xffffffff;
				break;
			}
			/* Check value */
			char*	Sp		= NULL;
			int		N_val	= (int)strtol(ANS.S_sel, &Sp, 10);
			if (!Sp || Sp[0] || !ANS.S_sel[0]) {
				printf("Error : [%s] is invalid integer value\n", ANS.S_sel);
				goto __ask;
			}
			ANS.N_sel = N_val;
			if (X_type == ASK_INTEGER || X_type == ASK_INTEGERNULL)
				break;
			if (ANS.N_sel < 1 || ANS.N_sel > N_numSel) {
				printf("Error : Select one of selection\n");
				goto __ask;
			}
		} break;
	case ASK_YN:
		/* Trim */
		for (char *Sp=ANS.S_sel+strlen(ANS.S_sel)-1 ; Sp>=ANS.S_sel &&
			(*Sp == '\t' || *Sp == '\n' || *Sp == '\r') ; Sp--)
			*Sp = '\0';
		if ((ANS.S_sel[0] == 'y' && !ANS.S_sel[1]) ||
			(ANS.S_sel[0] == 'Y' && !ANS.S_sel[1]) ||
			(ANS.S_sel[0] == 'n' && !ANS.S_sel[1]) ||
			(ANS.S_sel[0] == 'N' && !ANS.S_sel[1])) {
			ANS.B_yes = ANS.S_sel[0] == 'y' || ANS.S_sel[0] == 'Y';
		} else {
			printf("Error : The answer should be 'y/Y' or 'n/N'\n");
			goto __ask;
		}
		break;
	case ASK_STRING:
		/* Trim */
		for (char *Sp=ANS.S_sel+strlen(ANS.S_sel)-1 ; Sp>=ANS.S_sel &&
			(*Sp == '\t' || *Sp == '\n' || *Sp == '\r') ; Sp--)
			*Sp = '\0';
		break;
	default:
		break;
	}


	va_end(H_varList);
}

void ask(askType X_type, const char *S_question)
{
	ask(X_type, S_question, 0);
}

void _ask_gstest();
void _ask_pca();
void _ask_scoretest();
void _ask_gxgtest();
void _ask_recode();
void _ask_cor();
void _ask_filter();

wsStrCst getArgType(char u)
{
	switch (u) {
	case 'p': return "Positive integer";
	case 's': return "String with no whitespace";
	default: return "Unknown";
	}
}

static char S_itrtError[1024];
static int N_itrtEL = 0;
char	B_cmdExists	= 0;
void itrtSetError(int N_errLevel, wsStrCst S_errFormat, ...)
{
	if (N_itrtEL > N_errLevel) return;
	va_list	H_varList;
	va_start(H_varList, S_errFormat);
	int X = sprintf(S_itrtError, "[ERROR] ");
	vsprintf(S_itrtError+X, S_errFormat, H_varList);
	strcat(S_itrtError, "\n");
	pverbose("Error [%d] level [%s] set\n", N_errLevel, S_itrtError);
	//LOGnf(S_buf);
	va_end(H_varList);
}

void itrtClearError()
{
	N_itrtEL = 0;
	S_itrtError[0] = '\0';
	pverbose("Error removed\n");
}

void itrtSysErr(wsStrCst S_errFormat, ...)
{
	char S_buf[1024];
	va_list	H_varList;
	va_start(H_varList, S_errFormat);
	int X = sprintf(S_buf, "[SYSERR] ");
	vsprintf(S_buf+X, S_errFormat, H_varList);
	strcat(S_buf, "\n");
	LOGnf(S_buf);
	va_end(H_varList);
}

void itrtMsg(wsStrCst S_errFormat, ...)
{
	char S_buf[1024];
	va_list	H_varList;
	va_start(H_varList, S_errFormat);
	int X = sprintf(S_buf, "[MESSAGE] ");
	vsprintf(S_buf+X, S_errFormat, H_varList);
	strcat(S_buf, "\n");
	LOGnf(S_buf);
	va_end(H_varList);
}

bool extractArgs(char *S_params, xItrtCmd &X_cmd, vParam &X_params)
{
	char*	Sp_args = X_cmd.S_fmt;
	char*	Sp_par	= S_params;
	bool	B_ret	= true;
	int		i		= 1;			// Index of argument being processed (1-based), for logging purpose

	/* For all arguments */
	char	*Sp_end	= NULL;
	for (char *p=Sp_args ; *p ; p++,i++,Sp_par=Sp_end) {
		xParam	X_par	= { 0, };
		char	*Sp_chk	= NULL;
		if (Sp_par == NULL) {
			itrtSetError(1, "Command [%s] requires [%s] type for [%d]th argument, "
				"but not assigned!", X_cmd.S_cmd, getArgType(*p), i);
			return false;
		}
		getString(&Sp_par, &Sp_end);

		switch (*p) {
		case 'p':
			X_par.N_intVal = strtol(Sp_par, &Sp_chk, 10);
			if (Sp_chk[0] || X_par.N_intVal <= 0) {
				itrtSetError(2, "Command [%s] requires [%s] type for [%d]th argument, "
					"but [%s] is inappropriate!", X_cmd.S_cmd, getArgType(*p),
					i, Sp_par);
//				B_ret = false;
				return false;
			}
			break;
		case 's':
			X_par.S_strVal = strdup(Sp_par);
			break;
		}
		X_params.push_back(X_par);
	}

	/* #param check */
	if (X_params.size() != strlen(X_cmd.S_fmt)) {
		itrtSetError(1, "Command [%s] requires [%d] arguments, but [%d] assigned!",
			X_cmd.S_cmd, strlen(X_cmd.S_fmt), X_params.size());
		return false;
	}

	return B_ret;
}

void itrtHelp(xItrtCmd *Xa_cmds, wsUint N_cmds, char *S_arg, wsUint N_cmdPrt)
{
	char S_fmt[32];

	if (S_arg == NULL) {
		sprintf(S_fmt, "\t%%%ds - %%s\n", N_cmdPrt);
		/* List all available commands */
		for (wsUint i=0 ; i<N_cmds ; i++)
			printf(S_fmt, Xa_cmds[i].S_cmd, Xa_cmds[i].S_desc);
	} else {
		sprintf(S_fmt, " %%%ds - %%s\n", N_cmdPrt);
		/* Have specific query */
		char *Sp_shouldNull = NULL;
		getString(&S_arg, &Sp_shouldNull);
		if (Sp_shouldNull && Sp_shouldNull[0])
			itrtMsg("Unused part[%s] will be ignored");

		/* Find command */
		wsUint N_hit = 0;
		for (wsUint i=0 ; i<N_cmds ; i++)
			if (!stricmp(Xa_cmds[i].S_cmd, S_arg)) {
				char *Sp_desc = Xa_cmds[i].S_optDesc;
				char *q = strtok(Sp_desc, ",");
				printf(S_fmt, Xa_cmds[i].S_cmd, Xa_cmds[i].S_optDesc);
				printf("-----------------------------------------------\n");
				for (char *p=Xa_cmds[i].S_fmt ; *p ; p++) {
					if (q)
						printf("[%s] %s\n", getArgType(*p), q);
					else
						printf("[%s] <NO DESCRIPTION AVAILABLE>\n", getArgType(*p));
					q = strtok(NULL, ",");
				}
				printf("-----------------------------------------------\n");
				N_hit++;
			}
		/* Error - no such command */
		if (!N_hit)
			itrtSetError(1, "Command [%s] does not exist", S_arg);
	}
}

void interactiveShell(xItrtCmd *Xa_cmds, wsUint N_cmds, void *Vp_data)
{
	char *S_cline = new char[512];

	/* Longest command length */
	wsUint N_cmdPrt = 0;
	for (wsUint i=0 ; i<N_cmds ; i++) {
		wsUint N_len = (wsUint)strlen(Xa_cmds[i].S_cmd);
		if (N_len > N_cmdPrt) N_cmdPrt = N_len;
	}

	LOGnf("\n"
		" Welcome to the interactive shell of " TOOLSET_NAME ".\n"
		" Type 'help' <enter> to see the help, and\n"
		" type 'help <command>' <enter> to see the detailed help.\n"
		" Enter 'exit' to terminate the interactive mode.\n"
		"\n");

	fprintf(stderr, " > ");
	while (fgets(S_cline, 512, stdin)) {
		/* Clear error */
		S_itrtError[0] = '\0';

		/* Chop */
		char	*p = S_cline + strlen(S_cline) - 1;
		while (*p == '\r' || *p == '\n' || *p == ' ' || *p == '\t') *(p--) = '\0';
		if (!stricmp("exit", S_cline)) break;
		B_cmdExists	= 0;
		char	*S_arg1		= NULL;
		wsUint	i			= 0;

		if (!S_cline[0]) goto _next;

		/* Get the first command part */
		getString(&S_cline, &S_arg1);

		/* Special for 'help' */
		if (!stricmp("help", S_cline)) {
			itrtHelp(Xa_cmds, N_cmds, S_arg1, N_cmdPrt);
			goto _next;
		}

		/* Find command */
		S_itrtError[0] = '\0';
		for (i = 0; i<N_cmds; i++) {
			vParam	V_param;
			char*	Sp_prevErr = strdup(S_itrtError);
			if (!stricmp(Xa_cmds[i].S_cmd, S_cline)) {
				B_cmdExists = 1;
				itrtClearError();
				if (extractArgs(S_arg1, Xa_cmds[i], V_param)) {
					if (Xa_cmds[i].H_func(Xa_cmds[i], V_param, Vp_data) == 0) {
						/* Clear previous errors because it has done */
						if (strcmp(Sp_prevErr, S_itrtError) == 0)
							S_itrtError[0] = '\0';
						break;
					}
				}
			} else
				strcpy(S_itrtError, Sp_prevErr);
			free(Sp_prevErr);
		}
		if (!B_cmdExists)
			itrtSetError(1, "Command [%s] does not exist", S_cline);

_next:
		/* Print if error */
		if (S_itrtError[0])
			LOGnf(S_itrtError);

		fprintf(stderr, " > ");
	}
	delete [] S_cline;
}

void __queryVariantValueAdditive(xItrtData_EXPL *Xp_d, wsUint N_idx,
	wsUint N_nonZonly=0)
{
	cIO*		Cp_IO		= Xp_d->Cp_IO;

	vSampPtr&	Xa_samp	= Cp_IO->getSample();
	char**		Na_data	= Cp_IO->getGenotype();
	wsUint		N_samp	= Cp_IO->sizeSample();
	char		S_fmtPrint[32];
	char		S_fmtHeader[32];
	sprintf(S_fmtPrint, "%%%ds	%%d\n", (int)Xp_d->N_sampName);
	sprintf(S_fmtHeader, "%%%ds	%%s\n", (int)Xp_d->N_sampName);

	printf(S_fmtHeader, "IID", "ADDITIVE");
	wsUint		Na_cnts[4] = { 0, };
	for (wsUint i=0 ; i<N_samp ; i++) {
		char N_geno = Na_data[i][N_idx];
		if (!N_nonZonly || (N_nonZonly && N_geno)) {
			if (isMissingReal(N_geno)) {
				printf(S_fmtHeader, Xa_samp[i]->S_IID.c_str(), "<NA>");
				Na_cnts[3]++;
			} else {
				printf(S_fmtPrint, Xa_samp[i]->S_IID.c_str(), N_geno);
				Na_cnts[(int)N_geno]++;
			}
		}
	}

	printf("\n\n## Genotype summary\nNUM_0	%d\n", Na_cnts[0]);
	printf("NUM_1	%d\n", Na_cnts[1]);
	printf("NUM_2	%d\n", Na_cnts[2]);
	printf("NUM_NA	%d\n", Na_cnts[3]);
}

void __querySampleValueAdditive(xItrtData_EXPL *Xp_d, wsUint N_idx,
								wsUint N_nonZonly=0)
{
	cIO*		Cp_IO	= Xp_d->Cp_IO;

	vVariant&	Xa_vrt	= Cp_IO->getVariant();
	char**		Na_data	= Cp_IO->getGenotype();
	wsUint		N_vrt	= Cp_IO->sizeVariant();
	char		S_fmtPrint[32];
	char		S_fmtHeader[32];
	sprintf(S_fmtPrint, "%%%ds	%%d\n", (int)Xp_d->N_vrtName);
	sprintf(S_fmtHeader, "%%%ds	%%s\n", (int)Xp_d->N_vrtName);

	printf(S_fmtHeader, "VARIANT", Cp_IO->getSample()[N_idx]->S_IID.c_str());
	wsUint		Na_cnts[4] = { 0, };
	for (wsUint i=0 ; i<N_vrt ; i++) {
		char N_geno = Na_data[N_idx][i];
		if (!N_nonZonly || (N_nonZonly && N_geno)) {
			if (isMissingReal(N_geno)) {
				printf(S_fmtHeader, Xa_vrt[i].name, "<NA>");
				Na_cnts[3]++;
			} else {
				printf(S_fmtPrint, Xa_vrt[i].name, N_geno);
				Na_cnts[(int)N_geno]++;
			}
		}
	}

	printf("\n\n## Genotype summary\nNUM_0	%d\n", Na_cnts[0]);
	printf("NUM_1	%d\n", Na_cnts[1]);
	printf("NUM_2	%d\n", Na_cnts[2]);
	printf("NUM_NA	%d\n", Na_cnts[3]);
}

ITRTFUNC(func_dmn)
{
	xItrtData_EXPL
				*Xp_d		= (xItrtData_EXPL *)Vp_data;
	cIO*		Cp_IO		= Xp_d->Cp_IO;
	if (!Cp_IO) {
		itrtSysErr("SYSERR : Dataset is NULL");
		return 1;
	}
	/* Fetch first element */
	wsStrCst		Sp_vrtName	= V_param[0].S_strVal;
	mDataIdx_it	X_res		= Xp_d->Xm_vrt2idx.find(Sp_vrtName);
	/* Invalid variant name */
	if (X_res == Xp_d->Xm_vrt2idx.end()) {
		itrtSetError(2, "Variant [%s] does not exist in the final dataset",
			Sp_vrtName);
		return 1;
	}
	/* Do print */
	__queryVariantValueAdditive(Xp_d, X_res->second,
		!stricmp("dmna", X_cmd.S_cmd));

	/* Success */ return 0;
}

ITRTFUNC(func_dmi)
{
	xItrtData_EXPL
				*Xp_d		= (xItrtData_EXPL *)Vp_data;
	cIO*		Cp_IO		= Xp_d->Cp_IO;
	if (!Cp_IO) {
		itrtSysErr("SYSERR : Dataset is NULL");
		return 1;
	}
	/* Fetch first element */
	wsUint		N_vrt		= Cp_IO->sizeVariant();
	int			N_idx		= V_param[0].N_intVal;
	/* Invalid variant index */
	if (N_idx<=0 || (wsUint)N_idx>N_vrt) {
		itrtSetError(2, "Dataset have only [%d] variants, so requested index[%d] "
			"is invalid!", N_vrt, N_idx);
		return 1;
	}
	/* Do print */
	__queryVariantValueAdditive(Xp_d, N_idx-1,
		!stricmp("dmia", X_cmd.S_cmd));

	/* Success */ return 0;
}

ITRTFUNC(func_dsn)
{
	xItrtData_EXPL
				*Xp_d	= (xItrtData_EXPL *)Vp_data;
	cIO*		Cp_IO	= Xp_d->Cp_IO;
	if (!Cp_IO) {
		itrtSysErr("SYSERR : Dataset is NULL");
		return 1;
	}
	/* Fetch first element */
	wsStrCst		Sp_IID	= V_param[0].S_strVal;
	mSamp_it	X_res	= Cp_IO->getSampleData().find(Sp_IID);
	/* Invalid variant name */
	if (X_res == Cp_IO->getSampleData().end()) {
		itrtSetError(2, "Sample [%s] does not exist in the final dataset", Sp_IID);
		return 1;
	}
	/* Do print */
	__querySampleValueAdditive(Xp_d, X_res->second.N_idx,
		!stricmp("dsna", X_cmd.S_cmd));

	/* Success */ return 0;
}

ITRTFUNC(func_dsi)
{
	xItrtData_EXPL
				*Xp_d	= (xItrtData_EXPL *)Vp_data;
	cIO*		Cp_IO	= Xp_d->Cp_IO;
	if (!Cp_IO) {
		itrtSysErr("SYSERR : Dataset is NULL");
		return 1;
	}
	/* Fetch first element */
	wsUint		N_samp	= Cp_IO->sizeSample();
	int			N_idx	= V_param[0].N_intVal;
	/* Invalid variant index */
	if (N_idx<=0 || (wsUint)N_idx>N_samp) {
		itrtSetError(2, "Dataset have only [%d] samples, so requested index[%d] "
			"is invalid!", N_samp, N_idx);
		return 1;
	}
	/* Do print */
	__querySampleValueAdditive(Xp_d, N_idx-1, !stricmp("dsia", X_cmd.S_cmd));

	/* Success */ return 0;
}

ITRTFUNC(func_lm)
{
	xItrtData_EXPL
				*Xp_d		= (xItrtData_EXPL *)Vp_data;
	cIO*		Cp_IO		= Xp_d->Cp_IO;
	vVariant&	Xv_vrt		= Cp_IO->getVariant();
	if (!Cp_IO) {
		itrtSysErr("SYSERR : Dataset is NULL");
		return 1;
	}
	/* Fetch first element */
	wsUint		N_vrt		= Cp_IO->sizeVariant();
	int			N_from		= V_param[0].N_intVal;
	int			N_to		= N_from + V_param[1].N_intVal - 1;
	/* Invalid variant index */
	if (N_from<=0 || (wsUint)N_from>N_vrt) {
		itrtSetError(2, "Dataset have only [%d] variants, so requested staring index[%d] "
			"is invalid!", N_vrt, N_from);
		return 1;
	}
	/* Adjust size */
	if ((wsUint)N_to>N_vrt) N_to = N_vrt;
	/* Do print */
	for (int i=N_from ; i<=N_to ; i++)
		printf("%08d	%s\n", i, Xv_vrt[i-1].name);

	/* Success */ return 0;
}

ITRTFUNC(func_lm2)
{
	xItrtData_EXPL
		*Xp_d = (xItrtData_EXPL *)Vp_data;
	cIO*		Cp_IO = Xp_d->Cp_IO;
	vVariant&	Xv_vrt = Cp_IO->getVariant();
	if (!Cp_IO) {
		itrtSysErr("SYSERR : Dataset is NULL");
		return 1;
	}

	const TRexChar*	Sp_err = NULL;
	TRex*		X_regex = trex_compile(V_param[0].S_strVal, &Sp_err);
	if (!X_regex) {
		itrtSetError(2, "Parameter [%s] is error by [%s]!", V_param[0].S_strVal,
			Sp_err ? Sp_err : "???");
		return 1;
	}

	wsUint N_match = 0;
	const TRexChar*	s = NULL;
	const TRexChar*	e = NULL;
	wsUint I = 0;
	FOREACHDO(vVariant_it, Xv_vrt, i, I++) {
		if (!trex_search(X_regex, i->name, &s, &e)) continue;
		printf("%08d	%s\n", I, i->name);
		N_match++;
	}
	if (N_match == 0) {
		itrtSetError(2, "No match with parameter [%s]", V_param[0].S_strVal);
		return 1;
	}

	/* Success */ return 0;
}

ITRTFUNC(func_ls)
{
	xItrtData_EXPL
				*Xp_d		= (xItrtData_EXPL *)Vp_data;
	cIO*		Cp_IO		= Xp_d->Cp_IO;
	vSampPtr	&Xa_samp	= Cp_IO->getSample();
	if (!Cp_IO) {
		itrtSysErr("SYSERR : Dataset is NULL");
		return 1;
	}
	/* Fetch first element */
	wsUint		N_samp		= Cp_IO->sizeSample();
	int			N_from		= V_param[0].N_intVal;
	int			N_to		= N_from + V_param[1].N_intVal - 1;
	/* Invalid variant index */
	if (N_from<=0 || (wsUint)N_from>N_samp) {
		itrtSetError(2, "Dataset have only [%d] samples, so requested staring index[%d] "
			"is invalid!", N_samp, N_from);
		return 1;
	}
	/* Adjust size */
	if ((wsUint)N_to>N_samp) N_to = N_samp;
	/* Do print */
	for (int i=N_from ; i<=N_to ; i++)
		printf("%08d	%s\n", i, Xa_samp[i-1]->S_IID.c_str());

	/* Success */ return 0;
}

ITRTFUNC(func_ls2)
{
	xItrtData_EXPL
		*Xp_d		= (xItrtData_EXPL *)Vp_data;
	cIO*		Cp_IO		= Xp_d->Cp_IO;
	vSampPtr	&Xa_samp	= Cp_IO->getSample();
	if (!Cp_IO) {
		itrtSysErr("SYSERR : Dataset is NULL");
		return 1;
	}
	/* Fetch first element */
	const TRexChar*	Sp_err		= NULL;
	TRex*		X_regex		= trex_compile(V_param[0].S_strVal, &Sp_err);
	if (!X_regex) {
		itrtSetError(2, "Parameter [%s] is error by [%s]!", V_param[0].S_strVal,
			Sp_err?Sp_err:"???");
		return 1;
	}

	wsUint N_match = 0;
	const TRexChar*	s	= NULL;
	const TRexChar*	e	= NULL;
	FOREACH (vSampPtr_it, Xa_samp, i) {
		if (!trex_search(X_regex, (*i)->S_IID.c_str(), &s, &e)) continue;
		printf("%08d	%s\n", (*i)->N_idx, (*i)->S_IID.c_str());
		N_match++;
	}
	if (N_match == 0)
		itrtSetError(2, "No match with parameter [%s]", V_param[0].S_strVal);

	/* Success */ return 0;
}

ITRTFUNC(func_lf)
{
	xItrtData_EXPL
		*Xp_d		= (xItrtData_EXPL *)Vp_data;
	cIO*		Cp_IO		= Xp_d->Cp_IO;
	vSampPtr	&Xa_samp	= Cp_IO->getSample();
	if (!Cp_IO) {
		itrtSysErr("SYSERR : Dataset is NULL");
		return 1;
	}
	/* Fetch first element */
	const TRexChar*	Sp_err		= NULL;
	TRex*		X_regex		= trex_compile(V_param[0].S_strVal, &Sp_err);
	if (!X_regex) {
		itrtSetError(2, "Parameter [%s] is error by [%s]!", V_param[0].S_strVal,
			Sp_err?Sp_err:"???");
		return 1;
	}

	wsUint N_match = 0;
	const TRexChar*	s	= NULL;
	const TRexChar*	e	= NULL;
	FOREACH (vSampPtr_it, Xa_samp, i) {
		if (!trex_search(X_regex, (*i)->S_FID.c_str(), &s, &e)) continue;
		printf("%08d	%s\n", (*i)->N_idx, (*i)->S_IID.c_str());
		N_match++;
	}
	if (N_match == 0)
		itrtSetError(2, "No match with parameter [%s]", V_param[0].S_strVal);

	/* Success */ return 0;
}

ITRTFUNC(func_sz)
{
	xItrtData_EXPL
		*Xp_d		= (xItrtData_EXPL *)Vp_data;
	cIO*	Cp_IO	= Xp_d->Cp_IO;
	wsUint	N_samp	= Cp_IO->sizeSample();
	wsUint	N_fam	= (wsUint)Cp_IO->getFamilyData().size();
	wsUint	N_vrt	= Cp_IO->sizeVariant();
	wsUint	N_fnd	= Cp_IO->sizeFounder();
	printf("N_FAMILY	%d\n"
		"N_FOUNDER	%d\n"
		"N_SAMPLE	%d\n"
		"N_VARIANT	%d\n", N_fam, N_fnd, N_samp, N_vrt);

	/* Success */ return 0;
}

ITRTFUNC(func_si)
{
	xItrtData_EXPL
		*Xp_d = (xItrtData_EXPL *)Vp_data;
	cIO*		Cp_IO = Xp_d->Cp_IO;
	vSampPtr	&Xa_samp = Cp_IO->getSample();
	if (!Cp_IO) {
		itrtSysErr("SYSERR : Dataset is NULL");
		return 1;
	}
	/* Fetch first element */
	const TRexChar*	Sp_err = NULL;
	TRex*		X_regex = trex_compile(V_param[0].S_strVal, &Sp_err);
	if (!X_regex) {
		itrtSetError(2, "Parameter [%s] is error by [%s]!", V_param[0].S_strVal,
			Sp_err ? Sp_err : "???");
		return 1;
	}

	/* If there is pheno & covariates */
	wsMat	Ra_phe = Cp_IO->getPhenos();
	wsUint	N_phe = Cp_IO->sizePheno();
	wsMat	Ra_cov = Cp_IO->getCovariates();
	wsUint	N_cov = Cp_IO->sizeCovar();
	vPheno&	Xv_phe = Cp_IO->getPhenoInfo();
	vCovar&	Xv_cov = Cp_IO->getCovInfo();

	wsUint N_match = 0;
	const TRexChar*	s = NULL;
	const TRexChar*	e = NULL;
	FOREACH(vSampPtr_it, Xa_samp, i) {
		if (!trex_search(X_regex, (*i)->S_IID.c_str(), &s, &e)) continue;
		xSample& X_smp = *(*i);
		printf("FAMILY_ID     %s", X_smp.S_FID.c_str());
		printf("INDIVIDUAL_ID %s", X_smp.S_IID.c_str());
		printf("SEX           %s", X_smp.N_sex == 2 ? "FEMALE" : (X_smp.N_sex == 1 ? "MALE" : "UNKNOWN"));
		LOOP(j, N_phe)
			printf("%13s %g", Xv_phe[j].S_name.c_str(), Ra_phe[j][X_smp.N_idx]);
		LOOP(j, N_cov)
			printf("%13s %g", Xv_cov[j].Sp_varName, Ra_cov[j][X_smp.N_idx]);
		N_match++;
	}
	if (N_match == 0)
		itrtSetError(2, "No match with parameter [%s]", V_param[0].S_strVal);

	/* Success */ return 0;
}

void interactiveExplore(cIO *Cp_IO)
{
	mDataIdx	Xm_vrt2idx;
	size_t		N_vrtName	= 0;
	size_t		N_sampName	= 0;

	/* Build variant idx and print length of variant */
	vVariant&	Xv_vrt		= Cp_IO->getVariant();
	wsUint		I			= 0;
	FOREACHDO (vVariant_it, Xv_vrt, i, I++) {
		string s = i->name;
		Xm_vrt2idx[s] = I; /* Variant name to data index */
		if (N_vrtName < s.length())
			N_vrtName = s.length();
	}
	/* Build length of sample */
	vSampPtr&	Xv_samp		= Cp_IO->getSample();
	FOREACH (vSampPtr_it, Xv_samp, i) {
		if (N_sampName < (*i)->S_IID.length())
			N_sampName = (*i)->S_IID.length();
	}
	/* Build command mapping */
	xItrtCmd X_cmds[] = {
		{ "dmn",  "s",  CMD_EXPL_QUERYDATA_BY_SNPNAME,       func_dmn, "Variant name",
			"Show genotype by variant name" },
		{ "dmi",  "p",  CMD_EXPL_QUERYDATA_BY_SNPIDX,        func_dmi, "Index(1-base)",
			"Show genotype by stored index of variant" },
		{ "dmna", "s",  CMD_EXPL_QUERYDATA_BY_SNPNAME_AVAIL, func_dmn, "Variant name",
			"Show genotype by variant name, non-majorhomo only" },
		{ "dmia", "p",  CMD_EXPL_QUERYDATA_BY_SNPIDX_AVAIL,  func_dmi, "Index(1-base)",
			"Show genotype by stored index of variant, non-majorhomo only" },
		{ "dsn",  "s",  CMD_EXPL_QUERYDATA_BY_SAMPIID,       func_dsn, "Sample IID",
			"Show genotype by IID" },
		{ "dsi",  "p",  CMD_EXPL_QUERYDATA_BY_SAMPIDX,       func_dsi, "Index(1-base)",
			"Show genotype by stored index of sample" },
		{ "dsna", "s",  CMD_EXPL_QUERYDATA_BY_SAMPIID_AVAIL, func_dsn, "Sample IID",
			"Show genotype by IID, non-majorhomo only" },
		{ "dsia", "p",  CMD_EXPL_QUERYDATA_BY_SAMPIDX_AVAIL, func_dsi, "Index(1-base)",
			"Show genotype by stored index of sample, non-majorhomo only" },
		{ "lm",   "pp", CMD_EXPL_QUERY_SNPNAME,              func_lm,  "Index(1-base),count",
			"Show variant names by index" },
		{ "lm",   "pp", CMD_EXPL_QUERY_SNPNAME2,             func_lm2, "expression",
			"Search variant names" },
		{ "ls",   "pp", CMD_EXPL_QUERY_SAMPIID,              func_ls,  "Index(1-base),count",
			"Show sample IIDs by index" },
		{ "ls",   "s",  CMD_EXPL_QUERY_SAMPIID2,             func_ls2, "expression",
			"Search sample IIDs" },
		{ "lf",   "s",  CMD_EXPL_QUERY_SAMPFID,              func_lf,  "expression",
			"Search samples by FID" },
		{ "sz",   "",   CMD_EXPL_QUERY_SIZE,                 func_sz,  "<none>",
			"View data size" },
		{ "si",   "",   CMD_EXPL_QUERY_SAMPINFO,             func_si,  "Sample IID",
			"See sample information by IID" },
	};
	/* Build required data */
	xItrtData_EXPL X_data = { Cp_IO, Xm_vrt2idx, N_sampName, N_vrtName };
	interactiveShell(X_cmds, len(X_cmds, xItrtCmd), &X_data);
}

char isPathValid(wsStrCst S_path, wsStrCst S_posExt=NULL)
{
	/* Check the file is existing */
	FILE* H_fp = fopen(S_path, "r");
	if (H_fp == NULL) {
		/* If this not supports possible extension, just return 0 */
		if (S_posExt == NULL) goto _err;

		char S_extPath[MAX_PATH];
		sprintf(S_extPath, "%s.%s", S_path, S_posExt);
		/* Try with possible extension */
		H_fp = fopen(S_extPath, "r");
		if (H_fp == NULL) goto _err;
	}
	return 1;
_err:
	notice("Given path [%s] is not valid, re-enter path\n", S_path);
	return 0;
}

void interactiveExecution()
{
	/* Determining input */
	ask(ASK_SELECT, "Enter the type of data you will use", 10,
		"Trio simulation", "Family simulation", "BED file format",
		"PED/RAW file format", "LGEN file format", "Transposed PED format",
		"VCF format file", "GEN file format", "Binary GEN file", "Dosage format");
	switch (ANS.N_sel) {
	case 1:
		sprintf(ANS.S_sel, "%d", ANS.N_sel);
		OPTION().assign("simtrio", ANS.S_sel);
		break;
	case 2:
		OPTION().assign("simfam", "1");
		ask(ASK_SELECT, "Enter the type of family you will use", 4,
			"Independent samples",
			"Trio family",
			"Extended family",
			"Arbitrary type of family (via FAM file)");
		switch (ANS.N_sel) {
		case 1: OPTION().assign("indep", "1");
			ask(ASK_STRING, "Enter the number of samples being simulated");
			OPTION().assign("szfam", ANS.S_sel);
			break;
		case 2: OPTION().assign("trio", "1");
			ask(ASK_STRING, "Enter the number of families being simulated");
			OPTION().assign("szfam", ANS.S_sel);
			break;
		case 3: OPTION().assign("extfam", "1");
			ask(ASK_STRING, "Enter the number of families being simulated");
			OPTION().assign("szfam", ANS.S_sel);
			break;
		case 4:
			ask(ASK_STRING, "Enter the path of FAM file located");
			OPTION().assign("fam", ANS.S_sel);
			break;
		}
		ask(ASK_STRING, "Enter the number of variants being simulated");
		OPTION().assign("szvar", ANS.S_sel);
		break;
	case 3:
		do {
			ask(ASK_STRING, "Enter the path of BED file located");
		} while (!isPathValid(ANS.S_sel, "bed"));
		OPTION().assign("bed", ANS.S_sel);
		break;
	case 4:
		do {
			ask(ASK_STRING, "Enter the path of PED/RAW file located (with extension)");
		} while (!isPathValid(ANS.S_sel));
		OPTION().assign("ped", ANS.S_sel);
		break;
	case 5:
		do {
			ask(ASK_STRING, "Enter the path of LGEN file located (with extension)");
		} while (!isPathValid(ANS.S_sel));
		OPTION().assign("lgen", ANS.S_sel);
		break;
	case 6:
		do {
			ask(ASK_STRING, "Enter the path of transposed PED file located (with extension)");
		} while (!isPathValid(ANS.S_sel));
		OPTION().assign("tped", ANS.S_sel);
		break;
	case 7:
		do {
			ask(ASK_STRING, "Enter the path of VCF file located (with extension)");
		} while (!isPathValid(ANS.S_sel));
		OPTION().assign("vcf", ANS.S_sel);
		ask(ASK_YN, "Is this VCF file have corresponding FAM file?");
		if (ANS.B_yes) {
			ask(ASK_STRING, "Enter the path of FAM file located");
			OPTION().assign("fam", ANS.S_sel);
		}
		break;
	case 8:
		do {
			ask(ASK_STRING, "Enter the path of GEN file located (with extension)");
		} while (!isPathValid(ANS.S_sel));
		OPTION().assign("gen", ANS.S_sel);
		break;
	case 9:
		do {
			ask(ASK_STRING, "Enter the path of binary GEN file located (with extension)");
		} while (!isPathValid(ANS.S_sel));
		OPTION().assign("bgen", ANS.S_sel);
		break;
	case 10:
		do {
			ask(ASK_STRING, "Enter the path of dosage file located (with extension)");
		} while (!isPathValid(ANS.S_sel));
		OPTION().assign("in", ANS.S_sel);
		break;
	}

	/* Determining additional input */
	ask(ASK_YN, "Do you want to load covariate/alternative phenotype?");
	if (ANS.B_yes) {
		ask(ASK_STRING, "Enter the path of covariate/alternative phenotype file located");
		OPTION().assign("pheno", ANS.S_sel);

		ask(ASK_STRING, "Enter the name of columns to be included as covariates\n"
			"('*' mean all columns)");
		OPTION().assign("cname", ANS.S_sel);

		ask(ASK_STRING, "Enter the name of columns to be included as alternative phenotype\n"
			"('*' mean all columns, leave empty mean use default phenotype)");
		OPTION().assign("pname", ANS.S_sel);
	}

	/* Determining filters */
	ask(ASK_YN, "Do you want to apply some filters on samples/variants?");
	wsUint N_selIdx = 0xffffffff;
	if (ANS.B_yes) do {
		ask(ASK_SELECT, "Select the filter what you want to apply", 2,
			"Selecting/filtering subset of variants",
			"Selecting/filtering subset of samples");

#define NUM_OF_VARFILTERS 12
		switch (ANS.N_sel) {
		case 1:
			ask(ASK_SELECT, "What kind of selecting/filtering scheme do you want?", NUM_OF_VARFILTERS,
				"Selecting/filtering variants by the names of variant",
				"Selecting/filtering variants by MAF of variant",
				"Selecting/filtering variants by MAC of variant",
				"Selecting/filtering variants by the p-value of HWE test",
				"Selecting/filtering variants by genetic distance of variant",
				"Selecting/filtering variants by genotyping rate of variant",
				"Selecting indels only",
				"Selecting SNVs only",
				"Selecting autosome variants",
				"Selecting sex-chromosome variants",
				"Selecting variants of specific chromosome(s)",
				"END SELECT");
			N_selIdx = ANS.N_sel;
			/* If requires sel/fil */
			if (N_selIdx <= 6) {
				char B_yes = -1;
				ask(ASK_YN, "If you want apply this filter to select variants, select (y/Y)");
				B_yes = ANS.B_yes;
				switch (N_selIdx) {
				case 1:
					do {
						ask(ASK_STRING, "Give the location of file contains the list of variants to be filtered/selected");
					} while (isPathValid(ANS.S_sel));
					OPTION().assign(ANS.B_yes?"incvariant":"filvariant",
						ANS.S_sel);
					break;
				case 2:
					ask(ASK_STRING, "Give the range-formatted MAF to be filtered/selected");
					OPTION().assign(ANS.B_yes?"filmaf":"filmaf", ANS.S_sel);
					break;
				case 3:
					ask(ASK_STRING, "Give the range-formatted MAC to be filtered/selected");
					OPTION().assign(ANS.B_yes?"incmac":"filmac", ANS.S_sel);
					break;
				case 4:
					ask(ASK_STRING, "Give the range-formatted p-value from HWE to be filtered/selected");
					OPTION().assign(ANS.B_yes?"inchwe":"filhwe", ANS.S_sel);
					break;
				case 5:
					ask(ASK_STRING, "Give the range-formatted genetic distance to be filtered/selected");
					OPTION().assign(ANS.B_yes?"incgdist":"filgdist", ANS.S_sel);
					break;
				case 6:
					ask(ASK_STRING, "Give the range-formatted genotyping rate to be filtered/selected");
					OPTION().assign(ANS.B_yes?"incgvar":"filgvar", ANS.S_sel);
					break;
				}
			} else switch (N_selIdx) { /* It is only for fixed sel/fil */
			case 7: OPTION().assign("indelonly", "1"); break;
			case 8: OPTION().assign("snvonly", "1"); break;
			case 9: OPTION().assign("autoonly", "1"); break;
			case 10: OPTION().assign("sexonly", "1"); break;
			case 11:
				ask(ASK_STRING, "Give the list of chromosomes to be selected");
				OPTION().assign("chr", ANS.S_sel);
				break;
			}
			break;
		}
	} while (N_selIdx != NUM_OF_VARFILTERS);

	/* Determining analyses */
#define NUM_OF_ANALYSES 9
	wsUint N_tmp;
	do {
		ask(ASK_SELECT, "Select the analysis what you want to perform", NUM_OF_ANALYSES,
			"Association analysis",
			"Association analysis with population stratification",
			"Family-based association",
			"Gene-level association analysis",
			"Gene-level association analysis with population stratification",
			"Family-based gene-level association analysis",
			"Correlation-related setting",
			"Apply filter before analysis",
			"END SELECT");
		N_tmp = ANS.N_sel;

		switch (N_tmp) {
			/* FIXME */
		case 1: _ask_scoretest(); break;
		case 2: _ask_gstest(); break;
		case 3: _ask_gxgtest(); break;
		case 4: _ask_pca(); break;
		case 5: _ask_recode(); break;
		case 6: _ask_cor(); break;
		}
	} while (N_tmp != NUM_OF_ANALYSES);
}

void _ask_gstest()
{
	OPTION().assign("genetest", (char *)"1");
	ask(ASK_STRING, "Gene-set test requires prevalence parameter, please assign this\n"
		"(Range should be 0 < prevalence < 1)\n"
		"(One value to give same prevalence for male/female,\n"
		" or two values to give different prevalence [male and female]");
	OPTION().assign("prevalence", ANS.S_sel);

	ask(ASK_STRING, "Gene-set test requires gene-set definition file,\n"
		"Enter the path of gene-set definition file located");
	OPTION().assign("set", ANS.S_sel);

	ask(ASK_YN, "Do you want to export the summary of given gene-set definition file?");
	if (ANS.B_yes)
		OPTION().assign("genesummary", (char *)"1");

	ask(ASK_YN, "Do you want to perform family-based CMC test?");
	if (ANS.B_yes) {
		OPTION().assign("pedcmc", (char *)"1");
		ask(ASK_STRING, "Input MAF that categorizes variant into 'rare' in PEDCMC");
		OPTION().assign("raremaf", ANS.S_sel);
	} else {
		ask(ASK_YN, "Will you include samples having missing phenotype into gene-set analysis?");
		if (ANS.B_yes)
			OPTION().assign("imputepheno", (char *)"1");
	}

	ask(ASK_STRING, "Please input the size of gene-set to be included\n"
		"(To include gene-sets having smaller size than given value, `<SIZE`)\n"
		"(To include gene-sets having larger size than given value, `>SIZE`)\n"
		"(To include gene-sets having specific range, `SMALLEST-LARGEST`)\n"
		"(Leave it blank to use all available gene-set)");
	if (ANS.S_sel[0])
		OPTION().assign("gsrange", ANS.S_sel);
}

void _ask_gxgtest()
{
	ask(ASK_SELECT, "Select GxG test to perform", 3,
		"Multiple Dimensionality Reduction (MDR) analysis",
		"BOOST analysis",
		"PLINK fast-epistasis analysis");
	switch (ANS.N_sel) {
	case 1: OPTION().assign("--mdr", (char *)"1"); break;
	case 2: OPTION().assign("--boost", (char *)"1"); break;
	case 3: OPTION().assign("--quickepi", (char *)"1"); break;
	}
}

void _ask_pca()
{
	char S_askBuf[1024];

	OPTION().assign("pca", (char *)"1");
	ask(ASK_YN, "Do you want to extract all possible PCs? "
		"it requires long time");
	if (!ANS.B_yes) {
		sprintf(S_askBuf, "Give the number of PCs you will extract\n"
			"(Default value is %s, leave it blank then default will used)",
			OPTION().getDefVal("npc"));
		ask(ASK_INTEGERNULL, S_askBuf);

		/* Default value handling */
		if (ANS.N_sel == 0xffffffff)
			strcpy(ANS.S_sel, OPTION().getDefVal("npc"));

		OPTION().assign("npc", ANS.S_sel);
	}
}

void _ask_scoretest()
{
	OPTION().assign("scoretest", (char *)"1");

	ask(ASK_YN, "Do you want to use ML method instead of REML method?");
	if (ANS.B_yes)
		OPTION().assign("ml", (char *)"1");

	/* FIXME : Need to be added */
}

void _ask_recode()
{
	ask(ASK_SELECT, "Select the format of recode", 6,
		"PLINK PED format",
		"PLINK binary PED format",
		"PLINK RAW format",
		"VCF format",
		"PLINK long file format",
		"PLINK transposed PED format");
	switch (ANS.N_sel) {
	case 1: OPTION().assign("makeped", (char *)"1"); break;
	case 2: OPTION().assign("makebed", (char *)"1"); break;
	case 3: OPTION().assign("makeraw", (char *)"1"); break;
	case 4: OPTION().assign("makevcf", (char *)"1"); break;
	case 5: OPTION().assign("makelgen", (char *)"1"); break;
	case 6: OPTION().assign("maketped", (char *)"1"); break;
	}
}

void _ask_cor()
{
	ask(ASK_SELECT, "Select the task what you will apply to correlation calculation", 4,
		"Using theoretical correlation (fast calculation)",
		"Using IBS instead of correlation",
		"Assign alternative correlation matrix into analysis",
		"Using median in the calculation of empirical correlation");
	switch (ANS.N_sel) {
	case 1: OPTION().assign("kinship", (char *)"1"); break;
	case 2: OPTION().assign("ibs", (char *)"1"); break;
	case 3:
		ask(ASK_STRING, "Input the path of alternative correlation file");
		OPTION().assign("cor", ANS.S_sel);
		break;
	case 4: OPTION().assign("medcor", (char *)"1"); break;
	}
}

void _ask_filter()
{
	ask(ASK_SELECT, "Select the type of filter what you will apply", 5,
		"SELECT ONLY specified variants from dataset",
		"REMOVE specified variants from dataset",
		"SELECT ONLY specified samples from dataset",
		"REMOVE specified samples from dataset",
		"Return to previous menu");
	switch (ANS.N_sel) {
	case 1: case 2:
		ask(ASK_STRING, "Enter the path of variant filter file\n");
		OPTION().assign(ANS.N_sel==1?"selvariant":"remvariant", ANS.S_sel);
		break;
	case 3: case 4:
		ask(ASK_STRING, "Enter the path of sample filter file\n");
		OPTION().assign(ANS.N_sel==3?"selsamp":"remsamp", ANS.S_sel);
		break;
	}
}

#else

void interactiveExplore(cIO *Cp_IO) {}

#endif

} // End namespace ONETOOL
