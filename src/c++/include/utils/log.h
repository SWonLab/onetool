#pragma once
#ifndef __WISARD_LOG_H__
#define __WISARD_LOG_H__

#ifdef _WIN32
#	include <WinSock.h>
#endif
#include <stdio.h>

/* Print something and do halt */
#define halt(...)		__halt(__FUNC__, __VA_ARGS__)
#define halt_fmt(...)	__halt_format(__FUNC__, __VA_ARGS__)
#ifdef _WIN32
#	define LOG_COLOR_CYAN	0x0b
#	define LOG_COLOR_RED	0x0c
#	define LOG_COLOR_GREEN	0x0a
#else
#	define LOG_COLOR_CYAN	36
#	define LOG_COLOR_RED	31
#	define LOG_COLOR_GREEN	32
#endif

#ifdef WISARD_LOG_NO_FUNC
/* Print something */
#	define LOG(...)		__notice(1, "MESG", __VA_ARGS__)
/* Print output */
#	define LOGoutput(...)	__output("OUTPUT", LOG_COLOR_CYAN, 1, "FILE", __VA_ARGS__)
#	define LOGwarn(...)	__output("WARNING", LOG_COLOR_RED, 1, "WARN", __VA_ARGS__)
#	define LOGnote(...)	__output("NOTE", LOG_COLOR_GREEN, 1, "NOTE", __VA_ARGS__)
/* Print something with no function print*/
#	define LOGnf(...)		__notice(1, NULL, __VA_ARGS__)
/* Print something only on console */
#	define notice(...)		__notice(0, "MESG", __VA_ARGS__)
/* Print something only on console with no function print */
#	define noticenf(...)	__notice(0, NULL, __VA_ARGS__)
/* Print something only on --verbose */
#	define lverbose(...)	__verbose(1, "MESG", __VA_ARGS__)
/* Print something only on console */
#	define pverbose(...)	__verbose(0, "MESG", __VA_ARGS__)
#	define verbosenf(...)	__verbose(0, NULL, __VA_ARGS__)
#else
/* Print something */
#	define LOG(...)		__notice(1, __FUNC__, __VA_ARGS__)
/* Print output */
#	define LOGoutput(...)	__output("OUTPUT", LOG_COLOR_CYAN, 1, __FUNC__, __VA_ARGS__)
#	define LOGwarn(...)	__output("WARNING", LOG_COLOR_RED, 1, __FUNC__, __VA_ARGS__)
#	define LOGnote(...)	__output("NOTE", LOG_COLOR_GREEN, 1, __FUNC__, __VA_ARGS__)
/* Print something with no function print*/
#	define LOGnf(...)		__notice(1, NULL, __VA_ARGS__)
/* Print something only on console */
#	define notice(...)		__notice(0, __FUNC__, __VA_ARGS__)
/* Print something only on console with no function print */
#	define noticenf(...)	__notice(0, NULL, __VA_ARGS__)
/* Print something only on --verbose */
#	define lverbose(...)	__verbose(1, __FUNC__, __VA_ARGS__)
/* Print something only on console */
#	define pverbose(...)	__verbose(0, __FUNC__, __VA_ARGS__)
#	define verbosenf(...)	__verbose(0, NULL, __VA_ARGS__)
#endif

namespace ONETOOL {

typedef enum _xErrType {
	WISARD_SUCCESS,
	WISARD_DUPL_SNPINGENESET,
	WISARD_DUPL_OPTION,

	/**/WISARD_CANT_DO_W_SOMEOPT,
	/**/WISARD_CANT_DO_WO_SOMEOPT,
	/**/WISARD_CANT_OPT_W_OTHEROPT,
	/**/WISARD_CANT_OPT_W_SOMESTATE,

	WISARD_CANT_MULTILOCUS_WO_INDEL,
	/**/WISARD_CANT_ONEALLELE_MISSING,
	WISARD_CANT_NONAUTOCHR_W_AUTOONLY,
	WISARD_CANT_ONESIDE_PARENT,
	WISARD_CANT_MQLS_WO_BLUPORPREV,
	WISARD_CANT_MQLS_W_BOTHBLUPPREV,
	WISARD_CANT_BINPHENO_W_LONGIANA,
	WISARD_CANT_PEDCMC_W_MULTIPHENO,
	WISARD_CANT_PEDCMC_W_CONTPHENO,
	WISARD_CANT_MEDCOR_W_CORPEARSON,
	WISARD_CANT_FQLS_W_CONTPHENO,
	WISARD_CANT_FQLS_W_BLUP,
	WISARD_CANT_GETINDEL_WO_INDEL,
	WISARD_CANT_EXCL_OPT,

	WISARD_FAIL_OPEN_SAMPSEL,
	/**/WISARD_FAIL_OPEN_URLWITHDESC,
	/**/WISARD_FAIL_OPEN_URL,
	/**/WISARD_FAIL_OPEN_FILE,
	/**/WISARD_FAIL_OPEN_FILEWITHDESC,
	/**/WISARD_FAIL_OPEN_GZ,
	/**/WISARD_FAIL_OPEN_GZWITHDESC,
	WISARD_FAIL_OPEN_RUNAS,
	WISARD_FAIL_AUTOOPEN_TFAMEXTMISS,
	WISARD_FAIL_AUTOOPEN_TFAMNOTEXIST,
	/**/WISARD_FAIL_INIT_WINSOCK,
	/**/WISARD_FAIL_INIT_LOG,
	WISARD_FAIL_SOCKCREATE,
	WISARD_FAIL_SOCKCONN,
	WISARD_FAIL_OPEN_ANNOSNP,
	WISARD_FAIL_MEMALLOC_W_PROC,
	/**/WISARD_FAIL_READ_BEDSIG,
	WISARD_FAIL_ESTIMATE_PARAM,
	WISARD_FAIL_FIT_NULLMODEL,
	/**/WISARD_FAIL_MEMALLOC_W_SYSMSG,
	/**/WISARD_FAIL_MEMALLOC,
	WISARD_FAIL_OPEN_GSFILE_IDX,
	WISARD_FAIL_LAUNCH_R,

	WISARD_INVL_VARIANTPOS,
	WISARD_INVL_ALLELEFORMAT,
	WISARD_INVL_NVRTOFSAMPLE,
	WISARD_INVL_NSAMPLEOFSNP,
	WISARD_INVL_OPTNAME,
	WISARD_INVL_OPTVAL_SC,
	WISARD_INVL_OPTVAL_GENERAL_INT,
	WISARD_INVL_OPTVAL_GENERAL_STR,
	WISARD_INVL_LONGIVARDIM,
	WISARD_INVL_WSOCKVER,
	WISARD_INVL_URLFORMAT,
	WISARD_INVL_ANNOGENE_FORMAT,
	WISARD_INVL_SNPLIST_FORMAT,
	WISARD_INVL_RANGELIST_FORMAT,
	WISARD_INVL_SAMPLIST_FORMAT,
	WISARD_INVL_PEDCONTENTS,
	WISARD_INVL_CHROMFORMAT,
	WISARD_INVL_BEDVERSION,
	WISARD_INVL_BEDSIZE,
	WISARD_INVL_CORFORMAT,
	WISARD_INVL_ALTCORSIZE,
	WISARD_INVL_NPREV_NPHENO,
	WISARD_INVL_CORFILE_MORESZSAMP,
	WISARD_INVL_CORFILE_LESSSZSAMP,
	WISARD_INVL_GSFMT_NONPOS_GSRANGE,
	WISARD_INVL_GSFMT_REV_GSRANGE,
	WISARD_INVL_GSFMT_CHR,
	WISARD_INVL_GSFMT_DBLEND,
	WISARD_INVL_BIMFORMAT,
	WISARD_INVL_BEDSIGNATURE,
	WISARD_DUPL_GSID,
	WISARD_INVL_PLINK_INPUT,
	WISARD_INVL_SAMPNAME,
	WISARD_INVL_NVRTMAP_W_NVRTPED,
	WISARD_INVL_VCFMETA,
	WISARD_INVL_VCFMANDTRY_COLNAME,
	WISARD_INVL_PEDCHAR,
	WISARD_INVL_PEDMISCHAR,
	WISARD_INVL_RANGE_NUMBER,
	WISARD_INVL_RANGE_REAL,
	WISARD_INVL_RANGE_STRING,
	WISARD_INVL_FILE_GENERAL,
	WISARD_INVL_FILE_INCMPL,
	WISARD_INVL_FILE_INVALID,
	WISARD_INVL_FILE_INVALID_DESC,
	WISARD_INVL_EOF_LINE,
	WISARD_INVL_OPERATOR,

	WISARD_NULL_GENEDEF,
	WISARD_NULL_SAMPLE,
	WISARD_NULL_OPTVALUE,
	WISARD_NULL_ESSENTIALOPT,
	WISARD_NULL_ESSENTIALOPT_OR,
	WISARD_NULL_BLUP,
	WISARD_NULL_ANNOGENENAME,
	WISARD_NULL_VTCOVAR,
	WISARD_NULL_PEDCONTENTS,
	WISARD_NULL_IID_IN_DATASET,
	WISARD_NULL_FEATURE,
	WISARD_NULL_GENESET_CONTENTS,
	WISARD_NULL_FILEWITHDESC,

	WISARD_SYST_FAIL_ACTV_SIGHDLR,
	WISARD_SYST_FAIL_GET_FQLS,

	WISARD_SYST_NULL_MDR_RESULT,
	WISARD_SYST_NULL_PHENOTYPE,
	WISARD_SYST_NULL_PARENT,
	WISARD_SYST_NULL_CLPGENO,
	WISARD_SYST_NULL_IO,
	WISARD_SYST_NULL_FAMSTRUCT,
	WISARD_SYST_NULL_IMPLEMENT,
	WISARD_SYST_NULL_VRTDATA_OF_SAMP,
	WISARD_SYST_NULL_INPUT,
	WISARD_SYST_NULL_REGRGENO,
	WISARD_SYST_NULL_DATA,
	WISARD_SYST_NULL_CHILDGENO,

	WISARD_SYST_DUPL_LONGOPT,
	WISARD_SYST_DUPL_SHORTOPT,

	WISARD_SYST_INVL_MDR_CELLRANGE,
	WISARD_SYST_INVL_OPTCLASSINIT,
	WISARD_SYST_INVL_OPTTYPE,
	WISARD_SYST_INVL_BOOLOPTVAL,
	WISARD_SYST_INVL_NUMBEROPTVAL,
	WISARD_SYST_INVL_REALOPTVAL,
	WISARD_SYST_INVL_LOGSTREAM,
	WISARD_SYST_INVL_NOSUCHOPT,
	WISARD_SYST_INVL_COVTYPE,
	WISARD_SYST_INVL_PEDCMCIDX,
	WISARD_SYST_INVL_SNPDATA,
	WISARD_SYST_INVL_MEMTYPE,
	WISARD_SYST_INVL_PEDGENO,
	WISARD_SYST_INVL_ALTCORSIZE,
	WISARD_SYST_INVL_GS_SZSAMP,
	WISARD_SYST_INVL_DISTFUNC_CMD,
	WISARD_SYST_INVL_DIM,
	WISARD_SYST_INVL_DIM_MSG,
	WISARD_SYST_INVL_GSTYPE,
	WISARD_SYST_INVL_MISSFOUNDER,
	WISARD_SYST_INVL_EXPCOUNT,
	WISARD_SYST_INVL_SAMEORMULTSZ,
	WISARD_SYST_INVL_REPTYPE,
	WISARD_SYST_INVL_RANGEDEF,

	WISARD_SYST_CANT_RUNANALYSIS,
	WISARD_SYST_CANT_MISGENO,
	WISARD_SYST_CANT_REPREPLACE,
} xErrType;

typedef struct _xErrMsg {
	xErrType	X_errType;
	const char	*S_errCode;
	const char	*S_errMsg;
} xErrMsg;

typedef enum _xRepType {
	WISARD_REP_NONE,
	WISARD_REP_DATESTART,	/* Date begin */
	WISARD_REP_DATEEND,		/* Date ended */
	WISARD_REP_TOTEXECTIME, /* Total execution time */
	WISARD_REP_RANDSEED,	/* Random seed */
	WISARD_REP_CACT,		/* Case-control coding */
	WISARD_REP_CORSTR,		/* Correlation structure type */
	WISARD_REP_OUTPREFIX,	/* Output prefix */
	WISARD_REP_WARNING,		/* General warnings */
	WISARD_REP_EXECRES,		/* Execution result */
	WISARD_REP_TOTNCORE,	/* Total number of installed cores */
	WISARD_REP_USEDNCORE,	/* Total number of used cores */
} xRepType;

typedef struct _xRepData {
	const xRepType	X_type;
	const char		*name;
	const char		S_fmt;
	const char		B_replaceable;	/* Is the value can be changed after initialization? */
	timeval			X_time;
	struct _xRepData
					*Xa_subs;
	char			*Sp_val;
} xRepData;

/* Internal usage */
void __halt(const char *S_func, const char *S_fmt, ...);
void __halt_format(const char *S_func, xErrType X_errType, ...);
void __notice(int B_logThis, const char *S_func, const char *S_fmt, ...);
void __output(const char *S_prefix, unsigned char N_color, int B_logThis,
	const char *S_func, const char *S_prtFormat, ...);
void __verbose(int B_logThis, const char *S_func, const char *S_fmt, ...);

/* 로그 타입 */
enum LOGTYPE {
	/* stderr에 출력, 출력 스트림이 stdout이 아니면 그곳에도 출력 */
	L_FORCELOG,

	/* stderr에만 출력 */
	L_NOTICE,

	/* 출력 스트림이 어느쪽이든지 출력 */
	L_LOG,

	/* 출력 스트림이 stdout이 아닐 때만 출력 */
	L_LOGONLY,

	/* Ignore all log but halt */
	L_NOLOG,
};

/* 로깅 함수 */
int		__log(LOGTYPE X_logType, const char *S_prtFormat, ...);
int		logSet(LOGTYPE X_logType, const char *S_FN);
int		logSet(LOGTYPE X_logType, FILE *H_target);
void	logEnd();

/* Report functions */
void	report(xRepType X_rtype, const char *S_val);
void	report(xRepType X_rtype, wsUint N_val);
void	report(xRepType X_rtype);
void	reportEnd();

} // End namespace ONETOOL

#endif
