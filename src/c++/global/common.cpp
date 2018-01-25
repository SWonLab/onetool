#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sstream>
#include "global/common.h"
#include "utils/util.h"
#include "global/option.h"
#include "utils/stack.h"
#include "global/Rconn.h"

namespace ONETOOL {

const wsReal WISARD_NAN		= numeric_limits<wsReal>::quiet_NaN();
const wsReal WISARD_INF		= numeric_limits<wsReal>::infinity();
const wsReal WISARD_PI		= REAL_CONST(3.141592653589793238462643383279502884197169399375);
const wsReal WISARD_NA_REAL	= REAL_CONST(-9);
const wsReal W0				= REAL_CONST(0.0);
const wsReal W1				= REAL_CONST(1.0);
const wsReal W2				= REAL_CONST(2.0);


xErrMsg Xa_halt[] = {
	{ WISARD_DUPL_SNPINGENESET,          "DUPL_SNPINGENESET",          "Variant [%s] is already exists in same gene set [%s]" },
	{ WISARD_DUPL_GSID,                  "DUPL_GSID",                  "Set identifier [%s] is duplicated [%d] line" },
	{ WISARD_DUPL_OPTION,                "DUPL_OPTION",                "Option [%s] allocated twice, check the command line!" },

	{ WISARD_CANT_DO_W_SOMEOPT,          "CANT_DO_W_SOMEOPT",          "[%s] is not allowed with [%s] option" },
	{ WISARD_CANT_DO_WO_SOMEOPT,         "CANT_DO_WO_SOMEOPT",         "[%s] is not allowed without [%s] option" },
	{ WISARD_CANT_OPT_W_OTHEROPT,        "CANT_OPT_W_OTHEROPT",        "[%s] option is not allowed with [%s] option" },
	{ WISARD_CANT_OPT_W_SOMESTATE,       "CANT_OPT_W_SOMESTATE",       "[%s] option is not allowed when [%s]" },
/*
	
	
	"--autoonly", "the sex chromosome is included in --chr"
*/
	{ WISARD_CANT_MULTILOCUS_WO_INDEL,   "CANT_MULTILOCUS_WO_INDEL",   "Multiloci genotype found at line [%d] from input [%s], can't get without --indel" },
	{ WISARD_CANT_ONEALLELE_MISSING,     "CANT_ONEALLELE_MISSING",     "Genotype of variant [%s] in line [%d] have one-side missing, which is not allowed!" },
	{ WISARD_CANT_ONESIDE_PARENT,        "CANT_ONESIDE_PARENT",        "Sample [%s:%s] have only one parent [%s] and [%s], it is not allowed" },
	{ WISARD_CANT_NONAUTOCHR_W_AUTOONLY, "CANT_NONAUTOCHR_W_AUTOONLY", "--autoonly enabled but sex chromosome is activated with --chr" },
	{ WISARD_CANT_MQLS_WO_BLUPORPREV,    "CANT_MQLS_WO_BLUPORPREV",    "One of --blup/--prevalence is required in order perform --mqls" },
	{ WISARD_CANT_MQLS_W_BOTHBLUPPREV,   "CANT_MQLS_W_BOTHBLUPPREV",   "Both of --blup and --prevalence cannot be assigned concurrently" },
	{ WISARD_CANT_BINPHENO_W_LONGIANA,   "CANT_BINPHENO_W_LONGIANA",   "Dichotomous phenotype cannot be used in longitudinal analysis, at [%d]th phenotype" },
	{ WISARD_CANT_PEDCMC_W_MULTIPHENO,   "CANT_PEDCMC_W_MULTIPHENO",   "PEDCMC method only can be applied to single phenotype, [%d] phenotype assigned" },
/*D*/{ WISARD_CANT_PEDCMC_W_CONTPHENO,    "CANT_PEDCMC_W_CONTPHENO",    "PEDCMC method only can be applied to dichomotous phenotype" },
	{ WISARD_CANT_MEDCOR_W_CORPEARSON,   "CANT_MEDCOR_W_CORPEARSON",   "Cannot using median to get empirical correlation with Pearson's correlation" },
	{ WISARD_CANT_FQLS_W_BLUP,           "CANT_FQLS_W_BLUP",           "Cannot perform FQLS analysis with BLUP subtraction option" },
	{ WISARD_CANT_FQLS_W_CONTPHENO,      "CANT_FQLS_W_CONTPHENO",      "Cannot perform FQLS with continuous phenotype" },
	{ WISARD_CANT_GETINDEL_WO_INDEL,     "CANT_GETINDEL_WO_INDEL",     "Indel variant found at [%d] line, but cannot proceed without --indel" },
	{ WISARD_CANT_EXCL_OPT,              "CANT_EXCL_OPT",              "[%s] option and [%s] option are mutually exclusive!" },

	{ WISARD_FAIL_OPEN_URLWITHDESC,	     "FAIL_OPEN_URLWITHDESC",      "Failed to open [%s] URL [%s]" },
	{ WISARD_FAIL_OPEN_URL, 	         "FAIL_OPEN_URL",              "Failed to open URL [%s]" },
	{ WISARD_FAIL_OPEN_FILEWITHDESC,     "FAIL_OPEN_FILEWITHDESC",     "Failed to open [%s] file [%s]" },
	{ WISARD_FAIL_OPEN_FILE, 	         "FAIL_OPEN_FILE",             "Failed to open file [%s]" },
	{ WISARD_FAIL_OPEN_GZWITHDESC,       "FAIL_OPEN_GZWITHDESC",       "Failed to open [%s] gzipped file [%s], not supported" },
	{ WISARD_FAIL_OPEN_GZ, 	             "FAIL_OPEN_GZ",               "Failed to open gzipped file [%s], not supported" },
	{ WISARD_FAIL_AUTOOPEN_TFAMEXTMISS,  "FAIL_AUTOOPEN_TFAMEXTMISS",  "Failed to open TFAM file since TPED [%s] does not have extension" },
	{ WISARD_FAIL_AUTOOPEN_TFAMNOTEXIST, "FAIL_AUTOOPEN_TFAMNOTEXIST", "Failed to open corresponding TFAM file of TPED [%s]" },
	{ WISARD_FAIL_OPEN_FILE, 	   	   	 "FAIL_OPEN_FILE",             "Failed to open gzipped file [%s], not supported" },
	{ WISARD_FAIL_INIT_WINSOCK, 	   	 "FAIL_INIT_WINSOCK",          "Failed to initialize WinSock [%d.%d]" },
	{ WISARD_FAIL_INIT_LOG,	   	   	     "FAIL_INIT_LOG",              "Failed to open log file [%s]" },
	{ WISARD_FAIL_SOCKCREATE, 	   	     "FAIL_SOCKCREATE",            "Failed to create socket for URL [%s]" },
	{ WISARD_FAIL_SOCKCONN, 	   	   	 "FAIL_SOCKCONN",              "Failed to connect to URL [%s]" },
	{ WISARD_FAIL_OPEN_ANNOSNP, 	   	 "FAIL_OPEN_ANNOSNP",          "Failed to open variant annotation file [%s]" },
	{ WISARD_FAIL_MEMALLOC_W_PROC,       "FAIL_MEMALLOC_W_PROC",       "Failed to allocate memory [%d] byte(s)" },
	{ WISARD_FAIL_READ_BEDSIG,           "FAIL_READ_BEDSIG",           "Failed to read the signature of input BED file" },
	{ WISARD_FAIL_OPEN_RUNAS,            "FAIL_OPEN_RUNAS",            "Failed to open runas description file [%s]" },
	{ WISARD_FAIL_ESTIMATE_PARAM,        "FAIL_ESTIMATE_PARAM",        "Failed to estimate parameter [%s]" },
	{ WISARD_FAIL_FIT_NULLMODEL,         "FAIL_FIT_NULLMODEL",         "Failed to fitting NULL model" },
	{ WISARD_FAIL_MEMALLOC_W_SYSMSG,     "FAIL_MEMALLOC_W_SYSMSG",     "Failed to allocate [%d] bytes from [%s] at [%d] line\n System [%s]" },
	{ WISARD_FAIL_MEMALLOC,              "FAIL_MEMALLOC",              "Failed to allocate [%d] bytes from [%s] at [%d] line" },
	{ WISARD_FAIL_OPEN_GSFILE_IDX,       "FAIL_OPEN_GSFILE_IDX",       "Failed to open [%d] th gene-set definition file [%s]" },
	{ WISARD_FAIL_LAUNCH_R,              "FAIL_LAUNCH_R",              "Failed to launch [%s] [Return code : %d]" },

	{ WISARD_INVL_VARIANTPOS,	   	     "INVL_VARIANTPOS",            "Variant [%s] have invalid position [%s]" },
	{ WISARD_INVL_ALLELEFORMAT,	   	     "INVL_ALLELEFORMAT", 		   "[%d] th variant in [%d] th sample is incorrect allele [%s]" },
	{ WISARD_INVL_NVRTOFSAMPLE, 	   	 "INVL_NVRTOFSAMPLE",          "Sample [%s:%s] expected to have [%d/%d] variants in file/after initial filtering, but found [%d/%d]" },
	{ WISARD_INVL_NSAMPLEOFSNP, 	   	 "INVL_NSAMPLEOFSNP",          "Variant [%s] expected to have [%d] samples but found [%d] samples" },
	{ WISARD_INVL_OPTNAME, 	   	   	     "INVL_OPTNAME",               "Argument [%s] is not valid option name" },
	{ WISARD_INVL_OPTVAL_SC,	         "INVL_OPTVAL_SC",             "Argument [%s] is not valid option name, maybe there is a space between [%s] and [%s]?" },
	{ WISARD_INVL_OPTVAL_GENERAL_INT,    "INVL_OPTVAL_GENERAL_INT",    "Integer value [%d] for option [%s] is not valid" },
	{ WISARD_INVL_OPTVAL_GENERAL_STR,    "INVL_OPTVAL_GENERAL_STR",    "String value [%s] for option [%s] is not valid" },
	{ WISARD_INVL_LONGIVARDIM,	   	     "INVL_LONGIVARDIM",           "Invalid phenotype variance-covariance matrix dimension[%d*%d], should be [%d*%d]" },
	{ WISARD_INVL_WSOCKVER, 	   	   	 "INVL_WSOCKVER",              "WinSock initialized [%d.%d] but requested version was [%d.%d]" },
	{ WISARD_INVL_URLFORMAT,	   	     "INVL_URLFORMAT",             "Invalid URL format [%s]" },
	{ WISARD_INVL_ANNOGENE_FORMAT,	     "INVL_ANNOGENE_FORMAT",       "Invalid gene annotation definition format at line [%d], expected to have [%d] columns but [%d] columns found" },
	{ WISARD_INVL_SNPLIST_FORMAT,        "INVL_SNPLIST_FORMAT",        "Invalid format for variant list file [%s] at line [%d]" },
	{ WISARD_INVL_SAMPLIST_FORMAT,       "INVL_SAMPLIST_FORMAT",       "Invalid format for sample list file [%s] at line [%d]" },
	{ WISARD_INVL_PEDCONTENTS,           "INVL_PEDCONTENTS",           "Unexpected end of line found at line [%d] from input [%s]" },
	{ WISARD_INVL_CHROMFORMAT,           "INVL_CHROMFORMAT",           "Invalid chromosome definition [%s] found" },
	{ WISARD_INVL_BEDVERSION,            "INVL_BEDVERSION",            "Unsupported or invalid BED file" },
	{ WISARD_INVL_BEDSIZE,               "INVL_BEDSIZE",               "BED size is %s than expected, possibly error? [%d variants*%d samples]" },
	{ WISARD_INVL_NPREV_NPHENO,          "INVL_NPREV_NPHENO",          "Number of given prevalence[%d] should be matched with the number of phenotypes[%d]" },
	{ WISARD_INVL_CORFORMAT,             "INVL_CORFORMAT",             "Invalid alternative correlation format at line [%d], found [%d] columns but expected [%d] columns" },
	{ WISARD_INVL_ALTCORSIZE,            "INVL_ALTCORSIZE",            "Correlation alternation file [%s] have not match size [%d] with expected size [%d]" },
	{ WISARD_INVL_CORFILE_MORESZSAMP,    "INVL_CORFILE_MORESZSAMP",    "[%d]th line in correlation alternation file have more samples than expected ([%d] expected)" },
	{ WISARD_INVL_CORFILE_LESSSZSAMP,    "INVL_CORFILE_LESSSZSAMP",    "[%d]th line in correlation alternation file have less samples than expected ([%d] expected)" },
	{ WISARD_INVL_GSFMT_NONPOS_GSRANGE,  "INVL_GSFMT_NONPOS_GSRANGE",  "Range definition [%d~%d] have negative value at [%d] line" },
	{ WISARD_INVL_GSFMT_REV_GSRANGE,     "INVL_GSFMT_REV_GSRANGE",     "Range definition [%d~%d] have reverse range at [%d] line" },
	{ WISARD_INVL_GSFMT_CHR,             "INVL_GSFMT_CHR",             "Chromosome definition [%s] is not valid at [%d] line" },
	{ WISARD_INVL_GSFMT_DBLEND,          "INVL_GSFMT_DBLEND",          "END identifier used double times at [%d] line" },
	{ WISARD_INVL_BIMFORMAT,             "INVL_BIMFORMAT",             "Assigned BIM file have invalid format at [%d] line" },
	{ WISARD_INVL_BEDSIGNATURE,          "INVL_BEDSIGNATURE",          "Invalid BED file signature, cannot read" },
	{ WISARD_INVL_PLINK_INPUT,           "INVL_PLINK_INPUT",           "Invalid PLINK input [%s], please check the assigned file" },
	{ WISARD_INVL_SAMPNAME,              "INVL_SAMPNAME",              "Sample [%s::%s] at line [%d] exists in LGEN but not exists in FAM file" },
	{ WISARD_INVL_NVRTMAP_W_NVRTPED,     "INVL_NVRTMAP_W_NVRTPED",     "Variant count in MAP file [%d] is different in variant count in PED file [%d]" },
	{ WISARD_INVL_VCFMETA,               "INVL_VCFMETA",               "VCF meta data error in line [%d] : %s" },
	{ WISARD_INVL_VCFMANDTRY_COLNAME,    "INVL_VCFMANDTRY_COLNAME",    "[%d]th mandatory column of VCF file should be [%s], but [%s] found" },
	{ WISARD_INVL_PEDCHAR,               "INVL_PEDCHAR",               "PED file [%s] seems [%s] coding, but %s[%c] found" },
	{ WISARD_INVL_PEDMISCHAR,            "INVL_PEDMISCHAR",            "PED file [%s] seems [%s] coding, but %s[%s] found" },
	{ WISARD_INVL_RANGE_NUMBER,          "INVL_RANGE_NUMBER",          "Option [%s] requires range %s, but the parameter is [%d]" },
	{ WISARD_INVL_RANGE_REAL,            "INVL_RANGE_REAL",            "Option [%s] requires range %s, but the parameter is [%g]" },
	{ WISARD_INVL_RANGE_STRING,          "INVL_RANGE_STRING",          "Option [%s] requires range %s, but the parameter is [%s]" },
	{ WISARD_INVL_FILE_GENERAL,          "INVL_FILE_GENERAL",          "[%s] file [%s] error : %s" },
	{ WISARD_INVL_FILE_INCMPL,           "INVL_FILE_INCMPL",           "[%s] file [%s] have incomplete data in line [%d]" },
	{ WISARD_INVL_FILE_INVALID,          "INVL_FILE_INVALID",          "[%s] file [%s] have invalid data [%s] in line [%d]" },
	{ WISARD_INVL_FILE_INVALID_DESC,     "INVL_FILE_INVALID_DESC",     "[%s] file [%s] have invalid %s [%s] in line [%d]" },
	{ WISARD_INVL_EOF_LINE,              "INVL_EOF_LINE",              "[%s] file [%s] expected to have [%d] lines, but only [%d] lines found" },
	{ WISARD_INVL_OPERATOR,              "INVL_OPERATOR",              "Incompatible operator [%s] was assigned in %s action" },
	
/**/{ WISARD_NULL_GENEDEF,    	      	 "NULL_GENEDEF",               "No gene definition left after matching with dataset" },
/**/{ WISARD_NULL_SAMPLE,    	      	 "NULL_SAMPLE",                "No sample left to analysis in dataset" },
/**/{ WISARD_NULL_OPTVALUE,    	      	 "NULL_OPTVALUE",              "Option [%s] requires value, but no value found" },
	{ WISARD_NULL_ESSENTIALOPT,    	     "NULL_ESSENTIALOPT",          "Option [%s] is essential in order to continue, but not assigned" },
	{ WISARD_NULL_ESSENTIALOPT_OR, 	     "NULL_ESSENTIALOPT_OR",       "Option [%s] or [%s] is essential in order to continue, but not assigned" },
	{ WISARD_NULL_BLUP,    	      	     "NULL_BLUP", 				   "BLUP should be assigned since --prevalence is not assigned, but not given" },
	{ WISARD_NULL_ANNOGENENAME,   	     "NULL_ANNOGENENAME",          "Gene name of annotation definition is EMPTY at line [%d]" },
	{ WISARD_NULL_VTCOVAR,   	      	 "NULL_VTCOVAR",               "VT method does not allows covariate missing but missing found at %dth cov, sample [%s:%s]" },
	{ WISARD_NULL_PEDCONTENTS,           "NULL_PEDCONTENTS",           "PED file [%s] have no contents" },
	{ WISARD_NULL_IID_IN_DATASET,        "NULL_IID_IN_DATASET",        "IID [%s or %s] in line [%d] is not in the given dataset" },
	{ WISARD_NULL_FEATURE,               "NULL_FEATURE",               "Program tried to fetch feature [%s], which is not supported in this version" },
	{ WISARD_NULL_GENESET_CONTENTS,      "NULL_GENESET_CONTENTS",      "Given gene-set file have empty content!" },
	{ WISARD_NULL_FILEWITHDESC,          "NULL_FILEWITHDESC",          "[%s] file [%s] have no contents" },

	{ WISARD_SYST_FAIL_ACTV_SIGHDLR,     "SYST_FAIL_ACTV_SIGHDLR",     "Cannot activate Unix signal handler to control unexpected exception" },
	{ WISARD_SYST_FAIL_GET_FQLS,         "SYST_FAIL_GET_FQLS",         "Family [%s] error while getting FQLS %s statistics" },

	{ WISARD_SYST_NULL_MDR_RESULT,   	 "SYST_NULL_MDR_RESULT",       "SYSERR: Null result class" },
	{ WISARD_SYST_NULL_PHENOTYPE,    	 "SYST_NULL_PHENOTYPE", 	   "SYSERR : [%s:%s] included despite of its phenotype missing" },
	{ WISARD_SYST_NULL_PARENT,    	     "SYST_NULL_PARENT", 		   "SYSERR : %s of [%s:%s] is missing but %s is non-missing, it is not permitted" },
	{ WISARD_SYST_NULL_CLPGENO,   	     "SYST_NULL_CLPGENO",          "SYSERR : Collapsing genotype array is NULL" },
	{ WISARD_SYST_NULL_IO,    	      	 "SYST_NULL_IO",               "SYSERR : IO system is NULL" },
	{ WISARD_SYST_NULL_FAMSTRUCT,        "SYST_NULL_FAMSTRUCT",        "SYSERR : Family structure analysis required but not performed" },
	{ WISARD_SYST_NULL_IMPLEMENT,        "SYST_NULL_IMPLEMENT",        "SYSERR : [%s] feature is not implemented but called" },
	{ WISARD_SYST_NULL_VRTDATA_OF_SAMP,  "SYST_NULL_VRTDATA_OF_SAMP",  "SYSERR : Sample [%s:%s] is non-missing, but the data was not found" },
	{ WISARD_SYST_NULL_INPUT,            "SYST_NULL_INPUT",            "SYSERR : Give input [%s] is NULL pointer" },
	{ WISARD_SYST_NULL_REGRGENO,         "SYST_NULL_REGRGENO",         "SYSERR : Genotype of variant [%s] of [%s:%s] cannot be MISSING in regression" },
	{ WISARD_SYST_NULL_DATA,             "SYST_NULL_DATA",             "SYSERR : A pointer of data for [%s] should exists, but NULL given" },
	{ WISARD_SYST_NULL_CHILDGENO,        "SYST_NULL_CHILDGENO",        "SYSERR : Child [%s:%s] of parent [%s:%s] does not have data" },

	{ WISARD_SYST_DUPL_LONGOPT,   	     "SYST_DUPL_LONGOPT",          "SYSERR : Long option definition [%s] duplicated" },
	{ WISARD_SYST_DUPL_SHORTOPT,   	     "SYST_DUPL_SHORTOPT",         "SYSERR : Short option definition [%s] duplicated" },

	{ WISARD_SYST_INVL_MDR_CELLRANGE,    "SYST_INVL_MDR_CELLRANGE",    "SYSERR : MDR cell have value %d, but should not over %d" },
	{ WISARD_SYST_INVL_OPTCLASSINIT,     "SYST_INVL_OPTCLASSINIT", 	   "SYSERR : Initiation of option class called twice!" },
	{ WISARD_SYST_INVL_OPTTYPE,    	     "SYST_INVL_OPTTYPE",          "SYSERR : Option type [%d] is not valid for option [%s]" },
	{ WISARD_SYST_INVL_BOOLOPTVAL,       "SYST_INVL_BOOLOPTVAL",       "SYSERR : Invalid option parameter [%d] for ON/OFF option [%s]" },
	{ WISARD_SYST_INVL_NUMBEROPTVAL,     "SYST_INVL_NUMBEROPTVAL",     "SYSERR : Invalid option parameter [%s] for integer option [%s]" },
	{ WISARD_SYST_INVL_REALOPTVAL,       "SYST_INVL_REALOPTVAL",       "SYSERR : Invalid option parameter [%s] for real option [%s]" },
	{ WISARD_SYST_INVL_LOGSTREAM,        "SYST_INVL_LOGSTREAM",        "SYSERR : Invalid logging stream [%d] to print message [%s]" },
	{ WISARD_SYST_INVL_NOSUCHOPT,        "SYST_INVL_NOSUCHOPT",        "SYSERR : No such option [%s] is defined" },
	{ WISARD_SYST_INVL_COVTYPE,          "SYST_INVL_COVTYPE",          "SYSERR : Type of covariate [%s] is undefined" },
	{ WISARD_SYST_INVL_PEDCMCIDX,        "SYST_INVL_PEDCMCIDX",        "SYSERR : Too large colIdx on PEDCMC %d>%d" },
	{ WISARD_SYST_INVL_LOGSTREAM,        "SYST_INVL_LOGSTREAM"         "SYSERR : Invalid log stream [%d] given for printing message [%s]" },
	{ WISARD_SYST_INVL_SNPDATA,          "SYST_INVL_SNPDATA",          "SYSERR : Genotype data of [%s] for sample [%s:%s] have invalid value [%d]" },
	{ WISARD_SYST_INVL_MEMTYPE,          "SYST_INVL_MEMTYPE",          "SYSERR : Can't free memory, invalid memory type [%d]" },
	{ WISARD_SYST_INVL_PEDGENO,          "SYST_INVL_PEDGENO",          "SYSERR : Parsed genotype [%d] is invalid in variant [%s] of line [%d]" },
	{ WISARD_SYST_INVL_ALTCORSIZE,       "SYST_INVL_ALTCORSIZE",       "SYSERR : Retrieved number of samples [%d] in alternative correlation file [%s] is larger than expected [%d]" },
	{ WISARD_SYST_INVL_GS_SZSAMP,        "SYST_INVL_GS_SZSAMP",        "SYSERR : Retrieved # of samples is higher than expected number, logic error" },
	{ WISARD_SYST_INVL_DISTFUNC_CMD,     "SYST_INVL_DISTFUNC_CMD",     "SYSERR : Unimplemented distribution command [%d]" },
	{ WISARD_SYST_INVL_DIM,              "SYST_INVL_DIM",              "SYSERR : Invalid dimension match, source wants [%d] but target have [%d]" },
	{ WISARD_SYST_INVL_DIM_MSG,          "SYST_INVL_DIM_MSG",          "SYSERR : Dimension error on [%s], source wants [%d] but target have [%d]" },
	{ WISARD_SYST_INVL_GSTYPE,           "SYST_INVL_GSTYPE",           "SYSERR : Given gene-set [%s] classified invalid gene-set type [%d]" },
	{ WISARD_SYST_INVL_MISSFOUNDER,      "SYST_INVL_MISSFOUNDER",      "SYSERR : Sample [%s:%s] should be missing founder but not missing" },
	{ WISARD_SYST_INVL_EXPCOUNT,         "SYST_INVL_EXPCOUNT",         "SYSERR : %s[%d] must equal to %s[%d], but not equal!" },
	{ WISARD_SYST_INVL_SAMEORMULTSZ,     "SYST_INVL_SAMEORMULTSZ",     "SYSERR : %s[%d] must equal to %s[%d] or its exact multiply!" },
	{ WISARD_SYST_INVL_REPTYPE,          "SYST_INVL_REPTYPE",          "SYSERR : Report entry [%s] requires [%s], but string assgiend!" },
	{ WISARD_SYST_INVL_RANGEDEF,         "SYST_INVL_RANGEDEF",         "SYSERR : Option [%s] have incompatible range definition [%d]" },

	{ WISARD_SYST_CANT_RUNANALYSIS,      "SYST_CANT_RUNANALYSIS",      "SYSERR : This analysis NOT runnable" },
	{ WISARD_SYST_CANT_MISGENO,          "SYST_CANT_MISGENO",          "SYSERR : This analysis does not allow missing genotype, but variant [%s] of sample [%s] with missing accessed" },
	{ WISARD_SYST_CANT_REPREPLACE,       "SYST_CANT_REPREPLACE",       "SYSERR : Report entry [%s] does not allow replacement, a value [%s] was already assigned!" },
};

/* 출력 스트림 */
FILE	*H_logFile;
char	*S_buf	= NULL;
char	*Sp_buf	= NULL;
char	B_noLog	= 0;

void __halt_main(const char *S_errCode, const char *S_func, const char *S_prtFormat, va_list *args)
{
	extern cTimer t;
	FILE	*H_targetFile = stderr;
	char	S_buf[4096];

	vsprintf(S_buf, S_prtFormat, *args);

	fprintf(H_targetFile, "-------------------------------------------------------------------------------\n"
		" *** HALTED in [%s], after [%s] from start\n"
		"-------------------------------------------------------------------------------\n", S_func, t.getReadable());
	if (S_errCode)
		fprintf(H_targetFile, " ERROR CODE : %s\n"
			"-------------------------------------------------------------------------------\n", S_errCode);
	fwrite(S_buf, 1, strlen(S_buf), H_targetFile);
	fprintf(H_targetFile, "\n-------------------------------------------------------------------------------\n");

	if (H_logFile && H_logFile != stderr) {
		if (H_logFile == (FILE *)0xffff) {
			if (Sp_buf) {
				Sp_buf += sprintf(Sp_buf, "-------------------------------------------------------------------------------\n"
					" *** HALTED in [%s], after [%s] from start\n"
					"-------------------------------------------------------------------------------\n", S_func, t.getReadable());
				if (S_errCode)
					Sp_buf += sprintf(Sp_buf, " ERROR CODE : %s\n"
						"-------------------------------------------------------------------------------\n", S_errCode);
				//va_start(H_varList, S_prtFormat);
				Sp_buf += strlen(strcat(Sp_buf, S_buf));
				//Sp_buf += sprintf(Sp_buf, S_buf);
				Sp_buf += sprintf(Sp_buf, "\n-------------------------------------------------------------------------------\n");
			}// else {
			// 				printf("-------------------------------------------------------------------------------\n"
			// 					" *** HALTED in [%s]! ***\n"
			// 					"-------------------------------------------------------------------------------\n", S_func);
			// 				va_start(H_varList, S_prtFormat);
			// 				vprintf(S_prtFormat, H_varList);
			// 				printf("\n-------------------------------------------------------------------------------\n");
			// 			}
		} else {
			fprintf(H_logFile, "-------------------------------------------------------------------------------\n"
				" *** HALTED in [%s], after [%s] from start\n"
				"-------------------------------------------------------------------------------\n", S_func, t.getReadable());
			if (S_errCode)
				fprintf(H_logFile, " ERROR CODE : %s\n"
					"-------------------------------------------------------------------------------\n", S_errCode);
			//va_start(H_varList, S_prtFormat);
			fwrite(S_buf, 1, strlen(S_buf), H_logFile);
			//fprintf(H_logFile, S_buf);
			fprintf(H_logFile, "\n-------------------------------------------------------------------------------\n");
		}
	}
#ifdef _WIN32
//	MyStackWalker mw;
//	mw.ShowCallstack();
#else
#endif

	report(WISARD_REP_TOTEXECTIME, t.getReadable());
	report(WISARD_REP_EXECRES, "Failed");
	logEnd();
	//getc(stdin);
	exit(1);
}
	
void __halt_format(const char *S_func, xErrType X_errType, ...)
{
	va_list args;
	va_start (args, X_errType);
	wsUint N_sz = (wsUint)(len(Xa_halt, xErrMsg));
	for (wsUint i=0 ; i<N_sz ; i++)
		if (Xa_halt[i].X_errType == X_errType) {
			__halt_main(Xa_halt[i].S_errCode, S_func, Xa_halt[i].S_errMsg, &args);
			break;
		}
		va_end(args);
}

void __halt(const char *S_func, const char *S_prtFormat, ...)
{
	va_list	H_varList;
	va_start(H_varList, S_prtFormat);
	__halt_main(NULL, S_func, S_prtFormat, &H_varList);
	va_end(H_varList);
}

void __notice(int B_logThis, const char *S_func, const char *S_prtFormat, ...)
{
	va_list	H_varList;
	FILE*	H_targetFile = NULL;

	/* Set log type */
	LOGTYPE	X_logType = B_noLog ? L_NOLOG : (B_logThis ? L_FORCELOG : L_NOTICE);

	switch (X_logType) {
	case L_NOLOG:    H_targetFile = NULL;		break;
	case L_FORCELOG: H_targetFile = stderr;		break;
	case L_NOTICE:	 H_targetFile = stderr;		break;
	case L_LOG:		 H_targetFile = H_logFile;	break;
	case L_LOGONLY:	 H_targetFile = H_logFile;	break;
	default:
		halt_fmt(WISARD_SYST_INVL_LOGSTREAM, X_logType, S_prtFormat);
		break;
	}
	if (H_targetFile == NULL)
		return;

	va_start(H_varList, S_prtFormat);
	if ((X_logType == L_LOGONLY && H_logFile != stdout) ||
		X_logType != L_LOGONLY) {
		if (S_func) {
			if (H_targetFile == stdout || H_targetFile == stderr) {
#ifdef _WIN32
				HANDLE hstdout = GetStdHandle( STD_OUTPUT_HANDLE );
				CONSOLE_SCREEN_BUFFER_INFO csbi;
				GetConsoleScreenBufferInfo( hstdout, &csbi );
				SetConsoleTextAttribute( hstdout, 0x0f | csbi.wAttributes&0xf0 );
				fprintf(H_targetFile, "[%s] ", S_func);
				SetConsoleTextAttribute( hstdout, csbi.wAttributes );
#else
				fprintf(H_targetFile, "\033[4m[%s]\033[0m ", S_func);
#endif
			} else
				fprintf(H_targetFile, "[%s] ", S_func);
		}
		vfprintf(H_targetFile, S_prtFormat, H_varList);
	}
	if (X_logType == L_FORCELOG && H_logFile && H_logFile != stderr) {
		if (H_logFile == (FILE *)0xffff) {
			if (Sp_buf) {
				if (S_func)
					Sp_buf += sprintf(Sp_buf, "[%s] ", S_func);
				va_start(H_varList, S_prtFormat);
				Sp_buf += vsprintf(Sp_buf, S_prtFormat, H_varList);
			}// else {
// 				printf("[%s] ", S_func);
// 				va_start(H_varList, S_prtFormat);
// 				vprintf(S_prtFormat, H_varList);
// 			}
		} else {
			if (S_func)
				fprintf(H_logFile, "[%s] ", S_func);
			va_start(H_varList, S_prtFormat);
			vfprintf(H_logFile, S_prtFormat, H_varList);
		}
	}
	va_end(H_varList);
}

void __output(const char *S_prefix, unsigned char N_color, int B_logThis,
	const char *S_func, const char *S_prtFormat, ...)
{
	va_list	H_varList;
	FILE*	H_targetFile = NULL;

	/* Set log type */
	LOGTYPE	X_logType = B_noLog ? L_NOLOG : (B_logThis ? L_FORCELOG : L_NOTICE);

	switch (X_logType) {
	case L_NOLOG:    H_targetFile = NULL;		break;
	case L_FORCELOG: H_targetFile = stderr;		break;
	case L_NOTICE:	 H_targetFile = stderr;		break;
	case L_LOG:		 H_targetFile = H_logFile;	break;
	case L_LOGONLY:	 H_targetFile = H_logFile;	break;
	default:
		halt_fmt(WISARD_SYST_INVL_LOGSTREAM, X_logType, S_prtFormat);
		break;
	}
	if (H_targetFile == NULL)
		return;

	va_start(H_varList, S_prtFormat);
	if ((X_logType == L_LOGONLY && H_logFile != stdout) ||
		X_logType != L_LOGONLY) {
		if (S_func) {
			if (H_targetFile == stdout || H_targetFile == stderr) {
#ifdef _WIN32
				HANDLE hstdout = GetStdHandle( STD_OUTPUT_HANDLE );
				CONSOLE_SCREEN_BUFFER_INFO csbi;
				GetConsoleScreenBufferInfo( hstdout, &csbi );
				SetConsoleTextAttribute( hstdout, 0x0f | csbi.wAttributes&0xf0 );
				fprintf(H_targetFile, "[%s] ", S_func);
				SetConsoleTextAttribute( hstdout, N_color | csbi.wAttributes&0xf0 );
				fprintf(H_targetFile, "[%s] ", S_prefix);
				SetConsoleTextAttribute( hstdout, csbi.wAttributes );
#else
				fprintf(H_targetFile, "\033[4m[%s]", S_func);
				fprintf(H_targetFile, "\033[0;%dm [%s]\033[0m ", N_color, S_prefix);
#endif
			} else
				fprintf(H_targetFile, "[%s] [%s] ", S_func, S_prefix);
		}
		vfprintf(H_targetFile, S_prtFormat, H_varList);
	}
	if (X_logType == L_FORCELOG && H_logFile && H_logFile != stderr) {
		if (H_logFile == (FILE *)0xffff) {
			if (Sp_buf) {
				if (S_func)
					Sp_buf += sprintf(Sp_buf, "[%s] ", S_func);
				Sp_buf += sprintf(Sp_buf, "[%s] ", S_prefix);
				va_start(H_varList, S_prtFormat);
				Sp_buf += vsprintf(Sp_buf, S_prtFormat, H_varList);
			}// else {
			// 				printf("[%s] ", S_func);
			// 				va_start(H_varList, S_prtFormat);
			// 				vprintf(S_prtFormat, H_varList);
			// 			}
		} else {
			if (S_func)
				fprintf(H_logFile, "[%s] ", S_func);
			fprintf(H_logFile, "[%s] ", S_prefix);
			va_start(H_varList, S_prtFormat);
			vfprintf(H_logFile, S_prtFormat, H_varList);
		}
	}
	va_end(H_varList);
}

void __verbose(int B_logThis, const char *S_func, const char *S_prtFormat, ...)
{
	if (OPT_NUMBER(verbose) == 0)
		return;

	va_list	H_varList;
	FILE	*H_targetFile;

	LOGTYPE	X_logType = B_noLog ? L_NOLOG : (B_logThis ? L_FORCELOG : L_NOTICE);

	switch (X_logType) {
	case L_NOLOG:    H_targetFile = NULL;		break;
	case L_FORCELOG: H_targetFile = stderr;		break;
	case L_NOTICE:	 H_targetFile = stderr;		break;
	case L_LOG:		 H_targetFile = H_logFile;	break;
	case L_LOGONLY:	 H_targetFile = H_logFile;	break;
	default:
		notice("Invalid stream given! [%s]", S_prtFormat);
		exit(1);
	}
	if (H_targetFile == NULL)
		return;

	va_start(H_varList, S_prtFormat);
	if ((X_logType == L_LOGONLY && H_logFile != stdout) ||
		X_logType != L_LOGONLY) {
		if (S_func)
			fprintf(H_targetFile, "LOG [%s] ", S_func);
		vfprintf(H_targetFile, S_prtFormat, H_varList);
	}
	if (X_logType == L_FORCELOG && H_logFile && H_logFile != stderr) {
		if (H_logFile == (FILE *)0xffff) {
			if (Sp_buf) {
				if (S_func)
					Sp_buf += sprintf(Sp_buf, "LOG [%s] ", S_func);
// 				else
// 					Sp_buf += sprintf(Sp_buf, "LOG ");
				va_start(H_varList, S_prtFormat);
				Sp_buf += vsprintf(Sp_buf, S_prtFormat, H_varList);
			}// else {
// 				if (S_func)
// 					printf("LOG [%s] ", S_func);
// 				else
// 					printf("LOG ");
// 				va_start(H_varList, S_prtFormat);
// 				vprintf(S_prtFormat, H_varList);
// 			}
		} else {
			if (S_func)
				fprintf(H_logFile, "LOG [%s] ", S_func);
// 			else
// 				fprintf(H_logFile, "LOG ");
			va_start(H_varList, S_prtFormat);
			vfprintf(H_logFile, S_prtFormat, H_varList);
		}
	}
	va_end(H_varList);
}

/* 로깅 함수 */
int __log(LOGTYPE X_logType, const char *S_prtFormat, ...)
{
	va_list	H_varList;
	int		N_prtChrCount	= 0;
	FILE	*H_targetFile	= NULL;

	switch (X_logType) {
	case L_NOLOG:		H_targetFile = NULL;		break;
	case L_FORCELOG:	H_targetFile = stderr;		break;
	case L_NOTICE:		H_targetFile = stderr;		break;
	case L_LOG:			H_targetFile = H_logFile;	break;
	case L_LOGONLY:		H_targetFile = H_logFile;	break;
	default:
		halt_fmt(WISARD_SYST_INVL_LOGSTREAM, X_logType, S_prtFormat);
		break;
	}
	if (H_targetFile == NULL)
		return 0;

	va_start(H_varList, S_prtFormat);
	if ((X_logType == L_LOGONLY && H_logFile != stdout) ||
		X_logType != L_LOGONLY)
		N_prtChrCount = vfprintf(H_targetFile, S_prtFormat, H_varList);
	if (X_logType == L_FORCELOG && H_logFile && H_logFile != stderr) {
		if (H_logFile == (FILE *)0xffff) {
			va_start(H_varList, S_prtFormat);
			Sp_buf += vsprintf(Sp_buf, S_prtFormat, H_varList);
		} else {
			va_start(H_varList, S_prtFormat);
			vfprintf(H_logFile, S_prtFormat, H_varList);
		}
	}
	va_end(H_varList);

	return N_prtChrCount;
}

int logSet(LOGTYPE X_logType, const char *S_FN)
{
	switch (X_logType) {
		case L_NOLOG:
			H_logFile = NULL;
			B_noLog = 1;
			break;
		case L_LOG:
			if (!S_FN[0]) {
				H_logFile = (FILE *)0xffff;
				wsAlloc(S_buf, char, 1024*1024);
				S_buf[0] = '\0';
				Sp_buf = S_buf;
			} else {
				int B_save = H_logFile == (FILE *)0xffff;

				/* --version */
				if (OPT_ENABLED(version)) {
					S_buf[0] = '\0';
					B_save = 0;
				} else {
					H_logFile = fopen(S_FN, "w+");
					if (!H_logFile) {
						DEALLOC(S_buf);
						halt_fmt(WISARD_FAIL_INIT_LOG, S_FN);
						return 0;
					}
				}
				if (B_save) {
					fprintf(H_logFile, "%s", S_buf);
					DEALLOC(S_buf);
					S_buf = Sp_buf = NULL;
				}
				LOGoutput("Execution log is exported to [%s.log]\n",
					OPT_STRING(out));
			}
			break;
		default:
			fprintf(stderr, "Invalid stream assigned by '%s'", S_FN);
			exit(1);
	}

	return 1;
}

int logSet(LOGTYPE X_logType, FILE *H_target)
{
	switch (X_logType) {
		case L_LOG:
			H_logFile = H_target;
			break;
		default:
			fprintf(stderr, "Invalid stream assigned!!");
			exit(1);
	}

	return 1;
}

void logEnd()
{
	/* Export system report */
	reportEnd();
	if (H_logFile != stdout && H_logFile != (FILE *)0xffff &&
		H_logFile) fclose(H_logFile);
	H_logFile = (FILE *)0xffff;
}

#ifdef _M_ARM
sse_t neon_div(sse_t a, sse_t b)
{
	sse_t rec = vrecpeq_f32(b);
	rec = vmulq_f32(vrecpsq_f32(b, rec), rec);
	rec = vmulq_f32(vrecpsq_f32(b, rec), rec);
	return vmulq_f32(a , rec);
}
sse_t neon_and(sse_t a, sse_t b)
{
	uint32x4_t *va = (uint32x4_t *)&a;
	uint32x4_t *vb = (uint32x4_t *)&b;
	uint32x4_t vr = vandq_s32(*va, *vb);
	return *((sse_t *)&vr);
}
sse_t neon_xor(sse_t a, sse_t b)
{
	uint32x4_t *va = (uint32x4_t *)&a;
	uint32x4_t *vb = (uint32x4_t *)&b;
	uint32x4_t vr = veorq_s32(*va, *vb);
	return *((sse_t *)&vr);
}
int _mm_popcnt_u32(unsigned int v)
{
	v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
	v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
	return ((v + (v >> 4) & 0xF0F0F0F) * 0x1010101) >> 24; // count
}
#endif

xSysSpec G;

} // End namespace ONETOOL
