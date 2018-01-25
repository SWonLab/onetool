// 120119
// Sungyoung Lee, BIBS
//
// cOption.cpp
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits>
#include <string>
#include <map>
#include "global/common.h"
#include "global/option.h"
#include "global/interactive.h"
#include "utils/util.h"
#include "input/stream.h"
#if TOOLSET_TYPE == TOOLSET_ONETOOL
#	include "global/worker.h"
#endif
using namespace std;

#define OPTDEF(ocat, apply, allo, lname, lsym, ssym, sname, defv, otype, rng) \
	{ apply, allo, Z(lname), Z(lsym), Z(ssym), Z(sname), Z(defv), otype, rng, {0, }, }
#define FEATURE_ALL
#define FEATURE_FQL
#define FEATURE_FVT
#define FEATURE_WSD
#define FEATURE_QTS
#define FEATURE_HMN
#define FEATURE_GSC

namespace ONETOOL {

xOption Xa_defOption[] = {
	/*APPLY  ALLOC  LONGNAME            SRTN      DEFVAL        TYPE       {0, }, },*/

	/* OPTIONS FOR PROGRAM FLOW */
	OPTDEF(FEATURE_ALL, true,  false, "--verbose", 0, 0, "-vb", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--seed", 0, 0, "-S", "", OT_NUMBER, RT_NUM),
	OPTDEF(FEATURE_ALL, false, false, "--script", 0, 0, "-r", "0", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, true,  false, "--thread", 0, 0, "-t", "1", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_ALL, true,  false, "--out", 0, 0, "-o", "onetool", OT_STRING, RT_OTHER),

	/* OPTIONS FOR DATA INPUT */
	OPTDEF(FEATURE_ALL, false, false, "--lgen", 0, 0, "-l", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--tped", "--tfile", 0, "-tp", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--data", 0, 0, "-D", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--bed", "--bfile", 0, "-j", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--fam", 0, 0, "-F", "", OT_STRING, RT_OTHER),

	/* OPTIONS FOR ADDITIONAL INPUT */
	OPTDEF(FEATURE_ALL, false, false, "--bim", 0, 0, "-bi", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--map", 0, 0, "-M", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--pheno", "--sampvar", 0, "-h", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--sampvarflag", "", 0, "-ai", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--geneset", 0, 0, "-pw", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--set", 0, 0, "-s", "", OT_STRING, RT_OTHER),

	/* OPTIONS FOR HANDLING ADDITIONAL INPUT */
	OPTDEF(FEATURE_ALL, false, false, "--ignorefid", 0, 0, "-If", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--nofid", 0, 0, "-Nf", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--noparent", 0, 0, "-np", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--nosex", 0, 0, "-nX", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--nopos", 0, 0, "-no", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--nogdist", 0, 0, "-ng", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--acgt", 0, 0, "-ac", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--nomap", 0, 0, "-nm", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--1234", 0, 0, "-zz", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--sepallele", 0, 0, "-sa", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--consecallele", 0, 0, "-cl", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, true, false, "--1case", 0, 0, "-x1", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, true, false, "--mispheno", 0, 0, "-mh", "-9", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, true, false, "--misgeno", 0, 0, "-mg", "0", OT_STRING, RT_OTHER),

	/* ELEMENTARY ANALYSIS OPTIONS */
	OPTDEF(FEATURE_ALL, false, false, "--freq", 0, 0, "-f", "founder", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--hwe", 0, 0, "-H", "founder", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--heritability", 0, 0, "-he", "0", OT_ONOFF, RT_ONOFF),

	/* GENE-BASED TESTS */
	OPTDEF(FEATURE_ALL, false, false, "--noweight", "--equalweight", 0, "-ew", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--betaweight", 0, 0, "-bw", "1,25", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--weight", 0, 0, "-w", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--mafweight", 0, 0, "-Mw", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--gsmacthr", 0, 0, "-gt", "0", OT_NUMBER, RT_NONNEG),

	/* VARIANT-LEVEL TESTS */
	OPTDEF(FEATURE_ALL, false, false, "--prevalence", 0, 0, "-pv", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--scoretest", 0, 0, "-ct", "0", OT_ONOFF, RT_ONOFF),

	/* SAMPLE RELATEDNESS OPTIONS */
	OPTDEF(FEATURE_ALL, false, false, "--kinship", "--pddt", 0, "-k", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--ibs", 0, 0, "-bs", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--ktau", 0, 0, "-kt", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--empktau", 0, 0, "-ek", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--corpearson", 0, 0, "-cp", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--cordiag1", 0, 0, "-cd", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--medcor", 0, 0, "-Mc", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--indep", 0, 0, "-id", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--hybrid", 0, 0, "-hy", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--cor", 0, 0, "-c", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--corpair", 0, 0, "-ca", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--corgrm", 0, 0, "-cg", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--corepacts", 0, 0, "-ep", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--x", 0, 0, "-x", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--x2", 0, 0, "-x2", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--cormaf", 0, 0, "-C", "(0.05,1]", OT_RANGE, RT_01),

	/* EXPORT-RELATED OPTIONS */

	OPTDEF(FEATURE_ALL, false, false, "--quiet", 0, 0, "-q", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--version", 0, 0, "-v", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--ped", "--file", 0, "-i", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--variantvar", 0, 0, "-a", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--filvariant", 0, 0, "-fk", "", OT_EXPR, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--incvariant", 0, 0, "-ik", "", OT_EXPR, RT_OTHER),	
	OPTDEF(FEATURE_ALL, false, false, "--filgeno", 0, 0, "-fj", "", OT_EXPR, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--incgeno", 0, 0, "-ij", "", OT_EXPR, RT_OTHER),	
	OPTDEF(FEATURE_ALL, false, false, "--pname", 0, 0, "-n", "*", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--remsamp", "--indrem", 0, "-ir", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--selsamp", "--indsel", 0, "-is", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--nasamp", 0, 0, "-xs", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--randnasamp", 0, 0, "-rx", "", OT_REAL, RT_PROPNUM),
	OPTDEF(FEATURE_ALL, false, false, "--remvariant", "--snprem", 0, "-sr", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--selvariant", "--snpsel", 0, "-Ss", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--remfam", 0, 0, "-fF", "0", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--selfam", 0, 0, "-iF", "0", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--filrange", 0, 0, "-rr", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--incrange", 0, 0, "-rs", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--varresize", "--varesize", 0, "-sz", "0", OT_REAL, RT_PROPNUM),
	OPTDEF(FEATURE_ALL, false, false, "--varwindow", 0, 0, "-vw", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_ALL, false, false, "--sampresize", 0, 0, "-az", "0", OT_REAL, RT_PROPNUM),
	OPTDEF(FEATURE_ALL, false, false, "--emcount", 0, 0, "-ec", "5", OT_NUMBER, RT_NONNEG),
	OPTDEF(FEATURE_ALL, false, false, "--aithr", 0, 0, "-at", "1e-8", OT_REAL, RT_POS),
	OPTDEF(FEATURE_ALL, false, false, "--autoonly", 0, 0, "-U", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--sexonly", 0, 0, "-sx", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--chr", 0, 0, "-ch", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--makeblup", 0, 0, "-mu", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--filmaf", "--filfreq", 0, "-ff", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_ALL, false, false, "--filmac", 0, 0, "-fm", "", OT_RANGE, RT_NONNEG),
	OPTDEF(FEATURE_ALL, false, false, "--filhwe", 0, 0, "-fh", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_ALL, false, false, "--filgind", 0, 0, "-fg", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_ALL, false, false, "--filgvar", 0, 0, "-fp", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_ALL, false, false, "--filnf", 0, 0, "-fN", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--filmf", 0, 0, "-fs", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--filcase", 0, 0, "-fa", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--filcontrol", 0, 0, "-ft", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--filmispheno", 0, 0, "-fo", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--filsample", 0, 0, "-fn", "", OT_EXPR, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--filgdist", 0, 0, "-fi", "", OT_RANGE, RT_POS),
	OPTDEF(FEATURE_ALL, false, false, "--filnosex", 0, 0, "-fx", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--filmale", 0, 0, "-fl", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--filfemale", 0, 0, "-fe", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--snvonly", 0, 0, "-vo", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--indelonly", 0, 0, "-xo", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--incmaf", "--incfreq", 0, "-if", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_ALL, false, false, "--incmac", 0, 0, "-im", "", OT_RANGE, RT_NONNEG),
	OPTDEF(FEATURE_ALL, false, false, "--inchwe", 0, 0, "-ih", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_ALL, false, false, "--incgind", 0, 0, "-ig", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_ALL, false, false, "--incgvar", 0, 0, "-iv", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_ALL, false, false, "--incgdist", 0, 0, "-ii", "", OT_RANGE, RT_POS),
	OPTDEF(FEATURE_ALL, false, false, "--incsample", 0, 0, "-in", "", OT_EXPR, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--ml", 0, 0, "-ml", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--blup", 0, 0, "-b", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--imputepheno", 0, 0, "-ip", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--cname", 0, 0, "-cn", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--makecor", 0, 0, "-MC", "default", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--empiall", 0, 0, "-ea", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--cact", 0, 0, "-cc", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, true, false, "--1sex", 0, 0, "-s1", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--mafe", 0, 0, "-MF", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--ginv", 0, 0, "-gi", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--makecov", 0, 0, "-mc", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--makepheno", 0, 0, "-mP", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--founderonly", 0, 0, "-Fo", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--famsummary", 0, 0, "-Fs", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--ncsummary", 0, 0, "-nc", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--specdcmp", 0, 0, "-sd", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--nostop", 0, 0, "-xx", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--ignoreparent", 0, 0, "-Ip", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--sepid", 0, 0, "-Si", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--nopheno", 0, 0, "-nP", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--probandcol", 0, 0, "-pb", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--twincol", 0, 0, "-tw", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--chrwise", 0, 0, "-cw", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--est", 0, 0, "-E", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--makeclgeno", 0, 0, "-mz", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--indel", 0, 0, "-ID", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--logistic", 0, 0, "-lg", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--sortvariant", 0, 0, "-sv", "asc", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--sortpos", 0, 0, "-sp", "asc", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--sortsample", 0, 0, "-ss", "asc", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--sortiid", 0, 0, "-si", "asc", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--baseline", 0, 0, "-bl", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--gz", 0, 0, "-Z", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--makenrm", 0, 0, "-my", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--makeweight", 0, 0, "-mw", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--nophenohdr", "--nosampvarhdr", 0, "-nh", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--pvalrange", 0, 0, "-pr", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_ALL, false, false, "--time", 0, 0, "-ti", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--phenostdize", 0, 0, "-pz", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--outmispheno", 0, 0, "-op", "-9", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--outmisgeno", 0, 0, "-og", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--gxecovs", 0, 0, "-gc", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--passemptyline", 0, 0, "-pe", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--listvariant", 0, 0, "-lm", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--listsample", 0, 0, "-ls", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--listfounder", 0, 0, "-lf", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--nolmm", 0, 0, "-nl", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--gsetconsec", 0, 0, "-GC", "", OT_RANGE, RT_POS),
	OPTDEF(FEATURE_ALL, false, false, "--nosysfeat", 0, 0, "-ns", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--fname", 0, 0, "-fc", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--makeev", 0, 0, "-mV", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--ev", 0, 0, "-ev", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--natural", 0, 0, "-nt", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--dupnaming", 0, 0, "-dp", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--citation", "--cite", 0, "-ce", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--miss", 0, 0, "-mx", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--misparent", 0, 0, "-mq", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--sampmajor", 0, 0, "-sm", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--nospecdcmp", 0, 0, "-ne", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--species", 0, 0, "-sc", "human", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--nodata", 0, 0, "-nd", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--regex", 0, 0, "-er", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, true,  false, "--model", 0, 0, "-J", "additive", OT_STRING, RT_OTHER),

#if (TOOLSET_TYPE == TOOLSET_MFQLS) || ((TOOLSET_TYPE & 0x100) == 0x100)
	/* VARIANT-LEVEL TESTS */
	OPTDEF(FEATURE_FQL, false, false, "--mqls", 0, 0, "-qm", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FQL, false, false, "--fqls", 0, 0, "-qf", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FQL, false, false, "--heri", 0, 0, "-hi", "0", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_FQL, false, false, "--multifqls", "--mfarvat", 0, "-qx", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FQL, false, false, "--mqlsconsec", 0, 0, "-qp", "", OT_NUMBER, RT_POS),

	OPTDEF(FEATURE_FQL, false, false, "--fqlsnopddt", 0, 0, "-pf", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FQL, false, false, "--retestthr", 0, 0, "-rt", "0.05", OT_REAL, RT_01),
	OPTDEF(FEATURE_FQL, false, false, "--avail", 0, 0, "-av", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FQL, false, false, "--fastmqls", 0, 0, "-Fa", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FQL, false, false, "--fastfqls", 0, 0, "-FF", "0", OT_ONOFF, RT_ONOFF),
#endif

#if (TOOLSET_TYPE == TOOLSET_FARVAT) || (TOOLSET_TYPE == TOOLSET_MFARVAT) || \
	((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE == TOOLSET_FARVATX)
	/* GENE-BASED TESTS */
	OPTDEF(FEATURE_FVT, false, false, "--genesummary", 0, 0, "-sg", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--genetest", 0, 0, "-tg", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--raremaf", 0, 0, "-rm", "0.01", OT_REAL, RT_01),
	OPTDEF(FEATURE_FVT, false, false, "--genesize", "--gsrange", 0, "-gr", "", OT_RANGE, RT_NONNEG),
	OPTDEF(FEATURE_FVT, false, false, "--pedcmc", 0, 0, "-pc", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--wsum", 0, 0, "-ws", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--kbac", 0, 0, "-K", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--asum", 0, 0, "-au", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--farvat", 0, 0, "-Fv", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--pedgene", 0, 0, "-pg", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--skato", 0, 0, "-so", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--skat", 0, 0, "-sk", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--mfhom", 0, 0, "-mo", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--mfhet", 0, 0, "-mf", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--adjf1", 0, 0, "-f1", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--adjf2", 0, 0, "-f2", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--farvatx", 0, 0, "-Fx", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--farvatxd", 0, 0, "-Fd", "0.5", OT_REAL, RT_01),

	OPTDEF(FEATURE_FVT, false, false, "--gmapsummary", 0, 0, "-xg", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--genemiss", "--gsmissthr", 0, "-gm", "0.05", OT_REAL, RT_01),
	OPTDEF(FEATURE_FVT, false, false, "--skatondiv", 0, 0, "-ov", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_FVT, false, false, "--skatodivs", 0, 0, "-oz", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_FVT, false, false, "--kbacalpha", 0, 0, "-ka", "0.05", OT_REAL, RT_01),
	OPTDEF(FEATURE_FVT, false, false, "--kbac2side", 0, 0, "-k2", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--kbackernel", 0, 0, "-kk", "hypergeometric", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_FVT, false, false, "--makegeno", 0, 0, "-mZ", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--makefarvat", 0, 0, "-kf", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_FVT, false, false, "--farvatxndiv", 0, 0, "-Xv", "5", OT_NUMBER, RT_POS),
#endif

#if (TOOLSET_TYPE == TOOLSET_QTEST) || ((TOOLSET_TYPE & 0x100) == 0x100)
	OPTDEF(FEATURE_QTS, false, false, "--qtest", 0, 0, "-Q", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_QTS, false, false, "--qtestrange", 0, 0, "-Qr", "(0,0.05)", OT_RANGE, RT_01),
	OPTDEF(FEATURE_QTS, false, false, "--qtestclump", 0, 0, "-Qc", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_QTS, false, false, "--qteststt", 0, 0, "-Qs", "0.2", OT_REAL, RT_01),
	OPTDEF(FEATURE_QTS, false, false, "--qtestbetacov", 0, 0, "-Qb", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_QTS, false, false, "--makebeta", 0, 0, "-aa", "0", OT_ONOFF, RT_ONOFF),
#endif

#if (TOOLSET_TYPE == TOOLSET_HIMINI) || ((TOOLSET_TYPE & 0x100) == 0x100)
	OPTDEF(FEATURE_HMN, false, false, "--mdr", 0, 0, "-md", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_HMN, false, false, "--order", 0, 0, "-O", "1", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_HMN, false, false, "--top", 0, 0, "-T", "1000", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_HMN, false, false, "--hmdr", 0, 0, "-hm", "", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_HMN, false, false, "--hmdrall", 0, 0, "-ha", "1", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_HMN, false, false, "--hmdrprior", 0, 0, "-hp", "1", OT_STRING, RT_OTHER),
#endif

#if (TOOLSET_TYPE == TOOLSET_PHARAOH) || ((TOOLSET_TYPE & 0x100) == 0x100)
	OPTDEF(FEATURE_GSC, false, false, "--pharaoh", 0, 0, "-gs", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_GSC, false, false, "--proopt", 0, 0, "-ho", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_GSC, false, false, "--prolambda", 0, 0, "-Pl", "1000", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_GSC, false, false, "--prorange", 0, 0, "-Pr", "", OT_RANGE, RT_POS),
	OPTDEF(FEATURE_GSC, false, false, "--prothr", 0, 0, "-Pt", "0.0001", OT_REAL, RT_POS),
	OPTDEF(FEATURE_GSC, false, false, "--promaxiter", 0, 0, "-gq", "100", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_GSC, false, false, "--progenesize", 0, 0, "-Pg", "", OT_REAL, RT_POS),
	OPTDEF(FEATURE_GSC, false, false, "--progsetsize", 0, 0, "-pt", "", OT_REAL, RT_POS),
	OPTDEF(FEATURE_GSC, false, false, "--prosingle", 0, 0, "-Pe", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_GSC, false, false, "--propermcov", 0, 0, "-Pp", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--nperm", 0, 0, "-Ne", "1000", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--seqperm", 0, 0, "-SP", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--permfile", 0, 0, "-PF", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--cv", 0, 0, "-cv", "", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_GSC, false, false, "--gesca", 0, 0, "-ge", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--modeltype", 0, 0, "-MT", "formative", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--ggpath", 0, 0, "-gp", "", OT_STRING, RT_OTHER),
#endif

#if TOOLSET_TYPE == TOOLSET_ONETOOL
	OPTDEF(FEATURE_GSC, false, false, "--debug", 0, 0, "-db", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_GSC, false, false, "--help", 0, 0, "-hl", "0", OT_ONOFF, RT_ONOFF),
#endif

#if ((TOOLSET_TYPE & 0x100) == 0x100)
	/* OPTIONS FOR DATA INPUT */
	OPTDEF(FEATURE_WSD, false, false, "--dosage", 0, 0, "-d", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--expression", 0, 0, "-EX", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--genoprob", 0, 0, "-go", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--vcf", 0, 0, "-V", "", OT_STRING, RT_OTHER),

	/* ELEMENTARY ANALYSIS OPTIONS */
	OPTDEF(FEATURE_WSD, false, false, "--fst", 0, 0, "-Ft", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--mendel", 0, 0, "-m", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--famuniq", 0, 0, "-fU", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--tstv", 0, 0, "-tt", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--pca", 0, 0, "-P", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--npc", 0, 0, "-N", "5", OT_NUMBER, RT_POS),

	/* GENE-BASED TESTS */
	OPTDEF(FEATURE_WSD, false, false, "--genesplit", 0, 0, "-GS", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--rvtdt", 0, 0, "-rd", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--fbskat", 0, 0, "-fK", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--famvt", 0, 0, "-fv", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--ggemma", 0, 0, "-zg", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--longitudinal", 0, 0, "-L", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--vt", 0, 0, "-vt", "0", OT_ONOFF, RT_ONOFF),

	/* VARIANT-LEVEL TESTS */
	OPTDEF(FEATURE_WSD, false, false, "--fisher", 0, 0, "-fy", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--trend", 0, 0, "-CA", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--regression", 0, 0, "-rg", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--qls", 0, 0, "-ql", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--tdt", 0, 0, "-td", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--sdt", 0, 0, "-st", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--emmax", 0, 0, "-ex", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--gemma", "--femma", 0, "-G", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--powercalc", 0, 0, "-pl", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--powercalc2", 0, 0, "-px", "0", OT_ONOFF, RT_ONOFF),

	/* SAMPLE RELATEDNESS OPTIONS */

	/* EXPORT-RELATED OPTIONS */
	OPTDEF(FEATURE_WSD, false, false, "--makeped", 0, 0, "-mp", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--maketped", 0, 0, "-mT", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makebed", "--make-bed", 0, "-mb", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makeraw", "--recodeA", 0, "-mr", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makedom", 0, 0, "-mD", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makerec", 0, 0, "-mR", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makevcf", 0, 0, "-mv", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makebcf", 0, 0, "-mF", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makelgen", 0, 0, "-mn", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makegen", 0, 0, "-mG", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makebgen", 0, 0, "-mB", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makebeagle", 0, 0, "-mL", "", OT_ONOFF, RT_ONOFF),

	/* S.A.G.E. OPTIONS */
	OPTDEF(FEATURE_ALL, false, false, "--relpair", 0, 0, "-srel", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--fcor", 0, 0, "-sfco", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--segreg", 0, 0, "-sseg", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--lodlink", 0, 0, "-slod", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--typ", "--mfile", 0, "-sty", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--par", "--parfile", 0, "-spf", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_ALL, false, false, "--fcorStdErrOff", 0, 0, "-fseo", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--lodlinkLinkageTestOff", 0, 0, "-llto", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--lodlinkLinkageHomogOff", 0, 0, "-lho", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--lodlinkLinkageSexSpecific", 0, 0, "-lss", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--lodlinkSmithHomogTest", 0, 0, "-lsht", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--lodlinkGenotypes", 0, 0, "-lgt", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--merlin", 0, 0, "-mer", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_ALL, false, false, "--plot", 0, 0, "-rP", "0", OT_ONOFF, RT_ONOFF),

	OPTDEF(FEATURE_WSD, false, false, "--simtrio", 0, 0, "-xt", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--szfam", 0, 0, "-sf", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--szvar", 0, 0, "-z", "0", OT_NUMBER, RT_NONNEG),
	OPTDEF(FEATURE_WSD, false, false, "--simfam", 0, 0, "-sF", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--trio", 0, 0, "-sT", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--extfam", 0, 0, "-sE", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--nsig", 0, 0, "-nx", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--sigmaf", 0, 0, "-XM", "0.1", OT_REAL, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--simfreq", 0, 0, "-sq", "0", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--mafvar", 0, 0, "-ms", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--proppc", 0, 0, "-pp", "", OT_REAL, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--usemf", "--usepp", 0, "-u", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--fullpca", 0, 0, "-fu", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--rpath", 0, 0, "-rl", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--grm", 0, 0, "-gR", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--grmalpha", 0, 0, "-ga", "0.001", OT_REAL, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--makemdr", 0, 0, "-mM", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--zipbgen", 0, 0, "-zb", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--filqual", 0, 0, "-fq", "", OT_RANGE, RT_NONNEG),
	OPTDEF(FEATURE_WSD, false, false, "--incqual", 0, 0, "-iq", "", OT_RANGE, RT_NONNEG),
	OPTDEF(FEATURE_WSD, false, false, "--filmendelfam", 0, 0, "-rF", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--incmendelfam", 0, 0, "-xF", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--filmendelsamp", 0, 0, "-rS", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--incmendelsamp", 0, 0, "-xS", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--filmendelvar", 0, 0, "-ds", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--incmendelvar", 0, 0, "-xm", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--lrt", 0, 0, "-lr", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--vcfqc", 0, 0, "-vq", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--phasedonly", 0, 0, "-po", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--unphasedonly", 0, 0, "-uo", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--interactive", 0, 0, "-I", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--check", 0, 0, "-ck", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--gxe", 0, 0, "-gx", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--beta", 0, 0, "-be", "0", OT_REAL, RT_NUM),
	OPTDEF(FEATURE_WSD, false, false, "--rho", 0, 0, "-rh", "0", OT_REAL, RT_ABS1),
	OPTDEF(FEATURE_WSD, false, false, "--rhopheno", 0, 0, "-rp", "0", OT_REAL, RT_ABS1),
	OPTDEF(FEATURE_WSD, false, false, "--nsamp", 0, 0, "-Ns", "0", OT_RANGE, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--npheno", 0, 0, "-Np", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--noshuffle", 0, 0, "-Nh", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--shuffle", 0, 0, "-sh", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--split", 0, 0, "-sl", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--merge", 0, 0, "-me", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--mergemode", 0, 0, "-mm", "1", OT_NUMBER, RT_NONNEG), /**/
	OPTDEF(FEATURE_WSD, false, false, "--testmatrix", 0, 0, "-Xm", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--testmatfunc", 0, 0, "-Xt", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--testmatclass", 0, 0, "-Xc", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--testfunc", 0, 0, "-Xf", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--ld", 0, 0, "-ld", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--ldcor", 0, 0, "-lc", "0", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--annogene", 0, 0, "-ag", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--annorange", 0, 0, "-ar", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--annovar", 0, 0, "-as", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--density", 0, 0, "-de", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--sep", 0, 0, "-xp", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--dsgdist", 0, 0, "-dd", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--pc2cov", 0, 0, "-p2", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--donull", 0, 0, "-dn", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--explore", 0, 0, "-e", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--updvariant", 0, 0, "-um", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--updchr", 0, 0, "-uc", "*", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--updname", 0, 0, "-un", "*", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--updgdist", 0, 0, "-ud", "*", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--updpos", 0, 0, "-up", "*", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--updgeno", 0, 0, "-ug", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--ref", 0, 0, "-rf", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--ldbin", 0, 0, "-lb", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--ldsize", 0, 0, "-lS", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--ldvar", 0, 0, "-lM", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--remna", 0, 0, "-rn", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--boost", 0, 0, "-bo", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--thrboost", 0, 0, "-Bt", "30", OT_REAL, RT_NONNEG),
	OPTDEF(FEATURE_WSD, false, false, "--quickepi", 0, 0, "-qe", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--ext", 0, 0, "-X", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--gmdr", 0, 0, "-gd", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--impute", 0, 0, "-ie", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--lod", 0, 0, "-lo", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makeimpute", 0, 0, "-mi", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--sxa", 0, 0, "-xa", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--R", 0, 0, "-R", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--randbinpheno", 0, 0, "-rb", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--randpheno", 0, 0, "-re", "1", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--genoctrl", 0, 0, "-gl", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--usergc", 0, 0, "-eg", "", OT_REAL, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--adjust", 0, 0, "-A", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--popuniq", 0, 0, "-pU", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--monotone", 0, 0, "-tm", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--singleton", 0, 0, "-ts", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--doubleton", 0, 0, "-tb", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--genemdr", 0, 0, "-Gm", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--variant2cov", 0, 0, "-m2", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--inbreed", 0, 0, "-ib", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--group", 0, 0, "-g", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--filgenic", 0, 0, "-fG", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--filintergenic", 0, 0, "-fI", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--variantsummary", 0, 0, "-ks", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--sampleorder", 0, 0, "-sO", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--variantorder", 0, 0, "-sM", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--cosi", 0, 0, "-co", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--setconsec", 0, 0, "-xc", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--setoverlap", 0, 0, "-xv", "0", OT_NUMBER, RT_NONNEG),
	OPTDEF(FEATURE_WSD, false, false, "--setrandom", 0, 0, "-xr", "", OT_RANGE, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--makeset", 0, 0, "-mS", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--settype", 0, 0, "-et", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--outcact", 0, 0, "-oc", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--out1case", 0, 0, "-o1", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--mistest", 0, 0, "-mt", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--incmistest", 0, 0, "-iM", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--filmistest", 0, 0, "-fM", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--nageno", 0, 0, "-Ng", "", OT_REAL, RT_NONNEG),
	OPTDEF(FEATURE_WSD, false, false, "--napheno", 0, 0, "-NH", "", OT_REAL, RT_NONNEG),
	OPTDEF(FEATURE_WSD, false, false, "--filtreport", 0, 0, "-fr", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--genofield", 0, 0, "-gf", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--outphenoonly", 0, 0, "-oo", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--outnoheader", 0, 0, "-on", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--makemerlin", 0, 0, "-om", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--mds", 0, 0, "-MD", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--gxg", 0, 0, "-gg", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--invnorm", 0, 0, "-vn", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--forceconv", 0, 0, "-fV", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--out1234", 0, 0, "-t1", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--outacgt", 0, 0, "-oa", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--dfam", 0, 0, "-da", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--updallele", 0, 0, "-ua", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--window", 0, 0, "-W", "0", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--ci", 0, 0, "-ci", "0.05", OT_REAL, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--lasso", 0, 0, "-la", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--lassolambda", 0, 0, "-ll", "0.1", OT_REAL, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--lassoall", 0, 0, "-lL", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--pls", 0, 0, "-ps", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--sampleweight", 0, 0, "-sw", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--nskip", 0, 0, "-nk", "", OT_NUMBER, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--singleparent", 0, 0, "-sP", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--bcf", 0, 0, "-B", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--loocv", 0, 0, "-lv", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--mdrthr", 0, 0, "-dt", "", OT_RANGE, RT_01),
	OPTDEF(FEATURE_WSD, false, false, "--ldcontrast", 0, 0, "-lx", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--flip", 0, 0, "-fL", "TGCA", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--varsubset", 0, 0, "-fT", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--hethom", 0, 0, "-hh", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--markercheck", 0, 0, "-mk", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--meta", 0, 0, "-p", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--fid", 0, 0, "-xf", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--outformat", 0, 0, "-of", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--maf", 0, 0, "-ma", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--het", 0, 0, "-ht", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--variantblup", 0, 0, "-vp", "", OT_REAL, RT_POS),
	OPTDEF(FEATURE_WSD, false, false, "--famsplit", 0, 0, "-FS", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--setspan", 0, 0, "-sn", "0", OT_NUMBER, RT_NONNEG),

	OPTDEF(FEATURE_WSD, false, false, "--tridge", 0, 0, "-tr", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--hamming", 0, 0, "-hn", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--baldingnichols", "--bn", 0, "-bn", "", OT_ONOFF, RT_ONOFF),

	OPTDEF(FEATURE_WSD, false, false, "--fuzzymdr", 0, 0, "-zm", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--gxgall", 0, 0, "-gX", "", OT_ONOFF, RT_ONOFF),
	OPTDEF(FEATURE_WSD, false, false, "--gxglist", 0, 0, "-xl", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--gxglambda", 0, 0, "-xb", "", OT_STRING, RT_OTHER),

	OPTDEF(FEATURE_WSD, false, false, "--prunevif", 0, 0, "-rv", "", OT_STRING, RT_OTHER),
	OPTDEF(FEATURE_WSD, false, false, "--prunepw", 0, 0, "-rw", "", OT_STRING, RT_OTHER),
#endif

	/*** STEP 1 : Add the specification of new option here */
};

const char* Sa_descOption[] = {
	/* OPTIONS FOR PROGRAM FLOW */
	"Provide verbose output",
	"Random seed for program execution",
	"Run WISARD with the parameters specified in given path",
	"Number of threads used in analysis",
	"Output prefix",

	/* OPTIONS FOR DATA INPUT */
	"Long file format name",
	"TPED file name",
	"Shared data path",
	"PLINK Binary PED file name",
	"Alternative family definition file",

	/* OPTIONS FOR ADDITIONAL INPUT */
	"Assign alternative BIM file",
	"Use alternative MAP file when PED is given",
	"Alternative phenotype file",
	"Sample variable reading flags",
	"Path of gene-set definition",
	"List of gene-variant pairs that represents relationship",

	/* OPTIONS FOR HANDLING ADDITIONAL INPUT */
	"Ignores FID fields from input file",
	"Assume there is no FID definition in the pedigree info",
	"Assume there is no parent definition in the pedigree info",
	"Assume there is no sex definition in the pedigree info",
	"Assume there is no position definition in the variants info",
	"Assume there is no genetic distance definition in the variants info",
	"Redefine character for A,C,G,T",
	"Do not use map file when PED file is input",
	"Assume 1,2,3,4 are corresponding to A,C,G,T respectively",
	"Define a separator between two alleles",
	"Assume two alleles of a genotype is consecutive",
	"Treat 1=case/0=control for binary phenotype",
	"String for phenotype missingness",
	"String for genotype missingness",

	/* ELEMENTARY ANALYSIS OPTIONS */
	"Report frequency summary",
	"Report Hardy-Weinberg Equilibrium test summary",
	"Calculate heritability-related estimations",

	/* GENE-BASED TESTS */
	"Giving equal weight across all variants in gene-level test",
	"Giving beta-distribution based weight for each variant",
	"User-defined variant-wise weight file",
	"Giving inverse MAF-based weight for each variant",
	"Filtering out gene sets by its pooled minor allele count",

	/* VARIANT-LEVEL TESTS */
	"Prevalence of disease",
	"Use maximum-likelihood method for score test",

	/* SAMPLE RELATEDNESS OPTIONS */
	"Use 2*kinship coefficient instead of correlation",
	"Use IBS instead of correlation",
	"Set the sample correlation structure as Kendall's tau",
	"Set the sample correlation structure as Kendall's tau from empirical genotype data",
	"Calculate Pearson's correlation matrix",
	"Force the diagonals correlation matrix to 1",
	"Using median instead of mean in case of correlation",
	"Assume independency on samples",
	"Using hybrid correlation structure",
	"Assign correlation matrix",
	"Export sample relatedness with paired form instead of matrix form",
	"Export sample relatedness with GCTA form instead of matrix form",
	"Export sample relatedness matrix with EPACTS kin format",
	"Assume X chromosome inheritance on kinship computation",
	"Assume X chromosome inheritance on kinship computation with second method",
	"Range of MAF for relatedness estimation",

	/* EXPORT-RELATED OPTIONS */

	"Enabling quiet mode",
	"Print program version and stop",
	"PLINK PED file name",
	"Variant variable file", 
	"Filtering variants satisfying given condition(s)",
	"Including variants satisfying given condition(s)",
	"Filtering genotypes satisfying given condition(s)",
	"Including genotypes satisfying given condition(s)",
	"Phenotype(s) to be included",
	"List of samples to remove from analysis",
	"List of samples to select from analysis",
	"List of samples to nullify its genotype",
	"Number of samples that randomly nullified its genotype",
	"List of variants to remove from analysis",
	"List of variants to select from analysis",
	"List of families to remove from analysis",
	"List of families to select from analysis",
	"Range of variants to remove from analysis",
	"Range of variants to select from analysis",
	"Number of variants to be randomly resized",
	"Size of variant-thinning window",
	"Number of samples to be randomly resized",
	"Iteration count of EM",
	"AI convergence threshold",
	"Use only autosomes in analysis",
	"Use only sex chromosomes in analysis",
	"Select chromosomes to use in analysis",
	"Export BLUP value from AI/NR analysis",
	"Filtering out variants by minor allele frequency",
	"Filtering out variants by minor allele count",
	"Filtering out variants by Hardy-Weinberg Equilibrium test",
	"Filtering out samples by genotyping rate",
	"Filtering out variants by genotyping rate",
	"Filtering out non-founder samples",
	"Filtering out missing founder samples",
	"Filtering out case samples",
	"Filtering out control samples",
	"Filtering out missing phenotype samples",
	"Filtering out samples that does not satisfies given condition(s) for sample variables",
	"Filtering out variants by genetic distance",
	"Filtering out samples having no sex",
	"Filtering out male samples",
	"Filtering out female samples",
	"Filtering out all variants which are not SNV",
	"Filtering out all variants which are not indel",
	"Including variants by minor allele frequency",
	"Including variants by minor allele count",
	"Including variants by Hardy-Weinberg Equilibrium test",
	"Including samples by genotyping rate",
	"Including variants by genotyping rate",
	"Including variants by genetic distance",
	"Including samples that does not satisfies given condition(s) for sample variables",
	"Perform score test",
	"Adjust Y using blup",
	"Impute phenotype value if missing",
	"Covariates to use in analysis",
	"Export loaded/computed sample relatedness",
	"Compute empirical kinship based on all samples",
	"Assign values representing case and control",
	"Treat 1=female/0=male for sex encoding",
	"Assign values representing male and female",
	"Force to use generalized inverse matrix to get phiInv",
	"Make covariate file from final dataset and given covariates",
	"Make phenotype file from final dataset and given phenotypes",
	"Use founder only in the analysis",
	"Export family structure summary",
	"Export nuclear family structure summary",
	"Use spectral decomposition instead of EM-AI algorithm",
	"Do not stop iteration until converge",
	"Ignores parent fields from input file",
	"Assume there is no phenotype definition in the pedgree info",
	"Assume there is a separator between FID and IID rather than whitespace",
	"Column name for proband information in --pheno",
	"Column name for twin information in --pheno",
	"Perform analysis chromosome-wisely",
	"Alternative estimation of sig2/sig2g as exported file",
	"Export clumped genotype",
	"Allows loading indel of the dataset",
	"Perform gene-level test with logistic regression fitting",
	"Sort the variants in final dataset before analysis",
	"Sort the variants in final dataset by physical position before analysis",
	"Sort the samples in final dataset before analysis, by FID and IID",
	"Sort the samples in final dataset before analysis, by only IID",
	"Level of factors to be a baseline for that factor",
	"Deflate all outputs by gzip",
	"Export normalized genotype matrix",
	"Export computed/imported variant weights",
	"Assume there is no header in file from --pheno",
	"Apply p-value inclusion range to the all results",
	"Reporting computation time for each test to the result",
	"Standardize phenotype by dividing original phenotype value by its s.d.",
	"Missing string for NA when make output",
	"Covariate names for GxE interaction",
	"Pass empty line if exists",
	"Write out the list of remained variants in the final dataset",
	"Write out the list of remained samples in the final dataset",
	"Write out the list of remained founders in the final dataset",
	"Do not apply the result of LMM to the analyses ",
	"Automatically generates gene-set definition by consecutive manner",
	"Disable all system features",
	"Column names to be recognized as factor covariates",
	"Export eigendecomposition result",
	"Import eigendecomposition result",
	"Apply natural sort on sorting",
	"Automatically renaming duplicated sample record",
	"Print out citation information",
	"Summarize genotyping missing",
	"Set the missing indicator for parent",
	"Export BED file with sample-major manner",
	"Do not use spectral decomposition in any case",
	"Set the species to be analyzed",
	"Prevent program to stop from the lack of analyzable sample",
	"Use regular expression in the acceptable options",
	"Genetic model used in the dataset",

#if (TOOLSET_TYPE == TOOLSET_MFQLS) || ((TOOLSET_TYPE & 0x100) == 0x100)
	/* VARIANT-LEVEL TESTS */
	"Perform Modified QLS test",
	"Perform family QLS test",
	"Heritability value that given to prior",
	"Perform multiple FQLS test",
	"A pairing approach will be applied to MQLS procedure",

	"Do not force kinship coefficient to calculate offset in family QLS test",
	"Re-testing p-value threshold for further analysis",
	"Use available genotypes only",
	"Perform fast version of MQLS with genotype imputation",
	"Perform fast version of FQLS with genotype imputation",
#endif

#if (TOOLSET_TYPE == TOOLSET_FARVAT) || (TOOLSET_TYPE == TOOLSET_MFARVAT) || \
	((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE == TOOLSET_FARVATX)
	/* GENE-BASED TESTS */
	"Print summary of given gene definition",
	"Perform Gene-level test",
	"Threshold of MAF determines 'rare' variant",
	"Range of # variant in a gene to include analysis",
	"Perform family-based CMC test on gene-level analysis",
	"Perform weight sum test",
	"Perform KBAC test",
	"Perform adaptive-sum test",
	"Perform FARVAT forcely",
	"Perform pedgene analysis",
	"Perform SKAT-optimal test on --genetest",
	"Perform SKAT test on --genetest",
	"Perform FARVAT-multi with homogeneous case",
	"Perform FARVAT-multi with heterogeneous case",
	"Adjust phenotype with FQLS method 1",
	"Adjust phenotype with FQLS method 2",
	"Perform FARVATX analysis",
	"Assumed D in FARVATX analysis",

	"Print summary of gene-variants mapping result under single binary phenotype",
	"Acceptable genotyping missing rate in gene-level test",
	"Number of divisions for optimal test",
	"User-defined weights for optimal test",
	"Set the significance level of KBAC test",
	"Perform two-side KBAC test",
	"Kernel to do KBAC test",
	"Make per-gene genotype output",
	"Generate extra output of FARVAT for meta-analysis",
	"Number of divisions for FARVATX test",
#endif

#if (TOOLSET_TYPE == TOOLSET_QTEST) || ((TOOLSET_TYPE & 0x100) == 0x100)
	"Perform Qtest",
	"Range of inclusion for Qtest",
	"Perform genotype clumping by r^2 in Qtest",
	"Set STT in Qtest",
	"Export beta's variance-covariance matrix",
	"Export estimated beta(s)",
#endif

#if (TOOLSET_TYPE == TOOLSET_HIMINI) || ((TOOLSET_TYPE & 0x100) == 0x100)
	"Perform MDR analysis",
	"The order of combination to investigate",
	"Report top N combinations only when performing MDR analysis",	
	"Perform Hierarchical MDR analysis",
	"Reports all possible HMDR results",
	"Give priors to HMDR analysis",
#endif

#if (TOOLSET_TYPE == TOOLSET_PHARAOH) || ((TOOLSET_TYPE & 0x100) == 0x100)
	"Perform PHARAOH analysis",
	"Perform the detection of optimal lambda for PHARAOH only",
	"Set the fixed lambda penalty for PHARAOH",
	"Set the search range of optimal lambda",
	"Set the convergence threshold of finding optimal lambda",
	"Set the number of maximum iterations in PHARAOH",
	"Set the range of # variants in the gene to be included in the analysis",
	"Set the range of # genes in the pathway to be included in the analysis",
	"Use a single penalty on PHARAOH analysis",
	"Permute covariates in PHARAOH",
	"Number of permutations to perform",
	"Sequence of permutations to perform",
	"User-defined permutations to perform, by sample IIDs",
	"Number of cross-validation to perform",
	"Perform GeSCA",
	"Model type of GeSCA to perform (formative/reflective)",
	"Gene-gene path data for pathway",
#endif

#if TOOLSET_TYPE == TOOLSET_ONETOOL
	"Run on debug mode",
	"Print help",
#endif

#if ((TOOLSET_TYPE & 0x100) == 0x100)
	/* OPTIONS FOR DATA INPUT */
	"Dosage input file prefix",
	"Expression input file prefix",
	"VCF input file name",

	/* ELEMENTARY ANALYSIS OPTIONS */
	"Compute Fst using population information",
	"Computing Mendelian error",
	"Report family-specific variants",
	"Calculate window-based transition/transversion rate",
	"Adjust phenotype using PCA",
	"Number of PCs selected in PCA",

	/* GENE-BASED TESTS */
	"Split dataset by gene",
	"Perform rare-variant TDT analysis",
	"Perform FB-SKAT analysis",
	"Perform variable-threshold test for family in gene-level test",
	"Perform gene-level GEMMA test",
	"Apply longitudinal analysis with given variance-covariance matrix",
	"Perform variable-threshold test in gene-level test",

	/* VARIANT-LEVEL TESTS */
	"Perform variant-wise Fisher's exact test",
	"Perform variant-wise Cochran-Armitage trend test",
	"Perform regression analysis",
	"Perform QLS test for each variant",
	"Perform TDT",
	"Perform SDT",
	"Perform EMMAX",
	"Perform GEMMA",
	"Power calculation",
	"Power calculation version 2",

	/* SAMPLE RELATEDNESS OPTIONS */

	/* EXPORT-RELATED OPTIONS */
	"Make PED file from final dataset",
	"Make TPED file from final dataset",
	"Make BED file from final dataset",
	"Make RAW file from final dataset",
	"Make dominant-coded RAW file from final dataset",
	"Make recessive-coded RAW file from final dataset",
	"Make VCF file from final dataset",
	"Make Binary VCF file from final dataset",
	"Make LGEN file from final dataset",
	"Make GEN file from final dataset",
	"Make Binary GEN file from final dataset",
	"Make Beagle input file from final dataset",

	/* S.A.G.E. OPTIONS */
	"Generate xml file with relative pairs",
	"Run S.A.G.E. fcor",
	"Run S.A.G.E. segreg",
	"Run S.A.G.E. lodlink",
	"Type probability file from segreg",
	"Parameter file for segreg",
	"Specify the option not to compute the standard errors for correlations in fcor",
	"Specify the option not to perform linkage tests for lodscores in lodlink",
	"Specify the option not to assume linkage homogeneity in lodlink",
	"Specify the option to use sex-specific recombination fraction for linkage test in lodlink",
	"Specify the option to perform Smith's test for linkage homogeneity in lodlink",
	"Specify the option to calculate genotype probabilities in lodlink",
	"Perform MERLIN analysis for OneTool",
	"Plot pedigree(s)",

	"Use trio-simulation data instead of input",
	"Number of families in simulation",
	"Number of variants in simulation",
	"Use family-simulation using given pedigree",
	"Generate trio in family simulation",
	"Generate extended family in family simulation",
	"Set the number of statistically significant variants in simulation",
	"Set MAF of simulation",
	"Use frequency files",
	"Assign specific MAF to variants simulated",
	"Select PCs in PCA, with variance explained proportion",
	"Use information of missing founders",
	"Perform full eigenreduction in PCA",
	"R executable path",
	"GRM matrix from GCTA to compare with kinship coefficient",
	"Alpha level for GRM matrix comparison",
	"Make MDR input file from final dataset",
	"Compress variant block probability data while making Binary GEN file",
	"Filtering out variants by VCF quality",
	"Including variants by VCF quality",
	"Filtering family by Mendelian error rate",
	"Selecting family by Mendelian error rate",
	"Filtering sample by Mendelian error rate",
	"Selecting sample by Mendelian error rate",
	"Filtering variant by Mendelian error rate",
	"Selecting variant by Mendelian error rate",
	"Perform LRT while scoring test scheme",
	"Only use genotypes which are passed QC in VCF file",
	"Only remains phased genotype in VCF file",
	"Only remains unphased genotype in VCF file",
	"Execute WISARD interactively",
	"Check system integrity",
	"Perform regression analysis with interaction between genotype and covariate",
	"Beta value that given to prior",
	"Rho value that given to prior",
	"Phenotype rho value that given to prior",
	"Number of samples that given to prior",
	"Number of phenotypes that given to prior",
	"Do not shuffling sample indices on power calculation",
	"Shuffle case/control split in MDR analysis for train/test scheme",
	"Split dataset along with chromosomes",
	"Merge multiple files",
	"Set merge mode",
	"Test matrix calculation performance",
	"Test special matrix manipulation",
	"Test matrix class integrity",
	"Test functions integrity",
	"Calculate LD for each two-way combinations of variants",
	"Use correlation instead of r in LD computation",
	"Annotate gene information to result",
	"Set annotation range based on its proximity",
	"Annotate variant information to result",
	"Calculate variant density by given window size",
	"A non-whitespace separator from input dataset",
	"Report the distribution of values of given dosage dataset",
	"Incorporates computed PC score as covariates",
	"Do null test in regression",
	"Explore final dataset",
	"Update variant information",
	"Update variant's chromosome",
	"Update variant's name itself",
	"Update variant's genetic distance",
	"Update variant's position",
	"Update genotype",
	"Reference sequence for VCF reference/alternative",
	"Compute LD among adjacent variants",
	"Compute LD among the variants in specific range",
	"Compute LD for given variants",
	"Remove NA from result",
	"Perform BOOST analysis",
	"Set the reporting threshold of BOOST analysis",
	"Perform a quick epistasis test",
	"Perform external analyses",
	"Perform generalized MDR",
	"Impute genotypes",
	"Perform LoD analysis",
	"Export imputed dataset",
	"Do sample selection from given pedigree",
	"Enable R in WISARD",
	"Make random binary phenotype if there is no phenotype",
	"Make random continuous phenotype if there is no phenotype",
	"Compute and report genomic control",
	"User-defined genomic control",
	"Perform p-value adjustment",
	"Report population-specific variants",
	"Report monotone variants",
	"Report singleton variants",
	"Report doubleton variants",
	"Perform gene-MDR analysis",
	"Make specified variants to covariates",
	"Computing inbreeding coefficient",
	"Grouping samples with given condition",
	"Filtering out variants in genic region",
	"Filtering out variants in intergenic region",
	"Make the summary of variants",
	"Reordering samples as specified",
	"Reordering variants as specified",
	"Generate dataset using COSI simulation",
	"Automatically generate set with consecutive variants",
	"Set the size of overlap when generating set",
	"Automatically generate set with random variants",
	"Export retrieved gene-variant mapping",
	"Define exported gene-variant mapping format",
	"Assign case/control string for export dataset",
	"Set 1=case,0=control for export dataset",
	"Do missingness test for each variant by dichotomous trait",
	"Selecting variants by the p-value of missingness test",
	"Removing variants by the p-value of missingness test",
	"Make given portion or number of genotype as missing",
	"Make given portion or number of phenotype as missing",
	"Make variant-filtering report",
	"Set the name of field for hard-coded genotype",
	"Export only phenotype when make dataset",
	"Do not export header when make dataset",
	"Export dataset as Merlin format",
	"Generate MDS plot points and draw if possible",
	"Investigate GxG to available analyses",
	"Transform phenotype using inverse normalization",
	"Forcely converge LMM",
	"Export A/C/G/T coded genotype to 1/2/3/4",
	"Export A/C/G/T coded genotype to other coding",
	"Perform DFAM analysis",
	"Update allele information",
	"Set the window for various options",
	"Print out confidence interval if available",
	"Perform LASSO",
	"Cut lambda value for LASSO",
	"Report all results in LASSO",
	"Perform regression with partial least square (PLS)",
	"Provide sample-wise weight for LMM",
	"Set the number of lines to be skipped",
	"Allow single-side parent in the retrieval of pedigree structure",
	"Retrieve BCF file as an input",
	"Use LOOCV in training-testing scheme",
	"MDR reporting threshold",
	"Perform LD contrast method",
	"Flip the genotype or change with custom flip order",
	"Assign a list of variants to be flipped",
	"Compute sample-site het/hom ratio",
	"Validate marker information based on dbSNP",
	"Perform meta-analysis",
	"Force FID to fixed one", /* FIXME : Need to implement */
	"Output type of table type outputs",
	"Pre-defined MAF",
	"Compute heterozygosity",
	"Compute variant-level BLUP",
	"Export final dataset by FID",
	"Gene plus-minus span for mapping variant when gene definition is ranged form",

	"Perform truncated ridge analysis",
	"Calculate relatedness matrix using Hamming distance",
	"Calculate relatedness matrix using Balding-Nichols model",

	"Perform fuzzy MDR analysis",
	"Perform GxG analysis with all variables into a single model",
	"List of gene pairs to test gene-gene interaction",
	"List of lambdas of gene/pathway to try for gene-gene interaction",

	"Prune variants by Variance Inflation Factor",
	"Prune variants by pair-wise manner",
#endif

	/*** STEP 3 : Add the description of new option here */
};

cOption C_option;

void cOption::_loadConf()
{
	char S_buf[512];
	if (!getItselfPath(S_buf)) halt("Failed to get the path of binary!");

	/* Trying to open configuration file in default path */
	sprintf(S_buf + strlen(S_buf), WISARD_DIR_LETTER "%s.cnf", TOOLSET_NAME);
	cStrFile C_conf(S_buf, "Configuration file", 1);
	/* If there is no conf, just return */
	if (C_conf.isFailed()) {
		if (!IS_ASSIGNED(data)) {
			getItselfPath(S_buf);
			assign("data", S_buf, 1);
			FORCE_OPT_STRING(data);
		}
		return;
	}

	_runAs(S_buf);
	/* Set the datapath if not assigned */
	if (!IS_ASSIGNED(data)) {
		getItselfPath(S_buf);
		assign("data", S_buf, 1);
		FORCE_OPT_STRING(data);
	}
}

cOption& OPTION()
{
	return C_option;
}

cOption::cOption()
{
	if (this != &C_option)
		halt_fmt(WISARD_SYST_INVL_OPTCLASSINIT);
//	halt("This class cannot initiated by default constructor");
}

cOption::cOption(const int N_argc, char *Sa_argv[])
{
	if (this != &C_option)
		halt_fmt(WISARD_SYST_INVL_OPTCLASSINIT);
	init(N_argc, Sa_argv);
}

cOption::~cOption()
{
	clear();
}

void cOption::clear()
{
	/* Free if string value is assigned */
	for (wsUint i=0 ; i<L_defOpt ; i++) {
		if (Xa_defOption[i].E_type == OT_STRING &&
			Xa_defOption[i].S_strVal != NULL)
			MFREE(Xa_defOption[i].S_strVal);
	}
	if (OPT_STRING(cmdLine)) {
		char *tmp = OPT_STRING(cmdLine);
		delete [] tmp;
		OPT_STRING(cmdLine) = NULL;
	}
}

int cOption::__procOption(char *S_opt, char *Sp_prev)
{
	/* For all available options */
	for (wsUint i=0 ; i<L_defOpt ; i++) {
		xOption& X_currOpt = Xa_defOption[i];

		/* Check both of long and short name */
		char B_depUsed = 0;
		if (!strcmp(S_opt, X_currOpt.S_longName) ||
			!strcmp(S_opt, X_currOpt.S_shortName) ||
			(B_depUsed = X_currOpt.S_longSynonym && !strcmp(S_opt, X_currOpt.S_longSynonym))) {
			if (B_depUsed)
				LOGwarn("Option [%s] is deprecated, use [%s] instead\n",
					S_opt, X_currOpt.S_longName);


			/* Is this option requires value? */
			if (X_currOpt.E_type == OT_ONOFF) {
				__procValue(i, (char *)"1");
				return 0; /* 0 == NO MODE */
			} else
				return (i+1); /* Because 0 == NO MODE */
		}
	}

	/* Not valid option */
	if (Sp_prev == NULL)
		halt_fmt(WISARD_INVL_OPTNAME, S_opt);
	else if (Sp_prev[strlen(Sp_prev)-1] == ',') /* Maybe spaced comma */
		halt_fmt(WISARD_INVL_OPTVAL_SC, S_opt, Sp_prev, S_opt);
	else
		halt_fmt(WISARD_INVL_OPTNAME, S_opt);
	return -1;
}

int cOption::__procValue(int N_idxOpt, const char *S_oval, char B_over/*=0*/)
{
	int			N_ret	= 0;
	xOption&	X_opt	= Xa_defOption[N_idxOpt];
	char		*Sp_parseEnd = NULL;

	if (B_over==0 && X_opt.B_isAssigned == true)
		halt_fmt(WISARD_DUPL_OPTION, X_opt.S_longName);

	/* If value starts with --..., give default value */
	if (S_oval[0] && S_oval[0] == '-' && S_oval[1] && S_oval[1] == '-' && S_oval[2])
		N_ret = 1;
	wsStrCst S_val = N_ret ? getDefVal(X_opt.S_longName+2) : S_oval;

	switch (X_opt.E_type) {
	default:
		halt_fmt(WISARD_SYST_INVL_OPTTYPE, X_opt.E_type, X_opt.S_longName);
		break;
	case OT_ONOFF:
		/* If ONOFF option have INVALID value */
		if (S_val[1] || (S_val[0] != '0' && S_val[0] != '1'))
			halt_fmt(WISARD_SYST_INVL_BOOLOPTVAL, S_val, X_opt.S_longName);
//		LOG("[%s(`%s`)] enabled\n", Sa_descOption[N_idxOpt],
//			X_opt.S_longName);
	case OT_NUMBER:
		X_opt.N_intVal = (int)strtol(S_val, &Sp_parseEnd, 10);
		if (Sp_parseEnd[0] != '\0') {
			wsReal R_tmp = str2dbl(S_val, &Sp_parseEnd);
			if (Sp_parseEnd[0] == '\0' && (wsReal)(int)R_tmp == R_tmp) {
				X_opt.N_intVal = (int)R_tmp;
				break;
			}
			/* Parse special case : [0-9]+K [0-9]+M ... */
			char		S_st	= S_val[0];
			const char	*S_p	= NULL;

			/* Find the end of plain number */
			for (S_p=S_st=='-'?S_val+1:S_val ; *S_p>='0'&&*S_p<='9' ; S_p++)
				if (*S_p == '\0') break;

			/* Is the end of plain number have multiplier key? */
			if (*S_p == 'k' || *S_p == 'K' || *S_p == 'm' || *S_p == 'M' ||
				*S_p == 'g' || *S_p == 'G') {
				/* No extra character is allowed */
				if (*(S_p+1))
					halt_fmt(WISARD_SYST_INVL_NUMBEROPTVAL, S_val, X_opt.S_longName);

				/* Multiply given value as the multiplier */
				if (*S_p == 'k' || *S_p == 'K')
					X_opt.N_intVal *= 1000;
				else if (*S_p == 'm' || *S_p == 'M')
					X_opt.N_intVal *= 1000000;
				else if (*S_p == 'g' || *S_p == 'G')
					X_opt.N_intVal *= 1000000000;
			} else
				halt_fmt(WISARD_SYST_INVL_NUMBEROPTVAL, S_val, X_opt.S_longName);
		}
		/* Check range */
		switch (X_opt.X_rng) {
		case RT_ONOFF: break;
		case RT_PROPNUM:
			halt("SYSERR: This option has a bug (prop/num type cannot be an integer type"); break;
		case RT_ABS1: if (X_opt.N_intVal < -1 || X_opt.N_intVal > 1)
			halt_fmt(WISARD_INVL_RANGE_NUMBER, X_opt.S_longName, "[-1, 1]", X_opt.N_intVal); break;
		case RT_NONNEG: if (X_opt.N_intVal < 0)
			halt_fmt(WISARD_INVL_RANGE_NUMBER, X_opt.S_longName, "[0, inf)", X_opt.N_intVal); break;
		case RT_POS: if (X_opt.N_intVal <= 0)
			halt_fmt(WISARD_INVL_RANGE_NUMBER, X_opt.S_longName, "(0, inf)", X_opt.N_intVal); break;
		case RT_NUM: break;
		default:
			halt_fmt(WISARD_SYST_INVL_RANGEDEF, X_opt.S_longName, X_opt.X_rng); break;
		}
		break;
	case OT_STRING: {
			char *Sp_newVal = NULL;
			/* Looks like an option... then */
			if (S_val[2] && S_val[0] == '-' && S_val[1] == '-') {
				/* If not null default value, assign! */
				if (X_opt.S_defValue[0]) {
					wsAlloc(Sp_newVal, char, strlen(X_opt.S_defValue)+1);
					strcpy(Sp_newVal, X_opt.S_defValue);
					N_ret = 1;
				} else {
					halt_fmt(WISARD_NULL_OPTVALUE, X_opt.S_longName);
				}
			} else {
				/* If null default value, do error */
				if (!S_val || (S_val[2] && S_val[0] == '-' && S_val[1] == '-'))
					halt_fmt(WISARD_NULL_OPTVALUE, X_opt.S_longName);
				wsAlloc(Sp_newVal, char, strlen(S_val) + 1);
				strcpy(Sp_newVal, S_val);
			}

//			LOG("[%s(`%s`)] allocated value `%s`\n", Sa_descOption[N_idxOpt],
//				X_opt.S_longName, S_val);
			X_opt.S_strVal = Sp_newVal;
		} break;
	case OT_REAL:
		X_opt.R_realVal = (wsReal)str2dbl(S_val, &Sp_parseEnd);
		if (Sp_parseEnd[0] != '\0')
			halt_fmt(WISARD_SYST_INVL_REALOPTVAL, S_val, X_opt.S_longName);
		/* Check range */
		switch (X_opt.X_rng) {
		case RT_PROPNUM: {
				wsReal R_intVal = rint(X_opt.R_realVal);
				if ((X_opt.R_realVal > W0 && X_opt.R_realVal < W1) ||
					(X_opt.R_realVal >= W1 && R_intVal == X_opt.R_realVal)) {
				} else halt_fmt(WISARD_INVL_RANGE_REAL, X_opt.S_longName, "(0,1) or integer >= 1", X_opt.R_realVal);
			} break;
		case RT_ABS1: if (X_opt.R_realVal < REAL_CONST(-1.0) || X_opt.R_realVal > W1)
			halt_fmt(WISARD_INVL_RANGE_REAL, X_opt.S_longName, "[-1, 1]", X_opt.R_realVal); break;
		case RT_NONNEG: if (X_opt.R_realVal < W0)
			halt_fmt(WISARD_INVL_RANGE_REAL, X_opt.S_longName, "[0, inf)", X_opt.R_realVal); break;
		case RT_POS: if (X_opt.R_realVal <= W0)
			halt_fmt(WISARD_INVL_RANGE_REAL, X_opt.S_longName, "(0, inf)", X_opt.R_realVal); break;
		case RT_01: if (X_opt.R_realVal < W0 || X_opt.R_realVal > W1)
			halt_fmt(WISARD_INVL_RANGE_REAL, X_opt.S_longName, "[0, 1]", X_opt.R_realVal); break;
		case RT_M01: if (X_opt.R_realVal <= W0 || X_opt.R_realVal > W1)
			halt_fmt(WISARD_INVL_RANGE_REAL, X_opt.S_longName, "(0, 1]", X_opt.R_realVal); break;
		case RT_NUM: break;
		default:
			halt_fmt(WISARD_SYST_INVL_RANGEDEF, X_opt.S_longName, X_opt.X_rng); break;
		}
		break;
	case OT_EXPR: {
			/* Do parse */
			xOperation *X_op = parse(S_val);
			memcpy(&(X_opt.X_expVal), X_op, sizeof(xOperation));
			delete X_op;
		} break;		
	case OT_RANGE:
		procRange(&X_opt, S_val);
		switch (X_opt.X_rng) {
		case RT_POS: if ((X_opt.X_rngVal.R_s <= W0 && X_opt.X_rngVal.R_sEQ) ||
						 (X_opt.X_rngVal.R_e <= W0 && X_opt.X_rngVal.R_eEQ))
			halt_fmt(WISARD_INVL_RANGE_STRING, X_opt.S_longName, "(0, inf)", S_val);
		default:
			break;
		}
		break;
	}

	/* Make it assigned */
	//	notice("Option `%s` assigned", X_opt.S_longName);
	if (!X_opt.B_isAssigned)
		Xv_asgnOpts.push_back(&X_opt);
	X_opt.B_isAssigned = true;
	return N_ret;
}

void cOption::assign(const char *S_name, const char *S_val/*=NULL*/,
	char B_over/*=0*/)
{
	wsUint N_idx = findOption(S_name);
	if (S_val == NULL)
		S_val = getDefVal(S_name);
	__procValue(N_idx, S_val, B_over);
}

int cOption::init(const int N_argc, char *Sa_argv[])
{
	int N_mode = 0; /* 0 for optName, otherwise optValue for `N_mode`th option */

	L_defOpt = len(Xa_defOption, xOption);

	/* Check option integrity */
#ifdef USE_OPTCHECK
	_checkOptionIntegrity();
#endif

	/* Print machine name */ {
		char S_mname[256];
		getMachineName(S_mname);
		LOG("Working machine : %s\n", S_mname);
	}

	/* Print working directory */
	char S_curPath[FILENAME_MAX];
	if (GetCurrentDir(S_curPath, FILENAME_MAX))
		LOG("Working directory : %s\n", S_curPath);

	/* Print run command */ {
		char *p = Sa_argv[0] + strlen(Sa_argv[0]) - 1;
		while (p>Sa_argv[0] && *p!=WISARD_DIR_CHAR) p--;
		p = p==Sa_argv[0] ? Sa_argv[0] : p+1;

		LOG("Run %s as :", p);
	}

	string S_cmdLine;
	if (N_argc == 1) {
		S_cmdLine = "<NO ARGUMENT>";
		LOGnf(" <NO ARGUMENT>");
	} else for (int i = 1; i < N_argc; i++) {
		LOGnf(" %s", Sa_argv[i]);
		S_cmdLine += Sa_argv[i];
		if ((i+1) < N_argc) S_cmdLine += " ";
	}
	char* tmp = NULL;
	wsAlloc(tmp, char, strlen(S_cmdLine.c_str())+1);
	strcpy(tmp, S_cmdLine.c_str());
	OPT_STRING(cmdLine) = tmp;
	LOGnf("\n\n");

	/* For all arguments EXCEPT Sa_argv[0] */
	for (int i=1 ; i<N_argc ; i++) {
		/* Remove whitespace */
		for (char *p=Sa_argv[i] + strlen(Sa_argv[i]) - 1 ; *p == '\r' || *p =='\n' ; p--)
			*p = '\0';
		//LOG("Argument [%d] : [%s]\n", i, Sa_argv[i]);
		if (N_mode == 0)
			N_mode = __procOption(Sa_argv[i], Sa_argv[i-1]);
		else {
			if (__procValue(N_mode-1, Sa_argv[i])) /* N_mode(1-based) into 0-based */
				N_mode = __procOption(Sa_argv[i], Sa_argv[i-1]);
			else
				N_mode = 0;
		}
	}
	/* Check all of option values were appropriately assigned */
	if (N_mode != 0) {
		wsStrCst V = getDefVal(Sa_argv[N_argc-1] + 2);
		if (V[0])
			__procValue(N_mode-1, V);
		else
			halt_fmt(WISARD_NULL_OPTVALUE, Sa_argv[N_argc-1]);
	}

	/*** VERBOSE OPTION WILL WORK BELOW CODE ***/

	/* Configuration file process */
	_loadConf();

	/* If --script given, read file and alter it as */
	ASG_OPT_STRING(script);
	if (IS_ASSIGNED(script))
		_runAs(OPT_STRING(script));

	/* Assigns default value to NOT-assigned-by-user options */
	for (wsUint i=0 ; i<L_defOpt ; i++) if (Xa_defOption[i].B_isAssigned == false) {
		/* IF THIS IS ESSENTIAL, raise ERROR */
// 		if (Xa_defOption[i].B_isEssential == true)
// 			halt("Option `%s` should be assigned, but not given", Xa_defOption[i].S_longName);

		/* Assign default value only it have */
		if (Xa_defOption[i].B_isApplyDef == true && Xa_defOption[i].S_defValue != NULL) {
			LOG("[%s(`%s`)] allocated default value `%s`\n", Sa_descOption[i], Xa_defOption[i].S_longName, Xa_defOption[i].S_defValue);
			__procValue(i, Xa_defOption[i].S_defValue);
		}
	}

	/* Do interaction */
#if ((TOOLSET_TYPE & 0x100) == 0x100)
	ASG_OPT_NUMBER(interactive);
	if (OPT_ENABLED(interactive))
		interactiveExecution();
#endif

	/* Assigns value into class */

	/* OPTIONS FOR PROGRAM FLOW */
	ASG_OPT_NUMBER(verbose);		/**/
	ASG_OPT_NUMBER(seed);			/*%*/
	ASG_OPT_STRING(script);			/*%*/
	ASG_OPT_NUMBER(thread);			/**/// Number of threads utilizing
	ASG_OPT_STRING(out);			/**/// Input dataset file name

	/* OPTIONS FOR DATA INPUT */
	ASG_OPT_STRING(lgen);			/*%*/
	ASG_OPT_STRING(tped);			/*%*/
	ASG_OPT_STRING(data);
	ASG_OPT_STRING(bed);			/*%*/
	ASG_OPT_STRING(fam);			/*%*/

	/* OPTIONS FOR ADDITIONAL INPUT */
	ASG_OPT_STRING(bim);			/*%*/
	ASG_OPT_STRING(map);			/*%*/
	ASG_OPT_STRING(pheno);			/**/// Alternative phenotype file name
	ASG_OPT_STRING(sampvarflag);
	ASG_OPT_STRING(geneset);		/**/
	ASG_OPT_STRING(set);			/**/

	/* OPTIONS FOR HANDLING ADDITIONAL INPUT */
	ASG_OPT_NUMBER(ignorefid);		/*%*/
	ASG_OPT_NUMBER(nofid);			/*%*/
	ASG_OPT_NUMBER(noparent);		/*%*/
	ASG_OPT_NUMBER(nosex);			/*%*/
	ASG_OPT_NUMBER(nopos);			/**/
	ASG_OPT_NUMBER(nogdist);		/**/
	ASG_OPT_STRING(acgt);			/*%*/
	ASG_OPT_NUMBER(nomap);			/*%*/
	ASG_OPT_NUMBER(1234);			/*%*/
	ASG_OPT_STRING(sepallele);		/*%*/
	ASG_OPT_NUMBER(consecallele);	/*%*/
	ASG_OPT_NUMBER(1case);			/*%*/
	ASG_OPT_STRING(mispheno);		/*%*/
	ASG_OPT_STRING(misgeno);		/*%*/

	/* ELEMENTARY ANALYSIS OPTIONS */
	ASG_OPT_STRING(freq);			/*%*/
	ASG_OPT_STRING(hwe);			/*%*/
	ASG_OPT_NUMBER(heritability);	/*%*/

	/* GENE-BASED TESTS */
	ASG_OPT_NUMBER(noweight);		/*%*/
	ASG_OPT_STRING(betaweight);		/*%*/
	ASG_OPT_STRING(weight);			/*%*/
	ASG_OPT_NUMBER(mafweight);
	ASG_OPT_NUMBER(gsmacthr);		/**/

	/* VARIANT-LEVEL TESTS */
	ASG_OPT_STRING(prevalence);		/**/
	ASG_OPT_NUMBER(scoretest);		/*%*/

	/* SAMPLE RELATEDNESS OPTIONS */
	ASG_OPT_NUMBER(kinship);		/**/// Use PDDT instead of correlation calculation?
	ASG_OPT_NUMBER(ibs);			/*%*/
	ASG_OPT_NUMBER(ktau);			/*%*/
	ASG_OPT_NUMBER(empktau);		/*%*/
	ASG_OPT_NUMBER(corpearson);		/**/
	ASG_OPT_NUMBER(cordiag1);		/**/
	ASG_OPT_NUMBER(medcor);			/**/
	ASG_OPT_NUMBER(indep);			/**/
	ASG_OPT_NUMBER(hybrid);			/*%*/
	ASG_OPT_STRING(cor);			/**/
	ASG_OPT_NUMBER(corpair);		/*%*/
	ASG_OPT_NUMBER(corgrm);
	ASG_OPT_NUMBER(corepacts);		/*%*/
	ASG_OPT_NUMBER(x);				/**/
	ASG_OPT_NUMBER(x2);				/**/
	ASG_OPT_RANGE(cormaf);			/**/

	/* EXPORT-RELATED OPTIONS */

	ASG_OPT_NUMBER(quiet);
	ASG_OPT_NUMBER(version);		/*%*/
	ASG_OPT_STRING(ped);			/*%*/
	ASG_OPT_STRING(variantvar);		// Marker variable name
	ASG_OPT_EXPR(filvariant);
	ASG_OPT_EXPR(incvariant);
	ASG_OPT_EXPR(filgeno);			/*%*/
	ASG_OPT_EXPR(incgeno);			/*%*/
	ASG_OPT_STRING(pname);			/*%*/
	ASG_OPT_STRING(remsamp);		/*%*/
	ASG_OPT_STRING(selsamp);		/*%*/
	ASG_OPT_STRING(nasamp);			/*%*/
	ASG_OPT_REAL(randnasamp);		/*%*/
	ASG_OPT_STRING(remvariant);		/*%*/
	ASG_OPT_STRING(selvariant);		/*%*/
	ASG_OPT_STRING(remfam);			/*%*/
	ASG_OPT_STRING(selfam);			/*%*/
	ASG_OPT_STRING(filrange);		/*%*/
	ASG_OPT_STRING(incrange);		/*%*/
	ASG_OPT_REAL(varresize);		/*%*/
	ASG_OPT_NUMBER(varwindow);		/*%*/
	ASG_OPT_REAL(sampresize);		/*%*/
	ASG_OPT_NUMBER(emcount);		/**/
	ASG_OPT_REAL(aithr);			/**/
	ASG_OPT_NUMBER(autoonly);		/*%*/
	ASG_OPT_NUMBER(sexonly);		/*%*/
	ASG_OPT_STRING(chr);			/*%*/
	ASG_OPT_NUMBER(makeblup);		/**/
	ASG_OPT_RANGE(filmaf);			/*%*/
	ASG_OPT_RANGE(filmac);			/*%*/
	ASG_OPT_RANGE(filhwe);			/*%*/
	ASG_OPT_RANGE(filgind);			/*%*/
	ASG_OPT_RANGE(filgvar);			/*%*/
	ASG_OPT_NUMBER(filnf);			/*%*/
	ASG_OPT_NUMBER(filmf);			/*%*/
	ASG_OPT_NUMBER(filcase);		/*%*/
	ASG_OPT_NUMBER(filcontrol);		/*%*/
	ASG_OPT_NUMBER(filmispheno);	/*%*/
	ASG_OPT_EXPR(filsample);		/**/
	ASG_OPT_RANGE(filgdist);		/*%*/
	ASG_OPT_NUMBER(filnosex);		/*%*/
	ASG_OPT_NUMBER(filmale);		/*%*/
	ASG_OPT_NUMBER(filfemale);		/*%*/
	ASG_OPT_NUMBER(snvonly);		/*%*/
	ASG_OPT_NUMBER(indelonly);		/*%*/
	ASG_OPT_RANGE(incmaf);			/*%*/
	ASG_OPT_RANGE(incmac);			/*%*/
	ASG_OPT_RANGE(inchwe);			/*%*/
	ASG_OPT_RANGE(incgind);			/*%*/
	ASG_OPT_RANGE(incgvar);			/*%*/
	ASG_OPT_RANGE(incgdist);		/*%*/
	ASG_OPT_EXPR(incsample);
	ASG_OPT_NUMBER(ml);				/*%*/
	ASG_OPT_STRING(blup);			/**/
	ASG_OPT_NUMBER(imputepheno);	/**/
	ASG_OPT_STRING(cname);			/*%*/
	ASG_OPT_STRING(makecor);		/**/
	ASG_OPT_NUMBER(empiall);		/**/
	ASG_OPT_STRING(cact);			/*%*/
	ASG_OPT_NUMBER(1sex);			/*%*/
	ASG_OPT_STRING(mafe);			/*%*/
	ASG_OPT_NUMBER(ginv);			/**/
	ASG_OPT_NUMBER(makecov);		/*%*/
	ASG_OPT_NUMBER(makepheno);		/*%*/
	ASG_OPT_NUMBER(founderonly);	/**/
	ASG_OPT_NUMBER(famsummary);		/**/
	ASG_OPT_NUMBER(ncsummary);		/**/
	ASG_OPT_NUMBER(specdcmp);		/*%*/
	ASG_OPT_NUMBER(nostop);			/*%*/
	ASG_OPT_NUMBER(ignoreparent);	/*%*/
	ASG_OPT_STRING(sepid);			/*%*/
	ASG_OPT_NUMBER(nopheno);		/*%*/
	ASG_OPT_STRING(probandcol);		/*%*/
	ASG_OPT_STRING(twincol);		/*%*/
	ASG_OPT_NUMBER(chrwise);		/**/
	ASG_OPT_STRING(est);			/**/
	ASG_OPT_NUMBER(makeclgeno);
	ASG_OPT_NUMBER(indel);			/*%*/
	ASG_OPT_NUMBER(logistic);		/*%*/
	ASG_OPT_STRING(sortvariant);	/*%*/
	ASG_OPT_STRING(sortpos);		/*%*/
	ASG_OPT_STRING(sortsample);		/*%*/
	ASG_OPT_STRING(sortiid);		/*%*/
	ASG_OPT_STRING(baseline);		/*%*/
	ASG_OPT_NUMBER(gz);				/**/
	ASG_OPT_NUMBER(makenrm);
	ASG_OPT_NUMBER(makeweight);		/**/
	ASG_OPT_NUMBER(nophenohdr);	/*%*/
	ASG_OPT_RANGE(pvalrange);		/*%*/
	ASG_OPT_NUMBER(time);			/*%*/
	ASG_OPT_NUMBER(phenostdize);
	ASG_OPT_STRING(outmispheno);	/*%*/
	ASG_OPT_STRING(outmisgeno);		/*%*/
	ASG_OPT_STRING(gxecovs);		/*%*/
	ASG_OPT_NUMBER(passemptyline);	/**/
	ASG_OPT_NUMBER(listvariant);	/*%*/
	ASG_OPT_NUMBER(listsample);		/*%*/
	ASG_OPT_NUMBER(listfounder);	/*%*/
	ASG_OPT_NUMBER(nolmm);			/**/
	ASG_OPT_RANGE(gsetconsec);
	ASG_OPT_NUMBER(nosysfeat);		/**/
	ASG_OPT_STRING(fname);			/*%*/
	ASG_OPT_NUMBER(makeev);			/*%*/
	ASG_OPT_STRING(ev);				/*%*/
	ASG_OPT_NUMBER(natural);		/*%*/
	ASG_OPT_NUMBER(dupnaming);		/*#*/
	ASG_OPT_NUMBER(citation);		/**/
	ASG_OPT_NUMBER(miss);			/*%*/
	ASG_OPT_STRING(misparent);		/*%*/
	ASG_OPT_NUMBER(sampmajor);		/*%*/
	ASG_OPT_NUMBER(nospecdcmp);		/*%*/
	ASG_OPT_STRING(species);
	ASG_OPT_NUMBER(nodata);
	ASG_OPT_NUMBER(regex);
	ASG_OPT_STRING(model);

#if (TOOLSET_TYPE == TOOLSET_MFQLS) || ((TOOLSET_TYPE & 0x100) == 0x100)
	/* VARIANT-LEVEL TESTS */
	ASG_OPT_NUMBER(mqls);			/**/
	ASG_OPT_NUMBER(fqls);			/**/
	ASG_OPT_STRING(heri);			/**/
	ASG_OPT_NUMBER(multifqls);		/*%*/
	ASG_OPT_NUMBER(mqlsconsec);		/**/

	ASG_OPT_NUMBER(fqlsnopddt);		/**/
	ASG_OPT_REAL(retestthr);		/**/
	ASG_OPT_NUMBER(avail);			/*%*/
	ASG_OPT_NUMBER(fastmqls);
	ASG_OPT_NUMBER(fastfqls);
#endif

#if (TOOLSET_TYPE == TOOLSET_FARVAT) || (TOOLSET_TYPE == TOOLSET_MFARVAT) || \
	((TOOLSET_TYPE & 0x100) == 0x100) || (TOOLSET_TYPE == TOOLSET_FARVATX)
	/* GENE-BASED TESTS */
	ASG_OPT_NUMBER(genesummary);	/*%*/
	ASG_OPT_NUMBER(genetest);		/*%*/
	ASG_OPT_REAL(raremaf);			/**/
	ASG_OPT_RANGE(genesize);		/*%*/
	ASG_OPT_NUMBER(pedcmc);			/*%*/
	ASG_OPT_NUMBER(wsum);			/*%*/
	ASG_OPT_NUMBER(kbac);			/*%*/
	ASG_OPT_NUMBER(asum);
	ASG_OPT_NUMBER(farvat);
	ASG_OPT_NUMBER(pedgene);
	ASG_OPT_NUMBER(skato);			/*%*/
	ASG_OPT_NUMBER(skat);			/*%*/
	ASG_OPT_NUMBER(mfhom);			/**/
	ASG_OPT_NUMBER(mfhet);			/**/
	ASG_OPT_NUMBER(adjf1);
	ASG_OPT_NUMBER(adjf2);
	ASG_OPT_NUMBER(farvatx);
	ASG_OPT_REAL(farvatxd);

	ASG_OPT_NUMBER(gmapsummary);	/*%*/
	ASG_OPT_REAL(genemiss);			/*%*/
	ASG_OPT_NUMBER(skatondiv);		/*%*/
	ASG_OPT_STRING(skatodivs);		/*%*/
	ASG_OPT_REAL(kbacalpha);		/*%*/
	ASG_OPT_NUMBER(kbac2side);		/**/
	ASG_OPT_STRING(kbackernel);		/**/
	ASG_OPT_NUMBER(makegeno);
	ASG_OPT_NUMBER(makefarvat);
	ASG_OPT_NUMBER(farvatxndiv);
#endif

#if (TOOLSET_TYPE == TOOLSET_QTEST) || ((TOOLSET_TYPE & 0x100) == 0x100)
	ASG_OPT_NUMBER(qtest);			/*%*/
	ASG_OPT_RANGE(qtestrange);		/**/
	ASG_OPT_NUMBER(qtestclump);		/**/
	ASG_OPT_REAL(qteststt);			/**/
	ASG_OPT_NUMBER(qtestbetacov);	/**/
	ASG_OPT_NUMBER(makebeta);
#endif

#if (TOOLSET_TYPE == TOOLSET_HIMINI) || ((TOOLSET_TYPE & 0x100) == 0x100)
	ASG_OPT_NUMBER(mdr);			/**/
	ASG_OPT_NUMBER(order);			/**/
	ASG_OPT_NUMBER(top);			/**/
	ASG_OPT_NUMBER(hmdr);			/**/
	ASG_OPT_STRING(hmdrall);
	ASG_OPT_STRING(hmdrprior);
#endif

#if (TOOLSET_TYPE == TOOLSET_PHARAOH) || ((TOOLSET_TYPE & 0x100) == 0x100)
	ASG_OPT_NUMBER(pharaoh);
	ASG_OPT_NUMBER(proopt);
	ASG_OPT_STRING(prolambda);
	ASG_OPT_RANGE(prorange);
	ASG_OPT_REAL(prothr);
	ASG_OPT_NUMBER(promaxiter);
	ASG_OPT_RANGE(progenesize);
	ASG_OPT_RANGE(progsetsize);
	ASG_OPT_NUMBER(prosingle);
	ASG_OPT_NUMBER(propermcov);
	ASG_OPT_NUMBER(nperm);			/**/
	ASG_OPT_STRING(seqperm);
	ASG_OPT_STRING(permfile);
	ASG_OPT_NUMBER(cv);
	ASG_OPT_NUMBER(gesca);
	ASG_OPT_STRING(modeltype);
	ASG_OPT_STRING(ggpath);
#endif

#if TOOLSET_TYPE == TOOLSET_ONETOOL
	ASG_OPT_NUMBER(debug);
	ASG_OPT_NUMBER(help);
#endif

#if ((TOOLSET_TYPE & 0x100) == 0x100)
	/* OPTIONS FOR DATA INPUT */
	ASG_OPT_STRING(dosage);			/**/
	ASG_OPT_STRING(expression);
	ASG_OPT_STRING(genoprob);
	ASG_OPT_STRING(vcf);			/*%*/

	/* ELEMENTARY ANALYSIS OPTIONS */
	ASG_OPT_STRING(fst);			/**/
	ASG_OPT_NUMBER(mendel);			/*%*/
	ASG_OPT_NUMBER(famuniq);		/*%*/
	ASG_OPT_NUMBER(tstv);			/*%*/
	ASG_OPT_NUMBER(pca);			/*%*/
	ASG_OPT_NUMBER(npc);			/*%*/

	/* GENE-BASED TESTS */
	ASG_OPT_NUMBER(genesplit);
	ASG_OPT_NUMBER(rvtdt);
	ASG_OPT_NUMBER(fbskat);
	ASG_OPT_NUMBER(famvt);			/*%*/
	ASG_OPT_NUMBER(ggemma);
	ASG_OPT_STRING(longitudinal);	/**/
	ASG_OPT_NUMBER(vt);				/*%*/

	/* VARIANT-LEVEL TESTS */
	ASG_OPT_NUMBER(fisher);			/*%*/
	ASG_OPT_NUMBER(trend);			/*%*/
	ASG_OPT_NUMBER(regression);		/*%*/
	ASG_OPT_NUMBER(qls);			/*%*/
	ASG_OPT_NUMBER(tdt);			/*%*/
	ASG_OPT_NUMBER(sdt);			/*%*/
	ASG_OPT_NUMBER(emmax);			/*%*/
	ASG_OPT_NUMBER(gemma);			/*%*/
	ASG_OPT_NUMBER(powercalc);		/**/
	ASG_OPT_NUMBER(powercalc2);

	/* SAMPLE RELATEDNESS OPTIONS */

	/* EXPORT-RELATED OPTIONS */
	ASG_OPT_NUMBER(makeped);		/*%*/
	ASG_OPT_NUMBER(maketped);		/*%*/
	ASG_OPT_NUMBER(makebed);		/*%*/
	ASG_OPT_NUMBER(makeraw);		/*%*/
	ASG_OPT_NUMBER(makedom);		/*%*/
	ASG_OPT_NUMBER(makerec);		/*%*/
	ASG_OPT_NUMBER(makevcf);		/*%*/
	ASG_OPT_NUMBER(makebcf);		/*%*/
	ASG_OPT_NUMBER(makelgen);		/*%*/
	ASG_OPT_NUMBER(makegen);		/*%*/
	ASG_OPT_NUMBER(makebgen);		/*%*/
	ASG_OPT_NUMBER(makebeagle);		/**/

	/* S.A.G.E. OPTIONS */
	ASG_OPT_NUMBER(relpair);
	ASG_OPT_NUMBER(fcor);
	ASG_OPT_NUMBER(segreg);
	ASG_OPT_NUMBER(lodlink);
	ASG_OPT_STRING(typ);
	ASG_OPT_STRING(par);
	ASG_OPT_NUMBER(fcorStdErrOff);
	ASG_OPT_NUMBER(lodlinkLinkageTestOff);
	ASG_OPT_NUMBER(lodlinkLinkageHomogOff);
	ASG_OPT_NUMBER(lodlinkLinkageSexSpecific);
	ASG_OPT_NUMBER(lodlinkSmithHomogTest);
	ASG_OPT_NUMBER(lodlinkGenotypes);
	ASG_OPT_NUMBER(merlin);
	ASG_OPT_NUMBER(plot);

	ASG_OPT_NUMBER(simtrio);		/**/// Use trio simulation instead of data input?
	ASG_OPT_NUMBER(szfam);			/**/// # of families in the simulation
	ASG_OPT_NUMBER(szvar);			/**/// Use trio simulation instead of data input?
	ASG_OPT_NUMBER(simfam);			/*%*/
	ASG_OPT_NUMBER(trio);			/*%*/
	ASG_OPT_NUMBER(extfam);			/*%*/
	ASG_OPT_NUMBER(nsig);			/**/
	ASG_OPT_REAL(sigmaf);			/**/
	ASG_OPT_STRING(simfreq);		/**/
	ASG_OPT_STRING(mafvar);			/**/
	ASG_OPT_REAL(proppc);			/*%*/
	ASG_OPT_NUMBER(usemf);			/*%*/
	ASG_OPT_NUMBER(fullpca);		/**/// Use full eigenreduction in PCA?
	ASG_OPT_STRING(rpath);			/**/// A path to R
	ASG_OPT_STRING(grm);			/**/
	ASG_OPT_REAL(grmalpha);			/**/
	ASG_OPT_NUMBER(makemdr);		/**/
	ASG_OPT_NUMBER(zipbgen);		/*%*/
	ASG_OPT_RANGE(filqual);			/*%*/
	ASG_OPT_RANGE(incqual);			/*%*/
	ASG_OPT_RANGE(filmendelfam);	/*%*/
	ASG_OPT_RANGE(incmendelfam);	/*%*/
	ASG_OPT_RANGE(filmendelsamp);	/*%*/
	ASG_OPT_RANGE(incmendelsamp);	/*%*/
	ASG_OPT_RANGE(filmendelvar);	/*%*/
	ASG_OPT_RANGE(incmendelvar);	/*%*/
	ASG_OPT_NUMBER(lrt);			/**/
	ASG_OPT_NUMBER(vcfqc);			/*%*/
	ASG_OPT_NUMBER(phasedonly);		/*%*/
	ASG_OPT_NUMBER(unphasedonly);	/*%*/
	ASG_OPT_NUMBER(interactive);	/**/
	ASG_OPT_NUMBER(check);
	ASG_OPT_NUMBER(gxe);			/*%*/
	ASG_OPT_REAL(beta);				/**/
	ASG_OPT_REAL(rho);				/**/
	ASG_OPT_REAL(rhopheno);			/**/
	ASG_OPT_RANGE(nsamp);			/**/
	ASG_OPT_NUMBER(npheno);			/**/
	ASG_OPT_NUMBER(noshuffle);		/**/
	ASG_OPT_NUMBER(shuffle);		/**/
	ASG_OPT_NUMBER(split);			/**/
	ASG_OPT_STRING(merge);			/*%*/
	ASG_OPT_NUMBER(mergemode);		/*%*/
	ASG_OPT_NUMBER(testmatrix);
	ASG_OPT_NUMBER(testmatfunc);
	ASG_OPT_NUMBER(testmatclass);
	ASG_OPT_NUMBER(testfunc);
	ASG_OPT_NUMBER(ld);				/*%*/
	ASG_OPT_NUMBER(ldcor);			/*%*/
	ASG_OPT_STRING(annogene);		/**/
	ASG_OPT_NUMBER(annorange);		/**/
	ASG_OPT_STRING(annovar);
	ASG_OPT_NUMBER(density);		/**/
	ASG_OPT_STRING(sep);
	ASG_OPT_NUMBER(dsgdist);		/**/
	ASG_OPT_NUMBER(pc2cov);			/*%*/
	ASG_OPT_NUMBER(donull);
	ASG_OPT_NUMBER(explore);
	ASG_OPT_STRING(updvariant);		/**/
	ASG_OPT_STRING(updchr);			/**/
	ASG_OPT_STRING(updname);		/**/
	ASG_OPT_STRING(updgdist);		/**/
	ASG_OPT_STRING(updpos);			/**/
	ASG_OPT_STRING(updgeno);		/**/
	ASG_OPT_STRING(ref);			/**/
	ASG_OPT_NUMBER(ldbin);			/*%*/
	ASG_OPT_NUMBER(ldsize);			/*%*/
	ASG_OPT_STRING(ldvar);			/*%*/
	ASG_OPT_NUMBER(remna);			/*%*/
	ASG_OPT_NUMBER(boost);			/*%*/
	ASG_OPT_REAL(thrboost);			/*%*/
	ASG_OPT_NUMBER(quickepi);		/**/
	ASG_OPT_STRING(ext);
	ASG_OPT_NUMBER(gmdr);			/**/
	ASG_OPT_NUMBER(impute);			/**/
	ASG_OPT_NUMBER(lod);			/**/
	ASG_OPT_NUMBER(makeimpute);		/**/
	ASG_OPT_NUMBER(sxa);
	ASG_OPT_STRING(R);				/**/
	ASG_OPT_NUMBER(randbinpheno);	/**/
	ASG_OPT_NUMBER(randpheno);		/**/
	ASG_OPT_NUMBER(genoctrl);		/*%*/
	ASG_OPT_REAL(usergc);			/**/
	ASG_OPT_NUMBER(adjust);			/**/
	ASG_OPT_STRING(popuniq);		/*%*/
	ASG_OPT_NUMBER(monotone);		/*%*/
	ASG_OPT_NUMBER(singleton);		/*%*/
	ASG_OPT_NUMBER(doubleton);		/*%*/
	ASG_OPT_NUMBER(genemdr);
	ASG_OPT_STRING(variant2cov);	/*%*/
	ASG_OPT_NUMBER(inbreed);		/*%*/
	ASG_OPT_STRING(group);
	ASG_OPT_NUMBER(filgenic);
	ASG_OPT_NUMBER(filintergenic);
	ASG_OPT_NUMBER(variantsummary);	/*%*/
	ASG_OPT_STRING(sampleorder);	/*%*/
	ASG_OPT_STRING(variantorder);	/*%*/
	ASG_OPT_STRING(cosi);
	ASG_OPT_NUMBER(setconsec);		/**/
	ASG_OPT_NUMBER(setoverlap);		/**/
	ASG_OPT_RANGE(setrandom);
	ASG_OPT_NUMBER(makeset);		/*%*/
	ASG_OPT_STRING(outcact);		/*%*/
	ASG_OPT_NUMBER(out1case);		/*%*/
	ASG_OPT_NUMBER(settype);		/*%*/
	ASG_OPT_NUMBER(mistest);		/*%*/
	ASG_OPT_RANGE(incmistest);		/*%*/
	ASG_OPT_RANGE(filmistest);		/*%*/
	ASG_OPT_REAL(nageno);			/*%*/
	ASG_OPT_REAL(napheno);
	ASG_OPT_NUMBER(filtreport);
	ASG_OPT_STRING(genofield);
	ASG_OPT_NUMBER(outphenoonly);	/*%*/
	ASG_OPT_NUMBER(outnoheader);	/*%*/
	ASG_OPT_NUMBER(makemerlin);		/*%*/
	ASG_OPT_NUMBER(mds);			/*%*/
	ASG_OPT_NUMBER(gxg);
	ASG_OPT_NUMBER(invnorm);		/**/
	ASG_OPT_NUMBER(forceconv);
	ASG_OPT_NUMBER(out1234);		/**/
	ASG_OPT_STRING(outacgt);		/**/
	ASG_OPT_NUMBER(dfam);
	ASG_OPT_STRING(updallele);
	ASG_OPT_RANGE(window);
	ASG_OPT_REAL(ci);
	ASG_OPT_NUMBER(lasso);			/**/
	ASG_OPT_NUMBER(lassoall);
	ASG_OPT_REAL(lassolambda);		/**/
	ASG_OPT_NUMBER(pls);
	ASG_OPT_STRING(sampleweight);
	ASG_OPT_NUMBER(nskip);
	ASG_OPT_NUMBER(singleparent);
	ASG_OPT_STRING(bcf);			/**/
	ASG_OPT_NUMBER(loocv);
	ASG_OPT_RANGE(mdrthr);			/**/
	ASG_OPT_NUMBER(ldcontrast);
	ASG_OPT_STRING(flip);
	ASG_OPT_STRING(varsubset);
	ASG_OPT_NUMBER(hethom);			/**/
	ASG_OPT_NUMBER(markercheck);
	ASG_OPT_STRING(meta);
	ASG_OPT_STRING(fid);
	ASG_OPT_STRING(outformat);
	ASG_OPT_STRING(maf);
	ASG_OPT_NUMBER(het);
	ASG_OPT_REAL(variantblup);
	ASG_OPT_NUMBER(famsplit);
	ASG_OPT_NUMBER(setspan);

	ASG_OPT_NUMBER(tridge);
	ASG_OPT_NUMBER(hamming);
	ASG_OPT_NUMBER(bn);

	ASG_OPT_NUMBER(fuzzymdr);
	ASG_OPT_NUMBER(gxgall);
	ASG_OPT_STRING(gxglist);
	ASG_OPT_STRING(gxglambda);

	ASG_OPT_STRING(prunevif);
	ASG_OPT_STRING(prunepw);
#endif

	/*** STEP 2 : Add the assignment of new option here */

	/* Option post-processing */
	if (OPT_ENABLED(1case)) {
		if (IS_ASSIGNED(cact))
			halt_fmt(WISARD_CANT_EXCL_OPT, "--1case", "--cact");
		OPT_NUMBER(phenoCase) = 1;
		OPT_NUMBER(phenoCtrl) = 0;
	} else if (IS_ASSIGNED(cact)) {
		/* Load 2 integer value */
		wsUint N_v		= 0;
		wsReal *Ra_cact	= loadRealValues(OPT_STRING(cact), &N_v);

		/* Number of values, and format check */
		if (N_v != 2 || (int)(wsUint)Ra_cact[0] != Ra_cact[0] ||
			(int)(wsUint)Ra_cact[1] != Ra_cact[1])
			halt("Argument of --cact [%s] must be two integers divided by "
				"comma[,] with no whitespaces", OPT_STRING(cact));

		/* Assign values into case and control SEQUENTIALLY */
		OPT_NUMBER(phenoCase) = (int)Ra_cact[0];
		OPT_NUMBER(phenoCtrl) = (int)Ra_cact[1];

		sseFree(Ra_cact);
	} else {
		OPT_NUMBER(phenoCase) = 2;
		OPT_NUMBER(phenoCtrl) = 1;
	}
	LOG("Case phenotype is [%d], control phenotype is [%d]\n",
		OPT_NUMBER(phenoCase), OPT_NUMBER(phenoCtrl));

	/* --outformat */
	if (IS_ASSIGNED(outformat)) {
		wsStrCst S_fmt = OPT_STRING(outformat);
		if (!stricmp(S_fmt, "tabular")) {
			/* Do nothing */
		} else if (!stricmp(S_fmt, "json")) {
			LOGnote("Output will be provided as JSON format\n");
		} else if (!stricmp(S_fmt, "xml")) {
			LOGnote("Output will be provided as XML format\n");
		} else if (!stricmp(S_fmt, "csv")) {
			LOGnote("Output will be provided as CSV format\n");
		} else
			halt("Invalid --outformat value [%s]", S_fmt);
	}

	if (OPT_ENABLED(1sex)) {
		if (IS_ASSIGNED(mafe))
			halt_fmt(WISARD_CANT_EXCL_OPT, "--1sex", "--mafe");
		OPT_NUMBER(sexMale) = '0';
		OPT_NUMBER(sexFema) = '1';
	} else if (IS_ASSIGNED(mafe)) {
		/* Load 2 character value */
		wsUint	N_v		= 0;
		char**	Sa_mafe	= loadStringValues(OPT_STRING(mafe), &N_v);

		/* Number of values, and format check */
		if (N_v != 2 || Sa_mafe[0][1] || Sa_mafe[1][1])
			halt("Argument of --mafe [%s] must be two single characters divided by "
				"comma[,] with no whitespaces", OPT_STRING(mafe));

		/* Assign values into case and control SEQUENTIALLY */
		OPT_NUMBER(sexMale) = Sa_mafe[0][0];
		OPT_NUMBER(sexFema) = Sa_mafe[1][0];

		sseFree(Sa_mafe);
	} else {
		OPT_NUMBER(sexMale) = '1';
		OPT_NUMBER(sexFema) = '2';
	}
	LOG("Male coding is [%c], female coding is [%c]\n",
		OPT_NUMBER(sexMale), OPT_NUMBER(sexFema));

	/* species */
	if (!IS_ASSIGNED(species) || !stricmp(OPT_STRING(species), "human")) {
		OPT_NUMBER(maxNumChr)		= 26;
		OPT_NUMBER(maxNumAutoChr)	= 22; /* X, XY, Y, MT */
		OPT_STRING(nonAutoChrSeq)	= (char *)"XZYM";
	} else {
		if (!stricmp(OPT_STRING(species), "rice")) {
			OPT_NUMBER(maxNumChr)		= 14;
			OPT_NUMBER(maxNumAutoChr)	= 12; /* MT, CH */
			OPT_STRING(nonAutoChrSeq)	= (char *)"MC";
		} if (!stricmp(OPT_STRING(species), "mouse")) {
			OPT_NUMBER(maxNumChr)		= 22;
			OPT_NUMBER(maxNumAutoChr)	= 19; /* X, Y, MT */
			OPT_STRING(nonAutoChrSeq)	= (char *)"XYM";
		} else if (!stricmp(OPT_STRING(species), "rat")) {
			OPT_NUMBER(maxNumChr)		= 23;
			OPT_NUMBER(maxNumAutoChr)	= 20; /* X, Y, MT */
			OPT_STRING(nonAutoChrSeq)	= (char *)"XYM";
		} else if (!stricmp(OPT_STRING(species), "rabbit")) {
			OPT_NUMBER(maxNumChr)		= 24;
			OPT_NUMBER(maxNumAutoChr)	= 21; /* X, Y, MT */
			OPT_STRING(nonAutoChrSeq)	= (char *)"XYM";
		} else if (!stricmp(OPT_STRING(species), "sheep")) {
			OPT_NUMBER(maxNumChr)		= 29;
			OPT_NUMBER(maxNumAutoChr)	= 26; /* X, Y, MT */
			OPT_STRING(nonAutoChrSeq)	= (char *)"XYM";
		} else if (!stricmp(OPT_STRING(species), "cow")) {
			OPT_NUMBER(maxNumChr)		= 32;
			OPT_NUMBER(maxNumAutoChr)	= 29; /* X, Y, MT */
			OPT_STRING(nonAutoChrSeq)	= (char *)"XYM";
		} else if (!stricmp(OPT_STRING(species), "horse")) {
			OPT_NUMBER(maxNumChr)		= 34;
			OPT_NUMBER(maxNumAutoChr)	= 31; /* X, Y, MT */
			OPT_STRING(nonAutoChrSeq)	= (char *)"XYM";
		} else if (!stricmp(OPT_STRING(species), "dog")) {
			OPT_NUMBER(maxNumChr)		= 41;
			OPT_NUMBER(maxNumAutoChr)	= 38; /* X, Y, MT */
			OPT_STRING(nonAutoChrSeq)	= (char *)"XYM";
		} else if (!stricmp(OPT_STRING(species), "pig")) {
			OPT_NUMBER(maxNumChr)		= 21;
			OPT_NUMBER(maxNumAutoChr)	= 18; /* X, Y, MT */
			OPT_STRING(nonAutoChrSeq)	= (char *)"XYM";
		}
	}
	
#if TOOLSET_TYPE == TOOLSET_ONETOOL
	/* Check properly initialized */
	int N_actThread = WORKER().init(OPT_NUMBER(thread));
	if (N_actThread != OPT_NUMBER(thread)) {
		char S_thread[8] = {0, };
		sprintf(S_thread, "%d", N_actThread);
		OPTION().assign("thread", S_thread, 1);
		OPTION().FORCE_OPT_NUMBER(thread);
	}
#endif

	return 0;
}

wsUint cOption::findOption(const char *S_longName, char B_noError/*=0*/)
{
	for (wsUint i=0 ; i<L_defOpt ; i++) {
		if (!stricmp(S_longName, Xa_defOption[i].S_longName + 2) ||
			(Xa_defOption[i].S_longSynonym && !stricmp(S_longName, Xa_defOption[i].S_longSynonym + 2)))
			return i;
	}

	/* Return -1 if can't find such option */
	if (B_noError) return 0xffffffff;
	halt_fmt(WISARD_INVL_OPTNAME, S_longName);
//	halt("Can't find such option `%s`", S_longName);
	return 0;
}

char cOption::isAssigned(const char *S_longName)
{
	wsUint N_idx = findOption(S_longName, 1);
  
	if (N_idx == 0xffffffff) return 0;
	return Xa_defOption[N_idx].B_isAssigned;
}

void cOption::setUnassigned(const char *S_longName)
{
	wsUint N_idx = findOption(S_longName, 1);

	if (N_idx != 0xffffffff)
		Xa_defOption[N_idx].B_isAssigned = 0;
}

void cOption::assertOption(const char *S_longName)
{
	for (wsUint i=0 ; i<L_defOpt ; i++) {
		if (!strcmp(S_longName, Xa_defOption[i].S_longName+2)) {
			if (Xa_defOption[i].B_isAssigned == false)
				halt_fmt(WISARD_NULL_ESSENTIALOPT, S_longName);
			else
				return; /* Option is normally assigned */
		}
	}

	halt_fmt(WISARD_SYST_INVL_NOSUCHOPT, S_longName);
}

/* Assert if both false */
void cOption::assertOptionAnd(wsStrCst S_ln1, wsStrCst S_ln2)
{
	for (wsUint i=0 ; i<L_defOpt ; i++)
		if (!strcmp(S_ln1, Xa_defOption[i].S_longName+2)) {
			if (Xa_defOption[i].B_isAssigned == false) {
				for (wsUint ii=0 ; ii<L_defOpt ; ii++)
					if (!strcmp(S_ln2, Xa_defOption[ii].S_longName+2)) {
						if (Xa_defOption[ii].B_isAssigned == false)
							halt_fmt(WISARD_NULL_ESSENTIALOPT, S_ln1, S_ln2);
						else
							return; /* Option is normally assigned */
					}

				halt_fmt(WISARD_SYST_INVL_NOSUCHOPT, S_ln2);
			} else
				return; /* Option is normally assigned */
		}

	halt_fmt(WISARD_SYST_INVL_NOSUCHOPT, S_ln1);
}

void cOption::_checkOptionIntegrity()
{
	strMap Xm_long, Xm_short;

	/* Check option duplication check */
	for (wsUint i=0 ; i<L_defOpt ; i++) {
		string		S_lKey	= Xa_defOption[i].S_longName;
		string		S_sKey	= Xa_defOption[i].S_shortName;
		strMap_it	X_lFind	= Xm_long.find(S_lKey);
		strMap_it	X_sFind	= Xm_short.find(S_sKey);

		if (X_lFind != Xm_long.end())
			halt_fmt(WISARD_SYST_DUPL_LONGOPT, S_lKey.c_str());
		if (X_sFind != Xm_short.end())
			halt_fmt(WISARD_SYST_DUPL_SHORTOPT, S_lKey.c_str());

		Xm_long.insert(make_pair(S_lKey, 1));
		Xm_short.insert(make_pair(S_sKey, 1));
	}

	Xm_long.clear();
	Xm_short.clear();
}

char* _isWisardLog(FILE *H_fp)
{
	static char S_buf[1024];
	char *Sp_tmp = NULL;
	/* Read the second line */
	Sp_tmp = fgets(S_buf, 1024, H_fp);
	Sp_tmp = fgets(S_buf, 1024, H_fp);

	/* Check */
	if (!memcmp("WISARD", S_buf+1, 6)) {
		LOG("WISARD log file given, options will be inherited\n");
		char *S_ptr = NULL;
		do {
			Sp_tmp = fgets(S_buf, 1024, H_fp);
			if (Sp_tmp == NULL) break;
			S_ptr = strstr(S_buf, "as : ");
		} while (!S_ptr && !feof(H_fp));

		if (!S_ptr) {
			rewind(H_fp);
			return NULL;
		}

		LOG("Following options will inherited : %s\n", S_ptr+5);
		return S_ptr+5;
	}
	rewind(H_fp);
	return NULL;
}

void cOption::_runAs(const char *S_fn)
{
	FILE *H_fp = fopen(S_fn, "r");
	char *S_buf = NULL, *s1, *s2;
	int N_mode = 0;
	if (H_fp == NULL)
		halt_fmt(WISARD_FAIL_OPEN_RUNAS, S_fn);
		//halt("Failed to open option description file `%s`", S_fn);
	wsAlloc(S_buf, char, 4096);

	/* Is it the log file of WISARD? */
	char *Sp_args = NULL;
	if ((Sp_args = _isWisardLog(H_fp)) != NULL) {
		char *s1=Sp_args, *s2=NULL;
		char *Sp_prev = NULL;

		/* Skip the first entry */
		for ( ; s1 ; s1=s2) {
			getString(&s1, &s2);
			if (N_mode == 0)
				N_mode = __procOption(s1, Sp_prev);
			else {
				if (__procValue(N_mode-1, s1)) /* N_mode(1-based) into 0-based */
					/* Try to evaluate again */
					N_mode = __procOption(s1, Sp_prev);
				else
					N_mode = 0;
			}
			Sp_prev = s1;
		}
		return;
	}

	char *Sp_prev = NULL;
	for ( ; fgets(S_buf, 4096, H_fp) ; ) {
		s1 = S_buf;

		while (s1 && s1[0]) {
			getString(&s1, &s2);

			if (N_mode == 0)
				N_mode = __procOption(s1, Sp_prev);
			else {
				if (__procValue(N_mode-1, s1)) /* N_mode(1-based) into 0-based */
					N_mode = __procOption(s1, Sp_prev);
				else
					N_mode = 0;
			}
			Sp_prev = s1;
			s1 = s2;
		}
	}

	DEALLOC(S_buf);
	fclose(H_fp);
}

double getReal(wsStrCst S_v)
{
	char *p = NULL;
	double R_ret = str2dbl(S_v, &p);
	double R_mul = 1;
	if (!p || !p[0]) return R_ret;

	if (p[1]) return WISARD_NAN;
	if (p[0] == 'k' || p[0] == 'K') R_mul = 1000;
	else if (p[0] == 'm' || p[0] == 'M') R_mul = 1000000;
	else if (p[0] == 'm' || p[0] == 'M') R_mul = 1000000000;
	else return WISARD_NAN;

	return R_ret*R_mul;
}

const char* cOption::getDefVal( const char *S_longName )
{
	wsUint N_idx = findOption(S_longName);

	return Xa_defOption[N_idx].S_defValue;
}

int cOption::procRange(xOption *Xp_opt, const char *S_val)
{
	char*	S_valDup	= NULL;
	wsUint	N_range		= 0;
	char*	Sp_val2;

	wsAlloc(S_valDup, char, strlen(S_val)+1);
	strcpy(S_valDup, S_val);
	char*	Sp_delim	= strchr(S_valDup, ',');

	/* Search bogus backslashes */
	for (char* Sp=S_valDup ; *Sp ; Sp++)
		if (*Sp == '\\')
			halt("Bogus backslash found in value [%s] of option [%s]!",
				S_val, Xp_opt->S_longName);

	if (Sp_delim == NULL) {
		/* Greater/Lower */
		if (S_valDup[0] == '>' || S_valDup[0] == '<') {
			char B_eq = S_valDup[1] == '=';
			/* Greater have lower bound, but no upper bound */
			if (S_valDup[0] == '>') {
				Xp_opt->X_rngVal.R_s = (wsReal)getReal(S_valDup+1+B_eq);
				if (NA(Xp_opt->X_rngVal.R_s))
					halt("Invalid range value [%s] of option [%s]", 
						S_val, Xp_opt->S_longName);
				Xp_opt->X_rngVal.R_e = numeric_limits<wsReal>::infinity();

				Xp_opt->X_rngVal.R_sEQ = B_eq;
				Xp_opt->X_rngVal.R_eEQ = 0;
			} else {
				Xp_opt->X_rngVal.R_s = numeric_limits<wsReal>::infinity();
				Xp_opt->X_rngVal.R_e = (wsReal)getReal(S_valDup+1+B_eq);
				if (NA(Xp_opt->X_rngVal.R_e))
					halt("Invalid range value [%s] of option [%s]", 
						S_val, Xp_opt->S_longName);

				Xp_opt->X_rngVal.R_sEQ = 0;
				Xp_opt->X_rngVal.R_eEQ = B_eq;
			}

			/* Set return value to 'have range' */
			N_range = 2;
		} else {
//			char*	Sp		= NULL;
			wsReal	R_tmp	= (wsReal)getReal(S_valDup);
			if (NA(R_tmp))
				halt("Invalid range value [%s] of option [%s]", 
					S_val, Xp_opt->S_longName);
			/* Both are same */
			Xp_opt->X_rngVal.R_s = R_tmp;
			Xp_opt->X_rngVal.R_e = Xp_opt->X_rngVal.R_s;
			Xp_opt->X_rngVal.R_sEQ = 1;
			Xp_opt->X_rngVal.R_eEQ = 1;

			/* Set return value to 'have point' */
			N_range = 1;
		}

		/* Undeterministic? */
		//halt("Invalid RANGE value `%s` for option `%s`", S_val, Xp_opt->S_longName);
	} else {
		*Sp_delim = '\0';
		Sp_val2 = Sp_delim+1;
		size_t N_pos = strlen(Sp_val2)-1;

		/* If not valid range expression */
		if ((S_valDup[0] != '[' && S_valDup[0] != '(') ||
			(Sp_val2[N_pos] != ']' && Sp_val2[N_pos] != ')'))
			goto _rngRet;

		/* Set limit */
		Xp_opt->X_rngVal.R_sEQ = S_valDup[0] == '[';
		Xp_opt->X_rngVal.R_eEQ = Sp_val2[N_pos] == ']';
		Sp_val2[N_pos] = '\0';

		Xp_opt->X_rngVal.R_s = (wsReal)getReal(S_valDup+1);
		Xp_opt->X_rngVal.R_e = (wsReal)getReal(Sp_val2);
		if (NA(Xp_opt->X_rngVal.R_s) || NA(Xp_opt->X_rngVal.R_e))
			halt("Invalid range value [%s] of option [%s]", 
				S_val, Xp_opt->S_longName);

		/* Set return value to 'have range' */
		N_range = 2;
	}
_rngRet:
	Sp_delim = strchr(S_valDup, '~');
	if (Sp_delim != NULL) {
		*Sp_delim = '\0';
		Xp_opt->X_rngVal.R_s = (wsReal)getReal(S_valDup);
		Xp_opt->X_rngVal.R_e = (wsReal)getReal(Sp_delim+1);
		if (NA(Xp_opt->X_rngVal.R_s) || NA(Xp_opt->X_rngVal.R_e))
			halt("Invalid range value [%s] of option [%s]", 
				S_val, Xp_opt->S_longName);
		Xp_opt->X_rngVal.R_sEQ = Xp_opt->X_rngVal.R_eEQ = 1;

		/* Set return value to 'have range' */
		N_range = 2;
	}

	if (N_range == 0)
		halt("Invalid range definition [%s], maybe bogus character or invalid expression?",
			S_val);

// 	LOGnote("Option [%s] value [%g]%s val %s[%g]\n", Xp_opt->S_longName,
// 		Xp_opt->X_rngVal.R_s, Xp_opt->X_rngVal.R_sEQ?"<=":"<",
// 		Xp_opt->X_rngVal.R_eEQ?"<=":"<", Xp_opt->X_rngVal.R_e);

	DEALLOC(S_valDup);
	return N_range;
}

char isInRange(xOptRange& X_range, wsReal R_val)
{
	/* Have range */
	if (X_range.R_s != X_range.R_e) {
		/* Greater than R_s */
		if (X_range.R_e == numeric_limits<wsReal>::infinity()) {
			return X_range.R_sEQ ? X_range.R_s <= R_val : X_range.R_s < R_val;
		} else if (X_range.R_s == numeric_limits<wsReal>::infinity()) {
			return X_range.R_eEQ ? R_val <= X_range.R_e : R_val < X_range.R_e;
		} else {
			return X_range.R_sEQ ?
				(X_range.R_eEQ ?
					X_range.R_s <= R_val && R_val <= X_range.R_e
					: X_range.R_s <= R_val && R_val < X_range.R_e
				)
				: (X_range.R_eEQ ?
					X_range.R_s < R_val && R_val <= X_range.R_e
					: X_range.R_s < R_val && R_val < X_range.R_e
				);
		}
	} else
		return X_range.R_s == R_val;
}

char isInRangeOR(xOptRange& X_range, wsReal* Ra_val, wsUint N_val)
{
	bool B_ret = false;
	/* Have range */
	if (X_range.R_s != X_range.R_e) {
		/* Greater than R_s */
		if (X_range.R_e == numeric_limits<wsReal>::infinity()) {
			for (wsUint i=0 ; i<N_val ; i++)
				B_ret |= X_range.R_sEQ ? X_range.R_s <= Ra_val[i] : X_range.R_s < Ra_val[i];
		} else if (X_range.R_s == numeric_limits<wsReal>::infinity()) {
			for (wsUint i=0 ; i<N_val ; i++)
				B_ret |= X_range.R_eEQ ? Ra_val[i] <= X_range.R_e : Ra_val[i] < X_range.R_e;
		} else {
			for (wsUint i=0 ; i<N_val ; i++)
				B_ret |= X_range.R_sEQ ?
					(X_range.R_eEQ ?
					X_range.R_s <= Ra_val[i] && Ra_val[i] <= X_range.R_e
					: X_range.R_s <= Ra_val[i] && Ra_val[i] < X_range.R_e
					)
					: (X_range.R_eEQ ?
					X_range.R_s < Ra_val[i] && Ra_val[i] <= X_range.R_e
					: X_range.R_s < Ra_val[i] && Ra_val[i] < X_range.R_e
					);
		}
	} else for (wsUint i=0 ; i<N_val ; i++)
		B_ret |= X_range.R_s == Ra_val[i];
	return B_ret;
}

char isInRangeAND(xOptRange& X_range, wsReal* Ra_val, wsUint N_val)
{
	bool B_ret = true;
	/* Have range */
	if (X_range.R_s != X_range.R_e) {
		/* Greater than R_s */
		if (X_range.R_e == numeric_limits<wsReal>::infinity()) {
			for (wsUint i=0 ; i<N_val ; i++)
				B_ret &= X_range.R_sEQ ? X_range.R_s <= Ra_val[i] : X_range.R_s < Ra_val[i];
		} else if (X_range.R_s == numeric_limits<wsReal>::infinity()) {
			for (wsUint i=0 ; i<N_val ; i++)
				B_ret &= X_range.R_eEQ ? Ra_val[i] <= X_range.R_e : Ra_val[i] < X_range.R_e;
		} else {
			for (wsUint i=0 ; i<N_val ; i++)
				B_ret &= X_range.R_sEQ ?
					(X_range.R_eEQ ?
					X_range.R_s <= Ra_val[i] && Ra_val[i] <= X_range.R_e
					: X_range.R_s <= Ra_val[i] && Ra_val[i] < X_range.R_e
					)
					: (X_range.R_eEQ ?
					X_range.R_s < Ra_val[i] && Ra_val[i] <= X_range.R_e
					: X_range.R_s < Ra_val[i] && Ra_val[i] < X_range.R_e
					);
		}
	} else for (wsUint i=0 ; i<N_val ; i++)
		B_ret &= X_range.R_s == Ra_val[i];
	return B_ret;
}

} // End namespace ONETOOL
