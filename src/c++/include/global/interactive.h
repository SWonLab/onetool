#pragma once
#ifndef __WISARD_INTERACTIVE_H__
#define __WISARD_INTERACTIVE_H__

namespace ONETOOL {

typedef	vector<xParam>		vParam;
typedef	vParam::iterator	vParam_it;

typedef enum _xItrtCmdKey
{
	CMD_EXPL_QUERYDATA_BY_SNPNAME,
	CMD_EXPL_QUERYDATA_BY_SNPIDX,
	CMD_EXPL_QUERYDATA_BY_SNPNAME_AVAIL,
	CMD_EXPL_QUERYDATA_BY_SNPIDX_AVAIL,
	CMD_EXPL_QUERYDATA_BY_SAMPIID,
	CMD_EXPL_QUERYDATA_BY_SAMPIDX,
	CMD_EXPL_QUERYDATA_BY_SAMPIID_AVAIL,
	CMD_EXPL_QUERYDATA_BY_SAMPIDX_AVAIL,
	CMD_EXPL_QUERY_SNPNAME,
	CMD_EXPL_QUERY_SNPNAME2,
	CMD_EXPL_QUERY_SAMPIID,
	CMD_EXPL_QUERY_SAMPIID2,
	CMD_EXPL_QUERY_SAMPFID,
	CMD_EXPL_QUERY_SIZE,
	CMD_EXPL_QUERY_SAMPINFO,
} xItrtCmdKey;

struct _xItrtCmdDef;
typedef int (*hItrtFunc)(struct _xItrtCmdDef &X_cmd, vParam &V_param, void *Vp_data);

typedef struct _xItrtCmdDef
{
	char		S_cmd[64];
	char		S_fmt[8];
	xItrtCmdKey	X_cmd;
	hItrtFunc	H_func;
	char		S_optDesc[256];
	char		S_desc[256];
} xItrtCmd;

#define ITRTFUNC(a) int a(xItrtCmd &X_cmd, vParam &V_param, void *Vp_data)

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cIO;
void interactiveExplore(cIO *Cp_IO);
void interactiveExecution();

#else

class cIO;
void interactiveExplore(cIO *Cp_IO);

#endif

} // End namespace ONETOOL

#endif
