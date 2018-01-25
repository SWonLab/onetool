#ifndef __WISARD_MACH_H__
#define __WISARD_MACH_H__
#pragma once
#include "global/io.h"

namespace ONETOOL {

class cMachIO : public cIO
{
	char	_loadMain(char *Sp_buf, wsUint fiIdx, wsUint diIdx,
		wsReal *Rp_pheno, vBool& Ba_isSampFiltOrig, vBool& Ba_isVrtFiltOrig,
		cExporter **Cp_misInds);
public:
	cMachIO(char *S_fn, char B_inDry);
	void	init(char *S_fn);
	wsStrCst	getFormat() { return "io.mach"; }
};

} // End namespace ONETOOL

#endif
