#pragma once
#ifndef __WISARD_META_H__
#define __WISARD_META_H__
#include "global/io.h"

namespace ONETOOL {

/* !!! WISARD-specific analysis !!! */
#if (TOOLSET_TYPE & 0x100) == 0x100

class cMetaIO : public cIO
{
public:
	cMetaIO(wsStrCst S_fn, char B_inDry);
	void	init(wsStrCst S_fn);
	wsStrCst	getFormat() { return "io.meta"; }
};

#else

class cMetaIO : public cIO
{
public:
	cMetaIO(wsStrCst S_fn, char B_inDry) {}
	void	init(wsStrCst S_fn) {}
	wsStrCst	getFormat() { return "io.meta"; }
};

#endif

} // End namespace ONETOOL

#endif
