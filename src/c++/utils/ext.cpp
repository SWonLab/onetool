#ifdef WIN32
#	include <windows.h>
#else
#	include <dlfcn.h>
#endif
#include "global/common.h"
#include "utils/util.h"
#include "utils/ext.h"

namespace ONETOOL {

void doWISARDext(wsStrCst S_fnExt)
{
	WISARD_extIns	H_inst = loadExtLib(S_fnExt);

	if (!isValidExtLibHandle(H_inst))
		halt("External load failed [%s]", S_fnExt);

	/* Do something */

	unloadExtLib(H_inst);
}

} // End namespace ONETOOL
