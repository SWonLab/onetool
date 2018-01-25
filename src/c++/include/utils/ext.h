#pragma once
#ifndef __WISARD_EXT_H__
#define __WISARD_EXT_H__

#ifdef WIN32
#	define WISARD_extIns			HINSTANCE
#	define loadExtLib(f)			LoadLibrary(f)
#	define unloadExtLib(h)			FreeLibrary(h)
#	define loadFunc(h, n)			GetProcAddress(h, n)
#	define isValidExtLibHandle(a)	((a) == 0)
#else
#	define WISARD_extIns			void *
#	define loadExtLib(f)			dlopen(f, RTLD_LAZY)
#	define unloadExtLib(h)			dlclose(h)
#	define loadFunc(h, n)			dlsym(h, n)
#	define isValidExtLibHandle(a)	((a) == 0)
#endif

namespace ONETOOL {

void doWISARDext(wsStrCst S_fnExt);

} // End namespace ONETOOL

#endif
