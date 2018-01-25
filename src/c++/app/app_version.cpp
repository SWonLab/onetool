//==========================================================================
// File:      app_version.cpp
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                              yes Nov 14
//
// Notes:     This source file implements app_version.h.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "app/app_version.h"

#ifdef _WIN32
#	define ONETOOL_MAIN_VERSION 2
#	define ONETOOL_SUB_VERSION 0
#	define ONETOOL_MICRO_VERSION 0
#	define ONETOOL_BETA_VERSION 0
#else
#if !defined(ONETOOL_MAIN_VERSION) || !defined(ONETOOL_SUB_VERSION) || !defined(ONETOOL_MICRO_VERSION) || !defined(ONETOOL_BETA_VERSION)
#error You must declare all four components of this version for app to compile correctly (ie: gcc -DONETOOL_MAIN_VERSION=1 -DONETOOL_SUB_VERSION=2 -DONETOOL_MICRO_VERSION=3 -DONETOOL_BETA_VERSION=0)
#endif
#endif

namespace ONETOOL {
namespace APP    {

app_version version(ONETOOL_MAIN_VERSION, ONETOOL_SUB_VERSION, ONETOOL_MICRO_VERSION, ONETOOL_BETA_VERSION);

} // End namespace APP
} // End namespace ONETOOL
