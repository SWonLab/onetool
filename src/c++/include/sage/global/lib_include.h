#ifndef __LIB_INCLUDE_H
#define __LIB_INCLUDE_H

//==========================================================================
// File:      lib_include.h
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                            yes Dec 2014
//
// Notes:     This header file lists the #include directives for all globally
//            needed c/c++ library.
//            - pooled from all SAGE codes
//            - any new include should be added here
//            - both c++ and c style exist, but c++ style preferred
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include <algorithm>
#include <bitset>
#include <deque>
#include <cassert>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <climits>
#include <functional>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <list>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <assert.h>
#include <errno.h>
#include <math.h>
#ifndef _WIN32
#	include <pwd.h>
#endif
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#ifndef _WIN32
#   include <unistd.h>
#endif

#ifdef WIN32
#   include <windows.h>
#	pragma warning(disable:4267)
#	pragma warning(disable:4244)
#	pragma warning(disable:4996)
#	pragma warning(disable:4146)
#else
#   include <limits.h>
#endif

#endif  
