#ifndef __SAGE_DEFINITION_H
#define __SAGE_DEFINITION_H

//==========================================================================
// File:      definition.h
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                            yes Dec 2014
//
// Notes:     This header file provides globally used typedef, define ...
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "sage/global/lib_include.h"
#include "sage/global/clapack.h"
#include "sage/global/kahan.h"

using namespace std;

namespace SAGE {

using std::isnan;

const double QNAN = std::numeric_limits<double>::quiet_NaN();
const size_t NPOS = ((size_t) -1) / 2;

const double POSITIVE_INF        = std::numeric_limits<double>::infinity();
const double NEGATIVE_INF        = -POSITIVE_INF;
const double DEFAULT_LOWER_BOUND = NEGATIVE_INF;
const double DEFAULT_UPPER_BOUND = POSITIVE_INF;

template <class T>
class numeric_constants : public std::numeric_limits<T>
{};

template<>
class numeric_constants<double> : public std::numeric_limits<double>
{
  public:
    static double maxlog() { return (max_exponent-1.0)*log(2.0);    }
    static double minlog() { return (min_exponent-digits)*log(2.0); }
    static double pi()     { return 3.14159265358979323846;         }
};


#ifdef __WIN32__
const std::string path_delimiter = "\\";
#else
const std::string path_delimiter = "/";
#endif

#if defined(MINGW) || _WIN32
typedef unsigned int    uint;
typedef unsigned short  ushort;
#endif

typedef unsigned char   uchar;

#if defined(__MACOS__)
typedef unsigned int    uint;
typedef unsigned short  ushort;
#endif

// __U32 is an unsigned 32 bit integer.
#define __U32 unsigned int
#define __WORD_BIT (int(CHAR_BIT*sizeof(__U32)))

} // End namesapce SAGE

#endif
