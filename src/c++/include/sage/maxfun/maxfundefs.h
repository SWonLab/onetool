#ifndef __MAXFUN_DEFINITION_H
#define __MAXFUN_DEFINITION_H

//==========================================================================
// File:      def_maxfun.h
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                            yes Jun 2015
//
// Notes:     This header file provides globally used typedef, define ...
//            in MAXFUN & include all needed header files.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "sage/global/error.h"
#include "sage/global/get_mem.h"
#include "sage/global/functions.h"
#include "sage/global/log_double.h"
#include "sage/global/aparser.h"

namespace SAGE   {
namespace MAXFUN {

const double MF_INFINITY = std::numeric_limits<double>::infinity();

class SequenceCfg;
class Datatypes;
class DebugCfg;
class ParameterMgr;  
class APIMaxFunction;
class Maximizer;
class Parameter;
class OutputFormatter; 
class LogitTransformer;
class Transformer;
class ParameterIterator;
class ParameterConstIterator;
class ParamCalculator;
class Results;    
class Submodel;   
class NewSubmodel;

  
} // End namesapce MAXFUN
} // End namesapce SAGE

#endif  
