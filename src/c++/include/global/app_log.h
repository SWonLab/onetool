#ifndef __APP_LOG_H
#define __APP_LOG_H

//============================================================================
// File:      app_log_infoqc.h
//
// Author:    S.A.G.E. crew
//            Yeunjoo E. Song
//
// History:   Initial implementation in S.A.G.E.
//            Modified for INFOQC                                 yes Jan 2015
//
// Notes:     This file defines the global error function
//            - app_log.ipp : inline functions
//            - app_log.cpp :
//
// Copyright (c) 2015 Sungho Won
// All Rights Reserved                                                    
//============================================================================

#include "utils/util.h"
#include "sage/global/output_streams.h"

namespace ONETOOL {

typedef SAGE::Output_Streams app_log;


string to_upper(const string &str);
string to_lower(const string &str);
string strip_ws(const string& s, const char *ws=NULL);
} // End INFOQC namespace

#endif
