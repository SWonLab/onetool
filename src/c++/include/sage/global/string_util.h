#ifndef __SAGE_UTIL_H
#define __SAGE_UTIL_H

//==========================================================================
// File:      string_util.h
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                            yes Dec 2014
//
// Notes:     This header file provides globally used utility functions
//            for string.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "sage/global/error.h"

namespace SAGE {

//
// string utility
//
//int    kill_ws(FILE *i, const char *ws = " \t", const char *eol = "\r\n");
//string getString(FILE *i, const char *delim = "", const char *ws = " \t\n\r" );

int    kill_ws(   istream &i, const char *ws = " \t", const char *eol = "\r\n");
string getStringSAGE( istream &i, const char *delim = "", const char *ws = " \t\n\r" );
string getQStringSAGE(istream &i, const char *delim = "", const char *ws = " \t",
                              const char   *eol = "\n\r", int *lines = NULL );
int    printQuoted(ostream &o, string& s, int _pos = 0, int ind = 0, int _start = 0);

string to_upper(const string &s);
string to_lower(const string &s);

string toUpper(const string &s);
string toLower(const string &s);

string strip_ws(const string& s, const char *ws = NULL);

string doub2str(double d,      int width = 0, int pres = -1, long flags = 0, char fill = '\0');
string long2str(long   l,      int width = 0, long flags = 0, char fill = '\0');
string ptr2str (const void* v, int width = 0, long flags = 0, char fill = '\0');

double str2doub(const string &s);
long   str2long(const string &s);

//
// print utility
//
string fp(double d, size_t w, size_t p=4, char invalid='-');
string pval(double p, size_t w, int prec=-1, int max_stars=2);

string fp_scientific(double d, size_t w, size_t p=4, char invalid='-');
string pval_scientific(double p, size_t w, int prec=-1, int max_stars=2);

}


#endif  
