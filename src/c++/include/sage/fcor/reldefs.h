#ifndef RELATIONDEFS_H
#define RELATIONDEFS_H

//****************************************************************************
//* File:      reldefs.h                                                     *
//*                                                                          *
//* Author:    Yeunjoo Song                                                  *
//*                                                                          *
//* History:   Version 1.0                                                   *
//*                                                                          *
//* Notes:     This header file defines various types used in computing      *
//*            relationships between individuals in a reference pedigree.    *
//*                                                                          *
//* Copyright (c) 1999 R.C. Elston                                           *
//*   All Rights Reserved                                                    *
//****************************************************************************

#include "sage/rped/rped.h"

namespace SAGE {


extern const char* ordinal_names[];
extern const int   max_ordinal_names;
extern const char* ordinal_th_names[];
extern const int   max_ordinal_th_names;
extern const char* ordinal_tens_names[];
extern const int   max_ordinal_tens_names;
extern const char* ordinal_power_names[];
extern const int   max_ordinal_power_names;
extern const char* ordinal_st_names[];
extern const int   max_ordinal_st_names;
extern const char* relative_words[];
extern const int   max_relative_words;

typedef RPED::member_type                               pedigree_member_type;
typedef RPED::member_pointer                            pedigree_member_pointer;

typedef std::pair<pedigree_member_pointer,
                  pedigree_member_pointer>              pedigree_member_pair;

typedef std::pair<const pedigree_member_type*,
                  const pedigree_member_type*>          const_pedigree_member_pair;

typedef std::list<pedigree_member_pointer>              pedigree_member_list;

typedef std::pair<pedigree_member_list,
                  pedigree_member_list>                 pedigree_member_list_pair;

} // end of namespace SAGE

#endif
