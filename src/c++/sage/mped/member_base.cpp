//============================================================================
//  File:       mpbase.cpp
//
//  Author:     Bob Steagall
//
//  History:    Version 0.98
//
//  Copyright (c) 1998 R.C. Elston
//  All Rights Reserved
//============================================================================
//
#include "sage/mped/spbase.h"

namespace SAGE {
namespace MPED {

//============================================================================
//  IMPLEMENTATION: member_base
//============================================================================
//
member_base::member_base(const string& mname, SexCode G)
  : my_mpindex((uint) -1),
    my_index((uint) -1),
    my_subindex((uint) -1),
    my_name(mname), 
    my_subped(0),
    my_origin(0), 
    my_siblings(0), 
    my_gender(G), 
    my_offspring_count(0), 
    my_mate_count(0)
{}

}
}
