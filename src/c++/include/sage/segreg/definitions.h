#ifndef SEGREG_DEFINITIONS_H
#define SEGREG_DEFINITIONS_H
//============================================================================
// File:      definitions.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   3/8/2 created         -djb
//                                                                          
// Notes:     For definitions common to model and sub-models.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//#include "boost/smart_ptr.hpp"
#include "boost/bind.hpp"
#include "boost/array.hpp"
#include "boost/iterator/counting_iterator.hpp"
//#include "sage/global/string_tokenizer.h"
#include "sage/global/sinfo.h"
#include "sage/maxfun/Maximizer.h"
#include "sage/rped/peeler3.h"

//using namespace ONETOOL;

namespace SAGE   {
namespace SEGREG {

using std::vector;

using namespace RPED;
using namespace MAXFUN;

enum model_class { model_A, model_D, model_FPMM, model_MLM, model_INVALID };
enum primary_type { pt_NONE, pt_CONTINUOUS, pt_BINARY, pt_ONSET };

// moved from segreg_options.h
const bool  ONSET_AVAILABLE = true;
const bool  FPMM_AVAILABLE = true;
const bool  PREV_CONSTRAINT_AVAILABLE = true;
const bool  PREV_ESTIMATE_AVAILABLE = true;
const bool  TYPE_SUSCEPT_AVAILABLE = true;
const bool  INTERACTIONS_AVAILABLE = true;
const bool  TYPE_CORR_AVAILABLE = false;
const bool  EACH_PEDIGREE_AVAILABLE = false;
const bool  TYPE_PROBABILITIES_AVAILABLE = true;
const bool  PEN_FUNC_OUT_AVAILABLE = true;

// moved from sub_model_base.h
const int  DUMP_PRECISION = 12;
const int  NUM_OF_TYPES = 3;

enum gi_type { index_AA = 0, index_AB, index_BB, index_INVALID };
/*
class genotype_index
{
  public:

    inline genotype_index(unsigned int);
    inline genotype_index(gi_type = index_AA);   
    inline genotype_index(const genotype_index&);

    inline genotype_index& operator=(unsigned int);
    inline genotype_index& operator=(gi_type);
    inline genotype_index& operator=(const genotype_index&);

    inline operator unsigned int () const;

    inline genotype_index& operator ++();   
    inline genotype_index  operator ++(int);

    inline bool operator==(const genotype_index&) const;
    inline bool operator!=(const genotype_index&) const;

    inline bool operator==(gi_type) const;
    inline bool operator!=(gi_type) const;

  protected:

    gi_type my_type;
};
*/
enum genotype_info { no_geno = 0, 
                     AA = 1, AB = 2, BB = 4,
                     AA_AB = 3, AA_BB = 5, AB_BB = 6,   
                     all = 7
                   };

enum freq_index { index_freq_A = 3, index_corr };


// moved from segreg_datatypes.h
enum relationship_type  { undef = -1, no_rel = 0, sib = 1, parent_off = 2, mate = 3 };

class segreg_errors
{
  public:

    enum error_code { EVAL_OK = 0,
                      MCC_FAILED_TRANSFORM,
                      MCC_VARIANCE_INVALID,
                      BAD_INIT_PARAM,
                      ZERO_LIKELIHOOD,
                      BAD_THRESHOLDS,
                      BAD_LIKELIHOOD,
                      BAD_LEX_VALUE,
                      MLM_CORR_BOUNDARY
                    };
};

#define MAX_POLYGENOTYPE 11
#define NO_POLYGENIC_DATA 999
/*
struct genetic_info
{
  genetic_info(genotype_index genotype_param     = index_AA,
               size_t         polygenotype_param = NO_POLYGENIC_DATA)
  {
    genotype     = genotype_param;
    polygenotype = polygenotype_param;
  }

  genotype_index genotype;
  size_t         polygenotype;
};
*/

}
}

#endif

