#ifndef LODLINK_ANALYSIS_H
#define LODLINK_ANALYSIS_H
//============================================================================
// File:      analysis.h
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   12/20/2 - created.                                   djb
//                                                                          
// Notes:     declaration of a analysis class. 
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================


#include "sage/lodlink/linkage_tests.h"
#include "sage/lodlink/homogeneity_tests.h"
#include "sage/lodlink/lods.h"
#include "sage/lodlink/genotypes.h"

using std::vector;
using std::ostream;

namespace SAGE
{

namespace LODLINK
{

void  build_headers();

//----------------------------------------------------------------------------
//  Class:    analysis
//                                                                          
//  Purpose:  execute and write results of a set of user instructions 
//            corresponding to a parameter file lodlink analysis block.
//                                                                          
//----------------------------------------------------------------------------
//
class analysis
{
  public:

    analysis(const instructions& instr, const RPED::RefMultiPedigree& mped, cerrorstream& errors);
    ~analysis();

    void run_analysis(ostream& summary, ostream& detail);

    static void  clear();

  private:

    void  build();
    void  analyze();
    void  write(ostream& summary, ostream& detail);

    static void  write_results_label(ostream& out);

    void  build_filtered_mped();
    void  write_file_header(ostream& out, const string& file_type);
    bool  missing_sex_parents() const;
    void  assign_arbitrary_sexes();
  
    // Data members.
    const RPED::RefMultiPedigree&  my_original_mped;
    cerrorstream&                  my_errors;
    FPED::FilteredMultipedigree    my_mped;
    instructions                   my_instructions;
    
    typedef std::shared_ptr<task>  task_ptr;
    vector<task_ptr>  my_tasks;
    
    bool  aborted;
};

#include "sage/lodlink/analysis.ipp"
}
}

#endif

