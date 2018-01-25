//============================================================================
// File:      analysis.ipp
//                                                                          
// Author:    Dan Baechle
//                                                                          
// History:   12/20/2 created        -djb
//                                                                          
// Notes:     Inline implementation of analysis.
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

//============================================================================
// IMPLEMENTATION:  analysis
//============================================================================
//
inline
analysis::analysis(const instructions& instr,
                   const RPED::RefMultiPedigree& mped,
                   cerrorstream& errors)
      : my_original_mped(mped), my_errors(errors), my_mped(mped),
        aborted(false)
{
  my_instructions = instr;
  build_filtered_mped();
}

// - Deletion of tasks not necessary because boost smart_ptr's are 
//   used.
//
inline
analysis::~analysis()
{}
