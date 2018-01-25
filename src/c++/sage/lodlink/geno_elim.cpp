//============================================================================
// File:      geno_elim.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   3/14/3 - created.                         djb
//                                                                          
// Notes:     Implementation a class for creating and storing genotype-
//            eliminated penetrance models.   
//               
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "sage/lodlink/geno_elim.h"

namespace SAGE    {
namespace LODLINK {

//============================================================================
// IMPLEMENTATION:  ge_models
//============================================================================
//
pedigree_imodel_generator  ge_models::my_generator;
ge_models::model_map       ge_models::my_models;

/// Clears the current set of models from storage, preventing accidentally looking
/// up the wrong thing between analyses.
void ge_models::clear_models()
{
  my_models.clear();
}

const MLOCUS::penetrance_model&
ge_models::get_model(const subpedigree& subped, size_t locus_idx)
{
  //my_generator.set_genotype_elimination(false);
  model_map::const_iterator  iter = my_models.find(make_pair(&subped, locus_idx));
  if(iter != my_models.end())
  {
    return  iter->second;
  }
  else
  {

#if 0
  cout << "original subpedigree :" << endl;
  cout << " name = " << subped.name() << ", mem_count = " << subped.member_count() << endl;

  const RPED::RefMPedInfo& m_info = subped.multipedigree()->info();
  cout << "  trait count  = " << m_info.trait_count() << endl;
  cout << "  marker count = " << m_info.marker_count() << endl;

  cout << "Family Info : " << endl;
  FPED::FilteredMultipedigree::family_const_iterator fi = subped.family_begin();
  for( ; fi != subped.family_end(); ++fi )
  {
    cout << fi->index() << ": " << fi->name() << endl;
    cout << "p1 = " << fi->parent1()->name() << endl;
    cout << "p2 = " << fi->parent2()->name() << endl;

    FPED::FilteredMultipedigree::offspring_const_iterator oi = fi->offspring_begin();
    for( ; oi != fi->offspring_end(); ++oi )
      cout << "offspring " << oi->index() << " : " << oi->name() << endl;
  }

  cout << "Member Info : " << endl;
  //const RPED::RefPedInfo& p_info = subped.pedigree()->info();
  const FPED::FilteredPedigreeInfo& p_info = subped.pedigree()->info();
  cout << "  marker count = " << p_info.marker_count() << endl;
  cout << "  member count = " << p_info.member_count() << endl;
  for( size_t i = 0; i < p_info.member_count(); ++i )
  {
    cout << "    " << i << "  " << subped.member_index(i).index()
         << "," << subped.member_index(i).name() << " : ";
    for( size_t m = 0; m < p_info.marker_count(); ++m )
      cout << p_info.phenotype(i, m) << "    ";
    cout << endl;
  }
#endif

    my_models.insert(make_pair(make_pair(&subped, locus_idx), my_generator(subped, locus_idx)));
    return  get_model(subped, locus_idx);
  }
}

} // end of namespace LODLINK
} // end of namespace SAGE



