//============================================================================
// File:      analysis.cpp                      
//                                                                          
// Author:    Dan Baechle                                     
//                                                                          
// History:   12/20/2 - created.                         djb
//                                                                          
// Notes:     Implementation of analysis class.   
//               
//                                                                          
// Copyright (c) 2002 R.C. Elston                                           
// All Rights Reserved                                                    
//============================================================================

#include "sage/lodlink/analysis.h"

namespace SAGE
{

namespace LODLINK
{

void
build_headers()
{
  cleves_elston_result::build_headers();
  non_ss_lod_ratio_result::build_headers();
  ss_lod_ratio_result::build_headers();
  non_ss_smiths_result::build_headers();
  ss_smiths_result::build_headers();
  non_ss_faraways_result::build_headers();
  ss_faraways_result::build_headers();
  non_ss_mortons_result::build_headers();
  ss_mortons_result::build_headers();
  non_ss_genotype_result::build_headers();
  ss_genotype_result::build_headers();
}

//============================================================================
// IMPLEMENTATION:  analysis
//============================================================================
//
void
analysis::run_analysis(ostream& sum, ostream& det)
{
  build_headers();
  build();
  analyze();
  write(sum, det);

  return;
}

// - Create tasks per user instructions in the order in which output will be written.
//
void
analysis::build()
{
  if(my_instructions.valid)
  {
    // - Lod scores.
    //
    if(my_instructions.average_thetas.size())
    {
      my_tasks.push_back(task_ptr(new non_ss_lods(my_errors, my_mped, my_instructions)));
    }
    
    if(my_instructions.male_female_thetas.size())
    {
      my_tasks.push_back(task_ptr(new ss_lods(my_errors, my_mped, my_instructions)));
    }
    
    // - Linkage tests.
    //
    if(my_instructions.linkage_test)
    {
      if(my_instructions.linkage_sex_specific)
      {
        if(my_instructions.linkage_homog)
        {
          my_tasks.push_back(task_ptr(new ss_lod_ratio_test(my_errors, my_mped, my_instructions)));
          my_tasks.push_back(task_ptr(new cleves_elston_test(my_errors, my_mped, my_instructions)));
        }
        else
        {
          my_tasks.push_back(task_ptr(new ss_smiths_faraways_test(my_errors, my_mped, my_instructions, sf_LINKAGE)));  
        }
      }
      else
      {
        if(my_instructions.linkage_homog)
        {
          my_tasks.push_back(task_ptr(new non_ss_lod_ratio_test(my_errors, my_mped, my_instructions)));
        }
        else
        {
          my_tasks.push_back(task_ptr(new non_ss_smiths_faraways_test(my_errors, my_mped, my_instructions, sf_LINKAGE)));
        }
      }
    }
    
    // - Linkage homogeneity tests.
    //
    if(my_instructions.smiths_test)
    {
      if(my_instructions.smiths_sex_specific)
      {
        my_tasks.push_back(task_ptr(new ss_smiths_faraways_test(my_errors, my_mped, my_instructions, sf_HOMOGENEITY)));
      }
      else
      {
        my_tasks.push_back(task_ptr(new non_ss_smiths_faraways_test(my_errors, my_mped, my_instructions, sf_HOMOGENEITY)));
      }
    }
    
    if(my_instructions.mortons_test)
    {
      // - groups of subpedigrees must be mutually exclusive and comprehenive.
      //                                                        -RCE 2-3-3   
      //                                               
      // - No groups means that user tried unsuccessfully to specify a set of 
      //   groups that meet the above criteria.  Skip Morton's test.
      //
      if(! my_instructions.groups.empty())
      {
        if(my_instructions.mortons_sex_specific)
        {
          my_tasks.push_back(task_ptr(new ss_mortons_test(my_errors, my_mped, my_instructions)));
        }
        else
        {
          my_tasks.push_back(task_ptr(new non_ss_mortons_test(my_errors, my_mped, my_instructions)));
        }
      }
    }
      
    if(my_instructions.genotypes == true)
    {
      if(my_instructions.genotypes_sex_specific == true)
      {
        my_tasks.push_back(task_ptr(new ss_genotype_probs(my_errors, my_mped, my_instructions)));
      }
      else
      {
        my_tasks.push_back(task_ptr(new non_ss_genotype_probs(my_errors, my_mped, my_instructions)));
      }
    }
  }
}

// - Execute tasks.
//
void
analysis::analyze()
{
  if(my_instructions.sex_specific())
  {
    if(missing_sex_parents())
    {
      my_errors << "Sex specific analysis specified with missing information in sex field of some "
                << "parents.  Please supply the missing information and run again.  Skipping analysis ..." << endl;
      aborted = true;
      
      return;
    }
  }
  else
  {
    if(missing_sex_parents())
    {
      assign_arbitrary_sexes();
    }
  }
  
  for(size_t i = 0; i < my_tasks.size(); ++i)
  {
    my_tasks[i]->announce_start();
    my_tasks[i]->calculate();
  }
}

// - Write results of tasks.
//
void
//analysis::write(ostream& summary, ostream& detail, APP::SAGEapp& app)
analysis::write(ostream& summary, ostream& detail)
{
  if(! aborted)
  {
    //write_file_header(summary, "Summary", app);
    write_file_header(summary, "Summary");
    my_instructions.write(summary);
    analysis::write_results_label(summary);
    
    //write_file_header(detail, "Detail", app);
    write_file_header(summary, "Summary");
    my_instructions.write(detail);
    analysis::write_results_label(detail);

    for(size_t i = 0; i < my_tasks.size(); ++i)
    {
      my_tasks[i]->write(summary, detail);
    }
  }
}

void
analysis::write_results_label(ostream& out)
{
  const string  RESULTS = "Results";
  
  ios::fmtflags old_flags = out.flags();
  
  out << setfill('=') << left
      << RESULTS << "\n"
      << setw(RESULTS.size()) << "" << "\n" << endl;
  
  out.flags(old_flags);
}

void
//analysis::write_file_header(ostream& out, const string& file_type, APP::SAGEapp& app)
analysis::write_file_header(ostream& out, const string& file_type)
{
  //app.print_title(out);

  ios::fmtflags old_flags = out.flags();

  string  title = my_instructions.title + " " + toUpper(file_type) + " FILE";
  out << left << setfill('=')
      << setw(title.size()) << "" << "\n"
      << title << "\n"
      << setw(title.size()) << "" << "\n\n" << endl;
  
  out.flags(old_flags); 
}

// - Clear static elements from Smith's and Faraway's task classes.
//
void
analysis::clear()
{
  non_ss_smiths_faraways_test::clear();
  ss_smiths_faraways_test::clear();
  ge_models::clear_models();
}

void
analysis::build_filtered_mped()
{
  FPED::MPFilterer::add_multipedigree_filtered_by_members(my_mped, my_original_mped, FPED::always_keep());
  my_mped.construct();

#if 0
  cout << endl;
  cout << "original_multipedigree :" << endl;
  cout << "                  info :" << endl;
  const RPED::RefMPedInfo& m_info = my_original_mped.info();
  cout << "  trait count  = " << m_info.trait_count() << endl;
  cout << "  marker count = " << m_info.marker_count() << endl;

  for( size_t mi = 0; mi < m_info.marker_count(); ++mi )
  {
    cout << "          " << mi << " : "
         << m_info.marker_info(mi).name() << "  "
         << m_info.marker_info(mi).allele_count() << "  "
         << m_info.marker_info(mi).phased_genotype_count() << " "
         << m_info.marker_info(mi).unphased_genotype_count() << "       "
         << m_info.marker_info(mi).phenotype_count() << endl;
  }

  cout << "Member Info : " << endl;
  cout << "  member count = " << my_original_mped.member_count() << endl;

  for( size_t i = 0; i < my_original_mped.member_count(); ++i )
  {
    const RPED::member_type& mem = my_original_mped.member_index(i);

    const RPED::RefPedInfo& p_info = mem.pedigree()->info();

    cout << "    " << i << "  " << mem.mpindex() << "  " << mem.index()
         << "," << mem.name() << " : ";
    for( size_t m = 0; m < p_info.marker_count(); ++m )
      cout << p_info.phenotype(mem.index(), m) << "    ";
    cout << endl;
  }

  cout << endl;
  cout << "filtered_multipedigree :" << endl;
  cout << "                  info :" << endl;
  const RPED::RefMPedInfo& f_info = my_mped.info();
  cout << "  trait count  = " << f_info.trait_count() << endl;
  cout << "  marker count = " << f_info.marker_count() << endl;

  for( size_t mi = 0; mi < f_info.marker_count(); ++mi )
  {
    cout << "          " << mi << " : "
         << f_info.marker_info(mi).name() << "  "
         << f_info.marker_info(mi).allele_count() << "  "
         << f_info.marker_info(mi).phased_genotype_count() << " "
         << f_info.marker_info(mi).unphased_genotype_count() << "       "
         << f_info.marker_info(mi).phenotype_count() << endl;
  }

  cout << "Member Info : " << endl;
  cout << "  member count = " << my_mped.member_count() << endl;

  for( size_t i = 0; i < my_mped.member_count(); ++i )
  {
    const FPED::Member& mem = my_mped.member_index(i);

    const FPED::FilteredPedigreeInfo& p_info = mem.pedigree()->info();

    cout << "    " << i << "  " << mem.mpindex() << "  " << mem.index()
         << "," << mem.name() << " : ";
    for( size_t m = 0; m < p_info.marker_count(); ++m )
      cout << p_info.phenotype(mem.index(), m) << "    ";
    cout << endl;
  }
#endif
}

// - Are there any parents in my data who are missing sex information?
//
bool
analysis::missing_sex_parents() const
{
  size_t  member_count = my_mped.member_count();
  for(size_t m = 0; m < member_count; ++m)
  {
    const FPED::Member&  member = my_mped.member_index(m);
    if(member.mate_count() && member.is_sex_unknown())
    {
      return  true;
    }
  }
  
  return  false;
}

void
analysis::assign_arbitrary_sexes()
{
  size_t  member_count = my_mped.member_count();
  for(size_t m = 0; m < member_count; ++m)
  {
    FPED::Member&  member = my_mped.member_index(m);
    if(member.is_sex_unknown())
    {
      member.set_sex(MPED::SEX_XMALE);  // SEX_XMALE means arbitrarily assigned as male
      member.pedigree()->build();
    }
  }
}


}
}
