#include "sage/rped/rpfile.h"

#define DEBUG_RPEDFILE(x)

namespace SAGE {
namespace RPED {

//==================================================================
//
//  CONSTRUCTOR
//
//==================================================================
RefPedigreeFile::RefPedigreeFile(cerrorstream &err) : errors(err)
{
  set_treat_as_sibs          (false);
  set_verbose_output         (10);
  set_reject_partial_lineage ();
  set_require_record         (false);
  set_skip_traits            (false);
  set_skip_markers           (false);
  set_dynamic_markers        (false);
  set_sex_linked_exist       (false);
  set_no_sex_field           (false);
  set_no_sex_ok_option       (false);
  set_sex_code_trait         (false);
  set_pedigree_id_trait      (false);
  set_sex_field_name         ("");
  set_pedigree_id_name       ("");
  reset_counts               ();
  validate                   ();

  my_inds.clear();
}

//==================================================================   
//
//  DESTRUCTOR
//
//==================================================================
RefPedigreeFile::~RefPedigreeFile() 
{
  errors.flush();
}

//==================================================================   
//
//  reset_counts()
//
//==================================================================
void 
RefPedigreeFile::reset_counts()
{
  my_skip_count           = 0;
  my_study_id_count       = 0;
  my_pedigree_id_count    = 0;
  my_individual_id_count  = 0;
  my_parent_id_count      = 0;
  my_sex_count            = 0;
  my_marker_count         = 0;
  my_invalid_marker_count = 0;
  my_invalid_marker_cov_count = 0;
  my_trait_count          = 0;
  my_invalid_trait_count  = 0;
  my_string_count         = 0;
  my_invalid_string_count = 0;
}

//==================================================================   
//
//  input(...)
//
//==================================================================
bool 
RefPedigreeFile::input(RefMultiPedigree &p, const string &filename, ostream &messages)
{
  // Read in pedigree meta-information:
  if(!input_pedigree(p, filename, messages, false))
    return false;

  // Read in actual pedigree data:
  if(!input_data(p, filename, messages, true))
    return false;

  if(sex_linked_exist())
    update_sex_linked_marker_info(p);

  return do_no_sex_structural_test(p);
}

//==================================================================   
//
//  print_mped(...)
//
//==================================================================
void 
RefPedigreeFile::print_mped(
	const RefMultiPedigree & p, 
	const string           & filename, 
	      ostream          & messages,
	      bool               dump_trait, 
	      bool               dump_marker)
{
  const RefMPedInfo & mped_info = p.info();

  print_family_structure_header(messages, filename);

  // Loop through the first X individuals and print out their info:
  for(size_t i = 0; i < my_ind_list.size(); ++i)
  {
    std::string              ped_name = my_ind_list[i].first,
                             ind_name = my_ind_list[i].second;
    RPED::MemberConstPointer ind      = p.member_find(ped_name, ind_name);
    std::string              p1       = ind->parent1() ? ind->parent1()->name() : mped_info.individual_missing_code(),
                             p2       = ind->parent2() ? ind->parent2()->name() : mped_info.individual_missing_code(),
                             ind_sex  = "";

    if( ind->is_male() )
      ind_sex = mped_info.sex_code_male    ();
    else if( ind->is_female() )
      ind_sex = mped_info.sex_code_female  ();
    else
      ind_sex = mped_info.sex_code_unknown ();

    print_family_structure(messages, ped_name, ind_name, ind_sex, p1, p2);

  } // End loop across first X individuals

  print_family_structure_footer(messages);

  if( dump_trait && (mped_info.trait_count() || mped_info.string_count()) )
  {
    print_trait_header(messages, mped_info, filename);

    for( size_t i = 0; i < my_ind_list.size(); ++i )
    {
      std::string pn = my_ind_list[i].first,
                  id = my_ind_list[i].second;

      RefMultiPedigree::member_const_pointer mem = p.member_find(pn, id);

      const RefPedInfo &ped_info = mem->pedigree()->info();

      std::vector<std::pair<size_t,std::string> > trait_values;

      size_t sti = mped_info.trait_find("SEX_CODE");
      if( sti < mped_info.trait_count() )
      {
        std::string t_value = mped_info.trait_info(sti).string_missing_code();

        if( !ped_info.trait_missing(mem->index(), sti) )
          t_value = doub2str(ped_info.trait(mem->index(), sti));

        trait_values.push_back(make_pair(sti, t_value));
      }

      size_t t1 = mped_info.trait_find("FAMILIAL_INDICATOR");
      size_t t2 = mped_info.trait_find("FOUNDER_INDICATOR");
      size_t t3 = mped_info.trait_find("PEDIGREE_SIZE");

      if( t1 < mped_info.trait_count() )
      {
        std::string t_value = mped_info.trait_info(t1).string_missing_code();

        if( !ped_info.trait_missing(mem->index(), t1) )
          t_value = doub2str(ped_info.trait(mem->index(), t1));

        trait_values.push_back(make_pair(t1, t_value));
      }

      if( t2 < mped_info.trait_count() )
      {
        std::string t_value = mped_info.trait_info(t2).string_missing_code();

        if( !ped_info.trait_missing(mem->index(), t2) )
          t_value = doub2str(ped_info.trait(mem->index(), t2));

        trait_values.push_back(make_pair(t2, t_value));
      }

      if( t3 < mped_info.trait_count() )
      {
        std::string t_value = mped_info.trait_info(t3).string_missing_code();

        if( !ped_info.trait_missing(mem->index(), t3) )
          t_value = doub2str(ped_info.trait(mem->index(), t3));

        trait_values.push_back(make_pair(t3, t_value));
      }

      for(field_list_type::const_iterator field_info = my_fields.begin(); field_info != my_fields.end(); ++field_info )
      {
        if( field_info->type != trait )
          continue;

        if( field_info->index < mped_info.trait_count() )
        {
          double val     = ped_info.trait(mem->index(), field_info->index);

          trait_values.push_back(make_pair(field_info->index, doub2str(val)));
        }
      }

      std::vector<std::pair<size_t,std::string> > string_values;

      for( field_list_type::const_iterator field_info = my_fields.begin(); field_info != my_fields.end(); ++field_info )
      {
        if( field_info->type == string_field && field_info->index < mped_info.string_count() )
        {
          string_values.push_back(make_pair(field_info->index, ped_info.get_string(mem->index(), field_info->index)));
        }
      }

      print_trait(messages, mped_info, pn, id, trait_values, string_values);
    }

    print_trait_footer(messages);
  }

  if(dump_marker && mped_info.marker_count())
  {
    print_marker_header(messages, mped_info, filename);

    for( size_t i = 0; i < my_ind_list.size(); ++i )
    {
      std::string pn = my_ind_list[i].first;
      std::string id = my_ind_list[i].second;

      RefMultiPedigree::member_const_pointer mem = p.member_find(pn, id);

      const RefPedInfo &ped_info = mem->pedigree()->info();

      std::vector<std::pair<size_t, std::string> >  marker_values;
      set<size_t> marker_added;

      for(field_list_type::const_iterator field_info = my_fields.begin() ; field_info != my_fields.end(); ++field_info)
      {
        if( field_info->type != marker && field_info->type != allele )
          continue;

        if( marker_added.find(field_info->index) != marker_added.end() )
          continue;

        if( field_info->index < mped_info.marker_count())
        {
          marker_added.insert(field_info->index);
          size_t m_index = field_info->index;

          const RefMarkerInfo& minfo = mped_info.marker_info(m_index);

          uint pheno_id = minfo.get_missing_phenotype_id();

          if( !ped_info.phenotype_missing(mem->index(), m_index, minfo) )
            pheno_id = ped_info.phenotype(mem->index(), m_index);

          std::string pheno_name = minfo.get_phenotype(pheno_id).name();

          if( pheno_id == minfo.get_missing_phenotype_id() )
          {
            std::string a = minfo.missing_allele_name();
            if( a == "*missing" )
              a = "?";
            pheno_name = a + minfo.gmodel().unphased_separator() + a;
          }
          else if( minfo.codominant(pheno_id) )
          {
            std::string al1, al2;
            MLOCUS::Ordering order;

            std::string geno_name = minfo.unphased_penetrance_begin(pheno_id).unphased_geno().name();

            minfo.gmodel().parse_genotype_name(geno_name, al1, al2, order);

            if( !al1.size() || al1.substr(0,1) == "~" ) al1 = minfo.missing_allele_name();
            if( !al2.size() || al2.substr(0,1) == "~" ) al2 = minfo.missing_allele_name();

            if( al1 == "*missing" )
              al1 = "?";

            if( al2 == "*missing" )
              al2 = "?";

            pheno_name = al1 + minfo.gmodel().unphased_separator() + al2;
          }

          marker_values.push_back(make_pair(m_index, pheno_name));
        }
      }

      print_marker(messages, mped_info, pn, id, marker_values);
    }

    print_marker_footer(messages);
  }
}

//==================================================================   
//
//  print_family_structure_header(...)
//
//==================================================================
void 
RefPedigreeFile::print_family_structure_header(ostream &messages, const std::string &filename) const
{
  messages << std::endl 
           << "Family structure information on the first " 
           << verbose_output() 
           << " individuals read from file: " 
           << filename
           << std::endl 
           << std::endl
           << "     PED. ID       IND. ID       SEX       PARENT1       PARENT2     " << std::endl
           << "     ------------  ------------  --------  ------------  ------------" << std::endl;
}

//==================================================================   
//
//  print_family_structure(...)
//
//==================================================================
void 
RefPedigreeFile::print_family_structure(
	      ostream & messages,
        const std::string  & pn, 
	const std::string  & id,
        const std::string  & sex, 
        const std::string  & p1, 
	const std::string  & p2) const
{
  messages << left
           << "     " << setw(12) << pn
           << "  "    << setw(12) << id 
           << "  "    << setw(8)  << sex  
           << "  "    << setw(12) << p1  
           << "  "    << setw(12) << p2 
           << std::endl;
}

//==================================================================   
//
//  print_family_structure_footer(...)
//
//==================================================================
void 
RefPedigreeFile::print_family_structure_footer(ostream &messages) const
{
  messages << std::endl;
}

//==================================================================   
//
//  print_trait_header(...)
//
//==================================================================
void 
RefPedigreeFile::print_trait_header(ostream &messages, const RefMPedInfo &mped_info, const std::string &filename) const
{
  messages << std::endl
           << "Phenotypes for the first "
           << verbose_output()
           << " individuals read from file: " 
           << filename
           << std::endl
           << std::endl
           << "     PED. ID       IND. ID     " << left;

  size_t num_trait_printed  = 0;
  size_t num_string_printed = 0;

  size_t sti = mped_info.trait_find("SEX_CODE");
  if( sti < mped_info.trait_count() )
  {
    messages << "  " << setw(20) << mped_info.trait_info(sti).name();
    ++num_trait_printed;
  }

  size_t t1 = mped_info.trait_find("FAMILIAL_INDICATOR");
  size_t t2 = mped_info.trait_find("FOUNDER_INDICATOR");
  size_t t3 = mped_info.trait_find("PEDIGREE_SIZE");

  if( t1 < mped_info.trait_count() )
  {
    messages << "  " << setw(20) << mped_info.trait_info(t1).name();
    ++num_trait_printed;
  }

  if( t2 < mped_info.trait_count() )
  {
    messages << "  " << setw(20) << mped_info.trait_info(t2).name();
    ++num_trait_printed;
  }

  if( t3 < mped_info.trait_count() )
  {
    messages << "  " << setw(20) << mped_info.trait_info(t3).name();
    ++num_trait_printed;
  }

  for(field_list_type::const_iterator field_info = my_fields.begin() ; field_info != my_fields.end(); ++field_info)
  {
    if( field_info->type != trait )
      continue;

    if( field_info->index < mped_info.trait_count() )
    {
      messages << "  " << setw(20) << mped_info.trait_info(field_info->index).name();
      ++num_trait_printed;
    }
  }

  for(field_list_type::const_iterator field_info = my_fields.begin(); field_info != my_fields.end() ; ++field_info)
  {
    if(field_info->type == string_field && field_info->index < mped_info.string_count())
    {
      messages << "  " << setw(20) << mped_info.string_info(field_info->index).name();
      ++num_string_printed;
    }
  }

  messages << std::endl << "     ------------  ------------";

  for( size_t i = 0; i < num_trait_printed + num_string_printed; ++i )
    messages << "  --------------------";

  messages << std::endl;
}

//==================================================================   
//
//  print_trait(...)
//
//==================================================================
void 
RefPedigreeFile::print_trait(
        ostream                                     & messages, 
  const RefMPedInfo                                 & mped_info,
  const std::string                                 & pn, 
  const std::string                                 & id,
  const std::vector<std::pair<size_t,std::string> > & trait_values,
  const std::vector<std::pair<size_t,std::string> > & string_values) const
{
  messages << "     " << setw(12) << pn << "  " << setw(12) << id;

  for(size_t i=0; i < trait_values.size(); ++i)
  {
    if( trait_values[i].first < mped_info.trait_count())
    {
      messages << "  " << setw(20) << trait_values[i].second;
    }
  }

  for( size_t i=0; i < string_values.size(); ++i)
  {
    if( string_values[i].first < mped_info.string_count())
    {
      messages << "  " << setw(20) << string_values[i].second;
    }
  }

  messages << std::endl;
}

//==================================================================   
//
//  print_trait_footer(...)
//
//==================================================================
void 
RefPedigreeFile::print_trait_footer(ostream &messages) const
{
  messages << std::endl;
}

//==================================================================   
//
//  print_marker_header(...)
//
//==================================================================
void 
RefPedigreeFile::print_marker_header(ostream &messages, const RefMPedInfo &mped_info, const std::string &filename) const
{
  messages << std::endl
           << "Markers for the first " 
           << verbose_output()
           << " individuals read from file: " 
           << filename
           << std::endl 
           << std::endl
           << "     PED. ID       IND. ID     " << left;

  size_t num_marker_printed  = 0;

  set<size_t> marker_added;

  for(field_list_type::const_iterator field_info = my_fields.begin() ; field_info != my_fields.end(); ++field_info)
  {
    if((field_info->type == marker || field_info->type == allele)   &&
       (marker_added.find(field_info->index) == marker_added.end()) &&
       (field_info->index < mped_info.marker_count()))
    {
      marker_added.insert(field_info->index);

      messages << "  " << setw(12) << mped_info.marker_info(field_info->index).name();

      ++num_marker_printed;
    }
  }

  messages << std::endl << "     ------------  ------------";

  for( size_t i = 0; i < num_marker_printed; ++i )
    messages << "  ------------";

  messages << std::endl;
}

//==================================================================   
//
//  print_marker(...)
//
//==================================================================
void 
RefPedigreeFile::print_marker(
        std::ostream                                 & messages, 
  const RefMPedInfo                                  & mped_info,
  const std::string                                  & pn, 
  const std::string                                  & id,
  const std::vector<std::pair<size_t, std::string> > & marker_values) const
{
  messages << "     " << setw(12) << pn << "  " << setw(12) << id;

  for(size_t m = 0; m < marker_values.size(); ++m)
  {
    if(marker_values[m].first < mped_info.marker_count())
    {
      messages << "  " << setw(12) << marker_values[m].second;
    }
  }

  messages << std::endl;
}

//==================================================================   
//
//  print_marker_footer(...)
//
//==================================================================
void 
RefPedigreeFile::print_marker_footer(ostream &messages) const
{
  messages << std::endl;
}

//==================================================================   
//
//  parse_dynamic_markers_missing(...)
//
//==================================================================
void 
RefPedigreeFile::parse_dynamic_markers_missing(RefMarkerInfo& marker, const std::string& missing)
{
  char        sep0 = marker.gmodel().separators()[0];
  std::string mv   = missing;

  if( strchr(mv.c_str(), sep0) )
  {
    std::string mv2 = strip_ws(mv.substr(mv.find_first_of(sep0), mv.size()), "/ \t");

    mv = strip_ws(mv.substr(0, mv.find_first_of(sep0)), " \t");

    if( mv != mv2 )
    {
      errors << priority(warning) 
             << "The Missing Value code must use the "
             << "same value for both alleles for '"
             << marker.name() 
             << "'.  Will use " 
             << mv 
             << "/" 
             << mv
             << " for the missing value." 
             << std::endl;
    }
  }

  if(mv != marker.missing_allele_name())
  {
    marker.gmodel().set_missing_allele_name(mv);
    marker.set_missing_phenotype_name(mv);

    std::string sep = marker.gmodel().separators();

    marker.alias_phenotype(marker.get_missing_phenotype_id(), mv + sep[0] + mv);
    marker.alias_phenotype(marker.get_missing_phenotype_id(), mv + sep[1] + mv);
    marker.alias_phenotype(marker.get_missing_phenotype_id(), mv + sep[2] + mv);
  }
}

//==================================================================   
//
//  validate_fields(...)
//
//==================================================================
bool 
RefPedigreeFile::validate_fields(bool data_only, bool quiet)
{
  // NOTE: We treat invalid markers and traits as non-fatal states

  reset_counts();

  typedef std::map<std::string, size_t> name_count_map;

  name_count_map trait_map,
                 string_map,
                 marker_map,
                 marker_cov_map;

  for(field_list_type::iterator i = my_fields.begin(); i != my_fields.end(); ++i)
  {
    switch( i->type )
    {
      case skip          : ++my_skip_count;                    break;
      case study_id      : ++my_study_id_count;                break;
      case pedigree_id   : ++my_pedigree_id_count;             break;
      case individual_id : ++my_individual_id_count;           break;
      case parent_id     : ++my_parent_id_count;               break;
      case sex_code      : ++my_sex_count;                     break;

      case trait         : trait_map  [toUpper(i->name)]++;    break;
      case string_field  : string_map [toUpper(i->name)]++;    break;
      case marker        : marker_map [toUpper(i->name)] += 2; break;
      case allele        : marker_map [toUpper(i->name)]++;    break;
    }
  }

  // Set validity true to get started:
  bool valid_fields = true;

  // Check structural stuff:
  if(!data_only)
  {
    // Wrong count of any particular field:
    
    // Check study id count:
    if (my_study_id_count > 1 && !quiet)
    {
      errors << priority(warning) << "More than one study ID specified." << std::endl;
    }

    // Individual id:
    if(my_individual_id_count == 0)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "Individual ID field not found." << std::endl;
    }
    else if(my_individual_id_count > 1)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "Too many Individual ID fields found." << std::endl;
    }

    // Ped id:
    if(my_pedigree_id_count > 1)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "More than one pedigree id field indicated." << std::endl;
    }
  
    // Parent id:
    if(my_parent_id_count == 1)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "Only one parent id field indicated; either none or two are required." << std::endl;
    }
    else if(my_parent_id_count > 2)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "More than two parent id fields indicated." << std::endl;
    }
  
    // Parent id:
    if(my_sex_count > 1)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "More than one sex field indicated." << std::endl;
    }
  
    // Moving on to more complicated structural checking:  
    
    // Parents present, sex missing, but no_sex_field not set to true:
    if(my_parent_id_count == 2 && my_sex_count == 0 && my_no_sex_field == false)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "Parent ID's present, but sex field is missing and 'NO_SEX_FIELD' has not been specified." << std::endl;
    }
    
    // Parents present but treat_as_sibs specified:
    if(my_parent_id_count == 2 && get_treat_as_sibs() == true)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "User has specified 'TREAT_AS_SIBS', but still listed parent ids." << std::endl;
    }
    
    // Sex present but NO_SEX_FIELD also present:
    if(my_sex_count == 1 && my_no_sex_field == true)
    {
      valid_fields = false;
      if(!quiet)
        errors << priority(error) << "Both 'SEX_FIELD' and 'NO_SEX_FIELD' were specified in the parameter file." << std::endl;
    }
    
    // Sex missing, NO_SEX_FIELD specified, deliver an advance warning (but analysis still valid):
    if(my_sex_count == 0 && my_no_sex_field == true)
    {
      if(!quiet) 
        errors << priority(warning) << "Please note that no sex field is present in the data; "
                                    << "analyses that depend on sex-specific data may produce unpredictable results." 
                                    << std::endl;
    }
        
  } // End structural checking
  
  // Check traits:

  for(field_list_type::iterator i = my_fields.begin(); i != my_fields.end(); ++i)
  {
    if( i->type == trait )
    {
      size_t &count = trait_map[ toUpper(i->name) ];

      if( count == 1 )
      {
        ++my_trait_count;
      }
      else if( count != (size_t)-1 )  // Invalidate new bad trait
      {
        if(!quiet)
          errors << priority(warning) << "Trait assigned to more than one field.  Skipping trait '" << i->name << "'." << std::endl;

        // Count trait as bad and mark count to avoid further warnings
        ++my_invalid_trait_count;
        count = (size_t)-1;
        i->index = (size_t)-1;
      }
      else // invalidate known bad trait
      {
        i->index = (size_t)-1;
      }
    }
    else if( i->type == string_field )
    {
      size_t &count = string_map[ toUpper(i->name) ];

      if( count == 1 )
      {
        ++my_string_count;
      }
      else if( count != (size_t)-1 )  // Invalidate new bad string
      {
        if(!quiet)
          errors << priority(warning) << "String field assigned to more than one field.  Skipping field '" << i->name << "'." << std::endl;

        // Count field as bad and mark count to avoid further warnings
        ++my_invalid_string_count;
        count = (size_t)-1;
        i->index = (size_t)-1;
      }
      else // invalidate known bad field
      {
        i->index = (size_t)-1;
      }
    }
    else if( i->type == marker || i->type == allele )
    {
      size_t &count = marker_map[ toUpper(i->name) ];

      if( count == 2 )
      {
        ++my_marker_count;
      }
      else if( count != (size_t)-1 )  // Invalidate new bad marker
      {
        if(!quiet && count < 2)
        {
          errors << priority(warning) << "Marker assigned to too few allele fields.  Skipping marker '" << i->name << "'." << std::endl;
        }
        else if(!quiet && count > 2)
        {
          errors << priority(warning) << "Marker assigned to too many fields.  Skipping marker '" << i->name << "'." << std::endl;
        }

        // Count marker as bad and mark count to avoid further warnings
        ++my_invalid_marker_count;

        count    = (size_t)-1;
        i->index = (size_t)-1;
      }
      else // invalidate known bad marker
      {
        i->index = (size_t)-1;
      }
    }
  }

  // NOTE: We treat invalid markers and traits as non-fatal states

  return valid_fields;
}

//==================================================================   
//
//  build_pedigree(...)
//
//==================================================================
bool 
RefPedigreeFile::build_pedigree(RefMultiPedigree &p)
{
  const RefMPedInfo &mped_info = p.info();

  p.build();
  
  bool found_fatal_errors = false;

  // Construct the pedigrees and fill in index information
  // (much of this should/will be moved into the RefPedigree directly
  //  as a build action callback)
  RefMultiPedigree::pedigree_iterator i = p.pedigree_begin();
  for( int j=0; i != p.pedigree_end(); ++i, ++j)
  {
    if( i->error_count() )
    {
      if(!report_pedigree_build_errors(*i))
        found_fatal_errors = true;
    }
    
    if( i->warning_count() )
    {
      if(!report_pedigree_build_warnings(*i))
        found_fatal_errors = true;
    }

    // Create indices
    i->info().build(*i);
    i->info().resize_traits(  mped_info.trait_count() );
    i->info().resize_strings( mped_info.string_count() );
    i->info().resize_markers( mped_info.marker_count(), mped_info );
  }

  return !found_fatal_errors;
}

/// Reports structural errors detected when building the pedigree and a count
/// of the number of errors.  Returns a true/false based upon if any of the errors
/// are fatal (ie, unignorable).  Fatal errors currently include inconsistencies
/// in parental sex.
bool
RefPedigreeFile::report_pedigree_build_errors(const Pedigree &p) const
{
  bool found_fatal_error = false;

  // Create a buffer stream to hold the errors.  We want to report the error count
  // first, but some errors are suppressed in certain circumstances, so we need 
  // to count them first.  We can use a buffer to queue the messages.
  //SAGE::bufferederrorstream<char> buf(errors);
  stringstream buf;

  // Create a counter to count the reported errors
  size_t error_count = 0;
  
  RefMultiPedigree::pedigree_type::error_iterator err;
  for( err = p.error_begin(); err != p.error_end(); ++err )
  {
    switch( err->state )
    {
      case MPED::error_info::bad_sibship:
        buf //<< priority(warning)
            << "     Error in sibship due to individuals: '" 
            << err->name1 << "', '" << err->name2 << "'." << std::endl;
        ++error_count;
        break;

      case MPED::error_info::bad_marriage:
        buf //<< priority(warning)
            << "     Inconsistent mating due to individuals: '" 
            << err->name1 << "', '" << err->name2 << "'." << std::endl;
        ++error_count;
        break;

      case MPED::error_info::bad_lineage:
        buf //<< priority(warning)
            << "     Inconsistent lineage due to individuals: '" 
            << err->name2 << "', '" << err->name3 << "' -> '"
            << err->name1 << "'." << std::endl;
        ++error_count;
        break;

      case MPED::error_info::bad_gender:
        buf //<< priority(critical)
            << "     Individual '" << err->name1
            << "' is present as both male and female." << std::endl;
        found_fatal_error = true;
        ++error_count;
        break;

      case MPED::error_info::same_sex_marriage:
        buf //<< priority(critical)
            << "Both parents '" << err->name1 << "' and '" << err->name2
            << "' are labeled as or inferred to be " << err->name3 << "."
            << std::endl;
        found_fatal_error = true;
        ++error_count;
        break;
        
      case MPED::error_info::bad_marriage_loop:
        buf //<< priority(critical)
            << "Individual '" << err->name1 << "' belongs to a marriage loop "
            << "of odd length.  Such marriage loops cannot be assigned sex "
            << "consistently." << std::endl;
        found_fatal_error = true;
        ++error_count;
        break;
        
      case MPED::error_info::no_sex_parents:
        buf //<< priority(warning)
            << "Both parents '" << err->name1 << "' and '" << err->name2
            << "' have unknown sex and cannot be inferred." << std::endl;
        ++error_count;
        break;
        
      default: break;
    }
  }
  
  // If we have reported errors, report the number of errors first, then dump our
  // buffer
  if(error_count)
  {
    errors << priority(error) << "Error(s) building pedigree '" << p.name() << "': "
           << error_count << " errors." << std::endl;
    errors << buf.str() << std::endl;
    //buf.flush_buffer();
  }
  
  return !found_fatal_error;
}

/// Reports structural warnings detected when building the pedigree.
/// Returns a true/false based upon if any of the errors
/// are fatal (ie, unignorable).  In the current version, no warning is fatal,
/// but the boolean is included for consistency.
bool
RefPedigreeFile::report_pedigree_build_warnings(const Pedigree &p) const
{
  RefMultiPedigree::pedigree_type::error_iterator err;
  for( err = p.warning_begin(); err != p.warning_end(); ++err )
  {
    switch( err->state )
    {
      case MPED::error_info::gender_inferred:
        errors << priority(warning)
               << "The sex of individual '" << err->name1
               << "' in pedigree '" << p.name()
               << "' is believed to be " << err->name2
               << ".  Sex code will be reset." << std::endl;
        break;

      default: break;
    }
  }
  
  return true;
}    

//==================================================================   
//
//  build_data(...)
//
//==================================================================
bool 
RefPedigreeFile::build_data(RefMultiPedigree &p)
{
  RefMPedInfo &mped_info = p.info();

  // fill in index information
  RefMultiPedigree::pedigree_iterator i = p.pedigree_begin();
  for( int j=0; i != p.pedigree_end(); ++i, ++j)
  {
    // Create indices
    i->info().build(*i);
    i->info().resize_traits(  mped_info.trait_count() );
    i->info().resize_strings( mped_info.string_count() );
    i->info().resize_markers( mped_info.marker_count(), mped_info );
  }

  // Build implicit traits
  if( sex_code_trait() )
  {
    size_t t = mped_info.trait_find("SEX_CODE");

    if( t < mped_info.trait_count() )
    {
      RefMultiPedigree::member_const_iterator j;

      for( i=p.pedigree_begin(); i != p.pedigree_end(); ++i )
      {
        RefMultiPedigree::pedinfo_type &info = i->info();

        for( j=i->member_begin(); j != i->member_end(); ++j )
        {
          int ind_num = j->index();

          if(j->is_male())
            info.set_trait(ind_num, t, mped_info.sex_code_male(), mped_info.trait_info(t));
          else if(j->is_female())
            info.set_trait(ind_num, t, mped_info.sex_code_female(), mped_info.trait_info(t));
        }
      }
    }
    ++my_trait_count;
  }

  if( pedigree_id_trait() )
  {
    size_t t1 = mped_info.trait_find("FAMILIAL_INDICATOR");
    size_t t2 = mped_info.trait_find("FOUNDER_INDICATOR");
    size_t t3 = mped_info.trait_find("PEDIGREE_SIZE");

    if(    t1 < mped_info.trait_count()
        && t2 < mped_info.trait_count()
        && t3 < mped_info.trait_count() )
    {
      RefMultiPedigree::member_const_iterator j;

      for( i=p.pedigree_begin(); i != p.pedigree_end(); ++i )
      {
        RefMultiPedigree::pedinfo_type &info = i->info();

        for( j=i->member_begin(); j != i->member_end(); ++j )
        {
          int ind_num = j->index();

          if( j->pedigree()->member_count() > 1 )
            info.set_trait(ind_num, t1, "1", mped_info.trait_info(t1));
          else
            info.set_trait(ind_num, t1, "0", mped_info.trait_info(t1));

          if( j->is_founder() )
            info.set_trait(ind_num, t2, "1", mped_info.trait_info(t2));
          else
            info.set_trait(ind_num, t2, "0", mped_info.trait_info(t2));

          info.set_trait(ind_num, t3, j->pedigree()->member_count());

        }
      }
    }

    my_trait_count += 3;
  }

  return true;
}

//==================================================================   
//
//  add_member(...)
//
//==================================================================
void 
RefPedigreeFile::add_member(
        RefMultiPedigree & multipedigree, 
  const std::string      & ped_name,
  const std::string      & ind_name, 
  const std::string      & sex,
  const std::string      & parent1, 
  const std::string      & parent2,
        size_t             line,
        size_t             count)
{
  if(count == (size_t)-1)
  {
     count = line;
  }
  
  if(!ped_name.size())
  {
    errors << priority(warning) << "[" << line << "] record is missing pedigree ID.  Skipping..." << std::endl;

    return;
  }

  if(ind_name == "" || ind_name == multipedigree.info().individual_missing_code())
  {
    errors << priority(warning) << "[" << line
           << "] Found an individual ID that is missing or"
           << " the same as the missing individual/parent code.  Skipping..."
           << std::endl;

    return;
  }
  
  if(reject_partial_lineage() && ((parent1 == multipedigree.info().individual_missing_code()) ^ (parent2 == multipedigree.info().individual_missing_code())))
  {
    errors << priority(warning) << "[" << line << "] Individual '" << ind_name << "' in pedigree '" << ped_name << "' has only one parent specified and will be skipped." << std::endl;

    return;
  }

  // Add entries for the given parents, assuming require_record is set to false (the default).
  if(!require_record())
  { 
    if(parent1 != multipedigree.info().individual_missing_code())
      multipedigree.add_member(ped_name, parent1, MPED::SEX_ARB);
      
    if(parent2 != multipedigree.info().individual_missing_code())
      multipedigree.add_member(ped_name, parent2, MPED::SEX_ARB);
  }

  if(sex == multipedigree.info().sex_code_male())
  {
    multipedigree.add_member(ped_name, ind_name, MPED::SEX_MALE);
  }
  else if(sex == multipedigree.info().sex_code_female())
  {
    multipedigree.add_member(ped_name, ind_name, MPED::SEX_FEMALE);
  }
  else
  {
    if(sex != multipedigree.info().sex_code_unknown())
    {
      errors << priority(warning) << "[" << line
             << "] Unable to read sex of individual "  << ind_name
             << " in pedigree " << ped_name
             << ".  Found '" << sex
             << "'.  Setting sex to unknown." << std::endl;
    }
    
    multipedigree.add_member(ped_name, ind_name, MPED::SEX_MISSING);
  }

  if(parent1 != multipedigree.info().individual_missing_code())
  {
    if(parent2 != multipedigree.info().individual_missing_code())
    {
      multipedigree.add_lineage(ped_name, ind_name, parent1, parent2);
    }
    else
    {
      multipedigree.add_lineage(ped_name, ind_name, parent1);
    }
  }
  else if(parent2 != multipedigree.info().individual_missing_code())
  {
    multipedigree.add_lineage(ped_name, ind_name, parent2);
  }

  // Add person to list of individuals:
  my_inds.insert(make_pair(ped_name, ind_name));

  // Check to see if there's a duplicate record:
  if(my_inds.count(make_pair(ped_name, ind_name)) > 1)
  {
    errors << priority(warning) 
           << "[" 
           << line
           << "] Duplicate record of individual "  
           << ind_name
           << " in pedigree " 
           << ped_name
           << "." 
           << std::endl;
  }
  else if(count <= verbose_output()) // Ok, there's no duplicate record, so stick 'em in the initial list:
  {
    my_ind_list.push_back(make_pair(ped_name, ind_name));
  }
}

//==================================================================   
//
//  update_marker_delimeter_info(...)
//
//==================================================================
void 
RefPedigreeFile::update_marker_delimiter_info(RefMarkerInfo& marker, char sep0)
{
  if( sep0 != marker.gmodel().separators()[0] )
  {
    marker.gmodel().set_unphased_separator(sep0);

    RefMarkerInfo::phenotype_iterator phi = marker.phenotype_begin();
    for( ; phi != marker.phenotype_end(); ++phi )
    {
      uint     id = phi->id();
      std::string   name = phi->name();
      std::string   allele1, allele2;
      MLOCUS::Ordering order;

      marker.gmodel().parse_genotype_name(name, allele1, allele2, order);

      std::string sep = marker.gmodel().separators();

      marker.alias_phenotype(id, allele1 + sep[0] + allele2);
      marker.alias_phenotype(id, allele1 + sep[1] + allele2);
      marker.alias_phenotype(id, allele1 + sep[2] + allele2);
    }
  }
}

//==================================================================   
//
//  update_marker_missing_info(...)
//
//==================================================================
void 
RefPedigreeFile::update_marker_missing_info(RefMarkerInfo& marker, const std::string& missing)
{
  parse_dynamic_markers_missing(marker, missing);
}

//==================================================================   
//
//  update_sex_linked_marker_info(...)
//
//==================================================================
void
RefPedigreeFile::update_sex_linked_marker_info(RefMultiPedigree &mp)
{
}

//
//--------------------------------------------------------------------------------------
//

bool
RefPedigreeFile::do_no_sex_structural_test(const RefMultiPedigree& mp)
{
  if(get_no_sex_ok_option()) return true;
  
  for(PedigreeConstIterator ped = mp.pedigree_begin(); ped != mp.pedigree_end(); ++ped)
  {
    // If the pedigree has no families, there is no structure, so we  skip it
    if(!ped->family_count()) continue;

    // Iterate through the subpedigrees
    for(SubpedigreeConstIterator sped = ped->subpedigree_begin();
        sped != ped->subpedigree_end(); ++sped)
    {
      // Iterate through the members, looking for members without sex information.
      for(MemberConstIterator mem = sped->member_begin(); mem != sped->member_end();
          ++mem)
      {
        // If we encounter an individual without sex, print the warning and exit.
        if(mem->is_sex_unknown())
        {
          errors << priority(critical)
                 << "Many S.A.G.E. algorithms rely on sex for pedigree structure "
                 << "information. Even when sex does not affect the outcome of an "
                 << "analysis, the lack of sex information may cause the program to "
                 << "behave improperly or even crash. Individuals without sex "
                 << "information, either direct or inferred, but with pedigree "
                 << "structure information have been detected in your data set. If "
                 << "your analyses involve pedigree structure, we recommend that "
                 << "this be corrected before continuing your analyses. If you would "
                 << "like to continue the analysis without correcting this issue, "
                 << "you must place a \"no_sex_ok\" parameter into your pedigree "
                 << "block." << endl;
                 
          return false;
        }
      }
    }
  }
  
  return true;
}

//
//--------------------------------------------------------------------------------------
//

//==================================================================   
//
//  CONSTRUCTOR
//
//==================================================================
RefLSFPedigreeFile::RefLSFPedigreeFile(cerrorstream &errors) : RefPedigreeFile(errors)
{}

//==================================================================   
//
//  process_parameters(...)
//
//==================================================================
bool
RefLSFPedigreeFile::process_parameters(RefMPedInfo &mped_info, const LSFBase* params)
{
  if(!params)
  {
    errors << priority(critical) << "RefPedigreeFile Error: No parameters specified." << std::endl;

    invalidate();
    
    return false;
  }
  else // params is valid
  {
    for(LSFList::const_iterator i=params->List()->begin(); i != params->List()->end(); ++i)
    {
      if(!process_parameter(mped_info, *i)) return false;
    }

    // Create sex_code covariate always as long as sex_field exist.
    //
    if( sex_field_name() != "" )
    {
      set_sex_code_trait(true);

      size_t t = mped_info.trait_find("SEX_CODE");

      // If not previously available, add it.
      if( t >= mped_info.trait_count() )
      {
        t = mped_info.add_binary_trait("SEX_CODE", RefTraitInfo::trait_covariate);

        mped_info.trait_info(t).set_string_missing_code     ("");
        mped_info.trait_info(t).set_numeric_missing_code    (-1);
        mped_info.trait_info(t).set_string_affected_code    (mped_info.sex_code_female());
        mped_info.trait_info(t).set_string_unaffected_code  (mped_info.sex_code_male());
        mped_info.trait_info(t).set_numeric_affected_code   (1);
        mped_info.trait_info(t).set_numeric_unaffected_code (0);
        mped_info.trait_info(t).set_alias_name(sex_field_name());
      }
    }

    // Create indicator covariates always as long as pedigree_id exist.
    //
    if( pedigree_id_name() != "" )
    {
      set_pedigree_id_trait(true);

      size_t t = mped_info.trait_find("FAMILIAL_INDICATOR");

      // If not previously available, add it.
      if( t >= mped_info.trait_count() )
      {
        t = mped_info.add_binary_trait("FAMILIAL_INDICATOR", RefTraitInfo::trait_covariate);

        mped_info.trait_info(t).set_string_missing_code     ("");
        mped_info.trait_info(t).set_numeric_missing_code    (-1);
        mped_info.trait_info(t).set_string_affected_code    ("1");
        mped_info.trait_info(t).set_string_unaffected_code  ("0");
        mped_info.trait_info(t).set_numeric_affected_code   (1);
        mped_info.trait_info(t).set_numeric_unaffected_code (0);
        mped_info.trait_info(t).set_alias_name("is_familiar");
      }

      t = mped_info.trait_find("FOUNDER_INDICATOR");

      // If not previously available, add it.
      if( t >= mped_info.trait_count() )
      {
        t = mped_info.add_binary_trait("FOUNDER_INDICATOR", RefTraitInfo::trait_covariate);

        mped_info.trait_info(t).set_string_missing_code     ("");
        mped_info.trait_info(t).set_numeric_missing_code    (-1);
        mped_info.trait_info(t).set_string_affected_code    ("1");
        mped_info.trait_info(t).set_string_unaffected_code  ("0");
        mped_info.trait_info(t).set_numeric_affected_code   (1);
        mped_info.trait_info(t).set_numeric_unaffected_code (0);
        mped_info.trait_info(t).set_alias_name("is_founder");
      }

      t = mped_info.trait_find("PEDIGREE_SIZE");

      // If not previously available, add it.
      if( t >= mped_info.trait_count() )
      {
        t = mped_info.add_continuous_trait("PEDIGREE_SIZE", RefTraitInfo::trait_covariate);

        mped_info.trait_info(t).set_string_missing_code     ("");
        mped_info.trait_info(t).set_numeric_missing_code    (std::numeric_limits<double>::quiet_NaN());
        mped_info.trait_info(t).set_alias_name("ped_size");
      }
    }
  }

  return true;
}

//==================================================================   
//
//  process_parameter(...)
//
//==================================================================
bool
RefLSFPedigreeFile::process_parameter(RefMPedInfo &mped_info, const LSFBase* param)
{
  if(!param)
    return true;

  // Set up local variables:
  
  std::string              name = toUpper(param->name());

       if(name == "FORMAT")                   return process_format                   (mped_info, param);
  else if(name == "VERBOSE")                  return process_verbose                  (mped_info, param);
  else if(name == "REQUIRE_RECORD")           return process_require_record           (mped_info, param);
  else if(name == "SKIP_TRAITS")              return process_skip_traits              (mped_info, param);
  else if(name == "SKIP_MARKERS")             return process_skip_markers             (mped_info, param);
  else if(name == "DYNAMIC_MARKERS")          return process_dynamic_markers          (mped_info, param);
  else if(name == "REJECT_PARTIAL_LINEAGE")   return process_reject_partial_lineage   (mped_info, param);
  else if(name == "INDIVIDUAL_MISSING_VALUE") return process_individual_missing_value (mped_info, param);
  else if(name == "SEX_CODE")                 return process_sex_code                 (mped_info, param);
  else if(name == "NO_SEX_OK")                return process_no_sex_ok                (mped_info, param);

  return true;
}

/// Processes the "FORMAT" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_format(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if(val.has_value())
  {
    set_format(val.String());
  }

  return true;
}

/// Processes the "VERBOSE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_verbose(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if( val.has_value() && val.Int() > 0 )
  {
    set_verbose_output(val.Int());
  }

  return true;
}

/// Processes the "REQUIRE_RECORD" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_require_record(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if(val.has_value())
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_require_record(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_require_record(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'require_record'" << std::endl;
    }
  }

  return true;
}

/// Processes the "SKIP_TRAITS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_skip_traits(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if(val.has_value() && !force_skip_traits())
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_skip_traits(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_skip_traits(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'skip_traits'" << std::endl;
    }
  }

  return true;
}

/// Processes the "SKIP_MARKERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_skip_markers(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if(val.has_value() && !force_skip_markers())
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_skip_markers(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_skip_markers(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'skip_markers'" << std::endl;
    }
  }

  return true;
}

/// Processes the "DYNAMIC_MARKERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_dynamic_markers(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if(val.has_value() && !force_dynamic_markers())
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_dynamic_markers(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_dynamic_markers(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'dynamic_markers'" << std::endl;
    }
  }

  return true;
}

/// Processes the "REJECT_PARTIAL_LINEAGE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_reject_partial_lineage(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if( val.has_value() )
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_reject_partial_lineage(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_reject_partial_lineage(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'reject_partial_lineage'" << std::endl;
    }
  }

  return true;
}

/// Processes the "INDIVIDUAL_MISSING_VALUE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_individual_missing_value(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if( val.has_value() )
  {
    mped_info.set_individual_missing_code(val.String());

    DEBUG_RPEDFILE(errors << priority(debug) << "Found individual missing value = " << mped_info.individual_missing_code() << std::endl;)
  }

  return true;
}

/// Processes the "SEX_CODE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_sex_code(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrList::const_iterator a;
  
  if(param->attrs())
  {
    if( (a=param->attrs()->find("missing")) != param->attrs()->end() )
    {
      mped_info.set_sex_code_unknown(strip_ws(a->second.String()));

      DEBUG_RPEDFILE(errors << priority(debug) << "Found sex missing value code = " << mped_info.sex_code_unknown() << std::endl;)
    }
    if( (a=param->attrs()->find("male"))    != param->attrs()->end() )
    {
      mped_info.set_sex_code_male(strip_ws(a->second.String()));

      DEBUG_RPEDFILE(errors << priority(debug) << "Found sex male code = " << mped_info.sex_code_male() << std::endl;)
    }
    if( (a=param->attrs()->find("female"))  != param->attrs()->end() )
    {
      mped_info.set_sex_code_female(strip_ws(a->second.String()));

      DEBUG_RPEDFILE(errors << priority(debug) << "Found sex female code = " << mped_info.sex_code_female() << std::endl;)
    }
  }

  return true;
}

/// Processes the "NO_SEX_OK" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool
RefLSFPedigreeFile::process_no_sex_ok(RefMPedInfo &mped_info, const LSFBase* param)
{
  AttrVal val = attr_value(param, 0);

  if( val.has_value() )
  {
    if( toUpper(val.String()) == "TRUE")
    {
      set_no_sex_ok_option(true);
    }
    else if( toUpper(val.String()) == "FALSE")
    {
      set_no_sex_ok_option(false);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'no_sex_ok'" << std::endl;
    }
  }

  return true;
}

//=======================================================================
//=======================================================================
//=======================================================================
//
//  RefDelimitedPedigreeFile
//
//=======================================================================
//=======================================================================
//=======================================================================

//=======================================================================
//
//  Constructor
//
//=======================================================================
RefDelimitedPedigreeFile::RefDelimitedPedigreeFile(cerrorstream &err) : RefPedigreeFile(err)
{
  set_format_in_file              (true);
  set_whitespace                  (" \r\n");
  set_delimiters                  ("\t,");
  set_skip_consecutive_delimiters (false);
  set_skip_leading_delimiters     (false);
  set_skip_trailing_delimiters    (false);
  reset_counts                    ();
}

//=======================================================================
//
//  Destructor
//
//=======================================================================
RefDelimitedPedigreeFile::~RefDelimitedPedigreeFile()
{}

//=======================================================================
//
//  set_delimiters(...)
//
//=======================================================================
void 
RefDelimitedPedigreeFile::set_delimiters(const std::string &d)
{
  my_delimiters = d;

  if(my_delimiters.find("\n") == std::string::npos)
    my_delimiters += "\n";

  if(my_delimiters.find("\r") == std::string::npos)
    my_delimiters += "\r";
}

//=======================================================================
//
//  set_whitespace(...)
//
//=======================================================================
void 
RefDelimitedPedigreeFile::set_whitespace(const std::string &w)
{
  my_whitespace = w;

  if(my_whitespace.find("\n") == std::string::npos)
    my_whitespace += "\n";

  if(my_whitespace.find("\r") == std::string::npos)
    my_whitespace += "\r";
}


//=======================================================================
//
//  build_fields(...)
//
//=======================================================================
bool 
RefDelimitedPedigreeFile::build_fields(string_tokenizer&  header,
                                       const RefMPedInfo& mped_info,
                                       bool               quiet)
{
  my_fields.clear();

  for(string_tokenizer::iterator i = header.begin(); i != header.end(); ++i)
  {
    my_fields.push_back(field());

    if(i->empty())
       continue;

    field_map_type::const_iterator f = my_field_map.find( toUpper(*i) );
    
    if(f != my_field_map.end())
    {
      // If the field is in our map
      build_field(mped_info, f, quiet);
    }
    else
    {
      DEBUG_RPEDFILE(errors << priority(information) << "Field '" << *i << "' in pedigree file is skipped." << std::endl;)
    }
  }

  return true;
}

void 
RefDelimitedPedigreeFile::build_field(const RefMPedInfo&             mped_info,
                                      field_map_type::const_iterator f,
                                      bool                           quiet)
{
  // Verify that the field should not be skipped
  if(f->second.type == skip)
  {
    DEBUG_RPEDFILE(errors << priority(information) << "Field '" << *i << "' in pedigree file is skipped." << std::endl;)
    return;
  }

  field & current_field = my_fields.back();

  // This makes sense -- current field is a reference
  current_field = f->second;

  if( current_field.type == trait )
  {
    build_trait_field(mped_info, current_field);
  }
  else if( current_field.type == string_field )
  {
    build_string_field(mped_info, current_field);
  }
  else if( current_field.type == marker || current_field.type == allele )
  {
    build_marker_field(mped_info, current_field);
  }

  DEBUG_RPEDFILE(errors << priority(information) << "Using field '" << *i << "'." << std::endl;)
}

void
RefDelimitedPedigreeFile::build_trait_field(const RefMPedInfo& mped_info,
                                            field&             current_field)
{
  size_t t = mped_info.trait_find(current_field.name);

  if(t < mped_info.trait_count())
  {
    current_field.index = t;
  }
  else
  {
    current_field.index = (size_t)-1;

    // Can never happen (yet) since we always add the trait to mped_info
    // in the LSF parser.  But we have to look out for the future.
    errors << priority(warning) << "Invalid trait '" << current_field.name << "'.  Skipping..." << std::endl;
  }
}

void
RefDelimitedPedigreeFile::build_marker_field(const RefMPedInfo& mped_info,
                                             field&             current_field)
{
  size_t m = mped_info.marker_find(current_field.name);

  if(m < mped_info.marker_count())
  {
    current_field.index = m;
  }
  else
  {
    current_field.index = (size_t)-1;

    if(mped_info.marker_count())
    {
      errors << priority(warning) << "Invalid marker '" << current_field.name << "'.  Skipping..." << std::endl;
    }
  }
}
void
RefDelimitedPedigreeFile::build_string_field(const RefMPedInfo& mped_info,
                                             field&             current_field)
{
  size_t s = mped_info.string_find(current_field.name);

  if(s < mped_info.string_count())
  {
    current_field.index = s;
  }
  else
  {
    current_field.index = (size_t)-1;

    // Can never happen (yet) since we always add the trait to mped_info
    // in the LSF parser.  But we have to look out for the future.
    errors << priority(warning) << "Invalid string field '" << current_field.name
           << "'.  Skipping..." << std::endl;
  }
}

//=================================================
//
//  check_header_vs_field_map(...)
//
//=================================================
bool 
RefDelimitedPedigreeFile::check_header_vs_field_map(string_tokenizer & header)
{
  bool match = true;

  for(field_map_type::const_iterator f = my_field_map.begin(); f != my_field_map.end(); ++f)
  {
    const field& current_field = f->second;

    if(current_field.type == trait)
    {
      bool exist_in_header = false;

      for(string_tokenizer::iterator h = header.begin(); h != header.end(); ++h)
      {
        if(toUpper(f->first) == toUpper(*h))
        {
          exist_in_header = true;

          break;
        }
      }

      if(!exist_in_header)
      {
        match = false;

        errors << priority(warning) << "Trait '" << current_field.name << "' doesn't exist in the pedigree file." << std::endl;
      }
    }
    else if(current_field.type == string_field)
    {
      bool exist_in_header = false;

      for(string_tokenizer::iterator h = header.begin() ; h != header.end(); ++h )
      {
        if(toUpper(f->first) == toUpper(*h))
        {
          exist_in_header = true;
          break;
        }
      }

      if( !exist_in_header )
      {
        match = false;

        errors << priority(warning) << "String '" << current_field.name << "' doesn't exist in the pedigree file." << std::endl;
      }
    }
  }

  return match;
}

/// Validates the format and file, then reads in the data.
///
bool
RefDelimitedPedigreeFile::input(RefMultiPedigree&  p,
                                const std::string& filename,
                                ostream&           messages)
{
  if(!validate_file(filename) || !validate_format(filename))
    return false;

  return RefPedigreeFile::input(p, filename, messages);
}

/// Verifies that the filename given exists and is valid.
///
bool
RefDelimitedPedigreeFile::validate_file(const std::string &filename)
{
  if(!filename.size())
  {
    errors << priority(fatal) << "No Family Data file specified." << std::endl;

    return false;
  }

  std::ifstream infile(filename.c_str());

  if(!infile.good())
  {
    errors << priority(fatal) << "Unable to open Family Data file '" << filename << "'. Please check your file." << std::endl;

    return false;
  }

  return true;
}

/// Checks to see if the format is valid, or, if it's in the file, reads it
/// from the file and checks that.  In either case, the format is then available
/// for use.
bool
RefDelimitedPedigreeFile::validate_format(const std::string &filename)
{
  // If the format should be in the file, we have to go and get it before we
  // can validate it.
  if(format_in_file())
  {
    // We assume that validate_file() has already been done
    std::ifstream infile(filename.c_str());

    std::string line;

    getline(infile, line);
    set_format(line);
  }
  
  // Check the format
  if(!format().size())
  {
    errors << priority(fatal) << "No format specified to read headerless Family Data file." << std::endl;

    return false;
  }

  return true;
}

//====================================================
//
//  input_pedigree(...)
//
//=====================================================
bool 
RefDelimitedPedigreeFile::input_pedigree(RefMultiPedigree& p,
                                         const string&     filename,
                                         ostream&          messages,
                                         bool              quiet)
{
  std::ifstream infile(filename.c_str());

  string_tokenizer tokenizer( format() );
  
  setup_tokenizer(tokenizer);

  std::string line;   // Current line

  // Skip the header if it's in the file.  validate_format() will have already
  // set it.
  if( format_in_file() )
  {
    getline(infile, line);
  }

  const RefMPedInfo &mped_info = p.info();

  if( !build_fields(tokenizer, mped_info, quiet) )
  {
    errors << priority(error) << "Cannot build list of fields to read.  Aborting." << std::endl;

    invalidate();
    return false;
  }

  if( !validate_fields(false, quiet) )
  {
    errors << priority(error) << "Invalid list of fields to read.  Aborting." << std::endl;

    invalidate();
    return false;
  }

  // Figure out if there's supposed to be a pedigree id field:
//  bool has_pedigree_id_field = false;
  
//  for(field_list_type::const_iterator itr = my_fields.begin(); itr != my_fields.end(); ++itr)
//  {
//    if(itr->type == pedigree_id)
//      has_pedigree_id_field = true;
//  }

  for( size_t count = 1; infile.good(); ++count )
  {
    getline(infile, line);

    if(!line.size())
      continue;

    std::string ped_name = "",                           // Pedigree id
                ind_name = "",                           // Individual id
                parent1  = "",                           // Parent1
                parent2  = "",                           // Parent2
                sex      = mped_info.sex_code_unknown(); // Sex code

    tokenizer.set_str(line);

    field_list_type::const_iterator  field_info     = my_fields.begin (),
                                     field_info_end = my_fields.end   ();
    string_tokenizer::const_iterator field          = tokenizer.begin (),
                                     field_end      = tokenizer.end   ();

    // Loop across tokens and field_info's at the same time, processing each:
    for( ; field_info != field_info_end && field != field_end; ++field_info, ++field)
    {
      switch( field_info->type )
      {
        case pedigree_id   : ped_name = *field;    break;
        case individual_id : ind_name = *field;    break;
        case parent_id     : if(!parent1.size())
                               parent1 = *field;
                             else
                               parent2 = *field;  break;
        case sex_code      : sex = *field;        break;
        default:                                  break;
      }

    } // End loop-across-tokens
    
    // Check for treat_ped_id options:
    if( !pedigree_id_count() )
      ped_name = "0";
    else if( get_treat_as_sibs() == true ) // There's a pedigree_id field and treat_as_sibs is enabled:
    {
      parent1 = ped_name + "_parent1";
      parent2 = ped_name + "_parent2";
      
      p.add_member(ped_name, parent1, MPED::SEX_MALE);
      p.add_member(ped_name, parent2, MPED::SEX_FEMALE);
    }

    // If parents are still empty for some reason, set them to be missing:
    parent1 = parent1.size() ? parent1 : mped_info.individual_missing_code();
    parent2 = parent2.size() ? parent2 : mped_info.individual_missing_code();
    
    // Skip this person if any essential bits of info are missing:
    if(!ped_name.size() && !ind_name.size() && !parent1.size() && !parent2.size() && !sex.size())
      continue;

    DEBUG_RPEDFILE(errors << priority(debug) << "Found (" << pn  << "," << id  << "," << sex << "," << p1 << "," << p2  << ")" << std::endl; )

    add_member(p, ped_name, ind_name, sex, parent1, parent2, count + format_in_file(), count);
  }

  infile.close();

  if( !build_pedigree(p) )
  {
    errors << priority(error)
           << "Fatal error building pedigree data structure.  Aborting." << std::endl;
    return false;    
  }

  return true;
}

bool 
RefDelimitedPedigreeFile::input_data(RefMultiPedigree& p,
                                     const string&     filename,
                                     ostream&          messages,
                                     bool              quiet)
{
  // Check to see if a second pass is necessary
  // (traits and covariates were defered in the first pass)
  //
  RefMPedInfo &mped_info = p.info();

  if( trait_count()            == 0 &&
      !sex_code_trait()             &&
      !pedigree_id_trait()          &&
      marker_count()           == 0 && 
      string_count()           == 0 &&
      mped_info.trait_count()  == 0 &&
      mped_info.marker_count() == 0 &&
      mped_info.string_count() == 0 )
    return true;

  // We assume validate_file() has been called already
  std::ifstream infile( filename.c_str());

  string_tokenizer tokenizer( format() );
  
  setup_tokenizer(tokenizer);

  string line;   // Current line

  // Skip header if it is in the file.  We assume that it's already been 
  // set by validate_format()
  if( format_in_file() )
  {
    getline(infile, line);
  }

  if( !build_fields(tokenizer, mped_info, quiet) )
  {
    errors << priority(error)
           << "Cannot build list of fields to read.  Aborting." << std::endl;
    invalidate();
    return false;
  }

  if( !validate_fields(true, quiet) )
  {
    errors << priority(error)
           << "Invalid list of fields to read.  Aborting." << std::endl;
    invalidate();
    return false;
  }

  if( !check_header_vs_field_map(tokenizer) )
  {
    errors << priority(error)
           << "Pedigree block specifies misspelled or missing pedigree data fields."
           << " Please see inf file for details."<< std::endl;
    invalidate();
    return false;
  }

  if( !build_data(p) )
  {
    errors << priority(error)
           << "Error building data indices.  Aborting." << std::endl;
    invalidate();
    return false;
  }

  field_list_type::const_iterator field_info      = my_fields.begin();
  field_list_type::const_iterator field_info_end  = my_fields.end();

  vector<pair<size_t,string> >  trait_values( trait_count() );
  vector<pair<size_t,string> > string_values( string_count() );

  // FIXME: This approach to reading markers will certainly have serious
  //        performance problems in the large scale and is here only for the
  //        simplicity.  A hash_map implementation would be much better.  A
  //        separate index vector is the optimal solution in the long term.
  typedef pair<string, string>          string_pair;
  typedef std::map<size_t, string_pair> marker_value_map;

  for( size_t count=1; infile.good(); ++count )
  {
    getline(infile, line);

    if(!line.size())
      continue;

    string pn;     // Pedigree id
    string id;     // Individual id
    string p1, p2; // Parents
    string sex;    // Sex code

    tokenizer.set_str(line);

    string_tokenizer::const_iterator field     = tokenizer.begin();
    string_tokenizer::const_iterator field_end = tokenizer.end();

    marker_value_map marker_values;
    marker_value_map marker_covariate_values;

    for( size_t t = 0; t < trait_count(); ++t )
      trait_values[t].first = (size_t)-1;

    for( size_t s = 0; s < string_count(); ++s )
      string_values[s].first = (size_t)-1;

    size_t tfound = 0;
    size_t sfound = 0;

    for( field_info = my_fields.begin() ;
         field_info != field_info_end && field != field_end ;
         ++field_info, ++field)
    {
      switch( field_info->type )
      {
        case   pedigree_id: pn = *field;    break;
        case individual_id: id = *field;    break;
        case         trait:
                            if( field_info->index >= mped_info.trait_count() )
                              break;
                            trait_values[tfound].first  = field_info->index;
                            trait_values[tfound].second = *field;
                            ++tfound;
                            break;

        case  string_field:
                            if( field_info->index >= mped_info.string_count() )
                              break;
                            string_values[sfound].first  = field_info->index;
                            string_values[sfound].second = strip_ws(*field);
                            ++sfound;
                            break;

        case        allele:
        case        marker:
        {
                            size_t index = field_info->index;
                            if( index >= mped_info.marker_count() )
                              break;

                            string_pair &values = marker_values[index];
                            if( !values.first.size())
                              values.first = strip_ws(*field);
                            else
                              values.second = strip_ws(*field);
                            break;
        }

        default:                            break;
      }
    }

    if( !pedigree_id_count() )
      pn = "0";

    if( !pn.size() || !id.size() )
      continue;

    if( !id.size() || id == mped_info.individual_missing_code() )
    {
      errors << priority(warning) << "[" << count + format_in_file()
             << "] Found an individual ID that is missing or"
             << " the same as the missing individual/parent code.  Skipping..."
             << std::endl;
      continue;
    }

    RefMultiPedigree::member_pointer mem = p.member_find(pn,id);

    // This should only happen when an error has already been reported
    if( !mem ) continue;

    RefMultiPedigree::pedinfo_type &info = mem->pedigree()->info();
    int ind_num = mem->index();

    // Make the traits
    for( size_t tt=0; tt < tfound; ++tt )
    {
      size_t t = trait_values[tt].first;
      const string &value = trait_values[tt].second;
      int code = info.set_trait(ind_num, t, value, mped_info.trait_info(t));

      switch(code)
      {
        case 0:                       // trait ok
        case 1:                       // trait ok, but missing
        case 4:                       // trait not set legitimately
                 break;
        case 3:                       // invalid ind. or trait id
                 errors << priority(warning) << "["<<count+format_in_file()
                        << "] Cannot set trait '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'." << std::endl;
                 break;
        case 2:                       // bad trait value
                 errors << priority(warning) << "["<<count+format_in_file()
                        << "] Unrecognized value for trait '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'.  Found '" << value << "'." << std::endl;

                 break;
       default:
                 errors << priority(error) << "[" <<count+format_in_file()
                        << "] Unexpected error setting trait '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'." << std::endl;
                 break;
      }
    }

    for( size_t ss=0; ss < sfound; ++ss )
    {
      size_t s = string_values[ss].first;
      const string &value = string_values[ss].second;
      bool code = info.set_string(ind_num, s, value);

      if(!code)
        errors << priority(warning) << "["<<count+format_in_file() 
	       << "] Cannot set string field '"
               << mped_info.string_info(s).name()
               << "' for individual '" << id << "' in pedigree '"
               << pn << "'." << std::endl;
    }

    marker_value_map::const_iterator mm;
    for( mm = marker_values.begin(); mm != marker_values.end(); ++mm )
    {
      size_t m = mm->first;
      const string_pair &values = mm->second;

      int code;
      code = info.set_phenotype(ind_num, mem->get_effective_sex(), m, values.first, values.second,
                                mped_info.marker_info(m));

      switch(code)
      {
        case 0:                       // marker ok
        case 1:                       // marker ok, but missing
        case 3:                       // invalid ind. or marker id <-- Shouldn't ever happen
                 break;
        case 2:                       // bad marker value
                 errors << priority(warning) << "["<<count+format_in_file()
                        << "] Unrecognized value for marker '"
                        << mped_info.marker_info(m).name()
                        << "' of individual '" << id << "' in pedigree '"
                        << pn << "': Found '" << values.first << "'";
                 if(values.second.size())
                   errors << ", '" << values.second << "'";
                 errors << ". Marker will be set to missing for this individual." << std::endl;
                 break;
         case 4:
                 errors << priority(warning) << "[" << count+format_in_file()
                        << "] Marker '" << mped_info.marker_info(m).name()
                        << "' is sex-dependent, but the individual '" << id 
                        << "' in pedigree '" << pn 
                        << "' has unknown sex.  Marker will be set to missing for this individual." << std::endl;
                 break;
         case 5:
                 errors << priority(warning) << "[" << count+format_in_file()
                        << "] Phenotype '" << values.first << "'";
                 if(values.second.size())
                   errors << ", '" << values.second << "'";
                 errors << " at sex-dependent marker '" << mped_info.marker_info(m).name()
                        << "' is inconsistent with the sex of individual '" << id
                        << "' in pedigree '" << pn 
                        << "'.  Marker will be set to missing for this individual." << std::endl;
                 break;
      }
    }

    for( mm = marker_covariate_values.begin(); mm != marker_covariate_values.end(); ++mm )
    {
      size_t t = mm->first;
      const string_pair &values = mm->second;

      string mcov_name = mped_info.trait_info(t).name();

      const string &value = get_marker_covariate_value(mcov_name, values.first, values.second);

      int code = info.set_trait(ind_num, t, value, mped_info.trait_info(t));

      switch(code)
      {
        case 0:                       // trait ok
        case 1:                       // trait ok, but missing
        case 4:                       // trait not set legitimately
                 break;
        case 3:                       // invalid ind. or trait id
                 errors << priority(warning) << "["<<count+format_in_file()
                        << "] Cannot set marker covariate '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'." << std::endl;
                 break;
        case 2:                       // bad trait value
                 errors << priority(warning) << "["<<count+format_in_file()
                        << "] Unrecognized value for marker covariate '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'.  Found '" << value << "'." << std::endl;

                 break;
       default:
                 errors << priority(error) << "[" <<count+format_in_file()
                        << "] Unexpected error setting marker covariate '"
                        << mped_info.trait_info(t).name()
                        << "' for individual '" << id << "' in pedigree '"
                        << pn << "'." << std::endl;
                 break;
      }
    }
  }

  return true;
}

bool
RefDelimitedPedigreeFile::output(RefMultiPedigree& p,
                                 const string&     filename,
                                 ostream&          messages,
                                 bool              quiet)
{
  if( !filename.size() || !delimiters().size()) return false;

  const RefMPedInfo &mped_info = p.info();

  if( format().size() )
  {
    string_tokenizer tokenizer( format() );
    tokenizer.set_whitespace( whitespace() );
    tokenizer.set_delimiters( delimiters() );
    tokenizer.set_skip_consecutive_delimiters( skip_consecutive_delimiters() );
    tokenizer.set_skip_leading_delimiters( skip_leading_delimiters() );
    tokenizer.set_skip_trailing_delimiters( skip_trailing_delimiters() );

    if(!build_fields(tokenizer, mped_info, quiet))
    {
      errors << priority(error)
             << "Cannot build list of fields to write.  Aborting..." << std::endl;
      invalidate();
      return false;
    }

    if(!validate_fields(false, quiet))
    {
      errors << priority(error)
             << "Invalid list of fields to write.  Aborting." << std::endl;
      invalidate();
      return false;
    }
  }

  if( !my_fields.size() )
  {
    errors << priority(error)
           << "No fields to write.  Aborting..." << std::endl;
    invalidate();
    return false;
  }

  std::ofstream outfile(filename.c_str());

  if(!outfile.good())
  {
    errors << priority(fatal) << "Unable to open output Family Data file '" << filename
           << "'." << std::endl;
    return false;
  }

  string delimiter = delimiters().substr(0,1);

  field_list_type::const_iterator field_info;
  field_list_type::const_iterator field_info_end  = my_fields.end();

  bool first_field = true;
  for( field_info = my_fields.begin();
       field_info != field_info_end ;
     ++field_info )
  {
    if( !first_field )
      outfile << delimiter;
    first_field = false;
    outfile << field_info->field_name;
  }
  outfile << std::endl;

  RefPedigree::member_type* mid;
  RefPedigree::member_type* pid;
  size_t parent_count;

  vector<size_t> markers_written( mped_info.marker_count() );

  RefMultiPedigree::pedigree_iterator ped;
  for(ped = p.pedigree_begin(); ped != p.pedigree_end(); ++ped)
  {
    const RefPedInfo &ped_info = ped->info();

    for(unsigned int i = 0; i < ped->member_count(); ++i)
    {
      mid = &ped->member_index(i);
      parent_count = 0;
      markers_written.clear();

      first_field = true;

      for(field_info = my_fields.begin(); field_info != field_info_end; ++field_info)
      {
        if( !first_field )
          outfile << delimiter;
        first_field = false;

        switch( field_info->type )
        {
          case pedigree_id:   outfile << ped->name();
                              break;

          case individual_id: outfile << mid->name();
                              break;

          case parent_id:
          {
                              string value = mped_info.individual_missing_code();
                              pid = NULL;
                              if(parent_count == 0 && mid->parent1())
                              {
                                pid = mid->parent1();
                                parent_count = 1;
                              }
                              if(!pid && parent_count < 2 && mid->parent2())
                              {
                                pid = mid->parent2();
                                parent_count = 2;
                              }
                              if(pid)
                                value = pid->name();
                              else
                                parent_count = 2;
                              outfile << value;
                              break;
          }
          case      sex_code:
                              if( mid->is_male() )
                                outfile << mped_info.sex_code_male();
                              else if( mid->is_female() )
                                outfile << mped_info.sex_code_female();
                              else
                                outfile << mped_info.sex_code_unknown();
                              break;

          case         trait:
          {
                              if( field_info->index >= mped_info.trait_count() )
                                break;

                              if( ped_info.trait_missing(i, field_info->index) )
                              {
                                outfile << mped_info.trait_info(
                                      field_info->index).string_missing_code();
                                break;
                              }

                              outfile << ped_info.trait(i, field_info->index);
                              break;
          }

          case  string_field:
          {
                              if( field_info->index >= mped_info.string_count() )
                                break;
                              outfile << ped_info.get_string(i, field_info->index);
                              break;
          }

          case        allele:
          {
                              if( field_info->index >= mped_info.marker_count() )
                                break;

                              if(markers_written[field_info->index] < 2)
                              {
                                const RefMarkerInfo& minfo =
                                   mped_info.marker_info(field_info->index);

                                uint pheno_id = ped_info.phenotype(i, field_info->index);

                                string pheno_name = minfo.get_phenotype(pheno_id).name();

                                string al1, al2;
                                MLOCUS::Ordering order;

                                minfo.gmodel().parse_genotype_name(pheno_name, al1, al2, order);

                                if(!al1.size()) al1 = minfo.missing_allele_name();
                                if(!al2.size()) al2 = minfo.missing_allele_name();

                                if(markers_written[field_info->index] == 0)
                                  outfile << al1;
                                else
                                  outfile << al2;
                                markers_written[field_info->index]++;
                              }
                              break;
          }

          case        marker:
          {
                              if( field_info->index >= mped_info.marker_count() )
                                break;
                              if(markers_written[field_info->index] == 0)
                              {
                                const RefMarkerInfo& minfo =
                                   mped_info.marker_info(field_info->index);

                                uint pheno = ped_info.phenotype(i, field_info->index);

                                outfile << minfo.get_phenotype(pheno).name();
                                markers_written[field_info->index] = 2;
                              }
                              break;
          }

          default:
                              break;
        }
      }
      outfile << std::endl;
    }
  }
  outfile.close();
  return true;
}

void 
RefDelimitedPedigreeFile::setup_tokenizer(string_tokenizer& tokenizer)
{
  tokenizer.set_whitespace( whitespace() );
  tokenizer.set_delimiters( delimiters() );
  tokenizer.set_skip_consecutive_delimiters( skip_consecutive_delimiters() );
  tokenizer.set_skip_leading_delimiters( skip_leading_delimiters() );
  tokenizer.set_skip_trailing_delimiters( skip_trailing_delimiters() );
}

//============================================================
//============================================================
//============================================================
//
//  RefLSFDelimitedPedigreeFile
//
//============================================================
//============================================================
//============================================================

RefLSFDelimitedPedigreeFile::RefLSFDelimitedPedigreeFile(cerrorstream &errors)
{
  set_error_sink(errors);

  // Set defaults
  set_delimiters(",\t");
  set_whitespace(" ");
  set_skip_trailing_delimiters(false);
  set_skip_leading_delimiters(false);
  set_skip_consecutive_delimiters(false);
}

bool
RefLSFDelimitedPedigreeFile::process_parameters(RefMPedInfo &mped_info, const LSFBase *params)
{
  return RefLSFPedigreeFile::process_parameters(mped_info, params);
}
  
bool 
RefLSFDelimitedPedigreeFile::process_parameter(RefMPedInfo &mped_info, const LSFBase *param)
{
  if(!param)
    return true;

  // Allow generic super-class to handle parameters
  RefLSFPedigreeFile::process_parameter(mped_info, param);

  string name = toUpper( param->name() );

       if(name == "FORMAT")                      return process_format2                     (mped_info, param);
  else if(name == "WHITESPACE")                  return process_whitespace                  (mped_info, param);
  else if(name == "DELIMITER_MODE")              return process_delimiter_mode              (mped_info, param);
  else if(name == "DELIMITERS")                  return process_delimiters                  (mped_info, param);
  else if(name == "SKIP_LEADING_DELIMITERS")     return process_skip_leading_delimiters     (mped_info, param);
  else if(name == "SKIP_TRAILING_DELIMITERS")    return process_skip_trailing_delimiters    (mped_info, param);
  else if(name == "SKIP_CONSECUTIVE_DELIMITERS") return process_skip_consecutive_delimiters (mped_info, param);
  else if(name == "STUDY_ID")                    return process_study_id                    (mped_info, param);
  else if(name == "PEDIGREE_ID")                 return process_pedigree_id                 (mped_info, param);
  else if(name == "TREAT_AS_SIBS")               { set_treat_as_sibs(true); return true; }
  else if(name == "INDIVIDUAL_ID")               return process_individual_id               (mped_info, param);
  else if(name == "PARENT_ID")                   return process_parent_id                   (mped_info, param);
  else if(name == "SEX_FIELD")                   return process_sex_field                   (mped_info, param);
  else if(name == "NO_SEX_FIELD")                { set_no_sex_field(true); return true; }
  else if(name == "MARKER" || 
          name == "ALLELE" ||
          name == "TRAIT_MARKER" )               return process_marker                      (mped_info, param);
  else if(name == "PHENOTYPE" ||
          name == "TRAIT"     ||
          name == "COVARIATE")                   return process_phenotype                   (mped_info, param);
  else if(name == "STRING")                      return process_string                      (mped_info, param);
  else if(name == "MARKER_LIST")                 return process_marker_list                 (mped_info, param);
  else if(name == "COVARIATE_LIST")              return process_covariate_list              (mped_info, param);

  else return true;
}

/// Validates the format and file, then reads in the data.
///
bool
RefLSFDelimitedPedigreeFile::input(RefMultiPedigree&  p,
                                   const std::string& filename,
                                   ostream&           messages)
{
  if(!validate_file(filename) || !validate_format(filename))
    return false;

  if(!build_marker_list_parameters(p.info())) return false;
  if(!build_covariate_list_parameters(p.info())) return false;
  
  return RefPedigreeFile::input(p, filename, messages);
}

/// Searches the list of MarkerListElement for one that begins with start.
///
/// \param start The start condition of the marker list
std::list<RefLSFDelimitedPedigreeFile::MarkerListElement>::iterator 
RefLSFDelimitedPedigreeFile::find_marker_list(string start)
{
  std::list<MarkerListElement>::iterator i = my_marker_lists.begin();
  
  for( ; i != my_marker_lists.end() && i->start_marker != start; ++i);
  
  return i;
}

bool
RefLSFDelimitedPedigreeFile::build_marker_list_parameters(RefMPedInfo &mped_info)
{
  string_tokenizer header(format());
  setup_tokenizer(header);
  
  for(string_tokenizer::iterator i = header.begin(); i != header.end(); ++i)
  {
    if(i->empty())
       continue;

    // Check to see if the field begins a marker list
    std::list<MarkerListElement>::iterator mlist_index = find_marker_list(toUpper(*i));
      
    if(mlist_index != my_marker_lists.end())
    {
      // If the field begins a marker_list, we can now add the fields to the map.
      if(!build_marker_list(header, i, mped_info, mlist_index))
        return false;
      
      my_marker_lists.erase(mlist_index);
    }
  }
  
  // Verify that all marker lists have been found.  If any remain, then there
  // is a problem!
  
  if(!my_marker_lists.empty())
  {
    for(std::list<MarkerListElement>::iterator i = my_marker_lists.begin();
        i != my_marker_lists.end(); ++i)
    {
      errors << priority(fatal) << "Marker list with start marker '"
             << i->start_marker << "' and end marker '"
             << i->end_marker << "' not found in pedigree file.   Please fix this and re-run S.A.G.E."
             << std::endl;
    }
    
    return false;
  }
  
  return true;
}

bool
RefLSFDelimitedPedigreeFile::build_marker_list(string_tokenizer&                      header,
                                               string_tokenizer::iterator             i, 
                                               RefMPedInfo&                           mped_info,
                                               std::list<MarkerListElement>::iterator mlist_index)
{
  MarkerListElement& mlist = *mlist_index;
  
  // Create the list of markers we must build by iterating through the 
  // header, looking for the end marker in the marker list.
  
  std::list<string> mlist_fields;
  
  for( ; i != header.end(); ++i)
  {
    // Check to make sure it's not also another parameter.  If it is, we've
    // got problems.
    field_map_type::const_iterator f = field_map().find( toUpper(*i) );
    
    if(f != field_map().end())
    {
      errors << priority(fatal) << "Overlap detected with marker_list with "
             << "field " << *i << ".  Marker list fields must be exclusive with "
             << "regards to other fields in the pedigree block.  Please fix this and re-run S.A.G.E."
             << std::endl;
      return false;
    }
    mlist_fields.push_back(toUpper(*i));
    
    if(toUpper(*i) == mlist.end_marker) break;
  }
  
  
  // If we didn't find the end, this is a fatal error
  if(i == header.end())
  {
    errors << priority(fatal) << "End marker '" << mlist.end_marker
           << "' of marker list starting with '" << mlist.start_marker 
           << "' not found in pedigree file fields following the start marker. "
           << "Please check your parameters."
           << std::endl;
           
    return false;
  }
  
  // Add the found markers to the marker list

  for(std::list<string>::const_iterator m = mlist_fields.begin(); m != mlist_fields.end(); ++m)
  {
    mlist.marker_params->attrs(true)->set(0, *m);
    if(!process_marker(mped_info, &*mlist.marker_params))
      return false;
  }
  return true;
}

std::list<RefLSFDelimitedPedigreeFile::MarkerListElement>::iterator 
RefLSFDelimitedPedigreeFile::find_covariate_list(string start)
{
  std::list<MarkerListElement>::iterator i = my_covariate_lists.begin();
  
  for( ; i != my_covariate_lists.end() && i->start_marker != start; ++i );
  
  return i;
}

bool
RefLSFDelimitedPedigreeFile::build_covariate_list_parameters(RefMPedInfo &mped_info)
{
  string_tokenizer header(format());
  setup_tokenizer(header);
  
  for( string_tokenizer::iterator i = header.begin(); i != header.end(); ++i )
  {
    if( i->empty() )
       continue;

    // Check to see if the field begins a covariate list
    std::list<MarkerListElement>::iterator clist_index = find_covariate_list(toUpper(*i));
      
    if( clist_index != my_covariate_lists.end() )
    {
      // If the field begins a covariate_list, we can now add the fields to the map.
      if( !build_covariate_list(header, i, mped_info, clist_index) )
        return false;
      
      my_covariate_lists.erase(clist_index);
    }
  }
  
  // Verify that all covariate lists have been found.  If any remain, then there
  // is a problem!
  
  if( !my_covariate_lists.empty() )
  {
    for( std::list<MarkerListElement>::iterator i = my_covariate_lists.begin();
        i != my_covariate_lists.end(); ++i )
    {
      errors << priority(fatal) << "Covariate list with start covariate '"
             << i->start_marker << "' and end covariate '"
             << i->end_marker << "' not found in pedigree file.   Please fix this and re-run S.A.G.E."
             << std::endl;
    }
    
    return false;
  }
  
  return true;
}

bool
RefLSFDelimitedPedigreeFile::build_covariate_list(string_tokenizer&                      header,
                                                  string_tokenizer::iterator             i, 
                                                  RefMPedInfo&                           mped_info,
                                                  std::list<MarkerListElement>::iterator clist_index)
{
  MarkerListElement& clist = *clist_index;
  
  // Create the list of covariates we must build by iterating through the 
  // header, looking for the end covariate in the covariate list.
  
  std::list<string> clist_fields;
  
  for( ; i != header.end(); ++i)
  {
    // Check to make sure it's not also another parameter.  If it is, we've
    // got problems.
    field_map_type::const_iterator f = field_map().find( toUpper(*i) );
    
    if( f != field_map().end() )
    {
      errors << priority(fatal) << "Overlap detected with covariate_list with "
             << "field " << *i << ".  Covariate list fields must be exclusive with "
             << "regards to other fields in the pedigree block.  Please fix this and re-run S.A.G.E."
             << std::endl;
      return false;
    }
    clist_fields.push_back(toUpper(*i));
    
    if( toUpper(*i) == clist.end_marker ) break;
  }
  
  
  // If we didn't find the end, this is a fatal error
  if( i == header.end() )
  {
    errors << priority(fatal) << "End covariate '" << clist.end_marker
           << "' of covariate list starting with '" << clist.start_marker 
           << "' not found in pedigree file fields following the start covariate. "
           << "Please check your parameters."
           << std::endl;
           
    return false;
  }
  
  // Add the found covariates to the covariate list

  for( std::list<string>::const_iterator m = clist_fields.begin(); m != clist_fields.end(); ++m )
  {
    clist.marker_params->attrs(true)->set(0, *m);
    if( !process_phenotype(mped_info, &*clist.marker_params) )
      return false;
  }
  return true;
}

/// Does LSFDelimited level processing of the "FORMAT" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_format2(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
    set_format_in_file(false);

  return true;
}

/// Processes the "WHITESPACE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_whitespace(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    set_whitespace(v.String());
    DEBUG_RPEDFILE(errors << priority(debug) << "Found whitespace: '" << whitespace() << "'" << std::endl;)
  }

  return true;
}

/// Processes the "DELIMITER_MODE" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_delimiter_mode(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);

  if( v.has_value() )
  {
    if( toUpper(v.String()) == "SINGLE")
    {
      set_skip_trailing_delimiters    (false);
      set_skip_leading_delimiters     (false);
      set_skip_consecutive_delimiters (false);
    }
    else if( toUpper(v.String()) == "MULTIPLE")
    {
      set_skip_trailing_delimiters    (true);
      set_skip_leading_delimiters     (true);
      set_skip_consecutive_delimiters (true);
    }
    else
    {
      errors << priority(error) << "Unknown value for parameter 'delimiter_mode'" << std::endl;
    }

    DEBUG_RPEDFILE(errors << priority(debug) << "Found column delimited flag." << std::endl;)
  }

  return true;
}

/// Processes the "DELIMITERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_delimiters(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    set_delimiters(v.String());
    DEBUG_RPEDFILE(errors << priority(debug) << "Found delimiters: " << delimiters() << std::endl;)
  }

  return true;
}

/// Processes the "SKIP_LEADING_DELIMITERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_skip_leading_delimiters(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    if( toUpper(v.String()) == "TRUE")
      set_skip_leading_delimiters(true);
    else if( toUpper(v.String()) == "FALSE")
      set_skip_leading_delimiters(false);
    else
      errors << priority(error) << "Unknown value for parameter 'skip_leading_delimiters'" << std::endl;
  }

  return true;
}

/// Processes the "SKIP_TRAILING_DELIMITERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_skip_trailing_delimiters(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    if( toUpper(v.String()) == "TRUE")
      set_skip_trailing_delimiters(true);
    else if( toUpper(v.String()) == "FALSE")
      set_skip_trailing_delimiters(false);
    else
      errors << priority(error) << "Unknown value for parameter 'skip_trailing_delimiters'" << std::endl;
  }

  return true;
}

/// Processes the "SKIP_CONSECUTIVE_DELIMITERS" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_skip_consecutive_delimiters(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    if( toUpper(v.String()) == "TRUE")
      set_skip_consecutive_delimiters(true);
    else if( toUpper(v.String()) == "FALSE")
      set_skip_consecutive_delimiters(false);
    else
      errors << priority(error) << "Unknown value for parameter 'skip_consecutive_delimiters'" << std::endl;
  }

  return true;
}

/// Processes the "STUDY_ID" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_study_id(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    add_study_id_field(v.String());
    DEBUG_RPEDFILE(errors << priority(debug) << "Found study id field = " << v.String() << std::endl;)
  }

  return true;
}
/// Processes the "PEDIGREE_ID" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_pedigree_id(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    add_pedigree_id_field(v.String());
    set_pedigree_id_name(v.String());
    
    DEBUG_RPEDFILE(errors << priority(debug) << "Found pedigree id field = " << v.String() << std::endl;)
  }

  return true;
}

/// Processes the "INDIVIDUAL_ID" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_individual_id(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    add_individual_id_field(v.String());
    DEBUG_RPEDFILE(errors << priority(debug) << "Found individual id field = " << v.String() << std::endl;)
  }

  return true;
}

/// Processes the "PARENT_ID" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_parent_id(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    add_parent_id_field(v.String());
    DEBUG_RPEDFILE(errors << priority(debug) << "Found parent id field = " << v.String() << std::endl;)
  }

  return true;
}

/// Processes the "SEX_FIELD" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_sex_field(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrVal v=attr_value(param, 0);
  if( v.has_value() )
  {
    add_sex_field(v.String());
    set_sex_field_name(v.String());

    DEBUG_RPEDFILE(errors << priority(debug) << "Found sex field = " << v.String() << std::endl;)

    RefLSFPedigreeFile::process_sex_code(mped_info, param);
  }

  return true;
}

/// Processes the "MARKER", "TRAIT_MARKER", and "ALLELE" parameters
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_marker(RefMPedInfo &mped_info, const LSFBase *param)
{
  string name = toUpper(param->name());
  string field_name;

  if( param->attrs() )
    field_name = param->attrs()->StringAttr(0);

  if(!field_name.size())
  {
    errors << priority(warning) << "Marker with no field name specified.  Skipping." << std::endl;
    return true;
  }
  
  string marker_name = param->attrs()->StringAttr("name");

  if(!marker_name.size())
    marker_name = field_name;

  // Test to see if the marker should be skipped, and return if so.
  if(test_skip_marker(param)) 
    return true;

  // Find the marker's model, if it exists
  size_t marker_id = mped_info.marker_find(marker_name);

  // Does it exist?
  bool has_model = marker_id < mped_info.marker_count();
  
  if(!has_model)
  {
    // We can only continue if dynamic parsing is allowed, so test that.
    if(!test_allow_dynamic(param))
    {
        errors << priority(warning) << "Marker '" << marker_name
               << "' not found.  Skipping." << std::endl;
        return true;
    }

    marker_id = setup_dynamic_marker(mped_info, param, marker_name);
  }
  
  parse_marker_parameters(mped_info, param, marker_id);

  // Finally, we have set up everything we need to set up for adding the marker.
  if(name == "ALLELE")
    add_allele_field( field_name, marker_name );
  else
    add_marker_field( field_name, marker_name );

  return true;
}

/// Processes the "TRAIT", "COVARIATE", and "PHENOTYPE" parameters
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_phenotype(RefMPedInfo &mped_info, const LSFBase *param)
{
  AttrList::const_iterator a;

  std::string name       = toUpper(param->name()),
              field_name = "";

  if( param->attrs() )
  {
    field_name = param->attrs()->StringAttr(0);
  }

  if(!field_name.size())
  {
    errors << priority(warning) << "Trait with no field name specified.  Skipping." << std::endl;
    return true;
  }

  string trait_name;

  trait_name = param->attrs()->StringAttr("name");

  if(!trait_name.size())
    trait_name = field_name;

  bool skip = skip_traits();

  if( !force_skip_traits() && has_attr(param, "skip") )
  {
    skip = true;

    if( toUpper(attr_value(param, "skip").String()) == "FALSE" )
    {
      skip = false;
    }
  }

  if(skip)
  {
    return true;
  }

  add_trait_field( field_name, trait_name );

  size_t t         = mped_info.trait_find(trait_name);
  bool   has_model = t < mped_info.trait_count();

  if(has_model)
  {
    return true;
  }

  bool trvar = param->attrs()->has_attr("variate")   || name == "TRAIT",
       cov   = param->attrs()->has_attr("covariate") || name == "COVARIATE";

  RefTraitInfo::trait_use use = RefTraitInfo::unknown_use;

  if(trvar && cov)
  {
    errors << priority(warning) << "Trait '" << trait_name
           << "' has multiple types listed.  Setting to unknown." << std::endl;
  }
  else if(trvar)
  {
    use = RefTraitInfo::trait_variate;
  }
  else if(cov)
  {
    use = RefTraitInfo::trait_covariate;
  }

  // Binary!
  if( param->attrs()->has_attr("binary"))
  {
    std::string affected   = (a = param->attrs()->find("affected"))   == param->attrs()->end() ? "1" : strip_ws(a->second.String()),
                unaffected = (a = param->attrs()->find("unaffected")) == param->attrs()->end() ? "0" : strip_ws(a->second.String()),
                missing    = (a = param->attrs()->find("missing"))    == param->attrs()->end() ? ""  : strip_ws(a->second.String());
    double      threshold  = (a = param->attrs()->find("threshold"))  == param->attrs()->end() ? std::numeric_limits<double>::quiet_NaN() : a->second.Real();
    
    if( affected == unaffected )
    {
      errors << priority(error) << "Pedigree File: Affected and unaffected codes "
             << "must be specified and may not be the same." << std::endl;
      skip = true;
    }
    else if( affected == missing || unaffected == missing )
    {
      errors << priority(error) << "Pedigree File Error: Affected and unaffected codes "
             << "must not be the same as the missing value code." << std::endl;
      skip = true;
    }
    else
    {
      DEBUG_RPEDFILE(errors << priority(debug) << "Found binary trait = " << v.String() << std::endl;)

      size_t t = mped_info.add_binary_trait( trait_name, use );

      mped_info.trait_info(t).set_string_missing_code              (missing);
      mped_info.trait_info(t).set_numeric_missing_code    (str2doub(missing));
      mped_info.trait_info(t).set_string_affected_code             (affected);
      mped_info.trait_info(t).set_string_unaffected_code           (unaffected);
      mped_info.trait_info(t).set_numeric_affected_code   (str2doub(affected));
      mped_info.trait_info(t).set_numeric_unaffected_code (str2doub(unaffected));
      mped_info.trait_info(t).set_threshold                        (threshold);
    }
  }

  // Doesn't have binary or categorical attribute; must be continuous!
  else
  {
    string missing;

    if( (a=param->attrs()->find("missing"))  != param->attrs()->end() )
      missing = strip_ws(a->second.String());

    DEBUG_RPEDFILE(errors << priority(debug) << "Found continuous trait = " << trait_name << std::endl;)

    size_t t = mped_info.add_continuous_trait( trait_name, use );

    mped_info.trait_info(t).set_string_missing_code(missing);
    mped_info.trait_info(t).set_numeric_missing_code(str2doub(missing));
  }
  if(skip)
  {
    DEBUG_RPEDFILE(errors << priority(debug) << "Skipping trait = " << trait_name << std::endl;)
    mped_info.add_trait( trait_name, RefTraitInfo::invalid_trait );
  }

  return true;
}

/// Processes the "STRING" parameter
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_string(RefMPedInfo &mped_info, const LSFBase *param)
{
  string field_name;

  if( param->attrs() )
    field_name = param->attrs()->StringAttr(0);

  if(!field_name.size())
  {
    errors << priority(warning) << "String field with no field name specified.  Skipping." << std::endl;
    return true;
  }

  string string_name = param->attrs()->StringAttr("name");
  if(!string_name.size())
    string_name = field_name;

  bool skip = false;
  if( has_attr(param, "skip") )
  {
    skip = true;
    if( toUpper(attr_value(param, "skip").String()) == "FALSE" )
      skip = false;
  }

  if(skip)
    return true;

  add_string_field(field_name, string_name);
  mped_info.add_string_field(string_name);

  return true;
}

/// Processes the "MARKER_LIST" parameter
///
/// Note that because the format is generally in the file, the
/// markers cannot be added to the set of fields at this time.  This
/// isn't done until the input() method is called.
///
/// \param mped_info The multipedigree info used for any checking.
/// \param param     The parameter we're processing
bool 
RefLSFDelimitedPedigreeFile::process_marker_list(RefMPedInfo& mped_info, const LSFBase* param)
{
  // Test for required attributes
  if(!param->attrs() || !has_attr(param,"start") || !has_attr(param,"end"))
  {
    errors << priority(fatal) << "marker_list parameter requires both start "
           << "and end attributes specified.  Please check your file and re-run S.A.G.E."
           << std::endl;
    return false;
  }
  
  // Make sure the required attributes have content
  string start_marker = toUpper(param->attrs()->StringAttr("start"));
  string end_marker = toUpper(param->attrs()->StringAttr("end"));

  if(start_marker.empty() || end_marker.empty())
  {
    errors << priority(fatal) << "marker_list parameter requires both start "
           << "and end attributes specified.  Please check your file and re-run S.A.G.E."
           << std::endl;
    return false;
  }

  // Make sure there isn't a marker list with this start already.
  //
  // Note:  This may not be necessary.  Marker lists which overlap are detected
  //        when their markers are inserted into the set of fields.  There
  //        may not be any reason to do this check.  But it's a better, and more
  //        informative error message, so it is left for that reason.
  if(find_marker_list(start_marker) != my_marker_lists.end())
  {
    errors << priority(fatal) << "Overlapping marker lists detected with duplicate first field "
           << start_marker << ".  Marker lists may not overlap.  Please correct this "
           << "and re-run S.A.G.E." << endl;
    return false;
  }
  
  // Create a new MarkerListElement for storing the list information until it can
  // be used.
  my_marker_lists.push_back(MarkerListElement());
  
  MarkerListElement& new_mlist = my_marker_lists.back();
  
  new_mlist.start_marker  = start_marker;
  new_mlist.end_marker    = end_marker;
  
  // We copy the attributes into a new LSFBase object so that they can be
  // used to set up the markers when it is time.
  new_mlist.marker_params = new LSFBase("MARKER_LIST");
  
  (*new_mlist.marker_params->attrs(true)) = (*param->attrs());

  return true;
}

bool 
RefLSFDelimitedPedigreeFile::process_covariate_list(RefMPedInfo& mped_info, const LSFBase* param)
{
  // Test for required attributes
  if( !param->attrs() || !has_attr(param,"start") || !has_attr(param,"end" ))
  {
    errors << priority(fatal) << "covariate_list parameter requires both start "
           << "and end attributes specified.  Please check your file and re-run S.A.G.E."
           << std::endl;
    return false;
  }
  
  // Make sure the required attributes have content
  string start_covariate = toUpper(param->attrs()->StringAttr("start"));
  string end_covariate = toUpper(param->attrs()->StringAttr("end"));

  if( start_covariate.empty() || end_covariate.empty() )
  {
    errors << priority(fatal) << "covariate_list parameter requires both start "
           << "and end attributes specified.  Please check your file and re-run S.A.G.E."
           << std::endl;
    return false;
  }

  // Make sure there isn't a covariate list with this start already.
  //
  // Note:  This may not be necessary.  Marker lists which overlap are detected
  //        when their markers are inserted into the set of fields.  There
  //        may not be any reason to do this check.  But it's a better, and more
  //        informative error message, so it is left for that reason.
  if( find_covariate_list(start_covariate) != my_covariate_lists.end() )
  {
    errors << priority(fatal) << "Overlapping covariate lists detected with duplicate first field "
           << start_covariate << ".  Covariate lists may not overlap.  Please correct this "
           << "and re-run S.A.G.E." << endl;
    return false;
  }
  
  // Create a new MarkerListElement for storing the list information until it can
  // be used.
  my_covariate_lists.push_back(MarkerListElement());
  
  MarkerListElement& new_clist = my_covariate_lists.back();
  
  new_clist.start_marker  = start_covariate;
  new_clist.end_marker    = end_covariate;
  
  // We copy the attributes into a new LSFBase object so that they can be
  // used to set up the covariates when it is time.
  new_clist.marker_params = new LSFBase("COVARIATE_LIST");
  
  (*new_clist.marker_params->attrs(true)) = (*param->attrs());

  return true;
}

/// Tests to see if a particular marker should be skipped.
///
/// The test is as follows:
///  -# If skip_markers() is \c true, we skip *unless* force_skip_markers() is
///     \c false and the parameter has a "skip=FALSE" attribute.
///  -# If skip_markers() is \c false, we skip *only* if force_skip_markers() is
///     \c false and the parameter has a "skip" attribute which is *not* "FALSE".
///
/// \param param The parameter containing the marker information
/// \returns \c true if the marker should be skipped, \c false otherwise
bool
RefLSFDelimitedPedigreeFile::test_skip_marker(const LSFBase* param) const
{
  // Get Default skip status
  bool skip = skip_markers();

  // Determine if skipping markers should be done
  if( !force_skip_markers() && has_attr(param, "skip") )
  {
    skip = true;
    if( toUpper(attr_value(param, "skip").String()) == "FALSE" )
      skip = false;
  }

  return skip;
}

/// Tests to see if a particular marker can be dynamic.
///
/// The test is as follows:
///  -# If dynamic_markers() is \c true, we allow dynamic *unless* 
///     force_dynamic_markers() is
///     \c false and the parameter has a "dynamic=FALSE" attribute.
///  -# If dynamic_markers() is \c false, we allow dynamic *only* if 
///     force_skip_markers() is
///     \c false and the parameter has a "dynamic" attribute which is *not* "FALSE".
///
/// \param param The parameter containing the marker information
/// \returns \c true if the marker is allowed to be dynamic, \c false otherwise
bool
RefLSFDelimitedPedigreeFile::test_allow_dynamic(const LSFBase* param) const
{
  // Get default dynamic status
  bool dynamic = dynamic_markers();

  if( !force_dynamic_markers() && has_attr(param, "dynamic") )
  {
    dynamic = true;
    if( toUpper(attr_value(param, "dynamic").String()) == "FALSE" )
      dynamic = false;
  }
  
  return dynamic;
}

/// Sets up a dynamic marker when one is parsed
///
/// \param param The parameter containing the dynamic marker
/// \returns The id of the new dynamic marker
size_t
RefLSFDelimitedPedigreeFile::setup_dynamic_marker
    (RefMPedInfo&   mped_info,
     const LSFBase* param,
     const string&  marker_name)
{
  size_t marker_id = mped_info.add_marker(marker_name);
  mped_info.marker_info(marker_id).gmodel().set_dynamic_alleles(true);

  const PhenotypeReaderInfo& pr = mped_info.get_pheno_reader_info();
  
  char   sep     = pr.get_allele_delimiter();
  string missing = pr.get_allele_missing();
  mped_info.marker_info(marker_id).gmodel().set_unphased_separator(sep);

  parse_dynamic_markers_missing(mped_info.marker_info(marker_id), missing);

  return marker_id;
}

/// This function parses the basic marker parameters which don't require a lot
/// of special setup.  This includes the x/y linked, delimiters, and missing
/// parameters.
void
RefLSFDelimitedPedigreeFile::parse_marker_parameters
    (RefMPedInfo &mped_info,
     const LSFBase* param,
     size_t marker_id)
{
  AttrList::const_iterator a;
  
  if(    (a=param->attrs()->find("allele_missing")) != param->attrs()->end()
      || (a=param->attrs()->find("missing"))        != param->attrs()->end() )
  {
    string global_missing = mped_info.get_pheno_reader_info().get_allele_missing();
    string local_missing  = strip_ws(a->second.String());

    if( global_missing != local_missing )
      update_marker_missing_info(mped_info.marker_info(marker_id), local_missing);
  }

  if( param->attrs()->has_attr("x_linked") )
  {
    mped_info.marker_info(marker_id).set_model_type(MLOCUS::X_LINKED);
    set_sex_linked_exist(true);
  }
  else if( param->attrs()->has_attr("y_linked") )
  {
    mped_info.marker_info(marker_id).set_model_type(MLOCUS::Y_LINKED);
    set_sex_linked_exist(true);
  }

  if(    (a=param->attrs()->find("allele_delimiter")) != param->attrs()->end()
      || (a=param->attrs()->find("delimiter"))        != param->attrs()->end() )
  {
    char global_sep = mped_info.get_pheno_reader_info().get_allele_delimiter();
    char local_sep  = global_sep;

    if( a->second.String().size() == 1 )
      local_sep = a->second.String()[0];
    else
      local_sep = strip_ws(a->second.String())[0];

    if( global_sep != local_sep )
      update_marker_delimiter_info(mped_info.marker_info(marker_id), local_sep);
  }
}

} // End namespace RPED
} // End namespace SAGE
