// rped.cpp -- Implements the PedigreeSort function
#include "sage/rped/rped.h"

#ifdef _WIN32
#	ifdef max
#		undef max
#		undef min
#	endif
#	define finite std::isfinite
#endif

namespace SAGE{
namespace RPED {

#define DEBUG_PEDSORT(x)

RefMPedInfo::RefMPedInfo()
{
  set_sex_code_male           ("M");
  set_sex_code_female         ("F");
  set_sex_code_unknown        ("");
  set_individual_missing_code ("");

  my_trait_info.clear();
}

size_t 
RefMPedInfo::trait_find(const std::string &name) const
{
  std::string tn = toUpper(name);

  for( size_t t = 0; t < trait_count(); ++t )
    if( toUpper(trait_info(t).name()) == tn || toUpper(trait_info(t).alias_name()) == tn )
      return t;

  return (size_t)-1;
}

size_t 
RefMPedInfo::add_trait(const string &trait_name, RefTraitInfo::trait_t trait_type, RefTraitInfo::trait_use use)
{
  if(trait_exists(trait_name))
  {
    return trait_find(trait_name);
  }
  else
  {
    my_trait_info.push_back(RefTraitInfo(trait_name, trait_type, use));

    return trait_count() - 1;
  }
}


void
RefMPedInfo::remove_trait_info(size_t t_id)
{
  assert(t_id < my_trait_info.size());
  my_trait_info.erase(my_trait_info.begin() + t_id);
}

void
RefMPedInfo::remove_marker_info(size_t m_id)
{
  if( my_marker_info.name_find(my_marker_info.name(m_id)) != my_marker_info.name_end() )
    my_marker_info.erase(my_marker_info.name(m_id));
}

void
RefMPedInfo::remove_marker_info(string m_name)
{
  if( my_marker_info.name_find(m_name) != my_marker_info.name_end() )
    my_marker_info.erase(m_name);
}

void
RefPedInfo::remove_trait(size_t t_id)
{
  assert(t_id < my_traits.size());
  my_traits.erase(my_traits.begin() + t_id);

  --my_trait_count;
}

void
RefPedInfo::remove_marker(size_t m_id)
{
  assert(m_id < my_markers.size());
  my_markers[m_id] = my_markers[my_markers.size()-1];

  my_markers.pop_back();

  --my_marker_count;
}


// Returns: 0 - trait set ok
//          1 - trait set ok, but missing
//          2 - bad trait value, assumed missing
//          3 - invalid invididual or trait id
//          4 - trait not set due to trait settings
int
RefPedInfo::set_trait(size_t i, size_t t, const std::string & value, RefTraitInfo & trait_info)
{
  // Set up return value:
  int code = 0;

  // Trait # or member # out of range:
  if(t >= trait_count() || i >= member_count())
  {
    code = 3;
  }

  // Invalid trait type:
  else if(trait_info.type() == RefTraitInfo::invalid_trait)
  {
    code = 4;
  }

  // Continuous or binary trait
  else
  {
    double d = str2doub(value);

    if(!finite(d))
      code = 2;

    std::string smiss  = trait_info.string_missing_code  ();
    double      nmiss  = trait_info.numeric_missing_code ();

    if( value == smiss || (finite(nmiss) && d == nmiss))
    {
      code = 1;
      d    = numeric_limits<double>::quiet_NaN();
    }
    else if( trait_info.type() == RefTraitInfo::binary_trait )
    {
      double thresh = trait_info.threshold();
        
      if(value == trait_info.string_affected_code() )
      {    
        code = 0;
        d    = 1.0;
      }
      else if( value == trait_info.string_unaffected_code() )
      {
        code = 0;
        d    = 0.0;
      }
      else if( finite(d) && d == trait_info.numeric_affected_code() )
      {    
        code = 0;
        d = 1.0;
      }
      else if( finite(d) && d == trait_info.numeric_unaffected_code() )
      {
        code = 0;
        d = 0.0;
      }
      else if( finite(d) && finite(thresh) )
      {
        code = 0;
        d    = (d > thresh) ? 1.0 : 0.0;
      }
      else
      {
        code = 2;
        d = numeric_limits<double>::quiet_NaN();
      }
    }

    if(!set_trait(i, t, d))
    {
      code = 3;
    }
  }


  return code;
}

inline
std::string compose_genotype(const std::string& value1,
                             const std::string& value2,
                             RefMarkerInfo&     marker_info)
{
  std::string value = value1;

  if(value2.size())
    value  += marker_info.gmodel().unphased_separator() + value2;

  // Add any genotype information we need to parse the genotype.  Does
  // nothing if inappropriate
  marker_info.add_genotype_dynamically(value);
  
  // Convert to canonical form, stripping internal white space, etc
  string canonical_value = marker_info.gmodel().parse_genotype_name(value);
  
  return canonical_value.size() ? canonical_value : value;
}

inline bool phenotype_not_found(uint phenotype)
{
  return phenotype == MLOCUS::NPOS;
}
inline bool allele_invalid(MLOCUS::allele a1)
{
  return !a1.is_valid();
}

/// Given an allele string, adds the allele to the marker if it isn't already
/// present with an allele frequency of 0.0
inline void add_marker_allele(RefMarkerInfo& marker_info, std::string allele)
{
  const MLOCUS::genotype_model& gmodel = marker_info.gmodel();
  
  if( allele.size()                          &&
      allele != gmodel.missing_allele_name() && 
      allele_invalid(gmodel.get_allele(allele)) )
  {
    marker_info.add_allele(allele, 0.0, true, true);
  }
}

/// Returns \c true if the marker allows expansion alleles, \c false otherwise.
///
inline bool marker_is_expandable(RefMarkerInfo& marker_info)
{
  return marker_info.gmodel().dynamic_alleles() && marker_info.codominant();
}

/// If the marker is expandable (allows new alleles), adds alleles which may be
/// missing
inline void expand_marker(RefMarkerInfo& marker_info, std::string phenotype)
{
    uint pid = marker_info.get_phenotype_id(phenotype);
    
    if( phenotype_not_found(pid) && marker_is_expandable(marker_info))
    {
      string   allele1;
      string   allele2;
      MLOCUS::Ordering order;
 
      marker_info.gmodel().parse_genotype_name(phenotype, allele1, allele2, order);

      add_marker_allele(marker_info, allele1);
      add_marker_allele(marker_info, allele2);
    }
}
// Returns: 0 - marker set ok
//          1 - marker set ok, but missing
//          2 - bad marker value, assumed missing
//          3 - invalid invididual or marker id
int
RefPedInfo::set_autosomal_phenotype(size_t i, size_t m, const std::string& value1,
                        const std::string& value2, RefMarkerInfo& marker_info)
{
  int code = 0;

  std::string value = compose_genotype(value1, value2, marker_info);

  uint p = marker_info.get_phenotype_id(value);
  
  if( value.empty() || p == marker_info.get_missing_phenotype_id())                                  // If value empty, marker considered missing.
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());
    
    return 1;
  }

  if( phenotype_not_found(p) || !set_phenotype(i, m, p))
  {
    code = 2;
  }

  return code;
}

// Returns: 0 - marker set ok
//          1 - marker set ok, but missing
//          2 - bad marker value, assumed missing
//          3 - invalid invididual or marker id
//          4 - sex based marker, but individual missing sex
//          5 - phenotype incompatible with sex of individual
int
RefPedInfo::set_x_linked_phenotype
  (size_t i, MPED::SexCode s, size_t m, const std::string& value1,
   const std::string& value2, RefMarkerInfo& marker_info)
{
  std::string value = compose_genotype(value1, value2, marker_info);
  
  if(value.empty()) // If value empty, marker considered missing.
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());
    
    return 1;
  }
  
  // Look up the phenotype id, if we can
  uint p = MLOCUS::NPOS;
  
  if(MPED::is_male(s))
    p = marker_info.get_phenotype_id(value + "(male)");
  
  if(p == MLOCUS::NPOS)
    p = marker_info.get_phenotype_id(value);

  if(p == marker_info.get_missing_phenotype_id())
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());
    
    return 1;
  }

  if(MPED::is_sex_unknown(s))
    return 4;
  
  if( phenotype_not_found(p) || !set_phenotype(i, m, p))
  {
    return 2;
  }

  // Check validity of the phenotype given the sex
  bool is_valid = false;
  
  for(RefMarkerInfo::phased_penetrance_iterator it = marker_info.phased_penetrance_begin(p);
      !is_valid && it != marker_info.phased_penetrance_end(p); ++it)
  {
    if((MPED::is_male(s)   && it.phased_geno().is_male_compatible()) ||
       (MPED::is_female(s) && it.phased_geno().is_female_compatible()) )
      is_valid = true;
  }
  for(RefMarkerInfo::unphased_penetrance_iterator it = marker_info.unphased_penetrance_begin(p);
      !is_valid && it != marker_info.unphased_penetrance_end(p); ++it)
  {
    if((MPED::is_male(s)   && it.unphased_geno().is_male_compatible()) ||
       (MPED::is_female(s) && it.unphased_geno().is_female_compatible()) )
      is_valid = true;
  }
  
  if(!is_valid)
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());
    
    return 5;
  }
  
  return 0;
}

// Returns: 0 - marker set ok
//          1 - marker set ok, but missing
//          2 - bad marker value, assumed missing
//          3 - invalid invididual or marker id
//          4 - sex based marker, but individual missing sex
//          5 - phenotype incompatible with sex of individual

int
RefPedInfo::set_y_linked_phenotype
  (size_t i, MPED::SexCode s, size_t m, const std::string& value1,
   const std::string& value2, RefMarkerInfo& marker_info)
{
  std::string value = compose_genotype(value1, value2, marker_info);

  if(MPED::is_female(s))
  {
    uint female_id = marker_info.get_phenotype_id("~X/~X");

    set_phenotype(i, m, female_id);

    if(value.empty()) // If value empty, marker considered missing.
    {
      return 1;
    }
    
    uint p = marker_info.get_phenotype_id(value);

    if(p == marker_info.get_missing_phenotype_id())
    {
      return 1;
    }
    
    if(phenotype_not_found(p))
      return 2;
    
    if(p != female_id)
    {
      return 5;
    }
    return 0;
  }
  
  // Male
  uint p = marker_info.get_phenotype_id(value);
  
  if(value.empty()) // If value empty, marker considered missing.
  {
    return 1;
  }

  if(p == marker_info.get_missing_phenotype_id())
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());
    
    return 1;
  }

  if(MPED::is_sex_unknown(s))
    return 4;
  
  if( phenotype_not_found(p) || !set_phenotype(i, m, p))
  {
    return 2;
  }

  if(p == marker_info.get_phenotype_id("~X/~X"))
  {
    set_phenotype(i, m, marker_info.get_missing_phenotype_id());

    return 5;
  }
  
  return 0;
}


// Returns: 0 - marker set ok
//          1 - marker set ok, but missing
//          2 - bad marker value, assumed missing
//          3 - invalid invididual or marker id
//          4 - sex based marker, but individual missing sex
//          5 - phenotype incompatible with sex of individual
int
RefPedInfo::set_phenotype(size_t i, size_t m, const std::string& value1,
                          const std::string& value2, RefMarkerInfo& marker_info)
{
  if( m >= marker_count() || i >= member_count() )
    return 3;

  return set_autosomal_phenotype(i,m,value1,value2, marker_info);
}

// Added for X & Y-linkage  - yjs  Mar. 2002
// Adjust the data according to the sex for X-linked & Y-linked marker.
//
//   X-linked marker:
//     female - no adjustment needed, missing allele is missing allele.
//       male - if one is missing, the missing allele is replaced with pseudo-allele "~Y".
//
//   Y-linked marker:
//     female - if allele exist, should be reported as error.
//              if one missing, replaced with "~X" to be reported as error.
//       male - if one missing, the missing allele is replaced with pseudo-allele "~X".
//
// Returns: 0 - marker set ok
//          1 - marker set ok, but missing
//          2 - bad marker value, assumed missing
//          3 - invalid invididual or marker id
//          4 - sex based marker, but individual missing sex
//          5 - phenotype incompatible with sex of individual

int
RefPedInfo::set_phenotype(size_t i, MPED::SexCode s, size_t m, const std::string& value1,
                          const std::string& value2, RefMarkerInfo& marker_info)
{
  if( m >= marker_count() || i >= member_count() )
    return 3;
  
  switch(marker_info.get_model_type())
  {
    case MLOCUS::AUTOSOMAL : return set_autosomal_phenotype(i,m,value1,value2, marker_info);
    case MLOCUS::X_LINKED  : return set_x_linked_phenotype(i,s,m,value1,value2, marker_info);
    case MLOCUS::Y_LINKED  : return set_y_linked_phenotype(i,s,m,value1,value2, marker_info);
  }
  
  return 3;
}

// Sort a reference pedigree into topological order based on lineage
void PedigreeSort(RefPedigree& p)
{
  int p1, p2, sp1, sp2;
  const RefPedigree::member_type *mid;

  //RefPedInfo& info = p.info();

  for(size_t count = 0; count < p.member_count(); ++count)
  {
    mid = &p.member_index(count);
    p1 = p2 = -1;
    if(mid->parent1())
      p1 = mid->parent1()->index();
    if(mid->parent2())
      p2 = mid->parent2()->index();

    /*
    if(p1 == mid->index() || p2 == mid->index() )
    {
      cerr << "Your pedigree is broken.  Please fix it." << endl;
      cerr << "Indiviual " << mid->name() << " in pedigree " 
           << mid->pedigree()->name() << " is his/her own "
           << "parent (" << mid->parent1()->name() << "," 
                         << mid->parent2()->name() << ")" << endl;
      exit(1);
    }
    */
    
    // - Generalization of above check to insure that a person is not his/her own
    //   ancestor.  -djb  8/1/3
    //
    set<size_t>  ancestors;
    const RefPedigree::member_type*  member = ancestor_error(ancestors, mid);
    if(member)
    {
      cerr << "Problem detected in pedigree '" << mid->pedigree()->name() << "'."
           << "  One or more members is his/her own ancestor.\n"
           << "Member '" << member->name() << "' is involved.  Please "
           << "fix this problem and rerun program." << endl;
      exit(1);
    }

    DEBUG_PEDSORT(std::cout << "i=" << mid->index() << ", p=(" << p1 << "," << p2 << ")" << std::endl;)

    if(p1 == -1 && p2 == -1) continue;

    size_t maxp = max(p1,p2);

    if(mid->index() < maxp)
    {
      //info.swap_members(count, maxp);
      p.member_index_swap(count, maxp);
      --count;
    }
  }

  // Need to sort subpedigrees too.
  RefMultiPedigree::subpedigree_iterator si = p.subpedigree_begin();
  for( ; si != p.subpedigree_end(); ++si )
  {
    RefMultiPedigree::member_iterator mi = si->member_begin();
    for( ; mi != si->member_end(); ++mi )
    {
      sp1 = sp2 = -1;
      if( mi->parent1() )
        sp1 = mi->parent1()->subindex();
      if( mi->parent2() )
        sp2 = mi->parent2()->subindex();

      if( sp1 == -1 && sp2 == -1 ) continue;

      size_t maxsp = max(sp1,sp2);

      if( mi->subindex() < maxsp )
      {
        mi->subpedigree()->member_index_swap(mi->subindex(), maxsp);
        if (mi != si->member_begin()) --mi;
      }
    }
  }
}


// - If individual is in ancestor set, return pointer to that individual 
//   otherwise, return 0.  Check ancestry via parent 1.  Check ancestry via parent 2.
//
const RefPedigree::member_type*
ancestor_error(std::set<size_t>& ancestors, const RefPedigree::member_type* ind)
{
  if(! ind)
  {
    return  0;
  }
  
  if(ancestors.find(ind->index()) != ancestors.end())
  {
    return  ind;
  }
  
  ancestors.insert(ind->index());
  const RefPedigree::member_type*  p1_error = ancestor_error(ancestors, ind->parent1());
  
  if(p1_error)
  {
    ancestors.erase(ind->index());
    return  p1_error;
  }
  else
  {
    const RefPedigree::member_type*  p2_error = ancestor_error(ancestors, ind->parent2());
    ancestors.erase(ind->index());  
    return  p2_error;
  }
}

void print_mped(ostream &o, const RefMultiPedigree& mped)
{
  o << std::endl
    << "Family structure information :" //on the first "
    //<< "10"
    //<< " individuals"
    << std::endl
    << std::endl
    << "     PED. ID       IND. ID       SEX       PARENT1       PARENT2     " << std::endl
    << "     ------------  ------------  --------  ------------  ------------" << std::endl;

  for( size_t i = 0; i < mped.member_count(); ++i )
  {
    const member_type& mem = mped.member_index(i);

    string pn = mem.pedigree()->name();
    string id = mem.name();
    string sex = "";
    if( mem.is_male() ) sex = mped.info().sex_code_male();
    else if( mem.is_female() ) sex = mped.info().sex_code_female();
    string p1 = mem.parent1() ? mem.parent1()->name() : mped.info().individual_missing_code();
    string p2 = mem.parent2() ? mem.parent2()->name() : mped.info().individual_missing_code();

    o << left
      << "     " << std::setw(12) << pn
      << "  "    << std::setw(12) << id
      << "  "    << std::setw(8)  << sex
      << "  "    << std::setw(12) << p1
      << "  "    << std::setw(12) << p2
      << std::endl;
  }
  o << std::endl;

  if( mped.info().trait_count() > 0 )
  {
    o << std::endl
      << "Phenotypes :" // for the first "
      //<< "10"
      //<< " individuals"
      << std::endl
      << std::endl
      << "     PED. ID       IND. ID     " << left;

    for( size_t t = 0; t < mped.info().trait_count(); ++t )
    {
      o << "  " << std::setw(20) << mped.info().trait_info(t).name();
    }
    o << std::endl;

    o << "     ------------  ------------";

    for( size_t t = 0; t < mped.info().trait_count(); ++t )
    {
      o << "  --------------------";
    }
    o << std::endl;

    for( size_t i = 0; i < mped.member_count(); ++i )
    {
      const member_type& mem = mped.member_index(i);

      string pn = mem.pedigree()->name();
      string id = mem.name();

      const RefPedInfo& ped_info = mem.pedigree()->info();

      o << "     " << std::setw(12) << pn << "  " << std::setw(12) << id;

      for( size_t t = 0; t < mped.info().trait_count(); ++t )
      {
        if( SAGE::isnan(ped_info.trait(mem.index(), t)) )
          o << "  " << std::setw(20) << mped.info().trait_info(t).string_missing_code();
        else
          o << "  " << std::setw(20) << ped_info.trait(mem.index(), t);
      }
      o << std::endl;
    }
    o << std::endl;
  }

  if( mped.info().marker_count() > 0 )
  {
    o << std::endl
      << "Markers :" //for the first "
      //<< "10"
      //<< " individuals"
      << std::endl
      << std::endl
      << "     PED. ID       IND. ID     " << left;
    size_t m_cnt = 0;
    for( size_t m = 0; m < mped.info().marker_count() && m_cnt < 10; ++m )
    {
      const RefMarkerInfo& minfo = mped.info().marker_info(m);
      if( minfo.name().find("~~~") != std::string::npos )
        continue;

      o << "  " << std::setw(12) << mped.info().marker_info(m).name();
      ++m_cnt;
    }
    o << std::endl;

    o << "     ------------  ------------";

    m_cnt = 0;
    for( size_t m = 0; m < mped.info().marker_count() && m_cnt < 10; ++m )
    {
      const RefMarkerInfo& minfo = mped.info().marker_info(m);
      if( minfo.name().find("~~~") != std::string::npos )
        continue;

      o << "  ------------";
      ++m_cnt;
    }
    o << std::endl;

    for( size_t i = 0; i < mped.member_count(); ++i )
    {
      const member_type& mem = mped.member_index(i);

      string pn = mem.pedigree()->name();
      string id = mem.name();

      const RefPedInfo& ped_info = mem.pedigree()->info();

      o << "     " << std::setw(12) << pn << "  " << std::setw(12) << id;

      m_cnt = 0;
      for( size_t m = 0; m < mped.info().marker_count() && m_cnt < 10; ++m )
      {
        const RefMarkerInfo& minfo = mped.info().marker_info(m);

        if( minfo.name().find("~~~") != std::string::npos )
          continue;

        uint pheno_id = minfo.get_missing_phenotype_id();

        if( !ped_info.phenotype_missing(mem.index(), m, minfo) )
          pheno_id = ped_info.phenotype(mem.index(), m);

        string pheno_name = minfo.get_phenotype(pheno_id).name();

        if( pheno_id == minfo.get_missing_phenotype_id() )
          pheno_name = "?" + string(1, minfo.gmodel().unphased_separator()) + "?";
        else if( minfo.codominant(pheno_id) )
        {
          string al1, al2;
          MLOCUS::Ordering morder;
          string geno_name = minfo.unphased_penetrance_begin(pheno_id).unphased_geno().name();
          minfo.gmodel().parse_genotype_name(geno_name, al1, al2, morder);
          pheno_name = al1 + minfo.gmodel().unphased_separator() + al2;
        }
        //else
          //pheno_name = long2str(pheno_id, 2);

        o << "  " << std::setw(12) << pheno_name;
        ++m_cnt;
      }
      o << std::endl;
    }
    o << std::endl;
#if 0
    for( size_t m = 0; m < mped.info().marker_count(); ++m )
    {
      const RefMarkerInfo& minfo = mped.info().marker_info(m);
      MLOCUS::inheritance_model_test(o, minfo, minfo.name().c_str());
      MLOCUS::inheritance_model_print(o, minfo, minfo.name().c_str());
    }
#endif
  }
  o << std::endl;

  return;
}

} // end namespace RPED
} // end namespace SAGE

