#ifndef _RPED_H
#define _RPED_H

//============================================================================
//  Reference Pedigree structure -- Definition of General Reference pedigree
//                                  storage classes
//
//  Author: Kevin Jacobs (jacobs@darwin.cwru.edu)
//
//  History:   0.1  kbj Initial implementation                   May 07 98
//                  yjs Added a new class PhenotypeReaderInfo    Jan 2003
//                       & new read_..() function in RefMpedInfo.
//                  djb Added a check to PedigreeSort().          8/1/3
//
//  Copyright (c) 1998  R.C. Elston
//    All Rights Reserved
//============================================================================

#include "sage/mlocus/mfile.h"

namespace SAGE {
namespace RPED {

class RefTraitInfo
{
  public:
    /// Identifies trait category
    enum trait_t   { unknown_trait,    //< Trait is of unknown type.
                     invalid_trait,    //< Trait is an invalid type.
                     continuous_trait, //< Trait is continuous (-INF to +INF, eg. 0.55).
                     binary_trait      //< Trait is binary (0/1, true/false, etc.).
                   }; 

    /// Identifies the intended use of the trait
    enum trait_use { unknown_use,     //< Trait's use is unknown.
                     trait_variate,   //< Trait is intended for use as a primary trait.
                     trait_covariate  //< Trait is intended for use as a covariate.
                   };

    RefTraitInfo();

    explicit RefTraitInfo(string    trait_name, 
                          trait_t   trait_type = continuous_trait,
                          trait_use usage      = unknown_use);

    void set_name(const string &n);
    void set_alias_name(const string &n);

    void set_type(trait_t typ);
    void set_usage(trait_use tru);

    void set_threshold(double th);
    void set_numeric_missing_code(double nm);
    void set_numeric_affected_code(double ac);
    void set_numeric_unaffected_code(double uc);

    void set_string_missing_code(const string &sm);
    void set_string_affected_code(const string &ac);
    void set_string_unaffected_code(const string &uc);

    const string & name() const;
    const string & alias_name() const;

    trait_t   type() const;
    trait_use usage() const;

    double threshold() const;
    double numeric_missing_code() const;
    double numeric_affected_code() const;
    double numeric_unaffected_code() const;

    const string &string_missing_code() const;
    const string &string_affected_code() const;
    const string &string_unaffected_code() const;

  private:

    void init();

    // Basic info:
    string    my_name;
    string    my_alias_name;
    trait_t   my_type;
    trait_use my_use;

    double  my_threshold;

    // Data members for missing value stuff:
    double  my_numeric_missing_code;
    string  my_string_missing_code;

    // Data members for binary trait stuff:
    double  my_numeric_affected_code;
    double  my_numeric_unaffected_code;
    string  my_string_affected_code;
    string  my_string_unaffected_code;
};

class RefStringInfo
{
  public:

      RefStringInfo();

      explicit RefStringInfo(string name);

      void set_name(const string &n);

      const string &name()  const;
    
      void set_string_missing_code(const string &sm);

      const string &string_missing_code() const;
    
  private:

    void init();

    string  my_name;
    string  my_string_missing_code;
};

typedef MLOCUS::inheritance_model RefMarkerInfo;

class PhenotypeReaderInfo
{
  public:

    PhenotypeReaderInfo();

    void set_allele_delimiter(char s);
    void set_allele_missing(const string s);

    char get_allele_delimiter() const;
    const string get_allele_missing() const;

  private:

    char       my_allele_delimiter;
    string     my_allele_missing;
};

class RefMPedInfo
{
  public:

    RefMPedInfo();

    void set_individual_missing_code(const string &imc);
    void set_sex_code_male(const string &c);
    void set_sex_code_female(const string &c);
    void set_sex_code_unknown(const string &c);

    const string &individual_missing_code() const;
    const string &sex_code_male() const;
    const string &sex_code_female() const;
    const string &sex_code_unknown() const;

    size_t add_trait(const string&                 trait_name, 
                           RefTraitInfo::trait_t   trait_type = RefTraitInfo::continuous_trait,
                           RefTraitInfo::trait_use use        = RefTraitInfo::unknown_use);

    size_t add_continuous_trait(const string&                 trait_name, 
                                      RefTraitInfo::trait_use use = RefTraitInfo::unknown_use);

    size_t add_binary_trait(const string &trait_name,
                   RefTraitInfo::trait_use use = RefTraitInfo::unknown_use);

    size_t add_continuous_covariate(const string &covariate_name);

    size_t add_binary_covariate(const string &covariate_name);

    RefTraitInfo & trait_info(size_t t);

    void remove_last_trait();
    void remove_trait_info(size_t t_id);

    size_t add_marker(const string &marker_name);
    size_t add_marker(const string &marker_name, const MLOCUS::inheritance_model& model);

    MLOCUS::inheritance_model_map & markers();

    RefMarkerInfo & marker_info(size_t m);
    RefMarkerInfo & marker_info(const string& m);

    void remove_marker_info(size_t m_id);
    void remove_marker_info(string m_name);

    size_t add_string_field(const string &name);
    RefStringInfo &string_info(size_t s);

    const RefTraitInfo &trait_info(size_t t) const;
    size_t trait_count() const;
    bool trait_exists(const string& trait_name) const;
    RefTraitInfo::trait_t get_trait_type(const string& trait_name) const;
    size_t trait_find(const std::string &name) const;

    size_t marker_count() const;
    const MLOCUS::inheritance_model_map& markers() const;
    const RefMarkerInfo &marker_info(size_t m) const;
    bool marker_exists(const string& marker_name) const;
    size_t marker_find(const std::string &name) const;

    size_t string_count() const;
    const RefStringInfo & string_info(size_t s) const;
    size_t string_find(const std::string &name) const;

    PhenotypeReaderInfo& get_pheno_reader_info();
    const PhenotypeReaderInfo& get_pheno_reader_info() const;

  private:

    string my_individual_missing_code;
    string my_sex_code_male;
    string my_sex_code_female;
    string my_sex_code_unknown;

    typedef vector<RefTraitInfo>     trait_vector;
    typedef vector<RefStringInfo>    string_vector;  
    typedef MLOCUS::inheritance_model_map    marker_map;

    trait_vector         my_trait_info;
    string_vector        my_string_info;
    marker_map           my_marker_info;
    PhenotypeReaderInfo  my_pheno_reader_info;  
};

class RefPedInfo
{
  public:

    RefPedInfo();

    void build(MPED::pedigree_base &ped);

    void resize_traits(size_t traits);
    void resize_markers(size_t markers, const RefMPedInfo& mped_info);
    void resize_strings(size_t s);
    void add_marker(const MLOCUS::inheritance_model& m);

    void remove_trait(size_t t_id);
    void remove_marker(size_t m_id);

    bool set_trait(size_t i, size_t t, double d);

    int set_trait(size_t i,  size_t t, 
                  const std::string  & value, 
                  RefTraitInfo & trait_info);

    bool set_phenotype(size_t i, size_t m, uint   p);

    int set_phenotype(size_t i,  size_t m, 
                      const std::string & allele1,
                      const std::string & allele2,
                      RefMarkerInfo & marker);

    int set_phenotype(size_t ind_id, 
                      MPED::SexCode ind_sex, 
                      size_t marker_id, 
                      const std::string & allele1,
                      const std::string & allele2,
                      RefMarkerInfo & marker);

    bool set_string(size_t i,  size_t s,  const string& val);

    void swap_members(size_t m1, size_t m2);

    size_t trait_count() const;
    size_t string_count() const;
    size_t marker_count() const;
    size_t member_count() const;

    double trait(size_t i, size_t t) const;
    bool trait_missing(size_t i, size_t t) const;

    uint phenotype(size_t i, size_t m) const;
    bool phenotype_missing(size_t i, size_t m, const RefMarkerInfo& mi) const;

    string get_string(size_t i, size_t s) const;

  private:

    int set_autosomal_phenotype(size_t i,  size_t m, const std::string & allele1,
                                const std::string & allele2, RefMarkerInfo & marker);
    int set_x_linked_phenotype(size_t i,
                               MPED::SexCode ind_sex, 
                               size_t m, const std::string & allele1,
                               const std::string & allele2, RefMarkerInfo & marker);
    int set_y_linked_phenotype(size_t i,
                               MPED::SexCode ind_sex, 
                               size_t m, const std::string & allele1,
                               const std::string & allele2, RefMarkerInfo & marker);
                      
  typedef vector<double> dvector;
  typedef vector<string> svector;
  typedef vector<uint>   uvector;

  vector<dvector> my_traits;
  vector<svector> my_strings;
  vector<uvector> my_markers;

  size_t my_member_count;
  size_t my_marker_count;
  size_t my_trait_count;
  size_t my_string_count;
};

typedef MPED::multipedigree <MPED::no_info, MPED::no_info, MPED::no_info, RefPedInfo, RefMPedInfo> RefMultiPedigree;

typedef RefMultiPedigree::pedigree_type     RefPedigree;
typedef RefMultiPedigree::subpedigree_type  RefSubpedigree;
typedef RefMultiPedigree::family_type       RefFamily;
typedef RefMultiPedigree::member_type       RefMember;

typedef MPED::multipedigree <MPED::no_info, MPED::no_info, MPED::no_info, RefPedInfo, RefMPedInfo> MultiPedigree;

typedef MultiPedigree::member_const_pointer        member_const_pointer;  

typedef MultiPedigree::pedigree_type     Pedigree;
typedef MultiPedigree::subpedigree_type  Subpedigree;
typedef MultiPedigree::family_type       Family;
typedef MultiPedigree::member_type       Member;

typedef MultiPedigree::member_pointer        MemberPointer;
typedef MultiPedigree::family_pointer        FamilyPointer;
typedef MultiPedigree::subpedigree_pointer   SubpedigreePointer;
typedef MultiPedigree::pedigree_pointer      PedigreePointer;
typedef MultiPedigree::multipedigree_pointer MultipedigreePointer;

typedef MultiPedigree::member_const_pointer        MemberConstPointer;
typedef MultiPedigree::family_const_pointer        FamilyConstPointer;
typedef MultiPedigree::subpedigree_const_pointer   SubpedigreeConstPointer;
typedef MultiPedigree::pedigree_const_pointer      PedigreeConstPointer;
typedef MultiPedigree::multipedigree_const_pointer MultipedigreeConstPointer;

typedef MultiPedigree::family_iterator        FamilyIterator;
typedef MultiPedigree::mate_iterator          MateIterator;
typedef MultiPedigree::member_iterator        MemberIterator;
typedef MultiPedigree::offspring_iterator     OffspringIterator;
typedef MultiPedigree::parent_iterator        ParentIterator;
typedef MultiPedigree::pedigree_iterator      PedigreeIterator;
typedef MultiPedigree::progeny_iterator       ProgenyIterator;
typedef MultiPedigree::sibling_iterator       SiblingIterator;
typedef MultiPedigree::subpedigree_iterator   SubpedigreeIterator;

typedef MultiPedigree::family_const_iterator        FamilyConstIterator;
typedef MultiPedigree::mate_const_iterator          MateConstIterator;
typedef MultiPedigree::member_const_iterator        MemberConstIterator;
typedef MultiPedigree::offspring_const_iterator     OffspringConstIterator;
typedef MultiPedigree::parent_const_iterator        ParentConstIterator;
typedef MultiPedigree::pedigree_const_iterator      PedigreeConstIterator;
typedef MultiPedigree::progeny_const_iterator       ProgenyConstIterator;
typedef MultiPedigree::sibling_const_iterator       SiblingConstIterator;
typedef MultiPedigree::subpedigree_const_iterator   SubpedigreeConstIterator;

typedef MultiPedigree RefMultiPedigree;

typedef RefMultiPedigree                               multipedigree_type;
typedef RefMultiPedigree::pedigree_type                pedigree_type;
typedef RefMultiPedigree::subpedigree_type             subpedigree_type;
typedef RefMultiPedigree::family_type                  family_type;
typedef RefMultiPedigree::member_type                  member_type;

typedef member_type                                    RefMember;
typedef pedigree_type                                  RefPedigree;

typedef RefMultiPedigree::family_pointer               family_pointer;
typedef RefMultiPedigree::member_pointer               member_pointer;
typedef RefMultiPedigree::pedigree_pointer             pedigree_pointer;
typedef RefMultiPedigree::subpedigree_pointer          subpedigree_pointer;

typedef RefMultiPedigree::family_const_pointer         family_const_pointer;
typedef RefMultiPedigree::member_const_pointer         member_const_pointer;
typedef RefMultiPedigree::pedigree_const_pointer       pedigree_const_pointer;
typedef RefMultiPedigree::subpedigree_const_pointer    subpedigree_const_pointer;

typedef RefMultiPedigree::family_iterator              family_iterator;
typedef RefMultiPedigree::mate_iterator                mate_iterator;
typedef RefMultiPedigree::member_iterator              member_iterator;
typedef RefMultiPedigree::offspring_iterator           offspring_iterator;
typedef RefMultiPedigree::parent_iterator              parent_iterator;
typedef RefMultiPedigree::pedigree_iterator            pedigree_iterator;
typedef RefMultiPedigree::progeny_iterator             progeny_iterator;
typedef RefMultiPedigree::sibling_iterator             sibling_iterator;
typedef RefMultiPedigree::subpedigree_iterator         subpedigree_iterator;

typedef RefMultiPedigree::family_const_iterator        family_const_iterator;
typedef RefMultiPedigree::mate_const_iterator          mate_const_iterator;
typedef RefMultiPedigree::member_const_iterator        member_const_iterator;
typedef RefMultiPedigree::offspring_const_iterator     offspring_const_iterator;
typedef RefMultiPedigree::parent_const_iterator        parent_const_iterator;
typedef RefMultiPedigree::pedigree_const_iterator      pedigree_const_iterator;
typedef RefMultiPedigree::progeny_const_iterator       progeny_const_iterator;
typedef RefMultiPedigree::sibling_const_iterator       sibling_const_iterator;
typedef RefMultiPedigree::subpedigree_const_iterator   subpedigree_const_iterator;

//==================
// helper functions
//==================

void PedigreeSort(RefPedigree &p);

MemberConstPointer ancestor_error(std::set<size_t>& ancestors, const MemberConstPointer ind);

void print_mped(ostream &o, const RefMultiPedigree& p);

} // End namespace RPED
} // End namespace SAGE

#include "sage/rped/rped.ipp"

#endif
