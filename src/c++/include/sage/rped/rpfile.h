#ifndef REF_PED_FILE_H
#define REF_PED_FILE_H

#include "sage/global/LSFfile.h"
#include "sage/rped/rped.h"

namespace SAGE {
namespace RPED {

class RefPedigreeFile
{
  public:

  /// Uniquely identifies the type of a field.
    enum field_t 
    { 
      skip,          // Skipped over
      study_id,      // Indicates the study id number
      pedigree_id,   // Indicates the pedigree id number
      individual_id, // Indicates the individual id number
      parent_id,     // Indicates a parent id number
      sex_code,      // Indicates an individual's sex
      trait,         // Indicates a trait-type field
      marker,        // Indicates a marker
      allele,        // Indicates an allele
      string_field   // Indicates a string-type field
    };

    struct field
    {
      field(field_t t = skip, const std::string & f = "", const std::string & n  ="",  size_t i = (size_t)-1);

      field_t type;
      
      std::string field_name;
      std::string name;
      
      size_t index;

      void invalidate();
    };

    /// List of fields
    typedef std::list<field> field_list_type;

    /// Map of fields (mapping field name to field instance)
    typedef std::map<std::string, field> field_map_type;

    RefPedigreeFile(cerrorstream &err = sage_cerr);

    virtual ~RefPedigreeFile();

    void set_error_sink(cerrorstream &err);
    void set_verbose_output(size_t v = (size_t)-1);
    void set_reject_partial_lineage(bool r = true);
    void set_require_record(bool a = true);
    void set_skip_traits(bool s = true);
    void set_skip_markers(bool s = true);
    void set_dynamic_markers(bool d = true);
    void set_force_skip_traits(bool s=true);
    void set_force_skip_markers(bool s = true);
    void set_force_dynamic_markers(bool d = true);
    void set_sex_linked_exist(bool s = true);
    void set_build_incremental(bool b = true);
    void set_sex_code_trait(bool a = false);
    void set_sex_field_name(const std::string & s);
    void set_pedigree_id_trait(bool a = false);
    void set_pedigree_id_name(const std::string & s);
    void set_format(const std::string & f);
    void set_treat_as_sibs(bool t);
    void set_no_sex_field(bool s);
    void set_no_sex_ok_option(bool s);

    cerrorstream error_sink() const;

    size_t verbose_output() const;
    
    bool reject_partial_lineage() const;
    bool require_record() const;
    bool skip_traits() const;
    bool skip_markers() const;
    bool dynamic_markers() const;
    bool force_skip_traits() const;
    bool force_skip_markers() const;
    bool force_dynamic_markers() const;
    bool sex_linked_exist() const;
    bool build_incremental() const;
    bool sex_code_trait() const;

    const std::string & sex_field_name() const;

    bool pedigree_id_trait() const;

    const std::string & pedigree_id_name() const;
    const std::string & format() const;

    bool get_treat_as_sibs() const;
    bool get_no_sex_field() const;
    bool get_no_sex_ok_option() const;

    virtual bool input(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout);
    virtual bool input_pedigree(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false) = 0;
    virtual bool input_data(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false) = 0;

    void validate();
    void invalidate();
    
    bool valid() const;

    virtual void print_mped(const RefMultiPedigree& p, 
                            const std::string&      filename, 
                            ostream&                messages = cout,
                            bool                    dump_trait = true, 
                            bool                    dump_marker = false);
    virtual bool output(RefMultiPedigree&  p, 
                        const std::string& filename, 
                        ostream&            messages = cout, 
                        bool               quiet = false) = 0;
          
                              
    size_t field_count() const;
    size_t trait_count() const;
    size_t invalid_trait_count() const;
    size_t string_count() const;
    size_t invalid_string_count() const;
    size_t marker_count() const;
    size_t invalid_marker_count() const;
    size_t invalid_marker_cov_count() const;
    size_t sex_count() const;
    size_t skip_count() const;
    size_t study_id_count() const;
    size_t pedigree_id_count() const;
    size_t individual_id_count()  const;
    size_t parent_id_count() const;

    const field_list_type & field_list() const;
          field_list_type & field_list();

    const vector<pair<std::string, std::string> > & get_ind_list() const;

  protected:

    void reset_counts();

    void print_family_structure_header(ostream &messages, const std::string &filename) const;
    void print_family_structure_footer(ostream &messages) const;
    void print_family_structure       (ostream &messages,
                                       const std::string& pn,
                                       const std::string& in,
                                       const std::string& s,
                                       const std::string& p1,
                                       const std::string& p2) const;

    void print_trait_header (ostream &messages, const RefMPedInfo &mped_info,  const std::string &filename) const;
    void print_trait_footer (ostream &messages) const;
    void print_trait        (ostream &messages,
                             const RefMPedInfo &mped_info,
                             const std::string& pn,
                             const std::string& in,
                             const vector<pair<size_t,std::string> >& trait_values,
                             const vector<pair<size_t,std::string> >& string_values) const;
  
    void print_marker_header(ostream &messages, const RefMPedInfo &mped_info, const std::string &filename) const;
    void print_marker_footer(ostream &messages) const;
    void print_marker       (ostream &messages,
                             const RefMPedInfo &mped_info,
                             const std::string& pn,
                             const std::string& in,
                             const vector<pair<size_t, std::string> >& marker_values) const;

    void parse_dynamic_markers_missing(RefMarkerInfo& d_marker, const std::string& missing);

    bool validate_fields(bool data_only = false, bool quiet = false);

    bool build_pedigree(RefMultiPedigree &p);
    bool build_data(RefMultiPedigree &p);

    void add_member(RefMultiPedigree &p, const std::string &pn,
                    const std::string &id, const std::string &sex,
                    const std::string &parent1, const std::string &parent2,
                    size_t line, size_t count=(size_t)-1);

    void update_marker_delimiter_info(RefMarkerInfo& marker, char sep);
    void update_marker_missing_info  (RefMarkerInfo& marker, const std::string& missing);

    // Update sex-linked marker information for members according to sex.
    //
    void update_sex_linked_marker_info(RefMultiPedigree &mp);

    /// Test for missing sexes of individuals in pedigree structures.  If there
    /// are individuals with unknown sex after teh build is complete (inferred 
    /// sexes are ok) that are linked into a pedigree structure, this function
    /// produces a warning message that indicates this to the user.
    ///
    /// NOTE:  This is a temporary item, created specifically to deal with ticket
    /// #1601 dealing with the fact that many of our algorithms require sex,
    /// even though the analysis results are not affected by it.  This warning 
    /// should be removed post 5.3 when sex checking should be done by the 
    /// specific algorithms, which can take actions appropriate to thier needs, 
    /// rather than globally where we don't have any ability to determine what 
    /// the needs are.
    bool do_no_sex_structural_test(const RefMultiPedigree& mp);

    bool report_pedigree_build_errors   (const Pedigree &p) const;
    bool report_pedigree_build_warnings (const Pedigree &p) const;

    // Data members
    //
    field_list_type my_fields;

    size_t     my_skip_count;
    size_t     my_trait_count;
    size_t     my_invalid_trait_count;
    size_t     my_string_count;
    size_t     my_invalid_string_count;
    size_t     my_marker_count;
    size_t     my_invalid_marker_count;
    size_t     my_invalid_marker_cov_count;
    size_t     my_sex_count;
    size_t     my_study_id_count;
    size_t     my_pedigree_id_count;
    size_t     my_individual_id_count;
    size_t     my_parent_id_count;
   
    mutable cerrorstream errors;

    std::vector<std::pair<std::string, std::string> > my_ind_list;

  private:

    bool         my_treat_as_sibs;
    size_t       my_verbose_output;
    bool         my_reject_partial_lineage;
    bool         my_require_record;
    bool         my_skip_traits;
    bool         my_skip_markers;
    bool         my_dynamic_markers;
    bool         my_force_skip_traits;
    bool         my_force_skip_markers;
    bool         my_force_dynamic_markers;
    bool         my_sex_linked_exist;
    bool         my_build_incremental;
    bool         my_sex_code_trait;
    bool         my_pedigree_id_trait;
    bool         my_no_sex_field;
    bool         my_no_sex_ok_option;
    bool         my_valid;
    std::string  my_sex_field_name;
    std::string  my_pedigree_id_name;
    std::string  my_format;
    
    // This set is only used for checking for duplicate records.
    std::multiset<std::pair<std::string, std::string> > my_inds;
};

class RefLSFPedigreeFile : virtual public RefPedigreeFile
{
  public:

    RefLSFPedigreeFile(cerrorstream &errors = sage_cerr);

    virtual bool process_parameters(RefMPedInfo &mped_info, const LSFBase *params);
    virtual bool process_parameter(RefMPedInfo &mped_info, const LSFBase *param);

    bool process_sex_code                 (RefMPedInfo &mped_info, const LSFBase* param);

  private:

    bool process_format                   (RefMPedInfo &mped_info, const LSFBase* param);
    bool process_verbose                  (RefMPedInfo &mped_info, const LSFBase* param);
    bool process_require_record           (RefMPedInfo &mped_info, const LSFBase* param);
    bool process_skip_traits              (RefMPedInfo &mped_info, const LSFBase* param);
    bool process_skip_markers             (RefMPedInfo &mped_info, const LSFBase* param);
    bool process_dynamic_markers          (RefMPedInfo &mped_info, const LSFBase* param);
    bool process_reject_partial_lineage   (RefMPedInfo &mped_info, const LSFBase* param);
    bool process_individual_missing_value (RefMPedInfo &mped_info, const LSFBase* param);
    bool process_no_sex_ok                (RefMPedInfo &mped_info, const LSFBase* param);
};

class RefDelimitedPedigreeFile : virtual public RefPedigreeFile
{
  public:
    
    RefDelimitedPedigreeFile(cerrorstream &err = sage_cerr);

    ~RefDelimitedPedigreeFile();

    virtual bool input(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout);
    virtual bool input_pedigree(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false);
    virtual bool input_data(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false);
    virtual bool output(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout, bool quiet = false);

    bool format_in_file() const;

    const std::string &delimiters() const;
    const std::string &whitespace() const;

    bool skip_consecutive_delimiters() const;
    bool skip_leading_delimiters()     const;
    bool skip_trailing_delimiters()    const;

    void set_format_in_file(bool f);
    void set_delimiters(const std::string &d);
    void set_whitespace(const std::string &w);
    void set_skip_consecutive_delimiters(bool skip=true);
    void set_skip_leading_delimiters(bool skip=true);
    void set_skip_trailing_delimiters(bool skip=true);
    void add_study_id_field(const std::string &field_name);
    void add_pedigree_id_field(const std::string &field_name);
    void add_individual_id_field(const std::string &field_name);
    void add_parent_id_field(const std::string &field_name);
    void add_sex_field(const std::string &field_name);
    void add_trait_field(const std::string &field_name, const std::string &trait_name);
    void add_string_field(const std::string &field_name, const std::string &string_name);
    void add_allele_field(const std::string &field_name, const std::string &marker_name);
    void add_marker_field(const std::string &field_name, const std::string &marker_name);

    const field_map_type & field_map() const;
          field_map_type & field_map();

  protected:

    bool validate_file  (const std::string &filename);
    bool validate_format(const std::string &filename);
    
    void setup_tokenizer(string_tokenizer& tokenizer);


  private:
    bool build_fields(string_tokenizer &header, const RefMPedInfo &mped_info,
                      bool quiet = false);
                      
    void build_field(const RefMPedInfo& mped_info, field_map_type::const_iterator f,
                     bool quiet);
                     
    void build_trait_field  (const RefMPedInfo& mped_info, field& current_field);
    void build_marker_field (const RefMPedInfo& mped_info, field& current_field);
    void build_string_field (const RefMPedInfo& mped_info, field& current_field);

    bool check_header_vs_field_map(string_tokenizer &header);

    virtual string get_marker_covariate_value(string m, const string& v1, const string& v2) { return ""; }

    bool   my_format_in_file;
    std::string my_delimiters;
    std::string my_whitespace;
    
    bool my_skip_leading_delimiters;
    bool my_skip_trailing_delimiters;
    bool my_skip_consecutive_delimiters;

    field_map_type  my_field_map;
};

class RefLSFDelimitedPedigreeFile : virtual public RefDelimitedPedigreeFile, virtual public RefLSFPedigreeFile
{
  public:

    RefLSFDelimitedPedigreeFile(cerrorstream &errors = sage_cerr);
    
    virtual bool process_parameters(RefMPedInfo &mped_info, const LSFBase *params);
    virtual bool process_parameter(RefMPedInfo &mped_info, const LSFBase *param);

    virtual bool input(RefMultiPedigree &p, const std::string &filename, ostream &messages = cout);

  private:

    struct MarkerListElement
    {
      std::string start_marker;
      std::string end_marker;
      
      LSF_ptr<LSFBase> marker_params;
    };

    std::list<MarkerListElement> my_marker_lists;

    std::list<MarkerListElement>::iterator find_marker_list(string start);

    bool build_marker_list_parameters (RefMPedInfo &mped_info);
    bool build_marker_list ( string_tokenizer& header,
                             string_tokenizer::iterator i, 
                             RefMPedInfo &mped_info,
                             std::list<MarkerListElement>::iterator mlist_index);

    std::list<MarkerListElement> my_covariate_lists;

    std::list<MarkerListElement>::iterator find_covariate_list(string start);

    bool build_covariate_list_parameters (RefMPedInfo &mped_info);
    bool build_covariate_list ( string_tokenizer& header,
                                string_tokenizer::iterator i, 
                                RefMPedInfo &mped_info,
                                std::list<MarkerListElement>::iterator mlist_index);

    bool process_format2                     (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_whitespace                  (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_delimiter_mode              (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_delimiters                  (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_skip_leading_delimiters     (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_skip_trailing_delimiters    (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_skip_consecutive_delimiters (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_study_id                    (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_pedigree_id                 (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_individual_id               (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_parent_id                   (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_sex_field                   (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_marker                      (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_phenotype                   (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_string                      (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_marker_list                 (RefMPedInfo &mped_info, const LSFBase *param);
    bool process_covariate_list              (RefMPedInfo &mped_info, const LSFBase *param);

    bool test_skip_marker(const LSFBase* param) const;
    bool test_allow_dynamic(const LSFBase* param) const;

    size_t setup_dynamic_marker(RefMPedInfo &mped_info, const LSFBase* param, const string& marker_name);  

    void parse_marker_parameters(RefMPedInfo &mped_info,const LSFBase* param, size_t marker_id);

};

//============================
//  INLINE FUNCTIONS
//
//  RefPedigreeFile
//============================

inline cerrorstream RefPedigreeFile::error_sink() const        { return errors; }
inline void RefPedigreeFile::set_error_sink(cerrorstream &err) { errors = err;  }

inline size_t RefPedigreeFile::verbose_output         () const { return my_verbose_output;                              }
inline bool   RefPedigreeFile::reject_partial_lineage () const { return my_reject_partial_lineage;                      }
inline bool   RefPedigreeFile::require_record         () const { return my_require_record;                              }
inline bool   RefPedigreeFile::skip_traits            () const { return my_force_skip_traits     || my_skip_traits;     }
inline bool   RefPedigreeFile::skip_markers           () const { return my_force_skip_markers    || my_skip_markers;    }
inline bool   RefPedigreeFile::dynamic_markers        () const { return my_force_dynamic_markers || my_dynamic_markers; }
inline bool   RefPedigreeFile::force_skip_traits      () const { return my_force_skip_traits;                           }
inline bool   RefPedigreeFile::force_skip_markers     () const { return my_force_skip_markers;                          }
inline bool   RefPedigreeFile::force_dynamic_markers  () const { return my_force_dynamic_markers;                       }
inline bool   RefPedigreeFile::sex_linked_exist       () const { return my_sex_linked_exist;                            }
inline bool   RefPedigreeFile::build_incremental      () const { return my_build_incremental;                           }
inline bool   RefPedigreeFile::sex_code_trait         () const { return my_sex_code_trait;                              }
inline bool   RefPedigreeFile::pedigree_id_trait      () const { return my_pedigree_id_trait;                              }

inline const std::string & RefPedigreeFile::sex_field_name()   const { return my_sex_field_name;   }
inline const std::string & RefPedigreeFile::pedigree_id_name() const { return my_pedigree_id_name; }
inline const std::string & RefPedigreeFile::format()           const { return my_format;           }

inline void RefPedigreeFile::set_treat_as_sibs          (bool   t) { my_treat_as_sibs          = t; }
inline void RefPedigreeFile::set_no_sex_field           (bool   s) { my_no_sex_field           = s; }
inline void RefPedigreeFile::set_verbose_output         (size_t v) { my_verbose_output         = v; }
inline void RefPedigreeFile::set_reject_partial_lineage (bool   r) { my_reject_partial_lineage = r; }
inline void RefPedigreeFile::set_require_record         (bool   a) { my_require_record         = a; }
inline void RefPedigreeFile::set_skip_traits            (bool   s) { my_skip_traits            = s; }
inline void RefPedigreeFile::set_skip_markers           (bool   s) { my_skip_markers           = s; }
inline void RefPedigreeFile::set_dynamic_markers        (bool   d) { my_dynamic_markers        = d; }
inline void RefPedigreeFile::set_force_skip_traits      (bool   s) { my_force_skip_traits      = s; }
inline void RefPedigreeFile::set_force_skip_markers     (bool   s) { my_force_skip_markers     = s; }
inline void RefPedigreeFile::set_force_dynamic_markers  (bool   d) { my_force_dynamic_markers  = d; }
inline void RefPedigreeFile::set_sex_linked_exist       (bool   s) { my_sex_linked_exist       = s; }
inline void RefPedigreeFile::set_build_incremental      (bool   b) { my_build_incremental      = b; }
inline void RefPedigreeFile::set_sex_code_trait         (bool   a) { my_sex_code_trait         = a; }
inline void RefPedigreeFile::set_pedigree_id_trait      (bool   a) { my_pedigree_id_trait      = a; }
inline void RefPedigreeFile::set_no_sex_ok_option       (bool   a) { my_no_sex_ok_option       = a; }

inline void RefPedigreeFile::set_sex_field_name   (const std::string &f) { my_sex_field_name = f;   }
inline void RefPedigreeFile::set_pedigree_id_name (const std::string &f) { my_pedigree_id_name = f; }
inline void RefPedigreeFile::set_format           (const std::string &f) { my_format = f;           }

inline void RefPedigreeFile::validate()    { my_valid = true;  }
inline void RefPedigreeFile::invalidate()  { my_valid = false; }
inline bool RefPedigreeFile::valid() const { return my_valid;  }

inline RefPedigreeFile::field::field(field_t t, const std::string &f, const std::string &n, size_t i) : 
  type(t), field_name(f), name(n), index(i) 
{ }

inline void RefPedigreeFile::field::invalidate()
{
  index = (size_t)-1;
  type  = skip;
}

inline size_t RefPedigreeFile::field_count          () const { return my_fields.size();        }
inline size_t RefPedigreeFile::trait_count          () const { return my_trait_count;          }
inline size_t RefPedigreeFile::invalid_trait_count  () const { return my_invalid_trait_count;  }
inline size_t RefPedigreeFile::string_count         () const { return my_string_count;         }
inline size_t RefPedigreeFile::invalid_string_count () const { return my_invalid_string_count; }
inline size_t RefPedigreeFile::marker_count         () const { return my_marker_count;         }
inline size_t RefPedigreeFile::invalid_marker_count () const { return my_invalid_marker_count; }
inline size_t RefPedigreeFile::sex_count            () const { return my_sex_count;            }
inline size_t RefPedigreeFile::skip_count           () const { return my_skip_count;           }
inline size_t RefPedigreeFile::study_id_count       () const { return my_study_id_count;       }
inline size_t RefPedigreeFile::pedigree_id_count    () const { return my_pedigree_id_count;    }
inline size_t RefPedigreeFile::individual_id_count  () const { return my_individual_id_count;  }
inline size_t RefPedigreeFile::parent_id_count      () const { return my_parent_id_count;      }

inline bool RefPedigreeFile::get_treat_as_sibs    () const { return my_treat_as_sibs;    }
inline bool RefPedigreeFile::get_no_sex_field     () const { return my_no_sex_field;     }
inline bool RefPedigreeFile::get_no_sex_ok_option () const { return my_no_sex_ok_option; }

inline const RefPedigreeFile::field_list_type & RefPedigreeFile::field_list() const { return my_fields; }
inline       RefPedigreeFile::field_list_type & RefPedigreeFile::field_list()       { return my_fields; }

inline const vector<pair<std::string, std::string> > & RefPedigreeFile::get_ind_list() const { return my_ind_list; }

//==========================
//  INLINE FUNCTIONS
//
//  RefDelimitedPedigreeFile
//==========================

inline bool RefDelimitedPedigreeFile::format_in_file() const     { return my_format_in_file; }

inline const std::string & RefDelimitedPedigreeFile::delimiters() const { return my_delimiters; }
inline const std::string & RefDelimitedPedigreeFile::whitespace() const { return my_whitespace; }

inline bool RefDelimitedPedigreeFile::skip_consecutive_delimiters () const { return my_skip_consecutive_delimiters; }
inline bool RefDelimitedPedigreeFile::skip_leading_delimiters     () const { return my_skip_leading_delimiters;     }
inline bool RefDelimitedPedigreeFile::skip_trailing_delimiters    () const { return my_skip_trailing_delimiters;    }

inline void RefDelimitedPedigreeFile::set_format_in_file              (bool f)    { my_format_in_file              = f;    }
inline void RefDelimitedPedigreeFile::set_skip_consecutive_delimiters (bool skip) { my_skip_consecutive_delimiters = skip; }
inline void RefDelimitedPedigreeFile::set_skip_leading_delimiters     (bool skip) { my_skip_leading_delimiters     = skip; }
inline void RefDelimitedPedigreeFile::set_skip_trailing_delimiters    (bool skip) { my_skip_trailing_delimiters    = skip; }

inline const RefDelimitedPedigreeFile::field_map_type & RefDelimitedPedigreeFile::field_map() const { return my_field_map; }
inline       RefDelimitedPedigreeFile::field_map_type & RefDelimitedPedigreeFile::field_map()       { return my_field_map; }

inline void RefDelimitedPedigreeFile::add_study_id_field      (const std::string &field_name)                                  { my_field_map[toUpper(field_name)] = field(study_id,      field_name);              }
inline void RefDelimitedPedigreeFile::add_pedigree_id_field   (const std::string &field_name)                                  { my_field_map[toUpper(field_name)] = field(pedigree_id,   field_name);              }
inline void RefDelimitedPedigreeFile::add_individual_id_field (const std::string &field_name)                                  { my_field_map[toUpper(field_name)] = field(individual_id, field_name);              }
inline void RefDelimitedPedigreeFile::add_parent_id_field     (const std::string &field_name)                                  { my_field_map[toUpper(field_name)] = field(parent_id,     field_name);              }
inline void RefDelimitedPedigreeFile::add_sex_field           (const std::string &field_name)                                  { my_field_map[toUpper(field_name)] = field(sex_code,      field_name);              }
inline void RefDelimitedPedigreeFile::add_trait_field         (const std::string &field_name, const std::string & trait_name)  { my_field_map[toUpper(field_name)] = field(trait,         field_name, trait_name);  }
inline void RefDelimitedPedigreeFile::add_string_field        (const std::string &field_name, const std::string & string_name) { my_field_map[toUpper(field_name)] = field(string_field,  field_name, string_name); }
inline void RefDelimitedPedigreeFile::add_allele_field        (const std::string &field_name, const std::string & marker_name) { my_field_map[toUpper(field_name)] = field(allele,        field_name, marker_name); }
inline void RefDelimitedPedigreeFile::add_marker_field        (const std::string &field_name, const std::string & marker_name) { my_field_map[toUpper(field_name)] = field(marker,        field_name, marker_name); }

} // End namespace RPED
} // End namespace SAGE

#endif
