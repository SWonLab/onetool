#ifndef SEGREG_ANALYSIS_H
#include "sage/segreg/analysis_iterative_models.h"
#endif

namespace SAGE   {
namespace SEGREG {

//=========================
// primary_analysis inlines
//=========================

inline 
primary_analysis::primary_analysis
  (Output_Streams& o)
  : out(o),
    errors(out.errors()),
    my_quality(true),
    my_models(), 
    my_results()
{ }

inline
primary_analysis::~primary_analysis()
{ }

inline
const primary_analysis_results& primary_analysis::get_results() const
{
  return my_results;
}

inline
bool primary_analysis::is_quality() const
{
  return my_quality;
}

inline
void primary_analysis::set_quality(bool b)
{
  my_quality = b;
}

//=================================
// primary_analysis_results inlines
//=================================

inline 
primary_analysis_results::primary_analysis_results()
  : my_valid(false),
    my_ped_data(), 
    my_model(), 
    my_maxfun_results(),
    my_subpedigree_count(0),
    my_unconnected_count(0),
    my_type_probs(),
    my_pfunc_probs(),
    my_mlm_resid_corr()
{ }

//lint -e{1554}
inline 
primary_analysis_results::primary_analysis_results
  (const primary_analysis_results& p)
  : my_valid            (p.my_valid),
    my_ped_data         (p.my_ped_data), 
    my_model            (p.my_model), 
    my_maxfun_results   (p.my_maxfun_results),
    my_subpedigree_count(p.my_subpedigree_count),
    my_unconnected_count(p.my_unconnected_count),
    my_type_probs       (p.my_type_probs),
    my_pfunc_probs      (p.my_pfunc_probs),
    my_mlm_resid_corr   (p.my_mlm_resid_corr)
{ }

inline 
primary_analysis_results::~primary_analysis_results()
{ }

inline 
primary_analysis_results&
    primary_analysis_results::operator=(const primary_analysis_results& p)
{
  my_valid             = p.my_valid;

  //lint -e{1555}
  my_ped_data          = p.my_ped_data; 
  my_model             = p.my_model; 
  my_maxfun_results    = p.my_maxfun_results;
  my_subpedigree_count = p.my_subpedigree_count;
  my_unconnected_count = p.my_unconnected_count;
  my_type_probs        = p.my_type_probs;
  my_pfunc_probs       = p.my_pfunc_probs;
  my_mlm_resid_corr    = p.my_mlm_resid_corr;

  return *this;
}

inline
bool primary_analysis_results::is_valid() const
{
  return my_valid;
}

inline
const model& primary_analysis_results::get_final_model() const
{
  return my_model;
}

inline
const MAXFUN::Results& primary_analysis_results::get_maxfun_results() const
{
  return my_maxfun_results;
}

inline
uint primary_analysis_results::get_subpedigree_count() const
{
  return my_subpedigree_count;
}

inline
uint primary_analysis_results::get_unconnected_count() const
{
  return my_unconnected_count;
}

inline
const pg::post_geno_map& primary_analysis_results::get_type_probs() const
{
  return my_type_probs;
}

inline
const pf::pen_func_map& primary_analysis_results::get_pfunc_probs() const
{
  return my_pfunc_probs;
}

inline
const MlmResidCorrelationCalculator& primary_analysis_results::get_mlm_resid_corr() const
{
  return my_mlm_resid_corr;
}

//=================================
// Iterative_Models inlines
//=================================

inline
Iterative_Models::Iterative_Models
    (const PedigreeDataSet& ped_data,
     Output_Streams&   o,
     primary_analysis& a)
  : my_ped_data(ped_data),
    analyzer(a),
    out(o),
    output(true),
    output_mean(true),
    output_transm(true),
    one_mean_res(QNAN),
    one_susc_res(QNAN) 
{ }

inline
void Iterative_Models::set_target_model(const model& target)
{
  clear();

  my_target_model = target;

  // Determine the completeness of the three primary submodels

  target_mean_complete   = target.type_dependent_sub_model().is_complete();
  target_transm_complete = target.transm_sub_model          .is_complete();
  target_freq_complete   = target.freq_sub_model            .is_complete();

  // Initialize Sampler

  // The sample is initialized with the type dependent trait.  In general,
  // this is the primary trait, but in the case of age of onset, it may not
  // be, since the 'primary' trait is the binary component, which is not
  // always the type dependent component.  So if it's age of onset, and not
  // suceptibility dependent on type, we have to do a bit more work.
  
  //lint -e(534) since return value not needed

  if(my_target_model.get_primary_trait_type() != pt_ONSET ||
     my_target_model.ons_sub_model.t_option() == onset_sub_model::t_S)
    my_trait_sample(*my_ped_data.get_raw_data(), my_target_model.get_primary_trait());
  else
  {
    // We have to do this based on either the age_onset or age_exam (if
    // there aren't any age_onsets)
    
    my_trait_sample(*my_ped_data.get_raw_data(), my_target_model.ons_sub_model.age_of_onset());

    if(my_trait_sample.get_n() == 0)
      my_trait_sample(*my_ped_data.get_raw_data(), my_target_model.ons_sub_model.age_at_exam());
  }
}

inline
const model& Iterative_Models::get_target_model() const
{
  return my_target_model;
}

inline void Iterative_Models::set_output(bool o)
{
  output = o;
}

inline bool Iterative_Models::get_output() const
{
  return output;
}

inline void Iterative_Models::set_output_mean(bool o)
{
  output_mean = o;
}

inline bool Iterative_Models::get_output_mean() const
{
  return output_mean;
}

inline void Iterative_Models::set_output_transm(bool o)
{
  output_transm = o;
}

inline bool Iterative_Models::get_output_transm() const
{
  return output_transm;
}

// =========================================
// Iterative_Models::model_results functions
// =========================================

inline
Iterative_Models::model_results::model_results()
{
  clear();
}
          
inline void
Iterative_Models::model_results::clear()
{
  availability = empty;
  result       = primary_analysis_results();
}
              
inline bool
Iterative_Models::model_results::is_good() const
{
  return availability == good;
}

inline bool
Iterative_Models::model_results::is_bad() const
{
  return availability == bad;
}

inline bool
Iterative_Models::model_results::is_empty() const
{
  return availability == empty;
}

inline bool
Iterative_Models::model_results::is_available() const
{
  return availability != empty;
}

// ================================================
// Iterative_Models::model_classification functions
// ================================================

inline
Iterative_Models::model_classification::model_classification()
  : mean_opt  (genotype_specific_mean_sub_model::one),
    transm_opt(transmission_sub_model::homog_no_trans)
{ }

inline
Iterative_Models::model_classification::model_classification
    (mean_option m, transm_option t)
  : mean_opt  (m),
    transm_opt(t)
{ }

inline bool
Iterative_Models::model_classification::operator< (const model_classification& rhs) const
{
  return mean_opt < rhs.mean_opt ||
         (mean_opt == rhs.mean_opt && transm_opt < rhs.transm_opt);
}

// =====================================
// IM protected/private member functions
// =====================================

inline
const model& Iterative_Models::get_initial_model()
{
  if(!initial_model_built) create_initial_model(my_initial_model);

  return my_initial_model;
}

inline 
bool Iterative_Models::use_target_mean_model(mean_option mean_opt) const
{
  return target_mean_complete && my_target_model.type_dependent_sub_model().option() == mean_opt;
}

inline
bool Iterative_Models::use_target_freq_model(mean_option mean_opt) const
{
  if(target_freq_complete)
  {
    if(mean_opt ==  genotype_specific_sub_model::two     ||
       mean_opt ==  genotype_specific_sub_model::two_rec ||
       mean_opt ==  genotype_specific_sub_model::two_dom)
      return my_target_model.freq_sub_model.option() == genotype_frequency_sub_model::hwe;
  }

  return target_freq_complete;
}

inline
bool Iterative_Models::use_target_transm_model(transm_option transm_opt) const
{
  return target_transm_complete && my_target_model.transm_sub_model.option() == transm_opt;
}

inline
void Iterative_Models::update_means_to_target(model_vector& v) const
{
    for(size_t i = 0; i < v.size(); ++i)
      v[i].type_dependent_sub_model() = my_target_model.type_dependent_sub_model();
}

inline
void Iterative_Models::update_freqs_to_target(model_vector& v) const 
{
    for(size_t i = 0; i < v.size(); ++i)
      v[i].freq_sub_model = my_target_model.freq_sub_model;
}

inline
void Iterative_Models::set_one_mean_res(double meanval)
{
     one_mean_res = meanval;
}

inline
void Iterative_Models::set_one_susc_res(double suscval)
{
     one_susc_res = suscval;
}


/// Static functions for use by the iterative models.
///
/// These functions are to make the choice between discrete and continuous
/// models a little easier.
///
/// NOTE: They currently only distinguish between discrete (binary) and
///       continuous models.  They do not, at present, know anything
///       about age of onset variables which are both.

//@{

inline bool is_model_continuous(const model& m)
{
  return m.get_primary_trait_type() == pt_CONTINUOUS ||
         (m.get_primary_trait_type() == pt_ONSET &&
          m.ons_sub_model.t_option() == onset_sub_model::t_A);
}

inline bool is_model_binary(const model& m)
{
  return m.get_primary_trait_type() == pt_BINARY ||
         (m.get_primary_trait_type() == pt_ONSET &&
          m.ons_sub_model.t_option() == onset_sub_model::t_S);
}

inline bool is_model_onset(const model& m)
{
  return m.get_primary_trait_type() == pt_ONSET;
}

/// get_type_model() returns the mean or susceptibility sub_model as
/// necessary based on model type

inline
const genotype_specific_mean_susc_sub_model& get_type_model(const model& m)
{
  return m.type_dependent_sub_model();
}

/// get_type_model() returns the mean or susceptibility sub_model as
/// necessary based on model type

inline
genotype_specific_mean_susc_sub_model& get_type_model(model& m)
{
  return m.type_dependent_sub_model();
}

inline
string singular(const model& m)
{
  const string mean_singular("mean");
  const string susc_singular("susc.");

  if(is_model_continuous(m)) return mean_singular;
  else                       return susc_singular;
}

inline
string plural(const model& m)
{
  const string mean_plural("means");
  const string susc_plural("suscs.");

  if(is_model_continuous(m)) return mean_plural;
  else                       return susc_plural;
}

//@}


}
}
