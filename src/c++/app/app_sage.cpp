//==========================================================================
// File:      app_sage.cpp
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                              yes Nov 14
//
// Notes:     This source file implements S.A.G.E. app analyses.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "app/app_sage.h"

using namespace std;

namespace ONETOOL {

void
onetool::run_sage_analyses(app_data& ot_data, cIO* vin, SAGE::cerrormultistream& errors)
{
  string oname = OPT_STRING(out);

  const SAGE::RPED::RefMPedInfo& mped_info = ot_data.pedigrees().info();

  cout << endl << endl
       << "Performing PEDINFO analysis............................" << flush;

  SAGE::PEDINFO::MP_stats mps(errors);

  if( mped_info.trait_count() == 0 )
  {
    mps.compute(&ot_data.pedigrees());
  }
  else
  {
    mps.compute(&ot_data.pedigrees(), 0);
  }

  ofstream out_file;
  out_file.open(oname + ".pedinfo");
  print_title(out_file);

  SAGE::PEDINFO::mp_stats_viewer mp_view(out_file, mps);
  mp_view.view(ot_data.pedigrees().pedigree_count()>1);

  cout << "done." << endl << endl;

  if( mp_view.get_loop_count() )
  {
    errors << SAGE::priority(SAGE::error)
           << "Pedigree(s) with loop(s) detected. "
           << "FCOR, SEGREG and LODLINK analyses will be skipped!\n\n" << flush;

    return;
  }

  if( OPT_ENABLED(relpair) )
  {
    cout << endl
         << "Generating RELPAIR file................................" << flush;

    ofstream xml_file;
    xml_file.open(oname + ".relpair.xml");

    SAGE::FCOR::analysis_option_type fcor_opt;

    SAGE::FCOR::FcorAnalysis fcor_anal(fcor_opt, &ot_data.pedigrees(), errors);
    fcor_anal.view_pairs_xml(xml_file);

    cout << "done." << endl << endl;
  }

  if( OPT_ENABLED(fcor) )
  {
    if( mped_info.trait_count() == 0 )
    {
      errors << SAGE::priority(SAGE::error)
             << "No phenotype names are selected for further analyses. "
             << "FCOR analysis will be skipped!\n\n" << flush;
    }
    else
    {
      cout << endl
           << "Performing FCOR analysis..............................." << flush;

      ostream* null_o = NULL;
      
      ofstream out_file;
      out_file.open(oname + ".fcor");
      print_title(out_file);

      SAGE::FCOR::analysis_option_type fcor_opt;
      parse_fcor_options(fcor_opt);

      SAGE::FCOR::FcorAnalysis fcor_anal(fcor_opt, &ot_data.pedigrees(), errors);
      fcor_anal.run_analysis(out_file, *null_o, *null_o, *null_o, *null_o);

      cout << "done." << endl << endl;
    }
  }

  if( OPT_ENABLED(segreg) )
  {
    if( mped_info.trait_count() == 0 )
    {
      errors << SAGE::priority(SAGE::error)
             << "No phenotype names are selected for further analyses. "
             << "SEGREG analysis will be skipped!\n\n" << flush;
    }
    else
    {
      cout << endl
           << "Performing SEGREG analysis............................." << endl << flush;

      SAGE::SEGREG::AnalysisController segreg_anal(ot_data.get_ostreams());

      if( IS_ASSIGNED(par) )
      {
        segreg_anal.do_analysis_from_parameter_file(ot_data.pedigrees(), OPT_STRING(par), oname);
      }
      else
      {
        SAGE::SEGREG::model s_model;
        parse_segreg_options(s_model, ot_data.pedigrees(), errors);

        segreg_anal.do_analysis(ot_data.pedigrees(), s_model);
      }
    }
  }

  if( OPT_ENABLED(lodlink) )
  {
    if( IS_ASSIGNED(typ) )
    {
      cout << endl
           << "Performing LODLINK analysis............................" << endl << flush;

      ofstream  sum_file;
      sum_file.open(oname + ".lodlink.sum");
      print_title(sum_file);

      ofstream  det_file;
      det_file.open(oname + ".lodlink.det");
      print_title(det_file);

      SAGE::LODLINK::instructions lodlink_inst;
      parse_lodlink_options(lodlink_inst, ot_data.pedigrees(), errors);

      SAGE::LODLINK::analysis ld_anal(lodlink_inst, ot_data.pedigrees(), errors);
      ld_anal.run_analysis(sum_file, det_file);
    }
    else
    {
      errors << SAGE::priority(SAGE::error)
             << "The type probability file name is not specified. "
             << "LODLINK analysis will be skipped!\n\n" << flush;
    }
  }

  return;
}

void
parse_fcor_options(SAGE::FCOR::analysis_option_type& fcor_opt)
{
  if( OPT_ENABLED(fcorStdErrOff) )
    fcor_opt.standard_error = false;

  return;
}

void
parse_lodlink_options(SAGE::LODLINK::instructions& lodlink_inst,
                      const SAGE::RPED::RefMultiPedigree& mp, SAGE::cerrorstream& errors)
{
  lodlink_inst.title          = "LODLINK Analysis";
  lodlink_inst.file_name_root = "lodlink_analysis";
  lodlink_inst.linkage = SAGE::LODLINK::instructions::MARKER;
  lodlink_inst.trait = mp.info().marker_info(mp.info().marker_count()-1).name();

  if( OPT_ENABLED(lodlinkLinkageTestOff) )
    lodlink_inst.linkage_test = false;

  if( OPT_ENABLED(lodlinkLinkageHomogOff) )
    lodlink_inst.linkage_homog = false;

  if( OPT_ENABLED(lodlinkLinkageSexSpecific) )
    lodlink_inst.linkage_sex_specific = true;

  if( OPT_ENABLED(lodlinkSmithHomogTest) )
      lodlink_inst.smiths_test = true;

  if( OPT_ENABLED(lodlinkGenotypes) )
    lodlink_inst.genotypes = true;

  return;
}

void
parse_segreg_options(SAGE::SEGREG::model& s_model,
                     const SAGE::RPED::RefMultiPedigree& mp, SAGE::cerrorstream& errors)
{
  string oname = OPT_STRING(out);
  s_model.reset(errors);
  s_model.set_file_name_root(oname + ".segreg");

  // set primary_trait
  s_model.set_primary_trait(mp.info().trait_info(0).name());

  if( mp.info().trait_info(0).type() == SAGE::RPED::RefTraitInfo::continuous_trait )
  {
    s_model.set_primary_trait_type(SAGE::SEGREG::pt_CONTINUOUS);
  }
  else if( mp.info().trait_info(0).type() == SAGE::RPED::RefTraitInfo::binary_trait )
  {
    s_model.set_primary_trait_type(SAGE::SEGREG::pt_BINARY);
  }
  else
  {
    s_model.set_model_class(SAGE::SEGREG::model_INVALID);

    errors << SAGE::priority(SAGE::error)
           << "Invalid primary trait type detected. Segregation analysis will be skipped!"
           << endl;
  }

  // set output
  s_model.set_type_prob(true);

  // TYPE MEAN
  // TYPE SUSCEPT
  // COMPOSITE TRAIT

  // set covariate - MEAN COVARIATE, VARIANCE COVARIATE, SUSCEPTIBILITY COVARIATE
  SAGE::MAXFUN::model_input coeff = SAGE::MAXFUN::model_input(SAGE::QNAN, SAGE::SEGREG::COV_DEFAULT_FIXED);
  bool                interaction = SAGE::SEGREG::COV_DEFAULT_INTERACTION;

  for( size_t t = 0; t < mp.info().trait_count(); ++t )
  {
    if( mp.info().trait_info(t).usage() == SAGE::RPED::RefTraitInfo::trait_covariate )
    {
      bool set_valid = s_model.susc_cov_sub_model.add_covariate(&mp, mp.info().trait_info(t).name(),
                                                                coeff.value, interaction, coeff.fixed,
                                                                s_model.get_type_missing());
      if( !set_valid )
        s_model.set_model_class(SAGE::SEGREG::model_INVALID);
    }
  }

  // set class
  if( mp.info().trait_info(0).type() == SAGE::RPED::RefTraitInfo::binary_trait )
    s_model.set_model_class(SAGE::SEGREG::model_MLM);

  if( s_model.get_model_class() != SAGE::SEGREG::model_INVALID )
    s_model.resid_sub_model.set_as_default(s_model.get_model_class());

  // FPMM
  // RESID

  // set transformation
  if( mp.info().trait_info(0).type() == SAGE::RPED::RefTraitInfo::binary_trait )
    s_model.transf_sub_model.set_option(SAGE::MAXFUN::TransformationSubmodel::no_trans);

  // GENO FREQUENCY
  // TRANSMISSION
  // set variance
  // set ascertainment
  // set prevalence
  s_model.prev_sub_model.set_fpmm_option(s_model.get_model_class() == SAGE::SEGREG::model_FPMM);
  s_model.prev_sub_model.set_onset_option(s_model.get_primary_trait_type() == SAGE::SEGREG::pt_ONSET);

  // set prev_constraint
  // set prev_estimate

  // check_meta_constraints
  /*
  bitset<SAGE::SEGREG::i_num_i> incon = s_model.check_consistency();
  for( int i = 0; i < SAGE::SEGREG::i_num_i; ++i )
  {
    if( incon.test(i) )
    {
      errors << SAGE::priority(SAGE::error)
             << "incon detected!"
             << endl;
    }
  }
  */
  return;
}

} // end of namespace ONETOOL
