//==========================================================================
// File:      app_data.cpp
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                            yes Jan 2015
//
// Notes:     This header file defines the base structure for input data.
//
// Copyright (c) 2015 Sungho Won
//   All Rights Reserved
//==========================================================================
#include "input/app_data.h"

namespace ONETOOL {

void
print_fill(char* fname)
{
  // Print out a number of .'s equal to 35 - the filename size to line up the
  // columns
  char old_fill = cout.fill('.');

  if( string(fname).size() < 35 )
  {
    cout << std::setw(35-string(fname).size()) << '.';
  }

  cout << "done." << endl << endl;
  cout.fill(old_fill);

  return;
}

pair<string, string>
get_acgt_allele(char al, char a2, char N_geno)
{
  // Convert geno num to a1/a2
  char al1 = N_geno>1;
  char al2 = N_geno>0;

  string ia1;
  string ia2;

  if( (wsReal)al1 == 0 )
    ia1 = string(1, al);
  else
    ia1 = string(1, a2);

  if( (wsReal)al2 == 0 )
    ia2 = string(1, al);
  else
    ia2 = string(1, a2);

  return make_pair(ia1, ia2);
}

app_data::app_data(string name, bool debug)
        : my_program_name(name), my_output(name, debug)
{}

app_data::~app_data()
{}

cIO*
app_data::input_wisard()
{
  cOption& user_opt = OPTION();

  if (IS_ASSIGNED(vcf)) {
	  cIO* vin = cIO::summon(FT_VCF, OPT_STRING(vcf));
	  vin->_finalize();

	  return vin;
  } else if (IS_ASSIGNED(bed)) {
	  cIO* vin = cIO::summon(FT_BED, OPT_STRING(bed));
	  vin->_finalize();

	  return vin;
  } else if (IS_ASSIGNED(ped)) {
	  cIO* vin = cIO::summon(FT_PED, OPT_STRING(ped));
	  vin->_finalize();

	  return vin;
  } else if (IS_ASSIGNED(dosage)) {
	  cIO* vin = cIO::summon(FT_DOSAGE, OPT_STRING(dosage));
	  vin->_finalize();

	  return vin;
  } else if (IS_ASSIGNED(fam)) {
	  user_opt.N_szvar = 0;
	  user_opt.assign("szvar", "0");
	  user_opt.FORCE_OPT_NUMBER(szvar);
	  cIO* vin = cIO::summon(FT_SIM, OPT_STRING(fam));
	  vin->_finalize();

	  return vin;
  }

  LOGwarn("An input is required to run ONETOOL");
  return NULL;
}

void
app_data::input_sage(cIO* vin)
{
  cOption& user_opt = OPTION();

  read_pheno_file(vin, OPT_STRING(pheno));

  read_locus_description_file(vin, OPT_STRING(typ));

  read_pedigree_file(vin, OPT_STRING(fam));

  return;
}

bool
app_data::read_pedigree_file(cIO* vin, char* fname)
{
  //cout << "\nReading Pedigree File.................................." << endl;
  //cout << "               from " << fname << flush;

  if( IS_ASSIGNED(misparent) )
    my_pedigrees.info().set_individual_missing_code(OPT_STRING(misparent));
  else
    my_pedigrees.info().set_individual_missing_code("0");

  bool pedigree_loaded =  false;

  // 1. Set pedigree members
  //
  vSampPtr& Xv_samp = vin->getSample();
  mSamp& Xm_samp = vin->getSampleData();

  for( mSamp_it i = Xm_samp.begin() ; i != Xm_samp.end() ; i++ )
  {
    xSample& thisSample = i->second;
    //cout << "add mem " << thisSample.S_IID << endl;

    if( thisSample.N_sex == 1 )
      my_pedigrees.add_member(thisSample.S_FID, thisSample.S_IID, SAGE::MPED::SEX_MALE);
    else if( thisSample.N_sex == 2 )
      my_pedigrees.add_member(thisSample.S_FID, thisSample.S_IID, SAGE::MPED::SEX_FEMALE);
    else
      my_pedigrees.add_member(thisSample.S_FID, thisSample.S_IID, SAGE::MPED::SEX_MISSING);

    if( thisSample.Xp_mat != NULL) {
      if( thisSample.Xp_pat != NULL )
        my_pedigrees.add_lineage(thisSample.S_FID, thisSample.S_IID, thisSample.Xp_mat->S_IID, thisSample.Xp_pat->S_IID);
      else
        my_pedigrees.add_lineage(thisSample.S_FID, thisSample.S_IID, thisSample.Xp_mat->S_IID);
    }
    else if( thisSample.Xp_pat != NULL )
      my_pedigrees.add_lineage(thisSample.S_FID, thisSample.Xp_mat->S_IID, thisSample.Xp_pat->S_IID);
  }

  // 2. Build pedigrees & infos
  //
  my_pedigrees.build();

  SAGE::RPED::RefMPedInfo& mped_info = my_pedigrees.info();

  SAGE::RPED::RefMultiPedigree::pedigree_iterator j;
  for( j = my_pedigrees.pedigree_begin(); j != my_pedigrees.pedigree_end(); ++j )
  {
    j->info().build(*j);
    j->info().resize_traits(mped_info.trait_count());
    j->info().resize_markers(mped_info.marker_count(), mped_info);
  }

  //print_fill(fname);

  //cout << "\nSorting Pedigrees......................................" << flush;
  for( j = my_pedigrees.pedigree_begin(); j != my_pedigrees.pedigree_end(); ++j )
    PedigreeSort(*j);
  //cout << "done." << endl << endl;

#if 0
  cout << "trait cnt  = " << mped_info.trait_count() << endl;
  cout << "marker cnt = " << mped_info.marker_count() << endl;
#endif

  // 3. Set member's pheno, covariate, marker values
  vPheno&   phe_info = vin->getPhenoInfo();
  vCovar&   cov_info = vin->getCovInfo();
  vVariant& var_info = vin->getVariant();

  wsMat  phe_vals = vin->getPhenos();
  wsMat  cov_vals = vin->getCovariates();
  char** var_vals = vin->getGenotype();

  int N_idx = 0;
  for( vSampPtr_it i = Xv_samp.begin() ; i != Xv_samp.end() ; i++,N_idx++ )
  {
    xSample& thisSample = **i;

    if( thisSample.N_idx != N_idx) halt("Sample data indexing error");

    SAGE::RPED::RefMultiPedigree::member_pointer mem = my_pedigrees.member_find(thisSample.S_FID, thisSample.S_IID);

    if( !mem ) continue;

    SAGE::RPED::RefMultiPedigree::pedinfo_type& pinfo = mem->pedigree()->info();
    int ind_num = mem->index();
#if 0
    cout << "set trait/marker for mem " << thisSample.S_IID
         << ", " << thisSample.N_oriIdx << ", " << thisSample.N_idx << ", " << N_idx
         << "\t" << ind_num << ", " << mem->name() << endl;
#endif
//    if( IS_ASSIGNED(pname) )
//    {
      for( size_t p = 0; p < phe_info.size(); ++p )
      {
        double p_val = phe_vals[p][N_idx];
        size_t t     = mped_info.trait_find(phe_info[p].S_name);

        if( p_val == mped_info.trait_info(t).numeric_missing_code() )
          p_val = SAGE::QNAN;

        int code = pinfo.set_trait(ind_num, t, SAGE::doub2str(p_val), mped_info.trait_info(t));
#if 0
        cout << "  p = " << p << ", p_name = " << phe_info[p].S_name
             << "\tt = " << t
             << ", val = " << p_val << "\tcode = " << code
             << "\t" << phe_vals[p][N_idx]
             << "\t" << phe_vals[p][ind_num] << endl;
#endif
      }
//    }

    if( IS_ASSIGNED(cname) )
    {
      for( size_t c = 0; c < cov_info.size(); ++c )
      {
        double c_val = cov_vals[c][N_idx];
        //double c_val = cov_vals[c][thisSample.N_oriIdx];
        size_t t     = mped_info.trait_find(cov_info[c].Sp_varName);

        if( c_val == mped_info.trait_info(t).numeric_missing_code() )
          c_val = SAGE::QNAN;

        int code = pinfo.set_trait(ind_num, t, SAGE::doub2str(c_val), mped_info.trait_info(t));
#if 0
        cout << "  c = " << c << ", c_name = " << cov_info[c].Sp_varName
             << ", t = " << t
             << ", c_val = " << c_val << ", code = " << code << endl;
#endif
      }
    }

    if( OPT_ENABLED(lodlink) )
    {
      for( size_t v = 0; v < var_info.size(); ++v )
      {
        string v_name = var_info[v].name;
        char N_geno = var_vals[N_idx][v];

        if( !isMissing(N_geno) )
        {
          pair<string, string> als = get_acgt_allele(var_info[v].al1, var_info[v].al2, N_geno);

          char al1 = N_geno>1;
          char al2 = N_geno>0;
          //stringstream s1; s1 << (wsReal)al1;
          //stringstream s2; s2 << (wsReal)al2;

          size_t m = mped_info.marker_find(v_name);

          int code = pinfo.set_phenotype(ind_num, mem->get_effective_sex(), m,
                                         //s1.str(), s2.str(), mped_info.marker_info(m));
                                         als.first, als.second, mped_info.marker_info(m));
#if 0
          cout << "  v=" << v << ", v_name=" << v_name
               << ", a1=" << var_info[v].al1 << ", a2=" << var_info[v].al2
               << ", m=" << m << ", N_geno='" << (wsReal)N_geno
               << "'(" << (wsReal)al1 << "/" << (wsReal)al2 << ")"
               //<< " (" << s1.str() << "/" << s2.str() << ")"
               << " (" << als.first << "/" << als.second << ")"
               << ", code = " << code << endl;
#endif
        }
      }
    }
  }

  SAGE::RPED::print_mped(info(), my_pedigrees);

  // 4. check marker data
  if( IS_ASSIGNED(typ) )
  {
    for( size_t m = 0; m < mped_info.marker_count(); ++m )
    {
      MLOCUS::inheritance_model& model = mped_info.marker_info(m);

      if( model.allele_count() == 0 )
      {
        my_output.errors() << SAGE::priority(SAGE::warning)
                           << "No usable data found for marker '"
                           << model.name()
                           << "'. This marker will be considered unknown for all analyses."
                           << endl;
        model.add_allele("A", 1.0, true, true);
      }
      else if( toUpper(model.name()) == toUpper(mped_info.markers().name(m)) + "~~~" )
      {
        for( RPED::RefMultiPedigree::pedigree_iterator ped = my_pedigrees.pedigree_begin(); ped != my_pedigrees.pedigree_end(); ++ped )
        {
          for( RPED::RefMultiPedigree::member_const_iterator ind = ped->member_begin(); ind != ped->member_end(); ++ind )
          {
            if( ped->info().phenotype_missing(ind->index(), m, model) )
            {
              ped->info().set_phenotype((size_t) ind->index(), m, ped->name() + "~" + ind->name(), std::string(), model);
            }
          }
        }
        model.set_name(mped_info.markers().name(m));
      }
    }
  }

  return pedigree_loaded;
}

bool
app_data::read_vcf_file(cIO* vin, char* fname)
{
  //cout << "\nReading Variant Calling Format (VCF) File....." << flush;
  //vin.read_vcf_file(fname);
  //cout << "done." << endl;

  return true;
}

bool
app_data::read_pheno_file(cIO* vin, char* fname)
{
  //cout << "\nReading Phenotype File................................." << flush;

  // Already done in vin->_finalize()
  //loadSampVar(vin);

  // Add in pheno & covariate meta data
  SAGE::RPED::RefMPedInfo& mped_info = my_pedigrees.info();

  string p_miss = IS_ASSIGNED(mispheno) ? OPT_STRING(mispheno) : ".";

//  if( IS_ASSIGNED(pname) )
//  {
    vPheno& phe_info = vin->getPhenoInfo();
    for( size_t p = 0; p < phe_info.size(); ++p )
    {
      string p_name = phe_info[p].S_name;
//LOG("Phenotype [%s] found\n", p_name.c_str());
      if( vin->isContinuous(p) )
      {
        size_t t = mped_info.add_continuous_trait(p_name, SAGE::RPED::RefTraitInfo::trait_variate);

        mped_info.trait_info(t).set_string_missing_code(p_miss);
        mped_info.trait_info(t).set_numeric_missing_code(SAGE::str2doub(p_miss));
      }
      else
      {
        size_t t = mped_info.add_binary_trait(p_name, SAGE::RPED::RefTraitInfo::trait_variate);

        mped_info.trait_info(t).set_string_missing_code(p_miss);
        mped_info.trait_info(t).set_numeric_missing_code(SAGE::str2doub(p_miss));
        mped_info.trait_info(t).set_string_affected_code("1");
        mped_info.trait_info(t).set_string_unaffected_code("0");
        mped_info.trait_info(t).set_numeric_affected_code(1);
        mped_info.trait_info(t).set_numeric_unaffected_code(0);
        mped_info.trait_info(t).set_threshold(SAGE::QNAN);
      }
    }
//  }

  if( IS_ASSIGNED(cname) )
  {
    vCovar& cov_info = vin->getCovInfo();
    for( size_t c = 0; c < cov_info.size(); ++c )
    {
      string c_name = cov_info[c].Sp_varName;

      size_t t = mped_info.add_continuous_trait(c_name, SAGE::RPED::RefTraitInfo::trait_covariate);

      mped_info.trait_info(t).set_string_missing_code(p_miss);
      mped_info.trait_info(t).set_numeric_missing_code(SAGE::str2doub(p_miss));
    }
  }
#if 0
  cout << "phe cnt = " << vin->getPhenoInfo().size() << endl;
  cout << "cov cnt = " << vin->getCovInfo().size() << endl;
  cout << "trait cnt = " << mped_info.trait_count() << endl;
#endif
  //cout << "done." << endl;

  return true;
}

bool
app_data::read_locus_description_file(cIO* vin, char* fname)
{
  SAGE::RPED::RefMPedInfo& mped_info = my_pedigrees.info();

  char separator = '/';
  string a_miss  = ".";

  SAGE::MLOCUS::inheritance_model_map& marker_map = mped_info.markers();

  if( OPT_ENABLED(lodlink) )
  {
    // Add in SNP meta data from vcf file.
    vVariant& var_info = vin->getVariant();
    xMaf* maf_info = vin->getMAF();

    for( size_t v = 0; v < var_info.size(); ++v )
    {
      string v_name = var_info[v].name;
      if( marker_map.name_find(v_name) != marker_map.name_end() )
      {
        my_output.errors() << SAGE::priority(SAGE::warning)
                           << "Variant '" << v_name << "' is already present.  It will be overwritten."
                           << endl;
      }
#if 0
  cout << v << " " << v_name
       << ", maf = " << maf_info[v].R_maf
       << ", allmaf = " << maf_info[v].R_allMaf
       << ", N_allele = " << maf_info[v].N_allele
       << ", N_allMac = " << maf_info[v].N_allMac
       << ", N_mac = " << maf_info[v].N_mac
       << ", B_flipped = " << maf_info[v].B_flipped
       << endl;
#endif

      double d = maf_info[v].R_allMaf;
      stringstream s1; s1 << var_info[v].al1;
      stringstream s2; s2 << var_info[v].al2;
#if 0
  cout << v << " " << v_name
       << ", allele1 = " << s1.str() << ", d1 = " << 1.-d
       << ", allele2 = " << s2.str() << ", d2 = " << d
       << endl;
#endif
      SAGE::MLOCUS::genotype_model gm(v_name);
      gm.set_unphased_separator(separator);
      gm.set_missing_allele_name(a_miss);
      gm.add_allele(s1.str(), 1.-d);
      gm.add_allele(s2.str(), d);

      SAGE::MLOCUS::penetrance_model& pm = marker_map[v_name] = SAGE::MLOCUS::penetrance_model(gm, true, true);
      pm.set_name(v_name);
      pm.set_missing_phenotype_name(a_miss);

      string sep = gm.separators();
      pm.alias_phenotype(pm.get_missing_phenotype_id(), a_miss + sep[0] + a_miss);
      pm.alias_phenotype(pm.get_missing_phenotype_id(), a_miss + sep[1] + a_miss);
      pm.alias_phenotype(pm.get_missing_phenotype_id(), a_miss + sep[2] + a_miss);
    }

#if 0
  cout << "var cnt = " << var_info.size() << endl;
  cout << "marker cnt = " << mped_info.marker_count() << endl;
#endif
  }

  if( IS_ASSIGNED(typ) )
  {
    cout << "\n\nReading Type Probability File.........................." << endl;
    cout << "               from " << fname << flush;

    string mfile(fname);

    SAGE::MLOCUS::InheritanceModelFile marker_reader(my_output.errors());

    if( !marker_reader.input(marker_map, mfile, separator, a_miss, my_output.info()) )
    {
      std::cout << std::endl;
      my_output.errors() << SAGE::priority(SAGE::fatal)
                         << "Error reading marker locus description from file '" << mfile << "'." << std::endl;

      exit(EXIT_FAILURE);
    }

    print_fill(fname);
  }


  return true;
}

} // End namespace ONETOOL
