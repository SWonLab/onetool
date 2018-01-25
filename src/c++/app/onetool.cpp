//==========================================================================
// File:      onetool.h
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                              yes Nov 14
//
// Notes:     This source file implements onetool app. derived from app_base.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "app/onetool.h"
#include "app/app_wisard.h"
#include "app/app_sage.h"
#include "app/app_merlin.h"
#include "global/option.h"
#include "input/impute.h"

using namespace std;

namespace ONETOOL {

onetool::onetool(int argc, char** argv)
      : APP::app_base("onetool", argc, argv)
{
  logSet(L_LOG, LOG_BUFFER);
  cOption& user_opt = OPTION();
  user_opt.init(argc, argv);
  /* Init log AGAIN */ {
	  char S_fn[256];
	  sprintf(S_fn, "%s.log", OPT_STRING(out));
	  logSet(L_LOG, S_fn);
  }

  if( user_opt.N_debug )
  {
    print_debug(cerr);
  }

  if( user_opt.N_help )
  {
    print_help(cerr);
    exit(EXIT_SUCCESS);
  }

  // S60715 S_bim is not essential
  //if(   (user_opt.S_vcf == NULL || user_opt.S_fam == NULL)
  //      && user_opt.S_bed == NULL )
  if(    (user_opt.S_vcf == NULL || user_opt.S_fam == NULL)
	  && user_opt.S_fam == NULL
      && (user_opt.S_bed == NULL || user_opt.S_fam == NULL)
      && (user_opt.S_map == NULL || user_opt.S_ped == NULL)
	  && user_opt.S_dosage == NULL)
  {
    if (user_opt.S_fam == NULL) {
      print_help(cerr);
      exit(EXIT_FAILURE);
    }
  }
}

extern cTimer t;
int onetool::main()
{
  app_data ot_data(my_name, debug());
  t.start();

  print_title(ot_data.info());
  print_exe_info(ot_data.info());

  // Perform WISARD analyses
  cIO *io = ot_data.input_wisard();

  onetool_analysis(io);

  // Perform SAGE analyses
  ot_data.input_sage(io);

  SAGE::cerrormultistream errors;
  errors.insert(ot_data.errors());

  run_sage_analyses(ot_data, io, errors);
  
  // Perform Merlin analyses
  run_merlin_analyses(ot_data, io, errors);

  cout << endl << "Analysis complete!" << endl << endl;

  return EXIT_SUCCESS;
}

} // end of namespace ONETOOL

int main(int argc, char* argv[])
{
  free(malloc(1));

  ONETOOL::onetool onetool_inst(argc, argv);

  onetool_inst.main();

  exit(EXIT_SUCCESS);
}
