//==========================================================================
// File:      app_merlin.cpp
//
// Author:    Sungyoung Lee
//
// History:   Initial implementation                              syl Nov 22
//
// Notes:     This source file implements Merlin app analyses.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "app/onetool.h"
#include "app/app_merlin.h"

using namespace std;

namespace ONETOOL {

void
onetool::run_merlin_analyses(app_data& ot_data, cIO* vin, SAGE::cerrormultistream& errors)
{	
  if( !(OPT_ENABLED(merlin) && IS_ASSIGNED(map)) )
    return;

  cout << endl << endl
       << "Prepare MERLIN NPL analysis............................" << flush;

  int N_argc = 0;
  char** Sa_argv = NULL;
  if (OPT_ENABLED(merlin))
  {
    N_argc = 12;
    Sa_argv = new char*[12];
    const char* Sa_inp[] = { "__dummy__", "-d", "onetool", "-p", "onetool", "-m", "onetool",
                             "--pairs", "--npl", "--tabulate", "--prefix" };
    for (int i=0 ; i<11 ; i++)
      Sa_argv[i] = strdup(Sa_inp[i]);
    Sa_argv[11] = OPT_STRING(out);
  }
  ::main_ONETOOL(vin, N_argc, Sa_argv);

  // Cleanup
  if (Sa_argv)
  {
    for (int i=0 ; i<11 ; i++)
      free(Sa_argv[i]);
    delete [] Sa_argv;
  }

  return;
}
	
} // end of namespace ONETOOL
