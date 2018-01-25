#ifndef __APP_SAGE_H
#define __APP_SAGE_H

//==========================================================================
// File:      app_sage.h
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                              yes Nov 14
//
// Notes:     This header file defines S.A.G.E. app analyses.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "app/onetool.h"
#include "sage/rped/rpfile.h"
#include "sage/pedinfo/stats_view.h"
#include "sage/fcor/analysis.h"
#include "sage/segreg/AnalysisController.h"
#include "sage/lodlink/analysis.h"

namespace ONETOOL {

//void run_sage_analyses(app_data& ot_data, cIO* vin, SAGE::cerrormultistream& errors);
void parse_fcor_options(SAGE::FCOR::analysis_option_type& fcor_opt);
void parse_lodlink_options(SAGE::LODLINK::instructions& lodlink_inst,
                           const SAGE::RPED::RefMultiPedigree& mp, SAGE::cerrorstream& errors);
void parse_segreg_options(SAGE::SEGREG::model& segreg_mo,
                           const SAGE::RPED::RefMultiPedigree& mp, SAGE::cerrorstream& errors);

} // end of namespace ONETOOL

#endif
