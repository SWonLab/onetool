#ifndef SEGREG_ANALYSIS_CONTROLLER_H
#include "sage/segreg/AnalysisController.h"
#endif

namespace SAGE
{
namespace SEGREG
{

//===================================
// AnalysisController inlines
//===================================

/// Simple constructor
///
inline
AnalysisController::AnalysisController(Output_Streams& o)
  : my_out(o),
    my_analysis(my_out)
{ }


}
}
