#pragma once
#ifndef __APP_WISARD_H__
#define __APP_WISARD_H__

#include "global/io.h"

namespace ONETOOL {

cPPPAnalysisV2* getPPP(cIO *io);
int onetool_analysis(cIO *io);

} // End namespace ONETOOL

#endif
