#ifndef SUBMODEL_H
#include "sage/maxfun/Submodel.h"
#endif

namespace SAGE   {
namespace MAXFUN {

//=======================================================================
//
//  getParameterMgr()
//
//=======================================================================
inline       ParameterMgr * Submodel::getParameterMgr()       { return my_info; }
inline const ParameterMgr * Submodel::getParameterMgr() const { return my_info; }

} // End of namespace MAXFUN
} // End of namespace SAGE

