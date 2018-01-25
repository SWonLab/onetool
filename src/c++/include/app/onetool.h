#ifndef __ONETOOL_H
#define __ONETOOL_H

//==========================================================================
// File:      onetool.h
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                              yes Nov 14
//
// Notes:     This header file defines onetool app. derived from app_base.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "app/app_base.h"

namespace ONETOOL {

class onetool : public APP::app_base
{
  public:

    onetool(int argc=0, char **argv=NULL);

    virtual int main();

  protected:

    void run_sage_analyses(app_data& ot_data, cIO* vin, SAGE::cerrormultistream& errors);
	void run_merlin_analyses(app_data& ot_data, cIO* vin, SAGE::cerrormultistream& errors);
};

} // end of namespace ONETOOL

#endif
