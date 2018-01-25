#ifndef __APP_BASE_H
#define __APP_BASE_H

//==========================================================================
// File:      app_base.h
//
// Author:    Yeunjoo E. Song
//            Sungyoung Lee
//
// History:   Initial implementation                              yes Nov 14
//            Windows compatibility update                        syl Dec 12
//
// Notes:     This header file defines the base structure for pinfoqc program.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include "app/app_version.h"
#include "sage/global/output_streams.h"
#include "utils/util.h"
#include "input/app_data.h"

namespace ONETOOL {
namespace APP    {

class app_base
{
  public:

    app_base(string name, int argc=0, char **argv=NULL);

    virtual int main() = 0;

    virtual inline ~app_base();

    void   print_title(std::ostream& o);
    void   print_exe_info(std::ostream& o);
    void   print_debug(std::ostream& o);
    void   print_help(std::ostream& o);

    int    get_argc() { return argc; }
    char** get_argv() { return argv; }

    inline bool debug() const { return !!dbg; }

    static string get_release_string();

    static const char *build_date;
    static const char *cxxflags;
    static const char *ldflags;
    static const char *ldlibs;
    static const char *build;

  protected:

    int             argc;
    char**          argv;
    string          my_name;

  private:

    void calculate_release_string() const;

    static string release_string;

    int dbg;
};

char getItselfPath(char *Sp_path);

//==================================================================================
// INLINE FUNCTIONS
//==================================================================================

inline app_base::~app_base()
{
  argv = NULL;
}

} // End namespace APP
} // End namespace ONETOOL

#endif
