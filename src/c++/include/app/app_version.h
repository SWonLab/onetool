#ifndef __APP_VERSION_NUMBER_H
#define __APP_VERSION_NUMBER_H
//==========================================================================
// File:      app_version.h
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                              yes Nov 14
//
// Notes:     This header file defines the structure for app version number.
//            A version is composed of four numbers: main version,
//            sub version, micro version, and beta version (where
//            a beta version of 0 indicates that the version is not
//            a beta version).
//
//            Examples: 1.2.3.4, 5.0.1.0, 2.2.3.119, and so on.
//
// Copyright (c) 2014 Sungho Won
//   All Rights Reserved
//==========================================================================

#include <string>
#include <sstream>

namespace ONETOOL {
namespace APP    {

class app_version
{
  public:

    app_version(int main, int sub, int micro, int beta);

    ///
    /// Renders this version in pretty print.
    /// (eg: (non-beta) "INFOQC v1.2.3", (beta) "INFOQC 3.2.4 Beta 4")
    std::string print_pretty() const;

    /// Indicates whether or not this is a beta version (beta version != 0, that is).
    bool is_beta() const;

    /// Returns the main version.
    int get_main_version() const;

    /// Returns the sub version.
    int get_sub_version() const;

    /// Returns the micro version.
    int get_micro_version() const;

    /// Returns the beta version (0 = not beta)
    int get_beta_version() const;

  private:

    int my_main_version;
    int my_sub_version;
    int my_micro_version;
    int my_beta_version;
};

extern app_version version;

inline int app_version::get_main_version  () const { return my_main_version;  }
inline int app_version::get_sub_version   () const { return my_sub_version;   }
inline int app_version::get_micro_version () const { return my_micro_version; }
inline int app_version::get_beta_version  () const { return my_beta_version;  }

inline
app_version::app_version(int main, int sub, int micro, int beta)
{
  my_main_version  = main;
  my_sub_version   = sub;
  my_micro_version = micro;
  my_beta_version  = beta;
}

inline bool app_version::is_beta() const { return my_beta_version != 0; }

inline std::string
app_version::print_pretty() const
{
  std::ostringstream s;

  s << "v"
    << my_main_version << "."
    << my_sub_version << "."
    << my_micro_version;

  if( is_beta() )
    s << " Beta " << my_beta_version;

  return s.str();
}

} // End namespace APP
} // End namespace ONETOOL

#endif
