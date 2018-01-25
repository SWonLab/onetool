//==========================================================================
// File:      app_data.ipp
//
// Author:    Yeunjoo E. Song
//
// History:   Initial implementation                            yes Jan 2015
//
// Notes:     Inline function of app_data.
//
// Copyright (c) 2015 Sungho Won
//   All Rights Reserved
//==========================================================================

namespace ONETOOL {

// ================
// Inline Functions
// ================

inline const SAGE::RPED::RefMultiPedigree&
app_data::pedigrees() const
{
  return my_pedigrees;
}

inline SAGE::RPED::RefMultiPedigree&
app_data::pedigrees()
{
  return my_pedigrees;
}

inline SAGE::Output_Streams&
app_data::get_ostreams() const
{
  return my_output;
}

inline SAGE::cerrorstream&
app_data::errors() const
{
  return my_output.errors();
}

inline ostream&
app_data::info() const
{
  return my_output.info();
}

inline ostream&
app_data::screen() const
{
  return my_output.screen();
}

inline ostream&
app_data::messages() const
{
  return my_output.messages();
}

/*
inline std::vector<std::pair<string, string> >&
app_data::first_ten_ind()
{
  return my_first_ten_ind;
}

inline const std::vector<std::pair<string, string> >&
app_data::first_ten_ind() const
{
  return my_first_ten_ind;
}
*/
} // End namespace ONETOOL

